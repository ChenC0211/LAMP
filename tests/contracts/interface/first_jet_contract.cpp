#include "interface_first_jet_fixture.h"
#include "Core/Geometry/ActiveQuarticHermiteCompletion.h"
#include <cstdint>
#include <map>
#include <tuple>
namespace {

using MacroCompletion=LAMP::Geometry::ActiveQuarticHermiteCompletion;
using MacroSensitivity=LAMP::Geometry::ActiveQuarticHermiteSensitivity;
using MacroSource=std::tuple<int,int,int>;
std::array<double,6> macroJetChannels(const Snapshot::ScalarJet2& j)
{return {j.value,j.dx,j.dy,j.dxx,j.dxy,j.dyy};}
std::unique_ptr<QuadtreeMesh> macroMesh(int nx,int ny)
{
    DictionaryMesh d;d.Nx=nx;d.Ny=ny;d.maxCellNum=1024;d.maxInitialIterations=0;
    d.minCoord={.125,-.25};d.maxCoord={1.875,1.375};
    return std::make_unique<QuadtreeMesh>(d);
}
long double macroPolynomial(Vector2D q,int dx,int dy)
{
    long double result=0;
    for(int j=dy;j<5;++j)for(int i=dx;i<5;++i){
        long double factor=((i+j)%2?-1.L:1.L)*(1+(i+1)*(j+2))/200;
        for(int k=0;k<dx;++k)factor*=i-k;
        for(int k=0;k<dy;++k)factor*=j-k;
        result+=factor*std::pow(static_cast<long double>(q.x),i-dx)*std::pow(static_cast<long double>(q.y),j-dy);
    }
    return result;
}
Patch macroPatch(const QuadtreeMesh& mesh,bool arbitrary)
{
    Patch p;p.activeCells.assign(mesh.cellIndexNum(),1);p.nodes.resize(mesh.vertIndexNum());p.cellCenters.resize(mesh.cellIndexNum());
    std::uint32_t seed=137713;
    const auto random=[&]{seed=1664525u*seed+1013904223u;return static_cast<double>(seed)/4294967296.0-.5;};
    for(int v=0;v<mesh.vertIndexNum();++v){const auto q=point(mesh,v);
        p.nodes[v]={static_cast<double>(macroPolynomial(q,0,0)),static_cast<double>(macroPolynomial(q,1,0)),
            static_cast<double>(macroPolynomial(q,0,1)),static_cast<double>(macroPolynomial(q,1,1))};
        if(arbitrary)p.nodes[v]={random(),random(),random(),random()};
    }
    std::set<Key> edges;
    for(int c=0;c<mesh.cellIndexNum();++c){const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);
        p.cellCenters[c]=arbitrary?random():static_cast<double>(macroPolynomial(cell.centerCoord(),0,0));
        for(int k=0;k<4;++k){Key key{mesh.verts(cell.cellToVert[k]).index(),mesh.verts(cell.cellToVert[(k+1)%4]).index()};
            const auto a=point(mesh,key[0]),b=point(mesh,key[1]);if(a.y==b.y?a.x>b.x:a.y>b.y)std::swap(key[0],key[1]);edges.insert(key);}}
    for(const auto& key:edges){const auto a=point(mesh,key[0]),b=point(mesh,key[1]),q=a+(b-a)*.5;const bool horizontal=a.y==b.y;
        p.edgeMidpoints.push_back({key,arbitrary?random():static_cast<double>(macroPolynomial(q,0,0)),
            arbitrary?random():static_cast<double>(macroPolynomial(q,horizontal?0:1,horizontal?1:0))});}
    std::reverse(p.edgeMidpoints.begin(),p.edgeMidpoints.end());
    return p;
}
// Independent ten-condition polynomial solve. No production seven-basis,
// piece transform, evaluation helper or Bernstein conversion is reused.
struct MacroReference {
    using Basis=std::array<std::array<std::array<long double,7>,5>,2>;
    Basis basis{};
    std::array<std::array<long double,25>,4> coefficients{};
    long double hx,hy;
    static long double derivative(int power,int order,long double x)
    {if(power<order)return 0;long double factor=1;for(int k=0;k<order;++k)factor*=power-k;
        return factor*std::pow(x,power-order);}
    MacroReference(const std::array<double,49>& physical,Vector2D widths):hx(widths.x),hy(widths.y)
    {
        for(int datum=0;datum<7;++datum){
            long double a[10][11]{};
            const auto condition=[&](int row,int half,int order,long double x,long double scale){
                for(int k=0;k<5;++k)a[row][5*half+k]=scale*derivative(k,order,x);};
            condition(0,0,0,0,1);condition(1,0,0,1,1);condition(2,1,0,1,1);
            condition(3,0,1,0,2);condition(4,1,1,1,2);condition(5,0,2,0,4);condition(6,1,2,1,4);
            for(int order=0;order<3;++order){condition(7+order,0,order,1,1);condition(7+order,1,order,0,-1);}
            a[datum][10]=1;
            for(int column=0;column<10;++column){int pivot=column;
                for(int row=column+1;row<10;++row)if(std::abs(a[row][column])>std::abs(a[pivot][column]))pivot=row;
                check(a[pivot][column]!=0,"singular independent seven-datum macro system");
                for(int k=column;k<11;++k)std::swap(a[pivot][k],a[column][k]);
                const auto value=a[column][column];for(int k=column;k<11;++k)a[column][k]/=value;
                for(int row=0;row<10;++row)if(row!=column){const auto factor=a[row][column];
                    for(int k=column;k<11;++k)a[row][k]-=factor*a[column][k];}}
            for(int half=0;half<2;++half)for(int k=0;k<5;++k)basis[half][k][datum]=a[5*half+k][10];
        }
        constexpr int order[7]={0,0,0,1,1,2,2};
        for(int py=0;py<2;++py)for(int px=0;px<2;++px)for(int y=0;y<5;++y)for(int x=0;x<5;++x)
            for(int j=0;j<7;++j)for(int i=0;i<7;++i)
                coefficients[px+2*py][x+5*y]+=physical[i+7*j]*std::pow(hx,order[i])*std::pow(hy,order[j])
                    *basis[px][x][i]*basis[py][y][j];
    }
    std::array<long double,6> jet(int piece,Vector2D local) const
    {
        constexpr int dx[6]={0,1,0,2,1,0},dy[6]={0,0,1,0,1,2};std::array<long double,6> result{};
        for(int c=0;c<6;++c)for(int y=0;y<5;++y)for(int x=0;x<5;++x)
            result[c]+=coefficients[piece][x+5*y]*derivative(x,dx[c],local.x)*derivative(y,dy[c],local.y)
                *std::pow(2/hx,dx[c])*std::pow(2/hy,dy[c]);
        return result;
    }
};
MacroSource macroSourceKey(const MacroSensitivity& source)
{return {static_cast<int>(source.kind),source.index,source.component};}
std::vector<MacroSource> macroSources(const Patch& p)
{
    std::vector<MacroSource> result;
    for(int v=0;v<static_cast<int>(p.nodes.size());++v)for(int k=0;k<4;++k)
        result.emplace_back(static_cast<int>(MacroSensitivity::Kind::Vertex),v,k);
    for(int e=0;e<static_cast<int>(p.edgeMidpoints.size());++e)for(int k=0;k<2;++k)
        result.emplace_back(static_cast<int>(MacroSensitivity::Kind::EdgeMidpoint),e,k);
    for(int c=0;c<static_cast<int>(p.cellCenters.size());++c)result.emplace_back(static_cast<int>(MacroSensitivity::Kind::CellCenter),c,0);
    return result;
}
void macroSetSource(Patch& p,MacroSource source,double value)
{
    const auto [kind,index,component]=source;
    if(kind==static_cast<int>(MacroSensitivity::Kind::Vertex)){
        auto& node=p.nodes.at(index);if(component==0)node.value=value;else if(component==1)node.dx=value;
        else if(component==2)node.dy=value;else node.dxy=value;
    }else if(kind==static_cast<int>(MacroSensitivity::Kind::EdgeMidpoint)){
        auto& edge=p.edgeMidpoints.at(index);if(component==0)edge.value=value;else edge.normalDerivative=value;
    }else p.cellCenters.at(index)=value;
}
void c2MacroLocalContract()
{
    constexpr int dx[6]={0,1,0,2,1,0},dy[6]={0,0,1,0,1,2};
    int cells=0,interpolations=0,seams=0,reproductions=0,pullbacks=0,neighborCenters=0,guards=0;
    long double maximum=0;
    for(const auto dimensions:{std::array<int,2>{2,2},std::array<int,2>{4,3}}){
        auto mesh=macroMesh(dimensions[0],dimensions[1]);const auto patch=macroPatch(*mesh,true);
        const auto completion=MacroCompletion::fromUniformPatch(*mesh,patch);
        auto scalar=Snapshot::fromUniformC2Completion(*mesh,completion);
        check(scalar.usesC2Macro()&&scalar.reconstruction()==LevelSetHermiteGeometry::Reconstruction::C2MacroQuartic,
            "macro snapshot representation metadata");
        const auto compare=[&](const auto& a,const auto& b,Vector2D widths,const std::string& why){
            for(int c=0;c<6;++c){const long double scale=std::pow(static_cast<long double>(widths.x),dx[c])
                *std::pow(static_cast<long double>(widths.y),dy[c]);
                maximum=std::max(maximum,close(scale*a[c],scale*b[c],3e-11,why));}
        };
        for(int cell=0;cell<mesh->cellIndexNum();++cell){const auto box=scalar.nativeCellBounds(cell);const auto widths=box[1]-box[0];
            const auto data=scalar.nativeC2MacroPhysicalDofs(cell);MacroReference reference(data,widths);
            for(int piece=0;piece<4;++piece)for(Vector2D q:{Vector2D{0,0},Vector2D{1,1},Vector2D{.23,.71}}){
                compare(macroJetChannels(scalar.scalarJetAtNativeC2MacroPieceLocal(cell,piece,q)),reference.jet(piece,q),widths,
                    "independent arbitrary 49-datum piece jet");}
            for(const auto region:{std::array<Vector2D,2>{Vector2D{0,0},Vector2D{1,1}},
                    std::array<Vector2D,2>{Vector2D{.13,.21},Vector2D{.82,.94}},
                    std::array<Vector2D,2>{Vector2D{.5,.5},Vector2D{.5,.5}}}){
                for(int channel=0;channel<6;++channel){
                    const auto enclosure=scalar.nativeC2MacroScalarBounds(cell,region[0],region[1],dx[channel],dy[channel]);
                    check(std::isfinite(enclosure[0])&&std::isfinite(enclosure[1])&&enclosure[0]<=enclosure[1],"invalid macro physical derivative enclosure");
                    for(double tx:{0.,.23,.5,1.})for(double ty:{0.,.67,1.}){
                        const Vector2D q{region[0].x+tx*(region[1].x-region[0].x),region[0].y+ty*(region[1].y-region[0].y)};
                        const int px=q.x>.5,py=q.y>.5;const auto exact=reference.jet(px+2*py,{2*q.x-px,2*q.y-py})[channel];
                        const long double scale=std::pow(static_cast<long double>(widths.x),dx[channel])*std::pow(static_cast<long double>(widths.y),dy[channel]);
                        check(scale*(exact-enclosure[0])>=-3e-11L&&scale*(enclosure[1]-exact)>=-3e-11L,"independent macro derivative outside outward bound");
                        if(region[0].x==region[1].x&&region[0].y==region[1].y)
                            check(scale*(enclosure[1]-enclosure[0])<=3e-11L,"macro point derivative enclosure is unresolved");
                    }
                }
            }
            constexpr double locations[5]={0,.5,1,0,1};constexpr int orders[5]={0,0,0,1,1};
            for(int y=0;y<5;++y)for(int x=0;x<5;++x){const auto actual=macroJetChannels(scalar.scalarJetAtNativeCellChart(cell,
                    LAMP::Geometry::cellChartFromLocal({locations[x],locations[y]})));
                const int channel=orders[x]?orders[y]?4:1:orders[y]?2:0;
                const long double scale=std::pow(static_cast<long double>(widths.x),orders[x])
                    *std::pow(static_cast<long double>(widths.y),orders[y]);
                maximum=std::max(maximum,close(scale*actual[channel],scale*data[x+7*y],3e-11,"macro original 25 interpolation"));++interpolations;}
            for(double t:{0.,.19,.5,.83,1.})for(int transverse=0;transverse<2;++transverse){
                compare(macroJetChannels(scalar.scalarJetAtNativeC2MacroPieceLocal(cell,2*transverse,{1,t})),
                    macroJetChannels(scalar.scalarJetAtNativeC2MacroPieceLocal(cell,2*transverse+1,{0,t})),widths,"macro vertical half seam");
                compare(macroJetChannels(scalar.scalarJetAtNativeC2MacroPieceLocal(cell,transverse,{t,1})),
                    macroJetChannels(scalar.scalarJetAtNativeC2MacroPieceLocal(cell,transverse+2,{t,0})),widths,"macro horizontal half seam");seams+=2;}
            for(int other=cell+1;other<mesh->cellIndexNum();++other){const auto b=scalar.nativeCellBounds(other);
                if(box[1].x==b[0].x&&box[0].y==b[0].y&&box[1].y==b[1].y)
                    for(double t:{0.,.19,.5,.83,1.}){compare(macroJetChannels(scalar.scalarJetAtNativeCellChart(cell,LAMP::Geometry::cellChartFromLocal({1,t}))),
                        macroJetChannels(scalar.scalarJetAtNativeCellChart(other,LAMP::Geometry::cellChartFromLocal({0,t}))),widths,"macro native vertical seam");++seams;}
                if(box[1].y==b[0].y&&box[0].x==b[0].x&&box[1].x==b[1].x)
                    for(double t:{0.,.19,.5,.83,1.}){compare(macroJetChannels(scalar.scalarJetAtNativeCellChart(cell,LAMP::Geometry::cellChartFromLocal({t,1}))),
                        macroJetChannels(scalar.scalarJetAtNativeCellChart(other,LAMP::Geometry::cellChartFromLocal({t,0}))),widths,"macro native horizontal seam");++seams;}
            }++cells;
        }
        // Independent analytic tensor Q4 input, not merely a smooth circle.
        const auto polynomial=macroPatch(*mesh,false);const auto reproduced=MacroCompletion::fromUniformPatch(*mesh,polynomial);
        auto smooth=Snapshot::fromUniformC2Completion(*mesh,reproduced);
        for(int cell=0;cell<mesh->cellIndexNum();++cell){const auto box=smooth.nativeCellBounds(cell);const auto widths=box[1]-box[0];
            for(int piece=0;piece<4;++piece)for(Vector2D q:{Vector2D{0,0},Vector2D{1,1},Vector2D{.31,.67}}){
                const Vector2D global{box[0].x+(piece%2+q.x)*widths.x/2,box[0].y+(piece/2+q.y)*widths.y/2};
                std::array<long double,6> expected{};for(int c=0;c<6;++c)expected[c]=macroPolynomial(global,dx[c],dy[c]);
                compare(macroJetChannels(smooth.scalarJetAtNativeC2MacroPieceLocal(cell,piece,q)),expected,widths,"macro analytic Q4 reproduction");++reproductions;}}
        // All source columns of A/E at a fixed actual native chart. Unit
        // inputs remove subtraction/cancellation and retain nonlocal centers.
        const int selected=0;const auto box=scalar.nativeCellBounds(selected);const auto widths=box[1]-box[0];
        for(Vector2D q:{Vector2D{.23,.71},Vector2D{.5,.5},Vector2D{1,.37}}){
            const auto chart=LAMP::Geometry::cellChartFromLocal(q);std::array<std::map<MacroSource,long double>,6> weights;
            for(int channel=0;channel<6;++channel){std::array<long double,6> selector{};selector[channel]=1;
                for(const auto& w:scalar.nativeC2MacroDifferentialSensitivity(selected,chart,selector)){
                    check(std::isfinite(w.coefficient)&&w.coefficient!=0&&weights[channel].emplace(macroSourceKey(w),w.coefficient).second,
                        "macro pullback invalid or unmerged original source");
                    if(w.kind==MacroSensitivity::Kind::CellCenter&&w.index!=selected)++neighborCenters;}}
            const int px=q.x>.5,py=q.y>.5,piece=px+2*py;const Vector2D local{2*q.x-px,2*q.y-py};
            auto zero=patch;for(const auto& source:macroSources(zero))macroSetSource(zero,source,0);
            for(const auto& source:macroSources(patch)){auto unit=zero;macroSetSource(unit,source,1);
                const auto completed=MacroCompletion::fromUniformPatch(*mesh,unit);MacroReference reference(completed.physicalCellDofs(selected),widths);
                std::array<long double,6> actual{};for(int channel=0;channel<6;++channel)actual[channel]=weights[channel][source];
                compare(actual,reference.jet(piece,local),widths,"independent original-source E/A point transpose");++pullbacks;}
        }
        std::array<long double,25> legacy{};rejects([&]{scalar.nativeQuarticCellSensitivity(0,legacy);},"25-column native sensitivity");++guards;
        std::array<long double,6> differential{};rejects([&]{scalar.nativeQuarticDifferentialSensitivity(0,{},differential);},"Q4 differential sensitivity");++guards;
        differential[0]=std::numeric_limits<long double>::quiet_NaN();
        rejects([&]{scalar.nativeC2MacroDifferentialSensitivity(0,{},differential);},"nonfinite C2 macro differential weight");++guards;
        rejects([&]{scalar.scalarJetAtNativeC2MacroPieceLocal(0,4,{});},"piece index");++guards;
        rejects([&]{scalar.nativeC2MacroScalarBounds(0,{-.001,0},{1,1});},"leave the actual native cell");
        rejects([&]{scalar.nativeC2MacroScalarBounds(0,{0,0},{1,1},5,0);},"derivative order exceeds Q4");
        auto original=Snapshot::fromNativeLeaves(*mesh,patch);
        differential={};rejects([&]{original.nativeC2MacroDifferentialSensitivity(0,{},differential);},"requires a completed four-piece");++guards;
    }
    check(cells==16&&interpolations==400&&seams>300&&reproductions==192&&pullbacks>500&&neighborCenters>0&&guards==10,
        "macro mathematical coverage incomplete");
    std::cout<<"C2_MACRO_B,PASS,cells="<<cells<<",original25="<<interpolations<<",seams="<<seams
        <<",Q4="<<reproductions<<",source_transposes="<<pullbacks<<",nonlocal_center="<<neighborCenters
        <<",guards="<<guards<<",max_scaled_error="<<maximum<<'\n';
}


void c2MacroRootContract()
{
    int queries=0,variations=0,sideRanges=0;
    for(const auto size:{std::array<int,2>{2,2},std::array<int,2>{4,4}}){
        auto mesh=macroMesh(size[0],size[1]);
        for(int axis=0;axis<2;++axis){
            const auto lo=mesh->minCoord,hi=mesh->maxCoord;
            const double low=axis?lo.y:lo.x,high=axis?hi.y:hi.x;
            const double line=low+(high-low)/(2*size[axis]),fixed=axis?.71:.19;
            auto patch=macroPatch(*mesh,false);
            const auto value=[&](Vector2D q){return (axis?q.y:q.x)-line;};
            for(int v=0;v<mesh->vertIndexNum();++v)patch.nodes[v]={value(point(*mesh,v)),axis?0.:1.,axis?1.:0.,0};
            for(auto& e:patch.edgeMidpoints){const auto a=point(*mesh,e.vertices[0]),b=point(*mesh,e.vertices[1]);
                e.value=value(a+(b-a)*.5);e.normalDerivative=a.y==b.y?(axis?1.:0.):(axis?0.:1.);}
            for(int c=0;c<mesh->cellIndexNum();++c)patch.cellCenters[c]=value(mesh->cells(mesh->cellIndexToLabel[c]).centerCoord());
            const auto completed=MacroCompletion::fromUniformPatch(*mesh,patch);
            auto source=std::make_shared<const Snapshot>(Snapshot::fromUniformC2Completion(*mesh,completed));
            const Vector2D first=axis?Vector2D{fixed,low}:Vector2D{low,fixed},last=axis?Vector2D{fixed,high}:Vector2D{high,fixed};
            for(bool reverse:{false,true}){
                LevelSetMappedSegment segment(source,reverse?last:first,reverse?first:last);
                check(segment.usesSource(*source)&&segment.gridPieceCount()==static_cast<std::size_t>(2*size[axis]),"C2 segment misses actual half cuts");
                long double negative=0,total=0;for(const auto& range:segment.ranges()){total+=range.length;if(range.negative)negative+=range.length;}
                close(total,1,1e-12,"C2 segment physical range coverage");close(negative,(line-low)/(high-low),1e-12,"C2 half-seam root negative measure");
                const auto slopes=segment.nativeQuarticNegativePhaseVariation();check(slopes.size()==1&&!slopes.front().endpoint,"C2 root half incidence not uniquely accumulated");
                long double constantSlope=0;for(const auto& weight:slopes.front().independent)if(weight.component==0)constantSlope+=weight.coefficient;
                close(constantSlope,-1/(high-low),1e-11,"C2 full original-source root transpose");++variations;
                const long double sample=reverse?1-(line-low)/(high-low):(line-low)/(high-low);
                const auto phi=segment.nativeQuarticValueVariation(false,sample);close(phi.value,0,1e-12,"C2 root value variation actual source");
                long double constant=0;for(const auto& weight:phi.independent)if(weight.component==0)constant+=weight.coefficient;
                close(constant,1,1e-11,"C2 original-source point partition of unity");++variations;
                auto trace=LevelSetMappedSegment::nativeQuarticTrace(source,reverse?last:first,reverse?first:last,1-axis);
                check(trace.usesSource(*source)&&trace.coincidentRanges().empty()&&trace.gridPieceCount()==segment.gridPieceCount(),"C2 ordinary trace changed geometry");
                close(trace.ordinary().ranges().front().length,segment.ranges().front().length,1e-12,"C2 trace/query same range");++queries;
            }
            const Vector2D sideFirst=axis?Vector2D{lo.x,line}:Vector2D{line,lo.y},sideLast=axis?Vector2D{hi.x,line}:Vector2D{line,hi.y};
            auto coincident=LevelSetMappedSegment::nativeQuarticTrace(source,sideFirst,sideLast,axis);
            check(coincident.usesSource(*source)&&coincident.ordinaryRanges().empty()&&!coincident.coincidentRanges().empty(),"C2 exact halfline trace lost typed sides");
            for(const auto& range:coincident.coincidentRanges()){
                check(range.minus.available&&range.plus.available&&range.minus.normalOrder==1&&range.plus.normalOrder==1,
                    "C2 shared halfline normal Taylor source missing");
                for(const auto& r:range.minus.ranges)check(r.negative,"C2 minus halfline phase");
                for(const auto& r:range.plus.ranges)check(!r.negative,"C2 plus halfline phase");++sideRanges;}
            ++queries;
        }
    }
    check(queries==12&&variations==16&&sideRanges>0,"C2 restriction/root coverage incomplete");
    std::cout<<"C2_MACRO_C1,PASS,queries="<<queries<<",original_source_variations="<<variations<<",typed_halfline_sides="<<sideRanges<<'\n';
}

void firstJetPublicationContract()
{
    int cases=0,guards=0;
    for(bool radial:{false,true})for(bool amr:{false,true})for(bool nonpolynomial:{false,true}){
        auto mesh=makeMesh(8,amr?1:0);Case s;s.pureFirstJet=true;s.radial=radial;s.amr=amr;s.nonpolynomial=nonpolynomial;
        if(radial)s.radius=.27;
        const auto initialGeometry=initial(*mesh,s);auto p=initialGeometry.state().primary();
        check(p.transport==HermiteGalsTransport::FirstJet && p.phaseVolumes.empty(),"pure initial state contains M");
        HermiteGalsCandidate candidate(*mesh,{},p);const auto id=candidate.identity();
        const auto mixed=closeHermiteGalsMixedAmr(*mesh,candidate,s.symmetry(*mesh));
        HermiteGalsCenterPrediction predicted{id,centers(*mesh,s.profile(0))};
        const auto reject=[&](auto f,const std::string& message){rejects(f,message);++guards;};
        auto bad=p;bad.phaseVolumes.push_back({p.geometryCells.front(),0});
        reject([&]{HermiteGalsCandidate invalid(*mesh,{},bad);},"FirstJet state cannot carry conserved inventory");
        bad=p;bad.transport=static_cast<HermiteGalsTransport>(-1);
        reject([&]{HermiteGalsCandidate invalid(*mesh,{},bad);},"unknown transport mode");
        bad=p;bad.transport=HermiteGalsTransport::ConservativeVolume;
        reject([&]{HermiteGalsCandidate invalid(*mesh,{},bad);},"geometry cell lacks conserved inventory");
        HermiteGalsCandidate other(*mesh,{},p);const auto foreign=closeHermiteGalsMixedAmr(*mesh,other,s.symmetry(*mesh));
        auto wrong=predicted;wrong.primaryIdentity=other.identity();
        reject([&]{publishHermiteGalsFirstJetGeometry(*mesh,std::move(candidate),mixed,wrong);},"FirstJet primary revision mismatch");
        reject([&]{publishHermiteGalsFirstJetGeometry(*mesh,std::move(candidate),foreign,predicted);},"FirstJet primary revision mismatch");
        wrong=predicted;wrong.values.pop_back();
        reject([&]{publishHermiteGalsFirstJetGeometry(*mesh,std::move(candidate),mixed,wrong);},"FirstJet closure layout mismatch");
        wrong=predicted;wrong.values.front()=std::numeric_limits<double>::quiet_NaN();
        reject([&]{publishHermiteGalsFirstJetGeometry(*mesh,std::move(candidate),mixed,wrong);},"FirstJet nonfinite center prediction");
        HermiteGalsPublicationConfig invalidConfig;invalidConfig.unitNormalTolerance=0;
        reject([&]{publishHermiteGalsFirstJetGeometry(*mesh,std::move(candidate),mixed,predicted,invalidConfig);},"invalid publication tolerance");
        reject([&]{closeHermiteGalsBubble(*mesh,candidate,mixed,{id,predicted.values});},"inventory requires complete scalar support");
        check(candidate.identity()==id && !candidate.hasGeometry() && candidate.primary().phaseVolumes.empty(),"pure publication failure mutated candidate");
        const auto pub=publishHermiteGalsFirstJetGeometry(*mesh,std::move(candidate),mixed,predicted);
        reject([&]{candidate.identity();},"moved-from or consumed handle");
        reject([&]{pub.volumes().initializeCandidate(*mesh);},"FirstJet state cannot carry conserved inventory");
        check(pub.state().identity()==id && pub.state().timeLayer().time()==0 && pub.state().primary().phaseVolumes.empty(),"pure publication authority");
        check(pub.state().closure().cellBubbleCoefficients==predicted.values
            && pub.state().closure().vertexMixedDerivatives==mixed.vertexMixedDerivatives(),"pure publication adjusted closure");
        const auto& scalar=pub.state().geometry().snapshot();
        check(&pub.volumes().geometry().snapshot()==&scalar && &pub.facePartitions().source()==&scalar,"pure geometry has multiple owners");
        check(std::isnan(pub.report().maxRelativeInventoryResidual) && std::isnan(pub.volumes().report().totalSourceInventory)
            && std::isnan(pub.volumes().report().geometricInventoryDifference),"pure mode claims inventory residual");
        check(pub.volumes().cells().size()==p.geometryCells.size() && pub.report().quadratureNodes>0,"pure geometric diagnostics missing");
        for(std::size_t k=0;k<p.geometryCells.size();++k){const int cell=p.geometryCells[k];
            close(scalar.scalarJetInNativeCell(cell,mesh->cells(mesh->cellIndexToLabel[cell]).centerCoord()).value,
                predicted.values[k],2e-15,"pure center interpolation");}
        for(std::size_t k=0;k<p.vertices.size();++k){const auto& a=p.vertices[k];const auto& b=pub.state().primary().vertices[k];
            check(a.phi==b.phi && a.dx==b.dx && a.dy==b.dy,"pure publication changed master first jet");}
        const auto e=errors(pub,s);if(!nonpolynomial)check(std::max({e.position,e.gradient,e.curvature})<=1e-10,"pure Q2 geometry reproduction");
        std::cout<<"FIRST_JET_PUBLICATION,"<<radial<<','<<amr<<','<<nonpolynomial<<','<<e.position<<','<<e.gradient<<','<<e.curvature<<','<<e.volume<<",PASS\n";
        ++cases;
    }
    std::cout<<"I12p1 FirstJet publication PASS cases="<<cases<<" guards="<<guards<<'\n';
}
void firstJetStepContract()
{
    int cases=0,guards=0;
    const auto reject=[&](auto f,const std::string& why){rejects(f,why);++guards;};
    for(bool radial:{false,true})for(bool amr:{false,true})for(bool nonpolynomial:{false,true}){
        auto mesh=makeMesh(8,amr?1:0);Case s;s.pureFirstJet=true;s.radial=radial;s.amr=amr;s.nonpolynomial=nonpolynomial;
        s.motion=1;if(radial)s.radius=.27;
        const auto source=initial(*mesh,s);const auto id=source.state().identity();const auto& scalar=source.state().geometry().snapshot();
        const auto domain=traceDomain(*mesh,radial);HermiteGalsPublicationConfig config;config.volume.symmetryPlanes=s.symmetry(*mesh);
        const HermiteGalsTimeLayer time{1,0,.01,Position::End};HermiteGalsCandidate authority(*mesh,time,source.state().primary());
        const auto velocity=HermiteGalsTemporalVelocitySource::constant(source.state(),authority.identity(),time,s.velocity({},0).value);
        std::size_t queried=0;const HermiteGalsBoundaryFirstJetQuery boundary=[&](Vector2D q,HermiteGalsTimeLayer t,Boundary){
            check(q.x==mesh->minCoord.x || q.x==mesh->maxCoord.x || q.y==mesh->minCoord.y || q.y==mesh->maxCoord.y,"boundary query used for interior reset");
            check(t.step==1 && t.begin==0 && t.end==.01,"boundary query lost target clock");++queried;
            const auto j=s.exact(q,t.time());return HermiteGalsBoundaryFirstJet{j.v,{j.x,j.y}};};
        for(Position position:{Position::Midpoint,Position::End}){
            const HermiteGalsTimeLayer target{1,0,.01,position};const auto arrival=velocity.constantArrivalView(authority.identity(),target);
            queried=0;auto result=prepareHermiteGalsFirstJetStep(*mesh,source,domain,velocity.previous(),arrival,boundary,config);
            const auto& pub=result.publication();const auto& r=result.report();const auto& p=pub.state().primary();
            check(result.sourceIdentity()==id && result.arrivalVelocityIdentity()==authority.identity() && pub.state().identity()!=id
                && pub.state().identity()!=authority.identity() && pub.state().timeLayer().time()==target.time(),"step provenance/time authority");
            check(p.transport==HermiteGalsTransport::FirstJet && p.phaseVolumes.empty() && r.tracedCenters==p.geometryCells.size(),"step introduced inventory or lost centers");
            check(r.transportedVertices+r.prescribedVertices==p.vertices.size() && r.transportedEdges+r.prescribedEdges==p.edges.size()
                && r.prescribedVertices && r.prescribedEdges && queried==r.prescribedVertices+r.prescribedEdges,"step boundary/transport coverage");
            check(&pub.volumes().geometry().snapshot()==&pub.state().geometry().snapshot()
                && &pub.facePartitions().source()==&pub.state().geometry().snapshot(),"step has multiple scalar owners");
            const auto e=errors(pub,s);if(!nonpolynomial)check(std::max({e.position,e.gradient,e.curvature})<=1e-10 && e.cellVolume<=5e-12,"step Q2 accuracy gate");
            std::cout<<"FIRST_JET_STEP,"<<radial<<','<<amr<<','<<nonpolynomial<<','<<(position==Position::Midpoint?"midpoint":"end")
                <<','<<r.transportedVertices<<','<<r.transportedEdges<<','<<r.tracedCenters<<','<<r.prescribedVertices<<','<<r.prescribedEdges<<",PASS\n";
            if(position==Position::End){const auto retained=std::move(result);reject([&]{result.publication();},"moved-from result");
                check(retained.sourceIdentity()==id,"moved step lost provenance");}
            ++cases;
        }
        reject([&]{prepareHermiteGalsFirstJetStep(*mesh,source,domain,velocity.previous(),velocity.arrival(),{},config);},"requires prescribed first jet");
        const HermiteGalsBoundaryFirstJetQuery nan=[](Vector2D,HermiteGalsTimeLayer,Boundary){return HermiteGalsBoundaryFirstJet{std::numeric_limits<double>::quiet_NaN(),{0,0}};};
        reject([&]{prepareHermiteGalsFirstJetStep(*mesh,source,domain,velocity.previous(),velocity.arrival(),nan,config);},"nonfinite prescribed boundary");
        const HermiteGalsBoundaryFirstJetQuery throwing=[](Vector2D,HermiteGalsTimeLayer,Boundary)->HermiteGalsBoundaryFirstJet{throw std::runtime_error("fixture boundary failure");};
        reject([&]{prepareHermiteGalsFirstJetStep(*mesh,source,domain,velocity.previous(),velocity.arrival(),throwing,config);},"fixture boundary failure");
        const HermiteGalsVelocityView wrongPrevious(authority.identity(),source.state().timeLayer(),[&](Vector2D q){return s.velocity(q,0);});
        reject([&]{prepareHermiteGalsFirstJetStep(*mesh,source,domain,wrongPrevious,velocity.arrival(),boundary,config);},"previous velocity geometry differs");
        auto wrongMesh=makeMesh(8,amr?0:1);
        reject([&]{prepareHermiteGalsFirstJetStep(*wrongMesh,source,domain,velocity.previous(),velocity.arrival(),boundary,config);},"mesh");
        const HermiteGalsVelocityView wrongStep(authority.identity(),{2,0,.01,Position::End},[&](Vector2D q){return s.velocity(q,.01);});
        reject([&]{prepareHermiteGalsFirstJetStep(*mesh,source,domain,velocity.previous(),wrongStep,boundary,config);},"next step");
        const HermiteGalsVelocityView wrongBegin(authority.identity(),{1,.001,.01,Position::End},[&](Vector2D q){return s.velocity(q,.01);});
        reject([&]{prepareHermiteGalsFirstJetStep(*mesh,source,domain,velocity.previous(),wrongBegin,boundary,config);},"does not begin at source time");
        auto bad=config;bad.unitNormalTolerance=0;
        reject([&]{prepareHermiteGalsFirstJetStep(*mesh,source,domain,velocity.previous(),velocity.arrival(),boundary,bad);},"invalid publication tolerance");
        check(source.state().identity()==id && &source.state().geometry().snapshot()==&scalar && source.state().primary().phaseVolumes.empty(),"failed step changed source");
    }
    auto mesh=makeMesh(8);Case s;s.pureFirstJet=true;s.motion=1;const auto source=initial(*mesh,s);
    const HermiteGalsTimeLayer time{1,0,.01,Position::End};HermiteGalsCandidate authority(*mesh,time,source.state().primary());
    const auto velocity=HermiteGalsTemporalVelocitySource::constant(source.state(),authority.identity(),time,s.velocity({},0).value);
    const auto domain=traceDomain(*mesh,false);
    const HermiteGalsBoundaryFirstJetQuery boundary=[&](Vector2D q,HermiteGalsTimeLayer t,Boundary){
        check(q.x==0 || q.x==1 || q.y==0 || q.y==1,"interior boundary reset");const auto j=s.exact(q,t.time());return HermiteGalsBoundaryFirstJet{j.v,{j.x,j.y}};};
    auto legacyCase=s;legacyCase.pureFirstJet=false;const auto legacy=initial(*mesh,legacyCase);
    reject([&]{prepareHermiteGalsFirstJetStep(*mesh,legacy,domain,velocity.previous(),velocity.arrival(),boundary);},"requires FirstJet source without inventory");
    const HermiteGalsVelocityView initialTarget(authority.identity(),{1,0,0,Position::Initial},[&](Vector2D q){return s.velocity(q,0);});
    reject([&]{prepareHermiteGalsFirstJetStep(*mesh,source,domain,velocity.previous(),initialTarget,boundary);},"target must be Midpoint/End");
    auto wrongDomain=domain;wrongDomain.maxCoord.x=2;
    reject([&]{prepareHermiteGalsFirstJetStep(*mesh,source,wrongDomain,velocity.previous(),velocity.arrival(),boundary);},"domain must match");
    wrongDomain=domain;wrongDomain.cellSize={1,1};
    reject([&]{prepareHermiteGalsFirstJetStep(*mesh,source,wrongDomain,velocity.previous(),velocity.arrival(),boundary);},"cell scale");
    const auto mid=prepareHermiteGalsFirstJetStep(*mesh,source,domain,velocity.previous(),velocity.constantArrivalView(authority.identity(),{1,0,.01,Position::Midpoint}),boundary);
    reject([&]{prepareHermiteGalsFirstJetStep(*mesh,mid.publication(),domain,velocity.previous(),velocity.arrival(),boundary);},"source must be accepted Initial/End");
    const auto tooLarge=HermiteGalsTemporalVelocitySource::constant(source.state(),authority.identity(),{1,0,1,Position::End},s.velocity({},0).value);
    reject([&]{prepareHermiteGalsFirstJetStep(*mesh,source,domain,tooLarge.previous(),tooLarge.arrival(),boundary);},"sampled speed times full dt exceeds half a cell");
    const auto large=HermiteGalsTemporalVelocitySource::constant(source.state(),authority.identity(),{1,0,.0625,Position::End},{1,0});
    reject([&]{prepareHermiteGalsFirstJetStep(*mesh,source,domain,large.previous(),large.arrival(),boundary);},"interior Jacobian unavailable: boundary event");
    // This shear vanishes on horizontal master edges. Interior vertical edges
    // remain inside, while first-column centers actually hit the inflow wall.
    auto smallMesh=makeMesh(2);const auto smallSource=initial(*smallMesh,s);
    const auto shear=[](Vector2D q){const double y=q.y,f=y*(y-.5)*(y-1),df=3*y*y-3*y+.5;
        return LevelSetBoundaryTrace::VelocityJet{{256*f*f,0},{0,512*f*df,0,0}};};
    const HermiteGalsTemporalVelocitySource centerHit(smallSource.state(),authority.identity(),{1,0,4./9,Position::End},shear,shear);
    reject([&]{prepareHermiteGalsFirstJetStep(*smallMesh,smallSource,traceDomain(*smallMesh,false),centerHit.previous(),centerHit.arrival(),boundary);},
        "center prediction requires unresolved boundary data");
    std::cout<<"I12p3a FirstJet step PASS cases="<<cases<<" guards="<<guards<<'\n';
}
void localFirstJetContract()
{
    using Rule=HermiteGalsTimeProjection;using Jacobian=LevelSetBoundaryTrace::Jacobian;
    // One requested cell on the same original quartic closed map. The n8
    // support is large enough to distinguish a real local closure from full.
    auto mesh=makeMesh(8);Case s;s.pureFirstJet=true;const auto source=initial(*mesh,s);
    const auto sourceId=source.state().identity();const auto* sourceAddress=&source.state().geometry().snapshot();
    const HermiteGalsTimeLayer time{1,0,.01,Position::Midpoint};
    HermiteGalsCandidate authority(*mesh,time,source.state().primary());
    const auto& original=*sourceAddress;
    std::vector<std::array<double,25>> saved;
    for(int c=0;c<mesh->cellIndexNum();++c)saved.push_back(original.nativeQuarticLocalDofs(c));
    check(supportsHermiteGalsFirstJetLocalSample(*mesh,source),"uniform XY original source has no declared local FirstJet support");
    int requested=-1;
    for(int c=0;c<mesh->cellIndexNum();++c){const auto q=mesh->cells(mesh->cellIndexToLabel[c]).centerCoord();
        if(q.x==.4375&&q.y==.4375){requested=c;break;}}
    check(requested>=0,"existing n8 fixture has no actual requested native cell");
    const std::vector<int> cells{requested};const auto box=original.nativeCellBounds(requested);
    const auto map=[](Vector2D p){constexpr double a=.04;const double x=p.x*(1-p.x),y=p.y*(1-p.y);
        return std::make_pair(Vector2D{p.x+a*x*x*y*y,p.y},
            Jacobian{1+2*a*x*(1-2*p.x)*y*y,2*a*x*x*y*(1-2*p.y),0,1});};
    const auto motion=[&](std::size_t& calls){return HermiteGalsResolvedMotion(sourceId,source.state().timeLayer(),authority.identity(),time,
        [&,map](Vector2D q){++calls;const auto mapped=map(q);
            return LevelSetBoundaryTrace::JetResult{{mapped.first,time.time(),{},false,false},mapped.second,1};});};
    const auto sameTime=[](HermiteGalsTimeLayer a,HermiteGalsTimeLayer b){
        return a.step==b.step&&a.begin==b.begin&&a.end==b.end&&a.position==b.position;};
    std::optional<HermiteGalsFirstJetSample> fullPullback;
    for(Rule rule:{Rule::CharacteristicPullback,Rule::ClosedProjectionCompensated}){
        std::size_t fullQueries=0,localQueries=0;const auto fullMotion=motion(fullQueries),localMotion=motion(localQueries);
        const auto full=prepareHermiteGalsFirstJetSample(*mesh,source,fullMotion,{}, {},HermiteGalsBoundaryTreatment::InflowHistory,rule);
        const auto local=prepareHermiteGalsFirstJetLocalSample(*mesh,source,localMotion,cells,{}, {},HermiteGalsBoundaryTreatment::InflowHistory,rule);
        full.validateMesh(*mesh);local.validateMesh(*mesh);const auto fullView=full.geometry(),localView=local.geometry();
        const auto& complete=fullView.snapshot();const auto& partial=localView.snapshot();const auto& r=local.report();
        check(local.identity()!=sourceId&&local.identity()!=full.identity()&&local.identity()!=authority.identity()
            &&local.sourceIdentity()==sourceId&&local.arrivalVelocityIdentity()==authority.identity()
            &&sameTime(local.timeLayer(),time)&&sameTime(localView.timeLayer(),time)&&localView.identity()==local.identity()
            &&full.sourceIdentity()==sourceId&&full.arrivalVelocityIdentity()==authority.identity()&&sameTime(full.timeLayer(),time),
            "local/full sample lost its distinct scalar owner, actual source/arrival or complete Midpoint clock");
        check(local.nativeCells()==cells&&r.requestedNativeCells==cells.size()&&r.evaluatedNativeCells==cells.size(),
            "local work inventory is not the actual requested/output native support");
        std::size_t active=0;int missing=-1;
        for(int c=0;c<mesh->cellIndexNum();++c){if(partial.nativeCellActive(c))++active;else if(missing<0)missing=c;}
        check(active==cells.size()&&missing>=0,"local scalar silently prepared a full output geometry");
        check(partial.identicalNativeQuarticCell(requested,complete,requested),"local requested exact Q4 differs from the same-time full scalar");
        const auto localCopy=local;const auto copiedView=localCopy.geometry();
        check(localCopy.identity()==local.identity()&&&copiedView.snapshot()==&partial,"local copy changed its unique immutable scalar owner");
        check(r.evaluatedPointQueries>0&&fullQueries>localQueries&&localQueries>0,
            "local kernel did not reduce actual motion point queries");
        if(rule==Rule::CharacteristicPullback){check(r.stagedNativeCells==0,"pullback claims a compensated stage");fullPullback=full;}
        else check(r.stagedNativeCells>1&&r.stagedNativeCells<=static_cast<std::size_t>(mesh->cellIndexNum()),
            "compensated sample lost the actual multi-cell departure support");
        double correction=0,fullJetCorrection=0,localJetCorrection=0;
        for(double x:{0.,.5,1.})for(double y:{0.,.5,1.}){
            const Vector2D q{box[0].x+x*(box[1].x-box[0].x),box[0].y+y*(box[1].y-box[0].y)};
            const auto first=map(q);const auto one=original.scalarJetNative(first.first);auto expected=one.value;
            auto gradient=first.second.transposeTimes({one.dx,one.dy});
            if(rule==Rule::ClosedProjectionCompensated){check(fullPullback.has_value(),"missing actual full first-departure stage");
                const auto stageView=fullPullback->geometry();const auto staged=stageView.snapshot().scalarJetNative(first.first);
                const auto second=map(first.first);const auto twice=original.scalarJetNative(second.first);
                const auto twiceGradient=second.second.transposeTimes({twice.dx,twice.dy});
                expected+=.5*(staged.value-twice.value);
                gradient=first.second.transposeTimes({one.dx+.5*(staged.dx-twiceGradient.x),one.dy+.5*(staged.dy-twiceGradient.y)});
                correction=std::hypot(correction,.5*(staged.value-twice.value));}
            const auto actual=partial.scalarJetInNativeCell(requested,q);
            if(rule==Rule::ClosedProjectionCompensated){
                const auto pullbackView=fullPullback->geometry();
                const auto pullback=pullbackView.snapshot().scalarJetInNativeCell(requested,q);
                const auto oracle=complete.scalarJetInNativeCell(requested,q);
                const std::array<double,3> fullDelta{oracle.value-pullback.value,oracle.dx-pullback.dx,oracle.dy-pullback.dy};
                const std::array<double,3> localDelta{actual.value-pullback.value,actual.dx-pullback.dx,actual.dy-pullback.dy};
                for(std::size_t k=0;k<fullDelta.size();++k){
                    fullJetCorrection=std::hypot(fullJetCorrection,fullDelta[k]);
                    localJetCorrection=std::hypot(localJetCorrection,localDelta[k]);
                    close(localDelta[k],fullDelta[k],1e-12,"local/full same-point actual FirstJet correction");
                }
            }
            close(actual.value,expected,1e-12,"local actual second-departure/primary-center scalar formula");
            if(x!=.5&&y!=.5){close(actual.dx,gradient.x,1e-12,"local actual vertex pullback x gradient");
                close(actual.dy,gradient.y,1e-12,"local actual vertex pullback y gradient");}
        }
        if(rule==Rule::ClosedProjectionCompensated){
            // The old full n4 norm includes derivative and closure channels.
            // Use the actual same-cell full phi/dx/dy correction as the oracle.
            std::cout<<"LOCAL_FIRST_JET_CORRECTION_ORACLE,"<<correction<<','<<fullJetCorrection<<','<<localJetCorrection<<'\n';
            check(fullJetCorrection>0&&localJetCorrection>0,"actual full/local FirstJet compensation is inactive");
            close(localJetCorrection,fullJetCorrection,1e-12,"local/full actual FirstJet correction norm");
        }
        const Vector2D edge{box[1].x,(box[0].y+box[1].y)*.5};const auto departure=map(edge).first;
        check(departure.x>box[1].x,"actual first departure did not exercise neighboring staged native support");
        rejects([&]{partial.scalarJetInNativeCell(missing,mesh->cells(mesh->cellIndexToLabel[missing]).centerCoord());},"inactive leaf");
        std::cout<<"LOCAL_FIRST_JET_WORK,"<<int(rule)<<','<<r.requestedNativeCells<<','<<r.evaluatedNativeCells
            <<','<<r.stagedNativeCells<<','<<r.evaluatedPointQueries<<','<<fullQueries<<','<<localQueries<<','<<correction
            <<",distinct_output_and_stage_support,"<<full.report().timings.total<<','<<local.report().timings.total
            <<",single_kernel_observation_not_solver_speedup,PASS\n";
    }
    for(int c=0;c<mesh->cellIndexNum();++c)check(saved[c]==original.nativeQuarticLocalDofs(c),"local sample changed original Q4 fields");
    check(source.state().identity()==sourceId&&&source.state().geometry().snapshot()==sourceAddress,"local preparation replaced the original scalar authority");
    std::cout<<"LOCAL_FIRST_JET_PASS,2,exact_requested_Q4,second_departure,Midpoint_fullclock,missing_support\n";
}


void compatibleIncrementContract()
{
    using Jacobian=LevelSetBoundaryTrace::Jacobian;using Update=HermiteGalsGradientUpdate;
    check(HermiteGalsPublicationConfig{}.gradientUpdate==Update::CharacteristicPullback,"gradient default changed");
    auto mesh=makeMesh(8);
    const Profile profile=[](Vector2D q){return Jet{(q.x-.5)*(q.x-.5)+(q.y-.5)*(q.y-.5)-.04+.005*std::pow(q.x,4)+.003*std::pow(q.y,4),
        2*(q.x-.5)+.02*std::pow(q.x,3),2*(q.y-.5)+.012*std::pow(q.y,3),0};};
    const auto initialState=stateFor(*mesh,profile,false,{}, {},{},HermiteGalsTransport::FirstJet);
    HermiteGalsCandidate candidate(*mesh,{},initialState.primary());const auto mixed=closeHermiteGalsMixedAmr(*mesh,candidate);
    const auto source=publishHermiteGalsFirstJetGeometry(*mesh,std::move(candidate),mixed,{candidate.identity(),centers(*mesh,profile)});
    HermiteGalsPublicationConfig config;config.gradientUpdate=Update::CompatibleIncrement;
    auto prepared=HermiteGalsFirstJetLocalPreparation::prepare(*mesh,source);int cases=0;
    for(Position position:{Position::Midpoint,Position::End})for(int mapKind:{0,1,2}){
        const HermiteGalsTimeLayer time{1,0,.01,position};HermiteGalsCandidate authority(*mesh,time,source.state().primary());
        const auto map=[&](Vector2D q){
            if(mapKind==0)return std::make_pair(q,Jacobian::identity());
            if(mapKind==1)return std::make_pair(Vector2D{.99*q.x+.005,.99*q.y+.005},Jacobian{.99,0,0,.99});
            const double a=.04*(position==Position::End?1:.5),x=q.x*(1-q.x),y=q.y*(1-q.y);
            return std::make_pair(Vector2D{q.x+a*x*x*y*y,q.y},Jacobian{1+2*a*x*(1-2*q.x)*y*y,2*a*x*x*y*(1-2*q.y),0,1});};
        const HermiteGalsResolvedMotion motion(source.state().identity(),source.state().timeLayer(),authority.identity(),time,
            [&](Vector2D q){const auto actual=map(q);return LevelSetBoundaryTrace::JetResult{{actual.first,time.time(),{},false,false},actual.second,1};});
        const auto full=prepareHermiteGalsFirstJetStep(*mesh,source,motion,{},config);
        check(full.publication().state().primary().gradientUpdate==Update::CompatibleIncrement,"actual state mode tag missing");
        if(mapKind<2){
            for(const auto& v:full.publication().state().primary().vertices){const auto b=map(point(*mesh,v.vertex));const auto j=profile(b.first);
                const auto g=b.second.transposeTimes({j.dx,j.dy});close(v.phi,j.value,1e-12,"quartic value reproduction");
                close(v.dx,g.x,1e-11,"quartic gradient x / one-sided boundary");close(v.dy,g.y,1e-11,"quartic gradient y / one-sided boundary");}
            for(const auto& e:full.publication().state().primary().edges){const auto a=point(*mesh,e.vertices[0]),b=point(*mesh,e.vertices[1]);
                const auto mapped=map(a+(b-a)*.5);const auto j=profile(mapped.first);const auto g=mapped.second.transposeTimes({j.dx,j.dy});
                close(e.dn,a.y==b.y?g.y:g.x,1e-11,"quartic edge derivative");}
        }
        auto view=full.publication().state().geometry();const auto& scalar=view.snapshot();
        for(int cell:{0,27,63})for(auto execution:{HermiteGalsPointExecution::Serial,HermiteGalsPointExecution::Parallel}){
            const auto local=prepared.sample(motion,{cell},config,HermiteGalsTimeProjection::CharacteristicPullback,execution);
            check(local.geometry().snapshot().identicalNativeQuarticCell(cell,scalar,cell),"compatible increment full/local mismatch");
            const auto box=scalar.nativeCellBounds(cell);using Trace=HermiteGalsFirstJetTraceSample;
            const std::vector<Trace::Segment> segments{{box[0],{box[1].x,box[0].y},1},
                {{box[0].x,(box[0].y+box[1].y)*.5},{box[1].x,(box[0].y+box[1].y)*.5},1}};
            std::set<int> traceSupport{cell};
            for(const auto& segment:segments)for(int c:scalar.nativeCellCandidates(segment.begin,segment.end))traceSupport.insert(c);
            const auto packet=prepared.sampleTraces(motion,std::vector<int>(traceSupport.begin(),traceSupport.end()),segments,config,HermiteGalsTimeProjection::CharacteristicPullback,execution);
            const auto owner=std::make_shared<const HermiteGalsGeometryView>(view);
            const std::shared_ptr<const Snapshot> retained(owner,&owner->snapshot());
            for(std::size_t k=0;k<segments.size();++k){const auto& segment=segments[k];const auto ref=LevelSetMappedSegment::nativeQuarticTrace(retained,segment.begin,segment.end,segment.normalAxis);
                std::visit([&](const auto& trace){std::array<long double,2> a{},b{};
                    for(const auto& range:trace.ordinaryRanges())a[range.negative?0:1]+=range.length;
                    for(const auto& range:ref.ordinaryRanges())b[range.negative?0:1]+=range.length;
                    for(int phase=0;phase<2;++phase)close(a[phase],b[phase],1e-12,"compatible cardinal/full phase lengths");},packet.trace(k));}
        }
        if(mapKind==0){Case shape;shape.pureFirstJet=true;const auto sourceBefore=source.state().identity();
            rejects([&]{prepareHermiteGalsFirstJetStep(*mesh,source,motion,{},config,HermiteGalsBoundaryTreatment::InflowHistory,
                HermiteGalsTimeProjection::ClosedProjectionCompensated);},"requires characteristic temporal projection");
            auto invalid=config;invalid.scalarRepresentation=HermiteGalsScalarRepresentation::C2MacroQuartic;
            rejects([&]{prepareHermiteGalsFirstJetStep(*mesh,source,motion,{},invalid);},"requires SingleQuartic source");
            check(source.state().identity()==sourceBefore,"rejected update mutated original");}
        ++cases;
    }
    std::cout<<"COMPATIBLE_INCREMENT_PASS,cases="<<cases<<",Q4_zero_affine,one_sided,Midpoint_End,full_local_cardinal,serial_parallel\n";
}

void closedProjectionContract()
{
    using Rule=HermiteGalsTimeProjection;using Jacobian=LevelSetBoundaryTrace::Jacobian;
    auto mesh=makeMesh(4);Case s;s.pureFirstJet=true;const auto source=initial(*mesh,s);
    const HermiteGalsTimeLayer time{1,0,.01,Position::End};HermiteGalsCandidate authority(*mesh,time,source.state().primary());
    const auto fields=[](const HermiteGalsPublication& pub){std::vector<double> v;const auto& p=pub.state().primary();
        for(const auto& a:p.vertices){v.push_back(a.phi);v.push_back(a.dx);v.push_back(a.dy);}
        for(const auto& a:p.edges){v.push_back(a.phi);v.push_back(a.dn);}
        const auto& c=pub.state().closure();v.insert(v.end(),c.vertexMixedDerivatives.begin(),c.vertexMixedDerivatives.end());
        v.insert(v.end(),c.cellBubbleCoefficients.begin(),c.cellBubbleCoefficients.end());return v;};
    const auto saved=fields(source);const auto id=source.state().identity();const auto* scalarAddress=&source.state().geometry().snapshot();
    const auto motion=[&](auto map){return HermiteGalsResolvedMotion(id,source.state().timeLayer(),authority.identity(),time,
        [map](Vector2D p){const auto a=map(p);return LevelSetBoundaryTrace::JetResult{{a.first,.01,{},false,false},a.second,1};});};
    for(int kind:{0,1}){
        const auto actual=motion([kind](Vector2D p){return std::make_pair(kind==0?p:Vector2D{.99*p.x+.005,.99*p.y+.005},
            kind==0?Jacobian::identity():Jacobian{.99,0,0,.99});});
        const auto old=prepareHermiteGalsFirstJetStep(*mesh,source,actual);
        const auto next=prepareHermiteGalsFirstJetStep(*mesh,source,actual,{}, {},HermiteGalsBoundaryTreatment::InflowHistory,Rule::ClosedProjectionCompensated);
        const auto a=fields(old.publication()),b=fields(next.publication());double gap=0,na=0,nb=0;
        for(std::size_t i=0;i<a.size();++i){gap=std::hypot(gap,a[i]-b[i]);na=std::hypot(na,a[i]);nb=std::hypot(nb,b[i]);}
        const double scale=std::max(na,nb);check(scale==0?gap==0:gap/scale<=1e-12,"zero/affine projection changed exact representable field");
        check(next.sourceIdentity()==id&&next.arrivalVelocityIdentity()==authority.identity()&&next.publication().state().timeLayer().time()==time.time(),
            "compensated map lost original provenance/clock");
        std::cout<<"CLOSED_PROJECTION_EXACT,"<<kind<<','<<gap<<','<<scale<<"\n";
    }
    // A quadratic scalar pulled through the previous biquadratic map still
    // belongs to Q4, where compensation correctly vanishes. This closed
    // quartic map produces degree-eight terms outside the reconstruction.
    const auto map=[](Vector2D p){constexpr double a=.04;const double x=p.x*(1-p.x),y=p.y*(1-p.y);
        return std::make_pair(Vector2D{p.x+a*x*x*y*y,p.y},
            Jacobian{1+2*a*x*(1-2*p.x)*y*y,2*a*x*x*y*(1-2*p.y),0,1});};
    const auto actual=motion(map);const auto old=prepareHermiteGalsFirstJetStep(*mesh,source,actual);
    const auto next=prepareHermiteGalsFirstJetStep(*mesh,source,actual,{}, {},HermiteGalsBoundaryTreatment::InflowHistory,Rule::ClosedProjectionCompensated);
    const auto a=fields(old.publication()),b=fields(next.publication());double correction=0;
    for(std::size_t i=0;i<a.size();++i)correction=std::hypot(correction,a[i]-b[i]);check(correction>1e-10,"nonaffine projection compensation is inactive");
    const auto& original=source.state().geometry().snapshot();const auto& projected=old.publication().state().geometry().snapshot();
    const auto value=[&](Vector2D q){const auto d=map(q).first,d2=map(d).first;
        return original.scalarJetNative(d).value+.5*(projected.scalarJetNative(d).value-original.scalarJetNative(d2).value);};
    double fd=0;constexpr double delta=1e-5;
    for(const auto& v:next.publication().state().primary().vertices){const auto q=point(*mesh,v.vertex);
        close(v.phi,value(q),1e-12,"compensated primary pullback value");
        if(q.x==0||q.x==1||q.y==0||q.y==1)continue;
        const double gx=(value(q+Vector2D{delta,0})-value(q-Vector2D{delta,0}))/(2*delta);
        const double gy=(value(q+Vector2D{0,delta})-value(q-Vector2D{0,delta}))/(2*delta);
        fd=std::max({fd,std::abs(v.dx-gx),std::abs(v.dy-gy)});
    }
    check(fd<=5e-7,"compensated Jacobian differs from independent scalar finite difference");
    const HermiteGalsResolvedMotion secondHit(id,source.state().timeLayer(),authority.identity(),time,[](Vector2D p){
        if(p.x>0&&p.x<1&&std::abs(8*p.x-std::round(8*p.x))>1e-10)
            return LevelSetBoundaryTrace::JetResult{{{0,p.y},.005,Boundary::west,true,false},{},1};
        return LevelSetBoundaryTrace::JetResult{{{p.x==0||p.x==1?p.x:p.x-.01,p.y},.01,{},false,false},Jacobian::identity(),1};});
    rejects([&]{prepareHermiteGalsFirstJetStep(*mesh,source,secondHit,{}, {},HermiteGalsBoundaryTreatment::InflowHistory,Rule::ClosedProjectionCompensated);},
        "closed projection second query requires interior Jacobian");
    const HermiteGalsBoundaryFirstJetQuery supplied=[](Vector2D,HermiteGalsTimeLayer,Boundary){return HermiteGalsBoundaryFirstJet{};};
    rejects([&]{prepareHermiteGalsFirstJetStep(*mesh,source,actual,supplied,{},HermiteGalsBoundaryTreatment::InflowHistory,Rule::ClosedProjectionCompensated);},
        "closed projection requires no prescribed boundary");
    rejects([&]{prepareHermiteGalsFirstJetStep(*mesh,source,actual,{}, {},HermiteGalsBoundaryTreatment::PrescribedOuterFirstJet,Rule::ClosedProjectionCompensated);},
        "closed projection requires no prescribed boundary");
    rejects([&]{prepareHermiteGalsFirstJetStep(*mesh,source,actual,{}, {},HermiteGalsBoundaryTreatment::InflowHistory,static_cast<Rule>(-1));},
        "unknown temporal projection rule");
    check(source.state().identity()==id && &source.state().geometry().snapshot()==scalarAddress&&fields(source)==saved,
        "failed compensated step changed original publication");
    std::cout<<"CLOSED_PROJECTION_PASS,"<<correction<<','<<fd<<",guards=4\n";
    localFirstJetContract();
}
void firstJetMeshTransferContract()
{
    using Kind=HermiteGalsMeshTransferKind;
    int exactCases=0,projectionCases=0,guards=0;
    const auto reject=[&](auto f,const std::string& why){rejects(f,why);++guards;};
    for(bool radial:{false,true})for(bool amr:{false,true})for(bool nonpolynomial:{false,true}){
        auto sourceMesh=makeMesh(8,amr?1:0);Case s;s.pureFirstJet=true;s.radial=radial;s.amr=amr;s.nonpolynomial=nonpolynomial;if(radial)s.radius=.27;
        const auto source=initial(*sourceMesh,s);const auto id=source.state().identity();const auto& scalar=source.state().geometry().snapshot();
        std::optional<HermiteGalsPublication> refined;std::unique_ptr<QuadtreeMesh> refinedMesh;
        for(bool identical:{true,false}){
            auto target=makeMesh(8,identical?(amr?1:0):(amr?2:1));
            const auto result=transferHermiteGalsFirstJetMesh(*sourceMesh,source,*target,Kind::ExactRefinement);const auto& pub=result.publication();
            const auto& actual=pub.state().geometry().snapshot();const auto time=pub.state().timeLayer();
            check(result.sourceIdentity()==id && pub.state().identity()!=id && time.time()==0 && time.step==0
                && pub.state().primary().phaseVolumes.empty(),"exact regrid state/time authority");
            check(&pub.volumes().geometry().snapshot()==&actual && &pub.facePartitions().source()==&actual,"regrid has multiple scalar owners");
            double jet=0,trace=0,curvature=0;std::vector<long double> childVolumes(sourceMesh->cellIndexNum(),0);
            for(int c=0;c<target->cellIndexNum();++c){const auto b=actual.nativeCellBounds(c);const auto h=b[1]-b[0];const auto center=b[0]+h*.5;int parent=-1;
                for(int p:scalar.nativeCellCandidates(center,center)){const auto a=scalar.nativeCellBounds(p);
                    if(center.x>a[0].x && center.x<a[1].x && center.y>a[0].y && center.y<a[1].y){check(parent<0,"ambiguous exact parent");parent=p;}}
                check(parent>=0,"missing exact parent");childVolumes[parent]+=pub.volumes().phaseVolumes()[c].liquid;
                for(double x:{.25,.5,.75})for(double y:{.25,.5,.75}){const auto p=b[0]+Vector2D{x*h.x,y*h.y};
                    const auto a=scalar.scalarJetNative(p),z=actual.scalarJetAtNativeCellChart(c,{{x,false},{y,false}});
                    const double av[]={a.value,a.dx,a.dy,a.dxx,a.dxy,a.dyy},zv[]={z.value,z.dx,z.dy,z.dxx,z.dxy,z.dyy};
                    for(int k=0;k<6;++k)jet=std::max(jet,std::abs(zv[k]-av[k])/(1+std::abs(av[k])));}
                for(double t:{0.,.5,1.})for(Vector2D local:std::array<Vector2D,4>{{{0,t},{1,t},{t,0},{t,1}}}){
                    const auto a=scalar.scalarJetNative(b[0]+Vector2D{local.x*h.x,local.y*h.y}),z=actual.scalarJetAtNativeCellChart(c,
                        {{local.x,false},{local.y,false}});trace=std::max({trace,std::abs(a.value-z.value),std::abs(a.dx-z.dx),std::abs(a.dy-z.dy)});}
            }
            long double volume=0;for(const auto& c:source.volumes().cells())volume=std::max(volume,
                std::abs(childVolumes[c.cell]-source.volumes().phaseVolumes()[c.cell].liquid)/c.capacity);
            for(std::size_t c=0;c<pub.topology().curves.size();++c)for(std::size_t g=0;g<pub.topology().curves[c].graphs().size();++g)
                for(const auto& q:pub.quadrature(c,g)){const auto p=q.sample.frame.point;const auto j=scalar.scalarJetNative(p);
                    double reference=curvatureOf(j);if(radial)reference+=j.dx/(std::hypot(j.dx,j.dy)*p.x);
                    curvature=std::max(curvature,std::abs(q.sample.curvature-reference));}
            check(jet<=1e-10 && trace<=1e-12 && curvature<=1e-9 && volume<=1e-10,"exact regrid scalar/geometry gate");
            const auto& r=result.report();check(r.coarsenedCells==0 && r.refinedCells+r.unchangedCells==std::size_t(target->cellIndexNum()),"exact regrid leaf report");
            check((identical && r.refinedCells==0) || (!identical && r.refinedCells>0),"exact regrid fixture did not refine");
            std::cout<<"FIRST_JET_EXACT_REGRID,"<<radial<<','<<amr<<','<<nonpolynomial<<','<<identical<<','<<jet<<','<<trace<<','<<curvature<<','<<volume<<",PASS\n";
            ++exactCases;
            if(!identical){refined=pub;refinedMesh=std::move(target);}
        }
        reject([&]{transferHermiteGalsFirstJetMesh(*refinedMesh,*refined,*sourceMesh,Kind::ExactRefinement);},"ExactRefinement cannot coarsen");
        const auto projected=transferHermiteGalsFirstJetMesh(*refinedMesh,*refined,*sourceMesh,Kind::FirstJetProjection);
        const auto& pub=projected.publication();const auto e=errors(pub,s);
        check(projected.report().coarsenedCells>0 && pub.state().timeLayer().time()==0 && pub.state().primary().phaseVolumes.empty(),"projection did not coarsen with zero time");
        if(!nonpolynomial)check(std::max({e.position,e.gradient,e.curvature})<=1e-10 && e.cellVolume<=5e-12,"Q2 coarsening reproduction");
        std::cout<<"FIRST_JET_PROJECT_REGRID,"<<radial<<','<<amr<<','<<nonpolynomial<<','<<e.position<<','<<e.gradient<<','<<e.curvature
            <<','<<e.volume<<','<<projected.report().relativeDomainVolumeChange<<",PASS\n";++projectionCases;
        reject([&]{transferHermiteGalsFirstJetMesh(*sourceMesh,source,*sourceMesh,static_cast<Kind>(-1));},"unknown transfer kind");
        reject([&]{transferHermiteGalsFirstJetMesh(*refinedMesh,source,*sourceMesh,Kind::FirstJetProjection);},"mesh");
        HermiteGalsPublicationConfig bad;bad.unitNormalTolerance=0;
        reject([&]{transferHermiteGalsFirstJetMesh(*sourceMesh,source,*sourceMesh,Kind::ExactRefinement,bad);},"invalid publication tolerance");
        auto moved=projected;const auto retained=std::move(moved);
        reject([&]{moved.publication();},"moved-from result");
        check(retained.publication().state().identity()==pub.state().identity() && source.state().identity()==id
            && &source.state().geometry().snapshot()==&scalar,"regrid failure mutated source/owner");
    }
    auto mesh=makeMesh(8);Case s;s.pureFirstJet=true;const auto source=initial(*mesh,s);
    auto nonnested=makeMesh(7);reject([&]{transferHermiteGalsFirstJetMesh(*mesh,source,*nonnested,Kind::FirstJetProjection);},"nested source/target leaves");
    DictionaryMesh d;d.Nx=d.Ny=8;d.maxInitialIterations=0;d.maxCoord={2,1};QuadtreeMesh foreign(d);
    reject([&]{transferHermiteGalsFirstJetMesh(*mesh,source,foreign,Kind::ExactRefinement);},"changes domain");
    auto legacyCase=s;legacyCase.pureFirstJet=false;const auto legacy=initial(*mesh,legacyCase);
    reject([&]{transferHermiteGalsFirstJetMesh(*mesh,legacy,*mesh,Kind::FirstJetProjection);},"requires FirstJet source without inventory");
    // A partial planar publication needs its actual outer boundary contacts.
    auto partial=source.state().primary();partial.geometryCells.pop_back();std::set<int> vertices;std::set<Key> sides;
    for(int c:partial.geometryCells){const auto& cell=mesh->cells(mesh->cellIndexToLabel[c]);for(int k=0;k<4;++k){
        const int v=mesh->verts(cell.cellToVert[k]).index();vertices.insert(v);Key key{v,mesh->verts(cell.cellToVert[(k+1)%4]).index()};
        const auto a=point(*mesh,key[0]),b=point(*mesh,key[1]);if(a.y==b.y?a.x>b.x:a.y>b.y)std::swap(key[0],key[1]);sides.insert(key);}}
    partial.vertices.erase(std::remove_if(partial.vertices.begin(),partial.vertices.end(),[&](const auto& v){return !vertices.count(v.vertex);}),partial.vertices.end());
    partial.edges.erase(std::remove_if(partial.edges.begin(),partial.edges.end(),[&](const auto& e){return !sides.count(e.vertices);}),partial.edges.end());
    for(auto& v:partial.vertices){v.phi=point(*mesh,v.vertex).x-.53;v.dx=1;v.dy=0;}
    for(auto& e:partial.edges){const auto a=point(*mesh,e.vertices[0]),b=point(*mesh,e.vertices[1]);e.phi=(a.x+b.x)*.5-.53;e.dn=a.y==b.y?0:1;}
    HermiteGalsCandidate candidate(*mesh,{},partial);const auto mixed=closeHermiteGalsMixedAmr(*mesh,candidate);
    std::vector<double> center;for(int c:partial.geometryCells)center.push_back(mesh->cells(mesh->cellIndexToLabel[c]).centerCoord().x-.53);
    HermiteGalsPublicationConfig partialConfig;partialConfig.volume.boundaryContacts={{.53,0},{.53,1}};
    const auto local=publishHermiteGalsFirstJetGeometry(*mesh,std::move(candidate),mixed,{candidate.identity(),center},partialConfig);
    reject([&]{transferHermiteGalsFirstJetMesh(*mesh,local,*mesh,Kind::ExactRefinement);},"requires full scalar support");
    // Coarsening must not silently lose a resolved compact drop.
    auto fine=makeMesh(32),coarse=makeMesh(4);
    const Profile compact=[](Vector2D p){const D2 x{p.x-.4375,1,0,0,0,0},y{p.y-.4375,0,1,0,0,0};
        const auto j=1-2*exponential((-1/(.04*.04))*(x*x+y*y));return Jet{j.v,j.x,j.y,j.xy};};
    const auto seed=stateFor(*fine,compact,false,{},[](int,Vector2D,Vector2D){return 0.;},{},HermiteGalsTransport::FirstJet);
    HermiteGalsCandidate small(*fine,{},seed.primary());const auto smallMixed=closeHermiteGalsMixedAmr(*fine,small);
    const auto smallDrop=publishHermiteGalsFirstJetGeometry(*fine,std::move(small),smallMixed,{small.identity(),centers(*fine,compact)});
    reject([&]{transferHermiteGalsFirstJetMesh(*fine,smallDrop,*coarse,Kind::FirstJetProjection);},"changed components or endpoint kinds");
    std::cout<<"I12p2b FirstJet mesh transfer PASS exact="<<exactCases<<" projection="<<projectionCases<<" guards="<<guards<<'\n';
}
void firstJetCrossing()
{
    int cases=0;long double maximumVolume=0;double maximumJet=0,maximumCurvature=0;
    for(int configuration=0;configuration<3;++configuration)for(bool amr:{false,true})
        for(int direction:{-1,1})for(int sign:{-1,1})for(double fraction:{.25,.75}){
        Case s;s.pureFirstJet=true;s.planar=true;s.radial=configuration==2;s.axialTranslation=configuration!=0;
        s.amr=amr;s.motion=1;s.speed=direction;s.planeSign=sign;const double dt=.02;s.planeOffset=.5-s.speed*fraction*dt;
        std::string stage="initialize";
        try{
            auto mesh=makeMesh(4,amr?1:0);const auto old=initial(*mesh,s);const auto id=old.state().identity();
            const HermiteGalsTimeLayer endTime{1,0,dt,Position::End};
            HermiteGalsCandidate velocityGeometry(*mesh,endTime,old.state().primary());
            const auto velocity=HermiteGalsTemporalVelocitySource::constant(old.state(),velocityGeometry.identity(),endTime,s.velocity({},0).value);
            const HermiteGalsTemporalCharacteristicSource characteristic(*mesh,traceDomain(*mesh,s.radial));
            const auto mid=advanceFirstJet(*mesh,old,s,{1,0,dt,Position::Midpoint},stage,&velocity,&characteristic),
                end=advanceFirstJet(*mesh,old,s,endTime,stage,&velocity,&characteristic);
            check(old.state().identity()==id && old.state().primary().phaseVolumes.empty(),"pure plane mutated old state");
            long double volumeError=0;double jetError=0,curvature=0;
            for(const auto* pub:{&mid,&end}){const auto t=pub->state().timeLayer().time();const double threshold=s.planeOffset+s.speed*t;
                for(const auto& m:pub->volumes().phaseVolumes()){
                    const auto b=pub->state().geometry().snapshot().nativeCellBounds(m.cell);const double low=s.axial()?b[0].y:b[0].x,high=s.axial()?b[1].y:b[1].x;
                    const long double length=std::max(0.,std::min(threshold,high)-low),cap=pub->domain().cellVolumes()[m.cell];
                    const long double measure=s.radial?pi*(static_cast<long double>(b[1].x)-b[0].x)*(static_cast<long double>(b[1].x)+b[0].x)
                        :s.axial()?b[1].x-b[0].x:b[1].y-b[0].y;
                    const auto reference=sign>0?length*measure:cap-length*measure;
                    volumeError=std::max(volumeError,std::abs(m.liquid-reference)/cap);
                }
                for(const auto& v:pub->state().primary().vertices){const auto j=s.exact(point(*mesh,v.vertex),t);
                    jetError=std::max({jetError,std::abs(v.phi-j.v),std::abs(v.dx-j.x),std::abs(v.dy-j.y)});}
                for(const auto& e:pub->state().primary().edges){const auto a=point(*mesh,e.vertices[0]),b=point(*mesh,e.vertices[1]);const auto j=s.exact(a+(b-a)*.5,t);
                    jetError=std::max({jetError,std::abs(e.phi-j.v),std::abs(e.dn-(a.y==b.y?j.y:j.x))});}
                for(std::size_t c=0;c<pub->topology().curves.size();++c)for(std::size_t g=0;g<pub->topology().curves[c].graphs().size();++g)
                    for(const auto& q:pub->quadrature(c,g)){const auto normal=q.sample.frame.normal;
                        jetError=std::max({jetError,std::abs(normal.x-(s.axial()?0:sign)),std::abs(normal.y-(s.axial()?sign:0))});
                        curvature=std::max(curvature,std::abs(q.sample.curvature));}
            }
            check(mid.state().timeLayer().time()==dt/2 && end.state().timeLayer().time()==dt,"pure plane stage times");
            check(volumeError<=5e-12 && jetError<=1e-11 && curvature<=1e-10,"pure plane analytic gate");
            maximumVolume=std::max(maximumVolume,volumeError);maximumJet=std::max(maximumJet,jetError);maximumCurvature=std::max(maximumCurvature,curvature);++cases;
            std::cout<<"FIRST_JET_CROSSING,"<<configuration<<','<<amr<<','<<direction<<','<<sign<<','<<fraction
                <<','<<volumeError<<','<<jetError<<','<<curvature<<",PASS\n";
        }catch(const std::exception& e){++failures;std::cout<<"FIRST_JET_CROSSING_FAILURE,"<<configuration<<','<<amr<<','<<direction<<','<<sign<<','<<fraction<<','<<stage<<','<<e.what()<<'\n';}
    }
    std::cout<<"I12p2a FirstJet crossing "<<(failures?"FAIL":"PASS")<<" cases="<<cases<<" planned=48 failures="<<failures
        <<" max_volume="<<maximumVolume<<" max_jet="<<maximumJet<<" max_curvature="<<maximumCurvature<<'\n';
}
}
int main(int argc,char** argv) try {
    if(argc!=2) throw std::invalid_argument("publication | step | closed-projection | mesh-transfer | crossing");
    const std::string mode=argv[1];
    if(mode=="publication") {firstJetPublicationContract();c2MacroLocalContract();c2MacroRootContract();}
    else if(mode=="c2-macro") c2MacroLocalContract();
    else if(mode=="c2-roots") c2MacroRootContract();
    else if(mode=="step") {firstJetStepContract();closedProjectionContract();compatibleIncrementContract();}
    else if(mode=="closed-projection") closedProjectionContract();
    else if(mode=="compatible-increment") compatibleIncrementContract();
    else if(mode=="mesh-transfer") firstJetMeshTransferContract();
    else if(mode=="crossing") firstJetCrossing();
    else throw std::invalid_argument("unknown contract");
    return failures?1:0;
} catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 2;}
