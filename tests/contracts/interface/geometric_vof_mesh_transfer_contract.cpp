#include "Modules/Interface/GeometricVofMeshTransfer.h" // Public API stands alone.
#include "Core/Field/CellField.h"
#include "Core/Field/VertField.h"
#include "Core/Field/HEdgeField.h"
#include "Core/Field/VEdgeField.h"
#include "Core/Mesh/CartesianMeasure.h"
#include "Core/Mesh/FaceSegmentGraph.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/NativeMeshIdentity.h"
#include "Core/Mesh/DyadicGridCoordinates.h"
#include "Core/Geometry/Detail/FiniteVolumeClosure.h"
#include <algorithm>
#include <array>
#include <cstdint>
#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <tuple>
#include <vector>

namespace
{
using namespace LAMP::Interface;
using Publication=GeometricVofPublication;
using Transfer=GeometricVofMeshTransfer;
using Kind=GeometricVofTransferOverlap::Kind;
using Point=LAMP::Geometry::PlanarCellCut::Point;
using Clock=GeometricVofTimeLayer;
using Position=Clock::Position;
const Clock accepted{7,.125,.25,Position::End};
constexpr long double tolerance=1e-12L;
constexpr long double arithmetic=2048*std::numeric_limits<long double>::epsilon();
const long double twoPi=CartesianMeasure::twoPi;
unsigned assertions=0,scenarios=0,guards=0,overlapChecks=0,adaptations=0;
long double peakInventory=0,peakGas=0,peakRaw=0,peakGeometry=0,peakGreen=0,peakDonorClosure=0,peakCellClosure=0,peakMinorFactor=0;
std::string stage="setup";
static_assert(!std::is_default_constructible_v<Transfer>);
static_assert(std::is_same_v<decltype(std::declval<const Transfer&>().candidate()),const Publication&>);
static_assert(std::is_same_v<decltype(std::declval<const Transfer&>().overlaps()),const std::vector<GeometricVofTransferOverlap>&>);
void check(bool condition,const std::string& name)
{ ++assertions;if(!condition)throw std::runtime_error(stage+": "+name); }
void near(long double actual,long double expected,long double scale,const std::string& name,long double tol=arithmetic)
{
    check(std::isfinite(actual)&&std::isfinite(expected)&&std::isfinite(scale)&&scale>0,name+" invalid finite comparison");
    check(std::abs(actual-expected)<=scale*tol,name+" residual="+std::to_string(static_cast<double>((actual-expected)/scale)));
}
template<class F> void rejects(F&& f,const std::string& name)
{
    bool thrown=false;try{f();}catch(const std::exception&){thrown=true;}
    check(thrown,"missing explicit rejection: "+name);++guards;
}
struct Box {long double x0,x1,y0,y1;};
Box bounds(const QuadtreeMesh& mesh,int index)
{
    const auto& cell=mesh.cells(mesh.cellIndexToLabel[index]);
    const auto a=mesh.verts(cell.cellToVert[0]).coord(),b=mesh.verts(cell.cellToVert[2]).coord();
    return {a.x,b.x,a.y,b.y};
}
Box intersection(Box a,Box b)
{ return {std::max(a.x0,b.x0),std::min(a.x1,b.x1),std::max(a.y0,b.y0),std::min(a.y1,b.y1)}; }
long double capacity(Box b,bool rz)
{
    if(b.x1<=b.x0||b.y1<=b.y0)return 0;
    return (b.x1-b.x0)*(b.y1-b.y0)*(rz?twoPi*(b.x0+b.x1)/2:1);
}
NativeLeafDomain domain(const QuadtreeMesh& mesh,bool rz)
{ FaceSegmentGraph graph;graph.build(mesh);return NativeLeafDomain::build(mesh,graph,rz); }
std::unique_ptr<QuadtreeMesh> meshFor(int base=8,int level=0,bool local=false)
{
    DictionaryMesh d;d.Nx=d.Ny=base;d.maxCellNum=64*base*base;d.maxInitialIterations=0;d.maxLevel_local={2};
    auto mesh=std::make_unique<QuadtreeMesh>(d);
    if(level){DictionaryMesh::RefineBox r;r.level=level;r.min={local?.5:0,0};r.max={1,1};mesh->applyRefineBoxes({r});}
    return mesh;
}
bool sameClock(Clock a,Clock b)
{ return a.step==b.step&&a.begin==b.begin&&a.end==b.end&&a.position==b.position; }
Publication publish(const QuadtreeMesh& mesh,bool rz,const std::vector<double>& alpha,Clock clock=accepted)
{ return Publication::prepare(mesh,domain(mesh,rz),clock,alpha); }
// Independent strip antiderivatives for x+y<intercept. They provide ONLY
// cell-average input; no analytic normal or reconstructed geometry enters.
long double lineVolume(Box b,bool rz,long double intercept=.90625L)
{
    // Analytic corner classification makes full/outside states exact,
    // avoiding differently associated RZ moment arithmetic for pure cells.
    if(intercept>=b.x1+b.y1)return capacity(b,rz);
    if(intercept<=b.x0+b.y0)return 0;
    const auto full=std::clamp(intercept-b.y1,b.x0,b.x1),end=std::clamp(intercept-b.y0,b.x0,b.x1);
    long double area=(full-b.x0)*(b.y1-b.y0),radial=(full*full-b.x0*b.x0)*(b.y1-b.y0)/2;
    if(end>full){const auto h=intercept-b.y0,second=(end*end-full*full)/2,third=(end*end*end-full*full*full)/3;
        area+=h*(end-full)-second;radial+=h*second-third;}
    return rz?twoPi*radial:area;
}
using Integral=std::function<long double(Box,bool)>;
std::vector<double> input(const QuadtreeMesh& mesh,bool rz,const Integral& integral)
{
    std::vector<double> a(mesh.cellIndexNum());
    for(int i=0;i<mesh.cellIndexNum();++i){const auto b=bounds(mesh,i);const auto c=capacity(b,rz),v=integral(b,rz);
        std::ostringstream detail;detail<<std::setprecision(17)<<"analytic input exceeds true cell capacity i="<<i
            <<" box="<<b.x0<<","<<b.x1<<","<<b.y0<<","<<b.y1<<" v="<<v<<" capacity="<<c<<" normalized="<<(v-c)/c;
        check(v>=0&&v<=c,detail.str());
        a[i]=v==0?0:(v==c?1:static_cast<double>(v/c));check(a[i]>=0&&a[i]<=1,"analytic conversion outside bounds");}
    return a;
}
std::vector<double> lineInput(const QuadtreeMesh& mesh,bool rz)
{ return input(mesh,rz,[](Box b,bool radial){return lineVolume(b,radial);}); }
// Test-owned rectangle clipping and area/first-moment integration. Neither
// Core clipping nor Core measures are used to verify the transfer records.
std::vector<Point> clipSide(const std::vector<Point>& polygon,int axis,long double edge,bool greater)
{
    std::vector<Point> out;if(polygon.empty())return out;
    const auto coord=[axis](Point p){return axis==0?p.x:p.y;};
    const auto inside=[&](Point p){return greater?coord(p)>=edge:coord(p)<=edge;};
    auto a=polygon.back();bool ia=inside(a);
    for(auto b:polygon){const bool ib=inside(b);
        if(ia!=ib){const auto u=(edge-coord(a))/(coord(b)-coord(a));
            Point p{a.x+u*(b.x-a.x),a.y+u*(b.y-a.y)};if(axis==0)p.x=edge;else p.y=edge;out.push_back(p);}
        if(ib)out.push_back(b);a=b;ia=ib;}
    return out;
}
std::vector<Point> clip(const std::vector<Point>& polygon,Box b,Point origin)
{
    auto p=clipSide(polygon,0,b.x0-origin.x,true);p=clipSide(p,0,b.x1-origin.x,false);
    p=clipSide(p,1,b.y0-origin.y,true);return clipSide(p,1,b.y1-origin.y,false);
}
long double measure(const std::vector<Point>& polygon,Point origin,bool rz)
{
    long double twiceArea=0,sixMoment=0;if(polygon.size()<3)return 0;
    auto a=polygon.back();for(auto b:polygon){const auto c=a.x*b.y-b.x*a.y;
        twiceArea+=c;sixMoment+=(a.x+b.x)*c;a=b;}
    if(twiceArea<0){twiceArea=-twiceArea;sixMoment=-sixMoment;}
    return rz?twoPi*(origin.x*twiceArea/2+sixMoment/6):twiceArea/2;
}
long double phaseOverlap(const Publication& p,int sourceCell,Box target,bool liquid)
{
    long double total=0;for(const auto& polygon:p.phasePolygons(sourceCell,liquid))
        total+=measure(clip(polygon,target,p.frameOrigin()),p.frameOrigin(),p.domain().axisymmetric());
    return total;
}
Point global(Point p,Point origin){return {p.x+origin.x,p.y+origin.y};}
Point interpolate(Point a,Point b,long double t){return {a.x+t*(b.x-a.x),a.y+t*(b.y-a.y)};}
long double green(Point a,Point b,bool rz)
{
    return rz?twoPi*(b.y-a.y)*(a.x*a.x+a.x*b.x+b.x*b.x)/6:(a.x*b.y-b.x*a.y)/2;
}
void verifyGeometry(const QuadtreeMesh& mesh,const Publication& p)
{
    p.validateFor(mesh,p.identity(),accepted);const bool rz=p.domain().axisymmetric();const auto origin=p.frameOrigin();
    long double liquid=0,inventory=0,full=0,greenTotal=0;
    for(const auto& c:p.cells()){
        const auto cap=capacity(bounds(mesh,c.nativeCell),rz);long double l=0,g=0;
        for(const auto& polygon:p.phasePolygons(c.nativeCell,true))l+=measure(polygon,origin,rz);
        for(const auto& polygon:p.phasePolygons(c.nativeCell,false))g+=measure(polygon,origin,rz);
        near(c.capacity,cap,cap,"independent cell capacity");near(l+g,cap,cap,"independent complementary polygons");
        near(c.measured,l,cap,"actual geometric readback");near(l,c.inventory,cap,"independent inventory/geometry",tolerance);
        if(p.alpha()[c.nativeCell]>0&&p.alpha()[c.nativeCell]<1)check(l>0&&g>0,"mixed leaf lost an actual phase");
        peakGeometry=std::max(peakGeometry,std::abs(l-c.inventory)/cap);liquid+=l;inventory+=c.inventory;full+=cap;
    }
    for(const auto& f:p.facets())greenTotal+=green(global(p.vertices()[f.vertices[0]],origin),global(p.vertices()[f.vertices[1]],origin),rz);
    for(const auto& f:p.faces()){
        const auto native=p.domain().face(mesh,f.nativeRow);if(native.outwardSign==0)continue;
        const auto a=global(p.vertices()[f.vertices[0]],origin),b=global(p.vertices()[f.vertices[1]],origin);
        for(unsigned k=0;k<f.rangeCount;++k){const auto& r=f.ranges[k];
            if(f.negativeCell>=0?r.minusLiquid:r.plusLiquid){const auto u=interpolate(a,b,r.begin),v=interpolate(a,b,r.end);
                const bool reverse=native.horizontal?native.outwardSign>0:native.outwardSign<0;
                greenTotal+=reverse?green(v,u,rz):green(u,v,rz);}}
    }
    near(greenTotal,liquid,full,"independent Gamma plus wet-exterior Green",arithmetic);
    near(p.report().inventory,inventory,full,"inventory report independent sum");
    near(p.report().measured,liquid,full,"geometry report independent sum");
    peakGreen=std::max(peakGreen,std::abs(greenTotal-liquid)/full);
}
void verifyTransfer(const QuadtreeMesh& sourceMesh,const Publication& source,const QuadtreeMesh& targetMesh,const Transfer& result)
{
    const auto& target=result.candidate();const bool rz=source.domain().axisymmetric();
    check(result.sourceIdentity()==source.identity(),"source identity missing");
    check(target.identity()!=source.identity(),"regrid did not create new owner");
    check(sameClock(target.timeLayer(),source.timeLayer())&&sameClock(target.timeLayer(),accepted),"regrid changed complete accepted clock");
    check(target.domain().axisymmetric()==rz,"regrid changed physical metric");
    std::vector<long double> sourceCover(sourceMesh.cellIndexNum(),0),targetCover(targetMesh.cellIndexNum(),0),contributions(targetMesh.cellIndexNum(),0);
    std::vector<long double> donorStock(sourceMesh.cellIndexNum(),0),donorCorrections(sourceMesh.cellIndexNum(),0);std::set<std::pair<int,int>> pairs;
    std::size_t unchanged=0,coarsened=0,refined=0;std::set<int> refinedDonors;
    std::vector<Kind> targetKind(targetMesh.cellIndexNum(),Kind::Unchanged);std::vector<bool> seen(targetMesh.cellIndexNum(),false);
    long double sourceInventory=0,targetInventory=0,sourceGas=0,targetGas=0,full=0,correction=0;
    for(const auto& o:result.overlaps()){
        ++overlapChecks;check(o.sourceCell>=0&&o.sourceCell<sourceMesh.cellIndexNum()&&o.targetCell>=0&&o.targetCell<targetMesh.cellIndexNum(),"overlap foreign cell index");
        check(pairs.emplace(o.sourceCell,o.targetCell).second,"duplicate cell overlap");
        const auto b=intersection(bounds(sourceMesh,o.sourceCell),bounds(targetMesh,o.targetCell));const auto cap=capacity(b,rz);
        check(cap>0,"overlap is lower dimensional");near(o.capacity,cap,cap,"physical XY/RZ overlap capacity");
        const auto l=phaseOverlap(source,o.sourceCell,b,true),g=phaseOverlap(source,o.sourceCell,b,false);
        near(o.liquidMeasured,l,cap,"independent raw actual-fan liquid");near(o.gasMeasured,g,cap,"independent raw actual-fan gas");
        near(l+g,cap,cap,"raw physical complementarity");peakRaw=std::max(peakRaw,std::abs(o.liquidMeasured-l)/cap);
        check(o.inventoryContribution>=0&&o.inventoryContribution<=cap,"contribution outside phase bounds");
        if(l==0)check(o.inventoryContribution==0,"pure gas cut gained liquid");
        if(g==0)near(o.inventoryContribution,cap,cap,"pure liquid cut lost exact capacity");
        if(o.kind==Kind::Unchanged){check(source.alpha()[o.sourceCell]==target.alpha()[o.targetCell],"unchanged exact alpha changed");
            near(o.inventoryContribution,source.cells()[o.sourceCell].inventory,cap,"unchanged stock");check(o.closureCorrection==0,"unchanged cut received closure");}
        else if(o.kind==Kind::Coarsened){near(o.inventoryContribution,source.cells()[o.sourceCell].inventory,cap,"coarsening must sum original stock");
            check(o.closureCorrection==0,"restriction invented closure correction");}
        else {refinedDonors.insert(o.sourceCell);
            if(l>0&&g>0){const auto raw=cap*l/(l+g);near(o.inventoryContribution,raw+o.closureCorrection,cap,"raw normalization plus explicit local closure");
                const auto minor=std::min(raw,cap-raw);check(std::abs(o.closureCorrection)<minor,"refine closure can destroy a phase");}
            else check(o.closureCorrection==0,"pure refined cut was adjusted");}
        donorCorrections[o.sourceCell]+=o.closureCorrection;sourceCover[o.sourceCell]+=cap;targetCover[o.targetCell]+=cap;donorStock[o.sourceCell]+=o.inventoryContribution;
        contributions[o.targetCell]+=o.inventoryContribution;correction+=o.closureCorrection;
        if(seen[o.targetCell])check(targetKind[o.targetCell]==o.kind,"target mixes nonnested transfer kinds");
        seen[o.targetCell]=true;targetKind[o.targetCell]=o.kind;
    }
    for(const auto& c:source.cells()){near(sourceCover[c.nativeCell],c.capacity,c.capacity,"complete source volume coverage");
        const auto minor=std::min(c.inventory,c.capacity-c.inventory);
        if(minor>0)check(std::abs(donorCorrections[c.nativeCell])<=tolerance*minor,"donor correction exceeds original source geometry budget");
        else check(donorCorrections[c.nativeCell]==0,"pure donor received correction");
        const auto scale=c.inventory>0?c.inventory:c.capacity;
        near(donorStock[c.nativeCell],c.inventory,scale,"source donor inventory closure",tolerance);
        sourceInventory+=c.inventory;sourceGas+=c.capacity-c.inventory;full+=c.capacity;}
    for(const auto& c:target.cells()){check(seen[c.nativeCell],"target missing donor");
        near(targetCover[c.nativeCell],c.capacity,c.capacity,"complete target volume coverage");
        near(contributions[c.nativeCell],c.inventory,c.capacity,"double alpha conversion/target contribution",tolerance);
        targetInventory+=c.inventory;targetGas+=c.capacity-c.inventory;
        if(targetKind[c.nativeCell]==Kind::Unchanged)++unchanged;else if(targetKind[c.nativeCell]==Kind::Coarsened)++coarsened;else ++refined;}
    near(targetGas,sourceGas,sourceGas>0?sourceGas:full,"independent gas inventory conservation",tolerance);
    peakGas=std::max(peakGas,std::abs(targetGas-sourceGas)/(sourceGas>0?sourceGas:full));
    const auto scale=sourceInventory>0?sourceInventory:full;
    near(targetInventory,sourceInventory,scale,"total inventory conservation",tolerance);
    peakInventory=std::max(peakInventory,std::abs(targetInventory-sourceInventory)/scale);
    const auto& r=result.report();
    peakDonorClosure=std::max(peakDonorClosure,r.maximumRelativeDonorClosureCorrection);
    peakCellClosure=std::max(peakCellClosure,r.maximumRelativeClosureCorrection);
    peakMinorFactor=std::max(peakMinorFactor,r.maximumRelativeMinorPhaseCorrection);
    check(r.refinedDonors==refinedDonors.size()&&r.refinedTargets==refined&&r.coarsenedTargets==coarsened&&r.unchangedTargets==unchanged,"actual transfer-kind report");
    check(r.sourceComponents==source.topology().components&&r.targetComponents==target.topology().components,"component report detached from true Gamma");
    check(r.sourceJunctions==source.topology().junctionVertices.size()&&r.targetJunctions==target.topology().junctionVertices.size(),"junction report detached from true Gamma");
    near(r.sourceGasInventory,sourceGas,full,"source gas report");near(r.targetGasInventory,targetGas,full,"target gas report");
    check(r.maximumRelativeMinorPhaseCorrection<1,"minor factor can destroy a phase");
    near(r.sourceInventory,sourceInventory,full,"source inventory report");near(r.targetInventory,targetInventory,full,"target inventory report");
    near(r.sourceMeasured,source.report().measured,full,"source geometry report");near(r.targetMeasured,target.report().measured,full,"target geometry report");
    near(r.closureCorrection,correction,full,"explicit closure correction total");
    check(r.maximumRelativeCoverageResidual<=tolerance&&r.maximumRelativeInventoryResidual<=tolerance&&r.maximumRelativeGeometryResidual<=tolerance&&r.maximumRelativeDonorClosureCorrection<=tolerance,"reported frozen gate exceeded");
    verifyGeometry(targetMesh,target);
}
Transfer transfer(const QuadtreeMesh& sourceMesh,const Publication& source,const QuadtreeMesh& targetMesh)
{ return Transfer::prepare(sourceMesh,source,targetMesh,domain(targetMesh,source.domain().axisymmetric()),accepted); }
void staticTransfers()
{
    for(bool rz:{false,true}){
        stage=rz?"RZ native 2:1 transfers":"XY native 2:1 transfers";
        auto coarse=meshFor(),fine=meshFor(8,1,true);check(!domain(*fine,rz).uniformLeaves(),"fixture is not true native 2:1");
        const auto source=publish(*coarse,rz,lineInput(*coarse,rz));const auto old=source.alpha();const auto id=source.identity();
        auto identical=coarse->clone();const auto same=transfer(*coarse,source,*identical);verifyTransfer(*coarse,source,*identical,same);
        check(same.report().unchangedTargets==std::size_t(coarse->cellIndexNum()),"unchanged fixture unexpectedly adapted");
        const auto refined=transfer(*coarse,source,*fine);verifyTransfer(*coarse,source,*fine,refined);
        check(refined.report().refinedTargets>0&&refined.report().unchangedTargets>0,"local refine lacks real refined/unchanged leaves");
        const auto restricted=transfer(*fine,refined.candidate(),*coarse);verifyTransfer(*fine,refined.candidate(),*coarse,restricted);
        check(restricted.report().coarsenedTargets>0,"real restriction absent");
        check(source.identity()==id&&source.alpha()==old,"successful candidate mutated source");scenarios+=3;
    }
}
void adaptTo(QuadtreeMesh& mesh,int level)
{
    CellField<double> aim(mesh,0.),gradient(mesh,0.);
    for(int i=0;i<mesh.cellIndexNum();++i)aim.value[i]=mesh.cells(mesh.cellIndexToLabel[i]).centerCoord().x>=.5?level:0;
    check(mesh.mayAdaptToTargets(aim),"actual adapt preflight found no topology change");
    std::vector<VertField<double>> levelSets;std::vector<VertField<double>*> vertices;
    std::vector<HEdgeField<double>*> horizontal;std::vector<VEdgeField<double>*> vertical;std::vector<CellField<double>*> cells{&aim};
    mesh.adapt(levelSets,gradient,vertices,horizontal,vertical,cells,&aim);++adaptations;
}
void repeatedAdaptation()
{
    for(bool rz:{false,true}){
        stage=rz?"RZ repeated actual adapt":"XY repeated actual adapt";auto mesh=meshFor();
        // Freeze a plane away from finest-lattice corners. The aligned
        // 29/32 case creates a genuinely positive sub-double gas sliver after
        // reconstruction; its atomic refusal is tested separately below.
        auto state=publish(*mesh,rz,input(*mesh,rz,[](Box b,bool radial){return lineVolume(b,radial,.913L);}));
        const auto first=state.report().inventory;
        for(int desired:{1,0,2,0}){
            auto nextMesh=mesh->clone();const auto oldCount=mesh->cellIndexNum();const auto id=state.identity();const auto saved=state.alpha();
            adaptTo(*nextMesh,desired);check(nextMesh->cellIndexNum()!=oldCount,"actual adapt did not change leaf layout");
            const auto candidate=transfer(*mesh,state,*nextMesh);verifyTransfer(*mesh,state,*nextMesh,candidate);
            check(state.identity()==id&&state.alpha()==saved&&mesh->cellIndexNum()==oldCount,"candidate staging polluted accepted mesh/state");
            state=candidate.candidate();mesh=std::move(nextMesh);
            near(state.report().inventory,first,first,"multi-round stock closure",4*tolerance);++scenarios;
        }
        check(mesh->cellIndexNum()==64,"coarsening did not return actual base leaves");
    }
}
void pureAndTiny()
{
    for(bool rz:{false,true})for(double a:{0.,1.}){
        stage=std::string(rz?"RZ":"XY")+" exact pure transfer";auto coarse=meshFor(4),fine=meshFor(4,1,true);
        const auto source=publish(*coarse,rz,std::vector<double>(coarse->cellIndexNum(),a));const auto result=transfer(*coarse,source,*fine);
        verifyTransfer(*coarse,source,*fine,result);for(double q:result.candidate().alpha())check(q==a,"pure alpha changed");
        check(result.candidate().facets().empty(),"pure phase invented Gamma");const auto back=transfer(*fine,result.candidate(),*coarse);verifyTransfer(*fine,result.candidate(),*coarse,back);for(double q:back.candidate().alpha())check(q==a,"pure restriction alpha changed");scenarios+=2;
    }
    for(bool rz:{false,true}){
        stage=std::string(rz?"RZ":"XY")+" resolvable tiny phase";auto coarse=meshFor(4),fine=meshFor(4,1);
        constexpr long double width=1e-9L;
        const auto alpha=input(*coarse,rz,[width](Box b,bool radial){const auto end=std::clamp(width,b.x0,b.x1);
            return (b.y1-b.y0)*(radial?twoPi*(end*end-b.x0*b.x0)/2:end-b.x0);});
        const auto source=publish(*coarse,rz,alpha);const auto result=transfer(*coarse,source,*fine);verifyTransfer(*coarse,source,*fine,result);
        check(result.candidate().report().inventory>0&&result.candidate().report().measured>0,"tiny phase swallowed");
        const auto back=transfer(*fine,result.candidate(),*coarse);verifyTransfer(*fine,result.candidate(),*coarse,back);scenarios+=2;
    }
}
void smallPhysicalDomain()
{
    constexpr double extent=1e-9;
    for(bool rz:{false,true})for(double a:{0.,1.}){
        stage=std::string(rz?"RZ":"XY")+" small physical-domain transfer";
        DictionaryMesh d;d.Nx=d.Ny=4;d.maxCoord={extent,extent};d.maxInitialIterations=0;d.maxCellNum=1024;d.maxLevel_local={1};
        auto coarse=std::make_unique<QuadtreeMesh>(d);const Vector2D point{extent/8,extent/8};
        const auto exact=coarse->findCellLabel(point,0.);
        check(exact==coarse->cellIndexToLabel[0],"zero-tolerance lookup missed first small-domain leaf");
        check(coarse->findCellLabel(point)==coarse->findCellLabel(point,MINVAL),"old lookup changed default boundary tolerance");
        check(coarse->findCellLabel(point)!=exact,"small-domain fixture does not expose old absolute tolerance");
        check(coarse->findCellLabel({0,0},0.)==coarse->cellIndexToLabel[0],"exact domain minimum lookup failed");
        check(coarse->findCellLabel({extent,extent},0.)>=0,"exact domain maximum lookup failed");
        check(coarse->findCellLabel({-extent/8,extent/2},0.)==-1&&coarse->findCellLabel({extent+extent/8,extent/2},0.)==-1,"zero-tolerance lookup accepted outside point");
        check(coarse->findCellLabel({-extent/8,extent/2})>=0,"legacy boundary allowance unexpectedly changed");
        const auto nan=std::numeric_limits<double>::quiet_NaN(),inf=std::numeric_limits<double>::infinity();
        check(coarse->findCellLabel({nan,extent/2},0.)==-1&&coarse->findCellLabel({extent/2,inf},0.)==-1,"nonfinite coordinate entered lookup");
        check(coarse->findCellLabel(point,-1.)==-1&&coarse->findCellLabel(point,nan)==-1&&coarse->findCellLabel(point,inf)==-1,"invalid explicit tolerance entered lookup");
        auto fine=coarse->clone();DictionaryMesh::RefineBox region;region.level=1;region.min={extent/2,0};region.max={extent,extent};fine->applyRefineBoxes({region});
        check(!domain(*fine,rz).uniformLeaves(),"small-domain fixture lost true 2:1 mesh");
        for(int i=0;i<fine->cellIndexNum();++i){const auto& cell=fine->cells(fine->cellIndexToLabel[i]);
            check(fine->findCellLabel(cell.centerCoord(),0.)==cell.label(),"zero-tolerance lookup missed actual fine center");}
        const auto source=publish(*coarse,rz,std::vector<double>(coarse->cellIndexNum(),a));
        auto identical=coarse->clone();const auto same=transfer(*coarse,source,*identical);verifyTransfer(*coarse,source,*identical,same);
        const auto refined=transfer(*coarse,source,*fine);verifyTransfer(*coarse,source,*fine,refined);
        const auto back=transfer(*fine,refined.candidate(),*coarse);verifyTransfer(*fine,refined.candidate(),*coarse,back);
        for(double q:back.candidate().alpha())check(q==a,"small-domain pure phase changed");scenarios+=3;
    }
}
void syntheticTopology()
{
    for(int shape:{0,1}){
        stage=shape==0?"synthetic reconstruction merge":"synthetic reconstruction split";
        auto fine=meshFor(8,2),coarse=meshFor(8);std::vector<Box> liquids;
        if(shape==0)liquids={{13.L/32,15.L/32,15.L/32,17.L/32},{17.L/32,19.L/32,15.L/32,17.L/32}};
        else liquids={{8.L/32,15.L/32,8.L/32,24.L/32},{17.L/32,24.L/32,8.L/32,24.L/32},{15.L/32,17.L/32,15.L/32,16.L/32}};
        const auto alpha=input(*fine,false,[&](Box b,bool rz){long double v=0;for(auto q:liquids)v+=capacity(intersection(b,q),rz);return v;});
        for(double q:alpha)check(q==0||q==1,"synthetic source must be exact native pure inventory");
        const auto source=publish(*fine,false,alpha);check(source.topology().components==(shape==0?2u:1u),"synthetic source Gamma graph incorrect");
        const auto result=transfer(*fine,source,*coarse);verifyTransfer(*fine,source,*coarse,result);
        check(result.candidate().topology().components==(shape==0?1u:2u),"frozen synthetic reconstruction topology not observed");
        check(result.report().sourceComponents!=result.report().targetComponents,"transfer incorrectly enforces component bijection");++scenarios;
    }
}
void rejectionTransaction()
{
    stage="full-clock/domain/stale/failure transaction";auto mesh=meshFor(4),target=meshFor(4,1,true);
    const auto source=publish(*mesh,false,lineInput(*mesh,false));const auto id=source.identity();const auto saved=source.alpha();
    const auto* address=&source.alpha();const auto count=mesh->cellIndexNum(),targetCount=target->cellIndexNum();
    const auto targetDomain=domain(*target,false);
    for(Clock wrong:{Clock{8,.125,.25,Position::End},Clock{7,0,.25,Position::End},Clock{7,0,.5,Position::Midpoint}})
        rejects([&]{Transfer::prepare(*mesh,source,*target,targetDomain,wrong);},"full accepted clock mismatch");
    const auto initial=publish(*mesh,false,lineInput(*mesh,false),{});
    rejects([&]{Transfer::prepare(*mesh,initial,*target,targetDomain,{});},"not accepted End");
    rejects([&]{Transfer::prepare(*mesh,source,*target,domain(*target,true),accepted);},"XY/RZ metric change");
    FaceSegmentGraph partial;partial.build(*target,[](const QuadtreeMesh&,const Cell& c){return c.centerCoord().x<.5;});
    const auto omitted=NativeLeafDomain::build(*target,partial,false);
    rejects([&]{Transfer::prepare(*mesh,source,*target,omitted,accepted);},"partial target coverage");
    auto foreign=meshFor(8);rejects([&]{transfer(*mesh,source,*foreign);},"foreign background lattice");
    DictionaryMesh d;d.Nx=d.Ny=4;d.maxInitialIterations=0;d.maxCoord={2,1};QuadtreeMesh stretched(d);
    rejects([&]{transfer(*mesh,source,stretched);},"physical domain change");
    auto staleSource=mesh->clone();staleSource->verts(staleSource->vertIndexToLabel[0]).setCoord({-.125,0});
    rejects([&]{Transfer::prepare(*staleSource,source,*target,targetDomain,accepted);},"stale source geometry");
    auto staleTarget=mesh->clone();const auto staleDomain=domain(*staleTarget,false);adaptTo(*staleTarget,1);
    rejects([&]{Transfer::prepare(*mesh,source,*staleTarget,staleDomain,accepted);},"stale target domain");
    GeometricVofConfig bad;bad.maximumVolumeIterations=0;
    rejects([&]{Transfer::prepare(*mesh,source,*target,targetDomain,accepted,bad);},"invalid root budget");
    bad={};bad.relativeVolumeTolerance=2e-12L;
    rejects([&]{Transfer::prepare(*mesh,source,*target,targetDomain,accepted,bad);},"cannot loosen conservation gate");
    bad={};bad.maximumVolumeIterations=1;
    rejects([&]{Transfer::prepare(*mesh,source,*target,targetDomain,accepted,bad);},"failed target reconstruction budget");
    check(source.identity()==id&&source.alpha()==saved&&&source.alpha()==address&&mesh->cellIndexNum()==count&&target->cellIndexNum()==targetCount,"failed transfer replaced accepted state/mesh");
    source.validateFor(*mesh,id,accepted);verifyGeometry(*mesh,source);++scenarios;
    // Restriction of an exact native checkerboard gives .5 everywhere. The
    // target has no resolved alpha direction, so publication must fail.
    stage="post-inventory reconstruction failure keeps source";auto fine=meshFor(4,1),coarse=meshFor(4);
    std::vector<double> checker(fine->cellIndexNum());
    for(int i=0;i<fine->cellIndexNum();++i){const auto b=bounds(*fine,i);checker[i]=(static_cast<int>(8*b.x0)+static_cast<int>(8*b.y0))%2;}
    const auto acceptedChecker=publish(*fine,false,checker);const auto checkerId=acceptedChecker.identity();
    rejects([&]{transfer(*fine,acceptedChecker,*coarse);},"homogeneous mixed target cannot invent normal");
    check(acceptedChecker.identity()==checkerId&&acceptedChecker.alpha()==checker&&coarse->cellIndexNum()==16,"target reconstruction failure modified source");
    acceptedChecker.validateFor(*fine,checkerId,accepted);++scenarios;
    stage="aligned refinement cannot erase sub-double gas";
    auto alignedMesh=meshFor(8),alignedFine=alignedMesh->clone();adaptTo(*alignedFine,2);
    const auto aligned=publish(*alignedMesh,false,lineInput(*alignedMesh,false));const auto alignedId=aligned.identity();
    const auto alignedAlpha=aligned.alpha();
    std::string refusal;
    try{transfer(*alignedMesh,aligned,*alignedFine);}catch(const std::exception& error){refusal=error.what();}
    check(refusal.find("unrepresentable")!=std::string::npos,"aligned cut did not reject for the actual representation limit");++guards;
    std::cout<<"GEOMETRIC_VOF_EXPECTED_REFUSAL "<<refusal<<'\n';
    check(aligned.identity()==alignedId&&aligned.alpha()==alignedAlpha&&alignedMesh->cellIndexNum()==64,
        "unrepresentable cut failure changed accepted state");aligned.validateFor(*alignedMesh,alignedId,accepted);++scenarios;
}
namespace LocalPlicContract
{
using Stocks=std::vector<std::array<long double,2>>;
unsigned cases=0,piecesChecked=0;long double peakOwn=0,peakJoint=0,peakRawPlane=0;
static_assert(std::is_same_v<decltype(std::declval<const Transfer&>().nativeOverlaps()),const LAMP::Mesh::NativeLeafOverlap&>);
Stocks retained(const Publication& p)
{ Stocks result;for(const auto& c:p.plicCells())result.push_back(c.inventory);return result; }
Publication publishStocks(const QuadtreeMesh& mesh,bool rz,const Stocks& stocks,Clock clock=accepted)
{ return Publication::prepareLocalPlicInventory(mesh,domain(mesh,rz),clock,stocks,{}); }
Stocks independentPlane(const QuadtreeMesh& mesh,bool rz,long double x)
{
    Stocks stocks(mesh.cellIndexNum());
    for(int i=0;i<mesh.cellIndexNum();++i){const auto b=bounds(mesh,i);
        const Box l{b.x0,std::min(b.x1,x),b.y0,b.y1},g{std::max(b.x0,x),b.x1,b.y0,b.y1};
        stocks[i]={capacity(l,rz),capacity(g,rz)};}
    return stocks;
}
int namedCell(const QuadtreeMesh& mesh)
{
    for(int i=0;i<mesh.cellIndexNum();++i){const auto b=bounds(mesh,i);if(b.x0==.5L&&b.y0==.5L)return i;}
    throw std::runtime_error("missing named test source cell");
}
void ownNear(long double actual,long double expected,const std::string& name)
{
    if(expected==0){check(actual==0,name+" zero source phase gained stock");return;}
    check(actual>0,name+" lost positive phase");const auto residual=std::abs(actual-expected)/expected;
    peakOwn=std::max(peakOwn,residual);check(residual<=tolerance,name+" own-phase relative error");
}
long double independentCut(const GeometricVofPlicCell& source,Box box,bool liquid,bool rz)
{
    const Point origin{(static_cast<long double>(source.rectangle.lower.x)+source.rectangle.upper.x)/2,
        (static_cast<long double>(source.rectangle.lower.y)+source.rectangle.upper.y)/2};
    std::vector<Point> polygon{{box.x0-origin.x,box.y0-origin.y},{box.x1-origin.x,box.y0-origin.y},
        {box.x1-origin.x,box.y1-origin.y},{box.x0-origin.x,box.y1-origin.y}};
    const auto sign=liquid?1.L:-1.L,nx=sign*source.normal.x,ny=sign*source.normal.y;
    const auto limit=sign*source.offset*std::hypot(static_cast<long double>(source.normal.x),static_cast<long double>(source.normal.y));
    const auto value=[&](Point p){return nx*p.x+ny*p.y-limit;};
    std::vector<Point> out;auto a=polygon.back();auto va=value(a);bool ia=va<=0;
    for(auto b:polygon){const auto vb=value(b);const bool ib=vb<=0;
        if(ia!=ib){const auto t=va/(va-vb);out.push_back({a.x+t*(b.x-a.x),a.y+t*(b.y-a.y)});}
        if(ib)out.push_back(b);a=b;va=vb;ia=ib;}
    return measure(out,origin,rz);
}
// D17 test-owned integer oracle. All existing LocalPlic fixtures and the
// actual D16 rectangle lie exactly on this fixed 1/128 physical lattice.
// No production exact-capacity or partition helper supplies these weights.
std::uint64_t latticeWeight(Box b,bool rz)
{
    const auto integer=[](long double x){
        const auto q=x*128;const auto k=static_cast<std::int64_t>(q);
        check(q==k,"independent CDF oracle coordinate is not on actual fixture lattice");return k;
    };
    const auto x0=integer(b.x0),x1=integer(b.x1),y0=integer(b.y0),y1=integer(b.y1);
    check(x1>x0&&y1>y0&&(!rz||x0>=0),"independent CDF oracle physical rectangle");
    return static_cast<std::uint64_t>((rz?x1*x1-x0*x0:x1-x0)*(y1-y0));
}
std::uint64_t integerCdf(std::uint64_t m,std::uint64_t prefix,std::uint64_t whole)
{
    check(whole>0&&prefix<=whole,"independent CDF oracle prefix");
    // Divide before multiplying: m can fill the native significand, while
    // the actual fixture's small integer geometric weights fit uint64_t.
    const auto remainder=m%whole;
    check(prefix==0||remainder<=std::numeric_limits<std::uint64_t>::max()/prefix,"integer oracle product overflow");
    const auto product=remainder*prefix;
    auto count=(m/whole)*prefix+product/whole;const auto tail=product%whole;
    if(tail>whole-tail||(tail==whole-tail&&(count&1)))++count;
    check(count<=m,"independent CDF oracle count exceeds source");return count;
}
long double independentCdf(const QuadtreeMesh& mesh,const Publication& source,
    const Transfer& result,std::size_t selected,int phase)
{
    check(selected<result.overlaps().size()&&phase>=0&&phase<2,"CDF oracle selected piece/phase index");
    const auto si=result.overlaps()[selected].sourceCell;
    check(si>=0&&static_cast<std::size_t>(si)<source.plicCells().size()
        &&static_cast<std::size_t>(si)<result.nativeOverlaps().bySource().size(),"CDF oracle actual source index");
    const auto stock=source.plicCells()[si].inventory[phase];
    if(stock==0)return 0;
    check(std::isfinite(stock)&&stock>0,"independent CDF oracle source stock");
    int exponent=0;std::frexp(stock,&exponent);
    const auto quantum=std::fpclassify(stock)==FP_SUBNORMAL?std::numeric_limits<long double>::denorm_min():
        std::ldexp(1.L,exponent-std::numeric_limits<long double>::digits);
    const auto nativeCount=stock/quantum;
    check(quantum>0&&nativeCount>0&&nativeCount<std::ldexp(1.L,64)
        &&nativeCount==std::floor(nativeCount),"independent native quantum/count");
    const auto m=static_cast<std::uint64_t>(nativeCount);
    check(static_cast<long double>(m)*quantum==stock,"integer oracle source is not exactly represented");
    auto ids=result.nativeOverlaps().bySource()[si];
    const auto& rectangles=result.nativeOverlaps().targetRectangles();
    check(!ids.empty(),"CDF oracle empty actual source partition");
    for(auto id:ids){check(id<result.overlaps().size(),"CDF oracle grouped piece index");
        const auto& o=result.overlaps()[id];
        check(o.sourceCell==si&&o.kind==Kind::Refined&&o.targetCell>=0
            &&static_cast<std::size_t>(o.targetCell)<rectangles.size()
            &&o.targetCell<mesh.cellIndexNum(),"CDF oracle grouped source/target/kind mismatch");}
    std::sort(ids.begin(),ids.end(),[&](std::size_t a,std::size_t b){
        const auto& x=rectangles[result.overlaps()[a].targetCell];const auto& y=rectangles[result.overlaps()[b].targetCell];
        return std::tie(x.lower.x,x.lower.y,x.upper.x,x.upper.y)<
            std::tie(y.lower.x,y.lower.y,y.upper.x,y.upper.y);
    });
    const bool rz=source.domain().axisymmetric();std::uint64_t whole=0;
    for(auto id:ids)whole+=latticeWeight(bounds(mesh,result.overlaps()[id].targetCell),rz);
    // mesh is the TARGET mesh; source rectangle is independently retained.
    const auto& parent=source.plicCells()[si].rectangle;
    check(whole==latticeWeight({parent.lower.x,parent.upper.x,parent.lower.y,parent.upper.y},rz),
        "independent CDF integer geometric partition does not cover parent");
    std::uint64_t prefix=0,previous=0;
    for(auto id:ids){prefix+=latticeWeight(bounds(mesh,result.overlaps()[id].targetCell),rz);
        const auto current=integerCdf(m,prefix,whole);
        if(id==selected)return static_cast<long double>(current-previous)*quantum;
        previous=current;
    }
    throw std::runtime_error("missing selected actual CDF piece");
}
void verify(const QuadtreeMesh& sourceMesh,const Publication& source,const QuadtreeMesh& targetMesh,const Transfer& result)
{
    const auto& target=result.candidate();const auto& cover=result.nativeOverlaps();const bool rz=source.domain().axisymmetric();
    check(result.report().localPlic,"LocalPlic report unavailable");
    check(target.geometryKind()==GeometricVofGeometryKind::LocalPlic,"target changed VOF geometry kind");
    check(result.sourceIdentity()==source.identity()&&target.identity()!=source.identity(),"source/fresh target identity");
    check(sameClock(target.timeLayer(),accepted)&&sameClock(source.timeLayer(),accepted),"same accepted End clock");
    cover.validateFor(sourceMesh,targetMesh);
    check(cover.pieces().size()==result.overlaps().size(),"native coverage/piece record mismatch");
    Stocks bySource(source.plicCells().size()),byTarget(target.plicCells().size());
    for(std::size_t i=0;i<result.overlaps().size();++i){const auto& o=result.overlaps()[i];const auto& p=cover.pieces()[i];
        check(o.localPlic&&o.sourceCell==p.sourceCell&&o.targetCell==p.targetCell,"actual native overlap association");
        const auto b=o.kind==Kind::Refined?bounds(targetMesh,o.targetCell):bounds(sourceMesh,o.sourceCell);
        near(o.capacity,capacity(b,rz),o.capacity,"independent overlap capacity");
        const auto& donor=source.plicCells()[o.sourceCell];
        const bool finitePartition=o.kind==Kind::Refined&&
            (!donor.resolved||donor.inventory[0]==0||donor.inventory[1]==0);
        for(std::size_t phase=0;phase<2;++phase){const auto value=o.phaseContribution[phase];
            check(std::isfinite(value)&&value>=0,"nonnegative independent overlap stock");
            if(finitePartition)check(value==independentCdf(targetMesh,source,result,i,static_cast<int>(phase)),
                "declared finite pure/P0 differs from independent integer CDF");
            else if(o.kind!=Kind::Refined)check(o.phaseContribution==donor.inventory,
                "nonrefined certified stock was changed");
            else check(value<=o.capacity,"ordinary resolved overlap stock exceeds capacity");
            bySource[o.sourceCell][phase]+=value;byTarget[o.targetCell][phase]+=value;}
        const auto joint=std::abs((o.phaseContribution[0]+o.phaseContribution[1])-o.capacity)/o.capacity;
        peakJoint=std::max(peakJoint,joint);check(joint<=tolerance,"piece independent phase closure");
        if(o.kind==Kind::Refined&&donor.resolved){
            check(!o.usedP0&&o.geometricMeasureAvailable,"resolved source silently fell back to P0");
            const auto l=independentCut(donor,b,true,rz),g=independentCut(donor,b,false,rz);
            near(o.liquidMeasured,l,o.capacity,"independent parent plane liquid cut",2e-12L);
            near(o.gasMeasured,g,o.capacity,"independent parent plane gas cut",2e-12L);
            peakRawPlane=std::max({peakRawPlane,std::abs(o.liquidMeasured-l)/o.capacity,std::abs(o.gasMeasured-g)/o.capacity});
        }
        if(o.kind==Kind::Refined&&!donor.resolved){
            check(o.usedP0&&!o.geometricMeasureAvailable,"unresolved P0 tag lost");
            for(std::size_t phase=0;phase<2;++phase){
                const auto expected=independentCdf(targetMesh,source,result,i,static_cast<int>(phase));
                check(o.phaseContribution[phase]==expected,"declared independent P0 integer allocation");}
        }
        if(o.kind==Kind::Unchanged)check(o.phaseContribution==donor.inventory,"unchanged two stocks not bitwise retained");
        ++piecesChecked;
    }
    std::array<long double,2> sourceTotal{},targetTotal{};
    for(std::size_t i=0;i<bySource.size();++i)for(std::size_t phase=0;phase<2;++phase){
        ownNear(bySource[i][phase],source.plicCells()[i].inventory[phase],"independent donor sum");
        sourceTotal[phase]+=source.plicCells()[i].inventory[phase];}
    for(std::size_t i=0;i<byTarget.size();++i){
        const auto& stocks=target.plicCells()[i].inventory;
        if(target.domain().cellRowByIndex()[i]<0){
            check(stocks[0]==0&&stocks[1]==0&&byTarget[i][0]==0&&byTarget[i][1]==0,
                "inactive target acquired independent stock");continue;
        }
        const auto full=capacity(bounds(targetMesh,static_cast<int>(i)),rz);
        for(std::size_t phase=0;phase<2;++phase){ownNear(stocks[phase],byTarget[i][phase],"target retained independent sum");targetTotal[phase]+=stocks[phase];}
        check(std::abs((stocks[0]+stocks[1])-full)<=tolerance*full,"target capacity closure");}
    for(std::size_t phase=0;phase<2;++phase){ownNear(targetTotal[phase],sourceTotal[phase],"independent global sum");
        check(result.report().maximumRelativeDonorPhaseResidual[phase]<=tolerance&&result.report().relativeGlobalPhaseResidual[phase]<=tolerance,"own-phase report gate");}
    rejects([&]{target.validateFor(targetMesh,source.identity(),accepted);},"old identity cannot name new candidate");
    source.validateFor(sourceMesh,source.identity(),accepted);target.validateFor(targetMesh,target.identity(),accepted);
}
// Same physical rectangle as D16 target972; only background level labels
// become 0/1. Source refinement and coarsening are genuine public mesh calls.
void capacityCycles()
{
    namespace FV=LAMP::Geometry::Detail::FiniteVolume;
    stage="D17 actual RZ pure capacity and saturated integer partition";
    DictionaryMesh d;d.Nx=d.Ny=1;d.minCoord={11./64,5./16};d.maxCoord={3./16,21./64};
    d.maxCellNum=64;d.maxInitialIterations=0;d.maxLevel_local={1};
    auto original=std::make_unique<QuadtreeMesh>(d);auto fine=original->clone();
    DictionaryMesh::RefineBox box;box.level=1;box.min=d.minCoord;box.max=d.maxCoord;
    fine->applyRefineBoxes({box});check(original->cellIndexNum()==1&&fine->cellIndexNum()==4,"D17 actual 1-to-4 tree fixture");
    auto coarse=fine->clone();
    CellField<double> aim(*coarse,0.),gradient(*coarse,0.);
    check(coarse->mayAdaptToTargets(aim),"D17 actual coarsening preflight");
    std::vector<VertField<double>> levelSets;std::vector<VertField<double>*> vertices;
    std::vector<HEdgeField<double>*> horizontal;std::vector<VEdgeField<double>*> vertical;std::vector<CellField<double>*> cells{&aim};
    coarse->adapt(levelSets,gradient,vertices,horizontal,vertical,cells,&aim);
    check(coarse->cellIndexNum()==1,"D17 actual four siblings did not coarsen");
    for(int phase=0;phase<2;++phase){
        Stocks pure(fine->cellIndexNum());
        for(int c=0;c<fine->cellIndexNum();++c)pure[c][phase]=capacity(bounds(*fine,c),true);
        const auto source=publishStocks(*fine,true,pure);const auto saved=retained(source);const auto sourceId=source.identity();
        const auto down=transfer(*fine,source,*coarse);verify(*fine,source,*coarse,down);
        const auto parent=retained(down.candidate());const auto q=parent[0][phase],cap=down.candidate().cells()[0].capacity;
        check(parent[0][1-phase]==0&&down.candidate().alpha()[0]==(phase==0?1:0),"D17 pure phase/diagnostic alpha changed");
        if(std::numeric_limits<long double>::digits==53){
            const auto expectedC=std::ldexp(static_cast<long double>(0x12106ca4910069ULL),-64);
            const auto expectedQ=std::ldexp(static_cast<long double>(0x12106ca491006aULL),-64);
            check(cap==expectedC&&q==expectedQ&&q>cap,"D17 physical D16 capacity/stock bit witness");
            check(down.candidate().capacityDefects()[0].compare(FV::Dyadic(std::ldexp(1.L,-64)))==0,
                "D17 retained one-ULP certificate does not match independent hex oracle");
            rejects([&]{publishStocks(*coarse,true,parent);},"D17 strict raw root cannot self-sign certified excess");
        }
        auto unchanged=coarse->clone();const auto same=transfer(*coarse,down.candidate(),*unchanged);
        verify(*coarse,down.candidate(),*unchanged,same);
        check(retained(same.candidate())==parent,"D17 unchanged raw certified inventory differs");
        const auto up=transfer(*coarse,down.candidate(),*fine);verify(*coarse,down.candidate(),*fine,up);
        const auto first=retained(up.candidate());
        const auto back=transfer(*fine,up.candidate(),*coarse);verify(*fine,up.candidate(),*coarse,back);
        check(retained(back.candidate())==parent,"D17 saturated refine/coarsen loses raw parent phase");
        const auto again=transfer(*coarse,back.candidate(),*fine);verify(*coarse,back.candidate(),*fine,again);
        check(retained(again.candidate())==first,"D17 identical integer-CDF cycle is not idempotent");
        check(source.identity()==sourceId&&retained(source)==saved,"D17 capacity cycle changed source");
        std::cout<<std::setprecision(std::numeric_limits<long double>::max_digits10)
            <<"D17_CAPACITY_PHASE phase="<<phase<<" capacity="<<cap<<" stock="<<q
            <<" excess="<<q-cap<<" children="<<fine->cellIndexNum()<<'\n';
        cases+=5;
    }
    std::cout<<"LOCAL_PLIC_CAPACITY_PARTITION_PASS phases=2 same_actual_partition cycles=2"<<std::endl;
}
void group()
{
    for(bool rz:{false,true}){
        stage=rz?"LocalPlic RZ true refine/coarsen/unchanged":"LocalPlic XY true refine/coarsen/unchanged";
        auto coarse=meshFor(4);auto stocks=independentPlane(*coarse,rz,.734375L);const auto named=namedCell(*coarse);
        // A declared legal geometry/stock defect forces donor-local correction;
        // both inputs remain independent stocks and alpha is unchanged.
        stocks[named][1]+=4e-14L*capacity(bounds(*coarse,named),rz);
        const auto source=publishStocks(*coarse,rz,stocks);const auto saved=retained(source);const auto id=source.identity();
        auto same=coarse->clone();const auto unchanged=transfer(*coarse,source,*same);verify(*coarse,source,*same,unchanged);
        check(unchanged.report().unchangedTargets==static_cast<std::size_t>(same->cellIndexNum()),"unchanged coverage count");
        auto fine=coarse->clone();adaptTo(*fine,1);
        const auto refined=transfer(*coarse,source,*fine);verify(*coarse,source,*fine,refined);
        check(refined.report().resolvedRefinedDonors>0&&refined.report().refinedTargets>0&&refined.report().unchangedTargets>0,"actual refine and unchanged neighbors absent");
        check(refined.report().correctedRefinedDonors>0,"declared stock defect did not exercise conservative correction");
        auto back=fine->clone();adaptTo(*back,0);
        const auto restricted=transfer(*fine,refined.candidate(),*back);verify(*fine,refined.candidate(),*back,restricted);
        check(restricted.report().coarsenedTargets>0,"actual coarsen absent");
        check(source.identity()==id&&retained(source)==saved,"candidate changed source accepted stocks");cases+=3;
        stage=rz?"LocalPlic RZ independent micro phase P0":"LocalPlic XY independent micro phase P0";
        Stocks micro(coarse->cellIndexNum());for(int i=0;i<coarse->cellIndexNum();++i)micro[i]={capacity(bounds(*coarse,i),rz),0};
        micro[named][1]=micro[named][0]*1e-40L;
        const auto tiny=publishStocks(*coarse,rz,micro);const auto tinyId=tiny.identity();const auto tinySaved=retained(tiny);
        check(tiny.alpha()[named]==1&&tiny.plicCells()[named].inventory[1]>0&&!tiny.plicCells()[named].resolved,"independent sub-double gas fixture");
        const auto tinyFine=transfer(*coarse,tiny,*fine);verify(*coarse,tiny,*fine,tinyFine);
        check(tinyFine.report().unresolvedRefinedDonors>0&&tinyFine.report().p0RefinedPieces>0,"P0 counters absent");
        for(auto index:tinyFine.nativeOverlaps().bySource()[named]){
            const auto& o=tinyFine.overlaps()[index];check(o.usedP0&&o.phaseContribution[1]>0,"micro phase erased in child");}
        const auto tinyBack=transfer(*fine,tinyFine.candidate(),*back);verify(*fine,tinyFine.candidate(),*back,tinyBack);
        ownNear(tinyBack.candidate().plicCells()[namedCell(*back)].inventory[1],micro[named][1],"micro phase round trip");cases+=2;
        const auto targetDomain=domain(*fine,rz);Clock wrong=accepted;wrong.end=.26;
        rejects([&]{Transfer::prepare(*coarse,tiny,*fine,targetDomain,wrong);},"LocalPlic complete clock mismatch");
        auto stale=coarse->clone();const auto staleDomain=domain(*stale,rz);adaptTo(*stale,1);
        rejects([&]{Transfer::prepare(*coarse,tiny,*stale,staleDomain,accepted);},"LocalPlic stale target geometry");
        check(tiny.identity()==tinyId&&retained(tiny)==tinySaved,"LocalPlic rejected transaction changed source");
        tiny.validateFor(*coarse,tinyId,accepted);
    }
    stage="LocalPlic no mixed overlap slack refuses non-arithmetic repair";
    auto mesh=meshFor(4),fine=mesh->clone();adaptTo(*fine,1);auto stocks=independentPlane(*mesh,false,.625L);
    const auto named=namedCell(*mesh);stocks[named][1]+=4e-13L*capacity(bounds(*mesh,named),false);
    const auto aligned=publishStocks(*mesh,false,stocks);const auto alignedId=aligned.identity();const auto saved=retained(aligned);
    std::string reason;try{transfer(*mesh,aligned,*fine);}catch(const std::exception& e){reason=e.what();}
    check(reason.find("no mixed overlap slack")!=std::string::npos,"aligned pure children did not reject absent legal slack");
    check(aligned.identity()==alignedId&&retained(aligned)==saved,"no-slack failure modified source");
    const auto initial=publishStocks(*mesh,false,independentPlane(*mesh,false,.734375L),{});
    rejects([&]{Transfer::prepare(*mesh,initial,*fine,domain(*fine,false),{});},"LocalPlic initialization is not accepted End");
    capacityCycles();
    std::cout<<std::setprecision(17)<<"LOCAL_PLIC_MESH_TRANSFER_PASS cases="<<cases<<" pieces="<<piecesChecked
        <<" peak_own_phase="<<peakOwn<<" peak_joint="<<peakJoint<<" peak_independent_plane="<<peakRawPlane<<'\n';
}
}


namespace SelectedDomainAmrContract
{
using Stocks=LocalPlicContract::Stocks;
using Exact=LAMP::Geometry::Detail::FiniteVolume::Dyadic;
unsigned changes=0,activeChecks=0,solidChecks=0,coverageChecks=0,negativeChecks=0;
long double peakCoverage=0;
bool fluid(const Cell& cell)
{const auto p=cell.centerCoord();return !(p.x>1&&p.y<3);}
NativeLeafDomain selected(const QuadtreeMesh& mesh,bool rz,int mode=0,int omit=-1)
{
    FaceSegmentGraph graph;graph.build(mesh,[mode,omit](const QuadtreeMesh&,const Cell& cell){
        const auto p=cell.centerCoord();
        // mode=1 exchanges equal-volume top/bottom right leaves: both XY and
        // RZ totals are identical, while the actual physical mask differs.
        const bool active=mode==1?!(p.x>1&&p.y>3):fluid(cell);
        return active&&cell.index()!=omit;
    });return NativeLeafDomain::build(mesh,graph,rz);
}
std::unique_ptr<QuadtreeMesh> nozzle()
{
    DictionaryMesh d;d.Nx=d.Ny=2;d.minCoord={0,0};d.maxCoord={2,6};
    d.maxCellNum=512;d.maxInitialIterations=0;d.maxLevel_local={2};
    return std::make_unique<QuadtreeMesh>(d);
}
void adapt(QuadtreeMesh& mesh,int level)
{
    CellField<double> aim(mesh,0.),gradient(mesh,0.);
    for(int ci=0;ci<mesh.cellIndexNum();++ci)
        aim.value[ci]=fluid(mesh.cells(mesh.cellIndexToLabel[ci]))?level:0;
    check(mesh.mayAdaptToTargets(aim),"selected actual refine/coarsen preflight");
    std::vector<VertField<double>> levelSets;std::vector<VertField<double>*> vertices;
    std::vector<HEdgeField<double>*> horizontal;std::vector<VEdgeField<double>*> vertical;
    std::vector<CellField<double>*> cells{&aim};const auto old=mesh.cellIndexNum();
    mesh.adapt(levelSets,gradient,vertices,horizontal,vertical,cells,&aim);
    check(old!=mesh.cellIndexNum(),"selected actual tree did not change");++changes;
}
Stocks planeStocks(const QuadtreeMesh& mesh,const NativeLeafDomain& d)
{
    Stocks stocks(mesh.cellIndexNum());
    for(int ci:d.cellIndexByRow()){
        const auto b=bounds(mesh,ci);const long double y=4.234375L;
        stocks[ci]={capacity({b.x0,b.x1,b.y0,std::min(b.y1,y)},d.axisymmetric()),
            capacity({b.x0,b.x1,std::max(b.y0,y),b.y1},d.axisymmetric())};
    }return stocks;
}
long double physicalVolume(const QuadtreeMesh& mesh,const NativeLeafDomain& d)
{long double total=0;for(int ci:d.cellIndexByRow())total+=capacity(bounds(mesh,ci),d.axisymmetric());return total;}
void inspect(const QuadtreeMesh& mesh,const Publication& p)
{
    const bool rz=p.domain().axisymmetric();std::array<long double,2> sum{};
    for(int ci=0;ci<mesh.cellIndexNum();++ci){const auto& c=p.cells()[ci];const auto& q=p.plicCells()[ci];
        if(p.domain().cellRowByIndex()[ci]<0){
            check(c.nativeCell==-1&&c.capacity==0&&c.inventory==0&&c.measured==0&&p.alpha()[ci]==0,
                "selected solid native sentinel contains capacity or stock");
            check(q.inventory[0]==0&&q.inventory[1]==0&&q.measured[0]==0&&q.measured[1]==0&&
                q.rectangle.lower.x==0&&q.rectangle.lower.y==0&&q.rectangle.upper.x==0&&q.rectangle.upper.y==0,
                "selected solid contains PLIC geometry or phase inventory");
            check(p.capacityDefects()[ci].compare(Exact(0))==0,"selected solid contains closure evidence");++solidChecks;continue;
        }
        const auto cap=capacity(bounds(mesh,ci),rz);near(c.capacity,cap,cap,"selected independent active capacity",tolerance);
        near(q.inventory[0]+q.inventory[1],cap,cap,"selected independent active phase closure",tolerance);
        Exact defect(q.inventory[0]);defect.add(q.inventory[1]);defect.add(c.capacity,-1);
        check(defect.compare(p.capacityDefects()[ci])==0,"selected retained exact closure association");
        for(int phase=0;phase<2;++phase){check(q.inventory[phase]>=0,"selected negative stock");sum[phase]+=q.inventory[phase];}
        ++activeChecks;
    }
    check(sum[0]>0&&sum[1]>0,"selected transfer must retain two globally nonzero phases");
    const auto count=p.report().mixedCells+p.report().pureLiquidCells+p.report().pureGasCells;
    check(count==p.domain().cellIndexByRow().size(),"selected publication report counted solid as a phase");
}
void verify(const QuadtreeMesh& sourceMesh,const Publication& source,const QuadtreeMesh& targetMesh,const Transfer& result)
{
    LocalPlicContract::verify(sourceMesh,source,targetMesh,result);inspect(targetMesh,result.candidate());
    const auto& coverage=result.nativeOverlaps();const auto& target=result.candidate();
    const bool rz=source.domain().axisymmetric();std::size_t pairs=0;
    std::vector<long double> sourceCovered(sourceMesh.cellIndexNum()),targetCovered(targetMesh.cellIndexNum());
    // Independent test-owned all-pairs rectangle intersection. Production uses
    // the real two quadtree lookup sweeps and cannot borrow this oracle.
    for(int si=0;si<sourceMesh.cellIndexNum();++si){
        if(source.domain().cellRowByIndex()[si]<0){
            const auto& b=coverage.sourceRectangles()[si];check(coverage.bySource()[si].empty()&&
                b.lower.x==0&&b.lower.y==0&&b.upper.x==0&&b.upper.y==0,"solid source has overlap geometry");continue;
        }
        for(int ti=0;ti<targetMesh.cellIndexNum();++ti){
            if(target.domain().cellRowByIndex()[ti]<0)continue;
            const auto box=intersection(bounds(sourceMesh,si),bounds(targetMesh,ti));const auto expected=capacity(box,rz);
            if(expected==0)continue;
            std::size_t matches=0;
            for(const auto& piece:coverage.pieces())if(piece.sourceCell==si&&piece.targetCell==ti){
                near(piece.capacity,expected,expected,"selected independent positive intersection",tolerance);++matches;}
            check(matches==1,"selected actual overlap lacks unique independent rectangle intersection");
            sourceCovered[si]+=expected;targetCovered[ti]+=expected;++pairs;++coverageChecks;
        }
    }
    check(pairs==coverage.pieces().size(),"selected overlap includes a nonphysical or solid pair");
    const auto checkSide=[&](const QuadtreeMesh& mesh,const NativeLeafDomain& d,const auto& totals,const auto& groups,const auto& boxes){
        for(int ci=0;ci<mesh.cellIndexNum();++ci){
            if(d.cellRowByIndex()[ci]<0){const auto& b=boxes[ci];check(groups[ci].empty()&&totals[ci]==0&&
                b.lower.x==0&&b.lower.y==0&&b.upper.x==0&&b.upper.y==0,"selected solid coverage is not empty/zero");continue;}
            const auto cap=capacity(bounds(mesh,ci),rz),error=std::abs(totals[ci]-cap)/cap;
            peakCoverage=std::max(peakCoverage,error);check(error<=tolerance,"selected independent per-leaf coverage failed");
        }
    };
    checkSide(sourceMesh,source.domain(),sourceCovered,coverage.bySource(),coverage.sourceRectangles());
    checkSide(targetMesh,target.domain(),targetCovered,coverage.byTarget(),coverage.targetRectangles());
    const auto& r=result.report();check(r.maximumRelativeCoverageResidual<=tolerance&&
        r.maximumRelativeTargetPhaseClosure<=tolerance&&r.maximumRelativeGeometryResidual<=tolerance,
        "selected coverage/closure/geometry report gate");
    check(r.unchangedTargets+r.refinedTargets+r.coarsenedTargets==target.domain().cellIndexByRow().size(),
        "selected target transfer counters include inactive cells");
}
template<class F> void rejectOverlap(F&& f,const std::string& name)
{
    bool thrown=false;try{f();}catch(const std::invalid_argument& error){
        check(std::string(error.what()).find("[NativeLeafOverlap]")!=std::string::npos,
            "changed physical mask rejected outside actual overlap: "+std::string(error.what()));thrown=true;}
    check(thrown,"actual overlap accepted changed physical mask: "+name);++negativeChecks;
}
void run()
{
    for(bool rz:{false,true}){
        stage=rz?"selected nozzle L3 RZ actual repeated AMR":"selected nozzle L3 XY actual repeated AMR";
        auto root=nozzle();const auto d=selected(*root,rz);
        check(root->cellIndexNum()==4&&d.cellIndexByRow().size()==3,"actual nozzle whole-leaf L3 fixture");
        const auto stocks=planeStocks(*root,d);
        const auto source=Publication::prepareLocalPlicInventory(*root,d,accepted,stocks,{});
        const auto saved=LocalPlicContract::retained(source);const auto id=source.identity();inspect(*root,source);
        auto mesh=root->clone();auto unchanged=Transfer::prepare(*root,source,*mesh,selected(*mesh,rz),accepted);
        verify(*root,source,*mesh,unchanged);
        check(LocalPlicContract::retained(unchanged.candidate())==saved,"selected unchanged stocks not exact");
        Publication current=unchanged.candidate();unsigned up=0,down=0;bool sawGraded=false;
        for(int level:{1,0,1,2,1,0}){
            auto next=mesh->clone();adapt(*next,level);const auto nextDomain=selected(*next,rz);
            sawGraded|=!nextDomain.uniformLeaves();
            const auto oldStock=LocalPlicContract::retained(current);const auto oldId=current.identity();
            const auto result=Transfer::prepare(*mesh,current,*next,nextDomain,accepted);verify(*mesh,current,*next,result);
            up+=result.report().refinedTargets>0;down+=result.report().coarsenedTargets>0;
            check(current.identity()==oldId&&LocalPlicContract::retained(current)==oldStock,
                "selected successful candidate changed predecessor");
            current=result.candidate();mesh=std::move(next);
        }
        check(up>=3&&down>=3&&sawGraded,"selected AMR requires repeated real refinement/coarsening and graded fluid/solid tree");
        for(int phase=0;phase<2;++phase){long double first=0,last=0;
            for(const auto& q:source.plicCells())first+=q.inventory[phase];for(const auto& q:current.plicCells())last+=q.inventory[phase];
            LocalPlicContract::ownNear(last,first,"selected multiround independent phase conservation");}
        stage=rz?"selected RZ same-volume physical-mask rejection":"selected XY same-volume physical-mask rejection";
        const auto wrong=selected(*root,rz,1);
        check(wrong.cellRowByIndex()!=d.cellRowByIndex(),"equal-volume mask witness unchanged");
        near(physicalVolume(*root,wrong),physicalVolume(*root,d),physicalVolume(*root,d),"changed mask has equal physical volume",tolerance);
        rejectOverlap([&]{Transfer::prepare(*root,source,*root,wrong,accepted);},"same-size equal-volume solid exchange");
        rejectOverlap([&]{LAMP::Mesh::NativeLeafOverlap::build(*root,d,*root,wrong);},"direct public Core overlap equal-volume exchange");
        rejectOverlap([&]{Transfer::prepare(*root,source,*root,domain(*root,rz),accepted);},"selected to enlarged full domain");
        const auto full=Publication::prepareLocalPlicInventory(*root,domain(*root,rz),accepted,
            planeStocks(*root,domain(*root,rz)),{});
        rejectOverlap([&]{Transfer::prepare(*root,full,*root,d,accepted);},"full to selected physical stock loss");
        auto fine=root->clone();adapt(*fine,1);const auto fd=selected(*fine,rz);
        const auto missing=selected(*fine,rz,0,fd.cellIndexByRow().front());
        rejectOverlap([&]{Transfer::prepare(*root,source,*fine,missing,accepted);},"one missing fine fluid child");
        const auto missingSource=Publication::prepareLocalPlicInventory(*fine,missing,accepted,planeStocks(*fine,missing),{});
        rejectOverlap([&]{Transfer::prepare(*fine,missingSource,*root,d,accepted);},"coarse target partially covered by selected fine leaves");
        Clock clock=accepted;clock.end=.26;rejects([&]{Transfer::prepare(*root,source,*fine,fd,clock);},"selected End full-clock mismatch");
        clock=accepted;clock.position=Position::Midpoint;
        rejects([&]{Transfer::prepare(*root,source,*fine,fd,clock);},"selected transfer at Midpoint");
        rejects([&]{Transfer::prepare(*root,source,*fine,selected(*fine,!rz),accepted);},"selected physical measure mismatch");
        auto invalidStocks=saved;int solid=-1;for(int ci=0;ci<root->cellIndexNum();++ci)if(d.cellRowByIndex()[ci]<0)solid=ci;
        check(solid>=0,"selected solid negative-control witness missing");invalidStocks[solid][1]=1;
        rejects([&]{Publication::prepareLocalPlicInventory(*root,d,accepted,invalidStocks,{});},"raw solid phase stock forbidden");
        auto stale=root->clone();const auto staleDomain=selected(*stale,rz);adapt(*stale,1);
        rejects([&]{Transfer::prepare(*root,source,*stale,staleDomain,accepted);},"selected stale target domain");
        rejects([&]{Transfer::prepare(*root,source,*fine,d,accepted);},"selected target foreign mesh owner");
        check(source.identity()==id&&LocalPlicContract::retained(source)==saved,"selected rejection changed accepted source");
        source.validateFor(*root,id,accepted);
        stage=rz?"selected RZ independent unresolved microphase":"selected XY independent unresolved microphase";
        Stocks micro(root->cellIndexNum());for(int ci:d.cellIndexByRow()){
            const auto cap=capacity(bounds(*root,ci),rz);micro[ci]={cap,cap*1e-40L};}
        const auto tiny=Publication::prepareLocalPlicInventory(*root,d,accepted,micro,{});
        const auto tinySaved=LocalPlicContract::retained(tiny);
        const auto tinyUp=Transfer::prepare(*root,tiny,*fine,fd,accepted);verify(*root,tiny,*fine,tinyUp);
        check(tinyUp.report().unresolvedRefinedDonors==3&&tinyUp.report().p0RefinedPieces==12,
            "selected explicit P0 microphase refinement did not cover active donors only");
        for(int ci:fd.cellIndexByRow())check(tinyUp.candidate().plicCells()[ci].inventory[1]>0,
            "selected unresolved positive microphase erased");
        const auto tinyDown=Transfer::prepare(*fine,tinyUp.candidate(),*root,d,accepted);verify(*fine,tinyUp.candidate(),*root,tinyDown);
        check(LocalPlicContract::retained(tinyDown.candidate())==tinySaved,"selected exact P0 microphase round trip changed stocks");
        std::cout<<"SELECTED_DOMAIN_AMR_METRIC metric="<<(rz?"RZ":"XY")<<" root_fluid_leaves=3 refinements="<<up
            <<" coarsenings="<<down<<" final_fluid_leaves="<<current.domain().cellIndexByRow().size()<<'\n';
    }
    check(changes>=12&&solidChecks>0&&activeChecks>0&&coverageChecks>0&&negativeChecks==12,"selected bounded AMR coverage incomplete");
    std::cout<<std::setprecision(17)<<"SELECTED_DOMAIN_AMR_PASS actual_adaptations="<<changes<<" active_checks="<<activeChecks
        <<" solid_checks="<<solidChecks<<" independent_overlaps="<<coverageChecks<<" physical_mask_guards="<<negativeChecks
        <<" peak_coverage="<<peakCoverage<<" peak_own_phase="<<LocalPlicContract::peakOwn
        <<" peak_joint="<<LocalPlicContract::peakJoint<<'\n';
}
}


namespace NativeLeafOverlapContract
{
namespace Cut=LAMP::Geometry::PlanarCellCut;
using Coverage=LAMP::Mesh::NativeLeafOverlap;
std::size_t pieces=0,capacityChecks=0;
long double peakNativeVolume=0;
void adaptBand(QuadtreeMesh& mesh,bool refine)
{
    CellField<double> aim(mesh,0.),gradient(mesh,0.);
    for(int ci=0;ci<mesh.cellIndexNum();++ci){const auto p=mesh.cells(mesh.cellIndexToLabel[ci]).centerCoord();
        aim.value[ci]=refine&&p.y>3.6&&p.y<4.6?1:0;}
    check(mesh.mayAdaptToTargets(aim),"native overlap actual band tree preflight");
    std::vector<VertField<double>> levelSets;std::vector<VertField<double>*> vertices;
    std::vector<HEdgeField<double>*> horizontal;std::vector<VEdgeField<double>*> vertical;
    std::vector<CellField<double>*> cells{&aim};const auto before=mesh.cellIndexNum();
    mesh.adapt(levelSets,gradient,vertices,horizontal,vertical,cells,&aim);
    check(mesh.cellIndexNum()!=before,"native overlap actual adaptation missing");
}
void support(const QuadtreeMesh& mesh,const Publication& publication,
    const std::vector<Cut::Rectangle>& rectangles)
{
    const auto& d=publication.domain();
    check(rectangles.size()==static_cast<std::size_t>(mesh.cellIndexNum()),"native overlap rectangle map dimensions");
    long double volume=0;
    for(int ci=0;ci<mesh.cellIndexNum();++ci){const auto& box=rectangles[ci];const auto& q=publication.plicCells()[ci];
        if(d.cellRowByIndex()[ci]<0){check(box.lower.x==0&&box.lower.y==0&&box.upper.x==0&&box.upper.y==0
            &&q.inventory[0]==0&&q.inventory[1]==0,"native overlap solid acquired geometry/inventory");continue;}
        const auto& cell=mesh.cells(mesh.cellIndexToLabel[ci]);
        const auto a=mesh.verts(cell.cellToVert[0]).coord(),b=mesh.verts(cell.cellToVert[2]).coord();
        const auto address=LAMP::Mesh::cellAddress(mesh,cell);
        const LAMP::Mesh::DyadicGridCoordinates coordinates(mesh.minCoord,mesh.maxCoord,mesh.Nx,mesh.Ny);
        const auto lo=coordinates.point(address.i,address.j,address.level),hi=coordinates.point(address.i+1,address.j+1,address.level);
        check(box.lower.x==a.x&&box.lower.y==a.y&&box.upper.x==b.x&&box.upper.y==b.y,
            "Core overlap support differs from actual native vertices");
        check(a.x==lo.x&&a.y==lo.y&&b.x==hi.x&&b.y==hi.y,"native overlap vertices differ from canonical ancestry");
        check(q.rectangle.lower.x==box.lower.x&&q.rectangle.lower.y==box.lower.y
            &&q.rectangle.upper.x==box.upper.x&&q.rectangle.upper.y==box.upper.y,
            "Core overlap and Interface phase support differ");
        const auto full=Cut::capacity(box,true),independent=capacity(bounds(mesh,ci),true);
        check(full==publication.cells()[ci].capacity,"Core overlap/Publication exact capacity association");
        near(full,independent,full,"native overlap independently integrated RZ capacity",tolerance);
        const auto nativeVolume=d.cellVolumes()[d.cellRowByIndex()[ci]];
        peakNativeVolume=std::max(peakNativeVolume,std::abs(full-nativeVolume)/full);
        check(std::abs(full-nativeVolume)<=tolerance*full,"native overlap native cell-volume gate");
        volume+=full;++capacityChecks;
    }
    near(volume,328*std::acos(-1.L),volume,"native overlap selected nozzle total RZ volume",tolerance);
    for(const auto f:d.faces(mesh))if(f.negativeCell>=0&&f.positiveCell>=0)
        check(std::abs(mesh.cells(mesh.cellIndexToLabel[f.negativeCell]).level()
            -mesh.cells(mesh.cellIndexToLabel[f.positiveCell]).level())<=1,"native overlap target lost actual 2:1 balance");
}
void inspect(const QuadtreeMesh& sourceMesh,const Publication& source,
    const QuadtreeMesh& targetMesh,const Transfer& transfer)
{
    const auto& target=transfer.candidate();const auto& coverage=transfer.nativeOverlaps();
    coverage.validateFor(sourceMesh,targetMesh);support(sourceMesh,source,coverage.sourceRectangles());
    support(targetMesh,target,coverage.targetRectangles());
    check(source.timeLayer().position==Position::End&&sameClock(source.timeLayer(),accepted)
        &&sameClock(target.timeLayer(),accepted),"native overlap lost actual complete accepted-End clock");
    check(transfer.sourceIdentity()==source.identity()&&target.identity()!=source.identity(),"native overlap source/fresh target owner");
    check(coverage.maximumRelativeCoverageResidual()<=512*std::numeric_limits<long double>::epsilon(),
        "native overlap original strict physical coverage gate");
    std::array<long double,2> before{},after{};
    for(const auto& q:source.plicCells())for(int phase=0;phase<2;++phase)before[phase]+=q.inventory[phase];
    for(const auto& q:target.plicCells())for(int phase=0;phase<2;++phase)after[phase]+=q.inventory[phase];
    for(int phase=0;phase<2;++phase){check(before[phase]>0,"native overlap source missing actual phase");
        near(after[phase],before[phase],before[phase],"native overlap independent phase conservation",tolerance);}
    near(before[0],28*std::acos(-1.L),before[0],"native overlap actual liquid source volume",tolerance);
    near(before[1],300*std::acos(-1.L),before[1],"native overlap actual gas source volume",tolerance);
    check(coverage.pieces().size()==transfer.overlaps().size(),"native overlap Core/Interface piece dimensions");
    for(std::size_t i=0;i<coverage.pieces().size();++i){const auto& piece=coverage.pieces()[i];const auto& o=transfer.overlaps()[i];
        const auto& box=piece.kind==LAMP::Mesh::NativeLeafOverlapPiece::Kind::Refined?
            coverage.targetRectangles()[piece.targetCell]:coverage.sourceRectangles()[piece.sourceCell];
        check(piece.capacity==Cut::capacity(box,true)&&o.capacity==piece.capacity
            &&o.sourceCell==piece.sourceCell&&o.targetCell==piece.targetCell,"native overlap actual piece capacity/binding");++pieces;}
    source.validateFor(sourceMesh,source.identity(),accepted);target.validateFor(targetMesh,target.identity(),accepted);
    rejects([&]{target.validateFor(targetMesh,source.identity(),accepted);},"native overlap foreign Publication owner");
    auto wrongClock=accepted;wrongClock.begin=.126;
    rejects([&]{Transfer::prepare(sourceMesh,source,targetMesh,target.domain(),wrongClock);},"native overlap changed accepted-End clock");
}
void run()
{
    stage="native overlap original CSG base40x80 RZ actual sameEnd band refinement/coarsening";
    DictionaryMesh d;d.Nx=40;d.Ny=80;d.minCoord={0,0};d.maxCoord={5,16};
    d.levelSetCount=0;d.maxCellNum=16384;d.maxInitialIterations=0;d.maxLevel_local={1};
    auto root=std::make_unique<QuadtreeMesh>(d);const NativeMeshIdentity meshIdentity(*root);
    const auto selected=SelectedDomainAmrContract::selected(*root,true);
    check(root->cellIndexNum()==3200&&selected.cellIndexByRow().size()==2720,"native overlap original nozzle base counts");
    const auto& witness=root->cells(root->cellIndexToLabel[800]);const auto vertex=root->verts(witness.cellToVert[0]).coord();
    check(vertex.y==4&&witness.centerCoord().y-witness.size().y/2>vertex.y,
        "native overlap exact old center-minus-half failure witness missing");
    std::vector<double> alpha(root->cellIndexNum());
    for(int ci:selected.cellIndexByRow())alpha[ci]=root->cells(root->cellIndexToLabel[ci]).centerCoord().y<4?1:0;
    const auto source=Publication::prepareLocalPlic(*root,selected,accepted,alpha,{});
    const auto sourceId=source.identity();const auto saved=LocalPlicContract::retained(source);const auto savedAlpha=source.alpha();
    auto fine=root->clone();adaptBand(*fine,true);const auto fineDomain=SelectedDomainAmrContract::selected(*fine,true);
    check(!fineDomain.uniformLeaves(),"native overlap band did not create real graded leaves");
    const auto up=Transfer::prepare(*root,source,*fine,fineDomain,accepted);inspect(*root,source,*fine,up);
    check(up.report().refinedTargets>0&&up.nativeOverlaps().bySource()[800].size()==4,"native overlap failed actual named parent refinement");
    auto back=fine->clone();adaptBand(*back,false);const auto backDomain=SelectedDomainAmrContract::selected(*back,true);
    const auto down=Transfer::prepare(*fine,up.candidate(),*back,backDomain,accepted);inspect(*fine,up.candidate(),*back,down);
    check(back->cellIndexNum()==3200&&down.report().coarsenedTargets>0,"native overlap failed actual coarsening to original base");
    for(std::size_t ci=0;ci<saved.size();++ci)for(int phase=0;phase<2;++phase)
        near(down.candidate().plicCells()[ci].inventory[phase],saved[ci][phase],source.cells()[ci].capacity>0?source.cells()[ci].capacity:1,
            "native overlap local phase roundtrip",tolerance);
    check(source.identity()==sourceId&&LocalPlicContract::retained(source)==saved&&source.alpha()==savedAlpha,
        "native overlap candidate changed accepted source owner/stocks/alpha");meshIdentity.validate(*root);
    std::cout<<std::setprecision(17)<<"NATIVE_LEAF_OVERLAP_PASS base=40x80 native=3200 active=2720 solid=480 refinements=1 coarsenings=1 pieces="
        <<pieces<<" capacity_checks="<<capacityChecks<<" peak_native_volume="<<peakNativeVolume<<'\n';
}
}

}
int main(int argc,char** argv)
{
    try{
        if(argc==2&&std::string(argv[1])=="--native-leaf-overlap-only"){NativeLeafOverlapContract::run();return 0;}
        if(argc==2&&std::string(argv[1])=="--selected-domain-amr-only"){SelectedDomainAmrContract::run();return 0;}
        if(argc==2&&std::string(argv[1])=="--local-plic-capacity-partition-only"){LocalPlicContract::capacityCycles();return 0;}
        if(argc==2&&std::string(argv[1])=="--local-plic-only"){LocalPlicContract::group();return 0;}
        check(argc==1,"unknown contract argument");
        staticTransfers();repeatedAdaptation();pureAndTiny();smallPhysicalDomain();syntheticTopology();rejectionTransaction();
        LocalPlicContract::group();
        std::cout<<std::setprecision(17)<<"GEOMETRIC_VOF_MESH_TRANSFER_PASS scenarios="<<scenarios<<" assertions="<<assertions
            <<" overlaps="<<overlapChecks<<" actual_adaptations="<<adaptations<<" guards="<<guards
            <<" peak_inventory="<<peakInventory<<" peak_gas="<<peakGas<<" peak_donor_closure="<<peakDonorClosure
            <<" peak_cell_closure="<<peakCellClosure<<" peak_minor_factor="<<peakMinorFactor<<" peak_raw="<<peakRaw<<" peak_geometry="<<peakGeometry<<" peak_green="<<peakGreen<<'\n';
        return 0;
    }catch(const std::exception& e){std::cerr<<"GEOMETRIC_VOF_MESH_TRANSFER_FAIL stage="<<stage<<" reason="<<e.what()<<'\n';return 1;}
}




