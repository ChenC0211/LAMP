#include "Modules/Hydrodynamics/CchMomentumGeometry.h" // Public header stands alone.
#include "Modules/Hydrodynamics/CchSpatialAdvector.h"
#include "Modules/Hydrodynamics/Detail/CchRegularAdvectorGeometry.h"
#include "Modules/Hydrodynamics/CchSinglePhaseSaddle.h"
#include "Modules/Hydrodynamics/CchPhaseViscosity.h"
#include "Modules/Hydrodynamics/CchInterfaceNormalLoad.h"
#include "Modules/Hydrodynamics/CchMaterialTransport.h"
#include "Modules/Hydrodynamics/CchDualMaterial.h"
#include "Modules/Hydrodynamics/CchKepConvection.h"
#include "interface_first_jet_fixture.h"
#include "Core/Mesh/CartesianMeasure.h"
#include <map>
#include <tuple>
#include <Eigen/SVD>
#include <Eigen/LU>
#include <cfenv>
#include <chrono>
#include <xmmintrin.h>

namespace
{
using namespace Quantum::Hydrodynamics;
int geometryScenarios=0,massScenarios=0,overlapChecks=0,sharedPairs=0,segmentChecks=0,rejected=0,spaceRows=0,rebuilds=0;
std::array<int,2> amrCutDuals{};
long double peakGeometry=0,peakPhase=0,peakMass=0,peakPartition=0,averageGap=0;
static_assert(!std::is_default_constructible_v<CchMomentumGeometry>);
static_assert(!std::is_default_constructible_v<CchMomentumMass>);
static_assert(std::is_same_v<decltype(std::declval<CchMomentumGeometry&>().volumes()),const std::vector<CchMomentumVolume>&>);
static_assert(std::is_same_v<decltype(std::declval<CchMomentumMass&>().masses()),const std::vector<double>&>);
long double cap(Vector2D a,Vector2D b,bool rz)
{ return (static_cast<long double>(b.x)-a.x)*(static_cast<long double>(b.y)-a.y)*(rz?pi*(static_cast<long double>(a.x)+b.x):1.L); }
long double area(Vector2D a,Vector2D b,int normal,bool rz)
{ return (normal==0?static_cast<long double>(b.y)-a.y:static_cast<long double>(b.x)-a.x)*(rz?pi*(static_cast<long double>(a.x)+b.x):1.L); }
double component(Vector2D p,int axis) { return axis==0?p.x:p.y; }
int axis(const CchDofLayout& l,int u) { return l.normalAxis(u)==CchDofLayout::NormalAxis::First?0:1; }
bool pointEqual(Vector2D a,Vector2D b) { return a.x==b.x && a.y==b.y; }
long double norm(long double a,long double b,long double scale,long double limit,const std::string& why)
{
    const auto e=std::abs(a-b)/(scale>0?scale:1);
    check(std::isfinite(e) && e<=limit,why+" error="+std::to_string(static_cast<double>(e)));return e;
}
CchPhysicalDivergence metricsFor(const QuadtreeMesh& mesh,bool rz,int selection=0)
{
    FaceSegmentGraph graph;graph.build(mesh,[selection](const QuadtreeMesh&,const Cell& c) {
        const auto p=c.centerCoord();
        if(selection==1)return p.x<.5;
        if(selection==2)return !(p.x>.25 && p.x<.75 && p.y>.25 && p.y<.75);
        if(selection==3)return p.x<.25 || p.x>.75;
        return true;
    });
    return CchPhysicalDivergence::build(mesh,CchIntegralDivergence::build(mesh,CchDofLayout::build(mesh,NativeLeafDomain::build(mesh,graph,rz))));
}
// Independent ownership enumeration starts from a native face and its actual
// neighboring cells, not the production cell-side builder or reported overlaps.
struct ReferenceBox { Vector2D lower{0,0},upper{0,0}; std::vector<CchMomentumOverlap> pieces; };
ReferenceBox reference(const QuadtreeMesh& mesh,const CchDofLayout& layout,int u)
{
    ReferenceBox box;const auto face=layout.domain().face(mesh,layout.nativeFace(u));const int a=axis(layout,u);
    const auto& edge=face.horizontal?mesh.hEdges(face.edgeLabel):mesh.vEdges(face.edgeLabel);
    auto start=mesh.verts(edge.edgeToVert[0]).coord(),end=mesh.verts(edge.edgeToVert[1]).coord();
    if(component(start,1-a)>component(end,1-a))std::swap(start,end);
    for(int native:{face.negativeCell,face.positiveCell})if(native>=0) {
        const auto& cell=mesh.cells(mesh.cellIndexToLabel[native]);
        auto lo=mesh.verts(cell.cellToVert[0]).coord(),hi=mesh.verts(cell.cellToVert[2]).coord();
        const auto whole=cap(lo,hi,layout.domain().axisymmetric());const auto center=lo+(hi-lo)*.5;
        if(a==0) {lo.y=start.y;hi.y=end.y;if(native==face.negativeCell)lo.x=center.x;else hi.x=center.x;}
        else {lo.x=start.x;hi.x=end.x;if(native==face.negativeCell)lo.y=center.y;else hi.y=center.y;}
        const auto v=cap(lo,hi,layout.domain().axisymmetric());
        if(box.pieces.empty()){box.lower=lo;box.upper=hi;}
        else{box.lower={std::min(box.lower.x,lo.x),std::min(box.lower.y,lo.y)};box.upper={std::max(box.upper.x,hi.x),std::max(box.upper.y,hi.y)};}
        box.pieces.push_back({layout.pressureDof(native),lo,hi,v,v/whole});
    }
    return box;
}
void inspectGeometry(const QuadtreeMesh& mesh,const CchMomentumGeometry& geometry)
{
    geometry.validateMesh(mesh);const auto& l=geometry.layout();const bool rz=l.domain().axisymmetric();
    std::vector<ReferenceBox> boxes;std::array<std::vector<long double>,2> cellVolumes;
    for(auto& row:cellVolumes)row.resize(l.pressureCount());
    check(geometry.volumes().size()==static_cast<std::size_t>(l.velocityCount()),"new global velocity count");
    for(int u=0;u<l.velocityCount();++u) {
        boxes.push_back(reference(mesh,l,u));const auto& expected=boxes.back();const auto& actual=geometry.volumes()[u];
        check(pointEqual(expected.lower,actual.lower) && pointEqual(expected.upper,actual.upper),"dual uses virtual-fine centers");
        check(expected.pieces.size()==actual.overlaps.size(),"wrong dual owner count");
        peakGeometry=std::max(peakGeometry,norm(actual.volume,cap(expected.lower,expected.upper,rz),actual.volume,1e-12,"independent dual capacity"));
        for(const auto& ref:expected.pieces) {
            const auto found=std::find_if(actual.overlaps.begin(),actual.overlaps.end(),[&](const auto& p){return p.pressure==ref.pressure;});
            check(found!=actual.overlaps.end() && pointEqual(found->lower,ref.lower) && pointEqual(found->upper,ref.upper),"dual half-cell overlap differs");
            peakGeometry=std::max({peakGeometry,norm(found->volume,ref.volume,ref.volume,1e-12,"independent overlap capacity"),norm(found->cellFraction,ref.cellFraction,ref.cellFraction,1e-12,"true geometric transfer weight")});
            cellVolumes[axis(l,u)][ref.pressure]+=found->volume;++overlapChecks;
        }
    }
    for(const auto& row:cellVolumes)for(int p=0;p<l.pressureCount();++p)
        peakGeometry=std::max(peakGeometry,norm(row[p],geometry.metrics().cellVolumes()[p],row[p],1e-12,"independent cell partition"));
    // Each common side must equal the independent intersection of two dual
    // rectangles. Production can subdivide it at actual primal boundaries.
    using Key=std::tuple<int,int,int>;
    std::map<Key,std::vector<int>> byPair;
    std::vector<int> uses(geometry.segments().size());
    for(int u=0;u<l.velocityCount();++u)for(int side=0;side<4;++side) {
        std::vector<std::pair<double,double>> intervals;long double sum=0;
        for(int id:geometry.volumes()[u].sides[side]) {
            check(id>=0 && static_cast<std::size_t>(id)<uses.size(),"invalid side segment index");++uses[id];
            const auto& s=geometry.segments()[id];const int normal=side/2;
            check(s.normalAxis==normal && (side%2==0?s.positiveVelocity:s.negativeVelocity)==u,"incorrect side incidence");
            check(component(s.begin,normal)==component(side%2==0?boxes[u].lower:boxes[u].upper,normal) && component(s.end,normal)==component(s.begin,normal),"misplaced dual segment");
            intervals.emplace_back(component(s.begin,1-normal),component(s.end,1-normal));sum+=s.area;
        }
        std::sort(intervals.begin(),intervals.end());double cursor=component(boxes[u].lower,1-side/2);
        for(const auto& interval:intervals){check(interval.first==cursor && interval.second>interval.first,"side coverage gap/overlap");cursor=interval.second;}
        check(cursor==component(boxes[u].upper,1-side/2),"incomplete side coverage");
        Vector2D a=boxes[u].lower,b=boxes[u].upper;
        if(side/2==0)a.x=b.x=component(side%2==0?boxes[u].lower:boxes[u].upper,0);
        else a.y=b.y=component(side%2==0?boxes[u].lower:boxes[u].upper,1);
        peakGeometry=std::max(peakGeometry,norm(sum,area(a,b,side/2,rz),geometry.volumes()[u].volume,1e-12,"independent side area"));
    }
    for(std::size_t id=0;id<uses.size();++id) {
        const auto& s=geometry.segments()[id];const int expected=(s.negativeVelocity>=0?1:0)+(s.positiveVelocity>=0?1:0);
        check(uses[id]==expected,"shared segment not stored reciprocally");
        peakGeometry=std::max(peakGeometry,norm(s.area,area(s.begin,s.end,s.normalAxis,rz),1,1e-12,"independent segment area"));
        check((s.primalVelocity>=0)!=(s.ownerPressure>=0),"ambiguous segment source");
        if(s.primalVelocity>=0){
            const auto primal=l.domain().face(mesh,l.nativeFace(s.primalVelocity));const auto& e=primal.horizontal?mesh.hEdges(primal.edgeLabel):mesh.vEdges(primal.edgeLabel);
            const auto a=mesh.verts(e.edgeToVert[0]).coord(),b=mesh.verts(e.edgeToVert[1]).coord();
            check(s.normalAxis==(primal.horizontal?1:0) && component(s.begin,s.normalAxis)==component(primal.center,s.normalAxis),"wrong primal source plane");
            check(component(s.begin,1-s.normalAxis)>=std::min(component(a,1-s.normalAxis),component(b,1-s.normalAxis)) && component(s.end,1-s.normalAxis)<=std::max(component(a,1-s.normalAxis),component(b,1-s.normalAxis)),"segment leaves primal source");
        }else{
            const auto& c=mesh.cells(mesh.cellIndexToLabel[l.nativeCell(s.ownerPressure)]);
            const auto a=mesh.verts(c.cellToVert[0]).coord(),b=mesh.verts(c.cellToVert[2]).coord(),m=a+(b-a)*.5;
            check(component(s.begin,s.normalAxis)==component(m,s.normalAxis) && component(s.begin,1-s.normalAxis)>=component(a,1-s.normalAxis) && component(s.end,1-s.normalAxis)<=component(b,1-s.normalAxis),"wrong interior cell source");
        }
        if(expected==2){check(axis(l,s.negativeVelocity)==axis(l,s.positiveVelocity),"mixed momentum components");byPair[{s.negativeVelocity,s.positiveVelocity,s.normalAxis}].push_back(static_cast<int>(id));}
        if(rz && s.normalAxis==0 && s.begin.x==0)check(s.area==0,"nonzero axis area");++segmentChecks;
    }
    std::size_t pairCount=0;
    for(int u=0;u<l.velocityCount();++u)for(int v=u+1;v<l.velocityCount();++v)if(axis(l,u)==axis(l,v))for(int normal=0;normal<2;++normal) {
        int negative=u,positive=v;
        if(component(boxes[v].upper,normal)==component(boxes[u].lower,normal))std::swap(negative,positive);
        else if(component(boxes[u].upper,normal)!=component(boxes[v].lower,normal))continue;
        const double lo=std::max(component(boxes[u].lower,1-normal),component(boxes[v].lower,1-normal));
        const double hi=std::min(component(boxes[u].upper,1-normal),component(boxes[v].upper,1-normal));if(lo>=hi)continue;
        const auto found=byPair.find({negative,positive,normal});check(found!=byPair.end(),"missing independent shared rectangle side");
        auto ids=found->second;std::sort(ids.begin(),ids.end(),[&](int a,int b){return component(geometry.segments()[a].begin,1-normal)<component(geometry.segments()[b].begin,1-normal);});
        double cursor=lo;for(int id:ids){const auto& s=geometry.segments()[id];check(component(s.begin,1-normal)==cursor,"shared pair gap/overlap");cursor=component(s.end,1-normal);}check(cursor==hi,"shared pair coverage differs");++pairCount;++sharedPairs;
    }
    check(pairCount==byPair.size(),"unexpected shared rectangle side");++geometryScenarios;
}
long double liquid(Vector2D a,Vector2D b,const Case& s)
{
    if(!s.planar)return circleCell(a,b,{s.radial?0:s.cx,s.cy},s.radius,s.radial);
    if(s.axial()){if(s.planeSign>0)b.y=std::min(b.y,s.planeOffset);else a.y=std::max(a.y,s.planeOffset);}
    else{if(s.planeSign>0)b.x=std::min(b.x,s.planeOffset);else a.x=std::max(a.x,s.planeOffset);}
    return a.x>=b.x || a.y>=b.y?0:cap(a,b,s.radial);
}
long double inspectMass(const QuadtreeMesh& mesh,const Case& s,const CchMomentumMass& mass,bool exact=true)
{
    const auto& l=mass.layout();mass.validateFor(mesh,mass.geometryIdentity(),mass.timeLayer());
    check(mass.phaseGeometry().publication().state().primary().phaseVolumes.empty(),"invented material inventory");
    std::array<long double,2> totalMass{},volumeSum{},phaseSum{};
    std::array<std::vector<CchPhaseMeasure>,2> cellSums;for(auto& a:cellSums)a.resize(l.pressureCount());long double error=0;
    for(int u=0;u<l.velocityCount();++u) {
        const auto box=reference(mesh,l,u);const auto volume=cap(box.lower,box.upper,s.radial),vl=liquid(box.lower,box.upper,s);
        const long double expected=mass.liquidDensity()*vl+mass.gasDensity()*(volume-vl);
        const auto actual=mass.phaseVolumes()[u];const auto scale=(mass.liquidDensity()+mass.gasDensity())*volume;
        if(exact){peakPhase=std::max({peakPhase,norm(actual.liquid,vl,volume,1e-10,"analytic dual liquid"),norm(actual.gas,volume-vl,volume,1e-10,"analytic dual gas")});peakMass=std::max(peakMass,norm(mass.masses()[u],expected,scale,1e-10,"analytic dual mass"));}
        error+=std::abs(mass.masses()[u]-expected);const int a=axis(l,u);
        const auto face=l.domain().face(mesh,l.nativeFace(u));
        if(face.negativeCell>=0 && face.positiveCell>=0 && actual.liquid>0 && actual.gas>0) {
            const auto& negative=mesh.cells(mesh.cellIndexToLabel[face.negativeCell]);
            const auto& positive=mesh.cells(mesh.cellIndexToLabel[face.positiveCell]);
            if(negative.size().x!=positive.size().x)++amrCutDuals[s.radial?1:0];
        }
        totalMass[a]+=mass.masses()[u];volumeSum[a]+=volume;phaseSum[a]+=actual.liquid;
        long double average=0,overlapLiquid=0,overlapGas=0;
        for(std::size_t i=0;i<mass.geometry().volumes()[u].overlaps.size();++i) {
            const auto& p=mass.geometry().volumes()[u].overlaps[i];const auto phases=mass.overlapPhaseVolumes()[u][i];
            const auto exactLiquid=liquid(p.lower,p.upper,s);
            if(exact)peakPhase=std::max({peakPhase,norm(phases.liquid,exactLiquid,p.volume,1e-10,"analytic native overlap liquid"),norm(phases.gas,p.volume-exactLiquid,p.volume,1e-10,"analytic native overlap gas")});
            overlapLiquid+=phases.liquid;overlapGas+=phases.gas;
            auto& row=cellSums[a][p.pressure];row.liquid+=phases.liquid;row.gas+=phases.gas;
            const auto& cell=mass.phaseGeometry().cellVolumes()[p.pressure];average+=p.cellFraction*(mass.liquidDensity()*cell.liquid+mass.gasDensity()*cell.gas);
        }
        peakPartition=std::max({peakPartition,norm(overlapLiquid,actual.liquid,volume,1e-10,"overlap liquid assembly"),norm(overlapGas,actual.gas,volume,1e-10,"overlap gas assembly")});
        if(s.planar && mass.liquidDensity()==1e6 && s.planeOffset==.43)averageGap=std::max(averageGap,std::abs(average-mass.masses()[u])/scale);
    }
    long double cellLiquid=0,cellMass=0;
    for(int p=0;p<l.pressureCount();++p) {
        const auto& ref=mass.phaseGeometry().cellVolumes()[p];cellLiquid+=ref.liquid;cellMass+=mass.liquidDensity()*ref.liquid+mass.gasDensity()*ref.gas;
        for(const auto& row:cellSums)peakPartition=std::max({peakPartition,norm(row[p].liquid,ref.liquid,mass.geometry().metrics().cellVolumes()[p],1e-10,"cell liquid partition"),norm(row[p].gas,ref.gas,mass.geometry().metrics().cellVolumes()[p],1e-10,"cell gas partition")});
    }
    for(int a=0;a<2;++a){peakPartition=std::max({peakPartition,norm(totalMass[a],cellMass,cellMass,1e-10,"component mass partition"),norm(phaseSum[a],cellLiquid,volumeSum[a],1e-10,"component liquid partition")});norm(mass.report().totalMass[a],totalMass[a],cellMass,1e-12,"reported mass total");}
    ++massScenarios;return error/((mass.liquidDensity()+mass.gasDensity())*(volumeSum[0]+volumeSum[1]));
}
CchMomentumMass buildMass(const QuadtreeMesh& mesh,const CchMomentumGeometry& geometry,const Case& s,double rho=1000)
{ return CchMomentumMass::build(mesh,geometry,CchPhaseGeometry::build(mesh,geometry.metrics(),initial(mesh,s)),rho,1); }
void convergence(bool rz,int amr)
{
    double previous=0,order=0;
    for(int n:{4,8,16}) {
        auto mesh=makeMesh(n,amr);Case s;s.pureFirstJet=true;s.radial=rz;s.nonpolynomial=true;s.radius=rz?.34:.19;
        const auto geometry=CchMomentumGeometry::build(*mesh,metricsFor(*mesh,rz));const auto mass=buildMass(*mesh,geometry,s);
        const auto e=static_cast<double>(inspectMass(*mesh,s,mass,false));check(e>0,"zero nonpolynomial reference error");
        if(n>4){check(e<previous,"dual mass error did not decrease");order=std::log2(previous/e);}previous=e;++spaceRows;
        std::cout<<"SPACE,"<<rz<<','<<amr<<','<<n<<','<<e<<'\n';
    }
    check(order>=2,"dual geometry mass order below gate");std::cout<<"ORDER,"<<rz<<','<<amr<<','<<order<<'\n';
}
template<class F>void reject(F f,const std::string& why){rejects(f,why);++rejected;}
void declaredContacts()
{
    auto mesh=makeMesh(4,1);const auto geometry=CchMomentumGeometry::build(*mesh,metricsFor(*mesh,false));
    const Profile f=[](Vector2D p){return Jet{p.x+p.y-.7,1,1,0};};
    const auto seed=stateFor(*mesh,f,false,{},[](int,Vector2D,Vector2D){return 0.;},{},HermiteGalsTransport::FirstJet);
    HermiteGalsCandidate candidate(*mesh,{},seed.primary());const auto mixed=closeHermiteGalsMixedAmr(*mesh,candidate);
    HermiteGalsPublicationConfig config;config.volume.boundaryContacts={{0,.7},{.7,0}};
    const auto id=candidate.identity();const auto pub=publishHermiteGalsFirstJetGeometry(*mesh,std::move(candidate),mixed,{id,centers(*mesh,f)},config);
    const auto mass=CchMomentumMass::build(*mesh,geometry,CchPhaseGeometry::build(*mesh,geometry.metrics(),pub),1e6,1);
    std::array<long double,2> totals{};long double maximum=0;
    // Independent analytic vertical-slice integral of x+y<0.7.
    for(int u=0;u<mass.layout().velocityCount();++u) {
        const auto box=reference(*mesh,mass.layout(),u);const auto a=box.lower,b=box.upper;
        const long double fullEnd=std::max(a.x,std::min(b.x,.7-b.y));
        const long double lower=std::max(a.x,.7-b.y),upper=std::min(b.x,.7-a.y);
        long double vl=(fullEnd-a.x)*(b.y-a.y);
        if(upper>lower)vl+=(upper-lower)*(.7L-a.y-(upper+lower)/2);
        const auto volume=cap(a,b,false);
        maximum=std::max(maximum,norm(mass.phaseVolumes()[u].liquid,vl,volume,1e-10,"declared contact clipped triangle"));
        norm(mass.masses()[u],1e6*vl+volume-vl,1e6*volume,1e-10,"declared contact material mass");
        totals[axis(mass.layout(),u)]+=mass.phaseVolumes()[u].liquid;
    }
    for(auto total:totals)norm(total,.7L*.7L/2,1,1e-10,"declared contact triangle total");
    ++massScenarios;peakPhase=std::max(peakPhase,maximum);std::cout<<"CONTACT,"<<maximum<<'\n';
}
void overflow()
{
    DictionaryMesh d;d.Nx=d.Ny=4;d.maxCellNum=1024;d.maxInitialIterations=0;d.maxCoord={8,8};
    QuadtreeMesh mesh(d);Case s;s.pureFirstJet=true;s.planar=true;s.planeOffset=16;
    const auto geometry=CchMomentumGeometry::build(mesh,metricsFor(mesh,false));
    const auto phase=CchPhaseGeometry::build(mesh,geometry.metrics(),initial(mesh,s));
    reject([&]{CchMomentumMass::build(mesh,geometry,phase,std::numeric_limits<double>::max(),1);},
        std::numeric_limits<long double>::digits==std::numeric_limits<double>::digits?"material product overflow or underflow":"unrepresentable momentum mass");
}
void guards(bool rz)
{
    auto mesh=makeMesh(4);Case s;s.pureFirstJet=true;s.radial=rz;s.planar=true;s.planeOffset=.43;
    const auto metrics=metricsFor(*mesh,rz);const auto geometry=CchMomentumGeometry::build(*mesh,metrics);
    const auto phase=CchPhaseGeometry::build(*mesh,metrics,initial(*mesh,s));const auto mass=CchMomentumMass::build(*mesh,geometry,phase,1000,1);const auto saved=mass.masses();
    for(double bad:{0.,-1.,std::numeric_limits<double>::quiet_NaN(),std::numeric_limits<double>::infinity()}) {
        reject([&]{CchMomentumMass::build(*mesh,geometry,phase,bad,1);},"invalid material density");reject([&]{CchMomentumMass::build(*mesh,geometry,phase,1,bad);},"invalid material density");
    }
    const auto other=CchMomentumGeometry::build(*mesh,metricsFor(*mesh,rz));
    reject([&]{CchMomentumMass::build(*mesh,other,phase,1,1);},"different layout identity");
    const auto otherDi=CchMomentumGeometry::build(*mesh,CchPhysicalDivergence::build(*mesh,CchIntegralDivergence::build(*mesh,metrics.layout())));
    reject([&]{CchMomentumMass::build(*mesh,otherDi,phase,1,1);},"different retained D_I");
    reject([&]{mass.validateFor(*mesh,{},mass.timeLayer());},"different geometry identity");
    for(auto wrong:{HermiteGalsTimeLayer{1,0,0,Position::Initial},HermiteGalsTimeLayer{1,-.01,.01,Position::Midpoint},HermiteGalsTimeLayer{1,0,.1,Position::End}})
        reject([&]{mass.validateFor(*mesh,mass.geometryIdentity(),wrong);},"different full time layer");
    reject([&]{CchMomentumMass::build(*mesh,geometry,phase,std::numeric_limits<double>::denorm_min(),std::numeric_limits<double>::denorm_min());},
        std::numeric_limits<long double>::digits==std::numeric_limits<double>::digits?"material product overflow or underflow":"unrepresentable momentum mass");
    check(mass.masses()==saved,"failed build mutated retained masses");
    check(geometry.sameIdentity(mass.geometry()) && &mass.geometry().metrics().divergence().matrix()==&metrics.divergence().matrix(),"detached geometry or D_I");
    const auto copied=mass;check(&copied.masses()==&mass.masses() && &copied.geometry().segments()==&geometry.segments(),"copy rebuilt objects");
    auto equivalent=makeMesh(4);mesh.reset();copied.validateFor(*equivalent,copied.geometryIdentity(),copied.timeLayer());inspectMass(*equivalent,s,copied);
    auto& vertex=equivalent->verts(equivalent->vertIndexToLabel[0]);const auto old=vertex.coord();vertex.setCoord(old+Vector2D{.01,0});
    reject([&]{copied.validateFor(*equivalent,copied.geometryIdentity(),copied.timeLayer());},"NativeMeshIdentity");vertex.setCoord(old);
    DictionaryMesh::RefineBox box;box.level=1;box.min={.25,.25};box.max={.75,.75};equivalent->applyRefineBoxes({box});
    reject([&]{copied.validateFor(*equivalent,copied.geometryIdentity(),copied.timeLayer());},"NativeMeshIdentity");
    reject([&]{CchMomentumMass::build(*equivalent,geometry,phase,1,1);},"NativeMeshIdentity");
    const auto rebuilt=CchMomentumGeometry::build(*equivalent,metricsFor(*equivalent,rz));
    check(!rebuilt.sameIdentity(geometry),"regrid kept old geometry identity");inspectGeometry(*equivalent,rebuilt);inspectMass(*equivalent,s,buildMass(*equivalent,rebuilt,s));++rebuilds;
    std::cout<<"REGRID,"<<rz<<','<<rebuilt.layout().pressureCount()<<','<<rebuilt.layout().velocityCount()<<'\n';
}
}
void coincidentMass(){
    int scenarios=0,rows=0,overlaps=0;long double peak=0;
    for(bool rz:{false,true})for(int amr:{0,1})for(int sign:{-1,1})for(int selection:{0,1}){
        auto mesh=makeMesh(4,amr);Case s;s.pureFirstJet=true;s.radial=rz;s.planar=true;
        s.planeOffset=amr?.125:.375;s.planeSign=sign;
        const auto geometry=CchMomentumGeometry::build(*mesh,metricsFor(*mesh,rz,selection));
        const auto pub=initial(*mesh,s);const auto phase=CchPhaseGeometry::build(*mesh,geometry.metrics(),pub);
        const auto mass=CchMomentumMass::build(*mesh,geometry,phase,1000,1);
        check(pub.state().primary().phaseVolumes.empty(),"coincident mass introduced inventory");
        for(int u=0;u<mass.layout().velocityCount();++u){const auto box=reference(*mesh,mass.layout(),u);
            const auto capacity=cap(box.lower,box.upper,rz),vl=liquid(box.lower,box.upper,s),vg=capacity-vl;
            const auto m=mass.phaseVolumes()[u];check(m.liquid>=0 && m.gas>=0,"coincident dual negative phase");
            peak=std::max({peak,norm(m.liquid,vl,capacity,1e-12,"coincident dual liquid"),
                norm(m.gas,vg,capacity,1e-12,"coincident dual gas"),
                norm(mass.masses()[u],1000*vl+vg,1001*capacity,1e-12,"coincident true dual mass")});
            for(std::size_t k=0;k<box.pieces.size();++k){const auto& piece=box.pieces[k];
                const auto found=std::find_if(mass.geometry().volumes()[u].overlaps.begin(),mass.geometry().volumes()[u].overlaps.end(),
                    [&](const auto& p){return p.pressure==piece.pressure;});
                check(found!=mass.geometry().volumes()[u].overlaps.end(),"coincident actual overlap missing");
                const auto index=static_cast<std::size_t>(found-mass.geometry().volumes()[u].overlaps.begin());
                const auto actual=mass.overlapPhaseVolumes()[u][index];const auto expected=liquid(piece.lower,piece.upper,s);
                peak=std::max({peak,norm(actual.liquid,expected,piece.volume,1e-12,"coincident overlap liquid"),
                    norm(actual.gas,piece.volume-expected,piece.volume,1e-12,"coincident overlap gas")});++overlaps;
            }
            ++rows;
        }
        ++scenarios;
    }
    std::cout<<"COINCIDENT_MASS,"<<scenarios<<','<<rows<<','<<overlaps<<','<<peak<<'\n';
}
void nativeTraceMass()
{
    int publications=0,masses=0,duals=0,overlaps=0,queries=0;long double peak=0;std::size_t sideGraphs=0;
    const std::vector<HermiteGalsTimeLayer> times{{},{7,.2,.3,Position::Midpoint},{7,.2,.3,Position::End}};
    for(bool rz:{false,true})for(bool amr:{false,true})for(int normal:{0,1})for(int sign:{-1,1})for(const auto& time:times){
        auto mesh=makeMesh(4,amr?1:0);const double line=.5;
        const Profile f=[=](Vector2D p){const double n=normal==0?(rz?p.x*p.x-line*line:p.x-line):p.y-line;
            const double t=normal==0?p.y:p.x,a=rz && normal==1?1+t*t:1+t+t*t,da=rz && normal==1?2*t:1+2*t,dn=rz && normal==0?2*p.x:1;
            return Jet{sign*n*a,sign*(normal==0?dn*a:n*da),sign*(normal==1?a:n*da),sign*dn*da};};
        const auto seed=stateFor(*mesh,f,rz,time,[](int,Vector2D,Vector2D){return 0.;},{},HermiteGalsTransport::FirstJet);
        HermiteGalsCandidate candidate(*mesh,time,seed.primary());const auto id=candidate.identity();const auto mixed=closeHermiteGalsMixedAmr(*mesh,candidate);
        HermiteGalsPublicationConfig controls;controls.volume.boundaryContacts=normal==0?std::vector<Vector2D>{{line,0},{line,1}}:(rz?std::vector<Vector2D>{{1,line}}:std::vector<Vector2D>{{0,line},{1,line}});
        const auto pub=publishHermiteGalsFirstJetGeometry(*mesh,std::move(candidate),mixed,{id,centers(*mesh,f)},controls);
        const auto& scalar=pub.state().geometry().snapshot();const auto* address=&scalar;
        HermiteGalsPhaseRegionConfig config;config.boundaryContacts=controls.volume.boundaryContacts;
        const HermiteGalsPhaseRegionSource region(*mesh,pub.state(),config);
        check(region.nativeSideGraphs()==pub.report().graphs && region.nativeSideGraphs()>0 && &region.source().geometry().snapshot()==address,"native region source/topology changed");
        sideGraphs+=region.nativeSideGraphs();
        const auto referenceLiquid=[&](Vector2D lo,Vector2D hi){const auto capacity=cap(lo,hi,rz);auto clipped=hi;
            if(normal==0)clipped.x=std::max(lo.x,std::min(hi.x,line));else clipped.y=std::max(lo.y,std::min(hi.y,line));
            const auto prefix=cap(lo,clipped,rz);return sign>0?prefix:capacity-prefix;};
        for(const auto& box:std::vector<std::array<Vector2D,2>>{{{{0,0},{1,1}}},{{{.125,.125},{.875,.875}}},{{{.25,.25},{.5,.5}}}}){
            const auto m=region.measure(*mesh,box[0],box[1]);const auto expected=referenceLiquid(box[0],box[1]);
            peak=std::max({peak,norm(m.liquid,expected,m.capacity,1e-12,"native region liquid"),norm(m.gas,m.capacity-expected,m.capacity,1e-12,"native region gas"),
                norm(m.rawAxialLiquid,expected,m.capacity,1e-12,"native region axial"),norm(m.rawRadialLiquid,expected,m.capacity,1e-12,"native region radial")});
            check(m.sourceIdentity==id && m.integratedArcs==0,"native region passed side graph to regular integrator");++queries;
        }
        rejects([&]{region.measure(*mesh,{0,0},{1,2});},"active scalar support");
        const auto repeated=region.measure(*mesh,{0,0},{1,1});check(repeated.sourceIdentity==id && &region.source().geometry().snapshot()==address,"failed region query changed source");
        for(int selection:{0,1}) {
            const auto geometry=CchMomentumGeometry::build(*mesh,metricsFor(*mesh,rz,selection));const auto phase=CchPhaseGeometry::build(*mesh,geometry.metrics(),pub);
            for(double density:{1.,1000.,1e6}) {
                const auto mass=CchMomentumMass::build(*mesh,geometry,phase,density,1);
                check(mass.geometryIdentity()==id && mass.timeLayer().step==time.step && mass.timeLayer().time()==time.time()
                    && &mass.phaseGeometry().publication().state().geometry().snapshot()==address && pub.state().primary().phaseVolumes.empty(),"native dual material authority/inventory");
                for(int u=0;u<mass.layout().velocityCount();++u){const auto box=reference(*mesh,mass.layout(),u);const auto capacity=cap(box.lower,box.upper,rz);
                    const auto expected=referenceLiquid(box.lower,box.upper),gas=capacity-expected;const auto actual=mass.phaseVolumes()[u];
                    peak=std::max({peak,norm(actual.liquid,expected,capacity,1e-12,"native true dual liquid"),norm(actual.gas,gas,capacity,1e-12,"native true dual gas"),
                        norm(mass.masses()[u],density*expected+gas,(density+1)*capacity,1e-12,"native true dual mass")});
                    for(std::size_t k=0;k<geometry.volumes()[u].overlaps.size();++k){const auto& piece=geometry.volumes()[u].overlaps[k];const auto expectedPiece=referenceLiquid(piece.lower,piece.upper);
                        const auto m=mass.overlapPhaseVolumes()[u][k];peak=std::max({peak,norm(m.liquid,expectedPiece,piece.volume,1e-12,"native actual overlap liquid"),
                            norm(m.gas,piece.volume-expectedPiece,piece.volume,1e-12,"native actual overlap gas")});++overlaps;}
                    ++duals;
                }
                auto wrong=time;++wrong.step;rejects([&]{mass.validateFor(*mesh,id,wrong);},"different full time layer");++masses;
            }
        }
        ++publications;
    }
    check(publications==48 && masses==288 && queries==144 && sideGraphs>0,"native region/mass gate coverage");
    std::cout<<"NATIVE_TRACE_MASS,"<<publications<<','<<masses<<','<<duals<<','<<overlaps<<','<<queries<<','<<sideGraphs<<','<<peak<<",PASS\n";
}
void nativeTraceRegionExtras()
{
    int scenes=0;long double peak=0;
    for(bool amr:{false,true})for(int shape=0;shape<5;++shape)for(int sign:{-1,1}) {
        // Two physical exterior planes, their radial exterior counterpart,
        // an ordinary/side join, and two disconnected side components.
        const bool rz=shape==2;const double line=shape==0?0:1;
        auto mesh=makeMesh(4,amr?1:0);
        const Profile f=[=](Vector2D p){
            if(shape==3){const double v=std::max(0.,p.y-.5);return Jet{sign*(p.x-.5+v*v),static_cast<double>(sign),sign*2*v,0};}
            if(shape==4)return Jet{sign*(p.x-.25)*(p.x-.75),sign*(2*p.x-1),0,0};
            const double n=rz?p.x*p.x-line*line:p.x-line,a=1+p.y+p.y*p.y,da=1+2*p.y,dn=rz?2*p.x:1;
            return Jet{sign*n*a,sign*dn*a,sign*n*da,sign*dn*da};
        };
        const auto seed=stateFor(*mesh,f,rz,{},[](int,Vector2D,Vector2D){return 0.;},{},HermiteGalsTransport::FirstJet);
        HermiteGalsCandidate candidate(*mesh,{},seed.primary());const auto id=candidate.identity();const auto mixed=closeHermiteGalsMixedAmr(*mesh,candidate);
        HermiteGalsPublicationConfig controls;
        if(shape==3)controls.volume.boundaryContacts={{.5,0},{.25,1}};
        else if(shape==4)controls.volume.boundaryContacts={{.25,0},{.25,1},{.75,0},{.75,1}};
        else controls.volume.boundaryContacts={{line,0},{line,1}};
        const auto pub=publishHermiteGalsFirstJetGeometry(*mesh,std::move(candidate),mixed,{id,centers(*mesh,f)},controls);
        HermiteGalsPhaseRegionConfig config;config.boundaryContacts=controls.volume.boundaryContacts;
        const HermiteGalsPhaseRegionSource region(*mesh,pub.state(),config);
        const auto value=region.measure(*mesh,{0,0},{1,1});
        const long double prefix=shape==3?11.L/24:shape==4?.5L:line*(rz?pi:1.L);
        const auto expected=sign>0?prefix:value.capacity-prefix;
        peak=std::max({peak,norm(value.liquid,expected,value.capacity,1e-12,"exterior/join/component region liquid"),
            norm(value.gas,value.capacity-expected,value.capacity,1e-12,"exterior/join/component region gas"),
            norm(value.rawAxialLiquid,expected,value.capacity,1e-12,"exterior/join/component region axial"),
            norm(value.rawRadialLiquid,expected,value.capacity,1e-12,"exterior/join/component region radial")});
        check(value.sourceIdentity==id && region.nativeSideGraphs()>0 && region.components()==(shape==4?2:1)
            && &region.source().geometry().snapshot()==&pub.state().geometry().snapshot(),"exterior/join/component region provenance");
        if(shape<3){
            check(value.integratedArcs==0,"exterior native side entered region integration");
            const auto geometry=CchMomentumGeometry::build(*mesh,metricsFor(*mesh,rz));
            const auto mass=CchMomentumMass::build(*mesh,geometry,CchPhaseGeometry::build(*mesh,geometry.metrics(),pub),1000,1);
            for(int u=0;u<mass.layout().velocityCount();++u){const auto capacity=geometry.volumes()[u].volume,liquid=expected>0?capacity:0;
                peak=std::max(peak,norm(mass.masses()[u],1000*liquid+capacity-liquid,1001*capacity,1e-12,"missing exterior true dual mass"));}
        }
        ++scenes;
    }
    check(scenes==20,"native region extras coverage");
    std::cout<<"NATIVE_TRACE_REGION_EXTRAS,"<<scenes<<','<<peak<<",PASS\n";
}
namespace
{
static_assert(!std::is_default_constructible_v<CchAdvectorReconstruction>);
static_assert(!std::is_default_constructible_v<CchSpatialAdvector>);
static_assert(!std::is_constructible_v<CchAdvectorIdentity,std::uint64_t>);
// Physical divergence density: polynomial through degree five on each actual
// piece for both supported sources. GL3 is exact; the work scale retains each
// signed derivative contribution before cancellation. No volume-coordinate fit.
using PhysicalCuts=std::array<std::vector<double>,2>;
std::vector<PhysicalCuts> physicalPartitions(const CchAdvectorReconstruction& r)
{
    std::vector<PhysicalCuts> out(r.momentumGeometry().layout().pressureCount());
    for(const auto& p:r.pieces())for(int d=0;d<2;++d){out[p.pressure][d].push_back(d?p.lower.y:p.lower.x);out[p.pressure][d].push_back(d?p.upper.y:p.upper.x);}
    for(auto& cell:out)for(auto& line:cell){std::sort(line.begin(),line.end());line.erase(std::unique(line.begin(),line.end()),line.end());check(line.size()>=2,"missing physical partition");}
    return out;
}
template<class Query> std::array<long double,2> physicalVolumeIntegral(bool rz,int p,Vector2D lo,Vector2D hi,const PhysicalCuts& partition,Query query)
{
    PhysicalCuts cuts;
    for(int d=0;d<2;++d){const auto a=d?lo.y:lo.x,b=d?hi.y:hi.x;check(b>a,"empty physical integration region");cuts[d].push_back(a);for(auto k:partition[d])if(a<k&&k<b)cuts[d].push_back(k);cuts[d].push_back(b);}
    const long double g=std::sqrt(.6L),nodes[3]={-g,0,g},weights[3]={5.L/9,8.L/9,5.L/9};long double value=0,work=0;
    for(std::size_t x=1;x<cuts[0].size();++x)for(std::size_t y=1;y<cuts[1].size();++y){const long double hx=cuts[0][x]-cuts[0][x-1],hy=cuts[1][y]-cuts[1][y-1];
        for(int i=0;i<3;++i)for(int j=0;j<3;++j){const Vector2D point{static_cast<double>(cuts[0][x-1]+hx*(1+nodes[i])/2),static_cast<double>(cuts[1][y-1]+hy*(1+nodes[j])/2)};
            const auto a=query(p,point).jet;const long double xx=a.derivative.xx,yy=a.derivative.yy;
            const long double terms[3]={rz?point.x*xx:xx,rz?point.x*yy:yy,rz?a.value.x:0.L};const auto w=weights[i]*weights[j]*hx*hy/4*(rz?2*pi:1);
            for(auto term:terms){value+=w*term;work+=w*std::abs(term);}}
    }return {value,work};
}

// Six interior physical nodes differentiate all five-degree weighted values.
// Node coordinates use the ACTUAL rounded query positions. This reference uses
// values only, not production J, polynomial extraction, or reconstruction rows.
template<class Value> std::array<long double,4> measuredPhysicalJet(bool rz,int p,Vector2D point,Vector2D lower,Vector2D upper,Value value)
{
    const auto measured=[&](int direction,int c,int order){const long double a=component(lower,direction),h=component(upper-lower,direction),x=(component(point,direction)-a)/h;
        std::array<long double,6> nodes{},values{};
        for(int i=0;i<6;++i){auto q=point;const auto s=(1-std::cos((2*i+1)*pi/12))/2;(direction==0?q.x:q.y)=static_cast<double>(a+h*s);nodes[i]=(component(q,direction)-a)/h;
            const auto v=value(p,q);values[i]=component(v,c)*(rz&&c==0?static_cast<long double>(q.x):1.L);}
        long double result=0;
        for(int i=0;i<6;++i){long double denominator=1,weight=0;for(int j=0;j<6;++j)if(j!=i)denominator*=nodes[i]-nodes[j];
            for(int m=0;m<6;++m)if(m!=i){
                if(order==1){long double v=1;for(int j=0;j<6;++j)if(j!=i&&j!=m)v*=x-nodes[j];weight+=v;}
                else for(int n=0;n<6;++n)if(n!=i&&n!=m){long double v=1;for(int j=0;j<6;++j)if(j!=i&&j!=m&&j!=n)v*=x-nodes[j];weight+=v;}}
            result+=values[i]*weight/denominator;}
        return result/(order==1?h:h*h);
    };
    const long double radius=point.x;const auto v=value(p,point);
    if(rz&&radius==0)return {measured(0,0,2)/2,0,measured(0,1,1),measured(1,1,1)};
    const auto xx=measured(0,0,1),xy=measured(1,0,1),yx=measured(0,1,1),yy=measured(1,1,1);
    return rz?std::array<long double,4>{(xx-v.x)/radius,xy/radius,yx,yy}:std::array<long double,4>{xx,xy,yx,yy};
}

struct AdvectorReference
{
    const QuadtreeMesh& mesh;const CchDofLayout& l;bool rz;int kind;
    std::vector<Vector2D> lo,hi,mean;std::vector<double> vel;
    std::vector<std::vector<int>> faces;
    Vector2D analytic(Vector2D p) const {
        const double x=rz?.5*p.x*p.x:p.x,y=p.y;
        if(rz)return {p.x*(.13+kind*(.07*y*y+.05*x)),.2+.31*x+.19*y+kind*(.09*x*x-.07*x*y)};
        return {.3+.17*x-.23*y+kind*(.11*x*y+.08*y*y),-.2+.31*x+.19*y+kind*(.09*x*x-.07*x*y)};
    }
    template<class F> long double gauss(long double a,long double b,F f) const {
        const long double node=std::sqrt(.6L),nodes[3]={-node,0,node},weights[3]={5.L/9,8.L/9,5.L/9};long double v=0;
        for(int k=0;k<3;++k)v+=weights[k]*f(a+(b-a)*(1+nodes[k])/2);return v*(b-a)/2;
    }
    AdvectorReference(const QuadtreeMesh& m,const CchDofLayout& layout,int profile):mesh(m),l(layout),rz(l.domain().axisymmetric()),kind(profile),faces(l.pressureCount()) {
        for(int u=0;u<l.velocityCount();++u){const auto f=l.domain().face(mesh,l.nativeFace(u));const int ax=axis(l,u);const auto& e=f.horizontal?mesh.hEdges(f.edgeLabel):mesh.vEdges(f.edgeLabel);
            auto a=mesh.verts(e.edgeToVert[0]).coord(),b=mesh.verts(e.edgeToVert[1]).coord();if(component(a,1-ax)>component(b,1-ax))std::swap(a,b);lo.push_back(a);hi.push_back(b);
            const long double volume=gauss(component(a,1-ax),component(b,1-ax),[&](long double t){return rz?(ax==0?a.x:t):1.L;});
            const auto sample=gauss(component(a,1-ax),component(b,1-ax),[&](long double t){const auto p=ax==0?Vector2D{a.x,static_cast<double>(t)}:Vector2D{static_cast<double>(t),a.y};return component(analytic(p),ax)*(rz?p.x:1.L);});
            vel.push_back(volume>0?static_cast<double>(sample/volume):0);
            mean.push_back({rz?(ax==0?.5*a.x*a.x:.25*(a.x*a.x+b.x*b.x)):f.center.x,f.center.y});
            const int p=f.negativeCell<0?-1:l.pressureDof(f.negativeCell),q=f.positiveCell<0?-1:l.pressureDof(f.positiveCell);
            if(p>=0)faces[p].push_back(u);if(q>=0)faces[q].push_back(u);}
    }
};
void spatialAdvectorContract()
{
    int scenarios=0,updates=0,deficient=0,missing=0;std::size_t atoms=0,rows=0,views=0;long double peakMean=0,peakAffine=0,peakIndependent=0,peakDiv=0,dualAverageGap=0;
    for(bool rz:{false,true})for(int amr:{0,1})for(int n:{4,8,16,32}) {
        auto mesh=makeMesh(n,amr);Case shape;shape.pureFirstJet=true;shape.radial=rz;shape.planar=true;shape.planeOffset=2;const auto pub=initial(*mesh,shape);
        for(int selection=0;selection<4;++selection){const auto geometry=CchMomentumGeometry::build(*mesh,metricsFor(*mesh,rz,selection));const auto begin=std::chrono::steady_clock::now();
            const auto reconstruction=CchAdvectorReconstruction::prepare(*mesh,geometry);const auto copied=reconstruction;check(copied.sameIdentity(reconstruction) && geometry.sameIdentity(reconstruction.momentumGeometry()),"prepared geometry authority");
            const auto prepared=std::chrono::steady_clock::now();const auto& l=geometry.layout();const auto& report=reconstruction.report();
            std::size_t actualViews=0;for(const auto& v:reconstruction.primalAtoms())actualViews+=v.size();for(const auto& v:reconstruction.dualAtoms())actualViews+=v.size();
            check(report.atoms==reconstruction.atoms().size()&&report.viewEntries==actualViews,"actual prepared geometry inventory");
            deficient+=static_cast<int>(report.rankDeficient);missing+=static_cast<int>(report.tangentUnidentified);
            const auto partitions=physicalPartitions(reconstruction);
            const bool exactFamily=n>=8||selection<2; // Frozen original 56/64 analytical eligibility.
            std::optional<CchSpatialAdvector> first;
            for(int kind:{0,1}){AdvectorReference reference(*mesh,l,kind);const CchNativeValues values(l,reference.vel,std::vector<double>(l.pressureCount()));
                const auto advector=reconstruction.evaluate(*mesh,values,pub.state().identity(),pub.state().timeLayer());++updates;
                advector.validateFor(*mesh,geometry,pub.state().identity(),pub.state().timeLayer());
                check(advector.reconstruction().sameIdentity(reconstruction) && &advector.reconstruction().regularGeometry().rows()==&copied.regularGeometry().rows() && &advector.reconstruction().atoms()==&copied.atoms(),"update repeated geometry preparation");
                if(!first)first=advector;else check(first->identity()!=advector.identity() && first->values().normalVelocity()!=advector.values().normalVelocity(),"update failed to create actual revision");
                const auto oldDiv=geometry.metrics().apply(*mesh,values);
                for(int p=0;p<l.pressureCount();++p){long double gross=0;for(int u:reference.faces[p])gross+=std::abs(static_cast<long double>(geometry.metrics().faceMeasures()[u])*reference.vel[u]);
                    peakDiv=std::max(peakDiv,norm(advector.cellDivergence()[p]*geometry.metrics().cellVolumes()[p],oldDiv[p]*geometry.metrics().cellVolumes()[p],std::max(gross,1e-30L),1e-12,"original integral divergence"));
                    const auto& cell=mesh->cells(mesh->cellIndexToLabel[l.nativeCell(p)]);const auto a=mesh->verts(cell.cellToVert[0]).coord(),b=mesh->verts(cell.cellToVert[2]).coord();
                    const auto integral=physicalVolumeIntegral(rz,p,a,b,partitions[p],[&](int cell,Vector2D point){return advector.cellVelocityJet(cell,point);});
                    peakDiv=std::max(peakDiv,norm(integral[0],oldDiv[p]*geometry.metrics().cellVolumes()[p],std::max(gross+integral[1],1e-30L),1e-12,"independent physical cell Gauss integral"));
                    for(Vector2D f:std::vector<Vector2D>{{.21,.31},{.63,.77}}){const Vector2D point=a+Vector2D{f.x*(b.x-a.x),f.y*(b.y-a.y)};const auto actual=advector.cellVelocity(p,point),query=advector.cellVelocityJet(p,point).jet.value;
                        check(pointEqual(actual,query),"value/jet source mismatch");
                        if(kind==0 && exactFamily){const auto expected=reference.analytic(point);peakAffine=std::max({peakAffine,norm(actual.x,expected.x,1,1e-12,"cell affine x"),norm(actual.y,expected.y,1,1e-12,"cell affine y")});}++rows;}}
                for(std::size_t i=0;i<reconstruction.atoms().size();++i){const auto& atom=reconstruction.atoms()[i];const auto& poly=advector.atomPolynomials()[i];const auto ax=atom.normalAxis;
                    int ownerPressure=atom.ownerPressure;
                    if(atom.primalVelocity>=0){const auto face=l.domain().face(*mesh,l.nativeFace(atom.primalVelocity));ownerPressure=l.pressureDof(face.negativeCell>=0?face.negativeCell:face.positiveCell);}
                    check(ownerPressure>=0,"missing actual atom owner");
                    const auto referenceIntegral=reference.gauss(0,poly.length,[&](long double t){const Vector2D p=ax==0?Vector2D{atom.begin.x,static_cast<double>(atom.begin.y+t)}:Vector2D{static_cast<double>(atom.begin.x+t),atom.begin.y};
                        const auto value=component(advector.cellVelocity(ownerPressure,p),ax);
                        return value*(rz?2*pi*p.x:1.L);});
                    peakIndependent=std::max(peakIndependent,norm(poly.integral(false,0,poly.length),referenceIntegral,std::max(atom.area,1e-30L),1e-12,"independent atom integral"));
                    norm(poly.integral(true,0,poly.length),referenceIntegral,std::max(atom.area,1e-30L),1e-12,"end chart atom integral");
                    close(poly.integral(false,0,poly.length,-1),-poly.integral(false,0,poly.length),0,"normal reverse");
                    long double quarter=0;for(int k=0;k<4;++k)quarter+=poly.integral(false,poly.length*k/4,poly.length/4);
                    norm(quarter,referenceIntegral,std::max(atom.area,1e-30L),1e-12,"actual subinterval addback");}
                for(int u=0;u<l.velocityCount();++u){long double q=0;for(int atom:reconstruction.primalAtoms()[u])q+=advector.atomPolynomials()[atom].integral(false,0,advector.atomPolynomials()[atom].length);
                    peakMean=std::max(peakMean,norm(q,static_cast<long double>(geometry.metrics().faceMeasures()[u])*reference.vel[u],std::max(static_cast<long double>(geometry.metrics().faceMeasures()[u]),1e-30L),1e-12,"physical face mean"));
                    const auto v=advector.nativeVelocity(u);check(std::isfinite(v.x) && std::isfinite(v.y),"native point nonfinite");
                    long double dualFlux=0,dualGross=0,divIntegral=0,divWork=0,averageIntegral=0;
                    for(int side=0;side<4;++side)for(int segment:geometry.volumes()[u].sides[side])for(int atom:reconstruction.dualAtoms()[segment]){const auto& poly=advector.atomPolynomials()[atom];const auto q=poly.integral(false,0,poly.length);dualFlux+=(side%2?1:-1)*q;dualGross+=std::abs(q);}
                    for(const auto& overlap:geometry.volumes()[u].overlaps){const auto integral=physicalVolumeIntegral(rz,overlap.pressure,overlap.lower,overlap.upper,partitions[overlap.pressure],[&](int cell,Vector2D point){return advector.cellVelocityJet(cell,point);});divIntegral+=integral[0];divWork+=integral[1];averageIntegral+=overlap.volume*oldDiv[overlap.pressure];}
                    dualAverageGap=std::max(dualAverageGap,std::abs(dualFlux-averageIntegral));
                    norm(dualFlux,divIntegral,std::max(dualGross+divWork,1e-30L),1e-12,"actual dual boundary/divergence integral");
                    if(kind==0 && exactFamily){const auto expected=reference.analytic(l.domain().face(*mesh,l.nativeFace(u)).center);peakAffine=std::max({peakAffine,norm(v.x,expected.x,1,1e-12,"native affine x"),norm(v.y,expected.y,1,1e-12,"native affine y")});}}
                for(std::size_t g=0;g<geometry.segments().size();++g){const auto& s=geometry.segments()[g];long double measure=0;double cursor=component(s.begin,1-s.normalAxis);
                    for(int i:reconstruction.dualAtoms()[g]){const auto& atom=reconstruction.atoms()[i];check(component(atom.begin,1-s.normalAxis)==cursor && atom.primalVelocity==s.primalVelocity && atom.ownerPressure==s.ownerPressure,"dual atom actual origin/ordered coverage");cursor=component(atom.end,1-s.normalAxis);measure+=atom.area;++views;}
                    check(cursor==component(s.end,1-s.normalAxis),"dual atom gap");norm(measure,s.area,std::max(s.area,1e-30L),1e-12,"dual true area");}
            }
            atoms+=reconstruction.atoms().size();++scenarios;const auto done=std::chrono::steady_clock::now();
            std::cout<<"ADVECTOR_LAYOUT,"<<rz<<','<<amr<<','<<n<<','<<selection<<','<<l.velocityCount()<<','<<report.svdCount<<','<<report.weights<<','<<report.cellPieces<<','<<report.atoms<<','<<report.viewEntries<<','<<report.rankDeficient<<','<<report.tangentUnidentified<<','<<std::chrono::duration<double>(prepared-begin).count()<<','<<std::chrono::duration<double>(done-prepared).count()<<'\n';
        }
    }
    check(scenarios==64 && updates==128,"original full/selected coverage changed");
    std::cout<<"ADVECTOR_DUAL_AVERAGE_DIAGNOSTIC,"<<dualAverageGap<<'\n';
    std::cout<<"ADVECTOR_SUMMARY,"<<scenarios<<','<<updates<<','<<atoms<<','<<rows<<','<<views<<','<<deficient<<','<<missing<<','<<peakMean<<','<<peakAffine<<','<<peakIndependent<<','<<peakDiv<<",PASS\n";
}
void advectorGuards()
{
    for(bool rz:{false,true}){auto mesh=makeMesh(4,1);const auto geometry=CchMomentumGeometry::build(*mesh,metricsFor(*mesh,rz));const auto rec=CchAdvectorReconstruction::prepare(*mesh,geometry);
        Case shape;shape.pureFirstJet=true;shape.radial=rz;shape.planar=true;shape.planeOffset=2;const auto pub=initial(*mesh,shape);AdvectorReference ref(*mesh,geometry.layout(),0);const auto clock=pub.state().timeLayer();const auto id=pub.state().identity();
        const CchNativeValues values(geometry.layout(),ref.vel,std::vector<double>(geometry.layout().pressureCount()));const auto a=rec.evaluate(*mesh,values,id,clock);const auto copy=a;
        const auto atomBefore=a.atomPolynomials()[0].beginCoefficients;
        const CchNativeValues zeros(geometry.layout(),std::vector<double>(ref.vel.size()),std::vector<double>(geometry.layout().pressureCount()));const auto stationary=rec.evaluate(*mesh,zeros,id,clock);
        for(const auto& polynomial:stationary.atomPolynomials())check(std::all_of(polynomial.beginCoefficients.begin(),polynomial.beginCoefficients.end(),[](auto c){return c==0;}) && polynomial.integral(false,0,polynomial.length)==0,"zero native input produced moving trace");
        int zeroMeanFace=0;while(axis(geometry.layout(),zeroMeanFace)!=1)++zeroMeanFace;auto signedValues=ref.vel;const double target=ref.mean[zeroMeanFace].x;
        for(std::size_t u=0;u<signedValues.size();++u)signedValues[u]=axis(geometry.layout(),static_cast<int>(u))==1?ref.mean[u].x-target:0;
        const auto crossing=rec.evaluate(*mesh,CchNativeValues(geometry.layout(),signedValues,std::vector<double>(geometry.layout().pressureCount())),id,clock);
        long double net=0,gross=0;bool nonconstant=false;for(int i:rec.primalAtoms()[zeroMeanFace]){const auto& p=crossing.atomPolynomials()[i];net+=p.integral(false,0,p.length);gross+=p.length*(std::abs(p.value(false,0))+std::abs(p.value(true,0)));nonconstant=nonconstant || std::any_of(p.beginCoefficients.begin()+1,p.beginCoefficients.end(),[](auto c){return c!=0;});}
        check(signedValues[zeroMeanFace]==0 && nonconstant && gross>0,"zero face mean erased sign-changing normal polynomial");norm(net,0,gross,1e-12,"zero-net spatial trace");
        rejects([&]{rec.evaluate(*mesh,values,{},clock);},"unbound provenance");auto wrong=clock;wrong.step=-1;rejects([&]{rec.evaluate(*mesh,values,id,wrong);},"invalid step");
        wrong=clock;++wrong.step;rejects([&]{a.validateFor(*mesh,geometry,id,wrong);},"full time layer");wrong={7,.2,.3,Position::Midpoint};rejects([&]{a.validateFor(*mesh,geometry,id,wrong);},"full time layer");
        for(auto t:std::vector<HermiteGalsTimeLayer>{{7,.2,.3,Position::Midpoint},{7,.2,.3,Position::End}}){const auto revision=rec.evaluate(*mesh,values,id,t);revision.validateFor(*mesh,geometry,id,t);check(revision.identity()!=a.identity(),"clock revision identity");auto changed=t;changed.position=t.position==Position::End?Position::Midpoint:Position::End;rejects([&]{revision.validateFor(*mesh,geometry,id,changed);},"full time layer");}
        rejects([&]{a.validateFor(*mesh,CchMomentumGeometry::build(*mesh,geometry.metrics()),id,clock);},"momentum geometry instance");rejects([&]{a.validateFor(*mesh,geometry,{},clock);},"geometry identity");
        const auto fresh=CchMomentumGeometry::build(*mesh,metricsFor(*mesh,rz));rejects([&]{rec.evaluate(*mesh,CchNativeValues(fresh.layout(),ref.vel,std::vector<double>(fresh.layout().pressureCount())),id,clock);},"layout identity");
        auto bad=ref.vel;bad[0]=std::numeric_limits<double>::infinity();rejects([&]{CchNativeValues(geometry.layout(),bad,std::vector<double>(geometry.layout().pressureCount()));},"nonfinite native");
        bad=ref.vel;for(std::size_t u=0;u<bad.size();++u)bad[u]=(rz && geometry.metrics().faceMeasures()[u]==0)?0:std::numeric_limits<double>::max();rejects([&]{rec.evaluate(*mesh,CchNativeValues(geometry.layout(),bad,std::vector<double>(geometry.layout().pressureCount())),id,clock);},"nonfinite");
        if(rz){bad=ref.vel;bad[geometry.metrics().axisVelocityDofs()[0]]=1;rejects([&]{rec.evaluate(*mesh,CchNativeValues(geometry.layout(),bad,std::vector<double>(geometry.layout().pressureCount())),id,clock);},"nonzero radial axis");}
        const auto rounding=std::fegetround();std::fesetround(FE_UPWARD);try{rejects([&]{rec.evaluate(*mesh,values,id,clock);},"round-to-nearest");}catch(...){std::fesetround(rounding);throw;}std::fesetround(rounding);
        const auto control=_mm_getcsr();_mm_setcsr(control|(1U<<15));try{rejects([&]{a.nativeVelocity(0);},"gradual underflow");}catch(...){_mm_setcsr(control);throw;}_mm_setcsr(control);
        rejects([&]{a.cellVelocity(0,{-1,0});},"point outside");rejects([&]{a.atomPolynomials()[0].integral(false,0,-1);},"invalid shifted");
        // This width cannot be represented as a global endpoint increment.
        const auto found=std::find_if(a.atomPolynomials().begin(),a.atomPolynomials().end(),[](const auto& p){return p.value(false,.3L*p.length)>0 && (!p.radialMeasure || p.radiusBegin>0 || p.radialTangent);});check(found!=a.atomPolynomials().end(),"no positive tiny interval");const auto& poly=*found;const long double offset=.3L*poly.length,width=1e-25L;check(offset+width==offset && poly.integral(false,offset,width)>0,"tiny interval lost positive width");
        norm(poly.integral(false,offset,width),poly.value(false,offset)*width*(rz?2*pi*(poly.radiusBegin+(poly.radialTangent?offset:0)):1.L),width,1e-12,"tiny shifted polynomial");
        check(copy.identity()==a.identity() && copy.atomPolynomials()[0].beginCoefficients==atomBefore,"failed evaluation changed original owner");
        auto equivalent=makeMesh(4,1);mesh.reset();copy.validateFor(*equivalent,geometry,id,clock);check(std::isfinite(copy.nativeVelocity(0).y),"owner lifetime after mesh release");
        const auto retained=[&]{const auto g=CchMomentumGeometry::build(*equivalent,metricsFor(*equivalent,rz));const auto r=CchAdvectorReconstruction::prepare(*equivalent,g);const CchNativeValues v(g.layout(),ref.vel,std::vector<double>(g.layout().pressureCount()));return r.evaluate(*equivalent,v,id,clock);}();
        retained.validateFor(*equivalent,retained.momentumGeometry(),id,clock);check(retained.values().normalVelocity()==ref.vel && std::isfinite(retained.nativeVelocity(0).y),"owner lifetime after all original handles released");
        auto& vertex=equivalent->verts(equivalent->vertIndexToLabel[0]);const auto old=vertex.coord();vertex.setCoord(old+Vector2D{.01,0});rejects([&]{copy.validateFor(*equivalent,geometry,id,clock);},"NativeMeshIdentity");vertex.setCoord(old);
        DictionaryMesh::RefineBox box;box.level=2;box.min={.25,.25};box.max={.75,.75};equivalent->applyRefineBoxes({box});rejects([&]{copy.validateFor(*equivalent,geometry,id,clock);},"NativeMeshIdentity");
    }
    std::cout<<"ADVECTOR_GUARDS,PASS\n";
}
// Same physical-value derivative reference for every source representation.
struct JetAudit
{
    std::size_t rows=0,axisRows=0,knotRows=0,faces=0;
    long double derivative=0,affine=0,divergence=0,normal=0,tangentJump=0,knotJump=0;
};
void inspectCellJet(const CchSpatialAdvector& source,int p,Vector2D point,
    Vector2D lower,Vector2D upper,bool affine,JetAudit& report)
{
    const auto sample=source.cellVelocityJet(p,point);const auto& j=sample.jet.derivative;
    const auto value=source.cellVelocity(p,point);const bool rz=source.momentumGeometry().layout().domain().axisymmetric();
    check(pointEqual(sample.jet.value,value),"cell jet changed original value");
    check(pointEqual(sample.lower,lower) && pointEqual(sample.upper,upper),"cell jet selected wrong actual piece");
    const auto expected=measuredPhysicalJet(rz,p,point,lower,upper,[&](int cell,Vector2D x){return source.cellVelocity(cell,x);});
    const long double radius=point.x;
    const std::array<double,4> actual{j.xx,j.xy,j.yx,j.yy};
    for(int k=0;k<4;++k)report.derivative=std::max(report.derivative,
        norm(actual[k],expected[k],1+std::abs(expected[k]),2e-11,"measured same-branch Jacobian"));
    const long double div=j.xx+j.yy+(rz?(radius==0?j.xx:value.x/radius):0);
    const long double divScale=1+std::abs(j.xx)+std::abs(j.yy)+(rz?(radius==0?std::abs(j.xx):std::abs(value.x/radius)):0);
    report.divergence=std::max(report.divergence,std::abs(div-source.cellDivergence()[p])/divScale); // Diagnostic, not a pointwise constraint.
    if(affine) {
        const std::array<double,4> analytic=rz?std::array<double,4>{.13,0,.31*point.x,.19}:std::array<double,4>{.17,-.23,.31,.19};
        for(int k=0;k<4;++k)report.affine=std::max(report.affine,norm(actual[k],analytic[k],1,1e-12,"analytic affine Jacobian"));
    }
    if(rz && radius==0){check(value.x==0 && j.xy==0 && j.yx==0,"axis parity violated");++report.axisRows;}
    ++report.rows;
}
void spatialAdvectorJetContract()
{
    std::cout<<"CELL_JET_BEGIN\n";
    JetAudit total;int scenes=0;
    for(bool rz:{false,true})for(int amr:{0,1})for(int n:{4,8,16}) {
        auto mesh=makeMesh(n,amr);Case shape;shape.pureFirstJet=true;shape.radial=rz;shape.planar=true;shape.planeOffset=2;
        const auto publication=initial(*mesh,shape);
        for(int selection:{0,1,2,3}) {
            const auto geometry=CchMomentumGeometry::build(*mesh,metricsFor(*mesh,rz,selection));const auto& layout=geometry.layout();
            const auto reconstruction=CchAdvectorReconstruction::prepare(*mesh,geometry);const auto partitions=physicalPartitions(reconstruction);
            for(int kind:{0,1,2}) {
                AdvectorReference reference(*mesh,layout,kind==0?0:1);auto native=reference.vel;
                if(kind==2)for(int u=0;u<layout.velocityCount();++u)
                    native[u]=(rz && geometry.metrics().faceMeasures()[u]==0)?0:.1*std::sin(1.7*u)+.07*std::cos(.4*u);
                const auto source=reconstruction.evaluate(*mesh,CchNativeValues(layout,native,std::vector<double>(layout.pressureCount())),publication.state().identity(),publication.state().timeLayer());
                JetAudit local;
                for(int p=0;p<layout.pressureCount();++p) {
                    const auto& cell=mesh->cells(mesh->cellIndexToLabel[layout.nativeCell(p)]);
                    const auto lo=mesh->verts(cell.cellToVert[0]).coord(),hi=mesh->verts(cell.cellToVert[2]).coord();
                    const auto& xs=partitions[p][0];const auto& ys=partitions[p][1];
                    for(std::size_t ix=1;ix<xs.size();++ix)for(std::size_t iy=1;iy<ys.size();++iy) {
                        const Vector2D lower{xs[ix-1],ys[iy-1]},upper{xs[ix],ys[iy]};
                        // A true interior point, the lower knot, and each cell
                        // upper endpoint. No neighboring-side value averaging.
                        std::vector<Vector2D> points{lower,lower+Vector2D{.31*(upper.x-lower.x),.63*(upper.y-lower.y)}};
                        if(ix+1==xs.size())points.push_back({upper.x,(lower.y+upper.y)/2});
                        if(iy+1==ys.size())points.push_back({(lower.x+upper.x)/2,upper.y});
                        if(ix+1==xs.size() && iy+1==ys.size())points.push_back(upper);
                        for(const auto point:points)inspectCellJet(source,p,point,lower,upper,kind==0 && (n>=8||selection<2),local);
                        for(int direction=0;direction<2;++direction)if(component(lower,direction)>component(lo,direction)) {
                            auto right=lower+(upper-lower)*.5;auto left=right;
                            const double knot=component(lower,direction);
                            (direction==0?right.x:right.y)=knot;
                            (direction==0?left.x:left.y)=std::nextafter(knot,component(lo,direction));
                            const auto a=source.cellVelocityJet(p,left),b=source.cellVelocityJet(p,right);
                            check(component(a.upper,direction)==knot && component(b.lower,direction)==knot,"internal knot direction not retained");
                            local.normal=std::max(local.normal,norm(component(a.jet.value,direction),component(b.jet.value,direction),1,1e-12,"internal piece normal continuity"));
                            local.knotJump=std::max(local.knotJump,std::abs(static_cast<long double>(component(a.jet.value,1-direction))-component(b.jet.value,1-direction)));++local.knotRows;
                        }
                    }
                }
                // Independent actual-face neighbors, including both sides of
                // coarse/fine native segments, test the common normal trace.
                for(int u=0;u<layout.velocityCount();++u) {
                    const auto face=layout.domain().face(*mesh,layout.nativeFace(u));if(face.negativeCell<0 || face.positiveCell<0)continue;
                    const auto point=reference.lo[u]+(reference.hi[u]-reference.lo[u])*.37;const int normal=axis(layout,u);
                    const auto a=source.cellVelocityJet(layout.pressureDof(face.negativeCell),point),b=source.cellVelocityJet(layout.pressureDof(face.positiveCell),point);
                    local.normal=std::max(local.normal,norm(component(a.jet.value,normal),component(b.jet.value,normal),1,1e-12,"shared face normal trace"));
                    local.tangentJump=std::max(local.tangentJump,std::abs(static_cast<long double>(component(a.jet.value,1-normal))-component(b.jet.value,1-normal)));++local.faces;
                }
                ++scenes;total.rows+=local.rows;total.axisRows+=local.axisRows;total.knotRows+=local.knotRows;total.faces+=local.faces;
                total.derivative=std::max(total.derivative,local.derivative);total.affine=std::max(total.affine,local.affine);total.divergence=std::max(total.divergence,local.divergence);
                total.normal=std::max(total.normal,local.normal);total.tangentJump=std::max(total.tangentJump,local.tangentJump);total.knotJump=std::max(total.knotJump,local.knotJump);
                std::cout<<"CELL_JET,"<<rz<<','<<amr<<','<<n<<','<<selection<<','<<kind<<','<<local.rows<<','<<local.axisRows<<','<<local.knotRows<<','<<local.faces<<','<<local.derivative<<','<<local.affine<<','<<local.divergence<<','<<local.normal<<','<<local.tangentJump<<','<<local.knotJump<<'\n';
            }
        }
    }
    check(scenes==144 && total.axisRows>0 && total.knotRows>0 && total.faces>0,"missing actual branch/axis/jump coverage");
    std::cout<<"CELL_JET_SUMMARY,"<<scenes<<','<<total.rows<<','<<total.axisRows<<','<<total.knotRows<<','<<total.faces<<','<<total.derivative<<','<<total.affine<<','<<total.divergence<<','<<total.normal<<','<<total.tangentJump<<','<<total.knotJump<<",PASS\n";
}
void cellJetGuards()
{
    int guardRows=0;
    for(bool rz:{false,true}) {
        auto mesh=makeMesh(4,1);auto geometry=CchMomentumGeometry::build(*mesh,metricsFor(*mesh,rz));
        Case shape;shape.pureFirstJet=true;shape.radial=rz;shape.planar=true;shape.planeOffset=2;const auto pub=initial(*mesh,shape);const auto id=pub.state().identity();const auto clock=pub.state().timeLayer();
        const auto retained=[&] {const auto prepared=CchAdvectorReconstruction::prepare(*mesh,geometry);AdvectorReference ref(*mesh,geometry.layout(),1);
            return prepared.evaluate(*mesh,CchNativeValues(geometry.layout(),ref.vel,std::vector<double>(geometry.layout().pressureCount())),id,clock);}();
        const auto& cell=mesh->cells(mesh->cellIndexToLabel[geometry.layout().nativeCell(0)]);const auto point=cell.centerCoord();
        const auto before=retained.cellVelocityJet(0,point);const auto copy=retained;
        rejects([&]{retained.cellVelocityJet(-1,point);},"pressure index");rejects([&]{retained.cellVelocityJet(geometry.layout().pressureCount(),point);},"pressure index");
        rejects([&]{retained.cellVelocityJet(0,{-1,0});},"point outside");
        rejects([&]{retained.cellVelocityJet(0,{std::numeric_limits<double>::quiet_NaN(),0});},"point outside");
        rejects([&]{retained.cellVelocityJet(0,{point.x,std::numeric_limits<double>::infinity()});},"point outside");guardRows+=5;
        const auto mode=std::fegetround();std::fesetround(FE_UPWARD);try{rejects([&]{retained.cellVelocityJet(0,point);},"round-to-nearest");}catch(...){std::fesetround(mode);throw;}std::fesetround(mode);
        const auto control=_mm_getcsr();_mm_setcsr(control|(1U<<15));try{rejects([&]{retained.cellVelocityJet(0,point);},"gradual underflow");}catch(...){_mm_setcsr(control);throw;}_mm_setcsr(control);guardRows+=2;
        auto equivalent=makeMesh(4,1);mesh.reset();copy.validateFor(*equivalent,geometry,id,clock);const auto after=copy.cellVelocityJet(0,point);
        check(copy.identity()==retained.identity() && pointEqual(after.jet.value,before.jet.value) && after.jet.derivative.xx==before.jet.derivative.xx && after.jet.derivative.xy==before.jet.derivative.xy && after.jet.derivative.yx==before.jet.derivative.yx && after.jet.derivative.yy==before.jet.derivative.yy,"query failure or lifetime changed owner");
        auto changed=clock;++changed.step;rejects([&]{copy.validateFor(*equivalent,geometry,id,changed);},"full time layer");
        rejects([&]{copy.validateFor(*equivalent,geometry,{},clock);},"geometry identity");guardRows+=2;
        // Rebuild on a legal 2:1 mesh: the level-2 core needs a level-1 halo.
        DictionaryMesh::RefineBox halo;halo.level=1;halo.min={.125,.125};halo.max={.875,.875};equivalent->applyRefineBoxes({halo});
        DictionaryMesh::RefineBox box;box.level=2;box.min={.25,.25};box.max={.75,.75};equivalent->applyRefineBoxes({box});
        rejects([&]{copy.validateFor(*equivalent,geometry,id,clock);},"NativeMeshIdentity");++guardRows;
        // Regrid creates a NEW geometry/source; it does not retag old jets.
        const auto freshGeometry=CchMomentumGeometry::build(*equivalent,metricsFor(*equivalent,rz));const auto freshPub=initial(*equivalent,shape);
        const auto reconstruction=CchAdvectorReconstruction::prepare(*equivalent,freshGeometry);AdvectorReference ref(*equivalent,freshGeometry.layout(),0);
        const auto fresh=reconstruction.evaluate(*equivalent,CchNativeValues(freshGeometry.layout(),ref.vel,std::vector<double>(freshGeometry.layout().pressureCount())),freshPub.state().identity(),freshPub.state().timeLayer());
        fresh.validateFor(*equivalent,freshGeometry,freshPub.state().identity(),freshPub.state().timeLayer());check(fresh.identity()!=copy.identity(),"regrid reused advector source");
        const auto& freshCell=equivalent->cells(equivalent->cellIndexToLabel[freshGeometry.layout().nativeCell(0)]);const auto freshPoint=freshCell.centerCoord();
        check(pointEqual(fresh.cellVelocityJet(0,freshPoint).jet.value,fresh.cellVelocity(0,freshPoint)),"fresh regrid jet/value mismatch");
    }
    check(guardRows==20,"jet guard inventory changed");std::cout<<"CELL_JET_GUARDS,"<<guardRows<<",2,PASS\n";
}
// Preparation-only contract. The accepted advector remains covered above.
void regularPreparationContract()
{
    int guards=0;
    for(bool rz:{false,true}) {
        auto mesh=makeMesh(4,1);const auto geometry=CchMomentumGeometry::build(*mesh,metricsFor(*mesh,rz));
        const auto retained=[&]{const auto owner=CchAdvectorReconstruction::prepare(*mesh,geometry);return owner.regularGeometry();}();
        const auto copy=retained;check(copy.sameIdentity(retained)&&copy.momentumGeometry().sameIdentity(geometry),"regular geometry lost original owner");
        const auto& l=geometry.layout();std::vector<double> zero(l.velocityCount()),pressure(l.pressureCount());
        const CchNativeValues z(l,zero,pressure);const auto a=copy.apply(*mesh,z);
        for(const auto& c:a.cellDofs)for(const auto& v:c)for(double x:v)check(x==0,"regular zero field changed");
        for(double x:a.faceCorrections)check(x==0,"regular zero correction changed");
        check(copy.report().svdCount>0 && copy.report().weights>0 && copy.report().maxSupport<=24,"regular support/preparation inventory");
        // All physical P2 polynomials compatible with the RZ axis, summed.
        const auto fun=[&](Vector2D p,int axis){return rz?(axis==0?p.x+p.x*p.y:1+p.y+p.x*p.x+p.y*p.y):(axis==0?1+p.x+p.y+p.x*p.x+p.x*p.y+p.y*p.y:2-p.x+2*p.y);};
        const double q[5]={.046910077030668,.2307653449471585,.5,.7692346550528415,.953089922969332},w[5]={.1184634425280945,.2393143352496832,.2844444444444444,.2393143352496832,.1184634425280945};
        std::vector<double> u(l.velocityCount());
        for(int f=0;f<l.velocityCount();++f){const auto& face=copy.faces()[f];if(face.area==0)continue;for(int k=0;k<5;++k){const auto p=face.lower+(face.upper-face.lower)*q[k];u[f]+=w[k]*fun(p,face.normalAxis)*(rz?p.x/((face.lower.x+face.upper.x)/2):1);}}
        const CchNativeValues v(l,u,pressure);const auto coefficients=copy.apply(*mesh,v);
        for(std::size_t p=0;p<copy.cells().size();++p){const auto& c=copy.cells()[p];for(int axis=0;axis<2;++axis){const auto& raw=coefficients.cellDofs[p][axis];
            check(std::abs(raw[24]-fun((c.lower+c.upper)*.5,axis))<1e-11,"regular P2 center reproduction");}}
        for(double delta:coefficients.faceCorrections)check(std::abs(delta)<1e-11,"regular P2 face correction");
        const auto again=copy.apply(*mesh,v);check(again.cellDofs==coefficients.cellDofs&&again.faceCorrections==coefficients.faceCorrections,"regular repeated application changed");
        std::vector<int> walls;for(int f=0;f<l.velocityCount();++f)if(copy.faces()[f].external&&copy.faces()[f].area>0)walls.push_back(f);
        const auto law=CchAdvectorBoundary::impermeable(*mesh,geometry,walls);
        copy.apply(*mesh,z,&law);rejects([&]{copy.apply(*mesh,v,&law);},"nonzero prescribed wall");++guards;
        const auto other=CchMomentumGeometry::build(*mesh,metricsFor(*mesh,rz));
        const auto wrongLaw=CchAdvectorBoundary::impermeable(*mesh,other,walls);
        rejects([&]{copy.apply(*mesh,z,&wrongLaw);},"momentum geometry instance");
        rejects([&]{copy.validateFor(*mesh,other);},"momentum geometry instance");
        rejects([&]{copy.apply(*mesh,CchNativeValues(other.layout(),zero,pressure));},"layout identity");guards+=3;
        const auto mode=std::fegetround();std::fesetround(FE_UPWARD);try{rejects([&]{copy.apply(*mesh,z);},"round-to-nearest");}catch(...){std::fesetround(mode);throw;}std::fesetround(mode);
        const auto control=_mm_getcsr();_mm_setcsr(control|(1U<<15));try{rejects([&]{copy.apply(*mesh,z);},"gradual underflow");}catch(...){_mm_setcsr(control);throw;}_mm_setcsr(control);guards+=2;
        if(rz){auto bad=zero;for(int f=0;f<l.velocityCount();++f)if(copy.faces()[f].area==0){bad[f]=1;break;}rejects([&]{copy.apply(*mesh,CchNativeValues(l,bad,pressure));},"axis velocity");++guards;}
        mesh.reset();auto equivalent=makeMesh(4,1);const auto after=copy.apply(*equivalent,v);
        check(after.cellDofs==coefficients.cellDofs&&after.faceCorrections==coefficients.faceCorrections,"regular borrowed original mesh or caller values");
        DictionaryMesh::RefineBox halo;halo.level=1;halo.min={.125,.125};halo.max={.875,.875};equivalent->applyRefineBoxes({halo});
        DictionaryMesh::RefineBox box;box.level=2;box.min={.25,.25};box.max={.75,.75};equivalent->applyRefineBoxes({box});
        rejects([&]{copy.apply(*equivalent,v);},"NativeMeshIdentity");++guards;
        const auto fresh=CchMomentumGeometry::build(*equivalent,metricsFor(*equivalent,rz));const auto rebuilt=CchAdvectorReconstruction::prepare(*equivalent,fresh);
        check(!rebuilt.regularGeometry().sameIdentity(copy),"regular regrid reused geometry");
        rebuilt.regularGeometry().apply(*equivalent,CchNativeValues(fresh.layout(),std::vector<double>(fresh.layout().velocityCount()),std::vector<double>(fresh.layout().pressureCount())));
    }
    check(guards==15,"regular guard inventory");std::cout<<"REGULAR_PREPARATION,2,"<<guards<<",PASS\n";
}

// Complete public normal polynomial: reference is a physical-coordinate
// antiderivative, independent of production shifted/Horner integration.
void completeNormalPolynomialContract()
{
    using Coefficients=decltype(CchNormalPolynomial::beginCoefficients);
    static_assert(std::tuple_size<Coefficients>::value==5);
    const auto began=std::chrono::steady_clock::now();
    int ranges=0,tiny=0;long double peak=0;
    const auto reverse=[](const Coefficients& c,long double length) {
        Coefficients e{};
        for(int j=0;j<5;++j){int binomial=1;for(int k=0;k<=j;++k){e[k]+=c[j]*binomial*std::pow(length,j-k)*(k%2?-1:1);binomial=binomial*(j-k)/(k+1);}}
        return e;
    };
    const auto physicalIntegral=[](const Coefficients& c,long double lo,long double hi,int metric,long double radius,bool absolute) {
        long double result=0;
        for(int k=0;k<5;++k){const auto m=(std::pow(hi,k+1)-std::pow(lo,k+1))/(k+1);
            const auto rm=(std::pow(hi,k+2)-std::pow(lo,k+2))/(k+2);
            result+=(absolute?std::abs(c[k]):c[k])*(metric?CartesianMeasure::twoPi*(radius*m+(metric>=2?rm:0)):m);}
        return result;
    };
    std::cout<<"NORMAL_CONSTANT,"<<CartesianMeasure::twoPi<<'\n';
    for(int family=0;family<6;++family)for(long double length:{.5L,2.L})for(int metric=0;metric<4;++metric) {
        CchNormalPolynomial p;p.length=length;p.radialMeasure=metric>0;p.radialTangent=metric>=2;
        p.radiusBegin=metric==3?0:.25L;p.radiusEnd=p.radiusBegin+(p.radialTangent?length:0);
        if(family<5)p.beginCoefficients[family]=1;else p.beginCoefficients={1,-2,3,-4,5};
        p.endCoefficients=reverse(p.beginCoefficients,length);
        for(bool end:{false,true})for(auto range:std::array<std::array<long double,2>,3>{{{0,1},{.125L,.375L},{.75L,.125L}}})for(int sign:{-1,1}) {
            const auto off=length*range[0],width=length*range[1],lo=end?length-off-width:off,hi=lo+width;
            const auto expected=sign*physicalIntegral(p.beginCoefficients,lo,hi,metric,p.radiusBegin,false);
            const auto work=physicalIntegral(p.beginCoefficients,lo,hi,metric,p.radiusBegin,true);
            const auto actual=p.integral(end,off,width,sign);peak=std::max(peak,norm(actual,expected,work,1e-12,"complete normal physical antiderivative"));
            const auto x=end?length-off:off;long double value=0,valueWork=0;
            for(int k=0;k<5;++k){value+=p.beginCoefficients[k]*std::pow(x,k);valueWork+=std::abs(p.beginCoefficients[k])*std::pow(x,k);}
            const auto actualValue=p.value(end,off);peak=std::max(peak,norm(actualValue,value,valueWork,1e-12,"complete normal physical value"));
            std::cout<<"NORMAL_RANGE,"<<family<<','<<length<<','<<metric<<','<<end<<','<<off<<','<<width<<','<<sign<<','<<actual<<','<<expected<<','<<work<<','<<actualValue<<','<<value<<','<<valueWork<<'\n';++ranges;
        }
    }
    // At the exactly specified fourfold zero the retained integral is width^5,
    // and an axis-start radial segment additionally has a leading width^6.
    for(long double root:{0.L,.5L})for(int metric=0;metric<4;++metric)for(bool end:{false,true}) {
        if(root==0 && end)continue; // No positive interval beyond the physical lower endpoint.
        CchNormalPolynomial p;p.length=1;p.radialMeasure=metric>0;p.radialTangent=metric>=2;
        p.radiusBegin=metric==3?0:.25L;p.radiusEnd=p.radiusBegin+(p.radialTangent?1:0);
        p.beginCoefficients={std::pow(root,4),-4*std::pow(root,3),6*root*root,-4*root,1};p.endCoefficients=reverse(p.beginCoefficients,1);
        const auto offset=end?1-root:root;
        for(int exponent:{60,80})for(int sign:{-1,1}) {
            const auto width=std::ldexp(1.L,-exponent),base=std::pow(width,5)/5;
            const auto radius=p.radiusBegin+(p.radialTangent?root:0);
            const auto expected=sign*(metric?CartesianMeasure::twoPi*(radius*base+(p.radialTangent?(end?-1:1)*std::pow(width,6)/6:0)):base);
            const auto actual=p.integral(end,offset,width,sign);
            check(actual!=0 && (offset==0 || static_cast<double>(offset)+static_cast<double>(width)==static_cast<double>(offset)),"lost independent quartic width");
            check(p.value(end,offset)==0,"known fourfold zero not reproduced");
            peak=std::max(peak,norm(actual,expected,std::abs(expected),1e-12,"quartic tiny physical moment"));
            std::cout<<"NORMAL_TINY,"<<root<<','<<metric<<','<<end<<','<<exponent<<','<<sign<<','<<actual<<','<<expected<<'\n';++tiny;
        }
    }
    check(ranges==576 && tiny==48,"complete normal coverage");
    std::cout<<"NORMAL_COMPLETE_SUMMARY,"<<ranges<<','<<tiny<<','<<peak<<",PASS\n";
    std::cout<<"NORMAL_COMPLETE_SECONDS,"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-began).count()<<'\n';
}

// Independent Bernstein/de Casteljau reference for the same raw Q4 DOFs.
// It does not call production cardinal, derivative or polynomial extraction.
template<std::size_t N> long double regularBezier(std::array<long double,N> p,long double s)
{for(std::size_t n=N;n>1;--n)for(std::size_t i=0;i<n-1;++i)p[i]=(1-s)*p[i]+s*p[i+1];return p[0];}
std::array<long double,5> regularControls(const std::array<long double,5>& d)
{const auto b1=d[0]+d[2]/4,b3=d[1]-d[3]/4;return {d[0],b1,(16*d[4]-d[0]-4*b1-4*b3-d[1])/6,b3,d[1]};}
std::array<long double,6> regularJetValues(const CchAdvectorCellJet& a)
{const auto& j=a.jet;return {j.value.x,j.value.y,j.derivative.xx,j.derivative.xy,j.derivative.yx,j.derivative.yy};}
struct RegularReference
{
    const Quantum::Hydrodynamics::Detail::CchRegularAdvectorGeometry& geometry;
    const Quantum::Hydrodynamics::Detail::CchRegularCoefficients& data;
    std::vector<std::array<std::array<long double,25>,2>> controls;
    RegularReference(const Quantum::Hydrodynamics::Detail::CchRegularAdvectorGeometry& g,const Quantum::Hydrodynamics::Detail::CchRegularCoefficients& d):geometry(g),data(d),controls(d.cellDofs.size())
    {
        for(std::size_t p=0;p<controls.size();++p)for(int c=0;c<2;++c){std::array<long double,25> temp{};
            for(int y=0;y<5;++y){std::array<long double,5> row{};for(int x=0;x<5;++x)row[x]=d.cellDofs[p][c][5*y+x];row=regularControls(row);for(int x=0;x<5;++x)temp[5*y+x]=row[x];}
            for(int x=0;x<5;++x){std::array<long double,5> col{};for(int y=0;y<5;++y)col[y]=temp[5*y+x];col=regularControls(col);for(int y=0;y<5;++y)controls[p][c][5*y+x]=col[y];}}
    }
    std::array<long double,6> sample(int p,Vector2D point) const
    {
        const auto& cell=geometry.cells()[p];const auto h=cell.upper-cell.lower;const long double sx=(point.x-cell.lower.x)/h.x,sy=(point.y-cell.lower.y)/h.y;std::array<long double,6> out{};
        for(int c=0;c<2;++c){std::array<long double,5> value{},dx{};
            for(int y=0;y<5;++y){std::array<long double,5> row{};std::array<long double,4> derivative{};
                for(int x=0;x<5;++x)row[x]=controls[p][c][5*y+x];for(int x=0;x<4;++x)derivative[x]=4*(row[x+1]-row[x]);
                value[y]=regularBezier(row,sx);dx[y]=regularBezier(derivative,sx)/h.x;}
            std::array<long double,4> dy{};for(int y=0;y<4;++y)dy[y]=4*(value[y+1]-value[y]);
            out[c]=regularBezier(value,sy);out[2+2*c]=regularBezier(dx,sy);out[3+2*c]=regularBezier(dy,sy)/h.y;}
        for(int f:cell.faces){const auto& face=geometry.faces()[f];const int n=face.normalAxis,t=1-n;const auto delta=data.faceCorrections[f];if(delta==0)continue;
            const long double length=component(face.upper-face.lower,t),width=face.supportWidth[p==face.positivePressure?1:0];
            const auto s=(component(point,t)-component(face.lower,t))/length,xi=(component(point,n)-component(face.lower,n))/width,z=std::abs(xi);if(s<0||s>1||z>1)continue;
            const auto b=regularBezier(std::array<long double,5>{1,1,2.L/3,0,0},z),g=regularBezier(std::array<long double,5>{0,0,8.L/3,0,0},s);
            const auto db=regularBezier(std::array<long double,4>{0,-4.L/3,-8.L/3,0},z)*(xi<0?-1:1)/width;
            const auto dg=regularBezier(std::array<long double,4>{0,32.L/3,-32.L/3,0},s)/length;
            out[n]+=delta*b*g;out[2+2*n+n]+=delta*db*g;out[2+2*n+t]+=delta*b*dg;
        }
        return out;
    }
};
std::vector<double> regularCuts(const Quantum::Hydrodynamics::Detail::CchRegularAdvectorGeometry& g,int p,int direction,double lo,double hi,int other=-1)
{
    std::vector<double> cuts{lo,hi};for(int cell:{p,other})if(cell>=0)for(double knot:g.cells()[cell].supportKnots[direction])if(lo<knot&&knot<hi)cuts.push_back(knot);
    std::sort(cuts.begin(),cuts.end());cuts.erase(std::unique(cuts.begin(),cuts.end()),cuts.end());return cuts;
}
void localSourceStencilContract()
{
    std::size_t samples=0,units=0,maxSupport=0;long double peak=0;
    for(int rz=0;rz<2;++rz)for(int amr=0;amr<2;++amr)for(int selection=0;selection<4;++selection){
        auto mesh=makeMesh(4,amr);const auto momentum=CchMomentumGeometry::build(*mesh,metricsFor(*mesh,rz,selection));
        const auto owner=CchAdvectorReconstruction::prepare(*mesh,momentum);const auto& g=owner.regularGeometry();const auto& l=momentum.layout();
        std::vector<int> walls;for(int u=0;u<l.velocityCount();++u)if(g.faces()[u].external&&g.faces()[u].area>0)walls.push_back(u);
        const auto boundary=CchAdvectorBoundary::impermeable(*mesh,momentum,walls);
        for(int bound=0;bound<2;++bound){const auto* law=bound?&boundary:nullptr;std::vector<double> velocity(l.velocityCount());
            for(int u=0;u<l.velocityCount();++u)if(g.faces()[u].area>0&&(!bound||!std::binary_search(walls.begin(),walls.end(),u)))velocity[u]=std::sin(1.37*(u+1));
            const auto actual=g.apply(*mesh,CchNativeValues(l,velocity,std::vector<double>(l.pressureCount())),law);
            for(int p=0;p<l.pressureCount();++p){const auto stencil=g.cellStencil(*mesh,p,law);const auto& cell=g.cells()[p];
                for(Vector2D s:std::vector<Vector2D>{{0,0},{.31,.63},{1,1}}){const Vector2D x=cell.lower+Vector2D{s.x*(cell.upper.x-cell.lower.x),s.y*(cell.upper.y-cell.lower.y)};
                    const auto weights=g.cellJetWeights(stencil,x);maxSupport=std::max(maxSupport,weights.size());std::array<long double,6> recovered{};
                    for(const auto& w:weights){check(w.velocityDof>=0&&w.velocityDof<l.velocityCount(),"local stencil native index");for(int k=0;k<6;++k)recovered[k]+=w.jet[k]*velocity[w.velocityDof];}
                    const auto exact=regularJetValues(g.cellJet(actual,p,x));for(int k=0;k<6;++k)peak=std::max(peak,norm(recovered[k],exact[k],1+std::abs(exact[k]),1e-11,"local stencil actual source jet"));++samples;
                    if(p==0&&s.x==.31){for(int u=0;u<l.velocityCount();++u){if(g.faces()[u].area==0||(bound&&std::binary_search(walls.begin(),walls.end(),u)))continue;
                            std::vector<double> unit(l.velocityCount());unit[u]=1;const auto basis=g.apply(*mesh,CchNativeValues(l,unit,std::vector<double>(l.pressureCount())),law);
                            const auto expected=regularJetValues(g.cellJet(basis,p,x));std::array<double,6> observed{};for(const auto& w:weights)if(w.velocityDof==u)observed=w.jet;
                            for(int k=0;k<6;++k)peak=std::max(peak,norm(observed[k],expected[k],1+std::abs(expected[k]),1e-11,"local stencil independent native basis"));++units;}}
                }
            }
        }
    }
    check(samples>0&&units>0,"local source stencil coverage");
    std::cout<<"LOCAL_SOURCE_STENCIL,"<<samples<<','<<units<<','<<maxSupport<<','<<peak<<",PASS\n";
}
void normalFluxPortContract()
{
    int ports=0,traces=0;long double peak=0;
    constexpr long double node=.774596669241483377L;
    const long double quadrature[3]={(1-node)/2,.5L,(1+node)/2},weights[3]={5.L/18,4.L/9,5.L/18};
    for(int rz=0;rz<2;++rz)for(int amr=0;amr<2;++amr)for(int selection=0;selection<4;++selection){
        auto mesh=makeMesh(8,amr);const auto momentum=CchMomentumGeometry::build(*mesh,metricsFor(*mesh,rz,selection));
        const auto owner=CchAdvectorReconstruction::prepare(*mesh,momentum);const auto& g=owner.regularGeometry();const auto& l=momentum.layout();
        check(g.normalPorts().size()==g.report().normalPorts,"normal port preparation inventory");
        for(int c=0;c<2;++c)for(int degree=0;degree<=3;++degree){
            std::vector<double> velocity(l.velocityCount());
            for(int u=0;u<l.velocityCount();++u){const auto& face=g.faces()[u];if(face.normalAxis!=c||face.area==0)continue;
                for(int q=0;q<3;++q){const auto x=face.lower+(face.upper-face.lower)*static_cast<double>(quadrature[q]);
                    velocity[u]+=static_cast<double>(weights[q]*std::pow(component(x,1-c),degree)*(rz?x.x/((face.lower.x+face.upper.x)/2):1));}}
            const auto data=g.apply(*mesh,CchNativeValues(l,velocity,std::vector<double>(l.pressureCount())));
            for(const auto& port:g.normalPorts())if(port.component==c&&port.observedModes==4){
                const auto x=point(*mesh,port.vertex);int p=-1;
                for(int i=0;i<l.pressureCount();++i){const auto& cell=g.cells()[i];if(x.x>=cell.lower.x&&x.x<=cell.upper.x&&x.y>=cell.lower.y&&x.y<=cell.upper.y){p=i;break;}}
                check(p>=0,"normal-port vertex owner");const auto j=regularJetValues(g.cellJet(data,p,x));const long double t=component(x,1-c);
                peak=std::max(peak,norm(j[c],std::pow(t,degree),1,1e-11,"native line cubic value"));
                peak=std::max(peak,norm(j[2+2*c+1-c],degree?degree*std::pow(t,degree-1):0,1,1e-11,"native line cubic derivative"));++ports;}
        }
        // All native horizontal means on y=1/2 are zero; transverse
        // variation must not invent a normal trace or its tangent derivative.
        std::vector<double> velocity(l.velocityCount());
        for(int u=0;u<l.velocityCount();++u){const auto& f=g.faces()[u];if(f.normalAxis!=1||f.area==0)continue;
            for(int q=0;q<3;++q){const auto x=f.lower+(f.upper-f.lower)*static_cast<double>(quadrature[q]);
                velocity[u]+=static_cast<double>(weights[q]*(x.y-.5)*(1+x.x+x.x*x.x)*(rz?x.x/((f.lower.x+f.upper.x)/2):1));}}
        const auto data=g.apply(*mesh,CchNativeValues(l,velocity,std::vector<double>(l.pressureCount())));
        for(int p=0;p<l.pressureCount();++p){const auto& cell=g.cells()[p];if(cell.lower.y!=.5&&cell.upper.y!=.5)continue;
            for(double s:{0.,.17,.43,.81,1.}){const auto j=g.cellJet(data,p,{cell.lower.x+s*(cell.upper.x-cell.lower.x),.5}).jet;
                peak=std::max({peak,norm(j.value.y,0,1,1e-11,"zero normal face trace"),norm(j.derivative.yx,0,1,1e-11,"zero normal face tangent derivative")});++traces;}}
    }
    check(ports>0&&traces>0,"normal-port negative-control coverage");
    std::cout<<"NORMAL_FLUX_PORT,"<<ports<<','<<traces<<','<<peak<<",PASS\n";
}
void regularEvaluationContract(const std::vector<int>& sizes={4})
{
    using Clock=std::chrono::steady_clock;const auto began=Clock::now();double updateSeconds=0,querySeconds=0,kernelSeconds=0;
    int groups=0,fields=0,p2Layouts=0;std::size_t jets=0,atoms=0,c1=0,gauss=0,high=0,quartic=0,tiny=0,axisRows=0,wallRows=0,nativeAtoms=0,kernelQueries=0,maxScalarStorage=0;
    long double peakJet=0,peakAtom=0,peakC1=0,peakMean=0,peakGauss=0,peakP2=0,peakTiny=0,pointDivGap=0;
    const long double g=std::sqrt(3.L/5),nodes[3]={-g,0,g},weights[3]={5.L/9,8.L/9,5.L/9};
    Eigen::Matrix<double,5,5> interpolation;const int binomial[5]={1,4,6,4,1};
    for(int i=0;i<5;++i)for(int k=0;k<5;++k)interpolation(i,k)=binomial[k]*std::pow(i/4.,k)*std::pow(1-i/4.,4-k);
    const Eigen::Matrix<double,5,5> inverse=interpolation.inverse();
    for(bool rz:{false,true})for(int amr:{0,1})for(int n:sizes)for(int selected=0;selected<4;++selected){
        auto mesh=makeMesh(n,amr);const auto momentum=CchMomentumGeometry::build(*mesh,metricsFor(*mesh,rz,selected));
        const auto owner=CchAdvectorReconstruction::prepare(*mesh,momentum);const auto& geometry=owner.regularGeometry();const auto& l=momentum.layout();
        std::vector<int> walls;for(int u=0;u<l.velocityCount();++u)if(geometry.faces()[u].external&&geometry.faces()[u].area>0)walls.push_back(u);
        check(!walls.empty(),"missing actual regular wall fixture");
        const bool identifiable=std::all_of(geometry.rows().begin(),geometry.rows().end(),[](const auto& r){return r.zeroNormal||r.rank==static_cast<int>(r.columns.size());});
        if(identifiable)++p2Layouts;
        const auto polynomial=[&](Vector2D x){return rz?std::array<long double,6>{x.x*(1+x.y),1+x.y+x.x*x.x+x.y*x.y,1+x.y,x.x,2*x.x,1+2*x.y}
            :std::array<long double,6>{1+x.x+x.y+x.x*x.x+x.x*x.y+x.y*x.y,2-x.x+2*x.y,1+2*x.x+x.y,1+x.x+2*x.y,-1,2};};
        for(int bound=0;bound<3;++bound){auto chosen=bound==2?std::vector<int>{walls.front()}:walls;const auto law=CchAdvectorBoundary::impermeable(*mesh,momentum,chosen);
            const CchNativeValues zeros(l,std::vector<double>(l.velocityCount()),std::vector<double>(l.pressureCount()));const auto zero=geometry.apply(*mesh,zeros,bound?&law:nullptr);
            for(int p=0;p<l.pressureCount();++p){const auto& cell=geometry.cells()[p];const auto z=regularJetValues(geometry.cellJet(zero,p,(cell.lower+cell.upper)*.5));for(auto v:z)check(v==0,"regular zero jet not exact");}
            for(int family=0;family<2;++family){std::vector<double> velocity(l.velocityCount());
                for(int u=0;u<l.velocityCount();++u){const auto& f=geometry.faces()[u];if(f.area==0)continue;
                    for(int k=0;k<3;++k){const auto x=f.lower+(f.upper-f.lower)*static_cast<double>((1+nodes[k])/2);
                        const auto value=family==0?polynomial(x)[f.normalAxis]:std::sin(7*x.x+3*x.y)+.07*std::cos(2*x.x-5*x.y)+.03*std::sin(u+1.);
                        velocity[u]+=static_cast<double>(weights[k]/2*value*(rz?x.x/((f.lower.x+f.upper.x)/2):1));}}
                if(bound)for(int u:chosen)velocity[u]=0;
                const CchNativeValues values(l,velocity,std::vector<double>(l.pressureCount()));const auto updateStart=Clock::now();const auto data=geometry.apply(*mesh,values,bound?&law:nullptr);
                updateSeconds+=std::chrono::duration<double>(Clock::now()-updateStart).count();const RegularReference reference(geometry,data);
                const auto original=momentum.metrics().apply(*mesh,values);const auto queryStart=Clock::now();
                for(int p=0;p<l.pressureCount();++p){const auto& cell=geometry.cells()[p];const auto& xs=cell.supportKnots[0];const auto& ys=cell.supportKnots[1];long double volumeIntegral=0,volumeWork=0;
                    for(std::size_t ix=1;ix<xs.size();++ix)for(std::size_t iy=1;iy<ys.size();++iy){const Vector2D lo{xs[ix-1],ys[iy-1]},hi{xs[ix],ys[iy]};const auto h=hi-lo;
                        for(Vector2D fraction:std::vector<Vector2D>{{0,0},{.23,.61},{1,1}}){const auto point=lo+Vector2D{h.x*fraction.x,h.y*fraction.y};const auto actual=geometry.cellJet(data,p,point);const auto observed=regularJetValues(actual),expected=reference.sample(p,point);
                            check(point.x>=actual.lower.x&&point.x<=actual.upper.x&&point.y>=actual.lower.y&&point.y<=actual.upper.y,"regular piece closure mismatch");
                            for(int k=0;k<6;++k){peakJet=std::max(peakJet,norm(observed[k],expected[k],1+std::abs(expected[k]),1e-11,"independent regular Bernstein jet"));
                                if(family==0&&bound==0&&identifiable){const auto exact=polynomial(point);peakP2=std::max(peakP2,norm(observed[k],exact[k],1,1e-11,"regular identifiable P2 jet"));}}
                            if(rz&&point.x==0){check(observed[0]==0&&observed[3]==0&&observed[4]==0,"regular axis first-order parity");++axisRows;}++jets;}
                        for(int normal=0;normal<2;++normal){Vector2D a=lo,b=hi;(normal==0?a.x:a.y)=component(lo,normal)+.37*component(h,normal);(normal==0?b.x:b.y)=component(a,normal);
                            const auto poly=geometry.normalPolynomial(data,p,normal,a,b);long double integral=0,work=0;
                            for(int k=0;k<3;++k){const auto fraction=(1+nodes[k])/2;const auto point=a+(b-a)*static_cast<double>(fraction);const auto ref=reference.sample(p,point)[normal];
                                peakAtom=std::max({peakAtom,norm(poly.value(false,poly.length*fraction),ref,1+std::abs(ref),1e-11,"regular atom begin value"),norm(poly.value(true,poly.length*(1-fraction)),ref,1+std::abs(ref),1e-11,"regular atom end value")});
                                const auto measure=rz?CartesianMeasure::twoPi*point.x:1;integral+=weights[k]/2*poly.length*ref*measure;work+=weights[k]/2*poly.length*(1+std::abs(ref))*measure;}
                            peakAtom=std::max({peakAtom,norm(poly.integral(false,0,poly.length),integral,work,1e-11,"regular atom physical GL3"),norm(poly.integral(true,0,poly.length),integral,work,1e-11,"regular end physical GL3")});
                            check(poly.integral(false,0,poly.length,-1)==-poly.integral(false,0,poly.length),"regular reverse normal sharing");
                            if(poly.beginCoefficients[3]!=0||poly.beginCoefficients[4]!=0)++high;
                            if(poly.beginCoefficients[4]!=0)++quartic;
                            const auto offset=.375L*poly.length,width=std::ldexp(poly.length,-60);const auto point=a+(b-a)*.375;const auto ref=reference.sample(p,point)[normal];const auto measure=rz?CartesianMeasure::twoPi*point.x:1;
                            const auto small=poly.integral(false,offset,width);peakTiny=std::max(peakTiny,norm(small,ref*measure*width,(1+std::abs(ref))*measure*width,1e-11,"regular atom independent tiny width"));
                            check(static_cast<double>(offset)+static_cast<double>(width)==static_cast<double>(offset),"regular tiny fixture became global interval");if(std::abs(ref)>.05)check(small!=0,"regular tiny rate erased");++tiny;++atoms;}
                        for(int x=0;x<3;++x)for(int y=0;y<3;++y){const Vector2D point{lo.x+h.x*static_cast<double>((1+nodes[x])/2),lo.y+h.y*static_cast<double>((1+nodes[y])/2)};
                            const auto j=regularJetValues(geometry.cellJet(data,p,point));const auto density=rz?CartesianMeasure::twoPi*(point.x*(j[2]+j[5])+j[0]):j[2]+j[5];
                            const auto term=weights[x]*weights[y]*h.x*h.y*density/4;volumeIntegral+=term;volumeWork+=std::abs(term);++gauss;
                            pointDivGap=std::max(pointDivGap,std::abs((j[2]+j[5]+(rz?j[0]/point.x:0))-original[p]));}
                    }
                    long double gross=volumeWork;for(int u:cell.faces)gross+=std::abs(geometry.faces()[u].area*velocity[u]);
                    peakGauss=std::max(peakGauss,norm(volumeIntegral,original[p]*momentum.metrics().cellVolumes()[p],gross,1e-11,"regular Gauss / original DI"));
                }
                for(int u=0;u<l.velocityCount();++u){const auto& f=geometry.faces()[u];const int normal=f.normalAxis,t=1-normal,p=f.negativePressure>=0?f.negativePressure:f.positivePressure;
                    const auto cuts=regularCuts(geometry,p,t,component(f.lower,t),component(f.upper,t),f.negativePressure>=0&&f.positivePressure>=0?f.positivePressure:-1);long double integral=0,work=0;
                    for(std::size_t piece=1;piece<cuts.size();++piece){Vector2D a=f.lower,b=f.upper;(t==0?a.x:a.y)=cuts[piece-1];(t==0?b.x:b.y)=cuts[piece];
                        const auto nativePoly=geometry.normalPolynomial(data,p,f.normalAxis,a,b);long double pieceIntegral=0,pieceWork=0;
                        for(int k=0;k<3;++k){const auto point=a+(b-a)*static_cast<double>((1+nodes[k])/2);const auto v=geometry.cellJet(data,p,point).jet.value;const auto value=component(v,normal);const auto measure=rz?CartesianMeasure::twoPi*point.x:1;
                            pieceIntegral+=weights[k]/2*(cuts[piece]-cuts[piece-1])*value*measure;pieceWork+=weights[k]/2*(cuts[piece]-cuts[piece-1])*(1+std::abs(value))*measure;
                            if(bound&&std::find(chosen.begin(),chosen.end(),u)!=chosen.end()){check(value==0,"regular declared wall normal not zero");++wallRows;}}
                        integral+=pieceIntegral;work+=pieceWork;peakAtom=std::max({peakAtom,norm(nativePoly.integral(false,0,nativePoly.length),pieceIntegral,pieceWork,1e-11,"regular native polynomial query integral"),norm(nativePoly.integral(true,0,nativePoly.length),pieceIntegral,pieceWork,1e-11,"regular native end polynomial integral")});
                        if((rz&&f.area==0)||(bound&&std::find(chosen.begin(),chosen.end(),u)!=chosen.end()))for(const auto* coefficients:{&nativePoly.beginCoefficients,&nativePoly.endCoefficients})for(auto coefficient:*coefficients)check(coefficient==0,"regular axis/wall normal polynomial not exact zero");++nativeAtoms;
                        if(f.negativePressure>=0&&f.positivePressure>=0){std::array<std::array<long double,6>,5> diff{};long double scale=1;
                            for(int k=0;k<5;++k){const auto point=a+(b-a)*(k/4.);const auto left=regularJetValues(geometry.cellJet(data,f.negativePressure,point)),right=regularJetValues(geometry.cellJet(data,f.positivePressure,point));
                                for(int j=0;j<6;++j){diff[k][j]=left[j]-right[j];scale=std::max(scale,1+std::abs(left[j])+std::abs(right[j]));}}
                            for(int k=0;k<5;++k)for(int j=0;j<6;++j){long double coefficient=0;for(int i=0;i<5;++i)coefficient+=inverse(k,i)*diff[i][j];peakC1=std::max(peakC1,norm(coefficient,0,scale,1e-11,"regular full-degree internal C1 trace"));}++c1;}
                    }
                    peakMean=std::max(peakMean,norm(integral,f.area*velocity[u],work,1e-11,"regular true native face moment"));
                }
                querySeconds+=std::chrono::duration<double>(Clock::now()-queryStart).count();++fields;
                std::vector<CchAdvectorCellJet> measured(l.pressureCount());const auto kernelStart=Clock::now();
                for(int repeat=0;repeat<8;++repeat)for(int p=0;p<l.pressureCount();++p){const auto& cell=geometry.cells()[p];measured[p]=geometry.cellJet(data,p,(cell.lower+cell.upper)*.5);}
                kernelSeconds+=std::chrono::duration<double>(Clock::now()-kernelStart).count();kernelQueries+=8*measured.size();for(const auto& j:measured)check(std::isfinite(j.jet.value.x+j.jet.value.y),"regular timed query nonfinite");
                maxScalarStorage=std::max(maxScalarStorage,50*data.cellDofs.size()+data.faceCorrections.size());
                std::cout<<"REGULAR_EVALUATION_FIELD,"<<rz<<','<<amr<<','<<n<<','<<selected<<','<<bound<<','<<family<<','<<l.pressureCount()<<','<<l.velocityCount()<<','<<identifiable<<','<<geometry.report().rankDeficient<<','<<peakJet<<','<<peakAtom<<','<<peakC1<<','<<peakMean<<','<<peakGauss<<','<<peakP2<<','<<peakTiny<<'\n';
            }++groups;
        }
    }
    check(groups==48*static_cast<int>(sizes.size())&&fields==2*groups&&high>0&&quartic>0&&tiny==atoms&&axisRows>0&&wallRows>0&&nativeAtoms>0&&pointDivGap>1e-3,"regular evaluation coverage or constant-divergence overconstraint");
    std::cout<<"REGULAR_EVALUATION_SUMMARY,"<<groups<<','<<fields<<','<<p2Layouts<<','<<jets<<','<<atoms<<','<<c1<<','<<gauss<<','<<high<<','<<tiny<<','<<axisRows<<','<<wallRows<<','<<peakJet<<','<<peakAtom<<','<<peakC1<<','<<peakMean<<','<<peakGauss<<','<<peakP2<<','<<peakTiny<<','<<pointDivGap<<",PASS\n";
    std::cout<<"REGULAR_EVALUATION_SECONDS,"<<updateSeconds<<','<<querySeconds<<','<<std::chrono::duration<double>(Clock::now()-began).count()<<'\n';
    std::cout<<"REGULAR_EVALUATION_KERNEL,"<<kernelQueries<<','<<kernelSeconds<<','<<maxScalarStorage<<','<<sizeof(Quantum::Hydrodynamics::Detail::CchRegularCoefficients)<<'\n';
    std::cout<<"REGULAR_EVALUATION_NATIVE,"<<nativeAtoms<<','<<quartic<<",PASS\n";
}
void regularEvaluationGuards()
{
    int count=0;for(bool rz:{false,true}){
        auto mesh=makeMesh(4,1);const auto g=CchMomentumGeometry::build(*mesh,metricsFor(*mesh,rz));const auto original=CchAdvectorReconstruction::prepare(*mesh,g);
        const auto prepared=original.regularGeometry();const auto& l=g.layout();const CchNativeValues zeros(l,std::vector<double>(l.velocityCount()),std::vector<double>(l.pressureCount()));const auto data=prepared.apply(*mesh,zeros);
        const auto& cell=prepared.cells()[0];const auto point=(cell.lower+cell.upper)*.5;const auto checkReject=[&](auto call,const char* text){rejects(call,text);++count;};
        Quantum::Hydrodynamics::Detail::CchRegularCoefficients unbound;
        checkReject([&]{prepared.cellJet(unbound,0,point);},"coefficient geometry identity");
        const auto other=CchAdvectorReconstruction::prepare(*mesh,g);checkReject([&]{other.regularGeometry().cellJet(data,0,point);},"coefficient geometry identity");
        auto wrong=data;wrong.cellDofs.pop_back();checkReject([&]{prepared.cellJet(wrong,0,point);},"coefficient inventory");
        wrong=data;wrong.faceCorrections.pop_back();checkReject([&]{prepared.cellJet(wrong,0,point);},"coefficient inventory");
        checkReject([&]{prepared.cellJet(data,-1,point);},"pressure index");checkReject([&]{prepared.cellJet(data,l.pressureCount(),point);},"pressure index");
        checkReject([&]{prepared.cellJet(data,0,{-1,point.y});},"outside regular cell");checkReject([&]{prepared.cellJet(data,0,{std::numeric_limits<double>::quiet_NaN(),point.y});},"outside regular cell");
        wrong=data;wrong.cellDofs[0][0][0]=std::numeric_limits<double>::infinity();checkReject([&]{prepared.cellJet(wrong,0,point);},"nonfinite regular Q4 datum");
        wrong=data;wrong.faceCorrections[cell.faces[0]]=std::numeric_limits<double>::quiet_NaN();checkReject([&]{prepared.cellJet(wrong,0,point);},"nonfinite regular face correction");
        Vector2D a{point.x,cell.lower.y},b{point.x,cell.upper.y};
        checkReject([&]{prepared.normalPolynomial(data,0,2,a,b);},"normal axis");checkReject([&]{prepared.normalPolynomial(data,0,0,b,a);},"canonical regular normal interval");
        checkReject([&]{prepared.normalPolynomial(data,0,0,a,a);},"canonical regular normal interval");checkReject([&]{prepared.normalPolynomial(data,0,0,a,b);},"crosses regular support knot");
        const auto mode=std::fegetround();std::fesetround(FE_UPWARD);try{checkReject([&]{prepared.cellJet(data,0,point);},"round-to-nearest");}catch(...){std::fesetround(mode);throw;}std::fesetround(mode);
        const auto control=_mm_getcsr();_mm_setcsr(control|(1U<<15));try{checkReject([&]{prepared.cellJet(data,0,point);},"gradual underflow");}catch(...){_mm_setcsr(control);throw;}_mm_setcsr(control);
        auto temporary=CchAdvectorReconstruction::prepare(*mesh,g).regularGeometry();const auto retained=temporary.apply(*mesh,zeros);mesh.reset();
        const auto result=regularJetValues(temporary.cellJet(retained,0,point));for(auto value:result)check(value==0,"regular query borrowed released mesh");
        const auto again=regularJetValues(prepared.cellJet(data,0,point));for(auto value:again)check(value==0,"failed regular query changed original data");
    }check(count==32,"regular query guard inventory");std::cout<<"REGULAR_EVALUATION_GUARDS,"<<count<<",PASS\n";
}

// H08a2r10: declared physical low-space basis vectors and two independent
// mixtures. Eligibility is the ORIGINAL affine/RZ gate, not full P2 rank.
void regularLowerSpaceContract(const std::vector<int>& sizes={4,8})
{
    int layouts=0,eligible=0,fields=0,priorityMissing=0;std::size_t queries=0;
    long double peak=0,supplemental=0,deficit=0;double maxNorm=0;std::size_t maxSupport=0;
    for(bool rz:{false,true})for(int amr:{0,1})for(int n:sizes)for(int selection=0;selection<4;++selection){
        auto mesh=makeMesh(n,amr);const auto momentum=CchMomentumGeometry::build(*mesh,metricsFor(*mesh,rz,selection));
        const auto owner=CchAdvectorReconstruction::prepare(*mesh,momentum);const auto& g=owner.regularGeometry();const auto& l=momentum.layout();
        const bool gate=owner.report().tangentUnidentified==0;const int dimension=rz?4:6;
        maxNorm=std::max(maxNorm,g.report().maxCoefficientNorm);maxSupport=std::max(maxSupport,g.report().maxSupport);
        for(const auto& r:g.rows()){
            check(r.priorityRank>=0&&r.priorityRank<=r.priorityDimension&&std::isfinite(r.priorityDeficit),"invalid priority information report");
            if(r.priorityRank==r.priorityDimension){check(r.priorityDeficit<=1e-11,"identifiable lower span lost reproduction");deficit=std::max(deficit,static_cast<long double>(r.priorityDeficit));}
            else ++priorityMissing;
        }
        long double layoutPeak=0;
        for(int family=0;family<dimension+2;++family){
            const auto exact=[&](Vector2D x){std::array<long double,6> value{};
                for(int basis=0;basis<dimension;++basis){const long double weight=family<dimension?(family==basis?1:0):family==dimension?(basis%2?-1.L:1.L)*(basis+1)/7:(basis%2?1.L:-1.L)*(dimension-basis)/11;
                    std::array<long double,6> b{};
                    if(!rz){const int c=basis/3,k=basis%3;b[c]=k==0?1:k==1?x.x:x.y;if(k>0)b[2+2*c+k-1]=1;}
                    else if(basis==0){b[0]=x.x;b[2]=1;}else if(basis==1)b[1]=1;else if(basis==2){b[1]=x.x*x.x;b[4]=2*x.x;}else{b[1]=x.y;b[5]=1;}
                    for(int k=0;k<6;++k)value[k]+=weight*b[k];}
                return value;
            };
            std::vector<double> velocity(l.velocityCount());const long double node=std::sqrt(.6L),nodes[3]={-node,0,node},weights[3]={5.L/9,8.L/9,5.L/9};
            for(int u=0;u<l.velocityCount();++u){const auto& f=g.faces()[u];if(f.area==0)continue;long double sum=0;
                for(int k=0;k<3;++k){const auto x=f.lower+(f.upper-f.lower)*static_cast<double>((1+nodes[k])/2);sum+=weights[k]/2*exact(x)[f.normalAxis]*(rz?x.x/((f.lower.x+f.upper.x)/2):1);}velocity[u]=static_cast<double>(sum);}
            const auto data=g.apply(*mesh,CchNativeValues(l,velocity,std::vector<double>(l.pressureCount())));
            const auto inspect=[&](int p,Vector2D point,bool native){const auto a=regularJetValues(g.cellJet(data,p,point)),b=exact(point);++queries;
                for(int k=0;k<(native?2:6);++k){check(std::isfinite(a[k]),"nonfinite lower physical jet");if(!gate)continue;
                    const auto error=std::abs(a[k]-b[k]);if(k<2||n<=16){check(error<=1e-12,"original lower physical precision");layoutPeak=std::max(layoutPeak,error);}else supplemental=std::max(supplemental,error);}
            };
            for(int p=0;p<l.pressureCount();++p){const auto& c=g.cells()[p];const auto h=c.upper-c.lower;for(Vector2D f:std::vector<Vector2D>{{.21,.31},{.63,.77}})inspect(p,c.lower+Vector2D{f.x*h.x,f.y*h.y},false);}
            for(int u=0;u<l.velocityCount();++u){const auto& f=g.faces()[u];inspect(f.negativePressure>=0?f.negativePressure:f.positivePressure,(f.lower+f.upper)*.5,true);}
            ++fields;
        }
        ++layouts;if(gate)++eligible;peak=std::max(peak,layoutPeak);
        std::cout<<"LOWER_SPACE_LAYOUT,"<<rz<<','<<amr<<','<<n<<','<<selection<<','<<gate<<','<<g.report().rankDeficient<<','<<layoutPeak<<",PASS\n";
    }
    check(layouts==16*static_cast<int>(sizes.size())&&fields==112*static_cast<int>(sizes.size())&&priorityMissing>0,"lower physical coverage or missing information erased");
    std::cout<<"LOWER_SPACE_SUMMARY,"<<layouts<<','<<eligible<<','<<fields<<','<<queries<<','<<priorityMissing<<','<<peak<<','<<supplemental<<','<<deficit<<','<<maxNorm<<','<<maxSupport<<",PASS\n";
}

// A deficient whole affine span must not exempt its compatible constants or
// observable axial variation. References here are physical fields/integrals,
// independent of the prepared weights and of the older 56-layout eligibility.
void regularIdentifiableConstantContract(const std::vector<int>& sizes={4})
{
    int layouts=0,fields=0,missing=0;std::size_t queries=0,polynomials=0;
    long double peakValue=0,peakJac=0,peakCoefficient=0,peakIntegral=0;
    for(bool rz:{false,true})for(int amr:{0,1})for(int n:sizes)for(int selection=0;selection<4;++selection){
        auto mesh=makeMesh(n,amr);const auto momentum=CchMomentumGeometry::build(*mesh,metricsFor(*mesh,rz,selection));
        const auto owner=CchAdvectorReconstruction::prepare(*mesh,momentum);const auto& g=owner.regularGeometry();const auto& l=momentum.layout();
        for(const auto& row:g.rows())if(row.priorityRank<row.priorityDimension)++missing;
        for(int family=0;family<4;++family){
            const auto exact=[&](Vector2D x)->std::array<long double,6>{
                if(family==0)return rz?std::array<long double,6>{x.x,0,1,0,0,0}:std::array<long double,6>{1,0,0,0,0,0};
                if(family==1)return {0,1,0,0,0,0};
                if(family==2)return {0,x.y,0,0,0,1};
                return {rz?.3*x.x:.3+.2*x.x,-.1+.4*x.y,rz?.3:.2,0,0,.4};
            };
            std::vector<int> walls;
            for(int u=0;u<l.velocityCount();++u){const auto& f=g.faces()[u];
                if(f.external&&f.area>0&&exact(f.lower)[f.normalAxis]==0&&exact(f.upper)[f.normalAxis]==0)walls.push_back(u);}
            const auto law=CchAdvectorBoundary::impermeable(*mesh,momentum,walls);
            for(int bound:{0,1})for(double scale:(family==1&&n==4&&selection==3?std::vector<double>{0,1,1e-150,1e150}:std::vector<double>{1})){
                std::vector<double> velocity(l.velocityCount());
                for(int u=0;u<l.velocityCount();++u){const auto& f=g.faces()[u];if(f.area>0)velocity[u]=scale*static_cast<double>(exact((f.lower+f.upper)*.5)[f.normalAxis]);}
                const auto data=g.apply(*mesh,CchNativeValues(l,velocity,std::vector<double>(l.pressureCount())),bound?&law:nullptr);
                const auto error=[&](long double actual,long double expected,const char* message){
                    check(std::isfinite(actual),"nonfinite compatible low-space field");
                    if(scale==0){check(actual==0,"zero compatible low-space field");return 0.L;}
                    return norm(actual,scale*expected,std::abs(scale),1e-12,message);
                };
                for(int p=0;p<l.pressureCount();++p){const auto& c=g.cells()[p];const auto h=c.upper-c.lower;
                    for(Vector2D f:std::vector<Vector2D>{{.13,.17},{.37,.43},{.81,.79}}){const auto x=c.lower+Vector2D{f.x*h.x,f.y*h.y};
                        const auto a=regularJetValues(g.cellJet(data,p,x)),b=exact(x);++queries;
                        for(int k=0;k<6;++k){const auto e=error(a[k],b[k],"compatible constant/observable affine value or J");
                            (k<2?peakValue:peakJac)=std::max(k<2?peakValue:peakJac,e);}}
                }
                const auto integral=[&](int p,int normal,Vector2D a,Vector2D b){
                    const int t=1-normal;const auto cuts=regularCuts(g,p,t,component(a,t),component(b,t));long double sum=0;
                    for(std::size_t i=1;i<cuts.size();++i){auto lo=a,hi=b;(t==0?lo.x:lo.y)=cuts[i-1];(t==0?hi.x:hi.y)=cuts[i];
                        const auto poly=g.normalPolynomial(data,p,normal,lo,hi);const auto begin=exact(lo),end=exact(hi);const auto slope=begin[2+2*normal+t];
                        for(int k=0;k<5;++k){const auto dimensional=std::pow(poly.length,k);
                            peakCoefficient=std::max({peakCoefficient,error(poly.beginCoefficients[k]*dimensional,(k==0?begin[normal]:k==1?slope*poly.length:0),"compatible begin polynomial"),
                                error(poly.endCoefficients[k]*dimensional,(k==0?end[normal]:k==1?-slope*poly.length:0),"compatible end polynomial")});}
                        long double expected=0;const long double node=std::sqrt(.6L),nodes[3]={-node,0,node},weights[3]={5.L/9,8.L/9,5.L/9};
                        for(int k=0;k<3;++k){const auto x=lo+(hi-lo)*static_cast<double>((1+nodes[k])/2);
                            expected+=weights[k]*poly.length/2*exact(x)[normal]*(rz?2*pi*x.x:1);}
                        peakIntegral=std::max({peakIntegral,error(poly.integral(false,0,poly.length),expected,"compatible physical begin integral"),error(poly.integral(true,0,poly.length),expected,"compatible physical end integral")});
                        check(poly.integral(false,0,poly.length,-1)==-poly.integral(false,0,poly.length),"compatible opposite integral sharing");
                        sum+=poly.integral(false,0,poly.length);++polynomials;
                    }return sum;
                };
                for(int u=0;u<l.velocityCount();++u){const auto& f=g.faces()[u];const int p=f.negativePressure>=0?f.negativePressure:f.positivePressure;
                    peakIntegral=std::max(peakIntegral,error(integral(p,f.normalAxis,f.lower,f.upper),f.area*exact((f.lower+f.upper)*.5)[f.normalAxis],"compatible native mean"));}
                for(const auto& s:momentum.segments()){int p=s.ownerPressure;
                    if(s.primalVelocity>=0){const auto& f=g.faces()[s.primalVelocity];p=f.negativePressure>=0?f.negativePressure:f.positivePressure;}
                    integral(p,s.normalAxis,s.begin,s.end);}
                ++fields;
            }
        }++layouts;
    }
    check(layouts==16*static_cast<int>(sizes.size())&&fields==8*layouts+24&&missing>0,"compatible constant coverage or rank diagnostic lost");
    std::cout<<"IDENTIFIABLE_CONSTANT_SUMMARY,"<<layouts<<','<<fields<<','<<missing<<','<<queries<<','<<polynomials<<','<<peakValue<<','<<peakJac<<','<<peakCoefficient<<','<<peakIntegral<<",PASS\n";
}

template<class Query> std::array<long double,2> physicalTraceIntegral(bool rz,int p,int normal,Vector2D a,Vector2D b,const PhysicalCuts& partition,Query query)
{
    const int t=1-normal;const long double lo=component(a,t),hi=component(b,t);std::vector<long double> cuts{lo};
    for(auto k:partition[t])if(lo<k&&k<hi)cuts.push_back(k);cuts.push_back(hi);long double sum=0,work=0;
    const long double g=std::sqrt(.6L),nodes[3]={-g,0,g},weights[3]={5.L/9,8.L/9,5.L/9};
    for(std::size_t k=1;k<cuts.size();++k)for(int j=0;j<3;++j){const auto width=cuts[k]-cuts[k-1];auto x=a;(t==0?x.x:x.y)=static_cast<double>(cuts[k-1]+width*(1+nodes[j])/2);
        const long double value=component(query(p,x).jet.value,normal),w=weights[j]*width/2*(rz?2*pi*x.x:1);sum+=w*value;work+=w*std::abs(value);}
    return {sum,work};
}
// Identical mathematical consumer for the accepted source and passive regular
// data. No source-kind switch affects any formula, threshold or eligibility.
void matchingSpatialConsumers(const std::vector<int>& sizes={4})
{
    int layouts=0,fields=0,signMutation=0,degreeMutation=0,supportMutation=0;std::size_t derivatives=0,cells=0,duals=0,native=0;
    long double peakDerivative=0,peakCell=0,peakDual=0,peakMean=0,peakAnalytic=0;
    std::array<long double,2> pointGap{},dualGap{};std::array<long double,3> mutationGap{};
    for(bool rz:{false,true})for(int amr:{0,1})for(int n:sizes){auto mesh=makeMesh(n,amr);Case shape;shape.pureFirstJet=true;shape.radial=rz;shape.planar=true;shape.planeOffset=2;const auto pub=initial(*mesh,shape);
        for(int selection=0;selection<4;++selection){const auto momentum=CchMomentumGeometry::build(*mesh,metricsFor(*mesh,rz,selection));const auto owner=CchAdvectorReconstruction::prepare(*mesh,momentum);const auto& g=owner.regularGeometry();const auto& l=momentum.layout();
            const auto oldCuts=physicalPartitions(owner);std::vector<PhysicalCuts> regularPartitions;for(const auto& c:g.cells())regularPartitions.push_back(c.supportKnots);
            std::vector<int> walls;for(int u=0;u<l.velocityCount();++u)if(g.faces()[u].external&&g.faces()[u].area>0)walls.push_back(u);check(!walls.empty(),"missing matching boundary fixture");
            for(int boundary=0;boundary<3;++boundary){auto chosen=boundary==2?std::vector<int>{walls.front()}:walls;const auto law=CchAdvectorBoundary::impermeable(*mesh,momentum,chosen);
                for(int family=0;family<3;++family){AdvectorReference input(*mesh,l,family==0?0:1);auto velocity=input.vel;
                    if(family==2)for(int u=0;u<l.velocityCount();++u)velocity[u]=g.faces()[u].area==0?0:.1*std::sin(1.7*u)+.07*std::cos(.4*u);
                    if(boundary)for(int u:chosen)velocity[u]=0;
                    const CchNativeValues values(l,velocity,std::vector<double>(l.pressureCount()));const auto mean=momentum.metrics().apply(*mesh,values);
                    const auto accepted=boundary?owner.evaluate(*mesh,values,pub.state().identity(),pub.state().timeLayer(),law):owner.evaluate(*mesh,values,pub.state().identity(),pub.state().timeLayer());
                    const auto data=g.apply(*mesh,values,boundary?&law:nullptr);
                    const auto regularQuery=[&](int p,Vector2D x){return g.cellJet(data,p,x);};
                    const auto inspect=[&](int label,const auto& partition,auto query){
                        const bool exact=boundary==0&&family==0&&(n>=8||selection<2);
                        for(int p=0;p<l.pressureCount();++p){const auto& cell=g.cells()[p];const auto integral=physicalVolumeIntegral(rz,p,cell.lower,cell.upper,partition[p],query);long double gross=0;for(int u:cell.faces)gross+=std::abs(g.faces()[u].area*velocity[u]);
                            peakCell=std::max(peakCell,norm(integral[0],mean[p]*momentum.metrics().cellVolumes()[p],std::max(gross+integral[1],1e-30L),1e-12,"matching physical cell integral"));++cells;
                            const auto& xs=partition[p][0];const auto& ys=partition[p][1];
                            for(std::size_t ix=1;ix<xs.size();++ix)for(std::size_t iy=1;iy<ys.size();++iy){const Vector2D lo{xs[ix-1],ys[iy-1]},hi{xs[ix],ys[iy]};
                                for(Vector2D f:std::vector<Vector2D>{{0,0},{.31,.63}}){const auto x=lo+Vector2D{f.x*(hi.x-lo.x),f.y*(hi.y-lo.y)};const auto sample=query(p,x);const auto& j=sample.jet;const auto expected=measuredPhysicalJet(rz,p,x,lo,hi,[&](int q,Vector2D y){return query(q,y).jet.value;});
                                    const std::array<double,4> actual{j.derivative.xx,j.derivative.xy,j.derivative.yx,j.derivative.yy};for(int k=0;k<4;++k)peakDerivative=std::max(peakDerivative,norm(actual[k],expected[k],1+std::abs(expected[k]),2e-11,"matching value-derived physical J"));++derivatives;
                                    const long double div=j.derivative.xx+j.derivative.yy+(rz?(x.x==0?j.derivative.xx:j.value.x/x.x):0);
                                    pointGap[label]=std::max(pointGap[label],std::abs(div-mean[p]));
                                    if(exact){const auto a=input.analytic(x);const std::array<double,4> derivative=rz?std::array<double,4>{.13,0,.31*x.x,.19}:std::array<double,4>{.17,-.23,.31,.19};
                                        peakAnalytic=std::max({peakAnalytic,norm(j.value.x,a.x,1,1e-12,"matching analytic x"),norm(j.value.y,a.y,1,1e-12,"matching analytic y")});for(int k=0;k<4;++k)peakAnalytic=std::max(peakAnalytic,norm(actual[k],derivative[k],1,1e-12,"matching analytic J"));}
                                }}
                        }
                        for(int u=0;u<l.velocityCount();++u){const auto& face=g.faces()[u];const int p=face.negativePressure>=0?face.negativePressure:face.positivePressure;
                            const auto q=physicalTraceIntegral(rz,p,face.normalAxis,face.lower,face.upper,partition[p],query);const auto expected=face.area*velocity[u];peakMean=std::max(peakMean,norm(q[0],expected,std::max(q[1]+std::abs(expected),1e-30L),1e-12,"matching native physical mean"));++native;
                            if(exact){const auto x=(face.lower+face.upper)*.5,a=query(p,x).jet.value,b=input.analytic(x);peakAnalytic=std::max({peakAnalytic,norm(a.x,b.x,1,1e-12,"matching native x"),norm(a.y,b.y,1,1e-12,"matching native y")});}
                            long double flux=0,work=0,volume=0,average=0;
                            for(int side=0;side<4;++side)for(int s:momentum.volumes()[u].sides[side]){const auto& segment=momentum.segments()[s];int pressure=segment.ownerPressure;
                                if(segment.primalVelocity>=0){const auto& f=g.faces()[segment.primalVelocity];pressure=f.negativePressure>=0?f.negativePressure:f.positivePressure;}
                                const auto integral=physicalTraceIntegral(rz,pressure,segment.normalAxis,segment.begin,segment.end,partition[pressure],query);flux+=(side%2?1:-1)*integral[0];work+=integral[1];}
                            for(const auto& overlap:momentum.volumes()[u].overlaps){const auto integral=physicalVolumeIntegral(rz,overlap.pressure,overlap.lower,overlap.upper,partition[overlap.pressure],query);volume+=integral[0];work+=integral[1];average+=overlap.volume*mean[overlap.pressure];}
                            peakDual=std::max(peakDual,norm(flux,volume,std::max(work,1e-30L),1e-12,"matching physical dual integral"));dualGap[label]=std::max(dualGap[label],std::abs(flux-average));++duals;
                        }
                    };
                    inspect(0,oldCuts,[&](int p,Vector2D x){return accepted.cellVelocityJet(p,x);});inspect(1,regularPartitions,regularQuery);++fields;
                    // Perturb ONLY test-side consumers. Their independent physical
                    // references must see the wrong sign, missing degree/support.
                    if(!signMutation||!degreeMutation||!supportMutation)for(int p=0;p<l.pressureCount();++p){const auto& c=g.cells()[p];const auto& cuts=regularPartitions[p];
                        if(!supportMutation){const auto good=physicalVolumeIntegral(rz,p,c.lower,c.upper,cuts,regularQuery),bad=physicalVolumeIntegral(rz,p,c.lower,c.upper,PhysicalCuts{},regularQuery);const auto error=std::abs(good[0]-bad[0])/(1+good[1]);if(error>1e-8){supportMutation=1;mutationGap[2]=error;}}
                        if(!signMutation||!degreeMutation){const Vector2D a{(cuts[0][0]+cuts[0][1])*.5,cuts[1][0]},b{a.x,cuts[1][1]};const auto poly=g.normalPolynomial(data,p,0,a,b);const auto reference=physicalTraceIntegral(rz,p,0,a,b,cuts,regularQuery);
                            const auto signError=std::abs(-poly.integral(false,0,poly.length)-reference[0])/(1+reference[1]);if(signError>1e-8){signMutation=1;mutationGap[0]=signError;}
                            auto truncated=poly;truncated.beginCoefficients[3]=truncated.beginCoefficients[4]=0;const auto degreeError=std::abs(truncated.integral(false,0,poly.length)-reference[0])/(1+reference[1]);if(degreeError>1e-8){degreeMutation=1;mutationGap[1]=degreeError;}}
                    }
                }
            }
            ++layouts;std::cout<<"MATCHED_LAYOUT,"<<rz<<','<<amr<<','<<n<<','<<selection<<",PASS\n";
        }
    }
    check(layouts==16*static_cast<int>(sizes.size())&&fields==9*layouts&&signMutation&&degreeMutation&&supportMutation,"matching coverage or ineffective negative control");
    check(pointGap[1]>1e-3&&dualGap[1]>1e-5,"regular fixture did not distinguish released constraints");
    std::cout<<"MATCHED_SUMMARY,"<<layouts<<','<<fields<<','<<derivatives<<','<<cells<<','<<duals<<','<<native<<','<<peakDerivative<<','<<peakCell<<','<<peakDual<<','<<peakMean<<','<<peakAnalytic<<",PASS\n";
    std::cout<<"MATCHED_RELEASED_DIAGNOSTICS,"<<pointGap[0]<<','<<pointGap[1]<<','<<dualGap[0]<<','<<dualGap[1]<<'\n';
    std::cout<<"MATCHED_NEGATIVE_CONTROLS,"<<signMutation<<','<<degreeMutation<<','<<supportMutation<<','<<mutationGap[0]<<','<<mutationGap[1]<<','<<mutationGap[2]<<",PASS\n";
}

void phaseSourceSupportContract()
{
    for(bool rz:{false,true})for(int amr:{0,1}){
        auto mesh=makeMesh(4,amr);
        const auto geometry=CchMomentumGeometry::build(*mesh,metricsFor(*mesh,rz));
        const auto source=CchAdvectorReconstruction::prepare(*mesh,geometry);
        const auto& regular=source.regularGeometry();
        constexpr double plane=11./32.;
        int cut=-1;double centerDistance=2;
        for(std::size_t p=0;p<regular.cells().size();++p){const auto& cell=regular.cells()[p];
            const double distance=std::abs((cell.lower.y+cell.upper.y)/2-.5);
            if(cell.lower.x<plane&&plane<cell.upper.x&&distance<centerDistance){cut=static_cast<int>(p);centerDistance=distance;}}
        check(cut>=0,"phase support fixture has no cut cell");
        const auto& cell=regular.cells()[cut];
        const auto fraction=[&](double lo,double hi){const long double end=std::clamp(plane,lo,hi);
            return static_cast<double>(rz?(end*end-static_cast<long double>(lo)*lo)/
                (static_cast<long double>(hi)*hi-static_cast<long double>(lo)*lo):(end-lo)/(hi-lo));};
        const double alpha=fraction(cell.lower.x,cell.upper.x);
        const Eigen::Vector2d hidden(1-alpha,-alpha);
        Eigen::Matrix<double,2,2> local;local.row(0)<<alpha,1-alpha;local.row(1)=local.row(0);
        Eigen::JacobiSVD<Eigen::MatrixXd> localSvd(local,Eigen::ComputeThinU|Eigen::ComputeThinV);
        check(localSvd.rank()==1&&(local*hidden).norm()<=1e-15,"single-cell negative control is not hidden");
        const auto stencil=regular.cellStencil(*mesh,cut);std::set<int> support;
        for(const auto& dof:stencil.cellDofs[1])for(const auto& weight:dof)
            if(weight.weight!=0&&regular.faces()[weight.velocityDof].normalAxis==1)support.insert(weight.velocityDof);
        for(const auto& correction:stencil.faceCorrections)for(const auto& weight:correction)
            if(weight.weight!=0&&regular.faces()[weight.velocityDof].normalAxis==1)support.insert(weight.velocityDof);
        Eigen::MatrixXd observed(support.size(),2);Eigen::VectorXd native(support.size());
        std::size_t row=0,outsideNonzero=0;
        for(int u:support){const auto& face=regular.faces()[u];const double f=fraction(face.lower.x,face.upper.x);
            observed.row(row)<<f,1-f;native[row]=f*hidden[0]+(1-f)*hidden[1];
            if(std::find(cell.faces.begin(),cell.faces.end(),u)==cell.faces.end()&&std::abs(native[row])>1e-12)++outsideNonzero;
            ++row;}
        Eigen::JacobiSVD<Eigen::MatrixXd> patch(observed,Eigen::ComputeThinU|Eigen::ComputeThinV);
        patch.setThreshold(128*std::numeric_limits<double>::epsilon()*std::max<Eigen::Index>(observed.rows(),2));
        check(patch.rank()==2&&outsideNonzero>0&&native.norm()>1e-3,
            "actual source support did not observe the two phase constants");
        const auto recovered=patch.solve(native);
        const double coefficientError=(recovered-hidden).cwiseAbs().maxCoeff();
        check(coefficientError<=1e-12,"independent phase constant fit failed on actual source support");
        std::cout<<"PHASE_SOURCE_SUPPORT,"<<rz<<','<<amr<<','<<cut<<','<<alpha<<','<<support.size()
            <<','<<localSvd.rank()<<','<<patch.rank()<<','<<outsideNonzero<<','<<patch.singularValues()[1]
            <<','<<coefficientError<<",PASS\n";
    }
}

long double independentModalAtomRate(const Case& shape,const CchSpatialAdvector& state,const CchAdvectorAtom& atom)
{
    int p=atom.ownerPressure;
    if(atom.primalVelocity>=0){const auto& face=state.reconstruction().regularGeometry().faces()[atom.primalVelocity];
        p=face.negativePressure>=0?face.negativePressure:face.positivePressure;}
    const int tangent=1-atom.normalAxis;const double lo=component(atom.begin,tangent),hi=component(atom.end,tangent);
    std::vector<double> cuts{lo,hi};const auto add=[&](double x){if(x>lo&&x<hi)cuts.push_back(x);};
    if(shape.planar){if(atom.normalAxis!=(shape.axial()?1:0))add(shape.planeOffset);}
    else{const Vector2D center{shape.radial?0:shape.cx,shape.cy};const double delta=component(atom.begin,atom.normalAxis)-component(center,atom.normalAxis);
        const double squared=shape.radius*shape.radius-delta*delta;
        if(squared>0){add(component(center,tangent)-std::sqrt(squared));add(component(center,tangent)+std::sqrt(squared));}}
    std::sort(cuts.begin(),cuts.end());
    constexpr double nodes[8]={-.96028985649753623,-.79666647741362674,-.52553240991632899,-.18343464249564980,
        .18343464249564980,.52553240991632899,.79666647741362674,.96028985649753623};
    constexpr double weights[8]={.10122853629037626,.22238103445337447,.31370664587788729,.36268378337836198,
        .36268378337836198,.31370664587788729,.22238103445337447,.10122853629037626};
    long double sum=0;
    for(std::size_t c=1;c<cuts.size();++c){auto middle=atom.begin;(tangent==0?middle.x:middle.y)=(cuts[c-1]+cuts[c])/2;
        if(shape.exact(middle,0).v<0)continue;
        const double width=(cuts[c]-cuts[c-1])/8;
        for(int s=0;s<8;++s)for(int k=0;k<8;++k){auto point=atom.begin;
            (tangent==0?point.x:point.y)=cuts[c-1]+width*(s+(1+nodes[k])/2);
            const auto value=state.cellPhysicalVelocityJet(p,point,false).jet.value;
            sum+=width*weights[k]/2*component(value,atom.normalAxis)*(shape.radial?2*pi*point.x:1);}}
    return sum;
}
void tangentialSourceContract()
{
    for(int kind:{0,1})for(bool rz:{false,true})for(int amr:{0,1}){
        auto mesh=makeMesh(4,amr);Case shape;shape.pureFirstJet=true;shape.radial=rz;
        shape.planar=kind==0;shape.planeOffset=.43;shape.radius=rz?.34:.19;
        const auto geometry=CchMomentumGeometry::build(*mesh,metricsFor(*mesh,rz));
        const auto phase=CchPhaseGeometry::build(*mesh,geometry.metrics(),initial(*mesh,shape));
        const auto source=CchAdvectorReconstruction::prepare(*mesh,geometry).withTangentialModes(*mesh,phase);
        const CchNativeValues zero(geometry.layout(),std::vector<double>(geometry.layout().velocityCount()),
            std::vector<double>(geometry.layout().pressureCount()));
        check(!source.tangentialModes().empty(),"tangential source has no cut modes");
        rejects([&]{source.evaluate(*mesh,zero,phase.geometryIdentity(),phase.timeLayer());},"explicit tangential velocity state required");
        rejects([&]{source.evaluate(*mesh,zero,std::vector<double>{},phase.geometryIdentity(),phase.timeLayer());},"tangential state size mismatch");
        std::vector<double> initial(source.tangentialModes().size());
        for(std::size_t m=0;m<initial.size();++m)initial[m]=.125*(m+1);
        const auto state=source.evaluate(*mesh,zero,initial,phase.geometryIdentity(),phase.timeLayer());
        check(state.tangentialVelocity()==initial,"tangential initial state was not retained");
        rejects([&]{state.atomPolynomials();},"physical tangential phase flux coefficients required");
        double peakDerivative=0,peakDivergence=0,peakNormal=0,peakTangent=0;int axisModes=0;
        for(std::size_t m=0;m<initial.size();++m){const auto& mode=source.tangentialModes()[m];const auto h=mode.upper-mode.lower;
            check(source.tangentialModeForPressure(mode.pressure)==static_cast<int>(m),"modal pressure owner mismatch");
            const Vector2D point=mode.lower+Vector2D{.39*h.x,.61*h.y};
            const auto jet=source.tangentialModeJet(static_cast<int>(m),point,false);
            const auto physical=state.cellPhysicalVelocityJet(mode.pressure,point,false);
            close(physical.jet.value.x,initial[m]*jet.value.x,1e-12,"modal physical value x");
            close(physical.jet.value.y,initial[m]*jet.value.y,1e-12,"modal physical value y");
            check(state.cellVelocityJet(mode.pressure,point).jet.value.x==0&&state.cellVelocityJet(mode.pressure,point).jet.value.y==0,
                "zero-normal tangential mode changed geometry extension");
            for(int axis=0;axis<2;++axis){auto a=point,b=point;const double step=(axis==0?h.x:h.y)*1e-5;
                (axis==0?a.x:a.y)-=step;(axis==0?b.x:b.y)+=step;
                const auto lo=source.tangentialModeJet(static_cast<int>(m),a,false).value;
                const auto hi=source.tangentialModeJet(static_cast<int>(m),b,false).value;
                const auto expected=axis==0?Vector2D{jet.derivative.xx,jet.derivative.yx}:Vector2D{jet.derivative.xy,jet.derivative.yy};
                const auto measured=(hi-lo)*(1/(2*step));
                peakDerivative=std::max({peakDerivative,std::abs(measured.x-expected.x)/(1+std::abs(expected.x)),
                    std::abs(measured.y-expected.y)/(1+std::abs(expected.y))});}
            const double divergence=jet.derivative.xx+jet.derivative.yy+(rz?jet.value.x/point.x:0);
            peakDivergence=std::max(peakDivergence,std::abs(divergence)/(1+std::abs(jet.derivative.xx)+std::abs(jet.derivative.yy)));
            for(int side=0;side<4;++side){Vector2D boundary=mode.lower+Vector2D{.37*h.x,.63*h.y};
                (side<2?boundary.x:boundary.y)=component(side%2?mode.upper:mode.lower,side<2?0:1);
                const auto trace=source.tangentialModeJet(static_cast<int>(m),boundary,false);
                check(trace.value.x==0&&trace.value.y==0&&trace.derivative.xx==0&&trace.derivative.xy==0
                    &&trace.derivative.yx==0&&trace.derivative.yy==0,"modal native face value/jet leakage");}
            if(rz&&mode.lower.x==0){const auto axis=source.tangentialModeJet(static_cast<int>(m),{0,point.y},false);
                check(axis.value.x==0&&axis.derivative.xy==0&&axis.derivative.yx==0,"modal RZ axis parity");++axisModes;}
        }
        const auto pressure=CchPressureLoad::build(geometry.metrics().divergence());
        const auto load=CchInterfaceNormalLoad::capillary(*mesh,pressure,phase,1,source);
        check(load.tangentialBalanceLoads().size()==initial.size(),"normal stress modal row count");
        for(double force:load.tangentialBalanceLoads())check(force==0,"normal stress performed tangential modal work");
        const auto& publication=phase.publication();
        for(std::size_t c=0;c<publication.topology().curves.size();++c)
            for(std::size_t g=0;g<publication.topology().curves[c].graphs().size();++g)
                for(const auto& q:publication.quadrature(c,g)){const auto& frame=q.sample.frame;
                    const int p=geometry.layout().pressureDof(frame.nativeCellIndex);if(p<0)continue;
                    const int m=source.tangentialModeForPressure(p);if(m<0)continue;
                    const auto v=source.tangentialModeJet(m,frame.point,false).value;
                    const double normal=v.x*frame.normal.x+v.y*frame.normal.y;
                    peakNormal=std::max(peakNormal,std::abs(normal)/(1+std::hypot(v.x,v.y)));
                    peakTangent=std::max(peakTangent,std::abs(v.x*frame.normal.y-v.y*frame.normal.x));}
        check(peakDerivative<=2e-6&&peakDivergence<=2e-12&&peakNormal<=1e-10&&peakTangent>1e-4,
            "tangential first-jet/zero-normal/divergence identity");
        const auto rate=CchPhaseRate::prepareSpatial(*mesh,phase,state);
        const auto flux=CchMaterialFlux::prepare(*mesh,rate,1000,1);long double peakRate=0;
        for(std::size_t i=0;i<source.atoms().size();++i){const auto reference=independentModalAtomRate(shape,state,source.atoms()[i]);
            peakRate=std::max(peakRate,std::abs(rate.atomRates()[i].gas-reference)/(1+std::abs(reference)));
            check(rate.atomRates()[i].liquid==0&&flux.atomRates()[i]==rate.atomRates()[i].gas,"modal gas flux density/source changed");}
        for(const auto& face:rate.faceRates())check(face.liquid==0&&face.gas==0,"modal native face mean changed");
        std::vector<long double> outward(geometry.layout().velocityCount()),gross(outward.size());
        for(std::size_t f=0;f<geometry.segments().size();++f){const auto& segment=geometry.segments()[f];const auto flow=flux.dualFaceRates()[f];
            for(const auto [u,sign]:{std::pair<int,int>{segment.negativeVelocity,1},{segment.positiveVelocity,-1}})
                if(u>=0){outward[u]+=sign*flow;gross[u]+=std::abs(flow);}}
        for(std::size_t u=0;u<outward.size();++u)check(std::abs(outward[u])<=1e-12L*(1+gross[u]),"modal dual density flux curl closure");
        check(peakRate<=1e-12,"modal gas rate differs from independent physical trace integral");
        if(kind==1&&rz)check(axisModes>0,"missing modal axis fixture");
        std::cout<<"TANGENTIAL_SOURCE,"<<kind<<','<<rz<<','<<amr<<','<<initial.size()<<','<<axisModes<<','
            <<peakDerivative<<','<<peakDivergence<<','<<peakNormal<<','<<peakTangent<<','<<rate.spatialReport().tangentialAtoms
            <<','<<peakRate<<",PASS\n";
    }
}

template<class F>long double independentModeVolume(Vector2D lo,Vector2D hi,int subdivisions,F function)
{
    // Different partition from the dual/native/support cubature in production.
    constexpr long double nodes[8]={-.96028985649753623L,-.79666647741362674L,-.52553240991632899L,-.18343464249564980L,
        .18343464249564980L,.52553240991632899L,.79666647741362674L,.96028985649753623L};
    constexpr long double weights[8]={.10122853629037626L,.22238103445337447L,.31370664587788729L,.36268378337836198L,
        .36268378337836198L,.31370664587788729L,.22238103445337447L,.10122853629037626L};
    const long double dx=(static_cast<long double>(hi.x)-lo.x)/subdivisions,dy=(static_cast<long double>(hi.y)-lo.y)/subdivisions;
    long double sum=0;
    for(int x=0;x<subdivisions;++x)for(int y=0;y<subdivisions;++y)for(int i=0;i<8;++i)for(int j=0;j<8;++j){
        const Vector2D p{static_cast<double>(lo.x+dx*(x+(1+nodes[i])/2)),static_cast<double>(lo.y+dy*(y+(1+nodes[j])/2))};
        sum+=dx*dy*weights[i]*weights[j]/4*function(p);}
    return sum;
}
void tangentialKineticEvolutionContract()
{
    for(bool rz:{false,true})for(int amr:{0,1}){
        auto mesh=makeMesh(4,amr);Case shape;shape.pureFirstJet=true;shape.radial=rz;shape.planar=true;shape.planeOffset=.43;
        const auto geometry=CchMomentumGeometry::build(*mesh,metricsFor(*mesh,rz));const auto& layout=geometry.layout();
        const auto phase=CchPhaseGeometry::build(*mesh,geometry.metrics(),initial(*mesh,shape));
        const auto source=CchAdvectorReconstruction::prepare(*mesh,geometry).withTangentialModes(*mesh,phase);
        const auto hodge=CchKineticHodge::build(*mesh,CchMomentumMass::build(*mesh,geometry,phase,1000,1),source);
        const int nu=layout.velocityCount(),nm=static_cast<int>(source.tangentialModes().size()),nv=nu+nm;
        check(hodge.augmentedMatrix().rows()==nv&&hodge.report().tangentialModes==static_cast<std::size_t>(nm),"missing full modal Hodge blocks");
        std::vector<double> modal(nm),zero(nu),normal(nu);for(int m=0;m<nm;++m)modal[m]=.1*(m+1);
        const auto psi=[&](Vector2D p){const double radial=rz?p.x*p.x*(1-p.x*p.x)*(1-p.x*p.x):p.x*p.x*(1-p.x)*(1-p.x);
            return radial*p.y*p.y*(1-p.y)*(1-p.y);};
        for(int u=0;u<nu;++u){const auto& face=source.regularGeometry().faces()[u];if(face.area==0)continue;
            normal[u]=static_cast<double>((face.normalAxis==0?1:-1)*(rz?2*pi:1)*(psi(face.upper)-psi(face.lower))/face.area);}
        const CchNativeValues native(layout,normal,std::vector<double>(layout.pressureCount()));
        const CchNativeValues rest(layout,zero,std::vector<double>(layout.pressureCount()));
        const auto onlyMode=source.evaluate(*mesh,rest,modal,phase.geometryIdentity(),phase.timeLayer());
        const auto initial=source.evaluate(*mesh,native,modal,phase.geometryIdentity(),phase.timeLayer());
        const auto uniqueFlux=CchMaterialFlux::prepare(*mesh,CchPhaseRate::prepareSpatial(*mesh,phase,onlyMode),1000,1);
        const auto convection=CchKepConvection::modalBulkTransportMatrix(*mesh,hodge,uniqueFlux);
        const CchKineticHodge::Matrix symmetric=convection+CchKineticHodge::Matrix(convection.transpose());
        check(convection.norm()>0&&symmetric.norm()<=1e-12*convection.norm(),
            "closed curl carrier raw physical convection failed integration by parts");
        Eigen::VectorXd testVelocity=Eigen::VectorXd::Zero(nv),transported=Eigen::VectorXd::Zero(nv);
        for(int u=0;u<nu;++u){const auto& face=source.regularGeometry().faces()[u];if(face.normalAxis!=1)continue;
            const long double a=face.lower.x,b=face.upper.x;
            testVelocity[u]=static_cast<double>(rz?(a*a+b*b)/2:(a*a+a*b+b*b)/3);}
        for(int m=0;m<nm;++m)transported[nu+m]=modal[m];
        long double convectionReference=0;
        for(int m=0;m<nm;++m){const auto& mode=source.tangentialModes()[m];auto lo=mode.lower;
            (shape.axial()?lo.y:lo.x)=std::max(shape.axial()?lo.y:lo.x,shape.planeOffset);
            convectionReference-=independentModeVolume(lo,mode.upper,8,[&](Vector2D p){
                const auto v=source.tangentialModeJet(m,p,false).value*modal[m];
                return 2.L*p.x*v.x*v.y*(rz?2*pi*p.x:1);});}
        const double convectionWork=testVelocity.dot(convection*transported);
        check(std::abs(convectionWork-convectionReference)<=2e-10L*(1+std::abs(convectionReference)),
            "raw physical convection weak work differs from independent tensor integration");
        const long double energy=hodge.kineticEnergy(*mesh,onlyMode);long double reference[2]={0,0};
        std::vector<Vector2D> body;for(const auto& q:hodge.strainQuadrature())body.push_back({0,q.point.x});
        const auto force=hodge.bodyForceLoad(*mesh,body);double peakForceError=0;
        for(int m=0;m<nm;++m){const auto& mode=source.tangentialModes()[m];auto lo=mode.lower;
            (shape.axial()?lo.y:lo.x)=std::max(shape.axial()?lo.y:lo.x,shape.planeOffset);
            for(int level=0;level<2;++level)reference[level]+=.5L*modal[m]*modal[m]*independentModeVolume(lo,mode.upper,level?8:4,[&](Vector2D p){
                const auto v=source.tangentialModeJet(m,p,false).value;
                return (static_cast<long double>(v.x)*v.x+static_cast<long double>(v.y)*v.y)*(rz?2*pi*p.x:1);});
            const auto weak=independentModeVolume(lo,mode.upper,8,[&](Vector2D p){const double sx=rz&&mode.lower.x==0?
                    p.x*p.x/std::pow(mode.upper.x-mode.lower.x,2):(p.x-mode.lower.x)/(mode.upper.x-mode.lower.x),
                    sy=(p.y-mode.lower.y)/(mode.upper.y-mode.lower.y);
                const double bx=4*sx*(1-sx),by=4*sy*(1-sy),b=bx*bx*bx*by*by*by;
                const long double primitive=((shape.axial()?p.y:p.x)-shape.planeOffset)*b/mode.gradientScale*
                    (rz?p.x*p.x/(mode.upper.x-mode.lower.x):1);
                return primitive*(rz?2*pi:1);});
            peakForceError=std::max(peakForceError,static_cast<double>(std::abs(force[nu+m]-weak)/(1+std::abs(weak))));}
        check(energy>0&&std::abs(reference[1]-reference[0])<=1e-10L*energy&&std::abs(energy-reference[1])<=2e-9L*energy,
            "modal Hodge differs from independent physical cell integral");
        check(peakForceError<=1e-11,"modal body load failed independent streamfunction integration by parts");
        double cross=0;for(int u=0;u<nu;++u)for(CchKineticHodge::Matrix::InnerIterator it(hodge.augmentedMatrix(),u);it;++it)
            if(it.col()>=nu)cross=std::max(cross,std::abs(it.value()));
        check(cross>1e-10,"modal Hodge cross block disappeared");
        const auto pressure=CchPressureLoad::build(geometry.metrics().divergence());
        Eigen::VectorXd before(nv);for(int u=0;u<nu;++u)before[u]=normal[u];for(int m=0;m<nm;++m)before[nu+m]=modal[m];
        const double dt=.01;long double peakBalance=0;
        for(int scenario:{0,1,2,3}){
            const auto k=CchPhaseViscosity::modalBulkMatrix(*mesh,hodge,scenario==2?.1:0,scenario==2?.02:0);
            if(scenario!=2)check(k.nonZeros()==0,"exact mu=0 modal viscosity is not zero");
            CchSaddleVelocityBlock::Matrix a=hodge.augmentedMatrix()/dt+k*.5;
            Eigen::VectorXd rhs=hodge.augmentedMatrix()*before/dt-k*before*.5;
            if(scenario==3)rhs-=convection*before*.5;
            Eigen::VectorXd forcing=Eigen::VectorXd::Zero(nv);
            if(scenario==1){for(int i=0;i<nv;++i)forcing[i]=force[i];rhs+=forcing;}
            std::vector<double> nativeRhs(nu),modeRhs(nm);
            for(int u=0;u<nu;++u)nativeRhs[u]=rhs[u];for(int m=0;m<nm;++m)modeRhs[m]=rhs[nu+m];
            auto block=CchSaddleVelocityBlock::build(layout,a,nm);
            if(scenario==3)block=block.withSkew(CchSaddleVelocityBlock::Matrix(convection*.5));
            const auto saddle=CchSinglePhaseSaddle::buildClosed(pressure,block);
            const auto solved=saddle.solve(*mesh,CchSaddleRhs(layout,nativeRhs,zero,modeRhs));
            check(solved.report().freeTangentialCount==nm&&solved.tangentialVelocity().size()==modal.size(),"modal equation not solved");
            Eigen::VectorXd after(nv);for(int u=0;u<nu;++u)after[u]=solved.values().normalVelocity()[u];
            for(int m=0;m<nm;++m)after[nu+m]=solved.tangentialVelocity()[m];
            const auto state=source.evaluate(*mesh,solved.values(),solved.tangentialVelocity(),phase.geometryIdentity(),phase.timeLayer());
            const Eigen::VectorXd mid=(before+after)/2;
            const long double beginEnergy=hodge.kineticEnergy(*mesh,initial),endEnergy=hodge.kineticEnergy(*mesh,state);
            const long double work=dt*mid.dot(forcing),dissipation=dt*mid.dot(k*mid),balance=endEnergy-beginEnergy-work+dissipation;
            peakBalance=std::max(peakBalance,std::abs(balance)/(1+std::abs(beginEnergy)+std::abs(work)+std::abs(dissipation)));
            check(dissipation>=-1e-12&&std::abs(balance)<=2e-10L*(1+std::abs(beginEnergy)+std::abs(work)+std::abs(dissipation)),
                "full modal midpoint physical work/energy balance");
            if(scenario==0)check((after-before).cwiseAbs().maxCoeff()<=2e-10,"free inviscid modal initial state was not preserved");
            if(scenario==2)check(endEnergy<beginEnergy&&dissipation>0,"modal viscosity did not dissipate physical energy");
            if(scenario==3){check((after-before).norm()>1e-7,"physical convection did not evolve full modal velocity");
                check(std::abs(mid.dot(convection*mid))<=1e-11*(1+mid.squaredNorm()),"physical closed modal convection produced power");}
        }
        std::cout<<"TANGENTIAL_KINETIC_EVOLUTION,"<<rz<<','<<amr<<','<<nm<<','<<hodge.report().points<<','<<cross<<','
            <<static_cast<double>(std::abs(energy-reference[1])/energy)<<','<<peakForceError<<','<<peakBalance<<','
            <<symmetric.norm()/convection.norm()<<','<<convectionWork-convectionReference<<",PASS\n";
    }
}

void tangentialMovingMetricContract()
{
    for(bool rz:{false,true}){
        auto mesh=makeMesh(4);Case shape;shape.pureFirstJet=true;shape.radial=rz;shape.planar=true;shape.planeOffset=.43;
        shape.motion=1;shape.speed=.02;const double dt=.02;
        const auto geometry=CchMomentumGeometry::build(*mesh,metricsFor(*mesh,rz));const auto& layout=geometry.layout();
        const auto base=CchAdvectorReconstruction::prepare(*mesh,geometry);
        const auto at=[&](HermiteGalsTimeLayer clock){const auto f=shape.profile(clock.time());
            const auto seed=stateFor(*mesh,f,rz,clock,[](int,Vector2D,Vector2D){return 0.;},shape.symmetry(*mesh),HermiteGalsTransport::FirstJet);
            HermiteGalsCandidate candidate(*mesh,clock,seed.primary());const auto mixed=closeHermiteGalsMixedAmr(*mesh,candidate,shape.symmetry(*mesh));
            HermiteGalsPublicationConfig config;config.volume.symmetryPlanes=shape.symmetry(*mesh);const auto id=candidate.identity();
            return CchPhaseGeometry::build(*mesh,geometry.metrics(),publishHermiteGalsFirstJetGeometry(*mesh,std::move(candidate),mixed,{id,centers(*mesh,f)},config));};
        const auto g0=at({}),gm=at({1,0,dt,Position::Midpoint}),g1=at({1,0,dt,Position::End});
        const auto r0=base.withTangentialModes(*mesh,g0),rm=base.withTangentialModes(*mesh,gm),r1=base.withTangentialModes(*mesh,g1);
        const auto m0=CchMomentumMass::build(*mesh,geometry,g0,1000,1),mm=CchMomentumMass::build(*mesh,geometry,gm,1000,1),m1=CchMomentumMass::build(*mesh,geometry,g1,1000,1);
        const auto h0=CchKineticHodge::build(*mesh,m0,r0),hm=CchKineticHodge::build(*mesh,mm,rm),h1=CchKineticHodge::build(*mesh,m1,r1);
        const int nu=layout.velocityCount(),nm=static_cast<int>(rm.tangentialModes().size());
        std::vector<double> u0(nu),um(nu),u1(nu),q0(nm),qm(nm),q1(nm);
        for(int u=0;u<nu;++u){const auto& f=base.regularGeometry().faces()[u];if(f.normalAxis!=1)continue;
            const long double a=f.lower.x,b=f.upper.x;
            u0[u]=static_cast<double>(rz?(b*b+a*a)/2:(b*b+a*b+a*a)/3);u1[u]=1.05*u0[u];um[u]=(u0[u]+u1[u])/2;}
        for(int m=0;m<nm;++m){q0[m]=.2*(m+1);q1[m]=1.02*q0[m];qm[m]=(q0[m]+q1[m])/2;}
        const auto make=[&](const CchAdvectorReconstruction& r,const CchPhaseGeometry& g,const std::vector<double>& u,const std::vector<double>& q){
            return r.evaluate(*mesh,CchNativeValues(layout,u,std::vector<double>(layout.pressureCount())),q,g.geometryIdentity(),g.timeLayer());};
        const auto before=make(r0,g0,u0,q0),mid=make(rm,gm,um,qm),after=make(r1,g1,u1,q1);
        const auto transport=CchMaterialTransport::prepare(*mesh,g0,CchPhaseRate::prepareSpatial(*mesh,gm,mid),g1,1000,1);
        const auto material=CchDualMaterial::prepare(*mesh,transport,m0,mm,m1);
        const auto report=material.modalEnergyChange(*mesh,h0,hm,h1,before,after);
        const auto& change=report.kinetic;
        check(std::abs(change.identityResidual)<=2e-12L*(1+change.beginEnergy+change.endEnergy),"moving full modal metric polarization identity");
        check(std::abs(change.metricVariation)>1e-6&&material.report().components[1].l1IntegratedDefect>0,
            "moving modal test hid metric variation or material R");
        Eigen::VectorXd x(nu+nm);for(int u=0;u<nu;++u)x[u]=um[u];for(int m=0;m<nm;++m)x[nu+m]=qm[m];
        const CchKineticHodge::Matrix difference=h1.augmentedMatrix()-h0.augmentedMatrix();
        long double cross=0;
        for(int u=0;u<nu;++u)for(CchKineticHodge::Matrix::InnerIterator it(difference,u);it;++it)
            if(it.col()>=nu)cross+=static_cast<long double>(x[u])*it.value()*x[it.col()];
        check(std::abs(cross)>1e-10,"moving metric fixture failed to expose U-q cross change");
        std::cout<<"TANGENTIAL_MOVING_METRIC,"<<rz<<','<<nm<<','<<change.metricVariation<<','<<change.endpointQuadratureRemainder<<','
            <<cross<<','<<change.identityResidual<<','<<material.report().components[1].l1IntegratedDefect<<','
            <<report.legacyDiagonalMaterialChange[1]<<",PASS\n";
    }
}
void phaseCutKineticHodgeConstantEnergy(int selectedKind=-1,int selectedRz=-1,int selectedAmr=-1)
{
    for(int kind:{0,1})for(bool rz:{false,true})for(int amr:{0,1}){
        if((selectedKind>=0 && kind!=selectedKind)||(selectedRz>=0 && static_cast<int>(rz)!=selectedRz)
            ||(selectedAmr>=0 && amr!=selectedAmr))continue;
        auto mesh=makeMesh(4,amr);Case shape;shape.pureFirstJet=true;shape.radial=rz;
        shape.planar=kind==0;shape.planeOffset=.43;shape.radius=rz?.34:.19;
        const auto geometry=CchMomentumGeometry::build(*mesh,metricsFor(*mesh,rz));
        const auto phase=CchPhaseGeometry::build(*mesh,geometry.metrics(),initial(*mesh,shape));
        const auto mass=CchMomentumMass::build(*mesh,geometry,phase,1000,1);
        const auto source=CchAdvectorReconstruction::prepare(*mesh,geometry);
        const auto hodge=CchKineticHodge::build(*mesh,mass,source);
        hodge.validateFor(*mesh,phase.geometryIdentity(),phase.timeLayer());
        check(hodge.markedDuals().size()==static_cast<std::size_t>(geometry.layout().velocityCount())
            && hodge.report().markedDuals>0 && hodge.report().points>0 && hodge.report().nonZeros>0,
            "phase-cut kinetic Hodge missing marked quadrature");
        check(hodge.report().symmetryDefect<=1e-12 && hodge.report().maxRelativePhaseVolumeDifference<=1e-9,
            "phase-cut kinetic Hodge symmetry/phase measure");
        std::vector<double> velocity(geometry.layout().velocityCount(),0);
        for(int u=0;u<geometry.layout().velocityCount();++u)
            if(geometry.layout().normalAxis(u)==CchDofLayout::NormalAxis::Second)velocity[u]=1;
        const CchNativeValues values(geometry.layout(),velocity,std::vector<double>(geometry.layout().pressureCount()));
        const auto action=hodge.apply(*mesh,values,phase.geometryIdentity(),phase.timeLayer());
        long double matrixPair=0;
        for(int u=0;u<geometry.layout().velocityCount();++u){
            long double row=0;
            for(CchKineticHodge::Matrix::InnerIterator it(hodge.matrix(),u);it;++it)
                row+=static_cast<long double>(it.value())*velocity[it.col()];
            close(action[u],row,2e-12L,"kinetic Hodge apply differs from sparse matrix");
            matrixPair+=static_cast<long double>(velocity[u])*row;
        }
        const long double expected=.5L*mass.report().totalMass[1];
        const long double energy=hodge.kineticEnergy(*mesh,values,phase.geometryIdentity(),phase.timeLayer());
        check(expected>0 && std::abs(energy-expected)<=2e-9L*expected,
            "constant phase-cut field does not recover physical total mass");
        close(energy,matrixPair/2,2e-12L,"kinetic Hodge energy differs from independent quadratic form");
        std::cout<<"KINETIC_HODGE_CONSTANT_ENERGY,"<<kind<<','<<rz<<','<<amr<<','<<hodge.report().markedDuals<<','
            <<hodge.report().points<<','<<hodge.report().nonZeros<<','<<hodge.report().maxRelativePhaseVolumeDifference<<','
            <<energy<<",PASS\n";
    }
}
}
int main(int argc,char** argv)
{
    try {
        std::cout<<std::setprecision(17);
        if(argc==2 && std::string(argv[1])=="--kinetic-hodge-only"){
            phaseCutKineticHodgeConstantEnergy();return 0;
        }
        if(argc==2 && std::string(argv[1])=="--kinetic-hodge-circle-amr-only"){
            phaseCutKineticHodgeConstantEnergy(1,0,1);return 0;
        }
        if(argc==2 && std::string(argv[1])=="--phase-source-support-only"){
            phaseSourceSupportContract();return 0;
        }
        if(argc==2 && std::string(argv[1])=="--tangential-modes-only"){
            tangentialSourceContract();tangentialKineticEvolutionContract();tangentialMovingMetricContract();return 0;
        }
        check(argc==1,"unexpected H13 contract argument");
        for(bool rz:{false,true})for(int amr:{0,1,2})for(int selection:{0,1,2,3}) {
            auto mesh=makeMesh(4,amr);const auto geometry=CchMomentumGeometry::build(*mesh,metricsFor(*mesh,rz,selection));inspectGeometry(*mesh,geometry);
            std::cout<<"GEOMETRY,"<<rz<<','<<amr<<','<<selection<<','<<geometry.layout().pressureCount()<<','<<geometry.layout().velocityCount()<<','<<geometry.segments().size()<<'\n';
            if(amr>1 || selection>1)continue;
            for(int kind=0;kind<4;++kind) {
                Case s;s.pureFirstJet=true;s.radial=rz;s.radius=rz?.34:.19;s.planar=kind!=0;s.planeOffset=kind==1?.43:kind==2?2:-1;
                const auto phase=CchPhaseGeometry::build(*mesh,geometry.metrics(),initial(*mesh,s));
                for(double density:{1.,1000.,1e6})inspectMass(*mesh,s,CchMomentumMass::build(*mesh,geometry,phase,density,1));
            }
        }
        check(averageGap>1e-3,"fixture did not distinguish geometric dual mass from averaged cell density");std::cout<<"AVERAGE_GAP,"<<averageGap<<'\n';
        for(bool rz:{false,true})for(bool gas:{false,true}) {
            auto mesh=makeMesh(4);Case s;s.pureFirstJet=true;s.radial=rz;s.planar=true;s.planeOffset=gas?1-std::ldexp(1.,-40):std::ldexp(1.,-40);
            const auto geometry=CchMomentumGeometry::build(*mesh,metricsFor(*mesh,rz));
            const auto phase=CchPhaseGeometry::build(*mesh,geometry.metrics(),initial(*mesh,s));
            const auto mass=CchMomentumMass::build(*mesh,geometry,phase,gas?1:1e6,gas?1e6:1);inspectMass(*mesh,s,mass);int tiny=0;
            for(int u=0;u<mass.layout().velocityCount();++u){const auto phases=mass.phaseVolumes()[u];const auto small=gas?phases.gas:phases.liquid;
                if(small>0 && small<1e-10*geometry.volumes()[u].volume){check(std::isfinite(small*1e6) && small*1e6>0 && mass.masses()[u]!=(gas?phases.liquid:phases.gas),"lost tiny phase material contribution");++tiny;}}
            check(tiny>0,"missing tiny dual phase");std::cout<<"TINY,"<<rz<<','<<gas<<','<<tiny<<'\n';
        }
        declaredContacts();overflow();
        for(bool rz:{false,true}){for(int amr:{0,1})convergence(rz,amr);guards(rz);}
        for(int rz=0;rz<2;++rz){check(amrCutDuals[rz]>0,"missing cut dual at actual coarse-fine boundary");std::cout<<"AMR_CUT,"<<rz<<','<<amrCutDuals[rz]<<'\n';}
        std::cout<<"SUMMARY,"<<geometryScenarios<<','<<massScenarios<<','<<overlapChecks<<','<<sharedPairs<<','<<segmentChecks<<','<<rejected<<','<<spaceRows<<','<<rebuilds<<','<<peakGeometry<<','<<peakPhase<<','<<peakMass<<','<<peakPartition<<'\n';
        std::cout<<"H09a CCH native momentum geometry and same-Q4 material mass PASS\n";
        coincidentMass();
        nativeTraceMass();
        nativeTraceRegionExtras();
        spatialAdvectorContract();
        advectorGuards();
        spatialAdvectorJetContract();
        cellJetGuards();
        regularPreparationContract();
        normalFluxPortContract();
        localSourceStencilContract();
        completeNormalPolynomialContract();
        regularEvaluationContract();
        regularEvaluationGuards();
        regularLowerSpaceContract();
        regularIdentifiableConstantContract();
        matchingSpatialConsumers();
        phaseSourceSupportContract();
        tangentialSourceContract();tangentialKineticEvolutionContract();
        tangentialMovingMetricContract();
        phaseCutKineticHodgeConstantEnergy();
    } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
