#include "Modules/Hydrodynamics/HydroMomentumGeometry.h" // Public header stands alone.
#include "Modules/Hydrodynamics/HydroSpatialAdvector.h"
#include "Modules/Hydrodynamics/Detail/HydroRegularAdvectorGeometry.h"
#include "Modules/Hydrodynamics/HydroSinglePhaseSaddle.h"
#include "Modules/Hydrodynamics/HydroPhaseViscosity.h"
#include "Modules/Hydrodynamics/HydroInterfaceNormalLoad.h"
#include "Modules/Hydrodynamics/HydroMaterialTransport.h"
#include "Modules/Hydrodynamics/HydroDualMaterial.h"
#include "Modules/Hydrodynamics/HydroKepConvection.h"
#include "interface_first_jet_fixture.h"
#include "Modules/Interface/Detail/HermiteGalsRegionMeasure.h"
#include "Core/Geometry/Detail/PolynomialSegmentRoots.h"
#include "Core/Geometry/Detail/LevelSetHermiteQuarticLeafData.h"
#include "Core/Geometry/Detail/LevelSetMappedInterface.h"
#include "Core/Mesh/CartesianMeasure.h"
#include "Core/Mesh/NativeLeafOverlap.h"
#include "Core/Geometry/Detail/FiniteVolumeClosure.h"
#include <cerrno>
#include <cstring>
#if defined(_MSC_VER)
#include <float.h>
#endif
#include <map>
#include <tuple>
#include <Eigen/SVD>
#include <Eigen/LU>
#include <cfenv>
#include <chrono>
#include <exception>
#include <xmmintrin.h>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace
{
using namespace LAMP::Hydrodynamics;
int geometryScenarios=0,massScenarios=0,overlapChecks=0,sharedPairs=0,segmentChecks=0,rejected=0,spaceRows=0,rebuilds=0;
std::array<int,2> amrCutDuals{};
long double peakGeometry=0,peakPhase=0,peakMass=0,peakPartition=0,averageGap=0;
static_assert(!std::is_default_constructible_v<HydroMomentumGeometry>);
static_assert(!std::is_default_constructible_v<HydroMomentumMass>);
static_assert(std::is_same_v<decltype(std::declval<HydroMomentumGeometry&>().volumes()),const std::vector<HydroMomentumVolume>&>);
static_assert(std::is_same_v<decltype(std::declval<HydroMomentumMass&>().masses()),const std::vector<double>&>);
long double cap(Vector2D a,Vector2D b,bool rz)
{ return (static_cast<long double>(b.x)-a.x)*(static_cast<long double>(b.y)-a.y)*(rz?pi*(static_cast<long double>(a.x)+b.x):1.L); }
long double area(Vector2D a,Vector2D b,int normal,bool rz)
{ return (normal==0?static_cast<long double>(b.y)-a.y:static_cast<long double>(b.x)-a.x)*(rz?pi*(static_cast<long double>(a.x)+b.x):1.L); }
double component(Vector2D p,int axis) { return axis==0?p.x:p.y; }
int axis(const HydroDofLayout& l,int u) { return l.normalAxis(u)==HydroDofLayout::NormalAxis::First?0:1; }
bool pointEqual(Vector2D a,Vector2D b) { return a.x==b.x && a.y==b.y; }
long double norm(long double a,long double b,long double scale,long double limit,const std::string& why)
{
    const auto e=std::abs(a-b)/(scale>0?scale:1);
    check(std::isfinite(e) && e<=limit,why+" error="+std::to_string(static_cast<double>(e)));return e;
}
HydroPhysicalDivergence metricsFor(const QuadtreeMesh& mesh,bool rz,int selection=0)
{
    FaceSegmentGraph graph;graph.build(mesh,[selection](const QuadtreeMesh&,const Cell& c) {
        const auto p=c.centerCoord();
        if(selection==1)return p.x<.5;
        if(selection==2)return !(p.x>.25 && p.x<.75 && p.y>.25 && p.y<.75);
        if(selection==3)return p.x<.25 || p.x>.75;
        return true;
    });
    return HydroPhysicalDivergence::build(mesh,HydroIntegralDivergence::build(mesh,HydroDofLayout::build(mesh,NativeLeafDomain::build(mesh,graph,rz))));
}
// Independent ownership enumeration starts from a native face and its actual
// neighboring cells, not the production cell-side builder or reported overlaps.
struct ReferenceBox { Vector2D lower{0,0},upper{0,0}; std::vector<HydroMomentumOverlap> pieces; };
ReferenceBox reference(const QuadtreeMesh& mesh,const HydroDofLayout& layout,int u)
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
void inspectGeometry(const QuadtreeMesh& mesh,const HydroMomentumGeometry& geometry)
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
long double inspectMass(const QuadtreeMesh& mesh,const Case& s,const HydroMomentumMass& mass,bool exact=true)
{
    const auto& l=mass.layout();mass.validateFor(mesh,mass.geometryIdentity(),mass.timeLayer());
    check(mass.phaseGeometry().publication().state().primary().phaseVolumes.empty(),"invented material inventory");
    std::array<long double,2> totalMass{},volumeSum{},phaseSum{};
    std::array<std::vector<HydroPhaseMeasure>,2> cellSums;for(auto& a:cellSums)a.resize(l.pressureCount());long double error=0;
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
HydroMomentumMass buildMass(const QuadtreeMesh& mesh,const HydroMomentumGeometry& geometry,const Case& s,double rho=1000)
{ return HydroMomentumMass::build(mesh,geometry,HydroPhaseGeometry::build(mesh,geometry.metrics(),initial(mesh,s)),rho,1); }
void convergence(bool rz,int amr)
{
    double previous=0,order=0;
    for(int n:{4,8,16}) {
        auto mesh=makeMesh(n,amr);Case s;s.pureFirstJet=true;s.radial=rz;s.nonpolynomial=true;s.radius=rz?.34:.19;
        const auto geometry=HydroMomentumGeometry::build(*mesh,metricsFor(*mesh,rz));const auto mass=buildMass(*mesh,geometry,s);
        const auto e=static_cast<double>(inspectMass(*mesh,s,mass,false));check(e>0,"zero nonpolynomial reference error");
        if(n>4){check(e<previous,"dual mass error did not decrease");order=std::log2(previous/e);}previous=e;++spaceRows;
        std::cout<<"SPACE,"<<rz<<','<<amr<<','<<n<<','<<e<<'\n';
    }
    check(order>=2,"dual geometry mass order below gate");std::cout<<"ORDER,"<<rz<<','<<amr<<','<<order<<'\n';
}
template<class F>void reject(F f,const std::string& why){rejects(f,why);++rejected;}
void exactParallelMass(const HydroMomentumMass& a,const HydroMomentumMass& b)
{
    const auto ta=a.timeLayer(),tb=b.timeLayer();
    check(a.geometry().sameIdentity(b.geometry()) && a.layout().sameIdentity(b.layout())
        && a.geometryIdentity()==b.geometryIdentity() && ta.step==tb.step && ta.begin==tb.begin
        && ta.end==tb.end && ta.position==tb.position
        && &a.phaseGeometry().publication().state().geometry().snapshot()==&b.phaseGeometry().publication().state().geometry().snapshot(),
        "parallel mass detached geometry/layout/scalar/full clock");
    check(a.liquidDensity()==b.liquidDensity() && a.gasDensity()==b.gasDensity()
        && a.masses()==b.masses() && a.phaseVolumes().size()==b.phaseVolumes().size()
        && a.overlapPhaseVolumes().size()==b.overlapPhaseVolumes().size(),"parallel ordered masses changed");
    for(std::size_t u=0;u<a.phaseVolumes().size();++u) {
        const auto x=a.phaseVolumes()[u],y=b.phaseVolumes()[u];
        check(x.liquid==y.liquid && x.gas==y.gas,"parallel dual phase measure changed");
        const auto& p=a.overlapPhaseVolumes()[u];const auto& q=b.overlapPhaseVolumes()[u];
        check(p.size()==q.size(),"parallel native overlap ordering changed");
        for(std::size_t k=0;k<p.size();++k)
            check(p[k].liquid==q[k].liquid && p[k].gas==q[k].gas,"parallel native overlap measure changed");
    }
    const auto& x=a.report();const auto& y=b.report();
    for(int c=0;c<2;++c)check(x.totalPhaseVolumes[c].liquid==y.totalPhaseVolumes[c].liquid
        && x.totalPhaseVolumes[c].gas==y.totalPhaseVolumes[c].gas && x.totalMass[c]==y.totalMass[c],
        "parallel component totals changed the original addition order");
    check(x.maxRelativeCellPartition==y.maxRelativeCellPartition && x.maxRelativeLedgerDefect==y.maxRelativeLedgerDefect
        && x.maxRelativeMomentEstimate==y.maxRelativeMomentEstimate && x.cutOverlaps==y.cutOverlaps,
        "parallel mass report changed");
}
// Compare the local mutable accumulator with the untouched value-returning
// add, including components beyond inline storage and exact cancellation.
void exactRestrictionAccumulationWitness()
{
    namespace E=LAMP::Geometry::Detail::ExactExpansion;
    E::Expansion old,actual;std::size_t checks=0,peakComponents=0;
    const auto consume=[&](const E::Expansion& term) {
        const auto retained=term;
        old=E::add(old,term);E::accumulate(actual,term);
        check(actual.components()==old.components(),"restriction accumulation changed original expansion components");
        check(term.components()==retained.components(),"restriction accumulation changed its immutable term");
        peakComponents=std::max(peakComponents,actual.components().size());++checks;
    };
    consume(E::fromDouble(0));
    for(int k=0;k<6;++k)consume(E::fromDouble(std::ldexp(1.,-54*k)));
    check(peakComponents>4,"restriction accumulation did not exercise spilled expansion storage");
    const auto spilled=actual;
    const auto tail=E::multiply(E::fromDouble(1.+std::ldexp(1.,-27)),E::fromDouble(1.-std::ldexp(1.,-27)));
    check(tail.components().size()>1,"restriction accumulation witness lacks an actual product tail");
    consume(tail);consume(E::subtract(E::Expansion{},tail));
    consume(E::subtract(E::Expansion{},spilled));
    check(actual.components().empty(),"restriction accumulation did not preserve exact zero elimination");
    for(const auto& source:{E::fromDouble(.25),spilled}) {
        auto alias=source;const auto expected=E::add(source,source);
        E::accumulate(alias,alias);
        check(alias.components()==expected.components(),"restriction accumulation self-alias changed original components");++checks;
    }
    const auto maximum=E::fromDouble(std::numeric_limits<double>::max());
    rejects([&]{E::add(maximum,maximum);},"non-finite");
    rejects([&]{auto scratch=maximum;E::accumulate(scratch,maximum);},"non-finite");
    const auto term=E::fromDouble(.125);auto guarded=spilled;
    const auto rounding=std::fegetround();std::fesetround(FE_UPWARD);
    try {rejects([&]{E::accumulate(guarded,term);},"round-to-nearest");}
    catch(...){std::fesetround(rounding);throw;}std::fesetround(rounding);
    check(guarded.components()==spilled.components(),"restriction accumulation mutated scratch before FP rejection");
    const auto csr=_mm_getcsr();_mm_setcsr(csr|(1u<<15)|(1u<<6));
    try {rejects([&]{E::accumulate(guarded,term);},"gradual underflow");}
    catch(...){_mm_setcsr(csr);throw;}_mm_setcsr(csr);
    check(guarded.components()==spilled.components(),"restriction accumulation mutated scratch before FTZ rejection");
    const auto restored=E::add(spilled,term);E::accumulate(guarded,term);
    check(guarded.components()==restored.components(),"restriction accumulation did not recover after FP rejection");
    std::cout<<"EXACT_RESTRICTION_ACCUMULATION,"<<checks<<','<<peakComponents<<",PASS\n";
}
// Actual retained primary/closure inputs construct independent immutable
// owners for the filtered certificate and the original cold exact restriction.
void filteredTransverseSourceWitness(const QuadtreeMesh& mesh,const HermiteGalsPublication& publication,
    bool radial,std::size_t& successes,long double& peak)
{
    namespace E=LAMP::Geometry::Detail::ExactExpansion;
    namespace Q=LAMP::Geometry::Detail::LevelSetHermiteScalar;
    namespace I=LAMP::Geometry::Detail::MappedHermiteInterface;
    namespace P=LAMP::Geometry::Detail::PolynomialSegmentRoots;
    const auto& state=publication.state();const auto& primary=state.primary();const auto& closure=state.closure();
    const auto sourceId=state.identity();const auto* source=&state.geometry().snapshot();const auto clock=state.timeLayer();
    Patch patch;patch.activeCells.assign(mesh.cellIndexNum(),0);patch.nodes.resize(mesh.vertIndexNum());patch.cellCenters.resize(mesh.cellIndexNum());
    for(std::size_t k=0;k<primary.vertices.size();++k) {
        const auto& v=primary.vertices[k];patch.nodes[v.vertex]={v.phi,v.dx,v.dy,closure.vertexMixedDerivatives[k]};
    }
    for(const auto& edge:primary.edges)patch.edgeMidpoints.push_back({edge.vertices,edge.phi,edge.dn});
    for(std::size_t k=0;k<primary.geometryCells.size();++k) {
        patch.activeCells[primary.geometryCells[k]]=1;patch.cellCenters[primary.geometryCells[k]]=closure.cellBubbleCoefficients[k];
    }
    Q::QuarticLeafData filtered(mesh,patch,radial),original(mesh,patch,radial);
    std::array<bool,2> covered{};std::array<E::Interval,5> retained{};
    int retainedCell=-1,retainedAxis=-1,retainedSign=0;double retainedQ=0,retainedLower=0,retainedUpper=1,retainedParameter=0;
    for(const auto& chain:publication.topology().curves)for(const auto& graph:chain.graphs()) {
        if(graph.nativeSide || graph.nativeCellIndex<0)continue;
        const int axis=graph.coordinateIsX?0:1;if(covered[axis])continue;
        const double q=.5*component(graph.localLower,axis)+.5*component(graph.localUpper,axis);
        const double lo=component(graph.localLower,1-axis),hi=component(graph.localUpper,1-axis);
        std::array<E::Interval,5> controls;I::FilteredTransverseRoot actual;
        if(!filtered.transverseBernsteinControls(graph.nativeCellIndex,axis,{q,q},controls)
            || !I::filteredTransverseRoot(controls,lo,hi,actual))continue;
        check(filtered.identicalCell(graph.nativeCellIndex,original,graph.nativeCellIndex),
            "filtered transverse witness changed the actual exact native Q4");
        const Q::RationalSegmentCoordinate fixed{E::fromDouble(q),E::fromDouble(q),E::fromDouble(1.)};
        const Q::RationalSegmentCoordinate varying{E::fromDouble(lo),E::fromDouble(hi),E::fromDouble(1.)};
        const auto polynomial=original.restrictLocal(graph.nativeCellIndex,axis==0?fixed:varying,axis==0?varying:fixed);
        const auto reference=P::isolateEndpointRelative(polynomial.numerator);
        check(reference.roots.size()==1 && reference.enclosures.size()==1,
            "filtered transverse changed original complete exact root inventory");
        const double oldPoint=lo+(hi-lo)*reference.roots.front();const auto bracket=reference.enclosures.front();
        const double oldLower=std::nextafter(lo+(hi-lo)*bracket.lo,-std::numeric_limits<double>::infinity());
        const double oldUpper=std::nextafter(lo+(hi-lo)*bracket.hi,std::numeric_limits<double>::infinity());
        const double distance=std::min({actual.parameter,1.-actual.parameter,std::abs(actual.parameter-.5)});
        check(actual.lower>0 && actual.upper<1 && actual.lower<=actual.parameter && actual.parameter<=actual.upper
            && !(actual.lower<=.5 && actual.upper>=.5) && actual.upper-actual.lower<=1e-13*distance,
            "filtered transverse violates original native root width/distance gate");
        check(actual.lower<=oldUpper && oldLower<=actual.upper,
            "filtered transverse certified root misses independent exact root enclosure");
        const long double error=std::abs(static_cast<long double>(actual.parameter)-oldPoint);
        check(error<=1e-13*distance,"filtered transverse representative violates original root allowance");
        P::ExpandedQuarticPolynomial derivative{};
        for(std::size_t k=1;k<polynomial.numerator.size();++k)
            derivative[k-1]=E::multiply(polynomial.numerator[k],E::fromDouble(static_cast<double>(k)));
        check(actual.derivativeSign==E::sign(E::evaluate(derivative,reference.roots.front())),
            "filtered transverse derivative orientation differs from exact original Q4");
        I::FilteredTransverseJet jet;
        check(I::filteredTransverseGeometryJet(controls,actual.parameter,actual.derivativeSign,jet),
            "actual filtered geometry node did not produce its consumed jet certificate");
        check(std::isfinite(jet.value) && std::isfinite(jet.valueError) && jet.valueError>=0
            && std::isfinite(jet.derivative) && std::isfinite(jet.derivativeError) && jet.derivativeError>=0
            && std::abs(jet.derivative)>jet.derivativeError && (jet.derivative>0?1:-1)==actual.derivativeSign,
            "filtered geometry jet lost the original regular transverse derivative/roundoff gate");
        const auto parameter=E::fromDouble(actual.parameter);
        const Q::RationalSegmentCoordinate at{parameter,parameter,E::fromDouble(1.)};
        const auto value=original.restrictLocal(graph.nativeCellIndex,axis==0?fixed:at,axis==0?at:fixed);
        const auto physicalDerivative=original.restrictLocal(graph.nativeCellIndex,axis==0?fixed:at,axis==0?at:fixed,
            axis==0?0:1,axis==0?1:0);
        const auto bounds=original.leafBounds(graph.nativeCellIndex);
        const auto cellWidth=E::subtract(E::fromDouble(component(bounds[1],1-axis)),E::fromDouble(component(bounds[0],1-axis)));
        // Compare to exact numerator/positive-denominator authority directly;
        // do not round the independent source into a numerical jet oracle.
        const auto enclosed=[&](long double estimate,long double uncertainty,const auto& reference,const E::Expansion& factor) {
            const double lower=std::nextafter(static_cast<double>(estimate-uncertainty),-std::numeric_limits<double>::infinity());
            const double upper=std::nextafter(static_cast<double>(estimate+uncertainty),std::numeric_limits<double>::infinity());
            const auto numerator=E::multiply(reference.numerator[0],factor);
            check(E::sign(reference.denominator)>0 && E::sign(E::subtract(numerator,
                E::multiply(reference.denominator,E::fromDouble(lower))))>=0
                && E::sign(E::subtract(E::multiply(reference.denominator,E::fromDouble(upper)),numerator))>=0,
                "filtered geometry jet error bound misses the independent exact source value/residual or D1");
        };
        enclosed(jet.value,jet.valueError,value,E::fromDouble(1.));
        enclosed(jet.derivative,jet.derivativeError,physicalDerivative,cellWidth);
        check(E::sign(physicalDerivative.numerator[0])==actual.derivativeSign,
            "filtered geometry node D1 sign differs from independent exact point derivative");
        peak=std::max(peak,error/distance);covered[axis]=true;++successes;
        retained=controls;retainedCell=graph.nativeCellIndex;retainedAxis=axis;retainedQ=q;retainedLower=lo;retainedUpper=hi;
        retainedParameter=actual.parameter;retainedSign=actual.derivativeSign;
    }
    check(covered[0] && covered[1],"actual XY/RZ AMR filtered transverse did not succeed on both axes");
    const auto rounding=std::fegetround();std::fesetround(FE_UPWARD);
    try {
        rejects([&]{std::array<E::Interval,5> out;filtered.transverseBernsteinControls(retainedCell,retainedAxis,{retainedQ,retainedQ},out);},"round-to-nearest");
        rejects([&]{I::FilteredTransverseRoot out;I::filteredTransverseRoot(retained,retainedLower,retainedUpper,out);},"round-to-nearest");
        rejects([&]{I::FilteredTransverseJet out;I::filteredTransverseGeometryJet(retained,retainedParameter,retainedSign,out);},"round-to-nearest");
    } catch(...){std::fesetround(rounding);throw;}std::fesetround(rounding);
    const auto csr=_mm_getcsr();_mm_setcsr(csr|(1u<<15)|(1u<<6));
    try {
        rejects([&]{std::array<E::Interval,5> out;filtered.transverseBernsteinControls(retainedCell,retainedAxis,{retainedQ,retainedQ},out);},"gradual underflow");
        rejects([&]{I::FilteredTransverseRoot out;I::filteredTransverseRoot(retained,retainedLower,retainedUpper,out);},"gradual underflow");
        rejects([&]{I::FilteredTransverseJet out;I::filteredTransverseGeometryJet(retained,retainedParameter,retainedSign,out);},"gradual underflow");
    } catch(...){_mm_setcsr(csr);throw;}_mm_setcsr(csr);
    std::array<E::Interval,5> restored;I::FilteredTransverseRoot root;
    check(filtered.transverseBernsteinControls(retainedCell,retainedAxis,{retainedQ,retainedQ},restored)
        && I::filteredTransverseRoot(restored,retainedLower,retainedUpper,root),"filtered transverse did not recover after actual FP rejection");
    I::FilteredTransverseJet restoredJet;
    check(I::filteredTransverseGeometryJet(restored,root.parameter,root.derivativeSign,restoredJet),
        "filtered geometry jet did not recover after actual FP rejection");
    const auto restoredClock=state.timeLayer();
    check(state.identity()==sourceId && &state.geometry().snapshot()==source && restoredClock.step==clock.step
        && restoredClock.begin==clock.begin && restoredClock.end==clock.end && restoredClock.position==clock.position,
        "filtered transverse changed its original publication source/full clock");
}
std::size_t filteredTransverseFallbackWitness()
{
    namespace E=LAMP::Geometry::Detail::ExactExpansion;
    namespace I=LAMP::Geometry::Detail::MappedHermiteInterface;
    namespace P=LAMP::Geometry::Detail::PolynomialSegmentRoots;
    std::size_t fallbacks=0;
    for(int scene=0;scene<3;++scene) {
        const double root=scene==0?0.:scene==1?.5:std::ldexp(1.,-70);
        std::array<E::Interval,5> controls;
        for(int k=0;k<5;++k) {
            const double value=.25*k-root;controls[k]={value,value};
            if(scene==2)controls[k]={k==0?-2*root:std::nextafter(value,-std::numeric_limits<double>::infinity()),value};
        }
        if(scene==2)controls[0].hi=0.; // Outward near-endpoint uncertainty, never a guessed zero.
        I::FilteredTransverseRoot output;
        check(!I::filteredTransverseRoot(controls,0.,1.,output),"endpoint/half-root filter did not defer to exact authority");
        P::ExpandedQuarticPolynomial polynomial{};polynomial[0]=E::fromDouble(-root);polynomial[1]=E::fromDouble(1.);
        const auto reference=P::isolateEndpointRelative(polynomial);
        check(reference.roots.size()==1 && reference.enclosures.size()==1 && reference.enclosures[0].lo<=root
            && root<=reference.enclosures[0].hi,"filtered fallback lost the original exact endpoint/half inventory");++fallbacks;
    }
    const std::array<E::Interval,5> multiple{{{2.25,2.25},{-.75,-.75},{-1.75,-1.75},{-.75,-.75},{2.25,2.25}}};
    I::FilteredTransverseRoot output;
    check(!I::filteredTransverseRoot(multiple,0.,1.,output),"multiple transverse roots were falsely qualified");
    P::ExpandedQuarticPolynomial polynomial{};polynomial[0]=E::fromDouble(.1875);polynomial[1]=E::fromDouble(-1.);polynomial[2]=E::fromDouble(1.);
    const auto reference=P::isolateEndpointRelative(polynomial);
    check(reference.roots.size()==2 && reference.enclosures.size()==2 && reference.enclosures[0].lo<=.25
        && .25<=reference.enclosures[0].hi && reference.enclosures[1].lo<=.75 && .75<=reference.enclosures[1].hi,
        "multiple-root fallback changed original exact inventory");++fallbacks;
    const std::array<E::Interval,5> tangent{{{3.,3.},{0.,0.},{-1.,-1.},{0.,0.},{3.,3.}}};
    check(!I::filteredTransverseRoot(tangent,0.,1.,output),"near-critical transverse root was falsely qualified");
    I::FilteredTransverseJet rejectedJet{7.,11.,13.,17.};
    check(!I::filteredTransverseGeometryJet(tangent,.5,1,rejectedJet)
        && rejectedJet.value==7. && rejectedJet.valueError==11. && rejectedJet.derivative==13. && rejectedJet.derivativeError==17.,
        "uncertified transverse D1 did not preserve the fallback jet output");
    polynomial[0]=E::fromDouble(.25);rejects([&]{P::isolateEndpointRelative(polynomial);},"multiple");++fallbacks;
    return fallbacks;
}
// Direct independent engine calls do not share the source-owned initial
// root cache. The geometric all-channel witness below may legitimately reuse
// that same actual source certificate; independent process runs remain the
// full consumer comparison. This small check proves candidate engine work.
void guardedTransverseEngineWitness()
{
    namespace P=LAMP::Geometry::Detail::PolynomialSegmentRoots;
    namespace E=LAMP::Geometry::Detail::ExactExpansion;
    const auto polynomial=[](std::initializer_list<double> coefficients) {
        P::ExpandedQuarticPolynomial p{};std::size_t k=0;
        for(double value:coefficients)p[k++]=E::fromDouble(value);return p;
    };
    std::size_t qualified=0,suggestions=0,signedSteps=0,enclosed=0;long double peak=0;
    for(int scene=0;scene<2;++scene) {
        const double constant=scene==0?.3:std::ldexp(1.,-40);
        const auto p=scene==0?polynomial({-constant,1.,1.}):polynomial({-constant,1.});
        const long double exact=scene==0?2*static_cast<long double>(constant)
            /(1+std::sqrt(1+4*static_cast<long double>(constant))):constant;
        const auto old=P::isolateEndpointRelative(p);P::GuardedNewtonWork work;
        const auto actual=P::isolateEndpointRelativeGuardedNewton(p,&work);
        check(old.roots.size()==1 && old.enclosures.size()==1 && actual.roots.size()==1
            && actual.enclosures.size()==1,"guarded transverse engine changed complete simple-root inventory");
        const long double scale=std::min(exact,1-exact);
        const long double limit=128*std::numeric_limits<double>::epsilon()*scale;
        for(const auto* result:{&old,&actual}) {
            const auto bracket=result->enclosures.front();const long double point=result->roots.front();
            check(bracket.lo<=point && point<=bracket.hi && bracket.lo<=exact && exact<=bracket.hi,
                "guarded transverse engine analytic root outside original proven enclosure");
            check(std::abs(point-exact)<=limit,"guarded transverse engine original endpoint-relative root accuracy");
            peak=std::max(peak,std::abs(point-exact)/scale);
        }
        check(std::abs(static_cast<long double>(actual.roots.front())-old.roots.front())<=2*limit,
            "guarded transverse engine differs from independent original engine");
        check(work.qualifiedRoots==1 && work.suggestions>0 && work.signedSteps+work.enclosedRoots>0,
            "guarded transverse candidate was not actually used");
        qualified+=work.qualifiedRoots;suggestions+=work.suggestions;
        signedSteps+=work.signedSteps;enclosed+=work.enclosedRoots;
    }
    // A derivative hull touching zero is unqualified; no alternate stopping
    // allowance, endpoint snapping or guessed multiple root is introduced.
    const auto unqualified=polynomial({-.3,0.,1.});P::GuardedNewtonWork fallbackWork;
    const auto original=P::isolateEndpointRelative(unqualified);
    const auto fallback=P::isolateEndpointRelativeGuardedNewton(unqualified,&fallbackWork);
    check(fallbackWork.qualifiedRoots==0 && fallbackWork.suggestions==0 && fallback.roots==original.roots
        && fallback.enclosures.size()==original.enclosures.size(),"unqualified transverse engine did not retain original path");
    for(std::size_t k=0;k<fallback.enclosures.size();++k)
        check(fallback.enclosures[k].lo==original.enclosures[k].lo && fallback.enclosures[k].hi==original.enclosures[k].hi,
            "unqualified transverse engine changed original root enclosure");
    const auto multiple=polynomial({.25,-1.,1.});
    rejects([&]{P::isolateEndpointRelative(multiple);},"multiple");
    rejects([&]{P::isolateEndpointRelativeGuardedNewton(multiple);},"multiple");
    const auto zero=P::ExpandedQuarticPolynomial{};
    check(P::isolateEndpointRelativeGuardedNewton(zero).identicallyZero,
        "guarded transverse engine hid an identically-zero restriction");
    const auto valid=polynomial({-.3,1.,1.});
    const auto rounding=std::fegetround();std::fesetround(FE_UPWARD);
    try {rejects([&]{P::isolateEndpointRelativeGuardedNewton(valid);},"round-to-nearest");}
    catch(...){std::fesetround(rounding);throw;}std::fesetround(rounding);
    const auto csr=_mm_getcsr();_mm_setcsr(csr|(1u<<15)|(1u<<6));
    try {rejects([&]{P::isolateEndpointRelativeGuardedNewton(valid);},"gradual underflow");}
    catch(...){_mm_setcsr(csr);throw;}_mm_setcsr(csr);
    const auto restored=P::isolateEndpointRelativeGuardedNewton(valid);
    check(restored.roots.size()==1,"guarded transverse engine did not recover after actual FP rejection");
    std::cout<<"GUARDED_TRANSVERSE_ENGINE,"<<qualified<<','<<suggestions<<','<<signedSteps<<','<<enclosed<<','<<peak<<",PASS\n";
}
void independentRegionMassWitness(const QuadtreeMesh& mesh,const HydroMomentumMass& mass)
{
    namespace RM=LAMP::Interface::Detail::RegionMeasure;
    const auto& publication=mass.phaseGeometry().publication();const auto& state=publication.state();
    const auto view=std::make_shared<const HermiteGalsGeometryView>(state.geometry());
    const std::shared_ptr<const Snapshot> scalar(view,&view->snapshot());
    check(scalar.get()==&state.geometry().snapshot(),"batch witness uses a different retained scalar");
    const auto clock=state.timeLayer(),actualClock=mass.timeLayer();
    check(mass.geometryIdentity()==state.identity() && clock.step==actualClock.step
        && clock.begin==actualClock.begin && clock.end==actualClock.end && clock.position==actualClock.position,
        "batch witness detached source/full clock");
    std::vector<std::vector<LevelSetMappedInterface::Graph>> graphs(scalar->nativeCellCount());
    for(const auto& chain:publication.topology().curves){
        check(chain.usesSource(*scalar),"batch witness graph uses a different source");
        for(const auto& graph:chain.graphs())if(!graph.nativeSide)graphs.at(graph.nativeCellIndex).push_back(graph);}
    const auto& layout=mass.layout();const bool radial=layout.domain().axisymmetric();
    int independentCuts=0,cutHalfLeaves=0,overlaps=0,batchTargets=0;long double peak=0;
    // A direct successful batch call on ONE existing cut leaf proves this
    // fixture reaches the new Core API, even if a mass query could fall back.
    for(int native=0;native<scalar->nativeCellCount();++native)
        if(publication.topology().cellPhase.at(native)==0 && !graphs[native].empty()){
            const auto middle=mesh.cells(mesh.cellIndexToLabel[native]).centerCoord();
            const auto halves=RM::halfBoxes(*scalar,native,middle);
            const auto batch=RM::integrateHalfCells(scalar,halves,graphs[native],radial);
            for(std::size_t k=0;k<halves.size();++k){
                const auto& half=halves[k];const auto capacity=cap(half.lower,half.upper,radial);
                const auto single=RM::integrateCell(scalar,half,graphs[native],radial);
                const auto& actual=batch[k];
                check(actual.error>=0 && actual.error<=1e-10L*capacity
                    && single.error>=0 && single.error<=1e-10L*capacity,
                    "batch/single half-target moment estimate exceeds original gate");
                for(const auto* measured:{&actual,&single}){
                    check(std::max(std::abs(measured->fixedAxialAssemblyDefect),
                        std::abs(measured->fixedRadialAssemblyDefect))<=1e-10L*measured->assemblyCapacity,
                        "batch/single half-target fixed/coincident ledger");
                    for(const auto value:{measured->axialLiquid,measured->radialLiquid,
                        measured->axialGas,measured->radialGas})
                        check(value>=0 && value<=capacity,"batch/single half-target outside strict capacity");
                    peak=std::max({peak,
                        norm(measured->axialLiquid,measured->radialLiquid,capacity,1e-10,"batch/single half-target liquid ledger"),
                        norm(measured->axialGas,measured->radialGas,capacity,1e-10,"batch/single half-target gas ledger"),
                        norm(measured->axialLiquid+measured->axialGas,capacity,capacity,1e-10,"batch/single half-target axial closure"),
                        norm(measured->radialLiquid+measured->radialGas,capacity,capacity,1e-10,"batch/single half-target radial closure")});
                }
                peak=std::max({peak,
                    norm(actual.axialLiquid,single.axialLiquid,capacity,1e-10,"batch versus single raw axial liquid"),
                    norm(actual.radialLiquid,single.radialLiquid,capacity,1e-10,"batch versus single raw radial liquid"),
                    norm(actual.axialGas,single.axialGas,capacity,1e-10,"batch versus single raw axial gas"),
                    norm(actual.radialGas,single.radialGas,capacity,1e-10,"batch versus single raw radial gas")});
                ++batchTargets;
            }
            // Each orientation keeps its original two-half partition and
            // seven complete cuts. The old all-channel port is the oracle.
            using I=LevelSetMappedInterface;
            using Trace=LevelSetMappedSegment::NativeQuarticTrace;
            const auto fullBounds=scalar->nativeCellBounds(native);
            const auto fullCapacity=cap(fullBounds[0],fullBounds[1],radial);
            int selectedChecks=0;long double selectedPeak=0,fixedPeak=0;
            for(int group=0;group<2;++group) {
                std::vector<I::NativeMomentTarget> targets;
                std::vector<I::NativeSelectedMomentTarget> selectedTargets;
                for(int part=0;part<2;++part) {
                    const auto& requested=halves[2*group+part];I::NativeMomentTarget target;
                    for(int r=0;r<2;++r) {
                        const auto& half=halves[2*group+r];const auto a=half.lower,b=half.upper;
                        target.regions.push_back({r==part?0:1,{{a,{b.x,a.y},b,{a.x,b.y}}}});
                    }
                    // Actual original half-box gauge calculation, not a
                    // translated previously integrated global moment.
                    const auto center=(requested.lower+requested.upper)*.5;auto estimate=center;
                    const auto jet=scalar->scalarJetNative(estimate);const double n=jet.dx*jet.dx+jet.dy*jet.dy;
                    if(n>0 && std::isfinite(n)) {
                        estimate.x=std::clamp(estimate.x-jet.value*jet.dx/n,requested.lower.x,requested.upper.x);
                        estimate.y=std::clamp(estimate.y-jet.value*jet.dy/n,requested.lower.y,requested.upper.y);
                    }
                    target.greenOrigin={estimate.x<=center.x?requested.lower.x:requested.upper.x,
                        estimate.y<=center.y?requested.lower.y:requested.upper.y};
                    selectedTargets.push_back({target.regions,target.greenOrigin,{0}});targets.push_back(std::move(target));
                }
                const auto makeCuts=[&] {
                    std::vector<Trace> cuts;std::set<RM::CutKey> keys;
                    for(const auto& region:targets.front().regions)for(int e=0;e<4;++e) {
                        const auto key=RM::key(region.vertices[e],region.vertices[(e+1)%4]);
                        if(keys.insert(key).second)cuts.push_back(LevelSetMappedSegment::nativeQuarticTrace(
                            scalar,{key[0],key[1]},{key[2],key[3]},key[0]==key[2]?0:1));
                    }
                    check(cuts.size()==7,"selected pair changed original seven-cut partition");return cuts;
                };
                const auto retry=[&](auto invoke) {
                    auto cuts=makeCuts();
                    for(int attempt=0;attempt<=256;++attempt) {
                        std::vector<const Trace*> refs;for(const auto& cut:cuts)refs.push_back(&cut);
                        try{return invoke(refs);}
                        catch(const I::RootResolutionError& failure) {
                            check(attempt<256 && !failure.requests.empty(),"selected pair original face-root refinement limit");
                            std::set<std::size_t> replaced;
                            for(const auto& request:failure.requests) {
                                check(request.segmentIndex<cuts.size(),"selected pair foreign cut index");
                                if(replaced.insert(request.segmentIndex).second) {
                                    auto& cut=cuts[request.segmentIndex];
                                    check(request.rootIndex<cut.events().size(),"selected pair foreign root index");
                                    cut=cut.refined(request.rootIndex);
                                }
                            }
                        }
                    }
                    throw std::runtime_error("unreachable selected pair refinement exhaustion");
                };
                const auto selected=retry([&](const auto& refs){return I::integrateNativeQuarticSelectedRegionMomentsBatch(
                    *scalar,native,graphs[native],refs,selectedTargets,radial);});
                const auto original=retry([&](const auto& refs){return I::integrateNativeQuarticRegionMomentsBatch(
                    *scalar,native,graphs[native],refs,targets,radial);});
                check(selected.size()==2 && original.size()==2,"selected pair changed original target inventory");
                for(int part=0;part<2;++part) {
                    const auto& a=selected[part];const auto& old=original[part];
                    check(a.selectedRegions.size()==1 && a.selectedRegions.front().regionId==0
                        && a.fixedRegions.size()==old.regions.size(),"selected pair exposed an unconsumed free/phase region");
                    RM::Sum fixedAxial,fixedRadial,fixedCapacity;
                    for(std::size_t r=0;r<a.fixedRegions.size();++r) {
                        const auto& f=a.fixedRegions[r];const auto& b=old.regions[r];
                        check(f.regionId==b.regionId,"selected pair reordered complete fixed partition");
                        fixedPeak=std::max({fixedPeak,
                            norm(f.capacity,b.capacity,fullCapacity,1e-12,"selected fixed region original capacity"),
                            norm(f.fixedAxialLiquid,b.fixedAxialLiquid,fullCapacity,1e-12,"selected complete fixed axial"),
                            norm(f.fixedRadialLiquid,b.fixedRadialLiquid,fullCapacity,1e-12,"selected complete fixed radial")});
                        fixedAxial.add(f.fixedAxialLiquid);fixedRadial.add(f.fixedRadialLiquid);fixedCapacity.add(f.capacity);
                    }
                    const auto found=std::find_if(old.regions.begin(),old.regions.end(),[](const auto& r){return r.regionId==0;});
                    check(found!=old.regions.end(),"old all-channel oracle lacks actual selected owner");
                    const auto& x=a.selectedRegions.front();const auto& y=*found;
                    const auto halfCapacity=cap(halves[2*group+part].lower,halves[2*group+part].upper,radial);
                    const std::array<long double,10> xv{{x.capacity,x.axialLiquid,x.radialLiquid,x.axialGas,x.radialGas,
                        x.fixedAxialLiquid,x.fixedRadialLiquid,x.freeAxialLiquid,x.freeRadialLiquid,x.errorEstimate}},
                        yv{{y.capacity,y.axialLiquid,y.radialLiquid,y.axialGas,y.radialGas,
                        y.fixedAxialLiquid,y.fixedRadialLiquid,y.freeAxialLiquid,y.freeRadialLiquid,y.errorEstimate}};
                    for(std::size_t k=0;k<xv.size();++k)selectedPeak=std::max(selectedPeak,
                        norm(xv[k],yv[k],halfCapacity,1e-10,"selected versus original all-channel moment"));
                    for(const auto* m:{&x,&y}) {
                        check(m->errorEstimate>=0 && m->errorEstimate<=1e-10L*halfCapacity,
                            "selected pair original moment-estimate gate");
                        for(auto v:{m->axialLiquid,m->radialLiquid,m->axialGas,m->radialGas})
                            check(v>=0 && v<=halfCapacity,"selected pair original strict phase capacity");
                        norm(m->axialLiquid,m->radialLiquid,halfCapacity,1e-10,"selected pair liquid Green ledger");
                        norm(m->axialGas,m->radialGas,halfCapacity,1e-10,"selected pair gas Green ledger");
                        norm(m->axialLiquid+m->axialGas,halfCapacity,halfCapacity,1e-10,"selected pair axial capacity ledger");
                        norm(m->radialLiquid+m->radialGas,halfCapacity,halfCapacity,1e-10,"selected pair radial capacity ledger");
                    }
                    norm(fixedCapacity.value,fullCapacity,fullCapacity,1e-12,"selected complete fixed capacity");
                    norm(fixedAxial.value,a.exteriorFixedAxialLiquid+a.coincidentAxialLiquid,fullCapacity,1e-10,"selected fixed axial assembly");
                    norm(fixedRadial.value,a.exteriorFixedRadialLiquid+a.coincidentRadialLiquid,fullCapacity,1e-10,"selected fixed radial assembly");
                    fixedPeak=std::max({fixedPeak,
                        norm(a.exteriorFixedAxialLiquid,old.exteriorFixedAxialLiquid,fullCapacity,1e-12,"selected exterior axial"),
                        norm(a.exteriorFixedRadialLiquid,old.exteriorFixedRadialLiquid,fullCapacity,1e-12,"selected exterior radial"),
                        norm(a.coincidentAxialLiquid,old.coincidentAxialLiquid,fullCapacity,1e-12,"selected coincident axial"),
                        norm(a.coincidentRadialLiquid,old.coincidentRadialLiquid,fullCapacity,1e-12,"selected coincident radial")});
                    check(a.coincidentArcs==old.coincidentArcs && a.sourceRootIncidences==old.sourceRootIncidences,
                        "selected pair changed actual coincidence/root inventory");++selectedChecks;
                }
                if(group==0) {
                    auto cuts=makeCuts();std::vector<const Trace*> refs;for(const auto& cut:cuts)refs.push_back(&cut);
                    const auto rounding=std::fegetround();std::fesetround(FE_UPWARD);
                    try {rejects([&]{I::integrateNativeQuarticSelectedRegionMomentsBatch(
                        *scalar,native,graphs[native],refs,selectedTargets,radial);},"round-to-nearest");}
                    catch(...){std::fesetround(rounding);throw;}std::fesetround(rounding);
                }
            }
            check(selectedChecks==4,"selected pair witness missed an original half target");
            check(state.identity()==mass.geometryIdentity() && scalar.get()==&state.geometry().snapshot()
                && state.timeLayer().step==clock.step && state.timeLayer().begin==clock.begin
                && state.timeLayer().end==clock.end && state.timeLayer().position==clock.position,
                "selected pair changed actual source/full clock");
            std::cout<<"SELECTED_PAIR_OLD_ALL_CHANNEL,"<<radial<<','<<native<<','<<selectedChecks<<','<<selectedPeak
                <<','<<fixedPeak<<",PASS\n";
            break;
        }
    check(batchTargets==4,"batch witness did not execute four original half-targets");
    for(int u=0;u<layout.velocityCount();++u){
        const auto box=reference(mesh,layout,u);const auto capacity=cap(box.lower,box.upper,radial);
        check(box.pieces.size()==mass.overlapPhaseVolumes()[u].size(),"batch witness overlap count changed");
        RM::Sum liquidSum,gasSum;bool allNegative=true,allPositive=true;
        for(std::size_t k=0;k<box.pieces.size();++k){
            const auto& piece=box.pieces[k];const int native=layout.nativeCell(piece.pressure);
            const auto pieceCapacity=cap(piece.lower,piece.upper,radial);
            int phase=publication.topology().cellPhase.at(native);
            if(phase==0)phase=scalar->nativeQuarticBoxSign(native,piece.lower,piece.upper);
            RM::CellMeasure measured;
            if(phase!=0){
                check(phase==-1 || phase==1,"batch witness invalid phase certificate");
                measured.axialLiquid=measured.radialLiquid=phase<0?pieceCapacity:0;
                measured.axialGas=measured.radialGas=phase>0?pieceCapacity:0;
            }else{
                check(!graphs[native].empty(),"batch witness unresolved graph-free leaf");
                // This retained helper calls the ORIGINAL single-target Core
                // API with its own cuts/refinements and original box gauge.
                // It does not call region.measure or the batched mass path.
                measured=RM::integrateCell(scalar,{piece.lower,piece.upper,native},graphs[native],radial);
                ++independentCuts;
                const auto bounds=scalar->nativeCellBounds(native);
                const auto middle=mesh.cells(mesh.cellIndexToLabel[native]).centerCoord();
                const bool xHalf=piece.lower.y==bounds[0].y && piece.upper.y==bounds[1].y
                    && ((piece.lower.x==bounds[0].x && piece.upper.x==middle.x)
                        || (piece.lower.x==middle.x && piece.upper.x==bounds[1].x));
                const bool yHalf=piece.lower.x==bounds[0].x && piece.upper.x==bounds[1].x
                    && ((piece.lower.y==bounds[0].y && piece.upper.y==middle.y)
                        || (piece.lower.y==middle.y && piece.upper.y==bounds[1].y));
                if(xHalf || yHalf)++cutHalfLeaves;
                check(measured.error>=0 && measured.error<=1e-10L*pieceCapacity,
                    "independent single-target moment estimate exceeds original gate");
                check(std::max(std::abs(measured.fixedAxialAssemblyDefect),std::abs(measured.fixedRadialAssemblyDefect))
                    <=1e-10L*measured.assemblyCapacity,"independent single-target fixed/coincident ledger");
            }
            allNegative=allNegative && phase<0;allPositive=allPositive && phase>0;
            peak=std::max({peak,
                norm(measured.axialLiquid,measured.radialLiquid,pieceCapacity,1e-10,"single-target liquid Green ledger"),
                norm(measured.axialGas,measured.radialGas,pieceCapacity,1e-10,"single-target gas Green ledger"),
                norm(measured.axialLiquid+measured.axialGas,pieceCapacity,pieceCapacity,1e-10,"single-target axial capacity"),
                norm(measured.radialLiquid+measured.radialGas,pieceCapacity,pieceCapacity,1e-10,"single-target radial capacity")});
            const bool axial=piece.upper.y-piece.lower.y<=piece.upper.x-piece.lower.x;
            const auto liquid=axial?measured.axialLiquid:measured.radialLiquid,gas=axial?measured.axialGas:measured.radialGas;
            check(liquid>=0 && liquid<=pieceCapacity && gas>=0 && gas<=pieceCapacity,
                "independent single-target phase outside strict capacity");
            const auto& actualPieces=mass.geometry().volumes()[u].overlaps;
            const auto found=std::find_if(actualPieces.begin(),actualPieces.end(),[&](const auto& candidate){
                return candidate.pressure==piece.pressure && pointEqual(candidate.lower,piece.lower)
                    && pointEqual(candidate.upper,piece.upper);
            });
            check(found!=actualPieces.end(),"batch witness independent overlap has no actual owner/box match");
            const auto actualIndex=static_cast<std::size_t>(found-actualPieces.begin());
            const auto expected=RM::complement(pieceCapacity,liquid,gas);
            const auto actual=mass.overlapPhaseVolumes()[u][actualIndex];
            peak=std::max({peak,
                norm(actual.liquid,expected[0],pieceCapacity,1e-10,"actual batched overlap versus original single-target liquid"),
                norm(actual.gas,expected[1],pieceCapacity,1e-10,"actual batched overlap versus original single-target gas")});
            check((expected[0]==0 || actual.liquid>0) && (expected[1]==0 || actual.gas>0),
                "actual batched overlap erased independently resolved small phase");
            liquidSum.add(liquid);gasSum.add(gas);++overlaps;
        }
        auto expected=RM::complement(capacity,liquidSum.value,gasSum.value);
        if(allNegative)expected={capacity,0};else if(allPositive)expected={0,capacity};
        const auto actual=mass.phaseVolumes()[u];
        const long double expectedMass=mass.liquidDensity()*expected[0]+mass.gasDensity()*expected[1];
        peak=std::max({peak,
            norm(actual.liquid,expected[0],capacity,1e-10,"actual batched dual versus original single-target liquid"),
            norm(actual.gas,expected[1],capacity,1e-10,"actual batched dual versus original single-target gas"),
            norm(mass.masses()[u],expectedMass,(mass.liquidDensity()+mass.gasDensity())*capacity,1e-10,
                "actual batched mass versus original single-target mass")});
        norm(actual.liquid+actual.gas,capacity,capacity,1e-12,"actual batched dual capacity closure");
    }
    check(independentCuts>0 && cutHalfLeaves>0,"batch witness did not exercise cut native half-leaf targets");
    std::cout<<"BATCH_MASS_SINGLE_TARGET,"<<radial<<','<<layout.velocityCount()<<','<<overlaps<<','
        <<independentCuts<<','<<cutHalfLeaves<<','<<batchTargets<<','<<peak<<",PASS\n";
}
void retainedRegionBatchWitness(const QuadtreeMesh& mesh,const HydroMomentumMass& mass)
{
    const auto& phase=mass.phaseGeometry();const auto& publication=phase.publication();
    const auto& state=publication.state();HermiteGalsPhaseRegionConfig config;
    for(const auto& curve:publication.topology().curves){
        for(const auto& point:curve.boundaryContacts())
            if(std::none_of(config.boundaryContacts.begin(),config.boundaryContacts.end(),
                [point](Vector2D other){return pointEqual(point,other);}))config.boundaryContacts.push_back(point);
        for(const auto& plane:curve.symmetryPlanes())
            if(std::find(config.symmetryPlanes.begin(),config.symmetryPlanes.end(),plane)==config.symmetryPlanes.end())
                config.symmetryPlanes.push_back(plane);
    }
    const HermiteGalsPhaseRegionSource legacy(mesh,state,config),retained(mesh,publication,config);
    const auto clock=state.timeLayer(),actual=retained.source().timeLayer();
    check(retained.source().identity()==state.identity()
        && &retained.source().geometry().snapshot()==&state.geometry().snapshot()
        && actual.step==clock.step && actual.begin==clock.begin && actual.end==clock.end && actual.position==clock.position
        && retained.components()==legacy.components() && retained.nativeSideGraphs()==legacy.nativeSideGraphs(),
        "retained region detached actual scalar/full clock/topology");
    int cut=-1;
    for(int u=0;u<mass.layout().velocityCount();++u)
        if(mass.phaseVolumes()[u].liquid>0 && mass.phaseVolumes()[u].gas>0){cut=u;break;}
    check(cut>=0,"retained region witness lacks existing cut dual");
    std::vector<std::array<Vector2D,2>> boxes;
    for(int u:{0,cut}){const auto& volume=mass.geometry().volumes()[u];boxes.push_back({volume.lower,volume.upper});}
    const auto equal=[](const HermiteGalsPhaseRegionMeasure& a,const HermiteGalsPhaseRegionMeasure& b){
        check(a.sourceIdentity==b.sourceIdentity && pointEqual(a.lower,b.lower) && pointEqual(a.upper,b.upper)
            && std::tie(a.capacity,a.liquid,a.gas,a.rawAxialLiquid,a.rawRadialLiquid,a.rawAxialGas,a.rawRadialGas,
                a.momentErrorEstimate,a.maxRelativeLedgerDefect,a.intersectedLeaves,a.cutLeaves,a.integratedArcs,a.rootRefinements)
                ==std::tie(b.capacity,b.liquid,b.gas,b.rawAxialLiquid,b.rawRadialLiquid,b.rawAxialGas,b.rawRadialGas,
                    b.momentErrorEstimate,b.maxRelativeLedgerDefect,b.intersectedLeaves,b.cutLeaves,b.integratedArcs,b.rootRefinements)
            && a.leaves.size()==b.leaves.size(),"retained batch differs from original public region query");
        for(std::size_t k=0;k<a.leaves.size();++k){const auto& x=a.leaves[k];const auto& y=b.leaves[k];
            check(x.cell==y.cell && pointEqual(x.lower,y.lower) && pointEqual(x.upper,y.upper)
                && std::tie(x.capacity,x.liquid,x.gas)==std::tie(y.capacity,y.liquid,y.gas),
                "retained batch changed original leaf order/measure");}
    };
    const auto batch=retained.measureBatch(mesh,boxes);
    check(batch.size()==boxes.size(),"retained region batch changed query inventory");
    std::vector<HermiteGalsPhaseRegionMeasure> reference;
    for(std::size_t k=0;k<boxes.size();++k){
        check(!batch[k].failure,"valid retained region batch failed");
        reference.push_back(legacy.measure(mesh,boxes[k][0],boxes[k][1]));
        equal(reference.back(),retained.measure(mesh,boxes[k][0],boxes[k][1]));
        equal(reference.back(),batch[k].measure);
    }
    const std::array<Vector2D,2> invalid{boxes[0][0],Vector2D{boxes[0][0].x,boxes[0][1].y}};
    const std::array<Vector2D,2> outside{mesh.maxCoord+Vector2D{1,1},mesh.maxCoord+Vector2D{2,2}};
    const auto failed=retained.measureBatch(mesh,{boxes[0],invalid,outside,boxes[1]});
    check(failed.size()==4 && !failed[0].failure && failed[1].failure && failed[2].failure && !failed[3].failure,
        "retained region batch lost ordered success/failure slots");
    equal(reference[0],failed[0].measure);equal(reference[1],failed[3].measure);
    rejects([&]{for(const auto& slot:failed)if(slot.failure)std::rethrow_exception(slot.failure);},"invalid rectangle bounds");
    rejects([&]{retained.measureBatch(mesh,{invalid,boxes[0]});},"invalid rectangle bounds");
    auto foreign=makeMesh(2);
    rejects([&]{retained.measureBatch(*foreign,boxes);},"mesh");
    const auto rounding=std::fegetround();std::fesetround(FE_UPWARD);
    try {rejects([&]{retained.measureBatch(mesh,boxes);},"round-to-nearest");}
    catch(...){std::fesetround(rounding);throw;}std::fesetround(rounding);
    const auto warmed=retained.measureBatch(mesh,boxes);
    for(std::size_t k=0;k<boxes.size();++k){check(!warmed[k].failure,"restored FP batch failed");equal(reference[k],warmed[k].measure);}
    std::cout<<"RETAINED_REGION_BATCH,"<<mass.layout().domain().axisymmetric()<<','<<boxes.size()<<",EXACT,PASS\n";
}
void parallelMassContract(const QuadtreeMesh& mesh,const HydroMomentumGeometry& geometry,
    const HydroPhaseGeometry& phase,const Case& shape,bool requireFour=false)
{
#ifdef _OPENMP
    struct Restore {
        int threads=omp_get_max_threads(),dynamic=omp_get_dynamic();
        ~Restore(){omp_set_num_threads(threads);omp_set_dynamic(dynamic);}
    } restore;
    omp_set_dynamic(0);omp_set_num_threads(4);
    int team=0;
#pragma omp parallel
    {
#pragma omp single
        team=omp_get_num_threads();
    }
    check(team>=2 && (omp_get_thread_limit()<4 || team==4),"mass contract did not create the requested available team");
#else
    const int team=1;
#endif
    check(!requireFour || team==4,"four-thread mass contract requires an actual four-thread team");
    const auto& state=phase.publication().state();const auto primary=state.primary();const auto closure=state.closure();
    const auto id=state.identity();const auto clock=state.timeLayer();const auto* scalar=&state.geometry().snapshot();
    const auto sourceValues=[&] {
        std::vector<double> values;
        for(const auto& v:state.primary().vertices)values.insert(values.end(),{v.phi,v.dx,v.dy});
        for(const auto& e:state.primary().edges)values.insert(values.end(),{e.phi,e.dn});
        const auto& c=state.closure();values.insert(values.end(),c.vertexMixedDerivatives.begin(),c.vertexMixedDerivatives.end());
        values.insert(values.end(),c.cellBubbleCoefficients.begin(),c.cellBubbleCoefficients.end());return values;
    };
    const auto saved=sourceValues();
    // The first mass query on this newly constructed scalar uses the team.
    // Reusing its SAME authority afterward detects ordering/cache differences.
    const auto parallel=HydroMomentumMass::build(mesh,geometry,phase,1e6,1);
#ifdef _OPENMP
    omp_set_num_threads(1);
#endif
    const auto serial=HydroMomentumMass::build(mesh,geometry,phase,1e6,1);
    exactParallelMass(serial,parallel);inspectMass(mesh,shape,parallel);
    independentRegionMassWitness(mesh,parallel);
    retainedRegionBatchWitness(mesh,parallel);
#ifdef _OPENMP
    omp_set_num_threads(4);
    std::optional<HydroMomentumMass> nested;std::exception_ptr failure;
#pragma omp parallel num_threads(2)
    {
#pragma omp single
        {
            try {nested.emplace(HydroMomentumMass::build(mesh,geometry,phase,1e6,1));}
            catch(...) {failure=std::current_exception();}
        }
    }
    if(failure)std::rethrow_exception(failure);
    check(nested.has_value(),"nested caller did not finish mass preparation");exactParallelMass(serial,*nested);
#endif
    const auto rounding=std::fegetround();std::fesetround(FE_UPWARD);
    try {reject([&]{HydroMomentumMass::build(mesh,geometry,phase,1e6,1);},"round-to-nearest");}
    catch(...){std::fesetround(rounding);throw;}std::fesetround(rounding);
    const auto csr=_mm_getcsr();_mm_setcsr(csr|(1u<<15)|(1u<<6));
    try {reject([&]{HydroMomentumMass::build(mesh,geometry,phase,1e6,1);},"gradual underflow");}
    catch(...){_mm_setcsr(csr);throw;}_mm_setcsr(csr);
    const auto after=state.timeLayer();
    check(sourceValues()==saved && state.identity()==id && &state.geometry().snapshot()==scalar
        && after.step==clock.step && after.begin==clock.begin && after.end==clock.end && after.position==clock.position
        && state.primary().measure==primary.measure && state.primary().transport==primary.transport
        && state.primary().geometryCells==primary.geometryCells && state.primary().phaseVolumes.empty()
        && state.primary().vertices.size()==primary.vertices.size() && state.primary().edges.size()==primary.edges.size()
        && state.closure().vertexMixedDerivatives==closure.vertexMixedDerivatives
        && state.closure().cellBubbleCoefficients==closure.cellBubbleCoefficients,"parallel/failed mass changed retained FirstJet source");
    for(std::size_t i=0;i<primary.vertices.size();++i)check(state.primary().vertices[i].vertex==primary.vertices[i].vertex,"parallel source vertex association changed");
    for(std::size_t i=0;i<primary.edges.size();++i)check(state.primary().edges[i].vertices==primary.edges[i].vertices,"parallel source edge association changed");
    exactParallelMass(serial,parallel);
    std::cout<<"PARALLEL_MASS,"<<shape.radial<<','<<team<<','<<geometry.layout().velocityCount()<<",EXACT,PASS\n";
}
void declaredContacts()
{
    auto mesh=makeMesh(4,1);const auto geometry=HydroMomentumGeometry::build(*mesh,metricsFor(*mesh,false));
    const Profile f=[](Vector2D p){return Jet{p.x+p.y-.7,1,1,0};};
    const auto seed=stateFor(*mesh,f,false,{},[](int,Vector2D,Vector2D){return 0.;},{},HermiteGalsTransport::FirstJet);
    HermiteGalsCandidate candidate(*mesh,{},seed.primary());const auto mixed=closeHermiteGalsMixedAmr(*mesh,candidate);
    HermiteGalsPublicationConfig config;config.volume.boundaryContacts={{0,.7},{.7,0}};
    const auto id=candidate.identity();const auto pub=publishHermiteGalsFirstJetGeometry(*mesh,std::move(candidate),mixed,{id,centers(*mesh,f)},config);
    const auto mass=HydroMomentumMass::build(*mesh,geometry,HydroPhaseGeometry::build(*mesh,geometry.metrics(),pub),1e6,1);
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
    const auto geometry=HydroMomentumGeometry::build(mesh,metricsFor(mesh,false));
    const auto phase=HydroPhaseGeometry::build(mesh,geometry.metrics(),initial(mesh,s));
    const auto failure=[&] {
        std::string error;
        reject([&]{try {HydroMomentumMass::build(mesh,geometry,phase,std::numeric_limits<double>::max(),1);}
            catch(const std::exception& e){error=e.what();throw;}},
            std::numeric_limits<long double>::digits==std::numeric_limits<double>::digits?"material product overflow or underflow":"unrepresentable momentum mass");
        return error;
    };
#ifdef _OPENMP
    struct Restore {int count=omp_get_max_threads();~Restore(){omp_set_num_threads(count);}} restore;
    omp_set_num_threads(1);
#endif
    const auto serial=failure();
#ifdef _OPENMP
    omp_set_num_threads(4);
#endif
    check(serial==failure(),"parallel mass changed the original first failing exception");
    check(phase.publication().state().primary().phaseVolumes.empty(),"failed parallel mass invented source inventory");
}
void guards(bool rz)
{
    auto mesh=makeMesh(4);Case s;s.pureFirstJet=true;s.radial=rz;s.planar=true;s.planeOffset=.43;
    const auto metrics=metricsFor(*mesh,rz);const auto geometry=HydroMomentumGeometry::build(*mesh,metrics);
    const auto phase=HydroPhaseGeometry::build(*mesh,metrics,initial(*mesh,s));const auto mass=HydroMomentumMass::build(*mesh,geometry,phase,1000,1);const auto saved=mass.masses();
    for(double bad:{0.,-1.,std::numeric_limits<double>::quiet_NaN(),std::numeric_limits<double>::infinity()}) {
        reject([&]{HydroMomentumMass::build(*mesh,geometry,phase,bad,1);},"invalid material density");reject([&]{HydroMomentumMass::build(*mesh,geometry,phase,1,bad);},"invalid material density");
    }
    const auto other=HydroMomentumGeometry::build(*mesh,metricsFor(*mesh,rz));
    reject([&]{HydroMomentumMass::build(*mesh,other,phase,1,1);},"different layout identity");
    const auto otherDi=HydroMomentumGeometry::build(*mesh,HydroPhysicalDivergence::build(*mesh,HydroIntegralDivergence::build(*mesh,metrics.layout())));
    reject([&]{HydroMomentumMass::build(*mesh,otherDi,phase,1,1);},"different retained D_I");
    reject([&]{mass.validateFor(*mesh,{},mass.timeLayer());},"different geometry identity");
    for(auto wrong:{HermiteGalsTimeLayer{1,0,0,Position::Initial},HermiteGalsTimeLayer{1,-.01,.01,Position::Midpoint},HermiteGalsTimeLayer{1,0,.1,Position::End}})
        reject([&]{mass.validateFor(*mesh,mass.geometryIdentity(),wrong);},"different full time layer");
    reject([&]{HydroMomentumMass::build(*mesh,geometry,phase,std::numeric_limits<double>::denorm_min(),std::numeric_limits<double>::denorm_min());},
        std::numeric_limits<long double>::digits==std::numeric_limits<double>::digits?"material product overflow or underflow":"unrepresentable momentum mass");
    check(mass.masses()==saved,"failed build mutated retained masses");
    check(geometry.sameIdentity(mass.geometry()) && &mass.geometry().metrics().divergence().matrix()==&metrics.divergence().matrix(),"detached geometry or D_I");
    const auto copied=mass;check(&copied.masses()==&mass.masses() && &copied.geometry().segments()==&geometry.segments(),"copy rebuilt objects");
    auto equivalent=makeMesh(4);mesh.reset();copied.validateFor(*equivalent,copied.geometryIdentity(),copied.timeLayer());inspectMass(*equivalent,s,copied);
    auto& vertex=equivalent->verts(equivalent->vertIndexToLabel[0]);const auto old=vertex.coord();vertex.setCoord(old+Vector2D{.01,0});
    reject([&]{copied.validateFor(*equivalent,copied.geometryIdentity(),copied.timeLayer());},"NativeMeshIdentity");vertex.setCoord(old);
    DictionaryMesh::RefineBox box;box.level=1;box.min={.25,.25};box.max={.75,.75};equivalent->applyRefineBoxes({box});
    reject([&]{copied.validateFor(*equivalent,copied.geometryIdentity(),copied.timeLayer());},"NativeMeshIdentity");
    reject([&]{HydroMomentumMass::build(*equivalent,geometry,phase,1,1);},"NativeMeshIdentity");
    const auto rebuilt=HydroMomentumGeometry::build(*equivalent,metricsFor(*equivalent,rz));
    check(!rebuilt.sameIdentity(geometry),"regrid kept old geometry identity");inspectGeometry(*equivalent,rebuilt);inspectMass(*equivalent,s,buildMass(*equivalent,rebuilt,s));++rebuilds;
    std::cout<<"REGRID,"<<rz<<','<<rebuilt.layout().pressureCount()<<','<<rebuilt.layout().velocityCount()<<'\n';
}
}
void coincidentMass(){
    int scenarios=0,rows=0,overlaps=0;long double peak=0;
    for(bool rz:{false,true})for(int amr:{0,1})for(int sign:{-1,1})for(int selection:{0,1}){
        auto mesh=makeMesh(4,amr);Case s;s.pureFirstJet=true;s.radial=rz;s.planar=true;
        s.planeOffset=amr?.125:.375;s.planeSign=sign;
        const auto geometry=HydroMomentumGeometry::build(*mesh,metricsFor(*mesh,rz,selection));
        const auto pub=initial(*mesh,s);const auto phase=HydroPhaseGeometry::build(*mesh,geometry.metrics(),pub);
        const auto mass=HydroMomentumMass::build(*mesh,geometry,phase,1000,1);
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
            const auto geometry=HydroMomentumGeometry::build(*mesh,metricsFor(*mesh,rz,selection));const auto phase=HydroPhaseGeometry::build(*mesh,geometry.metrics(),pub);
            for(double density:{1.,1000.,1e6}) {
                const auto mass=HydroMomentumMass::build(*mesh,geometry,phase,density,1);
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
            const auto geometry=HydroMomentumGeometry::build(*mesh,metricsFor(*mesh,rz));
            const auto mass=HydroMomentumMass::build(*mesh,geometry,HydroPhaseGeometry::build(*mesh,geometry.metrics(),pub),1000,1);
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
static_assert(!std::is_default_constructible_v<HydroAdvectorReconstruction>);
static_assert(!std::is_default_constructible_v<HydroSpatialAdvector>);
static_assert(!std::is_constructible_v<HydroAdvectorIdentity,std::uint64_t>);
// Physical divergence density: polynomial through degree five on each actual
// piece for both supported sources. GL3 is exact; the work scale retains each
// signed derivative contribution before cancellation. No volume-coordinate fit.
using PhysicalCuts=std::array<std::vector<double>,2>;
std::vector<PhysicalCuts> physicalPartitions(const HydroAdvectorReconstruction& r)
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
    const QuadtreeMesh& mesh;const HydroDofLayout& l;bool rz;int kind;
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
    AdvectorReference(const QuadtreeMesh& m,const HydroDofLayout& layout,int profile):mesh(m),l(layout),rz(l.domain().axisymmetric()),kind(profile),faces(l.pressureCount()) {
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
        for(int selection=0;selection<4;++selection){const auto geometry=HydroMomentumGeometry::build(*mesh,metricsFor(*mesh,rz,selection));const auto begin=std::chrono::steady_clock::now();
            const auto reconstruction=HydroAdvectorReconstruction::prepare(*mesh,geometry);const auto copied=reconstruction;check(copied.sameIdentity(reconstruction) && geometry.sameIdentity(reconstruction.momentumGeometry()),"prepared geometry authority");
            const auto prepared=std::chrono::steady_clock::now();const auto& l=geometry.layout();const auto& report=reconstruction.report();
            std::size_t actualViews=0;for(const auto& v:reconstruction.primalAtoms())actualViews+=v.size();for(const auto& v:reconstruction.dualAtoms())actualViews+=v.size();
            check(report.atoms==reconstruction.atoms().size()&&report.viewEntries==actualViews,"actual prepared geometry inventory");
            deficient+=static_cast<int>(report.rankDeficient);missing+=static_cast<int>(report.tangentUnidentified);
            const auto partitions=physicalPartitions(reconstruction);
            const bool exactFamily=n>=8||selection<2; // Frozen original 56/64 analytical eligibility.
            std::optional<HydroSpatialAdvector> first;
            for(int kind:{0,1}){AdvectorReference reference(*mesh,l,kind);const HydroNativeValues values(l,reference.vel,std::vector<double>(l.pressureCount()));
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
    for(bool rz:{false,true}){auto mesh=makeMesh(4,1);const auto geometry=HydroMomentumGeometry::build(*mesh,metricsFor(*mesh,rz));const auto rec=HydroAdvectorReconstruction::prepare(*mesh,geometry);
        Case shape;shape.pureFirstJet=true;shape.radial=rz;shape.planar=true;shape.planeOffset=2;const auto pub=initial(*mesh,shape);AdvectorReference ref(*mesh,geometry.layout(),0);const auto clock=pub.state().timeLayer();const auto id=pub.state().identity();
        const HydroNativeValues values(geometry.layout(),ref.vel,std::vector<double>(geometry.layout().pressureCount()));const auto a=rec.evaluate(*mesh,values,id,clock);const auto copy=a;
        const auto atomBefore=a.atomPolynomials()[0].beginCoefficients;
        const HydroNativeValues zeros(geometry.layout(),std::vector<double>(ref.vel.size()),std::vector<double>(geometry.layout().pressureCount()));const auto stationary=rec.evaluate(*mesh,zeros,id,clock);
        for(const auto& polynomial:stationary.atomPolynomials())check(std::all_of(polynomial.beginCoefficients.begin(),polynomial.beginCoefficients.end(),[](auto c){return c==0;}) && polynomial.integral(false,0,polynomial.length)==0,"zero native input produced moving trace");
        int zeroMeanFace=0;while(axis(geometry.layout(),zeroMeanFace)!=1)++zeroMeanFace;auto signedValues=ref.vel;const double target=ref.mean[zeroMeanFace].x;
        for(std::size_t u=0;u<signedValues.size();++u)signedValues[u]=axis(geometry.layout(),static_cast<int>(u))==1?ref.mean[u].x-target:0;
        const auto crossing=rec.evaluate(*mesh,HydroNativeValues(geometry.layout(),signedValues,std::vector<double>(geometry.layout().pressureCount())),id,clock);
        long double net=0,gross=0;bool nonconstant=false;for(int i:rec.primalAtoms()[zeroMeanFace]){const auto& p=crossing.atomPolynomials()[i];net+=p.integral(false,0,p.length);gross+=p.length*(std::abs(p.value(false,0))+std::abs(p.value(true,0)));nonconstant=nonconstant || std::any_of(p.beginCoefficients.begin()+1,p.beginCoefficients.end(),[](auto c){return c!=0;});}
        check(signedValues[zeroMeanFace]==0 && nonconstant && gross>0,"zero face mean erased sign-changing normal polynomial");norm(net,0,gross,1e-12,"zero-net spatial trace");
        rejects([&]{rec.evaluate(*mesh,values,{},clock);},"unbound provenance");auto wrong=clock;wrong.step=-1;rejects([&]{rec.evaluate(*mesh,values,id,wrong);},"invalid step");
        wrong=clock;++wrong.step;rejects([&]{a.validateFor(*mesh,geometry,id,wrong);},"full time layer");wrong={7,.2,.3,Position::Midpoint};rejects([&]{a.validateFor(*mesh,geometry,id,wrong);},"full time layer");
        for(auto t:std::vector<HermiteGalsTimeLayer>{{7,.2,.3,Position::Midpoint},{7,.2,.3,Position::End}}){const auto revision=rec.evaluate(*mesh,values,id,t);revision.validateFor(*mesh,geometry,id,t);check(revision.identity()!=a.identity(),"clock revision identity");auto changed=t;changed.position=t.position==Position::End?Position::Midpoint:Position::End;rejects([&]{revision.validateFor(*mesh,geometry,id,changed);},"full time layer");}
        rejects([&]{a.validateFor(*mesh,HydroMomentumGeometry::build(*mesh,geometry.metrics()),id,clock);},"momentum geometry instance");rejects([&]{a.validateFor(*mesh,geometry,{},clock);},"geometry identity");
        const auto fresh=HydroMomentumGeometry::build(*mesh,metricsFor(*mesh,rz));rejects([&]{rec.evaluate(*mesh,HydroNativeValues(fresh.layout(),ref.vel,std::vector<double>(fresh.layout().pressureCount())),id,clock);},"layout identity");
        auto bad=ref.vel;bad[0]=std::numeric_limits<double>::infinity();rejects([&]{HydroNativeValues(geometry.layout(),bad,std::vector<double>(geometry.layout().pressureCount()));},"nonfinite native");
        bad=ref.vel;for(std::size_t u=0;u<bad.size();++u)bad[u]=(rz && geometry.metrics().faceMeasures()[u]==0)?0:std::numeric_limits<double>::max();rejects([&]{rec.evaluate(*mesh,HydroNativeValues(geometry.layout(),bad,std::vector<double>(geometry.layout().pressureCount())),id,clock);},"nonfinite");
        if(rz){bad=ref.vel;bad[geometry.metrics().axisVelocityDofs()[0]]=1;rejects([&]{rec.evaluate(*mesh,HydroNativeValues(geometry.layout(),bad,std::vector<double>(geometry.layout().pressureCount())),id,clock);},"nonzero radial axis");}
        const auto rounding=std::fegetround();std::fesetround(FE_UPWARD);try{rejects([&]{rec.evaluate(*mesh,values,id,clock);},"round-to-nearest");}catch(...){std::fesetround(rounding);throw;}std::fesetround(rounding);
        const auto control=_mm_getcsr();_mm_setcsr(control|(1U<<15));try{rejects([&]{a.nativeVelocity(0);},"gradual underflow");}catch(...){_mm_setcsr(control);throw;}_mm_setcsr(control);
        rejects([&]{a.cellVelocity(0,{-1,0});},"point outside");rejects([&]{a.atomPolynomials()[0].integral(false,0,-1);},"invalid shifted");
        // This width cannot be represented as a global endpoint increment.
        const auto found=std::find_if(a.atomPolynomials().begin(),a.atomPolynomials().end(),[](const auto& p){return p.value(false,.3L*p.length)>0 && (!p.radialMeasure || p.radiusBegin>0 || p.radialTangent);});check(found!=a.atomPolynomials().end(),"no positive tiny interval");const auto& poly=*found;const long double offset=.3L*poly.length,width=1e-25L;check(offset+width==offset && poly.integral(false,offset,width)>0,"tiny interval lost positive width");
        norm(poly.integral(false,offset,width),poly.value(false,offset)*width*(rz?2*pi*(poly.radiusBegin+(poly.radialTangent?offset:0)):1.L),width,1e-12,"tiny shifted polynomial");
        check(copy.identity()==a.identity() && copy.atomPolynomials()[0].beginCoefficients==atomBefore,"failed evaluation changed original owner");
        auto equivalent=makeMesh(4,1);mesh.reset();copy.validateFor(*equivalent,geometry,id,clock);check(std::isfinite(copy.nativeVelocity(0).y),"owner lifetime after mesh release");
        const auto retained=[&]{const auto g=HydroMomentumGeometry::build(*equivalent,metricsFor(*equivalent,rz));const auto r=HydroAdvectorReconstruction::prepare(*equivalent,g);const HydroNativeValues v(g.layout(),ref.vel,std::vector<double>(g.layout().pressureCount()));return r.evaluate(*equivalent,v,id,clock);}();
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
void inspectCellJet(const HydroSpatialAdvector& source,int p,Vector2D point,
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
            const auto geometry=HydroMomentumGeometry::build(*mesh,metricsFor(*mesh,rz,selection));const auto& layout=geometry.layout();
            const auto reconstruction=HydroAdvectorReconstruction::prepare(*mesh,geometry);const auto partitions=physicalPartitions(reconstruction);
            for(int kind:{0,1,2}) {
                AdvectorReference reference(*mesh,layout,kind==0?0:1);auto native=reference.vel;
                if(kind==2)for(int u=0;u<layout.velocityCount();++u)
                    native[u]=(rz && geometry.metrics().faceMeasures()[u]==0)?0:.1*std::sin(1.7*u)+.07*std::cos(.4*u);
                const auto source=reconstruction.evaluate(*mesh,HydroNativeValues(layout,native,std::vector<double>(layout.pressureCount())),publication.state().identity(),publication.state().timeLayer());
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
        auto mesh=makeMesh(4,1);auto geometry=HydroMomentumGeometry::build(*mesh,metricsFor(*mesh,rz));
        Case shape;shape.pureFirstJet=true;shape.radial=rz;shape.planar=true;shape.planeOffset=2;const auto pub=initial(*mesh,shape);const auto id=pub.state().identity();const auto clock=pub.state().timeLayer();
        const auto retained=[&] {const auto prepared=HydroAdvectorReconstruction::prepare(*mesh,geometry);AdvectorReference ref(*mesh,geometry.layout(),1);
            return prepared.evaluate(*mesh,HydroNativeValues(geometry.layout(),ref.vel,std::vector<double>(geometry.layout().pressureCount())),id,clock);}();
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
        const auto freshGeometry=HydroMomentumGeometry::build(*equivalent,metricsFor(*equivalent,rz));const auto freshPub=initial(*equivalent,shape);
        const auto reconstruction=HydroAdvectorReconstruction::prepare(*equivalent,freshGeometry);AdvectorReference ref(*equivalent,freshGeometry.layout(),0);
        const auto fresh=reconstruction.evaluate(*equivalent,HydroNativeValues(freshGeometry.layout(),ref.vel,std::vector<double>(freshGeometry.layout().pressureCount())),freshPub.state().identity(),freshPub.state().timeLayer());
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
        auto mesh=makeMesh(4,1);const auto geometry=HydroMomentumGeometry::build(*mesh,metricsFor(*mesh,rz));
        const auto retained=[&]{const auto owner=HydroAdvectorReconstruction::prepare(*mesh,geometry);return owner.regularGeometry();}();
        const auto copy=retained;check(copy.sameIdentity(retained)&&copy.momentumGeometry().sameIdentity(geometry),"regular geometry lost original owner");
        const auto& l=geometry.layout();std::vector<double> zero(l.velocityCount()),pressure(l.pressureCount());
        const HydroNativeValues z(l,zero,pressure);const auto a=copy.apply(*mesh,z);
        for(const auto& c:a.cellDofs)for(const auto& v:c)for(double x:v)check(x==0,"regular zero field changed");
        for(double x:a.faceCorrections)check(x==0,"regular zero correction changed");
        check(copy.report().svdCount>0 && copy.report().weights>0 && copy.report().maxSupport<=24,"regular support/preparation inventory");
        // All physical P2 polynomials compatible with the RZ axis, summed.
        const auto fun=[&](Vector2D p,int axis){return rz?(axis==0?p.x+p.x*p.y:1+p.y+p.x*p.x+p.y*p.y):(axis==0?1+p.x+p.y+p.x*p.x+p.x*p.y+p.y*p.y:2-p.x+2*p.y);};
        const double q[5]={.046910077030668,.2307653449471585,.5,.7692346550528415,.953089922969332},w[5]={.1184634425280945,.2393143352496832,.2844444444444444,.2393143352496832,.1184634425280945};
        std::vector<double> u(l.velocityCount());
        for(int f=0;f<l.velocityCount();++f){const auto& face=copy.faces()[f];if(face.area==0)continue;for(int k=0;k<5;++k){const auto p=face.lower+(face.upper-face.lower)*q[k];u[f]+=w[k]*fun(p,face.normalAxis)*(rz?p.x/((face.lower.x+face.upper.x)/2):1);}}
        const HydroNativeValues v(l,u,pressure);const auto coefficients=copy.apply(*mesh,v);
        for(std::size_t p=0;p<copy.cells().size();++p){const auto& c=copy.cells()[p];for(int axis=0;axis<2;++axis){const auto& raw=coefficients.cellDofs[p][axis];
            check(std::abs(raw[24]-fun((c.lower+c.upper)*.5,axis))<1e-11,"regular P2 center reproduction");}}
        for(double delta:coefficients.faceCorrections)check(std::abs(delta)<1e-11,"regular P2 face correction");
        const auto again=copy.apply(*mesh,v);check(again.cellDofs==coefficients.cellDofs&&again.faceCorrections==coefficients.faceCorrections,"regular repeated application changed");
        std::vector<int> walls;for(int f=0;f<l.velocityCount();++f)if(copy.faces()[f].external&&copy.faces()[f].area>0)walls.push_back(f);
        const auto law=HydroAdvectorBoundary::impermeable(*mesh,geometry,walls);
        copy.apply(*mesh,z,&law);rejects([&]{copy.apply(*mesh,v,&law);},"nonzero prescribed wall");++guards;
        const auto other=HydroMomentumGeometry::build(*mesh,metricsFor(*mesh,rz));
        const auto wrongLaw=HydroAdvectorBoundary::impermeable(*mesh,other,walls);
        rejects([&]{copy.apply(*mesh,z,&wrongLaw);},"momentum geometry instance");
        rejects([&]{copy.validateFor(*mesh,other);},"momentum geometry instance");
        rejects([&]{copy.apply(*mesh,HydroNativeValues(other.layout(),zero,pressure));},"layout identity");guards+=3;
        const auto mode=std::fegetround();std::fesetround(FE_UPWARD);try{rejects([&]{copy.apply(*mesh,z);},"round-to-nearest");}catch(...){std::fesetround(mode);throw;}std::fesetround(mode);
        const auto control=_mm_getcsr();_mm_setcsr(control|(1U<<15));try{rejects([&]{copy.apply(*mesh,z);},"gradual underflow");}catch(...){_mm_setcsr(control);throw;}_mm_setcsr(control);guards+=2;
        if(rz){auto bad=zero;for(int f=0;f<l.velocityCount();++f)if(copy.faces()[f].area==0){bad[f]=1;break;}rejects([&]{copy.apply(*mesh,HydroNativeValues(l,bad,pressure));},"axis velocity");++guards;}
        mesh.reset();auto equivalent=makeMesh(4,1);const auto after=copy.apply(*equivalent,v);
        check(after.cellDofs==coefficients.cellDofs&&after.faceCorrections==coefficients.faceCorrections,"regular borrowed original mesh or caller values");
        DictionaryMesh::RefineBox halo;halo.level=1;halo.min={.125,.125};halo.max={.875,.875};equivalent->applyRefineBoxes({halo});
        DictionaryMesh::RefineBox box;box.level=2;box.min={.25,.25};box.max={.75,.75};equivalent->applyRefineBoxes({box});
        rejects([&]{copy.apply(*equivalent,v);},"NativeMeshIdentity");++guards;
        const auto fresh=HydroMomentumGeometry::build(*equivalent,metricsFor(*equivalent,rz));const auto rebuilt=HydroAdvectorReconstruction::prepare(*equivalent,fresh);
        check(!rebuilt.regularGeometry().sameIdentity(copy),"regular regrid reused geometry");
        rebuilt.regularGeometry().apply(*equivalent,HydroNativeValues(fresh.layout(),std::vector<double>(fresh.layout().velocityCount()),std::vector<double>(fresh.layout().pressureCount())));
    }
    check(guards==15,"regular guard inventory");std::cout<<"REGULAR_PREPARATION,2,"<<guards<<",PASS\n";
}

// Complete public normal polynomial: reference is a physical-coordinate
// antiderivative, independent of production shifted/Horner integration.
void completeNormalPolynomialContract()
{
    using Coefficients=decltype(HydroNormalPolynomial::beginCoefficients);
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
        HydroNormalPolynomial p;p.length=length;p.radialMeasure=metric>0;p.radialTangent=metric>=2;
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
        HydroNormalPolynomial p;p.length=1;p.radialMeasure=metric>0;p.radialTangent=metric>=2;
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
// Algebraic power-to-Bernstein conversion for actual compatible quarter patches.
// Independent of production evaluation and of its polynomial extraction.
std::array<long double,5> regularPowerControls(const std::array<long double,5>& a)
{
    constexpr long double ratio[5][5]={{1,0,0,0,0},{1,.25L,0,0,0},
        {1,.5L,1.L/6,0,0},{1,.75L,.5L,.25L,0},{1,1,1,1,1}};
    std::array<long double,5> out{};
    for(int k=0;k<5;++k)for(int i=0;i<=k;++i)out[k]+=ratio[k][i]*a[i];
    return out;
}
std::array<long double,6> regularJetValues(const HydroAdvectorCellJet& a)
{const auto& j=a.jet;return {j.value.x,j.value.y,j.derivative.xx,j.derivative.xy,j.derivative.yx,j.derivative.yy};}
struct RegularReference
{
    const LAMP::Hydrodynamics::Detail::HydroRegularAdvectorGeometry& geometry;
    const LAMP::Hydrodynamics::Detail::HydroRegularCoefficients& data;
    std::vector<std::array<std::array<long double,25>,2>> controls;
    RegularReference(const LAMP::Hydrodynamics::Detail::HydroRegularAdvectorGeometry& g,const LAMP::Hydrodynamics::Detail::HydroRegularCoefficients& d):geometry(g),data(d),controls(d.cellDofs.size())
    {
        for(std::size_t p=0;p<controls.size();++p)for(int c=0;c<2;++c){std::array<long double,25> temp{};
            for(int y=0;y<5;++y){std::array<long double,5> row{};for(int x=0;x<5;++x)row[x]=d.cellDofs[p][c][5*y+x];row=g.usesCompatibleTensor()?regularPowerControls(row):regularControls(row);for(int x=0;x<5;++x)temp[5*y+x]=row[x];}
            for(int x=0;x<5;++x){std::array<long double,5> col{};for(int y=0;y<5;++y)col[y]=temp[5*y+x];col=g.usesCompatibleTensor()?regularPowerControls(col):regularControls(col);for(int y=0;y<5;++y)controls[p][c][5*y+x]=col[y];}}
    }
    std::array<long double,6> sample(int p,Vector2D point) const
    {
        const auto& cell=geometry.cells()[p];auto lower=cell.lower,h=cell.upper-cell.lower;int patch=p;
        if(geometry.usesCompatibleTensor()){
            const int quarterX=std::clamp(int(std::floor(4*(point.x-lower.x)/h.x)),0,3),
                      quarterY=std::clamp(int(std::floor(4*(point.y-lower.y)/h.y)),0,3);
            patch=16*p+4*quarterY+quarterX;h=h*.25;
            lower=lower+Vector2D{quarterX*h.x,quarterY*h.y};
        }
        const long double sx=(point.x-lower.x)/h.x,sy=(point.y-lower.y)/h.y;std::array<long double,6> out{};
        for(int c=0;c<2;++c){std::array<long double,5> value{},dx{};
            for(int y=0;y<5;++y){std::array<long double,5> row{};std::array<long double,4> derivative{};
                for(int x=0;x<5;++x)row[x]=controls[patch][c][5*y+x];for(int x=0;x<4;++x)derivative[x]=4*(row[x+1]-row[x]);
                value[y]=regularBezier(row,sx);dx[y]=regularBezier(derivative,sx)/h.x;}
            std::array<long double,4> dy{};for(int y=0;y<4;++y)dy[y]=4*(value[y+1]-value[y]);
            out[c]=regularBezier(value,sy);out[2+2*c]=regularBezier(dx,sy);out[3+2*c]=regularBezier(dy,sy)/h.y;}
        if(geometry.usesCompatibleTensor())return out;
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
std::vector<double> regularCuts(const LAMP::Hydrodynamics::Detail::HydroRegularAdvectorGeometry& g,int p,int direction,double lo,double hi,int other=-1)
{
    std::vector<double> cuts{lo,hi};for(int cell:{p,other})if(cell>=0)for(double knot:g.cells()[cell].supportKnots[direction])if(lo<knot&&knot<hi)cuts.push_back(knot);
    std::sort(cuts.begin(),cuts.end());cuts.erase(std::unique(cuts.begin(),cuts.end()),cuts.end());return cuts;
}
void localSourceStencilContract()
{
    std::size_t samples=0,units=0,maxSupport=0;long double peak=0;
    for(int rz=0;rz<2;++rz)for(int amr=0;amr<2;++amr)for(int selection=0;selection<4;++selection){
        auto mesh=makeMesh(4,amr);const auto momentum=HydroMomentumGeometry::build(*mesh,metricsFor(*mesh,rz,selection));
        const auto owner=HydroAdvectorReconstruction::prepare(*mesh,momentum);const auto& g=owner.regularGeometry();const auto& l=momentum.layout();
        std::vector<int> walls;for(int u=0;u<l.velocityCount();++u)if(g.faces()[u].external&&g.faces()[u].area>0)walls.push_back(u);
        const auto boundary=HydroAdvectorBoundary::impermeable(*mesh,momentum,walls);
        for(int bound=0;bound<2;++bound){const auto* law=bound?&boundary:nullptr;std::vector<double> velocity(l.velocityCount());
            for(int u=0;u<l.velocityCount();++u)if(g.faces()[u].area>0&&(!bound||!std::binary_search(walls.begin(),walls.end(),u)))velocity[u]=std::sin(1.37*(u+1));
            const auto actual=g.apply(*mesh,HydroNativeValues(l,velocity,std::vector<double>(l.pressureCount())),law);
            for(int p=0;p<l.pressureCount();++p){const auto stencil=g.cellStencil(*mesh,p,law);const auto& cell=g.cells()[p];
                for(Vector2D s:std::vector<Vector2D>{{0,0},{.31,.63},{1,1}}){const Vector2D x=cell.lower+Vector2D{s.x*(cell.upper.x-cell.lower.x),s.y*(cell.upper.y-cell.lower.y)};
                    const auto weights=g.cellJetWeights(stencil,x);maxSupport=std::max(maxSupport,weights.size());std::array<long double,6> recovered{};
                    for(const auto& w:weights){check(w.velocityDof>=0&&w.velocityDof<l.velocityCount(),"local stencil native index");for(int k=0;k<6;++k)recovered[k]+=w.jet[k]*velocity[w.velocityDof];}
                    const auto exact=regularJetValues(g.cellJet(actual,p,x));for(int k=0;k<6;++k)peak=std::max(peak,norm(recovered[k],exact[k],1+std::abs(exact[k]),1e-11,"local stencil actual source jet"));++samples;
                    if(p==0&&s.x==.31){for(int u=0;u<l.velocityCount();++u){if(g.faces()[u].area==0||(bound&&std::binary_search(walls.begin(),walls.end(),u)))continue;
                            std::vector<double> unit(l.velocityCount());unit[u]=1;const auto basis=g.apply(*mesh,HydroNativeValues(l,unit,std::vector<double>(l.pressureCount())),law);
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
        auto mesh=makeMesh(8,amr);const auto momentum=HydroMomentumGeometry::build(*mesh,metricsFor(*mesh,rz,selection));
        const auto owner=HydroAdvectorReconstruction::prepare(*mesh,momentum);const auto& g=owner.regularGeometry();const auto& l=momentum.layout();
        check(g.normalPorts().size()==g.report().normalPorts,"normal port preparation inventory");
        for(int c=0;c<2;++c)for(int degree=0;degree<=3;++degree){
            std::vector<double> velocity(l.velocityCount());
            for(int u=0;u<l.velocityCount();++u){const auto& face=g.faces()[u];if(face.normalAxis!=c||face.area==0)continue;
                for(int q=0;q<3;++q){const auto x=face.lower+(face.upper-face.lower)*static_cast<double>(quadrature[q]);
                    velocity[u]+=static_cast<double>(weights[q]*std::pow(component(x,1-c),degree)*(rz?x.x/((face.lower.x+face.upper.x)/2):1));}}
            const auto data=g.apply(*mesh,HydroNativeValues(l,velocity,std::vector<double>(l.pressureCount())));
            // Uniform tensor reconstruction has no line-port fitting. Test
            // the same independent line polynomials at actual shared vertices.
            if(g.usesCompatibleTensor())for(int p=0;p<l.pressureCount();++p){
                const auto& cell=g.cells()[p];
                for(Vector2D x:{cell.lower,cell.upper,Vector2D{cell.lower.x,cell.upper.y},Vector2D{cell.upper.x,cell.lower.y}}){
                    const auto j=regularJetValues(g.cellJet(data,p,x));const long double t=component(x,1-c);
                    peak=std::max(peak,norm(j[c],std::pow(t,degree),1,1e-11,"compatible native line cubic value"));
                    peak=std::max(peak,norm(j[2+2*c+1-c],degree?degree*std::pow(t,degree-1):0,1,1e-11,"compatible native line cubic derivative"));++ports;
                }
            }
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
        const auto data=g.apply(*mesh,HydroNativeValues(l,velocity,std::vector<double>(l.pressureCount())));
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
        auto mesh=makeMesh(n,amr);const auto momentum=HydroMomentumGeometry::build(*mesh,metricsFor(*mesh,rz,selected));
        const auto unrestricted=HydroAdvectorReconstruction::prepare(*mesh,momentum);const auto& initialGeometry=unrestricted.regularGeometry();const auto& l=momentum.layout();
        std::vector<int> walls;for(int u=0;u<l.velocityCount();++u)if(initialGeometry.faces()[u].external&&initialGeometry.faces()[u].area>0)walls.push_back(u);
        check(!walls.empty(),"missing actual regular wall fixture");
        const auto polynomial=[&](Vector2D x){return rz?std::array<long double,6>{x.x*(1+x.y),1+x.y+x.x*x.x+x.y*x.y,1+x.y,x.x,2*x.x,1+2*x.y}
            :std::array<long double,6>{1+x.x+x.y+x.x*x.x+x.x*x.y+x.y*x.y,2-x.x+2*x.y,1+2*x.x+x.y,1+x.x+2*x.y,-1,2};};
        for(int bound=0;bound<3;++bound){auto chosen=bound==2?std::vector<int>{walls.front()}:walls;const auto law=HydroAdvectorBoundary::impermeable(*mesh,momentum,chosen);
            const auto owner=bound?HydroAdvectorReconstruction::prepare(*mesh,momentum,law):unrestricted;
            const auto& geometry=owner.regularGeometry();
            const bool identifiable=std::all_of(geometry.rows().begin(),geometry.rows().end(),[](const auto& r){return r.zeroNormal||r.p2Rank==r.p2ColumnCount;});
            if(bound==0&&identifiable)++p2Layouts;
            const HydroNativeValues zeros(l,std::vector<double>(l.velocityCount()),std::vector<double>(l.pressureCount()));
            if(bound==2&&initialGeometry.usesCompatibleTensor()){
                check(!geometry.usesCompatibleTensor()&&!geometry.rows().empty(),"partial wall was not prepared with the actual constrained source");
                rejects([&]{initialGeometry.apply(*mesh,zeros,&law);},"partial impermeable boundary requires boundary-aware source preparation");
                rejects([&]{initialGeometry.cellStencil(*mesh,0,&law);},"partial impermeable boundary requires boundary-aware source preparation");
            }else check(geometry.usesCompatibleTensor()==initialGeometry.usesCompatibleTensor(),"complete boundary changed prepared method");
            const auto zero=geometry.apply(*mesh,zeros,bound?&law:nullptr);
            for(int p=0;p<l.pressureCount();++p){const auto& cell=geometry.cells()[p];const auto z=regularJetValues(geometry.cellJet(zero,p,(cell.lower+cell.upper)*.5));for(auto v:z)check(v==0,"regular zero jet not exact");}
            for(int family=0;family<2;++family){std::vector<double> velocity(l.velocityCount());
                for(int u=0;u<l.velocityCount();++u){const auto& f=geometry.faces()[u];if(f.area==0)continue;
                    for(int k=0;k<3;++k){const auto x=f.lower+(f.upper-f.lower)*static_cast<double>((1+nodes[k])/2);
                        const auto value=family==0?polynomial(x)[f.normalAxis]:std::sin(7*x.x+3*x.y)+.07*std::cos(2*x.x-5*x.y)+.03*std::sin(u+1.);
                        velocity[u]+=static_cast<double>(weights[k]/2*value*(rz?x.x/((f.lower.x+f.upper.x)/2):1));}}
                if(bound)for(int u:chosen)velocity[u]=0;
                const HydroNativeValues values(l,velocity,std::vector<double>(l.pressureCount()));const auto updateStart=Clock::now();const auto data=geometry.apply(*mesh,values,bound?&law:nullptr);
                updateSeconds+=std::chrono::duration<double>(Clock::now()-updateStart).count();const RegularReference reference(geometry,data);
                if(bound==2&&initialGeometry.usesCompatibleTensor())for(int p=0;p<l.pressureCount();++p){
                    const auto& cell=geometry.cells()[p];const auto x=(cell.lower+cell.upper)*.5;
                    const auto stencil=geometry.cellStencil(*mesh,p,&law);const auto weights=geometry.cellJetWeights(stencil,x);
                    std::array<long double,6> symbolic{};for(const auto& row:weights)for(int k=0;k<6;++k)symbolic[k]+=row.jet[k]*velocity[row.velocityDof];
                    const auto actual=regularJetValues(geometry.cellJet(data,p,x));
                    for(int k=0;k<6;++k)norm(actual[k],symbolic[k],1+std::abs(symbolic[k]),1e-11,"partial wall actual/symbolic source");
                }
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
                std::vector<HydroAdvectorCellJet> measured(l.pressureCount());const auto kernelStart=Clock::now();
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
    std::cout<<"REGULAR_EVALUATION_KERNEL,"<<kernelQueries<<','<<kernelSeconds<<','<<maxScalarStorage<<','<<sizeof(LAMP::Hydrodynamics::Detail::HydroRegularCoefficients)<<'\n';
    std::cout<<"REGULAR_EVALUATION_NATIVE,"<<nativeAtoms<<','<<quartic<<",PASS\n";
}
void regularEvaluationGuards()
{
    int count=0;for(bool rz:{false,true}){
        auto mesh=makeMesh(4,1);const auto g=HydroMomentumGeometry::build(*mesh,metricsFor(*mesh,rz));const auto original=HydroAdvectorReconstruction::prepare(*mesh,g);
        const auto prepared=original.regularGeometry();const auto& l=g.layout();const HydroNativeValues zeros(l,std::vector<double>(l.velocityCount()),std::vector<double>(l.pressureCount()));const auto data=prepared.apply(*mesh,zeros);
        const auto& cell=prepared.cells()[0];const auto point=(cell.lower+cell.upper)*.5;const auto checkReject=[&](auto call,const char* text){rejects(call,text);++count;};
        LAMP::Hydrodynamics::Detail::HydroRegularCoefficients unbound;
        checkReject([&]{prepared.cellJet(unbound,0,point);},"coefficient geometry identity");
        const auto other=HydroAdvectorReconstruction::prepare(*mesh,g);checkReject([&]{other.regularGeometry().cellJet(data,0,point);},"coefficient geometry identity");
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
        auto temporary=HydroAdvectorReconstruction::prepare(*mesh,g).regularGeometry();const auto retained=temporary.apply(*mesh,zeros);mesh.reset();
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
        auto mesh=makeMesh(n,amr);const auto momentum=HydroMomentumGeometry::build(*mesh,metricsFor(*mesh,rz,selection));
        const auto owner=HydroAdvectorReconstruction::prepare(*mesh,momentum);const auto& g=owner.regularGeometry();const auto& l=momentum.layout();
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
            const auto data=g.apply(*mesh,HydroNativeValues(l,velocity,std::vector<double>(l.pressureCount())));
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
        auto mesh=makeMesh(n,amr);const auto momentum=HydroMomentumGeometry::build(*mesh,metricsFor(*mesh,rz,selection));
        const auto owner=HydroAdvectorReconstruction::prepare(*mesh,momentum);const auto& g=owner.regularGeometry();const auto& l=momentum.layout();
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
            const auto law=HydroAdvectorBoundary::impermeable(*mesh,momentum,walls);
            for(int bound:{0,1})for(double scale:(family==1&&n==4&&selection==3?std::vector<double>{0,1,1e-150,1e150}:std::vector<double>{1})){
                std::vector<double> velocity(l.velocityCount());
                for(int u=0;u<l.velocityCount();++u){const auto& f=g.faces()[u];if(f.area>0)velocity[u]=scale*static_cast<double>(exact((f.lower+f.upper)*.5)[f.normalAxis]);}
                const auto data=g.apply(*mesh,HydroNativeValues(l,velocity,std::vector<double>(l.pressureCount())),bound?&law:nullptr);
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

// H14a3r2: analytic physical polynomials supply independent face means
// and jets. P3 observability NEVER narrows original P2 eligibility.
void regularCubicSourceContract()
{
    constexpr int powers[10][2]={{0,0},{1,0},{0,1},{2,0},{1,1},{0,2},{3,0},{2,1},{1,2},{0,3}};
    const long double node=std::sqrt(3.L/5),nodes[3]={-node,0,node},weights[3]={5.L/9,8.L/9,5.L/9};
    int fields=0,deficientP3=0,independentP2=0,wallFields=0;std::size_t jets=0,integrals=0;
    long double peak=0,peakMean=0,peakP2Deficit=0;
    for(int n:{3,4}){
        auto mesh=makeMesh(n);const auto momentum=HydroMomentumGeometry::build(*mesh,metricsFor(*mesh,false));
        const auto owner=HydroAdvectorReconstruction::prepare(*mesh,momentum);const auto& g=owner.regularGeometry();const auto& l=momentum.layout();
        check(g.usesCompatibleTensor()&&g.rows().empty()&&g.normalPorts().empty()&&g.report().svdCount==0,
            "uniform source retained obsolete moment fits or normal ports");
        // Independent information certificate of the actual native face means.
        // Three transverse intervals identify P2; four identify full P3.
        for(int c=0;c<2;++c){
            std::vector<int> observations;
            for(int u=0;u<l.velocityCount();++u)if(g.faces()[u].normalAxis==c)observations.push_back(u);
            Eigen::MatrixXd observed(observations.size(),10);observed.setZero();
            for(std::size_t i=0;i<observations.size();++i){const auto& face=g.faces()[observations[i]];
                for(int k=0;k<10;++k)for(int q=0;q<3;++q){
                    const auto x=face.lower+(face.upper-face.lower)*static_cast<double>((1+nodes[q])/2);
                    observed(i,k)+=static_cast<double>(weights[q]/2*std::pow(x.x,powers[k][0])*std::pow(x.y,powers[k][1]));
                }
            }
            Eigen::JacobiSVD<Eigen::MatrixXd> p2(observed.leftCols(6)),p3(observed);
            const double threshold=128*std::numeric_limits<double>::epsilon()*std::max<std::size_t>(observations.size(),10);
            p2.setThreshold(threshold);p3.setThreshold(threshold);
            check(p2.rank()==6,"native observations lost independent P2 information");
            if(n==3){++independentP2;deficientP3+=p3.rank()<10;}
            else check(p3.rank()==10,"P3 precision fixture lacks independent native information");
        }
        std::vector<int> walls;for(int u=0;u<l.velocityCount();++u)if(g.faces()[u].external&&g.faces()[u].area>0)walls.push_back(u);
        const auto law=HydroAdvectorBoundary::impermeable(*mesh,momentum,walls);
        Case shape;shape.pureFirstJet=true;shape.planar=true;shape.planeOffset=2;const auto publication=initial(*mesh,shape);
        const int degreeCount=n==3?6:10;
        for(int family=0;family<2*degreeCount+3;++family){const bool wall=family==2*degreeCount+2;
            const auto exact=[&](Vector2D x){std::array<long double,6> out{};
                if(wall){
                    // n3 deliberately lacks full P3 information: retain an
                    // independent old P2 wall field there. n4 checks true
                    // mixed-cubic wall compatibility of the actual source.
                    out=n==3?std::array<long double,6>{x.x*(1-x.x),x.y*(1-x.y),1-2*x.x,0,0,1-2*x.y}
                        :std::array<long double,6>{x.x*(1-x.x)*(1+x.y),x.y*(1-x.y)*(1+x.x),(1-2*x.x)*(1+x.y),x.x*(1-x.x),x.y*(1-x.y),(1-2*x.y)*(1+x.x)};return out;}
                for(int c=0;c<2;++c)for(int k=0;k<degreeCount;++k){
                    const long double factor=family<2*degreeCount?(family==c*degreeCount+k?1.L:0.L):family==2*degreeCount?((k+c)%2?-1.L:1.L)*(k+1)/11:((k+c)%2?1.L:-1.L)*(degreeCount-k)/13;
                    const int a=powers[k][0],b=powers[k][1];out[c]+=factor*std::pow(static_cast<long double>(x.x),a)*std::pow(static_cast<long double>(x.y),b);
                    if(a)out[2+2*c]+=factor*a*std::pow(static_cast<long double>(x.x),a-1)*std::pow(static_cast<long double>(x.y),b);
                    if(b)out[3+2*c]+=factor*b*std::pow(static_cast<long double>(x.x),a)*std::pow(static_cast<long double>(x.y),b-1);
                }return out;
            };
            std::vector<double> velocity(l.velocityCount());
            for(int u=0;u<l.velocityCount();++u){const auto& f=g.faces()[u];long double value=0;
                for(int q=0;q<3;++q){const auto x=f.lower+(f.upper-f.lower)*static_cast<double>((1+nodes[q])/2);value+=weights[q]/2*exact(x)[f.normalAxis];}
                velocity[u]=static_cast<double>(value);
            }
            if(wall)for(int u:walls)check(velocity[u]==0,"wall-compatible cubic input has nonzero normal mean");
            const HydroNativeValues values(l,velocity,std::vector<double>(l.pressureCount()));
            const auto source=wall?owner.evaluate(*mesh,values,publication.state().identity(),publication.state().timeLayer(),law)
                :owner.evaluate(*mesh,values,publication.state().identity(),publication.state().timeLayer());
            for(int p=0;p<l.pressureCount();++p){const auto& cell=g.cells()[p];const auto h=cell.upper-cell.lower;
                for(Vector2D f:std::vector<Vector2D>{{0,0},{.21,.31},{.63,.77},{1,1}}){const auto x=cell.lower+Vector2D{f.x*h.x,f.y*h.y};
                    const auto a=regularJetValues(source.cellVelocityJet(p,x)),b=exact(x);
                    for(int k=0;k<6;++k)peak=std::max(peak,norm(a[k],b[k],1+std::abs(b[k]),1e-12,"same-source P2/P3 analytic physical jet n="+std::to_string(n)+" family="+std::to_string(family)+" cell="+std::to_string(p)+" entry="+std::to_string(k)+" x="+std::to_string(x.x)+" y="+std::to_string(x.y)));++jets;
                }
            }
            const auto coefficients=g.apply(*mesh,HydroNativeValues(l,velocity,std::vector<double>(l.pressureCount())),wall?&law:nullptr);
            for(int u=0;u<l.velocityCount();++u){const auto& f=g.faces()[u];const int p=f.negativePressure>=0?f.negativePressure:f.positivePressure,t=1-f.normalAxis;
                const auto cuts=regularCuts(g,p,t,component(f.lower,t),component(f.upper,t));long double sum=0;
                for(std::size_t i=1;i<cuts.size();++i){auto a=f.lower,b=f.upper;(t==0?a.x:a.y)=cuts[i-1];(t==0?b.x:b.y)=cuts[i];
                    const auto poly=g.normalPolynomial(coefficients,p,f.normalAxis,a,b);long double expected=0;
                    for(int q=0;q<3;++q){const auto x=a+(b-a)*static_cast<double>((1+nodes[q])/2);expected+=weights[q]*poly.length/2*exact(x)[f.normalAxis];}
                    peakMean=std::max({peakMean,norm(poly.integral(false,0,poly.length),expected,1+std::abs(expected),1e-12,"P2/P3 begin physical normal integral"),norm(poly.integral(true,0,poly.length),expected,1+std::abs(expected),1e-12,"P2/P3 end physical normal integral")});
                    sum+=poly.integral(false,0,poly.length);++integrals;
                }
                peakMean=std::max(peakMean,norm(sum,f.area*velocity[u],1+std::abs(f.area*velocity[u]),1e-12,"P2/P3 actual native face mean"));
            }
            ++fields;wallFields+=wall;
        }
    }
    check(deficientP3>0&&independentP2>0&&wallFields==2&&fields==38,"P3/P2-alias or independent polynomial coverage changed");
    std::cout<<"CUBIC_SOURCE_SUMMARY,"<<fields<<','<<jets<<','<<integrals<<','<<independentP2<<','<<deficientP3<<','<<wallFields<<','<<peak<<','<<peakMean<<','<<peakP2Deficit<<",PASS\n";
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
        for(int selection=0;selection<4;++selection){const auto momentum=HydroMomentumGeometry::build(*mesh,metricsFor(*mesh,rz,selection));const auto unrestricted=HydroAdvectorReconstruction::prepare(*mesh,momentum);const auto& initialGeometry=unrestricted.regularGeometry();const auto& l=momentum.layout();
            std::vector<int> walls;for(int u=0;u<l.velocityCount();++u)if(initialGeometry.faces()[u].external&&initialGeometry.faces()[u].area>0)walls.push_back(u);check(!walls.empty(),"missing matching boundary fixture");
            for(int boundary=0;boundary<3;++boundary){auto chosen=boundary==2?std::vector<int>{walls.front()}:walls;const auto law=HydroAdvectorBoundary::impermeable(*mesh,momentum,chosen);
                const auto owner=boundary?HydroAdvectorReconstruction::prepare(*mesh,momentum,law):unrestricted;const auto& g=owner.regularGeometry();
                const auto oldCuts=physicalPartitions(owner);std::vector<PhysicalCuts> regularPartitions;for(const auto& c:g.cells())regularPartitions.push_back(c.supportKnots);
                for(int family=0;family<3;++family){AdvectorReference input(*mesh,l,family==0?0:1);auto velocity=input.vel;
                    if(family==2)for(int u=0;u<l.velocityCount();++u)velocity[u]=g.faces()[u].area==0?0:.1*std::sin(1.7*u)+.07*std::cos(.4*u);
                    if(boundary)for(int u:chosen)velocity[u]=0;
                    const HydroNativeValues values(l,velocity,std::vector<double>(l.pressureCount()));const auto mean=momentum.metrics().apply(*mesh,values);
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
        const auto geometry=HydroMomentumGeometry::build(*mesh,metricsFor(*mesh,rz));
        const auto source=HydroAdvectorReconstruction::prepare(*mesh,geometry);
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
        // Query the actual source API: compatible stencils are opaque handles.
        // Five interior points per direction observe every Q4 piece polynomial.
        constexpr double supportNodes[5]={.1,.3,.5,.7,.9};
        for(const auto& piece:source.pieces())if(piece.pressure==cut){const auto h=piece.upper-piece.lower;
            for(double sx:supportNodes)for(double sy:supportNodes){
                const auto point=piece.lower+Vector2D{sx*h.x,sy*h.y};
                for(const auto& weight:regular.cellJetWeights(stencil,point))
                    if(regular.faces()[weight.velocityDof].normalAxis==1&&
                        (weight.jet[1]!=0||weight.jet[4]!=0||weight.jet[5]!=0))support.insert(weight.velocityDof);
            }}
        check(support.size()>=2,"actual source support has fewer than two observations");
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

long double independentModalAtomRate(const Case& shape,const HydroSpatialAdvector& state,const HydroAdvectorAtom& atom)
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
        const auto geometry=HydroMomentumGeometry::build(*mesh,metricsFor(*mesh,rz));
        const auto phase=HydroPhaseGeometry::build(*mesh,geometry.metrics(),initial(*mesh,shape));
        const auto source=HydroAdvectorReconstruction::prepare(*mesh,geometry).withTangentialModes(*mesh,phase);
        const HydroNativeValues zero(geometry.layout(),std::vector<double>(geometry.layout().velocityCount()),
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
        const auto pressure=HydroPressureLoad::build(geometry.metrics().divergence());
        const auto load=HydroInterfaceNormalLoad::capillary(*mesh,pressure,phase,1,source);
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
        const auto rate=HydroPhaseRate::prepareSpatial(*mesh,phase,state);
        const auto flux=HydroMaterialFlux::prepare(*mesh,rate,1000,1);long double peakRate=0;
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
        const auto geometry=HydroMomentumGeometry::build(*mesh,metricsFor(*mesh,rz));const auto& layout=geometry.layout();
        const auto phase=HydroPhaseGeometry::build(*mesh,geometry.metrics(),initial(*mesh,shape));
        const auto source=HydroAdvectorReconstruction::prepare(*mesh,geometry).withTangentialModes(*mesh,phase);
        const auto hodge=HydroKineticHodge::build(*mesh,HydroMomentumMass::build(*mesh,geometry,phase,1000,1),source);
        const int nu=layout.velocityCount(),nm=static_cast<int>(source.tangentialModes().size()),nv=nu+nm;
        check(hodge.augmentedMatrix().rows()==nv&&hodge.report().tangentialModes==static_cast<std::size_t>(nm),"missing full modal Hodge blocks");
        std::vector<double> modal(nm),zero(nu),normal(nu);for(int m=0;m<nm;++m)modal[m]=.1*(m+1);
        const auto psi=[&](Vector2D p){const double radial=rz?p.x*p.x*(1-p.x*p.x)*(1-p.x*p.x):p.x*p.x*(1-p.x)*(1-p.x);
            return radial*p.y*p.y*(1-p.y)*(1-p.y);};
        for(int u=0;u<nu;++u){const auto& face=source.regularGeometry().faces()[u];if(face.area==0)continue;
            normal[u]=static_cast<double>((face.normalAxis==0?1:-1)*(rz?2*pi:1)*(psi(face.upper)-psi(face.lower))/face.area);}
        const HydroNativeValues native(layout,normal,std::vector<double>(layout.pressureCount()));
        const HydroNativeValues rest(layout,zero,std::vector<double>(layout.pressureCount()));
        const auto onlyMode=source.evaluate(*mesh,rest,modal,phase.geometryIdentity(),phase.timeLayer());
        const auto initial=source.evaluate(*mesh,native,modal,phase.geometryIdentity(),phase.timeLayer());
        const auto uniqueFlux=HydroMaterialFlux::prepare(*mesh,HydroPhaseRate::prepareSpatial(*mesh,phase,onlyMode),1000,1);
        const auto convection=HydroKepConvection::modalBulkTransportMatrix(*mesh,hodge,uniqueFlux);
        const HydroKineticHodge::Matrix symmetric=convection+HydroKineticHodge::Matrix(convection.transpose());
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
        double cross=0;for(int u=0;u<nu;++u)for(HydroKineticHodge::Matrix::InnerIterator it(hodge.augmentedMatrix(),u);it;++it)
            if(it.col()>=nu)cross=std::max(cross,std::abs(it.value()));
        check(cross>1e-10,"modal Hodge cross block disappeared");
        const auto pressure=HydroPressureLoad::build(geometry.metrics().divergence());
        Eigen::VectorXd before(nv);for(int u=0;u<nu;++u)before[u]=normal[u];for(int m=0;m<nm;++m)before[nu+m]=modal[m];
        const double dt=.01;long double peakBalance=0;
        for(int scenario:{0,1,2,3}){
            const auto k=HydroPhaseViscosity::modalBulkMatrix(*mesh,hodge,scenario==2?.1:0,scenario==2?.02:0);
            if(scenario!=2)check(k.nonZeros()==0,"exact mu=0 modal viscosity is not zero");
            HydroSaddleVelocityBlock::Matrix a=hodge.augmentedMatrix()/dt+k*.5;
            Eigen::VectorXd rhs=hodge.augmentedMatrix()*before/dt-k*before*.5;
            if(scenario==3)rhs-=convection*before*.5;
            Eigen::VectorXd forcing=Eigen::VectorXd::Zero(nv);
            if(scenario==1){for(int i=0;i<nv;++i)forcing[i]=force[i];rhs+=forcing;}
            std::vector<double> nativeRhs(nu),modeRhs(nm);
            for(int u=0;u<nu;++u)nativeRhs[u]=rhs[u];for(int m=0;m<nm;++m)modeRhs[m]=rhs[nu+m];
            auto block=HydroSaddleVelocityBlock::build(layout,a,nm);
            if(scenario==3)block=block.withSkew(HydroSaddleVelocityBlock::Matrix(convection*.5));
            const auto saddle=HydroSinglePhaseSaddle::buildClosed(pressure,block);
            const auto solved=saddle.solve(*mesh,HydroSaddleRhs(layout,nativeRhs,zero,modeRhs));
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
        const auto geometry=HydroMomentumGeometry::build(*mesh,metricsFor(*mesh,rz));const auto& layout=geometry.layout();
        const auto base=HydroAdvectorReconstruction::prepare(*mesh,geometry);
        const auto at=[&](HermiteGalsTimeLayer clock){const auto f=shape.profile(clock.time());
            const auto seed=stateFor(*mesh,f,rz,clock,[](int,Vector2D,Vector2D){return 0.;},shape.symmetry(*mesh),HermiteGalsTransport::FirstJet);
            HermiteGalsCandidate candidate(*mesh,clock,seed.primary());const auto mixed=closeHermiteGalsMixedAmr(*mesh,candidate,shape.symmetry(*mesh));
            HermiteGalsPublicationConfig config;config.volume.symmetryPlanes=shape.symmetry(*mesh);const auto id=candidate.identity();
            return HydroPhaseGeometry::build(*mesh,geometry.metrics(),publishHermiteGalsFirstJetGeometry(*mesh,std::move(candidate),mixed,{id,centers(*mesh,f)},config));};
        const auto g0=at({}),gm=at({1,0,dt,Position::Midpoint}),g1=at({1,0,dt,Position::End});
        const auto r0=base.withTangentialModes(*mesh,g0),rm=base.withTangentialModes(*mesh,gm),r1=base.withTangentialModes(*mesh,g1);
        const auto m0=HydroMomentumMass::build(*mesh,geometry,g0,1000,1),mm=HydroMomentumMass::build(*mesh,geometry,gm,1000,1),m1=HydroMomentumMass::build(*mesh,geometry,g1,1000,1);
        const auto h0=HydroKineticHodge::build(*mesh,m0,r0),hm=HydroKineticHodge::build(*mesh,mm,rm),h1=HydroKineticHodge::build(*mesh,m1,r1);
        const int nu=layout.velocityCount(),nm=static_cast<int>(rm.tangentialModes().size());
        std::vector<double> u0(nu),um(nu),u1(nu),q0(nm),qm(nm),q1(nm);
        for(int u=0;u<nu;++u){const auto& f=base.regularGeometry().faces()[u];if(f.normalAxis!=1)continue;
            const long double a=f.lower.x,b=f.upper.x;
            u0[u]=static_cast<double>(rz?(b*b+a*a)/2:(b*b+a*b+a*a)/3);u1[u]=1.05*u0[u];um[u]=(u0[u]+u1[u])/2;}
        for(int m=0;m<nm;++m){q0[m]=.2*(m+1);q1[m]=1.02*q0[m];qm[m]=(q0[m]+q1[m])/2;}
        const auto make=[&](const HydroAdvectorReconstruction& r,const HydroPhaseGeometry& g,const std::vector<double>& u,const std::vector<double>& q){
            return r.evaluate(*mesh,HydroNativeValues(layout,u,std::vector<double>(layout.pressureCount())),q,g.geometryIdentity(),g.timeLayer());};
        const auto before=make(r0,g0,u0,q0),mid=make(rm,gm,um,qm),after=make(r1,g1,u1,q1);
        const auto transport=HydroMaterialTransport::prepare(*mesh,g0,HydroPhaseRate::prepareSpatial(*mesh,gm,mid),g1,1000,1);
        const auto material=HydroDualMaterial::prepare(*mesh,transport,m0,mm,m1);
        const auto report=material.modalEnergyChange(*mesh,h0,hm,h1,before,after);
        const auto& change=report.kinetic;
        check(std::abs(change.identityResidual)<=2e-12L*(1+change.beginEnergy+change.endEnergy),"moving full modal metric polarization identity");
        check(std::abs(change.metricVariation)>1e-6&&material.report().components[1].l1IntegratedDefect>0,
            "moving modal test hid metric variation or material R");
        Eigen::VectorXd x(nu+nm);for(int u=0;u<nu;++u)x[u]=um[u];for(int m=0;m<nm;++m)x[nu+m]=qm[m];
        const HydroKineticHodge::Matrix difference=h1.augmentedMatrix()-h0.augmentedMatrix();
        long double cross=0;
        for(int u=0;u<nu;++u)for(HydroKineticHodge::Matrix::InnerIterator it(difference,u);it;++it)
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
        const auto geometry=HydroMomentumGeometry::build(*mesh,metricsFor(*mesh,rz));
        const auto phase=HydroPhaseGeometry::build(*mesh,geometry.metrics(),initial(*mesh,shape));
        const auto mass=HydroMomentumMass::build(*mesh,geometry,phase,1000,1);
        const auto source=HydroAdvectorReconstruction::prepare(*mesh,geometry);
        const auto hodge=HydroKineticHodge::build(*mesh,mass,source);
        hodge.validateFor(*mesh,phase.geometryIdentity(),phase.timeLayer());
        check(hodge.markedDuals().size()==static_cast<std::size_t>(geometry.layout().velocityCount())
            && hodge.report().markedDuals>0 && hodge.report().points>0 && hodge.report().nonZeros>0,
            "phase-cut kinetic Hodge missing marked quadrature");
        check(hodge.report().symmetryDefect<=1e-12 && hodge.report().maxRelativePhaseVolumeDifference<=1e-9,
            "phase-cut kinetic Hodge symmetry/phase measure");
        std::vector<double> velocity(geometry.layout().velocityCount(),0);
        for(int u=0;u<geometry.layout().velocityCount();++u)
            if(geometry.layout().normalAxis(u)==HydroDofLayout::NormalAxis::Second)velocity[u]=1;
        const HydroNativeValues values(geometry.layout(),velocity,std::vector<double>(geometry.layout().pressureCount()));
        const auto action=hodge.apply(*mesh,values,phase.geometryIdentity(),phase.timeLayer());
        long double matrixPair=0;
        for(int u=0;u<geometry.layout().velocityCount();++u){
            long double row=0;
            for(HydroKineticHodge::Matrix::InnerIterator it(hodge.matrix(),u);it;++it)
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
// Bounded structural-template contract. Mathematical references use actual
// native edges and physical polynomials, never the template key or its routing.
void regularTemplateBankContract()
{
    using Regular=LAMP::Hydrodynamics::Detail::HydroRegularAdvectorGeometry;
    std::size_t p2=0,ports=0,means=0,axes=0,walls=0,guards=0;
    std::array<std::size_t,3> kinds{};long double peak=0;
    const auto meshFor=[](double radial,double axial,bool expanded){
        DictionaryMesh d;d.Nx=d.Ny=8;d.minCoord={radial,axial};d.maxCoord={radial+1,axial+1};
        d.maxCellNum=16384;d.maxInitialIterations=0;auto mesh=std::make_unique<QuadtreeMesh>(d);
        DictionaryMesh::RefineBox box;box.level=1;const double a=expanded?.125:.25,b=expanded?.875:.75;
        box.min={radial+a,axial+a};box.max={radial+b,axial+b};mesh->applyRefineBoxes({box});
        double smallest=1,largest=0;for(int p=0;p<mesh->cellIndexNum();++p){const auto h=mesh->cells(mesh->cellIndexToLabel[p]).size().x;smallest=std::min(smallest,h);largest=std::max(largest,h);}
        check(largest==2*smallest,"template fixture is not actual mixed-level 2:1");return mesh;};
    const auto wallIds=[](const Regular& g){std::vector<int> result;
        for(std::size_t u=0;u<g.faces().size();++u)if(g.faces()[u].external&&g.faces()[u].area>0)result.push_back(static_cast<int>(u));return result;};
    const auto physical=[](Vector2D x){return std::array<long double,6>{
        x.x*(1.L+x.y),1.L+x.y+x.x*x.x+x.y*x.y,1.L+x.y,x.x,2.L*x.x,1.L+2.L*x.y};};
    const auto compare=[&](const QuadtreeMesh& mesh,const Regular& warm,const Regular& cold,
        const HydroAdvectorBoundary* warmWall,const HydroAdvectorBoundary* coldWall,int family){
        const auto& l=warm.momentumGeometry().layout();const bool bound=warmWall!=nullptr;
        check(l.velocityCount()==cold.momentumGeometry().layout().velocityCount(),"template target changed cold native inventory");
        const auto ends=[&](int u){const auto f=l.domain().face(mesh,l.nativeFace(u));
            const auto& edge=f.horizontal?mesh.hEdges(f.edgeLabel):mesh.vEdges(f.edgeLabel);
            std::array<Vector2D,2> e{mesh.verts(edge.edgeToVert[0]).coord(),mesh.verts(edge.edgeToVert[1]).coord()};
            if(component(e[0],1-axis(l,u))>component(e[1],1-axis(l,u)))std::swap(e[0],e[1]);return e;};
        std::vector<double> velocity(l.velocityCount());
        for(int u=0;u<l.velocityCount();++u){const auto e=ends(u);const auto& f=warm.faces()[u];
            check(l.nativeFace(u)==cold.momentumGeometry().layout().nativeFace(u),"template cold comparison changed actual native order");
            check(pointEqual(f.lower,e[0])&&pointEqual(f.upper,e[1]),"template target retained source face coordinates");
            if(f.area==0||(bound&&std::binary_search(warmWall->velocityDofs().begin(),warmWall->velocityDofs().end(),u)))continue;
            const long double a=e[0].x,b=e[1].x,z=e[0].y;
            if(family==0)velocity[u]=static_cast<double>(axis(l,u)==0?a*(1.L+(e[0].y+e[1].y)/2):1.L+z+z*z+(a*a+b*b)/2);
            else {const auto x=(e[0]+e[1])*.5;velocity[u]=std::sin(x.x+2*x.y)+.1*std::cos(3*x.x-x.y);}
        }
        const auto data=warm.apply(mesh,HydroNativeValues(l,velocity,std::vector<double>(l.pressureCount())),warmWall);
        const auto& cl=cold.momentumGeometry().layout();
        const auto referenceData=cold.apply(mesh,HydroNativeValues(cl,velocity,std::vector<double>(cl.pressureCount())),coldWall);
        check(data.cellDofs.size()==referenceData.cellDofs.size()&&data.faceCorrections.size()==referenceData.faceCorrections.size(),"template target cold coefficient inventory");
        for(std::size_t p=0;p<data.cellDofs.size();++p)for(int c=0;c<2;++c)for(int k=0;k<25;++k)
            peak=std::max(peak,norm(data.cellDofs[p][c][k],referenceData.cellDofs[p][c][k],1,1e-11,"template warm/cold native action"));
        for(std::size_t u=0;u<data.faceCorrections.size();++u)
            peak=std::max(peak,norm(data.faceCorrections[u],referenceData.faceCorrections[u],1,1e-11,"template warm/cold face correction"));
        const RegularReference bernstein(warm,data);
        for(int p=0;p<l.pressureCount();++p){const auto& cell=warm.cells()[p];bool identifiable=true;
            for(const auto& row:warm.rows())if(!row.zeroNormal){bool used=false;
                if(row.kind==decltype(row.kind)::Center)used=row.index==p;
                if(row.kind==decltype(row.kind)::Vertex)used=std::find(cell.vertices.begin(),cell.vertices.end(),row.index)!=cell.vertices.end();
                if(row.kind==decltype(row.kind)::Edge)used=std::find(cell.edges.begin(),cell.edges.end(),row.index)!=cell.edges.end();
                if(used&&row.p2Rank!=row.p2ColumnCount)identifiable=false;}
            for(Vector2D fraction:std::array<Vector2D,2>{{{.5,.5},{.23,.61}}}){
                const auto x=cell.lower+Vector2D{fraction.x*(cell.upper.x-cell.lower.x),fraction.y*(cell.upper.y-cell.lower.y)};
                const auto actual=regularJetValues(warm.cellJet(data,p,x)),expected=bernstein.sample(p,x);
                for(int k=0;k<6;++k)peak=std::max(peak,norm(actual[k],expected[k],1+std::abs(expected[k]),1e-11,"template independent Bernstein jet"));
                if(family==0&&!bound&&identifiable){const auto exact=physical(x);
                    for(int k=0;k<6;++k)peak=std::max(peak,norm(actual[k],exact[k],1,1e-11,"template independent identifiable physical P2"));++p2;}
            }
            if(cell.lower.x==0){const auto x=Vector2D{0,cell.lower.y+.37*(cell.upper.y-cell.lower.y)};
                const auto j=warm.cellJet(data,p,x).jet;check(j.value.x==0&&j.derivative.xy==0&&j.derivative.yx==0,"template RZ axis parity");++axes;}
        }
        // Integrate the actual normal trace against physical 2*pi*r. GL3 is
        // exact on each prepared Q4 support piece, including the radial weight.
        const long double node=std::sqrt(3.L/5),nodes[3]={-node,0,node},weights[3]={5.L/9,8.L/9,5.L/9};
        for(int u=0;u<l.velocityCount();++u){const auto& f=warm.faces()[u];if(f.area==0)continue;
            const auto e=ends(u);const int n=axis(l,u),t=1-n,p=f.negativePressure>=0?f.negativePressure:f.positivePressure;
            const auto cuts=regularCuts(warm,p,t,component(e[0],t),component(e[1],t));long double integral=0;
            for(std::size_t j=1;j<cuts.size();++j)for(int q=0;q<3;++q){auto x=e[0];(t==0?x.x:x.y)=static_cast<double>((cuts[j-1]+cuts[j])/2+(cuts[j]-cuts[j-1])*nodes[q]/2);
                const auto value=regularJetValues(warm.cellJet(data,p,x))[n];integral+=weights[q]/2*(cuts[j]-cuts[j-1])*CartesianMeasure::twoPi*x.x*value;}
            peak=std::max(peak,norm(integral,static_cast<long double>(velocity[u])*f.area,f.area*(1+std::abs(velocity[u])),1e-11,"template independent native normal mean"));++means;
            if(bound&&std::binary_search(warmWall->velocityDofs().begin(),warmWall->velocityDofs().end(),u)){
                check(velocity[u]==0,"template wall fixture is not impermeable");++walls;}
        }
    };
    const auto cubicPorts=[&](const QuadtreeMesh& mesh,const Regular& g){const auto& l=g.momentumGeometry().layout();
        const long double node=std::sqrt(3.L/5),q[3]={(1-node)/2,.5L,(1+node)/2},w[3]={5.L/18,4.L/9,5.L/18};
        for(int c=0;c<2;++c){std::vector<double> velocity(l.velocityCount());
            for(int u=0;u<l.velocityCount();++u){const auto& f=g.faces()[u];if(f.normalAxis!=c||f.area==0)continue;
                for(int k=0;k<3;++k){const auto x=f.lower+(f.upper-f.lower)*static_cast<double>(q[k]);
                    velocity[u]+=static_cast<double>(w[k]*std::pow(component(x,1-c),3)*x.x/((f.lower.x+f.upper.x)/2));}}
            const auto data=g.apply(mesh,HydroNativeValues(l,velocity,std::vector<double>(l.pressureCount())));
            for(const auto& port:g.normalPorts())if(port.component==c&&port.observedModes==4){const auto x=point(mesh,port.vertex);int owner=-1;
                for(int p=0;p<l.pressureCount();++p){const auto& cell=g.cells()[p];if(x.x>=cell.lower.x&&x.x<=cell.upper.x&&x.y>=cell.lower.y&&x.y<=cell.upper.y){owner=p;break;}}
                check(owner>=0,"template normal-port actual vertex owner");const auto actual=regularJetValues(g.cellJet(data,owner,x));const long double t=component(x,1-c);
                peak=std::max({peak,norm(actual[c],t*t*t,1,1e-11,"template independent line cubic value"),
                    norm(actual[2+2*c+1-c],3*t*t,1,1e-11,"template independent line cubic slope")});++ports;}
        }};
    auto sourceMesh=meshFor(0,0,false);const auto source=HydroMomentumGeometry::build(*sourceMesh,metricsFor(*sourceMesh,true));
    const auto sourceOwner=HydroAdvectorReconstruction::prepare(*sourceMesh,source);const auto& sourceRegular=sourceOwner.regularGeometry();
    const auto sourceLaw=HydroAdvectorBoundary::impermeable(*sourceMesh,source,wallIds(sourceRegular));
    const auto sourceBound=HydroAdvectorReconstruction::prepare(*sourceMesh,source,sourceLaw);
    compare(*sourceMesh,sourceRegular,sourceRegular,nullptr,nullptr,0);cubicPorts(*sourceMesh,sourceRegular);
    compare(*sourceMesh,sourceBound.regularGeometry(),sourceBound.regularGeometry(),&sourceLaw,&sourceLaw,1);
    check(sourceRegular.report().templateRowBuilds>0&&sourceRegular.report().templateEntryCount>0,"template source did not populate mathematical library");
    for(int scene=0;scene<3;++scene){auto mesh=meshFor(scene==1?.25:0,.25,scene==2);
        const auto warm=HydroMomentumGeometry::buildWithRegularTemplates(*mesh,metricsFor(*mesh,true),source);
        const auto cold=HydroMomentumGeometry::build(*mesh,metricsFor(*mesh,true));
        check(!warm.sameIdentity(source)&&!warm.sameIdentity(cold)&&!warm.layout().sameIdentity(source.layout())&&!warm.layout().sameIdentity(cold.layout()),"template library replaced fresh target geometry/layout");
        check(&warm.metrics().divergence().matrix()!=&source.metrics().divergence().matrix(),"template target retained source divergence owner");
        check((scene==2)==(mesh->cellIndexNum()!=sourceMesh->cellIndexNum()),"template fixture lost actual AMR topology change");
        const auto owner=HydroAdvectorReconstruction::prepare(*mesh,warm),referenceOwner=HydroAdvectorReconstruction::prepare(*mesh,cold);
        const auto& g=owner.regularGeometry();const auto& reference=referenceOwner.regularGeometry();
        const auto law=HydroAdvectorBoundary::impermeable(*mesh,warm,wallIds(g)),coldLaw=HydroAdvectorBoundary::impermeable(*mesh,cold,wallIds(reference));
        const auto bound=HydroAdvectorReconstruction::prepare(*mesh,warm,law),coldBound=HydroAdvectorReconstruction::prepare(*mesh,cold,coldLaw);
        compare(*mesh,g,reference,nullptr,nullptr,0);compare(*mesh,g,reference,nullptr,nullptr,1);
        compare(*mesh,bound.regularGeometry(),coldBound.regularGeometry(),&law,&coldLaw,1);cubicPorts(*mesh,g);
        for(const auto& row:g.rows())++kinds[static_cast<int>(row.kind)];
        const auto& report=g.report();check(report.templateRowQueries>0&&report.templatePortQueries>0,"template actual row/port requests missing");
        check(report.templateEntryCount<=report.templateEntryLimit&&report.templateRetainedBytes<=report.templateByteLimit,"template bounded library exceeded declared limits");
        if(scene==0){check(report.templateRowHits>0&&report.templatePortHits>0,"dyadic z translation did not reuse row/port mathematics");
            check(report.templateRowBuilds<reference.report().templateRowBuilds&&report.svdCount<reference.report().svdCount,"warm z translation did not avoid actual row construction/SVD");}
        if(scene==1)check(report.templateRowBuilds+report.templatePortBuilds>0,"different RZ radius reused all source templates without new construction");
        if(scene==2){check(report.templateRowHits>0,"AMR target failed to reuse unchanged mathematics");
            rejects([&]{sourceRegular.validateFor(*mesh,source);},"NativeMeshIdentity");++guards;}
        if(scene==0){const auto& l=warm.layout();const HydroNativeValues zero(l,std::vector<double>(l.velocityCount()),std::vector<double>(l.pressureCount()));
            const auto& cl=cold.layout();const HydroNativeValues foreign(cl,std::vector<double>(cl.velocityCount()),std::vector<double>(cl.pressureCount()));
            rejects([&]{g.apply(*mesh,foreign);},"layout identity");rejects([&]{g.apply(*mesh,zero,&coldLaw);},"momentum geometry instance");
            const auto old=reference.apply(*mesh,foreign);const auto x=(g.cells()[0].lower+g.cells()[0].upper)*.5;
            rejects([&]{g.cellJet(old,0,x);},"coefficient geometry identity");guards+=3;
            const auto rounding=std::fegetround();std::fesetround(FE_UPWARD);try{rejects([&]{g.apply(*mesh,zero);},"round-to-nearest");}catch(...){std::fesetround(rounding);throw;}std::fesetround(rounding);
            const auto csr=_mm_getcsr();_mm_setcsr(csr|(1U<<15));try{rejects([&]{g.apply(*mesh,zero);},"gradual underflow");}catch(...){_mm_setcsr(csr);throw;}_mm_setcsr(csr);guards+=2;}
        std::cout<<"REGULAR_TEMPLATE_BANK,"<<scene<<','<<warm.layout().pressureCount()<<','<<report.templateRowQueries<<','<<report.templateRowHits<<','<<report.templateRowBuilds<<','
            <<report.templatePortQueries<<','<<report.templatePortHits<<','<<report.templatePortBuilds<<','<<report.templateFallbacks<<','<<report.svdCount<<','<<reference.report().svdCount<<",PASS\n";
    }
    check(p2>0&&ports>0&&means>0&&axes>0&&walls>0&&guards==6&&std::all_of(kinds.begin(),kinds.end(),[](auto n){return n>0;}),"template structural/physical/guard coverage incomplete");
    std::cout<<"REGULAR_TEMPLATE_BANK_SUMMARY,"<<p2<<','<<ports<<','<<means<<','<<axes<<','<<walls<<','<<guards<<','<<kinds[0]<<','<<kinds[1]<<','<<kinds[2]<<','<<peak<<",PASS\n";
}
}

// Task52: independent ORIGINAL quadratic parent oracle on the existing six
// actual AMR fixtures. No private fast-path probe and no synthetic route map.
namespace ExactEdgeParentRoutingContract {
using Regular=LAMP::Hydrodynamics::Detail::HydroRegularAdvectorGeometry;
using Dyadic=LAMP::Geometry::Detail::FiniteVolume::Dyadic;
bool enabled=false;
std::size_t geometries=0,edgesChecked=0,parented=0,groupsChecked=0,quadraticVisits=0,largestGroup=0,neutralPairs=0,exactOpportunityGroups=0;
struct Environment
{
    std::fenv_t native{};unsigned csr=0;int error=errno;
    Environment(){check(std::fegetenv(&native)==0,"edge parent oracle cannot capture native environment");csr=_mm_getcsr();}
    void restore()const noexcept{if(std::fesetenv(&native)!=0)std::terminate();_mm_setcsr(csr);errno=error;}
    ~Environment(){restore();}
};
std::uint64_t bits(double x){std::uint64_t b=0;std::memcpy(&b,&x,sizeof b);return b;}
bool normalOrZero(double x)
{const auto b=bits(x)&UINT64_C(0x7fffffffffffffff),e=b&UINT64_C(0x7ff0000000000000);return b==0||(e!=0&&e!=UINT64_C(0x7ff0000000000000));}
bool exactDifference(double a,double b,double actual)
{
    if(!normalOrZero(a)||!normalOrZero(b)||!normalOrZero(actual))return false;
    Dyadic exact(a);exact.add(Dyadic(b),-1);return exact.compare(Dyadic(actual))==0;
}
void inspect(const Regular& g)
{
    if(!enabled)return;Environment caller;const auto& edges=g.edges();
    std::map<std::pair<int,double>,std::vector<int>> groups;
    for(std::size_t i=0;i<edges.size();++i){const auto& e=edges[i];check(e.normalAxis==0||e.normalAxis==1,"edge parent oracle invalid public axis");groups[{e.normalAxis,component(e.lower,e.normalAxis)}].push_back(static_cast<int>(i));}
    std::vector<int> expected(edges.size(),-1);bool allNeutral=true;
    for(const auto& group:groups){bool neutral=true;largestGroup=std::max(largestGroup,group.second.size());
        for(int i:group.second){const auto& e=edges[i];const auto difference=e.upper-e.lower;
            neutral=neutral&&exactDifference(e.upper.x,e.lower.x,difference.x)&&exactDifference(e.upper.y,e.lower.y,difference.y);
            const int t=1-e.normalAxis;double longest=component(e.upper-e.lower,t);
            for(int j:group.second){const auto& other=edges[j];const double length=component(other.upper-other.lower,t);
                if(length>longest&&component(other.lower,t)<=component(e.lower,t)&&component(other.upper,t)>=component(e.upper,t)){expected[i]=j;longest=length;}
                ++quadraticVisits;}}
        if(neutral)++exactOpportunityGroups;allNeutral=allNeutral&&neutral;++groupsChecked;
    }
    for(std::size_t i=0;i<edges.size();++i){check(edges[i].parent==expected[i],"edge parent differs from original ordered quadratic containment search");if(expected[i]>=0)++parented;++edgesChecked;}
#if defined(_MSC_VER) && defined(_M_X64) && !defined(__FAST_MATH__) && (!defined(_M_FP_FAST) || !_M_FP_FAST)
    unsigned control=0;const bool qualified=std::fegetround()==FE_TONEAREST&&_controlfp_s(&control,0,0)==0&&
        (control&_MCW_EM)==_MCW_EM&&(control&_MCW_RC)==_RC_NEAR&&(_mm_getcsr()&0x1f80U)==0x1f80U&&(_mm_getcsr()&((1U<<15)|(1U<<6)))==0;
    if(qualified&&allNeutral){
        std::vector<int> repeated(edges.size(),-1);caller.restore();check(std::feclearexcept(FE_ALL_EXCEPT)==0,"edge parent oracle cannot clear hardware status");
        _mm_setcsr(_mm_getcsr()&~0x3fU);errno=E2BIG;const Environment before;
        for(const auto& group:groups)for(int i:group.second){const auto& e=edges[i];const int t=1-e.normalAxis;double longest=component(e.upper-e.lower,t);
            for(int j:group.second){const auto& other=edges[j];const double length=component(other.upper-other.lower,t);
                if(length>longest&&component(other.lower,t)<=component(e.lower,t)&&component(other.upper,t)>=component(e.upper,t)){repeated[i]=j;longest=length;}}}
        const Environment after;
        check(std::memcmp(&before.native,&after.native,sizeof after.native)==0&&before.csr==after.csr&&before.error==after.error,
            "original exact edge search changed full native/raw FP status or errno");
        check(repeated==expected,"repeated original parent oracle changed routing");++neutralPairs;
    }
#endif
    ++geometries;
}
void finish()
{
    check(geometries==12&&edgesChecked>0&&parented>0&&groupsChecked>0&&largestGroup>1&&exactOpportunityGroups>0,
        "edge parent original six-scene topology/oracle coverage incomplete");
#if defined(_MSC_VER) && defined(_M_X64) && !defined(__FAST_MATH__) && (!defined(_M_FP_FAST) || !_M_FP_FAST)
    check(neutralPairs==geometries,"edge parent six dyadic AMR scenes lost neutral full-environment witness");
#endif
    std::cout<<"REGULAR_EDGE_PARENT_ROUTING_PASS,actual_AMR_scenes,6,public_geometries,"<<geometries<<",edges,"<<edgesChecked
        <<",nontrivial_parents,"<<parented<<",collinear_groups,"<<groupsChecked<<",largest_group,"<<largestGroup
        <<",original_quadratic_visits,"<<quadraticVisits<<",exact_difference_opportunity_groups,"<<exactOpportunityGroups
        <<",original_loop_full_FP_errno_neutral_pairs,"<<neutralPairs
        <<",original_P2_axis_ports_native_means_source_guards_retained,private_fastpath_hit_unobserved_no_speed_claim\n";
}
}

// Task35: actual overlap preparation reuse, independently checked on fresh target owners.
namespace Task35RegularRemappedPreparation {
using Regular=LAMP::Hydrodynamics::Detail::HydroRegularAdvectorGeometry;
std::size_t cases=0,rows=0,ports=0,jets=0,p2=0,means=0,axisSamples=0,guards=0;
std::array<std::size_t,3> rowKinds{};
long double peak=0;
void number(long double a,long double b,const char* why)
{peak=std::max(peak,norm(a,b,1+std::abs(b),1e-11L,why));}
std::unique_ptr<QuadtreeMesh> meshFor(int scene,bool refined)
{
    DictionaryMesh d;d.Nx=d.Ny=8;d.maxCellNum=16384;d.maxInitialIterations=0;
    const double r=scene==2?.25:0,height=scene==0?2:1;
    d.minCoord={r,0};d.maxCoord={r+1,height};auto mesh=std::make_unique<QuadtreeMesh>(d);
    DictionaryMesh::RefineBox center;center.level=1;center.min={r+.125,.125*height};center.max={r+.875,.875*height};
    std::vector<DictionaryMesh::RefineBox> boxes{center};
    if(refined){auto corner=center;corner.min={r,0};corner.max={r+.125,.125*height};boxes.push_back(corner);}
    mesh->applyRefineBoxes(boxes);double lo=1,hi=0;
    for(int p=0;p<mesh->cellIndexNum();++p){const auto h=mesh->cells(mesh->cellIndexToLabel[p]).size().x;lo=std::min(lo,h);hi=std::max(hi,h);}
    check(hi==2*lo,"Task35 fixture is not actual mixed-level 2:1");return mesh;
}
std::vector<int> walls(const Regular& g)
{std::vector<int> ids;for(std::size_t u=0;u<g.faces().size();++u)if(g.faces()[u].external&&g.faces()[u].area>0)ids.push_back(static_cast<int>(u));return ids;}
std::array<long double,6> polynomial(Vector2D x)
{return {x.x*(1.L+x.y),1.L+x.y+x.x*x.x+x.y*x.y,1.L+x.y,x.x,2.L*x.x,1.L+2.L*x.y};}
std::vector<double> values(const Regular& g,int family,const HydroAdvectorBoundary* wall)
{
    std::vector<double> v(g.faces().size());const bool rz=g.momentumGeometry().layout().domain().axisymmetric();
    const long double a=std::sqrt(3.L/5),q[3]={(1-a)/2,.5L,(1+a)/2},w[3]={5.L/18,4.L/9,5.L/18};
    for(std::size_t u=0;u<v.size();++u){const auto& f=g.faces()[u];
        if(f.area==0||(wall&&std::binary_search(wall->velocityDofs().begin(),wall->velocityDofs().end(),static_cast<int>(u))))continue;
        if(family){const auto x=(f.lower+f.upper)*.5;v[u]=std::sin(1.3*x.x+2.1*x.y)+.17*std::cos(3*x.x-x.y);continue;}
        long double integral=0;const int t=1-f.normalAxis;const long double length=component(f.upper-f.lower,t);
        for(int k=0;k<3;++k){const auto x=f.lower+(f.upper-f.lower)*static_cast<double>(q[k]);
            integral+=w[k]*length*(rz?CartesianMeasure::twoPi*x.x:1.L)*polynomial(x)[f.normalAxis];}
        v[u]=static_cast<double>(integral/f.area);
    }return v;
}
void topology(const Regular& a,const Regular& b)
{
    ExactEdgeParentRoutingContract::inspect(a);ExactEdgeParentRoutingContract::inspect(b);
    check(!a.sameIdentity(b),"Task35 warm/cold Regular owners alias");
    check(a.cells().size()==b.cells().size()&&a.faces().size()==b.faces().size()&&a.rows().size()==b.rows().size()&&a.normalPorts().size()==b.normalPorts().size(),"Task35 target topology inventory");
    for(std::size_t p=0;p<a.cells().size();++p){const auto& x=a.cells()[p];const auto& y=b.cells()[p];
        check(x.nativeCell==y.nativeCell&&pointEqual(x.lower,y.lower)&&pointEqual(x.upper,y.upper)&&x.vertices==y.vertices&&x.edges==y.edges&&x.faces==y.faces&&x.supportKnots==y.supportKnots,"Task35 fresh target cell/support topology");}
    for(std::size_t u=0;u<a.faces().size();++u){const auto& x=a.faces()[u];const auto& y=b.faces()[u];
        check(x.normalAxis==y.normalAxis&&x.negativePressure==y.negativePressure&&x.positivePressure==y.positivePressure&&pointEqual(x.lower,y.lower)&&pointEqual(x.upper,y.upper)&&x.area==y.area&&x.external==y.external&&x.supportWidth==y.supportWidth,"Task35 fresh target actual face geometry");}
    for(std::size_t i=0;i<a.rows().size();++i){const auto& x=a.rows()[i];const auto& y=b.rows()[i];
        check(x.kind==y.kind&&x.index==y.index&&x.component==y.component&&x.rank==y.rank&&x.rings==y.rings&&x.degree==y.degree&&x.p2Rank==y.p2Rank&&x.p2ColumnCount==y.p2ColumnCount&&x.zeroNormal==y.zeroNormal&&pointEqual(x.point,y.point)&&pointEqual(x.scale,y.scale)&&x.columns==y.columns&&x.support==y.support,"Task35 complete ordered cold row routing/rank");
        for(int k=0;k<4;++k){check(x.weights[k].size()==y.weights[k].size(),"Task35 row weights inventory");for(std::size_t j=0;j<x.weights[k].size();++j)number(x.weights[k][j],y.weights[k][j],"Task35 cold row coefficient");}++rows;++rowKinds[static_cast<int>(x.kind)];}
    for(std::size_t i=0;i<a.normalPorts().size();++i){const auto& x=a.normalPorts()[i];const auto& y=b.normalPorts()[i];
        check(x.vertex==y.vertex&&x.component==y.component&&x.observedModes==y.observedModes&&x.support==y.support,"Task35 complete ordered cold normal port");
        for(int k=0;k<2;++k){check(x.weights[k].size()==y.weights[k].size(),"Task35 port weights inventory");for(std::size_t j=0;j<x.weights[k].size();++j)number(x.weights[k][j],y.weights[k][j],"Task35 cold port coefficient");}++ports;}
}
void action(const QuadtreeMesh& mesh,const Regular& warm,const Regular& cold,int family,
    const HydroAdvectorBoundary* wall=nullptr,const HydroAdvectorBoundary* coldWall=nullptr)
{
    const auto& l=warm.momentumGeometry().layout();const auto& cl=cold.momentumGeometry().layout();const auto v=values(warm,family,wall);
    const auto a=warm.apply(mesh,HydroNativeValues(l,v,std::vector<double>(l.pressureCount())),wall);
    const auto b=cold.apply(mesh,HydroNativeValues(cl,v,std::vector<double>(cl.pressureCount())),coldWall);
    check(a.cellDofs.size()==b.cellDofs.size()&&a.faceCorrections.size()==b.faceCorrections.size(),"Task35 cold action inventory");
    for(std::size_t p=0;p<a.cellDofs.size();++p)for(int c=0;c<2;++c)for(int k=0;k<25;++k)number(a.cellDofs[p][c][k],b.cellDofs[p][c][k],"Task35 cold actual cell coefficient");
    for(std::size_t u=0;u<a.faceCorrections.size();++u)number(a.faceCorrections[u],b.faceCorrections[u],"Task35 cold actual face correction");
    const RegularReference independent(warm,a);const bool rz=l.domain().axisymmetric();
    for(int p=0;p<l.pressureCount();++p){const auto& cell=warm.cells()[p];bool identifiable=true;
        for(const auto& row:warm.rows())if(!row.zeroNormal){bool used=false;
            if(row.kind==decltype(row.kind)::Center)used=row.index==p;
            if(row.kind==decltype(row.kind)::Vertex)used=std::find(cell.vertices.begin(),cell.vertices.end(),row.index)!=cell.vertices.end();
            if(row.kind==decltype(row.kind)::Edge)used=std::find(cell.edges.begin(),cell.edges.end(),row.index)!=cell.edges.end();
            if(used&&row.p2Rank!=row.p2ColumnCount)identifiable=false;}
        for(Vector2D fraction:std::array<Vector2D,2>{{{.37,.43},{.71,.62}}}){const auto x=cell.lower+Vector2D{fraction.x*(cell.upper.x-cell.lower.x),fraction.y*(cell.upper.y-cell.lower.y)};
            const auto actual=regularJetValues(warm.cellJet(a,p,x)),ref=independent.sample(p,x),other=regularJetValues(cold.cellJet(b,p,x));
            for(int k=0;k<6;++k){number(actual[k],ref[k],"Task35 independent Bernstein jet");number(actual[k],other[k],"Task35 cold target jet");}
            ++jets;if(family==0&&!wall&&identifiable){const auto exact=polynomial(x);for(int k=0;k<6;++k)number(actual[k],exact[k],"Task35 independent identifiable physical P2");++p2;}}
        if(rz&&cell.lower.x==0){const auto j=warm.cellJet(a,p,{0,cell.lower.y+.37*(cell.upper.y-cell.lower.y)}).jet;
            check(j.value.x==0&&j.derivative.xy==0&&j.derivative.yx==0,"Task35 true-axis parity");++axisSamples;}
    }
    const long double node=std::sqrt(3.L/5),nodes[3]={-node,0,node},weights[3]={5.L/9,8.L/9,5.L/9};
    for(int u=0;u<l.velocityCount();++u){const auto& f=warm.faces()[u];if(f.area==0)continue;const int n=f.normalAxis,t=1-n,p=f.negativePressure>=0?f.negativePressure:f.positivePressure;
        const auto cuts=regularCuts(warm,p,t,component(f.lower,t),component(f.upper,t));long double integral=0;
        for(std::size_t j=1;j<cuts.size();++j)for(int q=0;q<3;++q){auto x=f.lower;(t==0?x.x:x.y)=static_cast<double>((cuts[j-1]+cuts[j])/2+(cuts[j]-cuts[j-1])*nodes[q]/2);
            integral+=weights[q]/2*(cuts[j]-cuts[j-1])*(rz?CartesianMeasure::twoPi*x.x:1.L)*regularJetValues(warm.cellJet(a,p,x))[n];}
        peak=std::max(peak,norm(integral,static_cast<long double>(v[u])*f.area,f.area*(1+std::abs(v[u])),1e-11L,"Task35 independent physical native normal mean"));++means;
    }
}
bool sameReport(const HydroRegularPreparationReuseReport& a,const HydroRegularPreparationReuseReport& b)
{return a.requested==b.requested&&a.enabled==b.enabled&&a.targetDatums==b.targetDatums&&a.reusedRows==b.reusedRows&&a.reusedPorts==b.reusedPorts&&a.recomputedRows==b.recomputedRows&&a.recomputedPorts==b.recomputedPorts&&a.stencilRejectedRows==b.stencilRejectedRows&&a.portSpanRejectedRows==b.portSpanRejectedRows&&a.columnRejectedRows==b.columnRejectedRows&&a.fpRejectedRows==b.fpRejectedRows&&a.retainedMetadataBytes==b.retainedMetadataBytes;}
void scenario(int scene,bool coarsen)
{
    const bool rz=scene!=0;auto sourceMesh=meshFor(scene,coarsen),targetMesh=meshFor(scene,!coarsen);
    const auto source=HydroMomentumGeometry::build(*sourceMesh,metricsFor(*sourceMesh,rz));
    const auto sourceOwner=HydroAdvectorReconstruction::prepare(*sourceMesh,source);const auto& old=sourceOwner.regularGeometry();
    check(!old.usesCompatibleTensor()&&!old.usesVirtualCompatibleTensor(),"Task35 positive source must exercise general Regular");
    const auto targetMetrics=metricsFor(*targetMesh,rz);
    const auto overlap=LAMP::Mesh::NativeLeafOverlap::build(*sourceMesh,source.layout().domain(),*targetMesh,targetMetrics.layout().domain());
    std::size_t changed=0;for(const auto& piece:overlap.pieces())if(piece.kind!=LAMP::Mesh::NativeLeafOverlapPiece::Kind::Unchanged)++changed;
    check(changed>0&&(coarsen?targetMesh->cellIndexNum()<sourceMesh->cellIndexNum():targetMesh->cellIndexNum()>sourceMesh->cellIndexNum()),"Task35 requires a genuine refine/coarsen overlap");
    const auto warm=HydroMomentumGeometry::buildWithRemappedRegularPreparation(*sourceMesh,old,*targetMesh,targetMetrics,overlap);
    const auto pending=warm.regularPreparationReuseReport();check(pending.requested&&!pending.enabled,"Task35 remap request must be pending before actual consumption");
    const auto cold=HydroMomentumGeometry::build(*targetMesh,targetMetrics);
    check(!warm.sameIdentity(source)&&!warm.sameIdentity(cold)&&warm.layout().sameIdentity(targetMetrics.layout())&& &warm.metrics().divergence().matrix()== &targetMetrics.divergence().matrix()&& &warm.metrics().divergence().matrix()!= &source.metrics().divergence().matrix(),"Task35 fresh authoritative target layout/D");
    const auto owner=HydroAdvectorReconstruction::prepare(*targetMesh,warm),reference=HydroAdvectorReconstruction::prepare(*targetMesh,cold);
    const auto& g=owner.regularGeometry();const auto& ref=reference.regularGeometry();
    check(!g.usesCompatibleTensor()&&!g.usesVirtualCompatibleTensor(),"Task35 positive target must exercise general Regular");
    const auto r=warm.regularPreparationReuseReport();check(r.requested&&r.enabled&&r.reusedRows>0&&r.recomputedRows>0&&r.reusedPorts>0,"Task35 positive overlap omitted no actual row/port work");
    check(r.reusedRows+r.recomputedRows==g.rows().size()&&r.reusedPorts+r.recomputedPorts==g.normalPorts().size()&&r.targetDatums>0&&r.targetDatums<=g.rows().size()&&r.retainedMetadataBytes>0,"Task35 actual construction counters/inventory");
    topology(g,ref);action(*targetMesh,g,ref,0);action(*targetMesh,g,ref,1);
    const auto wall=HydroAdvectorBoundary::impermeable(*targetMesh,warm,walls(g)),coldWall=HydroAdvectorBoundary::impermeable(*targetMesh,cold,walls(ref));
    const auto bound=HydroAdvectorReconstruction::prepare(*targetMesh,warm,wall),coldBound=HydroAdvectorReconstruction::prepare(*targetMesh,cold,coldWall);
    action(*targetMesh,bound.regularGeometry(),coldBound.regularGeometry(),1,&wall,&coldWall);
    const auto again=HydroAdvectorReconstruction::prepare(*targetMesh,warm,wall);(void)again;
    check(sameReport(r,warm.regularPreparationReuseReport()),"Task35 one-shot construction report changed on cached boundary prepare");
    if(scene==0&&!coarsen){
        rejects([&]{HydroMomentumGeometry::buildWithRemappedRegularPreparation(*targetMesh,old,*targetMesh,targetMetrics,overlap);},"NativeMeshIdentity");
        const auto wrongOverlap=LAMP::Mesh::NativeLeafOverlap::build(*sourceMesh,source.layout().domain(),*sourceMesh,source.layout().domain());
        rejects([&]{HydroMomentumGeometry::buildWithRemappedRegularPreparation(*sourceMesh,old,*targetMesh,targetMetrics,wrongOverlap);},"NativeMeshIdentity");
        const auto& l=warm.layout();const auto& cl=cold.layout();const HydroNativeValues zero(l,std::vector<double>(l.velocityCount()),std::vector<double>(l.pressureCount())),foreign(cl,std::vector<double>(cl.velocityCount()),std::vector<double>(cl.pressureCount()));
        // The cold oracle deliberately shares the actual target layout; a new actual layout is the negative.
        const auto other=HydroMomentumGeometry::build(*targetMesh,metricsFor(*targetMesh,rz));const auto& ol=other.layout();const HydroNativeValues foreignValues(ol,std::vector<double>(ol.velocityCount()),std::vector<double>(ol.pressureCount()));
        rejects([&]{g.apply(*targetMesh,foreignValues);},"layout identity");
        const auto wrongData=ref.apply(*targetMesh,foreign);const auto x=(g.cells()[0].lower+g.cells()[0].upper)*.5;
        rejects([&]{g.cellJet(wrongData,0,x);},"coefficient geometry identity");guards+=4;
        const int rounding=std::fegetround();check(std::fesetround(FE_UPWARD)==0,"Task35 cannot set unsupported rounding");
        try{rejects([&]{g.apply(*targetMesh,zero);},"round-to-nearest");rejects([&]{HydroAdvectorReconstruction::prepare(*targetMesh,warm);},"round-to-nearest");}
        catch(...){std::fesetround(rounding);throw;}check(std::fesetround(rounding)==0,"Task35 cannot restore rounding");guards+=2;
        const auto csr=_mm_getcsr();_mm_setcsr(csr|(1U<<15));try{rejects([&]{g.apply(*targetMesh,zero);},"gradual underflow");}catch(...){_mm_setcsr(csr);throw;}_mm_setcsr(csr);++guards;
    }
    ++cases;std::cout<<"REGULAR_REMAP_PREPARATION,"<<scene<<','<<coarsen<<','<<sourceMesh->cellIndexNum()<<','<<targetMesh->cellIndexNum()<<','<<r.targetDatums<<','<<r.reusedRows<<','<<r.recomputedRows<<','<<r.reusedPorts<<','<<r.recomputedPorts<<','<<r.stencilRejectedRows<<','<<r.portSpanRejectedRows<<','<<r.columnRejectedRows<<','<<r.fpRejectedRows<<','<<r.retainedMetadataBytes<<",PASS\n";
}
void run()
{
    for(int scene=0;scene<3;++scene)for(bool coarsen:{false,true})scenario(scene,coarsen);
    check(cases==6&&rows>0&&ports>0&&p2>0&&means>0&&axisSamples>0&&guards==7&&std::all_of(rowKinds.begin(),rowKinds.end(),[](auto n){return n>0;}),"Task35 missing physical/cold/guard coverage");
    std::cout<<"REGULAR_REMAP_PREPARATION_PASS,cases,"<<cases<<",rows,"<<rows<<",ports,"<<ports<<",jets,"<<jets<<",independent_P2,"<<p2<<",native_means,"<<means<<",axis_samples,"<<axisSamples<<",guards,"<<guards<<",peak,"<<peak<<"\n";
}
} // namespace Task35RegularRemappedPreparation

// F11 additive namespace: insert after Task35RegularRemappedPreparation and
// before the existing main. Add only the selector branch shown at the end.
namespace CompletedRegularDependencyContract {
using Regular=LAMP::Hydrodynamics::Detail::HydroRegularAdvectorGeometry;
using Environment=ExactEdgeParentRoutingContract::Environment;
struct Prepared {
    std::unique_ptr<QuadtreeMesh> mesh;
    HydroMomentumGeometry geometry;
    HydroAdvectorReconstruction owner;
};
std::unique_ptr<QuadtreeMesh> meshFor(int ordinal,bool target){
    DictionaryMesh d;d.Nx=d.Ny=8;d.maxCellNum=16384;d.maxInitialIterations=0;
    // Bounded R, unique odd dyadic offsets. Unlike a leaf-width translation,
    // these do not deliberately recycle an earlier radial rho signature.
    const double radial=.25+std::ldexp(static_cast<double>(2*ordinal+1),-20);
    d.minCoord={radial,0};d.maxCoord={radial+1,1};auto mesh=std::make_unique<QuadtreeMesh>(d);
    DictionaryMesh::RefineBox center;center.level=1;center.min={radial+.125,.125};center.max={radial+.875,.875};
    std::vector<DictionaryMesh::RefineBox> boxes{center};
    if(target){auto corner=center;corner.min={radial,0};corner.max={radial+.125,.125};boxes.push_back(corner);}
    mesh->applyRefineBoxes(boxes);double smallest=1,largest=0;
    for(int p=0;p<mesh->cellIndexNum();++p){const double h=mesh->cells(mesh->cellIndexToLabel[p]).size().x;
        smallest=std::min(smallest,h);largest=std::max(largest,h);}
    check(largest==2*smallest,"F11 fixture lost actual RZ2:1 leaves");return mesh;
}
Prepared prepare(int ordinal,const HydroMomentumGeometry& bankOwner){
    auto mesh=meshFor(ordinal,false);const auto geometry=HydroMomentumGeometry::buildWithRegularTemplates(*mesh,metricsFor(*mesh,true),bankOwner);
    const auto owner=HydroAdvectorReconstruction::prepare(*mesh,geometry);
    check(!owner.regularGeometry().usesCompatibleTensor()&&!owner.regularGeometry().usesVirtualCompatibleTensor(),"F11 source must execute general normalized Regular path");
    return {std::move(mesh),geometry,owner};
}
bool sameStock(const HydroRegularTemplateStatistics& a,const HydroRegularTemplateStatistics& b){
    return a.entries==b.entries&&a.rowEntries==b.rowEntries&&a.routingEntries==b.routingEntries&&a.retainedBytes==b.retainedBytes;
}
bool capacityMiss(const HydroRegularTemplateStatistics& before,const HydroRegularTemplateStatistics& after,
    const LAMP::Hydrodynamics::Detail::HydroRegularPreparationReport& report){
    return after.entries<after.entryLimit&&after.retainedBytes<after.byteLimit&&sameStock(before,after)
        &&report.templateRowBuilds>0&&report.templatePublishedEntries==0
        &&report.templateRowHits==0&&report.templatePortHits==0&&report.templateRoutingHits==0;
}
void prime(){check(std::feclearexcept(FE_ALL_EXCEPT)==0&&std::feraiseexcept(FE_INEXACT)==0,"F11 cannot prime actual native flags");
    _mm_setcsr(_mm_getcsr()|0x20u);errno=EDOM;}
void sameFp(const Environment& a,const Environment& b){
    check(std::memcmp(&a.native,&b.native,sizeof a.native)==0&&a.csr==b.csr&&a.error==b.error,"F11 changed complete native/raw FP/errno");
}
void run(){
#if defined(_MSC_VER) && defined(_M_X64) && !defined(__FAST_MATH__) && (!defined(_M_FP_FAST) || !_M_FP_FAST)
    Environment caller;check(std::fesetround(FE_TONEAREST)==0,"F11 requires nearest mode");_mm_setcsr(_mm_getcsr()&~0x8040u);
    unsigned control=0;check(_controlfp_s(&control,0,0)==0&&(control&_MCW_EM)==_MCW_EM&&(_mm_getcsr()&0x1f80u)==0x1f80u,"F11 requires masked original controls");
    auto anchorMesh=meshFor(0,false);const auto anchor=HydroMomentumGeometry::build(*anchorMesh,metricsFor(*anchorMesh,true));
    constexpr int maximumFillers=256;int stopped=0;std::size_t actualSvdRequests=0;
    for(int ordinal=1;ordinal<=maximumFillers;++ordinal){
        const auto before=anchor.regularTemplateStatistics();const auto filler=prepare(ordinal,anchor);
        const auto after=anchor.regularTemplateStatistics();const auto& r=filler.owner.regularGeometry().report();actualSvdRequests+=r.svdCount;
        check(after.entries<=after.entryLimit&&after.retainedBytes<=after.byteLimit,"F11 actual bank exceeded original declared limits");
        std::cout<<"F11_BANK_FILL,"<<ordinal<<",entries,"<<after.entries<<",bytes,"<<after.retainedBytes
            <<",row_builds,"<<r.templateRowBuilds<<",published,"<<r.templatePublishedEntries<<",svd_requests,"<<r.svdCount<<'\n';
        if(capacityMiss(before,after,r)){stopped=ordinal;break;}
        check(after.entries<after.entryLimit&&after.retainedBytes<after.byteLimit,
            "F11 witness unavailable: exact full cap/entry limit reached before still-room publish refusal");
    }
    check(stopped>0,"F11 witness unavailable: no real still-room publish refusal within 256 actual fillers");
    // Do not reuse a filler's rho or fabricate a saturated cache. A further
    // fresh source must independently show zero old-key hits and zero growth.
    const auto beforeSource=anchor.regularTemplateStatistics();auto source=prepare(stopped+1,anchor);
    const auto afterSource=anchor.regularTemplateStatistics();const auto& old=source.owner.regularGeometry();const auto& sr=old.report();
    check(capacityMiss(beforeSource,afterSource,sr),"F11 fresh-rho source did not independently witness uncached completed rows");
    std::cout<<"F11_SATURATED_SOURCE,ordinal,"<<stopped+1<<",entries,"<<afterSource.entries<<",entry_limit,"<<afterSource.entryLimit
        <<",bytes,"<<afterSource.retainedBytes<<",byte_limit,"<<afterSource.byteLimit<<",row_builds,"<<sr.templateRowBuilds
        <<",port_builds,"<<sr.templatePortBuilds<<",published,"<<sr.templatePublishedEntries<<",old_row_hits,"<<sr.templateRowHits
        <<",old_port_hits,"<<sr.templatePortHits<<",old_route_hits,"<<sr.templateRoutingHits<<std::endl;
    auto targetMesh=meshFor(stopped+1,true);const auto targetMetrics=metricsFor(*targetMesh,true);
    const auto overlap=LAMP::Mesh::NativeLeafOverlap::build(*source.mesh,source.geometry.layout().domain(),*targetMesh,targetMetrics.layout().domain());
    std::size_t changed=0,unchanged=0;for(const auto& p:overlap.pieces())
        if(p.kind==LAMP::Mesh::NativeLeafOverlapPiece::Kind::Unchanged)++unchanged;else ++changed;
    check(changed>0&&unchanged>0&&targetMesh->cellIndexNum()>source.mesh->cellIndexNum(),"F11 requires genuine local refine plus unchanged native cells");
    const auto warm=HydroMomentumGeometry::buildWithRemappedRegularPreparation(*source.mesh,old,*targetMesh,targetMetrics,overlap);
    const auto cold=HydroMomentumGeometry::buildWithRegularTemplates(*targetMesh,targetMetrics,source.geometry);
    const auto fresh=HydroMomentumGeometry::build(*targetMesh,targetMetrics);
    check(!warm.sameIdentity(source.geometry)&&!warm.sameIdentity(cold)&&!warm.sameIdentity(fresh)
        &&warm.layout().sameIdentity(targetMetrics.layout())&& &warm.metrics().divergence().matrix()== &targetMetrics.divergence().matrix(),"F11 must retain fresh target geometry/layout/D ownership");
    prime();const Environment entry;
    const auto warmOwner=HydroAdvectorReconstruction::prepare(*targetMesh,warm);const Environment afterWarm;
    entry.restore();const auto coldOwner=HydroAdvectorReconstruction::prepare(*targetMesh,cold);const Environment afterCold;sameFp(afterWarm,afterCold);
    entry.restore();const auto freshOwner=HydroAdvectorReconstruction::prepare(*targetMesh,fresh);const Environment afterFresh;sameFp(afterWarm,afterFresh);
    const auto& g=warmOwner.regularGeometry();const auto& c=coldOwner.regularGeometry();const auto& f=freshOwner.regularGeometry();
    const auto reuse=warm.regularPreparationReuseReport();
    std::cout<<"F11_COMPLETED_DEPENDENCY_OBSERVATION,enabled,"<<reuse.enabled<<",reused_rows,"<<reuse.reusedRows
        <<",recomputed_rows,"<<reuse.recomputedRows<<",reused_ports,"<<reuse.reusedPorts<<",recomputed_ports,"<<reuse.recomputedPorts<<std::endl;
    // Old publication-coupled code must fail here: the independently verified
    // fresh source had no retained-entry hits and no successful publications.
    check(reuse.requested&&reuse.enabled&&reuse.reusedRows>0&&reuse.reusedPorts>0&&reuse.recomputedRows>0,
        "F11 uncached completed source rows/ports were denied actual local-AMR reuse");
    check(reuse.reusedRows+reuse.recomputedRows==g.rows().size()&&reuse.reusedPorts+reuse.recomputedPorts==g.normalPorts().size()
        &&reuse.retainedMetadataBytes>0,"F11 reuse counters do not cover actual target inventory");
    check(!g.usesCompatibleTensor()&&!g.usesVirtualCompatibleTensor(),"F11 target bypassed general Regular path");
    Task35RegularRemappedPreparation::topology(g,c);Task35RegularRemappedPreparation::topology(g,f);
    for(int family:{0,1}){Task35RegularRemappedPreparation::action(*targetMesh,g,c,family);Task35RegularRemappedPreparation::action(*targetMesh,g,f,family);}
    const auto wall=HydroAdvectorBoundary::impermeable(*targetMesh,warm,Task35RegularRemappedPreparation::walls(g));
    const auto coldWall=HydroAdvectorBoundary::impermeable(*targetMesh,cold,Task35RegularRemappedPreparation::walls(c));
    const auto freshWall=HydroAdvectorBoundary::impermeable(*targetMesh,fresh,Task35RegularRemappedPreparation::walls(f));
    const auto bound=HydroAdvectorReconstruction::prepare(*targetMesh,warm,wall);
    const auto coldBound=HydroAdvectorReconstruction::prepare(*targetMesh,cold,coldWall);
    const auto freshBound=HydroAdvectorReconstruction::prepare(*targetMesh,fresh,freshWall);
    Task35RegularRemappedPreparation::action(*targetMesh,bound.regularGeometry(),coldBound.regularGeometry(),1,&wall,&coldWall);
    Task35RegularRemappedPreparation::action(*targetMesh,bound.regularGeometry(),freshBound.regularGeometry(),1,&wall,&freshWall);
    check(Task35RegularRemappedPreparation::sameReport(reuse,warm.regularPreparationReuseReport()),"F11 later cached boundary preparation rewrote one-shot reuse receipt");
    const auto finalStock=anchor.regularTemplateStatistics();
    check(finalStock.entries<=finalStock.entryLimit&&finalStock.retainedBytes<=finalStock.byteLimit,"F11 preparation exceeded original storage contract");
    check(Task35RegularRemappedPreparation::p2>0&&Task35RegularRemappedPreparation::means>0&&Task35RegularRemappedPreparation::ports>0,
        "F11 omitted original P2/native-mean/ordered-port oracle coverage");
    std::cout<<"REGULAR_COMPLETED_DEPENDENCY_PASS,actual_fillers,"<<stopped<<",fill_svd_requests,"<<actualSvdRequests
        <<",source_fresh_rho_uncached_completed,1,shared_bank_cold_and_independent_fresh,1,local_overlap_changed,"<<changed
        <<",unchanged,"<<unchanged<<",reused_rows,"<<reuse.reusedRows<<",reused_ports,"<<reuse.reusedPorts
        <<",recomputed_rows,"<<reuse.recomputedRows<<",FP_pairs,2,old_Task35_P2_mean_Bernstein_wall_gates,original"
        <<",final_bank_entries,"<<finalStock.entries<<",final_bank_bytes,"<<finalStock.retainedBytes
        <<",no_capacity_override_or_hook,no_time_or_speed_claim\n";
#else
    std::cout<<"REGULAR_COMPLETED_DEPENDENCY_UNAVAILABLE,requires_original_strict_MSVC_x64_template_path\n";
#endif
}
}

int main(int argc,char** argv)
{
    try {
        std::cout<<std::setprecision(17);
        if(argc==2 && std::string(argv[1])=="--regular-completed-dependency-only"){CompletedRegularDependencyContract::run();return 0;}
        if(argc==2 && std::string(argv[1])=="--regular-edge-parent-routing-only"){ExactEdgeParentRoutingContract::enabled=true;Task35RegularRemappedPreparation::run();ExactEdgeParentRoutingContract::finish();return 0;}
        if(argc==2 && std::string(argv[1])=="--regular-remapped-preparation-only"){Task35RegularRemappedPreparation::run();return 0;}
        if(argc==2 && std::string(argv[1])=="--parallel-mass-only"){
            exactRestrictionAccumulationWitness();
            guardedTransverseEngineWitness();
            const auto filteredFallbacks=filteredTransverseFallbackWitness();std::size_t filteredSuccesses=0;long double filteredPeak=0;
            for(bool rz:{false,true}){
                auto mesh=makeMesh(4,1);Case shape;shape.pureFirstJet=true;shape.radial=rz;shape.radius=rz?.34:.19;
                const auto geometry=HydroMomentumGeometry::build(*mesh,metricsFor(*mesh,rz));
                const auto phase=HydroPhaseGeometry::build(*mesh,geometry.metrics(),initial(*mesh,shape));
                filteredTransverseSourceWitness(*mesh,phase.publication(),rz,filteredSuccesses,filteredPeak);
                parallelMassContract(*mesh,geometry,phase,shape,true);
            }
            overflow();
            std::cout<<"FILTERED_TRANSVERSE_Q4,"<<filteredSuccesses<<','<<filteredFallbacks<<','<<filteredPeak<<",PASS\n";
            return 0;
        }
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
        if(argc==2 && std::string(argv[1])=="--regular-template-bank-only"){
            regularTemplateBankContract();return 0;
        }
        if(argc==2 && std::string(argv[1])=="--p3-source-only"){
            regularCubicSourceContract();return 0;
        }
        if(argc==2 && std::string(argv[1])=="--regular-source-only"){
            spatialAdvectorContract();advectorGuards();spatialAdvectorJetContract();cellJetGuards();
            regularPreparationContract();normalFluxPortContract();localSourceStencilContract();completeNormalPolynomialContract();
            regularEvaluationContract();regularEvaluationGuards();regularLowerSpaceContract();regularIdentifiableConstantContract();
            regularCubicSourceContract();return 0;
        }
        check(argc==1,"unexpected H13 contract argument");
        for(bool rz:{false,true})for(int amr:{0,1,2})for(int selection:{0,1,2,3}) {
            auto mesh=makeMesh(4,amr);const auto geometry=HydroMomentumGeometry::build(*mesh,metricsFor(*mesh,rz,selection));inspectGeometry(*mesh,geometry);
            std::cout<<"GEOMETRY,"<<rz<<','<<amr<<','<<selection<<','<<geometry.layout().pressureCount()<<','<<geometry.layout().velocityCount()<<','<<geometry.segments().size()<<'\n';
            if(amr>1 || selection>1)continue;
            for(int kind=0;kind<4;++kind) {
                Case s;s.pureFirstJet=true;s.radial=rz;s.radius=rz?.34:.19;s.planar=kind!=0;s.planeOffset=kind==1?.43:kind==2?2:-1;
                const auto phase=HydroPhaseGeometry::build(*mesh,geometry.metrics(),initial(*mesh,s));
                if(amr==1 && selection==0 && kind==0)parallelMassContract(*mesh,geometry,phase,s);
                for(double density:{1.,1000.,1e6})inspectMass(*mesh,s,HydroMomentumMass::build(*mesh,geometry,phase,density,1));
            }
        }
        check(averageGap>1e-3,"fixture did not distinguish geometric dual mass from averaged cell density");std::cout<<"AVERAGE_GAP,"<<averageGap<<'\n';
        for(bool rz:{false,true})for(bool gas:{false,true}) {
            auto mesh=makeMesh(4);Case s;s.pureFirstJet=true;s.radial=rz;s.planar=true;s.planeOffset=gas?1-std::ldexp(1.,-40):std::ldexp(1.,-40);
            const auto geometry=HydroMomentumGeometry::build(*mesh,metricsFor(*mesh,rz));
            const auto phase=HydroPhaseGeometry::build(*mesh,geometry.metrics(),initial(*mesh,s));
            const auto mass=HydroMomentumMass::build(*mesh,geometry,phase,gas?1:1e6,gas?1e6:1);inspectMass(*mesh,s,mass);int tiny=0;
            for(int u=0;u<mass.layout().velocityCount();++u){const auto phases=mass.phaseVolumes()[u];const auto small=gas?phases.gas:phases.liquid;
                if(small>0 && small<1e-10*geometry.volumes()[u].volume){check(std::isfinite(small*1e6) && small*1e6>0 && mass.masses()[u]!=(gas?phases.liquid:phases.gas),"lost tiny phase material contribution");++tiny;}}
            check(tiny>0,"missing tiny dual phase");std::cout<<"TINY,"<<rz<<','<<gas<<','<<tiny<<'\n';
        }
        declaredContacts();overflow();
        for(bool rz:{false,true}){for(int amr:{0,1})convergence(rz,amr);guards(rz);}
        for(int rz=0;rz<2;++rz){check(amrCutDuals[rz]>0,"missing cut dual at actual coarse-fine boundary");std::cout<<"AMR_CUT,"<<rz<<','<<amrCutDuals[rz]<<'\n';}
        std::cout<<"SUMMARY,"<<geometryScenarios<<','<<massScenarios<<','<<overlapChecks<<','<<sharedPairs<<','<<segmentChecks<<','<<rejected<<','<<spaceRows<<','<<rebuilds<<','<<peakGeometry<<','<<peakPhase<<','<<peakMass<<','<<peakPartition<<'\n';
        std::cout<<"H09a HYDRO native momentum geometry and same-Q4 material mass PASS\n";
        coincidentMass();
        nativeTraceMass();
        nativeTraceRegionExtras();
        spatialAdvectorContract();
        advectorGuards();
        spatialAdvectorJetContract();
        cellJetGuards();
        regularPreparationContract();
        regularCubicSourceContract();
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
