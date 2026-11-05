// BEGIN L1 COST DIRECT INCLUDES
#include <chrono>
#ifdef _OPENMP
#include <omp.h>
#endif
// END L1 COST DIRECT INCLUDES
#include "Core/Mesh/NativeMeshIdentity.h"
#include "Core/Geometry/Detail/FiniteVolumeClosure.h"
#include "Modules/Interface/GeometricVofTransport.h" // Public API stands alone.
#include "Modules/Interface/Detail/VofLocalFilteredDyadic.h"
#include "Modules/Electric/ElectricVofChargeTransport.h"
#include "Modules/Interface/VofSmoothGeometry.h"
#include "Modules/Interface/GeometricVofInitialization.h"
#include "Modules/Interface/GeometricVofPhaseBudget.h"
#include "Modules/Hydrodynamics/HydroVofAdvector.h"
#include "Modules/Hydrodynamics/HydroVofStep.h"
#include "Modules/Hydrodynamics/HydroVofMeshTransfer.h"
#include "Modules/Hydrodynamics/HydroMomentumGeometry.h"
#include "Modules/Hydrodynamics/HydroDualMaterial.h"
#include "Modules/Hydrodynamics/HydroKepConvection.h"
#include "Modules/Hydrodynamics/Detail/HydroFiniteContraction.h"
#include "Modules/Interface/GeometricVofMeshTransfer.h"
#include "Core/Field/CellField.h"
#include "Core/Field/VEdgeField.h"
#include "Core/Field/HEdgeField.h"
#include "Core/Mesh/CartesianMeasure.h"
#include "Core/Mesh/FaceSegmentGraph.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Geometry/PlanarFaceTimeIntegral.h"
#include "Core/Geometry/Detail/ExactRationalRound.h"
#include <cerrno>
#include <cfenv>
#if defined(_MSC_VER) && defined(_M_X64)
#include <float.h>
#include <xmmintrin.h>
#endif
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <exception>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <typeinfo>
#include <type_traits>
#include <vector>

namespace
{
using namespace LAMP::Interface;
using namespace LAMP::Hydrodynamics;
using Publication=GeometricVofPublication;
using Transport=GeometricVofTransport;
using Slab=GeometricVofAdvectorSlab;
using Clock=GeometricVofTimeLayer;
using Position=Clock::Position;
using Point=LAMP::Geometry::PlanarCellCut::Point;
constexpr long double conservationTolerance=1e-12L;
constexpr long double rawTolerance=2e-9L; // Frozen independent space/time quadrature gate.
constexpr long double arithmetic=2048*std::numeric_limits<long double>::epsilon();
const long double twoPi=CartesianMeasure::twoPi;
unsigned assertions=0,scenarios=0,guards=0,rawChecks=0,coarseChecks=0,counterflowFaces=0,limitedFaces=0;
long double peakConservation=0,peakRaw=0,peakInventory=0,peakGeometry=0,peakCorrection=0;
unsigned rationalRoundScenarios=0,rationalRoundRejections=0;
unsigned capturedBudgetScenarios=0,capturedBudgetGrossChecks=0,capturedBudgetPureChecks=0;
unsigned microScenarios=0;
long double peakMicroRelative=0;
unsigned faceTimeScenarios=0,faceTimeAttempts=0,faceTimePositiveChecks=0,faceTimeZeroChecks=0;
unsigned maximumFaceTimeRootRefinements=0;
long double peakFaceTimeRelative=0,peakFaceTimeUncertainty=0;
std::string stage="setup";
static_assert(!std::is_default_constructible_v<Transport>);
static_assert(std::is_base_of_v<std::invalid_argument,GeometricVofCourantExceeded>);
static_assert(!std::is_default_constructible_v<Slab>);
static_assert(std::is_same_v<decltype(std::declval<const Transport&>().endpoint()),const Publication&>);
void check(bool ok,const std::string& reason)
{ ++assertions;if(!ok)throw std::runtime_error(stage+": "+reason); }
void near(long double actual,long double expected,long double scale,const std::string& name,long double tolerance=arithmetic)
{
    check(std::isfinite(actual)&&std::isfinite(expected)&&std::isfinite(scale)&&scale>0,name+" invalid comparison");
    const auto error=std::abs(actual-expected)/scale;
    std::ostringstream message;message<<std::setprecision(17)<<name<<" actual="<<actual<<" expected="<<expected<<" residual="<<error;
    check(error<=tolerance,message.str());
}
template<class F> void rejects(F&& f,const std::string& name)
{
    bool thrown=false;try{f();}catch(const std::exception&){thrown=true;}
    check(thrown,"missing explicit rejection: "+name);++guards;
}
bool sameClock(Clock a,Clock b)
{ return a.step==b.step&&a.begin==b.begin&&a.end==b.end&&a.position==b.position; }
struct Box {long double x0,x1,y0,y1;};
Box bounds(const QuadtreeMesh& mesh,int ci)
{
    const auto& cell=mesh.cells(mesh.cellIndexToLabel[ci]);
    const auto a=mesh.verts(cell.cellToVert[0]).coord(),b=mesh.verts(cell.cellToVert[2]).coord();
    return {a.x,b.x,a.y,b.y};
}
long double capacity(Box b,bool rz)
{ return (b.x1-b.x0)*(b.y1-b.y0)*(rz?twoPi*(b.x0+b.x1)/2:1); }
NativeLeafDomain domain(const QuadtreeMesh& mesh,bool rz)
{ FaceSegmentGraph graph;graph.build(mesh);return NativeLeafDomain::build(mesh,graph,rz); }
// BEGIN V4a1b INDEPENDENT AP MATERIAL ORACLE (original B628 Wide/clip/Green/M math)
NativeLeafDomain domainFor(const QuadtreeMesh& mesh,bool rz){return domain(mesh,rz);}
HydroPhysicalDivergence metricsFor(const QuadtreeMesh& mesh,bool rz)
{return HydroPhysicalDivergence::build(mesh,HydroIntegralDivergence::build(mesh,HydroDofLayout::build(mesh,domainFor(mesh,rz))));}
// V4a: same retained VOF fan -> actual native dual material mass. The
// independent reference keeps twofold clipping intersections and integrates
// both phases, without calling the production Core clip/measure routines.
#if defined(_MSC_VER)
#pragma float_control(precise,on,push)
#endif
struct VofMaterialWide { long double high=0,low=0; };
VofMaterialWide vofMaterialSum(long double a,long double b)
{
    const auto s=a+b,p=s-a;return {s,(a-(s-p))+(b-p)};
}
VofMaterialWide vofMaterialAdd(VofMaterialWide a,VofMaterialWide b)
{
    const auto h=vofMaterialSum(a.high,b.high),l=vofMaterialSum(a.low,b.low);
    const auto m=vofMaterialSum(h.low,l.high),n=vofMaterialSum(h.high,m.high);
    return vofMaterialSum(n.high,n.low+m.low+l.low);
}
VofMaterialWide vofMaterialNeg(VofMaterialWide a){return {-a.high,-a.low};}
VofMaterialWide vofMaterialSub(VofMaterialWide a,VofMaterialWide b)
{return vofMaterialAdd(a,vofMaterialNeg(b));}
VofMaterialWide vofMaterialMul(VofMaterialWide a,VofMaterialWide b)
{
    const auto p=a.high*b.high;auto r=VofMaterialWide{p,std::fma(a.high,b.high,-p)};
    r=vofMaterialAdd(r,{a.high*b.low,0});r=vofMaterialAdd(r,{a.low*b.high,0});
    return vofMaterialAdd(r,{a.low*b.low,0});
}
VofMaterialWide vofMaterialDiv(VofMaterialWide a,VofMaterialWide b)
{
    check(b.high!=0,"VOF reference division by zero");const auto q=a.high/b.high;
    const auto r=vofMaterialSub(a,vofMaterialMul({q,0},b));
    const auto q1=(r.high+r.low)/b.high;
    const auto r1=vofMaterialSub(r,vofMaterialMul({q1,0},b));
    return vofMaterialAdd(vofMaterialAdd({q,0},{q1,0}),{(r1.high+r1.low)/b.high,0});
}
long double vofMaterialValue(VofMaterialWide a){return a.high+a.low;}
int vofMaterialSign(VofMaterialWide a)
{return a.high>0?1:a.high<0?-1:a.low>0?1:a.low<0?-1:0;}
using VofMaterialPoint=std::array<VofMaterialWide,2>;
using VofMaterialPolygon=std::vector<VofMaterialPoint>;
using VofMaterialBounds=std::array<VofMaterialWide,4>; // x0,x1,y0,y1
using VofPublication=LAMP::Interface::GeometricVofPublication;
using VofClock=LAMP::Interface::GeometricVofTimeLayer;
using VofTransport=LAMP::Interface::GeometricVofTransport;
unsigned vofMaterialLayers=0,vofMaterialComparisons=0,vofMaterialGuards=0,vofMaterialSlabs=0;
std::array<unsigned,2> vofMaterialAmrCuts{};
long double vofMaterialPeakPhase=0,vofMaterialPeakMass=0,vofMaterialPeakPartition=0,
    vofMaterialPeakInventory=0,vofMaterialPeakGasInventory=0,vofMaterialProjectionGap=0;
bool vofMaterialDifferentCandidates=false;
unsigned vofMaterialFaceChecks=0,vofMaterialCoincidentFaces=0,vofMaterialAxisFaces=0;
void vofMaterialNear(long double actual,VofMaterialWide expected,const std::string& what)
{
    ++vofMaterialComparisons;const auto ref=vofMaterialValue(expected);
    check(std::isfinite(actual)&&std::isfinite(ref)&&vofMaterialSign(expected)>=0,"VOF invalid "+what);
    if(vofMaterialSign(expected)==0){check(actual==0,"VOF true zero changed: "+what);return;}
    const auto e=std::abs(actual-ref)/ref;vofMaterialPeakPhase=std::max(vofMaterialPeakPhase,e);
    std::ostringstream message;message<<std::setprecision(17)<<"VOF "<<what<<" actual="<<actual<<" expected="<<ref<<" relative="<<e;
    check(actual>0&&e<=1e-12L,message.str());
}
VofMaterialBounds vofMaterialLocalBounds(Vector2D lo,Vector2D hi,LAMP::Geometry::PlanarCellCut::Point origin)
{return {vofMaterialSum(lo.x,-origin.x),vofMaterialSum(hi.x,-origin.x),vofMaterialSum(lo.y,-origin.y),vofMaterialSum(hi.y,-origin.y)};}
VofMaterialPolygon vofMaterialClip(VofMaterialPolygon polygon,const VofMaterialBounds& b)
{
    for(int side=0;side<4&&!polygon.empty();++side){
        const int axis=side/2;const auto bound=b[side];const int direction=side%2?-1:1;
        const auto inside=[&](const VofMaterialPoint& p){return direction*vofMaterialSign(vofMaterialSub(p[axis],bound))>=0;};
        VofMaterialPolygon result;auto a=polygon.back();bool ain=inside(a);
        for(const auto& c:polygon){const bool cin=inside(c);
            if(ain!=cin){const auto t=vofMaterialDiv(vofMaterialSub(bound,a[axis]),vofMaterialSub(c[axis],a[axis]));
                VofMaterialPoint cut;cut[axis]=bound;cut[1-axis]=vofMaterialAdd(a[1-axis],vofMaterialMul(t,vofMaterialSub(c[1-axis],a[1-axis])));result.push_back(cut);}
            if(cin)result.push_back(c);a=c;ain=cin;
        }polygon=std::move(result);
    }return polygon;
}
VofMaterialWide vofMaterialGreen(const VofMaterialPolygon& p,long double originX,bool rz)
{
    if(p.size()<3)return {};VofMaterialWide result;const auto& a=p.front();
    for(std::size_t i=1;i+1<p.size();++i){
        const auto bx=vofMaterialSub(p[i][0],a[0]),by=vofMaterialSub(p[i][1],a[1]);
        const auto cx=vofMaterialSub(p[i+1][0],a[0]),cy=vofMaterialSub(p[i+1][1],a[1]);
        const auto det=vofMaterialSub(vofMaterialMul(bx,cy),vofMaterialMul(by,cx));
        check(vofMaterialSign(det)>=0,"VOF reference polygon not convex CCW");
        auto measure=vofMaterialDiv(det,{2,0});
        if(rz){const auto centroidRadius=vofMaterialAdd(vofMaterialAdd({originX,0},a[0]),vofMaterialDiv(vofMaterialAdd(bx,cx),{3,0}));
            measure=vofMaterialMul(vofMaterialMul(measure,centroidRadius),{CartesianMeasure::twoPi,0});}
        result=vofMaterialAdd(result,measure);
    }return result;
}
VofMaterialWide vofMaterialReference(const VofPublication& p,int native,bool liquid,const VofMaterialBounds* b=nullptr)
{
    VofMaterialWide sum;for(const auto& polygon:p.anchoredPhasePolygons(native,liquid)){
        VofMaterialPolygon wide;for(const auto& q:polygon)wide.push_back({{vofMaterialSum(q.anchor.x,q.offset.x),vofMaterialSum(q.anchor.y,q.offset.y)}});
        if(b)wide=vofMaterialClip(std::move(wide),*b);
        sum=vofMaterialAdd(sum,vofMaterialGreen(wide,p.frameOrigin().x,p.domain().axisymmetric()));
    }return sum;
}
long double vofMaterialCapacity(Vector2D lo,Vector2D hi,bool rz)
{return (static_cast<long double>(hi.x)-lo.x)*(static_cast<long double>(hi.y)-lo.y)*(rz?CartesianMeasure::twoPi*(static_cast<long double>(lo.x)+hi.x)/2:1.L);}
std::unique_ptr<QuadtreeMesh> vofMaterialMesh(int n)
{
    DictionaryMesh d;d.Nx=d.Ny=n;d.maxCellNum=4096;d.maxInitialIterations=0;d.maxLevel_local={2};
    auto mesh=std::make_unique<QuadtreeMesh>(d);DictionaryMesh::RefineBox b;b.level=1;b.min={.5,0};b.max={1,1};mesh->applyRefineBoxes({b});
    check(mesh->cellIndexNum()==5*n*n/2,"VOF material fixture leaf count changed");return mesh;
}
long double vofMaterialLine(Vector2D lo,Vector2D hi,bool rz,long double c)
{
    if(c>=static_cast<long double>(hi.x)+hi.y)return vofMaterialCapacity(lo,hi,rz);
    if(c<=static_cast<long double>(lo.x)+lo.y)return 0;
    const auto full=std::clamp(c-hi.y,static_cast<long double>(lo.x),static_cast<long double>(hi.x));
    const auto end=std::clamp(c-lo.y,static_cast<long double>(lo.x),static_cast<long double>(hi.x));
    auto a=(full-lo.x)*(hi.y-lo.y),r=(full*full-lo.x*lo.x)*(hi.y-lo.y)/2;
    if(end>full){const auto h=c-lo.y,s=(end*end-full*full)/2,t=(end*end*end-full*full*full)/3;
        a+=h*(end-full)-s;r+=h*s-t;}
    return rz?CartesianMeasure::twoPi*r:a;
}
VofPublication vofMaterialPublish(const QuadtreeMesh& mesh,bool rz,long double c,int kind=-1)
{
    std::vector<double> alpha(mesh.cellIndexNum());
    for(int ci=0;ci<mesh.cellIndexNum();++ci){const auto& cell=mesh.cells(mesh.cellIndexToLabel[ci]);
        const auto lo=mesh.verts(cell.cellToVert[0]).coord(),hi=mesh.verts(cell.cellToVert[2]).coord();
        const auto cap=vofMaterialCapacity(lo,hi,rz);long double liquid=0;
        if(kind==0)liquid=0;else if(kind==1)liquid=cap;else if(kind==2){auto cut=hi;cut.x=std::min(cut.x,.5);liquid=cut.x<=lo.x?0:vofMaterialCapacity(lo,cut,rz);}
        else liquid=vofMaterialLine(lo,hi,rz,c);
        check(liquid>=0&&liquid<=cap,"VOF material analytic alpha invalid");alpha[ci]=liquid==0?0:liquid==cap?1:static_cast<double>(liquid/cap);
    }
    return VofPublication::prepare(mesh,domainFor(mesh,rz),{7,.125,.25,VofClock::Position::End},alpha);
}
HydroPhysicalDivergence vofMaterialMetrics(const QuadtreeMesh& mesh,const VofPublication& p)
{return HydroPhysicalDivergence::build(mesh,HydroIntegralDivergence::build(mesh,HydroDofLayout::build(mesh,p.domain())));}
bool vofMaterialSameClock(VofClock a,VofClock b)
{return a.step==b.step&&a.begin==b.begin&&a.end==b.end&&a.position==b.position;}

std::array<VofMaterialWide,2> vofMaterialTraceReference(const VofPublication& p,const LAMP::Interface::GeometricVofFace& face,bool minus)
{
    std::array<VofMaterialWide,2> result;
    for(unsigned i=0;i<face.rangeCount;++i){const auto& range=face.ranges[i];
        const auto a=p.anchoredVertices()[range.vertices[0]],b=p.anchoredVertices()[range.vertices[1]];
        const auto ax=vofMaterialSum(a.anchor.x,a.offset.x),ay=vofMaterialSum(a.anchor.y,a.offset.y);
        const auto dx=vofMaterialSub(vofMaterialSum(b.anchor.x,b.offset.x),ax),dy=vofMaterialSub(vofMaterialSum(b.anchor.y,b.offset.y),ay);
        const auto length=vofMaterialSign(dx)==0?dy:dx;
        check(vofMaterialSign(length)>0,"VOF retained native trace orientation invalid");auto physical=length;
        if(p.domain().axisymmetric()){
            const auto radius=vofMaterialAdd(vofMaterialAdd({p.frameOrigin().x,0},ax),vofMaterialDiv(dx,{2,0}));
            physical=vofMaterialMul(vofMaterialMul(physical,radius),{CartesianMeasure::twoPi,0});
        }
        const bool liquid=minus?range.minusLiquid:range.plusLiquid;
        result[liquid?0:1]=vofMaterialAdd(result[liquid?0:1],physical);
    }return result;
}
std::array<VofMaterialWide,2> vofMaterialAnalyticTrace(const QuadtreeMesh& mesh,const NativeLeafFace& f,int native,int kind,bool rz)
{
    std::array<VofMaterialWide,2> result;if(native<0)return result;
    const auto& edge=f.horizontal?mesh.hEdges(f.edgeLabel):mesh.vEdges(f.edgeLabel);
    auto a=mesh.verts(edge.edgeToVert[0]).coord(),b=mesh.verts(edge.edgeToVert[1]).coord();if(f.horizontal?a.x>b.x:a.y>b.y)std::swap(a,b);
    if(kind<2||a.x==b.x){bool liquid=kind==1;
        if(kind==2){const auto& cell=mesh.cells(mesh.cellIndexToLabel[native]);liquid=cell.centerCoord().x<.5;}
        result[liquid?0:1]={f.area,0};return result;
    }
    const auto measure=[&](long double lo,long double hi){
        if(!(hi>lo))return VofMaterialWide{};const auto width=vofMaterialSum(hi,-lo);
        return rz?vofMaterialMul(vofMaterialMul(width,vofMaterialDiv(vofMaterialSum(lo,hi),{2,0})),{CartesianMeasure::twoPi,0}):width;
    };
    result[0]=measure(a.x,std::min(static_cast<long double>(b.x),.5L));
    result[1]=measure(std::max(static_cast<long double>(a.x),.5L),b.x);return result;
}
void vofMaterialFaceContract(const QuadtreeMesh& mesh,const HydroPhaseGeometry& phase,const VofPublication& p,int analytic)
{
    const auto& layout=phase.layout();for(int u=0;u<layout.velocityCount();++u){
        const int row=layout.nativeFace(u);const auto& source=p.faces()[row];const auto& trace=phase.faceTraces()[u];const auto native=layout.domain().face(mesh,row);
        check(trace.negativeCell==source.negativeCell&&trace.positiveCell==source.positiveCell&&trace.coincidentHalfFaces==0,"VOF whole-face metadata fabricated/changed");
        const auto minus=vofMaterialTraceReference(p,source,true),plus=vofMaterialTraceReference(p,source,false);
        std::array<VofMaterialWide,2> spatial;
        if(trace.minusAvailable()){
            vofMaterialNear(trace.minus.liquid,minus[0],"minus face liquid row="+std::to_string(row));
            vofMaterialNear(trace.minus.gas,minus[1],"minus face gas row="+std::to_string(row));
            vofMaterialNear(trace.minus.liquid,{source.liquidMinus,0},"actual source minus liquid");}
        if(trace.plusAvailable()){
            vofMaterialNear(trace.plus.liquid,plus[0],"plus face liquid row="+std::to_string(row));
            vofMaterialNear(trace.plus.gas,plus[1],"plus face gas row="+std::to_string(row));
            vofMaterialNear(trace.plus.liquid,{source.liquidPlus,0},"actual source plus liquid");}
        const auto expectedKind=!source.coincident?HydroPhaseFaceTraceKind::Ordinary:
            trace.minusAvailable()&&trace.plusAvailable()?HydroPhaseFaceTraceKind::InteriorCoincident:HydroPhaseFaceTraceKind::BoundaryCoincident;
        check(trace.kind==expectedKind,"VOF actual whole-face kind changed");vofMaterialCoincidentFaces+=source.coincident;
        for(int k=0;k<2;++k){
            spatial[k]=trace.minusAvailable()&&trace.plusAvailable()?vofMaterialDiv(vofMaterialAdd(minus[k],plus[k]),{2,0}):trace.minusAvailable()?minus[k]:plus[k];
            vofMaterialNear(k==0?phase.faceMeasures()[u].liquid:phase.faceMeasures()[u].gas,spatial[k],"spatial same-face phase");}
        if(analytic>=0){const auto am=vofMaterialAnalyticTrace(mesh,native,source.negativeCell,analytic,p.domain().axisymmetric());
            const auto ap=vofMaterialAnalyticTrace(mesh,native,source.positiveCell,analytic,p.domain().axisymmetric());
            if(trace.minusAvailable()){vofMaterialNear(trace.minus.liquid,am[0],"analytic minus face liquid");vofMaterialNear(trace.minus.gas,am[1],"analytic minus face gas");}
            if(trace.plusAvailable()){vofMaterialNear(trace.plus.liquid,ap[0],"analytic plus face liquid");vofMaterialNear(trace.plus.gas,ap[1],"analytic plus face gas");}
        }
        if(p.domain().axisymmetric()&&native.area==0){++vofMaterialAxisFaces;check(phase.faceMeasures()[u].liquid==0&&phase.faceMeasures()[u].gas==0
            &&trace.minus.liquid==0&&trace.minus.gas==0&&trace.plus.liquid==0&&trace.plus.gas==0,"VOF radial axis physical phase not zero");}
        ++vofMaterialFaceChecks;
    }
}

HydroMomentumMass vofMaterialInspect(const QuadtreeMesh& mesh,const HydroMomentumGeometry& geometry,
    const VofPublication& publication,int analytic=-1)
{
    const auto phase=HydroPhaseGeometry::build(mesh,geometry.metrics(),publication);
    const auto mass=HydroMomentumMass::build(mesh,geometry,phase,100,1);const auto& layout=geometry.layout();const bool rz=publication.domain().axisymmetric();
    check(phase.model()==HydroPhaseGeometry::Model::GeometricVof&&mass.model()==phase.model(),"VOF material model lost");
    check(phase.vofIdentity()==publication.identity()&&mass.vofIdentity()==publication.identity(),"VOF material typed identity changed");
    check(&phase.vofPublication().alpha()==&publication.alpha()&&&mass.phaseGeometry().vofPublication().alpha()==&publication.alpha(),"VOF material detached actual alpha owner");
    check(vofMaterialSameClock(phase.timeLayer(),publication.timeLayer())&&vofMaterialSameClock(mass.timeLayer(),publication.timeLayer()),"VOF material full clock changed");
    phase.validateVofFor(mesh,publication.identity(),publication.timeLayer());mass.validateVofFor(mesh,publication.identity(),publication.timeLayer());
    check(mass.geometry().sameIdentity(geometry)&&mass.layout().sameIdentity(layout)&&phase.layout().sameIdentity(layout)
        &&&mass.geometry().metrics().divergence().matrix()==&phase.metrics().divergence().matrix(),"VOF material native geometry/D changed");
    check(mass.representation()==HydroMomentumMass::DensityRepresentation::SharpGeometry,"VOF actual fan mass projected");
    check(mass.liquidDensity()==100&&mass.gasDensity()==1,"VOF actual density ratio changed");
    vofMaterialFaceContract(mesh,phase,publication,analytic);
    std::array<std::vector<std::array<VofMaterialWide,2>>,2> partition;for(auto& row:partition)row.resize(layout.pressureCount());
    std::array<VofMaterialWide,2> totalMass;std::array<std::array<VofMaterialWide,2>,2> totalPhase{};
    for(int p=0;p<layout.pressureCount();++p){const auto native=layout.nativeCell(p);const auto& source=publication.cells()[native];
        check(phase.liquidFractions()[p]==publication.alpha()[native],"VOF material copied/rederived alpha changed bits");
        vofMaterialNear(phase.cellVolumes()[p].liquid,vofMaterialReference(publication,native,true),"native cell liquid");
        vofMaterialNear(phase.cellVolumes()[p].gas,vofMaterialReference(publication,native,false),"native cell gas");
        const auto minor=std::min(source.inventory,source.capacity-source.inventory);
        if(minor>0){const auto e=std::abs(phase.cellVolumes()[p].liquid-source.inventory)/minor;
            check(e<=1e-12L,"VOF measured liquid/inventory original minor gate");vofMaterialPeakInventory=std::max(vofMaterialPeakInventory,e);
            vofMaterialPeakGasInventory=std::max(vofMaterialPeakGasInventory,std::abs(phase.cellVolumes()[p].gas-(source.capacity-source.inventory))/minor);}
    }
    for(int u=0;u<layout.velocityCount();++u){const auto& volume=geometry.volumes()[u];const int axis=layout.normalAxis(u)==HydroDofLayout::NormalAxis::First?0:1;
        std::array<VofMaterialWide,2> dual;const auto wholeBounds=vofMaterialLocalBounds(volume.lower,volume.upper,publication.frameOrigin());
        long double projected=0;
        check(mass.overlapPhaseVolumes()[u].size()==volume.overlaps.size(),"VOF overlap ordering changed");
        for(std::size_t k=0;k<volume.overlaps.size();++k){const auto& overlap=volume.overlaps[k];const auto native=layout.nativeCell(overlap.pressure);
            const auto b=vofMaterialLocalBounds(overlap.lower,overlap.upper,publication.frameOrigin());
            const auto l=vofMaterialReference(publication,native,true,&b),g=vofMaterialReference(publication,native,false,&b);
            const auto actual=mass.overlapPhaseVolumes()[u][k];vofMaterialNear(actual.liquid,l,"dual overlap liquid");vofMaterialNear(actual.gas,g,"dual overlap gas");
            check(std::abs(actual.liquid+actual.gas-overlap.volume)<=1e-12L*overlap.volume,"VOF actual overlap closure");
            partition[axis][overlap.pressure][0]=vofMaterialAdd(partition[axis][overlap.pressure][0],l);
            partition[axis][overlap.pressure][1]=vofMaterialAdd(partition[axis][overlap.pressure][1],g);
            dual[0]=vofMaterialAdd(dual[0],vofMaterialReference(publication,native,true,&wholeBounds));
            dual[1]=vofMaterialAdd(dual[1],vofMaterialReference(publication,native,false,&wholeBounds));
            projected+=overlap.cellFraction*(100*phase.cellVolumes()[overlap.pressure].liquid+phase.cellVolumes()[overlap.pressure].gas);
        }
        vofMaterialNear(mass.phaseVolumes()[u].liquid,dual[0],"whole dual liquid");vofMaterialNear(mass.phaseVolumes()[u].gas,dual[1],"whole dual gas");
        const auto expected=vofMaterialAdd(vofMaterialMul({100,0},dual[0]),dual[1]);const auto expectedValue=vofMaterialValue(expected);
        const auto e=std::abs(mass.masses()[u]-expectedValue)/expectedValue;vofMaterialPeakMass=std::max(vofMaterialPeakMass,e);
        check(mass.masses()[u]>0&&std::isfinite(mass.masses()[u])&&e<=1e-12L,"VOF actual rho100:1 native mass");
        vofMaterialProjectionGap=std::max(vofMaterialProjectionGap,std::abs(projected-expectedValue)/(101*volume.volume));
        check(std::abs(mass.phaseVolumes()[u].liquid+mass.phaseVolumes()[u].gas-volume.volume)<=1e-12L*volume.volume,"VOF actual dual closure");
        const auto nativeFace=layout.domain().face(mesh,layout.nativeFace(u));
        if(nativeFace.negativeCell>=0&&nativeFace.positiveCell>=0&&vofMaterialSign(dual[0])>0&&vofMaterialSign(dual[1])>0){
            const auto& a=mesh.cells(mesh.cellIndexToLabel[nativeFace.negativeCell]);const auto& b=mesh.cells(mesh.cellIndexToLabel[nativeFace.positiveCell]);
            if(a.size().x!=b.size().x)++vofMaterialAmrCuts[rz?1:0];}
        totalMass[axis]=vofMaterialAdd(totalMass[axis],expected);for(int k=0;k<2;++k)totalPhase[axis][k]=vofMaterialAdd(totalPhase[axis][k],dual[k]);
        if(analytic>=0){auto leftHi=volume.upper,rightLo=volume.lower;long double l=0,g=0;
            if(analytic==0)g=vofMaterialCapacity(volume.lower,volume.upper,rz);
            else if(analytic==1)l=vofMaterialCapacity(volume.lower,volume.upper,rz);
            else {leftHi.x=std::min(leftHi.x,.5);rightLo.x=std::max(rightLo.x,.5);
                if(leftHi.x>volume.lower.x)l=vofMaterialCapacity(volume.lower,leftHi,rz);
                if(rightLo.x<volume.upper.x)g=vofMaterialCapacity(rightLo,volume.upper,rz);}
            vofMaterialNear(mass.phaseVolumes()[u].liquid,{l,0},"analytic pure/vertical dual liquid");
            vofMaterialNear(mass.phaseVolumes()[u].gas,{g,0},"analytic pure/vertical dual gas");}
    }
    for(int axis=0;axis<2;++axis){
        for(int p=0;p<layout.pressureCount();++p)for(int k=0;k<2;++k){const auto ref=vofMaterialReference(publication,layout.nativeCell(p),k==0),sum=partition[axis][p][k];
            vofMaterialNear(vofMaterialValue(sum),ref,"component dual/cell phase partition");
            if(vofMaterialSign(ref)>0)vofMaterialPeakPartition=std::max(vofMaterialPeakPartition,std::abs(vofMaterialValue(vofMaterialSub(sum,ref)))/vofMaterialValue(ref));}
        vofMaterialNear(mass.report().totalPhaseVolumes[axis].liquid,totalPhase[axis][0],"reported liquid total");
        vofMaterialNear(mass.report().totalPhaseVolumes[axis].gas,totalPhase[axis][1],"reported gas total");
        vofMaterialNear(mass.report().totalMass[axis],totalMass[axis],"reported native mass total");
    }
    ++vofMaterialLayers;std::cout<<"V4A_MASS_LAYER,"<<vofMaterialLayers<<','<<rz<<','<<layout.pressureCount()<<','<<layout.velocityCount()<<','
        <<publication.timeLayer().step<<','<<static_cast<int>(publication.timeLayer().position)<<','<<mass.report().cutOverlaps<<",PASS\n";return mass;
}
long double vofMaterialPsi(Vector2D p,int variant)
{
    const long double x=p.x,y=p.y;const auto base=x*x*(1-x)*(1-x)*y*y*(1-y)*(1-y);
    return base*(variant==0?1:variant==1?x+.4L*y+.2L:-2.L);
}
std::vector<double> vofMaterialNormals(const QuadtreeMesh& mesh,const VofPublication& p,int variant)
{
    std::vector<double> result(p.domain().nativeFaceCount());
    for(int row=0;row<p.domain().nativeFaceCount();++row){const auto f=p.domain().face(mesh,row);const auto& edge=f.horizontal?mesh.hEdges(f.edgeLabel):mesh.vEdges(f.edgeLabel);
        auto a=mesh.verts(edge.edgeToVert[0]).coord(),b=mesh.verts(edge.edgeToVert[1]).coord();if(f.horizontal?a.x>b.x:a.y>b.y)std::swap(a,b);
        const auto q=(f.horizontal?-1:1)*(vofMaterialPsi(b,variant)-vofMaterialPsi(a,variant));
        if(f.outwardSign!=0||f.area==0){check(q==0,"VOF material streamfunction boundary leakage");result[row]=0;}else result[row]=static_cast<double>(q/f.area);
    }return result;
}
long double vofMaterialOutflow(const QuadtreeMesh& mesh,const VofPublication& p,const std::vector<double>& a,const std::vector<double>& b)
{
    std::vector<long double> sums(p.cells().size());const auto positive=[](long double x,long double y){
        if(x<=0&&y<=0)return 0.L;if(x>=0&&y>=0)return (x+y)/2;const auto z=std::max(x,y);return z*z/(2*std::abs(y-x));};
    for(int row=0;row<p.domain().nativeFaceCount();++row){const auto f=p.domain().face(mesh,row);
        if(f.negativeCell>=0)sums[f.negativeCell]+=f.area*positive(a[row],b[row]);
        if(f.positiveCell>=0)sums[f.positiveCell]+=f.area*positive(-a[row],-b[row]);}
    long double maximum=0;for(const auto& c:p.cells())maximum=std::max(maximum,sums[c.nativeCell]/c.capacity);check(maximum>0,"VOF material slab zero motion");return maximum;
}
LAMP::Interface::GeometricVofAdvectorSlab vofMaterialSlab(const QuadtreeMesh& mesh,const VofPublication& p,const HydroDofLayout& layout,
    VofClock next,const std::vector<double>& a,const std::vector<double>& b)
{
    std::vector<double> av(layout.velocityCount()),bv(av.size());for(int u=0;u<layout.velocityCount();++u){av[u]=a[layout.nativeFace(u)];bv[u]=b[layout.nativeFace(u)];}
    const HydroNativeValues begin(layout,av,std::vector<double>(layout.pressureCount(),.125)),end(layout,bv,std::vector<double>(layout.pressureCount(),-.25));
    const auto result=prepareHydroVofAdvectorSlab(mesh,p,next,begin,end);
    check(result.beginNormal()==a&&result.endNormal()==b&&begin.normalVelocity()==av&&end.normalVelocity()==bv,"VOF material native slab input changed");return result;
}
template<class F>void vofMaterialReject(F action,const std::string& what)
{
    bool caught=false;try{action();}catch(const std::invalid_argument&){caught=true;}
    check(caught,"missing explicit VOF material rejection: "+what);++vofMaterialGuards;
}
void vofMaterialGuardContract(QuadtreeMesh& mesh,const HydroMomentumGeometry& geometry,const VofPublication& p,const HydroMomentumMass& mass,bool mutate=false)
{
    const auto saved=p.alpha(),savedMass=mass.masses();const auto identity=p.identity();const auto clock=p.timeLayer();const auto& phase=mass.phaseGeometry();
    const auto foreign=VofPublication::prepare(mesh,p.domain(),clock,p.alpha());
    vofMaterialReject([&]{phase.validateVofFor(mesh,foreign.identity(),clock);},"same clock foreign owner");
    vofMaterialReject([&]{mass.validateVofFor(mesh,foreign.identity(),clock);},"same clock foreign mass owner");
    auto wrong=clock;wrong.begin-=.0625;vofMaterialReject([&]{mass.validateVofFor(mesh,identity,wrong);},"full clock begin");
    wrong=clock;wrong.position=VofClock::Position::Midpoint;vofMaterialReject([&]{phase.validateVofFor(mesh,identity,wrong);},"full clock position");
    vofMaterialReject([&]{(void)phase.publication();},"VOF cannot pretend FirstJet publication");
    vofMaterialReject([&]{(void)phase.geometryIdentity();},"VOF cannot pretend FirstJet identity");
    vofMaterialReject([&]{(void)mass.geometryIdentity();},"VOF mass cannot pretend FirstJet identity");
    vofMaterialReject([&]{HydroMomentumMass::buildCellAverageProjection(mesh,geometry,phase,100,1);},"CompatibleFan VOF cannot use independent inventory projection");
    for(double rho:{0.,-1.,std::numeric_limits<double>::infinity(),std::numeric_limits<double>::quiet_NaN()})
        vofMaterialReject([&]{HydroMomentumMass::build(mesh,geometry,phase,rho,1);},"invalid density");
    vofMaterialReject([&]{HydroPhaseGeometry::build(mesh,metricsFor(mesh,!p.domain().axisymmetric()),p);},"XY/RZ measure");
    vofMaterialReject([&]{HydroPhaseGeometry::build(mesh,metricsFor(mesh,p.domain().axisymmetric()),p);},"foreign equivalent native mesh domain");
    const auto native=geometry.layout().nativeCell(0);const auto& cell=mesh.cells(mesh.cellIndexToLabel[native]);
    const auto lo=mesh.verts(cell.cellToVert[0]).coord(),hi=mesh.verts(cell.cellToVert[2]).coord();
    HydroVofBox bad{native,{lo,hi}};bad.rectangle.lower.x-=.125;
    vofMaterialReject([&]{phase.measureVofBoxes(mesh,{bad});},"outside selected native leaf box");
    bad={-1,{lo,hi}};vofMaterialReject([&]{phase.measureVofBoxes(mesh,{bad});},"invalid native leaf box");
    bad={native,{lo,hi}};bad.rectangle.lower.x=hi.x;vofMaterialReject([&]{phase.measureVofBoxes(mesh,{bad});},"zero width native leaf box");
    const auto otherD=HydroPhysicalDivergence::build(mesh,HydroIntegralDivergence::build(mesh,geometry.layout()));
    const auto otherPhase=HydroPhaseGeometry::build(mesh,otherD,p);
    vofMaterialReject([&]{HydroMomentumMass::build(mesh,geometry,otherPhase,100,1);},"different retained D_I");
    if(mutate){DictionaryMesh::RefineBox b;b.level=1;b.min={0,.5};b.max={.5,1};mesh.applyRefineBoxes({b});
        vofMaterialReject([&]{mass.validateVofFor(mesh,identity,clock);},"mesh revision");}
    check(p.alpha()==saved&&p.identity()==identity&&vofMaterialSameClock(p.timeLayer(),clock)&&mass.masses()==savedMass,"VOF material rejection polluted old source/mass");
}
#if defined(_MSC_VER)
#pragma float_control(pop)
#endif
using V4bAP=LAMP::Geometry::PlanarCellCut::AnchoredPoint;
VofMaterialPoint v4bWidePoint(V4bAP q)
{return {{vofMaterialSum(q.anchor.x,q.offset.x),vofMaterialSum(q.anchor.y,q.offset.y)}};}
long double measure(const std::vector<V4bAP>& polygon,Point origin,bool rz)
{VofMaterialPolygon wide;for(auto q:polygon)wide.push_back(v4bWidePoint(q));return vofMaterialValue(vofMaterialGreen(wide,origin.x,rz));}
// END V4a1b INDEPENDENT AP MATERIAL ORACLE

std::unique_ptr<QuadtreeMesh> meshFor(bool amr=true)
{
    DictionaryMesh d;d.Nx=d.Ny=8;d.maxCellNum=4096;d.maxInitialIterations=0;d.maxLevel_local={2};
    auto mesh=std::make_unique<QuadtreeMesh>(d);
    if(amr){DictionaryMesh::RefineBox r;r.level=1;r.min={.5,0};r.max={1,1};mesh->applyRefineBoxes({r});}
    return mesh;
}
std::array<Point,2> faceEndpoints(const QuadtreeMesh& mesh,const NativeLeafFace& f)
{
    const auto& edge=f.horizontal?mesh.hEdges(f.edgeLabel):mesh.vEdges(f.edgeLabel);
    auto a=mesh.verts(edge.edgeToVert[0]).coord(),b=mesh.verts(edge.edgeToVert[1]).coord();
    if(f.horizontal?a.x>b.x:a.y>b.y)std::swap(a,b);
    return {{{a.x,a.y},{b.x,b.y}}};
}
long double lineVolume(Box b,bool rz,long double intercept)
{
    if(intercept>=b.x1+b.y1)return capacity(b,rz);
    if(intercept<=b.x0+b.y0)return 0;
    const auto full=std::clamp(intercept-b.y1,b.x0,b.x1),end=std::clamp(intercept-b.y0,b.x0,b.x1);
    long double area=(full-b.x0)*(b.y1-b.y0),moment=(full*full-b.x0*b.x0)*(b.y1-b.y0)/2;
    if(end>full){const auto h=intercept-b.y0,second=(end*end-full*full)/2,third=(end*end*end-full*full*full)/3;
        area+=h*(end-full)-second;moment+=h*second-third;}
    return rz?twoPi*moment:area;
}
std::vector<double> lineInput(const QuadtreeMesh& mesh,bool rz,long double intercept=.913L)
{
    std::vector<double> alpha(mesh.cellIndexNum());
    for(int ci=0;ci<mesh.cellIndexNum();++ci){const auto b=bounds(mesh,ci);const auto cap=capacity(b,rz),liquid=lineVolume(b,rz,intercept);
        check(liquid>=0&&liquid<=cap,"independent input outside capacity");
        alpha[ci]=liquid==0?0:(liquid==cap?1:static_cast<double>(liquid/cap));}
    return alpha;
}
Publication publish(const QuadtreeMesh& mesh,bool rz,const std::vector<double>& alpha,Clock clock={7,.125,.25,Position::End})
{ return Publication::prepare(mesh,domain(mesh,rz),clock,alpha); }
// Endpoint differences give exactly the same telescoping discrete-divergence
// construction on uniform/2:1 leaves and XY/RZ. Only input means are prescribed.
long double psi(Point p,int variant)
{
    if(variant==3){
        const auto collar=[](long double x){
            if(x>=.125L&&x<=.875L)return 1.L;
            const auto t=(x<.125L?x:1-x)/.125L;
            return t*t*(3-2*t);
        };
        return (p.y-p.x)*collar(p.x)*collar(p.y);
    }
    const auto base=p.x*p.x*(1-p.x)*(1-p.x)*p.y*p.y*(1-p.y)*(1-p.y);
    return base*(variant==0?1:variant==1?(p.x+.4L*p.y+.2L):(1.3L-p.x+.6L*p.y));
}
std::vector<double> normals(const QuadtreeMesh& mesh,const Publication& p,int variant=0,long double multiplier=1)
{
    std::vector<double> result(p.domain().nativeFaceCount());
    for(int row=0;row<p.domain().nativeFaceCount();++row){const auto f=p.domain().face(mesh,row);const auto ends=faceEndpoints(mesh,f);
        const auto q=(f.horizontal?-1:1)*(psi(ends[1],variant)-psi(ends[0],variant))*multiplier;
        if(f.outwardSign!=0||f.area==0){check(q==0,"streamfunction boundary/axis leakage");result[row]=0;}
        else result[row]=static_cast<double>(q/f.area);}
    return result;
}
double deltaTime(const QuadtreeMesh& mesh,const Publication& p,const std::vector<double>& a,const std::vector<double>& b,long double cfl)
{
    std::vector<long double> out(p.cells().size());
    const auto positiveMean=[](long double begin,long double end){
        if(begin>=0&&end>=0)return (begin+end)/2;
        if(begin<=0&&end<=0)return 0.L;
        const auto positive=begin>0?begin:end;
        return positive*positive/(2*std::abs(end-begin));
    };
    for(int row=0;row<p.domain().nativeFaceCount();++row){const auto f=p.domain().face(mesh,row);
        const auto positive=f.area*positiveMean(a[row],b[row]),negative=f.area*positiveMean(-a[row],-b[row]);
        if(f.negativeCell>=0)out[f.negativeCell]+=positive;
        if(f.positiveCell>=0)out[f.positiveCell]+=negative;}
    long double maximum=0;for(const auto& c:p.cells())maximum=std::max(maximum,out[c.nativeCell]/c.capacity);
    return maximum>0?static_cast<double>(cfl/maximum):.01;
}
Clock nextClock(const Publication& p,double dt)
{ const auto old=p.timeLayer();return {old.step+1,old.time(),old.time()+dt,Position::End}; }
Slab directSlab(const QuadtreeMesh& mesh,const Publication& p,Clock next,const std::vector<double>& a,const std::vector<double>& b)
{ return Slab::prepare(mesh,p,next,a,b); }
Slab hydroSlab(const QuadtreeMesh& mesh,const Publication& p,Clock next,const std::vector<double>& a,const std::vector<double>& b)
{
    const auto layout=HydroDofLayout::build(mesh,p.domain());std::vector<double> av(layout.velocityCount()),bv(av.size());
    for(int row=0;row<layout.velocityCount();++row){av[row]=a[layout.nativeFace(row)];bv[row]=b[layout.nativeFace(row)];}
    const HydroNativeValues begin(layout,av,std::vector<double>(layout.pressureCount(),.125)),
        end(layout,bv,std::vector<double>(layout.pressureCount(),-.25));
    const auto beforeA=begin.normalVelocity(),beforeB=end.normalVelocity(),beforeP=begin.pressure(),beforeQ=end.pressure();
    const auto result=prepareHydroVofAdvectorSlab(mesh,p,next,begin,end);
    check(begin.normalVelocity()==beforeA&&end.normalVelocity()==beforeB&&begin.pressure()==beforeP&&end.pressure()==beforeQ,
        "adapter mutated Hydro native values");
    check(result.beginNormal()==a&&result.endNormal()==b,"adapter did not retain actual native means");
    return result;
}
long double measure(const std::vector<Point>& polygon,Point origin,bool rz)
{
    if(polygon.size()<3)return 0;long double twice=0,six=0;auto a=polygon.back();
    for(auto b:polygon){const auto cross=a.x*b.y-b.x*a.y;twice+=cross;six+=(a.x+b.x)*cross;a=b;}
    if(twice<0){twice=-twice;six=-six;}
    return rz?twoPi*(origin.x*twice/2+six/6):twice/2;
}
void geometryCheck(const QuadtreeMesh& mesh,const Publication& p)
{
    p.validateFor(mesh,p.identity(),p.timeLayer());long double inventory=0,measured=0,total=0;
    for(const auto& c:p.cells()){const auto cap=capacity(bounds(mesh,c.nativeCell),p.domain().axisymmetric());long double l=0,g=0;
        for(const auto& polygon:p.anchoredPhasePolygons(c.nativeCell,true))l+=measure(polygon,p.frameOrigin(),p.domain().axisymmetric());
        for(const auto& polygon:p.anchoredPhasePolygons(c.nativeCell,false))g+=measure(polygon,p.frameOrigin(),p.domain().axisymmetric());
        near(c.capacity,cap,cap,"independent physical capacity");near(l+g,cap,cap,"independent phase complement");
        near(l,c.inventory,cap,"candidate actual geometry/stock",conservationTolerance);
        check(p.alpha()[c.nativeCell]>=0&&p.alpha()[c.nativeCell]<=1,"candidate alpha outside strict bounds");
        if(p.alpha()[c.nativeCell]>0&&p.alpha()[c.nativeCell]<1)check(l>0&&g>0,"candidate deleted a representable phase");
        peakGeometry=std::max(peakGeometry,std::abs(l-c.inventory)/cap);inventory+=c.inventory;measured+=l;total+=cap;}
    near(p.report().inventory,inventory,total,"inventory independent report");
    near(p.report().measured,measured,total,"geometry independent report");
}
void verify(const QuadtreeMesh& mesh,const Publication& old,const Slab& slab,const Transport& t)
{
    check(t.sourceIdentity()==old.identity()&&slab.sourceIdentity()==old.identity(),"foreign transport source owner");
    check(sameClock(t.sourceTime(),old.timeLayer())&&sameClock(slab.sourceTime(),old.timeLayer()),"source complete clock changed");
    check(sameClock(t.endTime(),slab.endTime())&&sameClock(t.endpoint().timeLayer(),slab.endTime()),"End complete clock mismatch");
    auto mid=slab.endTime();mid.position=Position::Midpoint;
    check(sameClock(t.midpoint().timeLayer(),mid),"Midpoint detached from complete global interval");
    check(t.endpoint().identity()!=old.identity()&&t.midpoint().identity()!=old.identity()&&t.endpoint().identity()!=t.midpoint().identity(),
        "transport geometry did not create independent immutable candidates");
    check(t.faces().size()==old.faces().size(),"native face coverage changed");
    const auto count=old.cells().size();std::vector<long double> half(count),full(count),totalHalf(count),totalFull(count);
    std::map<std::pair<int,int>,std::vector<int>> sideRows;
    for(std::size_t row=0;row<t.faces().size();++row){const auto& f=t.faces()[row];const auto native=old.domain().face(mesh,static_cast<int>(row));
        check(f.nativeRow==static_cast<int>(row)&&f.negativeCell==native.negativeCell&&f.positiveCell==native.positiveCell,"native subface incidence mismatch");
        for(const auto* q:{&f.half,&f.full}){
            const auto scale=std::max(1.L,q->positiveTotal+q->negativeTotal);
            check(q->positiveTotal>=0&&q->negativeTotal>=0&&q->positiveLiquid>=0&&q->negativeLiquid>=0&&q->positiveGas>=0&&q->negativeGas>=0,"negative gross phase flux");
            std::ostringstream grossBounds;grossBounds<<std::setprecision(17)<<"phase gross flux exceeds actual total row="<<row
                <<" extent="<<(q==&f.half?"half":"full")<<" positiveTotal="<<q->positiveTotal<<" positiveLiquid="<<q->positiveLiquid
                <<" positiveGas="<<q->positiveGas<<" negativeTotal="<<q->negativeTotal<<" negativeLiquid="<<q->negativeLiquid
                <<" negativeGas="<<q->negativeGas<<" limiter="<<q->limiter;
            check(q->positiveLiquid<=q->positiveTotal&&q->negativeLiquid<=q->negativeTotal&&q->positiveGas<=q->positiveTotal
                &&q->negativeGas<=q->negativeTotal,grossBounds.str());
            near(q->positiveLiquid+q->positiveGas,q->positiveTotal,scale,"independent positive two-phase complement");
            near(q->negativeLiquid+q->negativeGas,q->negativeTotal,scale,"independent negative two-phase complement");
            near(q->gas,q->positiveGas-q->negativeGas,scale,"independent gas gross/signed closure");
            check(q->limiter>=0&&q->limiter<=1,"limiter outside algorithm interval");
            near(q->total,q->positiveTotal-q->negativeTotal,scale,"total gross/signed closure");
            near(q->liquid,q->positiveLiquid-q->negativeLiquid,scale,"liquid gross/signed closure");
            near(q->gas,q->total-q->liquid,scale,"shared gas complement");
            near(q->mass,100*q->liquid+q->gas,100*scale,"same accepted phase mass flux");
            if(native.outwardSign!=0)check(q->total==0&&q->liquid==0&&q->gas==0&&q->mass==0,"closed wall leaked inventory");
        }
        if(f.full.positiveTotal>0&&f.full.negativeTotal>0)++counterflowFaces;
        if(f.full.limiter<1||f.half.limiter<1)++limitedFaces;
        for(auto incidence:{std::pair<int,int>{f.negativeCell,-1},{f.positiveCell,1}})if(incidence.first>=0){
            half[incidence.first]+=incidence.second*f.half.liquid;full[incidence.first]+=incidence.second*f.full.liquid;
            totalHalf[incidence.first]+=incidence.second*f.half.total;totalFull[incidence.first]+=incidence.second*f.full.total;
            sideRows[{incidence.first,static_cast<int>(native.horizontal? (incidence.second>0?FaceSide::South:FaceSide::North)
                :(incidence.second>0?FaceSide::West:FaceSide::East))}].push_back(static_cast<int>(row));}
    }
    long double oldStock=0,halfStock=0,fullStock=0,total=0;
    for(const auto& c:old.cells()){const auto h=t.midpoint().cells()[c.nativeCell].inventory,e=t.endpoint().cells()[c.nativeCell].inventory;
        near(h,c.inventory+half[c.nativeCell],c.capacity,"half inventory consumes final shared flux",conservationTolerance);
        near(e,c.inventory+full[c.nativeCell],c.capacity,"full inventory consumes final shared flux",conservationTolerance);
        near(totalHalf[c.nativeCell],0,c.capacity,"half native divergence",conservationTolerance);
        near(totalFull[c.nativeCell],0,c.capacity,"full native divergence",conservationTolerance);
        peakInventory=std::max({peakInventory,std::abs(h-c.inventory-half[c.nativeCell])/c.capacity,std::abs(e-c.inventory-full[c.nativeCell])/c.capacity});
        oldStock+=c.inventory;halfStock+=h;fullStock+=e;total+=c.capacity;}
    near(halfStock,oldStock,total,"closed half inventory",conservationTolerance);
    near(fullStock,oldStock,total,"closed full inventory",conservationTolerance);
    peakConservation=std::max({peakConservation,std::abs(halfStock-oldStock)/total,std::abs(fullStock-oldStock)/total});
    for(const auto& entry:sideRows)if(entry.second.size()==2){
        const int cell=entry.first.first;const auto& coarse=mesh.cells(mesh.cellIndexToLabel[cell]);
        long double coarseLiquid=0,fineLiquid=0,coarseMass=0,fineMass=0;
        for(int row:entry.second){const auto& f=t.faces()[row];const int other=f.negativeCell==cell?f.positiveCell:f.negativeCell;
            check(other>=0&&mesh.cells(mesh.cellIndexToLabel[other]).level()==coarse.level()+1,"coarse side is not actual native 2:1");
            const int sign=f.negativeCell==cell?1:-1;coarseLiquid+=sign*f.full.liquid;fineLiquid-=sign*f.full.liquid;
            coarseMass+=sign*f.full.mass;fineMass-=sign*f.full.mass;}
        near(coarseLiquid+fineLiquid,0,old.cells()[cell].capacity,"coarse actual subface liquid sum");
        near(coarseMass+fineMass,0,1,"coarse actual subface mass sum");++coarseChecks;
    }
    check(t.report().maximumInventoryResidual<=conservationTolerance&&t.report().halfRelativeConservation<=conservationTolerance
        &&t.report().fullRelativeConservation<=conservationTolerance,"reported stock ledger exceeds frozen gate");
    check(t.report().maximumCapacityReconciliation<=512*std::numeric_limits<long double>::epsilon(),"capacity reconciliation exceeds roundoff budget");
    peakCorrection=std::max(peakCorrection,t.report().maximumFluxCorrection);
    geometryCheck(mesh,t.midpoint());geometryCheck(mesh,t.endpoint());++scenarios;
}
// Independent spatial reference: intersect an Eulerian face line translated
// back at each time with the OLD authoritative convex fan pieces. No production
// clipping, temporal roots, envelopes or quadrature are called.
long double lineFraction(const std::vector<Point>& polygon,Point a,Point d,long double& lo,long double& hi)
{
    lo=0;hi=1;if(polygon.size()<3)return 0;long double area=0;auto prev=polygon.back();
    for(auto p:polygon){area+=prev.x*p.y-p.x*prev.y;prev=p;}
    const int orientation=area>=0?1:-1;prev=polygon.back();
    for(auto p:polygon){const auto ex=p.x-prev.x,ey=p.y-prev.y;
        const auto c=orientation*(ex*(a.y-prev.y)-ey*(a.x-prev.x)),s=orientation*(ex*d.y-ey*d.x);
        if(s>0)lo=std::max(lo,-c/s);else if(s<0)hi=std::min(hi,-c/s);else if(c<0)return 0;
        if(hi<=lo)return 0;prev=p;}
    return hi-lo;
}
// BEGIN V4a1b ACTUAL AP READ FOR ORIGINAL SPATIAL REFERENCE
long double lineFraction(const std::vector<V4bAP>& polygon,Point a,Point d,long double& lo,long double& hi)
{
    lo=0;hi=1;if(polygon.size()<3)return 0;VofMaterialWide area;auto prev=v4bWidePoint(polygon.back());
    for(auto site:polygon){const auto p=v4bWidePoint(site);area=vofMaterialAdd(area,vofMaterialSub(vofMaterialMul(prev[0],p[1]),vofMaterialMul(p[0],prev[1])));prev=p;}
    const int orientation=vofMaterialSign(area)>=0?1:-1;prev=v4bWidePoint(polygon.back());
    for(auto site:polygon){const auto p=v4bWidePoint(site);const auto ex=vofMaterialSub(p[0],prev[0]),ey=vofMaterialSub(p[1],prev[1]);
        const auto c=vofMaterialValue(vofMaterialSub(vofMaterialMul(ex,vofMaterialSub({a.y,0},prev[1])),vofMaterialMul(ey,vofMaterialSub({a.x,0},prev[0]))))*orientation;
        const auto slope=vofMaterialValue(vofMaterialSub(vofMaterialMul(ex,{d.y,0}),vofMaterialMul(ey,{d.x,0})))*orientation;
        if(slope>0)lo=std::max(lo,-c/slope);else if(slope<0)hi=std::min(hi,-c/slope);else if(c<0)return 0;
        if(hi<=lo)return 0;prev=p;
    }return hi-lo;
}
// END V4a1b ACTUAL AP READ FOR ORIGINAL SPATIAL REFERENCE

using Function=std::function<long double(long double)>;
long double adaptive(const Function& f,long double a,long double b,long double fa,long double fm,long double fb,long double whole,long double tolerance,unsigned depth)
{
    const auto m=(a+b)/2,l=(a+m)/2,r=(m+b)/2,fl=f(l),fr=f(r);
    const auto left=(m-a)*(fa+4*fl+fm)/6,right=(b-m)*(fm+4*fr+fb)/6,delta=left+right-whole;
    if(depth==0||std::abs(delta)<=15*tolerance)return left+right+delta/15;
    return adaptive(f,a,m,fa,fl,fm,left,tolerance/2,depth-1)+adaptive(f,m,b,fm,fr,fb,right,tolerance/2,depth-1);
}
long double integrate(const Function& f,long double end,long double tolerance)
{
    long double sum=0;
    for(unsigned k=0;k<32;++k){const auto a=end*k/32,b=end*(k+1)/32,m=(a+b)/2,fa=f(a),fm=f(m),fb=f(b);
        sum+=adaptive(f,a,b,fa,fm,fb,(b-a)*(fa+4*fm+fb)/6,tolerance/32,18);}
    return sum;
}
struct Reference {long double value=0,normalOnly=0;};
Reference rawReference(const QuadtreeMesh& mesh,const Publication& p,const Slab& slab,int row,long double extent)
{
    const auto native=p.domain().face(mesh,row);const auto endpoints=faceEndpoints(mesh,native);const auto origin=p.frameOrigin();
    const Point a{endpoints[0].x-origin.x,endpoints[0].y-origin.y},d{endpoints[1].x-endpoints[0].x,endpoints[1].y-endpoints[0].y};
    const auto v0=slab.beginVelocity()[row],v1=slab.endVelocity()[row];const auto dt=slab.endTime().end-slab.endTime().begin;
    std::vector<std::vector<V4bAP>> polygons;
    // Conservative local support selection avoids an O(N^2) reference scan.
    const auto dx=dt*(std::abs(v0.x)+std::abs(v1.x)),dy=dt*(std::abs(v0.y)+std::abs(v1.y));
    for(const auto& c:p.cells()){const auto b=bounds(mesh,c.nativeCell);
        if(b.x1<endpoints[0].x-dx||b.x0>endpoints[1].x+dx||b.y1<endpoints[0].y-dy||b.y0>endpoints[1].y+dy)continue;
        for(auto polygon:p.anchoredPhasePolygons(c.nativeCell,true))polygons.push_back(std::move(polygon));}
    check(endpoints[0].x-dx>0&&endpoints[1].x+dx<1&&endpoints[0].y-dy>0&&endpoints[1].y+dy<1,"raw reference requires interior old-source support");
    const auto rate=[&](long double tau,bool normalOnly){Point shift{dt*(v0.x*tau+(v1.x-v0.x)*tau*tau/2),
        dt*(v0.y*tau+(v1.y-v0.y)*tau*tau/2)};
        if(normalOnly){if(native.horizontal)shift.x=0;else shift.y=0;}
        const Point departure{a.x-shift.x,a.y-shift.y};long double trace=0;
        std::vector<std::pair<long double,long double>> intervals;
        for(const auto& polygon:polygons){long double lo,hi;if(lineFraction(polygon,departure,d,lo,hi)>0)intervals.emplace_back(lo,hi);}
        std::sort(intervals.begin(),intervals.end());
        // A shared old cell/fan edge can coincide with the query at tau=0.
        // The physical trace is the union, never the sum of duplicate edges.
        std::vector<std::pair<long double,long double>> merged;
        for(auto interval:intervals){
            if(!merged.empty()&&interval.first<=merged.back().second)merged.back().second=std::max(merged.back().second,interval.second);
            else merged.push_back(interval);}
        for(auto interval:merged){const auto lo=interval.first,hi=interval.second;
            trace+=native.tangentSpacing*(hi-lo)*(p.domain().axisymmetric()?twoPi*(endpoints[0].x+d.x*(lo+hi)/2):1);}
        const auto speed=slab.beginNormal()[row]+tau*(slab.endNormal()[row]-slab.beginNormal()[row]);
        return dt*speed*trace;};
    const auto scale=std::max<long double>(1e-18L,dt*native.area*std::max(std::abs(slab.beginNormal()[row]),std::abs(slab.endNormal()[row])));
    return {integrate([&](long double tau){return rate(tau,false);},extent,scale*2e-12L),
        integrate([&](long double tau){return rate(tau,true);},extent,scale*2e-12L)};
}
void checkRaw(const QuadtreeMesh& mesh,const Publication& p,const Slab& slab,const Transport& t,bool requireTransverse)
{
    unsigned count=0;long double difference=0;
    for(std::size_t row=0;row<t.faces().size()&&count<6;++row){const auto native=p.domain().face(mesh,static_cast<int>(row));
        if(native.negativeCell<0||native.positiveCell<0||native.area<=0)continue;
        if(p.alpha()[native.negativeCell]==p.alpha()[native.positiveCell])continue;
        const auto ends=faceEndpoints(mesh,native);const auto v0=slab.beginVelocity()[row],v1=slab.endVelocity()[row];
        const auto dt=slab.endTime().end-slab.endTime().begin;
        const auto dx=dt*(std::abs(v0.x)+std::abs(v1.x)),dy=dt*(std::abs(v0.y)+std::abs(v1.y));
        if(ends[0].x-dx<=0||ends[1].x+dx>=1||ends[0].y-dy<=0||ends[1].y+dy>=1)continue;
        const auto scale=dt*native.area*std::max(std::abs(slab.beginNormal()[row]),std::abs(slab.endNormal()[row]));
        if(scale<=0)continue;
        for(long double extent:{.5L,1.L}){const auto reference=rawReference(mesh,p,slab,static_cast<int>(row),extent);
            const auto actual=extent==.5L?t.faces()[row].half.rawLiquid:t.faces()[row].full.rawLiquid;
            near(actual,reference.value,scale,"independent Eulerian face-time old-fan flux",rawTolerance);
            peakRaw=std::max(peakRaw,std::abs(actual-reference.value)/scale);
            difference=std::max(difference,std::abs(reference.value-reference.normalOnly)/scale);++rawChecks;}
        ++count;}
    check(count>=2,"independent raw reference lacks mixed interior faces");
    if(requireTransverse)check(difference>1e-5L,"raw fixture did not distinguish transverse transport from normal-only strips");
}

// This exercises the generic retained-rational conversion independently of
// the geometric owner and of the budget's demand/limiter construction.
void correctlyRoundedRationals()
{
    namespace X=LAMP::Geometry::Detail::ExactExpansion;
    namespace R=LAMP::Geometry::Detail::ExactRationalRound;
    using E=X::Expansion;
    const auto scalar=[](long double x){
        E e;constexpr unsigned count=(std::numeric_limits<long double>::digits+52)/53+1;
        for(unsigned i=0;x!=0&&i<count;++i){const double part=static_cast<double>(x);
            check(std::isfinite(part)&&part!=0,"test scalar outside exact input window");
            X::accumulate(e,X::fromDouble(part));x-=static_cast<long double>(part);}
        check(x==0,"test scalar did not retain every binary bit");return e;
    };
    const int p=std::numeric_limits<long double>::digits;
    const auto one=scalar(1),two=scalar(2),tail=scalar(std::ldexp(1.L,-p-20));
    const auto negative=[](const E& e){return X::subtract(E{},e);};
    const auto exact=[&](const std::string& name,const E& n,const E& d,long double expected){
        stage="correctly rounded rational "+name;
        const auto actual=R::quotient(n,d);
        std::ostringstream detail;detail<<std::setprecision(std::numeric_limits<long double>::max_digits10)
            <<name<<" actual="<<actual<<" expected="<<expected;
        check(std::isfinite(actual)&&actual==expected,detail.str());++rationalRoundScenarios;
    };
    const auto mirrors=[&](const std::string& name,const E& n,const E& d,long double expected){
        exact(name+" positive",n,d,expected);exact(name+" negative",negative(n),d,-expected);
    };
    exact("mathematical zero",E{},one,0);
    mirrors("unit identity",one,one,1);
    const long double adjacent=std::nextafter(1.L,2.L),second=std::nextafter(adjacent,2.L);
    const long double below=std::nextafter(1.L,0.L);
    const auto firstMid=X::add(one,scalar(std::ldexp(1.L,-p)));
    const auto secondMid=X::add(one,scalar(3*std::ldexp(1.L,-p)));
    const auto lowerMid=X::subtract(one,scalar(std::ldexp(1.L,-p-1)));
    long double third=0;
    for(int k=1;k<=(p+1)/2;++k)third+=std::ldexp(1.L,-2*k);
    if(p%2==0)third+=std::ldexp(1.L,-p-1);
    for(int shift:{-512,0,512}){
        const auto d=X::scalePow2(one,shift);
        const auto run=[&](const std::string& name,const E& n,long double expected){
            mirrors(name+" common_scale=2^"+std::to_string(shift),X::scalePow2(n,shift),d,expected);
        };
        run("first tie even",firstMid,1);
        run("first tie below",X::subtract(firstMid,tail),1);
        run("first tie above",X::add(firstMid,tail),adjacent);
        run("second tie even",secondMid,second);
        run("second tie below",X::subtract(secondMid,tail),adjacent);
        run("second tie above",X::add(secondMid,tail),second);
        run("power boundary tie",lowerMid,1);
        run("power boundary below",X::subtract(lowerMid,tail),below);
        run("power boundary above",X::add(lowerMid,tail),1);
        mirrors("one third common_scale=2^"+std::to_string(shift),
            X::scalePow2(one,shift),X::scalePow2(scalar(3),shift),third);
        mirrors("two thirds common_scale=2^"+std::to_string(shift),
            X::scalePow2(two,shift),X::scalePow2(scalar(3),shift),2*third);
        mirrors("retained micro common_scale=2^"+std::to_string(shift),
            X::scalePow2(scalar(std::ldexp(1.L,-306)),shift),
            X::scalePow2(scalar(std::ldexp(1.L,306)),shift),std::ldexp(1.L,-612));
    }
    // Original row190 total, encoded with its original binary64 bits.
    const long double capturedTotal=0x1.660a54aa27911p-11L;
    mirrors("captured binary total",scalar(capturedTotal),one,capturedTotal);
    for(const auto pair:{std::pair<int,int>{3,1},{7,5},{1023,1022}}){
        const auto d=scalar(pair.first),n=scalar(pair.second),t=scalar(capturedTotal);
        const auto blended=X::add(X::multiply(X::subtract(d,n),t),X::multiply(n,t));
        for(int shift:{-512,0,512})mirrors("pure pair blend D="+std::to_string(pair.first)+
            " N="+std::to_string(pair.second)+" common_scale=2^"+std::to_string(shift),
            X::scalePow2(blended,shift),X::scalePow2(d,shift),capturedTotal);
    }
    if(std::numeric_limits<long double>::digits==std::numeric_limits<double>::digits&&
        std::numeric_limits<long double>::min_exponent==std::numeric_limits<double>::min_exponent&&
        std::numeric_limits<long double>::max_exponent==std::numeric_limits<double>::max_exponent){
        const auto sub=std::numeric_limits<long double>::denorm_min();
        stage="correctly rounded rational positive-to-zero rejection";
        rejects([&]{R::quotient(scalar(sub),two);},"true positive rounds zero");++rationalRoundRejections;
        rejects([&]{R::quotient(negative(scalar(sub)),two);},"true negative rounds zero");++rationalRoundRejections;
        mirrors("subnormal three halves even",scalar(3*sub),two,2*sub);
        mirrors("subnormal five halves even",scalar(5*sub),two,2*sub);
        mirrors("normal/subnormal tie even",X::subtract(one,scalar(std::ldexp(1.L,-p))),
            scalar(std::ldexp(1.L,1-std::numeric_limits<long double>::min_exponent)),
            std::numeric_limits<long double>::min());
        const auto den=scalar(std::ldexp(1.L,1-std::numeric_limits<long double>::max_exponent));
        const auto maxMantissa=X::subtract(two,scalar(std::ldexp(1.L,1-p)));
        mirrors("max finite exact",maxMantissa,den,std::numeric_limits<long double>::max());
        mirrors("max finite quarter ulp",X::add(maxMantissa,scalar(std::ldexp(1.L,-p-1))),
            den,std::numeric_limits<long double>::max());
        const auto overflowTie=X::add(maxMantissa,scalar(std::ldexp(1.L,-p)));
        stage="correctly rounded rational overflow rejection";
        rejects([&]{R::quotient(overflowTie,den);},"positive overflow midpoint");++rationalRoundRejections;
        rejects([&]{R::quotient(negative(overflowTie),den);},"negative overflow midpoint");++rationalRoundRejections;
    }
    std::cout<<"correct_round retained_rational_cases="<<rationalRoundScenarios
        <<" explicit_rejections="<<rationalRoundRejections<<" precision="<<p<<" PASS\n";
}

void capturedBudgetReplay()
{
    using Cell=LAMP::Interface::Detail::GeometricVofPhaseBudgetCell;
    // Complete original arrays, all values encoded as their exact binary64
    // values. Source SHA256: 4DF94BEE38A9C060C95020D9F39D78BBDCC28E3370D943BA307764F753D593EF.
    const std::array<Cell,160> originalCells{{
        {0x1.0000000000000p-6L,0x1.0000000000000p+0}, // native cell 0
        {0x1.0000000000000p-6L,0x1.0000000000000p+0}, // native cell 1
        {0x1.0000000000000p-6L,0x1.0000000000000p+0}, // native cell 2
        {0x1.0000000000000p-6L,0x1.fc504816f0069p-1}, // native cell 3
        {0x1.0000000000000p-6L,0x1.0000000000000p+0}, // native cell 4
        {0x1.0000000000000p-6L,0x1.0000000000000p+0}, // native cell 5
        {0x1.0000000000000p-6L,0x1.fc504816f006ap-1}, // native cell 6
        {0x1.0000000000000p-6L,0x1.8c7e28240b780p-2}, // native cell 7
        {0x1.0000000000000p-6L,0x1.0000000000000p+0}, // native cell 8
        {0x1.0000000000000p-6L,0x1.fc504816f006ap-1}, // native cell 9
        {0x1.0000000000000p-6L,0x1.8c7e28240b780p-2}, // native cell 10
        {0x1.0000000000000p-6L,0x0.0p+0}, // native cell 11
        {0x1.0000000000000p-6L,0x1.fc504816f0069p-1}, // native cell 12
        {0x1.0000000000000p-6L,0x1.8c7e28240b77cp-2}, // native cell 13
        {0x1.0000000000000p-6L,0x0.0p+0}, // native cell 14
        {0x1.0000000000000p-6L,0x0.0p+0}, // native cell 15
        {0x1.0000000000000p-6L,0x1.8c7e28240b77fp-2}, // native cell 16
        {0x1.0000000000000p-6L,0x0.0p+0}, // native cell 17
        {0x1.0000000000000p-6L,0x0.0p+0}, // native cell 18
        {0x1.0000000000000p-6L,0x0.0p+0}, // native cell 19
        {0x1.0000000000000p-6L,0x0.0p+0}, // native cell 20
        {0x1.0000000000000p-6L,0x0.0p+0}, // native cell 21
        {0x1.0000000000000p-6L,0x0.0p+0}, // native cell 22
        {0x1.0000000000000p-6L,0x0.0p+0}, // native cell 23
        {0x1.0000000000000p-6L,0x0.0p+0}, // native cell 24
        {0x1.0000000000000p-6L,0x0.0p+0}, // native cell 25
        {0x1.0000000000000p-6L,0x0.0p+0}, // native cell 26
        {0x1.0000000000000p-6L,0x0.0p+0}, // native cell 27
        {0x1.0000000000000p-6L,0x0.0p+0}, // native cell 28
        {0x1.0000000000000p-6L,0x0.0p+0}, // native cell 29
        {0x1.0000000000000p-6L,0x0.0p+0}, // native cell 30
        {0x1.0000000000000p-6L,0x0.0p+0}, // native cell 31
        {0x1.0000000000000p-8L,0x1.f141205bc01b4p-1}, // native cell 32
        {0x1.0000000000000p-8L,0x1.27bb2fec56d50p-2}, // native cell 33
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 34
        {0x1.0000000000000p-8L,0x1.27bb2fec56d30p-2}, // native cell 35
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 36
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 37
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 38
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 39
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 40
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 41
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 42
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 43
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 44
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 45
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 46
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 47
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 48
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 49
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 50
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 51
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 52
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 53
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 54
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 55
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 56
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 57
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 58
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 59
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 60
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 61
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 62
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 63
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 64
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 65
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 66
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 67
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 68
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 69
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 70
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 71
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 72
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 73
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 74
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 75
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 76
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 77
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 78
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 79
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 80
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 81
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 82
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 83
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 84
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 85
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 86
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 87
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 88
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 89
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 90
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 91
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 92
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 93
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 94
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 95
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 96
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 97
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 98
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 99
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 100
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 101
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 102
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 103
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 104
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 105
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 106
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 107
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 108
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 109
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 110
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 111
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 112
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 113
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 114
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 115
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 116
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 117
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 118
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 119
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 120
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 121
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 122
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 123
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 124
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 125
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 126
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 127
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 128
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 129
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 130
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 131
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 132
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 133
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 134
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 135
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 136
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 137
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 138
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 139
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 140
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 141
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 142
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 143
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 144
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 145
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 146
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 147
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 148
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 149
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 150
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 151
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 152
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 153
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 154
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 155
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 156
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 157
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 158
        {0x1.0000000000000p-8L,0x0.0p+0}, // native cell 159
    }};
    struct CapturedFace {int row,negative,positive;std::array<long double,6> half,full;};
    const std::array<CapturedFace,348> originalFaces{{
        {0,-1,0,{{0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L}},{{0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L}}},
        {1,0,1,{{0x1.3488c7220e160p-13L,0x1.bc4a1482f5904p-17L,0x1.3488c7220e160p-13L,0x1.bc4a1482f5904p-17L,0x0.0p+0L,0x0.0p+0L}},{{0x1.3488c7220e160p-13L,0x1.8aec8402a163cp-12L,0x1.3488c7220e160p-13L,0x1.8aec8402a163cp-12L,0x0.0p+0L,0x0.0p+0L}}},
        {2,1,2,{{0x1.7de853404c156p-12L,0x1.ce0f9a73b7b4dp-15L,0x1.7de853404c156p-12L,0x1.ce0f9a73b7b4dp-15L,0x0.0p+0L,0x0.0p+0L}},{{0x1.7de853404c156p-12L,0x1.2dc10c0ce0760p-10L,0x1.7de853404c156p-12L,0x1.2dc10c0ce0760p-10L,0x0.0p+0L,0x0.0p+0L}}},
        {3,2,3,{{0x1.eb22def24aaadp-12L,0x1.eb22def24aaadp-14L,0x1.eb22def24aaadp-12L,0x1.eb22def24aaadp-14L,0x0.0p+0L,0x0.0p+0L}},{{0x1.eb22def24aaadp-12L,0x1.eb22def24aaadp-10L,0x1.eb22def24aaadp-12L,0x1.eb22def24aaadp-10L,0x0.0p+0L,0x0.0p+0L}}},
        {4,-1,4,{{0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L}},{{0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L}}},
        {5,4,5,{{0x1.57ef6d846cd55p-12L,0x1.2b8a2d97dddcdp-16L,0x1.57ef6d846cd55p-12L,0x1.2b8a2d97dddcdp-16L,0x0.0p+0L,0x0.0p+0L}},{{0x1.57ef6d846cd55p-12L,0x1.71e5413efdfadp-11L,0x1.57ef6d846cd55p-12L,0x1.71e5413efdfadp-11L,0x0.0p+0L,0x0.0p+0L}}},
        {6,5,6,{{0x1.b110b2e49d207p-11L,0x1.472f0d0a75697p-14L,0x1.b110b2e49d207p-11L,0x1.472f0d0a75697p-14L,0x0.0p+0L,0x0.0p+0L}},{{0x1.b110b2e49d207p-11L,0x1.1a3fb1b8dd111p-9L,0x1.b110b2e49d207p-11L,0x1.1a3fb1b8dd111p-9L,0x0.0p+0L,0x0.0p+0L}}},
        {7,6,7,{{0x1.1c78f6894177cp-10L,0x1.68594ecc4ef07p-13L,0x1.022d60dacf761p-10L,0x1.4a57706707fdfp-13L,0x1.a4b95ae7201a9p-14L,0x1.e01de6546f27dp-17L}},{{0x1.1c78f6894177cp-10L,0x1.cab7ea476b3b6p-9L,0x1.022d60dacf761p-10L,0x1.95b8e827dbaecp-9L,0x1.a4b95ae7201a9p-14L,0x1.a7f810fc7c651p-12L}}},
        {8,-1,8,{{0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L}},{{0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L}}},
        {9,8,9,{{0x1.60f36c3ee91c2p-12L,0x1.058b3f4654aa3p-17L,0x1.60f36c3ee91c2p-12L,0x1.058b3f4654aa3p-17L,0x0.0p+0L,0x0.0p+0L}},{{0x1.60f36c3ee91c2p-12L,0x1.2c3de410bc9e9p-11L,0x1.60f36c3ee91c2p-12L,0x1.2c3de410bc9e9p-11L,0x0.0p+0L,0x0.0p+0L}}},
        {10,9,10,{{0x1.c58e5abcac6efp-11L,0x1.4035a4550d6acp-15L,0x1.9070bcc4c26b0p-11L,0x1.1f2e178edf66dp-15L,0x1.a8ecefbf501fap-14L,0x1.083c6631701f9p-18L}},{{0x1.c58e5abcac6efp-11L,0x1.c959f8e1e3ef1p-10L,0x1.9070bcc4c26b0p-11L,0x1.a728c8ab02367p-10L,0x1.a8ecefbf501fap-14L,0x1.118981b70dc4fp-13L}}},
        {11,10,11,{{0x1.316ee13a1f649p-10L,0x1.7d017c171bceap-14L,0x0.0p+0L,0x0.0p+0L,0x1.316ee13a1f649p-10L,0x1.7d017c171bceap-14L}},{{0x1.316ee13a1f649p-10L,0x1.72e8a9898d03ep-9L,0x0.0p+0L,0x0.0p+0L,0x1.316ee13a1f649p-10L,0x1.72e8a9898d03ep-9L}}},
        {12,-1,12,{{0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L}},{{0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L}}},
        {13,12,13,{{0x1.9be64a5677490p-13L,0x0.0p+0L,0x1.68da16e46cfddp-13L,0x0.0p+0L,0x1.98619b9052599p-16L,0x0.0p+0L}},{{0x1.9be8b65066215p-13L,0x1.940d04aabaf2ap-13L,0x1.68dc253573530p-13L,0x1.7bce98a6daa8dp-13L,0x1.986488d796725p-16L,0x1.83e6c03e049cbp-17L}}},
        {14,13,14,{{0x1.12de7e7dddeaap-11L,0x1.b5b225be348fdp-22L,0x0.0p+0L,0x0.0p+0L,0x1.12de7e7dddeaap-11L,0x1.b5b225be348fdp-22L}},{{0x1.12de7e7dddeaap-11L,0x1.3261cc1e119a2p-11L,0x0.0p+0L,0x1.b4d5e325e0e30p-20L,0x1.12de7e7dddeaap-11L,0x1.3187612c7ea9bp-11L}}},
        {15,14,15,{{0x1.831942c060012p-11L,0x1.a4faa2ad2c86ap-19L,0x0.0p+0L,0x0.0p+0L,0x1.831942c060012p-11L,0x1.a4faa2ad2c86ap-19L}},{{0x1.831942c060012p-11L,0x1.ee98ed04d1c84p-11L,0x0.0p+0L,0x0.0p+0L,0x1.831942c060012p-11L,0x1.ee98ed04d1c84p-11L}}},
        {16,-1,16,{{0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L}},{{0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L}}},
        {17,16,17,{{0x1.3cd790dbd1157p-16L,0x1.dff397d32f464p-15L,0x0.0p+0L,0x0.0p+0L,0x1.3cd790dbd1157p-16L,0x1.dff397d32f464p-15L}},{{0x1.15137e7afeb82p-12L,0x1.dff397d32f464p-15L,0x1.b2a4fe727d021p-27L,0x0.0p+0L,0x1.1510193101d32p-12L,0x1.dff397d32f464p-15L}}},
        {18,17,18,{{0x1.2f523bee0089fp-14L,0x1.140740afbc9e7p-13L,0x0.0p+0L,0x0.0p+0L,0x1.2f523bee0089fp-14L,0x1.140740afbc9e7p-13L}},{{0x1.a94572856c3eap-11L,0x1.140740afbc9e7p-13L,0x0.0p+0L,0x0.0p+0L,0x1.a94572856c3eap-11L,0x1.140740afbc9e7p-13L}}},
        {19,18,19,{{0x1.2ea6c41725c3ap-13L,0x1.4232576968ff3p-13L,0x0.0p+0L,0x0.0p+0L,0x1.2ea6c41725c3ap-13L,0x1.4232576968ff3p-13L}},{{0x1.5bbc6240275c9p-10L,0x1.4232576968ff3p-13L,0x0.0p+0L,0x0.0p+0L,0x1.5bbc6240275c9p-10L,0x1.4232576968ff3p-13L}}},
        {20,-1,20,{{0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L}},{{0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L}}},
        {21,20,21,{{0x1.a69de18f4d7ccp-17L,0x1.37927016b5dcap-12L,0x0.0p+0L,0x0.0p+0L,0x1.a69de18f4d7ccp-17L,0x1.37927016b5dcap-12L}},{{0x1.367e7178f502dp-11L,0x1.37927016b5dcap-12L,0x0.0p+0L,0x0.0p+0L,0x1.367e7178f502dp-11L,0x1.37927016b5dcap-12L}}},
        {22,21,22,{{0x1.db34a80be34fap-15L,0x1.8b0ba9d130c97p-11L,0x0.0p+0L,0x0.0p+0L,0x1.db34a80be34fap-15L,0x1.8b0ba9d130c97p-11L}},{{0x1.d98fb68616c99p-10L,0x1.8b0ba9d130c97p-11L,0x0.0p+0L,0x0.0p+0L,0x1.d98fb68616c99p-10L,0x1.8b0ba9d130c97p-11L}}},
        {23,22,23,{{0x1.0b10b176b0937p-13L,0x1.05b951c11daf4p-10L,0x0.0p+0L,0x0.0p+0L,0x1.0b10b176b0937p-13L,0x1.05b951c11daf4p-10L}},{{0x1.8092eb1150307p-9L,0x1.05b951c11daf4p-10L,0x0.0p+0L,0x0.0p+0L,0x1.8092eb1150307p-9L,0x1.05b951c11daf4p-10L}}},
        {24,-1,24,{{0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L}},{{0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L}}},
        {25,24,25,{{0x1.dafa68fd082ebp-18L,0x1.b482c894d1dd0p-12L,0x0.0p+0L,0x0.0p+0L,0x1.dafa68fd082ebp-18L,0x1.b482c894d1dd0p-12L}},{{0x1.5aeee1433bfb9p-11L,0x1.b482c894d1dd0p-12L,0x0.0p+0L,0x0.0p+0L,0x1.5aeee1433bfb9p-11L,0x1.b482c894d1dd0p-12L}}},
        {26,25,26,{{0x1.328a638973963p-15L,0x1.19fac5697dafep-10L,0x0.0p+0L,0x0.0p+0L,0x1.328a638973963p-15L,0x1.19fac5697dafep-10L}},{{0x1.081836976d653p-9L,0x1.19fac5697dafep-10L,0x0.0p+0L,0x0.0p+0L,0x1.081836976d653p-9L,0x1.19fac5697dafep-10L}}},
        {27,26,27,{{0x1.7924b887604fap-14L,0x1.7e3eb7543692bp-10L,0x0.0p+0L,0x0.0p+0L,0x1.7924b887604fap-14L,0x1.7e3eb7543692bp-10L}},{{0x1.ac1bb5ef8a7b9p-9L,0x1.7e3eb7543692bp-10L,0x0.0p+0L,0x0.0p+0L,0x1.ac1bb5ef8a7b9p-9L,0x1.7e3eb7543692bp-10L}}},
        {28,-1,28,{{0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L}},{{0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L}}},
        {29,28,29,{{0x1.f3867d762e622p-20L,0x1.f542833b09369p-13L,0x0.0p+0L,0x0.0p+0L,0x1.f3867d762e622p-20L,0x1.f542833b09369p-13L}},{{0x1.5ae481caf58b8p-12L,0x1.f542833b09369p-13L,0x0.0p+0L,0x0.0p+0L,0x1.5ae481caf58b8p-12L,0x1.f542833b09369p-13L}}},
        {30,29,30,{{0x1.7d26309043cc6p-17L,0x1.47289615f5fcfp-11L,0x0.0p+0L,0x0.0p+0L,0x1.7d26309043cc6p-17L,0x1.47289615f5fcfp-11L}},{{0x1.07c5640b4185bp-10L,0x1.47289615f5fcfp-11L,0x0.0p+0L,0x0.0p+0L,0x1.07c5640b4185bp-10L,0x1.47289615f5fcfp-11L}}},
        {31,30,31,{{0x1.02615f0741fe6p-15L,0x1.c0f66a692c535p-11L,0x0.0p+0L,0x0.0p+0L,0x1.02615f0741fe6p-15L,0x1.c0f66a692c535p-11L}},{{0x1.ab13312bbb5cep-10L,0x1.c0f66a692c535p-11L,0x0.0p+0L,0x0.0p+0L,0x1.ab13312bbb5cep-10L,0x1.c0f66a692c535p-11L}}},
        {32,3,32,{{0x1.dbb55d27f9aa1p-14L,0x1.ca46ab6f7a61ap-15L,0x1.dbb55d27f9aa1p-14L,0x1.ca46ab6f7a61ap-15L,0x0.0p+0L,0x0.0p+0L}},{{0x1.dbb55d27f9aa1p-14L,0x1.531c49ef374e7p-11L,0x1.dbb55d27f9aa1p-14L,0x1.531c49ef374e7p-11L,0x0.0p+0L,0x0.0p+0L}}},
        {33,3,35,{{0x1.47ab6930d5c51p-12L,0x1.fd7dd09516655p-14L,0x1.02fcfebeb49cbp-12L,0x1.bf46c9c42d5b1p-14L,0x1.12b9a9c884a17p-14L,0x1.f1b8368748523p-17L}},{{0x1.47ab6930d5c51p-12L,0x1.9d95adfb381d2p-10L,0x1.02fcfebeb49cbp-12L,0x1.f4aa6c07e5d64p-11L,0x1.12b9a9c884a17p-14L,0x1.4680efee8a640p-11L}}},
        {34,33,36,{{0x1.3c8e2b1097946p-14L,0x1.0269bbe56b455p-14L,0x0.0p+0L,0x0.0p+0L,0x1.3c8e2b1097946p-14L,0x1.0269bbe56b455p-14L}},{{0x1.3c8e2b1097946p-14L,0x1.37c7fc666264ap-11L,0x0.0p+0L,0x0.0p+0L,0x1.3c8e2b1097946p-14L,0x1.37c7fc666264ap-11L}}},
        {35,34,39,{{0x1.bccdcb04db135p-13L,0x1.21b2eeff953adp-13L,0x0.0p+0L,0x0.0p+0L,0x1.bccdcb04db135p-13L,0x1.21b2eeff953adp-13L}},{{0x1.bccdcb04db135p-13L,0x1.7bef4ee199a99p-10L,0x0.0p+0L,0x0.0p+0L,0x1.bccdcb04db135p-13L,0x1.7bef4ee199a99p-10L}}},
        {36,32,33,{{0x1.939dc48abe1a2p-14L,0x1.f834920717d5dp-15L,0x1.5f0f0e345180fp-14L,0x1.b69035b9b9d4ap-15L,0x1.a475b2b364c9bp-17L,0x1.069171357804cp-17L}},{{0x1.939dc48abe1a2p-14L,0x1.4fff351022a63p-11L,0x1.5f0f0e345180fp-14L,0x1.0cba0123951a6p-11L,0x1.a475b2b364c9bp-17L,0x1.0d14cfb2362f3p-13L}}},
        {37,35,34,{{0x1.188d9912e8c41p-12L,0x1.1984b001acb89p-13L,0x0.0p+0L,0x0.0p+0L,0x1.188d9912e8c41p-12L,0x1.1984b001acb89p-13L}},{{0x1.188d9912e8c41p-12L,0x1.999e9108fb0b8p-10L,0x0.0p+0L,0x0.0p+0L,0x1.188d9912e8c41p-12L,0x1.999e9108fb0b8p-10L}}},
        {38,37,40,{{0x1.20cc6df8dedb2p-15L,0x1.a1588f8236942p-15L,0x0.0p+0L,0x0.0p+0L,0x1.20cc6df8dedb2p-15L,0x1.a1588f8236942p-15L}},{{0x1.20cc6df8dedb2p-15L,0x1.a25bf17934ed4p-12L,0x0.0p+0L,0x0.0p+0L,0x1.20cc6df8dedb2p-15L,0x1.a25bf17934ed4p-12L}}},
        {39,38,43,{{0x1.a160207f2f96dp-14L,0x1.d6f9786f8b687p-14L,0x0.0p+0L,0x0.0p+0L,0x1.a160207f2f96dp-14L,0x1.d6f9786f8b687p-14L}},{{0x1.a160207f2f96dp-14L,0x1.fd57921683a78p-11L,0x0.0p+0L,0x0.0p+0L,0x1.a160207f2f96dp-14L,0x1.fd57921683a78p-11L}}},
        {40,36,37,{{0x1.c53174c1c17bfp-15L,0x1.e8b97fcdeac8bp-15L,0x0.0p+0L,0x0.0p+0L,0x1.c53174c1c17bfp-15L,0x1.e8b97fcdeac8bp-15L}},{{0x1.c53174c1c17bfp-15L,0x1.0c2961de32ddcp-11L,0x0.0p+0L,0x0.0p+0L,0x1.c53174c1c17bfp-15L,0x1.0c2961de32ddcp-11L}}},
        {41,39,38,{{0x1.427a5071fa242p-13L,0x1.12e7b91ed999dp-13L,0x0.0p+0L,0x0.0p+0L,0x1.427a5071fa242p-13L,0x1.12e7b91ed999dp-13L}},{{0x1.427a5071fa242p-13L,0x1.46a23cae1d245p-10L,0x0.0p+0L,0x0.0p+0L,0x1.427a5071fa242p-13L,0x1.46a23cae1d245p-10L}}},
        {42,41,44,{{0x1.fff412a135030p-18L,0x1.619aa312cdf40p-16L,0x0.0p+0L,0x0.0p+0L,0x1.fff412a135030p-18L,0x1.619aa312cdf40p-16L}},{{0x1.fff412a135030p-18L,0x1.2b2b551fef1b6p-13L,0x0.0p+0L,0x0.0p+0L,0x1.fff412a135030p-18L,0x1.2b2b551fef1b6p-13L}}},
        {43,42,47,{{0x1.8229d6516040ep-16L,0x1.9126473592e7ep-15L,0x0.0p+0L,0x0.0p+0L,0x1.8229d6516040ep-16L,0x1.9126473592e7ep-15L}},{{0x1.8229d6516040ep-16L,0x1.6bdd0402d7fc1p-12L,0x0.0p+0L,0x0.0p+0L,0x1.8229d6516040ep-16L,0x1.6bdd0402d7fc1p-12L}}},
        {44,40,41,{{0x1.38072c04ae521p-16L,0x1.33df05a07b86ep-15L,0x0.0p+0L,0x0.0p+0L,0x1.38072c04ae521p-16L,0x1.33df05a07b86ep-15L}},{{0x1.38072c04ae521p-16L,0x1.1b04c37fe45acp-12L,0x0.0p+0L,0x0.0p+0L,0x1.38072c04ae521p-16L,0x1.1b04c37fe45acp-12L}}},
        {45,43,42,{{0x1.cb95e7d1d55d7p-15L,0x1.5c6672689daccp-14L,0x0.0p+0L,0x0.0p+0L,0x1.cb95e7d1d55d7p-15L,0x1.5c6672689daccp-14L}},{{0x1.cb95e7d1d55d7p-15L,0x1.5865ffd27c824p-11L,0x0.0p+0L,0x0.0p+0L,0x1.cb95e7d1d55d7p-15L,0x1.5865ffd27c824p-11L}}},
        {46,45,-1,{{0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L}},{{0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L}}},
        {47,46,-1,{{0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L}},{{0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L}}},
        {48,44,45,{{0x1.c30b2056ad0fdp-20L,0x1.c30b2056ad0fdp-18L,0x0.0p+0L,0x0.0p+0L,0x1.c30b2056ad0fdp-20L,0x1.c30b2056ad0fdp-18L}},{{0x1.c30b2056ad0fdp-20L,0x1.6060b143b7346p-15L,0x0.0p+0L,0x0.0p+0L,0x1.c30b2056ad0fdp-20L,0x1.6060b143b7346p-15L}}},
        {49,47,46,{{0x1.5ef878c6b352ep-18L,0x1.006a03165ae1cp-16L,0x0.0p+0L,0x0.0p+0L,0x1.5ef878c6b352ep-18L,0x1.006a03165ae1cp-16L}},{{0x1.5ef878c6b352ep-18L,0x1.ac5838afd9cf6p-14L,0x0.0p+0L,0x0.0p+0L,0x1.5ef878c6b352ep-18L,0x1.ac5838afd9cf6p-14L}}},
        {50,7,48,{{0x1.e42905a5a54a6p-12L,0x1.24e794fc803dap-13L,0x0.0p+0L,0x0.0p+0L,0x1.e42905a5a54a6p-12L,0x1.24e794fc803dap-13L}},{{0x1.e42905a5a54a6p-12L,0x1.0ae32fdc39f1bp-9L,0x0.0p+0L,0x0.0p+0L,0x1.e42905a5a54a6p-12L,0x1.0ae32fdc39f1bp-9L}}},
        {51,7,51,{{0x1.21bbffaaa32b0p-11L,0x1.086d97018553ap-13L,0x0.0p+0L,0x0.0p+0L,0x1.21bbffaaa32b0p-11L,0x1.086d97018553ap-13L}},{{0x1.21bbffaaa32b0p-11L,0x1.14efcb74abfe8p-9L,0x0.0p+0L,0x0.0p+0L,0x1.21bbffaaa32b0p-11L,0x1.14efcb74abfe8p-9L}}},
        {52,49,52,{{0x1.4f5c0f7546559p-12L,0x1.50cb5b5ace795p-13L,0x0.0p+0L,0x0.0p+0L,0x1.4f5c0f7546559p-12L,0x1.50cb5b5ace795p-13L}},{{0x1.4f5c0f7546559p-12L,0x1.e9e0fce80a24dp-10L,0x0.0p+0L,0x0.0p+0L,0x1.4f5c0f7546559p-12L,0x1.e9e0fce80a24dp-10L}}},
        {53,50,55,{{0x1.9959fa0f407cbp-12L,0x1.346cbefe5d627p-13L,0x0.0p+0L,0x0.0p+0L,0x1.9959fa0f407cbp-12L,0x1.346cbefe5d627p-13L}},{{0x1.9959fa0f407cbp-12L,0x1.fbcd204ba6652p-10L,0x0.0p+0L,0x0.0p+0L,0x1.9959fa0f407cbp-12L,0x1.fbcd204ba6652p-10L}}},
        {54,48,49,{{0x1.a273be5dfd32ep-12L,0x1.45975637ffa4bp-13L,0x0.0p+0L,0x0.0p+0L,0x1.a273be5dfd32ep-12L,0x1.45975637ffa4bp-13L}},{{0x1.a273be5dfd32ep-12L,0x1.083495eaf5c0dp-9L,0x0.0p+0L,0x0.0p+0L,0x1.a273be5dfd32ep-12L,0x1.083495eaf5c0dp-9L}}},
        {55,51,50,{{0x1.f96da96998b3fp-12L,0x1.282dc088af073p-13L,0x0.0p+0L,0x0.0p+0L,0x1.f96da96998b3fp-12L,0x1.282dc088af073p-13L}},{{0x1.f96da96998b3fp-12L,0x1.1203fb42ee174p-9L,0x0.0p+0L,0x0.0p+0L,0x1.f96da96998b3fp-12L,0x1.1203fb42ee174p-9L}}},
        {56,53,56,{{0x1.43aa0ed473141p-13L,0x1.142971f07549dp-13L,0x0.0p+0L,0x0.0p+0L,0x1.43aa0ed473141p-13L,0x1.142971f07549dp-13L}},{{0x1.43aa0ed473141p-13L,0x1.4806492dbcffbp-10L,0x0.0p+0L,0x0.0p+0L,0x1.43aa0ed473141p-13L,0x1.4806492dbcffbp-10L}}},
        {57,54,59,{{0x1.95bb1281d2e15p-13L,0x1.ff78f5787d277p-14L,0x0.0p+0L,0x0.0p+0L,0x1.95bb1281d2e15p-13L,0x1.ff78f5787d277p-14L}},{{0x1.95bb1281d2e15p-13L,0x1.53a4bd0d13446p-10L,0x0.0p+0L,0x0.0p+0L,0x1.95bb1281d2e15p-13L,0x1.53a4bd0d13446p-10L}}},
        {58,52,53,{{0x1.ec8c6d803ed96p-13L,0x1.4111e9ad849d6p-13L,0x0.0p+0L,0x0.0p+0L,0x1.ec8c6d803ed96p-13L,0x1.4111e9ad849d6p-13L}},{{0x1.ec8c6d803ed96p-13L,0x1.a4f0732f15f14p-10L,0x0.0p+0L,0x0.0p+0L,0x1.ec8c6d803ed96p-13L,0x1.a4f0732f15f14p-10L}}},
        {59,55,54,{{0x1.304f95eed3e05p-12L,0x1.27c482e99d63cp-13L,0x0.0p+0L,0x0.0p+0L,0x1.304f95eed3e05p-12L,0x1.27c482e99d63cp-13L}},{{0x1.304f95eed3e05p-12L,0x1.b4198dbdf8b84p-10L,0x0.0p+0L,0x0.0p+0L,0x1.304f95eed3e05p-12L,0x1.b4198dbdf8b84p-10L}}},
        {60,57,60,{{0x1.3821d48ec9aa6p-15L,0x1.d9b4aadcccf81p-15L,0x0.0p+0L,0x0.0p+0L,0x1.3821d48ec9aa6p-15L,0x1.d9b4aadcccf81p-15L}},{{0x1.3821d48ec9aa6p-15L,0x1.d421a5ec16559p-12L,0x0.0p+0L,0x0.0p+0L,0x1.3821d48ec9aa6p-15L,0x1.d421a5ec16559p-12L}}},
        {61,58,63,{{0x1.965ec9b37be15p-15L,0x1.ba8f94d9773d4p-15L,0x0.0p+0L,0x0.0p+0L,0x1.965ec9b37be15p-15L,0x1.ba8f94d9773d4p-15L}},{{0x1.965ec9b37be15p-15L,0x1.e41dd6a255380p-12L,0x0.0p+0L,0x0.0p+0L,0x1.965ec9b37be15p-15L,0x1.e41dd6a255380p-12L}}},
        {62,56,57,{{0x1.6b1163ed2bdaep-14L,0x1.9a12fd4636824p-14L,0x0.0p+0L,0x0.0p+0L,0x1.6b1163ed2bdaep-14L,0x1.9a12fd4636824p-14L}},{{0x1.6b1163ed2bdaep-14L,0x1.bb593993d0fbap-11L,0x0.0p+0L,0x0.0p+0L,0x1.6b1163ed2bdaep-14L,0x1.bb593993d0fbap-11L}}},
        {63,59,58,{{0x1.cefcde21fd549p-14L,0x1.7d8556f58236fp-14L,0x0.0p+0L,0x0.0p+0L,0x1.cefcde21fd549p-14L,0x1.7d8556f58236fp-14L}},{{0x1.cefcde21fd549p-14L,0x1.cac6c26ffbd55p-11L,0x0.0p+0L,0x0.0p+0L,0x1.cefcde21fd549p-14L,0x1.cac6c26ffbd55p-11L}}},
        {64,61,-1,{{0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L}},{{0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L}}},
        {65,62,-1,{{0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L}},{{0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L}}},
        {66,60,61,{{0x1.2405d5ded964ap-17L,0x1.2fad4b2defa8ap-16L,0x0.0p+0L,0x0.0p+0L,0x1.2405d5ded964ap-17L,0x1.2fad4b2defa8ap-16L}},{{0x1.2405d5ded964ap-17L,0x1.136024e45a35ep-13L,0x0.0p+0L,0x0.0p+0L,0x1.2405d5ded964ap-17L,0x1.136024e45a35ep-13L}}},
        {67,63,62,{{0x1.86090e8efd9f0p-17L,0x1.1cc4aeb1b3d87p-16L,0x0.0p+0L,0x0.0p+0L,0x1.86090e8efd9f0p-17L,0x1.1cc4aeb1b3d87p-16L}},{{0x1.86090e8efd9f0p-17L,0x1.1c972ccbe6cafp-13L,0x0.0p+0L,0x0.0p+0L,0x1.86090e8efd9f0p-17L,0x1.1c972ccbe6cafp-13L}}},
        {68,11,64,{{0x1.302f8fc56d2c5p-11L,0x1.8b5b9de616cedp-14L,0x0.0p+0L,0x0.0p+0L,0x1.302f8fc56d2c5p-11L,0x1.8b5b9de616cedp-14L}},{{0x1.302f8fc56d2c5p-11L,0x1.f0261ad1ae4f9p-10L,0x0.0p+0L,0x0.0p+0L,0x1.302f8fc56d2c5p-11L,0x1.f0261ad1ae4f9p-10L}}},
        {69,11,67,{{0x1.1d0640f97964dp-11L,0x1.cd933de9f2972p-15L,0x0.0p+0L,0x0.0p+0L,0x1.1d0640f97964dp-11L,0x1.cd933de9f2972p-15L}},{{0x1.1d0640f97964dp-11L,0x1.7d90b42fcb9eep-10L,0x0.0p+0L,0x0.0p+0L,0x1.1d0640f97964dp-11L,0x1.7d90b42fcb9eep-10L}}},
        {70,65,68,{{0x1.b66f16a3d3acbp-12L,0x1.d68a64ada015fp-14L,0x0.0p+0L,0x0.0p+0L,0x1.b66f16a3d3acbp-12L,0x1.d68a64ada015fp-14L}},{{0x1.b66f16a3d3acbp-12L,0x1.c658707c9739bp-10L,0x0.0p+0L,0x0.0p+0L,0x1.b66f16a3d3acbp-12L,0x1.c658707c9739bp-10L}}},
        {71,66,71,{{0x1.a41444db523bfp-12L,0x1.1c504e5659f7ep-14L,0x0.0p+0L,0x0.0p+0L,0x1.a41444db523bfp-12L,0x1.1c504e5659f7ep-14L}},{{0x1.a41444db523bfp-12L,0x1.5ce4fff3074c5p-10L,0x0.0p+0L,0x0.0p+0L,0x1.a41444db523bfp-12L,0x1.5ce4fff3074c5p-10L}}},
        {72,64,65,{{0x1.0bce8e16cb1b1p-11L,0x1.bfa052df1cfb2p-14L,0x0.0p+0L,0x0.0p+0L,0x1.0bce8e16cb1b1p-11L,0x1.bfa052df1cfb2p-14L}},{{0x1.0bce8e16cb1b1p-11L,0x1.eaa240123dd7dp-10L,0x0.0p+0L,0x0.0p+0L,0x1.0bce8e16cb1b1p-11L,0x1.eaa240123dd7dp-10L}}},
        {73,67,66,{{0x1.fb2a5443e9fc9p-12L,0x1.0a2ec140530bep-14L,0x0.0p+0L,0x0.0p+0L,0x1.fb2a5443e9fc9p-12L,0x1.0a2ec140530bep-14L}},{{0x1.fb2a5443e9fc9p-12L,0x1.790c3b5af3217p-10L,0x0.0p+0L,0x0.0p+0L,0x1.fb2a5443e9fc9p-12L,0x1.790c3b5af3217p-10L}}},
        {74,69,72,{{0x1.be2914fcc9fc9p-13L,0x1.8c47447056c19p-14L,0x0.0p+0L,0x0.0p+0L,0x1.be2914fcc9fc9p-13L,0x1.8c47447056c19p-14L}},{{0x1.be2914fcc9fc9p-13L,0x1.2f807dfdd4b22p-10L,0x0.0p+0L,0x0.0p+0L,0x1.be2914fcc9fc9p-13L,0x1.2f807dfdd4b22p-10L}}},
        {75,70,75,{{0x1.b7ed203c8410dp-13L,0x1.eb856efaf2344p-15L,0x0.0p+0L,0x0.0p+0L,0x1.b7ed203c8410dp-13L,0x1.eb856efaf2344p-15L}},{{0x1.b7ed203c8410dp-13L,0x1.d15dc1381da43p-11L,0x0.0p+0L,0x0.0p+0L,0x1.b7ed203c8410dp-13L,0x1.d15dc1381da43p-11L}}},
        {76,68,69,{{0x1.49eed66bd42ffp-12L,0x1.c6f9df9620799p-14L,0x0.0p+0L,0x0.0p+0L,0x1.49eed66bd42ffp-12L,0x1.c6f9df9620799p-14L}},{{0x1.49eed66bd42ffp-12L,0x1.85f2cd85e65fdp-10L,0x0.0p+0L,0x0.0p+0L,0x1.49eed66bd42ffp-12L,0x1.85f2cd85e65fdp-10L}}},
        {77,71,70,{{0x1.4061df172ecabp-12L,0x1.16c0ef3e667bdp-14L,0x0.0p+0L,0x0.0p+0L,0x1.4061df172ecabp-12L,0x1.16c0ef3e667bdp-14L}},{{0x1.4061df172ecabp-12L,0x1.2b34bbe61fb6ap-10L,0x0.0p+0L,0x0.0p+0L,0x1.4061df172ecabp-12L,0x1.2b34bbe61fb6ap-10L}}},
        {78,73,76,{{0x1.cf5f530b30856p-15L,0x1.5b1b28fc5c2b3p-15L,0x0.0p+0L,0x0.0p+0L,0x1.cf5f530b30856p-15L,0x1.5b1b28fc5c2b3p-15L}},{{0x1.cf5f530b30856p-15L,0x1.afff9e824f45ap-12L,0x0.0p+0L,0x0.0p+0L,0x1.cf5f530b30856p-15L,0x1.afff9e824f45ap-12L}}},
        {79,74,79,{{0x1.daaf17bb49c62p-15L,0x1.b74d6fec8bde9p-16L,0x0.0p+0L,0x0.0p+0L,0x1.daaf17bb49c62p-15L,0x1.b74d6fec8bde9p-16L}},{{0x1.daaf17bb49c62p-15L,0x1.4a9c8c1c0a105p-12L,0x0.0p+0L,0x0.0p+0L,0x1.daaf17bb49c62p-15L,0x1.4a9c8c1c0a105p-12L}}},
        {80,72,73,{{0x1.02d19329b99f2p-13L,0x1.29838dbbe2f8fp-14L,0x0.0p+0L,0x0.0p+0L,0x1.02d19329b99f2p-13L,0x1.29838dbbe2f8fp-14L}},{{0x1.02d19329b99f2p-13L,0x1.99adb23f8515dp-11L,0x0.0p+0L,0x0.0p+0L,0x1.02d19329b99f2p-13L,0x1.99adb23f8515dp-11L}}},
        {81,75,74,{{0x1.03be99f40b683p-13L,0x1.74fd7af424ba0p-15L,0x0.0p+0L,0x0.0p+0L,0x1.03be99f40b683p-13L,0x1.74fd7af424ba0p-15L}},{{0x1.03be99f40b683p-13L,0x1.39d02ba3be73cp-11L,0x0.0p+0L,0x0.0p+0L,0x1.03be99f40b683p-13L,0x1.39d02ba3be73cp-11L}}},
        {82,77,-1,{{0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L}},{{0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L}}},
        {83,78,-1,{{0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L}},{{0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L}}},
        {84,76,77,{{0x1.c76e2fe027c90p-17L,0x1.c0fe5ff95f2eap-17L,0x0.0p+0L,0x0.0p+0L,0x1.c76e2fe027c90p-17L,0x1.c0fe5ff95f2eap-17L}},{{0x1.c76e2fe027c90p-17L,0x1.fb86a2a6a9bafp-14L,0x0.0p+0L,0x0.0p+0L,0x1.c76e2fe027c90p-17L,0x1.fb86a2a6a9bafp-14L}}},
        {85,79,78,{{0x1.ddfba49fddcfep-17L,0x1.1e8640e5bdb00p-17L,0x0.0p+0L,0x0.0p+0L,0x1.ddfba49fddcfep-17L,0x1.1e8640e5bdb00p-17L}},{{0x1.ddfba49fddcfep-17L,0x1.840be31dbd06ap-14L,0x0.0p+0L,0x0.0p+0L,0x1.ddfba49fddcfep-17L,0x1.840be31dbd06ap-14L}}},
        {86,15,80,{{0x1.d58b6af087745p-12L,0x1.2cd44fd94881bp-16L,0x0.0p+0L,0x0.0p+0L,0x1.d58b6af087745p-12L,0x1.2cd44fd94881bp-16L}},{{0x1.d58b6af087745p-12L,0x1.cc4b45a39a418p-11L,0x0.0p+0L,0x0.0p+0L,0x1.d58b6af087745p-12L,0x1.cc4b45a39a418p-11L}}},
        {87,15,83,{{0x1.46f32053d3337p-12L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x1.46f32053d3337p-12L,0x0.0p+0L}},{{0x1.48200230aaba8p-12L,0x1.fc8c61d42c7adp-13L,0x0.0p+0L,0x0.0p+0L,0x1.48200230aaba8p-12L,0x1.fc8c61d42c7adp-13L}}},
        {88,81,84,{{0x1.649e27014be67p-12L,0x1.978889cdc577bp-16L,0x0.0p+0L,0x0.0p+0L,0x1.649e27014be67p-12L,0x1.978889cdc577bp-16L}},{{0x1.649e27014be67p-12L,0x1.a3dd228577a20p-11L,0x0.0p+0L,0x0.0p+0L,0x1.649e27014be67p-12L,0x1.a3dd228577a20p-11L}}},
        {89,82,87,{{0x1.08fb6ee8930a9p-12L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x1.08fb6ee8930a9p-12L,0x0.0p+0L}},{{0x1.094af98f6f475p-12L,0x1.cc7067ee6719dp-13L,0x0.0p+0L,0x0.0p+0L,0x1.094af98f6f475p-12L,0x1.cc7067ee6719dp-13L}}},
        {90,80,81,{{0x1.a7c3e015f3f57p-12L,0x1.6d63ab03415f6p-16L,0x0.0p+0L,0x0.0p+0L,0x1.a7c3e015f3f57p-12L,0x1.6d63ab03415f6p-16L}},{{0x1.a7c3e015f3f57p-12L,0x1.c64de8e34a465p-11L,0x0.0p+0L,0x0.0p+0L,0x1.a7c3e015f3f57p-12L,0x1.c64de8e34a465p-11L}}},
        {91,83,82,{{0x1.307a43983837bp-12L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x1.307a43983837bp-12L,0x0.0p+0L}},{{0x1.3129c5f74013ep-12L,0x1.f418d3039bf67p-13L,0x0.0p+0L,0x0.0p+0L,0x1.3129c5f74013ep-12L,0x1.f418d3039bf67p-13L}}},
        {92,85,88,{{0x1.83c95198f43b9p-13L,0x1.79162c34c48a7p-16L,0x0.0p+0L,0x0.0p+0L,0x1.83c95198f43b9p-13L,0x1.79162c34c48a7p-16L}},{{0x1.83c95198f43b9p-13L,0x1.1747ea5e44f03p-11L,0x0.0p+0L,0x0.0p+0L,0x1.83c95198f43b9p-13L,0x1.1747ea5e44f03p-11L}}},
        {93,86,91,{{0x1.36893128498a6p-13L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x1.36893128498a6p-13L,0x0.0p+0L}},{{0x1.368b88b57dcd0p-13L,0x1.2fd6ec38c6d4fp-13L,0x0.0p+0L,0x0.0p+0L,0x1.368b88b57dcd0p-13L,0x1.2fd6ec38c6d4fp-13L}}},
        {94,84,85,{{0x1.14d88d4677edcp-12L,0x1.9e7bf40943fbdp-16L,0x0.0p+0L,0x0.0p+0L,0x1.14d88d4677edcp-12L,0x1.9e7bf40943fbdp-16L}},{{0x1.14d88d4677edcp-12L,0x1.679b25022374ap-11L,0x0.0p+0L,0x0.0p+0L,0x1.14d88d4677edcp-12L,0x1.679b25022374ap-11L}}},
        {95,87,86,{{0x1.aa83e42fb52b3p-13L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x1.aa83e42fb52b3p-13L,0x0.0p+0L}},{{0x1.aab0b6c829e30p-13L,0x1.88d12721435c0p-13L,0x0.0p+0L,0x0.0p+0L,0x1.aab0b6c829e30p-13L,0x1.88d12721435c0p-13L}}},
        {96,89,92,{{0x1.b70b5920a496bp-15L,0x1.6298eb7a03084p-17L,0x0.0p+0L,0x0.0p+0L,0x1.b70b5920a496bp-15L,0x1.6298eb7a03084p-17L}},{{0x1.b70b5920a496bp-15L,0x1.8bb1c7474cc71p-13L,0x0.0p+0L,0x0.0p+0L,0x1.b70b5920a496bp-15L,0x1.8bb1c7474cc71p-13L}}},
        {97,90,95,{{0x1.7ff9257ebe7cep-15L,0x1.22178802d21b4p-25L,0x0.0p+0L,0x0.0p+0L,0x1.7ff9257ebe7cep-15L,0x1.22178802d21b4p-25L}},{{0x1.7ff9257ebe7cep-15L,0x1.aad32a3a27f09p-15L,0x0.0p+0L,0x0.0p+0L,0x1.7ff9257ebe7cep-15L,0x1.aad32a3a27f09p-15L}}},
        {98,88,89,{{0x1.d4560be16398ap-14L,0x1.2607a6b27a930p-16L,0x0.0p+0L,0x0.0p+0L,0x1.d4560be16398ap-14L,0x1.2607a6b27a930p-16L}},{{0x1.d4560be16398ap-14L,0x1.78226e8da5feap-12L,0x0.0p+0L,0x0.0p+0L,0x1.d4560be16398ap-14L,0x1.78226e8da5feap-12L}}},
        {99,91,90,{{0x1.87237c1ad6239p-14L,0x1.57155f2c2669cp-27L,0x0.0p+0L,0x0.0p+0L,0x1.87237c1ad6239p-14L,0x1.57155f2c2669cp-27L}},{{0x1.87237c1ad6239p-14L,0x1.977eda1d47adep-14L,0x0.0p+0L,0x0.0p+0L,0x1.87237c1ad6239p-14L,0x1.977eda1d47adep-14L}}},
        {100,93,-1,{{0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L}},{{0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L}}},
        {101,94,-1,{{0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L}},{{0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L}}},
        {102,92,93,{{0x1.c74a0b653135fp-17L,0x1.d8414dab2747bp-19L,0x0.0p+0L,0x0.0p+0L,0x1.c74a0b653135fp-17L,0x1.d8414dab2747bp-19L}},{{0x1.c74a0b653135fp-17L,0x1.cfbbbe22f08f9p-15L,0x0.0p+0L,0x0.0p+0L,0x1.c74a0b653135fp-17L,0x1.cfbbbe22f08f9p-15L}}},
        {103,95,94,{{0x1.a2eff8a999551p-17L,0x1.b4babbe65aa8cp-26L,0x0.0p+0L,0x0.0p+0L,0x1.a2eff8a999551p-17L,0x1.b4babbe65aa8cp-26L}},{{0x1.a2eff8a999551p-17L,0x1.f1f6c62b21c77p-17L,0x0.0p+0L,0x0.0p+0L,0x1.a2eff8a999551p-17L,0x1.f1f6c62b21c77p-17L}}},
        {104,19,96,{{0x1.2cde8f7917cf4p-13L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x1.2cde8f7917cf4p-13L,0x0.0p+0L}},{{0x1.1954b2124a7e3p-11L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x1.1954b2124a7e3p-11L,0x0.0p+0L}}},
        {105,19,99,{{0x1.ad895313ef55ap-14L,0x1.530401b4698fbp-13L,0x0.0p+0L,0x0.0p+0L,0x1.ad895313ef55ap-14L,0x1.530401b4698fbp-13L}},{{0x1.1cad70eb7673ap-10L,0x1.530401b4698fbp-13L,0x0.0p+0L,0x0.0p+0L,0x1.1cad70eb7673ap-10L,0x1.530401b4698fbp-13L}}},
        {106,97,100,{{0x1.35301d1ba117ep-13L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x1.35301d1ba117ep-13L,0x0.0p+0L}},{{0x1.062e462d79dc8p-11L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x1.062e462d79dc8p-11L,0x0.0p+0L}}},
        {107,98,103,{{0x1.dfe913400ced9p-14L,0x1.b64268624958cp-14L,0x0.0p+0L,0x0.0p+0L,0x1.dfe913400ced9p-14L,0x1.b64268624958cp-14L}},{{0x1.06059bbc82e90p-10L,0x1.b64268624958cp-14L,0x0.0p+0L,0x0.0p+0L,0x1.06059bbc82e90p-10L,0x1.b64268624958cp-14L}}},
        {108,96,97,{{0x1.3c35cebc49dfdp-13L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x1.3c35cebc49dfdp-13L,0x0.0p+0L}},{{0x1.18d4c3beda637p-11L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x1.18d4c3beda637p-11L,0x0.0p+0L}}},
        {109,99,98,{{0x1.d64159405b2c6p-14L,0x1.1bd8577ed94a3p-13L,0x0.0p+0L,0x0.0p+0L,0x1.d64159405b2c6p-14L,0x1.1bd8577ed94a3p-13L}},{{0x1.1a36f758d150ep-10L,0x1.1bd8577ed94a3p-13L,0x0.0p+0L,0x0.0p+0L,0x1.1a36f758d150ep-10L,0x1.1bd8577ed94a3p-13L}}},
        {110,101,104,{{0x1.c50b07d442e61p-14L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x1.c50b07d442e61p-14L,0x0.0p+0L}},{{0x1.62af7fcdebd79p-12L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x1.62af7fcdebd79p-12L,0x0.0p+0L}}},
        {111,102,107,{{0x1.80aafe441d89ap-14L,0x1.7f27086be4d76p-15L,0x0.0p+0L,0x0.0p+0L,0x1.80aafe441d89ap-14L,0x1.7f27086be4d76p-15L}},{{0x1.60037431ebec8p-11L,0x1.7f27086be4d76p-15L,0x0.0p+0L,0x0.0p+0L,0x1.60037431ebec8p-11L,0x1.7f27086be4d76p-15L}}},
        {112,100,101,{{0x1.169b15c5794e3p-13L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x1.169b15c5794e3p-13L,0x0.0p+0L}},{{0x1.c524ded0010d6p-12L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x1.c524ded0010d6p-12L,0x0.0p+0L}}},
        {113,103,102,{{0x1.c40a04120692cp-14L,0x1.33c814a854a7fp-14L,0x0.0p+0L,0x0.0p+0L,0x1.c40a04120692cp-14L,0x1.33c814a854a7fp-14L}},{{0x1.c2fe123427834p-11L,0x1.33c814a854a7fp-14L,0x0.0p+0L,0x0.0p+0L,0x1.c2fe123427834p-11L,0x1.33c814a854a7fp-14L}}},
        {114,105,108,{{0x1.5b4e0fbc53d7bp-15L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x1.5b4e0fbc53d7bp-15L,0x0.0p+0L}},{{0x1.fcc17743c079ep-14L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x1.fcc17743c079ep-14L,0x0.0p+0L}}},
        {115,106,111,{{0x1.43fd78f3786f6p-15L,0x1.3d3fb2685df4cp-17L,0x0.0p+0L,0x0.0p+0L,0x1.43fd78f3786f6p-15L,0x1.3d3fb2685df4cp-17L}},{{0x1.f81e7a59806dap-13L,0x1.3d3fb2685df4cp-17L,0x0.0p+0L,0x0.0p+0L,0x1.f81e7a59806dap-13L,0x1.3d3fb2685df4cp-17L}}},
        {116,104,105,{{0x1.3dfca7a27480ep-14L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x1.3dfca7a27480ep-14L,0x0.0p+0L}},{{0x1.e0e6a10d6ef20p-13L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x1.e0e6a10d6ef20p-13L,0x0.0p+0L}}},
        {117,107,106,{{0x1.1ae0ae059cc64p-14L,0x1.91e1aedb28f4bp-16L,0x0.0p+0L,0x0.0p+0L,0x1.1ae0ae059cc64p-14L,0x1.91e1aedb28f4bp-16L}},{{0x1.dc947d59877c2p-12L,0x1.91e1aedb28f4bp-16L,0x0.0p+0L,0x0.0p+0L,0x1.dc947d59877c2p-12L,0x1.91e1aedb28f4bp-16L}}},
        {118,109,-1,{{0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L}},{{0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L}}},
        {119,110,-1,{{0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L}},{{0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L}}},
        {120,108,109,{{0x1.a5115fa121055p-17L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x1.a5115fa121055p-17L,0x0.0p+0L}},{{0x1.2b795c59f9d65p-15L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x1.2b795c59f9d65p-15L,0x0.0p+0L}}},
        {121,111,110,{{0x1.9c3d32b439030p-17L,0x1.096161db7dfadp-19L,0x0.0p+0L,0x0.0p+0L,0x1.9c3d32b439030p-17L,0x1.096161db7dfadp-19L}},{{0x1.291a1667a38cep-14L,0x1.096161db7dfadp-19L,0x0.0p+0L,0x0.0p+0L,0x1.291a1667a38cep-14L,0x1.096161db7dfadp-19L}}},
        {122,23,112,{{0x1.af82924f7cb7cp-14L,0x1.8ab52e9e3d1e5p-12L,0x0.0p+0L,0x0.0p+0L,0x1.af82924f7cb7cp-14L,0x1.8ab52e9e3d1e5p-12L}},{{0x1.9ce764be65b63p-10L,0x1.8ab52e9e3d1e5p-12L,0x0.0p+0L,0x0.0p+0L,0x1.9ce764be65b63p-10L,0x1.8ab52e9e3d1e5p-12L}}},
        {123,23,115,{{0x1.9a64cc9a15fcdp-14L,0x1.2995ed6c6bdd4p-11L,0x0.0p+0L,0x0.0p+0L,0x1.9a64cc9a15fcdp-14L,0x1.2995ed6c6bdd4p-11L}},{{0x1.f28079c4d9209p-10L,0x1.2995ed6c6bdd4p-11L,0x0.0p+0L,0x0.0p+0L,0x1.f28079c4d9209p-10L,0x1.2995ed6c6bdd4p-11L}}},
        {124,113,116,{{0x1.f29b02234afa1p-14L,0x1.13708e18776fep-12L,0x0.0p+0L,0x0.0p+0L,0x1.f29b02234afa1p-14L,0x1.13708e18776fep-12L}},{{0x1.7ace337f02f2dp-10L,0x1.13708e18776fep-12L,0x0.0p+0L,0x0.0p+0L,0x1.7ace337f02f2dp-10L,0x1.13708e18776fep-12L}}},
        {125,114,119,{{0x1.e6984210f124ap-14L,0x1.ab8aba1926d16p-12L,0x0.0p+0L,0x0.0p+0L,0x1.e6984210f124ap-14L,0x1.ab8aba1926d16p-12L}},{{0x1.c8974aa5c1b7ep-10L,0x1.ab8aba1926d16p-12L,0x0.0p+0L,0x0.0p+0L,0x1.c8974aa5c1b7ep-10L,0x1.ab8aba1926d16p-12L}}},
        {126,112,113,{{0x1.e0ea2b6e4048fp-14L,0x1.5651b2cc66c37p-12L,0x0.0p+0L,0x0.0p+0L,0x1.e0ea2b6e4048fp-14L,0x1.5651b2cc66c37p-12L}},{{0x1.98ade6fe11d5ep-10L,0x1.5651b2cc66c37p-12L,0x0.0p+0L,0x0.0p+0L,0x1.98ade6fe11d5ep-10L,0x1.5651b2cc66c37p-12L}}},
        {127,115,114,{{0x1.cfb670bac870ep-14L,0x1.0599d837b8666p-11L,0x0.0p+0L,0x0.0p+0L,0x1.cfb670bac870ep-14L,0x1.0599d837b8666p-11L}},{{0x1.ed02363464951p-10L,0x1.0599d837b8666p-11L,0x0.0p+0L,0x0.0p+0L,0x1.ed02363464951p-10L,0x1.0599d837b8666p-11L}}},
        {128,117,120,{{0x1.9a6a21bc9cd7ep-14L,0x1.0c8a932b5496dp-13L,0x0.0p+0L,0x0.0p+0L,0x1.9a6a21bc9cd7ep-14L,0x1.0c8a932b5496dp-13L}},{{0x1.fb1734041238bp-11L,0x1.0c8a932b5496dp-13L,0x0.0p+0L,0x0.0p+0L,0x1.fb1734041238bp-11L,0x1.0c8a932b5496dp-13L}}},
        {129,118,123,{{0x1.989d443399d0bp-14L,0x1.b13fd27388cdep-13L,0x0.0p+0L,0x0.0p+0L,0x1.989d443399d0bp-14L,0x1.b13fd27388cdep-13L}},{{0x1.311160be7d000p-10L,0x1.b13fd27388cdep-13L,0x0.0p+0L,0x0.0p+0L,0x1.311160be7d000p-10L,0x1.b13fd27388cdep-13L}}},
        {130,116,117,{{0x1.dc4b18210d187p-14L,0x1.966d6b0dadb42p-13L,0x0.0p+0L,0x0.0p+0L,0x1.dc4b18210d187p-14L,0x1.966d6b0dadb42p-13L}},{{0x1.456e82c50a426p-10L,0x1.966d6b0dadb42p-13L,0x0.0p+0L,0x0.0p+0L,0x1.456e82c50a426p-10L,0x1.966d6b0dadb42p-13L}}},
        {131,119,118,{{0x1.d5c68f8ee5836p-14L,0x1.411ba77adaba9p-12L,0x0.0p+0L,0x0.0p+0L,0x1.d5c68f8ee5836p-14L,0x1.411ba77adaba9p-12L}},{{0x1.87eacac40b951p-10L,0x1.411ba77adaba9p-12L,0x0.0p+0L,0x0.0p+0L,0x1.87eacac40b951p-10L,0x1.411ba77adaba9p-12L}}},
        {132,121,124,{{0x1.6111c268e6697p-15L,0x1.06c59d6b8dc03p-15L,0x0.0p+0L,0x0.0p+0L,0x1.6111c268e6697p-15L,0x1.06c59d6b8dc03p-15L}},{{0x1.69ad72f2774f9p-12L,0x1.06c59d6b8dc03p-15L,0x0.0p+0L,0x0.0p+0L,0x1.69ad72f2774f9p-12L,0x1.06c59d6b8dc03p-15L}}},
        {133,122,127,{{0x1.651a04eccbfcdp-15L,0x1.bf6629ea24a61p-15L,0x0.0p+0L,0x0.0p+0L,0x1.651a04eccbfcdp-15L,0x1.bf6629ea24a61p-15L}},{{0x1.b2547a366d11ep-12L,0x1.bf6629ea24a61p-15L,0x0.0p+0L,0x0.0p+0L,0x1.b2547a366d11ep-12L,0x1.bf6629ea24a61p-15L}}},
        {134,120,121,{{0x1.3134360359862p-14L,0x1.2f38090cfeb21p-14L,0x0.0p+0L,0x0.0p+0L,0x1.3134360359862p-14L,0x1.2f38090cfeb21p-14L}},{{0x1.569bf6d915fbep-11L,0x1.2f38090cfeb21p-14L,0x0.0p+0L,0x0.0p+0L,0x1.569bf6d915fbep-11L,0x1.2f38090cfeb21p-14L}}},
        {135,123,122,{{0x1.3268568366222p-14L,0x1.f5572f69d7130p-14L,0x0.0p+0L,0x0.0p+0L,0x1.3268568366222p-14L,0x1.f5572f69d7130p-14L}},{{0x1.9bd6fdffd6c91p-11L,0x1.f5572f69d7130p-14L,0x0.0p+0L,0x0.0p+0L,0x1.9bd6fdffd6c91p-11L,0x1.f5572f69d7130p-14L}}},
        {136,125,-1,{{0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L}},{{0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L}}},
        {137,126,-1,{{0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L}},{{0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L}}},
        {138,124,125,{{0x1.c54362c108cf2p-17L,0x1.f097fe1df5a9fp-18L,0x0.0p+0L,0x0.0p+0L,0x1.c54362c108cf2p-17L,0x1.f097fe1df5a9fp-18L}},{{0x1.a9680fa85300fp-14L,0x1.f097fe1df5a9fp-18L,0x0.0p+0L,0x0.0p+0L,0x1.a9680fa85300fp-14L,0x1.f097fe1df5a9fp-18L}}},
        {139,127,126,{{0x1.cd7c341d228afp-17L,0x1.b61c5c3ce8e84p-17L,0x0.0p+0L,0x0.0p+0L,0x1.cd7c341d228afp-17L,0x1.b61c5c3ce8e84p-17L}},{{0x1.fe54597e1793bp-14L,0x1.b61c5c3ce8e84p-17L,0x0.0p+0L,0x0.0p+0L,0x1.fe54597e1793bp-14L,0x1.b61c5c3ce8e84p-17L}}},
        {140,27,128,{{0x1.64cc9598bccabp-14L,0x1.6d1c8354a0667p-11L,0x0.0p+0L,0x0.0p+0L,0x1.64cc9598bccabp-14L,0x1.6d1c8354a0667p-11L}},{{0x1.077c7a2213a5bp-9L,0x1.6d1c8354a0667p-11L,0x0.0p+0L,0x0.0p+0L,0x1.077c7a2213a5bp-9L,0x1.6d1c8354a0667p-11L}}},
        {141,27,131,{{0x1.12b2c9bba3663p-14L,0x1.79a7a604069f2p-11L,0x0.0p+0L,0x0.0p+0L,0x1.12b2c9bba3663p-14L,0x1.79a7a604069f2p-11L}},{{0x1.e540d1d5ce229p-10L,0x1.79a7a604069f2p-11L,0x0.0p+0L,0x0.0p+0L,0x1.e540d1d5ce229p-10L,0x1.79a7a604069f2p-11L}}},
        {142,129,132,{{0x1.b126029af88b5p-14L,0x1.0add7bd27d564p-11L,0x0.0p+0L,0x0.0p+0L,0x1.b126029af88b5p-14L,0x1.0add7bd27d564p-11L}},{{0x1.e220c0e1b15c5p-10L,0x1.0add7bd27d564p-11L,0x0.0p+0L,0x0.0p+0L,0x1.e220c0e1b15c5p-10L,0x1.0add7bd27d564p-11L}}},
        {143,130,135,{{0x1.5576938be6cb2p-14L,0x1.1780535668ee7p-11L,0x0.0p+0L,0x0.0p+0L,0x1.5576938be6cb2p-14L,0x1.1780535668ee7p-11L}},{{0x1.bb90980a9a00dp-10L,0x1.1780535668ee7p-11L,0x0.0p+0L,0x0.0p+0L,0x1.bb90980a9a00dp-10L,0x1.1780535668ee7p-11L}}},
        {144,128,129,{{0x1.9848b9878921dp-14L,0x1.4395cdea89f04p-11L,0x0.0p+0L,0x0.0p+0L,0x1.9848b9878921dp-14L,0x1.4395cdea89f04p-11L}},{{0x1.04709e2e0ac0ep-9L,0x1.4395cdea89f04p-11L,0x0.0p+0L,0x0.0p+0L,0x1.04709e2e0ac0ep-9L,0x1.4395cdea89f04p-11L}}},
        {145,131,130,{{0x1.3e5e0aacc9e3dp-14L,0x1.50adac20a41a6p-11L,0x0.0p+0L,0x0.0p+0L,0x1.3e5e0aacc9e3dp-14L,0x1.50adac20a41a6p-11L}},{{0x1.df6f26dd60ce5p-10L,0x1.50adac20a41a6p-11L,0x0.0p+0L,0x0.0p+0L,0x1.df6f26dd60ce5p-10L,0x1.50adac20a41a6p-11L}}},
        {146,133,136,{{0x1.7251e02725be1p-14L,0x1.1492657ac36c5p-12L,0x0.0p+0L,0x0.0p+0L,0x1.7251e02725be1p-14L,0x1.1492657ac36c5p-12L}},{{0x1.41bd1396cf26ap-10L,0x1.1492657ac36c5p-12L,0x0.0p+0L,0x0.0p+0L,0x1.41bd1396cf26ap-10L,0x1.1492657ac36c5p-12L}}},
        {147,134,139,{{0x1.2922a9e2af570p-14L,0x1.264e7e0aa64b1p-12L,0x0.0p+0L,0x0.0p+0L,0x1.2922a9e2af570p-14L,0x1.264e7e0aa64b1p-12L}},{{0x1.27b8252084554p-10L,0x1.264e7e0aa64b1p-12L,0x0.0p+0L,0x0.0p+0L,0x1.27b8252084554p-10L,0x1.264e7e0aa64b1p-12L}}},
        {148,132,133,{{0x1.a63122db5f576p-14L,0x1.9513ed8652b66p-12L,0x0.0p+0L,0x0.0p+0L,0x1.a63122db5f576p-14L,0x1.9513ed8652b66p-12L}},{{0x1.9d9733875bdd0p-10L,0x1.9513ed8652b66p-12L,0x0.0p+0L,0x0.0p+0L,0x1.9d9733875bdd0p-10L,0x1.9513ed8652b66p-12L}}},
        {149,135,134,{{0x1.4ffb680e5b56dp-14L,0x1.ab6d297cd53ccp-12L,0x0.0p+0L,0x0.0p+0L,0x1.4ffb680e5b56dp-14L,0x1.ab6d297cd53ccp-12L}},{{0x1.7c5480ae4da24p-10L,0x1.ab6d297cd53ccp-12L,0x0.0p+0L,0x0.0p+0L,0x1.7c5480ae4da24p-10L,0x1.ab6d297cd53ccp-12L}}},
        {150,137,140,{{0x1.482fe126ad29cp-15L,0x1.265a1fc8fd915p-14L,0x0.0p+0L,0x0.0p+0L,0x1.482fe126ad29cp-15L,0x1.265a1fc8fd915p-14L}},{{0x1.c974fe7098e76p-12L,0x1.265a1fc8fd915p-14L,0x0.0p+0L,0x0.0p+0L,0x1.c974fe7098e76p-12L,0x1.265a1fc8fd915p-14L}}},
        {151,138,143,{{0x1.0aed1c326ff74p-15L,0x1.3fd2df414cbcbp-14L,0x0.0p+0L,0x0.0p+0L,0x1.0aed1c326ff74p-15L,0x1.3fd2df414cbcbp-14L}},{{0x1.a4058f1033822p-12L,0x1.3fd2df414cbcbp-14L,0x0.0p+0L,0x0.0p+0L,0x1.a4058f1033822p-12L,0x1.3fd2df414cbcbp-14L}}},
        {152,136,137,{{0x1.17c17162903adp-14L,0x1.4484ad9540d4bp-13L,0x0.0p+0L,0x0.0p+0L,0x1.17c17162903adp-14L,0x1.4484ad9540d4bp-13L}},{{0x1.b21041a1cdfd2p-11L,0x1.4484ad9540d4bp-13L,0x0.0p+0L,0x0.0p+0L,0x1.b21041a1cdfd2p-11L,0x1.4484ad9540d4bp-13L}}},
        {153,139,138,{{0x1.c42fd0d7f36a4p-15L,0x1.5caef14db1ab3p-13L,0x0.0p+0L,0x0.0p+0L,0x1.c42fd0d7f36a4p-15L,0x1.5caef14db1ab3p-13L}},{{0x1.8ec186282f345p-11L,0x1.5caef14db1ab3p-13L,0x0.0p+0L,0x0.0p+0L,0x1.8ec186282f345p-11L,0x1.5caef14db1ab3p-13L}}},
        {154,141,-1,{{0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L}},{{0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L}}},
        {155,142,-1,{{0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L}},{{0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L}}},
        {156,140,141,{{0x1.aa99d565a7b00p-17L,0x1.25da9fb6c98cdp-16L,0x0.0p+0L,0x0.0p+0L,0x1.aa99d565a7b00p-17L,0x1.25da9fb6c98cdp-16L}},{{0x1.0c8fa0721e586p-13L,0x1.25da9fb6c98cdp-16L,0x0.0p+0L,0x0.0p+0L,0x1.0c8fa0721e586p-13L,0x1.25da9fb6c98cdp-16L}}},
        {157,143,142,{{0x1.5ceda73d0396ap-17L,0x1.43791361f7e71p-16L,0x0.0p+0L,0x0.0p+0L,0x1.5ceda73d0396ap-17L,0x1.43791361f7e71p-16L}},{{0x1.ece4493a0432dp-14L,0x1.43791361f7e71p-16L,0x0.0p+0L,0x0.0p+0L,0x1.ece4493a0432dp-14L,0x1.43791361f7e71p-16L}}},
        {158,31,144,{{0x1.57f5987b8ee20p-15L,0x1.35d236a388b4ap-11L,0x0.0p+0L,0x0.0p+0L,0x1.57f5987b8ee20p-15L,0x1.35d236a388b4ap-11L}},{{0x1.6920ad476c907p-10L,0x1.35d236a388b4ap-11L,0x0.0p+0L,0x0.0p+0L,0x1.6920ad476c907p-10L,0x1.35d236a388b4ap-11L}}},
        {159,31,147,{{0x1.d131f608771b7p-17L,0x1.099e42f868539p-12L,0x0.0p+0L,0x0.0p+0L,0x1.d131f608771b7p-17L,0x1.099e42f868539p-12L}},{{0x1.1e29f429089abp-11L,0x1.099e42f868539p-12L,0x0.0p+0L,0x0.0p+0L,0x1.1e29f429089abp-11L,0x1.099e42f868539p-12L}}},
        {160,145,148,{{0x1.b64bceafe43bfp-15L,0x1.cf1dce2458f83p-12L,0x0.0p+0L,0x0.0p+0L,0x1.b64bceafe43bfp-15L,0x1.cf1dce2458f83p-12L}},{{0x1.49dabc462c03bp-10L,0x1.cf1dce2458f83p-12L,0x0.0p+0L,0x0.0p+0L,0x1.49dabc462c03bp-10L,0x1.cf1dce2458f83p-12L}}},
        {161,146,151,{{0x1.2fc096e7d6adcp-16L,0x1.900981528ec1fp-13L,0x0.0p+0L,0x0.0p+0L,0x1.2fc096e7d6adcp-16L,0x1.900981528ec1fp-13L}},{{0x1.0538da639e8bap-11L,0x1.900981528ec1fp-13L,0x0.0p+0L,0x0.0p+0L,0x1.0538da639e8bap-11L,0x1.900981528ec1fp-13L}}},
        {162,144,145,{{0x1.940761b1d73ffp-15L,0x1.157d54f3a77b7p-11L,0x0.0p+0L,0x0.0p+0L,0x1.940761b1d73ffp-15L,0x1.157d54f3a77b7p-11L}},{{0x1.64aa565e69609p-10L,0x1.157d54f3a77b7p-11L,0x0.0p+0L,0x0.0p+0L,0x1.64aa565e69609p-10L,0x1.157d54f3a77b7p-11L}}},
        {163,147,146,{{0x1.14e459b140dc0p-16L,0x1.dd7f56c15856ep-13L,0x0.0p+0L,0x0.0p+0L,0x1.14e459b140dc0p-16L,0x1.dd7f56c15856ep-13L}},{{0x1.1a8af841ee975p-11L,0x1.dd7f56c15856ep-13L,0x0.0p+0L,0x0.0p+0L,0x1.1a8af841ee975p-11L,0x1.dd7f56c15856ep-13L}}},
        {164,149,152,{{0x1.8462c240ded45p-15L,0x1.edb102699b19cp-13L,0x0.0p+0L,0x0.0p+0L,0x1.8462c240ded45p-15L,0x1.edb102699b19cp-13L}},{{0x1.b77624b07eb6ap-11L,0x1.edb102699b19cp-13L,0x0.0p+0L,0x0.0p+0L,0x1.b77624b07eb6ap-11L,0x1.edb102699b19cp-13L}}},
        {165,150,155,{{0x1.11f0e3bdb5e8bp-16L,0x1.ae72975df1f5ap-14L,0x0.0p+0L,0x0.0p+0L,0x1.11f0e3bdb5e8bp-16L,0x1.ae72975df1f5ap-14L}},{{0x1.5bcae4808571fp-12L,0x1.ae72975df1f5ap-14L,0x0.0p+0L,0x0.0p+0L,0x1.5bcae4808571fp-12L,0x1.ae72975df1f5ap-14L}}},
        {166,148,149,{{0x1.b37b401b5b57fp-15L,0x1.642a7fdad5cbep-12L,0x0.0p+0L,0x0.0p+0L,0x1.b37b401b5b57fp-15L,0x1.642a7fdad5cbep-12L}},{{0x1.1ab7a3dc9ff43p-10L,0x1.642a7fdad5cbep-12L,0x0.0p+0L,0x0.0p+0L,0x1.1ab7a3dc9ff43p-10L,0x1.642a7fdad5cbep-12L}}},
        {167,151,150,{{0x1.30acc78c907ffp-16L,0x1.3502ed2e4cf93p-13L,0x0.0p+0L,0x0.0p+0L,0x1.30acc78c907ffp-16L,0x1.3502ed2e4cf93p-13L}},{{0x1.bfa3d1bb7b35cp-12L,0x1.3502ed2e4cf93p-13L,0x0.0p+0L,0x0.0p+0L,0x1.bfa3d1bb7b35cp-12L,0x1.3502ed2e4cf93p-13L}}},
        {168,153,156,{{0x1.61bb57b31e5fap-16L,0x1.108ae2dc11fa9p-14L,0x0.0p+0L,0x0.0p+0L,0x1.61bb57b31e5fap-16L,0x1.108ae2dc11fa9p-14L}},{{0x1.37d0e102dd349p-12L,0x1.108ae2dc11fa9p-14L,0x0.0p+0L,0x0.0p+0L,0x1.37d0e102dd349p-12L,0x1.108ae2dc11fa9p-14L}}},
        {169,154,159,{{0x1.f99917b86f6fdp-18L,0x1.e0f03ede9cf5bp-16L,0x0.0p+0L,0x0.0p+0L,0x1.f99917b86f6fdp-18L,0x1.e0f03ede9cf5bp-16L}},{{0x1.ed30f104718a5p-14L,0x1.e0f03ede9cf5bp-16L,0x0.0p+0L,0x0.0p+0L,0x1.ed30f104718a5p-14L,0x1.e0f03ede9cf5bp-16L}}},
        {170,152,153,{{0x1.29b3712db7766p-15L,0x1.26a0c89aed4acp-13L,0x0.0p+0L,0x0.0p+0L,0x1.29b3712db7766p-15L,0x1.26a0c89aed4acp-13L}},{{0x1.28299a4c52305p-11L,0x1.26a0c89aed4acp-13L,0x0.0p+0L,0x0.0p+0L,0x1.28299a4c52305p-11L,0x1.26a0c89aed4acp-13L}}},
        {171,155,154,{{0x1.a6ea66fd4b5bap-17L,0x1.02517a796226ap-14L,0x0.0p+0L,0x0.0p+0L,0x1.a6ea66fd4b5bap-17L,0x1.02517a796226ap-14L}},{{0x1.d49ac5cb5cb7dp-13L,0x1.02517a796226ap-14L,0x0.0p+0L,0x0.0p+0L,0x1.d49ac5cb5cb7dp-13L,0x1.02517a796226ap-14L}}},
        {172,157,-1,{{0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L}},{{0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L}}},
        {173,158,-1,{{0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L}},{{0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L}}},
        {174,156,157,{{0x1.d1070595df3e7p-18L,0x1.165cc2cc8346cp-16L,0x0.0p+0L,0x0.0p+0L,0x1.d1070595df3e7p-18L,0x1.165cc2cc8346cp-16L}},{{0x1.6dbd9512a6ea8p-14L,0x1.165cc2cc8346cp-16L,0x0.0p+0L,0x0.0p+0L,0x1.6dbd9512a6ea8p-14L,0x1.165cc2cc8346cp-16L}}},
        {175,159,158,{{0x1.4e2058437d88fp-19L,0x1.eec641a0b881fp-18L,0x0.0p+0L,0x0.0p+0L,0x1.4e2058437d88fp-19L,0x1.eec641a0b881fp-18L}},{{0x1.2121637fc4454p-15L,0x1.eec641a0b881fp-18L,0x0.0p+0L,0x0.0p+0L,0x1.2121637fc4454p-15L,0x1.eec641a0b881fp-18L}}},
        {176,-1,0,{{0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L}},{{0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L}}},
        {177,-1,1,{{0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L}},{{0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L}}},
        {178,-1,2,{{0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L}},{{0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L}}},
        {179,-1,3,{{0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L}},{{0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L}}},
        {180,0,4,{{0x1.bc4a1482f5904p-17L,0x1.3488c7220e160p-13L,0x1.bc4a1482f5904p-17L,0x1.3488c7220e160p-13L,0x0.0p+0L,0x0.0p+0L}},{{0x1.8aec8402a163cp-12L,0x1.3488c7220e160p-13L,0x1.8aec8402a163cp-12L,0x1.3488c7220e160p-13L,0x0.0p+0L,0x0.0p+0L}}},
        {181,1,5,{{0x1.66cff8059bf24p-15L,0x1.c93c980b327d4p-13L,0x1.66cff8059bf24p-15L,0x1.c93c980b327d4p-13L,0x0.0p+0L,0x0.0p+0L}},{{0x1.968904439a545p-11L,0x1.c93c980b327d4p-13L,0x1.968904439a545p-11L,0x1.c93c980b327d4p-13L,0x0.0p+0L,0x0.0p+0L}}},
        {182,2,6,{{0x1.1d407d642c2fap-14L,0x1.ce0f9a73b7b4dp-14L,0x1.1d407d642c2fap-14L,0x1.ce0f9a73b7b4dp-14L,0x0.0p+0L,0x0.0p+0L}},{{0x1.7de853404c156p-11L,0x1.ce0f9a73b7b4dp-14L,0x1.7de853404c156p-11L,0x1.ce0f9a73b7b4dp-14L,0x0.0p+0L,0x0.0p+0L}}},
        {183,3,7,{{0x1.a9a6c13862d85p-14L,0x0.0p+0L,0x1.7143606103247p-14L,0x0.0p+0L,0x1.c31b06bafd9edp-17L,0x0.0p+0L}},{{0x1.9c8dee799ae16p-12L,0x0.0p+0L,0x1.5963b70aed154p-12L,0x0.0p+0L,0x1.0ca8ddbab7309p-14L,0x0.0p+0L}}},
        {184,4,8,{{0x1.026f045d8df3bp-15L,0x1.f1e6be2390147p-12L,0x1.026f045d8df3bp-15L,0x1.f1e6be2390147p-12L,0x0.0p+0L,0x0.0p+0L}},{{0x1.1b9a7ce3ae637p-10L,0x1.f1e6be2390147p-12L,0x1.1b9a7ce3ae637p-10L,0x1.f1e6be2390147p-12L,0x0.0p+0L,0x0.0p+0L}}},
        {185,5,9,{{0x1.b3bc831fed805p-14L,0x1.77e922d44787cp-11L,0x1.b3bc831fed805p-14L,0x1.77e922d44787cp-11L,0x0.0p+0L,0x0.0p+0L}},{{0x1.2388e2a5c9341p-9L,0x1.77e922d44787cp-11L,0x1.2388e2a5c9341p-9L,0x1.77e922d44787cp-11L,0x0.0p+0L,0x0.0p+0L}}},
        {186,6,10,{{0x1.6793285afd614p-13L,0x1.8d5eeba9a3122p-12L,0x1.571c59867dce8p-13L,0x1.706b07c000315p-12L,0x1.076ced47f92c8p-17L,0x1.cf3e3e9a2e0ccp-16L}},{{0x1.11355f74be609p-9L,0x1.8d5eeba9a3122p-12L,0x1.f4edebfaad79bp-10L,0x1.706b07c000315p-12L,0x1.6be697767a3b6p-13L,0x1.cf3e3e9a2e0ccp-16L}}},
        {187,7,11,{{0x1.17b309e07779ap-12L,0x1.112b9ed2691a0p-21L,0x0.0p+0L,0x0.0p+0L,0x1.17b309e07779ap-12L,0x1.112b9ed2691a0p-21L}},{{0x1.240c786fbcfe9p-10L,0x1.112b9ed2691a0p-21L,0x0.0p+0L,0x0.0p+0L,0x1.240c786fbcfe9p-10L,0x1.112b9ed2691a0p-21L}}},
        {188,8,12,{{0x1.3804152a7dbe1p-15L,0x1.a8b03940f2425p-11L,0x1.3804152a7dbe1p-15L,0x1.a8b03940f2425p-11L,0x0.0p+0L,0x0.0p+0L}},{{0x1.b15b00f3e787bp-10L,0x1.a8b03940f2425p-11L,0x1.b15b00f3e787bp-10L,0x1.a8b03940f2425p-11L,0x0.0p+0L,0x0.0p+0L}}},
        {189,9,13,{{0x1.1547e39e33d54p-13L,0x1.45f581566b94ep-10L,0x1.c98a3c5bd73efp-14L,0x1.1340cf23e8684p-10L,0x1.84162b8241ae3p-16L,0x1.95a5919419653p-13L}},{{0x1.bce1b4e161f45p-9L,0x1.45f581566b94ep-10L,0x1.964467f9d01edp-9L,0x1.1340cf23e8684p-10L,0x1.34ea673c8eac1p-12L,0x1.95a5919419653p-13L}}},
        {190,10,14,{{0x1.de3459c856a5bp-13L,0x1.660a54aa27911p-11L,0x0.0p+0L,0x0.0p+0L,0x1.de3459c856a5bp-13L,0x1.660a54aa27911p-11L}},{{0x1.9ff3ea54ca588p-9L,0x1.660a54aa27911p-11L,0x0.0p+0L,0x0.0p+0L,0x1.9ff3ea54ca588p-9L,0x1.660a54aa27911p-11L}}},
        {191,11,15,{{0x1.85dd69ce3ae93p-12L,0x1.8688885fc85e5p-18L,0x0.0p+0L,0x0.0p+0L,0x1.85dd69ce3ae93p-12L,0x1.8688885fc85e5p-18L}},{{0x1.b82a504d4e7c4p-10L,0x1.8688885fc85e5p-18L,0x0.0p+0L,0x0.0p+0L,0x1.b82a504d4e7c4p-10L,0x1.8688885fc85e5p-18L}}},
        {192,12,16,{{0x1.09ef32378cb8fp-15L,0x1.06643ed3b0822p-10L,0x1.cd5a8b19dbcb9p-16L,0x1.d7798134d4a65p-11L,0x1.1a0f6554f6993p-18L,0x1.aa77e39462ef9p-14L}},{{0x1.e26bacf26982cp-10L,0x1.06643ed3b0822p-10L,0x1.b83ef4ab88cf8p-10L,0x1.d7798134d4a65p-11L,0x1.5165c237059a0p-13L,0x1.aa77e39462ef9p-14L}}},
        {193,13,17,{{0x1.f8fb5608d258ep-14L,0x1.98b35514f668cp-10L,0x0.0p+0L,0x0.0p+0L,0x1.f8fb5608d258ep-14L,0x1.98b35514f668cp-10L}},{{0x1.ee9f2d430efcfp-9L,0x1.98b35514f668cp-10L,0x1.43bdbe284b411p-15L,0x0.0p+0L,0x1.e990364a6dcffp-9L,0x1.98b35514f668cp-10L}}},
        {194,14,18,{{0x1.cae19922a7573p-13L,0x1.d0022465486e1p-11L,0x0.0p+0L,0x0.0p+0L,0x1.cae19922a7573p-13L,0x1.d0022465486e1p-11L}},{{0x1.cd70f56ca215ap-9L,0x1.d0022465486e1p-11L,0x0.0p+0L,0x0.0p+0L,0x1.cd70f56ca215ap-9L,0x1.d0022465486e1p-11L}}},
        {195,15,19,{{0x1.8af27f82565e9p-12L,0x1.1f7ede0fb0d46p-16L,0x0.0p+0L,0x0.0p+0L,0x1.8af27f82565e9p-12L,0x1.1f7ede0fb0d46p-16L}},{{0x1.e3ae41defab80p-10L,0x1.1f7ede0fb0d46p-16L,0x0.0p+0L,0x0.0p+0L,0x1.e3ae41defab80p-10L,0x1.1f7ede0fb0d46p-16L}}},
        {196,16,20,{{0x1.502faaa2f741cp-16L,0x1.f2928ae2ab870p-11L,0x0.0p+0L,0x0.0p+0L,0x1.502faaa2f741cp-16L,0x1.f2928ae2ab870p-11L}},{{0x1.9f0b70afe8905p-10L,0x1.f2928ae2ab870p-11L,0x0.0p+0L,0x0.0p+0L,0x1.9f0b70afe8905p-10L,0x1.f2928ae2ab870p-11L}}},
        {197,17,21,{{0x1.5b9b72f187475p-14L,0x1.895dd10794624p-10L,0x0.0p+0L,0x0.0p+0L,0x1.5b9b72f187475p-14L,0x1.895dd10794624p-10L}},{{0x1.a906241611ce0p-9L,0x1.895dd10794624p-10L,0x0.0p+0L,0x0.0p+0L,0x1.a906241611ce0p-9L,0x1.895dd10794624p-10L}}},
        {198,18,22,{{0x1.508625e395bb2p-13L,0x1.cb9feb6f224e9p-11L,0x0.0p+0L,0x0.0p+0L,0x1.508625e395bb2p-13L,0x1.cb9feb6f224e9p-11L}},{{0x1.8bae441bfab52p-9L,0x1.cb9feb6f224e9p-11L,0x0.0p+0L,0x0.0p+0L,0x1.8bae441bfab52p-9L,0x1.cb9feb6f224e9p-11L}}},
        {199,19,23,{{0x1.345c965d8bbf1p-12L,0x1.d776a0ca500cfp-16L,0x0.0p+0L,0x0.0p+0L,0x1.345c965d8bbf1p-12L,0x1.d776a0ca500cfp-16L}},{{0x1.9b0cee7ed0f0ep-10L,0x1.d776a0ca500cfp-16L,0x0.0p+0L,0x0.0p+0L,0x1.9b0cee7ed0f0ep-10L,0x1.d776a0ca500cfp-16L}}},
        {200,20,24,{{0x1.24a641c1742bep-17L,0x1.5774e60f7be54p-11L,0x0.0p+0L,0x0.0p+0L,0x1.24a641c1742bep-17L,0x1.5774e60f7be54p-11L}},{{0x1.0422018f83b53p-10L,0x1.5774e60f7be54p-11L,0x0.0p+0L,0x0.0p+0L,0x1.0422018f83b53p-10L,0x1.5774e60f7be54p-11L}}},
        {201,21,25,{{0x1.531f013b10392p-15L,0x1.1228427caa019p-10L,0x0.0p+0L,0x0.0p+0L,0x1.531f013b10392p-15L,0x1.1228427caa019p-10L}},{{0x1.0a13ba5d43f04p-9L,0x1.1228427caa019p-10L,0x0.0p+0L,0x0.0p+0L,0x1.0a13ba5d43f04p-9L,0x1.1228427caa019p-10L}}},
        {202,22,26,{{0x1.62e3c5eb2fe34p-14L,0x1.4884c2df86370p-11L,0x0.0p+0L,0x0.0p+0L,0x1.62e3c5eb2fe34p-14L,0x1.4884c2df86370p-11L}},{{0x1.ee6c512c23119p-10L,0x1.4884c2df86370p-11L,0x0.0p+0L,0x0.0p+0L,0x1.ee6c512c23119p-10L,0x1.4884c2df86370p-11L}}},
        {203,23,27,{{0x1.5db2fe9d38a53p-13L,0x1.dea0f8a23dbdcp-16L,0x0.0p+0L,0x0.0p+0L,0x1.5db2fe9d38a53p-13L,0x1.dea0f8a23dbdcp-16L}},{{0x1.fd4d003473720p-11L,0x1.dea0f8a23dbdcp-16L,0x0.0p+0L,0x0.0p+0L,0x1.fd4d003473720p-11L,0x1.dea0f8a23dbdcp-16L}}},
        {204,24,28,{{0x1.f3867d762e622p-20L,0x1.f542833b09369p-13L,0x0.0p+0L,0x0.0p+0L,0x1.f3867d762e622p-20L,0x1.f542833b09369p-13L}},{{0x1.5ae481caf58b8p-12L,0x1.f542833b09369p-13L,0x0.0p+0L,0x0.0p+0L,0x1.5ae481caf58b8p-12L,0x1.f542833b09369p-13L}}},
        {205,25,29,{{0x1.53a83a62b1307p-17L,0x1.9457815a70f80p-12L,0x0.0p+0L,0x0.0p+0L,0x1.53a83a62b1307p-17L,0x1.9457815a70f80p-12L}},{{0x1.626c52970d123p-11L,0x1.9457815a70f80p-12L,0x0.0p+0L,0x0.0p+0L,0x1.626c52970d123p-11L,0x1.9457815a70f80p-12L}}},
        {206,26,30,{{0x1.87fb84847c753p-16L,0x1.ef70cd249ca59p-13L,0x0.0p+0L,0x0.0p+0L,0x1.87fb84847c753p-16L,0x1.ef70cd249ca59p-13L}},{{0x1.48a9f936e4818p-11L,0x1.ef70cd249ca59p-13L,0x0.0p+0L,0x0.0p+0L,0x1.48a9f936e4818p-11L,0x1.ef70cd249ca59p-13L}}},
        {207,27,31,{{0x1.a469f602fe923p-15L,0x1.d4e069d672641p-17L,0x0.0p+0L,0x0.0p+0L,0x1.a469f602fe923p-15L,0x1.d4e069d672641p-17L}},{{0x1.4fdb00a268809p-12L,0x1.d4e069d672641p-17L,0x0.0p+0L,0x0.0p+0L,0x1.4fdb00a268809p-12L,0x1.d4e069d672641p-17L}}},
        {208,28,-1,{{0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L}},{{0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L}}},
        {209,29,-1,{{0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L}},{{0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L}}},
        {210,30,-1,{{0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L}},{{0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L}}},
        {211,31,-1,{{0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L}},{{0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L}}},
        {212,-1,32,{{0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L}},{{0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L}}},
        {213,-1,33,{{0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L}},{{0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L}}},
        {214,35,48,{{0x1.51115bcf0159ep-14L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x1.51115bcf0159ep-14L,0x0.0p+0L}},{{0x1.5873f0f2c69f1p-14L,0x1.587f1bf49acd5p-15L,0x0.0p+0L,0x0.0p+0L,0x1.5873f0f2c69f1p-14L,0x1.587f1bf49acd5p-15L}}},
        {215,34,49,{{0x1.65113f961cb8cp-14L,0x1.d75cdf8b1d529p-19L,0x0.0p+0L,0x0.0p+0L,0x1.65113f961cb8cp-14L,0x1.d75cdf8b1d529p-19L}},{{0x1.65113f961cb8cp-14L,0x1.610a5fcf10bb3p-13L,0x0.0p+0L,0x0.0p+0L,0x1.65113f961cb8cp-14L,0x1.610a5fcf10bb3p-13L}}},
        {216,32,35,{{0x1.7c3a2fa429282p-16L,0x0.0p+0L,0x1.4fe702e8fe94dp-16L,0x0.0p+0L,0x1.629965d9549a5p-19L,0x0.0p+0L}},{{0x1.8403e77a6ad7ap-16L,0x1.8e9041d02340ep-17L,0x1.56aa6cba259f4p-16L,0x1.1baa96e20294cp-17L,0x1.6acbd60229c2dp-19L,0x1.cb96abb882b07p-19L}}},
        {217,33,34,{{0x1.88df9c248166bp-16L,0x1.3636ab469e5d0p-20L,0x0.0p+0L,0x0.0p+0L,0x1.88df9c248166bp-16L,0x1.3636ab469e5d0p-20L}},{{0x1.88df9c248166bp-16L,0x1.99c425b9f7c0ep-15L,0x0.0p+0L,0x0.0p+0L,0x1.88df9c248166bp-16L,0x1.99c425b9f7c0ep-15L}}},
        {218,-1,36,{{0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L}},{{0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L}}},
        {219,-1,37,{{0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L}},{{0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L}}},
        {220,39,52,{{0x1.63d66cca864fep-14L,0x1.0375bad3eaf81p-16L,0x0.0p+0L,0x0.0p+0L,0x1.63d66cca864fep-14L,0x1.0375bad3eaf81p-16L}},{{0x1.63d66cca864fep-14L,0x1.31bfff9b94896p-12L,0x0.0p+0L,0x0.0p+0L,0x1.63d66cca864fep-14L,0x1.31bfff9b94896p-12L}}},
        {221,38,53,{{0x1.41cf9e5555b65p-14L,0x1.fa3c320086680p-16L,0x0.0p+0L,0x0.0p+0L,0x1.41cf9e5555b65p-14L,0x1.fa3c320086680p-16L}},{{0x1.41cf9e5555b65p-14L,0x1.98d2c771a5fd8p-12L,0x0.0p+0L,0x0.0p+0L,0x1.41cf9e5555b65p-14L,0x1.98d2c771a5fd8p-12L}}},
        {222,36,39,{{0x1.814651dbe8c11p-16L,0x1.4691fc5b93adbp-18L,0x0.0p+0L,0x0.0p+0L,0x1.814651dbe8c11p-16L,0x1.4691fc5b93adbp-18L}},{{0x1.814651dbe8c11p-16L,0x1.6350f808bf905p-14L,0x0.0p+0L,0x0.0p+0L,0x1.814651dbe8c11p-16L,0x1.6350f808bf905p-14L}}},
        {223,37,38,{{0x1.5750ec1999335p-16L,0x1.3a917e3e78b5bp-17L,0x0.0p+0L,0x0.0p+0L,0x1.5750ec1999335p-16L,0x1.3a917e3e78b5bp-17L}},{{0x1.5750ec1999335p-16L,0x1.db7d00aeb8359p-14L,0x0.0p+0L,0x0.0p+0L,0x1.5750ec1999335p-16L,0x1.db7d00aeb8359p-14L}}},
        {224,-1,40,{{0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L}},{{0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L}}},
        {225,-1,41,{{0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L}},{{0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L}}},
        {226,43,56,{{0x1.0408aff21ee1dp-14L,0x1.6ebdc4b0c2c63p-15L,0x0.0p+0L,0x0.0p+0L,0x1.0408aff21ee1dp-14L,0x1.6ebdc4b0c2c63p-15L}},{{0x1.0408aff21ee1dp-14L,0x1.d2be185984656p-12L,0x0.0p+0L,0x0.0p+0L,0x1.0408aff21ee1dp-14L,0x1.d2be185984656p-12L}}},
        {227,42,57,{{0x1.6b234d14324a3p-15L,0x1.af558e6f9982ep-15L,0x0.0p+0L,0x0.0p+0L,0x1.6b234d14324a3p-15L,0x1.af558e6f9982ep-15L}},{{0x1.6b234d14324a3p-15L,0x1.caf1bc29a9661p-12L,0x0.0p+0L,0x0.0p+0L,0x1.6b234d14324a3p-15L,0x1.caf1bc29a9661p-12L}}},
        {228,40,43,{{0x1.10feaebf2b0d5p-16L,0x1.c4c0252b23875p-17L,0x0.0p+0L,0x0.0p+0L,0x1.10feaebf2b0d5p-16L,0x1.c4c0252b23875p-17L}},{{0x1.10feaebf2b0d5p-16L,0x1.0f9bfbcce49a4p-13L,0x0.0p+0L,0x0.0p+0L,0x1.10feaebf2b0d5p-16L,0x1.0f9bfbcce49a4p-13L}}},
        {229,41,42,{{0x1.75f9fffc3de05p-17L,0x1.091640cfe6f87p-16L,0x0.0p+0L,0x0.0p+0L,0x1.75f9fffc3de05p-17L,0x1.091640cfe6f87p-16L}},{{0x1.75f9fffc3de05p-17L,0x1.0b3c8cf41155ep-13L,0x0.0p+0L,0x0.0p+0L,0x1.75f9fffc3de05p-17L,0x1.0b3c8cf41155ep-13L}}},
        {230,-1,44,{{0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L}},{{0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L}}},
        {231,-1,45,{{0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L}},{{0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L}}},
        {232,47,60,{{0x1.8fcc8d560d021p-16L,0x1.8a277472b2621p-15L,0x0.0p+0L,0x0.0p+0L,0x1.8fcc8d560d021p-16L,0x1.8a277472b2621p-15L}},{{0x1.8fcc8d560d021p-16L,0x1.6a69d26799a1bp-12L,0x0.0p+0L,0x0.0p+0L,0x1.8fcc8d560d021p-16L,0x1.6a69d26799a1bp-12L}}},
        {233,46,61,{{0x1.cecd2ab315212p-18L,0x1.70f145a1b3c84p-16L,0x0.0p+0L,0x0.0p+0L,0x1.cecd2ab315212p-18L,0x1.70f145a1b3c84p-16L}},{{0x1.cecd2ab315212p-18L,0x1.2e3cd7f790692p-13L,0x0.0p+0L,0x0.0p+0L,0x1.cecd2ab315212p-18L,0x1.2e3cd7f790692p-13L}}},
        {234,44,47,{{0x1.91b6b3a0423fep-18L,0x1.e2f26a84a1a09p-17L,0x0.0p+0L,0x0.0p+0L,0x1.91b6b3a0423fep-18L,0x1.e2f26a84a1a09p-17L}},{{0x1.91b6b3a0423fep-18L,0x1.a64ea82f4e24bp-14L,0x0.0p+0L,0x0.0p+0L,0x1.91b6b3a0423fep-18L,0x1.a64ea82f4e24bp-14L}}},
        {235,45,46,{{0x1.c30b2056ad0fdp-20L,0x1.c30b2056ad0fdp-18L,0x0.0p+0L,0x0.0p+0L,0x1.c30b2056ad0fdp-20L,0x1.c30b2056ad0fdp-18L}},{{0x1.c30b2056ad0fdp-20L,0x1.6060b143b7346p-15L,0x0.0p+0L,0x0.0p+0L,0x1.c30b2056ad0fdp-20L,0x1.6060b143b7346p-15L}}},
        {236,51,64,{{0x1.003be9886a9d5p-12L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x1.003be9886a9d5p-12L,0x0.0p+0L}},{{0x1.0696968947e6dp-12L,0x1.f274c6bf927fep-14L,0x0.0p+0L,0x0.0p+0L,0x1.0696968947e6dp-12L,0x1.f274c6bf927fep-14L}}},
        {237,50,65,{{0x1.1c7f52c1c8c4bp-12L,0x1.009b918144936p-17L,0x0.0p+0L,0x0.0p+0L,0x1.1c7f52c1c8c4bp-12L,0x1.009b918144936p-17L}},{{0x1.1c7f52c1c8c4bp-12L,0x1.fba0e9211778cp-12L,0x0.0p+0L,0x0.0p+0L,0x1.1c7f52c1c8c4bp-12L,0x1.fba0e9211778cp-12L}}},
        {238,48,51,{{0x1.4ca2fdb250431p-13L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x1.4ca2fdb250431p-13L,0x0.0p+0L}},{{0x1.5467496296693p-13L,0x1.4b9850d69957bp-14L,0x0.0p+0L,0x0.0p+0L,0x1.5467496296693p-13L,0x1.4b9850d69957bp-14L}}},
        {239,49,50,{{0x1.68ea538de15a2p-13L,0x1.8b7889985c774p-18L,0x0.0p+0L,0x0.0p+0L,0x1.68ea538de15a2p-13L,0x1.8b7889985c774p-18L}},{{0x1.68ea538de15a2p-13L,0x1.52bf1697c0729p-12L,0x0.0p+0L,0x0.0p+0L,0x1.68ea538de15a2p-13L,0x1.52bf1697c0729p-12L}}},
        {240,55,68,{{0x1.23ef1599fd483p-12L,0x1.3a0b81005e5f0p-15L,0x0.0p+0L,0x0.0p+0L,0x1.23ef1599fd483p-12L,0x1.3a0b81005e5f0p-15L}},{{0x1.23ef1599fd483p-12L,0x1.b694d723f21d2p-11L,0x0.0p+0L,0x0.0p+0L,0x1.23ef1599fd483p-12L,0x1.b694d723f21d2p-11L}}},
        {241,54,69,{{0x1.0f0dd282046cfp-12L,0x1.3b4ac7f878954p-14L,0x0.0p+0L,0x0.0p+0L,0x1.0f0dd282046cfp-12L,0x1.3b4ac7f878954p-14L}},{{0x1.0f0dd282046cfp-12L,0x1.24c15848df8fep-10L,0x0.0p+0L,0x0.0p+0L,0x1.0f0dd282046cfp-12L,0x1.24c15848df8fep-10L}}},
        {242,52,55,{{0x1.6d235617dc306p-13L,0x1.c9a4ba80938ecp-16L,0x0.0p+0L,0x0.0p+0L,0x1.6d235617dc306p-13L,0x1.c9a4ba80938ecp-16L}},{{0x1.6d235617dc306p-13L,0x1.25042ed1c5798p-11L,0x0.0p+0L,0x0.0p+0L,0x1.6d235617dc306p-13L,0x1.25042ed1c5798p-11L}}},
        {243,53,54,{{0x1.4ec3510347395p-13L,0x1.c4a484a7c2e5ap-15L,0x0.0p+0L,0x0.0p+0L,0x1.4ec3510347395p-13L,0x1.c4a484a7c2e5ap-15L}},{{0x1.4ec3510347395p-13L,0x1.877c0086b907fp-11L,0x0.0p+0L,0x0.0p+0L,0x1.4ec3510347395p-13L,0x1.877c0086b907fp-11L}}},
        {244,59,72,{{0x1.c26f07baff79dp-13L,0x1.cfabc65c8dbb5p-14L,0x0.0p+0L,0x0.0p+0L,0x1.c26f07baff79dp-13L,0x1.cfabc65c8dbb5p-14L}},{{0x1.c26f07baff79dp-13L,0x1.4dcc243db2ca6p-10L,0x0.0p+0L,0x0.0p+0L,0x1.c26f07baff79dp-13L,0x1.4dcc243db2ca6p-10L}}},
        {245,58,73,{{0x1.44f18de37ed49p-13L,0x1.13523ef3a0bd5p-13L,0x0.0p+0L,0x0.0p+0L,0x1.44f18de37ed49p-13L,0x1.13523ef3a0bd5p-13L}},{{0x1.44f18de37ed49p-13L,0x1.47d4c60da3c05p-10L,0x0.0p+0L,0x0.0p+0L,0x1.44f18de37ed49p-13L,0x1.47d4c60da3c05p-10L}}},
        {246,56,59,{{0x1.12794e732a04bp-13L,0x1.4a45fc2b904ebp-14L,0x0.0p+0L,0x0.0p+0L,0x1.12794e732a04bp-13L,0x1.4a45fc2b904ebp-14L}},{{0x1.12794e732a04bp-13L,0x1.bea74b5b7a91cp-11L,0x0.0p+0L,0x0.0p+0L,0x1.12794e732a04bp-13L,0x1.bea74b5b7a91cp-11L}}},
        {247,57,58,{{0x1.85f20693f7046p-14L,0x1.86435573b3a14p-14L,0x0.0p+0L,0x0.0p+0L,0x1.85f20693f7046p-14L,0x1.86435573b3a14p-14L}},{{0x1.85f20693f7046p-14L,0x1.b6ed417f1d5f2p-11L,0x0.0p+0L,0x0.0p+0L,0x1.85f20693f7046p-14L,0x1.b6ed417f1d5f2p-11L}}},
        {248,63,76,{{0x1.748da9f035616p-14L,0x1.fab56fa3dca70p-14L,0x0.0p+0L,0x0.0p+0L,0x1.748da9f035616p-14L,0x1.fab56fa3dca70p-14L}},{{0x1.748da9f035616p-14L,0x1.02952ac7edcf7p-10L,0x0.0p+0L,0x0.0p+0L,0x1.748da9f035616p-14L,0x1.02952ac7edcf7p-10L}}},
        {249,62,77,{{0x1.c4a8bd6ab3cbcp-16L,0x1.dca8a0042d257p-15L,0x0.0p+0L,0x0.0p+0L,0x1.c4a8bd6ab3cbcp-16L,0x1.dca8a0042d257p-15L}},{{0x1.c4a8bd6ab3cbcp-16L,0x1.aed8f4dc58e51p-12L,0x0.0p+0L,0x0.0p+0L,0x1.c4a8bd6ab3cbcp-16L,0x1.aed8f4dc58e51p-12L}}},
        {250,60,63,{{0x1.b6a9fb2384eb3p-15L,0x1.65d4678cf94f6p-14L,0x0.0p+0L,0x0.0p+0L,0x1.b6a9fb2384eb3p-15L,0x1.65d4678cf94f6p-14L}},{{0x1.b6a9fb2384eb3p-15L,0x1.5a67e846d81fcp-11L,0x0.0p+0L,0x0.0p+0L,0x1.b6a9fb2384eb3p-15L,0x1.5a67e846d81fcp-11L}}},
        {251,61,62,{{0x1.049599f69347fp-16L,0x1.4fbefa95025f2p-15L,0x0.0p+0L,0x0.0p+0L,0x1.049599f69347fp-16L,0x1.4fbefa95025f2p-15L}},{{0x1.049599f69347fp-16L,0x1.20bc74b39b645p-12L,0x0.0p+0L,0x0.0p+0L,0x1.049599f69347fp-16L,0x1.20bc74b39b645p-12L}}},
        {252,67,80,{{0x1.9dca9065cf883p-12L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x1.9dca9065cf883p-12L,0x0.0p+0L}},{{0x1.a942f1b785a75p-12L,0x1.7f8d5445f28aep-13L,0x0.0p+0L,0x0.0p+0L,0x1.a942f1b785a75p-12L,0x1.7f8d5445f28aep-13L}}},
        {253,66,81,{{0x1.dfe2a74d287d3p-12L,0x1.1c268d7da22d3p-17L,0x0.0p+0L,0x0.0p+0L,0x1.dfe2a74d287d3p-12L,0x1.1c268d7da22d3p-17L}},{{0x1.dfe2a74d287d3p-12L,0x1.84422a58335a2p-11L,0x0.0p+0L,0x0.0p+0L,0x1.dfe2a74d287d3p-12L,0x1.84422a58335a2p-11L}}},
        {254,64,67,{{0x1.560f1a23f04b0p-12L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x1.560f1a23f04b0p-12L,0x0.0p+0L}},{{0x1.5f0a4102ef854p-12L,0x1.44bc879413f51p-13L,0x0.0p+0L,0x0.0p+0L,0x1.5f0a4102ef854p-12L,0x1.44bc879413f51p-13L}}},
        {255,65,66,{{0x1.844eb261f3131p-12L,0x1.1d7645da245b1p-17L,0x0.0p+0L,0x0.0p+0L,0x1.844eb261f3131p-12L,0x1.1d7645da245b1p-17L}},{{0x1.844eb261f3131p-12L,0x1.49b4c0c70cdb3p-11L,0x0.0p+0L,0x0.0p+0L,0x1.844eb261f3131p-12L,0x1.49b4c0c70cdb3p-11L}}},
        {256,71,84,{{0x1.f95b87552ddfcp-12L,0x1.8e19a757b1dc0p-15L,0x0.0p+0L,0x0.0p+0L,0x1.f95b87552ddfcp-12L,0x1.8e19a757b1dc0p-15L}},{{0x1.f95b87552ddfcp-12L,0x1.4eaec58554e88p-10L,0x0.0p+0L,0x0.0p+0L,0x1.f95b87552ddfcp-12L,0x1.4eaec58554e88p-10L}}},
        {257,70,85,{{0x1.dff7169aa73c9p-12L,0x1.9d8a3fa84b1ddp-14L,0x0.0p+0L,0x0.0p+0L,0x1.dff7169aa73c9p-12L,0x1.9d8a3fa84b1ddp-14L}},{{0x1.dff7169aa73c9p-12L,0x1.be226d6793d90p-10L,0x0.0p+0L,0x0.0p+0L,0x1.dff7169aa73c9p-12L,0x1.be226d6793d90p-10L}}},
        {258,68,71,{{0x1.93cee4532eba6p-12L,0x1.7428ff38ed43fp-15L,0x0.0p+0L,0x0.0p+0L,0x1.93cee4532eba6p-12L,0x1.7428ff38ed43fp-15L}},{{0x1.93cee4532eba6p-12L,0x1.1c87f228f665dp-10L,0x0.0p+0L,0x0.0p+0L,0x1.93cee4532eba6p-12L,0x1.1c87f228f665dp-10L}}},
        {259,69,70,{{0x1.7b5a4d14d287fp-12L,0x1.7bc61db3bdf27p-14L,0x0.0p+0L,0x0.0p+0L,0x1.7b5a4d14d287fp-12L,0x1.7bc61db3bdf27p-14L}},{{0x1.7b5a4d14d287fp-12L,0x1.7b90337a48f7ep-10L,0x0.0p+0L,0x0.0p+0L,0x1.7b5a4d14d287fp-12L,0x1.7b90337a48f7ep-10L}}},
        {260,75,88,{{0x1.983eab11dc27fp-12L,0x1.356201d763f07p-13L,0x0.0p+0L,0x0.0p+0L,0x1.983eab11dc27fp-12L,0x1.356201d763f07p-13L}},{{0x1.983eab11dc27fp-12L,0x1.fc0d9caf9a83dp-10L,0x0.0p+0L,0x0.0p+0L,0x1.983eab11dc27fp-12L,0x1.fc0d9caf9a83dp-10L}}},
        {261,74,89,{{0x1.2e72165bd1abap-12L,0x1.7388f8b9f0effp-13L,0x0.0p+0L,0x0.0p+0L,0x1.2e72165bd1abap-12L,0x1.7388f8b9f0effp-13L}},{{0x1.2e72165bd1abap-12L,0x1.f2696be4aef15p-10L,0x0.0p+0L,0x0.0p+0L,0x1.2e72165bd1abap-12L,0x1.f2696be4aef15p-10L}}},
        {262,72,75,{{0x1.3f11ba1aab3cap-12L,0x1.1994a92fc763dp-13L,0x0.0p+0L,0x0.0p+0L,0x1.3f11ba1aab3cap-12L,0x1.1994a92fc763dp-13L}},{{0x1.3f11ba1aab3cap-12L,0x1.b0816670adc5dp-10L,0x0.0p+0L,0x0.0p+0L,0x1.3f11ba1aab3cap-12L,0x1.b0816670adc5dp-10L}}},
        {263,73,74,{{0x1.d2d1488abc0d2p-13L,0x1.503337d2cae9dp-13L,0x0.0p+0L,0x0.0p+0L,0x1.d2d1488abc0d2p-13L,0x1.503337d2cae9dp-13L}},{{0x1.d2d1488abc0d2p-13L,0x1.a8887714dc712p-10L,0x0.0p+0L,0x0.0p+0L,0x1.d2d1488abc0d2p-13L,0x1.a8887714dc712p-10L}}},
        {264,79,92,{{0x1.66056af89963cp-13L,0x1.58a47e978585ep-13L,0x0.0p+0L,0x0.0p+0L,0x1.66056af89963cp-13L,0x1.58a47e978585ep-13L}},{{0x1.66056af89963cp-13L,0x1.88b5402619c80p-10L,0x0.0p+0L,0x0.0p+0L,0x1.66056af89963cp-13L,0x1.88b5402619c80p-10L}}},
        {265,78,93,{{0x1.c33f6bd4dc366p-15L,0x1.460d301dbaa3cp-14L,0x0.0p+0L,0x0.0p+0L,0x1.c33f6bd4dc366p-15L,0x1.460d301dbaa3cp-14L}},{{0x1.c33f6bd4dc366p-15L,0x1.46d7d4a6b1596p-11L,0x0.0p+0L,0x0.0p+0L,0x1.c33f6bd4dc366p-15L,0x1.46d7d4a6b1596p-11L}}},
        {266,76,79,{{0x1.10296d88e1703p-13L,0x1.369342dd6c860p-13L,0x0.0p+0L,0x0.0p+0L,0x1.10296d88e1703p-13L,0x1.369342dd6c860p-13L}},{{0x1.10296d88e1703p-13L,0x1.4eacdd1796a88p-10L,0x0.0p+0L,0x0.0p+0L,0x1.10296d88e1703p-13L,0x1.4eacdd1796a88p-10L}}},
        {267,77,78,{{0x1.51333baaa10aap-15L,0x1.24f5c47fe111dp-14L,0x0.0p+0L,0x0.0p+0L,0x1.51333baaa10aap-15L,0x1.24f5c47fe111dp-14L}},{{0x1.51333baaa10aap-15L,0x1.16ad83d2d57d1p-11L,0x0.0p+0L,0x0.0p+0L,0x1.51333baaa10aap-15L,0x1.16ad83d2d57d1p-11L}}},
        {268,83,96,{{0x1.e613edae9d90cp-12L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x1.e613edae9d90cp-12L,0x0.0p+0L}},{{0x1.f4fb51a3bad30p-12L,0x1.adeb87c831086p-13L,0x0.0p+0L,0x0.0p+0L,0x1.f4fb51a3bad30p-12L,0x1.adeb87c831086p-13L}}},
        {269,82,97,{{0x1.25818805a8d63p-11L,0x1.aeb05e4830c37p-18L,0x0.0p+0L,0x0.0p+0L,0x1.25818805a8d63p-11L,0x1.aeb05e4830c37p-18L}},{{0x1.25818805a8d63p-11L,0x1.b0aafd1c2b5f8p-11L,0x0.0p+0L,0x0.0p+0L,0x1.25818805a8d63p-11L,0x1.b0aafd1c2b5f8p-11L}}},
        {270,80,83,{{0x1.cf9b10f30294fp-12L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x1.cf9b10f30294fp-12L,0x0.0p+0L}},{{0x1.dd2349c59d8ddp-12L,0x1.a3b461ae3b470p-13L,0x0.0p+0L,0x0.0p+0L,0x1.dd2349c59d8ddp-12L,0x1.a3b461ae3b470p-13L}}},
        {271,81,82,{{0x1.126ceba8c3363p-11L,0x1.020baddf5600ap-17L,0x0.0p+0L,0x0.0p+0L,0x1.126ceba8c3363p-11L,0x1.020baddf5600ap-17L}},{{0x1.126ceba8c3363p-11L,0x1.a79bac2de0ee9p-11L,0x0.0p+0L,0x0.0p+0L,0x1.126ceba8c3363p-11L,0x1.a79bac2de0ee9p-11L}}},
        {272,87,100,{{0x1.3c4f5ebda8215p-11L,0x1.68c0ef0ea243ep-15L,0x0.0p+0L,0x0.0p+0L,0x1.3c4f5ebda8215p-11L,0x1.68c0ef0ea243ep-15L}},{{0x1.3c4f5ebda8215p-11L,0x1.74266bf0a07c0p-10L,0x0.0p+0L,0x0.0p+0L,0x1.3c4f5ebda8215p-11L,0x1.74266bf0a07c0p-10L}}},
        {273,86,101,{{0x1.32654be091790p-11L,0x1.866cbe408b4e6p-14L,0x0.0p+0L,0x0.0p+0L,0x1.32654be091790p-11L,0x1.866cbe408b4e6p-14L}},{{0x1.32654be091790p-11L,0x1.ef5e937f5ca2bp-10L,0x0.0p+0L,0x0.0p+0L,0x1.32654be091790p-11L,0x1.ef5e937f5ca2bp-10L}}},
        {274,84,87,{{0x1.2471bcec992eep-11L,0x1.88b2b88376c41p-15L,0x0.0p+0L,0x0.0p+0L,0x1.2471bcec992eep-11L,0x1.88b2b88376c41p-15L}},{{0x1.2471bcec992eep-11L,0x1.6cc05a794b207p-10L,0x0.0p+0L,0x0.0p+0L,0x1.2471bcec992eep-11L,0x1.6cc05a794b207p-10L}}},
        {275,85,86,{{0x1.1897f9c23ff8fp-11L,0x1.9ff7935cd68f5p-14L,0x0.0p+0L,0x0.0p+0L,0x1.1897f9c23ff8fp-11L,0x1.9ff7935cd68f5p-14L}},{{0x1.1897f9c23ff8fp-11L,0x1.e5dd48d579d49p-10L,0x0.0p+0L,0x0.0p+0L,0x1.1897f9c23ff8fp-11L,0x1.e5dd48d579d49p-10L}}},
        {276,91,104,{{0x1.09d9b993603d0p-11L,0x1.2a124798d83b1p-13L,0x0.0p+0L,0x0.0p+0L,0x1.09d9b993603d0p-11L,0x1.2a124798d83b1p-13L}},{{0x1.09d9b993603d0p-11L,0x1.19bb0aade325ap-9L,0x0.0p+0L,0x0.0p+0L,0x1.09d9b993603d0p-11L,0x1.19bb0aade325ap-9L}}},
        {277,90,105,{{0x1.92b32d7d624a7p-12L,0x1.6a9a23a279b9dp-13L,0x0.0p+0L,0x0.0p+0L,0x1.92b32d7d624a7p-12L,0x1.6a9a23a279b9dp-13L}},{{0x1.92b32d7d624a7p-12L,0x1.1416fa1666dc0p-9L,0x0.0p+0L,0x0.0p+0L,0x1.92b32d7d624a7p-12L,0x1.1416fa1666dc0p-9L}}},
        {278,88,91,{{0x1.e24a8572acc89p-12L,0x1.3a3d3ba10c015p-13L,0x0.0p+0L,0x0.0p+0L,0x1.e24a8572acc89p-12L,0x1.3a3d3ba10c015p-13L}},{{0x1.e24a8572acc89p-12L,0x1.147c11af3fab4p-9L,0x0.0p+0L,0x0.0p+0L,0x1.e24a8572acc89p-12L,0x1.147c11af3fab4p-9L}}},
        {279,89,90,{{0x1.69721087688e0p-12L,0x1.7bb82387452fep-13L,0x0.0p+0L,0x0.0p+0L,0x1.69721087688e0p-12L,0x1.7bb82387452fep-13L}},{{0x1.69721087688e0p-12L,0x1.0f17639a81bdfp-9L,0x0.0p+0L,0x0.0p+0L,0x1.69721087688e0p-12L,0x1.0f17639a81bdfp-9L}}},
        {280,95,108,{{0x1.e9297d897d987p-13L,0x1.537816e5122c0p-13L,0x0.0p+0L,0x0.0p+0L,0x1.e9297d897d987p-13L,0x1.537816e5122c0p-13L}},{{0x1.e9297d897d987p-13L,0x1.b2a11948a5db7p-10L,0x0.0p+0L,0x0.0p+0L,0x1.e9297d897d987p-13L,0x1.b2a11948a5db7p-10L}}},
        {281,94,109,{{0x1.3d3daf937b90fp-14L,0x1.43413f5e06680p-14L,0x0.0p+0L,0x0.0p+0L,0x1.3d3daf937b90fp-14L,0x1.43413f5e06680p-14L}},{{0x1.3d3daf937b90fp-14L,0x1.696642c6736e0p-11L,0x0.0p+0L,0x0.0p+0L,0x1.3d3daf937b90fp-14L,0x1.696642c6736e0p-11L}}},
        {282,92,95,{{0x1.b1b42fc7e7467p-13L,0x1.61cd9755f0e9cp-13L,0x0.0p+0L,0x0.0p+0L,0x1.b1b42fc7e7467p-13L,0x1.61cd9755f0e9cp-13L}},{{0x1.b1b42fc7e7467p-13L,0x1.aaf9ad059ad68p-10L,0x0.0p+0L,0x0.0p+0L,0x1.b1b42fc7e7467p-13L,0x1.aaf9ad059ad68p-10L}}},
        {283,93,94,{{0x1.157f586553793p-14L,0x1.4fc59b9953154p-14L,0x0.0p+0L,0x0.0p+0L,0x1.157f586553793p-14L,0x1.4fc59b9953154p-14L}},{{0x1.157f586553793p-14L,0x1.63325caaa8494p-11L,0x0.0p+0L,0x0.0p+0L,0x1.157f586553793p-14L,0x1.63325caaa8494p-11L}}},
        {284,99,112,{{0x1.b8a4776721701p-12L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x1.b8a4776721701p-12L,0x0.0p+0L}},{{0x1.c77423f367e7ap-12L,0x1.745550ba31f69p-13L,0x0.0p+0L,0x0.0p+0L,0x1.c77423f367e7ap-12L,0x1.745550ba31f69p-13L}}},
        {285,98,113,{{0x1.144ebacda9833p-11L,0x1.cbe72fca3c15cp-19L,0x0.0p+0L,0x0.0p+0L,0x1.144ebacda9833p-11L,0x1.cbe72fca3c15cp-19L}},{{0x1.144ebacda9833p-11L,0x1.749cc3ed9bdd2p-11L,0x0.0p+0L,0x0.0p+0L,0x1.144ebacda9833p-11L,0x1.749cc3ed9bdd2p-11L}}},
        {286,96,99,{{0x1.de684e0d04887p-12L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x1.de684e0d04887p-12L,0x0.0p+0L}},{{0x1.edc8236ce6eebp-12L,0x1.9d85720cc8d4cp-13L,0x0.0p+0L,0x0.0p+0L,0x1.edc8236ce6eebp-12L,0x1.9d85720cc8d4cp-13L}}},
        {287,97,98,{{0x1.26796de1ce981p-11L,0x1.49ed1845f8b5ap-18L,0x0.0p+0L,0x0.0p+0L,0x1.26796de1ce981p-11L,0x1.49ed1845f8b5ap-18L}},{{0x1.26796de1ce981p-11L,0x1.9efc6566f09a7p-11L,0x0.0p+0L,0x0.0p+0L,0x1.26796de1ce981p-11L,0x1.9efc6566f09a7p-11L}}},
        {288,103,116,{{0x1.300272a65e48ep-11L,0x1.eda932b28484bp-16L,0x0.0p+0L,0x0.0p+0L,0x1.300272a65e48ep-11L,0x1.eda932b28484bp-16L}},{{0x1.300272a65e48ep-11L,0x1.3fd31d1b6f08ap-10L,0x0.0p+0L,0x0.0p+0L,0x1.300272a65e48ep-11L,0x1.3fd31d1b6f08ap-10L}}},
        {289,102,117,{{0x1.2b9d942db2199p-11L,0x1.18d9679759cf5p-14L,0x0.0p+0L,0x0.0p+0L,0x1.2b9d942db2199p-11L,0x1.18d9679759cf5p-14L}},{{0x1.2b9d942db2199p-11L,0x1.a92359eb6ffdep-10L,0x0.0p+0L,0x0.0p+0L,0x1.2b9d942db2199p-11L,0x1.a92359eb6ffdep-10L}}},
        {290,100,103,{{0x1.40bf292aeacccp-11L,0x1.3569788a2dd41p-15L,0x0.0p+0L,0x0.0p+0L,0x1.40bf292aeacccp-11L,0x1.3569788a2dd41p-15L}},{{0x1.40bf292aeacccp-11L,0x1.649065c48526dp-10L,0x0.0p+0L,0x0.0p+0L,0x1.40bf292aeacccp-11L,0x1.649065c48526dp-10L}}},
        {291,101,102,{{0x1.3979f56e48b41p-11L,0x1.56e6e6f79570ap-14L,0x0.0p+0L,0x0.0p+0L,0x1.3979f56e48b41p-11L,0x1.56e6e6f79570ap-14L}},{{0x1.3979f56e48b41p-11L,0x1.da4b9085b2f2dp-10L,0x0.0p+0L,0x0.0p+0L,0x1.3979f56e48b41p-11L,0x1.da4b9085b2f2dp-10L}}},
        {292,107,120,{{0x1.0874aad176631p-11L,0x1.b73e81313aac5p-14L,0x0.0p+0L,0x0.0p+0L,0x1.0874aad176631p-11L,0x1.b73e81313aac5p-14L}},{{0x1.0874aad176631p-11L,0x1.e30998ee62b23p-10L,0x0.0p+0L,0x0.0p+0L,0x1.0874aad176631p-11L,0x1.e30998ee62b23p-10L}}},
        {293,106,121,{{0x1.981bbd2f82f1cp-12L,0x1.0f37b51384b5dp-13L,0x0.0p+0L,0x0.0p+0L,0x1.981bbd2f82f1cp-12L,0x1.0f37b51384b5dp-13L}},{{0x1.981bbd2f82f1cp-12L,0x1.d8e312b74e21dp-10L,0x0.0p+0L,0x0.0p+0L,0x1.981bbd2f82f1cp-12L,0x1.d8e312b74e21dp-10L}}},
        {294,104,107,{{0x1.126adfe89f467p-11L,0x1.08cfb0d4ed2e3p-13L,0x0.0p+0L,0x0.0p+0L,0x1.126adfe89f467p-11L,0x1.08cfb0d4ed2e3p-13L}},{{0x1.126adfe89f467p-11L,0x1.0d97ce5a4c5c3p-9L,0x0.0p+0L,0x0.0p+0L,0x1.126adfe89f467p-11L,0x1.0d97ce5a4c5c3p-9L}}},
        {295,105,106,{{0x1.a3b8a621b9571p-12L,0x1.447a450902887p-13L,0x0.0p+0L,0x0.0p+0L,0x1.a3b8a621b9571p-12L,0x1.447a450902887p-13L}},{{0x1.a3b8a621b9571p-12L,0x1.080f4ad438d24p-9L,0x0.0p+0L,0x0.0p+0L,0x1.a3b8a621b9571p-12L,0x1.080f4ad438d24p-9L}}},
        {296,111,124,{{0x1.fa620c7ebe016p-13L,0x1.00a12272cd0b3p-13L,0x0.0p+0L,0x0.0p+0L,0x1.fa620c7ebe016p-13L,0x1.00a12272cd0b3p-13L}},{{0x1.fa620c7ebe016p-13L,0x1.73db7f650ed7ap-10L,0x0.0p+0L,0x0.0p+0L,0x1.fa620c7ebe016p-13L,0x1.73db7f650ed7ap-10L}}},
        {297,110,125,{{0x1.500fe80a39620p-14L,0x1.ec696131e9cfbp-15L,0x0.0p+0L,0x0.0p+0L,0x1.500fe80a39620p-14L,0x1.ec696131e9cfbp-15L}},{{0x1.500fe80a39620p-14L,0x1.34ef12a492972p-11L,0x0.0p+0L,0x0.0p+0L,0x1.500fe80a39620p-14L,0x1.34ef12a492972p-11L}}},
        {298,108,111,{{0x1.01c5122b04e1dp-12L,0x1.31564fbc9b71bp-13L,0x0.0p+0L,0x0.0p+0L,0x1.01c5122b04e1dp-12L,0x1.31564fbc9b71bp-13L}},{{0x1.01c5122b04e1dp-12L,0x1.9f7ce190cb27fp-10L,0x0.0p+0L,0x0.0p+0L,0x1.01c5122b04e1dp-12L,0x1.9f7ce190cb27fp-10L}}},
        {299,109,110,{{0x1.526c6d9fc10cfp-14L,0x1.23cdd17627c34p-14L,0x0.0p+0L,0x0.0p+0L,0x1.526c6d9fc10cfp-14L,0x1.23cdd17627c34p-14L}},{{0x1.526c6d9fc10cfp-14L,0x1.595484c25c802p-11L,0x0.0p+0L,0x0.0p+0L,0x1.526c6d9fc10cfp-14L,0x1.595484c25c802p-11L}}},
        {300,115,128,{{0x1.229a81dc06a62p-12L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x1.229a81dc06a62p-12L,0x0.0p+0L}},{{0x1.2d3abe51767cfp-12L,0x1.d5b78ecd0ca95p-14L,0x0.0p+0L,0x0.0p+0L,0x1.2d3abe51767cfp-12L,0x1.d5b78ecd0ca95p-14L}}},
        {301,114,129,{{0x1.798a5d7b0dce3p-12L,0x1.337b4c556833fp-20L,0x0.0p+0L,0x0.0p+0L,0x1.798a5d7b0dce3p-12L,0x1.337b4c556833fp-20L}},{{0x1.798a5d7b0dce3p-12L,0x1.d3860e8e64cfap-12L,0x0.0p+0L,0x0.0p+0L,0x1.798a5d7b0dce3p-12L,0x1.d3860e8e64cfap-12L}}},
        {302,112,115,{{0x1.77e7154d9a30fp-12L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x1.77e7154d9a30fp-12L,0x0.0p+0L}},{{0x1.8517c3584856ep-12L,0x1.3697992500889p-13L,0x0.0p+0L,0x0.0p+0L,0x1.8517c3584856ep-12L,0x1.3697992500889p-13L}}},
        {303,113,114,{{0x1.dff4cedcd7f5ap-12L,0x1.1e410125b3ad6p-19L,0x0.0p+0L,0x0.0p+0L,0x1.dff4cedcd7f5ap-12L,0x1.1e410125b3ad6p-19L}},{{0x1.dff4cedcd7f5ap-12L,0x1.35f99bea38387p-11L,0x0.0p+0L,0x0.0p+0L,0x1.dff4cedcd7f5ap-12L,0x1.35f99bea38387p-11L}}},
        {304,119,132,{{0x1.a739104ec74e2p-12L,0x1.dde84781aa136p-17L,0x0.0p+0L,0x0.0p+0L,0x1.a739104ec74e2p-12L,0x1.dde84781aa136p-17L}},{{0x1.a739104ec74e2p-12L,0x1.907c64187c1f5p-11L,0x0.0p+0L,0x0.0p+0L,0x1.a739104ec74e2p-12L,0x1.907c64187c1f5p-11L}}},
        {305,118,133,{{0x1.a786d01aba41bp-12L,0x1.2164db3c033f0p-15L,0x0.0p+0L,0x0.0p+0L,0x1.a786d01aba41bp-12L,0x1.2164db3c033f0p-15L}},{{0x1.a786d01aba41bp-12L,0x1.09d53b9f467f0p-10L,0x0.0p+0L,0x0.0p+0L,0x1.a786d01aba41bp-12L,0x1.09d53b9f467f0p-10L}}},
        {306,116,119,{{0x1.0a93ffab5fe0fp-11L,0x1.6a38b463c95c2p-16L,0x0.0p+0L,0x0.0p+0L,0x1.0a93ffab5fe0fp-11L,0x1.6a38b463c95c2p-16L}},{{0x1.0a93ffab5fe0fp-11L,0x1.09caa9085f49ap-10L,0x0.0p+0L,0x0.0p+0L,0x1.0a93ffab5fe0fp-11L,0x1.09caa9085f49ap-10L}}},
        {307,117,118,{{0x1.08ccddbe35d79p-11L,0x1.a870daf773739p-15L,0x0.0p+0L,0x0.0p+0L,0x1.08ccddbe35d79p-11L,0x1.a870daf773739p-15L}},{{0x1.08ccddbe35d79p-11L,0x1.611470ecfbda9p-10L,0x0.0p+0L,0x0.0p+0L,0x1.08ccddbe35d79p-11L,0x1.611470ecfbda9p-10L}}},
        {308,123,136,{{0x1.7b6cbd8e805a9p-12L,0x1.d1cdb479ee043p-15L,0x0.0p+0L,0x0.0p+0L,0x1.7b6cbd8e805a9p-12L,0x1.d1cdb479ee043p-15L}},{{0x1.7b6cbd8e805a9p-12L,0x1.2db74bc1a62ddp-10L,0x0.0p+0L,0x0.0p+0L,0x1.7b6cbd8e805a9p-12L,0x1.2db74bc1a62ddp-10L}}},
        {309,122,137,{{0x1.2972d39dd9683p-12L,0x1.24be6fd079622p-14L,0x0.0p+0L,0x0.0p+0L,0x1.2972d39dd9683p-12L,0x1.24be6fd079622p-14L}},{{0x1.2972d39dd9683p-12L,0x1.27176e7d663e6p-10L,0x0.0p+0L,0x0.0p+0L,0x1.2972d39dd9683p-12L,0x1.27176e7d663e6p-10L}}},
        {310,120,123,{{0x1.d70f16b952d33p-12L,0x1.507cb71b3a0aap-14L,0x0.0p+0L,0x0.0p+0L,0x1.d70f16b952d33p-12L,0x1.507cb71b3a0aap-14L}},{{0x1.d70f16b952d33p-12L,0x1.90f33c7318beep-10L,0x0.0p+0L,0x0.0p+0L,0x1.d70f16b952d33p-12L,0x1.90f33c7318beep-10L}}},
        {311,121,122,{{0x1.6e79fda54a1bfp-12L,0x1.a31251867794ep-14L,0x0.0p+0L,0x0.0p+0L,0x1.6e79fda54a1bfp-12L,0x1.a31251867794ep-14L}},{{0x1.6e79fda54a1bfp-12L,0x1.885557ca463f7p-10L,0x0.0p+0L,0x0.0p+0L,0x1.6e79fda54a1bfp-12L,0x1.885557ca463f7p-10L}}},
        {312,127,140,{{0x1.77b3b95e2322ep-13L,0x1.18708616ad4cep-14L,0x0.0p+0L,0x0.0p+0L,0x1.77b3b95e2322ep-13L,0x1.18708616ad4cep-14L}},{{0x1.77b3b95e2322ep-13L,0x1.cfab1b9a7db10p-11L,0x0.0p+0L,0x0.0p+0L,0x1.77b3b95e2322ep-13L,0x1.cfab1b9a7db10p-11L}}},
        {313,126,141,{{0x1.fc2cf65ad40cbp-15L,0x1.0f60b893b9313p-15L,0x0.0p+0L,0x0.0p+0L,0x1.fc2cf65ad40cbp-15L,0x1.0f60b893b9313p-15L}},{{0x1.fc2cf65ad40cbp-15L,0x1.80e3f2a2b5d01p-12L,0x0.0p+0L,0x0.0p+0L,0x1.fc2cf65ad40cbp-15L,0x1.80e3f2a2b5d01p-12L}}},
        {314,124,127,{{0x1.caf6ba80cb46fp-13L,0x1.8ee47ae14a0c1p-14L,0x0.0p+0L,0x0.0p+0L,0x1.caf6ba80cb46fp-13L,0x1.8ee47ae14a0c1p-14L}},{{0x1.caf6ba80cb46fp-13L,0x1.345ece507654ep-10L,0x0.0p+0L,0x0.0p+0L,0x1.caf6ba80cb46fp-13L,0x1.345ece507654ep-10L}}},
        {315,125,126,{{0x1.339e6d7bbd32fp-14L,0x1.804893286df31p-15L,0x0.0p+0L,0x0.0p+0L,0x1.339e6d7bbd32fp-14L,0x1.804893286df31p-15L}},{{0x1.339e6d7bbd32fp-14L,0x1.00151159f49c7p-11L,0x0.0p+0L,0x0.0p+0L,0x1.339e6d7bbd32fp-14L,0x1.00151159f49c7p-11L}}},
        {316,131,144,{{0x1.973d28245fec7p-14L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x1.973d28245fec7p-14L,0x0.0p+0L}},{{0x1.a75447e500344p-14L,0x1.3b2053a011c62p-15L,0x0.0p+0L,0x0.0p+0L,0x1.a75447e500344p-14L,0x1.3b2053a011c62p-15L}}},
        {317,130,145,{{0x1.116da70138d3dp-13L,0x1.45dc98b520aa9p-23L,0x0.0p+0L,0x0.0p+0L,0x1.116da70138d3dp-13L,0x1.45dc98b520aa9p-23L}},{{0x1.116da70138d3dp-13L,0x1.38036247ae573p-13L,0x0.0p+0L,0x0.0p+0L,0x1.116da70138d3dp-13L,0x1.38036247ae573p-13L}}},
        {318,128,131,{{0x1.855c1c184d480p-13L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x1.855c1c184d480p-13L,0x0.0p+0L}},{{0x1.942bb358e1f4fp-13L,0x1.33de29088db75p-14L,0x0.0p+0L,0x0.0p+0L,0x1.942bb358e1f4fp-13L,0x1.33de29088db75p-14L}}},
        {319,129,130,{{0x1.01347165e8fe0p-12L,0x1.0b0b584b4cf18p-21L,0x0.0p+0L,0x0.0p+0L,0x1.01347165e8fe0p-12L,0x1.0b0b584b4cf18p-21L}},{{0x1.01347165e8fe0p-12L,0x1.319ed8bfc89dcp-12L,0x0.0p+0L,0x0.0p+0L,0x1.01347165e8fe0p-12L,0x1.319ed8bfc89dcp-12L}}},
        {320,135,148,{{0x1.37bddaf6033e5p-13L,0x1.d9d1f3d186f14p-19L,0x0.0p+0L,0x0.0p+0L,0x1.37bddaf6033e5p-13L,0x1.d9d1f3d186f14p-19L}},{{0x1.37bddaf6033e5p-13L,0x1.0ac299d14fae2p-12L,0x0.0p+0L,0x0.0p+0L,0x1.37bddaf6033e5p-13L,0x1.0ac299d14fae2p-12L}}},
        {321,134,149,{{0x1.3c3bd6f558870p-13L,0x1.36bd41a6e7e5dp-17L,0x0.0p+0L,0x0.0p+0L,0x1.3c3bd6f558870p-13L,0x1.36bd41a6e7e5dp-17L}},{{0x1.3c3bd6f558870p-13L,0x1.61b24f61722c6p-12L,0x0.0p+0L,0x0.0p+0L,0x1.3c3bd6f558870p-13L,0x1.61b24f61722c6p-12L}}},
        {322,132,135,{{0x1.22d5be1010773p-12L,0x1.0eb845830458dp-17L,0x0.0p+0L,0x0.0p+0L,0x1.22d5be1010773p-12L,0x1.0eb845830458dp-17L}},{{0x1.22d5be1010773p-12L,0x1.058b2553c9b04p-11L,0x0.0p+0L,0x0.0p+0L,0x1.22d5be1010773p-12L,0x1.058b2553c9b04p-11L}}},
        {323,133,134,{{0x1.2515aa4292987p-12L,0x1.5452cedd9a26ap-16L,0x0.0p+0L,0x0.0p+0L,0x1.2515aa4292987p-12L,0x1.5452cedd9a26ap-16L}},{{0x1.2515aa4292987p-12L,0x1.5afe68772761cp-11L,0x0.0p+0L,0x0.0p+0L,0x1.2515aa4292987p-12L,0x1.5afe68772761cp-11L}}},
        {324,139,152,{{0x1.1f0491bc4a8c1p-13L,0x1.02f89e4005a3ap-16L,0x0.0p+0L,0x0.0p+0L,0x1.1f0491bc4a8c1p-13L,0x1.02f89e4005a3ap-16L}},{{0x1.1f0491bc4a8c1p-13L,0x1.91087b20da8d8p-12L,0x0.0p+0L,0x0.0p+0L,0x1.1f0491bc4a8c1p-13L,0x1.91087b20da8d8p-12L}}},
        {325,138,153,{{0x1.c837ad8fb3ed1p-14L,0x1.4c68ca099c02cp-16L,0x0.0p+0L,0x0.0p+0L,0x1.c837ad8fb3ed1p-14L,0x1.4c68ca099c02cp-16L}},{{0x1.c837ad8fb3ed1p-14L,0x1.87de5f57653fdp-12L,0x0.0p+0L,0x0.0p+0L,0x1.c837ad8fb3ed1p-14L,0x1.87de5f57653fdp-12L}}},
        {326,136,139,{{0x1.08560f6d24c5cp-12L,0x1.1677db66fe67dp-15L,0x0.0p+0L,0x0.0p+0L,0x1.08560f6d24c5cp-12L,0x1.1677db66fe67dp-15L}},{{0x1.08560f6d24c5cp-12L,0x1.89a1623edfc21p-11L,0x0.0p+0L,0x0.0p+0L,0x1.08560f6d24c5cp-12L,0x1.89a1623edfc21p-11L}}},
        {327,137,138,{{0x1.a16658b9b151cp-13L,0x1.618b1abd81adbp-15L,0x0.0p+0L,0x0.0p+0L,0x1.a16658b9b151cp-13L,0x1.618b1abd81adbp-15L}},{{0x1.a16658b9b151cp-13L,0x1.80cf2e5cfb16cp-11L,0x0.0p+0L,0x0.0p+0L,0x1.a16658b9b151cp-13L,0x1.80cf2e5cfb16cp-11L}}},
        {328,143,156,{{0x1.2480f23693d13p-14L,0x1.431222f36e871p-16L,0x0.0p+0L,0x0.0p+0L,0x1.2480f23693d13p-14L,0x1.431222f36e871p-16L}},{{0x1.2480f23693d13p-14L,0x1.3398f0ab59109p-12L,0x0.0p+0L,0x0.0p+0L,0x1.2480f23693d13p-14L,0x1.3398f0ab59109p-12L}}},
        {329,142,157,{{0x1.91eb1f96fc05ep-16L,0x1.3bc531a0643edp-17L,0x0.0p+0L,0x0.0p+0L,0x1.91eb1f96fc05ep-16L,0x1.3bc531a0643edp-17L}},{{0x1.91eb1f96fc05ep-16L,0x1.fe4579eb1bb4dp-14L,0x0.0p+0L,0x0.0p+0L,0x1.91eb1f96fc05ep-16L,0x1.fe4579eb1bb4dp-14L}}},
        {330,140,143,{{0x1.09ae9ae2b5265p-13L,0x1.550a1628f5c5ep-15L,0x0.0p+0L,0x0.0p+0L,0x1.09ae9ae2b5265p-13L,0x1.550a1628f5c5ep-15L}},{{0x1.09ae9ae2b5265p-13L,0x1.2e2fabdb46ba1p-11L,0x0.0p+0L,0x0.0p+0L,0x1.09ae9ae2b5265p-13L,0x1.2e2fabdb46ba1p-11L}}},
        {331,141,142,{{0x1.6a536d1f7f30cp-15L,0x1.4b9c13b4be5f5p-16L,0x0.0p+0L,0x0.0p+0L,0x1.6a536d1f7f30cp-15L,0x1.4b9c13b4be5f5p-16L}},{{0x1.6a536d1f7f30cp-15L,0x1.f57d367b51426p-13L,0x0.0p+0L,0x0.0p+0L,0x1.6a536d1f7f30cp-15L,0x1.f57d367b51426p-13L}}},
        {332,147,-1,{{0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L}},{{0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L}}},
        {333,146,-1,{{0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L}},{{0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L}}},
        {334,144,147,{{0x1.da34d828c7d0ep-16L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x1.da34d828c7d0ep-16L,0x0.0p+0L}},{{0x1.eda39a14e7c18p-16L,0x1.67333af8cb52cp-17L,0x0.0p+0L,0x0.0p+0L,0x1.eda39a14e7c18p-16L,0x1.67333af8cb52cp-17L}}},
        {335,145,146,{{0x1.4373adbcbb852p-15L,0x1.71cb3252445c5p-26L,0x0.0p+0L,0x0.0p+0L,0x1.4373adbcbb852p-15L,0x1.71cb3252445c5p-26L}},{{0x1.4373adbcbb852p-15L,0x1.62be35e695ed2p-15L,0x0.0p+0L,0x0.0p+0L,0x1.4373adbcbb852p-15L,0x1.62be35e695ed2p-15L}}},
        {336,151,-1,{{0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L}},{{0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L}}},
        {337,150,-1,{{0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L}},{{0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L}}},
        {338,148,151,{{0x1.73b03525692f0p-15L,0x1.c7f3108148bdcp-21L,0x0.0p+0L,0x0.0p+0L,0x1.73b03525692f0p-15L,0x1.c7f3108148bdcp-21L}},{{0x1.73b03525692f0p-15L,0x1.2f027e79388c2p-14L,0x0.0p+0L,0x0.0p+0L,0x1.73b03525692f0p-15L,0x1.2f027e79388c2p-14L}}},
        {339,149,150,{{0x1.7b68006be70d5p-15L,0x1.39f6c56045f62p-19L,0x0.0p+0L,0x0.0p+0L,0x1.7b68006be70d5p-15L,0x1.39f6c56045f62p-19L}},{{0x1.7b68006be70d5p-15L,0x1.918472232298fp-14L,0x0.0p+0L,0x0.0p+0L,0x1.7b68006be70d5p-15L,0x1.918472232298fp-14L}}},
        {340,155,-1,{{0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L}},{{0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L}}},
        {341,154,-1,{{0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L}},{{0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L}}},
        {342,152,155,{{0x1.5a662ffb46971p-15L,0x1.0b0e728d78b41p-18L,0x0.0p+0L,0x0.0p+0L,0x1.5a662ffb46971p-15L,0x1.0b0e728d78b41p-18L}},{{0x1.5a662ffb46971p-15L,0x1.c70801846fd4fp-14L,0x0.0p+0L,0x0.0p+0L,0x1.5a662ffb46971p-15L,0x1.c70801846fd4fp-14L}}},
        {343,153,154,{{0x1.1500fef01288dp-15L,0x1.5aed01a70cfaep-18L,0x0.0p+0L,0x0.0p+0L,0x1.1500fef01288dp-15L,0x1.5aed01a70cfaep-18L}},{{0x1.1500fef01288dp-15L,0x1.bc6faf4896407p-14L,0x0.0p+0L,0x0.0p+0L,0x1.1500fef01288dp-15L,0x1.bc6faf4896407p-14L}}},
        {344,159,-1,{{0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L}},{{0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L}}},
        {345,158,-1,{{0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L}},{{0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L,0x0.0p+0L}}},
        {346,156,159,{{0x1.6598ca7daa953p-16L,0x1.53f15bb39fab9p-18L,0x0.0p+0L,0x0.0p+0L,0x1.6598ca7daa953p-16L,0x1.53f15bb39fab9p-18L}},{{0x1.6598ca7daa953p-16L,0x1.5cb6c6465e57cp-14L,0x0.0p+0L,0x0.0p+0L,0x1.6598ca7daa953p-16L,0x1.5cb6c6465e57cp-14L}}},
        {347,157,158,{{0x1.eec641a0b881fp-18L,0x1.4e2058437d88fp-19L,0x0.0p+0L,0x0.0p+0L,0x1.eec641a0b881fp-18L,0x1.4e2058437d88fp-19L}},{{0x1.eec641a0b881fp-18L,0x1.2121637fc4454p-15L,0x0.0p+0L,0x0.0p+0L,0x1.eec641a0b881fp-18L,0x1.2121637fc4454p-15L}}},
    }};
    for(int exponent:{0,-512})for(bool half:{true,false}){
        stage="captured complete budget replay scale=2^"+std::to_string(exponent)+(half?" half":" full");
        const auto scale=std::ldexp(1.L,exponent);
        std::vector<Cell> cells(originalCells.begin(),originalCells.end());
        for(auto& cell:cells)cell.capacity*=scale;
        std::vector<GeometricVofTransportFace> faces;faces.reserve(originalFaces.size());
        const auto flux=[&](const std::array<long double,6>& a){
            GeometricVofIntegratedFlux f;f.positiveTotal=a[0]*scale;f.negativeTotal=a[1]*scale;
            f.positiveLiquid=a[2]*scale;f.negativeLiquid=a[3]*scale;
            f.positiveGas=a[4]*scale;f.negativeGas=a[5]*scale;
            f.total=f.positiveTotal-f.negativeTotal;f.rawLiquid=f.positiveLiquid-f.negativeLiquid;
            f.liquid=f.rawLiquid;f.gas=f.positiveGas-f.negativeGas;return f;
        };
        for(const auto& original:originalFaces){GeometricVofTransportFace f;
            f.nativeRow=original.row;f.negativeCell=original.negative;f.positiveCell=original.positive;
            f.half=flux(original.half);f.full=flux(original.full);faces.push_back(f);}
        const auto result=LAMP::Interface::Detail::prepareGeometricVofPhaseBudget(cells,faces,half);
        check(result.phaseInventory.size()==160&&result.flux.size()==348,"captured complete native coverage");
        check(result.maximumOutflowCfl<=.9L,"captured gross CFL exceeds frozen .9");
        check(result.maximumTotalDivergence<=512*std::numeric_limits<long double>::epsilon(),
            "captured total divergence exceeds original arithmetic gate");
        check(result.maximumCapacityReconciliation<=512*std::numeric_limits<long double>::epsilon(),
            "captured phase capacity exceeds original arithmetic gate");
        for(const auto& stock:result.phaseInventory)for(const auto amount:stock)
            check(std::isfinite(amount)&&amount>=0,"captured stock invalid");
        for(std::size_t row=0;row<faces.size();++row){const auto& accepted=result.flux[row];
            const auto& raw=half?faces[row].half:faces[row].full;
            const auto gross=[&](long double l,long double g,long double total,long double rawL,long double rawG,int upwind){
                check(std::isfinite(l)&&std::isfinite(g)&&l>=0&&g>=0&&l<=total&&g<=total,
                    "captured strict gross bounds row="+std::to_string(row));++capturedBudgetGrossChecks;
                if(total==0){check(l==0&&g==0,"captured zero-total gross acquired phase");return;}
                check(upwind>=0,"captured nonzero exterior gross");
                if(cells[upwind].alpha==0&&rawL==0&&rawG==total){
                    check(l==0&&g==total,"captured pure gas blend lost exact identity row="+std::to_string(row));
                    ++capturedBudgetPureChecks;
                }
                if(cells[upwind].alpha==1&&rawG==0&&rawL==total){
                    check(g==0&&l==total,"captured pure liquid blend lost exact identity row="+std::to_string(row));
                    ++capturedBudgetPureChecks;
                }
            };
            gross(accepted.positiveLiquid,accepted.positiveGas,raw.positiveTotal,raw.positiveLiquid,raw.positiveGas,faces[row].negativeCell);
            gross(accepted.negativeLiquid,accepted.negativeGas,raw.negativeTotal,raw.negativeLiquid,raw.negativeGas,faces[row].positiveCell);
        }
        const auto& proof=result.flux[190];const auto& original=half?faces[190].half:faces[190].full;
        check(proof.negativeLiquid==0&&proof.negativeGas==original.negativeTotal,
            "original row190 pure negative gas is not exactly original binary total");
        if(!half)check(proof.limiter>0&&proof.limiter<1,"captured full row190 must actually be limited");
        ++capturedBudgetScenarios;
        std::cout<<std::setprecision(17)<<"captured_budget extent="<<(half?"half":"full")<<" scale_exponent="<<exponent
            <<" row190_limiter="<<proof.limiter<<" row190_negative_total="<<original.negativeTotal
            <<" row190_negative_gas="<<proof.negativeGas<<" PASS\n";
    }
}

// Fixed binary inputs expose sub-ulp phase transfer directly, before the
// separate double-alpha and geometry publication constraints are applied.
// Expected values below are elementary phase ledgers, not a test copy of FCT.
void microscopicPhaseBudgets()
{
    using Cell=LAMP::Interface::Detail::GeometricVofPhaseBudgetCell;
    using Budget=LAMP::Interface::Detail::GeometricVofPhaseBudgetResult;
    const long double q=std::ldexp(1.L,-10),micro=std::ldexp(1.L,-100);
    long double fixtureScale=1;
    const auto face=[&](int row,int negative,int positive,long double hpL,long double hpG,
        long double hnL,long double hnG){
        GeometricVofTransportFace f;f.nativeRow=row;f.negativeCell=negative;f.positiveCell=positive;
        GeometricVofIntegratedFlux gross;gross.positiveTotal=gross.negativeTotal=q;gross.total=0;
        gross.positiveLiquid=hpL;gross.positiveGas=hpG;gross.negativeLiquid=hnL;gross.negativeGas=hnG;
        gross.rawLiquid=hpL-hnL;gross.liquid=gross.rawLiquid;gross.gas=hpG-hnG;
        f.half=f.full=gross;return f;
    };
    const auto phase=[&](long double actual,long double expected,const std::string& name){
        expected*=fixtureScale;
        if(expected==0){check(actual==0,name+" must remain theoretical exact zero");return;}
        near(actual,expected,std::abs(expected),name,conservationTolerance);
        peakMicroRelative=std::max(peakMicroRelative,std::abs(actual-expected)/std::abs(expected));
    };
    const auto run=[&](const std::string& name,const std::vector<Cell>& cells,
        const std::vector<GeometricVofTransportFace>& faces,const std::function<void(const Budget&)>& expected){
        // A common physical rescaling must leave the dimensionless limiter
        // and relative micro-phase accuracy unchanged. Scale expected ledgers
        // directly; do not normalize the kernel outputs before comparison.
        for(int scaleExponent:{0,-512}){
            fixtureScale=std::ldexp(1.L,scaleExponent);
            auto scaledCells=cells;auto scaledFaces=faces;
            for(auto& cell:scaledCells)cell.capacity*=fixtureScale;
            for(auto& physicalFace:scaledFaces)for(auto* flux:{&physicalFace.half,&physicalFace.full}){
                flux->total*=fixtureScale;flux->positiveTotal*=fixtureScale;flux->negativeTotal*=fixtureScale;
                flux->rawLiquid*=fixtureScale;flux->liquid*=fixtureScale;flux->gas*=fixtureScale;flux->mass*=fixtureScale;
                flux->positiveLiquid*=fixtureScale;flux->negativeLiquid*=fixtureScale;
                flux->positiveGas*=fixtureScale;flux->negativeGas*=fixtureScale;
            }
            for(bool half:{true,false}){
                stage="fixed microscopic phase budget "+name+" scale=2^"+std::to_string(scaleExponent)+(half?" half":" full");
                const auto result=LAMP::Interface::Detail::prepareGeometricVofPhaseBudget(scaledCells,scaledFaces,half);
                check(result.phaseInventory.size()==cells.size()&&result.flux.size()==faces.size(),"micro kernel output coverage");
                check(result.maximumOutflowCfl<=.9L,"micro fixture exceeded frozen gross CFL");
                check(result.maximumTotalDivergence==0,"balanced micro fixture acquired native divergence");
                check(result.maximumCapacityReconciliation<=512*std::numeric_limits<long double>::epsilon(),
                    "micro phase capacity exceeded unchanged arithmetic gate");
                for(const auto& inventory:result.phaseInventory)
                    check(std::isfinite(inventory[0])&&std::isfinite(inventory[1])&&inventory[0]>=0&&inventory[1]>=0,
                        "micro kernel published invalid phase inventory");
                expected(result);++microScenarios;
                std::cout<<"micro_case="<<name<<" extent="<<(half?"half":"full")<<" scale_exponent="<<scaleExponent<<" PASS\n";
            }
        }
    };
    run("birth-without-gas-outflow",{{1,.5},{1,1}},
        {face(0,0,1,q,micro,q,0)},[&](const Budget& r){
            check(r.flux[0].limiter==1&&r.limitedFaceExtents==0,"unbounded micro birth was unexpectedly limited");
            phase(r.flux[0].positiveGas,micro,"birth positive gross gas");
            phase(r.flux[0].negativeGas,0,"birth absent negative gross gas");
            phase(r.flux[0].liquid,-micro,"birth exact signed liquid");
            phase(r.flux[0].gas,micro,"birth exact signed gas");
            phase(r.phaseInventory[1][1],micro,"birth recipient gas inventory");
            phase(r.phaseInventory[0][0],.5L,"birth donor liquid inventory");
            phase(r.phaseInventory[0][1],.5L,"birth donor gas inventory");
            phase(r.phaseInventory[1][0],1,"birth recipient liquid inventory");
        });
    run("counterflow-net-cancellation",{{1,.5},{1,.5}},
        {face(0,0,1,q,micro,q,2*micro)},[&](const Budget& r){
            check(r.flux[0].limiter==1&&r.limitedFaceExtents==0,"micro counterflow was unexpectedly limited");
            check(r.flux[0].positiveLiquid==q*fixtureScale&&r.flux[0].negativeLiquid==q*fixtureScale,
                "counterflow fixture must retain coincident rounded dominant gross values");
            phase(r.flux[0].positiveGas,micro,"counterflow positive gross gas");
            phase(r.flux[0].negativeGas,2*micro,"counterflow negative gross gas");
            phase(r.flux[0].liquid,micro,"counterflow exact signed liquid");
            phase(r.flux[0].gas,-micro,"counterflow exact signed gas");
            for(const auto& inventory:r.phaseInventory){phase(inventory[0],.5L,"counterflow liquid stock");phase(inventory[1],.5L,"counterflow gas stock");}
        });
    // Face 1 transfers micro liquid from cell 2 with no geometric correction.
    // Face 0 alone has equal raw demands. Its donor now has micro additional
    // liquid budget: the two controlling ratios differ below a scalar ulp.
    // The receiver saturates first, leaving micro liquid in the donor.
    const double b=std::ldexp(1.,-90);
    run("near-equal-limiter-slack",{{1,0},{1,1},{1,b}},
        {face(0,0,1,q,0,0,q),face(1,0,2,0,q,micro,q)},[&](const Budget& r){
            check(r.flux[0].limiter==.5L&&r.flux[1].limiter==1&&r.limitedFaceExtents==1,
                "near-equal limiter fixture did not use one shared limiting face");
            phase(r.flux[0].liquid,0,"slack shared face net liquid");
            phase(r.flux[0].gas,0,"slack shared face net gas");
            phase(r.flux[1].liquid,-micro,"slack second face signed liquid");
            phase(r.flux[1].gas,micro,"slack second face signed gas");
            phase(r.phaseInventory[0][0],micro,"slack donor retained liquid inventory");
            phase(r.phaseInventory[0][1],1,"slack donor gas inventory");
            phase(r.phaseInventory[1][0],1,"slack receiver liquid inventory");
            phase(r.phaseInventory[1][1],0,"slack receiver saturated gas inventory");
            phase(r.phaseInventory[2][0],static_cast<long double>(b)*(1-q),"slack second donor retained liquid");
            phase(r.phaseInventory[2][1],1,"slack second donor gas inventory");
        });
    const double a=std::ldexp(1.,-12);
    const long double high=static_cast<long double>(a)*(2-q);
    run("exact-saturated-zero",{{1,a},{1,0}},
        {face(0,0,1,high,q-high,0,q)},[&](const Budget& r){
            check(r.flux[0].limiter==.5L&&r.limitedFaceExtents==1,"saturated limiter must be exactly one half");
            phase(r.flux[0].positiveLiquid,a,"saturation accepted positive liquid");
            phase(r.flux[0].negativeLiquid,0,"saturation absent negative liquid");
            phase(r.flux[0].liquid,a,"saturation signed liquid");
            phase(r.flux[0].gas,-static_cast<long double>(a),"saturation signed gas");
            phase(r.phaseInventory[0][0],0,"saturated donor exact zero liquid");
            phase(r.phaseInventory[0][1],1,"saturated donor gas inventory");
            phase(r.phaseInventory[1][0],a,"saturation recipient liquid inventory");
            phase(r.phaseInventory[1][1],1-static_cast<long double>(a),"saturation recipient gas inventory");
        });
}
// Analytic tests of the ORIGINAL public points and complete trajectory. These
// expected integrals come from trace lengths/change of variable, independently
// of the event roots, lambda bounds and integration arithmetic in Core.
void coreFaceTimeOriginalGeometry()
{
    namespace FaceTime=LAMP::Geometry::PlanarFaceTime;
    using Poly=std::vector<Point>;
    struct Fixture
    {
        const char* name;
        Poly polygon;
        Point a,b,v0,v1;
        long double duration;
        std::array<long double,4> expected;
        bool translationSafe;
    };
    const auto rectangle=[](long double x0,long double x1,long double y0,long double y1)
    {return Poly{{x0,y0},{x1,y0},{x1,y1},{x0,y1}};};
    const long double e=std::ldexp(1.L,-100),contactE=std::ldexp(1.L,-120),q=std::ldexp(1.L,-28);
    const Point cornerV0{-.0485L,.00428L},cornerV1{-2*cornerV0.x,-2*cornerV0.y};
    const std::vector<Fixture> fixtures{
        // f=t*(1-1.5*t): f>0 gives departure x>0; f<0 gives y>0.
        // At f=0 only one corner touches, hence all four exact zeros.
        {"zero-corner-original-json",{{0,0},{-1,0},{0,-1}},
            {0,0},{0,1},cornerV0,cornerV1,1.1788732394366197L,{0,0,0,0},true},
        // Wet length is 1-t/e on [0,e], so the positive integral is e/2.
        {"positive-micro-triangle",{{0,0},{0,1},{-e,0}},
            {0,0},{0,1},{1,0},{1,0},1,{e/2,0,e/2,0},true},
        // departure x=-(t-.5)^2. The right polygon touches at one instant.
        {"zero-isolated-half-contact",rectangle(0,1,0,1),
            {-.25L,0},{-.25L,1},{-1,0},{1,0},1,{0,0,0,0},true},
        // The left thin polygon has |t-.5|<=2^-60. Both true roots round
        // to the same binary64 time, yet each side transports contactE.
        {"positive-micro-half-window",rectangle(-contactE,0,0,1),
            {-.25L,0},{-.25L,1},{-1,0},{1,0},1,{0,contactE,contactE,contactE},true},
        // x=f, y=lambda-f. Wet physical length is e for 0<=f<=e and
        // e+f for -e<=f<=0. |normal|dt=|df|. The returning branch near
        // 2/3 has no binary64 interior time but transports 3*e^2/2.
        {"positive-micro-two-thirds-window",rectangle(-e,e,0,e),
            {0,0},{0,1},{-1,1},{2,-2},1,{0,e*e,1.5L*e*e,e*e},false},
        // X_reverse(t)=X_forward(1-t). Forming p-a in ordinary floating
        // arithmetic before accurate predicates would already erase e.
        {"positive-micro-time-reflection",rectangle(-e,e,0,e),
            {-.5L,.5L},{-.5L,1.5L},{-2,2},{1,-1},1,{0,1.5L*e*e,e*e,1.5L*e*e},false},
        // Wet length is f for 0<=f<=e, zero otherwise; each gross is e^2/2.
        {"positive-micro-shifted-corner",rectangle(-1,e,-1,0),
            {0,0},{0,1},{-1,1},{2,-2},1,{0,e*e/2,e*e/2,e*e/2},true},
        {"zero-shifted-corner-limit",rectangle(-1,0,-1,0),
            {0,0},{0,1},{-1,1},{2,-2},1,{0,0,0,0},true},
        // Area is det((1+q,1),(q,q))=q^2. Departure x=4t crosses the
        // entire actual parallelogram before half time, so flux equals area.
        {"positive-thin-parallelogram",{{0,0},{1+q,1},{1+2*q,1+q},{q,q}},
            {0,0},{0,2},{-4,0},{-4,0},1,{0,q*q,0,q*q},true}
    };
    unsigned failures=0;
    const auto verifyGross=[&](const FaceTime::Gross& actual,const std::array<long double,4>& expected)
    {
        const std::array<long double,4> values{actual.halfPositive,actual.halfNegative,actual.fullPositive,actual.fullNegative};
        const std::array<const char*,4> names{"half-positive","half-negative","full-positive","full-negative"};
        for(std::size_t i=0;i<values.size();++i){
            const auto value=values[i],target=expected[i];
            check(std::isfinite(value)&&value>=0,std::string(names[i])+" invalid gross integral");
            if(target==0){
                check(value==0,std::string(names[i])+" theoretical zero became geometric phase flux");
                ++faceTimeZeroChecks;
            }else{
                check(std::isfinite(target)&&target>0,"analytic positive target unrepresentable");
                check(value>0,std::string(names[i])+" deleted a strictly positive actual phase trace");
                const auto residual=std::abs(value-target)/target;
                peakFaceTimeRelative=std::max(peakFaceTimeRelative,residual);
                std::ostringstream detail;detail<<std::setprecision(17)<<names[i]<<" actual="<<value<<" expected="<<target<<" micro_relative="<<residual;
                check(residual<=conservationTolerance,detail.str());++faceTimePositiveChecks;
            }
        }
        check(std::isfinite(actual.maximumRelativeUncertainty)&&actual.maximumRelativeUncertainty>=0,
            "invalid face-time certified uncertainty diagnostic");
        check(actual.maximumRelativeUncertainty<=32*std::numeric_limits<long double>::epsilon(),
            "face-time relative uncertainty exceeded frozen positive-integral certificate");
        check(actual.maximumRootRefinements<=256,"face-time root refinement exceeded frozen work budget");
        maximumFaceTimeRootRefinements=std::max(maximumFaceTimeRootRefinements,actual.maximumRootRefinements);
        peakFaceTimeUncertainty=std::max(peakFaceTimeUncertainty,actual.maximumRelativeUncertainty);
    };
    for(const auto& fixture:fixtures)for(int scaleExponent:{0,-512})for(bool rz:{false,true})
        for(bool clockwise:{false,true})for(bool reverseFace:{false,true}){
            // A common axial frame translation of 16 is applied only where
            // every original public Point difference remains representable.
            // Do not translate y~2^-100 rectangles into y~16 or combine that
            // large frame translation with the 2^-512 physical-scale case.
            const unsigned translations=scaleExponent==0&&fixture.translationSafe?2:1;
            for(unsigned translated=0;translated<translations;++translated){
                const auto scale=std::ldexp(1.L,scaleExponent),shift=translated?16.L:0.L;
                auto polygon=fixture.polygon;auto a=fixture.a,b=fixture.b,v0=fixture.v0,v1=fixture.v1;
                for(auto& p:polygon){p.y=p.y*scale-shift;}
                a.y=a.y*scale-shift;b.y=b.y*scale-shift;v0.y*=scale;v1.y*=scale;
                if(clockwise)std::reverse(polygon.begin(),polygon.end());
                if(reverseFace)std::swap(a,b);
                const Point origin{1-a.x,shift}; // Actual Eulerian radial face coordinate is exactly 1.
                auto expected=fixture.expected;
                for(auto& value:expected)value*=scale*(rz?twoPi:1);
                stage=std::string("Core original face-time ")+fixture.name+(rz?" RZ":" XY")+
                    (clockwise?" CW":" CCW")+(reverseFace?" reversed-endpoints":" forward-endpoints")+
                    " axial-scale=2^"+std::to_string(scaleExponent)+(translated?" axial-frame-shift=16":" axial-frame-shift=0");
                ++faceTimeAttempts;
                try{
                    const auto actual=FaceTime::integratePolygon(polygon,a,b,v0,v1,0,fixture.duration,origin,rz,512,{});
                    verifyGross(actual,expected);++faceTimeScenarios;
                }catch(const std::exception& error){
                    ++failures;std::cerr<<"face_time_case FAIL stage="<<stage<<" "<<error.what()<<"\n";
                }
            }
        }
    // A source owner may supply one known closed-wall affine reflection.
    // Keep ORIGINAL x in [0,contactE] plus x'=.5-x: materializing x' first
    // would collapse the two distinct actual source endpoints to .5.
    // y'=y or scale-y both describe the same physical [0,scale] trace.
    for(int scaleExponent:{0,-512})for(int signY:{1,-1})for(bool rz:{false,true})
        for(bool clockwise:{false,true})for(bool reverseFace:{false,true}){
            const auto scale=std::ldexp(1.L,scaleExponent);
            auto polygon=rectangle(0,contactE,0,scale);
            Point a{.25L,0},b{.25L,scale};
            if(clockwise)std::reverse(polygon.begin(),polygon.end());
            if(reverseFace)std::swap(a,b);
            const FaceTime::Reflection reflection{-1,signY,{.5L,signY<0?scale:0}};
            const auto amount=contactE*scale*(rz?twoPi:1);
            const std::array<long double,4> expected{0,amount,amount,amount};
            stage=std::string("Core original reflected micro-half-window")+(rz?" RZ":" XY")+
                (clockwise?" source-CW":" source-CCW")+(reverseFace?" reversed-endpoints":" forward-endpoints")+
                " reflectionY="+std::to_string(signY)+" axial-scale=2^"+std::to_string(scaleExponent);
            ++faceTimeAttempts;
            try{
                const auto actual=FaceTime::integratePolygon(polygon,a,b,{-1,0},{1,0},0,1,{.75L,0},rz,512,reflection);
                verifyGross(actual,expected);++faceTimeScenarios;
            }catch(const std::exception& error){
                ++failures;std::cerr<<"face_time_case FAIL stage="<<stage<<" "<<error.what()<<"\n";
            }
        }
    stage="Core original face-time explicit input guards";
    const auto& fixture=fixtures[1];
    rejects([&]{FaceTime::integratePolygon(fixture.polygon,fixture.a,{.25L,1},fixture.v0,fixture.v1,0,1,{1,0},false);},
        "Core face-time non-axis-aligned actual face");
    rejects([&]{FaceTime::integratePolygon(fixture.polygon,fixture.a,fixture.b,fixture.v0,fixture.v1,1,1,{1,0},false);},
        "Core face-time nonpositive complete interval");
    rejects([&]{FaceTime::integratePolygon(fixture.polygon,fixture.a,fixture.b,fixture.v0,fixture.v1,0,1,{1,0},false,2);},
        "Core face-time invalid finite event budget");
    std::cout<<std::setprecision(17)<<"face_time_original_geometry cases="<<faceTimeAttempts<<" passed="<<faceTimeScenarios
        <<" failed="<<failures<<" positive_checks="<<faceTimePositiveChecks<<" exact_zero_checks="<<faceTimeZeroChecks
        <<" peak_micro_relative="<<peakFaceTimeRelative<<" peak_relative_uncertainty="<<peakFaceTimeUncertainty
        <<" maximum_root_refinements="<<maximumFaceTimeRootRefinements<<"\n";
    check(failures==0,"frozen original-geometry face-time case failures="+std::to_string(failures));
}
void genericMatrix()
{
    for(bool rz:{false,true})for(bool reverse:{false,true}){
        stage=std::string(rz?"RZ":"XY")+(reverse?" temporal reversal":" general nonconstant native flow");
        auto mesh=meshFor();const auto old=publish(*mesh,rz,lineInput(*mesh,rz));const auto saved=old.alpha();const auto identity=old.identity();
        auto a=normals(*mesh,old),b=reverse?normals(*mesh,old,0,-2):normals(*mesh,old,1);
        const auto next=nextClock(old,deltaTime(*mesh,old,a,b,.12L));
        const auto slab=hydroSlab(*mesh,old,next,a,b);const auto result=Transport::prepare(*mesh,old,slab,100,1);
        verify(*mesh,old,slab,result);checkRaw(*mesh,old,slab,result,!reverse);
        check(old.identity()==identity&&old.alpha()==saved,"successful prepare mutated accepted source");
        const auto cslab=hydroSlab(*mesh,old,next,normals(*mesh,old,2,.8L),normals(*mesh,old,0,.9L));
        const auto corrected=Transport::prepare(*mesh,old,cslab,100,1);verify(*mesh,old,cslab,corrected);
        check(corrected.sourceIdentity()==result.sourceIdentity(),"C advanced A candidate instead of same accepted old fan");
        check(sameClock(corrected.endTime(),result.endTime()),"C changed global clock");
        check(old.alpha()==saved&&old.identity()==identity,"C staging mutated accepted source");
        check(result.endpoint().alpha()!=corrected.endpoint().alpha(),"A/C fixture velocities did not alter transported inventory");
    }
    check(counterflowFaces>0,"normal temporal reversal did not retain positive/negative gross flux");
}
void constantDiagonalPatch()
{
    stage="constant diagonal XY interior patch with closed outer walls";auto mesh=meshFor(false);
    const auto old=publish(*mesh,false,lineInput(*mesh,false,1.03L));const auto a=normals(*mesh,old,3);
    const auto slab=hydroSlab(*mesh,old,nextClock(old,deltaTime(*mesh,old,a,a,.2L)),a,a);
    const auto result=Transport::prepare(*mesh,old,slab,100,1);verify(*mesh,old,slab,result);
    unsigned count=0;long double difference=0;
    for(std::size_t row=0;row<result.faces().size();++row){
        const auto f=old.domain().face(*mesh,static_cast<int>(row));if(f.negativeCell<0||f.positiveCell<0)continue;
        const auto velocity=slab.beginVelocity()[row];
        if(std::abs(velocity.x-1)>1e-12L||std::abs(velocity.y-1)>1e-12L)continue;
        if(old.alpha()[f.negativeCell]==old.alpha()[f.positiveCell])continue;
        const auto ends=faceEndpoints(*mesh,f);const auto dt=slab.endTime().end-slab.endTime().begin;
        if(ends[0].x-dt<=0||ends[1].x+dt>=1||ends[0].y-dt<=0||ends[1].y+dt>=1)continue;
        const auto reference=rawReference(*mesh,old,slab,static_cast<int>(row),1);
        const auto scale=dt*f.area;
        near(result.faces()[row].full.rawLiquid,reference.value,scale,"constant diagonal old-fan exact reference",rawTolerance);
        difference=std::max(difference,std::abs(reference.value-reference.normalOnly)/scale);++count;++rawChecks;
    }
    check(count>=2,"native local fit did not reproduce the constant diagonal interior patch");
    check(difference>1e-5L,"constant diagonal fixture failed to detect the missing transverse quadratic term");
}
void pureMatrix()
{
    for(bool rz:{false,true})for(double alpha:{0.,1.}){
        stage=std::string(rz?"RZ":"XY")+" pure phase transport";auto mesh=meshFor();
        const auto old=publish(*mesh,rz,std::vector<double>(mesh->cellIndexNum(),alpha));const auto a=normals(*mesh,old),b=normals(*mesh,old,1);
        const auto slab=hydroSlab(*mesh,old,nextClock(old,deltaTime(*mesh,old,a,b,.25L)),a,b);
        const auto result=Transport::prepare(*mesh,old,slab,100,1);verify(*mesh,old,slab,result);
        check(result.midpoint().alpha()==old.alpha()&&result.endpoint().alpha()==old.alpha(),"pure phase was not preserved exactly");
    }
}
void limiterMatrix()
{
    const auto previous=limitedFaces;
    for(bool rz:{false,true})for(long double intercept:{.61L,.913L,1.19L}){
        stage=std::string(rz?"RZ":"XY")+" coupled geometric flux limiter";auto mesh=meshFor();
        const auto old=publish(*mesh,rz,lineInput(*mesh,rz,intercept));const auto a=normals(*mesh,old,2),b=normals(*mesh,old,0,-2);
        const auto slab=directSlab(*mesh,old,nextClock(old,deltaTime(*mesh,old,a,b,.85L)),a,b);
        const auto result=Transport::prepare(*mesh,old,slab,100,1);verify(*mesh,old,slab,result);
    }
    check(limitedFaces>previous,"frozen generic large-CFL fixtures did not activate the coupled limiter");
}
void sequentialCandidates()
{
    for(bool rz:{false,true}){
        stage=std::string(rz?"RZ":"XY")+" actual multi-step old-fan rebuild";auto mesh=meshFor();
        auto state=publish(*mesh,rz,lineInput(*mesh,rz));const auto first=state.report().inventory;
        for(unsigned step=0;step<3;++step){const auto a=normals(*mesh,state,step%3),b=normals(*mesh,state,(step+1)%3);
            const auto slab=hydroSlab(*mesh,state,nextClock(state,deltaTime(*mesh,state,a,b,.08L)),a,b);
            const auto result=Transport::prepare(*mesh,state,slab,100,1);verify(*mesh,state,slab,result);
            state=result.endpoint();near(state.report().inventory,first,first,"multi-step closed stock",conservationTolerance);}
    }
}
void failuresAndFields()
{
    stage="full-clock binding and atomic failure";auto mesh=meshFor();const auto old=publish(*mesh,false,lineInput(*mesh,false));
    const auto identity=old.identity();const auto saved=old.alpha();const auto a=normals(*mesh,old),b=normals(*mesh,old,1);
    const auto next=nextClock(old,deltaTime(*mesh,old,a,b,.12L));const auto slab=hydroSlab(*mesh,old,next,a,b);
    CellField<double> alphaField(*mesh,.375),pressure(*mesh,.625);VEdgeField<double> u(*mesh,.125);HEdgeField<double> v(*mesh,-.25);
    const auto* alphaAddress=alphaField.value.get();const auto* uAddress=u.value.get();const auto* vAddress=v.value.get();
    const auto fieldUnchanged=[&]{check(alphaField.value.get()==alphaAddress&&u.value.get()==uAddress&&v.value.get()==vAddress,"candidate preparation replaced Field storage");
        for(int i=0;i<alphaField.size();++i)check(alphaField.value[i]==.375&&pressure.value[i]==.625,"candidate changed cell Field");
        for(int i=0;i<u.size();++i)check(u.value[i]==.125,"candidate changed native radial Field");
        for(int i=0;i<v.size();++i)check(v.value[i]==-.25,"candidate changed native axial Field");
        check(old.identity()==identity&&old.alpha()==saved,"failed candidate polluted old immutable VOF");};
    auto bad=next;bad.step+=1;rejects([&]{directSlab(*mesh,old,bad,a,b);},"skipped accepted step");
    bad=next;bad.begin=std::nextafter(bad.begin,1.);rejects([&]{directSlab(*mesh,old,bad,a,b);},"clock begin differs by one ULP");
    bad=next;bad.position=Position::Midpoint;rejects([&]{directSlab(*mesh,old,bad,a,b);},"slab end position is Midpoint");
    bad=next;bad.end=std::nextafter(bad.begin,std::numeric_limits<double>::infinity());
    rejects([&]{directSlab(*mesh,old,bad,a,b);},"full interval has no representable internal Midpoint");
    const auto layout1=HydroDofLayout::build(*mesh,old.domain()),layout2=HydroDofLayout::build(*mesh,old.domain());
    const HydroNativeValues begin1(layout1,a,std::vector<double>(layout1.pressureCount())),
        end2(layout2,b,std::vector<double>(layout2.pressureCount()));
    rejects([&]{prepareHydroVofAdvectorSlab(*mesh,old,next,begin1,end2);},"Hydro endpoint layout identity differs");
    const auto detachedLayout=HydroDofLayout::build(*mesh,domain(*mesh,false));
    const HydroNativeValues detachedBegin(detachedLayout,a,std::vector<double>(detachedLayout.pressureCount())),
        detachedEnd(detachedLayout,b,std::vector<double>(detachedLayout.pressureCount()));
    rejects([&]{prepareHydroVofAdvectorSlab(*mesh,old,next,detachedBegin,detachedEnd);},"rebuilt Hydro domain has different mesh binding");
    auto wrong=a;wrong.pop_back();rejects([&]{directSlab(*mesh,old,next,wrong,b);},"incomplete native coverage");
    wrong=a;wrong[0]=std::numeric_limits<double>::quiet_NaN();rejects([&]{directSlab(*mesh,old,next,wrong,b);},"nonfinite native velocity");
    wrong=a;for(int row=0;row<old.domain().nativeFaceCount();++row)if(old.domain().face(*mesh,row).outwardSign!=0){wrong[row]=.01;break;}
    rejects([&]{directSlab(*mesh,old,next,wrong,b);},"moving closed wall");
    const auto foreign=publish(*mesh,false,saved);rejects([&]{Transport::prepare(*mesh,foreign,slab,100,1);},"foreign same-clock source owner");
    rejects([&]{Transport::prepare(*mesh,old,slab,0,1);},"nonpositive liquid density");
    GeometricVofTransportConfig config;config.maximumOutflowCfl=.00001L;
    rejects([&]{Transport::prepare(*mesh,old,slab,100,1,config);},"CFL violation");
    const auto result=Transport::prepare(*mesh,old,slab,100,1);verify(*mesh,old,slab,result);fieldUnchanged();
    auto other=meshFor();const auto otherState=publish(*other,false,lineInput(*other,false));
    rejects([&]{slab.validateFor(*other,otherState);},"foreign mesh identity");fieldUnchanged();
    // The ordinary unit-domain fixture has |integrated flux|<1, so the
    // largest finite density need not overflow. A physically larger complete
    // domain tests the material product after all geometric/CFL guards pass.
    DictionaryMesh largeDictionary;largeDictionary.Nx=largeDictionary.Ny=8;
    largeDictionary.maxCoord={32,32};largeDictionary.maxCellNum=256;
    largeDictionary.maxInitialIterations=0;largeDictionary.maxLevel_local={0};
    QuadtreeMesh largeMesh(largeDictionary);
    const auto pure=publish(largeMesh,false,std::vector<double>(largeMesh.cellIndexNum(),1.));
    const auto pureIdentity=pure.identity();const auto pureSaved=pure.alpha();
    std::vector<double> largeNormals(pure.domain().nativeFaceCount());
    for(int row=0;row<pure.domain().nativeFaceCount();++row){const auto f=pure.domain().face(largeMesh,row);
        const auto ends=faceEndpoints(largeMesh,f);
        const Point first{ends[0].x/32,ends[0].y/32},last{ends[1].x/32,ends[1].y/32};
        const auto q=(f.horizontal?-1:1)*4096*(psi(last,0)-psi(first,0));
        largeNormals[row]=f.outwardSign!=0?0:static_cast<double>(q/f.area);}
    const auto largeSlab=hydroSlab(largeMesh,pure,nextClock(pure,deltaTime(largeMesh,pure,largeNormals,largeNormals,.3L)),largeNormals,largeNormals);
    const auto finiteMaterial=Transport::prepare(largeMesh,pure,largeSlab,1,1);
    long double largest=0;for(const auto& face:finiteMaterial.faces())largest=std::max(largest,std::abs(face.full.liquid));
    check(largest>1,"mass-overflow fixture does not have a material multiplier greater than one");
    bool massRejected=false;
    try{Transport::prepare(largeMesh,pure,largeSlab,std::numeric_limits<long double>::max(),1);}
    catch(const std::exception& error){
        massRejected=std::string(error.what()).find("unrepresentable shared material integral")!=std::string::npos;
        check(massRejected,std::string("overflow rejected by unrelated guard: ")+error.what());
    }
    check(massRejected,"finite density/material product overflow was not explicitly rejected");++guards;++scenarios;
    check(pure.identity()==pureIdentity&&pure.alpha()==pureSaved,"mass failure polluted its accepted source");fieldUnchanged();
}
}
// BEGIN V4a1a ADDED ANCHORED TEST BLOCK
namespace
{
// Independent twofold reference: supplied components remain separate until final scalar output.
struct AnchoredRefNumber {long double hi=0,lo=0;};
AnchoredRefNumber arSum(long double a,long double b)
{const auto s=a+b,v=s-a;return {s,(a-(s-v))+(b-v)};}
AnchoredRefNumber arAdd(AnchoredRefNumber a,AnchoredRefNumber b)
{const auto s=arSum(a.hi,b.hi),t=arSum(a.lo,b.lo),u=arSum(s.lo,t.hi),v=arSum(s.hi,u.hi);return arSum(v.hi,v.lo+u.lo+t.lo);}
AnchoredRefNumber arNeg(AnchoredRefNumber a){return {-a.hi,-a.lo};}
AnchoredRefNumber arSub(AnchoredRefNumber a,AnchoredRefNumber b){return arAdd(a,arNeg(b));}
AnchoredRefNumber arMul(AnchoredRefNumber a,AnchoredRefNumber b)
{const auto p=a.hi*b.hi;auto e=arAdd({std::fma(a.hi,b.hi,-p),0},{a.hi*b.lo,0});e=arAdd(e,{a.lo*b.hi,0});e=arAdd(e,{a.lo*b.lo,0});return arAdd({p,0},e);}
AnchoredRefNumber arDiv(AnchoredRefNumber a,AnchoredRefNumber b)
{const auto q=a.hi/b.hi;const auto r=arSub(a,arMul(b,{q,0}));return arAdd({q,0},{(r.hi+r.lo)/b.hi,0});}
long double arValue(AnchoredRefNumber a){return a.hi+a.lo;}
struct AnchoredRefPoint {AnchoredRefNumber x,y;};
using AnchoredSite=LAMP::Geometry::PlanarCellCut::AnchoredPoint;
using AnchoredSites=std::vector<AnchoredSite>;
using AnchoredRefPolygon=std::vector<AnchoredRefPoint>;
AnchoredRefPoint arPoint(AnchoredSite p){return {arSum(p.anchor.x,p.offset.x),arSum(p.anchor.y,p.offset.y)};}
AnchoredRefPoint arSubPoint(AnchoredRefPoint a,AnchoredRefPoint b){return {arSub(a.x,b.x),arSub(a.y,b.y)};}
AnchoredRefNumber arCross(AnchoredRefPoint a,AnchoredRefPoint b){return arSub(arMul(a.x,b.y),arMul(a.y,b.x));}
AnchoredRefPolygon arPoints(const AnchoredSites& p){AnchoredRefPolygon out;for(auto q:p)out.push_back(arPoint(q));return out;}
struct AnchoredRefMoments {AnchoredRefNumber area,x,y;};
AnchoredRefMoments arMoments(const AnchoredRefPolygon& p)
{
    if(p.size()<3)return {};const auto base=p.front();AnchoredRefNumber area,x,y;
    for(std::size_t i=0;i<p.size();++i){const auto a=arSubPoint(p[i],base),b=arSubPoint(p[(i+1)%p.size()],base);const auto c=arCross(a,b);
        area=arAdd(area,c);x=arAdd(x,arMul(arAdd(a.x,b.x),c));y=arAdd(y,arMul(arAdd(a.y,b.y),c));}
    if(arValue(area)<0){area=arNeg(area);x=arNeg(x);y=arNeg(y);}
    area=arDiv(area,{2,0});x=arAdd(arMul(base.x,area),arDiv(x,{6,0}));y=arAdd(arMul(base.y,area),arDiv(y,{6,0}));return {area,x,y};
}
AnchoredRefNumber arMeasure(const AnchoredRefPolygon& p,LAMP::Geometry::PlanarCellCut::Point origin,bool rz)
{const auto m=arMoments(p);return rz?arMul({static_cast<long double>(CartesianMeasure::twoPi),0},arAdd(arMul({origin.x,0},m.area),m.x)):m.area;}
AnchoredRefPoint arCrossing(AnchoredRefPoint a,AnchoredRefPoint b,long double va,long double vb)
{
    const auto wa=std::abs(va),wb=std::abs(vb);const auto t=arDiv({wa,0},arAdd({wa,0},{wb,0}));
    const auto d=arSubPoint(b,a);return {arAdd(a.x,arMul(t,d.x)),arAdd(a.y,arMul(t,d.y))};
}
AnchoredRefPolygon arTriangleNegative(const std::array<AnchoredSite,3>& triangle,const std::array<long double,3>& scalar)
{
    AnchoredRefPolygon p;for(unsigned i=0;i<3;++i){const auto j=(i+1)%3;const auto a=arPoint(triangle[i]),b=arPoint(triangle[j]);
        if(scalar[i]<=0)p.push_back(a);if((scalar[i]<0&&scalar[j]>0)||(scalar[i]>0&&scalar[j]<0))p.push_back(arCrossing(a,b,scalar[i],scalar[j]));}return p;
}
AnchoredRefPolygon arHalfPlane(const AnchoredRefPolygon& p,bool top)
{
    AnchoredRefPolygon out;if(p.empty())return out;
    for(std::size_t i=0;i<p.size();++i){const auto a=p[i],b=p[(i+1)%p.size()];auto va=arSub(a.y,{.5L,0}),vb=arSub(b.y,{.5L,0});if(top){va=arNeg(va);vb=arNeg(vb);}
        const bool ina=arValue(va)<=0,inb=arValue(vb)<=0;if(ina)out.push_back(a);
        if(ina!=inb){const auto t=arDiv(va,arSub(va,vb));const auto d=arSubPoint(b,a);out.push_back({arAdd(a.x,arMul(t,d.x)),arAdd(a.y,arMul(t,d.y))});}}return out;
}
bool arSame(AnchoredSite a,AnchoredSite b)
{const auto d=arSubPoint(arPoint(a),arPoint(b));return d.x.hi==0&&d.x.lo==0&&d.y.hi==0&&d.y.lo==0;}
unsigned anchoredFaceTimeCases=0,anchoredFaceTimeImages=0,anchoredFaceTimePositive=0,anchoredFaceTimeZeros=0,anchoredFaceTimeGuards=0;
long double anchoredFaceTimePeak=0;
void anchoredTimeCheck(bool ok,const char* reason){if(!ok)throw std::runtime_error(std::string("anchored FaceTime: ")+reason);}
void anchoredTimeNear(long double actual,long double expected)
{anchoredTimeCheck(std::isfinite(actual)&&expected>0,"invalid own-phase comparison");const auto error=std::abs(actual-expected)/expected;anchoredFaceTimePeak=std::max(anchoredFaceTimePeak,error);anchoredTimeCheck(error<=1e-12L,"actual-site own-phase gross discrepancy");++anchoredFaceTimePositive;}
template<class F>void anchoredTimeReject(F&& f,const char* reason)
{bool rejected=false;try{f();}catch(const std::exception&){rejected=true;}anchoredTimeCheck(rejected,reason);++anchoredFaceTimeGuards;}
void coreAnchoredFaceTime()
{
    namespace FaceTime=LAMP::Geometry::PlanarFaceTime;namespace Cut=LAMP::Geometry::PlanarCellCut;
    const auto delta=std::ldexp(1.L,-(std::numeric_limits<long double>::digits+16));
    const AnchoredSites band{{{1,0},{0,0}},{{1,0},{delta,0}},{{1,0},{delta,1}},{{1,0},{0,1}}};
    const AnchoredSites triangle{{{1,0},{0,0}},{{1,0},{delta,0}},{{1,0},{0,1}}};
    for(bool triangular:{false,true})for(bool rz:{false,true})for(bool reversePolygon:{false,true})for(bool reverseFace:{false,true})for(int sx:{1,-1})for(int sy:{1,-1}){
        auto polygon=triangular?triangle:band;if(reversePolygon)std::reverse(polygon.begin(),polygon.end());Point a{0,1},b{2,1};if(reverseFace)std::swap(a,b);
        const FaceTime::Reflection reflection{sx,sy,{sx<0?2.L:0.L,sy<0?1.L:0.L}};
        auto original=arPoints(polygon);AnchoredRefPolygon physical;
        for(auto p:original)physical.push_back({arAdd(arMul({static_cast<long double>(sx),0},p.x),{reflection.translation.x,0}),arAdd(arMul({static_cast<long double>(sy),0},p.y),{reflection.translation.y,0})});
        const auto full=arValue(arMeasure(physical,{},rz)),half=arValue(arMeasure(arHalfPlane(physical,true),{},rz));
        const auto actual=FaceTime::integrateAnchoredPolygon(polygon,a,b,{0,1},{0,1},0,1,{},rz,512,reflection);
        anchoredTimeNear(actual.halfPositive,half);anchoredTimeNear(actual.fullPositive,full);
        anchoredTimeCheck(actual.halfNegative==0&&actual.fullNegative==0,"static axial negative gross not strict zero");anchoredFaceTimeZeros+=2;++anchoredFaceTimeCases;
        if(!triangular&&!rz&&!reversePolygon&&!reverseFace&&sx==1&&sy==1){std::vector<Point> collapsed;for(auto p:polygon)collapsed.push_back({p.anchor.x+p.offset.x,p.anchor.y+p.offset.y});
            const auto control=FaceTime::integratePolygon(collapsed,a,b,{0,1},{0,1},0,1,{},false,512,{});
            anchoredTimeCheck(control.halfPositive==0&&control.fullPositive==0,"old Point collapse negative control retained tiny band");}
    }
    anchoredTimeCheck(anchoredFaceTimeCases==64,"frozen anchored static Fubini coverage changed");
    const auto empty=FaceTime::integrateAnchoredPolygon({}, {0,1},{2,1},{0,1},{0,1},0,1,{},false);
    anchoredTimeCheck(empty.halfPositive==0&&empty.halfNegative==0&&empty.fullPositive==0&&empty.fullNegative==0,"empty polygon strict zero");anchoredFaceTimeZeros+=4;
    const AnchoredSites axis{{{0,0},{}},{{0,0},{delta,0}},{{0,1},{delta,0}},{{0,1},{}}};
    const auto axisGross=FaceTime::integrateAnchoredPolygon(axis,{0,0},{0,1},{1,0},{1,0},0,1,{},true);
    anchoredTimeCheck(axisGross.halfPositive==0&&axisGross.halfNegative==0&&axisGross.fullPositive==0&&axisGross.fullNegative==0,"RZ axis actual flux strict zero");anchoredFaceTimeZeros+=4;
    anchoredTimeReject([&]{FaceTime::integrateAnchoredPolygon(band,{0,1},{2,2},{0,1},{0,1},0,1,{},false);},"non-axis-aligned actual query accepted");
    anchoredTimeReject([&]{FaceTime::integrateAnchoredPolygon(band,{0,1},{2,1},{0,1},{0,1},1,1,{},false);},"nonpositive complete physical clock accepted");
    anchoredTimeReject([&]{FaceTime::integrateAnchoredPolygon(band,{0,1},{2,1},{0,1},{0,1},0,1,{},false,2);},"unsupported event budget accepted");
    const auto nan=std::numeric_limits<long double>::quiet_NaN();auto invalid=band;invalid[0].offset.x=nan;
    anchoredTimeReject([&]{FaceTime::integrateAnchoredPolygon(invalid,{0,1},{2,1},{0,1},{0,1},0,1,{},false);},"nonfinite public site offset accepted");
    anchoredTimeReject([&]{FaceTime::integrateAnchoredPolygon(band,{0,1},{2,1},{0,1},{0,1},0,1,{},false,512,{0,1,{}});},"unsupported reflection sign accepted");
    // Frozen candidate02 correction: images may cross the axis; the physical face clips their support.
    const AnchoredSites crossAxisImage{{{0,0},{-delta,0}},{{0,0},{delta,0}},{{0,1},{}}};
    for(bool rz:{false,true}){
        const auto square=arMul({delta,0},{delta,0});
        const auto full=rz?arMul({static_cast<long double>(CartesianMeasure::twoPi),0},arDiv(square,{6,0})):arDiv({delta,0},{2,0});
        const auto half=rz?arMul({static_cast<long double>(CartesianMeasure::twoPi),0},arDiv(square,{48,0})):arDiv({delta,0},{8,0});
        const auto image=FaceTime::integrateAnchoredPolygon(crossAxisImage,{0,1},{2,1},{0,1},{0,1},0,1,{},rz,512,{});
        anchoredTimeNear(image.halfPositive,arValue(half));anchoredTimeNear(image.fullPositive,arValue(full));
        anchoredTimeCheck(image.halfNegative==0&&image.fullNegative==0,"cross-axis image static negative gross not strict zero");
        anchoredFaceTimeZeros+=2;++anchoredFaceTimeImages;
    }
    anchoredTimeCheck(anchoredFaceTimeImages==2,"frozen cross-axis image XY/RZ coverage changed");
    anchoredTimeReject([&]{FaceTime::integrateAnchoredPolygon(band,{-delta,1},{2,1},{0,1},{0,1},0,1,{},true);},"negative Eulerian face radius accepted");
    std::cout<<std::setprecision(17)<<"anchored_face_time cases="<<anchoredFaceTimeCases<<" image_cases="<<anchoredFaceTimeImages<<" positive_checks="<<anchoredFaceTimePositive<<" exact_zero_checks="<<anchoredFaceTimeZeros
        <<" guards="<<anchoredFaceTimeGuards<<" peak_own_relative="<<anchoredFaceTimePeak<<" delta="<<delta<<" own_gate=1e-12 PASS\n";
}
}
// END V4a1a ADDED ANCHORED TEST BLOCK
// BEGIN V4a1b ORIGINAL NATIVE PIPELINE AND TRANSFER TEST BLOCK
unsigned v4bRealSlabs=0,v4bTransferCases=0,v4bNativeFailures=0;
std::uint64_t v4bNativeBits(double value){std::uint64_t result=0;std::memcpy(&result,&value,sizeof(result));return result;}
template<std::size_t N>void v4bOriginalBits(const std::vector<double>& actual,const std::array<std::uint64_t,N>& expected,const char* reason)
{check(actual.size()==N,reason);for(std::size_t i=0;i<N;++i)check(v4bNativeBits(actual[i])==expected[i],reason);}
std::vector<double> v4bOriginalEndAlpha()
{
    static constexpr std::array<std::uint64_t,40> data{{4607182418800017408ULL,4605337116763832501ULL,4605491442021919829ULL,4589752509994103044ULL,4588296838971833504ULL,0ULL,0ULL,0ULL,
        4599697512656935662ULL,0ULL,0ULL,4277057994932860420ULL,0ULL,0ULL,0ULL,0ULL,
        0ULL,0ULL,0ULL,0ULL,0ULL,0ULL,0ULL,0ULL,
        0ULL,0ULL,0ULL,0ULL,0ULL,0ULL,0ULL,0ULL,
        0ULL,0ULL,0ULL,0ULL,0ULL,0ULL,0ULL,0ULL}};std::vector<double> result(data.size());
    for(std::size_t i=0;i<data.size();++i)std::memcpy(&result[i],&data[i],sizeof(double));return result;
}
void originalAnchoredNativePipeline()
{
    stage="V4a1b original fixture4 real A/C native pipeline";auto mesh=vofMaterialMesh(4);const auto old=vofMaterialPublish(*mesh,true,.61L);
    static constexpr std::array<std::uint64_t,40> oldAlpha{{4607182418800017408ULL,4605475070154084071ULL,4604885030548504164ULL,4589996554519582192ULL,4583841426862459738ULL,0ULL,0ULL,0ULL,
        4600326453645682503ULL,0ULL,0ULL,0ULL,0ULL,0ULL,0ULL,0ULL,
        0ULL,0ULL,0ULL,0ULL,0ULL,0ULL,0ULL,0ULL,
        0ULL,0ULL,0ULL,0ULL,0ULL,0ULL,0ULL,0ULL,
        0ULL,0ULL,0ULL,0ULL,0ULL,0ULL,0ULL,0ULL}};
    static constexpr std::array<std::uint64_t,94> beginBits{{0ULL,4569403720218141126ULL,0ULL,4567790987397734753ULL,0ULL,13791163024252510561ULL,0ULL,13792775757072916934ULL,
        4566423257143376648ULL,4570658067967346082ULL,4559724548954897478ULL,4564026556979717178ULL,4563816565727060098ULL,4568131173227039222ULL,0ULL,0ULL,
        4551814903980710823ULL,4556083313405105390ULL,4569403720218141126ULL,4563197791502563903ULL,4563085796167813460ULL,4556968076549633886ULL,4567249209965879487ULL,4560950897870414824ULL,
        0ULL,0ULL,4554985759124551053ULL,4548992621544999671ULL,13786569828357339711ULL,13792775757072916934ULL,13780340113404409694ULL,13786457833022589268ULL,
        13784322934725190632ULL,13790621246820655295ULL,0ULL,0ULL,13772364658399775479ULL,13778357795979326861ULL,13794030104822121890ULL,13789795293998152456ULL,
        13787398593834492986ULL,13783096585809673286ULL,13791503210081815030ULL,13787188602581835906ULL,0ULL,0ULL,13779455350259881198ULL,13775186940835486631ULL,
        0ULL,0ULL,13797279356700287430ULL,13788540946248947500ULL,13800976589917454739ULL,13792208314043514691ULL,13797279356700287430ULL,13788540946248947500ULL,
        0ULL,0ULL,0ULL,0ULL,4558694191875193407ULL,4563580611919529052ULL,4551692953505093675ULL,4556717614255555625ULL,
        0ULL,0ULL,4563535500511965237ULL,4558407483818232274ULL,4556656212617482653ULL,4551302711983118800ULL,4562580572768828130ULL,4567424457211278760ULL,
        4561818861714730754ULL,4566785320198646120ULL,4567384358182333146ULL,4562325721162640457ULL,4566750076911486889ULL,4561370880375728984ULL,4558694191875193407ULL,4563580611919529052ULL,
        4561818861714730754ULL,4566785320198646120ULL,4563535500511965237ULL,4558407483818232274ULL,4566750076911486889ULL,4561370880375728984ULL,0ULL,0ULL,
        4551692953505093675ULL,4556717614255555625ULL,0ULL,0ULL,4556656212617482653ULL,4551302711983118800ULL}};
    static constexpr std::array<std::uint64_t,94> AEndBits{{0ULL,4565625850359953498ULL,0ULL,4566432216770156685ULL,0ULL,13786901334548201022ULL,0ULL,13790961432649959764ULL,
        4564228148582267974ULL,4569175249735250223ULL,4559724548954897478ULL,4564510376825839089ULL,4563044847873545330ULL,4567836660495187668ULL,0ULL,0ULL,
        4552775263976195868ULL,4557578731112360676ULL,4568960218692529375ULL,4564658210667709665ULL,4564113913340822523ULL,4559442320711326362ULL,4567575221385629605ULL,4563063747086284459ULL,
        0ULL,0ULL,4556978996094772050ULL,4552093212387565674ULL,13782478371561850835ULL,13791767799060162951ULL,13778566107301962681ULL,13787062607830241658ULL,
        13781683344679291131ULL,13790550374772883531ULL,0ULL,0ULL,13771867131896428304ULL,13779857133523298411ULL,13793810593966011023ULL,13790128187359416472ULL,
        13788829894212603640ULL,13785072183514671094ULL,13792096960348702923ULL,13788269007576756581ULL,0ULL,0ULL,13781535510837420555ULL,13777481712473241524ULL,
        0ULL,0ULL,13793501486842099802ULL,13789454828180511112ULL,13798408269674571891ULL,13793689639004480542ULL,13795465032277330260ULL,13790491052244007905ULL,
        0ULL,0ULL,0ULL,0ULL,13767766982691214604ULL,4561485115617174437ULL,13763252930165934074ULL,4553911367515501922ULL,
        0ULL,0ULL,4563450211757039899ULL,4559329429413897918ULL,4556139095696836860ULL,4552199047978904842ULL,13754076991273502299ULL,4566117000809804527ULL,
        13766448530384146701ULL,4564626941061177651ULL,4567832339521334916ULL,4563613518185157545ULL,4566913627790960192ULL,4562684712875238438ULL,4543367681399722044ULL,4563003326875860862ULL,
        4541060577503862927ULL,4565554720595644382ULL,4564628747279644557ULL,4560383081523230081ULL,4567373993229477636ULL,4563096295730446315ULL,0ULL,0ULL,
        4539255318808413123ULL,4556335964698908319ULL,0ULL,0ULL,4558094047643881661ULL,4553744738820181400ULL}};
    static constexpr std::array<std::uint64_t,94> CEndBits{{0ULL,13797279356700287430ULL,0ULL,13795666623879881057ULL,0ULL,4572294587025105249ULL,0ULL,4573907319845511622ULL,
        13794298893625522952ULL,13798533704449492386ULL,13787600185437043782ULL,13791902193461863482ULL,13791692202209206402ULL,13796006809709185526ULL,0ULL,0ULL,
        13779690540462857127ULL,13783958949887251694ULL,13797279356700287430ULL,13791073427984710207ULL,13790961432649959764ULL,13784843713031780190ULL,13795124846448025791ULL,13788826534352561128ULL,
        0ULL,0ULL,13782861395606697357ULL,13776868258027145975ULL,4567701391129934399ULL,4573907319845511622ULL,4561471676177004382ULL,4567589395795183956ULL,
        4565454497497785320ULL,4571752809593249983ULL,0ULL,0ULL,4553496221172370167ULL,4559489358751921549ULL,4575161667594716578ULL,4570926856770747144ULL,
        4568530156607087674ULL,4564228148582267974ULL,4572634772854409718ULL,4568320165354430594ULL,0ULL,0ULL,4560586913032475886ULL,4556318503608081319ULL,
        0ULL,0ULL,4578410919472882118ULL,4569672509021542188ULL,4582108152690049427ULL,4573339876816109379ULL,4578410919472882118ULL,4569672509021542188ULL,
        0ULL,0ULL,0ULL,0ULL,13786569828357339711ULL,13791456248401675356ULL,13779568589987239979ULL,13784593250737701929ULL,
        0ULL,0ULL,13791411136994111541ULL,13786283120300378578ULL,13784531849099628957ULL,13779178348465265104ULL,13790456209250974434ULL,13795300093693425064ULL,
        13789694498196877058ULL,13794660956680792424ULL,13795259994664479450ULL,13790201357644786761ULL,13794625713393633193ULL,13789246516857875288ULL,13786569828357339711ULL,13791456248401675356ULL,
        13789694498196877058ULL,13794660956680792424ULL,13791411136994111541ULL,13786283120300378578ULL,13794625713393633193ULL,13789246516857875288ULL,0ULL,0ULL,
        13779568589987239979ULL,13784593250737701929ULL,0ULL,0ULL,13784531849099628957ULL,13779178348465265104ULL}};
    v4bOriginalBits(old.alpha(),oldAlpha,"original fixture4 old alpha bits changed");
    const auto saved=old.alpha();const auto id=old.identity();const auto clock=old.timeLayer();const auto* alphaStorage=old.alpha().data();const auto* siteStorage=old.anchoredVertices().data();
    const auto geometry=HydroMomentumGeometry::build(*mesh,vofMaterialMetrics(*mesh,old));const auto oldMass=vofMaterialInspect(*mesh,geometry,old);const auto savedMass=oldMass.masses();
    const auto a=vofMaterialNormals(*mesh,old,0),b=vofMaterialNormals(*mesh,old,1),c=vofMaterialNormals(*mesh,old,2);
    v4bOriginalBits(a,beginBits,"original fixture4 native begin normal bits changed");v4bOriginalBits(b,AEndBits,"original fixture4 A native end normal bits changed");v4bOriginalBits(c,CEndBits,"original fixture4 C native end normal bits changed");
    const auto out=std::max(vofMaterialOutflow(*mesh,old,a,b),vofMaterialOutflow(*mesh,old,a,c));const auto dt=static_cast<double>(.2L/out);
    const VofClock next{8,.25,.25+dt,VofClock::Position::End};check(v4bNativeBits(next.end)==4618004176285829495ULL,"original fixture4 CFL-derived fullclock changed");
    const auto slabA=vofMaterialSlab(*mesh,old,geometry.layout(),next,a,b),slabC=vofMaterialSlab(*mesh,old,geometry.layout(),next,a,c);
    for(const auto* slab:{&slabA,&slabC}){
        try{const auto transport=VofTransport::prepare(*mesh,old,*slab,100,1);++v4bRealSlabs;
            check(transport.sourceIdentity()==id&&vofMaterialSameClock(transport.sourceTime(),clock),"real AP candidate did not retain one old source");
            auto mid=next;mid.position=Position::Midpoint;check(vofMaterialSameClock(transport.midpoint().timeLayer(),mid)&&vofMaterialSameClock(transport.endpoint().timeLayer(),next),"real AP candidates changed fullclock");
            vofMaterialInspect(*mesh,geometry,transport.midpoint());vofMaterialInspect(*mesh,geometry,transport.endpoint());
        }catch(const std::exception& error){++v4bNativeFailures;std::cerr<<"V4B_REAL_SLAB,"<<(slab==&slabA?"A":"C")<<",FAIL,"<<error.what()<<'\n';}
    }
    try{stage="V4a1b authenticated original C_end actual AP transfer";const auto alpha=v4bOriginalEndAlpha();
        const auto source=Publication::prepare(*mesh,old.domain(),next,alpha);check(v4bNativeBits(source.alpha()[11])==4277057994932860420ULL,"original C_end microphase stock replaced");
        const auto sourceId=source.identity();const auto* sourceSites=source.anchoredVertices().data();auto target=mesh->clone();
        const auto moved=GeometricVofMeshTransfer::prepare(*mesh,source,*target,domain(*target,true),next);
        check(moved.sourceIdentity()==sourceId&&sameClock(moved.candidate().timeLayer(),next)&&moved.candidate().alpha()==alpha,"same-lattice actual AP transfer changed raw original stock/clock");
        check(moved.overlaps().size()==40,"actual native same-lattice transfer overlap coverage changed");
        for(const auto& overlap:moved.overlaps()){
            check(overlap.kind==GeometricVofTransferOverlap::Kind::Unchanged&&overlap.sourceCell==overlap.targetCell,"same-lattice transfer remapped original cells");
            vofMaterialNear(overlap.liquidMeasured,vofMaterialReference(source,overlap.sourceCell,true),"original C_end AP transfer liquid");
            vofMaterialNear(overlap.gasMeasured,vofMaterialReference(source,overlap.sourceCell,false),"original C_end AP transfer gas");
        }
        check(source.identity()==sourceId&&source.alpha()==alpha&&source.anchoredVertices().data()==sourceSites,"AP transfer modified accepted source storage");++v4bTransferCases;
    }catch(const std::exception& error){++v4bNativeFailures;std::cerr<<"V4B_ORIGINAL_END_TRANSFER,FAIL,"<<error.what()<<'\n';}
    stage="V4a1b actual AP failures keep old owner";vofMaterialGuardContract(*mesh,geometry,old,oldMass);
    const auto foreign=Publication::prepare(*mesh,old.domain(),clock,saved);
    rejects([&]{Transport::prepare(*mesh,foreign,slabA,100,1);},"actual AP sameclock foreign old owner");
    auto bad=next;bad.begin=std::nextafter(bad.begin,1.);rejects([&]{vofMaterialSlab(*mesh,old,geometry.layout(),bad,a,b);},"actual AP begin one ULP differs");
    GeometricVofTransportConfig config;config.maximumOutflowCfl=.00001L;rejects([&]{Transport::prepare(*mesh,old,slabA,100,1,config);},"actual AP CFL failure");
    check(old.identity()==id&&old.alpha()==saved&&old.alpha().data()==alphaStorage&&old.anchoredVertices().data()==siteStorage&&sameClock(old.timeLayer(),clock)
        &&oldMass.masses()==savedMass&&mesh->cellIndexNum()==40,"AP failure or candidate replaced old owner/M/native mesh");
    std::cout<<std::setprecision(17)<<"V4B_REAL_NATIVE_PIPELINE,"<<v4bRealSlabs<<",2,"<<vofMaterialLayers<<",5,"<<v4bTransferCases<<",1,"<<v4bNativeFailures
        <<','<<vofMaterialComparisons<<','<<vofMaterialGuards<<','<<vofMaterialPeakPhase<<','<<vofMaterialPeakMass<<','<<vofMaterialPeakPartition<<','<<vofMaterialPeakInventory<<'\n';
    check(v4bRealSlabs==2&&vofMaterialLayers==5&&v4bTransferCases==1&&v4bNativeFailures==0,"actual original A/C five-layer/transfer AP coverage FAIL");
}
// END V4a1b ORIGINAL NATIVE PIPELINE AND TRANSFER TEST BLOCK

// BEGIN V4b1 RETAINED NATIVE SUBFACE FORMAL TEST BLOCK
namespace {
unsigned v4b1Cases=0,v4b1BudgetCases=0,v4b1Queries=0,v4b1Whole=0,v4b1Zeros=0,v4b1Positive=0,v4b1Guards=0,v4b1RawReferences=0;
long double v4b1PeakPartition=0,v4b1PeakReference=0,v4b1PeakRational=0;
std::array<long double,12> v4b1Fields(const GeometricVofIntegratedFlux& f)
{return {{f.total,f.positiveTotal,f.negativeTotal,f.rawLiquid,f.liquid,f.gas,f.mass,f.positiveLiquid,f.negativeLiquid,f.positiveGas,f.negativeGas,f.limiter}};}
void v4b1Exact(const GeometricVofIntegratedFlux& a,const GeometricVofIntegratedFlux& b,const char* why)
{const auto x=v4b1Fields(a),y=v4b1Fields(b);for(std::size_t i=0;i<x.size();++i)check(x[i]==y[i]&&std::signbit(x[i])==std::signbit(y[i]),std::string(why)+" field="+std::to_string(i));++v4b1Whole;}
void v4b1Zero(const GeometricVofIntegratedFlux& f,const char* why)
{const auto x=v4b1Fields(f);for(std::size_t i=0;i+1<x.size();++i)check(x[i]==0,why);++v4b1Zeros;}
void v4b1Own(long double actual,long double expected,const char* why,long double& peak)
{check(std::isfinite(actual)&&std::isfinite(expected),std::string(why)+" nonfinite");if(expected==0){check(actual==0,std::string(why)+" theoreticalzero");++v4b1Zeros;return;}const auto e=std::abs(actual-expected)/std::abs(expected);peak=std::max(peak,e);near(actual,expected,std::abs(expected),why,1e-12L);++v4b1Positive;}
void v4b1Partition(const GeometricVofIntegratedFlux& parent,const GeometricVofIntegratedFlux& a,const GeometricVofIntegratedFlux& b)
{
    const auto sum=[](long double x,long double y){return vofMaterialValue(vofMaterialSum(x,y));};
    for(const auto& r:std::array<std::array<long double,3>,6>{{{parent.positiveTotal,a.positiveTotal,b.positiveTotal},{parent.negativeTotal,a.negativeTotal,b.negativeTotal},
        {parent.positiveLiquid,a.positiveLiquid,b.positiveLiquid},{parent.negativeLiquid,a.negativeLiquid,b.negativeLiquid},{parent.positiveGas,a.positiveGas,b.positiveGas},{parent.negativeGas,a.negativeGas,b.negativeGas}}})
        v4b1Own(sum(r[1],r[2]),r[0],"retained ownphase partition",v4b1PeakPartition);
    for(const auto& r:std::array<std::array<long double,4>,4>{{{parent.total,a.total,b.total,parent.positiveTotal+parent.negativeTotal},
        {parent.liquid,a.liquid,b.liquid,parent.positiveLiquid+parent.negativeLiquid},{parent.gas,a.gas,b.gas,parent.positiveGas+parent.negativeGas},
        {parent.mass,a.mass,b.mass,100*(parent.positiveLiquid+parent.negativeLiquid)+parent.positiveGas+parent.negativeGas}}}){
        if(r[3]==0)check(r[0]==0&&r[1]==0&&r[2]==0,"zero gross acquired signed restriction");else near(sum(r[1],r[2]),r[0],r[3],"retained signed gross-scaled partition",1e-12L);
    }
}
long double v4b1NormalPrimitive(long double a,long double b,long double dt,long double area,long double extent,bool positive)
{
    std::vector<long double> cuts{0,extent};if(b!=a){const auto root=-a/(b-a);if(root>0&&root<extent)cuts.push_back(root);}std::sort(cuts.begin(),cuts.end());long double result=0;
    for(std::size_t i=1;i<cuts.size();++i){const auto lo=cuts[i-1],hi=cuts[i];if((a+(b-a)*(lo+hi)/2>0)!=positive)continue;
        const auto q=dt*area*(a*(hi-lo)+(b-a)*(hi*hi-lo*lo)/2);result+=positive?q:-q;}return result;
}
long double v4b1WorldArea(Point a,Point b,bool rz)
{const auto length=a.x==b.x?std::abs(b.y-a.y):std::abs(b.x-a.x);return length*(rz?twoPi*(a.x==b.x?a.x:(a.x+b.x)/2):1.L);}
std::array<long double,4> v4b1IndependentGross(const QuadtreeMesh& mesh,const Publication& p,const Slab& slab,int row,Point worldA,Point worldB,bool liquid)
{
    using W=VofMaterialWide;using WP=std::array<W,2>;
    const auto add=vofMaterialAdd,sub=vofMaterialSub,mul=vofMaterialMul,div=vofMaterialDiv;
    const auto compare=[&](W a,W b){return vofMaterialSign(sub(a,b));};const auto half=[&](W a){return mul(a,{.5L,0});};
    const auto native=p.domain().face(mesh,row);if(native.horizontal?worldA.x>worldB.x:worldA.y>worldB.y)std::swap(worldA,worldB);
    const auto origin=p.frameOrigin();const Point d{worldB.x-worldA.x,worldB.y-worldA.y};const WP a{{sub({worldA.x,0},{origin.x,0}),sub({worldA.y,0},{origin.y,0})}};
    const auto v0=slab.beginVelocity()[row],v1=slab.endVelocity()[row];const auto normal=slab.beginNormal()[row];
    check(v1.x==-2*v0.x&&v1.y==-2*v0.y&&slab.endNormal()[row]==-2*normal,"fixed physical distance reference requires exact end=-2begin tuple");
    if(normal==0)return {};const W dt=sub({slab.endTime().end,0},{slab.endTime().begin,0});const auto dtValue=vofMaterialValue(dt);
    const auto dx=dtValue*(std::abs(v0.x)+std::abs(v1.x)),dy=dtValue*(std::abs(v0.y)+std::abs(v1.y));
    check(worldA.x-dx>mesh.minCoord.x&&worldB.x+dx<mesh.maxCoord.x&&worldA.y-dy>mesh.minCoord.y&&worldB.y+dy<mesh.maxCoord.y,"independent subface oracle requires interior sameold support");
    std::vector<std::vector<WP>> polygons;for(const auto& c:p.cells()){const auto box=bounds(mesh,c.nativeCell);if(box.x1<worldA.x-dx||box.x0>worldB.x+dx||box.y1<worldA.y-dy||box.y0>worldB.y+dy)continue;
        for(const auto& poly:p.anchoredPhasePolygons(c.nativeCell,liquid)){std::vector<WP> sites;for(auto site:poly)sites.push_back(v4bWidePoint(site));polygons.push_back(std::move(sites));}}
    // s=dt(tau-1.5tau^2): departure=a-v0*s, |normal(t)|dt d tau=|normal0| |ds|.
    // Each AP event is the linear rational N/D. Equality below compares the
    // exact supplied coefficient words, so no rounded-root or epsilon merging.
    struct Event{W n,d,x;};
    const auto event=[&](W n,W denominator){check(vofMaterialSign(denominator)!=0,"independent linear AP event zero denominator");if(vofMaterialSign(denominator)<0){n=vofMaterialNeg(n);denominator=vofMaterialNeg(denominator);}return Event{n,denominator,div(n,denominator)};};
    const auto order=[&](const Event& a,const Event& b){std::vector<long double> expansion;
        const auto term=[&](long double x){long double q=x;std::vector<long double> next;for(auto word:expansion){const auto sum=vofMaterialSum(q,word);if(sum.low!=0)next.push_back(sum.low);q=sum.high;}if(q!=0)next.push_back(q);expansion=std::move(next);};
        const auto product=[&](long double x,long double y,int sign){const auto value=x*y;check(std::isfinite(value)&&(value!=0||x==0||y==0),"linear event equality product out of reference range");term(sign*std::fma(x,y,-value));term(sign*value);};
        for(auto x:{a.n.high,a.n.low})for(auto y:{b.d.high,b.d.low})product(x,y,1);for(auto x:{b.n.high,b.n.low})for(auto y:{a.d.high,a.d.low})product(x,y,-1);
        return expansion.empty()?0:(expansion.back()>0?1:-1);};
    std::vector<Event> geometryEvents;const auto retain=[&](W numerator,W denominator){if(vofMaterialSign(denominator)!=0)geometryEvents.push_back(event(numerator,denominator));};
    const auto cross=[&](W ex,W ey,W x,W y){return sub(mul(ex,y),mul(ey,x));};
    for(const auto& poly:polygons){if(poly.empty())continue;auto previous=poly.back();for(auto current:poly){const int axis=native.horizontal?1:0;
            retain(sub(a[axis],current[axis]),{axis?v0.y:v0.x,0});
            auto low=previous,high=current;const auto xorder=compare(low[0],high[0]);if(xorder>0||(xorder==0&&compare(low[1],high[1])>0))std::swap(low,high);
            const auto ex=sub(high[0],low[0]),ey=sub(high[1],low[1]);
            // A polygon edge parallel to the face has exactly its vertex-normal
            // event already retained above, including both moving endpoints.
            if(vofMaterialSign(native.horizontal?ey:ex)!=0)for(int endpoint=0;endpoint<2;++endpoint){const auto x=sub(add(a[0],{endpoint*d.x,0}),low[0]),y=sub(add(a[1],{endpoint*d.y,0}),low[1]);
                retain(cross(ex,ey,x,y),cross(ex,ey,{v0.x,0},{v0.y,0}));}previous=current;
        }}
    const auto interval=[&](const std::vector<WP>& poly,const WP& departure,W& lo,W& hi){
        lo={0,0};hi={1,0};if(poly.size()<3)return false;W area{};auto previous=poly.back();for(auto current:poly){area=add(area,cross(previous[0],previous[1],current[0],current[1]));previous=current;}
        const int orientation=vofMaterialSign(area)>=0?1:-1;previous=poly.back();for(auto current:poly){const auto ex=sub(current[0],previous[0]),ey=sub(current[1],previous[1]);
            const auto c=mul({static_cast<long double>(orientation),0},cross(ex,ey,sub(departure[0],previous[0]),sub(departure[1],previous[1])));
            const auto slope=mul({static_cast<long double>(orientation),0},cross(ex,ey,{d.x,0},{d.y,0}));const auto sign=vofMaterialSign(slope);
            if(sign>0){const auto root=div(vofMaterialNeg(c),slope);if(compare(root,lo)>0)lo=root;}
            else if(sign<0){const auto root=div(vofMaterialNeg(c),slope);if(compare(root,hi)<0)hi=root;}
            else if(vofMaterialSign(c)<0)return false;if(compare(hi,lo)<=0)return false;previous=current;
        }return compare(hi,lo)>0;
    };
    const auto trace=[&](W distance){const WP departure{{sub(a[0],mul({v0.x,0},distance)),sub(a[1],mul({v0.y,0},distance))}};
        std::vector<std::pair<W,W>> intervals;for(const auto& poly:polygons){W lo,hi;if(interval(poly,departure,lo,hi))intervals.emplace_back(lo,hi);}
        std::sort(intervals.begin(),intervals.end(),[&](const auto& a,const auto& b){const auto c=compare(a.first,b.first);return c<0||(c==0&&compare(a.second,b.second)<0);});
        std::vector<std::pair<W,W>> merged;for(auto entry:intervals){if(!merged.empty()&&compare(entry.first,merged.back().second)<=0){if(compare(entry.second,merged.back().second)>0)merged.back().second=entry.second;}else merged.push_back(entry);}
        W trace{};for(auto entry:merged){auto length=mul({native.horizontal?d.x:d.y,0},sub(entry.second,entry.first));if(p.domain().axisymmetric()){
                const auto radius=add({worldA.x,0},mul({d.x,0},half(add(entry.first,entry.second))));length=mul(length,mul({twoPi,0},radius));}trace=add(trace,length);}
        return trace;
    };
    // On a fixed active-edge s panel, XY trace is linear and RZ trace quadratic.
    // GL3 integrates both directly; the original two-half-panel self gate stays.
    const W node{std::sqrt(3.L/5.L),0},edgeWeight=div({5,0},{9,0}),centerWeight=div({8,0},{9,0});
    std::array<long double,4> result{};for(int extentIndex=0;extentIndex<2;++extentIndex)for(int sign=0;sign<2;++sign){
        const bool forward=(normal>0)==(sign==0);const auto lo=forward?event({0,0},{1,0}):(extentIndex?event(vofMaterialNeg(dt),{2,0}):event(dt,{8,0}));const auto hi=event(dt,{6,0});
        std::vector<Event> events{lo,hi};for(const auto& e:geometryEvents)if(order(e,lo)>0&&order(e,hi)<0)events.push_back(e);
        std::sort(events.begin(),events.end(),[&](const Event& a,const Event& b){return order(a,b)<0;});events.erase(std::unique(events.begin(),events.end(),[&](const Event& a,const Event& b){return order(a,b)==0;}),events.end());
        const auto panel=[&](W lo,W hi){const auto m=half(add(lo,hi)),h=half(sub(hi,lo)),z=mul(node,h),left=sub(m,z),right=add(m,z);
            check(compare(m,lo)>0&&compare(m,hi)<0&&compare(left,lo)>0&&compare(right,hi)<0,"independent AP distance panel not representable; distinct rational tail cannot be merged");
            return mul(h,add(mul(edgeWeight,add(trace(left),trace(right))),mul(centerWeight,trace(m))));};
        W coarseSum{},fineSum{};for(std::size_t k=1;k<events.size();++k){const auto a=events[k-1].x,b=events[k].x,m=half(add(a,b));coarseSum=add(coarseSum,panel(a,b));fineSum=add(fineSum,add(panel(a,m),panel(m,b)));}
        const auto coarse=vofMaterialValue(mul({std::abs(normal),0},coarseSum)),fine=vofMaterialValue(mul({std::abs(normal),0},fineSum));if(fine>0)near(fine,coarse,fine,"independent rawphase ownscale refinement",2e-13L);else check(coarse==0,"independent rawphase unresolved zero");result[2*extentIndex+sign]=fine;
    }return result;
}
void v4b1SealedRationals()
{
    using Cell=LAMP::Interface::Detail::GeometricVofPhaseBudgetCell;
    const auto makeFace=[](long double q,long double pl,long double pg,long double nl,long double ng){GeometricVofTransportFace f;f.nativeRow=0;f.negativeCell=0;f.positiveCell=1;
        f.half.positiveTotal=f.half.negativeTotal=q;f.half.positiveLiquid=pl;f.half.positiveGas=pg;f.half.negativeLiquid=nl;f.half.negativeGas=ng;f.half.rawLiquid=pl-nl;f.full=f.half;return f;};
    for(int scale:{0,-512})for(bool half:{true,false}){
        const long double q=std::ldexp(1.L,-10+scale),micro=std::ldexp(1.L,-100+scale);
        {
            stage="V4b1 sealed unit microcounterflow";const std::vector<Cell> cells{{std::ldexp(1.L,scale),.5},{std::ldexp(1.L,scale),.5}};const std::vector<GeometricVofTransportFace> faces{makeFace(q,q,micro,q,2*micro)};
            const auto b=LAMP::Interface::Detail::prepareGeometricVofPhaseBudget(cells,faces,half);check(b.restriction&&b.flux[0].limiter==1&&b.limitedFaceExtents==0,"unit exactoperator binding");
            const auto first=makeFace(q/2,q/2,micro,q/2,2*micro),second=makeFace(q/2,q/2,0,q/2,0);const auto x=b.restriction->measure(0,half?first.half:first.full),y=b.restriction->measure(0,half?second.half:second.full);
            v4b1Own(x.positiveGas,micro,"unit child positivegas",v4b1PeakRational);v4b1Own(x.negativeGas,2*micro,"unit child negativegas",v4b1PeakRational);v4b1Own(x.gas,-micro,"unit child signedgas",v4b1PeakRational);
            v4b1Own(x.liquid,micro/2,"unit firstchild small signedliquid",v4b1PeakRational);v4b1Own(y.liquid,micro/2,"unit secondchild small signedliquid",v4b1PeakRational);
            v4b1Own(y.positiveGas,0,"unit high truezero positivegas",v4b1PeakRational);v4b1Own(y.negativeGas,0,"unit high truezero negativegas",v4b1PeakRational);
            v4b1Partition(b.flux[0],x,y);++v4b1BudgetCases;
        }
        {
            stage="V4b1 sealed nearunit diagnostic1";const long double epsilon=std::ldexp(1.L,-std::numeric_limits<long double>::digits-8),delta=q*epsilon;
            const std::vector<Cell> cells{{std::ldexp(1.L,scale),0},{std::ldexp(1.L,scale),1}};const std::vector<GeometricVofTransportFace> faces{makeFace(q,q,0,q,delta)};
            const auto b=LAMP::Interface::Detail::prepareGeometricVofPhaseBudget(cells,faces,half);check(b.restriction&&b.flux[0].limiter==1&&b.limitedFaceExtents==1,"nearunit factor lost exact complement despite diagnostic1");
            const auto first=makeFace(q/2,q/2,0,q/2,delta),second=makeFace(q/2,q/2,0,q/2,0);const auto x=b.restriction->measure(0,half?first.half:first.full),y=b.restriction->measure(0,half?second.half:second.full);
            const auto denominator=vofMaterialSum(1,epsilon);const auto minority=vofMaterialMul({q,0},vofMaterialDiv({epsilon,0},denominator));
            const auto tiny=vofMaterialMul({q/2,0},vofMaterialDiv(vofMaterialMul({epsilon,0},{epsilon,0}),denominator));
            const auto mg=vofMaterialValue(minority),sl=vofMaterialValue(tiny);check(mg>0&&sl>0,"independent nearunit positive rational unrepresentable");
            v4b1Own(x.positiveGas,mg/2,"nearunit complement positivegas",v4b1PeakRational);v4b1Own(x.negativeGas,mg,"nearunit negativegas",v4b1PeakRational);
            v4b1Own(x.gas,-mg/2,"nearunit signedgas",v4b1PeakRational);v4b1Own(x.liquid,-sl,"nearunit negative sub-ULP signedliquid",v4b1PeakRational);v4b1Own(y.liquid,sl,"nearunit positive sub-ULP signedliquid",v4b1PeakRational);
            check((1-b.flux[0].limiter)*(q/2)==0&&x.positiveGas>0,"rounded limiter negativecontrol did not erase truepositive minor");v4b1Partition(b.flux[0],x,y);++v4b1BudgetCases;
        }
        {
            stage="V4b1 sealed saturated zero";const std::vector<Cell> cells{{std::ldexp(1.L,scale),0},{std::ldexp(1.L,scale),0}};const std::vector<GeometricVofTransportFace> faces{makeFace(q,q,0,0,q)};
            const auto b=LAMP::Interface::Detail::prepareGeometricVofPhaseBudget(cells,faces,half);check(b.restriction&&b.flux[0].limiter==0,"saturated exactzero operator missing");const auto piece=makeFace(q/2,q/2,0,0,q/2);const auto x=b.restriction->measure(0,half?piece.half:piece.full);
            v4b1Own(x.positiveLiquid,0,"saturated accepted positive liquid zero",v4b1PeakRational);v4b1Own(x.negativeLiquid,0,"saturated accepted negative liquid zero",v4b1PeakRational);v4b1Own(x.liquid,0,"saturated signed liquid zero",v4b1PeakRational);v4b1Own(x.positiveGas,q/2,"saturated positive gas",v4b1PeakRational);v4b1Own(x.negativeGas,q/2,"saturated negative gas",v4b1PeakRational);++v4b1BudgetCases;
        }
    }
    check(v4b1BudgetCases==12,"V4b1 exact rational coverage changed");
}
void v4b1NativeRestriction()
{
    v4b1SealedRationals();
    for(int fixture=0;fixture<8;++fixture){const bool translated=fixture>=6,rz=translated||fixture/3!=0;const double alpha=translated?(fixture-6):std::array<double,3>{{0,1,.5}}[fixture%3];
        stage="V4b1 actual retained native subface "+std::to_string(fixture);
        auto mesh=vofMaterialMesh(4);if(translated){DictionaryMesh d;d.Nx=d.Ny=4;d.minCoord={2,3};d.maxCoord={3,4};d.maxCellNum=4096;d.maxInitialIterations=0;d.maxLevel_local={2};DictionaryMesh::RefineBox box;box.level=1;box.min={2.5,3};box.max={3,4};d.refineBoxes={box};mesh=std::make_unique<QuadtreeMesh>(d);}
        std::vector<double> inputAlpha(mesh->cellIndexNum(),alpha);
        if(alpha==.5){
            rejects([&]{Publication::prepare(*mesh,domain(*mesh,rz),{7,.125,.25,Position::End},inputAlpha);},"uniform mixed zero-gradient orientation remains unsupported");++v4b1Guards;
            for(int cell=0;cell<mesh->cellIndexNum();++cell)inputAlpha[cell]=.4+.2*(mesh->cells(mesh->cellIndexToLabel[cell]).centerCoord().x-mesh->minCoord.x);
            check(std::all_of(inputAlpha.begin(),inputAlpha.end(),[](double x){return x>=.4&&x<=.6;}),"bounded nonzero-gradient source profile changed");
        }
        const auto old=Publication::prepare(*mesh,domain(*mesh,rz),{7,.125,.25,Position::End},inputAlpha);
        check(old.frameOrigin().y==(translated?3.5L:.5L)&&old.frameOrigin().x==(translated?2.L:0.L),"actual source nonzero frame coverage changed");
        std::vector<double> a(old.domain().nativeFaceCount()),b(a.size());for(int row=0;row<old.domain().nativeFaceCount();++row){const auto f=old.domain().face(*mesh,row);auto e=faceEndpoints(*mesh,f);
            const Vector2D p{static_cast<double>(e[0].x-mesh->minCoord.x),static_cast<double>(e[0].y-mesh->minCoord.y)},q{static_cast<double>(e[1].x-mesh->minCoord.x),static_cast<double>(e[1].y-mesh->minCoord.y)};
            const auto integral=(f.horizontal?-1:1)*(vofMaterialPsi(q,0)-vofMaterialPsi(p,0));a[row]=f.area==0||f.outwardSign!=0?0:static_cast<double>(integral/f.area);b[row]=-2*a[row];}
        const auto dt=static_cast<double>(.05L/vofMaterialOutflow(*mesh,old,a,b));const Clock next{8,.25,.25+dt,Position::End};const auto slab=Slab::prepare(*mesh,old,next,a,b);const auto t=Transport::prepare(*mesh,old,slab,100,1),copy=t;
        if(alpha==.5){ // Independent exact-unit proof: low budgets >=.4(C-O), demand <=O+I.
            std::vector<long double> incoming(old.cells().size()),outgoing(old.cells().size());
            for(int row=0;row<old.domain().nativeFaceCount();++row){const auto f=old.domain().face(*mesh,row);const auto e=faceEndpoints(*mesh,f);const auto area=v4b1WorldArea(e[0],e[1],rz);
                const auto positive=v4b1NormalPrimitive(a[row],b[row],next.end-next.begin,area,1,true),negative=v4b1NormalPrimitive(a[row],b[row],next.end-next.begin,area,1,false);
                if(f.negativeCell>=0){outgoing[f.negativeCell]+=positive;incoming[f.negativeCell]+=negative;}if(f.positiveCell>=0){incoming[f.positiveCell]+=positive;outgoing[f.positiveCell]+=negative;}}
            const auto eps=std::numeric_limits<long double>::epsilon();for(const auto& c:old.cells()){
                check(old.alpha()[c.nativeCell]>=.4&&old.alpha()[c.nativeCell]<=.6,"actual source alpha leaves independent limiter bound");
                check(outgoing[c.nativeCell]<=.05L*c.capacity*(1+512*eps),"independent full gross CFL bound");
                near(incoming[c.nativeCell],outgoing[c.nativeCell],c.capacity,"independent closed discrete-divergence premise",512*eps);
                check(.4L*(c.capacity-outgoing[c.nativeCell])>outgoing[c.nativeCell]+incoming[c.nativeCell],"independent lowphase budget does not strictly dominate correction demand");}
        }
        const auto savedAlpha=old.alpha();const auto* savedStorage=old.alpha().data();const auto savedFaces=t.faces();const auto oldId=old.identity();
        check(t.sameInstance(copy)&&t.advector().sameInstance(slab)&&slab.sameInstance(t.advector())&&&t.source().alpha()==&old.alpha(),"retained exact source/slab storage detached");
        check(t.liquidDensity()==100&&t.gasDensity()==1&&sameClock(t.sourceTime(),old.timeLayer())&&sameClock(t.endTime(),next),"retained density/fullclock lost");t.validateFor(*mesh,old,slab,100,1);
        check(t.restrictNativeFaces(*mesh,{}).empty(),"empty queryvector changed coverage");
        for(int row:{9,19,60,70}){const auto native=old.domain().face(*mesh,row);const auto e=faceEndpoints(*mesh,native);const Point mid{(e[0].x+e[1].x)/2,(e[0].y+e[1].y)/2};
            const std::vector<GeometricVofNativeSubface> queries{{row,e[0],e[1]},{row,e[1],e[0]},{row,e[0],mid},{row,mid,e[1]},{row,e[0],e[0]}};GeometricVofTransportReport diagnostic;const auto r=t.restrictNativeFaces(*mesh,queries,&diagnostic);v4b1Queries+=static_cast<unsigned>(queries.size());
            check(r.size()==5,"subface batch output coverage");for(int i=0;i<5;++i)check(r[i].nativeRow==row&&r[i].negativeCell==native.negativeCell&&r[i].positiveCell==native.positiveCell,"subface changed native incidence");
            v4b1Exact(r[0].half,t.faces()[row].half,"cached whole half");v4b1Exact(r[0].full,t.faces()[row].full,"cached whole full");v4b1Exact(r[1].half,r[0].half,"reversed canonical whole half");v4b1Exact(r[1].full,r[0].full,"reversed canonical whole full");v4b1Zero(r[4].half,"empty half");v4b1Zero(r[4].full,"empty full");
            v4b1Partition(t.faces()[row].half,r[2].half,r[3].half);v4b1Partition(t.faces()[row].full,r[2].full,r[3].full);
            check(diagnostic.maximumCoverageResidual<=1e-12L,"subface original normal/coverage budget exceeded");
            if(alpha!=.5){for(int part=0;part<2;++part){const auto begin=part?mid:e[0],end=part?e[1]:mid;for(int extent=0;extent<2;++extent){const auto& f=extent?r[part+2].full:r[part+2].half;const auto area=v4b1WorldArea(begin,end,rz);
                const auto positive=v4b1NormalPrimitive(a[row],b[row],next.end-next.begin,area,extent?1.L:.5L,true),negative=v4b1NormalPrimitive(a[row],b[row],next.end-next.begin,area,extent?1.L:.5L,false);
                v4b1Own(f.positiveLiquid,alpha?positive:0,"pure child positive liquid",v4b1PeakReference);v4b1Own(f.negativeLiquid,alpha?negative:0,"pure child negative liquid",v4b1PeakReference);v4b1Own(f.positiveGas,alpha?0:positive,"pure child positive gas",v4b1PeakReference);v4b1Own(f.negativeGas,alpha?0:negative,"pure child negative gas",v4b1PeakReference);
                near(f.mass,100*f.liquid+f.gas,100*(f.positiveLiquid+f.negativeLiquid)+f.positiveGas+f.negativeGas,"restricted same density mass",1e-12L);}}}
            else if(row==70){check(t.report().maximumOutflowCfl<=.05L*(1+512*std::numeric_limits<long double>::epsilon()),"bounded nonzero-gradient exactunit a priori CFL proof failed");
                const auto parentL=v4b1IndependentGross(*mesh,old,slab,row,e[0],e[1],true),parentG=v4b1IndependentGross(*mesh,old,slab,row,e[0],e[1],false);
                for(int k=0;k<4;++k)check(parentL[k]>0&&parentG[k]>0,"fixed row70 independent parent lost genuine two-phase support");
                for(int part=0;part<2;++part){const auto childL=v4b1IndependentGross(*mesh,old,slab,row,part?mid:e[0],part?e[1]:mid,true),childG=v4b1IndependentGross(*mesh,old,slab,row,part?mid:e[0],part?e[1]:mid,false);
                    for(int extent=0;extent<2;++extent){const auto& parent=extent?t.faces()[row].full:t.faces()[row].half;const auto& f=extent?r[part+2].full:r[part+2].half;
                        for(int sign=0;sign<2;++sign){const auto k=2*extent+sign;const auto expected=[&](long double h,long double child,long double whole){if(whole==0){check(h==0&&child==0,"independent zero rawparent phase acquired child");return 0.L;}return vofMaterialValue(vofMaterialMul({h,0},vofMaterialDiv({child,0},{whole,0})));};
                            v4b1Own(sign?f.negativeLiquid:f.positiveLiquid,expected(sign?parent.negativeLiquid:parent.positiveLiquid,childL[k],parentL[k]),"independent actualAP calibrated liquid",v4b1PeakReference);
                            v4b1Own(sign?f.negativeGas:f.positiveGas,expected(sign?parent.negativeGas:parent.positiveGas,childG[k],parentG[k]),"independent actualAP calibrated gas",v4b1PeakReference);++v4b1RawReferences;}
                    }
                }
            }
        }
        if(rz&&!translated){const auto e=faceEndpoints(*mesh,old.domain().face(*mesh,0));const Point mid{e[0].x,(e[0].y+e[1].y)/2};const auto r=t.restrictNativeFaces(*mesh,{{0,e[0],mid}});v4b1Zero(r[0].half,"RZ zeroaxis half");v4b1Zero(r[0].full,"RZ zeroaxis full");++v4b1Queries;}
        if(fixture==0){const auto e=faceEndpoints(*mesh,old.domain().face(*mesh,19));const auto reject=[&](auto f,const char* why){rejects(f,why);++v4b1Guards;};
            reject([&]{t.restrictNativeFaces(*mesh,{{-1,e[0],e[1]}});},"restriction invalid native row");reject([&]{t.restrictNativeFaces(*mesh,{{old.domain().nativeFaceCount(),e[0],e[1]}});},"restriction past native row");
            auto bad=e[0];bad.x+=.125L;reject([&]{t.restrictNativeFaces(*mesh,{{19,bad,e[1]}});},"restriction off native plane");bad=e[1];bad.y+=.125L;reject([&]{t.restrictNativeFaces(*mesh,{{19,e[0],bad}});},"restriction outside native interval");bad=e[0];bad.x=std::numeric_limits<long double>::quiet_NaN();reject([&]{t.restrictNativeFaces(*mesh,{{19,bad,e[1]}});},"restriction nonfinite coordinate");
            const auto foreign=Publication::prepare(*mesh,old.domain(),old.timeLayer(),old.alpha());reject([&]{t.validateFor(*mesh,foreign,slab,100,1);},"restriction foreign sameclock source");const auto otherSlab=Slab::prepare(*mesh,old,next,a,b);reject([&]{t.validateFor(*mesh,old,otherSlab,100,1);},"restriction samevalue foreign slab");reject([&]{t.validateFor(*mesh,old,slab,101,1);},"restriction liquid density mismatch");reject([&]{t.validateFor(*mesh,old,slab,100,2);},"restriction gas density mismatch");
            const auto tailEndpoints=faceEndpoints(*mesh,old.domain().face(*mesh,8));
            check(tailEndpoints[0].x==.5L&&tailEndpoints[1].x==.5L&&tailEndpoints[0].y==0&&tailEndpoints[1].y==.125L,"fixed row8 representability guard geometry changed");
            auto tail=tailEndpoints[0];tail.y=std::ldexp(1.L,-80);check(tail.y>tailEndpoints[0].y&&tail.y<tailEndpoints[1].y,"world tail must remain within actual native segment");
            reject([&]{t.restrictNativeFaces(*mesh,{{8,tail,tailEndpoints[1]}});},"restriction exact world-to-frame tail cannot be represented");
            GeometricVofTransportReport sentinel;sentinel.maximumCoverageResidual=.375L;reject([&]{t.restrictNativeFaces(*mesh,{{19,e[0],e[1]},{-1,e[0],e[1]}},&sentinel);},"restriction partial batch failure");check(sentinel.maximumCoverageResidual==.375L,"failed restriction changed caller diagnostic");
        }
        check(old.alpha()==savedAlpha&&old.alpha().data()==savedStorage&&old.identity()==oldId&&sameClock(old.timeLayer(),{7,.125,.25,Position::End}),"native restriction changed old source/storage/clock");
        for(std::size_t i=0;i<savedFaces.size();++i){v4b1Exact(t.faces()[i].half,savedFaces[i].half,"query changed parent half payload");v4b1Exact(t.faces()[i].full,savedFaces[i].full,"query changed parent full payload");}++v4b1Cases;
    }
    check(v4b1Cases==8&&v4b1BudgetCases==12&&v4b1Guards==13&&v4b1RawReferences==16,"V4b1 frozen formal coverage changed");
    std::cout<<std::setprecision(17)<<"V4B1_NATIVE_SUBFACE_SUMMARY,"<<v4b1Cases<<','<<v4b1BudgetCases<<','<<v4b1Queries<<','<<v4b1Whole<<','<<v4b1Zeros<<','<<v4b1Positive<<','<<v4b1Guards<<','<<v4b1RawReferences<<','<<v4b1PeakPartition<<','<<v4b1PeakReference<<','<<v4b1PeakRational<<",PASS\n";
}
}
// END V4b1 RETAINED NATIVE SUBFACE FORMAL TEST BLOCK

// BEGIN V4b2 SMALL ANALYTIC INTERIOR PRODUCER TEST BLOCK
namespace {
unsigned v4b2ProducerCases=0,v4b2ProducerGuards=0,v4b2ProducerQueries=0,v4b2ProducerPositive=0,v4b2ProducerZeros=0;
long double v4b2ProducerPeak=0,v4b2ProducerVelocityPeak=0;
std::array<Point,2> v4b2IndependentP1(const QuadtreeMesh& mesh,const Publication& source,const Slab& slab,int nativeCell,Point query)
{
    const auto box=bounds(mesh,nativeCell);const long double cx=(box.x0+box.x1)/2,cy=(box.y0+box.y1)/2,h=std::min(box.x1-box.x0,box.y1-box.y0);
    std::array<std::array<long double,3>,3> gram{};std::array<std::array<long double,3>,4> rhs{};
    for(int row=0;row<source.domain().nativeFaceCount();++row){const auto f=source.domain().face(mesh,row);if(f.negativeCell!=nativeCell&&f.positiveCell!=nativeCell)continue;
        Point sample{f.center.x,f.center.y};if(source.domain().axisymmetric()&&f.horizontal){const auto e=faceEndpoints(mesh,f);const auto a=e[0].x,b=e[1].x;sample.x=2*(a*a+a*b+b*b)/(3*(a+b));}
        const long double x=(sample.x-cx)/h,y=(sample.y-cy)/h,weight=1/(1+x*x+y*y);const std::array<long double,3> design{{1,x,y}};
        const auto a=slab.beginVelocity()[row],b=slab.endVelocity()[row];const std::array<long double,4> values{{a.x,a.y,b.x,b.y}};
        for(int i=0;i<3;++i){for(int j=0;j<3;++j)gram[i][j]+=weight*design[i]*design[j];for(int k=0;k<4;++k)rhs[k][i]+=weight*design[i]*values[k];}
    }
    const auto determinant=[](const std::array<std::array<long double,3>,3>& a){return a[0][0]*(a[1][1]*a[2][2]-a[1][2]*a[2][1])-a[0][1]*(a[1][0]*a[2][2]-a[1][2]*a[2][0])+a[0][2]*(a[1][0]*a[2][1]-a[1][1]*a[2][0]);};
    const auto det=determinant(gram);check(det>0&&std::isfinite(det),"independent internal P1 moment matrix singular");const std::array<long double,3> at{{1,(query.x-cx)/h,(query.y-cy)/h}};std::array<long double,4> value{};
    for(int k=0;k<4;++k)for(int j=0;j<3;++j){auto replaced=gram;for(int i=0;i<3;++i)replaced[i][j]=rhs[k][i];value[k]+=at[j]*determinant(replaced)/det;}
    return {{{value[0],value[1]},{value[2],value[3]}}};
}
void v4b2ProducerNear(long double actual,long double expected,const char* why)
{
    check(std::isfinite(actual)&&std::isfinite(expected),why);if(expected==0){check(actual==0,why);++v4b2ProducerZeros;return;}
    const auto error=std::abs(actual-expected)/std::abs(expected);v4b2ProducerPeak=std::max(v4b2ProducerPeak,error);near(actual,expected,std::abs(expected),why,1e-12L);++v4b2ProducerPositive;
}
void v4b2InteriorProducer()
{
    for(bool rz:{false,true}){stage="V4b2 pure analytic interior producer "+std::to_string(rz);auto mesh=vofMaterialMesh(4);const double alpha=rz?0:1;
        const auto source=Publication::prepare(*mesh,domain(*mesh,rz),{7,.125,.25,Position::End},std::vector<double>(mesh->cellIndexNum(),alpha));
        const auto a=vofMaterialNormals(*mesh,source,0),b=vofMaterialNormals(*mesh,source,2);const auto dt=static_cast<double>(.05L/vofMaterialOutflow(*mesh,source,a,b));const Clock next{8,.25,.25+dt,Position::End};
        const auto slab=Slab::prepare(*mesh,source,next,a,b);const auto accepted=Transport::prepare(*mesh,source,slab,100,1);const auto saved=source.alpha();const auto* storage=source.alpha().data();
        const auto box=bounds(*mesh,0);const Point first{(box.x0+box.x1)/2,box.y0+(box.y1-box.y0)/4},last{first.x,box.y0+3*(box.y1-box.y0)/4},centroid{first.x,(first.y+last.y)/2};
        const auto expected=v4b2IndependentP1(*mesh,source,slab,0,centroid);GeometricVofTransportReport report;
        const auto result=accepted.measureInteriorFaces(*mesh,{{0,first,last},{0,last,first}},&report);v4b2ProducerQueries+=2;check(result.size()==2&&result[0].nativeCell==0&&result[1].nativeCell==0,"internal producer query incidence changed");
        for(int k=0;k<2;++k){const auto actual=k?result[0].endVelocity:result[0].beginVelocity;const auto scale=std::max(std::abs(expected[k].x),std::abs(expected[k].y));check(scale>0,"analytic internal velocity fixture is trivial");
            for(const auto values:std::array<std::array<long double,2>,2>{{{actual.x,expected[k].x},{actual.y,expected[k].y}}}){const auto error=std::abs(values[0]-values[1])/scale;v4b2ProducerVelocityPeak=std::max(v4b2ProducerVelocityPeak,error);near(values[0],values[1],scale,"independent P1 internal centroid velocity",1e-12L);}}
        v4b1Exact(result[0].half,result[1].half,"interior reversed canonical half");v4b1Exact(result[0].full,result[1].full,"interior reversed canonical full");
        for(int extent=0;extent<2;++extent){const auto& f=extent?result[0].full:result[0].half;const auto area=v4b1WorldArea(first,last,rz);const auto positive=v4b1NormalPrimitive(expected[0].x,expected[1].x,next.end-next.begin,area,extent?1.L:.5L,true),negative=v4b1NormalPrimitive(expected[0].x,expected[1].x,next.end-next.begin,area,extent?1.L:.5L,false),total=positive-negative;
            v4b2ProducerNear(f.positiveTotal,positive,"analytic internal positive normal");v4b2ProducerNear(f.negativeTotal,negative,"analytic internal negative normal");v4b2ProducerNear(f.total,total,"analytic internal signed normal");
            v4b2ProducerNear(f.positiveLiquid,alpha?positive:0,"analytic pure internal positive L");v4b2ProducerNear(f.negativeLiquid,alpha?negative:0,"analytic pure internal negative L");v4b2ProducerNear(f.positiveGas,alpha?0:positive,"analytic pure internal positive G");v4b2ProducerNear(f.negativeGas,alpha?0:negative,"analytic pure internal negative G");
            v4b2ProducerNear(f.liquid,alpha?total:0,"analytic pure internal signed L");v4b2ProducerNear(f.gas,alpha?0:total,"analytic pure internal signed G");v4b2ProducerNear(f.rawLiquid,f.liquid,"internal raw equals accepted high L");v4b2ProducerNear(f.mass,(alpha?100.L:1.L)*total,"analytic pure internal mass");check(f.limiter==1,"interior unit extension changed");}
        const auto reject=[&](auto action,const char* why){rejects(action,why);++v4b2ProducerGuards;};reject([&]{accepted.measureInteriorFaces(*mesh,{{-1,first,last}});},"internal invalid native cell");
        auto boundary=first;boundary.x=box.x0;reject([&]{accepted.measureInteriorFaces(*mesh,{{0,boundary,{boundary.x,last.y}}});},"internal boundary is not owner interior");auto diagonal=last;diagonal.x+=.03125L;reject([&]{accepted.measureInteriorFaces(*mesh,{{0,first,diagonal}});},"internal diagonal support");
        GeometricVofTransportReport sentinel;sentinel.maximumCoverageResidual=.375L;reject([&]{accepted.measureInteriorFaces(*mesh,{{0,first,first}},&sentinel);},"internal degenerate support");check(sentinel.maximumCoverageResidual==.375L,"invalid internal query changed caller diagnostic");
        check(source.alpha()==saved&&source.alpha().data()==storage&&accepted.sourceIdentity()==source.identity()&&sameClock(source.timeLayer(),{7,.125,.25,Position::End}),"internal query changed accepted old storage/clock");
        std::cout<<std::setprecision(17)<<"V4B2_INTERIOR_PRODUCER,"<<rz<<','<<alpha<<','<<result[0].beginVelocity.x<<','<<result[0].beginVelocity.y<<','<<v4b2ProducerPeak<<",PASS\n";++v4b2ProducerCases;
    }
    check(v4b2ProducerCases==2&&v4b2ProducerGuards==8&&v4b2ProducerQueries==4,"V4b2 producer frozen coverage changed");
    std::cout<<"V4B2_INTERIOR_PRODUCER_SUMMARY,"<<v4b2ProducerCases<<','<<v4b2ProducerQueries<<','<<v4b2ProducerGuards<<','<<v4b2ProducerPositive<<','<<v4b2ProducerZeros<<','<<v4b2ProducerVelocityPeak<<','<<v4b2ProducerPeak<<",PASS\n";
}
}
// END V4b2 SMALL ANALYTIC INTERIOR PRODUCER TEST BLOCK

// BEGIN L1 LOCAL PLIC SWEEP AND COST CONTRACT
namespace {
unsigned l1Cases=0,l1LedgerCells=0,l1Fluxes=0,l1Guards=0,l1Analytic=0,l1CoarseFine=0;long double l1PeakLedger=0,l1PeakConservation=0,l1PeakStrip=0;
long double l1Cap(Vector2D lo,Vector2D hi,bool rz){return vofMaterialCapacity(lo,hi,rz);}
std::vector<double> l1PlaneAlpha(const QuadtreeMesh& mesh,bool rz,int axis)
{std::vector<double> alpha(mesh.cellIndexNum());for(int cell=0;cell<mesh.cellIndexNum();++cell){const auto b=bounds(mesh,cell);Vector2D lo{static_cast<double>(b.x0),static_cast<double>(b.y0)},hi{static_cast<double>(b.x1),static_cast<double>(b.y1)};const auto c=l1Cap(lo,hi,rz);if(axis==0)hi.x=std::max(lo.x,std::min(hi.x,.73));else hi.y=std::max(lo.y,std::min(hi.y,.73));alpha[cell]=static_cast<double>(l1Cap(lo,hi,rz)/c);}return alpha;}
bool l1SymmetricHalo(const QuadtreeMesh& mesh,int cell)
{const auto b=bounds(mesh,cell);const auto dx=b.x1-b.x0,dy=b.y1-b.y0,x=(b.x0+b.x1)/2,y=(b.y0+b.y1)/2;for(int i=-1;i<=1;++i)for(int j=-1;j<=1;++j){bool found=false;for(int k=0;k<mesh.cellIndexNum();++k){const auto q=bounds(mesh,k);if(q.x1-q.x0==dx&&q.y1-q.y0==dy&&(q.x0+q.x1)/2==x+i*dx&&(q.y0+q.y1)/2==y+j*dy)found=true;}if(!found)return false;}return true;}
void l1Near(long double actual,long double expected,long double c,const char* why,long double& peak)
{near(actual,expected,c,why,1e-12L);peak=std::max(peak,std::abs(actual-expected)/c);}
void l1Inspect(const QuadtreeMesh& mesh,const Publication& source,const Slab& slab,const Transport& t,int axis)
{
    check(t.sourceIdentity()==source.identity()&&t.source().plicCells().data()==source.plicCells().data()&&t.advector().sameInstance(slab),"local sweep detached original owner/slab");const auto dt=slab.endTime().end-slab.endTime().begin;const bool rz=source.domain().axisymmetric();auto mid=slab.endTime();mid.position=Position::Midpoint;
    check(sameClock(t.midpoint().timeLayer(),mid)&&sameClock(t.endpoint().timeLayer(),slab.endTime()),"local half/full complete clock changed");check(t.report().localDirectionalSweeps==4&&t.report().maximumLocalColumnCfl<=.5L&&t.report().maximumIncomingCfl<=.5L&&t.report().maximumCompressionResidual<=1e-12L,"local sweep work/CFL/compression contract");
    for(int extent=0;extent<2;++extent){const auto& target=extent?t.endpoint():t.midpoint();check(target.geometryKind()==GeometricVofGeometryKind::LocalPlic,"local output changed representation");std::vector<std::array<long double,2>> net(source.plicCells().size());std::vector<long double> totalNet(net.size());long double oldL=0,newL=0,oldG=0,newG=0,totalCapacity=0;
        for(std::size_t row=0;row<t.faces().size();++row){const auto native=source.domain().face(mesh,static_cast<int>(row));const auto& f=extent?t.faces()[row].full:t.faces()[row].half;const auto a=slab.beginNormal()[row],b=slab.endNormal()[row];const auto mean=extent?(static_cast<long double>(a)+b)/2:(3*static_cast<long double>(a)+b)/4;const auto expected=mean*native.area*(extent?dt:dt/2);
            long double scale=0;for(int cell:{native.negativeCell,native.positiveCell})if(cell>=0)scale+=source.cells()[cell].capacity;check(scale>0,"local face has no native capacity");l1Near(f.total,expected,scale,"frozen local temporal mean",l1PeakLedger);l1Near(f.liquid+f.gas,f.total,scale,"shared local phase F closure",l1PeakLedger);check(f.mass==100*f.liquid+f.gas,"same rho-weighted local material F");
            if(native.area==0)check(f.total==0&&f.liquid==0&&f.gas==0&&f.mass==0,"local zero-axis flux is not strictzero");for(int cell:{native.negativeCell,native.positiveCell})if(cell>=0){const auto sign=cell==native.negativeCell?1.L:-1.L;net[cell][0]+=sign*f.liquid;net[cell][1]+=sign*f.gas;totalNet[cell]+=sign*f.total;}if(native.negativeCell>=0&&native.positiveCell>=0){check(f.liquid+(-f.liquid)==0&&f.gas+(-f.gas)==0,"actual shared local opposite incidence");if(mesh.cells(mesh.cellIndexToLabel[native.negativeCell]).level()!=mesh.cells(mesh.cellIndexToLabel[native.positiveCell]).level())++l1CoarseFine;}++l1Fluxes;
            // The first direction uses the actual OLD plane; qualify exact
            // symmetric native support before applying an axial strip oracle.
            if((native.horizontal?1:0)==axis&&f.total!=0){const auto e=faceEndpoints(mesh,native);const auto coordinate=axis?e[0].y:e[0].x;const auto tangentLo=axis?e[0].x:e[0].y,tangentHi=axis?e[1].x:e[1].y;if(coordinate!=.75L||tangentLo<.25L||tangentHi>.75L)continue;
                const int donor=f.total>0?native.negativeCell:native.positiveCell;check(donor>=0,"analytic local donor outside");const auto& p=source.plicCells()[donor];const bool mixed=p.inventory[0]>0&&p.inventory[1]>0;if(mixed&&!l1SymmetricHalo(mesh,donor))continue;if(mixed){check(p.resolved&&(axis?p.normal.y:p.normal.x)>0&&std::abs(axis?p.normal.x:p.normal.y)<=1e-10,"qualified strip source is not axial");}
                auto lo=p.rectangle.lower,hi=p.rectangle.upper;if(axis){lo.x=static_cast<double>(e[0].x);hi.x=static_cast<double>(e[1].x);const auto section=(static_cast<long double>(hi.x)-lo.x)*(rz?twoPi*(static_cast<long double>(hi.x)+lo.x)/2:1.L),depth=std::abs(f.total)/section;if(f.total>0)lo.y=static_cast<double>(static_cast<long double>(hi.y)-depth);else hi.y=static_cast<double>(static_cast<long double>(lo.y)+depth);}
                else{lo.y=static_cast<double>(e[0].y);hi.y=static_cast<double>(e[1].y);const auto length=static_cast<long double>(hi.y)-lo.y;if(rz){const auto radial=2*std::abs(f.total)/(twoPi*length);if(f.total>0)lo.x=static_cast<double>(std::sqrt(static_cast<long double>(hi.x)*hi.x-radial));else hi.x=static_cast<double>(std::sqrt(static_cast<long double>(lo.x)*lo.x+radial));}else{const auto depth=std::abs(f.total)/length;if(f.total>0)lo.x=static_cast<double>(static_cast<long double>(hi.x)-depth);else hi.x=static_cast<double>(static_cast<long double>(lo.x)+depth);}}
                auto lhi=hi,glo=lo;if(axis){lhi.y=std::max(lo.y,std::min(hi.y,.73));glo.y=std::min(hi.y,std::max(lo.y,.73));}else{lhi.x=std::max(lo.x,std::min(hi.x,.73));glo.x=std::min(hi.x,std::max(lo.x,.73));}const auto sign=f.total>0?1.L:-1.L;
                l1Near(f.liquid,sign*l1Cap(lo,lhi,rz),source.cells()[donor].capacity,"independent first-axis liquid strip",l1PeakStrip);l1Near(f.gas,sign*l1Cap(glo,hi,rz),source.cells()[donor].capacity,"independent first-axis gas strip",l1PeakStrip);++l1Analytic;
            }
        }
        for(std::size_t cell=0;cell<net.size();++cell){const auto& old=source.plicCells()[cell].inventory;const auto& now=target.plicCells()[cell].inventory;const auto c=source.cells()[cell].capacity,cc=old[0]>old[1]?1.L:0.L;for(unsigned phase=0;phase<2;++phase){check(now[phase]>=0&&now[phase]<=c,"local post-update stock out of bounds; no clipping permitted");const auto expected=old[phase]-net[cell][phase]+(phase?1-cc:cc)*totalNet[cell];l1Near(now[phase],expected,c,"actual local stock/phase-flux/compression ledger",l1PeakLedger);}l1Near(now[0]+now[1],c,c,"local physical pair capacity",l1PeakLedger);oldL+=old[0];oldG+=old[1];newL+=now[0];newG+=now[1];totalCapacity+=c;++l1LedgerCells;}
        l1Near(newL,oldL,totalCapacity,"local closed liquid conservation",l1PeakConservation);l1Near(newG,oldG,totalCapacity,"local closed gas conservation",l1PeakConservation);
    }
}
// Independent boundary-only Gauss references for LocalPlic Coons extension.
// Read the retained native P0 payloads; do not call its interior integrator.
std::array<long double,12> l1BoundaryPart(const QuadtreeMesh& mesh,const Publication& source,
    const Transport& transport,int cell,int axis,int side,long double lo,long double hi,int extent)
{
    std::array<long double,12> result{};const auto box=bounds(mesh,cell);const bool horizontal=axis==1;
    const auto coordinate=axis?(side?box.y1:box.y0):(side?box.x1:box.x0);
    for(int row=0;row<source.domain().nativeFaceCount();++row){const auto face=source.domain().face(mesh,row);
        if(face.horizontal!=horizontal||(face.negativeCell!=cell&&face.positiveCell!=cell)
            ||(axis?face.center.y:face.center.x)!=coordinate)continue;
        const auto ends=faceEndpoints(mesh,face);const auto a=horizontal?ends[0].x:ends[0].y,b=horizontal?ends[1].x:ends[1].y;
        const auto first=std::max(lo,std::min(a,b)),last=std::min(hi,std::max(a,b));if(last<=first)continue;
        const auto area=[&](long double x,long double y){return (y-x)*(source.domain().axisymmetric()&&horizontal?(x+y)/2:1.L);};
        const auto ratio=area(first,last)/area(std::min(a,b),std::max(a,b));
        const auto values=v4b1Fields(extent?transport.faces()[row].full:transport.faces()[row].half);
        for(std::size_t k=0;k+1<values.size();++k)result[k]+=ratio*values[k];
    }return result;
}
void l1CoonsPayload(const GeometricVofIntegratedFlux& flux)
{
    const auto values=v4b1Fields(flux);for(auto value:values)check(std::isfinite(value),"nonfinite Coons payload");
    for(std::size_t k:{1U,2U,7U,8U,9U,10U})check(values[k]>=0,"Coons discarded an original phase source by a negative gross");
    check(flux.mass==100*flux.liquid+flux.gas,"Coons mass detached from its own signed phase currents");
    near(flux.positiveTotal,flux.positiveLiquid+flux.positiveGas,std::max(1.L,flux.positiveTotal),"Coons positive phase gross addback",1e-12L);
    near(flux.negativeTotal,flux.negativeLiquid+flux.negativeGas,std::max(1.L,flux.negativeTotal),"Coons negative phase gross addback",1e-12L);
    const auto net=[&](long double actual,long double positive,long double negative){const auto gross=positive+negative;
        if(gross==0)check(actual==0,"Coons exact zero phase source became nonzero");
        else near(actual,positive-negative,gross,"Coons signed current and routed source gross",1e-12L);};
    net(flux.liquid,flux.positiveLiquid,flux.negativeLiquid);net(flux.gas,flux.positiveGas,flux.negativeGas);
    const auto gross=flux.positiveTotal+flux.negativeTotal;
    if(gross==0)check(flux.total==0,"Coons exact zero volume source became nonzero");
    else near(flux.total,flux.liquid+flux.gas,gross,"Coons signed physical phase closure",1e-12L);
}
void l1CoonsQuadrants(const QuadtreeMesh& mesh,const Publication& source,const Transport& transport)
{
    int cell=-1;
    for(int row=0;row<source.domain().nativeFaceCount()&&cell<0;++row){const auto face=source.domain().face(mesh,row);
        if(face.negativeCell<0||face.positiveCell<0)continue;
        const int a=mesh.cells(mesh.cellIndexToLabel[face.negativeCell]).level(),b=mesh.cells(mesh.cellIndexToLabel[face.positiveCell]).level();
        if(a!=b)cell=a<b?face.negativeCell:face.positiveCell;}
    check(cell>=0,"existing scene lacks its real coarse leaf beside a hanging subface");
    const auto box=bounds(mesh,cell);const bool rz=source.domain().axisymmetric();
    const std::array<long double,3> x{{box.x0,(box.x0+box.x1)/2,box.x1}},y{{box.y0,(box.y0+box.y1)/2,box.y1}};
    std::vector<GeometricVofInteriorSubface> queries;
    for(int j=0;j<2;++j)queries.push_back({cell,{x[1],y[j]},{x[1],y[j+1]}});
    for(int i=0;i<2;++i)queries.push_back({cell,{x[i],y[1]},{x[i+1],y[1]}});
    const auto count=queries.size();for(std::size_t i=0;i<count;++i)queries.push_back({cell,queries[i].end,queries[i].begin});
    const auto interior=transport.measureInteriorFaces(mesh,queries);long double peak=0,oldWitness=0;
    for(int extent=0;extent<2;++extent){std::array<std::array<long double,12>,4> face;
        for(std::size_t i=0;i<count;++i){const auto& f=extent?interior[i].full:interior[i].half;
            l1CoonsPayload(f);face[i]=v4b1Fields(f);const auto reversed=v4b1Fields(extent?interior[i+count].full:interior[i+count].half);
            for(std::size_t k=0;k<face[i].size();++k)check(face[i][k]==reversed[k],"Coons actual internal face reversal changed shared payload");}
        const auto west=l1BoundaryPart(mesh,source,transport,cell,0,0,y[0],y[2],extent),east=l1BoundaryPart(mesh,source,transport,cell,0,1,y[0],y[2],extent);
        const auto south=l1BoundaryPart(mesh,source,transport,cell,1,0,x[0],x[2],extent),north=l1BoundaryPart(mesh,source,transport,cell,1,1,x[0],x[2],extent);
        for(int i=0;i<2;++i)for(int j=0;j<2;++j){
            const auto left=i?face[j]:l1BoundaryPart(mesh,source,transport,cell,0,0,y[j],y[j+1],extent);
            const auto right=i?l1BoundaryPart(mesh,source,transport,cell,0,1,y[j],y[j+1],extent):face[j];
            const auto bottom=j?face[2+i]:l1BoundaryPart(mesh,source,transport,cell,1,0,x[i],x[i+1],extent);
            const auto top=j?l1BoundaryPart(mesh,source,transport,cell,1,1,x[i],x[i+1],extent):face[2+i];
            const auto fraction=capacity({x[i],x[i+1],y[j],y[j+1]},rz)/capacity(box,rz);
            for(std::size_t k:{0U,3U,4U,5U,6U}){
                const auto out=right[k]-left[k]+top[k]-bottom[k],target=fraction*(east[k]-west[k]+north[k]-south[k]);
                const auto gross=std::abs(right[k])+std::abs(left[k])+std::abs(top[k])+std::abs(bottom[k])+std::abs(target);
                if(gross==0)check(out==0,"Coons zero quadrant source became nonzero");
                else{near(out,target,gross,"independent physical quadrant Gauss balance",512*std::numeric_limits<double>::epsilon());peak=std::max(peak,std::abs(out-target)/gross);}
            }
        }
        // The retired LocalPlic opposite-side average has a real nonzero
        // quadrant defect on this SAME hanging-face leaf and source velocity.
        const auto w=l1BoundaryPart(mesh,source,transport,cell,0,0,y[0],y[1],extent),e=l1BoundaryPart(mesh,source,transport,cell,0,1,y[0],y[1],extent);
        const auto s=l1BoundaryPart(mesh,source,transport,cell,1,0,x[0],x[1],extent),n=l1BoundaryPart(mesh,source,transport,cell,1,1,x[0],x[1],extent);
        const auto xi=capacity({x[0],x[1],y[0],y[2]},rz)/capacity(box,rz),fraction=capacity({x[0],x[1],y[0],y[1]},rz)/capacity(box,rz);
        const auto oldVertical=(1-xi)*w[0]+xi*e[0],oldHorizontal=.5L*(s[0]+n[0]);
        const auto oldOut=oldVertical-w[0]+oldHorizontal-s[0],target=fraction*(east[0]-west[0]+north[0]-south[0]);
        const auto gross=std::abs(oldVertical)+std::abs(w[0])+std::abs(oldHorizontal)+std::abs(s[0]);
        if(gross>0)oldWitness=std::max(oldWitness,std::abs(oldOut-target)/gross);
    }
    check(oldWitness>1e-8L,"existing coarse-leaf source failed to distinguish old noncommuting interior extension");
    std::cout<<"LOCAL_PLIC_COONS_QUADRANTS,"<<rz<<','<<cell<<','<<peak<<','<<oldWitness<<",PASS\n";
}


// Independent conservative-VOF qualification. It only observes the accepted
// source, public shared F and actual endpoint masses; it never rewrites R/F/M.
namespace ConservativeVofMomentumContract {
namespace X=LAMP::Geometry::Detail::FiniteVolume;
constexpr long double gate=256.L*std::numeric_limits<double>::epsilon();
unsigned qualified=0,rows=0,coarseFine=0,axes=0,finiteRows=0,negativeControls=0,constantCases=0;
long double peakCanonical=0,peakStored=0,peakCommutation=0;
bool admissible(const X::Dyadic& residual,const X::Dyadic& scale)
{return residual.absolute().compare(scale.multiplied(gate))<=0;}
long double ratio(const X::Dyadic& residual,const X::Dyadic& scale)
{if(scale.compare(X::Dyadic{})==0){check(residual.compare(X::Dyadic{})==0,"conservative zero-scale nonzero residual");return 0;}
 return residual.absolute().rounded()/scale.rounded();}
void qualify(const QuadtreeMesh& mesh,const HydroDualMaterial& dual)
{
    const auto& m=dual.material();const auto& g=dual.geometry();const auto& l=dual.layout();const auto& d=l.domain();
    check(m.vofTransport().source().geometryKind()==GeometricVofGeometryKind::LocalPlic&&dual.beginMass().representation()==HydroMomentumMass::DensityRepresentation::CellAverageProjection,
        "conservative qualification did not use real LocalPlic P0");
    check(dual.usesConservativeVofMomentum(),"real LocalPlic P0 did not select conservative momentum");
    std::vector<X::Dyadic> primal(l.pressureCount()),primalGross(l.pressureCount());
    for(int u=0;u<l.velocityCount();++u){const auto f=d.face(mesh,l.nativeFace(u));const auto q=m.faceIntegrals()[u];
        if(f.negativeCell>=0){primal[l.pressureDof(f.negativeCell)].add(q);primalGross[l.pressureDof(f.negativeCell)].add(std::abs(q));}
        if(f.positiveCell>=0){primal[l.pressureDof(f.positiveCell)].add(q,-1);primalGross[l.pressureDof(f.positiveCell)].add(std::abs(q));}
        if(f.area==0)++axes;
        if(f.negativeCell>=0&&f.positiveCell>=0&&mesh.cells(mesh.cellIndexToLabel[f.negativeCell]).level()!=mesh.cells(mesh.cellIndexToLabel[f.positiveCell]).level())++coarseFine;
    }
    long double ownCanonical=0,ownStored=0;
    for(int u=0;u<l.velocityCount();++u){X::Dyadic outward,gross,mapped,mappedGross,scale;
        for(int side=0;side<4;++side)for(int id:g.volumes()[u].sides[side]){const auto q=dual.faceIntegrals()[id];outward.add(q,side%2?1:-1);gross.add(std::abs(q));}
        for(const auto& o:g.volumes()[u].overlaps){const auto cell=l.nativeCell(o.pressure);const auto box=bounds(mesh,cell);
            const auto physical=capacity({o.lower.x,o.upper.x,o.lower.y,o.upper.y},d.axisymmetric())/capacity(box,d.axisymmetric());
            mapped.add(primal[o.pressure].multiplied(physical));mappedGross.add(primalGross[o.pressure].multiplied(physical));
            const auto component=l.normalAxis(u)==HydroDofLayout::NormalAxis::First?0:1;
            if(dual.beginMass().usesFiniteP0Partition(o.pressure,component,0)||dual.beginMass().usesFiniteP0Partition(o.pressure,component,1)
                ||dual.endMass().usesFiniteP0Partition(o.pressure,component,0)||dual.endMass().usesFiniteP0Partition(o.pressure,component,1))++finiteRows;
        }
        for(const auto* mass:{&dual.beginMass(),&dual.endMass()}){const auto phase=mass->representedPhaseVolumes()[u];
            scale.add(X::Dyadic::product({m.liquidDensity(),phase.liquid}));scale.add(X::Dyadic::product({m.gasDensity(),phase.gas}));}
        scale.add(gross);
        X::Dyadic canonical(dual.cells()[u].massChange);canonical.add(outward);
        X::Dyadic stored(dual.endMass().masses()[u]);stored.add(dual.beginMass().masses()[u],-1);stored.add(outward);
        auto commutator=outward;commutator.add(mapped,-1);auto commutationScale=scale;commutationScale.add(mappedGross);
        check(admissible(canonical,scale),"LocalPlic P0 canonical DeltaM+same-F outflow exceeds256eps endpoint-mass budget");
        check(admissible(stored,scale),"actual stored M1-M0+same-F outflow exceeds256eps endpoint-mass budget");
        check(admissible(commutator,commutationScale),"Coons partial-dual incidence fails independently weighted primal incidence");
        const auto reported=dual.cells()[u].integratedDefect;
        near(reported,canonical.rounded(),scale.rounded(),"canonical R is the observed actual residual, not a forced zero",gate);
        ownCanonical=std::max(ownCanonical,ratio(canonical,scale));ownStored=std::max(ownStored,ratio(stored,scale));
        peakCommutation=std::max(peakCommutation,ratio(commutator,commutationScale));
        if(negativeControls==0&&scale.compare(X::Dyadic{})>0){auto invalid=stored;invalid.add(scale.multiplied(1e-6L));
            check(!admissible(invalid,scale),"local continuity gate accepted manufactured finite imbalance");++negativeControls;}
        ++rows;
    }
    const auto& r=dual.report();
    check(std::isfinite(r.maxRelativeContinuityDefect)&&r.maxRelativeContinuityDefect>=0&&r.maxRelativeContinuityDefect<=gate
        &&std::isfinite(r.maxRelativeStoredMassContinuityDefect)&&r.maxRelativeStoredMassContinuityDefect>=0&&r.maxRelativeStoredMassContinuityDefect<=gate,
        "conservative continuity qualification report missing or outside original roundoff gate");
    for(const auto value:{r.maxRelativePrimalContinuityDefect,r.maxRelativeMassTransferDefect,r.maxRelativeFluxCommutatorDefect})
        check(std::isfinite(value)&&value>=0&&value<=gate,"independent primal/mass-transfer/commutator qualification report missing");
    // The independent gates above remain authoritative even if a native
    // diagnostic ratio rounds a sufficiently tiny exact numerator to zero.
    peakCanonical=std::max(peakCanonical,ownCanonical);peakStored=std::max(peakStored,ownStored);++qualified;
}
}

void l1InventoryProjectionPair(QuadtreeMesh& mesh,const Publication& source,const Transport& transport,
    const HydroMomentumGeometry& geometry,const HydroMaterialTransport& material,
    const HydroPhaseGeometry& before,const HydroPhaseGeometry& middle,const HydroPhaseGeometry& after)
{
    const bool rz=source.domain().axisymmetric();const auto& layout=geometry.layout();
    const std::array<const HydroPhaseGeometry*,3> phases{{&before,&middle,&after}};
    const std::array<HydroMomentumMass,3> projected{{
        HydroMomentumMass::buildCellAverageProjection(mesh,geometry,before,100,1),
        HydroMomentumMass::buildCellAverageProjection(mesh,geometry,middle,100,1),
        HydroMomentumMass::buildCellAverageProjection(mesh,geometry,after,100,1)}};
    rejects([&]{(void)projected[0].phaseVolumes();},"P0 inventory mass cannot publish a sharp phase certificate");
    rejects([&]{(void)projected[0].overlapPhaseVolumes();},"P0 inventory mass cannot publish sharp overlap certificates");
    for(std::size_t layer=0;layer<phases.size();++layer){const auto& mass=projected[layer];const auto& phase=*phases[layer];
        check(mass.representation()==HydroMomentumMass::DensityRepresentation::CellAverageProjection
            &&mass.model()==HydroPhaseGeometry::Model::GeometricVof&&mass.vofIdentity()==phase.vofIdentity(),"P0 inventory mass changed representation or authentic VOF owner");
        for(int u=0;u<layout.velocityCount();++u){std::array<long double,2> expected{};
            const auto& volume=geometry.volumes()[u];const auto& represented=mass.representedOverlapPhaseVolumes()[u];
            check(represented.size()==volume.overlaps.size(),"P0 inventory overlap support changed");
            for(std::size_t k=0;k<volume.overlaps.size();++k){const auto& overlap=volume.overlaps[k];const auto native=layout.nativeCell(overlap.pressure);
                // Independent physical rectangle measure, rather than using
                // the factory's retained overlap.cellFraction as its oracle.
                const auto weight=vofMaterialCapacity(overlap.lower,overlap.upper,rz)/phase.metrics().cellVolumes()[overlap.pressure];
                const auto stock=phase.vofPublication().plicCells()[native].inventory;
                for(int component=0;component<2;++component){const auto reference=weight*stock[component];
                    const auto actual=component?represented[k].gas:represented[k].liquid;
                    if(reference==0)check(actual==0,"P0 projection invented a zero-stock phase");
                    else{check(actual>0,"P0 projection erased a positive independent phase stock");near(actual,reference,reference,"independent inventory overlap projection",1e-12L);}
                    expected[component]+=reference;}
            }
            const auto representedTotal=mass.representedPhaseVolumes()[u];
            for(int component=0;component<2;++component){const auto actual=component?representedTotal.gas:representedTotal.liquid;
                if(expected[component]==0)check(actual==0,"P0 zero represented phase changed");
                else near(actual,expected[component],expected[component],"independent whole-dual inventory projection",1e-12L);}
            const auto expectedMass=100*expected[0]+expected[1];near(mass.masses()[u],expectedMass,expectedMass,"actual P0 native momentum mass",1e-12L);
        }
    }
    check(material.report().usesPrimaryAlphaInventory,"LocalPlic primal material did not declare its independent inventory source");
    std::vector<long double> primalChange(layout.pressureCount()),primalOut(layout.pressureCount()),primalGross(layout.pressureCount());
    for(int u=0;u<layout.velocityCount();++u){const auto face=layout.domain().face(mesh,layout.nativeFace(u));const auto flux=material.faceIntegrals()[u];
        for(int native:{face.negativeCell,face.positiveCell})if(native>=0){const int pressure=layout.pressureDof(native);
            primalOut[pressure]+=(native==face.negativeCell?1.L:-1.L)*flux;primalGross[pressure]+=std::abs(flux);}}
    for(int pressure=0;pressure<layout.pressureCount();++pressure){const auto native=layout.nativeCell(pressure);
        const auto a=before.vofPublication().plicCells()[native].inventory,b=after.vofPublication().plicCells()[native].inventory;
        primalChange[pressure]=100*(b[0]-a[0])+(b[1]-a[1]);
        const auto scale=100*std::abs(b[0]-a[0])+std::abs(b[1]-a[1]);
        if(scale==0)check(material.cells()[pressure].massChange==0,"unchanged primal inventory acquired mass change");
        else near(material.cells()[pressure].massChange,primalChange[pressure],scale,"independent primal inventory mass change",1e-12L);
    }
    const auto dual=HydroDualMaterial::prepare(mesh,material,projected[0],projected[1],projected[2]);
    ConservativeVofMomentumContract::qualify(mesh,dual);
    long double peakTransfer=0,peakExchange=0,peakMappedR=0;
    for(int u=0;u<layout.velocityCount();++u){long double mappedChange=0,mappedR=0,mappedOut=0,mappedGross=0,changeGross=0,gross=0;
        for(const auto& overlap:geometry.volumes()[u].overlaps){const auto pressure=overlap.pressure;
            const auto weight=vofMaterialCapacity(overlap.lower,overlap.upper,rz)/before.metrics().cellVolumes()[pressure];
            mappedChange+=weight*primalChange[pressure];mappedOut+=weight*primalOut[pressure];
            mappedR+=weight*(primalChange[pressure]+primalOut[pressure]);mappedGross+=weight*primalGross[pressure];changeGross+=weight*std::abs(primalChange[pressure]);}
        for(int side=0;side<4;++side)for(int id:geometry.volumes()[u].sides[side])gross+=std::abs(dual.faceIntegrals()[id]);
        const auto& row=dual.cells()[u];const auto massScale=changeGross+std::abs(mappedChange),fluxScale=gross+mappedGross,totalScale=massScale+fluxScale;
        const auto zero=[&](long double actual,long double scale,const char* why,long double& peak){
            if(scale==0)check(actual==0,std::string(why)+" with zero source");
            else{near(actual,0,scale,why,512*std::numeric_limits<double>::epsilon());peak=std::max(peak,std::abs(actual)/scale);}};
        zero(row.geometryTransferChange,massScale,"P0 same-inventory geometry transfer commutator",peakTransfer);
        zero(row.outwardIntegral-mappedOut,fluxScale,"Coons same-F spatial exchange commutator",peakExchange);
        zero(row.integratedDefect-mappedR,totalScale,"actual P0 Coons dual R equals independently mapped primal R",peakMappedR);
        if(massScale>0)near(row.massChange,mappedChange,massScale,"canonical dual deltaM from independent stock differences",1e-12L);
    }
    check(dual.material().vofTransport().sameInstance(transport)&& &dual.faceRates()==&material.dualFaceRates(),"P0 paired dual detached its unique accepted material F owner");
    // Same existing XY/RZ P0 owners and mesh: each public entry must check
    // live coordinates on EVERY call, even after a successful validation.
    before.validateVofFor(mesh,before.vofIdentity(),before.timeLayer());
    projected[0].validateVofFor(mesh,before.vofIdentity(),before.timeLayer());
    dual.validateFor(mesh,dual.sources(),dual.timeLayer());
    auto& vertex=mesh.verts(mesh.vertIndexToLabel[0]);const auto saved=vertex.coord();
    vertex.setCoord(saved+Vector2D{1./1024,0});
    try{
        rejects([&]{before.validateVofFor(mesh,before.vofIdentity(),before.timeLayer());},"VOF phase live mesh mutation after prior valid call");
        rejects([&]{projected[0].validateVofFor(mesh,before.vofIdentity(),before.timeLayer());},"P0 VOF mass live mesh mutation after prior valid call");
        rejects([&]{dual.validateFor(mesh,dual.sources(),dual.timeLayer());},"VOF dual live mesh mutation after prior valid call");
    }catch(...){vertex.setCoord(saved);throw;}
    vertex.setCoord(saved);
    before.validateVofFor(mesh,before.vofIdentity(),before.timeLayer());
    projected[0].validateVofFor(mesh,before.vofIdentity(),before.timeLayer());
    dual.validateFor(mesh,dual.sources(),dual.timeLayer());
    std::cout<<"LOCAL_PLIC_LIVE_VALIDATION,"<<rz<<",phase_mass_dual_reject_and_restore,PASS\n";
    std::cout<<"LOCAL_PLIC_P0_COONS_PAIR,"<<rz<<','<<dual.cells().size()<<','<<peakTransfer<<','<<peakExchange<<','<<peakMappedR<<",PASS\n";
}

void l1MaterialConsumer(QuadtreeMesh& mesh,const Publication& source,const Slab& slab,const Transport& t)
{
    const bool rz=source.domain().axisymmetric();
    {
        // Reuse one real AMR leaf in each existing XY/RZ scene. The canonical
        // interior rectangle must retain its supplied vertices exactly; its
        // half-plane moment has an independent axial physical-measure oracle.
        namespace Cut=LAMP::Geometry::PlanarCellCut;
        const auto& cell=source.plicCells().front().rectangle;
        const Point origin{(static_cast<long double>(cell.lower.x)+cell.upper.x)/2,
            (static_cast<long double>(cell.lower.y)+cell.upper.y)/2};
        const auto hx=(static_cast<long double>(cell.upper.x)-cell.lower.x)/4;
        const auto hy=(static_cast<long double>(cell.upper.y)-cell.lower.y)/4;
        const std::vector<Point> rectangle{{-hx,-hy},{hx,-hy},{hx,hy},{-hx,hy}};
        const auto clipped=Cut::clipPolygonInFrame(cell,rectangle,origin);
        check(clipped.size()==rectangle.size(),"contained rectangle changed vertex count");
        for(std::size_t k=0;k<rectangle.size();++k)check(clipped[k].x==rectangle[k].x&&clipped[k].y==rectangle[k].y,
            "contained canonical rectangle changed an input coordinate");
        const auto value=Cut::overlapInFrame(cell,rectangle,origin,{1,0},0,rz);
        const auto total=4*hx*hy*(rz?twoPi*origin.x:1.L);
        const auto liquid=2*hx*hy*(rz?twoPi*(origin.x-hx/2):1.L);
        near(value.total,total,total,"contained rectangle actual XY/RZ capacity",1e-12L);
        near(value.liquid,liquid,total,"contained rectangle independent axial phase measure",1e-12L);
        // Cyclic clockwise input retains the general orientation path.
        auto clockwise=rectangle;std::reverse(clockwise.begin(),clockwise.end());
        const auto reversed=Cut::overlapInFrame(cell,clockwise,origin,{1,0},0,rz);
        near(reversed.total,total,total,"clockwise general rectangle capacity",1e-12L);
        near(reversed.liquid,liquid,total,"clockwise general rectangle liquid",1e-12L);
        auto invalid=rectangle;std::swap(invalid[1],invalid[2]);
        rejects([&]{Cut::clipPolygonInFrame(cell,invalid,origin);},"crossed rectangle cannot use contained fast path");
        invalid=rectangle;invalid[2].x=std::numeric_limits<long double>::infinity();
        rejects([&]{Cut::clipPolygonInFrame(cell,invalid,origin);},"nonfinite rectangle cannot use contained fast path");
        const std::vector<Point> empty{{-hx,-hy},{-hx,-hy},{-hx,hy},{-hx,hy}};
        check(Cut::clipPolygonInFrame(cell,empty,origin).empty(),"degenerate rectangle must retain zero-measure path");
    }

    // Same existing 2:1 scene: cached whole/native reversal, physical-measure
    // subface addback, and independently integrated opposite-side RT0 traces.
    for(int row=0;row<source.domain().nativeFaceCount();++row){const auto e=faceEndpoints(mesh,source.domain().face(mesh,row));
        const Point mid{e[0].x/2+e[1].x/2,e[0].y/2+e[1].y/2};
        const auto result=t.restrictNativeFaces(mesh,{{row,e[0],e[1]},{row,e[1],e[0]},{row,e[0],mid},{row,mid,e[1]},{row,e[0],e[0]}});
        for(int extent=0;extent<2;++extent){const auto& parent=extent?t.faces()[row].full:t.faces()[row].half;
            const auto fields=v4b1Fields(parent),whole=v4b1Fields(extent?result[0].full:result[0].half),reverse=v4b1Fields(extent?result[1].full:result[1].half);
            const auto a=v4b1Fields(extent?result[2].full:result[2].half),c=v4b1Fields(extent?result[3].full:result[3].half),empty=v4b1Fields(extent?result[4].full:result[4].half);
            const auto ratio=v4b1WorldArea(e[0],e[1],rz)>0?v4b1WorldArea(e[0],mid,rz)/v4b1WorldArea(e[0],e[1],rz):.5L;
            for(std::size_t k=0;k<fields.size();++k){check(fields[k]==whole[k]&&fields[k]==reverse[k],"local whole/reverse accepted-cache identity");if(k+1==fields.size())continue;
                check(empty[k]==0,"local empty accepted trace");const auto scale=std::abs(fields[k]);if(scale==0)check(a[k]==0&&c[k]==0,"local theoreticalzero P0 trace");
                else {near(a[k],ratio*fields[k],scale,"independent physical-measure P0 restriction",1e-12L);near(a[k]+c[k],fields[k],scale,"local accepted phase addback",1e-12L);}}
        }
    }
    for(int cell=0;cell<mesh.cellIndexNum();++cell){const auto box=bounds(mesh,cell);const Point center{box.x0/2+box.x1/2,box.y0/2+box.y1/2};
        for(int axis=0;axis<2;++axis){const bool horizontal=axis==1;const Point a=horizontal?Point{box.x0,center.y}:Point{center.x,box.y0},b=horizontal?Point{box.x1,center.y}:Point{center.x,box.y1};
            const auto result=t.measureInteriorFaces(mesh,{{cell,a,b},{cell,b,a}});
            const auto lower=axis?box.y0:box.x0,coordinate=axis?center.y:center.x;
            const auto tangentLo=axis?box.x0:box.y0,tangentHi=axis?box.x1:box.y1;
            const auto lowFraction=axis?.5L:capacity({box.x0,center.x,box.y0,box.y1},rz)/capacity(box,rz);
            for(int extent=0;extent<2;++extent){
                const auto actual=v4b1Fields(extent?result[0].full:result[0].half),reverse=v4b1Fields(extent?result[1].full:result[1].half);
                const auto low=l1BoundaryPart(mesh,source,t,cell,axis,0,tangentLo,tangentHi,extent),high=l1BoundaryPart(mesh,source,t,cell,axis,1,tangentLo,tangentHi,extent);
                const auto otherLo=axis?box.y0:box.x0,otherHi=axis?box.y1:box.x1;
                const auto negative=l1BoundaryPart(mesh,source,t,cell,1-axis,0,otherLo,otherHi,extent),positive=l1BoundaryPart(mesh,source,t,cell,1-axis,1,otherLo,otherHi,extent);
                const auto prefixNegative=l1BoundaryPart(mesh,source,t,cell,1-axis,0,lower,coordinate,extent),prefixPositive=l1BoundaryPart(mesh,source,t,cell,1-axis,1,lower,coordinate,extent);
                for(std::size_t k=0;k<actual.size();++k)check(actual[k]==reverse[k],"local interior canonical reversal");
                for(std::size_t k:{0U,3U,4U,5U,6U}){
                    const auto expected=low[k]+lowFraction*(high[k]-low[k]+positive[k]-negative[k])-(prefixPositive[k]-prefixNegative[k]);
                    const auto gross=std::abs(low[k])+std::abs(high[k])+std::abs(positive[k])+std::abs(negative[k]);
                    if(gross==0)check(actual[k]==0,"zero boundary source acquired an interior signed current");
                    else near(actual[k],expected,gross,"independent low-half-cell Gauss middle flux",1e-12L);
                }
                l1CoonsPayload(extent?result[0].full:result[0].half);
            }
        }
    }
    l1CoonsQuadrants(mesh,source,t);
    const auto geometry=HydroMomentumGeometry::build(mesh,vofMaterialMetrics(mesh,source));
    const auto before=HydroPhaseGeometry::build(mesh,geometry.metrics(),source),middle=HydroPhaseGeometry::build(mesh,geometry.metrics(),t.midpoint()),after=HydroPhaseGeometry::build(mesh,geometry.metrics(),t.endpoint());
#ifdef _OPENMP
    {
        // The existing real XY/RZ AMR dual overlaps provide every query;
        // compare concurrency choices without adding a layout or oracle.
        struct ThreadSettings {
            int threads=omp_get_max_threads(),dynamic=omp_get_dynamic();
            ~ThreadSettings(){omp_set_num_threads(threads);omp_set_dynamic(dynamic);}
        } restoreThreads;
        omp_set_dynamic(0);
        std::vector<HydroVofBox> queries;
        for(const auto& volume:geometry.volumes())for(const auto& overlap:volume.overlaps)
            queries.push_back({geometry.layout().nativeCell(overlap.pressure),{overlap.lower,overlap.upper}});
        check(queries.size()>=32,"existing AMR overlap scene lacks a parallel-sized box batch");
        omp_set_num_threads(1);const auto serial=before.measureVofBoxes(mesh,queries);
        omp_set_num_threads(8);int team=0;
#pragma omp parallel
        {
#pragma omp single
            team=omp_get_num_threads();
        }
        check(team>=2&&team<=8,"LocalPlic box comparison requires multiple actual workers within the requested budget");
        const auto parallel=before.measureVofBoxes(mesh,queries);
        const auto same=[&](const std::vector<HydroPhaseMeasure>& actual,const char* why){
            check(actual.size()==serial.size(),std::string(why)+" size changed");
            for(std::size_t i=0;i<serial.size();++i)check(actual[i].liquid==serial[i].liquid&&actual[i].gas==serial[i].gas,
                std::string(why)+" changed an actual phase measurement");
        };
        same(parallel,"requested-eight LocalPlic box batch");
        std::vector<HydroPhaseMeasure> nested;std::exception_ptr nestedFailure;
#pragma omp parallel num_threads(2)
        {
#pragma omp single
            {
                try{nested=before.measureVofBoxes(mesh,queries);}catch(...){nestedFailure=std::current_exception();}
            }
        }
        if(nestedFailure)std::rethrow_exception(nestedFailure);
        same(nested,"nested LocalPlic serial path");
        const std::vector<HydroVofBox> small(queries.begin(),queries.begin()+4);
        const auto smallResult=before.measureVofBoxes(mesh,small);
        for(std::size_t i=0;i<small.size();++i)check(smallResult[i].liquid==serial[i].liquid&&smallResult[i].gas==serial[i].gas,
            "small LocalPlic serial path changed an actual phase measurement");
        auto invalid=queries;
        const auto height=invalid[3].rectangle.upper.y-invalid[3].rectangle.lower.y;
        invalid[3].rectangle.lower.y-=2*height; // Valid capacity, then outside-owner bounds failure.
        invalid[invalid.size()-2].nativeCell=-1; // A later, distinct error must not win.
        const auto failure=[&](int threads){
            omp_set_num_threads(threads);std::string reason;
            try{(void)before.measureVofBoxes(mesh,invalid);}catch(const std::exception& error){reason=error.what();}
            return reason;
        };
        const auto firstSerial=failure(1),firstParallel=failure(8);
        check(firstSerial=="[HydroPhaseGeometry] VOF query box leaves its original native leaf"&&firstParallel==firstSerial,
            "parallel LocalPlic batch changed the first failing query or its diagnostic");
        check(before.vofIdentity()==source.identity()&&before.vofPublication().plicCells().data()==source.plicCells().data(),
            "LocalPlic box queries changed their immutable source owner");
        std::cout<<"LOCAL_PLIC_PARALLEL_BOXES,"<<rz<<','<<queries.size()<<",requested_threads=8,actual_threads="<<team<<",PASS\n";
    }
#endif
    const auto material=HydroMaterialTransport::prepareGeometricVof(mesh,geometry,before,middle,after,t);
    check(material.vofTransport().sameInstance(t)&&material.sources().vofStamp&&material.sources().vofStamp->sameInstance(t),"LocalPlic material lost its accepted flux owner");
    l1InventoryProjectionPair(mesh,source,t,geometry,material,before,middle,after);
    const auto mass0=HydroMomentumMass::build(mesh,geometry,before,100,1),massM=HydroMomentumMass::build(mesh,geometry,middle,100,1),mass1=HydroMomentumMass::build(mesh,geometry,after,100,1);
    const auto dual=HydroDualMaterial::prepare(mesh,material,mass0,massM,mass1);const auto& layout=geometry.layout();
    check(!dual.usesConservativeVofMomentum(),"LocalPlic Sharp reference must retain original nonconservative dispatch");
    {
        // Same native AMR topology, clocks, original solenoidal velocity and
        // sharp physical M. Only the typed publication is homogeneous gas.
        const auto gasSource=Publication::prepareLocalPlic(mesh,source.domain(),source.timeLayer(),std::vector<double>(source.alpha().size(),0));
        const auto gasSlab=Slab::prepare(mesh,gasSource,slab.endTime(),slab.beginNormal(),slab.endNormal());
        const auto gasTransport=Transport::prepareLocalPlicSweep(mesh,gasSource,gasSlab,100,1);
        const auto gasBegin=HydroPhaseGeometry::build(mesh,geometry.metrics(),gasSource),gasMiddle=HydroPhaseGeometry::build(mesh,geometry.metrics(),gasTransport.midpoint()),gasEnd=HydroPhaseGeometry::build(mesh,geometry.metrics(),gasTransport.endpoint());
        const auto gasMaterial=HydroMaterialTransport::prepareGeometricVof(mesh,geometry,gasBegin,gasMiddle,gasEnd,gasTransport);
        const auto gasDual=HydroDualMaterial::prepare(mesh,gasMaterial,
            HydroMomentumMass::build(mesh,geometry,gasBegin,100,1),HydroMomentumMass::build(mesh,geometry,gasMiddle,100,1),HydroMomentumMass::build(mesh,geometry,gasEnd,100,1));
        long double peakR=0,peakExchange=0;
        for(std::size_t u=0;u<gasDual.cells().size();++u){const auto& row=gasDual.cells()[u];long double gross=0;
            for(int side=0;side<4;++side)for(int id:geometry.volumes()[u].sides[side])gross+=std::abs(gasDual.faceIntegrals()[id]);
            check(row.massChange==0,"homogeneous gas sharp mass changed without physical motion of density");
            if(gross==0)check(row.integratedDefect==0,"homogeneous zero-flow gas dual acquired R");
            else{near(row.integratedDefect,0,gross,"actual homogeneous-gas dual R at roundoff",512*std::numeric_limits<double>::epsilon());peakR=std::max(peakR,std::abs(row.integratedDefect)/gross);}
            long double mapped=0,mappedGross=0;
            for(const auto& overlap:geometry.volumes()[u].overlaps){mapped+=overlap.cellFraction*gasMaterial.cells()[overlap.pressure].outwardIntegral;
                mappedGross+=overlap.cellFraction*std::abs(gasMaterial.cells()[overlap.pressure].outwardIntegral);}
            const auto exchange=row.outwardIntegral-mapped,scale=gross+mappedGross;
            if(scale==0)check(exchange==0,"zero gas source acquired spatial exchange");
            else{near(exchange,0,scale,"actual homogeneous-gas partial-dual spatial exchange",512*std::numeric_limits<double>::epsilon());peakExchange=std::max(peakExchange,std::abs(exchange)/scale);}
        }
        const auto interior=gasTransport.measureInteriorFaces(mesh,{{0,{(bounds(mesh,0).x0+bounds(mesh,0).x1)/2,bounds(mesh,0).y0},
            {(bounds(mesh,0).x0+bounds(mesh,0).x1)/2,bounds(mesh,0).y1}}});
        for(const auto* flux:{&interior[0].half,&interior[0].full})check(flux->liquid==0&&flux->positiveLiquid==0&&flux->negativeLiquid==0,
            "Coons gas source invented a liquid phase current");
        l1CoonsQuadrants(mesh,gasSource,gasTransport);
        check(gasSource.alpha()==std::vector<double>(source.alpha().size(),0)&&gasMaterial.vofTransport().sameInstance(gasTransport),"Coons pure-gas measurements changed accepted source or F owner");
        std::cout<<"LOCAL_PLIC_COONS_PURE_GAS_DUAL,"<<rz<<','<<gasDual.cells().size()<<','<<peakR<<','<<peakExchange<<",PASS\n";
    }
    std::vector<double> normal(layout.velocityCount());for(int u=0;u<layout.velocityCount();++u){normal[u]=slab.beginNormal()[layout.nativeFace(u)];const auto& parent=t.faces()[layout.nativeFace(u)].full;
        near(material.faceIntegrals()[u],parent.mass,std::max(1.L,std::abs(parent.mass)),"actual LocalPlic primal material F",1e-12L);}
    const HydroNativeValues values(layout,normal,std::vector<double>(layout.pressureCount()));const auto impulse=dual.consistencyImpulse(mesh,values,dual.sources(),dual.timeLayer());
    for(std::size_t u=0;u<impulse.size();++u)check(std::isfinite(dual.cells()[u].defectRate)&&impulse[u]==normal[u]*dual.cells()[u].integratedDefect,"genuine LocalPlic dual R/uR identity");
    {
        // Existing XY/RZ 2:1 streamfunction field: independent WHOLE-cell
        // endpoint differences, not the implementation's native-side loop.
        const auto nativeState=HydroTwoPhaseState::createVof(mesh,before,values);
        const auto sampled=nativeState.sampleCellCenterVelocities(mesh);
        long double componentScale=0;for(double value:normal)componentScale=std::max(componentScale,std::abs(static_cast<long double>(value)));
        for(int p=0;p<layout.pressureCount();++p){const auto box=bounds(mesh,layout.nativeCell(p));
            const long double x=(box.x0+box.x1)/2;
            const auto psi=[](long double x,long double y){return vofMaterialPsi({static_cast<double>(x),static_cast<double>(y)},0);};
            const auto lower=psi(box.x0,box.y1)-psi(box.x0,box.y0),upper=psi(box.x1,box.y1)-psi(box.x1,box.y0);
            const auto xi=rz?(x*x-box.x0*box.x0)/(box.x1*box.x1-box.x0*box.x0):.5L;
            const auto first=((1-xi)*lower+xi*upper)/((box.y1-box.y0)*(rz?twoPi*x:1.L));
            const auto second=(-psi(box.x1,box.y0)+psi(box.x0,box.y0)-psi(box.x1,box.y1)+psi(box.x0,box.y1))/2/((box.x1-box.x0)*(rz?twoPi*x:1.L));
            near(sampled[p].x,first,componentScale,"independent native RT0/Piola center first component",1e-12L);
            near(sampled[p].y,second,componentScale,"independent native RT0/Piola center second component",1e-12L);
            if(!rz)check(std::abs(sampled[p].x)<=componentScale&&std::abs(sampled[p].y)<=componentScale,"XY RT0 component exceeded global native component bound");
        }
        check(nativeState.values().normalVelocity()==normal,"VOF sampling modified authoritative native U");
    }
    {
        // One existing XY/RZ AMR trace, independently integrating its affine
        // parameter interval and physical radius; no new geometry fixture.
        unsigned checked=0;
        for(int u=0;u<layout.velocityCount()&&!checked;++u){const auto face=layout.domain().face(mesh,layout.nativeFace(u));if(face.area==0)continue;
            const auto ends=faceEndpoints(mesh,face);
            for(int native:{face.negativeCell,face.positiveCell})if(native>=0&&!checked){const auto& plane=source.plicCells()[native];
                const auto& rectangle=plane.rectangle;const long double cx=(static_cast<long double>(rectangle.lower.x)+rectangle.upper.x)/2,cy=(static_cast<long double>(rectangle.lower.y)+rectangle.upper.y)/2;
                const auto norm=std::hypot(static_cast<long double>(plane.normal.x),static_cast<long double>(plane.normal.y));
                const auto value=[&](Point p){return plane.normal.x/norm*(p.x-cx)+plane.normal.y/norm*(p.y-cy)-plane.offset;};
                const auto a=value(ends[0]),b=value(ends[1]);if(!((a<0&&b>0)||(a>0&&b<0)))continue;
                const auto root=std::abs(a)/(std::abs(a)+std::abs(b)),tail=std::abs(b)/(std::abs(a)+std::abs(b));
                const auto dx=ends[1].x-ends[0].x,dy=ends[1].y-ends[0].y,L=std::hypot(dx,dy);
                const auto first=L*root*(rz?twoPi*(ends[0].x+root*dx/2):1.L),last=L*tail*(rz?twoPi*(ends[1].x-tail*dx/2):1.L);
                const auto trace=face.negativeCell==native?before.faceTraces()[u].minus:before.faceTraces()[u].plus;
                const auto liquid=a<0?first:last,gas=a<0?last:first;
                near(trace.liquid,liquid,liquid,"independent XY/RZ PLIC parameter liquid trace",1e-12L);
                near(trace.gas,gas,gas,"independent XY/RZ PLIC parameter gas trace",1e-12L);++checked;
            }
        }
        check(checked==1,"existing XY/RZ AMR scene lacks a split physical trace");
    }
    if(!rz){
        // Reuse this real 2:1 mesh and its actual partial-box decomposition.
        // An oblique positive micro phase exercises whole-box versus polygon
        // arithmetic; zero phases and both positive stocks must remain intact.
        std::vector<std::array<long double,2>> microStocks(source.plicCells().size());int microCell=-1;
        for(std::size_t ci=0;ci<microStocks.size();++ci){const auto box=bounds(mesh,static_cast<int>(ci));
            const auto cx=(box.x0+box.x1)/2,cy=(box.y0+box.y1)/2,c=source.cells()[ci].capacity;
            // Positive affine background selects an upper-right oblique cut.
            // The selected fine leaf has a symmetric fine native one-ring.
            const auto f=std::clamp(.1L+.4L*(cx+.63L*cy),0.L,1.L);microStocks[ci]={f*c,(1-f)*c};
            if(cx<.25L&&cy<.25L)microStocks[ci]={0,c};
            if(cx<.25L&&cy>.75L)microStocks[ci]={c,0};
            if(box.x0==.75L&&box.y0==.25L)microCell=static_cast<int>(ci);}
        check(microCell>=0,"existing AMR scene lacks selected interior micro cell");
        // Retain a small real shape well ABOVE the explicit geometry resolution.
        const auto microCapacity=source.cells()[microCell].capacity,positiveGeometryStock=1e-10L*microCapacity;
        microStocks[microCell]={positiveGeometryStock,microCapacity-positiveGeometryStock};
        const auto retainedStocks=microStocks;
        const auto micro=Publication::prepareLocalPlicInventory(mesh,source.domain(),source.timeLayer(),microStocks,{});
        const auto microPhase=HydroPhaseGeometry::build(mesh,geometry.metrics(),micro);
        const auto microMass=HydroMomentumMass::build(mesh,geometry,microPhase,100,1);
        const auto pressure=layout.pressureDof(microCell);const auto actual=microPhase.cellVolumes()[pressure].liquid;
        check(positiveGeometryStock>64.L*std::numeric_limits<double>::epsilon()*microCapacity
            &&actual>0&&actual<1e-8L*microCapacity,"above-resolution actual PLIC geometry lost its positive small phase");
        std::array<long double,2> partition{};unsigned strictZeros=0;
        for(int u=0;u<layout.velocityCount();++u){const int component=layout.normalAxis(u)==HydroDofLayout::NormalAxis::First?0:1;
            const auto phases=microMass.phaseVolumes()[u];
            check(microMass.masses()[u]==static_cast<double>(100*phases.liquid+phases.gas),"micro material mass remapped measured phase geometry");
            for(std::size_t k=0;k<geometry.volumes()[u].overlaps.size();++k){const auto p=geometry.volumes()[u].overlaps[k].pressure;
                const auto q=microMass.overlapPhaseVolumes()[u][k];const auto ref=microPhase.cellVolumes()[p];
                check(std::isfinite(q.liquid)&&std::isfinite(q.gas)&&q.liquid>=0&&q.gas>=0,"invalid micro partial-box geometry");
                if(ref.liquid==0){check(q.liquid==0,"pure gas acquired a false micro liquid");++strictZeros;}
                if(ref.gas==0){check(q.gas==0,"pure liquid acquired a false micro gas");++strictZeros;}
                if(p==pressure)partition[component]+=q.liquid;}}
        long double witnessRelative=0;
        for(const auto measured:partition){check(measured>0,"positive micro phase vanished from a velocity partition");
            const auto error=std::abs(measured-actual);witnessRelative=std::max(witnessRelative,error/actual);
            check(error<=1e-12L*actual+64.L*std::numeric_limits<double>::epsilon()*source.cells()[microCell].capacity,"real micro partition exceeds declared allowance");}
        // Native binary half-boxes can take identical floating paths and give
        // exactly zero difference. Separately use the SAME actual PLIC owner
        // with one fixed non-dyadic 3/7 split in each axis: this named oblique
        // geometry witness does not require every native partition to differ.
        const auto& rectangle=micro.plicCells()[microCell].rectangle;
        const auto& plane=micro.plicCells()[microCell];
        check(plane.normal.x<0&&plane.normal.y<0,"named upper-right oblique witness changed quadrant");
        long double obliquePartialRelative=0;
        for(int axis=0;axis<2;++axis){auto lower=rectangle,upper=rectangle;
            if(axis==0){const double split=rectangle.lower.x+(rectangle.upper.x-rectangle.lower.x)*3/7;lower.upper.x=split;upper.lower.x=split;}
            else {const double split=rectangle.lower.y+(rectangle.upper.y-rectangle.lower.y)*3/7;lower.upper.y=split;upper.lower.y=split;}
            const auto partial=microPhase.measureVofBoxes(mesh,{{microCell,lower},{microCell,upper}});
            const auto measured=partial[0].liquid+partial[1].liquid;
            check(measured>0,"named oblique actual box partition lost positive micro liquid");
            const auto error=std::abs(measured-actual);obliquePartialRelative=std::max(obliquePartialRelative,error/actual);
            check(error<=1e-12L*actual+64.L*std::numeric_limits<double>::epsilon()*source.cells()[microCell].capacity,"named oblique actual box partition exceeds declared allowance");}
        check(strictZeros>0,"micro scene lacks strict pure-phase zero checks");
        // Exact agreement is valid; no artificial minimum roundoff defect.
        // The frozen physical run accepted the formerly rejected step 83.
        check(microMass.report().maxRelativeCellPartition>=witnessRelative
            &&std::isfinite(microMass.report().maxLocalPlicCapacityScaledCellPartition)
            &&microMass.report().maxLocalPlicPartitionBudgetFraction<=1,"micro diagnostics hid the actual partition error");
        for(std::size_t ci=0;ci<microStocks.size();++ci)check(micro.plicCells()[ci].inventory==retainedStocks[ci],"micro mass preparation changed retained independent phase stocks");
        std::cout<<"LOCAL_PLIC_MICRO_PARTITION,"<<actual<<','<<partition[0]<<','<<partition[1]<<','<<witnessRelative<<','
            <<obliquePartialRelative<<','<<microMass.report().maxLocalPlicCapacityScaledCellPartition<<','<<microMass.report().maxLocalPlicPartitionBudgetFraction<<",PASS\n";
        // Same existing cell, now BELOW the declared LocalPlic geometry
        // resolution: keep both independent stocks but derive strict pure gas.
        // This replaces the former below-resolution world-Point shape premise;
        // ordinary XY/RZ physical parameter traces and Core's own gates stay.
        auto traceStocks=retainedStocks;const auto eps=std::numeric_limits<long double>::epsilon();
        traceStocks[microCell]={source.cells()[microCell].capacity*eps*eps/2,source.cells()[microCell].capacity};
        const auto traceSource=Publication::prepareLocalPlicInventory(mesh,source.domain(),source.timeLayer(),traceStocks,{});
        const auto tracePhase=HydroPhaseGeometry::build(mesh,geometry.metrics(),traceSource);
        const auto& tiny=traceSource.plicCells()[microCell];const auto tinyCapacity=source.cells()[microCell].capacity;
        check(tiny.inventory[0]>0&&tiny.inventory[1]>0&&traceSource.alpha()[microCell]>0&&!tiny.resolved,
            "subresolution source discarded positive stock or claimed resolved geometry");
        check(tiny.measured[0]==0&&tiny.measured[1]==tinyCapacity,
            "subresolution liquid-minority source must derive exact pure-gas actual plane");
        const auto tinyPressure=layout.pressureDof(microCell);
        check(tracePhase.liquidFractions()[tinyPressure]>0&&tracePhase.cellVolumes()[tinyPressure].liquid==0
            &&tracePhase.cellVolumes()[tinyPressure].gas==tinyCapacity,"subresolution phase owner confused stock alpha with actual pure geometry");
        const auto traceProjection=HydroMomentumMass::buildCellAverageProjection(mesh,geometry,tracePhase,100,1);
        std::array<long double,2> retainedMicro{};unsigned positiveOverlaps=0;
        for(int u=0;u<layout.velocityCount();++u){const int component=layout.normalAxis(u)==HydroDofLayout::NormalAxis::First?0:1;
            for(std::size_t k=0;k<geometry.volumes()[u].overlaps.size();++k){const auto& overlap=geometry.volumes()[u].overlaps[k];if(overlap.pressure!=tinyPressure)continue;
                const auto represented=traceProjection.representedOverlapPhaseVolumes()[u][k];
                check(represented.liquid>0&&represented.gas>0,"P0 independent subresolution inventory phase was erased");
                retainedMicro[component]+=represented.liquid;++positiveOverlaps;}}
        for(auto actual:retainedMicro)near(actual,tiny.inventory[0],tiny.inventory[0],"P0 retains positive stock while actual sharp plane is pure gas",1e-12L);
        check(positiveOverlaps>0&&traceProjection.vofIdentity()==traceSource.identity(),"P0 micro inventory lost its actual owner");
        unsigned pureTraces=0;
        for(int u=0;u<layout.velocityCount();++u){const auto face=layout.domain().face(mesh,layout.nativeFace(u));if(face.negativeCell!=microCell&&face.positiveCell!=microCell)continue;
            const auto measured=face.negativeCell==microCell?tracePhase.faceTraces()[u].minus:tracePhase.faceTraces()[u].plus;
            check(measured.liquid==0&&measured.gas==face.area,"subresolution pure-plane trace manufactured a tiny interface");++pureTraces;}
        check(pureTraces==4,"existing interior subresolution cell lacks four real native traces");
        const auto traceSlab=Slab::prepare(mesh,traceSource,slab.endTime(),slab.beginNormal(),slab.endNormal());
        const auto traceTransport=Transport::prepareLocalPlicSweep(mesh,traceSource,traceSlab,100,1);
        check(traceTransport.sourceIdentity()==traceSource.identity()&&traceTransport.report().localUnresolvedDonors>0,
            "subresolution transport did not retain the real source and declared donor approximation");
        const bool firstHorizontal=(traceSlab.endTime().step%2)==0;unsigned approximatedFaces=0;
        for(const auto& accepted:traceTransport.faces()){const auto face=traceSource.domain().face(mesh,accepted.nativeRow);
            if(face.horizontal!=firstHorizontal)continue;
            for(const auto* flux:{&accepted.half,&accepted.full}){if(flux->total==0)continue;
                const int donor=flux->total>0?face.negativeCell:face.positiveCell;if(donor!=microCell)continue;
                const auto expected=flux->total*(tiny.inventory[0]/tinyCapacity);
                check(expected!=0&&flux->liquid!=0,"declared unresolved donor approximation erased a retained positive micro flux");
                near(flux->liquid,expected,std::abs(expected),"actual accepted subresolution donor-fraction shared liquid F",1e-12L);
                check(flux->mass==100*flux->liquid+flux->gas,"subresolution material flux lost its same accepted phase F");++approximatedFaces;}}
        check(approximatedFaces>0,"existing source velocity did not exercise subresolution donor outflow");
        std::array<long double,2> oldTotal{},newTotal{};long double totalCapacity=0;
        for(std::size_t ci=0;ci<traceStocks.size();++ci){check(traceSource.plicCells()[ci].inventory==traceStocks[ci],"subresolution trace/transport changed retained independent source stocks");
            const auto& candidate=traceTransport.endpoint().plicCells()[ci];totalCapacity+=source.cells()[ci].capacity;
            for(int phase=0;phase<2;++phase){oldTotal[phase]+=traceStocks[ci][phase];newTotal[phase]+=candidate.inventory[phase];}}
        for(int phase=0;phase<2;++phase)near(newTotal[phase],oldTotal[phase],totalCapacity,"subresolution closed-domain shared transport phase inventory",1e-12L);
        // BEGIN D11 EXISTING MICRO TRANSPORT MATERIAL CONSUMER
        const auto traceBegin=HydroPhaseGeometry::build(mesh,geometry.metrics(),traceSource);
        const auto traceMiddle=HydroPhaseGeometry::build(mesh,geometry.metrics(),traceTransport.midpoint());
        const auto traceEnd=HydroPhaseGeometry::build(mesh,geometry.metrics(),traceTransport.endpoint());
        const auto traceMaterial=HydroMaterialTransport::prepareGeometricVof(mesh,geometry,traceBegin,traceMiddle,traceEnd,traceTransport);
        check(traceMaterial.vofTransport().sameInstance(traceTransport),"existing micro transport lost its actual atom consumer");
        // END D11 EXISTING MICRO TRANSPORT MATERIAL CONSUMER
        std::cout<<"LOCAL_PLIC_SUBRESOLUTION_TRACE,"<<tiny.inventory[0]<<','<<tiny.measured[0]<<','<<pureTraces<<','<<approximatedFaces
            <<','<<traceTransport.report().localUnresolvedDonors<<",PASS\n";
        // Reuse this exact accepted LocalPlic/AMR source and native layout.
        // Zero flow isolates real State/Step plumbing without invoking HF.
        const HydroNativeValues zero(layout,std::vector<double>(layout.velocityCount()),std::vector<double>(layout.pressureCount()));
        auto accepted=HydroTwoPhaseState::createVof(mesh,before,zero);
        check(accepted.model()==HydroPhaseGeometry::Model::GeometricVof&&accepted.phaseGeometry().vofIdentity()==source.identity(),
            "pure VOF accepted state must retain the real LocalPlic owner");
        HydroVofAcConfig config;config.liquidDensity=100;config.gasDensity=1;config.liquidViscosity=config.gasViscosity=config.surfaceTension=0;
        check(config.momentumDensityRepresentation==HydroMomentumMass::DensityRepresentation::CellAverageProjection,"VOF Step default must be conservative P0");
        auto sharpConfig=config;sharpConfig.momentumDensityRepresentation=HydroMomentumMass::DensityRepresentation::SharpGeometry;
        rejects([&]{(void)HydroVofStep::prepare(mesh,accepted,slab.endTime(),sharpConfig);},"VOF production Step must reject explicit Sharp mass");
        const auto step=HydroVofStep::prepare(mesh,accepted,slab.endTime(),config);
        const auto& candidate=step.candidateState();const auto& frozen=step.frozenMidpoint();const auto& actualDual=frozen.material();
        const auto& solved=frozen.solveResult().values();
        check(step.report().linearSolves>=1&&step.report().linearSolves<=2&&!step.hasCapillaryLoad()
            &&candidate.phaseGeometry().vofPublication().alpha()==source.alpha()
            &&solved.normalVelocity()==candidate.values().normalVelocity()
            &&std::all_of(solved.normalVelocity().begin(),solved.normalVelocity().end(),[](double value){return value==0;}),
            "pure VOF zero candidate must consume the original frozen NS result with unchanged alpha and no force");
        const auto actualImpulse=actualDual.consistencyImpulse(mesh,step.solvedMidpointValues(),actualDual.sources(),actualDual.timeLayer());
        check(actualDual.sources().vofStamp&&actualDual.sources().vofStamp->sameInstance(step.transport())
            &&std::all_of(actualDual.material().faceIntegrals().begin(),actualDual.material().faceIntegrals().end(),[](long double value){return value==0;})
            &&std::all_of(actualDual.cells().begin(),actualDual.cells().end(),[](const auto& cell){return cell.integratedDefect==0&&cell.defectRate==0;})
            &&std::all_of(actualImpulse.begin(),actualImpulse.end(),[](long double value){return value==0;}),
            "pure VOF actual shared F, M/R and uR must retain the same zero-flow operator");
        auto projectedAccepted=HydroTwoPhaseState::createVof(mesh,before,zero);
        auto projectedConfig=config;projectedConfig.momentumDensityRepresentation=HydroMomentumMass::DensityRepresentation::CellAverageProjection;
        const auto projectedStep=HydroVofStep::prepare(mesh,projectedAccepted,slab.endTime(),projectedConfig);
        const auto& projectedFrozen=projectedStep.frozenMidpoint();const auto& projectedDual=projectedFrozen.material();
        for(const auto* mass:{&projectedDual.beginMass(),&projectedDual.midpointMass(),&projectedDual.endMass()})
            check(mass->representation()==HydroMomentumMass::DensityRepresentation::CellAverageProjection,"actual original NS did not consume explicit VOF P0 mass");
        check(projectedDual.material().vofTransport().sameInstance(projectedStep.transport())
            &&projectedStep.candidateState().retainedMass()!=nullptr
            &&projectedStep.candidateState().retainedMass()->masses().data()==projectedDual.endMass().masses().data(),
            "P0 frozen NS and accepted end candidate did not retain the same real mass owner");
        const auto& matrix=projectedFrozen.physicalVelocityMatrix();
        for(int u=0;u<layout.velocityCount();++u){const auto expected=projectedDual.endMass().masses()[u]/projectedDual.deltaTime();
            near(matrix.coeff(u,u),expected,expected,"original physical NS matrix consumes actual P0 endMass/dt",1e-12L);
            check(projectedDual.cells()[u].integratedDefect==0&&projectedStep.candidateState().values().normalVelocity()[u]==0,
                "P0 zero-flow original NS manufactured R or velocity");}
        for(int outer=0;outer<matrix.outerSize();++outer)for(HydroFrozenMidpoint::Matrix::InnerIterator entry(matrix,outer);entry;++entry)
            if(entry.row()!=entry.col())check(entry.value()==0,"zero-flow P0 physical NS matrix has unexpected off-diagonal load");
        std::cout<<"LOCAL_PLIC_P0_ORIGINAL_NS,"<<projectedDual.cells().size()<<",PASS\n";
        step.commit(mesh,accepted);
        check(accepted.sameIdentity(candidate)&&sameClock(accepted.timeLayer(),slab.endTime())
            &&accepted.phaseGeometry().vofPublication().alpha()==source.alpha(),"pure VOF immutable commit must adopt this exact end candidate");
        rejects([&]{step.commit(mesh,accepted);},"stale pure VOF commit must reject");
    }
    std::cout<<"LOCAL_PLIC_MATERIAL_CONSUMER,"<<rz<<','<<material.atomIntegrals().size()<<','<<dual.cells().size()<<",PASS\n";
}
// Read-only adaptive-rate prerequisites. The first field intentionally has
// nonzero cell divergence: rates are geometry measures, not a transport solve.
void localPlicStabilityRates()
{
    for(bool rz:{false,true}){
        stage=rz?"LocalPlic stability RZ measures":"LocalPlic stability XY 2:1 columns";
        auto mesh=vofMaterialMesh(4);
        const auto source=Publication::prepareLocalPlic(*mesh,domain(*mesh,rz),Clock{},l1PlaneAlpha(*mesh,rz,0));
        int coarse=-1;
        for(int cell=0;cell<mesh->cellIndexNum();++cell){const auto box=bounds(*mesh,cell);
            if(box.x0==.25L&&box.x1==.5L&&box.y0==.25L&&box.y1==.5L)coarse=cell;}
        check(coarse>=0,"rate witness must retain the original actual coarse leaf");
        std::vector<double> normal(source.domain().nativeFaceCount());unsigned radialSubfaces=0,axialFaces=0;
        for(int row=0;row<source.domain().nativeFaceCount();++row){const auto f=source.domain().face(*mesh,row);
            if(f.negativeCell!=coarse&&f.positiveCell!=coarse)continue;
            if(!f.horizontal&&f.center.x==.5){normal[row]=1;++radialSubfaces;}
            else if(f.horizontal&&f.center.y==.25){normal[row]=-1;++axialFaces;}
            else if(f.horizontal&&f.center.y==.5){normal[row]=1;++axialFaces;}
        }
        check(radialSubfaces==2&&axialFaces==2,"rate witness must exercise two actual hanging subfaces and both whole axial sides");
        const auto owner=source.identity();const auto* storage=source.plicCells().data();const auto oldAlpha=source.alpha();
        const auto rates=Transport::localPlicCourantRates(*mesh,source,normal);
        // XY: each donor column is .25*.125, fine receiver is .125^2;
        // outgoing coarse sums .25 + .25 + .25 over .25^2.
        // RZ: right-subface area=2*pi*.5*.125, column volume=
        // pi*(.5^2-.25^2)*.125; fine receiving volume=
        // pi*(.625^2-.5^2)*.125. Both axial sides each contribute 4/time.
        const long double column=rz?16.L/3:4.L,incoming=rz?64.L/9:8.L,outgoing=rz?40.L/3:12.L;
        near(rates.maximumLocalColumn,column,column,"independent actual donor subface-column rate",1e-12L);
        near(rates.maximumIncoming,incoming,incoming,"independent fine receiver rate",1e-12L);
        near(rates.maximumOutflow,outgoing,outgoing,"independent accumulated complete coarse-cell rate",1e-12L);
        for(auto& value:normal)value=-value;
        const auto reverse=Transport::localPlicCourantRates(*mesh,source,normal);
        near(reverse.maximumLocalColumn,incoming,incoming,"reversed actual fine donor-column rate",1e-12L);
        near(reverse.maximumIncoming,outgoing,outgoing,"reversed complete-cell incoming accumulation",1e-12L);
        near(reverse.maximumOutflow,incoming,incoming,"reversed complete-cell outflow",1e-12L);
        rejects([&]{Transport::localPlicCourantRates(*mesh,source,std::vector<double>(normal.size()-1));},"rate native dimensions");
        auto invalid=normal;invalid[0]=std::numeric_limits<double>::quiet_NaN();
        rejects([&]{Transport::localPlicCourantRates(*mesh,source,invalid);},"rate rejects a nonfinite actual native mean");
        // Use the existing closed discrete streamfunction for an actual
        // accepted Hydro state. Its compact velocity order is mapped back
        // into the original native rows without reconstructing point values.
        const auto layout=HydroDofLayout::build(*mesh,source.domain());
        const auto metrics=HydroPhysicalDivergence::build(*mesh,HydroIntegralDivergence::build(*mesh,layout));
        const auto native=vofMaterialNormals(*mesh,source,0);std::vector<double> compact(layout.velocityCount());
        for(int u=0;u<layout.velocityCount();++u)compact[u]=native[layout.nativeFace(u)];
        const auto accepted=HydroTwoPhaseState::createVof(*mesh,HydroPhaseGeometry::build(*mesh,metrics,source),
            HydroNativeValues(layout,compact,std::vector<double>(layout.pressureCount())));
        const auto oldTime=accepted.timeLayer();const auto acceptedVelocity=accepted.values().normalVelocity();
        HydroVofAcConfig config;config.liquidDensity=100;config.gasDensity=1;config.surfaceTension=0;
        const auto hydro=HydroVofStep::stabilityRates(*mesh,accepted,config);
        const auto direct=Transport::localPlicCourantRates(*mesh,source,native);
        check(hydro.transport.maximumLocalColumn==direct.maximumLocalColumn
            &&hydro.transport.maximumIncoming==direct.maximumIncoming&&hydro.transport.maximumOutflow==direct.maximumOutflow,
            "actual accepted Hydro means lost their canonical native-row mapping");
        check(hydro.flowRate==std::max({direct.maximumLocalColumn,direct.maximumIncoming,direct.maximumOutflow})&&hydro.flowRate>0,
            "accepted-state flow rate must retain all three actual rates");
        check(!hydro.capillaryAvailable&&hydro.capillaryRate==0&&hydro.interfaceMinimumSpacing==0,"sigma-zero declares no capillary constraint");
        if(!rz){
            config.surfaceTension=2;
            const auto capillary=HydroVofStep::stabilityRates(*mesh,accepted,config);
            // The x=.73 interface cuts the .125-spaced actual refined half.
            const long double h=.125L,rhoMean=50.5L;
            const auto expected=std::sqrt(std::acos(-1.L)*2/(rhoMean*h*h*h));
            check(capillary.capillaryAvailable&&capillary.interfaceMinimumSpacing==h,"capillary scale uses the actual mixed refined leaf spacing");
            near(capillary.capillaryRate,expected,expected,"independent interface-leaf capillary inverse time",1e-12L);
            check(capillary.flowRate==hydro.flowRate,"sigma must not change the retained actual flow rate");
            const auto pure=Publication::prepareLocalPlic(*mesh,source.domain(),Clock{},std::vector<double>(mesh->cellIndexNum()));
            const auto pureAccepted=HydroTwoPhaseState::createVof(*mesh,HydroPhaseGeometry::build(*mesh,metrics,pure),
                HydroNativeValues(layout,std::vector<double>(layout.velocityCount()),std::vector<double>(layout.pressureCount())));
            const auto pureRates=HydroVofStep::stabilityRates(*mesh,pureAccepted,config);
            check(pureRates.flowRate==0&&!pureRates.capillaryAvailable&&pureRates.capillaryRate==0&&pureRates.interfaceMinimumSpacing==0,
                "stationary pure phase declares neither flow nor capillary constraint");
        }
        check(source.identity()==owner&&source.plicCells().data()==storage&&source.alpha()==oldAlpha
            &&sameClock(accepted.timeLayer(),oldTime)&&accepted.values().normalVelocity()==acceptedVelocity,
            "read-only stability queries changed source identity, storage, clock or actual accepted values");
        std::cout<<"LOCAL_PLIC_STABILITY_RATES,"<<rz<<','<<rates.maximumLocalColumn<<','<<rates.maximumIncoming<<','<<rates.maximumOutflow<<",PASS\n";
    }
}
// BEGIN D11 ACTUAL NATIVE ATOM SUBNORMAL CONTRACT
// BEGIN D12 PURE VOF SMALL-VELOCITY RANGE PROBE
void localPlicConvectionRangeOnly(int exponent)
{
    check(exponent==120||exponent==160,"range probe exponent must be 120 or 160");
    const long double scale=exponent==120?1e-120L:1e-160L;
    const auto mark=[&](const char* operation){
        stage="D12 pure VOF range "+std::to_string(exponent)+" "+operation;
        std::cout<<"LOCAL_PLIC_KEP_RANGE_STAGE,exponent="<<exponent<<",operation="<<operation<<std::endl;
    };
    std::cout<<std::setprecision(std::numeric_limits<long double>::max_digits10)
        <<"LOCAL_PLIC_KEP_RANGE_FP,digits="<<std::numeric_limits<long double>::digits
        <<",min_exponent="<<std::numeric_limits<long double>::min_exponent<<",sizeof="<<sizeof(long double)
        <<",denorm_min="<<std::numeric_limits<long double>::denorm_min()<<",velocity_scale="<<scale<<std::endl;
    mark("actual native mesh");
    auto mesh=vofMaterialMesh(4);DictionaryMesh::RefineBox patch;
    patch.level=2;patch.min={.875,.875};patch.max={1,1};mesh->applyRefineBoxes({patch});
    const auto nativeDomain=domain(*mesh,false);std::vector<std::array<long double,2>> stocks(mesh->cellIndexNum());
    int cellIndex=-1;
    for(int ci=0;ci<mesh->cellIndexNum();++ci){const auto box=bounds(*mesh,ci);
        stocks[ci]={0,capacity(box,false)};if(box.x0==.75L&&box.y0==.25L)cellIndex=ci;}
    check(cellIndex>=0,"pure range probe lacks its original fine cell");
    mark("pure source publication");
    const auto publication=Publication::prepareLocalPlicInventory(*mesh,nativeDomain,Clock{},stocks,{});
    const auto sourceId=publication.identity();const auto* sourceStorage=publication.plicCells().data();
    const auto box=bounds(*mesh,cellIndex);const Point vertex{box.x1,box.y1};
    const long double amplitude=2*capacity(box,false)*scale;
    // Two adjacent REAL interior fine vertices keep endpoint-difference
    // circulation closed while overlapping transported U and material F.
    const Point neighbor{vertex.x,vertex.y+(box.y1-box.y0)};
    const auto psi=[&](Point p){return p.x==vertex.x&&p.y==vertex.y?amplitude:
        p.x==neighbor.x&&p.y==neighbor.y?amplitude/2:0.L;};
    std::vector<double> normal(publication.domain().nativeFaceCount());unsigned movingFaces=0;
    for(int row=0;row<publication.domain().nativeFaceCount();++row){const auto face=publication.domain().face(*mesh,row);
        const auto e=faceEndpoints(*mesh,face);const auto q=(face.horizontal?-1.L:1.L)*(psi(e[1])-psi(e[0]));
        if(face.area==0||face.outwardSign!=0){check(q==0,"pure range probe leaked through a physical wall");normal[row]=0;}
        else{normal[row]=static_cast<double>(q/face.area);movingFaces+=normal[row]!=0;}}
    check(movingFaces==7,"two-vertex range circulation lost a real nonzero native input");
    const auto originalNormal=normal;const Clock next{1,0,.125,Position::End};
    mark("actual native slab");const auto slab=Slab::prepare(*mesh,publication,next,normal,normal);
    mark("LocalPlic sweep");const auto transport=Transport::prepareLocalPlicSweep(*mesh,publication,slab,100,1);
    for(const auto& face:transport.faces())for(const auto* f:{&face.half,&face.full})
        check(f->liquid==0&&f->positiveLiquid==0&&f->negativeLiquid==0,"pure gas range probe manufactured a liquid current");
    mark("actual momentum geometry");const auto geometry=HydroMomentumGeometry::build(*mesh,vofMaterialMetrics(*mesh,publication));
    const auto before=HydroPhaseGeometry::build(*mesh,geometry.metrics(),publication);
    const auto middle=HydroPhaseGeometry::build(*mesh,geometry.metrics(),transport.midpoint());
    const auto after=HydroPhaseGeometry::build(*mesh,geometry.metrics(),transport.endpoint());
    mark("original material");const auto material=HydroMaterialTransport::prepareGeometricVof(*mesh,geometry,before,middle,after,transport);
    mark("P0 masses");
    const auto a=HydroMomentumMass::buildCellAverageProjection(*mesh,geometry,before,100,1);
    const auto m=HydroMomentumMass::buildCellAverageProjection(*mesh,geometry,middle,100,1);
    const auto b=HydroMomentumMass::buildCellAverageProjection(*mesh,geometry,after,100,1);
    mark("original dual material");const auto dual=HydroDualMaterial::prepare(*mesh,material,a,m,b);
    const auto& layout=geometry.layout();std::vector<int> walls;std::vector<double> u(layout.velocityCount());
    for(int i=0;i<layout.velocityCount();++i){u[i]=slab.beginNormal()[layout.nativeFace(i)];
        if(layout.domain().face(*mesh,layout.nativeFace(i)).outwardSign!=0)walls.push_back(i);}
    const HydroNativeValues values(layout,u,std::vector<double>(layout.pressureCount()));
    const auto sourceStamp=dual.sources();const auto clock=dual.timeLayer();
    // D13: authoritative receipts must remain intact when derived work RN0.
    const auto savedFaces=transport.faces();
    const auto savedDualCells=dual.cells();
    const auto savedNativeIntegrals=material.atomIntegrals(),savedDualIntegrals=dual.faceIntegrals(),savedDualRates=dual.faceRates();
    long double minF=std::numeric_limits<long double>::infinity(),maxF=0,minU=std::numeric_limits<long double>::infinity(),maxR=0;
    for(auto f:dual.faceRates())if(f!=0){minF=std::min(minF,std::abs(f));maxF=std::max(maxF,std::abs(f));}
    for(auto v:u)if(v!=0)minU=std::min(minU,std::abs(static_cast<long double>(v)));
    for(const auto& c:dual.cells())maxR=std::max(maxR,std::abs(c.integratedDefect));
    std::cout<<"LOCAL_PLIC_KEP_RANGE_INPUT,exponent="<<exponent<<",movingFaces="<<movingFaces
        <<",minimumNonzeroF="<<minF<<",maximumF="<<maxF<<",minimumNonzeroU="<<minU<<",maxAbsIntegratedR="<<maxR
        <<",dt="<<dual.deltaTime()<<std::endl;
    const auto law=HydroAdvectorBoundary::impermeable(*mesh,geometry,walls);
    const auto boundary=HydroConvectionBoundary::bind(*mesh,dual,std::vector<double>(dual.faceRates().size()),sourceStamp,clock);
    mark("buildImpermeable");const auto op=HydroKepConvection::buildImpermeable(*mesh,dual,law);
    unsigned supportedPairs=0;
    int coefficientMinExponent=std::numeric_limits<int>::max(),coefficientMaxExponent=std::numeric_limits<int>::min();
    int velocityMinExponent=std::numeric_limits<int>::max(),velocityMaxExponent=std::numeric_limits<int>::min();
    // The public frozen matrix is the actual same-owner closed physical
    // operator H. couplings() is a legacy correction table, not finalPairs.
    const auto physical=op.frozenMomentumOperator(*mesh,sourceStamp,clock);
    for(int outer=0;outer<physical.outerSize();++outer)
        for(HydroKineticHodge::Matrix::InnerIterator entry(physical,outer);entry;++entry)
            if(entry.row()<entry.col()&&entry.value()!=0&&u[entry.row()]!=0&&u[entry.col()]!=0){
        const int i=static_cast<int>(entry.row()),j=static_cast<int>(entry.col());
        const auto coefficient=static_cast<long double>(entry.value());
        const auto first=static_cast<long double>(u[i]),second=static_cast<long double>(u[j]);
        const int ce=std::ilogb(std::abs(coefficient)),ie=std::ilogb(std::abs(first)),je=std::ilogb(std::abs(second));
        coefficientMinExponent=std::min(coefficientMinExponent,ce);coefficientMaxExponent=std::max(coefficientMaxExponent,ce);
        velocityMinExponent=std::min({velocityMinExponent,ie,je});velocityMaxExponent=std::max({velocityMaxExponent,ie,je});
        if(supportedPairs<3)std::cout<<"LOCAL_PLIC_KEP_RANGE_SUPPORTED_PAIR,first="<<i<<",second="<<j
            <<",coefficient="<<coefficient<<",Ui="<<first<<",Uj="<<second<<std::endl;
        ++supportedPairs;
    }
    std::cout<<"LOCAL_PLIC_KEP_RANGE_SUPPORT,exponent="<<exponent<<",supportedPairs="<<supportedPairs
        <<",coefficient_min_binary_exponent="<<coefficientMinExponent<<",coefficient_max_binary_exponent="<<coefficientMaxExponent
        <<",velocity_min_binary_exponent="<<velocityMinExponent<<",velocity_max_binary_exponent="<<velocityMaxExponent<<std::endl;
    check(supportedPairs>0,"pure range probe lacks a public nonzero coefficient/Ui/Uj pair");
    // Independent integer-exponent absolute bound for COMPLETE work, not
    // the contraction helper or its multiplication algorithm. Nonzero
    // |x|<2^(ilogb|x|+1); no floating product is evaluated by this oracle.
    std::size_t powerTerms=0;
    int maximumPowerBoundExponent=std::numeric_limits<int>::min();
    const auto boundTerm=[&](int exponentBound){++powerTerms;
        maximumPowerBoundExponent=std::max(maximumPowerBoundExponent,exponentBound);};
    for(int outer=0;outer<physical.outerSize();++outer)
        for(HydroKineticHodge::Matrix::InnerIterator entry(physical,outer);entry;++entry)
            if(entry.value()!=0&&u[entry.row()]!=0&&u[entry.col()]!=0)
                boundTerm(std::ilogb(std::abs(entry.value()))+std::ilogb(std::abs(u[entry.row()]))
                    +std::ilogb(std::abs(u[entry.col()]))+3);
    const int dtExponent=std::ilogb(dual.deltaTime());
    for(int i=0;i<layout.velocityCount();++i)if(u[i]!=0){
        const int ue=std::ilogb(std::abs(u[i]));
        // H contains consistency diagonal: bound R and endpoint mass work
        // separately without changing either actual receipt.
        for(long double receipt:{dual.cells()[i].integratedDefect,dual.cells()[i].massChange})
            if(receipt!=0)boundTerm(2*ue+std::ilogb(std::abs(receipt))-dtExponent+3);
    }
    check(powerTerms>0,"actual complete-work oracle lacks physical factors");
    int countExponent=0;std::size_t countBound=1;
    while(countBound<powerTerms){countBound*=2;++countExponent;}
    const int absoluteSumBoundExponent=maximumPowerBoundExponent+countExponent;
    const int halfDenormExponent=std::numeric_limits<long double>::min_exponent
        -std::numeric_limits<long double>::digits-1;
    std::cout<<"LOCAL_PLIC_KEP_FINAL_POWER_BOUND,terms="<<powerTerms
        <<",absolute_sum_upper_binary_exponent="<<absoluteSumBoundExponent
        <<",half_denorm_binary_exponent="<<halfDenormExponent<<std::endl;
    check(absoluteSumBoundExponent<halfDenormExponent,
        "complete physical work is not certified below half a minimum subnormal");
    mark("apply");const auto result=op.apply(*mesh,values,boundary,sourceStamp,clock);
    for(const auto& c:result.report().components){
        check(std::isfinite(c.momentumIdentity)&&c.momentumIdentity<=1e-12L,"pure VOF range original momentum gate");
        check(std::isfinite(c.kineticIdentity)&&c.kineticIdentity<=1e-12L,"pure VOF range original kinetic power gate");
        // COMPLETE work bound certifies finite final RNE zero.
        check(std::isfinite(c.convectivePower)&&c.convectivePower==0,"actual final convective power must correctly round to zero");
        check(std::isfinite(c.materialKineticPower)&&c.materialKineticPower==0,"actual final material kinetic power must correctly round to zero");
        check(std::isfinite(c.kineticRate)&&c.kineticRate==0,"actual final load power must correctly round to zero");}
    check(material.atomIntegrals()==savedNativeIntegrals&&dual.faceIntegrals()==savedDualIntegrals&&dual.faceRates()==savedDualRates,
        "finite derived work changed authoritative material F");
    for(std::size_t i=0;i<savedFaces.size();++i)
        check(v4b1Fields(transport.faces()[i].half)==v4b1Fields(savedFaces[i].half)
            &&v4b1Fields(transport.faces()[i].full)==v4b1Fields(savedFaces[i].full),
            "finite derived work changed accepted parent phase net/gross F");
    for(std::size_t i=0;i<savedDualCells.size();++i)
        check(dual.cells()[i].integratedDefect==savedDualCells[i].integratedDefect
            &&dual.cells()[i].defectRate==savedDualCells[i].defectRate
            &&dual.cells()[i].massChange==savedDualCells[i].massChange
            &&dual.cells()[i].outwardIntegral==savedDualCells[i].outwardIntegral,
            "finite derived work changed actual material R/mass receipt");
    check(values.normalVelocity()==u,"finite derived work changed transported native U");
    check(publication.identity()==sourceId&&publication.plicCells().data()==sourceStorage&&normal==originalNormal,
        "pure VOF range changed accepted source owner or real native input");
    for(std::size_t ci=0;ci<stocks.size();++ci)
        check(publication.plicCells()[ci].inventory==stocks[ci],"pure VOF range changed independent source inventory");
    std::cout<<"LOCAL_PLIC_CONVECTION_RANGE_ONLY,exponent="<<exponent<<",PASS"<<std::endl;
}
// END D12 PURE VOF SMALL-VELOCITY RANGE PROBE


// BEGIN D14 INDEPENDENT INTEGER P0 PARTITION / ENDPOINT CONTRACT
struct D14P0Piece {int velocity=0;std::size_t overlap=0;std::uint64_t weight=0,prefix=0,whole=0;};
using D14P0Groups=std::map<std::pair<int,int>,std::vector<D14P0Piece>>;
unsigned d14FinitePartitions=0,d14ZeroPieces=0,d14PhaseChecks=0,d14EndpointDifferences=0,d14DualRows=0;
std::int64_t d14Coordinate(double x)
{
    const auto scaled=static_cast<long double>(x)*32;
    const auto integer=static_cast<std::int64_t>(scaled);
    check(scaled==integer,"D14 oracle requires its existing genuine 1/32 coordinate lattice");return integer;
}
std::uint64_t d14Weight(Vector2D lo,Vector2D hi,bool rz)
{
    const auto x0=d14Coordinate(lo.x),x1=d14Coordinate(hi.x),y0=d14Coordinate(lo.y),y1=d14Coordinate(hi.y);
    check(x0>=0&&x1>x0&&y1>y0,"D14 independent integer rectangle");
    // Common physical units cancel EXACTLY: XY 1/32^2, RZ pi/32^3.
    return static_cast<std::uint64_t>((rz?x1*x1-x0*x0:x1-x0)*(y1-y0));
}
std::uint64_t d14Units(long double stock)
{
    const auto unit=std::numeric_limits<long double>::denorm_min();
    check(stock>=0&&(stock==0||std::fpclassify(stock)==FP_SUBNORMAL),"D14 integer parent must be an actual subnormal stock");
    const auto n=static_cast<std::uint64_t>(stock/unit);
    check(n<=28&&static_cast<long double>(n)*unit==stock,"D14 integer parent units not exact");return n;
}
std::uint64_t d14Cumulative(std::uint64_t n,std::uint64_t prefix,std::uint64_t whole)
{
    check(whole>0&&prefix<=whole&&n<=28,"D14 integer CDF domain");
    // Independent integer quotient/remainder round-to-nearest, ties-to-even.
    const auto product=n*prefix,quotient=product/whole,remainder=product%whole;
    return quotient+((2*remainder>whole||(2*remainder==whole&&(quotient&1)))?1:0);
}
D14P0Groups d14Groups(const QuadtreeMesh& mesh,const HydroMomentumGeometry& geometry)
{
    const auto& layout=geometry.layout();const bool rz=layout.domain().axisymmetric();D14P0Groups groups;
    for(int u=0;u<layout.velocityCount();++u){const int component=layout.normalAxis(u)==HydroDofLayout::NormalAxis::First?0:1;
        const auto& overlaps=geometry.volumes()[u].overlaps;
        for(std::size_t k=0;k<overlaps.size();++k){const auto& o=overlaps[k];groups[{o.pressure,component}].push_back({u,k,d14Weight(o.lower,o.upper,rz),0,0});}}
    for(auto& item:groups){const int pressure=item.first.first,component=item.first.second;
        auto& pieces=item.second;
        const auto key=[&](const D14P0Piece& piece){const auto& o=geometry.volumes()[piece.velocity].overlaps[piece.overlap];
            return component==0?std::make_tuple(o.lower.x,o.upper.x,o.lower.y,o.upper.y):
                std::make_tuple(o.lower.y,o.upper.y,o.lower.x,o.upper.x);};
        std::sort(pieces.begin(),pieces.end(),[&](const auto& a,const auto& b){return key(a)<key(b);});
        const auto box=bounds(mesh,layout.nativeCell(pressure));
        const auto whole=d14Weight({static_cast<double>(box.x0),static_cast<double>(box.y0)},
            {static_cast<double>(box.x1),static_cast<double>(box.y1)},rz);
        std::uint64_t prefix=0;
        for(std::size_t k=0;k<pieces.size();++k){if(k)check(key(pieces[k-1])!=key(pieces[k]),"D14 duplicate actual rectangle");
            pieces[k].prefix=prefix;pieces[k].whole=whole;prefix+=pieces[k].weight;}
        check(prefix==whole&&pieces.size()>=2&&pieces.size()<=4,"D14 true component partition does not exactly cover whole cell");
    }
    return groups;
}
void d14InspectP0(const QuadtreeMesh& mesh,const HydroMomentumMass& mass,unsigned scene,int layer)
{
    const auto& layout=mass.layout();const auto& cells=mass.phaseGeometry().vofPublication().plicCells();
    const auto groups=d14Groups(mesh,mass.geometry());unsigned finite=0,zeros=0;
    for(const auto& item:groups){const int pressure=item.first.first,component=item.first.second;
        for(int phase=0;phase<2;++phase){const auto parent=cells[layout.nativeCell(pressure)].inventory[phase];
            const bool exceptional=mass.usesFiniteP0Partition(pressure,component,phase);long double sum=0;
            if(exceptional){check(parent>0,"D14 finite law flagged an empty original stock");++finite;}
            for(const auto& piece:item.second){const auto& o=mass.geometry().volumes()[piece.velocity].overlaps[piece.overlap];
                const auto m=mass.representedOverlapPhaseVolumes()[piece.velocity][piece.overlap];const auto actual=phase?m.gas:m.liquid;
                check(std::isfinite(actual)&&actual>=0,"D14 represented phase must be finite and nonnegative");
                if(exceptional){const auto n=d14Units(parent);
                    const auto before=d14Cumulative(n,piece.prefix,piece.whole),after=d14Cumulative(n,piece.prefix+piece.weight,piece.whole);
                    const auto expected=static_cast<long double>(after-before)*std::numeric_limits<long double>::denorm_min();
                    check(actual==expected,"D14 actual phase moment differs from independent integer prefix-CDF difference");
                    if(actual==0)++zeros;
                }else check(actual==o.cellFraction*parent,"D14 ordinary phase projection changed its original multiplication");
                sum+=actual;++d14PhaseChecks;
            }
            if(parent==0)check(sum==0,"D14 empty source acquired a projected phase");
            else near(sum,parent,parent,"D14 independent per-component own-phase addback",1e-12L);
        }
    }
    for(const auto& component:mass.report().maxRelativeP0OwnPhaseCellPartition)for(auto error:component)
        check(std::isfinite(error)&&error<=1e-12L,"D14 original own-phase projection gate");
    check(mass.report().finiteP0PhasePartitions==finite&&mass.report().finiteP0ZeroPhasePieces==zeros,"D14 finite partition diagnostic detached from actual public pieces");
    for(double value:mass.masses())check(std::isfinite(value)&&value>0,"D14 stored physical mass lost original positivity");
    d14FinitePartitions+=finite;d14ZeroPieces+=zeros;
    std::cout<<"LOCAL_PLIC_FINITE_P0_MASS,scene="<<scene<<",layer="<<layer<<",partitions="<<finite<<",zeros="<<zeros<<std::endl;
}
void d14InspectDual(const QuadtreeMesh& mesh,const HydroDualMaterial& dual,unsigned scene)
{
    ConservativeVofMomentumContract::qualify(mesh,dual);
    const auto& geometry=dual.geometry();const auto& layout=geometry.layout();
    const auto& begin=dual.beginMass();const auto& end=dual.endMass();
    const auto& a=begin.phaseGeometry().vofPublication().plicCells();const auto& b=end.phaseGeometry().vofPublication().plicCells();
    const auto groups=d14Groups(mesh,geometry);unsigned endpointDifferences=0;
    for(const auto& item:groups){const int pressure=item.first.first,component=item.first.second,native=layout.nativeCell(pressure);
        for(int phase=0;phase<2;++phase)if(begin.usesFiniteP0Partition(pressure,component,phase)||end.usesFiniteP0Partition(pressure,component,phase)){
            const auto na=d14Units(a[native].inventory[phase]),nb=d14Units(b[native].inventory[phase]);
            for(const auto& piece:item.second){
                const auto pieceUnits=[&](std::uint64_t n){return d14Cumulative(n,piece.prefix+piece.weight,piece.whole)-d14Cumulative(n,piece.prefix,piece.whole);};
                const auto endpointDelta=static_cast<std::int64_t>(pieceUnits(nb))-static_cast<std::int64_t>(pieceUnits(na));
                const auto deltaUnits=nb>=na?nb-na:na-nb;
                const auto deltaFirst=(nb>=na?1:-1)*static_cast<std::int64_t>(pieceUnits(deltaUnits));
                if(endpointDelta!=deltaFirst)++endpointDifferences;
                const auto ma=begin.representedOverlapPhaseVolumes()[piece.velocity][piece.overlap];
                const auto mb=end.representedOverlapPhaseVolumes()[piece.velocity][piece.overlap];
                const auto actual=(phase?mb.gas:mb.liquid)-(phase?ma.gas:ma.liquid);
                check(actual==static_cast<long double>(endpointDelta)*std::numeric_limits<long double>::denorm_min(),
                    "D14 real endpoint phase difference disagrees with independent integer oracle");
            }
        }
    }
    for(int u=0;u<layout.velocityCount();++u){const int component=layout.normalAxis(u)==HydroDofLayout::NormalAxis::First?0:1;
        long double expected=0,gross=0;
        for(std::size_t k=0;k<geometry.volumes()[u].overlaps.size();++k){const auto& o=geometry.volumes()[u].overlaps[k];
            const int native=layout.nativeCell(o.pressure);
            const auto ma=begin.representedOverlapPhaseVolumes()[u][k],mb=end.representedOverlapPhaseVolumes()[u][k];
            const std::array<long double,2> pa{{ma.liquid,ma.gas}},pb{{mb.liquid,mb.gas}};
            const std::array<bool,2> exceptional{{begin.usesFiniteP0Partition(o.pressure,component,0)||end.usesFiniteP0Partition(o.pressure,component,0),
                begin.usesFiniteP0Partition(o.pressure,component,1)||end.usesFiniteP0Partition(o.pressure,component,1)}};
            if(!exceptional[0]&&!exceptional[1]){
                const auto dl=dual.material().liquidDensity()*(b[native].inventory[0]-a[native].inventory[0]);
                const auto dg=dual.material().gasDensity()*(b[native].inventory[1]-a[native].inventory[1]);
                expected+=o.cellFraction*(dl+dg);gross+=o.cellFraction*(std::abs(dl)+std::abs(dg));
            }else for(int phase=0;phase<2;++phase){
                const auto rho=phase?dual.material().gasDensity():dual.material().liquidDensity();
                const auto term=exceptional[phase]?rho*(pb[phase]-pa[phase]):
                    o.cellFraction*(rho*(b[native].inventory[phase]-a[native].inventory[phase]));
                expected+=term;gross+=std::abs(term);}
        }
        const auto& row=dual.cells()[u];
        if(gross==0)check(row.massChange==0,"D14 unchanged independent endpoint receipts acquired canonical deltaM");
        else near(row.massChange,expected,gross,"D14 canonical deltaM consumes coherent real endpoint phase receipts",1e-12L);
        check(row.integratedDefect==row.massChange+row.outwardIntegral,"D14 canonical R no longer uses actual deltaM plus unchanged F");
        check(row.geometryTransferChange==row.massChange-row.mappedPrimalMassChange,"D14 represented metric difference hidden from actual receipt");
        ++d14DualRows;
    }
    check(dual.report().maxDefectDecomposition<=256.L*std::numeric_limits<double>::epsilon(),"D14 original dual decomposition gate");
    for(const auto& c:dual.report().components)check(c.balanceIdentity<=256.L*std::numeric_limits<double>::epsilon(),"D14 original component balance gate");
    d14EndpointDifferences+=endpointDifferences;
    std::cout<<"LOCAL_PLIC_FINITE_P0_ENDPOINTS,scene="<<scene<<",nonlinear_receipts="<<endpointDifferences<<std::endl;
}
// END D14 INDEPENDENT INTEGER P0 PARTITION / ENDPOINT CONTRACT

void localPlicAtomsOnly(bool convection=false,bool physicalDensity=false)
{
    const double rhoLiquid=physicalDensity?1000.:100.,rhoGas=physicalDensity?1.226:1.;
    if(physicalDensity)std::cout<<std::setprecision(std::numeric_limits<double>::max_digits10)
        <<"LOCAL_PLIC_PHYSICAL_DENSITIES,rhoL="<<rhoLiquid<<",rhoG="<<rhoGas<<std::endl;
    const long double unit=std::numeric_limits<long double>::denorm_min();
    std::cout<<std::setprecision(std::numeric_limits<long double>::max_digits10)
        <<"LOCAL_PLIC_ATOM_FP,digits="<<std::numeric_limits<long double>::digits
        <<",min_exponent="<<std::numeric_limits<long double>::min_exponent
        <<",sizeof="<<sizeof(long double)<<",denorm_min="<<unit<<std::endl;
    check(unit>0&&std::fpclassify(unit)==FP_SUBNORMAL,"actual C++ minimum subnormal is unavailable");
    struct Scene {bool rz,gas,reverse;int stockUnits,witnessExtent;};
    // Odd FULL parent 7 comes first: independent half-face rounding gives
    // two nonzero pieces 4+4. The ORIGINAL material consumer runs before
    // our addback assertions, proving its real unchanged production guard.
    const std::array<Scene,8> scenes{{
        {false,false,false,28,1},{true,false,false,28,1},
        {false,true,true,28,1},{true,true,true,28,1},
        {false,false,true,4,1},{true,true,false,4,1},
        {false,true,false,8,0},{true,false,true,8,0}}};
    unsigned sceneNumber=0,quantizedZeros=0,oddWitnesses=0,minimumWitnesses=0,phaseChecks=0,interiorChecks=0;
    unsigned convectionScenes=0;
    for(const auto scene:scenes){
        stage="D11 actual native atoms scene "+std::to_string(sceneNumber);
        auto mesh=vofMaterialMesh(4);
        // A real remote level-2 patch keeps this donor's fine one-ring but
        // makes the mesh 0/1/2, like the deep-AMR failure. The general native
        // chart then supplies its real midpoint atoms; the special two-level
        // XY virtual tensor otherwise subdivides this face into quarters.
        DictionaryMesh::RefineBox deepPatch;deepPatch.level=2;
        deepPatch.min={.875,.875};deepPatch.max={1,1};mesh->applyRefineBoxes({deepPatch});
        const auto nativeDomain=domain(*mesh,scene.rz);
        check(!nativeDomain.uniformLeaves(),"subnormal fixture lost its actual 2:1 mesh");
        unsigned level01=0,level12=0;
        for(const auto face:nativeDomain.faces(*mesh))if(face.negativeCell>=0&&face.positiveCell>=0){
            const int a=mesh->cells(mesh->cellIndexToLabel[face.negativeCell]).level();
            const int b=mesh->cells(mesh->cellIndexToLabel[face.positiveCell]).level();
            check(std::abs(a-b)<=1,"subnormal actual shared native face violates 2:1");
            if(std::min(a,b)==0&&std::max(a,b)==1)++level01;
            if(std::min(a,b)==1&&std::max(a,b)==2)++level12;
        }
        check(level01>0&&level12>0,"subnormal fixture lacks actual 0:1 and 1:2 shared native faces");
        std::cout<<"LOCAL_PLIC_ATOM_NATIVE_BALANCE,level01="<<level01<<",level12="<<level12<<std::endl;
        int microCell=-1;std::vector<std::array<long double,2>> stocks(mesh->cellIndexNum());
        for(int ci=0;ci<mesh->cellIndexNum();++ci){const auto b=bounds(*mesh,ci);const auto c=capacity(b,scene.rz);
            stocks[ci]=scene.gas?std::array<long double,2>{{c,0}}:std::array<long double,2>{{0,c}};
            if(b.x0==.75L&&b.y0==.25L)microCell=ci;}
        check(microCell>=0,"subnormal fixture lacks its real interior fine leaf");
        const auto microBounds=bounds(*mesh,microCell);const auto c=capacity(microBounds,scene.rz);
        const auto microscopic=static_cast<long double>(scene.stockUnits)*unit;
        stocks[microCell]=scene.gas?std::array<long double,2>{{c,microscopic}}:std::array<long double,2>{{microscopic,c}};
        const auto savedStocks=stocks;
        const auto source=Publication::prepareLocalPlicInventory(*mesh,nativeDomain,Clock{},stocks,{});
        const auto sourceId=source.identity();const auto* sourceStorage=source.plicCells().data();
        check(!source.plicCells()[microCell].resolved&&source.plicCells()[microCell].inventory[scene.gas?1:0]==microscopic,
            "subnormal fixture did not retain its independent unresolved micro inventory");
        // One REAL interior vertex streamfunction gives a closed circulation.
        // q=2*capacity, dt=1/8: target full sweep is a quarter cell. Choosing
        // its lower vertex reverses the micro donor's canonical x flux.
        const Point vertex=scene.reverse?Point{microBounds.x0,microBounds.y0}:Point{microBounds.x1,microBounds.y1};
        const long double amplitude=2*c;
        const auto psi=[&](Point q){return q.x==vertex.x&&q.y==vertex.y?amplitude:0.L;};
        std::vector<double> normal(source.domain().nativeFaceCount());int targetRow=-1;
        for(int row=0;row<source.domain().nativeFaceCount();++row){const auto face=source.domain().face(*mesh,row);
            const auto e=faceEndpoints(*mesh,face);const auto q=(face.horizontal?-1.L:1.L)*(psi(e[1])-psi(e[0]));
            if(face.area==0||face.outwardSign!=0){check(q==0,"subnormal closed circulation leaked through a physical boundary");normal[row]=0;}
            else normal[row]=static_cast<double>(q/face.area);
            if(!face.horizontal&&(scene.reverse?face.positiveCell:face.negativeCell)==microCell
                &&(scene.reverse?q<0:q>0))targetRow=row;}
        check(targetRow>=0,"subnormal donor lacks its actual first-direction outflow face");
        const Clock next{1,0,.125,Position::End};const auto savedNormal=normal;
        const auto slab=Slab::prepare(*mesh,source,next,normal,normal);
        const auto transport=Transport::prepareLocalPlicSweep(*mesh,source,slab,rhoLiquid,rhoGas);
        const auto& target=transport.faces()[targetRow];const auto& witness=scene.witnessExtent?target.full:target.half;
        const auto minority=[&](const GeometricVofIntegratedFlux& f){return scene.gas?f.gas:f.liquid;};
        const auto expectedUnits=scene.stockUnits==28?7.L:1.L;
        check(std::abs(minority(witness))==expectedUnits*unit&&std::fpclassify(minority(witness))==FP_SUBNORMAL,
            "actual transport did not produce the declared minimum/odd-subnormal parent");
        check(scene.reverse?minority(witness)<0:minority(witness)>0,"subnormal donor direction was not exercised");
        const auto geometry=HydroMomentumGeometry::build(*mesh,vofMaterialMetrics(*mesh,source));
        const auto reconstruction=HydroAdvectorReconstruction::prepare(*mesh,geometry);
        const auto& layout=geometry.layout();int targetVelocity=-1;
        for(int u=0;u<layout.velocityCount();++u)if(layout.nativeFace(u)==targetRow)targetVelocity=u;
        check(targetVelocity>=0,"subnormal face lost its actual velocity index");
        const auto& targetAtoms=reconstruction.primalAtoms()[targetVelocity];
        const auto targetEnds=faceEndpoints(*mesh,source.domain().face(*mesh,targetRow));
        std::cout<<"LOCAL_PLIC_ATOM_LAYOUT,row="<<targetRow<<",count="<<targetAtoms.size()
            <<",from="<<targetEnds[0].x<<','<<targetEnds[0].y<<",to="<<targetEnds[1].x<<','<<targetEnds[1].y;
        for(int id:targetAtoms){const auto& a=reconstruction.atoms()[id];
            std::cout<<",atom="<<id<<':'<<a.begin.x<<','<<a.begin.y<<':'<<a.end.x<<','<<a.end.y;}
        std::cout<<std::endl;
        check(targetAtoms.size()==2,"selected fine face no longer has two actual primal atoms");
        std::vector<GeometricVofNativeSubface> nativeQueries;
        std::vector<GeometricVofInteriorSubface> interiorQueries;
        std::vector<int> nativeQueryByAtom(reconstruction.atoms().size(),-1),interiorQueryByAtom(reconstruction.atoms().size(),-1);
        for(std::size_t id=0;id<reconstruction.atoms().size();++id){const auto& atom=reconstruction.atoms()[id];
            const Point a{atom.begin.x,atom.begin.y},b{atom.end.x,atom.end.y};
            if(atom.primalVelocity>=0){nativeQueryByAtom[id]=static_cast<int>(nativeQueries.size());
                nativeQueries.push_back({layout.nativeFace(atom.primalVelocity),a,b});}
            else{check(atom.ownerPressure>=0,"actual interior atom lacks its original pressure owner");
                interiorQueryByAtom[id]=static_cast<int>(interiorQueries.size());
                interiorQueries.push_back({layout.nativeCell(atom.ownerPressure),a,b});}}
        const auto printParent=[&](const char* extent,const GeometricVofIntegratedFlux& f){
            std::cout<<"LOCAL_PLIC_ATOM_PARENT,"<<sceneNumber<<",extent="<<extent
                <<",total="<<f.total<<",rawLiquid="<<f.rawLiquid<<",liquid="<<f.liquid<<",gas="<<f.gas
                <<",positiveLiquid="<<f.positiveLiquid<<",negativeLiquid="<<f.negativeLiquid
                <<",positiveGas="<<f.positiveGas<<",negativeGas="<<f.negativeGas
                <<",minority_units="<<minority(f)/unit<<std::endl;
        };
        printParent("half",target.half);printParent("full",target.full);
        const auto restricted=transport.restrictNativeFaces(*mesh,nativeQueries);
        std::cout<<"LOCAL_PLIC_ATOM_WITNESS,"<<sceneNumber<<",rz="<<scene.rz<<",minority="<<(scene.gas?"G":"L")
            <<",negative="<<scene.reverse<<",extent="<<scene.witnessExtent<<",parent_units="<<minority(witness)/unit;
        for(int id:targetAtoms){const auto qi=nativeQueryByAtom[id];check(qi>=0,"primal atom lacks its own native query");
            const auto& atom=reconstruction.atoms()[id];const auto e=faceEndpoints(*mesh,source.domain().face(*mesh,targetRow));
            const auto length=std::hypot(static_cast<long double>(atom.end.x)-atom.begin.x,static_cast<long double>(atom.end.y)-atom.begin.y);
            const auto whole=std::hypot(e[1].x-e[0].x,e[1].y-e[0].y);
            check(2*length==whole,"selected actual primal atom is not a half-face");
            std::cout<<",piece_units="<<minority(scene.witnessExtent?restricted[qi].full:restricted[qi].half)/unit;}
        std::cout<<std::endl;
        const auto before=HydroPhaseGeometry::build(*mesh,geometry.metrics(),source);
        const auto middle=HydroPhaseGeometry::build(*mesh,geometry.metrics(),transport.midpoint());
        const auto after=HydroPhaseGeometry::build(*mesh,geometry.metrics(),transport.endpoint());
        // The genuine old guard must be the first addback failure, not a
        // copied production formula or an assertion that all pieces survive.
        const auto material=HydroMaterialTransport::prepareGeometricVof(*mesh,geometry,before,middle,after,transport);
        check(material.vofTransport().sameInstance(transport),"subnormal material lost its accepted flux owner");
        check(material.liquidDensity()==rhoLiquid&&material.gasDensity()==rhoGas,"actual micro consumer densities differ from selected transport density");
        // BEGIN D12 REAL P0 DUAL / IMPERMEABLE KEP CONSUMER
        if(convection){
            const auto beginMass=HydroMomentumMass::buildCellAverageProjection(*mesh,geometry,before,rhoLiquid,rhoGas);
            const auto middleMass=HydroMomentumMass::buildCellAverageProjection(*mesh,geometry,middle,rhoLiquid,rhoGas);
            const auto endMass=HydroMomentumMass::buildCellAverageProjection(*mesh,geometry,after,rhoLiquid,rhoGas);
            d14InspectP0(*mesh,beginMass,sceneNumber,0);d14InspectP0(*mesh,middleMass,sceneNumber,1);d14InspectP0(*mesh,endMass,sceneNumber,2);
            const auto dual=HydroDualMaterial::prepare(*mesh,material,beginMass,middleMass,endMass);
            d14InspectDual(*mesh,dual,sceneNumber);
            const auto sourceStamp=dual.sources();const auto midpointClock=dual.timeLayer();
            const auto savedIntegrals=dual.faceIntegrals(),savedRates=dual.faceRates();
            std::vector<int> walls;std::vector<double> nativeU(layout.velocityCount());
            for(int u=0;u<layout.velocityCount();++u){nativeU[u]=slab.beginNormal()[layout.nativeFace(u)];
                if(layout.domain().face(*mesh,layout.nativeFace(u)).outwardSign!=0)walls.push_back(u);}
            const HydroNativeValues values(layout,nativeU,std::vector<double>(layout.pressureCount()));
            const auto law=HydroAdvectorBoundary::impermeable(*mesh,geometry,walls);
            unsigned subnormalIntegrals=0,subnormalRates=0,printed=0;
            long double minimumIntegral=std::numeric_limits<long double>::infinity();
            for(std::size_t f=0;f<savedIntegrals.size();++f){const auto integral=savedIntegrals[f],rate=savedRates[f];
                if(integral!=0)minimumIntegral=std::min(minimumIntegral,std::abs(integral));
                const bool tinyIntegral=std::fpclassify(integral)==FP_SUBNORMAL;
                const bool tinyRate=std::fpclassify(rate)==FP_SUBNORMAL;
                subnormalIntegrals+=tinyIntegral;subnormalRates+=tinyRate;
                if((tinyIntegral||tinyRate)&&printed<8){const auto& segment=geometry.segments()[f];
                    std::cout<<"LOCAL_PLIC_KEP_F,"<<sceneNumber<<",segment="<<f
                        <<",integral="<<integral<<",rate="<<rate
                        <<",negativeVelocity="<<segment.negativeVelocity<<",positiveVelocity="<<segment.positiveVelocity;
                    if(segment.negativeVelocity>=0)std::cout<<",negativeU="<<nativeU[segment.negativeVelocity];
                    if(segment.positiveVelocity>=0)std::cout<<",positiveU="<<nativeU[segment.positiveVelocity];
                    std::cout<<std::endl;++printed;}}
            std::cout<<"LOCAL_PLIC_KEP_PREPARE,"<<sceneNumber<<",subnormalIntegrals="<<subnormalIntegrals
                <<",subnormalRates="<<subnormalRates<<",minimumNonzeroIntegral="<<minimumIntegral
                <<",dt="<<dual.deltaTime()<<std::endl;
            const auto op=HydroKepConvection::buildImpermeable(*mesh,dual,law);
            check(op.isImpermeable()&&op.fixedVelocityDofs()==walls&&op.periodicPairs().empty(),
                "actual VOF convection lost its typed native impermeable law");
            check(&op.material().faceRates()==&dual.faceRates(),"actual VOF convection replaced unique material F");
            unsigned subnormalCouplings=0,observedSubnormalProducts=0,observedRoundedZeroProducts=0;
            for(const auto& pair:op.couplings()){
                if(std::fpclassify(pair.coefficient)==FP_SUBNORMAL)++subnormalCouplings;
                // Passive PUBLIC coefficient/sample diagnostics, not the
                // production checked-product implementation or its oracle.
                for(int u:{pair.first,pair.second}){const auto velocity=static_cast<long double>(nativeU[u]);
                    const auto sampleProduct=pair.coefficient*velocity;
                    if(std::fpclassify(sampleProduct)==FP_SUBNORMAL)++observedSubnormalProducts;
                    if(pair.coefficient!=0&&velocity!=0&&sampleProduct==0)++observedRoundedZeroProducts;}}
            std::cout<<"LOCAL_PLIC_KEP_APPLY,"<<sceneNumber<<",subnormalPublicCouplings="<<subnormalCouplings
                <<",observedSubnormalCoefficientTimesU="<<observedSubnormalProducts
                <<",observedRoundedZeroCoefficientTimesU="<<observedRoundedZeroProducts<<std::endl;
            const auto boundary=HydroConvectionBoundary::bind(*mesh,dual,
                std::vector<double>(dual.faceRates().size()),sourceStamp,midpointClock);
            const auto result=op.apply(*mesh,values,boundary,sourceStamp,midpointClock);
            check(&result.material().faceRates()==&dual.faceRates(),"actual VOF convection result replaced the original material owner");
            for(std::size_t u=0;u<result.momentumLoads().size();++u)check(result.momentumLoads()[u]==-result.momentumRates()[u],"finite-CDF conservative load retained applied uR");
            for(const auto& component:result.report().components){
                check(std::isfinite(component.momentumIdentity)&&component.momentumIdentity<=1e-12L,
                    "actual VOF convection original momentum identity failed");
                check(std::isfinite(component.kineticIdentity)&&component.kineticIdentity<=1e-12L,
                    "actual VOF convection original kinetic power identity failed");}
            check(dual.faceIntegrals()==savedIntegrals&&dual.faceRates()==savedRates
                &&values.normalVelocity()==nativeU,"actual VOF convection changed unique F or actual input U");
            dual.validateFor(*mesh,sourceStamp,midpointClock);++convectionScenes;
            std::cout<<"LOCAL_PLIC_KEP_CONSUMER,"<<sceneNumber<<",PASS"<<std::endl;
        }
        // END D12 REAL P0 DUAL / IMPERMEABLE KEP CONSUMER
        const auto interior=transport.measureInteriorFaces(*mesh,interiorQueries);
        const auto payload=[&](const GeometricVofIntegratedFlux& f){
            if(physicalDensity)check(f.mass==material.liquidDensity()*f.liquid+material.gasDensity()*f.gas,
                "physical-density atom mass detached from actual signed phases");
            else{
            check(f.mass==100*f.liquid+f.gas,"subnormal atom mass detached from signed phases");
            }
            for(int phase=0;phase<2;++phase){const auto pos=phase?f.positiveGas:f.positiveLiquid;
                const auto neg=phase?f.negativeGas:f.negativeLiquid,net=phase?f.gas:f.liquid;
                check(pos>=0&&neg>=0,"subnormal atom acquired negative gross");
                if(std::fpclassify(pos)==FP_SUBNORMAL||std::fpclassify(neg)==FP_SUBNORMAL||std::fpclassify(net)==FP_SUBNORMAL)
                    check(net==pos-neg,"subnormal signed phase is not its own routed gross difference");}};
        for(const auto& q:restricted){payload(q.half);payload(q.full);}
        for(const auto& q:interior){payload(q.half);payload(q.full);++interiorChecks;}
        for(int u=0;u<layout.velocityCount();++u)for(int extent=0;extent<2;++extent){
            const auto& parent=extent?transport.faces()[layout.nativeFace(u)].full:transport.faces()[layout.nativeFace(u)].half;
            std::array<long double,12> sum{};
            for(int id:reconstruction.primalAtoms()[u]){const auto qi=nativeQueryByAtom[id];check(qi>=0,"native partition lost a primal atom");
                const auto values=v4b1Fields(extent?restricted[qi].full:restricted[qi].half);
                for(std::size_t k=0;k+1<sum.size();++k)sum[k]+=values[k];}
            const auto values=v4b1Fields(parent);
            for(std::size_t k:{4U,5U,7U,8U,9U,10U}){const auto scale=k==4?parent.positiveLiquid+parent.negativeLiquid:
                    k==5?parent.positiveGas+parent.negativeGas:std::abs(values[k]);
                if(scale==0)check(sum[k]==0,"zero native phase acquired an atom contribution");
                else near(sum[k],values[k],scale,"actual native atom own-phase addback",1e-12L);++phaseChecks;}}
        for(std::size_t id=0;id<reconstruction.atoms().size();++id){const auto qi=nativeQueryByAtom[id];
            const auto expected=qi>=0?restricted[qi].full.mass:interior[interiorQueryByAtom[id]].full.mass;
            check(material.atomIntegrals()[id]==expected,"material atom table differs from original public measurements");}
        for(int id:targetAtoms){const auto qi=nativeQueryByAtom[id];const auto part=minority(scene.witnessExtent?restricted[qi].full:restricted[qi].half);
            if(expectedUnits==7)check(part!=0,"odd-subnormal witness requires two nonzero actual pieces");
            if(expectedUnits==1&&part==0)++quantizedZeros;}
        if(expectedUnits==7)++oddWitnesses;else ++minimumWitnesses;
        for(const auto* endpoint:{&transport.midpoint(),&transport.endpoint()}){long double minorTotal=0;
            for(const auto& cell:endpoint->plicCells())minorTotal+=cell.inventory[scene.gas?1:0];
            near(minorTotal,microscopic,microscopic,"closed subnormal inventory remains conserved",1e-12L);}
        check(source.identity()==sourceId&&source.plicCells().data()==sourceStorage&&normal==savedNormal,
            "subnormal prepare changed accepted source or velocity input");
        for(std::size_t ci=0;ci<savedStocks.size();++ci)
            check(source.plicCells()[ci].inventory==savedStocks[ci],"subnormal prepare changed independent source stocks");
        ++sceneNumber;
    }
    check(sceneNumber==8&&oddWitnesses==4&&minimumWitnesses==4&&quantizedZeros>0&&interiorChecks>0,"subnormal atom coverage incomplete");
    std::cout<<"LOCAL_PLIC_NATIVE_ATOMS_ONLY,scenes="<<sceneNumber<<",odd="<<oddWitnesses<<",minimum="<<minimumWitnesses
        <<",quantizedZeros="<<quantizedZeros<<",phaseChecks="<<phaseChecks<<",interior="<<interiorChecks
        <<",ownPhaseTolerance=1e-12,PASS"<<std::endl;
    if(convection){check(d14FinitePartitions>0&&d14ZeroPieces>0&&d14EndpointDifferences>0&&d14DualRows>0,"D14 actual finite partition/zero/nonlinear endpoint coverage incomplete");
        std::cout<<"LOCAL_PLIC_FINITE_P0_PASS,partitions="<<d14FinitePartitions<<",zeros="<<d14ZeroPieces<<",phaseChecks="<<d14PhaseChecks
            <<",nonlinear_receipts="<<d14EndpointDifferences<<",dualRows="<<d14DualRows<<std::endl;}
    if(convection){check(convectionScenes==scenes.size(),"actual VOF convection did not consume every original subnormal scene");
        std::cout<<"LOCAL_PLIC_CONVECTION_ONLY,scenes="<<convectionScenes<<",PASS"<<std::endl;}
}
// END D11 ACTUAL NATIVE ATOM SUBNORMAL CONTRACT


// BEGIN CONDITIONAL RZ C2 REAL CAPILLARY STEP CONTRACT
void localPlicRzCapillaryStepOnly()
{
    stage="C2 actual RZ sphere finite initialization";
    DictionaryMesh dictionary;dictionary.Nx=dictionary.Ny=32;dictionary.maxCellNum=16384;dictionary.maxInitialIterations=0;
    auto mesh=std::make_unique<QuadtreeMesh>(dictionary);
    GeometricVofInitializationConfig initial;initial.subcellsPerAxis=8;
    const auto source=GeometricVofInitialization::prepare(*mesh,true,[](Vector2D p){
        return p.x*p.x+(p.y-.5)*(p.y-.5)-.2*.2;},initial);
    const auto sourceId=source.identity();const auto* sourceStorage=source.plicCells().data();const auto saved=source.plicCells();
    const auto metrics=vofMaterialMetrics(*mesh,source);const auto& layout=metrics.layout();
    check(layout.domain().axisymmetric()&&layout.pressureCount()==1024,"C2 sphere lost its actual32x32 RZ native domain");
    const auto phase=HydroPhaseGeometry::build(*mesh,metrics,source);
    const HydroNativeValues zero(layout,std::vector<double>(layout.velocityCount()),std::vector<double>(layout.pressureCount()));
    auto accepted=HydroTwoPhaseState::createVof(*mesh,phase,zero);const auto oldState=accepted;
    HydroVofAcConfig config;config.liquidDensity=1000;config.gasDensity=1.226;
    config.liquidViscosity=config.gasViscosity=0;config.surfaceTension=.7;
    config.momentumDensityRepresentation=HydroMomentumMass::DensityRepresentation::CellAverageProjection;
    long double h=std::numeric_limits<long double>::infinity();std::size_t mixed=0,axisMixed=0;
    std::array<long double,2> oldTotal{};
    for(const auto& cell:saved){
        for(int phase=0;phase<2;++phase)oldTotal[phase]+=cell.inventory[phase];
        if(cell.inventory[0]>0&&cell.inventory[1]>0){++mixed;axisMixed+=cell.rectangle.lower.x==0;
            h=std::min(h,static_cast<long double>(std::min(cell.rectangle.upper.x-cell.rectangle.lower.x,cell.rectangle.upper.y-cell.rectangle.lower.y)));}
    }
    check(mixed>0&&axisMixed>=2&&h==1.L/32,"C2 genuine RZ source lacks mixed/axis native leaf scales");
    stage="C2 active RZ capillary stability rate";
    const auto rates=HydroVofStep::stabilityRates(*mesh,accepted,config);
    const auto mean=(static_cast<long double>(config.liquidDensity)+config.gasDensity)/2;
    const auto expectedRate=std::sqrt(std::acos(-1.L)*config.surfaceTension/(mean*h*h*h));
    check(rates.capillaryAvailable&&rates.interfaceMinimumSpacing==h&&rates.flowRate==0,
        "C2 RZ active capillary scale must use actual mixed source spacing, not effective HF spacing");
    near(rates.capillaryRate,expectedRate,expectedRate,"C2 independent RZ capillary inverse time",1e-12L);
    auto zeroSigma=config;zeroSigma.surfaceTension=0;
    const auto noCapillary=HydroVofStep::stabilityRates(*mesh,accepted,zeroSigma);
    check(noCapillary.flowRate==0&&!noCapillary.capillaryAvailable&&noCapillary.capillaryRate==0,
        "C2 sigma0 invented a capillary constraint");
    const double dt=static_cast<double>(std::min(1e-5L,.02L/expectedRate));const Clock end{1,0,dt,Position::End};
    auto radial=config;radial.bodyAcceleration.x=.01;
    rejects([&]{HydroVofStep::prepare(*mesh,accepted,end,radial);},"C2 RZ radial body acceleration remains forbidden");
    auto wrongEnd=end;wrongEnd.begin=dt/2;
    rejects([&]{HydroVofStep::prepare(*mesh,accepted,wrongEnd,config);},"C2 RZ actual full-clock continuity");
    bool capillaryRetry=false;const Clock tooLarge{1,0,static_cast<double>(2*config.maximumCapillaryCourant/expectedRate),Position::End};
    try{HydroVofStep::prepare(*mesh,accepted,tooLarge,config);}
    catch(const HydroTwoPhaseStepRetry& error){capillaryRetry=error.reason()==HydroTwoPhaseRetryReason::Capillary;}
    check(capillaryRetry,"C2 RZ retained hard capillary Courant retry gate");
    check(accepted.sameIdentity(oldState),"C2 read/rejection changed accepted state");
    stage="C2 real closed RZ sigma-positive A+C step";
    std::cout<<"LOCAL_PLIC_RZ_CAPILLARY_STEP_PREPARE,n=32,subcells=8,dt="<<dt<<",mixed="<<mixed
        <<",axisMixed="<<axisMixed<<",h="<<h<<",capillaryRate="<<rates.capillaryRate<<std::endl;
    const auto step=HydroVofStep::prepare(*mesh,accepted,end,config);
    const auto& frozen=step.frozenMidpoint();const auto& dual=frozen.material();const auto& candidate=step.candidateState();
    const auto& transport=step.transport();const auto midClock=transport.midpoint().timeLayer();const auto midId=transport.midpoint().identity();
    check(step.hasCapillaryLoad()&&frozen.hasInterfaceLoad()&&step.report().linearSolves==2&&step.report().startupCorrection,
        "C2 actual startup must consume RZ capillary load in its original two-stage NS");
    const auto& capillary=step.capillaryLoad();capillary.validateVofFor(*mesh,frozen.pressureLoad(),midId,midClock);
    const auto* actualForce=std::get_if<HydroBalancedCsfLoad>(&frozen.forces().interfaceLoad);
    check(actualForce&&actualForce->vofIdentity()==midId&&sameClock(capillary.timeLayer(),midClock)
        &&&actualForce->momentumRhsLoad()==&capillary.momentumRhsLoad(),
        "C2 original frozen NS does not retain the SAME actual midpoint RZ HF force");
    check(frozen.interfaceLoads()==capillary.momentumRhsLoad(),"C2 original NS changed actual capillary RHS sign/measure");
    check(dual.material().vofTransport().sameInstance(transport)&&dual.sources().vofStamp&&dual.sources().vofStamp->sameInstance(transport),
        "C2 actual NS lost its unique shared VOF F owner");
    for(const auto* mass:{&dual.beginMass(),&dual.midpointMass(),&dual.endMass()})
        check(mass->representation()==HydroMomentumMass::DensityRepresentation::CellAverageProjection
            &&mass->liquidDensity()==config.liquidDensity&&mass->gasDensity()==config.gasDensity,"C2 original NS lost actual P0 phase-density masses");
    for(int outer=0;outer<frozen.viscosityMatrix().outerSize();++outer)
        for(HydroFrozenMidpoint::Matrix::InnerIterator entry(frozen.viscosityMatrix(),outer);entry;++entry)
            check(entry.value()==0,"C2 zero physical viscosity acquired a matrix floor");
    const auto& solved=frozen.solveResult().values();const auto cp=frozen.pressureLoad().apply(*mesh,solved);
    const Eigen::Map<const Eigen::VectorXd> u(solved.normalVelocity().data(),layout.velocityCount());
    const Eigen::VectorXd au=frozen.velocityMatrix()*u;long double residual=0,scale=0;
    for(int row:frozen.saddle().freeVelocityDofs()){
        residual=std::hypot(residual,static_cast<long double>(au[row])+cp[row]-frozen.integralMomentumRhs()[row]);
        scale=std::hypot(scale,std::abs(static_cast<long double>(au[row]))+std::abs(static_cast<long double>(cp[row]))+std::abs(frozen.integralMomentumRhs()[row]));}
    check(scale>0&&residual/scale<=1e-11L,"C2 actual original A+C free-momentum equation failed");
    for(const auto& stageReport:step.report().stages){
        check(stageReport.momentumPowerRelativeDefect<=config.momentumPowerTolerance&&stageReport.saddle.algebraicRelativeResidual<=1e-11
            &&stageReport.saddle.freeMomentumDefect<=1e-11,"C2 original same-owner momentum/power/saddle gates");
    }
    for(const auto& component:frozen.appliedConvection().report().components)
        check(component.momentumIdentity<=1e-12L&&component.kineticIdentity<=1e-12L,"C2 original closed KEP momentum/kinetic gates");
    check(transport.report().halfRelativeConservation<=1e-12L&&transport.report().fullRelativeConservation<=1e-12L,
        "C2 original shared-face phase inventory gates");
    std::array<long double,2> total{};
    for(const auto& cell:transport.endpoint().plicCells())for(int phase=0;phase<2;++phase)total[phase]+=cell.inventory[phase];
    for(int phase=0;phase<2;++phase)near(total[phase],oldTotal[phase],oldTotal[phase],"C2 closed actual independent phase inventory",1e-12L);
    auto wrongMid=midClock;++wrongMid.step;
    rejects([&]{capillary.validateVofFor(*mesh,frozen.pressureLoad(),midId,wrongMid);},"C2 capillary full-clock source gate");
    rejects([&]{capillary.validateVofFor(*mesh,frozen.pressureLoad(),sourceId,midClock);},"C2 capillary distinct original source owner gate");
    check(accepted.sameIdentity(oldState)&&source.identity()==sourceId&&source.plicCells().data()==sourceStorage,
        "C2 prepare changed accepted source identity/storage");
    for(std::size_t ci=0;ci<saved.size();++ci)check(source.plicCells()[ci].inventory==saved[ci].inventory,"C2 prepare changed original independent source stocks");
    check(candidate.retainedMass()&&candidate.retainedMass()->masses().data()==dual.endMass().masses().data()
        &&candidate.values().normalVelocity()==solved.normalVelocity()&&sameClock(candidate.timeLayer(),end),
        "C2 candidate lost the real NS endpoint values/mass/full clock");
    step.commit(*mesh,accepted);
    check(accepted.sameIdentity(candidate)&&accepted.velocityHistory()&&sameClock(accepted.timeLayer(),end),"C2 transaction did not adopt exact real candidate/history");
    rejects([&]{step.commit(*mesh,accepted);},"C2 stale repeated accepted transaction");
    std::cout<<"LOCAL_PLIC_RZ_CAPILLARY_STEP_PASS,dt="<<dt<<",capillaryCourant="<<step.report().capillaryCourant
        <<",linearSolves="<<step.report().linearSolves<<",A_C_relative="<<residual/scale
        <<",module_step_only_no_dynamic_accuracy_claim"<<std::endl;
}
// END CONDITIONAL RZ C2 REAL CAPILLARY STEP CONTRACT
void vofViscosityCapillaryOverlapCase(int initialError)
{
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
    stage="v71 same-accepted 32x32 RZ sphere nonzero viscosity/capillary overlap";
    struct Runtime {
        int threads=omp_get_max_threads(),dynamic=omp_get_dynamic(),error=errno;
        std::fenv_t native{};unsigned csr=_mm_getcsr();
        Runtime(){check(std::fegetenv(&native)==0,"overlap cannot capture caller FP state");}
        ~Runtime(){omp_set_num_threads(threads);omp_set_dynamic(dynamic);std::fesetenv(&native);_mm_setcsr(csr);errno=error;}
    } restore;
    omp_set_dynamic(0);check(omp_get_thread_limit()>=8,"overlap qualification needs actual eight-thread capacity");
    check(std::fesetround(FE_TONEAREST)==0,"overlap cannot establish nearest controls");_mm_setcsr(_mm_getcsr()&~0x8040u);
    unsigned nativeControl=0;check(_controlfp_s(&nativeControl,0,0)==0&&(nativeControl&_MCW_EM)==_MCW_EM
        &&(_mm_getcsr()&0x1f80u)==0x1f80u,"overlap requires masked native/SSE controls");
    DictionaryMesh dictionary;dictionary.Nx=dictionary.Ny=32;dictionary.maxCellNum=16384;dictionary.maxInitialIterations=0;
    auto mesh=std::make_unique<QuadtreeMesh>(dictionary);const NativeMeshIdentity live(*mesh);
    GeometricVofInitializationConfig initial;initial.subcellsPerAxis=8;
    const auto source=GeometricVofInitialization::prepare(*mesh,true,[](Vector2D p){
        return p.x*p.x+(p.y-.5)*(p.y-.5)-.2*.2;},initial);
    const auto sourceId=source.identity();const auto* sourceStorage=source.plicCells().data();const auto saved=source.plicCells();
    const auto metrics=vofMaterialMetrics(*mesh,source);const auto& layout=metrics.layout();
    check(layout.domain().axisymmetric()&&layout.pressureCount()==1024,"overlap lost the original C2 RZ fixture");
    const auto phase=HydroPhaseGeometry::build(*mesh,metrics,source);
    const HydroNativeValues zero(layout,std::vector<double>(layout.velocityCount()),std::vector<double>(layout.pressureCount()));
    const auto accepted=HydroTwoPhaseState::createVof(*mesh,phase,zero);const auto acceptedBefore=accepted;
    HydroVofAcConfig config;config.liquidDensity=1000;config.gasDensity=1.226;
    config.liquidViscosity=.01;config.gasViscosity=.001;config.surfaceTension=.7;
    config.momentumDensityRepresentation=HydroMomentumMass::DensityRepresentation::CellAverageProjection;
    const auto rates=HydroVofStep::stabilityRates(*mesh,accepted,config);
    check(rates.capillaryAvailable&&rates.capillaryRate>0,"overlap has no actual positive-sigma interface");
    const double dt=static_cast<double>(std::min(1e-5L,.02L/rates.capillaryRate));const Clock end{1,0,dt,Position::End};
    const auto reset=[&]{check(std::feclearexcept(FE_ALL_EXCEPT)==0,"overlap cannot clear entry status");
        _mm_setcsr(_mm_getcsr()&~0x3fu);check(std::feraiseexcept(FE_INEXACT)==0,"overlap cannot prime actual sticky state");_mm_setcsr(_mm_getcsr()|0x20u);errno=initialError;};
    omp_set_num_threads(1);reset();const auto serial=HydroVofStep::prepare(*mesh,accepted,end,config);
    const int serialErrno=errno;const unsigned serialCsr=_mm_getcsr();
    omp_set_num_threads(8);reset();const auto concurrent=HydroVofStep::prepare(*mesh,accepted,end,config);
    const int concurrentErrno=errno;const unsigned concurrentCsr=_mm_getcsr();
    check(std::fegetround()==FE_TONEAREST&&(serialCsr&~0x3fu)==(concurrentCsr&~0x3fu)
        &&serialErrno==concurrentErrno&&(serialCsr&0x20u)&&(concurrentCsr&0x20u),
        "overlap changed controls/errno or lost preexisting sticky inexact");
    check(serial.beginState().sameIdentity(accepted)&&concurrent.beginState().sameIdentity(accepted)
        &&serial.report().linearSolves==concurrent.report().linearSolves&&serial.report().linearSolves==2
        &&serial.report().startupCorrection&&concurrent.report().startupCorrection,
        "overlap changed accepted owner or original actual AC solve count");
    check(serial.report().stages.size()==concurrent.report().stages.size(),"overlap changed actual stage count");
    for(const auto* step:{&serial,&concurrent}){
        const bool expected=step==&concurrent;
        for(const auto& report:step->report().stages){
            check(report.timing.concurrentViscosityCapillary==expected,"overlap actual producer receipt did not select requested route");
            check(report.momentumPowerRelativeDefect<=config.momentumPowerTolerance&&report.saddle.algebraicRelativeResidual<=1e-11
                &&report.saddle.freeMomentumDefect<=1e-11,"overlap failed original physical momentum/power/saddle gates");
        }
        const auto& frozen=step->frozenMidpoint();const auto& transport=step->transport();const auto& load=step->capillaryLoad();
        check(step->hasCapillaryLoad()&&frozen.hasInterfaceLoad()&&frozen.viscosityMatrix().norm()>0
            &&load.curvatureReport().heightFunctionCells+load.curvatureReport().heightFunctionFitCells>0
            &&load.curvatureReport().projectedFaces>0,"overlap did not prepare actual nonzero HF and viscosity");
        load.validateVofFor(*mesh,frozen.pressureLoad(),transport.midpoint().identity(),transport.midpoint().timeLayer());
        check(frozen.interfaceLoads()==load.momentumRhsLoad()&&sameClock(step->candidateState().timeLayer(),end)
            &&sameClock(transport.endTime(),end),"overlap lost actual same-midpoint HF owner or complete End clock");
        check(transport.report().halfRelativeConservation<=1e-12L&&transport.report().fullRelativeConservation<=1e-12L,
            "overlap failed original independent shared-face phase conservation");
        for(const auto& component:frozen.appliedConvection().report().components)
            check(component.momentumIdentity<=1e-12L&&component.kineticIdentity<=1e-12L,"overlap failed original closed KEP gates");
    }
    long double peak=0;std::size_t compared=0;
    const auto residual=[&](long double error,long double scale,const char* label){
        check(std::isfinite(error)&&std::isfinite(scale)&&(scale>0?error<=1e-11L*scale:error==0),label);
        if(scale>0)peak=std::max(peak,error/scale);++compared;};
    const auto values=[&](const auto& a,const auto& b,const char* label,bool unavailable=false){
        check(a.size()==b.size(),"overlap changed vector layout");long double error=0,scaleA=0,scaleB=0;
        for(std::size_t i=0;i<a.size();++i){if(unavailable&&(!std::isfinite(a[i])||!std::isfinite(b[i]))){
                check(std::isnan(a[i])&&std::isnan(b[i]),"overlap changed curvature availability");continue;}
            error=std::hypot(error,static_cast<long double>(a[i])-b[i]);scaleA=std::hypot(scaleA,static_cast<long double>(a[i]));scaleB=std::hypot(scaleB,static_cast<long double>(b[i]));}
        residual(error,std::max(scaleA,scaleB),label);};
    const auto matrix=[&](const auto& a,const auto& b,const char* label){
        check(a.rows()==b.rows()&&a.cols()==b.cols(),"overlap changed matrix layout");
        const HydroFrozenMidpoint::Matrix delta=a-b;residual(delta.norm(),std::max(a.norm(),b.norm()),label);};
    values(serial.candidateState().values().normalVelocity(),concurrent.candidateState().values().normalVelocity(),"overlap changed accepted native velocity");
    values(serial.candidateState().values().pressure(),concurrent.candidateState().values().pressure(),"overlap changed accepted physical pressure");
    const auto& a=serial.frozenMidpoint();const auto& b=concurrent.frozenMidpoint();
    matrix(a.viscosityMatrix(),b.viscosityMatrix(),"overlap changed actual viscosity matrix");
    matrix(a.velocityMatrix(),b.velocityMatrix(),"overlap changed actual Frozen velocity matrix");
    values(a.integralMomentumRhs(),b.integralMomentumRhs(),"overlap changed actual Frozen RHS");
    values(serial.capillaryLoad().momentumRhsLoad(),concurrent.capillaryLoad().momentumRhsLoad(),"overlap changed actual HF load");
    values(serial.capillaryLoad().heightFunctionCellEffectiveCurvatures(),concurrent.capillaryLoad().heightFunctionCellEffectiveCurvatures(),
        "overlap changed actual HF curvature",true);
    check(serial.capillaryLoad().heightFunctionCellCurvatureMethods()==concurrent.capillaryLoad().heightFunctionCellCurvatureMethods(),
        "overlap changed actual HF estimator methods");
    for(unsigned phaseIndex=0;phaseIndex<2;++phaseIndex){long double error=0,scaleA=0,scaleB=0;
        const auto& p=serial.candidateState().phaseGeometry().vofPublication().plicCells();
        const auto& q=concurrent.candidateState().phaseGeometry().vofPublication().plicCells();
        check(p.size()==q.size(),"overlap changed candidate stock layout");
        for(std::size_t i=0;i<p.size();++i){error+=std::abs(p[i].inventory[phaseIndex]-q[i].inventory[phaseIndex]);scaleA+=std::abs(p[i].inventory[phaseIndex]);scaleB+=std::abs(q[i].inventory[phaseIndex]);}
        residual(error,std::max(scaleA,scaleB),"overlap changed independent phase inventory");}
    const auto& p=serial.transport().faces();const auto& q=concurrent.transport().faces();check(p.size()==q.size(),"overlap changed native face count");
    for(auto member:{&GeometricVofIntegratedFlux::total,&GeometricVofIntegratedFlux::positiveTotal,&GeometricVofIntegratedFlux::negativeTotal,
            &GeometricVofIntegratedFlux::rawLiquid,&GeometricVofIntegratedFlux::liquid,&GeometricVofIntegratedFlux::gas,&GeometricVofIntegratedFlux::mass,
            &GeometricVofIntegratedFlux::positiveLiquid,&GeometricVofIntegratedFlux::negativeLiquid,&GeometricVofIntegratedFlux::positiveGas,
            &GeometricVofIntegratedFlux::negativeGas,&GeometricVofIntegratedFlux::limiter,&GeometricVofIntegratedFlux::exhaustionLiquidCorrection,
            &GeometricVofIntegratedFlux::exhaustionGasCorrection}){
        long double error=0,scaleA=0,scaleB=0;for(std::size_t i=0;i<p.size();++i){
            check(p[i].nativeRow==q[i].nativeRow&&p[i].negativeCell==q[i].negativeCell&&p[i].positiveCell==q[i].positiveCell,"overlap changed shared native face incidence");
            for(bool half:{true,false}){const auto x=(half?p[i].half:p[i].full).*member,y=(half?q[i].half:q[i].full).*member;
                error=std::hypot(error,x-y);scaleA=std::hypot(scaleA,x);scaleB=std::hypot(scaleB,y);}}
        residual(error,std::max(scaleA,scaleB),"overlap changed actual independent half/full flux");}
    check(accepted.sameIdentity(acceptedBefore)&&source.identity()==sourceId&&source.plicCells().data()==sourceStorage
        &&accepted.values().normalVelocity()==zero.normalVelocity()&&accepted.values().pressure()==zero.pressure(),"overlap mutated actual accepted/source state");
    for(std::size_t i=0;i<saved.size();++i)check(source.plicCells()[i].inventory==saved[i].inventory,"overlap mutated authoritative source inventory");live.validate(*mesh);
    std::cout<<std::setprecision(17)<<"VOF_VISCOSITY_CAPILLARY_OVERLAP_PASS,scenes,1,source_shared,1,serial_threads,1,candidate_threads,8"
        <<",entry_errno,"<<initialError<<",rho_l,1000,rho_g,1.226,mu_l,.01,mu_g,.001,sigma,.7,dt,"<<dt<<",compared_quantities,"<<compared
        <<",peak_relative,"<<peak<<",actual_route_receipts_checked,1,original_acceptance_gates_checked,1,no_speed_claim\n";
#else
    throw std::runtime_error("overlap qualification requires actual MSVC x64 OpenMP runtime");
#endif
}

void vofViscosityCapillaryOverlapOnly()
{vofViscosityCapillaryOverlapCase(0);vofViscosityCapillaryOverlapCase(ENOENT);}

// F11 actual native-CDF capacity witness, no injected/root defect.
void localPlicProjectedCapacityOnly()
{
    stage="F11 strict pure parent -> actual L4 integer-CDF -> P0 material/dual";
    if constexpr(std::numeric_limits<long double>::digits!=53){
        std::cout<<"LOCAL_PLIC_PROJECTED_CAPACITY_NOT_APPLICABLE,native_digits="
            <<std::numeric_limits<long double>::digits<<",actual_binary64_CDF_fixture"<<std::endl;return;
    }
    DictionaryMesh input;input.Nx=input.Ny=1;input.minCoord={0,9};input.maxCoord={1,10};
    input.maxCellNum=2048;input.maxInitialIterations=0;input.maxLevel_local={4};
    QuadtreeMesh coarse(input);auto fine=coarse.clone();
    DictionaryMesh::RefineBox box;box.level=4;box.min=input.minCoord;box.max=input.maxCoord;fine->applyRefineBoxes({box});
    check(coarse.cellIndexNum()==1&&fine->cellIndexNum()==256,"F11 real uniform L4 subdivision missing");
    const auto coarseDomain=domain(coarse,true),fineDomain=domain(*fine,true);
    const Clock oldClock{7,.125,.25,Position::End},next{8,.25,.25+1e-5,Position::End};
    constexpr long double parentC=0x1.921fb54442d18p+1L,childC=0x1.921fb54442d18p-11L;
    // Independent small integer RNE oracle: target is the ninth axial child
    // in the first radial ring; physical prefix weights are 8/4096,9/4096.
    // No production CDF/exact helper computes this expected inventory.
    const auto rnEven=[](std::uint64_t n,std::uint64_t d){const auto q=n/d,r=n%d;
        return q+std::uint64_t(2*r>d||(2*r==d&&(q&1)));};
    constexpr std::uint64_t units=7074237752028440ULL;
    const auto count=rnEven(units*9,4096)-rnEven(units*8,4096);
    const long double expected=static_cast<long double>(count)*std::ldexp(1.L,-51);
    check(expected==0x1.921fb54442000p-11L&&expected-childC==-0x1.a3p-52L,
        "F11 independent integer oracle lost actual failed native bits");
    check(coarseDomain.cellVolumes()[0]==parentC,"F11 exact strict native parent capacity changed");
    unsigned massChecks=0,dualChecks=0,phaseCases=0;
    for(unsigned phase:{1u,0u}){ // actual pure-G failure first, then symmetric pure-L
        std::vector<std::array<long double,2>> initial(1);initial[0][phase]=parentC;
        const auto root=Publication::prepareLocalPlicInventory(coarse,coarseDomain,oldClock,initial,{});
        const auto rootId=root.identity();const auto savedRoot=root.plicCells()[0].inventory;
        check(root.capacityDefects()[0].compare(LAMP::Geometry::Detail::FiniteVolume::Dyadic{})==0,
            "F11 strict root already contains a defect");
        const auto transfer=GeometricVofMeshTransfer::prepare(coarse,root,*fine,fineDomain,oldClock);
        const auto source=transfer.candidate();source.validateFor(*fine,source.identity(),oldClock);
        int named=-1;bool hasOvercapacity=false;long double sum=0;
        for(int ci=0;ci<fine->cellIndexNum();++ci){const auto& c=source.plicCells()[ci];const auto& r=c.rectangle;
            check(c.inventory[1-phase]==0,"F11 real refine created the other phase");sum+=c.inventory[phase];
            hasOvercapacity|=c.inventory[phase]>source.cells()[ci].capacity;
            if(r.lower.x==0&&r.upper.x==1./16&&r.lower.y==9.5&&r.upper.y==9.5625)named=ci;}
        check(named>=0&&hasOvercapacity,"F11 real CDF lacks named child or genuine quantization excess");
        const auto& namedCell=source.plicCells()[named];
        check(namedCell.inventory[phase]==expected&&source.cells()[named].capacity==childC,
            "F11 actual transfer did not reproduce observed Q/C bits");
        check(source.capacityDefects()[named].compare(LAMP::Geometry::Detail::FiniteVolume::Dyadic(-0x1.a3p-52L))==0,
            "F11 actual typed certificate differs from independent integer oracle");
        near(sum,parentC,parentC,"F11 true native-CDF phase conservation",1e-12L);
        check(std::abs(expected-childC)/childC>512.L*std::numeric_limits<double>::epsilon()
            &&std::abs(expected-childC)/childC<1e-12L,"F11 witness does not separate material/geometry budgets");
        std::vector<std::array<long double,2>> unproven(source.plicCells().size());
        for(std::size_t ci=0;ci<unproven.size();++ci)unproven[ci]=source.plicCells()[ci].inventory;
        rejects([&]{Publication::prepareLocalPlicInventory(*fine,fineDomain,oldClock,unproven,{});},
            "F11 overcapacity child cannot bypass producer lineage through raw root");
        const auto saved=source.plicCells();const auto id=source.identity();const auto* address=source.plicCells().data();
        std::cout<<std::setprecision(std::numeric_limits<long double>::max_digits10)
            <<"F11_REAL_CDF_BEFORE_MASS,phase="<<phase<<",parent="<<parentC<<",named="<<named
            <<",stock="<<namedCell.inventory[phase]<<",capacity="<<source.cells()[named].capacity
            <<",kappa="<<expected-childC<<",relative="<<std::abs(expected-childC)/childC
            <<",actual_refine_leaves="<<fine->cellIndexNum()<<std::endl;
        const std::vector<double> zero(fineDomain.nativeFaceCount(),0);
        const auto slab=Slab::prepare(*fine,source,next,zero,zero);
        const auto transport=Transport::prepareLocalPlicSweep(*fine,source,slab,1000,1.226L);
        const auto geometry=HydroMomentumGeometry::build(*fine,vofMaterialMetrics(*fine,source));
        const auto before=HydroPhaseGeometry::build(*fine,geometry.metrics(),source);
        const auto middle=HydroPhaseGeometry::build(*fine,geometry.metrics(),transport.midpoint());
        const auto after=HydroPhaseGeometry::build(*fine,geometry.metrics(),transport.endpoint());
        const auto material=HydroMaterialTransport::prepareGeometricVof(*fine,geometry,before,middle,after,transport);
        const auto m0=HydroMomentumMass::buildCellAverageProjection(*fine,geometry,before,1000,1.226L);
        const auto mm=HydroMomentumMass::buildCellAverageProjection(*fine,geometry,middle,1000,1.226L);
        const auto m1=HydroMomentumMass::buildCellAverageProjection(*fine,geometry,after,1000,1.226L);
        check(m0.masses()==mm.masses()&&m0.masses()==m1.masses(),"F11 zero-flow changed actual retained M");
        for(const auto* m:{&m0,&mm,&m1}){
            const auto& represented=m->representedOverlapPhaseVolumes();
            for(int u=0;u<m->layout().velocityCount();++u){const int component=m->layout().normalAxis(u)==HydroDofLayout::NormalAxis::First?0:1;
                const auto& pieces=geometry.volumes()[u].overlaps;
                for(std::size_t k=0;k<pieces.size();++k){const auto& p=pieces[k];const auto& q=source.plicCells()[m->layout().nativeCell(p.pressure)].inventory;
                    const auto actual=represented[u][k];check(actual.liquid==p.cellFraction*q[0]&&actual.gas==p.cellFraction*q[1],
                        "F11 capacity qualification changed ordinary independent P0 material bits");
                    check(!m->usesFiniteP0Partition(p.pressure,component,0)&&!m->usesFiniteP0Partition(p.pressure,component,1),
                        "F11 ordinary witness silently entered finite-CDF material branch");++massChecks;}}
            for(double mass:m->masses())check(std::isfinite(mass)&&mass>0,"F11 actual total M is not positive finite");
        }
        const auto dual=HydroDualMaterial::prepare(*fine,material,m0,mm,m1);
        for(const auto& c:dual.cells()){check(c.massChange==0&&c.outwardIntegral==0&&c.integratedDefect==0&&c.defectRate==0,
            "F11 zero-flow actual stocks/F acquired canonical DeltaM/R");++dualChecks;}
        for(std::size_t ci=0;ci<saved.size();++ci)for(unsigned p=0;p<2;++p){
            check(source.plicCells()[ci].inventory[p]==saved[ci].inventory[p]
                &&transport.midpoint().plicCells()[ci].inventory[p]==saved[ci].inventory[p]
                &&transport.endpoint().plicCells()[ci].inventory[p]==saved[ci].inventory[p],
                "F11 zero flow or P0 build changed a raw retained inventory");}
        check(source.identity()==id&&source.plicCells().data()==address&&root.identity()==rootId
            &&root.plicCells()[0].inventory==savedRoot,"F11 operation changed source owner/storage/strict root");
        auto wrong=oldClock;wrong.position=Position::Midpoint;
        rejects([&]{m0.validateVofFor(*fine,source.identity(),wrong);},"F11 fullclock mismatch");
        ++phaseCases;
    }
    check(phaseCases==2&&massChecks>0&&dualChecks>0,"F11 actual material/dual phases were not consumed");
    std::cout<<"LOCAL_PLIC_PROJECTED_CAPACITY_PASS,pure_phases="<<phaseCases<<",ordinary_overlap_checks="<<massChecks
        <<",canonical_dual_checks="<<dualChecks<<",rhoL=1000,rhoG=1.226,unchanged_512eps_and_ownphase_gates"<<std::endl;
}

// D17 actual D16 RZ rectangle; no new flow/parameter matrix.
void localPlicCapacityCarryOnly()
{
    stage="D17 true RZ pure-phase coarsen -> zero-flow independent inventory";
    if constexpr(std::numeric_limits<long double>::digits!=53){
        std::cout<<"LOCAL_PLIC_CAPACITY_CARRY_NOT_APPLICABLE,native_digits="
            <<std::numeric_limits<long double>::digits<<",actual_D16_binary64_fixture"<<std::endl;
        return;
    }
    DictionaryMesh dictionary;dictionary.Nx=dictionary.Ny=1;dictionary.maxCellNum=64;
    dictionary.maxInitialIterations=0;dictionary.maxLevel_local={1};
    dictionary.minCoord={11./64,5./16};dictionary.maxCoord={3./16,21./64};
    QuadtreeMesh coarse(dictionary);auto fine=coarse.clone();
    DictionaryMesh::RefineBox box;box.level=1;box.min=dictionary.minCoord;box.max=dictionary.maxCoord;
    fine->applyRefineBoxes({box});
    check(coarse.cellIndexNum()==1&&fine->cellIndexNum()==4,"D17 original rectangle partition changed");
    const auto coarseDomain=domain(coarse,true),fineDomain=domain(*fine,true);
    const long double nativeC=0x1.2106ca4910069p-12L,expectedQ=0x1.2106ca491006ap-12L;
    check(coarseDomain.cellVolumes()[0]==nativeC,"D17 actual native parent C bits changed");
    const Clock oldClock{7,.125,.25,Position::End},next{8,.25,.25+1e-5,Position::End};
    unsigned carryChecks=0,massChecks=0,dualChecks=0;
    const auto sameNative=[](long double a,long double b){return a==b&&std::signbit(a)==std::signbit(b);};
    for(unsigned phase=0;phase<2;++phase){
        std::vector<std::array<long double,2>> initial(fine->cellIndexNum());
        for(int ci=0;ci<fine->cellIndexNum();++ci)initial[ci][phase]=capacity(bounds(*fine,ci),true);
        const auto root=Publication::prepareLocalPlicInventory(*fine,fineDomain,oldClock,initial,{});
        const auto transfer=GeometricVofMeshTransfer::prepare(*fine,root,coarse,coarseDomain,oldClock);
        const auto source=transfer.candidate();source.validateFor(coarse,source.identity(),oldClock);
        check(source.plicCells()[0].inventory[phase]==expectedQ&&source.plicCells()[0].inventory[1-phase]==0,
            "D17 coarsen did not retain the actual independent donor phase sum");
        check(expectedQ>nativeC&&expectedQ-nativeC==0x1p-64L,"D17 one-native-ULP premise changed");
        check(source.alpha()[0]==(phase==0?1.:0.),"D17 pure diagnostic alpha must retain actual phase ownership");
        const std::vector<std::array<long double,2>> unproven{source.plicCells()[0].inventory};
        rejects([&]{Publication::prepareLocalPlicInventory(coarse,coarseDomain,oldClock,unproven,{});},
            "D17 arbitrary raw overcapacity stock cannot self-certify");
        const auto sourceId=source.identity();const auto* sourceAddress=source.plicCells().data();
        const auto saved=source.plicCells();
        const std::vector<double> zero(coarseDomain.nativeFaceCount(),0);
        const auto slab=Slab::prepare(coarse,source,next,zero,zero);
        const auto transport=Transport::prepareLocalPlicSweep(coarse,source,slab,1000,1.226L);
        auto midClock=next;midClock.position=Position::Midpoint;
        check(sameClock(transport.midpoint().timeLayer(),midClock)&&sameClock(transport.endpoint().timeLayer(),next),
            "D17 carried half/end actual clocks differ");
        for(const auto* extent:{&transport.midpoint(),&transport.endpoint()}){
            extent->validateFor(coarse,extent->identity(),extent->timeLayer());
            for(unsigned p=0;p<2;++p){
                check(sameNative(extent->plicCells()[0].inventory[p],saved[0].inventory[p]),
                    "D17 zero flow changed raw independent phase bits");++carryChecks;
            }
            check(extent->plicCells()[0].inventory[phase]>extent->cells()[0].capacity,
                "D17 zero flow erased the genuine finite capacity defect");
        }
        for(const auto& face:transport.faces())for(const auto* extent:{&face.half,&face.full})
            check(extent->total==0&&extent->liquid==0&&extent->gas==0&&extent->mass==0,
                "D17 zero flow acquired an exported phase/material flux");
        const auto geometry=HydroMomentumGeometry::build(coarse,vofMaterialMetrics(coarse,source));
        const auto before=HydroPhaseGeometry::build(coarse,geometry.metrics(),source);
        const auto middle=HydroPhaseGeometry::build(coarse,geometry.metrics(),transport.midpoint());
        const auto after=HydroPhaseGeometry::build(coarse,geometry.metrics(),transport.endpoint());
        const auto material=HydroMaterialTransport::prepareGeometricVof(coarse,geometry,before,middle,after,transport);
        check(material.vofTransport().sameInstance(transport)&&material.liquidDensity()==1000&&material.gasDensity()==1.226L,
            "D17 material lost actual flux owner or physical densities");
        for(const auto& cell:material.cells()){
            check(cell.massChange==0&&cell.outwardIntegral==0&&cell.integratedDefect==0,
                "D17 identical independent stocks acquired primal DeltaM or R");++massChecks;
        }
        const auto mass0=HydroMomentumMass::buildCellAverageProjection(coarse,geometry,before,1000,1.226L);
        const auto massM=HydroMomentumMass::buildCellAverageProjection(coarse,geometry,middle,1000,1.226L);
        const auto mass1=HydroMomentumMass::buildCellAverageProjection(coarse,geometry,after,1000,1.226L);
        check(mass0.masses()==massM.masses()&&mass0.masses()==mass1.masses(),
            "D17 unchanged independent material inventory changed P0 M");
        const auto dual=HydroDualMaterial::prepare(coarse,material,mass0,massM,mass1);
        check(dual.material().vofTransport().sameInstance(transport),"D17 dual replaced the accepted phase flux owner");
        for(const auto& cell:dual.cells()){
            check(cell.massChange==0&&cell.outwardIntegral==0&&cell.integratedDefect==0&&cell.defectRate==0,
                "D17 identical raw stocks and zero F acquired canonical dual DeltaM/R");++dualChecks;
        }
        check(source.identity()==sourceId&&source.plicCells().data()==sourceAddress&&sameClock(source.timeLayer(),oldClock),
            "D17 candidate changed original carried owner/storage/clock");
        for(unsigned p=0;p<2;++p)check(sameNative(source.plicCells()[0].inventory[p],saved[0].inventory[p]),
            "D17 prepare changed accepted carried inventory");
    }
    check(carryChecks==8&&massChecks>0&&dualChecks>0,"D17 actual carry/material/dual consumer coverage missing");
    std::cout<<std::setprecision(17)<<"LOCAL_PLIC_CAPACITY_CARRY_PASS,pure_phases=2,raw_checks="<<carryChecks
        <<",primal_checks="<<massChecks<<",dual_checks="<<dualChecks<<",native_C="<<nativeC
        <<",retained_Q="<<expectedQ<<",excess="<<expectedQ-nativeC<<",rhoL=1000,rhoG=1.226,zero_flow_only"<<std::endl;
}

void localSmoothOnly()
{
    stage="LocalSmooth explicit same-source split/WY";
    unsigned nonzeroFaces[2]{},positiveRemainders=0,changedStocks=0,sourceChecks=0,restrictionChecks=0;
    for(unsigned scene=0;scene<2;++scene){
        DictionaryMesh specification;specification.Nx=specification.Ny=scene?8:16;
        specification.maxCellNum=4096;specification.maxInitialIterations=0;QuadtreeMesh mesh(specification);
        const auto phi=[scene](Vector2D p){return scene?p.x-.43:std::hypot(p.x-.487,p.y-.513)-.18;};
        const auto source=GeometricVofInitialization::prepare(mesh,false,phi);
        const auto saved=source.plicCells();const auto sourceId=source.identity();const auto sourceClock=source.timeLayer();
        const auto* sourceStorage=source.plicCells().data();
        VofSmoothGeometryConfig smoothConfig;smoothConfig.allowBoundaryEndpoints=scene!=0;
        const auto smooth=VofSmoothGeometry::prepare(mesh,source,smoothConfig);
        std::vector<double> velocity(static_cast<std::size_t>(source.domain().nativeFaceCount()));
        if(scene){
            const auto psi=[](Point p){return (1.L/256)*p.x*(1-p.x)*p.y*(1-p.y);};
            for(const auto face:source.domain().faces(mesh)){const auto e=faceEndpoints(mesh,face);
                const auto integral=(face.horizontal?-1.L:1.L)*(psi(e[1])-psi(e[0]));
                velocity[static_cast<std::size_t>(face.nativeRow)]=static_cast<double>(integral/face.area);
                if(face.negativeCell<0||face.positiveCell<0)check(velocity[static_cast<std::size_t>(face.nativeRow)]==0,
                    "declared streamfunction is not an exactly closed native fixture");
            }
        }
        const Clock end{1,0,1./1024,Position::End};const auto slab=Slab::prepare(mesh,source,end,velocity,velocity);
        const GeometricVofTransportConfig config;
        const auto moved=Transport::prepareLocalSmoothSweep(mesh,source,smooth,slab,100,1,config,smoothConfig);
        // Retain the original physical ledgers/gates. Axis=-1 selects no old
        // axial PLIC oracle; that oracle is unrelated to this smooth method.
        l1Inspect(mesh,source,slab,moved,-1);
        check(moved.smoothSource()&&moved.smoothSource()->sameInstance(smooth)
            &&moved.smoothMidpoint()&&moved.smoothEndpoint(),"smooth transport omitted its actual immutable geometry owners");
        moved.smoothSource()->validateFor(mesh,source);
        for(bool half:{true,false}){
            const auto& target=half?moved.midpoint():moved.endpoint();
            const auto* geometry=half?moved.smoothMidpoint():moved.smoothEndpoint();
            geometry->validateFor(mesh,target);
            check(geometry->sourceIdentity()==target.identity()&&sameClock(geometry->sourceTime(),target.timeLayer())
                &&geometry->source().plicCells().data()==target.plicCells().data(),"smooth target is not the actual half/full publication instance");
            check(geometry->report().maximumOwnPhaseRelativeResidual<=smoothConfig.relativeVolumeTolerance,
                "smooth target failed the original own-phase geometry gate");
            for(std::size_t ci=0;ci<saved.size();++ci){
                const auto& old=saved[ci].inventory;const auto& now=target.plicCells()[ci].inventory;
                if(!scene)check(now==old,"zero smooth transport changed an independent phase stock");
                else if(old[0]>0&&old[1]>0){
                    check(now[0]>0&&now[1]>0,"partial smooth sweep lost a real positive phase remainder");
                    ++positiveRemainders;if(now!=old)++changedStocks;
                }
            }
        }
        if(!scene){
            const auto original=Transport::prepareLocalPlicSweep(mesh,source,slab,100,1,config);
            check(!original.smoothSource()&&!original.smoothMidpoint()&&!original.smoothEndpoint(),"original A path acquired a B1 geometry owner");
            for(std::size_t ci=0;ci<saved.size();++ci)check(original.midpoint().plicCells()[ci].inventory==moved.midpoint().plicCells()[ci].inventory
                &&original.endpoint().plicCells()[ci].inventory==moved.endpoint().plicCells()[ci].inventory,"zero smooth transport changed original A stocks");
            for(std::size_t row=0;row<moved.faces().size();++row){
                v4b1Exact(moved.faces()[row].half,original.faces()[row].half,"zero smooth/A cached half identity");
                v4b1Exact(moved.faces()[row].full,original.faces()[row].full,"zero smooth/A cached full identity");
                v4b1Zero(moved.faces()[row].half,"zero smooth half phase/material flux");
                v4b1Zero(moved.faces()[row].full,"zero smooth full phase/material flux");
            }
            const auto foreign=GeometricVofInitialization::prepare(mesh,false,phi);
            const auto foreignSmooth=VofSmoothGeometry::prepare(mesh,foreign,smoothConfig);
            rejects([&]{Transport::prepareLocalSmoothSweep(mesh,source,foreignSmooth,slab,100,1,config,smoothConfig);},
                "smooth sweep accepted a foreign same-clock source geometry");
        }
        for(int row=0;row<source.domain().nativeFaceCount();++row){
            const auto native=source.domain().face(mesh,row);const auto e=faceEndpoints(mesh,native);
            const Point mid{(e[0].x+e[1].x)/2,(e[0].y+e[1].y)/2};
            const auto parts=moved.restrictNativeFaces(mesh,{{row,e[0],e[1]},{row,e[1],e[0]},{row,e[0],mid},{row,mid,e[1]},{row,e[0],e[0]}});
            for(bool half:{true,false}){
                const auto& flux=half?moved.faces()[row].half:moved.faces()[row].full;
                v4b1Exact(half?parts[0].half:parts[0].full,flux,"smooth cached whole-face identity");
                v4b1Exact(half?parts[1].half:parts[1].full,flux,"smooth cached reversed-face identity");
                v4b1Zero(half?parts[4].half:parts[4].full,"smooth empty P0 restriction");
                const auto parent=v4b1Fields(flux),a=v4b1Fields(half?parts[2].half:parts[2].full),
                    b=v4b1Fields(half?parts[3].half:parts[3].full);
                for(std::size_t k=0;k+1<parent.size();++k){
                    if(parent[k]==0)check(a[k]==0&&b[k]==0,"zero smooth parent acquired a P0 subface value");
                    else {near(a[k],parent[k]/2,std::abs(parent[k]),"smooth retained numerical P0 half measure",1e-12L);
                        near(a[k]+b[k],parent[k],std::abs(parent[k]),"smooth P0 partition addback",1e-12L);}
                    ++restrictionChecks;
                }
                if(scene&&flux.total!=0)++nonzeroFaces[native.horizontal?1:0];
            }
        }
        check(source.identity()==sourceId&&sameClock(source.timeLayer(),sourceClock)&&source.plicCells().data()==sourceStorage,
            "smooth preparation/restriction/rejection replaced accepted source owner/clock/storage");
        for(std::size_t ci=0;ci<saved.size();++ci){
            const auto& old=saved[ci];const auto& now=source.plicCells()[ci];
            check(old.inventory==now.inventory&&old.measured==now.measured&&old.normal.x==now.normal.x
                &&old.normal.y==now.normal.y&&old.offset==now.offset&&old.resolved==now.resolved,"smooth method changed accepted phase/geometry bytes");
            ++sourceChecks;
        }
    }
    check(nonzeroFaces[0]>0&&nonzeroFaces[1]>0&&positiveRemainders>0&&changedStocks>0&&restrictionChecks>0,
        "declared smooth fixture missed real two-direction partial transport");
    std::cout<<"LOCAL_SMOOTH_PASS,scenes=2,zero_closed_circle=1,nonzero_plane=1,rhoL=100,rhoG=1,dt=0.0009765625"
        <<",x_nonzero="<<nonzeroFaces[0]<<",y_nonzero="<<nonzeroFaces[1]<<",positive_remainders="<<positiveRemainders
        <<",changed_stocks="<<changedStocks<<",source_checks="<<sourceChecks<<",P0_checks="<<restrictionChecks
        <<",original_gates_unchanged,charged_motion_unqualified\n";
}
void localPlicSingleOutgoingOnly()
{
    stage="LocalPlic exact single outgoing phase constraint";
    DictionaryMesh input;input.Nx=input.Ny=32;input.maxCellNum=4096;input.maxInitialIterations=0;
    QuadtreeMesh mesh(input);
    const auto taper=[](long double x){if(x<=0||x>=1)return 0.L;if(x>=.2L&&x<=.8L)return 1.L;
        const auto z=x<.2L?x/.2L:(1-x)/.2L;return z*z*(3-2*z);};
    const auto psi=[&](long double x,long double y){return .04L*y*taper(x)*taper(y);};
    const auto advance=[&](const Publication& source,double dt){
        std::vector<double> velocity(static_cast<std::size_t>(source.domain().nativeFaceCount()));
        for(const auto face:source.domain().faces(mesh)){const long double x=face.center.x,y=face.center.y,h=face.area/2;
            velocity[static_cast<std::size_t>(face.nativeRow)]=static_cast<double>(face.horizontal?
                -(psi(x+h,y)-psi(x-h,y))/face.area:(psi(x,y+h)-psi(x,y-h))/face.area);}
        const auto slab=Slab::prepare(mesh,source,{1,0,dt,Position::End},velocity,velocity);
        return Transport::prepareLocalPlicSweep(mesh,source,slab,100,1);
    };
    const auto rightFace=[&](const Publication& source,int cell){
        for(const auto face:source.domain().faces(mesh))if(!face.horizontal&&face.negativeCell==cell)return face.nativeRow;
        throw std::runtime_error("single outgoing fixture has no actual right face");
    };
    const auto retained=[&](const Publication& source,const auto& snapshot,GeometricVofIdentity identity,Clock clock){
        check(source.identity()==identity&&sameClock(source.timeLayer(),clock),"single outgoing preparation changed accepted owner/clock");
        for(std::size_t ci=0;ci<snapshot.size();++ci)check(source.plicCells()[ci].inventory==snapshot[ci].inventory,
            "single outgoing preparation changed independent accepted phase stock");
    };
    // The original D0 SDF initializer, not the distinct squared-circle App.
    const auto source=GeometricVofInitialization::prepare(mesh,false,[](Vector2D p){return std::hypot(p.x-.487,p.y-.513)-.18;});
    const auto saved=source.plicCells();const auto identity=source.identity();const auto clock=source.timeLayer();
    constexpr int donor=688;const auto stock=saved[donor].inventory[1];
    check(stock>0&&saved[donor].resolved,"D0 source lost its genuine resolved minority");
    const auto moved=advance(source,.025);const auto row=rightFace(source,donor);
    const auto& shared=moved.faces().at(static_cast<std::size_t>(row));const int receiver=shared.positiveCell;
    check(receiver>=0&&shared.negativeCell==donor,"D0 exhaust face is not actually shared");
    check(moved.report().limitedFaces>0&&moved.report().maximumFluxCorrection>0,"D0 did not exercise its geometric deficit restriction");
    for(const bool half:{true,false}){
        const auto& flux=half?shared.half:shared.full;const auto& target=half?moved.midpoint():moved.endpoint();
        check(flux.total>stock&&flux.gas==stock&&flux.liquid==flux.total-flux.gas,
            "unique outgoing phase flux did not saturate the independent stock exactly");
        check(flux.limiter>0&&flux.limiter<1,"D0 sole-face restriction diagnostic is inactive");
        check(target.plicCells()[donor].inventory[1]==0&&target.plicCells()[donor].resolved,
            "exact sole-face depletion left a positive numerical micro stock");
        for(const int cell:{donor,receiver}){
            std::array<long double,2> net{};long double total=0;
            for(const auto& face:moved.faces())if(face.negativeCell==cell||face.positiveCell==cell){
                const auto& f=half?face.half:face.full;const auto sign=face.negativeCell==cell?1.L:-1.L;
                net[0]+=sign*f.liquid;net[1]+=sign*f.gas;total+=sign*f.total;
                check(f.mass==100*f.liquid+f.gas,"single outgoing corrected face lost its same-phase material mass");
            }
            const auto& before=saved[static_cast<std::size_t>(cell)].inventory;
            const auto& after=target.plicCells()[static_cast<std::size_t>(cell)].inventory;
            const auto cc=before[0]>before[1]?1.L:0.L;
            for(unsigned phase=0;phase<2;++phase){const auto expected=before[phase]-net[phase]+(phase?1-cc:cc)*total;
                near(after[phase],expected,source.cells()[static_cast<std::size_t>(cell)].capacity,
                    "shared donor/receiver original WY phase ledger",1e-12L);}
        }
    }
    retained(source,saved,identity,clock);
    // No deficit: every pure-phase face retains its original raw flux.
    const auto pure=GeometricVofInitialization::prepare(mesh,false,[](Vector2D){return -1.;});
    const auto pureSaved=pure.plicCells();const auto pureId=pure.identity();const auto pureClock=pure.timeLayer();const auto pureMoved=advance(pure,.025);
    check(pureMoved.report().limitedFaces==0&&pureMoved.report().maximumFluxCorrection==0,"unlimited pure flow acquired a restriction");
    for(const auto& face:pureMoved.faces())for(const auto* f:{&face.half,&face.full})
        check(f->limiter==1&&f->gas==0&&f->liquid==f->total&&f->rawLiquid==f->liquid&&f->mass==100*f->liquid,
            "unlimited original pure-phase flux changed");
    retained(pure,pureSaved,pureId,pureClock);
    // A separately declared partial sweep must retain a real positive gas
    // remainder; unique incidence alone never authorizes exhaustion.
    const auto plane=GeometricVofInitialization::prepare(mesh,false,[](Vector2D p){return p.x-.423;});
    const auto planeSaved=plane.plicCells();const auto planeId=plane.identity();const auto planeClock=plane.timeLayer();const auto partial=advance(plane,.01);
    constexpr int partialCell=16*32+13;const auto partialRow=rightFace(plane,partialCell);
    for(const bool half:{true,false}){const auto& f=half?partial.faces()[partialRow].half:partial.faces()[partialRow].full;
        const auto& target=half?partial.midpoint():partial.endpoint();const auto old=planeSaved[partialCell].inventory[1];
        const auto actual=target.plicCells()[partialCell].inventory[1];
        check(f.limiter==1&&f.rawLiquid==f.liquid&&f.gas>0&&f.gas<old&&actual>0,
            "partial sole-face phase was incorrectly exhausted");
        near(actual,old-f.gas,old,"independent positive partial-sweep remainder",1e-12L);
    }
    retained(plane,planeSaved,planeId,planeClock);
    std::cout<<std::setprecision(17)<<"LOCAL_PLIC_SINGLE_OUTGOING_PASS,source=SDF,n=32,cell="<<donor
        <<",old_gas="<<stock<<",half_gas="<<moved.midpoint().plicCells()[donor].inventory[1]
        <<",full_gas="<<moved.endpoint().plicCells()[donor].inventory[1]<<",limited_faces="<<moved.report().limitedFaces
        <<",maximum_flux_correction="<<moved.report().maximumFluxCorrection<<",partial_remainder="
        <<partial.endpoint().plicCells()[partialCell].inventory[1]<<'\n';
}

void localPlicTransport()
{
    localPlicSingleOutgoingOnly();
    // BEGIN D11 SHORT MATRIX ENTRY
    localPlicAtomsOnly();
    // END D11 SHORT MATRIX ENTRY
    localPlicStabilityRates();
    for(bool rz:{false,true})for(bool amr:{false,true})for(int order:{0,1}){stage="local PLIC real sweep "+std::to_string(rz)+","+std::to_string(amr)+","+std::to_string(order);std::unique_ptr<QuadtreeMesh> mesh;if(amr)mesh=vofMaterialMesh(4);else{DictionaryMesh d;d.Nx=d.Ny=4;d.maxCellNum=4096;d.maxInitialIterations=0;mesh=std::make_unique<QuadtreeMesh>(d);}const Clock oldClock=order?Clock{1,0,.125,Position::End}:Clock{};
        const auto source=Publication::prepareLocalPlic(*mesh,domain(*mesh,rz),oldClock,l1PlaneAlpha(*mesh,rz,order));const auto a=vofMaterialNormals(*mesh,source,0),b=vofMaterialNormals(*mesh,source,order?2:0);const auto dt=static_cast<double>(.1L/vofMaterialOutflow(*mesh,source,a,b));const Clock next{oldClock.step+1,oldClock.time(),oldClock.time()+dt,Position::End};const auto slab=Slab::prepare(*mesh,source,next,a,b);const auto saved=source.plicCells();const auto* storage=source.plicCells().data();const auto id=source.identity();const auto alpha=source.alpha();const auto t=Transport::prepareLocalPlicSweep(*mesh,source,slab,100,1);const auto priorAnalytic=l1Analytic;l1Inspect(*mesh,source,slab,t,order);check(l1Analytic>priorAnalytic,"short scene missed qualified independent axial strip");
        if(amr&&order==0)l1MaterialConsumer(*mesh,source,slab,t);
        if(l1Cases==0){const auto reject=[&](auto action,const char* why){rejects(action,why);++l1Guards;};
            const auto courant=[&](auto action,const char* expected){bool caught=false;
                try{action();}catch(const GeometricVofCourantExceeded& error){caught=true;
                    check(std::string(error.what())==std::string("[Interface/GeometricVofTransport] ")+expected,
                        "typed Courant rejection changed its original complete diagnostic");}
                check(caught,"actual transport CFL failure lost its specific Courant exception type");};
            const auto hardFailure=[&](auto action,const char* why){bool ordinary=false;
                try{action();}catch(const GeometricVofCourantExceeded&){throw std::runtime_error(std::string(why)+" was incorrectly marked recoverable Courant");}
                catch(const HydroTwoPhaseStepRetry&){throw std::runtime_error(std::string(why)+" was incorrectly converted into StepRetry");}
                catch(const std::invalid_argument&){ordinary=true;}
                check(ordinary,std::string("missing ordinary invalid_argument: ")+why);};
            reject([&]{slab.beginVelocity();},"local has no tangential begin QR");reject([&]{slab.endVelocity();},"local has no tangential end QR");reject([&]{Transport::prepare(*mesh,source,slab,100,1);},"fan transport cannot consume LocalPlic");const auto e=faceEndpoints(*mesh,source.domain().face(*mesh,9));reject([&]{t.restrictNativeFaces(*mesh,{{-1,e[0],e[1]}});},"local invalid native restriction row");reject([&]{t.measureInteriorFaces(*mesh,{{0,{0,.0625L},{0,.1875L}}});},"local interior cannot lie on a cell boundary");
            auto wrong=next;++wrong.step;reject([&]{Slab::prepare(*mesh,source,wrong,a,b);},"local fullclock step mismatch");auto broken=a;int row=0;while(row<source.domain().nativeFaceCount()){const auto f=source.domain().face(*mesh,row);if(f.negativeCell>=0&&f.positiveCell>=0&&f.area>0)break;++row;}check(row<source.domain().nativeFaceCount(),"no actual shared face for bad divergence guard");broken[row]+=.1;const auto badDiv=Slab::prepare(*mesh,source,next,broken,b);reject([&]{Transport::prepareLocalPlicSweep(*mesh,source,badDiv,100,1);},"local rejects non-divergence-free endpoint");
            GeometricVofTransportConfig cfg;cfg.maximumLocalCourant=.51L;hardFailure([&]{Transport::prepareLocalPlicSweep(*mesh,source,slab,100,1,cfg);},"local Courant may not exceed frozen .5");++l1Guards;auto large=a;for(auto& v:large)v*=32;const auto badColumn=Slab::prepare(*mesh,source,next,large,large);courant([&]{Transport::prepareLocalPlicSweep(*mesh,source,badColumn,100,1);},"local donor subface-column Courant exceeded");++l1Guards;
            cfg={};cfg.maximumLocalCourant=.001L;reject([&]{Transport::prepareLocalPlicSweep(*mesh,source,slab,100,1,cfg);},"caller may tighten actual local CFL");const auto foreign=Publication::prepareLocalPlic(*mesh,source.domain(),source.timeLayer(),source.alpha());
            hardFailure([&]{slab.validateFor(*mesh,foreign);},"local sameclock foreign source binding");++l1Guards;
            hardFailure([&]{Transport::prepareLocalPlicSweep(*mesh,foreign,slab,100,1);},"transport foreign source owner");

            // Same original 4x4 XY scene and discrete streamfunction, no new
            // fixture. On its exact grid, max face/max total flow is 144/175.
            // Six times the original .1 total CFL gives .6 complete CFL,
            // while every donor column is at most 432/875, strictly below .5.
            check(!rz&&!amr&&order==0,"typed retry witness must reuse the existing first uniform XY scene");
            auto completeVelocity=a;for(auto& v:completeVelocity)v*=6;
            const auto completeSlab=Slab::prepare(*mesh,source,next,completeVelocity,completeVelocity);
            std::vector<long double> incoming(source.cells().size()),outgoing(incoming.size());long double maxColumn=0;
            for(int row=0;row<source.domain().nativeFaceCount();++row){const auto face=source.domain().face(*mesh,row);
                const auto q=static_cast<long double>(face.area)*completeVelocity[row]*dt;if(q==0)continue;
                const int donor=q>0?face.negativeCell:face.positiveCell,receiver=q>0?face.positiveCell:face.negativeCell;
                check(donor>=0&&receiver>=0,"existing closed retry streamfunction acquired exterior flux");
                const auto volume=std::abs(q);outgoing[donor]+=volume;incoming[receiver]+=volume;
                // Uniform full native subface: donor column is exactly its cell.
                maxColumn=std::max(maxColumn,volume/source.cells()[donor].capacity);}
            long double maxComplete=0;
            for(std::size_t ci=0;ci<incoming.size();++ci)maxComplete=std::max(maxComplete,
                std::max(incoming[ci],outgoing[ci])/source.cells()[ci].capacity);
            check(maxColumn<.5L&&maxComplete>.5L&&maxComplete/2<.5L,
                "independent complete-cell witness must pass individual columns and recover under half dt");
            near(maxColumn,432.L/875,1,"independent original streamfunction column CFL",1e-12L);
            near(maxComplete,.6L,1,"independent original streamfunction complete CFL",1e-12L);
            courant([&]{Transport::prepareLocalPlicSweep(*mesh,source,completeSlab,100,1);},
                "local complete-cell incoming/outgoing Courant exceeded");
            auto halfEnd=next;halfEnd.end=halfEnd.begin+(next.end-next.begin)/2;
            const auto halfSlab=Slab::prepare(*mesh,source,halfEnd,completeVelocity,completeVelocity);
            const auto halfTransport=Transport::prepareLocalPlicSweep(*mesh,source,halfSlab,100,1);
            check(halfTransport.sourceIdentity()==source.identity()&&sameClock(halfTransport.endpoint().timeLayer(),halfEnd)&&
                halfTransport.report().maximumLocalColumnCfl<=.5L&&halfTransport.report().maximumIncomingCfl<=.5L,
                "same source and original transport preparation must succeed after dt is halved");

            const auto retryMetrics=vofMaterialMetrics(*mesh,source);const auto& retryLayout=retryMetrics.layout();
            const auto retryPhase=HydroPhaseGeometry::build(*mesh,retryMetrics,source);
            std::vector<double> retryNormal(retryLayout.velocityCount());
            for(int u=0;u<retryLayout.velocityCount();++u)retryNormal[u]=completeVelocity[retryLayout.nativeFace(u)];
            const HydroNativeValues retryValues(retryLayout,retryNormal,std::vector<double>(retryLayout.pressureCount()));
            const auto retryAccepted=HydroTwoPhaseState::createVof(*mesh,retryPhase,retryValues),retainedAccepted=retryAccepted;
            const auto retainedClock=retryAccepted.timeLayer();const auto retainedNormal=retryAccepted.values().normalVelocity();
            const auto retainedPressure=retryAccepted.values().pressure();const auto retainedStocks=source.plicCells();
            // Homogeneous density and exact mu=sigma=0 isolate the retry
            // boundary, retaining the original NS/AC controls and real VOF.
            HydroVofAcConfig retryConfig;retryConfig.liquidDensity=retryConfig.gasDensity=1;
            retryConfig.momentumDensityRepresentation=HydroMomentumMass::DensityRepresentation::CellAverageProjection;
            bool transportRetry=false;
            try{HydroVofStep::prepare(*mesh,retryAccepted,next,retryConfig);}
            catch(const HydroTwoPhaseStepRetry& error){transportRetry=true;
                check(error.reason()==HydroTwoPhaseRetryReason::Transport&&
                    std::string(error.what()).find("local complete-cell incoming/outgoing Courant exceeded")!=std::string::npos,
                    "Hydro must translate only the specific transport CFL failure and retain its reason");}
            check(transportRetry,"actual Hydro VOF prepare did not expose the existing typed StepRetry boundary");
            const auto unchanged=[&]{
                check(retryAccepted.sameIdentity(retainedAccepted)&&sameClock(retryAccepted.timeLayer(),retainedClock)&&
                    retryAccepted.values().normalVelocity()==retainedNormal&&retryAccepted.values().pressure()==retainedPressure&&
                    retryAccepted.phaseGeometry().vofIdentity()==source.identity()&&
                    retryAccepted.phaseGeometry().vofPublication().plicCells().data()==storage,
                    "failed/reprepared Hydro candidate changed accepted clock, velocity, pressure or source owner");
                for(std::size_t ci=0;ci<retainedStocks.size();++ci)check(
                    retryAccepted.phaseGeometry().vofPublication().plicCells()[ci].inventory==retainedStocks[ci].inventory,
                    "typed retry changed an accepted independent phase stock pair");};
            unchanged();
            // Half dt recovery is established above through the unchanged
            // original Interface prepare; unrelated NS/AC acceptance is separate.
            auto invalidRetryConfig=retryConfig;invalidRetryConfig.transport.maximumLocalCourant=.51L;
            hardFailure([&]{HydroVofStep::prepare(*mesh,retryAccepted,halfEnd,invalidRetryConfig);},"Hydro transport configuration");
            auto wrongRetryClock=halfEnd;++wrongRetryClock.step;
            hardFailure([&]{HydroVofStep::prepare(*mesh,retryAccepted,wrongRetryClock,retryConfig);},"Hydro complete clock");
            unchanged();
            std::cout<<"LOCAL_PLIC_TYPED_COURANT,donor_typed,1,complete_typed,1,hydro_transport_retry,1,half_dt_transport_prepared,1"
                <<",max_column,"<<maxColumn<<",max_complete,"<<maxComplete<<",accepted_unchanged,1,PASS\n";
        }
        check(source.identity()==id&&source.plicCells().data()==storage&&source.alpha()==alpha&&sameClock(source.timeLayer(),oldClock),"local candidate/failure changed accepted source owner");for(std::size_t i=0;i<saved.size();++i)check(source.plicCells()[i].inventory==saved[i].inventory,"local candidate/failure changed source stock pair");++l1Cases;
    }
    check(l1Cases==8&&l1Guards==11&&l1CoarseFine>0&&l1Analytic>0,"local short coverage missing");std::cout<<std::setprecision(17)<<"LOCAL_PLIC_TRANSPORT_SUMMARY,"<<l1Cases<<','<<l1LedgerCells<<','<<l1Fluxes<<','<<l1Guards<<','<<l1Analytic<<','<<l1CoarseFine<<','<<l1PeakLedger<<','<<l1PeakConservation<<','<<l1PeakStrip<<",PASS\n";
}
long double l1EllipseArea(long double x0,long double x1,long double y0,long double y1)
{
    constexpr long double cx=.5L,cy=.5L,rx=.29L,ry=.22L;const auto left=std::max(x0,cx-rx),right=std::min(x1,cx+rx);if(right<=left)return 0;std::vector<long double> cuts{left,right};if(cx>left&&cx<right)cuts.push_back(cx);for(auto y:{y0,y1}){const auto q=std::abs((y-cy)/ry);if(q<1){const auto d=rx*std::sqrt(1-q*q);for(auto x:{cx-d,cx+d})if(x>left&&x<right)cuts.push_back(x);}}std::sort(cuts.begin(),cuts.end());cuts.erase(std::unique(cuts.begin(),cuts.end()),cuts.end());const auto primitive=[&](long double x){const auto t=std::clamp((x-cx)/rx,-1.L,1.L);return rx*(t*std::sqrt(std::max(0.L,1-t*t))+std::asin(t))/2;};long double area=0;for(std::size_t i=1;i<cuts.size();++i){const auto a=cuts[i-1],b=cuts[i],t=((a+b)/2-cx)/rx,q=std::sqrt(std::max(0.L,1-t*t));const auto upper=cy+ry*q,lower=cy-ry*q;if(std::min(y1,upper)<=std::max(y0,lower))continue;long double c=0,s=0;if(y1<upper)c+=y1;else{c+=cy;s+=ry;}if(y0>lower)c-=y0;else{c-=cy;s+=ry;}area+=c*(b-a)+s*(primitive(b)-primitive(a));}return area;
}
void localPlicCost(int n,int steps,double dt)
{
    check(n>=4&&n<=256&&steps>=1&&steps<=100&&std::isfinite(dt)&&dt>0,"invalid local cost dimensions/steps/dt");int threads=1;
#ifdef _OPENMP
    omp_set_dynamic(0);omp_set_num_threads(8);threads=omp_get_max_threads();check(omp_get_thread_limit()>=8,"local cost requires OMP_THREAD_LIMIT>=8");
#endif
    const auto start=std::chrono::steady_clock::now();DictionaryMesh d;d.Nx=d.Ny=n;d.maxCellNum=4*n*n;d.maxInitialIterations=0;auto mesh=std::make_unique<QuadtreeMesh>(d);std::vector<double> alpha(mesh->cellIndexNum());for(int i=0;i<mesh->cellIndexNum();++i){const auto b=bounds(*mesh,i);const auto c=(b.x1-b.x0)*(b.y1-b.y0);const auto area=l1EllipseArea(b.x0,b.x1,b.y0,b.y1);alpha[i]=static_cast<double>(area/c);check(alpha[i]>=0&&alpha[i]<=1,"independent initial ellipse fraction out of bounds");}
    auto state=Publication::prepareLocalPlic(*mesh,domain(*mesh,false),{},alpha);const auto initial=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();const auto normal=vofMaterialNormals(*mesh,state,0);long double initialL=0,initialG=0;for(auto q:state.plicCells()){initialL+=q.inventory[0];initialG+=q.inventory[1];}std::vector<double> elapsed;double slabSeconds=0,totalTransport=0;std::size_t sweeps=0,clips=0,unresolved=0;long double minStock=std::numeric_limits<long double>::infinity(),maxFraction=0,maximumCompression=0;
    std::cout<<std::setprecision(17)<<"LOCAL_PLIC_COST_STEP_HEADER,n,step,time,slab_seconds,transport_seconds,liquid,gas,mass,minimum_stock,maximum_phase_fraction,mixed_cells,directional_sweeps,strip_clips,unresolved_donors,compression_residual\n";
    for(int step=1;step<=steps;++step){const Clock next{step,(step-1)*dt,step*dt,Position::End};const auto slabStart=std::chrono::steady_clock::now();const auto slab=Slab::prepare(*mesh,state,next,normal,normal);const auto slabTime=std::chrono::duration<double>(std::chrono::steady_clock::now()-slabStart).count();const auto transportStart=std::chrono::steady_clock::now();const auto transport=Transport::prepareLocalPlicSweep(*mesh,state,slab,100,1);const auto seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-transportStart).count();const auto& candidate=transport.endpoint();long double liquid=0,gas=0,capacity=0;std::size_t mixed=0;for(std::size_t i=0;i<candidate.plicCells().size();++i){const auto& q=candidate.plicCells()[i];const auto c=candidate.cells()[i].capacity;for(auto value:q.inventory){check(value>=0&&value<=c,"local cost stock out of bounds");minStock=std::min(minStock,value);maxFraction=std::max(maxFraction,value/c);}liquid+=q.inventory[0];gas+=q.inventory[1];capacity+=c;if(q.inventory[0]>0&&q.inventory[1]>0)++mixed;}near(liquid,initialL,capacity,"real production-window liquid conservation",1e-12L);near(gas,initialG,capacity,"real production-window gas conservation",1e-12L);slabSeconds+=slabTime;totalTransport+=seconds;elapsed.push_back(seconds);sweeps+=transport.report().localDirectionalSweeps;clips+=transport.report().localStripClips;unresolved+=transport.report().localUnresolvedDonors;maximumCompression=std::max(maximumCompression,transport.report().maximumCompressionResidual);
        std::cout<<"LOCAL_PLIC_COST_STEP,"<<n<<','<<step<<','<<next.end<<','<<slabTime<<','<<seconds<<','<<liquid<<','<<gas<<','<<100*liquid+gas<<','<<minStock<<','<<maxFraction<<','<<mixed<<','<<transport.report().localDirectionalSweeps<<','<<transport.report().localStripClips<<','<<transport.report().localUnresolvedDonors<<','<<transport.report().maximumCompressionResidual<<'\n';state=candidate;
    }
    std::sort(elapsed.begin(),elapsed.end());const auto median=elapsed.size()%2?elapsed[elapsed.size()/2]:(elapsed[elapsed.size()/2-1]+elapsed[elapsed.size()/2])/2;const auto p95=elapsed[static_cast<std::size_t>(std::ceil(.95*elapsed.size()))-1];
    std::cout<<"LOCAL_PLIC_COST_SUMMARY_HEADER,n,steps,physical_time,OMP_threads,initial_seconds,slab_total_seconds,transport_median_seconds,transport_p95_seconds,transport_total_seconds,elapsed_total_seconds,cells,faces,directional_sweeps,strip_clips,unresolved_donors,minimum_stock,maximum_phase_fraction,maximum_compression_residual\n";
    std::cout<<"LOCAL_PLIC_COST_SUMMARY,"<<n<<','<<steps<<','<<steps*dt<<','<<threads<<','<<initial<<','<<slabSeconds<<','<<median<<','<<p95<<','<<totalTransport<<','<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<','<<mesh->cellIndexNum()<<','<<state.domain().nativeFaceCount()<<','<<sweeps<<','<<clips<<','<<unresolved<<','<<minStock<<','<<maxFraction<<','<<maximumCompression<<",PASS,transport_only_no_Hydro_AC_claim\n";
}
}
// END L1 LOCAL PLIC SWEEP AND COST CONTRACT

namespace
{

// BEGIN RZ C1 SAME-PLANE ORDINARY MERIDIONAL AREA CONTRACT
// Independent Fubini oracle of the ACTUAL published plane. Solve the
// principal coordinate, split at its two rectangle-edge intersections,
// and integrate each clipped affine strip exactly. No Core cut/measure,
// polygon clipping, inventory inversion or capacity-minus-phase is used.
std::array<long double,2> meridionalPlaneStripReference(const GeometricVofPlicCell& cell,
    const LAMP::Geometry::PlanarCellCut::Rectangle& box,bool liquid)
{
    const auto length=std::hypot(static_cast<long double>(cell.normal.x),static_cast<long double>(cell.normal.y));
    check(std::isfinite(length)&&length>0,"independent actual plane has invalid normal");
    const long double sign=liquid?1.L:-1.L,nx=sign*cell.normal.x/length,ny=sign*cell.normal.y/length,offset=sign*cell.offset;
    const long double cx=static_cast<long double>(cell.rectangle.lower.x)+(static_cast<long double>(cell.rectangle.upper.x)-cell.rectangle.lower.x)/2;
    const long double cy=static_cast<long double>(cell.rectangle.lower.y)+(static_cast<long double>(cell.rectangle.upper.y)-cell.rectangle.lower.y)/2;
    const bool solveRadial=std::abs(nx)>=std::abs(ny);
    const long double nu=solveRadial?nx:ny,nt=solveRadial?ny:nx;
    const long double cu=solveRadial?cx:cy,ct=solveRadial?cy:cx;
    const long double u0=solveRadial?box.lower.x:box.lower.y,u1=solveRadial?box.upper.x:box.upper.y;
    const long double t0=solveRadial?box.lower.y:box.lower.x,t1=solveRadial?box.upper.y:box.upper.x;
    check(nu!=0&&u1>u0&&t1>t0,"independent actual plane strip domain");
    std::vector<long double> cuts{t0,t1};
    if(nt!=0)for(auto bound:{u0,u1}){const auto t=ct+(offset-nu*(bound-cu))/nt;
        if(t>t0&&t<t1)cuts.push_back(t);}
    std::sort(cuts.begin(),cuts.end());cuts.erase(std::unique(cuts.begin(),cuts.end()),cuts.end());
    std::array<long double,2> result{}; // dA, integral(r dA)
    for(std::size_t i=1;i<cuts.size();++i){
        const auto lower=cuts[i-1],upper=cuts[i],span=upper-lower;
        const auto crossing=[&](long double t){return cu+(offset-nt*(t-ct))/nu;};
        // These min/max operations intersect the analytic strip with the
        // query rectangle; they never alter a published phase measure.
        const auto c0=std::max(u0,std::min(u1,crossing(lower))),c1=std::max(u0,std::min(u1,crossing(upper)));
        const long double a0=nu>0?u0:c0,b0=nu>0?c0:u1,a1=nu>0?u0:c1,b1=nu>0?c1:u1;
        const auto d0=b0-a0,d1=b1-a1;
        result[0]+=span*(d0+d1)/2;
        if(solveRadial)result[1]+=span*(d0*(b0+a0)+d0*b1+a0*d1+d1*(b1+a1))/6;
        else result[1]+=span*(lower*(d0+d1)/2+span*(d0+2*d1)/6);
    }
    check(std::isfinite(result[0])&&std::isfinite(result[1])&&result[0]>=0&&result[1]>=0,
        "independent actual plane strip integral invalid");return result;
}
void localPlicMeridionalOnly()
{
    stage="RZ C1 independent ordinary meridional area";
    std::size_t checked=0,mixed=0,axisBoxes=0,weightedDifference=0;
    for(bool rz:{false,true}){
        auto mesh=vofMaterialMesh(4);const auto native=domain(*mesh,rz);
        constexpr long double targetRadius=.587L;
        std::vector<std::array<long double,2>> stocks(mesh->cellIndexNum());
        for(int ci=0;ci<mesh->cellIndexNum();++ci){
            const auto& cell=mesh->cells(mesh->cellIndexToLabel[ci]);
            const auto lo=mesh->verts(cell.cellToVert[0]).coord(),hi=mesh->verts(cell.cellToVert[2]).coord();
            const auto left=std::min(static_cast<long double>(hi.x),targetRadius);
            const auto right=std::max(static_cast<long double>(lo.x),targetRadius);
            // Two independent intervals: never capacity minus the liquid stock.
            const auto measure=[&](long double a,long double b){
                return b<=a?0.L:(hi.y-lo.y)*(rz?twoPi*(b*b-a*a)/2:b-a);};
            stocks[ci]={measure(lo.x,left),measure(right,hi.x)};
        }
        const Clock clock{3,.125,.25,Position::Midpoint};
        const auto source=Publication::prepareLocalPlicInventory(*mesh,native,clock,stocks,{});
        const auto metrics=vofMaterialMetrics(*mesh,source);
        const auto phase=HydroPhaseGeometry::build(*mesh,metrics,source);
        std::vector<HydroVofBox> boxes;
        for(int ci=0;ci<mesh->cellIndexNum();++ci){
            const auto& b=source.plicCells()[ci].rectangle;
            boxes.push_back({ci,b});
            const auto dx=b.upper.x-b.lower.x,dy=b.upper.y-b.lower.y;
            boxes.push_back({ci,{{b.lower.x+.125*dx,b.lower.y+.25*dy},
                {b.lower.x+.8125*dx,b.lower.y+.875*dy}}});
        }
        const auto areas=phase.measureVofMeridionalBoxes(*mesh,boxes);
        const auto physical=phase.measureVofBoxes(*mesh,boxes);
        check(areas.size()==boxes.size()&&physical.size()==boxes.size(),"meridional batch result count");
        for(std::size_t q=0;q<boxes.size();++q){
            const auto& cell=source.plicCells()[boxes[q].nativeCell];const auto& box=boxes[q].rectangle;
            const auto area=static_cast<long double>(box.upper.x-box.lower.x)*(box.upper.y-box.lower.y);
            long double liquid=0,gas=0,liquidPhysical=0,gasPhysical=0;
            if(cell.measured[0]==0){gas=area;gasPhysical=vofMaterialCapacity(box.lower,box.upper,rz);}
            else if(cell.measured[1]==0){liquid=area;liquidPhysical=vofMaterialCapacity(box.lower,box.upper,rz);}
            else{
                if(!(cell.resolved&&cell.normal.x==1&&cell.normal.y==0))std::cout<<"MERIDIONAL_PLANE_DIAGNOSTIC,rz="<<rz<<",native="<<boxes[q].nativeCell<<",resolved="<<cell.resolved<<",nx="<<cell.normal.x<<",ny="<<cell.normal.y<<",offset="<<cell.offset<<",r0="<<cell.rectangle.lower.x<<",r1="<<cell.rectangle.upper.x<<",z0="<<cell.rectangle.lower.y<<",z1="<<cell.rectangle.upper.y<<std::endl;
                check(cell.resolved,"ordinary-area oracle requires its actual resolved mixed plane");
                const auto l=meridionalPlaneStripReference(cell,box,true),g=meridionalPlaneStripReference(cell,box,false);
                liquid=l[0];gas=g[0];liquidPhysical=rz?twoPi*l[1]:liquid;gasPhysical=rz?twoPi*g[1]:gas;
                ++mixed;
                if(rz&&liquid>0&&gas>0){
                    const auto fake=physical[q].liquid/(twoPi*(box.lower.x+box.upper.x)/2);
                    check(std::abs(fake-liquid)>1e-6L*area,"ordinary area must differ from volume divided by centre radius");
                    ++weightedDifference;
                }
            }
            near(areas[q].liquid,liquid,area,"same-plane ordinary liquid",1e-12L);
            near(areas[q].gas,gas,area,"same-plane ordinary gas",1e-12L);
            near(areas[q].liquid+areas[q].gas,area,area,"independent ordinary phase closure",1e-12L);
            const auto cap=vofMaterialCapacity(box.lower,box.upper,rz);
            near(physical[q].liquid,liquidPhysical,cap,"retained physical liquid",1e-12L);
            near(physical[q].gas,gasPhysical,cap,"retained physical gas",1e-12L);
            axisBoxes+=rz&&box.lower.x==0;++checked;
        }
        for(std::size_t ci=0;ci<stocks.size();++ci)
            check(source.plicCells()[ci].inventory==stocks[ci],"meridional queries preserve both exact source stocks");
        phase.validateVofFor(*mesh,source.identity(),clock);
        auto wrong=clock;wrong.begin=0;
        rejects([&]{phase.validateVofFor(*mesh,source.identity(),wrong);},"clock");
        auto bad=boxes.front();bad.rectangle.lower.x-=.01;
        rejects([&]{phase.measureVofMeridionalBoxes(*mesh,{bad});},"original native leaf");
    }
    check(mixed>0&&axisBoxes>0&&weightedDifference>0,"meridional matrix lacks mixed/axis/physical-weight witnesses");
    std::cout<<"LOCAL_PLIC_MERIDIONAL_PASS boxes="<<checked<<" mixed="<<mixed
        <<" axis="<<axisBoxes<<" weighted_difference="<<weightedDifference<<'\n';
}
// END RZ C1 SAME-PLANE ORDINARY MERIDIONAL AREA CONTRACT

// F1 Interface-only open transport; no Hydro boundary solver is exercised.
void localPlicOpenBoundaryOnly()
{
    using Boundary=Slab::Boundary;using Mode=Boundary::Mode;
    unsigned scenes=0,boundaryFaces=0,negativeControls=0;
    const auto rulesFor=[](const QuadtreeMesh& mesh,const Publication& source,Mode low,Mode high){
        std::vector<Mode> rules(source.domain().nativeFaceCount(),Mode::Interior);
        for(int row=0;row<source.domain().nativeFaceCount();++row){const auto f=source.domain().face(mesh,row);
            if(f.negativeCell>=0&&f.positiveCell>=0)continue;
            rules[row]=f.horizontal||f.area==0?Mode::Closed:(f.positiveCell>=0?low:high);}
        return rules;
    };
    const auto normals=[](const QuadtreeMesh& mesh,const Publication& source,double amplitude){
        std::vector<double> result(source.domain().nativeFaceCount());
        for(int row=0;row<source.domain().nativeFaceCount();++row){const auto f=source.domain().face(mesh,row);
            if(!f.horizontal&&f.area>0)result[row]=source.domain().axisymmetric()?amplitude/f.center.x:amplitude;}
        return result;
    };
    const auto fields=[](const GeometricVofIntegratedFlux& f){
        return std::array<long double,12>{f.total,f.positiveTotal,f.negativeTotal,f.rawLiquid,f.liquid,f.gas,f.mass,
            f.positiveLiquid,f.negativeLiquid,f.positiveGas,f.negativeGas,f.limiter};
    };
    const auto inspect=[&](const QuadtreeMesh& mesh,const Publication& source,const Transport& transport){
        check(transport.advector().boundary()!=nullptr,"open transport lost its retained Boundary");
        transport.validateFor(mesh,source,transport.advector(),100,1);
        for(bool half:{true,false}){
            const auto& ledger=half?transport.report().boundaryHalf:transport.report().boundaryFull;
            const auto& candidate=half?transport.midpoint():transport.endpoint();
            check(ledger.enabled&&ledger.sourceIdentity==source.identity()&&sameClock(ledger.sourceTime,source.timeLayer())
                &&sameClock(ledger.targetTime,candidate.timeLayer()),"boundary ledger owner/fullclock");
            std::array<long double,2> old{},now{},in{},out{};
            for(std::size_t i=0;i<source.plicCells().size();++i)for(unsigned p=0;p<2;++p){
                old[p]+=source.plicCells()[i].inventory[p];now[p]+=candidate.plicCells()[i].inventory[p];}
            for(const auto& f:transport.faces()){
                const auto& q=half?f.half:f.full;
                near(q.mass,100*q.liquid+q.gas,std::max(std::abs(q.mass),1e-30L),"open same shared density-weighted flux",1e-12L);
                if(f.negativeCell>=0&&f.positiveCell>=0)continue;
                ++boundaryFaces;const auto sign=f.negativeCell>=0?1.L:-1.L;
                for(unsigned p=0;p<2;++p){const auto positive=p==0?q.positiveLiquid:q.positiveGas;
                    const auto negative=p==0?q.negativeLiquid:q.negativeGas;
                    out[p]+=sign>0?positive:negative;in[p]+=sign>0?negative:positive;}
            }
            long double c=0;for(const auto& cell:source.cells())c+=cell.capacity;
            for(unsigned p=0;p<2;++p){
                near(now[p]-old[p]+out[p]-in[p],0,c,"independent physical phase boundary balance",1e-12L);
                near(ledger.source[p],old[p],c,"ledger independent source",1e-12L);
                near(ledger.target[p],now[p],c,"ledger independent target",1e-12L);
                near(ledger.inflow[p],in[p],c,"ledger independent inflow",1e-12L);
                near(ledger.outflow[p],out[p],c,"ledger independent outflow",1e-12L);
                near(ledger.physicalResidual[p]-ledger.compression[p]-ledger.roundingError[p],
                    ledger.algebraicResidual[p],c,"WY/RN algebra is distinct from physical boundary gate",1e-12L);
            }
            near(100*(now[0]-old[0]+out[0]-in[0])+now[1]-old[1]+out[1]-in[1],0,100*c,
                "independent physical mass boundary balance",1e-12L);
            check(ledger.physicalRelativeConservation<=1e-12L&&ledger.algebraicRelativeConservation<=1e-12L
                &&ledger.massRelativeConservation<=1e-12L,"open physical/algebraic gates");
        }
    };
    stage="F1 XY actual 2:1 liquid reservoir into gas";
    auto mesh=vofMaterialMesh(4);const auto source=Publication::prepareLocalPlic(*mesh,domain(*mesh,false),{},
        std::vector<double>(mesh->cellIndexNum(),0));
    const auto saved=source.plicCells();const auto* storage=source.plicCells().data();const auto oldAlpha=source.alpha();
    auto rules=rulesFor(*mesh,source,Mode::PureLiquidReservoir,Mode::PureGasReservoir);
    const auto boundary=Boundary::prepare(*mesh,source,rules),copy=boundary;
    const auto another=Boundary::prepare(*mesh,source,rules);
    check(copy.sameInstance(boundary)&&!another.sameInstance(boundary),"Boundary sameInstance is owner identity");
    const auto a=normals(*mesh,source,.125),b=normals(*mesh,source,.25);
    const auto rates=Transport::localPlicBoundaryCourantRates(*mesh,source,boundary,a);
    check(rates.maximumIncoming>0&&rates.maximumOutflow>0&&rates.maximumLocalColumn>0,"pre-dt open rates unavailable");
    const double dt=1./32;const Clock end{1,0,dt,Position::End};
    const auto slab=Slab::prepareBoundary(*mesh,source,boundary,end,a,b);
    check(slab.boundary()&&slab.boundary()->sameInstance(boundary),"slab substituted boundary owner");
    const auto transport=Transport::prepareLocalPlicBoundarySweep(*mesh,source,slab,100,1);
    inspect(*mesh,source,transport);++scenes;
    for(const auto& f:transport.faces())if(f.negativeCell<0&&!source.domain().face(*mesh,f.nativeRow).horizontal){
        const auto face=source.domain().face(*mesh,f.nativeRow);
        near(f.half.total,face.area*dt*(3.L*.125L+.25L)/8,face.area*dt,
            "analytic half boundary linear endpoint integral",1e-12L);
        near(f.full.total,face.area*dt*(.125L+.25L)/2,face.area*dt,
            "analytic full boundary linear endpoint integral",1e-12L);
        check(f.full.gas==0&&f.full.liquid==f.full.total&&f.full.negativeTotal==0,"pure reservoir independent phase pair");}
    check(transport.report().boundaryFull.inflow[0]>0&&transport.report().boundaryFull.inflow[1]==0
        &&transport.report().boundaryFull.outflow[0]==0&&transport.report().boundaryFull.outflow[1]>0,
        "initial gas outflow must use actual inside donor, not declared gas inlet law");
    const auto reject=[&](auto&& action,const std::string& why){rejects(action,why);++negativeControls;};
    reject([&]{Transport::prepareLocalPlicSweep(*mesh,source,slab,100,1);},"closed sweep cannot consume open slab");
    reject([&]{Slab::prepare(*mesh,source,end,a,b);},"old slab still rejects nonzero exterior");
    reject([&]{Transport::localPlicCourantRates(*mesh,source,a);},"old rates still reject nonzero exterior");
    auto reversed=b;for(auto& value:reversed)value=-value;
    reject([&]{const auto reversedSlab=Slab::prepareBoundary(*mesh,source,boundary,end,a,reversed);
        Transport::prepareLocalPlicBoundarySweep(*mesh,source,reversedSlab,100,1);},"same-slab mismatched pure-reservoir reversal rejected");
    const auto zero=normals(*mesh,source,0);
    const auto oneZero=Slab::prepareBoundary(*mesh,source,boundary,end,a,zero);
    check(oneZero.boundary()!=nullptr,"one-zero endpoint incorrectly rejected");
    auto wrong=end;++wrong.step;
    reject([&]{Slab::prepareBoundary(*mesh,source,boundary,wrong,a,b);},"next fullclock mismatch");
    const auto foreign=Publication::prepareLocalPlic(*mesh,source.domain(),source.timeLayer(),source.alpha());
    reject([&]{boundary.validateFor(*mesh,foreign);},"same-clock foreign publication");
    auto otherMesh=vofMaterialMesh(6);
    reject([&]{boundary.validateFor(*otherMesh,source);},"changed mesh content");
    auto badRules=rules;int inside=-1,exterior=-1,wall=-1;
    for(int row=0;row<source.domain().nativeFaceCount();++row){const auto f=source.domain().face(*mesh,row);
        if(f.negativeCell>=0&&f.positiveCell>=0)inside=row;
        else {exterior=row;if(f.horizontal)wall=row;}}
    check(inside>=0&&exterior>=0&&wall>=0,"missing boundary negative witnesses");
    badRules[inside]=Mode::Closed;
    reject([&]{Boundary::prepare(*mesh,source,badRules);},"interior rule cannot masquerade as Closed");
    badRules=rules;badRules[exterior]=Mode::Interior;
    reject([&]{Boundary::prepare(*mesh,source,badRules);},"exterior needs explicit reservoir or Closed");
    auto wallNormal=a;wallNormal[wall]=.125;
    reject([&]{Slab::prepareBoundary(*mesh,source,boundary,end,wallNormal,b);},"nonzero Closed face");
    reject([&]{Boundary::prepare(*mesh,transport.midpoint(),rules);},"midpoint is not accepted source for rates");
    auto large=end;large.end=1;
    const auto largeSlab=Slab::prepareBoundary(*mesh,source,boundary,large,a,b);
    bool courant=false;
    try{Transport::prepareLocalPlicBoundarySweep(*mesh,source,largeSlab,100,1);}
    catch(const GeometricVofCourantExceeded&){courant=true;}
    check(courant,"incoming/outgoing/column excess lost typed Courant gate");++negativeControls;
    // A zero net divergence must not cancel two real incoming/outgoing sides.
    auto cornerRules=rules;
    for(auto& rule:cornerRules)if(rule!=Mode::Interior)rule=Mode::PureGasReservoir;
    const auto cornerBoundary=Boundary::prepare(*mesh,source,cornerRules);
    std::vector<double> diagonal(source.domain().nativeFaceCount(),.125);
    const Clock cornerEnd{1,0,.375,Position::End};
    const auto cornerSlab=Slab::prepareBoundary(*mesh,source,cornerBoundary,cornerEnd,diagonal,diagonal);
    bool grossCourant=false;
    try{Transport::prepareLocalPlicBoundarySweep(*mesh,source,cornerSlab,100,1);}
    catch(const GeometricVofCourantExceeded& error){grossCourant=
        std::string(error.what()).find("local complete-cell incoming/outgoing Courant exceeded")!=std::string::npos;}
    check(grossCourant,"gross complete-cell CFL must reject although per-column <=.5 and net divergence is zero");++negativeControls;
    stage="F1 source End and cross-step backflow";
    const auto& accepted=transport.endpoint();
    const auto backRules=rulesFor(*mesh,accepted,Mode::PureGasReservoir,Mode::PureGasReservoir);
    const auto back=Boundary::prepare(*mesh,accepted,backRules);
    const auto reverseA=normals(*mesh,accepted,-.125),reverseB=normals(*mesh,accepted,-.25);
    const auto endRates=Transport::localPlicBoundaryCourantRates(*mesh,accepted,back,reverseA);
    check(endRates.maximumIncoming>0,"accepted End pre-dt rate query");
    const Clock next{2,dt,2*dt,Position::End};
    const auto reverseSlab=Slab::prepareBoundary(*mesh,accepted,back,next,reverseA,reverseB);
    const auto backwards=Transport::prepareLocalPlicBoundarySweep(*mesh,accepted,reverseSlab,100,1);
    inspect(*mesh,accepted,backwards);++scenes;
    check(backwards.report().boundaryFull.inflow[0]==0&&backwards.report().boundaryFull.inflow[1]>0
        &&backwards.report().boundaryFull.outflow[0]>0,"backflow gas reservoir must not overwrite inside liquid donor");
    stage="F1 mixed real donor with pure reservoir";
    const auto mixed=Publication::prepareLocalPlic(*mesh,domain(*mesh,false),{},
        std::vector<double>(mesh->cellIndexNum(),.5));
    const auto mixedBoundary=Boundary::prepare(*mesh,mixed,rulesFor(*mesh,mixed,Mode::PureLiquidReservoir,Mode::PureGasReservoir));
    const auto mixedSlab=Slab::prepareBoundary(*mesh,mixed,mixedBoundary,end,a,b);
    const auto mixedTransport=Transport::prepareLocalPlicBoundarySweep(*mesh,mixed,mixedSlab,100,1);
    inspect(*mesh,mixed,mixedTransport);++scenes;
    check(mixedTransport.report().boundaryFull.outflow[0]>0&&mixedTransport.report().boundaryFull.outflow[1]>0,
        "mixed inside outflow must keep both independent phases");
    stage="F1 exact zero input phase";
    const auto gasBoundary=Boundary::prepare(*mesh,source,rulesFor(*mesh,source,Mode::PureGasReservoir,Mode::PureGasReservoir));
    const auto gasSlab=Slab::prepareBoundary(*mesh,source,gasBoundary,end,a,b);
    const auto gasTransport=Transport::prepareLocalPlicBoundarySweep(*mesh,source,gasSlab,100,1);
    inspect(*mesh,source,gasTransport);++scenes;
    for(const auto& q:gasTransport.endpoint().plicCells())check(q.inventory[0]==0,"zero liquid source and inflow must stay exact zero");
    stage="F1 AllClosed delegates original payload";
    auto closedRules=rules;for(auto& rule:closedRules)if(rule!=Mode::Interior)rule=Mode::Closed;
    const auto closed=Boundary::prepare(*mesh,source,closedRules);
    const auto ca=vofMaterialNormals(*mesh,source,0);
    const auto oldSlab=Slab::prepare(*mesh,source,end,ca,ca);
    const auto closedSlab=Slab::prepareBoundary(*mesh,source,closed,end,ca,ca);
    check(closed.allClosed()&&closedSlab.boundary()==nullptr,"AllClosed must delegate unchanged old slab");
    const auto oldRates=Transport::localPlicCourantRates(*mesh,source,ca);
    const auto closedRates=Transport::localPlicBoundaryCourantRates(*mesh,source,closed,ca);
    check(oldRates.maximumIncoming==closedRates.maximumIncoming&&oldRates.maximumOutflow==closedRates.maximumOutflow
        &&oldRates.maximumLocalColumn==closedRates.maximumLocalColumn,"AllClosed exact rates");
    const auto oldTransport=Transport::prepareLocalPlicSweep(*mesh,source,oldSlab,100,1);
    const auto closedTransport=Transport::prepareLocalPlicBoundarySweep(*mesh,source,closedSlab,100,1);
    for(std::size_t i=0;i<source.plicCells().size();++i){
        check(oldTransport.midpoint().plicCells()[i].inventory==closedTransport.midpoint().plicCells()[i].inventory
            &&oldTransport.endpoint().plicCells()[i].inventory==closedTransport.endpoint().plicCells()[i].inventory,
            "AllClosed exact raw independent stock pair");}
    for(std::size_t row=0;row<oldTransport.faces().size();++row)
        check(fields(oldTransport.faces()[row].half)==fields(closedTransport.faces()[row].half)
            &&fields(oldTransport.faces()[row].full)==fields(closedTransport.faces()[row].full),"AllClosed exact accepted flux payload");
    check(!closedTransport.report().boundaryHalf.enabled&&!closedTransport.report().boundaryFull.enabled,"AllClosed must keep old report path");
    ++scenes;
    stage="F1 RZ annular actual 2:1 physical throughflow";
    DictionaryMesh dictionary;dictionary.Nx=dictionary.Ny=4;dictionary.minCoord={1,0};dictionary.maxCoord={2,1};
    dictionary.maxCellNum=4096;dictionary.maxInitialIterations=0;dictionary.maxLevel_local={1};
    QuadtreeMesh annulus(dictionary);DictionaryMesh::RefineBox refine;refine.level=1;refine.min={1.5,0};refine.max={2,1};
    annulus.applyRefineBoxes({refine});
    const auto radial=Publication::prepareLocalPlic(annulus,domain(annulus,true),{},std::vector<double>(annulus.cellIndexNum(),0));
    const auto radialBoundary=Boundary::prepare(annulus,radial,rulesFor(annulus,radial,Mode::PureLiquidReservoir,Mode::PureGasReservoir));
    const auto ra=normals(annulus,radial,.125),rb=normals(annulus,radial,.25);
    const auto radialRates=Transport::localPlicBoundaryCourantRates(annulus,radial,radialBoundary,ra);
    check(radialRates.maximumIncoming>0&&radialRates.maximumLocalColumn>0,"RZ rate unavailable");
    const auto radialSlab=Slab::prepareBoundary(annulus,radial,radialBoundary,end,ra,rb);
    const auto radialTransport=Transport::prepareLocalPlicBoundarySweep(annulus,radial,radialSlab,100,1);
    inspect(annulus,radial,radialTransport);++scenes;
    near(radialTransport.report().boundaryFull.inflow[0],twoPi*dt*(.125L+.25L)/2,twoPi*dt,
        "RZ annulus true 2pi*r input measure",1e-12L);
    near(radialTransport.report().boundaryFull.outflow[1],twoPi*dt*(.125L+.25L)/2,twoPi*dt,
        "RZ annulus true 2pi*r output measure",1e-12L);
    stage="F1 true RZ axis cannot be reservoir or nonzero normal";
    const auto axisSource=Publication::prepareLocalPlic(*mesh,domain(*mesh,true),{},std::vector<double>(mesh->cellIndexNum(),0));
    auto axisRules=rulesFor(*mesh,axisSource,Mode::PureLiquidReservoir,Mode::PureGasReservoir);
    int axis=-1;for(int row=0;row<axisSource.domain().nativeFaceCount();++row)if(axisSource.domain().face(*mesh,row).area==0){axis=row;break;}
    check(axis>=0,"RZ witness lacks true zero-area axis");
    auto brokenAxis=axisRules;brokenAxis[axis]=Mode::PureLiquidReservoir;
    reject([&]{Boundary::prepare(*mesh,axisSource,brokenAxis);},"zero-area axis cannot use reservoir");
    const auto axisBoundary=Boundary::prepare(*mesh,axisSource,axisRules);
    auto axisNormal=normals(*mesh,axisSource,0);axisNormal[axis]=.125;
    reject([&]{Transport::localPlicBoundaryCourantRates(*mesh,axisSource,axisBoundary,axisNormal);},"true axis rate requires exact zero normal");
    reject([&]{Slab::prepareBoundary(*mesh,axisSource,axisBoundary,end,axisNormal,axisNormal);},"true axis endpoint requires exact zero normal");
    check(source.identity()==boundary.sourceIdentity()&&source.plicCells().data()==storage&&source.alpha()==oldAlpha,
        "open candidates/rejections changed accepted source owner/storage/alpha");
    for(std::size_t i=0;i<saved.size();++i)check(source.plicCells()[i].inventory==saved[i].inventory,"open candidate/rejection changed raw source phase stocks");
    std::cout<<std::setprecision(17)<<"LOCAL_PLIC_OPEN_BOUNDARY_SUMMARY,scenes="<<scenes
        <<",boundary_checks="<<boundaryFaces<<",negative_controls="<<negativeControls
        <<",assertions="<<assertions<<",INTERFACE_ONLY,PASS\n";
}

// F2: actual F1 open transport -> same material/P0 dual -> physical KEP.
// No pressure solve, App transaction, Electric field or new transport is built.
void localPlicOpenKepOnly()
{
    using Mode=Slab::Boundary::Mode;
    unsigned scenes=0,guardsChecked=0,tangentSegments=0,uRChecks=0,externalSegments=0;
    long double peakMomentum=0,peakPower=0,maximumActualR=0,maximumInletDissipation=0;
    const auto fail=[&](auto&& action,const char* reason){rejects(action,reason);++guardsChecked;};
    for(bool rz:{false,true}){
        stage=rz?"F2 true annular RZ 2:1":"F2 true XY 2:1";
        std::unique_ptr<QuadtreeMesh> mesh;
        if(!rz)mesh=vofMaterialMesh(4);
        else {
            DictionaryMesh dictionary;dictionary.Nx=dictionary.Ny=4;dictionary.minCoord={1,0};dictionary.maxCoord={2,1};
            dictionary.maxCellNum=4096;dictionary.maxInitialIterations=0;dictionary.maxLevel_local={1};
            mesh=std::make_unique<QuadtreeMesh>(dictionary);
            DictionaryMesh::RefineBox box;box.level=1;box.min={1.5,0};box.max={2,1};mesh->applyRefineBoxes({box});
        }
        const auto source=Publication::prepareLocalPlic(*mesh,domain(*mesh,rz),{},std::vector<double>(mesh->cellIndexNum(),0));
        const auto saved=source.plicCells();const auto alpha=source.alpha();const auto sourceId=source.identity();
        const auto* sourceStorage=source.plicCells().data();const auto& native=source.domain();
        check(!native.uniformLeaves(),"F2 fixture must retain actual2:1 leaves");
        std::vector<Mode> rules(native.nativeFaceCount(),Mode::Interior);
        std::vector<double> a(native.nativeFaceCount()),b(a.size());
        for(int row=0;row<native.nativeFaceCount();++row){const auto face=native.face(*mesh,row);
            if(!face.horizontal&&face.area>0){a[row]=rz?.125/face.center.x:.125;b[row]=rz?.25/face.center.x:.25;}
            if(face.negativeCell<0||face.positiveCell<0)
                rules[row]=face.horizontal||face.area==0?Mode::Closed:
                    (face.positiveCell>=0?Mode::PureLiquidReservoir:Mode::PureGasReservoir);
        }
        const auto boundary=Slab::Boundary::prepare(*mesh,source,rules);
        const double dt=1./64;const Clock end{1,0,dt,Position::End};
        const auto slab=Slab::prepareBoundary(*mesh,source,boundary,end,a,b);
        const auto volume=Transport::prepareLocalPlicBoundarySweep(*mesh,source,slab,100,1);
        check(volume.report().boundaryFull.enabled,"F2 did not consume actual F1 open ledger");
        const auto geometry=HydroMomentumGeometry::build(*mesh,vofMaterialMetrics(*mesh,source));
        const auto before=HydroPhaseGeometry::build(*mesh,geometry.metrics(),source);
        const auto middle=HydroPhaseGeometry::build(*mesh,geometry.metrics(),volume.midpoint());
        const auto after=HydroPhaseGeometry::build(*mesh,geometry.metrics(),volume.endpoint());
        const auto material=HydroMaterialTransport::prepareGeometricVof(*mesh,geometry,before,middle,after,volume);
        const auto m0=HydroMomentumMass::buildCellAverageProjection(*mesh,geometry,before,100,1);
        const auto mh=HydroMomentumMass::buildCellAverageProjection(*mesh,geometry,middle,100,1);
        const auto m1=HydroMomentumMass::buildCellAverageProjection(*mesh,geometry,after,100,1);
        const auto dual=HydroDualMaterial::prepare(*mesh,material,m0,mh,m1);
        ConservativeVofMomentumContract::qualify(*mesh,dual);
        const auto sources=dual.sources();const auto clock=dual.timeLayer();const auto& layout=geometry.layout();
        check(material.vofTransport().sameInstance(volume)&&dual.material().vofTransport().sameInstance(volume)
            && &dual.faceRates()==&material.dualFaceRates(),"F2 introduced a second material F owner");
        check(sameClock(clock,{1,0,dt,Position::Midpoint}),"F2 lost complete actual Midpoint interval");
        const auto ratesSaved=dual.faceRates(),integralsSaved=dual.faceIntegrals();const auto cellsSaved=dual.cells();
        std::vector<int> walls,freeOutflow;int incomingDof=-1;
        std::vector<double> u(layout.velocityCount());
        for(int i=0;i<layout.velocityCount();++i){const int row=layout.nativeFace(i);const auto face=native.face(*mesh,row);
            const bool first=layout.normalAxis(i)==HydroDofLayout::NormalAxis::First;
            u[i]=first?.5*a[row]+.5*b[row]:.375*face.center.y*(1-face.center.y);
            if(face.negativeCell<0||face.positiveCell<0){
                if(rules[row]==Mode::Closed){walls.push_back(i);u[i]=0;}
                else if(face.outwardSign*(.5*a[row]+.5*b[row])>0)freeOutflow.push_back(i);
                else incomingDof=i;
            }
        }
        check(!walls.empty()&&!freeOutflow.empty()&&incomingDof>=0,"F2 fixture lacks actual closed/inlet/outlet ports");
        const auto wallLaw=HydroAdvectorBoundary::impermeable(*mesh,geometry,walls);
        const HydroNativeValues values(layout,u,std::vector<double>(layout.pressureCount()));
        const auto& segments=geometry.segments();
        std::vector<Vector2D> begin(segments.size(),Vector2D(0,0)),finish(begin);
        int firstIn=-1,unused=-1;
        for(std::size_t id=0;id<segments.size();++id){const auto& g=segments[id];
            if(g.negativeVelocity>=0&&g.positiveVelocity>=0){unused=static_cast<int>(id);continue;}
            const int owner=g.negativeVelocity>=0?g.negativeVelocity:g.positiveVelocity;
            const auto q=(g.negativeVelocity>=0?1.L:-1.L)*dual.faceRates()[id];
            if(q>=0){unused=static_cast<int>(id);continue;}
            check(g.primalVelocity>=0,"F2 actual external segment has no native source");
            const int row=layout.nativeFace(g.primalVelocity);
            // The normal is exactly the real slab endpoint; tangent is explicit.
            begin[id]=g.normalAxis==0?Vector2D(a[row],.375):Vector2D(.375,a[row]);
            finish[id]=g.normalAxis==0?Vector2D(b[row],.625):Vector2D(.625,b[row]);
            const int component=layout.normalAxis(owner)==HydroDofLayout::NormalAxis::First?0:1;
            if(component!=g.normalAxis)++tangentSegments;
            firstIn=static_cast<int>(id);
        }
        check(firstIn>=0&&unused>=0,"F2 fixture lacks actual vector positive/negative witness");
        const auto inlet=HydroConvectionBoundary::bindVofEndpointVelocity(*mesh,dual,begin,finish,sources,clock);
        check(inlet.hasVofEndpointVelocity(),"F2 explicit inlet lost its endpoint full vectors");
        for(std::size_t i=0;i<begin.size();++i)check(inlet.beginVelocity()[i].x==begin[i].x&&inlet.beginVelocity()[i].y==begin[i].y
            &&inlet.endVelocity()[i].x==finish[i].x&&inlet.endVelocity()[i].y==finish[i].y,"F2 endpoint vector payload changed");
        std::vector<double> traces(segments.size());
        for(std::size_t id=0;id<segments.size();++id){const auto& g=segments[id];
            if(g.negativeVelocity>=0&&g.positiveVelocity>=0)continue;
            const int owner=g.negativeVelocity>=0?g.negativeVelocity:g.positiveVelocity;
            const auto q=(g.negativeVelocity>=0?1.L:-1.L)*dual.faceRates()[id];
            traces[id]=q<0?inlet.traces()[id]:q>0?u[owner]:0;
        }
        const auto appliedBoundary=HydroConvectionBoundary::bind(*mesh,dual,traces,sources,clock);
        // Independent raw-stock projection + actual signed dual flux scatter.
        // This keeps the true small R; it never assigns R=0 from primal closure.
        const auto impulse=dual.consistencyImpulse(*mesh,values,sources,clock);
        for(int i=0;i<layout.velocityCount();++i){
            VofMaterialWide change{},out{};long double scale=0;
            for(const auto& overlap:geometry.volumes()[i].overlaps){const int cell=layout.nativeCell(overlap.pressure);
                const auto weight=capacity({overlap.lower.x,overlap.upper.x,overlap.lower.y,overlap.upper.y},rz)/source.cells()[cell].capacity;
                const auto& old=source.plicCells()[cell].inventory;const auto& now=volume.endpoint().plicCells()[cell].inventory;
                auto local=vofMaterialAdd(vofMaterialMul({100,0},vofMaterialSub({now[0],0},{old[0],0})),
                    vofMaterialSub({now[1],0},{old[1],0}));
                change=vofMaterialAdd(change,vofMaterialMul({weight,0},local));
                scale+=weight*(100*std::abs(now[0]-old[0])+std::abs(now[1]-old[1]));
            }
            for(std::size_t g=0;g<segments.size();++g){const auto& segment=segments[g];
                if(segment.negativeVelocity==i){out=vofMaterialAdd(out,{dual.faceIntegrals()[g],0});scale+=std::abs(dual.faceIntegrals()[g]);}
                if(segment.positiveVelocity==i){out=vofMaterialSub(out,{dual.faceIntegrals()[g],0});scale+=std::abs(dual.faceIntegrals()[g]);}}
            if(scale==0){check(dual.cells()[i].massChange==0&&dual.cells()[i].integratedDefect==0,"F2 zero material row acquired R");}
            else {
                near(dual.cells()[i].massChange,vofMaterialValue(change),scale,"independent raw-stock P0 dual mass change",1e-12L);
                near(dual.cells()[i].integratedDefect,vofMaterialValue(vofMaterialAdd(change,out)),scale,"independent dual DeltaM+F R",1e-12L);
                near(dual.cells()[i].defectRate*dt,dual.cells()[i].integratedDefect,scale,"actual dual R/dt",1e-12L);
            }
            maximumActualR=std::max(maximumActualR,std::abs(dual.cells()[i].defectRate));
            const auto expected=vofMaterialValue(vofMaterialMul({u[i],0},{dual.cells()[i].integratedDefect,0}));
            if(expected==0)check(impulse[i]==0,"F2 zero uR acquired impulse");
            else {near(impulse[i],expected,std::abs(expected),"actual uR own-scale consistency impulse",1e-12L);++uRChecks;}
        }
        const auto verify=[&](const HydroKepConvection& op){
            check(op.isPrescribedFlow()&&!op.isImpermeable(),"F2 invoked a closed/legacy operator");
            const auto result=op.apply(*mesh,values,appliedBoundary,sources,clock);
            std::array<long double,2> boundaryMomentum{},physicalKE{},inletDiss{},boundaryWork{},momentumScale{},powerScale{};
            for(std::size_t id=0;id<segments.size();++id){const auto& g=segments[id];
                if(g.negativeVelocity>=0&&g.positiveVelocity>=0)continue;
                const int owner=g.negativeVelocity>=0?g.negativeVelocity:g.positiveVelocity;
                const int component=layout.normalAxis(owner)==HydroDofLayout::NormalAxis::First?0:1;
                const auto q=(g.negativeVelocity>=0?1.L:-1.L)*dual.faceRates()[id];
                const auto trace=static_cast<long double>(traces[id]),value=static_cast<long double>(u[owner]);
                boundaryMomentum[component]+=q*trace;
                physicalKE[component]+=.5L*q*trace*trace;
                inletDiss[component]+=q<0?.5L*std::abs(q)*(value-trace)*(value-trace):0;
                boundaryWork[component]+=value*q*trace-.5L*q*value*value;
                momentumScale[component]+=std::abs(q*trace);
                powerScale[component]+=.5L*std::abs(q)*(trace*trace+value*value)+(q<0?.5L*std::abs(q)*(value-trace)*(value-trace):0);
                ++externalSegments;
            }
            std::array<long double,2> momentum{},materialPower{},kinetic{};
            for(int i=0;i<layout.velocityCount();++i){
                const int c=layout.normalAxis(i)==HydroDofLayout::NormalAxis::First?0:1;const auto value=static_cast<long double>(u[i]);
                const auto r=dual.cells()[i].defectRate;
                const auto load=-result.momentumRates()[i];
                const auto scale=std::abs(result.momentumRates()[i]);
                check(result.momentumLoads()[i]==load,"F2 conservative applied load must be exactly -N without diagnostic uR");
                momentum[c]+=result.momentumRates()[i];materialPower[c]-=.5L*value*value*r;
                kinetic[c]+=value*result.momentumLoads()[i]-.5L*value*value*dual.cells()[i].massChange/dt;
                momentumScale[c]+=std::abs(result.momentumRates()[i]);
                powerScale[c]+=std::abs(value*result.momentumRates()[i])+std::abs(value*result.momentumLoads()[i])
                    +std::abs(.5L*value*value*r)+std::abs(.5L*value*value*dual.cells()[i].massChange/dt);
            }
            for(int c=0;c<2;++c){const auto& report=result.report().components[c];
                const auto ms=momentumScale[c],ps=powerScale[c];
                check(ms>0&&ps>0,"F2 must exercise both transported components");
                near(report.boundaryMomentumRate,boundaryMomentum[c],ms,"actual external momentum",1e-12L);
                near(momentum[c],boundaryMomentum[c],ms,"independent compatible momentum balance",1e-12L);
                near(report.boundaryPhysicalKineticPower,physicalKE[c],ps,"actual outward physical KE",1e-12L);
                near(report.inletDissipationPower,inletDiss[c],ps,"actual one-sided inlet upwind dissipation",1e-12L);
                near(report.boundaryKineticPower,boundaryWork[c],ps,"actual boundary kinetic work",1e-12L);
                near(boundaryWork[c],physicalKE[c]+inletDiss[c],ps,"physical KE plus inlet dissipation",1e-12L);
                near(report.materialKineticPower,materialPower[c],ps,"actual negative half-u-squared R power",1e-12L);
                near(kinetic[c]-materialPower[c]+physicalKE[c]+inletDiss[c],0,ps,"independent open conservative kinetic identity",1e-12L);
                check(report.momentumIdentity<=1e-12L&&report.kineticIdentity<=1e-12L&&report.inletDissipationPower>=0,
                    "F2 relaxed old momentum/kinetic gates");
                peakMomentum=std::max(peakMomentum,std::abs(momentum[c]-boundaryMomentum[c])/ms);
                peakPower=std::max(peakPower,std::abs(kinetic[c]-materialPower[c]+physicalKE[c]+inletDiss[c])/ps);
                maximumInletDissipation=std::max(maximumInletDissipation,inletDiss[c]);
            }
            // Same frozen affine load, without a second material flux or solve.
            const auto matrix=op.frozenMomentumOperator(*mesh,sources,clock);
            const auto load=op.frozenBoundaryMomentumLoad(*mesh,sources,clock);
            std::vector<long double> frozen(load.begin(),load.end()),frozenGross(load.size());
            for(std::size_t i=0;i<load.size();++i)frozenGross[i]=std::abs(load[i]);
            for(int column=0;column<matrix.outerSize();++column)
                for(HydroKineticHodge::Matrix::InnerIterator entry(matrix,column);entry;++entry){
                    const auto term=static_cast<long double>(entry.value())*u[entry.col()];
                    frozen[entry.row()]-=term;frozenGross[entry.row()]+=std::abs(term);}
            for(std::size_t i=0;i<frozen.size();++i){const auto scale=frozenGross[i]+std::abs(result.momentumLoads()[i]);
                if(scale==0)check(frozen[i]==result.momentumLoads()[i],"F2 zero frozen operator changed");
                else near(frozen[i],result.momentumLoads()[i],scale,"same physical frozen -H*u+B",1e-12L);}
            ++scenes;
        };
        const auto prescribed=HydroKepConvection::buildVofPrescribedFlow(*mesh,dual,wallLaw,inlet);
        const auto outlet=HydroKepConvection::buildVofPrescribedFlow(*mesh,dual,wallLaw,inlet,freeOutflow);
        verify(prescribed);verify(outlet);
        for(int dof:freeOutflow)check(std::find(outlet.fixedVelocityDofs().begin(),outlet.fixedVelocityDofs().end(),dof)==outlet.fixedVelocityDofs().end(),
            "F2 pressure-outlet normal remained fixed");
        auto wrongClock=clock;++wrongClock.step;
        fail([&]{HydroConvectionBoundary::bindVofEndpointVelocity(*mesh,dual,begin,finish,sources,wrongClock);},"F2 equal-time foreign fullclock");
        auto wrongSources=sources;wrongSources.vofIds[1]=source.identity();
        fail([&]{HydroConvectionBoundary::bindVofEndpointVelocity(*mesh,dual,begin,finish,wrongSources,clock);},"F2 foreign VOF source IDs");
        auto wrong=begin;const auto axis=segments[firstIn].normalAxis;
        if(axis==0)wrong[firstIn].x+=.03125;else wrong[firstIn].y+=.03125;
        fail([&]{HydroConvectionBoundary::bindVofEndpointVelocity(*mesh,dual,wrong,finish,sources,clock);},"F2 wrong actual F1 endpoint normal");
        wrong=begin;wrong[unused]=Vector2D(.125,.25);
        fail([&]{HydroConvectionBoundary::bindVofEndpointVelocity(*mesh,dual,wrong,finish,sources,clock);},"F2 nonzero unused internal/zero/outflow vector");
        wrong=begin;wrong[firstIn].x=std::numeric_limits<double>::quiet_NaN();
        fail([&]{HydroConvectionBoundary::bindVofEndpointVelocity(*mesh,dual,wrong,finish,sources,clock);},"F2 nonfinite actual full vector");
        const auto untyped=HydroConvectionBoundary::bind(*mesh,dual,inlet.traces(),sources,clock);
        fail([&]{HydroKepConvection::buildVofPrescribedFlow(*mesh,dual,wallLaw,untyped);},"F2 scalar trace cannot forge endpoint full-vector owner");
        fail([&]{HydroKepConvection::buildVofPrescribedFlow(*mesh,dual,wallLaw,inlet,{incomingDof});},"F2 free pressure outlet rejects actual backflow");
        fail([&]{outlet.apply(*mesh,values,appliedBoundary,sources,wrongClock);},"F2 apply foreign fullclock");
        const auto otherMaterial=HydroMaterialTransport::prepareGeometricVof(*mesh,geometry,before,middle,after,volume);
        const auto otherDual=HydroDualMaterial::prepare(*mesh,otherMaterial,m0,mh,m1);
        fail([&]{inlet.validateFor(*mesh,otherDual,otherDual.sources(),clock);},"F2 same-numeric foreign material owner");
        check(source.identity()==sourceId&&source.plicCells().data()==sourceStorage&&source.alpha()==alpha,
            "F2 prepared operators or rejection changed original VOF source");
        for(std::size_t i=0;i<saved.size();++i)check(source.plicCells()[i].inventory==saved[i].inventory,"F2 changed independent raw source stocks");
        check(dual.faceRates()==ratesSaved&&dual.faceIntegrals()==integralsSaved,"F2 changed the unique material F");
        for(std::size_t i=0;i<cellsSaved.size();++i)check(dual.cells()[i].integratedDefect==cellsSaved[i].integratedDefect
            &&dual.cells()[i].defectRate==cellsSaved[i].defectRate,"F2 repaired or replaced real R");
        std::cout<<std::setprecision(17)<<"LOCAL_PLIC_OPEN_KEP_SCENE,rz="<<rz<<",cells="<<mesh->cellIndexNum()
            <<",segments="<<segments.size()<<",freeOutflow="<<freeOutflow.size()<<",maxR="<<maximumActualR<<",PASS\n";
    }
    check(tangentSegments>0&&maximumInletDissipation>0,"F2 did not exercise explicit nonzero tangential inlet work");
    check(uRChecks>0&&maximumActualR>0,"F2 actual native material must exercise nonzero R/uR");
    std::cout<<std::setprecision(17)<<"LOCAL_PLIC_OPEN_KEP_SUMMARY,scenes="<<scenes<<",guards="<<guardsChecked
        <<",tangent_segments="<<tangentSegments<<",uR_checks="<<uRChecks<<",external_segments="<<externalSegments
        <<",momentum_peak="<<peakMomentum<<",kinetic_peak="<<peakPower<<",maximum_actual_R="<<maximumActualR
        <<",inlet_dissipation_max="<<maximumInletDissipation<<",NO_SADDLE_NO_APP,PASS\n";
}

// F3: explicit fixed-native-mesh open VOF -> SAME original A/C and saddle.
// This selector does not create an App, an open remap or a capillary boundary.
void localPlicOpenStepOnly()
{
    using Mode=Slab::Boundary::Mode;
    unsigned scenes=0,acceptedSteps=0,guardCount=0,endpointVectorChecks=0,phaseLedgerChecks=0,uRChecks=0;
    long double peakPhase=0,peakMomentum=0,peakEnergy=0,maximumReaction=0,maximumOutlet=0;
    const auto reject=[&](auto&& action,const char* label,const char* expected=nullptr){
        bool failed=false;std::string message;try{action();}catch(const std::exception& e){failed=true;message=e.what();}
        check(failed,std::string("F3 missing rejection: ")+label);
        if(expected)check(message.find(expected)!=std::string::npos,std::string("F3 wrong rejecting owner: ")+label+" "+message);
        ++guardCount;
    };
    for(bool rz:{false,true})for(bool pressureOutlet:{false,true}){
        stage=std::string("F3 ")+(rz?"RZ":"XY")+(pressureOutlet?" free pressure outlet":" prescribed normals");
        std::unique_ptr<QuadtreeMesh> mesh;
        if(!rz)mesh=vofMaterialMesh(4);
        else {DictionaryMesh cfg;cfg.Nx=cfg.Ny=4;cfg.minCoord={1,0};cfg.maxCoord={2,1};
            cfg.maxCellNum=4096;cfg.maxInitialIterations=0;cfg.maxLevel_local={1};mesh=std::make_unique<QuadtreeMesh>(cfg);
            DictionaryMesh::RefineBox box;box.level=1;box.min={1.5,0};box.max={2,1};mesh->applyRefineBoxes({box});}
        const auto initial=Publication::prepareLocalPlic(*mesh,domain(*mesh,rz),{},std::vector<double>(mesh->cellIndexNum(),0));
        const auto initialId=initial.identity();const auto originalStocks=initial.plicCells();
        const auto metrics=vofMaterialMetrics(*mesh,initial);const auto& layout=metrics.layout();const auto native=initial.domain();
        check(!native.uniformLeaves()&&mesh->cellIndexNum()==40,"F3 fixture must retain actual native 2:1 leaves");
        HydroVofPrescribedFlow flow;
        flow.rulesByNativeRow.assign(native.nativeFaceCount(),Mode::Interior);
        flow.normalVelocity.assign(layout.velocityCount(),0);
        std::vector<double> velocity(layout.velocityCount());
        int inlet=-1,outlet=-1,inside=-1;
        for(int u=0;u<layout.velocityCount();++u){
            const int row=layout.nativeFace(u);const auto f=native.face(*mesh,row);
            if(!f.horizontal&&f.area>0)velocity[u]=rz?.125/f.center.x:.125;
            if(f.outwardSign==0){inside=u;continue;}
            if(f.horizontal||f.area==0)flow.rulesByNativeRow[row]=Mode::Closed;
            else {flow.rulesByNativeRow[row]=f.outwardSign<0?Mode::PureLiquidReservoir:Mode::PureGasReservoir;
                if(f.outwardSign>0&&pressureOutlet){flow.pressureOutlets.push_back({u,.125});outlet=u;}
                else {flow.prescribedVelocityDofs.push_back(u);flow.normalVelocity[u]=velocity[u];
                    if(f.outwardSign<0)inlet=u;else outlet=u;}}
        }
        check(inlet>=0&&outlet>=0&&inside>=0,"F3 fixture lacks authentic inlet/outlet/interior");
        const auto* meshPtr=mesh.get();
        flow.endpointVelocity=[native,meshPtr,rz](int row,Clock endpoint){
            const auto f=native.face(*meshPtr,row);const double a=rz?.125/f.center.x:.125;
            return Vector2D(a,.03125*(1+endpoint.end));
        };
        const HydroNativeValues initialValues(layout,velocity,std::vector<double>(layout.pressureCount(),pressureOutlet?.125:0));
        const auto initialPhase=HydroPhaseGeometry::build(*mesh,metrics,initial);
        auto accepted=HydroTwoPhaseState::createVofPrescribedFlow(*mesh,initialPhase,initialValues,flow);
        const auto initialState=accepted;
        check(accepted.vofPrescribedFlow()&&!accepted.prescribedFlow(),"F3 crossed into old FirstJet or closed state");
        HydroVofAcConfig config;config.liquidDensity=2;config.gasDensity=1;
        config.momentumDensityRepresentation=HydroMomentumMass::DensityRepresentation::CellAverageProjection;
        config.surfaceTension=config.liquidViscosity=config.gasViscosity=0;
        constexpr double dt=1./16384;
        // Reject bad boundary declarations before any time advancement.
        auto badRules=flow;badRules.rulesByNativeRow[layout.nativeFace(inside)]=Mode::Closed;
        reject([&]{HydroTwoPhaseState::createVofPrescribedFlow(*mesh,initialPhase,initialValues,badRules);},"interior rule");
        auto badNormal=flow;badNormal.normalVelocity[inlet]*=2;
        reject([&]{HydroTwoPhaseState::createVofPrescribedFlow(*mesh,initialPhase,initialValues,badNormal);},"accepted steady normal mismatch");
        auto badVector=flow;const auto actualQuery=flow.endpointVelocity;
        badVector.endpointVelocity=[actualQuery](int row,Clock t){auto p=actualQuery(row,t);p.x*=2;return p;};
        reject([&]{HydroTwoPhaseState::createVofPrescribedFlow(*mesh,initialPhase,initialValues,badVector);},"source full-vector normal mismatch");
        // Source normal is valid; only the actual NEXT endpoint changes it.
        auto temporal=flow;
        temporal.endpointVelocity=[actualQuery](int row,Clock t){auto p=actualQuery(row,t);if(t.time()>0)p.x*=2;return p;};
        auto temporalState=HydroTwoPhaseState::createVofPrescribedFlow(*mesh,initialPhase,initialValues,temporal);
        const auto temporalSaved=temporalState;
        reject([&]{HydroVofStep::preparePrescribedFlow(*mesh,temporalState,{1,0,dt,Position::End},config);},
            "time-variable normal is not the declared steady law");
        check(temporalState.sameIdentity(temporalSaved),"F3 failed future-vector preparation changed accepted state");
        if(pressureOutlet){
            auto reversed=velocity;reversed[outlet]=-reversed[outlet];
            const HydroNativeValues backflow(layout,reversed,initialValues.pressure());
            reject([&]{HydroTwoPhaseState::createVofPrescribedFlow(*mesh,initialPhase,backflow,flow);},
                "pressure outlet backflow","pressure outlet reverses");
            auto fixedOutlet=flow;fixedOutlet.normalVelocity[outlet]=velocity[outlet];
            reject([&]{HydroTwoPhaseState::createVofPrescribedFlow(*mesh,initialPhase,initialValues,fixedOutlet);},
                "pressure outlet cannot also have prescribed payload");
        }
        for(int index=1;index<=2;++index){
            const auto old=accepted;const auto& source=old.phaseGeometry().vofPublication();const auto sourceId=source.identity();
            const auto sourceStocks=source.plicCells();const auto* sourceStorage=source.plicCells().data();
            const auto oldVelocity=old.values().normalVelocity(),oldPressure=old.values().pressure();
            const Clock next{index,(index-1)*dt,index*dt,Position::End};
            reject([&]{HydroVofStep::stabilityRates(*mesh,accepted,config);},"closed rate entry rejects open");
            reject([&]{HydroVofStep::prepare(*mesh,accepted,next,config);},"closed step entry rejects open");
            auto wrong=next;++wrong.step;
            reject([&]{HydroVofStep::preparePrescribedFlow(*mesh,accepted,wrong,config);},"noncontiguous full End clock");
            const auto rate=HydroVofStep::stabilityRatesPrescribedFlow(*mesh,accepted,config);
            std::vector<double> nativeNormal(native.nativeFaceCount());
            for(int u=0;u<layout.velocityCount();++u)nativeNormal[layout.nativeFace(u)]=accepted.values().normalVelocity()[u];
            const auto rateBoundary=Slab::Boundary::prepare(*mesh,source,flow.rulesByNativeRow);
            const auto directRate=Transport::localPlicBoundaryCourantRates(*mesh,source,rateBoundary,nativeNormal);
            const auto expectedRate=std::max({directRate.maximumIncoming,directRate.maximumOutflow,directRate.maximumLocalColumn});
            check(rate.flowRate==expectedRate&&rate.flowRate>0&&!rate.capillaryAvailable&&rate.capillaryRate==0,
                "F3 active rates did not consume actual F1 native accepted source");
            std::cout<<"LOCAL_PLIC_OPEN_STEP_PREPARE,rz="<<rz<<",pressureOutlet="<<pressureOutlet<<",step="<<index<<std::endl;
            const auto step=HydroVofStep::preparePrescribedFlow(*mesh,accepted,next,config);
            const auto& frozen=step.frozenMidpoint();const auto& dual=frozen.material();
            const auto& volume=step.transport();const auto& candidate=step.candidateState();
            const auto& um=step.solvedMidpointValues();const auto& u1=frozen.solveResult().values();
            const auto clock=dual.timeLayer();const auto sources=dual.sources();
            check(volume.sourceIdentity()==sourceId&&sameClock(volume.sourceTime(),old.timeLayer())
                &&sameClock(volume.endTime(),next)&&volume.advector().boundary(),"F3 lost actual accepted F1 boundary source/clock");
            check(dual.material().vofTransport().sameInstance(volume)&&frozen.convection().isPrescribedFlow()
                &&!frozen.hasInterfaceLoad()&&!frozen.hasMaxwellLoad()&&!step.hasCapillaryLoad(),
                "F3 original NS did not retain same unique actual F1/F2 owner");
            check(step.report().linearSolves>=1&&step.report().linearSolves<=2&&frozen.solveCount()==1,
                "F3 introduced an extra solver or bypassed original AC");
            check(step.report().usedVelocityHistory==(index==2)&&step.report().startupCorrection==(index==1),
                "F3 genuine predecessor velocity history was not retained");
            for(const auto& r:step.report().stages)check(r.saddle.algebraicRelativeResidual<=1e-11
                &&r.saddle.freeMomentumDefect<=1e-11&&r.momentumPowerRelativeDefect<=config.momentumPowerTolerance,
                "F3 relaxed original saddle/momentum-power gates");
            // Two independent stock sums and actual gross exterior phase fluxes.
            for(bool half:{true,false}){
                const auto& endpoint=half?volume.midpoint():volume.endpoint();
                const auto& ledger=half?volume.report().boundaryHalf:volume.report().boundaryFull;
                std::array<VofMaterialWide,2> change{};std::array<long double,2> in{},out{},now{},begin{};long double totalCapacity=0;
                for(std::size_t ci=0;ci<source.plicCells().size();++ci){totalCapacity+=source.cells()[ci].capacity;
                    for(unsigned p=0;p<2;++p){const auto a=source.plicCells()[ci].inventory[p],b=endpoint.plicCells()[ci].inventory[p];
                        change[p]=vofMaterialAdd(change[p],vofMaterialSub({b,0},{a,0}));begin[p]+=a;now[p]+=b;}}
                for(const auto& f:volume.faces()){
                    const auto& q=half?f.half:f.full;
                    const auto weighted=2*q.liquid+q.gas,faceMassScale=std::abs(q.mass)+2*std::abs(q.liquid)+std::abs(q.gas);
                    if(faceMassScale==0)check(q.mass==0&&weighted==0,"F3 zero independent phase flux acquired mass");
                    else near(q.mass,weighted,faceMassScale,"F3 SAME density-weighted accepted F",1e-12L);
                    if(f.negativeCell>=0&&f.positiveCell>=0)continue;
                    const auto sign=f.negativeCell>=0?1.L:-1.L;
                    for(unsigned p=0;p<2;++p){const auto positive=p==0?q.positiveLiquid:q.positiveGas;
                        const auto negative=p==0?q.negativeLiquid:q.negativeGas;
                        out[p]+=sign>0?positive:negative;in[p]+=sign>0?negative:positive;}
                    const auto face=native.face(*mesh,f.nativeRow);
                    if(face.outwardSign<0&&!face.horizontal){
                        const auto extent=half?dt/2:dt;const auto expected=static_cast<long double>(face.area)*extent*nativeNormal[f.nativeRow];
                        near(q.total,expected,std::abs(expected),"independent actual steady inlet volume",1e-12L);
                        check(q.liquid==q.total&&q.gas==0,"pure liquid reservoir was replaced by a donor/rounded alpha source");}
                }
                for(unsigned p=0;p<2;++p){const auto residual=vofMaterialValue(change[p])+out[p]-in[p];
                    near(residual,0,totalCapacity,"independent open DeltaQ+Out-In",1e-12L);
                    near(ledger.source[p],begin[p],totalCapacity,"independent phase ledger source",1e-12L);
                    near(ledger.target[p],now[p],totalCapacity,"independent phase ledger target",1e-12L);
                    near(ledger.inflow[p],in[p],totalCapacity,"independent phase ledger inflow",1e-12L);
                    near(ledger.outflow[p],out[p],totalCapacity,"independent phase ledger outflow",1e-12L);
                    peakPhase=std::max(peakPhase,std::abs(residual)/totalCapacity);++phaseLedgerChecks;}
                near(2*(vofMaterialValue(change[0])+out[0]-in[0])+vofMaterialValue(change[1])+out[1]-in[1],0,
                    3*totalCapacity,"independent physical boundary mass",1e-12L);
                check(ledger.enabled&&ledger.sourceIdentity==sourceId&&ledger.physicalRelativeConservation<=1e-12L
                    &&ledger.massRelativeConservation<=1e-12L&&ledger.algebraicRelativeConservation<=1e-12L,
                    "F3 changed original physical/WY boundary gates");
                check(in[0]>0&&in[1]==0&&out[0]==0&&out[1]>0,"F3 actual short front must inject liquid and emit original inside gas");
            }
            // Independent F2 incoming component from complete physical endpoint vectors.
            const auto& segments=dual.geometry().segments();std::vector<double> traces(segments.size());
            std::array<long double,2> boundaryMomentum{};long double physicalKE=0,inletDiss=0,boundaryWork=0,boundaryScale=0;
            for(std::size_t g=0;g<segments.size();++g){const auto& s=segments[g];
                if(s.negativeVelocity>=0&&s.positiveVelocity>=0)continue;
                const int owner=s.negativeVelocity>=0?s.negativeVelocity:s.positiveVelocity;
                const auto q=(s.negativeVelocity>=0?1.L:-1.L)*dual.faceRates()[g],v=static_cast<long double>(um.normalVelocity()[owner]);
                if(q<0){
                    check(s.primalVelocity>=0,"F3 incoming dual segment lacks actual native source");
                    const int row=layout.nativeFace(s.primalVelocity);const auto a=flow.endpointVelocity(row,old.timeLayer()),b=flow.endpointVelocity(row,next);
                    const int c=layout.normalAxis(owner)==HydroDofLayout::NormalAxis::First?0:1;
                    traces[g]=c==0?.5*a.x+.5*b.x:.5*a.y+.5*b.y;
                    check(traces[g]==frozen.appliedConvection().boundary().traces()[g],"F3 changed actual full-vector inlet component");
                    const auto normalA=s.normalAxis==0?a.x:a.y,normalB=s.normalAxis==0?b.x:b.y;
                    check(normalA==volume.advector().beginNormal()[row]&&normalB==volume.advector().endNormal()[row],
                        "F3 did not retain steady normal at both actual endpoint clocks");++endpointVectorChecks;
                }else if(q>0)traces[g]=um.normalVelocity()[owner];
                const auto trace=static_cast<long double>(traces[g]);
                const int c=layout.normalAxis(owner)==HydroDofLayout::NormalAxis::First?0:1;
                boundaryMomentum[c]+=q*trace;physicalKE+=.5L*q*trace*trace;
                inletDiss+=q<0?.5L*std::abs(q)*(v-trace)*(v-trace):0;
                boundaryWork+=v*q*trace-.5L*q*v*v;
                boundaryScale+=std::abs(q*trace)+std::abs(q)*(.5L*trace*trace+.5L*v*v)+
                    (q<0?.5L*std::abs(q)*(v-trace)*(v-trace):0);
            }
            check(boundaryScale>0&&inletDiss>0,"F3 did not exercise real tangential inlet dissipation");
            const auto cb=HydroConvectionBoundary::bind(*mesh,dual,traces,sources,clock);
            const auto convection=frozen.convection().apply(*mesh,um,cb,sources,clock);
            const auto& report=step.report().stages.back();
            for(int c=0;c<2;++c){near(report.boundaryMomentumRate[c],boundaryMomentum[c],boundaryScale,"F3 independent boundary momentum rate",1e-12L);
                check(convection.report().components[c].momentumIdentity<=1e-12L&&convection.report().components[c].kineticIdentity<=1e-12L,
                    "F3 changed F2 original momentum/KEP gates");}
            near(report.boundaryPhysicalKineticWork,dt*physicalKE,dt*boundaryScale,"F3 physical boundary KE work",1e-12L);
            near(report.inletDissipationWork,dt*inletDiss,dt*boundaryScale,"F3 independent inlet dissipation work",1e-12L);
            near(report.boundaryKineticWork,dt*boundaryWork,dt*boundaryScale,"F3 actual boundary work",1e-12L);
            near(boundaryWork,physicalKE+inletDiss,boundaryScale,"F3 boundary kinetic split",1e-12L);
            // Actual outlet traction is -p_out*n_out*Af, not a pressure-cell constraint.
            std::vector<long double> traction(layout.velocityCount());
            for(const auto& port:flow.pressureOutlets){const auto f=native.face(*mesh,layout.nativeFace(port.velocityDof));
                traction[port.velocityDof]=-static_cast<long double>(f.outwardSign)*f.area*port.pressure;
                check(f.outwardSign*u1.normalVelocity()[port.velocityDof]>0,"accepted pressure outlet reversed");
                near(frozen.saddle().pressureOutletLoads()[port.velocityDof],traction[port.velocityDof],
                    std::abs(traction[port.velocityDof]),"independent total-normal outlet traction",1e-12L);}
            check(frozen.saddle().pressureGaugeCount()==(pressureOutlet?0:1),"F3 imposed a gauge on an outlet-anchored component");
            const auto cp=frozen.pressureLoad().apply(*mesh,u1);
            std::vector<bool> free(layout.velocityCount());for(int u:frozen.saddle().freeVelocityDofs())free[u]=true;
            const Eigen::Map<const Eigen::VectorXd> endpointU(u1.normalVelocity().data(),layout.velocityCount());
            const Eigen::VectorXd au=frozen.velocityMatrix()*endpointU;
            long double equation=0,equationScale=0,physical=0,physicalScale=0,kinetic=0,pressureWork=0,outletWork=0,reactionWork=0,convectiveWork=0,massCross=0,energyScale=0;
            ConservativeVofMomentumContract::qualify(*mesh,dual);
            check(report.consistencyPower==0,"conservative actual Step applied a diagnostic uR power");
            const auto impulse=dual.consistencyImpulse(*mesh,um,sources,clock);
            for(int u=0;u<layout.velocityCount();++u){
                const long double v=um.normalVelocity()[u],a=old.values().normalVelocity()[u],b=u1.normalVelocity()[u];
                const long double m0=dual.beginMass().masses()[u],m1=dual.endMass().masses()[u];
                const long double rateP=(m1*b-m0*a)/dt,load=convection.momentumLoads()[u],residual=rateP+cp[u]-load-traction[u];
                check(load==-convection.momentumRates()[u],"F3 actual conservative Step load retained uR");
                const long double rowScale=std::abs(rateP)+std::abs(static_cast<long double>(cp[u]))+std::abs(load)+std::abs(traction[u]);
                if(free[u]){
                    physical=std::hypot(physical,residual);physicalScale=std::hypot(physicalScale,rowScale);
                    const auto rhs=frozen.integralMomentumRhs()[u]+frozen.saddle().pressureOutletLoads()[u];
                    equation=std::hypot(equation,static_cast<long double>(au[u])+cp[u]-rhs);
                    equationScale=std::hypot(equationScale,std::abs(static_cast<long double>(au[u]))+std::abs(static_cast<long double>(cp[u]))+std::abs(static_cast<long double>(rhs)));
                }else reactionWork+=dt*v*residual;
                const auto ku=.5L*(m1*b*b-m0*a*a),pw=dt*v*cp[u],ow=dt*v*traction[u],cw=dt*v*load,mc=.5L*(m1-m0)*a*b;
                kinetic+=ku;pressureWork+=pw;outletWork+=ow;convectiveWork+=cw;massCross+=mc;
                energyScale+=std::abs(ku)+std::abs(pw)+std::abs(ow)+std::abs(cw)+std::abs(mc)+dt*std::abs(v)*rowScale;
                const auto expectedImpulse=vofMaterialValue(vofMaterialMul({v,0},{dual.cells()[u].integratedDefect,0}));
                if(expectedImpulse==0)check(impulse[u]==0,"F3 zero uR acquired impulse");
                else near(impulse[u],expectedImpulse,std::abs(expectedImpulse),"F3 actual R/uR own-scale",1e-12L);
                ++uRChecks;
            }
            check(equationScale>0&&physicalScale>0&&energyScale>0,"F3 physical equation reference is inactive");
            near(equation,0,equationScale,"F3 original A+C free equation including outlet T once",1e-11L);
            near(physical,0,physicalScale,"F3 independent actual material/convection/pressure momentum",1e-11L);
            near(report.kineticChange,kinetic,energyScale,"F3 independent actual endpoint kinetic change",1e-12L);
            near(report.energy.pressureWork,pressureWork,energyScale,"F3 actual pressure work",1e-12L);
            near(report.pressureOutletWork,outletWork,energyScale,"F3 independent outlet pressure work",1e-12L);
            near(report.prescribedReactionWork,reactionWork,energyScale,"F3 fixed normal physical reaction work",1e-12L);
            const auto energyResidual=kinetic+pressureWork-convectiveWork-outletWork-reactionWork+massCross;
            near(energyResidual,0,energyScale,"F3 independent complete endpoint KE/mass-cross/boundary work",1e-11L);
            check(std::isfinite(report.energy.closureRelativeResidual)&&report.energy.closureRelativeResidual<=1e-11L,
                "F3 accepted energy ledger did not close");
            peakMomentum=std::max(peakMomentum,physical/physicalScale);peakEnergy=std::max(peakEnergy,std::abs(energyResidual)/energyScale);
            maximumReaction=std::max(maximumReaction,std::abs(reactionWork));maximumOutlet=std::max(maximumOutlet,std::abs(outletWork));
            auto wrongClock=clock;++wrongClock.step;
            reject([&]{frozen.validateFor(*mesh,sources,wrongClock);},"frozen fullclock provenance");
            auto wrongSources=sources;wrongSources.vofIds[1]=sourceId;
            reject([&]{frozen.validateFor(*mesh,wrongSources,clock);},"foreign midpoint VOF identity");
            reject([&]{candidate.validateVofFor(*mesh,sourceId,next);},"accepted candidate cannot reuse begin source identity");
            check(accepted.sameIdentity(old)&&source.plicCells().data()==sourceStorage&&accepted.values().normalVelocity()==oldVelocity
                &&accepted.values().pressure()==oldPressure,"F3 preparation/failure changed accepted source handle/U/p");
            for(std::size_t ci=0;ci<sourceStocks.size();++ci)check(source.plicCells()[ci].inventory==sourceStocks[ci].inventory
                &&source.plicCells()[ci].normal.x==sourceStocks[ci].normal.x&&source.plicCells()[ci].normal.y==sourceStocks[ci].normal.y
                &&source.plicCells()[ci].offset==sourceStocks[ci].offset,"F3 changed source stocks/normal/offset");
            check(candidate.vofPrescribedFlow()&&candidate.retainedMass()&&candidate.retainedMass()->masses().data()==dual.endMass().masses().data()
                &&candidate.values().normalVelocity()==u1.normalVelocity()&&sameClock(candidate.timeLayer(),next),
                "F3 candidate lost actual endpoint U/M/clock/flow");
            step.commit(*mesh,accepted);
            check(accepted.sameIdentity(candidate)&&accepted.velocityHistory()&&sameClock(accepted.timeLayer(),next),
                "F3 commit did not atomically adopt actual candidate and real predecessor history");
            reject([&]{step.commit(*mesh,accepted);},"stale repeated physical step commit");
            ++acceptedSteps;
            std::cout<<"LOCAL_PLIC_OPEN_STEP_ACCEPTED,rz="<<rz<<",pressureOutlet="<<pressureOutlet<<",step="<<index
                <<",actualT="<<next.end<<",solves="<<step.report().linearSolves<<",physical_momentum="<<physical/physicalScale
                <<",energy="<<std::abs(energyResidual)/energyScale<<",outletWork="<<outletWork<<",reactionWork="<<reactionWork<<std::endl;
        }
        // An authentic no-op Interface receipt reaches the explicit OPEN guard
        // before the closed fluid remapper's changed-tree work. No open AMR.
        const auto& finalPublication=accepted.phaseGeometry().vofPublication();
        const auto receipt=GeometricVofMeshTransfer::prepare(*mesh,finalPublication,*mesh,native,accepted.timeLayer());
        reject([&]{HydroVofMeshTransfer::prepare(*mesh,accepted,*mesh,receipt,config.liquidDensity,config.gasDensity);},
            "unsupported open remap","explicit physical boundary projection");
        check(initial.identity()==initialId&&initialState.phaseGeometry().vofIdentity()==initialId,
            "F3 two committed steps rewrote original publication/state");
        for(std::size_t ci=0;ci<originalStocks.size();++ci)check(initial.plicCells()[ci].inventory==originalStocks[ci].inventory,
            "F3 initial source was not immutable across accepted history");
        ++scenes;
    }
    check(scenes==4&&acceptedSteps==8&&endpointVectorChecks>0&&phaseLedgerChecks==32&&uRChecks>0&&maximumReaction>0&&maximumOutlet>0,
        "F3 required physical port/history/work evidence incomplete");
    std::cout<<"LOCAL_PLIC_OPEN_STEP_PASS,scenes="<<scenes<<",accepted_steps="<<acceptedSteps<<",guards="<<guardCount
        <<",endpoint_vectors="<<endpointVectorChecks<<",phase_ledger_checks="<<phaseLedgerChecks<<",uR_checks="<<uRChecks
        <<",phase_peak="<<peakPhase<<",momentum_peak="<<peakMomentum<<",energy_peak="<<peakEnergy
        <<",maximum_reaction_work="<<maximumReaction<<",maximum_outlet_work="<<maximumOutlet
        <<",FIXED_NATIVE_NO_APP_NO_OPEN_AMR_NO_HF_BOUNDARY_CLAIM\n";
}

// BEGIN F15 STRICT PURE-RESERVOIR TEMPORAL GROSS CONTRACT
namespace F15BoundaryGrossContract
{
using Exact=LAMP::Geometry::Detail::FiniteVolume::Dyadic;
using Boundary=Slab::Boundary;
using Mode=Boundary::Mode;
unsigned scenes=0,negativeControls=0,regressions=0,grossChecks=0,zeroNetFaces=0;
unsigned halfRootFaces=0,lateRootFaces=0,roundingFaces=0,incomingOtherWitnesses=0;
long double peakPhysical=0,peakAlgebraic=0,peakAnalytic=0;
std::unique_ptr<QuadtreeMesh> makeMesh(bool rz)
{
    (void)rz;
    DictionaryMesh d;d.Nx=d.Ny=4;d.maxCellNum=4096;d.maxInitialIterations=0;d.maxLevel_local={2};
    auto mesh=std::make_unique<QuadtreeMesh>(d);
    DictionaryMesh::RefineBox box;box.level=1;box.min={.5,.5};box.max={1,1};mesh->applyRefineBoxes({box});
    check(mesh->cellIndexNum()==28,"F15 actual asymmetric 2:1 fixture leaves");
    return mesh;
}
std::vector<double> normals(const QuadtreeMesh& mesh,const Publication& source,double value,unsigned axis=1)
{
    std::vector<double> result(source.domain().nativeFaceCount());
    for(int row=0;row<source.domain().nativeFaceCount();++row){const auto f=source.domain().face(mesh,row);
        if(f.horizontal==(axis==1)&&f.area>0)result[row]=value;}
    return result;
}
std::vector<Mode> rulesFor(const QuadtreeMesh& mesh,const Publication& source,unsigned phase,unsigned axis=1)
{
    check(phase<2,"F15 invalid test phase");std::vector<Mode> result(source.domain().nativeFaceCount(),Mode::Interior);
    for(int row=0;row<source.domain().nativeFaceCount();++row){const auto f=source.domain().face(mesh,row);
        if(f.negativeCell>=0&&f.positiveCell>=0)continue;
        result[row]=f.horizontal==(axis==1)&&f.area>0?(phase==0?Mode::PureLiquidReservoir:Mode::PureGasReservoir):Mode::Closed;}
    return result;
}
// Independent piecewise trapezoid of a linear endpoint law on t/dt in [0,s].
// This deliberately does not call the production FaceTime integrator.
std::array<long double,2> analyticGross(long double a,long double b,long double dt,long double area,long double s)
{
    std::array<long double,3> nodes{0,s,s};unsigned count=2;
    if(a*b<0){const auto root=-a/(b-a);if(root>0&&root<s){nodes[1]=root;nodes[2]=s;count=3;}}
    std::array<long double,2> result{};
    for(unsigned j=1;j<count;++j){const auto lo=nodes[j-1],hi=nodes[j];
        const auto first=a+(b-a)*lo,last=a+(b-a)*hi,mean=a+(b-a)*(lo+hi)/2;
        const auto amount=area*dt*(hi-lo)*(std::abs(first)+std::abs(last))/2;
        result[mean>=0?0:1]+=amount;}
    return result;
}
std::array<long double,12> fields(const GeometricVofIntegratedFlux& f)
{return {f.total,f.positiveTotal,f.negativeTotal,f.rawLiquid,f.liquid,f.gas,f.mass,
    f.positiveLiquid,f.negativeLiquid,f.positiveGas,f.negativeGas,f.limiter};}
bool sameFields(const GeometricVofIntegratedFlux& a,const GeometricVofIntegratedFlux& b)
{
    const auto x=fields(a),y=fields(b);
    for(unsigned i=0;i<x.size();++i)if(x[i]!=y[i]||std::signbit(x[i])!=std::signbit(y[i]))return false;
    return true;
}
template<class F> void rejectLiteral(F&& action,const char* literal)
{
    bool caught=false;try{action();}catch(const std::invalid_argument& error){
        caught=std::string(error.what()).find(literal)!=std::string::npos;
        check(caught,std::string("F15 wrong guard instead of ")+literal+": "+error.what());}
    check(caught,std::string("F15 missing specific guard: ")+literal);++negativeControls;
}
void inspect(const QuadtreeMesh& mesh,const Publication& source,const Transport& transport,unsigned phase,
    double a,double b,double dt,unsigned pattern)
{
    transport.validateFor(mesh,source,transport.advector(),100,1);
    const auto& report=transport.report();
    long double expectedIncoming=0,expectedOutgoing=0,expectedColumn=0;
    for(bool half:{true,false}){
        const auto& candidate=half?transport.midpoint():transport.endpoint();
        const auto& ledger=half?report.boundaryHalf:report.boundaryFull;
        check(ledger.enabled&&ledger.sourceIdentity==source.identity()&&sameClock(ledger.sourceTime,source.timeLayer())
            &&sameClock(ledger.targetTime,candidate.timeLayer()),"F15 actual boundary ledger owner/fullclock");
        std::array<Exact,2> before,after,in,out,rounding;
        Exact compression;
        std::vector<long double> incoming(source.plicCells().size()),outgoing(incoming.size());
        for(std::size_t i=0;i<source.plicCells().size();++i){
            check(candidate.plicCells()[i].inventory==source.plicCells()[i].inventory,
                "F15 matched pure transport must retain both stock values exactly");
            for(unsigned q=0;q<2;++q){before[q].add(source.plicCells()[i].inventory[q]);after[q].add(candidate.plicCells()[i].inventory[q]);}}
        for(const auto& face:transport.faces()){
            const auto native=source.domain().face(mesh,face.nativeRow);const auto& flux=half?face.half:face.full;
            check(flux.limiter==1,"F15 strict pure counterflow must not invoke a limiter");
            near(flux.mass,100*flux.liquid+flux.gas,std::max(std::abs(flux.mass),1e-30L),"F15 same phase-weighted mass flux",1e-12L);
            const bool outside=face.negativeCell<0||face.positiveCell<0;
            if(outside){
                const auto expected=native.horizontal&&native.area>0?analyticGross(a,b,dt,native.area,half?.5L:1.L):std::array<long double,2>{};
                const auto scale=std::max(expected[0]+expected[1],1e-30L);
                near(flux.positiveTotal,expected[0],scale,"F15 independent positive temporal area",1e-12L);
                near(flux.negativeTotal,expected[1],scale,"F15 independent negative temporal area",1e-12L);
                peakAnalytic=std::max(peakAnalytic,std::max(std::abs(flux.positiveTotal-expected[0]),std::abs(flux.negativeTotal-expected[1]))/scale);
                if(expected[0]+expected[1]>0){++grossChecks;
                    check(phase==0?(flux.positiveLiquid==flux.positiveTotal&&flux.negativeLiquid==flux.negativeTotal&&flux.positiveGas==0&&flux.negativeGas==0):
                        (flux.positiveGas==flux.positiveTotal&&flux.negativeGas==flux.negativeTotal&&flux.positiveLiquid==0&&flux.negativeLiquid==0),
                        "F15 gross phases must be the actual matching pure reservoir, not net or complementary recovery");
                    check((phase==0?flux.liquid:flux.gas)==flux.total&&(phase==0?flux.gas:flux.liquid)==0,"F15 canonical pure phase net");
                    if(!half&&pattern<2){check(flux.total==0&&flux.positiveTotal>0&&flux.negativeTotal>0,"F15 net zero must retain two actual nonzero gross branches");++zeroNetFaces;}
                    if(half&&pattern==2){check(flux.total==0&&flux.positiveTotal>0&&flux.negativeTotal>0,"F15 root before Midpoint must be represented in half extent");++halfRootFaces;}
                    if(half&&pattern==3){check(flux.positiveTotal==0&&flux.negativeTotal>0,"F15 root after Midpoint must not invent half positive flow");++lateRootFaces;}
                }
                const int sign=face.negativeCell>=0?1:-1;
                compression.add(flux.total,sign);
                for(unsigned q=0;q<2;++q){const auto positive=q==0?flux.positiveLiquid:flux.positiveGas;
                    const auto negative=q==0?flux.negativeLiquid:flux.negativeGas,stored=q==0?flux.liquid:flux.gas;
                    out[q].add(sign>0?positive:negative);in[q].add(sign>0?negative:positive);
                    Exact error(stored);error.add(positive,-1);error.add(negative);
                    if(error.compare(Exact{})!=0)++roundingFaces;
                    rounding[q].add(error,-sign);}
            }
            if(face.negativeCell>=0){outgoing[face.negativeCell]+=flux.positiveTotal;incoming[face.negativeCell]+=flux.negativeTotal;}
            if(face.positiveCell>=0){incoming[face.positiveCell]+=flux.positiveTotal;outgoing[face.positiveCell]+=flux.negativeTotal;}
            if(native.area>0&&native.horizontal){
                const int positiveDonor=face.negativeCell>=0?face.negativeCell:face.positiveCell;
                const int negativeDonor=face.positiveCell>=0?face.positiveCell:face.negativeCell;
                const auto& positiveBox=source.plicCells()[positiveDonor].rectangle;
                const auto& negativeBox=source.plicCells()[negativeDonor].rectangle;
                const auto positiveColumn=static_cast<long double>(native.area)*(positiveBox.upper.y-positiveBox.lower.y);
                const auto negativeColumn=static_cast<long double>(native.area)*(negativeBox.upper.y-negativeBox.lower.y);
                expectedColumn=std::max({expectedColumn,flux.positiveTotal/positiveColumn,flux.negativeTotal/negativeColumn});}
        }
        Exact capacitySum;for(const auto& cell:source.cells())capacitySum.add(cell.capacity);const auto scale=capacitySum.rounded();
        for(std::size_t i=0;i<incoming.size();++i){expectedIncoming=std::max(expectedIncoming,incoming[i]/source.cells()[i].capacity);
            expectedOutgoing=std::max(expectedOutgoing,outgoing[i]/source.cells()[i].capacity);}
        for(unsigned q=0;q<2;++q){
            check(ledger.source[q]==before[q].rounded()&&ledger.target[q]==after[q].rounded()
                &&ledger.inflow[q]==in[q].rounded()&&ledger.outflow[q]==out[q].rounded(),"F15 independently accumulated physical source/target/In/Out");
            check(ledger.netFluxRounding[q]==rounding[q].rounded(),"F15 signed stored-net minus exact-gross rounding receipt");
            check(ledger.compression[q]==(q==phase?compression.rounded():0)&&ledger.roundingError[q]==0,
                "F15 pure-stock WY source and final stock RNE receipts");
            auto physical=after[q];physical.add(before[q],-1);physical.add(out[q]);physical.add(in[q],-1);
            check(ledger.physicalResidual[q]==physical.rounded(),"F15 physical residual must use gross boundary flux");
            auto algebraic=physical;if(q==phase)algebraic.add(compression,-1);algebraic.add(rounding[q],-1);
            check(ledger.algebraicResidual[q]==algebraic.rounded(),"F15 exact signed gross/net/WY algebraic identity");
            check(physical.absolute().compare(Exact::product({1e-12L,scale}))<=0
                &&algebraic.absolute().compare(Exact::product({1e-12L,scale}))<=0,"F15 original physical and algebraic gates");
            peakPhysical=std::max(peakPhysical,std::abs(ledger.physicalResidual[q])/scale);
            peakAlgebraic=std::max(peakAlgebraic,std::abs(ledger.algebraicResidual[q])/scale);
        }
        check(ledger.physicalRelativeConservation<=1e-12L&&ledger.algebraicRelativeConservation<=1e-12L
            &&ledger.massRelativeConservation<=1e-12L,"F15 original physical/mass thresholds");
    }
    near(report.maximumIncomingCfl,expectedIncoming,std::max(expectedIncoming,1e-30L),"F15 real gross receiver Courant",1e-12L);
    near(report.maximumOutflowCfl,expectedOutgoing,std::max(expectedOutgoing,1e-30L),"F15 real gross donor Courant",1e-12L);
    near(report.maximumLocalColumnCfl,expectedColumn,std::max(expectedColumn,1e-30L),"F15 real gross subface-column Courant",1e-12L);
    check(expectedIncoming>0&&expectedOutgoing>0&&expectedColumn>0,"F15 net cancellation cannot erase gross CFL");
}
void checkSteadyPureGas(const QuadtreeMesh& mesh,const Publication& pure,double dt)
{
    const auto& native=pure.domain();const Clock end{1,0,dt,Position::End};
    const auto boundary=Boundary::prepare(mesh,pure,rulesFor(mesh,pure,1));
        const auto steady=normals(mesh,pure,.125);const auto steadySlab=Slab::prepareBoundary(mesh,pure,boundary,end,steady,steady);
        const auto steadyResult=Transport::prepareLocalPlicBoundarySweep(mesh,pure,steadySlab,100,1);
        for(const auto& f:steadyResult.faces()){
            const auto face=native.face(mesh,f.nativeRow);
            for(bool half:{true,false}){const auto& q=half?f.half:f.full;
                // Original F1 exterior authority is the existing Core temporal
                // event integrator; ideal Af*mean*duration can differ by 1 ULP.
                // Internal faces retain the original frozen-mean arithmetic.
                const auto canonical=LAMP::Geometry::PlanarFaceTime::integrateNormal(
                    steady[f.nativeRow],steady[f.nativeRow],end.begin,end.end,face.area);
                const bool exterior=f.negativeCell<0||f.positiveCell<0;
                const auto total=exterior?(half?canonical.halfPositive-canonical.halfNegative:canonical.fullPositive-canonical.fullNegative)
                    :static_cast<long double>(face.area)*.125L*(half?dt/2:dt)*(face.horizontal?1:0);
                const std::array<long double,12> expected{total,total,0,0,0,total,total,0,0,total,0,1};
                const auto actual=fields(q);
                const std::array<const char*,12> names{"total","positiveTotal","negativeTotal","rawLiquid","liquid","gas","mass",
                    "positiveLiquid","negativeLiquid","positiveGas","negativeGas","limiter"};
                for(unsigned field=0;field<actual.size();++field)if(actual[field]!=expected[field]){
                    std::cout<<std::setprecision(std::numeric_limits<long double>::max_digits10)
                        <<"F15_NONREVERSE_MISMATCH,rz="<<pure.domain().axisymmetric()<<",row="<<f.nativeRow
                        <<",negative="<<f.negativeCell<<",positive="<<f.positiveCell<<",horizontal="<<face.horizontal
                        <<",half="<<half<<",field="<<names[field]<<",actual="<<actual[field]<<",expected="<<expected[field]
                        <<std::hexfloat<<",actual_hex="<<actual[field]<<",expected_hex="<<expected[field]<<",native_area_hex="<<face.area
                        <<std::defaultfloat<<std::endl;}
                check(actual==expected,"F15 original nonreverse pure-gas accepted flux bits");}}
        check(steadyResult.report().boundaryHalf.netFluxRounding==std::array<long double,2>{}
            &&steadyResult.report().boundaryFull.netFluxRounding==std::array<long double,2>{},"F15 nonreverse rounding receipt remains exact zero");
}
void nonreverseOnly()
{
    for(bool rz:{false,true}){stage="F15 isolated exact nonreversing rz="+std::to_string(rz);
        auto mesh=makeMesh(rz);const auto pure=Publication::prepareLocalPlic(*mesh,domain(*mesh,rz),{},std::vector<double>(mesh->cellIndexNum(),0));
        checkSteadyPureGas(*mesh,pure,1./32);}
    std::cout<<"F15_NONREVERSE_COMPATIBILITY_PASS,XY_RZ_SAME_FIXTURE_ORIGINAL_EXACT_ASSERTION_NO_REVERSAL\n";
}
void run()
{
    const std::array<std::array<double,2>,4> patterns{{{.125,-.125},{-.125,.125},{.125,-.375},{-.375,.125}}};
    constexpr double dt=1./32;
    for(bool rz:{false,true}){
        auto mesh=makeMesh(rz);const auto native=domain(*mesh,rz);
        check(!native.uniformLeaves(),"F15 requires genuine 2:1 leaves");unsigned hanging=0,axisFaces=0;
        for(int row=0;row<native.nativeFaceCount();++row){const auto f=native.face(*mesh,row);
            if(f.area==0)++axisFaces;
            if(f.negativeCell>=0&&f.positiveCell>=0){const auto a=bounds(*mesh,f.negativeCell),b=bounds(*mesh,f.positiveCell);
                if(a.x1-a.x0!=b.x1-b.x0){const auto ratio=std::max(a.x1-a.x0,b.x1-b.x0)/std::min(a.x1-a.x0,b.x1-b.x0);check(ratio==2,"F15 actual 2:1 size ratio");++hanging;}}}
        check(hanging>0&&(!rz||axisFaces>0),"F15 missing actual hanging/true-axis witnesses");
        for(unsigned pattern=0;pattern<patterns.size();++pattern){
            stage="F15 strict reservoir reversal rz="+std::to_string(rz)+" pattern="+std::to_string(pattern);
            const unsigned phase=pattern%2;const auto source=Publication::prepareLocalPlic(*mesh,native,{},std::vector<double>(mesh->cellIndexNum(),phase==0?1:0));
            const auto saved=source.plicCells();const auto* storage=source.plicCells().data();const auto alpha=source.alpha();const auto id=source.identity();
            for(const auto& cell:source.plicCells())check(cell.resolved&&cell.inventory[1-phase]==0&&cell.measured[1-phase]==0&&cell.measured[phase]>0,
                "F15 producer must genuinely certify strict inside purity");
            const auto boundary=Boundary::prepare(*mesh,source,rulesFor(*mesh,source,phase));
            const auto a=normals(*mesh,source,patterns[pattern][0]),b=normals(*mesh,source,patterns[pattern][1]);
            const Clock end{pattern%2?2:1,pattern%2?.125:0,(pattern%2?.125:0)+dt,Position::End};
            // An accepted End source is producer-created for the even split order.
            const auto acceptedSource=pattern%2?Publication::prepareLocalPlic(*mesh,native,{1,0,.125,Position::End},alpha):source;
            const auto actualBoundary=pattern%2?Boundary::prepare(*mesh,acceptedSource,rulesFor(*mesh,acceptedSource,phase)):boundary;
            const auto slab=Slab::prepareBoundary(*mesh,acceptedSource,actualBoundary,end,a,b);
            const auto result=Transport::prepareLocalPlicBoundarySweep(*mesh,acceptedSource,slab,100,1);
            inspect(*mesh,acceptedSource,result,phase,patterns[pattern][0],patterns[pattern][1],dt,pattern);
            check(acceptedSource.identity()==slab.sourceIdentity()&&sameClock(acceptedSource.timeLayer(),slab.sourceTime()),"F15 retained actual accepted End/fullclock");
            if(pattern%2)for(std::size_t i=0;i<saved.size();++i)check(acceptedSource.plicCells()[i].inventory==saved[i].inventory,"F15 End clock preparation preserves source stock values");
            check(source.identity()==id&&source.plicCells().data()==storage&&source.alpha()==alpha,"F15 candidate changed original source owner/storage/alpha");
            for(std::size_t i=0;i<saved.size();++i)check(source.plicCells()[i].inventory==saved[i].inventory,"F15 source raw stocks changed");
            ++scenes;
        }
        stage="F15 exact nonreversing/AllClosed rz="+std::to_string(rz);
        const auto pure=Publication::prepareLocalPlic(*mesh,native,{},std::vector<double>(mesh->cellIndexNum(),0));
        const Clock end{1,0,dt,Position::End};const auto boundary=Boundary::prepare(*mesh,pure,rulesFor(*mesh,pure,1));
        checkSteadyPureGas(*mesh,pure,dt);++regressions;
        auto closedRules=rulesFor(*mesh,pure,1);for(auto& rule:closedRules)if(rule!=Mode::Interior)rule=Mode::Closed;
        const auto closed=Boundary::prepare(*mesh,pure,closedRules);const auto ca=vofMaterialNormals(*mesh,pure,0),cb=vofMaterialNormals(*mesh,pure,2);
        const auto oldSlab=Slab::prepare(*mesh,pure,end,ca,cb),delegatedSlab=Slab::prepareBoundary(*mesh,pure,closed,end,ca,cb);
        const auto oldResult=Transport::prepareLocalPlicSweep(*mesh,pure,oldSlab,100,1),delegated=Transport::prepareLocalPlicBoundarySweep(*mesh,pure,delegatedSlab,100,1);
        check(!delegatedSlab.boundary()&&!delegated.report().boundaryHalf.enabled&&!delegated.report().boundaryFull.enabled,"F15 AllClosed delegates original path");
        for(std::size_t i=0;i<pure.plicCells().size();++i)check(oldResult.midpoint().plicCells()[i].inventory==delegated.midpoint().plicCells()[i].inventory
            &&oldResult.endpoint().plicCells()[i].inventory==delegated.endpoint().plicCells()[i].inventory,"F15 old closed raw stock pair bits");
        for(std::size_t row=0;row<oldResult.faces().size();++row)check(sameFields(oldResult.faces()[row].half,delegated.faces()[row].half)
            &&sameFields(oldResult.faces()[row].full,delegated.faces()[row].full),"F15 old closed accepted flux field bits");++regressions;
        stage="F15 targeted purity/incoming-other/gross-CFL guards rz="+std::to_string(rz);
        const auto reverseA=normals(*mesh,pure,.125),reverseB=normals(*mesh,pure,-.125);
        const auto wrongPhase=Boundary::prepare(*mesh,pure,rulesFor(*mesh,pure,0));
        rejectLiteral([&]{const auto slab=Slab::prepareBoundary(*mesh,pure,wrongPhase,end,reverseA,reverseB);
            Transport::prepareLocalPlicBoundarySweep(*mesh,pure,slab,100,1);},"counterflow boundary requires a strict scratch-pure inside cell matching its reservoir");
        // F18 makes a producer-unresolved constant-alpha donor a legal P0
        // class. Keep this guard on a genuinely resolved mixed plane instead.
        constexpr long double mixedCut=.53125L;
        std::vector<double> resolvedAlpha(mesh->cellIndexNum());
        for(int ci=0;ci<mesh->cellIndexNum();++ci){const auto box=bounds(*mesh,ci);
            resolvedAlpha[ci]=box.x1<=mixedCut?1:box.x0>=mixedCut?0:static_cast<double>(rz?
                ((mixedCut-box.x0)*(mixedCut+box.x0))/((box.x1-box.x0)*(box.x1+box.x0)):
                (mixedCut-box.x0)/(box.x1-box.x0));}
        const auto mixed=Publication::prepareLocalPlic(*mesh,native,{},resolvedAlpha);
        auto mixedRules=rulesFor(*mesh,mixed,0);unsigned resolvedBoundaryCells=0;
        for(int row=0;row<native.nativeFaceCount();++row){const auto f=native.face(*mesh,row);
            if(!f.horizontal||(f.negativeCell>=0&&f.positiveCell>=0))continue;
            const int inside=f.negativeCell>=0?f.negativeCell:f.positiveCell;const auto& cell=mixed.plicCells()[inside];
            mixedRules[row]=cell.inventory[0]==0?Mode::PureGasReservoir:Mode::PureLiquidReservoir;
            if(cell.inventory[0]>0&&cell.inventory[1]>0){check(cell.resolved,
                "F18 migrated negative requires actual public-producer resolved mixed geometry");++resolvedBoundaryCells;}}
        check(resolvedBoundaryCells>0,"F18 resolved mixed exterior witness is required before rejection");
        const auto mixedBoundary=Boundary::prepare(*mesh,mixed,mixedRules);
        rejectLiteral([&]{const auto slab=Slab::prepareBoundary(*mesh,mixed,mixedBoundary,end,reverseA,reverseB);
            Transport::prepareLocalPlicBoundarySweep(*mesh,mixed,slab,100,1);},"counterflow boundary requires a strict scratch-pure inside cell matching its reservoir");
        std::vector<double> layeredAlpha(mesh->cellIndexNum(),1);
        for(int ci=0;ci<mesh->cellIndexNum();++ci)if(bounds(*mesh,ci).y1==1)layeredAlpha[ci]=0;
        const auto layered=Publication::prepareLocalPlic(*mesh,native,{},layeredAlpha);auto layeredRules=rulesFor(*mesh,layered,0);
        for(int row=0;row<native.nativeFaceCount();++row){const auto f=native.face(*mesh,row);
            if(f.horizontal&&f.negativeCell>=0&&f.positiveCell<0)layeredRules[row]=Mode::PureGasReservoir;
            if(f.horizontal&&f.negativeCell>=0&&f.positiveCell>=0&&layered.plicCells()[f.positiveCell].inventory[0]==0
                &&layered.plicCells()[f.negativeCell].inventory[1]==0){check(layered.plicCells()[f.negativeCell].inventory[0]>0
                    &&layered.plicCells()[f.positiveCell].inventory[1]>0,"F15 incoming-other witness phase sources");++incomingOtherWitnesses;}}
        const auto layeredSaved=layered.plicCells();const auto layeredBoundary=Boundary::prepare(*mesh,layered,layeredRules);
        const auto otherA=normals(*mesh,layered,-.125),otherB=normals(*mesh,layered,.375);
        rejectLiteral([&]{const auto slab=Slab::prepareBoundary(*mesh,layered,layeredBoundary,end,otherA,otherB);
            Transport::prepareLocalPlicBoundarySweep(*mesh,layered,slab,100,1);},"counterflow scratch-pure cell receives another phase in the actual direction");
        for(std::size_t ci=0;ci<layeredSaved.size();++ci)check(layered.plicCells()[ci].inventory==layeredSaved[ci].inventory,"F15 actual incoming-other rejection preserves source stocks");
        bool typed=false;try{const Clock large{1,0,4,Position::End};const auto slab=Slab::prepareBoundary(*mesh,pure,boundary,large,reverseA,reverseB);
            Transport::prepareLocalPlicBoundarySweep(*mesh,pure,slab,100,1);}
        catch(const GeometricVofCourantExceeded& error){typed=std::string(error.what()).find("Courant exceeded")!=std::string::npos;}
        check(typed,"F15 gross nonzero even when full net is zero must retain typed Courant rejection");++negativeControls;
    }
    check(scenes==8&&regressions==4&&negativeControls==8&&zeroNetFaces>0&&halfRootFaces>0&&lateRootFaces>0
        &&incomingOtherWitnesses>0&&roundingFaces>0,"F15 required distinct reversal/gross/rounding/phase witnesses incomplete");
    std::cout<<std::setprecision(17)<<"LOCAL_PLIC_BOUNDARY_GROSS_PASS,scenes="<<scenes<<",regressions="<<regressions
        <<",guards="<<negativeControls<<",gross_checks="<<grossChecks<<",zero_net_faces="<<zeroNetFaces
        <<",early_half_root_faces="<<halfRootFaces<<",late_root_half_faces="<<lateRootFaces
        <<",rounding_faces="<<roundingFaces<<",incoming_other_witnesses="<<incomingOtherWitnesses
        <<",physical_peak="<<peakPhysical<<",algebraic_peak="<<peakAlgebraic<<",analytic_peak="<<peakAnalytic
        <<",STRICT_PURE_RESERVOIR_INTERFACE_ONLY_NO_HYDRO_APP_CLAIM\n";
}
}
// END F15 STRICT PURE-RESERVOIR TEMPORAL GROSS CONTRACT

}

// BEGIN F16 EXPLICIT GROSS OPEN MOMENTUM CONTRACT
namespace F16OpenMomentumContract
{
using E=LAMP::Hydrodynamics::Detail::HydroFiniteContraction;
unsigned scenes=0,guards=0,netZeroActive=0,nonzeroAmbient=0,checkedRows=0;
long double peakMatrix=0,peakMomentum=0,peakPower=0,maximumRepresentation=0;
enum class ApplyRequest { Complete,WithoutCenter,MomentumLoads };
unsigned parallelPairs=0,literalFallbackPairs=0,omittedPhysicalPairs=0,omittedOverflowChecks=0;
unsigned momentumPhysicalPairs=0,momentumHugePositiveScenes=0;
unsigned productionPhysicalPairs=0,productionUnderflowWitnesses=0,productionZeroWitnesses=0,productionRouteTransitions=0,productionOwnerRebuilds=0;
unsigned productionClosedZeroWitnesses=0,productionClosedPositiveZeros=0,productionClosedNegativeZeros=0;
void equalOpenApply(const HydroKepConvectionResult& first,const HydroKepConvectionResult& second,bool compareCenter=true)
{
    check(first.hasComponentReport()&&second.hasComponentReport()&&first.hasPowerReport()&&second.hasPowerReport(),"F16 original full-power comparator received missing report");
    const auto equal=[](long double a,long double b){return a==b&&(a!=0||std::signbit(a)==std::signbit(b));};
    const auto rows=[&](const auto& a,const auto& b){check(a.size()==b.size(),"F16 serial/parallel dimension");
        for(std::size_t i=0;i<a.size();++i)check(equal(a[i],b[i]),"F16 serial/parallel original native value differs");};
    const auto centerState=[](const HydroKepConvectionResult& result){
        if(result.hasCenterMomentumRates())check(result.centerMomentumRates().size()==result.openApplyReport().segments,
            "F16 requested center result has wrong segment dimensions");
        else check(result.centerMomentumRates().empty(),"F16 unavailable center result must be empty, not fake zeros");};
    centerState(first);centerState(second);
    if(compareCenter){check(first.hasCenterMomentumRates()==second.hasCenterMomentumRates(),"F16 same-mode center availability differs");
        rows(first.centerMomentumRates(),second.centerMomentumRates());}
    else check(first.hasCenterMomentumRates()&&!second.hasCenterMomentumRates(),"F16 complete/omitted comparison has wrong availability");
    rows(first.momentumRates(),second.momentumRates());
    rows(first.momentumLoads(),second.momentumLoads());rows(first.dissipationRates(),second.dissipationRates());
    for(unsigned axis=0;axis<2;++axis){const auto values=[](const HydroKepComponentReport& r){return std::array<long double,17>{{
        r.momentumRate,r.boundaryMomentumRate,r.convectivePower,r.boundaryKineticPower,r.materialKineticPower,r.kineticRate,
        r.boundaryPhysicalKineticPower,r.inletDissipationPower,r.boundaryFluxRoundingPower,r.maximumBoundaryCanonicalGrossRateDefect,
        r.boundaryDiagonalPower,r.momentumIdentity,r.kineticIdentity,r.dissipationPower,r.maximumSensor,
        static_cast<long double>(r.dissipationConnections),static_cast<long double>(r.activeDissipationConnections)}};};
        rows(values(first.report().components[axis]),values(second.report().components[axis]));}
    const auto& a=first.openApplyReport();const auto& b=second.openApplyReport();
    check(a.available&&b.available&&a.vofSources==b.vofSources&&sameClock(a.timeLayer,b.timeLayer)
        &&a.nativeRows==b.nativeRows&&a.segments==b.segments&&a.pairs==b.pairs,"F16 actual apply receipt provenance changed");
    const std::array<double,9> times{{b.validationSeconds,b.preparationSeconds,b.segmentSeconds,b.diagonalSeconds,b.pairSeconds,
        b.rowSeconds,b.ledgerSeconds,b.fallbackSeconds,b.totalSeconds}};
    for(double value:times)check(std::isfinite(value)&&value>=0,"F16 passive actual apply time is invalid");
}
// The minimal request retains actual momentum/eta and known sensor-free zeros.
// It must never masquerade as a computed zero power report.
void checkMomentumOnlyReport(const HydroKepConvectionResult& result)
{
    check(result.hasComponentReport()&&!result.hasPowerReport()&&!result.hasCenterMomentumRates()&&result.centerMomentumRates().empty(),
        "F16 momentum request exposed unavailable center/power as available");
    for(const auto& c:result.report().components){
        const std::array<long double,9> missing{{c.convectivePower,c.boundaryKineticPower,c.materialKineticPower,c.kineticRate,
            c.boundaryPhysicalKineticPower,c.inletDissipationPower,c.boundaryFluxRoundingPower,c.boundaryDiagonalPower,c.kineticIdentity}};
        for(auto value:missing)check(std::isnan(value),"F16 unrequested power field is a fake finite value");
        for(auto value:{c.momentumRate,c.boundaryMomentumRate,c.momentumIdentity,c.maximumBoundaryCanonicalGrossRateDefect})
            check(std::isfinite(value),"F16 requested momentum/eta result is missing or nonfinite");
        check(c.momentumIdentity>=0&&c.momentumIdentity<=1e-12L,"F16 minimal request omitted original momentum gate");
        check(c.dissipationPower==0&&c.maximumSensor==0&&c.dissipationConnections==0&&c.activeDissipationConnections==0,
            "F16 known sensor-free stabilization values changed");
    }
}
void equalMomentumApply(const HydroKepConvectionResult& first,const HydroKepConvectionResult& second)
{
    checkMomentumOnlyReport(second);if(!first.hasPowerReport())checkMomentumOnlyReport(first);
    const auto equal=[](long double a,long double b){return a==b&&(a!=0||std::signbit(a)==std::signbit(b));};
    const auto rows=[&](const auto& a,const auto& b){check(a.size()==b.size(),"F16 minimal momentum dimension differs");
        for(std::size_t i=0;i<a.size();++i)check(equal(a[i],b[i]),"F16 minimal momentum original native value differs");};
    rows(first.momentumRates(),second.momentumRates());rows(first.momentumLoads(),second.momentumLoads());
    rows(first.dissipationRates(),second.dissipationRates());
    for(unsigned axis=0;axis<2;++axis){const auto values=[](const HydroKepComponentReport& c){return std::array<long double,8>{{
        c.momentumRate,c.boundaryMomentumRate,c.momentumIdentity,c.maximumBoundaryCanonicalGrossRateDefect,
        c.dissipationPower,c.maximumSensor,static_cast<long double>(c.dissipationConnections),static_cast<long double>(c.activeDissipationConnections)}};};
        rows(values(first.report().components[axis]),values(second.report().components[axis]));}
    const auto& a=first.openApplyReport();const auto& b=second.openApplyReport();
    check(a.available&&b.available&&a.vofSources==b.vofSources&&sameClock(a.timeLayer,b.timeLayer)
        &&a.nativeRows==b.nativeRows&&a.segments==b.segments&&a.pairs==b.pairs,"F16 minimal apply receipt provenance changed");
    const std::array<double,9> times{{b.validationSeconds,b.preparationSeconds,b.segmentSeconds,b.diagonalSeconds,b.pairSeconds,
        b.rowSeconds,b.ledgerSeconds,b.fallbackSeconds,b.totalSeconds}};
    for(double value:times)check(std::isfinite(value)&&value>=0,"F16 minimal actual stage time is invalid");
}
HydroKepConvectionResult pairedApply(const QuadtreeMesh& mesh,const HydroKepConvection& op,const HydroNativeValues& values,
    HydroMaterialSources sources,LAMP::TimeLayer clock,ApplyRequest request=ApplyRequest::Complete)
{
    const bool includeCenterDiagnostic=request==ApplyRequest::Complete;
    const bool includePowerReport=request!=ApplyRequest::MomentumLoads;
    const auto apply=[&](LAMP::TimeLayer actualClock){return includeCenterDiagnostic
        ?op.applyVofOpen(mesh,values,sources,actualClock)
        :includePowerReport?op.applyVofOpenWithoutCenterDiagnostic(mesh,values,sources,actualClock)
        :op.applyVofOpenMomentumLoads(mesh,values,sources,actualClock);};
    const auto equalRequested=[&](const HydroKepConvectionResult& a,const HydroKepConvectionResult& b){
        if(includePowerReport)equalOpenApply(a,b);else equalMomentumApply(a,b);};
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
    struct Runtime {int threads=omp_get_max_threads(),dynamic=omp_get_dynamic(),error=errno;std::fenv_t env{};unsigned csr=_mm_getcsr();
        Runtime(){check(std::fegetenv(&env)==0,"F16 paired original FP capture");}
        ~Runtime(){omp_set_dynamic(dynamic);omp_set_num_threads(threads);std::fesetenv(&env);_mm_setcsr(csr);errno=error;}} restore;
    omp_set_dynamic(0);
    std::fenv_t entry{};check(std::feraiseexcept(FE_INEXACT)==0&&std::fegetenv(&entry)==0,"F16 prime public FP receipt");const unsigned csr=_mm_getcsr();
    const auto reset=[&]{check(std::fesetenv(&entry)==0,"F16 common FP restore");_mm_setcsr(csr);errno=173;};
    struct Final {std::fenv_t env{};unsigned csr=0;int error=0;
        void capture(){error=errno;csr=_mm_getcsr();check(std::fegetenv(&env)==0,"F16 completed actual FP capture");}
        bool same(const Final& other)const{return error==other.error&&csr==other.csr&&std::memcmp(&env,&other.env,sizeof env)==0;}} sf,pf;
    omp_set_num_threads(1);reset();const auto serial=apply(clock);sf.capture();
    check(serial.hasCenterMomentumRates()==includeCenterDiagnostic,"F16 literal result ignored center request");
    check(serial.hasPowerReport()==includePowerReport,"F16 literal result ignored power request");
    check(!serial.openApplyReport().parallelUsed,"F16 one-worker apply used optional route");
    omp_set_num_threads(2);reset();const auto actual=apply(clock);pf.capture();
    equalRequested(serial,actual);check(sf.same(pf),"F16 actual serial/parallel raw native FP/CSR/errno differs");
    const auto& r=actual.openApplyReport();
    check(r.parallelRequested&&r.parallelEligible&&r.parallelUsed&&!r.parallelFallback&&r.contributionWorkers>=2
        &&r.gatherWorkers>=2&&r.scratchBytes>0&&r.scratchBytes<=64u*1024u*1024u,"F16 actual contribution/gather was not consumed");++parallelPairs;
    // The same six physical scenes (and their changed-U calls) also consume
    // an eight-worker request. Observe any real runtime team limit separately.
    unsigned expectedWorkers=0;omp_set_num_threads(8);
#pragma omp parallel num_threads(8)
    {
#pragma omp master
        {expectedWorkers=static_cast<unsigned>(omp_get_num_threads());}
    }
    reset();const auto eight=apply(clock);Final ef;ef.capture();
    equalRequested(serial,eight);check(sf.same(ef),"F04 actual eight-worker full FP/CSR/errno differs from literal");
    const auto& e=eight.openApplyReport();check(e.parallelRequested&&e.parallelEligible&&e.parallelUsed&&!e.parallelFallback
        &&expectedWorkers>=2&&e.contributionWorkers==expectedWorkers&&e.gatherWorkers==expectedWorkers
        &&e.scratchBytes>0&&e.scratchBytes<=64u*1024u*1024u,"F04 requested/actual team partial route was not consumed");
    // A real nested caller is ineligible and must use the literal operation.
    std::optional<HydroKepConvectionResult> nested;std::exception_ptr error;
#pragma omp parallel num_threads(2)
    {
#pragma omp single
        {try{nested.emplace(apply(clock));}catch(...){error=std::current_exception();}}
    }
    if(error)std::rethrow_exception(error);check(nested.has_value(),"F16 nested literal result missing");equalRequested(serial,*nested);
    check(!nested->openApplyReport().parallelUsed&&!nested->openApplyReport().parallelEligible,"F16 nested caller escaped literal fallback");++literalFallbackPairs;
    // Both actual calls retain the same original first validation error.
    auto stale=clock;stale.begin=.001;std::string first,second;
    for(int workers:{1,2}){omp_set_num_threads(workers);reset();try{(void)apply(stale);}
        catch(const std::invalid_argument& e){(workers==1?first:second)=e.what();}}
    check(!first.empty()&&first==second,"F16 serial/parallel first full-clock guard differs");
    return actual;
#else
    const auto actual=apply(clock);check(actual.openApplyReport().available,"F16 literal actual receipt missing");
    check(actual.hasCenterMomentumRates()==includeCenterDiagnostic,"F16 literal result ignored center request");
    check(actual.hasPowerReport()==includePowerReport,"F16 literal result ignored power request");return actual;
#endif
}
// Compare the actual typed producer against the retained full binding path.
// Same-batch support is shared; immutable-v40 Main remains the cross-version oracle.
void equalOpenBinding(const HydroConvectionBoundary& reference,const HydroConvectionBoundary& actual)
{
    check(reference.hasVofOpenMomentum()&&actual.hasVofOpenMomentum()&&reference.hasVofEndpointVelocity()&&actual.hasVofEndpointVelocity(),
        "F16 production binder lost typed gross or full endpoint payload");
    check(&reference.material().cells()==&actual.material().cells()&&reference.freePressureOutletVelocityDofs()==actual.freePressureOutletVelocityDofs(),
        "F16 production binder changed genuine material owner or ordered ports");
    const auto equal=[](long double a,long double b){return a==b&&(a!=0||std::signbit(a)==std::signbit(b));};
    const auto rows=[&](const auto& a,const auto& b){check(a.size()==b.size(),"F16 production binding dimensions differ");
        for(std::size_t i=0;i<a.size();++i)check(equal(a[i],b[i]),"F16 production binding scalar or zero sign differs");};
    rows(reference.traces(),actual.traces());
    const auto& a=reference.vofOpenMomentumFluxes();const auto& b=actual.vofOpenMomentumFluxes();check(a.size()==b.size(),"F16 production gross table dimensions differ");
    for(std::size_t i=0;i<a.size();++i){check(a[i].nativeRow==b[i].nativeRow&&a[i].velocityDof==b[i].velocityDof&&a[i].incomingModel==b[i].incomingModel,
            "F16 production binding changed gross row/orientation/model");
        const auto values=[](const auto& f){return std::array<long double,8>{{f.outgoingPhaseIntegral[0],f.outgoingPhaseIntegral[1],
            f.incomingPhaseIntegral[0],f.incomingPhaseIntegral[1],f.canonicalOutwardMassRate,f.outgoingMassRate,f.incomingMassRate,f.incomingMomentum}};};
        rows(values(a[i]),values(b[i]));
        const auto& r0=reference.beginVelocity()[i];const auto& r1=reference.endVelocity()[i];
        const auto& p0=actual.beginVelocity()[i];const auto& p1=actual.endVelocity()[i];
        rows(std::array<long double,4>{{r0.x,r0.y,r1.x,r1.y}},std::array<long double,4>{{p0.x,p0.y,p1.x,p1.y}});
    }
}
void equalPortCertificate(const HydroKepConvection& reference,const HydroKepConvection& production)
{
    const auto a=reference.momentReport().maxNativeMeanDefect,b=production.momentReport().maxNativeMeanDefect;
    check(a==b&&(a!=0||std::signbit(a)==std::signbit(b)),"F16 actual atom-partition/native-mean certificate differs from retained full certification");
}
// Production contracts compare physical values with the retained complete
// reference, but compare native FP/first-error only within the same request.
void checkProductionReport(const HydroKepConvectionResult& result,bool boundaryBudget)
{
    check(!result.hasComponentReport()&&!result.hasPowerReport()&&!result.hasCenterMomentumRates()&&result.centerMomentumRates().empty(),
        "F16 production result claimed unrequested complete diagnostics");
    check((result.openBoundaryBudget()!=nullptr)==boundaryBudget,"F16 production boundary capability differs from request");
    for(const auto& c:result.report().components){
        const std::array<long double,13> missing{{c.momentumRate,c.boundaryMomentumRate,c.convectivePower,c.boundaryKineticPower,
            c.materialKineticPower,c.kineticRate,c.boundaryPhysicalKineticPower,c.inletDissipationPower,c.boundaryFluxRoundingPower,
            c.maximumBoundaryCanonicalGrossRateDefect,c.boundaryDiagonalPower,c.momentumIdentity,c.kineticIdentity}};
        for(auto value:missing)check(std::isnan(value),"F16 production missing component field is fake finite data");
        check(c.dissipationPower==0&&c.maximumSensor==0&&c.dissipationConnections==0&&c.activeDissipationConnections==0,
            "F16 production sensor-free stabilization changed");
    }
}
void equalProductionApply(const HydroKepConvectionResult& reference,const HydroKepConvectionResult& actual,bool boundaryBudget,bool sameSource=true)
{
    checkProductionReport(actual,boundaryBudget);if(!reference.hasComponentReport())checkProductionReport(reference,boundaryBudget);
    const auto equal=[](long double a,long double b){return a==b&&(a!=0||std::signbit(a)==std::signbit(b));};
    const auto rows=[&](const auto& a,const auto& b){check(a.size()==b.size(),"F16 production native dimensions differ");
        for(std::size_t i=0;i<a.size();++i)check(equal(a[i],b[i]),"F16 production native value or zero sign differs");};
    rows(reference.momentumRates(),actual.momentumRates());rows(reference.momentumLoads(),actual.momentumLoads());rows(reference.dissipationRates(),actual.dissipationRates());
    if(boundaryBudget)for(unsigned axis=0;axis<2;++axis){
        const auto values=[](const auto& c){return std::array<long double,6>{{c.boundaryMomentumRate,c.boundaryKineticPower,
            c.boundaryPhysicalKineticPower,c.inletDissipationPower,c.boundaryFluxRoundingPower,c.boundaryDiagonalPower}};};
        const auto expected=reference.hasComponentReport()?values(reference.report().components[axis]):values(reference.openBoundaryBudget()->components[axis]);
        rows(expected,values(actual.openBoundaryBudget()->components[axis]));
    }
    const auto& a=reference.openApplyReport();const auto& b=actual.openApplyReport();
    check(a.available&&b.available&&(!sameSource||a.vofSources==b.vofSources)&&sameClock(a.timeLayer,b.timeLayer)
        &&a.nativeRows==b.nativeRows&&a.segments==b.segments&&a.pairs==b.pairs,"F16 production receipt lost actual source/clock/dimensions");
    for(double value:{b.validationSeconds,b.preparationSeconds,b.segmentSeconds,b.diagonalSeconds,b.pairSeconds,b.rowSeconds,b.ledgerSeconds,b.fallbackSeconds,b.totalSeconds})
        check(std::isfinite(value)&&value>=0,"F16 production passive time invalid");
    check(b.gatherWorkers==0&&b.diagonalSeconds==0&&b.pairSeconds==0&&b.ledgerSeconds==0,
        "F16 fused production receipt invented removed gather/stages");
    if(!boundaryBudget)check(b.segmentSeconds==0,"F16 prediction requested boundary power work");
}
HydroKepConvectionResult pairedProductionApply(const QuadtreeMesh& mesh,const HydroKepConvection& op,const HydroNativeValues& values,
    HydroMaterialSources sources,LAMP::TimeLayer clock,bool boundaryBudget)
{
    const auto apply=[&](const HydroNativeValues& state,HydroMaterialSources actualSources,LAMP::TimeLayer actualClock){return boundaryBudget
        ?op.applyVofOpenBoundaryBudget(mesh,state,actualSources,actualClock):op.applyVofOpenProductionLoads(mesh,state,actualSources,actualClock);};
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
    struct Runtime {int threads=omp_get_max_threads(),dynamic=omp_get_dynamic(),error=errno;std::fenv_t env{};unsigned csr=_mm_getcsr();
        Runtime(){check(std::fegetenv(&env)==0,"F16 production runtime FP capture");}
        ~Runtime(){omp_set_dynamic(dynamic);omp_set_num_threads(threads);std::fesetenv(&env);_mm_setcsr(csr);errno=error;}} restore;
    omp_set_dynamic(0);std::fenv_t entry{};check(std::feraiseexcept(FE_INEXACT)==0&&std::fegetenv(&entry)==0,"F16 production prime FP");const unsigned csr=_mm_getcsr();
    const auto reset=[&]{check(std::fesetenv(&entry)==0,"F16 production common FP restore");_mm_setcsr(csr);errno=173;};
    struct Final {std::fenv_t env{};unsigned csr=0;int error=0;
        void capture(){error=errno;csr=_mm_getcsr();check(std::fegetenv(&env)==0,"F16 production final FP capture");}
        bool same(const Final& other)const{return error==other.error&&csr==other.csr&&std::memcmp(&env,&other.env,sizeof env)==0;}} sf,pf;
    omp_set_num_threads(1);reset();const auto serial=apply(values,sources,clock);sf.capture();checkProductionReport(serial,boundaryBudget);
    check(!serial.openApplyReport().parallelUsed,"F16 production one-worker path used optional team");
    omp_set_num_threads(2);reset();const auto actual=apply(values,sources,clock);pf.capture();equalProductionApply(serial,actual,boundaryBudget);
    check(sf.same(pf),"F16 same production request serial/parallel FP/CSR/errno differs");
    const auto& r=actual.openApplyReport();check(r.parallelRequested&&r.parallelEligible&&r.parallelUsed&&!r.parallelFallback
        &&r.contributionWorkers>=2&&r.gatherWorkers==0&&r.scratchBytes>0&&r.scratchBytes<=64u*1024u*1024u,
        "F16 production native-row actual team was not consumed");
    unsigned expectedWorkers=0;omp_set_num_threads(8);
#pragma omp parallel num_threads(8)
    {
#pragma omp master
        {expectedWorkers=static_cast<unsigned>(omp_get_num_threads());}
    }
    reset();const auto eight=apply(values,sources,clock);Final ef;ef.capture();equalProductionApply(serial,eight,boundaryBudget);
    check(sf.same(ef),"F16 production eight-worker own FP/CSR/errno differs");
    const auto& e=eight.openApplyReport();check(e.parallelRequested&&e.parallelEligible&&e.parallelUsed&&!e.parallelFallback
        &&expectedWorkers>=2&&e.contributionWorkers==expectedWorkers&&e.gatherWorkers==0&&e.scratchBytes>0&&e.scratchBytes<=64u*1024u*1024u,
        "F16 production requested/actual native-row team differs");
    std::optional<HydroKepConvectionResult> nested;std::exception_ptr error;
#pragma omp parallel num_threads(2)
    {
#pragma omp single
        {try{nested.emplace(apply(values,sources,clock));}catch(...){error=std::current_exception();}}
    }
    if(error)std::rethrow_exception(error);check(nested.has_value(),"F16 production nested result missing");equalProductionApply(serial,*nested,boundaryBudget);
    check(!nested->openApplyReport().parallelUsed&&!nested->openApplyReport().parallelEligible,"F16 production nested request lost literal fallback");
    const auto sameFirstError=[&](const HydroNativeValues& state,HydroMaterialSources actualSources,LAMP::TimeLayer actualClock){
        std::string first;Final rejected;
        for(int workers:{1,2,8}){omp_set_num_threads(workers);reset();std::string message;
            try{(void)apply(state,actualSources,actualClock);}catch(const std::invalid_argument& ex){message=ex.what();}
            Final finished;finished.capture();check(!message.empty(),"F16 production stale/Closed request unexpectedly accepted");
            if(workers==1){first=message;rejected=finished;}else check(message==first&&rejected.same(finished),"F16 production common-input first error or FP differs");}
    };
    auto stale=clock;stale.begin=.001;sameFirstError(values,sources,stale);
    auto wrong=sources;std::swap(wrong.vofIds[0],wrong.vofIds[1]);sameFirstError(values,wrong,clock);
    auto invalidU=values.normalVelocity();check(!op.fixedVelocityDofs().empty(),"F16 production Closed guard fixture missing");invalidU[op.fixedVelocityDofs().front()]=1;
    const HydroNativeValues invalid(values.layout(),std::move(invalidU),values.pressure());sameFirstError(invalid,sources,clock);
    reset();const auto recovered=apply(values,sources,clock);equalProductionApply(serial,recovered,boundaryBudget);
    return actual;
#else
    const auto actual=apply(values,sources,clock);checkProductionReport(actual,boundaryBudget);return actual;
#endif
}
// Same existing large RZ fixture: one-worker build selects the original
// per-port checker, eight-worker qualified build selects the bounded batch.
HydroKepConvection pairedOutflowBuild(const QuadtreeMesh& mesh,const HydroDualMaterial& material,
    const HydroAdvectorBoundary& walls,const HydroConvectionBoundary& incoming,
    std::optional<HydroKepConvection>& literal)
{
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
    struct Runtime {int threads=omp_get_max_threads(),dynamic=omp_get_dynamic(),error=errno;std::fenv_t env{};unsigned csr=_mm_getcsr();
        Runtime(){check(std::fegetenv(&env)==0,"F06 original builder FP capture");}
        ~Runtime(){omp_set_dynamic(dynamic);omp_set_num_threads(threads);std::fesetenv(&env);_mm_setcsr(csr);errno=error;}} restore;
    omp_set_dynamic(0);std::fenv_t entry{};
    check(incoming.freePressureOutletVelocityDofs().size()>1,"F06 actual scene needs multiple free ports");
    check(std::feraiseexcept(FE_INEXACT)==0&&std::fegetenv(&entry)==0,"F06 prime builder FP receipt");const unsigned csr=_mm_getcsr();
    const auto reset=[&]{check(std::fesetenv(&entry)==0,"F06 common builder FP restore");_mm_setcsr(csr);errno=173;};
    struct Final {std::fenv_t env{};unsigned csr=0;int error=0;
        void capture(){error=errno;csr=_mm_getcsr();check(std::fegetenv(&env)==0,"F06 completed builder FP capture");}
        bool same(const Final& other)const{return error==other.error&&csr==other.csr&&std::memcmp(&env,&other.env,sizeof env)==0;}} sf,bf;
    omp_set_num_threads(1);reset();literal.emplace(HydroKepConvection::buildVofOpenFlow(mesh,material,walls,incoming));sf.capture();
    omp_set_num_threads(8);reset();auto actual=HydroKepConvection::buildVofOpenFlow(mesh,material,walls,incoming);bf.capture();
    check(sf.same(bf),"F06 original/batched builder full FP/CSR/errno differs");
    const auto a=literal->momentReport().maxNativeMeanDefect,b=actual.momentReport().maxNativeMeanDefect;
    check(a==b&&(a!=0||std::signbit(a)==std::signbit(b)),"F06 original/batched native mean certificate differs");
    check(literal->fixedVelocityDofs()==actual.fixedVelocityDofs(),"F06 batch changed fixed/free ownership");return actual;
#else
    (void)literal;return HydroKepConvection::buildVofOpenFlow(mesh,material,walls,incoming);
#endif
}
void run()
{
    const auto fail=[](auto&& action,const char* why){rejects(action,why);++guards;};
    for(bool rz:{false,true})for(unsigned pattern:{0u,1u,2u}){
        // One existing RZ ambient scene crosses all fixed256 chunk boundaries;
        // the other five original28-leaf scenes and every physical gate stay.
        auto mesh=rz&&pattern==1?[]{DictionaryMesh d;d.Nx=d.Ny=16;d.maxCellNum=4096;d.maxInitialIterations=0;d.maxLevel_local={2};
            auto result=std::make_unique<QuadtreeMesh>(d);DictionaryMesh::RefineBox box;box.level=1;box.min={.5,.5};box.max={1,1};
            result->applyRefineBoxes({box});return result;}():F15BoundaryGrossContract::makeMesh(rz);const unsigned phase=pattern%2;
        const auto source=Publication::prepareLocalPlic(*mesh,domain(*mesh,rz),{},std::vector<double>(mesh->cellIndexNum(),phase==0?1.:0.));
        const auto saved=source.plicCells();const auto sourceId=source.identity();const auto alpha=source.alpha();
        const auto a=F15BoundaryGrossContract::normals(*mesh,source,-.125);
        const auto b=F15BoundaryGrossContract::normals(*mesh,source,pattern?.375:.125);
        const auto rules=F15BoundaryGrossContract::rulesFor(*mesh,source,phase);
        const auto boundary=Slab::Boundary::prepare(*mesh,source,rules);constexpr double dt=1./32;
        const auto slab=Slab::prepareBoundary(*mesh,source,boundary,{1,0,dt,Position::End},a,b);
        constexpr double rhoL=100,rhoG=1.226;const double rho=phase==0?rhoL:rhoG;
        const auto volume=Transport::prepareLocalPlicBoundarySweep(*mesh,source,slab,rhoL,rhoG);
        const auto geometry=HydroMomentumGeometry::build(*mesh,vofMaterialMetrics(*mesh,source));
        const auto g0=HydroPhaseGeometry::build(*mesh,geometry.metrics(),source);
        const auto gm=HydroPhaseGeometry::build(*mesh,geometry.metrics(),volume.midpoint());
        const auto g1=HydroPhaseGeometry::build(*mesh,geometry.metrics(),volume.endpoint());
        const auto material=HydroMaterialTransport::prepareGeometricVof(*mesh,geometry,g0,gm,g1,volume);
        const auto m0=HydroMomentumMass::buildCellAverageProjection(*mesh,geometry,g0,rhoL,rhoG);
        const auto mh=HydroMomentumMass::buildCellAverageProjection(*mesh,geometry,gm,rhoL,rhoG);
        const auto m1=HydroMomentumMass::buildCellAverageProjection(*mesh,geometry,g1,rhoL,rhoG);
        const auto dual=HydroDualMaterial::prepare(*mesh,material,m0,mh,m1);const auto& layout=geometry.layout();
        ConservativeVofMomentumContract::qualify(*mesh,dual);
        const auto clock=dual.timeLayer();const auto sources=dual.sources();const auto fluxBefore=dual.faceRates(),integralBefore=dual.faceIntegrals();const auto cellsBefore=dual.cells();
        std::vector<int> walls,free;std::vector<double> u(layout.velocityCount());std::vector<Vector2D> ambient(source.domain().nativeFaceCount());
        const Vector2D declared=pattern==1?Vector2D(.0625,-.03125):Vector2D(0,0);
        const auto model=pattern==2?HydroVofPressureOutletMomentumModel::InteriorContinuationZeroShear:HydroVofPressureOutletMomentumModel::Ambient;
        for(int i=0;i<layout.velocityCount();++i){const int row=layout.nativeFace(i);const auto f=source.domain().face(*mesh,row);
            u[i]=layout.normalAxis(i)==HydroDofLayout::NormalAxis::First?.375+.125*f.center.x:.25+.125*f.center.x-.0625*f.center.y;
            if(f.negativeCell<0||f.positiveCell<0){if(rules[row]==Slab::Boundary::Mode::Closed){walls.push_back(i);u[i]=0;}
                else {free.push_back(i);ambient[row]=declared;}}}
        check(!free.empty()&&!walls.empty(),"F16 requires genuine free reservoirs and Closed walls");
        const auto wallLaw=HydroAdvectorBoundary::impermeable(*mesh,geometry,walls);
        std::vector<Vector2D> first(geometry.segments().size()),last(first);
        const auto incoming=HydroConvectionBoundary::bindVofOpenMomentum(*mesh,dual,first,last,free,ambient,sources,clock,model);
        check(incoming.hasVofOpenMomentum()&&incoming.freePressureOutletVelocityDofs()==free,"F16 lost explicit open momentum ports");
        std::optional<HydroKepConvection> literalOutflow;
        const auto op=rz&&pattern==1?pairedOutflowBuild(*mesh,dual,wallLaw,incoming,literalOutflow):
            HydroKepConvection::buildVofOpenFlow(*mesh,dual,wallLaw,incoming);
        check(op.isVofOpenFlow()&&op.fixedVelocityDofs()==walls,"F16 fixed the free native pressure normal");
        unsigned unusedEndpointCalls=0;
        const auto productionIncoming=HydroConvectionBoundary::bindVofOpenProductionMomentum(*mesh,dual,free,ambient,
            [&](int,LAMP::TimeLayer){++unusedEndpointCalls;return Vector2D{};},sources,clock,model);
        check(unusedEndpointCalls==0,"F16 all-free/Closed source called an unused endpoint callback");equalOpenBinding(incoming,productionIncoming);
        // Genuine accepted Closed zero measures on this same actual fixture,
        // excluding the axis which never needed a positive-area subface query.
        for(std::size_t id=0;id<geometry.segments().size();++id){const auto& g=geometry.segments()[id];
            if(g.negativeVelocity>=0&&g.positiveVelocity>=0)continue;
            const int row=layout.nativeFace(g.primalVelocity);const auto face=source.domain().face(*mesh,row);
            if(rules[row]!=Slab::Boundary::Mode::Closed||face.area==0)continue;
            check(face.area>0&&(face.negativeCell>=0)!=(face.positiveCell>=0)
                &&row>=0&&static_cast<std::size_t>(row)<volume.faces().size()&&volume.faces()[row].nativeRow==row,
                "F16 Closed-zero witness lacks its actual bounded active native face");
            const auto& q=volume.faces()[row].full;
            const std::array<long double,13> acceptedZero{{q.total,q.positiveTotal,q.negativeTotal,q.rawLiquid,q.liquid,q.gas,q.mass,
                q.positiveLiquid,q.negativeLiquid,q.positiveGas,q.negativeGas,q.exhaustionLiquidCorrection,q.exhaustionGasCorrection}};
            for(auto value:acceptedZero)check(value==0,"F16 Closed-zero witness has nonzero accepted full gross or correction");
            check(q.limiter==1&&slab.beginNormal()[row]==0&&slab.endNormal()[row]==0&&dual.faceRates()[id]==0&&dual.faceIntegrals()[id]==0,
                "F16 Closed-zero witness uses a limiter, normal or segment-flux approximation");
            const auto& actual=productionIncoming.vofOpenMomentumFluxes()[id];
            for(auto value:std::array<long double,7>{{actual.outgoingPhaseIntegral[0],actual.outgoingPhaseIntegral[1],
                actual.incomingPhaseIntegral[0],actual.incomingPhaseIntegral[1],actual.outgoingMassRate,actual.incomingMassRate,actual.incomingMomentum}})
                check(value==0,"F16 actual Closed-zero production payload is nonzero");
            check(actual.canonicalOutwardMassRate==0,"F16 actual Closed-zero production canonical rate is nonzero");
            ++productionClosedZeroWitnesses;if(std::signbit(actual.canonicalOutwardMassRate))++productionClosedNegativeZeros;else ++productionClosedPositiveZeros;
        }
        if(rz&&pattern==1){
            const auto sameSourceFailure=[&](HydroMaterialSources bad,LAMP::TimeLayer when,const char* expected){
                std::string reference,actual;
                try{dual.validateFor(*mesh,bad,when);}catch(const std::invalid_argument& error){reference=error.what();}
                try{(void)HydroConvectionBoundary::bindVofOpenProductionMomentum(*mesh,dual,free,ambient,
                    [&](int,LAMP::TimeLayer){++unusedEndpointCalls;return Vector2D{};},bad,when);}
                catch(const std::invalid_argument& error){actual=error.what();}
                check(!reference.empty()&&actual==reference&&actual.find(expected)!=std::string::npos,
                    "F16 production binder changed original source/full-clock guard or accepted stale input");++guards;
            };
            auto stale=clock;stale.begin=.001;sameSourceFailure(sources,stale,"full clock mismatch");
            auto bad=sources;bad.vofIds[0]=volume.endpoint().identity();sameSourceFailure(bad,clock,"foreign VOF owner");
            bad=sources;bad.vofStamp.reset();sameSourceFailure(bad,clock,"different accepted VOF transport instance");
            int axis=-1;for(int i:walls)if(geometry.metrics().faceMeasures()[i]==0){axis=i;break;}
            check(axis>=0,"F16 existing RZ fixture lacks an actual zero-measure axis port");auto badPorts=free;badPorts.push_back(axis);
            std::string reference,actual;
            try{(void)HydroConvectionBoundary::bindVofOpenMomentum(*mesh,dual,first,last,badPorts,ambient,sources,clock,model);}
            catch(const std::invalid_argument& error){reference=error.what();}
            try{(void)HydroConvectionBoundary::bindVofOpenProductionMomentum(*mesh,dual,badPorts,ambient,
                [&](int,LAMP::TimeLayer){++unusedEndpointCalls;return Vector2D{};},sources,clock,model);}
            catch(const std::invalid_argument& error){actual=error.what();}
            check(!reference.empty()&&actual==reference&&actual.find("ambient pressure port is internal or axis")!=std::string::npos,
                "F16 production binder queried a degenerate axis instead of retaining its real port guard");++guards;
            check(unusedEndpointCalls==0,"F16 rejected source/clock/axis binder invoked an unused endpoint provider");
        }
        const auto productionOp=HydroKepConvection::buildVofOpenProductionFlow(*mesh,dual,wallLaw,productionIncoming);equalPortCertificate(op,productionOp);
        check(productionOp.points().empty()&&productionOp.weights().empty()&&productionOp.fixedVelocityDofs()==op.fixedVelocityDofs(),
            "F16 production factory retained diagnostic tables or changed fixed rows");
        const HydroNativeValues values(layout,u,std::vector<double>(layout.pressureCount()));
        const auto applied=pairedApply(*mesh,op,values,sources,clock);
        const auto omitted=pairedApply(*mesh,op,values,sources,clock,ApplyRequest::WithoutCenter);
        equalOpenApply(applied,omitted,false);++omittedPhysicalPairs;
        const auto momentum=pairedApply(*mesh,op,values,sources,clock,ApplyRequest::MomentumLoads);
        equalMomentumApply(applied,momentum);++momentumPhysicalPairs;
        for(bool budget:{false,true}){const auto production=pairedProductionApply(*mesh,productionOp,values,sources,clock,budget);
            equalProductionApply(applied,production,budget);++productionPhysicalPairs;}
        const auto H=op.frozenMomentumOperator(*mesh,sources,clock);
        const auto productionH=productionOp.frozenMomentumOperator(*mesh,sources,clock);
        check(productionH.rows()==H.rows()&&productionH.cols()==H.cols()&&productionH.nonZeros()==H.nonZeros()
            &&productionOp.momentReport().physicalPairs==op.momentReport().physicalPairs,"F16 production factory shared matrix dimensions differ");
        for(int outer=0;outer<H.outerSize();++outer){HydroKineticHodge::Matrix::InnerIterator a(H,outer),b(productionH,outer);
            for(;a&&b;++a,++b)check(a.row()==b.row()&&a.col()==b.col()&&a.value()==b.value(),"F16 production factory changed actual shared coefficient");
            check(!a&&!b,"F16 production factory changed sparse column length");}
        equalOpenApply(applied,productionOp.applyVofOpen(*mesh,values,sources,clock));
        if(literalOutflow){
            equalOpenApply(literalOutflow->applyVofOpen(*mesh,values,sources,clock),applied);
            const auto originalH=literalOutflow->frozenMomentumOperator(*mesh,sources,clock);
            check(originalH.rows()==H.rows()&&originalH.cols()==H.cols()&&originalH.nonZeros()==H.nonZeros(),"F06 frozen matrix structure differs");
            for(int outer=0;outer<H.outerSize();++outer){HydroKineticHodge::Matrix::InnerIterator a(originalH,outer),b(H,outer);
                for(;a&&b;++a,++b)check(a.row()==b.row()&&a.col()==b.col()&&a.value()==b.value()
                    &&(a.value()!=0||std::signbit(a.value())==std::signbit(b.value())),"F06 frozen original matrix entry differs");
                check(!a&&!b,"F06 frozen matrix incidence differs");}
            const auto originalB=literalOutflow->frozenBoundaryMomentumLoad(*mesh,sources,clock),actualB=op.frozenBoundaryMomentumLoad(*mesh,sources,clock);
            check(originalB.size()==actualB.size(),"F06 frozen boundary dimensions differ");
            for(std::size_t i=0;i<originalB.size();++i)check(originalB[i]==actualB[i]
                &&(originalB[i]!=0||std::signbit(originalB[i])==std::signbit(actualB[i])),"F06 frozen original boundary value differs");
        }
        // The SAME immutable operator must consume each actual U afresh.
        // Compare every native/scalar output with the independent literal
        // route, retaining its FP/errno, nested fallback and stale-clock checks.
        auto changedU=u;for(std::size_t i=0;i<changedU.size();++i)if(u[i]!=0)
            changedU[i]=-.5*u[i]+(layout.normalAxis(static_cast<int>(i))==HydroDofLayout::NormalAxis::First?.03125:-.0625);
        const HydroNativeValues changedValues(layout,changedU,std::vector<double>(layout.pressureCount()));
        const auto changedApplied=pairedApply(*mesh,op,changedValues,sources,clock);
        const auto changedOmitted=pairedApply(*mesh,op,changedValues,sources,clock,ApplyRequest::WithoutCenter);
        equalOpenApply(changedApplied,changedOmitted,false);++omittedPhysicalPairs;
        const auto changedMomentum=pairedApply(*mesh,op,changedValues,sources,clock,ApplyRequest::MomentumLoads);
        equalMomentumApply(changedApplied,changedMomentum);++momentumPhysicalPairs;
        for(bool budget:{false,true}){const auto production=pairedProductionApply(*mesh,productionOp,changedValues,sources,clock,budget);
            equalProductionApply(changedApplied,production,budget);++productionPhysicalPairs;}
        bool changedRate=false;for(std::size_t i=0;i<u.size();++i)changedRate=changedRate||changedApplied.momentumRates()[i]!=applied.momentumRates()[i];
        check(changedRate,"F16 same operator reused old-U ledger values");
        const auto repeated=op.applyVofOpen(*mesh,values,sources,clock);equalOpenApply(applied,repeated);
        const auto repeatedOmitted=op.applyVofOpenWithoutCenterDiagnostic(*mesh,values,sources,clock);equalOpenApply(omitted,repeatedOmitted);
        const auto repeatedMomentum=op.applyVofOpenMomentumLoads(*mesh,values,sources,clock);equalMomentumApply(momentum,repeatedMomentum);
        equalProductionApply(applied,productionOp.applyVofOpenProductionLoads(*mesh,values,sources,clock),false);
        equalProductionApply(applied,productionOp.applyVofOpenBoundaryBudget(*mesh,values,sources,clock),true);
        if(rz&&pattern==1){const auto& actual=applied.openApplyReport();check(actual.nativeRows>256&&actual.segments>256&&actual.pairs>256,
            "F16 RZ ambient must exercise actual multi-chunk native rows/segments/pairs");
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
            // A finite public U reaches the original work-range guard AFTER
            // contribution/gather phases. No injected worker or solver failure.
            struct Runtime {int threads=omp_get_max_threads(),dynamic=omp_get_dynamic(),error=errno;std::fenv_t env{};unsigned csr=_mm_getcsr();
                Runtime(){check(std::fegetenv(&env)==0,"F02 original runtime FP capture");}
                ~Runtime(){omp_set_dynamic(dynamic);omp_set_num_threads(threads);std::fesetenv(&env);_mm_setcsr(csr);errno=error;}} restore;
            omp_set_dynamic(0);std::fenv_t entry{};
            check(std::feraiseexcept(FE_INEXACT)==0&&std::fegetenv(&entry)==0,"F02 prime original inexact");const unsigned entryCsr=_mm_getcsr();
            const auto reset=[&]{check(std::fesetenv(&entry)==0,"F02 common FP restore");_mm_setcsr(entryCsr);errno=173;};
            auto hugeU=u;for(auto& value:hugeU)value*=1e200;
            const HydroNativeValues hugeValues(layout,std::move(hugeU),std::vector<double>(layout.pressureCount()));
            for(bool includeCenterDiagnostic:{true,false}){
            std::string literalError;std::fenv_t literalFp{};unsigned literalCsr=0;int literalErrno=0;
            for(int workers:{1,2,8}){omp_set_num_threads(workers);reset();std::string actualError;
                try{if(includeCenterDiagnostic)(void)op.applyVofOpen(*mesh,hugeValues,sources,clock);
                    else (void)op.applyVofOpenWithoutCenterDiagnostic(*mesh,hugeValues,sources,clock);}
                catch(const std::invalid_argument& error){actualError=error.what();}
                check(actualError.find("[HydroFiniteContraction] final contraction overflow")!=std::string::npos,
                    "F02 finite large-U witness did not reach the original exact work-range rejection");
                const int actualErrno=errno;const unsigned actualCsr=_mm_getcsr();std::fenv_t actualFp{};
                check(std::fegetenv(&actualFp)==0,"F02 rejected final FP capture");
                if(workers==1){literalError=actualError;literalFp=actualFp;literalCsr=actualCsr;literalErrno=actualErrno;}
                else check(actualError==literalError&&actualErrno==literalErrno&&actualCsr==literalCsr
                    &&std::memcmp(&actualFp,&literalFp,sizeof actualFp)==0,"F02 partial-apply replay changed original error or full FP/CSR/errno");
                reset();const auto recovered=includeCenterDiagnostic?op.applyVofOpen(*mesh,values,sources,clock)
                    :op.applyVofOpenWithoutCenterDiagnostic(*mesh,values,sources,clock);
                equalOpenApply(includeCenterDiagnostic?applied:omitted,recovered);
                if(!includeCenterDiagnostic)++omittedOverflowChecks;
                if(workers>1){const auto& receipt=recovered.openApplyReport();check(receipt.parallelUsed&&!receipt.parallelFallback
                    &&receipt.contributionWorkers>=2&&receipt.gatherWorkers>=2,"F02 recovery apply lost actual team consumption");}}
            }
            // Same existing finite huge-U fixture: powers overflow in BOTH
            // old capabilities above, but the actually requested loads are finite.
            reset();const auto hugeMomentum=pairedApply(*mesh,op,hugeValues,sources,clock,ApplyRequest::MomentumLoads);
            checkMomentumOnlyReport(hugeMomentum);
            for(const auto* rows:{&hugeMomentum.momentumRates(),&hugeMomentum.momentumLoads(),&hugeMomentum.dissipationRates()})
                for(auto value:*rows)check(std::isfinite(value),"F16 minimal huge-U finite load positive control failed");
            ++momentumHugePositiveScenes;
            const auto recoveredMomentum=op.applyVofOpenMomentumLoads(*mesh,values,sources,clock);equalMomentumApply(momentum,recoveredMomentum);
            for(bool budget:{false,true}){
                if(!budget){const auto production=pairedProductionApply(*mesh,productionOp,hugeValues,sources,clock,false);equalProductionApply(hugeMomentum,production,false);}
                else {std::string first;std::fenv_t firstFp{};unsigned firstCsr=0;int firstErrno=0;
                    for(int workers:{1,2,8}){omp_set_num_threads(workers);reset();std::string message;
                        try{(void)productionOp.applyVofOpenBoundaryBudget(*mesh,hugeValues,sources,clock);}catch(const std::invalid_argument& ex){message=ex.what();}
                        const int error=errno;const unsigned csr=_mm_getcsr();std::fenv_t fp{};check(std::fegetenv(&fp)==0,"F16 production rejected-budget FP capture");
                        check(message.find("[HydroFiniteContraction] final contraction overflow")!=std::string::npos,"F16 production requested huge boundary power did not reject");
                        if(workers==1){first=message;firstFp=fp;firstCsr=csr;firstErrno=error;}
                        else check(message==first&&error==firstErrno&&csr==firstCsr&&std::memcmp(&fp,&firstFp,sizeof fp)==0,"F16 production budget rejection/replay differs");
                        reset();equalProductionApply(applied,productionOp.applyVofOpenBoundaryBudget(*mesh,values,sources,clock),true);
                    }
                }
            }

#endif
        }
        if(rz&&pattern==1){
            // Same actual fixture/operator: two nonzero interior contributions
            // cancel exactly before terminal RN. Coefficients are taken from
            // the retained full operator; no rounded native product is the oracle.
            int cancelledRow=-1,firstColumn=-1,secondColumn=-1;double firstCoefficient=0,secondCoefficient=0;
            for(int row=0;row<H.rows()&&cancelledRow<0;++row){
                bool inside=true;for(const auto& side:geometry.volumes()[row].sides)for(int id:side){const auto& g=geometry.segments()[id];
                    inside=inside&&g.negativeVelocity>=0&&g.positiveVelocity>=0;}
                if(!inside)continue;int first=-1;double coefficient=0;
                for(int column=0;column<H.cols();++column){if(column==row||std::binary_search(walls.begin(),walls.end(),column))continue;
                    const double value=H.coeff(row,column);if(value==0)continue;
                    if(first<0){first=column;coefficient=value;}
                    else {cancelledRow=row;firstColumn=first;secondColumn=column;firstCoefficient=coefficient;secondCoefficient=value;break;}}
            }
            check(cancelledRow>=0&&firstCoefficient!=0&&secondCoefficient!=0,"F16 existing RZ fixture lacks two real nonzero interior terms");
            std::vector<double> cancelledU(u.size());cancelledU[firstColumn]=secondCoefficient;cancelledU[secondColumn]=-firstCoefficient;
            const HydroNativeValues cancelledValues(layout,std::move(cancelledU),std::vector<double>(layout.pressureCount()));
            const auto cancelledReference=op.applyVofOpen(*mesh,cancelledValues,sources,clock);
            check(cancelledReference.momentumRates()[cancelledRow]==0&&!std::signbit(cancelledReference.momentumRates()[cancelledRow])
                &&cancelledReference.momentumLoads()[cancelledRow]==0&&!std::signbit(cancelledReference.momentumLoads()[cancelledRow]),
                "F16 two nonzero interior products did not cancel to exact positive zero in the original Full route");
            for(bool budget:{false,true})equalProductionApply(cancelledReference,
                pairedProductionApply(*mesh,productionOp,cancelledValues,sources,clock,budget),budget);
            // Reuse this existing real 448-leaf RZ fixture. Choose a nonzero
            // off-diagonal binary64 coefficient and a genuinely interior target.
            // Then one denorm_min seed has exact nonzero N below half a quantum.
            if constexpr(std::numeric_limits<long double>::digits==std::numeric_limits<double>::digits
                &&std::numeric_limits<long double>::min_exponent==std::numeric_limits<double>::min_exponent){
                int target=-1,seed=-1;
                for(int outer=0;outer<H.outerSize()&&target<0;++outer)for(HydroKineticHodge::Matrix::InnerIterator entry(H,outer);entry;++entry){
                    if(entry.row()==entry.col()||entry.value()==0||std::abs(entry.value())>=.5
                        ||std::binary_search(walls.begin(),walls.end(),static_cast<int>(entry.col())))continue;
                    bool inside=true;for(const auto& side:geometry.volumes()[entry.row()].sides)for(int id:side){const auto& g=geometry.segments()[id];
                        inside=inside&&g.negativeVelocity>=0&&g.positiveVelocity>=0;}
                    if(inside){target=static_cast<int>(entry.row());seed=static_cast<int>(entry.col());break;}
                }
                check(target>=0&&seed>=0,"F16 existing RZ fixture lacks selected small interior coefficient");
                std::vector<double> tinyU(u.size());const HydroNativeValues zeroValues(layout,tinyU,std::vector<double>(layout.pressureCount()));
                const auto zeroReference=op.applyVofOpen(*mesh,zeroValues,sources,clock);
                check(zeroReference.momentumRates()[target]==0&&!std::signbit(zeroReference.momentumRates()[target])
                    &&zeroReference.momentumLoads()[target]==0&&!std::signbit(zeroReference.momentumLoads()[target]),"F16 selected exact-zero row lost positive zero");
                for(bool budget:{false,true})equalProductionApply(zeroReference,pairedProductionApply(*mesh,productionOp,zeroValues,sources,clock,budget),budget);
                ++productionZeroWitnesses;
                for(double sign:{-1.,1.}){tinyU[seed]=std::copysign(std::numeric_limits<double>::denorm_min(),sign);
                    const HydroNativeValues tinyValues(layout,tinyU,std::vector<double>(layout.pressureCount()));
                    const auto tinyReference=op.applyVofOpen(*mesh,tinyValues,sources,clock);
                    check(tinyReference.momentumRates()[target]==0&&tinyReference.momentumLoads()[target]==0
                        &&std::signbit(tinyReference.momentumRates()[target])!=std::signbit(tinyReference.momentumLoads()[target]),
                        "F16 selected exact nonzero subquantum row did not witness opposite signed RN zeros");
                    for(bool budget:{false,true})equalProductionApply(tinyReference,pairedProductionApply(*mesh,productionOp,tinyValues,sources,clock,budget),budget);
                    ++productionUnderflowWitnesses;
                }
            }
            // Exercise the production-prepared shortcut on a publicly mutable
            // mesh. Immutable source tags cannot authenticate changed live words.
            auto& movedVertex=mesh->verts(mesh->vertIndexToLabel[0]);const auto savedCoordinate=movedVertex.coord();
            movedVertex.setCoord(savedCoordinate+Vector2D{1./1024,0});
            try{for(bool budget:{false,true}){bool rejected=false;
                // This new binder has no retained incoming receipt yet: its
                // public Full restriction must reject stale mesh before queries.
                if(!budget){bool rejectedBinding=false;
                    try{(void)HydroConvectionBoundary::bindVofOpenProductionMomentum(*mesh,dual,free,ambient,
                        [&](int,LAMP::TimeLayer){++unusedEndpointCalls;return Vector2D{};},sources,clock,model);}
                    catch(const std::invalid_argument& error){rejectedBinding=std::string(error.what()).find("mesh geometry/topology changed since identity construction")!=std::string::npos;
                        check(rejectedBinding,"F16 production pre-query binder selected a different stale-mesh guard");}
                    check(rejectedBinding,"F16 production binder used immutable segments without checking the actual live mesh");++guards;
                }
                try{if(budget)(void)productionOp.applyVofOpenBoundaryBudget(*mesh,values,sources,clock);
                    else (void)productionOp.applyVofOpenProductionLoads(*mesh,values,sources,clock);}
                catch(const std::invalid_argument& error){rejected=std::string(error.what()).find("mesh geometry/topology changed since identity construction")!=std::string::npos;
                    check(rejected,"F16 production changed-mesh request selected a different guard");}
                check(rejected,"F16 production prepared source bypassed live mutable mesh validation");++guards;
            }
                bool buildRejected=false;try{(void)HydroKepConvection::buildVofOpenProductionFlow(*mesh,dual,wallLaw,productionIncoming);}
                catch(const std::invalid_argument& error){buildRejected=std::string(error.what()).find("mesh geometry/topology changed since identity construction")!=std::string::npos;
                    check(buildRejected,"F16 production certified bind/build selected a different mesh guard");}
                check(buildRejected,"F16 production factory reused a bind-time live-mesh success");++guards;
            }catch(...){movedVertex.setCoord(savedCoordinate);throw;}
            movedVertex.setCoord(savedCoordinate);
            equalProductionApply(applied,productionOp.applyVofOpenProductionLoads(*mesh,values,sources,clock),false);
            const auto rebound=HydroConvectionBoundary::bindVofOpenProductionMomentum(*mesh,dual,free,ambient,
                [&](int,LAMP::TimeLayer){++unusedEndpointCalls;return Vector2D{};},sources,clock,model);
            equalOpenBinding(incoming,rebound);check(unusedEndpointCalls==0,"F16 recovered binder called unused endpoint data");
            equalProductionApply(applied,productionOp.applyVofOpenBoundaryBudget(*mesh,values,sources,clock),true);
            // Same source and geometry owner; only genuine slab/transport data
            // change. A cached active-pair list would fail 0 -> + -> 0 -> +.
            const auto route=[&](const HydroMomentumGeometry& routeGeometry,double firstNormal,double lastNormal){
                const auto r0=F15BoundaryGrossContract::normals(*mesh,source,firstNormal),r1=F15BoundaryGrossContract::normals(*mesh,source,lastNormal);
                const auto routeSlab=Slab::prepareBoundary(*mesh,source,boundary,{1,0,dt,Position::End},r0,r1);
                const auto routeVolume=Transport::prepareLocalPlicBoundarySweep(*mesh,source,routeSlab,rhoL,rhoG);
                const auto p0=HydroPhaseGeometry::build(*mesh,routeGeometry.metrics(),source);
                const auto pm=HydroPhaseGeometry::build(*mesh,routeGeometry.metrics(),routeVolume.midpoint());
                const auto p1=HydroPhaseGeometry::build(*mesh,routeGeometry.metrics(),routeVolume.endpoint());
                const auto transport=HydroMaterialTransport::prepareGeometricVof(*mesh,routeGeometry,p0,pm,p1,routeVolume);
                const auto mass=[&](const HydroPhaseGeometry& phaseGeometry){return HydroMomentumMass::buildCellAverageProjection(*mesh,routeGeometry,phaseGeometry,rhoL,rhoG);};
                const auto current=HydroDualMaterial::prepare(*mesh,transport,mass(p0),mass(pm),mass(p1));
                const auto routeWalls=HydroAdvectorBoundary::impermeable(*mesh,routeGeometry,walls);
                const std::vector<Vector2D> endpoints(routeGeometry.segments().size());
                const auto routeBoundary=HydroConvectionBoundary::bindVofOpenMomentum(*mesh,current,endpoints,endpoints,free,ambient,current.sources(),current.timeLayer(),model);
                const auto referenceOp=HydroKepConvection::buildVofOpenFlow(*mesh,current,routeWalls,routeBoundary);
                unsigned callbacks=0;const auto productionBoundary=HydroConvectionBoundary::bindVofOpenProductionMomentum(*mesh,current,free,ambient,
                    [&](int,LAMP::TimeLayer){++callbacks;return Vector2D{};},current.sources(),current.timeLayer(),model);
                check(callbacks==0,"F16 route-only free-port source called unused endpoint provider");equalOpenBinding(routeBoundary,productionBoundary);
                const auto currentOp=HydroKepConvection::buildVofOpenProductionFlow(*mesh,current,routeWalls,productionBoundary);equalPortCertificate(referenceOp,currentOp);
                const HydroNativeValues currentValues(routeGeometry.layout(),u,std::vector<double>(routeGeometry.layout().pressureCount()));
                const auto reference=referenceOp.applyVofOpen(*mesh,currentValues,current.sources(),current.timeLayer());
                const auto predicted=pairedProductionApply(*mesh,currentOp,currentValues,current.sources(),current.timeLayer(),false);
                auto actual=pairedProductionApply(*mesh,currentOp,currentValues,current.sources(),current.timeLayer(),true);
                equalProductionApply(reference,predicted,false);equalProductionApply(reference,actual,true);++productionRouteTransitions;
                return std::pair<int,HydroKepConvectionResult>{currentOp.momentReport().physicalPairs,std::move(actual)};
            };
            const auto zero=route(geometry,0,0),nonzero=route(geometry,.125,.125),cancelled=route(geometry,-.125,.125),again=route(geometry,.125,.125);
            check(zero.first==0&&nonzero.first>0&&cancelled.first==0&&again.first==nonzero.first,"F16 same-owner dynamic physical pair support was cached stale");
            // Independently rebuilt transports have fresh phase IDs; each was
            // already checked against its own reference/source above.
            equalProductionApply(nonzero.second,again.second,true,false);
            bool foreignRejected=false;try{(void)HydroKepConvection::buildVofOpenProductionFlow(*mesh,nonzero.second.material(),wallLaw,productionIncoming);}
            catch(const std::invalid_argument& error){foreignRejected=std::string(error.what()).find("different prepared material owner")!=std::string::npos;
                check(foreignRejected,"F16 foreign incoming certificate did not retain original owner guard");}
            check(foreignRejected,"F16 production factory reused a foreign material's port certificate");++guards;
            const auto rebuilt=HydroMomentumGeometry::build(*mesh,geometry.metrics());check(!rebuilt.sameIdentity(geometry),"F16 distinct geometry owner fixture missing");
            const auto fresh=route(rebuilt,.125,.125);check(fresh.first==nonzero.first,"F16 geometry owner rebuild lost physical pair support");
            equalProductionApply(nonzero.second,fresh.second,true,false);++productionOwnerRebuilds;
        }
        const auto B=op.frozenBoundaryMomentumLoad(*mesh,sources,clock);const Eigen::Map<const Eigen::VectorXd> U(u.data(),u.size());const Eigen::VectorXd hu=H*U;
        // The norm must retain each actual matrix term BEFORE cancellation.
        std::vector<E> matrixRowScale(u.size());
        for(int outer=0;outer<H.outerSize();++outer)
            for(HydroKineticHodge::Matrix::InnerIterator entry(H,outer);entry;++entry){
                E term;term.addProduct({static_cast<long double>(entry.value()),u[entry.col()]});
                matrixRowScale[entry.row()].addAbsolute(term);
            }
        const auto& receipt=incoming.vofOpenMomentumFluxes();const auto& segments=geometry.segments();
        std::array<E,2> n,p,d,w,rep,nativeRepresentation,scale,diagonalExpected;std::vector<E> externalExpected(u.size());std::vector<long double> canonicalExterior(u.size());
        for(std::size_t id=0;id<segments.size();++id){const auto& g=segments[id];if(g.negativeVelocity>=0&&g.positiveVelocity>=0)continue;
            const int i=g.negativeVelocity>=0?g.negativeVelocity:g.positiveVelocity;const int c=layout.normalAxis(i)==HydroDofLayout::NormalAxis::First?0:1;
            const auto& q=receipt[id];const auto native=source.domain().face(*mesh,q.nativeRow);const int sign=g.negativeVelocity>=0?1:-1;
            canonicalExterior[i]+=sign*dual.faceRates()[id];
            const long double area=g.normalAxis==0?std::abs(static_cast<long double>(g.end.y)-g.begin.y)*(rz?CartesianMeasure::twoPi*g.begin.x:1.L):
                std::abs(static_cast<long double>(g.end.x)-g.begin.x)*(rz?.5L*CartesianMeasure::twoPi*(static_cast<long double>(g.begin.x)+g.end.x):1.L);
            const auto analytical=native.horizontal&&native.area>0?F15BoundaryGrossContract::analyticGross(-.125,pattern?.375:.125,dt,area,1):std::array<long double,2>{};
            const auto out=analytical[sign>0?0:1],in=analytical[sign>0?1:0];
            const long double geometricScale=out+in;
            for(unsigned ph=0;ph<2;++ph){const auto expectedOut=ph==phase?out:0.L,expectedIn=ph==phase?in:0.L;
                if(geometricScale==0)check(q.outgoingPhaseIntegral[ph]==0&&q.incomingPhaseIntegral[ph]==0,"F16 zero exterior gross acquired phase");
                else {near(q.outgoingPhaseIntegral[ph],expectedOut,geometricScale,"F16 independent actual exterior phase outgross",1e-12L);
                    near(q.incomingPhaseIntegral[ph],expectedIn,geometricScale,"F16 independent actual exterior phase ingross",1e-12L);}}
            check(q.canonicalOutwardMassRate==sign*dual.faceRates()[id],"F16 changed canonical exterior signed F");
            const long double value=u[i],ambientComponent=q.incomingModel==HydroVofPressureOutletMomentumModel::InteriorContinuationZeroShear?value:(c==0?declared.x:declared.y);
            E exactNativeWork;exactNativeWork.addProduct({-.5L,q.canonicalOutwardMassRate,value,value});
            for(unsigned ph=0;ph<2;++ph){const long double density=ph==0?rhoL:rhoG;
                exactNativeWork.addQuotientProduct({.5L,density,q.outgoingPhaseIntegral[ph],value,value},dt);
                exactNativeWork.addQuotientProduct({-.5L,density,q.incomingPhaseIntegral[ph],value,value},dt);}
            nativeRepresentation[c].add(exactNativeWork);
            E ni,pi,di,wi,ri;ni.addQuotientProduct({rho,out,value},dt);ni.addQuotientProduct({-rho,in,ambientComponent},dt);
            pi.addQuotientProduct({.5L,rho,out,value,value},dt);pi.addQuotientProduct({-.5L,rho,in,ambientComponent,ambientComponent},dt);
            di.addQuotientProduct({.5L,rho,in,value,value},dt);di.addQuotientProduct({-rho,in,value,ambientComponent},dt);di.addQuotientProduct({.5L,rho,in,ambientComponent,ambientComponent},dt);
            wi.addQuotientProduct({rho,out,value,value},dt);wi.addQuotientProduct({-rho,in,value,ambientComponent},dt);wi.addProduct({-.5L,q.canonicalOutwardMassRate,value,value});
            ri.addQuotientProduct({.5L,rho,out,value,value},dt);ri.addQuotientProduct({-.5L,rho,in,value,value},dt);ri.addProduct({-.5L,q.canonicalOutwardMassRate,value,value});
            diagonalExpected[c].addQuotientProduct({rho,out,value,value},dt);
            if(q.incomingModel==HydroVofPressureOutletMomentumModel::InteriorContinuationZeroShear)diagonalExpected[c].addQuotientProduct({-rho,in,value,value},dt);
            check(di.compare(E{})>=0,"F16 independent dissipation must be nonnegative");n[c].add(ni);p[c].add(pi);d[c].add(di);w[c].add(wi);rep[c].add(ri);externalExpected[i].add(ni);
            scale[c].addAbsolute(pi);scale[c].addAbsolute(di);scale[c].addAbsolute(wi);scale[c].addAbsolute(ri);
            if(pattern==0&&out>0&&in>0){check(q.canonicalOutwardMassRate==0,"F16 true symmetric reversal must have native net0");if(ni.compare(E{})!=0&&di.compare(E{})>0)++netZeroActive;}
            if(pattern==1&&in>0&&ambientComponent!=0)++nonzeroAmbient;
        }
        std::array<E,2> sumN,kinetic,materialPower;std::array<long double,2> rowScale{};
        for(int i=0;i<layout.velocityCount();++i){const int c=layout.normalAxis(i)==HydroDofLayout::NormalAxis::First?0:1;
            const auto actual=applied.momentumLoads()[i],expected=-static_cast<long double>(hu[i])+B[i];
            check(dual.usesConservativeVofMomentum()&&actual==-applied.momentumRates()[i],"F16 gross-open applied load retained uR");
            E matrixComparisonScale=matrixRowScale[i];matrixComparisonScale.add(std::abs(actual));matrixComparisonScale.add(std::abs(B[i]));
            const auto s=matrixComparisonScale.rounded();
            if(s==0)check(actual==expected,"F16 zero H/B row");else {near(actual,expected,s,"F16 same actual H/B gives applied load",1e-12L);peakMatrix=std::max(peakMatrix,std::abs(actual-expected)/s);}
            check(dual.cells()[i].integratedDefect==cellsBefore[i].integratedDefect&&dual.cells()[i].massChange==cellsBefore[i].massChange,"F16 changed actual R or deltaM");
            diagonalExpected[c].addProduct({-.5L,canonicalExterior[i],u[i],u[i]});
            E row;row.addProduct({u[i],applied.momentumLoads()[i]});row.addQuotientProduct({-.5L,u[i],u[i],dual.cells()[i].massChange},dt);kinetic[c].add(row);
            materialPower[c].addQuotientProduct({-.5L,u[i],u[i],dual.cells()[i].integratedDefect},dt);sumN[c].add(applied.momentumRates()[i]);rowScale[c]+=std::abs(applied.momentumRates()[i]);++checkedRows;
        }
        for(int c=0;c<2;++c){const auto& r=applied.report().components[c];
            if(pattern==2){check(r.inletDissipationPower==0,"F16 interior continuation must have zero incoming upwind dissipation");
                check(std::all_of(B.begin(),B.end(),[](long double q){return q==0;}),"F16 implicit continuation must not freeze an affine incoming load");}
            const auto s=scale[c].rounded();check(s>0,"F16 must exercise both component powers");
            near(r.boundaryMomentumRate,n[c].rounded(),std::max(rowScale[c],std::abs(n[c].rounded())),"F16 independently integrated exterior momentum",1e-12L);
            near(sumN[c].rounded(),n[c].rounded(),std::max(rowScale[c],std::abs(n[c].rounded())),"F16 momentum conservation to actual boundary",1e-12L);
            near(r.boundaryPhysicalKineticPower,p[c].rounded(),s,"F16 independent outward physical power",1e-12L);
            near(r.inletDissipationPower,d[c].rounded(),s,"F16 independent nonnegative incoming dissipation",1e-12L);
            near(r.boundaryKineticPower,w[c].rounded(),s,"F16 independent actual boundary power",1e-12L);
            near(r.boundaryFluxRoundingPower,rep[c].rounded(),s,"F16 independent canonical-gross signed work",1e-12L);
            const auto exactRepresentation=nativeRepresentation[c].rounded(),own=std::abs(exactRepresentation)+std::abs(r.boundaryFluxRoundingPower);
            if(own==0)check(r.boundaryFluxRoundingPower==0,"F16 zero native representation work");
            else near(r.boundaryFluxRoundingPower,exactRepresentation,own,"F16 native quantum representation work own-scale",1e-12L);
            maximumRepresentation=std::max(maximumRepresentation,std::abs(r.boundaryFluxRoundingPower));
            near(r.boundaryDiagonalPower,diagonalExpected[c].rounded(),s,"F16 independent actual exterior diagonal power",1e-12L);
            E decomposition;decomposition.add(w[c]);decomposition.add(p[c],-1);decomposition.add(d[c],-1);decomposition.add(rep[c],-1);
            near(decomposition.rounded(),0,s,"F16 independently derived P+D+Wrep identity",1e-12L);
            E balance;balance.add(kinetic[c]);balance.add(materialPower[c],-1);balance.add(w[c]);near(balance.rounded(),0,s,"F16 actual kinetic/material/boundary identity",1e-12L);
            check(r.momentumIdentity<=1e-12L&&r.kineticIdentity<=1e-12L&&r.inletDissipationPower>=0,"F16 original component gates");
            peakMomentum=std::max(peakMomentum,r.momentumIdentity);peakPower=std::max(peakPower,r.kineticIdentity);
        }
        if(rz&&pattern==1){const auto plain=HydroConvectionBoundary::bind(*mesh,dual,std::vector<double>(geometry.segments().size()),sources,clock);
            fail([&]{(void)HydroKepConvection::buildVofOpenProductionFlow(*mesh,dual,wallLaw,plain);},"F16 production factory accepted untyped boundary");}
        fail([&]{(void)HydroKepConvection::buildVofPrescribedFlow(*mesh,dual,wallLaw,incoming,free);},"F16 new law cannot silently enter original one-way factory");
        fail([&]{(void)op.apply(*mesh,values,incoming,sources,clock);},"F16 new law cannot silently use a single applied upwind trace");
        auto stale=clock;stale.begin=.001;fail([&]{(void)op.applyVofOpen(*mesh,values,sources,stale);},"F16 stale full clock accepted");
        auto wrongAmbient=ambient;const int wallRow=layout.nativeFace(walls.front());wrongAmbient[wallRow]={1,0};
        fail([&]{(void)HydroConvectionBoundary::bindVofOpenMomentum(*mesh,dual,first,last,free,wrongAmbient,sources,clock);},"F16 ambient outside free port accepted");
        check(dual.faceRates()==fluxBefore&&dual.faceIntegrals()==integralBefore&&source.identity()==sourceId&&source.alpha()==alpha,"F16 changed actual source F/identity/alpha");
        for(std::size_t ci=0;ci<saved.size();++ci)check(source.plicCells()[ci].inventory==saved[ci].inventory,"F16 changed independent source stocks");
        std::cout<<"LOCAL_PLIC_OPEN_MOMENTUM_SCENE rz="<<rz<<" pattern="<<pattern<<" free="<<free.size()<<" cells="<<mesh->cellIndexNum()<<'\n';++scenes;
    }
    check(netZeroActive>0&&nonzeroAmbient>0,"F16 lacks active net0 and nonzero ambient witnesses");
    check(omittedPhysicalPairs==12,"F16 complete/omitted original and changed-U coverage incomplete");
    check(momentumPhysicalPairs==12,"F16 complete/minimal original and changed-U coverage incomplete");
    check(productionPhysicalPairs==24&&productionRouteTransitions==5&&productionOwnerRebuilds==1,"F16 production physical/route/owner coverage incomplete");
    check(productionClosedZeroWitnesses>0&&productionClosedPositiveZeros>0&&productionClosedNegativeZeros>0,"F16 genuine Closed-zero/signed canonical witnesses incomplete");
    if constexpr(std::numeric_limits<long double>::digits==std::numeric_limits<double>::digits
        &&std::numeric_limits<long double>::min_exponent==std::numeric_limits<double>::min_exponent)
        check(productionZeroWitnesses==1&&productionUnderflowWitnesses==2,"F16 production exact/subquantum-zero witnesses incomplete");
    std::cout<<"LOCAL_PLIC_OPEN_MOMENTUM_PASS scenes="<<scenes<<" netZeroActive="<<netZeroActive<<" nonzeroAmbient="<<nonzeroAmbient<<" rows="<<checkedRows<<" guards="<<guards<<" matrix="<<peakMatrix<<" momentum="<<peakMomentum<<" power="<<peakPower<<" representation="<<maximumRepresentation<<" omitted_physical_pairs="<<omittedPhysicalPairs
        <<" omitted_overflow_checks="<<omittedOverflowChecks<<" momentum_physical_pairs="<<momentumPhysicalPairs
        <<" momentum_huge_positive_scenes="<<momentumHugePositiveScenes<<" production_physical_pairs="<<productionPhysicalPairs
        <<" production_zero_witnesses="<<productionZeroWitnesses<<" production_underflow_witnesses="<<productionUnderflowWitnesses
        <<" production_closed_zero_witnesses="<<productionClosedZeroWitnesses<<" production_closed_positive_zeros="<<productionClosedPositiveZeros
        <<" production_closed_negative_zeros="<<productionClosedNegativeZeros<<" production_route_transitions="<<productionRouteTransitions<<" production_owner_rebuilds="<<productionOwnerRebuilds<<'\n';
}
}
// END F16 EXPLICIT GROSS OPEN MOMENTUM CONTRACT
// BEGIN F17 ACTUAL OPT-IN FREE-PRESSURE MOMENTUM STEP CONTRACT
namespace F17OpenMomentumStepContract
{
using Mode=Slab::Boundary::Mode;
using Model=HydroVofPressureOutletMomentumModel;
unsigned scenes=0,acceptedSteps=0,guardChecks=0,freeIncoming=0,zeroInitialPorts=0,actualContinuation=0,ambientDifferent=0,phaseChecks=0,uRChecks=0;
long double peakMomentum=0,peakEnergy=0,peakBoundary=0;
template<class F> void guard(F&& action,const char* literal)
{
    bool caught=false;try{action();}catch(const std::exception& error){caught=std::string(error.what()).find(literal)!=std::string::npos;
        check(caught,std::string("F17 wrong rejection: ")+error.what()+" instead of "+literal);}
    check(caught,std::string("F17 missing rejection: ")+literal);++guardChecks;
}
void run()
{
    for(bool rz:{false,true})for(Model model:{Model::Ambient,Model::InteriorContinuationZeroShear}){
        stage=std::string("F17 actual Step ")+(rz?"RZ":"XY")+(model==Model::Ambient?" ambient":" current-U continuation");
        std::unique_ptr<QuadtreeMesh> mesh;
        if(!rz)mesh=vofMaterialMesh(4);
        else{DictionaryMesh d;d.Nx=d.Ny=4;d.minCoord={1,0};d.maxCoord={2,1};d.maxCellNum=4096;d.maxInitialIterations=0;d.maxLevel_local={1};
            mesh=std::make_unique<QuadtreeMesh>(d);DictionaryMesh::RefineBox b;b.level=1;b.min={1.5,0};b.max={2,1};mesh->applyRefineBoxes({b});}
        const auto initial=Publication::prepareLocalPlic(*mesh,domain(*mesh,rz),{},std::vector<double>(mesh->cellIndexNum(),0));
        const auto metrics=vofMaterialMetrics(*mesh,initial);const auto& layout=metrics.layout();const auto native=initial.domain();
        check(!native.uniformLeaves()&&mesh->cellIndexNum()==40,"F17 must use actual 2:1 source leaves");
        const auto profile=[rz](double x){const double q=x-(rz?1:0);return q<.25?-.0625:q<.5?0.:.0625;};
        HydroVofPrescribedFlow flow;flow.rulesByNativeRow.assign(native.nativeFaceCount(),Mode::Interior);
        flow.normalVelocity.assign(layout.velocityCount(),0);flow.pressureOutletBackflowModel=model;
        std::vector<double> velocity(layout.velocityCount());std::vector<bool> freePort(layout.velocityCount());
        unsigned initialNegative=0,initialZero=0;int positiveFeed=-1;
        for(int u=0;u<layout.velocityCount();++u){const int row=layout.nativeFace(u);const auto f=native.face(*mesh,row);
            if(f.horizontal)velocity[u]=profile(f.center.x);
            if(f.outwardSign==0)continue;
            if(!f.horizontal||f.area==0){flow.rulesByNativeRow[row]=Mode::Closed;continue;}
            if(f.outwardSign>0){flow.rulesByNativeRow[row]=Mode::PureGasReservoir;flow.pressureOutlets.push_back({u,.125});freePort[u]=true;
                initialNegative+=velocity[u]<0;initialZero+=velocity[u]==0;}
            else if(velocity[u]==0)flow.rulesByNativeRow[row]=Mode::Closed;
            else{flow.rulesByNativeRow[row]=Mode::PureGasReservoir;flow.prescribedVelocityDofs.push_back(u);flow.normalVelocity[u]=velocity[u];
                if(velocity[u]>0)positiveFeed=u;}}
        check(initialNegative>0&&initialZero>0&&positiveFeed>=0,"F17 requires genuine prescribed feed, negative and zero free normals");zeroInitialPorts+=initialZero;
        const auto* meshPtr=mesh.get();flow.endpointVelocity=[native,meshPtr,profile](int row,Clock endpoint){endpoint.time();const auto f=native.face(*meshPtr,row);
            return Vector2D(.015625,profile(f.center.x));};
        struct Query {int row;Clock clock;};const auto queries=std::make_shared<std::vector<Query>>();
        if(model==Model::Ambient)flow.pressureOutletBackflowMomentum=[queries](int row,Clock clock){clock.time();queries->push_back({row,clock});return Vector2D(.03125,.09375);};
        const HydroNativeValues initialValues(layout,velocity,std::vector<double>(layout.pressureCount(),.125));
        const auto phase=HydroPhaseGeometry::build(*mesh,metrics,initial);const auto initialStocks=initial.plicCells();
        const auto initialDiv=metrics.divergence().apply(*mesh,initialValues);const auto& di=metrics.divergence().matrix();
        for(int row=0;row<layout.pressureCount();++row){long double scale=0;
            for(HydroIntegralDivergence::Matrix::InnerIterator it(di,row);it;++it)scale+=std::abs(static_cast<long double>(it.value())*velocity[it.col()]);
            check(scale>0?std::abs(initialDiv[row])/scale<=1e-12L:initialDiv[row]==0,"F17 actual source DI is not divergence free");}
        auto legacy=flow;legacy.pressureOutletBackflowModel.reset();legacy.pressureOutletBackflowMomentum={};
        guard([&]{HydroTwoPhaseState::createVofPrescribedFlow(*mesh,phase,initialValues,legacy);},"accepted VOF pressure outlet reverses or has zero flow");
        auto payload=flow;payload.normalVelocity[flow.pressureOutlets.front().velocityDof]=.0625;
        guard([&]{HydroTwoPhaseState::createVofPrescribedFlow(*mesh,phase,initialValues,payload);},"free VOF pressure outlet has prescribed normal payload");
        if(model==Model::Ambient){auto missing=flow;missing.pressureOutletBackflowMomentum={};
            guard([&]{HydroTwoPhaseState::createVofPrescribedFlow(*mesh,phase,initialValues,missing);},"VOF ambient backflow momentum requires an explicit query");}
        auto unused=flow;unused.pressureOutletBackflowModel.reset();unused.pressureOutletBackflowMomentum=[](int,Clock){return Vector2D{};};
        guard([&]{HydroTwoPhaseState::createVofPrescribedFlow(*mesh,phase,initialValues,unused);},"VOF backflow callback requires an explicit momentum model");
        if(model==Model::InteriorContinuationZeroShear){auto wrongCallback=flow;wrongCallback.pressureOutletBackflowMomentum=[](int,Clock){return Vector2D{};};
            guard([&]{HydroTwoPhaseState::createVofPrescribedFlow(*mesh,phase,initialValues,wrongCallback);},"VOF interior continuation cannot consume ambient callback data");}
        auto accepted=HydroTwoPhaseState::createVofPrescribedFlow(*mesh,phase,initialValues,flow);accepted.validateFor(*mesh);
        HydroVofAcConfig config;config.liquidDensity=2;config.gasDensity=1;config.surfaceTension=config.liquidViscosity=config.gasViscosity=0;
        config.momentumDensityRepresentation=HydroMomentumMass::DensityRepresentation::CellAverageProjection;
        constexpr double dt=1./16384;std::size_t totalStages=0;
        for(int index=1;index<=2;++index){const auto old=accepted;const auto& source=old.phaseGeometry().vofPublication();const auto stocks=source.plicCells();
            const auto oldU=old.values().normalVelocity(),oldP=old.values().pressure();const auto* storage=source.plicCells().data();const auto id=source.identity();
            const Clock end{index,(index-1)*dt,index*dt,Position::End};auto wrong=end;++wrong.step;
            guard([&]{HydroVofStep::preparePrescribedFlow(*mesh,accepted,wrong,config);},"requires contiguous complete End clock");
            const auto step=HydroVofStep::preparePrescribedFlow(*mesh,accepted,end,config);totalStages+=step.report().stages.size();
            for(const auto& ac:step.report().stages){const auto& c=ac.openConvection;
                check(c.available&&c.advector.available&&c.solvedMidpoint.available&&c.advector.vofSources==c.solvedMidpoint.vofSources
                    &&sameClock(c.advector.timeLayer,c.solvedMidpoint.timeLayer)&&c.advector.timeLayer.step==end.step
                    &&c.advector.timeLayer.begin==end.begin&&c.advector.timeLayer.end==end.end,
                    "F17 lost actual advector or necessary solved-midpoint apply receipt");
                check(c.advector.gatherWorkers==0&&c.solvedMidpoint.gatherWorkers==0&&c.advector.diagonalSeconds==0&&c.solvedMidpoint.pairSeconds==0,
                    "F17 actual Step still consumed phased reference diagnostic kernel");
                check(c.advectorApplySeconds>=c.advector.totalSeconds&&c.solvedMidpointApplySeconds>=c.solvedMidpoint.totalSeconds,
                    "F17 actual call time does not contain apply subdivision");}
            const auto& frozen=step.frozenMidpoint();const auto& dual=frozen.material();const auto& volume=step.transport();const auto& um=step.solvedMidpointValues();
            const auto& endpoint=frozen.solveResult().values();const auto& report=step.report().stages.back();const auto clock=dual.timeLayer();const auto sources=dual.sources();
            check(frozen.convection().isVofOpenFlow()&&frozen.convection().isPrescribedFlow()&&frozen.appliedConvection().boundary().hasVofOpenMomentum(),
                "F17 actual original NS must consume explicit F16 rather than legacy trace KEP");
            check(!frozen.appliedConvection().hasCenterMomentumRates()&&frozen.appliedConvection().centerMomentumRates().empty(),
                "F17 actual HydroVofStep predictor requested unused center diagnostic");
            F16OpenMomentumContract::checkProductionReport(frozen.appliedConvection(),false);
            check(frozen.convection().points().empty()&&frozen.convection().weights().empty(),"F17 actual Step used full diagnostic factory");
            check(dual.material().vofTransport().sameInstance(volume)&&volume.sourceIdentity()==id&&sameClock(volume.endTime(),end),"F17 unique actual F1/material source");
            check(frozen.viscosityMatrix().nonZeros()==0&&!step.hasCapillaryLoad()&&!frozen.hasMaxwellLoad(),"F17 mu0 exact zero/no unrelated forces");
            const auto reassessed=frozen.convection().applyVofOpen(*mesh,um,sources,clock);
            check(reassessed.hasCenterMomentumRates()&&!reassessed.centerMomentumRates().empty(),
                "F17 explicit original full API lost center diagnostic");
            check(reassessed.hasComponentReport()&&reassessed.hasPowerReport(),"F17 independent explicit full reassessment lost original power report");
            F16OpenMomentumContract::equalProductionApply(reassessed,frozen.convection().applyVofOpenBoundaryBudget(*mesh,um,sources,clock),true);
            const auto& typed=frozen.appliedConvection().boundary();const auto& gross=typed.vofOpenMomentumFluxes();const auto& segments=dual.geometry().segments();
            std::vector<GeometricVofNativeSubface> nativeQueries;std::vector<std::size_t> segmentIds;
            for(std::size_t g=0;g<segments.size();++g){const auto& s=segments[g];if(s.negativeVelocity>=0&&s.positiveVelocity>=0)continue;
                const auto f=native.face(*mesh,layout.nativeFace(s.primalVelocity));if(f.area==0)continue;
                nativeQueries.push_back({f.nativeRow,{s.begin.x,s.begin.y},{s.end.x,s.end.y}});segmentIds.push_back(g);}
            const auto actual=volume.restrictNativeFaces(*mesh,nativeQueries);
            if(rz&&model==Model::Ambient&&index==1){
                const auto& ports=typed.freePressureOutletVelocityDofs();std::vector<Vector2D> ambient(native.nativeFaceCount());
                for(int u:ports)ambient[layout.nativeFace(u)]={.03125,.09375};
                const auto referenceBinding=HydroConvectionBoundary::bindVofOpenMomentum(*mesh,dual,typed.beginVelocity(),typed.endVelocity(),ports,ambient,sources,clock,model);
                std::vector<Query> endpointCalls,expectedCalls;std::vector<bool> seen(native.nativeFaceCount());
                const auto& refGross=referenceBinding.vofOpenMomentumFluxes();
                for(std::size_t g=0;g<segments.size();++g){const auto& s=segments[g];
                    if(s.negativeVelocity>=0&&s.positiveVelocity>=0)continue;
                    const auto& q=refGross[g];if(q.nativeRow<0||freePort[s.primalVelocity]||seen[q.nativeRow]
                        ||!(q.incomingPhaseIntegral[0]>0||q.incomingPhaseIntegral[1]>0)||native.face(*mesh,q.nativeRow).area==0)continue;
                    seen[q.nativeRow]=true;expectedCalls.push_back({q.nativeRow,volume.sourceTime()});expectedCalls.push_back({q.nativeRow,volume.endTime()});
                }
                check(!expectedCalls.empty(),"F17 existing RZ feed lacks actual incoming nonfree callback witness");
                const auto productionBinding=HydroConvectionBoundary::bindVofOpenProductionMomentum(*mesh,dual,ports,ambient,
                    [&](int row,Clock endpoint){endpointCalls.push_back({row,endpoint});return flow.endpointVelocity(row,endpoint);},sources,clock,model);
                F16OpenMomentumContract::equalOpenBinding(referenceBinding,productionBinding);
                check(endpointCalls.size()==expectedCalls.size(),"F17 production endpoint callback count differs from actual incoming native rows");
                for(std::size_t q=0;q<endpointCalls.size();++q)check(endpointCalls[q].row==expectedCalls[q].row&&sameClock(endpointCalls[q].clock,expectedCalls[q].clock),
                    "F17 production callback did not retain segment-first source/end native-row order");
                const auto referenceOp=HydroKepConvection::buildVofOpenFlow(*mesh,dual,old.boundary(),referenceBinding);
                const auto productionOp=HydroKepConvection::buildVofOpenProductionFlow(*mesh,dual,old.boundary(),productionBinding);
                F16OpenMomentumContract::equalPortCertificate(referenceOp,productionOp);
                const auto full=referenceOp.applyVofOpen(*mesh,um,sources,clock);
                F16OpenMomentumContract::equalProductionApply(full,productionOp.applyVofOpenProductionLoads(*mesh,um,sources,clock),false);
                F16OpenMomentumContract::equalProductionApply(full,productionOp.applyVofOpenBoundaryBudget(*mesh,um,sources,clock),true);
                unsigned throwCalls=0;
                guard([&]{(void)HydroConvectionBoundary::bindVofOpenProductionMomentum(*mesh,dual,ports,ambient,
                    [&](int,Clock)->Vector2D{++throwCalls;throw std::runtime_error("F17 endpoint exception witness");},sources,clock,model);},"F17 endpoint exception witness");
                check(throwCalls==1,"F17 production binder replayed an endpoint callback exception");
                auto& moved=mesh->verts(mesh->vertIndexToLabel[0]);const auto coordinate=moved.coord();unsigned changedCalls=0;bool changed=false;
                try{guard([&]{(void)HydroConvectionBoundary::bindVofOpenProductionMomentum(*mesh,dual,ports,ambient,
                    [&](int row,Clock endpoint){endpoint.time();++changedCalls;if(!changed){moved.setCoord(coordinate+Vector2D{1./1024,0});changed=true;}
                        return Vector2D(.015625,flow.normalVelocity[layout.velocityDof(row)]);},sources,clock,model);},"mesh geometry/topology changed since identity construction");}
                catch(...){moved.setCoord(coordinate);throw;}
                moved.setCoord(coordinate);check(changed&&changedCalls>0,"F17 callback mutation control never entered actual incoming provider");
                const auto recovered=HydroConvectionBoundary::bindVofOpenProductionMomentum(*mesh,dual,ports,ambient,flow.endpointVelocity,sources,clock,model);
                F16OpenMomentumContract::equalOpenBinding(referenceBinding,recovered);
                F16OpenMomentumContract::equalPortCertificate(referenceOp,HydroKepConvection::buildVofOpenProductionFlow(*mesh,dual,old.boundary(),recovered));
            }
            std::array<long double,2> momentum{};long double physical=0,dissipation=0,representation=0,boundaryWork=0,diagonal=0,boundaryScale=0;
            std::vector<long double> exterior(layout.velocityCount());
            for(std::size_t k=0;k<segmentIds.size();++k){const auto g=segmentIds[k];const auto& s=segments[g];const auto& a=actual[k].full;const int sign=s.negativeVelocity>=0?1:-1;
                const int owner=s.negativeVelocity>=0?s.negativeVelocity:s.positiveVelocity,c=layout.normalAxis(owner)==HydroDofLayout::NormalAxis::First?0:1;
                const auto incoming=sign>0?std::array<long double,2>{a.negativeLiquid,a.negativeGas}:std::array<long double,2>{a.positiveLiquid,a.positiveGas};
                const auto outgoing=sign>0?std::array<long double,2>{a.positiveLiquid,a.positiveGas}:std::array<long double,2>{a.negativeLiquid,a.negativeGas};
                check(gross[g].incomingPhaseIntegral==incoming&&gross[g].outgoingPhaseIntegral==outgoing,"F17 typed momentum must retain same actual two-phase subface gross");
                const long double pin=(2*incoming[0]+incoming[1])/dt,pout=(2*outgoing[0]+outgoing[1])/dt,q=sign*dual.faceRates()[g],v=um.normalVelocity()[owner];
                long double b=0;
                if(freePort[s.primalVelocity]){check(typed.beginVelocity()[g].x==0&&typed.beginVelocity()[g].y==0
                    &&typed.endVelocity()[g].x==0&&typed.endVelocity()[g].y==0,"F17 free pressure momentum must not prescribe endpoint velocity");
                    if(pin>0){++freeIncoming;
                        if(model==Model::InteriorContinuationZeroShear){b=v;++actualContinuation;check(gross[g].incomingMomentum==0,"F17 continuation must not freeze old/guessed inside U");}
                        else{b=c==0?.03125L:.09375L;ambientDifferent+=v!=b;check(gross[g].incomingMomentum==b,"F17 actual nonzero ambient component");}}
                    check(std::find(frozen.saddle().freeVelocityDofs().begin(),frozen.saddle().freeVelocityDofs().end(),s.primalVelocity)!=frozen.saddle().freeVelocityDofs().end(),
                        "F17 pressure normal was accidentally fixed to ambient momentum");
                }else if(pin>0){b=c==0?.015625L:flow.normalVelocity[s.primalVelocity];
                    check(gross[g].incomingMomentum==b&&typed.beginVelocity()[g].x==typed.endVelocity()[g].x
                        &&typed.beginVelocity()[g].y==typed.endVelocity()[g].y,"F17 actual steady full-vector prescribed feed");}
                const auto n=pout*v-pin*b,p=.5L*pout*v*v-.5L*pin*b*b,d=.5L*pin*(v-b)*(v-b),w=.5L*(pout-pin-q)*v*v;
                momentum[c]+=n;physical+=p;dissipation+=d;representation+=w;boundaryWork+=v*n-.5L*q*v*v;
                diagonal+=(pout-(freePort[s.primalVelocity]&&model==Model::InteriorContinuationZeroShear?pin:0))*v*v;exterior[owner]+=q;
                boundaryScale+=std::abs(n)+std::abs(p)+std::abs(d)+std::abs(w)+std::abs(v*n)+std::abs(q*v*v);
            }
            for(int u=0;u<layout.velocityCount();++u)diagonal-=.5L*exterior[u]*um.normalVelocity()[u]*um.normalVelocity()[u];
            check(boundaryScale>0&&freeIncoming>0,"F17 did not exercise real free-port backflow work");
            for(int c=0;c<2;++c){near(report.boundaryMomentumRate[c],momentum[c],boundaryScale,"F17 independent gross boundary momentum",1e-12L);
                check(reassessed.report().components[c].momentumIdentity<=1e-12L&&reassessed.report().components[c].kineticIdentity<=1e-12L,"F17 original F16 momentum/kinetic identities");}
            near(report.boundaryPhysicalKineticWork,dt*physical,dt*boundaryScale,"F17 physical gross KE",1e-12L);
            near(report.inletDissipationWork,dt*dissipation,dt*boundaryScale,"F17 real inlet upwind dissipation",1e-12L);
            near(report.boundaryFluxRoundingWork,dt*representation,dt*boundaryScale,"F17 signed canonical/gross work split",1e-12L);
            near(report.boundaryKineticWork,dt*boundaryWork,dt*boundaryScale,"F17 P+D+Wrep counted once",1e-12L);
            near(boundaryWork,physical+dissipation+representation,boundaryScale,"F17 independent boundary decomposition",1e-12L);
            std::vector<long double> traction(layout.velocityCount());
            for(const auto& port:flow.pressureOutlets){const auto f=native.face(*mesh,layout.nativeFace(port.velocityDof));traction[port.velocityDof]=-static_cast<long double>(f.outwardSign)*f.area*port.pressure;
                near(frozen.saddle().pressureOutletLoads()[port.velocityDof],traction[port.velocityDof],std::abs(traction[port.velocityDof]),"F17 original pressure traction",1e-12L);}
            check(frozen.saddle().pressureGaugeCount()==0,"F17 must retain physical pressure reference");
            const auto cp=frozen.pressureLoad().apply(*mesh,endpoint);std::vector<bool> free(layout.velocityCount());for(int u:frozen.saddle().freeVelocityDofs())free[u]=true;
            long double kinetic=0,pwork=0,outwork=0,reaction=0,convection=0,massCross=0,scale=0,momentumError=0,momentumScale=0,rawOperator=0,operatorScale=0;
            ConservativeVofMomentumContract::qualify(*mesh,dual);
            check(report.consistencyPower==0,"conservative actual Step applied a diagnostic uR power");
            const auto impulse=dual.consistencyImpulse(*mesh,um,sources,clock);
            for(int u=0;u<layout.velocityCount();++u){const long double v=um.normalVelocity()[u],a=oldU[u],b=endpoint.normalVelocity()[u],m0=dual.beginMass().masses()[u],m1=dual.endMass().masses()[u];
                const auto rate=(m1*b-m0*a)/dt,load=reassessed.momentumLoads()[u],residual=rate+cp[u]-load-traction[u];
                check(load==-reassessed.momentumRates()[u],"F17 actual conservative Step load retained uR");
                const auto rowScale=std::abs(rate)+std::abs(static_cast<long double>(cp[u]))+std::abs(load)+std::abs(traction[u]);
                if(free[u]){momentumError=std::hypot(momentumError,residual);momentumScale=std::hypot(momentumScale,rowScale);}else reaction+=dt*v*residual;
                const auto k=.5L*(m1*b*b-m0*a*a),pw=dt*v*cp[u],ow=dt*v*traction[u],cw=dt*v*load,mc=.5L*(m1-m0)*a*b;
                kinetic+=k;pwork+=pw;outwork+=ow;convection+=cw;massCross+=mc;scale+=std::abs(k)+std::abs(pw)+std::abs(ow)+std::abs(cw)+std::abs(mc)+dt*std::abs(v)*rowScale;
                const auto op=(.5L*dual.cells()[u].outwardIntegral-static_cast<long double>(dt)*frozen.frozenConvectionMatrix().coeff(u,u))*v*v;
                rawOperator+=op;operatorScale+=std::abs(op);
                const auto expected=LAMP::Geometry::Detail::FiniteVolume::Dyadic::product({v,dual.cells()[u].integratedDefect}).rounded();
                if(expected==0)check(impulse[u]==0,"F17 zero actual uR must remain zero");else near(impulse[u],expected,std::abs(expected),"F17 actual uR own scale",1e-12L);++uRChecks;
            }
            check(scale>0&&momentumScale>0&&operatorScale+std::abs(dt*diagonal)>0,"F17 independent physical references inactive");
            near(momentumError,0,momentumScale,"F17 physical free momentum including pressure traction",1e-11L);
            const auto energy=kinetic+pwork-convection-outwork-reaction+massCross;near(energy,0,scale,"F17 independent endpoint KE/mass/pressure/reaction closure",1e-11L);
            near(report.kineticChange,kinetic,scale,"F17 native endpoint kinetic change",1e-12L);
            near(report.pressureOutletWork,outwork,scale,"F17 pressure work exactly once",1e-12L);near(report.prescribedReactionWork,reaction,scale,"F17 fixed feed reaction",1e-12L);
            near(report.energy.operatorRateRoundoffWork,rawOperator+dt*diagonal,operatorScale+std::abs(dt*diagonal),"F17 actual exterior diagonal replaces legacy abs-net correction",1e-12L);
            check(report.momentumPowerRelativeDefect<=config.momentumPowerTolerance&&std::isfinite(report.energy.closureRelativeResidual)
                &&report.energy.closureRelativeResidual<=1e-11L,"F17 original momentum/energy closure gates");
            for(bool half:{true,false}){const auto& ledger=half?volume.report().boundaryHalf:volume.report().boundaryFull;const auto& now=half?volume.midpoint():volume.endpoint();
                LAMP::Geometry::Detail::FiniteVolume::Dyadic change[2],out[2],in[2];long double capacitySum=0;
                for(std::size_t ci=0;ci<stocks.size();++ci){capacitySum+=source.cells()[ci].capacity;for(unsigned phase=0;phase<2;++phase){change[phase].add(now.plicCells()[ci].inventory[phase]);change[phase].add(stocks[ci].inventory[phase],-1);}}
                for(const auto& f:volume.faces())if(f.negativeCell<0||f.positiveCell<0){const auto& q=half?f.half:f.full;const int sign=f.negativeCell>=0?1:-1;
                    for(unsigned phase=0;phase<2;++phase){const auto p=phase==0?q.positiveLiquid:q.positiveGas,n=phase==0?q.negativeLiquid:q.negativeGas;out[phase].add(sign>0?p:n);in[phase].add(sign>0?n:p);}}
                for(unsigned phase=0;phase<2;++phase){change[phase].add(out[phase]);change[phase].add(in[phase],-1);near(change[phase].rounded(),0,capacitySum,"F17 real gross physical phase balance",1e-12L);++phaseChecks;}
                check(ledger.enabled&&ledger.sourceIdentity==id&&ledger.physicalRelativeConservation<=1e-12L&&ledger.algebraicRelativeConservation<=1e-12L,"F17 original phase/WY/rounding gates");}
            check(accepted.sameIdentity(old)&&source.identity()==id&&source.plicCells().data()==storage&&accepted.values().normalVelocity()==oldU&&accepted.values().pressure()==oldP,
                "F17 prepare or rejected prefix changed accepted source/U/p");
            for(std::size_t ci=0;ci<stocks.size();++ci)check(source.plicCells()[ci].inventory==stocks[ci].inventory,"F17 prepare altered raw accepted stocks");
            peakMomentum=std::max(peakMomentum,momentumError/momentumScale);peakEnergy=std::max(peakEnergy,std::abs(energy)/scale);
            peakBoundary=std::max(peakBoundary,std::abs(boundaryWork-physical-dissipation-representation)/boundaryScale);
            step.commit(*mesh,accepted);check(sameClock(accepted.timeLayer(),end)&&accepted.vofPrescribedFlow()->pressureOutletBackflowModel==model,"F17 accepted actual Step/model clock");++acceptedSteps;
        }
        if(model==Model::InteriorContinuationZeroShear)check(queries->empty(),"F17 continuation unexpectedly froze ambient/inside callback data");
        else{check(queries->size()==totalStages*flow.pressureOutlets.size(),"F17 ambient not prepared once per actual stage/native pressure port");
            for(const auto& query:*queries){const int u=layout.velocityDof(query.row);check(u>=0&&freePort[u]&&query.clock.position==Position::Midpoint
                &&query.clock.step>=1&&query.clock.step<=2,"F17 ambient query outside actual fullMid pressure owner");}}
        for(std::size_t ci=0;ci<initialStocks.size();++ci)check(initial.plicCells()[ci].inventory==initialStocks[ci].inventory,"F17 changed original source across commits");++scenes;
    }
    check(scenes==4&&acceptedSteps==8&&freeIncoming>0&&zeroInitialPorts>0&&actualContinuation>0&&ambientDifferent>0&&phaseChecks==32&&uRChecks>0,
        "F17 missing actual free-backflow/zero/ambient/current-U/full-Step evidence");
    std::cout<<std::setprecision(17)<<"LOCAL_PLIC_OPEN_MOMENTUM_STEP_PASS,scenes="<<scenes<<",accepted="<<acceptedSteps<<",guards="<<guardChecks
        <<",free_incoming="<<freeIncoming<<",initial_zero_ports="<<zeroInitialPorts<<",actual_U_continuation="<<actualContinuation
        <<",ambient_differs_from_unknown="<<ambientDifferent<<",phase_checks="<<phaseChecks<<",uR_checks="<<uRChecks
        <<",momentum_peak="<<peakMomentum<<",energy_peak="<<peakEnergy<<",boundary_split_peak="<<peakBoundary
        <<",EXPLICIT_STEP_ONLY_NO_APP_NO_NOZZLE_NO_CONTACT_PIN_CLAIM\n";
}
}
// END F17 ACTUAL OPT-IN FREE-PRESSURE MOMENTUM STEP CONTRACT

// BEGIN F18 PRODUCER-UNRESOLVED P0 BOUNDARY GROSS CONTRACT
namespace F18P0BoundaryGrossContract
{
using Exact=F15BoundaryGrossContract::Exact;
using Boundary=Slab::Boundary;
using Mode=Boundary::Mode;
unsigned scenes=0,grossChecks=0,minorOutChecks=0,zeroNetFaces=0,phaseLedgers=0;
long double peakP0=0,peakPhase=0,peakAlgebraic=0;
void inspect(const QuadtreeMesh& mesh,const Publication& source,const Transport& result,
    unsigned minor,double a,double b,const Clock& end)
{
    result.validateFor(mesh,source,result.advector(),100,1);
    const auto& report=result.report();
    check(report.localUnresolvedDonors>0,"F18 must consume actual producer-unresolved P0 donors");
    check(report.maximumIncomingCfl>0&&report.maximumOutflowCfl>0&&report.maximumLocalColumnCfl>0,
        "F18 actual gross donor/receiver/column Courant cannot vanish with net cancellation");
    for(bool half:{true,false}){
        const auto& candidate=half?result.midpoint():result.endpoint();
        const auto& ledger=half?report.boundaryHalf:report.boundaryFull;
        check(ledger.enabled&&ledger.sourceIdentity==source.identity()&&sameClock(ledger.sourceTime,source.timeLayer())
            &&sameClock(ledger.targetTime,candidate.timeLayer()),"F18 actual owner/fullclock half/full receipt");
        std::array<Exact,2> before,after,in,out,netRound;
        for(std::size_t ci=0;ci<source.plicCells().size();++ci){
            const auto& old=source.plicCells()[ci];const auto& now=candidate.plicCells()[ci];
            check(!old.resolved&&old.inventory[0]>0&&old.inventory[1]>0,
                "F18 public producer must supply the existing double-positive unresolved P0 classification");
            check(now.inventory[minor]>0,"F18 may not erase an independently positive microphase");
            Exact defect(now.inventory[0]);defect.add(now.inventory[1]);defect.add(candidate.cells()[ci].capacity,-1);
            check(defect.compare(candidate.capacityDefects()[ci])==0,"F18 exact certified closure must follow actual updated stocks");
            for(unsigned p=0;p<2;++p){before[p].add(old.inventory[p]);after[p].add(now.inventory[p]);}
        }
        for(const auto& f:result.faces()){
            if(f.negativeCell>=0&&f.positiveCell>=0)continue;
            const auto native=source.domain().face(mesh,f.nativeRow);
            const auto& flux=half?f.half:f.full;
            const auto expected=native.horizontal&&native.area>0?
                F15BoundaryGrossContract::analyticGross(a,b,end.end-end.begin,native.area,half?.5L:1.L):std::array<long double,2>{};
            const auto temporalScale=expected[0]+expected[1];
            if(temporalScale>0){
                near(flux.positiveTotal,expected[0],temporalScale,"F18 independent actual positive temporal gross",1e-12L);
                near(flux.negativeTotal,expected[1],temporalScale,"F18 independent actual negative temporal gross",1e-12L);
                const bool insideNegative=f.negativeCell>=0;const int inside=insideNegative?f.negativeCell:f.positiveCell;
                const auto outgoing=insideNegative?flux.positiveTotal:flux.negativeTotal;
                const auto incoming=insideNegative?flux.negativeTotal:flux.positiveTotal;
                const auto& old=source.plicCells()[inside];
                const auto expectedMinor=outgoing*(old.inventory[minor]/source.cells()[inside].capacity);
                const auto expectedMajor=outgoing-expectedMinor;
                const auto positive=std::array<long double,2>{flux.positiveLiquid,flux.positiveGas};
                const auto negative=std::array<long double,2>{flux.negativeLiquid,flux.negativeGas};
                const auto& actualOut=insideNegative?positive:negative;const auto& actualIn=insideNegative?negative:positive;
                check(flux.limiter==1,"F18 low-CFL fixture must not rely on donor correction");
                if(outgoing>0){check(expectedMinor>0&&actualOut[minor]>0,"F18 nonzero microphase outgoing must survive P0 splitting");
                    near(actualOut[minor],expectedMinor,expectedMinor,"F18 independent original P0 minor ratio",1e-12L);
                    near(actualOut[1-minor],expectedMajor,expectedMajor,"F18 independent P0 major outgoing",1e-12L);
                    peakP0=std::max(peakP0,std::abs(actualOut[minor]-expectedMinor)/expectedMinor);++minorOutChecks;}
                check(actualIn[minor]==0&&actualIn[1-minor]==incoming,"F18 only declared pure-major reservoir may supply incoming gross");
                check(flux.liquid==flux.positiveLiquid-flux.negativeLiquid&&flux.gas==flux.positiveGas-flux.negativeGas,
                    "F18 phase canonical net must be one native difference of its independent gross branches");
                if(!half&&a==-b){check(flux.total==0&&flux.positiveTotal>0&&flux.negativeTotal>0,
                    "F18 net-zero total must retain positive incoming and outgoing gross");++zeroNetFaces;}
                ++grossChecks;
            }else check(flux.positiveTotal==0&&flux.negativeTotal==0,"F18 closed radial/wall trace must stay zero");
            const int sign=f.negativeCell>=0?1:-1;
            for(unsigned p=0;p<2;++p){const auto positive=p==0?flux.positiveLiquid:flux.positiveGas;
                const auto negative=p==0?flux.negativeLiquid:flux.negativeGas,net=p==0?flux.liquid:flux.gas;
                out[p].add(sign>0?positive:negative);in[p].add(sign>0?negative:positive);
                Exact error(net);error.add(positive,-1);error.add(negative);netRound[p].add(error,-sign);}
        }
        for(unsigned p=0;p<2;++p){
            check(ledger.source[p]==before[p].rounded()&&ledger.target[p]==after[p].rounded()
                &&ledger.inflow[p]==in[p].rounded()&&ledger.outflow[p]==out[p].rounded(),"F18 independent physical phase boundary ledger");
            check(ledger.netFluxRounding[p]==netRound[p].rounded(),"F18 per-phase canonical-net rounding receipt");
            auto physical=after[p];physical.add(before[p],-1);physical.add(out[p]);physical.add(in[p],-1);
            check(ledger.physicalResidual[p]==physical.rounded(),"F18 physical boundary residual is not a WY source");
            auto algebraic=physical;algebraic.add(ledger.compression[p],-1);algebraic.add(ledger.roundingError[p],-1);
            algebraic.add(netRound[p],-1);
            // Report compression/stock-rounding totals are native scalars;
            // their final contraction rounding is not silently an exact oracle.
            const auto phaseScale=before[p].absolute();auto actualScale=phaseScale;actualScale.add(in[p].absolute());actualScale.add(out[p].absolute());
            check(actualScale.compare(Exact{})>0,"F18 actual phase scale must be positive without a capacity floor");
            check(physical.absolute().compare(actualScale.multiplied(1e-12L))<=0,
                "F18 raw phase-scale physical conservation (including the tiny phase)");
            check(algebraic.absolute().compare(actualScale.multiplied(1e-12L))<=0,
                "F18 phase-scale reported WY/RN/canonical-net algebraic conservation");
            check(std::abs(ledger.algebraicResidual[p])<=1e-12L*actualScale.rounded(),"F18 original exact owner algebraic residual");
            peakPhase=std::max(peakPhase,std::abs(ledger.physicalResidual[p])/actualScale.rounded());
            peakAlgebraic=std::max(peakAlgebraic,std::abs(ledger.algebraicResidual[p])/actualScale.rounded());++phaseLedgers;
        }
        check(ledger.physicalRelativeConservation<=1e-12L&&ledger.algebraicRelativeConservation<=1e-12L
            &&ledger.massRelativeConservation<=1e-12L,"F18 original physical/algebraic/mass gates unchanged");
    }
}
void scene(bool rz,unsigned pattern,bool actualMagnitude=false)
{
    stage="F18 producer P0 actual boundary gross rz="+std::to_string(rz)+" pattern="+std::to_string(pattern)+" magnitude="+std::to_string(actualMagnitude);
    auto mesh=F15BoundaryGrossContract::makeMesh(rz);
    Clock sourceClock{};Clock end{pattern%2?2:1,pattern%2?.125:0,(pattern%2?.125:0)+1./32,Position::End};
    const std::array<std::array<double,2>,4> patterns{{{.125,-.125},{-.125,.125},{.125,-.375},{-.375,.125}}};
    double a=patterns[pattern][0],b=patterns[pattern][1];
    if(actualMagnitude){
        DictionaryMesh d;d.Nx=d.Ny=1;d.minCoord={4,60};d.maxCoord={8,64};d.maxCellNum=64;d.maxInitialIterations=0;d.maxLevel_local={1};
        mesh=std::make_unique<QuadtreeMesh>(d);sourceClock={44,.10814720710750886,.11095367219757830,Position::End};
        end={45,.11095367219757830,.11376013728764775,Position::End};a=5.3449021203856062e-6;b=-2.616277548365469e-6;
    }else if(pattern%2)sourceClock={1,0,.125,Position::End};
    const unsigned minor=actualMagnitude||rz?0:1;
    const auto native=domain(*mesh,rz);
    check(actualMagnitude||!native.uniformLeaves(),"F18 standard fixture must retain real 2:1 native leaves");
    const auto seed=Publication::prepareLocalPlic(*mesh,native,sourceClock,std::vector<double>(mesh->cellIndexNum(),minor==0?0:1));
    std::vector<std::array<long double,2>> inventory(seed.plicCells().size());
    for(std::size_t ci=0;ci<inventory.size();++ci){inventory[ci][1-minor]=seed.cells()[ci].capacity;
        inventory[ci][minor]=actualMagnitude?4.1058512395836025e-197L:std::ldexp(seed.cells()[ci].capacity,-80);}
    const auto source=Publication::prepareLocalPlicInventory(*mesh,native,sourceClock,inventory,{});
    const auto saved=source.plicCells();const auto* storage=source.plicCells().data();const auto id=source.identity();const auto alpha=source.alpha();
    for(const auto& cell:saved)check(!cell.resolved&&cell.inventory[0]>0&&cell.inventory[1]>0,
        "F18 witness classification must be issued by public producer, never a hand-set flag");
    const auto boundary=Boundary::prepare(*mesh,source,F15BoundaryGrossContract::rulesFor(*mesh,source,1-minor));
    const auto ua=F15BoundaryGrossContract::normals(*mesh,source,a),ub=F15BoundaryGrossContract::normals(*mesh,source,b);
    const auto slab=Slab::prepareBoundary(*mesh,source,boundary,end,ua,ub);
    const auto result=Transport::prepareLocalPlicBoundarySweep(*mesh,source,slab,100,1);
    inspect(*mesh,source,result,minor,a,b,end);
    check(source.identity()==id&&source.plicCells().data()==storage&&source.alpha()==alpha,"F18 may not replace accepted source identity/storage/alpha");
    for(std::size_t ci=0;ci<saved.size();++ci)check(source.plicCells()[ci].inventory==saved[ci].inventory
        &&source.plicCells()[ci].normal.x==saved[ci].normal.x&&source.plicCells()[ci].normal.y==saved[ci].normal.y
        &&source.plicCells()[ci].offset==saved[ci].offset&&source.plicCells()[ci].resolved==saved[ci].resolved,
        "F18 candidate altered accepted independent stocks or geometry");
    if(actualMagnitude)std::cout<<std::setprecision(20)<<"F18_ACTUAL_MAGNITUDE_MECHANISM,source_L="<<source.plicCells()[0].inventory[0]
        <<",source_G="<<source.plicCells()[0].inventory[1]<<",source_resolved="<<source.plicCells()[0].resolved
        <<",FRESH_PUBLIC_PRODUCER_NOT_SOURCE397_RESTART_NOT_ACTUAL_DIRECTIONAL_SCRATCH\n";
    ++scenes;
}
void run()
{
    for(bool rz:{false,true})for(unsigned pattern=0;pattern<4;++pattern)scene(rz,pattern);
    scene(true,0,true);
    check(scenes==9&&grossChecks>0&&minorOutChecks>0&&zeroNetFaces>0&&phaseLedgers==36,"F18 necessary P0/gross/half/full witnesses incomplete");
    std::cout<<std::setprecision(17)<<"LOCAL_PLIC_P0_BOUNDARY_GROSS_PASS,scenes="<<scenes<<",phase_ledgers="<<phaseLedgers
        <<",gross_checks="<<grossChecks<<",nonzero_minor_out="<<minorOutChecks<<",net0_gross_faces="<<zeroNetFaces
        <<",p0_peak="<<peakP0<<",phase_scale_peak="<<peakPhase<<",algebraic_peak="<<peakAlgebraic
        <<",PRODUCER_UNRESOLVED_P0_ONLY_NO_RESOLVED_MIXED_CHRONOLOGY_NO_APP_REPLAY_CLAIM\n";
}
}
// END F18 PRODUCER-UNRESOLVED P0 BOUNDARY GROSS CONTRACT


namespace {
namespace TransportValidationBatchContract {
using Boundary=Slab::Boundary;using Mode=Boundary::Mode;
unsigned cases=0,rejections=0,changedMeshes=0,acceptedCopies=0;std::size_t physicalFaces=0;

// P0-2 compares against the original public source/domain order. These literal
// references deliberately use only public getters and never the batch API.
unsigned literalValidationPairs=0,literalValidationFailures=0,latePreparePairs=0;
void literalBoundary(const QuadtreeMesh& mesh,const Publication& source,const Boundary& boundary) {
    const auto owner=boundary.sourceIdentity();const auto clock=boundary.sourceTime();
    source.validateFor(mesh,owner,clock);boundary.domain().validateMeshIdentity(mesh);
}
void literalSlab(const QuadtreeMesh& mesh,const Publication& source,const Slab& slab) {
    const auto owner=slab.sourceIdentity();const auto clock=slab.sourceTime();
    source.validateFor(mesh,owner,clock);slab.domain().validateMeshIdentity(mesh);
    if(slab.boundary())literalBoundary(mesh,source,*slab.boundary());
}
struct ValidationNative {
    std::fenv_t native{};int error=0;
#if defined(_MSC_VER) && defined(_M_X64)
    unsigned csr=0;
#endif
    void capture(){error=errno;
#if defined(_MSC_VER) && defined(_M_X64)
        csr=_mm_getcsr();
#endif
        check(std::fegetenv(&native)==0,"P0-2 validation native FP capture");}
    void restore()const{std::fesetenv(&native);
#if defined(_MSC_VER) && defined(_M_X64)
        _mm_setcsr(csr);
#endif
        errno=error;}
    bool same(const ValidationNative& other)const{return error==other.error
#if defined(_MSC_VER) && defined(_M_X64)
        &&csr==other.csr
#endif
        &&std::memcmp(&native,&other.native,sizeof native)==0;}
};
template<class Literal,class Actual>void validationPair(Literal literal,Actual actual,const char* expected=nullptr) {
    ValidationNative saved;saved.capture();struct Restore {const ValidationNative& state;~Restore(){state.restore();}} restore{saved};
    check(std::feraiseexcept(FE_INEXACT)==0,"P0-2 validation prime native flags");errno=173;
    ValidationNative entry;entry.capture();
    struct Outcome {ValidationNative exit;bool failed=false;std::string what,type;};
    const auto observe=[&](auto action){entry.restore();Outcome result;
        try{action();result.exit.capture();}
        catch(const std::exception& error){result.exit.capture();result.failed=true;result.what=error.what();result.type=typeid(error).name();}
        return result;};
    const auto reference=observe(literal),candidate=observe(actual);
    check(reference.failed==candidate.failed&&reference.type==candidate.type&&reference.what==candidate.what,
        "P0-2 old literal validator exception type/first what changed");
    check(reference.exit.same(candidate.exit),"P0-2 old literal validator complete FP/rawCSR/errno changed");
    if(expected)check(reference.failed&&reference.what.find(expected)!=std::string::npos,
        "P0-2 competing later guard displaced original first error");
    ++literalValidationPairs;if(reference.failed)++literalValidationFailures;
}
void pairedValidators(const QuadtreeMesh& mesh,const Publication& source,const Boundary& boundary,const Slab& slab) {
    validationPair([&]{literalBoundary(mesh,source,boundary);},[&]{boundary.validateFor(mesh,source);});
    validationPair([&]{literalSlab(mesh,source,slab);},[&]{slab.validateFor(mesh,source);});
}
void competingPrepareErrors(QuadtreeMesh& mesh,const Publication& source,const Publication& foreign,
    const Boundary&,const Slab& slab,bool open) {
    auto& leaf=mesh.cells(mesh.cellIndexToLabel[0]);const auto center=leaf.centerCoord();
    struct Restore {decltype(leaf) cell;Vector2D center;~Restore(){cell.setCenterCoord(center);}} restore{leaf,center};
    leaf.setCenterCoord(center+Vector2D{.00390625,0});
    const auto original=[&](const Publication& input){
        if(open)literalBoundary(mesh,input,*slab.boundary());literalSlab(mesh,input,slab);};
    const auto actual=[&](const Publication& input){
        if(open)(void)Transport::prepareLocalPlicBoundarySweep(mesh,input,slab,-1,1);
        else (void)Transport::prepareLocalPlicSweep(mesh,input,slab,-1,1);};
    validationPair([&]{original(source);},[&]{actual(source);},"mesh geometry/topology changed since identity construction");++latePreparePairs;
    validationPair([&]{original(foreign);},[&]{actual(foreign);},"foreign VOF owner");++latePreparePairs;
}

template<class F>void invalid(F action,const char* reason) {
    bool caught=false;try{action();}catch(const std::invalid_argument&){caught=true;}
    check(caught,std::string("private mesh validation batch bypassed public guard: ")+reason);++rejections;
}
void validate(const QuadtreeMesh& mesh,const Publication& source,const Boundary& boundary,const Slab& slab,const Transport& transport) {
    source.validateFor(mesh,source.identity(),source.timeLayer());boundary.validateFor(mesh,source);slab.validateFor(mesh,source);
    transport.validateFor(mesh,source,slab,100,1);
}
void invalidMesh(const QuadtreeMesh& mesh,const Publication& source,const Boundary& boundary,const Slab& slab,const Transport& transport) {
    pairedValidators(mesh,source,boundary,slab);
    invalid([&]{source.validateFor(mesh,source.identity(),source.timeLayer());},"standalone Publication live identity");
    invalid([&]{boundary.validateFor(mesh,source);},"standalone Boundary live identity");
    invalid([&]{slab.validateFor(mesh,source);},"standalone Slab live identity");
    invalid([&]{transport.validateFor(mesh,source,slab,100,1);},"Transport after an earlier successful complete validation");++changedMeshes;
}
void mutations(QuadtreeMesh& mesh,const Publication& source,const Boundary& boundary,const Slab& slab,const Transport& transport) {
    validate(mesh,source,boundary,slab,transport);
    auto& leaf=mesh.cells(mesh.cellIndexToLabel[0]);const auto center=leaf.centerCoord();
    leaf.setCenterCoord(center+Vector2D{.00390625,0});invalidMesh(mesh,source,boundary,slab,transport);leaf.setCenterCoord(center);
    validate(mesh,source,boundary,slab,transport);
    int parent=-1;for(int i=0;i<mesh.cellIndexNum()&&parent<0;++i){const int candidate=mesh.cells(mesh.cellIndexToLabel[i]).parent;
        if(candidate>=0&&std::find(mesh.cellIndexToLabel,mesh.cellIndexToLabel+mesh.cellIndexNum(),candidate)==mesh.cellIndexToLabel+mesh.cellIndexNum())parent=candidate;}
    check(parent>=0,"validation batch ancestor mutation lacks an actual nonindexed ancestor");
    auto& ancestor=mesh.cells(parent);check(ancestor.son[0]>=0&&ancestor.son[1]>=0&&ancestor.son[0]!=ancestor.son[1],"ancestor-only mutation has no two actual distinct children");
    std::swap(ancestor.son[0],ancestor.son[1]);invalidMesh(mesh,source,boundary,slab,transport);std::swap(ancestor.son[0],ancestor.son[1]);
    validate(mesh,source,boundary,slab,transport);
    auto& edge=mesh.vEdges(mesh.vEdgeIndexToLabel[0]);check(edge.edgeToVert[0]!=edge.edgeToVert[1],"connectivity guard has a degenerate original edge");
    std::swap(edge.edgeToVert[0],edge.edgeToVert[1]);invalidMesh(mesh,source,boundary,slab,transport);std::swap(edge.edgeToVert[0],edge.edgeToVert[1]);
    validate(mesh,source,boundary,slab,transport);
    auto& vertex=mesh.verts(mesh.vertIndexToLabel[0]);const auto coordinate=vertex.coord();
    vertex.setCoord(coordinate+Vector2D{0,.00390625});invalidMesh(mesh,source,boundary,slab,transport);vertex.setCoord(coordinate);
    validate(mesh,source,boundary,slab,transport);
}
void scopedLiveProofs(QuadtreeMesh& mesh,const Publication& source)
{
    using Batch=NativeMeshIdentity::ValidationBatch;
    const auto valid=[&]{source.validateFor(mesh,source.identity(),source.timeLayer());};
    const auto rejected=[&](auto action){bool caught=false;try{action();}catch(const std::invalid_argument&){caught=true;}
        check(caught,"scoped live proof concealed a foreign binding, callback mutation or stale clock");};
    auto foreign=mesh.clone();auto& other=foreign->cells(foreign->cellIndexToLabel[0]);
    other.setCenterCoord(other.centerCoord()+Vector2D{.00390625,0});
    auto& leaf=mesh.cells(mesh.cellIndexToLabel[0]);const auto center=leaf.centerCoord();
    struct Restore {decltype(leaf) cell;Vector2D center;~Restore(){cell.setCenterCoord(center);}} restore{leaf,center};
    {
        Batch outer(mesh);valid();
        NativeMeshIdentity retained(mesh);retained.validate(mesh);const auto words=retained.checkpointWords();
        auto& self=retained;retained=std::move(self);
        check(retained.checkpointWords()==words,"identity self move detached its checked geometry words");retained.validate(mesh);
        auto moved=std::move(retained);moved.validate(mesh);
        check(moved.checkpointWords()==words,"identity move changed checked geometry words");
        rejected([&]{source.validateFor(*foreign,source.identity(),source.timeLayer());});
        {Batch inner(mesh);valid();auto wrong=source.timeLayer();wrong.begin=std::nextafter(wrong.begin,1.);
            rejected([&]{source.validateFor(mesh,source.identity(),wrong);});
            // Recreate a different immutable identity at exactly the same address.
            // The strong lifetime token must prevent address-reuse proof hits.
            std::optional<NativeMeshIdentity> slot;slot.emplace(mesh);slot->validate(mesh);
            const auto* address=&*slot;slot.reset();slot.emplace(*foreign);
            check(address==&*slot,"optional identity did not reuse its actual storage");
            rejected([&]{slot->validate(mesh);});
            struct CallbackFailure {};
            try{Batch::Pause callback;leaf.setCenterCoord(center+Vector2D{.00390625,0});
                rejected(valid);throw CallbackFailure{};}
            catch(const CallbackFailure&){}
            rejected(valid); // Exception must resume an empty proof chain.
            leaf.setCenterCoord(center);valid();
        }
        valid();
        struct SegmentFailure {};
        try{Batch inner(mesh);valid();leaf.setCenterCoord(center+Vector2D{.00390625,0});throw SegmentFailure{};}
        catch(const SegmentFailure&){}
        rejected(valid); // Unwinding a nested segment clears outer proofs too.
        leaf.setCenterCoord(center);valid();
    }
    // A public call after the lexical batch must perform its own live scan.
    leaf.setCenterCoord(center+Vector2D{.00390625,0});rejected(valid);
    leaf.setCenterCoord(center);valid();
}
void scene(bool rz,bool open) {
    stage="private transport validation batch "+std::to_string(rz)+","+std::to_string(open);
    auto mesh=F15BoundaryGrossContract::makeMesh(rz);const auto native=domain(*mesh,rz);
    check(!native.uniformLeaves(),"validation batch fixture lost real nonuniform 2:1 leaves");
    const Clock oldClock{7,.125,.25,Position::End};
    const auto source=Publication::prepareLocalPlic(*mesh,native,oldClock,std::vector<double>(mesh->cellIndexNum(),0.));
    const auto saved=source.alpha();const auto stocks=source.plicCells();const auto* storage=source.plicCells().data();
    auto rules=F15BoundaryGrossContract::rulesFor(*mesh,source,1);
    if(!open)for(auto& rule:rules)if(rule!=Mode::Interior)rule=Mode::Closed;
    const auto boundary=Boundary::prepare(*mesh,source,rules);const auto a=open?F15BoundaryGrossContract::normals(*mesh,source,.125):vofMaterialNormals(*mesh,source,0);
    const Clock next{8,.25,.28125,Position::End};
    const auto slab=open?Slab::prepareBoundary(*mesh,source,boundary,next,a,a):Slab::prepare(*mesh,source,next,a,a);
    const auto transport=open?Transport::prepareLocalPlicBoundarySweep(*mesh,source,slab,100,1):Transport::prepareLocalPlicSweep(*mesh,source,slab,100,1);
    validate(*mesh,source,boundary,slab,transport);pairedValidators(*mesh,source,boundary,slab);
    check(open?(slab.boundary()!=nullptr):(slab.boundary()==nullptr),"open/closed fixture lost original boundary owner distinction");
    bool moving=false;for(std::size_t row=0;row<transport.faces().size();++row){const auto f=native.face(*mesh,static_cast<int>(row));
        for(bool half:{true,false}){const auto& flow=half?transport.faces()[row].half:transport.faces()[row].full;
            const auto expected=static_cast<long double>(a[row])*f.area*(next.end-next.begin)*(half?.5L:1.L);
            if(expected==0)check(flow.total==0,"validation batch produced flux on exact zero boundary/axis");
            else {near(flow.total,expected,std::abs(expected),"independent steady normal swept volume",1e-12L);moving=true;}
            check(flow.liquid==0&&flow.gas==flow.total&&flow.mass==flow.gas,"validation-only change altered pure-gas physical/material flux");++physicalFaces;}}
    check(moving,"validation batch fixture has no actual nonzero transport");
    for(const auto* candidate:{&transport.midpoint(),&transport.endpoint()})for(std::size_t i=0;i<stocks.size();++i)
        check(candidate->plicCells()[i].inventory==stocks[i].inventory,"solenoidal strict pure gas changed independent phase inventory");
    if(open)check(transport.report().boundaryFull.enabled&&transport.report().boundaryFull.physicalRelativeConservation<=1e-12L
        &&transport.report().boundaryFull.algebraicRelativeConservation<=1e-12L,"open source-bound physical/algebraic budget changed");
    scopedLiveProofs(*mesh,source);
    auto equivalent=mesh->clone();validate(*equivalent,source,boundary,slab,transport);
    pairedValidators(*equivalent,source,boundary,slab);++acceptedCopies;
    // The same source domain is intentionally accepted on a content-identical
    // mesh copy. C++ object-address inequality is not a foreign topology.
    mutations(*mesh,source,boundary,slab,transport);mutations(*equivalent,source,boundary,slab,transport);
    const auto foreign=Publication::prepareLocalPlic(*mesh,domain(*mesh,rz),oldClock,saved);
    invalid([&]{source.validateFor(*mesh,foreign.identity(),oldClock);},"foreign publication identity");
    invalid([&]{boundary.validateFor(*mesh,foreign);},"Boundary foreign publication owner");
    invalid([&]{slab.validateFor(*mesh,foreign);},"Slab foreign publication owner");
    invalid([&]{transport.validateFor(*mesh,foreign,slab,100,1);},"Transport foreign publication owner/domain");
    pairedValidators(*mesh,foreign,boundary,slab);
    competingPrepareErrors(*mesh,source,foreign,boundary,slab,open);
    const auto sameValues=open?Slab::prepareBoundary(*mesh,source,boundary,next,a,a):Slab::prepare(*mesh,source,next,a,a);
    invalid([&]{transport.validateFor(*mesh,source,sameValues,100,1);},"same-value foreign actual slab");
    auto wrong=oldClock;wrong.begin=std::nextafter(wrong.begin,0.);
    invalid([&]{source.validateFor(*mesh,source.identity(),wrong);},"same endpoint but different complete source clock");
    invalid([&]{transport.validateFor(*mesh,source,slab,101,1);},"different actual liquid density");
    invalid([&]{transport.validateFor(*mesh,source,slab,100,2);},"different actual gas density");
    validate(*mesh,source,boundary,slab,transport);
    check(source.alpha()==saved&&source.plicCells().data()==storage&&sameClock(source.timeLayer(),oldClock),"validation/rejection replaced source storage or clock");
    for(std::size_t i=0;i<stocks.size();++i)check(source.plicCells()[i].inventory==stocks[i].inventory,"validation/rejection modified raw source stock");++cases;
    std::cout<<"TRANSPORT_VALIDATION_BATCH,"<<(rz?"RZ":"XY")<<','<<(open?"open":"closed")<<",cells,"<<mesh->cellIndexNum()
        <<",content_identical_clone,accepted,leaf_ancestor_connectivity_coordinate,independent_public_rejections\n";
}
void run(){for(bool rz:{false,true})for(bool open:{false,true})scene(rz,open);
    check(cases==4&&acceptedCopies==4&&changedMeshes==32&&rejections==160,"validation batch small fixture/guard coverage incomplete");
    check(literalValidationPairs==96&&literalValidationFailures==80&&latePreparePairs==8,
        "P0-2 independent literal/FP/late-first-error coverage incomplete");
    std::cout<<"TRANSPORT_PREPARE_VALIDATION_BATCH,literal_FP_pairs,"<<literalValidationPairs
        <<",rejected_pairs,"<<literalValidationFailures<<",competing_prepare_errors,"<<latePreparePairs
        <<",same_actual_identity_only,no_scan_counter_or_speed_claim\n";
    std::cout<<"TRANSPORT_VALIDATION_BATCH_PASS,cases,"<<cases<<",same_content_mesh_copies,"<<acceptedCopies<<",mutations,"<<changedMeshes
        <<",explicit_rejections,"<<rejections<<",independent_physical_face_checks,"<<physicalFaces
        <<",different_identity_binding_private_fallback_static_only,no_scan_count_or_speed_claim\n";
}
}
}
// BEGIN TASK51 LITERAL/PARALLEL FACE VALUES AND ORDERED CELL GATHER
namespace F51LocalFaceParallelContract
{
using Exact=LAMP::Geometry::Detail::FiniteVolume::Dyadic;
using Boundary=Slab::Boundary;
using Mode=Boundary::Mode;
unsigned pairs=0,parallelDirections=0,firstErrorPairs=0,zeroNetGross=0,physicalLedgers=0,pureRoundingWitnesses=0;
template<class T> void bits(T a,T b,const char* what){check(std::memcmp(&a,&b,sizeof(T))==0,what);}
void sameLedger(const GeometricVofBoundaryLedger& a,const GeometricVofBoundaryLedger& b)
{
    check(a.enabled==b.enabled&&a.sourceIdentity==b.sourceIdentity&&sameClock(a.sourceTime,b.sourceTime)
        &&sameClock(a.targetTime,b.targetTime),"Task51 boundary owner/fullclock");
    for(auto member:{&GeometricVofBoundaryLedger::source,&GeometricVofBoundaryLedger::target,
        &GeometricVofBoundaryLedger::inflow,&GeometricVofBoundaryLedger::outflow,
        &GeometricVofBoundaryLedger::compression,&GeometricVofBoundaryLedger::roundingError,
        &GeometricVofBoundaryLedger::netFluxRounding,&GeometricVofBoundaryLedger::physicalResidual,
        &GeometricVofBoundaryLedger::algebraicResidual})
        for(unsigned p=0;p<2;++p)bits((a.*member)[p],(b.*member)[p],"Task51 original phase ledger bits");
    for(auto member:{&GeometricVofBoundaryLedger::physicalRelativeConservation,
        &GeometricVofBoundaryLedger::algebraicRelativeConservation,&GeometricVofBoundaryLedger::massResidual,
        &GeometricVofBoundaryLedger::massRelativeConservation})bits(a.*member,b.*member,"Task51 original boundary aggregate bits");
}
void sameReport(const GeometricVofTransportReport& a,const GeometricVofTransportReport& b)
{
    for(auto member:{&GeometricVofTransportReport::donorQueries,&GeometricVofTransportReport::polygonIntegrals,
        &GeometricVofTransportReport::temporalIntervals,&GeometricVofTransportReport::limitedFaces,
        &GeometricVofTransportReport::temporalEvents,&GeometricVofTransportReport::localDirectionalSweeps,
        &GeometricVofTransportReport::localStripClips,&GeometricVofTransportReport::localUnresolvedDonors,
        &GeometricVofTransportReport::localPinnedSweeps,&GeometricVofTransportReport::localWholePhaseExhaustions,
        &GeometricVofTransportReport::localPhaseUpdateCalls,&GeometricVofTransportReport::localPhaseUpdateFallback,
        &GeometricVofTransportReport::localPhaseUpdateFast,&GeometricVofTransportReport::localPhaseUpdatePromoted,
        &GeometricVofTransportReport::localPhaseUpdateNoop,
        &GeometricVofTransportReport::localNativeP0Restrictions,&GeometricVofTransportReport::localInteriorP0Extensions})
        check(a.*member==b.*member,"Task51 original report count");
    check(a.maximumTemporalRefinements==b.maximumTemporalRefinements,"Task51 original temporal refinements");
    for(auto member:{&GeometricVofTransportReport::maximumLocalColumnCfl,&GeometricVofTransportReport::maximumIncomingCfl,
        &GeometricVofTransportReport::maximumCompressionResidual,&GeometricVofTransportReport::maximumLocalStripResidual,
        &GeometricVofTransportReport::maximumFaceTimeRelativeUncertainty,&GeometricVofTransportReport::maximumOutflowCfl,
        &GeometricVofTransportReport::maximumTotalDivergence,&GeometricVofTransportReport::maximumCoverageResidual,
        &GeometricVofTransportReport::maximumFluxCorrection,&GeometricVofTransportReport::maximumCapacityReconciliation,
        &GeometricVofTransportReport::maximumInventoryResidual,&GeometricVofTransportReport::halfRelativeConservation,
        &GeometricVofTransportReport::fullRelativeConservation})bits(a.*member,b.*member,"Task51 original report scalar bits");
    sameLedger(a.boundaryHalf,b.boundaryHalf);sameLedger(a.boundaryFull,b.boundaryFull);
}
void samePublication(const Publication& a,const Publication& b)
{
    check(sameClock(a.timeLayer(),b.timeLayer())&&a.geometryKind()==b.geometryKind(),"Task51 target geometry/fullclock");
    check(a.alpha().size()==b.alpha().size()&&a.cells().size()==b.cells().size()
        &&a.plicCells().size()==b.plicCells().size()&&a.capacityDefects().size()==b.capacityDefects().size(),"Task51 target shape");
    for(std::size_t i=0;i<a.alpha().size();++i)bits(a.alpha()[i],b.alpha()[i],"Task51 alpha bits");
    for(std::size_t i=0;i<a.cells().size();++i){const auto& x=a.cells()[i];const auto& y=b.cells()[i];
        check(x.nativeCell==y.nativeCell&&x.boundaryVertices==y.boundaryVertices&&x.triangles.size()==y.triangles.size(),
            "Task51 native cell geometry");
        for(auto member:{&GeometricVofCell::capacity,&GeometricVofCell::inventory,&GeometricVofCell::measured,
            &GeometricVofCell::centerScalar})bits(x.*member,y.*member,"Task51 native cell scalar bits");
        const auto& p=a.plicCells()[i];const auto& q=b.plicCells()[i];
        bits(p.rectangle.lower.x,q.rectangle.lower.x,"Task51 rectangle lower x");bits(p.rectangle.lower.y,q.rectangle.lower.y,"Task51 rectangle lower y");
        bits(p.rectangle.upper.x,q.rectangle.upper.x,"Task51 rectangle upper x");bits(p.rectangle.upper.y,q.rectangle.upper.y,"Task51 rectangle upper y");
        bits(p.normal.x,q.normal.x,"Task51 normal x");bits(p.normal.y,q.normal.y,"Task51 normal y");bits(p.offset,q.offset,"Task51 PLIC offset");
        check(p.resolved==q.resolved,"Task51 producer resolution");
        for(unsigned phase=0;phase<2;++phase){bits(p.inventory[phase],q.inventory[phase],"Task51 independent stock bits");
            bits(p.measured[phase],q.measured[phase],"Task51 measured geometry bits");}
        check(a.capacityDefects()[i].compare(b.capacityDefects()[i])==0,"Task51 exact capacity defect");
    }
    const auto& x=a.report();const auto& y=b.report();
    check(x.mixedCells==y.mixedCells&&x.pureLiquidCells==y.pureLiquidCells&&x.pureGasCells==y.pureGasCells
        &&x.zeroCollars==y.zeroCollars&&x.coincidentFaces==y.coincidentFaces
        &&x.maximumVolumeIterations==y.maximumVolumeIterations,"Task51 publication report counts");
    for(auto member:{&GeometricVofReport::maximumRelativeVolumeResidual,&GeometricVofReport::maximumRelativeCapacityResidual,
        &GeometricVofReport::inventory,&GeometricVofReport::measured})bits(x.*member,y.*member,"Task51 publication report bits");
}
void samePayload(const Transport& a,const Transport& b)
{
    check(a.source().identity()==b.source().identity()&&a.advector().sameInstance(b.advector()),"Task51 actual source/slab owner");
    check(a.faces().size()==b.faces().size(),"Task51 face shape");
    for(std::size_t i=0;i<a.faces().size();++i){const auto& x=a.faces()[i];const auto& y=b.faces()[i];
        check(x.nativeRow==y.nativeRow&&x.negativeCell==y.negativeCell&&x.positiveCell==y.positiveCell,"Task51 face incidence");
        check(F15BoundaryGrossContract::sameFields(x.half,y.half)&&F15BoundaryGrossContract::sameFields(x.full,y.full),
            "Task51 all twelve half/full native flux scalar bits/signs");}
    samePublication(a.midpoint(),b.midpoint());samePublication(a.endpoint(),b.endpoint());sameReport(a.report(),b.report());
}
void independent(const QuadtreeMesh& mesh,const Publication& source,const Transport& actual)
{
    actual.validateFor(mesh,source,actual.advector(),100,1);
    for(bool half:{true,false}){const auto& candidate=half?actual.midpoint():actual.endpoint();
        const auto& ledger=half?actual.report().boundaryHalf:actual.report().boundaryFull;
        std::array<Exact,2> before,after,in,out;Exact capacity;
        for(std::size_t i=0;i<source.cells().size();++i){capacity.add(source.cells()[i].capacity);
            Exact closure(candidate.plicCells()[i].inventory[0]);closure.add(candidate.plicCells()[i].inventory[1]);
            closure.add(candidate.cells()[i].capacity,-1);
            check(closure.compare(candidate.capacityDefects()[i])==0,"Task51 independent exact phase-pair/defect closure");
            for(unsigned p=0;p<2;++p){before[p].add(source.plicCells()[i].inventory[p]);after[p].add(candidate.plicCells()[i].inventory[p]);}}
        for(const auto& face:actual.faces()){const auto& flux=half?face.half:face.full;
            if(!half&&flux.total==0&&flux.positiveTotal>0&&flux.negativeTotal>0)++zeroNetGross;
            if(face.negativeCell>=0&&face.positiveCell>=0)continue;
            const bool insideNegative=face.negativeCell>=0;
            for(unsigned p=0;p<2;++p){const auto positive=p==0?flux.positiveLiquid:flux.positiveGas;
                const auto negative=p==0?flux.negativeLiquid:flux.negativeGas;
                out[p].add(insideNegative?positive:negative);in[p].add(insideNegative?negative:positive);}}
        for(unsigned p=0;p<2;++p){auto residual=after[p];residual.add(before[p],-1);residual.add(out[p]);residual.add(in[p],-1);
            if(ledger.enabled){bits(ledger.source[p],before[p].rounded(),"Task51 independent source boundary sum");
                bits(ledger.target[p],after[p].rounded(),"Task51 independent target boundary sum");
                bits(ledger.inflow[p],in[p].rounded(),"Task51 independent incoming gross");
                bits(ledger.outflow[p],out[p].rounded(),"Task51 independent outgoing gross");
                bits(ledger.physicalResidual[p],residual.rounded(),"Task51 independent physical residual");
                check(ledger.physicalRelativeConservation<=1e-12L&&ledger.algebraicRelativeConservation<=1e-12L
                    &&ledger.massRelativeConservation<=1e-12L,"Task51 original physical/algebraic/mass gates");++physicalLedgers;
            }else check(residual.absolute().compare(capacity.multiplied(1e-12L))<=0,"Task51 original closed capacity phase gate");}
    }
}
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
struct Runtime
{
    std::fenv_t native{};unsigned csr=0;int error=errno,threads=omp_get_max_threads(),dynamic=omp_get_dynamic();
    Runtime(){check(std::fegetenv(&native)==0,"Task51 original native FP capture");csr=_mm_getcsr();}
    ~Runtime(){omp_set_dynamic(dynamic);omp_set_num_threads(threads);std::fesetenv(&native);_mm_setcsr(csr);errno=error;}
};
struct Receipt
{
    std::fenv_t native{};unsigned csr=0;int error=0;
    void capture(){error=errno;csr=_mm_getcsr();check(std::fegetenv(&native)==0,"Task51 final native FP capture");}
    bool same(const Receipt& b)const{return error==b.error&&csr==b.csr&&std::memcmp(&native,&b.native,sizeof native)==0;}
};
struct Entry
{
    std::fenv_t native{};unsigned csr=0;
    Entry(){check(std::feraiseexcept(FE_INEXACT)==0,"Task51 prime actual public inexact");
        csr=_mm_getcsr();check(std::fegetenv(&native)==0,"Task51 common native FP capture");}
    void restore()const{check(std::fesetenv(&native)==0,"Task51 common native FP restore");_mm_setcsr(csr);errno=173;}
};
// P0-3: actual half/full producer data, independent ORIGINAL mass products.
// This is a test-only mass-block replay; it does not replace phase/WY gates.
unsigned massOracleExtents=0,massOracleNonzero=0,massOracleMicro=0,massOracleOddDensity=0;
void massOracle(const Publication& source,const Transport& actual)
{
    using Exact=LAMP::Geometry::Detail::FiniteVolume::Dyadic;
    struct Result {Exact exact;long double rounded=0,relative=0,scale=0;};
    Runtime restore;
    auto clean=restore.native;clean._Fe_stat=0;
    check(std::fesetenv(&clean)==0,"P0-3 mass oracle clean native entry");
    _mm_setcsr(restore.csr&~0x3fu);
    Entry entry; // Same nonzero sticky state and errno=173 for both replays.
    bool micro=false;
    for(std::size_t i=0;i<source.plicCells().size();++i)for(auto q:source.plicCells()[i].inventory)
        if(q>0&&q<std::ldexp(source.cells()[i].capacity,-60))micro=true;
    for(bool half:{true,false}){
        const auto& candidate=half?actual.midpoint():actual.endpoint();
        const auto& ledger=half?actual.report().boundaryHalf:actual.report().boundaryFull;
        check(ledger.enabled,"P0-3 actual open half/full ledger required");
        const auto evaluate=[&](bool grouped){
            std::array<Exact,2> before,after,in,out;
            for(std::size_t i=0;i<source.plicCells().size();++i)for(unsigned p=0;p<2;++p){
                before[p].add(source.plicCells()[i].inventory[p]);after[p].add(candidate.plicCells()[i].inventory[p]);}
            for(const auto& f:actual.faces())if((f.negativeCell>=0)!=(f.positiveCell>=0)){
                const auto& q=half?f.half:f.full;const bool forward=f.negativeCell>=0;
                for(unsigned p=0;p<2;++p){const auto positive=p==0?q.positiveLiquid:q.positiveGas;
                    const auto negative=p==0?q.negativeLiquid:q.negativeGas;
                    out[p].add(forward?positive:negative);in[p].add(forward?negative:positive);}}
            Result r;
            for(unsigned p=0;p<2;++p){
                // Retain the original per-phase RNE prefix and RN scale order.
                const auto initial=before[p].rounded();const auto target=after[p].rounded();
                const auto incoming=in[p].rounded();const auto outgoing=out[p].rounded();
                (void)target;(void)outgoing;
                auto physical=after[p];physical.add(before[p],-1);physical.add(out[p]);physical.add(in[p],-1);
                const auto density=p==0?actual.liquidDensity():actual.gasDensity();
                if(grouped)r.exact.add(physical.multiplied(density));
                else {
                    // Literal old product ordering; no grouped residual oracle.
                    for(std::size_t i=0;i<source.plicCells().size();++i){
                        r.exact.add(Exact::product({density,candidate.plicCells()[i].inventory[p]}));
                        r.exact.add(Exact::product({density,source.plicCells()[i].inventory[p]}),-1);}
                    for(const auto& f:actual.faces())if((f.negativeCell>=0)!=(f.positiveCell>=0)){
                        const auto& q=half?f.half:f.full;const auto positive=p==0?q.positiveLiquid:q.positiveGas;
                        const auto negative=p==0?q.negativeLiquid:q.negativeGas;const bool forward=f.negativeCell>=0;
                        r.exact.add(Exact::product({density,forward?positive:negative}));
                        r.exact.add(Exact::product({density,forward?negative:positive}),-1);}
                }
                r.scale+=density*(initial+incoming);
            }
            r.rounded=r.exact.rounded();r.relative=std::abs(r.rounded)/r.scale;return r;
        };
        entry.restore();const auto literal=evaluate(false);Receipt literalExit;literalExit.capture();
        entry.restore();const auto grouped=evaluate(true);Receipt groupedExit;groupedExit.capture();
        check(literalExit.same(groupedExit),"P0-3 old/grouped full native FP/rawCSR/errno");
        check(literal.exact.compare(grouped.exact)==0,"P0-3 actual exact old-products/grouped mass numerator");
        bits(literal.rounded,grouped.rounded,"P0-3 actual old/grouped final RNE mass bits/sign");
        bits(literal.scale,grouped.scale,"P0-3 original RN mass scale order");
        bits(literal.relative,grouped.relative,"P0-3 old/grouped normalized mass bits");
        bits(ledger.massResidual,literal.rounded,"P0-3 actual producer versus independent original mass products");
        bits(ledger.massRelativeConservation,literal.relative,"P0-3 actual producer original RN normalized mass");
        check(literal.scale>0&&std::isfinite(literal.scale)&&literal.relative<=1e-12L,"P0-3 original actual mass gate");
        if(literal.exact.compare(Exact{})!=0)++massOracleNonzero;
        if(micro)++massOracleMicro;
        if(actual.liquidDensity()==7.1L)++massOracleOddDensity;
        ++massOracleExtents;
    }
}

template<class Expected=std::invalid_argument,class Action> void rejectedPair(const Action& action,const char* firstLiteral=nullptr)
{
    Entry entry;std::string serialError,actualError;Receipt sf,pf;bool serialType=false,actualType=false;
    for(int threads:{1,8}){omp_set_num_threads(threads);entry.restore();
        try{action();}catch(const std::exception& error){
            (threads==1?serialError:actualError)=error.what();
            (threads==1?serialType:actualType)=typeid(error)==typeid(Expected);}
        if(threads==1)sf.capture();else pf.capture();}
    check(serialType&&actualType,"Task51 original literal exception type changed");
    check(!serialError.empty()&&actualError==serialError,"Task51 literal first exception changed");
    if(firstLiteral)check(serialError.find(firstLiteral)!=std::string::npos,"Task51 original first guard order changed");
    check(sf.same(pf),"Task51 rejected native/rawCSR/errno differs");++firstErrorPairs;
}
void scene(bool rz,unsigned kind)
{
    stage="Task51 face/gather literal pairing "+std::to_string(rz)+","+std::to_string(kind);omp_set_num_threads(1);
    auto mesh=F15BoundaryGrossContract::makeMesh(rz);const auto native=domain(*mesh,rz);
    check(!native.uniformLeaves()&&mesh->cellIndexNum()==28,"Task51 fixture must retain true asymmetric native 2:1");
    const Clock old=kind==0?Clock{7,.125,.25,Position::End}:kind==3?Clock{1,0,.125,Position::End}:Clock{};
    const Clock end{old.step+1,old.time(),old.time()+1./32,Position::End};
    const unsigned minor=rz?0:1;std::vector<double> alpha(mesh->cellIndexNum(),minor==0?0.:1.);
    if(kind==0)for(int ci=0;ci<mesh->cellIndexNum();++ci)alpha[ci]=.3+.4*mesh->cells(mesh->cellIndexToLabel[ci]).centerCoord().x;
    const auto seed=Publication::prepareLocalPlic(*mesh,native,old,alpha);
    std::vector<std::array<long double,2>> stocks;stocks.reserve(seed.plicCells().size());
    for(const auto& cell:seed.plicCells())stocks.push_back(cell.inventory);
    if(kind==3)for(std::size_t i=0;i<stocks.size();++i){stocks[i][1-minor]=seed.cells()[i].capacity;
        stocks[i][minor]=std::ldexp(seed.cells()[i].capacity,-80);}
    const auto source=kind==3?Publication::prepareLocalPlicInventory(*mesh,native,old,stocks,{}):seed;
    const auto sourceId=source.identity();const auto* sourceStorage=source.plicCells().data();const auto sourceCells=source.plicCells();
    auto rules=F15BoundaryGrossContract::rulesFor(*mesh,source,1-minor);
    if(kind==0)for(auto& rule:rules)if(rule!=Mode::Interior)rule=Mode::Closed;
    const auto boundary=Boundary::prepare(*mesh,source,rules);
    const double begin=.125,finish=kind==1?.125:(kind==3||(kind==2&&rz))?-.375:-.125;
    const auto a=kind==0?vofMaterialNormals(*mesh,source,0):F15BoundaryGrossContract::normals(*mesh,source,begin);
    const auto b=kind==0?vofMaterialNormals(*mesh,source,2):F15BoundaryGrossContract::normals(*mesh,source,finish);
    const auto slab=kind==0?Slab::prepare(*mesh,source,end,a,b):Slab::prepareBoundary(*mesh,source,boundary,end,a,b);
    const auto prepare=[&](long double rhoL,const GeometricVofTransportConfig& config){
        return kind==0?Transport::prepareLocalPlicSweep(*mesh,source,slab,rhoL,1,config):
            Transport::prepareLocalPlicBoundarySweep(*mesh,source,slab,rhoL,1,config);};
    Entry entry;omp_set_num_threads(1);entry.restore();const auto serial=prepare(100,{});Receipt sf;sf.capture();
    omp_set_num_threads(8);entry.restore();const auto actual=prepare(100,{});Receipt pf;pf.capture();
    check(sf.same(pf),"Task51 complete success native/rawCSR/errno differs from literal");
    samePayload(serial,actual);independent(*mesh,source,actual);
    check(actual.report().localPhaseUpdateCalls==8*source.plicCells().size()
        &&actual.report().localPhaseUpdateFallback==actual.report().localPhaseUpdateCalls,
        "v79 full donor-to-stock pipeline lost original half/full two-direction Core phase updates");
    if(kind!=0)massOracle(source,actual);
    if(kind==3){const auto oddDensity=prepare(7.1L,{});massOracle(source,oddDensity);}
    const auto& sr=serial.report().localFaceExecution;const auto& r=actual.report().localFaceExecution;
    check(!sr.requested&&!sr.faceValuesParallelUsed&&!sr.cellGatherParallelUsed&&!sr.fallback
        &&sr.directionsUsed==0&&sr.maximumWorkers==0,"Task51 OMP1 must execute original literal");
    check(r.requested&&r.eligible&&r.faceValuesParallelUsed&&r.cellGatherParallelUsed&&!r.fallback
        &&r.maximumWorkers>=2&&r.directionsUsed==4,"Task51 actual four directional parallel passes missing");
    check(r.nativeFaces==actual.faces().size()&&r.faceEvaluations==2*actual.faces().size()
        &&r.cellGathers==4*source.cells().size()&&r.incidences>0,"Task51 actual face/gather work counters");
    check(r.scratchBytes>0&&r.scratchBytes<=r.scratchLimitBytes&&r.scratchLimitBytes==64u*1024u*1024u,"Task51 actual scratch cap");
    if(kind==3){check(actual.report().localUnresolvedDonors>0,"Task51 actual unresolved-P0 branch missing");
        F18P0BoundaryGrossContract::inspect(*mesh,source,actual,minor,begin,finish,end);}
    if(kind==2&&rz){check(serial.report().boundaryFull.netFluxRounding[1-minor]!=0,
        "Task51 literal pure-counterflow net-RN ledger must be nonzero to catch its sign");++pureRoundingWitnesses;}
    GeometricVofTransportConfig bad;bad.relativeConservationTolerance=2e-12L;
    rejectedPair([&]{const auto unexpected=prepare(-1,bad);(void)unexpected;},"densities must be finite and positive");
    const auto foreign=Publication::prepareLocalPlicInventory(*mesh,native,old,stocks,{});
    rejectedPair([&]{actual.validateFor(*mesh,foreign,slab,100,1);});
    check(source.identity()==sourceId&&source.plicCells().data()==sourceStorage,"Task51 original source owner/storage changed");
    for(std::size_t i=0;i<sourceCells.size();++i)check(source.plicCells()[i].inventory==sourceCells[i].inventory
        &&source.plicCells()[i].measured==sourceCells[i].measured&&source.plicCells()[i].resolved==sourceCells[i].resolved,
        "Task51 original source stocks/geometry changed");
    if(kind==2){const unsigned csr=_mm_getcsr();_mm_setcsr(csr|0x8000u);
        rejectedPair<std::runtime_error>([&]{const auto unexpected=prepare(100,{});(void)unexpected;},
            "[ExactExpansionCandidate] exact arithmetic requires gradual underflow with FTZ/DAZ disabled; the floating-point mode was not changed");_mm_setcsr(csr);}
    parallelDirections+=static_cast<unsigned>(r.directionsUsed);++pairs;
    std::cout<<"LOCAL_PLIC_FACE_PARALLEL_PAIR,"<<(rz?"RZ":"XY")<<",kind,"<<kind<<",directions,"<<r.directionsUsed
        <<",workers,"<<r.maximumWorkers<<",scratch_bytes,"<<r.scratchBytes<<",original_payload_bits,exact,FP_errno,exact\n";
}
void run()
{
    Runtime restore;omp_set_dynamic(0);
    check(sizeof(long double)==sizeof(double)&&std::numeric_limits<long double>::digits==53,"Task51 MSVC binary64 native contract required");
    std::cout<<"LOCAL_PLIC_FACE_PARALLEL_ENV,round,"<<std::fegetround()<<",mxcsr,"<<_mm_getcsr()<<"\n";
    for(bool rz:{false,true})for(unsigned kind=0;kind<4;++kind)scene(rz,kind);
    check(pairs==8&&parallelDirections==32&&firstErrorPairs==18&&zeroNetGross>0&&physicalLedgers==24&&pureRoundingWitnesses==1,"Task51 minimum finite coverage incomplete");
    check(massOracleExtents==16&&massOracleNonzero>0&&massOracleMicro==8&&massOracleOddDensity==4,
        "P0-3 actual half/full, nonzero mass, microphase and non-power-of-two density witnesses");
#if !defined(__FAST_MATH__) && (!defined(_M_FP_FAST) || !_M_FP_FAST)
    static_assert(LAMP::Geometry::Detail::FiniteVolume::NativeSmallExact::supported,
        "P0-3 production exact-scaling compile format must be qualified");
    std::cout<<"LOCAL_PLIC_OPEN_MASS_CONTRACTION_PASS,actual_extents,"<<massOracleExtents<<",nonzero_exact_mass,"<<massOracleNonzero
        <<",microphase_extents,"<<massOracleMicro<<",non_power_of_two_density_extents,"<<massOracleOddDensity
        <<",old_product_oracle,exact,full_native_rawCSR_errno,exact,compile_format,native_x64_binary64,no_speed_or_resource_fault_claim\n";
#endif
    std::cout<<"LOCAL_PLIC_FACE_PARALLEL_PASS,paired_factories,"<<pairs<<",actual_parallel_directions,"<<parallelDirections
        <<",first_error_pairs,"<<firstErrorPairs<<",zero_net_gross_faces,"<<zeroNetGross<<",phase_ledgers,"<<physicalLedgers<<",nonzero_literal_pure_rounding,"<<pureRoundingWitnesses
        <<",threads,1_vs_8,full_donor_to_stock,exact,all_original_bits,exact,full_native_rawCSR_errno,exact,no_speed_or_resource_fault_claim\n";
}
#else
void run(){throw std::runtime_error("Task51 parallel contract requires actual MSVC x64 OpenMP");}
#endif
}
// END TASK51 LITERAL/PARALLEL FACE VALUES AND ORDERED CELL GATHER


namespace {
namespace ConservativeVofMomentumContract {
void physicalFamilies()
{
    for(bool rz:{false,true})for(bool amr:{false,true}){
        std::unique_ptr<QuadtreeMesh> mesh;if(amr)mesh=vofMaterialMesh(4);else{DictionaryMesh q;q.Nx=q.Ny=4;q.maxCellNum=4096;q.maxInitialIterations=0;mesh=std::make_unique<QuadtreeMesh>(q);}
        const auto source=Publication::prepareLocalPlic(*mesh,domain(*mesh,rz),Clock{},l1PlaneAlpha(*mesh,rz,0));
        const auto a=vofMaterialNormals(*mesh,source,0);const auto dt=static_cast<double>(.1L/vofMaterialOutflow(*mesh,source,a,a));
        const Clock end{1,0,dt,Position::End};const auto slab=Slab::prepare(*mesh,source,end,a,a);
        const auto transport=Transport::prepareLocalPlicSweep(*mesh,source,slab,100,1);
        const auto geometry=HydroMomentumGeometry::build(*mesh,vofMaterialMetrics(*mesh,source));
        const auto before=HydroPhaseGeometry::build(*mesh,geometry.metrics(),source),middle=HydroPhaseGeometry::build(*mesh,geometry.metrics(),transport.midpoint()),after=HydroPhaseGeometry::build(*mesh,geometry.metrics(),transport.endpoint());
        const auto material=HydroMaterialTransport::prepareGeometricVof(*mesh,geometry,before,middle,after,transport);
        l1InventoryProjectionPair(*mesh,source,transport,geometry,material,before,middle,after);
        if(amr)l1CoonsQuadrants(*mesh,source,transport);
        // A different actual source provides a translating density jump. RZ
        // uses only axial velocity so the true axis is exactly impermeable.
        std::vector<double> alpha(mesh->cellIndexNum());for(int ci=0;ci<mesh->cellIndexNum();++ci){const auto box=bounds(*mesh,ci);alpha[ci]=static_cast<double>(std::clamp((.43L-box.y0)/(box.y1-box.y0),0.L,1.L));}
        const auto front=Publication::prepareLocalPlic(*mesh,domain(*mesh,rz),Clock{},alpha);using Mode=Slab::Boundary::Mode;
        std::vector<Mode> rules(front.domain().nativeFaceCount(),Mode::Interior);std::vector<double> normal(rules.size());
        for(int row=0;row<front.domain().nativeFaceCount();++row){const auto f=front.domain().face(*mesh,row);normal[row]=f.horizontal?.125:0;
            if(f.negativeCell<0||f.positiveCell<0)rules[row]=f.horizontal?(f.outwardSign<0?Mode::PureLiquidReservoir:Mode::PureGasReservoir):Mode::Closed;}
        const auto law=Slab::Boundary::prepare(*mesh,front,rules);const Clock target{1,0,.03125,Position::End};
        const auto movingSlab=Slab::prepareBoundary(*mesh,front,law,target,normal,normal);
        const auto moving=Transport::prepareLocalPlicBoundarySweep(*mesh,front,movingSlab,1000,1.226);
        const auto geom=HydroMomentumGeometry::build(*mesh,vofMaterialMetrics(*mesh,front));const auto& layout=geom.layout();
        const auto p0=HydroPhaseGeometry::build(*mesh,geom.metrics(),front),pm=HydroPhaseGeometry::build(*mesh,geom.metrics(),moving.midpoint()),p1=HydroPhaseGeometry::build(*mesh,geom.metrics(),moving.endpoint());
        const auto mt=HydroMaterialTransport::prepareGeometricVof(*mesh,geom,p0,pm,p1,moving);
        const auto m0=HydroMomentumMass::buildCellAverageProjection(*mesh,geom,p0,1000,1.226),mm=HydroMomentumMass::buildCellAverageProjection(*mesh,geom,pm,1000,1.226),m1=HydroMomentumMass::buildCellAverageProjection(*mesh,geom,p1,1000,1.226);
        const auto dual=HydroDualMaterial::prepare(*mesh,mt,m0,mm,m1);qualify(*mesh,dual);
        std::vector<int> walls;std::vector<double> u(layout.velocityCount());for(int i=0;i<layout.velocityCount();++i){const auto f=front.domain().face(*mesh,layout.nativeFace(i));u[i]=f.horizontal?.125:0;if(rules[layout.nativeFace(i)]==Mode::Closed)walls.push_back(i);}
        const HydroNativeValues values(layout,u,std::vector<double>(layout.pressureCount()));const auto wallLaw=HydroAdvectorBoundary::impermeable(*mesh,geom,walls);
        std::vector<Vector2D> first(geom.segments().size()),last(first);for(std::size_t k=0;k<first.size();++k){const auto& seg=geom.segments()[k];if(seg.negativeVelocity>=0&&seg.positiveVelocity>=0)continue;
            const int owner=seg.negativeVelocity>=0?seg.negativeVelocity:seg.positiveVelocity;const auto outward=(seg.negativeVelocity>=0?1.L:-1.L)*dual.faceRates()[k];if(outward<0)first[k]=last[k]={0,.125};(void)owner;}
        const auto inlet=HydroConvectionBoundary::bindVofEndpointVelocity(*mesh,dual,first,last,dual.sources(),dual.timeLayer());
        const auto op=HydroKepConvection::buildVofPrescribedFlow(*mesh,dual,wallLaw,inlet);
        std::vector<double> traces(geom.segments().size());for(std::size_t k=0;k<traces.size();++k){const auto& seg=geom.segments()[k];if(seg.negativeVelocity>=0&&seg.positiveVelocity>=0)continue;
            const int owner=seg.negativeVelocity>=0?seg.negativeVelocity:seg.positiveVelocity;const auto out=(seg.negativeVelocity>=0?1.L:-1.L)*dual.faceRates()[k];
            traces[k]=out<0?inlet.traces()[k]:out>0?u[owner]:0;}
        const auto actualBoundary=HydroConvectionBoundary::bind(*mesh,dual,traces,dual.sources(),dual.timeLayer());
        const auto result=op.apply(*mesh,values,actualBoundary,dual.sources(),dual.timeLayer());
        for(int i=0;i<layout.velocityCount();++i){check(result.momentumLoads()[i]==-result.momentumRates()[i],"constant density-jump applied load retained uR");
            X::Dyadic residual=X::Dyadic::product({dual.endMass().masses()[i],u[i]});residual.add(X::Dyadic::product({dual.beginMass().masses()[i],u[i]}),-1);
            residual.add(X::Dyadic::product({dual.deltaTime(),result.momentumRates()[i]}));
            X::Dyadic scale(std::abs(dual.beginMass().masses()[i]*u[i]));scale.add(std::abs(dual.endMass().masses()[i]*u[i]));
            scale.add(X::Dyadic::product({dual.deltaTime(),std::abs(result.momentumRates()[i])}));
            check(admissible(residual,scale),"constant velocity through actual moving density jump is not retained to roundoff");}
        ++constantCases;
    }
}
void run()
{
    physicalFamilies();localPlicAtomsOnly(true,true);localPlicOpenKepOnly();F16OpenMomentumContract::run();F17OpenMomentumStepContract::run();localPlicOpenStepOnly();
    check(qualified>0&&rows>0&&coarseFine>0&&axes>0&&finiteRows>0&&negativeControls==1&&constantCases==4,"conservative VOF finite required coverage incomplete");
    std::cout<<"LOCAL_PLIC_CONSERVATIVE_MOMENTUM_PASS,qualified,"<<qualified<<",rows,"<<rows<<",coarse_fine,"<<coarseFine<<",axis_faces,"<<axes<<",finite_overlap_rows,"<<finiteRows
        <<",constant_jump_cases,"<<constantCases<<",negative_controls,"<<negativeControls<<",max_relative_canonical_R,"<<peakCanonical<<",max_relative_stored_M_R,"<<peakStored
        <<",max_relative_commutation,"<<peakCommutation<<",endpoint_mass_and_gross_scale,256doubleeps,actual_load_minus_N,negative_half_R_energy,no_cell_or_flux_repair\n";
}
}
} // same unnamed namespace as the independent qualification oracle

namespace SelectedFluidTransportContract
{
using Mode=Slab::Boundary::Mode;
unsigned steps=0,inactiveChecks=0,wallChecks=0,changedCells=0;
long double peakLedger=0;
std::unique_ptr<QuadtreeMesh> nozzleMesh(bool amr)
{
    DictionaryMesh d;d.Nx=10;d.Ny=32;d.minCoord={0,0};d.maxCoord={5,16};
    d.maxCellNum=8192;d.maxInitialIterations=0;d.maxLevel_local={2};
    auto mesh=std::make_unique<QuadtreeMesh>(d);
    if(amr){DictionaryMesh::RefineBox b;b.level=1;b.min={0,3};b.max={1,5};mesh->applyRefineBoxes({b});}
    return mesh;
}
NativeLeafDomain selected(const QuadtreeMesh& mesh,bool rz)
{
    FaceSegmentGraph graph;graph.build(mesh,[](const QuadtreeMesh&,const Cell& c){
        const auto p=c.centerCoord();return !(p.x>1&&p.y<3);});
    return NativeLeafDomain::build(mesh,graph,rz);
}
std::array<long double,2> inventories(const Publication& p)
{
    std::array<long double,2> total{};for(int i:p.domain().cellIndexByRow())
        for(unsigned k=0;k<2;++k)total[k]+=p.plicCells()[i].inventory[k];
    return total;
}
void inspect(const QuadtreeMesh& mesh,const Publication& before,const Transport& moved)
{
    for(const auto* after:{&moved.midpoint(),&moved.endpoint()}){
        check(after->domain().cellRowByIndex()==before.domain().cellRowByIndex(),"selected transport changed active map");
        for(std::size_t i=0;i<after->cells().size();++i){
            if(after->domain().cellRowByIndex()[i]<0){
                check(after->cells()[i].nativeCell<0&&after->cells()[i].capacity==0&&after->alpha()[i]==0
                    &&after->plicCells()[i].inventory==std::array<long double,2>{0,0},"solid gained stock/capacity/phase identity");
                ++inactiveChecks;
            }else if(after->plicCells()[i].inventory!=before.plicCells()[i].inventory)++changedCells;
        }
    }
    for(const auto f:before.domain().faces(mesh)){
        if(f.negativeCell>=0&&f.positiveCell>=0)continue;
        const auto& trace=moved.faces()[f.nativeRow];
        const bool excluded=f.negativeCell<0&&f.positiveCell<0;
        const bool wall=(f.horizontal&&f.center.y==3&&f.center.x>1)
            ||(!f.horizontal&&f.center.x==1&&f.center.y<3);
        if(excluded||wall){
            for(const auto* flux:{&trace.half,&trace.full})check(flux->total==0&&flux->liquid==0&&flux->gas==0
                &&flux->mass==0&&flux->positiveTotal==0&&flux->negativeTotal==0
                &&flux->positiveLiquid==0&&flux->negativeLiquid==0&&flux->positiveGas==0&&flux->negativeGas==0,
                "solid or electrode wall acquired shared flux");
            if(wall)++wallChecks;
        }
    }
    check(moved.report().halfRelativeConservation<=conservationTolerance&&moved.report().fullRelativeConservation<=conservationTolerance,
        "selected shared transport conservation gate");
    ++steps;
}
void runScene(bool rz,bool amr)
{
    stage=std::string("selected nozzle transport ")+(rz?"RZ":"XY")+(amr?" 2:1":" uniform");
    auto mesh=nozzleMesh(amr);const auto fluid=selected(*mesh,rz);
    auto source=GeometricVofInitialization::prepare(*mesh,fluid,[](Vector2D p){
        return std::max(p.x-1,p.y-4.125);});
    const auto n=static_cast<std::size_t>(fluid.nativeFaceCount());
    std::vector<Mode> rules(n,Mode::Closed);std::vector<double> through(n),vortex(n);
    int excludedRow=-1,wallRow=-1;unsigned hanging=0;
    const auto psi=[](Point p){
        if(p.x>=1||p.y<=3||p.y>=5)return 0.L;
        return p.x*p.x*(1-p.x)*(1-p.x)*(p.y-3)*(p.y-3)*(5-p.y)*(5-p.y)/256;
    };
    for(const auto f:fluid.faces(*mesh)){
        if(f.negativeCell<0&&f.positiveCell<0){excludedRow=f.nativeRow;continue;}
        const auto e=faceEndpoints(*mesh,f);const bool interior=f.negativeCell>=0&&f.positiveCell>=0;
        if(interior)rules[f.nativeRow]=Mode::Interior;
        else if(f.horizontal&&f.center.y==0)rules[f.nativeRow]=Mode::PureLiquidReservoir;
        else if(f.horizontal&&f.center.y==16)rules[f.nativeRow]=Mode::PureGasReservoir;
        else if((f.horizontal&&f.center.y==3&&f.center.x>1)||(!f.horizontal&&f.center.x==1&&f.center.y<3))wallRow=f.nativeRow;
        if(f.horizontal&&e[0].x<1&&e[1].x<=1)through[f.nativeRow]=1;
        const auto q=(f.horizontal?-1.L:1.L)*(psi(e[1])-psi(e[0]));
        if(f.area>0)vortex[f.nativeRow]=static_cast<double>(q/f.area);
        else check(q==0,"selected streamfunction leaks at true axis");
        if(!interior)check(vortex[f.nativeRow]==0,"selected streamfunction leaks at wall/open domain");
        if(interior&&mesh->cells(mesh->cellIndexToLabel[f.negativeCell]).level()!=
            mesh->cells(mesh->cellIndexToLabel[f.positiveCell]).level())++hanging;
    }
    check(excludedRow>=0&&wallRow>=0,"selected fixture lacks excluded row or actual electrode wall");
    if(amr)check(hanging>0,"selected AMR fixture lacks actual coarse/fine subfaces");
    const auto initial=source;
    const std::vector<double> zero(n);
    auto boundary=Slab::Boundary::prepare(*mesh,source,rules);
    auto badRules=rules;badRules[excludedRow]=Mode::PureGasReservoir;
    rejects([&]{Slab::Boundary::prepare(*mesh,source,badRules);},"excluded row cannot be gas reservoir");
    badRules[excludedRow]=Mode::Interior;
    rejects([&]{Slab::Boundary::prepare(*mesh,source,badRules);},"excluded row cannot be fluid interior");
    auto badVelocity=zero;badVelocity[excludedRow]=1;
    const Clock firstEnd{1,0,.01,Position::End};
    rejects([&]{Slab::prepareBoundary(*mesh,source,boundary,firstEnd,badVelocity,zero);},"excluded slab normal");
    rejects([&]{Transport::localPlicBoundaryCourantRates(*mesh,source,boundary,badVelocity);},"excluded rate normal");
    badVelocity=zero;badVelocity[wallRow]=1;
    rejects([&]{Slab::prepareBoundary(*mesh,source,boundary,firstEnd,badVelocity,zero);},"electrode wall normal");
    // Closed circulation changes the real interface without boundary stock.
    for(unsigned k=0;k<2;++k){
        const auto t=source.timeLayer();const Clock end{t.step+1,t.time(),t.time()+.01,Position::End};
        const auto slab=Slab::prepare(*mesh,source,end,vortex,vortex);
        const auto moved=Transport::prepareLocalPlicSweep(*mesh,source,slab,1240,1.184L);
        inspect(*mesh,source,moved);source=moved.endpoint();
    }
    source=initial;
    const long double section=rz?CartesianMeasure::twoPi/2:1.L;
    for(unsigned k=0;k<4;++k){
        const bool reverse=k==3;auto last=through;if(reverse)for(auto& v:last)v*=-.5;
        const auto t=source.timeLayer();const Clock end{t.step+1,t.time(),t.time()+.01,Position::End};
        const auto b=Slab::Boundary::prepare(*mesh,source,rules);
        const auto rates=Transport::localPlicBoundaryCourantRates(*mesh,source,b,through);
        check(rates.maximumIncoming>0&&rates.maximumOutflow>0&&rates.maximumLocalColumn>0,"selected rate lost real flux");
        const auto slab=Slab::prepareBoundary(*mesh,source,b,end,through,last);
        const auto moved=Transport::prepareLocalPlicBoundarySweep(*mesh,source,slab,1240,1.184L);
        inspect(*mesh,source,moved);
        for(bool half:{true,false}){
            const auto& target=half?moved.midpoint():moved.endpoint();const auto a=inventories(source),c=inventories(target);
            const auto duration=static_cast<long double>(end.end)-end.begin;
            const auto fraction=reverse?(half?5.L/16:1.L/4):(half?.5L:1.L);
            const auto net=section*duration*fraction;
            near(c[0]-a[0],net,rz?328*CartesianMeasure::twoPi/2:68.L,"independent selected liquid supply balance",conservationTolerance);
            near(c[1]-a[1],-net,rz?328*CartesianMeasure::twoPi/2:68.L,"independent selected gas displacement balance",conservationTolerance);
            const auto& ledger=half?moved.report().boundaryHalf:moved.report().boundaryFull;
            const auto positive=section*duration*(reverse?(half?5.L/16:1.L/3):(half?.5L:1.L));
            const auto negative=reverse&&!half?section*duration/12:0.L;
            near(ledger.inflow[0],positive,section*duration,"independent gross liquid inflow",conservationTolerance);
            near(ledger.outflow[1],positive,section*duration,"independent gross gas outflow",conservationTolerance);
            near(ledger.outflow[0],negative,section*duration,"independent reverse liquid outflow",conservationTolerance);
            near(ledger.inflow[1],negative,section*duration,"independent reverse gas inflow",conservationTolerance);
            peakLedger=std::max(peakLedger,ledger.physicalRelativeConservation);
        }
        const auto endpoints=faceEndpoints(*mesh,fluid.face(*mesh,excludedRow));
        rejects([&]{moved.restrictNativeFaces(*mesh,{{excludedRow,endpoints[0],endpoints[1]}});},"excluded accepted trace query");
        source=moved.endpoint();
    }
    rejects([&]{boundary.validateFor(*mesh,source);},"old selected boundary cannot bind new accepted source");
}
void run()
{
    for(bool rz:{false,true})for(bool amr:{false,true})runScene(rz,amr);
    check(steps==24&&inactiveChecks>0&&wallChecks>0&&changedCells>0,"selected fluid transport coverage missing");
    std::cout<<std::setprecision(17)<<"SELECTED_FLUID_TRANSPORT_PASS steps="<<steps<<" inactive_checks="<<inactiveChecks
        <<" wall_checks="<<wallChecks<<" changed_cells="<<changedCells<<" peak_boundary_conservation="<<peakLedger<<'\n';
}
}
namespace PinnedCornerTransportContract
{
using Pin=Publication::PinCertificate;
namespace Cut=LAMP::Geometry::PlanarCellCut;
unsigned steps=0,changingSteps=0,layerChecks=0,wallChecks=0,solidChecks=0,refines=0,coarsens=0,coincidentLayers=0;
long double peakPhaseResidual=0,peakGeometryResidual=0,peakDivergence=0,peakPinnedCfl=0;
unsigned pinnedSweeps=0,positivePinFaces=0,negativePinFaces=0,pinPeakRejections=0;
unsigned wholePhaseExhaustions=0,wholePhaseControls=0;long double peakExhaustionCorrection=0;
struct Support { int vertex=-1,contact=-1,wetCell=-1,dryCell=-1,wetRow=-1,dryRow=-1; };
Support support(const QuadtreeMesh& mesh,const Publication& p)
{
    Support s;unsigned vertices=0,active=0,solid=0,wetRows=0,dryRows=0;
    for(int v=0;v<mesh.vertIndexNum();++v){const auto q=mesh.verts(mesh.vertIndexToLabel[v]).coord();
        if(q.x==1&&q.y==3){s.vertex=v;++vertices;}}
    for(int ci=0;ci<mesh.cellIndexNum();++ci){const auto b=bounds(mesh,ci);
        if(!((b.x0==1||b.x1==1)&&(b.y0==3||b.y1==3)))continue;
        if(p.domain().cellRowByIndex()[ci]<0){++solid;continue;}++active;
        if(b.x1==1&&b.y0==3)s.contact=ci;
        if(b.x1==1&&b.y1==3)s.wetCell=ci;
        if(b.x0==1&&b.y0==3)s.dryCell=ci;
    }
    check(vertices==1&&active==3&&solid==1&&s.contact>=0&&s.wetCell>=0&&s.dryCell>=0,
        "pin fixture lacks one actual three-fluid/one-solid native corner");
    for(const auto f:p.domain().faces(mesh)){const auto e=faceEndpoints(mesh,f);
        if(!f.horizontal&&e[0].x==1&&e[1].x==1&&e[0].y<3&&e[1].y==3){
            s.wetRow=f.nativeRow;++wetRows;
            check((f.negativeCell==s.wetCell&&f.positiveCell<0)||(f.positiveCell==s.wetCell&&f.negativeCell<0),
                "wet leg is not the actual selected solid wall");}
        if(f.horizontal&&e[0].y==3&&e[1].y==3&&e[0].x==1&&e[1].x>1){
            s.dryRow=f.nativeRow;++dryRows;
            check((f.negativeCell==s.dryCell&&f.positiveCell<0)||(f.positiveCell==s.dryCell&&f.negativeCell<0),
                "dry leg is not the actual selected solid wall");}
    }
    check(wetRows==1&&dryRows==1,"pin fixture lacks unique actual wall legs");return s;
}
void inspect(const QuadtreeMesh& mesh,const Publication& p,Clock expected)
{
    p.validateFor(mesh,p.identity(),expected);const auto s=support(mesh,p);const auto pin=p.pinnedCorner();
    check(pin.has_value(),"typed descendant dropped the corner pin");pin->validateFor(mesh,p);
    check(pin->sourceIdentity()==p.identity()&&sameClock(pin->sourceTime(),expected)&&sameClock(p.timeLayer(),expected),
        "pin certificate is not bound to actual owner and complete stage clock");
    const auto point=pin->physicalPoint();
    check(point.x==1&&point.y==3&&pin->nativeVertex()==s.vertex&&pin->nativeCell()==s.contact
        &&pin->corner()==Cut::RectangleCorner::LowerRight&&pin->wetNativeRow()==s.wetRow&&pin->dryNativeRow()==s.dryRow,
        "pin certificate lost actual physical point/native support");
    const auto& wet=p.plicCells()[s.wetCell];const auto& dry=p.plicCells()[s.dryCell];
    check(wet.inventory[0]>0&&wet.inventory[1]==0&&dry.inventory[0]==0&&dry.inventory[1]>0,
        "closed pin circulation contaminated wet-liquid or dry-gas wall support");
    const auto& c=p.plicCells()[s.contact];const bool mixed=c.inventory[0]>0&&c.inventory[1]>0;
    check(c.resolved&&pin->kind()==(mixed?Pin::Kind::ResolvedCorner:Pin::Kind::NativeCoincident),
        "pin kind does not describe actual retained contact stocks");
    const auto segment=Cut::interfaceSegment(c.rectangle,c.normal,c.offset);
    check(segment.has_value(),"actual pinned PLIC has no representable segment");
    const auto cx=(static_cast<long double>(c.rectangle.lower.x)+c.rectangle.upper.x)/2;
    const auto cy=(static_cast<long double>(c.rectangle.lower.y)+c.rectangle.upper.y)/2;
    const Point corner{1-cx,3-cy};const auto equal=[](Point a,Point b){return a.x==b.x&&a.y==b.y;};
    check(equal((*segment)[0],corner)||equal((*segment)[1],corner),
        "actual PLIC segment endpoint moved off exact native physical corner");
    for(unsigned phase=0;phase<2;++phase){
        if(c.inventory[phase]==0)check(c.measured[phase]==0,"pure contact invented a measured microphase");
        else{const auto residual=std::abs(c.measured[phase]-c.inventory[phase])/c.inventory[phase];
            peakGeometryResidual=std::max(peakGeometryResidual,residual);
            check(residual<=1e-12L,"actual pinned phase geometry lost independent stock accuracy");}}
    if(!mixed){
        unsigned matched=0;const int opposite=c.inventory[1]==0?s.dryCell:s.wetCell;
        for(const auto f:p.domain().faces(mesh)){
            if(!((f.negativeCell==s.contact&&f.positiveCell==opposite)||(f.positiveCell==s.contact&&f.negativeCell==opposite)))continue;
            const auto e=faceEndpoints(mesh,f);const Point a{e[0].x-cx,e[0].y-cy},b{e[1].x-cx,e[1].y-cy};
            if((equal((*segment)[0],a)&&equal((*segment)[1],b))||(equal((*segment)[0],b)&&equal((*segment)[1],a)))++matched;
        }
        check(matched==1,"pure contact is not one actual complete opposite-phase native face");++coincidentLayers;
    }
    for(int ci=0;ci<mesh.cellIndexNum();++ci)if(p.domain().cellRowByIndex()[ci]<0){
        check(p.cells()[ci].nativeCell<0&&p.cells()[ci].capacity==0&&p.alpha()[ci]==0
            &&p.plicCells()[ci].inventory==std::array<long double,2>{0,0},"pinned solid acquired fluid stock or capacity");++solidChecks;}
    ++layerChecks;
}
std::vector<double> circulation(const QuadtreeMesh& mesh,const Publication& p)
{
    // q=delta(psi) is a discrete solenoidal physical-volume flux in XY/RZ.
    // Both solid wall legs have psi=0, but the two INTERNAL continuation
    // rays have nonzero native means. Point velocity vanishes only at pin.
    const auto psi=[](Point q){if(q.x<=0||q.x>=2||q.y<=0||q.y>=5)return 0.L;
        const auto left=std::min(0.L,q.x-1),up=std::max(0.L,q.y-3);
        return (left*left+up*up)*q.x*q.x*(2-q.x)*(2-q.x)*q.y*q.y*(5-q.y)*(5-q.y)/256;};
    std::vector<double> velocity(p.domain().nativeFaceCount());
    std::vector<long double> divergence(mesh.cellIndexNum()),gross(divergence.size());unsigned movingFaces=0,hanging=0;
    for(const auto f:p.domain().faces(mesh)){const auto e=faceEndpoints(mesh,f);
        const auto exact=(f.horizontal?-1.L:1.L)*(psi(e[1])-psi(e[0]));
        if(f.area>0)velocity[f.nativeRow]=static_cast<double>(exact/f.area);
        else check(exact==0,"pin streamfunction leaks at the actual RZ axis");
        if(f.negativeCell<0||f.positiveCell<0)check(velocity[f.nativeRow]==0,"pin streamfunction leaks through selected wall/excluded face");
        const auto q=static_cast<long double>(f.area)*velocity[f.nativeRow];movingFaces+=q!=0;
        if(f.negativeCell>=0){divergence[f.negativeCell]+=q;gross[f.negativeCell]+=std::abs(q);}
        if(f.positiveCell>=0){divergence[f.positiveCell]-=q;gross[f.positiveCell]+=std::abs(q);}
        if(f.negativeCell>=0&&f.positiveCell>=0&&mesh.cells(mesh.cellIndexToLabel[f.negativeCell]).level()!=
            mesh.cells(mesh.cellIndexToLabel[f.positiveCell]).level())++hanging;
    }
    check(movingFaces>0,"pin circulation is not an actual moving native field");
    if(!p.domain().uniformLeaves())check(hanging>0,"pinned refined fixture lacks actual coarse/fine native faces");
    for(int ci:p.domain().cellIndexByRow()){
        if(gross[ci]==0)check(divergence[ci]==0,"zero-flow cell has nonzero discrete divergence");
        else{const auto relative=std::abs(divergence[ci])/gross[ci];peakDivergence=std::max(peakDivergence,relative);
            check(relative<=64.L*std::numeric_limits<double>::epsilon(),"pin psi field lost discrete volume conservation");}}
    return velocity;
}
void phaseConservation(const Publication& before,const Publication& after)
{
    const auto a=SelectedFluidTransportContract::inventories(before),b=SelectedFluidTransportContract::inventories(after);
    for(unsigned phase=0;phase<2;++phase){const auto residual=std::abs(b[phase]-a[phase])/a[phase];
        peakPhaseResidual=std::max(peakPhaseResidual,residual);
        near(b[phase],a[phase],a[phase],"pinned independent closed phase conservation",1e-12L);}
}
void wallFluxes(const QuadtreeMesh& mesh,const Publication& before,const Transport& moved)
{
    for(const auto f:before.domain().faces(mesh))if(f.negativeCell<0||f.positiveCell<0){
        for(const auto* flux:{&moved.faces()[f.nativeRow].half,&moved.faces()[f.nativeRow].full})
            check(flux->total==0&&flux->liquid==0&&flux->gas==0&&flux->mass==0&&flux->positiveTotal==0&&flux->negativeTotal==0
                &&flux->positiveLiquid==0&&flux->negativeLiquid==0&&flux->positiveGas==0&&flux->negativeGas==0,
                "pinned actual wall/excluded native row acquired a gross phase flux");++wallChecks;}
}
void pinFluxes(const QuadtreeMesh& mesh,const Publication& before,const Transport& moved,const std::vector<double>& velocity)
{
    const auto pin=before.pinnedCorner();unsigned faces=0;
    const auto physical=pin->physicalPoint();const Point point{physical.x,physical.y};
    for(const auto f:before.domain().faces(mesh)){
        if(f.negativeCell<0||f.positiveCell<0||f.area<=0)continue;
        const auto e=faceEndpoints(mesh,f);const auto same=[](Point a,Point b){return a.x==b.x&&a.y==b.y;};
        if(!same(e[0],point)&&!same(e[1],point))continue;
        check(velocity[f.nativeRow]!=0,"pinned sweep test has zero mean on an internal pin continuation");++faces;
        for(const auto* flux:{&moved.faces()[f.nativeRow].half,&moved.faces()[f.nativeRow].full}){
            check(flux->total!=0,"nonzero internal pin face was not consumed");
            if(flux->total>0)++positivePinFaces;else ++negativePinFaces;
            // A strict zero, not a tolerance, detects the former 1.89e-23 wedge.
            if(f.horizontal)check(flux->gas==0&&flux->positiveGas==0&&flux->negativeGas==0&&flux->liquid==flux->total,
                "pinned horizontal sweep transported a gas micro-wedge to wet support");
            else check(flux->liquid==0&&flux->positiveLiquid==0&&flux->negativeLiquid==0&&flux->gas==flux->total,
                "pinned vertical sweep transported a liquid micro-wedge to dry support");
            check(flux->positiveLiquid+flux->positiveGas==flux->positiveTotal&&
                flux->negativeLiquid+flux->negativeGas==flux->negativeTotal,"pin polygon changed canonical gross total FT");
        }
    }
    check(faces==2&&moved.report().localPinnedSweeps==4,"both half/full actual pin-incident polygons were not consumed");
    check(moved.report().maximumPinnedPeakCfl>0&&moved.report().maximumPinnedPeakCfl<=.5L,
        "actual pin peak-depth Courant is missing or unsafe");
    const auto rates=Transport::localPlicCourantRates(mesh,before,velocity);
    check(rates.maximumLocalColumn>0,"sealed-pin instantaneous rate omitted actual positive flow");
    peakPinnedCfl=std::max(peakPinnedCfl,moved.report().maximumPinnedPeakCfl);
    pinnedSweeps+=static_cast<unsigned>(moved.report().localPinnedSweeps);
}
void peakCflReject(const QuadtreeMesh& mesh,const Publication& source)
{
    const auto s=support(mesh,source);const auto b=bounds(mesh,s.contact);
    const Point site{1,b.y1};
    // One actual native-vertex circulation has old rectangular/whole-cell
    // Courant .375, while the pin triangle's peak depth exceeds .5 widths.
    const long double amplitude=.375L*source.cells()[s.contact].capacity;
    const auto psi=[&](Point q){return q.x==site.x&&q.y==site.y?amplitude:0.L;};
    std::vector<double> velocity(source.domain().nativeFaceCount());
    long double originalMaximum=0;
    for(const auto f:source.domain().faces(mesh)){
        const auto e=faceEndpoints(mesh,f);const auto q=(f.horizontal?-1.L:1.L)*(psi(e[1])-psi(e[0]));
        if(f.area>0)velocity[f.nativeRow]=static_cast<double>(q/f.area);else check(q==0,"CFL negative control leaks on axis");
        if(q==0)continue;
        check(f.negativeCell>=0&&f.positiveCell>=0,"CFL negative control leaks through selected wall");
        const auto donor=q>0?f.negativeCell:f.positiveCell;const auto d=bounds(mesh,donor);
        const auto tangent=f.horizontal?std::abs(e[1].x-e[0].x):std::abs(e[1].y-e[0].y);
        const auto wholeTangent=f.horizontal?d.x1-d.x0:d.y1-d.y0;
        const auto column=source.cells()[donor].capacity*tangent/wholeTangent;
        originalMaximum=std::max(originalMaximum,std::abs(q)/column);
    }
    check(originalMaximum<.5L,"negative control already violates old mean-column Courant");
    const auto rates=Transport::localPlicCourantRates(mesh,source,velocity);
    check(rates.maximumLocalColumn>.5L,"pin-aware instantaneous rate missed an unsafe profile peak");
    const Clock end{source.timeLayer().step+1,source.timeLayer().time(),source.timeLayer().time()+1,Position::End};
    const auto slab=Slab::prepare(mesh,source,end,velocity,velocity);
    bool rejected=false;
    try{(void)Transport::prepareLocalPlicSweep(mesh,source,slab,1240,1.184L);}
    catch(const GeometricVofCourantExceeded& error){rejected=std::string(error.what()).find("pinned actual peak")!=std::string::npos;}
    check(rejected,"prepare did not reject actual pin peak depth with old mean Courant still legal");++pinPeakRejections;
    inspect(mesh,source,source.timeLayer());
}
void wholePhaseReceipts(const Publication& before,const Transport& moved)
{
    unsigned corrections=0;
    for(const auto& face:moved.faces())for(const auto* flux:{&face.half,&face.full}){
        if(flux->exhaustionLiquidCorrection==0&&flux->exhaustionGasCorrection==0)continue;
        ++corrections;check(face.negativeCell>=0&&face.positiveCell>=0&&flux->total!=0,
            "whole-phase receipt is not on an actual nonzero INTERNAL shared face");
        check(flux->limiter==1,"whole-phase receipt changed the meaning of the original limiter");
        check(flux->liquid==flux->rawLiquid+flux->exhaustionLiquidCorrection,
            "whole-phase accepted liquid lost its separately retained raw geometric flux and correction");
        check(flux->positiveLiquid+flux->positiveGas==flux->positiveTotal&&
            flux->negativeLiquid+flux->negativeGas==flux->negativeTotal,"whole-phase correction changed canonical FT");
    }
    check(corrections==moved.report().localWholePhaseExhaustions,"actual native shared corrections do not match strict geometry receipts");
    if(corrections){
        check(before.pinnedCorner().has_value()&&moved.report().maximumWholePhaseExhaustionRelativeCorrection>0&&
            moved.report().maximumWholePhaseExhaustionRelativeCorrection<=1e-12L,"whole-phase correction exceeded original own-phase budget");
        wholePhaseExhaustions+=corrections;
        peakExhaustionCorrection=std::max(peakExhaustionCorrection,moved.report().maximumWholePhaseExhaustionRelativeCorrection);
    }
}
void wholePhaseStepControls(const QuadtreeMesh& mesh,const Publication& before,const Transport& moved)
{
    const bool fullControl=(before.domain().axisymmetric()&&moved.endTime().step==6)||
        (!before.domain().axisymmetric()&&moved.endTime().step==8);
    if(!fullControl)return;
    const auto star=support(mesh,before);std::vector<int> neighbours;
    for(const auto f:before.domain().faces(mesh)){
        if(f.horizontal||f.negativeCell<0||f.positiveCell<0||
            (f.negativeCell!=star.dryCell&&f.positiveCell!=star.dryCell))continue;
        const int ci=f.negativeCell==star.dryCell?f.positiveCell:f.negativeCell;
        if(ci==star.contact||std::find(neighbours.begin(),neighbours.end(),ci)!=neighbours.end())continue;
        neighbours.push_back(ci);
    }
    for(int ci:neighbours){
        const auto& old=before.plicCells()[ci];const auto& mid=moved.midpoint().plicCells()[ci];const auto& end=moved.endpoint().plicCells()[ci];
        if(old.inventory[0]==0)continue;
        check(mid.inventory[0]>0&&mid.resolved,"half sweep with actual remaining phase was falsely exhausted");
        check(end.inventory[0]==0&&end.resolved,"actual full whole-phase sweep retained a nongeometric stock tail");
        bool receipt=false;LAMP::Geometry::Detail::FiniteVolume::Dyadic ledger(old.inventory[0]);
        for(const auto f:before.domain().faces(mesh)){
            if(f.negativeCell!=ci&&f.positiveCell!=ci)continue;const auto& flux=moved.faces()[f.nativeRow].full;
            ledger.add(flux.liquid,f.negativeCell==ci?-1:1);
            const bool outgoing=flux.total>0?f.negativeCell==ci:f.positiveCell==ci;
            receipt|=outgoing&&flux.exhaustionLiquidCorrection!=0;
        }
        check(receipt,"full XY/RZ swept minority did not consume a real strict-geometry shared-flux certificate");
        check(ledger.compare(LAMP::Geometry::Detail::FiniteVolume::Dyadic(end.inventory[0]))==0,
            "whole-phase exhausted stock does not come from the exact same shared face ledger");
        ++wholePhaseControls;
    }
}
void advance(const QuadtreeMesh& mesh,Publication& accepted,unsigned count,bool requireChange)
{
    const auto velocity=requireChange?circulation(mesh,accepted):std::vector<double>(accepted.domain().nativeFaceCount());
    for(unsigned k=0;k<count;++k){
        const auto saved=accepted.plicCells();const auto id=accepted.identity();const auto oldClock=accepted.timeLayer();
        const auto oldPin=accepted.pinnedCorner();const auto contact=support(mesh,accepted).contact;
        const Clock end{oldClock.step+1,oldClock.time(),oldClock.time()+1./64,Position::End};
        auto mid=end;mid.position=Position::Midpoint;
        auto actualVelocity=velocity;if(k%2)for(auto& value:actualVelocity)value=-value;
        const auto slab=Slab::prepare(mesh,accepted,end,actualVelocity,actualVelocity);
        const auto moved=Transport::prepareLocalPlicSweep(mesh,accepted,slab,1240,1.184L);
        moved.validateFor(mesh,accepted,slab,1240,1.184L);
        check(moved.sourceIdentity()==id&&sameClock(moved.sourceTime(),oldClock)&&sameClock(moved.endTime(),end)
            &&moved.midpoint().identity()!=id&&moved.endpoint().identity()!=id
            &&moved.midpoint().identity()!=moved.endpoint().identity(),"pinned transport lost actual source/stage owners or complete clocks");
        wholePhaseReceipts(accepted,moved);
        wholePhaseStepControls(mesh,accepted,moved);
        inspect(mesh,moved.midpoint(),mid);inspect(mesh,moved.endpoint(),end);
        phaseConservation(accepted,moved.midpoint());phaseConservation(accepted,moved.endpoint());wallFluxes(mesh,accepted,moved);
        check(moved.report().halfRelativeConservation<=1e-12L&&moved.report().fullRelativeConservation<=1e-12L,
            "pinned transport weakened original half/full conservation gates");
        if(requireChange){
            pinFluxes(mesh,accepted,moved,actualVelocity);
            check(moved.midpoint().plicCells()[contact].inventory!=saved[contact].inventory
                &&moved.endpoint().plicCells()[contact].inventory!=saved[contact].inventory,
                "pin persistence used a stationary contact instead of actual changing stocks");
            check(moved.midpoint().plicCells()[contact].normal.x!=saved[contact].normal.x
                ||moved.midpoint().plicCells()[contact].normal.y!=saved[contact].normal.y,
                "pinned midpoint did not reconstruct geometry for its changed contact inventory");++changingSteps;
        }else for(std::size_t ci=0;ci<saved.size();++ci)
            check(moved.midpoint().plicCells()[ci].inventory==saved[ci].inventory&&moved.endpoint().plicCells()[ci].inventory==saved[ci].inventory,
                "zero-flow coincident pin changed independent stocks");
        rejects([&]{oldPin->validateFor(mesh,moved.midpoint());},"old pin cannot bind actual midpoint owner/clock");
        rejects([&]{oldPin->validateFor(mesh,moved.endpoint());},"old pin cannot bind actual endpoint owner/clock");
        check(accepted.identity()==id&&sameClock(accepted.timeLayer(),oldClock),"preparation changed accepted pin owner/clock");
        oldPin->validateFor(mesh,accepted);
        for(std::size_t ci=0;ci<saved.size();++ci)check(accepted.plicCells()[ci].inventory==saved[ci].inventory,
            "preparation or rejected stale certificate changed accepted pin stocks");
        accepted=moved.endpoint();++steps;
    }
}
Publication regrid(const QuadtreeMesh& sourceMesh,const Publication& source,const QuadtreeMesh& targetMesh,bool refine)
{
    const auto saved=source.plicCells();const auto oldPin=source.pinnedCorner();const auto time=source.timeLayer();
    check(time.position==Position::End,"pinned regrid is not at accepted End");
    const auto targetDomain=SelectedFluidTransportContract::selected(targetMesh,source.domain().axisymmetric());
    const auto transfer=GeometricVofMeshTransfer::prepare(sourceMesh,source,targetMesh,targetDomain,time);
    const auto candidate=transfer.candidate();inspect(targetMesh,candidate,time);phaseConservation(source,candidate);
    check(transfer.sourceIdentity()==source.identity()&&candidate.identity()!=source.identity()&&sameClock(candidate.timeLayer(),time)
        &&transfer.report().localPlic,"pinned actual regrid lost source/target/same-End ownership");
    if(refine){check(transfer.report().refinedDonors>0&&targetMesh.cellIndexNum()>sourceMesh.cellIndexNum(),"pinned refinement is a no-op");++refines;}
    else{check(transfer.report().coarsenedTargets>0&&targetMesh.cellIndexNum()<sourceMesh.cellIndexNum(),"pinned coarsening is a no-op");++coarsens;}
    for(unsigned phase=0;phase<2;++phase)check(transfer.report().relativeGlobalPhaseResidual[phase]<=1e-12L,
        "pinned mesh transfer weakened independent phase conservation");
    rejects([&]{oldPin->validateFor(targetMesh,candidate);},"old native pin certificate cannot bind regridded target");
    oldPin->validateFor(sourceMesh,source);
    for(std::size_t ci=0;ci<saved.size();++ci)check(source.plicCells()[ci].inventory==saved[ci].inventory,
        "same-End pin transfer or stale-certificate rejection changed source stocks");
    return candidate;
}
void runScene(bool rz)
{
    stage=std::string("pinned actual nozzle transport ")+(rz?"RZ":"XY");
    auto coarse=SelectedFluidTransportContract::nozzleMesh(false),fine=SelectedFluidTransportContract::nozzleMesh(true);
    const auto initial=GeometricVofInitialization::prepare(*coarse,SelectedFluidTransportContract::selected(*coarse,rz),
        [](Vector2D q){return q.x+q.y-4;});
    const auto s=support(*coarse,initial);const auto saved=initial.plicCells();
    auto accepted=Publication::withPinnedCorner(*coarse,initial,{1,3},s.wetRow,s.dryRow);
    check(!initial.pinnedCorner()&&accepted.identity()!=initial.identity()&&sameClock(accepted.timeLayer(),initial.timeLayer()),
        "adding pin mutated source or advanced its physical clock");
    for(std::size_t ci=0;ci<saved.size();++ci)check(accepted.plicCells()[ci].inventory==saved[ci].inventory,
        "adding a pin changed authoritative independent phase stocks");
    inspect(*coarse,accepted,initial.timeLayer());peakCflReject(*coarse,accepted);advance(*coarse,accepted,4,true);
    accepted=regrid(*coarse,accepted,*fine,true);advance(*fine,accepted,4,true);
    accepted=regrid(*fine,accepted,*coarse,false);advance(*coarse,accepted,4,true);
    check(accepted.timeLayer().step==12&&accepted.timeLayer().time()==12./64,"pinned regrid advanced physical time or step");
    const auto pure=GeometricVofInitialization::prepare(*coarse,SelectedFluidTransportContract::selected(*coarse,rz),
        [](Vector2D q){return q.x-1;});
    const auto ps=support(*coarse,pure);auto coincident=Publication::withPinnedCorner(*coarse,pure,{1,3},ps.wetRow,ps.dryRow);
    inspect(*coarse,coincident,pure.timeLayer());
    check(coincident.pinnedCorner()->kind()==Pin::Kind::NativeCoincident,"pure contact lost actual native coincident classification");
    advance(*coarse,coincident,1,false);
}
unsigned prescribedScenes=0,prescribedLayers=0,prescribedGuards=0,prescribedRemaps=0;
unsigned prescribedIncomingFaces=0,prescribedZeroWallRows=0;
template<class F>void prescribedReject(F action,const char* expected)
{
    bool caught=false;std::string actual;
    try{action();}catch(const std::invalid_argument& error){caught=true;actual=error.what();}
    catch(const std::exception& error){throw std::runtime_error(std::string("prescribed pin wrong exception: ")+error.what());}
    check(caught,"prescribed pin negative did not throw");
    if(actual.find(expected)==std::string::npos)throw std::runtime_error(std::string("prescribed pin expected {")+expected+"}, actual {"+actual+"}");
    ++prescribedGuards;
}
void prescribedInspect(const QuadtreeMesh& mesh,const Publication& p,Clock clock)
{
    p.validateFor(mesh,p.identity(),clock);const auto star=support(mesh,p);const auto pin=p.pinnedCorner();
    check(pin&&pin->wallTraceLaw()==Pin::WallTraceLaw::PrescribedWetDry&&pin->prescribedWetWallTrace()&&pin->prescribedDryWallTrace(),
        "current prescribed wall receipts/law are missing");pin->validateFor(mesh,p);
    check(pin->sourceIdentity()==p.identity()&&sameClock(pin->sourceTime(),clock)&&
        pin->nativeCell()==star.contact&&pin->nativeVertex()==star.vertex&&
        pin->wetNativeRow()==star.wetRow&&pin->dryNativeRow()==star.dryRow,
        "prescribed wall data lost its actual owner/full clock/target support");
    const auto& contact=p.plicCells()[star.contact];
    check(pin->kind()==Pin::Kind::ResolvedCorner&&contact.resolved&&contact.inventory[0]>0&&contact.inventory[1]>0,
        "prescribed wall law replaced the actual resolved mixed contact");
    const auto segment=Cut::interfaceSegment(contact.rectangle,contact.normal,contact.offset);
    const Point corner{1-(static_cast<long double>(contact.rectangle.lower.x)+contact.rectangle.upper.x)/2,
        3-(static_cast<long double>(contact.rectangle.lower.y)+contact.rectangle.upper.y)/2};
    const auto equal=[](Point a,Point b){return a.x==b.x&&a.y==b.y;};
    check(segment&&(equal((*segment)[0],corner)||equal((*segment)[1],corner)),
        "prescribed pin contact lost its exact native corner endpoint");
    for(unsigned phase=0;phase<2;++phase)near(contact.measured[phase],contact.inventory[phase],contact.inventory[phase],
        "prescribed contact independent actual phase measure",1e-12L);
    for(unsigned expected=0;expected<2;++expected){const int native=expected==0?star.wetCell:star.dryCell;const auto& cell=p.plicCells()[native];
        check(cell.inventory[expected]>0&&cell.measured[expected]==p.cells()[native].capacity&&cell.measured[1-expected]==0,
            "prescribed neighbour lacks the actual expected-majority geometry");
        if(cell.inventory[1-expected]>0)check(!cell.resolved&&cell.inventory[expected]>cell.inventory[1-expected]&&
            cell.normal.x==1&&cell.normal.y==0,"prescribed wall data falsely resolved or changed minority stocks");
    }
    for(int native=0;native<mesh.cellIndexNum();++native)if(p.domain().cellRowByIndex()[native]<0)
        check(p.plicCells()[native].inventory==std::array<long double,2>{0,0}&&p.cells()[native].capacity==0,
            "prescribed pin law filled excluded solid stock");
    ++prescribedLayers;
}
void prescribedConservation(const Publication& a,const Publication& b)
{
    const auto before=SelectedFluidTransportContract::inventories(a),after=SelectedFluidTransportContract::inventories(b);
    for(unsigned phase=0;phase<2;++phase)near(after[phase],before[phase],before[phase],
        "prescribed independent closed phase conservation",1e-12L);
}
void prescribedScene(bool rz)
{
    stage=std::string("explicit prescribed pin wall ")+(rz?"RZ":"XY");
    auto mesh=SelectedFluidTransportContract::nozzleMesh(false),fine=SelectedFluidTransportContract::nozzleMesh(true);
    const auto domain=SelectedFluidTransportContract::selected(*mesh,rz);
    const auto initial=GeometricVofInitialization::prepare(*mesh,domain,[](Vector2D q){return q.x+q.y-4;});
    const auto star=support(*mesh,initial);std::vector<std::array<long double,2>> stocks;
    for(const auto& cell:initial.plicCells())stocks.push_back(cell.inventory);
    const long double tiny=3.7684908769063313e-119L;
    auto micro=stocks;micro[star.wetCell][1]=tiny;micro[star.dryCell][0]=tiny;
    const auto raw=Publication::prepareLocalPlicInventory(*mesh,domain,initial.timeLayer(),micro,{});
    check(!raw.plicCells()[star.wetCell].resolved&&!raw.plicCells()[star.dryCell].resolved,
        "prescribed fixture did not publish actual unresolved neighbours");
    prescribedReject([&]{Publication::withPinnedCorner(*mesh,raw,{1,3},star.wetRow,star.dryRow);},
        "unambiguous pure wet-liquid and dry-gas wall neighbours");
    prescribedReject([&]{Publication::withPinnedCorner(*mesh,raw,{1,3},star.wetRow,star.dryRow,static_cast<Pin::WallTraceLaw>(99));},
        "invalid pinned wall trace law");
    const auto saved=raw.plicCells();const auto declared=Publication::withPinnedCorner(*mesh,raw,{1,3},star.wetRow,star.dryRow,Pin::WallTraceLaw::PrescribedWetDry);
    prescribedInspect(*mesh,declared,raw.timeLayer());
    check(!raw.pinnedCorner()&&declared.identity()!=raw.identity()&&sameClock(declared.timeLayer(),raw.timeLayer()),
        "prescribed cold entry changed original source or time");
    for(std::size_t i=0;i<saved.size();++i){check(declared.plicCells()[i].inventory==saved[i].inventory&&raw.plicCells()[i].inventory==saved[i].inventory,
        "prescribed wall declaration changed an independent stock");}
    for(int i:{star.wetCell,star.dryCell})check(declared.plicCells()[i].measured==saved[i].measured&&
        declared.plicCells()[i].resolved==saved[i].resolved&&declared.plicCells()[i].normal.x==saved[i].normal.x&&
        declared.plicCells()[i].normal.y==saved[i].normal.y&&declared.plicCells()[i].offset==saved[i].offset,
        "prescribed law fabricated neighbouring minor geometry");
    for(bool wrongMajor:{false,true}){
        auto bad=stocks;const auto capacity=initial.cells()[star.dryCell].capacity;
        bad[star.dryCell]=wrongMajor?std::array<long double,2>{capacity,tiny}:
            std::array<long double,2>{capacity/4,3*capacity/4};
        const auto p=Publication::prepareLocalPlicInventory(*mesh,domain,initial.timeLayer(),bad,{});
        prescribedReject([&]{Publication::withPinnedCorner(*mesh,p,{1,3},star.wetRow,star.dryRow,Pin::WallTraceLaw::PrescribedWetDry);},
            "prescribed pin wall requires current pure or unresolved expected-majority boundary-placeholder neighbour");
    }
    const auto pure=GeometricVofInitialization::prepare(*mesh,domain,[](Vector2D q){return q.x-1;});
    const auto ps=support(*mesh,pure);
    const auto purePin=Publication::withPinnedCorner(*mesh,pure,{1,3},ps.wetRow,ps.dryRow,Pin::WallTraceLaw::PrescribedWetDry);
    check(purePin.pinnedCorner()->kind()==Pin::Kind::NativeCoincident&&
        purePin.pinnedCorner()->wallTraceLaw()==Pin::WallTraceLaw::PrescribedWetDry&&
        !purePin.pinnedCorner()->prescribedWetWallTrace()&&!purePin.pinnedCorner()->prescribedDryWallTrace(),
        "NativeCoincident borrowed the prescribed-wall receipt");
    std::vector<std::array<long double,2>> pureMicro;for(const auto& cell:pure.plicCells())pureMicro.push_back(cell.inventory);
    pureMicro[ps.dryCell][0]=tiny;
    const auto notCoincident=Publication::prepareLocalPlicInventory(*mesh,domain,pure.timeLayer(),pureMicro,{});
    prescribedReject([&]{Publication::withPinnedCorner(*mesh,notCoincident,{1,3},ps.wetRow,ps.dryRow,Pin::WallTraceLaw::PrescribedWetDry);},
        "unambiguous pure wet-liquid and dry-gas wall neighbours");
    // A local vertex circulation on the far edge of the dry neighbour.
    // Its positive minor arrival is away from the pin's incident faces and
    // comes from the unchanged unresolved/P0 donor law, as in actual step54.
    int donor=-1,entry=-1;
    for(const auto face:domain.faces(*mesh))if(!face.horizontal&&face.negativeCell>=0&&face.positiveCell>=0&&
        (face.negativeCell==star.dryCell||face.positiveCell==star.dryCell)){
        const int other=face.negativeCell==star.dryCell?face.positiveCell:face.negativeCell;
        if(other!=star.contact){check(donor<0,"prescribed fixture has ambiguous dry far face");donor=other;entry=face.nativeRow;}}
    check(donor>=0&&stocks[donor][0]==0,"prescribed fixture has no independent pure-G far donor");
    stocks[donor][0]=tiny;const auto source=Publication::prepareLocalPlicInventory(*mesh,domain,initial.timeLayer(),stocks,{});
    const auto strict=Publication::withPinnedCorner(*mesh,source,{1,3},star.wetRow,star.dryRow);
    check(strict.pinnedCorner()->wallTraceLaw()==Pin::WallTraceLaw::StrictWholeCell&&
        !strict.pinnedCorner()->prescribedWetWallTrace()&&!strict.pinnedCorner()->prescribedDryWallTrace(),
        "default pin silently enabled prescribed wall data");
    const auto accepted=Publication::withPinnedCorner(*mesh,source,{1,3},star.wetRow,star.dryRow,Pin::WallTraceLaw::PrescribedWetDry);
    prescribedInspect(*mesh,accepted,source.timeLayer());const auto oldPin=accepted.pinnedCorner();const auto savedSource=accepted.plicCells();
    const auto face=domain.face(*mesh,entry);const auto endpoints=faceEndpoints(*mesh,face);const Point site=endpoints[1];
    const long double amplitude=source.cells()[donor].capacity/16;
    const auto psi=[&](Point q){return q.x==site.x&&q.y==site.y?-amplitude:0.L;};
    std::vector<double> velocity(domain.nativeFaceCount());
    for(const auto f:domain.faces(*mesh)){const auto e=faceEndpoints(*mesh,f);const auto q=(f.horizontal?-1.L:1.L)*(psi(e[1])-psi(e[0]));
        velocity[f.nativeRow]=f.area>0?static_cast<double>(q/f.area):0;
        if(f.area==0)check(q==0,"prescribed local circulation leaks at axis");
        if(f.negativeCell<0||f.positiveCell<0)check(velocity[f.nativeRow]==0,"prescribed local circulation leaks through a solid wall");}
    const Clock end{1,0,1,Position::End};auto mid=end;mid.position=Position::Midpoint;
    const auto strictSlab=Slab::prepare(*mesh,strict,end,velocity,velocity);
    prescribedReject([&]{Transport::prepareLocalPlicSweep(*mesh,strict,strictSlab,1240,1.184L);},
        "unambiguous pure wet-liquid and dry-gas wall neighbours");
    const auto slab=Slab::prepare(*mesh,accepted,end,velocity,velocity);
    const auto moved=Transport::prepareLocalPlicSweep(*mesh,accepted,slab,1240,1.184L);
    moved.validateFor(*mesh,accepted,slab,1240,1.184L);
    prescribedInspect(*mesh,moved.midpoint(),mid);prescribedInspect(*mesh,moved.endpoint(),end);
    prescribedConservation(accepted,moved.midpoint());prescribedConservation(accepted,moved.endpoint());
    check(moved.report().localUnresolvedDonors>0&&moved.report().halfRelativeConservation<=1e-12L&&
        moved.report().fullRelativeConservation<=1e-12L,"prescribed transport changed original P0/conservation gates");
    for(const auto* next:{&moved.midpoint(),&moved.endpoint()}){
        const auto& dry=next->plicCells()[star.dryCell];
        check(dry.inventory[0]>0&&!dry.resolved&&dry.measured[0]==0,
            "prescribed actual dry-neighbour minor was erased or falsely resolved");}
    const auto& flux=moved.faces()[entry].full;const auto incoming=face.positiveCell==star.dryCell?flux.liquid:-flux.liquid;
    const auto rawIncoming=face.positiveCell==star.dryCell?flux.rawLiquid:-flux.rawLiquid;
    check(incoming>0&&rawIncoming>0&&flux.exhaustionLiquidCorrection==0&&flux.limiter==1,
        "prescribed trace fixture did not consume a positive unchanged P0 shared transfer");++prescribedIncomingFaces;
    LAMP::Geometry::Detail::FiniteVolume::Dyadic dryLedger(accepted.plicCells()[star.dryCell].inventory[0]);
    for(const auto f:domain.faces(*mesh)){
        if(f.negativeCell==star.dryCell)dryLedger.add(moved.faces()[f.nativeRow].full.liquid,-1);
        if(f.positiveCell==star.dryCell)dryLedger.add(moved.faces()[f.nativeRow].full.liquid,1);
        if(f.negativeCell<0||f.positiveCell<0){for(const auto* value:{&moved.faces()[f.nativeRow].half,&moved.faces()[f.nativeRow].full})
            check(value->total==0&&value->liquid==0&&value->gas==0&&value->mass==0&&value->positiveTotal==0&&value->negativeTotal==0&&
                value->positiveLiquid==0&&value->negativeLiquid==0&&value->positiveGas==0&&value->negativeGas==0,
                "prescribed wall data injected a phase flux through solid");++prescribedZeroWallRows;}
    }
    check(dryLedger.rounded()==moved.endpoint().plicCells()[star.dryCell].inventory[0],
        "positive dry microphase does not come from the unchanged shared ledger");
    prescribedReject([&]{oldPin->validateFor(*mesh,moved.midpoint());},"");
    prescribedReject([&]{oldPin->validateFor(*mesh,moved.endpoint());},"");
    oldPin->validateFor(*mesh,accepted);
    for(std::size_t i=0;i<savedSource.size();++i)check(accepted.plicCells()[i].inventory==savedSource[i].inventory,
        "prescribed preparation/rejected stale certificate mutated accepted stocks");
    const auto refined=GeometricVofMeshTransfer::prepare(*mesh,moved.endpoint(),*fine,SelectedFluidTransportContract::selected(*fine,rz),end);
    check(refined.report().refinedDonors>0,"prescribed same-End refinement is not actual");
    prescribedInspect(*fine,refined.candidate(),end);prescribedConservation(moved.endpoint(),refined.candidate());++prescribedRemaps;
    prescribedReject([&]{moved.endpoint().pinnedCorner()->validateFor(*fine,refined.candidate());},"");
    const auto coarse=GeometricVofMeshTransfer::prepare(*fine,refined.candidate(),*mesh,domain,end);
    check(coarse.report().coarsenedTargets>0,"prescribed same-End coarsening is not actual");
    prescribedInspect(*mesh,coarse.candidate(),end);prescribedConservation(refined.candidate(),coarse.candidate());++prescribedRemaps;
    prescribedReject([&]{refined.candidate().pinnedCorner()->validateFor(*mesh,coarse.candidate());},"");
    check(coarse.candidate().plicCells()[star.dryCell].inventory[0]>0&&!coarse.candidate().plicCells()[star.dryCell].resolved,
        "same-End AMR erased the prescribed-neighbour independent microphase");
    ++prescribedScenes;
}
void run()
{
    for(bool rz:{false,true})runScene(rz);
    for(bool rz:{false,true})prescribedScene(rz);
    check(prescribedScenes==2&&prescribedLayers==12&&prescribedRemaps==4&&prescribedIncomingFaces==2&&prescribedZeroWallRows>0&&prescribedGuards>=20,
        "prescribed pin wall stock/transport/target requalification coverage incomplete");
    std::cout<<"PRESCRIBED_PIN_WALL_PASS scenes="<<prescribedScenes<<" current_layers="<<prescribedLayers
        <<" same_End_remaps="<<prescribedRemaps<<" positive_minor_entry_faces="<<prescribedIncomingFaces
        <<" zero_wall_rows="<<prescribedZeroWallRows<<" guards="<<prescribedGuards
        <<" stocks_unchanged_by_law=1 unresolved_geometry_retained=1 no_internal_trace_claim=1\n";
    check(steps==26&&changingSteps==24&&refines==2&&coarsens==2&&coincidentLayers==6&&wallChecks>0&&solidChecks>0
        &&pinnedSweeps==96&&positivePinFaces>0&&negativePinFaces>0&&pinPeakRejections==2
        &&wholePhaseExhaustions>0&&wholePhaseControls==2,
        "pinned XY/RZ changing transport, same-End AMR or pure-coincident coverage incomplete");
    std::cout<<std::setprecision(17)<<"PINNED_CORNER_TRANSPORT_PASS steps="<<steps<<" changing_steps="<<changingSteps
        <<" stage_pin_checks="<<layerChecks<<" same_End_refines="<<refines<<" same_End_coarsens="<<coarsens
        <<" native_coincident_layers="<<coincidentLayers<<" wall_checks="<<wallChecks<<" solid_checks="<<solidChecks
        <<" pin_polygon_sweeps="<<pinnedSweeps<<" positive_pin_faces="<<positivePinFaces<<" negative_pin_faces="<<negativePinFaces
        <<" peak_pin_Courant="<<peakPinnedCfl<<" peak_Courant_rejections="<<pinPeakRejections
        <<" whole_phase_exhaustions="<<wholePhaseExhaustions<<" whole_phase_controls="<<wholePhaseControls
        <<" own_phase_exhaustion_correction="<<peakExhaustionCorrection
        <<" phase_conservation="<<peakPhaseResidual<<" contact_geometry="<<peakGeometryResidual<<" psi_divergence="<<peakDivergence
        <<" actual_Interface_transport_no_Hydro_HF_Electric_or_App_claim\n";
}
}
namespace SelectedFluidHydroContract
{
using namespace SelectedFluidTransportContract;
void run()
{
    unsigned acceptedSteps=0,compactFaces=0,internalWalls=0;long double peakStoredR=0,peakEnergy=0;
    for(bool rz:{false,true})for(bool amr:{false,true})for(bool open:{false,true}){
        stage=std::string("selected Hydro ")+(rz?"RZ":"XY")+(amr?" 2:1":" uniform")+(open?" reservoir":" closed");
        DictionaryMesh meshSpec;meshSpec.Nx=20;meshSpec.Ny=64;meshSpec.minCoord={0,0};meshSpec.maxCoord={5,16};
        meshSpec.maxCellNum=16384;meshSpec.maxInitialIterations=0;meshSpec.maxLevel_local={2};
        auto mesh=std::make_unique<QuadtreeMesh>(meshSpec);
        if(amr){DictionaryMesh::RefineBox b;b.level=1;b.min={0,3};b.max={1,5};mesh->applyRefineBoxes({b});}
        const auto fluid=selected(*mesh,rz);
        const auto initial=GeometricVofInitialization::prepare(*mesh,fluid,[](Vector2D p){return std::max(p.x-1,p.y-4.125);});
        const auto metrics=HydroPhysicalDivergence::build(*mesh,HydroIntegralDivergence::build(*mesh,HydroDofLayout::build(*mesh,fluid)));
        const auto& layout=metrics.layout();
        check(static_cast<std::size_t>(layout.pressureCount())==fluid.cellIndexByRow().size()&&layout.pressureCount()<mesh->cellIndexNum()
            &&layout.velocityCount()<fluid.nativeFaceCount(),"selected Hydro retained a solid pressure/velocity unknown");
        const auto psi=[](Point p){if(p.x>=1||p.y<=3||p.y>=5)return 0.L;
            return p.x*p.x*(1-p.x)*(1-p.x)*(p.y-3)*(p.y-3)*(5-p.y)*(5-p.y)/256;};
        std::vector<double> velocity(layout.velocityCount());HydroVofPrescribedFlow flow;
        flow.rulesByNativeRow.assign(fluid.nativeFaceCount(),Mode::Closed);flow.normalVelocity.resize(layout.velocityCount());
        flow.endpointVelocity=[](int,Clock){return Vector2D(0,1);};
        flow.pressureOutletBackflowModel=HydroVofPressureOutletMomentumModel::InteriorContinuationZeroShear;
        for(const auto f:fluid.faces(*mesh)){
            const int u=layout.velocityDof(f.nativeRow);
            if(f.negativeCell<0&&f.positiveCell<0){check(u<0,"solid face has native Hydro velocity");++compactFaces;continue;}
            check(u>=0&&layout.nativeFace(u)==f.nativeRow,"selected Hydro scatter is not invertible");
            const auto ends=faceEndpoints(*mesh,f);const bool interior=f.negativeCell>=0&&f.positiveCell>=0;
            if(interior)flow.rulesByNativeRow[f.nativeRow]=Mode::Interior;
            if(open){
                if(f.horizontal&&ends[0].x<1&&ends[1].x<=1)velocity[u]=1;
                if(!interior&&f.horizontal&&f.center.y==0){flow.rulesByNativeRow[f.nativeRow]=Mode::PureLiquidReservoir;
                    flow.prescribedVelocityDofs.push_back(u);flow.normalVelocity[u]=1;}
                if(!interior&&f.horizontal&&f.center.y==16){flow.rulesByNativeRow[f.nativeRow]=Mode::PureGasReservoir;
                    flow.pressureOutlets.push_back({u,0});}
            }else{
                const auto q=(f.horizontal?-1.L:1.L)*(psi(ends[1])-psi(ends[0]));
                velocity[u]=f.area>0?static_cast<double>(q/f.area):0;
            }
        }
        const HydroNativeValues values(layout,velocity,std::vector<double>(layout.pressureCount()));
        const auto phase=HydroPhaseGeometry::build(*mesh,metrics,initial);
        auto accepted=open?HydroTwoPhaseState::createVofPrescribedFlow(*mesh,phase,values,flow):HydroTwoPhaseState::createVof(*mesh,phase,values);
        for(const auto f:fluid.faces(*mesh))if((f.horizontal&&f.center.y==3&&f.center.x>1)||(!f.horizontal&&f.center.x==1&&f.center.y<3)){
            const auto u=layout.velocityDof(f.nativeRow);if(u<0)continue;
            check(std::find(accepted.boundary().velocityDofs().begin(),accepted.boundary().velocityDofs().end(),u)
                !=accepted.boundary().velocityDofs().end(),"actual electrode wall missing from impermeable Hydro boundary");++internalWalls;
        }
        HydroVofAcConfig config;config.liquidDensity=2;config.gasDensity=1;
        config.liquidViscosity=.01;config.gasViscosity=.005;config.surfaceTension=0;
        for(unsigned k=0;k<2;++k){
            const auto rates=open?HydroVofStep::stabilityRatesPrescribedFlow(*mesh,accepted,config):HydroVofStep::stabilityRates(*mesh,accepted,config);
            check(rates.flowRate>0&&!rates.capillaryAvailable,"selected Hydro stability lost actual flow or added sigma0 capillarity");
            const auto before=accepted;const auto t=accepted.timeLayer();const Clock end{t.step+1,t.time(),t.time()+1./1024,Position::End};
            const auto step=open?HydroVofStep::preparePrescribedFlow(*mesh,accepted,end,config):HydroVofStep::prepare(*mesh,accepted,end,config);
            check(accepted.sameIdentity(before)&&step.report().linearSolves>=1&&step.report().linearSolves<=2,
                "selected Hydro changed accepted owner during preparation or added solves");
            inspect(*mesh,accepted.phaseGeometry().vofPublication(),step.transport());
            ConservativeVofMomentumContract::qualify(*mesh,step.frozenMidpoint().material());
            peakStoredR=std::max(peakStoredR,step.frozenMidpoint().material().report().maxRelativeStoredMassContinuityDefect);
            for(const auto& report:step.report().stages){
                check(report.consistencyPower==0&&report.momentumPowerRelativeDefect<=config.momentumPowerTolerance
                    &&report.energy.closureRelativeResidual<=1e-11L,"selected Hydro weakened original conservative momentum/power/KE gates");
                peakEnergy=std::max(peakEnergy,report.energy.closureRelativeResidual);
            }
            check(step.candidateState().sampleCellCenterVelocities(*mesh).size()==fluid.cellIndexByRow().size(),
                "selected Hydro cell sample is not compact pressure order");
            step.commit(*mesh,accepted);check(sameClock(accepted.timeLayer(),end),"selected Hydro commit lost actual full clock");++acceptedSteps;
        }
    }
    check(acceptedSteps==16&&compactFaces>0&&internalWalls>0,"selected Hydro scene/compact/wall coverage incomplete");
    std::cout<<std::setprecision(17)<<"SELECTED_FLUID_HYDRO_PASS accepted_steps="<<acceptedSteps<<" excluded_faces="<<compactFaces
        <<" actual_walls="<<internalWalls<<" max_stored_dual_R="<<peakStoredR<<" max_energy_closure="<<peakEnergy<<'\n';
}
}
namespace SelectedFluidEhdStepContract
{
namespace E=LAMP::Electric;
void run()
{
    unsigned steps=0,stages=0;long double maximumForce=0,maximumPower=0,maximumGauss=0;
    for(bool rz:{false,true}){
        stage=std::string("selected actual EHD step ")+(rz?"RZ":"XY");
        DictionaryMesh d;d.Nx=20;d.Ny=64;d.minCoord={0,0};d.maxCoord={5,16};
        d.maxCellNum=16384;d.maxInitialIterations=0;d.maxLevel_local={2};QuadtreeMesh mesh(d);
        DictionaryMesh::RefineBox box;box.level=1;box.min={0,3};box.max={1,5};mesh.applyRefineBoxes({box});
        const auto fluid=SelectedFluidTransportContract::selected(mesh,rz);
        const auto source=GeometricVofInitialization::prepare(mesh,fluid,[](Vector2D p){
            const auto z=std::max(p.y-3.,0.);return p.x*p.x+z*z-1.;});
        const auto metrics=HydroPhysicalDivergence::build(mesh,HydroIntegralDivergence::build(mesh,HydroDofLayout::build(mesh,fluid)));
        const auto& layout=metrics.layout();const HydroNativeValues rest(layout,std::vector<double>(layout.velocityCount()),std::vector<double>(layout.pressureCount()));
        auto accepted=HydroTwoPhaseState::createVof(mesh,HydroPhaseGeometry::build(mesh,metrics,source),rest);
        auto charge=E::ElectricVofCharge::fromIntegrated(mesh,source,std::vector<double>(mesh.cellIndexNum()));
        const E::ElectricVofMaterials materials{2,1,.3,0};const E::ElectricVofControls controls;
        std::vector<E::ElectricVofBoundary> bc(fluid.nativeFaceCount());
        for(const auto f:fluid.faces(mesh)){
            if((f.negativeCell<0&&f.positiveCell<0)||(f.negativeCell>=0&&f.positiveCell>=0))continue;
            const bool electrode=(f.horizontal&&f.center.y==3&&f.center.x>1)||(!f.horizontal&&f.center.x==1&&f.center.y<3);
            if(electrode)bc[f.nativeRow]={E::ElectricVofBoundary::Type::Dirichlet,6};
            else if(f.horizontal&&f.center.y==16)bc[f.nativeRow]={E::ElectricVofBoundary::Type::Dirichlet,0};
            else bc[f.nativeRow]={E::ElectricVofBoundary::Type::HomogeneousNeumann,0};
        }
        const auto initialField=E::ElectricVofField::prepareElectrostaticSnapshot(mesh,source,charge,materials,bc,controls);
        check(initialField.chargeCandidate().sameInstance(charge),"selected coupled Initial snapshot advanced Q");
        HydroVofAcConfig config;config.liquidDensity=2;config.gasDensity=1;
        config.liquidViscosity=.01;config.gasViscosity=.005;config.surfaceTension=0;config.requireMaxwellCorrection=true;
        for(unsigned k=0;k<2;++k){
            const auto acceptedCharge=charge;const auto before=accepted;
            config.maxwellLoad=[acceptedCharge,materials,controls,bc](const QuadtreeMesh& grid,const Transport& transport,
                const HydroPhaseGeometry& mid,const HydroPressureLoad& pressure){
                const auto moved=E::ElectricVofChargeTransport::prepareClosed(grid,acceptedCharge,transport,controls);
                const auto field=E::ElectricVofField::prepareOhmic(grid,mid.vofPublication(),moved.midpoint(),materials,bc,controls);
                return HydroMaxwellLoad::prepare(grid,pressure,mid,E::ElectricVofMaxwellStress::prepare(grid,field));
            };
            const auto t=accepted.timeLayer();const Clock end{t.step+1,t.time(),t.time()+1./4096,Position::End};
            const auto step=HydroVofStep::prepare(mesh,accepted,end,config);
            const auto moved=E::ElectricVofChargeTransport::prepareClosed(mesh,acceptedCharge,step.transport(),controls);
            const auto field=E::ElectricVofField::prepareOhmic(mesh,step.transport().endpoint(),moved.endpoint(),materials,bc,controls);
            check(accepted.sameIdentity(before)&&charge.sameInstance(acceptedCharge),"selected EHD prepare changed accepted fluid/Q separately");
            check(step.frozenMidpoint().hasMaxwellLoad()&&step.report().maxwellCorrectionChecked,"selected EHD omitted actual Maxwell correction");
            for(double f:step.frozenMidpoint().maxwellLoads())maximumForce=std::max(maximumForce,static_cast<long double>(std::abs(f)));
            for(const auto& report:step.report().stages){
                check(report.maxwellLoadReportAvailable&&report.maxwellLoadReport.volumePoints>0&&report.maxwellLoadReport.boundaryPoints>0,
                    "selected EHD lacks actual volume and wall stress traction");
                check(report.consistencyPower==0&&report.momentumPowerRelativeDefect<=config.momentumPowerTolerance
                    &&report.energy.closureRelativeResidual<=1e-11L,
                    "selected EHD violated retained physical momentum/power/energy gate");
                maximumPower=std::max(maximumPower,static_cast<long double>(report.momentumPowerRelativeDefect));++stages;
            }
            maximumGauss=std::max(maximumGauss,static_cast<long double>(field.report().maximumRelativeGaussResidual));
            check(field.report().maximumRelativeGaussResidual<=1e-10&&sameClock(field.chargeCandidate().timeLayer(),end),
                "selected EHD retained actual Gauss/full clock");
            for(std::size_t ci=0;ci<field.potential().size();++ci)if(fluid.cellRowByIndex()[ci]<0)
                check(field.potential()[ci]==0&&field.chargeCandidate().integratedCharge()[ci]==0
                    &&step.transport().endpoint().plicCells()[ci].inventory==std::array<long double,2>{0,0},"selected EHD solid gained field/Q/material");
            ConservativeVofMomentumContract::qualify(mesh,step.frozenMidpoint().material());
            step.commit(mesh,accepted);charge=field.chargeCandidate();charge.validateFor(mesh,accepted.timeLayer());++steps;
        }
    }
    check(steps==4&&stages>=steps&&maximumForce>0,"selected EHD real force/stage coverage missing");
    std::cout<<std::setprecision(17)<<"SELECTED_FLUID_EHD_STEP_PASS accepted_steps="<<steps<<" actual_stages="<<stages
        <<" maximum_force="<<maximumForce<<" momentum_power_defect="<<maximumPower<<" Gauss="<<maximumGauss
        <<" sigma0_positive_mu_K_actual_internal_electrode_no_pinned_or_AMR_claim\n";
}
}

namespace VofFilteredLocalUpdateContract
{
namespace FV=LAMP::Geometry::Detail::FiniteVolume;
using Filter=LAMP::Interface::Detail::VofLocalFilteredDyadic;
using Path=LAMP::Interface::Detail::VofLocalPhaseUpdatePath;
#if defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
struct State {
    std::fenv_t native{};unsigned csr=0;int error=0;
    void capture(){error=errno;csr=_mm_getcsr();check(std::fegetenv(&native)==0,"filtered FP capture");errno=error;}
    void restore()const{check(std::fesetenv(&native)==0,"filtered FP restore");_mm_setcsr(csr);errno=error;}
    bool same(const State& x)const{return error==x.error&&csr==x.csr&&std::memcmp(&native,&x.native,sizeof native)==0;}
};
struct Restore {State state;Restore(){state.capture();}~Restore(){state.restore();}};
struct Result {std::optional<FV::Rounded> rounded;std::string failure;State exit;Path path=Path::Noop;};
unsigned oracleCases=0,originalFailures=0,fast=0,promoted=0,fallback=0,noops=0;
void one(long double old,const std::vector<long double>& terms,const std::vector<long double>& compression,
    std::optional<Path> expected=std::nullopt)
{
    State entry;entry.capture();
    const auto reference=[&]{entry.restore();Result r;
        try{FV::Dyadic n,c;for(auto x:terms)n.add(x);for(auto x:compression)c.add(x);
            r.rounded=FV::phaseUpdate(old,n,c);r.exit.capture();}
        catch(const std::exception& error){r.exit.capture();r.failure=error.what();}return r;};
    const auto actual=[&]{entry.restore();Result r;
        try{Filter n,c;for(auto x:terms)n.add(x);for(auto x:compression)c.add(x);
            r.rounded=LAMP::Interface::Detail::vofLocalPhaseUpdate(old,n,c,r.path);r.exit.capture();}
        catch(const std::exception& error){r.exit.capture();r.failure=error.what();}return r;};
    const auto ref=reference(),got=actual();++oracleCases;
    check(ref.exit.same(got.exit),"filtered versus Core complete native FP/rawCSR/errno");
    check(ref.failure==got.failure&&ref.rounded.has_value()==got.rounded.has_value(),"filtered changed original first what");
    if(ref.rounded){
        check(std::memcmp(&ref.rounded->value,&got.rounded->value,sizeof(long double))==0,"filtered Core phaseUpdate value bits/sign");
        check(ref.rounded->error.compare(got.rounded->error)==0,"filtered Core exact rounding error");
        if(expected)check(got.path==*expected,"filtered actual phaseUpdate path");
        switch(got.path){case Path::Fast:++fast;break;case Path::Promoted:++promoted;break;
            case Path::Fallback:++fallback;break;case Path::Noop:++noops;break;}
    }else ++originalFailures;
    entry.restore();
}
void replay(const QuadtreeMesh& mesh,const Publication& source,const Transport& actual)
{
    const auto n=source.plicCells().size();
    for(bool half:{true,false}){
        std::vector<std::array<long double,2>> stocks;for(const auto& c:source.plicCells())stocks.push_back(c.inventory);
        auto defects=source.capacityDefects();const bool firstHorizontal=(actual.endTime().step%2)==0;
        for(unsigned direction=0;direction<2;++direction){
            const bool horizontal=direction==0?firstHorizontal:!firstHorizontal;
            std::vector<std::array<FV::Dyadic,2>> net(n);std::vector<FV::Dyadic> total(n);
            for(const auto& face:actual.faces()){
                if(source.domain().face(mesh,face.nativeRow).horizontal!=horizontal)continue;
                const auto& flux=half?face.half:face.full;
                for(int ci:{face.negativeCell,face.positiveCell})if(ci>=0){const int sign=ci==face.negativeCell?1:-1;
                    net[ci][0].add(flux.liquid,sign);net[ci][1].add(flux.gas,sign);total[ci].add(flux.total,sign);}
            }
            for(std::size_t i=0;i<n;++i){if(source.domain().cellRowByIndex()[i]<0)continue;
                const auto tracked=source.plicCells()[i].inventory[0]>source.plicCells()[i].inventory[1]?1u:0u;
                defects[i].add(net[i][0],-1);defects[i].add(net[i][1],-1);defects[i].add(total[i]);
                for(unsigned phase=0;phase<2;++phase){const auto c=phase==tracked?FV::Dyadic{}:total[i];
                    const auto update=FV::phaseUpdate(stocks[i][phase],net[i][phase],c);
                    stocks[i][phase]=update.value;defects[i].add(update.error);}
            }
        }
        const auto& candidate=half?actual.midpoint():actual.endpoint();
        for(std::size_t i=0;i<n;++i){for(unsigned phase=0;phase<2;++phase)
            check(std::memcmp(&stocks[i][phase],&candidate.plicCells()[i].inventory[phase],sizeof(long double))==0,
                "filtered actual shared-face Core replay independent phase stock bits");
            check(defects[i].compare(candidate.capacityDefects()[i])==0,"filtered actual independent predecessor exact defect");}
    }
}
void actualTransport()
{
    unsigned scenes=0;std::size_t calls=0,coreCalls=0;
    for(bool rz:{false,true})for(bool adaptive:{false,true}){
        std::unique_ptr<QuadtreeMesh> mesh;
        if(adaptive)mesh=vofMaterialMesh(4);else{DictionaryMesh d;d.Nx=d.Ny=4;d.maxCellNum=4096;d.maxInitialIterations=0;
            mesh=std::make_unique<QuadtreeMesh>(d);}
        const Clock initial{0,0,0,Position::Initial};
        const auto source=Publication::prepareLocalPlic(*mesh,domain(*mesh,rz),initial,l1PlaneAlpha(*mesh,rz,0));
        const auto a=vofMaterialNormals(*mesh,source,0),b=a;
        const double dt=static_cast<double>(.1L/vofMaterialOutflow(*mesh,source,a,b));
        const Clock end{1,0,dt,Position::End};const auto slab=Slab::prepare(*mesh,source,end,a,b);
        const auto oldStocks=source.plicCells();const auto id=source.identity();
        const auto t=Transport::prepareLocalPlicSweep(*mesh,source,slab,100,1);replay(*mesh,source,t);
        const auto& r=t.report();
        check(r.localPhaseUpdateCalls==r.localPhaseUpdateFast+r.localPhaseUpdatePromoted+r.localPhaseUpdateFallback+r.localPhaseUpdateNoop,
            "Core actual returned-call receipt");
        check(r.localPhaseUpdateCalls==8*source.plicCells().size()&&r.localPhaseUpdateFallback==r.localPhaseUpdateCalls&&
            r.localPhaseUpdateFast==0&&r.localPhaseUpdatePromoted==0&&r.localPhaseUpdateNoop==0,
            "actual production did not return through original Core phaseUpdate exclusively");
        check(source.identity()==id&&source.plicCells().size()==oldStocks.size(),"filtered actual changed original owner");
        for(std::size_t i=0;i<oldStocks.size();++i)check(source.plicCells()[i].inventory==oldStocks[i].inventory,"filtered actual modified source inventory");
        calls+=r.localPhaseUpdateCalls;coreCalls+=r.localPhaseUpdateFallback;++scenes;
    }
    check(scenes==4,"filtered actual XY/RZ regular/2:1 witness count");
    std::cout<<"LOCAL_PLIC_CORE_ACTUAL scenes="<<scenes<<" phase_calls="<<calls<<" core_calls="<<coreCalls
        <<" filtered_production=disabled shared_face_Core_stock_bits=exact predecessor_defects=exact\n";
}
void run()
{
    Restore restore;
    check(sizeof(long double)==8&&std::numeric_limits<long double>::digits==53,"filtered contract requires actual MSVC binary64");
    // Clear complete x87 + raw SSE flags for a sensitive independent suffix oracle.
    auto entry=restore.state;entry.native._Fe_stat=0;entry.csr&=~0x3fu;entry.error=173;entry.restore();
    const auto p=[](int e){return std::scalbn(1.L,e);};
    one(-0.L,{}, {},Path::Noop);
    one(-0.L,{1.L},{1.L},Path::Fallback); // nonnoop exact zero must run original environment and return +0.
    one(.75L,{p(-5),-p(-7),p(-9)},{p(-10)},Path::Fast);
    one(1.L,{p(-54)},{},Path::Fallback); // lower-binade half tie: original quantum/suffix.
    one(1.5L,{-p(-53)},{},Path::Fast); // same-binade genuine RN-even tie.
    one(1.L,{-p(-200),-p(-400)},{},Path::Promoted); // same rounded stock; two-component error must survive.
    one(1.L,{-p(-1000)},{},Path::Fast);
    one(.75L,{1.L,-1.L},{},Path::Noop);
    one(0.L,{1.L,-1.L},{p(-54)},Path::Fast);
    const auto tiny=std::numeric_limits<long double>::min();one(tiny,{tiny/2},{},Path::Fallback);
    const auto maximum=std::numeric_limits<long double>::max();
    one(maximum,{maximum,maximum,-maximum},{maximum},Path::Fallback); // overflowed EFT intermediate, exact final remains finite.
    one(maximum,{-maximum},{}); // original final-overflow what/FP, no new helper overflow.
    one(-1.L,{},{});one(std::numeric_limits<long double>::quiet_NaN(),{},{});
    one(1.L,{std::numeric_limits<long double>::infinity()},{});
    // Lexicographic fast comparison is checked against the existing exact oracle.
    Filter x,y;FV::Dyadic dx,dy;
    for(auto term:{1.L,p(-200)}){x.add(term);dx.add(term);}
    for(auto term:{1.L,p(-199)}){y.add(term);dy.add(term);}
    check(x.isFast()&&y.isFast()&&x.compare(y)==dx.compare(dy),"filtered equal-head exact tail comparison");
    Filter third=x;third.add(p(-400));dx.add(p(-400));
    check(!third.isFast()&&third.dyadic().compare(dx)==0,"filtered third nonzero tail must permanently promote");
    third.add(-1.L);dx.add(-1.L);check(!third.isFast()&&third.dyadic().compare(dx)==0,"filtered promoted object returned to fast or lost tail");
    const int oldRound=std::fegetround();check(std::fesetround(FE_UPWARD)==0,"filtered directed rounding setup");
    one(.75L,{p(-5)},{p(-10)});check(std::fesetround(oldRound)==0,"filtered directed rounding restore");
    const unsigned csr=_mm_getcsr();_mm_setcsr(csr|0x8000u);one(.75L,{p(-5)},{p(-10)});_mm_setcsr(csr);
    _mm_setcsr(csr|0x40u);one(.75L,{p(-5)},{p(-10)});_mm_setcsr(csr);
    check(fast>=4&&promoted>=1&&fallback>=4&&noops==2&&originalFailures>=4,"filtered original oracle path coverage");
    actualTransport();
    std::cout<<"LOCAL_PLIC_FILTERED_HELPER_REFERENCE_PASS oracle_cases="<<oracleCases<<" fast="<<fast<<" promoted="<<promoted
        <<" fallback="<<fallback<<" noop="<<noops<<" original_failures="<<originalFailures
        <<" historical_helper_only=1 production_uses_helper=0 Core_value_bits_and_exact_error=exact full_native_rawCSR_errno=exact no_speed_claim\n";
}
#else
void run(){throw std::runtime_error("filtered contract requires actual MSVC x64");}
#endif
}

namespace MaterialFullOnlyContract
{
using Mode=Slab::Boundary::Mode;
unsigned scenes=0,nativeQueries=0,interiorQueries=0,materialAtoms=0,negativeControls=0,limiterDifferences=0;
void sameFlux(const GeometricVofIntegratedFlux& a,const GeometricVofIntegratedFlux& b)
{
    check(F15BoundaryGrossContract::sameFields(a,b),"full-only changed original full scalar bits or zero signs");
    for(auto field:{&GeometricVofIntegratedFlux::exhaustionLiquidCorrection,&GeometricVofIntegratedFlux::exhaustionGasCorrection})
        check(a.*field==b.*field&&std::signbit(a.*field)==std::signbit(b.*field),"full-only changed original exhaustion receipt");
}
template<class Old,class Full>void sameFailure(const Old& old,const Full& full)
{
    std::string before,after;bool oldType=false,fullType=false;
    try{old();}catch(const std::invalid_argument& error){oldType=true;before=error.what();}
    try{full();}catch(const std::invalid_argument& error){fullType=true;after=error.what();}
    check(oldType&&fullType&&!before.empty()&&before==after,"full-only changed common original query guard/first what");++negativeControls;
}
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
void parallelInterior(QuadtreeMesh& mesh,const Transport& t,const std::vector<GeometricVofInteriorSubface>& original)
{
    using namespace F51LocalFaceParallelContract;
    Runtime restore;omp_set_dynamic(0);unsigned control=0;
    check(_controlfp_s(&control,_MCW_EM|_RC_NEAR,_MCW_EM|_MCW_RC)==0,"Coons masked RNE setup");
    check(std::fesetround(FE_TONEAREST)==0,"Coons native RNE setup");_mm_setcsr((_mm_getcsr()|0x1f80u)&~0xe040u);
    Entry entry;auto queries=original;check(!queries.empty(),"Coons original owner queries missing");
    // Exercise interleaved/repeated owners and preserve caller result order.
    std::reverse(queries.begin(),queries.end());
    while(queries.size()<128){const auto copy=queries;queries.insert(queries.end(),copy.begin(),copy.end());}
    omp_set_num_threads(1);entry.restore();const auto serial=t.measureInteriorFaces(mesh,queries);Receipt pairedSerial;pairedSerial.capture();
    omp_set_num_threads(8);entry.restore();const auto parallel=t.measureInteriorFaces(mesh,queries);Receipt pairedParallel;pairedParallel.capture();
    check(pairedSerial.same(pairedParallel),"Coons paired native FP/rawCSR/errno changed");
    omp_set_num_threads(1);entry.restore();const auto serialFull=t.measureInteriorFacesFull(mesh,queries);Receipt fullSerial;fullSerial.capture();
    omp_set_num_threads(8);entry.restore();const auto parallelFull=t.measureInteriorFacesFull(mesh,queries);Receipt fullParallel;fullParallel.capture();
    check(fullSerial.same(fullParallel),"Coons full native FP/rawCSR/errno changed");
    check(serial.size()==parallel.size()&&serialFull.size()==parallelFull.size(),"Coons owner gather result extent");
    for(std::size_t i=0;i<queries.size();++i){
        check(serial[i].nativeCell==parallel[i].nativeCell&&serialFull[i].nativeCell==parallelFull[i].nativeCell,
            "Coons owner gather changed caller order");
        sameFlux(serial[i].half,parallel[i].half);sameFlux(serial[i].full,parallel[i].full);sameFlux(serialFull[i].full,parallelFull[i].full);
        for(const auto pair:{std::pair<long double,long double>{serial[i].beginVelocity.x,parallel[i].beginVelocity.x},
            {serial[i].beginVelocity.y,parallel[i].beginVelocity.y},{serial[i].endVelocity.x,parallel[i].endVelocity.x},
            {serial[i].endVelocity.y,parallel[i].endVelocity.y}})
            check(pair.first==pair.second&&std::signbit(pair.first)==std::signbit(pair.second),"Coons paired normal bits/sign changed");
    }
    auto bad=queries;bad[1].end=bad[1].begin;bad[2].nativeCell=-1;
    rejectedPair([&]{(void)t.measureInteriorFacesFull(mesh,bad);},"local interior face must be nondegenerate");
    // Directed rounding is ineligible: preserve the literal route and status.
    check(std::fesetround(FE_UPWARD)==0,"Coons fallback directed setup");Entry directed;
    std::vector<GeometricVofInteriorFullFlux> directedSerial,directedParallel;std::string se,pe;
    omp_set_num_threads(1);directed.restore();try{directedSerial=t.measureInteriorFacesFull(mesh,queries);}
    catch(const std::invalid_argument& e){se=e.what();}Receipt ds;ds.capture();
    omp_set_num_threads(8);directed.restore();try{directedParallel=t.measureInteriorFacesFull(mesh,queries);}
    catch(const std::invalid_argument& e){pe=e.what();}Receipt dp;dp.capture();
    check(se==pe&&ds.same(dp),"Coons directed fallback guard/FP/rawCSR/errno changed");
    check(directedSerial.size()==directedParallel.size(),"Coons directed fallback extent changed");
    for(std::size_t i=0;i<directedSerial.size();++i)sameFlux(directedSerial[i].full,directedParallel[i].full);
}
#endif
void inspect(QuadtreeMesh& mesh,const Publication& source,const Transport& t,bool material=true)
{
    const auto saved=source.plicCells();const auto id=source.identity(),target=t.endpoint().identity();
    const auto clock=source.timeLayer();const auto* storage=source.plicCells().data();
    std::vector<GeometricVofNativeSubface> nq;
    std::vector<GeometricVofInteriorSubface> iq;
    for(const auto f:source.domain().faces(mesh)){
        if(f.negativeCell<0&&f.positiveCell<0)continue;
        const auto e=faceEndpoints(mesh,f);const Point mid{e[0].x/2+e[1].x/2,e[0].y/2+e[1].y/2};
        nq.push_back({f.nativeRow,e[0],e[1]});nq.push_back({f.nativeRow,e[1],e[0]});
        nq.push_back({f.nativeRow,e[0],mid});nq.push_back({f.nativeRow,mid,e[1]});nq.push_back({f.nativeRow,e[0],e[0]});
        if(t.faces()[f.nativeRow].half.limiter!=t.faces()[f.nativeRow].full.limiter)++limiterDifferences;
    }
    for(int ci:source.domain().cellIndexByRow()){
        const auto box=bounds(mesh,ci);const Point center{box.x0/2+box.x1/2,box.y0/2+box.y1/2};
        for(int axis=0;axis<2;++axis){const Point a=axis?Point{box.x0,center.y}:Point{center.x,box.y0};
            const Point b=axis?Point{box.x1,center.y}:Point{center.x,box.y1};
            iq.push_back({ci,a,b});iq.push_back({ci,b,a});iq.push_back({ci,a,center});iq.push_back({ci,center,b});}
    }
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
    parallelInterior(mesh,t,iq);
#endif
    GeometricVofTransportReport nr,nf,ir,ifull;
    const auto oldNative=t.restrictNativeFaces(mesh,nq,&nr);const auto fullNative=t.restrictNativeFacesFull(mesh,nq,&nf);
    const auto oldInterior=t.measureInteriorFaces(mesh,iq,&ir);const auto fullInterior=t.measureInteriorFacesFull(mesh,iq,&ifull);
    check(oldNative.size()==fullNative.size()&&oldInterior.size()==fullInterior.size(),"full-only result shape changed");
    check(nr.localNativeP0Restrictions==nf.localNativeP0Restrictions&&ir.localInteriorP0Extensions==ifull.localInteriorP0Extensions,
        "full-only actual query receipt changed");
    for(std::size_t i=0;i<nq.size();++i){const auto& a=oldNative[i];const auto& b=fullNative[i];
        check(a.nativeRow==b.nativeRow&&a.negativeCell==b.negativeCell&&a.positiveCell==b.positiveCell,"full-only native incidence changed");sameFlux(a.full,b.full);}
    for(std::size_t i=0;i<iq.size();++i){check(oldInterior[i].nativeCell==fullInterior[i].nativeCell,"full-only actual interior owner changed");
        sameFlux(oldInterior[i].full,fullInterior[i].full);}
    // A whole tangent section consumes all four actual cell sides. Its
    // original full limiter includes BOTH time extents on every such side.
    std::vector<long double> minimum(source.plicCells().size(),1.L);
    for(const auto& f:t.faces())for(int ci:{f.negativeCell,f.positiveCell})if(ci>=0){
        minimum[ci]=std::min(minimum[ci],f.half.limiter);minimum[ci]=std::min(minimum[ci],f.full.limiter);}
    for(std::size_t i=0;i<iq.size();i+=4)check(fullInterior[i].full.limiter==minimum[iq[i].nativeCell],
        "full-only lost half limiter contribution to actual full interior metadata");
    sameFailure([&]{t.restrictNativeFaces(mesh,{{-1,{0,0},{0,0}}});},
        [&]{t.restrictNativeFacesFull(mesh,{{-1,{0,0},{0,0}}});});
    sameFailure([&]{t.measureInteriorFaces(mesh,{{-1,{0,0},{1,0}}});},
        [&]{t.measureInteriorFacesFull(mesh,{{-1,{0,0},{1,0}}});});
    auto outside=nq.front();outside.begin.x=std::numeric_limits<long double>::infinity();
    sameFailure([&]{t.restrictNativeFaces(mesh,{outside});},[&]{t.restrictNativeFacesFull(mesh,{outside});});
    auto badInterior=iq.front();badInterior.end=badInterior.begin;
    sameFailure([&]{t.measureInteriorFaces(mesh,{badInterior});},[&]{t.measureInteriorFacesFull(mesh,{badInterior});});
    if(material){
        const auto geometry=HydroMomentumGeometry::build(mesh,vofMaterialMetrics(mesh,source));
        const auto before=HydroPhaseGeometry::build(mesh,geometry.metrics(),source),mid=HydroPhaseGeometry::build(mesh,geometry.metrics(),t.midpoint()),
            after=HydroPhaseGeometry::build(mesh,geometry.metrics(),t.endpoint());
        const auto mt=HydroMaterialTransport::prepareGeometricVof(mesh,geometry,before,mid,after,t);
        mt.validateFor(mesh,mt.sources(),mt.timeLayer());
        nq.clear();iq.clear();std::vector<std::pair<bool,std::size_t>> origins;
        for(const auto& atom:mt.fluxReconstruction().atoms()){
            const Point a{atom.begin.x,atom.begin.y},b{atom.end.x,atom.end.y};
            if(atom.primalVelocity>=0){origins.emplace_back(true,nq.size());nq.push_back({geometry.layout().nativeFace(atom.primalVelocity),a,b});}
            else {origins.emplace_back(false,iq.size());iq.push_back({geometry.layout().nativeCell(atom.ownerPressure),a,b});}}
        const auto oldN=t.restrictNativeFaces(mesh,nq);const auto oldI=t.measureInteriorFaces(mesh,iq);
        check(mt.atomIntegrals().size()==origins.size()&&mt.phaseAtomIntegrals().size()==origins.size(),"full-only material atom shape changed");
        for(std::size_t k=0;k<origins.size();++k){const auto& origin=origins[k];const auto& expected=origin.first?oldN[origin.second].full:oldI[origin.second].full;
            const auto& phases=mt.phaseAtomIntegrals()[k];
            for(const auto pair:{std::pair<long double,long double>{mt.atomIntegrals()[k],expected.mass},{phases.liquid,expected.liquid},{phases.gas,expected.gas}})
                check(pair.first==pair.second&&std::signbit(pair.first)==std::signbit(pair.second),"material did not consume original full atom bits");}
        materialAtoms+=static_cast<unsigned>(origins.size());
    }
    nativeQueries+=static_cast<unsigned>(oldNative.size());interiorQueries+=static_cast<unsigned>(oldInterior.size());
    check(source.identity()==id&&t.endpoint().identity()==target&&sameClock(source.timeLayer(),clock)&&source.plicCells().data()==storage,
        "full-only query replaced source/target owner/fullclock/storage");
    for(std::size_t i=0;i<saved.size();++i)check(source.plicCells()[i].inventory==saved[i].inventory
        &&source.plicCells()[i].measured==saved[i].measured&&source.plicCells()[i].resolved==saved[i].resolved,
        "full-only material query modified source stocks/measured/unknown label");++scenes;
}
void standard(bool rz,bool micro)
{
    stage="full-only original RZ/XY 2:1 "+std::to_string(rz)+","+std::to_string(micro);
    auto mesh=F15BoundaryGrossContract::makeMesh(rz);const auto native=domain(*mesh,rz);
    check(!native.uniformLeaves(),"full-only fixture lost actual 2:1 faces");
    const unsigned minor=rz?0:1;
    const auto seed=Publication::prepareLocalPlic(*mesh,native,Clock{},micro?std::vector<double>(mesh->cellIndexNum(),minor?1.:0.):l1PlaneAlpha(*mesh,rz,0));
    auto source=seed;
    if(micro){std::vector<std::array<long double,2>> stock;for(std::size_t ci=0;ci<seed.plicCells().size();++ci){
            auto inventory=seed.plicCells()[ci].inventory;inventory[1-minor]=seed.cells()[ci].capacity;inventory[minor]=std::ldexp(inventory[1-minor],-80);stock.push_back(inventory);}
        source=Publication::prepareLocalPlicInventory(*mesh,native,Clock{},stock,{});}
    const auto a=micro?F15BoundaryGrossContract::normals(*mesh,source,.125):vofMaterialNormals(*mesh,source,0);
    const auto b=micro?F15BoundaryGrossContract::normals(*mesh,source,-.375):a;
    const double dt=micro?1./32:static_cast<double>(.1L/vofMaterialOutflow(*mesh,source,a,b));const Clock end{1,0,dt,Position::End};
    if(micro){const auto boundary=Slab::Boundary::prepare(*mesh,source,F15BoundaryGrossContract::rulesFor(*mesh,source,1-minor));
        const auto slab=Slab::prepareBoundary(*mesh,source,boundary,end,a,b);const auto t=Transport::prepareLocalPlicBoundarySweep(*mesh,source,slab,100,1);
        check(t.report().localUnresolvedDonors>0&&t.report().boundaryFull.enabled,"full-only missing actual open unknown-micro consumption");
        bool both=false,positiveMinor=false;for(const auto& f:t.faces()){
            both|=f.full.positiveTotal>0&&f.full.negativeTotal>0;
            positiveMinor|=minor?(f.full.positiveGas>0||f.full.negativeGas>0):(f.full.positiveLiquid>0||f.full.negativeLiquid>0);}
        check(both&&positiveMinor,"full-only open fixture lacks two gross signs or positive minor");inspect(*mesh,source,t);
    }else {const auto slab=Slab::prepare(*mesh,source,end,a,b);const auto t=Transport::prepareLocalPlicSweep(*mesh,source,slab,100,1);inspect(*mesh,source,t);}
}
void selected(bool rz)
{
    stage="full-only actual selected nozzle "+std::to_string(rz);
    auto mesh=SelectedFluidTransportContract::nozzleMesh(true);const auto fluid=SelectedFluidTransportContract::selected(*mesh,rz);
    const auto source=GeometricVofInitialization::prepare(*mesh,fluid,[](Vector2D p){return std::max(p.x-1,p.y-4.125);});
    std::vector<double> velocity(fluid.nativeFaceCount());const auto psi=[](Point p){
        if(p.x>=1||p.y<=3||p.y>=5)return 0.L;
        return p.x*p.x*(1-p.x)*(1-p.x)*(p.y-3)*(p.y-3)*(5-p.y)*(5-p.y)/256;};
    bool hanging=false;for(const auto f:fluid.faces(*mesh)){
        if(f.negativeCell<0||f.positiveCell<0)continue;
        const auto e=faceEndpoints(*mesh,f);const auto q=(f.horizontal?-1.L:1.L)*(psi(e[1])-psi(e[0]));
        if(f.area>0)velocity[f.nativeRow]=static_cast<double>(q/f.area);else check(q==0,"full-only selected axis leak");
        hanging|=mesh->cells(mesh->cellIndexToLabel[f.negativeCell]).level()!=mesh->cells(mesh->cellIndexToLabel[f.positiveCell]).level();}
    check(hanging&&fluid.cellIndexByRow().size()<static_cast<std::size_t>(mesh->cellIndexNum()),"full-only selected fixture lost hanging/excluded leaves");
    const auto slab=Slab::prepare(*mesh,source,{1,0,.01,Position::End},velocity,velocity);
    const auto t=Transport::prepareLocalPlicSweep(*mesh,source,slab,100,1);inspect(*mesh,source,t);
}
void differingLimiter()
{
    // Reuse the original D0 32x32 sole-outgoing resolved-minor fixture. The
    // same full sweep is retained; reversing the first endpoint makes the
    // half mean nearly zero instead of consuming that tiny corner phase.
    stage="full-only actual half/full limiter distinction";
    DictionaryMesh d;d.Nx=d.Ny=32;d.maxCellNum=4096;d.maxInitialIterations=0;QuadtreeMesh mesh(d);
    const auto source=GeometricVofInitialization::prepare(mesh,false,[](Vector2D p){return std::hypot(p.x-.487,p.y-.513)-.18;});
    const auto taper=[](long double x){if(x<=0||x>=1)return 0.L;if(x>=.2L&&x<=.8L)return 1.L;
        const auto z=x<.2L?x/.2L:(1-x)/.2L;return z*z*(3-2*z);};
    const auto psi=[&](long double x,long double y){return .04L*y*taper(x)*taper(y);};
    std::vector<double> a(source.domain().nativeFaceCount()),b(a.size());
    for(const auto f:source.domain().faces(mesh)){const long double x=f.center.x,y=f.center.y,h=f.area/2;
        const auto v=static_cast<double>(f.horizontal?-(psi(x+h,y)-psi(x-h,y))/f.area:(psi(x,y+h)-psi(x,y-h))/f.area);
        a[f.nativeRow]=-v;b[f.nativeRow]=3*v;}
    const auto slab=Slab::prepare(mesh,source,{1,0,.025,Position::End},a,b);
    const auto t=Transport::prepareLocalPlicSweep(mesh,source,slab,100,1);const auto before=limiterDifferences;
    inspect(mesh,source,t,false);
    check(t.report().limitedFaces>0&&limiterDifferences>before,"full-only fixture did not distinguish actual half/full retained limiter");
}
void run()
{
    for(bool rz:{false,true}){standard(rz,false);standard(rz,true);selected(rz);}differingLimiter();
    check(scenes==7&&nativeQueries>0&&interiorQueries>0&&materialAtoms>0&&negativeControls==28&&limiterDifferences>0,
        "full-only finite required coverage incomplete");
    std::cout<<"LOCAL_PLIC_FULL_MATERIAL_PASS scenes="<<scenes<<" native_queries="<<nativeQueries<<" interior_queries="<<interiorQueries
        <<" material_atoms="<<materialAtoms<<" common_guards="<<negativeControls<<" differing_half_full_limiters="<<limiterDifferences
        <<" original_full_bits=exact source_stock_owner_clock=unchanged no_speed_claim\n";
}
}

namespace StockHeightTransportContract
{
using Interval=std::array<long double,2>;
std::optional<Interval> trace(const Publication& p,const QuadtreeMesh& mesh,int ci,bool horizontal,long double fixed,unsigned phase)
{
    const auto& q=p.plicCells()[ci];const auto b=bounds(mesh,ci);const long double cx=b.x0+(b.x1-b.x0)/2,cy=b.y0+(b.y1-b.y0)/2;
    const long double hx=(b.x1-b.x0)/2,hy=(b.y1-b.y0)/2;
    check(horizontal?(fixed==b.y0||fixed==b.y1):(fixed==b.x0||fixed==b.x1),"HEIGHT_TRANSPORT native face does not lie on the actual whole cell side");
    const auto side=horizontal?(fixed==b.y0?-hy:hy):(fixed==b.x0?-hx:hx);
    const Point a=horizontal?Point{-hx,side}:Point{side,-hy},z=horizontal?Point{hx,side}:Point{side,hy};
    const long double lower=horizontal?b.x0:b.y0,upper=horizontal?b.x1:b.y1;
    if(q.measured[phase]==0)return {};if(q.measured[1-phase]==0)return Interval{lower,upper};
    if(!q.resolved)return {}; // Unknown geometry is reported, not given a fabricated trace.
    const auto length=std::hypot(static_cast<long double>(q.normal.x),static_cast<long double>(q.normal.y));
    const long double sign=phase?-1.L:1.L,nx=sign*q.normal.x/length,ny=sign*q.normal.y/length,offset=sign*q.offset;
    const auto level=[&](Point r){return nx*r.x+ny*r.y-offset;};const auto fa=level(a),fz=level(z);
    if(fa<=0&&fz<=0)return Interval{lower,upper};if(fa>0&&fz>0)return {};
    // Use the actual public finite PLIC segment, rather than an ideal plane
    // crossing, to decide which part of this whole native cell side is wet.
    const auto segment=LAMP::Geometry::PlanarCellCut::interfaceSegment(q.rectangle,q.normal,q.offset);
    check(segment.has_value(),"HEIGHT_TRANSPORT mixed side lost its native finite PLIC segment");
    std::optional<long double> crossing;
    for(auto site:*segment)if(horizontal?site.y==a.y:site.x==a.x)crossing=horizontal?site.x+cx:site.y+cy;
    check(crossing.has_value(),"HEIGHT_TRANSPORT actual finite segment misses a signed crossing side");
    return fa<=0?Interval{lower,*crossing}:Interval{*crossing,upper};
}
unsigned geometryBodies(const Publication& p,const QuadtreeMesh& mesh,unsigned phase,unsigned& mixedWitnesses,unsigned& splitWitnesses,unsigned& singletons,unsigned& unknown)
{
    const auto count=p.plicCells().size();std::vector<std::vector<int>> neighbours(count);
    for(int row=0;row<p.domain().nativeFaceCount();++row){const auto f=p.domain().face(mesh,row);
        if(f.negativeCell<0||f.positiveCell<0)continue;
        const auto& edge=f.horizontal?mesh.hEdges(f.edgeLabel):mesh.vEdges(f.edgeLabel);
        const auto a=mesh.verts(edge.edgeToVert[0]).coord(),z=mesh.verts(edge.edgeToVert[1]).coord();
        const long double fixed=f.horizontal?a.y:a.x,lo=f.horizontal?std::min(a.x,z.x):std::min(a.y,z.y),hi=f.horizontal?std::max(a.x,z.x):std::max(a.y,z.y);
        const auto u=trace(p,mesh,f.negativeCell,f.horizontal,fixed,phase),v=trace(p,mesh,f.positiveCell,f.horizontal,fixed,phase);
        if(!u||!v||std::min({hi,(*u)[1],(*v)[1]})<=std::max({lo,(*u)[0],(*v)[0]}))continue;
        neighbours[f.negativeCell].push_back(f.positiveCell);neighbours[f.positiveCell].push_back(f.negativeCell);
        if(mesh.cells(mesh.cellIndexToLabel[f.negativeCell]).level()!=mesh.cells(mesh.cellIndexToLabel[f.positiveCell]).level())++splitWitnesses;
    }
    std::vector<bool> seen(count);unsigned bodies=0,mixed=0;singletons=unknown=0;
    for(std::size_t seed=0;seed<count;++seed){const auto& q=p.plicCells()[seed];
        if(q.measured[phase]>0&&q.measured[1-phase]>0&&!q.resolved){++unknown;continue;}
        if(seen[seed]||q.measured[phase]==0)continue;
        ++bodies;std::vector<int> queue{static_cast<int>(seed)};seen[seed]=true;
        for(std::size_t head=0;head<queue.size();++head){const int ci=queue[head];if(p.plicCells()[ci].measured[1-phase]>0)++mixed;
            for(int next:neighbours[ci])if(!seen[next]){seen[next]=true;queue.push_back(next);}}
        if(queue.size()==1)++singletons;
    }
    mixedWitnesses+=mixed;return bodies;
}
std::array<unsigned,2> stockBodies(const Publication& p,const std::vector<std::vector<int>>& neighbours,unsigned phase)
{
    std::array<unsigned,2> counts{};
    for(unsigned kind=0;kind<2;++kind){std::vector<bool> seen(p.plicCells().size());
        const auto present=[&](int ci){const auto& q=p.plicCells()[ci].inventory;
            return q[phase]>0&&(kind==0||q[phase]>=q[1-phase]);};
        for(std::size_t seed=0;seed<seen.size();++seed){if(seen[seed]||!present(static_cast<int>(seed)))continue;
            ++counts[kind];std::vector<int> queue{static_cast<int>(seed)};seen[seed]=true;
            for(std::size_t head=0;head<queue.size();++head)for(int next:neighbours[queue[head]])if(!seen[next]&&present(next)){
                seen[next]=true;queue.push_back(next);}}
    }return counts;
}
void run()
{
    unsigned scenes=0,states=0,mixedWitnesses=0,splitWitnesses=0,forward=0,reverse=0;long double peakVolume=0;
    for(bool adaptive:{false,true}){stage="stock-height real split transport "+std::to_string(adaptive);
        std::unique_ptr<QuadtreeMesh> mesh;
        // Dyadic leaves make equal-size qualification a literal mesh fact,
        // including each side of the static 2:1 interface.
        if(adaptive)mesh=vofMaterialMesh(16);else{DictionaryMesh d;d.Nx=d.Ny=32;d.maxCellNum=4096;d.maxInitialIterations=0;mesh=std::make_unique<QuadtreeMesh>(d);}
        std::vector<std::array<long double,2>> stocks(mesh->cellIndexNum());std::array<long double,2> initial{};
        for(int ci=0;ci<mesh->cellIndexNum();++ci){const auto b=bounds(*mesh,ci);const auto c=capacity(b,false);
            const auto liquid=l1EllipseArea(b.x0,b.x1,b.y0,b.y1);check(liquid>=0&&liquid<=c,"HEIGHT_TRANSPORT independent ellipse outside capacity");
            stocks[ci]={liquid,c-liquid};for(unsigned phase=0;phase<2;++phase)initial[phase]+=stocks[ci][phase];}
        auto accepted=Publication::prepareLocalPlicInventory(*mesh,domain(*mesh,false),{},stocks,{});
        std::vector<std::vector<int>> stockNeighbours(stocks.size());
        for(const auto f:accepted.domain().faces(*mesh))if(f.area>0&&f.negativeCell>=0&&f.positiveCell>=0){
            stockNeighbours[f.negativeCell].push_back(f.positiveCell);stockNeighbours[f.positiveCell].push_back(f.negativeCell);}
        std::array<std::array<unsigned,2>,2> initialStockBodies{};std::array<unsigned,2> initialGeometry{},peakGeometry{},initialSingletons{},peakSingletons{};
        unsigned unknownGeometry=0;bool firstState=true;
        const auto velocity=vofMaterialNormals(*mesh,accepted,0);const double dt=static_cast<double>(.1L/vofMaterialOutflow(*mesh,accepted,velocity,velocity));
        const auto inspect=[&](const Publication& p){std::array<long double,2> total{};
            for(std::size_t ci=0;ci<p.plicCells().size();++ci){const auto& q=p.plicCells()[ci];const auto c=p.cells()[ci].capacity;
                for(unsigned phase=0;phase<2;++phase){check(q.inventory[phase]>=0&&q.inventory[phase]<=c,"HEIGHT_TRANSPORT invalid independent phase stock");
                    total[phase]+=q.inventory[phase];near(q.measured[phase],q.inventory[phase],c,"HEIGHT_TRANSPORT retained stock inverse volume",1e-12L);}}
            for(unsigned phase=0;phase<2;++phase){peakVolume=std::max(peakVolume,std::abs(total[phase]-initial[phase]));
                near(total[phase],initial[phase],1,"HEIGHT_TRANSPORT real split independent phase conservation",1e-12L);
                const auto current=stockBodies(p,stockNeighbours,phase);
                if(firstState){initialStockBodies[phase]=current;check(current[0]==1&&current[1]==1,"HEIGHT_TRANSPORT analytic initial stock-support/majority body is not connected");}
                std::ostringstream reason;reason<<"HEIGHT_TRANSPORT new inventory-supported phase component phase="<<phase<<" step="<<p.timeLayer().step
                    <<" position="<<static_cast<int>(p.timeLayer().position)<<" strict_positive="<<current[0]<<" majority="<<current[1]
                    <<" initial="<<initialStockBodies[phase][0]<<'/'<<initialStockBodies[phase][1];
                check(current==initialStockBodies[phase],reason.str());
                unsigned singles=0,unknown=0;const auto geometry=geometryBodies(p,*mesh,phase,mixedWitnesses,splitWitnesses,singles,unknown);
                if(firstState){initialGeometry[phase]=geometry;initialSingletons[phase]=singles;}
                peakGeometry[phase]=std::max(peakGeometry[phase],geometry);peakSingletons[phase]=std::max(peakSingletons[phase],singles);unknownGeometry+=unknown;
            }firstState=false;++states;};
        inspect(accepted);
        for(int step=1;step<=8;++step){auto u=velocity;if(step>4){for(auto& v:u)v=-v;++reverse;}else ++forward;
            const auto old=accepted.plicCells();const auto identity=accepted.identity();const auto oldClock=accepted.timeLayer();const auto* storage=accepted.plicCells().data();
            const Clock next{oldClock.step+1,oldClock.time(),oldClock.time()+dt,Position::End};
            const auto slab=Slab::prepare(*mesh,accepted,next,u,u);const auto candidate=Transport::prepareLocalPlicSweep(*mesh,accepted,slab,100,1);
            candidate.validateFor(*mesh,accepted,slab,100,1);inspect(candidate.midpoint());inspect(candidate.endpoint());
            check(accepted.identity()==identity&&accepted.plicCells().data()==storage&&sameClock(accepted.timeLayer(),oldClock),"HEIGHT_TRANSPORT preparation mutated accepted owner/clock/storage");
            for(std::size_t ci=0;ci<old.size();++ci){const auto& q=accepted.plicCells()[ci];check(q.inventory==old[ci].inventory&&q.measured==old[ci].measured
                &&q.normal.x==old[ci].normal.x&&q.normal.y==old[ci].normal.y&&q.offset==old[ci].offset&&q.resolved==old[ci].resolved,
                "HEIGHT_TRANSPORT candidate changed accepted stock or actual plane");}
            accepted=candidate.endpoint();
        }
        // Native PLIC independently fits each inventory. It promises neither
        // C0 continuity nor one connected finite polygon union: the analytic
        // initial ellipse itself has four tiny separated gas cut tips. Keep
        // that diagnostic separate from genuinely new inventory support.
        std::cout<<"STOCK_HEIGHT_TRANSPORT_TOPOLOGY,adaptive,"<<adaptive<<",initial_stock_positive_bodies,"<<initialStockBodies[0][0]<<'/'<<initialStockBodies[1][0]
            <<",initial_stock_majority_bodies,"<<initialStockBodies[0][1]<<'/'<<initialStockBodies[1][1]<<",initial_finite_geometry_bodies,"<<initialGeometry[0]<<'/'<<initialGeometry[1]
            <<",max_finite_geometry_bodies,"<<peakGeometry[0]<<'/'<<peakGeometry[1]<<",initial_finite_singletons,"<<initialSingletons[0]<<'/'<<initialSingletons[1]
            <<",max_finite_singletons,"<<peakSingletons[0]<<'/'<<peakSingletons[1]<<",unknown_mixed_geometry,"<<unknownGeometry<<"\n";++scenes;
    }
    check(scenes==2&&states==34&&forward==8&&reverse==8&&splitWitnesses>0,"HEIGHT_TRANSPORT uniform/AMR forward/reverse coverage missing");
    std::cout<<std::setprecision(17)<<"STOCK_HEIGHT_TRANSPORT_PASS,scenes,"<<scenes<<",half_end_states,"<<states<<",positive_steps,"<<forward<<",reverse_steps,"<<reverse
        <<",actual_mixed_component_witnesses,"<<mixedWitnesses<<",actual_2to1_trace_witnesses,"<<splitWitnesses<<",max_phase_volume_error,"<<peakVolume
        <<",stock_support_and_majority_components_unchanged,1,stocks_planes_accepted_unchanged,NO_C0_OR_REVERSIBILITY_ACCURACY_CLAIM\n";
}
}

namespace NativeInterfaceTransportContract
{
void opposedFilmForce()
{
    unsigned scenes=0,closedCells=0;
    for(unsigned liquidSide:{0u,1u}){
        stage="native interface support one-cell-thick sheet force phase "+std::to_string(liquidSide);
        DictionaryMesh dict;dict.Nx=dict.Ny=16;dict.maxCellNum=2048;dict.maxInitialIterations=0;QuadtreeMesh mesh(dict);
        constexpr long double h=1.L/16;std::map<std::pair<int,int>,int> grid;std::vector<std::array<long double,2>> stocks(mesh.cellIndexNum());
        // Exact rectangle [4h,12h] x [7.25h,8.25h], entirely inside the real
        // fluid domain. Neither phase lacks a local majority witness here.
        for(int ci=0;ci<mesh.cellIndexNum();++ci){const auto b=bounds(mesh,ci);const int x=static_cast<int>(std::llround(b.x0/h)),y=static_cast<int>(std::llround(b.y0/h));
            grid[{x,y}]=ci;const auto c=capacity(b,false);const long double fraction=x>=4&&x<12?(y==7?.75L:y==8?.25L:0.L):0.L;
            stocks[ci][liquidSide]=fraction*c;stocks[ci][1-liquidSide]=(1-fraction)*c;}
        GeometricVofConfig config;config.localResolution=GeometricVofLocalResolution::NativeInterfaceSupport;
        const Clock clock{39,.2,.2002,Position::Midpoint};const auto publication=Publication::prepareLocalPlicInventory(mesh,domain(mesh,false),clock,stocks,config);
        const auto original=publication.plicCells();const auto mask=publication.localPlicSubgridMask(mesh);const int lower=grid.at({7,7}),upper=grid.at({7,8});
        check(stocks[lower][liquidSide]>stocks[lower][1-liquidSide]&&stocks[upper][1-liquidSide]>stocks[upper][liquidSide],
            "NATIVE_INTERFACE_SUPPORT sheet force fixture lacks its independent witnesses for both phases");
        const auto& a=original[lower];const auto& b=original[upper];
        check(a.resolved&&b.resolved&&a.measured[0]>0&&a.measured[1]>0&&b.measured[0]>0&&b.measured[1]>0
            &&static_cast<long double>(a.normal.x)*b.normal.x+static_cast<long double>(a.normal.y)*b.normal.y<=0,
            "NATIVE_INTERFACE_SUPPORT sheet force lacks independently reconstructed opposed actual planes");
        unsigned mixed=0;for(std::size_t ci=0;ci<stocks.size();++ci)if(stocks[ci][0]>0&&stocks[ci][1]>0){
            check(mask[ci],"NATIVE_INTERFACE_SUPPORT complete interior thin-sheet member did not qualify");++mixed;}
        check(mixed==16,"NATIVE_INTERFACE_SUPPORT exact rectangular-sheet raw-stock coverage changed");
        const auto metrics=HydroPhysicalDivergence::build(mesh,HydroIntegralDivergence::build(mesh,HydroDofLayout::build(mesh,publication.domain())));
        const auto phase=HydroPhaseGeometry::build(mesh,metrics,publication);const auto pressure=HydroPressureLoad::build(metrics.divergence());
        HydroCsfCurvatureConfig curvature;curvature.forceMethod=HydroCsfForceMethod::HeightFunction;
        const auto load=HydroBalancedCsfLoad::localPlicHeightFunctionCapillary(mesh,pressure,phase,.1,curvature);
        load.validateVofFor(mesh,pressure,publication.identity(),clock);
        for(std::size_t ci=0;ci<stocks.size();++ci){const auto& q=publication.plicCells()[ci];const auto& saved=original[ci];
            check(q.inventory==stocks[ci]&&q.inventory==saved.inventory&&q.measured==saved.measured&&q.resolved==saved.resolved
                &&q.normal.x==saved.normal.x&&q.normal.y==saved.normal.y&&q.offset==saved.offset,
                "NATIVE_INTERFACE_SUPPORT thin-sheet force closure changed stock or actual PLIC geometry");
            if(!mask[ci])continue;const int row=phase.layout().pressureDof(static_cast<int>(ci));
            check(load.heightFunctionCellCurvatureMethods()[row]==HydroCsfCellCurvatureMethod::Unavailable
                &&std::isnan(load.heightFunctionCellCurvatures()[row])&&load.heightFunctionCellEffectiveCurvatures()[row]==0
                &&!load.heightFunctionCellJets()[row],"NATIVE_INTERFACE_SUPPORT thin-sheet closure reported a physical/HF/fit/isolated curvature");++closedCells;}
        check(load.curvatureReport().heightFunctionCells==0&&load.curvatureReport().heightFunctionFitCells==0&&load.curvatureReport().unresolvedIsolatedFragmentCells==0,
            "NATIVE_INTERFACE_SUPPORT thin sheet was counted as a successful or isolated geometric estimator");
        check(std::all_of(load.momentumRhsLoad().begin(),load.momentumRhsLoad().end(),[](double q){return q==0;})
            &&std::any_of(load.supportReports().begin(),load.supportReports().end(),[](const auto& q){return q.underresolvedFragmentClosure;}),
            "NATIVE_INTERFACE_SUPPORT thin-sheet load lost its explicit zero effective/unavailable physical curvature closure");++scenes;
    }
    check(scenes==2&&closedCells==32,"NATIVE_INTERFACE_SUPPORT thin-sheet phase/force coverage incomplete");
    std::cout<<"NATIVE_INTERFACE_SUPPORT_SHEET_FORCE_PASS,phase_scenes,"<<scenes<<",unavailable_physical_curvature_cells,"<<closedCells
        <<",both_phase_majorities_present,1,opposed_actual_PLIC_normals,1,raw_stocks_and_planes_unchanged,1\n";
}
void run()
{
    unsigned scenes=0,stepsTaken=0,geometricFaces=0,strictZeroFaces=0,triangleFaces=0,sharedLedgers=0,closures=0;
    long double maxConservation=0,maximumRetained=0,peakTriangleRelative=0;
    for(unsigned minor:{0u,1u}){
        stage="native interface support resolved geometric transport and unavailable physical curvature phase "+std::to_string(minor);
        DictionaryMesh dict;dict.Nx=dict.Ny=16;dict.maxCellNum=2048;dict.maxInitialIterations=0;QuadtreeMesh mesh(dict);
        constexpr long double h=1.L/16;std::map<std::pair<int,int>,int> grid;const auto native=domain(mesh,false);
        std::vector<std::array<long double,2>> stocks(mesh.cellIndexNum());std::array<long double,2> initial{};
        for(int ci=0;ci<mesh.cellIndexNum();++ci){const auto b=bounds(mesh,ci);grid[{static_cast<int>(std::llround(b.x0/h)),static_cast<int>(std::llround(b.y0/h))}]=ci;
            stocks[ci][1-minor]=capacity(b,false);}
        const int source=grid.at({7,7}),east=grid.at({8,7}),north=grid.at({7,8});
        const auto put=[&](int ci,long double fraction){const auto c=capacity(bounds(mesh,ci),false);stocks[ci][minor]=fraction*c;stocks[ci][1-minor]=(1-fraction)*c;};
        put(source,1.L/64);put(east,1.L/1024);put(north,1.L/512);
        for(const auto& q:stocks)for(unsigned phase=0;phase<2;++phase)initial[phase]+=q[phase];
        GeometricVofConfig resolution;resolution.localResolution=GeometricVofLocalResolution::NativeInterfaceSupport;
        auto accepted=Publication::prepareLocalPlicInventory(mesh,native,{},stocks,resolution);const auto initialSource=stocks[source][minor],initialEast=stocks[east][minor];
        const auto modelMask=accepted.localPlicSubgridMask(mesh);for(int ci:{source,east,north})check(modelMask[ci]&&accepted.plicCells()[ci].resolved,
            "NATIVE_INTERFACE_SUPPORT geometric source/tails require actual resolved planes and strict owner majority-deficit classification");
        const auto original=accepted.plicCells();const auto oldId=accepted.identity();const auto oldClock=accepted.timeLayer();
        const auto& triangle=original[source];const long double orientation=minor==0?1.L:-1.L;
        const auto norm=std::hypot(static_cast<long double>(triangle.normal.x),static_cast<long double>(triangle.normal.y));
        const auto nx=orientation*triangle.normal.x/norm,ny=orientation*triangle.normal.y/norm,offset=orientation*triangle.offset;
        const auto cx=(static_cast<long double>(triangle.rectangle.lower.x)+triangle.rectangle.upper.x)/2,
            cy=(static_cast<long double>(triangle.rectangle.lower.y)+triangle.rectangle.upper.y)/2;
        const auto hx=(static_cast<long double>(triangle.rectangle.upper.x)-triangle.rectangle.lower.x)/2,
            hy=(static_cast<long double>(triangle.rectangle.upper.y)-triangle.rectangle.lower.y)/2;
        check(nx<0&&ny<0&&nx*hx+ny*hy-offset<0&&-nx*hx+ny*hy-offset>0&&nx*hx-ny*hy-offset>0&&-nx*hx-ny*hy-offset>0,
            "NATIVE_INTERFACE_SUPPORT geometric fixture must independently contain only its strict NE minor-phase corner");
        // Solve the actual plane's top and right edge intercepts. This is an
        // independent triangle, not a polygon returned by the Core cutter.
        const auto legX=triangle.rectangle.upper.x-(cx+(offset-ny*hy)/nx),legY=triangle.rectangle.upper.y-(cy+(offset-nx*hx)/ny);
        const auto triangleArea=legX*legY/2;
        check(legX>0&&legX<h&&legY>0&&legY<h,"NATIVE_INTERFACE_SUPPORT actual corner triangle escaped its native rectangle");
        near(triangle.measured[minor],triangleArea,triangleArea,"NATIVE_INTERFACE_SUPPORT independent actual triangle area",1e-12L);
        const auto metrics=HydroPhysicalDivergence::build(mesh,HydroIntegralDivergence::build(mesh,
            HydroDofLayout::build(mesh,accepted.domain())));const auto phase=HydroPhaseGeometry::build(mesh,metrics,accepted);
        const auto pressure=HydroPressureLoad::build(metrics.divergence());HydroCsfCurvatureConfig curvature;curvature.forceMethod=HydroCsfForceMethod::HeightFunction;
        const auto load=HydroBalancedCsfLoad::localPlicHeightFunctionCapillary(mesh,pressure,phase,.1,curvature);
        load.validateVofFor(mesh,pressure,oldId,oldClock);
        for(int ci:{source,east,north}){const int row=phase.layout().pressureDof(ci);
            check(row>=0&&load.heightFunctionCellCurvatureMethods()[row]==HydroCsfCellCurvatureMethod::Unavailable
                &&std::isnan(load.heightFunctionCellCurvatures()[row])&&load.heightFunctionCellEffectiveCurvatures()[row]==0
                &&!load.heightFunctionCellJets()[row],"NATIVE_INTERFACE_SUPPORT subgrid force invented a physical curvature/HF/fit/isolated estimate");++closures;}
        check(load.curvatureReport().heightFunctionCells==0&&load.curvatureReport().heightFunctionFitCells==0
            &&load.curvatureReport().unresolvedIsolatedFragmentCells==0,"NATIVE_INTERFACE_SUPPORT declared closure entered a successful or isolated geometric estimator counter");
        check(std::all_of(load.momentumRhsLoad().begin(),load.momentumRhsLoad().end(),[](double q){return q==0;}),"NATIVE_INTERFACE_SUPPORT subgrid-only source retained a capillary force");
        check(std::any_of(load.supportReports().begin(),load.supportReports().end(),[](const auto& q){return q.underresolvedFragmentClosure;}),
            "NATIVE_INTERFACE_SUPPORT face report lost the explicit unavailable-geometry force closure");
        check(accepted.identity()==oldId&&sameClock(accepted.timeLayer(),oldClock),"NATIVE_INTERFACE_SUPPORT force preparation changed source owner/full clock");
        for(std::size_t ci=0;ci<stocks.size();++ci)check(accepted.plicCells()[ci].inventory==stocks[ci]
            &&accepted.plicCells()[ci].normal.x==original[ci].normal.x&&accepted.plicCells()[ci].normal.y==original[ci].normal.y
            &&accepted.plicCells()[ci].offset==original[ci].offset&&accepted.plicCells()[ci].measured==original[ci].measured
            &&accepted.plicCells()[ci].resolved==original[ci].resolved,"NATIVE_INTERFACE_SUPPORT force closure changed source inventory/PLIC geometry");
        // Keep the historical carrier and four-step duration. The archived
        // forced-P0 ratio/release acceptance belongs to that old approximation;
        // the new contract instead checks the real resolved swept geometry.
        // This closed solenoidal carrier is exactly u=(.5,0) on the patch.
        const auto psi=[](Vector2D p){const auto r=std::hypot(static_cast<long double>(p.x)-.5L,static_cast<long double>(p.y)-.5L);
            const auto taper=r<=.25L?1.L:r>=.5L?0.L:(.5L-r)/.25L;return .5L*(p.y-.5L)*taper;};
        std::vector<double> velocity(native.nativeFaceCount());for(const auto f:native.faces(mesh)){
            const auto& edge=f.horizontal?mesh.hEdges(f.edgeLabel):mesh.vEdges(f.edgeLabel);auto a=mesh.verts(edge.edgeToVert[0]).coord(),b=mesh.verts(edge.edgeToVert[1]).coord();
            if(f.horizontal?a.x>b.x:a.y>b.y)std::swap(a,b);const auto flux=(f.horizontal?-1.L:1.L)*(psi(b)-psi(a));
            if(f.outwardSign!=0)check(flux==0,"NATIVE_INTERFACE_SUPPORT geometric carrier leaked through a wall");velocity[f.nativeRow]=static_cast<double>(flux/f.area);}
        GeometricVofTransportConfig config;config.geometry=resolution;constexpr double dt=1./128;
        const auto inspect=[&](const Publication& p){std::array<long double,2> total{};
            for(std::size_t ci=0;ci<p.plicCells().size();++ci){const auto& q=p.plicCells()[ci];const auto c=p.cells()[ci].capacity;
                using Exact=LAMP::Geometry::Detail::FiniteVolume::Dyadic;Exact closure(q.inventory[0]);closure.add(q.inventory[1]);closure.add(c,-1);const auto& receipt=p.capacityDefects()[ci];
                check(closure.compare(receipt)==0&&receipt.absolute().compare(Exact::product({1e-12L,c}))<=0,
                    "NATIVE_INTERFACE_SUPPORT geometric exact capacity receipt failed original budget");
                Exact bound(c);bound.add(receipt.absolute());for(unsigned phase=0;phase<2;++phase){
                    check(std::isfinite(q.inventory[phase])&&q.inventory[phase]>=0&&Exact(q.inventory[phase]).compare(bound)<=0,
                        "NATIVE_INTERFACE_SUPPORT geometric stock escaped original certified typed capacity bound");total[phase]+=q.inventory[phase];}}
            for(unsigned phase=0;phase<2;++phase){maxConservation=std::max(maxConservation,std::abs(total[phase]-initial[phase]));
                near(total[phase],initial[phase],1,"NATIVE_INTERFACE_SUPPORT geometric conserved both original phase stocks",1e-12L);}
            near(100*total[0]+total[1],100*initial[0]+initial[1],100*initial[0]+initial[1],"NATIVE_INTERFACE_SUPPORT geometric original two-density mass conservation",1e-12L);};
        const auto unchanged=[&](const Publication& p,const std::vector<GeometricVofPlicCell>& saved,GeometricVofIdentity id,Clock clock){
            check(p.identity()==id&&sameClock(p.timeLayer(),clock),"NATIVE_INTERFACE_SUPPORT geometric candidate changed accepted owner/full clock");
            for(std::size_t ci=0;ci<saved.size();++ci)check(p.plicCells()[ci].inventory==saved[ci].inventory&&p.plicCells()[ci].measured==saved[ci].measured
                &&p.plicCells()[ci].resolved==saved[ci].resolved&&p.plicCells()[ci].normal.x==saved[ci].normal.x&&p.plicCells()[ci].normal.y==saved[ci].normal.y
                &&p.plicCells()[ci].offset==saved[ci].offset,"NATIVE_INTERFACE_SUPPORT geometric candidate changed accepted stocks/planes/measurement/resolution");};
        // Independent probes share the very same immutable public source.
        // With u=.5, duration=legX gives d=legX/2 at full and legX/4 at half.
        // Reversing u consumes a strictly disjoint west strip of that triangle.
        for(int direction:{1,-1}){auto probeVelocity=velocity;for(auto& value:probeVelocity)value*=direction;
            const double duration=static_cast<double>(legX);check(duration>0,"NATIVE_INTERFACE_SUPPORT geometric probe duration vanished");
            const auto slab=Slab::prepare(mesh,accepted,{1,oldClock.time(),oldClock.time()+duration,Position::End},probeVelocity,probeVelocity);
            const auto moved=Transport::prepareLocalPlicSweep(mesh,accepted,slab,100,1,config);moved.validateFor(mesh,accepted,slab,100,1);
            inspect(moved.midpoint());inspect(moved.endpoint());unchanged(accepted,original,oldId,oldClock);
            const int receiver=direction>0?east:grid.at({6,7});unsigned named=0;
            for(const auto& face:moved.faces()){
                if(direction>0?!(face.negativeCell==source&&face.positiveCell==receiver):!(face.positiveCell==source&&face.negativeCell==receiver))continue;
                ++named;for(bool half:{true,false}){const auto& flux=half?face.half:face.full;const auto total=direction>0?flux.positiveTotal:flux.negativeTotal;
                    const auto actualMinor=minor?(direction>0?flux.positiveGas:flux.negativeGas):(direction>0?flux.positiveLiquid:flux.negativeLiquid);
                    const auto signedMinor=minor?flux.gas:flux.liquid;
                    const auto depth=total/h;auto strip=triangle.rectangle;
                    if(direction>0)strip.lower.x=static_cast<double>(static_cast<long double>(strip.upper.x)-depth);
                    else strip.upper.x=static_cast<double>(static_cast<long double>(strip.lower.x)+depth);
                    const auto liquid=meridionalPlaneStripReference(triangle,strip,true)[0],gas=meridionalPlaneStripReference(triangle,strip,false)[0];
                    const auto independentMinor=minor?gas:liquid;
                    check(total>0&&flux.limiter==1&&flux.exhaustionLiquidCorrection==0&&flux.exhaustionGasCorrection==0,
                        "NATIVE_INTERFACE_SUPPORT geometric probe was changed by a limiter or whole-phase correction");
                    near(flux.liquid,direction*liquid,total,"NATIVE_INTERFACE_SUPPORT independent signed actual-plane liquid strip",1e-12L);
                    near(flux.gas,direction*gas,total,"NATIVE_INTERFACE_SUPPORT independent signed actual-plane gas strip",1e-12L);
                    near(liquid+gas,total,total,"NATIVE_INTERFACE_SUPPORT independently measured strip total",1e-12L);
                    near(flux.mass,100*flux.liquid+flux.gas,100*std::abs(flux.liquid)+std::abs(flux.gas),"NATIVE_INTERFACE_SUPPORT shared geometric mass flux",1e-12L);
                    check(direction>0?flux.negativeTotal==0&&flux.negativeLiquid==0&&flux.negativeGas==0:
                        flux.positiveTotal==0&&flux.positiveLiquid==0&&flux.positiveGas==0,
                        "NATIVE_INTERFACE_SUPPORT signed shared strip acquired an opposite-direction gross flux");
                    if(direction>0){const auto expectedFraction=half?7.L/16:3.L/4;const auto expected=expectedFraction*triangleArea;
                        near(depth,legX*(half?.25L:.5L),legX,"NATIVE_INTERFACE_SUPPORT declared half/full triangle penetration",1e-12L);
                        near(independentMinor,expected,triangleArea,"NATIVE_INTERFACE_SUPPORT independent analytic triangle swept fraction",1e-12L);
                        near(actualMinor,expected,triangleArea,"NATIVE_INTERFACE_SUPPORT actual resolved triangle outgoing phase",1e-12L);
                        peakTriangleRelative=std::max(peakTriangleRelative,std::abs(actualMinor-expected)/triangleArea);++triangleFaces;
                    }else{check(strip.upper.x<triangle.rectangle.upper.x-legX&&independentMinor==0&&actualMinor==0&&signedMinor==0,
                            "NATIVE_INTERFACE_SUPPORT noncontacting resolved strip must retain strictly zero minor-phase flux");
                        const auto historicalP0=total*(triangle.inventory[minor]/(triangle.inventory[0]+triangle.inventory[1]));
                        check(historicalP0>0,"NATIVE_INTERFACE_SUPPORT strict-zero witness cannot distinguish the retired forced-P0 outgoing tail");++strictZeroFaces;}
                    ++geometricFaces;
                }
            }
            check(named==1,"NATIVE_INTERFACE_SUPPORT geometric probe lacks its unique actual shared native face");
            for(bool half:{true,false}){const auto& target=half?moved.midpoint():moved.endpoint();
                for(int ci:{source,receiver}){using Exact=LAMP::Geometry::Detail::FiniteVolume::Dyadic;Exact expected(original[ci].inventory[minor]);
                    // This minority is the frozen WY tracked phase: its
                    // compression coefficient is zero in both directions.
                    for(const auto& face:moved.faces())if(face.negativeCell==ci||face.positiveCell==ci){const auto& flux=half?face.half:face.full;
                        expected.add(minor?flux.gas:flux.liquid,face.negativeCell==ci?-1:1);}
                    check(target.plicCells()[ci].inventory[minor]==expected.rounded(),
                        "NATIVE_INTERFACE_SUPPORT shared donor/receiver stock differs from its independent exact signed minority ledger");++sharedLedgers;}}
        }
        for(int step=1;step<=4;++step){const auto before=accepted.plicCells();const auto id=accepted.identity();const auto clock=accepted.timeLayer();
            const auto slab=Slab::prepare(mesh,accepted,{step,clock.time(),step*dt,Position::End},velocity,velocity);
            const auto transport=Transport::prepareLocalPlicSweep(mesh,accepted,slab,100,1,config);transport.validateFor(mesh,accepted,slab,100,1);
            inspect(transport.midpoint());inspect(transport.endpoint());unchanged(accepted,before,id,clock);
            accepted=transport.endpoint();++stepsTaken;}
        const auto retained=accepted.plicCells()[source].inventory[minor]/initialSource;
        maximumRetained=std::max(maximumRetained,retained);
        std::cout<<std::setprecision(17)<<"NATIVE_INTERFACE_SUPPORT_CORNER_LOCK_LIMITATION,minor_phase,"<<minor<<",steps,4,dt,"<<dt
            <<",source_minor_retained_fraction,"<<retained<<",east_minor_stock_change,"<<accepted.plicCells()[east].inventory[minor]-initialEast
            <<",diagnostic_only_no_four_step_release_or_motion_accuracy_acceptance,1\n";++scenes;
    }
    check(scenes==2&&stepsTaken==8&&geometricFaces==8&&strictZeroFaces==4&&triangleFaces==4&&sharedLedgers==16&&closures==6,
        "NATIVE_INTERFACE_SUPPORT resolved geometric/strict-zero/shared-ledger/force symmetric coverage incomplete");
    std::cout<<std::setprecision(17)<<"NATIVE_INTERFACE_SUPPORT_RESOLVED_GEOMETRY_FORCE_PASS,scenes,"<<scenes<<",diagnostic_steps,"<<stepsTaken
        <<",half_full_geometric_faces,"<<geometricFaces<<",strict_noncontact_zero_faces,"<<strictZeroFaces<<",analytic_triangle_fraction_faces,"<<triangleFaces<<",shared_phase_ledgers,"<<sharedLedgers
        <<",unavailable_physical_curvature_closures,"<<closures<<",maximum_source_minority_retained_fraction,"<<maximumRetained
        <<",peak_triangle_relative_error,"<<peakTriangleRelative<<",maximum_two_phase_volume_error,"<<maxConservation
        <<",no_stock_deletion_or_plane_change,1,no_P0_or_motion_accuracy_claim,1\n";
    opposedFilmForce();
}
}

namespace TranslatedCircleShapeContract
{
long double area(Box box,long double cx,long double cy,long double radius)
{
    const auto capacity=(box.x1-box.x0)*(box.y1-box.y0);
    const auto farX=std::max(std::abs(box.x0-cx),std::abs(box.x1-cx)),farY=std::max(std::abs(box.y0-cy),std::abs(box.y1-cy));
    if(farX*farX+farY*farY<=radius*radius)return capacity;
    const auto nearX=std::max({box.x0-cx,0.L,cx-box.x1}),nearY=std::max({box.y0-cy,0.L,cy-box.y1});
    if(nearX*nearX+nearY*nearY>=radius*radius)return 0;
    const auto left=std::max(box.x0,cx-radius),right=std::min(box.x1,cx+radius);if(right<=left)return 0;
    std::vector<long double> cuts{left,right};
    for(auto y:{box.y0,box.y1})if(std::abs(y-cy)<radius){const auto dx=std::sqrt(radius*radius-(y-cy)*(y-cy));
        for(auto x:{cx-dx,cx+dx})if(x>left&&x<right)cuts.push_back(x);}
    std::sort(cuts.begin(),cuts.end());cuts.erase(std::unique(cuts.begin(),cuts.end()),cuts.end());
    const auto primitive=[&](long double x){const auto t=std::clamp((x-cx)/radius,-1.L,1.L);
        return radius*radius*(t*std::sqrt(std::max(0.L,1-t*t))+std::asin(t))/2;};
    long double result=0;
    for(std::size_t i=1;i<cuts.size();++i){const auto a=cuts[i-1],b=cuts[i],mid=(a+b)/2;
        const auto dy=std::sqrt(std::max(0.L,radius*radius-(mid-cx)*(mid-cx))),lower=cy-dy,upper=cy+dy;
        if(std::min(box.y1,upper)<=std::max(box.y0,lower))continue;
        long double constant=0,coefficient=0;
        if(box.y1<upper)constant+=box.y1;else{constant+=cy;++coefficient;}
        if(box.y0>lower)constant-=box.y0;else{constant-=cy;++coefficient;}
        result+=constant*(b-a)+coefficient*(primitive(b)-primitive(a));
    }return result;
}
void run()
{
    constexpr long double cx=.487L,cy=.513L,radius=.23L,ux=.7L,uy=.4L,finish=.02L;
    const auto totalLiquid=std::acos(-1.L)*radius*radius;std::array<long double,2> errors{};unsigned scenes=0,acceptedSteps=0;long double maxConservation=0;
    for(int n:{32,64}){stage="NATIVE_INTERFACE_SUPPORT independent translating circle "+std::to_string(n);const unsigned scene=n==32?0u:1u;
        const double dt=n==32?.005:.0025;const int steps=n==32?4:8;
        const auto pathRadius=std::max(std::hypot(cx-.5L,cy-.5L),std::hypot(cx+ux*finish-.5L,cy+uy*finish-.5L));
        check(pathRadius+radius+std::sqrt(2.L)/n<.3125L,"NATIVE_INTERFACE_SUPPORT_SHAPE circle-path intersected cell faces leave the exactly constant-velocity core");
        DictionaryMesh dict;dict.Nx=dict.Ny=n;dict.maxInitialIterations=0;dict.maxCellNum=4*n*n;QuadtreeMesh mesh(dict);
        std::vector<std::array<long double,2>> stocks(mesh.cellIndexNum());std::array<long double,2> initial{};
        for(int ci=0;ci<mesh.cellIndexNum();++ci){const auto b=bounds(mesh,ci);const auto c=capacity(b,false),liquid=area(b,cx,cy,radius);
            check(liquid>=0&&liquid<=c,"NATIVE_INTERFACE_SUPPORT_SHAPE independent initial circle integral outside capacity");stocks[ci]={liquid,c-liquid};
            for(unsigned phase=0;phase<2;++phase)initial[phase]+=stocks[ci][phase];}
        near(initial[0],totalLiquid,1,"NATIVE_INTERFACE_SUPPORT_SHAPE independent initial complete circle area",1e-12L);
        GeometricVofConfig resolution;resolution.localResolution=GeometricVofLocalResolution::NativeInterfaceSupport;
        auto accepted=Publication::prepareLocalPlicInventory(mesh,domain(mesh,false),{},stocks,resolution);
        // The gauge constant changes no physical velocity. Compact radial
        // taper closes every wall; all circle-intersected native face endpoints
        // stay in its exactly linear core for the entire declared translation.
        const auto psi=[](Vector2D p){const auto r=std::hypot(static_cast<long double>(p.x)-.5L,static_cast<long double>(p.y)-.5L);
            const auto factor=r<=.3125L?1.L:r>=.5L?0.L:(.5L-r)/.1875L;
            return (.7L*(p.y-.5L)-.4L*(p.x-.5L))*factor;};
        std::vector<double> velocity(accepted.domain().nativeFaceCount());
        for(const auto f:accepted.domain().faces(mesh)){const auto& edge=f.horizontal?mesh.hEdges(f.edgeLabel):mesh.vEdges(f.edgeLabel);
            auto a=mesh.verts(edge.edgeToVert[0]).coord(),b=mesh.verts(edge.edgeToVert[1]).coord();if(f.horizontal?a.x>b.x:a.y>b.y)std::swap(a,b);
            const auto flux=(f.horizontal?-1.L:1.L)*(psi(b)-psi(a));
            if(f.outwardSign!=0){check(flux==0,"NATIVE_INTERFACE_SUPPORT_SHAPE closed streamfunction leaked through a wall");velocity[f.nativeRow]=0;}
            else velocity[f.nativeRow]=static_cast<double>(flux/f.area);}
        const auto inspect=[&](const Publication& p){std::array<long double,2> total{};
            for(std::size_t ci=0;ci<p.plicCells().size();++ci){const auto& q=p.plicCells()[ci];const auto c=p.cells()[ci].capacity;
                using Exact=LAMP::Geometry::Detail::FiniteVolume::Dyadic;
                const auto& receipt=p.capacityDefects()[ci];Exact closure(q.inventory[0]);closure.add(q.inventory[1]);closure.add(c,-1);
                auto bound=Exact(c);bound.add(receipt.absolute());const auto ceiling=Exact::product({1e-12L,c});
                // Typed successors carry the producer's exact closure receipt.
                // A certified pure stock can exceed C by that receipt; preserve
                // the raw stock and verify the original capacity budget exactly.
                check(closure.compare(receipt)==0,"NATIVE_INTERFACE_SUPPORT_SHAPE independent exact closure differs from issued capacity receipt");
                check(receipt.absolute().compare(ceiling)<=0,"NATIVE_INTERFACE_SUPPORT_SHAPE issued capacity receipt exceeds original 1e-12 C budget");
                for(unsigned phase=0;phase<2;++phase){std::ostringstream diagnostic;
                    diagnostic<<std::setprecision(std::numeric_limits<long double>::max_digits10)
                        <<"NATIVE_INTERFACE_SUPPORT_SHAPE independent stock escaped original bounds native="<<ci<<" step="<<p.timeLayer().step
                        <<" position="<<static_cast<int>(p.timeLayer().position)<<" begin="<<p.timeLayer().begin<<" end="<<p.timeLayer().end
                        <<" phase="<<phase<<" Q="<<q.inventory[phase]<<" C="<<c<<" Q_minus_C="<<q.inventory[phase]-c
                        <<" pair="<<q.inventory[0]<<'/'<<q.inventory[1]<<" capacityDefect="<<receipt.rounded()<<" independent_exact_closure="<<closure.rounded()
                        <<" exact_closure_vs_receipt="<<closure.compare(receipt)<<" abs_receipt_vs_original_budget="<<receipt.absolute().compare(ceiling)
                        <<" Q_vs_C_plus_abs_receipt="<<Exact(q.inventory[phase]).compare(bound)<<" Q_hex="<<std::hexfloat<<q.inventory[phase]
                        <<" C_hex="<<c<<" difference_hex="<<q.inventory[phase]-c<<" receipt_hex="<<receipt.rounded();
                    check(std::isfinite(q.inventory[phase])&&q.inventory[phase]>=0&&Exact(q.inventory[phase]).compare(bound)<=0,diagnostic.str());
                    total[phase]+=q.inventory[phase];near(q.measured[phase],q.inventory[phase],c,"NATIVE_INTERFACE_SUPPORT_SHAPE inverse PLIC volume",1e-12L);}}
            for(unsigned phase=0;phase<2;++phase){maxConservation=std::max(maxConservation,std::abs(total[phase]-initial[phase]));
                near(total[phase],initial[phase],1,"NATIVE_INTERFACE_SUPPORT_SHAPE both-phase volume conservation",1e-12L);}};
        GeometricVofTransportConfig curved;curved.geometry=resolution;
        for(int step=1;step<=steps;++step){const auto saved=accepted.plicCells();const auto identity=accepted.identity();const auto clock=accepted.timeLayer();
            const Clock target{step,clock.time(),step*dt,Position::End};const auto slab=Slab::prepare(mesh,accepted,target,velocity,velocity);
            const auto transport=Transport::prepareLocalPlicSweep(mesh,accepted,slab,100,1,curved);transport.validateFor(mesh,accepted,slab,100,1);
            inspect(transport.midpoint());inspect(transport.endpoint());
            check(accepted.identity()==identity&&sameClock(accepted.timeLayer(),clock),"NATIVE_INTERFACE_SUPPORT_SHAPE candidate changed accepted owner or full clock");
            for(std::size_t ci=0;ci<saved.size();++ci)check(accepted.plicCells()[ci].inventory==saved[ci].inventory&&accepted.plicCells()[ci].normal.x==saved[ci].normal.x
                &&accepted.plicCells()[ci].normal.y==saved[ci].normal.y&&accepted.plicCells()[ci].offset==saved[ci].offset,"NATIVE_INTERFACE_SUPPORT_SHAPE candidate changed accepted inventories/planes");
            accepted=transport.endpoint();++acceptedSteps;}
        near(accepted.timeLayer().time(),finish,1,"NATIVE_INTERFACE_SUPPORT_SHAPE fixed final physical time");
        long double error=0;for(int ci=0;ci<mesh.cellIndexNum();++ci){const auto exact=area(bounds(mesh,ci),cx+ux*finish,cy+uy*finish,radius);
            error+=std::abs(accepted.plicCells()[ci].inventory[0]-exact);}
        errors[scene]=error/totalLiquid;const long double relativeSpacing=1.L/(n*radius);
        check(errors[scene]<=3*relativeSpacing*relativeSpacing,"NATIVE_INTERFACE_SUPPORT_SHAPE translating circle exceeds frozen independent L1/area gate");++scenes;
        std::cout<<std::setprecision(17)<<"NATIVE_INTERFACE_SUPPORT_TRANSLATED_CIRCLE_GRID,n,"<<n<<",steps,"<<steps<<",T,"<<finish<<",phase_L1_over_circle_area,"<<errors[scene]
            <<",frozen_gate,"<<3*relativeSpacing*relativeSpacing<<",global_complete_cfl,"<<dt*vofMaterialOutflow(mesh,accepted,velocity,velocity)<<'\n';
    }
    check(errors[1]<.6L*errors[0],"NATIVE_INTERFACE_SUPPORT_SHAPE independent translating-circle refinement ratio exceeds frozen .6");
    std::cout<<"NATIVE_INTERFACE_SUPPORT_TRANSLATED_CIRCLE_PASS,scenes,"<<scenes<<",accepted_steps,"<<acceptedSteps<<",L1_refinement_ratio,"<<errors[1]/errors[0]
        <<",maximum_phase_volume_error,"<<maxConservation<<",original_numerical_gates_unchanged,1\n";
}
}

int main(int argc,char** argv)
{
    try{
        if(argc==2&&std::string(argv[1])=="--local-plic-native-interface-support-only"){NativeInterfaceTransportContract::run();return 0;}
        if(argc==2&&std::string(argv[1])=="--local-plic-translated-circle-shape-only"){TranslatedCircleShapeContract::run();return 0;}
        if(argc==2&&std::string(argv[1])=="--local-plic-height-normal-transport-only"){StockHeightTransportContract::run();return 0;}
        if(argc==2&&std::string(argv[1])=="--local-plic-full-material-only"){MaterialFullOnlyContract::run();return 0;}
        if(argc==2&&std::string(argv[1])=="--local-plic-filtered-update-only"){VofFilteredLocalUpdateContract::run();return 0;}
        if(argc==2&&std::string(argv[1])=="--pinned-corner-transport-only"){PinnedCornerTransportContract::run();return 0;}
        if(argc==2&&std::string(argv[1])=="--selected-fluid-ehd-step-only"){SelectedFluidEhdStepContract::run();return 0;}
        if(argc==2&&std::string(argv[1])=="--selected-fluid-hydro-only"){SelectedFluidHydroContract::run();return 0;}
        if(argc==2&&std::string(argv[1])=="--selected-fluid-transport-only"){SelectedFluidTransportContract::run();return 0;}
        if(argc==2&&std::string(argv[1])=="--local-plic-conservative-momentum-only"){ConservativeVofMomentumContract::run();return 0;}
        if(argc==2&&std::string(argv[1])=="--local-plic-face-parallel-only"){F51LocalFaceParallelContract::run();return 0;}
        if(argc==2&&std::string(argv[1])=="--transport-validation-batch-only"){TransportValidationBatchContract::run();return 0;}
        if(argc==2&&std::string(argv[1])=="--local-plic-open-momentum-step-only"){F17OpenMomentumStepContract::run();return 0;}
                if(argc==2&&std::string(argv[1])=="--local-plic-open-momentum-only"){F16OpenMomentumContract::run();return 0;}
        if(argc==2&&std::string(argv[1])=="--local-plic-p0-boundary-gross-only"){F18P0BoundaryGrossContract::run();return 0;}
        if(argc==2&&std::string(argv[1])=="--local-plic-boundary-nonreverse-only"){F15BoundaryGrossContract::nonreverseOnly();return 0;}
        if(argc==2&&std::string(argv[1])=="--local-plic-boundary-gross-only"){F15BoundaryGrossContract::run();return 0;}
        if(argc==2&&std::string(argv[1])=="--local-plic-open-step-only"){localPlicOpenStepOnly();return 0;}
        if(argc==2&&std::string(argv[1])=="--local-plic-open-kep-only"){localPlicOpenKepOnly();return 0;}
        if(argc==2&&std::string(argv[1])=="--local-plic-open-boundary-only"){localPlicOpenBoundaryOnly();return 0;}
        if(argc==2&&std::string(argv[1])=="--local-plic-projected-capacity-only"){localPlicProjectedCapacityOnly();return 0;}
        if(argc==2&&std::string(argv[1])=="--local-plic-capacity-carry-only"){localPlicCapacityCarryOnly();return 0;}
        if(argc==2&&std::string(argv[1])=="--local-plic-meridional-only"){localPlicMeridionalOnly();return 0;}
        // BEGIN L1 LOCAL CLI ENTRIES
        // BEGIN D11 DEDICATED CLI ENTRY
        if((argc==2||argc==3)&&std::string(argv[1])=="--local-plic-convection-range-only"){
            localPlicConvectionRangeOnly(argc==3?std::stoi(argv[2]):120);return 0;}
        if(argc==2&&std::string(argv[1])=="--vof-viscosity-capillary-overlap-only"){vofViscosityCapillaryOverlapOnly();return 0;}
        if(argc==2&&std::string(argv[1])=="--local-plic-rz-capillary-step-only"){localPlicRzCapillaryStepOnly();return 0;}
        if(argc==2&&std::string(argv[1])=="--local-plic-finite-p0-physical-density-only"){localPlicAtomsOnly(true,true);return 0;}
        if(argc==2&&std::string(argv[1])=="--local-plic-finite-p0-only"){localPlicAtomsOnly(true);return 0;}
        if(argc==2&&std::string(argv[1])=="--local-plic-convection-only"){localPlicAtomsOnly(true);return 0;}
        if(argc==2&&std::string(argv[1])=="--local-plic-atoms-only"){localPlicAtomsOnly();return 0;}
        // END D11 DEDICATED CLI ENTRY
        if(argc==2&&std::string(argv[1])=="--local-smooth-only"){localSmoothOnly();return 0;}
        if(argc==2&&std::string(argv[1])=="--local-plic-single-outflow-only"){localPlicSingleOutgoingOnly();return 0;}
        if(argc==2&&std::string(argv[1])=="--local-plic-only"){localPlicTransport();return 0;}
        if(argc==5&&std::string(argv[1])=="--local-plic-cost"){localPlicCost(std::stoi(argv[2]),std::stoi(argv[3]),std::stod(argv[4]));return 0;}
        check(argc==1,"unknown geometric VOF transport contract argument");
        // END L1 LOCAL CLI ENTRIES
        unsigned failures=0;
        const std::array<std::pair<const char*,void(*)()>,10> groups{{
            {"correctly-rounded-rationals",correctlyRoundedRationals},
            {"captured-complete-budget",capturedBudgetReplay},
            {"core-face-time-original-geometry",coreFaceTimeOriginalGeometry},
            {"microscopic-phase-budget",microscopicPhaseBudgets},
            {"general-native-xy-rz",genericMatrix},{"constant-diagonal",constantDiagonalPatch},
            {"pure-phases",pureMatrix},{"coupled-limiter",limiterMatrix},
            {"sequential-rebuild",sequentialCandidates},{"guards-and-fields",failuresAndFields}}};
        for(const auto& group:groups){
            try{group.second();std::cout<<"group="<<group.first<<" PASS\n";}
            catch(const std::exception& error){++failures;std::cerr<<"group="<<group.first<<" FAIL stage="<<stage<<" "<<error.what()<<"\n";}
        }
        // BEGIN V4a1a ADDED FACETIME MAIN CALL
        try{coreAnchoredFaceTime();std::cout<<"group=core-anchored-static-fubini PASS\n";}
        catch(const std::exception& error){++failures;std::cerr<<"group=core-anchored-static-fubini FAIL "<<error.what()<<"\n";}
        // END V4a1a ADDED FACETIME MAIN CALL
        // BEGIN V4a1b ADDED MAIN CALL
        try{originalAnchoredNativePipeline();std::cout<<"group=actual-anchored-native-pipeline PASS\n";}
        catch(const std::exception& error){++failures;std::cerr<<"group=actual-anchored-native-pipeline FAIL "<<error.what()<<"\n";}
        // END V4a1b ADDED MAIN CALL
        // BEGIN V4b1 ADDED MAIN CALL
        try{v4b1NativeRestriction();std::cout<<"group=retained-native-subface PASS\n";}
        catch(const std::exception& error){++failures;std::cerr<<"group=retained-native-subface FAIL "<<error.what()<<'\n';}
        // END V4b1 ADDED MAIN CALL
        // BEGIN V4b2 ADDED PRODUCER MAIN CALL
        try{v4b2InteriorProducer();std::cout<<"group=analytic-vof-interior-producer PASS\n";}
        catch(const std::exception& error){++failures;std::cerr<<"group=analytic-vof-interior-producer FAIL "<<error.what()<<'\n';}
        // END V4b2 ADDED PRODUCER MAIN CALL
        // BEGIN L1 DEFAULT TRANSPORT GROUP
        try{localPlicTransport();std::cout<<"group=local-plic-sweep PASS\n";}
        catch(const std::exception& error){++failures;std::cerr<<"group=local-plic-sweep FAIL "<<error.what()<<'\n';}
        // END L1 DEFAULT TRANSPORT GROUP
        try{StockHeightTransportContract::run();std::cout<<"group=stock-height-real-split PASS\n";}
        catch(const std::exception& error){++failures;std::cerr<<"group=stock-height-real-split FAIL "<<error.what()<<'\n';}
        try{NativeInterfaceTransportContract::run();std::cout<<"group=native-interface-support-resolved-geometry-force PASS\n";}
        catch(const std::exception& error){++failures;std::cerr<<"group=native-interface-support-resolved-geometry-force FAIL "<<error.what()<<'\n';}
        check(coarseChecks>0,"no actual 2:1 coarse/fine flux sum checks");
        std::cout<<std::setprecision(17)<<"geometric_vof_transport_contract "<<(failures==0?"PASS":"FAIL")
            <<" failed_groups="<<failures<<" scenarios="<<scenarios
            <<" assertions="<<assertions<<" guards="<<guards<<" raw_checks="<<rawChecks<<" coarse_checks="<<coarseChecks
            <<" counterflow_faces="<<counterflowFaces<<" limited_faces="<<limitedFaces
            <<" rational_round_scenarios="<<rationalRoundScenarios<<" rational_round_rejections="<<rationalRoundRejections
            <<" captured_budget_scenarios="<<capturedBudgetScenarios<<" captured_budget_gross_checks="<<capturedBudgetGrossChecks
            <<" captured_budget_pure_checks="<<capturedBudgetPureChecks
            <<" micro_scenarios="<<microScenarios<<" peak_micro_relative="<<peakMicroRelative
            <<" face_time_attempts="<<faceTimeAttempts<<" face_time_scenarios="<<faceTimeScenarios
            <<" face_time_positive_checks="<<faceTimePositiveChecks<<" face_time_exact_zero_checks="<<faceTimeZeroChecks
            <<" peak_face_time_micro_relative="<<peakFaceTimeRelative<<" peak_face_time_uncertainty="<<peakFaceTimeUncertainty
            <<" maximum_face_time_root_refinements="<<maximumFaceTimeRootRefinements
            <<" peak_stock="<<peakInventory<<" peak_conservation="<<peakConservation<<" peak_raw="<<peakRaw
            <<" peak_geometry="<<peakGeometry<<" peak_flux_correction="<<peakCorrection<<"\n";
        return failures==0?0:1;
    }catch(const std::exception& e){std::cerr<<"geometric_vof_transport_contract FAIL stage="<<stage<<" "<<e.what()<<"\n";return 1;}
}
