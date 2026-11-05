#include "Modules/Hydrodynamics/HydroTwoPhaseStep.h"
#include "Modules/Hydrodynamics/HydroVofStep.h"
#include "Modules/Hydrodynamics/HydroVofMeshTransfer.h"
#include "Core/Mesh/CartesianMeasure.h"
#include "interface_first_jet_fixture.h"
#include <algorithm>
#include <sstream>
#include <filesystem>
#include <fstream>
#include <map>
#include <chrono>
#include <cfenv>
#include <cerrno>
#include <cstring>
#if defined(_M_X64) || defined(__SSE2__)
#include <xmmintrin.h>
#endif
#ifdef _OPENMP
#include <omp.h>
#endif

namespace
{
using namespace LAMP::Hydrodynamics;
using Matrix=HydroTwoPhaseStep::Matrix;
int movingCases=0,acceptedSteps=0,guardCases=0;
int capillarySteps=0,capillaryRejects=0;
int relaxedComparisons=0,updateToleranceRejects=0;
int linearCases=0,linearRejects=0;
long double peakEquation=0,peakEnergy=0,peakMaterial=0;

HydroPhysicalDivergence stepMetrics(const QuadtreeMesh& mesh,bool rz=false,bool subset=false)
{
    FaceSegmentGraph graph;
    graph.build(mesh,[subset](const QuadtreeMesh&,const Cell& c){return !subset||c.centerCoord().x<.5;});
    return HydroPhysicalDivergence::build(mesh,HydroIntegralDivergence::build(mesh,
        HydroDofLayout::build(mesh,NativeLeafDomain::build(mesh,graph,rz))));
}

// The integral curl uses only physical endpoints, so every cell has a
// telescoping zero outward flux. The normal wall data are exactly zero.
double stepPsi(Vector2D p)
{
    const double x=p.x*(1-p.x),y=p.y*(1-p.y);
    // Break the square's reflection symmetries: on n=4 the symmetric
    // vortex can be a discrete steady Euler mode and hide a frozen update.
    return 16*x*x*y*y*(1+.3*p.x+.7*p.y);
}
HydroNativeValues stepVelocity(const QuadtreeMesh& mesh,const HydroDofLayout& layout,bool stationary=false)
{
    std::vector<double> u(layout.velocityCount());
    for(int f=0;f<layout.velocityCount();++f){
        const auto face=layout.domain().face(mesh,layout.nativeFace(f));
        const auto& edge=face.horizontal?mesh.hEdges(face.edgeLabel):mesh.vEdges(face.edgeLabel);
        auto a=mesh.verts(edge.edgeToVert[0]).coord(),b=mesh.verts(edge.edgeToVert[1]).coord();
        if(face.horizontal?a.x>b.x:a.y>b.y)std::swap(a,b);
        if(!stationary)u[f]=(face.horizontal?stepPsi(a)-stepPsi(b):stepPsi(b)-stepPsi(a))/face.tangentSpacing;
    }
    return HydroNativeValues(layout,u,std::vector<double>(layout.pressureCount()));
}

// Independent physical geometry: these formulas do not consume the metrics
// under test. XY retains its original Cartesian measure and all old gates.
long double stepBoxMeasure(Vector2D lo,Vector2D hi,bool rz)
{
    if(!rz)return static_cast<long double>(hi.x-lo.x)*(hi.y-lo.y);
    return pi*(static_cast<long double>(hi.x)*hi.x-static_cast<long double>(lo.x)*lo.x)*(hi.y-lo.y);
}
long double stepFaceMeasure(const QuadtreeMesh& mesh,const HydroDofLayout& l,int u)
{
    const auto face=l.domain().face(mesh,l.nativeFace(u));
    if(!l.domain().axisymmetric())return face.tangentSpacing;
    const auto& edge=face.horizontal?mesh.hEdges(face.edgeLabel):mesh.vEdges(face.edgeLabel);
    auto a=mesh.verts(edge.edgeToVert[0]).coord(),b=mesh.verts(edge.edgeToVert[1]).coord();
    if(face.horizontal?a.x>b.x:a.y>b.y)std::swap(a,b);
    return face.horizontal?pi*(static_cast<long double>(b.x)*b.x-static_cast<long double>(a.x)*a.x)
        :2*pi*a.x*(b.y-a.y);
}
HydroNativeValues rzStepVelocity(const QuadtreeMesh& mesh,const HydroDofLayout& l)
{
    const auto psi=[](Vector2D p){const double r2=p.x*p.x,z=p.y*(1-p.y);
        return 16*r2*(1-r2)*(1-r2)*z*z*(1+.7*p.y);};
    std::vector<double> u(l.velocityCount());
    for(int i=0;i<l.velocityCount();++i){const auto face=l.domain().face(mesh,l.nativeFace(i));
        const auto& edge=face.horizontal?mesh.hEdges(face.edgeLabel):mesh.vEdges(face.edgeLabel);
        auto a=mesh.verts(edge.edgeToVert[0]).coord(),b=mesh.verts(edge.edgeToVert[1]).coord();
        if(face.horizontal?a.x>b.x:a.y>b.y)std::swap(a,b);
        const auto area=stepFaceMeasure(mesh,l,i);
        u[i]=area>0?static_cast<double>(2*pi*(face.horizontal?psi(a)-psi(b):psi(b)-psi(a))/area):0;
    }
    return {l,u,std::vector<double>(l.pressureCount())};
}

long double scaled(long double actual,long double expected,long double scale,const char* why,double tolerance=1e-11)
{
    const long double gap=std::abs(actual-expected),relative=scale==0?(gap==0?0:std::numeric_limits<long double>::infinity()):gap/scale;
    if(!std::isfinite(relative)||relative>tolerance){std::ostringstream message;message<<std::setprecision(17)<<why
        <<" actual="<<actual<<" expected="<<expected<<" scale="<<scale<<" relative="<<relative;throw std::runtime_error(message.str());}
    return relative;
}
std::vector<double> matvec(const Matrix& a,const std::vector<double>& u)
{
    const Eigen::Map<const Eigen::VectorXd> x(u.data(),u.size());const Eigen::VectorXd y=a*x;
    return {y.data(),y.data()+y.size()};
}

// Copies every primary/closure/native payload. Addresses additionally ensure
// rejection preserves the retained authority rather than an equivalent clone.
struct StateCopy
{
    HermiteGalsIdentity identity;
    HermiteGalsTimeLayer clock;
    HermiteGalsPrimaryState primary;
    HermiteGalsClosure closure;
    std::vector<double> u,p;
    const HermiteGalsPrimaryState* primaryAddress=nullptr;
    const LevelSetHermiteGeometry::Snapshot* geometryAddress=nullptr;
};
StateCopy capture(const HermiteGalsPublication& pub,const HydroNativeValues& native)
{
    const auto& s=pub.state();return {s.identity(),s.timeLayer(),s.primary(),s.closure(),native.normalVelocity(),native.pressure(),
        &s.primary(),&s.geometry().snapshot()};
}
bool sameClock(HermiteGalsTimeLayer a,HermiteGalsTimeLayer b)
{return a.step==b.step&&a.begin==b.begin&&a.end==b.end&&a.position==b.position;}
void unchanged(const StateCopy& saved,const HermiteGalsPublication& pub,const HydroNativeValues& native)
{
    const auto& s=pub.state();const auto& p=s.primary();const auto& c=s.closure();
    check(saved.identity==s.identity()&&sameClock(saved.clock,s.timeLayer()),"failed step changed accepted identity or complete clock");
    check(saved.primaryAddress==&p&&saved.geometryAddress==&s.geometry().snapshot(),"failed step replaced old geometry authority");
    check(saved.u==native.normalVelocity()&&saved.p==native.pressure(),"failed step changed native values");
    check(saved.primary.measure==p.measure&&saved.primary.transport==p.transport&&saved.primary.geometryCells==p.geometryCells,
        "failed step changed primary topology/transport");
    check(saved.primary.vertices.size()==p.vertices.size()&&saved.primary.edges.size()==p.edges.size()
        &&saved.primary.phaseVolumes.size()==p.phaseVolumes.size(),"failed step changed primary dimensions");
    for(std::size_t i=0;i<p.vertices.size();++i){const auto& a=saved.primary.vertices[i];const auto& b=p.vertices[i];
        check(a.vertex==b.vertex&&a.phi==b.phi&&a.dx==b.dx&&a.dy==b.dy,"failed step changed vertex FirstJet");}
    for(std::size_t i=0;i<p.edges.size();++i){const auto& a=saved.primary.edges[i];const auto& b=p.edges[i];
        check(a.vertices==b.vertices&&a.phi==b.phi&&a.dn==b.dn,"failed step changed edge FirstJet");}
    for(std::size_t i=0;i<p.phaseVolumes.size();++i){const auto& a=saved.primary.phaseVolumes[i];const auto& b=p.phaseVolumes[i];
        check(a.cell==b.cell&&a.liquid==b.liquid,"failed step changed forbidden inventory");}
    check(saved.closure.vertexMixedDerivatives==c.vertexMixedDerivatives&&saved.closure.cellBubbleCoefficients==c.cellBubbleCoefficients,
        "failed step changed closure fields");
}

void inspectIteration(const HydroTwoPhaseIterationReport& r,int expected,HermiteGalsTimeLayer endClock,std::size_t expectedInitialPanels)
{
    check(r.iteration==expected&&r.geometryCompared==(expected>1),"iteration callback lost the actual ordered Picard history");
    const HermiteGalsTimeLayer midpointClock{endClock.step,endClock.begin,endClock.end,Position::Midpoint};
    check(sameClock(r.timeLayer,midpointClock),"iteration callback lost the actual complete transport clock");
    check(r.reusedTimePlan==(expected>1)&&r.initialPlanPanels==expectedInitialPanels,
        "Picard time plan was not cold first, then the actual preceding accepted partition");
    check(std::isfinite(r.interfaceMismatch)&&r.interfaceMismatch>=0
        &&std::isfinite(r.linearMomentumDefect)&&r.linearMomentumDefect>=0,"iteration callback has invalid linear stage diagnostics");
    check(std::isfinite(r.frozenOwnerMomentumDefect)&&r.frozenOwnerMomentumDefect>=0,"iteration callback has an invalid actual nonlinear momentum defect");
    check(std::isfinite(r.momentumPowerRelativeDefect)&&r.momentumPowerRelativeDefect>=0,"iteration callback has an invalid actual physical momentum power defect");
    for(double seconds:r.stageSeconds)check(std::isfinite(seconds)&&seconds>=0,"iteration callback has an invalid stage time");
    check(std::isfinite(r.csfColorSeconds)&&r.csfColorSeconds>=0
        &&std::isfinite(r.csfProjectionSeconds)&&r.csfProjectionSeconds>=0,"iteration callback has invalid CSF split times");
    check(r.geometrySamples>=2&&r.acceptedPanels>0&&r.deepestPanel>=0,"iteration callback has no actual transport work inventory");
}
void inspectIterationObserver(const HydroTwoPhaseStep& step,const std::vector<HydroTwoPhaseIterationReport>& observed)
{
    const auto& history=step.report().history;
    check(observed.size()==history.size(),"iteration callback skipped or repeated a completed Picard iteration");
    for(std::size_t i=0;i<history.size();++i){
        const auto& a=observed[i];const auto& b=history[i];inspectIteration(a,static_cast<int>(i)+1,step.candidateState().timeLayer(),
            i==0?1:history[i-1].acceptedPanels);
        check(a.iteration==b.iteration&&a.timeScheme==b.timeScheme&&a.timeScheme==step.report().timeScheme
            &&a.interfaceMismatch==b.interfaceMismatch&&a.linearMomentumDefect==b.linearMomentumDefect
            &&a.geometryCompared==b.geometryCompared&&sameClock(a.timeLayer,b.timeLayer)
            &&a.velocityUpdateRelative==b.velocityUpdateRelative&&a.geometryUpdateRelative==b.geometryUpdateRelative
            &&a.midpointGeometryUpdateRelative==b.midpointGeometryUpdateRelative&&a.endGeometryUpdateRelative==b.endGeometryUpdateRelative
            &&a.stageSeconds==b.stageSeconds&&a.csfColorSeconds==b.csfColorSeconds
            &&a.csfProjectionSeconds==b.csfProjectionSeconds&&a.geometrySamples==b.geometrySamples
            &&a.acceptedPanels==b.acceptedPanels&&a.deepestPanel==b.deepestPanel
            &&a.initialPlanPanels==b.initialPlanPanels&&a.reusedTimePlan==b.reusedTimePlan
            &&a.endpointSupportPanels==b.endpointSupportPanels&&a.frozenOwnerMomentumDefect==b.frozenOwnerMomentumDefect
            &&a.momentumPowerRelativeDefect==b.momentumPowerRelativeDefect,
            "iteration callback differs from the retained actual history");
        const auto& x=a.saddle;const auto& y=b.saddle;
        check(x.algebraicRelativeResidual==y.algebraicRelativeResidual&&x.freeMomentumDefect==y.freeMomentumDefect
            &&x.maxCellContinuityDefect==y.maxCellContinuityDefect&&x.maxPressureMeanDefect==y.maxPressureMeanDefect
            &&x.maxBoundaryFluxDefect==y.maxBoundaryFluxDefect&&x.maxAbsGaugeMultiplier==y.maxAbsGaugeMultiplier
            &&x.freeVelocityCount==y.freeVelocityCount&&x.pressureGaugeCount==y.pressureGaugeCount
            &&x.freeTangentialCount==y.freeTangentialCount&&x.linearRefinements==y.linearRefinements,"iteration callback changed the actual saddle report");
    }
    check(!history.empty(),"successful step has no completed iteration");
    const auto& actual=step.materialTransport().spaceTimeReport();const auto& last=history.back();
    const auto& cost=actual.samplingCost;const auto& transport=step.materialTransport();
    const auto& motion=transport.characteristicSource();
    check(cost.atomLocal==!transport.layout().domain().axisymmetric()
        &&cost.sourceIdentity==transport.beginGeometry().geometryIdentity()
        &&cost.previousAdvectorIdentity==motion.previous().identity()&&cost.arrivalAdvectorIdentity==motion.arrival().identity(),
        "coupled local time integration detached its actual source/endpoint owners");
    check(cost.fullSamples>=3&&cost.fullPreparedNativeCells>0
        &&cost.fullAtomRateEvaluations+cost.localAtomRateEvaluations>=3*transport.phaseAtomIntegrals().size(),
        "coupled time integration did not report its actually evaluated full geometry/rates");
    std::cout<<"STEP_PHASE_SAMPLING,"<<cost.localSamples<<','<<cost.fullSamples<<','<<cost.localPreparedNativeCells
        <<','<<cost.fullPreparedNativeCells<<','<<cost.localAtomRateEvaluations<<','<<cost.fullAtomRateEvaluations<<",PASS\n";
    check(sameClock(last.timeLayer,step.materialTransport().timeLayer())
        &&last.geometrySamples==actual.geometrySamples&&last.acceptedPanels==actual.acceptedPanels
        &&last.deepestPanel==actual.deepestPanel&&last.initialPlanPanels==actual.initialPlanPanels
        &&last.reusedTimePlan==actual.reusedTimePlan&&last.endpointSupportPanels==actual.endpointSupportPanels,
        "iteration counters do not belong to the retained actual transport");
    check(last.frozenOwnerMomentumDefect==step.report().frozenOwnerMomentumDefect
        &&last.momentumPowerRelativeDefect==step.report().momentumPowerRelativeDefect
        &&last.linearMomentumDefect==step.report().linearMomentumDefect,
        "final iteration nonlinear momentum/power defects differ from the actual step reassessment");
    std::cout<<"STEP_ITERATION_OBSERVER,"<<observed.size()<<','<<last.geometrySamples<<','<<last.acceptedPanels<<','<<last.deepestPanel<<",PASS\n";
}

template<class F>void stepRejects(F action,const char* expected=nullptr)
{
    try{action();}catch(const std::exception& e){
        if(expected)check(std::string(e.what()).find(expected)!=std::string::npos,
            std::string("rejection did not reach the required stage: expected=")+expected+", actual="+e.what());
        ++guardCases;std::cout<<"STEP_REJECT,"<<e.what()<<'\n';return;}
    throw std::runtime_error("invalid two-phase input accepted");
}

// Additional capillary guards have their own inventory; the original forty-six
// rejection cases and six sigma-zero accepted steps remain unchanged.
template<class F>void capillaryRejectsAction(F action,const char* expected)
{
    try{action();}catch(const std::exception& e){check(std::string(e.what()).find(expected)!=std::string::npos,
        std::string("capillary rejection reached a different stage: ")+e.what());++capillaryRejects;return;}
    throw std::runtime_error("invalid capillary step/getter accepted");
}
long double actualInterfaceMeasure(const HydroPhaseGeometry& phase)
{
    const auto& publication=phase.publication();long double area=0;
    for(std::size_t c=0;c<publication.topology().curves.size();++c)
        for(std::size_t g=0;g<publication.topology().curves[c].graphs().size();++g)
            for(const auto& q:publication.quadrature(c,g)){
                check(std::isfinite(q.physicalWeight)&&q.physicalWeight>0,"capillary area has nonpositive actual physical weight");
                check(q.sample.geometryIdentity==phase.geometryIdentity()&&sameClock(q.sample.timeLayer,phase.timeLayer()),
                    "capillary area sample lost publication owner/full clock");area+=q.physicalWeight;}
    check(std::isfinite(area)&&area>0,"capillary step has no actual interface measure");return area;
}
void nativeConstraints(const QuadtreeMesh& mesh,const HydroNativeValues& values)
{
    const auto& l=values.layout();const auto& u=values.normalVelocity();const auto& p=values.pressure();
    std::vector<long double> outward(l.pressureCount()),gross(l.pressureCount());long double mean=0,meanScale=0;
    for(int f=0;f<l.velocityCount();++f){const auto face=l.domain().face(mesh,l.nativeFace(f));const long double flow=stepFaceMeasure(mesh,l,f)*u[f];
        if(face.negativeCell>=0){const int cell=l.pressureDof(face.negativeCell);outward[cell]+=flow;gross[cell]+=std::abs(flow);}
        if(face.positiveCell>=0){const int cell=l.pressureDof(face.positiveCell);outward[cell]-=flow;gross[cell]+=std::abs(flow);}
        if(face.negativeCell<0||face.positiveCell<0)check(u[f]==0,"impermeable normal wall is not exact zero");}
    for(int cell=0;cell<l.pressureCount();++cell){peakEquation=std::max(peakEquation,scaled(outward[cell],0,gross[cell],"independent physical cell continuity"));
        const auto& actual=mesh.cells(mesh.cellIndexToLabel[l.nativeCell(cell)]);
        const auto lo=mesh.verts(actual.cellToVert[0]).coord(),hi=mesh.verts(actual.cellToVert[2]).coord();
        const long double volume=stepBoxMeasure(lo,hi,l.domain().axisymmetric());
        mean+=volume*p[cell];meanScale+=volume*std::abs(p[cell]);}
    peakEquation=std::max(peakEquation,scaled(mean,0,meanScale,"independent volume pressure gauge"));
}

// Independent actual-source face integral; the normal mean still observes U.
// Split public support-piece closures and use positive GL8; do not consume
// production mean preparation or its observer. It is not a KEP P1/P2 truth.
// The original rho/mu cases, counts, tolerances and source checks remain.
struct StepNativeMeanReference
{
    std::array<long double,2> mean{};
    std::array<long double,2> gross{}; // Positive physical value-work scale, including cancellation.
};
StepNativeMeanReference stepNativeMeanReference(const QuadtreeMesh& mesh,const HydroSpatialAdvector& actual,int u)
{
    const auto& l=actual.momentumGeometry().layout();check(!l.domain().axisymmetric(),"step mean reference requires XY");
    const auto face=l.domain().face(mesh,l.nativeFace(u));const int normal=face.horizontal?1:0,tangent=1-normal;
    const auto& edge=face.horizontal?mesh.hEdges(face.edgeLabel):mesh.vEdges(face.edgeLabel);
    auto a=mesh.verts(edge.edgeToVert[0]).coord(),b=mesh.verts(edge.edgeToVert[1]).coord();
    const auto coordinate=[](Vector2D p,int c){return c==0?p.x:p.y;};
    if(coordinate(a,tangent)>coordinate(b,tangent))std::swap(a,b);
    const double lo=coordinate(a,tangent),hi=coordinate(b,tangent),plane=coordinate(a,normal);
    check(hi>lo && coordinate(b,normal)==plane,"step mean reference native face is not canonical");
    const int native=face.negativeCell>=0?face.negativeCell:face.positiveCell,p=l.pressureDof(native);
    check(p>=0,"step mean reference lacks actual pressure owner");std::vector<double> cuts{lo,hi};
    for(const auto& piece:actual.reconstruction().pieces())if(piece.pressure==p &&
        plane>=coordinate(piece.lower,normal) && plane<=coordinate(piece.upper,normal)){
        const double x=std::max(lo,coordinate(piece.lower,tangent)),y=std::min(hi,coordinate(piece.upper,tangent));
        if(x<y){cuts.push_back(x);cuts.push_back(y);}}
    std::sort(cuts.begin(),cuts.end());cuts.erase(std::unique(cuts.begin(),cuts.end()),cuts.end());
    constexpr std::array<long double,8> nodes{{-.96028985649753623168L,-.79666647741362673959L,-.52553240991632898582L,-.18343464249564980494L,
        .18343464249564980494L,.52553240991632898582L,.79666647741362673959L,.96028985649753623168L}};
    constexpr std::array<long double,8> weights{{.10122853629037625915L,.22238103445337447054L,.31370664587788728734L,.36268378337836198297L,
        .36268378337836198297L,.31370664587788728734L,.22238103445337447054L,.10122853629037625915L}};
    StepNativeMeanReference result;long double measured=0;
    for(std::size_t k=1;k<cuts.size();++k){const long double x=cuts[k-1],y=cuts[k],half=(y-x)/2,center=(x+y)/2;
        check(half>0,"empty independent step source piece");
        for(std::size_t q=0;q<nodes.size();++q){const double t=static_cast<double>(center+half*nodes[q]);
            const auto point=normal==0?Vector2D{plane,t}:Vector2D{t,plane};const auto v=actual.cellVelocityJet(p,point).jet.value;
            const long double w=half*weights[q];check(w>0,"independent step mean has nonpositive physical weight");measured+=w;
            result.mean[0]+=w*v.x;result.mean[1]+=w*v.y;
            result.gross[0]+=w*std::abs(v.x);result.gross[1]+=w*std::abs(v.y);}}
    const long double length=static_cast<long double>(hi)-lo;
    scaled(measured,length,length,"independent step native line capacity",1e-12);
    scaled(actual.momentumGeometry().metrics().faceMeasures()[u],length,length,"independent step true XY face measure",1e-12);
    for(auto& value:result.mean)value/=length;for(auto& value:result.gross)value/=length;return result;
}

// Independent unit-basis reconstruction in every original moving rho/mu
// contract. The physical pair truth is assembled ONLY from unique F and its
// original negative/positive dual incidence; reaction rows do not enlarge
// the FREE error scale. No old QP P1/P2 moment defines this operator.
void inspectNativeConvection(const QuadtreeMesh& mesh,const HydroTwoPhaseStep& step)
{
    const auto& op=step.convection();const auto& m=op.material();const auto& l=m.layout();const int n=l.velocityCount();
    std::vector<bool> fixed(n);for(int i:op.fixedVelocityDofs())fixed[i]=true;
    const auto entries=static_cast<std::size_t>(n)*n;
    std::vector<long double> gross(n),k(entries),expected(entries),pairGross(entries);
    for(int i=0;i<n;++i)for(int side=0;side<4;++side)for(int g:m.geometry().volumes()[i].sides[side])
        gross[i]+=std::abs(m.faceRates()[g]);
    int physicalSegments=0,physicalPairs=0;
    const auto& segments=m.geometry().segments();
    check(m.faceRates().size()==segments.size(),"physical pair truth lost unique F inventory");
    for(std::size_t g=0;g<segments.size();++g){const auto& segment=segments[g];
        const int i=segment.negativeVelocity,j=segment.positiveVelocity;
        if(i<0||j<0){check(m.faceRates()[g]==0,"impermeable physical pair has a nonzero exterior F");continue;}
        ++physicalSegments;check(i!=j,"physical dual segment has identical incident native rows");
        const auto ij=static_cast<std::size_t>(i)*n+j,ji=static_cast<std::size_t>(j)*n+i;
        const long double coefficient=.5L*m.faceRates()[g],magnitude=.5L*std::abs(m.faceRates()[g]);
        expected[ij]+=coefficient;expected[ji]-=coefficient;pairGross[ij]+=magnitude;pairGross[ji]+=magnitude;
    }
    for(int i=0;i<n;++i)for(int j=i+1;j<n;++j)if(expected[static_cast<std::size_t>(i)*n+j]!=0)++physicalPairs;
    for(int j=0;j<n;++j){std::vector<double> u(n);u[j]=1;
        const auto value=op.apply(mesh,HydroNativeValues(l,u,std::vector<double>(l.pressureCount())),step.convectionResult().boundary(),m.sources(),m.timeLayer());
        for(int i=0;i<n;++i)k[static_cast<std::size_t>(i)*n+j]=value.momentumRates()[i]-value.dissipationRates()[i]-(i==j?.5L*m.cells()[i].outwardIntegral/m.deltaTime():0);}
    long double peak=0,reactionNorm=0,reactionConstantRaw=0,reactionConstantScale=0,coefficientRatio=0;
    for(int i=0;i<n;++i){long double rowSum=0,norm=0;
        for(int j=0;j<n;++j){const auto ij=static_cast<std::size_t>(i)*n+j;
            const auto a=k[ij],b=k[static_cast<std::size_t>(j)*n+i];
            check(std::isfinite(a),"nonfinite actual free/reaction coefficient");
            peak=std::max(peak,scaled(a,expected[ij],pairGross[ij]+std::abs(a),"actual unique-F physical pair coefficient",1e-12));
            peak=std::max(peak,scaled(a,-b,std::abs(a)+std::abs(b),"actual rho-dependent native skew",1e-12));
            rowSum+=a;norm+=std::abs(a);}
        const auto outward=m.cells()[i].outwardIntegral/m.deltaTime();
        peak=std::max(peak,scaled(rowSum,.5L*outward,norm+std::abs(outward),"actual rho-dependent canonical constant moment",1e-12));
        if(fixed[i]){check(step.solvedMidpointValues().normalVelocity()[i]==0,"reaction row acquired moving velocity");
            reactionNorm=std::max(reactionNorm,norm);reactionConstantRaw=std::max(reactionConstantRaw,std::abs(rowSum-.5L*outward));
            reactionConstantScale=std::max(reactionConstantScale,norm+std::abs(outward));continue;}
        check(norm<=2*gross[i]+1e-12*gross[i],"actual rho-dependent free complete-row L1 bound");
        if(gross[i]>0)coefficientRatio=std::max(coefficientRatio,norm/(2*gross[i]));
        const auto sourceMean=stepNativeMeanReference(mesh,m.material().midpointRate().spatialAdvector(),i);
        const int normal=l.normalAxis(i)==HydroDofLayout::NormalAxis::First?0:1;
        const long double sourceNative=m.material().midpointRate().spatialAdvector().values().normalVelocity()[i];
        peak=std::max(peak,scaled(sourceMean.mean[normal],sourceNative,sourceMean.gross[normal]+std::abs(sourceNative),"independent actual-source native mean",1e-12));
    }
    const auto& report=op.momentReport();
    check(report.physicalPairing==1&&report.physicalSegments==physicalSegments&&report.physicalPairs==physicalPairs,
        "physical pair report is not the unique-F incidence inventory");
    check(report.constraints==0&&report.iterations==0&&report.projectionCalls==0&&report.meanP1Rows==0&&report.meanP2Rows==0
        &&report.affineError==0&&report.quadraticError==0&&report.quadraticBudgetRatio==0&&report.quadraticBoundError==0
        &&report.maxQuadraticTargetRatio==0,"physical pair operator retained an old QP P1/P2 certificate");
    check(report.boundedRows==n-static_cast<int>(op.fixedVelocityDofs().size())&&report.freeRows==report.boundedRows
        &&report.fixedRows==static_cast<int>(op.fixedVelocityDofs().size()),"reaction rows counted as free stability rows");
    peak=std::max(peak,scaled(report.coefficientRatio,coefficientRatio,1,"actual free coefficient bound ratio",1e-12));
    peak=std::max(peak,scaled(report.maxReactionCoefficientL1,reactionNorm,reactionNorm,"actual reaction coefficient norm",1e-12));
    peak=std::max(peak,scaled(report.maxReactionConstantRawError,reactionConstantRaw,reactionConstantScale,"actual reaction constant raw defect",1e-12));
    std::cout<<"REACTION_NATIVE,"<<m.material().liquidDensity()<<','<<report.coefficientRatio<<','<<reactionNorm<<','<<report.maxReactionConstantRawError<<','<<peak<<'\n';
}

// The RZ short witness checks the actual final action against unique F/2
// incidence. It needs no axis-incompatible unit basis or new operator matrix.
void inspectRzConvection(const QuadtreeMesh& mesh,const HydroTwoPhaseStep& step)
{
    const auto& dual=step.dualMaterial();const auto& l=dual.layout();
    const auto& u=step.solvedMidpointValues().normalVelocity();
    const auto& result=step.reassessedConvectionResult();
    std::vector<long double> expected(u.size()),gross(u.size());
    long double power=0,expectedPower=0,powerScale=0;
    for(std::size_t i=0;i<u.size();++i){const auto value=.5L*dual.cells()[i].outwardIntegral/dual.deltaTime()*u[i];
        expected[i]=value;gross[i]=std::abs(value);expectedPower+=u[i]*value;}
    for(std::size_t g=0;g<dual.geometry().segments().size();++g){const auto& face=dual.geometry().segments()[g];
        const int i=face.negativeVelocity,j=face.positiveVelocity;const auto f=dual.faceRates()[g];
        if(i<0||j<0){check(f==0,"RZ closed exterior dual face has nonzero unique F");continue;}
        const long double left=.5L*f*u[j],right=.5L*f*u[i];
        expected[i]+=left;expected[j]-=right;gross[i]+=std::abs(left);gross[j]+=std::abs(right);
    }
    for(std::size_t i=0;i<u.size();++i){const auto actual=result.momentumRates()[i]-result.dissipationRates()[i];
        peakEquation=std::max(peakEquation,scaled(actual,expected[i],gross[i]+std::abs(actual),"RZ actual unique-F paired momentum",1e-12));
        power+=u[i]*actual;powerScale+=std::abs(u[i])*(gross[i]+std::abs(actual));}
    peakEnergy=std::max(peakEnergy,scaled(power,expectedPower,powerScale,"RZ actual paired KEP power",1e-12));
    for(const auto& component:result.report().components)
        check(component.kineticIdentity<=1e-12&&component.momentumIdentity<=1e-12
            &&component.boundaryKineticPower==0,"RZ actual KEP identity or closed-boundary work failed");
    const auto& report=step.convection().momentReport();
    check(report.physicalPairing==1&&report.constraints==0&&report.projectionCalls==0
        &&report.freeRows==static_cast<int>(step.saddle().freeVelocityDofs().size()),"RZ did not consume the physical native pairing");
    for(int axis:step.beginState().phaseGeometry().metrics().axisVelocityDofs()){
        check(u[axis]==0&&stepFaceMeasure(mesh,l,axis)==0
            &&!std::binary_search(step.saddle().freeVelocityDofs().begin(),step.saddle().freeVelocityDofs().end(),axis),
            "RZ zero-area radial axis was not eliminated");}
}

// The force color is a read-only view of the actual material midpoint. Check
// its native C action from physical incidence, independently of C's matrix.
void inspectCsfOwner(const QuadtreeMesh& mesh,const HydroTwoPhaseStep& step)
{
    const auto& load=step.capillaryLoad();const auto& midpoint=step.materialTransport().midpointRate().geometry();
    const auto& source=step.materialTransport().midpointRate().spatialAdvector();
    const auto& phase=load.phaseGeometry();const auto& l=midpoint.layout();
    const bool projected=load.discretization()==HydroCsfDiscretization::GeometricProjectedAdjoint;
    if(projected)load.validateFor(mesh,step.saddle().pressureLoad(),midpoint.geometryIdentity(),midpoint.timeLayer(),
        source,load.surfaceTension()!=0?&load.mixedSensitivity():nullptr);
    else load.validateFor(mesh,step.saddle().pressureLoad(),midpoint.geometryIdentity(),midpoint.timeLayer(),
        source.reconstruction(),source.boundary());
    check(projected||load.discretization()==HydroCsfDiscretization::GeometricAdjoint,
        "production step did not use a retained motion adjoint CSF discretization");
    if(projected)check(load.projectedMotionAdvector().identity()==source.identity(),
        "projected CSF lost the actual directional motion owner");
    check(load.isCapillary()&&load.surfaceTension()==step.report().surfaceTension
        &&load.layout().sameIdentity(l)&&load.reconstruction().sameIdentity(source.reconstruction())
        &&&load.pressureLoad().matrix()==&step.saddle().pressureLoad().matrix(),
        "CSF step lost the retained pressure/reconstruction owner");
    check(&phase.publication().state().geometry().snapshot()==&midpoint.publication().state().geometry().snapshot()
        &&&phase.cellVolumes()==&midpoint.cellVolumes(),
        "CSF replaced actual midpoint publication or material geometry");
    check(phase.geometryIdentity()==midpoint.geometryIdentity()&&sameClock(phase.timeLayer(),midpoint.timeLayer()),
        "CSF phase fraction lost actual midpoint identity/full clock");
    const auto& c=phase.liquidFractions();const auto& g=load.colorGradient();const auto& f=load.momentumRhsLoad();
    const auto& k0=load.referenceCurvatures();const auto& residual=load.adjointResidualLoad();
    const auto& reference=load.referencePressureLoad();
    check(c.size()==static_cast<std::size_t>(l.pressureCount())&&k0.size()==c.size()
        &&g.size()==static_cast<std::size_t>(l.velocityCount())&&f.size()==g.size()
        &&reference.size()==g.size()&&residual.size()==g.size(),"CSF native component/load inventory");
    for(double value:c)check(std::isfinite(value)&&value>=0&&value<=1,"CSF cell fraction is unbounded");
    for(int u=0;u<l.velocityCount();++u){const auto face=l.domain().face(mesh,l.nativeFace(u));
        const int pm=face.negativeCell>=0?l.pressureDof(face.negativeCell):-1;
        const int pp=face.positiveCell>=0?l.pressureDof(face.positiveCell):-1;
        const long double cm=pm>=0?c[pm]:0,cp=pp>=0?c[pp]:0;
        const long double area=stepFaceMeasure(mesh,l,u),expected=area*(cp-cm);
        scaled(g[u],expected,std::abs(area*cp)+std::abs(area*cm),"CSF actual native C alpha action",1e-12);
        const long double km=pm>=0&&std::isfinite(k0[pm])?k0[pm]:0;
        const long double kp=pp>=0&&std::isfinite(k0[pp])?k0[pp]:0;
        const long double pressurePart=load.surfaceTension()*area*(kp*cp-km*cm);
        scaled(reference[u],pressurePart,std::abs(load.surfaceTension()*area*kp*cp)
            +std::abs(load.surfaceTension()*area*km*cm),"CSF component reference C action",1e-12);
        check(std::isfinite(f[u])&&std::isfinite(residual[u]),"CSF adjoint load is nonfinite");
        scaled(f[u],pressurePart+residual[u],std::abs(pressurePart)+std::abs(residual[u]),
            "CSF actual pressure plus unsmoothed traction RHS",1e-12);
    }
}

// Sharp uses the same actual midpoint publication and reconstruction owner.
void inspectSharpOwner(const QuadtreeMesh& mesh,const HydroTwoPhaseStep& step)
{
    const auto& load=step.sharpNormalStressLoad();const auto& midpoint=step.materialTransport().midpointRate().geometry();
    const auto& source=step.materialTransport().midpointRate().spatialAdvector();const auto& l=midpoint.layout();
    const auto& physical=step.capillaryMomentumRhsLoad();
    load.validateFor(mesh,step.saddle().pressureLoad(),midpoint.geometryIdentity(),midpoint.timeLayer(),
        source.reconstruction(),source.boundary());
    check(&load.geometry().cellVolumes()==&midpoint.cellVolumes()
        &&&load.geometry().publication().state().geometry().snapshot()==&midpoint.publication().state().geometry().snapshot(),
        "sharp replaced actual midpoint material/publication owner");
    check(load.reconstruction().sameIdentity(source.reconstruction())&&load.surfaceTension()==step.report().surfaceTension,
        "sharp lost actual midpoint source/sigma");
    const auto& chi=load.phaseIndicators();const auto& scalar=midpoint.publication().state().geometry().snapshot();
    for(int p=0;p<l.pressureCount();++p){const int native=l.nativeCell(p);
        const double phi=scalar.scalarJetInNativeCell(native,mesh.cells(mesh.cellIndexToLabel[native]).centerCoord()).value;
        check(chi[p]==(phi<0?1:phi>0?0:.5),"sharp pressure-center phase is not same-scalar");
    }
    for(int u=0;u<l.velocityCount();++u){const auto face=l.domain().face(mesh,l.nativeFace(u));
        const bool internal=face.negativeCell>=0&&face.positiveCell>=0&&stepFaceMeasure(mesh,l,u)>0;
        const long double jump=internal?chi[l.pressureDof(face.positiveCell)]-chi[l.pressureDof(face.negativeCell)]:0;
        const long double expected=static_cast<long double>(load.surfaceTension())*load.linkCurvatures()[u]*stepFaceMeasure(mesh,l,u)*jump;
        scaled(physical[u],expected,std::abs(expected),"sharp same-native-incidence physical RHS",1e-12);
        check(physical[u]==-load.balanceLoads()[u],"sharp physical RHS did not negate balance exactly once");
        if(!internal||jump==0)check(physical[u]==0&&load.linkCurvatures()[u]==0&&load.linkLengths()[u]==0,
            "sharp loaded a same-phase or outer-wall link");
    }
}

void inspectStep(const QuadtreeMesh& mesh,const HydroTwoPhaseStep& step,bool moving,double velocityUpdateTolerance=1e-11,double geometryUpdateTolerance=1e-11,double powerBalanceTolerance=1e-11)
{
    if(step.beginState().phaseGeometry().layout().domain().axisymmetric())inspectRzConvection(mesh,step);
    else inspectNativeConvection(mesh,step);
    const bool linear=step.report().timeScheme!=HydroTwoPhaseTimeScheme::FullyCoupledMidpoint;
    const auto& material=step.materialTransport();const auto& dual=step.dualMaterial();const auto& l=dual.layout();
    const auto& m0=dual.beginMass().masses();const auto& m1=dual.endMass().masses();
    const auto& body=step.bodyLoads();check(body.size()==static_cast<std::size_t>(l.velocityCount()),"body native load inventory");
    const auto& u0=step.beginState().values().normalVelocity();const auto& u1=step.candidateState().values().normalVelocity();
    const auto& guessed=material.midpointRate().spatialAdvector().values().normalVelocity();
    const auto& bar=step.solvedMidpointValues().normalVelocity();const auto& pressure=step.candidateState().values().pressure();
    const auto& convection=step.convectionResult();const auto& reassessed=step.reassessedConvectionResult();
    const Matrix k=step.viscosity().viscosity().matrix();const auto ku0=matvec(k,u0),kbar=matvec(k,bar);
    Matrix h(l.velocityCount(),l.velocityCount());
    if(linear)h=step.convection().frozenMomentumOperator(mesh,dual.sources(),dual.timeLayer());
    const auto hu0=matvec(h,u0);
    const double dt=step.report().deltaTime;const auto sources=material.sources();
    check(material.isSpaceTime()&&material.isSpatial(),"step did not consume actual spacetime material transport");
    for(const auto* phase:std::array<const HydroPhaseGeometry*,3>{{&material.beginGeometry(),&material.midpointRate().geometry(),&material.endGeometry()}}){
        const auto& primary=phase->publication().state().primary();
        check(primary.transport==HermiteGalsTransport::FirstJet&&primary.phaseVolumes.empty(),"candidate family is not PureFirstJet without inventory");}
    check(sources.begin==step.beginState().geometryIdentity()&&sources.midpoint==dual.midpointMass().geometryIdentity()
        &&sources.end==step.candidateState().geometryIdentity(),"candidate three geometry identities are not bound");
    check(sources.begin!=sources.midpoint&&sources.midpoint!=sources.end&&sources.begin!=sources.end,"candidate clock family collapsed identities");
    check(material.midpointFlux().advectorIdentity()==material.midpointRate().spatialAdvector().identity(),"unique actual material source changed");
    check(&dual.material().faceIntegrals()==&material.faceIntegrals(),"dual material copied or reintegrated primal F");
    check(&convection.material().faceRates()==&dual.faceRates() && &reassessed.material().faceRates()==&dual.faceRates(),"KEP replaced unique dual F");
    check(step.convection().isImpermeable()&&step.convection().periodicPairs().empty()
        &&step.convection().fixedVelocityDofs()==step.beginState().boundary().velocityDofs(),
        "owner did not consume the typed impermeable native convection path");
    std::vector<int> free;for(int i=0;i<l.velocityCount();++i)
        if(!std::binary_search(step.convection().fixedVelocityDofs().begin(),step.convection().fixedVelocityDofs().end(),i))free.push_back(i);
    check(free==step.saddle().freeVelocityDofs(),"closed convection/saddle normal constraints differ");
    check(step.beginState().timeLayer().time()==material.beginGeometry().timeLayer().time()
        &&step.candidateState().timeLayer().time()==material.endGeometry().timeLayer().time(),"material interval clocks differ from state clocks");
    const std::array<const HydroMomentumMass*,3> masses{{&dual.beginMass(),&dual.midpointMass(),&dual.endMass()}};
    for(const auto* mass:masses){
        check(mass->representation()==step.report().momentumDensityRepresentation,"three masses differ from selected density representation");
        for(int i=0;i<l.velocityCount();++i){
            const auto phase=mass->representedPhaseVolumes()[i];const auto& box=dual.geometry().volumes()[i];
            const bool rz=l.domain().axisymmetric();
            const long double volume=stepBoxMeasure(box.lower,box.upper,rz);
            check(phase.liquid>=0&&phase.gas>=0,"three-layer represented phase moment is negative");
            peakMaterial=std::max(peakMaterial,scaled(phase.liquid+phase.gas,volume,volume,"three-layer physical dual rectangle measure",1e-12));
            long double expectedLiquid=phase.liquid,expectedGas=phase.gas;
            if(mass->representation()==HydroMomentumMass::DensityRepresentation::SharpGeometry){
                check(&mass->representedPhaseVolumes()==&mass->phaseVolumes(),"sharp represented moments replaced original actual phase geometry");
            }else{
                expectedLiquid=expectedGas=0;
                const auto& actual=mass->phaseGeometry();
                for(const auto& overlap:box.overlaps){
                    const auto cell=actual.publication().state().geometry().snapshot().nativeCellBounds(l.nativeCell(overlap.pressure));
                    const long double cellVolume=stepBoxMeasure(cell[0],cell[1],rz);
                    const long double weight=stepBoxMeasure(overlap.lower,overlap.upper,rz)/cellVolume;
                    expectedLiquid+=weight*actual.cellVolumes()[overlap.pressure].liquid;
                    expectedGas+=weight*actual.cellVolumes()[overlap.pressure].gas;
                }
                scaled(phase.liquid,expectedLiquid,volume,"independent projected dual liquid moment",1e-12);
                scaled(phase.gas,expectedGas,volume,"independent projected dual gas moment",1e-12);
            }
            const long double expected=material.liquidDensity()*expectedLiquid+material.gasDensity()*expectedGas;
            peakMaterial=std::max(peakMaterial,scaled(mass->masses()[i],expected,expected,"three-layer selected density mass",1e-12));
        }
    }
    const auto& report=step.report();
    const std::vector<double>* capillary=step.hasCapillaryLoad()?&step.capillaryMomentumRhsLoad():nullptr;
    check((capillary!=nullptr)==(report.surfaceTension>0),"capillary report/active owner mismatch");
    if(capillary){if(step.capillaryModel()==HydroCapillaryModel::BalancedCsf)inspectCsfOwner(mesh,step);
        else if(step.capillaryModel()==HydroCapillaryModel::SharpNormalStress)inspectSharpOwner(mesh,step);
        }
    else check(report.surfaceTension==0&&report.capillaryPower==0&&report.interfaceMeasureBegin==0
        &&report.interfaceMeasureMidpoint==0&&report.interfaceMeasureEnd==0&&report.surfaceEnergyBegin==0
        &&report.surfaceEnergyEnd==0&&report.surfaceEnergyChange==0&&report.capillaryEnergyDefect==0,
        "sigma-zero step did not retain exact unsampled capillary fast path");
    check(report.history.size()==static_cast<std::size_t>(report.iterations),"incorrect retained stage history");
    if(!capillary)for(const auto& iteration:report.history)
        check(iteration.csfColorSeconds==0&&iteration.csfProjectionSeconds==0,"sigma-zero step sampled CSF factories");
    if(!linear){
        check(report.iterations>=2&&report.iterations<=12,"incorrect corrector history");
        check(report.velocityUpdateRelative<=velocityUpdateTolerance&&report.geometryUpdateRelative<=geometryUpdateTolerance
            &&report.advectorGuessSolveRelative<=velocityUpdateTolerance&&report.frozenOwnerMomentumDefect<=1e-11,"returned unaccepted corrector");
    }else{
        check(report.iterations==(report.corrected?2:1)&&report.frozenOwnerMomentumDefect<=1e-11
            &&report.linearMomentumDefect<=1e-11,"linear midpoint changed its stage budget or momentum gate");
        check(std::isfinite(report.interfaceMismatch)&&report.interfaceMismatch>=0
            &&std::isfinite(report.finalInterfaceMismatch)&&report.finalInterfaceMismatch>=0
            &&std::isfinite(report.correctionRelative)&&report.correctionRelative>=0,"invalid actual linear stage diagnostics");
    }
    const auto& saddle=step.solveResult().report();
    check(step.materialTransport().spaceTimeReport().firstJetTimeProjection==HermiteGalsTimeProjection::CharacteristicPullback,
        "step did not use the declared characteristic pullback time projection");
    check(saddle.algebraicRelativeResidual<=1e-11&&saddle.freeMomentumDefect<=1e-11&&saddle.maxCellContinuityDefect<=1e-11
        &&saddle.maxPressureMeanDefect<=1e-11&&saddle.maxBoundaryFluxDefect<=1e-12,"original saddle gate changed");
    nativeConstraints(mesh,step.candidateState().values());

    // A and RHS use the retained assembly guess; the independent physical
    // equation below uses the FINAL midpoint and reassessed KEP output.
    Matrix expectedA=.5*k;if(linear)expectedA+=.5*h;
    for(int i=0;i<l.velocityCount();++i)expectedA.coeffRef(i,i)+=m1[i]/dt;
    const auto au1=matvec(expectedA,u1);
    const Matrix da=step.velocityMatrix()-expectedA;
    peakEquation=std::max(peakEquation,scaled(da.norm(),0,expectedA.norm(),"actual variable-mass velocity matrix",1e-12));
    long double rhsGap=0,rhsScale=0,equationGap=0,equationScale=0,ke0=0,ke1=0,inertial=0,cross=0;
    long double pressurePower=0,viscousPower=0,convectionPower=0,consistencyPower=0,dissipationPower=0,nonzeroR=0,massMovement=0;
    long double temporalScale=0,powerScale=0,capillaryPower=0,capillaryWork=0,bodyPower=0,bodyWork=0;
    long double physicalPowerBalance=0,physicalPowerWork=0,linearGap=0,linearWork=0;
    for(int i=0;i<l.velocityCount();++i){
        scaled(bar[i],(.5L*u0[i]+.5L*u1[i]),std::abs(.5L*u0[i])+std::abs(.5L*u1[i]),"actual solved midpoint",1e-12);
        const long double force=capillary?(*capillary)[i]:0;
        const long double a=m0[i]*u0[i]/dt,b=-.5L*ku0[i],
            c=linear?-.5L*hu0[i]:convection.momentumLoads()[i],rhs=a+b+c+force+body[i];
        rhsGap=std::hypot(rhsGap,step.integralMomentumRhs()[i]-rhs);rhsScale=std::hypot(rhsScale,std::abs(a)+std::abs(b)+std::abs(c)+std::abs(force)+std::abs(body[i]));
        const auto begin=dual.beginMass().representedPhaseVolumes()[i],end=dual.endMass().representedPhaseVolumes()[i];
        const long double dm=material.liquidDensity()*(end.liquid-begin.liquid)+material.gasDensity()*(end.gas-begin.gas);
        long double outward=0,gross=0;
        for(int side=0;side<4;++side)for(int segment:dual.geometry().volumes()[i].sides[side]){
            const auto f=dual.faceRates()[segment];outward+=(side%2?1:-1)*f;gross+=std::abs(f);}
        const long double mdot=dm/dt,r=mdot+outward,rs=std::abs(mdot)+gross;
        peakMaterial=std::max(peakMaterial,scaled(dual.cells()[i].defectRate,r,rs,"actual R from selected represented mass and original F",1e-12));
        const long double source=bar[i]*r,load=-reassessed.momentumRates()[i]+source;
        peakMaterial=std::max(peakMaterial,scaled(reassessed.momentumLoads()[i],load,std::abs(reassessed.momentumRates()[i])+std::abs(bar[i])*rs,
            "one final uR load",1e-12));
        peakMaterial=std::max(peakMaterial,scaled(convection.momentumLoads()[i],-convection.momentumRates()[i]+guessed[i]*r,
            std::abs(convection.momentumRates()[i])+std::abs(guessed[i])*rs,"one original guess uR load",1e-12));
        nonzeroR+=std::abs(r);massMovement+=std::abs(m1[i]-m0[i]);
        const auto face=l.domain().face(mesh,l.nativeFace(i));
        const long double pm=face.negativeCell>=0?pressure[l.pressureDof(face.negativeCell)]:0;
        const long double pp=face.positiveCell>=0?pressure[l.pressureDof(face.positiveCell)]:0;
        const long double cp=stepFaceMeasure(mesh,l,i)*(pp-pm),impulse=m1[i]*u1[i]-m0[i]*u0[i],rate=impulse/dt;
        if(face.negativeCell>=0&&face.positiveCell>=0){equationGap=std::hypot(equationGap,rate+kbar[i]-load+cp-force-body[i]);
            equationScale=std::hypot(equationScale,(std::abs(m1[i]*u1[i])+std::abs(m0[i]*u0[i]))/dt+std::abs(kbar[i])+std::abs(load)+std::abs(cp)+std::abs(force)+std::abs(body[i]));}
        if(linear&&face.negativeCell>=0&&face.positiveCell>=0){
            linearGap=std::hypot(linearGap,static_cast<long double>(au1[i])+cp-rhs);
            linearWork=std::hypot(linearWork,std::abs(au1[i])+std::abs(cp)+std::abs(rhs));}
        ke0+=.5L*m0[i]*u0[i]*u0[i];ke1+=.5L*m1[i]*u1[i]*u1[i];inertial+=bar[i]*impulse;
        cross+=.5L*u0[i]*u1[i]*(m1[i]-m0[i]);
        temporalScale+=std::abs(bar[i]*impulse)+.5L*(std::abs(m0[i]*u0[i]*u0[i])+std::abs(m1[i]*u1[i]*u1[i]));
        pressurePower+=bar[i]*cp;viscousPower+=bar[i]*kbar[i];
        capillaryPower+=bar[i]*force;capillaryWork+=std::abs(bar[i]*force);
        bodyPower+=bar[i]*body[i];bodyWork+=std::abs(bar[i]*body[i]);
        convectionPower+=bar[i]*reassessed.momentumRates()[i];consistencyPower+=bar[i]*bar[i]*r;
        dissipationPower+=bar[i]*reassessed.dissipationRates()[i];
        powerScale+=std::abs(bar[i]*rate)+std::abs(bar[i]*kbar[i])+std::abs(bar[i]*load)+std::abs(bar[i]*cp)+std::abs(bar[i]*force)+std::abs(bar[i]*body[i]);
        // Independent physical power residual: true endpoint masses, actual
        // final loads and native-face pressure, before any report aggregate.
        const long double physicalRate=(static_cast<long double>(m1[i])*u1[i]-static_cast<long double>(m0[i])*u0[i])/dt;
        const long double velocity=bar[i],actualLoad=reassessed.momentumLoads()[i];
        physicalPowerBalance+=velocity*(physicalRate+kbar[i]+cp-actualLoad-force-body[i]);
        physicalPowerWork+=std::abs(velocity*physicalRate)+std::abs(velocity*kbar[i])
            +std::abs(velocity*cp)+std::abs(velocity*actualLoad)+std::abs(velocity*force)+std::abs(velocity*body[i]);
    }
    peakEquation=std::max(peakEquation,scaled(rhsGap,0,rhsScale,"actual variable-mass momentum RHS",1e-12));
    if(linear)peakEquation=std::max(peakEquation,scaled(linearGap,0,linearWork,"actual implicit frozen-H free momentum"));
    peakEquation=std::max(peakEquation,scaled(equationGap,0,equationScale,"final nonlinear physical free momentum"));
    peakEnergy=std::max(peakEnergy,scaled(inertial-(ke1-ke0),cross,temporalScale,"independent nonzero mass time cross term",1e-12));
    peakEnergy=std::max(peakEnergy,scaled(inertial/dt,-convectionPower+consistencyPower-viscousPower-pressurePower+capillaryPower+bodyPower,
        powerScale,"final physical momentum power",powerBalanceTolerance));
    const auto physicalRelative=scaled(physicalPowerBalance,0,physicalPowerWork,
        "independent final physical momentum power residual",powerBalanceTolerance);
    // Compare saved defect on the original physical work scale. This can
    // detect a fabricated zero without dividing by a tiny residual itself.
    scaled(static_cast<long double>(report.momentumPowerRelativeDefect)*physicalPowerWork,
        std::abs(physicalPowerBalance),physicalPowerWork,"saved actual physical momentum power defect",1e-12);
    check(report.momentumPowerRelativeDefect<=powerBalanceTolerance,"reported physical momentum power exceeds its explicit target");
    std::cout<<"PHYSICAL_MOMENTUM_POWER_TARGET,"<<physicalRelative<<','<<report.momentumPowerRelativeDefect
        <<','<<powerBalanceTolerance<<'\n';
    peakEnergy=std::max(peakEnergy,scaled(report.kineticBegin,ke0,std::abs(ke0),"begin kinetic ledger",1e-12));
    peakEnergy=std::max(peakEnergy,scaled(report.kineticEnd,ke1,std::abs(ke1),"end kinetic ledger",1e-12));
    peakEnergy=std::max(peakEnergy,scaled(report.temporalMassCrossTerm,cross,temporalScale,"retained temporal cross ledger",1e-12));
    peakEnergy=std::max(peakEnergy,scaled(report.inertialImpulseWork,inertial,temporalScale,"retained inertial ledger",1e-12));
    peakEnergy=std::max(peakEnergy,scaled(report.pressurePower,pressurePower,powerScale,"pressure power ledger",1e-12));
    peakEnergy=std::max(peakEnergy,scaled(report.viscousPower,viscousPower,powerScale,"viscous power ledger",1e-12));
    peakEnergy=std::max(peakEnergy,scaled(report.consistencyPower,consistencyPower,powerScale,"actual uR power ledger",1e-12));
    peakEnergy=std::max(peakEnergy,scaled(report.numericalDissipationPower,dissipationPower,powerScale,"separate numerical dissipation ledger",1e-12));
    peakEnergy=std::max(peakEnergy,scaled(report.interfacePower,capillaryPower,capillaryWork,"actual selected total interface force power ledger",1e-12));
    peakEnergy=std::max(peakEnergy,scaled(report.capillaryPower,capillaryPower,capillaryWork,"actual pure capillary physical force power ledger",1e-12));
    peakEnergy=std::max(peakEnergy,scaled(report.bodyPower,bodyPower,bodyWork,"actual retained midpoint body power ledger",1e-12));
    if(capillary){const long double beginArea=actualInterfaceMeasure(material.beginGeometry()),
            midpointArea=actualInterfaceMeasure(material.midpointRate().geometry()),endArea=actualInterfaceMeasure(material.endGeometry());
        scaled(report.interfaceMeasureBegin,beginArea,beginArea,"actual begin interface measure",1e-12);
        scaled(report.interfaceMeasureMidpoint,midpointArea,midpointArea,"actual midpoint interface measure",1e-12);
        scaled(report.interfaceMeasureEnd,endArea,endArea,"actual end interface measure",1e-12);
        const long double beginEnergy=report.surfaceTension*beginArea,endEnergy=report.surfaceTension*endArea,change=endEnergy-beginEnergy;
        scaled(report.surfaceEnergyBegin,beginEnergy,beginEnergy,"actual begin surface energy",1e-12);
        scaled(report.surfaceEnergyEnd,endEnergy,endEnergy,"actual end surface energy",1e-12);
        scaled(report.surfaceEnergyChange,change,std::abs(beginEnergy)+std::abs(endEnergy),"actual surface energy change",1e-12);
        scaled(report.capillaryEnergyDefect,change+dt*capillaryPower,std::abs(change)+std::abs(dt*capillaryPower),
            "observed surface-energy/capillary-work defect",1e-12);
        std::cout<<"CAPILLARY_STEP_OBSERVATION,"<<report.surfaceTension<<','<<beginArea<<','<<midpointArea<<','<<endArea<<','<<change<<','<<capillaryPower<<','<<report.capillaryEnergyDefect<<'\n';}
    check(report.wallReactionPower==0,"stationary impermeable wall did work");
    if(step.viscosity().liquidViscosity()==0&&step.viscosity().gasViscosity()==0)
        check(k.nonZeros()==0&&viscousPower==0&&!step.viscosity().hasQuadrature()
            &&step.viscosity().viscosity().strain().points().empty()
            &&step.viscosity().viscosity().strain().matrix().rows()==0
            &&step.viscosity().viscosity().strain().matrix().cols()==l.velocityCount(),"mu0 did not use the exact unsampled zero form");
    else check(k.norm()>0&&viscousPower>0&&step.viscosity().hasQuadrature()
        &&!step.viscosity().viscosity().strain().points().empty(),"finite physical viscosity did not retain actual quadrature");
    if(moving){check(u1!=u0,"nonzero two-phase velocity was frozen");
        check(step.beginState().publication().state().closure().cellBubbleCoefficients!=step.candidateState().publication().state().closure().cellBubbleCoefficients,
            "actual moving FirstJet centers were frozen");
        long double phaseMovement=0;
        for(int p=0;p<l.pressureCount();++p)phaseMovement+=std::abs(material.endGeometry().cellVolumes()[p].liquid-material.beginGeometry().cellVolumes()[p].liquid);
        check(phaseMovement>0,"actual two-phase geometry did not move");
        if(material.liquidDensity()!=material.gasDensity())check(massMovement>0&&nonzeroR>0,"moving contrast material mass/R was frozen");}
    std::cout<<"TWO_PHASE_STEP,"<<material.liquidDensity()<<','<<step.viscosity().liquidViscosity()<<','<<step.candidateState().timeLayer().step
        <<','<<report.iterations<<','<<massMovement<<','<<nonzeroR<<','<<equationGap<<','<<cross<<','<<viscousPower<<'\n';
}

void compareRelaxedEndpoint(const HydroTwoPhaseStep& strict,const HydroTwoPhaseStep& relaxed)
{
    check(strict.beginState().sameIdentity(relaxed.beginState())
        &&sameClock(strict.candidateState().timeLayer(),relaxed.candidateState().timeLayer()),
        "strict/relaxed comparison changed the actual begin owner or complete endpoint clock");
    check(relaxed.report().iterations<strict.report().iterations,
        "relaxed stopping did not reduce actual Picard iterations: strict="+std::to_string(strict.report().iterations)
        +", relaxed="+std::to_string(relaxed.report().iterations)
        +", nonlinear="+std::to_string(relaxed.report().frozenOwnerMomentumDefect));
    // Each physical/native channel has its own work norm. No absolute floor,
    // increment-only scale, rounded reference field or fitted geometry is used.
    std::array<long double,7> gap{},strictNorm{},relaxedNorm{},beginNorm{};
    const std::array<const char*,7> names{{"endpoint U","endpoint p","vertex phi","vertex dx","vertex dy","edge phi","edge dn"}};
    const auto add=[&](std::size_t k,double a,double b,double initial){
        gap[k]=std::hypot(gap[k],static_cast<long double>(b)-a);
        strictNorm[k]=std::hypot(strictNorm[k],static_cast<long double>(a));
        relaxedNorm[k]=std::hypot(relaxedNorm[k],static_cast<long double>(b));
        beginNorm[k]=std::hypot(beginNorm[k],static_cast<long double>(initial));};
    const auto& a=strict.candidateState().values();const auto& b=relaxed.candidateState().values();const auto& initial=strict.beginState().values();
    check(a.layout().sameIdentity(b.layout())&&a.layout().sameIdentity(initial.layout()),"relaxed endpoint comparison changed native layout");
    for(std::size_t i=0;i<a.normalVelocity().size();++i)add(0,a.normalVelocity()[i],b.normalVelocity()[i],initial.normalVelocity()[i]);
    for(std::size_t i=0;i<a.pressure().size();++i)add(1,a.pressure()[i],b.pressure()[i],initial.pressure()[i]);
    const auto& ap=strict.candidateState().publication().state().primary();
    const auto& bp=relaxed.candidateState().publication().state().primary();
    const auto& ip=strict.beginState().publication().state().primary();
    check(ap.vertices.size()==bp.vertices.size()&&ap.vertices.size()==ip.vertices.size()
        &&ap.edges.size()==bp.edges.size()&&ap.edges.size()==ip.edges.size()
        &&ap.geometryCells==bp.geometryCells&&ap.geometryCells==ip.geometryCells,
        "relaxed endpoint changed authoritative FirstJet channel topology");
    for(std::size_t i=0;i<ap.vertices.size();++i){const auto& x=ap.vertices[i];const auto& y=bp.vertices[i];const auto& z=ip.vertices[i];
        check(x.vertex==y.vertex&&x.vertex==z.vertex,"relaxed endpoint changed native vertex identity");
        add(2,x.phi,y.phi,z.phi);add(3,x.dx,y.dx,z.dx);add(4,x.dy,y.dy,z.dy);}
    for(std::size_t i=0;i<ap.edges.size();++i){const auto& x=ap.edges[i];const auto& y=bp.edges[i];const auto& z=ip.edges[i];
        check(x.vertices==y.vertices&&x.vertices==z.vertices,"relaxed endpoint changed native edge identity");
        add(5,x.phi,y.phi,z.phi);add(6,x.dn,y.dn,z.dn);}
    for(std::size_t k=0;k<gap.size();++k){const auto work=std::max({strictNorm[k],relaxedNorm[k],beginNorm[k]});
        const auto relative=scaled(gap[k],0,work,names[k],1e-7);
        std::cout<<"RELAXED_ENDPOINT_CHANNEL,"<<names[k]<<','<<gap[k]<<','<<work<<','<<relative<<'\n';}
    ++relaxedComparisons;
    std::cout<<"RELAXED_PICARD_COMPARISON,"<<strict.report().iterations<<','<<relaxed.report().iterations
        <<','<<relaxed.report().velocityUpdateRelative<<','<<relaxed.report().geometryUpdateRelative
        <<','<<relaxed.report().frozenOwnerMomentumDefect<<",PASS\n";
}

// Actual Frozen inputs use the original raw A and RHS directly. The public
// S/W roundtrip remains an independent legacy capability, not this solver oracle.
namespace FrozenPrepareContract {
unsigned pressureReferences=0,checkpointReferences=0,laggedReferences=0,threadGroups=0,threadPairs=0,coldCheckpointReferences=0;
struct Floating {
    std::fenv_t native{};int error=errno;
#if defined(_M_X64) || defined(__SSE2__)
    unsigned csr=_mm_getcsr();
#endif
    Floating(){check(std::fegetenv(&native)==0,"Frozen reference FP capture");}
    void restore()const{check(std::fesetenv(&native)==0,"Frozen reference FP restore");
#if defined(_M_X64) || defined(__SSE2__)
        _mm_setcsr(csr);
#endif
        errno=error;}
    bool same(const Floating& b)const{return error==b.error&&std::memcmp(&native,&b.native,sizeof native)==0
#if defined(_M_X64) || defined(__SSE2__)
        &&csr==b.csr
#endif
        ;}
};
struct Restore {
    Floating original;
#ifdef _OPENMP
    int threads=omp_get_max_threads(),dynamic=omp_get_dynamic();
#endif
    ~Restore(){
#ifdef _OPENMP
        omp_set_num_threads(threads);omp_set_dynamic(dynamic);
#endif
        original.restore();}
};
template<class T>void bits(const std::vector<T>& a,const std::vector<T>& b,const char* message){
    check(a.size()==b.size()&&(a.empty()||std::memcmp(a.data(),b.data(),a.size()*sizeof(T))==0),message);
}
void matrix(const Matrix& a,const Matrix& b,const char* message){
    check(a.isCompressed()&&b.isCompressed()&&a.rows()==b.rows()&&a.cols()==b.cols()&&a.nonZeros()==b.nonZeros(),message);
    check(std::memcmp(a.outerIndexPtr(),b.outerIndexPtr(),(a.outerSize()+1)*sizeof(Matrix::StorageIndex))==0,message);
    if(a.nonZeros())check(std::memcmp(a.innerIndexPtr(),b.innerIndexPtr(),a.nonZeros()*sizeof(Matrix::StorageIndex))==0
        &&std::memcmp(a.valuePtr(),b.valuePtr(),a.nonZeros()*sizeof(double))==0,message);
}
void vector(const Eigen::VectorXd& a,const Eigen::VectorXd& b,const char* message){
    check(a.size()==b.size()&&(a.size()==0||std::memcmp(a.data(),b.data(),a.size()*sizeof(double))==0),message);
}
void payload(const HydroFrozenMidpoint& a,const HydroFrozenMidpoint& b){
    matrix(a.viscosityMatrix(),b.viscosityMatrix(),"Frozen own-request K words/CSR differ");
    matrix(a.physicalVelocityMatrix(),b.physicalVelocityMatrix(),"Frozen own-request physical A words/CSR differ");
    matrix(a.velocityMatrix(),b.velocityMatrix(),"Frozen own-request raw A words/CSR differ");
    matrix(a.frozenConvectionMatrix(),b.frozenConvectionMatrix(),"Frozen own-request H words/CSR differ");
    matrix(a.saddle().velocityBlock().matrix(),b.saddle().velocityBlock().matrix(),"Frozen own-request actual velocity owner words/CSR differ");
    vector(a.beginViscousAction(),b.beginViscousAction(),"Frozen own-request K*u0 words differ");
    vector(a.beginFrozenAction(),b.beginFrozenAction(),"Frozen own-request H*u0 words differ");
    bits(a.integralMomentumRhs(),b.integralMomentumRhs(),"Frozen own-request RHS words differ");
    bits(a.solveResult().values().normalVelocity(),b.solveResult().values().normalVelocity(),"Frozen own-request solved U words differ");
    bits(a.solveResult().values().pressure(),b.solveResult().values().pressure(),"Frozen own-request solved p words differ");
    check(std::isnan(a.saddle().velocityBlock().symmetryDefect())&&std::isnan(b.saddle().velocityBlock().symmetryDefect()),
        "Frozen production reported a requested combined symmetry diagnostic");
}
// Independent original mixed equations, using raw physical inputs rather than
// another invocation of the new private factory or its assembly/index plan.
void mixedOriginalEquation(const QuadtreeMesh& mesh,const HydroFrozenMidpoint& actual,
    const Matrix& raw,const std::vector<double>& rhs,const HydroSaddleConfig& config){
    const auto& saddle=actual.saddle();const auto& layout=actual.material().layout();
    const auto& free=saddle.freeVelocityDofs();const auto& c=actual.pressureLoad().matrix();
    const int n=layout.velocityCount(),np=layout.pressureCount(),nf=static_cast<int>(free.size()),ng=saddle.pressureGaugeCount();
    std::vector<int> row(n,-1);std::vector<bool> outlet(n,false);
    for(const auto& port:actual.boundary().pressureOutlets)outlet.at(port.velocityDof)=true;
    for(int i=0;i<nf;++i){check(free[i]>=0&&free[i]<n&&row[free[i]]<0,"raw mixed duplicate/invalid free face");row[free[i]]=i;}
    for(int u=0;u<n;++u){const auto f=layout.domain().face(mesh,layout.nativeFace(u));
        check((row[u]>=0)==(stepFaceMeasure(mesh,layout,u)>0&&(f.outwardSign==0||outlet[u])),
            "raw mixed free rows disagree with actual interior/outlet/axis geometry");}
    std::vector<Eigen::Triplet<double,Eigen::Index>> entries;
    for(int u:free){for(Matrix::InnerIterator it(raw,u);it;++it)if(row[it.col()]>=0)
            entries.emplace_back(row[u],row[it.col()],it.value());
        for(Matrix::InnerIterator it(c,u);it;++it){entries.emplace_back(row[u],nf+it.col(),it.value());
            entries.emplace_back(nf+it.col(),row[u],it.value());}}
    const auto& component=saddle.pressureComponents();const auto& gauge=saddle.pressureGaugeByComponent();
    const auto& volumes=layout.domain().cellVolumes();
    for(int p=0;p<np;++p){const int group=gauge.at(component.at(p));if(group<0)continue;
        entries.emplace_back(nf+p,nf+np+group,volumes[p]);entries.emplace_back(nf+np+group,nf+p,volumes[p]);}
    Matrix mixed(nf+np+ng,nf+np+ng);mixed.setFromTriplets(entries.begin(),entries.end());mixed.makeCompressed();
    const Matrix actualMixed=saddle.matrix();matrix(mixed,actualMixed,"actual system differs from independent raw A/C/volume mixed equations");
    const auto& solved=actual.solveResult().values();const auto& known=actual.boundary().normalVelocity;
    const Eigen::Map<const Eigen::VectorXd> prescribed(known.data(),n),u(solved.normalVelocity().data(),n),p(solved.pressure().data(),np);
    const Eigen::VectorXd knownAction=raw*prescribed,aU=raw*u,cP=c*p,knownContinuity=c.transpose()*prescribed;
    std::vector<double> traction(n);for(const auto& port:actual.boundary().pressureOutlets)
        for(Matrix::InnerIterator it(c,port.velocityDof);it;++it)traction[port.velocityDof]+=it.value()*port.pressure;
    Eigen::VectorXd b=Eigen::VectorXd::Zero(mixed.rows()),x=Eigen::VectorXd::Zero(mixed.rows());
    double momentumError=0,momentumScale=0;
    for(int i=0;i<nf;++i){const int v=free[i];const double momentum=rhs[v]+traction[v];
        check(known[v]==0,"raw mixed free face has a prescribed value");b[i]=momentum-knownAction[v];x[i]=u[v];
        momentumError=std::hypot(momentumError,aU[v]+cP[v]-momentum);
        momentumScale=std::hypot(momentumScale,std::abs(aU[v])+std::abs(cP[v])+std::abs(momentum));}
    b.segment(nf,np)=-knownContinuity;x.segment(nf,np)=p;
    check(b.allFinite()&&x.allFinite(),"independent eliminated raw mixed RHS is nonfinite");
    const Eigen::VectorXd residual=mixed*x-b;
    scaled(momentumError,0,momentumScale,"independent raw momentum equation",config.momentumRelativeTolerance);
    std::vector<double> continuityScale(np);for(int v=0;v<n;++v)for(Matrix::InnerIterator it(c,v);it;++it)
        continuityScale[it.col()]+=std::abs(it.value()*u[v]);
    for(int q=0;q<np;++q)scaled(residual[nf+q],0,continuityScale[q],"independent raw mixed continuity",config.continuityRelativeTolerance);
    std::vector<double> mean(ng),gross(ng);for(int q=0;q<np;++q){const int g=gauge.at(component.at(q));if(g>=0){
        const double term=volumes[q]*p[q];mean[g]+=term;gross[g]+=std::abs(term);}}
    for(int g=0;g<ng;++g)scaled(mean[g],0,gross[g],"independent raw mixed pressure gauge",config.pressureMeanRelativeTolerance);
    const auto& gates=actual.solveResult().report();
    check(gates.algebraicRelativeResidual<=config.algebraicRelativeTolerance&&gates.freeMomentumDefect<=config.momentumRelativeTolerance
        &&gates.maxCellContinuityDefect<=config.continuityRelativeTolerance&&gates.maxPressureMeanDefect<=config.pressureMeanRelativeTolerance
        &&gates.maxBoundaryFluxDefect<=config.boundaryFluxRelativeTolerance,"raw mixed solve weakened an original physical/algebraic gate");
}
void verify(const QuadtreeMesh& mesh,const HydroFrozenMidpoint& actual,const HydroSaddleConfig& config,
    const char* context,bool ownThreads=false){
    Restore restore;check(std::fesetround(FE_TONEAREST)==0,"Frozen reference cannot set nearest controls");
#if defined(_M_X64) || defined(__SSE2__)
    _mm_setcsr(_mm_getcsr()&~0x8040u);
#endif
    const auto& dual=actual.material();const auto& layout=dual.layout();const int n=layout.velocityCount();const double dt=dual.deltaTime();
    const bool linear=actual.scheme()==HydroFrozenMidpointScheme::FrozenConvection;
    check(std::isnan(actual.saddle().velocityBlock().symmetryDefect()),"actual Frozen combined symmetry diagnostic was not marked unrequested");
    if(ownThreads){
#ifdef _OPENMP
        omp_set_dynamic(0);
#endif
        check(std::feclearexcept(FE_ALL_EXCEPT)==0&&std::feraiseexcept(FE_INEXACT)==0,"Frozen reference cannot prime actual sticky inexact");
        errno=173;const Floating entry;std::optional<HydroFrozenMidpoint> first;std::optional<Floating> firstFp;
        for(int threads:{1,2,8}){
#ifdef _OPENMP
            omp_set_num_threads(threads);
#endif
            entry.restore();
            const auto prepared=HydroFrozenMidpoint::prepare(mesh,dual,actual.convection(),actual.appliedConvection(),actual.viscosity(),
                actual.pressureLoad(),actual.beginValues(),actual.forces(),actual.boundary(),actual.scheme(),dual.sources(),dual.timeLayer(),config);
            const Floating completed; // Before any floating comparisons or reference arithmetic.
            payload(actual,prepared);
            if(first){check(firstFp->same(completed),"Frozen own-request 1/2/8 complete FP/errno differ");payload(*first,prepared);++threadPairs;}
            else{first=prepared;firstFp=completed;}
        }
        ++threadGroups;
    }
    // Independent old public H and B calls, not the new private joint producer.
    const Matrix& k=actual.viscosity().viscosity().matrix();
    Matrix physical=.5*k;for(int u=0;u<n;++u)physical.coeffRef(u,u)+=dual.endMass().masses()[u]/dt;physical.makeCompressed();
    Matrix h(n,n);std::vector<long double> boundary(n);Matrix raw=physical;
    if(linear){h=actual.convection().frozenMomentumOperator(mesh,dual.sources(),dual.timeLayer());
        boundary=actual.convection().frozenBoundaryMomentumLoad(mesh,dual.sources(),dual.timeLayer());raw+=.5*h;}
    h.makeCompressed();raw.makeCompressed();
    matrix(k,actual.viscosityMatrix(),"Frozen actual K differs from its genuine phase producer");
    matrix(physical,actual.physicalVelocityMatrix(),"Frozen actual physical A differs from independent old literal arithmetic");
    matrix(h,actual.frozenConvectionMatrix(),"Frozen actual H differs from independent old public entry");
    matrix(raw,actual.velocityMatrix(),"Frozen actual raw A differs from independent old literal arithmetic");
    const auto reference=[&]{
        if(!linear)return HydroSaddleVelocityBlock::build(layout,raw);
        const Matrix transposed=raw.transpose();Matrix symmetric=.5*(raw+transposed),skew=.5*(raw-transposed);
        symmetric.makeCompressed();skew.makeCompressed();
        return HydroSaddleVelocityBlock::build(layout,symmetric).withSkew(skew);
    }();
    check(std::isfinite(reference.symmetryDefect())&&reference.spdReport().tier!=HydroVelocitySpdTier::None,
        "retained public literal velocity path lost its complete requested report/SPD capability");
    const Matrix transpose=raw.transpose(),symmetric=.5*(raw+transpose),skew=.5*(raw-transpose);
    const Matrix legacy=linear?Matrix(symmetric+skew):raw;
    matrix(reference.matrix(),legacy,"retained public S/W path changed its independent original roundtrip");
    check(&actual.velocityMatrix()==&actual.saddle().velocityBlock().matrix(),
        "actual Frozen raw A is duplicated instead of sharing the solver's unique velocity owner");
    matrix(raw,actual.saddle().velocityBlock().matrix(),"actual solver did not consume the independent raw A directly");
    const Eigen::Map<const Eigen::VectorXd> u0(actual.beginValues().normalVelocity().data(),n);
    const Eigen::VectorXd ku0=k*u0,hu0=h*u0;
    vector(ku0,actual.beginViscousAction(),"Frozen actual K*u0 differs from independent Eigen action");
    vector(hu0,actual.beginFrozenAction(),"Frozen actual H*u0 differs from independent public H action");
    const auto acceleration=actual.forces().bodyAcceleration;const bool bodyActive=acceleration.x!=0||acceleration.y!=0;
    std::vector<long double> body(n);std::vector<double> rhs(n);
    for(int u=0;u<n;++u){const double g=layout.normalAxis(u)==HydroDofLayout::NormalAxis::First?acceleration.x:acceleration.y;
        if(g!=0)body[u]=static_cast<long double>(dual.midpointMass().masses()[u])*g;
        check(body[u]==actual.bodyLoads()[u]&&(body[u]!=0||std::signbit(body[u])==std::signbit(actual.bodyLoads()[u])),
            "Frozen actual body load differs from independent midpoint mass/acceleration");
        long double value=static_cast<long double>(dual.beginMass().masses()[u])*u0[u]/dt-.5L*ku0[u];
        if(linear){value-=.5L*hu0[u];value+=boundary[u];}else value+=actual.appliedConvection().momentumLoads()[u];
        if(actual.hasInterfaceLoad())value+=actual.interfaceLoads()[u];if(actual.hasMaxwellLoad())value+=actual.maxwellLoads()[u];
        if(bodyActive)value+=body[u];rhs[u]=static_cast<double>(value);
    }
    bits(rhs,actual.integralMomentumRhs(),"Frozen actual RHS differs from independent old public H/B and literal arithmetic");
    const auto& plan=actual.matrixAssemblyPlanReport();
    check(!plan.physicalColdUsed||(plan.physicalColdCompiled&&plan.physicalValuesUsed&&!plan.physicalPlanHit),
        "Frozen physical cold current-use receipt contradicts its actual route");
    check(!plan.unionColdUsed||(plan.unionColdCompiled&&plan.unionValuesUsed&&!plan.unionPlanHit),
        "Frozen union cold current-use receipt contradicts its actual route");
    check(plan.valuesMaximumWorkers>=0&&plan.valuesParallelUsed==(plan.valuesMaximumWorkers>1),
        "Frozen values receipt fabricates a successful numeric team");
    check(plan.valuesScratchBytes<=128u*1024u*1024u,"Frozen optional worker scratch exceeds the combined128MiB scope");
    const auto& block=actual.saddle().velocityBlock();const auto& split=block.splitPlanReport();
    check(!split.requested&&!split.eligible&&!split.planHit&&!split.symmetricValuesUsed&&!split.skewValuesUsed
        &&split.inputNonZeros==0&&split.unionNonZeros==0&&split.retainedBytes==0&&!split.parallelUsed&&!split.parallelFallback&&split.maximumWorkers==0,
        "direct raw A production unexpectedly requested an S/W roundtrip");
    const auto& spd=block.spdReport();const auto& validation=block.sparseValidationReport();
    check(spd.tier==HydroVelocitySpdTier::None&&!spd.gramAttempted&&!spd.gramBoundsAvailable&&spd.certificateSeconds==0&&spd.lltSeconds==0
        &&!validation.requested&&!validation.eligible&&!validation.exactSymmetryUsed&&!validation.exactSkewUsed
        &&!validation.appendPlanUsed&&validation.retainedBytes==0&&validation.appendScratchBytes==0,
        "direct raw A production requested obsolete SPD/symmetry/skew certification");
    mixedOriginalEquation(mesh,actual,raw,rhs,config);
    const auto& execution=actual.solveResult().report();
    check(execution.orthogonalizationParallelCalls>=0&&execution.orthogonalizationMaximumWorkers>=0
        &&(execution.orthogonalizationParallelCalls>0)==(execution.orthogonalizationMaximumWorkers>1),"Frozen CGS2 actual-team receipt is inconsistent");
    if(execution.solveBreakdownTimingsAvailable)for(double t:{execution.krylovPreconditionApplySeconds,execution.krylovSpmvSeconds,
        execution.krylovOrthogonalizationSeconds,execution.krylovTrueResidualSeconds,execution.physicalAssessmentSeconds})
        check(std::isfinite(t)&&t>=0,"Frozen CGS2 available solve breakdown is invalid");
    std::cout<<"FROZEN_PREPARE_REFERENCE_PASS,context,"<<context<<",scheme,"<<(linear?"linear":"lagged")
        <<",RZ,"<<layout.domain().axisymmetric()<<",graded,"<<!layout.domain().uniformLeaves()<<",matrix_words,"<<raw.nonZeros()
        <<",observed_physical_hit,"<<plan.physicalPlanHit<<",observed_union_hit,"<<plan.unionPlanHit
        <<",physical_cold_compiled,"<<plan.physicalColdCompiled<<",physical_cold_used,"<<plan.physicalColdUsed
        <<",union_cold_compiled,"<<plan.unionColdCompiled<<",union_cold_used,"<<plan.unionColdUsed
        <<",numeric_workers,"<<plan.valuesMaximumWorkers<<",numeric_fallback,"<<plan.valuesParallelFallback
        <<",own_1_2_8,"<<ownThreads<<",own_pairs,"<<(ownThreads?2:0)<<",unique_rawA_solver_owner,1,independent_raw_mixed_equations,1,legacy_public_SW_reference,1,cross_request_FP_comparison,0\n";
}
} // namespace FrozenPrepareContract

void movingCase(const QuadtreeMesh& mesh,const HydroPhaseGeometry& geometry,double rho,double mu,int count)
{
    auto accepted=HydroTwoPhaseState::create(mesh,geometry,stepVelocity(mesh,geometry.layout()));
    nativeConstraints(mesh,accepted.values());
    HydroTwoPhaseStepConfig config;config.timeScheme=HydroTwoPhaseTimeScheme::FullyCoupledMidpoint;
    config.liquidDensity=rho;config.gasDensity=1;config.liquidViscosity=mu;config.gasViscosity=.01*mu;
    config.spaceTime.characteristic.maxStep=.00025;
    config.spaceTime.initialPanels=1;
    const auto initialHandle=accepted;const auto initialCopy=capture(accepted.publication(),accepted.values());
    for(int i=0;i<count;++i){const auto before=accepted;const auto saved=capture(accepted.publication(),accepted.values());const auto oldClock=accepted.timeLayer();
        const HermiteGalsTimeLayer clock{oldClock.step+1,oldClock.time(),oldClock.time()+.001,Position::End};
        std::vector<HydroTwoPhaseIterationReport> observed;
        config.onIteration=[&](const HydroTwoPhaseIterationReport& r){unchanged(saved,accepted.publication(),accepted.values());
            inspectIteration(r,static_cast<int>(observed.size())+1,clock,
                observed.empty()?static_cast<std::size_t>(config.spaceTime.initialPanels):observed.back().acceptedPanels);observed.push_back(r);};
        const auto step=HydroTwoPhaseStep::prepare(mesh,accepted,clock,config);unchanged(saved,accepted.publication(),accepted.values());
        const auto& unrequestedSpd=step.saddle().velocityBlock().spdReport();
        check(unrequestedSpd.tier==HydroVelocitySpdTier::None&&!unrequestedSpd.gramAttempted&&!unrequestedSpd.gramBoundsAvailable
            &&unrequestedSpd.certificateSeconds==0&&unrequestedSpd.lltSeconds==0,
            "lagged Frozen production still requested a global SPD audit");
        inspectIterationObserver(step,observed);
        check(step.beginState().sameIdentity(accepted),"step replaced begin owner");inspectStep(mesh,step,true);
        if(rho==1e6&&mu==0&&i==0){
            // This owner does not expose its stored Frozen object. Reprepare
            // from the actual final stage owners, then compare to that Step.
            HydroFrozenMidpointForces forces;forces.bodyAcceleration=config.bodyAcceleration;
            const auto& dual=step.dualMaterial();
            const auto frozen=HydroFrozenMidpoint::prepare(mesh,dual,step.convection(),step.convectionResult(),step.viscosity(),
                HydroPressureLoad::build(dual.geometry().metrics().divergence()),step.beginState().values(),forces,
                HydroFrozenMidpointBoundary::closed(step.beginState().boundary()),HydroFrozenMidpointScheme::LaggedConvection,
                dual.sources(),dual.timeLayer());
            FrozenPrepareContract::matrix(frozen.velocityMatrix(),step.velocityMatrix(),"lagged reprepare changed actual Step A words/CSR");
            FrozenPrepareContract::matrix(frozen.saddle().velocityBlock().matrix(),step.saddle().velocityBlock().matrix(),"lagged reprepare changed actual Step velocity block");
            FrozenPrepareContract::bits(frozen.integralMomentumRhs(),step.integralMomentumRhs(),"lagged reprepare changed actual Step RHS words");
            FrozenPrepareContract::verify(mesh,frozen,{},"actual_FirstJet_highRho",true);++FrozenPrepareContract::laggedReferences;
        }
        if(rho==1e3&&mu==.01&&i==0){
            auto relaxedConfig=config;relaxedConfig.velocityUpdateTolerance=1e-8;relaxedConfig.geometryUpdateTolerance=1e-7;
            relaxedConfig.momentumPowerTolerance=1e-10;
            std::vector<HydroTwoPhaseIterationReport> relaxedObserved;
            relaxedConfig.onIteration=[&](const HydroTwoPhaseIterationReport& r){unchanged(saved,accepted.publication(),accepted.values());
                inspectIteration(r,static_cast<int>(relaxedObserved.size())+1,clock,
                    relaxedObserved.empty()?static_cast<std::size_t>(relaxedConfig.spaceTime.initialPanels):relaxedObserved.back().acceptedPanels);
                relaxedObserved.push_back(r);};
            const auto relaxed=HydroTwoPhaseStep::prepare(mesh,accepted,clock,relaxedConfig);
            unchanged(saved,accepted.publication(),accepted.values());inspectIterationObserver(relaxed,relaxedObserved);
            inspectStep(mesh,relaxed,true,relaxedConfig.velocityUpdateTolerance,relaxedConfig.geometryUpdateTolerance,relaxedConfig.momentumPowerTolerance);
            compareRelaxedEndpoint(step,relaxed);
        }
        step.commit(mesh,accepted);check(accepted.sameIdentity(step.candidateState())&&!accepted.sameIdentity(before),"successful step did not atomically replace entire handle");
        check(sameClock(clock,accepted.timeLayer()),"committed complete clock changed");unchanged(saved,before.publication(),before.values());
        const auto committedCopy=capture(accepted.publication(),accepted.values());stepRejects([&]{step.commit(mesh,accepted);});
        unchanged(committedCopy,accepted.publication(),accepted.values());++acceptedSteps;}
    unchanged(initialCopy,initialHandle.publication(),initialHandle.values());++movingCases;
}

// Reuse the original moving circle/curl with the adopted default AC scheme.
void linearTimeCase(const QuadtreeMesh& mesh,const HydroPhaseGeometry& geometry,bool inviscid=false)
{
    auto accepted=HydroTwoPhaseState::create(mesh,geometry,stepVelocity(mesh,geometry.layout()));
    const auto original=accepted;const auto originalCopy=capture(accepted.publication(),accepted.values());
    HydroTwoPhaseStepConfig config;config.liquidDensity=1000;config.gasDensity=1;
    config.liquidViscosity=inviscid?0:.01;config.gasViscosity=inviscid?0:.0001;config.surfaceTension=1;
    check(config.timeScheme==HydroTwoPhaseTimeScheme::AdaptiveLinearMidpoint,"linear family did not inherit default AC");
    check(config.spaceTime.integrationRule==HydroSpaceTimeIntegrationRule::EmbeddedSelectiveBoole,"linear family did not inherit endpoint-reusing quadrature");
    config.interfaceMismatchTolerance=1e6;
    config.spaceTime.characteristic.maxStep=.00025;config.spaceTime.initialPanels=1;
    const auto prepareObserved=[&](HermiteGalsTimeLayer clock,const HydroTwoPhaseStepConfig& settings){
        auto actualConfig=settings;const auto saved=capture(accepted.publication(),accepted.values());
        std::vector<HydroTwoPhaseIterationReport> observed;
        actualConfig.onIteration=[&](const HydroTwoPhaseIterationReport& r){
            unchanged(saved,accepted.publication(),accepted.values());
            inspectIteration(r,static_cast<int>(observed.size())+1,clock,
                observed.empty()?static_cast<std::size_t>(actualConfig.spaceTime.initialPanels):observed.back().acceptedPanels);
            observed.push_back(r);};
        const auto result=HydroTwoPhaseStep::prepare(mesh,accepted,clock,actualConfig);
        unchanged(saved,accepted.publication(),accepted.values());
        check(result.beginState().sameIdentity(accepted)&&result.report().timeScheme==settings.timeScheme,
            "linear stage changed its begin owner or selected scheme");
        check(sameClock(result.candidateState().timeLayer(),clock)
            &&sameClock(result.materialTransport().endGeometry().timeLayer(),clock),
            "linear stage lost its actual endpoint full clock");
        check(result.materialTransport().spaceTimeReport().integrationRule==settings.spaceTime.integrationRule,"linear stage did not use the actual requested integration rule");
        inspectIterationObserver(result,observed);inspectStep(mesh,result,true);
        check(result.report().interfaceMismatch==result.report().history.front().interfaceMismatch
            &&result.report().finalInterfaceMismatch==result.report().history.back().interfaceMismatch,
            "linear mismatch diagnostics replaced the first/final actual stages");
        ++linearCases;return result;
    };
    const HermiteGalsTimeLayer firstClock{1,0,.001,Position::End};
    const auto startup=prepareObserved(firstClock,config);
    check(!startup.report().usedVelocityHistory&&startup.report().startupCorrection&&startup.report().corrected
        &&startup.report().iterations==2,"history-free startup did not force exactly two stages");
    startup.commit(mesh,accepted);check(accepted.sameIdentity(startup.candidateState()),"linear startup did not commit");
    unchanged(originalCopy,original.publication(),original.values());

    // Variable dt distinguishes the stored accepted history from a copied guess.
    const auto beforeSecond=accepted;const auto beforeSecondCopy=capture(accepted.publication(),accepted.values());
    const HermiteGalsTimeLayer secondClock{2,.001,.0025,Position::End};
    const auto one=prepareObserved(secondClock,config);
    check(one.report().usedVelocityHistory&&!one.report().startupCorrection&&!one.report().corrected
        &&one.report().iterations==1&&one.report().interfaceMismatch<=config.interfaceMismatchTolerance,
        "accepted history did not select the one-stage extrapolated solve");
    long double guessGap=0,guessWork=0,midGap=0,midWork=0;
    const auto& previous=original.values().normalVelocity();const auto& current=beforeSecond.values().normalVelocity();
    const auto& endpoint=one.endpointGuess().normalVelocity();
    const auto& midpoint=one.materialTransport().midpointRate().spatialAdvector().values().normalVelocity();
    const long double ratio=(secondClock.end-secondClock.begin)/(firstClock.end-firstClock.begin);
    for(std::size_t i=0;i<current.size();++i){
        const long double expected=current[i]+ratio*(static_cast<long double>(current[i])-previous[i]);
        guessGap=std::hypot(guessGap,static_cast<long double>(endpoint[i])-expected);
        guessWork=std::hypot(guessWork,std::abs(current[i])+std::abs(ratio*(static_cast<long double>(current[i])-previous[i])));
        const long double expectedMid=.5L*(current[i]+expected);
        midGap=std::hypot(midGap,static_cast<long double>(midpoint[i])-expectedMid);
        midWork=std::hypot(midWork,.5L*(std::abs(current[i])+std::abs(expected)));
    }
    scaled(guessGap,0,guessWork,"actual accepted velocity-history endpoint extrapolation",1e-12);
    scaled(midGap,0,midWork,"actual history-extrapolated material midpoint",1e-12);
    check(one.report().capillaryCourant>0&&!one.report().capillaryForcedCorrection,
        "weak capillary history solve incorrectly forced a correction");
    one.commit(mesh,accepted);check(accepted.sameIdentity(one.candidateState()),"history one-stage solve did not commit");
    unchanged(beforeSecondCopy,beforeSecond.publication(),beforeSecond.values());

    const HermiteGalsTimeLayer thirdClock{3,.0025,.0035,Position::End};
    auto forcedConfig=config;forcedConfig.interfaceMismatchTolerance=1e-14;
    const auto corrected=prepareObserved(thirdClock,forcedConfig);
    check(corrected.report().usedVelocityHistory&&!corrected.report().startupCorrection&&corrected.report().corrected
        &&corrected.report().iterations==2&&corrected.report().interfaceMismatch>forcedConfig.interfaceMismatchTolerance
        &&corrected.report().correctionRelative>0,"small mismatch trigger did not perform exactly one real correction");

    const auto spacing=geometry.layout().domain().spacing();const double h=std::min(spacing.x,spacing.y);
    const double expectedCourant=(thirdClock.end-thirdClock.begin)
        /std::sqrt(.5*(config.liquidDensity+config.gasDensity)*h*h*h/(std::acos(-1.)*config.surfaceTension));
    auto capillaryConfig=config;capillaryConfig.capillaryCorrectionCourant=.5*expectedCourant;
    const auto capillaryCorrected=prepareObserved(thirdClock,capillaryConfig);
    check(capillaryCorrected.report().capillaryForcedCorrection&&capillaryCorrected.report().iterations==2
        &&capillaryCorrected.report().corrected&&!capillaryCorrected.report().startupCorrection
        &&capillaryCorrected.report().interfaceMismatch<=config.interfaceMismatchTolerance,
        "capillary Courant did not independently trigger exactly one AC correction");
    scaled(capillaryCorrected.report().capillaryCourant,expectedCourant,expectedCourant,
        "actual capillary Courant",1e-12);

    auto resetConfig=config;resetConfig.timeScheme=HydroTwoPhaseTimeScheme::ExtrapolatedMidpoint;resetConfig.resetVelocityHistory=true;
    const auto reset=prepareObserved(thirdClock,resetConfig);
    check(!reset.report().usedVelocityHistory&&reset.report().startupCorrection&&reset.report().corrected
        &&reset.report().iterations==2,"explicit velocity-history reset did not restore two-stage startup");

    const auto retained=accepted;const auto saved=capture(accepted.publication(),accepted.values());
    const auto reject=[&](const HydroTwoPhaseStepConfig& settings,const char* expected,bool correction){
        try{HydroTwoPhaseStep::prepare(mesh,accepted,thirdClock,settings);}
        catch(const std::exception& e){
            check(std::string(e.what()).find(expected)!=std::string::npos,
                std::string("linear rejection reached a different stage: ")+e.what());
            if(correction){const auto* retry=dynamic_cast<const HydroTwoPhaseStepRetry*>(&e);
                check(retry&&retry->reason()==HydroTwoPhaseRetryReason::Correction,"correction retry lost its typed reason");}
            check(accepted.sameIdentity(retained),"linear failure changed the accepted owner/history handle");
            unchanged(saved,accepted.publication(),accepted.values());++linearRejects;return;
        }
        throw std::runtime_error("invalid linear midpoint preparation accepted");
    };
    auto sensorConfig=config;sensorConfig.dissipation=HydroConvectionDissipationConfig{};
    reject(sensorConfig,"linear midpoint does not support nonlinear sensor dissipation",false);
    auto boundedConfig=forcedConfig;
    boundedConfig.maximumCorrectionRelative=std::min(1e-14,.5*corrected.report().correctionRelative);
    reject(boundedConfig,"linear midpoint correction exceeds configured bound",true);
    auto capillaryBoundConfig=capillaryConfig;capillaryBoundConfig.maximumCapillaryCourant=.5*expectedCourant;
    capillaryBoundConfig.capillaryCorrectionCourant=.25*expectedCourant;
    bool capillaryRejected=false;
    try{HydroTwoPhaseStep::prepare(mesh,accepted,thirdClock,capillaryBoundConfig);}
    catch(const HydroTwoPhaseStepRetry& e){
        check(e.reason()==HydroTwoPhaseRetryReason::Capillary,"capillary retry lost its typed reason");
        check(accepted.sameIdentity(retained),"capillary retry changed accepted owner/history");
        unchanged(saved,accepted.publication(),accepted.values());capillaryRejected=true;++linearRejects;
    }
    check(capillaryRejected,"over-budget capillary step was accepted");
    std::cout<<"LINEAR_TIME_SCHEME_PASS,"<<startup.report().iterations<<','<<one.report().iterations<<','
        <<corrected.report().iterations<<','<<reset.report().iterations<<','<<corrected.report().interfaceMismatch
        <<','<<corrected.report().correctionRelative<<','<<guessGap<<','<<linearRejects<<",PASS\n";
}

// One additional step on the original n4 moving field/geometry, with sigma=1.
void capillaryCase(const QuadtreeMesh& mesh,const HydroPhaseGeometry& geometry)
{
    auto accepted=HydroTwoPhaseState::create(mesh,geometry,stepVelocity(mesh,geometry.layout()));const auto retained=accepted;
    const auto saved=capture(accepted.publication(),accepted.values());const HermiteGalsTimeLayer clock{1,0,.001,Position::End};
    HydroTwoPhaseStepConfig config;config.timeScheme=HydroTwoPhaseTimeScheme::FullyCoupledMidpoint;
    config.liquidDensity=1000;config.gasDensity=1;config.surfaceTension=1;
    config.spaceTime.characteristic.maxStep=.00025;config.spaceTime.initialPanels=1;
    auto limited=config;limited.maximumIterations=1;
    capillaryRejectsAction([&]{HydroTwoPhaseStep::prepare(mesh,accepted,clock,limited);},"corrector did not converge");
    check(accepted.sameIdentity(retained),"failed capillary corrector changed accepted handle");unchanged(saved,accepted.publication(),accepted.values());
    std::vector<HydroTwoPhaseIterationReport> observed;
    config.onIteration=[&](const HydroTwoPhaseIterationReport& r){unchanged(saved,accepted.publication(),accepted.values());
        inspectIteration(r,static_cast<int>(observed.size())+1,clock,
                observed.empty()?static_cast<std::size_t>(config.spaceTime.initialPanels):observed.back().acceptedPanels);observed.push_back(r);};
    const auto step=HydroTwoPhaseStep::prepare(mesh,accepted,clock,config);unchanged(saved,accepted.publication(),accepted.values());
    inspectIterationObserver(step,observed);
    check(step.beginState().sameIdentity(accepted)&&step.hasCapillaryLoad()&&step.report().surfaceTension==config.surfaceTension,"moving capillary step lost begin/load owner or physical sigma");
    inspectStep(mesh,step,true);
    long double forceNorm=0,work=0;const auto& force=step.capillaryLoad().momentumRhsLoad();const auto& bar=step.solvedMidpointValues().normalVelocity();
    for(std::size_t i=0;i<force.size();++i){forceNorm=std::hypot(forceNorm,static_cast<long double>(force[i]));work+=static_cast<long double>(bar[i])*force[i];}
    check(forceNorm>0&&std::isfinite(work)&&work!=0,"moving capillary step erased actual force/work");
    step.commit(mesh,accepted);check(accepted.sameIdentity(step.candidateState())&&!accepted.sameIdentity(retained)
        &&sameClock(accepted.timeLayer(),clock),"successful capillary step did not atomically commit complete state/full clock");
    unchanged(saved,retained.publication(),retained.values());++capillarySteps;
    std::cout<<"CAPILLARY_STEP_PASS,"<<forceNorm<<','<<work<<','<<step.report().iterations<<'\n';
}
void inspectRzPublication(const QuadtreeMesh& mesh,const HydroPhaseGeometry& phase,bool analytic=false)
{
    const auto& pub=phase.publication();const auto& state=pub.state();const auto& primary=state.primary();
    check(primary.measure==HermiteGalsMeasure::AxisymmetricVolume&&primary.transport==HermiteGalsTransport::FirstJet
        &&primary.phaseVolumes.empty(),"RZ publication lost true measure or acquired inventory");
    for(const auto& v:primary.vertices)if(point(mesh,v.vertex).x==0)check(v.dx==0,"RZ scalar FirstJet lost exact even axis derivative");
    for(const auto& e:primary.edges)if(point(mesh,e.vertices[0]).x==0&&point(mesh,e.vertices[1]).x==0)
        check(e.dn==0,"RZ scalar edge FirstJet lost exact even axis derivative");
    std::size_t axes=0;
    const auto sample=[&](const HermiteGalsSurfaceSample& s){
        check(s.geometryIdentity==phase.geometryIdentity()&&sameClock(s.timeLayer,phase.timeLayer()),"RZ interface sample lost actual source/full clock");
        const auto& f=s.frame;const auto jet=state.geometry().snapshot().scalarJetAtNativeCellChart(f.nativeCellIndex,f.cellChartPoint);
        const auto norm=std::hypot(jet.dx,jet.dy);check(norm>0&&std::isfinite(s.curvature)&&std::isfinite(s.azimuthalCurvature),"RZ physical curvature is not finite");
        close(f.normal.x,jet.dx/norm,2e-13,"RZ normal is detached from actual scalar dx");
        close(f.normal.y,jet.dy/norm,2e-13,"RZ normal is detached from actual scalar dy");
        if(f.point.x==0){check(f.normal.x==0&&jet.dx==0&&s.azimuthalCurvature==f.curvature
                &&s.curvature==2*f.curvature,"RZ physical curvature lost the actual axis even limit");++axes;}
        if(analytic){close(f.curvature,1/.27,2e-8,"RZ initial sphere meridional curvature");
            close(s.curvature,2/.27,2e-8,"RZ initial sphere true physical curvature");}
    };
    for(std::size_t c=0;c<pub.topology().curves.size();++c)for(std::size_t g=0;g<pub.topology().curves[c].graphs().size();++g){
        sample(pub.sample(c,g,0));sample(pub.sample(c,g,1));
        for(const auto& q:pub.quadrature(c,g)){sample(q.sample);
            scaled(q.physicalWeight,2*pi*q.sample.frame.point.x*q.weight,q.physicalWeight,"RZ actual revolved interface area weight",1e-12);}}
    check(axes==2&&axes==pub.report().axisEndpoints,"RZ moving interface lost its two actual axis poles");
    long double liquid=0,capacity=0;
    for(int p=0;p<phase.layout().pressureCount();++p){const auto& cell=mesh.cells(mesh.cellIndexToLabel[phase.layout().nativeCell(p)]);
        const auto lo=mesh.verts(cell.cellToVert[0]).coord(),hi=mesh.verts(cell.cellToVert[2]).coord();
        const auto volume=stepBoxMeasure(lo,hi,true);const auto split=phase.cellVolumes()[p];
        check(split.liquid>=0&&split.gas>=0,"RZ actual cell has negative sharp phase measure");
        scaled(phase.metrics().cellVolumes()[p],volume,volume,"RZ actual native revolved cell volume",1e-12);
        scaled(split.liquid+split.gas,volume,volume,"RZ true two-phase cell capacity",1e-12);
        liquid+=split.liquid;capacity+=volume;}
    scaled(capacity,pi,pi,"RZ actual complete cylinder volume",1e-12);
    check(liquid>0&&liquid<capacity,"RZ fixture is not an actual two-phase fluid");
    for(int i=0;i<phase.layout().velocityCount();++i){const auto area=stepFaceMeasure(mesh,phase.layout(),i);
        scaled(phase.metrics().faceMeasures()[i],area,area,"RZ independent actual native face area",1e-12);}
    if(analytic){const long double volume=4*pi*.27L*.27L*.27L/3,area=4*pi*.27L*.27L;
        scaled(liquid,volume,capacity,"RZ initial exact sphere true volume",1e-10);
        // Fixed publication GL8 has no adaptive arclength error contract.
        // Verify analytic area through the existing controlled consumer of
        // these SAME retained chains; keep production GL8 as an observation.
        long double controlledArea=0;
        for(const auto& curve:pub.topology().curves){const LevelSetInterfacePanels material(curve);
            for(const auto& panel:material.panels())controlledArea+=panel.physicalMeasure;}
        scaled(controlledArea,area,area,"RZ initial exact sphere true surface area",1e-10);
        const auto fixedArea=actualInterfaceMeasure(phase);
        std::cout<<"RZ_INITIAL_SPHERE_AREA,"<<area<<','<<controlledArea<<','<<fixedArea<<','
            <<std::abs(fixedArea-area)/area<<",FIXED_GL8_OBSERVATION\n";}
}
// One actual open RZ step: independent unique-F incidence, physical inlet
// and outlet energy, and prescribed-normal reaction work. Closed inspectStep
// deliberately remains unchanged and is not the oracle for open boundaries.
void inspectPrescribedFlowStep(const QuadtreeMesh& mesh,const HydroTwoPhaseStep& step)
{
    const auto& material=step.materialTransport();const auto& dual=step.dualMaterial();const auto& l=dual.layout();
    const auto& op=step.convection();const auto& result=step.reassessedConvectionResult();const auto& report=step.report();
    const auto* flow=step.beginState().prescribedFlow();const auto sources=dual.sources();const double dt=dual.deltaTime();
    const bool pressureOutlet=flow&&!flow->pressureOutlets.empty();
    std::vector<bool> outlet(l.velocityCount());std::vector<long double> traction(l.velocityCount());
    if(pressureOutlet)for(const auto& port:flow->pressureOutlets){const auto face=l.domain().face(mesh,l.nativeFace(port.velocityDof));
        outlet[port.velocityDof]=true;traction[port.velocityDof]=-face.outwardSign*stepFaceMeasure(mesh,l,port.velocityDof)*port.pressure;}
    const auto& u0=step.beginState().values().normalVelocity();const auto& u1=step.candidateState().values().normalVelocity();
    const auto& bar=step.solvedMidpointValues().normalVelocity();const auto& pressure=step.candidateState().values().pressure();
    check(flow&&step.candidateState().prescribedFlow()&&op.isPrescribedFlow()&&!op.isImpermeable()
        &&op.periodicPairs().empty()&&l.domain().axisymmetric(),"open step lost the retained prescribed RZ flow");
    check(material.liquidDensity()==1000&&material.gasDensity()==1&&report.surfaceTension==1
        &&step.hasCapillaryLoad(),"open fixture lost its actual density contrast or capillary load");
    check(material.isSpaceTime()&&material.isSpatial()
        &&material.spaceTimeReport().firstJetBoundaryTreatment==HermiteGalsBoundaryTreatment::InflowHistory
        &&material.spaceTimeReport().firstJetTimeProjection==HermiteGalsTimeProjection::CharacteristicPullback,
        "open material did not consume explicit inlet-only FirstJet pullback");
    check(sources.begin==step.beginState().geometryIdentity()&&sources.midpoint==dual.midpointMass().geometryIdentity()
        &&sources.end==step.candidateState().geometryIdentity()&&sources.begin!=sources.midpoint&&sources.midpoint!=sources.end,
        "open step collapsed its three actual geometry owners");
    check(&dual.material().faceIntegrals()==&material.faceIntegrals()
        && &result.material().faceRates()==&dual.faceRates()
        && &step.convectionResult().material().faceRates()==&dual.faceRates(),"open consumers replaced the unique original F");
    check(material.midpointFlux().advectorIdentity()==material.midpointRate().spatialAdvector().identity()
        &&material.characteristicSource().sameIdentity(step.characteristicSource()),"open material detached its actual advector");
    const HermiteGalsTimeLayer mid{report.history.back().timeLayer.step,report.history.back().timeLayer.begin,
        report.history.back().timeLayer.end,Position::Midpoint};
    check(sameClock(material.beginGeometry().timeLayer(),step.beginState().timeLayer())
        &&sameClock(material.midpointRate().geometry().timeLayer(),mid)&&sameClock(dual.timeLayer(),mid)
        &&sameClock(material.endGeometry().timeLayer(),step.candidateState().timeLayer()),"open material full clocks differ");
    const std::array<const HydroPhaseGeometry*,3> phases{{&material.beginGeometry(),&material.midpointRate().geometry(),&material.endGeometry()}};
    const std::array<const HydroMomentumMass*,3> masses{{&dual.beginMass(),&dual.midpointMass(),&dual.endMass()}};
    for(std::size_t k=0;k<3;++k){inspectRzPublication(mesh,*phases[k]);
        check(masses[k]->geometryIdentity()==phases[k]->geometryIdentity()&&sameClock(masses[k]->timeLayer(),phases[k]->timeLayer()),
            "open three masses lost the actual scalar source/full clock");
        for(int u=0;u<l.velocityCount();++u){const auto v=masses[k]->phaseVolumes()[u];const auto& box=dual.geometry().volumes()[u];
            const auto capacity=stepBoxMeasure(box.lower,box.upper,true);check(v.liquid>=0&&v.gas>=0,"open sharp dual phase became negative");
            scaled(v.liquid+v.gas,capacity,capacity,"open true revolved dual capacity",1e-12);
            scaled(masses[k]->masses()[u],1000*v.liquid+v.gas,1000*v.liquid+v.gas,"open actual sharp three-mass formula",1e-12);}}
    const auto& motion=material.characteristicSource();
    check(motion.previous().geometryIdentity()==step.beginState().geometryIdentity()
        &&motion.arrival().geometryIdentity()==step.beginState().geometryIdentity()
        &&sameClock(motion.previous().timeLayer(),step.beginState().timeLayer())
        &&sameClock(motion.arrival().timeLayer(),step.candidateState().timeLayer()),"open endpoint motion source/full clock mismatch");
    // Reconstruct both integral views from the same actual phase atom data.
    const auto& reconstruction=step.beginState().reconstruction();const auto& atom=material.phaseAtomIntegrals();
    const auto sumAtoms=[&](const std::vector<int>& ids){std::array<long double,2> sum{};
        for(int id:ids){const auto term=1000*atom[id].liquid+atom[id].gas;sum[0]+=term;sum[1]+=std::abs(term);}return sum;};
    std::vector<long double> cellOut(l.pressureCount()),cellGross(l.pressureCount());int in=0,out=0;
    for(int u=0;u<l.velocityCount();++u){const auto a=sumAtoms(reconstruction.primalAtoms()[u]);
        scaled(material.faceIntegrals()[u],a[0],a[1],"open original primal atom integral",1e-12);
        const auto face=l.domain().face(mesh,l.nativeFace(u));const auto area=stepFaceMeasure(mesh,l,u);
        scaled(material.faceRates()[u],material.faceIntegrals()[u]/dt,a[1]/dt,"open original primal I/dt",1e-12);
        if(face.negativeCell>=0){const int p=l.pressureDof(face.negativeCell);cellOut[p]+=material.faceIntegrals()[u];cellGross[p]+=std::abs(material.faceIntegrals()[u]);}
        if(face.positiveCell>=0){const int p=l.pressureDof(face.positiveCell);cellOut[p]-=material.faceIntegrals()[u];cellGross[p]+=std::abs(material.faceIntegrals()[u]);}
        if(face.outwardSign!=0){
            if(!outlet[u])check(u1[u]==flow->normalVelocity[u],"open solve changed a prescribed normal mean");
            if(!pressureOutlet){const auto expected=face.horizontal?1000.L*area*.2L*dt:0.L;
                scaled(material.faceIntegrals()[u],expected,std::abs(expected),"open independent liquid inlet/outlet integral",1e-12);}
            if(face.horizontal){check(material.faceIntegrals()[u]>0,"open port has no genuine liquid transport");
                if(face.outwardSign<0)++in;else ++out;}
            else check(material.faceIntegrals()[u]==0,"open impermeable radial face has material flux");
        }
    }
    check(in>0&&out>0,"open step did not exercise both physical ports");long double phaseMovement=0;
    for(int p=0;p<l.pressureCount();++p){const auto a=phases[0]->cellVolumes()[p],b=phases[2]->cellVolumes()[p];
        const auto dm=1000*(b.liquid-a.liquid)+(b.gas-a.gas),r=(dm+cellOut[p])/dt;
        scaled(material.cells()[p].defectRate,r,(std::abs(dm)+cellGross[p])/dt,"open primal actual continuity R",1e-12);
        phaseMovement+=std::abs(b.liquid-a.liquid);}
    check(phaseMovement>0,"open true two-phase interface did not move");
    std::vector<long double> expected(bar.size()),gross(bar.size()),external(bar.size());
    std::array<long double,2> boundaryMomentum{},boundaryEnergy{},physicalEnergy{},inletDissipation{},boundaryGross{};
    int incoming=0,outgoing=0;
    for(std::size_t g=0;g<dual.geometry().segments().size();++g){const auto& s=dual.geometry().segments()[g];
        const auto a=sumAtoms(reconstruction.dualAtoms()[g]);
        scaled(dual.faceIntegrals()[g],a[0],a[1],"open original dual atom integral",1e-12);
        scaled(dual.faceRates()[g],dual.faceIntegrals()[g]/dt,a[1]/dt,"open original dual I/dt",1e-12);
        const int i=s.negativeVelocity,j=s.positiveVelocity;const auto f=dual.faceRates()[g];
        if(i>=0&&j>=0){const auto left=.5L*f*bar[j],right=.5L*f*bar[i];
            expected[i]+=left;expected[j]-=right;gross[i]+=std::abs(left);gross[j]+=std::abs(right);continue;}
        const int u=i>=0?i:j,axis=l.normalAxis(u)==HydroDofLayout::NormalAxis::First?0:1;
        const auto q=(i>=0?1.L:-1.L)*f;external[u]+=q;
        const double trace=q<0?(axis==0?0:.2):q>0?bar[u]:0.;
        check(result.boundary().traces()[g]==trace,"open KEP did not retain inlet data or natural outlet native trace");
        const auto n=q*trace,energy=q*(bar[u]*trace-.5L*bar[u]*bar[u]);
        const auto physical=.5L*q*trace*trace,dissipation=q<0?.5L*(-q)*(bar[u]-trace)*(bar[u]-trace):0.L;
        expected[u]+=n;gross[u]+=std::abs(n);boundaryMomentum[axis]+=n;boundaryEnergy[axis]+=energy;
        physicalEnergy[axis]+=physical;inletDissipation[axis]+=dissipation;
        boundaryGross[axis]+=std::abs(energy)+std::abs(physical)+dissipation;
        if(q<0)++incoming;else if(q>0)++outgoing;
    }
    check(incoming>0&&outgoing>0,"open KEP had no actual incoming or outgoing unique F");
    const Matrix h=op.frozenMomentumOperator(mesh,dual.sources(),dual.timeLayer());
    const auto b=op.frozenBoundaryMomentumLoad(mesh,dual.sources(),dual.timeLayer());
    const auto hbar=matvec(h,bar),hu0=matvec(h,u0);const auto& force=step.capillaryLoad().momentumRhsLoad();
    std::vector<long double> hProductWork(bar.size());
    for(int outer=0;outer<h.outerSize();++outer)for(Matrix::InnerIterator it(h,outer);it;++it)
        hProductWork[it.row()]+=std::abs(static_cast<long double>(it.value())*bar[it.col()]);
    check(step.viscosity().liquidViscosity()==0&&step.viscosity().gasViscosity()==0
        &&step.viscosity().viscosity().matrix().nonZeros()==0&&!step.viscosity().hasQuadrature(),"open exact mu0 sampled a viscous form");
    inspectCsfOwner(mesh,step);
    Matrix expectedA=.5*h;for(int u=0;u<l.velocityCount();++u)expectedA.coeffRef(u,u)+=dual.endMass().masses()[u]/dt;
    const Matrix matrixGap=step.velocityMatrix()-expectedA;
    scaled(matrixGap.norm(),0,expectedA.norm(),"open implicit variable-mass H matrix",1e-12);
    const auto& free=step.saddle().freeVelocityDofs();std::vector<int> independentFree;
    std::vector<long double> divergence(l.pressureCount()),divergenceGross(l.pressureCount());
    std::array<long double,2> componentPower{},componentWork{},outwardKinetic{},componentMomentum{},momentumWork{};
    long double rhsGap=0,rhsWork=0,equationGap=0,equationWork=0,totalPower=0,reactionPower=0,powerWork=0,tractionPower=0;
    long double netBoundary=0,grossBoundary=0,outletChange=0,tractionNorm=0;
    long double ke0=0,ke1=0,inertial=0,cross=0,timeWork=0,pressurePower=0,capillaryPower=0,forceNorm=0,massMovement=0;
    for(int u=0;u<l.velocityCount();++u){const auto face=l.domain().face(mesh,l.nativeFace(u));const int axis=face.horizontal?1:0;
        if(face.outwardSign==0||outlet[u])independentFree.push_back(u);
        const auto begin=dual.beginMass().phaseVolumes()[u],end=dual.endMass().phaseVolumes()[u];
        const auto dm=1000*(end.liquid-begin.liquid)+(end.gas-begin.gas);long double integral=0,integralGross=0;
        for(int side=0;side<4;++side)for(int g:dual.geometry().volumes()[u].sides[side]){
            integral+=(side%2?1:-1)*dual.faceIntegrals()[g];integralGross+=std::abs(dual.faceIntegrals()[g]);}
        const auto r=(dm+integral)/dt,rs=(std::abs(dm)+integralGross)/dt;
        scaled(dual.cells()[u].defectRate,r,rs,"open dual actual sharp geometry continuity R",1e-12);
        const auto diagonal=.5L*(integral/dt-external[u])*bar[u];expected[u]+=diagonal;gross[u]+=std::abs(diagonal);
        scaled(result.momentumRates()[u],expected[u],gross[u]+std::abs(result.momentumRates()[u]),"open independent original-F centered KEP",1e-12);
        scaled(result.momentumLoads()[u],-result.momentumRates()[u]+bar[u]*r,
            std::abs(result.momentumRates()[u])+std::abs(bar[u])*rs,"open one actual uR momentum load",1e-12);
        scaled(result.momentumLoads()[u],-static_cast<long double>(hbar[u])+b[u],
            hProductWork[u]+std::abs(b[u])+std::abs(bar[u])*rs,"open frozen H plus retained inlet B action",1e-12);
        check(result.dissipationRates()[u]==0,"open sensor-free AC manufactured numerical dissipation");
        const auto m0=static_cast<long double>(dual.beginMass().masses()[u]),m1=static_cast<long double>(dual.endMass().masses()[u]);
        const auto rhs=m0*u0[u]/dt-.5L*hu0[u]+b[u]+force[u];
        rhsGap=std::hypot(rhsGap,step.integralMomentumRhs()[u]-rhs);
        rhsWork=std::hypot(rhsWork,std::abs(m0*u0[u]/dt)+.5L*std::abs(hu0[u])+std::abs(b[u])+std::abs(force[u]));
        scaled(bar[u],.5L*u0[u]+.5L*u1[u],.5L*(std::abs(u0[u])+std::abs(u1[u])),"open actual solved midpoint",1e-12);
        const auto area=stepFaceMeasure(mesh,l,u);const auto volumeFlow=area*u1[u];
        if(face.negativeCell>=0){const int p=l.pressureDof(face.negativeCell);divergence[p]+=volumeFlow;divergenceGross[p]+=std::abs(volumeFlow);}
        if(face.positiveCell>=0){const int p=l.pressureDof(face.positiveCell);divergence[p]-=volumeFlow;divergenceGross[p]+=std::abs(volumeFlow);}
        const auto pm=face.negativeCell>=0?pressure[l.pressureDof(face.negativeCell)]:0;
        const auto pp=face.positiveCell>=0?pressure[l.pressureDof(face.positiveCell)]:0;
        const auto cp=area*(pp-pm),rate=(m1*u1[u]-m0*u0[u])/dt,residual=rate+cp-result.momentumLoads()[u]-force[u];
        const auto work=std::abs(bar[u]*rate)+std::abs(bar[u]*cp)+std::abs(bar[u]*result.momentumLoads()[u])+std::abs(bar[u]*force[u])+std::abs(bar[u]*traction[u]);
        totalPower+=bar[u]*residual;powerWork+=work;tractionPower+=bar[u]*traction[u];tractionNorm+=std::abs(traction[u]);
        scaled(step.saddle().pressureOutletLoads()[u],traction[u],std::abs(traction[u]),"open original native traction T",1e-12);
        if(face.outwardSign!=0){const auto q=face.outwardSign*volumeFlow;netBoundary+=q;grossBoundary+=std::abs(q);}
        if(outlet[u]){check(flow->normalVelocity[u]==0,"free pressure outlet has a prescribed velocity payload");
            check(face.outwardSign*u1[u]>0,"actual pressure outlet is not outgoing");outletChange=std::max(outletChange,std::abs(static_cast<long double>(u1[u])-u0[u]));}
        if(face.outwardSign!=0&&!outlet[u])reactionPower+=bar[u]*residual;
        else {equationGap=std::hypot(equationGap,residual-traction[u]);equationWork=std::hypot(equationWork,
            (std::abs(m1*u1[u])+std::abs(m0*u0[u]))/dt+std::abs(cp)+std::abs(result.momentumLoads()[u])+std::abs(force[u])+std::abs(traction[u]));}
        componentMomentum[axis]+=result.momentumRates()[u];momentumWork[axis]+=gross[u]+std::abs(result.momentumRates()[u]);
        componentPower[axis]+=bar[u]*result.momentumRates()[u];componentWork[axis]+=std::abs(bar[u])*(gross[u]+std::abs(result.momentumRates()[u]));
        outwardKinetic[axis]+=.5L*bar[u]*bar[u]*integral/dt;
        const auto impulse=m1*u1[u]-m0*u0[u];ke0+=.5L*m0*u0[u]*u0[u];ke1+=.5L*m1*u1[u]*u1[u];
        inertial+=bar[u]*impulse;cross+=.5L*u0[u]*u1[u]*(m1-m0);
        timeWork+=std::abs(bar[u]*impulse)+.5L*(std::abs(m0*u0[u]*u0[u])+std::abs(m1*u1[u]*u1[u]));
        pressurePower+=bar[u]*cp;capillaryPower+=bar[u]*force[u];forceNorm=std::hypot(forceNorm,force[u]);massMovement+=std::abs(m1-m0);
    }
    check(free==independentFree,"open saddle failed to eliminate physical prescribed normals");
    check(forceNorm>0&&massMovement>0,"open moving capillary contrast fixture erased force or mass movement");
    scaled(rhsGap,0,rhsWork,"open original H/B momentum RHS",1e-12);
    scaled(equationGap,0,equationWork,"open final actual free momentum",1e-11);
    long double pressureMean=0,pressureMeanWork=0;
    for(int p=0;p<l.pressureCount();++p){scaled(divergence[p],0,divergenceGross[p],"open true volume continuity",1e-11);
        const auto& cell=mesh.cells(mesh.cellIndexToLabel[l.nativeCell(p)]);const auto lo=mesh.verts(cell.cellToVert[0]).coord(),hi=mesh.verts(cell.cellToVert[2]).coord();
        const auto term=stepBoxMeasure(lo,hi,true)*pressure[p];pressureMean+=term;pressureMeanWork+=std::abs(term);}
    if(pressureOutlet){check(step.saddle().pressureGaugeCount()==0&&step.solveResult().report().pressureGaugeCount==0
            &&std::all_of(step.saddle().pressureGaugeByComponent().begin(),step.saddle().pressureGaugeByComponent().end(),[](int g){return g==-1;}),
            "pressure outlet component retained a pressure gauge");
        check(pressureMean!=0&&tractionNorm>0&&outletChange>0,"pressure outlet lost absolute pressure, actual T, or its freely solved velocity");
        check(step.saddle().pressureOutlets().size()==flow->pressureOutlets.size(),"Step dropped a pressure traction port");
        scaled(netBoundary,0,grossBoundary,"actual pressure outlet solved volume flux compatibility",1e-12);
    }else scaled(pressureMean,0,pressureMeanWork,"open true-volume pressure gauge",1e-11);
    for(int axis:step.beginState().phaseGeometry().metrics().axisVelocityDofs())
        check(u1[axis]==0&&material.faceIntegrals()[axis]==0&&stepFaceMeasure(mesh,l,axis)==0
            &&!std::binary_search(free.begin(),free.end(),axis),"open RZ zero-area radial axis carries flux or free velocity");
    long double totalBoundary=0,totalPhysical=0,totalDissipation=0,totalBoundaryWork=0;
    for(int a=0;a<2;++a){const auto& c=result.report().components[a];const auto work=componentWork[a]+boundaryGross[a];
        scaled(componentMomentum[a],boundaryMomentum[a],momentumWork[a],"open boundary momentum incidence identity",1e-12);
        scaled(componentPower[a],outwardKinetic[a]+boundaryEnergy[a],work,"open physical KEP boundary identity",1e-12);
        scaled(boundaryEnergy[a],physicalEnergy[a]+inletDissipation[a],boundaryGross[a],"open physical supply plus inlet upwind dissipation",1e-12);
        scaled(c.boundaryMomentumRate,boundaryMomentum[a],momentumWork[a],"open retained component boundary momentum",1e-12);
        scaled(report.boundaryMomentumRate[a],boundaryMomentum[a],momentumWork[a],"open step boundary momentum",1e-12);
        scaled(c.boundaryKineticPower,boundaryEnergy[a],work,"open component boundary kinetic power",1e-12);
        scaled(c.boundaryPhysicalKineticPower,physicalEnergy[a],work,"open component physical kinetic supply",1e-12);
        scaled(c.inletDissipationPower,inletDissipation[a],work,"open component inlet upwind dissipation",1e-12);
        check(c.kineticIdentity<=1e-12&&c.momentumIdentity<=1e-12&&c.inletDissipationPower>=0,"open original KEP identity failed");
        totalBoundary+=boundaryEnergy[a];totalPhysical+=physicalEnergy[a];totalDissipation+=inletDissipation[a];totalBoundaryWork+=work;}
    scaled(report.boundaryKineticPower,totalBoundary,totalBoundaryWork,"open step boundary kinetic power",1e-12);
    scaled(report.boundaryPhysicalKineticPower,totalPhysical,totalBoundaryWork,"open step physical kinetic supply",1e-12);
    scaled(report.inletDissipationPower,totalDissipation,totalBoundaryWork,"open step inlet upwind dissipation",1e-12);
    powerWork+=std::abs(reactionPower);
    scaled(totalPower-reactionPower-tractionPower,0,powerWork,"open physical momentum power after reaction and traction T",1e-11);
    scaled(report.wallReactionPower,reactionPower,powerWork,"open actual prescribed-normal reaction power",1e-12);
    scaled(report.pressureOutletPower,tractionPower,powerWork,"open independently assembled pressure outlet power",1e-12);
    scaled(report.momentumPowerDefect,totalPower-reactionPower-tractionPower,powerWork,"open retained actual physical power defect",1e-12);
    scaled(report.momentumPowerRelativeDefect*powerWork,std::abs(totalPower-reactionPower-tractionPower),powerWork,"open physical power relative report",1e-12);
    check(report.momentumPowerRelativeDefect<=1e-11&&report.frozenOwnerMomentumDefect<=1e-11&&report.linearMomentumDefect<=1e-11,
        "open AC bypassed an original physical/free momentum gate");
    scaled(inertial-(ke1-ke0),cross,timeWork,"open actual temporal mass cross identity",1e-12);
    scaled(report.kineticBegin,ke0,ke0,"open begin kinetic energy",1e-12);scaled(report.kineticEnd,ke1,ke1,"open end kinetic energy",1e-12);
    scaled(report.inertialImpulseWork,inertial,timeWork,"open actual inertia work",1e-12);
    scaled(report.temporalMassCrossTerm,cross,timeWork,"open retained temporal cross term",1e-12);
    scaled(report.pressurePower,pressurePower,powerWork,"open actual pressure power",1e-12);
    scaled(report.capillaryPower,capillaryPower,powerWork,"open actual nonzero-sigma power",1e-12);
    const std::array<long double,3> area{{actualInterfaceMeasure(*phases[0]),actualInterfaceMeasure(*phases[1]),actualInterfaceMeasure(*phases[2])}};
    scaled(report.interfaceMeasureBegin,area[0],area[0],"open begin actual surface measure",1e-12);
    scaled(report.interfaceMeasureMidpoint,area[1],area[1],"open midpoint actual surface measure",1e-12);
    scaled(report.interfaceMeasureEnd,area[2],area[2],"open end actual surface measure",1e-12);
    std::size_t outletJets=0;const auto& scalar=material.beginGeometry().publication().state().geometry().snapshot();
    HydroCharacteristicConfig characteristic;characteristic.maxStep=.00025;
    for(const auto* phase:std::array<const HydroPhaseGeometry*,2>{{phases[1],phases[2]}})for(const auto& v:phase->publication().state().primary().vertices){
        const auto q=point(mesh,v.vertex);if(q.y!=mesh.maxCoord.y)continue;
        const auto trace=motion.traceBackward(q,phase->timeLayer().time(),step.beginState().timeLayer().time(),characteristic);
        check(!trace.boundary()&&trace.jacobian(),"open outlet did not use the retained original interior motion");
        const auto old=scalar.scalarJetNative(trace.point());const auto gradient=trace.jacobian()->transposeTimes({old.dx,old.dy});
        scaled(v.phi,old.value,std::abs(old.value),"open natural outlet original scalar",1e-12);
        scaled(v.dx,gradient.x,std::hypot(gradient.x,gradient.y),"open natural outlet original dx",1e-12);
        scaled(v.dy,gradient.y,std::hypot(gradient.x,gradient.y),"open natural outlet original dy",1e-12);++outletJets;}
    check(outletJets>0,"open step never observed natural outlet FirstJets");
    const auto& saddle=step.solveResult().report();
    check(saddle.algebraicRelativeResidual<=1e-11&&saddle.freeMomentumDefect<=1e-11&&saddle.maxCellContinuityDefect<=1e-11
        &&saddle.maxPressureMeanDefect<=1e-11&&saddle.maxBoundaryFluxDefect<=1e-12,"open original saddle gate changed");
    check(report.iterations==(report.corrected?2:1)&&report.history.size()==static_cast<std::size_t>(report.iterations),"open AC stage history changed");
    std::cout<<(pressureOutlet?"PRESSURE_OUTLET_STEP,":"PRESCRIBED_FLOW_STEP,")<<step.candidateState().timeLayer().step<<','<<phaseMovement<<','<<massMovement<<','
        <<incoming<<','<<outgoing<<','<<totalBoundary<<','<<totalPhysical<<','<<totalDissipation<<','<<reactionPower
        <<','<<tractionPower<<','<<outletChange<<','<<pressureMean<<','<<step.saddle().pressureGaugeCount()
        <<','<<scaled(equationGap,0,equationWork,"reported actual free momentum",1e-11)
        <<','<<scaled(totalPower-reactionPower-tractionPower,0,powerWork,"reported actual power",1e-11)<<",PASS\n";
}
void prescribedFlowOnly(bool pressureOutlet=false)
{
    auto mesh=makeMesh(4);Case shape;shape.pureFirstJet=true;shape.radial=true;shape.radius=.27;shape.cy=.51;shape.motion=1;shape.speed=.2;
    const auto profile=[shape](double time)->Profile{return [shape,time](Vector2D q){const auto j=shape.profile(time)(q);return Jet{-j.value,-j.dx,-j.dy,-j.dxy};};};
    const auto f=profile(0);const auto seed=stateFor(*mesh,f,true,{},[](int,Vector2D,Vector2D){return 0.;},{},HermiteGalsTransport::FirstJet);
    HermiteGalsCandidate candidate(*mesh,{},seed.primary());const auto mixed=closeHermiteGalsMixedAmr(*mesh,candidate);
    const auto initialPub=publishHermiteGalsFirstJetGeometry(*mesh,std::move(candidate),mixed,{candidate.identity(),centers(*mesh,f)});
    const auto geometry=HydroPhaseGeometry::build(*mesh,stepMetrics(*mesh,true),initialPub);const auto& l=geometry.layout();
    inspectRzPublication(*mesh,geometry);std::vector<double> normal(l.velocityCount());HydroTwoPhasePrescribedFlow flow;
    flow.normalVelocity.resize(l.velocityCount());
    for(int u=0;u<l.velocityCount();++u){const auto face=l.domain().face(*mesh,l.nativeFace(u));if(face.horizontal)normal[u]=.2;
        if(face.horizontal&&face.outwardSign!=0){
            if(pressureOutlet&&face.outwardSign>0)flow.pressureOutlets.push_back({u,2.5+.025*face.center.x*face.center.x});
            else {flow.normalVelocity[u]=.2;(face.outwardSign<0?flow.inletVelocityDofs:flow.outletVelocityDofs).push_back(u);}}}
    const HydroNativeValues values(l,normal,std::vector<double>(l.pressureCount()));std::size_t inletJets=0,inletVelocities=0;
    const auto inletPoint=[&](Vector2D q,HermiteGalsTimeLayer t,Boundary side){
        check(q.y==mesh->minCoord.y&&(side==Boundary::south||side==Boundary::southWest||side==Boundary::southEast),
            "prescribed inlet callback queried an outlet, wall, or interior point");
        check(t.step>=1&&t.step<=2&&t.end>t.begin&&t.time()>t.begin&&t.time()<=t.end,"prescribed inlet callback lost actual arrival full clock");};
    flow.inletVelocity=[&](Vector2D q,HermiteGalsTimeLayer t,Boundary side){inletPoint(q,t,side);++inletVelocities;return Vector2D{0,.2};};
    flow.inletFirstJet=[&](Vector2D q,HermiteGalsTimeLayer t,Boundary side){inletPoint(q,t,side);++inletJets;
        const auto j=profile(t.time())(q);return HermiteGalsBoundaryFirstJet{j.value,{j.dx,j.dy}};};
    auto accepted=HydroTwoPhaseState::createPrescribedFlow(*mesh,geometry,values,flow);const auto initialHandle=accepted;
    const auto initialCopy=capture(accepted.publication(),accepted.values());check(accepted.prescribedFlow(),"open state forgot its immutable flow data");
    HydroTwoPhaseStepConfig config;config.liquidDensity=1000;config.gasDensity=1;config.surfaceTension=1;
    config.liquidViscosity=0;config.gasViscosity=0;config.spaceTime.initialPanels=1;config.spaceTime.characteristic.maxStep=.00025;
    config.spaceTime.firstJetBoundaryTreatment=HermiteGalsBoundaryTreatment::InflowHistory;
    check(config.timeScheme==HydroTwoPhaseTimeScheme::AdaptiveLinearMidpoint&&config.interfaceMismatchTolerance==.01
        &&config.momentumPowerTolerance==1e-11,"open fixture changed original AC or power targets");
    const HermiteGalsTimeLayer first{1,0,.00025,Position::End};int rejected=0;
    const auto refuse=[&](const std::function<void()>& action,const char* why){rejects(action,why);++rejected;
        check(accepted.sameIdentity(initialHandle),"open refusal changed accepted handle");unchanged(initialCopy,accepted.publication(),accepted.values());};
    auto bad=flow;
    if(pressureOutlet){
        bad.normalVelocity[flow.pressureOutlets.front().velocityDof]=.3;
        refuse([&]{HydroTwoPhaseState::createPrescribedFlow(*mesh,geometry,values,bad);},"pressure outlet prescribed normal payload");
        bad=flow;bad.pressureOutlets.push_back(flow.pressureOutlets.front());
        refuse([&]{HydroTwoPhaseState::createPrescribedFlow(*mesh,geometry,values,bad);},"duplicate, overlapping or invalid pressure outlet");
        bad=flow;bad.pressureOutlets.front().velocityDof=flow.inletVelocityDofs.front();
        refuse([&]{HydroTwoPhaseState::createPrescribedFlow(*mesh,geometry,values,bad);},"duplicate, overlapping or invalid pressure outlet");
        bad=flow;bad.pressureOutlets.front().pressure=std::numeric_limits<double>::quiet_NaN();
        refuse([&]{HydroTwoPhaseState::createPrescribedFlow(*mesh,geometry,values,bad);},"nonfinite pressure outlet");
        bad=flow;bad.pressureOutlets.front().velocityDof=geometry.metrics().axisVelocityDofs().front();
        refuse([&]{HydroTwoPhaseState::createPrescribedFlow(*mesh,geometry,values,bad);},"internal or zero-area axis");
        bad=flow;for(int u=0;u<l.velocityCount();++u){const auto face=l.domain().face(*mesh,l.nativeFace(u));
            if(face.outwardSign==0){bad.pressureOutlets.front().velocityDof=u;break;}}
        refuse([&]{HydroTwoPhaseState::createPrescribedFlow(*mesh,geometry,values,bad);},"internal or zero-area axis");
    }else{
        bad.normalVelocity[flow.outletVelocityDofs.front()]=.3;
        refuse([&]{HydroTwoPhaseState::createPrescribedFlow(*mesh,geometry,values,bad);},"accepted prescribed-flow normal mean differs");
    }
    bad=flow;bad.inletVelocityDofs.push_back(geometry.metrics().axisVelocityDofs().front());
    refuse([&]{HydroTwoPhaseState::createPrescribedFlow(*mesh,geometry,values,bad);},"zero-area axis");
    bad=flow;bad.inletVelocityDofs.push_back(flow.inletVelocityDofs.front());
    refuse([&]{HydroTwoPhaseState::createPrescribedFlow(*mesh,geometry,values,bad);},"duplicate or invalid prescribed-flow port");
    auto wrong=config;wrong.spaceTime.firstJetBoundaryTreatment=HermiteGalsBoundaryTreatment::PrescribedOuterFirstJet;
    refuse([&]{HydroTwoPhaseStep::prepare(*mesh,accepted,first,wrong);},"explicit inlet-only FirstJet boundary treatment");
    auto failing=flow;failing.inletFirstJet=[](Vector2D,HermiteGalsTimeLayer,Boundary)->HermiteGalsBoundaryFirstJet{
        throw std::runtime_error("PRESCRIBED_INLET_FIRSTJET_THROW");};
    auto failingState=HydroTwoPhaseState::createPrescribedFlow(*mesh,geometry,values,failing);const auto failingHandle=failingState;
    const auto failingCopy=capture(failingState.publication(),failingState.values());
    refuse([&]{HydroTwoPhaseStep::prepare(*mesh,failingState,first,config);},"PRESCRIBED_INLET_FIRSTJET_THROW");
    check(failingState.sameIdentity(failingHandle),"inlet callback failure changed its actual accepted handle");
    unchanged(failingCopy,failingState.publication(),failingState.values());
    auto late=config;int lateCalls=0;late.onIteration=[&](const HydroTwoPhaseIterationReport& r){++lateCalls;
        check(r.geometrySamples>=2&&r.acceptedPanels>0&&r.saddle.freeVelocityCount>0,"open late failure missed actual material/saddle work");
        unchanged(initialCopy,accepted.publication(),accepted.values());throw std::runtime_error("PRESCRIBED_COMPLETED_STAGE_THROW");};
    refuse([&]{HydroTwoPhaseStep::prepare(*mesh,accepted,first,late);},"PRESCRIBED_COMPLETED_STAGE_THROW");
    check(lateCalls==1,"open late observer exception was swallowed or repeated");
    for(int k=1;k<=2;++k){const auto retained=accepted;const auto saved=capture(accepted.publication(),accepted.values());
        const HermiteGalsTimeLayer clock{k,(k-1)*.00025,k*.00025,Position::End};auto observedConfig=config;
        std::vector<HydroTwoPhaseIterationReport> observed;observedConfig.onIteration=[&](const HydroTwoPhaseIterationReport& r){
            unchanged(saved,accepted.publication(),accepted.values());observed.push_back(r);};
        const auto step=HydroTwoPhaseStep::prepare(*mesh,accepted,clock,observedConfig);
        unchanged(saved,accepted.publication(),accepted.values());inspectIterationObserver(step,observed);inspectPrescribedFlowStep(*mesh,step);
        check(step.beginState().sameIdentity(retained)&&sameClock(step.candidateState().timeLayer(),clock)
            &&step.report().usedVelocityHistory==(k==2)&&step.report().startupCorrection==(k==1),"open state/AC history detached actual accepted source");
        step.commit(*mesh,accepted);check(accepted.sameIdentity(step.candidateState())&&!accepted.sameIdentity(retained),"open successful step did not commit actual candidate");
        unchanged(saved,retained.publication(),retained.values());const auto committed=capture(accepted.publication(),accepted.values());
        rejects([&]{step.commit(*mesh,accepted);},"stale candidate");++rejected;unchanged(committed,accepted.publication(),accepted.values());
    }
    check(inletJets>0&&inletVelocities>0&&rejected==(pressureOutlet?13:8),"open inlet/commit/refusal coverage incomplete");
    std::cout<<(pressureOutlet?"PRESSURE_OUTLET_ONLY,PASS":"PRESCRIBED_FLOW_ONLY,PASS")<<" steps=2 rejected="<<rejected
        <<" inlet_jets="<<inletJets<<" inlet_velocities="<<inletVelocities<<" rho=1000:1 mu=0 sigma=1"<<'\n';
}

void rzOnly()
{
    auto mesh=makeMesh(4);Case shape;shape.pureFirstJet=true;shape.radial=true;shape.radius=.27;shape.cy=.51;
    const auto geometry=HydroPhaseGeometry::build(*mesh,stepMetrics(*mesh,true),initial(*mesh,shape));
    inspectRzPublication(*mesh,geometry,true);
    auto accepted=HydroTwoPhaseState::create(*mesh,geometry,rzStepVelocity(*mesh,geometry.layout()));
    nativeConstraints(*mesh,accepted.values());
    auto badAxis=accepted.values().normalVelocity();badAxis[geometry.metrics().axisVelocityDofs().front()]=1;
    bool axisRejected=false;try{HydroTwoPhaseState::create(*mesh,geometry,HydroNativeValues(geometry.layout(),badAxis,accepted.values().pressure()));}
    catch(const std::exception&){axisRejected=true;}check(axisRejected,"nonzero actual radial axis velocity accepted");
    HydroTwoPhaseStepConfig config;config.liquidDensity=1000;config.gasDensity=1;config.surfaceTension=1;
    config.liquidViscosity=0;config.gasViscosity=0;config.spaceTime.initialPanels=1;config.spaceTime.characteristic.maxStep=.00025;
    check(config.timeScheme==HydroTwoPhaseTimeScheme::AdaptiveLinearMidpoint&&config.interfaceMismatchTolerance==.01
        &&config.momentumPowerTolerance==1e-11,"RZ full step changed default AC sensor or physical power target");
    const auto initialHandle=accepted;const auto initialCopy=capture(accepted.publication(),accepted.values());
    auto rejectedConfig=config;int lateCallbacks=0;
    rejectedConfig.onIteration=[&](const HydroTwoPhaseIterationReport& r){++lateCallbacks;
        check(r.geometrySamples>=2&&r.acceptedPanels>0&&r.saddle.freeVelocityCount>0,"RZ failure did not reach actual completed material/saddle stage");
        unchanged(initialCopy,accepted.publication(),accepted.values());throw std::runtime_error("RZ_COMPLETED_STAGE_OBSERVER_THROW");};
    rejects([&]{HydroTwoPhaseStep::prepare(*mesh,accepted,{1,0,.00025,Position::End},rejectedConfig);},"RZ_COMPLETED_STAGE_OBSERVER_THROW");
    check(lateCallbacks==1&&accepted.sameIdentity(initialHandle),"RZ late failure changed accepted handle or repeated callback");
    unchanged(initialCopy,accepted.publication(),accepted.values());
    for(int k=1;k<=2;++k){const auto retained=accepted;const auto saved=capture(accepted.publication(),accepted.values());
        const HermiteGalsTimeLayer clock{k,(k-1)*.00025,k*.00025,Position::End},mid{k,clock.begin,clock.end,Position::Midpoint};
        auto observedConfig=config;std::vector<HydroTwoPhaseIterationReport> observed;
        observedConfig.onIteration=[&](const HydroTwoPhaseIterationReport& r){unchanged(saved,accepted.publication(),accepted.values());observed.push_back(r);};
        const auto step=HydroTwoPhaseStep::prepare(*mesh,accepted,clock,observedConfig);
        unchanged(saved,accepted.publication(),accepted.values());inspectIterationObserver(step,observed);
        inspectStep(*mesh,step,true);const auto& material=step.materialTransport();
        check(step.beginState().sameIdentity(retained)&&sameClock(material.beginGeometry().timeLayer(),retained.timeLayer())
            &&sameClock(material.midpointRate().geometry().timeLayer(),mid)&&sameClock(material.endGeometry().timeLayer(),clock)
            &&sameClock(step.candidateState().timeLayer(),clock),"RZ actual three-layer source/full clock mismatch");
        const auto& motion=material.characteristicSource();
        check(motion.sameIdentity(step.characteristicSource())&&motion.previous().geometryIdentity()==retained.geometryIdentity()
            &&motion.arrival().geometryIdentity()==retained.geometryIdentity()&&sameClock(motion.previous().timeLayer(),retained.timeLayer())
            &&sameClock(motion.arrival().timeLayer(),clock),"RZ material detached its actual accepted characteristic source/clock");
        const auto& dual=step.dualMaterial();
        check(sameClock(dual.beginMass().timeLayer(),retained.timeLayer())&&sameClock(dual.midpointMass().timeLayer(),mid)
            &&sameClock(dual.endMass().timeLayer(),clock),"RZ true three-mass family lost actual full clocks");
        inspectRzPublication(*mesh,material.beginGeometry());inspectRzPublication(*mesh,material.midpointRate().geometry());inspectRzPublication(*mesh,material.endGeometry());
        for(int axis:geometry.metrics().axisVelocityDofs()){
            check(material.faceIntegrals()[axis]==0&&step.candidateState().values().normalVelocity()[axis]==0,
                "RZ axis carries physical material flux or velocity");
            const auto face=geometry.layout().domain().face(*mesh,geometry.layout().nativeFace(axis));
            const int owner=geometry.layout().pressureDof(face.positiveCell>=0?face.positiveCell:face.negativeCell);
            const auto center=mesh->cells(mesh->cellIndexToLabel[geometry.layout().nativeCell(owner)]).centerCoord();
            for(const auto* endpoint:std::array<const HydroSpatialAdvector*,2>{{&motion.previous(),&motion.arrival()}}){
                const auto jet=endpoint->cellVelocityJet(owner,{0,center.y}).jet;
                check(jet.value.x==0&&jet.derivative.xy==0&&jet.derivative.yx==0,
                    "RZ actual endpoint source lost radial odd/axial even axis constraints");}
        }
        const auto& force=step.capillaryLoad().momentumRhsLoad();long double forceNorm=0,work=0;
        for(std::size_t i=0;i<force.size();++i){forceNorm=std::hypot(forceNorm,force[i]);work+=step.solvedMidpointValues().normalVelocity()[i]*force[i];}
        check(forceNorm>0&&std::isfinite(work)&&work!=0,"RZ moving nonzero-sigma step erased actual capillary force/work");
        check(step.report().usedVelocityHistory==(k==2)&&step.report().startupCorrection==(k==1),"RZ AC lost startup/accepted velocity history");
        step.commit(*mesh,accepted);check(accepted.sameIdentity(step.candidateState())&&!accepted.sameIdentity(retained)
            &&sameClock(accepted.timeLayer(),clock),"RZ successful complete step did not atomically commit");
        unchanged(saved,retained.publication(),retained.values());const auto committed=capture(accepted.publication(),accepted.values());
        rejects([&]{step.commit(*mesh,accepted);},"stale candidate");unchanged(committed,accepted.publication(),accepted.values());
        std::cout<<"RZ_TWO_PHASE_STEP,"<<k<<",rho=1000:1,mu=0,sigma=1,"<<step.report().iterations<<','<<forceNorm<<','<<work<<",PASS\n";
    }
    unchanged(initialCopy,initialHandle.publication(),initialHandle.values());
    check(movingCases==0&&acceptedSteps==0&&guardCases==0&&capillarySteps==0&&capillaryRejects==0
        &&linearCases==0&&linearRejects==0,"RZ short entry changed original XY inventories");
    std::cout<<"RZ_ONLY,PASS steps=2 axis_rejected=1 atomic_prepare=1 stale_commit=2 equation="<<peakEquation
        <<" energy="<<peakEnergy<<" material="<<peakMaterial<<'\n';
}

// Constant acceleration is an external load, never a replacement mass or flux.
void inspectBodyAcceleration(const QuadtreeMesh& mesh,const HydroTwoPhaseStep& step,Vector2D g)
{
    inspectStep(mesh,step,true);const auto& dual=step.dualMaterial();const auto& material=step.materialTransport();
    const auto& l=dual.layout();const auto& bar=step.solvedMidpointValues().normalVelocity();const auto& load=step.bodyLoads();
    std::array<long double,2> norms{};long double work=0,gross=0;
    for(int u=0;u<l.velocityCount();++u){const auto face=l.domain().face(mesh,l.nativeFace(u));const int axis=face.horizontal?1:0;
        const auto phases=dual.midpointMass().phaseVolumes()[u];
        const long double mass=material.liquidDensity()*phases.liquid+material.gasDensity()*phases.gas;
        scaled(dual.midpointMass().masses()[u],mass,mass,"body independent midpoint sharp mass",1e-12);
        const long double expected=mass*(axis==0?g.x:g.y);
        scaled(load[u],expected,std::abs(expected),"body actual midpoint phase mass times physical axis",1e-12);
        norms[axis]=std::hypot(norms[axis],load[u]);work+=static_cast<long double>(bar[u])*load[u];gross+=std::abs(static_cast<long double>(bar[u])*load[u]);
        if(l.domain().axisymmetric()&&axis==0)check(load[u]==0,"axial RZ gravity acquired a radial load");
    }
    check((g.x==0||norms[0]>0)&&(g.y==0||norms[1]>0)&&std::isfinite(work)&&work!=0,
        "moving contrast gravity fixture erased its actual component or external work");
    scaled(step.report().bodyPower,work,gross,"body solved-midpoint external work",1e-12);
    const auto& r=step.report();const long double dt=r.deltaTime;
    const long double balance=r.inertialImpulseWork+dt*(r.pressurePower+r.viscousPower+r.convectionPower-r.consistencyPower-r.capillaryPower-r.bodyPower);
    const long double scale=std::abs(r.inertialImpulseWork)+dt*(std::abs(r.pressurePower)+std::abs(r.viscousPower)+std::abs(r.convectionPower)+std::abs(r.consistencyPower)+std::abs(r.capillaryPower)+gross);
    scaled(balance,0,scale,"body original kinetic ledger includes external work",1e-11);
    std::cout<<"BODY_ACCELERATION_STEP,"<<(l.domain().axisymmetric()?"RZ":"XY")<<",rho=1000:1,mu=0,sigma=1,"
        <<g.x<<','<<g.y<<','<<norms[0]<<','<<norms[1]<<','<<work<<",PASS\n";
}
void bodyAccelerationOnly()
{
    auto mesh=makeMesh(4);Case shape;shape.pureFirstJet=true;
    const auto geometry=HydroPhaseGeometry::build(*mesh,stepMetrics(*mesh),initial(*mesh,shape));
    auto accepted=HydroTwoPhaseState::create(*mesh,geometry,stepVelocity(*mesh,geometry.layout()));
    const auto initialHandle=accepted;const auto initialCopy=capture(accepted.publication(),accepted.values());
    HydroTwoPhaseStepConfig config;config.liquidDensity=1000;config.gasDensity=1;config.surfaceTension=1;
    config.spaceTime.initialPanels=1;config.spaceTime.characteristic.maxStep=.00025;
    const HermiteGalsTimeLayer clock{1,0,.00025,Position::End};int rejected=0;
    const auto refuse=[&](const std::function<void()>& action,const HydroTwoPhaseState& state,const char* why){
        const auto handle=state;const auto saved=capture(state.publication(),state.values());rejects(action,why);++rejected;
        check(state.sameIdentity(handle),"body failure changed the actual accepted handle");unchanged(saved,state.publication(),state.values());};
    auto bad=config;bad.bodyAcceleration.x=std::numeric_limits<double>::quiet_NaN();
    refuse([&]{HydroTwoPhaseStep::prepare(*mesh,accepted,clock,bad);},accepted,"invalid body acceleration");
    bad=config;bad.bodyAcceleration.y=std::numeric_limits<double>::infinity();
    refuse([&]{HydroTwoPhaseStep::prepare(*mesh,accepted,clock,bad);},accepted,"invalid body acceleration");
    bad=config;bad.bodyAcceleration.x=std::numeric_limits<double>::max();
    refuse([&]{HydroTwoPhaseStep::prepare(*mesh,accepted,clock,bad);},accepted,"unrepresentable body load");
    const auto implicitZero=HydroTwoPhaseStep::prepare(*mesh,accepted,clock,config);
    auto explicitConfig=config;explicitConfig.bodyAcceleration={0,0};
    const auto explicitZero=HydroTwoPhaseStep::prepare(*mesh,accepted,clock,explicitConfig);
    inspectStep(*mesh,implicitZero,true);
    check(implicitZero.integralMomentumRhs()==explicitZero.integralMomentumRhs()
        &&implicitZero.candidateState().values().normalVelocity()==explicitZero.candidateState().values().normalVelocity()
        &&implicitZero.candidateState().values().pressure()==explicitZero.candidateState().values().pressure()
        &&implicitZero.report().iterations==explicitZero.report().iterations,"explicit zero body acceleration changed original arithmetic");
    const Matrix zeroMatrixGap=implicitZero.velocityMatrix()-explicitZero.velocityMatrix();check(zeroMatrixGap.norm()==0,"zero body acceleration changed the original matrix");
    check(implicitZero.report().bodyPower==0&&explicitZero.report().bodyPower==0
        &&std::all_of(implicitZero.bodyLoads().begin(),implicitZero.bodyLoads().end(),[](long double f){return f==0;})
        &&implicitZero.bodyLoads()==explicitZero.bodyLoads(),"zero body acceleration is not exactly zero");
    config.bodyAcceleration={.4,-.7};auto late=config;int lateCalls=0;
    late.onIteration=[&](const HydroTwoPhaseIterationReport& r){++lateCalls;check(r.geometrySamples>=2&&r.acceptedPanels>0&&r.saddle.freeVelocityCount>0,
        "body late failure did not reach the actual completed material/saddle stage");throw std::runtime_error("BODY_COMPLETED_STAGE_THROW");};
    refuse([&]{HydroTwoPhaseStep::prepare(*mesh,accepted,clock,late);},accepted,"BODY_COMPLETED_STAGE_THROW");check(lateCalls==1,"body late exception was swallowed or repeated");
    const auto xy=HydroTwoPhaseStep::prepare(*mesh,accepted,clock,config);inspectBodyAcceleration(*mesh,xy,config.bodyAcceleration);
    check(accepted.sameIdentity(initialHandle),"body successful prepare changed accepted handle");unchanged(initialCopy,accepted.publication(),accepted.values());
    xy.commit(*mesh,accepted);check(accepted.sameIdentity(xy.candidateState())&&sameClock(accepted.timeLayer(),clock),"body commit lost actual complete candidate");
    refuse([&]{xy.commit(*mesh,accepted);},accepted,"stale candidate");unchanged(initialCopy,initialHandle.publication(),initialHandle.values());
    auto rzMesh=makeMesh(4);Case sphere;sphere.pureFirstJet=true;sphere.radial=true;sphere.radius=.27;sphere.cy=.51;
    const auto rzGeometry=HydroPhaseGeometry::build(*rzMesh,stepMetrics(*rzMesh,true),initial(*rzMesh,sphere));
    auto rz=HydroTwoPhaseState::create(*rzMesh,rzGeometry,rzStepVelocity(*rzMesh,rzGeometry.layout()));
    const auto rzHandle=rz;const auto rzCopy=capture(rz.publication(),rz.values());auto axial=config;axial.bodyAcceleration={0,-.7};
    auto radial=axial;radial.bodyAcceleration.x=.1;
    refuse([&]{HydroTwoPhaseStep::prepare(*rzMesh,rz,clock,radial);},rz,"RZ body acceleration must be axial");
    const auto z=HydroTwoPhaseStep::prepare(*rzMesh,rz,clock,axial);inspectBodyAcceleration(*rzMesh,z,axial.bodyAcceleration);
    unchanged(rzCopy,rz.publication(),rz.values());check(rz.sameIdentity(rzHandle),"RZ body prepare changed accepted handle");
    for(int u:rzGeometry.metrics().axisVelocityDofs())check(z.bodyLoads()[u]==0&&z.candidateState().values().normalVelocity()[u]==0,
        "RZ axial gravity lost exact zero radial axis load/velocity");
    z.commit(*rzMesh,rz);check(rz.sameIdentity(z.candidateState())&&sameClock(rz.timeLayer(),clock),"RZ body step did not atomically commit");
    unchanged(rzCopy,rzHandle.publication(),rzHandle.values());check(rejected==6,"body refusal/atomic coverage incomplete");
    std::cout<<"BODY_ACCELERATION_ONLY,PASS forced_steps=2 zero_preparations=2 rejected="<<rejected<<" exact_zero_path=1\n";
}

// Side-wall reflection constrains only the normal FirstJet derivative. The
// scalar on an impermeable wall still moves under its actual tangential flow.
void symmetryBoundaryOnly()
{
    DictionaryMesh dictionary;dictionary.Nx=4;dictionary.Ny=16;
    dictionary.minCoord={0,0};dictionary.maxCoord={1,4};dictionary.maxInitialIterations=0;dictionary.maxCellNum=16384;
    QuadtreeMesh mesh(dictionary);
    const Profile profile=[](Vector2D q){const double k=2*double(pi),a=.02;
        const double dx=q.x==0||q.x==1?0:-a*k*std::sin(k*q.x);
        return Jet{2+a*std::cos(k*q.x)-q.y,dx,-1,0};};
    HermiteGalsPublicationConfig publication;publication.volume.symmetryPlanes={{0,0},{0,1}};
    const auto seed=stateFor(mesh,profile,false,{},[](int,Vector2D,Vector2D){return 0.;},
        publication.volume.symmetryPlanes,HermiteGalsTransport::FirstJet);
    HermiteGalsCandidate candidate(mesh,{},seed.primary());
    const auto mixed=closeHermiteGalsMixedAmr(mesh,candidate,publication.volume.symmetryPlanes);
    const auto identity=candidate.identity();
    const auto initialPublication=publishHermiteGalsFirstJetGeometry(mesh,std::move(candidate),mixed,{identity,centers(mesh,profile)},publication);
    const auto geometry=HydroPhaseGeometry::build(mesh,stepMetrics(mesh),initialPublication);
    auto accepted=HydroTwoPhaseState::create(mesh,geometry,HydroNativeValues(geometry.layout(),
        std::vector<double>(geometry.layout().velocityCount()),std::vector<double>(geometry.layout().pressureCount())));
    const auto initialHandle=accepted;const auto initialCopy=capture(accepted.publication(),accepted.values());
    HydroTwoPhaseStepConfig config;config.liquidDensity=3;config.gasDensity=1;config.surfaceTension=0;
    config.bodyAcceleration={0,-1};config.publication=publication;config.spaceTime.initialPanels=1;
    config.spaceTime.characteristic.maxStep=.001;
    config.spaceTime.firstJetBoundaryTreatment=HermiteGalsBoundaryTreatment::SymmetryFirstJet;
    const auto side=[](Vector2D q){return q.x==0||q.x==1;};
    const auto inspectWall=[&](const HermiteGalsPublication& pub){
        std::array<int,2> vertices{},edges{};const auto& primary=pub.state().primary();
        for(const auto& v:primary.vertices){const auto q=point(mesh,v.vertex);if(!side(q))continue;
            check(v.dx==0,"symmetry wall primary vertex normal derivative is not exact zero");
            check(std::isfinite(v.phi)&&std::isfinite(v.dy),"symmetry wall vertex lost its finite scalar/tangent jet");++vertices[q.x==0?0:1];}
        for(const auto& e:primary.edges){const auto a=point(mesh,e.vertices[0]),b=point(mesh,e.vertices[1]);
            if(a.x!=b.x||!side(a))continue;
            check(e.dn==0,"symmetry wall primary edge normal derivative is not exact zero");
            check(std::isfinite(e.phi),"symmetry wall edge scalar is nonfinite");++edges[a.x==0?0:1];}
        check(vertices==std::array<int,2>{{17,17}}&&edges==std::array<int,2>{{16,16}},"symmetry wall did not inspect every actual side vertex/edge");
    };
    inspectWall(accepted.publication());int oldRejected=0,twoStages=0;long double wallMovement=0,totalBodyWork=0;
    // Historical InflowHistory keeps its original semantics. A coincidental
    // exact zero on this coarse mesh is recorded, never turned into a fallback.
    auto historical=config;historical.spaceTime.firstJetBoundaryTreatment=HermiteGalsBoundaryTreatment::InflowHistory;
    try{const auto old=HydroTwoPhaseStep::prepare(mesh,accepted,{1,0,.001,Position::End},historical);(void)old;
        std::cout<<"SYMMETRY_OLD_INFLOW,exact_coincidence_accepted\n";}
    catch(const std::exception& e){check(std::string(e.what()).find("symmetry requires exactly zero primary")!=std::string::npos,
        std::string("historical Inflow refusal reached an unrelated stage: ")+e.what());++oldRejected;std::cout<<"SYMMETRY_OLD_INFLOW,original_rejection,"<<e.what()<<'\n';}
    check(accepted.sameIdentity(initialHandle),"historical Inflow attempt changed accepted handle");unchanged(initialCopy,accepted.publication(),accepted.values());
    for(int k=1;k<=2;++k){const auto retained=accepted;const auto saved=capture(accepted.publication(),accepted.values());
        const HermiteGalsTimeLayer clock{k,(k-1)*.001,k*.001,Position::End};
        const auto step=HydroTwoPhaseStep::prepare(mesh,accepted,clock,config);const auto& r=step.report();
        unchanged(saved,accepted.publication(),accepted.values());check(accepted.sameIdentity(retained),"symmetry prepare changed accepted handle");
        check(r.timeScheme==HydroTwoPhaseTimeScheme::AdaptiveLinearMidpoint&&r.iterations==(r.corrected?2:1)
            &&r.linearMomentumDefect<=1e-11&&r.frozenOwnerMomentumDefect<=1e-11&&r.momentumPowerRelativeDefect<=1e-11,
            "symmetry step changed original AC/physical momentum gates");
        const auto& saddle=step.solveResult().report();
        check(saddle.algebraicRelativeResidual<=1e-11&&saddle.freeMomentumDefect<=1e-11&&saddle.maxCellContinuityDefect<=1e-11
            &&saddle.maxPressureMeanDefect<=1e-11&&saddle.maxBoundaryFluxDefect<=1e-12,"symmetry changed original saddle gates");
        check(step.materialTransport().spaceTimeReport().firstJetBoundaryTreatment==HermiteGalsBoundaryTreatment::SymmetryFirstJet,
            "symmetry material transport lost the explicitly selected boundary treatment");
        if(r.iterations==2)++twoStages;
        if(k==1)check(r.startupCorrection&&r.iterations==2,"symmetry gravity startup did not execute actual two-stage AC");
        inspectWall(step.materialTransport().midpointRate().geometry().publication());inspectWall(step.candidateState().publication());
        nativeConstraints(mesh,step.candidateState().values());
        long double volume=0;for(const auto& v:step.candidateState().phaseGeometry().cellVolumes()){
            check(std::isfinite(v.liquid)&&std::isfinite(v.gas)&&v.liquid>=0&&v.gas>=0,"symmetry actual phase volume is invalid");volume+=v.liquid;}
        check(volume>0&&volume<4&&std::isfinite(volume)&&std::isfinite(r.bodyPower)&&std::isfinite(r.kineticEnd),"symmetry volume/body work/kinetic energy is invalid");
        long double power=0,gross=0;const auto& bar=step.solvedMidpointValues().normalVelocity();
        for(std::size_t u=0;u<bar.size();++u){const long double work=static_cast<long double>(bar[u])*step.bodyLoads()[u];power+=work;gross+=std::abs(work);}
        scaled(r.bodyPower,power,gross,"symmetry actual retained midpoint external body work",1e-12);totalBodyWork+=r.deltaTime*r.bodyPower;
        const auto& begin=step.beginState().publication().state().primary();const auto& end=step.candidateState().publication().state().primary();
        check(begin.vertices.size()==end.vertices.size()&&begin.edges.size()==end.edges.size(),"symmetry changed primary inventory");
        for(std::size_t i=0;i<end.vertices.size();++i)if(side(point(mesh,end.vertices[i].vertex)))wallMovement=std::max(wallMovement,std::abs(static_cast<long double>(end.vertices[i].phi)-begin.vertices[i].phi));
        for(std::size_t i=0;i<end.edges.size();++i){const auto a=point(mesh,end.edges[i].vertices[0]),b=point(mesh,end.edges[i].vertices[1]);
            if(a.x==b.x&&side(a))wallMovement=std::max(wallMovement,std::abs(static_cast<long double>(end.edges[i].phi)-begin.edges[i].phi));}
        step.commit(mesh,accepted);check(accepted.sameIdentity(step.candidateState())&&sameClock(accepted.timeLayer(),clock),"symmetry candidate did not atomically commit");
        unchanged(saved,retained.publication(),retained.values());
        std::cout<<"SYMMETRY_BOUNDARY_STEP,"<<k<<",rho=3:1,g_y=-1,sigma=0,"<<r.iterations<<','<<volume<<','<<r.bodyPower<<','<<wallMovement<<",PASS\n";
    }
    check(twoStages>0&&std::isfinite(wallMovement)&&wallMovement>0&&std::isfinite(totalBodyWork)&&totalBodyWork!=0
        &&accepted.values().normalVelocity()!=initialHandle.values().normalVelocity(),"symmetry gravity fixture froze wall scalar or actual motion/work");
    unchanged(initialCopy,initialHandle.publication(),initialHandle.values());
    std::cout<<"SYMMETRY_BOUNDARY_ONLY,PASS steps=2 side_vertices=34 side_edges=32 exact_normal_jet=1 wall_phi_movement="
        <<wallMovement<<" body_work="<<totalBodyWork<<" old_inflow_rejected="<<oldRejected<<'\n';
}

void guards(const QuadtreeMesh& mesh,const HydroPhaseGeometry& geometry)
{
    auto accepted=HydroTwoPhaseState::create(mesh,geometry,stepVelocity(mesh,geometry.layout()));const auto retained=accepted;
    const auto saved=capture(accepted.publication(),accepted.values());const HermiteGalsTimeLayer clock{1,0,.001,Position::End};
    HydroTwoPhaseStepConfig config;config.timeScheme=HydroTwoPhaseTimeScheme::FullyCoupledMidpoint;
    config.liquidDensity=1e3;config.spaceTime.characteristic.maxStep=.00025;config.spaceTime.initialPanels=1;
    const auto fail=[&](const HydroTwoPhaseStepConfig& c,HermiteGalsTimeLayer t,const char* expected=nullptr){
        try{stepRejects([&]{HydroTwoPhaseStep::prepare(mesh,accepted,t,c);},expected);}
        catch(...){check(accepted.sameIdentity(retained),"wrong-stage failure changed complete accepted handle");
            unchanged(saved,accepted.publication(),accepted.values());
            std::cout<<"STEP_WRONG_STAGE_ATOMIC_STATE,PASS\n";throw;}
        check(accepted.sameIdentity(retained),"failure changed complete accepted handle");unchanged(saved,accepted.publication(),accepted.values());};
    for(double bad:{0.,-1.,std::numeric_limits<double>::infinity(),std::numeric_limits<double>::quiet_NaN()}){
        auto c=config;c.liquidDensity=bad;fail(c,clock);c=config;c.gasDensity=bad;fail(c,clock);}
    for(double bad:{-1.,std::numeric_limits<double>::infinity(),std::numeric_limits<double>::quiet_NaN()}){
        auto c=config;c.liquidViscosity=bad;fail(c,clock);c=config;c.gasViscosity=bad;fail(c,clock);}
    for(int bad:{0,13}){auto c=config;c.maximumIterations=bad;fail(c,clock);}
    for(double bad:{0.,-1.,std::numeric_limits<double>::infinity(),std::numeric_limits<double>::quiet_NaN()})
        for(int target=0;target<3;++target){
            auto c=config;if(target==0)c.velocityUpdateTolerance=bad;else if(target==1)c.geometryUpdateTolerance=bad;else c.momentumPowerTolerance=bad;
            bool rejected=false;
            try{HydroTwoPhaseStep::prepare(mesh,accepted,clock,c);}
            catch(const std::exception& e){check(std::string(e.what()).find("tolerance")!=std::string::npos,
                std::string("invalid update tolerance rejected at a different stage: ")+e.what());rejected=true;++updateToleranceRejects;}
            check(rejected,"nonpositive/nonfinite corrector update tolerance accepted");
            check(accepted.sameIdentity(retained),"update tolerance rejection changed the accepted handle");
            unchanged(saved,accepted.publication(),accepted.values());
        }
    std::cout<<"STEP_UPDATE_TOLERANCE_REJECTS,"<<updateToleranceRejects<<",ATOMIC_PASS\n";
    {auto c=config;c.spaceTime.firstJetTimeScheme=static_cast<HydroFirstJetTimeScheme>(999);fail(c,clock,"invalid FirstJet time scheme");}
    {
        auto c=config.spaceTime;c.firstJetTimeScheme=HydroFirstJetTimeScheme::ClosedProjectionCompensated;
        const auto& reconstruction=accepted.reconstruction();
        // The native data and geometry are unchanged; omitting the declared
        // normal law must reject BEFORE preparing any sampled FirstJet.
        stepRejects([&]{
            const auto previous=reconstruction.evaluate(mesh,accepted.values(),accepted.geometryIdentity(),accepted.timeLayer());
            const auto arrival=reconstruction.evaluate(mesh,accepted.values(),accepted.geometryIdentity(),clock);
            const auto missingLaw=accepted.characteristicGeometry().linear(mesh,previous,arrival);
            HydroMaterialTransport::prepareSpaceTime(mesh,accepted.publication(),missingLaw,config.liquidDensity,config.gasDensity,c);
        },"compensated FirstJet time scheme requires complete uniform XY impermeable base-native sources");
        check(accepted.sameIdentity(retained),"missing-law rejection changed accepted handle");unchanged(saved,accepted.publication(),accepted.values());
        int prescribedQueries=0;
        const HermiteGalsBoundaryFirstJetQuery prescribed=[&](Vector2D,HermiteGalsTimeLayer,Boundary)->HermiteGalsBoundaryFirstJet{
            ++prescribedQueries;throw std::runtime_error("prescribed query executed before time-scheme rejection");};
        stepRejects([&]{
            const auto previous=reconstruction.evaluate(mesh,accepted.values(),accepted.geometryIdentity(),accepted.timeLayer(),accepted.boundary());
            const auto arrival=reconstruction.evaluate(mesh,accepted.values(),accepted.geometryIdentity(),clock,accepted.boundary());
            const auto completeLaw=accepted.characteristicGeometry().linear(mesh,previous,arrival);
            HydroMaterialTransport::prepareSpaceTime(mesh,accepted.publication(),completeLaw,config.liquidDensity,config.gasDensity,c,prescribed);
        },"compensated FirstJet time scheme requires complete uniform XY impermeable base-native sources");
        check(prescribedQueries==0,"unsupported projection reached a sampled boundary query");
        check(accepted.sameIdentity(retained),"prescribed-law rejection changed accepted handle");unchanged(saved,accepted.publication(),accepted.values());
        std::cout<<"HYDRO_TIME_RULE_GUARDS,2,prescribed_queries="<<prescribedQueries<<"\n";
    }
    for(double bad:{-1.,std::numeric_limits<double>::infinity(),std::numeric_limits<double>::quiet_NaN()}){auto c=config;c.surfaceTension=bad;fail(c,clock);}
    for(auto bad:std::vector<HermiteGalsTimeLayer>{{2,0,.001,Position::End},{1,.0001,.001,Position::End},{1,0,0,Position::End},
            {1,0,-.001,Position::End},{1,0,.001,Position::Midpoint},{1,0,.001,Position::Initial},
            {1,0,std::numeric_limits<double>::infinity(),Position::End},{1,0,std::numeric_limits<double>::quiet_NaN(),Position::End}})fail(config,bad);
    auto limited=config;limited.maximumIterations=1;int diagnosticCalls=0;
    limited.onIteration=[&](const HydroTwoPhaseIterationReport& r){++diagnosticCalls;inspectIteration(r,1,clock,static_cast<std::size_t>(limited.spaceTime.initialPanels));
        unchanged(saved,accepted.publication(),accepted.values());
        throw std::runtime_error("STEP_ITERATION_DIAGNOSTIC_THROW");};
    fail(limited,clock,"STEP_ITERATION_DIAGNOSTIC_THROW");
    check(diagnosticCalls==1,"diagnostic exception was swallowed or the completed iteration callback was repeated");
    std::cout<<"STEP_ITERATION_DIAGNOSTIC_ATOMIC_REJECTION,PASS\n";
    const auto fresh=stepMetrics(mesh);stepRejects([&]{HydroTwoPhaseState::create(mesh,geometry,stepVelocity(mesh,fresh.layout()));});
    auto badVelocity=accepted.values().normalVelocity();
    for(int i=0;i<geometry.layout().velocityCount();++i){const auto f=geometry.layout().domain().face(mesh,geometry.layout().nativeFace(i));
        if(f.negativeCell<0||f.positiveCell<0){badVelocity[i]=1;break;}}
    stepRejects([&]{HydroTwoPhaseState::create(mesh,geometry,HydroNativeValues(geometry.layout(),badVelocity,accepted.values().pressure()));});
    stepRejects([&]{HydroTwoPhaseState::create(mesh,geometry,HydroNativeValues(geometry.layout(),accepted.values().normalVelocity(),
        std::vector<double>(geometry.layout().pressureCount(),1)));});
    badVelocity=accepted.values().normalVelocity();
    for(int i=0;i<geometry.layout().velocityCount();++i){const auto f=geometry.layout().domain().face(mesh,geometry.layout().nativeFace(i));
        if(f.negativeCell>=0&&f.positiveCell>=0){badVelocity[i]+=.1;break;}}
    stepRejects([&]{HydroTwoPhaseState::create(mesh,geometry,HydroNativeValues(geometry.layout(),badVelocity,accepted.values().pressure()));});
    auto& mutableMesh=const_cast<QuadtreeMesh&>(mesh);auto& v=mutableMesh.verts(mutableMesh.vertIndexToLabel[0]);const auto old=v.coord();v.setCoord(old+Vector2D{.01,0});
    fail(config,clock);v.setCoord(old);
    Case radial;radial.pureFirstJet=true;radial.radial=true;auto rzMesh=makeMesh(4);const auto rzGeometry=HydroPhaseGeometry::build(*rzMesh,stepMetrics(*rzMesh,true),initial(*rzMesh,radial));
    auto badAxis=std::vector<double>(rzGeometry.layout().velocityCount());
    check(!rzGeometry.metrics().axisVelocityDofs().empty(),"RZ axis constraint fixture has no radial axis");
    badAxis[rzGeometry.metrics().axisVelocityDofs().front()]=1;
    stepRejects([&]{HydroTwoPhaseState::create(*rzMesh,rzGeometry,HydroNativeValues(rzGeometry.layout(),badAxis,std::vector<double>(rzGeometry.layout().pressureCount())));});
    const auto subset=HydroPhaseGeometry::build(mesh,stepMetrics(mesh,false,true),geometry.publication());
    stepRejects([&]{HydroTwoPhaseState::create(mesh,subset,stepVelocity(mesh,subset.layout(),true));});
    auto amrMesh=makeMesh(4,1);Case amr;amr.pureFirstJet=true;const auto amrGeometry=HydroPhaseGeometry::build(*amrMesh,stepMetrics(*amrMesh),initial(*amrMesh,amr));
    stepRejects([&]{HydroTwoPhaseState::create(*amrMesh,amrGeometry,stepVelocity(*amrMesh,amrGeometry.layout(),true));});
    unchanged(saved,accepted.publication(),accepted.values());
}


// Explicit diagnosis only: one production AC step, no artificial power gate.
// Both sharp surface integrals use the SAME retained midpoint quadrature.
struct CsfCouplingSurface
{
    long double area=0,kappaIntegral=0,kappaSquare=0;
    long double characteristicPower=0,solvedPower=0;
    long double characteristicFlux=0,solvedFlux=0,characteristicAbsFlux=0,solvedAbsFlux=0;
    long double characteristicDivSquare=0,solvedDivSquare=0,velocityDifferenceSquare=0;
    double kappaMin=std::numeric_limits<double>::infinity(),kappaMax=-std::numeric_limits<double>::infinity();
    double maxCharacteristicDiv=0,maxSolvedDiv=0,maxInterpolationGap=0;
    std::size_t nodes=0;
};
const char* csfClockPosition(Position p)
{return p==Position::Initial?"initial":p==Position::Midpoint?"midpoint":"end";}
void csfSourceRow(std::ostream& out,const char* name,HermiteGalsIdentity geometry,
    HydroAdvectorIdentity velocity,HermiteGalsTimeLayer clock)
{
    out<<name<<','<<geometry.value()<<','<<velocity.value()<<','<<clock.step<<','
        <<clock.begin<<','<<clock.end<<','<<csfClockPosition(clock.position)<<','<<clock.time()<<'\n';
}
double csfNativeD(const QuadtreeMesh& mesh,const HydroPhaseGeometry& phase,const HydroNativeValues& values)
{
    double result=0;for(double d:phase.metrics().apply(mesh,values))result=std::max(result,std::abs(d));return result;
}
CsfCouplingSurface csfSurfaceObservation(const QuadtreeMesh& mesh,const HydroPhaseGeometry& phase,
    const HydroCharacteristicSource& motion,const HydroSpatialAdvector& actualMidpoint,
    const HydroSpatialAdvector& solvedMidpoint,double sigma,std::ostream* rows,int graphSubdivisions=0)
{
    const auto& publication=phase.publication();const auto& layout=phase.layout();
    actualMidpoint.validateFor(mesh,actualMidpoint.momentumGeometry(),phase.geometryIdentity(),phase.timeLayer());
    solvedMidpoint.validateFor(mesh,actualMidpoint.momentumGeometry(),phase.geometryIdentity(),phase.timeLayer());
    check(actualMidpoint.reconstruction().sameIdentity(motion.geometry().reconstruction())
        &&solvedMidpoint.reconstruction().sameIdentity(actualMidpoint.reconstruction()),
        "CSF diagnosis compared different reconstruction owners");
    if(rows)*rows<<"component,graph,node,native_cell,pressure_dof,x,y,nx,ny,kappa,physical_weight,"
        "characteristic_ux,characteristic_uy,characteristic_Jxx,characteristic_Jxy,characteristic_Jyx,characteristic_Jyy,"
        "solved_midpoint_ux,solved_midpoint_uy,solved_midpoint_Jxx,solved_midpoint_Jxy,solved_midpoint_Jyx,solved_midpoint_Jyy\n";
    CsfCouplingSurface out;
    for(std::size_t c=0;c<publication.topology().curves.size();++c)
        for(std::size_t g=0;g<publication.topology().curves[c].graphs().size();++g){
            std::size_t node=0;
            for(int part=0;part<std::max(1,graphSubdivisions);++part){
            const auto quadrature=graphSubdivisions>0
                ?publication.quadrature(c,g,double(part)/graphSubdivisions,double(part+1)/graphSubdivisions)
                :publication.quadrature(c,g);
            for(const auto& q:quadrature){
                const auto& f=q.sample.frame;const auto point=f.point;const int p=layout.pressureDof(f.nativeCellIndex);
                check(p>=0&&q.sample.geometryIdentity==phase.geometryIdentity()&&sameClock(q.sample.timeLayer,phase.timeLayer())
                    &&q.physicalWeight>0&&std::isfinite(q.physicalWeight),"CSF diagnosis lost midpoint surface owner");
                int piece=-1;
                for(int candidate:motion.geometry().locatePieces(point))
                    if(motion.geometry().reconstruction().pieces()[candidate].pressure==p){piece=candidate;break;}
                check(piece>=0,"CSF diagnosis could not locate actual native characteristic piece");
                const auto characteristic=motion.sample(piece,point,phase.timeLayer().time()).spatial;
                const auto actual=actualMidpoint.cellVelocityJet(p,point).jet;
                const auto solved=solvedMidpoint.cellVelocityJet(p,point).jet;
                const auto gap=std::hypot(characteristic.value.x-actual.value.x,characteristic.value.y-actual.value.y);
                const auto scale=std::max({1.,std::hypot(characteristic.value.x,characteristic.value.y),
                    std::hypot(actual.value.x,actual.value.y)});
                // This checks source identity/interpolation, not equality of force and sharp-area powers.
                check(gap<=1e-12*scale,"CSF actual characteristic and retained midpoint field disagree");
                out.maxInterpolationGap=std::max(out.maxInterpolationGap,gap);
                const long double w=q.physicalWeight,k=q.sample.curvature;
                const long double nc=f.normal.x*characteristic.value.x+f.normal.y*characteristic.value.y;
                const long double ns=f.normal.x*solved.value.x+f.normal.y*solved.value.y;
                const long double dc=characteristic.derivative.xx+characteristic.derivative.yy;
                const long double ds=solved.derivative.xx+solved.derivative.yy;
                const long double dux=solved.value.x-characteristic.value.x,duy=solved.value.y-characteristic.value.y;
                check(std::isfinite(k)&&std::isfinite(nc)&&std::isfinite(ns)&&std::isfinite(dc)&&std::isfinite(ds),
                    "CSF diagnosis encountered nonfinite actual surface data");
                out.area+=w;out.kappaIntegral+=w*k;out.kappaSquare+=w*k*k;
                out.kappaMin=std::min(out.kappaMin,double(k));out.kappaMax=std::max(out.kappaMax,double(k));
                out.characteristicPower-=sigma*w*k*nc;out.solvedPower-=sigma*w*k*ns;
                out.characteristicFlux+=w*nc;out.solvedFlux+=w*ns;
                out.characteristicAbsFlux+=w*std::abs(nc);out.solvedAbsFlux+=w*std::abs(ns);
                out.characteristicDivSquare+=w*dc*dc;out.solvedDivSquare+=w*ds*ds;
                out.velocityDifferenceSquare+=w*(dux*dux+duy*duy);
                out.maxCharacteristicDiv=std::max(out.maxCharacteristicDiv,double(std::abs(dc)));
                out.maxSolvedDiv=std::max(out.maxSolvedDiv,double(std::abs(ds)));
                ++out.nodes;
                if(rows)*rows<<c<<','<<g<<','<<node<<','<<f.nativeCellIndex<<','<<p<<','<<point.x<<','<<point.y<<','
                    <<f.normal.x<<','<<f.normal.y<<','<<k<<','<<w<<','<<characteristic.value.x<<','<<characteristic.value.y<<','
                    <<characteristic.derivative.xx<<','<<characteristic.derivative.xy<<','
                    <<characteristic.derivative.yx<<','<<characteristic.derivative.yy<<','
                    <<solved.value.x<<','<<solved.value.y<<','<<solved.derivative.xx<<','<<solved.derivative.xy<<','
                    <<solved.derivative.yx<<','<<solved.derivative.yy<<'\n';
                ++node;
            }
            }
        }
    check(out.nodes>0&&out.area>0,"CSF diagnosis has no actual midpoint surface quadrature");return out;
}
HermiteGalsPublication csfInitialEllipse(const QuadtreeMesh& mesh,double a,double b,
    HermiteGalsGradientUpdate update=HermiteGalsGradientUpdate::CharacteristicPullback)
{
    const Profile profile=[a,b](Vector2D p){const double x=p.x-.5,y=p.y-.5;
        return Jet{x*x/(a*a)+y*y/(b*b)-1,2*x/(a*a),2*y/(b*b),0};};
    const auto seed=stateFor(mesh,profile,false,{},[](int,Vector2D,Vector2D){return 0.;},{},HermiteGalsTransport::FirstJet);
    auto primary=seed.primary();primary.gradientUpdate=update;
    HermiteGalsCandidate candidate(mesh,{},primary);const auto mixed=closeHermiteGalsMixedAmr(mesh,candidate);
    const auto identity=candidate.identity();const auto targetCenters=centers(mesh,profile);
    HermiteGalsPublicationConfig config;config.gradientUpdate=update;
    return publishHermiteGalsFirstJetGeometry(mesh,std::move(candidate),mixed,{identity,targetCenters},config);
}
void csfLoadObservation(const char* shape,const HydroBalancedCsfLoad& load)
{
    const auto& fractions=load.phaseGeometry().liquidFractions();
    const auto range=std::minmax_element(fractions.begin(),fractions.end());
    double lower=std::numeric_limits<double>::infinity(),upper=-lower;std::size_t defined=0;
    if(load.discretization()==HydroCsfDiscretization::GeometricAdjoint
        ||load.discretization()==HydroCsfDiscretization::GeometricProjectedAdjoint){
        for(double k:load.referenceCurvatures())if(std::isfinite(k)){lower=std::min(lower,k);upper=std::max(upper,k);++defined;}
    }else for(std::size_t u=0;u<load.curvatureAverages().size();++u)if(load.supportReports()[u].curvatureDefined){
        const double k=load.curvatureAverages()[u];lower=std::min(lower,k);upper=std::max(upper,k);++defined;}
    check(defined>0&&*range.first>=0&&*range.second<=1,"CSF diagnosis has no local curvature or bounded phase fraction");
    std::cout<<"CSF_PHASE_CURVATURE,"<<shape<<",phase_fraction_min,"<<*range.first<<",phase_fraction_max,"<<*range.second
        <<",local_kappa_min,"<<lower<<",local_kappa_max,"<<upper<<",loaded_faces,"<<defined<<'\n';
}

// Fixed actual-source observation only: no new accepted geometry, trajectory
// truth, Hessian state, or FD curvature accuracy/stability acceptance gate.
void firstJetReplayOnly(const std::string& directory)
{
    namespace fs=std::filesystem;const fs::path output(directory);
    check(output.is_absolute()&&(output.root_name()==fs::path("D:")||output.root_name()==fs::path("d:")),"FirstJet replay output must be absolute on D:");
    fs::create_directories(output);
    for(const char* name:{"primary_replay.csv","seam_hessian.csv","source_clock.csv"})check(!fs::exists(output/name),"FirstJet replay refuses to overwrite observation");
#ifdef _OPENMP
    omp_set_dynamic(0);omp_set_num_threads(8);
#endif
    constexpr double a=.168,b=.15238095238095239,dt=.0002;
    auto mesh=makeMesh(64);const auto phase=HydroPhaseGeometry::build(*mesh,stepMetrics(*mesh),csfInitialEllipse(*mesh,a,b));
    auto accepted=HydroTwoPhaseState::create(*mesh,phase,stepVelocity(*mesh,phase.layout(),true));
    HydroTwoPhaseStepConfig cfg;cfg.liquidDensity=100;cfg.gasDensity=1;cfg.liquidViscosity=cfg.gasViscosity=0;cfg.surfaceTension=341.642;
    cfg.momentumDensityRepresentation=HydroMomentumMass::DensityRepresentation::CellAverageProjection;
    cfg.csfCurvature.forceMethod=HydroCsfForceMethod::HeightFunction;cfg.spaceTime.characteristic.maxStep=dt;cfg.spaceTime.initialPanels=1;
    cfg.onIteration=[](const HydroTwoPhaseIterationReport& r){std::cout<<"FIRST_JET_REPLAY_STAGE,"<<r.iteration<<",transport,"<<r.stageSeconds[1]<<",capillary,"<<r.csfProjectionSeconds<<'\n'<<std::flush;};
    std::cout<<"FIRST_JET_REPLAY_BEGIN,64x64,rho100:1,mu0,sigma341.642,dt.0002,P0,C3,AC,HF,one_actual_step\n"<<std::flush;
    const auto saved=capture(accepted.publication(),accepted.values());const auto step=HydroTwoPhaseStep::prepare(*mesh,accepted,{1,0,dt,Position::End},cfg);
    unchanged(saved,accepted.publication(),accepted.values());const auto& motion=step.characteristicSource();const auto& end=step.materialTransport().endGeometry().publication();
    const auto& scalar=end.state().geometry().snapshot();
    check(motion.sameIdentity(step.materialTransport().characteristicSource())&&end.state().identity()==step.candidateState().geometryIdentity(),"replay detached actual motion/End owner");
    check(step.materialTransport().spaceTimeReport().firstJetTimeProjection==HermiteGalsTimeProjection::CharacteristicPullback
        &&scalar.reconstruction()==LevelSetHermiteGeometry::Reconstruction::EnrichedQuartic,"replay requires actual single pullback/SingleQ4 path");
    std::ofstream provenance(output/"source_clock.csv"),primary(output/"primary_replay.csv"),hessian(output/"seam_hessian.csv");
    check(provenance.good()&&primary.good()&&hessian.good(),"replay output open failed");provenance<<std::setprecision(17);primary<<std::setprecision(17);hessian<<std::setprecision(17);
    provenance<<"field,geometry_identity,advector_identity,step,begin,end,position,physical_time\n";
    const auto sourceRow=[&](const char* name,HermiteGalsIdentity id,HydroAdvectorIdentity velocity,HermiteGalsTimeLayer time){
        std::cout<<"FIRST_JET_REPLAY_SOURCE,";csfSourceRow(std::cout,name,id,velocity,time);csfSourceRow(provenance,name,id,velocity,time);};
    sourceRow("begin_scalar",accepted.geometryIdentity(),{},accepted.timeLayer());sourceRow("candidate_End_scalar",end.state().identity(),{},end.state().timeLayer());
    sourceRow("actual_previous",motion.previous().geometryIdentity(),motion.previous().identity(),motion.previous().timeLayer());
    sourceRow("actual_arrival",motion.arrival().geometryIdentity(),motion.arrival().identity(),motion.arrival().timeLayer());
    const auto& midpoint=step.materialTransport().midpointRate().spatialAdvector();sourceRow("actual_midpoint",midpoint.geometryIdentity(),midpoint.identity(),midpoint.timeLayer());
    struct Seam {HermiteGalsSurfaceSample left,right;std::size_t component,graph,next;};std::optional<Seam> seam;double jump=-1;
    for(std::size_t c=0;c<end.topology().curves.size();++c){const auto& curve=end.topology().curves[c];const auto count=curve.graphs().size();
        for(std::size_t g=0;g<count;++g){if(g+1==count&&!curve.closed())continue;const auto next=(g+1)%count;
            const auto left=end.sample(c,g,1),right=end.sample(c,next,0);if(left.frame.nativeCellIndex==right.frame.nativeCellIndex)continue;
            const auto l=scalar.nativeCellBounds(left.frame.nativeCellIndex),r=scalar.nativeCellBounds(right.frame.nativeCellIndex);
            const bool adjacent=l[1].x==r[0].x||r[1].x==l[0].x||l[1].y==r[0].y||r[1].y==l[0].y;
            if(!adjacent||std::hypot(left.frame.point.x-right.frame.point.x,left.frame.point.y-right.frame.point.y)>1e-11)continue;
            const double candidate=std::abs(left.curvature-right.curvature);if(candidate>jump){jump=candidate;seam=Seam{left,right,c,g,next};}}
    }
    check(seam.has_value(),"replay lacks an observed genuine native-cell seam");
    std::set<int> vertices;std::set<Key> edges;
    for(int cell:{seam->left.frame.nativeCellIndex,seam->right.frame.nativeCellIndex}){const auto& native=mesh->cells(mesh->cellIndexToLabel[cell]);
        for(int k=0;k<4;++k){vertices.insert(mesh->verts(native.cellToVert[k]).index());Key key{{mesh->verts(native.cellToVert[k]).index(),mesh->verts(native.cellToVert[(k+1)%4]).index()}};
            const auto x=point(*mesh,key[0]),y=point(*mesh,key[1]);if(x.y==y.y?x.x>y.x:x.y>y.y)std::swap(key[0],key[1]);edges.insert(key);}}
    const auto starts=motion.prepareNativeStarts(*mesh,accepted.publication());int queries=0;
    struct Pull {Jet jet;HydroCharacteristicResult trace;};
    const auto pull=[&](Vector2D p){++queries;const auto trace=motion.traceBackward(p,dt,0,cfg.spaceTime.characteristic,starts);
        check(!trace.boundary()&&trace.jacobian().has_value(),"replay interior direct query lacks actual source Jacobian");const auto q=trace.point();const double x=q.x-.5,y=q.y-.5;
        const auto gradient=trace.jacobian()->transposeTimes({2*x/(a*a),2*y/(b*b)});
        return Pull{Jet{x*x/(a*a)+y*y/(b*b)-1,gradient.x,gradient.y,0},trace};};
    primary<<"kind,id,x,y,actual_phi,pulled_phi,value_error,actual_dx,actual_dy,pulled_dx,pulled_dy,gradient_error,departure_x,departure_y,Jxx,Jxy,Jyx,Jyy,initial_piece,final_piece,one_sided\n";
    double maxValue=0,maxGradient=0;int primaryPoints=0;const double unavailable=std::numeric_limits<double>::quiet_NaN();
    const auto primaryRow=[&](const char* kind,int id,Vector2D p,double value,double dx,double dy){const auto direct=pull(p);const auto& d=direct.jet;const auto& t=direct.trace;
        const double ve=std::abs(value-d.value),ge=std::hypot(std::isfinite(dx)?dx-d.dx:0,std::isfinite(dy)?dy-d.dy:0);
        maxValue=std::max(maxValue,ve);maxGradient=std::max(maxGradient,ge);++primaryPoints;const auto& J=*t.jacobian();
        primary<<kind<<','<<id<<','<<p.x<<','<<p.y<<','<<value<<','<<d.value<<','<<ve<<','<<dx<<','<<dy<<','<<d.dx<<','<<d.dy<<','<<ge<<','<<t.point().x<<','<<t.point().y<<','<<J.xx<<','<<J.xy<<','<<J.yx<<','<<J.yy<<','<<t.initialPiece()<<','<<t.finalPiece()<<','<<t.oneSided()<<'\n';};
    for(const auto& v:end.state().primary().vertices)if(vertices.count(v.vertex))primaryRow("vertex",v.vertex,point(*mesh,v.vertex),v.phi,v.dx,v.dy);
    for(std::size_t id=0;id<end.state().primary().edges.size();++id){const auto& e=end.state().primary().edges[id];if(!edges.count(e.vertices))continue;
        const auto x=point(*mesh,e.vertices[0]),y=point(*mesh,e.vertices[1]);const bool horizontal=x.y==y.y;
        primaryRow("edge",static_cast<int>(id),x+(y-x)*.5,e.phi,horizontal?unavailable:e.dn,horizontal?e.dn:unavailable);}
    hessian<<"point,owner,graph,x,y,delta,stencil_x,stencil_y,scalar_dxx,scalar_dxy,scalar_dyy,scalar_kappa,direct_fd_dxx,direct_fd_dxy,direct_fd_dyx,direct_fd_dyy,direct_fd_kappa,fd_hessian_change,fd_kappa_change,direct_phi_at_scalar_curve,direct_geometry_residual,scalar_phi,initial_piece,final_piece,one_sided\n";
    const auto curvature=[](Vector2D g,const std::array<double,4>& H){const double n=std::hypot(g.x,g.y);return (g.y*g.y*H[0]-g.x*g.y*(H[1]+H[2])+g.x*g.x*H[3])/(n*n*n);};
    const std::array<HermiteGalsSurfaceSample,4> points{{seam->left,seam->right,end.sample(seam->component,seam->graph,.99),end.sample(seam->component,seam->next,.01)}};
    for(int index=0;index<4;++index){const auto& q=points[index];const auto p=q.frame.point;const int cell=q.frame.nativeCellIndex;const auto box=scalar.nativeCellBounds(cell);
        const double width=std::min(box[1].x-box[0].x,box[1].y-box[0].y);const auto s=scalar.scalarJetInNativeCell(cell,p);const auto direct=pull(p);const Vector2D g{direct.jet.dx,direct.jet.dy};
        const double sk=curvature({s.dx,s.dy},{{s.dxx,s.dxy,s.dxy,s.dyy}});std::optional<std::array<double,4>> coarse;double coarseK=0;
        for(double scale:{1e-3,2.5e-4}){const double delta=scale*width;std::array<int,2> stencil{};std::array<Vector2D,2> derivative;
            for(int axis=0;axis<2;++axis){const double coordinate=axis==0?p.x:p.y,lo=axis==0?box[0].x:box[0].y,hi=axis==0?box[1].x:box[1].y;
                const int sign=coordinate-lo<2*delta?1:hi-coordinate<2*delta?-1:0;stencil[axis]=sign;const Vector2D shift{axis==0?delta:0,axis==1?delta:0};
                if(sign){const auto one=pull(p+shift*sign).jet,two=pull(p+shift*(2*sign)).jet;
                    derivative[axis]={sign*(-3*g.x+4*one.dx-two.dx)/(2*delta),sign*(-3*g.y+4*one.dy-two.dy)/(2*delta)};}
                else{const auto plus=pull(p+shift).jet,minus=pull(p-shift).jet;derivative[axis]={(plus.dx-minus.dx)/(2*delta),(plus.dy-minus.dy)/(2*delta)};}}
            const std::array<double,4> H{{derivative[0].x,derivative[1].x,derivative[0].y,derivative[1].y}};const double dk=curvature(g,H);double change=unavailable;
            if(coarse){change=0;for(int k=0;k<4;++k)change=std::max(change,std::abs(H[k]-(*coarse)[k]));}
            hessian<<index<<','<<cell<<','<<(index%2?seam->next:seam->graph)<<','<<p.x<<','<<p.y<<','<<delta<<','<<stencil[0]<<','<<stencil[1]<<','<<s.dxx<<','<<s.dxy<<','<<s.dyy<<','<<sk<<','<<H[0]<<','<<H[1]<<','<<H[2]<<','<<H[3]<<','<<dk<<','<<change<<','<<(coarse?std::abs(dk-coarseK):unavailable)<<','<<direct.jet.value<<','<<direct.jet.value/std::hypot(g.x,g.y)<<','<<s.value<<','<<direct.trace.initialPiece()<<','<<direct.trace.finalPiece()<<','<<direct.trace.oneSided()<<'\n';
            coarse=H;coarseK=dk;}}

    // Non-accepted spatial prototype. Delta-phi and all center predictions
    // remain the actual first step; only its gradient increments are projected.
    auto projectedPrimary=end.state().primary();const auto& oldPrimary=accepted.publication().state().primary();
    using PointKey=std::pair<double,double>;std::map<PointKey,double> deltaPhi;
    std::map<int,const HermiteGalsVertexJet*> oldVertices;
    for(const auto& v:oldPrimary.vertices)oldVertices.emplace(v.vertex,&v);
    std::map<Key,const HermiteGalsEdgeJet*> oldEdges;
    for(const auto& e:oldPrimary.edges)oldEdges.emplace(e.vertices,&e);
    for(const auto& v:projectedPrimary.vertices){const auto p=point(*mesh,v.vertex);deltaPhi[{p.x,p.y}]=v.phi-oldVertices.at(v.vertex)->phi;}
    for(const auto& e:projectedPrimary.edges){const auto p=(point(*mesh,e.vertices[0])+point(*mesh,e.vertices[1]))*.5;deltaPhi[{p.x,p.y}]=e.phi-oldEdges.at(e.vertices)->phi;}
    std::vector<double> actualCenters;const auto& beforeScalar=accepted.publication().state().geometry().snapshot();
    for(int cell:projectedPrimary.geometryCells){const auto box=scalar.nativeCellBounds(cell);const auto p=(box[0]+box[1])*.5;
        const double actual=scalar.scalarJetInNativeCell(cell,p).value;actualCenters.push_back(actual);
        deltaPhi[{p.x,p.y}]=actual-beforeScalar.scalarJetInNativeCell(cell,p).value;}
    const auto incrementDerivative=[&](Vector2D p,int axis)->std::optional<double>{
        constexpr double h=1./64;const Vector2D half{axis==0?h/2:0,axis==1?h/2:0};std::array<double,4> f;
        const int offsets[4]={-2,-1,1,2};for(int k=0;k<4;++k){const auto q=p+half*offsets[k];const auto found=deltaPhi.find({q.x,q.y});
            if(found==deltaPhi.end())return {};f[k]=found->second;}
        const long double derivative=(static_cast<long double>(f[0])-8.L*f[1]+8.L*f[2]-f[3])/(6.L*h);
        check(std::isfinite(derivative),"gradient-increment prototype has nonfinite derivative");return static_cast<double>(derivative);};
    std::size_t projectedComponents=0;double largestProjection=0;
    for(auto& v:projectedPrimary.vertices){const auto p=point(*mesh,v.vertex);const auto* old=oldVertices.at(v.vertex);
        if(const auto d=incrementDerivative(p,0)){const double value=old->dx+*d;largestProjection=std::max(largestProjection,std::abs(value-v.dx));v.dx=value;++projectedComponents;}
        if(const auto d=incrementDerivative(p,1)){const double value=old->dy+*d;largestProjection=std::max(largestProjection,std::abs(value-v.dy));v.dy=value;++projectedComponents;}}
    for(auto& e:projectedPrimary.edges){const auto x=point(*mesh,e.vertices[0]),y=point(*mesh,e.vertices[1]);const auto p=(x+y)*.5;
        if(const auto d=incrementDerivative(p,x.y==y.y?1:0)){const double value=oldEdges.at(e.vertices)->dn+*d;largestProjection=std::max(largestProjection,std::abs(value-e.dn));e.dn=value;++projectedComponents;}}
    HermiteGalsCandidate referenceCandidate(*mesh,end.state().timeLayer(),projectedPrimary);
    const auto referenceMixed=closeHermiteGalsMixedAmr(*mesh,referenceCandidate);
    const auto referenceId=referenceCandidate.identity();const auto reference=publishHermiteGalsFirstJetGeometry(*mesh,
        std::move(referenceCandidate),referenceMixed,{referenceId,actualCenters},{});
    double referenceJump=0;for(std::size_t c=0;c<reference.topology().curves.size();++c){const auto& curve=reference.topology().curves[c];const auto count=curve.graphs().size();
        for(std::size_t g=0;g<count;++g){if(g+1==count&&!curve.closed())continue;const auto next=(g+1)%count;const auto l=reference.sample(c,g,1),r=reference.sample(c,next,0);
            if(l.frame.nativeCellIndex==r.frame.nativeCellIndex)continue;referenceJump=std::max(referenceJump,std::abs(l.curvature-r.curvature));}}
    long double originalVolume=0,referenceVolume=0;for(const auto& v:end.volumes().phaseVolumes())originalVolume+=v.liquid;
    for(const auto& v:reference.volumes().phaseVolumes())referenceVolume+=v.liquid;
    std::ofstream projection(output/"increment_projection.csv");check(projection.good(),"projection observation open failed");projection<<std::setprecision(17);
    projection<<"scope,actual_End_identity,unaccepted_reference_identity,projected_components,max_gradient_change,original_max_seam_jump,reference_max_seam_jump,relative_volume_change\n"
        <<"unaccepted_spatial_prototype,"<<end.state().identity().value()<<','<<reference.state().identity().value()<<','<<projectedComponents<<','<<largestProjection<<','<<jump<<','<<referenceJump<<','<<(referenceVolume-originalVolume)/originalVolume<<'\n';
    projection.flush();check(projection.good(),"projection observation write failed");unchanged(saved,accepted.publication(),accepted.values());
    std::cout<<"GRADIENT_INCREMENT_PROTOTYPE,UNACCEPTED,components,"<<projectedComponents<<",maximum_gradient_change,"<<largestProjection<<",raw_seam_jump,"<<jump<<",reference_seam_jump,"<<referenceJump<<",relative_volume_change,"<<(referenceVolume-originalVolume)/originalVolume<<'\n';

    primary.flush();hessian.flush();provenance.flush();check(primary.good()&&hessian.good()&&provenance.good(),"replay diagnostic write failed");
    std::cout<<"FIRST_JET_REPLAY_COMPLETED,primary_points,"<<primaryPoints<<",primary_value_error_max,"<<maxValue<<",primary_gradient_error_max,"<<maxGradient<<",native_seam_kappa_jump,"<<jump<<",left_owner,"<<seam->left.frame.nativeCellIndex<<",right_owner,"<<seam->right.frame.nativeCellIndex<<",direct_queries,"<<queries<<",FD_observation_only,no_new_accuracy_or_stability_PASS\n";
}

void csfCircleStatic(const QuadtreeMesh& mesh,const HydroPhysicalDivergence& metrics,const HydroTwoPhaseStepConfig& config)
{
    constexpr double radius=.16;
    const auto phase=HydroPhaseGeometry::build(mesh,metrics,csfInitialEllipse(mesh,radius,radius));
    const auto pressure=HydroPressureLoad::build(metrics.divergence());
    const auto load=HydroBalancedCsfLoad::geometricCapillary(mesh,pressure,phase,config.surfaceTension,config.csfCurvature);
    csfLoadObservation("circle_static",load);
    long double forceError=0,forceScale=0;double kappaError=0;
    for(int u=0;u<metrics.layout().velocityCount();++u){
        const long double expected=config.surfaceTension/radius*load.colorGradient()[u];
        forceError=std::hypot(forceError,static_cast<long double>(load.momentumRhsLoad()[u])-expected);
        forceScale=std::hypot(forceScale,expected);
        if(load.supportReports()[u].curvatureDefined)kappaError=std::max(kappaError,std::abs(load.curvatureAverages()[u]-1/radius));
    }
    const auto area=actualInterfaceMeasure(phase);
    std::cout<<"CSF_CIRCLE_STATIC,reference_kappa,"<<1/radius<<",local_kappa_max_error,"<<kappaError
        <<",force_vs_constant_kappa_Calpha_L2,"<<forceError<<",constant_force_L2,"<<forceScale
        <<",relative_force_error,"<<(forceScale>0?forceError/forceScale:0)
        <<",sharp_area,"<<area<<",reference_area,"<<2*pi*radius
        <<",observation_only,no_time_window\n";
}
void csfCouplingOnly(const std::string& outputDirectory)
{
    namespace fs=std::filesystem;
    fs::path output;
    if(!outputDirectory.empty()){
        output=fs::path(outputDirectory);
        check(output.is_absolute()&&(output.root_name()==fs::path("D:")||output.root_name()==fs::path("d:")),
            "CSF diagnostic output directory must be absolute on D:");
        fs::create_directories(output);
    }
    int threads=1;
#ifdef _OPENMP
    omp_set_dynamic(0);omp_set_num_threads(8);threads=omp_get_max_threads();
#endif
    std::cout<<"CSF_COUPLING_CONFIG,nx,64,ny,64,a,0.168,b,0.15238095238095239,rho_liquid,100,rho_gas,1,"
        "mu_liquid,0,mu_gas,0,sigma,341.642,dt,0.0002,threads,"<<threads
        <<",scheme,adaptiveLinear,steps,1,diagnostic_only,begin\n"<<std::flush;
    auto mesh=makeMesh(64);const auto metrics=stepMetrics(*mesh);
    const auto phase=HydroPhaseGeometry::build(*mesh,metrics,csfInitialEllipse(*mesh,.168,.15238095238095239));
    auto accepted=HydroTwoPhaseState::create(*mesh,phase,stepVelocity(*mesh,phase.layout(),true));
    HydroTwoPhaseStepConfig config;config.liquidDensity=100;config.gasDensity=1;config.surfaceTension=341.642;
    config.liquidViscosity=config.gasViscosity=0;config.spaceTime.characteristic.maxStep=.0002;
    config.spaceTime.initialPanels=1;
    config.onIteration=[](const HydroTwoPhaseIterationReport& r){
        std::cout<<"CSF_COUPLING_STAGE,"<<r.iteration<<",step,"<<r.timeLayer.step
            <<",transport_seconds,"<<r.stageSeconds[1]<<",color_seconds,"<<r.csfColorSeconds
            <<",projection_seconds,"<<r.csfProjectionSeconds<<",saddle_seconds,"<<r.stageSeconds[6]
            <<",momentum_power_relative_defect,"<<r.momentumPowerRelativeDefect<<'\n'<<std::flush;
    };
    const auto saved=capture(accepted.publication(),accepted.values());
    const auto step=HydroTwoPhaseStep::prepare(*mesh,accepted,{1,0,.0002,Position::End},config);
    unchanged(saved,accepted.publication(),accepted.values());
    const auto& midpoint=step.materialTransport().midpointRate().geometry();
    const auto& actual=step.materialTransport().midpointRate().spatialAdvector();
    const auto& motion=step.characteristicSource();const auto& load=step.capillaryLoad();
    check(motion.sameIdentity(step.materialTransport().characteristicSource())
        &&motion.arrival().values().normalVelocity()==step.endpointGuess().normalVelocity()
        &&motion.previous().values().normalVelocity()==step.beginState().values().normalVelocity(),
        "CSF diagnosis detached actual transport endpoint guess from characteristic source");
    inspectCsfOwner(*mesh,step);
    const auto solved=step.beginState().reconstruction().evaluate(*mesh,step.solvedMidpointValues(),
        midpoint.geometryIdentity(),midpoint.timeLayer(),step.beginState().boundary());
    std::ofstream surface,native,provenance;
    if(!output.empty()){
        surface.open(output/"surface_midpoint_step_000001.csv");
        native.open(output/"native_faces_step_000001.csv");
        provenance.open(output/"source_provenance_step_000001.csv");
        check(surface.good()&&native.good()&&provenance.good(),"CSF diagnosis could not open actual source outputs");
        surface<<std::setprecision(17);native<<std::setprecision(17);provenance<<std::setprecision(17);
    }
    std::cout<<"CSF_SOURCE_HEADER,field,geometry_identity,advector_identity,step,begin,end,position,physical_time\n";
    if(provenance.is_open())provenance<<"field,geometry_identity,advector_identity,step,begin,end,position,physical_time\n";
    const auto sourceRow=[&](const char* name,HermiteGalsIdentity geometry,HydroAdvectorIdentity velocity,HermiteGalsTimeLayer clock){
        std::cout<<"CSF_SOURCE,";csfSourceRow(std::cout,name,geometry,velocity,clock);
        if(provenance.is_open())csfSourceRow(provenance,name,geometry,velocity,clock);
    };
    sourceRow("characteristic_previous",motion.previous().geometryIdentity(),motion.previous().identity(),motion.previous().timeLayer());
    sourceRow("characteristic_arrival_endpointGuess",motion.arrival().geometryIdentity(),motion.arrival().identity(),motion.arrival().timeLayer());
    sourceRow("actual_midpoint_transport_interpolation",actual.geometryIdentity(),actual.identity(),actual.timeLayer());
    sourceRow("solved_midpoint_reconstructed_same_source",solved.geometryIdentity(),solved.identity(),solved.timeLayer());
    sourceRow("CSF_load_and_surface_midpoint",load.geometryIdentity(),{},load.timeLayer());
    sourceRow("candidate_solved_endpoint",step.candidateState().geometryIdentity(),{},step.candidateState().timeLayer());
    const auto stats=csfSurfaceObservation(*mesh,midpoint,motion,actual,solved,config.surfaceTension,surface.is_open()?&surface:nullptr);
    long double solvedPower=0,actualPower=0,nativeGap=0,guessNorm=0,solvedNorm=0;
    const auto& l=midpoint.layout();const auto& before=step.beginState().values().normalVelocity();
    const auto& guess=step.endpointGuess().normalVelocity();const auto& bar=step.solvedMidpointValues().normalVelocity();
    const auto& end=step.candidateState().values().normalVelocity();const auto& um=actual.values().normalVelocity();
    const auto& force=load.momentumRhsLoad();
    if(native.is_open())native<<"velocity_dof,native_face,normal_axis,negative_cell,positive_cell,ax,ay,bx,by,cx,cy,physical_measure,"
        "begin_u,endpoint_guess_u,actual_midpoint_u,solved_midpoint_u,candidate_solved_endpoint_u,CSF_force,reference_component_kappa,"
        "curvature_defined\n";
    for(int u=0;u<l.velocityCount();++u){
        solvedPower+=static_cast<long double>(bar[u])*force[u];actualPower+=static_cast<long double>(um[u])*force[u];
        nativeGap=std::hypot(nativeGap,static_cast<long double>(bar[u])-um[u]);
        guessNorm=std::hypot(guessNorm,static_cast<long double>(um[u]));solvedNorm=std::hypot(solvedNorm,static_cast<long double>(bar[u]));
        if(native.is_open()){
            const auto face=l.domain().face(*mesh,l.nativeFace(u));
            const auto& edge=face.horizontal?mesh->hEdges(face.edgeLabel):mesh->vEdges(face.edgeLabel);
            const auto a=mesh->verts(edge.edgeToVert[0]).coord(),b=mesh->verts(edge.edgeToVert[1]).coord(),center=(a+b)*.5;
            native<<u<<','<<l.nativeFace(u)<<','<<(face.horizontal?1:0)<<','<<face.negativeCell<<','<<face.positiveCell
                <<','<<a.x<<','<<a.y<<','<<b.x<<','<<b.y<<','<<center.x<<','<<center.y<<','<<stepFaceMeasure(*mesh,l,u)
                <<','<<before[u]<<','<<guess[u]<<','<<um[u]<<','<<bar[u]<<','<<end[u]<<','<<force[u]<<',';
            const int p=face.negativeCell>=0?l.pressureDof(face.negativeCell):l.pressureDof(face.positiveCell);
            if(std::isfinite(load.referenceCurvatures()[p]))native<<load.referenceCurvatures()[p];else native<<"NA";
            native<<','<<load.supportReports()[u].curvatureDefined<<'\n';
        }
    }
    scaled(solvedPower,step.report().capillaryPower,std::abs(solvedPower)+std::abs(step.report().capillaryPower),
        "CSF diagnosis actual solved midpoint native force work",1e-12);
    check(step.report().frozenOwnerMomentumDefect<=1e-11&&step.report().linearMomentumDefect<=1e-11
        &&step.report().momentumPowerRelativeDefect<=config.momentumPowerTolerance,
        "CSF diagnosis production momentum gate failed");
    nativeConstraints(*mesh,step.candidateState().values());
    csfLoadObservation("ellipse_midpoint",load);
    const auto mean=stats.kappaIntegral/stats.area;
    std::cout<<"CSF_COUPLING_POWER,step,1,stages,"<<step.report().iterations<<",solved_native_uS,"<<solvedPower
        <<",actual_transport_native_uS,"<<actualPower<<",characteristic_sharp_minus_sigma_kappa_nu,"<<stats.characteristicPower
        <<",solved_sharp_minus_sigma_kappa_nu,"<<stats.solvedPower
        <<",same_solved_velocity_CSF_minus_sharp,"<<solvedPower-stats.solvedPower
        <<",same_actual_velocity_CSF_minus_sharp,"<<actualPower-stats.characteristicPower
        <<",sharp_solved_minus_characteristic,"<<stats.solvedPower-stats.characteristicPower
        <<",native_solved_minus_actual_work,"<<solvedPower-actualPower
        <<",area_change_rate,"<<(step.report().interfaceMeasureEnd-step.report().interfaceMeasureBegin)/step.report().deltaTime
        <<",surface_energy_change_rate,"<<step.report().surfaceEnergyChange/step.report().deltaTime
        <<",observation_only,no_power_equality_gate\n";
    std::cout<<"CSF_COUPLING_FIELD,step,1,surface_area,"<<stats.area<<",surface_nodes,"<<stats.nodes
        <<",raw_kappa_min,"<<stats.kappaMin<<",raw_kappa_max,"<<stats.kappaMax<<",raw_kappa_mean,"<<mean
        <<",raw_kappa_stddev,"<<std::sqrt(std::max(0.L,stats.kappaSquare/stats.area-mean*mean))
        <<",characteristic_closed_surface_flux,"<<stats.characteristicFlux
        <<",solved_closed_surface_flux,"<<stats.solvedFlux
        <<",characteristic_abs_normal_flux,"<<stats.characteristicAbsFlux<<",solved_abs_normal_flux,"<<stats.solvedAbsFlux
        <<",actual_native_Dmax,"<<csfNativeD(*mesh,midpoint,actual.values())
        <<",solved_native_Dmax,"<<csfNativeD(*mesh,midpoint,step.solvedMidpointValues())
        <<",actual_point_div_rms,"<<std::sqrt(stats.characteristicDivSquare/stats.area)
        <<",solved_point_div_rms,"<<std::sqrt(stats.solvedDivSquare/stats.area)
        <<",actual_point_div_max,"<<stats.maxCharacteristicDiv<<",solved_point_div_max,"<<stats.maxSolvedDiv
        <<",surface_velocity_guess_solve_rms,"<<std::sqrt(stats.velocityDifferenceSquare/stats.area)
        <<",native_velocity_guess_solve_L2,"<<nativeGap
        <<",native_actual_L2,"<<guessNorm<<",native_solved_L2,"<<solvedNorm
        <<",actual_interpolation_max_gap,"<<stats.maxInterpolationGap<<",capillary_courant,"<<step.report().capillaryCourant<<'\n';
    // Refine only this observational integral of the retained publication/source.
    // The production publication rule, step, load and all original gates stay fixed.
    CsfCouplingSurface preceding=stats;
    for(int split:{1,2,4}){
        const auto refined=csfSurfaceObservation(*mesh,midpoint,motion,actual,solved,config.surfaceTension,nullptr,split);
        std::cout<<"CSF_COUPLING_QUADRATURE,graph_subdivisions,"<<split<<",surface_nodes,"<<refined.nodes
            <<",area,"<<refined.area<<",characteristic_sharp_power,"<<refined.characteristicPower
            <<",solved_sharp_power,"<<refined.solvedPower
            <<",CSF_minus_solved_sharp,"<<solvedPower-refined.solvedPower
            <<",characteristic_closed_surface_flux,"<<refined.characteristicFlux
            <<",solved_closed_surface_flux,"<<refined.solvedFlux
            <<",characteristic_abs_normal_flux,"<<refined.characteristicAbsFlux
            <<",solved_abs_normal_flux,"<<refined.solvedAbsFlux
            <<",solved_sharp_change_from_canonical,"<<refined.solvedPower-stats.solvedPower
            <<",solved_sharp_change_from_preceding,"<<refined.solvedPower-preceding.solvedPower
            <<",characteristic_sharp_change_from_preceding,"<<refined.characteristicPower-preceding.characteristicPower
            <<",area_change_from_preceding,"<<refined.area-preceding.area
            <<",observation_only,no_new_gate\n";
        preceding=refined;
    }
    if(!output.empty()){
        surface.flush();native.flush();provenance.flush();
        check(surface.good()&&native.good()&&provenance.good(),"CSF diagnosis source output write failed");
    }
    std::cout<<"CSF_COUPLING_STATIC_CIRCLE,begin\n"<<std::flush;
    csfCircleStatic(*mesh,metrics,config);
    std::cout<<"CSF_COUPLING_ONLY,actual_step_and_source_checks,PASS,power_and_delta_errors,OBSERVED";
    if(!output.empty())std::cout<<",output_directory,"<<output.string();
    std::cout<<'\n';
}


// Opt-in spatial diagnosis: compare the actually projected FirstJet scalar
// motion against its retained midpoint source. No equality/physics gate is
// manufactured, and no work is reconstructed from the production energy defect.
void projectedMotionWorkOnly(const std::string& outputDirectory)
{
    namespace fs=std::filesystem;
    fs::path output;
    if(!outputDirectory.empty()){
        output=fs::path(outputDirectory);
        check(output.is_absolute()&&(output.root_name()==fs::path("D:")||output.root_name()==fs::path("d:")),
            "projected-motion output directory must be absolute on D:");
        fs::create_directories(output);
    }
    int threads=1;
#ifdef _OPENMP
    omp_set_dynamic(0);omp_set_num_threads(8);threads=omp_get_max_threads();
#endif
    constexpr double dt=.0002,sigma=341.642;
    std::cout<<"PROJECTED_MOTION_CONFIG,nx,64,ny,64,a,.168,b,.15238095238095239,rho_liquid,100,rho_gas,1,"
        "mu_liquid,0,mu_gas,0,sigma,"<<sigma<<",dt,"<<dt<<",threads,"<<threads
        <<",density,cellAverageProjection,capillary,geometricAdjoint,time,adaptiveLinear,steps,1\n"<<std::flush;
    auto mesh=makeMesh(64);
    const auto phase=HydroPhaseGeometry::build(*mesh,stepMetrics(*mesh),csfInitialEllipse(*mesh,.168,.15238095238095239));
    auto accepted=HydroTwoPhaseState::create(*mesh,phase,stepVelocity(*mesh,phase.layout(),true));
    HydroTwoPhaseStepConfig config;config.liquidDensity=100;config.gasDensity=1;
    config.liquidViscosity=config.gasViscosity=0;config.surfaceTension=sigma;
    config.momentumDensityRepresentation=HydroMomentumMass::DensityRepresentation::CellAverageProjection;
    config.spaceTime.characteristic.maxStep=dt;config.spaceTime.initialPanels=1;
    config.onIteration=[](const HydroTwoPhaseIterationReport& r){
        std::cout<<"PROJECTED_MOTION_STAGE,"<<r.iteration<<",transport_seconds,"<<r.stageSeconds[1]
            <<",csf_seconds,"<<r.csfProjectionSeconds<<",saddle_seconds,"<<r.stageSeconds[6]<<'\n'<<std::flush;
    };
    const auto saved=capture(accepted.publication(),accepted.values());
    const auto step=HydroTwoPhaseStep::prepare(*mesh,accepted,{1,0,dt,Position::End},config);
    unchanged(saved,accepted.publication(),accepted.values());
    const auto& material=step.materialTransport();const auto& begin=material.beginGeometry();
    const auto& mid=material.midpointRate().geometry();const auto& end=material.endGeometry();
    const auto& actual=material.midpointRate().spatialAdvector();const auto& load=step.capillaryLoad();
    check(load.discretization()==HydroCsfDiscretization::GeometricAdjoint
        &&step.report().iterations<=2&&step.report().momentumDensityRepresentation==config.momentumDensityRepresentation,
        "projected-motion diagnosis changed current production method");
    load.validateFor(*mesh,step.saddle().pressureLoad(),mid.geometryIdentity(),mid.timeLayer(),
        actual.reconstruction(),actual.boundary());
    check(begin.geometryIdentity()==accepted.geometryIdentity()
        &&end.geometryIdentity()==step.candidateState().geometryIdentity()
        &&begin.timeLayer().time()==0&&mid.timeLayer().time()==dt/2&&end.timeLayer().time()==dt,
        "projected-motion diagnosis detached begin/mid/end owners");
    actual.validateFor(*mesh,actual.momentumGeometry(),mid.geometryIdentity(),mid.timeLayer());
    const auto solved=actual.reconstruction().evaluate(*mesh,step.solvedMidpointValues(),
        mid.geometryIdentity(),mid.timeLayer(),step.beginState().boundary());
    for(const auto* mass:std::array<const HydroMomentumMass*,3>{
        &step.dualMaterial().beginMass(),&step.dualMaterial().midpointMass(),&step.dualMaterial().endMass()})
        check(mass->representation()==config.momentumDensityRepresentation,"projected-motion diagnosis lost actual P0 mass layer");

    std::ofstream nodes,provenance;
    if(!output.empty()){
        nodes.open(output/"projected_motion_midpoint_nodes.csv");
        provenance.open(output/"projected_motion_sources.csv");
        check(nodes.good()&&provenance.good(),"projected-motion diagnostic output open failed");
        nodes<<std::setprecision(17);provenance<<std::setprecision(17);
        nodes<<"component,graph,node,native_cell,pressure_dof,x,y,nx,ny,kappa,kappa0,physical_weight,"
            "phi_begin,phi_mid,phi_end,gradient_mid_norm,v_projected_secant,n_dot_transport_midpoint,"
            "n_dot_solved_midpoint,secant_minus_transport,secant_minus_solved\n";
        provenance<<"field,geometry_identity,advector_identity,step,begin,end,position,physical_time\n";
    }
    const auto sourceRow=[&](const char* name,HermiteGalsIdentity geometry,HydroAdvectorIdentity velocity,HermiteGalsTimeLayer clock){
        std::cout<<"PROJECTED_MOTION_SOURCE,";csfSourceRow(std::cout,name,geometry,velocity,clock);
        if(provenance.is_open())csfSourceRow(provenance,name,geometry,velocity,clock);
    };
    sourceRow("begin_scalar",begin.geometryIdentity(),{},begin.timeLayer());
    sourceRow("mid_scalar_and_force",mid.geometryIdentity(),{},mid.timeLayer());
    sourceRow("end_scalar",end.geometryIdentity(),{},end.timeLayer());
    sourceRow("actual_transport_midpoint",actual.geometryIdentity(),actual.identity(),actual.timeLayer());
    sourceRow("solved_midpoint_same_reconstruction",solved.geometryIdentity(),solved.identity(),solved.timeLayer());

    // Replay the accepted adjoint GL8 subarcs exactly once per real graph.
    // Native-row ledgers repeat a panel; deduplicate instead of double counting.
    using GraphKey=std::pair<std::size_t,std::size_t>;
    std::map<GraphKey,std::set<std::pair<double,double>>> panels;
    for(int u=0;u<mid.layout().velocityCount();++u)
        for(const auto& panel:load.quadraturePanels(u)){
            check(panel.coarseSubdivisions==1,"projected-motion diagnosis requires actual two-subarc adjoint rule");
            panels[{panel.component,panel.graph}].insert({panel.begin,panel.end});
        }
    std::map<std::string,long double> m;
    m["nx"]=m["ny"]=64;m["rho_liquid"]=100;m["rho_gas"]=1;m["mu_liquid"]=m["mu_gas"]=0;
    m["sigma"]=sigma;m["dt"]=dt;m["threads"]=threads;m["accepted_ac_stages"]=step.report().iterations;
    m["raw_kappa_min"]=std::numeric_limits<double>::infinity();
    m["raw_kappa_max"]=-std::numeric_limits<double>::infinity();
    const auto& sb=begin.publication().state().geometry().snapshot();
    const auto& sm=mid.publication().state().geometry().snapshot();
    const auto& se=end.publication().state().geometry().snapshot();
    const auto& publication=mid.publication();
    for(std::size_t c=0;c<publication.topology().curves.size();++c)
        for(std::size_t g=0;g<publication.topology().curves[c].graphs().size();++g){
            const auto found=panels.find({c,g});
            check(found!=panels.end()&&!found->second.empty(),"midpoint graph missing actual accepted force panels");
            std::vector<std::pair<double,double>> intervals(found->second.begin(),found->second.end());
            double covered=0;
            for(const auto& interval:intervals){
                check(interval.first==covered&&interval.second>interval.first,"accepted force graph partition is incomplete/overlapping");
                covered=interval.second;
            }
            check(covered==1,"accepted force graph partition does not cover actual graph");
            std::size_t node=0;
            for(const auto& interval:intervals){
                const double middle=interval.first+(interval.second-interval.first)/2;
                for(const auto& arc:std::array<std::pair<double,double>,2>{{{interval.first,middle},{middle,interval.second}}})
                    for(const auto& q:publication.quadrature(c,g,arc.first,arc.second)){
                        const auto& f=q.sample.frame;const int cell=f.nativeCellIndex,p=mid.layout().pressureDof(cell);
                        check(p>=0&&q.sample.geometryIdentity==mid.geometryIdentity()
                            &&sameClock(q.sample.timeLayer,mid.timeLayer())&&q.physicalWeight>0,
                            "projected-motion node lost actual midpoint surface owner");
                        const auto jb=sb.scalarJetInNativeCell(cell,f.point),jm=sm.scalarJetInNativeCell(cell,f.point),
                            je=se.scalarJetInNativeCell(cell,f.point);
                        const long double gradient=std::hypot(jm.dx,jm.dy),w=q.physicalWeight,k=q.sample.curvature,
                            k0=load.referenceCurvatures()[p];
                        check(gradient>0&&std::isfinite(gradient)&&std::isfinite(k)&&std::isfinite(k0),
                            "projected-motion node has singular/nonfinite scalar/curvature");
                        const long double vp=-(static_cast<long double>(je.value)-jb.value)/(dt*gradient);
                        const auto uc=actual.cellVelocity(p,f.point),us=solved.cellVelocity(p,f.point);
                        const long double vc=f.normal.x*uc.x+f.normal.y*uc.y,vs=f.normal.x*us.x+f.normal.y*us.y;
                        check(std::isfinite(vp)&&std::isfinite(vc)&&std::isfinite(vs),"nonfinite projected motion or source velocity");
                        m["surface_area_mid"]+=w;m["surface_nodes"]+=1;
                        m["raw_kappa_min"]=std::min(m["raw_kappa_min"],k);m["raw_kappa_max"]=std::max(m["raw_kappa_max"],k);
                        m["projected_surface_energy_rate"]+=sigma*w*k*vp;
                        m["transport_surface_energy_rate"]+=sigma*w*k*vc;
                        m["solved_surface_energy_rate"]+=sigma*w*k*vs;
                        m["projected_volume_flux"]+=w*vp;m["transport_volume_flux"]+=w*vc;m["solved_volume_flux"]+=w*vs;
                        m["projected_kappa0_volume_power"]+=sigma*w*k0*vp;
                        m["transport_kappa0_volume_power"]+=sigma*w*k0*vc;
                        m["solved_kappa0_volume_power"]+=sigma*w*k0*vs;
                        m["projected_speed_square_integral"]+=w*vp*vp;
                        m["transport_speed_square_integral"]+=w*vc*vc;m["solved_speed_square_integral"]+=w*vs*vs;
                        m["projected_minus_transport_square_integral"]+=w*(vp-vc)*(vp-vc);
                        m["projected_minus_solved_square_integral"]+=w*(vp-vs)*(vp-vs);
                        m["projected_minus_transport_maxabs"]=std::max(m["projected_minus_transport_maxabs"],std::abs(vp-vc));
                        m["projected_minus_solved_maxabs"]=std::max(m["projected_minus_solved_maxabs"],std::abs(vp-vs));
                        m["midpoint_scalar_maxabs"]=std::max(m["midpoint_scalar_maxabs"],std::abs(static_cast<long double>(jm.value)));
                        if(nodes.is_open())nodes<<c<<','<<g<<','<<node<<','<<cell<<','<<p<<','<<f.point.x<<','<<f.point.y<<','
                            <<f.normal.x<<','<<f.normal.y<<','<<k<<','<<k0<<','<<w<<','<<jb.value<<','<<jm.value<<','<<je.value<<','
                            <<gradient<<','<<vp<<','<<vc<<','<<vs<<','<<vp-vc<<','<<vp-vs<<'\n';
                        ++node;
                    }
            }
        }
    check(m["surface_nodes"]>0&&m["surface_area_mid"]>0,"projected-motion diagnosis has no real surface samples");
    m["surface_area_begin"]=actualInterfaceMeasure(begin);m["surface_area_end"]=actualInterfaceMeasure(end);
    m["measured_surface_energy_rate"]=sigma*(m["surface_area_end"]-m["surface_area_begin"])/dt;
    long double volumeBefore=0,volumeAfter=0;
    for(int p=0;p<mid.layout().pressureCount();++p){
        volumeBefore+=begin.cellVolumes()[p].liquid;volumeAfter+=end.cellVolumes()[p].liquid;
        m["measured_kappa0_volume_power"]+=sigma*load.referenceCurvatures()[p]
            *(end.cellVolumes()[p].liquid-begin.cellVolumes()[p].liquid)/dt;
    }
    m["volume_begin"]=volumeBefore;m["volume_end"]=volumeAfter;m["measured_volume_change_rate"]=(volumeAfter-volumeBefore)/dt;
    for(int u=0;u<mid.layout().velocityCount();++u){
        const long double vc=actual.values().normalVelocity()[u],vs=step.solvedMidpointValues().normalVelocity()[u];
        m["native_transport_capillary_work"]+=vc*load.momentumRhsLoad()[u];
        m["native_solved_capillary_work"]+=vs*load.momentumRhsLoad()[u];
        m["native_transport_reference_work"]+=vc*load.referencePressureLoad()[u];
        m["native_solved_reference_work"]+=vs*load.referencePressureLoad()[u];
        m["native_transport_residual_work"]+=vc*load.adjointResidualLoad()[u];
        m["native_solved_residual_work"]+=vs*load.adjointResidualLoad()[u];
    }
    m["native_transport_D_max"]=csfNativeD(*mesh,mid,actual.values());
    m["native_solved_D_max"]=csfNativeD(*mesh,mid,step.solvedMidpointValues());
    m["kinetic_begin"]=step.report().kineticBegin;m["kinetic_end"]=step.report().kineticEnd;
    m["kinetic_change_rate"]=(m["kinetic_end"]-m["kinetic_begin"])/dt;
    m["surface_rate_minus_projected_work"]=m["measured_surface_energy_rate"]-m["projected_surface_energy_rate"];
    m["projected_minus_transport_surface_work"]=m["projected_surface_energy_rate"]-m["transport_surface_energy_rate"];
    m["projected_minus_solved_surface_work"]=m["projected_surface_energy_rate"]-m["solved_surface_energy_rate"];
    m["transport_adjoint_power_residual"]=m["native_transport_capillary_work"]+m["transport_surface_energy_rate"]
        -m["native_transport_reference_work"]-m["transport_kappa0_volume_power"];
    m["solved_adjoint_power_residual"]=m["native_solved_capillary_work"]+m["solved_surface_energy_rate"]
        -m["native_solved_reference_work"]-m["solved_kappa0_volume_power"];
    m["surface_rate_plus_native_solved_work"]=m["measured_surface_energy_rate"]+m["native_solved_capillary_work"];
    m["kinetic_rate_minus_native_solved_work"]=m["kinetic_change_rate"]-m["native_solved_capillary_work"];
    m["total_energy_rate"]=m["kinetic_change_rate"]+m["measured_surface_energy_rate"];
    for(const auto* key:{"projected_minus_transport","projected_minus_solved"}){
        const auto integral=m[std::string(key)+"_square_integral"];
        m[std::string(key)+"_L2"]=std::sqrt(integral);
        m[std::string(key)+"_rms"]=std::sqrt(integral/m["surface_area_mid"]);
    }
    for(const auto* key:{"projected","transport","solved"})
        m[std::string(key)+"_normal_speed_rms"]=std::sqrt(m[std::string(key)+"_speed_square_integral"]/m["surface_area_mid"]);
    for(const auto& row:m){check(std::isfinite(row.second),"projected-motion diagnostic generated nonfinite summary");
        std::cout<<"PROJECTED_MOTION_METRIC,"<<row.first<<','<<row.second<<'\n';}
    if(!output.empty()){
        std::ofstream json(output/"projected_motion_work.json");check(json.good(),"projected-motion JSON open failed");
        json<<std::setprecision(17)<<"{\n  \"status\": \"OBSERVED_SINGLE_STEP\",\n"
            "  \"surface_energy_power_sign\": \"+sigma integral(kappa v ds); physical capillary power has opposite sign\",\n"
            "  \"surface_rule\": \"exact replay of unique actual adjoint accepted panels; two GL8 subarcs per panel\",\n"
            "  \"endpoint_measure_rule\": \"independently sum positive original begin/end publication GL8 physical weights\",\n"
            "  \"temporal_limit\": \"central endpoint secant at t_mid has O(dt^2) truncation for a temporally smooth scalar family; spatial projection and interval stage changes need not be smooth\",\n"
            "  \"cause_limit\": \"single-step direct observation does not establish a unique cause or full-cycle validity\",\n"
            "  \"metrics\": {\n";
        std::size_t i=0;for(const auto& row:m)json<<"    \""<<row.first<<"\": "<<row.second<<(++i<m.size()?",":"")<<'\n';
        json<<"  }\n}\n";nodes.flush();provenance.flush();json.flush();
        check(nodes.good()&&provenance.good()&&json.good(),"projected-motion diagnostic output write failed");
    }
    std::cout<<"PROJECTED_MOTION_WORK_ONLY,source_checks,PASS,power_and_motion_errors,OBSERVED,single_step,1\n"<<std::flush;
}

void heightFunctionIncrementWindow(const std::string& directory,int steps,double requestedDt=.0002)
{
    namespace fs=std::filesystem;const fs::path output(directory);
    check(output.is_absolute()&&(output.root_name()==fs::path("D:")||output.root_name()==fs::path("d:")),"capillary research output must be absolute on D:");
    check(steps>=1&&steps<=20,"capillary research window requires 1..20 fixed steps");
    constexpr const char* model="heightFunctionIncrement";
    fs::create_directories(output);check(!fs::exists(output/"metrics.csv"),"capillary research refuses to overwrite observation");
#ifdef _OPENMP
    omp_set_dynamic(0);omp_set_num_threads(8);
#endif
    constexpr double a=.168,b=.15238095238095239,sigma=341.642;const double dt=requestedDt;
    check(std::isfinite(dt)&&dt>0,"capillary research requires finite positive explicit dt");
    auto mesh=makeMesh(64);const auto phase=HydroPhaseGeometry::build(*mesh,stepMetrics(*mesh),
        csfInitialEllipse(*mesh,a,b,HermiteGalsGradientUpdate::CompatibleIncrement));
    auto accepted=HydroTwoPhaseState::create(*mesh,phase,stepVelocity(*mesh,phase.layout(),true));
    HydroTwoPhaseStepConfig config;config.liquidDensity=100;config.gasDensity=1;config.liquidViscosity=config.gasViscosity=0;
    config.surfaceTension=sigma;config.momentumDensityRepresentation=HydroMomentumMass::DensityRepresentation::CellAverageProjection;
    config.publication.gradientUpdate=HermiteGalsGradientUpdate::CompatibleIncrement;
    config.publication.scalarRepresentation=HermiteGalsScalarRepresentation::SingleQuartic;
    config.spaceTime.characteristic.maxStep=dt;config.spaceTime.initialPanels=1;
    config.capillaryModel=HydroCapillaryModel::BalancedCsf;
    config.csfCurvature.forceMethod=HydroCsfForceMethod::HeightFunctionIncrement;
    config.onIteration=[](const HydroTwoPhaseIterationReport& r){std::cout<<"CAPILLARY_RESEARCH_STAGE,step,"<<r.timeLayer.step
        <<",iteration,"<<r.iteration<<",transport,"<<r.stageSeconds[1]<<",capillary,"<<r.stageSeconds[4]
        <<",saddle,"<<r.stageSeconds[6]<<'\n'<<std::flush;};
    std::ofstream metrics(output/"metrics.csv");check(metrics.good(),"capillary research metrics open failed");metrics<<std::setprecision(17);
    metrics<<"step,time,geometry_identity,midpoint_identity,model,AC_stages,liquid_volume,relative_volume_drift,surface_measure,kinetic,total_energy,total_energy_change,umax,divergence_max,raw_kappa_abs_max,raw_kappa_endpoint_jump,raw_kappa_graph_d2_rms,raw_kappa_graph_d2_max,capillary_power,capillary_energy_defect,frozen_momentum_defect,momentum_power_relative,step_seconds,endpoint_seconds,transport_seconds,mass_seconds,convection_seconds,capillary_seconds,viscosity_seconds,saddle_seconds,reassessment_seconds,output_seconds\n";
    long double initialVolume=0,initialEnergy=0;for(const auto& q:phase.cellVolumes())initialVolume+=q.liquid;
    initialEnergy=sigma*actualInterfaceMeasure(phase);
    const auto record=[&](int index,const HydroTwoPhaseStep* step,double stepSeconds){
        const auto outputStarted=std::chrono::steady_clock::now();const auto& publication=accepted.publication();const auto& geometry=accepted.phaseGeometry();
        const auto& layout=geometry.layout();const auto clock=accepted.timeLayer();const auto identity=accepted.geometryIdentity();
        std::ostringstream suffix;suffix<<std::setfill('0')<<std::setw(6)<<index;
        std::vector<HermiteGalsSurfaceSample> nodes;std::vector<std::array<std::size_t,2>> segments;
        long double secondSquare=0;std::size_t secondCount=0;double maxCurvature=0,maxJump=0,maxSecond=0;
        for(std::size_t c=0;c<publication.topology().curves.size();++c){const auto& curve=publication.topology().curves[c];
            check(curve.closed(),"oscillating-drop research source unexpectedly open");
            for(std::size_t g=0;g<curve.graphs().size();++g){const auto first=nodes.size();
                for(int j=0;j<=8;++j){auto sample=publication.sample(c,g,j/8.);check(sample.geometryIdentity==identity&&sameClock(sample.timeLayer,clock),"research sample lost actual accepted clock");
                    maxCurvature=std::max(maxCurvature,std::abs(sample.curvature));nodes.push_back(std::move(sample));}
                for(int j=0;j<8;++j)segments.push_back({{first+j+1,first+j+2}});
                for(int j=1;j<8;++j){const double d=nodes[first+j-1].curvature-2*nodes[first+j].curvature+nodes[first+j+1].curvature;
                    secondSquare+=static_cast<long double>(d)*d;++secondCount;maxSecond=std::max(maxSecond,std::abs(d));}
                const auto next=(g+1)%curve.graphs().size();maxJump=std::max(maxJump,std::abs(nodes.back().curvature-publication.sample(c,next,0).curvature));
            }
        }
        std::ofstream surface(output/("interface_"+suffix.str()+".dat"));check(surface.good(),"research interface open failed");surface<<std::setprecision(17);
        surface<<"TITLE=\"Actual accepted FirstJet raw interface\"\nVARIABLES=\"X\" \"Y\" \"Nx\" \"Ny\" \"Kappa\"\nZONE T=\""<<model<<" step "<<index<<"\", NODES="<<nodes.size()
            <<", ELEMENTS="<<segments.size()<<", DATAPACKING=BLOCK, ZONETYPE=FELINESEG, SOLUTIONTIME="<<clock.time()
            <<"\nAUXDATA geometry_identity=\""<<identity.value()<<"\"\nAUXDATA full_clock=\""<<clock.step<<','<<clock.begin<<','<<clock.end<<','<<static_cast<int>(clock.position)<<"\"\n";
        for(int field=0;field<5;++field){for(const auto& sample:nodes){const auto& f=sample.frame;surface<<(field==0?f.point.x:field==1?f.point.y:field==2?f.normal.x:field==3?f.normal.y:sample.curvature)<<'\n';}}
        for(const auto& segment:segments)surface<<segment[0]<<' '<<segment[1]<<'\n';surface.flush();check(surface.good(),"research interface write failed");
        std::ofstream native(output/("native_"+suffix.str()+".csv"));check(native.good(),"research native fields open failed");native<<std::setprecision(17);
        native<<"step,time,geometry_identity,u,native_face,axis,x,y,normal_velocity,negative_pressure,positive_pressure,midpoint_identity,midpoint_time,midpoint_capillary_load\n";
        double umax=0;for(int u=0;u<layout.velocityCount();++u){const auto face=layout.domain().face(*mesh,layout.nativeFace(u));
            const auto& edge=face.horizontal?mesh->hEdges(face.edgeLabel):mesh->vEdges(face.edgeLabel);
            const auto p=(mesh->verts(edge.edgeToVert[0]).coord()+mesh->verts(edge.edgeToVert[1]).coord())*.5;
            const auto velocity=accepted.values().normalVelocity()[u];umax=std::max(umax,std::abs(velocity));
            native<<index<<','<<clock.time()<<','<<identity.value()<<','<<u<<','<<layout.nativeFace(u)<<','<<(face.horizontal?1:0)<<','<<p.x<<','<<p.y<<','<<velocity<<',';
            if(face.negativeCell>=0)native<<accepted.values().pressure()[layout.pressureDof(face.negativeCell)];else native<<"NA";native<<',';
            if(face.positiveCell>=0)native<<accepted.values().pressure()[layout.pressureDof(face.positiveCell)];else native<<"NA";native<<',';
            if(step){const auto& midpoint=step->materialTransport().midpointRate().geometry();const double combined=step->capillaryMomentumRhsLoad()[u];
                native<<midpoint.geometryIdentity().value()<<','<<midpoint.timeLayer().time()<<','<<combined;}
            else native<<"NA,NA,NA";native<<'\n';
        }
        native.flush();check(native.good(),"research native fields write failed");
        long double liquid=0;for(const auto& q:geometry.cellVolumes())liquid+=q.liquid;
        const auto area=actualInterfaceMeasure(geometry);const auto kinetic=step?step->report().kineticEnd:0;
        const auto energy=kinetic+sigma*area;double divergence=0;for(double d:geometry.metrics().apply(*mesh,accepted.values()))divergence=std::max(divergence,std::abs(d));
        std::array<double,8> stageSeconds{};if(step)for(const auto& stage:step->report().history)for(int k=0;k<8;++k)stageSeconds[k]+=stage.stageSeconds[k];
        const double outputSeconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-outputStarted).count();
        metrics<<index<<','<<clock.time()<<','<<identity.value()<<','<<(step?step->materialTransport().midpointRate().geometry().geometryIdentity().value():0)<<','<<model<<','<<(step?step->report().iterations:0)
            <<','<<liquid<<','<<(liquid-initialVolume)/initialVolume<<','<<area<<','<<kinetic<<','<<energy<<','<<energy-initialEnergy<<','<<umax<<','<<divergence<<','<<maxCurvature<<','<<maxJump
            <<','<<(secondCount?std::sqrt(secondSquare/secondCount):0)<<','<<maxSecond<<','<<(step?step->report().capillaryPower:0)<<','<<(step?step->report().capillaryEnergyDefect:0)
            <<','<<(step?step->report().frozenOwnerMomentumDefect:0)<<','<<(step?step->report().momentumPowerRelativeDefect:0)<<','<<stepSeconds;
        for(double seconds:stageSeconds)metrics<<','<<seconds;metrics<<','<<outputSeconds<<'\n';metrics.flush();check(metrics.good(),"research metrics write failed");
        // Created only after all fields of this accepted layer are complete.
        std::ofstream marker(output/("accepted_"+suffix.str()+".marker"));marker<<index<<','<<clock.time()<<','<<identity.value()<<'\n';marker.flush();check(marker.good(),"research accepted marker failed");
        std::cout<<"CAPILLARY_RESEARCH_ACCEPTED,model,"<<model<<",step,"<<index<<",time,"<<clock.time()<<",AC,"<<(step?step->report().iterations:0)
            <<",raw_kappa_max,"<<maxCurvature<<",seam_jump,"<<maxJump<<",dV/V,"<<(liquid-initialVolume)/initialVolume<<",dE,"<<energy-initialEnergy
            <<",umax,"<<umax<<",wall,"<<stepSeconds<<",output,"<<outputSeconds<<'\n'<<std::flush;
    };
    std::cout<<"CAPILLARY_RESEARCH_BEGIN,n,64,rho,100:1,mu,0,sigma,"<<sigma<<",dt,"<<dt<<",threads,8,gradient,compatibleIncrement,scalar,SingleQuartic,density,cellAverageProjection,AC,adaptiveLinear,model,"<<model<<",steps,"<<steps<<'\n'<<std::flush;
    record(0,nullptr,0);
    for(int index=1;index<=steps;++index){const auto started=std::chrono::steady_clock::now();const auto retained=accepted;const auto saved=capture(accepted.publication(),accepted.values());
        try{const auto step=HydroTwoPhaseStep::prepare(*mesh,accepted,{index,(index-1)*dt,index*dt,Position::End},config);
            unchanged(saved,accepted.publication(),accepted.values());check(accepted.sameIdentity(retained),"capillary research prepare changed accepted state");
            check(step.report().iterations<=2&&step.report().timeScheme==HydroTwoPhaseTimeScheme::AdaptiveLinearMidpoint
                &&step.report().capillaryModel==config.capillaryModel&&step.report().momentumDensityRepresentation==config.momentumDensityRepresentation
                &&step.candidateState().publication().state().primary().gradientUpdate==HermiteGalsGradientUpdate::CompatibleIncrement,
                "capillary research changed common AC/material/gradient route");
            step.commit(*mesh,accepted);record(index,&step,std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count());
        }catch(const std::exception& e){unchanged(saved,retained.publication(),retained.values());
            std::ofstream failure(output/"failure.txt");failure<<"attempted_step="<<index<<" accepted_step="<<accepted.timeLayer().step<<" accepted_time="<<accepted.timeLayer().time()<<" reason="<<e.what()<<'\n';
            throw;}
    }
    std::ofstream completed(output/"completed.txt");completed<<"model="<<model<<" steps="<<steps<<" end_time="<<accepted.timeLayer().time()<<" observation_only=1\n";
    std::cout<<"CAPILLARY_RESEARCH_COMPLETED,model,"<<model<<",accepted_steps,"<<steps<<",fixed_dt,"<<dt<<",raw_geometry_observation_only,no_long_time_stability_claim\n"<<std::flush;
}

void stationary(const QuadtreeMesh& mesh,const HydroPhaseGeometry& geometry)
{
    auto accepted=HydroTwoPhaseState::create(mesh,geometry,stepVelocity(mesh,geometry.layout(),true));
    HydroTwoPhaseStepConfig config;config.timeScheme=HydroTwoPhaseTimeScheme::FullyCoupledMidpoint;
    config.liquidDensity=1e6;config.spaceTime.characteristic.maxStep=.00025;config.spaceTime.initialPanels=1;
    const HermiteGalsTimeLayer clock{1,0,.001,Position::End};
    const auto saved=capture(accepted.publication(),accepted.values());std::vector<HydroTwoPhaseIterationReport> observed;
    config.onIteration=[&](const HydroTwoPhaseIterationReport& r){unchanged(saved,accepted.publication(),accepted.values());
        inspectIteration(r,static_cast<int>(observed.size())+1,clock,
                observed.empty()?static_cast<std::size_t>(config.spaceTime.initialPanels):observed.back().acceptedPanels);observed.push_back(r);};
    const auto step=HydroTwoPhaseStep::prepare(mesh,accepted,clock,config);
    inspectIterationObserver(step,observed);inspectStep(mesh,step,false);
    capillaryRejectsAction([&]{step.capillaryLoad();},"step has no interface load");
    check(step.candidateState().values().normalVelocity()==accepted.values().normalVelocity()
        &&step.candidateState().values().pressure()==accepted.values().pressure(),"stationary sigma0 native state drifted");
    stepRejects([&]{HydroTwoPhaseState::create(mesh,step.materialTransport().midpointRate().geometry(),step.solvedMidpointValues());});
    const auto& a=accepted.publication().state();const auto& b=step.candidateState().publication().state();
    check(a.primary().vertices.size()==b.primary().vertices.size()&&a.primary().edges.size()==b.primary().edges.size()
        &&a.primary().geometryCells==b.primary().geometryCells
        &&a.closure().vertexMixedDerivatives.size()==b.closure().vertexMixedDerivatives.size()
        &&a.closure().cellBubbleCoefficients.size()==b.closure().cellBubbleCoefficients.size(),"stationary FirstJet layout changed");
    long double geometryGap=0,beginNorm=0,endNorm=0,vertexGap=0,edgeGap=0,closureGap=0;
    const auto compare=[&](double av,double bv,long double& maximum){const long double delta=static_cast<long double>(bv)-av;
        geometryGap=std::hypot(geometryGap,delta);beginNorm=std::hypot(beginNorm,static_cast<long double>(av));
        endNorm=std::hypot(endNorm,static_cast<long double>(bv));maximum=std::max(maximum,std::abs(delta));};
    for(std::size_t i=0;i<a.primary().vertices.size();++i){const auto& av=a.primary().vertices[i];const auto& bv=b.primary().vertices[i];
        check(av.vertex==bv.vertex,"stationary FirstJet vertex identity changed");
        compare(av.phi,bv.phi,vertexGap);compare(av.dx,bv.dx,vertexGap);compare(av.dy,bv.dy,vertexGap);}
    for(std::size_t i=0;i<a.primary().edges.size();++i){const auto& av=a.primary().edges[i];const auto& bv=b.primary().edges[i];
        check(av.vertices==bv.vertices,"stationary FirstJet edge identity changed");compare(av.phi,bv.phi,edgeGap);compare(av.dn,bv.dn,edgeGap);}
    for(std::size_t i=0;i<a.closure().vertexMixedDerivatives.size();++i)
        compare(a.closure().vertexMixedDerivatives[i],b.closure().vertexMixedDerivatives[i],closureGap);
    for(std::size_t i=0;i<a.closure().cellBubbleCoefficients.size();++i)
        compare(a.closure().cellBubbleCoefficients[i],b.closure().cellBubbleCoefficients[i],closureGap);
    const auto geometryScale=std::max(beginNorm,endNorm);
    const auto relative=geometryScale==0?(geometryGap==0?0:std::numeric_limits<long double>::infinity()):geometryGap/geometryScale;
    std::cout<<"STATIONARY_GEOMETRY,"<<vertexGap<<','<<edgeGap<<','<<closureGap<<','<<geometryGap<<','<<geometryScale<<','<<relative<<'\n';
    peakEquation=std::max(peakEquation,scaled(geometryGap,0,geometryScale,"stationary FirstJet geometry norm",1e-12));
    auto differentHandle=HydroTwoPhaseState::create(mesh,geometry,accepted.values());const auto differentCopy=capture(differentHandle.publication(),differentHandle.values());
    stepRejects([&]{step.commit(mesh,differentHandle);});unchanged(differentCopy,differentHandle.publication(),differentHandle.values());
    step.commit(mesh,accepted);++acceptedSteps;
}
void capillaryModelOnly()
{
    auto mesh=makeMesh(4);Case shape;shape.pureFirstJet=true;
    const auto phase=HydroPhaseGeometry::build(*mesh,stepMetrics(*mesh),initial(*mesh,shape));
    const auto accepted=HydroTwoPhaseState::create(*mesh,phase,stepVelocity(*mesh,phase.layout()));
    const auto saved=capture(accepted.publication(),accepted.values());
    HydroTwoPhaseStepConfig cfg;cfg.liquidDensity=100;cfg.gasDensity=1;cfg.surfaceTension=1;
    cfg.spaceTime.initialPanels=1;cfg.spaceTime.characteristic.maxStep=.00025;
    const HermiteGalsTimeLayer clock{1,0,.00025,Position::End};
    int models=0,zeros=0,refusals=0;std::optional<HydroTwoPhaseStep> zeroReference;
    for(auto model:{HydroCapillaryModel::BalancedCsf,HydroCapillaryModel::SharpNormalStress}){
        cfg.capillaryModel=model;
        const auto step=HydroTwoPhaseStep::prepare(*mesh,accepted,clock,cfg);inspectStep(*mesh,step,true);
        check(step.capillaryModel()==model&&step.report().capillaryModel==model
            &&step.hasCapillaryLoad()&&step.report().iterations<=2,"selected force model changed AC path");
        long double norm=0;for(double f:step.capillaryMomentumRhsLoad())norm=std::hypot(norm,f);
        check(norm>0,"selected force model erased its actual nonzero-sigma RHS");++models;
        if(model==HydroCapillaryModel::BalancedCsf){
            check(&step.capillaryMomentumRhsLoad()==&step.capillaryLoad().momentumRhsLoad(),"CSF common getter replaced physical force");
            rejects([&]{step.sharpNormalStressLoad();},"no sharpNormalStress load");
        }else rejects([&]{step.capillaryLoad();},"no balancedCsf load");++refusals;
        auto zeroCfg=cfg;zeroCfg.surfaceTension=0;
        const auto zero=HydroTwoPhaseStep::prepare(*mesh,accepted,clock,zeroCfg);inspectStep(*mesh,zero,true);
        check(!zero.hasCapillaryLoad()&&zero.report().capillaryPower==0
            &&zero.report().interfaceMeasureBegin==0&&zero.report().interfaceMeasureEnd==0,
            "zero sigma selected a force factory or sampled surface energy");
        rejects([&]{zero.capillaryMomentumRhsLoad();},"step has no interface load");++refusals;++zeros;
        if(zeroReference){
            check(zero.integralMomentumRhs()==zeroReference->integralMomentumRhs()
                &&zero.candidateState().values().normalVelocity()==zeroReference->candidateState().values().normalVelocity()
                &&zero.candidateState().values().pressure()==zeroReference->candidateState().values().pressure(),
                "zero sigma method selection changed shared momentum arithmetic");
            const Matrix gap=zero.velocityMatrix()-zeroReference->velocityMatrix();
            check(gap.norm()==0,"zero sigma method selection changed shared matrix");
        }else zeroReference=zero;
        unchanged(saved,accepted.publication(),accepted.values());
    }
    cfg.capillaryModel=static_cast<HydroCapillaryModel>(99);
    rejects([&]{HydroTwoPhaseStep::prepare(*mesh,accepted,clock,cfg);},"invalid capillary model");++refusals;
    unchanged(saved,accepted.publication(),accepted.values());
    check(models==2&&zeros==2&&refusals==5,"incomplete independent method selection coverage");
    std::cout<<"CAPILLARY_MODEL_ONLY,PASS models="<<models<<" sigma0="<<zeros<<" rejected="<<refusals
        <<" equation="<<peakEquation<<" energy="<<peakEnergy<<" material="<<peakMaterial<<'\n';
}

void projectedDensityOnly()
{
    auto mesh=makeMesh(4);Case shape;shape.pureFirstJet=true;
    const auto phase=HydroPhaseGeometry::build(*mesh,stepMetrics(*mesh),initial(*mesh,shape));
    auto accepted=HydroTwoPhaseState::create(*mesh,phase,stepVelocity(*mesh,phase.layout()));
    HydroTwoPhaseStepConfig cfg;cfg.liquidDensity=100;cfg.gasDensity=1;cfg.surfaceTension=0;
    cfg.spaceTime.initialPanels=1;cfg.spaceTime.characteristic.maxStep=.00025;
    const HermiteGalsTimeLayer first{1,0,.00025,Position::End};
    const auto prior=HydroTwoPhaseStep::prepare(*mesh,accepted,first,cfg);
    check(prior.dualMaterial().endMass().representation()==HydroMomentumMass::DensityRepresentation::SharpGeometry,
        "default density no longer uses actual sharp geometry");
    prior.commit(*mesh,accepted);const auto retained=accepted;const auto saved=capture(accepted.publication(),accepted.values());
    cfg.surfaceTension=1;cfg.momentumDensityRepresentation=HydroMomentumMass::DensityRepresentation::CellAverageProjection;
    const HermiteGalsTimeLayer second{2,.00025,.0005,Position::End};
    const auto step=HydroTwoPhaseStep::prepare(*mesh,accepted,second,cfg);inspectStep(*mesh,step,true);
    const auto& dual=step.dualMaterial();
    check(step.report().momentumDensityRepresentation==cfg.momentumDensityRepresentation
        &&step.report().iterations<=2&&step.hasCapillaryLoad()
        &&step.materialTransport().liquidDensity()==100&&step.materialTransport().gasDensity()==1,
        "projected step changed physical density or shared AC/capillary path");
    for(const auto* mass:std::array<const HydroMomentumMass*,3>{{&dual.beginMass(),&dual.midpointMass(),&dual.endMass()}}){
        check(mass->representation()==HydroMomentumMass::DensityRepresentation::CellAverageProjection,
            "projected three-layer transaction reused a sharp mass cache");
        mass->validateFor(*mesh,mass->phaseGeometry().geometryIdentity(),mass->phaseGeometry().timeLayer());
        rejects([&]{mass->phaseVolumes();},"actual sharp geometry");
        rejects([&]{mass->overlapPhaseVolumes();},"actual sharp geometry");
    }
    check(&dual.beginMass().phaseGeometry().cellVolumes()==&accepted.phaseGeometry().cellVolumes(),
        "representation change replaced accepted begin sharp material geometry");
    long double norm=0;for(double f:step.capillaryMomentumRhsLoad())norm=std::hypot(norm,f);
    check(norm>0,"projected-density nonzero-sigma fixture erased physical force");
    check(accepted.sameIdentity(retained),"projected preparation changed accepted state");unchanged(saved,accepted.publication(),accepted.values());
    auto invalid=cfg;invalid.momentumDensityRepresentation=static_cast<HydroMomentumMass::DensityRepresentation>(99);
    rejects([&]{HydroTwoPhaseStep::prepare(*mesh,accepted,second,invalid);},"invalid momentum density representation");
    unchanged(saved,accepted.publication(),accepted.values());
    std::cout<<"PROJECTED_DENSITY_ONLY,PASS rho=100:1 mu=0 sigma=1 steps=2 cache_rep_checked=1"
        <<" equation="<<peakEquation<<" energy="<<peakEnergy<<" material="<<peakMaterial<<'\n';
}

void vofTransportProviderOnly()
{
    namespace I=LAMP::Interface;using Time=LAMP::TimeLayer;
    auto mesh=makeMesh(4);const auto domain=domainFor(*mesh,false);std::vector<double> alpha(mesh->cellIndexNum());
    for(int c=0;c<mesh->cellIndexNum();++c){const auto& cell=mesh->cells(mesh->cellIndexToLabel[c]);
        const auto lo=mesh->verts(cell.cellToVert[0]).coord(),hi=mesh->verts(cell.cellToVert[2]).coord();
        alpha[c]=std::clamp((.37-lo.y)/(hi.y-lo.y),0.,1.);}
    const auto source=I::GeometricVofPublication::prepareLocalPlic(*mesh,domain,{},alpha);
    const auto layout=HydroDofLayout::build(*mesh,source.domain());
    const auto metrics=HydroPhysicalDivergence::build(*mesh,HydroIntegralDivergence::build(*mesh,layout));
    const auto phase=HydroPhaseGeometry::build(*mesh,metrics,source);const auto zero=stepVelocity(*mesh,layout,true);
    const auto accepted=HydroTwoPhaseState::createVof(*mesh,phase,zero),saved=accepted;
    const auto oldCells=source.plicCells();const Time end{1,0,1e-5,Time::Position::End};
    const auto unchangedAccepted=[&]{check(accepted.sameIdentity(saved)&&accepted.phaseGeometry().vofIdentity()==source.identity()
        &&sameClock(accepted.timeLayer(),source.timeLayer())&&accepted.values().normalVelocity()==zero.normalVelocity()
        &&accepted.values().pressure()==zero.pressure(),"VOF provider attempt changed accepted source or native state");
        for(std::size_t c=0;c<oldCells.size();++c)check(source.plicCells()[c].inventory==oldCells[c].inventory,"VOF provider attempt changed independent phase stocks");};
    const auto samePayload=[&](const I::GeometricVofTransport& a,const I::GeometricVofTransport& b){
        check(a.faces().size()==b.faces().size(),"VOF provider changed actual native face count");
        for(std::size_t c=0;c<oldCells.size();++c)check(a.midpoint().plicCells()[c].inventory==b.midpoint().plicCells()[c].inventory
            &&a.endpoint().plicCells()[c].inventory==b.endpoint().plicCells()[c].inventory,"VOF provider changed the default half/full stocks");
        const auto flux=[&](const I::GeometricVofIntegratedFlux& x,const I::GeometricVofIntegratedFlux& y){
            check(x.total==y.total&&x.positiveTotal==y.positiveTotal&&x.negativeTotal==y.negativeTotal&&x.rawLiquid==y.rawLiquid
                &&x.liquid==y.liquid&&x.gas==y.gas&&x.mass==y.mass&&x.positiveLiquid==y.positiveLiquid&&x.negativeLiquid==y.negativeLiquid
                &&x.positiveGas==y.positiveGas&&x.negativeGas==y.negativeGas&&x.limiter==y.limiter,"VOF provider changed default shared phase/mass flux");};
        for(std::size_t f=0;f<a.faces().size();++f){const auto& x=a.faces()[f];const auto& y=b.faces()[f];
            check(x.nativeRow==y.nativeRow&&x.negativeCell==y.negativeCell&&x.positiveCell==y.positiveCell,"VOF provider changed shared face incidence");flux(x.half,y.half);flux(x.full,y.full);}};
    HydroVofAcConfig config;config.liquidDensity=10;config.gasDensity=1;
    config.momentumDensityRepresentation=HydroMomentumMass::DensityRepresentation::CellAverageProjection;
    check(!config.transportProvider,"default VOF transport provider must be empty");
    const auto baseline=HydroVofStep::prepare(*mesh,accepted,end,config);unchangedAccepted();
    const auto direct=I::GeometricVofTransport::prepareLocalPlicSweep(*mesh,source,baseline.transport().advector(),10,1,config.transport);
    samePayload(baseline.transport(),direct);
    auto configured=config;std::vector<I::GeometricVofTransport> receipts;
    configured.transportProvider=[&](const QuadtreeMesh& live,const I::GeometricVofPublication& old,
        const I::GeometricVofAdvectorSlab& slab,long double rhoL,long double rhoG,const I::GeometricVofTransportConfig& controls){
        check(&live==mesh.get()&&old.identity()==source.identity()&&sameClock(old.timeLayer(),source.timeLayer()),"VOF provider did not receive the actual accepted mesh/source");
        check(sameClock(slab.endTime(),end)&&rhoL==10&&rhoG==1&& &controls==&configured.transport,"VOF provider lost actual clock/densities/controls");
        slab.validateFor(live,old);auto result=I::GeometricVofTransport::prepareLocalPlicSweep(live,old,slab,rhoL,rhoG,controls);
        receipts.push_back(result);return result;};
    const auto provided=HydroVofStep::prepare(*mesh,accepted,end,configured);unchangedAccepted();
    check(receipts.size()==2&&receipts.size()==provided.report().stages.size()
        &&provided.transport().sameInstance(receipts.back())&&!receipts.front().advector().sameInstance(receipts.back().advector()),
        "VOF provider was skipped, reused, or lost the actual per-stage returned receipt");
    check(provided.candidateState().values().normalVelocity()==baseline.candidateState().values().normalVelocity()
        &&provided.candidateState().values().pressure()==baseline.candidateState().values().pressure(),"equivalent VOF provider changed native solve output");
    samePayload(provided.transport(),baseline.transport());
    const auto foreign=I::GeometricVofPublication::prepareLocalPlic(*mesh,source.domain(),source.timeLayer(),alpha);
    for(unsigned kind=0;kind<4;++kind){auto invalid=config;unsigned calls=0;
        invalid.transportProvider=[&,kind](const QuadtreeMesh& live,const I::GeometricVofPublication& old,
            const I::GeometricVofAdvectorSlab& slab,long double rhoL,long double rhoG,const I::GeometricVofTransportConfig& controls){
            ++calls;if(kind==3)return I::GeometricVofTransport::prepareLocalPlicSweep(live,old,slab,rhoL+1,rhoG,controls);
            const auto& actualSource=kind==0?foreign:old;auto next=slab.endTime();if(kind==1)next.end*=2;
            const auto other=I::GeometricVofAdvectorSlab::prepare(live,actualSource,next,slab.beginNormal(),slab.endNormal());
            return I::GeometricVofTransport::prepareLocalPlicSweep(live,actualSource,other,rhoL,rhoG,controls);};
        bool rejected=false;try{HydroVofStep::prepare(*mesh,accepted,end,invalid);}catch(const std::invalid_argument& error){
            const std::string reason=error.what(),expected=kind==0?"foreign VOF owner":kind==3?"different material densities":"different advector slab instance";
            check(reason.find(expected)!=std::string::npos,"VOF provider negative control failed for an unrelated reason");rejected=true;}
        check(rejected&&calls==1,"foreign source, wrong clock/slab or densities were accepted from a VOF provider");unchangedAccepted();}
    HydroVofPrescribedFlow flow;using Mode=I::GeometricVofAdvectorSlab::Boundary::Mode;
    flow.rulesByNativeRow.assign(source.domain().nativeFaceCount(),Mode::Interior);flow.normalVelocity.resize(layout.velocityCount());
    std::vector<double> openVelocity(layout.velocityCount());
    for(int u=0;u<layout.velocityCount();++u){const auto face=layout.domain().face(*mesh,layout.nativeFace(u));
        if(face.horizontal)openVelocity[u]=.02;
        if(face.outwardSign==0)continue;
        flow.rulesByNativeRow[layout.nativeFace(u)]=Mode::Closed;
        if(face.horizontal){flow.rulesByNativeRow[layout.nativeFace(u)]=face.negativeCell<0?Mode::PureLiquidReservoir:Mode::PureGasReservoir;
            flow.prescribedVelocityDofs.push_back(u);flow.normalVelocity[u]=.02;}}
    flow.endpointVelocity=[](int,LAMP::TimeLayer){return Vector2D{0,.02};};
    const HydroNativeValues openValues(layout,openVelocity,std::vector<double>(layout.pressureCount()));
    const auto opened=HydroTwoPhaseState::createVofPrescribedFlow(*mesh,phase,openValues,std::move(flow));
    unsigned openCalls=0;auto openConfig=config;openConfig.transportProvider=[&](const QuadtreeMesh&,
        const I::GeometricVofPublication&,const I::GeometricVofAdvectorSlab&,long double,long double,const I::GeometricVofTransportConfig&)->I::GeometricVofTransport{
        ++openCalls;throw std::runtime_error("open VOF transport provider must not run");};
    rejects([&]{HydroVofStep::preparePrescribedFlow(*mesh,opened,end,openConfig);},"explicit transport provider supports closed VOF flow only");
    rejects([&]{HydroVofStep::stabilityRatesPrescribedFlow(*mesh,opened,openConfig);},"explicit transport provider supports closed VOF flow only");
    check(openCalls==0&&opened.values().normalVelocity()==openVelocity&&opened.phaseGeometry().vofIdentity()==source.identity(),"open VOF provider was invoked or changed accepted state");
    unchangedAccepted();std::cout<<"VOF_TRANSPORT_PROVIDER_ONLY,PASS actual_stages="<<receipts.size()<<" negative_owners_clocks_densities=4 open_rejections=2 default_payload_exact=1\n";
}

}

namespace VofPressureDrivenContract {
using namespace LAMP::Hydrodynamics;
namespace I=LAMP::Interface;
using Clock=LAMP::TimeLayer;
using Mode=I::GeometricVofAdvectorSlab::Boundary::Mode;
using Model=HydroVofPressureOutletMomentumModel;
unsigned scenes=0,steps=0,guards=0,excludedChecks=0,wallChecks=0,freeChecks=0,noSlipOwnerChecks=0,noSlipReuseChecks=0,noSlipZeroChecks=0;
long double peakMomentum=0,peakEnergy=0,peakVolume=0,peakPressureAdjoint=0;
template<class Action>void reject(Action&& action,const char* expected) {
    bool refused=false;try{action();}catch(const std::invalid_argument& error){
        check(std::string(error.what()).find(expected)!=std::string::npos,
            std::string("pressure-driven wrong rejection: ")+error.what());refused=true;}
    check(refused,std::string("pressure-driven missing rejection: ")+expected);++guards;
}
// The actual fluid is a fitted channel cut from a larger native mesh. The
// omitted right-hand cells and their faces must never become fluid unknowns.
long double scene(bool rz,double pressure,Model model,bool negativeControls,bool noSlip=false,bool horizontal=false,bool wallComparison=false) {
    DictionaryMesh spec;spec.Nx=6;spec.Ny=8;spec.minCoord={0,0};spec.maxCoord={1.5,2};
    spec.maxCellNum=4096;spec.maxInitialIterations=0;
    if(horizontal){spec.Nx=8;spec.Ny=6;spec.minCoord={.5,0};spec.maxCoord={2.5,1.5};}
    auto mesh=std::make_unique<QuadtreeMesh>(spec);FaceSegmentGraph graph;
    graph.build(*mesh,[horizontal](const QuadtreeMesh&,const Cell& cell){return horizontal?cell.centerCoord().y<1:cell.centerCoord().x<1;});
    const auto native=NativeLeafDomain::build(*mesh,graph,rz);
    const auto layout=HydroDofLayout::build(*mesh,native);
    const auto divergence=HydroIntegralDivergence::build(*mesh,layout);
    const auto metrics=HydroPhysicalDivergence::build(*mesh,divergence);
    check(layout.pressureCount()<mesh->cellIndexNum()&&layout.velocityCount()<native.nativeFaceCount(),
        "pressure-driven channel must retain the real selected-domain scatter");
    std::vector<double> alpha(mesh->cellIndexNum());
    for(int c:native.cellIndexByRow()){
        const auto& cell=mesh->cells(mesh->cellIndexToLabel[c]);
        const auto lo=mesh->verts(cell.cellToVert[0]).coord(),hi=mesh->verts(cell.cellToVert[2]).coord();
        if(horizontal){const double front=std::clamp(1.5625,lo.x,hi.x);
            alpha[c]=rz?(front*front-lo.x*lo.x)/(hi.x*hi.x-lo.x*lo.x):(front-lo.x)/(hi.x-lo.x);}
        else alpha[c]=std::clamp((1.0625-lo.y)/(hi.y-lo.y),0.,1.);
    }
    const auto initial=I::GeometricVofPublication::prepareLocalPlic(*mesh,native,{},alpha);
    const auto phase=HydroPhaseGeometry::build(*mesh,metrics,initial);
    const HydroNativeValues zero(layout,std::vector<double>(layout.velocityCount()),std::vector<double>(layout.pressureCount()));
    HydroVofPrescribedFlow flow;check(!flow.pressureDrivenReservoirs,"pressure-driven capability changed its default");
    flow.pressureDrivenReservoirs=true;flow.pressureOutletBackflowModel=model;
    flow.normalVelocity.assign(layout.velocityCount(),0);flow.rulesByNativeRow.assign(native.nativeFaceCount(),Mode::Closed);
    std::vector<bool> port(layout.velocityCount());unsigned bottomPorts=0,topPorts=0,solidFaces=0;
    const auto lowPort=[&](const auto& face){return horizontal?(!face.horizontal&&face.center.x==.5):(face.horizontal&&face.center.y==0);};
    const auto highPort=[&](const auto& face){return horizontal?(!face.horizontal&&face.center.x==2.5):(face.horizontal&&face.center.y==2);};
    for(const auto f:native.faces(*mesh)){
        const int u=layout.velocityDof(f.nativeRow);
        if(f.negativeCell<0&&f.positiveCell<0){check(u<0,"pressure-driven solid face has a Hydro unknown");++solidFaces;continue;}
        check(u>=0&&layout.nativeFace(u)==f.nativeRow,"pressure-driven native scatter is not invertible");
        if(f.negativeCell>=0&&f.positiveCell>=0)flow.rulesByNativeRow[f.nativeRow]=Mode::Interior;
        else if(lowPort(f)||highPort(f)){
            const bool bottom=lowPort(f);flow.rulesByNativeRow[f.nativeRow]=bottom?Mode::PureLiquidReservoir:Mode::PureGasReservoir;
            flow.pressureOutlets.push_back({u,bottom?pressure:0});port[u]=true;
            if(bottom)++bottomPorts;else ++topPorts;
        }
    }
    check(solidFaces>0&&bottomPorts>0&&topPorts>0&&flow.prescribedVelocityDofs.empty()&&!flow.endpointVelocity,
        "pressure-driven channel lacks actual free reservoirs or has a fixed feed");
    std::vector<int> expectedNoSlip;int interior=-1,axis=-1;
    for(int u=0;u<layout.velocityCount();++u){const auto f=native.face(*mesh,layout.nativeFace(u));
        // Independent physical predicates, rather than the generated F1 rules.
        const bool wall=horizontal?(f.horizontal&&(f.center.y==0||f.center.y==1))
            :(!f.horizontal&&(f.center.x==0||f.center.x==1));
        if(wall&&stepFaceMeasure(*mesh,layout,u)>0)expectedNoSlip.push_back(u);
        if(f.outwardSign==0)interior=u;if(stepFaceMeasure(*mesh,layout,u)==0)axis=u;
    }
    check(!expectedNoSlip.empty(),"no-slip pressure fixture has no physical wall support");
    if(noSlip)flow.stationaryNoSlipVelocityDofs=expectedNoSlip;
    else check(flow.stationaryNoSlipVelocityDofs.empty(),"stationary no-slip changed the legacy default");
    for(int c=0;c<mesh->cellIndexNum();++c){const auto center=mesh->cells(mesh->cellIndexToLabel[c]).centerCoord();
        check((native.cellRowByIndex()[c]>=0)==(horizontal?center.y<1:center.x<1),"no-slip selected-domain membership differs from physical fixture");}
    unsigned ambientQueries=0;
    if(model==Model::Ambient)flow.pressureOutletBackflowMomentum=[&](int row,Clock clock){
        const int u=layout.velocityDof(row);check(u>=0&&port[u]&&clock.position==Clock::Position::Midpoint
            &&clock.step>=1&&clock.step<=3,"pressure-driven ambient query lost the actual free-port midpoint");
        ++ambientQueries;return Vector2D{};};
    for(double u:zero.normalVelocity())check(u==0&&!std::signbit(u),"pressure-driven initial U must be exact positive zero");
    if(negativeControls){
        const auto invalid=[&](HydroVofPrescribedFlow bad,const char* why){
            reject([&]{HydroTwoPhaseState::createVofPrescribedFlow(*mesh,phase,zero,std::move(bad));},why);};
        auto bad=flow;bad.pressureDrivenReservoirs=false;
        invalid(bad,"VOF prescribed flow requires real exterior velocity data");
        bad=flow;bad.pressureOutletBackflowModel.reset();bad.pressureOutletBackflowMomentum={};
        invalid(bad,"pressure-driven VOF reservoirs require an explicit backflow momentum model");
        bad=flow;bad.prescribedVelocityDofs.push_back(flow.pressureOutlets.front().velocityDof);
        invalid(bad,"pressure-driven VOF reservoirs cannot prescribe endpoint velocity");
        unsigned endpointQueries=0;bad=flow;bad.endpointVelocity=[&](int,Clock){++endpointQueries;return Vector2D{};};
        invalid(bad,"pressure-driven VOF reservoirs cannot prescribe endpoint velocity");
        check(endpointQueries==0,"pressure-driven endpoint rejection sampled a forbidden fixed feed");
        bad=flow;bad.normalVelocity[flow.pressureOutlets.front().velocityDof]=.03125;
        invalid(bad,"free VOF pressure outlet has prescribed normal payload");
        bad=flow;bad.normalVelocity.pop_back();invalid(bad,"different VOF prescribed normal dimensions");
        bad=flow;bad.pressureOutlets.resize(1);
        invalid(bad,"pressure-driven VOF reservoirs require multiple free pressure ports");
        bad=flow;for(auto& p:bad.pressureOutlets)p.pressure=pressure;
        invalid(bad,"pressure-driven VOF reservoirs require distinct physical pressures");
        if(model==Model::InteriorContinuationZeroShear){bad=flow;bad.pressureOutletBackflowMomentum=[](int,Clock){return Vector2D{};};
            invalid(bad,"VOF interior continuation cannot consume ambient callback data");}
        else{bad=flow;bad.pressureOutletBackflowMomentum={};invalid(bad,"VOF ambient backflow momentum requires an explicit query");}
        bad=flow;bad.pressureOutletBackflowModel=static_cast<Model>(99);
        invalid(bad,"unknown VOF pressure-outlet momentum model");
    }
    if(noSlip){
        const auto invalid=[&](std::vector<int> rows,const char* why){auto bad=flow;bad.stationaryNoSlipVelocityDofs=std::move(rows);
            reject([&]{HydroTwoPhaseState::createVofPrescribedFlow(*mesh,phase,zero,bad);},why);};
        invalid({expectedNoSlip[0],expectedNoSlip[0]},"wall rows must be sorted and unique");
        auto unsorted=expectedNoSlip;std::reverse(unsorted.begin(),unsorted.end());invalid(unsorted,"wall rows must be sorted and unique");
        invalid({-1},"positive-area Closed wall");invalid({layout.velocityCount()},"positive-area Closed wall");
        invalid({flow.pressureOutlets.front().velocityDof},"positive-area Closed wall");invalid({interior},"positive-area Closed wall");
        if(axis>=0)invalid({axis},"positive-area Closed wall");
        const auto freeMaterial=HydroPhaseViscosity::buildVofCellAverage(*mesh,phase,.3,.15);
        const auto wallMaterial=HydroPhaseViscosity::buildVofCellAverage(*mesh,phase,.3,.15,expectedNoSlip);
        reject([&]{HydroPhaseViscosity::buildVofCellAverage(*mesh,phase,.3,.15,expectedNoSlip,freeMaterial);},"different stationary no-slip boundary");
        reject([&]{HydroPhaseViscosity::buildVofCellAverage(*mesh,phase,.3,.15,wallMaterial);},"different stationary no-slip boundary");
        auto changed=expectedNoSlip;changed.pop_back();
        reject([&]{HydroPhaseViscosity::buildVofCellAverage(*mesh,phase,.3,.15,changed,wallMaterial);},"different stationary no-slip boundary");
        const auto reused=HydroPhaseViscosity::buildVofCellAverage(*mesh,phase,.6,.3,expectedNoSlip,wallMaterial);
        check(&reused.viscosity().strain().matrix()==&wallMaterial.viscosity().strain().matrix()
            &&reused.viscosity().viscosity()!=wallMaterial.viscosity().viscosity(),"same-wall reuse copied field coefficients or rebuilt the raw owner");
        ++noSlipReuseChecks;
    }
    auto accepted=HydroTwoPhaseState::createVofPrescribedFlow(*mesh,phase,zero,flow);
    HydroVofAcConfig config;config.liquidDensity=2;config.gasDensity=1;
    config.liquidViscosity=wallComparison?.3:.01;config.gasViscosity=wallComparison?.15:.005;config.surfaceTension=0;
    check(HydroVofStep::stabilityRatesPrescribedFlow(*mesh,accepted,config).flowRate==0,
        "pressure-driven resting state was seeded with artificial transport");
    constexpr double dt=1./16384;long double supplied=0,solvedSupplyRate=0;unsigned totalStages=0;
    if(noSlip){auto exactZero=config;exactZero.liquidViscosity=0;exactZero.gasViscosity=0;
        const auto prepared=HydroVofStep::preparePrescribedFlow(*mesh,accepted,Clock{1,0,dt,Clock::Position::End},exactZero);
        const auto& inactive=prepared.frozenMidpoint().viscosity();
        check(inactive.materialModel()==HydroPhaseViscosity::MaterialModel::ExactZero&&inactive.viscosity().matrix().nonZeros()==0
            &&inactive.viscosity().strain().points().empty()&&prepared.report().stages.back().energy.viscousWork==0,
            "no-slip mask inserted strain/penalty/viscosity into the exact zero path");++noSlipZeroChecks;
    }
    const auto initialStocks=initial.plicCells();
    for(int k=1;k<=3;++k){
        const auto before=accepted;const auto& source=before.phaseGeometry().vofPublication();const auto stocks=source.plicCells();
        const auto oldU=before.values().normalVelocity(),oldP=before.values().pressure();
        const Clock end{k,(k-1)*dt,k*dt,Clock::Position::End};
        const auto step=HydroVofStep::preparePrescribedFlow(*mesh,accepted,end,config);
        const auto& frozen=step.frozenMidpoint();const auto& dual=frozen.material();const auto& volume=step.transport();
        if(wallComparison&&noSlip&&!horizontal&&k<=2){
            FrozenPrepareContract::verify(*mesh,frozen,config.saddle,k==1?"pressure_startup":"pressure_retained_history");
            ++FrozenPrepareContract::pressureReferences;
        }
        const auto& mid=step.solvedMidpointValues();const auto& endpoint=frozen.solveResult().values();
        const auto& report=step.report().stages.back();totalStages+=static_cast<unsigned>(step.report().stages.size());
        check(step.report().linearSolves>=1&&step.report().linearSolves<=2&&step.report().usedVelocityHistory==(k>1)
            &&dual.usesConservativeVofMomentum()&&dual.material().vofTransport().sameInstance(volume)
            &&volume.sourceIdentity()==source.identity()&&sameClock(volume.endTime(),end),
            "pressure-driven step lost original AC/history/unique transport source");
        check(frozen.viscosityMatrix().nonZeros()>0&&!step.hasCapillaryLoad()&&!frozen.hasMaxwellLoad()
            &&frozen.convection().isVofOpenFlow()&&frozen.saddle().pressureGaugeCount()==0,
            "pressure-driven step lost positive viscosity, explicit open momentum or physical pressure reference");
        for(const auto& stage:step.report().stages){
            check(stage.saddle.velocitySpdTier==static_cast<int>(HydroVelocitySpdTier::None)
                &&!stage.saddle.velocitySpdGramAttempted&&!stage.saddle.velocitySpdGramBoundsAvailable
                &&stage.saddle.velocitySpdCertificateSeconds==0&&stage.saddle.velocitySpdLltSeconds==0,
                "pressure-driven production stage falsely reported a requested SPD audit");
            check(stage.consistencyPower==0&&stage.momentumPowerRelativeDefect<=config.momentumPowerTolerance
                &&std::isfinite(stage.energy.closureRelativeResidual)&&stage.energy.closureRelativeResidual<=1e-11L,
                "pressure-driven step violated original conservative momentum/power/energy gate");
        }
        const auto& productionVelocity=frozen.saddle().velocityBlock();const auto& unrequestedSpd=productionVelocity.spdReport();
        check(unrequestedSpd.tier==HydroVelocitySpdTier::None&&!unrequestedSpd.gramAttempted&&!unrequestedSpd.gramBoundsAvailable
            &&unrequestedSpd.certificateSeconds==0&&unrequestedSpd.lltSeconds==0,
            "actual Frozen velocity falsely reported a global SPD certificate");
        if(k==1&&rz&&negativeControls&&model==Model::InteriorContinuationZeroShear){
            // One retained public audit/solve of the SAME actual producer inputs.
            // Omitted diagnostic FP history is deliberately not compared.
            const auto referenceVelocity=HydroSaddleVelocityBlock::buildWithEnergyTransportSplit(*mesh,layout,
                frozen.velocityMatrix(),dual.endMass(),frozen.viscosity().viscosity(),dt,dual.endMass().geometry());
            check(referenceVelocity.spdReport().tier!=HydroVelocitySpdTier::None,
                "public arbitrary-matrix factory lost its original SPD certification");
            const auto& actualMatrix=productionVelocity.matrix();const auto& referenceMatrix=referenceVelocity.matrix();
            const Matrix at=actualMatrix.transpose(),symmetric=.5*(actualMatrix+at),skew=.5*(actualMatrix-at),legacy=symmetric+skew;
            FrozenPrepareContract::matrix(referenceMatrix,legacy,"public SPD reference lost its original S/W roundtrip");
            FrozenPrepareContract::verify(*mesh,frozen,config.saddle,"pressure_public_audit_raw_equations");
            const auto referenceSaddle=HydroSinglePhaseSaddle::buildWithPressureOutlets(frozen.pressureLoad(),referenceVelocity,
                frozen.boundary().pressureOutlets);
            const auto referenceSolved=referenceSaddle.solve(*mesh,HydroSaddleRhs(layout,frozen.integralMomentumRhs(),frozen.boundary().normalVelocity));
            // The old public roundtrip solves its own rounded equation. Its original
            // acceptance gates remain mandatory; cross-method byte equality is not claimed.
            const auto& solvedReport=referenceSolved.report();
            check(solvedReport.algebraicRelativeResidual<=1e-11&&solvedReport.freeMomentumDefect<=1e-11
                &&solvedReport.maxCellContinuityDefect<=1e-11&&solvedReport.maxPressureMeanDefect<=1e-11
                &&solvedReport.maxBoundaryFluxDefect<=1e-12,
                "retained public SPD reference solve failed the original actual residual gates");
            std::cout<<"PRODUCTION_SPD_UNREQUESTED_REFERENCE_PASS,RZ,actual_matrix_words,"<<actualMatrix.nonZeros()
                <<",public_audit_tier,"<<static_cast<int>(referenceVelocity.spdReport().tier)<<",independent_raw_equations,1,legacy_public_gates,1,cross_method_byte_identity,0\n";
        }
        constexpr long double continuityTolerance=256*std::numeric_limits<double>::epsilon();
        check(dual.report().maxRelativeContinuityDefect<=continuityTolerance
            &&dual.report().maxRelativeStoredMassContinuityDefect<=continuityTolerance,
            "pressure-driven step weakened original actual endpoint mass continuity");
        const auto applied=frozen.convection().applyVofOpen(*mesh,mid,dual.sources(),dual.timeLayer());
        const auto cp=frozen.pressureLoad().apply(*mesh,endpoint),div=divergence.apply(*mesh,mid);
        const Eigen::Map<const Eigen::VectorXd> vm(mid.normalVelocity().data(),layout.velocityCount());
        const Eigen::VectorXd viscous=frozen.viscosityMatrix()*vm;
        std::vector<bool> free(layout.velocityCount());for(int u:frozen.saddle().freeVelocityDofs())free[u]=true;
        long double momentumError=0,momentumScale=0,kinetic=0,massCross=0,pressureWork=0,tractionWork=0,viscousWork=0,convectionWork=0,energyScale=0;
        long double adjoint=0,adjointScale=0,bottomRate=0,topRate=0;
        for(int u=0;u<layout.velocityCount();++u){
            const auto f=native.face(*mesh,layout.nativeFace(u));const long double a=oldU[u],b=endpoint.normalVelocity()[u],v=mid.normalVelocity()[u];
            const long double m0=dual.beginMass().masses()[u],m1=dual.endMass().masses()[u];
            const long double traction=port[u]?-static_cast<long double>(f.outwardSign)*stepFaceMeasure(*mesh,layout,u)*(lowPort(f)?pressure:0):0;
            scaled(frozen.saddle().pressureOutletLoads()[u],traction,std::abs(traction),"pressure-driven original Saddle traction",1e-12);
            if(port[u]){check(free[u]&&std::find(frozen.convection().fixedVelocityDofs().begin(),frozen.convection().fixedVelocityDofs().end(),u)
                ==frozen.convection().fixedVelocityDofs().end(),"pressure-driven reservoir normal was fixed");++freeChecks;
                check(flow.normalVelocity[u]==0,"pressure-driven free reservoir acquired fixed normal payload");
                if(lowPort(f))bottomRate+=stepFaceMeasure(*mesh,layout,u)*v;else topRate+=stepFaceMeasure(*mesh,layout,u)*v;}
            if(f.outwardSign!=0&&!port[u]){check(a==0&&b==0&&v==0,"pressure-driven fitted wall leaked");++wallChecks;}
            check(applied.momentumLoads()[u]==-applied.momentumRates()[u],"pressure-driven conservative momentum acquired uR");
            const long double rate=(m1*b-m0*a)/dt,load=applied.momentumLoads()[u];
            const auto residual=rate+viscous[u]+cp[u]-load-traction;
            const auto rowScale=std::abs(rate)+std::abs(static_cast<long double>(viscous[u]))+std::abs(static_cast<long double>(cp[u]))+std::abs(load)+std::abs(traction);
            if(free[u]){momentumError=std::hypot(momentumError,residual);momentumScale=std::hypot(momentumScale,rowScale);}
            const auto ke=.5L*(m1*b*b-m0*a*a),mc=.5L*(m1-m0)*a*b,pw=dt*v*cp[u],tw=dt*v*traction,vw=dt*v*viscous[u],cw=dt*v*load;
            kinetic+=ke;massCross+=mc;pressureWork+=pw;tractionWork+=tw;viscousWork+=vw;convectionWork+=cw;
            energyScale+=std::abs(ke)+std::abs(mc)+std::abs(pw)+std::abs(tw)+std::abs(vw)+std::abs(cw)+dt*std::abs(v)*rowScale;
            adjoint+=v*cp[u];adjointScale+=std::abs(v*cp[u]);
        }
        for(int p=0;p<layout.pressureCount();++p){const long double term=endpoint.pressure()[p]*static_cast<long double>(div[p]);adjoint+=term;adjointScale+=std::abs(term);}
        const auto& factor=frozen.viscosity().viscosity();const auto strain=factor.strain().apply(*mesh,mid);
        const double retainedDissipation=factor.dissipation(*mesh,mid);
        if(noSlip){
            check(accepted.vofPrescribedFlow()->stationaryNoSlipVelocityDofs==expectedNoSlip
                &&step.candidateState().vofPrescribedFlow()->stationaryNoSlipVelocityDofs==expectedNoSlip
                &&factor.strain().stationaryNoSlipVelocityDofs()==expectedNoSlip&&!factor.strain().isStressCompatible()
                &&frozen.viscosity().geometry().vofIdentity()==volume.midpoint().identity()
                &&sameClock(frozen.viscosity().geometry().timeLayer(),volume.midpoint().timeLayer())
                && &factor.strain().metrics().divergence().matrix()==&metrics.divergence().matrix(),
                "actual no-slip consumer lost the exact wall/phase/clock/divergence owner");
            check(std::isfinite(retainedDissipation)&&retainedDissipation>0,"actual no-slip pressure flow has no positive physical dissipation");
            for(std::size_t q=0;q<factor.strain().points().size();++q){const auto& point=factor.strain().points()[q];
                const int c=layout.nativeCell(point.cellDof);const auto& inventory=volume.midpoint().plicCells()[c].inventory;
                const auto capacity=volume.midpoint().cells()[c].capacity;
                const double expected=static_cast<double>((inventory[0]/capacity)*config.liquidViscosity+(inventory[1]/capacity)*config.gasViscosity);
                check(factor.viscosity()[q]==expected,"no-slip stage reused previous viscosity coefficients instead of current independent stocks");}
            ++noSlipOwnerChecks;
        }else check(factor.strain().stationaryNoSlipVelocityDofs().empty(),"empty-mask branch acquired no-slip strain");
        long double matrixGross=0,matrixSigned=0,asymmetryWork=0,gramWork=0,gramGross=0,maxAbsStrain=0;
        const auto& matrix=frozen.viscosityMatrix();
        for(int row=0;row<matrix.outerSize();++row)for(Matrix::InnerIterator entry(matrix,row);entry;++entry){
            const long double term=static_cast<long double>(dt)*vm[entry.row()]*entry.value()*vm[entry.col()];
            matrixSigned+=term;matrixGross+=std::abs(term);
            if(entry.row()<entry.col())asymmetryWork+=std::abs(static_cast<long double>(dt)*vm[entry.row()]
                *(static_cast<long double>(entry.value())-matrix.coeff(entry.col(),entry.row()))*vm[entry.col()]);
        }
        const auto& e=factor.strain().matrix();
        for(int row=0;row<e.outerSize();++row){
            long double absolute=0;for(Matrix::InnerIterator entry(e,row);entry;++entry)absolute+=std::abs(static_cast<long double>(entry.value())*vm[entry.col()]);
            const long double weight=static_cast<long double>(dt)*factor.rowWeights()[row],value=strain[row];
            gramWork+=weight*value*value;gramGross+=weight*absolute*absolute;maxAbsStrain=std::max(maxAbsStrain,std::abs(value));
        }
        // The original public source-qualified certificate bounds each row of
        // diag(s)*(L(K)-E^T W E)*diag(s), with 0<s_i<=max|U|. Summing n
        // rows bounds its quadratic work; retain full-K asymmetry separately.
        const auto gramAssemblyBound=factor.coarseScaledGramLowerErrorBound(vm.cwiseAbs().maxCoeff());
        check(gramAssemblyBound.has_value(),"pressure-driven original stored Gram error certificate unavailable");
        // A single conservative gamma counts all stored sparse products and
        // reductions in K*U, E*U, and the two work sums. Both gross scales use
        // these original stored operators; no physical norm or absolute floor.
        const long double operations=2.L*(matrix.nonZeros()+e.nonZeros()+matrix.rows()+e.rows())+16;
        const long double epsilon=std::numeric_limits<double>::epsilon(),neps=operations*epsilon;
        check(neps<.5L,"pressure-driven viscosity roundoff count is not bounded");
        const long double viscousRoundoffBound=static_cast<long double>(dt)*layout.velocityCount()*(*gramAssemblyBound)
            +asymmetryWork+neps/(1-neps)*(matrixGross+gramGross+std::abs(gramWork));
        const long double nativeGramGap=std::abs(viscousWork-static_cast<long double>(dt)*retainedDissipation);
        const long double storedGramGap=std::abs(matrixSigned-gramWork);
        const bool viscosityCertified=std::isfinite(retainedDissipation)&&retainedDissipation>=0
            &&std::isfinite(viscousRoundoffBound)&&nativeGramGap<=viscousRoundoffBound&&storedGramGap<=viscousRoundoffBound;
        if(!(bottomRate>0&&topRate>0&&tractionWork>0&&viscosityCertified)){
            std::ostringstream diagnostic;diagnostic<<std::setprecision(std::numeric_limits<long double>::max_digits10)
                <<"pressure-driven actual pressure failed to supply solved liquid flow/work scene="<<(rz?"RZ":"XY")
                <<" model="<<(model==Model::Ambient?"ambient":"interiorContinuationZeroShear")<<" pressure="<<pressure<<" k="<<k
                <<" bottomRate="<<bottomRate<<" topRate="<<topRate<<" tractionWork="<<tractionWork<<" viscousWork="<<viscousWork
                <<" matrixSignedWork="<<matrixSigned<<" matrixProductGross="<<matrixGross
                <<" retainedGramWork="<<gramWork<<" maxAbsRetainedStrain="<<maxAbsStrain<<" retainedDissipation="<<retainedDissipation
                <<" viscousRoundoffBound="<<viscousRoundoffBound<<" nativeGramGap="<<nativeGramGap<<" storedGramGap="<<storedGramGap
                <<" matrixWorkOverEpsGross="<<(epsilon*matrixGross>0?viscousWork/(epsilon*matrixGross):0)
                <<" failedDirection="<<(!(bottomRate>0&&topRate>0))<<" failedPressureWork="<<(!(tractionWork>0))
                <<" failedViscosityCertificate="<<(!viscosityCertified);
            check(false,diagnostic.str());
        }
        scaled(bottomRate,topRate,std::abs(bottomRate)+std::abs(topRate),"pressure-driven solved physical volume flux",1e-12);solvedSupplyRate=bottomRate;
        peakPressureAdjoint=std::max(peakPressureAdjoint,scaled(adjoint,0,adjointScale,"pressure-driven original C=-DI transpose work",1e-12));
        peakMomentum=std::max(peakMomentum,scaled(momentumError,0,momentumScale,"pressure-driven independent free momentum",1e-11));
        const auto energy=kinetic+massCross+pressureWork+viscousWork-convectionWork-tractionWork;
        peakEnergy=std::max(peakEnergy,scaled(energy,0,energyScale,"pressure-driven independent native KE and pressure work",1e-11));
        scaled(report.pressureOutletWork,tractionWork,energyScale,"pressure-driven boundary pressure work counted once",1e-12);
        check(report.prescribedReactionWork==0,"pressure-driven channel invented fixed feed reaction work");
        for(bool half:{true,false}){
            const auto& now=half?volume.midpoint():volume.endpoint();const auto& ledger=half?volume.report().boundaryHalf:volume.report().boundaryFull;
            std::array<long double,2> oldStock{},newStock{},outward{},gross{};
            long double liquidSupply=0;
            for(int c:native.cellIndexByRow())for(unsigned phaseIndex=0;phaseIndex<2;++phaseIndex){oldStock[phaseIndex]+=stocks[c].inventory[phaseIndex];newStock[phaseIndex]+=now.plicCells()[c].inventory[phaseIndex];}
            for(const auto f:native.faces(*mesh)){
                const auto& q=half?volume.faces()[f.nativeRow].half:volume.faces()[f.nativeRow].full;
                if(f.negativeCell<0&&f.positiveCell<0){check(q.total==0&&q.liquid==0&&q.gas==0,"pressure-driven excluded solid face acquired flux");++excludedChecks;continue;}
                if(f.outwardSign==0)continue;
                outward[0]+=f.outwardSign*q.liquid;outward[1]+=f.outwardSign*q.gas;
                gross[0]+=q.positiveLiquid+q.negativeLiquid;gross[1]+=q.positiveGas+q.negativeGas;
                if(lowPort(f))liquidSupply-=f.outwardSign*q.liquid;
            }
            for(unsigned phaseIndex=0;phaseIndex<2;++phaseIndex){const auto scale=std::abs(oldStock[phaseIndex])+std::abs(newStock[phaseIndex])+gross[phaseIndex];
                peakVolume=std::max(peakVolume,scaled(newStock[phaseIndex]-oldStock[phaseIndex]+outward[phaseIndex],0,scale,"pressure-driven actual independent phase stock balance",1e-12));}
            check(ledger.enabled&&ledger.sourceIdentity==source.identity()&&ledger.physicalRelativeConservation<=1e-12L
                &&ledger.algebraicRelativeConservation<=1e-12L&&liquidSupply>0,"pressure-driven phase/WY/gross supply gate failed");
            if(!half)supplied+=liquidSupply;
            for(std::size_t c=0;c<now.plicCells().size();++c)if(native.cellRowByIndex()[c]<0){
                check(now.plicCells()[c].inventory==std::array<long double,2>{0,0},"pressure-driven omitted solid gained stock");++excludedChecks;}
        }
        check(accepted.sameIdentity(before)&&accepted.values().normalVelocity()==oldU&&accepted.values().pressure()==oldP,
            "pressure-driven preparation changed accepted native state");
        for(std::size_t c=0;c<stocks.size();++c)check(source.plicCells()[c].inventory==stocks[c].inventory,"pressure-driven preparation changed accepted phase stocks");
        step.commit(*mesh,accepted);check(sameClock(accepted.timeLayer(),end)&&accepted.vofPrescribedFlow()->pressureDrivenReservoirs
            &&accepted.vofPrescribedFlow()->prescribedVelocityDofs.empty()&&!accepted.vofPrescribedFlow()->endpointVelocity,
            "pressure-driven commit lost full clock or restored a fixed feed");++steps;
    }
    for(std::size_t c=0;c<initialStocks.size();++c)check(initial.plicCells()[c].inventory==initialStocks[c].inventory,"pressure-driven commits mutated original publication");
    check(model==Model::Ambient?ambientQueries==totalStages*flow.pressureOutlets.size():ambientQueries==0,
        "pressure-driven backflow model did not consume exactly its actual stage data");
    const long double rate=supplied/(3*dt);check(rate>0,"pressure-driven actual integrated supply is zero");++scenes;
    std::cout<<"VOF_PRESSURE_DRIVEN_SCENE,"<<(rz?"RZ":"XY")<<",model,"<<(model==Model::Ambient?"ambient":"interiorContinuationZeroShear")
        <<",wall_model,"<<(noSlip?"stationaryNoSlip":"normalImpermeable")<<",orientation,"<<(horizontal?"horizontal":"vertical")<<",pressure,"<<pressure<<",actual_solved_supply_rate,"<<solvedSupplyRate<<",transport_average_supply_rate,"<<rate<<",steps,3,initial_U,0,fixed_feed,0\n";
    return solvedSupplyRate;
}
void run(){
    std::cout<<std::setprecision(17);
    for(bool rz:{false,true}){
        const auto high=scene(rz,45.2,Model::InteriorContinuationZeroShear,true);
        const auto low=scene(rz,22.6,Model::InteriorContinuationZeroShear,false);
        check(high>low&&high-low>1e-6L*high,"pressure-driven physical pressure change did not change solved supply rate");
        scene(rz,45.2,Model::Ambient,true);
    }
    check(scenes==6&&steps==18&&guards==40&&excludedChecks>0&&wallChecks>0&&freeChecks>0,"pressure-driven selected XY/RZ coverage incomplete");
    std::cout<<"VOF_PRESSURE_DRIVEN_ONLY,PASS scenes="<<scenes<<" accepted_steps="<<steps<<" negative_controls="<<guards
        <<" momentum="<<peakMomentum<<" energy="<<peakEnergy<<" volume="<<peakVolume<<" pressure_adjoint="<<peakPressureAdjoint<<'\n';
}
} // namespace VofPressureDrivenContract
namespace VofNoSlipPressureContract {
void run(){
    using namespace VofPressureDrivenContract;
    std::cout<<std::setprecision(17);
    long double smallestResponse=std::numeric_limits<long double>::infinity();
    for(bool rz:{false,true})for(bool horizontal:{false,true}){
        const auto resisted=scene(rz,45.2,Model::InteriorContinuationZeroShear,false,true,horizontal,true);
        const auto original=scene(rz,45.2,Model::InteriorContinuationZeroShear,false,false,horizontal,true);
        const long double response=(original-resisted)/original;
        check(original>resisted&&resisted>0&&response>1e-6L,
            "actual stationary no-slip viscosity did not reduce solved supply against identical normal-impermeable pressure flow");
        smallestResponse=std::min(smallestResponse,response);
    }
    check(FrozenPrepareContract::pressureReferences==4,"Frozen pressure startup/history XY/RZ reference coverage incomplete");
    check(scenes==8&&steps==24&&noSlipOwnerChecks==12&&noSlipReuseChecks==4&&noSlipZeroChecks==4
        &&VofPressureDrivenContract::guards>=36&&wallChecks>0&&freeChecks>0&&excludedChecks>0,"no-slip pressure selected XY/RZ/orientation coverage incomplete");
    std::cout<<"VOF_NO_SLIP_PRESSURE_ONLY,PASS scenes="<<scenes<<" accepted_steps="<<steps<<" no_slip_steps="<<noSlipOwnerChecks
        <<" guards="<<VofPressureDrivenContract::guards<<" same_boundary_reuse="<<noSlipReuseChecks<<" exact_zero_checks="<<noSlipZeroChecks
        <<" minimum_supply_response="<<smallestResponse<<" momentum="<<peakMomentum<<" energy="<<VofPressureDrivenContract::peakEnergy
        <<" volume="<<peakVolume<<" pressure_adjoint="<<peakPressureAdjoint<<'\n';
}
}
namespace Task38StateValidationContract {
using namespace LAMP::Hydrodynamics;
namespace I=LAMP::Interface;
using Clock=LAMP::TimeLayer;
void require38(bool ok,const char* message){if(!ok)throw std::runtime_error(message);}
template<class Action>std::string rejected38(Action&& action) {
    try{action();}catch(const std::invalid_argument& e){return e.what();}
    throw std::runtime_error("Task38 expected original public invalid_argument");
}
struct Floating38 {
    std::fenv_t native{};int error=errno;
#if defined(_M_X64) || defined(__SSE2__)
    unsigned csr=_mm_getcsr();
#endif
    Floating38(){require38(std::fegetenv(&native)==0,"Task38 FP capture");}
    void restore()const{require38(std::fesetenv(&native)==0,"Task38 FP restore");
#if defined(_M_X64) || defined(__SSE2__)
        _mm_setcsr(csr);
#endif
        errno=error;}
    bool same(const Floating38& b)const{return error==b.error&&std::memcmp(&native,&b.native,sizeof native)==0
#if defined(_M_X64) || defined(__SSE2__)
        &&csr==b.csr
#endif
        ;}
};
struct FloatingRestore38 {Floating38 original;~FloatingRestore38(){original.restore();}};
template<class T>void bits38(const std::vector<T>& a,const std::vector<T>& b,const char* message) {
    require38(a.size()==b.size()&&(a.empty()||std::memcmp(a.data(),b.data(),a.size()*sizeof(T))==0),message);
}
void run() {
    FloatingRestore38 restore;std::size_t scenes=0,mutations=0,guards=0,callbackGuards=0;
    for(bool rz:{false,true}){
        auto mesh=::makeMesh(4,1);const auto domain=::domainFor(*mesh,rz);
        require38(!domain.uniformLeaves(),"Task38 needs actual native 2:1");
        const auto layout=HydroDofLayout::build(*mesh,domain);
        const auto divergence=HydroIntegralDivergence::build(*mesh,layout);
        const auto metrics=HydroPhysicalDivergence::build(*mesh,divergence);
        const Clock initial{0,0,0,Clock::Position::Initial},end{1,0,.03125,Clock::Position::End};
        std::vector<std::array<long double,2>> stocks(mesh->cellIndexNum());
        for(std::size_t i=0;i<stocks.size();++i)stocks[i]={domain.cellVolumes()[domain.cellRowByIndex()[i]],0};
        const auto pub=I::GeometricVofPublication::prepareLocalPlicInventory(*mesh,domain,initial,stocks,{});
        const auto endPub=I::GeometricVofPublication::prepareLocalPlicInventory(*mesh,domain,end,stocks,{});
        const auto phase=HydroPhaseGeometry::build(*mesh,metrics,pub),endPhase=HydroPhaseGeometry::build(*mesh,metrics,endPub);
        const HydroNativeValues zero(layout,std::vector<double>(layout.velocityCount()),std::vector<double>(layout.pressureCount()));
        const auto accepted=HydroTwoPhaseState::createVof(*mesh,phase,zero);
        const auto mass=HydroMomentumMass::buildCellAverageProjection(*mesh,accepted.momentumGeometry(),endPhase,10,1);
        const auto closed=HydroTwoPhaseState::createVofAccepted(*mesh,endPhase,zero,accepted,mass);
        require38(closed.velocityHistory()&&closed.retainedMass(),"Task38 actual End/history/mass missing");
        const auto oldU=closed.values().normalVelocity(),oldP=closed.values().pressure(),oldMass=closed.retainedMass()->masses();
        const auto oldD=divergence.apply(*mesh,closed.values());
        require38(std::feraiseexcept(FE_INEXACT)==0,"Task38 prime actual diagnostic inexact");
        errno=173;const Floating38 before;
        closed.validateFor(*mesh);closed.validateVofFor(*mesh,endPub.identity(),end);
        const Floating38 after;require38(before.same(after),"Task38 successful State validation changed complete FP/errno");
        auto clone=mesh->clone();closed.validateFor(*clone);closed.validateVofFor(*clone,endPub.identity(),end);
        bits38(oldD,divergence.apply(*clone,closed.values()),"Task38 public divergence/clone bits");
        const auto foreign=I::GeometricVofPublication::prepareLocalPlicInventory(*mesh,domain,end,stocks,{});
        const auto foreignError=rejected38([&]{endPub.validateFor(*mesh,foreign.identity(),end);});
        require38(rejected38([&]{closed.validateVofFor(*mesh,foreign.identity(),end);})==foreignError,"Task38 original foreign-owner first error");++guards;
        auto wrong=end;++wrong.step;
        const auto clockError=rejected38([&]{endPub.validateFor(*mesh,endPub.identity(),wrong);});
        require38(rejected38([&]{closed.validateVofFor(*mesh,endPub.identity(),wrong);})==clockError,"Task38 original complete-clock first error");++guards;
        const auto freshLayout=HydroDofLayout::build(*mesh,::domainFor(*mesh,rz));
        const HydroNativeValues foreignValues(freshLayout,oldU,oldP);
        const auto layoutError=rejected38([&]{foreignValues.validateFor(layout,*mesh);});
        require38(rejected38([&]{divergence.apply(*mesh,foreignValues);})==layoutError,"Task38 public DI intrinsic layout first error");++guards;
        int parent=-1;
        for(int i=0;i<mesh->cellIndexNum()&&parent<0;++i){const int p=mesh->cells(mesh->cellIndexToLabel[i]).parent;
            if(p>=0&&std::find(mesh->cellIndexToLabel,mesh->cellIndexToLabel+mesh->cellIndexNum(),p)==mesh->cellIndexToLabel+mesh->cellIndexNum())parent=p;}
        require38(parent>=0,"Task38 nonindexed actual ancestor missing");
        auto& ancestor=mesh->cells(parent);
        const auto badMesh=[&]{
            const auto expected=rejected38([&]{layout.validateMesh(*mesh);});
            require38(rejected38([&]{closed.validateFor(*mesh);})==expected,"Task38 public State stale-mesh first error");++guards;
            require38(rejected38([&]{closed.validateVofFor(*mesh,endPub.identity(),end);})==expected,"Task38 independent State wrapper reused prior proof");++guards;
            require38(rejected38([&]{divergence.apply(*mesh,closed.values());})==expected,"Task38 independent public DI reused prior proof");++guards;
            require38(rejected38([&]{closed.validateVofFor(*mesh,foreign.identity(),end);})==foreignError,"Task38 owner must precede bad mesh");++guards;
            ++mutations;
        };
        auto& leaf=mesh->cells(mesh->cellIndexToLabel[0]);const auto center=leaf.centerCoord();
        leaf.setCenterCoord(center+Vector2D{.0009765625,0});badMesh();leaf.setCenterCoord(center);closed.validateFor(*mesh);
        std::swap(ancestor.son[0],ancestor.son[1]);badMesh();std::swap(ancestor.son[0],ancestor.son[1]);closed.validateFor(*mesh);
        auto& edge=mesh->vEdges(mesh->vEdgeIndexToLabel[0]);
        std::swap(edge.edgeToVert[0],edge.edgeToVert[1]);badMesh();std::swap(edge.edgeToVert[0],edge.edgeToVert[1]);closed.validateFor(*mesh);
        auto& vertex=mesh->verts(mesh->vertIndexToLabel[0]);const auto coord=vertex.coord();
        vertex.setCoord(coord+Vector2D{0,.0009765625});badMesh();vertex.setCoord(coord);closed.validateFor(*mesh);
        bool armed=false,changed=false;
        HydroVofPrescribedFlow flow;using Mode=I::GeometricVofAdvectorSlab::Boundary::Mode;
        flow.rulesByNativeRow.assign(domain.nativeFaceCount(),Mode::Interior);
        flow.normalVelocity.resize(layout.velocityCount());std::vector<double> openU(layout.velocityCount());
        for(int u=0;u<layout.velocityCount();++u){const auto f=domain.face(*mesh,layout.nativeFace(u));
            if(f.horizontal)openU[u]=.03125;
            if(f.outwardSign==0)continue;
            flow.rulesByNativeRow[layout.nativeFace(u)]=Mode::Closed;
            if(f.horizontal){flow.rulesByNativeRow[layout.nativeFace(u)]=Mode::PureLiquidReservoir;
                flow.prescribedVelocityDofs.push_back(u);flow.normalVelocity[u]=.03125;}}
        flow.endpointVelocity=[&](int,Clock){
            if(armed&&!changed){std::swap(ancestor.son[0],ancestor.son[1]);changed=true;}
            return Vector2D{0,.03125};};
        const HydroNativeValues openValues(layout,openU,std::vector<double>(layout.pressureCount()));
        const auto opened=HydroTwoPhaseState::createVofPrescribedFlow(*mesh,phase,openValues,flow);
        const auto openMass=HydroMomentumMass::buildCellAverageProjection(*mesh,opened.momentumGeometry(),endPhase,10,1);
        const auto openEnd=HydroTwoPhaseState::createVofAccepted(*mesh,endPhase,openValues,opened,openMass);
        require38(openEnd.velocityHistory()&&openEnd.retainedMass(),"Task38 open actual history/endMass absent");
        openEnd.validateFor(*mesh);openEnd.validateVofFor(*mesh,endPub.identity(),end);
        openEnd.validateFor(*clone);openEnd.validateVofFor(*clone,endPub.identity(),end);
        require38(std::feraiseexcept(FE_INEXACT)==0,"Task38 open FP prime");errno=173;const Floating38 openBefore;
        openEnd.validateFor(*mesh);openEnd.validateVofFor(*mesh,endPub.identity(),end);
        const Floating38 openAfter;require38(openBefore.same(openAfter),"Task38 open successful FP/errno changed");
        armed=true;
        const auto callbackError=rejected38([&]{openEnd.validateFor(*mesh);});
        require38(changed,"Task38 actual endpoint callback was not invoked");
        const auto staleError=rejected38([&]{layout.validateMesh(*mesh);});
        require38(callbackError==staleError,"Task38 pre-callback proof leaked across ancestor-only mutation");
        std::swap(ancestor.son[0],ancestor.son[1]);armed=false;changed=false;++callbackGuards;
        openEnd.validateFor(*mesh);openEnd.validateVofFor(*mesh,endPub.identity(),end);
        // The original Boundary FP gate must remain before any later physical
        // validation, without this optimization normalizing caller controls.
        const Floating38 nearest;require38(std::fesetround(FE_UPWARD)==0,"Task38 set real upward");
#if defined(_M_X64) || defined(__SSE2__)
        _mm_setcsr((_mm_getcsr()&~0x6000u)|0x4000u);
#endif
        const auto fpError=rejected38([&]{closed.boundary().validateFor(*mesh,closed.momentumGeometry());});
        require38(rejected38([&]{closed.validateFor(*mesh);})==fpError,"Task38 original Boundary FP first error");
        nearest.restore();++guards;
        bits38(oldU,closed.values().normalVelocity(),"Task38 closed U changed");bits38(oldP,closed.values().pressure(),"Task38 closed p changed");
        bits38(oldMass,closed.retainedMass()->masses(),"Task38 retained M changed");bits38(oldD,divergence.apply(*mesh,closed.values()),"Task38 original D output changed");
        bits38(openU,openEnd.values().normalVelocity(),"Task38 open U changed");
        for(std::size_t i=0;i<stocks.size();++i)require38(pub.plicCells()[i].inventory==stocks[i]&&endPub.plicCells()[i].inventory==stocks[i],"Task38 immutable stocks changed");
        ++scenes;
    }
    require38(scenes==2&&mutations==8&&callbackGuards==2,"Task38 finite public coverage incomplete");
    std::cout<<"HYDRO_STATE_VALIDATION_BATCH_PASS,XY_RZ_true2to1,"<<scenes<<",mutations,"<<mutations<<",guards,"<<guards
        <<",actual_callback_guards,"<<callbackGuards<<",End_history_mass,1,content_clone,1,public_DI,1,no_PDE\n";
}
} // namespace Task38StateValidationContract

namespace VofCheckpointContract {
using namespace LAMP::Hydrodynamics;
namespace I=LAMP::Interface;
using Clock=LAMP::TimeLayer;
unsigned comparisons=0,guards=0,physicalSteps=0;
template<class Action>void refuse(Action&& action,const char* message){
    bool caught=false;try{action();}catch(const std::exception&){caught=true;}
    check(caught,message);++guards;
}
template<class T>void bits(const std::vector<T>& a,const std::vector<T>& b,const char* message){
    check(a.size()==b.size()&&(a.empty()||std::memcmp(a.data(),b.data(),a.size()*sizeof(T))==0),message);
}
void compare(const HydroTwoPhaseState& expected,const HydroTwoPhaseState& actual){
    check(sameClock(expected.timeLayer(),actual.timeLayer()),"checkpoint changed the complete accepted clock");
    bits(expected.values().normalVelocity(),actual.values().normalVelocity(),"checkpoint changed native U bits");
    bits(expected.values().pressure(),actual.values().pressure(),"checkpoint changed physical p bits");
    check(expected.phaseGeometry().vofPublication().checkpoint()==actual.phaseGeometry().vofPublication().checkpoint(),
        "checkpoint changed independent stock, exact capacity certificate or actual PLIC geometry");
    const auto* a=expected.velocityHistory();const auto* b=actual.velocityHistory();
    check(bool(a)==bool(b),"checkpoint changed the actual velocity-history presence");
    if(a){check(sameClock(a->clock,b->clock)&&a->model==b->model,"checkpoint changed predecessor clock/model");
        bits(a->values.normalVelocity(),b->values.normalVelocity(),"checkpoint changed predecessor U bits");
        bits(a->values.pressure(),b->values.pressure(),"checkpoint changed predecessor p bits");}
    const auto* m=expected.retainedMass();const auto* n=actual.retainedMass();
    check(bool(m)==bool(n),"checkpoint changed retained mass presence");
    if(m){check(m->representation()==n->representation()&&m->liquidDensity()==n->liquidDensity()
        &&m->gasDensity()==n->gasDensity(),"checkpoint changed retained mass representation/materials");
        bits(m->masses(),n->masses(),"checkpoint changed retained M bits");
        const auto phase=[&](const std::vector<HydroPhaseMeasure>& x,const std::vector<HydroPhaseMeasure>& y){
            check(x.size()==y.size(),"checkpoint changed represented moment count");
            for(std::size_t k=0;k<x.size();++k)check(x[k].liquid==y[k].liquid&&x[k].gas==y[k].gas,
                "checkpoint changed an independent represented phase moment");};
        phase(m->representedPhaseVolumes(),n->representedPhaseVolumes());
        const auto& x=m->representedOverlapPhaseVolumes();const auto& y=n->representedOverlapPhaseVolumes();
        check(x.size()==y.size(),"checkpoint changed represented overlap count");
        for(std::size_t k=0;k<x.size();++k)phase(x[k],y[k]);
        for(int p=0;p<m->layout().pressureCount();++p)for(int component=0;component<2;++component)
            for(int q=0;q<2;++q)check(m->usesFiniteP0Partition(p,component,q)==n->usesFiniteP0Partition(p,component,q),
                "checkpoint changed canonical finite-P0 partition flags");
    }
    ++comparisons;
}
HydroNativeValues vortex(const QuadtreeMesh& mesh,const HydroDofLayout& layout,bool rz){
    std::vector<double> u(layout.velocityCount());
    const auto psi=[&](Vector2D p)->long double{
        const long double r=p.x,z=p.y;return .01L*r*r*(1-r)*(1-r)*z*z*(1-z)*(1-z)*(1+.3L*r+.7L*z);};
    for(int k=0;k<layout.velocityCount();++k){const auto f=layout.domain().face(mesh,layout.nativeFace(k));
        if(f.negativeCell<0||f.positiveCell<0)continue;
        const auto& edge=f.horizontal?mesh.hEdges(f.edgeLabel):mesh.vEdges(f.edgeLabel);
        auto a=mesh.verts(edge.edgeToVert[0]).coord(),b=mesh.verts(edge.edgeToVert[1]).coord();
        if(f.horizontal?a.x>b.x:a.y>b.y)std::swap(a,b);
        u[k]=static_cast<double>((f.horizontal?-1:1)*(psi(b)-psi(a))*(rz?CartesianMeasure::twoPi:1)/f.area);
    }
    return HydroNativeValues(layout,u,std::vector<double>(layout.pressureCount()));
}
void run(){
    const double dt=.000244140625; // Exact binary clock interval.
    for(bool rz:{false,true}){
        auto mesh=makeMesh(4,1);const auto domain=domainFor(*mesh,rz);
        check(!domain.uniformLeaves(),"checkpoint qualification needs actual graded AMR leaves");
        std::vector<double> alpha(mesh->cellIndexNum());
        for(int ci=0;ci<mesh->cellIndexNum();++ci){const auto& cell=mesh->cells(mesh->cellIndexToLabel[ci]);
            const auto lo=mesh->verts(cell.cellToVert[0]).coord(),hi=mesh->verts(cell.cellToVert[2]).coord();
            alpha[ci]=std::clamp((.43-lo.y)/(hi.y-lo.y),0.,1.);}
        const auto publication=I::GeometricVofPublication::prepareLocalPlic(*mesh,domain,{},alpha);
        const auto layout=HydroDofLayout::build(*mesh,domain);
        const auto metrics=HydroPhysicalDivergence::build(*mesh,HydroIntegralDivergence::build(*mesh,layout));
        const auto values=vortex(*mesh,layout,rz);
        const auto initial=HydroTwoPhaseState::createVof(*mesh,HydroPhaseGeometry::build(*mesh,metrics,publication),values);
        auto restoredMesh=mesh->clone();
        const auto initialBytes=initial.checkpoint(*mesh);
        const auto restoredInitial=HydroTwoPhaseState::restoreCheckpoint(*restoredMesh,initialBytes);
        compare(initial,restoredInitial);check(restoredInitial.checkpoint(*restoredMesh)==initialBytes,
            "history-free initial checkpoint roundtrip changed bytes");
        HydroVofAcConfig config;config.liquidDensity=10;config.gasDensity=1;
        config.liquidViscosity=.01;config.gasViscosity=.001;
        const auto firstStep=HydroVofStep::prepare(*mesh,initial,{1,0,dt,Clock::Position::End},config);++physicalSteps;
        FrozenPrepareContract::verify(*mesh,firstStep.frozenMidpoint(),config.saddle,"checkpoint_graded_startup",rz);
        ++FrozenPrepareContract::checkpointReferences;
        if(rz){
            // Same original checkpoint specimen, genuine independent restored owner.
            // This complete Step is compared numerically, not given a new whole-Step FP promise.
            auto coldMesh=mesh->clone();const auto coldInitial=HydroTwoPhaseState::restoreCheckpoint(*coldMesh,initialBytes);
            check(!coldInitial.momentumGeometry().sameIdentity(initial.momentumGeometry()),"cold checkpoint restored an old matrix-plan owner");
            const auto coldStep=HydroVofStep::prepare(*coldMesh,coldInitial,{1,0,dt,Clock::Position::End},config);
            FrozenPrepareContract::payload(firstStep.frozenMidpoint(),coldStep.frozenMidpoint());
            FrozenPrepareContract::verify(*coldMesh,coldStep.frozenMidpoint(),config.saddle,"checkpoint_real_cloned_cold_startup");
            ++FrozenPrepareContract::coldCheckpointReferences;
        }
        auto original=firstStep.candidateState();check(original.velocityHistory()&&original.retainedMass(),
            "checkpoint scene lacks real physical End history and retained M");
        const auto bytes=original.checkpoint(*mesh,&initial);
        auto resumed=HydroTwoPhaseState::restoreCheckpoint(*restoredMesh,bytes);
        compare(original,resumed);check(!resumed.sameIdentity(original)
            &&resumed.phaseGeometry().vofIdentity()!=original.phaseGeometry().vofIdentity()
            &&resumed.velocityHistory()->vofSource!=original.velocityHistory()->vofSource,
            "checkpoint reused old physical geometry/history identities");
        check(resumed.checkpoint(*restoredMesh)==bytes,
            "restored history lost its genuine predecessor for repeated checkpoint capture");
        const auto repeated=HydroTwoPhaseState::restoreCheckpoint(*restoredMesh,resumed.checkpoint(*restoredMesh));
        compare(resumed,repeated);
        refuse([&]{(void)original.checkpoint(*mesh);},"fresh state silently dropped history without its actual predecessor");
        const auto foreignPublication=I::GeometricVofPublication::prepareLocalPlic(*mesh,domain,{},alpha);
        const auto foreignPrevious=HydroTwoPhaseState::createVof(*mesh,HydroPhaseGeometry::build(*mesh,metrics,foreignPublication),values);
        refuse([&]{(void)original.checkpoint(*mesh,&foreignPrevious);},"checkpoint accepted a foreign predecessor publication owner");
        auto changedU=values.normalVelocity();for(auto& value:changedU)value*=2;
        const auto changedPrevious=HydroTwoPhaseState::createVof(*mesh,initial.phaseGeometry(),HydroNativeValues(layout,changedU,values.pressure()));
        refuse([&]{(void)original.checkpoint(*mesh,&changedPrevious);},"checkpoint accepted different actual predecessor U bits");
        auto corrupted=bytes;corrupted[corrupted.size()/2]^=1;
        refuse([&]{(void)HydroTwoPhaseState::restoreCheckpoint(*restoredMesh,corrupted);},"checkpoint accepted a corrupted payload");
        refuse([&]{(void)HydroTwoPhaseState::restoreCheckpoint(*restoredMesh,bytes.substr(0,bytes.size()-1));},"checkpoint accepted truncated native data");
        auto foreignMesh=mesh->clone();auto& vertex=foreignMesh->verts(foreignMesh->vertIndexToLabel[0]);
        vertex.setCoord(vertex.coord()+Vector2D{0,.0009765625});
        refuse([&]{(void)HydroTwoPhaseState::restoreCheckpoint(*foreignMesh,bytes);},"checkpoint accepted equal-size foreign native geometry");
        const auto& currentFrozen=firstStep.frozenMidpoint();
        const auto reprepare=[&](const QuadtreeMesh& actualMesh,LAMP::TimeLayer clock,HydroFrozenMidpointForces forces){
            return HydroFrozenMidpoint::prepare(actualMesh,currentFrozen.material(),currentFrozen.convection(),currentFrozen.appliedConvection(),
                currentFrozen.viscosity(),currentFrozen.pressureLoad(),currentFrozen.beginValues(),std::move(forces),currentFrozen.boundary(),
                currentFrozen.scheme(),currentFrozen.material().sources(),clock,config.saddle);};
        refuse([&]{(void)reprepare(*foreignMesh,currentFrozen.material().timeLayer(),currentFrozen.forces());},
            "Frozen fused preparation accepted equal-size stale actual native geometry");
        auto wrongClock=currentFrozen.material().timeLayer();wrongClock.end+=dt;
        refuse([&]{(void)reprepare(*mesh,wrongClock,currentFrozen.forces());},"Frozen fused preparation accepted a foreign complete midpoint clock");
        for(double bad:{std::numeric_limits<double>::infinity(),std::numeric_limits<double>::quiet_NaN()}){
            auto forces=currentFrozen.forces();forces.bodyAcceleration.x=bad;
            refuse([&]{(void)reprepare(*mesh,currentFrozen.material().timeLayer(),forces);},"Frozen fused preparation accepted nonfinite actual force data");}

        refuse([&]{(void)HydroTwoPhaseState::restoreCheckpoint(*restoredMesh,bytes,HydroVofPrescribedFlow{});},"checkpoint accepted a different physical boundary capability");
        // Continue both the uninterrupted and independently restored owners.
        for(int step=2;step<=3;++step){const Clock target{step,(step-1)*dt,step*dt,Clock::Position::End};
            const auto next=HydroVofStep::prepare(*mesh,original,target,config);
            if(step==2){FrozenPrepareContract::verify(*mesh,next.frozenMidpoint(),config.saddle,"checkpoint_graded_history");
                ++FrozenPrepareContract::checkpointReferences;}
            const auto split=HydroVofStep::prepare(*restoredMesh,resumed,target,config);physicalSteps+=2;
            check(next.report().usedVelocityHistory&&split.report().usedVelocityHistory
                &&next.report().stages.size()==split.report().stages.size(),"restart changed actual AC history/correction decisions");
            compare(next.candidateState(),split.candidateState());original=next.candidateState();resumed=split.candidateState();
        }
        // A real changed-tree remap intentionally resets incompatible history.
        auto targetMesh=mesh->clone();DictionaryMesh::RefineBox box;box.min={.375,.375};box.max={.625,.625};box.level=2;
        targetMesh->applyRefineBoxes({box});check(targetMesh->cellIndexNum()!=mesh->cellIndexNum(),"checkpoint AMR scene did not actually change topology");
        auto interfaceTransfer=I::GeometricVofMeshTransfer::prepare(*mesh,original.phaseGeometry().vofPublication(),
            *targetMesh,domainFor(*targetMesh,rz),original.timeLayer());
        HydroVofMeshTransferConfig remapConfig;remapConfig.momentumDensityRepresentation=HydroMomentumMass::DensityRepresentation::CellAverageProjection;
        const auto remap=HydroVofMeshTransfer::prepare(*mesh,original,*targetMesh,std::move(interfaceTransfer),10,1,remapConfig);
        const auto remapped=remap.candidateState();check(!remapped.velocityHistory()&&remapped.retainedMass(),"real AMR did not reset history and retain target M");
        const auto changedBytes=remapped.checkpoint(*targetMesh,&original); // Old-mesh predecessor is deliberately not consumed.
        auto targetClone=targetMesh->clone();const auto remapRestored=HydroTwoPhaseState::restoreCheckpoint(*targetClone,changedBytes);
        compare(remapped,remapRestored);check(!remapRestored.velocityHistory()&&remapRestored.checkpoint(*targetClone)==changedBytes,
            "AMR restart invented history or changed accepted bytes");
        const Clock afterRemap{4,3*dt,4*dt,Clock::Position::End};
        const auto next=HydroVofStep::prepare(*targetMesh,remapped,afterRemap,config);
        FrozenPrepareContract::verify(*targetMesh,next.frozenMidpoint(),config.saddle,"checkpoint_changed_AMR_startup");
        ++FrozenPrepareContract::checkpointReferences;
        const auto split=HydroVofStep::prepare(*targetClone,remapRestored,afterRemap,config);physicalSteps+=2;
        check(!next.report().usedVelocityHistory&&!split.report().usedVelocityHistory,"AMR restart changed the required startup correction");
        compare(next.candidateState(),split.candidateState());
        // An independently stored minor stock remains present even though its
        // display alpha rounds to one. The retained finite-P0 M must also agree.
        std::vector<std::array<long double,2>> tiny(mesh->cellIndexNum());
        for(int ci=0;ci<mesh->cellIndexNum();++ci)tiny[ci]={domain.cellVolumes()[domain.cellRowByIndex()[ci]],std::numeric_limits<long double>::denorm_min()};
        const auto tinyPublication=I::GeometricVofPublication::prepareLocalPlicInventory(*mesh,domain,
            {1,0,dt,Clock::Position::End},tiny,{});
        const auto tinyPhase=HydroPhaseGeometry::build(*mesh,metrics,tinyPublication);
        const auto tinyMass=HydroMomentumMass::buildCellAverageProjection(*mesh,initial.momentumGeometry(),tinyPhase,10,1);
        const HydroNativeValues zero(layout,std::vector<double>(layout.velocityCount()),std::vector<double>(layout.pressureCount()));
        const auto tinyState=HydroTwoPhaseState::createVofRemapped(*mesh,tinyPhase,zero,tinyMass);
        const auto tinyBytes=tinyState.checkpoint(*mesh);const auto tinyRestored=HydroTwoPhaseState::restoreCheckpoint(*restoredMesh,tinyBytes);
        compare(tinyState,tinyRestored);check(tinyMass.report().finiteP0PhasePartitions>0,"checkpoint tiny mass did not exercise canonical finite-P0 moments");
        for(const auto& cell:tinyRestored.phaseGeometry().vofPublication().plicCells())
            check(cell.inventory[1]==std::numeric_limits<long double>::denorm_min()&&cell.inventory[1]>0,
                "checkpoint reconstructed a tiny independent gas stock from rounded alpha");
        std::cout<<"VOF_CHECKPOINT_SCENE,"<<(rz?"RZ":"XY")<<",graded_leaves,"<<mesh->cellIndexNum()
            <<",changed_leaves,"<<targetMesh->cellIndexNum()<<",history_and_mass_exact,1\n";
    }
    check(FrozenPrepareContract::checkpointReferences==6&&FrozenPrepareContract::threadGroups==1&&FrozenPrepareContract::threadPairs==2&&
        FrozenPrepareContract::coldCheckpointReferences==1,
        "Frozen checkpoint graded/history/changed-AMR XY/RZ reference and own-request FP coverage incomplete");
    std::cout<<"VOF_CHECKPOINT_PASS,XY_RZ_graded_and_changed_AMR,comparisons,"<<comparisons
        <<",physical_steps,"<<physicalSteps<<",guards,"<<guards<<",native_exact,1\n";
}
}

// F13 additive namespace: insert after existing namespaces, before main.
// Add only the selector branch shown below; retain every original selector/default.
namespace NoSlipRawOwnerContract {
using namespace LAMP::Hydrodynamics;
namespace I=LAMP::Interface;
using Floating=Task38StateValidationContract::Floating38;
using Restore=Task38StateValidationContract::FloatingRestore38;
using Matrix=HydroEnergyViscosity::Matrix;
std::size_t pairs=0,coldGuards=0,errorGuards=0,muChecks=0,zeroPlanRawReuse=0;
std::unique_ptr<QuadtreeMesh> meshFor(bool changed){
    // Original horizontal no-slip pressure channel: selected fluid y<1,
    // complete planar walls at y=0/1, exterior pressure ports at R=.5/2.5.
    DictionaryMesh d;d.Nx=8;d.Ny=6;d.minCoord={.5,0};d.maxCoord={2.5,1.5};d.maxCellNum=4096;d.maxInitialIterations=0;
    auto mesh=std::make_unique<QuadtreeMesh>(d);DictionaryMesh::RefineBox b;
    b.level=1;b.min={1,.25};b.max={2,.75};std::vector<DictionaryMesh::RefineBox> boxes{b};
    if(changed){b.min={.5,.25};b.max={.75,.5};boxes.push_back(b);}mesh->applyRefineBoxes(boxes);
    return mesh;
}
HydroPhysicalDivergence metricsFor(const QuadtreeMesh& mesh){
    FaceSegmentGraph graph;graph.build(mesh,[](const QuadtreeMesh&,const Cell& c){return c.centerCoord().y<1;});
    const auto domain=NativeLeafDomain::build(mesh,graph,true);
    check(!domain.uniformLeaves()&&domain.cellIndexByRow().size()<static_cast<std::size_t>(mesh.cellIndexNum()),"F13 needs actual selected RZ2:1 channel");
    return HydroPhysicalDivergence::build(mesh,HydroIntegralDivergence::build(mesh,HydroDofLayout::build(mesh,domain)));
}
std::vector<int> wallsFor(const QuadtreeMesh& mesh,const HydroPhysicalDivergence& metrics,bool bottomOnly=false){
    std::vector<int> rows;const auto& l=metrics.layout();for(int u=0;u<l.velocityCount();++u){
        const auto face=l.domain().face(mesh,l.nativeFace(u));
        if(face.horizontal&&(face.center.y==0||(!bottomOnly&&face.center.y==1))&&metrics.faceMeasures()[u]>0)rows.push_back(u);
    }check(!rows.empty(),"F13 missing physical complete no-slip plane");return rows;
}
HydroPhaseGeometry phaseFor(const QuadtreeMesh& mesh,const HydroPhysicalDivergence& metrics,double front){
    std::vector<double> alpha(mesh.cellIndexNum());for(int c:metrics.layout().domain().cellIndexByRow()){
        const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);const auto lo=mesh.verts(cell.cellToVert[0]).coord(),hi=mesh.verts(cell.cellToVert[2]).coord();
        alpha[c]=std::clamp((front-lo.y)/(hi.y-lo.y),0.,1.);
    }const auto publication=I::GeometricVofPublication::prepareLocalPlic(mesh,metrics.layout().domain(),{},alpha);
    return HydroPhaseGeometry::build(mesh,metrics,publication);
}
template<class T>void bits(const std::vector<T>& a,const std::vector<T>& b,const char* message){
    check(a.size()==b.size()&&(a.empty()||std::memcmp(a.data(),b.data(),a.size()*sizeof(T))==0),message);
}
void matrixBits(const Matrix& a,const Matrix& b){
    check(a.isCompressed()&&b.isCompressed()&&a.rows()==b.rows()&&a.cols()==b.cols()&&a.nonZeros()==b.nonZeros(),"F13 sparse inventory differs from original cold build");
    check(std::memcmp(a.outerIndexPtr(),b.outerIndexPtr(),(a.outerSize()+1)*sizeof(Matrix::StorageIndex))==0,"F13 sparse outer order changed");
    if(a.nonZeros())check(std::memcmp(a.innerIndexPtr(),b.innerIndexPtr(),a.nonZeros()*sizeof(Matrix::StorageIndex))==0
        &&std::memcmp(a.valuePtr(),b.valuePtr(),a.nonZeros()*sizeof(double))==0,"F13 sparse values/order differ from original cold build");
}
void prime(){
    check(std::feclearexcept(FE_ALL_EXCEPT)==0&&std::feraiseexcept(FE_INEXACT)==0,"F13 cannot establish actual sticky native inexact");
#if defined(_M_X64) || defined(__SSE2__)
    _mm_setcsr(_mm_getcsr()|0x20u);
#endif
    errno=173;
}
void currentMaterial(const HydroPhaseViscosity& a,const HydroPhaseGeometry& phase,double liquid,double gas){
    check(a.materialModel()==HydroPhaseViscosity::MaterialModel::VofCellAverage&&a.geometry().vofIdentity()==phase.vofIdentity()
        &&sameClock(a.geometry().timeLayer(),phase.timeLayer())&&!a.hasQuadrature(),"F13 material lost current genuine LocalPlic owner/clock");
    const auto& factor=a.viscosity();const auto& raw=factor.strain();const auto& source=phase.vofPublication();
    check(!raw.isStressCompatible()&&raw.components()==4&&!raw.points().empty(),"F13 bypassed actual RZ wall raw strain");
    for(std::size_t q=0;q<raw.points().size();++q){const auto& point=raw.points()[q];const int c=phase.layout().nativeCell(point.cellDof);
        const auto& stock=source.plicCells()[c].inventory;const auto capacity=source.cells()[c].capacity;
        const double expected=static_cast<double>((stock[0]/capacity)*liquid+(stock[1]/capacity)*gas);
        check(factor.viscosity()[q]==expected,"F13 cached old inventory or phase coefficient");
        const double weight=static_cast<double>(2*static_cast<long double>(expected)*point.physicalWeight);
        for(int k=0;k<4;++k)check(factor.rowWeights()[4*q+k]==weight,"F13 stale physical RZ material row weight");++muChecks;
    }
}
HydroPhaseViscosity pair(const QuadtreeMesh& mesh,const HydroPhaseGeometry& phase,double liquid,double gas,
    const std::vector<int>& walls,const HydroMomentumGeometry& owner){
    const Floating entry;const auto actual=HydroPhaseViscosity::buildVofCellAverage(mesh,phase,liquid,gas,walls,owner);const Floating after;
    entry.restore();const auto cold=HydroPhaseViscosity::buildVofCellAverage(mesh,phase,liquid,gas,walls);const Floating coldAfter;
    check(after.same(coldAfter),"F13 owner reuse changed full native/raw FP/errno from identical original cold entry");
    check(actual.viscosity().strain().stationaryNoSlipVelocityDofs()==walls,"F13 raw strain retained wrong physical wall rows");
    matrixBits(actual.viscosity().strain().matrix(),cold.viscosity().strain().matrix());matrixBits(actual.viscosity().matrix(),cold.viscosity().matrix());
    bits(actual.viscosity().viscosity(),cold.viscosity().viscosity(),"F13 mu bits differ from original cold material");
    bits(actual.viscosity().rowWeights(),cold.viscosity().rowWeights(),"F13 W bits differ from original cold material");
    const auto& p=actual.viscosity().strain().points();const auto& q=cold.viscosity().strain().points();check(p.size()==q.size(),"F13 stale volume quadrature inventory");
    for(std::size_t i=0;i<p.size();++i)check(p[i].cellDof==q[i].cellDof&&std::memcmp(&p[i].point.x,&q[i].point.x,sizeof(double))==0
        &&std::memcmp(&p[i].point.y,&q[i].point.y,sizeof(double))==0&&std::memcmp(&p[i].physicalWeight,&q[i].physicalWeight,sizeof(double))==0,"F13 changed original ordered volume quadrature");
    currentMaterial(actual,phase,liquid,gas);after.restore();++pairs;return actual;
}
void firstProduct(const HydroPhaseViscosity& a){const auto& r=a.viscosity().productReport();
    check(r.requested&&r.eligible&&r.currentFpQualified&&!r.sourceAvailable&&!r.sourceFpReceipt&&!r.localUsed
        &&r.fullAssembly&&r.memoPublished&&r.fallbackReason==0&&r.outputRows>0
        &&r.recomputedRows==r.outputRows&&r.reusedRows==0,
        "F13 cold raw owner did not compute and publish its actual first numeric row receipt");
}
void allReused(const HydroPhaseViscosity& a){const auto& r=a.viscosity().productReport();
    check(r.requested&&r.eligible&&r.sourceAvailable&&r.sourceFpReceipt&&r.currentFpQualified&&r.localUsed&&!r.fullAssembly&&r.memoPublished
        &&r.changedWeightRows==0&&r.recomputedRows==0&&r.reusedRows==r.outputRows&&r.outputRows>0,"F13 actual same-weight call did not reuse complete K row inventory");
}
void rejectPair(const QuadtreeMesh& mesh,const HydroPhaseGeometry& phase,double liquid,double gas,
    const std::vector<int>& walls,const HydroMomentumGeometry& owner){
    const auto error=[&](bool typed){try{if(typed)HydroPhaseViscosity::buildVofCellAverage(mesh,phase,liquid,gas,walls,owner);
        else HydroPhaseViscosity::buildVofCellAverage(mesh,phase,liquid,gas,walls);}catch(const std::invalid_argument& e){return std::string(e.what());}
        throw std::runtime_error("F13 missing original material rejection");};
    const Floating entry;const auto actual=error(true);const Floating after;entry.restore();const auto cold=error(false);const Floating coldAfter;
    check(actual==cold&&after.same(coldAfter),"F13 cache bypass changed original material rejection order/message/full FP");after.restore();++errorGuards;
}
void run(){
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(__FAST_MATH__) && (!defined(_M_FP_FAST) || !_M_FP_FAST) && !defined(__AVX__) && !defined(__AVX2__) && !defined(__AVX512F__) && !defined(__FMA__)
    Restore original;check(std::fesetround(FE_TONEAREST)==0,"F13 cannot establish actual RNE");_mm_setcsr(_mm_getcsr()&~0x8040u);
    auto mesh=meshFor(false);const auto metrics=metricsFor(*mesh);const auto walls=wallsFor(*mesh,metrics);
    const auto phase=phaseFor(*mesh,metrics,.5625),changed=phaseFor(*mesh,metrics,.59375);
    check(phase.vofIdentity()!=changed.vofIdentity(),"F13 changed phase reused publication identity");
    std::size_t changedStocks=0;for(std::size_t c=0;c<phase.vofPublication().plicCells().size();++c)
        changedStocks+=phase.vofPublication().plicCells()[c].inventory!=changed.vofPublication().plicCells()[c].inventory;
    check(changedStocks>0,"F13 has no actually changed independent L/G stocks");
    const auto owner=HydroMomentumGeometry::build(*mesh,metrics);
    prime();const auto first=pair(*mesh,phase,.3,.15,walls,owner);firstProduct(first);
    prime();const auto seed=pair(*mesh,phase,.3,.15,walls,owner);const auto& seedReport=seed.viscosity().productReport();
    check(&first.viscosity().strain().matrix()==&seed.viscosity().strain().matrix(),"F13 quiet same-owner second call rebuilt original raw E");
    allReused(seed);check(&seed.viscosity().matrix()==&first.viscosity().matrix(),"F13 all-clean owned K did not share actual immutable matrix");
    prime();const auto stock=pair(*mesh,changed,.3,.15,walls,owner);
    check(&stock.viscosity().strain().matrix()==&first.viscosity().strain().matrix()
        &&stock.viscosity().viscosity()!=first.viscosity().viscosity(),"F13 changing actual stocks rebuilt raw geometry or retained old mu");
    const auto& stockReport=stock.viscosity().productReport();
    check(stockReport.eligible&&stockReport.sourceAvailable&&stockReport.sourceFpReceipt&&stockReport.changedWeightRows>0
        &&stockReport.recomputedRows>0&&stockReport.reusedRows+stockReport.recomputedRows==stockReport.outputRows&&stockReport.memoPublished,"F13 changed-stock K report lost actual dirty/output inventory");
    prime();const auto coefficients=pair(*mesh,changed,.6,.3,walls,owner);
    check(&coefficients.viscosity().strain().matrix()==&first.viscosity().strain().matrix()
        &&coefficients.viscosity().viscosity()!=stock.viscosity().viscosity(),"F13 changing phase coefficients reused old material or rebuilt E");
    prime();const auto repeated=pair(*mesh,changed,.6,.3,walls,owner);allReused(repeated);check(&repeated.viscosity().matrix()==&coefficients.viscosity().matrix(),"F13 changed-material repeat copied actual all-clean K stock");
    prime();rejectPair(*mesh,phase,std::numeric_limits<double>::quiet_NaN(),.15,walls,owner);
    prime();rejectPair(*mesh,phase,std::numeric_limits<double>::denorm_min(),0,walls,owner);
    const auto independent=HydroMomentumGeometry::build(*mesh,metrics);check(!independent.sameIdentity(owner),"F13 needs fresh Momentum Data with identical metrics");
    prime();const auto ownerCold=pair(*mesh,phase,.3,.15,walls,independent);firstProduct(ownerCold);
    check(&ownerCold.viscosity().strain().matrix()!=&first.viscosity().strain().matrix(),"F13 new Momentum Data inherited another owner's raw slot");++coldGuards;
    // Actual all-gas stocks with a nonzero liquid coefficient traverse the
    // original raw wall factory, not the liquid=gas=0 inactive shortcut.
    const auto gasPhase=phaseFor(*mesh,metrics,0.);
    const auto zeroOwner=HydroMomentumGeometry::build(*mesh,metrics);
    check(!zeroOwner.sameIdentity(owner),"F13 zero-weight raw witness requires an independent actual owner");
    prime();const auto zeroFirst=pair(*mesh,gasPhase,.3,0.,walls,zeroOwner);
    const auto& zeroFactor=zeroFirst.viscosity();const auto& zeroReport=zeroFactor.productReport();
    check(zeroFactor.strain().matrix().nonZeros()>0&&!zeroFactor.strain().points().empty()
        &&std::all_of(zeroFactor.viscosity().begin(),zeroFactor.viscosity().end(),[](double x){return x==0;})
        &&std::all_of(zeroFactor.rowWeights().begin(),zeroFactor.rowWeights().end(),[](double x){return x==0;})
        &&zeroFactor.matrix().nonZeros()==0&&!zeroReport.requested&&!zeroReport.memoPublished&&zeroReport.symbolicStoredValues==0,
        "F13 all-gas zero weights did not retain a genuine raw factory without a numeric K plan");
    prime();const auto positive=pair(*mesh,phase,.3,0.,walls,zeroOwner);firstProduct(positive);
    check(&positive.viscosity().strain().matrix()==&zeroFactor.strain().matrix()
        &&std::any_of(positive.viscosity().viscosity().begin(),positive.viscosity().viscosity().end(),[](double x){return x>0;}),
        "F13 zero-weight raw geometry was rebuilt or its later actual positive material was lost");++zeroPlanRawReuse;
    const auto bottom=wallsFor(*mesh,metrics,true);check(bottom!=walls,"F13 changed wall plane is ineffective");
    prime();const auto wallCold=pair(*mesh,phase,.3,.15,bottom,owner);firstProduct(wallCold);
    check(&wallCold.viscosity().strain().matrix()!=&first.viscosity().strain().matrix(),"F13 different actual walls hit original raw slot");++coldGuards;
    auto target=meshFor(true);const auto targetMetrics=metricsFor(*target);check(target->cellIndexNum()>mesh->cellIndexNum(),"F13 AMR guard lacks actual changed leaves");
    const auto targetPhase=phaseFor(*target,targetMetrics,.5625);const auto targetOwner=HydroMomentumGeometry::build(*target,targetMetrics);
    prime();const auto amrCold=pair(*target,targetPhase,.3,.15,wallsFor(*target,targetMetrics),targetOwner);firstProduct(amrCold);
    check(&amrCold.viscosity().strain().matrix()!=&first.viscosity().strain().matrix(),"F13 new AMR Data hit original raw slot");++coldGuards;
    // Fresh guard owner isolates a complete-entry mismatch from prior wall eviction.
    const auto guardOwner=HydroMomentumGeometry::build(*mesh,metrics);
    prime();const auto guardSource=pair(*mesh,phase,.3,.15,walls,guardOwner);
    prime();const auto guardHit=pair(*mesh,phase,.3,.15,walls,guardOwner);
    check(&guardSource.viscosity().strain().matrix()==&guardHit.viscosity().strain().matrix(),"F13 guard source never obtained actual quiet raw reuse");
    prime();_mm_setcsr(_mm_getcsr()|0x02u);const auto csrCold=pair(*mesh,phase,.3,.15,walls,guardOwner);firstProduct(csrCold);
    check(&csrCold.viscosity().strain().matrix()!=&guardHit.viscosity().strain().matrix(),"F13 different actual raw-CSR sticky status hit original raw slot");++coldGuards;
    prime();_mm_setcsr(_mm_getcsr()|0x02u);errno=174;const auto errnoCold=pair(*mesh,phase,.3,.15,walls,guardOwner);firstProduct(errnoCold);
    check(&errnoCold.viscosity().strain().matrix()!=&csrCold.viscosity().strain().matrix(),"F13 different complete errno entry hit old raw slot");++coldGuards;
    prime();check(std::fesetround(FE_UPWARD)==0&&std::fegetround()==FE_UPWARD,"F13 cannot establish actual unsupported directed controls");
    const auto directed=pair(*mesh,phase,.3,.15,walls,guardOwner);
    check(&directed.viscosity().strain().matrix()!=&errnoCold.viscosity().strain().matrix()&&!directed.viscosity().productReport().memoPublished,
        "F13 unsupported directed full entry reused qualified raw/product state");++coldGuards;
    check(pairs==15&&coldGuards==6&&errorGuards==2&&muChecks>0&&zeroPlanRawReuse==1,"F13 missing actual owner/material/cold/FP gates");
    std::cout<<"NO_SLIP_RAW_OWNER_PASS,cold_pairs,"<<pairs<<",changed_stock_cells,"<<changedStocks<<",material_samples,"<<muChecks
        <<",cold_guards,"<<coldGuards<<",original_error_pairs,"<<errorGuards<<",seeded_rows,"<<seedReport.outputRows
        <<",repeat_reused_rows,"<<repeated.viscosity().productReport().reusedRows<<",zero_weight_no_plan_then_positive_same_raw,"<<zeroPlanRawReuse
        <<",full_native_raw_FP_errno_pairs,"<<pairs+errorGuards<<",same_true_RZ_2to1_noSlip,1,original_three_step_noSlip_selector,unchanged\n";
#else
    std::cout<<"NO_SLIP_RAW_OWNER_UNAVAILABLE,requires_original_strict_MSVC_x64_OpenMP_product_path\n";
#endif
}
}

int main(int argc,char** argv)try
{
    if(argc==2&&std::string(argv[1])=="--vof-no-slip-raw-owner-only"){NoSlipRawOwnerContract::run();return 0;}
    static_assert(!std::is_default_constructible_v<HydroTwoPhaseState>);
    static_assert(std::is_same_v<decltype(std::declval<const HydroTwoPhaseState&>().values()),const HydroNativeValues&>);
    static_assert(std::is_same_v<decltype(std::declval<const HydroTwoPhaseState&>().publication()),const HermiteGalsPublication&>);
    if(argc==2&&std::string(argv[1])=="--vof-checkpoint-only"){VofCheckpointContract::run();return 0;}
    if(argc==2&&std::string(argv[1])=="--vof-no-slip-pressure-only"){VofNoSlipPressureContract::run();return 0;}
    if(argc==2&&std::string(argv[1])=="--vof-pressure-driven-only"){VofPressureDrivenContract::run();return 0;}
    if(argc==2&&std::string(argv[1])=="--state-validation-batch-only"){Task38StateValidationContract::run();return 0;}
    if(argc==2&&std::string(argv[1])=="--vof-transport-provider-only"){vofTransportProviderOnly();return 0;}
    if(argc>1&&(std::string(argv[1])=="--normal-difference-window"||std::string(argv[1])=="--local-pressure-transport-probe"||std::string(argv[1])=="--local-pressure-only"))
        throw std::invalid_argument("NormalDifference capillary research is retired; this option is unavailable");
    const bool heightWindowEntry=(argc==4||argc==5)&&std::string(argv[1])=="--height-function-increment-window";
    const bool replayEntry=argc==3&&std::string(argv[1])=="--first-jet-replay-only";
    const bool csfCouplingEntry=(argc==2||argc==3)&&std::string(argv[1])=="--csf-coupling-only";
    const bool projectedMotionEntry=(argc==2||argc==3)&&std::string(argv[1])=="--projected-motion-work-only";
    const bool projectedDensityEntry=argc==2&&std::string(argv[1])=="--projected-density-only";
    const bool capillaryModelEntry=argc==2&&std::string(argv[1])=="--capillary-model-only";
    const bool adaptiveOnly=argc==2&&std::string(argv[1])=="--adaptive-linear-only";
    const bool highRhoOnly=argc==2&&std::string(argv[1])=="--high-rho-only";
    const bool inviscidAcOnly=argc==2&&std::string(argv[1])=="--inviscid-ac-only";
    const bool rzEntry=argc==2&&std::string(argv[1])=="--rz-only";
    const bool prescribedEntry=argc==2&&std::string(argv[1])=="--prescribed-flow-only";
    const bool pressureOutletEntry=argc==2&&std::string(argv[1])=="--pressure-outlet-only";
    const bool bodyEntry=argc==2&&std::string(argv[1])=="--body-acceleration-only";
    const bool symmetryEntry=argc==2&&std::string(argv[1])=="--symmetry-boundary-only";
    check(argc==1||heightWindowEntry||replayEntry||projectedMotionEntry||csfCouplingEntry||projectedDensityEntry||capillaryModelEntry||adaptiveOnly||inviscidAcOnly||highRhoOnly||rzEntry||prescribedEntry||pressureOutletEntry||bodyEntry||symmetryEntry,"unknown two-phase step contract argument");
    std::cout<<std::setprecision(17);
    static_assert(std::is_same_v<decltype(std::declval<const HydroTwoPhaseStep&>().bodyLoads()),const std::vector<long double>&>);
    static_assert(std::is_same_v<decltype(std::declval<const HydroTwoPhaseStep&>().capillaryLoad()),const HydroBalancedCsfLoad&>);
    const HydroTwoPhaseStepConfig defaults;
    check(defaults.momentumDensityRepresentation==HydroMomentumMass::DensityRepresentation::SharpGeometry,"default density representation changed");
    check(defaults.capillaryModel==HydroCapillaryModel::BalancedCsf,"default force model changed");
    check(defaults.csfCurvature.relativeQuadratureTolerance==1e-7,"default local CSF integration budget changed");
    check(defaults.bodyAcceleration.x==0&&defaults.bodyAcceleration.y==0&&HydroTwoPhaseStepReport{}.bodyPower==0,"default body acceleration/power is not exactly zero");
    check(defaults.spaceTime.enableAtomLocalTimeIntegration,"default atom-local time integration changed");
    check(defaults.spaceTime.firstJetTimeScheme==HydroFirstJetTimeScheme::CharacteristicPullback
        &&HydroSpaceTimeConfig{}.firstJetTimeScheme==HydroFirstJetTimeScheme::BoundaryLaw,
        "default step pullback or standalone boundary-law time scheme changed");
    check(defaults.timeScheme==HydroTwoPhaseTimeScheme::AdaptiveLinearMidpoint
        &&HydroTwoPhaseIterationReport{}.timeScheme==defaults.timeScheme
        &&HydroTwoPhaseStepReport{}.timeScheme==defaults.timeScheme
        &&defaults.interfaceMismatchTolerance==.01&&defaults.maximumCorrectionRelative==.25
        &&defaults.correctionReferenceVelocity==1e-8&&defaults.capillaryCorrectionCourant==.1
        &&defaults.maximumCapillaryCourant==.5&&!defaults.resetVelocityHistory,"default time scheme or adaptive controls changed");
    check(defaults.velocityUpdateTolerance==1e-11&&defaults.geometryUpdateTolerance==1e-11
        &&defaults.momentumPowerTolerance==1e-11,"default strict corrector targets changed");
    if(heightWindowEntry){heightFunctionIncrementWindow(argv[2],std::stoi(argv[3]),argc==5?std::stod(argv[4]):.0002);return 0;}
    if(replayEntry){firstJetReplayOnly(argv[2]);return 0;}
    if(projectedMotionEntry){projectedMotionWorkOnly(argc==3?argv[2]:"");return 0;}
    if(projectedDensityEntry){projectedDensityOnly();return 0;}
    if(capillaryModelEntry){capillaryModelOnly();return 0;}
    if(csfCouplingEntry){csfCouplingOnly(argc==3?argv[2]:"");return 0;}
    if(symmetryEntry){symmetryBoundaryOnly();return 0;}
    if(bodyEntry){bodyAccelerationOnly();return 0;}
    if(rzEntry){rzOnly();return 0;}
    if(prescribedEntry){prescribedFlowOnly();return 0;}
    if(pressureOutletEntry){prescribedFlowOnly(true);return 0;}
    auto mesh=makeMesh(4);Case shape;shape.pureFirstJet=true;
    const auto geometry=HydroPhaseGeometry::build(*mesh,stepMetrics(*mesh),initial(*mesh,shape));
    if(highRhoOnly){movingCase(*mesh,geometry,1e6,0,1);
        check(FrozenPrepareContract::laggedReferences==1&&FrozenPrepareContract::threadGroups==1&&FrozenPrepareContract::threadPairs==2,
            "Frozen actual lagged high-rho reference and own-request FP coverage incomplete");return 0;}
    if(adaptiveOnly||inviscidAcOnly){
        linearTimeCase(*mesh,geometry,inviscidAcOnly);
        check(linearCases==5&&linearRejects==3&&movingCases==0&&acceptedSteps==0&&guardCases==0
            &&capillarySteps==0&&capillaryRejects==0,"incomplete direct default AC contract coverage");
        std::cout<<(inviscidAcOnly?"INVISCID_AC_ONLY,PASS cases=":"ADAPTIVE_LINEAR_ONLY,PASS cases=")<<linearCases<<" rejected="<<linearRejects
            <<" equation="<<peakEquation<<" energy="<<peakEnergy<<" material="<<peakMaterial<<'\n';return 0;
    }
    prescribedFlowOnly();prescribedFlowOnly(true);
    guards(*mesh,geometry);stationary(*mesh,geometry);
    for(double rho:{1.,1e3,1e6})movingCase(*mesh,geometry,rho,0,rho==1e3?2:1);
    movingCase(*mesh,geometry,1e3,.01,1);
    check(movingCases==4&&acceptedSteps==6&&guardCases==46,"incomplete two-phase step contract coverage");
    check(relaxedComparisons==1&&updateToleranceRejects==12,"incomplete single relaxed comparison or separate tolerance rejection coverage");
    capillaryCase(*mesh,geometry);
    check(capillarySteps==1&&capillaryRejects==2&&movingCases==4&&acceptedSteps==6&&guardCases==46,"incomplete independent capillary step/guard coverage");
    std::cout<<"CAPILLARY_SUMMARY,"<<capillarySteps<<','<<capillaryRejects<<",PASS\n";
    linearTimeCase(*mesh,geometry);
    check(linearCases==5&&linearRejects==3&&movingCases==4&&acceptedSteps==6&&guardCases==46
        &&capillarySteps==1&&capillaryRejects==2,"linear scheme changed original family or refusal coverage");
    bodyAccelerationOnly();symmetryBoundaryOnly();capillaryModelOnly();projectedDensityOnly();
    std::cout<<"H14A1 TWO_PHASE_STEP PASS moving="<<movingCases<<" steps="<<acceptedSteps<<" rejected="<<guardCases
        <<" equation="<<peakEquation<<" energy="<<peakEnergy<<" material="<<peakMaterial<<'\n';return 0;
}
catch(const std::exception& error){std::cerr<<"H14A1 TWO_PHASE_STEP FAIL: "<<error.what()<<'\n';return 1;}
