#include <chrono> // F23 passive small-fixture timing only.
#include <cfenv>
#include <cerrno>
#include <cstring>
#include <optional>
#include <xmmintrin.h>
#if defined(_OPENMP)
#include <omp.h>
#endif
#include "Modules/Hydrodynamics/HydroPhaseViscosity.h" // F24 explicit public typed material.
#include "Modules/Hydrodynamics/HydroVofMeshTransfer.h" // Public API stands alone.
#include "Modules/Hydrodynamics/HydroVofRefinement.h"
#include "Core/Mesh/CartesianLeafLocator.h"
#include "Modules/Hydrodynamics/Detail/HydroRegularAdvectorGeometry.h"
#include "Core/Geometry/Detail/FiniteVolumeClosure.h"
#include "Electrohydrodynamics/IO/VofElectrohydrodynamicsOutput.h"
#include "Hydrodynamics/HydrodynamicsRun.h"
#include <filesystem>
#include <fstream>
#include <sstream>
#include <limits>
#include "Modules/Hydrodynamics/HydroVofStep.h"
#include "Modules/Electric/ElectricVofChargeTransport.h"
#include "Modules/Interface/GeometricVofInitialization.h"
#include <iomanip>
#include "Core/Field/CellField.h"
#include "Core/Field/VertField.h"
#include "Core/Field/HEdgeField.h"
#include "Core/Field/VEdgeField.h"
#include "Core/Mesh/CartesianMeasure.h"
#include "Core/Mesh/FaceSegmentGraph.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>

namespace
{
using namespace LAMP::Hydrodynamics;
using namespace LAMP::Interface;
using Density=HydroMomentumMass::DensityRepresentation;
using Clock=LAMP::TimeLayer;
unsigned assertions=0,remaps=0,steps=0;
void check(bool ok,const char* why){++assertions;if(!ok)throw std::runtime_error(why);}
void near(long double a,long double b,long double scale,const char* why)
{check(std::isfinite(a)&&std::isfinite(b)&&std::abs(a-b)<=2e-11L*std::max(1.L,scale),why);}
template<class F> void rejects(F&& f,const char* why)
{bool caught=false;try{f();}catch(const std::exception&){caught=true;}check(caught,why);}
bool sameClock(Clock a,Clock b)
{return a.step==b.step&&a.begin==b.begin&&a.end==b.end&&a.position==b.position;}
// F12: pair the real stage consumer on the existing Maxwell/no-slip fixtures.
namespace F12DiagnosticRateRows {
bool active=false;
std::size_t pairs=0,stages=0,rows=0,maxwellPairs=0,noSlipPairs=0,lateFailurePairs=0,lateProviderCalls=0;
int minimumWorkers=8,maximumWorkers=0;
void run();
template<class T>void exact(T a,T b,const char* why) {
    check(std::memcmp(&a,&b,sizeof(T))==0,why);
}
template<class Values>void sameValues(const Values& a,const Values& b,const char* why) {
    check(a.size()==b.size(),why);for(std::size_t i=0;i<a.size();++i)exact(a[i],b[i],why);
}
void sameNative(const HydroNativeValues& a,const HydroNativeValues& b) {
    sameValues(a.normalVelocity(),b.normalVelocity(),"F12 candidate/midpoint native U bits");
    sameValues(a.pressure(),b.pressure(),"F12 candidate/midpoint native p bits");
}
void samePublication(const GeometricVofPublication& a,const GeometricVofPublication& b) {
    check(sameClock(a.timeLayer(),b.timeLayer()),"F12 actual publication full clock");
    sameValues(a.alpha(),b.alpha(),"F12 actual phase alpha bits");
    check(a.plicCells().size()==b.plicCells().size()&&a.cells().size()==b.cells().size(),"F12 actual phase stock dimensions");
    for(std::size_t i=0;i<a.plicCells().size();++i){const auto& x=a.plicCells()[i];const auto& y=b.plicCells()[i];
        sameValues(x.inventory,y.inventory,"F12 actual independent liquid/gas stock bits");
        sameValues(x.measured,y.measured,"F12 actual measured phase bits");
        exact(x.normal.x,y.normal.x,"F12 actual PLIC normal x bits");exact(x.normal.y,y.normal.y,"F12 actual PLIC normal y bits");
        exact(x.offset,y.offset,"F12 actual PLIC offset bits");check(x.resolved==y.resolved,"F12 actual PLIC resolution");
        exact(a.cells()[i].capacity,b.cells()[i].capacity,"F12 actual phase capacity bits");}
}
void sameStage(const HydroVofStageReport& a,const HydroVofStageReport& b) {
    for(auto member:{&HydroVofStageReport::interfaceMismatch,&HydroVofStageReport::momentumPowerRelativeDefect})
        exact(a.*member,b.*member,"F12 actual stage physical double bits");
    for(auto member:{&HydroVofStageReport::maximumAbsDualDefect,&HydroVofStageReport::consistencyPower,
        &HydroVofStageReport::kineticChange,&HydroVofStageReport::capillaryWork,&HydroVofStageReport::maxwellWork,
        &HydroVofStageReport::boundaryKineticWork,
        &HydroVofStageReport::boundaryPhysicalKineticWork,&HydroVofStageReport::inletDissipationWork,
        &HydroVofStageReport::pressureOutletWork,&HydroVofStageReport::prescribedReactionWork,&HydroVofStageReport::boundaryFluxRoundingWork})
        exact(a.*member,b.*member,"F12 actual complete stage physical work bits");
    sameValues(a.boundaryMomentumRate,b.boundaryMomentumRate,"F12 actual boundary component momentum bits");
    using Energy=HydroVofStageReport::EnergyBudget;
    for(auto member:{&Energy::halfRWork,&Energy::mappedPrimalHalfRWork,&Energy::geometryTransferHalfRWork,&Energy::spatialExchangeHalfRWork,
        &Energy::timeMassWork,&Energy::roundedPhaseMassWork,&Energy::operatorRateRoundoffWork,&Energy::rawMinusCanonicalOutflowHalfWork,
        &Energy::halfRDecompositionResidual,&Energy::bodyWork,&Energy::viscousWork,&Energy::pressureWork,&Energy::momentumEquationWork,
        &Energy::kineticBudgetResidual,&Energy::closureResidual,&Energy::closureRelativeResidual,&Energy::stableGasHalfRWork,
        &Energy::stableGasSpatialExchangeHalfRWork,&Energy::stableGasMaxAbsIntegratedDefectOverMidMass,&Energy::peakGasAbsEndNormal,
        &Energy::peakGasIntegratedDefectOverMidMass,&Energy::peakGasMassChangeOverMidMass,&Energy::peakGasOutwardIntegralOverMidMass,
        &Energy::peakGasGeometryTransferOverMidMass,&Energy::peakGasSpatialExchangeOverMidMass})
        exact(a.energy.*member,b.energy.*member,"F12 actual complete endpoint energy ledger bits");
    check(a.energy.stableGasDofs==b.energy.stableGasDofs&&a.energy.peakGasNativeFace==b.energy.peakGasNativeFace,"F12 actual gas diagnostic row");
    for(auto member:{&HydroSaddleReport::algebraicRelativeResidual,&HydroSaddleReport::freeMomentumDefect,
        &HydroSaddleReport::maxCellContinuityDefect,&HydroSaddleReport::maxPressureMeanDefect,&HydroSaddleReport::maxBoundaryFluxDefect,
        &HydroSaddleReport::maxAbsGaugeMultiplier,&HydroSaddleReport::krylovRelativeResidual,&HydroSaddleReport::krylovScaledRelativeResidual,
        &HydroSaddleReport::velocitySpdGramBoundUpper,&HydroSaddleReport::velocitySpdScaledMassLower})
        exact(a.saddle.*member,b.saddle.*member,"F12 original stage Saddle physical gate bits");
    check(a.saddle.freeVelocityCount==b.saddle.freeVelocityCount&&a.saddle.pressureGaugeCount==b.saddle.pressureGaugeCount
        &&a.saddle.freeTangentialCount==b.saddle.freeTangentialCount&&a.saddle.linearRefinements==b.saddle.linearRefinements,
        "F12 original stage Saddle physical dimensions/refinements");
    check(a.maxwellLoadReportAvailable==b.maxwellLoadReportAvailable&&a.viscosityProductReportAvailable==b.viscosityProductReportAvailable
        &&sameClock(a.viscosityProductTime,b.viscosityProductTime)&&a.openConvection.available==b.openConvection.available,"F12 actual stage physical availability/clock");
    const auto& x=a.maxwellLoadReport;const auto& y=b.maxwellLoadReport;
    check(x.pairing==y.pairing&&x.cellSidePoints==y.cellSidePoints&&x.volumePoints==y.volumePoints
        &&x.boundaryPoints==y.boundaryPoints&&x.axisPoints==y.axisPoints,"F12 actual Maxwell physical sampling report");
}
void sameStep(const HydroVofStep& literal,const HydroVofStep& actual) {
    const auto& a=literal.report();const auto& b=actual.report();
    check(a.linearSolves==b.linearSolves&&a.usedVelocityHistory==b.usedVelocityHistory&&a.startupCorrection==b.startupCorrection
        &&a.corrected==b.corrected&&a.capillaryForcedCorrection==b.capillaryForcedCorrection&&a.maxwellCorrectionChecked==b.maxwellCorrectionChecked
        &&a.stages.size()==b.stages.size(),"F12 actual AC stages/correction decisions");
    for(auto member:{&HydroVofStepReport::deltaTime,&HydroVofStepReport::capillaryCourant,&HydroVofStepReport::interfaceMismatch,
        &HydroVofStepReport::finalInterfaceMismatch,&HydroVofStepReport::correctionRelative,&HydroVofStepReport::maxwellLoadChange,&HydroVofStepReport::maxwellWorkChange})
        exact(a.*member,b.*member,"F12 actual complete step physical report bits");
    for(std::size_t i=0;i<a.stages.size();++i)sameStage(a.stages[i],b.stages[i]);
    sameNative(literal.candidateState().values(),actual.candidateState().values());sameNative(literal.solvedMidpointValues(),actual.solvedMidpointValues());
    samePublication(literal.transport().midpoint(),actual.transport().midpoint());samePublication(literal.transport().endpoint(),actual.transport().endpoint());
    sameValues(literal.frozenMidpoint().maxwellLoads(),actual.frozenMidpoint().maxwellLoads(),"F12 actual retained Maxwell load bits");
    sameValues(literal.frozenMidpoint().bodyLoads(),actual.frozenMidpoint().bodyLoads(),"F12 actual retained body load bits");
}
HydroVofStep prepare(const QuadtreeMesh& mesh,const HydroTwoPhaseState& accepted,Clock end,const HydroVofAcConfig& config,bool open) {
    const auto call=[&]{return open?HydroVofStep::preparePrescribedFlow(mesh,accepted,end,config):HydroVofStep::prepare(mesh,accepted,end,config);};
    if(!active)return call();
#if defined(_MSC_VER) && defined(_M_X64) && defined(_OPENMP)
    struct Floating {
        std::fenv_t native{};unsigned csr=0;int error=errno;
        Floating(){csr=_mm_getcsr();check(std::fegetenv(&native)==0,"F12 capture complete native FP");}
        void restore()const{check(std::fesetenv(&native)==0,"F12 restore complete native FP");_mm_setcsr(csr);errno=error;}
        bool same(const Floating& b)const{return std::memcmp(&native,&b.native,sizeof native)==0&&csr==b.csr&&error==b.error;}
    };
    struct Runtime {Floating original;int threads=omp_get_max_threads(),dynamic=omp_get_dynamic();
        ~Runtime(){omp_set_dynamic(dynamic);omp_set_num_threads(threads);original.restore();}} runtime;
    const auto sourceU=accepted.values().normalVelocity(),sourceP=accepted.values().pressure();
    const auto source=accepted.phaseGeometry().vofPublication();const auto sourceAlpha=source.alpha();
    std::vector<std::array<long double,2>> sourceStocks;for(const auto& cell:source.plicCells())sourceStocks.push_back(cell.inventory);
    omp_set_dynamic(0);check(std::fesetenv(FE_DFL_ENV)==0,"F12 default nearest gradual entry");_mm_setcsr(0x1f80u);errno=173;const Floating entry;
    omp_set_num_threads(1);entry.restore();const auto literal=call();const Floating literalExit;
    // The real late momentum gate throws after the complete diagnostic fold.
    // Tighten only its tolerance using this same stage's nonzero physical witness.
    // Run negatives before final actual so the existing provider .back() owners
    // still denote the successful candidate consumed by the original oracle.
    if(config.maxwellLoad&&maxwellPairs==0){
        double witness=0;for(const auto& stage:literal.report().stages)if(stage.momentumPowerRelativeDefect>0){witness=stage.momentumPowerRelativeDefect;break;}
        check(std::isfinite(witness)&&witness>0,"F12 real momentum-power residual witness unavailable");
        auto bad=config;bad.momentumPowerTolerance=witness/2;
        check(bad.momentumPowerTolerance>0,"F12 real momentum-power residual witness cannot supply positive tightened gate");
        bad.maxwellLoad=[provider=config.maxwellLoad](const QuadtreeMesh& live,const GeometricVofTransport& volume,
            const HydroPhaseGeometry& phase,const HydroPressureLoad& pressure){++lateProviderCalls;return provider(live,volume,phase,pressure);};
        const auto badCall=[&]{return open?HydroVofStep::preparePrescribedFlow(mesh,accepted,end,bad):HydroVofStep::prepare(mesh,accepted,end,bad);};
        const auto preserved=[&]{
            sameValues(accepted.values().normalVelocity(),sourceU,"F12 late gate changed accepted U bits");
            sameValues(accepted.values().pressure(),sourceP,"F12 late gate changed accepted p bits");
            sameValues(source.alpha(),sourceAlpha,"F12 late gate changed accepted alpha bits");
            for(std::size_t i=0;i<sourceStocks.size();++i)sameValues(source.plicCells()[i].inventory,sourceStocks[i],"F12 late gate changed accepted independent stocks");
        };
        const auto failAtGate=[&](int threads,std::string& message,std::optional<Floating>& exit){
            omp_set_num_threads(threads);entry.restore();
            try{(void)badCall();}catch(const std::exception& error){exit.emplace();message=error.what();return;}
            check(false,"F12 tightened real late momentum gate did not reject");
        };
        std::string literalError,actualError;std::optional<Floating> literalFailure,actualFailure;
        const auto oldCalls=lateProviderCalls;failAtGate(1,literalError,literalFailure);preserved();const auto literalCalls=lateProviderCalls-oldCalls;
        failAtGate(8,actualError,actualFailure);preserved();const auto actualCalls=lateProviderCalls-oldCalls-literalCalls;
        check(literalError.find("[HydroVofStep] same-owner physical momentum power gate failed")==0&&literalError==actualError,
            "F12 original late physical momentum gate full error differs");
        check(literalFailure&&actualFailure&&literalFailure->same(*actualFailure),"F12 original late physical momentum gate native FP/raw MXCSR/errno differs");
        check(literalCalls>0&&literalCalls==actualCalls,"F12 late failure must retain each real stage provider callback");++lateFailurePairs;
    }
    omp_set_num_threads(8);entry.restore();const auto actual=call();const Floating actualExit;
    check(literalExit.same(actualExit),"F12 whole actual stage native FP/raw MXCSR/errno differs from literal");
    sameStep(literal,actual);const auto count=static_cast<std::size_t>(accepted.values().layout().velocityCount());
    for(std::size_t i=0;i<actual.report().stages.size();++i){const auto& s=literal.report().stages[i].momentumRateRows;const auto& p=actual.report().stages[i].momentumRateRows;
        check(!s.used,"F12 one-worker literal route unexpectedly used row team");
        check(p.requested&&p.eligible&&p.used&&!p.fallback&&p.rows==count&&p.workers>1&&p.workers<=8&&p.scratchBytes>0&&p.scratchBytes<=64u*1024u*1024u,
            "F12 actual stage did not consume all independently contracted rate rows with multiple workers");
        check(std::isfinite(p.seconds)&&p.seconds>=0,"F12 actual rate interval is invalid");++stages;rows+=p.rows;minimumWorkers=std::min(minimumWorkers,p.workers);maximumWorkers=std::max(maximumWorkers,p.workers);}
    check(accepted.values().normalVelocity()==sourceU&&accepted.values().pressure()==sourceP&&source.alpha()==sourceAlpha,"F12 paired preparation changed accepted U/p/alpha");
    for(std::size_t i=0;i<sourceStocks.size();++i)check(source.plicCells()[i].inventory==sourceStocks[i],"F12 paired preparation changed accepted independent phase stocks");
    ++pairs;if(config.maxwellLoad)++maxwellPairs;if(open&&!accepted.vofPrescribedFlow()->stationaryNoSlipVelocityDofs.empty())++noSlipPairs;
    return actual;
#else
    throw std::runtime_error("F12 diagnostic rate rows requires actual MSVC x64 OpenMP");
#endif
}
}
std::unique_ptr<QuadtreeMesh> meshFor(bool fine)
{
    DictionaryMesh d;d.Nx=d.Ny=4;d.maxCellNum=4096;d.maxInitialIterations=0;d.maxLevel_local={1};
    auto m=std::make_unique<QuadtreeMesh>(d);DictionaryMesh::RefineBox upper;
    upper.level=1;upper.min={0,.5};upper.max={1,1};std::vector<DictionaryMesh::RefineBox> boxes{upper};
    if(fine){auto lower=upper;lower.min={0,0};lower.max={.5,.5};boxes.push_back(lower);}
    m->applyRefineBoxes(boxes);return m;
}
NativeLeafDomain domain(const QuadtreeMesh& mesh,bool rz)
{FaceSegmentGraph graph;graph.build(mesh);return NativeLeafDomain::build(mesh,graph,rz);}
GeometricVofPublication initial(const QuadtreeMesh& mesh,bool rz)
{
    std::vector<double> a(mesh.cellIndexNum());
    for(int ci=0;ci<mesh.cellIndexNum();++ci){const auto& c=mesh.cells(mesh.cellIndexToLabel[ci]);
        const auto lo=mesh.verts(c.cellToVert[0]).coord(),hi=mesh.verts(c.cellToVert[2]).coord();
        a[ci]=std::clamp((11./32-lo.y)/(hi.y-lo.y),0.,1.);}
    return GeometricVofPublication::prepareLocalPlic(mesh,domain(mesh,rz),{0,0,0,Clock::Position::Initial},a);
}
HydroPhaseGeometry phaseFor(const QuadtreeMesh& mesh,const GeometricVofPublication& p)
{
    const auto l=HydroDofLayout::build(mesh,p.domain());
    return HydroPhaseGeometry::build(mesh,HydroPhysicalDivergence::build(mesh,HydroIntegralDivergence::build(mesh,l)),p);
}
long double psi(Vector2D p)
{const long double x=p.x,y=p.y;return 16*x*x*(1-x)*(1-x)*y*y*(1-y)*(1-y)*(1+.3L*x+.7L*y);}
HydroNativeValues valuesFor(const QuadtreeMesh& mesh,const HydroPhaseGeometry& phase)
{
    const auto& l=phase.layout();const bool rz=l.domain().axisymmetric();std::vector<double> u(l.velocityCount()),p(l.pressureCount());
    for(int i=0;i<l.velocityCount();++i){const auto f=l.domain().face(mesh,l.nativeFace(i));
        const auto& e=f.horizontal?mesh.hEdges(f.edgeLabel):mesh.vEdges(f.edgeLabel);
        auto a=mesh.verts(e.edgeToVert[0]).coord(),b=mesh.verts(e.edgeToVert[1]).coord();
        if(f.horizontal?a.x>b.x:a.y>b.y)std::swap(a,b);
        const auto flux=(f.horizontal?-1.L:1.L)*(rz?CartesianMeasure::twoPi:1.L)*(psi(b)-psi(a));
        u[i]=f.area>0?static_cast<double>(flux/f.area):0;}
    for(int i=0;i<l.pressureCount();++i)p[i]=mesh.cells(mesh.cellIndexToLabel[l.nativeCell(i)]).centerCoord().y-.5;
    return {l,u,p};
}
HydroMomentumMass massFor(const QuadtreeMesh& mesh,const HydroTwoPhaseState& state,Density r)
{return r==Density::CellAverageProjection
    ?HydroMomentumMass::buildCellAverageProjection(mesh,state.momentumGeometry(),state.phaseGeometry(),10,1)
    :HydroMomentumMass::build(mesh,state.momentumGeometry(),state.phaseGeometry(),10,1);}
HydroVofRemapNativeBudget budget(const HydroMomentumMass& mass,const HydroNativeValues& v)
{
    HydroVofRemapNativeBudget b;for(int u=0;u<mass.layout().velocityCount();++u){
        const int c=mass.layout().normalAxis(u)==HydroDofLayout::NormalAxis::First?0:1;
        const long double m=mass.masses()[u],x=v.normalVelocity()[u];b.mass[c]+=m;b.momentum[c]+=m*x;b.kinetic+=.5L*m*x*x;}return b;
}
void sameBudget(const HydroVofRemapNativeBudget& actual,const HydroVofRemapNativeBudget& expected)
{
    for(int c=0;c<2;++c){near(actual.mass[c],expected.mass[c],expected.mass[c],"native mass ledger");
        near(actual.momentum[c],expected.momentum[c],std::abs(expected.momentum[c]),"native momentum ledger");}
    near(actual.kinetic,expected.kinetic,expected.kinetic,"native kinetic ledger");
}
std::array<long double,2> stocks(const GeometricVofPublication& p)
{std::array<long double,2> v{};for(const auto& c:p.plicCells())for(int k=0;k<2;++k)v[k]+=c.inventory[k];return v;}
void verify(const QuadtreeMesh& a,const HydroTwoPhaseState& old,const QuadtreeMesh& b,
    const HydroVofMeshTransfer& t,Density representation)
{
    ++remaps;const auto& next=t.candidateState();const auto& r=t.report();const auto& l=next.values().layout();next.validateFor(b);
    check(t.sourceState().sameIdentity(old),"authentic source handle");
    check(!old.values().layout().domain().uniformLeaves()&&!l.domain().uniformLeaves(),"actual hanging 2:1 leaves");
    check(sameClock(next.timeLayer(),old.timeLayer()),"same complete accepted End clock");
    check(next.phaseGeometry().vofIdentity()==t.interfaceTransfer().candidate().identity(),"one actual VOF target owner");
    check(r.remapSolves==1&&r.historyReset&&!next.velocityHistory(),"one remap solve and explicit history reset");
    check(r.polynomialIntegrals>0&&r.polynomialIntegrals==r.faceSpans,"actual source polynomial flux spans");
    check(next.retainedMass()&&next.retainedMass()->representation()==representation,"retained same selected mass representation");
    check(next.retainedMass()->geometry().sameIdentity(next.momentumGeometry()),"target mass/geometry reused");
    check(t.projectionMultiplier().size()==next.values().pressure().size()
        &&t.projectionMultiplier()!=next.values().pressure(),"remap multiplier is not physical pressure");
    const auto& overlap=t.interfaceTransfer().nativeOverlaps();long double pressureMean=0,pressureGross=0;
    for(int p=0;p<l.pressureCount();++p){long double integral=0;
        for(auto id:overlap.byTarget()[l.nativeCell(p)]){const auto& o=overlap.pieces()[id];
            integral+=o.capacity*old.values().pressure()[old.values().layout().pressureDof(o.sourceCell)];}
        const auto expected=static_cast<double>(static_cast<double>(integral/next.phaseGeometry().metrics().cellVolumes()[p])-r.pressureGaugeShift);
        near(next.values().pressure()[p],expected,1,"physical overlap-average pressure");
        const auto x=static_cast<long double>(next.phaseGeometry().metrics().cellVolumes()[p])*next.values().pressure()[p];
        pressureMean+=x;pressureGross+=std::abs(x);}
    near(pressureMean,0,pressureGross,"complete physical pressure volume gauge");
    const auto flow=next.phaseGeometry().metrics().divergence().apply(b,next.values());
    const auto& d=next.phaseGeometry().metrics().divergence().matrix();
    for(int p=0;p<l.pressureCount();++p){double scale=0;for(HydroIntegralDivergence::Matrix::InnerIterator it(d,p);it;++it)
        scale+=std::abs(it.value()*next.values().normalVelocity()[it.col()]);
        check(scale>0?std::abs(flow[p])/scale<=1e-11:flow[p]==0,"original native divergence gate");}
    for(int u:next.boundary().velocityDofs())check(next.values().normalVelocity()[u]==0&&t.rawValues().normalVelocity()[u]==0,"closed native wall mean exactly zero");
    const auto ma=massFor(a,old,representation),mb=massFor(b,next,representation);
    const auto sa=budget(ma,old.values()),sr=budget(mb,t.rawValues()),sb=budget(mb,next.values());
    sameBudget(r.source,sa);sameBudget(r.rawTarget,sr);sameBudget(r.target,sb);
    for(int c=0;c<2;++c){near(r.momentumChange[c],sb.momentum[c]-sa.momentum[c],1,"reported full DeltaP");
        near(r.projectionMomentumChange[c],sb.momentum[c]-sr.momentum[c],1,"reported projection DeltaP");}
    near(r.kineticChange,sb.kinetic-sa.kinetic,1,"reported full DeltaKE");
    near(r.projectionKineticChange,sb.kinetic-sr.kinetic,1,"reported projection DeltaKE");
    long double correction=0;for(int u=0;u<l.velocityCount();++u){const long double delta=next.values().normalVelocity()[u]-static_cast<long double>(t.rawValues().normalVelocity()[u]);
        correction+=.5L*mb.masses()[u]*delta*delta;}
    near(r.projectionCorrectionEnergy,correction,sr.kinetic,"actual M correction norm");
    near(r.projectionEnergyIdentityResidual,(sb.kinetic-sr.kinetic)+correction,sr.kinetic,"reported M-orthogonal identity");
    near((sb.kinetic-sr.kinetic)+correction,0,sr.kinetic,"independent M-orthogonal energy identity");
    check(sb.kinetic<=sr.kinetic+2e-11L*std::max(1.L,sr.kinetic),"projection does not invent KE");
    const auto oldStock=stocks(old.phaseGeometry().vofPublication()),newStock=stocks(next.phaseGeometry().vofPublication());
    for(int k=0;k<2;++k)near(oldStock[k],newStock[k],oldStock[k],"independent liquid/gas stocks");
    if(l.domain().axisymmetric()){
        check(r.axisNormals>0,"true zero-area RZ axis");unsigned matched=0;
        // Actual unchanged native-face means must remain cardinal after the
        // regular RZ integral: catches missing OR duplicated radial weight.
        const auto& ol=old.values().layout();for(int u=0;u<l.velocityCount();++u){const auto f=l.domain().face(b,l.nativeFace(u));if(f.area==0)continue;
            for(int v=0;v<ol.velocityCount();++v){const auto q=ol.domain().face(a,ol.nativeFace(v));
                if(f.horizontal==q.horizontal&&f.center.x==q.center.x&&f.center.y==q.center.y&&f.area==q.area){
                    near(t.rawValues().normalVelocity()[u],old.values().normalVelocity()[v],1,"actual RZ cardinal physical flux");++matched;break;}}}
        check(matched>0,"RZ cardinal face supports checked");}
}
void scenario(bool rz,Density representation)
{
    auto coarse=meshFor(false),fine=meshFor(true);const auto start=initial(*coarse,rz);const auto p0=phaseFor(*coarse,start);
    const auto v=valuesFor(*coarse,p0);const auto first=HydroTwoPhaseState::createVof(*coarse,p0,v);
    std::vector<std::array<long double,2>> inventory(start.alpha().size());for(std::size_t ci=0;ci<start.plicCells().size();++ci)inventory[ci]=start.plicCells()[ci].inventory;
    const auto end=GeometricVofPublication::prepareLocalPlicInventory(*coarse,start.domain(),{1,0,.125,Clock::Position::End},inventory,{});
    const auto phase=HydroPhaseGeometry::build(*coarse,p0.metrics(),end);
    const auto retained=representation==Density::CellAverageProjection
        ?HydroMomentumMass::buildCellAverageProjection(*coarse,first.momentumGeometry(),phase,10,1)
        :HydroMomentumMass::build(*coarse,first.momentumGeometry(),phase,10,1);
    const auto old=HydroTwoPhaseState::createVofAccepted(*coarse,phase,v,first,retained);
    const auto oldU=old.values().normalVelocity(),oldP=old.values().pressure(),oldAlpha=end.alpha();const auto oldStock=stocks(end);
    HydroVofMeshTransferConfig cfg;cfg.momentumDensityRepresentation=representation;
    const auto transfer=GeometricVofMeshTransfer::prepare(*coarse,end,*fine,domain(*fine,rz),old.timeLayer());
    check(transfer.report().refinedTargets>0,"real refinement coverage");
    const auto refined=HydroVofMeshTransfer::prepare(*coarse,old,*fine,transfer,10,1,cfg);verify(*coarse,old,*fine,refined,representation);
    const auto back=GeometricVofMeshTransfer::prepare(*fine,refined.candidateState().phaseGeometry().vofPublication(),*coarse,start.domain(),old.timeLayer());
    check(back.report().coarsenedTargets>0,"real coarsening coverage");
    const auto coarsened=HydroVofMeshTransfer::prepare(*fine,refined.candidateState(),*coarse,back,10,1,cfg);
    verify(*fine,refined.candidateState(),*coarse,coarsened,representation);
    if(!rz){HydroVofAcConfig stepCfg;stepCfg.liquidDensity=10;stepCfg.gasDensity=1;stepCfg.momentumDensityRepresentation=representation;
        stepCfg.scheme=HydroTwoPhaseTimeScheme::ExtrapolatedMidpoint;const auto begin=refined.candidateState().timeLayer();
        const auto step=HydroVofStep::prepare(*fine,refined.candidateState(),{begin.step+1,begin.time(),begin.time()+1e-5,Clock::Position::End},stepCfg);
        step.candidateState().validateFor(*fine);check(!step.report().usedVelocityHistory&&step.candidateState().velocityHistory(),"actual VOF step after explicit history reset");++steps;}
    if(!rz){const auto foreign=GeometricVofPublication::prepareLocalPlicInventory(*coarse,start.domain(),old.timeLayer(),inventory,{});
        const auto foreignTransfer=GeometricVofMeshTransfer::prepare(*coarse,foreign,*fine,domain(*fine,rz),old.timeLayer());
        rejects([&]{HydroVofMeshTransfer::prepare(*coarse,old,*fine,foreignTransfer,10,1,cfg);},"reject foreign actual VOF source identity");}
    rejects([&]{HydroVofMeshTransfer::prepare(*coarse,first,*fine,transfer,10,1,cfg);},"reject Initial remap source");
    auto bad=cfg;bad.momentumDensityRepresentation=representation==Density::SharpGeometry?Density::CellAverageProjection:Density::SharpGeometry;
    rejects([&]{HydroVofMeshTransfer::prepare(*coarse,old,*fine,transfer,10,1,bad);},"reject changing actual accepted M representation");
    check(old.values().normalVelocity()==oldU&&old.values().pressure()==oldP&&end.alpha()==oldAlpha&&stocks(end)==oldStock
        &&old.velocityHistory()&&old.phaseGeometry().vofIdentity()==end.identity(),"all remaps preserve actual source state and history");
}

namespace MaxwellCoupledContract {
namespace E=LAMP::Electric;
void bounded(long double error,long double scale,const char* why,long double gate=1e-12L) {
    check(std::isfinite(error)&&std::isfinite(scale)&&scale>=0,why);check(scale>0?std::abs(error)<=gate*scale:error==0,why);
}
std::vector<std::array<long double,8>> sourceSnapshot(const GeometricVofPublication& p) {
    std::vector<std::array<long double,8>> result;
    for(const auto& c:p.plicCells())result.push_back({c.inventory[0],c.inventory[1],c.measured[0],c.measured[1],c.normal.x,c.normal.y,c.offset,c.resolved?1.L:0.L});
    return result;
}
void run() {
    auto mesh=meshFor(false);GeometricVofInitializationConfig initialization;initialization.subcellsPerAxis=8;
    const auto pub=GeometricVofInitialization::prepare(*mesh,false,[](Vector2D p){const double x=p.x-.47,y=p.y-.48;return x*x+y*y-.23*.23;},initialization);
    const auto phase=phaseFor(*mesh,pub);const auto& layout=phase.layout();check(!layout.domain().uniformLeaves(),"coupled Maxwell fixture lacks real 2:1 AMR");
    const HydroNativeValues zero(layout,std::vector<double>(layout.velocityCount()),std::vector<double>(layout.pressureCount()));
    const auto accepted=HydroTwoPhaseState::createVof(*mesh,phase,zero);const auto oldState=accepted;
    const auto sourceQ=E::ElectricVofCharge::fromIntegrated(*mesh,pub,std::vector<double>(mesh->cellIndexNum()));const auto oldQ=sourceQ;
    const auto oldU=accepted.values().normalVelocity(),oldP=accepted.values().pressure(),oldCharge=sourceQ.integratedCharge();
    const auto oldGeometry=sourceSnapshot(pub);const auto oldStock=stocks(pub);const auto oldClock=accepted.timeLayer();
    const Clock end{oldClock.step+1,oldClock.time(),oldClock.time()+1e-5,Clock::Position::End};
    const Clock middle{end.step,end.begin,end.end,Clock::Position::Midpoint};
    HydroVofAcConfig config;config.liquidDensity=10;config.gasDensity=1;config.liquidViscosity=config.gasViscosity=0;
    config.surfaceTension=0;config.momentumDensityRepresentation=Density::CellAverageProjection;
    config.scheme=HydroTwoPhaseTimeScheme::AdaptiveLinearMidpoint;
    const E::ElectricVofMaterials materials{2,5,0,0};std::size_t callbacks=0;
    std::vector<GeometricVofTransport> volumes;std::vector<E::ElectricVofChargeTransport> advected;
    std::vector<E::ElectricVofField> fields;std::vector<HydroMaxwellLoad> loads;
    std::vector<const void*> pressureOwners,phaseOwners;volumes.reserve(2);advected.reserve(2);fields.reserve(2);loads.reserve(2);
    config.maxwellLoad=[&](const QuadtreeMesh& live,const GeometricVofTransport& volume,
        const HydroPhaseGeometry& actualMid,const HydroPressureLoad& actualPressure)->HydroMaxwellLoad {
        ++callbacks;check(&live==mesh.get(),"coupled provider received another live mesh");
        check(volume.source().identity()==pub.identity()&&actualMid.vofIdentity()==volume.midpoint().identity()&&
            sameClock(actualMid.timeLayer(),middle)&&sameClock(volume.midpoint().timeLayer(),middle),"coupled provider did not receive actual coeval stage VOF");
        check(&actualPressure.divergence().matrix()==&actualMid.metrics().divergence().matrix(),"coupled provider pressure is not actual stage C owner");
        const auto moved=E::ElectricVofChargeTransport::prepareClosed(live,sourceQ,volume);
        check(moved.transport().sameInstance(volume)&&moved.sourceCharge().sameInstance(sourceQ)&&
            moved.midpoint().geometry().identity()==volume.midpoint().identity()&&sameClock(moved.midpoint().timeLayer(),middle),"coupled actual midpoint Qadv receipt lost");
        std::vector<E::ElectricVofBoundary> boundary(actualMid.vofPublication().domain().nativeFaceCount());
        for(int row=0;row<static_cast<int>(boundary.size());++row){const auto f=actualMid.vofPublication().domain().face(live,row);
            boundary[row]=f.negativeCell>=0&&f.positiveCell>=0?E::ElectricVofBoundary{E::ElectricVofBoundary::Type::Interior,0}:
                E::ElectricVofBoundary{E::ElectricVofBoundary::Type::Dirichlet,.01*f.center.x};}
        const auto field=E::ElectricVofField::prepareOhmic(live,actualMid.vofPublication(),moved.midpoint(),materials,boundary);
        check(field.advection()&&field.advection()->sameInstance(moved.midpoint())&&field.sourceCharge().sameInstance(sourceQ)&&
            field.geometry().identity()==actualMid.vofIdentity()&&sameClock(field.geometryTime(),middle)&&sameClock(field.candidateTime(),middle)&&
            sameClock(field.chargeSourceTime(),oldClock),"coupled field fabricated a Begin/geometry/candidate clock");
        check(field.chargeCandidate().integratedCharge()==oldCharge&&field.report().relativeLinearResidual<=1e-11&&
            field.report().maximumRelativeGaussResidual<=1e-10,"coupled actual Ohmic solve changed retained Q/linear/Gauss gates");
        const auto stress=E::ElectricVofMaxwellStress::prepare(live,field);
        const auto load=HydroMaxwellLoad::prepare(live,actualPressure,actualMid,stress);
        load.validateFor(live,actualPressure,actualMid,stress,middle);
        volumes.push_back(volume);advected.push_back(moved);fields.push_back(field);loads.push_back(load);
        pressureOwners.push_back(&actualPressure.matrix());phaseOwners.push_back(&actualMid.cellVolumes());return load;
    };
    const auto step=F12DiagnosticRateRows::prepare(*mesh,accepted,end,config,false);step.candidateState().validateFor(*mesh);
    const auto& report=step.report();const auto& frozen=step.frozenMidpoint();const auto& dual=frozen.material();
    check(callbacks==(F12DiagnosticRateRows::active?2u:1u)*report.stages.size()+F12DiagnosticRateRows::lateProviderCalls&&callbacks==(F12DiagnosticRateRows::active?2u:1u)*report.linearSolves+F12DiagnosticRateRows::lateProviderCalls&&report.stages.size()==2&&report.startupCorrection&&report.corrected,
        "coupled provider must run once for each real startup AC stage");
    check(frozen.hasMaxwellLoad()&&frozen.forces().maxwellLoad&&!frozen.hasInterfaceLoad()&&!step.hasCapillaryLoad(),"coupled NS did not retain actual standalone Maxwell owner");
    check(step.transport().sameInstance(volumes.back())&&dual.material().vofTransport().sameInstance(step.transport())&&
        fields.back().geometry().identity()==step.transport().midpoint().identity()&&
        frozen.forces().maxwellLoad->phaseGeometry().vofIdentity()==step.transport().midpoint().identity()&&
        &frozen.pressureLoad().matrix()==pressureOwners.back()&&&dual.material().midpointGeometry().cellVolumes()==phaseOwners.back()&&
        &frozen.forces().maxwellLoad->stress().cells()==&loads.back().stress().cells(),"coupled final NS/field/load owners do not match actual final stage");
    check(frozen.maxwellLoads()==loads.back().momentumRhsLoad()&&frozen.viscosityMatrix().nonZeros()==0&&
        frozen.viscosity().liquidViscosity()==0&&frozen.viscosity().gasViscosity()==0,"coupled Maxwell replaced real load or inserted viscosity");
    long double forceL1=0,maxResponse=0,maxPressure=0,maxR=0,maxMomentumDefect=0,maxClosure=0;
    for(int u:frozen.saddle().freeVelocityDofs())forceL1+=std::abs(frozen.maxwellLoads()[u]);
    for(std::size_t u=0;u<oldU.size();++u)maxResponse=std::max(maxResponse,std::abs(static_cast<long double>(step.candidateState().values().normalVelocity()[u])-oldU[u]));
    for(double p:step.candidateState().values().pressure())maxPressure=std::max(maxPressure,std::abs(static_cast<long double>(p)));
    check(std::isfinite(forceL1)&&forceL1>0&&std::isfinite(maxResponse)&&maxResponse>0&&maxPressure>0,"coupled Maxwell did not produce a real free native force/NS response");
    for(const auto& stage:report.stages){check(stage.saddle.algebraicRelativeResidual<=1e-11&&stage.saddle.freeMomentumDefect<=1e-11&&
            stage.saddle.maxCellContinuityDefect<=1e-11&&stage.saddle.maxPressureMeanDefect<=1e-11&&stage.saddle.maxBoundaryFluxDefect<=1e-12,
            "coupled native saddle original gates failed");
        check(stage.momentumPowerRelativeDefect<=config.momentumPowerTolerance,
            "coupled retained physical momentum/power gate failed");
        check(std::isfinite(stage.energy.closureRelativeResidual)&&stage.energy.closureRelativeResidual<=1e-11L&&stage.energy.viscousWork==0&&stage.energy.bodyWork==0&&stage.capillaryWork==0,
            "coupled actual KE/Maxwell budget or exact-zero non-electric work failed");
        maxMomentumDefect=std::max(maxMomentumDefect,static_cast<long double>(stage.momentumPowerRelativeDefect));
        maxClosure=std::max(maxClosure,stage.energy.closureRelativeResidual);}
    LAMP::Geometry::Detail::FiniteVolume::Dyadic actualPower;
    const auto& solved=step.solvedMidpointValues().normalVelocity();const auto& rhs=frozen.maxwellLoads();
    for(std::size_t i=0;i<rhs.size();++i)actualPower.add(LAMP::Geometry::Detail::FiniteVolume::Dyadic::product({solved[i],rhs[i]}));
    const long double expectedWork=report.deltaTime*actualPower.rounded();
    bounded(report.stages.back().maxwellWork-expectedWork,std::abs(report.stages.back().maxwellWork)+std::abs(expectedWork),"coupled final Maxwell work is not actual dt*U_mid.F");
    check(report.stages.back().maxwellWork>0,"coupled electric work did not drive the actual accepted-candidate KE");
    dual.validateFor(*mesh,dual.sources(),middle);
    for(const auto& c:dual.cells()){const long double residual=c.integratedDefect-c.massChange-c.outwardIntegral;
        bounded(residual,std::abs(c.integratedDefect)+std::abs(c.massChange)+std::abs(c.outwardIntegral),"coupled actual canonical R identity");
        maxR=std::max(maxR,std::abs(c.integratedDefect));}
    const auto newStock=stocks(step.candidateState().phaseGeometry().vofPublication());
    for(int k=0;k<2;++k)bounded(newStock[k]-oldStock[k],std::abs(oldStock[k])+std::abs(newStock[k]),"coupled independently retained closed phase conservation");
    check(sameClock(step.candidateState().timeLayer(),end)&&step.candidateState().velocityHistory()&&
        sourceQ.sameInstance(oldQ)&&sameClock(sourceQ.timeLayer(),oldClock)&&sourceQ.integratedCharge()==oldCharge,
        "coupled candidate/reusable accepted Q clocks were relabeled");
    const auto preserved=[&]{check(accepted.sameIdentity(oldState)&&sameClock(accepted.timeLayer(),oldClock)&&
            accepted.values().normalVelocity()==oldU&&accepted.values().pressure()==oldP&&!accepted.velocityHistory()&&
            sourceSnapshot(accepted.phaseGeometry().vofPublication())==oldGeometry&&sourceQ.sameInstance(oldQ)&&sourceQ.integratedCharge()==oldCharge&&
            sameClock(sourceQ.timeLayer(),oldClock),"coupled preparation/rejection changed accepted VOF/U/p/Q without commit");};
    preserved();
    auto throwing=config;std::size_t throwCalls=0;throwing.maxwellLoad=[&](const QuadtreeMesh&,const GeometricVofTransport&,
        const HydroPhaseGeometry&,const HydroPressureLoad&)->HydroMaxwellLoad {++throwCalls;throw std::runtime_error("COUPLED_PROVIDER_EXPECTED_THROW");};
    bool caught=false;try{HydroVofStep::prepare(*mesh,accepted,end,throwing);}catch(const std::exception& e){caught=std::string(e.what()).find("COUPLED_PROVIDER_EXPECTED_THROW")!=std::string::npos;}
    check(caught&&throwCalls==1,"coupled throwing provider was not the real rejection source");preserved();
    auto stale=config;std::size_t staleCalls=0;stale.maxwellLoad=[&](const QuadtreeMesh&,const GeometricVofTransport& volume,
        const HydroPhaseGeometry& actualMid,const HydroPressureLoad&)->HydroMaxwellLoad {++staleCalls;
        check(volume.midpoint().identity()==actualMid.vofIdentity()&&actualMid.vofIdentity()!=loads.front().phaseGeometry().vofIdentity()&&
            sameClock(actualMid.timeLayer(),loads.front().phaseGeometry().timeLayer()),"stale-load negative needs a genuine new owner at same fullclock");return loads.front();};
    caught=false;try{HydroVofStep::prepare(*mesh,accepted,end,stale);}catch(const std::exception& e){caught=std::string(e.what()).find("[HydroMaxwellLoad]")!=std::string::npos;}
    check(caught&&staleCalls==1,"coupled foreign same-clock Maxwell load was not rejected by actual typed consumer");preserved();
    std::cout<<std::setprecision(17)<<"MAXWELL_COUPLED_PASS,cells,"<<mesh->cellIndexNum()<<",AC_stages,"<<report.stages.size()<<",provider_calls,"<<callbacks
        <<",free_force_L1,"<<forceL1<<",NS_max_deltaU,"<<maxResponse<<",max_pressure,"<<maxPressure<<",Maxwell_work,"<<report.stages.back().maxwellWork
        <<",momentum_power,"<<maxMomentumDefect<<",KE_closure,"<<maxClosure<<",actual_R,"<<maxR
        <<",throw_rejection_calls,"<<throwCalls<<",stale_rejection_calls,"<<staleCalls<<",mu0,exact0,accepted_VOF_U_p_Q_preserved,no_App_commit\n";
}
} // namespace MaxwellCoupledContract

// BEGIN F23 ACTUAL MAXWELL GEOMETRY REUSE CONTRACT
namespace F23MaxwellGeometryReuseContract {
namespace E=LAMP::Electric;
using Timer=std::chrono::steady_clock;
unsigned scenes=0,loadsCompared=0,nativeValuesCompared=0,workValuesCompared=0,guards=0;
template<class Number> void exact(Number a,Number b,const char* why) {
    check(std::isfinite(a)&&std::isfinite(b)&&a==b&&(a!=0||std::signbit(a)==std::signbit(b)),why);
}
void sameValues(const std::vector<double>& a,const std::vector<double>& b,const char* why) {
    check(a.size()==b.size(),why);for(std::size_t k=0;k<a.size();++k){exact(a[k],b[k],why);++nativeValuesCompared;}
}
template<class F> void rejectsWith(F&& f,const char* token,const char* why) {
    bool caught=false;try{f();}catch(const std::exception& e){caught=std::string(e.what()).find(token)!=std::string::npos;}
    check(caught,why);++guards;
}
std::vector<E::ElectricVofBoundary> electricalBoundary(const QuadtreeMesh& mesh,const GeometricVofPublication& pub) {
    std::vector<E::ElectricVofBoundary> result(pub.domain().nativeFaceCount());
    for(int row=0;row<static_cast<int>(result.size());++row){const auto f=pub.domain().face(mesh,row);
        result[row]=f.negativeCell>=0&&f.positiveCell>=0?E::ElectricVofBoundary{E::ElectricVofBoundary::Type::Interior,0}:
            f.area==0?E::ElectricVofBoundary{E::ElectricVofBoundary::Type::HomogeneousNeumann,0}:
            E::ElectricVofBoundary{E::ElectricVofBoundary::Type::Dirichlet,.01*f.center.y};}
    return result;
}
void compare(const QuadtreeMesh& mesh,const HydroMaxwellLoad& fresh,const HydroMaxwellLoad& reused,
    const HydroStrainOperator& geometry,const HydroNativeValues& velocity) {
    check(&fresh.geometry().matrix()!=&geometry.matrix()&&&reused.geometry().matrix()==&geometry.matrix(),
        "F23 reused geometry must retain the prepared handle while fresh really builds another geometry");
    check(&reused.geometry().stressFaceVelocityMatrix()==&geometry.stressFaceVelocityMatrix()&&
        &reused.geometry().stressVolumeVelocityMatrix()==&geometry.stressVolumeVelocityMatrix()&&
        &reused.geometry().stressComponentIntegralMatrix()==&geometry.stressComponentIntegralMatrix(),
        "F23 all actual prepared geometry matrices must retain the same owner");
    sameValues(fresh.momentumRhsLoad(),reused.momentumRhsLoad(),"F23 exact native Maxwell momentum RHS differs");
    const auto& a=fresh.report();const auto& b=reused.report();
    check(a.volumePoints==b.volumePoints&&a.boundaryPoints==b.boundaryPoints&&a.axisPoints==b.axisPoints,
        "F23 Maxwell load report quadrature counts differ");
    const auto power=[&](const HydroMaxwellLoad& load){LAMP::Geometry::Detail::FiniteVolume::Dyadic sum;
        const auto& u=velocity.normalVelocity();for(std::size_t k=0;k<u.size();++k)
            sum.add(LAMP::Geometry::Detail::FiniteVolume::Dyadic::product({u[k],load.momentumRhsLoad()[k]}));return sum.rounded();};
    exact(power(fresh),power(reused),"F23 actual native RHS exact-dot power differs");++workValuesCompared;
    check(std::any_of(velocity.normalVelocity().begin(),velocity.normalVelocity().end(),[](double u){return u!=0;})
        &&std::any_of(fresh.stress().field().faces().begin(),fresh.stress().field().faces().end(),[](const auto& f){return f.integratedField!=0;}),
        "F23 comparison requires actual nonzero field and native velocity, not a nonzero net force");++loadsCompared;
}
void scene(bool rz) {
    ++scenes;auto mesh=meshFor(false);GeometricVofInitializationConfig init;init.subcellsPerAxis=8;
    auto pub=GeometricVofInitialization::prepare(*mesh,rz,[](Vector2D p){const double x=p.x-.47,y=p.y-.48;return x*x+y*y-.23*.23;},init);
    const auto first=phaseFor(*mesh,pub);const auto metrics=first.metrics();
    check(!first.layout().domain().uniformLeaves(),"F23 fixture needs actual hanging 2:1 leaves");
    const auto geometry=HydroStrainOperator::buildStressCompatible(*mesh,metrics);
    check(geometry.isStressCompatible()&&geometry.layout().sameIdentity(first.layout())&&
        &geometry.metrics().divergence().matrix()==&metrics.divergence().matrix(),"F23 prepared geometry actual native D owner");
    const auto pressure=HydroPressureLoad::build(metrics.divergence());
    auto charge=E::ElectricVofCharge::fromIntegrated(*mesh,pub,std::vector<double>(mesh->cellIndexNum()));
    const auto initialPub=pub;const auto initialCharge=charge;const auto savedGeometry=MaxwellCoupledContract::sourceSnapshot(pub);
    const auto savedQ=charge.integratedCharge();const E::ElectricVofMaterials materials{2,5,0,0};
    const auto initialField=E::ElectricVofField::prepareElectrostaticSnapshot(*mesh,pub,charge,materials,electricalBoundary(*mesh,pub));
    const auto initialStress=E::ElectricVofMaxwellStress::prepare(*mesh,initialField);
    rejectsWith([&]{HydroMaxwellLoad::prepare(*mesh,pressure,first,initialStress);},"requires a complete Midpoint clock","F23 old qualification must reject Initial stress");
    rejectsWith([&]{HydroMaxwellLoad::prepare(*mesh,pressure,first,initialStress,geometry);},"requires a complete Midpoint clock","F23 reuse must preserve original Initial/full-clock qualification");
    auto changedMesh=meshFor(true);const auto changedPub=initial(*changedMesh,rz);const auto changedPhase=phaseFor(*changedMesh,changedPub);
    check(changedMesh->cellIndexNum()!=mesh->cellIndexNum(),"F23 changed-mesh guard needs real refinement");
    const auto changedGeometry=HydroStrainOperator::buildStressCompatible(*changedMesh,changedPhase.metrics());
    const auto raw=HydroStrainOperator::build(*mesh,metrics,HydroStrainOperator::volumeQuadrature(*mesh,metrics));
    check(!raw.isStressCompatible(),"F23 raw strain negative must be genuinely non-compatible");
    const auto foreignMetrics=HydroPhysicalDivergence::build(*mesh,HydroIntegralDivergence::build(*mesh,first.layout()));
    const auto foreignD=HydroStrainOperator::buildStressCompatible(*mesh,foreignMetrics);
    check(foreignD.layout().sameIdentity(first.layout())&&&foreignD.metrics().divergence().matrix()!=&metrics.divergence().matrix(),
        "F23 foreign-D negative needs same actual layout but a rebuilt divergence owner");
    const auto foreignLayout=HydroDofLayout::build(*mesh,pub.domain());
    const auto foreignLayoutMetrics=HydroPhysicalDivergence::build(*mesh,HydroIntegralDivergence::build(*mesh,foreignLayout));
    const auto foreignLayoutGeometry=HydroStrainOperator::buildStressCompatible(*mesh,foreignLayoutMetrics);
    check(!foreignLayout.sameIdentity(first.layout()),"F23 foreign-layout negative must have another actual owner");
    long double freshSeconds=0,reuseSeconds=0,geometrySeconds=0;auto previousMid=initialPub.identity();
    for(int stage=0;stage<2;++stage){const auto old=pub.timeLayer();const Clock end{old.step+1,old.time(),old.time()+1e-5,Clock::Position::End};
        const std::vector<double> zero(pub.domain().nativeFaceCount());
        const auto slab=GeometricVofAdvectorSlab::prepare(*mesh,pub,end,zero,zero);
        const auto volume=GeometricVofTransport::prepareLocalPlicSweep(*mesh,pub,slab,10,1);
        const auto moved=E::ElectricVofChargeTransport::prepareClosed(*mesh,charge,volume);
        const auto phase=HydroPhaseGeometry::build(*mesh,metrics,volume.midpoint());
        const auto field=E::ElectricVofField::prepareOhmic(*mesh,volume.midpoint(),moved.midpoint(),materials,electricalBoundary(*mesh,volume.midpoint()));
        const auto stress=E::ElectricVofMaxwellStress::prepare(*mesh,field);
        const auto velocity=valuesFor(*mesh,phase);const auto middle=phase.timeLayer();
        check(sameClock(field.geometryTime(),middle)&&sameClock(field.candidateTime(),middle)&&
            phase.vofIdentity()==volume.midpoint().identity()&&field.geometry().identity()==phase.vofIdentity()&&
            field.advection()&&field.advection()->sameInstance(moved.midpoint()),"F23 field/stress must use actual coeval VOF/charge transport receipt");
        if(stage)check(previousMid!=phase.vofIdentity(),"F23 second phase must retain a distinct real midpoint owner");
        previousMid=phase.vofIdentity();
        const auto fresh=HydroMaxwellLoad::prepare(*mesh,pressure,phase,stress);
        const auto reused=HydroMaxwellLoad::prepare(*mesh,pressure,phase,stress,geometry);
        fresh.validateFor(*mesh,pressure,phase,stress,middle);reused.validateFor(*mesh,pressure,phase,stress,middle);
        compare(*mesh,fresh,reused,geometry,velocity);
        rejectsWith([&]{HydroMaxwellLoad::prepare(*mesh,pressure,phase,stress,raw);},"prepared Maxwell geometry requires actual RT0 Green capability","F23 raw prepared strain must fail");
        rejectsWith([&]{HydroMaxwellLoad::prepare(*mesh,pressure,phase,stress,foreignD);},"prepared Maxwell geometry must retain the actual phase divergence owner","F23 rebuilt same-layout D must fail");
        rejectsWith([&]{HydroMaxwellLoad::prepare(*mesh,pressure,phase,stress,foreignLayoutGeometry);},"prepared Maxwell geometry has a different actual native layout","F23 independently rebuilt same-mesh layout must fail");
        rejectsWith([&]{HydroMaxwellLoad::prepare(*mesh,pressure,phase,stress,changedGeometry);},"","F23 genuinely changed mesh geometry must fail");
        const auto badPressure=HydroPressureLoad::build(foreignMetrics.divergence());
        rejectsWith([&]{HydroMaxwellLoad::prepare(*mesh,badPressure,phase,stress);},"pressure must retain the actual phase divergence owner","F23 original qualified pressure owner must fail");
        rejectsWith([&]{HydroMaxwellLoad::prepare(*mesh,badPressure,phase,stress,geometry);},"pressure must retain the actual phase divergence owner","F23 reused path must retain original pressure owner guard");
        if(stage==0)for(int sample=0;sample<3;++sample){
            auto begin=Timer::now();const auto built=HydroStrainOperator::buildStressCompatible(*mesh,metrics);
            geometrySeconds+=std::chrono::duration<long double>(Timer::now()-begin).count();check(built.isStressCompatible(),"F23 timed genuine geometry build");
            begin=Timer::now();const auto a=HydroMaxwellLoad::prepare(*mesh,pressure,phase,stress);
            freshSeconds+=std::chrono::duration<long double>(Timer::now()-begin).count();
            begin=Timer::now();const auto b=HydroMaxwellLoad::prepare(*mesh,pressure,phase,stress,geometry);
            reuseSeconds+=std::chrono::duration<long double>(Timer::now()-begin).count();compare(*mesh,a,b,geometry,velocity);
        }
        const auto endField=E::ElectricVofField::prepareOhmic(*mesh,volume.endpoint(),moved.endpoint(),materials,electricalBoundary(*mesh,volume.endpoint()));
        pub=volume.endpoint();charge=endField.chargeCandidate();
    }
    check(initialCharge.sameInstance(initialField.sourceCharge())&&initialCharge.integratedCharge()==savedQ&&
        MaxwellCoupledContract::sourceSnapshot(initialPub)==savedGeometry&&sameClock(initialPub.timeLayer(),first.timeLayer()),
        "F23 fresh/reused preparation must preserve actual initial source VOF/Q/fullclock");
    std::cout<<std::setprecision(17)<<"F23_MAXWELL_GEOMETRY_REUSE_SCENE,"<<(rz?"RZ":"XY")<<",cells,"<<mesh->cellIndexNum()
        <<",actual_phase_stress_clocks,2,timing_samples,3,geometry_build_seconds,"<<geometrySeconds
        <<",fresh_load_seconds,"<<freshSeconds<<",reused_load_seconds,"<<reuseSeconds
        <<",timing_scope,small_public_prepare_only_no_gate_no_whole_App_claim\n";
}
void run(){scene(false);scene(true);std::cout<<"F23_MAXWELL_GEOMETRY_REUSE_PASS,scenes,"<<scenes<<",loads,"<<loadsCompared
    <<",native_values,"<<nativeValuesCompared<<",actual_RHS_power_checks,"<<workValuesCompared<<",guards,"<<guards
    <<",comparison,exact_value_and_signed_zero,real_transport_fields,source_preserved,no_App_commit\n";}
} // namespace F23MaxwellGeometryReuseContract
// END F23 ACTUAL MAXWELL GEOMETRY REUSE CONTRACT

namespace VofEhdOutputContract {
namespace E=LAMP::Electric;
namespace fs=std::filesystem;
using Writer=VofElectrohydrodynamicsOutput;
unsigned scenes=0,frames=0,ownerGuards=0,ioGuards=0;
std::string read(const fs::path& path) {
    std::ifstream in(path,std::ios::binary);check(bool(in),"output contract file missing");
    std::string s((std::istreambuf_iterator<char>(in)),{});check(!in.bad(),"output contract file read failed");return s;
}
bool has(const std::string& value,const std::string& token){return value.find(token)!=std::string::npos;}
struct Block {std::string text;std::size_t nodes=0,elements=0,variables=0;std::vector<double> values;};
Block block(const fs::path& file) {
    Block b;b.text=read(file);
    const auto number=[&](const char* key){const auto p=b.text.find(key);check(p!=std::string::npos,"native DAT dimension missing");return std::stoull(b.text.substr(p+std::string(key).size()));};
    b.nodes=number("NODES=");b.elements=number("ELEMENTS=");
    const auto zone=b.text.find("ZONE T=");check(zone!=std::string::npos,"native DAT zone missing");
    for(std::size_t i=0;i<zone;++i)if(b.text[i]=='"')++b.variables;b.variables/=2;
    const auto aux=b.text.rfind("AUXDATA ");check(aux!=std::string::npos,"native DAT owner metadata missing");
    const auto start=b.text.find('\n',aux);check(start!=std::string::npos,"native DAT body missing");
    std::istringstream numbers(b.text.substr(start+1));double x;while(numbers>>x)b.values.push_back(x);
    check(numbers.eof(),"native DAT numeric payload not fully readable");
    check(has(b.text,"AUXDATA App=\"VofElectrohydrodynamics\""),"wrong actual App name");
    return b;
}
std::vector<std::vector<long double>> rows(const fs::path& file) {
    const auto s=read(file);check(has(s,"VARIABLES=")&&has(s,"Available"),"series lacks explicit numeric availability");
    std::istringstream in(s);std::string line;std::vector<std::vector<long double>> result;
    while(std::getline(in,line)){std::istringstream values(line);long double x;if(!(values>>x))continue;
        std::vector<long double> r{x};while(values>>x)r.push_back(x);check(values.eof(),"series incomplete numeric row");result.push_back(std::move(r));}
    return result;
}
void verifyLocalPlicExecution(const fs::path& output,const HydroTwoPhaseState& accepted,
    const HydroVofStep& step,std::size_t count) {
    const auto file=output/"series/local_plic_execution.dat";const auto actual=rows(file);
    check(actual.size()==count&&actual.back().size()==56,"final-AC LocalPlic execution numeric/Available schema");
    const auto& v=step.transport();const auto s=v.sourceTime(),t=accepted.timeLayer();const auto& q=v.report().localFaceExecution;
    const std::array<long double,28> expected{{
        static_cast<long double>(t.step),
        static_cast<long double>(t.time()),
        static_cast<long double>(t.begin),
        static_cast<long double>(t.end),
        static_cast<long double>(static_cast<int>(t.position)),
        static_cast<long double>(step.candidateState().phaseGeometry().vofIdentity().value()),
        static_cast<long double>(accepted.phaseGeometry().vofIdentity().value()),
        static_cast<long double>(step.report().stages.size()),
        static_cast<long double>(v.sourceIdentity().value()),
        static_cast<long double>(s.step),
        static_cast<long double>(s.begin),
        static_cast<long double>(s.end),
        static_cast<long double>(static_cast<int>(s.position)),
        static_cast<long double>(v.midpoint().identity().value()),
        static_cast<long double>(v.endpoint().identity().value()),
        static_cast<long double>(q.requested),
        static_cast<long double>(q.eligible),
        static_cast<long double>(q.faceValuesParallelUsed),
        static_cast<long double>(q.cellGatherParallelUsed),
        static_cast<long double>(q.fallback),
        static_cast<long double>(q.nativeFaces),
        static_cast<long double>(q.incidences),
        static_cast<long double>(q.directionsUsed),
        static_cast<long double>(q.faceEvaluations),
        static_cast<long double>(q.cellGathers),
        static_cast<long double>(q.scratchBytes),
        static_cast<long double>(q.scratchLimitBytes),
        static_cast<long double>(q.maximumWorkers)
    }};
    const auto& r=actual.back();
    for(std::size_t i=0;i<expected.size();++i){
        const bool available=(i>=20&&i<=25)?q.eligible:i==27?q.directionsUsed!=0:true;
        check(r[expected.size()+i]==(available?1:0),"LocalPlic execution availability is not the actual receipt scope");
        check(r[i]==(available?expected[i]:0),"LocalPlic execution field differs from retained final transport report");
    }
    check(q.scratchLimitBytes==64u*1024u*1024u,"LocalPlic execution must retain the 64MiB optional limit");
    check(has(read(file),"source_geometry_identity")&&has(read(file),"accepted_geometry_identity"),"LocalPlic execution owner columns missing");
}
std::vector<E::ElectricVofBoundary> boundary(const QuadtreeMesh& mesh,const GeometricVofPublication& p,bool rz) {
    std::vector<E::ElectricVofBoundary> b(p.domain().nativeFaceCount());
    for(int n=0;n<static_cast<int>(b.size());++n){const auto f=p.domain().face(mesh,n);
        if(f.negativeCell>=0&&f.positiveCell>=0)b[n]={E::ElectricVofBoundary::Type::Interior,0};
        else if(f.area==0)b[n]={E::ElectricVofBoundary::Type::HomogeneousNeumann,0};
        else b[n]={E::ElectricVofBoundary::Type::Dirichlet,rz?0:.01*f.center.x};}
    return b;
}
void verifyFields(const fs::path& output,int frame,const QuadtreeMesh& mesh,const HydroTwoPhaseState& state,
    const E::ElectricVofCharge& charge,const E::ElectricVofField& field) {
    std::ostringstream suffix;suffix<<std::setw(6)<<std::setfill('0')<<frame;
    const auto volume=block(output/"fields"/("volume_"+suffix.str()+".dat"));
    const auto face=block(output/"fields"/("electric_faces_"+suffix.str()+".dat"));
    for(const auto& entry:fs::directory_iterator(output/"fields"))
        check(entry.is_regular_file(),"accepted VOF fields must be flat, with no category directories");
    const auto& p=state.phaseGeometry().vofPublication();const auto& q=charge.integratedCharge();
    check(volume.nodes==static_cast<std::size_t>(mesh.vertIndexNum())&&volume.elements==p.cells().size()
        &&volume.variables==25&&volume.values.size()==2*volume.nodes+27*volume.elements,"volume 19 borrowed plus 4 electric cell fields/order");
    check(face.elements==field.faces().size()&&face.nodes==2*face.elements&&face.variables==9
        &&face.values.size()==2*face.nodes+9*face.elements,"actual native electric face dimensions/location/connectivity");
    check(has(volume.text,field.isElectrostaticSnapshot()?"ElectricFieldKind=\"electrostaticSnapshot\"":"ElectricFieldKind=\"ohmicEnd\""),"snapshot and real End Ohmic diagnostic conflated");
    if(p.domain().axisymmetric())check((has(volume.text,"\"R\"\n\"Z\"")||has(volume.text,"\"R\"\r\n\"Z\""))&&(has(face.text,"\"R\"\n\"Z\"")||has(face.text,"\"R\"\r\n\"Z\"")),"RZ real coordinate names missing");
    for(std::size_t row=0;row<volume.elements;++row){
        const int native=static_cast<int>(volume.values[2*volume.nodes+row]);check(native==state.values().layout().nativeCell(static_cast<int>(row)),"volume native owner row");
        const std::array<double,4> expected{{q[native],static_cast<double>(static_cast<long double>(q[native])/p.cells()[native].capacity),field.potential()[native],field.epsilonByCell()[native]}};
        for(std::size_t c=0;c<4;++c)check(volume.values[2*volume.nodes+(19+c)*volume.elements+row]==expected[c],"actual Q/rhoQ/phi/epsilon ASCII round trip");
        check(volume.values[2*volume.nodes+volume.elements+row]==static_cast<double>(p.cells()[native].capacity),"true XY/RZ native capacity");
    }
    for(std::size_t n=0;n<field.faces().size();++n){const auto& f=field.faces()[n];const std::array<double,7> expected{{
            double(f.nativeRow),double(f.negativeCell),double(f.positiveCell),f.physicalMeasure,f.integratedField,f.displacement,f.current}};
        for(std::size_t c=0;c<7;++c)check(face.values[2*face.nodes+c*face.elements+n]==expected[c],"actual canonical native H/D/J incidence and measure round trip");}
    check(has(volume.text,"ChargeIdentity=\""+std::to_string(charge.identity())+"\"")&&has(face.text,"ElectricGeometryIdentity=\""+std::to_string(p.identity().value())+"\""),"actual electric output identities missing");
    ++frames;
}
void scene(const fs::path& work,bool rz) {
    ++scenes;auto mesh=meshFor(false),fine=meshFor(true);
    GeometricVofInitializationConfig init;init.subcellsPerAxis=8;
    const auto pub=rz?initial(*mesh,true):GeometricVofInitialization::prepare(*mesh,false,
        [](Vector2D p){const double x=p.x-.47,y=p.y-.48;return x*x+y*y-.23*.23;},init);
    const auto phase=phaseFor(*mesh,pub);const auto& layout=phase.layout();
    HydroTwoPhaseState accepted=HydroTwoPhaseState::createVof(*mesh,phase,HydroNativeValues(layout,
        std::vector<double>(layout.velocityCount()),std::vector<double>(layout.pressureCount())));
    std::vector<double> qi(pub.cells().size());
    for(std::size_t ci=0;ci<qi.size();++ci)qi[ci]=static_cast<double>((ci%2?1e-8L:-1e-8L)*pub.cells()[ci].capacity);
    auto charge=E::ElectricVofCharge::fromIntegrated(*mesh,pub,qi);
    const auto oldU=accepted.values().normalVelocity(),oldP=accepted.values().pressure();const auto oldStocks=stocks(pub);
    VofElectrohydrodynamicsInput input;input.fluid.radial=rz;input.fluid.run.endTime=2e-5;input.fluid.run.writeFields=true;
    input.fluid.run.writeInitial=input.fluid.run.writeFirst=input.fluid.run.writeFinal=true;input.fluid.run.writeEverySteps=1;
    input.fluid.run.diagnosticSeries=true; // This fixture verifies the actual LocalPlic execution receipt.
    input.fluid.vof.liquidDensity=10;input.fluid.vof.gasDensity=1;input.fluid.vof.surfaceTension=0;
    input.fluid.vof.momentumDensityRepresentation=Density::CellAverageProjection;input.fluid.vof.scheme=HydroTwoPhaseTimeScheme::AdaptiveLinearMidpoint;
    input.materials=rz?E::ElectricVofMaterials{2,2,0,0}:E::ElectricVofMaterials{2,5,0,0};
    for(const char* name:{"app.dict","mesh.dict","initial.dict","physics.dict","time.dict","output.dict","solver.dict","electric.dict"})
        input.fluid.documents.push_back({name,"// actual output fixture document\n","// actual output fixture effective document\n"});
    const auto root=work/(rz?"rz":"xy");Writer writer(input,{root});writer.begin();
    const auto initialField=E::ElectricVofField::prepareElectrostaticSnapshot(*mesh,pub,charge,input.materials,boundary(*mesh,pub,rz));
    writer.accepted(*mesh,accepted,charge,initialField,nullptr,Writer::Clock::now(),0,false);verifyFields(root/"output",0,*mesh,accepted,charge,initialField);
    const auto initialIndex=read(root/"output/meta/frames.csv");
    const auto foreign=E::ElectricVofCharge::fromIntegrated(*mesh,pub,qi);
    rejects([&]{writer.accepted(*mesh,accepted,foreign,initialField,nullptr,Writer::Clock::now(),0,false);},"reject same-clock foreign accepted Q owner");
    ++ownerGuards;check(read(root/"output/meta/frames.csv")==initialIndex,"owner rejection modified frame index");
    auto config=input.fluid.vof;
    const auto prepare=[&](const HydroTwoPhaseState& source,const E::ElectricVofCharge& sourceQ,const Clock& target){
        auto cfg=config;cfg.maxwellLoad=[&](const QuadtreeMesh& grid,const GeometricVofTransport& volume,const HydroPhaseGeometry& mid,const HydroPressureLoad& pressure){
            const auto moved=E::ElectricVofChargeTransport::prepareClosed(grid,sourceQ,volume);
            const auto field=E::ElectricVofField::prepareOhmic(grid,mid.vofPublication(),moved.midpoint(),input.materials,boundary(grid,mid.vofPublication(),rz));
            return HydroMaxwellLoad::prepare(grid,pressure,mid,E::ElectricVofMaxwellStress::prepare(grid,field));};
        return HydroVofStep::prepare(*mesh,source,target,cfg);
    };
    const Clock end1{1,0,1e-5,Clock::Position::End};const auto step1=prepare(accepted,charge,end1);
    const auto moved1=E::ElectricVofChargeTransport::prepareClosed(*mesh,charge,step1.transport());
    const auto field1=E::ElectricVofField::prepareOhmic(*mesh,step1.candidateState().phaseGeometry().vofPublication(),moved1.endpoint(),input.materials,
        boundary(*mesh,step1.candidateState().phaseGeometry().vofPublication(),rz));
    const auto midField=E::ElectricVofField::prepareOhmic(*mesh,step1.transport().midpoint(),moved1.midpoint(),input.materials,boundary(*mesh,step1.transport().midpoint(),rz));
    rejects([&]{writer.accepted(*mesh,step1.candidateState(),field1.chargeCandidate(),midField,&step1,Writer::Clock::now(),0,false);},"reject actual Mid field presented as accepted End");
    ++ownerGuards;check(read(root/"output/meta/frames.csv")==initialIndex,"clock rejection modified frame index");
    writer.accepted(*mesh,step1.candidateState(),field1.chargeCandidate(),field1,&step1,Writer::Clock::now(),0,false);
    verifyFields(root/"output",1,*mesh,step1.candidateState(),field1.chargeCandidate(),field1);
    verifyLocalPlicExecution(root/"output",step1.candidateState(),step1,1);
    accepted=step1.candidateState();charge=field1.chargeCandidate();
    check(!field1.isElectrostaticSnapshot()&&field1.sourceCharge().sameInstance(moved1.sourceCharge())&&!field1.sourceCharge().sameInstance(charge),"ordinary End must retain actual old source charge");
    const Clock end2{2,end1.end,2e-5,Clock::Position::End};const auto step2=prepare(accepted,charge,end2);
    const auto moved2=E::ElectricVofChargeTransport::prepareClosed(*mesh,charge,step2.transport());
    const auto& physical=step2.candidateState();const auto& endPub=physical.phaseGeometry().vofPublication();
    const auto field2=E::ElectricVofField::prepareOhmic(*mesh,endPub,moved2.endpoint(),input.materials,boundary(*mesh,endPub,rz));
    auto transfer=GeometricVofMeshTransfer::prepare(*mesh,endPub,*fine,domain(*fine,rz),end2);
    HydroVofMeshTransferConfig remapConfig;remapConfig.momentumDensityRepresentation=Density::CellAverageProjection;
    const auto fluidRemap=HydroVofMeshTransfer::prepare(*mesh,physical,*fine,std::move(transfer),10,1,remapConfig);
    const auto& actualTransfer=fluidRemap.interfaceTransfer();const auto& target=fluidRemap.candidateState();
    const auto qremap=E::ElectricVofChargeMeshTransfer::prepare(*mesh,endPub,field2.chargeCandidate(),*fine,actualTransfer);
    const auto& targetQ=qremap.chargeCandidate();
    const auto snapshot=E::ElectricVofField::prepareElectrostaticSnapshot(*fine,target.phaseGeometry().vofPublication(),targetQ,input.materials,boundary(*fine,target.phaseGeometry().vofPublication(),rz));
    HydrodynamicsVofRegridReport decision;decision.sourceClock=end2;decision.sourceLeaves=mesh->cellIndexNum();decision.targetLeaves=fine->cellIndexNum();
    decision.requested=decision.changed=true;
    const HydrodynamicsVofRemapReceipt receipt{mesh.get(),&physical,&step2,&actualTransfer,&fluidRemap,&decision};
    const auto prefix=read(root/"output/meta/frames.csv");
    rejects([&]{writer.accepted(*fine,target,targetQ,snapshot,&step2,Writer::Clock::now(),0,true,&receipt,nullptr);},"reject fluid remap without actual signed Q transfer receipt");
    ++ownerGuards;check(read(root/"output/meta/frames.csv")==prefix,"missing Q transfer modified complete index");
    const auto targetStocks=stocks(target.phaseGeometry().vofPublication());const auto targetValues=target.values().normalVelocity();const auto targetCharge=targetQ.integratedCharge();
    writer.accepted(*fine,target,targetQ,snapshot,&step2,Writer::Clock::now(),0,true,&receipt,&qremap);writer.complete();
    verifyFields(root/"output",2,*fine,target,targetQ,snapshot);
    verifyLocalPlicExecution(root/"output",target,step2,2);
    const auto manifest=read(root/"output/meta/manifest.json"),index=read(root/"output/meta/frames.csv");
    check(has(manifest,"\"app\":\"VofElectrohydrodynamics\"")&&has(manifest,"\"status\":\"COMPLETED\"")&&has(manifest,"\"complete_frames\":3"),"actual App manifest/completed three frames");
    check(std::count(index.begin(),index.end(),'\n')==4&&has(index,"ohmicEnd")&&has(index,"electrostaticSnapshot"),"payload-first index complete frame sources");
    const auto maxwell=rows(root/"output/series/maxwell_work.dat");check(maxwell.size()==3&&maxwell.back().size()==16,"last-stage Mid work explicit schema/order");
    const auto& actualWork=step2.report().stages.back().maxwellWork;
    std::ostringstream rendered;rendered<<std::scientific<<std::setprecision(16)<<actualWork;
    check(maxwell.back()[4]==std::stold(rendered.str()),"actual final AC Maxwell work serialization not accumulated across alternative stages");
    check(static_cast<double>(maxwell.back()[2])==step2.frozenMidpoint().forces().maxwellLoad->phaseGeometry().timeLayer().time()&&maxwell.back()[7]==1,"actual Mid force clock and remapped bit");
    check(rows(root/"output/series/regrid.dat").size()==1&&rows(root/"output/series/regrid_decision.dat").size()==1,"only real changed remap has independent phase and signed Q ledger");
    check(read(root/"output/series/energy.dat").find("surfaceAvailable")!=std::string::npos,"undefined surface energy availability missing");
    check(has(manifest,"\"series_profile\":\"diagnostic\""),"diagnostic series profile missing from manifest");
    // Replay the same immutable accepted physical owners with the default
    // series selection, including the genuine target-mesh remap receipt.
    auto compactInput=input;compactInput.fluid.run.diagnosticSeries=RunControls{}.diagnosticSeries;
    check(!compactInput.fluid.run.diagnosticSeries,"default run series profile must be physical");
    const auto compactRoot=work/(rz?"rz_physical":"xy_physical");Writer compact(compactInput,{compactRoot});compact.begin();
    compact.accepted(*mesh,step1.beginState(),moved1.sourceCharge(),initialField,nullptr,Writer::Clock::now(),0,false);
    compact.accepted(*mesh,step1.candidateState(),field1.chargeCandidate(),field1,&step1,Writer::Clock::now(),0,false);
    compact.accepted(*fine,target,targetQ,snapshot,&step2,Writer::Clock::now(),0,true,&receipt,&qremap);compact.complete();
    check(has(read(compactRoot/"output/meta/manifest.json"),"\"series_profile\":\"physical\""),"physical series profile missing from manifest");
    check(read(compactRoot/"output/meta/frames.csv")==index,"series profile changed the accepted electric frame index");
    for(const auto& entry:fs::directory_iterator(root/"output/fields"))
        check(read(compactRoot/"output/fields"/entry.path().filename())==read(entry.path()),"series profile changed an accepted physical field byte");
    const std::array<std::string,14> allowed{{"conservation","energy","solver","vof_csf","electric","maxwell_work","charge_balance",
        "phase_boundary","vof_boundary_work","regrid_projection_work","time_step","regrid","regrid_decision","events"}};
    for(const auto& entry:fs::directory_iterator(compactRoot/"output/series")){
        const auto name=entry.path().stem().string();
        check(entry.path().extension()==".dat"&&std::find(allowed.begin(),allowed.end(),name)!=allowed.end(),"physical profile emitted developer series");
        check(read(entry.path())==read(root/"output/series"/entry.path().filename()),"series profile changed a retained physical history byte");
    }
    for(const auto* name:{"conservation","energy","solver","electric","maxwell_work","regrid","regrid_decision"})
        check(fs::is_regular_file(compactRoot/"output/series"/(std::string(name)+".dat")),"physical profile omitted a process/solver/remap history");
    check(stocks(pub)==oldStocks&&qi==moved1.sourceCharge().integratedCharge()&&step1.beginState().values().normalVelocity()==oldU&&step1.beginState().values().pressure()==oldP
        &&stocks(target.phaseGeometry().vofPublication())==targetStocks&&target.values().normalVelocity()==targetValues&&targetQ.integratedCharge()==targetCharge,"passive writer modified a source or target physical owner");
    // Real later-payload write failure must not publish its already-written volume.
    const auto failedRoot=work/(rz?"rz_io_failure":"xy_io_failure");Writer failing(input,{failedRoot});failing.begin();
    failing.accepted(*mesh,step1.beginState(),moved1.sourceCharge(),initialField,nullptr,Writer::Clock::now(),0,false);
    const auto failurePrefix=read(failedRoot/"output/meta/frames.csv");
    fs::create_directory(failedRoot/"output/fields/electric_faces_000001.dat");
    rejects([&]{failing.accepted(*mesh,step1.candidateState(),field1.chargeCandidate(),field1,&step1,Writer::Clock::now(),0,false);},"actual later frame payload I/O failure");
    ++ioGuards;check(read(failedRoot/"output/meta/frames.csv")==failurePrefix,"failed electric-face payload published a partial frame");
    rejects([&]{failing.complete();},"poisoned output cannot complete");++ioGuards;
    failing.failure(2,1,end1.end,"output","D2 expected electric-face path is directory");
    check(has(read(failedRoot/"output/meta/manifest.json"),"\"status\":\"FAILED\"")&&has(read(failedRoot/"output/meta/failure.json"),"\"accepted_state_preserved\":true"),"native failure metadata missing after payload error");
    // Model the actual commit -> adaptive-history I/O -> accepted-output order:
    // the real first-step state exists, but it has no successful writer witness.
    const auto clockRoot=work/(rz?"rz_unobserved_clock":"xy_unobserved_clock");Writer unobserved(input,{clockRoot});unobserved.begin();
    unobserved.accepted(*mesh,step1.beginState(),moved1.sourceCharge(),initialField,nullptr,Writer::Clock::now(),0,false);
    const auto unobservedPrefix=read(clockRoot/"output/meta/frames.csv");
    unobserved.failure(2,1,end1.end,"adaptive history","D2 expected committed clock precedes accepted field output");
    const auto unavailable=read(clockRoot/"output/meta/manifest.json");
    check(has(unavailable,"\"geometry_identity\":0")&&has(unavailable,"\"charge_identity\":0")&&has(unavailable,"\"geometry_identity_available\":false")
        &&has(unavailable,"\"charge_identity_available\":false")&&read(clockRoot/"output/meta/frames.csv")==unobservedPrefix,"old geometry/Q IDs mislabeled at unobserved new accepted clock");
    ++ioGuards;
}
void run(const fs::path& work) {
    check(!fs::exists(work),"output contract requires unique fresh D work path");fs::create_directories(work);
    scene(work,false);scene(work,true);
    std::cout<<"VOF_EHD_OUTPUT_PASS,scenes,"<<scenes<<",verified_frames,"<<frames<<",owner_guards,"<<ownerGuards<<",io_guards,"<<ioGuards
        <<",actual_End_Ohmic_and_target_snapshot,actual_native_H_D_J,passive_Midwork,independent_remap_ledgers,no_registered_App_claim\n";
}
} // namespace VofEhdOutputContract

namespace VofOpenRemapContract {
namespace E=LAMP::Electric;
namespace FV=LAMP::Geometry::Detail::FiniteVolume;
using Mode=GeometricVofAdvectorSlab::Boundary::Mode;
unsigned sceneCount=0,remapCount=0,stepCount=0,guardCount=0,electricCount=0;
unsigned pressureShiftRemaps=0,gaugeRejections=0;
long double peakEnergy=0,peakDiv=0,maxBoundaryWork=0,maxCorrection=0;
void bound(long double residual,long double scale,long double tolerance,const char* why) {
    check(std::isfinite(residual)&&std::isfinite(scale)&&scale>=0,why);
    check(scale>0?std::abs(residual)<=tolerance*scale:residual==0,why);
}
template<class F>void reject(F&& f,const char* why){rejects(std::forward<F>(f),why);++guardCount;}
std::unique_ptr<QuadtreeMesh> meshForOpen(bool fine,bool rz) {
    DictionaryMesh d;d.Nx=d.Ny=4;d.minCoord={rz?1.:0.,0};d.maxCoord={rz?2.:1.,1};
    d.maxInitialIterations=0;d.maxCellNum=4096;d.maxLevel_local={1};auto mesh=std::make_unique<QuadtreeMesh>(d);
    const double left=rz?1.:0.;DictionaryMesh::RefineBox upper;upper.level=1;upper.min={left,.5};upper.max={left+1,1};
    std::vector<DictionaryMesh::RefineBox> boxes{upper};if(fine){auto lower=upper;lower.min={left,0};lower.max={left+.5,.5};boxes.push_back(lower);}
    mesh->applyRefineBoxes(boxes);return mesh;
}
HydroVofPrescribedFlow flowFor(const QuadtreeMesh& mesh,const HydroDofLayout& l,bool outlet,
    const std::vector<double>& normalByNative) {
    HydroVofPrescribedFlow f;f.rulesByNativeRow.assign(l.domain().nativeFaceCount(),Mode::Interior);f.normalVelocity.assign(l.velocityCount(),0);
    for(int u=0;u<l.velocityCount();++u){const auto face=l.domain().face(mesh,l.nativeFace(u));if(face.outwardSign==0)continue;
        if(face.horizontal||face.area==0){f.rulesByNativeRow[face.nativeRow]=Mode::Closed;continue;}
        f.rulesByNativeRow[face.nativeRow]=face.outwardSign<0?Mode::PureLiquidReservoir:Mode::PureGasReservoir;
        if(outlet&&face.outwardSign>0)f.pressureOutlets.push_back({u,.125});
        else{f.prescribedVelocityDofs.push_back(u);f.normalVelocity[u]=normalByNative.at(face.nativeRow);}}
    const auto native=l.domain();const auto* grid=&mesh;const auto normals=normalByNative;
    // Retains only THIS target tree/domain, never a borrowed old native-row map.
    f.endpointVelocity=[native,grid,normals](int row,Clock endpoint){const auto face=native.face(*grid,row);
        return face.outwardSign!=0&&!face.horizontal?Vector2D(normals.at(row),.03125*(1+endpoint.end)):Vector2D(0,0);};
    return f;
}
std::vector<E::ElectricVofBoundary> electricBoundary(const QuadtreeMesh& mesh,const GeometricVofPublication& pub) {
    std::vector<E::ElectricVofBoundary> b(pub.domain().nativeFaceCount());
    for(int n=0;n<static_cast<int>(b.size());++n){const auto f=pub.domain().face(mesh,n);
        b[n]=f.negativeCell>=0&&f.positiveCell>=0?E::ElectricVofBoundary{E::ElectricVofBoundary::Type::Interior,0}:
            (f.area==0?E::ElectricVofBoundary{E::ElectricVofBoundary::Type::HomogeneousNeumann,0}:
            E::ElectricVofBoundary{E::ElectricVofBoundary::Type::Dirichlet,.01*f.center.x+.02*f.center.y});}return b;
}
void physicalStep(const HydroVofStep& step,const HydroTwoPhaseState& old,const QuadtreeMesh& mesh) {
    ++stepCount;const auto& volume=step.transport();const auto& next=step.candidateState();const auto& pub=old.phaseGeometry().vofPublication();
    check(volume.advector().boundary()&&next.vofPrescribedFlow()&&sameClock(next.timeLayer(),volume.endTime()),"F5 real F3 step/transport/flow owner");
    for(const auto& r:step.report().stages)check(r.saddle.algebraicRelativeResidual<=1e-11&&r.saddle.freeMomentumDefect<=1e-11
        &&r.saddle.maxCellContinuityDefect<=1e-11&&r.momentumPowerRelativeDefect<=1e-11&&r.energy.closureRelativeResidual<=1e-11L,"F5 next original F3 gates");
    for(bool half:{true,false}){const auto& endpoint=half?volume.midpoint():volume.endpoint();const auto& ledger=half?volume.report().boundaryHalf:volume.report().boundaryFull;
        for(int p=0;p<2;++p){FV::Dyadic change,in,out;long double capacity=0;
            for(std::size_t ci=0;ci<pub.plicCells().size();++ci){change.add(endpoint.plicCells()[ci].inventory[p]);change.add(pub.plicCells()[ci].inventory[p],-1);capacity+=pub.cells()[ci].capacity;}
            for(const auto& face:volume.faces()){if(face.negativeCell>=0&&face.positiveCell>=0)continue;const auto& q=half?face.half:face.full;
                const auto positive=p==0?q.positiveLiquid:q.positiveGas,negative=p==0?q.negativeLiquid:q.negativeGas;
                out.add(face.negativeCell>=0?positive:negative);in.add(face.negativeCell>=0?negative:positive);}
            FV::Dyadic residual=change;residual.add(out);residual.add(in,-1);
            bound(residual.rounded(),capacity,1e-12L,"F5 independent actual open phase Delta+Out-In");
            bound(ledger.inflow[p]-in.rounded(),capacity,1e-12L,"F5 actual open phase inflow ledger");
            bound(ledger.outflow[p]-out.rounded(),capacity,1e-12L,"F5 actual open phase outflow ledger");}
        check(ledger.enabled&&ledger.physicalRelativeConservation<=1e-12L&&ledger.massRelativeConservation<=1e-12L
            &&ledger.algebraicRelativeConservation<=1e-12L,"F5 original boundary WY/physical/mass gates");}
    for(const auto& p:next.vofPrescribedFlow()->pressureOutlets){const auto f=next.values().layout().domain().face(mesh,next.values().layout().nativeFace(p.velocityDof));
        check(f.outwardSign*next.values().normalVelocity()[p.velocityDof]>0,"F5 physical F3 outlet reversed");}
}
void verify(const QuadtreeMesh& sourceMesh,const HydroTwoPhaseState& source,const QuadtreeMesh& targetMesh,
    const HydroVofMeshTransfer& remap,bool outlet) {
    ++remapCount;const auto& target=remap.candidateState();const auto& r=remap.report();const auto& l=target.values().layout();
    const auto& flow=*target.vofPrescribedFlow();const auto& mass=*target.retainedMass();const auto& u=target.values().normalVelocity();const auto& v=remap.rawValues().normalVelocity();
    const auto& lambda=remap.projectionMultiplier();const auto& b=remap.projectionMomentumRhs();const auto& metrics=target.phaseGeometry().metrics();
    const auto pressure=HydroPressureLoad::build(metrics.divergence());const auto& C=pressure.matrix();const auto& D=pressure.divergence().matrix();
    check(remap.sourceState().sameIdentity(source)&&target.vofPrescribedFlow()&&r.prescribedFlow&&r.historyReset&&!target.velocityHistory()&&r.remapSolves==1,"F5 actual open owner/remap solve/history reset");
    check(sameClock(target.timeLayer(),source.timeLayer())&&target.phaseGeometry().vofIdentity()==remap.interfaceTransfer().candidate().identity()
        &&!source.values().layout().domain().uniformLeaves()&&!l.domain().uniformLeaves(),"F5 same-End true native 2:1 source/target");
    check(mass.representation()==Density::CellAverageProjection&&mass.geometry().sameIdentity(target.momentumGeometry()),"F5 retained actual target M/geometry");
    check(r.projectionSaddle.pressureGaugeCount==(outlet?0:1)&&r.pressureGaugeApplied==!outlet,"F5 physical/projection pressure anchoring");
    check(r.projectionSaddle.algebraicRelativeResidual<=1e-11&&r.projectionSaddle.freeMomentumDefect<=1e-11&&r.projectionSaddle.maxCellContinuityDefect<=1e-11
        &&r.maximumFinalRelativeDivergence<=1e-11&&r.projectionEnergyIdentityRelativeResidual<=2e-11L,"F5 original projection gates");
    std::vector<bool> fixed(u.size());for(int k:target.boundary().velocityDofs()){fixed[k]=true;check(u[k]==0,"F5 closed/axis wall changed");}
    for(int k:flow.prescribedVelocityDofs){fixed[k]=true;check(u[k]==flow.normalVelocity[k]&&v[k]==flow.normalVelocity[k],"F5 prescribed normal was clamped from different raw seed");}
    check(r.fixedNormals==flow.prescribedVelocityDofs.size()&&r.pressureOutletNormals==flow.pressureOutlets.size()&&r.boundaryCoverageSpans>0,"F5 exact physical boundary coverage roles");
    for(const auto& port:flow.pressureOutlets){const auto f=l.domain().face(targetMesh,l.nativeFace(port.velocityDof));
        check(!fixed[port.velocityDof]&&flow.normalVelocity[port.velocityDof]==0&&f.outwardSign*u[port.velocityDof]>0,"F5 free outlet not outward/free");}
    FV::Dyadic delta,correction,WB,WRHS,freeWork,continuity,scale;
    for(int k=0;k<l.velocityCount();++k){const long double M=mass.masses()[k],a=u[k],raw=v[k];
        // Independent active kinetic scale, never a unit/mass floor or a
        // residual-sized denominator that mistakes roundoff for lost energy.
        scale.add(FV::Dyadic::product({.5L,M,a,a}));scale.add(FV::Dyadic::product({.5L,M,raw,raw}));
        check(FV::roundDouble(FV::Dyadic::product({M,raw})).value==b[k],"F5 projection RHS not actual independent RN(Mv)");
        FV::Dyadic dk=FV::Dyadic::product({.5L,M,a,a});dk.add(FV::Dyadic::product({.5L,M,raw,raw}),-1);delta.add(dk);scale.add(dk.absolute());
        FV::Dyadic ec=FV::Dyadic::product({.5L,M,a,a});ec.add(FV::Dyadic::product({.5L,M,raw,raw}));ec.add(FV::Dyadic::product({M,a,raw}),-1);
        check(ec.compare(FV::Dyadic{})>=0,"F5 exact positive M correction norm");correction.add(ec);scale.add(ec.absolute());
        FV::Dyadic wr=FV::Dyadic::product({a,b[k]});wr.add(FV::Dyadic::product({a,M,raw}),-1);WRHS.add(wr);scale.add(wr.absolute());
        FV::Dyadic work=FV::Dyadic::product({a,M,a});work.add(FV::Dyadic::product({a,b[k]}),-1);
        for(HydroPressureLoad::Matrix::InnerIterator it(C,k);it;++it)work.add(FV::Dyadic::product({a,it.value(),lambda[it.col()]}));
        if(fixed[k]){WB.add(work);scale.add(work.absolute());}else freeWork.add(work);
    }
    for(int p=0;p<l.pressureCount();++p){FV::Dyadic div;long double gross=0;
        for(HydroIntegralDivergence::Matrix::InnerIterator it(D,p);it;++it){div.add(FV::Dyadic::product({it.value(),u[it.col()]}));gross+=std::abs(static_cast<long double>(it.value())*u[it.col()]);
            continuity.add(FV::Dyadic::product({lambda[p],it.value(),u[it.col()]}));}
        bound(div.rounded(),gross,1e-11L,"F5 independent actual D u");if(gross>0)peakDiv=std::max(peakDiv,std::abs(div.rounded())/gross);}
    FV::Dyadic left=delta;left.add(correction);left.add(WB,-1);left.add(WRHS,-1);FV::Dyadic right=freeWork;right.add(continuity);
    check(left.compare(right)==0,"F5 exact C=-D^T projection energy algebra identity");
    bound(left.rounded(),scale.rounded(),2e-11L,"F5 independent uncompensated projection energy/no floor");
    const auto compare=[&](long double report,const FV::Dyadic& expected,const char* why){bound(report-expected.rounded(),r.projectionEnergyIdentityScale,1e-12L,why);};
    compare(r.projectionKineticChange,delta,"F5 actual independent projection DeltaK");compare(r.projectionCorrectionEnergy,correction,"F5 actual independent Ecor");
    compare(r.projectionBoundaryWork,WB,"F5 fixed reaction WB");compare(r.projectionRhsRoundingWork,WRHS,"F5 actual RHS rounding WRHS");
    compare(r.projectionFreeResidualWork,freeWork,"F5 free residual work");compare(r.projectionContinuityWork,continuity,"F5 actual lambda D u");
    compare(r.projectionEnergyIdentityResidual,left,"F5 independent uncorrected energy report");
    if(scale.rounded()>0)peakEnergy=std::max(peakEnergy,std::abs(left.rounded())/scale.rounded());
    maxBoundaryWork=std::max(maxBoundaryWork,std::abs(WB.rounded()));maxCorrection=std::max(maxCorrection,correction.rounded());
    // Independent source/raw/target native M/P/KE budgets; no total-P
    // conservation is asserted for this constrained remap law.
    const std::array<const HydroMomentumMass*,3> masses{{source.retainedMass(),&mass,&mass}};
    const std::array<const HydroNativeValues*,3> values{{&source.values(),&remap.rawValues(),&target.values()}};
    const std::array<const HydroVofRemapNativeBudget*,3> reports{{&r.source,&r.rawTarget,&r.target}};
    std::array<std::array<long double,2>,3> momentumActivity{};
    std::array<HydroVofRemapNativeBudget,3> nativeBudgets;
    for(int j=0;j<3;++j){check(masses[j]!=nullptr,"F5 real source/target mass owner missing");nativeBudgets[j]=budget(*masses[j],*values[j]);
        for(int k=0;k<masses[j]->layout().velocityCount();++k){const int c=masses[j]->layout().normalAxis(k)==HydroDofLayout::NormalAxis::First?0:1;
            momentumActivity[j][c]+=masses[j]->masses()[k]*std::abs(static_cast<long double>(values[j]->normalVelocity()[k]));}
        for(int c=0;c<2;++c){bound(reports[j]->mass[c]-nativeBudgets[j].mass[c],nativeBudgets[j].mass[c],1e-12L,"F5 actual native mass budget");
            bound(reports[j]->momentum[c]-nativeBudgets[j].momentum[c],momentumActivity[j][c],1e-12L,"F5 actual native momentum budget");}
        bound(reports[j]->kinetic-nativeBudgets[j].kinetic,nativeBudgets[j].kinetic,1e-12L,"F5 actual native kinetic budget");}
    for(int c=0;c<2;++c){bound(r.momentumChange[c]-(nativeBudgets[2].momentum[c]-nativeBudgets[0].momentum[c]),
            momentumActivity[0][c]+momentumActivity[2][c],1e-12L,"F5 source-to-target native DeltaP report");
        bound(r.projectionMomentumChange[c]-(nativeBudgets[2].momentum[c]-nativeBudgets[1].momentum[c]),
            momentumActivity[1][c]+momentumActivity[2][c],1e-12L,"F5 raw-to-target native DeltaP report");}
    // Physical pressure is independently overlap averaged; lambda is unrelated.
    const auto& overlap=remap.interfaceTransfer().nativeOverlaps();long double pressureMean=0,pressureL1=0;
    for(int p=0;p<l.pressureCount();++p){FV::Dyadic integral,activity;
        for(auto id:overlap.byTarget()[l.nativeCell(p)]){const auto& piece=overlap.pieces()[id];const auto oldP=source.values().pressure()[source.values().layout().pressureDof(piece.sourceCell)];
            const auto term=FV::Dyadic::product({piece.capacity,oldP});integral.add(term);activity.add(term.absolute());}
        const auto volume=static_cast<long double>(metrics.cellVolumes()[p]);const auto expected=static_cast<double>(static_cast<double>(integral.rounded()/volume)-r.pressureGaugeShift);
        bound(target.values().pressure()[p]-expected,activity.rounded()/volume+std::abs(r.pressureGaugeShift),1e-12L,"F5 independently mapped physical p");
        pressureMean+=volume*target.values().pressure()[p];pressureL1+=volume*std::abs(target.values().pressure()[p]);}
    if(outlet)check(r.pressureGaugeShift==0,"F5 imposed an artificial gauge on physical outlet pressure");
    else bound(pressureMean,pressureL1,1e-11L,"F5 all-fixed physical pressure volume gauge");
    const auto a=stocks(source.phaseGeometry().vofPublication()),z=stocks(target.phaseGeometry().vofPublication());
    for(int p=0;p<2;++p)bound(z[p]-a[p],std::abs(z[p])+std::abs(a[p]),1e-12L,"F5 independent retained two-phase stocks");
    // Independent port-side totals verify physical coverage, including true RZ Af.
    std::array<long double,4> oldFlux{},rawFlux{};
    for(int k=0;k<source.values().layout().velocityCount();++k){const auto f=source.values().layout().domain().face(sourceMesh,source.values().layout().nativeFace(k));if(f.outwardSign==0)continue;
        const int side=f.horizontal?(f.outwardSign<0?2:3):(f.outwardSign<0?0:1);oldFlux[side]+=static_cast<long double>(f.area)*source.values().normalVelocity()[k];}
    for(int k=0;k<l.velocityCount();++k){const auto f=l.domain().face(targetMesh,l.nativeFace(k));if(f.outwardSign==0)continue;
        const int side=f.horizontal?(f.outwardSign<0?2:3):(f.outwardSign<0?0:1);rawFlux[side]+=static_cast<long double>(f.area)*v[k];}
    for(int side=0;side<4;++side)bound(rawFlux[side]-oldFlux[side],std::abs(rawFlux[side])+std::abs(oldFlux[side]),1e-12L,"F5 raw seed preserves actual covered physical side flux");
    std::cout<<std::setprecision(17)<<"VOF_OPEN_REMAP,"<<(l.domain().axisymmetric()?"RZ":"XY")<<",pressure_outlet,"<<outlet<<",step,"<<target.timeLayer().step
        <<",source_cells,"<<sourceMesh.cellIndexNum()<<",target_cells,"<<targetMesh.cellIndexNum()<<",projection_DeltaK,"<<delta.rounded()<<",Ecor,"<<correction.rounded()
        <<",WB,"<<WB.rounded()<<",WRHS,"<<WRHS.rounded()<<",uncompensated_energy,"<<left.rounded()<<",free_residual_work,"<<freeWork.rounded()
        <<",lambda_Du,"<<continuity.rounded()<<",history_reset,"<<r.historyReset<<",pressure_gauge,"<<r.pressureGaugeApplied<<"\n";
}
void electric(const QuadtreeMesh& oldMesh,const HydroTwoPhaseState& old,const QuadtreeMesh& newMesh,const HydroVofMeshTransfer& remap) {
    const auto& pub=old.phaseGeometry().vofPublication();std::vector<double> qi(pub.cells().size());
    for(std::size_t ci=0;ci<qi.size();++ci)qi[ci]=static_cast<double>((ci%2?.125L:-.0625L)*pub.cells()[ci].capacity);
    const auto charge=E::ElectricVofCharge::fromIntegrated(oldMesh,pub,qi);const auto savedQ=charge.integratedCharge();const auto oldGeometry=MaxwellCoupledContract::sourceSnapshot(pub);
    const auto transfer=E::ElectricVofChargeMeshTransfer::prepare(oldMesh,pub,charge,newMesh,remap.interfaceTransfer());
    const auto& candidate=transfer.chargeCandidate();const auto& target=remap.candidateState().phaseGeometry().vofPublication();
    transfer.validateFor(oldMesh,pub,charge,newMesh,remap.interfaceTransfer());
    const E::ElectricVofMaterials material{2,1,.01,.02};
    const auto field=E::ElectricVofField::prepareElectrostaticSnapshot(newMesh,target,candidate,material,electricBoundary(newMesh,target));
    field.validateFor(newMesh,target,candidate,old.timeLayer());
    check(field.isElectrostaticSnapshot()&&field.effectiveDt()==0&&field.chargeCandidate().sameInstance(candidate)
        &&sameClock(field.candidateTime(),old.timeLayer())&&field.geometry().identity()==target.identity(),"F5 E2 target D0 actual sameEnd owner/clock");
    FV::Dyadic a,b,activity;for(double q:qi){a.add(q);activity.add(std::abs(q));}for(double q:candidate.integratedCharge()){b.add(q);activity.add(std::abs(q));}
    FV::Dyadic residual=b;residual.add(a,-1);bound(residual.rounded(),activity.rounded(),1e-10L,"F5 independent signed sameEnd Q conservation");
    check(transfer.report().relativeChargeResidual<=1e-10&&field.report().relativeLinearResidual<=1e-11&&field.report().maximumRelativeGaussResidual<=1e-10
        &&charge.integratedCharge()==savedQ&&MaxwellCoupledContract::sourceSnapshot(pub)==oldGeometry,"F5 E2/D0 original gates and immutable source");
    const auto foreign=E::ElectricVofCharge::fromIntegrated(oldMesh,pub,qi);
    reject([&]{transfer.validateFor(oldMesh,pub,foreign,newMesh,remap.interfaceTransfer());},"F5 foreign samevalue Q remap owner");
    ++electricCount;std::cout<<"VOF_OPEN_REMAP_E2_D0,same_End_signed_Q_source,explicitly_not_F4_time_advance,source_step,"<<old.timeLayer().step
        <<",target_Q_identity,"<<candidate.identity()<<",relative_Q,"<<transfer.report().relativeChargeResidual<<",Gauss,"<<field.report().maximumRelativeGaussResidual<<"\n";
}
void scene(bool rz,bool outlet) {
    ++sceneCount;auto coarse=meshForOpen(false,rz),fine=meshForOpen(true,rz);const auto pub=GeometricVofPublication::prepareLocalPlic(*coarse,domain(*coarse,rz),{},
        std::vector<double>(coarse->cellIndexNum(),0));const auto phase=phaseFor(*coarse,pub);const auto& l=phase.layout();std::vector<double> u(l.velocityCount()),normal(l.domain().nativeFaceCount());
    for(int k=0;k<l.velocityCount();++k){const auto f=l.domain().face(*coarse,l.nativeFace(k));if(!f.horizontal)u[k]=rz?.125/f.center.x:.125;
        const auto& edge=f.horizontal?coarse->hEdges(f.edgeLabel):coarse->vEdges(f.edgeLabel);auto a=coarse->verts(edge.edgeToVert[0]).coord(),b=coarse->verts(edge.edgeToVert[1]).coord();
        if(f.horizontal?a.x>b.x:a.y>b.y)std::swap(a,b);if(rz){a.x-=1;b.x-=1;}
        const auto flux=.03125L*(f.horizontal?-1.L:1.L)*(rz?CartesianMeasure::twoPi:1.L)*(psi(b)-psi(a));
        if(f.area>0)u[k]+=static_cast<double>(flux/f.area);normal[f.nativeRow]=u[k];}
    auto flow=flowFor(*coarse,l,outlet,normal);auto accepted=HydroTwoPhaseState::createVofPrescribedFlow(*coarse,phase,
        HydroNativeValues(l,u,std::vector<double>(l.pressureCount(),outlet?.125:0)),flow);
    HydroVofAcConfig config;config.liquidDensity=2;config.gasDensity=1;config.momentumDensityRepresentation=Density::CellAverageProjection;
    config.surfaceTension=config.liquidViscosity=config.gasViscosity=0;constexpr double dt=1./16384;
    const auto first=HydroVofStep::preparePrescribedFlow(*coarse,accepted,{1,0,dt,Clock::Position::End},config);physicalStep(first,accepted,*coarse);first.commit(*coarse,accepted);
    HydroVofMeshTransferConfig transferConfig;transferConfig.momentumDensityRepresentation=Density::CellAverageProjection;
    for(int cycle=0;cycle<2;++cycle){const QuadtreeMesh& sourceMesh=cycle?*fine:*coarse;const QuadtreeMesh& targetMesh=cycle?*coarse:*fine;
        const auto old=accepted;const auto oldU=old.values().normalVelocity(),oldP=old.values().pressure();const auto oldGeometry=MaxwellCoupledContract::sourceSnapshot(old.phaseGeometry().vofPublication());
        const auto targetDomain=domain(targetMesh,rz);const auto targetLayout=HydroDofLayout::build(targetMesh,targetDomain);
        const auto means=HydroVofMeshTransfer::boundaryNormalMeans(sourceMesh,old,targetMesh,targetDomain);auto targetFlow=flowFor(targetMesh,targetLayout,outlet,means);
        const auto volumes=GeometricVofMeshTransfer::prepare(sourceMesh,old.phaseGeometry().vofPublication(),targetMesh,targetDomain,old.timeLayer());
        check(cycle?volumes.report().coarsenedTargets>0:volumes.report().refinedTargets>0,"F5 true refine/coarsen coverage missing");
        const auto remap=HydroVofMeshTransfer::preparePrescribedFlow(sourceMesh,old,targetMesh,volumes,2,1,targetFlow,transferConfig);
        verify(sourceMesh,old,targetMesh,remap,outlet);electric(sourceMesh,old,targetMesh,remap);
        const auto& target=remap.candidateState();const auto& tl=target.values().layout();const auto saved=target.values().normalVelocity();
        reject([&]{HydroVofMeshTransfer::prepare(sourceMesh,old,targetMesh,volumes,2,1,transferConfig);},"F5 closed entry consumed open source");
        auto bad=targetFlow;const auto found=std::find_if(bad.prescribedVelocityDofs.begin(),bad.prescribedVelocityDofs.end(),[&](int u){const auto f=tl.domain().face(targetMesh,tl.nativeFace(u));return f.outwardSign*bad.normalVelocity[u]<0;});
        check(found!=bad.prescribedVelocityDofs.end(),"F5 target has no true incoming prescribed port");const int inlet=*found;bad.normalVelocity[inlet]*=1.001;
        reject([&]{HydroVofMeshTransfer::preparePrescribedFlow(sourceMesh,old,targetMesh,volumes,2,1,bad,transferConfig);},"F5 changed actual covered inlet g");
        bad=targetFlow;bad.rulesByNativeRow[tl.nativeFace(inlet)]=Mode::PureGasReservoir;
        reject([&]{HydroVofMeshTransfer::preparePrescribedFlow(sourceMesh,old,targetMesh,volumes,2,1,bad,transferConfig);},"F5 changed actual source reservoir phase");
        bad=targetFlow;const auto query=bad.endpointVelocity;bad.endpointVelocity=[query](int row,Clock t){auto q=query(row,t);q.y+=.01;return q;};
        reject([&]{HydroVofMeshTransfer::preparePrescribedFlow(sourceMesh,old,targetMesh,volumes,2,1,bad,transferConfig);},"F5 changed sameEnd source tangent");
        const auto& oldPub=old.phaseGeometry().vofPublication();std::vector<std::array<long double,2>> inventory(oldPub.plicCells().size());
        for(std::size_t ci=0;ci<inventory.size();++ci)inventory[ci]=oldPub.plicCells()[ci].inventory;
        const auto foreign=GeometricVofPublication::prepareLocalPlicInventory(sourceMesh,oldPub.domain(),old.timeLayer(),inventory,{});
        const auto foreignTransfer=GeometricVofMeshTransfer::prepare(sourceMesh,foreign,targetMesh,targetDomain,old.timeLayer());
        reject([&]{HydroVofMeshTransfer::preparePrescribedFlow(sourceMesh,old,targetMesh,foreignTransfer,2,1,targetFlow,transferConfig);},"F5 foreign same-clock VOF source");
        auto wrong=old.timeLayer();++wrong.step;
        reject([&]{GeometricVofMeshTransfer::prepare(sourceMesh,oldPub,targetMesh,targetDomain,wrong);},"F5 changed complete sameEnd clock");
        if(outlet){auto reversed=target.values().normalVelocity();reversed[targetFlow.pressureOutlets.front().velocityDof]*=-1;
            reject([&]{HydroTwoPhaseState::createVofPrescribedFlowRemapped(targetMesh,target.phaseGeometry(),HydroNativeValues(tl,reversed,target.values().pressure()),*target.retainedMass(),targetFlow);},"F5 target pressure outlet backflow");}
        // A physical outlet reference shift is a legal accepted state; the
        // all-fixed representation instead requires its canonical zero mean.
        constexpr double shift=2;auto shiftedFlow=*old.vofPrescribedFlow();for(auto& p:shiftedFlow.pressureOutlets)p.pressure+=shift;
        auto shiftedP=oldP;for(auto& p:shiftedP)p+=shift;
        if(outlet){
            const auto shiftedSource=HydroTwoPhaseState::createVofPrescribedFlowRemapped(sourceMesh,old.phaseGeometry(),HydroNativeValues(old.values().layout(),oldU,shiftedP),*old.retainedMass(),shiftedFlow);
            auto shiftedTargetFlow=targetFlow;for(auto& p:shiftedTargetFlow.pressureOutlets)p.pressure+=shift;
            const auto shifted=HydroVofMeshTransfer::preparePrescribedFlow(sourceMesh,shiftedSource,targetMesh,volumes,2,1,shiftedTargetFlow,transferConfig);
            check(shifted.candidateState().values().normalVelocity()==saved&&shifted.projectionMultiplier()==remap.projectionMultiplier(),"F5 physical pressure offset contaminated U/lambda projection");
            for(std::size_t p=0;p<target.values().pressure().size();++p){const auto delta=shifted.candidateState().values().pressure()[p]-target.values().pressure()[p];
                bound(delta-(outlet?shift:0),shift+std::abs(target.values().pressure()[p])+std::abs(shifted.candidateState().values().pressure()[p]),1e-12L,"F5 physical mapped pressure reference shift/gauge");}
            ++pressureShiftRemaps;
        }else{
            bool rejectedGauge=false;
            try{HydroTwoPhaseState::createVofPrescribedFlowRemapped(sourceMesh,old.phaseGeometry(),HydroNativeValues(old.values().layout(),oldU,shiftedP),*old.retainedMass(),shiftedFlow);}
            catch(const std::invalid_argument& e){rejectedGauge=std::string(e.what()).find("accepted pressure gauge is not mean zero")!=std::string::npos;}
            check(rejectedGauge,"F5 all-fixed accepted offset bypassed the original zero-mean pressure gauge");++guardCount;++gaugeRejections;
            check(old.sameIdentity(accepted)&&old.values().normalVelocity()==oldU&&old.values().pressure()==oldP
                &&MaxwellCoupledContract::sourceSnapshot(oldPub)==oldGeometry,"F5 rejected pressure gauge changed the actual accepted source");
        }
        const auto nextClock=Clock{old.timeLayer().step+1,old.timeLayer().time(),old.timeLayer().time()+dt,Clock::Position::End};
        const auto next=HydroVofStep::preparePrescribedFlow(targetMesh,target,nextClock,config);physicalStep(next,target,targetMesh);
        check(!next.report().usedVelocityHistory&&next.report().startupCorrection&&next.candidateState().velocityHistory(),"F5 true next F3 after reset did not restart history");
        check(old.sameIdentity(accepted)&&old.values().normalVelocity()==oldU&&old.values().pressure()==oldP
            &&MaxwellCoupledContract::sourceSnapshot(oldPub)==oldGeometry&&target.values().normalVelocity()==saved,"F5 remap/failures/nextStep changed accepted source or target");
        accepted=target;next.commit(targetMesh,accepted);
    }
}
void run() {
    for(bool rz:{false,true})for(bool outlet:{false,true})scene(rz,outlet);
    check(sceneCount==4&&remapCount==8&&stepCount==12&&electricCount==8,"F5 minimal real consumers incomplete");
    check(pressureShiftRemaps==4&&gaugeRejections==4,"F5 physical-outlet shift and all-fixed gauge witnesses incomplete");
    check(maxCorrection>0&&maxBoundaryWork>0,"F5 must exercise nonzero genuine projection correction/reaction; do not alter inlet g to manufacture work");
    std::cout<<std::setprecision(17)<<"VOF_OPEN_REMAP_PASS,scenes,"<<sceneCount<<",remaps,"<<remapCount<<",real_F3_steps,"<<stepCount
        <<",pressure_shift_remaps,"<<pressureShiftRemaps<<",all_fixed_gauge_rejections,"<<gaugeRejections<<",guards,"<<guardCount<<",E2_D0_sameEnd,"<<electricCount<<",energy_relative,"<<peakEnergy<<",div_relative,"<<peakDiv
        <<",maximum_WB,"<<maxBoundaryWork<<",maximum_Ecor,"<<maxCorrection<<",true_XY_RZ_2to1,fixed_and_free_outlet,"
        <<"physical_p_shift_only,history_reset,next_actual_F3,not_App_not_F4_advance_not_HF_contact\n";
}
} // namespace VofOpenRemapContract

namespace VofEhdBoundaryOutputContract {
namespace E=LAMP::Electric;
namespace fs=std::filesystem;
namespace Old=VofEhdOutputContract;
namespace Open=VofOpenRemapContract;
using Writer=VofElectrohydrodynamicsOutput;
unsigned verified=0,guards=0;
struct Table {
    std::vector<std::string> names;std::vector<std::vector<long double>> values;
    long double at(std::size_t row,const char* key)const {
        const auto i=std::find(names.begin(),names.end(),key);check(i!=names.end(),"F6C named table column missing");
        const auto k=static_cast<std::size_t>(i-names.begin());check(row<values.size()&&values[row].size()==2*names.size(),"F6C all-numeric table row/availability size");
        check(values[row][names.size()+k]==1,"F6C actual receipt value incorrectly unavailable");return values[row][k];
    }
};
Table table(const fs::path& path) {
    Table result;const auto text=Old::read(path);result.values=Old::rows(path);
    for(std::size_t i=0;;++i){const auto key="DATASETAUXDATA Column"+std::to_string(i)+"Name=\"";const auto start=text.find(key);
        if(start==std::string::npos)break;const auto first=start+key.size(),last=text.find('"',first);check(last!=std::string::npos,"F6C column metadata truncated");
        result.names.push_back(text.substr(first,last-first));}
    check(!result.names.empty(),"F6C dynamic TableHistory schema missing");return result;
}
long double rendered(long double value) {std::ostringstream out;out<<std::scientific<<std::setprecision(16)<<value;return std::stold(out.str());}
void equal(const Table& t,std::size_t row,const char* key,long double value) {
    check(t.at(row,key)==rendered(value),"F6C retained report value did not round-trip in its named column");++verified;
}
std::vector<double> inflow(const QuadtreeMesh& mesh,const GeometricVofTransport& volume) {
    std::vector<double> rho(volume.source().domain().nativeFaceCount());
    for(const auto& f:volume.faces()){if(f.negativeCell>=0&&f.positiveCell>=0)continue;
        const auto actual=volume.source().domain().face(mesh,f.nativeRow);if(actual.outwardSign*f.full.total<0)rho[f.nativeRow]=-.0625;}
    return rho;
}
void scene(const fs::path& work,bool rz) {
    auto mesh=Open::meshForOpen(false,rz),fine=Open::meshForOpen(true,rz);
    const auto pub=GeometricVofPublication::prepareLocalPlic(*mesh,domain(*mesh,rz),{},std::vector<double>(mesh->cellIndexNum(),0));
    const auto phase=phaseFor(*mesh,pub);const auto& l=phase.layout();std::vector<double> u(l.velocityCount()),normal(l.domain().nativeFaceCount());
    for(int k=0;k<l.velocityCount();++k){const auto f=l.domain().face(*mesh,l.nativeFace(k));if(!f.horizontal)u[k]=rz?.125/f.center.x:.125;
        const auto& edge=f.horizontal?mesh->hEdges(f.edgeLabel):mesh->vEdges(f.edgeLabel);auto a=mesh->verts(edge.edgeToVert[0]).coord(),b=mesh->verts(edge.edgeToVert[1]).coord();
        if(f.horizontal?a.x>b.x:a.y>b.y)std::swap(a,b);if(rz){a.x-=1;b.x-=1;}
        const auto flux=.03125L*(f.horizontal?-1.L:1.L)*(rz?CartesianMeasure::twoPi:1.L)*(psi(b)-psi(a));
        if(f.area>0)u[k]+=static_cast<double>(flux/f.area);normal[f.nativeRow]=u[k];}
    auto flow=Open::flowFor(*mesh,l,true,normal);
    const auto initial=HydroTwoPhaseState::createVofPrescribedFlow(*mesh,phase,HydroNativeValues(l,u,std::vector<double>(l.pressureCount(),.125)),flow);
    std::vector<double> qi(pub.cells().size());for(std::size_t ci=0;ci<qi.size();++ci)qi[ci]=static_cast<double>((ci%2?.125L:-.0625L)*pub.cells()[ci].capacity);
    const auto q0=E::ElectricVofCharge::fromIntegrated(*mesh,pub,qi);
    VofElectrohydrodynamicsInput input;constexpr double dt=1./16384;input.fluid.radial=rz;
    input.fluid.run.endTime=2*dt;input.fluid.run.writeFields=input.fluid.run.writeInitial=input.fluid.run.writeFirst=input.fluid.run.writeFinal=true;input.fluid.run.writeEverySteps=1;
    input.fluid.vof.liquidDensity=2;input.fluid.vof.gasDensity=1;input.fluid.vof.surfaceTension=0;
    input.fluid.vof.momentumDensityRepresentation=Density::CellAverageProjection;
    input.materials={2,5,.03,.07};
    const auto initField=E::ElectricVofField::prepareElectrostaticSnapshot(*mesh,pub,q0,input.materials,Open::electricBoundary(*mesh,pub));
    const auto root=work/(rz?"rz":"xy");Writer writer(input,{root});writer.begin();
    writer.accepted(*mesh,initial,q0,initField,nullptr,Writer::Clock::now(),0,false);
    check(!fs::exists(root/"output/series/charge_balance.dat")&&!fs::exists(root/"output/series/phase_boundary.dat"),"F6C Initial fabricated a physical timestep ledger");
    const auto oldU=initial.values().normalVelocity(),oldP=initial.values().pressure(),oldQ=q0.integratedCharge();const auto oldStocks=stocks(pub);
    const Clock end1{1,0,dt,Clock::Position::End};const auto step1=HydroVofStep::preparePrescribedFlow(*mesh,initial,end1,input.fluid.vof);
    const auto bc1=E::ElectricVofChargeBoundary::prepare(*mesh,q0,step1.transport(),inflow(*mesh,step1.transport()));
    const auto moved1=E::ElectricVofChargeTransport::prepareBoundary(*mesh,q0,step1.transport(),bc1);
    const auto& physical=step1.candidateState();const auto& physicalPub=physical.phaseGeometry().vofPublication();
    const auto field1=E::ElectricVofField::prepareOhmic(*mesh,physicalPub,moved1.endpoint(),input.materials,Open::electricBoundary(*mesh,physicalPub));
    const auto targetDomain=domain(*fine,rz);const auto targetLayout=HydroDofLayout::build(*fine,targetDomain);
    const auto means=HydroVofMeshTransfer::boundaryNormalMeans(*mesh,physical,*fine,targetDomain);
    const auto targetFlow=Open::flowFor(*fine,targetLayout,true,means);
    const auto volume=GeometricVofMeshTransfer::prepare(*mesh,physicalPub,*fine,targetDomain,end1);
    HydroVofMeshTransferConfig cfg;cfg.momentumDensityRepresentation=Density::CellAverageProjection;
    const auto remap=HydroVofMeshTransfer::preparePrescribedFlow(*mesh,physical,*fine,volume,2,1,targetFlow,cfg);
    const auto& target=remap.candidateState();const auto& actual=remap.interfaceTransfer();
    const auto qremap=E::ElectricVofChargeMeshTransfer::prepare(*mesh,physicalPub,field1.chargeCandidate(),*fine,actual);
    const auto snapshot=E::ElectricVofField::prepareElectrostaticSnapshot(*fine,target.phaseGeometry().vofPublication(),qremap.chargeCandidate(),input.materials,
        Open::electricBoundary(*fine,target.phaseGeometry().vofPublication()));
    HydrodynamicsVofRegridReport decision;decision.sourceClock=end1;decision.sourceLeaves=mesh->cellIndexNum();decision.targetLeaves=fine->cellIndexNum();decision.requested=decision.changed=true;
    const HydrodynamicsVofRemapReceipt receipt{mesh.get(),&physical,&step1,&actual,&remap,&decision};
    const auto prefix=Old::read(root/"output/meta/frames.csv");
    const auto rejected=[&](auto action,const char* why){rejects(action,why);++guards;check(Old::read(root/"output/meta/frames.csv")==prefix,"F6C bad receipt changed complete frame prefix");};
    rejected([&]{writer.accepted(*fine,target,qremap.chargeCandidate(),snapshot,&step1,Writer::Clock::now(),0,false,&receipt,&qremap,&moved1,nullptr);},"F6C unpaired physical receipt accepted");
    rejected([&]{writer.accepted(*fine,target,qremap.chargeCandidate(),snapshot,&step1,Writer::Clock::now(),0,false,&receipt,&qremap,&moved1,&snapshot);},"F6C target snapshot accepted as old-tree physical End Ohmic");
    const auto foreign=E::ElectricVofChargeTransport::prepareBoundary(*mesh,q0,step1.transport(),bc1);
    rejected([&]{writer.accepted(*fine,target,qremap.chargeCandidate(),snapshot,&step1,Writer::Clock::now(),0,false,&receipt,&qremap,&foreign,&field1);},"F6C equivalent new advection receipt accepted in place of actual field input");
    writer.accepted(*fine,target,qremap.chargeCandidate(),snapshot,&step1,Writer::Clock::now(),0,false,&receipt,&qremap,&moved1,&field1);
    const Clock end2{2,dt,2*dt,Clock::Position::End};const auto step2=HydroVofStep::preparePrescribedFlow(*fine,target,end2,input.fluid.vof);
    const auto bc2=E::ElectricVofChargeBoundary::prepare(*fine,qremap.chargeCandidate(),step2.transport(),inflow(*fine,step2.transport()));
    const auto moved2=E::ElectricVofChargeTransport::prepareBoundary(*fine,qremap.chargeCandidate(),step2.transport(),bc2);
    const auto& last=step2.candidateState();const auto field2=E::ElectricVofField::prepareOhmic(*fine,last.phaseGeometry().vofPublication(),moved2.endpoint(),input.materials,
        Open::electricBoundary(*fine,last.phaseGeometry().vofPublication()));
    writer.accepted(*fine,last,field2.chargeCandidate(),field2,&step2,Writer::Clock::now(),0,true,nullptr,nullptr,&moved2,&field2);writer.complete();
    const auto charges=table(root/"output/series/charge_balance.dat"),phases=table(root/"output/series/phase_boundary.dat"),works=table(root/"output/series/vof_boundary_work.dat"),projections=table(root/"output/series/regrid_projection_work.dat");
    check(charges.values.size()==2&&phases.values.size()==2&&works.values.size()==2&&projections.values.size()==1,"F6C real End/phase/work rows or unique changed projection row missing");
    check(!fs::exists(root/"output/series/electric_ohmic_finite.dat")&&charges.names.size()==33
        &&charges.names.back()=="dt_boundary_current_from_report_total","F6C retired Field audit schema remains or physical observations disappeared");
    for(const auto* name:{"ohmic_charge_residual","ohmic_relative_residual","combined_charge_residual","combined_relative_residual"})
        check(std::find(charges.names.begin(),charges.names.end(),name)==charges.names.end(),"F6C retired audit column was fabricated");
    const std::array<const HydroVofStep*,2> steps{{&step1,&step2}};
    const std::array<const E::ElectricVofChargeTransport*,2> moves{{&moved1,&moved2}};const std::array<const E::ElectricVofField*,2> fields{{&field1,&field2}};
    for(std::size_t k=0;k<2;++k){const auto& ar=moves[k]->endpoint().report();const auto& pl=steps[k]->transport().report().boundaryFull;const auto& w=steps[k]->report().stages.back();const auto t=steps[k]->candidateState().timeLayer();
        equal(charges,k,"source_charge_identity",static_cast<long double>(moves[k]->sourceCharge().identity()));equal(charges,k,"physical_end_charge_identity",static_cast<long double>(fields[k]->chargeCandidate().identity()));
        equal(charges,k,"source_geometry_identity",static_cast<long double>(steps[k]->beginState().phaseGeometry().vofIdentity().value()));equal(charges,k,"physical_end_geometry_identity",static_cast<long double>(fields[k]->geometry().identity().value()));
        equal(charges,k,"accepted_geometry_identity",static_cast<long double>(k?last.phaseGeometry().vofIdentity().value():target.phaseGeometry().vofIdentity().value()));
        equal(charges,k,"begin",t.begin);equal(charges,k,"end",t.end);equal(charges,k,"position",static_cast<int>(t.position));
        equal(charges,k,"advective_charge_in",ar.boundaryIncomingCharge);equal(charges,k,"advective_charge_out",ar.boundaryOutgoingCharge);equal(charges,k,"advective_charge_net_out",ar.boundaryNetCharge);
        LAMP::Geometry::Detail::FiniteVolume::Dyadic currentSum;
        for(const auto& face:fields[k]->faces())if((face.negativeCell>=0)!=(face.positiveCell>=0))
            currentSum.add(face.current,face.negativeCell>=0?1:-1);
        const long double boundaryCurrent=static_cast<double>(currentSum.rounded());
        equal(charges,k,"boundary_current",boundaryCurrent);equal(charges,k,"dt_boundary_current_from_report_total",static_cast<long double>(fields[k]->effectiveDt())*boundaryCurrent);
        const auto total=[](const std::vector<double>& q){LAMP::Geometry::Detail::FiniteVolume::Dyadic sum;for(double value:q)sum.add(value);return static_cast<double>(sum.rounded());};
        const long double sourceTotal=total(fields[k]->sourceCharge().integratedCharge()),physicalTotal=total(fields[k]->chargeCandidate().integratedCharge());
        equal(charges,k,"source_total_charge",sourceTotal);equal(charges,k,"advected_total_charge",total(fields[k]->implicitInputCharge()));
        equal(charges,k,"physical_end_total_charge",physicalTotal);equal(charges,k,"accepted_total_charge",k?physicalTotal:total(qremap.chargeCandidate().integratedCharge()));
        equal(charges,k,"physical_delta_charge",physicalTotal-sourceTotal);
        equal(charges,k,"advection_rounding_residual",ar.roundingResidual);equal(charges,k,"inlet_product_rounding_residual",ar.boundaryProductRoundingResidual);
        equal(phases,k,"inflow_liquid",pl.inflow[0]);equal(phases,k,"outflow_gas",pl.outflow[1]);equal(phases,k,"physical_residual_liquid",pl.physicalResidual[0]);
        equal(phases,k,"WY_compression_gas",pl.compression[1]);equal(phases,k,"rounding_gas",pl.roundingError[1]);equal(phases,k,"algebraic_residual_gas",pl.algebraicResidual[1]);
        equal(works,k,"pressure_outlet_work",w.pressureOutletWork);equal(works,k,"prescribed_reaction_work",w.prescribedReactionWork);equal(works,k,"boundary_kinetic_work",w.boundaryKineticWork);
        equal(works,k,"boundary_physical_kinetic_work",w.boundaryPhysicalKineticWork);equal(works,k,"inlet_dissipation_work",w.inletDissipationWork);
        check(!steps[k]->beginState().vofPrescribedFlow()->pressureOutletBackflowModel.has_value(),"F20 legacy fixture must remain one-way");
        check(works.names.back()=="boundary_flux_rounding_work","F20 new passive field must be appended after all legacy columns");
        const auto roundingColumn=works.names.size()-1;
        check(works.values[k][roundingColumn]==0&&works.values[k][works.names.size()+roundingColumn]==0,
            "F20 legacy representation work must serialize NA/Available0, not a valid zero");
        check(phases.at(k,"physical_relative_residual")<=1e-12L&&works.at(k,"closure_relative_residual")<=1e-11L,"F6C changed a retained physical gate/report");
    }
    equal(charges,0,"accepted_charge_identity",static_cast<long double>(qremap.chargeCandidate().identity()));equal(charges,0,"remapped",1);equal(charges,1,"remapped",0);
    equal(projections,0,"projection_boundary_work",remap.report().projectionBoundaryWork);equal(projections,0,"projection_rhs_rounding_work",remap.report().projectionRhsRoundingWork);
    equal(projections,0,"projection_kinetic_change",remap.report().projectionKineticChange);equal(projections,0,"source_to_target_kinetic_change",remap.report().kineticChange);
    check(moved1.endpoint().report().boundaryEnabled&&moved1.endpoint().report().boundaryIncomingCharge<0&&step1.transport().report().boundaryFull.inflow[0]>0
        &&std::any_of(field1.faces().begin(),field1.faces().end(),[](const auto& f){return f.current!=0;}),"F6C fixture lacks actual signed inlet/phase flux/nonzero Ohmic current");
    const auto index=Old::read(root/"output/meta/frames.csv");check(std::count(index.begin(),index.end(),'\n')==4&&Old::has(Old::read(root/"output/meta/manifest.json"),"\"status\":\"COMPLETED\""),"F6C complete three-frame prefix absent");
    // New-history I/O failure happens after payloads but before the only index.
    const auto failRoot=work/(rz?"rz_io_failure":"xy_io_failure");Writer failing(input,{failRoot});failing.begin();
    failing.accepted(*mesh,initial,q0,initField,nullptr,Writer::Clock::now(),0,false);const auto failurePrefix=Old::read(failRoot/"output/meta/frames.csv");
    fs::create_directory(failRoot/"output/series/charge_balance.dat");
    rejects([&]{failing.accepted(*mesh,physical,field1.chargeCandidate(),field1,&step1,Writer::Clock::now(),0,false,nullptr,nullptr,&moved1,&field1);},"F6C failed boundary-history payload accepted");++guards;
    check(Old::read(failRoot/"output/meta/frames.csv")==failurePrefix,"F6C new history I/O failure published a partial frame");
    rejects([&]{failing.complete();},"F6C poisoned new-history output completed");++guards;
    check(initial.values().normalVelocity()==oldU&&initial.values().pressure()==oldP&&q0.integratedCharge()==oldQ&&stocks(pub)==oldStocks,"F6C passive output changed source inventory/U/p/Q");
}
void run(const fs::path& work) {
    check(!fs::exists(work),"F6C output requires unique fresh D work directory");fs::create_directories(work);scene(work,false);scene(work,true);
    std::cout<<"VOF_EHD_BOUNDARY_OUTPUT_PASS,scenes,2,accepted,4,changed,2,frames,6,named_report_values,"<<verified<<",guards,"<<guards
        <<",actual_F3_F4_EndOhmic_F5_Qremap_D0,physical_vs_target_source,TableHistory_dynamic_schema,payload_flush_index,no_App_claim\n";
}
} // namespace VofEhdBoundaryOutputContract

}
// BEGIN F20 EXPLICIT BOUNDARY REPRESENTATION WORK OUTPUT CONTRACT
namespace F20BoundaryRoundingOutputContract {
namespace E=LAMP::Electric;
namespace fs=std::filesystem;
namespace Old=VofEhdOutputContract;
namespace Open=VofOpenRemapContract;
namespace BoundaryOutput=VofEhdBoundaryOutputContract;
using Writer=VofElectrohydrodynamicsOutput;
void run(const fs::path& work) {
    check(!fs::exists(work),"F20 output requires a unique fresh D work directory");fs::create_directories(work);
    auto mesh=Open::meshForOpen(false,false);
    const auto pub=GeometricVofPublication::prepareLocalPlic(*mesh,domain(*mesh,false),{},std::vector<double>(mesh->cellIndexNum(),0));
    const auto phase=phaseFor(*mesh,pub);const auto& layout=phase.layout();
    check(!layout.domain().uniformLeaves(),"F20 requires genuine native 2:1 geometry");
    std::vector<double> u(layout.velocityCount()),normal(layout.domain().nativeFaceCount());
    for(int k=0;k<layout.velocityCount();++k){const auto f=layout.domain().face(*mesh,layout.nativeFace(k));
        if(!f.horizontal)u[k]=.125;normal[f.nativeRow]=u[k];}
    auto flow=Open::flowFor(*mesh,layout,true,normal);
    // The opt-in F16 scope requires steady prescribed full vectors. Freeze
    // actual native-row data by value; never borrow the helper's tree/domain.
    std::vector<Vector2D> steadyEndpoints(normal.size());
    for(std::size_t row=0;row<steadyEndpoints.size();++row){const auto mode=flow.rulesByNativeRow[row];
        if(mode==GeometricVofAdvectorSlab::Boundary::Mode::PureLiquidReservoir||
            mode==GeometricVofAdvectorSlab::Boundary::Mode::PureGasReservoir)
            steadyEndpoints[row]=Vector2D(normal[row],.03125);}
    flow.endpointVelocity=[steadyEndpoints](int row,Clock endpoint){endpoint.time();return steadyEndpoints.at(row);};
    flow.pressureOutletBackflowModel=HydroVofPressureOutletMomentumModel::Ambient;
    const auto queryClocks=std::make_shared<std::vector<Clock>>();
    flow.pressureOutletBackflowMomentum=[queryClocks](int row,Clock fullMid){check(row>=0,"F20 real ambient native row");
        queryClocks->push_back(fullMid);return Vector2D(.03125,.09375);};
    const auto initial=HydroTwoPhaseState::createVofPrescribedFlow(*mesh,phase,
        HydroNativeValues(layout,u,std::vector<double>(layout.pressureCount(),.125)),flow);
    check(initial.vofPrescribedFlow()&&initial.vofPrescribedFlow()->pressureOutletBackflowModel.has_value()
        &&bool(initial.vofPrescribedFlow()->pressureOutletBackflowMomentum),"F20 requires actual F17 opt-in source/model/callback");
    std::vector<double> qi(pub.cells().size());for(std::size_t ci=0;ci<qi.size();++ci)qi[ci]=static_cast<double>((ci%2?.125L:-.0625L)*pub.cells()[ci].capacity);
    const auto charge=E::ElectricVofCharge::fromIntegrated(*mesh,pub,qi);
    VofElectrohydrodynamicsInput input;constexpr double dt=1./16384;
    input.fluid.run.endTime=dt;input.fluid.run.writeFields=input.fluid.run.writeInitial=input.fluid.run.writeFirst=input.fluid.run.writeFinal=true;
    input.fluid.run.writeEverySteps=1;input.fluid.vof.liquidDensity=2;input.fluid.vof.gasDensity=1;input.fluid.vof.surfaceTension=0;
    input.fluid.vof.momentumDensityRepresentation=Density::CellAverageProjection;input.materials={2,5,.03,.07};
    const auto initialField=E::ElectricVofField::prepareElectrostaticSnapshot(*mesh,pub,charge,input.materials,Open::electricBoundary(*mesh,pub));
    const auto sourceU=initial.values().normalVelocity(),sourceP=initial.values().pressure(),sourceQ=charge.integratedCharge();const auto sourceStocks=stocks(pub);
    Writer writer(input,{work/"optin"});writer.begin();writer.accepted(*mesh,initial,charge,initialField,nullptr,Writer::Clock::now(),0,false);
    check(!fs::exists(work/"optin/output/series/vof_boundary_work.dat"),"F20 Initial must not fabricate final-stage work");
    const Clock end{1,0,dt,Clock::Position::End};const auto step=HydroVofStep::preparePrescribedFlow(*mesh,initial,end,input.fluid.vof);
    check(!queryClocks->empty(),"F20 actual stage never consumed the explicit ambient callback");
    for(const auto clock:*queryClocks)check(clock.step==end.step&&clock.begin==end.begin&&clock.end==end.end
        &&clock.position==Clock::Position::Midpoint,"F20 ambient callback lacks full actual stage Midpoint clock");
    const auto boundary=E::ElectricVofChargeBoundary::prepare(*mesh,charge,step.transport(),BoundaryOutput::inflow(*mesh,step.transport()));
    const auto moved=E::ElectricVofChargeTransport::prepareBoundary(*mesh,charge,step.transport(),boundary);
    const auto& candidate=step.candidateState();const auto& candidatePub=candidate.phaseGeometry().vofPublication();
    const auto field=E::ElectricVofField::prepareOhmic(*mesh,candidatePub,moved.endpoint(),input.materials,Open::electricBoundary(*mesh,candidatePub));
    writer.accepted(*mesh,candidate,field.chargeCandidate(),field,&step,Writer::Clock::now(),0,true,nullptr,nullptr,&moved,&field);writer.complete();
    const auto table=BoundaryOutput::table(work/"optin/output/series/vof_boundary_work.dat");
    check(table.values.size()==1&&table.names.back()=="boundary_flux_rounding_work","F20 one real final-stage row and exact appended schema column");
    const auto& actual=step.report().stages.back();
    BoundaryOutput::equal(table,0,"boundary_flux_rounding_work",actual.boundaryFluxRoundingWork);
    BoundaryOutput::equal(table,0,"boundary_kinetic_work",actual.boundaryKineticWork);
    BoundaryOutput::equal(table,0,"boundary_physical_kinetic_work",actual.boundaryPhysicalKineticWork);
    BoundaryOutput::equal(table,0,"inlet_dissipation_work",actual.inletDissipationWork);
    const auto scale=std::abs(actual.boundaryKineticWork)+std::abs(actual.boundaryPhysicalKineticWork)
        +std::abs(actual.inletDissipationWork)+std::abs(actual.boundaryFluxRoundingWork);
    Open::bound(actual.boundaryKineticWork-actual.boundaryPhysicalKineticWork-actual.inletDissipationWork-actual.boundaryFluxRoundingWork,
        scale,1e-11L,"F20 passive P+D+Wrep split; Wrep must not be counted twice");
    check(table.at(0,"closure_relative_residual")<=1e-11L,"F20 writer changed actual final-stage closure gate");
    const auto index=Old::read(work/"optin/output/meta/frames.csv");
    check(std::count(index.begin(),index.end(),'\n')==3,"F20 only complete Initial/End frames may be published");
    check(initial.values().normalVelocity()==sourceU&&initial.values().pressure()==sourceP&&charge.integratedCharge()==sourceQ
        &&stocks(pub)==sourceStocks,"F20 passive serialization changed actual source stocks/U/p/Q");
    std::cout<<std::setprecision(17)<<"VOF_EHD_BOUNDARY_ROUNDING_OUTPUT_PASS,scenes,1,accepted,1,actual_final_Wrep,"<<actual.boundaryFluxRoundingWork
        <<",available,1,callback_queries,"<<queryClocks->size()<<",optin_zero_is_available,no_F19_remap_dependency,no_R8_App_claim\n";
}
} // namespace F20BoundaryRoundingOutputContract
// END F20 EXPLICIT BOUNDARY REPRESENTATION WORK OUTPUT CONTRACT

// BEGIN F19 ACTUAL BACKFLOW MODEL REMAP CONSUMER
namespace F19BackflowRemapContract {
namespace FV=LAMP::Geometry::Detail::FiniteVolume;
using Model=HydroVofPressureOutletMomentumModel;
using Mode=GeometricVofAdvectorSlab::Boundary::Mode;
unsigned sceneCount=0,remapCount=0,stepCount=0,guardCount=0,initialNegative=0,initialZero=0;
unsigned targetNegative=0,postStepNegative=0,ambientSpans=0,continuationIncoming=0,ambientIncoming=0,distinctCoarsens=0;
long double peakEnergy=0,peakDiv=0,peakPhase=0,maximumWB=0,maximumWRHS=0;
void bound(long double r,long double s,long double tol,const char* why)
{VofOpenRemapContract::bound(r,s,tol,why);}
template<class F> void guard(F&& f,const char* literal) {
    bool caught=false;try{f();}catch(const std::exception& e){const auto message=std::string(e.what());
        check(message.find(literal)!=std::string::npos,(std::string("F19 wrong rejection: ")+message+" instead of "+literal).c_str());caught=true;}
    check(caught,literal);++guardCount;
}
std::unique_ptr<QuadtreeMesh> meshForBackflow(bool fine,bool rz) {
    DictionaryMesh d;d.Nx=d.Ny=4;const double left=rz?1.:0.;d.minCoord={left,0};d.maxCoord={left+1,1};
    d.maxCellNum=4096;d.maxInitialIterations=0;d.maxLevel_local={1};auto mesh=std::make_unique<QuadtreeMesh>(d);
    DictionaryMesh::RefineBox right;right.level=1;right.min={left+.5,0};right.max={left+1,1};
    std::vector<DictionaryMesh::RefineBox> boxes{right};if(fine){auto top=right;top.min={left,.5};top.max={left+.5,1};boxes.push_back(top);}
    mesh->applyRefineBoxes(boxes);return mesh;
}
struct Query {int row;Clock clock;Vector2D value;};
using QueryLog=std::shared_ptr<std::vector<Query>>;
const Vector2D ambient{.03125,.09375};
HydroVofPrescribedFlow flowForBackflow(const QuadtreeMesh& mesh,const HydroDofLayout& layout,
    const std::vector<double>& normals,Model model,const QueryLog& log) {
    const auto native=layout.domain();check(normals.size()==static_cast<std::size_t>(native.nativeFaceCount()),"F19 actual target native normal map");
    HydroVofPrescribedFlow flow;flow.rulesByNativeRow.assign(native.nativeFaceCount(),Mode::Interior);
    flow.normalVelocity.assign(layout.velocityCount(),0);flow.pressureOutletBackflowModel=model;
    std::vector<bool> fixed(native.nativeFaceCount()),free(native.nativeFaceCount());
    for(int u=0;u<layout.velocityCount();++u){const auto face=native.face(mesh,layout.nativeFace(u));if(face.outwardSign==0)continue;
        if(!face.horizontal||face.area==0){flow.rulesByNativeRow[face.nativeRow]=Mode::Closed;continue;}
        if(face.outwardSign>0){flow.rulesByNativeRow[face.nativeRow]=Mode::PureGasReservoir;flow.pressureOutlets.push_back({u,.125});free[face.nativeRow]=true;}
        else if(normals[face.nativeRow]==0)flow.rulesByNativeRow[face.nativeRow]=Mode::Closed;
        else{flow.rulesByNativeRow[face.nativeRow]=Mode::PureGasReservoir;flow.prescribedVelocityDofs.push_back(u);
            flow.normalVelocity[u]=normals[face.nativeRow];fixed[face.nativeRow]=true;}}
    // CURRENT row map/P0 data are owned by value, never a source mesh pointer.
    flow.endpointVelocity=[fixed,normals](int row,Clock t){t.time();check(row>=0&&row<static_cast<int>(fixed.size())&&fixed[row],"F19 endpoint query is not a current fixed native row");
        check(t.position==Clock::Position::Initial||t.position==Clock::Position::End,"F19 endpoint clock is not actual Initial/End");return Vector2D(.015625,normals[row]);};
    if(model==Model::Ambient){std::vector<Vector2D> values(free.size(),Vector2D{});for(std::size_t row=0;row<free.size();++row)if(free[row])values[row]=ambient;
        flow.pressureOutletBackflowMomentum=[free,values,log](int row,Clock t){t.time();
            check(row>=0&&row<static_cast<int>(free.size())&&free[row],"F19 ambient query is not a current free native row");
            check(t.position==Clock::Position::Midpoint&&t.step>0&&t.begin<t.end,"F19 ambient query lacks actual complete Midpoint interval");
            const auto value=values[row];log->push_back({row,t,value});return value;};}
    return flow;
}
struct Span {double normal,lower,upper;int axis;};
Span span(const QuadtreeMesh& mesh,const NativeLeafFace& face) {
    const auto& edge=face.horizontal?mesh.hEdges(face.edgeLabel):mesh.vEdges(face.edgeLabel);
    const auto a=mesh.verts(edge.edgeToVert[0]).coord(),b=mesh.verts(edge.edgeToVert[1]).coord();
    return face.horizontal?Span{face.center.y,std::min(a.x,b.x),std::max(a.x,b.x),1}:
        Span{face.center.x,std::min(a.y,b.y),std::max(a.y,b.y),0};
}
std::vector<int> coveredSourceRows(const QuadtreeMesh& sourceMesh,const HydroTwoPhaseState& source,
    const QuadtreeMesh& targetMesh,const HydroDofLayout& tl,int targetVelocity) {
    const auto tf=tl.domain().face(targetMesh,tl.nativeFace(targetVelocity));const auto t=span(targetMesh,tf);
    const auto& sl=source.values().layout();std::vector<int> rows;FV::Dyadic total;
    for(const auto& port:source.vofPrescribedFlow()->pressureOutlets){const auto sf=sl.domain().face(sourceMesh,sl.nativeFace(port.velocityDof));const auto s=span(sourceMesh,sf);
        if(t.axis!=s.axis||t.normal!=s.normal||tf.outwardSign!=sf.outwardSign)continue;
        const auto lo=std::max(t.lower,s.lower),hi=std::min(t.upper,s.upper);if(!(lo<hi))continue;
        FV::Dyadic length(hi);length.add(lo,-1);total.add(length);rows.push_back(sf.nativeRow);}
    FV::Dyadic whole(t.upper);whole.add(t.lower,-1);check(total.compare(whole)==0&&!rows.empty(),"F19 independent exact free-port source coverage");return rows;
}
void queryCoverage(const QuadtreeMesh& sourceMesh,const HydroTwoPhaseState& source,const QuadtreeMesh& targetMesh,
    const HydroTwoPhaseState& target,Model model,const QueryLog& sourceLog,std::size_t sourceBegin,const QueryLog& targetLog,std::size_t targetBegin,bool coarsen) {
    const auto& sl=source.values().layout();const auto& tl=target.values().layout();const auto& sf=*source.vofPrescribedFlow();const auto& tf=*target.vofPrescribedFlow();
    check(sf.pressureOutletBackflowModel==model&&tf.pressureOutletBackflowModel==model,"F19 changed target lost explicit model");
    unsigned pairs=0,multiple=0;std::vector<unsigned> oldHits(sl.domain().nativeFaceCount()),newHits(tl.domain().nativeFaceCount());
    for(const auto& port:tf.pressureOutlets){const auto row=tl.nativeFace(port.velocityDof);const auto covered=coveredSourceRows(sourceMesh,source,targetMesh,tl,port.velocityDof);
        multiple+=covered.size()>1;for(int oldRow:covered){++oldHits[oldRow];++newHits[row];++pairs;}}
    check(pairs>0,"F19 missing actual physical pressure coverage");
    if(coarsen)check(multiple>0,"F19 coarsen did not merge any actual free face");
    else check(std::any_of(oldHits.begin(),oldHits.end(),[](unsigned n){return n>1;}),"F19 refine did not split any actual free face");
    if(model==Model::InteriorContinuationZeroShear){check(!sf.pressureOutletBackflowMomentum&&!tf.pressureOutletBackflowMomentum&&sourceLog->empty()&&targetLog->empty(),
        "F19 continuation carries/froze unused ambient callback");return;}
    auto mid=source.timeLayer();check(mid.position==Clock::Position::End,"F19 remap source is not actual End");mid.position=Clock::Position::Midpoint;
    check(sourceLog->size()-sourceBegin==2*pairs&&targetLog->size()-targetBegin==2*pairs,"F19 each real overlap must repeat-query both actual owners");
    const auto verify=[&](const QueryLog& log,std::size_t start,const std::vector<unsigned>& hits){for(std::size_t k=start;k<log->size();++k){const auto& q=(*log)[k];
        check(q.row>=0&&q.row<static_cast<int>(hits.size())&&hits[q.row]>0&&sameClock(q.clock,mid),"F19 remap query used foreign row/End/fabricated clock");
        check(q.value.x==ambient.x&&q.value.y==ambient.y,"F19 split/coarsen changed finite nonzero full-vector P0 data");}};
    verify(sourceLog,sourceBegin,oldHits);verify(targetLog,targetBegin,newHits);ambientSpans+=pairs;
}
void auditProjection(const QuadtreeMesh& sourceMesh,const HydroTwoPhaseState& source,const QuadtreeMesh& targetMesh,const HydroVofMeshTransfer& remap) {
    ++remapCount;const auto& target=remap.candidateState();const auto& l=target.values().layout();const auto& flow=*target.vofPrescribedFlow();const auto& r=remap.report();
    const auto& mass=*target.retainedMass();const auto& u=target.values().normalVelocity();const auto& v=remap.rawValues().normalVelocity();
    const auto& b=remap.projectionMomentumRhs();const auto& lambda=remap.projectionMultiplier();const auto& metrics=target.phaseGeometry().metrics();
    const auto pressure=HydroPressureLoad::build(metrics.divergence());const auto& C=pressure.matrix();const auto& D=pressure.divergence().matrix();target.validateFor(targetMesh);
    check(remap.sourceState().sameIdentity(source)&&sameClock(target.timeLayer(),source.timeLayer())&&r.prescribedFlow&&r.historyReset&&!target.velocityHistory()&&r.remapSolves==1,
        "F19 real source/target same-End owner/history reset");
    check(!source.values().layout().domain().uniformLeaves()&&!l.domain().uniformLeaves()&&r.boundaryCoverageSpans>0,"F19 requires real source/target 2:1 physical coverage");
    check(target.phaseGeometry().vofIdentity()==remap.interfaceTransfer().candidate().identity()&&mass.geometry().sameIdentity(target.momentumGeometry())&&mass.representation()==Density::CellAverageProjection,
        "F19 unique target VOF/geometry/retained P0 mass owner");
    check(!r.pressureGaugeApplied&&r.pressureGaugeShift==0&&r.projectionSaddle.pressureGaugeCount==0,"F19 physical outlet pressure was re-gauged");
    check(r.maximumFinalRelativeDivergence<=1e-11&&r.projectionSaddle.algebraicRelativeResidual<=1e-11&&r.projectionSaddle.freeMomentumDefect<=1e-11
        &&r.projectionSaddle.maxCellContinuityDefect<=1e-11&&r.projectionEnergyIdentityRelativeResidual<=2e-11L,"F19 original open projection gates");
    std::vector<bool> fixed(u.size());for(int k:target.boundary().velocityDofs()){fixed[k]=true;check(u[k]==0&&v[k]==0,"F19 closed wall changed its actual zero normal");}
    for(int k:flow.prescribedVelocityDofs){fixed[k]=true;check(u[k]==flow.normalVelocity[k]&&v[k]==flow.normalVelocity[k],"F19 fixed g is not actual source polynomial mean");}
    unsigned negatives=0;for(const auto& port:flow.pressureOutlets){const auto f=l.domain().face(targetMesh,l.nativeFace(port.velocityDof));
        check(!fixed[port.velocityDof]&&flow.normalVelocity[port.velocityDof]==0,"F19 incoming momentum accidentally fixed free normal");negatives+=f.outwardSign*u[port.velocityDof]<0;}
    check(negatives>0,"F19 changed target lost actual free-port backflow witness");targetNegative+=negatives;
    FV::Dyadic delta,correction,WB,WRHS,freeWork,continuity,scale;
    for(int k=0;k<l.velocityCount();++k){const long double m=mass.masses()[k],a=u[k],raw=v[k];
        check(FV::roundDouble(FV::Dyadic::product({m,raw})).value==b[k],"F19 projection RHS must be actual RN(Mv)");
        scale.add(FV::Dyadic::product({.5L,m,a,a}).absolute());scale.add(FV::Dyadic::product({.5L,m,raw,raw}).absolute());
        FV::Dyadic dk=FV::Dyadic::product({.5L,m,a,a});dk.add(FV::Dyadic::product({.5L,m,raw,raw}),-1);delta.add(dk);scale.add(dk.absolute());
        FV::Dyadic ec=FV::Dyadic::product({.5L,m,a,a});ec.add(FV::Dyadic::product({.5L,m,raw,raw}));ec.add(FV::Dyadic::product({m,a,raw}),-1);
        check(ec.compare(FV::Dyadic{})>=0,"F19 exact M correction norm must be positive");correction.add(ec);scale.add(ec.absolute());
        FV::Dyadic wr=FV::Dyadic::product({a,b[k]});wr.add(FV::Dyadic::product({a,m,raw}),-1);WRHS.add(wr);scale.add(wr.absolute());
        FV::Dyadic row=FV::Dyadic::product({a,m,a});row.add(FV::Dyadic::product({a,b[k]}),-1);
        for(HydroPressureLoad::Matrix::InnerIterator it(C,k);it;++it)row.add(FV::Dyadic::product({a,it.value(),lambda[it.col()]}));
        if(fixed[k]){WB.add(row);scale.add(row.absolute());}else freeWork.add(row);
    }
    for(int p=0;p<l.pressureCount();++p){FV::Dyadic div;long double activity=0;
        for(HydroIntegralDivergence::Matrix::InnerIterator it(D,p);it;++it){div.add(FV::Dyadic::product({it.value(),u[it.col()]}));
            activity+=std::abs(static_cast<long double>(it.value())*u[it.col()]);continuity.add(FV::Dyadic::product({lambda[p],it.value(),u[it.col()]}));}
        bound(div.rounded(),activity,1e-11L,"F19 independent actual D u");if(activity>0)peakDiv=std::max(peakDiv,std::abs(div.rounded())/activity);}
    FV::Dyadic left=delta;left.add(correction);left.add(WB,-1);left.add(WRHS,-1);FV::Dyadic right=freeWork;right.add(continuity);
    check(left.compare(right)==0,"F19 exact C=-D^T projection identity");bound(left.rounded(),scale.rounded(),2e-11L,"F19 uncompensated independent energy/original no-floor ceiling");
    const auto compare=[&](long double observed,const FV::Dyadic& value,const char* why){bound(observed-value.rounded(),r.projectionEnergyIdentityScale,1e-12L,why);};
    compare(r.projectionKineticChange,delta,"F19 actual projection DeltaK");compare(r.projectionCorrectionEnergy,correction,"F19 actual projection correction energy");
    compare(r.projectionBoundaryWork,WB,"F19 independent fixed reaction WB");compare(r.projectionRhsRoundingWork,WRHS,"F19 independent real RHS rounding WRHS");
    compare(r.projectionFreeResidualWork,freeWork,"F19 independent free residual work");compare(r.projectionContinuityWork,continuity,"F19 independent lambda D u");
    compare(r.projectionEnergyIdentityResidual,left,"F19 original uncorrected energy receipt");
    if(scale.rounded()>0)peakEnergy=std::max(peakEnergy,std::abs(left.rounded())/scale.rounded());maximumWB=std::max(maximumWB,std::abs(WB.rounded()));maximumWRHS=std::max(maximumWRHS,std::abs(WRHS.rounded()));
    const auto& overlap=remap.interfaceTransfer().nativeOverlaps();for(int p=0;p<l.pressureCount();++p){FV::Dyadic integral,activity;
        for(auto id:overlap.byTarget()[l.nativeCell(p)]){const auto& piece=overlap.pieces()[id];const auto oldP=source.values().pressure()[source.values().layout().pressureDof(piece.sourceCell)];
            const auto term=FV::Dyadic::product({piece.capacity,oldP});integral.add(term);activity.add(term.absolute());}
        const long double volume=metrics.cellVolumes()[p];const auto expected=static_cast<double>(integral.rounded()/volume);
        bound(target.values().pressure()[p]-expected,activity.rounded()/volume,1e-12L,"F19 actual overlap-averaged physical p, not lambda");}
    const auto a=stocks(source.phaseGeometry().vofPublication()),z=stocks(target.phaseGeometry().vofPublication());
    for(int k=0;k<2;++k)bound(z[k]-a[k],std::abs(z[k])+std::abs(a[k]),1e-12L,"F19 two independent remap stocks");
    const auto sourceBudget=budget(*source.retainedMass(),source.values()),rawBudget=budget(mass,remap.rawValues()),targetBudget=budget(mass,target.values());
    const std::array<HydroVofRemapNativeBudget,3> actual{{sourceBudget,rawBudget,targetBudget}};
    const std::array<const HydroVofRemapNativeBudget*,3> reported{{&r.source,&r.rawTarget,&r.target}};
    for(int j=0;j<3;++j){for(int k=0;k<2;++k)bound(reported[j]->mass[k]-actual[j].mass[k],actual[j].mass[k],1e-12L,"F19 independent actual native component mass");
        bound(reported[j]->kinetic-actual[j].kinetic,actual[j].kinetic,1e-12L,"F19 independent actual native KE, not physical remap work");}
}
void auditStep(const QuadtreeMesh& mesh,const HydroTwoPhaseState& old,const HydroVofStep& step,Model model) {
    ++stepCount;const auto& volume=step.transport();const auto& next=step.candidateState();const auto& flow=*old.vofPrescribedFlow();const auto& l=old.values().layout();
    const auto& source=old.phaseGeometry().vofPublication();const auto& frozen=step.frozenMidpoint();const auto& dual=frozen.material();const auto& mid=step.solvedMidpointValues();
    check(step.beginState().sameIdentity(old)&&volume.sourceIdentity()==source.identity()&&sameClock(next.timeLayer(),volume.endTime())
        &&next.vofPrescribedFlow()->pressureOutletBackflowModel==model,"F19 actual target-following Step/transport/model owner");
    check(frozen.convection().isVofOpenFlow()&&frozen.appliedConvection().boundary().hasVofOpenMomentum(),"F19 actual Step silently used legacy trace KEP");
    check(frozen.viscosityMatrix().nonZeros()==0&&!step.hasCapillaryLoad()&&!frozen.hasMaxwellLoad(),"F19 fixture must remain mu=sigma=0/no unrelated force");
    for(const auto& report:step.report().stages)check(report.saddle.algebraicRelativeResidual<=1e-11&&report.saddle.freeMomentumDefect<=1e-11
        &&report.saddle.maxCellContinuityDefect<=1e-11&&report.momentumPowerRelativeDefect<=1e-11&&report.energy.closureRelativeResidual<=1e-11L,"F19 original actual Step gates");
    std::vector<bool> free(l.velocityCount());for(const auto& port:flow.pressureOutlets)free[port.velocityDof]=true;
    const auto& typed=frozen.appliedConvection().boundary();const auto& segments=dual.geometry().segments();std::vector<GeometricVofNativeSubface> queries;std::vector<std::size_t> ids;
    for(std::size_t g=0;g<segments.size();++g){const auto& s=segments[g];if(s.negativeVelocity>=0&&s.positiveVelocity>=0)continue;
        const auto f=l.domain().face(mesh,l.nativeFace(s.primalVelocity));if(f.area==0)continue;queries.push_back({f.nativeRow,{s.begin.x,s.begin.y},{s.end.x,s.end.y}});ids.push_back(g);}
    const auto measured=volume.restrictNativeFaces(mesh,queries);const long double dt=step.report().deltaTime;long double physical=0,dissipation=0,representation=0,BKE=0,activity=0;
    for(std::size_t k=0;k<ids.size();++k){const auto g=ids[k];const auto& s=segments[g];const auto& f=measured[k].full;const int sign=s.negativeVelocity>=0?1:-1;
        const int owner=s.negativeVelocity>=0?s.negativeVelocity:s.positiveVelocity,axis=l.normalAxis(owner)==HydroDofLayout::NormalAxis::First?0:1;
        const auto incoming=sign>0?std::array<long double,2>{f.negativeLiquid,f.negativeGas}:std::array<long double,2>{f.positiveLiquid,f.positiveGas};
        const auto outgoing=sign>0?std::array<long double,2>{f.positiveLiquid,f.positiveGas}:std::array<long double,2>{f.negativeLiquid,f.negativeGas};
        const long double pin=(2*incoming[0]+incoming[1])/dt,pout=(2*outgoing[0]+outgoing[1])/dt,u=mid.normalVelocity()[owner],q=sign*dual.faceRates()[g];long double b=0;
        const auto& flux=typed.vofOpenMomentumFluxes()[g];check(flux.incomingPhaseIntegral==incoming&&flux.outgoingPhaseIntegral==outgoing,"F19 Step must consume actual two-phase subface gross");
        if(pin>0){if(free[s.primalVelocity]){if(model==Model::Ambient){b=axis==0?ambient.x:ambient.y;++ambientIncoming;check(flux.incomingMomentum==b,"F19 target Step lost nonzero ambient component");}
                else{b=u;++continuationIncoming;check(flux.incomingMomentum==0,"F19 continuation froze old U into affine incoming data");}}
            else b=axis==0?.015625L:flow.normalVelocity[s.primalVelocity];}
        const auto n=pout*u-pin*b,p=.5L*pout*u*u-.5L*pin*b*b,d=.5L*pin*(u-b)*(u-b),w=.5L*(pout-pin-q)*u*u;
        physical+=p;dissipation+=d;representation+=w;BKE+=u*n-.5L*q*u*u;activity+=std::abs(p)+std::abs(d)+std::abs(w)+std::abs(u*n)+std::abs(q*u*u);}
    const auto& report=step.report().stages.back();bound(report.boundaryPhysicalKineticWork-dt*physical,dt*activity,1e-12L,"F19 actual physical boundary KE");
    bound(report.inletDissipationWork-dt*dissipation,dt*activity,1e-12L,"F19 actual inlet dissipation");bound(report.boundaryFluxRoundingWork-dt*representation,dt*activity,1e-12L,"F19 passive Wrep from actual current U");
    bound(report.boundaryKineticWork-dt*BKE,dt*activity,1e-12L,"F19 P+D+Wrep in BKE once");bound(BKE-physical-dissipation-representation,activity,1e-12L,"F19 independent boundary power split");
    for(const auto& port:flow.pressureOutlets){const auto f=l.domain().face(mesh,l.nativeFace(port.velocityDof));const long double traction=-static_cast<long double>(f.outwardSign)*f.area*port.pressure;
        bound(frozen.saddle().pressureOutletLoads()[port.velocityDof]-traction,std::abs(traction),1e-12L,"F19 actual pressure traction survives remap");
        postStepNegative+=f.outwardSign*next.values().normalVelocity()[port.velocityDof]<0;}
    for(bool half:{true,false}){const auto& endpoint=half?volume.midpoint():volume.endpoint();const auto& ledger=half?volume.report().boundaryHalf:volume.report().boundaryFull;
        long double capacity=0;for(const auto& c:source.cells())capacity+=c.capacity;
        for(int phase=0;phase<2;++phase){FV::Dyadic change,out,in;
            for(std::size_t ci=0;ci<source.plicCells().size();++ci){change.add(endpoint.plicCells()[ci].inventory[phase]);change.add(source.plicCells()[ci].inventory[phase],-1);}
            for(const auto& face:volume.faces())if(face.negativeCell<0||face.positiveCell<0){const auto& q=half?face.half:face.full;
                const auto positive=phase==0?q.positiveLiquid:q.positiveGas,negative=phase==0?q.negativeLiquid:q.negativeGas;out.add(face.negativeCell>=0?positive:negative);in.add(face.negativeCell>=0?negative:positive);}
            change.add(out);change.add(in,-1);bound(change.rounded(),capacity,1e-12L,"F19 actual independent physical phase Delta+Out-In");
            bound(ledger.inflow[phase]-in.rounded(),capacity,1e-12L,"F19 actual gross phase inflow receipt");bound(ledger.outflow[phase]-out.rounded(),capacity,1e-12L,"F19 actual gross phase outflow receipt");
            if(capacity>0)peakPhase=std::max(peakPhase,std::abs(change.rounded())/capacity);}
        check(ledger.enabled&&ledger.sourceIdentity==source.identity()&&ledger.physicalRelativeConservation<=1e-12L&&ledger.algebraicRelativeConservation<=1e-12L,"F19 original physical/WY/rounding gates");}
}
void scene(bool rz,Model model) {
    ++sceneCount;auto coarse=meshForBackflow(false,rz),fine=meshForBackflow(true,rz);const auto pub=GeometricVofPublication::prepareLocalPlic(*coarse,domain(*coarse,rz),{},std::vector<double>(coarse->cellIndexNum(),0));
    const auto phase=phaseFor(*coarse,pub);const auto& l=phase.layout();const auto native=l.domain();check(coarse->cellIndexNum()==40&&!native.uniformLeaves(),"F19 genuine 2:1 initial source");
    std::vector<double> u(l.velocityCount()),normals(native.nativeFaceCount());const double left=rz?1.:0.;
    for(int k=0;k<l.velocityCount();++k){const auto f=native.face(*coarse,l.nativeFace(k));const auto x=f.center.x-left;
        if(f.horizontal)u[k]=x<.25?-.0625:x<.5?0.:.0625;normals[f.nativeRow]=u[k];}
    auto sourceLog=std::make_shared<std::vector<Query>>();auto flow=flowForBackflow(*coarse,l,normals,model,sourceLog);const HydroNativeValues initialValues(l,u,std::vector<double>(l.pressureCount(),.125));
    const auto div=phase.metrics().divergence().apply(*coarse,initialValues);const auto& D=phase.metrics().divergence().matrix();
    for(int p=0;p<l.pressureCount();++p){long double activity=0;for(HydroIntegralDivergence::Matrix::InnerIterator it(D,p);it;++it)activity+=std::abs(static_cast<long double>(it.value())*u[it.col()]);bound(div[p],activity,1e-12L,"F19 actual divergence-free Initial source");}
    auto accepted=HydroTwoPhaseState::createVofPrescribedFlow(*coarse,phase,initialValues,flow);unsigned negative=0,zero=0;
    for(const auto& port:accepted.vofPrescribedFlow()->pressureOutlets){const auto f=native.face(*coarse,l.nativeFace(port.velocityDof));negative+=f.outwardSign*accepted.values().normalVelocity()[port.velocityDof]<0;zero+=accepted.values().normalVelocity()[port.velocityDof]==0;}
    check(negative>0&&zero>0,"F19 actual accepted Initial must witness negative and exact-zero free normal");initialNegative+=negative;initialZero+=zero;
    auto unused=flow;unused.pressureOutletBackflowModel.reset();unused.pressureOutletBackflowMomentum=[](int,Clock){return ambient;};
    guard([&]{HydroTwoPhaseState::createVofPrescribedFlow(*coarse,phase,initialValues,unused);},"VOF backflow callback requires an explicit momentum model");
    HydroVofAcConfig config;config.liquidDensity=2;config.gasDensity=1;config.momentumDensityRepresentation=Density::CellAverageProjection;config.surfaceTension=config.liquidViscosity=config.gasViscosity=0;
    constexpr double dt=1./16384;const auto first=HydroVofStep::preparePrescribedFlow(*coarse,accepted,{1,0,dt,Clock::Position::End},config);auditStep(*coarse,accepted,first,model);first.commit(*coarse,accepted);
    HydroVofMeshTransferConfig remapConfig;remapConfig.momentumDensityRepresentation=Density::CellAverageProjection;
    for(int cycle=0;cycle<2;++cycle){const QuadtreeMesh& sourceMesh=cycle?*fine:*coarse;const QuadtreeMesh& targetMesh=cycle?*coarse:*fine;const auto old=accepted;
        const auto oldU=old.values().normalVelocity(),oldP=old.values().pressure();const auto& oldPub=old.phaseGeometry().vofPublication();const auto rawStocks=oldPub.plicCells();const auto oldIdentity=oldPub.identity();
        const auto targetDomain=domain(targetMesh,rz);const auto targetLayout=HydroDofLayout::build(targetMesh,targetDomain);const auto means=HydroVofMeshTransfer::boundaryNormalMeans(sourceMesh,old,targetMesh,targetDomain);
        auto targetLog=std::make_shared<std::vector<Query>>();auto targetFlow=flowForBackflow(targetMesh,targetLayout,means,model,targetLog);
        const auto volume=GeometricVofMeshTransfer::prepare(sourceMesh,oldPub,targetMesh,targetDomain,old.timeLayer());check(cycle?volume.report().coarsenedTargets>0:volume.report().refinedTargets>0,"F19 actual refine/coarsen native overlap");
        const auto sourceStart=sourceLog->size(),targetStart=targetLog->size();const auto remap=HydroVofMeshTransfer::preparePrescribedFlow(sourceMesh,old,targetMesh,volume,2,1,targetFlow,remapConfig);
        const auto& target=remap.candidateState();auditProjection(sourceMesh,old,targetMesh,remap);queryCoverage(sourceMesh,old,targetMesh,target,model,sourceLog,sourceStart,targetLog,targetStart,cycle!=0);
        auto changed=targetFlow;changed.pressureOutletBackflowModel=model==Model::Ambient?Model::InteriorContinuationZeroShear:Model::Ambient;
        changed.pressureOutletBackflowMomentum=model==Model::Ambient?HydroVofPrescribedFlow::PressureOutletBackflowMomentumQuery{}:HydroVofPrescribedFlow::PressureOutletBackflowMomentumQuery([](int,Clock){return ambient;});
        guard([&]{HydroVofMeshTransfer::preparePrescribedFlow(sourceMesh,old,targetMesh,volume,2,1,changed,remapConfig);},"target changes the actual pressure-outlet backflow model");
        changed=targetFlow;changed.pressureOutletBackflowModel.reset();changed.pressureOutletBackflowMomentum={};
        guard([&]{HydroVofMeshTransfer::preparePrescribedFlow(sourceMesh,old,targetMesh,volume,2,1,changed,remapConfig);},"target changes the actual pressure-outlet backflow model");
        if(model==Model::InteriorContinuationZeroShear){changed=targetFlow;changed.pressureOutletBackflowMomentum=[](int,Clock){return ambient;};
            guard([&]{HydroVofMeshTransfer::preparePrescribedFlow(sourceMesh,old,targetMesh,volume,2,1,changed,remapConfig);},"interior continuation cannot consume ambient remap data");}
        else{const auto originalQuery=targetFlow.pressureOutletBackflowMomentum;
            changed=targetFlow;changed.pressureOutletBackflowMomentum={};guard([&]{HydroVofMeshTransfer::preparePrescribedFlow(sourceMesh,old,targetMesh,volume,2,1,changed,remapConfig);},"ambient remap requires actual source/target queries");
            changed=targetFlow;changed.pressureOutletBackflowMomentum=[originalQuery](int row,Clock t){auto v=originalQuery(row,t);v.x+=.0625;return v;};
            guard([&]{HydroVofMeshTransfer::preparePrescribedFlow(sourceMesh,old,targetMesh,volume,2,1,changed,remapConfig);},"target changes/merges distinct actual P0 ambient momentum vectors");
            changed=targetFlow;changed.pressureOutletBackflowMomentum=[originalQuery](int row,Clock t){auto v=originalQuery(row,t);v.y=std::numeric_limits<double>::quiet_NaN();return v;};
            guard([&]{HydroVofMeshTransfer::preparePrescribedFlow(sourceMesh,old,targetMesh,volume,2,1,changed,remapConfig);},"nonfinite/nonrepeatable actual ambient momentum on remap");
            const auto calls=std::make_shared<unsigned>(0);changed=targetFlow;changed.pressureOutletBackflowMomentum=[originalQuery,calls](int row,Clock t){auto v=originalQuery(row,t);if((*calls)++%2)v.x+=.0625;return v;};
            guard([&]{HydroVofMeshTransfer::preparePrescribedFlow(sourceMesh,old,targetMesh,volume,2,1,changed,remapConfig);},"nonfinite/nonrepeatable actual ambient momentum on remap");
            if(cycle){int alteredRow=-1;for(const auto& port:targetFlow.pressureOutlets){const auto rows=coveredSourceRows(sourceMesh,old,targetMesh,targetLayout,port.velocityDof);if(rows.size()>1){alteredRow=rows.back();break;}}
                check(alteredRow>=0,"F19 coarsen guard lacks an actual multi-source pressure face");auto varied=*old.vofPrescribedFlow();const auto oldQuery=varied.pressureOutletBackflowMomentum;
                varied.pressureOutletBackflowMomentum=[oldQuery,alteredRow](int row,Clock t){auto v=oldQuery(row,t);if(row==alteredRow)v.x+=.0625;return v;};
                check(old.retainedMass()!=nullptr,"F19 actual source missing retained M");const auto discontinuous=HydroTwoPhaseState::createVofPrescribedFlowRemapped(sourceMesh,old.phaseGeometry(),old.values(),*old.retainedMass(),varied);
                guard([&]{HydroVofMeshTransfer::preparePrescribedFlow(sourceMesh,discontinuous,targetMesh,volume,2,1,targetFlow,remapConfig);},"target changes/merges distinct actual P0 ambient momentum vectors");++distinctCoarsens;}}
        const auto begin=target.timeLayer();const Clock end{begin.step+1,begin.time(),begin.time()+dt,Clock::Position::End};const auto targetU=target.values().normalVelocity(),targetP=target.values().pressure();
        const auto targetGeometry=MaxwellCoupledContract::sourceSnapshot(target.phaseGeometry().vofPublication());
        const auto next=HydroVofStep::preparePrescribedFlow(targetMesh,target,end,config);auditStep(targetMesh,target,next,model);
        check(!next.report().usedVelocityHistory&&next.report().startupCorrection&&next.candidateState().velocityHistory(),"F19 actual next Step failed to restart remapped history");
        check(old.sameIdentity(accepted)&&old.values().normalVelocity()==oldU&&old.values().pressure()==oldP&&oldPub.identity()==oldIdentity
            &&target.values().normalVelocity()==targetU&&target.values().pressure()==targetP
            &&MaxwellCoupledContract::sourceSnapshot(target.phaseGeometry().vofPublication())==targetGeometry,"F19 remap/guards/prepare mutated accepted source or target");
        for(std::size_t ci=0;ci<rawStocks.size();++ci)check(oldPub.plicCells()[ci].inventory==rawStocks[ci].inventory,"F19 remap changed raw source L/G stocks");
        accepted=target;next.commit(targetMesh,accepted);check(sameClock(accepted.timeLayer(),end)&&accepted.vofPrescribedFlow()->pressureOutletBackflowModel==model,"F19 actual changed target after Step lost model/clock");sourceLog=targetLog;
    }
}
void run() {
    for(bool rz:{false,true})for(Model model:{Model::Ambient,Model::InteriorContinuationZeroShear})scene(rz,model);
    check(sceneCount==4&&remapCount==8&&stepCount==12&&initialNegative>0&&initialZero>0&&targetNegative>0&&postStepNegative>0
        &&ambientSpans>0&&ambientIncoming>0&&continuationIncoming>0&&distinctCoarsens==2,"F19 missing actual model/refine/coarsen/next-Step/backflow witnesses");
    std::cout<<std::setprecision(17)<<"F19_BACKFLOW_REMAP_PASS,scenes="<<sceneCount<<",remaps="<<remapCount<<",accepted_steps="<<stepCount<<",guards="<<guardCount
        <<",initial_negative="<<initialNegative<<",initial_zero="<<initialZero<<",target_negative="<<targetNegative<<",post_step_negative="<<postStepNegative
        <<",ambient_overlap_spans="<<ambientSpans<<",ambient_incoming="<<ambientIncoming<<",current_U_continuation="<<continuationIncoming<<",distinct_coarsen_rejections="<<distinctCoarsens
        <<",projection_energy_peak="<<peakEnergy<<",div_peak="<<peakDiv<<",phase_peak="<<peakPhase<<",maximum_WB="<<maximumWB<<",maximum_WRHS="<<maximumWRHS
        <<",COMPLETE_NATIVE_FEED_PORTS_ONLY_NO_APP_LONGRUN_NOZZLE_CONTACT_PIN_CLAIM\n";
}
} // namespace F19BackflowRemapContract
// END F19 ACTUAL BACKFLOW MODEL REMAP CONSUMER


// BEGIN F24 ACTUAL VOF VISCOSITY GEOMETRY REUSE CONTRACT
namespace F24VofViscosityGeometryReuseContract {
using Timer=std::chrono::steady_clock;
unsigned scenes=0,pairs=0,scalars=0,matrixEntries=0,points=0,guards=0,phaseFresh=0,coefficientFresh=0,zeroPairs=0;
template<class Number> void exact(Number a,Number b,const char* why) {
    check(std::isfinite(a)&&std::isfinite(b)&&a==b&&(a!=0||std::signbit(a)==std::signbit(b)),why);++scalars;
}
void sameValues(const std::vector<double>& a,const std::vector<double>& b,const char* why) {
    check(a.size()==b.size(),why);for(std::size_t k=0;k<a.size();++k)exact(a[k],b[k],why);
}
void sameMatrix(const HydroEnergyViscosity::Matrix& a,const HydroEnergyViscosity::Matrix& b,const char* why) {
    check(a.rows()==b.rows()&&a.cols()==b.cols()&&a.nonZeros()==b.nonZeros(),why);
    for(int row=0;row<a.outerSize();++row) {
        HydroEnergyViscosity::Matrix::InnerIterator ia(a,row),ib(b,row);
        for(;ia&&ib;++ia,++ib){check(ia.row()==ib.row()&&ia.col()==ib.col(),why);exact(ia.value(),ib.value(),why);++matrixEntries;}
        check(!ia&&!ib,why);
    }
}
bool differentValues(const std::vector<double>& a,const std::vector<double>& b) {
    check(a.size()==b.size(),"F24 freshness vectors need the same actual quadrature");
    for(std::size_t k=0;k<a.size();++k)if(a[k]!=b[k]||(a[k]==0&&std::signbit(a[k])!=std::signbit(b[k])))return true;
    return false;
}
template<class F> void rejectsWith(F&& f,const char* token,const char* why) {
    bool caught=false;try{f();}catch(const std::exception& e){caught=std::string(e.what()).find(token)!=std::string::npos;}
    check(caught,why);++guards;
}
void sameStrain(const HydroStrainOperator& a,const HydroStrainOperator& b) {
    check(!a.isStressCompatible()&&!b.isStressCompatible()&&a.components()==b.components(),
        "F24 viscosity must retain the original raw strain convention");
    check(a.points().size()==b.points().size(),"F24 raw volume point count differs");
    for(std::size_t k=0;k<a.points().size();++k){const auto& x=a.points()[k];const auto& y=b.points()[k];
        check(x.cellDof==y.cellDof,"F24 raw volume point native owner differs");
        exact(x.point.x,y.point.x,"F24 raw volume point x differs");exact(x.point.y,y.point.y,"F24 raw volume point y differs");
        exact(x.physicalWeight,y.physicalWeight,"F24 raw volume quadrature physical weight differs");++points;}
    sameMatrix(a.matrix(),b.matrix(),"F24 exact raw strain matrix differs");
}
void compare(const QuadtreeMesh& mesh,const HydroPhaseGeometry& phase,const HydroPhaseViscosity& fresh,
    const HydroPhaseViscosity& reused,const HydroPhaseViscosity& previous,const HydroNativeValues& velocity) {
    const auto& a=fresh.viscosity();const auto& b=reused.viscosity();const auto& old=previous.viscosity();
    check(fresh.materialModel()==HydroPhaseViscosity::MaterialModel::VofCellAverage&&reused.materialModel()==fresh.materialModel()
        &&!fresh.hasQuadrature()&&!reused.hasQuadrature(),"F24 nonzero public VOF material model changed");
    check(fresh.geometry().vofIdentity()==phase.vofIdentity()&&reused.geometry().vofIdentity()==phase.vofIdentity()
        &&sameClock(fresh.geometry().timeLayer(),phase.timeLayer())&&sameClock(reused.geometry().timeLayer(),phase.timeLayer()),
        "F24 current immutable phase/full clock was frozen to the previous source");
    check(&a.strain().matrix()!=&old.strain().matrix()&&&b.strain().matrix()==&old.strain().matrix()
        &&&b.strain().points()==&old.strain().points(),"F24 fresh must rebuild while reused must retain the actual prepared raw owner");
    check(&a.matrix()!=&old.matrix()&&&b.matrix()!=&old.matrix()&&&a.matrix()!=&b.matrix(),
        "F24 current viscosity K must be freshly assembled in both paths");
    check(a.strain().layout().sameIdentity(phase.layout())&&b.strain().layout().sameIdentity(phase.layout())
        &&&a.strain().metrics().divergence().matrix()==&phase.metrics().divergence().matrix()
        &&&b.strain().metrics().divergence().matrix()==&phase.metrics().divergence().matrix(),
        "F24 current matrix lost the actual native layout/divergence owner");
    exact(fresh.liquidViscosity(),reused.liquidViscosity(),"F24 liquid coefficient differs");
    exact(fresh.gasViscosity(),reused.gasViscosity(),"F24 gas coefficient differs");
    sameValues(a.viscosity(),b.viscosity(),"F24 exact quadrature material vector differs");
    sameValues(a.rowWeights(),b.rowWeights(),"F24 exact row weights differ");sameStrain(a.strain(),b.strain());
    sameMatrix(a.matrix(),b.matrix(),"F24 exact current viscosity K differs");
    sameValues(a.strain().apply(mesh,velocity),b.strain().apply(mesh,velocity),"F24 exact strain application differs");
    sameValues(a.apply(mesh,velocity),b.apply(mesh,velocity),"F24 exact K application differs");
    const auto da=a.dissipation(mesh,velocity),db=b.dissipation(mesh,velocity);
    exact(da,db,"F24 exact viscous dissipation differs");check(da>0,"F24 nonzero viscosity oracle needs actual nonzero strain");++pairs;
}
void compareZero(const QuadtreeMesh& mesh,const HydroPhaseGeometry& phase,const HydroPhaseViscosity& previous,
    const HydroNativeValues& velocity) {
    const auto fresh=HydroPhaseViscosity::buildVofCellAverage(mesh,phase,0,0);
    const auto reused=HydroPhaseViscosity::buildVofCellAverage(mesh,phase,0,0,previous);
    for(const auto* value:{&fresh,&reused}) {
        const auto& visc=value->viscosity();
        check(value->materialModel()==HydroPhaseViscosity::MaterialModel::ExactZero&&!value->hasQuadrature()
            &&value->geometry().vofIdentity()==phase.vofIdentity()&&sameClock(value->geometry().timeLayer(),phase.timeLayer()),
            "F24 zero path must retain the actual current phase and full clock");
        check(visc.strain().points().empty()&&visc.strain().matrix().rows()==0&&visc.viscosity().empty()&&visc.rowWeights().empty()
            &&visc.matrix().rows()==phase.layout().velocityCount()&&visc.matrix().cols()==phase.layout().velocityCount()
            &&visc.matrix().nonZeros()==0,"F24 exact zero must have no sampled strain/material or hidden viscosity");
    }
    sameMatrix(fresh.viscosity().matrix(),reused.viscosity().matrix(),"F24 exact zero K differs");
    sameValues(fresh.viscosity().apply(mesh,velocity),reused.viscosity().apply(mesh,velocity),"F24 exact zero application differs");
    exact(fresh.viscosity().dissipation(mesh,velocity),0.,"F24 original zero dissipation changed");
    exact(reused.viscosity().dissipation(mesh,velocity),0.,"F24 reused zero dissipation changed");++zeroPairs;
}
void scene(bool rz) {
    ++scenes;auto mesh=meshFor(false);GeometricVofInitializationConfig init;init.subcellsPerAxis=8;
    auto pub=GeometricVofInitialization::prepare(*mesh,rz,[](Vector2D p){const double x=p.x-.47,y=p.y-.48;return x*x+y*y-.23*.23;},init);
    const auto first=phaseFor(*mesh,pub);const auto metrics=first.metrics();
    check(mesh->cellIndexNum()==40&&!first.layout().domain().uniformLeaves(),"F24 needs genuine 40-leaf hanging 2:1 XY/RZ scene");
    auto previous=HydroPhaseViscosity::buildVofCellAverage(*mesh,first,.0625,.015625);
    const auto initialPub=pub;const auto saved=MaxwellCoupledContract::sourceSnapshot(pub);
    const auto savedMu=previous.viscosity().viscosity();const auto sourceVelocity=valuesFor(*mesh,first);
    std::vector<double> normals(pub.domain().nativeFaceCount());
    for(int u=0;u<first.layout().velocityCount();++u)normals[first.layout().nativeFace(u)]=sourceVelocity.normalVelocity()[u];
    check(std::any_of(normals.begin(),normals.end(),[](double x){return x!=0;}),"F24 transport cannot be stationary for phase freshness");

    const auto exactZero=HydroPhaseViscosity::zero(*mesh,first);
    auto changedMesh=meshFor(true);const auto changedPhase=phaseFor(*changedMesh,initial(*changedMesh,rz));
    const auto changedPrevious=HydroPhaseViscosity::buildVofCellAverage(*changedMesh,changedPhase,.0625,.015625);
    check(changedMesh->cellIndexNum()!=mesh->cellIndexNum(),"F24 foreign mesh guard needs actual changed leaves");
    const auto foreignMetrics=HydroPhysicalDivergence::build(*mesh,HydroIntegralDivergence::build(*mesh,first.layout()));
    const auto foreignPhase=HydroPhaseGeometry::build(*mesh,foreignMetrics,pub);
    const auto foreignPrevious=HydroPhaseViscosity::buildVofCellAverage(*mesh,foreignPhase,.0625,.015625);
    check(foreignPhase.layout().sameIdentity(first.layout())&&
        &foreignMetrics.divergence().matrix()!=&metrics.divergence().matrix(),"F24 foreign D must share actual layout but have another owner");
    sameMatrix(foreignMetrics.divergence().matrix(),metrics.divergence().matrix(),"F24 foreign D must be numerically identical");
    const auto foreignLayout=HydroDofLayout::build(*mesh,pub.domain());
    const auto foreignLayoutMetrics=HydroPhysicalDivergence::build(*mesh,HydroIntegralDivergence::build(*mesh,foreignLayout));
    const auto foreignLayoutPhase=HydroPhaseGeometry::build(*mesh,foreignLayoutMetrics,pub);
    const auto foreignLayoutPrevious=HydroPhaseViscosity::buildVofCellAverage(*mesh,foreignLayoutPhase,.0625,.015625);
    check(!foreignLayout.sameIdentity(first.layout()),"F24 foreign layout needs genuinely distinct owner");
    const auto fan=GeometricVofPublication::prepare(*mesh,pub.domain(),pub.timeLayer(),std::vector<double>(mesh->cellIndexNum(),0));
    const auto fanPhase=HydroPhaseGeometry::build(*mesh,metrics,fan);
    check(fan.geometryKind()==GeometricVofGeometryKind::CompatibleFan,"F24 wrong current phase needs actual incompatible geometric kind");
    rejectsWith([&]{HydroPhaseViscosity::buildVofCellAverage(*mesh,fanPhase,.125,.03125);},
        "requires actual LocalPlic inventories","F24 original wrong current geometric source must fail");
    rejectsWith([&]{HydroPhaseViscosity::buildVofCellAverage(*mesh,fanPhase,.125,.03125,previous);},
        "requires actual LocalPlic inventories","F24 reuse must retain wrong current geometric source qualification");
    long double freshSeconds=0,reuseSeconds=0;auto lastMid=pub.identity();
    for(int stage=0;stage<2;++stage) {
        const auto old=pub.timeLayer();const Clock end{old.step+1,old.time(),old.time()+1e-5,Clock::Position::End};
        const auto slab=GeometricVofAdvectorSlab::prepare(*mesh,pub,end,normals,normals);
        const auto volume=GeometricVofTransport::prepareLocalPlicSweep(*mesh,pub,slab,10,1);
        const auto phase=HydroPhaseGeometry::build(*mesh,metrics,volume.midpoint());const auto velocity=valuesFor(*mesh,phase);
        check(phase.vofIdentity()==volume.midpoint().identity()&&phase.vofIdentity()!=lastMid
            &&phase.timeLayer().position==Clock::Position::Midpoint&&sameClock(phase.timeLayer(),volume.midpoint().timeLayer()),
            "F24 each current phase must be a distinct actual transported Midpoint");
        lastMid=phase.vofIdentity();
        bool stocksChanged=false;const auto& before=previous.geometry().vofPublication().plicCells();const auto& after=phase.vofPublication().plicCells();
        for(std::size_t k=0;k<after.size();++k)stocksChanged=stocksChanged||before[k].inventory!=after[k].inventory;
        check(stocksChanged,"F24 phase freshness needs actual changed independent transported stocks");
        const auto sameMuFresh=HydroPhaseViscosity::buildVofCellAverage(*mesh,phase,previous.liquidViscosity(),previous.gasViscosity());
        const auto sameMuReused=HydroPhaseViscosity::buildVofCellAverage(*mesh,phase,previous.liquidViscosity(),previous.gasViscosity(),previous);
        compare(*mesh,phase,sameMuFresh,sameMuReused,previous,velocity);
        check(differentValues(sameMuReused.viscosity().viscosity(),previous.viscosity().viscosity()),
            "F24 same mu must read changed current stocks rather than freeze prior material");++phaseFresh;
        const double liquid=stage?.09375:.125,gas=stage?.0078125:.03125;
        const auto fresh=HydroPhaseViscosity::buildVofCellAverage(*mesh,phase,liquid,gas);
        const auto reused=HydroPhaseViscosity::buildVofCellAverage(*mesh,phase,liquid,gas,previous);
        compare(*mesh,phase,fresh,reused,previous,velocity);
        exact(reused.liquidViscosity(),liquid,"F24 current liquid coefficient frozen");exact(reused.gasViscosity(),gas,"F24 current gas coefficient frozen");
        check(differentValues(reused.viscosity().viscosity(),sameMuReused.viscosity().viscosity()),
            "F24 same current phase must read new viscosity coefficients");++coefficientFresh;
        compareZero(*mesh,phase,previous,velocity);
        rejectsWith([&]{HydroPhaseViscosity::buildVofCellAverage(*mesh,phase,liquid,gas,exactZero);},
            "previous VofCellAverage source","F24 ExactZero cannot authenticate raw volume quadrature");
        rejectsWith([&]{HydroPhaseViscosity::buildVofCellAverage(*mesh,phase,0,0,exactZero);},
            "previous VofCellAverage source","F24 zero request must retain prior owner qualification");
        rejectsWith([&]{HydroPhaseViscosity::buildVofCellAverage(*mesh,phase,liquid,gas,changedPrevious);},
            "","F24 actual changed-mesh prior must fail");
        rejectsWith([&]{HydroPhaseViscosity::buildVofCellAverage(*mesh,phase,liquid,gas,foreignLayoutPrevious);},
            "different actual native layout","F24 numerically equivalent rebuilt layout must fail");
        rejectsWith([&]{HydroPhaseViscosity::buildVofCellAverage(*mesh,phase,liquid,gas,foreignPrevious);},
            "actual phase divergence owner","F24 numerically identical same-layout foreign D must fail");
        for(double bad:{-1.,std::numeric_limits<double>::infinity(),std::numeric_limits<double>::quiet_NaN()}) {
            rejectsWith([&]{HydroPhaseViscosity::buildVofCellAverage(*mesh,phase,bad,gas);},
                "negative/nonfinite phase viscosity","F24 original current coefficient guard lost");
            rejectsWith([&]{HydroPhaseViscosity::buildVofCellAverage(*mesh,phase,bad,gas,previous);},
                "negative/nonfinite phase viscosity","F24 reuse current coefficient guard lost");
        }
        if(stage==0)for(int sample=0;sample<3;++sample) {
            auto begin=Timer::now();const auto a=HydroPhaseViscosity::buildVofCellAverage(*mesh,phase,liquid,gas);
            freshSeconds+=std::chrono::duration<long double>(Timer::now()-begin).count();
            begin=Timer::now();const auto b=HydroPhaseViscosity::buildVofCellAverage(*mesh,phase,liquid,gas,previous);
            reuseSeconds+=std::chrono::duration<long double>(Timer::now()-begin).count();
            compare(*mesh,phase,a,b,previous,velocity);
        }
        pub=volume.endpoint();previous=reused;
    }
    check(MaxwellCoupledContract::sourceSnapshot(initialPub)==saved&&sameClock(initialPub.timeLayer(),first.timeLayer()),
        "F24 current material preparation mutated retained Initial source");
    sameValues(savedMu,HydroPhaseViscosity::buildVofCellAverage(*mesh,first,.0625,.015625).viscosity().viscosity(),
        "F24 original source material changed after transported current sources");
    std::cout<<std::setprecision(17)<<"F24_VOF_VISCOSITY_GEOMETRY_REUSE_SCENE,"<<(rz?"RZ":"XY")<<",cells,"<<mesh->cellIndexNum()
        <<",actual_transported_midpoints,2,timing_samples,3,fresh_prepare_seconds,"<<freshSeconds<<",reused_prepare_seconds,"<<reuseSeconds
        <<",timing_scope,small_public_prepare_only_no_gate_no_whole_App_claim\n";
}
void run() {
    scene(false);scene(true);
    check(scenes==2&&pairs==14&&phaseFresh==4&&coefficientFresh==4&&zeroPairs==4&&guards==48,
        "F24 missing actual scene/source/coefficient/zero/rejection coverage");
    std::cout<<"F24_VOF_VISCOSITY_GEOMETRY_REUSE_PASS,scenes,"<<scenes<<",pairs,"<<pairs<<",scalar_values,"<<scalars
        <<",stored_matrix_entries,"<<matrixEntries<<",raw_points,"<<points<<",guards,"<<guards<<",phase_fresh,"<<phaseFresh
        <<",coefficient_fresh,"<<coefficientFresh<<",zero_pairs,"<<zeroPairs
        <<",comparison,exact_value_and_signed_zero,actual_transported_stocks,no_Hermite_prior_fixture,no_App_commit\n";
}
} // namespace F24VofViscosityGeometryReuseContract
// END F24 ACTUAL VOF VISCOSITY GEOMETRY REUSE CONTRACT


// Task35: real VOF remap and real target State consumption, with a separate cold preparation.
namespace Task35VofRegularPreparationRemap {
using Regular=LAMP::Hydrodynamics::Detail::HydroRegularAdvectorGeometry;
std::size_t scenes=0,comparedRows=0,comparedPorts=0,comparedCoefficients=0,comparedJets=0,guards=0;
long double peak=0;
std::unique_ptr<QuadtreeMesh> meshFor(bool rz,bool refined)
{
    DictionaryMesh d;d.Nx=d.Ny=8;d.maxCellNum=16384;d.maxInitialIterations=0;
    const double height=rz?1:2;d.minCoord={0,0};d.maxCoord={1,height};auto m=std::make_unique<QuadtreeMesh>(d);
    DictionaryMesh::RefineBox center;center.level=1;center.min={.125,.125*height};center.max={.875,.875*height};
    std::vector<DictionaryMesh::RefineBox> boxes{center};if(refined){auto corner=center;corner.min={0,0};corner.max={.125,.125*height};boxes.push_back(corner);}
    m->applyRefineBoxes(boxes);return m;
}
HydroNativeValues sourceValues(const QuadtreeMesh& mesh,const HydroPhaseGeometry& phase)
{
    const auto potential=[&](Vector2D p){const long double x=(p.x-mesh.minCoord.x)/(mesh.maxCoord.x-mesh.minCoord.x),y=(p.y-mesh.minCoord.y)/(mesh.maxCoord.y-mesh.minCoord.y);
        return 16*x*x*(1-x)*(1-x)*y*y*(1-y)*(1-y)*(1+.3L*x+.7L*y);};
    const auto& l=phase.layout();const bool rz=l.domain().axisymmetric();std::vector<double> u(l.velocityCount()),p(l.pressureCount());
    for(int i=0;i<l.velocityCount();++i){const auto f=l.domain().face(mesh,l.nativeFace(i));const auto& e=f.horizontal?mesh.hEdges(f.edgeLabel):mesh.vEdges(f.edgeLabel);
        auto a=mesh.verts(e.edgeToVert[0]).coord(),b=mesh.verts(e.edgeToVert[1]).coord();if(f.horizontal?a.x>b.x:a.y>b.y)std::swap(a,b);
        const auto flux=(f.horizontal?-1.L:1.L)*(rz?CartesianMeasure::twoPi:1.L)*(potential(b)-potential(a));u[i]=f.area>0?static_cast<double>(flux/f.area):0;}
    for(int i=0;i<l.pressureCount();++i)p[i]=mesh.cells(mesh.cellIndexToLabel[l.nativeCell(i)]).centerCoord().y-(mesh.minCoord.y+mesh.maxCoord.y)/2;
    return {l,u,p};
}
void number(long double a,long double b,const char* why)
{near(a,b,1+std::abs(b),why);peak=std::max(peak,std::abs(a-b)/(1+std::abs(b)));++comparedCoefficients;}
void compare(const QuadtreeMesh& mesh,const HydroTwoPhaseState& warm)
{
    // This is a cold target reconstruction, not a second independent remap solve.
    // The original verify() below independently checks the actual projection and stocks.
    const auto coldGeometry=HydroMomentumGeometry::build(mesh,warm.phaseGeometry().metrics());
    const auto coldMass=HydroMomentumMass::buildCellAverageProjection(mesh,coldGeometry,warm.phaseGeometry(),10,1);
    const auto cold=HydroTwoPhaseState::createVofRemapped(mesh,warm.phaseGeometry(),warm.values(),coldMass);
    check(!cold.sameIdentity(warm)&&!cold.momentumGeometry().sameIdentity(warm.momentumGeometry()),"Task35 State cold preparation aliases reused target");
    check(!cold.momentumGeometry().regularPreparationReuseReport().requested,"Task35 independent cold State claimed a remap request");
    check(cold.phaseGeometry().vofIdentity()==warm.phaseGeometry().vofIdentity()&&sameClock(cold.timeLayer(),warm.timeLayer())&&cold.values().layout().sameIdentity(warm.values().layout())&&!cold.velocityHistory(),"Task35 cold State changed physical owner/fullclock/layout");
    const auto& a=warm.reconstruction().regularGeometry();const auto& b=cold.reconstruction().regularGeometry();
    check(!a.sameIdentity(b)&&!a.usesCompatibleTensor()&&!a.usesVirtualCompatibleTensor()&&!b.usesCompatibleTensor()&&!b.usesVirtualCompatibleTensor(),"Task35 State positive must use independent general preparations");
    check(a.cells().size()==b.cells().size()&&a.faces().size()==b.faces().size()&&a.rows().size()==b.rows().size()&&a.normalPorts().size()==b.normalPorts().size(),"Task35 cold State topology inventory");
    for(std::size_t p=0;p<a.cells().size();++p){const auto& x=a.cells()[p];const auto& y=b.cells()[p];
        check(x.nativeCell==y.nativeCell&&x.lower.x==y.lower.x&&x.lower.y==y.lower.y&&x.upper.x==y.upper.x&&x.upper.y==y.upper.y&&x.vertices==y.vertices&&x.edges==y.edges&&x.faces==y.faces&&x.supportKnots==y.supportKnots,"Task35 State fresh target cell topology");}
    for(std::size_t i=0;i<a.rows().size();++i){const auto& x=a.rows()[i];const auto& y=b.rows()[i];
        check(x.kind==y.kind&&x.index==y.index&&x.component==y.component&&x.rank==y.rank&&x.rings==y.rings&&x.degree==y.degree&&x.p2Rank==y.p2Rank&&x.p2ColumnCount==y.p2ColumnCount&&x.zeroNormal==y.zeroNormal&&x.columns==y.columns&&x.support==y.support,"Task35 State full ordered row routing/rank");
        for(int k=0;k<4;++k){check(x.weights[k].size()==y.weights[k].size(),"Task35 State row coefficient count");for(std::size_t j=0;j<x.weights[k].size();++j)number(x.weights[k][j],y.weights[k][j],"Task35 State cold row coefficient");}++comparedRows;}
    for(std::size_t i=0;i<a.normalPorts().size();++i){const auto& x=a.normalPorts()[i];const auto& y=b.normalPorts()[i];
        check(x.vertex==y.vertex&&x.component==y.component&&x.observedModes==y.observedModes&&x.support==y.support,"Task35 State ordered normal port routing");
        for(int k=0;k<2;++k){check(x.weights[k].size()==y.weights[k].size(),"Task35 State port coefficient count");for(std::size_t j=0;j<x.weights[k].size();++j)number(x.weights[k][j],y.weights[k][j],"Task35 State cold normal port coefficient");}++comparedPorts;}
    const auto av=a.apply(mesh,warm.values(),&warm.boundary()),bv=b.apply(mesh,cold.values(),&cold.boundary());
    check(av.cellDofs.size()==bv.cellDofs.size()&&av.faceCorrections.size()==bv.faceCorrections.size(),"Task35 State cold applied inventory");
    for(std::size_t p=0;p<av.cellDofs.size();++p)for(int c=0;c<2;++c)for(int k=0;k<25;++k)number(av.cellDofs[p][c][k],bv.cellDofs[p][c][k],"Task35 State cold applied cell coefficient");
    for(std::size_t u=0;u<av.faceCorrections.size();++u)number(av.faceCorrections[u],bv.faceCorrections[u],"Task35 State cold face correction");
    for(std::size_t p=0;p<a.cells().size();++p){const auto& c=a.cells()[p];for(Vector2D f:std::array<Vector2D,2>{{{.31,.42},{.67,.73}}}){
        const Vector2D x=c.lower+Vector2D{f.x*(c.upper.x-c.lower.x),f.y*(c.upper.y-c.lower.y)};const auto u=a.cellJet(av,static_cast<int>(p),x).jet,v=b.cellJet(bv,static_cast<int>(p),x).jet;
        const std::array<long double,6> actual{u.value.x,u.value.y,u.derivative.xx,u.derivative.xy,u.derivative.yx,u.derivative.yy},ref{v.value.x,v.value.y,v.derivative.xx,v.derivative.xy,v.derivative.yx,v.derivative.yy};
        for(int k=0;k<6;++k)number(actual[k],ref[k],"Task35 State cold actual jet");++comparedJets;}}
    const auto r=warm.momentumGeometry().regularPreparationReuseReport();
    check(r.requested&&r.enabled&&r.reusedRows>0&&r.reusedPorts>0&&r.recomputedRows>0&&r.reusedRows+r.recomputedRows==a.rows().size()&&r.reusedPorts+r.recomputedPorts==a.normalPorts().size()&&r.targetDatums>0&&r.retainedMetadataBytes>0,"Task35 actual State did not consume qualified overlap preparation");
    check(cold.values().normalVelocity()==warm.values().normalVelocity()&&cold.values().pressure()==warm.values().pressure(),"Task35 cold reconstruction altered actual remap U/p");
}
void scenario(bool rz,bool coarsen)
{
    auto targetMesh=meshFor(rz,!coarsen);
    const auto target=[&]{auto sourceMesh=meshFor(rz,coarsen);const auto pub=initial(*sourceMesh,rz);const auto phase0=phaseFor(*sourceMesh,pub);const auto values=sourceValues(*sourceMesh,phase0);
        const auto first=HydroTwoPhaseState::createVof(*sourceMesh,phase0,values);
        check(!first.reconstruction().regularGeometry().usesCompatibleTensor()&&!first.reconstruction().regularGeometry().usesVirtualCompatibleTensor(),"Task35 State source positive must use general reconstruction");
        std::vector<std::array<long double,2>> inventory(pub.plicCells().size());for(std::size_t i=0;i<inventory.size();++i)inventory[i]=pub.plicCells()[i].inventory;
        const auto end=GeometricVofPublication::prepareLocalPlicInventory(*sourceMesh,pub.domain(),{1,0,.125,Clock::Position::End},inventory,{});
        const auto phase=HydroPhaseGeometry::build(*sourceMesh,phase0.metrics(),end);
        const auto mass=HydroMomentumMass::buildCellAverageProjection(*sourceMesh,first.momentumGeometry(),phase,10,1);
        const auto old=HydroTwoPhaseState::createVofAccepted(*sourceMesh,phase,values,first,mass);
        const auto oldU=old.values().normalVelocity(),oldP=old.values().pressure(),oldAlpha=end.alpha();const auto oldStocks=stocks(end);
        const auto interfaceTransfer=GeometricVofMeshTransfer::prepare(*sourceMesh,end,*targetMesh,domain(*targetMesh,rz),old.timeLayer());
        check(coarsen?interfaceTransfer.report().coarsenedTargets>0:interfaceTransfer.report().refinedTargets>0,"Task35 State missing real refine/coarsen targets");
        HydroVofMeshTransferConfig cfg;cfg.momentumDensityRepresentation=Density::CellAverageProjection;
        const auto actual=HydroVofMeshTransfer::prepare(*sourceMesh,old,*targetMesh,interfaceTransfer,10,1,cfg);
        verify(*sourceMesh,old,*targetMesh,actual,Density::CellAverageProjection);compare(*targetMesh,actual.candidateState());
        check(old.values().normalVelocity()==oldU&&old.values().pressure()==oldP&&end.alpha()==oldAlpha&&stocks(end)==oldStocks&&old.phaseGeometry().vofIdentity()==end.identity()&&old.velocityHistory(),"Task35 preparation changed accepted source U/p/stocks/history");
        if(!rz&&!coarsen){const auto foreign=GeometricVofPublication::prepareLocalPlicInventory(*sourceMesh,pub.domain(),old.timeLayer(),inventory,{});
            const auto wrong=GeometricVofMeshTransfer::prepare(*sourceMesh,foreign,*targetMesh,domain(*targetMesh,rz),old.timeLayer());
            rejects([&]{HydroVofMeshTransfer::prepare(*sourceMesh,old,*targetMesh,wrong,10,1,cfg);},"Task35 foreign VOF owner accepted");
            rejects([&]{HydroVofMeshTransfer::prepare(*sourceMesh,first,*targetMesh,interfaceTransfer,10,1,cfg);},"Task35 Initial accepted as End remap source");guards+=2;}
        return actual.candidateState();
    }();
    // All source State, mesh, transfer, and source preparation owners above are gone.
    target.validateFor(*targetMesh);const auto samples=target.sampleCellCenterVelocities(*targetMesh);
    check(samples.size()==static_cast<std::size_t>(target.values().layout().pressureCount()),"Task35 target retained a required source lifetime");
    for(auto v:samples)check(std::isfinite(v.x)&&std::isfinite(v.y),"Task35 target source-lifetime sample nonfinite");
    const auto r=target.momentumGeometry().regularPreparationReuseReport();++scenes;
    std::cout<<"VOF_REGULAR_PREPARATION_REMAP,"<<rz<<','<<coarsen<<','<<targetMesh->cellIndexNum()<<','<<r.targetDatums<<','<<r.reusedRows<<','<<r.recomputedRows<<','<<r.reusedPorts<<','<<r.recomputedPorts<<','<<r.retainedMetadataBytes<<",PASS\n";
}
void run()
{
    for(bool rz:{false,true})for(bool coarsen:{false,true})scenario(rz,coarsen);
    check(scenes==4&&comparedRows>0&&comparedPorts>0&&comparedCoefficients>0&&comparedJets>0&&guards==2,"Task35 State cold/remap/guard coverage incomplete");
    std::cout<<"VOF_REGULAR_PREPARATION_REMAP_PASS,scenes,"<<scenes<<",rows,"<<comparedRows<<",ports,"<<comparedPorts<<",scalar_comparisons,"<<comparedCoefficients<<",jets,"<<comparedJets<<",guards,"<<guards<<",peak,"<<peak<<",real_remap_original_conservation_projection_oracle,cold_preparation_not_second_projection,source_owners_destroyed\n";
}
} // namespace Task35VofRegularPreparationRemap


// Task37: original matrix arithmetic versus the actual Frozen factory, on real VOF owners.
namespace Task37FrozenMatrixAssemblyPlan {
using Matrix=HydroFrozenMidpoint::Matrix;
std::size_t scenes=0,preparations=0,storedEntries=0,scalars=0,solves=0,fpPairs=0,guards=0,kChanges=0,hChanges=0;
std::size_t physicalUses=0,unionUses=0;
std::size_t compiledMissUses=0,equalNnzChanges=0,hOnlyInsertions=0,hOnlyRemovals=0;
struct Floating {
    std::fenv_t native{};unsigned csr=0;int error=errno;
    Floating(){csr=_mm_getcsr();check(std::fegetenv(&native)==0,"Task37 cannot capture actual FP state");}
    void restore()const{check(std::fesetenv(&native)==0,"Task37 cannot restore actual FP state");_mm_setcsr(csr);errno=error;}
    bool same(const Floating& b)const{return std::memcmp(&native,&b.native,sizeof native)==0&&csr==b.csr&&error==b.error;}
};
struct Inputs {
    HydroDualMaterial dual;HydroKepConvection convection;HydroKepConvectionResult applied;
    HydroPhaseViscosity viscosity;HydroPressureLoad pressure;HydroNativeValues before;
    HydroFrozenMidpointForces forces;HydroFrozenMidpointBoundary boundary;
};
void exact(double a,double b,const char* why)
{check(std::isfinite(a)&&std::isfinite(b)&&std::memcmp(&a,&b,sizeof a)==0,why);++scalars;}
template<class Sparse>void sameMatrix(const Sparse& a,const Sparse& b,const char* why)
{
    check(a.rows()==b.rows()&&a.cols()==b.cols()&&a.isCompressed()&&b.isCompressed()&&a.nonZeros()==b.nonZeros(),why);
    for(Eigen::Index i=0;i<=a.outerSize();++i)check(a.outerIndexPtr()[i]==b.outerIndexPtr()[i],why);
    for(Eigen::Index i=0;i<a.nonZeros();++i){check(a.innerIndexPtr()[i]==b.innerIndexPtr()[i],why);exact(a.valuePtr()[i],b.valuePtr()[i],why);++storedEntries;}
}
bool samePattern(const Matrix& a,const Matrix& b)
{
    if(a.rows()!=b.rows()||a.cols()!=b.cols()||a.nonZeros()!=b.nonZeros()||!a.isCompressed()||!b.isCompressed())return false;
    for(Eigen::Index i=0;i<=a.outerSize();++i)if(a.outerIndexPtr()[i]!=b.outerIndexPtr()[i])return false;
    for(Eigen::Index i=0;i<a.nonZeros();++i)if(a.innerIndexPtr()[i]!=b.innerIndexPtr()[i])return false;return true;
}
bool different(const Matrix& a,const Matrix& b)
{
    if(!samePattern(a,b))return true;
    for(Eigen::Index i=0;i<a.nonZeros();++i)if(std::memcmp(a.valuePtr()+i,b.valuePtr()+i,sizeof(double))!=0)return true;return false;
}
Inputs inputs(const QuadtreeMesh& mesh,const HydroMomentumGeometry& geometry,const GeometricVofPublication& pub,
    double dt,double rho,double mu,double speed,int curlVertex=0)
{
    const auto& metrics=geometry.metrics();const auto a=HydroPhaseGeometry::build(mesh,metrics,pub);const auto raw=valuesFor(mesh,a);
    auto velocity=raw.normalVelocity();for(double& v:velocity)v*=speed;
    const auto& l=geometry.layout();
    if(curlVertex){
        check(curlVertex>=1&&curlVertex<=15,"Task42 local streamfunction vertex fixture out of range");
        const Vector2D node{.25+((curlVertex-1)%5)/8.,.625+((curlVertex-1)/5)/8.};
        const auto localPsi=[&](Vector2D p)->long double{return p.x==node.x&&p.y==node.y?1.L/256:0.L;};
        bool actualVertex=false;for(int i=0;i<mesh.vertIndexNum();++i){const auto p=mesh.verts(mesh.vertIndexToLabel[i]).coord();actualVertex=actualVertex||(p.x==node.x&&p.y==node.y);}
        check(actualVertex,"Task42 local curl must use a real interior native vertex");
        for(int u=0;u<l.velocityCount();++u){const auto f=l.domain().face(mesh,l.nativeFace(u));
            const auto& e=f.horizontal?mesh.hEdges(f.edgeLabel):mesh.vEdges(f.edgeLabel);
            auto x=mesh.verts(e.edgeToVert[0]).coord(),y=mesh.verts(e.edgeToVert[1]).coord();
            if(f.horizontal?x.x>y.x:x.y>y.y)std::swap(x,y);
            const long double flux=(f.horizontal?-1.L:1.L)*(l.domain().axisymmetric()?CartesianMeasure::twoPi:1.L)*(localPsi(y)-localPsi(x));
            velocity[u]=f.area>0?static_cast<double>(flux/f.area)*speed:0.;}
    }
    const HydroNativeValues before(l,velocity,std::vector<double>(l.pressureCount()));
    std::vector<double> normals(pub.domain().nativeFaceCount());for(int u=0;u<l.velocityCount();++u)normals[l.nativeFace(u)]=velocity[u];
    const Clock end{pub.timeLayer().step+1,pub.timeLayer().time(),pub.timeLayer().time()+dt,Clock::Position::End};
    const auto slab=GeometricVofAdvectorSlab::prepare(mesh,pub,end,normals,normals);
    const auto moved=GeometricVofTransport::prepareLocalPlicSweep(mesh,pub,slab,10*rho,rho);
    const auto middle=HydroPhaseGeometry::build(mesh,metrics,moved.midpoint()),b=HydroPhaseGeometry::build(mesh,metrics,moved.endpoint());
    std::vector<int> fixed;for(int u=0;u<l.velocityCount();++u){const auto face=l.domain().face(mesh,l.nativeFace(u));if(face.negativeCell<0||face.positiveCell<0)fixed.push_back(u);}
    const auto walls=HydroAdvectorBoundary::impermeable(mesh,geometry,fixed);
    const auto material=HydroMaterialTransport::prepareGeometricVof(mesh,geometry,a,middle,b,moved,walls);
    const auto ma=HydroMomentumMass::buildCellAverageProjection(mesh,geometry,a,10*rho,rho),
        mm=HydroMomentumMass::buildCellAverageProjection(mesh,geometry,middle,10*rho,rho),mb=HydroMomentumMass::buildCellAverageProjection(mesh,geometry,b,10*rho,rho);
    const auto dual=HydroDualMaterial::prepare(mesh,material,ma,mm,mb);
    const auto convection=HydroKepConvection::buildImpermeable(mesh,dual,walls);
    const auto boundary=HydroConvectionBoundary::bind(mesh,dual,std::vector<double>(dual.faceRates().size()),dual.sources(),dual.timeLayer());
    const auto applied=convection.apply(mesh,before,boundary,dual.sources(),dual.timeLayer());
    const auto viscosity=HydroPhaseViscosity::buildVofCellAverage(mesh,middle,.0625*mu,.015625*mu);
    HydroFrozenMidpointForces force;force.bodyAcceleration={0,.125};
    return {dual,convection,applied,viscosity,HydroPressureLoad::build(metrics.divergence()),before,force,HydroFrozenMidpointBoundary::closed(walls)};
}
using Coordinate=std::pair<Eigen::Index,Eigen::Index>;
std::vector<Coordinate> hOnlyCoordinates(const Matrix& h)
{
    // These fixtures use genuine exact-zero K. Positive actual M/dt makes P
    // the complete diagonal, so off-diagonal H coordinates are exactly H-only.
    std::vector<Coordinate> result;for(Eigen::Index row=0;row<h.outerSize();++row)
        for(Matrix::InnerIterator e(h,row);e;++e)if(e.col()!=row)result.emplace_back(row,e.col());return result;
}
std::size_t coordinateDifference(const std::vector<Coordinate>& a,const std::vector<Coordinate>& b)
{std::size_t count=0;for(const auto& q:a)count+=!std::binary_search(b.begin(),b.end(),q);return count;}
std::pair<Inputs,Inputs> shiftedEqualCountH(const QuadtreeMesh& mesh,const HydroMomentumGeometry& geometry,const GeometricVofPublication& pub)
{
    std::vector<Inputs> flows;std::vector<Matrix> matrices;
    // Bounded deterministic fixture selection uses ONLY real public transport/
    // material/KEP factories. No Frozen prepare or solve is run in this search,
    // no private H is injected, and no cache report chooses an expected result.
    for(int vertex=1;vertex<=15;++vertex){
        auto candidate=inputs(mesh,geometry,pub,1./32768,2,0,.25,vertex);
        check(candidate.viscosity.viscosity().matrix().nonZeros()==0,"Task42 local H fixture lost original exact-zero K");
        const Matrix actual=candidate.convection.frozenMomentumOperator(mesh,candidate.dual.sources(),candidate.dual.timeLayer());
        if(actual.nonZeros()>0){const auto next=hOnlyCoordinates(actual);
            for(std::size_t j=0;j<matrices.size();++j)if(actual.nonZeros()==matrices[j].nonZeros()&&!samePattern(actual,matrices[j])){
                const auto previous=hOnlyCoordinates(matrices[j]);const auto added=coordinateDifference(next,previous),removed=coordinateDifference(previous,next);
                if(added&&removed){hOnlyInsertions+=added;hOnlyRemovals+=removed;++equalNnzChanges;
                    return {flows[j],std::move(candidate)};}}
            flows.push_back(std::move(candidate));matrices.push_back(actual);}
    }
    throw std::runtime_error("Task42 bounded real curl fixtures lack equal-nnz/different-coordinate H with both H-only insertion/removal");
}
HydroFrozenMidpoint prepare(const QuadtreeMesh& mesh,const Inputs& in,HydroFrozenMidpointScheme scheme=HydroFrozenMidpointScheme::FrozenConvection)
{return HydroFrozenMidpoint::prepare(mesh,in.dual,in.convection,in.applied,in.viscosity,in.pressure,in.before,in.forces,in.boundary,scheme,in.dual.sources(),in.dual.timeLayer());}
void report(const HydroSaddleReport& r)
{
    check(r.algebraicRelativeResidual<=1e-11&&r.freeMomentumDefect<=1e-11&&r.maxCellContinuityDefect<=1e-11&&r.maxPressureMeanDefect<=1e-11&&r.maxBoundaryFluxDefect<=1e-12,"Task37 original accepted Saddle gates changed");
}
void oracle(const QuadtreeMesh& mesh,const Inputs& in,const HydroFrozenMidpoint& actual)
{
    actual.validateFor(mesh,in.dual.sources(),in.dual.timeLayer());check(actual.solveCount()==1,"Task37 real Frozen factory did not solve exactly once");
    check(&actual.viscosityMatrix()==&actual.viscosity().viscosity().matrix(),"Task37 K getter is not the retained actual Energy matrix");
    const auto& l=in.dual.layout();const auto& k=in.viscosity.viscosity().matrix();const double dt=in.dual.deltaTime();const bool linear=actual.scheme()==HydroFrozenMidpointScheme::FrozenConvection;
    // Independent original literal sparse expressions, never the compiled routes.
    Matrix physical=.5*k;for(int u=0;u<l.velocityCount();++u)physical.coeffRef(u,u)+=in.dual.endMass().masses()[u]/dt;physical.makeCompressed();
    Matrix h(l.velocityCount(),l.velocityCount());if(linear)h=in.convection.frozenMomentumOperator(mesh,in.dual.sources(),in.dual.timeLayer());
    Matrix a=physical;if(linear)a+=.5*h;a.makeCompressed();
    sameMatrix(actual.viscosityMatrix(),k,"Task37 complete original K CSR bits");sameMatrix(actual.physicalVelocityMatrix(),physical,"Task37 complete original halfK+mass/dt CSR bits");
    sameMatrix(actual.frozenConvectionMatrix(),h,"Task37 independently produced actual H CSR bits");sameMatrix(actual.velocityMatrix(),a,"Task37 complete original physical+halfH CSR bits");
    const Eigen::Map<const Eigen::VectorXd> u0(in.before.normalVelocity().data(),l.velocityCount());
    const Eigen::VectorXd ku=k*u0;Eigen::VectorXd hu=Eigen::VectorXd::Zero(l.velocityCount());if(linear)hu=h*u0;
    for(int u=0;u<l.velocityCount();++u){exact(actual.beginViscousAction()[u],ku[u],"Task37 original K*u0 action bits");exact(actual.beginFrozenAction()[u],hu[u],"Task37 original H*u0 action bits");}
    const auto exterior=linear?in.convection.frozenBoundaryMomentumLoad(mesh,in.dual.sources(),in.dual.timeLayer()):std::vector<long double>(l.velocityCount());
    std::vector<double> rhs(l.velocityCount());for(int u=0;u<l.velocityCount();++u){long double value=static_cast<long double>(in.dual.beginMass().masses()[u])*u0[u]/dt-.5L*ku[u];
        if(linear){value-=.5L*hu[u];value+=exterior[u];}else value+=in.applied.momentumLoads()[u];
        const double acceleration=l.normalAxis(u)==HydroDofLayout::NormalAxis::First?in.forces.bodyAcceleration.x:in.forces.bodyAcceleration.y;
        const long double body=static_cast<long double>(in.dual.midpointMass().masses()[u])*acceleration;
        if(in.forces.bodyAcceleration.x!=0||in.forces.bodyAcceleration.y!=0)value+=body;rhs[u]=static_cast<double>(value);
        exact(actual.integralMomentumRhs()[u],rhs[u],"Task37 independent original full momentum RHS bits");exact(static_cast<double>(actual.bodyLoads()[u]),static_cast<double>(body),"Task37 actual midpoint body mass is fresh");}
    const auto velocity=linear?HydroSaddleVelocityBlock::buildWithTransportSplit(l,a):HydroSaddleVelocityBlock::build(l,a);
    const auto saddle=HydroSinglePhaseSaddle::buildClosed(in.pressure,velocity);
    sameMatrix(actual.saddle().matrix(),saddle.matrix(),"Task37 independent original monolithic system CSR bits");
    const auto solved=saddle.solve(mesh,HydroSaddleRhs(l,rhs,in.boundary.normalVelocity));report(solved.report());report(actual.solveResult().report());
    const auto& x=actual.solveResult().values();const auto& y=solved.values();
    for(int u=0;u<l.velocityCount();++u)near(x.normalVelocity()[u],y.normalVelocity()[u],1+std::abs(y.normalVelocity()[u]),"Task37 original literal accepted velocity solution");
    for(int p=0;p<l.pressureCount();++p)near(x.pressure()[p],y.pressure()[p],1+std::abs(y.pressure()[p]),"Task37 original literal accepted pressure solution");
    const auto divergence=in.dual.geometry().metrics().divergence().apply(mesh,x);const auto& di=in.dual.geometry().metrics().divergence().matrix();
    for(int p=0;p<l.pressureCount();++p){long double scale=0;for(HydroIntegralDivergence::Matrix::InnerIterator e(di,p);e;++e)scale+=std::abs(static_cast<long double>(e.value())*x.normalVelocity()[e.col()]);
        check(scale>0?std::abs(divergence[p])/scale<=1e-11L:divergence[p]==0,"Task37 independent original native continuity gate");}
    const auto& r=actual.matrixAssemblyPlanReport();check(r.viscosityNonZeros==static_cast<std::size_t>(k.nonZeros())&&r.convectionNonZeros==static_cast<std::size_t>(h.nonZeros())&&r.physicalNonZeros==static_cast<std::size_t>(physical.nonZeros())&&r.velocityNonZeros==static_cast<std::size_t>(a.nonZeros()),"Task37 actual per-preparation nnz report");
    check(r.retainedBytes<=128u*1024u*1024u,"Task37 pure index stock exceeds actual cap");
    if(r.physicalValuesUsed)check(r.requested&&r.eligible&&r.physicalPlanHit,"Task37 physical use lacks current qualification/binding");
    if(r.unionValuesUsed)check(r.requested&&r.eligible&&r.physicalPlanHit&&r.physicalValuesUsed,"Task42 union use lacks qualified actual fixed K/P owner");
    physicalUses+=r.physicalValuesUsed;unionUses+=r.unionValuesUsed;++preparations;++solves;
}
void paired(const QuadtreeMesh& mesh,const Inputs& in,const char* label,bool physicalHit,bool unionHit)
{
    // Both complete factories start from the SAME genuinely captured native/raw
    // entry, not fabricated _Fe_stat encodings. Original arithmetic remains live.
    check(std::feclearexcept(FE_ALL_EXCEPT)==0,"Task37 cannot clear actual entry flags");_mm_setcsr(_mm_getcsr()&~0x3fu);errno=173;const Floating entry;
    const auto first=prepare(mesh,in);const Floating firstExit;entry.restore();const auto warm=prepare(mesh,in);const Floating warmExit;
    check(firstExit.same(warmExit),"Task37 whole cold/warm factory FP controls/status/errno differ");++fpPairs;
    const auto& a=first.matrixAssemblyPlanReport();const auto& b=warm.matrixAssemblyPlanReport();
    check(a.requested&&a.eligible&&a.physicalPlanHit==physicalHit&&a.physicalValuesUsed==physicalHit&&a.unionPlanHit==unionHit&&a.unionValuesUsed==physicalHit,"Task42 complete actual H-route hit and newly compiled route use differ from real fixture");
    if(!a.unionPlanHit&&a.unionValuesUsed)++compiledMissUses;
    check(b.requested&&b.eligible&&b.physicalPlanHit&&b.physicalValuesUsed&&b.unionPlanHit&&b.unionValuesUsed&&b.retainedBytes>0,"Task37 repeated actual preparation did not use both fresh value stages");
    oracle(mesh,in,first);oracle(mesh,in,warm);
    std::cout<<"FROZEN_MATRIX_PLAN,"<<label<<','<<a.physicalPlanHit<<','<<a.unionPlanHit<<','<<a.viscosityNonZeros<<','<<a.convectionNonZeros<<','<<a.physicalNonZeros<<','<<a.velocityNonZeros<<','<<b.retainedBytes<<",PASS\n";
}
template<class Action>void rejectWith(Action action,const char* token)
{bool caught=false;try{action();}catch(const std::exception& e){caught=std::string(e.what()).find(token)!=std::string::npos;}check(caught,"Task37 original rejection missing or wrong reason");++guards;}
void scenario(bool rz)
{
    auto mesh=meshFor(false);const auto pub=initial(*mesh,rz);const auto phase=phaseFor(*mesh,pub);const auto geometry=HydroMomentumGeometry::build(*mesh,phase.metrics());
    check(mesh->cellIndexNum()==40&&!geometry.layout().domain().uniformLeaves(),"Task37 requires actual 40-leaf XY/RZ 2:1");
    const auto raw=pub.plicCells();const auto alpha=pub.alpha();
    auto base=inputs(*mesh,geometry,pub,1./16384,1,1,.25);paired(*mesh,base,rz?"RZ_cold_warm":"XY_cold_warm",false,false);
    auto freshK=base;freshK.viscosity=HydroPhaseViscosity::buildVofCellAverage(*mesh,base.dual.material().midpointGeometry(),.125,.03125);
    check(samePattern(base.viscosity.viscosity().matrix(),freshK.viscosity.viscosity().matrix())&&different(base.viscosity.viscosity().matrix(),freshK.viscosity.viscosity().matrix()),"Task37 fresh K must retain pattern and change real values");
    paired(*mesh,freshK,"fresh_K_same_pattern",true,true);
    auto freshMass=inputs(*mesh,geometry,pub,1./16384,2,2,.25);
    check(freshMass.dual.endMass().masses()!=base.dual.endMass().masses(),"Task37 fresh actual mass did not change");
    const auto oldH=base.convection.frozenMomentumOperator(*mesh,base.dual.sources(),base.dual.timeLayer()),newH=freshMass.convection.frozenMomentumOperator(*mesh,freshMass.dual.sources(),freshMass.dual.timeLayer());
    check(different(oldH,newH),"Task37 fresh actual H values did not change");
    paired(*mesh,freshMass,"fresh_M_and_H",samePattern(freshK.viscosity.viscosity().matrix(),freshMass.viscosity.viscosity().matrix()),samePattern(freshK.viscosity.viscosity().matrix(),freshMass.viscosity.viscosity().matrix())&&samePattern(oldH,newH));
    auto freshDt=inputs(*mesh,geometry,pub,1./32768,2,2,.25);check(freshDt.dual.deltaTime()!=freshMass.dual.deltaTime(),"Task37 genuine full-clock dt did not change");
    const auto dtH=freshDt.convection.frozenMomentumOperator(*mesh,freshDt.dual.sources(),freshDt.dual.timeLayer());
    const bool dtKsame=samePattern(freshMass.viscosity.viscosity().matrix(),freshDt.viscosity.viscosity().matrix());
    paired(*mesh,freshDt,"fresh_dt",dtKsame,dtKsame&&samePattern(newH,dtH));
    auto zero=inputs(*mesh,geometry,pub,1./32768,2,0,0);check(zero.viscosity.viscosity().matrix().nonZeros()==0,"Task37 zero viscosity is not original exact zero");
    const auto zeroH=zero.convection.frozenMomentumOperator(*mesh,zero.dual.sources(),zero.dual.timeLayer());check(zeroH.nonZeros()==0,"Task37 stationary physical H must be actual empty pattern");
    paired(*mesh,zero,"changed_K_to_exact_zero",false,false);++kChanges;
    const auto shifted=shiftedEqualCountH(*mesh,geometry,pub);const auto& movingZero=shifted.first;const auto& movedZero=shifted.second;
    const auto movingH=movingZero.convection.frozenMomentumOperator(*mesh,movingZero.dual.sources(),movingZero.dual.timeLayer());
    const auto movedH=movedZero.convection.frozenMomentumOperator(*mesh,movedZero.dual.sources(),movedZero.dual.timeLayer());
    check(movingH.nonZeros()>0&&!samePattern(zeroH,movingH),"Task37 H-change fixture requires actual changed complete pattern");
    check(movingH.nonZeros()==movedH.nonZeros()&&!samePattern(movingH,movedH),"Task42 same stored H count must still change complete coordinates");
    const auto oldOnly=hOnlyCoordinates(movingH),newOnly=hOnlyCoordinates(movedH);
    check(!oldOnly.empty()&&coordinateDifference(newOnly,oldOnly)>0&&coordinateDifference(oldOnly,newOnly)>0,"Task42 true H-only coordinate insertion/removal missing");
    paired(*mesh,movingZero,"changed_H_only",true,false);++hChanges;
    paired(*mesh,movedZero,"equal_nnz_moved_H_only_columns",true,false);++hChanges;
    paired(*mesh,base,"return_nonzero_K_pattern",false,false);++kChanges;
    const auto otherGeometry=HydroMomentumGeometry::build(*mesh,phase.metrics());check(!otherGeometry.sameIdentity(geometry)&&otherGeometry.layout().sameIdentity(geometry.layout())&& &otherGeometry.metrics().divergence().matrix()== &geometry.metrics().divergence().matrix(),"Task37 independent owner must be genuinely new with same actual layout/D");
    const auto other=inputs(*mesh,otherGeometry,pub,1./16384,1,1,.25);paired(*mesh,other,"fresh_actual_Momentum_owner",false,false);
    const auto lagged=prepare(*mesh,base,HydroFrozenMidpointScheme::LaggedConvection);const auto& lr=lagged.matrixAssemblyPlanReport();
    check(!lr.requested&&!lr.eligible&&!lr.physicalPlanHit&&!lr.physicalValuesUsed&&!lr.unionPlanHit&&!lr.unionValuesUsed&&lr.retainedBytes==0,"Task37 Lagged path claimed optional plan participation");oracle(*mesh,base,lagged);
    auto wrongClock=base.dual.timeLayer();wrongClock.end=std::nextafter(wrongClock.end,1.);
    rejectWith([&]{lagged.validateFor(*mesh,base.dual.sources(),wrongClock);},"full clock mismatch");
    auto foreignSources=base.dual.sources();foreignSources.vofIds[1]=pub.identity();rejectWith([&]{lagged.validateFor(*mesh,foreignSources,base.dual.timeLayer());},"foreign VOF owner");
    auto foreignPressure=base;const auto alien=HydroPhysicalDivergence::build(*mesh,HydroIntegralDivergence::build(*mesh,geometry.layout()));foreignPressure.pressure=HydroPressureLoad::build(alien.divergence());
    rejectWith([&]{prepare(*mesh,foreignPressure);},"different authoritative viscosity/pressure D_I");
    auto foreignValues=base;const auto alienLayout=HydroDofLayout::build(*mesh,pub.domain());foreignValues.before=HydroNativeValues(alienLayout,base.before.normalVelocity(),base.before.pressure());
    rejectWith([&]{prepare(*mesh,foreignValues);},"layout identity");
    auto foreignWall=base;foreignWall.boundary=other.boundary;rejectWith([&]{prepare(*mesh,foreignWall);},"different impermeable momentum geometry instance");
    const Floating saved;check(std::fesetround(FE_UPWARD)==0,"Task37 cannot install unsupported rounding");
    try{rejectWith([&]{base.boundary.walls.validateFor(*mesh,geometry);},"round-to-nearest");rejectWith([&]{prepare(*mesh,base);},"round-to-nearest");}catch(...){saved.restore();throw;}saved.restore();
    const unsigned csr=_mm_getcsr();_mm_setcsr(csr|(1U<<15));
    try{rejectWith([&]{base.boundary.walls.validateFor(*mesh,geometry);},"gradual underflow");rejectWith([&]{prepare(*mesh,base);},"gradual underflow");}catch(...){_mm_setcsr(csr);throw;}_mm_setcsr(csr);
    check(pub.alpha()==alpha&&pub.plicCells().size()==raw.size(),"Task37 source primary inventory storage changed");
    for(std::size_t i=0;i<raw.size();++i)check(pub.plicCells()[i].inventory==raw[i].inventory,"Task37 original independent source stocks changed");++scenes;
}
void run()
{
    const Floating caller;try{scenario(false);scenario(true);}catch(...){caller.restore();throw;}caller.restore();
    check(scenes==2&&preparations==38&&solves==38&&fpPairs==18&&guards==18&&kChanges==4&&hChanges==4&&physicalUses>0&&unionUses>0,"Task37 missing actual cold/warm/fresh/pattern/owner/original guard coverage");
    check(equalNnzChanges==2&&hOnlyInsertions>0&&hOnlyRemovals>0&&compiledMissUses>=4,"Task42 missing actual equal-count H-only insertion/removal and miss-used coverage");
    std::cout<<"FROZEN_DYNAMIC_H_ROUTE_PASS,equal_nnz_different_pattern_pairs,"<<equalNnzChanges<<",H_only_inserted_coordinates,"<<hOnlyInsertions<<",H_only_removed_coordinates,"<<hOnlyRemovals<<",new_complete_H_route_miss_used,"<<compiledMissUses<<",literal_CSR_RHS_solve_FP_oracle_preserved,1\n";
    std::cout<<"FROZEN_MATRIX_ASSEMBLY_PLAN_PASS,scenes,"<<scenes<<",preparations,"<<preparations<<",independent_reference_solves,"<<solves<<",whole_factory_FP_errno_pairs,"<<fpPairs<<",guards,"<<guards<<",K_pattern_changes,"<<kChanges<<",H_pattern_changes,"<<hChanges<<",stored_CSR_entries,"<<storedEntries<<",scalar_bits,"<<scalars<<",physical_used,"<<physicalUses<<",union_used,"<<unionUses<<",original_equations_and_local_gates,1,stock_cap,134217728\n";
}
} // namespace Task37FrozenMatrixAssemblyPlan


namespace SelectedDomainHydroAmrContract {
using Mode=GeometricVofAdvectorSlab::Boundary::Mode;
unsigned scenes=0,adaptations=0,physicalSteps=0,solidChecks=0,wallChecks=0,firstSolidIncidents=0,guards=0;
long double peakR=0,peakEnergy=0,peakProjection=0;
std::unique_ptr<QuadtreeMesh> nozzle(bool fine) {
    // Four actual base cells across x<1 retain the original P2 strain rank.
    DictionaryMesh d;d.Nx=8;d.Ny=24;d.minCoord={0,0};d.maxCoord={2,6};
    d.maxInitialIterations=0;d.maxCellNum=4096;d.maxLevel_local={1};
    auto mesh=std::make_unique<QuadtreeMesh>(d);DictionaryMesh::RefineBox top;
    top.level=1;top.min={0,5};top.max={2,6};std::vector<DictionaryMesh::RefineBox> boxes{top};
    if(fine){auto neck=top;neck.min={0,3};neck.max={1,4};boxes.push_back(neck);}
    mesh->applyRefineBoxes(boxes);return mesh;
}
NativeLeafDomain selected(const QuadtreeMesh& mesh,bool rz,bool exchanged=false) {
    FaceSegmentGraph graph;graph.build(mesh,[exchanged](const QuadtreeMesh&,const Cell& cell){
        const auto p=cell.centerCoord();return exchanged?!(p.x>1&&p.y>3):!(p.x>1&&p.y<3);});
    return NativeLeafDomain::build(mesh,graph,rz);
}
std::array<Vector2D,2> ends(const QuadtreeMesh& mesh,const NativeLeafFace& f) {
    const auto& e=f.horizontal?mesh.hEdges(f.edgeLabel):mesh.vEdges(f.edgeLabel);
    auto a=mesh.verts(e.edgeToVert[0]).coord(),b=mesh.verts(e.edgeToVert[1]).coord();
    if(f.horizontal?a.x>b.x:a.y>b.y)std::swap(a,b);return {a,b};
}
HydroVofPrescribedFlow pressureFlow(const QuadtreeMesh& mesh,const HydroDofLayout& l) {
    HydroVofPrescribedFlow flow;flow.pressureDrivenReservoirs=true;
    flow.pressureOutletBackflowModel=HydroVofPressureOutletMomentumModel::InteriorContinuationZeroShear;
    flow.normalVelocity.assign(l.velocityCount(),0);flow.rulesByNativeRow.assign(l.domain().nativeFaceCount(),Mode::Closed);
    for(int u=0;u<l.velocityCount();++u){const auto f=l.domain().face(mesh,l.nativeFace(u));
        if(f.outwardSign==0)flow.rulesByNativeRow[f.nativeRow]=Mode::Interior;
        else if(f.horizontal&&(f.center.y==0||f.center.y==6)){
            const bool bottom=f.center.y==0;flow.rulesByNativeRow[f.nativeRow]=bottom?Mode::PureLiquidReservoir:Mode::PureGasReservoir;
            flow.pressureOutlets.push_back({u,bottom?.125:0});}}
    return flow;
}
void inspect(const QuadtreeMesh& mesh,const HydroTwoPhaseState& state) {
    const auto& l=state.values().layout();const auto& pub=state.phaseGeometry().vofPublication();const auto& d=pub.domain();
    check(l.pressureCount()==static_cast<int>(d.cellIndexByRow().size())&&l.pressureCount()<mesh.cellIndexNum()
        &&l.velocityCount()<d.nativeFaceCount(),"selected AMR added solid pressure/velocity DOFs");
    const auto indicators=HydroVofRefinement::sample(mesh,state);
    check(indicators.vorticity.size()==static_cast<std::size_t>(mesh.cellIndexNum())&&indicators.hAbsVorticity.size()==indicators.vorticity.size(),
        "selected passive AMR indicator is not dense native order");
    double maximum=0;for(int ci=0;ci<mesh.cellIndexNum();++ci){const int p=l.pressureDof(ci);
        if(d.cellRowByIndex()[ci]<0){check(p<0&&pub.alpha()[ci]==0&&pub.cells()[ci].capacity==0
            &&pub.plicCells()[ci].inventory==std::array<long double,2>{0,0}&&indicators.vorticity[ci]==0&&indicators.hAbsVorticity[ci]==0,
            "selected solid gained DOF, inventory, capacity or refinement indicator");++solidChecks;
        }else{check(p==d.cellRowByIndex()[ci]&&l.nativeCell(p)==ci,"selected pressure map lost exact inverse");
            maximum=std::max(maximum,indicators.hAbsVorticity[ci]);}}
    check(maximum==indicators.maximumHAbsVorticity,"selected indicator maximum includes inactive cells");
    const CartesianLeafLocator locator(mesh);for(const auto f:d.faces(mesh)){
        const int u=l.velocityDof(f.nativeRow);if(f.negativeCell<0&&f.positiveCell<0){check(u<0,"solid native face gained a velocity DOF");++solidChecks;continue;}
        check(u>=0&&l.nativeFace(u)==f.nativeRow,"selected velocity map lost exact inverse");
        const bool electrode=(f.horizontal&&f.center.y==3&&f.center.x>1)||(!f.horizontal&&f.center.x==1&&f.center.y<3);
        if(electrode){check(state.values().normalVelocity()[u]==0&&std::find(state.boundary().velocityDofs().begin(),state.boundary().velocityDofs().end(),u)
            !=state.boundary().velocityDofs().end(),"selected internal solid wall lost exact impermeability");++wallChecks;
            const auto incident=locator.incident(f.center);bool active=false,solid=false;
            for(int k=0;k<incident.count;++k){active=active||l.pressureDof(incident.indices[k])>=0;solid=solid||l.pressureDof(incident.indices[k])<0;}
            if(active&&solid&&l.pressureDof(incident.indices[0])<0)++firstSolidIncidents;
        }
    }
}
void step(const QuadtreeMesh& mesh,HydroTwoPhaseState& accepted,bool open,const HydroVofAcConfig& config) {
    const auto before=accepted;const auto oldU=before.values().normalVelocity(),oldP=before.values().pressure();
    const auto t=before.timeLayer();const Clock end{t.step+1,t.time(),t.time()+1./16384,Clock::Position::End};
    std::cerr<<"SELECTED_HYDRO_AMR_STEP_BEGIN rz="<<before.values().layout().domain().axisymmetric()<<" pressure_only="<<open
        <<" accepted_step="<<t.step<<" next_step="<<end.step<<" native_cells="<<mesh.cellIndexNum()
        <<" active_cells="<<before.values().layout().pressureCount()<<'\n'<<std::flush;
    const auto actual=open?HydroVofStep::preparePrescribedFlow(mesh,accepted,end,config):HydroVofStep::prepare(mesh,accepted,end,config);
    const auto& dual=actual.frozenMidpoint().material();const auto& report=dual.report();const auto& transport=actual.transport();
    constexpr long double massGate=256*std::numeric_limits<double>::epsilon();
    check(dual.usesConservativeVofMomentum()&&report.maxRelativeContinuityDefect<=massGate&&report.maxRelativeStoredMassContinuityDefect<=massGate
        &&report.maxRelativePrimalContinuityDefect<=massGate&&report.maxRelativeMassTransferDefect<=massGate&&report.maxRelativeFluxCommutatorDefect<=massGate,
        "selected post-AMR original actual M/F/R gates changed");
    peakR=std::max(peakR,report.maxRelativeStoredMassContinuityDefect);
    for(const auto& r:actual.report().stages){check(r.consistencyPower==0&&r.momentumPowerRelativeDefect<=1e-11&&r.energy.closureRelativeResidual<=1e-11L
        &&r.saddle.algebraicRelativeResidual<=1e-11&&r.saddle.freeMomentumDefect<=1e-11&&r.saddle.maxCellContinuityDefect<=1e-11,
        "selected post-AMR original Saddle/momentum/energy gates changed");peakEnergy=std::max(peakEnergy,r.energy.closureRelativeResidual);}
    check(accepted.sameIdentity(before)&&accepted.values().normalVelocity()==oldU&&accepted.values().pressure()==oldP
        &&transport.sourceIdentity()==before.phaseGeometry().vofIdentity(),"selected AMR physical prepare mutated accepted source");
    for(const auto f:before.phaseGeometry().vofPublication().domain().faces(mesh)){
        const bool solid=f.negativeCell<0&&f.positiveCell<0;
        const auto u=before.values().layout().velocityDof(f.nativeRow);
        const bool wall=u>=0&&std::find(before.boundary().velocityDofs().begin(),before.boundary().velocityDofs().end(),u)!=before.boundary().velocityDofs().end();
        if(solid||wall){const auto& q=transport.faces()[f.nativeRow].full;
            check(q.total==0&&q.liquid==0&&q.gas==0&&q.positiveTotal==0&&q.negativeTotal==0,
                "selected solid/wall acquired physical shared flux");}}
    if(open){const auto& ledger=transport.report().boundaryFull;
        check(ledger.enabled&&ledger.physicalRelativeConservation<=1e-12L&&ledger.massRelativeConservation<=1e-12L,
            "selected pressure supply original phase/mass boundary gates changed");
        check(actual.candidateState().vofPrescribedFlow()->pressureDrivenReservoirs&&actual.frozenMidpoint().saddle().pressureGaugeCount()==0,
            "selected pressure supply lost free pressure reference after AMR");
        long double in=0,out=0;const auto& l=accepted.values().layout();
        for(const auto& port:accepted.vofPrescribedFlow()->pressureOutlets){const auto f=l.domain().face(mesh,l.nativeFace(port.velocityDof));
            const auto q=f.outwardSign*f.area*actual.candidateState().values().normalVelocity()[port.velocityDof];
            in+=std::max(0.,-q);out+=std::max(0.,q);
            check(actual.frozenMidpoint().saddle().pressureOutletLoads()[port.velocityDof]==-f.outwardSign*f.area*port.pressure,
                "selected remap changed actual physical pressure traction");}
        check(in>0&&out>0,"selected pressure-only supply lacks actual solved incoming/outgoing flow");
    }
    actual.commit(mesh,accepted);check(sameClock(accepted.timeLayer(),end),"selected physical commit lost actual full clock");++physicalSteps;inspect(mesh,accepted);
}
template<class Action>void reject(Action&& action,const char* why) {rejects(std::forward<Action>(action),why);++guards;}
// A current selected RZ 2:1 pressure-remap witness. Ambient zero callbacks
// force the untouched generic authentication loop for the SAME instantaneous
// source u/p/alpha/mass and pressure/phase/wall data. This does not assert the
// two backflow models have the same later physical time-step law.
void pairedBoundarySourceRoles(const QuadtreeMesh& sourceMesh,const HydroTwoPhaseState& source,
    const QuadtreeMesh& targetMesh,const GeometricVofMeshTransfer& volumes,const HydroVofPrescribedFlow& targetFlow,
    const HydroVofMeshTransferConfig& config,const HydroVofMeshTransfer& observed)
{
    using Floating=Task37FrozenMatrixAssemblyPlan::Floating;
    using Model=HydroVofPressureOutletMomentumModel;
    const Floating caller;
    try{
        check(source.vofPrescribedFlow()&&source.vofPrescribedFlow()->pressureDrivenReservoirs
            &&targetFlow.pressureDrivenReservoirs&&source.retainedMass(),"F08 needs the actual retained pressure-remap source");
        auto literalFlow=*source.vofPrescribedFlow();literalFlow.pressureOutletBackflowModel=Model::Ambient;
        unsigned sourceCallbacks=0,targetCallbacks=0;
        literalFlow.pressureOutletBackflowMomentum=[&](int,Clock){++sourceCallbacks;return Vector2D{};};
        const auto literalSource=HydroTwoPhaseState::createVofPrescribedFlowRemapped(sourceMesh,
            source.phaseGeometry(),source.values(),*source.retainedMass(),literalFlow);
        auto literalTarget=targetFlow;literalTarget.pressureOutletBackflowModel=Model::Ambient;
        literalTarget.pressureOutletBackflowMomentum=[&](int,Clock){++targetCallbacks;return Vector2D{};};
        check(std::feraiseexcept(FE_INEXACT)==0,"F08 cannot prime common FP entry");errno=173;const Floating entry;
        const auto sourceCallsBefore=sourceCallbacks,targetCallsBefore=targetCallbacks;
        const auto literal=HydroVofMeshTransfer::preparePrescribedFlow(sourceMesh,literalSource,targetMesh,volumes,2,1,literalTarget,config);
        const Floating literalExit;entry.restore();
        const auto filtered=HydroVofMeshTransfer::preparePrescribedFlow(sourceMesh,source,targetMesh,volumes,2,1,targetFlow,config);
        const Floating filteredExit;check(literalExit.same(filteredExit),"F08 full remap FP/CSR/errno differs from generic literal witness");
        check(sourceCallbacks>sourceCallsBefore&&targetCallbacks>targetCallsBefore,"F08 generic literal witness did not consume actual ambient callbacks");
        const auto equal=[](const auto& a,const auto& b,const char* why){check(a.size()==b.size(),why);
            for(std::size_t i=0;i<a.size();++i)check(a[i]==b[i]&&(a[i]!=0||std::signbit(a[i])==std::signbit(b[i])),why);};
        const auto compare=[&](const HydroVofMeshTransfer& a,const HydroVofMeshTransfer& b){
            equal(a.rawValues().normalVelocity(),b.rawValues().normalVelocity(),"F08 raw original polynomial means differ");
            equal(a.rawValues().pressure(),b.rawValues().pressure(),"F08 raw physical pressure differs");
            equal(a.candidateState().values().normalVelocity(),b.candidateState().values().normalVelocity(),"F08 projected original velocity differs");
            equal(a.candidateState().values().pressure(),b.candidateState().values().pressure(),"F08 candidate physical pressure differs");
            equal(a.projectionMomentumRhs(),b.projectionMomentumRhs(),"F08 original projection RHS differs");
            equal(a.projectionMultiplier(),b.projectionMultiplier(),"F08 original projection multiplier differs");
            equal(stocks(a.candidateState().phaseGeometry().vofPublication()),stocks(b.candidateState().phaseGeometry().vofPublication()),"F08 independent target phase stocks differ");
            const auto& ac=a.candidateState().phaseGeometry().vofPublication().plicCells();
            const auto& bc=b.candidateState().phaseGeometry().vofPublication().plicCells();
            check(ac.size()==bc.size(),"F08 target inventory table differs");
            for(std::size_t i=0;i<ac.size();++i)equal(ac[i].inventory,bc[i].inventory,"F08 independent cell phase stocks differ");
            equal(a.candidateState().phaseGeometry().vofPublication().alpha(),b.candidateState().phaseGeometry().vofPublication().alpha(),"F08 target alpha differs");
            const auto& x=a.report();const auto& y=b.report();
            check(x.boundaryCoverageSpans==y.boundaryCoverageSpans&&x.fixedNormals==y.fixedNormals
                &&x.pressureOutletNormals==y.pressureOutletNormals&&x.wallNormals==y.wallNormals
                &&x.axisNormals==y.axisNormals&&x.faceSpans==y.faceSpans&&x.polynomialIntegrals==y.polynomialIntegrals
                &&x.sourceLeafCandidates==y.sourceLeafCandidates,"F08 original coverage/measure work receipts differ");
            const std::array<long double,8> ax{{x.projectionBoundaryWork,x.projectionRhsRoundingWork,x.projectionFreeResidualWork,
                x.projectionContinuityWork,x.projectionEnergyIdentityResidual,x.projectionEnergyIdentityScale,
                x.projectionEnergyIdentityRelativeResidual,x.maximumFinalRelativeDivergence}};
            const std::array<long double,8> by{{y.projectionBoundaryWork,y.projectionRhsRoundingWork,y.projectionFreeResidualWork,
                y.projectionContinuityWork,y.projectionEnergyIdentityResidual,y.projectionEnergyIdentityScale,
                y.projectionEnergyIdentityRelativeResidual,y.maximumFinalRelativeDivergence}};
            equal(ax,by,"F08 original same-End projection/energy receipt differs");
            check(sameClock(a.candidateState().timeLayer(),b.candidateState().timeLayer())
                &&a.candidateState().phaseGeometry().vofIdentity()==b.candidateState().phaseGeometry().vofIdentity(),
                "F08 target alpha owner/fullclock differs");
        };
        compare(literal,filtered);compare(observed,filtered);
        auto badLiteral=literalTarget,badFiltered=targetFlow;
        badLiteral.pressureOutlets.front().pressure+=.0625;badFiltered.pressureOutlets.front().pressure+=.0625;
        const auto rejected=[&](const HydroTwoPhaseState& old,const HydroVofPrescribedFlow& flow){
            try{(void)HydroVofMeshTransfer::preparePrescribedFlow(sourceMesh,old,targetMesh,volumes,2,1,flow,config);}
            catch(const std::exception& error){return std::string(error.what());}return std::string{};
        };
        entry.restore();const auto literalError=rejected(literalSource,badLiteral);const Floating literalFailure;
        entry.restore();const auto filteredError=rejected(source,badFiltered);const Floating filteredFailure;
        check(!literalError.empty()&&literalError==filteredError
            &&literalError.find("target merges/changes physical outlet pressure")!=std::string::npos,
            "F08 original first pressure rejection differs");
        check(literalFailure.same(filteredFailure),"F08 first rejection FP/CSR/errno differs");
        check(source.sameIdentity(observed.sourceState())&&sameClock(source.timeLayer(),observed.sourceState().timeLayer()),
            "F08 paired witnesses mutated the actual source owner/fullclock");
        std::cout<<"F08_BOUNDARY_SOURCE_ROLES_PAIR_PASS,selected_RZ_2to1,1,generic_literal_callbacks,"
            <<sourceCallbacks+targetCallbacks<<",full_remap_FP_errno,1,original_values_and_receipts,1,first_error_FP,1\n";
    }catch(...){caller.restore();throw;}
    caller.restore();
}
void scene(bool rz,bool open) {
    std::cerr<<"SELECTED_HYDRO_AMR_SCENE_BEGIN rz="<<rz<<" pressure_only="<<open<<'\n'<<std::flush;
    ++scenes;auto coarse=nozzle(false),fine=nozzle(true);const auto d=selected(*coarse,rz);
    const auto pub=GeometricVofInitialization::prepare(*coarse,d,[](Vector2D p){return std::max(p.x-1,p.y-4.125);});
    const auto phase=phaseFor(*coarse,pub);const auto& l=phase.layout();std::vector<double> u(l.velocityCount());
    if(!open){const auto psi=[](Vector2D p){if(p.x>=1||p.y<=3||p.y>=5)return 0.L;
        return static_cast<long double>(p.x)*p.x*(1-p.x)*(1-p.x)*(p.y-3)*(p.y-3)*(5-p.y)*(5-p.y)/256;};
        for(int k=0;k<l.velocityCount();++k){const auto f=d.face(*coarse,l.nativeFace(k));const auto e=ends(*coarse,f);
            const auto flux=(f.horizontal?-1.L:1.L)*(rz?CartesianMeasure::twoPi:1.L)*(psi(e[1])-psi(e[0]));
            u[k]=f.area>0?static_cast<double>(flux/f.area):0;}}
    const HydroNativeValues values(l,u,std::vector<double>(l.pressureCount()));
    auto accepted=open?HydroTwoPhaseState::createVofPrescribedFlow(*coarse,phase,values,pressureFlow(*coarse,l)):
        HydroTwoPhaseState::createVof(*coarse,phase,values);
    HydroVofAcConfig config;config.liquidDensity=open?2:10;config.gasDensity=1;config.liquidViscosity=.01;config.gasViscosity=.005;config.surfaceTension=0;
    if(open)check(HydroVofStep::stabilityRatesPrescribedFlow(*coarse,accepted,config).flowRate==0,"selected pressure AMR initial state invented a velocity seed");
    step(*coarse,accepted,open,config);HydroVofMeshTransferConfig remapConfig;remapConfig.momentumDensityRepresentation=Density::CellAverageProjection;
    for(int cycle=0;cycle<4;++cycle){const auto& sourceMesh=cycle%2?*fine:*coarse;const auto& targetMesh=cycle%2?*coarse:*fine;
        std::cerr<<"SELECTED_HYDRO_AMR_REMAP_BEGIN rz="<<rz<<" pressure_only="<<open<<" cycle="<<cycle
            <<" coarsen="<<(cycle%2)<<" source_cells="<<sourceMesh.cellIndexNum()<<" target_cells="<<targetMesh.cellIndexNum()<<'\n'<<std::flush;
        const auto old=accepted;const auto sourceU=old.values().normalVelocity(),sourceP=old.values().pressure();const auto sourceStocks=stocks(old.phaseGeometry().vofPublication());
        const auto td=selected(targetMesh,rz);const auto volumes=GeometricVofMeshTransfer::prepare(sourceMesh,old.phaseGeometry().vofPublication(),targetMesh,td,old.timeLayer());
        check(cycle%2?volumes.report().coarsenedTargets>0:volumes.report().refinedTargets>0,"selected Hydro AMR lacks actual Interface refinement/coarsening");
        const auto tl=HydroDofLayout::build(targetMesh,td);const auto flow=pressureFlow(targetMesh,tl);
        const auto remap=open?HydroVofMeshTransfer::preparePrescribedFlow(sourceMesh,old,targetMesh,volumes,config.liquidDensity,1,flow,remapConfig):
            HydroVofMeshTransfer::prepare(sourceMesh,old,targetMesh,volumes,config.liquidDensity,1,remapConfig);
        if(open&&rz&&cycle==0)pairedBoundarySourceRoles(sourceMesh,old,targetMesh,volumes,flow,remapConfig,remap);
        if(open){F19BackflowRemapContract::auditProjection(sourceMesh,old,targetMesh,remap);
            check(remap.candidateState().vofPrescribedFlow()->pressureDrivenReservoirs&&remap.report().fixedNormals==0
                &&!remap.candidateState().vofPrescribedFlow()->endpointVelocity&&remap.candidateState().vofPrescribedFlow()->prescribedVelocityDofs.empty(),
                "selected pressure AMR invented a fixed velocity or dropped mode");
            check(remap.projectionMultiplier()!=remap.candidateState().values().pressure(),"selected pressure remap substituted projection lambda for physical p");
            peakProjection=std::max(peakProjection,remap.report().projectionEnergyIdentityRelativeResidual);
            if(cycle==0){auto bad=flow;bad.pressureDrivenReservoirs=false;
                reject([&]{HydroVofMeshTransfer::preparePrescribedFlow(sourceMesh,old,targetMesh,volumes,2,1,bad,remapConfig);},"remap changed pressure-driven mode");
                bad=flow;bad.pressureOutlets.front().pressure+=.0625;
                reject([&]{HydroVofMeshTransfer::preparePrescribedFlow(sourceMesh,old,targetMesh,volumes,2,1,bad,remapConfig);},"remap changed physical pressure");
                bad=flow;const int portRow=tl.nativeFace(bad.pressureOutlets.front().velocityDof);
                bad.rulesByNativeRow[portRow]=bad.rulesByNativeRow[portRow]==Mode::PureLiquidReservoir?Mode::PureGasReservoir:Mode::PureLiquidReservoir;
                reject([&]{HydroVofMeshTransfer::preparePrescribedFlow(sourceMesh,old,targetMesh,volumes,2,1,bad,remapConfig);},"remap changed physical phase reservoir");
                bad=flow;unsigned endpointCalls=0;bad.endpointVelocity=[&](int,Clock){++endpointCalls;return Vector2D{};};
                reject([&]{HydroVofMeshTransfer::preparePrescribedFlow(sourceMesh,old,targetMesh,volumes,2,1,bad,remapConfig);},"pressure remap accepted endpoint velocity");
                check(endpointCalls==0,"pressure remap evaluated forbidden endpoint velocity");
                reject([&]{HydroVofMeshTransfer::boundaryNormalMeans(sourceMesh,old,targetMesh,selected(targetMesh,rz,true));},"remap accepted equal-volume exchanged solid mask");
            }
        }else verify(sourceMesh,old,targetMesh,remap,Density::CellAverageProjection);
        check(old.values().normalVelocity()==sourceU&&old.values().pressure()==sourceP&&stocks(old.phaseGeometry().vofPublication())==sourceStocks
            &&sameClock(old.timeLayer(),accepted.timeLayer()),"selected remap mutated actual accepted source");
        inspect(targetMesh,remap.candidateState());accepted=remap.candidateState();++adaptations;step(targetMesh,accepted,open,config);
    }
    std::cout<<"SELECTED_HYDRO_AMR_SCENE,"<<(rz?"RZ":"XY")<<",pressure_only,"<<open<<",actual_remaps,4,physical_steps,5,P2_rank_retained,1\n";
}
void runPressureOnly() {
    // Explicit bounded subset: the full selector and its closed-flow failure
    // remain separate. These are the identical pressure-driven scenes.
    for(bool rz:{false,true})scene(rz,true);
    check(scenes==2&&adaptations==8&&physicalSteps==10&&solidChecks>0&&wallChecks>0&&firstSolidIncidents>0&&guards==10,
        "selected pressure-only Hydro AMR missing actual subset coverage");
    std::cout<<std::setprecision(17)<<"SELECTED_DOMAIN_HYDRO_AMR_PRESSURE_ONLY_PASS scenes="<<scenes<<" remaps="<<adaptations<<" physical_steps="<<physicalSteps
        <<" solid_checks="<<solidChecks<<" wall_checks="<<wallChecks<<" solid_first_incidents="<<firstSolidIncidents<<" guards="<<guards
        <<" peak_stored_R="<<peakR<<" peak_physical_energy="<<peakEnergy<<" peak_open_projection_energy="<<peakProjection
        <<" full_selector_closed_power_issue_pending=1"<<'\n';
}
void run() {
    for(bool rz:{false,true})for(bool open:{false,true})scene(rz,open);
    check(scenes==4&&adaptations==16&&physicalSteps==20&&solidChecks>0&&wallChecks>0&&firstSolidIncidents>0&&guards==10,
        "selected Hydro AMR missing actual selected/refine/coarsen/pressure/wall/solid-first coverage");
    std::cout<<std::setprecision(17)<<"SELECTED_DOMAIN_HYDRO_AMR_PASS scenes="<<scenes<<" remaps="<<adaptations<<" physical_steps="<<physicalSteps
        <<" solid_checks="<<solidChecks<<" wall_checks="<<wallChecks<<" solid_first_incidents="<<firstSolidIncidents<<" guards="<<guards
        <<" peak_stored_R="<<peakR<<" peak_physical_energy="<<peakEnergy<<" peak_open_projection_energy="<<peakProjection<<'\n';
}
}


namespace SelectedDomainNoSlipAmrContract {
namespace Base=SelectedDomainHydroAmrContract;
namespace FV=LAMP::Geometry::Detail::FiniteVolume;
using Mode=GeometricVofAdvectorSlab::Boundary::Mode;
unsigned scenes=0,changes=0,physicalSteps=0,ownerChecks=0,guards=0,wallSpans=0;
unsigned radialSplits=0,axialSplits=0,radialMerges=0,axialMerges=0;
long double peakR=0,peakEnergy=0,peakProjection=0,peakWallArea=0;
bool internalWall(const NativeLeafFace& f) {
    // Test-owned geometry only; production consumes explicit physical rows.
    return (!f.horizontal&&f.center.x==1&&f.center.y<3)||(f.horizontal&&f.center.y==3&&f.center.x>1);
}
std::unique_ptr<QuadtreeMesh> nozzle(bool fine) {
    auto mesh=Base::nozzle(fine);
    if(fine){DictionaryMesh::RefineBox stem;stem.level=1;stem.min={.75,2};stem.max={1.25,3.5};
        auto shoulder=stem;shoulder.min={1,3};shoulder.max={2,3.5};mesh->applyRefineBoxes({stem,shoulder});}
    return mesh;
}
HydroVofPrescribedFlow flowFor(const QuadtreeMesh& mesh,const HydroDofLayout& l) {
    auto flow=Base::pressureFlow(mesh,l);
    for(int u=0;u<l.velocityCount();++u)if(internalWall(l.domain().face(mesh,l.nativeFace(u))))flow.stationaryNoSlipVelocityDofs.push_back(u);
    check(!flow.stationaryNoSlipVelocityDofs.empty(),"no-slip AMR lacks actual internal solid walls");return flow;
}
void inspect(const QuadtreeMesh& mesh,const HydroTwoPhaseState& state) {
    Base::inspect(mesh,state);const auto& l=state.values().layout();const auto& flow=*state.vofPrescribedFlow();const auto expected=flowFor(mesh,l);
    check(flow.stationaryNoSlipVelocityDofs==expected.stationaryNoSlipVelocityDofs&&flow.rulesByNativeRow==expected.rulesByNativeRow
        &&flow.pressureDrivenReservoirs&&flow.prescribedVelocityDofs.empty()&&!flow.endpointVelocity
        &&flow.pressureOutletBackflowModel==expected.pressureOutletBackflowModel&&flow.pressureOutlets.size()==expected.pressureOutlets.size(),
        "no-slip AMR lost explicit phase/wall/pressure/backflow metadata");
    for(std::size_t i=0;i<flow.pressureOutlets.size();++i)check(flow.pressureOutlets[i].velocityDof==expected.pressureOutlets[i].velocityDof
        &&flow.pressureOutlets[i].pressure==expected.pressureOutlets[i].pressure,"no-slip AMR changed physical pressure traction data");
    for(int u=0;u<l.velocityCount();++u){const auto f=l.domain().face(mesh,l.nativeFace(u));const bool marked=std::binary_search(flow.stationaryNoSlipVelocityDofs.begin(),flow.stationaryNoSlipVelocityDofs.end(),u);
        check(marked==internalWall(f)&&flow.normalVelocity[u]==0,"no-slip AMR changed physical wall subset or invented a fixed normal");
        if(marked)check(f.area>0&&flow.rulesByNativeRow[f.nativeRow]==Mode::Closed&&state.values().normalVelocity()[u]==0,
            "no-slip AMR mask includes port, axis, interior or nonzero normal");}
}
void step(const QuadtreeMesh& mesh,HydroTwoPhaseState& accepted,const HydroVofAcConfig& config) {
    const auto before=accepted;const auto oldU=before.values().normalVelocity(),oldP=before.values().pressure();const auto& source=before.phaseGeometry().vofPublication();
    const auto t=before.timeLayer();const Clock end{t.step+1,t.time(),t.time()+1./16384,Clock::Position::End};
    std::cerr<<"NO_SLIP_AMR_STEP_BEGIN rz="<<source.domain().axisymmetric()<<" next_step="<<end.step<<" native_cells="<<mesh.cellIndexNum()<<'\n'<<std::flush;
    const auto actual=F12DiagnosticRateRows::prepare(mesh,accepted,end,config,true);const auto& frozen=actual.frozenMidpoint();
    const auto& dual=frozen.material();const auto& volume=actual.transport();const auto& factor=frozen.viscosity().viscosity();
    constexpr long double massGate=256*std::numeric_limits<double>::epsilon();const auto& balance=dual.report();
    check(dual.usesConservativeVofMomentum()&&balance.maxRelativeContinuityDefect<=massGate&&balance.maxRelativeStoredMassContinuityDefect<=massGate
        &&balance.maxRelativePrimalContinuityDefect<=massGate&&balance.maxRelativeMassTransferDefect<=massGate&&balance.maxRelativeFluxCommutatorDefect<=massGate,
        "no-slip AMR changed original actual M/F/R gates");peakR=std::max(peakR,balance.maxRelativeStoredMassContinuityDefect);
    for(const auto& r:actual.report().stages){check(r.consistencyPower==0&&r.momentumPowerRelativeDefect<=1e-11&&r.energy.closureRelativeResidual<=1e-11L
        &&r.saddle.algebraicRelativeResidual<=1e-11&&r.saddle.freeMomentumDefect<=1e-11&&r.saddle.maxCellContinuityDefect<=1e-11,
        "no-slip AMR changed original physical Saddle/power/KE gates");peakEnergy=std::max(peakEnergy,r.energy.closureRelativeResidual);}
    const auto& walls=before.vofPrescribedFlow()->stationaryNoSlipVelocityDofs;
    check(factor.strain().stationaryNoSlipVelocityDofs()==walls&&!factor.strain().isStressCompatible()&&factor.matrix().nonZeros()>0
        &&actual.candidateState().vofPrescribedFlow()->stationaryNoSlipVelocityDofs==walls
        &&frozen.viscosity().geometry().vofIdentity()==volume.midpoint().identity()
        &&sameClock(frozen.viscosity().geometry().timeLayer(),volume.midpoint().timeLayer())
        && &factor.strain().metrics().divergence().matrix()==&dual.geometry().metrics().divergence().matrix(),
        "no-slip AMR physical Step lost actual wall/phase/full-clock/native D strain owner");++ownerChecks;
    const double dissipation=factor.dissipation(mesh,actual.solvedMidpointValues());
    check(std::isfinite(dissipation)&&dissipation>0,"no-slip AMR actual solved flow has no positive retained Gram dissipation");
    for(std::size_t q=0;q<factor.strain().points().size();++q){const auto& point=factor.strain().points()[q];const int ci=dual.layout().nativeCell(point.cellDof);
        const auto& inventory=volume.midpoint().plicCells()[ci].inventory;const auto capacity=volume.midpoint().cells()[ci].capacity;
        const double expected=static_cast<double>((inventory[0]/capacity)*config.liquidViscosity+(inventory[1]/capacity)*config.gasViscosity);
        check(factor.viscosity()[q]==expected,"no-slip AMR viscosity reused old stocks/phase coefficients");}
    for(bool half:{true,false}){const auto& endpoint=half?volume.midpoint():volume.endpoint();const auto& ledger=half?volume.report().boundaryHalf:volume.report().boundaryFull;
        long double capacity=0;for(int ci:source.domain().cellIndexByRow())capacity+=source.cells()[ci].capacity;
        for(unsigned phase=0;phase<2;++phase){FV::Dyadic change,in,out;
            for(std::size_t ci=0;ci<source.plicCells().size();++ci){change.add(endpoint.plicCells()[ci].inventory[phase]);change.add(source.plicCells()[ci].inventory[phase],-1);
                if(source.domain().cellRowByIndex()[ci]<0)check(endpoint.plicCells()[ci].inventory==std::array<long double,2>{0,0},"no-slip AMR physical step filled a solid");}
            for(const auto f:source.domain().faces(mesh)){if(f.outwardSign==0)continue;const auto& q=half?volume.faces()[f.nativeRow].half:volume.faces()[f.nativeRow].full;
                const auto positive=phase==0?q.positiveLiquid:q.positiveGas,negative=phase==0?q.negativeLiquid:q.negativeGas;
                out.add(f.negativeCell>=0?positive:negative);in.add(f.negativeCell>=0?negative:positive);
                if(internalWall(f))check(q.total==0&&q.liquid==0&&q.gas==0&&positive==0&&negative==0,"no-slip AMR internal wall acquired shared flux");}
            change.add(out);change.add(in,-1);F19BackflowRemapContract::bound(change.rounded(),capacity,1e-12L,"no-slip AMR independent phase Delta+Out-In");
            F19BackflowRemapContract::bound(ledger.inflow[phase]-in.rounded(),capacity,1e-12L,"no-slip AMR actual phase incoming ledger");
            F19BackflowRemapContract::bound(ledger.outflow[phase]-out.rounded(),capacity,1e-12L,"no-slip AMR actual phase outgoing ledger");}
        check(ledger.enabled&&ledger.sourceIdentity==source.identity()&&ledger.physicalRelativeConservation<=1e-12L
            &&ledger.massRelativeConservation<=1e-12L&&ledger.algebraicRelativeConservation<=1e-12L,"no-slip AMR original phase/mass/WY gates");}
    long double incoming=0,outgoing=0;const auto& l=before.values().layout();
    for(const auto& port:before.vofPrescribedFlow()->pressureOutlets){const auto f=l.domain().face(mesh,l.nativeFace(port.velocityDof));
        const double q=f.outwardSign*f.area*actual.candidateState().values().normalVelocity()[port.velocityDof];incoming+=std::max(0.,-q);outgoing+=std::max(0.,q);
        check(frozen.saddle().pressureOutletLoads()[port.velocityDof]==-f.outwardSign*f.area*port.pressure,"no-slip AMR changed actual physical pressure load");}
    check(incoming>0&&outgoing>0&&frozen.saddle().pressureGaugeCount()==0,"no-slip AMR lost real free-pressure supply or absolute reference");
    check(accepted.sameIdentity(before)&&accepted.values().normalVelocity()==oldU&&accepted.values().pressure()==oldP,
        "no-slip AMR physical preparation mutated accepted source");
    actual.commit(mesh,accepted);check(sameClock(accepted.timeLayer(),end),"no-slip AMR physical commit lost full clock");++physicalSteps;inspect(mesh,accepted);
}
void coverage(const QuadtreeMesh& sm,const HydroTwoPhaseState& source,const QuadtreeMesh& tm,const HydroTwoPhaseState& target,bool coarsen) {
    const auto& sl=source.values().layout();const auto& tl=target.values().layout();const auto& old=source.vofPrescribedFlow()->stationaryNoSlipVelocityDofs;
    const auto& next=target.vofPrescribedFlow()->stationaryNoSlipVelocityDofs;std::vector<FV::Dyadic> oldCovered(old.size()),newCovered(next.size());
    std::vector<long double> oldAreas(old.size()),newAreas(next.size());std::vector<unsigned> oldHits(old.size()),newHits(next.size());
    check(coarsen?old.size()>next.size():old.size()<next.size(),"no-slip AMR did not really split/coarsen physical wall DOFs");
    for(std::size_t i=0;i<old.size();++i){const auto sf=sl.domain().face(sm,sl.nativeFace(old[i]));const auto a=Base::ends(sm,sf);
        for(std::size_t j=0;j<next.size();++j){const auto tf=tl.domain().face(tm,tl.nativeFace(next[j]));const auto b=Base::ends(tm,tf);
            if(sf.horizontal!=tf.horizontal||sf.outwardSign!=tf.outwardSign||(sf.horizontal?sf.center.y!=tf.center.y:sf.center.x!=tf.center.x))continue;
            const double left=std::max(sf.horizontal?a[0].x:a[0].y,tf.horizontal?b[0].x:b[0].y),right=std::min(sf.horizontal?a[1].x:a[1].y,tf.horizontal?b[1].x:b[1].y);
            if(!(left<right))continue;FV::Dyadic length(right);length.add(left,-1);oldCovered[i].add(length);newCovered[j].add(length);
            const auto area=CartesianMeasure::face(right-left,sf.horizontal?left+(right-left)/2:sf.center.x,sl.domain().axisymmetric());
            check(area>0,"no-slip AMR positive wall span has zero actual physical measure");oldAreas[i]+=area;newAreas[j]+=area;++oldHits[i];++newHits[j];++wallSpans;}}
    const auto side=[&](const QuadtreeMesh& mesh,const HydroDofLayout& l,const auto& rows,const auto& lengths,const auto& areas,const auto& hits,bool oldSide){
        for(std::size_t i=0;i<rows.size();++i){const auto f=l.domain().face(mesh,l.nativeFace(rows[i]));const auto e=Base::ends(mesh,f);
            FV::Dyadic whole(f.horizontal?e[1].x:e[1].y);whole.add(f.horizontal?e[0].x:e[0].y,-1);
            check(lengths[i].compare(whole)==0,"no-slip AMR lost or duplicated exact physical wall span");
            const auto defect=std::abs(areas[i]-f.area)/f.area;check(defect<=1e-12L,"no-slip AMR wrong XY/RZ wall overlap measure");peakWallArea=std::max(peakWallArea,defect);
            if(hits[i]>1){if(oldSide){if(f.horizontal)++axialSplits;else ++radialSplits;}else{if(f.horizontal)++axialMerges;else ++radialMerges;}}}
    };
    side(sm,sl,old,oldCovered,oldAreas,oldHits,true);side(tm,tl,next,newCovered,newAreas,newHits,false);
}
template<class Action>void reject(Action&& action,const char* token) {bool caught=false;try{action();}catch(const std::invalid_argument& error){
    check(std::string(error.what()).find(token)!=std::string::npos,"no-slip AMR mask rejected outside intended physical metadata check");caught=true;}
    check(caught,"no-slip AMR accepted changed or invalid physical wall mask");++guards;}
void scene(bool rz) {
    ++scenes;auto coarse=nozzle(false),fine=nozzle(true);const auto d=Base::selected(*coarse,rz);
    const auto pub=GeometricVofInitialization::prepare(*coarse,d,[](Vector2D p){return std::max(p.x-1,p.y-4.125);});const auto phase=phaseFor(*coarse,pub);const auto& l=phase.layout();
    auto accepted=HydroTwoPhaseState::createVofPrescribedFlow(*coarse,phase,HydroNativeValues(l,std::vector<double>(l.velocityCount()),std::vector<double>(l.pressureCount())),flowFor(*coarse,l));
    HydroVofAcConfig config;config.liquidDensity=2;config.gasDensity=1;config.liquidViscosity=.01;config.gasViscosity=.005;config.surfaceTension=0;
    check(HydroVofStep::stabilityRatesPrescribedFlow(*coarse,accepted,config).flowRate==0,"no-slip AMR invented initial velocity");step(*coarse,accepted,config);
    HydroVofMeshTransferConfig remapConfig;remapConfig.momentumDensityRepresentation=Density::CellAverageProjection;
    for(int cycle=0;cycle<4;++cycle){const auto& sm=cycle%2?*fine:*coarse;const auto& tm=cycle%2?*coarse:*fine;const auto old=accepted;
        std::cerr<<"NO_SLIP_AMR_REMAP_BEGIN rz="<<rz<<" cycle="<<cycle<<" coarsen="<<(cycle%2)<<'\n'<<std::flush;
        const auto td=Base::selected(tm,rz);const auto transfer=GeometricVofMeshTransfer::prepare(sm,old.phaseGeometry().vofPublication(),tm,td,old.timeLayer());
        check(cycle%2?transfer.report().coarsenedTargets>0:transfer.report().refinedTargets>0,"no-slip AMR lacks actual Interface refinement/coarsening");
        const auto tl=HydroDofLayout::build(tm,td);const auto flow=flowFor(tm,tl);const auto oldMask=old.vofPrescribedFlow()->stationaryNoSlipVelocityDofs;
        const auto remap=HydroVofMeshTransfer::preparePrescribedFlow(sm,old,tm,transfer,2,1,flow,remapConfig);
        F19BackflowRemapContract::auditProjection(sm,old,tm,remap);coverage(sm,old,tm,remap.candidateState(),cycle%2);
        check(remap.projectionMultiplier()!=remap.candidateState().values().pressure(),"no-slip AMR substituted artificial lambda for physical pressure");
        peakProjection=std::max(peakProjection,remap.report().projectionEnergyIdentityRelativeResidual);inspect(tm,remap.candidateState());
        auto bad=flow;bad.stationaryNoSlipVelocityDofs.erase(bad.stationaryNoSlipVelocityDofs.begin());
        reject([&]{HydroVofMeshTransfer::preparePrescribedFlow(sm,old,tm,transfer,2,1,bad,remapConfig);},"adds/removes a physical stationary no-slip wall span");
        bad=flow;int extra=-1;for(int u=0;u<tl.velocityCount();++u){const auto f=td.face(tm,tl.nativeFace(u));
            if(f.outwardSign!=0&&f.area>0&&flow.rulesByNativeRow[f.nativeRow]==Mode::Closed&&!internalWall(f)&&f.center.x==tm.maxCoord.x&&f.center.y>4){extra=u;break;}}
        check(extra>=0,"no-slip AMR add-mask control lacks another genuine positive-area Closed wall");
        bad.stationaryNoSlipVelocityDofs.insert(std::lower_bound(bad.stationaryNoSlipVelocityDofs.begin(),bad.stationaryNoSlipVelocityDofs.end(),extra),extra);
        reject([&]{HydroVofMeshTransfer::preparePrescribedFlow(sm,old,tm,transfer,2,1,bad,remapConfig);},"adds/removes a physical stationary no-slip wall span");
        if(cycle==0){bad=flow;bad.stationaryNoSlipVelocityDofs.insert(bad.stationaryNoSlipVelocityDofs.begin(),bad.stationaryNoSlipVelocityDofs.front());
            reject([&]{HydroVofMeshTransfer::preparePrescribedFlow(sm,old,tm,transfer,2,1,bad,remapConfig);},"sorted and unique");
            bad=flow;const int port=flow.pressureOutlets.front().velocityDof;bad.stationaryNoSlipVelocityDofs.insert(std::lower_bound(bad.stationaryNoSlipVelocityDofs.begin(),bad.stationaryNoSlipVelocityDofs.end(),port),port);
            reject([&]{HydroVofMeshTransfer::preparePrescribedFlow(sm,old,tm,transfer,2,1,bad,remapConfig);},"positive-area Closed wall");}
        check(old.vofPrescribedFlow()->stationaryNoSlipVelocityDofs==oldMask&&accepted.sameIdentity(old),"no-slip AMR rejected candidate mutated old wall owner");
        accepted=remap.candidateState();++changes;step(tm,accepted,config);
    }
    std::cout<<"NO_SLIP_AMR_SCENE,"<<(rz?"RZ":"XY")<<",actual_remaps,4,physical_steps,5,physical_internal_wall_split_merge,1\n";
}
void run() {
    for(bool rz:{false,true})scene(rz);
    check(scenes==2&&changes==8&&physicalSteps==10&&ownerChecks==10&&guards==20&&wallSpans>0
        &&radialSplits>0&&axialSplits>0&&radialMerges>0&&axialMerges>0,"no-slip AMR coverage incomplete");
    std::cout<<std::setprecision(17)<<"SELECTED_DOMAIN_NO_SLIP_AMR_PASS scenes="<<scenes<<" remaps="<<changes<<" physical_steps="<<physicalSteps<<" actual_viscous_owners="<<ownerChecks
        <<" guards="<<guards<<" positive_wall_spans="<<wallSpans<<" radial_splits="<<radialSplits<<" axial_splits="<<axialSplits
        <<" radial_merges="<<radialMerges<<" axial_merges="<<axialMerges<<" peak_stored_R="<<peakR<<" peak_physical_energy="<<peakEnergy
        <<" peak_artificial_projection_energy="<<peakProjection<<" peak_actual_wall_area="<<peakWallArea<<" full_selected_closed_failure_pending=1\n";
}
}

namespace { namespace F12DiagnosticRateRows {
void run() {
    active=true;MaxwellCoupledContract::run();SelectedDomainNoSlipAmrContract::run();active=false;
    check(pairs==11&&maxwellPairs==1&&noSlipPairs==10&&lateFailurePairs==1&&lateProviderCalls>=2&&stages>=pairs&&rows>0,"F12 actual Maxwell/open/no-slip consumer coverage incomplete");
    std::cout<<"F12_DIAGNOSTIC_RATE_ROWS_PASS,paired_preparations,"<<pairs<<",actual_stages,"<<stages<<",used_rate_rows,"<<rows
        <<",Maxwell_pairs,"<<maxwellPairs<<",open_no_slip_pairs,"<<noSlipPairs<<",late_failure_pairs,"<<lateFailurePairs<<",late_provider_calls,"<<lateProviderCalls
        <<",requested_threads,8,actual_workers_min,"<<minimumWorkers<<",actual_workers_max,"<<maximumWorkers
        <<",complete_stage_physical_report_bits,exact,candidate_U_p_PLIC_stocks,exact,full_native_MXCSR_errno,exact,original_physical_oracles,retained\n";
}
} } // Reopen the same anonymous namespace as the F12 pair helper.
int main(int argc,char** argv)
{
    try{if(argc==2&&std::string(argv[1])=="--diagnostic-rate-rows-only"){F12DiagnosticRateRows::run();return 0;}
        if(argc==2&&std::string(argv[1])=="--selected-domain-amr-only"){SelectedDomainHydroAmrContract::run();return 0;}
        if(argc==2&&std::string(argv[1])=="--selected-domain-amr-pressure-only"){SelectedDomainHydroAmrContract::runPressureOnly();return 0;}
        if(argc==2&&std::string(argv[1])=="--selected-domain-no-slip-amr-only"){SelectedDomainNoSlipAmrContract::run();return 0;}
        if(argc==2&&std::string(argv[1])=="--frozen-matrix-assembly-plan-only"){Task37FrozenMatrixAssemblyPlan::run();return 0;}
        if(argc==2&&std::string(argv[1])=="--vof-regular-preparation-remap-only"){Task35VofRegularPreparationRemap::run();return 0;}
        if(argc>=2&&std::string(argv[1])=="--vof-ehd-boundary-output-only"){check(argc==3,"F6C output selector requires a unique D work directory");VofEhdBoundaryOutputContract::run(argv[2]);return 0;}
        if(argc>=2&&std::string(argv[1])=="--vof-ehd-boundary-rounding-output-only"){check(argc==3,"F20 output selector requires a unique D work directory");F20BoundaryRoundingOutputContract::run(argv[2]);return 0;}
        if(argc==2&&std::string(argv[1])=="--vof-backflow-remap-only"){F19BackflowRemapContract::run();return 0;}
        if(argc==2&&std::string(argv[1])=="--vof-open-remap-only"){VofOpenRemapContract::run();return 0;}
        if(argc>=2&&std::string(argv[1])=="--vof-ehd-output-only"){check(argc==3,"output selector requires a unique D work directory");VofEhdOutputContract::run(argv[2]);return 0;}
        if(argc==2&&std::string(argv[1])=="--maxwell-coupled-only"){MaxwellCoupledContract::run();return 0;}
        if(argc==2&&std::string(argv[1])=="--maxwell-geometry-reuse-only"){F23MaxwellGeometryReuseContract::run();return 0;}
        if(argc==2&&std::string(argv[1])=="--vof-viscosity-geometry-reuse-only"){F24VofViscosityGeometryReuseContract::run();return 0;}
        scenario(false,Density::CellAverageProjection);scenario(true,Density::SharpGeometry);
        std::cout<<"VOF fluid mesh transfer PASS remaps="<<remaps<<" steps="<<steps<<" assertions="<<assertions<<'\n';return 0;}
    catch(const std::exception& e){std::cerr<<"VOF fluid mesh transfer FAIL "<<e.what()<<'\n';return 1;}
}
