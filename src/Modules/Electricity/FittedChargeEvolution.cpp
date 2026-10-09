#include "Modules/Electricity/FittedChargeEvolution.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/NativeMeshIdentity.h"
#include "Core/Geometry/CellChartPoint.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace Quantum::Electricity {
struct FittedChargeState::Data
{
    NativeMeshIdentity meshIdentity;
    struct Component {
        std::uint64_t surfaceId = 0;
        std::shared_ptr<SurfaceChargeField> charge;
        std::vector<double> previousMaterialCharge;
    };
    std::vector<Component> components;
    FittedSharpDielectricResult electric;
    double previousChargeDt = 0.0, time = 0.0;
    int completedStep = 0;
    explicit Data(const QuadtreeMesh& mesh) : meshIdentity(mesh), components(1) {}
};
struct FittedChargeAccess
{
    static const FittedChargeState::Data& read(const FittedChargeState& state) {
        if (!state.data_) throw std::logic_error("fitted electrical time level has not been published");
        return *state.data_;
    }
    static auto create(const QuadtreeMesh& mesh) { return std::make_shared<FittedChargeState::Data>(mesh); }
    static SurfaceChargeSource source(const std::shared_ptr<FittedChargeState::Data>& data) {
        FittedChargeState temporary;temporary.data_=data;return temporary.source();
    }
    static FittedChargeState publish(std::shared_ptr<FittedChargeState::Data> data) {
        if (data->components.empty() || !data->electric.report.success || !data->electric.report.conductionCurrent
            || !data->electric.potential.source())
            throw std::logic_error("incomplete fitted electrical publication");
        for (const auto& component:data->components)
          if (!component.charge || !component.charge->geometry().usesSource(*data->electric.potential.source())
            || (component.previousMaterialCharge.empty() ? data->previousChargeDt != 0.0
                : !std::isfinite(data->previousChargeDt) || data->previousChargeDt <= 0.0
                    || component.previousMaterialCharge.size() != component.charge->transport().panelCharge().size()))
            throw std::logic_error("incomplete or mismatched fitted electrical publication");
        FittedChargeState result; result.data_=std::move(data); return result;
    }
    static FittedChargeMeshTransfer transfer(std::shared_ptr<FittedChargeState::Data> data) {
        FittedChargeMeshTransfer result; result.data_=std::move(data); return result;
    }
    static auto complete(const FittedChargeMeshTransfer& transfer) {
        if (!transfer.data_) throw std::logic_error("fitted charge mesh transfer has not been prepared");
        return std::make_shared<FittedChargeState::Data>(*transfer.data_);
    }
};
std::size_t FittedChargeState::componentCount() const { return FittedChargeAccess::read(*this).components.size(); }
std::uint64_t FittedChargeState::surfaceId(std::size_t c) const { return FittedChargeAccess::read(*this).components.at(c).surfaceId; }
std::shared_ptr<const SurfaceChargeField> FittedChargeState::charge(std::size_t c) const { return FittedChargeAccess::read(*this).components.at(c).charge; }
const std::vector<double>& FittedChargeState::previousMaterialCharge(std::size_t c) const { return FittedChargeAccess::read(*this).components.at(c).previousMaterialCharge; }
std::shared_ptr<const SurfaceChargeField> FittedChargeState::charge() const {
    if(componentCount()!=1)throw std::logic_error("single-surface charge accessor cannot select from a collection");
    return charge(0);
}
const FittedSharpDielectricResult& FittedChargeState::electric() const { return FittedChargeAccess::read(*this).electric; }
const std::vector<double>& FittedChargeState::previousMaterialCharge() const {
    if(componentCount()!=1)throw std::logic_error("single-surface history accessor cannot select from a collection");
    return previousMaterialCharge(0);
}
double FittedChargeState::previousChargeDt() const { return FittedChargeAccess::read(*this).previousChargeDt; }
double FittedChargeState::time() const { return FittedChargeAccess::read(*this).time; }
int FittedChargeState::completedStep() const { return FittedChargeAccess::read(*this).completedStep; }
void FittedChargeState::validateMesh(const QuadtreeMesh& mesh) const { FittedChargeAccess::read(*this).meshIdentity.validate(mesh); }
SurfaceChargeSource FittedChargeState::source() const {
    if(componentCount()!=1) {
        std::vector<SurfaceChargeSource> sources;
        std::vector<int> owners;
        SurfaceChargeSource result;result.axisymmetric=true;
        for(std::size_t s=0;s<componentCount();++s) {
            const auto field=charge(s);
            if(!field->geometry().usesNativeLeaves())throw std::logic_error("charge collection requires native geometry");
            sources.push_back(field->source());
            if(s==0) {
                owners.assign(sources.back().cellIntegratedCharge.size(),-1);
                result.cellIntegratedCharge.assign(owners.size(),0);
            }
            for(const auto& panel:field->geometry().panels()) {
                auto& owner=owners.at(panel.cellIndex);
                if(owner>=0 && owner!=static_cast<int>(s))
                    throw std::runtime_error("multiple charged surfaces share a cut cell; refine before advancing charge");
                owner=static_cast<int>(s);
            }
            for(std::size_t c=0;c<owners.size();++c)result.cellIntegratedCharge[c]+=sources.back().cellIntegratedCharge.at(c);
        }
        result.evaluate=[owner=data_,sources=std::move(sources),owners=std::move(owners)](const auto& frame) {
            const int s=owners.at(frame.nativeCellIndex);
            if(s<0)throw std::runtime_error("charge query has no owning surface");
            return sources.at(s).evaluate(frame);
        };
        return result;
    }
    auto result=charge()->source();
    result.evaluate=[owner=data_, evaluate=std::move(result.evaluate)](const LevelSetHermiteGeometry::Frame& frame) {
        return evaluate(frame);
    };
    return result;
}
std::shared_ptr<const SurfaceChargeField> FittedChargeMeshTransfer::charge() const {
    if (!data_) throw std::logic_error("fitted charge mesh transfer has not been prepared");
    if(data_->components.size()!=1)throw std::logic_error("single-surface mesh transfer cannot select from a collection");
    return data_->components[0].charge;
}
namespace {
double totalPanelCharge(const Quantum::Electricity::SurfaceChargeField& field)
{
    long double total=0.0L;
    for (double value : field.transport().panelCharge()) total+=value;
    const double result=static_cast<double>(total);
    if (!std::isfinite(result)) throw std::runtime_error("integrated surface charge overflow");
    return result;
}

void require(bool condition, const char* reason)
{
    if (!condition) throw std::runtime_error(reason);
}

void appendSectionReport(SurfaceChargeFieldStep& physical,const SurfaceChargeFieldStep& event)
{
    require(physical.transport.success && event.transport.success && event.sectionDetached
        && event.transport.integratedNormalCurrent==0 && event.transport.netBoundaryCurrent==0,
        "section closure must be a successful zero-current material operation");
    physical.transport.totalChargeAfter=event.transport.totalChargeAfter;
    physical.transport.conservationResidual+=event.transport.conservationResidual;
    physical.transport.relativeConservationResidual=physical.transport.conservationScale==0?0:
        std::abs(physical.transport.conservationResidual)/physical.transport.conservationScale;
    physical.transport.maxAbsDensity=event.transport.maxAbsDensity;
    physical.remapConservationResidual+=event.remapConservationResidual;
    physical.sectionDetached=true;
    physical.sectionRemapConservationResidual=event.remapConservationResidual;
}

double relaxationRate(double conductivity, double permittivity)
{
    require(std::isfinite(permittivity) && permittivity > 0.0
        && std::isfinite(conductivity) && conductivity >= 0.0,
        "charge relaxation requires finite positive permittivity and nonnegative conductivity");
    const double rate = conductivity / permittivity;
    require(std::isfinite(rate) && rate >= 0.0 && (conductivity == 0.0 || rate > 0.0),
        "charge relaxation rate overflowed or lost a nonzero rate");
    return rate;
}
void requireSameCandidatePanels(const SurfaceChargeField& charge,
    const FittedChargeMotion& motion)
{
    require(motion.geometry && motion.panels
        && motion.panels->usesSource(*motion.geometry)
        && charge.geometry().usesSource(*motion.geometry),
        "charge and supplied geometry must retain the same final scalar");
    const auto& a = charge.geometry().panels();
    const auto& b = motion.panels->panels();
    require(a.size() == b.size(),"charge and supplied geometry panel inventories disagree");
    if(!motion.geometry->usesNativeLeaves())
        require(charge.geometry().cellToPanel() == motion.panels->cellToPanel(),
            "charge and supplied legacy geometry panel maps disagree");
    for (std::size_t i = 0; i < a.size(); ++i) {
        require(a[i].cellIndex == b[i].cellIndex && a[i].physicalMeasure == b[i].physicalMeasure,
            "charge and supplied geometry have different native panel identities or RZ areas");
        if(motion.geometry->usesNativeLeaves())
            require(a[i].coordinateIsX==b[i].coordinateIsX && a[i].measure==b[i].measure
                && a[i].beginFrame.hasCellChartPoint && b[i].beginFrame.hasCellChartPoint
                && a[i].endFrame.hasCellChartPoint && b[i].endFrame.hasCellChartPoint
                && Geometry::cellChartEqual(a[i].beginFrame.cellChartPoint,b[i].beginFrame.cellChartPoint)
                && Geometry::cellChartEqual(a[i].endFrame.cellChartPoint,b[i].endFrame.cellChartPoint),
                "charge and supplied native graph intervals disagree");
    }
}
}

FittedChargeStepReport preflightFittedChargeStep(const FittedChargeParameters& parameters, double dt)
{
    require(parameters.field.axisymmetric && parameters.field.model == "leakyDielectric"
        && parameters.field.interfaceDiscretization == "sharpHermiteP3",
        "fitted charge evolution requires RZ leakyDielectric/sharpHermiteP3");
    require(std::isfinite(dt) && dt > 0.0, "conducting time step must be finite and positive");
    const double bulkTolerance = parameters.closureTolerance;
    require(std::isfinite(bulkTolerance) && bulkTolerance > 0.0,
        "conducting closure requires a finite positive chargeClosure tolerance");
    FittedChargeStepReport result;
    result.dt = dt;
    result.relaxationRate = std::max(relaxationRate(parameters.conductivity.negative,parameters.field.permittivityNegative),
        relaxationRate(parameters.conductivity.positive,parameters.field.permittivityPositive));
    result.relaxationCfl = dt * result.relaxationRate;
    require(std::isfinite(result.relaxationCfl)
        && (result.relaxationRate == 0.0 || result.relaxationCfl > 0.0),
        "explicit charge relaxation CFL overflowed or underflowed");
    if (result.relaxationCfl > 0.25)
    {
        std::ostringstream message;
        message << std::setprecision(17) << "explicit charge relaxation requires dt*max(conductivity/permittivity)<=0.25"
            << "; dt=" << dt << ", rate=" << result.relaxationRate
            << ", CFL=" << result.relaxationCfl << "; reduce dt";
        throw std::runtime_error(message.str());
    }
    return result;
}

FittedSharpDielectricResult solveFittedChargeField(
    const QuadtreeMesh& mesh, const FittedChargeProblem& problem,
    std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> snapshot,
    const LevelSetInterfacePanels& panels, const SurfaceChargeSource& charge,
    const FittedPotentialInput& previousPotential, FittedSharpDielectricSolver* repeatedSolver)
{
    const auto& electric = problem.parameters.field;
    const auto& conductivity = problem.parameters.conductivity;
    if (snapshot && snapshot->usesNativeLeaves())
    {
        FittedSharpDielectricSolver initialSolver;
        return (repeatedSolver ? repeatedSolver : &initialSolver)->solveNative(mesh,
            problem.graph,problem.boundaries,std::move(snapshot),panels,electric,previousPotential,charge,&conductivity);
    }
    if (repeatedSolver)
        return repeatedSolver->solve(mesh,problem.graph,problem.boundaries,
            std::move(snapshot),panels,CartesianGridScope::UniformLeaves,electric,previousPotential,charge,&conductivity);
    return solveFittedSharpDielectricElectrostatics(mesh,problem.graph,problem.boundaries,
        std::move(snapshot),panels,CartesianGridScope::UniformLeaves,electric,previousPotential,charge,&conductivity);
}

FittedChargeState initializeFittedChargeState(const QuadtreeMesh& mesh, const FittedChargeProblem& problem,
    std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> geometry, const LevelSetInterfacePanels& panels,
    const std::vector<Vector2D>& fixedContacts, const FittedPotentialInput& previousPotential, double initialDensity)
{
    require(problem.parameters.field.axisymmetric && problem.parameters.field.model=="leakyDielectric"
        && problem.parameters.field.interfaceDiscretization=="sharpHermiteP3",
        "fitted charge initialization requires RZ leakyDielectric/sharpHermiteP3");
    require(geometry && panels.usesSource(*geometry), "initial charge must use the supplied interface panels");
    auto next=FittedChargeAccess::create(mesh);
    SurfaceChargeGeometryOptions options;
    options.axisymmetric=problem.parameters.field.axisymmetric;
    options.boundaryContacts=fixedContacts;
    options.scope=CartesianGridScope::UniformLeaves;
    next->components[0].charge=geometry->usesNativeLeaves()
        ? std::make_shared<SurfaceChargeField>(mesh,panels.nativeChain(),initialDensity)
        : std::make_shared<SurfaceChargeField>(mesh,geometry,initialDensity,std::move(options));
    FittedChargeMotion motion; motion.geometry=geometry; motion.panels=&panels;
    requireSameCandidatePanels(*next->components[0].charge,motion);
    next->electric=solveFittedChargeField(mesh,problem,geometry,next->components[0].charge->geometry(),next->components[0].charge->source(),previousPotential);
    return FittedChargeAccess::publish(std::move(next));
}

FittedChargeState initializeFittedChargeState(const QuadtreeMesh& mesh, const FittedChargeProblem& problem,
    std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> geometry,
    const std::vector<FittedChargeComponentInitial>& components, const FittedPotentialInput& previousPotential)
{
    require(geometry && geometry->usesNativeLeaves() && !components.empty()
        && problem.parameters.field.axisymmetric && problem.parameters.field.model=="leakyDielectric"
        && problem.parameters.field.interfaceDiscretization=="sharpHermiteP3",
        "charge collection initialization requires native RZ leakyDielectric/sharpHermiteP3");
    auto next=FittedChargeAccess::create(mesh);next->components.resize(components.size());
    std::vector<const LevelSetInterfacePanels*> surfaces;
    for(std::size_t s=0;s<components.size();++s) {
        const auto& input=components[s];
        require(input.surfaceId!=0 && input.panels && input.panels->usesSource(*geometry),
            "charge collection requires nonzero identities and the supplied scalar");
        for(std::size_t j=0;j<s;++j)require(components[j].surfaceId!=input.surfaceId,"duplicate charged surface identity");
        auto& component=next->components[s];component.surfaceId=input.surfaceId;
        component.charge=std::make_shared<SurfaceChargeField>(mesh,input.panels->nativeChain(),input.initialDensity);
        FittedChargeMotion binding;binding.geometry=geometry;binding.panels=input.panels;
        requireSameCandidatePanels(*component.charge,binding);surfaces.push_back(&component.charge->geometry());
    }
    const auto source=FittedChargeAccess::source(next);
    FittedSharpDielectricSolver solver;
    next->electric=solver.solveNative(mesh,problem.graph,problem.boundaries,geometry,surfaces,
        problem.parameters.field,previousPotential,source,&problem.parameters.conductivity);
    return FittedChargeAccess::publish(std::move(next));
}

FittedChargeCollectionCandidate prepareFittedChargeCollectionStep(
    const QuadtreeMesh& mesh, const FittedChargeProblem& problem, const FittedChargeState& oldState,
    const std::vector<FittedChargeComponentMotion>& motions)
{
    oldState.validateMesh(mesh);
    require(!motions.empty() && motions.size()==oldState.componentCount(),
        "charge collection continuation cannot create or remove surfaces");
    const auto& first=motions.front().motion;const double dt=first.dt;
    FittedChargeCollectionCandidate result;result.report=preflightFittedChargeStep(problem.parameters,dt);
    const double timeScale=std::max({std::abs(first.sourceTime),std::abs(first.targetTime),dt});
    require(first.sourceTime==oldState.time() && std::isfinite(first.targetTime) && first.targetTime>first.sourceTime
        && std::abs(first.targetTime-first.sourceTime-dt)<=16*std::numeric_limits<double>::epsilon()*timeScale
        && oldState.completedStep()<std::numeric_limits<int>::max(),"invalid collection material clock");
    const auto& oldReport=oldState.electric().report;
    require(oldReport.success && oldReport.usedSuperLU && oldReport.usedSharpInterface
        && oldReport.usedSurfaceCharge && oldReport.conductionCurrent,"collection lacks the accepted global conduction current");
    const auto& current=*oldReport.conductionCurrent;
    const std::size_t cells=mesh.cellIndexNum();
    require(current.cellIntegratedCurrent.size()==cells && std::isfinite(current.absoluteFaceCurrent)
        && current.absoluteFaceCurrent>=0 && std::isfinite(current.netBoundaryCurrent),"invalid global phase current");
    const double tolerance=problem.parameters.closureTolerance;
    std::vector<int> owner(cells,-1),sourceIndex(motions.size(),-1);
    std::vector<std::vector<double>> currents(motions.size(),std::vector<double>(cells,0));
    for(std::size_t s=0;s<motions.size();++s) {
        for(std::size_t j=0;j<s;++j)require(motions[j].surfaceId!=motions[s].surfaceId,"duplicate target surface identity");
        for(std::size_t old=0;old<oldState.componentCount();++old)
            if(oldState.surfaceId(old)==motions[s].surfaceId)sourceIndex[s]=static_cast<int>(old);
        require(motions[s].surfaceId!=0 && sourceIndex[s]>=0,"material continuation changed surface identity");
        const auto& motion=motions[s].motion;
        require(motion.geometry && motion.geometry==first.geometry && motion.geometry->usesNativeLeaves()
            && motion.panels && motion.panels->usesSource(*motion.geometry) && motion.dt==dt
            && motion.sourceTime==first.sourceTime && motion.targetTime==first.targetTime
            && !motion.sectionDetachment && (bool(motion.forwardMaterialPoint)!=bool(motion.preparedTransfer)),
            "collection requires one native target/time and unchanged-topology material maps");
        if(motion.preparedTransfer) {
            const auto& transfer=*motion.preparedTransfer;
            require(transfer.target==motion.panels && motion.fixedContacts.empty()
                && transfer.step==oldState.completedStep()+1 && transfer.sourceTime==motion.sourceTime
                && transfer.targetTime==motion.targetTime && transfer.dt==dt,"collection material map has another time or target");
        }
        const auto charge=oldState.charge(sourceIndex[s]);
        require(charge && charge->geometry().usesNativeLeaves() && charge->geometry().axisymmetric(),"collection lost native RZ charge");
        for(const auto& panel:charge->geometry().panels()) {
            auto& cellOwner=owner.at(panel.cellIndex);
            require(cellOwner<0 || cellOwner==static_cast<int>(s),"different charged surfaces share an unresolved current cell");
            cellOwner=static_cast<int>(s);
        }
    }
    long double bulk=0,bulkAbs=0;
    for(std::size_t c=0;c<cells;++c) {
        const double value=current.cellIntegratedCurrent[c];require(std::isfinite(value),"nonfinite global cell current");
        if(owner[c]>=0)currents[owner[c]][c]=value;
        else {bulk+=value;bulkAbs+=std::abs(value);}
    }
    require(bulkAbs<=static_cast<long double>(tolerance)*current.absoluteFaceCurrent,
        "global non-interface conduction current exceeds bulk roundoff allowance");
    result.unassignedBulkCurrent=static_cast<double>(bulk);result.absoluteBulkCurrent=static_cast<double>(bulkAbs);
    auto next=FittedChargeAccess::create(mesh);next->components.resize(motions.size());
    result.chargeSteps.resize(motions.size());std::vector<const LevelSetInterfacePanels*> surfaces;
    const auto started=std::chrono::steady_clock::now();long double residual=0,scale=0;
    for(std::size_t s=0;s<motions.size();++s) {
        const auto& motion=motions[s].motion;const auto old=oldState.charge(sourceIndex[s]);
        auto& component=next->components[s];component.surfaceId=motions[s].surfaceId;
        auto& report=result.chargeSteps[s];
        component.charge=motion.preparedTransfer
            ? old->advancePrepared(mesh,*motion.preparedTransfer,currents[s],current.absoluteFaceCurrent,tolerance,dt,report)
            : old->advanceMoving(mesh,motion.geometry,motion.fixedContacts,motion.forwardMaterialPoint,
                currents[s],current.absoluteFaceCurrent,tolerance,dt,report,motion.panels);
        require(component.charge && report.transport.success && report.remapped,"child charge material continuation failed");
        requireSameCandidatePanels(*component.charge,motion);surfaces.push_back(&component.charge->geometry());
        residual+=static_cast<long double>(report.transport.totalChargeAfter)-report.transport.totalChargeBefore;
        scale+=report.transport.conservationScale;
    }
    result.report.materialTransfer_s=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
    residual+=static_cast<long double>(dt)*current.netBoundaryCurrent;
    scale+=static_cast<long double>(dt)*current.absoluteFaceCurrent;
    require(std::isfinite(residual) && std::isfinite(scale) && scale>=0 && std::abs(residual)<=tolerance*scale,
        "collection charge increment does not close with the single old electrode current");
    result.report.chargeClosureResidual=static_cast<double>(residual);
    result.report.chargeClosureScale=static_cast<double>(scale);
    result.report.relativeChargeClosure=scale>0?static_cast<double>(std::abs(residual)/scale):0;
    const auto source=FittedChargeAccess::source(next);FittedSharpDielectricSolver solver;
    next->electric=solver.solveNative(mesh,problem.graph,problem.boundaries,first.geometry,surfaces,
        problem.parameters.field,oldState.electric().potential,source,&problem.parameters.conductivity);
    result.report.electricTimings=next->electric.timings;
    result.report.electricMatrix_s=next->electric.linearSolve.timeAssemble_s;
    result.report.electricFactor_s=next->electric.linearSolve.timeFactor_s;
    result.report.electricBacksolve_s=next->electric.linearSolve.timeSolve_s;
    next->completedStep=oldState.completedStep()+1;next->time=first.targetTime;
    result.state=FittedChargeAccess::publish(std::move(next));return result;
}

FittedChargeCollectionCandidate prepareFittedChargeInteriorSplitStep(
    const QuadtreeMesh& mesh,const FittedChargeProblem& problem,const FittedChargeState& oldState,
    const FittedChargeMotion& physical,const FittedChargeInteriorSplit& split)
{
    oldState.validateMesh(mesh);const double dt=physical.dt;
    FittedChargeCollectionCandidate result;result.report=preflightFittedChargeStep(problem.parameters,dt);
    require(oldState.componentCount()==1 && physical.preparedTransfer && !physical.forwardMaterialPoint
        && !physical.sectionDetachment && physical.fixedContacts.empty() && physical.geometry && physical.panels
        && physical.panels->usesSource(*physical.geometry) && physical.geometry->usesNativeLeaves()
        && split.geometry && split.geometry->usesNativeLeaves() && split.targets.size()==2,
        "interior charge event requires one connected native physical leg and two children");
    const auto& transfer=*physical.preparedTransfer;
    require(oldState.completedStep()<std::numeric_limits<int>::max() && transfer.step==oldState.completedStep()+1
        && transfer.sourceTime==oldState.time() && physical.sourceTime==transfer.sourceTime
        && transfer.targetTime==physical.targetTime && transfer.dt==dt && transfer.target==physical.panels,
        "interior charge event has another material source, target or clock");
    require(oldState.electric().report.success && oldState.electric().report.usedSuperLU
        && oldState.electric().report.conductionCurrent,"interior charge event has no old global current");
    const auto& current=*oldState.electric().report.conductionCurrent;
    const double tolerance=problem.parameters.closureTolerance;
    const auto started=std::chrono::steady_clock::now();SurfaceChargeFieldStep physicalReport;
    auto connected=oldState.charge()->advancePrepared(mesh,transfer,current.cellIntegratedCurrent,
        current.absoluteFaceCurrent,tolerance,dt,physicalReport);
    requireSameCandidatePanels(*connected,physical);
    auto next=FittedChargeAccess::create(mesh);next->components.resize(2);
    std::vector<const LevelSetInterfacePanels*> surfaces;std::vector<double> measures;
    std::array<std::size_t,3> offset{};
    for(std::size_t s=0;s<2;++s) {
        const auto& target=split.targets[s];
        require(target.surfaceId!=0 && target.surfaceId!=oldState.surfaceId(0) && target.panels
            && target.panels->usesSource(*split.geometry) && target.panels->usesNativeLeaves()
            && !target.panels->closed() && target.panels->axisymmetric()
            && (s==0 || target.surfaceId!=split.targets[0].surfaceId),"invalid interior child charge identity/geometry");
        auto& child=next->components[s];child.surfaceId=target.surfaceId;
        child.charge=std::make_shared<SurfaceChargeField>(mesh,*target.panels,0);
        offset[s]=measures.size();const auto& area=child.charge->transport().panelMeasure();
        measures.insert(measures.end(),area.begin(),area.end());surfaces.push_back(&child.charge->geometry());
    }
    offset[2]=measures.size();const auto& before=connected->geometry().panels();
    require(split.weights.size()==before.size() && std::isfinite(split.section),"incomplete interior material rows");
    for(std::size_t i=0;i<split.cells.size();++i) {
        require(split.cells[i]>=0 && split.cells[i]<mesh.cellIndexNum(),"invalid interior event cell");
        for(std::size_t j=0;j<i;++j)require(split.cells[i]!=split.cells[j],"duplicate interior event cell");
        const auto box=physical.geometry->nativeCellBounds(split.cells[i]);
        require(box[0].x==0 && (i<2?box[1].y<=split.section:box[0].y>=split.section),
            "interior event cells do not follow the axis/section partition");
    }
    const auto local=[&](int c){return std::find(split.cells.begin(),split.cells.end(),c)!=split.cells.end();};
    std::vector<std::vector<SurfaceChargeTransport::PanelWeight>> weights(before.size());
    for(std::size_t i=0;i<before.size();++i) {
        const auto& a=before[i];const auto box=physical.geometry->nativeCellBounds(a.cellIndex);
        require(box[1].y<=split.section || box[0].y>=split.section,"source material panel crosses the event section");
        const bool above=box[0].y>=split.section;
        for(const auto& weight:split.weights[i]) {
            require(weight.targetSurface<2 && weight.targetPanel<surfaces[weight.targetSurface]->panels().size(),
                "interior material weight has no child panel");
            const auto& b=surfaces[weight.targetSurface]->panels()[weight.targetPanel];
            const auto targetBox=split.geometry->nativeCellBounds(b.cellIndex);
            require(above?targetBox[0].y>=split.section:targetBox[1].y<=split.section,
                "interior charge crosses its declared material section");
            if(local(a.cellIndex))require(local(b.cellIndex),"local neck charge escaped its event region");
            else require(!local(b.cellIndex) && split.weights[i].size()==1 && weight.share==1
                && a.cellIndex==b.cellIndex && a.coordinateIsX==b.coordinateIsX
                && a.measure==b.measure && a.physicalMeasure==b.physicalMeasure
                && Geometry::cellChartEqual(a.beginFrame.cellChartPoint,b.beginFrame.cellChartPoint)
                && Geometry::cellChartEqual(a.endFrame.cellChartPoint,b.endFrame.cellChartPoint),
                "interior event changed remote material charge identity");
            weights[i].push_back({offset[weight.targetSurface]+weight.targetPanel,weight.share});
        }
    }
    auto material=SurfaceChargeTransport::fromIntegratedCharge(connected->transport().panelMeasure(),
        connected->transport().panelCharge(),false);
    const auto event=material.remapWeights(weights,measures);require(event.success,event.reason.c_str());
    result.chargeSteps.resize(2);
    for(std::size_t s=0;s<2;++s) {
        next->components[s].charge->assignIntegratedCharge(std::vector<double>(material.panelCharge().begin()+offset[s],
            material.panelCharge().begin()+offset[s+1]));
        auto& report=result.chargeSteps[s];report.remapped=true;report.transport.success=true;
        report.transport.stage="interior-split-child";report.transport.totalChargeAfter=totalPanelCharge(*next->components[s].charge);
    }
    result.interiorSplit=true;result.eventRemapResidual=event.conservationResidual;
    result.unassignedBulkCurrent=physicalReport.unassignedBulkCurrent;result.absoluteBulkCurrent=physicalReport.absoluteBulkCurrent;
    const long double residual=static_cast<long double>(event.totalChargeAfter)-physicalReport.transport.totalChargeBefore
        +static_cast<long double>(dt)*current.netBoundaryCurrent;
    const long double scale=physicalReport.transport.conservationScale+static_cast<long double>(dt)*current.absoluteFaceCurrent;
    require(std::isfinite(residual) && std::isfinite(scale) && scale>=0 && std::abs(residual)<=tolerance*scale,
        "interior split charge does not close with the one physical electrode current");
    result.report.chargeClosureResidual=static_cast<double>(residual);result.report.chargeClosureScale=static_cast<double>(scale);
    result.report.relativeChargeClosure=scale>0?static_cast<double>(std::abs(residual)/scale):0;
    result.report.materialTransfer_s=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
    const auto source=FittedChargeAccess::source(next);FittedSharpDielectricSolver solver;
    next->electric=solver.solveNative(mesh,problem.graph,problem.boundaries,split.geometry,surfaces,
        problem.parameters.field,oldState.electric().potential,source,&problem.parameters.conductivity);
    result.report.electricTimings=next->electric.timings;result.report.electricMatrix_s=next->electric.linearSolve.timeAssemble_s;
    result.report.electricFactor_s=next->electric.linearSolve.timeFactor_s;result.report.electricBacksolve_s=next->electric.linearSolve.timeSolve_s;
    next->completedStep=oldState.completedStep()+1;next->time=physical.targetTime;
    result.state=FittedChargeAccess::publish(std::move(next));return result;
}

FittedChargeMeshTransfer prepareFittedChargeMeshTransfer(const QuadtreeMesh& sourceMesh,
    const FittedChargeState& state, const QuadtreeMesh& target,
    std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> geometry, const std::vector<int>& sourceToTargetCell,
    const LevelSetInterfacePanels* transferredPanels)
{
    state.validateMesh(sourceMesh);
    std::vector<std::vector<double>> histories;
    if (!state.previousMaterialCharge().empty()) histories.push_back(state.previousMaterialCharge());
    auto mapped=state.charge()->rebindMesh(target,std::move(geometry),sourceToTargetCell,histories,transferredPanels);
    auto next=FittedChargeAccess::create(target);
    next->components[0].surfaceId=state.surfaceId(0);
    next->components[0].charge=std::move(mapped.field);
    if (!histories.empty()) next->components[0].previousMaterialCharge=std::move(mapped.history.front());
    next->previousChargeDt=state.previousChargeDt();
    next->time=state.time(); next->completedStep=state.completedStep();
    require(state.charge()->transport().panelCharge()==next->components[0].charge->transport().panelCharge()
        && next->components[0].previousMaterialCharge==state.previousMaterialCharge(),
        "open-chain AMR changed current charge or material history");
    return FittedChargeAccess::transfer(std::move(next));
}

FittedChargeState solveFittedChargeMeshTransfer(const QuadtreeMesh& mesh, const FittedChargeProblem& problem,
    const FittedChargeMeshTransfer& transfer, std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> geometry,
    const FittedPotentialInput& previousPotential)
{
    auto next=FittedChargeAccess::complete(transfer);
    next->meshIdentity.validate(mesh);
    require(geometry && next->components[0].charge->geometry().usesSource(*geometry), "transferred electric source mismatch");
    next->electric=solveFittedChargeField(mesh,problem,geometry,next->components[0].charge->geometry(),next->components[0].charge->source(),previousPotential);
    return FittedChargeAccess::publish(std::move(next));
}

FittedChargeCandidate prepareFittedChargeStep(
    const QuadtreeMesh& mesh, const FittedChargeProblem& problem, const FittedChargeState& oldState,
    const FittedChargeMotion& motion, const SurfaceChargeRelaxationOptions* relaxation,
    const FittedChargeObserver& observer)
{
    const double dt = motion.dt;
    const char* stage = "preflight";
    try
    {
        const auto preflight = preflightFittedChargeStep(problem.parameters,dt);
        stage = "time-level-input";
        oldState.validateMesh(mesh);
        const double timeScale=std::max({std::abs(motion.sourceTime),std::abs(motion.targetTime),dt});
        require(motion.sourceTime==oldState.time() && std::isfinite(motion.targetTime)
            && motion.targetTime>motion.sourceTime
            && std::abs((motion.targetTime-motion.sourceTime)-dt)<=16.0*std::numeric_limits<double>::epsilon()*timeScale
            && oldState.completedStep()<std::numeric_limits<int>::max(),
            "charge motion times do not match the retained source and numerical dt");
        require(oldState.charge()
            && oldState.charge()->geometry().axisymmetric() && oldState.charge()->geometry().hasRetainedSource(),
            "old conducting charge must own its retained active RZ panels");
        require(motion.targetSurfaceId==0 || motion.targetSurfaceId==oldState.surfaceId(0) || motion.sectionDetachment,
            "single-surface charge identity changed without a material event");
        const auto& oldReport = oldState.electric().report;
        require(oldReport.success && oldReport.usedSuperLU && oldReport.usedSharpInterface
            && oldReport.usedSurfaceCharge && oldReport.conductionCurrent,
            "old conducting time level lacks its successful actual P3/SuperLU phase current");
        require(motion.geometry && motion.panels && (bool(motion.forwardMaterialPoint)!=bool(motion.preparedTransfer)),
            "conducting candidate requires the matching final interface and its retained material trajectory");
        if (motion.preparedTransfer) {
            const auto& transfer=*motion.preparedTransfer;
            require(transfer.target==(motion.sectionDetachment?motion.sectionDetachment->source:motion.panels)
                && motion.panels->usesSource(*motion.geometry)
                && motion.fixedContacts.empty() && transfer.step==oldState.completedStep()+1
                && transfer.sourceTime==motion.sourceTime && transfer.targetTime==motion.targetTime && transfer.dt==dt,
                "prepared charge transfer has another target, boundary law or physical step");
        }
        const auto& oldCurrent = *oldReport.conductionCurrent;
        if(motion.sectionDetachment)
            require(motion.preparedTransfer && motion.sectionDetachment->target==motion.panels
                && motion.sectionDetachment->step==oldState.completedStep()+1
                && motion.sectionDetachment->time==motion.targetTime,
                "section detachment lacks its preceding physical transfer or matching completed clock");
        require(oldCurrent.cellIntegratedCurrent.size() == static_cast<std::size_t>(mesh.cellIndexNum())
            && std::isfinite(oldCurrent.absoluteFaceCurrent) && oldCurrent.absoluteFaceCurrent >= 0.0
            && std::isfinite(oldCurrent.netBoundaryCurrent),
            "old phase-current report has invalid native cell storage or boundary totals");
        const double closureTolerance = problem.parameters.closureTolerance;
        // Interface's unaccepted geometric candidate already contains the complete
        // native panel integration. Electricity checks its scalar identity
        // and owns a copy; no second reconstruction of this time level.
        const auto* preparedPanels=motion.geometry->usesNativeLeaves()
            ? motion.panels : nullptr;
        FittedChargeCandidate result;
        auto next=FittedChargeAccess::create(mesh);
        next->components[0].surfaceId=motion.targetSurfaceId?motion.targetSurfaceId:oldState.surfaceId(0);
        result.report = preflight;
        if (relaxation)
        {
            stage="implicit-material-history-remap";
            const bool hasHistory=!oldState.previousMaterialCharge().empty();
            require(hasHistory ? oldState.previousChargeDt()>0.0 : oldState.previousChargeDt()==0.0,
                "implicit charge history and previous dt must be retained together");
            const auto materialStarted=std::chrono::steady_clock::now();
            auto remapped=motion.preparedTransfer ? oldState.charge()->remapPrepared(mesh,*motion.preparedTransfer,
                hasHistory ? std::vector<std::vector<double>>{oldState.previousMaterialCharge()}
                           : std::vector<std::vector<double>>{}) : oldState.charge()->remapMoving(mesh,motion.geometry,
                motion.fixedContacts,motion.forwardMaterialPoint,
                hasHistory ? std::vector<std::vector<double>>{oldState.previousMaterialCharge()}
                           : std::vector<std::vector<double>>{},preparedPanels);
            if(motion.sectionDetachment) {
                auto closed=remapped.field->remapSectionDetachment(mesh,*motion.sectionDetachment,remapped.history);
                appendSectionReport(remapped.report,closed.report);
                remapped.field=std::move(closed.field);remapped.history=std::move(closed.history);
                remapped.historyReports=std::move(closed.historyReports);
            }
            result.report.materialTransfer_s=std::chrono::duration<double>(
                std::chrono::steady_clock::now()-materialStarted).count();
            next->components[0].charge=std::move(remapped.field);
            result.chargeStep=std::move(remapped.report);
            requireSameCandidatePanels(*next->components[0].charge,motion);
            const auto& history=next->components[0].charge->transport();
            next->components[0].previousMaterialCharge=history.panelCharge();
            next->previousChargeDt=dt;
            FittedSharpDielectricSolver electricSolver; // One owner, one candidate geometry.
            int evaluations=0;
            stage="implicit-electric-charge-iteration";
            const auto relaxed=Quantum::Electricity::relaxSurfaceCharge(history,
                hasHistory ? &remapped.history.front() : nullptr,dt,oldState.previousChargeDt(),
                [&](const std::vector<double>& trial) {
                    const auto trialSource=next->components[0].charge->source(trial);
                    auto electric=solveFittedChargeField(mesh,problem,motion.geometry,
                        next->components[0].charge->geometry(),trialSource,
                        evaluations ? next->electric.potential : oldState.electric().potential,&electricSolver);
                    require(electric.report.success && electric.report.conductionCurrent,
                        "implicit charge requires the actual successful P3 phase current");
                    const auto& current=*electric.report.conductionCurrent;
                    auto panelCurrent=next->components[0].charge->mapNormalCurrent(current.cellIntegratedCurrent,
                        current.absoluteFaceCurrent,closureTolerance,result.chargeStep);
                    result.report.electricFactorReuses+=electric.linearSolve.fittedFactorizationReused ? 1 : 0;
                    result.report.electricGeometryReuses+=electric.geometryAssemblyReused ? 1 : 0;
                    result.report.electricPreparedRhsSolves += electric.linearSolve.fittedRhsOnlySolve ? 1 : 0;
                    auto& timing = result.report.electricTimings;
                    timing.input_s += electric.timings.input_s;
                    timing.geometry_s += electric.timings.geometry_s;
                    timing.linear_s += electric.timings.linear_s;
                    timing.flux_s += electric.timings.flux_s;
                    timing.current_s += electric.timings.current_s;
                    result.report.electricMatrix_s += electric.linearSolve.timeAssemble_s;
                    result.report.electricFactor_s += electric.linearSolve.timeFactor_s;
                    result.report.electricBacksolve_s += electric.linearSolve.timeSolve_s;
                    ++evaluations;
                    if (observer) observer(evaluations,relaxation->maxIterations,
                        electric.geometryAssemblyReused,electric.linearSolve.fittedFactorizationReused);
                    next->electric=std::move(electric);
                    return panelCurrent;
                },*relaxation);
            result.report.relaxation=relaxed.report;
            if (!relaxed.report.success)
            {
                std::ostringstream message;
                message << "electric/charge iteration failed at " << relaxed.report.stage
                    << "; iterations=" << relaxed.report.iterations
                    << "; relative_equation_residual=" << relaxed.report.relativeEquationResidual
                    << "; reason=" << relaxed.report.reason;
                throw std::runtime_error(message.str());
            }
            stage="implicit-electrode-charge-closure";
            double totalCurrent=0.0, scale=0.0;
            const double historyWeight=hasHistory ? relaxed.report.beta*dt/oldState.previousChargeDt() : 0.0;
            for (std::size_t i=0; i<relaxed.integratedCharge.size(); ++i)
            {
                totalCurrent+=relaxed.integratedNormalCurrent[i];
                scale+=relaxed.report.alpha*std::max(std::abs(relaxed.integratedCharge[i]),std::abs(history.panelCharge()[i]));
                if (hasHistory) scale+=std::abs(historyWeight*(history.panelCharge()[i]-remapped.history.front()[i]));
            }
            const auto& current=*next->electric.report.conductionCurrent;
            scale+=dt*current.absoluteFaceCurrent;
            // The kernel's sum is the actual BDF derivative minus panel I.
            // Add the same-iterate panel and electrode currents to test the
            // complete domain equation; do not substitute old-current Euler.
            const double residual=relaxed.report.integratedEquationResidual
                +dt*(totalCurrent+current.netBoundaryCurrent);
            require(std::isfinite(residual) && std::isfinite(scale) && scale>=0.0
                && std::abs(residual)<=closureTolerance*scale,
                "implicit charge time derivative does not close with new electrode conduction current");
            result.report.chargeClosureResidual=residual;
            result.report.chargeClosureScale=scale;
            result.report.relativeChargeClosure=scale>0.0 ? std::abs(residual)/scale : 0.0;
            auto& transport=result.chargeStep.transport;
            transport.totalChargeBefore=totalPanelCharge(*oldState.charge());
            transport.integratedNormalCurrent=totalCurrent;
            transport.conservationResidual=relaxed.report.integratedEquationResidual;
            transport.conservationScale=scale;
            transport.relativeConservationResidual=scale>0.0 ? std::abs(transport.conservationResidual)/scale : 0.0;
            next->components[0].charge->assignIntegratedCharge(relaxed.integratedCharge);
            transport.totalChargeAfter=totalPanelCharge(*next->components[0].charge);
            transport.maxAbsDensity=0.0;
            for (double density : next->components[0].charge->transport().density())
                transport.maxAbsDensity=std::max(transport.maxAbsDensity,std::abs(density));
            transport.success=true; transport.stage="implicit-material-charge";
            transport.reason="material histories and coeval P3 current satisfy the implicit charge equation";
        }
        else
        {
            stage = "old-current-and-material-remap";
            const auto materialStarted=std::chrono::steady_clock::now();
            next->components[0].charge = motion.preparedTransfer ? oldState.charge()->advancePrepared(mesh,*motion.preparedTransfer,
                oldCurrent.cellIntegratedCurrent,oldCurrent.absoluteFaceCurrent,closureTolerance,dt,result.chargeStep)
                : oldState.charge()->advanceMoving(mesh,motion.geometry,
                motion.fixedContacts,motion.forwardMaterialPoint,oldCurrent.cellIntegratedCurrent,
                oldCurrent.absoluteFaceCurrent,closureTolerance,dt,result.chargeStep,preparedPanels);
            if(motion.sectionDetachment) {
                auto closed=next->components[0].charge->remapSectionDetachment(mesh,*motion.sectionDetachment);
                appendSectionReport(result.chargeStep,closed.report);
                next->components[0].charge=std::move(closed.field);
            }
            result.report.materialTransfer_s=std::chrono::duration<double>(
                std::chrono::steady_clock::now()-materialStarted).count();
            require(next->components[0].charge && result.chargeStep.transport.success && result.chargeStep.remapped,
                "surface charge material candidate did not complete its one current update and conservative remap");
            requireSameCandidatePanels(*next->components[0].charge,motion);

            stage = "electrode-charge-closure";
            const auto& charge = result.chargeStep.transport;
            // Preserve the existing EHD whole-domain accounting and arithmetic.
            // Net Q may vanish for signed charge, so the denominator includes
            // absolute old charge/current, not |totalChargeBefore| alone.
            const double residual = charge.totalChargeAfter - charge.totalChargeBefore
                + dt * oldCurrent.netBoundaryCurrent;
            const double scale = charge.conservationScale + dt * oldCurrent.absoluteFaceCurrent;
            const double allowed = closureTolerance * scale;
            if (!std::isfinite(residual) || !std::isfinite(scale) || scale < 0.0
                || !std::isfinite(allowed) || std::abs(residual) > allowed)
            {
                std::ostringstream message;
                message << std::setprecision(17) << "surface charge increment does not close with old electrode conduction current"
                    << "; Q_before=" << charge.totalChargeBefore << ", Q_after=" << charge.totalChargeAfter
                    << ", dt=" << dt << ", old_boundary_J=" << oldCurrent.netBoundaryCurrent
                    << ", residual=" << residual << ", scale=" << scale << ", tolerance=" << closureTolerance;
                throw std::runtime_error(message.str());
            }
            result.report.chargeClosureResidual = residual;
            result.report.chargeClosureScale = scale;
            result.report.relativeChargeClosure = scale > 0.0 ? std::abs(residual) / scale : 0.0;
            require(std::isfinite(result.report.relativeChargeClosure), "nonfinite relative electrode-charge closure");
        }

        // source() borrows the candidate Field. Its owner stays in result
        // through BOTH synchronous consumers and is then returned intact;
        // never save this callback beside a movable/destroyed local Field.
        stage = "candidate-charge-reconstruction";
        const auto source = next->components[0].charge->source();
        stage = "candidate-charged-electrostatics";
        if (!relaxation)
        {
            next->electric = solveFittedChargeField(mesh,problem,motion.geometry,
                next->components[0].charge->geometry(),source,oldState.electric().potential);
            result.report.electricTimings=next->electric.timings;
            result.report.electricMatrix_s=next->electric.linearSolve.timeAssemble_s;
            result.report.electricFactor_s=next->electric.linearSolve.timeFactor_s;
            result.report.electricBacksolve_s=next->electric.linearSolve.timeSolve_s;
        }
        next->completedStep=oldState.completedStep()+1;
        next->time=motion.targetTime;
        result.state=FittedChargeAccess::publish(std::move(next));
        return result;
    }
    catch (const std::exception& error)
    {
        throw std::runtime_error(std::string("[FittedChargeEvolution] ") + stage + ": " + error.what());
    }
}
}
