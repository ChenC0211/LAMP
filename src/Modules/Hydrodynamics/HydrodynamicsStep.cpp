#include "Modules/Hydrodynamics/Hydrodynamics.h"
#include "Modules/Hydrodynamics/HydroConservativeTransport.h"
#include "Modules/Hydrodynamics/HydroFittedDomain.h"
#include "Modules/Interface/Compatibility/InterfaceParametricMotion.h"
#include "Core/Mesh/NativeMeshIdentity.h"
#include "Core/Mesh/NativeLeafDomain.h"
#include "Core/Mesh/QuadtreeMesh.h"

#include <cmath>
#include <limits>
#include <stdexcept>

namespace
{
bool sameTime(double a, double b, double dt)
{
    return std::abs(a - b) <= 64.0 * std::numeric_limits<double>::epsilon()
        * std::max({std::abs(a), std::abs(b), dt});
}

template<class FieldType, class LabelAt>
void checkLabels(const FieldType& field, int count, LabelAt labelAt)
{
    if (field.size() != count || !field.value)
        throw std::invalid_argument("Hydro step field does not match mesh");
    for (int i = 0; i < count; ++i)
        if (field.label(i) != labelAt(i))
            throw std::invalid_argument("Hydro step field has stale native labels");
}
}

HydroForceSource::HydroForceSource(HydroDomainKind kind,
    std::shared_ptr<const NativeMeshIdentity> meshIdentity, double sourceTime,
    HydroForceEvolution evolution,
    std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> interfaceGeometry,
    std::shared_ptr<const Quantum::Geometry::Parametric::PolynomialCurve<5>> parametricGeometry)
    : kind_(kind), meshIdentity_(std::move(meshIdentity)), sourceTime_(sourceTime),
      evolution_(evolution), interfaceGeometry_(std::move(interfaceGeometry)),parametricGeometry_(std::move(parametricGeometry))
{
    if (!meshIdentity_ || !std::isfinite(sourceTime_) || sourceTime_ < 0.0
        || (evolution_ != HydroForceEvolution::Frozen && evolution_ != HydroForceEvolution::TimeDependent)
        || (interfaceGeometry_ && parametricGeometry_))
        throw std::invalid_argument("invalid Hydro force source identity, time or evolution");
}

std::shared_ptr<const HydroForceSource> HydroForceSource::fitted(const HydroFittedDomain& domain,
    double sourceTime, HydroForceEvolution evolution,
    std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> interfaceGeometry)
{
    if (!domain.axisymmetric()) throw std::invalid_argument("fitted Hydro force source requires RZ");
    return std::shared_ptr<const HydroForceSource>(new HydroForceSource(HydroDomainKind::FittedRZ,
        domain.geometry().meshIdentity(), sourceTime, evolution, std::move(interfaceGeometry)));
}

void HydroForceSource::validate(HydroDomainKind kind,
    const std::shared_ptr<const NativeMeshIdentity>& expectedMesh,
    double stepBegin, double dt, const LevelSetHermiteGeometry::Snapshot* expectedInterface,
    const InterfaceState* expectedParametric) const
{
    if (kind_ != kind || !expectedMesh || meshIdentity_ != expectedMesh)
        throw std::invalid_argument("Hydro force source belongs to another mesh/domain incarnation");
    if (interfaceGeometry_.get() != expectedInterface)
        throw std::invalid_argument("Hydro force source does not match the current retained interface");
    if((expectedParametric && !expectedParametric->hasParametricGeometry()) ||
        parametricGeometry_.get()!=(expectedParametric?expectedParametric->parametricGeometry().geometryHandle().get():nullptr))
        throw std::invalid_argument("Hydro force source does not match the retained parameter geometry");
    if ((evolution_ == HydroForceEvolution::TimeDependent && !sameTime(sourceTime_, stepBegin, dt))
        || (sourceTime_ > stepBegin && !sameTime(sourceTime_, stepBegin, dt)))
        throw std::invalid_argument("Hydro force source is stale or from a future time; relabelling evaluationTime cannot refresh it");
    // The shared identity is checked against the actual mesh by advance()
    // (complete domain) or the existing fitted kernel, before any field update.
    // Provenance does not add a second whole-mesh scan per loaded step.
}

std::shared_ptr<const HydroForceSource> Hydrodynamics::captureParametricForceSource(const QuadtreeMesh& mesh,
    const NativeLeafDomain& domain,const InterfaceState& state) const
{
    InterfaceParametricMotion::validateDomain(mesh,domain,state);
    if(!initialStateSet_ || stepProgress_.stopped || stepProgress_.geometry!=HydroDomainKind::CompleteXY || !stepMeshIdentity_)
        throw std::logic_error("parameter force provenance requires an initialized complete XY Hydro state");
    stepMeshIdentity_->validate(mesh);
    return std::shared_ptr<const HydroForceSource>(new HydroForceSource(HydroDomainKind::CompleteXY,
        stepMeshIdentity_,state.progress().endTime,HydroForceEvolution::TimeDependent,{},state.parametricGeometry().geometryHandle()));
}

std::shared_ptr<const HydroForceSource> Hydrodynamics::captureForceSource(double sourceTime,
    HydroForceEvolution evolution,
    std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> interfaceGeometry) const
{
    if (!initialStateSet_ || stepProgress_.stopped || !stepProgress_.geometry)
        throw std::logic_error("Hydro force provenance requires an initialized, non-stopped managed state");
    const auto identity = fittedDomain_ ? fittedDomain_->geometry().meshIdentity() : stepMeshIdentity_;
    return std::shared_ptr<const HydroForceSource>(new HydroForceSource(*stepProgress_.geometry,
        identity, sourceTime, evolution, std::move(interfaceGeometry)));
}

HydroSolveReport Hydrodynamics::advance(QuadtreeMesh& mesh, const HydroStepInput& input)
{
    return advanceStep(mesh, input, false);
}

HydroSolveReport Hydrodynamics::advanceStep(QuadtreeMesh& mesh, const HydroStepInput& input, bool deferAcceptance)
{
    HydroSolveReport rejected;
    rejected.stage = "step-input";
    try
    {
        const auto& time = input.time;
        const auto& force = input.force;
        validateStepTransactionInput();
        const bool fitted = input.geometry.kind == HydroDomainKind::FittedRZ;
        if(time.policy!=HydroTimePolicy::Fixed && time.policy!=HydroTimePolicy::ConservativeOneStep)
            throw std::invalid_argument("invalid Hydro time policy");
        if(time.policy==HydroTimePolicy::ConservativeOneStep && (fitted || !input.conservativeTransport))
            throw std::invalid_argument("variable Hydro time policy requires complete conservative mass transport");
        if (stepProgress_.stopped)
        {
            rejected.requiresReinitialization = true;
            throw std::invalid_argument("Hydro stopped after a failed substep; rebuild/restore all state before retrying");
        }
        if (stepProgress_.pendingAcceptance)
            throw std::invalid_argument("Hydro has a pending substep; accept it before another step");
        if (input.geometry.kind != HydroDomainKind::CompleteXY
            && input.geometry.kind != HydroDomainKind::CompleteRZ && !fitted)
            throw std::invalid_argument("unsupported Hydro domain capability");
        if (fitted != (input.geometry.fitted != nullptr))
            throw std::invalid_argument("Hydro domain kind and fitted geometry disagree");
        if (fitted && input.options.completeConservation)
            throw std::invalid_argument("fitted Hydro does not consume completeConservation");
        if(input.options.compatibleNativePressure && (fitted || !input.conservativeTransport || !input.options.boundaryPatches.empty()))
            throw std::invalid_argument("compatible native pressure requires closed complete conservative mass transport");
        if(input.options.compatibleNativePressure && static_cast<long long>(mesh.Nx)*mesh.Ny!=mesh.cellIndexNum())
            throw std::invalid_argument("compatible native pressure requires the uniform root grid until AMR Hodge is implemented");
        if(stepProgress_.completedSteps>0 && input.options.compatibleNativePressure!=activeOptions_.compatibleNativePressure)
            throw std::invalid_argument("cannot change pressure discretization across accepted Hydro history");
        if(bool(input.midpoint)!=input.options.conservativeMidpoint || (input.midpoint &&
            (!input.options.compatibleNativePressure || fitted ||
             time.policy!=HydroTimePolicy::ConservativeOneStep || !input.midpoint->domain ||
             !input.midpoint->velocityX || !input.midpoint->velocityY)))
            throw std::invalid_argument("midpoint momentum requires its complete compatible physical iterate");
        if(stepProgress_.completedSteps>0 && input.options.conservativeMidpoint!=activeOptions_.conservativeMidpoint)
            throw std::invalid_argument("cannot change momentum time scheme across accepted Hydro history");
        if(input.midpoint){
            input.midpoint->domain->validateMeshIdentity(mesh);
            if(input.midpoint->domain->axisymmetric()!=(input.geometry.kind==HydroDomainKind::CompleteRZ))
                throw std::invalid_argument("midpoint domain has a different coordinate measure");
            if(input.midpoint->interface)InterfaceParametricMotion::validateDomain(mesh,*input.midpoint->domain,*input.midpoint->interface);
            else if(!input.conservativeTransport || !input.conservativeTransport->isSinglePhase() || !materialState_
                || materialState_->evolution!=HydroMaterialEvolution::Frozen || input.forceInterface
                || input.midpoint->viscosityNegative!=0 || input.midpoint->viscosityPositive!=0)
                throw std::invalid_argument("midpoint without an interface requires single-phase mass and frozen owned viscosity");
        }
        if (!std::isfinite(time.begin) || time.begin < 0.0 || !std::isfinite(time.dt) || time.dt <= 0.0
            || !std::isfinite(time.begin + time.dt) || time.begin + time.dt <= time.begin
            || time.index <= 0 || stepProgress_.completedSteps == std::numeric_limits<int>::max()
            || time.index != stepProgress_.completedSteps + 1
            || !sameTime(time.begin, stepProgress_.endTime, time.dt))
            throw std::invalid_argument("Hydro step index/start time does not continue its state history");
        if (stepProgress_.geometry && *stepProgress_.geometry != input.geometry.kind)
            throw std::invalid_argument("cannot change Hydro coordinate/domain capability across time history");
        if (stepEntry_ == StepEntry::Legacy
            || (stepProgress_.completedSteps == 0 && (fixedDt_ > 0.0 || pressureMode_ != PressureMode::Unset)))
            throw std::invalid_argument("legacy Hydro history has no managed step identity; explicit restoration is required");
        if(force.timeLevel!=HydroForceTimeLevel::StepBegin && force.timeLevel!=HydroForceTimeLevel::StepEnd && force.timeLevel!=HydroForceTimeLevel::StepMidpoint)
            throw std::invalid_argument("invalid Hydro force quadrature time level");
        const bool targetForce=force.timeLevel==HydroForceTimeLevel::StepEnd;
        const bool midpointForce=force.timeLevel==HydroForceTimeLevel::StepMidpoint;
        const double forceTime=time.begin+(midpointForce?.5:(targetForce?1:0))*time.dt;
        if(midpointForce!=bool(input.midpoint) || (input.midpoint && input.midpoint->interface!=input.forceParametricInterface))
            throw std::invalid_argument("midpoint momentum and force require the same half-step geometry");
        if(input.forceInterface && input.forceParametricInterface)
            throw std::invalid_argument("Hydro force cannot bind two interface representations");
        if((targetForce || (midpointForce && input.midpoint->interface)) && (fitted || !input.conservativeTransport || (!input.forceInterface && !input.forceParametricInterface) || !force.source
            || force.source->evolution()!=HydroForceEvolution::TimeDependent))
            throw std::invalid_argument("end-step Hydro force requires a prepared conservative interface and time-dependent source");
        if (!std::isfinite(force.evaluationTime) || !sameTime(force.evaluationTime, forceTime, time.dt))
            throw std::invalid_argument("Hydro force density does not belong to its declared quadrature time");
        if (!initialStateSet_)
            throw std::invalid_argument("Hydro advance requires explicit initializeState before its first step");
        if (input.conservativeTransport) {
            if(fitted || (stepProgress_.completedSteps>0 && !conservativeTransportActive_))
                throw std::invalid_argument("conservative momentum requires complete geometry and selection at the first step");
            input.conservativeTransport->validateFor(mesh,input.geometry.kind==HydroDomainKind::CompleteRZ,time,rho);
        } else if(conservativeTransportActive_)
            throw std::invalid_argument("conservative momentum history requires its phase transport at every step");
        if(input.forceParametricInterface) {
            if(fitted || !input.conservativeTransport)
                throw std::invalid_argument("parameter force requires its conservative mass step");
            input.forceParametricInterface->validateMesh(mesh);
            input.conservativeTransport->validateParametricForce(*input.forceParametricInterface,force.timeLevel);
        }
        if (!materialState_
            || (materialState_->evolution == HydroMaterialEvolution::TimeDependent
                && !sameTime(materialState_->sourceTime, time.begin, time.dt)))
            throw std::invalid_argument("Hydro materials are missing or stale for the step source time");
        if ((force.cellX == nullptr) != (force.cellY == nullptr)
            || (force.faceX == nullptr) != (force.faceY == nullptr))
            throw std::invalid_argument("Hydro force density requires both components at each location");
        if (!input.options.closedDomainProjection || input.options.surfaceForces
            || input.options.evaporation || input.options.directForce)
            throw std::invalid_argument("Hydro advance requires incremental pressure and explicit external loads");

        // Reuse the existing exact identity; no second fitted-domain snapshot.
        // Complete domains capture once and reject remeshing without transfer.
        auto identity = stepMeshIdentity_;
        if (!fitted)
        {
            if (identity) identity->validate(mesh);
            else identity = std::make_shared<const NativeMeshIdentity>(mesh);
        }
        if ((force.cellX || force.faceX) && !force.source)
            throw std::invalid_argument("Hydro native force fields require their assembly-time source record");
        if (force.source)
            force.source->validate(input.geometry.kind,
                fitted ? input.geometry.fitted->geometry().meshIdentity() : identity,
                forceTime, time.dt, input.forceInterface,input.forceParametricInterface);
        else if (input.forceInterface || input.forceParametricInterface)
            throw std::invalid_argument("Hydro interface load context requires a force source record");
        const auto cell = [&](int i) { return mesh.cellIndexToLabel[i]; };
        const auto vertex = [&](int i) { return mesh.vertIndexToLabel[i]; };
        const auto vx = [&](int i) { return mesh.vEdgeIndexToLabel[i]; };
        const auto vy = [&](int i) { return mesh.hEdgeIndexToLabel[i]; };
        // Fitted kernels already validate all state/force labels. Complete
        // kernels retain their finite/positive/alias checks below this port.
        if (!fitted)
        {
            for (const auto* f : {&rho, &mu, &p}) checkLabels(*f, mesh.cellIndexNum(), cell);
            for (const auto* f : {&vertRho, &vertMu}) checkLabels(*f, mesh.vertIndexNum(), vertex);
            for (const auto* f : {&ux, &uxPre, &uxMid}) checkLabels(*f, mesh.vEdgeIndexNum(), vx);
            for (const auto* f : {&uy, &uyPre, &uyMid}) checkLabels(*f, mesh.hEdgeIndexNum(), vy);
            if (force.cellX) { checkLabels(*force.cellX, mesh.cellIndexNum(), cell); checkLabels(*force.cellY, mesh.cellIndexNum(), cell); }
            if (force.faceX) { checkLabels(*force.faceX, mesh.vEdgeIndexNum(), vx); checkLabels(*force.faceY, mesh.hEdgeIndexNum(), vy); }
            if(input.midpoint){checkLabels(*input.midpoint->velocityX,mesh.vEdgeIndexNum(),vx);checkLabels(*input.midpoint->velocityY,mesh.hEdgeIndexNum(),vy);
                if(!input.midpoint->interface){std::vector<double> velocity(mesh.vEdgeIndexNum()+mesh.hEdgeIndexNum());
                    for(int row=0;row<static_cast<int>(velocity.size());++row)velocity[row]=row<mesh.vEdgeIndexNum()
                        ?input.midpoint->velocityX->value[row]:input.midpoint->velocityY->value[row-mesh.vEdgeIndexNum()];
                    input.conservativeTransport->validateSinglePhaseMidpoint(velocity);}}
        }
        auto report = fitted
            ? solveFittedStep(mesh, time.dt, input.options, *input.geometry.fitted,
                force.cellX, force.cellY, force.faceX, force.faceY)
            : solveStep<QuadtreeMesh>(mesh, time.dt, input.geometry.kind == HydroDomainKind::CompleteRZ,
                input.options, nullptr, nullptr, true, force.cellX, force.cellY, force.faceX, force.faceY, StepEntry::Managed,
                input.conservativeTransport,time.policy==HydroTimePolicy::ConservativeOneStep,input.midpoint);
        if (report.success)
        {
            conservativeTransportActive_=input.conservativeTransport!=nullptr;
            stepProgress_.completedSteps = time.index;
            stepProgress_.endTime = time.begin + time.dt;
            stepProgress_.geometry = input.geometry.kind;
            stepMeshIdentity_ = std::move(identity);
            stepProgress_.pendingAcceptance = deferAcceptance;
            if (!deferAcceptance)
            {
                stepProgress_.acceptedSteps = time.index;
                stepProgress_.acceptedEndTime = stepProgress_.endTime;
            }
        }
        else
        {
            stepProgress_.stopped = true;
            report.requiresReinitialization = true;
        }
        return report;
    }
    catch (const std::exception& error)
    {
        rejected.reason = error.what();
        return rejected;
    }
}

HydroStateView Hydrodynamics::state() const
{
    return {HydroReadOnlyField<CellField<double>>(p), HydroReadOnlyField<CellField<double>>(rho),
        HydroReadOnlyField<CellField<double>>(mu), HydroReadOnlyField<VEdgeField<double>>(ux),
        HydroReadOnlyField<VEdgeField<double>>(uxPre), HydroReadOnlyField<VEdgeField<double>>(uxMid),
        HydroReadOnlyField<HEdgeField<double>>(uy), HydroReadOnlyField<HEdgeField<double>>(uyPre),
        HydroReadOnlyField<HEdgeField<double>>(uyMid), HydroReadOnlyField<CellField<double>>(divU),
        HydroReadOnlyField<CellField<double>>(uMag), HydroReadOnlyField<VertField<double>>(vertRho),
        HydroReadOnlyField<VertField<double>>(vertMu), HydroReadOnlyField<CellField<double>>(externalFx),
        HydroReadOnlyField<CellField<double>>(externalFy), stepProgress_, initialStateSet_, materialState_};
}
