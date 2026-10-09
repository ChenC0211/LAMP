#include "Modules/Electricity/DielectricEvolution.h"
#include "Core/Geometry/InterfaceGeometry.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace Quantum::Electricity
{
namespace {
std::shared_ptr<const DielectricChargeSnapshot> captureCharge(const QuadtreeMesh& mesh,
    const std::shared_ptr<SurfaceChargeField>& charge, int completedStep, double dt)
{
    auto result = std::make_shared<DielectricChargeSnapshot>(mesh);
    result->completedStep = completedStep;
    result->time = completedStep * dt;
    if (!std::isfinite(result->time)) throw std::runtime_error("charge publication time overflow");
    result->geometry = std::shared_ptr<const LevelSetInterfacePanels>(charge, &charge->geometry());
    result->integratedCharge = charge->transport().panelCharge();
    result->measure = charge->transport().panelMeasure();
    return result;
}
}
std::vector<double> DielectricChargeSnapshot::density() const
{
    std::vector<double> values(integratedCharge.size());
    for (std::size_t i = 0; i < values.size(); ++i) values[i] = integratedCharge[i] / measure[i];
    return values;
}

DielectricEvolution::DielectricEvolution(const QuadtreeMesh& mesh, DielectricConfig config, double thickness)
    : mesh_(&mesh), config_(std::move(config)), thickness_(thickness), electric_(mesh)
{
    if (config_.interfaceDiscretization != "sharpHermiteP3" && (!std::isfinite(thickness_) || !(thickness_ > 0.0)))
        throw std::invalid_argument("electric interface thickness must be finite and positive");
}

void DielectricEvolution::requireUsable() const
{
    if (failed_) throw std::logic_error("failed electric/charge state requires reinitialization");
}
std::shared_ptr<const DielectricEvaluation> DielectricEvolution::evaluated() const
{
    requireUsable();
    if (!publication_) throw std::logic_error("electric fields have not been solved");
    return publication_;
}
std::shared_ptr<const DielectricChargeSnapshot> DielectricEvolution::charge() const
{ requireUsable(); return chargePublication_; }

double DielectricEvolution::explicitRelaxationCfl(const DielectricConfig& config, double dt)
{
    if (!std::isfinite(dt) || !(dt > 0.0) ||
        !std::isfinite(config.permittivityNegative) || !(config.permittivityNegative > 0.0) ||
        !std::isfinite(config.permittivityPositive) || !(config.permittivityPositive > 0.0) ||
        !std::isfinite(config.conductivityNegative) || config.conductivityNegative < 0.0 ||
        !std::isfinite(config.conductivityPositive) || config.conductivityPositive < 0.0)
        throw std::invalid_argument("invalid explicit charge material or time input");
    const double rate = std::max(config.conductivityNegative / config.permittivityNegative,
        config.conductivityPositive / config.permittivityPositive);
    const double cfl = dt * rate;
    if (!std::isfinite(rate) || !std::isfinite(cfl) || cfl > .25)
        throw std::invalid_argument("explicit charge relaxation requires dt*max(conductivity/permittivity)<=0.25; reduce dt");
    return cfl;
}

std::shared_ptr<const DielectricEvaluation> DielectricEvolution::evaluate(const QuadtreeMesh& mesh, const VertField<double>& scalar,
    DielectricStepTime time, const InterfaceGeometry* geometry,
    DielectricSharpAssemblyDiagnostics* diagnostics, std::shared_ptr<const InterfaceTraceGeometry> traceQueries)
{
    requireUsable();
    if (pendingCurrent_) throw std::logic_error("previous charge step has not completed");
    if (&mesh != mesh_ || time.step < 1 || time.step != evaluatedTime_.step + 1 ||
        !std::isfinite(time.dt) || !(time.dt > 0.0) || !std::isfinite(time.forceTime) ||
        time.forceTime != (time.step-1)*time.dt || (evaluatedTime_.step && time.dt != evaluatedTime_.dt))
        throw std::invalid_argument("electric step requires the next fixed-dt layer on the original mesh");
    const bool charging = config_.model == "leakyDielectric";
    if (charging)
    {
        explicitRelaxationCfl(config_, time.dt);
        if (!geometry || config_.interfaceDiscretization != "sharpHermiteP3")
            throw std::invalid_argument("charging requires shared sharpHermiteP3 interface geometry");
    }
    if (geometry) geometry->validateSource(mesh, scalar);
    if (traceQueries && (config_.interfaceDiscretization != "sharpHermiteP3" || !geometry
        || traceQueries->source() != geometry->source()))
        throw std::invalid_argument("trace queries must use the same retained sharp scalar");

    failed_ = true; // Constructor/solve can mutate owned caches and fields.
    SurfaceChargeSource source;
    const DielectricConductivity conductivity{config_.conductivityNegative, config_.conductivityPositive};
    if (charging)
    {
        if (!charge_)
        {
            charge_ = std::make_shared<SurfaceChargeField>(mesh, geometry->snapshot(),
                config_.initialSurfaceChargeDensity, config_.axisymmetric);
            chargePublication_ = captureCharge(mesh, charge_, 0, time.dt);
        }
        if (!chargePublication_ || chargePublication_->completedStep != time.step - 1
            || chargePublication_->time != time.forceTime)
            throw std::logic_error("charge publication does not match the requested field time");
        source = charge_->source();
    }
    auto solution = electric_.solve(mesh, scalar, thickness_, config_, geometry, diagnostics,
        charging ? &source : nullptr, charging ? &conductivity : nullptr, std::move(traceQueries));
    const auto& report = solution.report;
    if (!report.success || !solution.fields) throw std::runtime_error("electric " + report.stage + ": " + report.reason);
    auto next = std::make_shared<DielectricEvaluation>(DielectricEvaluation{time, std::move(solution), chargePublication_});
    if (charging)
    {
        if (!next->electric.report.conductionCurrent)
            throw std::runtime_error("completed electric solve did not return this step's conduction current");
        pendingCurrent_ = next->electric.report.conductionCurrent;
    }
    publication_ = std::move(next);
    evaluatedTime_ = time;
    failed_ = false;
    return publication_;
}

DielectricChargeUpdate DielectricEvolution::completeCharge(DielectricStepTime time,
    const StaggeredVelocitySampler* velocity, const DielectricChargeMotion* motion)
{
    requireUsable();
    if (!pendingCurrent_ || time.step != evaluatedTime_.step || time.forceTime != evaluatedTime_.forceTime || time.dt != evaluatedTime_.dt)
        throw std::invalid_argument("charge completion does not match its staged electric current");
    if ((velocity && motion) || (config_.chargeConvection != bool(velocity || motion)) ||
        (motion && (&motion->mesh != mesh_ || !motion->forwardMaterialPoint)))
        throw std::invalid_argument("charge transport velocity/material-map input is inconsistent");

    failed_ = true;
    DielectricChargeUpdate result;
    result.sourceTime = time;
    result.current = std::move(pendingCurrent_); // One shot, including failure.
    const auto& current = *result.current;
    const double closureTolerance = config_.chargeClosureTolerance.value_or(16e-9);
    std::shared_ptr<SurfaceChargeField> nextCharge;
    if (motion)
        nextCharge = charge_->advanceMoving(motion->mesh, motion->snapshot, motion->forwardMaterialPoint,
            current.cellIntegratedCurrent, current.absoluteFaceCurrent,
            closureTolerance, time.dt, result.transport);
    else if (velocity)
        result.transport = charge_->advance(current.cellIntegratedCurrent, charge_->advectiveMeasureFlux(*velocity),
            current.absoluteFaceCurrent, closureTolerance, time.dt);
    else
        result.transport = charge_->advanceNormalCurrent(current.cellIntegratedCurrent,
            current.absoluteFaceCurrent, closureTolerance, time.dt);
    if (!result.transport.transport.success)
        throw std::runtime_error("surface charge " + result.transport.transport.stage + ": " + result.transport.transport.reason);
    const auto& q = result.transport.transport;
    result.electrodeResidual = q.totalChargeAfter - q.totalChargeBefore + time.dt*current.netBoundaryCurrent;
    const double scale = q.conservationScale + time.dt*current.absoluteFaceCurrent;
    if (!std::isfinite(result.electrodeResidual) || !std::isfinite(scale) ||
        std::abs(result.electrodeResidual) > closureTolerance*scale)
        throw std::runtime_error("surface charge increment does not close with electrode conduction current");
    result.charge = captureCharge(*mesh_, nextCharge ? nextCharge : charge_, time.step, time.dt);
    if (nextCharge) charge_ = std::move(nextCharge);
    chargePublication_ = result.charge;
    failed_ = false;
    return result;
}
}
