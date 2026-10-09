#include "Modules/Plasma/Physics/PicPoissonFieldSolve.h"
#include "Modules/Plasma/Physics/PicFieldState.h"
#include "Modules/Plasma/ExplicitPIC/Field/ElectricFieldEvaluator.h"
#include "Modules/Plasma/Physics/PicChargeDeposition.h"
#include "Modules/Plasma/Boundary/PicFieldBoundaryState.h"
#include "Modules/Plasma/Boundary/PicEmbeddedGeometry.h"
#include "Modules/Plasma/Diagnostics/PicStepDiagnostics.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <chrono>
#include <iostream>
#include <stdexcept>
#include <utility>

PicPoissonFieldSolve::PicPoissonFieldSolve(PicFieldState& fields, Options options)
    : fields_(fields), options_(std::move(options))
{
    if (!options_.cacheEnabled) solver_.clearCache();
}

void PicPoissonFieldSolve::clearCache() { solver_.clearCache(); }
void PicPoissonFieldSolve::resetAfterMeshChange()
{
    setupPerformed_ = false;
    solver_.clearCache();
}

void PicPoissonFieldSolve::solveExplicit(const Input& input, int stepForPrint, Report& rep)
{
    const auto& mesh = fields_.mesh_;
    const auto& faceGraph = input.graph;
    const auto& domainBC = input.domain;
    const auto* cutMetricsPtr = input.cut;
    const auto* embeddedBCPtr = input.embedded;
    const auto* embeddedBCByIdPtr = input.embeddedById;
    const auto& rhs = fields_.rhs_;
    auto& phi = fields_.phi_;
    rep = {};
    Report setupRep, solveRep;
    if (!options_.cacheEnabled)
    {
        // No-cache mode: call the wrapper solve(), which will rebuild each call.
        const bool ok = solver_.solve(mesh, faceGraph, domainBC, rhs, phi, &solveRep, options_,
                                      cutMetricsPtr, embeddedBCPtr, embeddedBCByIdPtr);
        if (!ok)
        {
            std::cerr << "[ExplicitPIC] Poisson solve failed at step " << stepForPrint
                      << ": " << solveRep.message << std::endl;
            throw std::runtime_error("Poisson solve failed");
        }
        rep = solveRep;
    }
    else
    {
        const bool needSetup = !setupPerformed_;
        if (needSetup)
        {
            setupRep = solver_.setup(mesh, faceGraph, domainBC, options_, cutMetricsPtr,
                                     embeddedBCPtr, embeddedBCByIdPtr);
            if (!setupRep.success)
            {
                std::cerr << "[ExplicitPIC] Poisson setup failed at step " << stepForPrint
                          << ": " << setupRep.message << std::endl;
                throw std::runtime_error("Poisson setup failed");
            }
            setupPerformed_ = true;
        }

        bool ok = solver_.solveWithCache(mesh, faceGraph, domainBC, rhs, phi, &solveRep, options_,
                                         cutMetricsPtr, embeddedBCPtr, embeddedBCByIdPtr);
        if (!ok)
        {
            // Safety fallback: wrapper solve() will refresh the cache if needed.
            solveRep = {};
            ok = solver_.solve(mesh, faceGraph, domainBC, rhs, phi, &solveRep, options_,
                               cutMetricsPtr, embeddedBCPtr, embeddedBCByIdPtr);
        }
        if (!ok)
        {
            std::cerr << "[ExplicitPIC] Poisson solve failed at step " << stepForPrint
                      << ": " << solveRep.message << std::endl;
            throw std::runtime_error("Poisson solve failed");
        }

        rep = solveRep;
        if (needSetup)
        {
            // Combine setup (A assemble + factor) and cached solve (RHS + solve).
            rep.timeAssemble_s = setupRep.timeAssemble_s + solveRep.timeAssemble_s;
            rep.timeFactor_s = setupRep.timeFactor_s;
            rep.timeSolve_s = solveRep.timeSolve_s;
            rep.timeTotal_s = setupRep.timeTotal_s + solveRep.timeTotal_s;
            rep.N = setupRep.N;
            rep.nnz = setupRep.nnz;
            rep.usedSuperLU = setupRep.usedSuperLU;
            if (!setupRep.message.empty())
            {
                if (!rep.message.empty())
                    rep.message = setupRep.message + " | " + rep.message;
                else
                    rep.message = setupRep.message;
            }
        }
    }

}

void PicPoissonFieldSolve::advanceExplicit(const FaceSegmentGraph& graph,
    const ExplicitStageInput& input, PicCurrentWeights& weights,
    ExplicitPicChargeDeposition& deposition, PicFieldBoundaryState& boundaries,
    Report& report, PicExplicitStepTiming* timing)
{
    if (&deposition.fields_ != &fields_ || &weights.mesh_ != &fields_.mesh_
        || &input.coupling.mesh_ != &fields_.mesh_)
        throw std::invalid_argument("Explicit PIC field stage owners belong to different fields/meshes");
    if (input.coupling.mode_ != PicCurrentWeights::Policy::Mode::Explicit)
        throw std::invalid_argument("Explicit PIC field stage requires explicit coupling policy");
    const int cells = fields_.mesh_.cellIndexNum();
    if (fields_.phi_.size() != cells || fields_.rhs_.size() != cells
        || fields_.q_.size() != cells || fields_.rho_.size() != cells
        || fields_.electric_.size() != cells || deposition.amrNeInstant.size() != cells)
        throw std::logic_error("Explicit PIC field stage storage must be reset after mesh adaptation");

    using Clock = std::chrono::steady_clock;
    const auto secondsBetween = [](const Clock::time_point& a, const Clock::time_point& b) {
        return std::chrono::duration<double>(b - a).count();
    };
    if (timing) *timing = PicExplicitStepTiming{};
    const auto tBuild0 = Clock::now();
    weights.rebuild(input.species, input.coupling, input.embedded);
    const auto& caches = weights.caches();
    const auto tBuild1 = Clock::now();
    if (timing) timing->buildWeights_s = secondsBetween(tBuild0, tBuild1);

    const auto tDeposit0 = Clock::now();
    deposition.deposit(input.species, input.macroWeight, caches, input.normalization,
        input.embedded, input.geometry, input.collectElectronAmount, input.electronIndex);
    if (timing) timing->deposit_s = secondsBetween(tDeposit0, Clock::now());

    const auto tRhs0 = Clock::now();
    fields_.copyChargeDensityToRhs();
    if (timing) timing->rhs_s = secondsBetween(tRhs0, Clock::now());

    boundaries.update(input.timeSI);
    const auto* cut = input.embedded.metricsReady() ? &input.embedded.metrics() : nullptr;
    const Input fieldInput{graph, boundaries.domain(), cut,
        (cut && boundaries.embedded().enabled) ? &boundaries.embedded() : nullptr,
        (cut && !boundaries.embeddedById().empty()) ? &boundaries.embeddedById() : nullptr};
    solveExplicit(fieldInput, input.step, report);

    const auto tE0 = Clock::now();
    evaluateCurrent(fieldInput, BoundaryPolicy::Configured, input.normalFieldBySegment);
    if (timing) timing->evalE_s = secondsBetween(tE0, Clock::now());
}

void PicPoissonFieldSolve::solveForRebuild(const Input& input)
{
    Report rep;
    const bool ok = solver_.solve(fields_.mesh_, input.graph, input.domain, fields_.rhs_, fields_.phi_,
        &rep, options_, input.cut, input.embedded, input.embeddedById);
    if (!ok) throw std::runtime_error(std::string("[ImplicitPIC] Poisson solve failed: ") + rep.message);
    evaluateCurrent(input, BoundaryPolicy::LegacyDefault);
}

void PicPoissonFieldSolve::requireImplicitRebuildState(const ImplicitPicFieldHistory* history) const
{
    if (history && &history->fields_ != &fields_)
        throw std::invalid_argument("PIC electric history belongs to a different field state");
    const int count = fields_.mesh_.cellIndexNum();
    if (fields_.phi_.size() != count || fields_.rhs_.size() != count
        || fields_.rho_.size() != count || fields_.q_.size() != count || fields_.electric_.size() != count)
        throw std::logic_error("Implicit PIC rebuild fields must be reset after mesh adaptation");
    if (history) history->requireCurrentSizes();
}

void PicPoissonFieldSolve::rebuildImplicit(const Input& input, const ImplicitRebuildInput& particles,
    PicCurrentWeights& weights, ImplicitPicFieldHistory* history)
{
    if (&weights.mesh_ != &fields_.mesh_ || &particles.coupling.mesh_ != &fields_.mesh_)
        throw std::invalid_argument("Implicit PIC rebuild weights/policy belong to a different mesh");
    if (particles.coupling.mode_ != PicCurrentWeights::Policy::Mode::Implicit)
        throw std::invalid_argument("Implicit PIC rebuild requires implicit coupling policy");
    requireImplicitRebuildState(history);
    weights.rebuild(particles.species, particles.coupling, particles.embedded, particles.useCutMetrics);
    fields_.depositFullCellCharge(particles.species, particles.macroWeight, weights.caches(), particles.normalization);
    fields_.copyChargeDensityToRhs();
    solveForRebuild(input);
    if (history) history->seedFromCurrent();
}

void PicPoissonFieldSolve::reseedImplicitFromPotential(const Input& input, ImplicitPicFieldHistory& history)
{
    requireImplicitRebuildState(&history);
    evaluateCurrent(input, BoundaryPolicy::LegacyDefault);
    history.seedFromCurrent();
}

bool PicPoissonFieldSolve::restoreImplicit(const std::filesystem::path& directory, const Input& input,
    const ImplicitRebuildInput& particles, PicCurrentWeights& weights, ImplicitPicFieldHistory& history)
{
    if (&history.fields_ != &fields_)
        throw std::invalid_argument("PIC electric history belongs to a different field state");
    if (history.readCheckpoint(directory)) return true;
    rebuildImplicit(input, particles, weights, &history);
    return false;
}

void PicPoissonFieldSolve::evaluate(const Input& input, CellField<Vector2D>& target,
    BoundaryPolicy policy, std::vector<double>* normalFieldBySegment)
{
    ExplicitPIC::Field::ElectricFieldEvaluator::Options eOpt;
    if (policy == BoundaryPolicy::Configured && options_.embeddedBCMissingPolicy == Options::EmbeddedBCMissingPolicy::Skip)
        eOpt.embeddedBCMissingPolicy = ExplicitPIC::Field::ElectricFieldEvaluator::Options::EmbeddedBCMissingPolicy::Skip;
    ExplicitPIC::Field::ElectricFieldEvaluator::evaluateCellCenteredField(
        fields_.mesh_, input.graph, input.cut, input.domain, input.embedded, input.embeddedById,
        fields_.phi_, target, normalFieldBySegment, nullptr, eOpt);
}

void PicPoissonFieldSolve::evaluateCurrent(const Input& input, BoundaryPolicy policy,
    std::vector<double>* normalFieldBySegment)
{
    evaluate(input, fields_.electric_, policy, normalFieldBySegment);
}

void PicPoissonFieldSolve::evaluateNext(const Input& input, ImplicitPicFieldHistory& history,
    std::vector<double>* normalFieldBySegment)
{
    if (&history.fields_ != &fields_)
        throw std::invalid_argument("PIC electric history belongs to a different field state");
    evaluate(input, history.next_, BoundaryPolicy::Configured, normalFieldBySegment);
}
