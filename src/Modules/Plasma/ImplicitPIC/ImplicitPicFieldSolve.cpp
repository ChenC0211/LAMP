#include "Modules/Plasma/ImplicitPIC/ImplicitPicFieldSolve.h"
#include "Modules/Plasma/ImplicitPIC/ImplicitPicPolarization.h"
#include "Modules/Plasma/Physics/PicFieldState.h"
#include "Modules/Plasma/Diagnostics/PicStepDiagnostics.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/PeriodicFaceUtils.h"
#include "Solvers/Poisson/CutCellMetrics.h"
#include <algorithm>
#include <cctype>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <unordered_set>

using Solvers::Poisson::CutCellMetrics;
namespace {
std::string toLowerCopy(std::string s)
{
    for (auto& c : s) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return s;
}
static void enforcePeriodicCellSymmetry(const QuadtreeMesh& mesh,
                                        const FaceSegmentGraph& faceGraph,
                                        const CutCellMetrics* cutMetricsPtr,
                                        const bool periodicX,
                                        const bool periodicY,
                                        CellField<double>& phi)
{
    if (!(periodicX || periodicY))
        return;

    std::unordered_set<std::uint64_t> seenPairs;
    seenPairs.reserve(faceGraph.segments().size());
    auto packPair = [](const int a, const int b) -> std::uint64_t {
        const std::uint32_t lo = static_cast<std::uint32_t>(std::min(a, b));
        const std::uint32_t hi = static_cast<std::uint32_t>(std::max(a, b));
        return (static_cast<std::uint64_t>(lo) << 32) | static_cast<std::uint64_t>(hi);
    };

    const bool hasCutMask =
        (cutMetricsPtr != nullptr) &&
        (cutMetricsPtr->betaBySegment().size() == faceGraph.segments().size()) &&
        (cutMetricsPtr->cellMetricsByIndex().size() == static_cast<std::size_t>(mesh.cellIndexNum()));
    constexpr double periodicOpenTol = 1e-6;

    for (std::size_t sid = 0; sid < faceGraph.segments().size(); ++sid)
    {
        const auto& seg = faceGraph.segments()[sid];
        if (seg.boundary == Boundary::non)
            continue;
        if (!((periodicX && (seg.boundary == Boundary::west || seg.boundary == Boundary::east)) ||
              (periodicY && (seg.boundary == Boundary::south || seg.boundary == Boundary::north))))
        {
            continue;
        }

        const int kLab = seg.K;
        if (kLab < 0)
            continue;
        const int kIdx = mesh.cells(kLab).index();
        if (kIdx < 0)
            continue;

        const Vector2D shift = PeriodicFaceUtils::shiftToMate(mesh, seg.boundary);
        const int lLab = PeriodicFaceUtils::findMateCellLabel(mesh, seg, shift);
        if (lLab < 0)
            continue;
        const int lIdx = mesh.cells(lLab).index();
        if (lIdx < 0)
            continue;
        if (!PeriodicFaceUtils::isCleanCellMate(mesh, seg, lLab, shift))
            continue;

        if (hasCutMask)
        {
            double beta = cutMetricsPtr->betaBySegment()[sid];
            if (!std::isfinite(beta))
                beta = 1.0;
            beta = std::clamp(beta, 0.0, 1.0);
            if (!(beta > 0.0))
                continue;

            const auto& cm = cutMetricsPtr->cellMetricsByIndex();
            const auto& cmK = cm[static_cast<std::size_t>(kIdx)];
            const auto& cmL = cm[static_cast<std::size_t>(lIdx)];
            double alphaK = cmK.alpha;
            double alphaL = cmL.alpha;
            if (!std::isfinite(alphaK))
                alphaK = 1.0;
            if (!std::isfinite(alphaL))
                alphaL = 1.0;
            alphaK = std::clamp(alphaK, 0.0, 1.0);
            alphaL = std::clamp(alphaL, 0.0, 1.0);
            if (!(alphaK > 0.0) || !(alphaL > 0.0))
                continue;
            if ((beta < (1.0 - periodicOpenTol)) || cmK.isCut || cmL.isCut)
                continue;
        }

        const std::uint64_t pairKey = packPair(kIdx, lIdx);
        if (!seenPairs.insert(pairKey).second)
            continue;

        const double avg = 0.5 * (phi.value[kIdx] + phi.value[lIdx]);
        phi.value[kIdx] = avg;
        phi.value[lIdx] = avg;
    }
}
} // namespace

ImplicitPicFieldSolve::ImplicitPicFieldSolve(PicFieldState& fields,
    const PicPoissonFieldSolve::Options& poissonOptions,
    const PicConfig::PolarizationPoissonConfig& config, bool scheduledMotion, bool consoleInfo,
    const ImplicitPicOverrides::Polarization& overrides)
    : fields_(fields), consoleInfo_(consoleInfo)
{
    const std::string polModeCfgRaw = toLowerCopy(config.mode);
    const std::string polModeRaw = toLowerCopy(overrides.mode.value_or(polModeCfgRaw));
    mode_ =
        (polModeRaw == "direct" || polModeRaw == "cg" || polModeRaw == "ic")
            ? polModeRaw
            : std::string("auto");
    options_.useSuperLU = poissonOptions.useSuperLU;
    if (mode_ == "direct")
    {
        options_.useCG = false;
        options_.useIC = false;
    }
    else if (mode_ == "cg")
    {
        options_.useCG = true;
        options_.useIC = false;
    }
    else if (mode_ == "ic")
    {
        options_.useCG = false;
        options_.useIC = true;
    }
    else
    {
        options_.useCG = poissonOptions.useCG;
        options_.useIC = poissonOptions.useIC;
        if (!options_.useCG && !options_.useIC)
        {
            options_.useCG = false;
            options_.useIC = true;
        }
    }
    options_.maxIter = poissonOptions.maxIter;
    options_.relTol = poissonOptions.relTol;
    options_.verbose = poissonOptions.verbose;
    // Prescribed moving geometry must solve the current coefficients. A matrix
    // frozen by sparsity alone is neither a geometry identity nor restart state.
    options_.freezeDirectMatrixUntilPatternChange =
        config.freezeDirectMatrixUntilPatternChange && !scheduledMotion;
    options_.reuseIterativePattern = config.reuseIterativePattern;
    options_.useTemporalExtrapolatedGuess = overrides.temporalGuess.value_or(config.temporalGuess);
    options_.temporalExtrapolationOmega = overrides.temporalOmega.value_or(std::max(0.0, config.temporalOmega));
    options_.reuseIterativePattern = overrides.reuseIterativePattern.value_or(options_.reuseIterativePattern);
    if (consoleInfo_)
    {
        std::cout << "[ImplicitPIC][polPoisson] mode=" << mode_
                  << " direct=" << (options_.useSuperLU ? "superlu/sparselu" : "sparselu")
                  << " iterative=" << ((options_.useIC || options_.useCG)
                                           ? (options_.useIC ? "ic-cg" : "cg")
                                           : "off")
                  << " freezeDirect=" << (options_.freezeDirectMatrixUntilPatternChange ? 1 : 0)
                  << " temporalGuess=" << (options_.useTemporalExtrapolatedGuess ? 1 : 0)
                  << " omega=" << options_.temporalExtrapolationOmega
                  << std::endl;
    }
    options_.embeddedBCMissingPolicy = (poissonOptions.embeddedBCMissingPolicy ==
                                      PoissonFVSolverCellCentered::Options::EmbeddedBCMissingPolicy::Skip)
                                         ? PolarizationPoissonFVSolverCellCentered::Options::EmbeddedBCMissingPolicy::Skip
                                         : PolarizationPoissonFVSolverCellCentered::Options::EmbeddedBCMissingPolicy::Error;
}

ImplicitPicFieldSolve::Result ImplicitPicFieldSolve::solve(
    const Input& input, const ImplicitPicPolarization& polarization, int step)
{
    using Clock = std::chrono::steady_clock;
    const auto& mesh = fields_.mesh_;
    const auto& faceGraph = input.graph;
    const auto& domainBCUsed = input.domain;
    const auto* cutMetricsPtr = input.cut;
    const auto* embeddedBCPtr = input.embedded;
    const auto* embeddedBCByIdPtr = input.embeddedById;
    const auto& rhoStar = polarization.chargeDensity();
    const auto& kappaBySeg = polarization.faceCoefficients();
    const auto& kappaByCell = polarization.cellCoefficients();
    auto& phi = fields_.phi_;
    Result result;
    auto& polRep = result.report;
    auto polSolveOpt = options_;
    if ((mode_ == "auto") && iterativeSuppressed_)
    {
        polSolveOpt.useCG = false;
        polSolveOpt.useIC = false;
    }
    const auto tPolPoisson0 = Clock::now();
    const bool okPol = solver_.solve(mesh, faceGraph, domainBCUsed, rhoStar, phi,
                                        &kappaBySeg, &kappaByCell, &polRep, polSolveOpt,
                                        /*cutMetrics=*/cutMetricsPtr,
                                        /*embeddedBC=*/embeddedBCPtr,
                                        /*embeddedBCById=*/embeddedBCByIdPtr);
    result.solveSeconds = std::chrono::duration<double>(Clock::now() - tPolPoisson0).count();
    if (!okPol)
        throw std::runtime_error(std::string("[ImplicitPIC] polarization Poisson failed: ") + polRep.message);
    if ((mode_ == "auto"))
    {
        const bool attemptedIterative = polSolveOpt.useCG || polSolveOpt.useIC;
        if (attemptedIterative && polRep.usedCG)
        {
            iterativeSuppressed_ = false;
        }
        else if (attemptedIterative && (polRep.cgFallbackToDirect || !polRep.usedCG))
        {
            if (!iterativeSuppressed_ && consoleInfo_)
            {
                std::cout << "[ImplicitPIC][polPoisson] suppress iterative path after fallback at step="
                          << step << "; retry after next AMR" << std::endl;
            }
            iterativeSuppressed_ = true;
        }
    }
    if (consoleInfo_ && polRep.cgFallbackToDirect)
    {
        std::cout << "[ImplicitPIC][Event][POL] CG fallback to direct at step=" << step
                  << " iter=" << polRep.cgIterations
                  << " residual=" << std::scientific << std::setprecision(3) << polRep.cgResidual
                  << std::defaultfloat << std::setprecision(6)
                  << std::endl;
    }

    return result;
}

void ImplicitPicFieldSolve::enforcePeriodicSymmetry(const Input& input, bool periodicX, bool periodicY)
{
    enforcePeriodicCellSymmetry(fields_.mesh_, input.graph, input.cut, periodicX, periodicY, fields_.phi_);
}

void ImplicitPicFieldSolve::clearDirectCache() { solver_.clearDirectCache(); }
void ImplicitPicFieldSolve::resetAfterMeshChange(int step)
{
    // Keep direct polarization solver cache consistent with AMR topology updates.
    solver_.clearDirectCache();
    if ((mode_ == "auto") && iterativeSuppressed_)
    {
        iterativeSuppressed_ = false;
        if (consoleInfo_)
        {
            std::cout << "[ImplicitPIC][polPoisson] re-enable iterative path after AMR at step="
                      << step << std::endl;
        }
    }
}

ImplicitPicFieldSolve::Result ImplicitPicFieldSolve::advance(const Input& input, const StageInput& stage,
    ImplicitPicPolarization& polarization, PicPoissonFieldSolve& poisson,
    ImplicitPicFieldHistory& fieldHistory, PicStepParticleBalance& particleBalance,
    PicImplicitStepTiming& stepTiming)
{
    // Reject mixed owners and stale mesh-sized storage before any deposition.
    if (&poisson.fields_ != &fields_ || &fieldHistory.fields_ != &fields_
        || &polarization.mesh != &fields_.mesh_)
        throw std::invalid_argument("Implicit PIC field stage requires one field state and mesh");
    polarization.requireCurrentSizes();
    fieldHistory.requireCurrentSizes();

    using Clock = std::chrono::steady_clock;
    const auto secondsBetween = [](const Clock::time_point& a, const Clock::time_point& b) -> double {
        return std::chrono::duration_cast<std::chrono::duration<double>>(b - a).count();
    };
    const auto& species = stage.species;
    const auto& macroWeight = stage.macroWeight;
    const auto& cachesBar = stage.averagedWeights;
    const auto& norm = stage.normalization;
    const double dt_si = stage.dtSI;
    const int step = stage.step;
    const int depositFillOmpMinParticles = stage.depositFillOmpMinParticles;
    const bool periodicX = stage.periodicX, periodicY = stage.periodicY, consoleDebug = stage.consoleDebug;
    const auto& faceGraph = input.graph;
    const auto& domainBCUsed = input.domain;
    const auto* cutMetricsPtr = input.cut;
    const auto* embeddedBCPtr = input.embedded;
    const auto* embeddedBCByIdPtr = input.embeddedById;
    auto& polPoisson = *this;

    auto hasAnyDirichletBC = [&]() -> bool {
        auto isDirichlet = [](const auto& bc) -> bool {
            return bc.type == PicBCType::Dirichlet;
        };
        return isDirichlet(domainBCUsed.left) ||
               isDirichlet(domainBCUsed.right) ||
               isDirichlet(domainBCUsed.bottom) ||
               isDirichlet(domainBCUsed.top);
    };

    // (5) Deposit Q*/N* using orbit-avg cache.
    const auto tDepositOrbit0 = Clock::now();
    const auto deposited = polarization.deposit(species, macroWeight, cachesBar,
        depositFillOmpMinParticles, particleBalance);
    stepTiming.depositOrbit_qp_size_total += deposited.particles;
    stepTiming.depositOrbit_nnz_total += deposited.entries;
    stepTiming.depositOrbit_s += secondsBetween(tDepositOrbit0, Clock::now());

    // Convert Q* -> rho* (code), and assemble mu/kappa (SI).
    const auto tAssemblePol0 = Clock::now();
    polarization.assembleCells(species, norm, dt_si);

    // Neutrality correction: for pure Neumann / fully periodic closures the elliptic operator
    // has a constant nullspace and requires net charge to be (close to) zero. Apply a
    // minimal "remove mean rho" correction to satisfy compatibility.
    if (!hasAnyDirichletBC())
    {
        const auto correction = polarization.removeMeanCharge(norm);
        const double totalQ = correction.totalQ;
        const double meanRho_SI = correction.meanRhoSI;
        if (correction.applied && consoleDebug && step == 0)
        {
            std::cout << "[ImplicitPIC] neutrality correction applied (no Dirichlet BC): "
                      << "totalQ[C]=" << std::setprecision(6) << std::scientific << totalQ
                      << " meanRho[C/m^3]=" << meanRho_SI
                      << std::defaultfloat << std::setprecision(6)
                      << std::endl;
        }
    }

    // Stage-matched charge consistency at deposition time level:
    // compare Qgrid from Qstar and Qmacro built from the same particle stage.
    particleBalance.finishStage(polarization.totalCharge());

    // (6) kappa by face-segment (harmonic mean).
    polarization.assembleFaces(faceGraph);
    stepTiming.assemblePolarization_s += secondsBetween(tAssemblePol0, Clock::now());

    // (7) Solve polarization Poisson for phi^{n+1}.
    const auto polResult = polPoisson.solve(
        {faceGraph, domainBCUsed, cutMetricsPtr, embeddedBCPtr, embeddedBCByIdPtr}, polarization, step + 1);
    const auto& polRep = polResult.report;
    stepTiming.polPoisson_total_s += polResult.solveSeconds;
    stepTiming.polPoisson_assemble_s += polRep.timeAssemble_s;
    stepTiming.polPoisson_factor_s += polRep.timeFactor_s;
    stepTiming.polPoisson_solve_s += polRep.timeSolve_s;
    stepTiming.pol_cg_iterations = polRep.cgIterations;
    stepTiming.pol_cg_residual = polRep.cgResidual;

    // Project the solved potential onto the periodic-equality subspace on
    // clean periodic mates. This removes small north/south or west/east
    // asymmetries from the linear solve before E reconstruction/output.
    polPoisson.enforcePeriodicSymmetry(
        {faceGraph, domainBCUsed, cutMetricsPtr, embeddedBCPtr, embeddedBCByIdPtr}, periodicX, periodicY);

    // (8) E^{n+1} and E^{n+1/2}.
    const auto tEvalE0 = Clock::now();
    poisson.evaluateNext({faceGraph, domainBCUsed, cutMetricsPtr, embeddedBCPtr, embeddedBCByIdPtr},
        fieldHistory, stage.normalFieldBySegment);

    fieldHistory.formHalfStep();
    stepTiming.evalE_s += secondsBetween(tEvalE0, Clock::now());
    return polResult;
}
