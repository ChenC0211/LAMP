#pragma once
#include "Modules/Plasma/ExplicitPIC/PicConfig.h"
#include "Modules/Plasma/ImplicitPIC/ImplicitPicOverrides.h"
#include "Modules/Plasma/Physics/PicPoissonFieldSolve.h"
#include "Solvers/Poisson/PolarizationPoissonFVSolverCellCentered.h"

class ImplicitPicPolarization;
class ParticleType;
class CouplingWeightsCache;
class NormalizationContext;
class PicStepParticleBalance;
struct PicImplicitStepTiming;

// Owns the implicit polarization solver, prepared options and auto-mode retry
// state. Bound fields/mesh outlive this service; stage inputs are borrowed.
class ImplicitPicFieldSolve
{
public:
    using Input = PicPoissonFieldSolve::Input;
    using Options = PolarizationPoissonFVSolverCellCentered::Options;
    using Report = PolarizationPoissonFVSolverCellCentered::SolveReport;
    struct Result { Report report; double solveSeconds = 0.; };

    // One field substep, sampled after pre-collision and orbit construction.
    // step is the caller's zero-based macro-step; no clock or acceptance moves.
    // All references, including the optional normal-field destination, are
    // borrowed only for advance(). Particle data and caches remain unchanged.
    struct StageInput {
        const std::vector<ParticleType>& species;
        const std::vector<double>& macroWeight;
        const std::vector<CouplingWeightsCache>& averagedWeights;
        const NormalizationContext& normalization;
        double dtSI;
        int step, depositFillOmpMinParticles;
        bool periodicX, periodicY, consoleDebug;
        std::vector<double>* normalFieldBySegment = nullptr;
    };

    ImplicitPicFieldSolve(PicFieldState& fields, const PicPoissonFieldSolve::Options& poissonOptions,
        const PicConfig::PolarizationPoissonConfig& config, bool scheduledMotion, bool consoleInfo,
        const ImplicitPicOverrides::Polarization& overrides = {});
    ImplicitPicFieldSolve(const ImplicitPicFieldSolve&) = delete;
    ImplicitPicFieldSolve& operator=(const ImplicitPicFieldSolve&) = delete;
    const Options& options() const { return options_; }
    const std::string& mode() const { return mode_; }
    bool iterativeSuppressed() const { return iterativeSuppressed_; }

    Result solve(const Input& input, const ImplicitPicPolarization& polarization, int step);
    // Deposit/assemble, original domain-only neutrality correction, stage
    // charge accounting, solve/projection, next E and half E, in that order.
    // Services must share this solver's field state and mesh. The current E
    // stays unaccepted. Numerical failure propagates; this is not a rollback.
    Result advance(const Input& input, const StageInput& stage,
        ImplicitPicPolarization& polarization, PicPoissonFieldSolve& poisson,
        ImplicitPicFieldHistory& history, PicStepParticleBalance& balance,
        PicImplicitStepTiming& timing);
    // Kept after solve/report accounting and before electric reconstruction.
    void enforcePeriodicSymmetry(const Input& input, bool periodicX, bool periodicY);
    // Geometry invalidation retains retry suppression; AMR permits a retry.
    void clearDirectCache();
    void resetAfterMeshChange(int step);

private:
    PicFieldState& fields_;
    PolarizationPoissonFVSolverCellCentered solver_;
    Options options_;
    std::string mode_;
    const bool consoleInfo_;
    bool iterativeSuppressed_ = false;
};
