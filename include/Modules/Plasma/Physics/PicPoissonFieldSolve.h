#pragma once

#include "Solvers/Poisson/PoissonFVSolverCellCentered.h"
#include "Modules/Plasma/Physics/PicCurrentWeights.h"
#include <filesystem>

class PicFieldState;
class ImplicitPicFieldHistory;
class ExplicitPicChargeDeposition;
class PicFieldBoundaryState;
class NormalizationContext;
struct PicExplicitStepTiming;
template<class Mesh> class ImmersedGeometry;

// Owns the ordinary PIC Poisson solver, options and explicit setup bookkeeping.
// Bound fields/mesh outlive this nonmovable service. Geometry and boundary
// inputs are borrowed only during each solve/evaluation.
class PicPoissonFieldSolve
{
public:
    using Options = PoissonFVSolverCellCentered::Options;
    using Report = PoissonFVSolverCellCentered::SolveReport;
    struct Input {
        const FaceSegmentGraph& graph;
        const PicDomainBC& domain;
        const Solvers::Poisson::CutCellMetrics* cut = nullptr;
        const PicEmbeddedBC* embedded = nullptr;
        const PicEmbeddedBCMap* embeddedById = nullptr;
    };
    enum class BoundaryPolicy { Configured, LegacyDefault };

    PicPoissonFieldSolve(PicFieldState& fields, Options options);
    PicPoissonFieldSolve(const PicPoissonFieldSolve&) = delete;
    PicPoissonFieldSolve& operator=(const PicPoissonFieldSolve&) = delete;
    const Options& options() const { return options_; }
    bool hasCache() const { return solver_.hasCache(); }

    void solveExplicit(const Input& input, int step, Report& report);
    struct ExplicitStageInput {
        const std::vector<ParticleType>& species;
        const std::vector<double>& macroWeight;
        const NormalizationContext& normalization;
        const PicCurrentWeights::Policy& coupling;
        const PicEmbeddedGeometry& embedded;
        const ImmersedGeometry<QuadtreeMesh>* geometry;
        double timeSI;
        int step;
        bool collectElectronAmount;
        int electronIndex;
        std::vector<double>* normalFieldBySegment = nullptr;
    };
    // Current-position weights -> cut-volume charge/RHS -> driven boundaries
    // -> ordinary Poisson -> current E. All owners must share this field/mesh;
    // AMR resets remain explicit. No particle update, clock or output callback.
    // Boundaries keep the evaluated values for same-stage diagnostics. Failure
    // after preflight stops the stage without rolling back completed operations.
    void advanceExplicit(const FaceSegmentGraph& graph, const ExplicitStageInput& input,
        PicCurrentWeights& weights, ExplicitPicChargeDeposition& deposition,
        PicFieldBoundaryState& boundaries, Report& report, PicExplicitStepTiming* timing = nullptr);
    // Original implicit initialization/rebuild solve and default-policy E.
    void solveForRebuild(const Input& input);
    struct ImplicitRebuildInput {
        const std::vector<ParticleType>& species;
        const std::vector<double>& macroWeight;
        const NormalizationContext& normalization;
        const PicCurrentWeights::Policy& coupling;
        const PicEmbeddedGeometry& embedded;
        bool useCutMetrics;
    };
    // Rebuild from current particles using the original full-cell-volume
    // closure and default boundary policy. Optional history is seeded only
    // after successful E evaluation; initial startup leaves it uncreated.
    void rebuildImplicit(const Input& input, const ImplicitRebuildInput& particles,
        PicCurrentWeights& weights, ImplicitPicFieldHistory* history = nullptr);
    // AMR fast path: preserve transferred phi and deposited fields, evaluate E,
    // then seed next/half. The caller performs the original mesh/field resets.
    void reseedImplicitFromPotential(const Input& input, ImplicitPicFieldHistory& history);
    // True preserves all IPF1 layers exactly, without rebuilding weights or
    // solving. False runs the same particle rebuild after the legacy reader
    // rejects the file. The outer checkpoint completeness gate stays in IO.
    bool restoreImplicit(const std::filesystem::path& directory, const Input& input,
        const ImplicitRebuildInput& particles, PicCurrentWeights& weights,
        ImplicitPicFieldHistory& history);
    void evaluateCurrent(const Input& input, BoundaryPolicy policy,
        std::vector<double>* normalFieldBySegment = nullptr);
    // Writes the next layer only; half-step formation/acceptance stay separate.
    void evaluateNext(const Input& input, ImplicitPicFieldHistory& history,
        std::vector<double>* normalFieldBySegment = nullptr);

    // Geometry invalidation preserves the old cached-solve retry path. AMR
    // explicitly restarts setup and its combined timing/report bookkeeping.
    void clearCache();
    void resetAfterMeshChange();

private:
    friend class ImplicitPicRefinement;
    friend class ExplicitPicRefinement;
    friend class ImplicitPicFieldSolve;
    void requireImplicitRebuildState(const ImplicitPicFieldHistory* history) const;
    void evaluate(const Input& input, CellField<Vector2D>& target, BoundaryPolicy policy,
        std::vector<double>* normalFieldBySegment);
    PicFieldState& fields_;
    const Options options_;
    PoissonFVSolverCellCentered solver_;
    bool setupPerformed_ = false;
};
