#pragma once

#include "Core/Field/CellField.h"
#include "Core/Base/Vector2D.h"
#include "Modules/Plasma/Physics/PicRefinementIndicators.h"
#include "Modules/Plasma/Physics/PicCurrentWeights.h"
#include "Modules/Plasma/Physics/PicPoissonFieldSolve.h"
#include <functional>

class QuadtreeMesh;
class QuadtreeMeshRefiner;
class UnitsContext;
class NormalizationContext;
class PicDomain;
class PicAmrGeometry;
class PicMovingGeometry;
class PicFieldState;
class ExplicitPicChargeDeposition;
class PicSpeciesState;
class PicPoissonFieldSolve;
class ImplicitPicPolarization;
class ImplicitPicOrbit;
struct AmrConfig;
template<class T> class VertField;

// The mesh outlives these owners. Samples/configuration are borrowed only for
// update(); const results describe the latest update, not a historical snapshot.
// Adaptation delegates to the existing Core transfer in its original order.
class ExplicitPicRefinement
{
public:
    struct Sample {
        const CellField<double>& numberAmount;
        const CellField<Vector2D>& electricField;
        const std::vector<double>* alpha;
        const UnitsContext& units;
        const NormalizationContext& normalization;
    };
    struct Report {
        double lambdaD_min_m = 0.0, h_min_m = 0.0, ratio_max = 0.0, ratio_mean = 0.0;
        int refine_cells = 0, coarsen_cells = 0;
        int debyeLowNeCount = 0, debyeBadCount = 0;
        double neFloor = 0.0;
    };
    explicit ExplicitPicRefinement(QuadtreeMesh& mesh);
    ExplicitPicRefinement(const ExplicitPicRefinement&) = delete;
    ExplicitPicRefinement& operator=(const ExplicitPicRefinement&) = delete;

    // Caller checks activation, electron availability and units. EMA updates
    // every call; true means the original heavy indicator/report stage ran.
    bool update(const AmrConfig& config, const Sample& sample, double temperatureEV,
        int step, int finalStep, QuadtreeMeshRefiner& refiner);
    void adaptMesh(std::vector<VertField<double>*>& vertexFields);
    struct Adaptation {
        bool enabled, embeddedEnabled, geometryRefinementActive, reportLevels;
        PicDomain& domain;
        PicAmrGeometry& geometry;
        PicMovingGeometry& moving;
        const std::vector<int>& lineToBoundaryId;
        PicFieldState& fields;
        ExplicitPicChargeDeposition& deposition;
        PicSpeciesState& particles;
        PicPoissonFieldSolve& poisson;
        PicCurrentWeights::Policy& coupling;
    };
    struct AdaptationReport {
        double seconds = 0.0;
        std::vector<int> levelsBefore, levelsAfter;
    };
    // Own the explicit AMR transition, including geometry, fields, particle
    // hosts and solver/coupling invalidation. Skips retain the last report.
    // Due calls reject inconsistent owners before adapting. Failures after
    // adaptation begins stop the session without rollback; no field solve or
    // particle cull is performed here. The original integrator does those later.
    bool adaptIfNeeded(const AmrConfig& config, int step, const Adaptation& state);
    const AdaptationReport& adaptationReport() const { return adaptationReport_; }
    void requestAdaptation() { amrPending = true; }
    void adaptationCompleted() { amrPending = false; amrPendingCount = 0; }
    bool pending() const { return amrPending; }
    const Report& report() const { return amrStats; }
    const CellField<double>& targetLevels() const { return amrTargetLevel; }
    const std::vector<RefinementIndicatorField>& indicators() const { return refinementIndicators.fields(); }

private:
    QuadtreeMesh& mesh;
    CellField<double> amrNeAvg, amrEAvg, amrIndicator, amrTargetLevel;
    bool amrNeAvgInit = false, amrEAvgInit = false, amrErefAvgInit = false;
    double amrErefAvg = 0.0;
    int amrPendingCount = 0;
    bool amrPending = false;
    PicRefinementIndicators refinementIndicators;
    Report amrStats;
    AdaptationReport adaptationReport_;
};

class ImplicitPicRefinement
{
public:
    struct Options {
        bool enabled = false, wantDebye = false, wantEField = false;
        int indicatorEvery = 1, persistenceSteps = 1, electronIndex = -1, maxChanges = 0;
        double temperatureEV = 2.0, debyeScale = 0.5, electricPercentile = 95.0;
        double gateElectricRelative = 0.0, gateDebyeRelative = 0.0;
    };
    struct Sample {
        const std::vector<CellField<double>>& numberAmounts;
        const CellField<Vector2D>& electricField;
        const std::vector<double>* alpha;
        const UnitsContext* units;
        const NormalizationContext& normalization;
    };
    explicit ImplicitPicRefinement(QuadtreeMesh& mesh);
    ImplicitPicRefinement(const ImplicitPicRefinement&) = delete;
    ImplicitPicRefinement& operator=(const ImplicitPicRefinement&) = delete;

    void update(const Options& options, const Sample& sample, int step, QuadtreeMeshRefiner& refiner);
    void adaptMesh(std::vector<VertField<double>*>& vertexFields, CellField<double>* potentialToTransfer);
    struct Adaptation {
        bool enabled, fastRemap, embeddedEnabled, geometryRefinementActive, resetTransients;
        int every;
        PicDomain& domain;
        PicAmrGeometry& geometry;
        PicMovingGeometry& moving;
        const std::vector<int>& lineToBoundaryId;
        PicFieldState& fields;
        ImplicitPicFieldHistory& history;
        ImplicitPicPolarization& polarization;
        PicSpeciesState& particles;
        PicPoissonFieldSolve& poisson;
        PicCurrentWeights::Policy& coupling;
        PicCurrentWeights& currentWeights;
        ImplicitPicOrbit& orbit;
        PicPoissonFieldSolve::Input fieldInput;
        const PicPoissonFieldSolve::ImplicitRebuildInput& rebuild;
        // Synchronous output-only hook, after successful field rebuild and
        // before completion/transient reset. Must not advance physical state.
        std::function<void(const QuadtreeMesh&)> refreshOutputCoordinates;
    };
    struct AdaptationReport {
        int cellsBefore = 0, cellsAfter = 0;
        bool useFastRemap = false;
        double meshAdaptS = 0.0, movingEbMeshChangedS = 0.0, faceGraphS = 0.0;
        double embeddedMetricsS = 0.0, fieldResetS = 0.0, relabelParticlesS = 0.0;
        double fieldRebuildS = 0.0, coordRefreshS = 0.0, seconds = 0.0;
    };
    // Full implicit adaptation/rebuild transition. Input views stay attached
    // to the same owners while their storage is rebuilt. Skips preserve the
    // previous report; invalid bindings reject before mesh mutation. Later
    // failures stop the session without rollback. Solver-cache invalidation
    // remains after the caller's existing AMR diagnostic output.
    bool adaptIfNeeded(int step, const Adaptation& state);
    const AdaptationReport& adaptationReport() const { return adaptationReport_; }
    // Called at the original AMR field-reset stage; change-gate history stays.
    void resetFieldsAfterMeshChange();
    void requestAdaptation();
    void adaptationCompleted() { amrPending = false; amrPendingCount = 0; }
    bool pending() const { return amrPending; }
    bool skippedByGate() const { return amrSkippedByGateLast; }
    const CellField<double>& targetLevels() const { return amrTargetLevel; }
    const std::vector<RefinementIndicatorField>& indicators() const { return refinementIndicators.fields(); }

private:
    QuadtreeMesh& mesh;
    CellField<double> amrTargetLevel, amrIndicator;
    bool amrPending = false, amrSkippedByGateLast = false;
    int amrPendingCount = 0;
    PicRefinementIndicators refinementIndicators;
    std::vector<double> prevAmrEFieldIndicator, prevAmrDebyeIndicator;
    AdaptationReport adaptationReport_;
};
