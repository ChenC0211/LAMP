#pragma once
#include "Modules/Plasma/Boundary/PicEmbeddedGeometry.h"
#include "Core/Geometry/ImmersedGeometry.h"

// Owns AMR geometry refresh policy/history. Capture before mesh.adapt(), then
// refresh after the existing MovingEB onMeshChanged notification. The report
// is read-only and valid until the next capture/refresh. Mesh/geometry/embedded
// state are borrowed; face graph, cut metrics and solver invalidation follow
// in their original driver order after this operation.
class PicAmrGeometry
{
public:
    using Mode = PicEmbeddedGeometry::Mode;
    using Geometry = PicEmbeddedGeometry::Geometry;
    enum class RefreshMode { Full, Adaptive, Skip };
    struct Options
    {
        RefreshMode mode = RefreshMode::Full;
        double amrGeomRefreshChangeRelTh = 0.2;
        bool amrGeomRefreshVerbose = false;
        bool amrGeomRefreshAllowWithMovingEb = false;
        bool amrReinitActiveMovingOnly = false;
        bool amrReinitHybridSensitiveInactive = false;
        bool amrReinitHybridUseDelta = true;
        int amrReinitHybridNearZeroMin = 1;
        int amrReinitHybridNarrowBandMin = 1;
        int amrReinitHybridNearZeroDeltaMin = 1;
        int amrReinitHybridNarrowBandDeltaMin = 1;
        bool amrReinitHybridIncludeMaxIter = true;
        int amrReinitHybridMaxIterThreshold = 500;
        int amrReinitActiveMovingMaxStreak = 0;
    };
    struct Adaptation
    {
        int step = 0;
        double cellCountChange = 0.0;
        bool embeddedEnabled = false;
        bool refinementActive = false;
        bool phaseActive = false;
        bool movingUpdated = false;
    };
    struct Report
    {
        double geometryRefreshS = 0.0;
        double geomInitializeLevelSetsS = 0.0;
        double geomPdeReinitS = 0.0;
        double geomRefreshDerivedS = 0.0;
        double geomCheckCutInfoS = 0.0;
        double geomRebuildCompositeS = 0.0;
        double geomMarkPhasesS = 0.0;
        bool geometryReinitSkipped = false;
        bool geometryRefreshUsedExactMovingEb = false;
        std::vector<int> amrRawActiveMovingLevelSets;
        std::vector<int> amrSelectiveMovingLevelSets;
        std::vector<int> amrHybridSensitiveInactiveLevelSets;
    };
    PicAmrGeometry(Mode mode, Options options);
    PicAmrGeometry(const PicAmrGeometry&) = delete;
    PicAmrGeometry& operator=(const PicAmrGeometry&) = delete;
    void captureBeforeAdapt(const Geometry* geometry);
    void refresh(QuadtreeMesh& mesh, Geometry* geometry, PicEmbeddedGeometry& embedded,
        const MovingEBController* moving, const Adaptation& adaptation);
    const Report& report() const { return report_; }
private:
    friend class ImplicitPicRefinement;
    friend class ExplicitPicRefinement;
    const Mode mode_;
    const Options options_;
    int selectiveStreak_ = 0;
    std::vector<Geometry::LevelSetUpdateDiagnostics> preAdaptDiagnostics_;
    Report report_;
};
