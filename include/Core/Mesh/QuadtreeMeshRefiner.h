#pragma once

#include <string>
#include <cstdint>
#include <functional>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "Core/Field/CellField.h"
#include "Core/IO/DictionaryGeometry.h"
#include "Core/IO/DictionaryMesh.h"
#include "Core/Mesh/MeshRefinementIndicator.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Geometry/ImmersedGeometry.h"

class QuadtreeMeshRefiner
{
public:
    struct BudgetOptions
    {
        int cellBudget = 0; // 0 = no budget
        const std::vector<unsigned char>* excludeMask = nullptr;
    };

    // Explicit current-time, current-leaf inputs; never cached across adapt().
    // phaseByCell uses 0-based phase IDs, -1 for unresolved mixed cells.
    // Flow indicators use finite nonnegative values and positive references in
    // matching units. The producing module owns their physical interpretation.
    struct CurrentFields
    {
        std::vector<int> cellLabels;
        std::vector<int> phaseByCell;
        std::vector<RefinementIndicatorField> indicators;
    };

    // Synchronous source for repeated runtime refinement sweeps. The source
    // mesh, indicators and already topology-closed protection mask stay fixed.
    // Phase classification belongs to the caller; receives target center and
    // the containing source leaf index. An absent callback means unresolved.
    struct CurrentSource
    {
        const QuadtreeMesh& mesh;
        const std::vector<RefinementIndicatorField>& indicators;
        const std::vector<std::uint8_t>& protectedCells;
        std::function<int(Vector2D, int)> phaseAt;
    };

    explicit QuadtreeMeshRefiner(const DictionaryMesh& dictMesh);

    bool attachGeometry(const DictionaryGeometry* dictGeometry,
                        ImmersedGeometry<QuadtreeMesh>* geometry);
    void clearGeometry();

    // Re-discretize boundary polylines from the current live geometry lines.
    // Call this after geometry has moved to keep geometry-driven AMR in sync.
    void refreshBoundaryPolylines();

    // Replace ONE configured boundary's refinement geometry with an owned
    // current-time polyline set (code coordinates). For transported interfaces
    // whose authoritative scalar is no longer the initial geometry lines.
    // Uses the existing distance/target-level algorithm, not a second refiner.
    // Requires an attached, resolved boundary rule; rejects empty/nonfinite or
    // degenerate segments before changing any rule. Other boundaries retain
    // their data. refresh/clear/attach intentionally replace this override.
    // These polylines select resolution only, never replace the physical scalar.
    void replaceBoundaryPolylines(int boundaryId,
        const std::vector<std::vector<Vector2D>>& polylines);

    // Clear cached boundary polylines, forcing geometry-band AMR to fall back
    // to level-set |psi| distances. Use this when iGeometry.lines are NOT in
    // the current geometry frame (e.g. after advected-path movingEB updates).
    void clearBoundaryPolylines();

    bool geometryActive() const { return geometryActive_; }
    bool phaseActive() const { return phaseActive_; }

    // The optional runtime policy permits one level of
    // coarsening outside all protected geometry bands/regions/phase minima.
    // Boundary coarsening uses the prospective parent's band (twice the
    // child's width) so a stationary interface does not refine/coarsen-cycle.
    // False retains the historical initialization/PIC policy exactly. Runtime
    // flow physics requires explicit CurrentFields, never PIC indicators or
    // initial vertex phase labels. Rebuild inputs after each adapt() sweep.
    // PIC indicators are borrowed only for this call; no process-wide lookup
    // or retained input pointer. Missing/wrong-sized indicators retain the
    // original warning/ignore policy. Do not combine them with CurrentFields.
    bool buildTargetLevels(const QuadtreeMesh& mesh,
                           CellField<double>& targetLevel,
                           CellField<double>* indicatorOut = nullptr,
                           const BudgetOptions* budget = nullptr,
                           bool coarsenOutsideGeometry = false,
                           const CurrentFields* current = nullptr,
                           const std::vector<RefinementIndicatorField>* picIndicators = nullptr);

    // Map current inputs to this sweep's leaves, retaining every covered old
    // indicator maximum on coarsening and parent values on refinement. Apply
    // the existing runtime policy, lock protected levels and ask the original
    // topology rules whether any adaptation can start. No mesh/state clone,
    // physical-field transfer, retained callback or acceptance occurs here.
    // Output fields are work buffers and may be changed before a rejection.
    bool buildMappedTargetLevels(const QuadtreeMesh& target, const CurrentSource& source,
        CellField<double>& targetLevel, CellField<double>* indicatorOut = nullptr);

private:
    struct GeoBoundaryRuleRuntime
    {
        int boundaryId = -1;
        int levelSetIndex = -1;
        double cGamma = 1.0;
        int maxLevel = -1;
        std::vector<int> lineIds;
        std::vector<std::vector<Vector2D>> linePolylines;
    };

    struct GeoRegionRuleRuntime
    {
        std::string name;
        DictionaryMesh::GeometryRefinementConfig::RegionShape shape =
            DictionaryMesh::GeometryRefinementConfig::RegionShape::kUnknown;
        Vector2D min;
        Vector2D max;
        int maxLevel = -1;
    };

    const DictionaryMesh& dictMesh_;
    const DictionaryGeometry* dictGeometry_ = nullptr;
    ImmersedGeometry<QuadtreeMesh>* geometry_ = nullptr;

    bool geometryActive_ = false;
    bool regionActive_ = false;
    bool boundaryActive_ = false;
    bool phaseActive_ = false;

    std::vector<GeoBoundaryRuleRuntime> geoBoundaryRules_;
    std::vector<GeoRegionRuleRuntime> geoRegionRules_;
    std::unordered_map<int, DictionaryMesh::GeometryRefinementConfig::PhaseRule> geoPhaseRules_;

    std::unordered_set<std::string> warnedMissingIndicators_;
    bool warnedMissingGeometry_ = false;
};
