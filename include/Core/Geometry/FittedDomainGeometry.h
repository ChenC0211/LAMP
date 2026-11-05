#pragma once

#include <vector>
#include <memory>
#include "Core/IO/DictionaryGeometry.h"

class QuadtreeMesh;
struct FaceSegment;
template<class Mesh> class ImmersedGeometry;

// Fixed-domain geometry binding, with no electric/fluid algorithm or file IO.
namespace FittedDomainGeometry
{
struct InitialPhaseSamplerData;
struct FixedSolidSelection;
// Immutable t=0 CSG values from the existing owned primitive samplers. No
// evolving level-set initialization, zero displacement or interface state.
// Level-set indices are zero based; public physical phase IDs are one based.
class InitialPhaseSampler
{
public:
    InitialPhaseSampler(QuadtreeMesh&, const DictionaryGeometry&);
    double levelSetValue(int levelSetIndex, Vector2D point) const;
    double phaseValue(int phaseId, Vector2D point) const;
    int phaseAt(Vector2D point) const; // -1 on uncovered exact interfaces
    int phaseCount() const;
private:
    std::shared_ptr<const InitialPhaseSamplerData> data_;
    friend FixedSolidSelection selectFixedSolids(const QuadtreeMesh&,
        const DictionaryGeometry&, const InitialPhaseSampler&,
        const std::vector<int>&, double);
};
struct FixedSolidSelection
{
    // Native leaf-index order: 0 participates, otherwise selected phase ID.
    std::vector<int> excludedPhaseByCell;
    // Same ordering as FaceSegmentGraph::build(mesh) with every leaf active.
    // 0 for nonexposed segments; declared ID on each participating domain or
    // selected-solid face. Both directions of internal fluid faces remain 0.
    std::vector<int> boundaryIdsBySegment;
};
// Sufficient whole-open-leaf certificate, deliberately conservative. Bounds
// of ALL original sign primitives and closure connectors enclose every possible
// sign change; interval Boolean evaluation must prove each selected phase
// wholly present/absent. Unsupported/cutting geometry is rejected before a
// mask is returned, including a closed solid hidden between cell samples.
// Exposed segments then reuse bindFace's full straight coverage and adjacency.
// Moving fluid/fluid curves may cancel from the selected solid expression.
// Recompute on each actual mesh; the immutable initial sampler may be reused.
FixedSolidSelection selectFixedSolids(const QuadtreeMesh&, const DictionaryGeometry&,
    const InitialPhaseSampler&, const std::vector<int>& excludedPhaseIds,
    double tolerance = 1e-12);
// Original cell-center interpolation and phase-expression order. Returns
// 1-based phase IDs or -1 when no phase covers the center; not a whole-cell
// certificate. A fitted domain must separately verify full solid coverage.
std::vector<int> cellPhases(const QuadtreeMesh&, const ImmersedGeometry<QuadtreeMesh>&);
// Original full straight-segment coverage and ambiguity checks. The excluded
// phase is a selected solid, regardless of its electrical material model.
// Historical failure reason strings are retained during ownership migration.
const DictionaryGeometry::BoundaryEntity& bindFace(const DictionaryGeometry&,
    const FaceSegment&, int ownerPhase, int excludedPhase, double tolerance);
}
