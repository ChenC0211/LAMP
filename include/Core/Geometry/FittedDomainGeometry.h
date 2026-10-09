#pragma once

#include <vector>
#include "Core/IO/DictionaryGeometry.h"

class QuadtreeMesh;
struct FaceSegment;
template<class Mesh> class ImmersedGeometry;

// Fixed-domain geometry binding, with no electric/fluid algorithm or file IO.
namespace FittedDomainGeometry
{
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
