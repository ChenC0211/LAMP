#pragma once
#include "Core/Geometry/ActiveHermitePatch.h"
#include <cstddef>
#include <vector>

class QuadtreeMesh;
namespace ActiveHermiteTransfer
{
struct Result
{
    LAMP::Geometry::ActiveHermitePatch patch;
    std::vector<int> sourceToTargetVertex;
    // All old leaves with positive-area overlap, including inactive leaves.
    std::vector<std::vector<int>> sourceCellsByTarget;
    // Absolute scalar L-infinity bound over each active target leaf. Computed
    // AFTER native hanging constraints from exact rational corner-jet defects
    // on every old/new overlap rectangle. Inactive entries are NaN.
    std::vector<double> scalarDifferenceBound;
    double maxScalarDifferenceBound = 0;
    // Refined descendants of a declared pure corner-sector cell may need a
    // new scalar extension: one shared fine jet cannot retain both old normal
    // traces. Nonzero entries certify SAME strict phase before/after, except
    // the retained declared zero corner. They do not claim scalar equality.
    std::vector<int> pureSectorExtension;
    double maxNonExtensionDifferenceBound = 0;
    std::size_t copiedVertices = 0, sampledVertices = 0;
};
// Zero-time bicubic scalar remap on the same Cartesian background lattice and
// fixed active domain. Interface leaves MAY refine/coarsen. Retained regular
// jets are copied; new independent jets sample the old constrained polynomial.
// Dependent target payload stays NaN and is resolved by the existing authority.
// Coarsening is interpolation, not an assertion of polynomial equality.
// The returned bound is for the scalar, NOT distance, volume, curvature or
// topology. Callers must accept shared geometry/topology and independently
// remap conserved quantities before publishing any coupled state.
// Corner sectors retain their original regular zero vertex and locally
// same-level reentrant topology. New seam jets prefer the non-sector chart;
// only proven pure descendants may reextend their scalar. All scalar defects
// remain reported. Lost contacts, invalid phase or unsupported layouts reject.
Result prepare(const QuadtreeMesh& sourceMesh,
    const LAMP::Geometry::ActiveHermitePatch& source,
    const QuadtreeMesh& targetMesh, const std::vector<std::uint8_t>& targetActive,
    bool axisymmetric);
}
