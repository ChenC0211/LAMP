#pragma once

#include "Core/Geometry/ActiveQuarticHermitePatch.h"
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <vector>

class QuadtreeMesh;

namespace ActiveQuarticHermiteTransfer
{
// A phase-preserving bulk polynomial changed through new shared constraints.
// The caller may rebuild a DIFFERENT topology retaining its native support;
// the rejected transfer itself is never accepted or relaxed.
class ProtectedPolynomialChange : public std::runtime_error
{
public:
    ProtectedPolynomialChange(int sourceCell, int targetCell);
};

// Retain the geometric carriers and vertex incidence of requested preserved
// Q4 polynomials. This augments topology requests, NOT the scalar acceptance
// mask. Uses existing native hanging identities; creates no constraint graph.
std::vector<std::uint8_t> topologyProtection(const QuadtreeMesh& mesh,
    const std::vector<std::uint8_t>& preservedCells);

struct Report
{
    std::size_t identicalActiveCells = 0;
    std::size_t certifiedChangedCells = 0;
    std::size_t protectedInterfaceCells = 0;
    std::size_t copiedVertices = 0, copiedMidpoints = 0, copiedCenters = 0;
};

struct Result
{
    LAMP::Geometry::ActiveQuarticHermitePatch patch;
    // Full SOURCE native index arrays, not labels. A cell entry is nonnegative
    // only when its bounds AND authoritative constrained Q4 polynomial are
    // exactly unchanged. Its retained local chart can then be reused directly.
    // Vertex entries identify equal stored coordinates, not scalar/DOF equality.
    std::vector<int> sourceToTargetCell, sourceToTargetVertex;
    Report report;
};

// Zero-time remeshing of the existing explicit Q4 scalar, not material motion,
// reinitialization, a curvature fit or a new contact law. Both meshes and the
// source patch are borrowed read-only for this call; the result owns its data.
// Reuses the native Q4 constraints, polynomial restriction and exact Bernstein
// certificates. Independent retained data are copied. Newly required data are
// sampled from the old polynomial at the EXACT target vertex/midpoint/center,
// then represented as doubles; bulk scalar interpolation is not claimed exact.
//
// All source cells without a strict CLOSED single-phase certificate must keep
// their bounds and their complete constrained polynomial exactly. The optional
// full native-cell binary mask adds a caller-selected support halo to this rule.
// Every other changed target cell must be certified to have the same strict
// phase as ALL overlapping source cells. Thus accepted remeshing changes no
// zero set or phase volume; no volume correction is used to hide transfer error.
// Refinement/coarsening of the interface itself is deliberately rejected: the
// caller must maintain a finest-level band ahead of the moving interface.
// Fixed active/solid coverage, radial symmetry and corner sectors are retained.
// Failure publishes nothing and changes neither mesh nor source patch. This
// does not transfer velocity/pressure/charge or accept an App AMR time step.
Result prepare(const QuadtreeMesh& sourceMesh,
    const LAMP::Geometry::ActiveQuarticHermitePatch& sourcePatch,
    const QuadtreeMesh& targetMesh, const std::vector<std::uint8_t>& targetActiveCells,
    bool radialAxisSymmetry = false,
    const std::vector<std::uint8_t>& protectedSourceCells = {});
}
