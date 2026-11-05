#pragma once

#include <cstdint>
#include <vector>

namespace LAMP::Geometry
{
// Explicit scalar authority on a subset of a complete Cartesian leaf layout.
// Entries use the mesh's native indices, not labels or Cartesian positions.
// Each active cell consumes its four nodal DOFs, shared except for explicitly
// declared corner sectors below. No solid extension, derivative reconstruction,
// interface selection or boundary law is implied.
// A snapshot copies these values and scales each consumed physical derivative
// once. Inactive-cell and unused-node payloads are never scalar query results.
struct ActiveHermitePatch
{
    struct NodalDof
    {
        double value = 0.0;
        double dx = 0.0, dy = 0.0, dxy = 0.0;
    };

    // A cell-local scalar extension at one declared zero boundary corner.
    // IDs are native indices. The vertex must be internal to the Cartesian
    // rectangle, with exactly three active incident cells and one inactive;
    // the selected active cell must share an edge with that inactive cell.
    // Its value remains the exact base zero. Shared active edges retain the
    // same value/tangential derivative traces (C0); normal/mixed derivatives
    // may differ, so no cross-sector C1 assertion is made. The selected cell
    // must be certified phaseSign everywhere except this one zero corner.
    struct CornerSector
    {
        int cellIndex = -1;
        int vertexIndex = -1;
        NodalDof dof;
        int phaseSign = 0;
    };

    // Full cellIndexNum length; every entry must be binary.
    std::vector<std::uint8_t> activeCells;
    // Full vertIndexNum length. Only active-cell corners must be finite;
    // unused entries may contain NaN and are not read during construction.
    std::vector<NodalDof> nodes;
    // At most one sector per native vertex and per active cell. Empty keeps
    // the original shared-DOF scalar. This is not a contact-angle condition.
    std::vector<CornerSector> cornerSectors;
};
} // namespace LAMP::Geometry
