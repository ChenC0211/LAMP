#pragma once

#include <array>

namespace Quantum::Geometry
{
// A CornerSector owns a separate replacement jet, not a regular-node DOF.
// HoldFixed differentiates shared patch data while ALL declared sector jets stay
// unchanged. The caller must still impose pin/shared-trace admissibility.
// No variation of the independently declared sector data is implemented.
enum class CornerSectorVariation
{
    RequireAbsent,
    HoldFixed
};

// A floating first variation with respect to one native REGULAR vertex's
// physical [value, dx, dy, dxy] data. Hanging-node chain rules have already
// been applied. No mesh, source patch or mutable cache is borrowed.
// Pinning, axis symmetry and the choice of free derivatives are separate
// admissibility constraints, imposed by the caller before solving/updating.
struct ActiveHermiteNodalSensitivity
{
    int vertexIndex = -1;
    std::array<long double, 4> coefficients{};
};

// Quartic shared-data pullback. Different independent identities and layout
// from the cubic regular-node result above; never convert between them.
struct ActiveQuarticHermiteSensitivity
{
    enum class Kind { Vertex, EdgeMidpoint, CellCenter };
    Kind kind=Kind::Vertex;
    // Native vertex/cell index, or source ActiveQuarticHermitePatch edge index.
    int index=-1;
    // Vertex: value/dx/dy/dxy = 0..3; edge: value/normal = 0..1; center: 0.
    int component=0;
    long double coefficient=0.0L;
};
} // namespace Quantum::Geometry
