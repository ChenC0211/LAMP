#pragma once

#include <cstdint>
#include "Core/Base/Vector2D.h"

class BackgroundMesh;
class Cell;

namespace Quantum::Mesh
{
struct DyadicCellAddress
{
    std::uint64_t i=0,j=0;
    unsigned level=0;
};

// One physical coordinate for one logical Cartesian address. Refinement
// history does not participate in evaluation: (i,j,L) and (2i,2j,L+1)
// produce identical coordinates. Centers use odd addresses one level finer.
// The domain endpoints are retained verbatim. This does not assert that
// subtracting two rounded coordinates reproduces the nominal spacing exactly.
// Implemented out of line with strict FP/no contraction, not in mesh callers
// with potentially different compiler arithmetic options.
class DyadicGridCoordinates
{
public:
    DyadicGridCoordinates(Vector2D minimum,Vector2D maximum,int nx,int ny);
    Vector2D point(std::uint64_t i,std::uint64_t j,unsigned level) const;
    Vector2D cellSize(unsigned level) const;
private:
    Vector2D minimum_,maximum_,spacing_;
    std::uint64_t nx_,ny_;
    double axis(std::uint64_t index,unsigned level,int direction) const;
};

// Recover the address from existing root/parent/son topology. No new mesh
// identity or serialized solver state is introduced. Noncanonical roots or
// malformed/cyclic ancestry reject before constructing any child geometry.
DyadicCellAddress cellAddress(const BackgroundMesh& mesh,const Cell& cell);
}
