#pragma once
#include "Core/Base/Vector2D.h"
#include <array>
#include <vector>

namespace LAMP::Hydrodynamics::Detail
{
// Pure prepared rows on the ORIGINAL native U. No virtual field/state owner.
struct HydroVirtualWeight {int velocityDof=-1;double weight=0;};
using HydroVirtualRow=std::vector<HydroVirtualWeight>;
struct HydroVirtualCellInput {
    std::array<int,2> lower{},upper{}; // finest logical vertices, extent1 or2
    std::array<HydroVirtualRow,2> value;
    std::array<std::array<HydroVirtualRow,2>,2> gradient; // [component][physical derivative axis]
};
struct HydroVirtualFaceInput {
    int velocityDof=-1,normalAxis=0,negativePressure=-1,positivePressure=-1;
    std::array<int,2> lower{},upper{}; // constant normal coordinate, positive tangent interval
    HydroVirtualRow tangentGradient; // physical derivative; caller authenticates boundary law
};
struct HydroVirtualFaceProlongation {
    int nx=0,ny=0;
    // x-face y*(nx+1)+x; y-face y*nx+x, canonical +x/+y means.
    std::array<std::vector<HydroVirtualRow>,2> rows;
    std::vector<int> cellPressure; // fine cell y*nx+x -> original pressure
};
// Complete square XY domain tiled by actual1x1/2x2 leaves only. Cell input
// order is the original pressure order. On a length2h native face, the two
// means are U +/- h/2*d_t U. Four coarse internal faces retain the predicted
// circulation and have one fourth of the ORIGINAL integrated divergence.
// Boundary tangent rows come from the same caller source: zero wall normals
// must imply zero boundary slopes; unbound affine data must not be zeroed.
HydroVirtualFaceProlongation prepareHydroVirtualFaceProlongation(int nx,int ny,Vector2D spacing,
    int nativeVelocityCount,const std::vector<HydroVirtualCellInput>& cells,
    const std::vector<HydroVirtualFaceInput>& faces);
}
