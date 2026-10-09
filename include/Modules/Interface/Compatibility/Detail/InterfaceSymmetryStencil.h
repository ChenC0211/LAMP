#pragma once
#include "Core/Geometry/LevelSetInterfacePanels.h"

namespace InterfaceSymmetryStencil
{
struct Reflection { bool even=false,axis=false;double origin=0; };
// Select three ACTUAL columns and an even moment basis about the declared
// section. No ghost volumes/cells are created; physical RZ weights remain in
// InterfaceHeightMoments. Simultaneous opposite-boundary support is unresolved.
inline Reflection select(Vector2D lower,Vector2D upper,int tangent,int extent,
    int& first,int& count,bool radialAxis,
    const std::vector<LevelSetInterfacePanels::SymmetryPlane>& planes)
{
    const bool low=first<0,high=first+count>extent;
    if(low==high || extent<3)return {};
    const double coordinate=tangent?(low?lower.y:upper.y):(low?lower.x:upper.x);
    const bool axis=radialAxis && tangent==0 && low && coordinate==0;
    bool admitted=axis;
    for(const auto& p:planes)admitted|=p.normalAxis==tangent && p.coordinate==coordinate;
    if(!admitted)return {};
    first=low?0:extent-3;count=3;
    return {true,axis,coordinate};
}
}
