#pragma once
#include "Core/Base/Vector2D.h"
#include <array>
#include <utility>
#include <vector>

class QuadtreeMesh;
class NativeLeafDomain;
namespace HydroFittedDetail
{
using MomentumTrace=std::vector<std::pair<int,double>>;
struct NativeVelocityTrace
{
    MomentumTrace value;
    std::array<MomentumTrace,2> gradient;
};
// Existing first-order native face reconstruction, shared by fitted-boundary
// and complete conservative stress assembly. No velocity DOFs or ghost mesh.
std::vector<int> nativeTraceSupport(const QuadtreeMesh&,const NativeLeafDomain&,const std::vector<int>& seed);
NativeVelocityTrace nativeVelocityTrace(const QuadtreeMesh&,const NativeLeafDomain&,
    const std::vector<int>& cells,int component,Vector2D point);
}
