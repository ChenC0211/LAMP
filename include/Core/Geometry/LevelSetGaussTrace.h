#pragma once
#include "Core/Base/Vector2D.h"
#include <Eigen/Core>
#include <functional>

// Frozen autonomous first-jet map. Patch exits, event durations and the total
// duration are differentiated together. A variable event duration does not
// confer a general exact volume-preservation guarantee. Only bracketed,
// transverse exits and the declared floating endpoint incidences are admitted.
namespace LevelSetGaussTrace {
struct VelocityJet {Vector2D value,dx,dy;};
struct TimeVelocityJet {Vector2D value,dx,dy,dt;};
struct SmoothPatch {Vector2D lower,upper;};
using Query=std::function<VelocityJet(Vector2D)>;
using PatchQuery=std::function<SmoothPatch(Vector2D,Vector2D)>;
inline bool contains(const SmoothPatch& p,Vector2D x) {
    return x.x>=p.lower.x && x.x<=p.upper.x && x.y>=p.lower.y && x.y<=p.upper.y;
}
struct Map {
    Vector2D point,timeDerivative,displacement;
    Eigen::Matrix2d derivative,naiveDerivative;
    int events=0,endpointIncidences=0,cornerIncidences=0;
    double eventResidual=0;
    // Coordinate incidence scale already admitted by the native branch logic.
    // This is not a certified ODE truncation-error bound.
    double coordinateResolution=0;
};
Map advance(const Query&,const PatchQuery&,Vector2D begin,double duration,int blocks=8);
// C1 spatial field on one physical domain. The query time is elapsed from
// the fixed beginning of this step. Differentiate the actual collocation map
// in both starting position and duration, including movement of stage times.
// Physical boundary exits are unsupported; no numerical-knot event restarts.
using TimeQuery=std::function<TimeVelocityJet(Vector2D,double)>;
Map advanceTimeDependent(const TimeQuery&,const SmoothPatch&,Vector2D begin,double duration,int blocks=8);
} // namespace LevelSetGaussTrace
