#pragma once

#include "Modules/Interface/HermiteGalsTrace.h"
#include "Core/Geometry/LevelSetHermiteGeometry.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

// Private implementation shared by the vertex and edge stages. No new state,
// interpolation, clock or geometry authority is introduced here.
namespace LAMP::Interface::Detail
{
inline void requireFirstJet(bool ok,const char* message)
{ if(!ok)throw std::invalid_argument(message); }

inline void validateFirstJetSource(const QuadtreeMesh& mesh,const HermiteGalsState& source,
    const LevelSetHermiteGeometry::Snapshot& scalar,const LevelSetBoundaryTrace::Domain& domain,
    const HermiteGalsVelocityView& previous,const HermiteGalsVelocityView& arrival)
{
    requireFirstJet(previous.geometryIdentity()==source.identity(),"previous velocity geometry differs from source");
    const auto sourceTime=source.timeLayer(),previousTime=previous.timeLayer(),arrivalTime=arrival.timeLayer();
    requireFirstJet(previousTime.step==sourceTime.step && previousTime.begin==sourceTime.begin &&
        previousTime.end==sourceTime.end && previousTime.position==sourceTime.position,
        "previous velocity time layer differs from source");
    const double oldTime=sourceTime.time(),newTime=arrivalTime.time();
    requireFirstJet(arrivalTime.begin==oldTime,"arrival layer does not begin at source time");
    requireFirstJet(newTime>=oldTime && (newTime>oldTime || arrivalTime.position==HermiteGalsTimeLayer::Position::Initial),
        "arrival time interval is unresolved or reversed");
    const auto same=[](Vector2D a,Vector2D b){return a.x==b.x && a.y==b.y;};
    requireFirstJet(same(domain.minCoord,mesh.minCoord) && same(domain.maxCoord,mesh.maxCoord) &&
        same(scalar.minCoord(),mesh.minCoord) && same(scalar.maxCoord(),mesh.maxCoord),
        "trace domain must match the native mesh and source snapshot bounds");
    requireFirstJet(domain.radialAxisSymmetry==scalar.radialAxisSymmetry(),"trace/source axis semantics differ");
    Vector2D minimum{std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity()};
    for(int c=0;c<mesh.cellIndexNum();++c){
        const auto size=mesh.cells(mesh.cellIndexToLabel[c]).size();
        minimum.x=std::min(minimum.x,size.x);minimum.y=std::min(minimum.y,size.y);
    }
    requireFirstJet(std::isfinite(domain.cellSize.x) && std::isfinite(domain.cellSize.y) &&
        domain.cellSize.x>0 && domain.cellSize.y>0 && domain.cellSize.x<=minimum.x && domain.cellSize.y<=minimum.y,
        "trace cell scale must not exceed the minimum native leaf size");
}

struct HermiteGalsFirstJet { double phi; Vector2D gradient; };

inline HermiteGalsFirstJet pullbackFirstJet(const LevelSetHermiteGeometry::Snapshot& scalar,
    const LevelSetBoundaryTrace::Domain& domain,Vector2D point,
    const HermiteGalsVelocityView& previous,const HermiteGalsVelocityView& arrival,const char*& stage)
{
    stage="trace";
    const auto traced=traceHermiteGalsWithJacobian(domain,point,previous,arrival);
    const auto jacobian=traced.interiorJacobian();
    stage="source sample";
    const auto old=scalar.scalarJetNative(traced.departurePoint());
    requireFirstJet(std::isfinite(old.value) && std::isfinite(old.dx) && std::isfinite(old.dy),
        "nonfinite source first jet");
    stage="first-jet pullback";
    const auto gradient=jacobian.transposeTimes({old.dx,old.dy});
    requireFirstJet(std::isfinite(gradient.x) && std::isfinite(gradient.y),"nonfinite transported gradient");
    return {old.value,gradient};
}
} // namespace LAMP::Interface::Detail
