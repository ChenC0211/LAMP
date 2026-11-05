#include "Modules/Interface/HermiteGalsTrace.h"

#include <cmath>
#include <stdexcept>
#include <string>
#include <utility>

namespace LAMP::Interface
{
namespace
{
std::string identities(HermiteGalsIdentity previous, HermiteGalsIdentity arrival)
{
    return "[Interface/HermiteGalsTrace previous=" + std::to_string(previous.value())
        + " arrival=" + std::to_string(arrival.value()) + "] ";
}
void require(bool ok, const char* message)
{ if (!ok) throw std::invalid_argument(message); }
}

HermiteGalsVelocityView::HermiteGalsVelocityView(HermiteGalsIdentity geometry,
    HermiteGalsTimeLayer time, Query query)
    : geometry_(geometry), time_(time), query_(std::move(query))
{
    try {
        require(geometry_.value() != 0, "unbound geometry identity");
        time_.time();
        require(static_cast<bool>(query_), "missing paired velocity/derivative query");
    } catch (const std::exception& error) {
        throw std::invalid_argument("[Interface/HermiteGalsTrace geometry="
            + std::to_string(geometry_.value()) + "] velocity view: " + error.what());
    }
}

HermiteGalsTraceResult::HermiteGalsTraceResult(Vector2D point,
    const HermiteGalsVelocityView& previous, const HermiteGalsVelocityView& arrival,
    LevelSetBoundaryTrace::JetResult trace)
    : arrivalPoint_(point), previousGeometry_(previous.geometryIdentity()),
      arrivalGeometry_(arrival.geometryIdentity()), previousTime_(previous.timeLayer()),
      arrivalTime_(arrival.timeLayer()), trace_(std::move(trace)) {}

double HermiteGalsTraceResult::departureTime() const
{
    // A completed full trace retains the supplied old physical time exactly;
    // avoid recovering it through a potentially rounded (new - dt) subtraction.
    return boundaryHit() ? arrivalTime_.time() - elapsed() : previousTime_.time();
}

LevelSetBoundaryTrace::Jacobian HermiteGalsTraceResult::interiorJacobian() const
{
    if (!trace_.derivative)
        throw std::logic_error(identities(previousGeometry_, arrivalGeometry_)
            + "interior Jacobian unavailable: boundary event requires a boundary-history tangent");
    return *trace_.derivative;
}

HermiteGalsTraceResult traceHermiteGalsWithJacobian(
    const LevelSetBoundaryTrace::Domain& domain, Vector2D arrivalPoint,
    const HermiteGalsVelocityView& previous, const HermiteGalsVelocityView& arrival)
{
    try {
        const double oldTime = previous.time_.time(), newTime = arrival.time_.time();
        require(arrival.time_.begin == oldTime,
            "previous velocity time does not equal arrival layer begin");
        const double dt = newTime - oldTime;
        require(std::isfinite(dt) && dt >= 0, "invalid resolved trace interval");
        require(dt > 0 || arrival.time_.position == HermiteGalsTimeLayer::Position::Initial,
            "advancing/midpoint time layer is unresolved in floating point");
        // Core's primal and tangent specializations share a single RK2 control
        // flow, time interpolation, domain guard and discrete event policy.
        auto result = LevelSetBoundaryTrace::traceWithJacobian(
            domain, arrivalPoint, dt, arrival.query_, previous.query_);
        return HermiteGalsTraceResult(arrivalPoint, previous, arrival, std::move(result));
    } catch (const std::exception& error) {
        throw std::runtime_error(identities(previous.geometry_, arrival.geometry_)
            + "trace: " + error.what());
    }
}
} // namespace LAMP::Interface
