#pragma once

#include "Core/Base/Enumeration.h"
#include "Core/Base/Vector2D.h"
#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <type_traits>

// Boundary events for the existing two-time-level explicit midpoint backtrace.
// This is a bounded numerical characteristic, not a certified first-hit oracle
// for an arbitrary oscillatory velocity callback or its continuous ODE. Callers
// must supply resolved, continuous interior/boundary velocity traces. No field,
// boundary value, exterior velocity or contact-line model is owned here.
namespace LevelSetBoundaryTrace
{
struct Domain
{
    Vector2D minCoord, maxCoord, cellSize;
    bool radialAxisSymmetry = false;
    // Optional selected-fluid coverage. Rejects internal wall crossings;
    // it does not invent a boundary history or reflect an RK stage.
    std::function<bool(Vector2D,Vector2D)> activeSegment;
};
struct Result
{
    Vector2D point;
    double elapsed = 0.0; // Backward from newTime; boundary time is newTime-elapsed.
    std::optional<Boundary> side;
    bool boundaryHit = false;
    // Boundary resolution from an unsampled midpoint inside the existing
    // roundoff-sized space/time bracket, never a completed interior RK stage.
    bool roundoffMidpointEvent = false;
};
using Velocity = std::function<Vector2D(const Vector2D&)>;

// Rows are output components, columns are physical input coordinates:
// [du_x/dx du_x/dy; du_y/dx du_y/dy]. Velocity values and derivatives must
// come from the SAME interpolation branch; no independent gradient fit.
struct Jacobian
{
    double xx=0.0, xy=0.0, yx=0.0, yy=0.0;
    static Jacobian identity() { return {1,0,0,1}; }
    Vector2D transposeTimes(Vector2D gradient) const
    { return {xx*gradient.x+yx*gradient.y,xy*gradient.x+yy*gradient.y}; }
};
struct VelocityJet { Vector2D value; Jacobian derivative; };
using DifferentiableVelocity=std::function<VelocityJet(const Vector2D&)>;
struct JetResult
{
    Result trajectory;
    // Absent at physical boundary hits: differentiating the interior RK map
    // does not differentiate the hit time or prescribed boundary history.
    std::optional<Jacobian> derivative;
    int acceptedSubsteps=0;
};

namespace Detail
{
inline void require(bool condition, const char* reason)
{
    if (!condition) throw std::runtime_error(std::string("[LevelSetBoundaryTrace] ") + reason);
}
inline bool finite(Vector2D p) { return std::isfinite(p.x) && std::isfinite(p.y); }
inline bool finite(const Jacobian& a)
{ return std::isfinite(a.xx) && std::isfinite(a.xy) && std::isfinite(a.yx) && std::isfinite(a.yy); }
template<class Matrix> inline Matrix multiply(const Matrix& a,const Matrix& b)
{ return {a.xx*b.xx+a.xy*b.yx,a.xx*b.xy+a.xy*b.yy,
          a.yx*b.xx+a.yy*b.yx,a.yx*b.xy+a.yy*b.yy}; }
template<class Scalar,class Matrix> inline Matrix identityMinus(Scalar factor,const Matrix& a)
{ return {1-factor*a.xx,-factor*a.xy,-factor*a.yx,1-factor*a.yy}; }

// Shared algebra for the actual double tracer and outward interval evaluation.
// Preserve these operation orders; this is not an alternate integrator.
template<class Point,class Scalar> inline Point rk2Midpoint(Point point,Point speed,Scalar duration)
{return point-speed*0.5*duration;}
template<class Point,class Scalar> inline Point rk2Target(Point point,Point speed,Scalar duration)
{return point-speed*duration;}
template<class Scalar,class Matrix> inline Matrix rk2Jacobian(Scalar duration,const Matrix& start,const Matrix& middle)
{return identityMinus(duration,multiply(middle,identityMinus(0.5*duration,start)));}

template<bool Differentiate>
class Tracer
{
    using Sample=std::conditional_t<Differentiate,VelocityJet,Vector2D>;
    using Callback=std::conditional_t<Differentiate,DifferentiableVelocity,Velocity>;
    using Output=std::conditional_t<Differentiate,JetResult,Result>;
    static Vector2D value(const Sample& sample)
    {
        if constexpr (Differentiate) return sample.value;
        else return sample;
    }
    static Output completed(Result result,const Jacobian& derivative,int substeps)
    {
        if constexpr (Differentiate)
            return {result,result.boundaryHit?std::optional<Jacobian>{}:derivative,substeps};
        else return result;
    }
public:
    Tracer(const Domain& domain, double dt, const Callback& current, const Callback& previous)
        : domain_(domain), dt_(dt), current_(current), previous_(previous)
    {
        require(finite(domain.minCoord) && finite(domain.maxCoord) && finite(domain.cellSize)
            && domain.minCoord.x < domain.maxCoord.x && domain.minCoord.y < domain.maxCoord.y
            && domain.cellSize.x > 0.0 && domain.cellSize.y > 0.0,
            "invalid finite rectangular domain or cell size");
        const Vector2D extent = domain.maxCoord - domain.minCoord;
        require(finite(extent) && domain.cellSize.x <= extent.x && domain.cellSize.y <= extent.y,
            "cell size exceeds or cannot resolve the domain");
        require(std::isfinite(dt_) && dt_ >= 0.0 && current_ && previous_,
            "finite nonnegative dt and both velocity callbacks are required");
        require(!domain.radialAxisSymmetry || domain.minCoord.x == 0.0,
            "radial symmetry requires the actual r=0 west boundary");
        h_ = std::min(domain.cellSize.x, domain.cellSize.y);
        const double scale = std::max({h_, std::abs(domain.minCoord.x), std::abs(domain.minCoord.y),
            std::abs(domain.maxCoord.x), std::abs(domain.maxCoord.y)});
        coordinateTolerance_ = 128.0 * std::numeric_limits<double>::epsilon() * scale;
        timeTolerance_ = 128.0 * std::numeric_limits<double>::epsilon() * dt_;
        require(std::isfinite(coordinateTolerance_) && coordinateTolerance_ > 0.0
            && coordinateTolerance_ <= 1e-8 * h_, "coordinate scale is unresolved relative to the cell size");
    }

    Output run(Vector2D point) const
    {
        require(inside(point), "initial point is outside the closed physical domain");
        checkActive(point,point);
        Jacobian derivative=Jacobian::identity();
        if (dt_ == 0.0) return completed({point, 0.0, std::nullopt, false},derivative,0);
        double elapsed = 0.0;
        for (int step = 0; step < 128; ++step)
        {
            const Sample startVelocity = velocity(point, elapsed);
            if (const auto side = outwardBacktraceSide(point, value(startVelocity)))
                return completed({point, elapsed, side, true},derivative,step);
            const double remaining = dt_ - elapsed;
            require(remaining > 0.0 && std::isfinite(remaining), "unresolved remaining trace time");
            double duration = remaining;
            Vector2D half;
            bool safeHalf = false;
            for (int reduction = 0; reduction < 128; ++reduction)
            {
                // Original arithmetic order when elapsed=0 and duration=dt.
                half = rk2Midpoint(point,value(startVelocity),duration);
                checkRadius(half);
                if (inside(half)) { safeHalf = true; break; }
                const double halfDuration=duration*0.5;
                const double spread=std::max(std::abs(half.x-point.x),std::abs(half.y-point.y));
                if (halfDuration<=timeTolerance_ && spread<=coordinateTolerance_)
                {
                    // This stage is NOT sampled or clipped. At the same
                    // resolution used by locateEvent, intersect the physical
                    // wall and validate its actual boundary velocity/coverage.
                    // Event time may round to elapsed; no zero-time interior
                    // substep or interior Jacobian is accepted by this branch.
                    auto hit=intersectEventBracket(point,elapsed,0.0,halfDuration,point,half);
                    hit.roundoffMidpointEvent=true;
                    return completed(hit,derivative,step);
                }
                duration *= 0.5;
                require(duration > 0.0 && elapsed + duration > elapsed,
                    "cannot shorten an out-of-domain midpoint to a resolved substep");
            }
            require(safeHalf, "midpoint remained outside after bounded substep shortening");
            Jacobian localDerivative;
            const Vector2D target = trial(point, elapsed, duration, startVelocity,&localDerivative);
            const double endTime = duration == remaining ? dt_ : elapsed + duration;
            require(endTime > elapsed && endTime <= dt_, "substep time did not advance within the old/new interval");
            if (!inside(target)) return completed(locateEvent(point, elapsed, duration, startVelocity, target),derivative,step);
            if (touchesNewWall(point, target)) return completed(event(target, endTime),derivative,step);
            if constexpr (Differentiate) {
                derivative=multiply(localDerivative,derivative);
                require(finite(derivative),"nonfinite composite numerical-map Jacobian");
            }
            point = target;
            elapsed = endTime;
            if (elapsed == dt_) return completed({point, dt_, std::nullopt, false},derivative,step+1);
        }
        throw std::runtime_error("[LevelSetBoundaryTrace] bounded RK2 substeps did not finish or resolve a boundary event");
    }

private:
    const Domain& domain_;
    double dt_, h_ = 0.0, coordinateTolerance_ = 0.0, timeTolerance_ = 0.0;
    const Callback& current_;
    const Callback& previous_;

    bool inside(Vector2D point) const
    {
        return finite(point) && point.x >= domain_.minCoord.x && point.x <= domain_.maxCoord.x
            && point.y >= domain_.minCoord.y && point.y <= domain_.maxCoord.y;
    }
    void checkRadius(Vector2D point) const
    {
        require(finite(point), "nonfinite numerical trajectory point");
        require(!domain_.radialAxisSymmetry || point.x >= 0.0,
            "radial backtrace reached negative r; clipping or reflection is forbidden");
    }
    void checkActive(Vector2D a,Vector2D b) const
    {
        require(!domain_.activeSegment || domain_.activeSegment(a,b),
            "numerical characteristic crosses outside the selected fluid domain");
    }
    Sample sampled(const Callback& callback, Vector2D point) const
    {
        require(inside(point), "velocity sampling outside the physical domain is forbidden");
        checkActive(point,point);
        const Sample sample = callback(point);
        const Vector2D speed=value(sample);
        require(finite(speed), "nonfinite sampled velocity");
        if constexpr (Differentiate) require(finite(sample.derivative),"nonfinite sampled velocity Jacobian");
        const double travel = std::hypot(speed.x, speed.y) * dt_;
        require(std::isfinite(travel) && travel <= 0.5 * h_, "sampled speed times full dt exceeds half a cell");
        require(!domain_.radialAxisSymmetry || point.x != 0.0 || speed.x == 0.0,
            "sampled radial velocity on the symmetry axis must be exactly zero");
        if constexpr (Differentiate)
            require(!domain_.radialAxisSymmetry || point.x!=0.0 || sample.derivative.xy==0.0,
                "radial velocity tangential derivative on the axis must be exactly zero");
        return sample;
    }
    Sample velocity(Vector2D point, double elapsed) const
    {
        require(std::isfinite(elapsed) && elapsed >= 0.0 && elapsed <= dt_, "velocity time is outside the old/new interval");
        if (elapsed == 0.0) return sampled(current_, point);
        if (elapsed == dt_) return sampled(previous_, point);
        const Sample current = sampled(current_, point), previous = sampled(previous_, point);
        // Preserve the old full-step midpoint average exactly. For a shortened
        // step starting a seconds back, its midpoint is newTime-a-duration/2,
        // NOT the old full-step midpoint or the complete old time.
        const Vector2D result = elapsed == 0.5 * dt_ ? (value(current) + value(previous)) * 0.5 :
            value(current) * (1.0 - elapsed / dt_) + value(previous) * (elapsed / dt_);
        require(finite(result), "interpolated velocity overflow");
        const double travel = std::hypot(result.x, result.y) * dt_;
        require(std::isfinite(travel) && travel <= 0.5 * h_, "interpolated speed times full dt exceeds half a cell");
        if constexpr (Differentiate) {
            const auto mix=[&](double a,double b) {
                return elapsed==0.5*dt_?(a+b)*0.5:a*(1.0-elapsed/dt_)+b*(elapsed/dt_);
            };
            const auto& a=current.derivative;const auto& b=previous.derivative;
            Jacobian derivative{mix(a.xx,b.xx),mix(a.xy,b.xy),mix(a.yx,b.yx),mix(a.yy,b.yy)};
            require(finite(derivative),"interpolated velocity Jacobian overflow");
            return {result,derivative};
        } else return result;
    }
    Vector2D trial(Vector2D point, double elapsed, double duration,const Sample& startVelocity,
        Jacobian* derivative=nullptr) const
    {
        const Vector2D half = rk2Midpoint(point,value(startVelocity),duration);
        checkRadius(half);
        require(inside(half), "event refinement produced an unsampleable midpoint");
        checkActive(point,half);
        const Sample middleVelocity = velocity(half, elapsed + 0.5 * duration);
        const Vector2D target = rk2Target(point,value(middleVelocity),duration);
        if constexpr (Differentiate) if (derivative) {
            *derivative=rk2Jacobian(duration,startVelocity.derivative,middleVelocity.derivative);
            require(finite(*derivative),"nonfinite numerical midpoint-map Jacobian");
        }
        checkRadius(target);
        if (inside(target)) checkActive(point,target);
        return target;
    }
    std::optional<Boundary> outwardBacktraceSide(Vector2D point, Vector2D value) const
    {
        std::optional<Boundary> result;
        const auto append = [&](bool hit, Boundary side)
        {
            if (!hit) return;
            require(!result, "simultaneous backward crossings of two walls are unresolved");
            result = side;
        };
        append(point.x == domain_.minCoord.x && value.x > 0.0, Boundary::west);
        append(point.y == domain_.minCoord.y && value.y > 0.0, Boundary::south);
        append(point.x == domain_.maxCoord.x && value.x < 0.0, Boundary::east);
        append(point.y == domain_.maxCoord.y && value.y < 0.0, Boundary::north);
        return result;
    }
    bool touchesNewWall(Vector2D from, Vector2D target) const
    {
        return (target.x == domain_.minCoord.x && from.x > domain_.minCoord.x)
            || (target.x == domain_.maxCoord.x && from.x < domain_.maxCoord.x)
            || (target.y == domain_.minCoord.y && from.y > domain_.minCoord.y)
            || (target.y == domain_.maxCoord.y && from.y < domain_.maxCoord.y);
    }
    Result event(Vector2D point, double elapsed, std::optional<Boundary> expected = std::nullopt) const
    {
        require(inside(point), "boundary event is outside the closed domain");
        const auto side = outwardBacktraceSide(point, value(velocity(point, elapsed)));
        require(side.has_value() && (!expected || side == expected),
            "boundary event is tangential, non-inflow, or inconsistent with its crossed wall");
        return {point, elapsed, side, true};
    }
    Result locateEvent(Vector2D point, double elapsed, double duration,
        const Sample& startVelocity, Vector2D outside) const
    {
        double lo = 0.0, hi = duration;
        Vector2D insidePoint = point, outsidePoint = outside;
        bool resolved = false;
        for (int iteration = 0; iteration < 96; ++iteration)
        {
            const double spread = std::max(std::abs(outsidePoint.x - insidePoint.x),
                std::abs(outsidePoint.y - insidePoint.y));
            if (hi-lo <= timeTolerance_ && spread <= coordinateTolerance_)
            { resolved = true; break; }
            const double middle = lo + 0.5 * (hi-lo);
            require(middle > lo && middle < hi, "boundary-event time bracket stagnated above roundoff resolution");
            const Vector2D candidate = trial(point, elapsed, middle, startVelocity);
            if (inside(candidate)) { lo = middle; insidePoint = candidate; }
            else { hi = middle; outsidePoint = candidate; }
        }
        require(resolved, "boundary-event isolation exceeded its bounded refinement");
        return intersectEventBracket(point,elapsed,lo,hi,insidePoint,outsidePoint);
    }
    Result intersectEventBracket(Vector2D point,double elapsed,double lo,double hi,
        Vector2D insidePoint,Vector2D outsidePoint) const
    {
        std::optional<Boundary> side;
        double boundaryCoordinate = 0.0;
        bool x = false;
        const auto crossing = [&](bool crossed, Boundary which, bool radial, double coordinate)
        {
            if (!crossed) return;
            require(!side, "multiple wall crossings cannot be separated at the event resolution");
            side = which; x = radial; boundaryCoordinate = coordinate;
        };
        crossing(outsidePoint.x < domain_.minCoord.x, Boundary::west, true, domain_.minCoord.x);
        crossing(outsidePoint.y < domain_.minCoord.y, Boundary::south, false, domain_.minCoord.y);
        crossing(outsidePoint.x > domain_.maxCoord.x, Boundary::east, true, domain_.maxCoord.x);
        crossing(outsidePoint.y > domain_.maxCoord.y, Boundary::north, false, domain_.maxCoord.y);
        require(side.has_value(), "event isolation lost its outside bracket");
        const double a = x ? insidePoint.x : insidePoint.y;
        const double b = x ? outsidePoint.x : outsidePoint.y;
        const double fraction = (boundaryCoordinate-a)/(b-a);
        require(std::isfinite(fraction) && fraction >= 0.0 && fraction <= 1.0,
            "unrepresentable boundary intersection in the resolved event bracket");
        // Geometric intersection of the final roundoff-sized bracket, not a
        // clipped RK stage. No velocity is evaluated at outsidePoint.
        Vector2D hit = insidePoint + (outsidePoint-insidePoint) * fraction;
        if (x) hit.x = boundaryCoordinate; else hit.y = boundaryCoordinate;
        const double hitTime = elapsed + (lo + (hi-lo)*fraction);
        require(std::isfinite(hitTime) && hitTime >= elapsed && hitTime <= dt_, "invalid boundary event time");
        checkActive(point,hit);
        return event(hit, hitTime, side);
    }
};
}

// Successful no-hit full steps retain the original RK2 expression. Shortening
// occurs only to make a midpoint sample legal; every accepted substep consumes
// its actual elapsed time. All new/old samples obey dt*|u| <= min(hx,hy)/2.
// Exceptions leave caller data untouched (callbacks themselves must be pure).
inline Result trace(const Domain& domain, Vector2D point, double dt,
    const Velocity& newVelocity, const Velocity& oldVelocity)
{
    return Detail::Tracer<false>(domain, dt, newVelocity, oldVelocity).run(point);
}

// Same numerical map, sampling order, physical events and half-cell bound as
// trace(). Differentiate every accepted RK2 substep and compose in trajectory
// order. This is the piecewise derivative on the selected interpolation and
// subdivision branch, not a derivative across a branch change. Callbacks must
// be fixed, repeatable fields; derivatives at a non-smooth seam require an
// explicit caller-side convention. No spatial finite differences are used.
inline JetResult traceWithJacobian(const Domain& domain,Vector2D point,double dt,
    const DifferentiableVelocity& newVelocity,const DifferentiableVelocity& oldVelocity)
{
    return Detail::Tracer<true>(domain,dt,newVelocity,oldVelocity).run(point);
}
}
