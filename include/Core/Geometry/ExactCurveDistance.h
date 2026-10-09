#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>
#include "Core/Base/Vector2D.h"

// Shared geometric primitives for level-set initialization. Angles in the
// geometry dictionary are multiples of pi, not radians. No mesh size or
// particular circle/level-set configuration enters these calculations.
namespace ExactCurveDistance
{
constexpr double pi = 3.141592653589793238462643383279502884;

inline void requireFinitePoint(const Vector2D& point)
{
    if (!std::isfinite(point.x) || !std::isfinite(point.y))
        throw std::invalid_argument("[ExactCurveDistance] point coordinates must be finite");
}

inline double pointDistance(const Vector2D& a, const Vector2D& b)
{
    return std::hypot(a.x - b.x, a.y - b.y);
}

class CircularArc
{
public:
    // center/radius, start's direction from center, and signed sweep define
    // the arc. A rounded dictionary start need not lie exactly at radius.
    // Positive sweep is CCW; reversed changes traversal, not the point set.
    // |sweepInPi| <= 2 supports minor/major arcs and one complete circle.
    // A zero radius or zero sweep is a point, not an implicit full circle.
    CircularArc(const Vector2D& center, const Vector2D& start,
                const double radius, const double sweepInPi,
                const bool reversed = false)
        : center_(center), radius_(radius)
    {
        requireFinitePoint(center);
        requireFinitePoint(start);
        if (!std::isfinite(radius) || radius < 0.0 ||
            !std::isfinite(sweepInPi) || std::abs(sweepInPi) > 2.0)
            throw std::invalid_argument("[ExactCurveDistance] arc requires finite radius >= 0 and sweep in [-2,2] pi");
        if (!std::isfinite(center.x - radius) || !std::isfinite(center.x + radius) ||
            !std::isfinite(center.y - radius) || !std::isfinite(center.y + radius))
            throw std::invalid_argument("[ExactCurveDistance] arc extent exceeds finite coordinates");

        const double dx = start.x - center.x;
        const double dy = start.y - center.y;
        if (!std::isfinite(dx) || !std::isfinite(dy) ||
            (radius > 0.0 && dx == 0.0 && dy == 0.0))
            throw std::invalid_argument("[ExactCurveDistance] nonzero-radius arc requires a finite noncentral start direction");
        angle_ = radius > 0.0 ? std::atan2(dy, dx) : 0.0;
        sweep_ = sweepInPi * pi;
        fullCircle_ = std::abs(sweepInPi) == 2.0;
        start_ = atAngle(angle_);
        end_ = fullCircle_ ? start_ : atAngle(angle_ + sweep_);
        if (reversed)
        {
            angle_ += sweep_;
            sweep_ = -sweep_;
            std::swap(start_, end_);
        }

        // Split at horizontal tangencies. The usual half-open segment test
        // then counts shared endpoints once and tangent contacts zero/twice.
        // Each piece has a unique circle/ray intersection on its x branch.
        if (radius_ > 0.0 && sweep_ != 0.0)
        {
            std::vector<double> cuts{0.0, 1.0};
            for (int k = -4; k <= 4; ++k)
            {
                const double t = (0.5 * pi + k * pi - angle_) / sweep_;
                if (t > 0.0 && t < 1.0)
                    cuts.push_back(t);
            }
            std::sort(cuts.begin(), cuts.end());
            for (std::size_t i = 1; i < cuts.size(); ++i)
            {
                const double t0 = cuts[i - 1];
                const double t1 = cuts[i];
                const double midAngle = angle_ + 0.5 * (t0 + t1) * sweep_;
                rayPieces_.push_back({point(t0).y, point(t1).y,
                                      std::cos(midAngle) >= 0.0 ? 1.0 : -1.0});
            }
        }
    }

    Vector2D start() const { return start_; }
    Vector2D end() const { return end_; }

    // Unit tangent from the selected endpoint INTO this oriented arc. This
    // describes the initial primitive, not a prescribed dynamic contact angle.
    Vector2D inwardEndpointTangent(bool atEnd) const
    {
        if (!(radius_ > 0.0) || sweep_ == 0.0)
            throw std::invalid_argument("[ExactCurveDistance] a point arc has no endpoint tangent");
        const Vector2D radial = (atEnd ? end_ : start_) - center_;
        const double magnitude = std::hypot(radial.x, radial.y);
        if (!(magnitude > 0.0) || !std::isfinite(magnitude))
            throw std::invalid_argument("[ExactCurveDistance] unrepresentable endpoint tangent");
        const double sign = (sweep_ > 0.0 ? 1.0 : -1.0) * (atEnd ? -1.0 : 1.0);
        return Vector2D(-sign * radial.y / magnitude, sign * radial.x / magnitude);
    }

    Vector2D point(const double t) const
    {
        if (!std::isfinite(t) || t < 0.0 || t > 1.0)
            throw std::invalid_argument("[ExactCurveDistance] arc parameter must be in [0,1]");
        if (t == 0.0) return start_;
        if (t == 1.0) return end_;
        return atAngle(angle_ + t * sweep_);
    }

    double distance(const Vector2D& query) const
    {
        requireFinitePoint(query);
        if (radius_ == 0.0 || sweep_ == 0.0)
            return pointDistance(query, start_);
        const double dx = query.x - center_.x;
        const double dy = query.y - center_.y;
        const double radial = std::hypot(dx, dy);
        if (radial == 0.0)
            return radius_; // every point on a nondegenerate arc is equidistant
        if (fullCircle_)
            return std::abs(radial - radius_);

        double delta = (sweep_ > 0.0 ? 1.0 : -1.0) *
                       (std::atan2(dy, dx) - angle_);
        delta = std::fmod(delta, 2.0 * pi);
        if (delta < 0.0) delta += 2.0 * pi;
        double result = std::min(pointDistance(query, start_), pointDistance(query, end_));
        if (delta <= std::abs(sweep_))
            result = std::min(result, std::abs(radial - radius_));
        return result;
    }

    // Parity of intersections with the ray (query.x,+infinity) at query.y.
    // As with polygon parity, values exactly on the curve have no unique
    // inside/outside meaning; initializeLevelSets retains its zero safeguard.
    bool rayCrossingParity(const Vector2D& query) const
    {
        requireFinitePoint(query);
        bool parity = false;
        for (const auto& piece : rayPieces_)
        {
            if ((piece.y0 > query.y) == (piece.y1 > query.y))
                continue;
            const double ratio = (query.y - center_.y) / radius_;
            const double root = std::sqrt(std::max(0.0, (1.0 - std::abs(ratio)) *
                                                        (1.0 + std::abs(ratio))));
            const double intersectionX = center_.x + piece.xSign * radius_ * root;
            if (query.x < intersectionX)
                parity = !parity;
        }
        return parity;
    }

private:
    struct RayPiece { double y0, y1, xSign; };
    Vector2D center_, start_, end_;
    double radius_ = 0.0, angle_ = 0.0, sweep_ = 0.0;
    bool fullCircle_ = false;
    std::vector<RayPiece> rayPieces_;

    Vector2D atAngle(const double angle) const
    {
        const double reduced = std::remainder(angle, 2.0 * pi);
        double cosine = std::cos(reduced);
        double sine = std::sin(reduced);
        // Avoid artificial gaps at cardinal endpoints shared by adjacent arcs.
        const double tolerance = 8.0 * std::numeric_limits<double>::epsilon();
        if (std::abs(cosine) < tolerance) cosine = 0.0;
        if (std::abs(sine) < tolerance) sine = 0.0;
        return Vector2D(center_.x + radius_ * cosine, center_.y + radius_ * sine);
    }
};

inline bool segmentRayCrossingParity(const Vector2D& query,
                                     const Vector2D& start, const Vector2D& end)
{
    if ((start.y > query.y) == (end.y > query.y))
        return false;
    // Preserve an exact vertical intersection. Weighted interpolation of two
    // identical x values can round one ulp inward and misclassify a domain
    // wall's one-sided query; no tolerance or coordinate displacement fixes
    // that geometric identity as directly as retaining the original x.
    if (start.x == end.x)
        return query.x < start.x;
    const double t = (query.y - start.y) / (end.y - start.y);
    return query.x < (1.0 - t) * start.x + t * end.x;
}
} // namespace ExactCurveDistance
