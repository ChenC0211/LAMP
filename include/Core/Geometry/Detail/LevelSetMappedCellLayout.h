#pragma once

#include "Core/Geometry/Detail/LevelSetHermiteGeometry.h"
#include "Core/Geometry/Detail/LevelSetHermiteQuarticGridData.h"
#include "Core/Geometry/Detail/LevelSetHermiteQuarticLeafData.h"
#include "Core/Geometry/Detail/LevelSetLeafHermiteScalar.h"
#include "Core/Geometry/Detail/LevelSetMappedCoordinates.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <ostream>
#include <tuple>
#include <utility>
#include <vector>

namespace Quantum::Geometry::Detail
{
namespace
{
namespace MappedHermiteSegment
{
using Polynomial = LevelSetHermiteScalar::Snapshot::ExpandedPolynomial;
using Coordinate = LevelSetHermiteScalar::Snapshot::RationalSegmentCoordinate;

// Logical Cartesian cells and native leaves are distinct identities. Native
// leaves never masquerade as an nx-by-one logical grid.
struct CellAddress
{
    int i = -1, j = -1, native = -1;
    bool operator==(const CellAddress &other) const
    {
        return i == other.i && j == other.j && native == other.native;
    }
};
inline std::ostream &operator<<(std::ostream &stream, const CellAddress &cell)
{
    if (cell.native >= 0)
        return stream << "native=" << cell.native;
    return stream << cell.i << ',' << cell.j;
}

// Original exact physical-to-local composition, shared by both layouts.
// Positive common scaling changes neither signs nor roots.
inline std::array<Coordinate, 2> localCoordinates(Vector2D lower, Vector2D upper, Vector2D anchorBegin,
                                                  Vector2D anchorEnd, const Ratio &begin, const Ratio &end)
{
    const bool common = sameDenominator(begin, end);
    const auto denominator = common ? begin.denominator : Exact::multiply(begin.denominator, end.denominator);
    const auto first = common ? begin.numerator : Exact::multiply(begin.numerator, end.denominator);
    const auto last = common ? end.numerator : Exact::multiply(end.numerator, begin.denominator);
    std::array<Coordinate, 2> local;
    for (int axis = 0; axis < 2; ++axis)
    {
        const double origin = axis == 0 ? lower.x : lower.y;
        const double a = axis == 0 ? anchorBegin.x : anchorBegin.y, b = axis == 0 ? anchorEnd.x : anchorEnd.y;
        const auto base = Exact::multiply(difference(a, origin), denominator), slope = difference(b, a);
        local[axis] = {Exact::add(base, Exact::multiply(slope, first)),
                       Exact::add(base, Exact::multiply(slope, last)),
                       Exact::multiply(denominator, difference(axis == 0 ? upper.x : upper.y, origin))};
        // The clearing factor is arbitrary and positive. Keep it near
        // unity before taking its third power: a tiny *chart denominator*
        // must not manufacture underflow in an otherwise resolved local
        // polynomial. Scale ALL three expansions by the same exact power
        // of two; no coordinate quotient, coefficient or root is rounded.
        require(Exact::sign(local[axis].denominator) > 0, "nonpositive shared cell width");
        int exponent = 0;
        std::frexp(Exact::estimate(local[axis].denominator), &exponent);
        local[axis].numeratorBegin = Exact::scalePow2(local[axis].numeratorBegin, -exponent);
        local[axis].numeratorEnd = Exact::scalePow2(local[axis].numeratorEnd, -exponent);
        local[axis].denominator = Exact::scalePow2(local[axis].denominator, -exponent);
    }

    return local;
}

template<class Geometry> class BasicUniformCellLayout
{
  public:
    using Source = Geometry;
    using Polynomial = typename Source::ExpandedPolynomial;
    explicit BasicUniformCellLayout(const Source &source) : source_(source)
    {
    }
    std::vector<Ratio> cuts(Vector2D a, Vector2D b)
    {
        a_ = a;
        b_ = b;
        std::vector<Ratio> cuts{constant(0.0), constant(1.0)};
        for (int axis = 0; axis < 2; ++axis)
        {
            const double first = axis == 0 ? a_.x : a_.y, last = axis == 0 ? b_.x : b_.y;
            if (first == last)
                continue;
            const double minimum = std::min(first, last), maximum = std::max(first, last);
            // Shared canonical grid lines are monotone. Locate the first
            // STRICTLY interior line, then retain the original increasing-k
            // cut order and exact ratios. No division-to-cell-index rounding,
            // snapping, root or endpoint-incidence changes are introduced.
            int lo = 1, hi = count(axis);
            while (lo < hi)
            {
                const int middle = lo + (hi - lo) / 2;
                if (coordinate(middle, axis) <= minimum)
                    lo = middle + 1;
                else
                    hi = middle;
            }
            for (int k = lo; k < count(axis); ++k)
            {
                const double line = coordinate(k, axis);
                if (line >= maximum)
                    break;
                cuts.push_back(ratio(difference(line, first), difference(last, first)));
            }
        }

        return cuts;
    }
    int strictSign(const CellAddress &cell) const
    {
        return source_.strictCellSign(cell.i, cell.j);
    }
    int pieceSign(const CellAddress& cell,const Ratio&,const Ratio&) const {return strictSign(cell);}
    Polynomial polynomial(const CellAddress &cell, const Ratio &begin, const Ratio &end) const
    {
        const auto local =
            localCoordinates({coordinate(cell.i, 0), coordinate(cell.j, 1)},
                             {coordinate(cell.i + 1, 0), coordinate(cell.j + 1, 1)}, a_, b_, begin, end);
        return source_.clearedCellSegment(cell.i, cell.j, local[0], local[1]);
    }
    bool sameTrace(const Polynomial &p, const Polynomial &q, int pivot) const
    {
        for (std::size_t k = 0; k < p.size(); ++k)
            if (Exact::sign(
                    Exact::subtract(Exact::multiply(p[k], q[pivot]), Exact::multiply(q[k], p[pivot]))) != 0)
                return false;
        return true;
    }
    int count(int axis) const
    {
        return axis == 0 ? source_.nx() : source_.ny();
    }

    double coordinate(int k, int axis) const
    {
        const auto lower = source_.minCoord(), upper = source_.maxCoord(), h = source_.cellSize();
        if (k == count(axis))
            return axis == 0 ? upper.x : upper.y;
        return (axis == 0 ? lower.x : lower.y) + k * (axis == 0 ? h.x : h.y);
    }

    std::vector<int> incidentIndicesAtInputPoint(Vector2D point, int axis) const
    {
        // Input endpoints are binary64 coordinates already validated against
        // the rectangle, not approximations to internal rational/root points.
        // Ordered double comparisons therefore give EXACT endpoint incidence.
        const double value = axis == 0 ? point.x : point.y;
        int lo = 0, hi = count(axis);
        while (lo < hi)
        {
            const int middle = lo + (hi - lo) / 2;
            if (coordinate(middle, axis) < value)
                lo = middle + 1;
            else
                hi = middle;
        }
        std::vector<int> result;
        if (lo > 0)
            result.push_back(lo - 1);
        if (lo < count(axis) && coordinate(lo, axis) == value)
            result.push_back(lo);
        return result;
    }

    std::vector<CellAddress> incidentCellsAtInputPoint(Vector2D point) const
    {
        const auto xs = incidentIndicesAtInputPoint(point, 0), ys = incidentIndicesAtInputPoint(point, 1);
        std::vector<CellAddress> result;
        for (int j : ys)
            for (int i : xs)
                if (source_.cellActive(i, j))
                    result.push_back({i, j});
        return result;
    }

    Ratio physical(const Ratio &parameter, int axis) const
    {
        const double a = axis == 0 ? a_.x : a_.y, b = axis == 0 ? b_.x : b_.y;
        return ratio(Exact::add(Exact::multiply(number(a), parameter.denominator),
                                Exact::multiply(difference(b, a), parameter.numerator)),
                     parameter.denominator);
    }

    std::vector<int> incidentIndices(const Ratio &parameter, int axis) const
    {
        const auto p = physical(parameter, axis);
        int lower = 0, upper = count(axis);
        while (lower < upper)
        {
            const int middle = lower + (upper - lower) / 2;
            if (compare(constant(coordinate(middle, axis)), p) < 0)
                lower = middle + 1;
            else
                upper = middle;
        }
        require(compare(p, constant(coordinate(0, axis))) >= 0 &&
                    compare(p, constant(coordinate(count(axis), axis))) <= 0,
                "exact segment parameter left the snapshot rectangle");
        std::vector<int> result;
        if (lower > 0)
            result.push_back(lower - 1);
        if (lower < count(axis) && compare(p, constant(coordinate(lower, axis))) == 0)
            result.push_back(lower);
        return result;
    }

    std::vector<CellAddress> incidentCells(const Ratio &parameter) const
    {
        const auto xs = incidentIndices(parameter, 0), ys = incidentIndices(parameter, 1);
        std::vector<CellAddress> result;
        for (int j : ys)
            for (int i : xs)
                if (source_.cellActive(i, j))
                    result.push_back({i, j});
        return result;
    }

  private:
    const Source &source_;
    Vector2D a_, b_;
};

using UniformCellLayout = BasicUniformCellLayout<LevelSetHermiteGeometry::Snapshot>;
using QuarticGeometry = LevelSetHermiteGeometry::BasicSnapshot<LevelSetHermiteScalar::QuarticGridData>;
using QuarticCellLayout = BasicUniformCellLayout<QuarticGeometry>;

// Only geometric addressing changes. Segment clipping uses the same retained
// rational parameter as roots and phase ranges, including point-only contacts.
// The shared query below retains the original isolator and acceptance checks.
template<class Scalar> class BasicNativeCellLayout
{
    struct Clip
    {
        int native;
        Ratio begin, end;
    };

  public:
    using Source = Scalar;
    using Polynomial = decltype(std::declval<typename Source::Polynomial>().numerator);
    explicit BasicNativeCellLayout(const Source &source) : source_(source)
    {
    }
    // Internal chart for a SAME-source normal Taylor coefficient. The query
    // keeps the original physical anchors; all local parameters, including
    // dyadic/refined subcharts, compose through this exact retained chart.
    BasicNativeCellLayout(const Source& source,CellAddress owner,Chart chart,int axis,int order)
        : source_(source),restrictedOwner_(owner.native),restrictedChart_(std::move(chart)),
          normalAxis_(axis),normalOrder_(order)
    {
        require(owner.native>=0 && source.cellActive(owner.native),"normal trace owner is inactive");
        require(axis==0 || axis==1,"invalid normal trace axis");
        require(order>=1 && order<=4,"invalid normal Taylor order");
    }
    std::vector<Ratio> cuts(Vector2D a, Vector2D b)
    {
        a_ = a;
        b_ = b;
        clips_.clear();
        std::vector<Ratio> result{constant(0.0), constant(1.0)};
        if(restrictedOwner_>=0)return result;
        const auto candidates = source_.leafCandidates({std::min(a.x, b.x), std::min(a.y, b.y)},
                                                       {std::max(a.x, b.x), std::max(a.y, b.y)});
        for (int index : candidates)
        {
            const auto bounds = source_.leafBounds(index);
            Ratio lo = constant(0.0), hi = constant(1.0);
            bool intersects = true;
            for (int axis = 0; axis < 2; ++axis)
            {
                const double first = axis == 0 ? a.x : a.y, last = axis == 0 ? b.x : b.y;
                const double lower = axis == 0 ? bounds[0].x : bounds[0].y;
                const double upper = axis == 0 ? bounds[1].x : bounds[1].y;
                if (first == last)
                {
                    if (first < lower || first > upper)
                    {
                        intersects = false;
                        break;
                    }
                    continue;
                }
                auto begin = ratio(difference(lower, first), difference(last, first));
                auto end = ratio(difference(upper, first), difference(last, first));
                if (compare(begin, end) > 0)
                    std::swap(begin, end);
                if (compare(begin, lo) > 0)
                    lo = std::move(begin);
                if (compare(end, hi) < 0)
                    hi = std::move(end);
                if (compare(lo, hi) > 0)
                {
                    intersects = false;
                    break;
                }
            }
            if (!intersects)
                continue;
            clips_.push_back({index, lo, hi});
            result.push_back(std::move(lo));
            result.push_back(std::move(hi));
        }
        require(!clips_.empty(), "segment has no native leaf coverage");
        return result;
    }
    std::vector<CellAddress> incidentCellsAtInputPoint(Vector2D point) const
    {
        if(restrictedOwner_>=0)return {{-1,-1,restrictedOwner_}};
        std::vector<CellAddress> result;
        for (const auto &clip : clips_)
        {
            if (!source_.cellActive(clip.native))
                continue;
            const auto box = source_.leafBounds(clip.native);
            if (point.x >= box[0].x && point.x <= box[1].x && point.y >= box[0].y && point.y <= box[1].y)
                result.push_back({-1, -1, clip.native});
        }
        sortCells(result);
        return result;
    }
    std::vector<CellAddress> incidentCells(const Ratio &parameter) const
    {
        if(restrictedOwner_>=0)return {{-1,-1,restrictedOwner_}};
        std::vector<CellAddress> result;
        for (const auto &clip : clips_)
            if (source_.cellActive(clip.native) && compare(clip.begin, parameter) <= 0 &&
                compare(parameter, clip.end) <= 0)
                result.push_back({-1, -1, clip.native});
        sortCells(result);
        return result;
    }
    int strictSign(const CellAddress &cell) const
    {
        if(restrictedOwner_>=0)return 0;
        return source_.cellSignExceptVertices(cell.native, 0u);
    }
    int pieceSign(const CellAddress& cell,const Ratio& begin,const Ratio& end) const
    {
        if(restrictedOwner_>=0)return 0;
        const int whole=strictSign(cell);if(whole)return whole;
        if constexpr(std::tuple_size<Polynomial>::value==9) {
            // A short grid-clipped part of a mixed leaf can be wholly pure.
            // Certify its exact-rational local box BEFORE expanding the Q4
            // line polynomial; no high powers of irrelevant tiny tails are
            // needed and no phase/root is discarded by a length threshold.
            const auto box=source_.leafBounds(cell.native);
            const auto local=localCoordinates(box[0],box[1],a_,b_,begin,end);
            std::array<double,2> lo{},hi{};
            for(int axis=0;axis<2;++axis) {
                const auto a=Exact::enclosure(local[axis].numeratorBegin),b=Exact::enclosure(local[axis].numeratorEnd);
                const auto d=Exact::enclosure(local[axis].denominator);
                require(d.lo>0 && std::isfinite(d.hi),"unresolved local phase-box denominator");
                lo[axis]=std::clamp(std::nextafter(std::min(a.lo,b.lo)/d.hi,-std::numeric_limits<double>::infinity()),0.,1.);
                hi[axis]=std::clamp(std::nextafter(std::max(a.hi,b.hi)/d.lo,std::numeric_limits<double>::infinity()),0.,1.);
            }
            return source_.localBoxSign(cell.native,{lo[0],lo[1]},{hi[0],hi[1]});
        }
        return 0;
    }
    Polynomial polynomial(const CellAddress &cell, const Ratio &begin, const Ratio &end) const
    {
        const auto box = source_.leafBounds(cell.native);
        const auto compose=[&](const Ratio& t){
            if(restrictedOwner_<0)return t;
            // between() takes a binary64 local value; this overload retains
            // the full ratio created by nested root/refinement charts.
            const auto width=subtract(restrictedChart_.end,restrictedChart_.begin);
            return add(restrictedChart_.begin,ratio(Exact::multiply(width.numerator,t.numerator),
                Exact::multiply(width.denominator,t.denominator)));
        };
        const auto local = localCoordinates(box[0], box[1], a_, b_, compose(begin), compose(end));
        if constexpr(std::tuple_size<Polynomial>::value==9)
            if(restrictedOwner_>=0)return source_.restrictLocal(cell.native,local[0],local[1],
                normalAxis_==0?normalOrder_:0,normalAxis_==1?normalOrder_:0).numerator;
        return source_.restrictLocal(cell.native, local[0], local[1]).numerator;
    }
    std::array<Vector2D,2> bounds(const CellAddress& cell) const{return source_.leafBounds(cell.native);}
    bool sameTrace(const Polynomial &p, const Polynomial &q, int pivot) const
    {
        // Native hanging constraints introduce positive clearing factors.
        // Their low bits may underflow when two near-unit polynomials are
        // cross-multiplied, even when the identity fits the expansion range.
        // Center each WHOLE polynomial's component exponents first. These
        // exact positive power-of-two factors cancel in every cross-product;
        // no coefficient is discarded or divided and root inputs are unchanged.
        const auto a = centeredTrace(p), b = centeredTrace(q);
        for (std::size_t k = 0; k < p.size(); ++k)
            if (Exact::sign(
                    Exact::subtract(Exact::multiply(a[k], b[pivot]), Exact::multiply(b[k], a[pivot]))) != 0)
                return false;
        return true;
    }

  private:
    const Source &source_;
    Vector2D a_, b_;
    std::vector<Clip> clips_;
    int restrictedOwner_=-1;
    Chart restrictedChart_;
    int normalAxis_=-1,normalOrder_=0;
    static Polynomial centeredTrace(const Polynomial &polynomial)
    {
        Exact::ExponentRange range;
        for (const auto &coefficient : polynomial)range.include(coefficient);
        require(range.lowest != std::numeric_limits<int>::max(), "cannot compare an identically zero trace");
        const int shift = range.centeredShift();
        Polynomial result;
        for (std::size_t k = 0; k < result.size(); ++k)
            result[k] = Exact::scalePow2(polynomial[k], shift);
        return result;
    }
    void sortCells(std::vector<CellAddress> &cells) const
    {
        std::sort(cells.begin(), cells.end(), [&](const CellAddress &a, const CellAddress &b) {
            const auto x = source_.leafBounds(a.native)[0], y = source_.leafBounds(b.native)[0];
            return std::make_tuple(x.y, x.x, a.native) < std::make_tuple(y.y, y.x, b.native);
        });
    }
};
using NativeCellLayout=BasicNativeCellLayout<LevelSetHermiteScalar::LeafSnapshot>;
using NativeQuarticCellLayout=BasicNativeCellLayout<LevelSetHermiteScalar::QuarticLeafData>;
} // namespace MappedHermiteSegment
} // namespace
} // namespace Quantum::Geometry::Detail
