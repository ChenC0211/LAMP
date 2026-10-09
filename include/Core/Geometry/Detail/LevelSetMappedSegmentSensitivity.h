#pragma once

#include "Core/Geometry/Detail/LevelSetMappedSegment.h"
#include <type_traits>

namespace Quantum::Geometry::Detail
{
namespace
{
namespace MappedHermiteSegment
{
// Terminal differential evaluation at the SAME retained rational root. Root
// identity/deduplication/endpoint decisions never use these floating values.
inline CellChartCoordinate variationCoordinate(const Coordinate &coordinate)
{
    const auto value = ratio(coordinate.numeratorBegin, coordinate.denominator);
    require(compare(value, constant(0)) >= 0 && compare(value, constant(1)) <= 0,
            "variation root lies outside its native cell");
    const bool upper = compare(value, constant(.5)) > 0;
    const auto distance = upper ? subtract(constant(1), value) : value;
    const double offset = static_cast<double>(approximate(distance));
    require(std::isfinite(offset) && offset >= 0 && offset <= .5 &&
                (compare(distance, constant(0)) == 0 || offset > 0),
            "variation root offset is not representable");
    return {offset, upper};
}

template<class Scalar>
auto valueVariation(const BasicQuery<BasicNativeCellLayout<Scalar>>& query,
    const Scalar& scalar,bool fromEnd,long double offset,CornerSectorVariation sectors)
{
    require(std::isfinite(offset) && offset>=0 && offset<=1 && static_cast<long double>(static_cast<double>(offset))==offset,
        "unrepresentable native trace sample parameter");
    Vector2D a=query.inputBegin,b=query.inputEnd;
    const bool reversed=std::make_pair(b.x,b.y)<std::make_pair(a.x,a.y);
    if (reversed) std::swap(a,b);
    auto parameter=constant(static_cast<double>(offset));
    if (fromEnd!=reversed) parameter=subtract(constant(1),parameter);
    const int cell=query.cellAtCanonicalParameter(parameter).native;
    const auto box=scalar.leafBounds(cell);
    const auto local=localCoordinates(box[0],box[1],a,b,parameter,parameter);
    constexpr bool quartic=Scalar::dofCount==25;
    using Variation=std::conditional_t<quartic,LevelSetMappedSegment::QuarticValueVariation,
        LevelSetMappedSegment::CubicValueVariation>;
    const long double hx=static_cast<long double>(box[1].x)-box[0].x,hy=static_cast<long double>(box[1].y)-box[0].y;
    const CellChartPoint point{variationCoordinate(local[0]),variationCoordinate(local[1])};
    std::array<long double,Scalar::dofCount> coefficients{};
    if constexpr (quartic) {
      const auto basis=[](const Coordinate& coordinate) {
        const auto exact=LevelSetHermiteScalar::QuarticCell::rationalBasisJet(coordinate.numeratorBegin,coordinate.denominator);
        std::array<long double,5> values{};
        for (int k=0;k<5;++k) values[k]=approximate(ratio(exact.value[k],exact.denominator));
        return values;
    };
    const auto bx=basis(local[0]),by=basis(local[1]);
    for (int y=0;y<5;++y) for (int x=0;x<5;++x)
        coefficients[5*y+x]=bx[x]*by[y]*(x==2 || x==3?hx:1.L)*(y==2 || y==3?hy:1.L);
    } else {
        const LevelSetHermiteScalar::Detail::Basis bx(point.x),by(point.y);
        constexpr int xs[4]={0,1,1,0},ys[4]={0,0,1,1};
        for (int c=0;c<4;++c) {
            const int x=xs[c],y=ys[c];
            coefficients[4*c]=bx.value[x]*by.value[y];
            coefficients[4*c+1]=hx*bx.value[x+2]*by.value[y];
            coefficients[4*c+2]=hy*bx.value[x]*by.value[y+2];
            coefficients[4*c+3]=hx*hy*bx.value[x+2]*by.value[y+2];
        }
    }
    return Variation{scalar.evaluateCellChart(cell,point).value,
        scalar.pullbackCellSensitivity(cell,coefficients,sectors)};
}

template<class Scalar>
auto negativePhaseVariation(
    const BasicQuery<BasicNativeCellLayout<Scalar>> &query, const Scalar &scalar,
    CornerSectorVariation sectors)
{
    constexpr bool quartic = Scalar::dofCount == 25;
    using Variation = std::conditional_t<quartic,
        LevelSetMappedSegment::QuarticNegativePhaseVariation, LevelSetMappedSegment::NegativePhaseVariation>;
    require(!query.hasMultipleRootCharts,"negative-phase variation requires transverse roots; exact multiple-root chart retained");
    std::vector<Variation> result;
    Vector2D a = query.inputBegin, b = query.inputEnd;
    const bool reversed = std::make_pair(b.x, b.y) < std::make_pair(a.x, a.y);
    if (reversed)
        std::swap(a, b);
    const long double tx = static_cast<long double>(b.x) - a.x, ty = static_cast<long double>(b.y) - a.y;
    for (std::size_t first = 0; first < query.roots.size();)
    {
        const auto &root = query.roots[first];
        std::size_t next = first + 1;
        while (next < query.roots.size() && compare(query.roots[next].parameter, root.parameter) == 0)
            ++next;
        // All equal-parameter closed-cell incidences describe ONE geometric
        // root. Prefer a positive-length trace owner, not a point-only cell.
        const RootRecord *owner = &root;
        for (std::size_t k = first; k < next; ++k)
            if (!query.roots[k].info.pointOnly)
            {
                owner = &query.roots[k];
                break;
            }
        const int cell = owner->info.nativeCellIndex;
        const auto box = scalar.leafBounds(cell);
        const auto local = localCoordinates(box[0], box[1], a, b, root.parameter, root.parameter);
        const CellChartPoint point{variationCoordinate(local[0]), variationCoordinate(local[1])};
        const auto jet = scalar.evaluateCellChart(cell, point);
        const long double derivative = tx * jet.dx + ty * jet.dy;
        const long double floor = 512.0L * std::numeric_limits<double>::epsilon() *
                                  scalar.derivativeResolutionScale(cell) * std::hypot(tx, ty);
        require(std::isfinite(derivative) && std::isfinite(floor) && std::abs(derivative) > floor,
                "negative-phase variation requires a resolved transverse root");
        const bool endpoint =
            compare(root.parameter, constant(0)) == 0 || compare(root.parameter, constant(1)) == 0;
        const long double weight = -(endpoint ? .5L : 1.0L) / std::abs(derivative);
        const long double hx = static_cast<long double>(box[1].x) - box[0].x;
        const long double hy = static_cast<long double>(box[1].y) - box[0].y;
        std::array<long double, Scalar::dofCount> coefficients{};
        if constexpr (quartic)
        {
            // Same quartic basis kernel as the scalar. Retain the rational
            // root until terminal basis evaluation; do not recover it from
            // a rounded physical point or discard the midpoint columns.
            const auto basis = [](const Coordinate &coordinate) {
                const auto exact = LevelSetHermiteScalar::QuarticCell::rationalBasisJet(
                    coordinate.numeratorBegin, coordinate.denominator);
                std::array<long double, 5> values{};
                for (std::size_t k = 0; k < values.size(); ++k)
                    values[k] = approximate(ratio(exact.value[k], exact.denominator));
                return values;
            };
            const auto bx = basis(local[0]), by = basis(local[1]);
            for (int y = 0; y < 5; ++y)
                for (int x = 0; x < 5; ++x)
                    coefficients[5*y+x] = weight * bx[x] * by[y]
                        * (x == 2 || x == 3 ? hx : 1.0L) * (y == 2 || y == 3 ? hy : 1.0L);
        }
        else
        {
            // Preserve the original cubic operation order and corner order.
            const LevelSetHermiteScalar::Detail::Basis bx(point.x), by(point.y);
            constexpr int xs[4] = {0, 1, 1, 0}, ys[4] = {0, 0, 1, 1};
            for (int corner = 0; corner < 4; ++corner)
            {
                const int x = xs[corner], y = ys[corner];
                coefficients[4 * corner] = weight * bx.value[x] * by.value[y];
                coefficients[4 * corner + 1] = weight * hx * bx.value[x + 2] * by.value[y];
                coefficients[4 * corner + 2] = weight * hy * bx.value[x] * by.value[y + 2];
                coefficients[4 * corner + 3] = weight * hx * hy * bx.value[x + 2] * by.value[y + 2];
            }
        }
        const bool fromEnd = compare(root.parameter, constant(.5)) > 0;
        const auto offset = fromEnd ? subtract(constant(1), root.parameter) : root.parameter;
        result.push_back({fromEnd != reversed, approximate(offset), endpoint,
                          scalar.pullbackCellSensitivity(cell, coefficients, sectors)});
        first = next;
    }
    if (reversed)
        std::reverse(result.begin(), result.end());
    return result;
}
} // namespace MappedHermiteSegment
} // namespace
} // namespace Quantum::Geometry::Detail
