#pragma once

#include "Core/Geometry/Detail/LevelSetMappedSegment.h"
#include "Core/Geometry/LevelSetMappedInterface.h"
#include <array>
#include <cstddef>
#include <stdexcept>

namespace Quantum::Geometry::Detail
{
namespace
{
namespace MappedHermiteInterface
{
// Only cell ownership, coefficients and physical embedding depend on layout.
// Both adapters use the existing scalar restriction and the same integrator.
struct CubicDofLayout
{
    static constexpr std::size_t dofCount=16;
    static std::array<int,2> tensorIndex(std::size_t k)
    {
        constexpr int xs[4]={0,1,1,0},ys[4]={0,0,1,1};
        return {{xs[k/4]+((k%4&1)?2:0),ys[k/4]+((k%4&2)?2:0)}};
    }
};
struct QuarticDofLayout
{
    static constexpr std::size_t dofCount=25;
    static std::array<int,2> tensorIndex(std::size_t k)
    { return {{static_cast<int>(k%5),static_cast<int>(k/5)}}; }
};
template<class Geometry,class DofLayout=CubicDofLayout>
class BasicUniformCurveCell : public DofLayout
{
  public:
    using Source = Geometry;
    using Query = MappedHermiteSegment::BasicQuery<MappedHermiteSegment::BasicUniformCellLayout<Source>>;
    using Polynomial = typename Query::Polynomial;
    BasicUniformCurveCell(const Source &source, const ::LevelSetMappedInterface::Graph &graph)
        : source_(source), cell_{graph.cellI, graph.cellJ, -1}
    {
        ExactExpansion::validateEnvironment();
        if (graph.nativeCellIndex != -1 || graph.cellI < 0 || graph.cellI >= source.nx() || graph.cellJ < 0 ||
            graph.cellJ >= source.ny() || !source.cellActive(graph.cellI, graph.cellJ))
            throw std::invalid_argument(
                "[LevelSetMappedInterface] graph requires a valid active uniform cell");
        nominal_ = source.cellSize();
        const auto low = source.minCoord(), high = source.maxCoord();
        const auto coordinate = [&](int k, int axis) {
            const int count = axis == 0 ? source.nx() : source.ny();
            if (k == count)
                return axis == 0 ? high.x : high.y;
            return (axis == 0 ? low.x : low.y) + k * (axis == 0 ? nominal_.x : nominal_.y);
        };
        bounds_ = {{{coordinate(cell_.i, 0), coordinate(cell_.j, 1)},
                    {coordinate(cell_.i + 1, 0), coordinate(cell_.j + 1, 1)}}};
    }
    MappedHermiteSegment::CellAddress address() const
    {
        return cell_;
    }
    std::array<Vector2D, 2> bounds() const
    {
        return bounds_;
    }
    Vector2D derivativeScale() const
    {
        return nominal_;
    }
    std::array<double, 4> dofs(int x, int y) const
    {
        return source_.mappedCellNodeDofs(cell_.i, cell_.j, x, y);
    }
    std::array<double,DofLayout::dofCount> coefficients() const
    {
        if constexpr (DofLayout::dofCount==25) return source_.mappedCellDofs(cell_.i,cell_.j);
        else
        {
            std::array<double,16> result;
            constexpr int xs[4]={0,1,1,0},ys[4]={0,0,1,1};
            for (int c=0;c<4;++c)
            {
                const auto node=dofs(xs[c],ys[c]);
                for (int d=0;d<4;++d) result[4*c+d]=node[d];
            }
            return result;
        }
    }
    Polynomial polynomial(const MappedHermiteSegment::Coordinate &x,
                                                const MappedHermiteSegment::Coordinate &y) const
    {
        return source_.clearedCellSegment(cell_.i, cell_.j, x, y);
    }
    void finalizeSensitivity(::LevelSetMappedInterface::CellArc<DofLayout::dofCount>&) const {}

  private:
    const Source &source_;
    MappedHermiteSegment::CellAddress cell_;
    Vector2D nominal_;
    std::array<Vector2D, 2> bounds_;
};

using UniformCurveCell = BasicUniformCurveCell<LevelSetHermiteGeometry::Snapshot>;
using QuarticCurveCell = BasicUniformCurveCell<MappedHermiteSegment::QuarticGeometry,QuarticDofLayout>;

template<class Scalar,class DofLayout=CubicDofLayout>
class BasicNativeCurveCell : public DofLayout
{
  public:
    using Source = Scalar;
    using Query = MappedHermiteSegment::BasicQuery<MappedHermiteSegment::BasicNativeCellLayout<Source>>;
    using Polynomial = typename Query::Polynomial;
    BasicNativeCurveCell(const Source &source, const ::LevelSetMappedInterface::Graph &graph)
        : source_(source), cell_{-1, -1, graph.nativeCellIndex}
    {
        if (graph.cellI != -1 || graph.cellJ != -1 || graph.nativeCellIndex < 0 ||
            !source.cellActive(graph.nativeCellIndex))
            throw std::invalid_argument(
                "[LevelSetMappedInterface] graph requires a valid active native leaf");
        bounds_ = source.leafBounds(cell_.native);
    }
    MappedHermiteSegment::CellAddress address() const
    {
        return cell_;
    }
    std::array<Vector2D, 2> bounds() const
    {
        return bounds_;
    }
    Vector2D derivativeScale() const
    {
        // This is a terminal sensitivity-column scale, not the root map.
        // Exact restriction uses the retained physical width expansions.
        return {bounds_[1].x - bounds_[0].x, bounds_[1].y - bounds_[0].y};
    }
    std::array<double, 4> dofs(int x, int y) const
    {
        return source_.mappedCellNodeDofs(cell_.native, x, y);
    }
    std::array<double,DofLayout::dofCount> coefficients() const
    {
        if constexpr(DofLayout::dofCount==25) return source_.mappedCellLocalDofs(cell_.native);
        else
        {
            std::array<double,16> result;
            constexpr int xs[4]={0,1,1,0},ys[4]={0,0,1,1};
            for (int c=0;c<4;++c)
            {
                const auto node=dofs(xs[c],ys[c]);
                for (int d=0;d<4;++d) result[4*c+d]=node[d];
            }
            return result;
        }
    }
    Polynomial polynomial(const MappedHermiteSegment::Coordinate &x,
                                                const MappedHermiteSegment::Coordinate &y) const
    {
        return source_.restrictLocal(cell_.native, x, y).numerator;
    }
    void finalizeSensitivity(::LevelSetMappedInterface::CellArc<DofLayout::dofCount>& arc) const
    {
        if constexpr(DofLayout::dofCount==25)
        {
            arc.independentSensitivity=source_.pullbackCellSensitivity(cell_.native,arc.sensitivity);
            arc.requiresNativeDofConstraints=false;
        }
        // Legacy cubic caller still chooses its explicit corner-sector policy
        // before the original 16-column pullback. Do not change that path.
    }

  private:
    const Source &source_;
    MappedHermiteSegment::CellAddress cell_;
    std::array<Vector2D, 2> bounds_;
};
using NativeCurveCell=BasicNativeCurveCell<LevelSetHermiteScalar::LeafSnapshot>;
using NativeQuarticCurveCell=BasicNativeCurveCell<LevelSetHermiteScalar::QuarticLeafData,QuarticDofLayout>;
} // namespace MappedHermiteInterface
} // namespace
} // namespace Quantum::Geometry::Detail
