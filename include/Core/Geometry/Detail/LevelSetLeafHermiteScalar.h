#pragma once

#include "Core/Geometry/ActiveHermitePatch.h"
#include "Core/Geometry/ActiveHermiteSensitivity.h"
#include "Core/Geometry/Detail/LevelSetHermiteCell.h"
#include "Core/Geometry/Detail/LevelSetHermiteQuarticCell.h"
#include "Core/Geometry/Detail/LevelSetNativeCornerSector.h"
#include "Core/Mesh/CartesianLeafLocator.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <iomanip>
#include <map>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>

namespace Quantum::Geometry::Detail
{
namespace
{
namespace LevelSetHermiteScalar
{
// Native rectangular 2:1 leaf ownership for the SAME bicubic cell kernel.
// Independent physical jets are copied at regular vertices. A used hanging
// vertex is a constraint, not an independent DOF: both its value trace and
// normal-derivative trace come from the coarse edge's Hermite polynomial.
// Keep exact numerators/positive denominators through these restrictions;
// independently rounded hanging doubles must never become phase authority.
// All data are owned. No uniform background grid, solid extension, field
// reconstruction, interface selection, force law or IO is performed here.
// The public opaque geometry Snapshot owns this implementation when explicitly
// constructed for native leaves; implementation types never cross a TU ABI.
class LeafSnapshot
{
    using Exact = ExactExpansion::Expansion;
    template <std::size_t N> struct Homogeneous
    {
        std::array<Exact, N> numerator{};
        Exact denominator = ExactExpansion::fromDouble(1.0);
    };
    using Jet = Homogeneous<4>; // physical phi, phi_x, phi_y, phi_xy
    struct Leaf
    {
        Vector2D lower, upper;
        std::array<int, 4> vertices{}; // SW, SE, NW, NE
        std::array<Node, 4> approximateNodes{};
        Homogeneous<16> dofs; // row-major Hermite basis order, dimensionless jets
        CellData cache;
        bool active = false;
        int sectorIndex = -1, sectorCorner = -1, sectorPhaseSign = 0;
    };
    struct Constraint
    {
        std::array<int, 2> parents;
        // Physical output jets by the two parents' physical input jets.
        std::array<std::array<long double, 8>, 4> derivative{};
    };

  public:
    static constexpr std::size_t dofCount=16;
    using Polynomial = Homogeneous<7>;
    LeafSnapshot(const QuadtreeMesh &mesh, const ActiveHermitePatch &patch, bool radialAxisSymmetry = false)
        : locator_(mesh), cells_(mesh.cellIndexNum()), jets_(mesh.vertIndexNum()),
          constraintByVertex_(mesh.vertIndexNum(), -1),
          lower_(mesh.minCoord), upper_(mesh.maxCoord), radialAxisSymmetry_(radialAxisSymmetry)
    {
        ExactExpansion::validateEnvironment();
        Detail::require(patch.activeCells.size() == cells_.size() && patch.nodes.size() == jets_.size(),
                        "native Hermite patch requires full native cell/vertex arrays");
        Detail::require(!radialAxisSymmetry || mesh.minCoord.x == 0.0,
                        "native radial symmetry requires an actual r=0 boundary");
        std::vector<std::uint8_t> used(jets_.size(), 0), state(jets_.size(), 0);
        bool anyActive = false;
        for (int i = 0; i < mesh.cellIndexNum(); ++i)
        {
            Detail::require(patch.activeCells[i] <= 1, "native Hermite activity mask must be binary");
            const auto &cell = mesh.cells(mesh.cellIndexToLabel[i]);
            auto &leaf = cells_[i];
            leaf.active = patch.activeCells[i] != 0;
            leaf.lower = mesh.verts(cell.cellToVert[0]).coord();
            leaf.upper = mesh.verts(cell.cellToVert[2]).coord();
            constexpr int corners[4] = {0, 1, 3, 2};
            for (int k = 0; k < 4; ++k)
            {
                const int vertex = mesh.verts(cell.cellToVert[corners[k]]).index();
                leaf.vertices[k] = vertex;
                if (leaf.active)
                    used[vertex] = 1;
            }
            anyActive = anyActive || leaf.active;
        }
        Detail::require(anyActive, "native Hermite patch must contain an active leaf");
        std::function<void(int)> resolve = [&](int index) {
            if (state[index] == 2)
                return;
            Detail::require(state[index] == 0, "cyclic native hanging-node constraint");
            state[index] = 1;
            const auto &vertex = mesh.verts(mesh.vertIndexToLabel[index]);
            if (vertex.supCell >= 0)
            {
                const int coarse = mesh.cells(vertex.supCell).index();
                Detail::require(cells_[coarse].active,
                                "active hanging scalar support crosses an inactive coarse leaf");
                // supVert belongs to the old interpolation stencil on the
                // OPPOSITE coarse side; it is not this hanging point's edge.
                // Locate the actual shared side from retained leaf corners.
                const auto &owner = cells_[coarse];
                const auto point = vertex.coord();
                int first = -1, last = -1;
                if (point.x == owner.lower.x)
                {
                    first = 0;
                    last = 2;
                }
                else if (point.x == owner.upper.x)
                {
                    first = 1;
                    last = 3;
                }
                else if (point.y == owner.lower.y)
                {
                    first = 0;
                    last = 1;
                }
                else if (point.y == owner.upper.y)
                {
                    first = 2;
                    last = 3;
                }
                Detail::require(first >= 0, "hanging vertex is not on its retained coarse leaf boundary");
                const auto &a = mesh.verts(mesh.vertIndexToLabel[owner.vertices[first]]);
                const auto &b = mesh.verts(mesh.vertIndexToLabel[owner.vertices[last]]);
                used[a.index()] = used[b.index()] = 1;
                resolve(a.index());
                resolve(b.index());
                jets_[index] =
                    constrain(jets_[a.index()], jets_[b.index()], a.coord(), b.coord(), vertex.coord());
                // Differentiate the SAME linear coarse-edge restriction by
                // applying it to its eight unit inputs, not a second set of
                // hand-written midpoint weights. Decimal edges need not be
                // exact half charts. Only the terminal Jacobian is floating;
                // the scalar/root authority above is unchanged.
                Constraint constraint{{{a.index(), b.index()}}};
                for (int input = 0; input < 8; ++input)
                {
                    Jet first, last;
                    (input < 4 ? first : last).numerator[input % 4] = number(1.0);
                    const auto column = constrain(first, last, a.coord(), b.coord(), vertex.coord());
                    for (int output = 0; output < 4; ++output)
                        constraint.derivative[output][input] =
                            approximate(column.numerator[output], column.denominator);
                }
                constraintByVertex_[index] = static_cast<int>(constraints_.size());
                constraints_.push_back(std::move(constraint));
                ++hangingCount_;
            }
            else
            {
                const auto &dof = patch.nodes[index];
                if (radialAxisSymmetry && vertex.coord().x == 0.0)
                    Detail::require(dof.dx == 0.0 && dof.dxy == 0.0,
                                    "native radial-axis nodal dx and dxy must be exactly zero");
                jets_[index].numerator = {
                    {number(dof.value), number(dof.dx), number(dof.dy), number(dof.dxy)}};
            }
            state[index] = 2;
        };
        for (int i = 0; i < mesh.vertIndexNum(); ++i)
            if (used[i])
                resolve(i);
        usedNodes_ = std::move(used);
        const auto sectors = prepareSectors(mesh, patch);
        for (auto &cell : cells_)
        {
            if (!cell.active)
                continue;
            std::array<const Jet *, 4> corners{{&jets_[cell.vertices[0]], &jets_[cell.vertices[1]],
                                                &jets_[cell.vertices[2]], &jets_[cell.vertices[3]]}};
            if (cell.sectorIndex >= 0)
                corners[cell.sectorCorner] = &sectors[cell.sectorIndex];
            const auto common = commonJets(corners);
            cell.dofs.denominator = common.denominator;
            const auto hx = difference(cell.upper.x, cell.lower.x);
            const auto hy = difference(cell.upper.y, cell.lower.y);
            for (int y = 0; y < 2; ++y)
                for (int x = 0; x < 2; ++x)
                {
                    const auto &jet = common.numerator[2 * y + x];
                    cell.dofs.numerator[4 * y + x] = jet[0];
                    cell.dofs.numerator[4 * y + x + 2] = ExactExpansion::multiply(jet[1], hx);
                    cell.dofs.numerator[4 * (y + 2) + x] = ExactExpansion::multiply(jet[2], hy);
                    cell.dofs.numerator[4 * (y + 2) + x + 2] =
                        ExactExpansion::multiply(ExactExpansion::multiply(jet[3], hx), hy);
                }
            normalize(cell.dofs);
            for (int y = 0; y < 2; ++y)
                for (int x = 0; x < 2; ++x)
                {
                    auto value = [&](int k) {
                        return approximate(cell.dofs.numerator[k], cell.dofs.denominator);
                    };
                    cell.approximateNodes[2 * y + x] = {value(4 * y + x), value(4 * y + x + 2),
                                                        value(4 * (y + 2) + x), value(4 * (y + 2) + x + 2)};
                }
            cell.cache = CellView::makeCell(nodePointers(cell));
            if (cell.sectorIndex >= 0)
            {
                constexpr unsigned masks[4] = {1u, 2u, 8u, 4u};
                Detail::require(CellView::expandedCellSignExceptVertices(
                                    cell.dofs.numerator, masks[cell.sectorCorner]) == cell.sectorPhaseSign,
                                "native corner sector lacks its exact single-phase certificate");
            }
        }
    }

    std::size_t cellCount() const
    {
        return cells_.size();
    }
    int hangingNodeCount() const
    {
        return hangingCount_;
    }
    std::size_t vertexCount() const { return jets_.size(); }
    Vector2D minCoord() const { return lower_; }
    Vector2D maxCoord() const { return upper_; }
    bool radialAxisSymmetry() const { return radialAxisSymmetry_; }
    std::array<Vector2D, 2> leafBounds(int index) const
    {
        cellActive(index); // Validate the ID; inactive bounds are not a scalar query.
        return {{cells_[index].lower, cells_[index].upper}};
    }
    std::vector<int> leafCandidates(Vector2D lower, Vector2D upper) const
    { return locator_.intersectingBox(lower, upper); }
    bool hasCornerSectors() const
    {
        return sectorCount_ != 0;
    }
    int cellSignExceptVertices(int index, unsigned allowedCornerMask) const
    {
        return CellView::expandedCellSignExceptVertices(activeCell(index).dofs.numerator, allowedCornerMask);
    }
    bool closedCellOutsideBand(int index, double halfWidth) const
    {
        Detail::require(std::isfinite(halfWidth) && halfWidth >= 0.0,
                        "native closed-cell band width must be finite and nonnegative");
        // Subtract a CONSTANT from the same homogeneous Hermite polynomial.
        // Only the four value DOFs change; all derivative DOFs and the positive
        // denominator are retained. The existing exact Bernstein sign test
        // proves the whole closed leaf, not a sample or a rounded min/max.
        const auto &dofs = activeCell(index).dofs;
        const auto shift = ExactExpansion::multiply(ExactExpansion::fromDouble(halfWidth), dofs.denominator);
        auto below = dofs.numerator, above = dofs.numerator;
        for (int y = 0; y < 2; ++y)
            for (int x = 0; x < 2; ++x)
            {
                above[4 * y + x] = ExactExpansion::subtract(above[4 * y + x], shift);
                below[4 * y + x] = ExactExpansion::add(below[4 * y + x], shift);
            }
        return CellView::expandedCellSignExceptVertices(above, 0u) == 1 ||
               CellView::expandedCellSignExceptVertices(below, 0u) == -1;
    }
    bool closedCellInsideBand(int index, double halfWidth) const
    {
        Detail::require(std::isfinite(halfWidth) && halfWidth >= 0.0,
                        "native closed-cell band width must be finite and nonnegative");
        const auto &dofs = activeCell(index).dofs;
        const auto shift = ExactExpansion::multiply(ExactExpansion::fromDouble(halfWidth), dofs.denominator);
        auto below = dofs.numerator, above = dofs.numerator;
        for (int y = 0; y < 2; ++y)
            for (int x = 0; x < 2; ++x)
            {
                above[4 * y + x] = ExactExpansion::subtract(above[4 * y + x], shift);
                below[4 * y + x] = ExactExpansion::add(below[4 * y + x], shift);
            }
        // Inclusive convex-hull bound: zero controls/edges/constants at the
        // layer limit are admitted. A mixed sign is inconclusive, never pass.
        const auto high = CellView::expandedBernsteinSigns(above),
                   low = CellView::expandedBernsteinSigns(below);
        return std::all_of(high.begin(), high.end(), [](int sign) { return sign <= 0; }) &&
               std::all_of(low.begin(), low.end(), [](int sign) { return sign >= 0; });
    }
    std::array<int, 4> cellVertices(int index) const
    {
        const auto& v=activeCell(index).vertices;
        return {{v[0],v[1],v[3],v[2]}}; // SW SE NE NW, native vertex indices.
    }
    std::array<int,2> zeroCornerGradientSigns(int index,int corner) const
    {
        Detail::require(corner>=0 && corner<4,"invalid native corner");
        constexpr int x[4]={0,1,1,0},y[4]={0,0,1,1};
        const auto& d=activeCell(index).dofs.numerator;
        Detail::require(ExactExpansion::sign(d[4*y[corner]+x[corner]])==0,
            "native corner-gradient query requires an exact zero corner");
        return {{ExactExpansion::sign(d[4*y[corner]+x[corner]+2]),
            ExactExpansion::sign(d[4*(y[corner]+2)+x[corner]])}};
    }
    bool resolvedCartesianPartial(int index, int axis = -1) const
    {
        const auto& cell=activeCell(index);
        const auto b=floatingBernstein(cell);
        long double scale=0;
        for(auto v:b) scale=std::max(scale,std::abs(v));
        floatingVariationScale(cell,b);
        return LevelSetHermiteScalar::resolvedCartesianPartial(b,
            512.0L*std::numeric_limits<double>::epsilon()*scale, axis);
    }
    long double derivativeResolutionScale(int index) const
    {
        const auto& cell=activeCell(index);
        return floatingVariationScale(cell,floatingBernstein(cell));
    }
    bool resolvedZeroSetPartial(int index,int axis=-1,
        Vector2D lower={0,0},Vector2D upper={1,1}) const
    {
        ExactExpansion::validateEnvironment();
        Detail::require(axis>=-1 && axis<2,"invalid zero-set partial axis");
        const auto partials=QuarticCell::zeroSetCartesianPartials(elevatedControls(index),lower,upper);
        return axis<0?partials[0] || partials[1]:partials[axis];
    }
    std::vector<QuarticCell::RegularChartBox> regularChartPartition(int index) const
    {
        ExactExpansion::validateEnvironment();
        return QuarticCell::regularChartPartition(elevatedControls(index));
    }
    std::array<int,2> localGradientSigns(int index,Vector2D point) const
    {
        const RationalSegmentCoordinate x{number(point.x),number(point.x),number(1)},
            y{number(point.y),number(point.y),number(1)};
        return {{ExactExpansion::sign(restrictLocal(index,x,y,1,0).numerator[0]),
            ExactExpansion::sign(restrictLocal(index,x,y,0,1).numerator[0])}};
    }
    std::array<Exact,25> elevatedControls(int index) const
    {
        const auto cubic=CellView::bernsteinControlsCleared(activeCell(index).dofs.numerator);
        // Exact degree elevation of the SAME cubic, with the positive common
        // denominator 4*4 left cleared. No Q4 samples, refit, or scalar-family
        // conversion. This reuses the Q4 outward zero-set partial certificate.
        constexpr int weights[5][4]={{4,0,0,0},{1,3,0,0},{0,2,2,0},{0,0,3,1},{0,0,0,4}};
        std::array<Exact,25> elevated;
        for (int y=0;y<5;++y) for (int x=0;x<5;++x)
            for (int j=0;j<4;++j) for (int i=0;i<4;++i) if (weights[y][j] && weights[x][i])
                elevated[5*y+x]=ExactExpansion::add(elevated[5*y+x],ExactExpansion::multiply(cubic[4*j+i],
                    ExactExpansion::fromDouble(weights[y][j]*weights[x][i])));
        return elevated;
    }
    bool cellActive(int index) const
    {
        Detail::require(index >= 0 && static_cast<std::size_t>(index) < cells_.size(),
                        "native Hermite cell index is out of range");
        return cells_[index].active;
    }
    Vector2D lower(int index) const
    {
        return activeCell(index).lower;
    }
    Vector2D upper(int index) const
    {
        return activeCell(index).upper;
    }
    // The floating approximation uses the original normalized CellData kernel.
    // Its independent rounding is not used to decide exact roots or signs.
    Sample evaluateCell(int index, Vector2D point) const
    {
        return evaluateCellChart(index,chartAtPoint(index,point));
    }
    CellChartPoint chartAtPoint(int index,Vector2D point) const
    {
        ExactExpansion::validateEnvironment();
        const auto &cell = activeCell(index);
        contains(cell, point);
        auto local = [](double point, double low, double high) {
            const bool upper = point - low > high - point;
            const double offset = static_cast<double>(
                upper ? (static_cast<long double>(high) - point) / (static_cast<long double>(high) - low)
                      : (static_cast<long double>(point) - low) / (static_cast<long double>(high) - low));
            return CellChartCoordinate{offset, upper};
        };
        return {local(point.x,cell.lower.x,cell.upper.x),local(point.y,cell.lower.y,cell.upper.y)};
    }
    Sample evaluateCellChart(int index,const CellChartPoint& point) const
    {
        validateCellChart(point);
        return view(activeCell(index)).evaluateCellBasis(Detail::Basis(point.x),Detail::Basis(point.y));
    }
    ThirdDerivatives thirdDerivativesCellChart(int index,const CellChartPoint& point) const
    {
        validateCellChart(point);
        return view(activeCell(index)).thirdDerivativesCellBasis(Detail::Basis(point.x),Detail::Basis(point.y));
    }
    Vector2D chartDisplayPoint(int index,const CellChartPoint& point) const
    {
        validateCellChart(point);
        const auto& cell=activeCell(index);
        const auto physical=[](CellChartCoordinate coordinate,double low,double high) {
            auto offset=ExactExpansion::multiply(ExactExpansion::subtract(ExactExpansion::fromDouble(high),
                ExactExpansion::fromDouble(low)),ExactExpansion::fromDouble(coordinate.offset));
            return ExactExpansion::estimate(coordinate.fromUpper
                ? ExactExpansion::subtract(ExactExpansion::fromDouble(high),offset)
                : ExactExpansion::add(ExactExpansion::fromDouble(low),offset));
        };
        return {physical(point.x,cell.lower.x,cell.upper.x),physical(point.y,cell.lower.y,cell.upper.y)};
    }
    Sample evaluate(Vector2D point) const
    { return evaluateCell(selectedCell(point), point); }
    int selectedCell(Vector2D point) const
    {
        const auto incident = locator_.incident(point);
        int chosen = -1;
        for (int k = 0; k < incident.count; ++k)
            if (cellActive(incident.indices[k]) && (chosen < 0 || incident.indices[k] < chosen))
                chosen = incident.indices[k];
        Detail::require(chosen >= 0, "native Hermite query is outside the active leaf domain");
        return chosen;
    }
    double valueInCell(int index, Vector2D point) const
    {
        const auto value = restrictPhysical(index, point, point);
        return approximate(value.numerator[0], value.denominator);
    }
    // P(t) / denominator is the requested physical derivative on the exact
    // represented physical segment. P has the same roots/signs when order=0.
    // Retain the denominator for cross-cell identity and physical magnitudes.
    Polynomial restrictPhysical(int index, Vector2D begin, Vector2D end, unsigned orderX = 0,
                                unsigned orderY = 0) const
    {
        const auto &cell = activeCell(index);
        contains(cell, begin);
        contains(cell, end);
        return restrictLocal(index, chart(begin.x, end.x, cell.lower.x, cell.upper.x),
                             chart(begin.y, end.y, cell.lower.y, cell.upper.y), orderX, orderY);
    }

    // Retain rational/endpoint-relative local charts; no rounded global refit.
    Polynomial restrictLocal(int index, const RationalSegmentCoordinate &x,
                             const RationalSegmentCoordinate &y,
                             unsigned orderX = 0, unsigned orderY = 0) const
    {
        const auto &cell = activeCell(index);
        for (const auto *coordinate : {&x, &y})
        {
            Detail::require(ExactExpansion::sign(coordinate->denominator) > 0,
                            "native restriction requires a positive chart denominator");
            for (const auto *numerator : {&coordinate->numeratorBegin, &coordinate->numeratorEnd})
                Detail::require(ExactExpansion::sign(*numerator) >= 0 &&
                                    ExactExpansion::sign(ExactExpansion::subtract(coordinate->denominator,
                                                                                *numerator)) >= 0,
                                "native restriction chart leaves its closed cell");
        }
        const auto bx = CellView::restrictedBasisCleared(x, orderX);
        const auto by = CellView::restrictedBasisCleared(y, orderY);
        Polynomial result;
        result.denominator = cell.dofs.denominator;
        for (int k = 0; k < 3; ++k)
            result.denominator = ExactExpansion::multiply(
                result.denominator, ExactExpansion::multiply(x.denominator, y.denominator));
        for (unsigned k = 0; k < orderX; ++k)
            result.denominator =
                ExactExpansion::multiply(result.denominator, difference(cell.upper.x, cell.lower.x));
        for (unsigned k = 0; k < orderY; ++k)
            result.denominator =
                ExactExpansion::multiply(result.denominator, difference(cell.upper.y, cell.lower.y));
        for (int row = 0; row < 4; ++row)
            for (int column = 0; column < 4; ++column)
            {
                const auto &dof = cell.dofs.numerator[4 * row + column];
                if (ExactExpansion::sign(dof) == 0)
                    continue;
                for (int v = 0; v < 4; ++v)
                    for (int u = 0; u < 4; ++u)
                    {
                        if (ExactExpansion::sign(bx[column][u]) == 0 || ExactExpansion::sign(by[row][v]) == 0)
                            continue;
                        const auto term = ExactExpansion::multiply(
                            dof, ExactExpansion::multiply(bx[column][u], by[row][v]));
                        result.numerator[u + v] = ExactExpansion::add(result.numerator[u + v], term);
                    }
            }
        normalize(result);
        return result;
    }
    std::array<double, 4> nodeJet(int index) const
    {
        Detail::require(index >= 0 && static_cast<std::size_t>(index) < jets_.size() && usedNodes_[index],
                        "native Hermite nodal query requires a used vertex");
        const auto &jet = jets_[index];
        return {
            {approximate(jet.numerator[0], jet.denominator), approximate(jet.numerator[1], jet.denominator),
             approximate(jet.numerator[2], jet.denominator), approximate(jet.numerator[3], jet.denominator)}};
    }
    // Terminal floating quadrature only. The root authority remains the
    // homogeneous cell polynomial, not these rounded normalized coefficients.
    std::array<double, 4> mappedCellNodeDofs(int index, int cornerX, int cornerY) const
    {
        Detail::require((cornerX == 0 || cornerX == 1) && (cornerY == 0 || cornerY == 1),
                        "native mapped nodal query requires a cell corner");
        const auto &node = activeCell(index).approximateNodes[2 * cornerY + cornerX];
        return {{node.value, node.x, node.y, node.xy}};
    }

    // Apply C^T to formal physical cell-corner sensitivities, SW/SE/NE/NW.
    // This is the differential of the retained scalar, not a value-only
    // interpolation: all four hanging components depend on coarse values
    // AND derivatives. Recursive support is resolved to regular vertices.
    std::vector<ActiveHermiteNodalSensitivity>
    pullbackCellSensitivity(int index, const std::array<long double, 16> &coefficients,
                            CornerSectorVariation sectors = CornerSectorVariation::RequireAbsent) const
    {
        ExactExpansion::validateEnvironment();
        const auto &cell = activeCell(index);
        Detail::require(sectors == CornerSectorVariation::RequireAbsent ||
                            sectors == CornerSectorVariation::HoldFixed,
                        "invalid native corner-sector variation policy");
        Detail::require(cell.sectorIndex < 0 || sectors == CornerSectorVariation::HoldFixed,
                        "native sector sensitivity requires an explicit sector-DOF variation model");
        for (auto value : coefficients)
            Detail::require(std::isfinite(value), "nonfinite native cell sensitivity");
        std::map<int, std::array<long double, 4>> regular;
        const auto add = [](long double &target, long double coefficient, long double weight) {
            const long double term = coefficient * weight;
            Detail::require(std::isfinite(term) && (coefficient == 0 || weight == 0 || term != 0),
                            "native sensitivity chain-rule overflow or underflow");
            target += term;
            Detail::require(std::isfinite(target), "nonfinite native sensitivity accumulation");
        };
        const auto accumulate = [&](const auto &self, int vertex,
                                    const std::array<long double, 4> &local) -> void {
            const int id = constraintByVertex_[vertex];
            if (id < 0)
            {
                auto &result = regular[vertex];
                for (int component = 0; component < 4; ++component)
                    add(result[component], local[component], 1.0L);
                return;
            }
            const auto &constraint = constraints_[id];
            for (int parent = 0; parent < 2; ++parent)
            {
                std::array<long double, 4> next{};
                for (int input = 0; input < 4; ++input)
                    for (int output = 0; output < 4; ++output)
                        add(next[input], local[output], constraint.derivative[output][4 * parent + input]);
                self(self, constraint.parents[parent], next);
            }
        };
        constexpr int corners[4] = {0, 1, 3, 2};
        for (int corner = 0; corner < 4; ++corner)
        {
            // Only this corner was replaced by the independent fixed jet.
            // The other three corners retain their full hanging-node chain rule.
            if (cell.sectorIndex >= 0 && corners[corner] == cell.sectorCorner)
                continue;
            std::array<long double, 4> local{};
            std::copy_n(coefficients.begin() + 4 * corner, 4, local.begin());
            accumulate(accumulate, cell.vertices[corners[corner]], local);
        }
        std::vector<ActiveHermiteNodalSensitivity> result;
        for (const auto &entry : regular)
            if (std::any_of(entry.second.begin(), entry.second.end(), [](long double value) { return value != 0; }))
                result.push_back({entry.first, entry.second});
        return result;
    }

  private:
    std::array<long double,16> floatingBernstein(const Leaf& cell) const
    {
        return sampledBernstein(cell.upper-cell.lower,[&](int x,int y) {
            return view(cell).evaluateCellBasis(Detail::Basis(static_cast<double>(x)),
                Detail::Basis(static_cast<double>(y)));
        });
    }
    static long double floatingVariationScale(const Leaf& cell,const std::array<long double,16>& b)
    {
        const auto bounds=std::minmax_element(b.begin(),b.end());
        const long double result=(*bounds.second-*bounds.first)/
            std::min(cell.upper.x-cell.lower.x,cell.upper.y-cell.lower.y);
        Detail::require(std::isfinite(result),"unrepresentable local derivative-resolution scale");
        return result;
    }
    CartesianLeafLocator locator_;
    std::vector<Leaf> cells_;
    std::vector<Jet> jets_;
    std::vector<int> constraintByVertex_;
    std::vector<Constraint> constraints_;
    std::vector<std::uint8_t> usedNodes_;
    int hangingCount_ = 0;
    int sectorCount_ = 0;
    Vector2D lower_, upper_;
    bool radialAxisSymmetry_ = false;
    std::vector<Jet> prepareSectors(const QuadtreeMesh &mesh, const ActiveHermitePatch &patch)
    {
        std::vector<Jet> sectors;
        std::vector<std::uint8_t> declared(jets_.size(), 0);
        for (const auto &declaration : patch.cornerSectors)
        {
            const int cellIndex = declaration.cellIndex, vertexIndex = declaration.vertexIndex;
            Detail::require(cellActive(cellIndex) && vertexIndex >= 0 &&
                                static_cast<std::size_t>(vertexIndex) < jets_.size(),
                            "native corner sector requires active cell and valid native vertex");
            auto &cell = cells_[cellIndex];
            const auto corner = std::find(cell.vertices.begin(), cell.vertices.end(), vertexIndex);
            Detail::require(corner != cell.vertices.end() && !declared[vertexIndex] && cell.sectorIndex < 0,
                            "native corner sector is not a unique incident vertex/cell");
            const auto shared = nativeCornerSharedComponents(mesh, locator_, cellIndex, vertexIndex,
                [&](int index) { return cells_[index].active; });
            const auto &dof = declaration.dof;
            Jet replacement;
            replacement.numerator = {{number(dof.value), number(dof.dx), number(dof.dy), number(dof.dxy)}};
            Detail::require(dof.value == 0.0 && ExactExpansion::sign(jets_[vertexIndex].numerator[0]) == 0,
                            "native corner sector must preserve the exact base zero");
            Detail::require(declaration.phaseSign == -1 || declaration.phaseSign == 1,
                            "native corner sector phase must be minus or plus one");
            for (int tangent = 1; tangent <= 2; ++tangent)
                if (shared[tangent])
                {
                    const auto &base = jets_[vertexIndex];
                    Detail::require(
                        ExactExpansion::sign(ExactExpansion::subtract(
                            ExactExpansion::multiply(replacement.numerator[tangent], base.denominator),
                            base.numerator[tangent])) == 0,
                        "native corner sector violates the exact active-edge C0 trace");
                }
            cell.sectorIndex = static_cast<int>(sectors.size());
            cell.sectorCorner = static_cast<int>(corner - cell.vertices.begin());
            cell.sectorPhaseSign = declaration.phaseSign;
            declared[vertexIndex] = 1;
            sectors.push_back(std::move(replacement));
        }
        sectorCount_ = static_cast<int>(sectors.size());
        return sectors;
    }
    static Exact number(double value)
    {
        return ExactExpansion::fromDouble(value);
    }
    static Exact difference(double a, double b)
    {
        return ExactExpansion::subtract(number(a), number(b));
    }
    static double approximate(const Exact &numerator, const Exact &denominator)
    {
        const long double value = static_cast<long double>(ExactExpansion::estimate(numerator)) /
                                  ExactExpansion::estimate(denominator);
        const double result = Detail::checked(value, "unrepresentable native Hermite ratio");
        Detail::require(ExactExpansion::sign(numerator) == 0 || result != 0.0,
                        "native Hermite ratio conversion underflow");
        return result;
    }
    template <std::size_t N> static void normalize(Homogeneous<N> &value)
    {
        Detail::require(ExactExpansion::sign(value.denominator) > 0,
                        "nonpositive native Hermite denominator");
        const int exponent = -std::ilogb(value.denominator.components().back());
        value.denominator = ExactExpansion::scalePow2(value.denominator, exponent);
        for (auto &entry : value.numerator)
            entry = ExactExpansion::scalePow2(entry, exponent);
    }
    template <std::size_t N> struct CommonJets
    {
        std::array<std::array<Exact, 4>, N> numerator;
        Exact denominator;
    };
    template <std::size_t N> static CommonJets<N> commonJets(const std::array<const Jet *, N> &jets)
    {
        std::vector<Exact> denominators;
        std::array<std::size_t, N> position{};
        for (std::size_t i = 0; i < N; ++i)
        {
            auto found =
                std::find_if(denominators.begin(), denominators.end(), [&](const Exact &denominator) {
                    return ExactExpansion::sign(
                               ExactExpansion::subtract(denominator, jets[i]->denominator)) == 0;
                });
            position[i] = static_cast<std::size_t>(found - denominators.begin());
            if (found == denominators.end())
                denominators.push_back(jets[i]->denominator);
        }
        CommonJets<N> result;
        result.denominator = number(1.0);
        for (const auto &denominator : denominators)
            result.denominator = ExactExpansion::multiply(result.denominator, denominator);
        for (std::size_t i = 0; i < N; ++i)
        {
            auto factor = number(1.0);
            for (std::size_t k = 0; k < denominators.size(); ++k)
                if (k != position[i])
                    factor = ExactExpansion::multiply(factor, denominators[k]);
            for (int k = 0; k < 4; ++k)
                result.numerator[i][k] = ExactExpansion::multiply(jets[i]->numerator[k], factor);
        }
        return result;
    }
    static RationalSegmentCoordinate chart(double begin, double end, double low, double high)
    {
        RationalSegmentCoordinate result{difference(begin, low), difference(end, low), difference(high, low)};
        const int exponent = -std::ilogb(result.denominator.components().back());
        result.numeratorBegin = ExactExpansion::scalePow2(result.numeratorBegin, exponent);
        result.numeratorEnd = ExactExpansion::scalePow2(result.numeratorEnd, exponent);
        result.denominator = ExactExpansion::scalePow2(result.denominator, exponent);
        return result;
    }
    static Jet constrain(const Jet &first, const Jet &last, Vector2D a, Vector2D b, Vector2D point)
    {
        const bool horizontal = a.y == b.y;
        if (!(horizontal ? (a.y == point.y && a.x < point.x && point.x < b.x)
                         : (a.x == b.x && a.x == point.x && a.y < point.y && point.y < b.y)))
        {
            std::ostringstream message;
            message << std::setprecision(17) << "[LevelSetLeafHermiteScalar] hanging node (" << point.x << ','
                    << point.y << ") is not inside its coarse edge (" << a.x << ',' << a.y << ") -> (" << b.x
                    << ',' << b.y << ')';
            throw std::runtime_error(message.str());
        }
        const double low = horizontal ? a.x : a.y, high = horizontal ? b.x : b.y;
        const auto coordinate =
            chart(horizontal ? point.x : point.y, horizontal ? point.x : point.y, low, high);
        const auto valueBasis = CellView::restrictedBasisCleared(coordinate);
        const auto slopeBasis = CellView::restrictedBasisCleared(coordinate, 1);
        const auto width = difference(high, low);
        const auto common = commonJets(std::array<const Jet *, 2>{{&first, &last}});
        Jet result;
        result.denominator = ExactExpansion::multiply(common.denominator, width);
        for (int k = 0; k < 3; ++k)
            result.denominator = ExactExpansion::multiply(result.denominator, coordinate.denominator);
        const int tangent = horizontal ? 1 : 2, normal = horizontal ? 2 : 1;
        for (int component = 0; component < 2; ++component)
        {
            const int valueIndex = component ? normal : 0, derivativeIndex = component ? 3 : tangent;
            const std::array<Exact, 4> data{
                {common.numerator[0][valueIndex], common.numerator[1][valueIndex],
                 ExactExpansion::multiply(width, common.numerator[0][derivativeIndex]),
                 ExactExpansion::multiply(width, common.numerator[1][derivativeIndex])}};
            Exact value, slope;
            for (int k = 0; k < 4; ++k)
            {
                value = ExactExpansion::add(value, ExactExpansion::multiply(data[k], valueBasis[k][0]));
                slope = ExactExpansion::add(slope, ExactExpansion::multiply(data[k], slopeBasis[k][0]));
            }
            result.numerator[valueIndex] = ExactExpansion::multiply(value, width);
            result.numerator[derivativeIndex] = std::move(slope);
        }
        normalize(result);
        return result;
    }
    const Leaf &activeCell(int index) const
    {
        Detail::require(cellActive(index), "native Hermite cell is inactive");
        return cells_[index];
    }
    static void contains(const Leaf &cell, Vector2D point)
    {
        Detail::require(std::isfinite(point.x) && std::isfinite(point.y) && point.x >= cell.lower.x &&
                            point.x <= cell.upper.x && point.y >= cell.lower.y && point.y <= cell.upper.y,
                        "native Hermite query must be in its actual closed leaf");
    }
    static std::array<const Node *, 4> nodePointers(const Leaf &cell)
    {
        return {{&cell.approximateNodes[0], &cell.approximateNodes[1], &cell.approximateNodes[2],
                 &cell.approximateNodes[3]}};
    }
    static CellView view(const Leaf &cell)
    {
        return CellView(nodePointers(cell), cell.cache, cell.upper.x - cell.lower.x,
                        cell.upper.y - cell.lower.y);
    }
};
} // namespace LevelSetHermiteScalar
} // namespace
} // namespace Quantum::Geometry::Detail
