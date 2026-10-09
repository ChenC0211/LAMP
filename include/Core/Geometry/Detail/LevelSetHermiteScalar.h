#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

#include "Core/Field/VertField.h"
#include "Core/Geometry/ActiveHermitePatch.h"
#include "Core/Geometry/Detail/LevelSetHermiteCell.h"
#include "Core/Geometry/LevelSetCurvature.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/UniformLeafGrid.h"

// Implementation-only engine. Shared nodal phi and Q4 nodal derivatives define a
// tensor-product bicubic Hermite scalar. Adjacent cells use the SAME retained
// double nodal DOFs: phi, hx*phi_x, hy*phi_y and hx*hy*phi_xy. These DOFs,
// including their original scaling roundoff, are the topology authority;
// exact shared-DOF restrictions below are C1 in exact arithmetic, except at
// explicitly declared corner-sector seams, which are only guaranteed C0.
// Every cell-local path uses that cell's SAME retained DOFs. The older
// independently normalized cell cache remains a floating approximation and
// its evaluate/evaluateCell/segmentPolynomial arithmetic is kept unchanged.
// It is generally NOT C2. No zero finding, contour, normal, curvature,
// distance-property claim or production-method selection is provided here.
// Implementation-only: the anonymous ancestor deliberately gives this engine
// TU-local linkage. No type or inline definition here is a cross-TU public ABI.
// Keep every include outside these namespaces; public callers use the opaque
// out-of-line snapshot, while legacy tests may inspect this local engine.
namespace Quantum::Geometry::Detail
{
namespace
{
namespace LevelSetHermiteScalar
{
class Snapshot
{
public:
    using ExpandedPolynomial = std::array<ExactExpansion::Expansion, 7>;
    using Interval = ExactExpansion::Interval;
    using BernsteinControls = std::array<long double, 16>;

    // Original sampled-corner Bernstein arithmetic, exposed to the common
    // geometry engine without assuming that every scalar has cubic degree.
    BernsteinControls cellBernsteinApproximate(int i, int j) const
    {
        return sampledBernstein(cellSize(), [&](int x,int y) { return evaluateCell(i,j,x,y); });
    }

    // A coordinate on one CLOSED local cell, retained without division:
    // u(t) = ((1-t)*numeratorBegin + t*numeratorEnd) / denominator.
    // The denominator is constant and strictly positive. Physical callers
    // must supply their declared map's exact differences, not rounded local
    // quotients. This type itself does not select a physical cell embedding.
    using RationalSegmentCoordinate = LevelSetHermiteScalar::RationalSegmentCoordinate;

    // Immutable and reference-free after construction. The shared Q4
    // constructor validates the complete uniform Cartesian mesh and phi.
    // Its derivative-only API also permits constant fields and stationary
    // points; a scalar interpolation must not require a nonzero gradient.
    Snapshot(const QuadtreeMesh& mesh, const VertField<double>& phi, bool radialAxisSymmetry = false,
        CartesianGridScope scope = CartesianGridScope::RootOnly)
        : radialAxisSymmetry_(radialAxisSymmetry)
    {
        const LevelSetCurvature::Reconstruction<QuadtreeMesh> derivatives(
            mesh, phi, 4, LevelSetCurvature::Basis::TensorProduct, radialAxisSymmetry_, scope);
        const auto& grid = derivatives.grid();
        nx_ = grid.nx(); ny_ = grid.ny(); lower_ = grid.minCoord(); upper_ = grid.maxCoord();
        hx_ = grid.spacing().x; hy_ = grid.spacing().y;

        std::vector<Node> nodes(vertexCount());
        for (std::size_t position = 0; position < grid.vertices().size(); ++position)
        {
            const int label = grid.vertices()[position];
            const int k = mesh.verts(label).index();
            const auto point = mesh.verts(label).coord();
            const auto sample = derivatives.evaluateDerivatives(point);
            auto& node = nodes[position];
            node.value = phi.value[k];
            node.x = Detail::checked(static_cast<long double>(sample.dx) * hx_,
                "unrepresentable scaled nodal x derivative");
            node.y = Detail::checked(static_cast<long double>(sample.dy) * hy_,
                "unrepresentable scaled nodal y derivative");
            node.xy = Detail::checked(static_cast<long double>(sample.dxy) * hx_ * hy_,
                "unrepresentable scaled nodal mixed derivative");
            Detail::require((sample.dx == 0.0 || node.x != 0.0)
                && (sample.dy == 0.0 || node.y != 0.0)
                && (sample.dxy == 0.0 || node.xy != 0.0),
                "scaled nodal derivative underflow; unsupported numeric scale");
        }

        cells_.resize(cellCount());
        for (int j = 0; j < ny_; ++j)
            for (int i = 0; i < nx_; ++i)
            {
                const std::array<const Node*, 4> corners{{&nodes[vertexIndex(i, j)],
                    &nodes[vertexIndex(i + 1, j)], &nodes[vertexIndex(i, j + 1)],
                    &nodes[vertexIndex(i + 1, j + 1)]}};
                cells_[cellIndex(i, j)] = makeCell(corners);
            }
        // Retain the SAME once-computed nodal doubles, after constructing the
        // legacy cache with its original operation sequence. No second Q4
        // reconstruction or refitting defines a competing scalar here.
        nodes_.swap(nodes);
    }

    // Explicit shared physical nodal DOFs: no Q4 refit and no values or
    // derivatives read outside the corners of active cells. The complete
    // native mesh still obeys the selected uniform-layout contract.
    Snapshot(const QuadtreeMesh& mesh, const ::Quantum::Geometry::ActiveHermitePatch& patch,
        bool radialAxisSymmetry = false, CartesianGridScope scope = CartesianGridScope::RootOnly)
        : radialAxisSymmetry_(radialAxisSymmetry), restrictedDomain_(true)
    {
        ExactExpansion::validateEnvironment();
        const UniformLeafGrid grid(mesh, scope);
        nx_ = grid.nx(); ny_ = grid.ny(); lower_ = grid.minCoord(); upper_ = grid.maxCoord();
        hx_ = grid.spacing().x; hy_ = grid.spacing().y;
        Detail::require(!radialAxisSymmetry_ || lower_.x == 0.0,
            "radial even extension requires the actual r=0 domain boundary");
        Detail::require(patch.activeCells.size() == cellCount(),
            "active-cell mask must have the full native cell-index size");
        Detail::require(patch.nodes.size() == vertexCount(),
            "active patch nodal DOFs must have the full native vertex-index size");

        activeCells_.resize(cellCount());
        usedNodes_.assign(vertexCount(), 0);
        for (int j = 0; j < ny_; ++j)
            for (int i = 0; i < nx_; ++i)
            {
                const auto position = cellIndex(i, j);
                const int index = mesh.cells(grid.cells()[position]).index();
                const auto active = patch.activeCells[index];
                Detail::require(active <= 1, "active-cell mask must be binary");
                activeCells_[position] = active;
                if (!active) continue;
                for (int y = 0; y < 2; ++y)
                    for (int x = 0; x < 2; ++x)
                        usedNodes_[vertexIndex(i + x, j + y)] = 1;
            }
        Detail::require(std::find(activeCells_.begin(), activeCells_.end(), std::uint8_t{1}) != activeCells_.end(),
            "active Hermite patch must contain at least one active cell");

        std::vector<Node> nodes(vertexCount());
        for (std::size_t position = 0; position < grid.vertices().size(); ++position)
        {
            if (!usedNodes_[position]) continue;
            const int index = mesh.verts(grid.vertices()[position]).index();
            const auto& dof = patch.nodes[index];
            Detail::require(std::isfinite(dof.value) && std::isfinite(dof.dx)
                && std::isfinite(dof.dy) && std::isfinite(dof.dxy),
                "active-cell nodal DOFs must be finite");
            if (radialAxisSymmetry_ && position % (static_cast<std::size_t>(nx_) + 1) == 0)
                Detail::require(dof.dx == 0.0 && dof.dxy == 0.0,
                    "active radial-axis nodal dx and dxy must be exactly zero");
            auto& node = nodes[position];
            node.value = dof.value;
            node.x = Detail::checked(static_cast<long double>(dof.dx) * hx_,
                "unrepresentable scaled nodal x derivative");
            node.y = Detail::checked(static_cast<long double>(dof.dy) * hy_,
                "unrepresentable scaled nodal y derivative");
            node.xy = Detail::checked(static_cast<long double>(dof.dxy) * hx_ * hy_,
                "unrepresentable scaled nodal mixed derivative");
            Detail::require((dof.dx == 0.0 || node.x != 0.0)
                && (dof.dy == 0.0 || node.y != 0.0)
                && (dof.dxy == 0.0 || node.xy != 0.0),
                "scaled nodal derivative underflow; unsupported numeric scale");
        }

        if (!patch.cornerSectors.empty())
        {
            // Install the retained base nodes before resolving cell-local
            // overrides. All caches and exact certificates use cellNode;
            // the no-sector construction below keeps its original order.
            nodes_.swap(nodes);
            initializeCornerSectors(mesh, grid, patch);
            cells_.resize(cellCount());
            for (int j = 0; j < ny_; ++j)
                for (int i = 0; i < nx_; ++i)
                {
                    if (!activeCells_[cellIndex(i, j)]) continue;
                    const std::array<const Node*, 4> corners{{&cellNode(i, j, 0, 0),
                        &cellNode(i, j, 1, 0), &cellNode(i, j, 0, 1), &cellNode(i, j, 1, 1)}};
                    cells_[cellIndex(i, j)] = makeCell(corners);
                }
            validateCornerSectors();
            return;
        }

        cells_.resize(cellCount());
        for (int j = 0; j < ny_; ++j)
            for (int i = 0; i < nx_; ++i)
            {
                if (!activeCells_[cellIndex(i, j)]) continue;
                const std::array<const Node*, 4> corners{{&nodes[vertexIndex(i, j)],
                    &nodes[vertexIndex(i + 1, j)], &nodes[vertexIndex(i, j + 1)],
                    &nodes[vertexIndex(i + 1, j + 1)]}};
                cells_[cellIndex(i, j)] = makeCell(corners);
            }
        nodes_.swap(nodes);
    }

    int nx() const { return nx_; }
    int ny() const { return ny_; }
    bool radialAxisSymmetry() const { return radialAxisSymmetry_; }
    // True for the explicit patch constructor, including an all-active mask.
    bool restrictedDomain() const { return restrictedDomain_; }
    bool cellActive(int i, int j) const
    {
        Detail::require(i >= 0 && i < nx_ && j >= 0 && j < ny_, "cell index is outside the snapshot");
        return !restrictedDomain_ || activeCells_[cellIndex(i, j)] != 0;
    }
    Vector2D minCoord() const { return lower_; }
    Vector2D maxCoord() const { return upper_; }
    Vector2D cellSize() const { return Vector2D(hx_, hy_); }
    std::size_t vertexCount() const
    { return (static_cast<std::size_t>(nx_) + 1) * (static_cast<std::size_t>(ny_) + 1); }
    std::size_t cellCount() const
    { return static_cast<std::size_t>(nx_) * static_cast<std::size_t>(ny_); }

    bool hasCornerSectors() const { return !cornerSectors_.empty(); }
    bool cellHasCornerSector(int i, int j) const
    {
        Detail::require(i >= 0 && i < nx_ && j >= 0 && j < ny_, "cell index is outside the snapshot");
        requireActiveCell(i, j);
        return !sectorByCell_.empty() && sectorByCell_[cellIndex(i, j)] >= 0;
    }

    // Read-only copy of the retained BASE nodal doubles, in this order:
    // {phi, hx*phi_x, hy*phi_y, hx*hy*phi_xy}. The last three are derivatives
    // already scaled to the dimensionless cell coordinates, NOT unscaled
    // physical derivatives. A declared cell-local sector does not change this
    // base/cap jet: use cellNodeDofs for a particular incident cell instead.
    // No normalized-cell-cache round trip is used.
    std::array<double, 4> nodeDofs(int i, int j) const
    {
        Detail::require(i >= 0 && i <= nx_ && j >= 0 && j <= ny_, "vertex index is outside the snapshot");
        Detail::require(!restrictedDomain_ || usedNodes_[vertexIndex(i, j)] != 0,
            "vertex is not consumed by an active cell");
        const auto& node = nodes_[vertexIndex(i, j)];
        return {{node.value, node.x, node.y, node.xy}};
    }

    // Authoritative DOFs for one active cell corner, with the same scaling as
    // nodeDofs. Both corner coordinates must be exactly integer 0 or 1.
    std::array<double, 4> cellNodeDofs(int i, int j, int cornerX, int cornerY) const
    {
        Detail::require(i >= 0 && i < nx_ && j >= 0 && j < ny_, "cell index is outside the snapshot");
        Detail::require((cornerX == 0 || cornerX == 1) && (cornerY == 0 || cornerY == 1),
            "cell corner coordinates must be zero or one");
        requireActiveCell(i, j);
        const auto& node = cellNode(i, j, cornerX, cornerY);
        return {{node.value, node.x, node.y, node.xy}};
    }

    Sample evaluate(const Vector2D& point) const
    {
        Detail::require(std::isfinite(point.x) && std::isfinite(point.y)
            && point.x >= lower_.x && point.x <= upper_.x
            && point.y >= lower_.y && point.y <= upper_.y,
            "query must lie in the finite physical domain; clipping is unsupported");
        const double x = gridCoordinate(point.x, lower_.x, upper_.x, hx_, nx_);
        const double y = gridCoordinate(point.y, lower_.y, upper_.y, hy_, ny_);
        const int i = std::min(static_cast<int>(x), nx_ - 1);
        const int j = std::min(static_cast<int>(y), ny_ - 1);
        return evaluateCell(i, j, x - i, y - j);
    }

    // localX/localY are dimensionless coordinates in the closed unit cell.
    // This entry point permits independent two-sided edge/vertex tests.
    // All returned derivatives are PHYSICAL x/y derivatives. On shared
    // boundaries second derivatives may differ between the chosen cells;
    // at declared sector seams first normal derivatives may also differ.
    Sample evaluateCell(int i, int j, double localX, double localY) const
    {
        Detail::require(i >= 0 && i < nx_ && j >= 0 && j < ny_, "cell index is outside the snapshot");
        Detail::require(std::isfinite(localX) && std::isfinite(localY)
            && localX >= 0.0 && localX <= 1.0 && localY >= 0.0 && localY <= 1.0,
            "local cell coordinates must be finite and in [0,1]");
        requireActiveCell(i, j);
        return evaluateCellBasis(i, j, Detail::Basis(localX), Detail::Basis(localY));
    }

    Sample evaluateCellChart(int i, int j, const CellChartPoint& point) const
    {
        Detail::require(i >= 0 && i < nx_ && j >= 0 && j < ny_, "cell index is outside the snapshot");
        validateCellChart(point);
        requireActiveCell(i, j);
        return evaluateCellBasis(i, j, Detail::Basis(point.x), Detail::Basis(point.y));
    }

private:
    // The legacy and endpoint-relative paths share the same accumulation and
    // normalized CellData. Only basis evaluation changes in the new API.
    Sample evaluateCellBasis(int i, int j, const Detail::Basis& bx, const Detail::Basis& by) const
    {
        return cellView(i, j).evaluateCellBasis(bx, by);
    }

public:
    // Separate opt-in path. Select the SAME cell/local coordinates as evaluate;
    // no finite differences, projection, neighboring-cell average or new Q4
    // reconstruction is used. At an exact grid boundary this is the selected
    // one-sided polynomial derivative, not a claim of global C3 regularity.
    ThirdDerivatives thirdDerivatives(const Vector2D& point) const
    {
        Detail::require(std::isfinite(point.x) && std::isfinite(point.y)
            && point.x >= lower_.x && point.x <= upper_.x
            && point.y >= lower_.y && point.y <= upper_.y,
            "query must lie in the finite physical domain; clipping is unsupported");
        const double x = gridCoordinate(point.x, lower_.x, upper_.x, hx_, nx_);
        const double y = gridCoordinate(point.y, lower_.y, upper_.y, hy_, ny_);
        const int i = std::min(static_cast<int>(x), nx_ - 1);
        const int j = std::min(static_cast<int>(y), ny_ - 1);
        return thirdDerivativesCell(i, j, x - i, y - j);
    }

    // Algebraic third differentiation of the SAME normalized CellData used by
    // evaluateCell. That cache is the existing floating approximation of the
    // authoritative shared nodal DOFs, not a second scalar fit. Constant cells
    // and stationary queries are valid here; only geometric normals need a
    // resolved nonzero gradient. Explicit cell selection enables edge tests.
    ThirdDerivatives thirdDerivativesCell(int i, int j, double localX, double localY) const
    {
        Detail::require(i >= 0 && i < nx_ && j >= 0 && j < ny_, "cell index is outside the snapshot");
        Detail::require(std::isfinite(localX) && std::isfinite(localY)
            && localX >= 0.0 && localX <= 1.0 && localY >= 0.0 && localY <= 1.0,
            "local cell coordinates must be finite and in [0,1]");
        requireActiveCell(i, j);
        return thirdDerivativesCellBasis(i, j, Detail::Basis(localX), Detail::Basis(localY));
    }

    ThirdDerivatives thirdDerivativesCellChart(int i, int j, const CellChartPoint& point) const
    {
        Detail::require(i >= 0 && i < nx_ && j >= 0 && j < ny_, "cell index is outside the snapshot");
        validateCellChart(point);
        requireActiveCell(i, j);
        return thirdDerivativesCellBasis(i, j, Detail::Basis(point.x), Detail::Basis(point.y));
    }

private:
    ThirdDerivatives thirdDerivativesCellBasis(int i, int j,
        const Detail::Basis& bx, const Detail::Basis& by) const
    {
        return cellView(i, j).thirdDerivativesCellBasis(bx, by);
    }

public:
    // Continuous D3(phi) reconstructed from retained nodal values, not the
    // third derivative of the piecewise C1 Hermite polynomial. Five-point
    // nodal differences are exact through degree four; bilinear interpolation
    // is second-order consistent for smooth resolved data. On-demand work
    // leaves legacy snapshot construction and its failure behavior unchanged.
    // Explicit patches do not declare the surrounding five-node stencils;
    // this reconstruction is unsupported there, even for an all-active mask.
    // The same-cell thirdDerivativesCell path remains available.
    ThirdDerivatives reconstructedThirdDerivativesCell(int i, int j, double localX, double localY) const
    {
        Detail::require(nx_ >= 4 && ny_ >= 4, "third derivative reconstruction requires five nodes per axis");
        Detail::require(i >= 0 && i < nx_ && j >= 0 && j < ny_, "cell index is outside the snapshot");
        Detail::require(std::isfinite(localX) && std::isfinite(localY)
            && localX >= 0.0 && localX <= 1.0 && localY >= 0.0 && localY <= 1.0,
            "local cell coordinates must be finite and in [0,1]");
        requireActiveCell(i, j);
        Detail::require(!restrictedDomain_,
            "nodal third derivative reconstruction is unsupported for an explicit active Hermite patch");
        std::array<long double, 4> sum{};
        for (int y = 0; y < 2; ++y)
            for (int x = 0; x < 2; ++x)
            {
                const long double weight = (x ? localX : 1.0L - localX)
                    * (y ? localY : 1.0L - localY);
                if (weight == 0.0L) continue;
                const auto node = reconstructedThirdDerivativesNode(i + x, j + y);
                sum[0] += weight * node.dxxx;
                sum[1] += weight * node.dxxy;
                sum[2] += weight * node.dxyy;
                sum[3] += weight * node.dyyy;
            }
        auto converted = [](long double value)
        {
            const double result = Detail::checked(value, "non-finite interpolated third derivative");
            Detail::require(value == 0.0L || result != 0.0, "interpolated third derivative conversion underflow");
            return result;
        };
        return {converted(sum[0]), converted(sum[1]), converted(sum[2]), converted(sum[3])};
    }

    // Physical scalar coefficients c[k] of phi(localBegin+t*(localEnd-localBegin)),
    // 0<=t<=1, in ascending powers through degree six. This algebraically
    // restricts the ACTUAL CellData/Hermite polynomial; it does not fit seven
    // evaluations. Coefficient conversion has its own floating-point rounding
    // and does not change evaluate/evaluateCell or promise bitwise agreement.
    std::array<double, 7> segmentPolynomial(int i, int j, const Vector2D& localBegin,
        const Vector2D& localEnd) const
    {
        return cellView(i, j).segmentPolynomial(localBegin, localEnd);
    }

    // Exact ascending coefficients of the polynomial defined by the retained
    // cell-local DOFs, restricted to the represented local endpoints. Expansion
    // arithmetic preserves ALL basis-conversion tails: coefficients are not
    // independently rounded to seven doubles. In particular a shared edge
    // or endpoint has the same value/zero in either adjacent cell. This does
    // not make physical-to-local coordinate conversion exact; callers must
    // provide a consistent local representation of their physical query.
    // Unsupported numeric range/environment is an explicit exception, never
    // a fallback to the legacy rounded restriction.
    ExpandedPolynomial segmentPolynomialExact(int i, int j, const Vector2D& localBegin,
        const Vector2D& localEnd) const
    {
        return cellView(i, j).segmentPolynomialExact(localBegin, localEnd);
    }

    // Positive-denominator-cleared restriction of the SAME retained cell
    // polynomial P: Q(t) = Dx^3 * Dy^3 * P(Nx(t)/Dx, Ny(t)/Dy).
    // Q has degree <= 6 and exactly the same phase signs and roots as P on
    // this segment. Its magnitude/derivatives are NOT physical phi/grad(phi).
    // No division, nodal refit, coordinate snapping or tolerance is used.
    // In particular positive widths formed from shared canonical bounds are
    // admissible even when the local endpoint quotient is not binary64.
    // This algebraic port does not make the old nominal-h physical embedding
    // equivalent to a shared-bound embedding; callers must keep a single map
    // through intersections, integration and geometric derivatives.
    ExpandedPolynomial segmentPolynomialCleared(int i, int j,
        const RationalSegmentCoordinate& x, const RationalSegmentCoordinate& y) const
    {
        return cellView(i, j).segmentPolynomialCleared(x, y);
    }

    ExactExpansion::Expansion evaluateExactCell(int i, int j, double localX, double localY) const
    {
        const Vector2D point(localX, localY);
        return segmentPolynomialExact(i, j, point, point)[0];
    }

    int evaluateSignCell(int i, int j, double localX, double localY) const
    { return ExactExpansion::sign(evaluateExactCell(i, j, localX, localY)); }

    // Outward convex-hull bound of ALL 16 bicubic Bernstein controls formed
    // directly from cell-local nodal DOFs, in original phi units. In particular
    // 1/3 and 1/9 are enclosed divisions, not rounded point coefficients.
    // A bound containing zero is INCONCLUSIVE and requires the exact path.
    // Infinite endpoints conservatively report an unsupported filter scale;
    // they do not change the authoritative scalar or drop any coefficient.
    // Optional CLOSED local subbox, including a degenerate line or point.
    // Tensor-product blossom restriction encloses the same controls with
    // outward interval arithmetic; it does not sample phi or isolate roots.
    // Default whole-cell queries retain their original control arithmetic.
    Interval cellBounds(int i, int j, Vector2D localLower = {0.0, 0.0},
        Vector2D localUpper = {1.0, 1.0}) const
    {
        return cellView(i, j).cellBounds(localLower, localUpper);
    }

    // Exact sufficient sign certificate permitting zeros only at declared
    // corners. Multiplying all Bernstein controls by nine removes the 1/3
    // and 1/9 divisions without rounding away an exact endpoint zero.
    int cellSignExceptVertices(int i, int j, unsigned allowedCornerMask) const
    {
        return cellView(i, j).cellSignExceptVertices(allowedCornerMask);
    }

    // Fast, strictly enclosing evaluation of the SAME cell-local DOF scalar at
    // represented local coordinates. No Expansion objects/heap allocation or
    // repeated floating-environment checks occur inside the arithmetic loop.
    Interval evaluateEnclosureCell(int i, int j, double localX, double localY) const
    {
        return cellView(i, j).evaluateEnclosureCell(localX, localY);
    }

private:
    int nx_, ny_;
    Vector2D lower_, upper_;
    double hx_, hy_;
    bool radialAxisSymmetry_ = false;
    bool restrictedDomain_ = false;
    std::vector<std::uint8_t> activeCells_, usedNodes_;
    std::vector<CellData> cells_;
    std::vector<Node> nodes_;
    struct SectorNode
    {
        int i = -1, j = -1, cornerX = -1, cornerY = -1;
        int phaseSign = 0;
        Node node;
    };
    std::vector<SectorNode> cornerSectors_;
    std::vector<int> sectorByCell_;

    CellView cellView(int i, int j) const
    {
        Detail::require(i >= 0 && i < nx_ && j >= 0 && j < ny_, "cell index is outside the snapshot");
        requireActiveCell(i, j);
        return CellView({{&cellNode(i, j, 0, 0), &cellNode(i, j, 1, 0),
                         &cellNode(i, j, 0, 1), &cellNode(i, j, 1, 1)}},
            cells_[cellIndex(i, j)], hx_, hy_);
    }

    // Arguments are already checked by each cell entry point. The empty
    // override table selects exactly the original shared nodal object.
    const Node& cellNode(int i, int j, int cornerX, int cornerY) const
    {
        if (!sectorByCell_.empty())
        {
            const int id = sectorByCell_[cellIndex(i, j)];
            if (id >= 0)
            {
                const auto& sector = cornerSectors_[id];
                if (sector.cornerX == cornerX && sector.cornerY == cornerY) return sector.node;
            }
        }
        return nodes_[vertexIndex(i + cornerX, j + cornerY)];
    }

    void initializeCornerSectors(const QuadtreeMesh& mesh, const UniformLeafGrid& grid,
        const ::Quantum::Geometry::ActiveHermitePatch& patch)
    {
        std::vector<int> cellPosition(mesh.cellIndexNum(), -1), vertexPosition(mesh.vertIndexNum(), -1);
        for (std::size_t p = 0; p < grid.cells().size(); ++p)
            cellPosition[mesh.cells(grid.cells()[p]).index()] = static_cast<int>(p);
        for (std::size_t p = 0; p < grid.vertices().size(); ++p)
            vertexPosition[mesh.verts(grid.vertices()[p]).index()] = static_cast<int>(p);
        sectorByCell_.assign(cellCount(), -1);
        std::vector<std::uint8_t> usedVertices(vertexCount(), 0);
        for (const auto& declaration : patch.cornerSectors)
        {
            Detail::require(declaration.cellIndex >= 0 && declaration.cellIndex < mesh.cellIndexNum()
                && declaration.vertexIndex >= 0 && declaration.vertexIndex < mesh.vertIndexNum(),
                "corner sector requires valid native cell and vertex indices");
            const int cp = cellPosition[declaration.cellIndex], vp = vertexPosition[declaration.vertexIndex];
            Detail::require(cp >= 0 && vp >= 0, "corner sector indices are not in the uniform leaf layout");
            const int i = cp % nx_, j = cp / nx_, vi = vp % (nx_ + 1), vj = vp / (nx_ + 1);
            Detail::require(vi > 0 && vi < nx_ && vj > 0 && vj < ny_,
                "corner sector requires an internal Cartesian boundary vertex, not a domain wall or axis");
            Detail::require(activeCells_[cp] && (vi == i || vi == i + 1) && (vj == j || vj == j + 1),
                "corner sector cell must be active and incident to its declared vertex");
            Detail::require(!usedVertices[vp] && sectorByCell_[cp] < 0,
                "duplicate corner sector vertex or cell");
            int activeCount = 0, inactiveI = -1, inactiveJ = -1;
            for (int y = vj - 1; y <= vj; ++y)
                for (int x = vi - 1; x <= vi; ++x)
                    if (activeCells_[cellIndex(x, y)]) ++activeCount;
                    else { inactiveI = x; inactiveJ = y; }
            Detail::require(activeCount == 3,
                "corner sector requires exactly three active incident cells and one inactive cell");
            Detail::require(std::abs(i - inactiveI) + std::abs(j - inactiveJ) == 1,
                "corner sector cell must share an edge with the inactive incident cell");
            const auto& dof = declaration.dof;
            Detail::require(std::isfinite(dof.value) && std::isfinite(dof.dx)
                && std::isfinite(dof.dy) && std::isfinite(dof.dxy), "corner sector DOFs must be finite");
            Detail::require(dof.value == 0.0 && dof.value == nodes_[vp].value,
                "corner sector must preserve the exact base zero value");
            Detail::require(declaration.phaseSign == -1 || declaration.phaseSign == 1,
                "corner sector phase sign must be minus or plus one");
            SectorNode sector;
            sector.i = i; sector.j = j; sector.cornerX = vi - i; sector.cornerY = vj - j;
            sector.phaseSign = declaration.phaseSign;
            sector.node.value = dof.value;
            sector.node.x = Detail::checked(static_cast<long double>(dof.dx) * hx_,
                "unrepresentable scaled sector x derivative");
            sector.node.y = Detail::checked(static_cast<long double>(dof.dy) * hy_,
                "unrepresentable scaled sector y derivative");
            sector.node.xy = Detail::checked(static_cast<long double>(dof.dxy) * hx_ * hy_,
                "unrepresentable scaled sector mixed derivative");
            Detail::require((dof.dx == 0.0 || sector.node.x != 0.0)
                && (dof.dy == 0.0 || sector.node.y != 0.0)
                && (dof.dxy == 0.0 || sector.node.xy != 0.0),
                "scaled sector derivative underflow; unsupported numeric scale");
            usedVertices[vp] = 1;
            sectorByCell_[cp] = static_cast<int>(cornerSectors_.size());
            cornerSectors_.push_back(sector);
        }
    }

    void validateCornerSectors() const
    {
        for (const auto& sector : cornerSectors_)
        {
            const int i = sector.i, j = sector.j;
            for (int side = 0; side < 4; ++side)
            {
                const bool vertical = side < 2;
                const int ni = i + (side == 0 ? -1 : side == 1 ? 1 : 0);
                const int nj = j + (side == 2 ? -1 : side == 3 ? 1 : 0);
                if (ni < 0 || ni >= nx_ || nj < 0 || nj >= ny_
                    || !activeCells_[cellIndex(ni, nj)]) continue;
                const int edge = side % 2;
                for (int t = 0; t < 2; ++t)
                {
                    const auto& a = cellNode(i, j, vertical ? edge : t, vertical ? t : edge);
                    const auto& b = cellNode(ni, nj, vertical ? 1 - edge : t, vertical ? t : 1 - edge);
                    // A cubic edge trace is defined by its two endpoint
                    // values and tangential derivatives. Exact equality of
                    // these SAME scaled doubles gives identical exact edge
                    // polynomials, without asserting normal/mixed continuity.
                    Detail::require(a.value == b.value && (vertical ? a.y == b.y : a.x == b.x),
                        "corner sector violates exact C0 active-edge trace continuity");
                }
            }
            const unsigned corner = sector.cornerY ? (sector.cornerX ? 2u : 3u) : (sector.cornerX ? 1u : 0u);
            Detail::require(cellSignExceptVertices(i, j, 1u << corner) == sector.phaseSign,
                "corner sector lacks its exact single-phase certificate except at the declared zero corner");
        }
    }

    // Public cell entry points retain their original range/argument checks
    // before this activity guard, and never expose inactive default storage.
    void requireActiveCell(int i, int j) const
    {
        Detail::require(!restrictedDomain_ || activeCells_[cellIndex(i, j)] != 0,
            "cell is inactive in this Hermite patch");
    }

    ThirdDerivatives reconstructedThirdDerivativesNode(int i, int j) const
    {
        // Lagrange derivative numerators at offsets 0..4 in a five-node
        // stencil. Clamping the stencil (not the query) gives reflected
        // one-sided formulas at physical boundaries, except the explicitly
        // selected radial axis where actual nodal values are extended evenly.
        static constexpr int first[5][5] = {
            {-25,48,-36,16,-3}, {-3,-10,18,-6,1}, {1,-8,0,8,-1},
            {-1,6,-18,10,3}, {3,-16,36,-48,25}};
        static constexpr int second[5][5] = {
            {35,-104,114,-56,11}, {11,-20,6,4,-1}, {-1,16,-30,16,-1},
            {-1,4,6,-20,11}, {11,-56,114,-104,35}};
        static constexpr int third[5][5] = {
            {-5,18,-24,14,-3}, {-3,10,-12,6,-1}, {-1,2,0,-2,1},
            {1,-6,12,-10,3}, {3,-14,24,-18,5}};
        const int beginX = std::clamp(i - 2, radialAxisSymmetry_ ? -2 : 0, nx_ - 4);
        const int beginY = std::clamp(j - 2, 0, ny_ - 4);
        const int offsetX = i - beginX, offsetY = j - beginY;
        const long double center = nodes_[vertexIndex(i, j)].value;
        auto difference = [&](int x, int y)
        { return static_cast<long double>(nodes_[vertexIndex(radialAxisSymmetry_ ? std::abs(x) : x, y)].value) - center; };
        long double dxxx = 0.0L, dxxy = 0.0L, dxyy = 0.0L, dyyy = 0.0L;
        for (int k = 0; k < 5; ++k)
        {
            dxxx += third[offsetX][k] * difference(beginX + k, j);
            dyyy += third[offsetY][k] * difference(i, beginY + k);
        }
        for (int y = 0; y < 5; ++y)
            for (int x = 0; x < 5; ++x)
            {
                const long double value = difference(beginX + x, beginY + y);
                dxxy += second[offsetX][x] * first[offsetY][y] * value;
                dxyy += first[offsetX][x] * second[offsetY][y] * value;
            }
        auto physical = [](long double value, double denominator, double a, double b, double c)
        {
            const long double normalized = value / denominator;
            Detail::require(std::isfinite(normalized) && (value == 0.0L || normalized != 0.0L),
                "nodal third derivative normalization overflow or underflow");
            value = normalized;
            for (double length : {a, b, c})
            {
                const long double divided = value / length;
                Detail::require(std::isfinite(divided) && (value == 0.0L || divided != 0.0L),
                    "nodal third derivative physical scaling overflow or underflow");
                value = divided;
            }
            const double result = Detail::checked(value, "non-finite nodal third derivative");
            Detail::require(value == 0.0L || result != 0.0, "nodal third derivative conversion underflow");
            return result;
        };
        if (radialAxisSymmetry_ && i == 0) { dxxx = 0.0L; dxyy = 0.0L; }
        return {physical(dxxx, 2.0, hx_, hx_, hx_), physical(dxxy, 144.0, hx_, hx_, hy_),
            physical(dxyy, 144.0, hx_, hy_, hy_), physical(dyyy, 2.0, hy_, hy_, hy_)};
    }
    std::size_t vertexIndex(int i, int j) const
    { return static_cast<std::size_t>(j) * (static_cast<std::size_t>(nx_) + 1) + i; }
    std::size_t cellIndex(int i, int j) const
    { return static_cast<std::size_t>(j) * static_cast<std::size_t>(nx_) + i; }

    static CellData makeCell(const std::array<const Node*, 4>& corners)
    { return CellView::makeCell(corners); }

    static double gridCoordinate(double point, double lower, double upper, double h, int count)
    {
        if (point == upper) return count;
        double value = Detail::checked((static_cast<long double>(point) - lower) / h,
            "unrepresentable normalized query coordinate");
        const double tolerance = 128.0 * std::numeric_limits<double>::epsilon() * count;
        Detail::require(value >= -tolerance && value <= count + tolerance,
            "physical domain and Cartesian coordinate conversion are inconsistent");
        // The physical query has ALREADY passed the strict domain check.
        // Only arithmetic roundoff of its normalized coordinate is bounded;
        // no outside physical point is brought into the domain.
        value = std::clamp(value, 0.0, static_cast<double>(count));
        return value;
    }
};
}
} // TU-local engine
} // namespace Quantum::Geometry::Detail
