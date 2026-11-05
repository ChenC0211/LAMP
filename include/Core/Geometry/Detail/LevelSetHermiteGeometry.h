#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <iomanip>
#include <map>
#include <sstream>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "Core/Geometry/Detail/LevelSetHermiteScalar.h"
#include "Core/Geometry/Detail/LevelSetHermiteDifferential.h"
#include "Core/Geometry/Detail/PolynomialSegmentRoots.h"

// Implementation-only engine. Shared nodal Hermite DOFs define the authoritative C1
// scalar. Phase/root arithmetic retains expansion tails; the legacy scaled
// cell cache supplies floating approximations to derivatives of that scalar.
// No reference shape or PL substitute is used.
// This is deliberately restricted to regular-resolved cells: a possibly
// interfacial cell must pass a floating-point sufficient sign test for one
// nonzero Cartesian partial derivative throughout its interior. Bernstein
// convex-hull bounds use a fixed arithmetic guard; this is NOT an exact-
// rational or outward-interval topology certification. Ambiguous cells are
// rejected, rather than accepting a possibly hidden internal zero loop.
// A restricted-patch cell may instead have an exact Bernstein sign
// certificate proving that its only zeros are exact nodal corners. Those
// corners still enter the complete edge-root inventory; they do not give
// the cell a nonzero closed-cell sign or declare an interface endpoint.
// The scalar is C1, NOT C2; curvature at a grid boundary uses
// the canonical (positive-coordinate, except at the domain maximum) cell.
// Restricted patches use another incident active cell only at an EXACT
// physical Cartesian edge/vertex when that canonical cell is inactive.
// projectedFrame is a MULTISTART LOCAL stationary projection, not a certified
// globally closest point, and must not be advertised as a distance oracle.
// Implementation-only: the anonymous ancestor deliberately gives this engine
// TU-local linkage. No type or inline definition here is a cross-TU public ABI.
// Keep every include outside these namespaces; public callers use the opaque
// out-of-line snapshot, while legacy tests may inspect this local engine.
namespace LAMP::Geometry::Detail
{
namespace
{
namespace LevelSetHermiteGeometry
{
struct Piece { Vector2D begin, end; bool negative = false; };
struct ParametricPiece
{
    Vector2D anchorBegin, anchorEnd;
    double tBegin = 0.0, tEnd = 1.0;
    bool negative = false;
};
struct SegmentRootIncidence
{
    int cellI = -1, cellJ = -1;
    Vector2D anchorBegin, anchorEnd, localBegin, localEnd, localPoint;
    double parameter = 0.0;
    bool pointOnly = false;
};
struct GridEdgeZeros { std::vector<Vector2D> points; bool coincident = false; std::vector<double> parameters; };
struct Frame
{
    Vector2D point, normal;
    double curvature = 0.0; // +div(n), n points from negative to positive.
    double levelSetResidual = 0.0; // |phi|/|grad phi|, physical length.
    double stationarityResidual = 0.0; // |t dot (point-query)|, physical length.
    int iterations = 0, attemptedSeeds = 0, convergedSeeds = 0;
    Vector2D cellLocalPoint;
    int cellI = -1, cellJ = -1;
    bool hasCellLocalPoint = false;
    CellChartPoint cellChartPoint;
    bool hasCellChartPoint = false;
};
struct BoundaryProjection
{
    Frame frame;
    bool endpoint = false;
    double inwardDerivative = 0.0, kktResidual = 0.0;
};
struct PatchEndpoint
{
    Vector2D point;
    int cellI = -1, cellJ = -1, inwardTangentSign = 0;
};

namespace Detail
{
inline void require(bool condition, const char* message)
{
    if (!condition)
        throw std::runtime_error(std::string("[LevelSetHermiteGeometryCandidate] ") + message);
}
inline double checked(long double value, const char* message)
{
    const double result = static_cast<double>(value);
    require(std::isfinite(value) && std::isfinite(result), message);
    return result;
}
inline Vector2D between(Vector2D a, Vector2D b, double t)
{
    if (t == 0.0) return a;
    if (t == 1.0) return b;
    return Vector2D(checked(static_cast<long double>(a.x) + t * (static_cast<long double>(b.x) - a.x),
        "unrepresentable segment x"),
        checked(static_cast<long double>(a.y) + t * (static_cast<long double>(b.y) - a.y),
        "unrepresentable segment y"));
}
}

template<class ScalarSnapshot>
class BasicSnapshot
{
public:
    using BernsteinControls = typename ScalarSnapshot::BernsteinControls;
    using ExpandedPolynomial = typename ScalarSnapshot::ExpandedPolynomial;

    explicit BasicSnapshot(ScalarSnapshot scalar)
        : scalar_(std::move(scalar)), nx_(scalar_.nx()), ny_(scalar_.ny()),
          lower_(scalar_.minCoord()), upper_(scalar_.maxCoord()), h_(scalar_.cellSize())
    { buildSeeds(); }

    BasicSnapshot(const QuadtreeMesh& mesh, const VertField<double>& phi, bool radialAxisSymmetry = false,
        CartesianGridScope scope = CartesianGridScope::RootOnly)
        : scalar_(mesh, phi, radialAxisSymmetry, scope), nx_(scalar_.nx()), ny_(scalar_.ny()),
          lower_(scalar_.minCoord()), upper_(scalar_.maxCoord()), h_(scalar_.cellSize())
    { buildSeeds(); }

    BasicSnapshot(const QuadtreeMesh& mesh, const LAMP::Geometry::ActiveHermitePatch& patch,
        bool radialAxisSymmetry = false, CartesianGridScope scope = CartesianGridScope::RootOnly)
        : scalar_(mesh, patch, radialAxisSymmetry, scope), nx_(scalar_.nx()), ny_(scalar_.ny()),
          lower_(scalar_.minCoord()), upper_(scalar_.maxCoord()), h_(scalar_.cellSize())
    { buildSeeds(); }

    int nx() const { return nx_; }
    int ny() const { return ny_; }
    bool radialAxisSymmetry() const { return scalar_.radialAxisSymmetry(); }
    bool restrictedDomain() const { return scalar_.restrictedDomain(); }
    bool hasCornerSectors() const { return scalar_.hasCornerSectors(); }
    bool cellActive(int i, int j) const { return scalar_.cellActive(i, j); }
    Vector2D minCoord() const { return lower_; }
    Vector2D maxCoord() const { return upper_; }
    Vector2D cellSize() const { return h_; }
    std::size_t candidateCellCount() const { return candidateCells_; }
    std::size_t seedCount() const { return seeds_.size(); }
    std::size_t vertexCount() const { return scalar_.vertexCount(); }

    // Internal same-source access for explicit shared-bound segment queries.
    // These do not change the nominal-h point/jet or old partition APIs.
    int strictCellSign(int i, int j) const
    {
        ExactExpansion::validateEnvironment();
        Detail::require(cellActive(i, j), "cell sign requires an active cell");
        return cellSigns_[static_cast<std::size_t>(j) * nx_ + i];
    }
    ExpandedPolynomial clearedCellSegment(int i, int j,
        const LevelSetHermiteScalar::RationalSegmentCoordinate& x,
        const LevelSetHermiteScalar::RationalSegmentCoordinate& y) const
    { return scalar_.segmentPolynomialCleared(i, j, x, y); }
    std::array<double,4> mappedCellNodeDofs(int i, int j, int x, int y) const
    { return scalar_.cellNodeDofs(i, j, x, y); }
    auto mappedCellDofs(int i,int j) const { return scalar_.cellMappedDofs(i,j); }

    double value(Vector2D point) const { return sample(point).value; }
    int closedCellSign(Vector2D point) const
    {
        const auto g = gridPoint(point);
        const auto cell = selectCell(point, g);
        const int i = cell[0], j = cell[1];
        return cellSigns_[static_cast<std::size_t>(j) * nx_ + i];
    }
    int closedCellSignExceptVertices(Vector2D point, unsigned allowedCornerMask) const
    {
        const auto g = gridPoint(point);
        const auto cell = selectCell(point, g);
        const int i = cell[0], j = cell[1];
        return scalar_.cellSignExceptVertices(i, j, allowedCornerMask);
    }
    std::array<int, 2> zeroVertexGradientSigns(Vector2D point) const
    {
        const auto g = gridPoint(point);
        selectCell(point, g);
        Detail::require(g[0] == std::round(g[0]) && g[1] == std::round(g[1])
            && point.x == coordinate(static_cast<int>(g[0]), true)
            && point.y == coordinate(static_cast<int>(g[1]), false),
            "zero-vertex gradient signs require an exact Cartesian vertex");
        const auto dofs = scalar_.nodeDofs(static_cast<int>(g[0]), static_cast<int>(g[1]));
        Detail::require(dofs[0] == 0.0, "zero-vertex gradient signs require exactly zero authoritative nodal phi");
        return {{(dofs[1] > 0.0) - (dofs[1] < 0.0), (dofs[2] > 0.0) - (dofs[2] < 0.0)}};
    }
    std::array<int, 2> zeroVertexGradientSignsInCell(int i, int j, Vector2D point) const
    {
        gridPoint(point);
        Detail::require(cellActive(i, j), "corner signs require an active incident cell");
        const int x = point.x == coordinate(i, true) ? 0 : point.x == coordinate(i+1, true) ? 1 : -1;
        const int y = point.y == coordinate(j, false) ? 0 : point.y == coordinate(j+1, false) ? 1 : -1;
        Detail::require(x >= 0 && y >= 0, "corner signs require an exact corner of the specified cell");
        const auto dofs = scalar_.cellNodeDofs(i, j, x, y);
        Detail::require(dofs[0] == 0.0, "corner signs require exactly zero authoritative phi");
        return {{(dofs[1] > 0.0) - (dofs[1] < 0.0), (dofs[2] > 0.0) - (dofs[2] < 0.0)}};
    }
    bool closedCellOutsideBand(Vector2D point, double halfWidth) const
    {
        const auto g = gridPoint(point);
        Detail::require(std::isfinite(halfWidth) && halfWidth >= 0.0,
            "band half-width must be finite and nonnegative");
        const auto cell = selectCell(point, g);
        const int i = cell[0], j = cell[1];
        const auto bound = scalar_.cellBounds(i, j);
        return bound.lo > halfWidth || bound.hi < -halfWidth;
    }
    bool closedCellInsideBand(Vector2D point, double halfWidth) const
    {
        const auto g = gridPoint(point);
        Detail::require(std::isfinite(halfWidth) && halfWidth >= 0.0,
            "band half-width must be finite and nonnegative");
        const auto cell = selectCell(point, g);
        const auto bound = scalar_.cellBounds(cell[0], cell[1]);
        return bound.lo >= -halfWidth && bound.hi <= halfWidth;
    }
    bool negative(Vector2D point) const
    {
        const auto g = gridPoint(point);
        const auto cell = selectCell(point, g);
        const int i = cell[0], j = cell[1];
        const int cellSign = cellSigns_[static_cast<std::size_t>(j) * nx_ + i];
        if (cellSign != 0) return cellSign < 0;
        const auto bound = scalar_.evaluateEnclosureCell(i, j, g[0] - i, g[1] - j);
        if (bound.lo > 0.0) return false;
        if (bound.hi < 0.0) return true;
        return scalar_.evaluateSignCell(i, j, g[0] - i, g[1] - j) <= 0;
    }
    Frame jet(Vector2D point) const;
    Frame jetInCell(int i, int j, Vector2D point) const;
    // Opt-in local derivative of +div(n) along t=(-ny,nx), in length^-2.
    // Uses the SAME canonical cell as jet; the scalar is C1, not C2/C3, so
    // this need not match the other cell's limit on a shared grid boundary.
    // The point is evaluated as given, NOT projected or averaged across cells.
    double curvatureTangentialDerivative(Vector2D point) const;
    // Original same-cell D3, including its one-sided sector. Retains the
    // caller's local chart rather than reselecting a cell from a display point.
    double curvatureTangentialDerivativeAtCellLocal(int i, int j, Vector2D localPoint) const;
    double curvatureTangentialDerivativeAtCellChart(int i, int j, const CellChartPoint& point) const;
    double reconstructedCurvatureTangentialDerivative(Vector2D point) const;
    Frame projectedFrame(Vector2D point) const;
    Frame projectedCellLocalFrame(Vector2D point) const;
    std::vector<Frame> stationaryCellLocalFrames(Vector2D point, Frame& diagnostics) const;
    BoundaryProjection boundaryProjectedFrame(Vector2D query,
        const std::vector<Vector2D>& boundaryEndpoints) const;
    BoundaryProjection patchBoundaryProjectedFrame(Vector2D query,
        const std::vector<PatchEndpoint>& boundaryEndpoints) const;
    bool segmentInsideActiveDomain(Vector2D a, Vector2D b) const;
    int certifiedSegmentSign(Vector2D a, Vector2D b) const;
    std::vector<Piece> partitionSegment(Vector2D a, Vector2D b) const;
    std::vector<ParametricPiece> partitionSegmentParametric(Vector2D a, Vector2D b, bool certifyGridSlivers = false) const;
    std::vector<SegmentRootIncidence> intersectSegmentWithProvenance(Vector2D a, Vector2D b,
        bool certifyGridSlivers = true) const;
    GridEdgeZeros intersectGridEdge(Vector2D a, Vector2D b) const;
    std::vector<Vector2D> intersectCellSegment(int i, int j, Vector2D localA, Vector2D localB) const;
    Frame jetAtCellLocal(int i, int j, Vector2D localPoint) const;
    Frame jetAtCellChart(int i, int j, const CellChartPoint& point) const;
    std::vector<CellChartPoint> intersectCellSection(int i, int j,
        bool coordinateIsX, CellChartCoordinate fixedCoordinate) const;
    LevelSetHermiteScalar::Sample scalarJetAtCellLocal(int i, int j, Vector2D localPoint) const
    {
        return scalar_.evaluateCell(i, j, localPoint.x, localPoint.y);
    }
    LevelSetHermiteScalar::Sample scalarJetAtCellChart(int i, int j, const CellChartPoint& point) const
    {
        return scalar_.evaluateCellChart(i, j, point);
    }
    double negativeFraction(Vector2D a, Vector2D b) const;

private:
    bool segmentCoveredByCells(Vector2D a, Vector2D b, int requiredSign) const;
    std::vector<Piece> partitionSegmentImpl(Vector2D a, Vector2D b,
        std::vector<ParametricPiece>* parametric, bool certifyGridSlivers = false,
        std::vector<SegmentRootIncidence>* rootIncidences = nullptr) const;
    double evaluateCurvatureDerivative(Vector2D point, bool reconstructThird) const;
    static double curvatureDerivativeFromSample(const LevelSetHermiteScalar::Sample& sample,
        const LevelSetHermiteScalar::ThirdDerivatives& third, long double gradient);
    ScalarSnapshot scalar_;
    int nx_, ny_;
    Vector2D lower_, upper_, h_;
    std::size_t candidateCells_ = 0;
    std::vector<Vector2D> seeds_;
    std::vector<long double> derivativeResolutionScales_;
    std::vector<int> cellSigns_; // Nonzero only after a rigorous closed-cell sign bound.
    // Restricted cells certified to have zeros ONLY at these exact corners.
    // Kept separate from cellSigns_: their closed cells are not strictly signed.
    std::vector<unsigned> cornerOnlyMasks_;

    bool inRectangle(Vector2D point) const
    {
        return std::isfinite(point.x) && std::isfinite(point.y)
            && point.x >= lower_.x && point.x <= upper_.x
            && point.y >= lower_.y && point.y <= upper_.y;
    }
    bool inDomain(Vector2D point) const
    {
        if (!inRectangle(point)) return false;
        if (!restrictedDomain()) return true;
        return findCell(point, gridPoint(point))[0] >= 0;
    }
    double coordinate(int k, bool x) const
    {
        const int count = x ? nx_ : ny_;
        if (k == count) return x ? upper_.x : upper_.y;
        return (x ? lower_.x : lower_.y) + k * (x ? h_.x : h_.y);
    }
    std::array<double, 2> gridPoint(Vector2D point) const
    {
        ExactExpansion::validateEnvironment(); // Cached sign filters cannot bypass this contract.
        Detail::require(inRectangle(point), "query must lie in the finite physical domain; clipping is unsupported");
        std::array<double, 2> result;
        for (int axis = 0; axis < 2; ++axis)
        {
            const bool x = axis == 0;
            const double p = x ? point.x : point.y, upper = x ? upper_.x : upper_.y;
            const int count = x ? nx_ : ny_;
            double normalized = p == upper ? count : Detail::checked(
                (static_cast<long double>(p) - (x ? lower_.x : lower_.y)) / (x ? h_.x : h_.y),
                "unrepresentable Cartesian query coordinate");
            const double guard = 128.0 * std::numeric_limits<double>::epsilon() * count;
            Detail::require(normalized >= -guard && normalized <= count + guard,
                "physical domain and Cartesian coordinates disagree");
            normalized = std::clamp(normalized, 0.0, static_cast<double>(count));
            const int nearest = static_cast<int>(std::llround(normalized));
            // Only an EXACT physical grid coordinate is canonicalized; nearby
            // physical points are not snapped or reclassified by a tolerance.
            if (p == coordinate(nearest, x)) normalized = nearest;
            result[axis] = normalized;
        }
        return result;
    }
    bool inClosedCell(Vector2D point, int i, int j) const
    {
        return point.x >= coordinate(i, true) && point.x <= coordinate(i + 1, true)
            && point.y >= coordinate(j, false) && point.y <= coordinate(j + 1, false);
    }
    std::array<int, 2> findCell(Vector2D point, const std::array<double, 2>& g) const
    {
        const int i = std::min(static_cast<int>(g[0]), nx_ - 1);
        const int j = std::min(static_cast<int>(g[1]), ny_ - 1);
        // The unrestricted selection and arithmetic remain unchanged. In a
        // restricted patch, division roundoff cannot authorize an inactive
        // physical point merely by rounding it onto an active-side gridline.
        if (!restrictedDomain()) return {{i, j}};
        if (cellActive(i, j) && inClosedCell(point, i, j)) return {{i, j}};
        // Normalization may round an interior physical point across an
        // integer index even though it is not exactly on that gridline.
        // Recover its actual containing cell, not a snapped representative:
        // every alternative must contain the unchanged point exactly.
        const std::array<int, 3> xs{{i, i - 1, i + 1}}, ys{{j, j - 1, j + 1}};
        for (int b : ys) for (int a : xs)
            if (a >= 0 && a < nx_ && b >= 0 && b < ny_
                && cellActive(a, b) && inClosedCell(point, a, b)) return {{a, b}};
        return {{-1, -1}};
    }
    std::array<int, 2> selectCell(Vector2D point, const std::array<double, 2>& g) const
    {
        const auto cell = findCell(point, g);
        Detail::require(cell[0] >= 0, "query has no active incident cell; inactive-domain clipping is unsupported");
        return cell;
    }
    std::array<int, 2> selectEdgeCell(int x, int y, bool vertical) const
    {
        if (!restrictedDomain()) return {{std::min(x, nx_ - 1), std::min(y, ny_ - 1)}};
        // The midpoint disambiguates an edge from other cells meeting only
        // its endpoint. Its fixed coordinate is an EXACT physical gridline.
        const Vector2D a(coordinate(x, true), coordinate(y, false));
        const Vector2D b(coordinate(x + (vertical ? 0 : 1), true),
            coordinate(y + (vertical ? 1 : 0), false));
        return selectCell(Detail::between(a, b, 0.5),
            {{static_cast<double>(x) + (vertical ? 0.0 : 0.5),
              static_cast<double>(y) + (vertical ? 0.5 : 0.0)}});
    }
    LevelSetHermiteScalar::Sample sample(Vector2D point) const
    {
        const auto g = gridPoint(point);
        const auto cell = selectCell(point, g);
        const int i = cell[0], j = cell[1];
        return sampleInCell(i, j, g);
    }
    LevelSetHermiteScalar::Sample sampleInCell(int i, int j, const std::array<double, 2>& g) const
    {
        auto result = scalar_.evaluateCell(i, j, g[0] - i, g[1] - j);
        const int cellSign = cellSigns_[static_cast<std::size_t>(j) * nx_ + i];
        const int approximateSign = (result.value > 0.0) - (result.value < 0.0);
        // The old cache is an approximation to the same shared-DOF scalar,
        // not a second topology model. It can supply a finite value only when
        // its sign agrees with an authoritative interval certificate. Near
        // zero or on any disagreement, retain the full exact-expansion value.
        if (cellSign != 0 && cellSign == approximateSign) return result;
        const auto bound = scalar_.evaluateEnclosureCell(i, j, g[0] - i, g[1] - j);
        const int boundedSign = bound.lo > 0.0 ? 1 : bound.hi < 0.0 ? -1 : 0;
        if (boundedSign != 0 && boundedSign == approximateSign
            && result.value >= bound.lo && result.value <= bound.hi) return result;
        result.value = ExactExpansion::estimate(
            scalar_.evaluateExactCell(i, j, g[0] - i, g[1] - j));
        return result;
    }
    ExactExpansion::Expansion exactValue(Vector2D point) const
    {
        const auto g = gridPoint(point);
        const auto cell = selectCell(point, g);
        const int i = cell[0], j = cell[1];
        return scalar_.evaluateExactCell(i, j, g[0] - i, g[1] - j);
    }
    bool usableGradient(const LevelSetHermiteScalar::Sample& s, long double gradient, Vector2D point) const
    {
        const auto g = gridPoint(point);
        const auto cell = selectCell(point, g);
        const int i = cell[0], j = cell[1];
        return usableGradientInCell(s, gradient, i, j);
    }
    bool usableGradientInCell(const LevelSetHermiteScalar::Sample& s, long double gradient, int i, int j) const
    {
        const long double localVariationScale = derivativeResolutionScales_[static_cast<std::size_t>(j) * nx_ + i];
        return HermiteDifferential::usableGradient(s,gradient,localVariationScale,h_);
    }
    BernsteinControls bernstein(int i, int j) const;
    void buildSeeds();
    bool project(Vector2D query, Vector2D seed, Frame& result, bool retainCellLocal = false) const;
    bool stationaryFrame(Vector2D query, Frame& best, bool retainCellLocal = false,
        std::vector<Frame>* candidates = nullptr) const;
    struct Endpoint { Frame frame; Vector2D inward; };
    BoundaryProjection selectBoundaryProjection(Vector2D query, const std::vector<Endpoint>& endpoints,
        long double length, long double tolerance) const;
    static Frame frameFromSample(Vector2D point, const LevelSetHermiteScalar::Sample& s, long double gradient);
};

// Existing callers retain the cubic specialization and public ABI. The same
// engine can consume retained quartic data without copying its geometry logic.
using Snapshot = BasicSnapshot<LevelSetHermiteScalar::Snapshot>;

template<class ScalarSnapshot>
inline typename BasicSnapshot<ScalarSnapshot>::BernsteinControls BasicSnapshot<ScalarSnapshot>::bernstein(int i, int j) const
{
    return scalar_.cellBernsteinApproximate(i,j);
}

template<class ScalarSnapshot>
inline void BasicSnapshot<ScalarSnapshot>::buildSeeds()
{
    // A canonical side is used ONCE for each physical edge. Neighboring cells
    // never independently manufacture slightly separated copies of its root.
    std::set<std::pair<std::size_t, int>> edges; // horizontal=(j*(nx+1)+i,0), vertical=(...,1)
    derivativeResolutionScales_.resize(scalar_.cellCount());
    cellSigns_.resize(scalar_.cellCount(), 0);
    if (restrictedDomain()) cornerOnlyMasks_.resize(scalar_.cellCount(), 0);
    auto vertexIndex = [&](int i, int j)
    { return static_cast<std::size_t>(j) * (static_cast<std::size_t>(nx_) + 1) + i; };
    for (int j = 0; j < ny_; ++j)
        for (int i = 0; i < nx_; ++i)
        {
            if (!cellActive(i, j)) continue;
            const auto b = bernstein(i, j);
            long double scale = 0.0L;
            for (auto v : b) scale = std::max(scale, std::abs(v));
            const long double guard = 512.0L * std::numeric_limits<double>::epsilon() * scale;
            const auto bounds = std::minmax_element(b.begin(), b.end());
            const long double variationScale = (*bounds.second - *bounds.first) / std::min(h_.x, h_.y);
            Detail::require(std::isfinite(variationScale), "unrepresentable local derivative-resolution scale");
            derivativeResolutionScales_[static_cast<std::size_t>(j) * nx_ + i] = variationScale;
            const auto authoritativeBound = scalar_.cellBounds(i, j);
            const int cellSign = authoritativeBound.lo > 0.0 ? 1 : authoritativeBound.hi < 0.0 ? -1 : 0;
            cellSigns_[static_cast<std::size_t>(j) * nx_ + i] = cellSign;
            if (cellSign != 0) continue;
            bool onlyZeroCorners = false;
            if (restrictedDomain())
            {
                unsigned zeroCorners = 0;
                const int corners[4][2] = {{0,0}, {1,0}, {1,1}, {0,1}}; // SW, SE, NE, NW.
                for (unsigned k = 0; k < 4; ++k)
                    if (scalar_.nodeDofs(i + corners[k][0], j + corners[k][1])[0] == 0.0)
                        zeroCorners |= 1u << k;
                // This certificate uses the authoritative scaled DOFs and
                // exact control signs, not near-zero nodal values. Common
                // weak sign plus a strict control in the cell and on every
                // edge excludes every zero except the listed corners.
                onlyZeroCorners = zeroCorners != 0
                    && scalar_.cellSignExceptVertices(i, j, zeroCorners) != 0;
                if (onlyZeroCorners)
                    cornerOnlyMasks_[static_cast<std::size_t>(j) * nx_ + i] = zeroCorners;
            }
            bool regular = false;
            if (!onlyZeroCorners) regular=LevelSetHermiteScalar::resolvedCartesianPartial(b,guard);
            if (restrictedDomain())
                Detail::require(regular || onlyZeroCorners,
                    ("unresolved zero topology: restricted candidate cell (" + std::to_string(i) + ","
                    + std::to_string(j) + ") fails the floating-point Cartesian-partial sign test and has no exact corner-only sign certificate; hidden loops/critical zeros are unsupported").c_str());
            else Detail::require(regular,
                "unresolved zero topology: candidate cell fails the floating-point Cartesian-partial sign test; hidden loops/critical zeros are unsupported");
            // Keep cellSigns_ == 0 and collect all four edges even for a
            // corner-only cell. Panels must still prove local incidence and
            // enforce its declared-endpoint/complete-chain contract.
            ++candidateCells_;
            edges.insert({vertexIndex(i, j), 0});
            edges.insert({vertexIndex(i, j + 1), 0});
            edges.insert({vertexIndex(i, j), 1});
            edges.insert({vertexIndex(i + 1, j), 1});
        }
    std::set<std::pair<double, double>> unique;
    for (const auto& edge : edges)
    {
        const int x = static_cast<int>(edge.first % (static_cast<std::size_t>(nx_) + 1));
        const int y = static_cast<int>(edge.first / (static_cast<std::size_t>(nx_) + 1));
        const bool vertical = edge.second != 0;
        const auto cell = selectEdgeCell(x, y, vertical);
        const int i = cell[0], j = cell[1];
        const Vector2D localA(x - i, y - j);
        const Vector2D localB(x - i + (vertical ? 0 : 1), y - j + (vertical ? 1 : 0));
        const auto coefficients = scalar_.segmentPolynomialExact(i, j, localA, localB);
        auto roots = hasCornerSectors() ? PolynomialSegmentRoots::isolateEndpointRelative(coefficients)
            : PolynomialSegmentRoots::isolate(coefficients);
        if (roots.identicallyZero) roots.roots = {0.0, 0.5, 1.0};
        const Vector2D a(coordinate(x, true), coordinate(y, false));
        const Vector2D b(coordinate(x + (vertical ? 0 : 1), true), coordinate(y + (vertical ? 1 : 0), false));
        for (double root : roots.roots)
        {
            const auto point = Detail::between(a, b, root);
            if (unique.insert({point.x, point.y}).second) seeds_.push_back(point);
        }
    }
}

template<class ScalarSnapshot>
inline GridEdgeZeros BasicSnapshot<ScalarSnapshot>::intersectGridEdge(Vector2D a, Vector2D b) const
{
    const auto ga = gridPoint(a), gb = gridPoint(b);
    const auto vertex = [&](Vector2D p, const std::array<double, 2>& g)
    {
        Detail::require(g[0] == std::round(g[0]) && g[1] == std::round(g[1])
            && p.x == coordinate(static_cast<int>(g[0]), true)
            && p.y == coordinate(static_cast<int>(g[1]), false),
            "grid-edge intersections require exact Cartesian vertices");
    };
    vertex(a, ga); vertex(b, gb);
    Detail::require((ga[0] == gb[0] && std::abs(ga[1]-gb[1]) == 1.0)
        || (ga[1] == gb[1] && std::abs(ga[0]-gb[0]) == 1.0),
        "grid-edge intersections require one complete adjacent-vertex edge");
    if (std::make_pair(b.x, b.y) < std::make_pair(a.x, a.y))
    {
        auto result = intersectGridEdge(b, a);
        std::reverse(result.points.begin(), result.points.end());
        std::reverse(result.parameters.begin(), result.parameters.end());
        for (auto& parameter : result.parameters) parameter = 1.0 - parameter;
        return result;
    }
    const auto cell = selectEdgeCell(static_cast<int>(ga[0]), static_cast<int>(ga[1]), a.x == b.x);
    const int i = cell[0], j = cell[1];
    const auto polynomial = scalar_.segmentPolynomialExact(i, j,
        {ga[0]-i, ga[1]-j}, {gb[0]-i, gb[1]-j});
    namespace Exact = ExactExpansion;
    Detail::require(Exact::sign(Exact::evaluate(polynomial, 0.0)) == Exact::sign(exactValue(a))
        && Exact::sign(Exact::evaluate(polynomial, 1.0)) == Exact::sign(exactValue(b)),
        "grid-edge restriction and canonical endpoint values disagree");
    const auto isolated = hasCornerSectors() ? PolynomialSegmentRoots::isolateEndpointRelative(polynomial)
        : PolynomialSegmentRoots::isolate(polynomial);
    GridEdgeZeros result;
    result.coincident = isolated.identicallyZero;
    if (result.coincident) return result;
    const bool horizontal = a.y == b.y;
    for (double root : isolated.roots)
    {
        Detail::require(std::isfinite(root) && root >= 0.0 && root <= 1.0, "invalid grid-edge root parameter");
        const auto point = Detail::between(a, b, root);
        const double position = horizontal ? point.x : point.y;
        Detail::require((root == 0.0 || position > (horizontal ? a.x : a.y))
            && (root == 1.0 || position < (horizontal ? b.x : b.y)),
            "interior grid-edge root collapsed onto an endpoint");
        Detail::require(result.points.empty() || position > (horizontal ? result.points.back().x : result.points.back().y),
            "distinct grid-edge roots collapse to one represented point");
        result.points.push_back(point);
        result.parameters.push_back(root);
    }
    return result;
}

template<class ScalarSnapshot>
inline Frame BasicSnapshot<ScalarSnapshot>::jet(Vector2D point) const
{
    const auto s = sample(point);
    const long double gradient = std::hypot(static_cast<long double>(s.dx), static_cast<long double>(s.dy));
    Detail::require(usableGradient(s, gradient, point), "zero, unresolved, or non-finite Hermite gradient cannot define an interface jet");
    return frameFromSample(point, s, gradient);
}

template<class ScalarSnapshot>
inline std::vector<Vector2D> BasicSnapshot<ScalarSnapshot>::intersectCellSegment(int i, int j, Vector2D a, Vector2D b) const
{
    Detail::require(cellActive(i, j), "local segment requires an active cell");
    Detail::require(a.x != b.x || a.y != b.y, "local segment must have positive length");
    const auto polynomial = scalar_.segmentPolynomialExact(i, j, a, b);
    const auto roots = PolynomialSegmentRoots::isolateEndpointRelative(polynomial);
    Detail::require(!roots.identicallyZero, "coincident local zero segment is unsupported");
    std::vector<Vector2D> result;
    for (double root : roots.roots)
    {
        const auto point = Detail::between(a, b, root);
        Detail::require((root == 0.0 || point.x != a.x || point.y != a.y)
            && (root == 1.0 || point.x != b.x || point.y != b.y),
            "interior local root collapsed to an endpoint");
        Detail::require(result.empty() || point.x != result.back().x || point.y != result.back().y,
            "distinct local roots collapse to one represented chart point");
        result.push_back(point);
    }
    return result;
}

template<class ScalarSnapshot>
inline Frame BasicSnapshot<ScalarSnapshot>::jetAtCellLocal(int i, int j, Vector2D local) const
{
    Detail::require(cellActive(i, j), "local jet requires an active cell");
    const auto sample = scalar_.evaluateCell(i, j, local.x, local.y);
    const long double gradient = std::hypot(static_cast<long double>(sample.dx), static_cast<long double>(sample.dy));
    Detail::require(usableGradientInCell(sample, gradient, i, j), "unresolved cell-local interface gradient");
    const Vector2D point(local.x == 1.0 ? coordinate(i+1, true) : coordinate(i, true) + h_.x*local.x,
        local.y == 1.0 ? coordinate(j+1, false) : coordinate(j, false) + h_.y*local.y);
    auto result = frameFromSample(point, sample, gradient);
    result.cellLocalPoint = local; result.cellI = i; result.cellJ = j; result.hasCellLocalPoint = true;
    return result;
}

template<class ScalarSnapshot>
inline std::vector<CellChartPoint> BasicSnapshot<ScalarSnapshot>::intersectCellSection(int i, int j,
    bool coordinateIsX, CellChartCoordinate fixedCoordinate) const
{
    namespace Exact = ExactExpansion;
    using Coordinate = LevelSetHermiteScalar::RationalSegmentCoordinate;
    validateCellChart({fixedCoordinate, fixedCoordinate});
    Detail::require(cellActive(i, j), "cell section requires an active cell");
    const auto one = Exact::fromDouble(1.0), half = Exact::fromDouble(0.5);
    const auto offset = Exact::fromDouble(fixedCoordinate.offset);
    const auto fixedValue = fixedCoordinate.fromUpper ? Exact::subtract(one, offset) : offset;
    const Coordinate fixed{fixedValue, fixedValue, one};
    std::array<std::vector<CellChartPoint>, 2> halves;
    for (int side = 0; side < 2; ++side)
    {
        // Both restrictions start at an original cell edge and end at the
        // same EXACT midpoint. Upper-half roots retain their distance from
        // one, rather than materializing 1-root/2 in a binary64 coordinate.
        const Coordinate variable{Exact::fromDouble(side ? 1.0 : 0.0), half, one};
        const auto polynomial = scalar_.segmentPolynomialCleared(i, j,
            coordinateIsX ? fixed : variable, coordinateIsX ? variable : fixed);
        const auto isolated = PolynomialSegmentRoots::isolateEndpointRelative(polynomial);
        Detail::require(!isolated.identicallyZero, "coincident cell section is unsupported");
        Detail::require(isolated.roots.size() == isolated.enclosures.size(),
            "cell section root/enclosure count mismatch");
        for (std::size_t k = 0; k < isolated.roots.size(); ++k)
        {
            const double root = isolated.roots[k];
            const auto bracket = isolated.enclosures[k];
            Detail::require(root >= 0.0 && root <= 1.0 && bracket.lo >= 0.0 && bracket.hi <= 1.0
                && bracket.lo <= root && root <= bracket.hi,
                "invalid cell section root enclosure");
            Detail::require((root != 0.0 || Exact::sign(polynomial[0]) == 0)
                && (root != 1.0 || Exact::sign(Exact::evaluate(polynomial, 1.0)) == 0),
                "interior cell section root collapsed to a nonzero chart endpoint");
            if (root == 1.0)
            {
                // The half-chart endpoint is INTERNAL to the original
                // section. Splitting must not turn an unsupported internal
                // multiple root into two otherwise admissible endpoint roots.
                auto slope = Exact::fromDouble(0.0);
                for (std::size_t degree = 1; degree < polynomial.size(); ++degree)
                    slope = Exact::add(slope,Exact::multiply(polynomial[degree],
                        Exact::fromDouble(static_cast<double>(degree))));
                Detail::require(Exact::sign(slope) != 0,
                    "internal multiple root at cell-section midpoint is unsupported");
            }
            const double tail = std::scalbn(root, -1);
            Detail::require(std::isfinite(tail) && std::scalbn(tail, 1) == root,
                "cell section root tail cannot be preserved by binary scaling");
            const CellChartCoordinate transverse{tail, side != 0};
            halves[side].push_back(coordinateIsX
                ? CellChartPoint{fixedCoordinate, transverse} : CellChartPoint{transverse, fixedCoordinate});
        }
    }
    std::vector<CellChartPoint> result = std::move(halves[0]);
    for (auto it = halves[1].rbegin(); it != halves[1].rend(); ++it)
    {
        if (!result.empty())
        {
            const auto previous = coordinateIsX ? result.back().y : result.back().x;
            const auto current = coordinateIsX ? it->y : it->x;
            const int order = compareCellChartCoordinate(previous, current);
            // The only shared identity between these disjoint half charts is
            // an exact root at their common midpoint, checked above at t=1.
            if (order == 0)
            {
                Detail::require(previous.offset == 0.5 && current.offset == 0.5
                    && !previous.fromUpper && current.fromUpper,
                    "distinct cell section roots have one represented chart point");
                continue;
            }
            Detail::require(order < 0, "cell section roots are not strictly ordered");
        }
        result.push_back(*it);
    }
    return result;
}

template<class ScalarSnapshot>
inline Frame BasicSnapshot<ScalarSnapshot>::jetAtCellChart(int i, int j, const CellChartPoint& chart) const
{
    namespace Exact = ExactExpansion;
    Detail::require(cellActive(i, j), "chart jet requires an active cell");
    const auto sample = scalar_.evaluateCellChart(i, j, chart);
    const long double gradient = std::hypot(static_cast<long double>(sample.dx), static_cast<long double>(sample.dy));
    Detail::require(usableGradientInCell(sample, gradient, i, j), "unresolved endpoint-relative interface gradient");
    const auto physical = [&](CellChartCoordinate component, int index, bool x)
    {
        if (component.offset == (component.fromUpper ? 1.0 : 0.0)) return coordinate(index, x);
        if (component.offset == (component.fromUpper ? 0.0 : 1.0)) return coordinate(index+1, x);
        const auto h = Exact::fromDouble(x ? h_.x : h_.y);
        auto displacement = Exact::multiply(h, Exact::fromDouble(component.offset));
        if (component.fromUpper) displacement = Exact::subtract(h, displacement);
        // Same nominal-h interior embedding as jetAtCellLocal, with the
        // actual coordinate(k+1) retained at an exact endpoint. This rounded
        // display point is not the coordinate/root authority.
        return Exact::estimate(Exact::add(Exact::fromDouble(coordinate(index, x)), displacement));
    };
    const Vector2D point(physical(chart.x, i, true), physical(chart.y, j, false));
    auto result = frameFromSample(point, sample, gradient);
    result.cellChartPoint = chart; result.cellI = i; result.cellJ = j; result.hasCellChartPoint = true;
    return result;
}

template<class ScalarSnapshot>
inline Frame BasicSnapshot<ScalarSnapshot>::jetInCell(int i, int j, Vector2D point) const
{
    const auto g = gridPoint(point);
    Detail::require(cellActive(i, j) && inClosedCell(point, i, j),
        "same-cell jet requires a point in the specified active closed cell");
    const auto s = sampleInCell(i, j, g);
    const long double gradient = std::hypot(static_cast<long double>(s.dx), static_cast<long double>(s.dy));
    Detail::require(usableGradientInCell(s, gradient, i, j),
        "zero, unresolved, or non-finite Hermite gradient cannot define an interface jet");
    return frameFromSample(point, s, gradient);
}

template<class ScalarSnapshot>
inline Frame BasicSnapshot<ScalarSnapshot>::frameFromSample(Vector2D point, const LevelSetHermiteScalar::Sample& s, long double gradient)
{
    return HermiteDifferential::frame<Frame>(point,s,gradient,Detail::checked);
}

template<class ScalarSnapshot>
inline double BasicSnapshot<ScalarSnapshot>::curvatureTangentialDerivative(Vector2D point) const
{
    return evaluateCurvatureDerivative(point, false);
}

template<class ScalarSnapshot>
inline double BasicSnapshot<ScalarSnapshot>::curvatureTangentialDerivativeAtCellLocal(int i, int j, Vector2D local) const
{
    Detail::require(cellActive(i, j), "local curvature derivative requires an active cell");
    const auto sample = scalar_.evaluateCell(i, j, local.x, local.y);
    const long double gradient = std::hypot(static_cast<long double>(sample.dx), static_cast<long double>(sample.dy));
    Detail::require(usableGradientInCell(sample, gradient, i, j), "unresolved cell-local interface gradient");
    // Same CellData as the local jet, including the declared sector. This
    // never asks the reconstructed-D3 path for a surrounding nodal stencil.
    const auto third = scalar_.thirdDerivativesCell(i, j, local.x, local.y);
    return curvatureDerivativeFromSample(sample, third, gradient);
}

template<class ScalarSnapshot>
inline double BasicSnapshot<ScalarSnapshot>::curvatureTangentialDerivativeAtCellChart(int i, int j, const CellChartPoint& chart) const
{
    Detail::require(cellActive(i, j), "chart curvature derivative requires an active cell");
    const auto sample = scalar_.evaluateCellChart(i, j, chart);
    const long double gradient = std::hypot(static_cast<long double>(sample.dx), static_cast<long double>(sample.dy));
    Detail::require(usableGradientInCell(sample, gradient, i, j), "unresolved endpoint-relative interface gradient");
    return curvatureDerivativeFromSample(sample, scalar_.thirdDerivativesCellChart(i, j, chart), gradient);
}

template<class ScalarSnapshot>
inline double BasicSnapshot<ScalarSnapshot>::reconstructedCurvatureTangentialDerivative(Vector2D point) const
{
    return evaluateCurvatureDerivative(point, true);
}

template<class ScalarSnapshot>
inline double BasicSnapshot<ScalarSnapshot>::evaluateCurvatureDerivative(Vector2D point, bool reconstructThird) const
{
    const auto s = sample(point);
    const long double gradient = std::hypot(static_cast<long double>(s.dx), static_cast<long double>(s.dy));
    Detail::require(usableGradient(s, gradient, point),
        "zero, unresolved, or non-finite Hermite gradient cannot define a curvature derivative");
    // Geometry's physical-to-local canonicalization is also used by sample
    // and jet. Do not independently ask scalar.thirdDerivatives(point) to
    // select a potentially different one-sided cell at a represented gridline.
    const auto g = gridPoint(point);
    const auto cell = selectCell(point, g);
    const int i = cell[0], j = cell[1];
    const auto third = reconstructThird
        ? scalar_.reconstructedThirdDerivativesCell(i, j, g[0] - i, g[1] - j)
        : scalar_.thirdDerivativesCell(i, j, g[0] - i, g[1] - j);
    return curvatureDerivativeFromSample(s, third, gradient);
}

template<class ScalarSnapshot>
inline double BasicSnapshot<ScalarSnapshot>::curvatureDerivativeFromSample(const LevelSetHermiteScalar::Sample& s,
    const LevelSetHermiteScalar::ThirdDerivatives& third, long double gradient)
{
    return HermiteDifferential::curvatureDerivative(s,third,gradient,Detail::checked);
}

template<class ScalarSnapshot>
inline bool BasicSnapshot<ScalarSnapshot>::project(Vector2D query, Vector2D seed, Frame& result, bool retainCellLocal) const
{
    // Two equations: phi(x)=0 and grad(phi) cross (x-query)=0. Row scaling
    // by |grad(phi)| gives two residuals with units of physical length. The
    // row scale is held fixed in each Newton linearization (no new equation).
    const long double length = std::max(h_.x, h_.y);
    const long double coordinateScale = std::max({length,
        std::abs(static_cast<long double>(lower_.x)), std::abs(static_cast<long double>(lower_.y)),
        std::abs(static_cast<long double>(upper_.x)), std::abs(static_cast<long double>(upper_.y))});
    const long double tolerance = 2048.0L * std::numeric_limits<double>::epsilon() * coordinateScale;
    Detail::require(tolerance <= std::sqrt(static_cast<long double>(std::numeric_limits<double>::epsilon()))
        * std::min(h_.x, h_.y),
        "coordinate precision cannot resolve the requested local Hermite projection; translate/rescale the input domain");
    Vector2D point = seed;
    for (int iteration = 0; iteration < 48; ++iteration)
    {
        const auto s = sample(point);
        const long double gx = s.dx, gy = s.dy, gradient = std::hypot(gx, gy);
        if (!usableGradient(s, gradient, point)) return false;
        const long double dx = static_cast<long double>(point.x) - query.x;
        const long double dy = static_cast<long double>(point.y) - query.y;
        const long double f = s.value / gradient, stationarity = (gx * dy - gy * dx) / gradient;
        const long double merit = std::hypot(f, stationarity);
        if (!std::isfinite(merit)) return false;
        if (std::abs(f) <= tolerance && std::abs(stationarity) <= tolerance)
        {
            if (restrictedDomain())
            {
                const auto cell = selectCell(point, gridPoint(point));
                const int i = cell[0], j = cell[1];
                const unsigned corners = cornerOnlyMasks_[static_cast<std::size_t>(j) * nx_ + i];
                if (corners != 0)
                {
                    // Small Newton residuals cannot contradict an exact
                    // single-phase certificate. Near a tangent zero corner,
                    // phi may be arbitrarily small on the wet extension but
                    // there is NO interior/edge zero whose jet can be used.
                    // Preserve exact corner candidates without snapping;
                    // the boundary entry point still applies its original
                    // endpoint KKT test and specified cap-side jet.
                    const int offsets[4][2] = {{0,0}, {1,0}, {1,1}, {0,1}};
                    bool atAllowedCorner = false;
                    for (unsigned k = 0; k < 4; ++k)
                        if ((corners & (1u << k)) != 0
                            && point.x == coordinate(i + offsets[k][0], true)
                            && point.y == coordinate(j + offsets[k][1], false))
                            atAllowedCorner = true;
                    if (!atAllowedCorner) return false;
                }
            }
            if (retainCellLocal)
            {
                // Preserve the SAME selected-cell sample parameters that
                // produced this Newton iterate's gradient and Hessian. Do
                // not attach a local flag to the old global-display Frame:
                // its physical-to-local-to-physical round trip need not be
                // bitwise exact on a nonbinary mesh or near a cell edge.
                const auto g = gridPoint(point);
                const auto cell = selectCell(point, g);
                result = jetAtCellLocal(cell[0], cell[1], {g[0]-cell[0], g[1]-cell[1]});
                const long double localDx = static_cast<long double>(result.point.x)-query.x;
                const long double localDy = static_cast<long double>(result.point.y)-query.y;
                const long double localStationarity = (gx*localDy-gy*localDx)/gradient;
                if (!std::isfinite(localStationarity) || result.levelSetResidual > tolerance
                    || std::abs(localStationarity) > tolerance) return false;
                result.stationarityResidual = Detail::checked(std::abs(localStationarity),
                    "non-finite cell-local projected stationarity residual");
            }
            else
            {
                result = jet(point);
                result.levelSetResidual = Detail::checked(std::abs(f), "non-finite projected scalar residual");
                result.stationarityResidual = Detail::checked(std::abs(stationarity), "non-finite projected stationarity residual");
            }
            result.iterations = iteration;
            return true;
        }
        const long double a = gx / gradient, b = gy / gradient;
        const long double c = (s.dxx * dy - s.dxy * dx - gy) / gradient;
        const long double d = (s.dxy * dy - s.dyy * dx + gx) / gradient;
        const long double determinant = a * d - b * c;
        const long double matrixScale = std::hypot(a, b) * std::hypot(c, d);
        if (!std::isfinite(determinant) || !std::isfinite(matrixScale)
            || std::abs(determinant) <= 512.0L * std::numeric_limits<double>::epsilon() * matrixScale) return false;
        long double stepX = (-d * f + b * stationarity) / determinant;
        long double stepY = (c * f - a * stationarity) / determinant;
        const long double stepLength = std::hypot(stepX, stepY);
        if (!std::isfinite(stepLength)) return false;
        if (stepLength > 2.0L * length)
        {
            stepX *= 2.0L * length / stepLength;
            stepY *= 2.0L * length / stepLength;
        }
        bool accepted = false;
        for (int backtrack = 0; backtrack < 24; ++backtrack)
        {
            const long double factor = std::ldexp(1.0L, -backtrack);
            const long double trialX = point.x + factor * stepX, trialY = point.y + factor * stepY;
            if (!std::isfinite(trialX) || !std::isfinite(trialY)) continue;
            const Vector2D trial(static_cast<double>(trialX), static_cast<double>(trialY));
            if (!inDomain(trial) || (trial.x == point.x && trial.y == point.y)) continue;
            const auto v = sample(trial);
            const long double trialGradient = std::hypot(static_cast<long double>(v.dx), static_cast<long double>(v.dy));
            if (!usableGradient(v, trialGradient, trial)) continue;
            const long double fTrial = v.value / trialGradient;
            const long double sTrial = (v.dx * (static_cast<long double>(trial.y) - query.y)
                - v.dy * (static_cast<long double>(trial.x) - query.x)) / trialGradient;
            const long double trialMerit = std::hypot(fTrial, sTrial);
            if (std::isfinite(trialMerit) && trialMerit < merit)
            { point = trial; accepted = true; break; }
        }
        if (!accepted) return false;
    }
    return false;
}

template<class ScalarSnapshot>
inline Frame BasicSnapshot<ScalarSnapshot>::projectedFrame(Vector2D query) const
{
    selectCell(query, gridPoint(query));
    Detail::require(!seeds_.empty(), "no resolved Hermite zero set exists for projectedFrame");
    Frame best;
    const bool converged = stationaryFrame(query, best);
    Detail::require(converged, "local multistart Hermite projection did not converge; no fallback is provided");
    return best;
}

template<class ScalarSnapshot>
inline Frame BasicSnapshot<ScalarSnapshot>::projectedCellLocalFrame(Vector2D query) const
{
    selectCell(query, gridPoint(query));
    Detail::require(!seeds_.empty(), "no resolved Hermite zero set exists for projectedCellLocalFrame");
    Frame best;
    const bool converged = stationaryFrame(query, best, true);
    Detail::require(converged, "local multistart Hermite projection did not converge with cell-local provenance; no endpoint fallback is provided");
    return best;
}

template<class ScalarSnapshot>
inline std::vector<Frame> BasicSnapshot<ScalarSnapshot>::stationaryCellLocalFrames(Vector2D query, Frame& diagnostics) const
{
    selectCell(query, gridPoint(query));
    Detail::require(!seeds_.empty(), "no resolved Hermite zero set exists for local stationary candidates");
    std::vector<Frame> candidates;
    stationaryFrame(query, diagnostics, true, &candidates);
    return candidates;
}

template<class ScalarSnapshot>
inline bool BasicSnapshot<ScalarSnapshot>::stationaryFrame(Vector2D query, Frame& best, bool retainCellLocal,
    std::vector<Frame>* candidates) const
{
    // Eight is a fixed algorithmic multistart budget, not a shape-dependent
    // or measured-error-dependent switch. Searching all roots for the eight
    // closest seeds is cheap; Newton is NOT run for every contour root.
    std::vector<std::pair<long double, std::size_t>> nearest;
    const long double length = std::max(h_.x, h_.y);
    for (std::size_t k = 0; k < seeds_.size(); ++k)
    {
        const long double distance = std::hypot((static_cast<long double>(seeds_[k].x) - query.x) / length,
            (static_cast<long double>(seeds_[k].y) - query.y) / length);
        Detail::require(std::isfinite(distance), "unrepresentable projection seed distance");
        const auto position = std::lower_bound(nearest.begin(), nearest.end(), std::make_pair(distance, k));
        if (position != nearest.end() || nearest.size() < 8)
        {
            nearest.insert(position, {distance, k});
            if (nearest.size() > 8) nearest.pop_back();
        }
    }
    long double bestDistance = std::numeric_limits<long double>::infinity();
    int converged = 0;
    for (const auto& entry : nearest)
    {
        Frame candidate;
        if (!project(query, seeds_[entry.second], candidate, retainCellLocal)) continue;
        ++converged;
        if (candidates) candidates->push_back(candidate);
        const long double distance = std::hypot((static_cast<long double>(candidate.point.x) - query.x) / length,
            (static_cast<long double>(candidate.point.y) - query.y) / length);
        if (distance < bestDistance) { bestDistance = distance; best = candidate; }
    }
    best.attemptedSeeds = static_cast<int>(nearest.size());
    best.convergedSeeds = converged;
    if (candidates) for (auto& candidate : *candidates)
    {
        candidate.attemptedSeeds = best.attemptedSeeds;
        candidate.convergedSeeds = best.convergedSeeds;
    }
    return converged > 0;
}

template<class ScalarSnapshot>
inline BoundaryProjection BasicSnapshot<ScalarSnapshot>::boundaryProjectedFrame(Vector2D query,
    const std::vector<Vector2D>& boundaryEndpoints) const
{
    if (boundaryEndpoints.empty())
    {
        const auto frame = projectedFrame(query);
        return {frame, false, 0.0, frame.stationarityResidual};
    }
    selectCell(query, gridPoint(query));
    const long double length = std::max(h_.x, h_.y);
    const long double coordinateScale = std::max({length,
        std::abs(static_cast<long double>(lower_.x)), std::abs(static_cast<long double>(lower_.y)),
        std::abs(static_cast<long double>(upper_.x)), std::abs(static_cast<long double>(upper_.y))});
    const long double tolerance = 2048.0L * std::numeric_limits<double>::epsilon() * coordinateScale;
    Detail::require(tolerance <= std::sqrt(static_cast<long double>(std::numeric_limits<double>::epsilon()))
        * std::min(h_.x, h_.y), "coordinate precision cannot resolve the requested boundary projection");
    std::vector<Endpoint> endpoints;
    endpoints.reserve(boundaryEndpoints.size());
    namespace Exact = ExactExpansion;
    for (const auto point : boundaryEndpoints)
    {
        const auto g = gridPoint(point);
        const bool west = point.x == lower_.x, east = point.x == upper_.x;
        const bool south = point.y == lower_.y, north = point.y == upper_.y;
        Detail::require(static_cast<int>(west) + east + south + north == 1,
            "projection endpoint must lie on exactly one noncorner physical wall");
        for (const auto& old : endpoints)
            Detail::require(point.x != old.frame.point.x || point.y != old.frame.point.y,
                "duplicate declared projection endpoint");
        const bool vertical = west || east;
        const auto cell = selectCell(point, g);
        const int i = cell[0], j = cell[1];
        const Vector2D localA(vertical ? g[0]-i : 0.0, vertical ? 0.0 : g[1]-j);
        const Vector2D localB(vertical ? g[0]-i : 1.0, vertical ? 1.0 : g[1]-j);
        const auto polynomial = scalar_.segmentPolynomialExact(i, j, localA, localB);
        // An exact zero takes the fast path. A general binary64 root position
        // is accepted only by identity with this scalar's canonical full-edge
        // root inventory, not by an epsilon test on its sampled value.
        if (Exact::sign(exactValue(point)) != 0)
        {
            const int x = vertical ? (west ? 0 : nx_) : i;
            const int y = vertical ? j : (south ? 0 : ny_);
            const Vector2D a(coordinate(x, true), coordinate(y, false));
            const Vector2D b(coordinate(x + (vertical ? 0 : 1), true),
                coordinate(y + (vertical ? 1 : 0), false));
            const auto zeros = intersectGridEdge(a, b);
            bool matched = false;
            for (const auto& zero : zeros.points) matched = matched || (zero.x == point.x && zero.y == point.y);
            Detail::require(!zeros.coincident && matched,
                "declared projection endpoint is not an authoritative wall zero");
        }
        ExpandedPolynomial derivative{};
        for (std::size_t k = 1; k < polynomial.size(); ++k)
            derivative[k-1] = Exact::multiply(polynomial[k], Exact::fromDouble(static_cast<double>(k)));
        const double parameter = vertical ? g[1]-j : g[0]-i;
        Detail::require(Exact::sign(Exact::evaluate(derivative, parameter)) != 0,
            "wall-coincident or tangential zero cannot define a one-sided projection endpoint");
        auto frame = jet(point);
        Detail::require(frame.levelSetResidual <= tolerance, "wall-root representation has an unresolved scalar residual");
        const Vector2D tangent(-frame.normal.y, frame.normal.x);
        const double intoDomain = west ? tangent.x : east ? -tangent.x : south ? tangent.y : -tangent.y;
        Detail::require(std::isfinite(intoDomain)
            && std::abs(intoDomain) > 512.0 * std::numeric_limits<double>::epsilon(),
            "wall crossing is too tangential to resolve an inward curve direction");
        endpoints.push_back({frame, intoDomain > 0.0 ? tangent : tangent * -1.0});
    }
    return selectBoundaryProjection(query, endpoints, length, tolerance);
}

template<class ScalarSnapshot>
inline BoundaryProjection BasicSnapshot<ScalarSnapshot>::patchBoundaryProjectedFrame(Vector2D query,
    const std::vector<PatchEndpoint>& boundaryEndpoints) const
{
    Detail::require(restrictedDomain(), "patch endpoint projection requires an explicitly restricted scalar");
    if (boundaryEndpoints.empty())
    {
        const auto frame = projectedFrame(query);
        return {frame, false, 0.0, frame.stationarityResidual};
    }
    selectCell(query, gridPoint(query));
    const long double length = std::max(h_.x, h_.y);
    const long double coordinateScale = std::max({length,
        std::abs(static_cast<long double>(lower_.x)), std::abs(static_cast<long double>(lower_.y)),
        std::abs(static_cast<long double>(upper_.x)), std::abs(static_cast<long double>(upper_.y))});
    const long double tolerance = 2048.0L * std::numeric_limits<double>::epsilon() * coordinateScale;
    Detail::require(tolerance <= std::sqrt(static_cast<long double>(std::numeric_limits<double>::epsilon()))
        * std::min(h_.x, h_.y), "coordinate precision cannot resolve the requested boundary projection");
    std::vector<Endpoint> endpoints;
    endpoints.reserve(boundaryEndpoints.size());
    for (const auto& endpoint : boundaryEndpoints)
    {
        const auto point = endpoint.point;
        const int i = endpoint.cellI, j = endpoint.cellJ;
        Detail::require(endpoint.inwardTangentSign == -1 || endpoint.inwardTangentSign == 1,
            "patch endpoint inward tangent sign must be exactly -1 or +1");
        auto frame = jetInCell(i, j, point);
        const auto g = gridPoint(point);
        for (const auto& old : endpoints)
            Detail::require(point.x != old.frame.point.x || point.y != old.frame.point.y,
                "duplicate declared projection endpoint");
        const bool exactZero = ExactExpansion::sign(
            scalar_.evaluateExactCell(i, j, g[0]-i, g[1]-j)) == 0;
        bool onBoundary = false, canonicalRoot = false;
        const auto active = [&](int x, int y)
        { return x >= 0 && x < nx_ && y >= 0 && y < ny_ && cellActive(x, y); };
        const auto inspectEdge = [&](int x, int y, bool vertical)
        {
            const auto a = Vector2D(coordinate(x, true), coordinate(y, false));
            const auto b = Vector2D(coordinate(x + (vertical ? 0 : 1), true),
                coordinate(y + (vertical ? 1 : 0), false));
            if (point.x < a.x || point.x > b.x || point.y < a.y || point.y > b.y) return;
            const bool low = vertical ? active(x-1, y) : active(x, y-1);
            const bool high = active(x, y);
            if (low == high) return;
            onBoundary = true;
            if (!exactZero)
            {
                const auto zeros = intersectGridEdge(a, b);
                if (!zeros.coincident) for (const auto& zero : zeros.points)
                    canonicalRoot = canonicalRoot || (zero.x == point.x && zero.y == point.y);
            }
        };
        // A concave active-domain corner may touch the supplied cap cell
        // only at its vertex. Inspect every incident boundary edge, not just
        // that cell's own sides. Exact physical equality alone selects them.
        const int x = static_cast<int>(std::llround(g[0]));
        const int y = static_cast<int>(std::llround(g[1]));
        if (point.x == coordinate(x, true))
        {
            const int first = std::min(static_cast<int>(g[1]), ny_-1);
            inspectEdge(x, first, true);
            if (y > 0 && y < ny_ && point.y == coordinate(y, false) && first != y-1)
                inspectEdge(x, y-1, true);
        }
        if (point.y == coordinate(y, false))
        {
            const int first = std::min(static_cast<int>(g[0]), nx_-1);
            inspectEdge(first, y, false);
            if (x > 0 && x < nx_ && point.x == coordinate(x, true) && first != x-1)
                inspectEdge(x-1, y, false);
        }
        Detail::require(onBoundary, "patch endpoint is not on an actual active-domain boundary edge");
        Detail::require(exactZero || canonicalRoot,
            "declared patch endpoint is not an authoritative boundary zero");
        Detail::require(frame.levelSetResidual <= tolerance, "wall-root representation has an unresolved scalar residual");
        const Vector2D inward = Vector2D(-frame.normal.y, frame.normal.x)
            * static_cast<double>(endpoint.inwardTangentSign);
        Detail::require((point.x != coordinate(i, true) || inward.x >= 0.0)
            && (point.x != coordinate(i+1, true) || inward.x <= 0.0)
            && (point.y != coordinate(j, false) || inward.y >= 0.0)
            && (point.y != coordinate(j+1, false) || inward.y <= 0.0),
            "patch endpoint tangent points outside its specified incident cell");
        // A zero component is valid (e.g. the initial lip or an axis pole).
        // The supplied, validated panel selects the curve's incident cell;
        // no wall-normal clipping or analytic tangent is manufactured here.
        endpoints.push_back({frame, inward});
    }
    return selectBoundaryProjection(query, endpoints, length, tolerance);
}

template<class ScalarSnapshot>
inline BoundaryProjection BasicSnapshot<ScalarSnapshot>::selectBoundaryProjection(Vector2D query,
    const std::vector<Endpoint>& endpoints, long double length, long double tolerance) const
{
    Detail::require(!seeds_.empty(), "no resolved Hermite zero set exists for projectedFrame");
    Frame stationary;
    bool found = stationaryFrame(query, stationary);
    BoundaryProjection best{stationary, false, 0.0, stationary.stationarityResidual};
    long double bestDistance = found ? std::hypot((static_cast<long double>(stationary.point.x)-query.x)/length,
        (static_cast<long double>(stationary.point.y)-query.y)/length) : std::numeric_limits<long double>::infinity();
    for (auto candidate : endpoints)
    {
        const long double dx = static_cast<long double>(candidate.frame.point.x)-query.x;
        const long double dy = static_cast<long double>(candidate.frame.point.y)-query.y;
        const long double inward = dx*candidate.inward.x + dy*candidate.inward.y;
        Detail::require(std::isfinite(inward), "nonfinite one-sided projection derivative");
        if (inward < -tolerance) continue; // Valid endpoint, but not a candidate for THIS query.
        const long double distance = std::hypot(dx/length, dy/length);
        Detail::require(std::isfinite(distance), "unrepresentable endpoint projection distance");
        const bool samePoint = candidate.frame.point.x == best.frame.point.x && candidate.frame.point.y == best.frame.point.y;
        if (!found || distance < bestDistance || (distance == bestDistance && samePoint))
        {
            candidate.frame.stationarityResidual = Detail::checked(std::abs(inward), "nonfinite endpoint stationarity residual");
            candidate.frame.attemptedSeeds = stationary.attemptedSeeds;
            candidate.frame.convergedSeeds = stationary.convergedSeeds;
            best = {candidate.frame, true, Detail::checked(inward, "nonfinite endpoint inward derivative"),
                Detail::checked(std::max(0.0L, -inward), "nonfinite endpoint KKT residual")};
            bestDistance = distance;
            found = true;
        }
    }
    Detail::require(found, "local boundary projection has no converged stationary or admissible endpoint candidate");
    return best;
}

template<class ScalarSnapshot>
inline bool BasicSnapshot<ScalarSnapshot>::segmentInsideActiveDomain(Vector2D a, Vector2D b) const
{
    return segmentCoveredByCells(a, b, 0);
}

template<class ScalarSnapshot>
inline int BasicSnapshot<ScalarSnapshot>::certifiedSegmentSign(Vector2D a, Vector2D b) const
{
    if (segmentCoveredByCells(a, b, -1)) return -1;
    if (segmentCoveredByCells(a, b, +1)) return +1;
    return 0;
}

template<class ScalarSnapshot>
inline bool BasicSnapshot<ScalarSnapshot>::segmentCoveredByCells(Vector2D a, Vector2D b, int requiredSign) const
{
    const auto ga = gridPoint(a), gb = gridPoint(b); // Finite physical rectangle / FP guards only.
    if (!requiredSign && !restrictedDomain()) return true; // A closed rectangle is convex.
    if (!requiredSign && a.x == b.x && a.y == b.y) return findCell(a, ga)[0] >= 0;

    namespace Exact = ExactExpansion;
    struct Rational { Exact::Expansion numerator, denominator; }; // Denominator strictly positive.
    struct Interval { Rational lower, upper; };
    const auto zero = Exact::fromDouble(0.0), one = Exact::fromDouble(1.0);
    const Rational begin{zero, one}, end{one, one};
    const auto compare = [](const Rational& x, const Rational& y) {
        return Exact::sign(Exact::subtract(Exact::multiply(x.numerator, y.denominator),
            Exact::multiply(y.numerator, x.denominator)));
    };
    const std::array<double, 2> av{{a.x, a.y}}, bv{{b.x, b.y}};
    const std::array<Exact::Expansion, 2> origin{{Exact::fromDouble(a.x), Exact::fromDouble(a.y)}};
    const std::array<Exact::Expansion, 2> finish{{Exact::fromDouble(b.x), Exact::fromDouble(b.y)}};
    std::array<Exact::Expansion, 2> length;
    std::array<int, 2> direction;
    for (int axis = 0; axis < 2; ++axis)
    {
        direction[axis] = (bv[axis] > av[axis]) - (bv[axis] < av[axis]);
        length[axis] = direction[axis] >= 0 ? Exact::subtract(finish[axis], origin[axis])
                                          : Exact::subtract(origin[axis], finish[axis]);
    }
    // Normalized coordinates restrict the search only. Each candidate is
    // clipped against the unchanged canonical binary64 cell coordinates.
    const int imin = std::max(0, static_cast<int>(std::min(ga[0], gb[0])) - 1);
    const int imax = std::min(nx_ - 1, static_cast<int>(std::max(ga[0], gb[0])) + 1);
    const int jmin = std::max(0, static_cast<int>(std::min(ga[1], gb[1])) - 1);
    const int jmax = std::min(ny_ - 1, static_cast<int>(std::max(ga[1], gb[1])) + 1);
    std::vector<Interval> intervals;
    for (int j = jmin; j <= jmax; ++j) for (int i = imin; i <= imax; ++i)
    {
        if (!cellActive(i, j)) continue;
        if (requiredSign && cellSigns_[static_cast<std::size_t>(j)*nx_+i] != 0
            && cellSigns_[static_cast<std::size_t>(j)*nx_+i] != requiredSign) continue;
        if (requiredSign && cellSigns_[static_cast<std::size_t>(j)*nx_+i] != requiredSign)
        {
            // A cut cell can still certify a single-phase half-face far from
            // its zero curve. Enclose the segment/cell intersection by a box;
            // this does not construct a rounded grid-cut point or midpoint.
            // The range is widened by outward coordinate arithmetic and by
            // exact side identities before intersecting the scalar unit cell.
            using Arithmetic = LevelSetHermiteScalar::Detail::Outward;
            Vector2D localLower, localUpper;
            bool boxExists = true;
            for (int axis = 0; axis < 2; ++axis)
            {
                const bool x = axis == 0;
                const int index = x ? i : j;
                const double first = coordinate(index, x), last = coordinate(index+1, x);
                const double lo = std::max(std::min(av[axis],bv[axis]), first);
                const double hi = std::min(std::max(av[axis],bv[axis]), last);
                if (lo > hi) { boxExists = false; break; }
                const auto normalized = [&](double p) {
                    auto range = Arithmetic::add(Arithmetic::point(p),
                        Arithmetic::negate(Arithmetic::point(x ? lower_.x : lower_.y)));
                    range = Arithmetic::dividePositive(range, x ? h_.x : h_.y);
                    range = Arithmetic::add(range, Arithmetic::point(-static_cast<double>(index)));
                    // gridPoint uses exact physical side identities as well
                    // as division. Include both representations in the range.
                    if (p == first) { range.lo = std::min(range.lo,0.0); range.hi = std::max(range.hi,0.0); }
                    if (p == last) { range.lo = std::min(range.lo,1.0); range.hi = std::max(range.hi,1.0); }
                    return range;
                };
                const auto aRange = normalized(lo), bRange = normalized(hi);
                if (!Arithmetic::finite(aRange) || !Arithmetic::finite(bRange)) { boxExists = false; break; }
                const double aLocal = std::max(0.0,aRange.lo), bLocal = std::min(1.0,bRange.hi);
                if (aLocal > bLocal) { boxExists = false; break; }
                (x ? localLower.x : localLower.y) = aLocal;
                (x ? localUpper.x : localUpper.y) = bLocal;
            }
            if (!boxExists) continue;
            const auto bound = scalar_.cellBounds(i,j,localLower,localUpper);
            const int sign = bound.lo > 0.0 ? 1 : bound.hi < 0.0 ? -1 : 0;
            if (sign != requiredSign) continue;
        }
        if (inClosedCell(a, i, j) && inClosedCell(b, i, j)) return true;
        const std::array<double, 2> lower{{coordinate(i, true), coordinate(j, false)}};
        const std::array<double, 2> upper{{coordinate(i + 1, true), coordinate(j + 1, false)}};
        Interval interval{begin, end};
        bool intersects = true;
        for (int axis = 0; axis < 2 && intersects; ++axis)
        {
            if (!direction[axis])
            {
                intersects = av[axis] >= lower[axis] && av[axis] <= upper[axis];
                continue;
            }
            const auto lo = Exact::fromDouble(lower[axis]), hi = Exact::fromDouble(upper[axis]);
            const Rational enter{direction[axis] > 0 ? Exact::subtract(lo, origin[axis])
                                                     : Exact::subtract(origin[axis], hi), length[axis]};
            const Rational leave{direction[axis] > 0 ? Exact::subtract(hi, origin[axis])
                                                     : Exact::subtract(origin[axis], lo), length[axis]};
            if (compare(enter, interval.lower) > 0) interval.lower = enter;
            if (compare(leave, interval.upper) < 0) interval.upper = leave;
            intersects = compare(interval.lower, interval.upper) <= 0;
        }
        if (intersects) intervals.push_back(std::move(interval));
    }
    std::sort(intervals.begin(), intervals.end(), [&](const Interval& x, const Interval& y) {
        const int first = compare(x.lower, y.lower);
        return first < 0 || (first == 0 && compare(x.upper, y.upper) < 0);
    });
    Rational covered = begin;
    for (const auto& interval : intervals)
    {
        if (compare(interval.lower, covered) > 0) return false; // Even a nonrepresentable gap counts.
        if (compare(interval.upper, covered) > 0) covered = interval.upper;
        if (compare(covered, end) == 0) return true;
    }
    return false;
}

template<class ScalarSnapshot>
inline std::vector<Piece> BasicSnapshot<ScalarSnapshot>::partitionSegment(Vector2D a, Vector2D b) const
{
    return partitionSegmentImpl(a, b, nullptr);
}

template<class ScalarSnapshot>
inline std::vector<ParametricPiece> BasicSnapshot<ScalarSnapshot>::partitionSegmentParametric(Vector2D a, Vector2D b, bool certifyGridSlivers) const
{
    std::vector<ParametricPiece> result;
    partitionSegmentImpl(a, b, &result, certifyGridSlivers);
    return result;
}

template<class ScalarSnapshot>
inline std::vector<SegmentRootIncidence> BasicSnapshot<ScalarSnapshot>::intersectSegmentWithProvenance(
    Vector2D a, Vector2D b, bool certifyGridSlivers) const
{
    std::vector<SegmentRootIncidence> result;
    partitionSegmentImpl(a, b, nullptr, certifyGridSlivers, &result);
    return result;
}

template<class ScalarSnapshot>
inline std::vector<Piece> BasicSnapshot<ScalarSnapshot>::partitionSegmentImpl(Vector2D a, Vector2D b,
    std::vector<ParametricPiece>* parametric, bool certifyGridSlivers,
    std::vector<SegmentRootIncidence>* rootIncidences) const
{
    selectCell(a, gridPoint(a)); selectCell(b, gridPoint(b));
    Detail::require(a.x != b.x || a.y != b.y, "zero-length segment has no phase fraction");
    // Evaluate each unoriented segment exactly once. Reversal only reverses
    // the returned pieces; it does not invoke another floating root search.
    if (std::make_pair(b.x, b.y) < std::make_pair(a.x, a.y))
    {
        auto reverse = partitionSegmentImpl(b, a, parametric, certifyGridSlivers, rootIncidences);
        std::reverse(reverse.begin(), reverse.end());
        for (auto& piece : reverse) std::swap(piece.begin, piece.end);
        if (parametric)
        {
            // Retain the same chart: 1-t would erase a tiny endpoint interval.
            std::reverse(parametric->begin(), parametric->end());
            for (auto& piece : *parametric) std::swap(piece.tBegin, piece.tEnd);
        }
        if (rootIncidences) std::reverse(rootIncidences->begin(), rootIncidences->end());
        return reverse;
    }
    struct Cut { double t; int x = -1, y = -1; };
    std::vector<Cut> gridCuts{{0.0, -1, -1}, {1.0, -1, -1}};
    for (int axis = 0; axis < 2; ++axis)
    {
        const bool x = axis == 0;
        const long double begin = x ? a.x : a.y, end = x ? b.x : b.y;
        if (begin == end) continue;
        for (int k = 1; k < (x ? nx_ : ny_); ++k)
        {
            const long double line = coordinate(k, x);
            if (!(line > std::min(begin, end) && line < std::max(begin, end))) continue;
            const double t = Detail::checked((line - begin) / (end - begin), "unrepresentable Cartesian segment split");
            Detail::require(t > 0.0 && t < 1.0, "Cartesian segment split collapsed to an endpoint");
            gridCuts.push_back({t, x ? k : -1, x ? -1 : k});
        }
    }
    std::sort(gridCuts.begin(), gridCuts.end(), [](const Cut& x, const Cut& y) { return x.t < y.t; });
    std::vector<Cut> cuts;
    for (const auto& cut : gridCuts)
    {
        if (!cuts.empty() && cuts.back().t == cut.t)
        {
            Detail::require((cut.x < 0 || cuts.back().x < 0 || cut.x == cuts.back().x)
                && (cut.y < 0 || cuts.back().y < 0 || cut.y == cuts.back().y),
                "distinct Cartesian crossings collapse to the same segment parameter");
            if (cut.x >= 0) cuts.back().x = cut.x;
            if (cut.y >= 0) cuts.back().y = cut.y;
            if (cuts.back().x >= 0 && cuts.back().y >= 0)
            {
                const auto direct = Detail::between(a, b, cut.t);
                Detail::require(direct.x == coordinate(cuts.back().x, true)
                    && direct.y == coordinate(cuts.back().y, false),
                    "cannot distinguish near-coincident Cartesian crossings without removing a subsegment");
            }
        }
        else cuts.push_back(cut);
    }
    auto pointAt = [&](const Cut& cut)
    {
        auto point = Detail::between(a, b, cut.t);
        if (cut.x >= 0) point.x = coordinate(cut.x, true);
        if (cut.y >= 0) point.y = coordinate(cut.y, false);
        return point;
    };
    std::vector<Piece> pieces;
    const auto appendPiece = [&](Vector2D first, Vector2D last, bool isNegative)
    {
        if (!pieces.empty() && pieces.back().negative == isNegative)
            pieces.back().end = last;
        else pieces.push_back({first, last, isNegative});
    };
    for (std::size_t k = 1; k < cuts.size(); ++k)
    {
        const auto begin = pointAt(cuts[k - 1]), end = pointAt(cuts[k]);
        const auto middle = Detail::between(begin, end, 0.5);
        const bool hasInterior = (middle.x != begin.x || middle.y != begin.y)
            && (middle.x != end.x || middle.y != end.y);
        if (!hasInterior && rootIncidences)
            Detail::require(begin.x != end.x || begin.y != end.y,
                "root-incidence grid piece has collapsed physical anchors");
        if (!hasInterior && parametric && certifyGridSlivers)
        {
            Detail::require(begin.x != end.x || begin.y != end.y,
                "unrepresentable Cartesian grid piece has collapsed anchors");
            const int sign = certifiedSegmentSign(begin, end);
            Detail::require(sign != 0, "unrepresentable Cartesian grid piece lacks an entire closed-phase certificate");
            // Keep its real nonzero anchor interval; no root/phase averaging
            // or global endpoint deletion is permitted by this policy.
            if (parametric) parametric->push_back({begin, end, 0.0, 1.0, sign < 0});
            continue;
        }
        // Root incidences use an explicit cell-local polynomial, not a
        // physical midpoint phase sample. A globally sub-ULP grid piece may
        // still have distinct local endpoints and an exact endpoint zero.
        // Physical/phase-piece APIs retain their original stricter policy.
        Detail::require(hasInterior || (rootIncidences && certifyGridSlivers),
            "unrepresentable Cartesian subsegment interior");
        const auto g = gridPoint(middle);
        auto cell = selectCell(middle, g);
        if (rootIncidences && !hasInterior
            && !(this->inClosedCell(begin, cell[0], cell[1]) && this->inClosedCell(end, cell[0], cell[1])))
        {
            // The rounded midpoint may select the cell on the OTHER side of
            // a gridline. Select only a real active closed cell containing
            // BOTH represented anchors; never extend its physical bounds.
            bool found = false;
            for (int j = std::max(0, cell[1]-1); j <= std::min(ny_-1, cell[1]+1) && !found; ++j)
                for (int i = std::max(0, cell[0]-1); i <= std::min(nx_-1, cell[0]+1); ++i)
                    if (cellActive(i, j) && this->inClosedCell(begin, i, j) && this->inClosedCell(end, i, j))
                    { cell = {i, j}; found = true; break; }
            Detail::require(found, "root-incidence grid piece has no active cell containing both anchors");
        }
        const int i = cell[0], j = cell[1];
        // Even if a rounded midpoint lands exactly on an active boundary,
        // both physical endpoints must bound a subsegment in that SAME
        // active closed cell. Never clamp away positive inactive length.
        if (restrictedDomain())
            Detail::require(inClosedCell(begin, i, j) && inClosedCell(end, i, j),
                "positive-length segment is not contained in an active cell; inactive traversal is unsupported");
        const int cellSign = cellSigns_[static_cast<std::size_t>(j) * nx_ + i];
        const auto globalBegin = gridPoint(begin), globalEnd = gridPoint(end);
        const auto inClosedCell = [&](const std::array<double, 2>& point)
        {
            return point[0] >= i && point[0] <= i + 1
                && point[1] >= j && point[1] <= j + 1;
        };
        if (!rootIncidences && cellSign != 0 && inClosedCell(globalBegin) && inClosedCell(globalEnd))
        {
            // Convex-hull bounds of authoritative Bernstein controls prove
            // no root anywhere in this closed cell. This is not an endpoint
            // sign heuristic and cannot discard a hidden opposite phase.
            // A rounded midpoint alone does not establish cell containment:
            // both represented endpoints must be in the certified cell,
            // without expanding its domain by the coordinate-roundoff guard.
            if (parametric) parametric->push_back({begin, end, 0.0, 1.0, cellSign < 0});
            else appendPiece(begin, end, cellSign < 0);
            continue;
        }
        auto localInCell = [&](Vector2D point, int cellI, int cellJ)
        {
            // Use the SAME global coordinate conversion as point queries.
            // Subtracting a separately rounded physical cell origin can give
            // a different local coordinate exactly where the phase is zero.
            const auto global = gridPoint(point);
            Vector2D q;
            for (int axis = 0; axis < 2; ++axis)
            {
                const bool x = axis == 0;
                const int cell = x ? cellI : cellJ;
                const double normalized = global[axis] - cell;
                const double guard = 128.0 * std::numeric_limits<double>::epsilon() * std::max(nx_, ny_);
                Detail::require(normalized >= -guard && normalized <= 1.0 + guard,
                    "segment endpoint is not in its selected Cartesian cell");
                // Strict physical-domain validation precedes this conversion.
                // This only bounds arithmetic coordinate roundoff at cell sides.
                (x ? q.x : q.y) = std::clamp(normalized, 0.0, 1.0);
            }
            return q;
        };
        const auto local = [&](Vector2D point) { return localInCell(point, i, j); };
        if (rootIncidences)
        {
            // Include each actual closed-cell incidence, not merely the
            // canonical midpoint owner. A curve/panel may belong to the other
            // active side of a grid edge, or touch this segment at one corner.
            std::vector<SegmentRootIncidence> partRoots;
            for (int cellJ = std::max(0, j-1); cellJ <= std::min(ny_-1, j+1); ++cellJ)
                for (int cellI = std::max(0, i-1); cellI <= std::min(nx_-1, i+1); ++cellI)
                {
                    if (!cellActive(cellI, cellJ)) continue;
                    const bool containsBegin = this->inClosedCell(begin, cellI, cellJ);
                    const bool containsEnd = this->inClosedCell(end, cellI, cellJ);
                    if (!containsBegin && !containsEnd) continue;
                    // A rigorous strict whole-cell sign cannot contain even
                    // an endpoint root. This is not a sample-sign filter.
                    if (cellSigns_[static_cast<std::size_t>(cellJ)*nx_+cellI] != 0) continue;
                    if (!(containsBegin && containsEnd))
                    {
                        const double parameter = containsBegin ? 0.0 : 1.0;
                        const auto q = localInCell(containsBegin ? begin : end, cellI, cellJ);
                        if (ExactExpansion::sign(scalar_.evaluateExactCell(cellI, cellJ, q.x, q.y)) == 0)
                            partRoots.push_back({cellI, cellJ, begin, end, q, q, q, parameter, true});
                        continue;
                    }
                    const auto localBegin = localInCell(begin, cellI, cellJ);
                    const auto localEnd = localInCell(end, cellI, cellJ);
                    Detail::require(localBegin.x != localEnd.x || localBegin.y != localEnd.y,
                        "nonzero physical subsegment has collapsed cell-local endpoints");
                    const auto polynomial = scalar_.segmentPolynomialExact(cellI, cellJ, localBegin, localEnd);
                    const auto roots = PolynomialSegmentRoots::isolateEndpointRelative(polynomial);
                    Detail::require(!roots.identicallyZero,
                        "coincident physical-segment zero restriction is unsupported");
                    Detail::require(ExactExpansion::sign(ExactExpansion::evaluate(polynomial, 0.0))
                            == ExactExpansion::sign(exactValue(begin))
                        && ExactExpansion::sign(ExactExpansion::evaluate(polynomial, 1.0))
                            == ExactExpansion::sign(exactValue(end)),
                        "root incidence restriction and canonical endpoint signs disagree");
                    Detail::require(roots.roots.size() == roots.enclosures.size(),
                        "root incidence representatives and enclosures disagree in count");
                    Vector2D previousPoint;
                    for (std::size_t rootIndex = 0; rootIndex < roots.roots.size(); ++rootIndex)
                    {
                        const double parameter = roots.roots[rootIndex];
                        const auto bracket = roots.enclosures[rootIndex];
                        Detail::require(std::isfinite(parameter) && parameter >= 0.0 && parameter <= 1.0
                            && bracket.lo >= 0.0 && bracket.hi <= 1.0
                            && bracket.lo <= parameter && parameter <= bracket.hi
                            && (rootIndex == 0 || (roots.roots[rootIndex-1] < parameter
                                && roots.enclosures[rootIndex-1].hi <= bracket.lo)),
                            "root incidence isolator returned unsorted or inconsistent enclosures");
                        const auto q = Detail::between(localBegin, localEnd, parameter);
                        Detail::require((parameter == 0.0 || q.x != localBegin.x || q.y != localBegin.y)
                            && (parameter == 1.0 || q.x != localEnd.x || q.y != localEnd.y),
                            "interior physical-segment root collapsed to a cell-local endpoint");
                        Detail::require(rootIndex == 0 || q.x != previousPoint.x || q.y != previousPoint.y,
                            "distinct physical-segment roots collapse to one cell-local point");
                        partRoots.push_back({cellI, cellJ, begin, end, localBegin, localEnd, q, parameter, false});
                        previousPoint = q;
                    }
                }
            // Equal representatives deliberately remain separate records.
            std::stable_sort(partRoots.begin(), partRoots.end(),
                [](const SegmentRootIncidence& x, const SegmentRootIncidence& y)
                { return x.parameter < y.parameter; });
            rootIncidences->insert(rootIncidences->end(), partRoots.begin(), partRoots.end());
            continue;
        }
        namespace Exact = ExactExpansion;
        auto chartBegin = begin, chartEnd = end;
        auto coefficients = scalar_.segmentPolynomialExact(i, j, local(chartBegin), local(chartEnd));
        bool reverseChart = false;
        const std::size_t parametricStart = parametric ? parametric->size() : 0;
        if (parametric)
        {
            // A real phase interval next to t=1 may be smaller than one ULP
            // there, while its distance from that endpoint is representable.
            // Retain a chart originating at the endpoint with the smaller
            // exact |phi|. Recompute the SAME restriction from its reversed
            // anchors; never form 1-root or round the root onto an endpoint.
            const auto absolute = [](const Exact::Expansion& value) {
                return Exact::sign(value) < 0 ? Exact::subtract(Exact::Expansion{}, value) : value;
            };
            reverseChart = Exact::sign(Exact::subtract(absolute(Exact::evaluate(coefficients, 1.0)),
                absolute(Exact::evaluate(coefficients, 0.0)))) < 0;
            if (reverseChart)
            {
                std::swap(chartBegin, chartEnd);
                coefficients = scalar_.segmentPolynomialExact(i, j, local(chartBegin), local(chartEnd));
            }
        }
        PolynomialSegmentRoots::Result isolated;
        try
        {
            isolated = parametric ? PolynomialSegmentRoots::isolateEndpointRelative(coefficients)
                : PolynomialSegmentRoots::isolate(coefficients);
        }
        catch (const std::exception& error)
        {
            const auto localBegin = local(chartBegin), localEnd = local(chartEnd);
            std::ostringstream message;
            message << std::setprecision(17) << "[LevelSetHermiteGeometryCandidate] segment root isolation"
                << "; cell=(" << i << ',' << j << "), begin=(" << chartBegin.x << ',' << chartBegin.y
                << "), end=(" << chartEnd.x << ',' << chartEnd.y << "), localBegin=(" << localBegin.x << ',' << localBegin.y
                << "), localEnd=(" << localEnd.x << ',' << localEnd.y << "); coefficients=";
            for (const auto& coefficient : coefficients)
            {
                message << '[';
                for (double component : coefficient.components()) message << component << ',';
                message << ']';
            }
            message << ": " << error.what();
            throw std::runtime_error(message.str());
        }
        const int firstSign = Exact::sign(Exact::evaluate(coefficients, 0.0));
        const int lastSign = Exact::sign(Exact::evaluate(coefficients, 1.0));
        if (firstSign != Exact::sign(exactValue(chartBegin))
            || lastSign != Exact::sign(exactValue(chartEnd)))
        {
            const auto localBegin = local(chartBegin), localEnd = local(chartEnd);
            std::ostringstream message;
            message << std::setprecision(17)
                << "[LevelSetHermiteGeometryCandidate] authoritative restriction and point query disagree at a shared endpoint"
                << "; cell=(" << i << ',' << j << "), begin=(" << chartBegin.x << ',' << chartBegin.y
                << "), end=(" << chartEnd.x << ',' << chartEnd.y << "), localBegin=(" << localBegin.x << ',' << localBegin.y
                << "), localEnd=(" << localEnd.x << ',' << localEnd.y << ")"
                << "; restriction endpoint signs=(" << firstSign << ',' << lastSign
                << "), canonical endpoint signs=(" << Exact::sign(exactValue(chartBegin))
                << ',' << Exact::sign(exactValue(chartEnd)) << ')';
            throw std::runtime_error(message.str());
        }
        struct Boundary { double point, lo, hi; };
        std::vector<Boundary> roots{{0.0, 0.0, 0.0}};
        Detail::require(isolated.roots.size() == isolated.enclosures.size(),
            "root representatives and enclosures disagree in count");
        for (std::size_t rootIndex = 0; rootIndex < isolated.roots.size(); ++rootIndex)
        {
            const double root = isolated.roots[rootIndex];
            const auto bracket = isolated.enclosures[rootIndex];
            Detail::require(std::isfinite(root) && root >= 0.0 && root <= 1.0,
                "root isolator returned an invalid segment parameter");
            if (parametric)
                Detail::require((root != 0.0 || firstSign == 0) && (root != 1.0 || lastSign == 0),
                    "interior root collapsed to a nonzero chart endpoint; phase interval cannot be discarded");
            if (root > 0.0 && root < 1.0)
            {
                Detail::require(root > roots.back().point && bracket.lo >= roots.back().hi
                    && bracket.lo <= root && root <= bracket.hi,
                    "root isolator returned unsorted or inconsistent enclosures");
                roots.push_back({root, bracket.lo, bracket.hi});
            }
        }
        roots.push_back({1.0, 1.0, 1.0});
        if (parametric && firstSign != 0 && lastSign != 0)
            Detail::require(lastSign == ((roots.size()-2)%2 == 0 ? firstSign : -firstSign),
                "simple-root inventory disagrees with exact endpoint phase signs");
        for (std::size_t r = 1; r < roots.size(); ++r)
        {
            const auto& left = roots[r - 1];
            const auto& right = roots[r];
            Detail::require(left.hi <= right.lo, "overlapping Hermite root enclosures leave no proven phase interval");
            const auto first = parametric ? Vector2D{} : Detail::between(begin, end, left.point);
            const auto last = parametric ? Vector2D{} : Detail::between(begin, end, right.point);
            // Root representatives describe geometry to their enclosure
            // accuracy. Their midpoint need not lie outside both enclosures!
            // Establish phase in the certified root-free gap instead. For a
            // near-domain-end root whose enclosure reaches the endpoint, the
            // exact nonzero endpoint sign supplies the one-sided phase.
            int phaseSign = 0;
            if (parametric && firstSign != 0)
            {
                // The isolator admits only SIMPLE interior roots: each has
                // a strict nonzero derivative enclosure; unresolved critical
                // zeros reject. Thus phase flips at every interior root.
                // Propagate the exact endpoint sign instead of evaluating
                // high powers of a sub-ULP endpoint parameter, whose exact
                // product tails can underflow even when the sign is proven.
                // Endpoint zeros themselves are not counted as crossings.
                phaseSign = (r-1)%2 == 0 ? firstSign : -firstSign;
            }
            else if (parametric && lastSign != 0)
                phaseSign = (roots.size()-1-r)%2 == 0 ? lastSign : -lastSign;
            else if (r == 1 && right.lo == 0.0) phaseSign = firstSign;
            else if (r + 1 == roots.size() && left.hi == 1.0) phaseSign = lastSign;
            else
            {
                const double middle = 0.5 * left.hi + 0.5 * right.lo;
                for (double t : {middle, left.hi, right.lo})
                {
                    phaseSign = Exact::sign(Exact::evaluate(coefficients, t));
                    if (phaseSign != 0) break;
                }
            }
            Detail::require(isolated.identicallyZero || phaseSign != 0,
                "unresolved scalar zero inside a purported root-free phase interval");
            const bool isNegative = phaseSign <= 0;
            // A grid cut is not an interface. Merge adjacent SAME-phase
            // parameter intervals before asking whether their physical image
            // has an interior representable double point. Opposite phases
            // are never removed or merged, including sub-ULP intervals.
            if (parametric)
            {
                Detail::require(left.point >= 0.0 && left.point < right.point && right.point <= 1.0,
                    "unrepresentable or unordered parametric phase interval");
                parametric->push_back({chartBegin, chartEnd, left.point, right.point, isNegative});
            }
            else appendPiece(first, last, isNegative);
        }
        if (parametric && reverseChart)
        {
            // Preserve traversal order without converting the retained chart
            // coordinate back to 1-t, which would erase its endpoint phase.
            std::reverse(parametric->begin() + parametricStart, parametric->end());
            for (std::size_t r = parametricStart; r < parametric->size(); ++r)
                std::swap((*parametric)[r].tBegin, (*parametric)[r].tEnd);
        }
    }
    if (rootIncidences) return pieces; // An empty inventory is a valid no-zero result.
    if (parametric)
    {
        Detail::require(!parametric->empty(), "empty parametric Hermite segment partition");
        return pieces; // No global root round-trip or physical-interior claim.
    }
    Detail::require(!pieces.empty(), "empty Hermite segment partition");
    for (const auto& piece : pieces)
    {
        Detail::require(piece.begin.x != piece.end.x || piece.begin.y != piece.end.y,
            "distinct Hermite crossings define an unrepresentable physical phase interval");
        const auto probe = Detail::between(piece.begin, piece.end, 0.5);
        if (!((probe.x != piece.begin.x || probe.y != piece.begin.y)
            && (probe.x != piece.end.x || probe.y != piece.end.y)))
        {
            std::ostringstream message;
            message << std::setprecision(17)
                << "[LevelSetHermiteGeometryCandidate] Hermite phase interval has no representable physical interior sample"
                << "; segment=(" << a.x << ',' << a.y << ")->(" << b.x << ',' << b.y
                << "), piece=(" << piece.begin.x << ',' << piece.begin.y << ")->("
                << piece.end.x << ',' << piece.end.y << "), negative=" << piece.negative;
            throw std::runtime_error(message.str());
        }
        Detail::require(negative(probe) == piece.negative,
            "authoritative scalar and rounded physical phase interval disagree");
    }
    return pieces;
}

template<class ScalarSnapshot>
inline double BasicSnapshot<ScalarSnapshot>::negativeFraction(Vector2D a, Vector2D b) const
{
    const auto pieces = partitionSegment(a, b);
    const long double length = std::hypot(static_cast<long double>(b.x) - a.x, static_cast<long double>(b.y) - a.y);
    long double negativeLength = 0.0L;
    for (const auto& piece : pieces)
        if (piece.negative) negativeLength += std::hypot(static_cast<long double>(piece.end.x) - piece.begin.x,
            static_cast<long double>(piece.end.y) - piece.begin.y);
    const double fraction = Detail::checked(negativeLength / length, "unrepresentable negative segment fraction");
    Detail::require(fraction >= 0.0 && fraction <= 1.0 + 128.0 * std::numeric_limits<double>::epsilon(),
        "Hermite negative segment fraction is outside [0,1]");
    return std::min(fraction, 1.0); // Final length-sum roundoff only; no phase interval is removed.
}
}
} // TU-local engine
} // namespace LAMP::Geometry::Detail
