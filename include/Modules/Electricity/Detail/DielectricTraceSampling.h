#pragma once

// Shared implementation-only lookup/window/value gathering, not an electric
// solver or a reconstruction algorithm. The two legacy policies retain their
// arithmetic and validation contracts. GlobalFaceSymmetric is an explicit
// physical-face policy, not a silent replacement for either legacy window.
// RZAxisOneSided is a separate load-trace policy using only real positive-r
// cells near the symmetry axis; it supplies no mirrored electric values.
// FittedActive explicitly combines actual uniform-leaf indexing with the SAME
// retained Snapshot's activity and phase. Physical-face windows keep the 6x5
// extent; old-local or explicit boundary-chart Frames keep the incident
// cell's SAME 5x5 extent without conversion between those coordinates. Only the
// declared radial symmetry axis retains the original 0..4 radial shift.
// Other walls only truncate that fixed range and remove inactive samples;
// no solid extension, mirrored potential, larger window or lower-order fit
// is supplied. The existing P3 rank/conditioning checks own sufficiency.
// Coordinates are copied from actual active mesh cells.
// NativeBoundaryWindow::OneSided additionally translates a fixed native box
// into the rectangular physical domain before removing inactive cells. This
// is phase-independent boundary support, not adaptive stencil growth or ghost
// values. Clipped retains the original boundary truncation as a control.
#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include <array>
#include <Eigen/Core>
#include "Core/Geometry/LevelSetHermiteGeometry.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/UniformLeafGrid.h"
#include "Core/Mesh/FaceSegmentGraph.h"
#include "Core/Mesh/DyadicGridCoordinates.h"
#include "Modules/Electricity/Detail/DielectricInterfaceTrace.h"

// Implementation-only: the anonymous ancestor gives every definition and
// type TU-local linkage, so different floating-point modes cannot share these
// engines through COMDAT. All includes stay outside these namespaces; no
// Detail type is a cross-TU public API.
namespace Quantum::Electricity::Detail
{
namespace
{
namespace DielectricTraceSampling
{
enum class Policy { GlobalClipped, TraceStrict, GlobalFaceSymmetric, RZAxisOneSided, FittedActive, FittedNative };
enum class NativeBoundaryWindow { Clipped, OneSided };

struct Window
{
    std::vector<DielectricInterfaceTrace::SamplePoint> samples;
    // Active-cell indexes, NOT labels. Both vectors have identical row-major
    // order: y/row outside, x/column inside. Phase comes only from the caller.
    std::vector<int> indices;
};

namespace Detail
{
inline void require(bool condition, const char* message)
{
    if (!condition)
        throw std::runtime_error(std::string("[DielectricTraceSampling] ") + message);
}
inline double finite(long double value)
{
    const double result = static_cast<double>(value);
    require(std::isfinite(value) && std::isfinite(result), "non-finite or unrepresentable measurement");
    return result;
}
inline void validate(const Window& window)
{
    require(window.samples.size() == window.indices.size() && !window.indices.empty(),
        "sample/index window sizes are inconsistent or empty");
}
}

// Immutable index/coordinate snapshot. Mesh and geometry must not change
// between construction, window construction and use of the sampled potentials.
// Callbacks must be deterministic, read-only actual-geometry/value queries.
class Lookup
{
public:
    Lookup(const QuadtreeMesh& mesh, Policy policy) : policy_(policy), nx_(mesh.Nx), ny_(mesh.Ny),
        min_(mesh.minCoord), hx_(mesh.cellSize.x), hy_(mesh.cellSize.y)
    {
        Detail::require(policy == Policy::GlobalClipped || policy == Policy::TraceStrict
            || policy == Policy::GlobalFaceSymmetric || policy == Policy::RZAxisOneSided,
            "unknown sampling policy");
        if (policy == Policy::GlobalClipped) initializeGlobal(mesh);
        else initializeStrict(mesh);
        if (policy == Policy::RZAxisOneSided)
            Detail::require(min_.x == 0.0, "RZ one-sided sampling requires the actual radial axis at r=0");
    }

    // Explicit fitted-active mode. Keep the referent immutable through every
    // alias; this lookup owns its actual cell/face coordinate maps and source.
    // The old policy constructor and all four old policies are unchanged.
    Lookup(const QuadtreeMesh& mesh, std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> source,
        CartesianGridScope scope)
        : policy_(Policy::FittedActive), nx_(0), ny_(0), min_(0.0,0.0), hx_(0.0), hy_(0.0),
          source_(std::move(source))
    {
        Detail::require(static_cast<bool>(source_), "fitted-active sampling requires a retained source");
        initializeFitted(mesh, scope);
    }

    // Explicit native-leaf policy. Use actual active centers in a fixed
    // Cartesian box around the anchor CELL center, +/-2 of that leaf's widths.
    // Membership uses retained dyadic addresses, not subtraction of rounded
    // physical coordinates at the box boundary.
    // Only the axis shifts the radial box to the original one-sided extent.
    // Sample count can vary across levels; do not fabricate uniform values,
    // grow a deficient stencil, lower P3 order, or consume solid payload.
    Lookup(const QuadtreeMesh& mesh, std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> source,
        NativeBoundaryWindow boundaryWindow=NativeBoundaryWindow::Clipped)
        : policy_(Policy::FittedNative), nx_(0), ny_(0), min_(0,0), hx_(0), hy_(0), source_(std::move(source)),
          nativeBoundaryWindow_(boundaryWindow)
    {
        Detail::require(boundaryWindow==NativeBoundaryWindow::Clipped || boundaryWindow==NativeBoundaryWindow::OneSided,
            "unknown native electric boundary-window policy");
        initializeNative(mesh);
    }

    Window window(const Vector2D& anchor,
        const std::function<bool(const Vector2D&)>& negative) const
    {
        using Detail::require;
        require(policy_ != Policy::FittedActive && policy_ != Policy::FittedNative,
            "fitted-active trace sampling requires the actual retained-source local Frame");
        require(policy_ != Policy::GlobalFaceSymmetric,
            "symmetric physical-face policy requires faceWindow and its orientation");
        require(static_cast<bool>(negative), "actual-phase query is empty");
        int ix, iy, firstRow, lastRow, firstColumn, lastColumn;
        if (policy_ == Policy::GlobalClipped)
        {
            // Preserve the original faceCenter DOUBLE arithmetic, truncation
            // toward zero, then clamp. This is NOT rounding or strict xGamma.
            // As in the old global fixture, this policy assumes finite valid
            // mesh/face coordinates supplied by the checked fixture/graph.
            ix = std::clamp(static_cast<int>((anchor.x-min_.x)/hx_),0,nx_-1);
            iy = std::clamp(static_cast<int>((anchor.y-min_.y)/hy_),0,ny_-1);
            firstRow=std::max(0,iy-2); lastRow=std::min(ny_-1,iy+2);
            firstColumn=std::max(0,ix-2); lastColumn=std::min(nx_-1,ix+2);
        }
        else if (policy_ == Policy::RZAxisOneSided)
        {
            // Explicit positive-radius data policy, not clipped legacy data
            // or parity-generated electric samples. Only the radial axis may
            // shift a window; axial and outer-radial boundaries stay strict.
            const double x=Detail::finite((static_cast<long double>(anchor.x)-min_.x)/hx_);
            const double y=Detail::finite((static_cast<long double>(anchor.y)-min_.y)/hy_);
            require(x>0. && x<nx_-2 && y>=2. && y<ny_-2,
                "RZ 5x5 trace requires positive anchor radius and separation from non-axis walls");
            ix=static_cast<int>(x); iy=static_cast<int>(y);
            firstRow=iy-2; lastRow=iy+2;
            firstColumn=ix<2 ? 0 : ix-2; lastColumn=firstColumn+4;
        }
        else
        {
            // Preserve LONG DOUBLE normalization, conversion through finite
            // double, range gate, and only then the original integer cast.
            const double x=Detail::finite((static_cast<long double>(anchor.x)-min_.x)/hx_);
            const double y=Detail::finite((static_cast<long double>(anchor.y)-min_.y)/hy_);
            require(x>=2. && x<nx_-2 && y>=2. && y<ny_-2,
                "fixed centered 5x5 trace window crosses domain boundary; no clipping or adaptive replacement is allowed");
            ix=static_cast<int>(x); iy=static_cast<int>(y);
            firstRow=iy-2; lastRow=iy+2; firstColumn=ix-2; lastColumn=ix+2;
        }
        return makeWindow(firstRow, lastRow, firstColumn, lastColumn, negative);
    }

    // Anchor the new fixed local/chart window by its actual incident cell, never by
    // inverting a rounded display point. Recheck the geometric Frame against
    // this source; projection diagnostics such as stationarity/iterations are
    // query-dependent and are not scalar-source identity fields.
    Window window(const LevelSetHermiteGeometry::Frame& anchor) const
    {
        using Detail::require;
        if (policy_ == Policy::FittedNative) return nativeWindow(anchor);
        require(policy_ == Policy::FittedActive && source_,
            "local-Frame window requires fitted-active retained-source sampling");
        require(anchor.hasCellLocalPoint != anchor.hasCellChartPoint,
            "fitted-active trace requires exactly one local or boundary-chart provenance");
        require(anchor.cellI >= 0 && anchor.cellI < nx_
            && anchor.cellJ >= 0 && anchor.cellJ < ny_
            && active(anchor.cellI, anchor.cellJ),
            "fitted-active trace anchor must retain an actual active cell");
        // Both public queries validate their own closed-cell representation.
        // In particular, never form 1-tail or invert a rounded display point
        // to turn the explicit boundary chart into the old local coordinate.
        const auto actual = anchor.hasCellChartPoint
            ? source_->jetAtCellChart(anchor.cellI, anchor.cellJ, anchor.cellChartPoint)
            : source_->jetAtCellLocal(anchor.cellI, anchor.cellJ, anchor.cellLocalPoint);
        const bool sameCoordinate = anchor.hasCellChartPoint
            ? actual.cellChartPoint.x.offset == anchor.cellChartPoint.x.offset
                && actual.cellChartPoint.x.fromUpper == anchor.cellChartPoint.x.fromUpper
                && actual.cellChartPoint.y.offset == anchor.cellChartPoint.y.offset
                && actual.cellChartPoint.y.fromUpper == anchor.cellChartPoint.y.fromUpper
            : actual.cellLocalPoint.x == anchor.cellLocalPoint.x && actual.cellLocalPoint.y == anchor.cellLocalPoint.y;
        require(actual.hasCellLocalPoint == anchor.hasCellLocalPoint && actual.hasCellChartPoint == anchor.hasCellChartPoint
            && sameCoordinate && actual.cellI == anchor.cellI && actual.cellJ == anchor.cellJ
            && actual.point.x == anchor.point.x && actual.point.y == anchor.point.y
            && actual.normal.x == anchor.normal.x && actual.normal.y == anchor.normal.y
            && actual.curvature == anchor.curvature && actual.levelSetResidual == anchor.levelSetResidual,
            "fitted-active trace Frame does not match the retained source's actual local/chart jet");
        const int ix = anchor.cellI, iy = anchor.cellJ;
        int firstColumn = ix - 2, lastColumn = ix + 2;
        if (source_->radialAxisSymmetry() && ix < 2)
        {
            require(min_.x == 0.0, "radial one-sided fitted sampling requires the actual axis at r=0");
            firstColumn = 0; lastColumn = 4; // SAME explicit RZAxisOneSided radial shift.
        }
        return makeFittedWindow(std::max(0, iy - 2), std::min(ny_ - 1, iy + 2),
            std::max(0, firstColumn), std::min(nx_ - 1, lastColumn));
    }

    // An interior physical face lies BETWEEN two cell centers in its normal
    // direction. A five-cell window selected by floor therefore has a 2/3
    // directional bias. Include all six centers at normal distances <=2.5h
    // and five at tangential distances <=2h instead. Physical-box clipping is
    // reflection covariant. This does not change the local x_Gamma policy.
    Window faceWindow(const Vector2D& anchor, bool horizontal,
        const std::function<bool(const Vector2D&)>& negative) const
    {
        using Detail::require;
        require(policy_ == Policy::GlobalFaceSymmetric,
            "faceWindow requires the explicit symmetric physical-face policy");
        require(static_cast<bool>(negative), "actual-phase query is empty");
        const double x = Detail::finite((static_cast<long double>(anchor.x) - min_.x) / hx_);
        const double y = Detail::finite((static_cast<long double>(anchor.y) - min_.y) / hy_);
        const double column = std::round(x - (horizontal ? 0.5 : 0.0));
        const double row = std::round(y - (horizontal ? 0.0 : 0.5));
        require(horizontal ? (column >= 0. && column < nx_ && row >= 1. && row < ny_)
            : (column >= 1. && column < nx_ && row >= 0. && row < ny_),
            "symmetric sampling anchor must be an interior physical face");
        const double tolerance = 256. * std::numeric_limits<double>::epsilon()
            * std::max({1., static_cast<double>(nx_), static_cast<double>(ny_)});
        require(std::abs(x - (column + (horizontal ? 0.5 : 0.0))) <= tolerance
            && std::abs(y - (row + (horizontal ? 0.0 : 0.5))) <= tolerance,
            "symmetric sampling anchor/orientation is not on the physical-face lattice");
        const int i = static_cast<int>(column), j = static_cast<int>(row);
        const int firstRow = std::max(0, j - (horizontal ? 3 : 2));
        const int lastRow = std::min(ny_ - 1, j + 2);
        const int firstColumn = std::max(0, i - (horizontal ? 2 : 3));
        const int lastColumn = std::min(nx_ - 1, i + 2);
        return makeWindow(firstRow, lastRow, firstColumn, lastColumn, negative);
    }

    // A full ACTUAL native face center, including outer and active/solid
    // boundaries. This does not accept a quarter-cell half-link midpoint or
    // reconstruct a face identity from a merely nearby coordinate.
    Window faceWindow(const Vector2D& anchor, bool horizontal) const
    {
        using Detail::require;
        require(policy_ == Policy::FittedActive && source_,
            "fitted faceWindow requires fitted-active retained-source sampling");
        const double x = Detail::finite((static_cast<long double>(anchor.x) - min_.x) / hx_);
        const double y = Detail::finite((static_cast<long double>(anchor.y) - min_.y) / hy_);
        const double column = std::round(x - (horizontal ? 0.5 : 0.0));
        const double row = std::round(y - (horizontal ? 0.0 : 0.5));
        require(horizontal ? (column >= 0.0 && column < nx_ && row >= 0.0 && row <= ny_)
            : (column >= 0.0 && column <= nx_ && row >= 0.0 && row < ny_),
            "fitted sampling anchor must be on an actual complete physical face");
        const int i = static_cast<int>(column), j = static_cast<int>(row);
        const auto& faces = horizontal ? horizontalFaces_ : verticalFaces_;
        const auto slot = static_cast<std::size_t>(j) * (horizontal ? nx_ : nx_ + 1) + i;
        require(slot < faces.size() && faces[slot].x == anchor.x && faces[slot].y == anchor.y,
            "fitted sampling anchor/orientation is not the actual native face center");
        require(horizontal ? (active(i, j - 1) || active(i, j)) : (active(i - 1, j) || active(i, j)),
            "fitted physical face has no active incident sample cell");
        return makeFittedWindow(std::max(0, j - (horizontal ? 3 : 2)), std::min(ny_ - 1, j + 2),
            std::max(0, i - (horizontal ? 2 : 3)), std::min(nx_ - 1, i + 2));
    }

    // Actual split-face identity, not a uniform lattice reconstructed from a
    // nearby point. The fixed 6-by-5 physical footprint uses the finer incident
    // leaf; clipping removes centers and never widens a deficient P3 stencil.
    Window faceWindow(const QuadtreeMesh& mesh, const FaceSegment& face) const
    {
        using Detail::require;
        require(policy_ == Policy::FittedNative && source_ && face.K >= 0 && face.L >= 0 &&
                    face.K < mesh.cellLabelCapacity() && face.L < mesh.cellLabelCapacity(),
                "native face sampling requires two actual incident leaves");
        const auto &edge = face.isHorizontal ? mesh.hEdges(face.edgeLabel) : mesh.vEdges(face.edgeLabel);
        require(edge.label() == face.edgeLabel && edge.index() >= 0 &&
                    ((edge.masterCell == face.K && edge.slaveCell == face.L) ||
                     (edge.masterCell == face.L && edge.slaveCell == face.K)) &&
                    edge.midpointCoord().x == face.faceCenter.x &&
                    edge.midpointCoord().y == face.faceCenter.y && edge.length() == face.Af,
                "native electric window has stale physical face identity");
        const auto &a = mesh.cells(face.K);
        const auto &b = mesh.cells(face.L);
        require(a.index() >= 0 && b.index() >= 0 && source_->nativeCellActive(a.index()) &&
                    source_->nativeCellActive(b.index()),
                "native electric face contains an inactive sample cell");
        const int fine = a.level() >= b.level() ? a.index() : b.index();
        const auto bounds = nativeBounds_.at(fine);
        const auto address = nativeAddresses_.at(fine);
        const bool upper = face.isHorizontal ? face.faceCenter.y == bounds[1].y
                                             : face.faceCenter.x == bounds[1].x;
        require(upper || (face.isHorizontal ? face.faceCenter.y == bounds[0].y
                                           : face.faceCenter.x == bounds[0].x),
                "native electric face is not on its fine incident leaf boundary");
        const auto i = 2 * address.i + (face.isHorizontal ? 1 : upper ? 2 : 0);
        const auto j = 2 * address.j + (face.isHorizontal ? (upper ? 2 : 0) : 1);
        const auto center = nativeCoordinates_->point(i, j, address.level + 1);
        require(center.x == face.faceCenter.x && center.y == face.faceCenter.y,
                "native electric split face does not match its dyadic fine-leaf address");
        return addressWindow(i, j, address.level + 1, face.isHorizontal ? 4 : 5,
                             face.isHorizontal ? 5 : 4);
    }

  private:
    Policy policy_;
    int nx_, ny_;
    Vector2D min_;
    double hx_, hy_;
    std::vector<int> grid_;
    std::vector<Vector2D> coordinates_;
    std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> source_;
    std::vector<Vector2D> horizontalFaces_, verticalFaces_;
    std::vector<std::array<Vector2D, 2>> nativeBounds_;
    std::vector<Quantum::Mesh::DyadicCellAddress> nativeAddresses_;
    std::optional<Quantum::Mesh::DyadicGridCoordinates> nativeCoordinates_;
    NativeBoundaryWindow nativeBoundaryWindow_=NativeBoundaryWindow::Clipped;

    void initializeNative(const QuadtreeMesh &mesh)
    {
        using Detail::require;
        require(source_ && source_->usesNativeLeaves() && mesh.cellIndexToLabel &&
                    source_->nativeCellCount() == static_cast<std::size_t>(mesh.cellIndexNum()),
                "native electric sampling requires the actual complete leaf layout");
        require(source_->minCoord().x == mesh.minCoord.x && source_->minCoord().y == mesh.minCoord.y &&
                    source_->maxCoord().x == mesh.maxCoord.x && source_->maxCoord().y == mesh.maxCoord.y,
                "native electric source and mesh bounds disagree");
        min_ = mesh.minCoord;
        nx_ = mesh.Nx; ny_ = mesh.Ny;
        nativeCoordinates_.emplace(mesh.minCoord,mesh.maxCoord,mesh.Nx,mesh.Ny);
        coordinates_.resize(mesh.cellIndexNum());
        nativeBounds_.resize(mesh.cellIndexNum());
        nativeAddresses_.resize(mesh.cellIndexNum());
        for (int index = 0; index < mesh.cellIndexNum(); ++index)
        {
            const int label = mesh.cellIndexToLabel[index];
            require(label >= 0 && label < mesh.cellLabelCapacity(), "native electric sample label out of range");
            const auto& cell = mesh.cells(label);
            require(cell.index() == index && cell.label() == label &&
                std::all_of(cell.son, cell.son+4, [](int son) { return son < 0; }),
                "native electric sample is not its actual leaf");
            const auto bounds = source_->nativeCellBounds(index);
            Vector2D lo(std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity());
            Vector2D hi(-lo.x, -lo.y);
            for (int vertex : cell.cellToVert)
            {
                const auto point = mesh.verts(vertex).coord();
                lo.x = std::min(lo.x,point.x); lo.y = std::min(lo.y,point.y);
                hi.x = std::max(hi.x,point.x); hi.y = std::max(hi.y,point.y);
            }
            require(bounds[0].x == lo.x && bounds[0].y == lo.y && bounds[1].x == hi.x && bounds[1].y == hi.y,
                "native electric source has stale actual leaf bounds");
            nativeBounds_[index] = bounds;
            coordinates_[index] = cell.centerCoord();
            const auto address = Quantum::Mesh::cellAddress(mesh,cell);
            nativeAddresses_[index] = address;
            const auto center = nativeCoordinates_->point(2*address.i+1,2*address.j+1,address.level+1);
            const auto lower = nativeCoordinates_->point(address.i,address.j,address.level);
            const auto upper = nativeCoordinates_->point(address.i+1,address.j+1,address.level);
            require(center.x == coordinates_[index].x && center.y == coordinates_[index].y
                && lower.x == lo.x && lower.y == lo.y && upper.x == hi.x && upper.y == hi.y,
                "native electric sample has inconsistent dyadic coordinates");
        }
    }

    Window nativeWindow(const LevelSetHermiteGeometry::Frame& anchor) const
    {
        using Detail::require;
        const int cell = anchor.nativeCellIndex;
        require(source_ && anchor.hasCellChartPoint && !anchor.hasCellLocalPoint &&
            anchor.cellI == -1 && anchor.cellJ == -1 && cell >= 0 &&
            cell < static_cast<int>(nativeBounds_.size()) && source_->nativeCellActive(cell),
            "native electric trace requires an actual active leaf/chart frame");
        const auto actual = source_->jetAtNativeCellChart(cell, anchor.cellChartPoint);
        require(actual.point.x == anchor.point.x && actual.point.y == anchor.point.y &&
            actual.normal.x == anchor.normal.x && actual.normal.y == anchor.normal.y &&
            actual.curvature == anchor.curvature && actual.levelSetResidual == anchor.levelSetResidual,
            "native electric frame does not match its retained scalar jet");
        const auto address = nativeAddresses_[cell];
        auto i = 2*address.i+1;
        if (source_->radialAxisSymmetry() && i < 5)
        {
            require(min_.x == 0, "native one-sided electric window requires the actual r=0 axis");
            i = 5;
        }
        return addressWindow(i,2*address.j+1,address.level+1,4,4);
    }

    Window addressWindow(std::uint64_t i,std::uint64_t j,unsigned level,
                         std::int64_t halfI,std::int64_t halfJ) const
    {
        // Canonical physical endpoints retrieve a superset. Exact dyadic
        // center addresses decide membership, including a center ON an edge
        // of the fixed stencil. Rounded physical widths can otherwise remove
        // an entire outer row on a non-binary domain, losing the only liquid
        // sample before a resolved-to-subcell neck transition.
        const auto bound = [](std::uint64_t value,std::int64_t offset,std::uint64_t count) {
            return static_cast<std::uint64_t>(std::clamp(static_cast<std::int64_t>(value)+offset,
                std::int64_t{0},static_cast<std::int64_t>(count)));
        };
        Detail::require(level<=53 && static_cast<std::uint64_t>(nx_)<=((std::uint64_t{1}<<53)>>level)
            && static_cast<std::uint64_t>(ny_)<=((std::uint64_t{1}<<53)>>level),
            "native electric window exceeds exact dyadic address range");
        const auto countI=static_cast<std::uint64_t>(nx_)<<level,countJ=static_cast<std::uint64_t>(ny_)<<level;
        if(nativeBoundaryWindow_==NativeBoundaryWindow::OneSided) {
            // At this half-leaf address scale actual same-level centers are
            // odd integers. Shift by whole cells, preserving the nominal
            // five/six-center footprint and its parity. This does not depend
            // on phase, rank or a previous failed reconstruction.
            const auto interiorCenter=[](std::uint64_t center,std::uint64_t count,std::int64_t half) {
                const auto margin=static_cast<std::uint64_t>(half+1);
                return count>=2*margin?std::clamp(center,margin,count-margin):center;
            };
            i=interiorCenter(i,countI,halfI);j=interiorCenter(j,countJ,halfJ);
        }
        const auto lower=nativeCoordinates_->point(bound(i,-halfI,countI),bound(j,-halfJ,countJ),level);
        const auto upper=nativeCoordinates_->point(bound(i,halfI,countI),bound(j,halfJ,countJ),level);
        const auto candidates = source_->nativeCellCandidates(lower,upper);
        const auto logical=[&](std::int64_t value){return std::ldexp(static_cast<double>(value),-static_cast<int>(level));};
        const double lowI=logical(static_cast<std::int64_t>(i)-halfI),highI=logical(static_cast<std::int64_t>(i)+halfI);
        const double lowJ=logical(static_cast<std::int64_t>(j)-halfJ),highJ=logical(static_cast<std::int64_t>(j)+halfJ);
        Window result;
        for (int index : candidates)
        {
            if (!source_->nativeCellActive(index))
                continue;
            const auto a=nativeAddresses_[index];
            const double x=std::ldexp(static_cast<double>(2*a.i+1),-static_cast<int>(a.level+1));
            const double y=std::ldexp(static_cast<double>(2*a.j+1),-static_cast<int>(a.level+1));
            if (x>=lowI && x<=highI && y>=lowJ && y<=highJ)
                result.indices.push_back(index);
        }
        std::sort(result.indices.begin(), result.indices.end(), [&](int a, int b) {
            const auto p = coordinates_[a], q = coordinates_[b];
            return p.y != q.y ? p.y < q.y : p.x != q.x ? p.x < q.x : a < b;
        });
        for (int index : result.indices)
        {
            const auto point = coordinates_[index];
            result.samples.push_back({point, source_->negative(point)
                                                 ? DielectricInterfaceTrace::Phase::Negative
                                                 : DielectricInterfaceTrace::Phase::Positive});
        }
        Detail::validate(result);
        return result;
    }

    bool active(int i, int j) const
    {
        return i >= 0 && i < nx_ && j >= 0 && j < ny_ && grid_[static_cast<std::size_t>(j) * nx_ + i] >= 0;
    }

    Window makeFittedWindow(int firstRow, int lastRow, int firstColumn, int lastColumn) const
    {
        Window result;
        for (int row = firstRow; row <= lastRow; ++row)
            for (int column = firstColumn; column <= lastColumn; ++column)
            {
                const int index = grid_[static_cast<std::size_t>(row) * nx_ + column];
                if (index < 0)
                    continue; // Inactive payload is never queried or gathered.
                const auto sample = coordinates_[index];
                result.samples.push_back({sample, source_->negative(sample)
                                                      ? DielectricInterfaceTrace::Phase::Negative
                                                      : DielectricInterfaceTrace::Phase::Positive});
                result.indices.push_back(index);
            }
        Detail::validate(result);
        return result;
    }

    void initializeFitted(const QuadtreeMesh& mesh, CartesianGridScope scope)
    {
        using Detail::require;
        const UniformLeafGrid layout(mesh, scope);
        nx_ = layout.nx(); ny_ = layout.ny(); min_ = layout.minCoord();
        const auto spacing = layout.spacing(), upper = layout.maxCoord();
        hx_ = spacing.x; hy_ = spacing.y;
        const auto sourceMin = source_->minCoord(), sourceMax = source_->maxCoord(), sourceSize = source_->cellSize();
        require(nx_ >= 5 && ny_ >= 5 && source_->nx() == nx_ && source_->ny() == ny_
            && sourceMin.x == min_.x && sourceMin.y == min_.y && sourceMax.x == upper.x && sourceMax.y == upper.y
            && sourceSize.x == hx_ && sourceSize.y == hy_,
            "fitted-active source and actual uniform-leaf lattice must agree and contain at least 5x5 cells");
        grid_.assign(layout.cells().size(), -1);
        coordinates_.resize(mesh.cellIndexNum());
        for (int j = 0; j < ny_; ++j) for (int i = 0; i < nx_; ++i)
        {
            const auto slot = static_cast<std::size_t>(j) * nx_ + i;
            const int label = layout.cells()[slot];
            const auto& cell = mesh.cells(label);
            const int index = cell.index();
            require(index >= 0 && index < mesh.cellIndexNum() && mesh.cellIndexToLabel[index] == label,
                "fitted sample has a stale native cell index/label");
            if (!source_->cellActive(i, j)) continue;
            grid_[slot] = index;
            coordinates_[index] = cell.centerCoord();
        }
        horizontalFaces_.reserve(layout.horizontal().size());
        for (int label : layout.horizontal()) horizontalFaces_.push_back(mesh.hEdges(label).midpointCoord());
        verticalFaces_.reserve(layout.vertical().size());
        for (int label : layout.vertical()) verticalFaces_.push_back(mesh.vEdges(label).midpointCoord());
    }

    Window makeWindow(int firstRow, int lastRow, int firstColumn, int lastColumn,
        const std::function<bool(const Vector2D&)>& negative) const
    {
        Window result;
        for(int row=firstRow;row<=lastRow;++row)
            for(int column=firstColumn;column<=lastColumn;++column)
            {
                const int k=grid_[static_cast<std::size_t>(row)*nx_+column];
                Detail::require(k>=0,"missing Cartesian sample cell");
                const auto sample=coordinates_[k];
                result.samples.push_back({sample,negative(sample)
                    ? DielectricInterfaceTrace::Phase::Negative
                    : DielectricInterfaceTrace::Phase::Positive});
                result.indices.push_back(k);
            }
        return result;
    }

    void initializeGlobal(const QuadtreeMesh& mesh)
    {
        // This is deliberately the old GLOBAL fixture's weaker lookup, not a
        // general mesh validator. Do not import the strict trace checks here:
        // validated fixed-mesh fixtures/FaceSegmentGraph own that precondition.
        grid_.assign(static_cast<std::size_t>(mesh.Nx)*mesh.Ny,-1);
        coordinates_.resize(mesh.cellIndexNum());
        const double hx=mesh.cellSize.x,hy=mesh.cellSize.y;
        for(int k=0;k<mesh.cellIndexNum();++k)
        {
            const auto x=mesh.cells(mesh.cellIndexToLabel[k]).centerCoord();
            const int i=static_cast<int>((x.x-mesh.minCoord.x)/hx),j=static_cast<int>((x.y-mesh.minCoord.y)/hy);
            Detail::require(i>=0 && i<mesh.Nx && j>=0 && j<mesh.Ny,"cell outside test Cartesian lookup");
            grid_[static_cast<std::size_t>(j)*mesh.Nx+i]=k;
            coordinates_[k]=x;
        }
    }

    void initializeStrict(const QuadtreeMesh& mesh)
    {
        using Detail::require;
        using Detail::finite;
        require(mesh.Nx>=5 && mesh.Ny>=5 && mesh.cellIndexToLabel,"complete uniform mesh with at least 5x5 cells required");
        const long long total=static_cast<long long>(mesh.Nx)*mesh.Ny;
        require(total==mesh.cellIndexNum(),"complete uniform-cell count contract failed");
        const double hx=mesh.cellSize.x,hy=mesh.cellSize.y,h=std::max(hx,hy);
        require(std::isfinite(hx) && hx>0. && std::isfinite(hy) && hy>0.,"uniform cell sizes must be finite and positive");
        require(std::isfinite(mesh.minCoord.x) && std::isfinite(mesh.minCoord.y)
            && std::isfinite(mesh.maxCoord.x) && std::isfinite(mesh.maxCoord.y)
            && mesh.maxCoord.x>mesh.minCoord.x && mesh.maxCoord.y>mesh.minCoord.y,"invalid Cartesian domain");
        const double coordinateScale=std::max({h,std::abs(mesh.minCoord.x),std::abs(mesh.minCoord.y),
            std::abs(mesh.maxCoord.x),std::abs(mesh.maxCoord.y)});
        const double tolerance=256.*std::numeric_limits<double>::epsilon()*coordinateScale;
        require(std::abs(finite(static_cast<long double>(hx)*mesh.Nx)-(mesh.maxCoord.x-mesh.minCoord.x))<=tolerance
            && std::abs(finite(static_cast<long double>(hy)*mesh.Ny)-(mesh.maxCoord.y-mesh.minCoord.y))<=tolerance,
            "domain extents disagree with uniform cell sizes");
        grid_.assign(static_cast<std::size_t>(total),-1);
        coordinates_.resize(mesh.cellIndexNum());
        for(int k=0;k<mesh.cellIndexNum();++k)
        {
            const int label=mesh.cellIndexToLabel[k];
            require(label>=0 && label<mesh.cellLabelCapacity(),"cell map label is outside allocated storage");
            const auto& cell=mesh.cells(label); const auto point=cell.centerCoord();
            require(cell.label()==label && cell.index()==k && cell.level()==0 && cell.parent==-1
                && std::all_of(cell.son,cell.son+4,[](int son){return son==-1;}),"uniform active leaf cell mapping required");
            require(std::isfinite(point.x) && std::isfinite(point.y)
                && std::isfinite(cell.size().x) && std::isfinite(cell.size().y)
                && std::abs(cell.size().x-hx)<=tolerance && std::abs(cell.size().y-hy)<=tolerance,
                "non-finite or non-uniform cell geometry");
            const double x=finite((static_cast<long double>(point.x)-mesh.minCoord.x)/hx);
            const double y=finite((static_cast<long double>(point.y)-mesh.minCoord.y)/hy);
            require(x>=0. && x<mesh.Nx && y>=0. && y<mesh.Ny,"cell center is outside Cartesian grid");
            const int ix=static_cast<int>(x),iy=static_cast<int>(y);
            require(std::abs(point.x-(mesh.minCoord.x+(ix+.5)*hx))<=tolerance
                && std::abs(point.y-(mesh.minCoord.y+(iy+.5)*hy))<=tolerance,"cell center is not on the uniform lattice");
            auto& slot=grid_[static_cast<std::size_t>(iy)*mesh.Nx+ix]; require(slot==-1,"duplicate Cartesian cell center"); slot=k;
            coordinates_[k]=point;
        }
        require(std::all_of(grid_.begin(),grid_.end(),[](int index){return index>=0;}),"missing Cartesian cell");
    }
};

inline std::vector<double> gather(const Window& window, const Eigen::VectorXd& values)
{
    Detail::validate(window);
    std::vector<double> result;
    for (const int k : window.indices)
    {
        Detail::require(k>=0 && k<values.size(),"potential index is outside the supplied active-cell array");
        Detail::require(std::isfinite(values[k]),"sampled potential is non-finite");
        result.push_back(values[k]);
    }
    return result;
}

inline std::vector<double> evaluate(const Window& window,
    const std::function<double(const Vector2D&)>& potential)
{
    Detail::validate(window);
    Detail::require(static_cast<bool>(potential),"potential query is empty");
    std::vector<double> result;
    for (const auto& sample : window.samples)
        result.push_back(Detail::finite(potential(sample.point)));
    return result;
}
}
} // TU-local engine
} // namespace Quantum::Electricity::Detail
