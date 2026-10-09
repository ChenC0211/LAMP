#include "Modules/Interface/InterfacePhaseGeometry.h"
#include "Core/Mesh/NativeLeafDomain.h"
#include "Core/Geometry/LevelSetHermiteGeometry.h"
#include "Core/Geometry/LevelSetMappedSegment.h"
#include "Core/Geometry/LevelSetInterfacePanels.h"
#include "Core/Geometry/AxisymmetricBoundaryIntegral.h"
#include "Core/Geometry/CartesianBoundaryIntegral.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/UniformLeafGrid.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <exception>
#include <limits>
#include <map>
#include <sstream>
#include <iomanip>
#include <stdexcept>
#include <string>
#include <tuple>
#include <type_traits>

namespace InterfacePhaseGeometry
{
namespace
{
void require(bool valid, const char* reason)
{
    if (!valid) throw std::runtime_error(std::string("[InterfacePhaseGeometry] ") + reason);
}
void validate(const QuadtreeMesh& mesh, const NativeLeafDomain& domain,
    const LevelSetHermiteGeometry::Snapshot& snapshot)
{
    domain.validateMeshIdentity(mesh);
    if (snapshot.usesNativeLeaves())
    {
        require(domain.axisymmetric() == snapshot.radialAxisSymmetry() &&
            snapshot.nativeCellCount() == static_cast<std::size_t>(mesh.cellIndexNum()) &&
            snapshot.minCoord().x == mesh.minCoord.x && snapshot.minCoord().y == mesh.minCoord.y &&
            snapshot.maxCoord().x == mesh.maxCoord.x && snapshot.maxCoord().y == mesh.maxCoord.y,
            "native scalar and selected XY/RZ leaf layout differ");
        for (int k = 0; k < mesh.cellIndexNum(); ++k)
        {
            const auto& cell = mesh.cells(mesh.cellIndexToLabel[k]);
            const auto low = mesh.verts(cell.cellToVert[0]).coord();
            const auto high = mesh.verts(cell.cellToVert[2]).coord();
            const auto bounds = snapshot.nativeCellBounds(k);
            require(bounds[0].x == low.x && bounds[0].y == low.y &&
                bounds[1].x == high.x && bounds[1].y == high.y &&
                snapshot.nativeCellActive(k) == (domain.cellRowByIndex()[k] >= 0),
                "native scalar leaf bounds/activity differ from the fitted fluid domain");
        }
        return;
    }
    const UniformLeafGrid grid(mesh, CartesianGridScope::UniformLeaves);
    require(domain.axisymmetric() && snapshot.radialAxisSymmetry()
        && snapshot.nx() == grid.nx() && snapshot.ny() == grid.ny()
        && snapshot.cellSize().x == grid.spacing().x && snapshot.cellSize().y == grid.spacing().y
        && snapshot.minCoord().x == grid.minCoord().x && snapshot.minCoord().y == grid.minCoord().y
        && snapshot.maxCoord().x == grid.maxCoord().x && snapshot.maxCoord().y == grid.maxCoord().y,
        "shared scalar and fitted RZ leaf layout differ");
    for (int j = 0; j < grid.ny(); ++j) for (int i = 0; i < grid.nx(); ++i)
    {
        const int k = mesh.cells(grid.cells()[static_cast<std::size_t>(j) * grid.nx() + i]).index();
        require(snapshot.cellActive(i, j) == (domain.cellRowByIndex()[k] >= 0),
            "shared scalar activity differs from the fitted fluid domain");
    }
}
}

void validateLayout(const QuadtreeMesh& mesh, const NativeLeafDomain& domain,
    const LevelSetHermiteGeometry::Snapshot& snapshot)
{ validate(mesh, domain, snapshot); }

struct MappedFacePartitions::Impl
{
    std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> source;
    const NativeLeafDomain* domain = nullptr;
    std::vector<std::array<std::unique_ptr<const LevelSetMappedSegment>,2>> halfFaces;
    std::vector<std::array<std::unique_ptr<const LevelSetMappedSegment::NativeQuarticTrace>,2>> traceHalfFaces;
    bool nativeTrace = false;
    std::size_t rootCount = 0;
    std::size_t eventCount = 0;
};
MappedFacePartitions::MappedFacePartitions(const QuadtreeMesh& mesh, const NativeLeafDomain& domain,
    std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> source)
    :MappedFacePartitions(mesh,domain,std::move(source),false){}
MappedFacePartitions MappedFacePartitions::prepareNativeTrace(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,
    std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> source)
{return MappedFacePartitions(mesh,domain,std::move(source),true);}
MappedFacePartitions::MappedFacePartitions(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,
    std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> source,bool nativeTrace)
{
    require(static_cast<bool>(source), "null mapped partition source");
    InterfacePhaseGeometry::validate(mesh,domain,*source);
    require(!nativeTrace || (source->usesNativeLeaves() && source->reconstruction()==LevelSetHermiteGeometry::Reconstruction::EnrichedQuartic),
        "native trace cache requires a native Q4 source");
    auto result = std::make_unique<Impl>();
    result->source = std::move(source); result->domain = &domain;result->nativeTrace=nativeTrace;
    const auto& faces = domain.faces(mesh);
    if(nativeTrace)result->traceHalfFaces.resize(faces.size());else result->halfFaces.resize(faces.size());
    const auto buildFace = [&](std::size_t index)
    {
        const auto& face = faces[static_cast<int>(index)];
        if (face.negativeCell < 0 && face.positiveCell < 0) return;
        const auto& edge = face.horizontal ? mesh.hEdges(face.edgeLabel) : mesh.vEdges(face.edgeLabel);
        auto a = mesh.verts(edge.edgeToVert[0]).coord(), b = mesh.verts(edge.edgeToVert[1]).coord();
        if (face.horizontal ? a.x > b.x : a.y > b.y) std::swap(a,b);
        for (int half = 0; half < 2; ++half)
        {
            const auto begin = half == 0 ? a : face.center, end = half == 0 ? face.center : b;
            try
            {
                if(nativeTrace)result->traceHalfFaces[face.nativeRow][half]=
                    std::make_unique<const LevelSetMappedSegment::NativeQuarticTrace>(
                        LevelSetMappedSegment::nativeQuarticTrace(result->source,begin,end,face.horizontal?1:0));
                else {
                    auto query = std::make_unique<const LevelSetMappedSegment>(result->source,begin,end);
                    result->halfFaces[face.nativeRow][half] = std::move(query);
                }
            }
            catch (const std::exception& error)
            {
                std::ostringstream message;
                message << std::setprecision(17) << "[InterfacePhaseGeometry] mapped face=" << face.nativeRow
                    << ", half=" << half << ", begin=(" << begin.x << ',' << begin.y << "), end=("
                    << end.x << ',' << end.y << "): " << error.what();
                throw std::runtime_error(message.str());
            }
        }
    };
    // Keep the first active query on the caller thread: it also enforces the
    // original strict floating-point environment contract before any workers.
    // Snapshot and mesh are immutable; each task owns a distinct face slot.
    std::size_t first = 0;
    while (first < faces.size() && faces[static_cast<int>(first)].negativeCell < 0
        && faces[static_cast<int>(first)].positiveCell < 0) ++first;
    if (first < faces.size()) buildFace(first);
    std::vector<std::exception_ptr> failures(faces.size());
    // Never throw through an OpenMP region. Preserve the original first-face
    // failure order after joining, and publish no partially constructed cache.
    // No floating-point sum is parallelized here or in the downstream ledger.
#pragma omp parallel for schedule(dynamic, 32) if(faces.size() > 128)
    for (std::ptrdiff_t index = static_cast<std::ptrdiff_t>(first)+1;
        index < static_cast<std::ptrdiff_t>(faces.size()); ++index)
    {
        try { buildFace(static_cast<std::size_t>(index)); }
        catch (...) { failures[index] = std::current_exception(); }
    }
    for (std::size_t index = first; index < faces.size(); ++index)
    {
        if (failures[index]) std::rethrow_exception(failures[index]);
        const auto row=faces[static_cast<int>(index)].nativeRow;
        if(nativeTrace){for(const auto& half:result->traceHalfFaces[row])if(half){
            result->eventCount+=half->events().size();
            for(const auto& e:half->events())if(e.normalOrder==0)++result->rootCount;}}
        else for (const auto& half : result->halfFaces[row])
            if (half) result->rootCount += half->roots().size();
    }
    impl_ = std::move(result);
}
MappedFacePartitions::~MappedFacePartitions() = default;
MappedFacePartitions::MappedFacePartitions(MappedFacePartitions&&) noexcept = default;
MappedFacePartitions& MappedFacePartitions::operator=(MappedFacePartitions&&) noexcept = default;
const MappedFacePartitions::Impl& MappedFacePartitions::checked() const
{
    require(static_cast<bool>(impl_), "query on a moved-from mapped partition cache");
    return *impl_;
}
void MappedFacePartitions::validate(const QuadtreeMesh& mesh, const NativeLeafDomain& domain) const
{
    const auto& data = checked();
    require(data.domain == &domain, "mapped partition cache belongs to a different fitted domain");
    InterfacePhaseGeometry::validate(mesh,domain,*data.source);
}
const LevelSetHermiteGeometry::Snapshot& MappedFacePartitions::source() const { return *checked().source; }
const LevelSetMappedSegment& MappedFacePartitions::halfFace(int nativeRow, int half) const
{
    const auto& data = checked();
    if(data.nativeTrace)return traceHalfFace(nativeRow,half).ordinary();
    require(nativeRow >= 0 && static_cast<std::size_t>(nativeRow) < data.halfFaces.size()
        && (half == 0 || half == 1), "invalid mapped native half-face index");
    const auto& query = data.halfFaces[nativeRow][half];
    require(static_cast<bool>(query), "solid-only face has no mapped phase partition");
    return *query;
}
const LevelSetMappedSegment::NativeQuarticTrace& MappedFacePartitions::traceHalfFace(int nativeRow,int half) const
{
    const auto& data=checked();require(data.nativeTrace,"typed half-face trace requires prepareNativeTrace");
    require(nativeRow>=0 && static_cast<std::size_t>(nativeRow)<data.traceHalfFaces.size() && (half==0 || half==1),
        "invalid mapped native half-face index");
    const auto& query=data.traceHalfFaces[nativeRow][half];
    require(static_cast<bool>(query),"solid-only face has no mapped phase partition");return *query;
}
std::size_t MappedFacePartitions::rootIncidenceCount() const { return checked().rootCount; }
std::size_t MappedFacePartitions::traceEventIncidenceCount() const
{const auto& data=checked();require(data.nativeTrace,"typed event count requires prepareNativeTrace");return data.eventCount;}
void MappedFacePartitions::refineRoot(int nativeRow, int half, std::size_t rootIndex)
{
    const auto& previous = halfFace(nativeRow, half);
    require(rootIndex < previous.roots().size(), "invalid mapped half-face root index");
    if(checked().nativeTrace){const auto& trace=traceHalfFace(nativeRow,half);const auto& event=trace.events().at(rootIndex);
        const auto& a=event.root;const auto& b=previous.roots()[rootIndex];
        require(trace.events().size()==previous.roots().size() && event.normalOrder==0
            && a.chart==b.chart && a.nativeCellIndex==b.nativeCellIndex && a.parameter==b.parameter
            && a.lower==b.lower && a.upper==b.upper && a.pointOnly==b.pointOnly,
            "ordinary root differs from its shared typed value event");
        refineTraceEvent(nativeRow,half,rootIndex);return;
    }
    try
    {
        auto replacement = std::make_unique<const LevelSetMappedSegment>(previous.refined(rootIndex));
        const auto same = [](Vector2D a, Vector2D b) { return a.x == b.x && a.y == b.y; };
        require(replacement->usesSource(*impl_->source)
            && same(replacement->begin(), previous.begin()) && same(replacement->end(), previous.end()),
            "root refinement changed the mapped source or actual half-face");
        require(replacement->roots().size() == previous.roots().size()
            && replacement->gridPieceCount() == previous.gridPieceCount(),
            "root refinement changed the mapped incidence or grid-piece count");
        // Keep the old query alive until every check succeeds. rootCount is
        // unchanged because the replacement has the same incidence count.
        impl_->halfFaces[nativeRow][half] = std::move(replacement);
    }
    catch (const std::exception& error)
    {
        std::ostringstream message;
        message << "[InterfacePhaseGeometry] refine mapped face=" << nativeRow
            << ", half=" << half << ", root=" << rootIndex << ": " << error.what();
        throw std::runtime_error(message.str());
    }
}
void MappedFacePartitions::refineTraceEvent(int nativeRow,int half,std::size_t eventIndex)
{
    const auto& previous=traceHalfFace(nativeRow,half);
    require(eventIndex<previous.events().size(),"invalid mapped half-face trace event index");
    try {
        auto replacement=std::make_unique<const LevelSetMappedSegment::NativeQuarticTrace>(previous.refined(eventIndex));
        const auto same=[](Vector2D a,Vector2D b){return a.x==b.x && a.y==b.y;};
        require(replacement->usesSource(*impl_->source) && same(replacement->begin(),previous.begin())
            && same(replacement->end(),previous.end()) && replacement->normalAxis()==previous.normalAxis(),
            "trace refinement changed the mapped source or actual half-face");
        require(replacement->events().size()==previous.events().size() && replacement->gridPieceCount()==previous.gridPieceCount(),
            "trace refinement changed the mapped incidence or grid-piece count");
        for(std::size_t k=0;k<previous.events().size();++k){const auto& a=previous.events()[k];const auto& b=replacement->events()[k];
            require(a.normalOrder==b.normalOrder && a.normalSide==b.normalSide && a.coincidentChart==b.coincidentChart
                && a.root.nativeCellIndex==b.root.nativeCellIndex && a.root.pointOnly==b.root.pointOnly,
                "trace refinement changed the mapped event identity");}
        impl_->traceHalfFaces[nativeRow][half]=std::move(replacement);
    } catch(const std::exception& error){std::ostringstream message;
        message<<"[InterfacePhaseGeometry] refine mapped trace face="<<nativeRow<<", half="<<half<<", event="<<eventIndex<<": "<<error.what();
        throw std::runtime_error(message.str());
    }
}

long double volume(const QuadtreeMesh& mesh, const NativeLeafDomain& domain,
    const LevelSetHermiteGeometry::Snapshot& snapshot, const LevelSetInterfacePanels& panels)
{
    validate(mesh, domain, snapshot);
    require(panels.axisymmetric(), "free panels must use physical RZ measures");
    require(!snapshot.usesNativeLeaves(), "native phase volume requires native free-curve assembly, not legacy panels");
    // Preserve the free-node then boundary-piece accumulation and arithmetic
    // order of the original App diagnostic. F=z e_z has unit divergence;
    // vertical and r=0 boundaries make no contribution to this volume formula.
    long double result = 0.0L;
    for (const auto& panel : panels.panels())
    {
        require(panel.cellIndex >= 0 && panel.cellIndex < mesh.cellIndexNum()
            && domain.cellRowByIndex()[panel.cellIndex] >= 0, "free panel is outside the active fluid domain");
        for (const auto& node : panel.quadratureNodes)
            result += static_cast<long double>(node.physicalWeight) * node.frame.point.y * node.frame.normal.y;
    }
    constexpr long double pi = 3.141592653589793238462643383279502884L;
    for (const auto& face : domain.faces(mesh))
    {
        if (!face.horizontal || face.outwardSign == 0.0) continue;
        const auto& edge = mesh.hEdges(face.edgeLabel);
        const auto a = mesh.verts(edge.edgeToVert[0]).coord(), b = mesh.verts(edge.edgeToVert[1]).coord();
        if (snapshot.hasCornerSectors())
        {
            // A moving contact can leave a real wet interval shorter than a
            // global-coordinate ULP. Integrate the certified restriction's
            // retained parameter interval, not two rounded global roots.
            // Keep the actual edge endpoints and the fitted face's z measure.
            for (const auto& piece : snapshot.partitionSegmentParametric(a, b))
                if (piece.negative)
                {
                    const long double dr = static_cast<long double>(piece.anchorEnd.x) - piece.anchorBegin.x;
                    const long double dt = static_cast<long double>(piece.tEnd) - piece.tBegin;
                    const long double radii = 2.0L * piece.anchorBegin.x
                        + dr * (static_cast<long double>(piece.tBegin) + piece.tEnd);
                    result += face.outwardSign * static_cast<long double>(face.center.y) * pi * std::abs(dr * dt * radii);
                }
            continue;
        }
        for (const auto& piece : snapshot.partitionSegment(a, b))
            if (piece.negative)
                result += face.outwardSign * static_cast<long double>(face.center.y) * pi
                    * std::abs((static_cast<long double>(piece.end.x)-piece.begin.x) * (piece.end.x+piece.begin.x));
    }
    require(std::isfinite(result) && result > 0.0L, "invalid complete fitted negative-phase volume");
    return result;
}

static BoundaryMoments boundaryMomentsImpl(const QuadtreeMesh& mesh, const NativeLeafDomain& domain,
    const LevelSetHermiteGeometry::Snapshot& snapshot, const MappedFacePartitions* mapped)
{
    validate(mesh, domain, snapshot);
    require(mapped || !snapshot.usesNativeLeaves(), "native boundary moments require retained mapped face partitions");
    BoundaryMoments result;
    const long double nan = std::numeric_limits<long double>::quiet_NaN();
    result.negativePhaseByNativeFace.assign(domain.faces(mesh).size(), nan);
    result.negativePhaseByCell.assign(mesh.cellIndexNum(), nan);
    result.capacityByCell.assign(mesh.cellIndexNum(), nan);
    for (int k : domain.cellIndexByRow())
        result.negativePhaseByCell[k] = result.capacityByCell[k] = 0.0L;
    for (const auto& face : domain.faces(mesh))
    {
        if (face.negativeCell < 0 && face.positiveCell < 0) continue;
        try
        {
            const auto& edge = face.horizontal ? mesh.hEdges(face.edgeLabel) : mesh.vEdges(face.edgeLabel);
            auto a = mesh.verts(edge.edgeToVert[0]).coord(), b = mesh.verts(edge.edgeToVert[1]).coord();
            if (face.horizontal ? a.x > b.x : a.y > b.y) std::swap(a, b);
            long double wet = 0.0L, whole = 0.0L;
            // Increasing z has right co-normal +r; increasing r has -z.
            const double orientation = face.horizontal ? -1.0 : 1.0;
            for (int half = 0; half < 2; ++half)
            {
                const auto begin = half == 0 ? a : face.center;
                const auto end = half == 0 ? face.center : b;
                const LevelSetHermiteGeometry::ParametricPiece complete{begin, end, 0.0, 1.0, false};
                const long double moment = orientation*CartesianBoundaryIntegral::yVolume(complete,domain.axisymmetric());
                whole += moment;
                // dr=0 makes this F's entire trace zero, regardless of phase.
                if (begin.x == end.x) continue;
                if (mapped)
                {
                    wet += orientation*CartesianBoundaryIntegral::yVolume(mapped->halfFace(face.nativeRow,half),domain.axisymmetric());
                    continue;
                }
                const int sign = snapshot.certifiedSegmentSign(begin, end);
                if (sign != 0) { if (sign < 0) wet += moment; continue; }
                for (const auto& piece : snapshot.partitionSegmentParametric(begin, end,
                    LevelSetHermiteGeometry::SegmentPartitionPolicy::CertifyUnresolvedGridPieces))
                    if (piece.negative) wet += orientation*CartesianBoundaryIntegral::yVolume(piece,domain.axisymmetric());
            }
            require(std::isfinite(wet) && std::isfinite(whole), "nonfinite native boundary moment");
            result.negativePhaseByNativeFace[face.nativeRow] = wet;
            if (face.negativeCell >= 0)
            {
                result.negativePhaseByCell[face.negativeCell] += wet;
                result.capacityByCell[face.negativeCell] += whole;
            }
            if (face.positiveCell >= 0)
            {
                result.negativePhaseByCell[face.positiveCell] -= wet;
                result.capacityByCell[face.positiveCell] -= whole;
            }
            result.exteriorNegativeMoment += face.outwardSign*wet;
        }
        catch (const std::exception& error)
        {
            std::ostringstream message;
            message << std::setprecision(17) << "[InterfacePhaseGeometry] boundary moment face=" << face.nativeRow
                << ", center=(" << face.center.x << ',' << face.center.y << "): " << error.what();
            throw std::runtime_error(message.str());
        }
    }
    for (int k : domain.cellIndexByRow())
    {
        require(std::isfinite(result.negativePhaseByCell[k]) && std::isfinite(result.capacityByCell[k])
            && result.capacityByCell[k] > 0.0L, "nonfinite moment or nonpositive actual control-volume capacity");
        result.summedCellNegativeMoment += result.negativePhaseByCell[k];
    }
    result.assemblyClosure = result.summedCellNegativeMoment-result.exteriorNegativeMoment;
    require(std::isfinite(result.exteriorNegativeMoment) && std::isfinite(result.summedCellNegativeMoment)
        && std::isfinite(result.assemblyClosure), "nonfinite assembled boundary moment");
    return result;
}

BoundaryMoments boundaryMoments(const QuadtreeMesh& mesh, const NativeLeafDomain& domain,
    const LevelSetHermiteGeometry::Snapshot& snapshot)
{ return boundaryMomentsImpl(mesh,domain,snapshot,nullptr); }
BoundaryMoments boundaryMoments(const QuadtreeMesh& mesh, const NativeLeafDomain& domain,
    const MappedFacePartitions& partitions)
{
    partitions.validate(mesh,domain);
    return boundaryMomentsImpl(mesh,domain,partitions.source(),&partitions);
}
}
