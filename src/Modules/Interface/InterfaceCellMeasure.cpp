#include "Modules/Interface/InterfaceCellMeasure.h"
#include "Modules/Interface/Detail/InterfaceNativeMeasure.h"
#include "Modules/Interface/Detail/HermiteGalsRegionMeasure.h"
#include "Core/Mesh/NativeLeafDomain.h"
#include "Core/Geometry/ActiveHermiteVolumeSensitivity.h"
#include "Core/Geometry/Detail/ExactExpansion.h"
#include "Core/Geometry/LevelSetMappedInterface.h"
#include "Core/Geometry/LevelSetMappedSegment.h"
#include "Core/Geometry/CartesianBoundaryIntegral.h"
#include "Core/Mesh/CartesianMeasure.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/UniformLeafGrid.h"

#include <algorithm>
#include <cmath>
#include <exception>
#include <iomanip>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <tuple>

namespace InterfaceCellMeasure
{
namespace
{
namespace Exact = Quantum::Geometry::Detail::ExactExpansion;
using Snapshot = LevelSetHermiteGeometry::Snapshot;
using Root = LevelSetHermiteGeometry::SegmentRootIncidence;
using Quad = std::array<Vector2D, 4>;
void require(bool ok, const char* reason)
{ if (!ok) throw std::runtime_error(std::string("[InterfaceCellMeasure] ")+reason); }
struct FaceResolutionError : std::runtime_error
{
    struct Request { int nativeRow,half; std::size_t rootIndex; };
    std::vector<Request> requests;
    explicit FaceResolutionError(const std::string& message) : std::runtime_error(message) {}
};
struct Position { Exact::Expansion x, y; };
Exact::Expansion difference(double a, double b)
{ return Exact::subtract(Exact::fromDouble(a), Exact::fromDouble(b)); }
int side(Vector2D a, Vector2D b, const Position& p)
{
    return Exact::sign(Exact::subtract(
        Exact::multiply(difference(b.x,a.x), Exact::subtract(p.y,Exact::fromDouble(a.y))),
        Exact::multiply(difference(b.y,a.y), Exact::subtract(p.x,Exact::fromDouble(a.x)))));
}
Position actual(Vector2D p) { return {Exact::fromDouble(p.x),Exact::fromDouble(p.y)}; }
bool contains(const Quad& q, const Position& point)
{
    for (int k=0;k<4;++k) if (side(q[k],q[(k+1)%4],point)<0) return false;
    return true;
}
struct Cell { int index; std::array<Quad,4> parts; };
struct LocalControlVolumes
{
    std::vector<Cell> cells;
    std::set<int> faces;
};

template<bool Quartic=false>
auto integrateWithRefinementRequests(const Snapshot& source,
    const LevelSetMappedInterface::Graph& graph, const std::vector<const LevelSetMappedSegment*>& cuts,
    const std::vector<int>& faces, const std::vector<LevelSetMappedInterface::Region>& regions,
    bool axisymmetric = true)
{
    require(cuts.size()==2*faces.size(), "half-face inventory size mismatch");
    try
    {
        if constexpr(Quartic) return LevelSetMappedInterface::integrateQuartic(source,graph,cuts,regions,axisymmetric);
        else return LevelSetMappedInterface::integrate(source,graph,cuts,regions,axisymmetric);
    }
    catch(const LevelSetMappedInterface::RootResolutionError& error)
    {
        FaceResolutionError request(std::string("[InterfaceCellMeasure] ")+error.what());
        // Multiple closed-cell incidences may name one root. Deduplicate its
        // original chart/representative before replacing any half-face query.
        std::set<std::tuple<int,int,std::size_t,double>> seen;
        for(const auto& item:error.requests)
        {
            require(item.segmentIndex<cuts.size(),"invalid Core face refinement request");
            const auto& roots=cuts[item.segmentIndex]->roots();
            require(item.rootIndex<roots.size(),"invalid Core root refinement request");
            const int face=faces[item.segmentIndex/2],half=static_cast<int>(item.segmentIndex%2);
            const auto& root=roots[item.rootIndex];
            if(seen.emplace(face,half,root.chart,root.parameter).second)
                request.requests.push_back({face,half,item.rootIndex});
        }
        require(!request.requests.empty(),"empty face refinement request");
        throw request;
    }
}

// Native CV geometry is shared by both scalar embeddings. Neither path
// replaces actual vertices, edge midpoints or cell centers by rectangles.
LocalControlVolumes localControlVolumes(const QuadtreeMesh& mesh, const NativeLeafDomain& domain,
    const UniformLeafGrid& grid, int cellI, int cellJ)
{
    const int nx=grid.nx(),ny=grid.ny();
    const auto vertex=[&](int i,int j) { return mesh.verts(grid.vertices()[i+(nx+1)*j]).coord(); };
    const auto faceRow=[&](int i,int j,bool horizontal) {
        return horizontal ? mesh.vEdgeIndexNum()+mesh.hEdges(grid.horizontal()[i+nx*j]).index()
            : mesh.vEdges(grid.vertical()[i+(nx+1)*j]).index();
    };
    LocalControlVolumes result;
    for(int j=std::max(0,cellJ-1);j<=std::min(ny-1,cellJ+1);++j)
        for(int i=std::max(0,cellI-1);i<=std::min(nx-1,cellI+1);++i)
        {
            const auto& native=mesh.cells(grid.cells()[i+nx*j]);
            if(domain.cellRowByIndex()[native.index()]<0) continue;
            const int west=faceRow(i,j,false),east=faceRow(i+1,j,false);
            const int south=faceRow(i,j,true),north=faceRow(i,j+1,true);
            for(int f:{west,east,south,north}) result.faces.insert(f);
            const auto sw=vertex(i,j),se=vertex(i+1,j),nw=vertex(i,j+1),ne=vertex(i+1,j+1);
            const auto w=domain.faces(mesh)[west].center,e=domain.faces(mesh)[east].center;
            const auto s=domain.faces(mesh)[south].center,n=domain.faces(mesh)[north].center,c=native.centerCoord();
            Cell cell{native.index(),{{Quad{{sw,s,c,w}},Quad{{s,se,e,c}},Quad{{w,c,n,nw}},Quad{{c,e,ne,n}}}}};
            for(const auto& q:cell.parts) for(int k=0;k<4;++k)
                require(side(q[k],q[(k+1)%4],actual(q[(k+2)%4]))>0,
                    "nonconvex/unresolved actual velocity subquad");
            result.cells.push_back(std::move(cell));
        }
    return result;
}
}

Measure evaluate(const QuadtreeMesh& mesh, const NativeLeafDomain& domain,
    const Snapshot& snapshot, const LevelSetInterfacePanels& panels)
{
    require(panels.usesSource(snapshot) && panels.axisymmetric(),
        "actual CV integration requires the same retained RZ panel source");
    Measure result;
    result.fixedBoundary = InterfacePhaseGeometry::boundaryMoments(mesh,domain,snapshot);
    const UniformLeafGrid grid(mesh,CartesianGridScope::UniformLeaves);
    const int nx=grid.nx();
    const auto h=grid.spacing(), lower=grid.minCoord();
    const long double nan=std::numeric_limits<long double>::quiet_NaN();
    result.freeBoundaryMomentByCell.assign(mesh.cellIndexNum(),nan);
    result.negativeVolumeByCell.assign(mesh.cellIndexNum(),nan);
    for (int k:domain.cellIndexByRow()) result.freeBoundaryMomentByCell[k]=0.0L;
    std::map<int,std::vector<Root>> inventories;
    for (int panelId=0;panelId<static_cast<int>(panels.panels().size());++panelId)
    {
        try
        {
            const auto& chart=panels.panelGraph(panelId);
            const double first=chart.coordinateIsX?chart.localBegin.x:chart.localBegin.y;
            const double last=chart.coordinateIsX?chart.localEnd.x:chart.localEnd.y;
            const auto coordinate=[&](Vector2D p) { return chart.coordinateIsX?p.x:p.y; };
            const auto overlay=localControlVolumes(mesh,domain,grid,chart.cellI,chart.cellJ);
            std::map<double,Vector2D> cuts;
            cuts.emplace(first,chart.localBegin); cuts.emplace(last,chart.localEnd);
            for(int row:overlay.faces)
            {
                const auto& face=domain.faces(mesh)[row];
                const auto& edge=face.horizontal?mesh.hEdges(face.edgeLabel):mesh.vEdges(face.edgeLabel);
                auto a=mesh.verts(edge.edgeToVert[0]).coord(),b=mesh.verts(edge.edgeToVert[1]).coord();
                if(face.horizontal?a.x>b.x:a.y>b.y) std::swap(a,b);
                for(int half=0;half<2;++half)
                {
                    const auto begin=half==0?a:face.center,end=half==0?face.center:b;
                    const int key=2*row+half;
                    auto found=inventories.find(key);
                    if(found==inventories.end())
                    {
                        std::vector<Root> roots;
                        if(snapshot.certifiedSegmentSign(begin,end)==0)
                        {
                            try { roots=snapshot.intersectSegmentWithProvenance(begin,end); }
                            catch(const std::exception& error)
                            {
                                std::ostringstream message;
                                message<<std::setprecision(17)<<"root inventory face="<<row<<", half="<<half
                                    <<", begin=("<<begin.x<<','<<begin.y<<"), end=("<<end.x<<','<<end.y
                                    <<"): "<<error.what();
                                throw std::runtime_error(message.str());
                            }
                        }
                        found=inventories.emplace(key,std::move(roots)).first;
                    }
                    for(const auto& root:found->second)
                    {
                        if(root.cellI!=chart.cellI || root.cellJ!=chart.cellJ) continue;
                        const double c=coordinate(root.localPoint);
                        // A cell may have a pure-phase zero corner outside its
                        // regular graph. It is not an interior graph cut.
                        if(c<std::min(first,last) || c>std::max(first,last)) continue;
                        const auto inserted=cuts.emplace(c,root.localPoint);
                        if(!inserted.second)
                            require(inserted.first->second.x==root.localPoint.x
                                && inserted.first->second.y==root.localPoint.y,
                                "distinct face roots collapse to one graph coordinate");
                    }
                }
            }
            std::vector<double> ordered;
            for(const auto& cut:cuts) ordered.push_back(cut.first);
            if(last<first) std::reverse(ordered.begin(),ordered.end());
            result.crossedPanels+=ordered.size()>2;
            const Vector2D origin{lower.x+chart.cellI*h.x,lower.y+chart.cellJ*h.y};
            for(std::size_t part=1;part<ordered.size();++part)
            {
                const auto samples=panels.subarcQuadrature(panelId,ordered[part-1],ordered[part]);
                require(!samples.empty(),"empty positive-length interface subarc");
                Arc arc;
                arc.panelIndex=panelId; arc.coordinateBegin=ordered[part-1]; arc.coordinateEnd=ordered[part];
                for(const auto& sample:samples)
                {
                    const auto local=sample.frame.cellLocalPoint;
                    // Same nominal chart as jetAtCellLocal, evaluated without
                    // rounding away a tiny offset before physical containment.
                    const Position point{
                        Exact::add(Exact::fromDouble(origin.x),Exact::multiply(Exact::fromDouble(h.x),Exact::fromDouble(local.x))),
                        Exact::add(Exact::fromDouble(origin.y),Exact::multiply(Exact::fromDouble(h.y),Exact::fromDouble(local.y)))};
                    int owner=-1;
                    for(const auto& cell:overlay.cells)
                    {
                        bool inside=false;
                        for(const auto& quad:cell.parts) inside|=contains(quad,point);
                        if(inside)
                        {
                            require(owner<0,"interface sample has multiple actual control-volume owners");
                            owner=cell.index;
                        }
                    }
                    require(owner>=0,"interface sample has no active actual control-volume owner");
                    require(arc.cellIndex<0 || arc.cellIndex==owner,
                        "subarc crosses an actual control-volume boundary missing from its root cuts");
                    arc.cellIndex=owner;
                    arc.freeBoundaryMoment+=static_cast<long double>(sample.physicalWeight)
                        *sample.frame.point.y*sample.frame.normal.y;
                }
                const int i=chart.cellI,j=chart.cellJ;
                const int positions[4]={i+(nx+1)*j,i+1+(nx+1)*j,i+1+(nx+1)*(j+1),i+(nx+1)*(j+1)};
                for(int corner=0;corner<4;++corner) arc.vertexIndices[corner]=mesh.verts(grid.vertices()[positions[corner]]).index();
                arc.sensitivity=ActiveHermiteVolumeSensitivity::assemble(i,j,h,samples);
                require(std::isfinite(arc.freeBoundaryMoment),"nonfinite free-boundary moment");
                result.freeBoundaryMomentByCell[arc.cellIndex]+=arc.freeBoundaryMoment;
                result.arcs.push_back(std::move(arc));
            }
        }
        catch(const std::exception& error)
        {
            std::ostringstream message;
            message<<"[InterfaceCellMeasure] panel="<<panelId<<", sourceCell="
                <<panels.panels()[panelId].cellIndex<<": "<<error.what();
            throw std::runtime_error(message.str());
        }
    }
    for(int k:domain.cellIndexByRow())
    {
        result.negativeVolumeByCell[k]=result.fixedBoundary.negativePhaseByCell[k]+result.freeBoundaryMomentByCell[k];
        require(std::isfinite(result.negativeVolumeByCell[k]),"nonfinite assembled negative-phase volume");
        result.summedNegativeVolume+=result.negativeVolumeByCell[k];
    }
    require(std::isfinite(result.summedNegativeVolume),"nonfinite total negative-phase volume");
    return result;
}

Measure evaluate(const QuadtreeMesh& mesh, const NativeLeafDomain& domain,
    const InterfacePhaseGeometry::MappedFacePartitions& partitions,
    const LevelSetInterfacePanels& panels)
{
    partitions.validate(mesh,domain);
    const auto& source=partitions.source();
    require(panels.usesSource(source) && panels.axisymmetric(),
        "mapped CV integration requires the same retained RZ panel source");
    Measure result;
    result.fixedBoundary=InterfacePhaseGeometry::boundaryMoments(mesh,domain,partitions);
    const UniformLeafGrid grid(mesh,CartesianGridScope::UniformLeaves);
    const int nx=grid.nx();
    const long double nan=std::numeric_limits<long double>::quiet_NaN();
    result.freeBoundaryMomentByCell.assign(mesh.cellIndexNum(),nan);
    result.negativeVolumeByCell.assign(mesh.cellIndexNum(),nan);
    std::vector<long double> compensation(mesh.cellIndexNum(),0.0L);
    for(int cell:domain.cellIndexByRow()) result.freeBoundaryMomentByCell[cell]=0.0L;
    for(int panelId=0;panelId<static_cast<int>(panels.panels().size());++panelId)
    {
        try
        {
            const auto& chart=panels.panelGraph(panelId);
            const auto overlay=localControlVolumes(mesh,domain,grid,chart.cellI,chart.cellJ);
            std::vector<LevelSetMappedInterface::Region> regions;
            for(const auto& cell:overlay.cells)
                for(const auto& quad:cell.parts) regions.push_back({cell.index,quad});
            std::vector<const LevelSetMappedSegment*> cuts;
            for(int face:overlay.faces)
                for(int half=0;half<2;++half) cuts.push_back(&partitions.halfFace(face,half));
            const LevelSetMappedInterface::Graph graph{chart.cellI,chart.cellJ,chart.coordinateIsX,
                {chart.beginRoot.side,chart.beginRoot.rootOrdinal},
                {chart.endRoot.side,chart.endRoot.rootOrdinal}};
            const auto integrated=integrateWithRefinementRequests(source,graph,cuts,
                std::vector<int>(overlay.faces.begin(),overlay.faces.end()),regions);
            require(!integrated.arcs.empty(),"mapped regular graph returned no positive subarc");
            result.crossedPanels+=integrated.arcs.size()>1;
            for(const auto& part:integrated.arcs)
            {
                Arc arc;
                arc.panelIndex=panelId; arc.cellIndex=part.regionId;
                require(arc.cellIndex>=0 && arc.cellIndex<mesh.cellIndexNum()
                    && domain.cellRowByIndex()[arc.cellIndex]>=0,"mapped arc owner is not an active native CV");
                arc.mappedCoordinate=true; arc.fromEnd=part.fromEnd;
                arc.coordinateOffset=part.offset; arc.coordinateLength=part.width;
                // Terminal display only. The integration and correction use
                // retained coordinates/assembled coefficients, never these.
                arc.coordinateBegin=static_cast<double>(part.fromEnd?1.0L-part.offset:part.offset);
                arc.coordinateEnd=static_cast<double>(part.fromEnd?1.0L-part.offset-part.width:part.offset+part.width);
                require(!part.requiresNativeDofConstraints && part.nativeCellIndex<0,
                    "native curve sensitivity requires hanging-node constraint assembly before volume correction");
                arc.freeBoundaryMoment=part.axialMoment; arc.sensitivity=part.sensitivity;
                const int i=chart.cellI,j=chart.cellJ;
                const int positions[4]={i+(nx+1)*j,i+1+(nx+1)*j,i+1+(nx+1)*(j+1),i+(nx+1)*(j+1)};
                for(int corner=0;corner<4;++corner) arc.vertexIndices[corner]=mesh.verts(grid.vertices()[positions[corner]]).index();
                require(std::isfinite(arc.freeBoundaryMoment),"nonfinite mapped free-boundary moment");
                const long double increment=arc.freeBoundaryMoment-compensation[arc.cellIndex];
                auto& sum=result.freeBoundaryMomentByCell[arc.cellIndex];
                const long double next=sum+increment;
                compensation[arc.cellIndex]=(next-sum)-increment;
                sum=next;
                result.arcs.push_back(std::move(arc));
            }
        }
        catch(const FaceResolutionError&) { throw; }
        catch(const std::exception& error)
        {
            std::ostringstream message;
            message<<"[InterfaceCellMeasure] mapped panel="<<panelId<<", sourceCell="
                <<panels.panels()[panelId].cellIndex<<": "<<error.what();
            throw std::runtime_error(message.str());
        }
    }
    long double totalCompensation=0.0L;
    for(int cell:domain.cellIndexByRow())
    {
        auto& volume=result.negativeVolumeByCell[cell];
        volume=result.fixedBoundary.negativePhaseByCell[cell]+result.freeBoundaryMomentByCell[cell];
        require(std::isfinite(volume),"nonfinite mapped negative-phase volume");
        const long double increment=volume-totalCompensation;
        const long double next=result.summedNegativeVolume+increment;
        totalCompensation=(next-result.summedNegativeVolume)-increment;
        result.summedNegativeVolume=next;
    }
    require(std::isfinite(result.summedNegativeVolume),"nonfinite mapped total negative-phase volume");
    return result;
}

template<bool Quartic,bool MomentsOnly=false>
static auto evaluateNative(const QuadtreeMesh& mesh, const NativeLeafDomain& domain,
    const InterfacePhaseGeometry::MappedFacePartitions& partitions,
    const std::vector<LevelSetMappedInterface::Graph>& nativeGraphs,
    const std::vector<int>* certifiedPhase = nullptr,
    const std::shared_ptr<const Snapshot>& localVolumeSource = {},
    const InterfacePhaseGeometry::BoundaryMoments* traceFixedBoundary = nullptr)
{
    partitions.validate(mesh,domain);
    const auto& source=partitions.source();
    require(!localVolumeSource || (Quartic && MomentsOnly && certifiedPhase && localVolumeSource.get()==&source),
        "local volume requires the same complete native Q4 source");
    require(source.usesNativeLeaves() && (certifiedPhase || domain.axisymmetric()),
        "native XY volume requires a complete state topology; legacy entries retain RZ scope");
    require(!certifiedPhase || certifiedPhase->size()==static_cast<std::size_t>(mesh.cellIndexNum()),
        "complete topology phase certificate has the wrong mesh size");
    require((source.reconstruction()==LevelSetHermiteGeometry::Reconstruction::EnrichedQuartic)==Quartic,
        "native volume sensitivity port does not match scalar data layout");
    std::conditional_t<MomentsOnly,Volumes,std::conditional_t<Quartic,QuarticMeasure,Measure>> result;
    result.fixedBoundary=traceFixedBoundary ? *traceFixedBoundary : InterfacePhaseGeometry::boundaryMoments(mesh,domain,partitions);
    const auto nan=std::numeric_limits<long double>::quiet_NaN();
    if constexpr(MomentsOnly) result.momentErrorEstimateByCell.assign(mesh.cellIndexNum(),nan);
    result.freeBoundaryMomentByCell.assign(mesh.cellIndexNum(),nan);
    result.negativeVolumeByCell.assign(mesh.cellIndexNum(),nan);
    std::vector<std::vector<int>> graphByCell(mesh.cellIndexNum());
    for(std::size_t i=0;i<nativeGraphs.size();++i)
    {
        const auto& graph=nativeGraphs[i];
        require(graph.cellI==-1 && graph.cellJ==-1 && graph.nativeCellIndex>=0
            && graph.nativeCellIndex<mesh.cellIndexNum()
            && domain.cellRowByIndex()[graph.nativeCellIndex]>=0, "native graph does not name an active source leaf");
                auto& inventory=graphByCell[graph.nativeCellIndex];
        require(certifiedPhase || inventory.empty(), "legacy entry accepts one native graph per scalar leaf");
        if (certifiedPhase) require((*certifiedPhase)[graph.nativeCellIndex]==0,
            "a cut leaf cannot also carry a pure-phase certificate");
        for (int old:inventory) {
            const auto& other=nativeGraphs[old];
            require(graph.localLower.x!=other.localLower.x || graph.localLower.y!=other.localLower.y ||
                graph.localUpper.x!=other.localUpper.x || graph.localUpper.y!=other.localUpper.y ||
                graph.begin.side!=other.begin.side || graph.begin.rootOrdinal!=other.begin.rootOrdinal ||
                graph.end.side!=other.end.side || graph.end.rootOrdinal!=other.end.rootOrdinal,
                "duplicate native graph identity");
        }
        inventory.push_back(static_cast<int>(i));
    }
    // Each cell owns its integrator and reads immutable scalar/half-face data.
    // Keep arcs, refinement requests and summation in the original cell order;
    // only the independent native integrations run on worker threads.
    using NativeArc=std::conditional_t<Quartic,QuarticArc,Arc>;
    struct CellResult {
        std::vector<NativeArc> arcs;
        std::vector<long double> momentErrors;
        int crossedPanels=0;
        std::exception_ptr failure;
    };
    const auto& cells=domain.cellIndexByRow();
    std::vector<CellResult> prepared(cells.size());
    const auto prepare=[&](std::size_t row) {
      auto& local=prepared[row];
      try {
        const int cell=cells[row];
        const auto& native=mesh.cells(mesh.cellIndexToLabel[cell]);
        const auto bounds=source.nativeCellBounds(cell);
        const auto lo=bounds[0],hi=bounds[1];
        auto& moment=result.freeBoundaryMomentByCell[cell];
        moment=0.0L;
        const auto& graphIds=graphByCell[cell];
        if(graphIds.empty())
        {
            // The same exact Bernstein certificate excludes all interior and
            // edge zeros, allowing only isolated zero vertices. It is not a
            // root-count heuristic and also handles declared pure-phase sectors.
            require(certifiedPhase ? std::abs((*certifiedPhase)[cell])==1 :
                source.closedCellSignExceptVertices(native.centerCoord(),15u)!=0,
                "native curve inventory omits an uncertified active leaf");
        }
        else
        {
            std::set<int> rows;
            for(int direction=0;direction<2;++direction)
                for(int part=0;part<2;++part)
                {
                    const int v=native.cellToVEdge[direction][part],h=native.cellToHEdge[direction][part];
                    if(v>=0) rows.insert(mesh.vEdges(v).index());
                    if(h>=0) rows.insert(mesh.vEdgeIndexNum()+mesh.hEdges(h).index());
                }
            std::vector<const LevelSetMappedSegment*> cuts;
            for(int row:rows)
            {
                const auto& face=domain.faces(mesh)[row];
                require(face.negativeCell==cell || face.positiveCell==cell,"native CV has an unrelated face");
                const auto& edge=face.horizontal?mesh.hEdges(face.edgeLabel):mesh.vEdges(face.edgeLabel);
                const auto a=mesh.verts(edge.edgeToVert[0]).coord(),b=mesh.verts(edge.edgeToVert[1]).coord();
                require(face.horizontal
                    ? (a.y==b.y && a.y==face.center.y && (a.y==lo.y || a.y==hi.y))
                    : (a.x==b.x && a.x==face.center.x && (a.x==lo.x || a.x==hi.x)),
                    "native rectangular CV has a non-collinear actual boundary");
                for(int half=0;half<2;++half) cuts.push_back(&partitions.halfFace(row,half));
            }
            const LevelSetMappedInterface::Region region{cell,{{lo,{hi.x,lo.y},hi,{lo.x,hi.y}}}};
            long double compensation=0.0L;
            for (int graphId:graphIds)
            {
                if constexpr(MomentsOnly)
                {
                    // Each graph lies in this actual native CV; artificial
                    // subchart boundaries are not additional control volumes.
                    // Fixed physical half-face moments still use partitions.
                    const auto integrated=LevelSetMappedInterface::nativeMoments(source,nativeGraphs[graphId],domain.axisymmetric());
                    require(std::isfinite(integrated.axial) && std::isfinite(integrated.errorEstimate),
                        "nonfinite native geometry-only moment");
                    const auto increment=integrated.axial-compensation,next=moment+increment;
                    compensation=(next-moment)-increment;moment=next;
                    local.momentErrors.push_back(integrated.errorEstimate);
                }
                else
                {
                    const auto integrated=integrateWithRefinementRequests<Quartic>(source,nativeGraphs[graphId],cuts,
                        std::vector<int>(rows.begin(),rows.end()),{region},domain.axisymmetric());
                    require(!integrated.arcs.empty(),"native regular graph returned no free arc");
                    local.crossedPanels+=integrated.arcs.size()>1;
                    for(const auto& part:integrated.arcs)
                    {
                        require(part.regionId==cell && part.nativeCellIndex==cell && !part.requiresNativeDofConstraints
                            && !part.independentSensitivity.empty(),"native CV curve lacks constrained scalar sensitivity");
                        std::conditional_t<Quartic,QuarticArc,Arc> arc;
                        arc.panelIndex=graphId; arc.cellIndex=cell; arc.nativeSourceCell=cell;
                        arc.mappedCoordinate=true; arc.fromEnd=part.fromEnd;
                        arc.coordinateOffset=part.offset; arc.coordinateLength=part.width;
                        arc.coordinateBegin=static_cast<double>(part.fromEnd?1.0L-part.offset:part.offset);
                        arc.coordinateEnd=static_cast<double>(part.fromEnd?1.0L-part.offset-part.width:part.offset+part.width);
                        arc.freeBoundaryMoment=part.axialMoment;
                        arc.independentSensitivity=part.independentSensitivity;
                        require(std::isfinite(arc.freeBoundaryMoment),"nonfinite native free-boundary moment");
                        const long double increment=arc.freeBoundaryMoment-compensation,next=moment+increment;
                        compensation=(next-moment)-increment; moment=next;
                        local.arcs.push_back(std::move(arc));
                    }
                }
            }
        }
        auto& volume=result.negativeVolumeByCell[cell];
        volume=result.fixedBoundary.negativePhaseByCell[cell]+moment;
        if constexpr(Quartic && MomentsOnly) if(localVolumeSource && !graphIds.empty()) {
            namespace Region=Quantum::Interface::Detail::RegionMeasure;
            std::vector<LevelSetMappedInterface::Graph> graphs;
            for(int graphId:graphIds)graphs.push_back(nativeGraphs[graphId]);
            const Region::Box box{lo,hi,cell};
            const auto measured=Region::integrateCell(localVolumeSource,box,graphs,domain.axisymmetric());
            const auto capacity=Region::capacity(lo,hi,domain.axisymmetric());
            for(long double value:{measured.axialLiquid,measured.radialLiquid,measured.axialGas,measured.radialGas,measured.error})
                require(std::isfinite(value),"nonfinite local native moment");
            const auto defect=std::max({std::abs(measured.axialLiquid-measured.radialLiquid),std::abs(measured.axialGas-measured.radialGas),
                std::abs(measured.axialLiquid+measured.axialGas-capacity),std::abs(measured.radialLiquid+measured.radialGas-capacity)})/capacity;
            // Same local Green/capacity check as the existing region consumer.
            require(defect<=1e-10L,"local native Green/capacity ledger defect");
            const bool axial=hi.y-lo.y<=hi.x-lo.x;
            const auto liquid=axial?measured.axialLiquid:measured.radialLiquid,gas=axial?measured.axialGas:measured.radialGas;
            require(liquid>=0 && liquid<=capacity && gas>=0 && gas<=capacity,"local native volume outside strict capacity; no clipping");
            require(measured.error>=0 && measured.error<=1e-10L*capacity,"local native moment estimate exceeds tolerance");
            volume=Region::complement(capacity,liquid,gas)[0];
            local.momentErrors={measured.error};
            // Keep global fixed/free moments untouched for their shared-face
            // ledger. Their rounded sum can differ from this local volume.
        }
        require(std::isfinite(volume),"nonfinite native negative-phase volume");
      } catch(...) {local.failure=std::current_exception();}
    };
    Exact::validateEnvironment(); // Check the caller before entering an OpenMP team.
    if(!cells.empty())prepare(0);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 8) if(nativeGraphs.size()>16)
#endif
    for(std::ptrdiff_t row=1;row<static_cast<std::ptrdiff_t>(cells.size());++row)
        prepare(static_cast<std::size_t>(row));
    long double totalCompensation=0.0L;
    for(std::size_t row=0;row<cells.size();++row) {
        auto& local=prepared[row];
        if(local.failure)std::rethrow_exception(local.failure);
        if constexpr(MomentsOnly) {
            for(long double error:local.momentErrors)result.momentErrorEstimate+=error;
            auto& cellError=result.momentErrorEstimateByCell[cells[row]];cellError=0;
            for(long double error:local.momentErrors)cellError+=error;
        } else {
            result.crossedPanels+=local.crossedPanels;
            for(auto& arc:local.arcs)result.arcs.push_back(std::move(arc));
        }
        const auto volume=result.negativeVolumeByCell[cells[row]];
        const long double increment=volume-totalCompensation,next=result.summedNegativeVolume+increment;
        totalCompensation=(next-result.summedNegativeVolume)-increment; result.summedNegativeVolume=next;
    }
    require(std::isfinite(result.summedNegativeVolume),"nonfinite total native negative-phase volume");
    return result;
}

Measure evaluate(const QuadtreeMesh& mesh, const NativeLeafDomain& domain,
    const InterfacePhaseGeometry::MappedFacePartitions& partitions,
    const std::vector<LevelSetMappedInterface::Graph>& nativeGraphs)
{ return evaluateNative<false>(mesh,domain,partitions,nativeGraphs); }

QuarticMeasure evaluateQuartic(const QuadtreeMesh& mesh, const NativeLeafDomain& domain,
    const InterfacePhaseGeometry::MappedFacePartitions& partitions,
    const std::vector<LevelSetMappedInterface::Graph>& nativeGraphs)
{ return evaluateNative<true>(mesh,domain,partitions,nativeGraphs); }

namespace
{
using CompleteNativeInventory = Detail::NativeInventory;
Measure evaluate(const QuadtreeMesh& mesh, const NativeLeafDomain& domain,
    const InterfacePhaseGeometry::MappedFacePartitions& partitions, const CompleteNativeInventory& inventory)
{ return evaluateNative<false>(mesh,domain,partitions,inventory.graphs,inventory.phase); }
QuarticMeasure evaluateQuartic(const QuadtreeMesh& mesh, const NativeLeafDomain& domain,
    const InterfacePhaseGeometry::MappedFacePartitions& partitions, const CompleteNativeInventory& inventory)
{ return evaluateNative<true>(mesh,domain,partitions,inventory.graphs,inventory.phase); }
}

namespace
{
using InterfaceCellMeasure::evaluate;
using InterfaceCellMeasure::evaluateQuartic;
template<bool Quartic=false,class Graphs>
auto evaluateRefinedImpl(const QuadtreeMesh& mesh, const NativeLeafDomain& domain,
    InterfacePhaseGeometry::MappedFacePartitions& partitions,
    const Graphs& panels, std::size_t maxRootRefinements)
{
    std::size_t refinements=0;
    for(;;)
    {
        try
        {
            auto result=[&] {
                if constexpr(Quartic) return evaluateQuartic(mesh,domain,partitions,panels);
                else return evaluate(mesh,domain,partitions,panels);
            }();
            result.rootRefinements=refinements;
            return result;
        }
        catch(const FaceResolutionError& error)
        {
            for(const auto& request:error.requests)
            {
                if(refinements==maxRootRefinements)
                    throw std::runtime_error(std::string("[InterfaceCellMeasure] root refinement work limit exhausted: ")+error.what());
                partitions.refineRoot(request.nativeRow,request.half,request.rootIndex);
                ++refinements;
            }
            // No partly computed measure escapes. Fixed-boundary moments,
            // free arcs and their sensitivities are all rebuilt together.
        }
    }
}
}

Measure evaluateRefined(const QuadtreeMesh& mesh, const NativeLeafDomain& domain,
    InterfacePhaseGeometry::MappedFacePartitions& partitions,
    const LevelSetInterfacePanels& panels, std::size_t maxRootRefinements)
{ return evaluateRefinedImpl(mesh,domain,partitions,panels,maxRootRefinements); }

Measure evaluateRefined(const QuadtreeMesh& mesh, const NativeLeafDomain& domain,
    InterfacePhaseGeometry::MappedFacePartitions& partitions,
    const std::vector<LevelSetMappedInterface::Graph>& nativeGraphs, std::size_t maxRootRefinements)
{ return evaluateRefinedImpl(mesh,domain,partitions,nativeGraphs,maxRootRefinements); }

Measure evaluateRefined(const QuadtreeMesh& mesh, const NativeLeafDomain& domain,
    InterfacePhaseGeometry::MappedFacePartitions& partitions,
    const LevelSetInterfacePanels::NativeGraphChain& chain, std::size_t maxRootRefinements)
{
    require(chain.usesSource(partitions.source()) && chain.axisymmetric(),
        "native RZ volume requires the same retained chain source and axisymmetric measure");
    return evaluateRefinedImpl(mesh,domain,partitions,chain.graphs(),maxRootRefinements);
}

QuarticMeasure evaluateQuarticRefined(const QuadtreeMesh& mesh, const NativeLeafDomain& domain,
    InterfacePhaseGeometry::MappedFacePartitions& partitions,
    const LevelSetInterfacePanels::NativeGraphChain& chain, std::size_t maxRootRefinements)
{
    require(chain.usesSource(partitions.source()) && chain.axisymmetric(),
        "native quartic RZ volume requires the same retained chain source and axisymmetric measure");
    return evaluateRefinedImpl<true>(mesh,domain,partitions,chain.graphs(),maxRootRefinements);
}
namespace Detail {
Measure evaluateRefined(const QuadtreeMesh& mesh, const NativeLeafDomain& domain,
    InterfacePhaseGeometry::MappedFacePartitions& partitions, const NativeInventory& inventory,
    std::size_t maxRootRefinements)
{ return evaluateRefinedImpl(mesh,domain,partitions,inventory,maxRootRefinements); }
QuarticMeasure evaluateQuarticRefined(const QuadtreeMesh& mesh, const NativeLeafDomain& domain,
    InterfacePhaseGeometry::MappedFacePartitions& partitions, const NativeInventory& inventory,
    std::size_t maxRootRefinements)
{ return evaluateRefinedImpl<true>(mesh,domain,partitions,inventory,maxRootRefinements); }
Volumes evaluateVolumes(const QuadtreeMesh& mesh, const NativeLeafDomain& domain,
    const InterfacePhaseGeometry::MappedFacePartitions& partitions, const NativeInventory& inventory)
{
    if(partitions.source().reconstruction()==LevelSetHermiteGeometry::Reconstruction::EnrichedQuartic)
        return evaluateNative<true,true>(mesh,domain,partitions,inventory.graphs,inventory.phase);
    return evaluateNative<false,true>(mesh,domain,partitions,inventory.graphs,inventory.phase);
}
}
DiscoveredNativeTraceVolumes measureNativeTraceVolumes(const QuadtreeMesh& mesh, const NativeLeafDomain& domain,
    const InterfacePhaseGeometry::MappedFacePartitions& partitions,
    const LevelSetInterfacePanels::NativeComponents& topology,
    const std::vector<Vector2D>& boundaryContacts,
    const std::vector<LevelSetInterfacePanels::SymmetryPlane>& symmetryPlanes)
{
    partitions.validate(mesh,domain);
    const auto& source=partitions.source();
    require(source.usesNativeLeaves() && source.reconstruction()==LevelSetHermiteGeometry::Reconstruction::EnrichedQuartic,
        "native trace volume requires original native Q4 source");
    // Synchronous borrowed alias only: the cache already owns this exact
    // immutable allocation. Re-discovery certifies completeness of the public
    // topology value; callers cannot supply a modified phase mask/arc list.
    const std::shared_ptr<const Snapshot> borrowed(&source,[](const Snapshot*){});
    const auto certified=LevelSetInterfacePanels::discoverNativeTraceComponents(mesh,borrowed,domain.axisymmetric(),boundaryContacts,symmetryPlanes);
    require(topology.cellPhase==certified.cellPhase && topology.curves.size()==certified.curves.size(),
        "trace topology differs from complete same-source certificate");
    std::vector<LevelSetMappedInterface::Graph> ordinary, sides;
    for(std::size_t c=0;c<topology.curves.size();++c) {
        const auto& chain=topology.curves[c];const auto& reference=certified.curves[c];
        require(chain.usesSource(source) && chain.axisymmetric()==domain.axisymmetric()
            && chain.graphs().size()==reference.graphs().size() && chain.closed()==reference.closed()
            && chain.beginEndpointKind()==reference.beginEndpointKind() && chain.endEndpointKind()==reference.endEndpointKind(),
            "native trace topology changed source/measure/component");
        for(std::size_t g=0;g<chain.graphs().size();++g) {
            const auto& a=chain.graphs()[g];const auto& b=reference.graphs()[g];
            require(a.nativeCellIndex==b.nativeCellIndex && a.cellI==b.cellI && a.cellJ==b.cellJ
                && a.coordinateIsX==b.coordinateIsX && a.controlledIntegration==b.controlledIntegration
                && a.localLower.x==b.localLower.x && a.localLower.y==b.localLower.y
                && a.localUpper.x==b.localUpper.x && a.localUpper.y==b.localUpper.y
                && a.begin.side==b.begin.side && a.begin.rootOrdinal==b.begin.rootOrdinal
                && a.end.side==b.end.side && a.end.rootOrdinal==b.end.rootOrdinal
                && static_cast<bool>(a.nativeSide)==static_cast<bool>(b.nativeSide),"incomplete/foreign native trace graph");
            if(a.nativeSide) {
                const auto& s=*a.nativeSide;const auto& r=*b.nativeSide;
                require(s.usesSource(source) && s.ownerCell()==r.ownerCell() && s.side()==r.side()
                    && s.edgeLabel()==r.edgeLabel() && s.horizontal()==r.horizontal()
                    && s.beginVertex()==r.beginVertex() && s.endVertex()==r.endVertex()
                    && s.begin().x==r.begin().x && s.begin().y==r.begin().y
                    && s.end().x==r.end().x && s.end().y==r.end().y
                    && s.minus().nativeCellIndex==r.minus().nativeCellIndex && s.minus().phase==r.minus().phase
                    && s.plus().nativeCellIndex==r.plus().nativeCellIndex && s.plus().phase==r.plus().phase,
                    "foreign native trace side support");
                sides.push_back(a);
            } else ordinary.push_back(a);
        }
    }
    DiscoveredNativeTraceVolumes result;auto& ledger=result.ledger;
    InterfacePhaseGeometry::BoundaryMoments fixed;
    const auto nan=std::numeric_limits<long double>::quiet_NaN();
    fixed.negativePhaseByNativeFace.assign(domain.faces(mesh).size(),nan);
    fixed.negativePhaseByCell.assign(mesh.cellIndexNum(),nan);fixed.capacityByCell.assign(mesh.cellIndexNum(),nan);
    for(auto* v:{&ledger.fixedRadialByCell,&ledger.freeRadialByCell,&ledger.axialVolumeByCell,&ledger.radialVolumeByCell,&ledger.radialCapacityByCell})
        v->assign(mesh.cellIndexNum(),nan);
    for(int cell:domain.cellIndexByRow()) {
        fixed.negativePhaseByCell[cell]=fixed.capacityByCell[cell]=0;
        ledger.fixedRadialByCell[cell]=ledger.freeRadialByCell[cell]=ledger.radialCapacityByCell[cell]=0;
    }
    for(const auto& face:domain.faces(mesh)) {
        if(face.negativeCell<0 && face.positiveCell<0)continue;
        const long double orientation=face.horizontal?-1:1;
        long double axial=0,radial=0,wholeAxial=0,wholeRadial=0;
        for(int half=0;half<2;++half) {
            const auto& trace=partitions.traceHalfFace(face.nativeRow,half);
            require(trace.usesSource(source) && trace.normalAxis()==(face.horizontal?1:0),"trace half changed source/normal");
            const auto a=trace.begin(),b=trace.end();
            const long double dx=static_cast<long double>(b.x)-a.x,dy=static_cast<long double>(b.y)-a.y;
            require(face.horizontal ? (dy==0 && dx>0) : (dx==0 && dy>0),"trace half is not an ordered actual face");
            const auto radialMoment=[&](bool fromEnd,long double offset,long double length) {
                const long double delta=fromEnd?-dx:dx;
                const auto x=std::fma(delta,offset,static_cast<long double>(fromEnd?b.x:a.x));
                const auto mean=domain.axisymmetric()
                    ? static_cast<long double>(CartesianMeasure::twoPi)*.5L*(x*x+delta*x*length+delta*delta*length*length/3)
                    : x+delta*length*.5L;
                const auto value=orientation*dy*length*mean;
                require(std::isfinite(value),"nonfinite retained radial boundary moment");return value;
            };
            wholeRadial+=radialMoment(false,0,1);
            wholeAxial+=orientation*CartesianBoundaryIntegral::yVolume({a,b,0,1,false},domain.axisymmetric());
            if(trace.coincidentRanges().empty()) {
                if(dx!=0)axial+=orientation*CartesianBoundaryIntegral::yVolume(trace.ordinary(),domain.axisymmetric());
                for(const auto& range:trace.ordinaryRanges())if(range.negative)radial+=radialMoment(range.fromEnd,range.offset,range.length);
            } else {
                // Actual halves lie in one subedge. Complete native topology
                // proves a regular whole zero edge, not mixed value support.
                require(trace.ordinaryRanges().empty(),"regular zero half has mixed ordinary support");
                ++ledger.zeroHalfFaces;
            }
        }
        fixed.negativePhaseByNativeFace[face.nativeRow]=axial;
        for(int side=0;side<2;++side) {
            const int cell=side?face.positiveCell:face.negativeCell;if(cell<0)continue;
            const int sign=side?-1:1;
            fixed.negativePhaseByCell[cell]+=sign*axial;fixed.capacityByCell[cell]+=sign*wholeAxial;
            ledger.fixedRadialByCell[cell]+=sign*radial;ledger.radialCapacityByCell[cell]+=sign*wholeRadial;
        }
        fixed.exteriorNegativeMoment+=face.outwardSign*axial;ledger.exteriorFixedRadial+=face.outwardSign*radial;
    }
    // Preserve the original ordinary axial arithmetic (also on typed caches).
    if(ledger.zeroHalfFaces==0)fixed=InterfacePhaseGeometry::boundaryMoments(mesh,domain,partitions);
    else {
        for(int cell:domain.cellIndexByRow())fixed.summedCellNegativeMoment+=fixed.negativePhaseByCell[cell];
        fixed.assemblyClosure=fixed.summedCellNegativeMoment-fixed.exteriorNegativeMoment;
    }
    result.measure=evaluateNative<true,true>(mesh,domain,partitions,ordinary,&topology.cellPhase,borrowed,&fixed);
    for(const auto& graph:ordinary)ledger.freeRadialByCell[graph.nativeCellIndex]+=LevelSetMappedInterface::nativeMoments(source,graph,domain.axisymmetric()).radial;
    for(const auto& graph:sides) {
        ++ledger.sideGraphs;const auto& support=*graph.nativeSide;int liquid=-1;
        for(const auto* side:{&support.minus(),&support.plus()})if(side->available() && side->phase<0) {
            require(liquid<0,"regular native side has two liquid owners");liquid=side->nativeCellIndex;
        }
        // A physical/selected boundary may expose only gas. Keep its geometry
        // but do not invent the absent exterior liquid CV or a volume moment.
        if(liquid<0)continue;
        require(domain.cellRowByIndex().at(liquid)>=0 && topology.cellPhase.at(liquid)==-1,
            "native side lacks a certified actual liquid owner");
        const auto moment=LevelSetMappedInterface::nativeMoments(source,graph,domain.axisymmetric());
        require(std::isfinite(moment.axial) && std::isfinite(moment.radial) && std::isfinite(moment.errorEstimate),"nonfinite native side moment");
        result.measure.freeBoundaryMomentByCell[liquid]+=moment.axial;ledger.freeRadialByCell[liquid]+=moment.radial;
        result.measure.momentErrorEstimateByCell[liquid]+=moment.errorEstimate;result.measure.momentErrorEstimate+=moment.errorEstimate;
        ++ledger.sideLiquidOwners;
    }
    result.measure.summedNegativeVolume=0;long double radialFixedSum=0,totalCompensation=0;
    for(int cell:domain.cellIndexByRow()) {
        const long double cap=domain.cellVolumes()[domain.cellRowByIndex()[cell]];
        require(std::isfinite(cap) && cap>0,"invalid native trace capacity");
        const auto capacityDefect=std::max(std::abs(fixed.capacityByCell[cell]-cap),std::abs(ledger.radialCapacityByCell[cell]-cap))/cap;
        require(std::isfinite(capacityDefect) && capacityDefect<=1e-10L,"trace two-direction capacity ledger defect");
        ledger.maxRelativeCapacityDefect=std::max(ledger.maxRelativeCapacityDefect,capacityDefect);
        const auto axial=fixed.negativePhaseByCell[cell]+result.measure.freeBoundaryMomentByCell[cell];
        const auto radial=ledger.fixedRadialByCell[cell]+ledger.freeRadialByCell[cell];
        const int phase=topology.cellPhase[cell];auto& volume=result.measure.negativeVolumeByCell[cell];
        if(phase!=0)volume=phase<0?cap:0;
        const auto defect=std::max({std::abs(axial-radial),std::abs(axial-volume),std::abs(radial-volume)})/cap;
        require(std::isfinite(volume) && volume>=0 && volume<=cap,"native trace volume outside strict capacity; no clipping");
        require(std::isfinite(defect) && defect<=1e-10L,"trace fixed/free Green ledger defect");
        ledger.maxRelativeGreenDefect=std::max(ledger.maxRelativeGreenDefect,defect);
        ledger.axialVolumeByCell[cell]=axial;ledger.radialVolumeByCell[cell]=radial;
        radialFixedSum+=ledger.fixedRadialByCell[cell];
        const auto increment=volume-totalCompensation,next=result.measure.summedNegativeVolume+increment;
        totalCompensation=(next-result.measure.summedNegativeVolume)-increment;result.measure.summedNegativeVolume=next;
    }
    ledger.fixedRadialClosure=radialFixedSum-ledger.exteriorFixedRadial;
    result.cellPhase=topology.cellPhase;result.components=topology.curves.size();result.graphs=ordinary.size()+sides.size();
    result.faceRootIncidences=partitions.rootIncidenceCount();return result;
}

DiscoveredNativeTraceVolumes discoverNativeTraceVolumes(const QuadtreeMesh& mesh, const NativeLeafDomain& domain,
    std::shared_ptr<const Snapshot> source, const std::vector<Vector2D>& contacts,
    const std::vector<LevelSetInterfacePanels::SymmetryPlane>& symmetry)
{
    require(static_cast<bool>(source),"null native trace volume source");
    InterfacePhaseGeometry::validateLayout(mesh,domain,*source);
    const auto topology=LevelSetInterfacePanels::discoverNativeTraceComponents(mesh,source,domain.axisymmetric(),contacts,symmetry);
    const auto partitions=InterfacePhaseGeometry::MappedFacePartitions::prepareNativeTrace(mesh,domain,source);
    return measureNativeTraceVolumes(mesh,domain,partitions,topology,contacts,symmetry);
}

template<bool Bubble>
static auto discoverCompleteVolumes(const QuadtreeMesh& mesh, const NativeLeafDomain& domain,
    std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> source,
    const std::vector<Vector2D>& boundaryContacts,
    const std::vector<LevelSetInterfacePanels::SymmetryPlane>& symmetryPlanes)
{
    require(static_cast<bool>(source), "null complete volume source");
    InterfacePhaseGeometry::validateLayout(mesh,domain,*source);
    require(source->usesNativeLeaves(), "complete volume discovery requires native leaves");
    auto topology=LevelSetInterfacePanels::discoverNativeComponents(mesh,source,
        domain.axisymmetric(),boundaryContacts,symmetryPlanes);
    InterfacePhaseGeometry::MappedFacePartitions partitions(mesh,domain,source);
    std::vector<LevelSetMappedInterface::Graph> graphs;
    for(const auto& chain:topology.curves) {
        require(chain.usesSource(*source) && chain.axisymmetric()==domain.axisymmetric(),
            "discovered volume chain changed source or measure");
        graphs.insert(graphs.end(),chain.graphs().begin(),chain.graphs().end());
    }
    std::conditional_t<Bubble,DiscoveredBubbleVolumes,DiscoveredVolumes> result;
    if constexpr(Bubble) {
        require(source->reconstruction()==LevelSetHermiteGeometry::Reconstruction::EnrichedQuartic,
            "bubble sensitivity requires native quartic scalar");
        const CompleteNativeInventory inventory{graphs,&topology.cellPhase};
        const auto sensitivity=evaluateRefinedImpl<true>(mesh,domain,partitions,inventory,256);
        result.centerDerivative.assign(mesh.cellIndexNum(),std::numeric_limits<long double>::quiet_NaN());
        for(int c:domain.cellIndexByRow())result.centerDerivative[c]=0;
        using Sensitivity=Quantum::Geometry::ActiveQuarticHermiteSensitivity;
        for(const auto& arc:sensitivity.arcs)for(const auto& entry:arc.independentSensitivity)
            if(entry.kind==Sensitivity::Kind::CellCenter) {
                require(entry.index==arc.cellIndex && entry.component==0,"bubble derivative escaped its own cell");
                result.centerDerivative[entry.index]+=entry.coefficient;
            }
    }
    if(source->reconstruction()==LevelSetHermiteGeometry::Reconstruction::EnrichedQuartic)
        result.measure=evaluateNative<true,true>(mesh,domain,partitions,graphs,&topology.cellPhase,source);
    else result.measure=evaluateNative<false,true>(mesh,domain,partitions,graphs,&topology.cellPhase);
    result.cellPhase=std::move(topology.cellPhase);
    result.components=topology.curves.size(); result.graphs=graphs.size();
    result.faceRootIncidences=partitions.rootIncidenceCount();
    return result;
}

DiscoveredVolumes discoverVolumes(const QuadtreeMesh& mesh, const NativeLeafDomain& domain,
    std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> source,
    const std::vector<Vector2D>& boundaryContacts,
    const std::vector<LevelSetInterfacePanels::SymmetryPlane>& symmetryPlanes)
{ return discoverCompleteVolumes<false>(mesh,domain,std::move(source),boundaryContacts,symmetryPlanes); }

DiscoveredBubbleVolumes discoverBubbleVolumes(const QuadtreeMesh& mesh, const NativeLeafDomain& domain,
    std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> source,
    const std::vector<Vector2D>& boundaryContacts,
    const std::vector<LevelSetInterfacePanels::SymmetryPlane>& symmetryPlanes)
{ return discoverCompleteVolumes<true>(mesh,domain,std::move(source),boundaryContacts,symmetryPlanes); }

}
