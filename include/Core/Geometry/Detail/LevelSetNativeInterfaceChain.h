#pragma once

#include "Core/Geometry/Detail/LevelSetInterfaceChainOrder.h"
#include "Core/Geometry/Detail/LevelSetMappedInterface.h"
#include "Core/Geometry/LevelSetInterfacePanels.h"
#include "Core/Mesh/CartesianLeafMesh.h"
#include "Core/Geometry/InterfaceGeometryResolution.h"
#include <array>
#include <limits>
#include <map>
#include <memory>
#include <iomanip>
#include <sstream>
#include <type_traits>

namespace LAMP::Geometry::Detail
{
namespace
{
namespace NativeInterfaceChain
{
namespace M = MappedHermiteSegment;
namespace I = MappedHermiteInterface;
using Kind = LevelSetInterfacePanels::EndpointKind;
using Graph = LevelSetMappedInterface::Graph;
struct SideArc
{
    int owner=-1,side=-1,edgeLabel=-1,beginVertex=-1,endVertex=-1;
    bool horizontal=false;
    Vector2D begin,end;
    LevelSetMappedInterface::NativeSideSupport::Side minus,plus;
};
inline void require(bool valid, const char *reason)
{
    if (!valid)
        throw std::runtime_error(std::string("[LevelSetNativeInterfaceChain] ") + reason);
}
struct Result
{
    std::vector<Graph> graphs;
    std::vector<SideArc> sideSupports;
    std::vector<int> mapping;
    std::vector<std::vector<int>> cellGraphs;
    std::vector<std::array<int,2>> rootNodes;
    std::vector<Vector2D> boundaryContacts;
    std::array<int,2> endpointSymmetryPlanes{{-1,-1}};
    Vector2D minimumSpacing;
    int negativeRegion = -1;
    bool closed = true;
    Kind begin = Kind::None, end = Kind::None;
};
struct Discovery
{
    std::vector<Result> curves;
    std::vector<LevelSetInterfacePanels::NativeNegativeRegion> negativeRegions;
    std::vector<int> cellPhase;
};
template<class Query> struct EdgeData
{
    Vector2D begin, end;
    int vertices[2];
    std::vector<int> cells, nodes;
    std::unique_ptr<Query> query;
    bool zero=false;
    int normalSign=0,graph=-1,edgeLabel=-1;
    LevelSetMappedInterface::NativeSideSupport::Side minus,plus;
};
struct Node
{
    I::Point point;
    int vertex = -1;
    Kind kind = Kind::None;
    std::vector<int> cells;
    I::Point lower=point,upper=point;
    int symmetryPlane = -1;
};
inline bool same(Vector2D a, Vector2D b)
{
    return a.x == b.x && a.y == b.y;
}
template<class Query> inline bool onEdge(Vector2D p, const EdgeData<Query> &edge)
{
    return I::orientation(edge.begin, edge.end, I::point(p)) == 0 && p.x >= edge.begin.x &&
           p.x <= edge.end.x && p.y >= edge.begin.y && p.y <= edge.end.y;
}
// A terminal polynomial layout may still describe a piecewise native
// scalar. Its common signed partial is a regularity certificate, never
// permission to integrate all native pieces as one polynomial.
template<class Scalar,class=void> struct RequiresPieceCharts:std::false_type {};
template<class Scalar> struct RequiresPieceCharts<Scalar,
    std::void_t<decltype(Scalar::requiresPieceCharts)>>
    :std::integral_constant<bool,Scalar::requiresPieceCharts> {};
inline void validateMacroCharts(
    const std::vector<LevelSetHermiteScalar::QuarticCell::RegularChartBox>& boxes,int cell)
{
    if(boxes.empty())throw InterfaceGeometryResolution(InterfaceGeometryResolution::Reason::CriticalCell,
        cell,"C2 macro lacks a certified complete piece-chart partition");
    // Validate complete tensor coverage even when every piece is pure and
    // therefore no graph construction follows to perform its own check.
    require(same(boxes.front().lower,{0,0}),"invalid macro partition origin");
    std::vector<double> xCuts{0},yCuts{0};
    for(const auto& box:boxes){
        if(box.lower.y!=0)break;
        require(box.lower.x==xCuts.back() && box.upper.x>box.lower.x,"nonconforming macro chart columns");
        xCuts.push_back(box.upper.x);
    }
    const std::size_t nx=xCuts.size()-1;
    require(nx>=2 && xCuts.back()==1 && boxes.size()%nx==0,"incomplete macro chart columns");
    const std::size_t ny=boxes.size()/nx;
    for(std::size_t y=0;y<ny;++y){
        const auto& first=boxes[y*nx];
        require(first.lower.y==yCuts.back() && first.upper.y>first.lower.y,"nonconforming macro chart rows");
        yCuts.push_back(first.upper.y);
        for(std::size_t x=0;x<nx;++x){
            const auto& box=boxes[y*nx+x];
            require(same(box.lower,{xCuts[x],yCuts[y]}) && same(box.upper,{xCuts[x+1],yCuts[y+1]}),
                "macro charts are not a complete conforming tensor product");
        }
    }
    require(ny>=2 && yCuts.back()==1,"incomplete macro chart rows");
    bool xHalf=false,yHalf=false;
    for(const auto& box:boxes){
        require(box.lower.x>=0 && box.lower.y>=0 && box.upper.x<=1 && box.upper.y<=1
            && box.lower.x<box.upper.x && box.lower.y<box.upper.y,"invalid macro chart extent");
        require(!(box.lower.x<.5 && box.upper.x>.5)
            && !(box.lower.y<.5 && box.upper.y>.5),"macro chart crosses an actual half-piece seam");
        require(box.phase==-1 || box.phase==0 || box.phase==1,"invalid macro chart phase certificate");
        require(box.phase!=0 || box.partial[0]!=0 || box.partial[1]!=0,
            "macro zero-containing chart lacks a strict Cartesian partial");
        xHalf=xHalf || box.lower.x==.5 || box.upper.x==.5;
        yHalf=yHalf || box.lower.y==.5 || box.upper.y==.5;
    }
    require(xHalf && yHalf,"macro chart partition lost a native half-piece seam");
}
template<class Scalar> inline bool resolvedPartial(const Scalar& source,int cell,bool legacy,int axis=-1)
{
    if (source.resolvedCartesianPartial(cell,axis)) return true;
    return !legacy && source.resolvedZeroSetPartial(cell,axis);
}

// A regular even contact has two incidences on its inward side and none
// on the other side. Keep its one root identity in both edge inventories.
template<class Scalar> inline void applyNativeIncidences(const Scalar& source,const Graph& chart,
    const std::array<std::vector<int>,4>& inventories,std::vector<int>& local,
    std::vector<int>* inwardContacts=nullptr)
{
    using Integrator=I::NativeIntegratorFor<Scalar>;
    std::map<int,int> counts;
    for(int side=0;side<4;++side) {
        if(inventories[side].empty())continue;
        const auto incidence=Integrator::chartEdgeIncidences(source,chart,side);
        require(incidence.size()==inventories[side].size(),"contact certificate lost complete edge-root ordinals");
        for(std::size_t k=0;k<incidence.size();++k) {
            if(incidence[k]<0)continue; // a corner's other side can prove its inward direction
            const int id=inventories[side][k];
            const auto found=counts.emplace(id,incidence[k]);
            require(found.second || found.first->second==incidence[k],"corner contact has inconsistent inward certificates");
        }
    }
    std::vector<int> retained;
    for(int id:local) {
        const auto found=counts.find(id);
        if(found==counts.end())throw InterfaceGeometryResolution(InterfaceGeometryResolution::Reason::RootIncidence,
            chart.nativeCellIndex,"contact lacks a nonzero original normal derivative certificate");
        for(int k=0;k<found->second;++k)retained.push_back(id);
        if(inwardContacts && found->second==2)inwardContacts->push_back(id);
    }
    local=std::move(retained);
}

// Partition geometry inside ONE native scalar cell. Macro boxes must also
// retain every real half-piece seam. Internal chart edges have
// one root inventory shared by their two chart owners; actual cell boundaries
// retain the existing mesh-edge node identities. Only zero arcs are emitted:
// these internal edges never enter phase-volume or physical-flux ledgers.
template<class Scalar> inline void appendChartGraphs(const Scalar& source,int cell,
    const std::vector<LevelSetHermiteScalar::QuarticCell::RegularChartBox>& boxes,
    const std::array<std::vector<int>,4>& boundaryNodes,std::vector<Node>& nodes,
    std::vector<Graph>& graphs,std::vector<std::array<int,2>>& connections)
{
    using Integrator=I::NativeIntegratorFor<Scalar>;
    const std::size_t firstGraph=graphs.size();
    const auto bounds=source.leafBounds(cell);
    const auto origin=I::point(bounds[0]);
    const I::Point width{{M::ratio(M::difference(bounds[1].x,bounds[0].x),M::number(1)),
        M::ratio(M::difference(bounds[1].y,bounds[0].y),M::number(1))}};
    const auto physical=[&](const I::Point& p) {
        return I::Point{{I::sum(origin[0],I::product(width[0],p[0])),
            I::sum(origin[1],I::product(width[1],p[1]))}};
    };
    require(!boxes.empty() && same(boxes.front().lower,{0,0}),"invalid chart partition origin");
    std::vector<double> xCuts{0},yCuts{0};
    for(const auto& box:boxes) {
        if(box.lower.y!=0) break;
        require(box.lower.x==xCuts.back() && box.upper.x>box.lower.x,"nonconforming chart columns");
        xCuts.push_back(box.upper.x);
    }
    const int nx=static_cast<int>(xCuts.size())-1;
    require(nx>=2 && xCuts.back()==1 && boxes.size()%nx==0,"invalid complete chart partition");
    const int ny=static_cast<int>(boxes.size()/nx);
    for(int y=0;y<ny;++y) {
        const auto& first=boxes[y*nx];
        require(first.lower.y==yCuts.back() && first.upper.y>first.lower.y,"nonconforming chart rows");
        yCuts.push_back(first.upper.y);
        for(int x=0;x<nx;++x) {
            const auto& box=boxes[y*nx+x];
            require(same(box.lower,{xCuts[x],yCuts[y]}) && same(box.upper,{xCuts[x+1],yCuts[y+1]}),
                "chart partition is not a conforming tensor product");
        }
    }
    require(ny>=2 && yCuts.back()==1,"incomplete chart partition extent");
    std::map<std::array<int,3>,std::vector<int>> edgeNodes;
    std::map<std::array<int,2>,int> vertexNodes;
    std::map<int,std::array<int,2>> nodeVertices;
    const auto boundaryIdentity=[&](int side,const I::Point& lower,const I::Point& upper) {
        int found=-1;
        for (int id:boundaryNodes[side]) {
            bool intersects=true;
            for (int a=0;a<2;++a) intersects=intersects &&
                M::compare(lower[a],nodes[id].upper[a])<=0 && M::compare(nodes[id].lower[a],upper[a])<=0;
            if (!intersects) continue;
            if (found>=0 && found!=id) throw InterfaceGeometryResolution(
                InterfaceGeometryResolution::Reason::RootIncidence,cell,
                "chart edge root has ambiguous actual-boundary enclosure identity");
            found=id;
        }
        if (found<0) throw InterfaceGeometryResolution(InterfaceGeometryResolution::Reason::RootIncidence,
            cell,"chart edge root is missing from its actual native boundary inventory");
        return found;
    };
    for (std::size_t b=0;b<boxes.size();++b) {
        const auto& box=boxes[b];
        if (box.phase) continue;
        const int ix=static_cast<int>(b)%nx,iy=static_cast<int>(b)/nx;
        Graph chart;chart.nativeCellIndex=cell;chart.localLower=box.lower;chart.localUpper=box.upper;
        chart.controlledIntegration=true;
        std::array<std::vector<int>,4> inventories;
        std::vector<int> local;
        for (int side=0;side<4;++side) {
            const bool vertical=side==0 || side==2;
            const int fixed=vertical?ix+(side==2):iy+(side==3),start=vertical?iy:ix;
            const std::array<int,3> key{{vertical?1:0,fixed,start}};
            auto found=edgeNodes.find(key);
            if (found==edgeNodes.end()) {
                std::vector<int> inventory;
                const auto roots=Integrator::chartEdgeRoots(source,chart,side);
                for (const auto& root:roots) {
                    I::Point p,lo,hi;
                    p[vertical?0:1]=lo[vertical?0:1]=hi[vertical?0:1]=M::constant(vertical?xCuts[fixed]:yCuts[fixed]);
                    p[vertical?1:0]=root.parameter;lo[vertical?1:0]=root.lower;hi[vertical?1:0]=root.upper;
                    const auto physicalPoint=physical(p),physicalLo=physical(lo),physicalHi=physical(hi);
                    int outer=-1;
                    if (M::compare(p[0],M::constant(0))==0) outer=0;
                    else if (M::compare(p[0],M::constant(1))==0) outer=2;
                    else if (M::compare(p[1],M::constant(0))==0) outer=1;
                    else if (M::compare(p[1],M::constant(1))==0) outer=3;
                    const auto& tangentCuts=vertical?yCuts:xCuts;
                    const int endpoint=M::compare(root.parameter,M::constant(tangentCuts[start]))==0?0:
                        M::compare(root.parameter,M::constant(tangentCuts[start+1]))==0?1:-1;
                    const std::array<int,2> vertex{{vertical?fixed:start+endpoint,vertical?start+endpoint:fixed}};
                    int id=outer>=0?boundaryIdentity(outer,physicalLo,physicalHi):-1;
                    if (endpoint>=0) {
                        const auto existing=vertexNodes.find(vertex);
                        if (existing!=vertexNodes.end()) {
                            require(id<0 || id==existing->second,"chart vertex disagrees with actual boundary identity");
                            id=existing->second;
                        }
                    }
                    if (id<0) {
                        id=static_cast<int>(nodes.size());nodes.push_back({physicalPoint});
                        nodes.back().lower=physicalLo;nodes.back().upper=physicalHi;
                    }
                    if (endpoint>=0) {vertexNodes.emplace(vertex,id);nodeVertices.emplace(id,vertex);}
                    inventory.push_back(id);
                }
                found=edgeNodes.emplace(key,std::move(inventory)).first;
            }
            inventories[side]=found->second;
            for (int id:inventories[side])
                if (std::find(local.begin(),local.end(),id)==local.end()) local.push_back(id);
        }
        if constexpr(Scalar::dofCount==25) applyNativeIncidences(source,chart,inventories,local);
        else local.erase(std::remove_if(local.begin(),local.end(),[&](int id) {
            const auto found=nodeVertices.find(id);
            if (found==nodeVertices.end()) return false;
            const auto vertex=found->second;
            if ((vertex[0]!=ix && vertex[0]!=ix+1) || (vertex[1]!=iy && vertex[1]!=iy+1)) return false;
            const auto gradient=source.localGradientSigns(cell,{xCuts[vertex[0]],yCuts[vertex[1]]});
            const int xIn=(vertex[0]==ix?1:-1)*gradient[0],yIn=(vertex[1]==iy?1:-1)*gradient[1];
            return xIn!=0 && xIn==yIn;
        }),local.end());
        if (local.empty()) continue; // strict partial excludes an interior loop
        if (local.size()%2) throw InterfaceGeometryResolution(InterfaceGeometryResolution::Reason::RootIncidence,
            cell,"unpaired chart crossings or unresolved tangency");
        const int chosen=box.partial[1]?0:1;
        chart.coordinateIsX=chosen==0;
        std::sort(local.begin(),local.end(),[&](int a,int b) {
            return M::compare(nodes[a].point[chosen],nodes[b].point[chosen])<0;
        });
        const auto endpoint=[&](int id) {
            for (int side=0;side<4;++side) {
                const auto found=std::find(inventories[side].begin(),inventories[side].end(),id);
                if (found!=inventories[side].end()) return LevelSetMappedInterface::Endpoint{
                    side,static_cast<int>(found-inventories[side].begin())};
            }
            throw std::runtime_error("[LevelSetNativeInterfaceChain] missing chart endpoint");
        };
        for (std::size_t k=0;k<local.size();k+=2) {
            std::array<int,2> ends{{local[k],local[k+1]}};
            Graph graph=chart;graph.begin=endpoint(ends[0]);graph.end=endpoint(ends[1]);
            if (Integrator::graphOrientation(source,graph)<0) {
                std::swap(graph.begin,graph.end);std::swap(ends[0],ends[1]);
            }
            const int index=static_cast<int>(graphs.size());graphs.push_back(graph);connections.push_back(ends);
            for (int id:ends) nodes[id].cells.push_back(index);
        }
    }
    // The original-cell negative-fragment builder joins curves to actual
    // boundary intervals. A wholly interior closed cycle additionally needs
    // hole/nesting ownership; keep that case explicit rather than assigning
    // its two boundary cycles to different physical liquids.
    std::vector<bool> visited(graphs.size()-firstGraph,false);
    for (std::size_t seed=firstGraph;seed<graphs.size();++seed) {
        if (visited[seed-firstGraph]) continue;
        bool boundary=false;std::vector<int> pending{static_cast<int>(seed)};
        while (!pending.empty()) {
            const int g=pending.back();pending.pop_back();
            if (visited[g-firstGraph]) continue;
            visited[g-firstGraph]=true;
            for (int id:connections[g]) {
                for (const auto& side:boundaryNodes)
                    boundary=boundary || std::find(side.begin(),side.end(),id)!=side.end();
                for (int next:nodes[id].cells)
                    if (next>=static_cast<int>(firstGraph)) pending.push_back(next);
            }
        }
        if (!boundary) throw InterfaceGeometryResolution(InterfaceGeometryResolution::Reason::CriticalCell,cell,
            "closed subcell cycle requires physical refinement for negative-fragment ownership");
    }
}

// Native address adapter for the existing panel topology rules. Edge roots
// come from one authoritative actual subedge query. Coarse sides concatenate
// that inventory; only exact shared VERTEX identity is deduplicated. Complete
// side ordinals retain roots excluded from a particular cell's incidence.
template<class Scalar> inline Discovery discoverAll(const QuadtreeMesh &mesh, const Scalar &source, bool axisymmetric,
                       const std::vector<Vector2D> &contacts, bool legacy = false,
                       const std::vector<LevelSetInterfacePanels::SymmetryPlane>& symmetryPlanes={},bool nativeTrace=false)
{
    validateCartesianLeafMesh(mesh);
    if constexpr(RequiresPieceCharts<Scalar>::value)
        require(!legacy && !axisymmetric,"C2 macro discovery requires complete nonlegacy XY geometry");
    require(!nativeTrace || (!legacy && Scalar::dofCount==25),"native trace discovery requires complete Q4 topology");
    require(source.cellCount() == static_cast<std::size_t>(mesh.cellIndexNum()) &&
                source.vertexCount() == static_cast<std::size_t>(mesh.vertIndexNum()),
            "mesh/source size mismatch");
    require(same(source.minCoord(), mesh.minCoord) && same(source.maxCoord(), mesh.maxCoord),
            "mesh/source domain mismatch");
    require(!legacy || contacts.size() <= 2, "one open chain accepts at most two declared endpoints");
    const bool axisContact = axisymmetric && source.radialAxisSymmetry() && source.minCoord().x == 0;
    for (std::size_t k=0;k<symmetryPlanes.size();++k) {
        const auto& plane=symmetryPlanes[k];
        require((plane.normalAxis==0 || plane.normalAxis==1) && std::isfinite(plane.coordinate),
            "invalid symmetry plane");
        const auto coordinate=[&](Vector2D p){return plane.normalAxis==0?p.x:p.y;};
        require(plane.coordinate==coordinate(mesh.minCoord) || plane.coordinate==coordinate(mesh.maxCoord),
            "symmetry plane must be an actual outer domain section");
        require(!(axisContact && plane.normalAxis==0 && plane.coordinate==0),
            "radial axis already has its own symmetry identity");
        for (std::size_t j=0;j<k;++j) require(!(plane==symmetryPlanes[j]),"duplicate symmetry plane");
        for (auto p:contacts) require(coordinate(p)!=plane.coordinate,
            "fixed contact and symmetry section cannot share an endpoint declaration");
    }
    using Query=M::BasicQuery<M::BasicNativeCellLayout<Scalar>>;
    std::vector<EdgeData<Query>> edges;
    std::map<std::pair<bool, int>, int> edgeIds;
    std::vector<std::array<std::vector<int>, 4>> sides(source.cellCount());
    std::vector<int> signs(source.cellCount(), 0), cellPhase(source.cellCount(), 0);
    std::vector<int> sidePhase(source.cellCount(),0);
    using ChartBox=LevelSetHermiteScalar::QuarticCell::RegularChartBox;
    std::vector<std::vector<ChartBox>> chartBoxes(source.cellCount());
    bool allActive = true;
    Vector2D minimumSpacing{std::numeric_limits<double>::infinity(),
                            std::numeric_limits<double>::infinity()};
    for (int c = 0; c < mesh.cellIndexNum(); ++c)
    {
        if (!source.cellActive(c))
        {
            allActive = false;
            continue;
        }
        const auto &cell = mesh.cells(mesh.cellIndexToLabel[c]);
        const auto bounds = source.leafBounds(c);
        minimumSpacing.x = std::min(minimumSpacing.x,bounds[1].x-bounds[0].x);
        minimumSpacing.y = std::min(minimumSpacing.y,bounds[1].y-bounds[0].y);
        require(same(bounds[0], mesh.verts(cell.cellToVert[0]).coord()) &&
                    same(bounds[1], mesh.verts(cell.cellToVert[2]).coord()),
                "mesh/source leaf bounds mismatch");
        const auto vertices = source.cellVertices(c);
        for (int k = 0; k < 4; ++k)
            require(vertices[k] == mesh.verts(cell.cellToVert[k]).index(),
                    "mesh/source native vertex identity mismatch");
        signs[c] = source.cellSignExceptVertices(c, 0);
        const bool uncertainCell =
            signs[c] == 0 && source.cellSignExceptVertices(c, 15) == 0;
        if (uncertainCell)
        {
            // One common signed partial on the complete zero set also
            // excludes two crossings of a transverse section (their
            // crossing directions would alternate). Off-interface extrema
            // do not require chart subdivision. Preserve this certificate
            // before partitioning: a regular straight zero can coincide
            // with an artificial dyadic edge.
            if constexpr(RequiresPieceCharts<Scalar>::value){
                // Even a full-native signed partial cannot erase a real
                // polynomial knot. Graphs and their quadrature stay on one
                // actual terminal piece; no whole-cell polynomial fallback.
                chartBoxes[c]=source.regularChartPartition(c);
                validateMacroCharts(chartBoxes[c],c);
                const int phase=chartBoxes[c].front().phase;
                if(phase!=0 && std::all_of(chartBoxes[c].begin(),chartBoxes[c].end(),
                    [phase](const ChartBox& box){return box.phase==phase;})){
                    // Complete outward piece bounds certify the original
                    // native cell, including internal seams. Share this
                    // phase with edge queries and the negative-region ledger.
                    signs[c]=cellPhase[c]=phase;
                }
            }else if (!legacy && !resolvedPartial(source,c,false))
                chartBoxes[c]=source.regularChartPartition(c);
            if (chartBoxes[c].empty() && !resolvedPartial(source,c,legacy)) {
            std::ostringstream detail;
            detail<<std::setprecision(17)<<"cell lacks a regular Cartesian chart (hidden loop or critical zero)"
                <<"; cell="<<c<<", bounds=("<<bounds[0].x<<','<<bounds[0].y<<")-("
                <<bounds[1].x<<','<<bounds[1].y<<')';
            if constexpr (RequiresPieceCharts<Scalar>::value) {
                // Actual completed physical data; a macro has no whole-native
                // 25-column polynomial. Failure diagnostics cannot invent one.
                detail<<", completed_macro_physical_49_data=[";
                for(double value:source.physicalCellDofs(c))detail<<value<<',';
                detail<<']';
            } else if constexpr (Scalar::dofCount==25) {
                // Failure-only local data for numerical diagnosis. At AMR
                // constraints these are rounded display data, not an exact
                // replacement scalar or permission to accept the cell.
                detail<<", approximate_local_Q4_data=[";
                for (double value:source.mappedCellLocalDofs(c)) detail<<value<<',';
                detail<<']';
            } else {
                detail<<", approximate_corner_dofs_SW_SE_NW_NE=[";
                for (int y=0;y<2;++y) for (int x=0;x<2;++x)
                    for (double value:source.mappedCellNodeDofs(c,x,y)) detail<<value<<',';
                detail<<']';
            }
            throw InterfaceGeometryResolution(InterfaceGeometryResolution::Reason::CriticalCell, c,
                detail.str());
            }
        }
        for (int side = 0; side < 4; ++side)
            for (int part = 0; part < 2; ++part)
            {
                const bool horizontal = side == 1 || side == 3;
                const int label =
                    horizontal ? cell.cellToHEdge[side / 2][part] : cell.cellToVEdge[side / 2][part];
                if (label < 0)
                    continue;
                const auto key = std::make_pair(horizontal, label);
                auto found = edgeIds.find(key);
                if (found == edgeIds.end())
                {
                    const auto &edge = horizontal ? mesh.hEdges(label) : mesh.vEdges(label);
                    const auto &a = mesh.verts(edge.edgeToVert[0]);
                    const auto &b = mesh.verts(edge.edgeToVert[1]);
                    EdgeData<Query> data{a.coord(), b.coord(), {a.index(), b.index()}};
                    data.edgeLabel=label;
                    if (std::make_pair(data.end.x, data.end.y) < std::make_pair(data.begin.x, data.begin.y))
                    {
                        std::swap(data.begin, data.end);
                        std::swap(data.vertices[0], data.vertices[1]);
                    }
                    const int id = static_cast<int>(edges.size());
                    edges.push_back(std::move(data));
                    found = edgeIds.emplace(key, id).first;
                }
                sides[c][side].push_back(found->second);
                edges[found->second].cells.push_back(c);
            }
    }
    require(!contacts.empty() || allActive,
            "an active-subdomain interface requires explicit boundary contacts");
    for (const auto& plane:symmetryPlanes) {
        bool found=false;
        for (const auto& edge:edges) {
            const auto coordinate=[&](Vector2D p){return plane.normalAxis==0?p.x:p.y;};
            if (edge.cells.size()!=1 || coordinate(edge.begin)!=plane.coordinate || coordinate(edge.end)!=plane.coordinate) continue;
            found=true;
            if constexpr(RequiresPieceCharts<Scalar>::value) {
                require(source.normalTraceIdenticallyZero(edge.cells.front(),edge.begin,edge.end,plane.normalAxis),
                    "symmetry section requires identically zero scalar normal derivative");
            } else {
                const auto normal=source.restrictPhysical(edge.cells.front(),edge.begin,edge.end,
                    plane.normalAxis==0?1:0,plane.normalAxis==1?1:0);
                for (const auto& coefficient:normal.numerator)
                    require(ExactExpansion::sign(coefficient)==0,
                        "symmetry section requires identically zero scalar normal derivative");
            }
        }
        require(found,"symmetry plane has no active boundary trace");
    }
    for (std::size_t k = 0; k < contacts.size(); ++k)
    {
        const auto p = contacts[k];
        require(std::isfinite(p.x) && std::isfinite(p.y) && (!axisymmetric || p.x > 0),
                "invalid non-axis boundary declaration");
        bool boundary = false;
        for (const auto &edge : edges)
            boundary = boundary || (edge.cells.size() == 1 && onEdge(p, edge));
        require(boundary, "declared endpoint is not on an actual active boundary");
        for (std::size_t j = 0; j < k; ++j)
            require(!same(p, contacts[j]), "duplicate boundary declaration");
    }
    std::vector<Node> nodes;
    std::vector<int> vertexNodes(source.vertexCount(), -1), matched(contacts.size(), -1);
    for (auto &edge : edges)
    {
        require(edge.cells.size() == 1 || edge.cells.size() == 2, "invalid actual subedge ownership");
        // Keep all four inventories on a cut leaf, including pure-phase sides;
        // the graph's complete side ordinals refer to these actual subedges.
        if (std::all_of(edge.cells.begin(), edge.cells.end(), [&](int c) { return signs[c] != 0; }))
            continue;
        if(nativeTrace)edge.query=std::make_unique<Query>(source,edge.begin,edge.end,M::BasicNativeCellLayout<Scalar>(source),true);
        else edge.query = std::make_unique<Query>(source, edge.begin, edge.end);
        if(nativeTrace && !edge.query->zeroTraces.empty()) {
            require(edge.query->phases.empty(),"actual zero subedge has mixed positive-length value support");
            const int normal=edge.begin.x==edge.end.x?0:1;
            if constexpr(Scalar::dofCount==25)for(int c:edge.cells) {
                const auto bounds=source.leafBounds(c);
                const double position=normal?edge.begin.y:edge.begin.x,lo=normal?bounds[0].y:bounds[0].x,hi=normal?bounds[1].y:bounds[1].x;
                require(position==lo || position==hi,"zero edge is not an actual incident leaf side");
                if(!source.resolvedZeroSetPartial(c,normal))throw InterfaceGeometryResolution(
                    InterfaceGeometryResolution::Reason::CriticalCell,c,"native zero side lacks a complete signed normal zero-set certificate");
                int derivative=0;
                if constexpr(RequiresPieceCharts<Scalar>::value) {
                    // The complete zero-set certificate above proves this
                    // continuous normal derivative is nonzero on the WHOLE
                    // actual side. Its exact endpoint sign therefore fixes
                    // both pieces without inventing a whole-face polynomial.
                    const auto p=source.restrictPhysical(c,edge.begin,edge.begin,normal==0,normal==1).numerator;
                    derivative=ExactExpansion::sign(p[0]);
                } else {
                    const auto p=source.restrictPhysical(c,edge.begin,edge.end,normal==0,normal==1).numerator;
                    const auto roots=PolynomialSegmentRoots::isolateEndpointRelative(p);
                    require(!roots.identicallyZero && roots.roots.empty(),"native zero side lacks a strict first normal derivative");
                    derivative=ExactExpansion::sign(p[0]);
                }
                require(derivative!=0 && (edge.normalSign==0 || edge.normalSign==derivative),"native zero side has inconsistent incident normals");
                edge.normalSign=derivative;
                const int inward=position==lo?1:-1,phase=inward*derivative;
                require(sidePhase[c]==0 || sidePhase[c]==phase,"one native leaf has conflicting zero-side phase certificates");
                sidePhase[c]=phase;
                auto& side=position==hi?edge.minus:edge.plus;
                require(!side.available(),"native zero side has overlapping incident support");side={c,phase};
            }
            for(int endpoint=0;endpoint<2;++endpoint) {
                const int vertex=edge.vertices[endpoint];int id=vertexNodes[vertex];const auto p=I::point(endpoint?edge.end:edge.begin);
                if(id<0){id=static_cast<int>(nodes.size());nodes.push_back({p,vertex});vertexNodes[vertex]=id;}
                else require(I::same(nodes[id].point,p),"zero subedge vertex disagrees with ordinary endpoint incidence");
                edge.nodes.push_back(id);
            }
            edge.zero=true;continue; // No fictitious finite value-root inventory.
        }
        const auto &roots = edge.query->roots;
        for (std::size_t r = 0; r < roots.size();)
        {
            const auto &root = roots[r];
            std::size_t next = r + 1;
            while (next < roots.size() && M::compare(root.parameter, roots[next].parameter) == 0)
                ++next;
            const int endpoint = M::compare(root.parameter, M::constant(0)) == 0   ? 0
                                 : M::compare(root.parameter, M::constant(1)) == 0 ? 1
                                                                                   : -1;
            const int vertex = endpoint >= 0 ? edge.vertices[endpoint] : -1;
            require(!legacy || !contacts.empty() || vertex < 0,
                    "default interface path does not allow grid-vertex contact");
            int id = vertex >= 0 ? vertexNodes[vertex] : -1;
            if (id < 0)
            {
                I::Point point;
                for (int a = 0; a < 2; ++a)
                    point[a] = I::sum(M::constant(a == 0 ? edge.begin.x : edge.begin.y),
                                      I::product(M::ratio(M::difference(a == 0 ? edge.end.x : edge.end.y,
                                                                        a == 0 ? edge.begin.x : edge.begin.y),
                                                          M::number(1)),
                                                 root.parameter));
                id = static_cast<int>(nodes.size());
                nodes.push_back({point, vertex});
                for (int a=0;a<2;++a) {
                    const auto origin=M::constant(a==0?edge.begin.x:edge.begin.y);
                    const auto width=M::ratio(M::difference(a==0?edge.end.x:edge.end.y,
                        a==0?edge.begin.x:edge.begin.y),M::number(1));
                    nodes.back().lower[a]=I::sum(origin,I::product(width,root.lower));
                    nodes.back().upper[a]=I::sum(origin,I::product(width,root.upper));
                }
                if (vertex >= 0)
                    vertexNodes[vertex] = id;
            }
            edge.nodes.push_back(id);
            if (edge.cells.size() == 1)
            {
                auto &node = nodes[id];
                bool declared = false;
                for (std::size_t k = 0; k < contacts.size(); ++k)
                {
                    if (!onEdge(contacts[k], edge))
                        continue;
                    const bool x = edge.begin.x != edge.end.x;
                    const auto parameter = M::ratio(
                        M::difference(x ? contacts[k].x : contacts[k].y, x ? edge.begin.x : edge.begin.y),
                        M::difference(x ? edge.end.x : edge.end.y, x ? edge.begin.x : edge.begin.y));
                    if (M::compare(root.lower, parameter) > 0 || M::compare(parameter, root.upper) > 0)
                        continue;
                    const auto value = source.restrictPhysical(edge.cells.front(), contacts[k], contacts[k]);
                    if (ExactExpansion::sign(value.numerator[0]) != 0)
                        continue;
                    require(matched[k] < 0 || matched[k] == id,
                            "declaration matches distinct retained roots");
                    matched[k] = id;
                    declared = true;
                }
                const bool axis = axisContact && M::compare(node.point[0], M::constant(0)) == 0 &&
                                  M::compare(node.point[1], M::constant(mesh.minCoord.y)) > 0 &&
                                  M::compare(node.point[1], M::constant(mesh.maxCoord.y)) < 0;
                int symmetry=-1;
                for (std::size_t k=0;k<symmetryPlanes.size();++k) {
                    const auto& plane=symmetryPlanes[k];const int normal=plane.normalAxis,tangent=1-normal;
                    if (M::compare(node.point[normal],M::constant(plane.coordinate))!=0) continue;
                    require(M::compare(node.point[tangent],M::constant(tangent==0?mesh.minCoord.x:mesh.minCoord.y))>0 &&
                        M::compare(node.point[tangent],M::constant(tangent==0?mesh.maxCoord.x:mesh.maxCoord.y))<0,
                        "symmetry endpoint at a domain corner is ambiguous");
                    require(symmetry<0 && !axis && !declared,"ambiguous symmetry endpoint identity");
                    symmetry=static_cast<int>(k);
                }
                require(axis || declared || symmetry>=0, "interface reaches an undeclared active boundary");
                const auto kind = axis ? Kind::Axis : symmetry>=0?Kind::SymmetryPlane:Kind::DeclaredBoundary;
                require(node.kind == Kind::None || node.kind == kind, "inconsistent boundary root identity");
                require(node.symmetryPlane<0 || node.symmetryPlane==symmetry,"inconsistent symmetry plane identity");
                node.kind = kind;
                node.symmetryPlane=symmetry;
            }
            r = next;
        }
        require(!legacy || !contacts.empty() || edge.nodes.size() <= 1,
                "default interface path does not allow multiple roots on one actual subedge");
    }
    if(!nativeTrace)for (int id : matched)
        require(id >= 0, "declaration is not an exact retained boundary zero");
    std::vector<std::array<std::vector<int>, 4>> sideNodes(source.cellCount());
    std::vector<std::array<int, 2>> connections;
    std::vector<Graph> graphs;
    std::vector<std::vector<int>> inwardContacts(source.cellCount());
    for (int c = 0; c < mesh.cellIndexNum(); ++c)
    {
        if (!source.cellActive(c))
            continue;
        if(sidePhase[c]!=0){cellPhase[c]=sidePhase[c];signs[c]=sidePhase[c];continue;}
        if constexpr(RequiresPieceCharts<Scalar>::value)
            if(signs[c]!=0){cellPhase[c]=signs[c];continue;}
        std::vector<int> local;
        for (int side = 0; side < 4; ++side)
        {
            auto &inventory = sideNodes[c][side];
            for (int edge : sides[c][side])
                inventory.insert(inventory.end(), edges[edge].nodes.begin(), edges[edge].nodes.end());
            const int axis = (side == 0 || side == 2) ? 1 : 0;
            std::sort(inventory.begin(), inventory.end(), [&](int a, int b) {
                return M::compare(nodes[a].point[axis], nodes[b].point[axis]) < 0;
            });
            inventory.erase(std::unique(inventory.begin(), inventory.end()), inventory.end());
            require(!legacy || !contacts.empty() || inventory.size() <= 1,
                    "default interface path does not allow multiple roots on one complete leaf side");
            for (int node : inventory)
                if (std::find(local.begin(), local.end(), node) == local.end())
                    local.push_back(node);
        }
        if (!chartBoxes[c].empty()) {
            appendChartGraphs(source,c,chartBoxes[c],sideNodes[c],nodes,graphs,connections);
            if constexpr(Scalar::dofCount==25) {
                // Subcharts already certified the incidences and matched the
                // actual boundary roots. Do not demand another full-coarse-
                // side isolation: its one floating parameter can be less
                // resolving than the retained actual subedge/chart ratios.
                for(int id:local)if(std::count_if(nodes[id].cells.begin(),nodes[id].cells.end(),
                    [&](int g){return graphs[g].nativeCellIndex==c;})==2)inwardContacts[c].push_back(id);
            }
            continue;
        }
        if constexpr(Scalar::dofCount==25) {
            if(!legacy) {
                Graph whole;whole.nativeCellIndex=c;
                applyNativeIncidences(source,whole,sideNodes[c],local,&inwardContacts[c]);
            }
        }
        const auto corners = source.cellVertices(c);
        const auto cornerOf = [&](int node) {
            if (nodes[node].vertex < 0)
                return -1;
            const auto found = std::find(corners.begin(), corners.end(), nodes[node].vertex);
            return found == corners.end() ? -1 : static_cast<int>(found - corners.begin());
        };
        unsigned excluded = 0;
        if(legacy || Scalar::dofCount!=25)local.erase(std::remove_if(local.begin(), local.end(),
                                   [&](int node) {
                                       const int corner = cornerOf(node);
                                       // A hanging vertex is an interior point of the COARSE side, not
                                       // a coarse corner. Keep that incidence; the fine corner rules
                                       // below must not erase the coarse graph's endpoint.
                                       if (corner < 0)
                                           return false;
                                       const auto d = source.zeroCornerGradientSigns(c, corner);
                                       const int xIn = (corner == 0 || corner == 3 ? 1 : -1) * d[0];
                                       const int yIn = (corner == 0 || corner == 1 ? 1 : -1) * d[1];
                                       if (xIn == 0 || xIn != yIn)
                                           return false;
                                       excluded |= 1u << corner;
                                       return true;
                                   }),
                    local.end());
        if (local.size() < 2)
        {
            if (!local.empty())
            {
                const int corner = cornerOf(local.front());
                require(corner >= 0, "isolated noncorner contact is unsupported");
                excluded |= 1u << corner;
            }
            const int phase = source.cellSignExceptVertices(c, excluded);
            if (phase == 0)
            {
                // Bernstein control signs are sufficient, not necessary.
                // A resolved Cartesian partial rules out an interior loop;
                // one common phase on ALL boundary intervals then certifies
                // the entire interior. This uses the retained edge inventory,
                // not corner samples or a relaxed residual threshold.
                int boundarySign = 0;
                bool certified = !legacy && local.empty() && resolvedPartial(source,c,legacy);
                for (const auto& side : sides[c]) for (int e : side)
                {
                    require(edges[e].query != nullptr, "candidate leaf lost its edge phase inventory");
                    for (const auto& phase : edges[e].query->phases)
                    {
                        const int sign = phase.negative ? -1 : 1;
                        certified = certified && (boundarySign == 0 || boundarySign == sign);
                        boundarySign = sign;
                    }
                }
                if (!certified || boundarySign == 0)
                    throw InterfaceGeometryResolution(InterfaceGeometryResolution::Reason::CriticalCell, c,
                        "no-graph leaf lacks a complete phase certificate");
                signs[c] = boundarySign;
            }
            cellPhase[c] = phase != 0 ? phase : signs[c];
            continue;
        }
        if (signs[c] != 0 || local.size() % 2 != 0 || (legacy && local.size() != 2))
            throw InterfaceGeometryResolution(InterfaceGeometryResolution::Reason::RootIncidence, c,
                "unpaired or inconsistent cell crossings");
        int multipleCoordinate = -1;
        if (local.size() > 2)
        {
            // A strict transverse partial makes the zero set single-valued.
            // Sorted boundary crossings then delimit disjoint graph intervals.
            // Multiple transverse branches are deliberately not guessed.
            if (resolvedPartial(source,c,legacy,1)) multipleCoordinate = 0;
            else if (resolvedPartial(source,c,legacy,0)) multipleCoordinate = 1;
            else throw InterfaceGeometryResolution(InterfaceGeometryResolution::Reason::GraphPairing, c,
                "several arcs need a certified common transverse partial");
            std::sort(local.begin(), local.end(), [&](int a, int b) {
                return M::compare(nodes[a].point[multipleCoordinate], nodes[b].point[multipleCoordinate]) < 0;
            });
        }
        const auto endpointIdentity = [&](int node) {
            for (int side = 0; side < 4; ++side)
            {
                const auto &inventory = sideNodes[c][side];
                const auto found = std::find(inventory.begin(), inventory.end(), node);
                if (found != inventory.end())
                    return LevelSetMappedInterface::Endpoint{side,
                                                             static_cast<int>(found - inventory.begin())};
            }
            throw std::runtime_error("[LevelSetNativeInterfaceChain] missing complete edge-root identity");
        };
        for (std::size_t pair = 0; pair < local.size(); pair += 2)
        {
            std::array<int, 2> ends{{local[pair], local[pair+1]}};
            Graph graph;
            graph.nativeCellIndex = c;
            // All charts of a complete native reconstruction use the same
            // controlled moments, including unsplit cubic neighbour cells.
            graph.controlledIntegration=!legacy;
            const auto bounds = source.leafBounds(c);
            I::Ratio span[2];
            for (int axis = 0; axis < 2; ++axis)
            {
                span[axis] = I::quotient(I::difference(nodes[ends[1]].point[axis], nodes[ends[0]].point[axis]),
                    M::ratio(M::difference(axis == 0 ? bounds[1].x : bounds[1].y,
                                          axis == 0 ? bounds[0].x : bounds[0].y), M::number(1)));
                if (ExactExpansion::sign(span[axis].numerator) < 0)
                    span[axis] = I::difference(M::constant(0), span[axis]);
            }
            graph.coordinateIsX = multipleCoordinate < 0 ? M::compare(span[0], span[1]) >= 0 : multipleCoordinate == 0;
            const bool coordinateConflict=!legacy && multipleCoordinate<0
                && !resolvedPartial(source,c,false,graph.coordinateIsX?1:0)
                && resolvedPartial(source,c,false,graph.coordinateIsX?0:1);
            if(coordinateConflict)
                graph.coordinateIsX=!graph.coordinateIsX;
            const int coordinate = graph.coordinateIsX ? 0 : 1;
            if (M::compare(span[coordinate], M::constant(0)) == 0)
                throw InterfaceGeometryResolution(InterfaceGeometryResolution::Reason::GraphPairing, c,
                    "zero graph extent or touching boundary roots");
            graph.begin = endpointIdentity(ends[0]);
            graph.end = endpointIdentity(ends[1]);
            if(!legacy && multipleCoordinate<0) {
                // End tangents are only a request to improve integration
                // coordinates, never a topology/regularity certificate.
                // Every resulting chart still needs the existing outward
                // polynomial proof and exact outer edge-root identities.
                bool splitCoordinate=coordinateConflict;
                if(!splitCoordinate) {
                    using Integrator=I::NativeIntegratorFor<Scalar>;
                    const auto sampling=Integrator::prepareGraphSampling(source,graph);
                    for(double endpoint:{0.,1.}) {
                        const auto sample=Integrator::graphSample(source,sampling,0,1,endpoint);
                        const auto jet=source.evaluateCellChart(c,sample.chart);
                        const double transverse=graph.coordinateIsX?jet.dy:jet.dx;
                        const double along=graph.coordinateIsX?jet.dx:jet.dy;
                        splitCoordinate|=std::abs(transverse)<.5*std::abs(along);
                    }
                }
                if(splitCoordinate) {
                    const auto partition=source.regularChartPartition(c);
                    if(!partition.empty()) {
                        appendChartGraphs(source,c,partition,sideNodes[c],nodes,graphs,connections);
                        break; // multipleCoordinate<0 implies this sole pair
                    }
                }
            }
            const int direction = I::NativeIntegratorFor<Scalar>::graphOrientation(source, graph);
            if (direction < 0) { std::swap(graph.begin, graph.end); std::swap(ends[0], ends[1]); }
            const int segment = static_cast<int>(graphs.size());
            graphs.push_back(graph);
            connections.push_back(ends);
            for (int node : ends) nodes[node].cells.push_back(segment);
        }
    }
    std::vector<SideArc> attached(graphs.size());
    if(nativeTrace)for(std::size_t e=0;e<edges.size();++e){auto& edge=edges[e];if(!edge.zero)continue;
        const int owner=*std::min_element(edge.cells.begin(),edge.cells.end());int side=-1;
        for(int s=0;s<4;++s)if(std::find(sides[owner][s].begin(),sides[owner][s].end(),static_cast<int>(e))!=sides[owner][s].end())side=s;
        require(side>=0 && edge.nodes.size()==2,"native zero edge lost actual side/vertex support");
        Graph graph;graph.nativeCellIndex=owner;graph.coordinateIsX=edge.begin.y==edge.end.y;graph.controlledIntegration=true;
        SideArc support;support.owner=owner;support.side=side;support.horizontal=graph.coordinateIsX;
        support.edgeLabel=edge.edgeLabel;
        support.begin=edge.begin;support.end=edge.end;support.beginVertex=edge.vertices[0];support.endVertex=edge.vertices[1];
        support.minus=edge.minus;support.plus=edge.plus;
        std::array<int,2> ends{{edge.nodes[0],edge.nodes[1]}};
        const int direction=graph.coordinateIsX?-edge.normalSign:edge.normalSign;
        if(direction<0){std::swap(ends[0],ends[1]);std::swap(support.begin,support.end);std::swap(support.beginVertex,support.endVertex);}
        edge.graph=static_cast<int>(graphs.size());graphs.push_back(graph);attached.push_back(support);connections.push_back(ends);
        for(int n:ends)nodes[n].cells.push_back(edge.graph);
    }
    // Only actual open endpoints need a boundary declaration. Intermediate
    // vertices of a zero boundary curve keep their ordinary two-valence.
    if(nativeTrace)for(auto& node:nodes)if(node.kind==Kind::None && node.cells.size()==1) {
        bool boundary=false;for(const auto& edge:edges)if(edge.cells.size()==1 && I::orientation(edge.begin,edge.end,node.point)==0) {
            const int tangent=edge.begin.x==edge.end.x?1:0;
            const auto lo=M::constant(tangent?edge.begin.y:edge.begin.x),hi=M::constant(tangent?edge.end.y:edge.end.x);
            boundary=boundary || (M::compare(node.point[tangent],lo)>=0 && M::compare(node.point[tangent],hi)<=0);
        }
        require(boundary,"native zero chain has a nonboundary open endpoint");
        bool declared=false;for(std::size_t k=0;k<contacts.size();++k)if(I::same(node.point,I::point(contacts[k]))) {
            require(matched[k]<0 || &node==&nodes[matched[k]],"zero endpoint declaration matches distinct nodes");
            matched[k]=static_cast<int>(&node-nodes.data());declared=true;
        }
        const bool axis=axisContact && M::compare(node.point[0],M::constant(0))==0
            && M::compare(node.point[1],M::constant(mesh.minCoord.y))>0 && M::compare(node.point[1],M::constant(mesh.maxCoord.y))<0;
        int symmetry=-1;for(std::size_t k=0;k<symmetryPlanes.size();++k){const auto& plane=symmetryPlanes[k];const int n=plane.normalAxis,t=1-n;
            if(M::compare(node.point[n],M::constant(plane.coordinate))!=0)continue;
            require(M::compare(node.point[t],M::constant(t?mesh.minCoord.y:mesh.minCoord.x))>0
                && M::compare(node.point[t],M::constant(t?mesh.maxCoord.y:mesh.maxCoord.x))<0,"symmetry zero endpoint at a domain corner is ambiguous");
            require(symmetry<0 && !axis && !declared,"ambiguous zero endpoint boundary identity");symmetry=static_cast<int>(k);
        }
        require(axis || declared || symmetry>=0,"native zero interface reaches an undeclared active boundary");
        node.kind=axis?Kind::Axis:symmetry>=0?Kind::SymmetryPlane:Kind::DeclaredBoundary;node.symmetryPlane=symmetry;
    }
    for(int id:matched)require(id>=0,"declaration is not an exact retained boundary zero");
    std::vector<std::vector<int>> nodeCells;
    for (const auto& node : nodes)
    {
        nodeCells.push_back(node.cells);
        if (node.cells.size() != (node.kind == Kind::None ? 2u : 1u))
            throw InterfaceGeometryResolution(InterfaceGeometryResolution::Reason::RootIncidence,
                node.cells.empty() ? -1 : graphs[node.cells.front()].nativeCellIndex,
                "root is branched, isolated, or missing its neighboring incidence");
    }
    Discovery discovery;
    discovery.cellPhase = std::move(cellPhase);
    auto& results = discovery.curves;
    std::vector<int> graphRegion(graphs.size(), -1);
    if (!legacy)
    {
        // Boundaries of NEGATIVE cell fragments: interface arcs and negative
        // actual-edge intervals, all oriented with the cell fragment on the left.
        // Join only across positive-length intervals, never a touching vertex.
        struct Arc {
            int begin, end, graph = -1, interval = -1;
            bool boundary = false;
            long double boundaryVolume = 0;
        };
        std::vector<std::vector<Arc>> arcs(source.cellCount());
        for (std::size_t g = 0; g < graphs.size(); ++g)if(attached[g].owner<0)
            arcs[graphs[g].nativeCellIndex].push_back({connections[g][0], connections[g][1], static_cast<int>(g)});
        for(const auto& edge:edges)if(edge.zero){int liquid=-1;
            for(const auto* s:{&edge.minus,&edge.plus})if(s->available() && s->phase<0)liquid=s->nativeCellIndex;
            if(liquid>=0)arcs[liquid].push_back({connections[edge.graph][0],connections[edge.graph][1],edge.graph,-1,
                edge.cells.size()==1 && !(axisContact && edge.begin.x==0),0});
        }
        const auto vertexNode = [&](int vertex) {
            return vertexNodes[vertex] >= 0 ? vertexNodes[vertex] : static_cast<int>(nodes.size()) + vertex;
        };
        int intervalCount = 0;
        for (std::size_t e = 0; e < edges.size(); ++e)
        {
            const auto& edge = edges[e];
            if(edge.zero)continue; // Interface boundary, not a shared negative interval.
            std::vector<int> points{vertexNode(edge.vertices[0])};
            std::vector<M::Ratio> parameters{M::constant(0)};
            for (int n : edge.nodes)
            {
                const int axis = edge.begin.x != edge.end.x ? 0 : 1;
                const auto parameter = I::quotient(I::difference(nodes[n].point[axis],
                    M::constant(axis == 0 ? edge.begin.x : edge.begin.y)),
                    M::ratio(M::difference(axis == 0 ? edge.end.x : edge.end.y,
                                          axis == 0 ? edge.begin.x : edge.begin.y), M::number(1)));
                if (M::compare(parameter, M::constant(0)) > 0 && M::compare(parameter, M::constant(1)) < 0)
                { points.push_back(n); parameters.push_back(parameter); }
            }
            points.push_back(vertexNode(edge.vertices[1])); parameters.push_back(M::constant(1));
            for (std::size_t k = 0; k + 1 < points.size(); ++k)
            {
                int sign = 0;
                if (!edge.query) sign = signs[edge.cells.front()];
                else for (const auto& phase : edge.query->phases)
                {
                    if (M::compare(phase.end, parameters[k]) <= 0 ||
                        M::compare(phase.begin, parameters[k+1]) >= 0) continue;
                    const int partSign = phase.negative ? -1 : 1;
                    require(sign == 0 || sign == partSign, "edge phase changes without a retained root");
                    sign = partSign;
                }
                require(sign != 0, "missing edge phase interval");
                if (sign > 0) continue;
                const int interval = intervalCount++;
                const bool boundary = edge.cells.size() == 1 &&
                    !(axisContact && edge.begin.x == 0 && edge.end.x == 0);
                for (int c : edge.cells)
                {
                    int side = -1;
                    for (int s = 0; s < 4; ++s)
                        if (std::find(sides[c][s].begin(), sides[c][s].end(), static_cast<int>(e)) != sides[c][s].end()) side = s;
                    require(side >= 0, "negative edge lost native side ownership");
                    const bool forward = side == 1 || side == 2; // W/S/E/N: CW boundary normals, CCW curve.
                    long double moment = 0;
                    if (boundary && edge.begin.x == edge.end.x)
                    {
                        const long double r = edge.begin.x;
                        const long double length = M::approximate(I::difference(parameters[k+1], parameters[k])) *
                            (static_cast<long double>(edge.end.y) - edge.begin.y);
                        moment = (side == 2 ? 1 : -1) * length *
                            (axisymmetric ? std::acos(-1.0L)*r*r : r);
                    }
                    arcs[c].push_back({points[k + (forward ? 0 : 1)], points[k + (forward ? 1 : 0)],
                        -1, interval, boundary, moment});
                }
            }
        }
        struct Fragment { int cell; bool boundary; long double boundaryVolume = 0; };
        std::vector<Fragment> fragments;
        std::vector<int> parent, intervalOwner(intervalCount, -1), graphFragment(graphs.size(), -1);
        const auto root = [&](int n) { while (parent[n] != n) n = parent[n]; return n; };
        for (int c = 0; c < mesh.cellIndexNum(); ++c)
        {
            const auto& local = arcs[c];
            std::map<int, std::vector<int>> outgoing;
            for (std::size_t a = 0; a < local.size(); ++a) {
                auto& next=outgoing[local[a].begin];next.push_back(static_cast<int>(a));
                if(next.size()>1 && (next.size()!=2 || std::find(inwardContacts[c].begin(),inwardContacts[c].end(),
                    local[a].begin)==inwardContacts[c].end()))
                    throw InterfaceGeometryResolution(InterfaceGeometryResolution::Reason::RootIncidence, c,
                        "negative fragment boundary touches or branches");
            }
            std::vector<bool> visited(local.size(), false);
            for (std::size_t seed = 0; seed < local.size(); ++seed)
            {
                if (visited[seed]) continue;
                const int f = static_cast<int>(fragments.size());
                parent.push_back(f); fragments.push_back({c, false});
                int a = static_cast<int>(seed);
                for (;;)
                {
                    require(!visited[a], "negative fragment failed to close"); visited[a] = true;
                    const auto& arc = local[a];
                    fragments[f].boundary = fragments[f].boundary || arc.boundary;
                    fragments[f].boundaryVolume += arc.boundaryVolume;
                    if (arc.graph >= 0) graphFragment[arc.graph] = f;
                    if (arc.interval >= 0)
                    {
                        auto& owner = intervalOwner[arc.interval];
                        if (owner >= 0) parent[root(f)] = root(owner); else owner = f;
                    }
                    const auto next = outgoing.find(arc.end);
                    require(next != outgoing.end(), "negative fragment has an open cell boundary");
                    int successor=-1;
                    for(int candidate:next->second) {
                        // At an inward tangency with liquid on the boundary,
                        // each fragment turns between the interface and the
                        // cell edge. Joining graph to graph would connect two
                        // fragments through a zero-measure point.
                        if(next->second.size()==2 && (arc.graph>=0)==(local[candidate].graph>=0))continue;
                        require(successor<0,"contact has ambiguous negative-fragment continuation");
                        successor=candidate;
                    }
                    require(successor>=0,"contact lost its negative-fragment continuation");
                    if(successor==static_cast<int>(seed))break;
                    a=successor;
                }
            }
        }
        std::map<int, int> regionIds;
        std::vector<int> fragmentRegion(fragments.size(), -1);
        for (std::size_t f = 0; f < fragments.size(); ++f)
        {
            const int representative = root(static_cast<int>(f));
            auto found = regionIds.find(representative);
            if (found == regionIds.end())
            {
                const int id = static_cast<int>(discovery.negativeRegions.size());
                found = regionIds.emplace(representative, id).first;
                discovery.negativeRegions.emplace_back();
            }
            const int id = found->second; fragmentRegion[f] = id;
            auto& region = discovery.negativeRegions[id];
            region.cells.push_back(fragments[f].cell);
            if (signs[fragments[f].cell] < 0) region.pureCells.push_back(fragments[f].cell);
            region.touchesNonAxisBoundary = region.touchesNonAxisBoundary || fragments[f].boundary;
            region.fixedBoundaryVolume += fragments[f].boundaryVolume;
        }
        for (auto& region : discovery.negativeRegions)
        {
            std::sort(region.cells.begin(), region.cells.end());
            region.cells.erase(std::unique(region.cells.begin(), region.cells.end()), region.cells.end());
            std::sort(region.pureCells.begin(), region.pureCells.end());
            region.pureCells.erase(std::unique(region.pureCells.begin(), region.pureCells.end()), region.pureCells.end());
        }
        for (std::size_t g = 0; g < graphs.size(); ++g)
        {
            if(graphFragment[g]<0 && attached[g].owner>=0
                && !(attached[g].minus.available() && attached[g].minus.phase<0)
                && !(attached[g].plus.available() && attached[g].plus.phase<0))continue;
            require(graphFragment[g] >= 0, "interface arc lost its negative region");
            graphRegion[g] = fragmentRegion[graphFragment[g]];
        }
    }

    std::vector<bool> used(graphs.size(), false);
    // Begin open chains at their oriented boundary root; then visit closed ones.
    std::vector<int> seeds;
    for (std::size_t g = 0; g < graphs.size(); ++g)
        if (nodes[connections[g][0]].kind != Kind::None) seeds.push_back(static_cast<int>(g));
    for (std::size_t g = 0; g < graphs.size(); ++g) seeds.push_back(static_cast<int>(g));
    for (int seed : seeds)
    {
        if (used[seed]) continue;
        Result result;
        result.negativeRegion = graphRegion[seed];
        result.mapping.assign(source.cellCount(), -1);
        result.cellGraphs.resize(source.cellCount());
        result.minimumSpacing = minimumSpacing;
        const int start = connections[seed][0];
        result.begin = nodes[start].kind;
        result.endpointSymmetryPlanes[0]=nodes[start].symmetryPlane;
        result.closed = result.begin == Kind::None;
        int segment = seed;
        for (;;)
        {
            require(!used[segment], "component traversal revisits an arc before closure");
            used[segment] = true;
            require(graphRegion[segment] == result.negativeRegion, "one surface borders inconsistent negative regions");
            const int c = graphs[segment].nativeCellIndex;
            auto& cellGraphs = result.cellGraphs[c];
            cellGraphs.push_back(static_cast<int>(result.graphs.size()));
            result.mapping[c] = cellGraphs.size() == 1 ? cellGraphs.front() : -2;
            result.graphs.push_back(graphs[segment]);
            result.sideSupports.push_back(attached[segment]);
            result.rootNodes.push_back(connections[segment]);
            const int end = connections[segment][1];
            if (end == start)
            {
                require(result.closed, "open component returned to its starting boundary");
                result.end = Kind::None;
                break;
            }
            if (nodes[end].kind != Kind::None)
            {
                require(!result.closed, "closed component encountered an unmatched endpoint");
                result.end = nodes[end].kind;
                result.endpointSymmetryPlanes[1]=nodes[end].symmetryPlane;
                break;
            }
            const auto& neighbors = nodeCells[end];
            const int next = neighbors[0] == segment ? neighbors[1] : neighbors[0];
            require(connections[next][0] == end, "neighboring arcs disagree on interface orientation");
            segment = next;
        }
        if (result.graphs.size() < 3)
            throw InterfaceGeometryResolution(InterfaceGeometryResolution::Reason::SmallComponent,
                graphs[seed].nativeCellIndex, "component has fewer than three retained graphs");
        for (std::size_t k = 0; k < contacts.size(); ++k)
            if (matched[k] == start || matched[k] == result.rootNodes.back()[1])
                result.boundaryContacts.push_back(contacts[k]);
        results.push_back(std::move(result));
    }
    if (legacy) require(results.size() == 1, "legacy consumer requires exactly one interface component");
    return discovery;
}
} // namespace NativeInterfaceChain
} // namespace
} // namespace LAMP::Geometry::Detail
