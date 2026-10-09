#pragma once

#include "Core/Geometry/ActiveQuarticHermitePatch.h"
#include "Core/Geometry/ActiveHermiteSensitivity.h"
#include "Core/Geometry/LevelSetHermiteGeometry.h"
#include "Core/Geometry/InterfaceGeometryResolution.h"
#include "Core/Geometry/PolynomialVelocitySlab.h"
#include "Core/Geometry/Detail/LevelSetHermiteQuarticCell.h"
#include "Core/Geometry/Detail/LevelSetNativeCornerSector.h"
#include "Core/Mesh/CartesianLeafLocator.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstring>
#include <exception>
#include <iterator>
#include <limits>
#include <map>
#include <memory>
#include <mutex>
#include <tuple>
#include <utility>
#include <vector>
#include "Core/Geometry/Detail/NativeHermiteTransferSampling.h"

namespace Quantum::Geometry::Detail
{
namespace
{
namespace LevelSetHermiteScalar
{
// Native data/constraint owner for the SAME QuarticCell tensor kernel.
// Maximal shared sides own two midpoint data. Hanging jets and fine-side
// midpoints are exact restrictions of their coarse value/normal traces.
// A retained linear dependency graph supplies the matching transpose; there
// is no separately fitted hanging scalar or hand-coded midpoint Jacobian.
// No geometry projection, root isolation, motion, boundary law or IO here.
class QuarticLeafData
{
    using Exact=ExactExpansion::Expansion;
    using Patch=ActiveQuarticHermitePatch;
    using Sensitivity=ActiveQuarticHermiteSensitivity;
    using EdgeKey=std::array<int,2>;
    template<std::size_t N> struct Homogeneous
    {
        std::array<Exact,N> numerator{};
        Exact denominator=ExactExpansion::fromDouble(1.0);
    };
    struct Dependency { int source=-1; long double weight=0; };
    struct Datum
    {
        Homogeneous<1> value;
        std::vector<Dependency> dependencies;
        Sensitivity identity;
        bool fixed=false;
    };
    struct VertexData
    {
        Vector2D point;
        int coarseCell=-1, coarseEdge=-1, state=0;
        std::array<int,4> jet{{-1,-1,-1,-1}};
        bool used=false;
    };
    struct EdgeData
    {
        EdgeKey vertices;
        int parent=-1, input=-1, state=0;
        std::array<int,2> midpoint{{-1,-1}};
        bool horizontal=false;
    };
    struct BandClassification
    {
        double halfWidth=0;
        bool inside=false,outside=false,ready=false;
    };
    struct Leaf
    {
        Vector2D lower,upper;
        std::array<int,4> vertices{},edges{}; // SW SE NW NE; south north west east
        std::array<int,25> physical{};
        Homogeneous<25> local;
        std::unique_ptr<QuarticCell> cache;
        // The mathematical data never change after construction. Cache only
        // their exact Bernstein signs, owned by this leaf in this snapshot.
        // call_once preserves the original on-query arithmetic failure and
        // safely publishes the same result to concurrent face queries.
        mutable std::once_flag phaseSignsReady;
        mutable std::array<int,25> phaseSigns{};
        mutable std::once_flag resolvedPhaseReady;
        mutable int resolvedPhase=0;
        mutable std::once_flag zeroPartialsReady;
        mutable std::array<bool,2> zeroPartials{{false,false}};
        // Immutable outward controls reused by all physical/local subbox
        // queries. Environment checks remain per call; failed initialization
        // remains retryable, as with the other exact-data caches.
        mutable std::once_flag boxControlsReady;
        mutable std::array<ExactExpansion::Interval,25> boxControls{};
        mutable int boxControlSign=0;
        // A bounded cache for repeated corner/midpoint/center point queries.
        // Values come from the ORIGINAL exact physical restriction, not from
        // approximate nodal jets or assumed cardinal interpolation. Each new
        // immutable scalar owns fresh slots; failures remain retryable.
        mutable std::array<std::once_flag,9> pointValuesReady;
        mutable std::array<double,9> pointValues{};
        // Whole-cell band certificates depend on this immutable scalar and
        // the exact requested width. Retain two widths, not an unbounded map.
        // A mutex protects replacement as well as same-leaf cold queries.
        mutable std::mutex bandMutex;
        mutable std::array<BandClassification,2> bands{};
        mutable std::size_t nextBand=0;
        bool active=false;
        int sectorIndex=-1,sectorCorner=-1,sectorPhaseSign=0;
        std::array<int,4> sectorData{{-1,-1,-1,-1}};
    };

  public:
    static constexpr std::size_t dofCount=25;
    using Polynomial=Homogeneous<9>;
    enum class Build { Geometry, Skeleton };
    QuarticLeafData(const QuadtreeMesh& mesh,const Patch& patch,bool radialAxisSymmetry=false,
        const std::vector<int>* queryCells=nullptr,Build build=Build::Geometry,const QuarticLeafData* restrictionSource=nullptr)
        : locator_(mesh), vertices_(mesh.vertIndexNum()), cells_(mesh.cellIndexNum()),
          lower_(mesh.minCoord),upper_(mesh.maxCoord),radial_(radialAxisSymmetry),
          restrictionSource_(restrictionSource),restrictionPatch_(restrictionSource?&patch:nullptr)
    {
        ExactExpansion::validateEnvironment();
        Detail::require(!restrictionSource || (build==Build::Geometry && !queryCells && patch.cornerSectors.empty()),
            "exact quartic refinement requires complete geometry without sectors");
        Detail::require(patch.nodes.size()==vertices_.size() && patch.activeCells.size()==cells_.size()
            && (build==Build::Skeleton || patch.cellCenters.size()==cells_.size()),
            "native quartic patch requires full vertex/cell arrays");
        Detail::require(build!=Build::Skeleton || (patch.cellCenters.empty() && patch.cornerSectors.empty() && !queryCells),
            "quartic skeleton input must not contain cell centers, corner sectors or query-cell geometry");
        Detail::require(!radial_ || lower_.x==0.0,"native quartic radial symmetry requires r=0");
        std::vector<std::uint8_t> selected;
        if (queryCells) {
            Detail::require(!queryCells->empty() && patch.cornerSectors.empty(),
                "quartic query subset requires cells and no corner sectors");
            selected.assign(cells_.size(),0);
            std::vector<int> queue=*queryCells;
            for (std::size_t next=0;next<queue.size();++next) {
                const int c=queue[next];
                Detail::require(c>=0 && static_cast<std::size_t>(c)<cells_.size() && patch.activeCells[c],
                    "quartic query subset selected an inactive or invalid leaf");
                if (selected[c]) continue;
                selected[c]=1;
                // A hanging vertex still belongs to its actual coarse trace.
                // Include that source and recursively retain its own support.
                for (int label:mesh.cells(mesh.cellIndexToLabel[c]).cellToVert) {
                    const int coarse=mesh.verts(label).supCell;
                    if (coarse>=0) queue.push_back(mesh.cells(coarse).index());
                }
            }
        }
        prepareTopology(mesh,queryCells?selected:patch.activeCells);
        for(std::size_t i=0;i<patch.edgeMidpoints.size();++i)
        {
            const auto found=edgeIndex_.find(patch.edgeMidpoints[i].vertices);
            if (queryCells && found==edgeIndex_.end()) continue;
            Detail::require(found!=edgeIndex_.end(),"quartic midpoint is not an ordered active leaf side");
            auto& edge=edges_[found->second];
            Detail::require(edge.input<0,"duplicate quartic midpoint side identity");
            edge.input=static_cast<int>(i);
        }
        for(std::size_t i=0;i<vertices_.size();++i)
            if(vertices_[i].used) resolveVertex(static_cast<int>(i),patch);
        for(std::size_t i=0;i<edges_.size();++i) resolveMidpoint(static_cast<int>(i),patch);
        if(build==Build::Skeleton) return; // immutable constraint graph only
        prepareSectors(mesh,patch);
        for(std::size_t i=0;i<cells_.size();++i)
        {
            auto& cell=cells_[i];
            if(!cell.active) continue;
            for(int y=0;y<2;++y) for(int x=0;x<2;++x)
            {
                const auto& jet=cell.sectorIndex>=0 && cell.sectorCorner==2*y+x
                    ? cell.sectorData:vertices_[cell.vertices[2*y+x]].jet;
                cell.physical[5*y+x]=jet[0];
                cell.physical[5*y+x+2]=jet[1];
                cell.physical[5*(y+2)+x]=jet[2];
                cell.physical[5*(y+2)+x+2]=jet[3];
            }
            for(int side=0;side<2;++side)
            {
                const auto& h=edges_[cell.edges[side]].midpoint;
                const auto& v=edges_[cell.edges[side+2]].midpoint;
                cell.physical[5*side+4]=h[0]; cell.physical[5*(side+2)+4]=h[1];
                cell.physical[20+side]=v[0]; cell.physical[22+side]=v[1];
            }
            cell.physical[24]=independent(patch.cellCenters[i],Sensitivity::Kind::CellCenter,
                static_cast<int>(i),0);
        }
        // Freeze the shared dependency graph and all independent IDs first.
        // Each remaining transform owns one leaf; common() only reads data_.
        // No parallel registration, change of constraint order or reduction.
        const auto prepareLeaf=[&](std::size_t i) {
            ExactExpansion::validateEnvironment();
            auto& cell=cells_[i];
            if (!cell.active) return;
            if(restrictionSource_)cell.local=restrictionSource_->restrictionData(cell.lower,cell.upper);
            else cell.local=common(cell.physical);
            const auto hx=difference(cell.upper.x,cell.lower.x),hy=difference(cell.upper.y,cell.lower.y);
            QuarticCell::Dofs approximate{};
            for(int y=0;y<5;++y) for(int x=0;x<5;++x)
            {
                auto& value=cell.local.numerator[5*y+x];
                if(!restrictionSource_ && (x==2 || x==3)) value=ExactExpansion::multiply(value,hx);
                if(!restrictionSource_ && (y==2 || y==3)) value=ExactExpansion::multiply(value,hy);
            }
            normalize(cell.local);
            for(std::size_t k=0;k<25;++k)
                approximate[k]=ratio(cell.local.numerator[k],cell.local.denominator);
            cell.cache=std::make_unique<QuarticCell>(approximate,cell.upper.x-cell.lower.x,cell.upper.y-cell.lower.y);
            if(cell.sectorIndex>=0)
            {
                constexpr unsigned masks[4]={1u,2u,8u,4u};
                Detail::require(cellSignExceptVertices(static_cast<int>(i),masks[cell.sectorCorner])==cell.sectorPhaseSign,
                    "native quartic corner sector lacks its exact single-phase certificate");
            }
        };
        if (!cells_.empty()) prepareLeaf(0);
        std::vector<std::exception_ptr> errors(cells_.size());
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 32) if(cells_.size()>64)
#endif
        for (std::ptrdiff_t i=1;i<static_cast<std::ptrdiff_t>(cells_.size());++i)
        {
            try {prepareLeaf(static_cast<std::size_t>(i));}
            catch (...) {errors[i]=std::current_exception();}
        }
        // Construction has not published a snapshot. Preserve native-cell
        // failure order and unwind all completed private leaves on failure.
        for (const auto& error:errors) if (error) std::rethrow_exception(error);
        // Borrowed construction inputs are never retained by the immutable result.
        restrictionSource_=nullptr;restrictionPatch_=nullptr;
    }

    static std::unique_ptr<QuarticLeafData> refined(const QuadtreeMesh& mesh,const QuarticLeafData& source)
    {
        ExactExpansion::validateEnvironment();
        Detail::require(source.sectorCount_==0,"exact quartic refinement does not support corner sectors");
        Detail::require(mesh.minCoord.x==source.lower_.x && mesh.minCoord.y==source.lower_.y
            && mesh.maxCoord.x==source.upper_.x && mesh.maxCoord.y==source.upper_.y,"exact quartic refinement changes domain");
        for(const auto& cell:source.cells_)Detail::require(cell.active,"exact quartic refinement requires full active source");
        // This patch declares only the target topology. Its dummy floating
        // payload is replaced by retained exact source data at registration.
        Patch topology;topology.activeCells.assign(mesh.cellIndexNum(),1);topology.nodes.resize(mesh.vertIndexNum());
        topology.cellCenters.assign(mesh.cellIndexNum(),0);std::map<EdgeKey,bool> sides;
        for(int c=0;c<mesh.cellIndexNum();++c){const auto& leaf=mesh.cells(mesh.cellIndexToLabel[c]);
            const auto lo=mesh.verts(leaf.cellToVert[0]).coord(),hi=mesh.verts(leaf.cellToVert[2]).coord();int parent=-1;
            for(int p:source.leafCandidates(lo,hi)){const auto b=source.leafBounds(p);
                if(std::max(lo.x,b[0].x)>=std::min(hi.x,b[1].x) || std::max(lo.y,b[0].y)>=std::min(hi.y,b[1].y))continue;
                Detail::require(parent<0 && lo.x>=b[0].x && lo.y>=b[0].y && hi.x<=b[1].x && hi.y<=b[1].y,
                    "exact quartic refinement cannot coarsen or span source leaves");parent=p;}
            Detail::require(parent>=0,"exact quartic refinement lacks source coverage");
            for(int k=0;k<4;++k){EdgeKey key{{mesh.verts(leaf.cellToVert[k]).index(),mesh.verts(leaf.cellToVert[(k+1)%4]).index()}};
                const auto a=mesh.verts(mesh.vertIndexToLabel[key[0]]).coord(),b=mesh.verts(mesh.vertIndexToLabel[key[1]]).coord();
                if(a.y==b.y?a.x>b.x:a.y>b.y)std::swap(key[0],key[1]);sides.emplace(key,true);}}
        for(const auto& side:sides)topology.edgeMidpoints.push_back({side.first,0,0});
        std::unique_ptr<QuarticLeafData> result;
        try {result=std::make_unique<QuarticLeafData>(mesh,topology,source.radial_,nullptr,Build::Geometry,&source);}
        catch(const std::exception& error){throw std::runtime_error(std::string("exact refinement construction: ")+error.what());}
        const NativeHermiteTransfer::DatumSampler<QuarticLeafData> sampler(source);
        // Hermite unisolvence: compare ALL 25 normalized data exactly, not a
        // floating sampling test. Shared constraints are checked by this too.
        for(std::size_t c=0;c<result->cells_.size();++c){const auto& cell=result->cells_[c];
            for(int y=0;y<5;++y)for(int x=0;x<5;++x){try {const auto coordinate=[](int i,double lo,double hi){
                    return std::array<double,2>{{i==1 || i==3?hi:lo,i==0 || i==2?lo:hi}};};
                const auto px=coordinate(x,cell.lower.x,cell.upper.x),py=coordinate(y,cell.lower.y,cell.upper.y);
                const unsigned dx=x==2 || x==3?1:0,dy=y==2 || y==3?1:0;
                const auto expected=source.restrictionDatum(sampler.locate({px[0],py[0]},{px[1],py[1]}),dx,dy);
                const auto& physical=result->data_.at(cell.physical[5*y+x]).value;
                Detail::require(equalRatios(physical.numerator[0],physical.denominator,expected.numerator[0],expected.denominator),
                    "exact refinement shared constraint mismatch");
                auto numerator=expected.numerator[0];if(dx)numerator=ExactExpansion::multiply(numerator,difference(cell.upper.x,cell.lower.x));
                if(dy)numerator=ExactExpansion::multiply(numerator,difference(cell.upper.y,cell.lower.y));
                const auto message="exact quartic refinement polynomial mismatch cell="+std::to_string(c);
                Detail::require(equalRatios(numerator,expected.denominator,cell.local.numerator[5*y+x],cell.local.denominator),message.c_str());}
                catch(const std::exception& error){throw std::runtime_error("exact refinement certificate cell="+std::to_string(c)
                    +" datum="+std::to_string(5*y+x)+": "+error.what());}}}
        return result;
    }

    std::size_t cellCount() const noexcept { return cells_.size(); }
    Vector2D minCoord() const noexcept { return lower_; }
    Vector2D maxCoord() const noexcept { return upper_; }
    std::size_t vertexCount() const noexcept { return vertices_.size(); }
    bool radialAxisSymmetry() const noexcept { return radial_; }
    bool hasCornerSectors() const noexcept { return sectorCount_!=0; }
    std::vector<int> leafCandidates(Vector2D lo,Vector2D hi) const
    { return locator_.intersectingBox(lo,hi); }
    bool cellActive(int index) const
    {
        Detail::require(index>=0 && static_cast<std::size_t>(index)<cells_.size(),"quartic leaf index out of range");
        return cells_[index].active;
    }
    std::array<Vector2D,2> leafBounds(int index) const
    {
        cellActive(index);
        return {{cells_[index].lower,cells_[index].upper}};
    }
    std::array<int,4> cellVertices(int index) const
    {
        const auto& corners=activeCell(index).vertices;
        return {{corners[0],corners[1],corners[3],corners[2]}}; // SW SE NE NW
    }
    int selectedCell(Vector2D point) const
    {
        const auto incident=locator_.incident(point);
        int selected=-1;
        for(int k=0;k<incident.count;++k)
            if(cellActive(incident.indices[k]) && (selected<0 || incident.indices[k]<selected))
                selected=incident.indices[k];
        Detail::require(selected>=0,"native quartic query is outside the active domain");
        return selected;
    }
    CellChartPoint chartAtPoint(int index,Vector2D point) const
    {
        ExactExpansion::validateEnvironment();
        const auto& cell=activeCell(index);
        const auto local=[](double point,double low,double high)
        {
            Detail::require(std::isfinite(point) && point>=low && point<=high,
                "native quartic point leaves its actual closed leaf");
            const bool upper=point-low>high-point;
            const double offset=static_cast<double>(upper
                ? (static_cast<long double>(high)-point)/(static_cast<long double>(high)-low)
                : (static_cast<long double>(point)-low)/(static_cast<long double>(high)-low));
            return CellChartCoordinate{offset,upper};
        };
        return {local(point.x,cell.lower.x,cell.upper.x),local(point.y,cell.lower.y,cell.upper.y)};
    }
    Vector2D chartDisplayPoint(int index,CellChartPoint point) const
    {
        validateCellChart(point);
        const auto& cell=activeCell(index);
        const auto physical=[](CellChartCoordinate coordinate,double low,double high)
        {
            const auto tail=ExactExpansion::multiply(difference(high,low),number(coordinate.offset));
            return ExactExpansion::estimate(coordinate.fromUpper
                ? ExactExpansion::subtract(number(high),tail):ExactExpansion::add(number(low),tail));
        };
        return {physical(point.x,cell.lower.x,cell.upper.x),physical(point.y,cell.lower.y,cell.upper.y)};
    }
    Sample evaluateCell(int index,Vector2D point) const
    { return evaluateCellChart(index,chartAtPoint(index,point)); }
    double valueInCell(int index,Vector2D point) const
    {
        ExactExpansion::validateEnvironment();
        const auto& cell=activeCell(index);
        const auto evaluate=[&] {
            const auto polynomial=restrictPhysical(index,point,point);
            return ratio(polynomial.numerator[0],polynomial.denominator);
        };
        const auto slot=[](double coordinate,double low,double high) {
            const std::array<double,3> points{{low,low+(high-low)*.5,high}};
            for(int k=0;k<3;++k)
                if(std::memcmp(&coordinate,&points[k],sizeof(double))==0) return k;
            return -1;
        };
        const int x=slot(point.x,cell.lower.x,cell.upper.x);
        const int y=slot(point.y,cell.lower.y,cell.upper.y);
        if(x<0 || y<0) return evaluate();
        const int position=3*y+x;
        std::call_once(cell.pointValuesReady[position],[&] {cell.pointValues[position]=evaluate();});
        return cell.pointValues[position];
    }
    int hangingNodeCount() const noexcept { return hangingCount_; }
    int constrainedEdgeCount() const noexcept { return constrainedEdgeCount_; }
    std::size_t independentDofCount() const noexcept { return independentCount_; }
    std::array<double,4> nodeJet(int index) const
    {
        Detail::require(index>=0 && static_cast<std::size_t>(index)<vertices_.size() && vertices_[index].used,
            "quartic nodal query requires a used native vertex");
        std::array<double,4> result;
        for(int k=0;k<4;++k) result[k]=value(vertices_[index].jet[k]);
        return result;
    }
    std::vector<int> skeletonVertices() const
    {
        std::vector<int> result;
        for(std::size_t i=0;i<vertices_.size();++i)if(vertices_[i].used)result.push_back(static_cast<int>(i));
        return result;
    }
    std::vector<EdgeKey> skeletonEdges() const
    {
        std::vector<EdgeKey> result;
        for(const auto& e:edgeIndex_)result.push_back(e.first);
        return result;
    }
    EdgeKey vertexMasterSide(int index) const
    {
        Detail::require(index>=0 && static_cast<std::size_t>(index)<vertices_.size() && vertices_[index].used,
            "quartic vertex ownership requires a used native vertex");
        const int edge=vertices_[index].coarseEdge;
        return edge<0?EdgeKey{{-1,-1}}:edges_[edge].vertices;
    }
    EdgeKey edgeParentSide(const EdgeKey& key) const
    {
        const auto& edge=skeletonEdge(key);
        return edge.parent<0?EdgeKey{{-1,-1}}:edges_[edge.parent].vertices;
    }
    std::array<double,2> skeletonMidpoint(const EdgeKey& key) const
    {
        const auto& edge=skeletonEdge(key);
        return {{value(edge.midpoint[0]),value(edge.midpoint[1])}};
    }
    std::array<double,4> skeletonTrace(const EdgeKey& key,double coordinate) const
    {
        ExactExpansion::validateEnvironment();
        const auto& edge=skeletonEdge(key);
        const auto a=vertices_[key[0]].point,b=vertices_[key[1]].point;
        const auto t=chart(coordinate,coordinate,edge.horizontal?a.x:a.y,edge.horizontal?b.x:b.y);
        std::array<double,4> result;
        visitTrace(edge,t.numeratorBegin,t.denominator,[&](int component,const auto& inputs,const auto& weights,const auto& denominator) {
            auto combined=linearValue(inputs,weights,denominator,[](std::size_t){});
            normalize(combined);
            result[component]=ratio(combined.numerator[0],combined.denominator);
        });
        return result;
    }
    std::array<double,25> mappedCellPhysicalDofs(int index) const
    {
        const auto& cell=activeCell(index);
        std::array<double,25> result;
        for(std::size_t k=0;k<25;++k) result[k]=value(cell.physical[k]);
        return result;
    }
    QuarticCell::Dofs mappedCellLocalDofs(int index) const { return activeCell(index).cache->localDofs(); }
    // Certify the retained polynomial, including exact hanging constraints.
    // Coefficients are physical {value at origin, dx, dy}/denominator.
    // Every closed intersecting leaf must carry the same affine polynomial.
    Homogeneous<3> affineOnBox(Vector2D lower,Vector2D upper,Vector2D origin) const
    {
        ExactExpansion::validateEnvironment();
        Detail::require(std::isfinite(lower.x) && std::isfinite(lower.y) &&
            std::isfinite(upper.x) && std::isfinite(upper.y) &&
            lower.x>=lower_.x && lower.y>=lower_.y && upper.x<=upper_.x && upper.y<=upper_.y &&
            lower.x<=upper.x && lower.y<=upper.y,"affine swept box leaves the retained domain");
        const auto candidates=leafCandidates(lower,upper);
        Detail::require(!candidates.empty(),"affine swept box has no retained cells");
        Homogeneous<3> result; bool first=true;
        for(int index:candidates) {
            const auto& cell=activeCell(index);
            Detail::require(cell.sectorIndex<0,"affine swept box contains a corner sector");
            const auto& n=cell.local.numerator;
            const std::array<double,5> coordinate{{0,1,0,0,.5}};
            for(int y=0;y<5;++y)for(int x=0;x<5;++x) {
                const bool dx=x==2 || x==3,dy=y==2 || y==3;
                const auto expected=dx?(dy?number(0):n[2]):dy?n[10]:ExactExpansion::add(n[0],
                    ExactExpansion::add(ExactExpansion::multiply(n[2],number(coordinate[x])),
                        ExactExpansion::multiply(n[10],number(coordinate[y]))));
                Detail::require(ExactExpansion::sign(ExactExpansion::subtract(n[5*y+x],expected))==0,
                    "affine swept box contains a non-affine native Q4");
            }
            const auto hx=difference(cell.upper.x,cell.lower.x),hy=difference(cell.upper.y,cell.lower.y);
            const auto area=ExactExpansion::multiply(hx,hy);
            Homogeneous<3> current;
            current.denominator=ExactExpansion::multiply(cell.local.denominator,area);
            current.numerator[1]=ExactExpansion::multiply(n[2],hy);
            current.numerator[2]=ExactExpansion::multiply(n[10],hx);
            current.numerator[0]=ExactExpansion::add(ExactExpansion::multiply(n[0],area),
                ExactExpansion::add(ExactExpansion::multiply(current.numerator[1],difference(origin.x,cell.lower.x)),
                    ExactExpansion::multiply(current.numerator[2],difference(origin.y,cell.lower.y))));
            ExactExpansion::centerHomogeneous(current.numerator,current.denominator);
            if(first){result=current;first=false;}
            else for(int k=0;k<3;++k)
                Detail::require(ExactExpansion::sign(ExactExpansion::subtract(
                    ExactExpansion::multiply(current.numerator[k],result.denominator),
                    ExactExpansion::multiply(result.numerator[k],current.denominator)))==0,
                    "affine swept box contains different affine polynomials");
        }
        return result;
    }
    ::LevelSetHermiteGeometry::Snapshot::QuarticCenterInterval centerInterval(int index,
        double displacementCells,double retainedFraction,
        ::LevelSetHermiteGeometry::Snapshot::QuarticCenterDirection direction) const
    {
        ExactExpansion::validateEnvironment();
        using Direction=::LevelSetHermiteGeometry::Snapshot::QuarticCenterDirection;
        Detail::require(std::isfinite(displacementCells) && displacementCells>0 &&
            std::isfinite(retainedFraction) && retainedFraction>0 && retainedFraction<1 &&
            (direction==Direction::Both || direction==Direction::Decrease || direction==Direction::Increase),
            "invalid Q4 center interval controls");
        const auto& cell=activeCell(index);
        Detail::require(cell.sectorIndex<0,"center interval does not alter a corner sector");
        using A=Detail::Outward;using Net=std::array<ExactExpansion::Interval,25>;
        const auto controls=QuarticCell::bernsteinControlsCleared(cell.local.numerator);
        const std::array<Exact,2> widths{{difference(cell.upper.x,cell.lower.x),difference(cell.upper.y,cell.lower.y)}};
        const double h=std::min(ExactExpansion::enclosure(widths[0]).lo,ExactExpansion::enclosure(widths[1]).lo);
        const double center=value(cell.physical[24]);
        const auto slope=ExactExpansion::multiply(cell.local.denominator,number(1024));
        const auto scalarDenominator=ExactExpansion::multiply(cell.local.denominator,number(144));
        const auto hull=[](const Net& net) {
            ExactExpansion::Interval r{std::numeric_limits<double>::infinity(),-std::numeric_limits<double>::infinity()};
            for(const auto& x:net){if(!A::finite(x))return A::whole();r.lo=std::min(r.lo,x.lo);r.hi=std::max(r.hi,x.hi);}
            return r;
        };
        const auto split=[](const Net& net,int axis) {
            std::array<Net,2> halves;
            for(int line=0;line<5;++line) {
                const auto at=[&](int k){return axis?5*k+line:5*line+k;};
                std::array<ExactExpansion::Interval,5> work,left,right;
                for(int k=0;k<5;++k)work[k]=net[at(k)];
                left[0]=work[0];right[4]=work[4];
                for(int stage=1;stage<5;++stage) {
                    for(int k=0;k<5-stage;++k)work[k]=A::multiply(A::add(work[k],work[k+1]),A::point(.5));
                    left[stage]=work[0];right[4-stage]=work[4-stage];
                }
                for(int k=0;k<5;++k){halves[0][at(k)]=left[k];halves[1][at(k)]=right[k];}
            }
            return halves;
        };
        // Differentiate degree four and elevate degree three back to four.
        // These integer weights retain the exact original-cell derivative.
        constexpr int weights[5][4]={{4,0,0,0},{1,3,0,0},{0,2,2,0},{0,0,3,1},{0,0,0,4}};
        const auto derivativeControls=[&](const auto& values,int axis) {
            std::array<Exact,25> d;
            for(int line=0;line<5;++line) {
                const auto at=[&](int k){return axis?5*k+line:5*line+k;};
                for(int k=0;k<5;++k)for(int j=0;j<4;++j)if(weights[k][j])
                    d[at(k)]=ExactExpansion::add(d[at(k)],ExactExpansion::multiply(
                        ExactExpansion::subtract(values[at(j+1)],values[at(j)]),number(weights[k][j])));
            }
            return d;
        };
        const double permitted=A::multiply(A::point(displacementCells),A::point(h)).lo;
        ::LevelSetHermiteGeometry::Snapshot::QuarticCenterInterval best{center,center,-1,0,0,0};
        for(int axis=0;axis<2;++axis) {
            if(!resolvedZeroSetPartial(index,axis))continue;
            const auto denominator=ExactExpansion::multiply(scalarDenominator,widths[axis]);
            const auto den=ExactExpansion::enclosure(denominator);
            Detail::require(den.lo>0 && std::isfinite(den.hi),"unresolved Q4 physical derivative scale");
            double maximum=0;
            for(const auto& x:derivativeControls(controls,axis)) {
                const auto e=ExactExpansion::enclosure(x);maximum=std::max({maximum,std::abs(e.lo),std::abs(e.hi)});
            }
            double minimum=retainedFraction*maximum/den.hi;
            for(int attempt=0;attempt<16;++attempt,minimum*=.5) {
                minimum=std::nextafter(minimum,0.);
                double extent=std::nextafter(displacementCells*h*minimum,0.);
                if(!(std::isfinite(extent) && extent>0 && std::isfinite(minimum) && minimum>0))continue;
                bool found=false;
                // Preserve the legacy symmetric work cap. A one-sided
                // regular tube may remain safe far beyond the symmetric
                // interval, especially when psi is tiny near a corner.
                // Additional trials spend work only; every extent still
                // passes the identical derivative/displacement proof.
                const int growthLimit=direction==Direction::Both?9:16;
                for(int growth=0;growth<growthLimit;++growth,extent*=2) {
                const double lower=direction==Direction::Increase?center:std::min(center,std::nextafter(center-extent,center));
                const double upper=direction==Direction::Decrease?center:std::max(center,std::nextafter(center+extent,center));
                if(!std::isfinite(lower) || !std::isfinite(upper))break;
                std::array<QuarticCell::ExpandedControls,2> endpoints{{controls,controls}};
                double deltaBound=0;
                for(int end=0;end<2;++end) {
                    const auto delta=ExactExpansion::subtract(number(end?upper:lower),number(center));
                    const auto absolute=ExactExpansion::sign(delta)<0?ExactExpansion::multiply(delta,number(-1)):delta;
                    deltaBound=std::max(deltaBound,ExactExpansion::enclosure(absolute).hi);
                    endpoints[end][12]=ExactExpansion::add(controls[12],ExactExpansion::multiply(slope,delta));
                }
                Net values,partial,bubble{};bubble[12]=ExactExpansion::enclosure(slope);
                const auto da=derivativeControls(endpoints[0],axis),db=derivativeControls(endpoints[1],axis);
                for(int k=0;k<25;++k) {
                    const auto a=ExactExpansion::enclosure(endpoints[0][k]),b=ExactExpansion::enclosure(endpoints[1][k]);
                    values[k]={std::min(a.lo,b.lo),std::max(a.hi,b.hi)};
                    const auto p=ExactExpansion::enclosure(da[k]),q=ExactExpansion::enclosure(db[k]);
                    partial[k]={std::min(p.lo,q.lo),std::max(p.hi,q.hi)};
                }
                const double required=ExactExpansion::enclosure(ExactExpansion::multiply(denominator,number(minimum))).hi;
                bool enlarged=false;
                for(int sign:{1,-1}) {
                    // Interval controls cover EVERY center between the two
                    // endpoints. Only boxes capable of containing a zero for
                    // ANY of those centers need a signed derivative bound.
                    // Pure parts may contain stationary scalar points.
                    struct Box{Net value,partial,bubble;int depth;};
                    std::vector<Box> pending{{values,partial,bubble,0}};int visited=0;bool certified=true;double displacementRate=0;
                    while(!pending.empty()) {
                        if(++visited>16384){certified=false;break;}
                        auto box=std::move(pending.back());pending.pop_back();
                        const auto f=hull(box.value);if(f.lo>0 || f.hi<0)continue;
                        const auto d=hull(box.partial);
                        const double signedMinimum=sign>0?d.lo:-d.hi;
                        if(signedMinimum>=required) {
                            // A bubble vanishes quadratically on the skeleton.
                            // Its actual zero-tube size is essential for small
                            // corner cuts; using psi<=1 everywhere spuriously
                            // rejects an otherwise tiny interface displacement.
                            const auto b=hull(box.bubble);
                            const auto numerator=A::multiply(A::point(std::max(0.,b.hi)),ExactExpansion::enclosure(widths[axis]));
                            const double rate=std::nextafter(numerator.hi/signedMinimum,std::numeric_limits<double>::infinity());
                            if(std::isfinite(rate) && A::multiply(A::point(deltaBound),A::point(rate)).hi<=permitted) {
                                displacementRate=std::max(displacementRate,rate);continue;
                            }
                        }
                        if(box.depth==24){certified=false;break;}
                        const auto fv=split(box.value,box.depth&1),dv=split(box.partial,box.depth&1),bv=split(box.bubble,box.depth&1);
                        pending.push_back({fv[1],dv[1],bv[1],box.depth+1});pending.push_back({fv[0],dv[0],bv[0],box.depth+1});
                    }
                    if(certified) {
                        if(upper-lower>best.upper-best.lower)best={lower,upper,axis,sign,minimum,displacementRate};
                        found=true;enlarged=true;break;
                    }
                }
                if(!enlarged)break;
                }
                if(found)break;
            }
        }
        if(best.partialAxis<0)throw InterfaceGeometryResolution(InterfaceGeometryResolution::Reason::CriticalCell,index,
            "bubble center interval has no resolved signed partial on its possible zero-set tube");
        return best;
    }
    ::LevelSetHermiteGeometry::Snapshot::QuarticChangeCertificate changeCertificate(int index,
        const QuarticLeafData& other,double budget) const
    {
        ExactExpansion::validateEnvironment();
        Detail::require(std::isfinite(budget) && budget>0,"invalid Q4 change displacement budget");
        const auto& a=activeCell(index);const auto& b=other.activeCell(index);
        Detail::require(a.lower.x==b.lower.x && a.lower.y==b.lower.y && a.upper.x==b.upper.x && a.upper.y==b.upper.y &&
            a.sectorIndex<0 && b.sectorIndex<0,"Q4 change certificate requires identical sector-free cells");
        using A=Detail::Outward;using Net=std::array<ExactExpansion::Interval,25>;
        auto ca=QuarticCell::bernsteinControlsCleared(a.local.numerator);
        auto cb=QuarticCell::bernsteinControlsCleared(b.local.numerator);
        auto denominator=a.local.denominator;
        if(ExactExpansion::sign(ExactExpansion::subtract(a.local.denominator,b.local.denominator))!=0) {
            for(int k=0;k<25;++k) {
                ca[k]=ExactExpansion::multiply(ca[k],b.local.denominator);
                cb[k]=ExactExpansion::multiply(cb[k],a.local.denominator);
            }
            denominator=ExactExpansion::multiply(a.local.denominator,b.local.denominator);
        }
        const auto den=ExactExpansion::enclosure(ExactExpansion::multiply(denominator,number(144)));
        Detail::require(den.lo>0 && A::finite(den),"unresolved Q4 change denominator");
        std::array<ExactExpansion::Interval,2> width{{ExactExpansion::enclosure(difference(a.upper.x,a.lower.x)),
            ExactExpansion::enclosure(difference(a.upper.y,a.lower.y))}};
        const double h=std::min(width[0].lo,width[1].lo);
        const auto hull=[](const Net& n){ExactExpansion::Interval r{INFINITY,-INFINITY};
            for(const auto x:n){if(!A::finite(x))return A::whole();r.lo=std::min(r.lo,x.lo);r.hi=std::max(r.hi,x.hi);}return r;};
        const auto split=[](const Net& n,int axis) {
            std::array<Net,2> out;
            for(int line=0;line<5;++line) {
                const auto at=[&](int j){return axis?5*j+line:5*line+j;};std::array<ExactExpansion::Interval,5> work;
                for(int j=0;j<5;++j)work[j]=n[at(j)];out[0][at(0)]=work[0];out[1][at(4)]=work[4];
                for(int stage=1;stage<5;++stage) {
                    for(int j=0;j<5-stage;++j)work[j]=A::multiply(A::add(work[j],work[j+1]),A::point(.5));
                    out[0][at(stage)]=work[0];out[1][at(4-stage)]=work[4-stage];
                }
            }return out;
        };
        constexpr int weights[5][4]={{4,0,0,0},{1,3,0,0},{0,2,2,0},{0,0,3,1},{0,0,0,4}};
        const auto derivative=[&](const auto& n,int axis) {
            std::array<Exact,25> out;
            for(int line=0;line<5;++line)for(int k=0;k<5;++k)for(int j=0;j<4;++j)if(weights[k][j]) {
                const auto at=[&](int p){return axis?5*p+line:5*line+p;};
                out[at(k)]=ExactExpansion::add(out[at(k)],ExactExpansion::multiply(
                    ExactExpansion::subtract(n[at(j+1)],n[at(j)]),number(weights[k][j])));
            }return out;
        };
        struct Box{Net value,delta,dx,dy;int depth=0;};Box initial;
        const auto ax=derivative(ca,0),ay=derivative(ca,1),bx=derivative(cb,0),by=derivative(cb,1);
        for(int k=0;k<25;++k) {
            const auto span=[](const Exact& x,const Exact& y){const auto p=ExactExpansion::enclosure(x),q=ExactExpansion::enclosure(y);
                return ExactExpansion::Interval{std::min(p.lo,q.lo),std::max(p.hi,q.hi)};};
            initial.value[k]=span(ca[k],cb[k]);initial.delta[k]=ExactExpansion::enclosure(ExactExpansion::subtract(cb[k],ca[k]));
            initial.dx[k]=span(ax[k],bx[k]);initial.dy[k]=span(ay[k],by[k]);
        }
        ::LevelSetHermiteGeometry::Snapshot::QuarticChangeCertificate result;
        result.minimumPhysicalPartial=INFINITY;std::vector<Box> pending{initial};
        while(!pending.empty()) {
            if(++result.visitedBoxes>16384)return result;
            auto box=std::move(pending.back());pending.pop_back();const auto value=hull(box.value);
            if(value.lo>0 || value.hi<0)continue;
            const auto delta=hull(box.delta);const double absolute=std::max(std::abs(delta.lo),std::abs(delta.hi));
            double bestRate=INFINITY,partial=0;
            for(int axis=0;axis<2;++axis) {
                const auto d=hull(axis?box.dy:box.dx);const double lower=d.lo>0?d.lo:d.hi<0?-d.hi:0;
                if(lower<=0)continue;
                const double rate=A::up(A::multiply(A::point(absolute),width[axis]).hi/lower);
                if(rate<bestRate){bestRate=rate;partial=A::down(lower/A::multiply(den,width[axis]).hi);}
            }
            if(std::isfinite(bestRate) && partial>0 && bestRate<=A::multiply(A::point(budget),A::point(h)).lo) {
                result.maxDisplacementCells=std::max(result.maxDisplacementCells,A::up(bestRate/h));
                result.minimumPhysicalPartial=std::min(result.minimumPhysicalPartial,partial);continue;
            }
            if(box.depth==24)return result;
            const int axis=box.depth&1;const auto f=split(box.value,axis),d=split(box.delta,axis),x=split(box.dx,axis),y=split(box.dy,axis);
            for(int half=1;half>=0;--half)pending.push_back({f[half],d[half],x[half],y[half],box.depth+1});
        }
        result.resolved=true;if(!std::isfinite(result.minimumPhysicalPartial))result.minimumPhysicalPartial=0;
        return result;
    }
    int boxSign(int index,Vector2D lower,Vector2D upper) const
    {
        ExactExpansion::validateEnvironment();const auto& cell=activeCell(index);
        const auto chart=[](double p,double lo,double hi,bool upperBound) {
            Detail::require(std::isfinite(p) && p>=lo && p<=hi,"Q4 box leaves its native cell");
            const auto n=ExactExpansion::enclosure(difference(p,lo));
            const auto d=ExactExpansion::enclosure(difference(hi,lo));
            Detail::require(d.lo>0 && std::isfinite(d.hi),"Q4 box has unresolved cell width");
            const double value=upperBound?std::nextafter(n.hi/d.lo,std::numeric_limits<double>::infinity()):
                std::nextafter(n.lo/d.hi,-std::numeric_limits<double>::infinity());
            return std::clamp(value,0.,1.);
        };
        Detail::require(lower.x<=upper.x && lower.y<=upper.y,"unordered Q4 physical box");
        const Vector2D a{chart(lower.x,cell.lower.x,cell.upper.x,false),chart(lower.y,cell.lower.y,cell.upper.y,false)};
        const Vector2D b{chart(upper.x,cell.lower.x,cell.upper.x,true),chart(upper.y,cell.lower.y,cell.upper.y,true)};
        return localBoxSign(index,a,b);
    }
    PolynomialAffineCoordinate affineEnclosure(int index,const std::array<PolynomialAffineCoordinate,2>& physical) const
    {
        ExactExpansion::validateEnvironment();const auto& cell=activeCell(index);
        // These are the original homogeneous constrained data, never the
        // rounded cache returned by mappedCellLocalDofs.
        auto powers=QuarticCell::bernsteinControlsCleared(cell.local.numerator);
        constexpr int binomial[5][5]={{1,0,0,0,0},{1,1,0,0,0},{1,2,1,0,0},{1,3,3,1,0},{1,4,6,4,1}};
        for(int axis=0;axis<2;++axis){const auto before=powers;
            for(int line=0;line<5;++line)for(int k=0;k<5;++k){Exact coefficient;
                const auto entry=[&](int q){return axis?5*q+line:5*line+q;};
                for(int j=0;j<=k;++j){const int factor=binomial[4][k]*binomial[k][j]*((k-j)%2?-1:1);
                    coefficient=ExactExpansion::add(coefficient,ExactExpansion::multiply(before[entry(j)],ExactExpansion::fromDouble(factor)));}
                powers[entry(k)]=std::move(coefficient);
            }
        }
        auto denominator=ExactExpansion::multiply(cell.local.denominator,ExactExpansion::fromDouble(144));
        ExactExpansion::centerHomogeneous(powers,denominator);
        const auto d=ExactExpansion::enclosure(denominator);
        Detail::require(d.lo>0 && std::isfinite(d.hi),"unresolved authoritative Q4 affine denominator");
        std::array<std::array<double,2>,25> coefficients{};
        for(std::size_t k=0;k<25;++k){
            if(ExactExpansion::sign(powers[k])==0)continue;
            const auto numerator=ExactExpansion::enclosure(powers[k]);
            // Positive denominator: extrema occur at endpoint quotients.
            // One adjacent binary64 value bounds each correctly rounded
            // division, including gradual underflow.
            const std::array<double,4> quotients{{numerator.lo/d.lo,numerator.lo/d.hi,
                numerator.hi/d.lo,numerator.hi/d.hi}};
            coefficients[k]={std::nextafter(*std::min_element(quotients.begin(),quotients.end()),-std::numeric_limits<double>::infinity()),
                std::nextafter(*std::max_element(quotients.begin(),quotients.end()),std::numeric_limits<double>::infinity())};
        }
        return encloseTensorQuarticAffine(coefficients,cell.lower,cell.upper,physical);
    }
    int localBoxSign(int index,Vector2D a,Vector2D b) const
    {
        ExactExpansion::validateEnvironment();
        Detail::require(a.x>=0 && a.y>=0 && b.x<=1 && b.y<=1 && a.x<=b.x && a.y<=b.y,"invalid native Q4 local sign box");
        const auto& cell=activeCell(index);
        std::call_once(cell.boxControlsReady,[&] {
            const auto exact=QuarticCell::bernsteinControlsCleared(cell.local.numerator);
            double lo=std::numeric_limits<double>::infinity(),hi=-lo;
            for(std::size_t k=0;k<25;++k) {
                cell.boxControls[k]=ExactExpansion::enclosure(exact[k]);
                if(!Detail::Outward::finite(cell.boxControls[k])) {lo=-std::numeric_limits<double>::infinity();hi=-lo;}
                else {lo=std::min(lo,cell.boxControls[k].lo);hi=std::max(hi,cell.boxControls[k].hi);}
            }
            cell.boxControlSign=lo>0?1:hi<0?-1:0;
        });
        if(cell.boxControlSign)return cell.boxControlSign;
        const auto hull=Detail::tensorBernsteinBounds(cell.boxControls,a,b);
        return hull.lo>0?1:hull.hi<0?-1:0;
    }
    // Remeshing identity is a polynomial identity, not agreement of sampled
    // values or rounded jet caches. Bounds must match so local charts retain
    // their physical meaning. This includes hanging constraints and sectors.
    bool identicalCell(int index,const QuarticLeafData& other,int otherIndex) const
    {
        ExactExpansion::validateEnvironment();
        const auto& a=activeCell(index);const auto& b=other.activeCell(otherIndex);
        if(a.lower.x!=b.lower.x || a.lower.y!=b.lower.y || a.upper.x!=b.upper.x
            || a.upper.y!=b.upper.y || radial_!=other.radial_)return false;
        for(std::size_t k=0;k<25;++k)
            if(ExactExpansion::sign(ExactExpansion::subtract(
                ExactExpansion::multiply(a.local.numerator[k],b.local.denominator),
                ExactExpansion::multiply(b.local.numerator[k],a.local.denominator)))!=0)return false;
        return true;
    }
    bool independentMidpoint(const std::array<int,2>& vertices) const
    {
        const auto found=edgeIndex_.find(vertices);
        Detail::require(found!=edgeIndex_.end(),"midpoint ownership requires an active native side");
        return edges_[found->second].parent<0;
    }
    bool identicalCellRepresentation(int index,const QuarticLeafData& other,int otherIndex) const
    {
        ExactExpansion::validateEnvironment();
        const auto& a=activeCell(index); const auto& b=other.activeCell(otherIndex);
        const std::array<double,4> boundsA{{a.lower.x,a.lower.y,a.upper.x,a.upper.y}};
        const std::array<double,4> boundsB{{b.lower.x,b.lower.y,b.upper.x,b.upper.y}};
        if(radial_!=other.radial_ || std::memcmp(boundsA.data(),boundsB.data(),sizeof(double)*4)!=0 ||
            a.local.denominator.components()!=b.local.denominator.components()) return false;
        for(std::size_t k=0;k<25;++k)
            if(a.local.numerator[k].components()!=b.local.numerator[k].components()) return false;
        const auto& rawA=a.cache->localDofs(); const auto& rawB=b.cache->localDofs();
        return std::memcmp(rawA.data(),rawB.data(),sizeof(double)*rawA.size())==0;
    }
    Sample evaluateCellChart(int index,CellChartPoint point) const
    { return activeCell(index).cache->approximateJetAtChart(point); }
    ThirdDerivatives thirdDerivativesCellChart(int index,CellChartPoint point) const
    { return activeCell(index).cache->approximateThirdDerivativesAtChart(point); }
    int cellSignExceptVertices(int index,unsigned allowedMask) const
    {
        Detail::require((allowedMask&~15u)==0,"invalid quartic corner mask");
        const auto& cell=activeCell(index);
        // A previously cached integer is not permission to query in an
        // invalid FP environment. Retain this check on EVERY public use.
        ExactExpansion::validateEnvironment();
        std::call_once(cell.phaseSignsReady,[&cell] {
            const auto controls=QuarticCell::bernsteinControlsCleared(cell.local.numerator);
            std::array<int,25> signs;
            for(std::size_t k=0;k<25;++k) signs[k]=ExactExpansion::sign(controls[k]);
            cell.phaseSigns=signs;
        });
        const int phase=CellView::certifyBernsteinSigns(cell.phaseSigns,allowedMask);
        if (phase!=0) return phase;
        std::call_once(cell.resolvedPhaseReady,[&] {
            // A monotone cell can be classified by its authoritative edge
            // inventory. Only the nonmonotone, inconclusive cells need this
            // bounded pure-phase certificate. Never infer phase from samples.
            if (resolvedCartesianPartial(index)) return;
            const auto exact=QuarticCell::bernsteinControlsCleared(cell.local.numerator);
            std::array<ExactExpansion::Interval,25> controls;
            for (std::size_t k=0;k<25;++k) controls[k]=ExactExpansion::enclosure(exact[k]);
            struct Box { Vector2D lower,upper; int depth; };
            std::vector<Box> pending{{{0,0},{1,1},0}};
            int commonSign=0,visited=0;
            while (!pending.empty()) {
                if (++visited>4096) return;
                const auto box=pending.back();pending.pop_back();
                const auto hull=Detail::tensorBernsteinBounds(controls,box.lower,box.upper);
                const int sign=hull.lo>0?1:hull.hi<0?-1:0;
                if (sign) {
                    if (commonSign && commonSign!=sign) return;
                    commonSign=sign;continue;
                }
                if (box.depth==12) return;
                auto first=box,second=box;first.depth=second.depth=box.depth+1;
                if ((box.depth&1)==0) {
                    const double mid=.5*(box.lower.x+box.upper.x);
                    first.upper.x=second.lower.x=mid;
                } else {
                    const double mid=.5*(box.lower.y+box.upper.y);
                    first.upper.y=second.lower.y=mid;
                }
                pending.push_back(second);pending.push_back(first);
            }
            // Strict same-sign hulls cover the CLOSED original cell. The
            // original positive homogeneous denominator is not approximated;
            // subdividing its numerator's Bernstein representation changes
            // neither the scalar, the zero set, nor the physical mesh.
            cell.resolvedPhase=commonSign;
        });
        return cell.resolvedPhase;
    }
    std::array<int,2> zeroCornerGradientSigns(int index,int corner) const
    {
        Detail::require(corner>=0 && corner<4,"invalid native quartic corner");
        constexpr int x[4]={0,1,1,0},y[4]={0,0,1,1};
        const auto& dofs=activeCell(index).local.numerator;
        Detail::require(ExactExpansion::sign(dofs[5*y[corner]+x[corner]])==0,
            "native quartic corner-gradient query requires an exact zero corner");
        return {{ExactExpansion::sign(dofs[5*y[corner]+x[corner]+2]),
            ExactExpansion::sign(dofs[5*(y[corner]+2)+x[corner]])}};
    }
    bool resolvedCartesianPartial(int index, int axis = -1) const
    {
        const auto controls=activeCell(index).cache->bernsteinControlsApproximate();
        long double scale=0;
        for(const auto value:controls) scale=std::max(scale,std::abs(value));
        return LevelSetHermiteScalar::resolvedCartesianPartial(controls,
            512.0L*std::numeric_limits<double>::epsilon()*scale, axis);
    }
    bool resolvedZeroSetPartial(int index,int axis=-1,
        Vector2D lower={0,0},Vector2D upper={1,1}) const
    {
        ExactExpansion::validateEnvironment();
        Detail::require(axis>=-1 && axis<2,"invalid zero-set partial axis");
        const auto& cell=activeCell(index);
        if(lower.x!=0 || lower.y!=0 || upper.x!=1 || upper.y!=1) {
            const auto partials=QuarticCell::zeroSetCartesianPartials(
                QuarticCell::bernsteinControlsCleared(cell.local.numerator),lower,upper);
            return axis<0?partials[0] || partials[1]:partials[axis];
        }
        std::call_once(cell.zeroPartialsReady,[&] {
            cell.zeroPartials=QuarticCell::zeroSetCartesianPartials(
                QuarticCell::bernsteinControlsCleared(cell.local.numerator));
        });
        return axis<0?cell.zeroPartials[0] || cell.zeroPartials[1]:cell.zeroPartials[axis];
    }
    std::vector<QuarticCell::RegularChartBox> regularChartPartition(int index) const
    {
        ExactExpansion::validateEnvironment();
        return QuarticCell::regularChartPartition(
            QuarticCell::bernsteinControlsCleared(activeCell(index).local.numerator));
    }
    std::array<int,2> localGradientSigns(int index,Vector2D point) const
    {
        const RationalSegmentCoordinate x{number(point.x),number(point.x),number(1)},
            y{number(point.y),number(point.y),number(1)};
        return {{ExactExpansion::sign(restrictLocal(index,x,y,1,0).numerator[0]),
            ExactExpansion::sign(restrictLocal(index,x,y,0,1).numerator[0])}};
    }
    long double derivativeResolutionScale(int index) const
    {
        const auto& cell=activeCell(index);
        const auto controls=cell.cache->bernsteinControlsApproximate();
        const auto limits=std::minmax_element(controls.begin(),controls.end());
        const auto scale=(*limits.second-*limits.first)/std::min(cell.upper.x-cell.lower.x,cell.upper.y-cell.lower.y);
        Detail::require(std::isfinite(scale),"unrepresentable quartic derivative-resolution scale");
        return scale;
    }
    bool closedCellOutsideBand(int index,double halfWidth) const
    {
        return classifyBand(index,halfWidth).outside;
    }
    bool closedCellInsideBand(int index,double halfWidth) const
    {
        return classifyBand(index,halfWidth).inside;
    }
    Polynomial restrictPhysical(int index,Vector2D a,Vector2D b,unsigned dx=0,unsigned dy=0) const
    {
        const auto& cell=activeCell(index);
        return restrictLocal(index,chart(a.x,b.x,cell.lower.x,cell.upper.x),
            chart(a.y,b.y,cell.lower.y,cell.upper.y),dx,dy);
    }
    Polynomial restrictLocal(int index,const RationalSegmentCoordinate& x,const RationalSegmentCoordinate& y,
        unsigned dx=0,unsigned dy=0) const
    {
        const auto& cell=activeCell(index);
        auto restricted=QuarticCell::homogeneousRestriction(cell.local.numerator,cell.local.denominator,x,y,dx,dy);
        Polynomial result{std::move(restricted.numerator),std::move(restricted.denominator)};
        for(unsigned k=0;k<dx;++k)
            result.denominator=ExactExpansion::multiply(result.denominator,difference(cell.upper.x,cell.lower.x));
        for(unsigned k=0;k<dy;++k)
            result.denominator=ExactExpansion::multiply(result.denominator,difference(cell.upper.y,cell.lower.y));
        ExactExpansion::centerHomogeneous(result.numerator,result.denominator);
        return result;
    }

    // Physical 5*y+x cell columns -> independent vertex/edge/center data.
    // The DAG is visited once per touched datum in reverse construction order;
    // cost follows local support, not all mesh vertices. Pinning/axis and the
    // choice of free input data remain the caller's admissibility constraints.
    std::vector<Sensitivity> pullbackNodeSensitivity(int index,int component) const
    {
        ExactExpansion::validateEnvironment();
        Detail::require(index>=0 && static_cast<std::size_t>(index)<vertices_.size() && vertices_[index].used
            && component>=0 && component<4,"quartic nodal sensitivity requires a used vertex and component 0..3");
        return pullbackDataSensitivity({{vertices_[index].jet[component],1.0L}});
    }
    std::vector<Sensitivity> pullbackCellSensitivity(int index,const std::array<long double,25>& columns,
        CornerSectorVariation sectors=CornerSectorVariation::RequireAbsent) const
    {
        ExactExpansion::validateEnvironment();
        const auto& cell=activeCell(index);
        Detail::require(sectors==CornerSectorVariation::RequireAbsent || sectors==CornerSectorVariation::HoldFixed,
            "invalid native quartic corner-sector variation policy");
        Detail::require(cell.sectorIndex<0 || sectors==CornerSectorVariation::HoldFixed,
            "native quartic sector sensitivity requires an explicit sector-DOF variation model");
        std::map<int,long double> pending;
        for(std::size_t k=0;k<25;++k) add(pending[cell.physical[k]],columns[k],1.0L);
        return pullbackDataSensitivity(std::move(pending));
    }

  private:
    std::vector<Sensitivity> pullbackDataSensitivity(std::map<int,long double> pending) const
    {
        std::map<std::tuple<int,int,int>,long double> independentColumns;
        while(!pending.empty())
        {
            const auto last=std::prev(pending.end());
            const int id=last->first;
            const auto coefficient=last->second;
            pending.erase(last);
            if(coefficient==0) continue;
            const auto& datum=data_[id];
            if(datum.fixed)continue;
            if(datum.dependencies.empty())
            {
                const auto& identity=datum.identity;
                add(independentColumns[{static_cast<int>(identity.kind),identity.index,identity.component}],
                    coefficient,1.0L);
            }
            else for(const auto& dependency:datum.dependencies)
                add(pending[dependency.source],coefficient,dependency.weight);
        }
        std::vector<Sensitivity> result;
        for(const auto& entry:independentColumns) if(entry.second!=0)
            result.push_back({static_cast<Sensitivity::Kind>(std::get<0>(entry.first)),std::get<1>(entry.first),
                std::get<2>(entry.first),entry.second});
        return result;
    }

    BandClassification classifyBand(int index,double halfWidth) const
    {
        ExactExpansion::validateEnvironment();
        Detail::require(std::isfinite(halfWidth) && halfWidth>=0,
            "native quartic band width must be nonnegative and finite");
        const auto& cell=activeCell(index);
        const std::lock_guard<std::mutex> lock(cell.bandMutex);
        for(const auto& band:cell.bands)
            if(band.ready && std::memcmp(&band.halfWidth,&halfWidth,sizeof(double))==0) return band;
        // Preserve both original shifted-control constructions and their
        // sign predicates. No approximate bounds or different certificate.
        const auto high=shiftedControlSigns(index,-halfWidth),low=shiftedControlSigns(index,halfWidth);
        BandClassification result;
        result.halfWidth=halfWidth;
        result.inside=std::all_of(high.begin(),high.end(),[](int sign){return sign<=0;})
            && std::all_of(low.begin(),low.end(),[](int sign){return sign>=0;});
        result.outside=CellView::certifyBernsteinSigns(high,0)==1 || CellView::certifyBernsteinSigns(low,0)==-1;
        result.ready=true;
        cell.bands[cell.nextBand]=result;
        cell.nextBand=(cell.nextBand+1)%cell.bands.size();
        return result;
    }

    std::array<int,25> shiftedControlSigns(int index,double shift) const
    {
        Detail::require(std::isfinite(shift),"native quartic band shift must be finite");
        const auto& local=activeCell(index).local;
        auto data=local.numerator;
        const auto amount=ExactExpansion::multiply(number(shift),local.denominator);
        for(int y:{0,1,4}) for(int x:{0,1,4})
            data[5*y+x]=ExactExpansion::add(data[5*y+x],amount);
        const auto controls=QuarticCell::bernsteinControlsCleared(data);
        std::array<int,25> result;
        for(std::size_t k=0;k<25;++k) result[k]=ExactExpansion::sign(controls[k]);
        return result;
    }
    CartesianLeafLocator locator_; // validates and owns actual leaf coordinates
    std::vector<VertexData> vertices_;
    std::vector<EdgeData> edges_;
    std::map<EdgeKey,int> edgeIndex_;
    std::vector<Leaf> cells_;
    std::vector<Datum> data_;
    Vector2D lower_,upper_;
    bool radial_=false;
    int hangingCount_=0,constrainedEdgeCount_=0,sectorCount_=0;
    std::size_t independentCount_=0;
    // Construction-only borrows, cleared before successful construction returns.
    const QuarticLeafData* restrictionSource_=nullptr;
    const Patch* restrictionPatch_=nullptr;

    void prepareTopology(const QuadtreeMesh& mesh,const std::vector<std::uint8_t>& active)
    {
        for(std::size_t i=0;i<vertices_.size();++i)
        {
            const auto& source=mesh.verts(mesh.vertIndexToLabel[i]);
            vertices_[i].point=source.coord();
            if(source.supCell>=0) vertices_[i].coarseCell=mesh.cells(source.supCell).index();
        }
        bool any=false;
        constexpr int cornerOrder[4]={0,1,3,2};
        constexpr int sideCorners[4][2]={{0,1},{2,3},{0,2},{1,3}};
        for(std::size_t i=0;i<cells_.size();++i)
        {
            Detail::require(active[i]<=1,"native quartic activity mask must be binary");
            auto& cell=cells_[i];
            const auto& source=mesh.cells(mesh.cellIndexToLabel[i]);
            cell.active=active[i]!=0;
            for(int k=0;k<4;++k)
                cell.vertices[k]=mesh.verts(source.cellToVert[cornerOrder[k]]).index();
            cell.lower=vertices_[cell.vertices[0]].point;
            cell.upper=vertices_[cell.vertices[3]].point;
            if(!cell.active) continue;
            any=true;
            for(int vertex:cell.vertices) vertices_[vertex].used=true;
            for(int side=0;side<4;++side)
            {
                const EdgeKey key{{cell.vertices[sideCorners[side][0]],cell.vertices[sideCorners[side][1]]}};
                const auto inserted=edgeIndex_.emplace(key,static_cast<int>(edges_.size()));
                if(inserted.second)
                {
                    EdgeData edge;
                    edge.vertices=key; edge.horizontal=side<2;
                    edges_.push_back(edge);
                }
                cell.edges[side]=inserted.first->second;
            }
        }
        Detail::require(any,"native quartic patch has no active leaves");
        for(auto& vertex:vertices_) if(vertex.used && vertex.coarseCell>=0)
        {
            const auto& coarse=activeCell(vertex.coarseCell);
            const auto p=vertex.point;
            int side=-1;
            if(p.y==coarse.lower.y) side=0;
            else if(p.y==coarse.upper.y) side=1;
            else if(p.x==coarse.lower.x) side=2;
            else if(p.x==coarse.upper.x) side=3;
            Detail::require(side>=0,"quartic hanging vertex is not on its actual coarse side");
            vertex.coarseEdge=coarse.edges[side];
            const auto& edge=edges_[vertex.coarseEdge];
            const auto a=vertices_[edge.vertices[0]].point,b=vertices_[edge.vertices[1]].point;
            Detail::require(edge.horizontal ? (a.x<p.x && p.x<b.x) : (a.y<p.y && p.y<b.y),
                "quartic hanging vertex is not strictly inside its coarse side");
            ++hangingCount_;
        }
        for(auto& edge:edges_)
        {
            const auto a=vertices_[edge.vertices[0]].point,b=vertices_[edge.vertices[1]].point;
            for(int vertex:edge.vertices)
            {
                const int parent=vertices_[vertex].coarseEdge;
                if(parent<0) continue;
                const auto& coarse=edges_[parent];
                if(coarse.horizontal!=edge.horizontal) continue;
                const auto lo=vertices_[coarse.vertices[0]].point,hi=vertices_[coarse.vertices[1]].point;
                const bool contained=edge.horizontal
                    ? (a.y==lo.y && b.y==lo.y && lo.x<=a.x && b.x<=hi.x && b.x-a.x<hi.x-lo.x)
                    : (a.x==lo.x && b.x==lo.x && lo.y<=a.y && b.y<=hi.y && b.y-a.y<hi.y-lo.y);
                if(!contained) continue;
                Detail::require(edge.parent<0 || edge.parent==parent,"quartic fine side has conflicting coarse owners");
                edge.parent=parent;
            }
            if(edge.parent>=0) ++constrainedEdgeCount_;
        }
    }

    void prepareSectors(const QuadtreeMesh& mesh,const Patch& patch)
    {
        std::vector<bool> declared(vertices_.size(),false);
        for(const auto& declaration:patch.cornerSectors)
        {
            const int cellIndex=declaration.cellIndex,vertexIndex=declaration.vertexIndex;
            Detail::require(cellActive(cellIndex) && vertexIndex>=0 && static_cast<std::size_t>(vertexIndex)<vertices_.size(),
                "native quartic corner sector requires active cell and valid native vertex");
            auto& cell=cells_[cellIndex];
            const auto corner=std::find(cell.vertices.begin(),cell.vertices.end(),vertexIndex);
            Detail::require(corner!=cell.vertices.end() && !declared[vertexIndex] && cell.sectorIndex<0,
                "native quartic corner sector is not a unique incident vertex/cell");
            const auto shared=nativeCornerSharedComponents(mesh,locator_,cellIndex,vertexIndex,
                [&](int index){return cells_[index].active;});
            const auto& d=declaration.dof;
            const std::array<double,4> values{{d.value,d.dx,d.dy,d.dxy}};
            const auto& base=vertices_[vertexIndex].jet;
            Detail::require(d.value==0 && ExactExpansion::sign(data_[base[0]].value.numerator[0])==0,
                "native quartic corner sector must preserve the exact base zero");
            Detail::require(declaration.phaseSign==-1 || declaration.phaseSign==1,
                "native quartic corner sector phase must be minus or plus one");
            for(int k=0;k<4;++k)
            {
                Datum fixed;fixed.fixed=true;fixed.value.numerator[0]=number(values[k]);
                if(shared[k])
                {
                    const auto& original=data_[base[k]].value;
                    Detail::require(ExactExpansion::sign(ExactExpansion::subtract(
                        ExactExpansion::multiply(fixed.value.numerator[0],original.denominator),original.numerator[0]))==0,
                        "native quartic corner sector violates the exact active-edge C0 trace");
                }
                cell.sectorData[k]=static_cast<int>(data_.size());data_.push_back(std::move(fixed));
            }
            cell.sectorIndex=sectorCount_++;
            cell.sectorCorner=static_cast<int>(corner-cell.vertices.begin());
            cell.sectorPhaseSign=declaration.phaseSign;
            declared[vertexIndex]=true;
        }
    }
    void resolveVertex(int index,const Patch& patch)
    {
        auto& vertex=vertices_[index];
        if(vertex.state==2) return;
        Detail::require(vertex.state==0,"cyclic native quartic vertex constraint");
        vertex.state=1;
        if(vertex.coarseEdge>=0)
        {
            const auto& edge=edges_[vertex.coarseEdge];
            const auto a=vertices_[edge.vertices[0]].point,b=vertices_[edge.vertices[1]].point;
            const double p=edge.horizontal?vertex.point.x:vertex.point.y;
            const auto coordinate=chart(p,p,edge.horizontal?a.x:a.y,edge.horizontal?b.x:b.y);
            vertex.jet=traceJet(vertex.coarseEdge,coordinate.numeratorBegin,coordinate.denominator,patch);
        }
        else
        {
            const auto& node=patch.nodes[index];
            Detail::require(!radial_ || vertex.point.x!=0.0 || (node.dx==0.0 && node.dxy==0.0),
                "native quartic axis dx and dxy must be exactly zero");
            const std::array<double,4> values{{node.value,node.dx,node.dy,node.dxy}};
            for(int k=0;k<4;++k) vertex.jet[k]=independent(values[k],Sensitivity::Kind::Vertex,index,k);
        }
        vertex.state=2;
    }
    void resolveMidpoint(int index,const Patch& patch)
    {
        auto& edge=edges_[index];
        if(edge.state==2) return;
        Detail::require(edge.state==0,"cyclic native quartic edge constraint");
        edge.state=1;
        if(edge.parent<0)
        {
            Detail::require(edge.input>=0,"missing independent coarse/shared quartic midpoint");
            const auto& input=patch.edgeMidpoints[edge.input];
            Detail::require(!radial_ || edge.horizontal || vertices_[edge.vertices[0]].point.x!=0.0
                || input.normalDerivative==0.0,"native quartic axis midpoint normal derivative must be zero");
            edge.midpoint={{independent(input.value,Sensitivity::Kind::EdgeMidpoint,edge.input,0),
                independent(input.normalDerivative,Sensitivity::Kind::EdgeMidpoint,edge.input,1)}};
        }
        else
        {
            const auto& parent=edges_[edge.parent];
            const auto a=vertices_[edge.vertices[0]].point,b=vertices_[edge.vertices[1]].point;
            const auto lo=vertices_[parent.vertices[0]].point,hi=vertices_[parent.vertices[1]].point;
            const double low=edge.horizontal?lo.x:lo.y,high=edge.horizontal?hi.x:hi.y;
            // Exact midpoint of actual fine endpoints, not rounded (a+b)/2.
            const auto numerator=ExactExpansion::add(difference(edge.horizontal?a.x:a.y,low),
                difference(edge.horizontal?b.x:b.y,low));
            const auto denominator=ExactExpansion::multiply(number(2),difference(high,low));
            const auto jet=traceJet(edge.parent,numerator,denominator,patch);
            edge.midpoint={{jet[0],jet[edge.horizontal?2:1]}};
        }
        edge.state=2;
    }
    std::array<int,4> traceJet(int edgeIndex,const Exact& numerator,const Exact& denominator,const Patch& patch)
    {
        const auto& edge=edges_[edgeIndex];
        resolveVertex(edge.vertices[0],patch); resolveVertex(edge.vertices[1],patch);
        resolveMidpoint(edgeIndex,patch);
        std::array<int,4> result;
        visitTrace(edge,numerator,denominator,[&](int component,const auto& inputs,const auto& weights,const auto& divisor) {
            result[component]=linear(inputs,weights,divisor);
        });
        return result;
    }
    template<class Apply> void visitTrace(const EdgeData& edge,const Exact& numerator,const Exact& denominator,const Apply& apply) const
    {
        const auto& first=vertices_[edge.vertices[0]].jet;
        const auto& last=vertices_[edge.vertices[1]].jet;
        const auto a=vertices_[edge.vertices[0]].point,b=vertices_[edge.vertices[1]].point;
        const auto width=difference(edge.horizontal?b.x:b.y,edge.horizontal?a.x:a.y);
        const auto basis=QuarticCell::rationalBasisJet(numerator,denominator);
        const auto commonDenominator=ExactExpansion::multiply(basis.denominator,width);
        const int tangent=edge.horizontal?1:2,normal=edge.horizontal?2:1;
        for(int group=0;group<2;++group)
        {
            const int valueIndex=group?normal:0,derivativeIndex=group?3:tangent;
            const std::array<int,5> inputs{{first[valueIndex],last[valueIndex],first[derivativeIndex],
                last[derivativeIndex],edge.midpoint[group]}};
            std::array<Exact,5> values=basis.value,derivatives=basis.derivative;
            for(int k=0;k<5;++k)
            {
                if(k==2 || k==3)
                {
                    values[k]=ExactExpansion::multiply(values[k],width);
                    derivatives[k]=ExactExpansion::multiply(derivatives[k],width);
                }
                values[k]=ExactExpansion::multiply(values[k],width);
            }
            apply(valueIndex,inputs,values,commonDenominator);
            apply(derivativeIndex,inputs,derivatives,commonDenominator);
        }
    }

    int independent(double value,Sensitivity::Kind kind,int index,int component)
    {
        Datum datum;
        if(restrictionSource_){Vector2D a,b;unsigned dx=0,dy=0;
            if(kind==Sensitivity::Kind::Vertex){a=b=vertices_.at(index).point;
                dx=component==1 || component==3?1:0;dy=component==2 || component==3?1:0;}
            else if(kind==Sensitivity::Kind::EdgeMidpoint){const auto& key=restrictionPatch_->edgeMidpoints.at(index).vertices;
                a=vertices_.at(key[0]).point;b=vertices_.at(key[1]).point;
                if(component){dx=a.y==b.y?0:1;dy=a.y==b.y?1:0;}}
            else {Detail::require(kind==Sensitivity::Kind::CellCenter,"unsupported exact refinement datum");
                a=cells_.at(index).lower;b=cells_.at(index).upper;}
            const NativeHermiteTransfer::DatumSampler<QuarticLeafData> sampler(*restrictionSource_);
            const auto exact=restrictionSource_->restrictionDatum(sampler.locate(a,b),dx,dy);
            datum.value.numerator[0]=exact.numerator[0];datum.value.denominator=exact.denominator;
            normalize(datum.value); // Same exact power-of-two gauge as constrained data.
        }else datum.value.numerator[0]=number(value);
        datum.identity={kind,index,component,0.0L};
        data_.push_back(std::move(datum));
        ++independentCount_;
        return static_cast<int>(data_.size()-1);
    }
    int linear(const std::array<int,5>& inputs,const std::array<Exact,5>& weights,const Exact& denominator)
    {
        Datum datum;
        datum.value=linearValue(inputs,weights,denominator,[&](std::size_t k) {
            datum.dependencies.push_back({inputs[k],ratio(weights[k],denominator)});
        });
        Detail::require(!datum.dependencies.empty(),"zero quartic constraint has no source identity");
        normalize(datum.value);
        data_.push_back(std::move(datum));
        return static_cast<int>(data_.size()-1);
    }
    template<class Record> Homogeneous<1> linearValue(const std::array<int,5>& inputs,
        const std::array<Exact,5>& weights,const Exact& denominator,const Record& record) const
    {
        const auto values=common(inputs);
        Homogeneous<1> result;
        result.denominator=ExactExpansion::multiply(values.denominator,denominator);
        for(std::size_t k=0;k<5;++k)if(ExactExpansion::sign(weights[k])!=0) {
            result.numerator[0]=ExactExpansion::add(result.numerator[0],ExactExpansion::multiply(weights[k],values.numerator[k]));
            record(k);
        }
        return result;
    }
    const EdgeData& skeletonEdge(const EdgeKey& key) const
    {
        const auto found=edgeIndex_.find(key);
        Detail::require(found!=edgeIndex_.end(),"quartic skeleton query requires an ordered active native side");
        return edges_[found->second];
    }
    template<std::size_t N> Homogeneous<N> common(const std::array<int,N>& ids) const
    {
        std::vector<Exact> denominators;
        std::array<std::size_t,N> positions;
        for(std::size_t k=0;k<N;++k)
        {
            Detail::require(ids[k]>=0 && static_cast<std::size_t>(ids[k])<data_.size(),"unresolved quartic datum");
            const auto& denominator=data_[ids[k]].value.denominator;
            const auto found=std::find_if(denominators.begin(),denominators.end(),[&](const Exact& other)
                { return ExactExpansion::sign(ExactExpansion::subtract(denominator,other))==0; });
            positions[k]=static_cast<std::size_t>(found-denominators.begin());
            if(found==denominators.end()) denominators.push_back(denominator);
        }
        Homogeneous<N> result;
        for(const auto& denominator:denominators)
            result.denominator=ExactExpansion::multiply(result.denominator,denominator);
        for(std::size_t k=0;k<N;++k)
        {
            auto factor=number(1);
            for(std::size_t j=0;j<denominators.size();++j) if(j!=positions[k])
                factor=ExactExpansion::multiply(factor,denominators[j]);
            result.numerator[k]=ExactExpansion::multiply(data_[ids[k]].value.numerator[0],factor);
        }
        return result;
    }
    static Exact number(double value) { return ExactExpansion::fromDouble(value); }
    static Exact difference(double a,double b) { return ExactExpansion::subtract(number(a),number(b)); }
    static bool equalRatios(const Exact& an,const Exact& ad,const Exact& bn,const Exact& bd)
    {
        if(an.components()==bn.components() && ad.components()==bd.components())return true;
        std::array<Exact,1> a{{an}},b{{bn}};auto da=ad,db=bd;
        // Equality needs cross-products, so balance complete homogeneous
        // objects first. No discarded tails or approximate division.
        ExactExpansion::centerHomogeneous(a,da);ExactExpansion::centerHomogeneous(b,db);
        return ExactExpansion::sign(ExactExpansion::subtract(
            ExactExpansion::multiply(a[0],db),ExactExpansion::multiply(b[0],da)))==0;
    }
    Homogeneous<1> restrictionDatum(const NativeHermiteTransfer::DatumSampler<QuarticLeafData>::Location& point,
        unsigned dx,unsigned dy) const
    {
        const auto slot=[](const RationalSegmentCoordinate& t,unsigned derivative){
            if(ExactExpansion::sign(t.numeratorBegin)==0)return derivative?2:0;
            if(ExactExpansion::sign(ExactExpansion::subtract(t.numeratorBegin,t.denominator))==0)return derivative?3:1;
            if(!derivative && ExactExpansion::sign(ExactExpansion::subtract(
                ExactExpansion::scalePow2(t.numeratorBegin,1),t.denominator))==0)return 4;
            return -1;};
        const int x=slot(point.x,dx),y=slot(point.y,dy);
        // Cardinal data already have their exact physical value. Preserve
        // this representation on unchanged leaves instead of repeatedly
        // introducing and clearing the same coordinate denominator.
        if(x>=0 && y>=0)return data_.at(activeCell(point.cell).physical[5*y+x]).value;
        const auto exact=restrictLocal(point.cell,point.x,point.y,dx,dy);
        Homogeneous<1> result;result.numerator[0]=exact.numerator[0];result.denominator=exact.denominator;return result;
    }
    Homogeneous<25> restrictionData(Vector2D lower,Vector2D upper) const
    {
        const NativeHermiteTransfer::DatumSampler<QuarticLeafData> sampler(*this);
        const auto& parent=activeCell(sampler.locate(lower,upper).cell);
        if(lower.x==parent.lower.x && lower.y==parent.lower.y && upper.x==parent.upper.x && upper.y==parent.upper.y)
            return parent.local;
        // One parent polynomial has ONE clearing denominator for all 25
        // data. Multiplying 25 separately cleared point denominators would
        // introduce redundant powers of its nonbinary physical widths.
        // The tensor restriction kernel and the shared constraints are the
        // existing ones; the factory checks their exact agreement afterwards.
        const auto hx=difference(parent.upper.x,parent.lower.x),hy=difference(parent.upper.y,parent.lower.y);
        const auto cx=difference(upper.x,lower.x),cy=difference(upper.y,lower.y);
        const auto dx=ExactExpansion::scalePow2(hx,1),dy=ExactExpansion::scalePow2(hy,1);
        const auto chart=[](int i,double lo,double hi,double oldLo,const Exact& denominator){
            const auto point=i==4?ExactExpansion::add(number(lo),number(hi)):
                ExactExpansion::scalePow2(number(i==1 || i==3?hi:lo),1);
            const auto n=ExactExpansion::subtract(point,ExactExpansion::scalePow2(number(oldLo),1));
            return RationalSegmentCoordinate{n,n,denominator};};
        Homogeneous<25> result;result.denominator=parent.local.denominator;
        for(int k=0;k<4;++k){result.denominator=ExactExpansion::multiply(result.denominator,dx);
            result.denominator=ExactExpansion::multiply(result.denominator,dy);}
        result.denominator=ExactExpansion::multiply(ExactExpansion::multiply(result.denominator,hx),hy);
        for(int y=0;y<5;++y)for(int x=0;x<5;++x){const bool derivativeX=x==2 || x==3,derivativeY=y==2 || y==3;
            const auto polynomial=QuarticCell::segmentPolynomialCleared(parent.local.numerator,
                chart(x,lower.x,upper.x,parent.lower.x,dx),chart(y,lower.y,upper.y,parent.lower.y,dy),derivativeX,derivativeY);
            result.numerator[5*y+x]=ExactExpansion::multiply(ExactExpansion::multiply(polynomial[0],derivativeX?cx:hx),derivativeY?cy:hy);}
        normalize(result);return result;
    }
    static double ratio(const Exact& numerator,const Exact& denominator)
    {
        Detail::require(ExactExpansion::sign(denominator)>0,"nonpositive quartic data denominator");
        const long double value=static_cast<long double>(ExactExpansion::estimate(numerator))
            /ExactExpansion::estimate(denominator);
        const auto result=Detail::checked(value,"unrepresentable quartic data ratio");
        Detail::require(ExactExpansion::sign(numerator)==0 || result!=0,"quartic data ratio underflow");
        return result;
    }
    double value(int id) const { return ratio(data_[id].value.numerator[0],data_[id].value.denominator); }
    template<std::size_t N> static void normalize(Homogeneous<N>& value)
    {
        Detail::require(ExactExpansion::sign(value.denominator)>0,"nonpositive quartic homogeneous denominator");
        const int exponent=-std::ilogb(value.denominator.components().back());
        value.denominator=ExactExpansion::scalePow2(value.denominator,exponent);
        for(auto& numerator:value.numerator) numerator=ExactExpansion::scalePow2(numerator,exponent);
    }
    static RationalSegmentCoordinate chart(double first,double last,double low,double high)
    {
        Detail::require(std::isfinite(first) && std::isfinite(last) && first>=low && first<=high
            && last>=low && last<=high,"quartic physical segment leaves its actual cell");
        RationalSegmentCoordinate result{difference(first,low),difference(last,low),difference(high,low)};
        const int exponent=-std::ilogb(result.denominator.components().back());
        result.numeratorBegin=ExactExpansion::scalePow2(result.numeratorBegin,exponent);
        result.numeratorEnd=ExactExpansion::scalePow2(result.numeratorEnd,exponent);
        result.denominator=ExactExpansion::scalePow2(result.denominator,exponent);
        return result;
    }
    static void add(long double& target,long double coefficient,long double weight)
    {
        const long double term=coefficient*weight;
        Detail::require(std::isfinite(term) && (coefficient==0 || weight==0 || term!=0),
            "quartic sensitivity chain-rule overflow or underflow");
        target+=term;
        Detail::require(std::isfinite(target),"nonfinite quartic sensitivity accumulation");
    }
    const Leaf& activeCell(int index) const
    {
        Detail::require(cellActive(index),"native quartic query/support crosses an inactive leaf");
        return cells_[index];
    }
};
} // namespace LevelSetHermiteScalar
} // namespace
} // namespace Quantum::Geometry::Detail
