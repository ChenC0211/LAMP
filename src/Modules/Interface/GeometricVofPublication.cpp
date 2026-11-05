#include "Modules/Interface/GeometricVofPublication.h"
#include "Core/IO/BinaryArchive.h"
#include "Core/Mesh/NativeLeafOverlap.h"
#include "Core/Geometry/Detail/FiniteVolumeClosure.h"
#include "Core/Mesh/CartesianMeasure.h"
#include "Core/Mesh/FaceSegmentGraph.h"
#include "Core/Mesh/NativeMeshIdentity.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <algorithm>
#include <atomic>
#include <cmath>
#include <cfenv>
#include <cerrno>
#include <cstring>
#if defined(_MSC_VER) && defined(_M_X64)
#include <float.h>
#if defined(_OPENMP)
#include <omp.h>
#endif
#endif
#include <exception>
#if defined(_M_X64) || defined(__SSE2__)
#include <xmmintrin.h>
#endif
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>

namespace LAMP::Interface
{
namespace
{
namespace Cut=Geometry::PlanarCellCut;
namespace FiniteVolume=Geometry::Detail::FiniteVolume;
using Point=Cut::Point;
using AnchoredPoint=Cut::AnchoredPoint;
using EdgeKey=std::pair<std::size_t,std::size_t>;
void require(bool ok,const std::string& message)
{ if(!ok) throw std::invalid_argument("[Interface/GeometricVofPublication] "+message); }
// Literal success checks need no temporary error string. The original string
// overload constructs the identical failure only after the condition fails.
void require(bool ok,const char* message)
{if(!ok)require(false,std::string(message));}
Point minus(Point a,Point b){ return {a.x-b.x,a.y-b.y}; }
AnchoredPoint coordinate(AnchoredPoint p,bool x)
{ return {{x?p.anchor.x:p.anchor.y,0},{x?p.offset.x:p.offset.y,0}}; }
Cut::LinearTriangleCut diagnosticCut(const Cut::AnchoredTriangleCut& actual)
{
    Cut::LinearTriangleCut result;result.kind=actual.kind;result.negativeCount=actual.negativeCount;
    result.zeroVertexMask=actual.zeroVertexMask;result.zeroEdgeMask=actual.zeroEdgeMask;result.unitNormal=actual.unitNormal;
    for(std::size_t i=0;i<actual.negativeCount;++i)result.negativeVertices[i]=Cut::approximateAnchoredPoint(actual.negativeVertices[i]);
    if(actual.interface)result.interface=std::array<Point,2>{Cut::approximateAnchoredPoint((*actual.interface)[0]),Cut::approximateAnchoredPoint((*actual.interface)[1])};
    return result;
}
std::array<AnchoredPoint,2> negativeLeftEdge(AnchoredPoint a,AnchoredPoint b,AnchoredPoint interior)
{
    // This exact source-side certificate stays in the existing strict Core TU.
    const auto cut=Cut::cutAnchoredTriangle({a,b,interior},{0,0,-1});
    require(static_cast<bool>(cut.interface),"zero-edge phase side unrepresentable");return *cut.interface;
}
EdgeKey edgeKey(std::size_t a,std::size_t b){ return std::minmax(a,b); }
bool mixed(double a){ return a>0 && a<1; }
struct PositiveSum
{
    long double sum=0,correction=0;
    void add(long double x){const auto y=x-correction,t=sum+y;correction=(t-sum)-y;sum=t;}
};
bool sameClock(GeometricVofTimeLayer a,GeometricVofTimeLayer b)
{ return a.step==b.step && a.begin==b.begin && a.end==b.end && a.position==b.position; }
struct RawCell
{
    Cut::Rectangle rectangle;
    Point center{},planeCenter{};
    long double width=0;
    Cut::Reconstruction plane;
    std::size_t centerId=0;
};
struct ZeroSide { int cell=-1; bool liquid=false; AnchoredPoint interior{}; };
// One independent mesh-only block, never a predecessor Publication/State.
// Public/raw construction is cold. Only a typed same-native derived producer
// may offer this already successful pure-integer incidence preparation.
struct LocalPlicNativeSupport
{
    std::shared_ptr<const NativeMeshIdentity> identity;
    bool radial=false;
    std::size_t cells=0,vertices=0,faces=0,bytes=0;
    std::vector<std::vector<std::size_t>> cellNodes;
    std::vector<std::vector<int>> nodeCells;
};
// Retained definition contains no borrowed tree/row references. Actual target
// support is rederived from the physical vertex and oriented wall legs.
struct PinnedCornerSupport
{
    Vector2D point;
    std::array<int,2> wetDirection{},dryDirection{};
    int vertex=-1,contact=-1,wet=-1,dry=-1,wetCell=-1,dryCell=-1;
    Cut::RectangleCorner corner=Cut::RectangleCorner::LowerLeft;
    bool negativeAlongHorizontal=true;
    GeometricVofPublication::PinCertificate::WallTraceLaw wallTraceLaw=
        GeometricVofPublication::PinCertificate::WallTraceLaw::StrictWholeCell;
    bool prescribedWet=false,prescribedDry=false;
    GeometricVofPublication::PinCertificate::Kind kind=GeometricVofPublication::PinCertificate::Kind::ResolvedCorner;
};
PinnedCornerSupport qualifyPinnedCorner(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,
    Vector2D point,int requestedWet,int requestedDry,const PinnedCornerSupport* retained=nullptr)
{
    require(std::isfinite(point.x)&&std::isfinite(point.y),"pinned corner point must be finite");
    PinnedCornerSupport pin;pin.point=point;
    if(retained)pin.wallTraceLaw=retained->wallTraceLaw; // New actual support, no old qualification flags.
    for(int v=0;v<mesh.vertIndexNum();++v){const auto p=mesh.verts(mesh.vertIndexToLabel[v]).coord();
        if(p.x==point.x&&p.y==point.y){require(pin.vertex<0,"ambiguous pinned native vertex");pin.vertex=v;}}
    require(pin.vertex>=0,"pin point must equal an actual native vertex");
    std::array<int,4> quadrant{{-1,-1,-1,-1}};int solid=-1,solidQuadrant=-1;unsigned active=0;
    for(int ci=0;ci<mesh.cellIndexNum();++ci){const auto& cell=mesh.cells(mesh.cellIndexToLabel[ci]);
        const auto lo=mesh.verts(cell.cellToVert[0]).coord(),hi=mesh.verts(cell.cellToVert[2]).coord();
        if(point.x<lo.x||point.x>hi.x||point.y<lo.y||point.y>hi.y)continue;
        require((point.x==lo.x||point.x==hi.x)&&(point.y==lo.y||point.y==hi.y),
            "pin requires actual contact-star rectangle corners, not a hanging point");
        const unsigned q=(point.x==lo.x?2u:0u)+(point.y==lo.y?1u:0u);
        require(quadrant[q]<0,"ambiguous overlapping pinned corner star");quadrant[q]=ci;
        if(domain.cellRowByIndex()[ci]>=0)++active;
        else {require(solid<0,"pin requires exactly one excluded solid corner leaf");solid=ci;solidQuadrant=static_cast<int>(q);}}
    require(active==3&&solid>=0&&std::all_of(quadrant.begin(),quadrant.end(),[](int ci){return ci>=0;}),
        "pin requires an internal three-fluid one-solid reentrant vertex");
    pin.contact=quadrant[static_cast<unsigned>(solidQuadrant)^3u];
    const auto& contact=mesh.cells(mesh.cellIndexToLabel[pin.contact]);
    bool contactCorner=false;for(unsigned k=0;k<4;++k)if(mesh.verts(contact.cellToVert[k]).index()==pin.vertex){
        pin.corner=static_cast<Cut::RectangleCorner>(k);contactCorner=true;}
    require(contactCorner,"pin is not an actual contact-cell corner");
    struct Leg {int row,cell;bool horizontal;std::array<int,2> direction;};std::vector<Leg> legs;
    for(int row=0;row<domain.nativeFaceCount();++row){const auto f=domain.face(mesh,row);
        const auto& edge=f.horizontal?mesh.hEdges(f.edgeLabel):mesh.vEdges(f.edgeLabel);
        const int a=mesh.verts(edge.edgeToVert[0]).index(),b=mesh.verts(edge.edgeToVert[1]).index();
        if(a!=pin.vertex&&b!=pin.vertex)continue;
        if((f.negativeCell>=0)==(f.positiveCell>=0))continue;
        require(edge.slaveCell>=0&&edge.masterCell>=0,"pin wall leg must be an internal selected/solid face");
        const int sa=mesh.cells(edge.slaveCell).index(),sb=mesh.cells(edge.masterCell).index();
        require(sa==solid||sb==solid,"ambiguous additional boundary at pin vertex");
        const int fluid=f.negativeCell>=0?f.negativeCell:f.positiveCell;
        require(fluid!=pin.contact,"wall leg cannot replace the opposite contact cell");
        const auto other=mesh.verts(mesh.vertIndexToLabel[a==pin.vertex?b:a]).coord();
        const std::array<int,2> direction{{other.x>point.x?1:other.x<point.x?-1:0,
            other.y>point.y?1:other.y<point.y?-1:0}};
        require((direction[0]==0)!=(direction[1]==0),"pin wall leg must be actual Cartesian geometry");
        legs.push_back({row,fluid,f.horizontal,direction});}
    require(legs.size()==2&&legs[0].horizontal!=legs[1].horizontal,"pin requires exactly two orthogonal actual wall legs");
    const auto choose=[&](bool wet)->const Leg&{
        for(const auto& leg:legs)if(retained?leg.direction==(wet?retained->wetDirection:retained->dryDirection):
            leg.row==(wet?requestedWet:requestedDry))return leg;
        require(false,"pinned wall row or retained physical leg does not match actual support");return legs.front();};
    const auto& wet=choose(true);const auto& dry=choose(false);
    require(wet.row!=dry.row&&wet.horizontal!=dry.horizontal,"wet/dry pin rows must be different orthogonal wall legs");
    pin.wet=wet.row;pin.dry=dry.row;pin.wetCell=wet.cell;pin.dryCell=dry.cell;
    pin.wetDirection=wet.direction;pin.dryDirection=dry.direction;
    pin.negativeAlongHorizontal=!wet.horizontal;
    return pin;
}
// Native ABI receipt only: no environment reset, flag union or changed control.
// Geometry/Plic results are copied only when their original call was quiet and
// the exact same complete qualified entry is present again.
struct LocalPlicQuietState
{
    std::fenv_t native{};unsigned csr=0;int error=0;bool valid=false,quiet=false;
    static LocalPlicQuietState read() noexcept;
    bool same(const LocalPlicQuietState&) const noexcept;
    bool matches() const noexcept;
    void completed(const LocalPlicQuietState&) noexcept;
};
struct LocalPlicGeometry
{
    std::shared_ptr<const NativeMeshIdentity> identity;
    bool radial=false;Point origin{};std::size_t bytes=0;
    std::vector<int> selectedRows,selectedCells;
    std::vector<Cut::Rectangle> rectangles;
    std::vector<Point> centers,vertices;
    std::vector<GeometricVofFace> faces;
    std::vector<std::vector<int>> neighbors,reflected;
    // Pure geometry only: unique equal-level whole-side neighbours for the
    // fixed inventory-height stencil. Missing/split/physical sides are -1.
    std::vector<std::array<int,4>> heightSides;
    LocalPlicQuietState vertexState,faceState;
};
struct BuildData
{
    std::optional<PinnedCornerSupport> pin;
    GeometricVofConfig localConfig;
    std::shared_ptr<const LocalPlicNativeSupport> localNativeSupport;
    // Ephemeral siblings: checkpoints never retain the mixed-plane reuse cache.
    std::shared_ptr<const LocalPlicGeometry> localGeometry;
    std::shared_ptr<std::vector<LocalPlicQuietState>> localMixedStates;
    GeometricVofIdentity identity;
    GeometricVofGeometryKind kind=GeometricVofGeometryKind::CompatibleFan;
    std::vector<GeometricVofPlicCell> plicCells;
    std::vector<FiniteVolume::Dyadic> capacityDefects;
    GeometricVofTimeLayer time;
    NativeLeafDomain domain;
    std::vector<double> alpha;
    Point origin{};
    std::vector<AnchoredPoint> anchoredVertices;
    std::vector<Point> vertices; // One-time diagnostic projections only.
    std::vector<long double> scalar;
    std::vector<GeometricVofCell> cells;
    std::vector<GeometricVofFace> faces;
    std::vector<GeometricVofFacet> facets;
    GeometricVofTopology topology;
    GeometricVofReport report;
};
// Optional integer storage must not leak allocation errno/FP state. This
// scope never encloses original capacities, planes, inventory or stock gates.
struct LocalPlicSupportEnvironment
{
    std::fenv_t native{};int error=errno;
#if defined(_M_X64) || defined(__SSE2__)
    unsigned csr=_mm_getcsr();
#endif
    bool captured=std::fegetenv(&native)==0;
    ~LocalPlicSupportEnvironment() noexcept {
        if(captured&&std::fesetenv(&native)!=0)std::terminate();
#if defined(_M_X64) || defined(__SSE2__)
        _mm_setcsr(csr);
#endif
        errno=error;
    }
};
#if defined(_MSC_VER) && !defined(__FAST_MATH__) && (!defined(_M_FP_FAST) || !_M_FP_FAST)
#pragma float_control(precise,on,push)
#pragma float_control(except,on)
#pragma fenv_access(on)
#endif
LocalPlicQuietState LocalPlicQuietState::read() noexcept
{
    LocalPlicQuietState out;out.error=errno;
#if defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE) && !defined(__FAST_MATH__) && (!defined(_M_FP_FAST) || !_M_FP_FAST) && !defined(__AVX__) && !defined(__AVX2__) && !defined(__AVX512F__) && !defined(__FMA__)
    if(sizeof(long double)!=8||std::numeric_limits<long double>::digits!=53
        ||out.error==EDOM||out.error==ERANGE||out.error==EINVAL)return out;
    unsigned control=0;out.csr=_mm_getcsr();
    const bool captured=std::fegetenv(&out.native)==0;
    const bool queried=_controlfp_s(&control,0,0)==0;errno=out.error;
    out.valid=captured&&queried&&std::fegetround()==FE_TONEAREST
        &&(control&_MCW_EM)==_MCW_EM&&(control&_MCW_RC)==_RC_NEAR
        &&(out.csr&0x1f80u)==0x1f80u&&(out.csr&0xe040u)==0;
#ifdef _EM_AMBIGUOUS
    if(control&_EM_AMBIGUOUS)out.valid=false;
#endif
#endif
    errno=out.error;return out;
}
bool LocalPlicQuietState::same(const LocalPlicQuietState& other) const noexcept
{return valid&&other.valid&&error==other.error&&csr==other.csr&&std::memcmp(&native,&other.native,sizeof native)==0;}
bool LocalPlicQuietState::matches() const noexcept
{return quiet&&same(read());}
void LocalPlicQuietState::completed(const LocalPlicQuietState& before) noexcept
{*this=before;quiet=before.same(read());}
#if defined(_MSC_VER) && !defined(__FAST_MATH__) && (!defined(_M_FP_FAST) || !_M_FP_FAST)
#if !defined(_M_FP_STRICT)
#pragma fenv_access(off)
#endif
#pragma float_control(pop)
#endif
template<class T> bool sameLocalPlicBits(const T& a,const T& b) noexcept
{return std::memcmp(&a,&b,sizeof(T))==0;}
std::size_t localPlicRecordBytes(const std::shared_ptr<std::vector<LocalPlicQuietState>>& records) noexcept
{return records?sizeof(std::vector<LocalPlicQuietState>)+512u+records->capacity()*sizeof(LocalPlicQuietState):0;}
void prepareLocalPlicRecords(BuildData& out,std::size_t count,const BuildData* predecessor) noexcept
{
    LocalPlicSupportEnvironment environment;if(!environment.captured||!LocalPlicQuietState::read().valid)return;
    try {
        constexpr std::size_t perArray=16u*1024u*1024u,limit=64u*1024u*1024u;
        if(count>perArray/sizeof(LocalPlicQuietState))return;
        const auto borrowed=predecessor?localPlicRecordBytes(predecessor->localMixedStates):0;
        const auto base=(out.localNativeSupport?out.localNativeSupport->bytes:0)
            +(out.localGeometry?out.localGeometry->bytes:0)+borrowed;
        const auto requested=sizeof(std::vector<LocalPlicQuietState>)+512u+count*sizeof(LocalPlicQuietState);
        if(base>limit||requested>limit-base)return;
        auto records=std::shared_ptr<std::vector<LocalPlicQuietState>>(new std::vector<LocalPlicQuietState>,
            [](std::vector<LocalPlicQuietState>* p) noexcept {LocalPlicSupportEnvironment released;delete p;});
        records->resize(count);
        if(localPlicRecordBytes(records)>limit-base)return;
        out.localMixedStates=std::move(records);
    }catch(...){out.localMixedStates.reset();}
}
const LocalPlicGeometry* qualifiedLocalPlicGeometry(const QuadtreeMesh& mesh,const BuildData& out) noexcept
{
#if defined(_MSC_VER) && defined(_M_X64)
    if(sizeof(long double)!=8||std::numeric_limits<long double>::digits!=53)return nullptr;
    const auto* g=out.localGeometry.get();const auto n=static_cast<std::size_t>(mesh.cellIndexNum());
    if(g&&g->identity==out.domain.meshIdentity()&&g->radial==out.domain.axisymmetric()
        &&sameLocalPlicBits(g->origin.x,out.origin.x)&&sameLocalPlicBits(g->origin.y,out.origin.y)
        &&g->selectedRows==out.domain.cellRowByIndex()&&g->selectedCells==out.domain.cellIndexByRow()
        &&g->rectangles.size()==n&&g->centers.size()==n
        &&g->neighbors.size()==n&&g->reflected.size()==n&&(g->radial||g->heightSides.size()==n)
        &&g->vertices.size()==static_cast<std::size_t>(mesh.vertIndexNum())
        &&g->faces.size()==static_cast<std::size_t>(out.domain.nativeFaceCount()))return g;
#else
    (void)mesh;(void)out;
#endif
    return nullptr;
}

void retainLocalPlicNativeSupport(BuildData& out,
    std::vector<std::vector<std::size_t>>& cellNodes,std::vector<std::vector<int>>& nodeCells) noexcept
{
    LocalPlicSupportEnvironment environment;if(!environment.captured)return;
    try {
        constexpr std::size_t limit=32u*1024u*1024u;
        std::size_t bytes=sizeof(LocalPlicNativeSupport)+512u;
        const auto account=[&](std::size_t count,std::size_t width){
            if(bytes>limit||(width&&count>(limit-bytes)/width))return false;
            bytes+=count*width;return true;
        };
        if(!account(cellNodes.capacity(),sizeof(std::vector<std::size_t>))
            ||!account(nodeCells.capacity(),sizeof(std::vector<int>)))return;
        for(const auto& list:cellNodes)if(!account(list.capacity(),sizeof(std::size_t))
            ||(!list.empty()&&!account(1,32u)))return;
        for(const auto& list:nodeCells)if(!account(list.capacity(),sizeof(int))
            ||(!list.empty()&&!account(1,32u)))return;
        auto support=std::shared_ptr<LocalPlicNativeSupport>(new LocalPlicNativeSupport,
            [](LocalPlicNativeSupport* value) noexcept {
                LocalPlicSupportEnvironment released;
                delete value;
            });
        support->identity=out.domain.meshIdentity();support->radial=out.domain.axisymmetric();
        support->cells=cellNodes.size();support->vertices=nodeCells.size();
        support->faces=static_cast<std::size_t>(out.domain.nativeFaceCount());support->bytes=bytes;
        // Both arrays have passed every original incidence/cardinality check.
        // Their final numeric consumer has completed; these swaps allocate
        // nothing and cannot leave a partially retained source preparation.
        support->cellNodes.swap(cellNodes);support->nodeCells.swap(nodeCells);
        out.localNativeSupport=std::move(support);
    }catch(...){/* Optional retention only; all physical work already passed. */}
}
struct Builder
{
    const QuadtreeMesh& mesh;
    BuildData& out;
    GeometricVofConfig config;
    std::vector<RawCell> raw;
    std::vector<std::vector<int>> incidence;
    std::map<EdgeKey,std::size_t> roots;
    std::map<EdgeKey,std::vector<ZeroSide>> zeroSides;
    std::map<EdgeKey,int> nativeRows;
    bool rz;
    std::size_t addPoint(AnchoredPoint p,long double scalar)
    {
        require(std::isfinite(p.anchor.x)&&std::isfinite(p.anchor.y)&&std::isfinite(p.offset.x)&&std::isfinite(p.offset.y)&&std::isfinite(scalar),"nonfinite derived site");
        const auto id=out.anchoredVertices.size();out.anchoredVertices.push_back(p);
        out.vertices.push_back(Cut::approximateAnchoredPoint(p));out.scalar.push_back(scalar);return id;
    }
    std::size_t addPoint(Point p,long double scalar){return addPoint(AnchoredPoint{p,{}},scalar);}
    long double planeValue(int cell,Point p) const
    {
        const auto& r=raw[cell]; const auto d=minus(p,r.planeCenter);
        return r.plane.normal.x*d.x+r.plane.normal.y*d.y-r.plane.offset;
    }
    long double planeValue(int cell,AnchoredPoint p) const
    {
        // Approximate reconstruction seed, not an authoritative plane predicate.
        // Native offset-zero sites retain the original floating operation path.
        if(p.offset.x==0&&p.offset.y==0)return planeValue(cell,p.anchor);
        const auto& r=raw[cell];const auto d=Cut::anchoredDifference(p,{r.planeCenter,{}});
        return r.plane.normal.x*d.x+r.plane.normal.y*d.y-r.plane.offset;
    }
    long double radius(AnchoredPoint p) const
    {
        const auto value=Cut::anchoredDifference(coordinate(p,true),{{-out.origin.x,0},{}}).x;
        require(std::isfinite(value)&&value>=0,"negative/nonfinite actual radius");return value;
    }
    bool onAxis(AnchoredPoint p) const
    {return Cut::sameAnchoredPoint(coordinate(p,true),{{-out.origin.x,0},{}});}
    void layout()
    {
        const int n=mesh.cellIndexNum();
        out.cells.resize(n);raw.resize(n);incidence.resize(mesh.vertIndexNum());
        out.anchoredVertices.reserve(static_cast<std::size_t>(mesh.vertIndexNum())+n);
        out.vertices.reserve(out.anchoredVertices.capacity());
        out.scalar.reserve(out.anchoredVertices.capacity());
        for(int vi=0;vi<mesh.vertIndexNum();++vi){
            const auto p=mesh.verts(mesh.vertIndexToLabel[vi]).coord();
            addPoint(Point{static_cast<long double>(p.x)-out.origin.x,static_cast<long double>(p.y)-out.origin.y},0);
        }
        FaceSegmentGraph graph;graph.build(mesh);
        for(int ci=0;ci<n;++ci){
            auto& cell=out.cells[ci];const auto& native=mesh.cells(mesh.cellIndexToLabel[ci]);
            const auto lower=mesh.verts(native.cellToVert[0]).coord(),upper=mesh.verts(native.cellToVert[2]).coord();
            auto& r=raw[ci];r.rectangle={lower,upper};
            const auto cc=native.centerCoord();r.center={static_cast<long double>(cc.x)-out.origin.x,static_cast<long double>(cc.y)-out.origin.y};
            r.planeCenter={static_cast<long double>(lower.x)+(static_cast<long double>(upper.x)-lower.x)/2-out.origin.x,
                static_cast<long double>(lower.y)+(static_cast<long double>(upper.y)-lower.y)/2-out.origin.y};
            r.width=std::min(static_cast<long double>(upper.x)-lower.x,static_cast<long double>(upper.y)-lower.y);
            cell.nativeCell=ci;cell.capacity=Cut::capacity(r.rectangle,rz);cell.inventory=out.alpha[ci]*cell.capacity;
            const auto row=out.domain.cellRowByIndex()[ci];
            require(row>=0,"selected domain omits native cell="+std::to_string(ci));
            const auto nativeVolume=out.domain.cellVolumes()[row];
            const auto capacityDefect=std::abs(cell.capacity-nativeVolume)/cell.capacity;
            require(capacityDefect<=config.relativeVolumeTolerance,"Core/native capacity mismatch cell="+std::to_string(ci));
            out.report.maximumRelativeCapacityResidual=std::max(out.report.maximumRelativeCapacityResidual,capacityDefect);
            std::map<std::size_t,std::size_t> next;
            const auto range=graph.segmentRangeForCellIndex(ci);
            for(auto si=range.first;si<range.second;++si){
                const auto& s=graph.segments()[si];const auto& e=s.isHorizontal?mesh.hEdges(s.edgeLabel):mesh.vEdges(s.edgeLabel);
                auto a=static_cast<std::size_t>(mesh.verts(e.edgeToVert[0]).index());
                auto b=static_cast<std::size_t>(mesh.verts(e.edgeToVert[1]).index());
                if(s.side==FaceSide::North || s.side==FaceSide::West)std::swap(a,b);
                require(a<incidence.size()&&b<incidence.size()&&next.emplace(a,b).second,"invalid leaf boundary incidence");
            }
            require(next.size()>=4&&next.size()<=8,"leaf fan requires 4--8 real 2:1 subfaces");
            const auto start=static_cast<std::size_t>(mesh.verts(native.cellToVert[0]).index());auto id=start;
            do {
                require(cell.boundaryVertices.size()<next.size(),"leaf boundary is not a single cycle");
                cell.boundaryVertices.push_back(id);incidence[id].push_back(ci);
                const auto found=next.find(id);require(found!=next.end(),"broken real subface boundary");id=found->second;
            } while(id!=start);
            require(cell.boundaryVertices.size()==next.size(),"disconnected leaf boundary");
            r.centerId=addPoint(r.center,0);
        }
        for(auto& cells:incidence){std::sort(cells.begin(),cells.end());cells.erase(std::unique(cells.begin(),cells.end()),cells.end());}
        for(int row=0;row<out.domain.nativeFaceCount();++row){
            const auto f=out.domain.face(mesh,row);const auto& e=f.horizontal?mesh.hEdges(f.edgeLabel):mesh.vEdges(f.edgeLabel);
            const auto a=static_cast<std::size_t>(mesh.verts(e.edgeToVert[0]).index());
            const auto b=static_cast<std::size_t>(mesh.verts(e.edgeToVert[1]).index());
            require(nativeRows.emplace(edgeKey(a,b),row).second,"duplicate real shared face");
        }
    }
    void planes()
    {
        for(std::size_t ci=0;ci<out.cells.size();++ci)if(mixed(out.alpha[ci])){
            std::set<int> neighbors;
            for(auto v:out.cells[ci].boundaryVertices)for(int k:incidence[v])if(k!=static_cast<int>(ci))neighbors.insert(k);
            // One physical node ring; no uniform-grid derivative or supplied normal.
            long double xx=0,xy=0,yy=0,bx=0,by=0;
            const auto h=raw[ci].width;
            for(int k:neighbors){
                const auto d=minus(raw[k].center,raw[ci].center);const auto x=d.x/h,y=d.y/h;
                const auto d2=x*x+y*y;require(d2>0&&std::isfinite(d2),"invalid alpha normal neighbor");
                const auto w=1/d2,da=static_cast<long double>(out.alpha[k])-out.alpha[ci];
                xx+=w*x*x;xy+=w*x*y;yy+=w*y*y;bx+=w*x*da;by+=w*y*da;
            }
            const auto determinant=xx*yy-xy*xy;
            require(std::isfinite(determinant)&&determinant>32*std::numeric_limits<long double>::epsilon()*(xx+yy)*(xx+yy),
                "rank-deficient alpha normal support cell="+std::to_string(ci));
            const auto gx=(yy*bx-xy*by)/determinant,gy=(xx*by-xy*bx)/determinant;
            const auto norm=std::hypot(gx,gy);
            require(std::isfinite(norm)&&norm>0,"zero alpha gradient; interface orientation unresolved cell="+std::to_string(ci));
            raw[ci].plane=Cut::reconstruct(raw[ci].rectangle,Vector2D(static_cast<double>(-gx/norm),static_cast<double>(-gy/norm)),out.alpha[ci],rz);
        }
    }
    void sharedScalars()
    {
        for(std::size_t v=0;v<incidence.size();++v){
            bool liquid=false,gas=false;long double total=0,weight=0,width=std::numeric_limits<long double>::max();
            for(int ci:incidence[v]){
                width=std::min(width,raw[ci].width);
                liquid|=out.alpha[ci]==1;gas|=out.alpha[ci]==0;
                if(mixed(out.alpha[ci])){const auto w=1/(raw[ci].width*raw[ci].width);total+=w*planeValue(ci,out.anchoredVertices[v]);weight+=w;}
            }
            require(!incidence[v].empty(),"native vertex lacks real leaf incidence");
            long double value=weight>0?total/weight:(liquid?-width:width);
            if(liquid&&gas)value=0;
            else if(liquid)value=-std::abs(value);
            else if(gas)value=std::abs(value);
            out.scalar[v]=value;
        }
    }
    GeometricVofTriangle triangle(std::size_t a,std::size_t b,std::size_t c) const
    {
        GeometricVofTriangle t;t.vertexIds={a,b,c};
        for(unsigned j=0;j<3;++j){t.anchoredVertices[j]=out.anchoredVertices[t.vertexIds[j]];t.vertices[j]=Cut::approximateAnchoredPoint(t.anchoredVertices[j]);t.scalar[j]=out.scalar[t.vertexIds[j]];}
        const auto edge=negativeLeftEdge(t.anchoredVertices[0],t.anchoredVertices[1],t.anchoredVertices[2]);
        require(Cut::sameAnchoredPoint(edge[0],t.anchoredVertices[0])&&Cut::sameAnchoredPoint(edge[1],t.anchoredVertices[1]),"non-CCW fan triangle");return t;
    }
    long double measured(const GeometricVofCell& cell,std::size_t center,long double value) const
    {
        PositiveSum volume;
        for(const auto& t:cell.triangles){auto s=t.scalar;for(unsigned j=0;j<3;++j)if(t.vertexIds[j]==center)s[j]=value;
                        try {volume.add(Cut::negativeAnchoredMeasure(t.anchoredVertices,s,out.origin,rz));}
            catch(const std::exception& error){std::ostringstream context;context<<"fan volume query cell="<<cell.nativeCell<<std::hexfloat<<" centre="<<value<<" scalars=("<<s[0]<<","<<s[1]<<","<<s[2]<<"): "<<error.what();throw std::runtime_error("[Interface/GeometricVofPublication] "+context.str());}}
        require(std::isfinite(volume.sum),"nonfinite fan volume");return volume.sum;
    }
    void reconstruct()
    {
        for(std::size_t ci=0;ci<out.cells.size();++ci){
            auto& cell=out.cells[ci];auto& r=raw[ci];const auto a=out.alpha[ci];
            if(!mixed(a)){
                cell.centerScalar=a==1?-r.width:r.width;out.scalar[r.centerId]=cell.centerScalar;
                cell.measured=a==1?cell.capacity:0;
                if(a==1)++out.report.pureLiquidCells;else++out.report.pureGasCells;
                for(auto v:cell.boundaryVertices)require(a==1?out.scalar[v]<=0:out.scalar[v]>=0,"pure phase node constraint violated");
            }else{
                ++out.report.mixedCells;long double collarVolume=0;
                for(std::size_t j=0;j<cell.boundaryVertices.size();++j){
                    const auto va=cell.boundaryVertices[j],vb=cell.boundaryVertices[(j+1)%cell.boundaryVertices.size()];
                    if(out.scalar[va]==0&&out.scalar[vb]==0){
                        const auto p=out.anchoredVertices[va].anchor,q=out.anchoredVertices[vb].anchor;
                        const Point midpoint{p.x+(q.x-p.x)/2,p.y+(q.y-p.y)/2};
                        const auto fraction=std::min(static_cast<long double>(a),1-static_cast<long double>(a))/12;
                        const auto d=Cut::anchoredDifference(out.anchoredVertices[r.centerId],{midpoint,{}});
                        const AnchoredPoint support{midpoint,{fraction*d.x,fraction*d.y}};
                        require(!Cut::sameAnchoredPoint(support,{midpoint,{}}),"collapsed fixed zero-edge support cell="+std::to_string(ci));
                        auto value=planeValue(static_cast<int>(ci),support);
                        if(value==0)value=(a>=.5?-1.L:1.L)*r.width*fraction;
                        require(value!=0&&std::isfinite(value),"unrepresentable fixed zero-edge support cell="+std::to_string(ci));
                        const auto vp=addPoint(support,value);
                        auto outer=triangle(va,vb,vp);
                        const std::vector<AnchoredPoint> outerPolygon(outer.anchoredVertices.begin(),outer.anchoredVertices.end());
                        require(Cut::anchoredPolygonInsideBox(outerPolygon,out.origin,r.rectangle),"fixed zero-edge support outside actual leaf");
                        collarVolume+=Cut::anchoredPolygonMeasureInFrame(outerPolygon,out.origin,rz);
                        cell.triangles.push_back(outer);cell.triangles.push_back(triangle(va,vp,r.centerId));cell.triangles.push_back(triangle(vp,vb,r.centerId));
                        ++out.report.zeroCollars;
                    }else cell.triangles.push_back(triangle(va,vb,r.centerId));
                }
                require(collarVolume<std::min(cell.inventory,cell.capacity-cell.inventory)/4,
                    "zero-edge collar exceeds actual minor phase capacity cell="+std::to_string(ci));
                long double span=r.width;
                for(const auto& t:cell.triangles)for(unsigned j=0;j<3;++j)if(t.vertexIds[j]!=r.centerId)span=std::max(span,std::abs(t.scalar[j]));
                long double lo=-span,hi=span,vl=measured(cell,r.centerId,lo),vh=measured(cell,r.centerId,hi);
                unsigned expansions=0;
                while((vl<cell.inventory||vh>cell.inventory)&&expansions<config.maximumBracketExpansions){
                    if(vl<cell.inventory){lo*=2;require(std::isfinite(lo),"negative volume bracket overflow");vl=measured(cell,r.centerId,lo);}
                    if(vh>cell.inventory){hi*=2;require(std::isfinite(hi),"positive volume bracket overflow");vh=measured(cell,r.centerId,hi);}++expansions;
                }
                require(vl>=cell.inventory&&vh<=cell.inventory,"finite fan volume bracket unavailable cell="+std::to_string(ci));
                long double best=std::abs(vl-cell.inventory)<std::abs(vh-cell.inventory)?lo:hi;
                auto actual=measured(cell,r.centerId,best);unsigned iterations=0;
                const auto tolerance=std::min(config.relativeVolumeTolerance*cell.capacity,
                    config.relativeVolumeTolerance*std::min(cell.inventory,cell.capacity-cell.inventory));
                require(std::isfinite(tolerance)&&tolerance>0,"minor-phase volume tolerance is unrepresentable cell="+std::to_string(ci));
                const auto resolved=[&]{return actual>0&&actual<cell.capacity&&std::abs(actual-cell.inventory)<=tolerance;};
                for(unsigned trial=0;!resolved()&&trial<config.maximumVolumeIterations;++trial){
                    const auto mid=lo+(hi-lo)/2;if(mid==lo||mid==hi)break;
                    ++iterations;const auto vm=measured(cell,r.centerId,mid);
                    if(std::abs(vm-cell.inventory)<std::abs(actual-cell.inventory)){best=mid;actual=vm;}
                    if(vm==cell.inventory)break;
                    if(vm>cell.inventory)lo=mid;else hi=mid;
                    // Stop at the stricter minor-phase residual once an actual representable two-phase query passes.
                }
                if(!resolved()){
                    std::ostringstream detail;detail.precision(std::numeric_limits<long double>::max_digits10);
                    detail<<"fan inventory residual or minor phase unresolved cell="<<ci<<" target="<<cell.inventory
                        <<" actual="<<actual<<" tolerance="<<tolerance<<" iterations="<<iterations
                        <<" centre="<<best<<" lo="<<lo<<" vlo="<<measured(cell,r.centerId,lo)
                        <<" hi="<<hi<<" vhi="<<measured(cell,r.centerId,hi)<<" frame="<<out.origin.x<<","<<out.origin.y;
                    for(const auto& t:cell.triangles){detail<<" triangle=";
                        for(unsigned j=0;j<3;++j){const auto p=t.anchoredVertices[j];
                            detail<<"(anchor="<<p.anchor.x<<","<<p.anchor.y<<" offset="<<p.offset.x<<","<<p.offset.y
                                <<" scalar="<<(t.vertexIds[j]==r.centerId?best:t.scalar[j])<<" id="<<t.vertexIds[j]<<")";}}
                    require(false,detail.str());}
                cell.centerScalar=best;out.scalar[r.centerId]=best;cell.measured=actual;
                out.report.maximumVolumeIterations=std::max(out.report.maximumVolumeIterations,iterations);
                for(auto& t:cell.triangles){
                    for(unsigned j=0;j<3;++j)if(t.vertexIds[j]==r.centerId)t.scalar[j]=best;
                    t.anchoredCut=Cut::cutAnchoredTriangle(t.anchoredVertices,t.scalar);t.cut=diagnosticCut(t.anchoredCut);
                    require(t.anchoredCut.kind!=Cut::TriangleRegion::AllZeroScalar,"positive-area all-zero fan scalar cell="+std::to_string(ci));
                }
                PositiveSum published;
                for(const auto& t:cell.triangles)if(t.anchoredCut.negativeCount>=3)
                    published.add(Cut::anchoredPolygonMeasureInFrame({t.anchoredCut.negativeVertices.begin(),t.anchoredCut.negativeVertices.begin()+t.anchoredCut.negativeCount},out.origin,rz));
                require(published.sum>0&&published.sum<cell.capacity&&std::abs(published.sum-cell.inventory)<=tolerance,
                    "published fan inventory residual or minor phase unresolved cell="+std::to_string(ci));
                cell.measured=published.sum;
            }
            const auto residual=std::abs(cell.measured-cell.inventory)/cell.capacity;
            out.report.maximumRelativeVolumeResidual=std::max(out.report.maximumRelativeVolumeResidual,residual);
            out.report.inventory+=cell.inventory;out.report.measured+=cell.measured;
        }
    }
    std::size_t zeroRoot(std::size_t a,std::size_t b)
    {
        if(out.scalar[a]==0)return a;if(out.scalar[b]==0)return b;
        const auto key=edgeKey(a,b);const auto found=roots.find(key);if(found!=roots.end())return found->second;
        const auto p=Cut::anchoredZeroCrossing(out.anchoredVertices[key.first],out.anchoredVertices[key.second],out.scalar[key.first],out.scalar[key.second]);
        const auto id=addPoint(p,0);roots.emplace(key,id);return id;
    }
    void facet(std::size_t a,std::size_t b,int liquid,int gas,int row,GeometricVofFacet::Kind kind)
    {
        // Callers supply the exact actual-site negative-left order.
        const auto pa=out.anchoredVertices[a],pb=out.anchoredVertices[b];const auto d=Cut::anchoredDifference(pb,pa);
        const auto length=std::hypot(d.x,d.y);require(std::isfinite(length)&&length>0,"unrepresentable actual Gamma segment");
        const Point outward{d.y/length,-d.x/length};
        const bool axis=rz&&onAxis(pa)&&onAxis(pb);
        const auto meanRadius=rz?radius(pa)/2+radius(pb)/2:0;
        const auto measure=rz?length*CartesianMeasure::twoPi*meanRadius:length;
        require(std::isfinite(measure)&&(axis?measure==0:measure>0),"unrepresentable physical Gamma measure");
        if(axis)return; // Exact actual radial-axis sites are domain boundary, never charged Gamma.
        out.facets.push_back({kind,{a,b},outward,liquid,gas,row,length,measure});
    }
    void gamma()
    {
        for(const auto& cell:out.cells){
            if(cell.triangles.empty()){
                for(std::size_t j=0;j<cell.boundaryVertices.size();++j){const auto a=cell.boundaryVertices[j],b=cell.boundaryVertices[(j+1)%cell.boundaryVertices.size()];
                    if(out.scalar[a]==0&&out.scalar[b]==0)zeroSides[edgeKey(a,b)].push_back({cell.nativeCell,out.alpha[cell.nativeCell]==1,{raw[cell.nativeCell].center,{}}});}
                continue;
            }
            for(const auto& t:cell.triangles){
                for(unsigned j=0;j<3;++j)if(t.scalar[j]==0&&t.scalar[(j+1)%3]==0){
                    const auto opposite=(j+2)%3;require(t.scalar[opposite]!=0,"ambiguous zero-edge phase");
                    zeroSides[edgeKey(t.vertexIds[j],t.vertexIds[(j+1)%3])].push_back({cell.nativeCell,t.scalar[opposite]<0,t.anchoredVertices[opposite]});
                }
                if(!t.anchoredCut.interface||t.anchoredCut.zeroEdgeMask)continue;
                std::vector<std::size_t> ends;
                for(unsigned j=0;j<3;++j){
                    if(t.scalar[j]==0)ends.push_back(t.vertexIds[j]);
                    const auto k=(j+1)%3;
                    if((t.scalar[j]<0&&t.scalar[k]>0)||(t.scalar[j]>0&&t.scalar[k]<0))ends.push_back(zeroRoot(t.vertexIds[j],t.vertexIds[k]));
                }
                std::sort(ends.begin(),ends.end());ends.erase(std::unique(ends.begin(),ends.end()),ends.end());
                require(ends.size()==2,"actual noncoincident triangle Gamma lacks two canonical roots");
                if(Cut::sameAnchoredPoint(out.anchoredVertices[ends[1]],(*t.anchoredCut.interface)[0]))std::swap(ends[0],ends[1]);
                require(Cut::sameAnchoredPoint(out.anchoredVertices[ends[0]],(*t.anchoredCut.interface)[0])&&
                    Cut::sameAnchoredPoint(out.anchoredVertices[ends[1]],(*t.anchoredCut.interface)[1]),"canonical Gamma roots disagree with actual fan cut");
                facet(ends[0],ends[1],cell.nativeCell,cell.nativeCell,-1,GeometricVofFacet::Kind::CellFan);
            }
        }
        for(const auto& pair:zeroSides){
            const auto& sides=pair.second;require(sides.size()<=2,"more than two actual phases incident to one fan edge");
            if(sides.size()!=2||sides[0].liquid==sides[1].liquid)continue;
            const auto& liquid=sides[sides[0].liquid?0:1];const auto& gas=sides[sides[0].liquid?1:0];
            auto a=pair.first.first,b=pair.first.second;
            const auto edge=negativeLeftEdge(out.anchoredVertices[a],out.anchoredVertices[b],liquid.interior);
            if(Cut::sameAnchoredPoint(edge[0],out.anchoredVertices[b]))std::swap(a,b);
            require(Cut::sameAnchoredPoint(edge[0],out.anchoredVertices[a])&&Cut::sameAnchoredPoint(edge[1],out.anchoredVertices[b]),"zero-edge phase side unrepresentable");
            const auto native=nativeRows.find(pair.first);const int row=native==nativeRows.end()?-1:native->second;
            facet(a,b,liquid.cell,gas.cell,row,row<0?GeometricVofFacet::Kind::CellFan:GeometricVofFacet::Kind::NativeCoincident);
        }
    }
    bool zeroFacePhase(EdgeKey key,int cell) const
    {
        const auto found=zeroSides.find(key);require(found!=zeroSides.end(),"zero face missing same-owner phase incidence");
        for(const auto& side:found->second)if(side.cell==cell)return side.liquid;
        require(false,"zero face lacks actual adjacent phase");return false;
    }
    long double rangeMeasure(std::size_t a,std::size_t b) const
    {
        const auto pa=out.anchoredVertices[a],pb=out.anchoredVertices[b];const auto d=Cut::anchoredDifference(pb,pa);
        const auto length=std::hypot(d.x,d.y);require(std::isfinite(length)&&length>0,"unrepresentable actual face trace");
        if(!rz)return length;
        const bool axis=onAxis(pa)&&onAxis(pb);const auto meanRadius=radius(pa)/2+radius(pb)/2;
        const auto measure=length*CartesianMeasure::twoPi*meanRadius;
        require(std::isfinite(measure)&&(axis?measure==0:measure>0),"unrepresentable physical face trace");return measure;
    }
    void faceTraces()
    {
        out.faces.reserve(out.domain.nativeFaceCount());
        for(int row=0;row<out.domain.nativeFaceCount();++row){
            const auto native=out.domain.face(mesh,row);const auto& edge=native.horizontal?mesh.hEdges(native.edgeLabel):mesh.vEdges(native.edgeLabel);
            GeometricVofFace f;f.nativeRow=row;f.negativeCell=native.negativeCell;f.positiveCell=native.positiveCell;
            f.vertices={static_cast<std::size_t>(mesh.verts(edge.edgeToVert[0]).index()),static_cast<std::size_t>(mesh.verts(edge.edgeToVert[1]).index())};
            const auto a=f.vertices[0],b=f.vertices[1];const auto va=out.scalar[a],vb=out.scalar[b];
            if(va==0&&vb==0){
                const bool minus=native.negativeCell>=0&&zeroFacePhase(edgeKey(a,b),native.negativeCell);
                const bool plus=native.positiveCell>=0&&zeroFacePhase(edgeKey(a,b),native.positiveCell);
                f.ranges[0]={0,1,minus,plus,{a,b}};f.rangeCount=1;
                f.coincident=native.negativeCell>=0&&native.positiveCell>=0&&minus!=plus;
                if(f.coincident)++out.report.coincidentFaces;
            }else if((va<0&&vb>0)||(va>0&&vb<0)){
                const auto root=zeroRoot(a,b);const auto pa=out.anchoredVertices[a],pb=out.anchoredVertices[b],pr=out.anchoredVertices[root];
                require(Cut::sameAnchoredPoint(coordinate(pa,!native.horizontal),coordinate(pr,!native.horizontal))&&
                    Cut::sameAnchoredPoint(coordinate(pb,!native.horizontal),coordinate(pr,!native.horizontal)),"shared face root off actual native segment");
                const auto firstChord=Cut::anchoredDifference(pr,pa),secondChord=Cut::anchoredDifference(pb,pr),wholeChord=Cut::anchoredDifference(pb,pa);
                const auto firstLength=native.horizontal?firstChord.x:firstChord.y,secondLength=native.horizontal?secondChord.x:secondChord.y;
                const auto wholeLength=native.horizontal?wholeChord.x:wholeChord.y;
                require(firstLength>0&&secondLength>0&&wholeLength>0,"shared face root collapsed onto native vertex");
                const auto parameter=firstLength/wholeLength; // Diagnostic only; may round to an endpoint.
                const bool first=va<0,second=vb<0;
                f.ranges[0]={0,parameter,native.negativeCell>=0&&first,native.positiveCell>=0&&first,{a,root}};
                f.ranges[1]={parameter,1,native.negativeCell>=0&&second,native.positiveCell>=0&&second,{root,b}};f.rangeCount=2;
            }else{
                const bool liquid=va<0||vb<0;f.ranges[0]={0,1,native.negativeCell>=0&&liquid,native.positiveCell>=0&&liquid,{a,b}};f.rangeCount=1;
            }
            for(unsigned j=0;j<f.rangeCount;++j){const auto& range=f.ranges[j];const auto measure=rangeMeasure(range.vertices[0],range.vertices[1]);
                if(range.minusLiquid)f.liquidMinus+=measure;if(range.plusLiquid)f.liquidPlus+=measure;}
            out.faces.push_back(f);
        }
    }
    void connections()
    {
        std::vector<std::vector<std::size_t>> adjacent(out.anchoredVertices.size());
        for(std::size_t i=0;i<out.facets.size();++i)for(auto v:out.facets[i].vertices)adjacent[v].push_back(i);
        const auto xmin=static_cast<long double>(mesh.minCoord.x)-out.origin.x,ymin=static_cast<long double>(mesh.minCoord.y)-out.origin.y;
        const auto xmax=static_cast<long double>(mesh.maxCoord.x)-out.origin.x,ymax=static_cast<long double>(mesh.maxCoord.y)-out.origin.y;
        for(std::size_t v=0;v<adjacent.size();++v)if(!adjacent[v].empty()&&adjacent[v].size()!=2){
            const auto p=out.anchoredVertices[v];
            if(adjacent[v].size()==1&&rz&&onAxis(p))out.topology.axisEndpoints.push_back(v);
            else if(adjacent[v].size()==1&&(Cut::sameAnchoredPoint(coordinate(p,true),{{xmin,0},{}})||
                Cut::sameAnchoredPoint(coordinate(p,false),{{ymin,0},{}})||Cut::sameAnchoredPoint(coordinate(p,true),{{xmax,0},{}})||
                Cut::sameAnchoredPoint(coordinate(p,false),{{ymax,0},{}})))out.topology.exteriorEndpoints.push_back(v);
            else out.topology.junctionVertices.push_back(v);
        }
        out.topology.facetComponent.assign(out.facets.size(),std::numeric_limits<std::size_t>::max());
        for(std::size_t i=0;i<out.facets.size();++i)if(out.topology.facetComponent[i]==std::numeric_limits<std::size_t>::max()){
            const auto component=out.topology.components++;std::vector<std::size_t> pending{i};out.topology.facetComponent[i]=component;
            for(std::size_t j=0;j<pending.size();++j)for(auto v:out.facets[pending[j]].vertices)for(auto k:adjacent[v])
                if(out.topology.facetComponent[k]==std::numeric_limits<std::size_t>::max()){out.topology.facetComponent[k]=component;pending.push_back(k);}
        }
    }
};
// This validates current-layer neighbour geometry only after it exists.
// The receipt is an explicit exterior boundary datum; no internal face trace,
// minority plane, whole pure star or stock change is inferred.
void qualifyPrescribedPinWalls(BuildData& out)
{
    auto& pin=*out.pin;pin.prescribedWet=false;pin.prescribedDry=false;
    using Certificate=GeometricVofPublication::PinCertificate;
    if(pin.wallTraceLaw!=Certificate::WallTraceLaw::PrescribedWetDry||pin.kind!=Certificate::Kind::ResolvedCorner)return;
    const auto& contact=out.plicCells[pin.contact];
    require(contact.resolved&&contact.inventory[0]>0&&contact.inventory[1]>0&&contact.measured[0]>0&&contact.measured[1]>0,
        "prescribed pin wall data require the actual resolved mixed contact");
    const auto qualify=[&](int native,unsigned expected,const char* wall){
        const auto& cell=out.plicCells[native];const auto capacity=out.cells[native].capacity;
        const unsigned other=1-expected;
        const bool pure=cell.inventory[expected]>0&&cell.inventory[other]==0&&cell.resolved;
        const auto hx=(static_cast<long double>(cell.rectangle.upper.x)-cell.rectangle.lower.x)/2;
        const bool placeholder=!cell.resolved&&cell.inventory[expected]>cell.inventory[other]&&cell.inventory[other]>0&&
            cell.normal.x==1&&cell.normal.y==0&&cell.offset==(expected==0?hx:-hx);
        if((pure||placeholder)&&cell.measured[expected]==capacity&&cell.measured[other]==0)return;
        std::ostringstream detail;detail.precision(std::numeric_limits<long double>::max_digits10);
        detail<<"prescribed pin wall requires current pure or unresolved expected-majority boundary-placeholder neighbour"
            <<" [wall="<<wall<<",native_cell="<<native<<",expected_phase="<<expected
            <<",L="<<cell.inventory[0]<<",G="<<cell.inventory[1]<<",measured_L="<<cell.measured[0]<<",measured_G="<<cell.measured[1]
            <<",capacity="<<capacity<<",resolved="<<cell.resolved<<",normal=("<<cell.normal.x<<','<<cell.normal.y<<"),offset="<<cell.offset
            <<",clock_position="<<(out.time.position==GeometricVofTimeLayer::Position::Initial?"Initial":
                out.time.position==GeometricVofTimeLayer::Position::Midpoint?"Midpoint":"End")
            <<",clock_step="<<out.time.step<<",clock_begin="<<out.time.begin<<",clock_end="<<out.time.end
            <<",pin_contact_cell="<<pin.contact<<"]";require(false,detail.str());
    };
    qualify(pin.wetCell,0,"wet");qualify(pin.dryCell,1,"dry");
    pin.prescribedWet=true;pin.prescribedDry=true;
}
void reconstructPinnedContact(const QuadtreeMesh& mesh,BuildData& out,
    const std::vector<std::array<long double,2>>& inventory,const GeometricVofConfig& config,bool wallGeometryReady=true)
{
    auto& pin=*out.pin;auto& plic=out.plicCells[pin.contact];auto& cell=out.cells[pin.contact];
    pin.prescribedWet=false;pin.prescribedDry=false;
    const bool prescribed=pin.wallTraceLaw==GeometricVofPublication::PinCertificate::WallTraceLaw::PrescribedWetDry&&
        plic.inventory[0]>0&&plic.inventory[1]>0;
    if(!prescribed){
    if(!(inventory[pin.wetCell][0]>0&&inventory[pin.wetCell][1]==0&&
        inventory[pin.dryCell][0]==0&&inventory[pin.dryCell][1]>0)){
        std::ostringstream detail;detail.precision(std::numeric_limits<long double>::max_digits10);
        detail<<"pinned corner requires unambiguous pure wet-liquid and dry-gas wall neighbours"
            <<" [wet_native_cell="<<pin.wetCell<<",wet_L="<<inventory[pin.wetCell][0]<<",wet_G="<<inventory[pin.wetCell][1]
            <<",dry_native_cell="<<pin.dryCell<<",dry_L="<<inventory[pin.dryCell][0]<<",dry_G="<<inventory[pin.dryCell][1]
            <<",clock_position="<<(out.time.position==GeometricVofTimeLayer::Position::Initial?"Initial":
                out.time.position==GeometricVofTimeLayer::Position::Midpoint?"Midpoint":"End")
            <<",clock_step="<<out.time.step<<",clock_begin="<<out.time.begin<<",clock_end="<<out.time.end
            <<",pin_contact_cell="<<pin.contact<<"]";
        require(false,detail.str());
    }
    }
    const bool mixed=plic.inventory[0]>0&&plic.inventory[1]>0;
    const auto resolution=std::min(config.relativeVolumeTolerance,64.L*std::numeric_limits<double>::epsilon());
    require(!mixed||std::min(plic.inventory[0],plic.inventory[1])>resolution*cell.capacity,
        "unresolved positive pinned contact microphase cannot be clipped or made pure");
    const bool solveLiquid=plic.inventory[0]<=plic.inventory[1];
    const auto plane=Cut::reconstructThroughCorner(plic.rectangle,pin.corner,
        plic.inventory[solveLiquid?0:1]/cell.capacity,out.domain.axisymmetric(),
        solveLiquid?pin.negativeAlongHorizontal:!pin.negativeAlongHorizontal);
    plic.normal=solveLiquid?plane.plane.normal:Vector2D(-plane.plane.normal.x,-plane.plane.normal.y);
    plic.offset=solveLiquid?plane.plane.offset:-plane.plane.offset;
    plic.measured={Cut::volume(plic.rectangle,plic.normal,plic.offset,out.domain.axisymmetric()),
        Cut::volume(plic.rectangle,Vector2D(-plic.normal.x,-plic.normal.y),-plic.offset,out.domain.axisymmetric())};
    const auto segment=Cut::interfaceSegment(plic.rectangle,plic.normal,plic.offset);
    require(segment&&(((*segment)[0].x==plane.localCorner.x&&(*segment)[0].y==plane.localCorner.y)||
        ((*segment)[1].x==plane.localCorner.x&&(*segment)[1].y==plane.localCorner.y)),
        "pinned actual PLIC segment must retain its exact local native corner");
    if(mixed){require(plic.measured[0]>0&&plic.measured[1]>0,"pinned mixed stocks require two positive actual phase volumes");
        pin.kind=GeometricVofPublication::PinCertificate::Kind::ResolvedCorner;
    }else{
        // A pure contact cell has no artificial mixed facet. Only the unique
        // actual opposite-phase shared native face beginning at this vertex
        // can qualify its coincident interface.
        const int opposite=plic.inventory[1]==0?pin.dryCell:pin.wetCell;
        unsigned candidates=0;bool matched=false;
        const auto cx=(static_cast<long double>(plic.rectangle.lower.x)+plic.rectangle.upper.x)/2;
        const auto cy=(static_cast<long double>(plic.rectangle.lower.y)+plic.rectangle.upper.y)/2;
        for(int row=0;row<out.domain.nativeFaceCount();++row){const auto f=out.domain.face(mesh,row);
            if(!((f.negativeCell==pin.contact&&f.positiveCell==opposite)||(f.positiveCell==pin.contact&&f.negativeCell==opposite)))continue;
            const auto& edge=f.horizontal?mesh.hEdges(f.edgeLabel):mesh.vEdges(f.edgeLabel);
            if(mesh.verts(edge.edgeToVert[0]).index()!=pin.vertex&&mesh.verts(edge.edgeToVert[1]).index()!=pin.vertex)continue;
            ++candidates;const auto a=mesh.verts(edge.edgeToVert[0]).coord(),b=mesh.verts(edge.edgeToVert[1]).coord();
            const Point x{static_cast<long double>(a.x)-cx,static_cast<long double>(a.y)-cy},
                y{static_cast<long double>(b.x)-cx,static_cast<long double>(b.y)-cy};
            const auto equal=[](Point a,Point b){return a.x==b.x&&a.y==b.y;};
            matched=(equal((*segment)[0],x)&&equal((*segment)[1],y))||(equal((*segment)[0],y)&&equal((*segment)[1],x));}
        require(candidates==1&&matched,"pure pinned corner requires one actual complete native coincident opposite-phase face");
        pin.kind=GeometricVofPublication::PinCertificate::Kind::NativeCoincident;
    }
    plic.resolved=true;
    for(unsigned phase=0;phase<2;++phase)require(std::isfinite(plic.measured[phase])&&plic.measured[phase]>=0&&
        std::abs(plic.measured[phase]-plic.inventory[phase])/cell.capacity<=config.relativeVolumeTolerance,
        "pinned actual measured/stock gate failed");
    require(std::abs(plic.measured[0]+plic.measured[1]-cell.capacity)/cell.capacity<=config.relativeVolumeTolerance,
        "pinned actual phase closure failed");
    cell.measured=plic.measured[0];
    if(wallGeometryReady)qualifyPrescribedPinWalls(out);
}
// Preserve the original TU arithmetic policy on the literal bodies. Only the
// state/eligibility wrappers below use a local strict scope; no whole-TU or
// CMake FP policy changes, and no inlining across that boundary on MSVC.
#if defined(_MSC_VER)
__declspec(noinline)
#endif
void literalLocalPlicVertices(const QuadtreeMesh& mesh,BuildData& out)
{
    for(int v=0;v<mesh.vertIndexNum();++v){const auto p=mesh.verts(mesh.vertIndexToLabel[v]).coord();
        out.vertices.push_back({static_cast<long double>(p.x)-out.origin.x,static_cast<long double>(p.y)-out.origin.y});}
}
#if defined(_MSC_VER)
__declspec(noinline)
#endif
void literalLocalPlicFaces(const QuadtreeMesh& mesh,BuildData& out,bool reuseSupport,
    std::vector<std::vector<std::size_t>>& cellNodes,std::vector<std::vector<int>>& nodeCells)
{
    const auto n=out.cells.size();
    for(int row=0;row<out.domain.nativeFaceCount();++row){
        const auto f=out.domain.face(mesh,row);const auto& edge=f.horizontal?mesh.hEdges(f.edgeLabel):mesh.vEdges(f.edgeLabel);
        GeometricVofFace face;face.nativeRow=row;face.negativeCell=f.negativeCell;face.positiveCell=f.positiveCell;
        face.vertices={static_cast<std::size_t>(mesh.verts(edge.edgeToVert[0]).index()),
            static_cast<std::size_t>(mesh.verts(edge.edgeToVert[1]).index())};
        const auto x=out.vertices[face.vertices[0]],y=out.vertices[face.vertices[1]];
        if(f.horizontal?x.x>y.x:x.y>y.y)std::swap(face.vertices[0],face.vertices[1]);
        out.faces.push_back(face);
        if(!reuseSupport)for(int ci:{f.negativeCell,f.positiveCell})if(ci>=0)for(auto v:face.vertices){
            require(static_cast<std::size_t>(ci)<n&&v<nodeCells.size(),"local PLIC native incidence outside source");
            cellNodes[ci].push_back(v);nodeCells[v].push_back(ci);}
    }
}
#if defined(_MSC_VER)
__declspec(noinline)
#endif
int localPlicHeightNeighbor(const QuadtreeMesh& mesh,const BuildData& out,int native,unsigned side,
    const LocalPlicGeometry* geometry,const std::vector<std::vector<std::size_t>>& cellNodes,
    const std::vector<std::vector<int>>& nodeCells)
{
    if(native<0)return -1;
    if(geometry&&geometry->heightSides.size()==out.plicCells.size())return geometry->heightSides[native][side];
    const auto& box=out.plicCells[native].rectangle;
    const auto level=mesh.cells(mesh.cellIndexToLabel[native]).level();int result=-1;
    for(auto vertex:cellNodes[native])for(int other:nodeCells[vertex]){
        if(other==native||other<0||out.domain.cellRowByIndex()[other]<0)continue;
        const auto& b=out.plicCells[other].rectangle;
        if(mesh.cells(mesh.cellIndexToLabel[other]).level()!=level||
            static_cast<long double>(b.upper.x)-b.lower.x!=static_cast<long double>(box.upper.x)-box.lower.x||
            static_cast<long double>(b.upper.y)-b.lower.y!=static_cast<long double>(box.upper.y)-box.lower.y)continue;
        const bool match=side<2?b.lower.y==box.lower.y&&b.upper.y==box.upper.y&&
            (side==0?b.upper.x==box.lower.x:b.lower.x==box.upper.x):
            b.lower.x==box.lower.x&&b.upper.x==box.upper.x&&
            (side==2?b.upper.y==box.lower.y:b.lower.y==box.upper.y);
        if(match){if(result>=0&&result!=other)return -1;result=other;}
    }
    return result;
}
using LocalPlicHeightStencil=std::array<int,21>;
LocalPlicHeightStencil localPlicHeightStencil(const QuadtreeMesh& mesh,const BuildData& out,int native,unsigned axis,
    const LocalPlicGeometry* geometry,const std::vector<std::vector<std::size_t>>& cellNodes,
    const std::vector<std::vector<int>>& nodeCells)
{
    LocalPlicHeightStencil result;result.fill(-1);
    const unsigned tangent=1-axis;
    for(int column=0;column<3;++column){const int center=column==1?native:
            localPlicHeightNeighbor(mesh,out,native,2*tangent+(column==2),geometry,cellNodes,nodeCells);
        result[7*column+3]=center;
        for(int direction:{-1,1}){int cell=center;
            for(int row=1;row<=3;++row){cell=localPlicHeightNeighbor(mesh,out,cell,2*axis+(direction>0),geometry,cellNodes,nodeCells);
                result[7*column+3+direction*row]=cell;}}
    }
    return result;
}
std::optional<Vector2D> localPlicInventoryHeightNormal(const QuadtreeMesh& mesh,const BuildData& out,int native,
    Vector2D seed,const LocalPlicGeometry* geometry,const std::vector<std::vector<std::size_t>>& cellNodes,
    const std::vector<std::vector<int>>& nodeCells)
{
    if(out.domain.axisymmetric())return std::nullopt;
    const auto& source=out.plicCells[native].rectangle;
    const unsigned first=std::abs(seed.x)>=std::abs(seed.y)?0u:1u;
    for(unsigned attempt=0;attempt<2;++attempt){const unsigned axis=attempt?1-first:first;
        const auto stencil=localPlicHeightStencil(mesh,out,native,axis,geometry,cellNodes,nodeCells);
        const long double spacing=axis?static_cast<long double>(source.upper.x)-source.lower.x:
            static_cast<long double>(source.upper.y)-source.lower.y;
        const long double origin=axis?(static_cast<long double>(source.lower.y)+source.upper.y)/2:
            (static_cast<long double>(source.lower.x)+source.upper.x)/2;
        if(!(std::isfinite(spacing)&&spacing>0))continue;
        // Qualify the chart with native finite phase-fraction pairs. This is
        // only a stencil-selection law, not a claim that a small inventory
        // cannot affect the complete height. Every raw stock, including tails
        // absorbed by the endpoint sum, still enters both height integrals.
        for(int width=3;width<=7;++width)for(int low=-3;low<=0;++low){const int high=low+width-1;
            if(high<0||high>3)continue;
            std::array<long double,3> heights{};int lowerPhase=-1;bool qualified=true;
            for(int column=0;column<3&&qualified;++column){
                const int lower=stencil[7*column+low+3],upper=stencil[7*column+high+3];
                if(lower<0||upper<0){qualified=false;break;}
                const auto phase=[](const GeometricVofPlicCell& c){const auto total=c.inventory[0]+c.inventory[1];
                    return total==c.inventory[0]?0:total==c.inventory[1]?1:-1;};
                const int lp=phase(out.plicCells[lower]),up=phase(out.plicCells[upper]);
                if(lp<0||up<0||lp==up||(lowerPhase>=0&&lp!=lowerPhase)){qualified=false;break;}
                lowerPhase=lp;std::array<long double,2> stock{},compensation{};
                std::array<long double,2> previous=lp==0?std::array<long double,2>{1.L,0.L}:std::array<long double,2>{0.L,1.L};
                for(int row=low;row<=high;++row){const int cell=stencil[7*column+row+3];
                    if(cell<0){qualified=false;break;}const auto& q=out.plicCells[cell].inventory;
                    const auto total=q[0]+q[1];if(!(std::isfinite(total)&&total>0)){qualified=false;break;}
                    const std::array<long double,2> fraction{q[0]/total,q[1]/total};
                    // A reversal must be distinguishable in BOTH phases.
                    // A tie in either fraction has no resolved ordering at
                    // this storage boundary; it never changes a raw stock.
                    const bool reversed=lp==0?fraction[0]>previous[0]&&fraction[1]<previous[1]:
                        fraction[0]<previous[0]&&fraction[1]>previous[1];
                    if(reversed){qualified=false;break;}previous=fraction;
                    for(unsigned p=0;p<2;++p){const auto term=q[p]-compensation[p],next=stock[p]+term;
                        compensation[p]=(next-stock[p])-term;stock[p]=next;}
                }
                if(!qualified)break;
                const auto& a=out.plicCells[lower].rectangle;const auto& b=out.plicCells[upper].rectangle;
                const long double lo=(axis?static_cast<long double>(a.lower.y):a.lower.x)-origin;
                const long double hi=(axis?static_cast<long double>(b.upper.y):b.upper.x)-origin;
                const auto liquid=lp==0?lo+stock[0]/spacing:hi-stock[0]/spacing;
                const auto gas=lp==0?hi-stock[1]/spacing:lo+stock[1]/spacing;
                const auto roundoff=64*std::numeric_limits<long double>::epsilon()*(hi-lo);
                if(!(std::isfinite(liquid)&&std::isfinite(gas)&&std::abs(liquid-gas)<=roundoff)){
                    qualified=false;break;}
                heights[column]=(liquid+gas)/2;
            }
            if(!qualified)continue;
            const auto centerSlope=(heights[2]-heights[0])/(2*spacing);
            const auto second=(heights[2]-2*heights[1]+heights[0])/(spacing*spacing);
            if(!(std::isfinite(centerSlope)&&std::isfinite(second)))continue;
            const long double sign=lowerPhase==0?1.L:-1.L;
            const auto normalAt=[&](long double slope){const auto length=std::hypot(1.L,slope);
                return axis?Vector2D(static_cast<double>(-sign*slope/length),static_cast<double>(sign/length)):
                    Vector2D(static_cast<double>(sign/length),static_cast<double>(-sign*slope/length));};
            auto normal=normalAt(centerSlope);
            // Column heights are tangential cell averages. A corner cut's
            // interface is away from that column centre. Evaluate the same
            // quadratic height derivative at its actual PLIC segment midpoint,
            // with two fixed updates; each offset still solves the OWN stock.
            const auto& inventory=out.plicCells[native].inventory;
            const bool solveLiquid=inventory[0]<=inventory[1];
            for(unsigned update=0;update<2&&second!=0;++update){
                const auto cut=Cut::reconstruct(source,solveLiquid?normal:Vector2D(-normal.x,-normal.y),
                    inventory[solveLiquid?0:1]/out.cells[native].capacity,false);
                const auto segment=Cut::interfaceSegment(source,cut.normal,cut.offset);if(!segment)break;
                const auto tangent=axis?((*segment)[0].x+(*segment)[1].x)/2:((*segment)[0].y+(*segment)[1].y)/2;
                const auto slope=centerSlope+second*tangent;if(!std::isfinite(slope))break;
                const auto candidate=normalAt(slope);
                if(!(std::isfinite(candidate.x)&&std::isfinite(candidate.y)&&
                    static_cast<long double>(candidate.x)*seed.x+static_cast<long double>(candidate.y)*seed.y>0))break;
                normal=candidate;
            }
            if(static_cast<long double>(normal.x)*seed.x+static_cast<long double>(normal.y)*seed.y>0)return normal;
        }
    }
    return std::nullopt;
}
// Keep the original LS and minor-stock offset arithmetic in its literal scope.
#if defined(_MSC_VER)
__declspec(noinline)
#endif
void literalLocalPlicMixed(const QuadtreeMesh& mesh,BuildData& out,std::size_t ci,const std::vector<Point>& centers,
    const std::vector<long double>& widths,const std::vector<std::vector<std::size_t>>& actualCellNodes,
    const std::vector<std::vector<int>>& actualNodeCells,const LocalPlicGeometry* geometry,GeometricVofReport& report)
{
    auto& plic=out.plicCells[ci];const auto& cell=out.cells[ci];const bool rz=out.domain.axisymmetric();
    std::vector<int> coldNeighbors;
    if(!geometry){coldNeighbors.reserve(32);
        for(auto v:actualCellNodes[ci])for(int k:actualNodeCells[v])if(k!=static_cast<int>(ci))coldNeighbors.push_back(k);
        std::sort(coldNeighbors.begin(),coldNeighbors.end());coldNeighbors.erase(std::unique(coldNeighbors.begin(),coldNeighbors.end()),coldNeighbors.end());}
    const auto& neighbors=geometry?geometry->neighbors[ci]:coldNeighbors;
    if(geometry)++report.localPlicReuse.nodeRings;
    require(neighbors.size()<=32,"local PLIC node-ring support is not bounded");
    long double xx=0,xy=0,yy=0,bx=0,by=0;
    for(int k:neighbors){const auto d=minus(centers[k],centers[ci]);const auto x=d.x/widths[ci],y=d.y/widths[ci];
        const auto d2=x*x+y*y;require(std::isfinite(d2)&&d2>0,"invalid local PLIC node-ring geometry");
        const auto w=1/d2,da=static_cast<long double>(out.alpha[k])-out.alpha[ci];
        xx+=w*x*x;xy+=w*x*y;yy+=w*y*y;bx+=w*x*da;by+=w*y*da;}
    if(rz&&mesh.minCoord.x==0&&plic.rectangle.lower.x==0){
        // Scalar alpha is even at the TRUE rotation axis. Complete
        // only this node ring with read-only reflected sample sites,
        // including the source's own mirror; no negative-r physical
        // measure, ghost inventory or prescribed normal is created.
        std::vector<int> coldReflected;
        if(!geometry){for(auto v:actualCellNodes[ci])if(mesh.verts(mesh.vertIndexToLabel[v]).coord().x==0)
            for(int k:actualNodeCells[v])coldReflected.push_back(k);
            std::sort(coldReflected.begin(),coldReflected.end());coldReflected.erase(std::unique(coldReflected.begin(),coldReflected.end()),coldReflected.end());}
        const auto& reflected=geometry?geometry->reflected[ci]:coldReflected;
        require(reflected.size()<=neighbors.size()+1,"RZ axis reflection leaves the original node ring");
        for(int k:reflected){const Point image{-centers[k].x,centers[k].y};const auto d=minus(image,centers[ci]);
            const auto x=d.x/widths[ci],y=d.y/widths[ci],d2=x*x+y*y;
            require(std::isfinite(d2)&&d2>0,"invalid true-axis local PLIC reflected sample");
            const auto w=1/d2,da=static_cast<long double>(out.alpha[k])-out.alpha[ci];
            xx+=w*x*x;xy+=w*x*y;yy+=w*y*y;bx+=w*x*da;by+=w*y*da;}
    }
    const auto determinant=xx*yy-xy*xy;
    require(std::isfinite(determinant)&&std::isfinite(xx)&&std::isfinite(yy)&&std::isfinite(bx)&&std::isfinite(by),
        "nonfinite local PLIC normal reconstruction");
    bool normalResolved=false;
    if(determinant>32*std::numeric_limits<long double>::epsilon()*(xx+yy)*(xx+yy)){
        const auto gx=(yy*bx-xy*by)/determinant,gy=(xx*by-xy*bx)/determinant;const auto length=std::hypot(gx,gy);
        require(std::isfinite(length),"nonfinite local PLIC gradient");
        if(length>0){plic.normal={static_cast<double>(-gx/length),static_cast<double>(-gy/length)};normalResolved=true;}
    }
    if(normalResolved)if(const auto height=localPlicInventoryHeightNormal(mesh,out,static_cast<int>(ci),plic.normal,
        geometry,actualCellNodes,actualNodeCells))plic.normal=*height;
    // Solve the smaller retained stock directly. A tiny gas inventory
    // is never reconstructed by rounding liquid alpha to one first.
    const bool solveLiquid=plic.inventory[0]<=plic.inventory[1];
    const auto normal=solveLiquid?plic.normal:Vector2D(-plic.normal.x,-plic.normal.y);
    const auto fraction=plic.inventory[solveLiquid?0:1]/cell.capacity;
    const auto reconstruction=Cut::reconstruct(plic.rectangle,normal,fraction,rz);
    plic.normal=solveLiquid?reconstruction.normal:Vector2D(-reconstruction.normal.x,-reconstruction.normal.y);
    plic.offset=solveLiquid?reconstruction.offset:-reconstruction.offset;
    plic.measured={Cut::volume(plic.rectangle,plic.normal,plic.offset,rz),
        Cut::volume(plic.rectangle,Vector2D(-plic.normal.x,-plic.normal.y),-plic.offset,rz)};
    plic.resolved=normalResolved&&plic.measured[0]>0&&plic.measured[1]>0;
}
#if defined(_MSC_VER) && !defined(__FAST_MATH__) && (!defined(_M_FP_FAST) || !_M_FP_FAST)
#pragma float_control(precise,on,push)
#pragma float_control(except,on)
#pragma fenv_access(on)
#endif
void localPlicVertices(const QuadtreeMesh& mesh,BuildData& out,const LocalPlicGeometry* geometry,LocalPlicQuietState& record)
{
    if(geometry&&geometry->vertexState.matches()){
        for(const auto& p:geometry->vertices)out.vertices.push_back(p);
        ++out.report.localPlicReuse.vertexSets;record=geometry->vertexState;return;}
    const auto before=LocalPlicQuietState::read();
    literalLocalPlicVertices(mesh,out);
    record.completed(before);
}
void localPlicFaces(const QuadtreeMesh& mesh,BuildData& out,bool reuseSupport,
    std::vector<std::vector<std::size_t>>& cellNodes,std::vector<std::vector<int>>& nodeCells,
    const LocalPlicGeometry* geometry,LocalPlicQuietState& record)
{
    if(geometry&&reuseSupport&&geometry->faceState.matches()){
        for(const auto& face:geometry->faces)out.faces.push_back(face);
        ++out.report.localPlicReuse.faceSets;record=geometry->faceState;return;}
    const auto before=LocalPlicQuietState::read();
    literalLocalPlicFaces(mesh,out,reuseSupport,cellNodes,nodeCells);
    record.completed(before);
}
bool unchangedLocalPlicMixed(const QuadtreeMesh& mesh,std::size_t ci,const BuildData& out,const BuildData* source,
    const LocalPlicGeometry* geometry,const GeometricVofConfig& config,
    const std::vector<std::vector<std::size_t>>& cellNodes,const std::vector<std::vector<int>>& nodeCells) noexcept
{
    if(!source||!geometry||!source->plicCells[ci].resolved
        ||!sameLocalPlicBits(config.relativeVolumeTolerance,source->localConfig.relativeVolumeTolerance)
        ||config.maximumVolumeIterations!=source->localConfig.maximumVolumeIterations
        ||config.maximumBracketExpansions!=source->localConfig.maximumBracketExpansions
        ||config.localResolution!=source->localConfig.localResolution)return false;
    const auto& old=source->plicCells[ci];const auto& now=out.plicCells[ci];
    if(!sameLocalPlicBits(old.inventory[0],now.inventory[0])||!sameLocalPlicBits(old.inventory[1],now.inventory[1])
        ||!sameLocalPlicBits(source->alpha[ci],out.alpha[ci]))return false;
    for(int k:geometry->neighbors[ci])if(!sameLocalPlicBits(source->alpha[k],out.alpha[k]))return false;
    for(int k:geometry->reflected[ci])if(!sameLocalPlicBits(source->alpha[k],out.alpha[k]))return false;
    // Qualification can change outside the LS node ring, even when rounded
    // alpha is identical. Both possible charts and both independent stocks
    // therefore belong to the cached normal's complete dependency set.
    if(!out.domain.axisymmetric())for(unsigned axis=0;axis<2;++axis){
        const auto stencil=localPlicHeightStencil(mesh,out,static_cast<int>(ci),axis,geometry,cellNodes,nodeCells);
        for(int k:stencil)if(k>=0)for(unsigned phase=0;phase<2;++phase)
            if(!sameLocalPlicBits(source->plicCells[k].inventory[phase],out.plicCells[k].inventory[phase]))return false;
    }
    return true;
}
void localPlicMixed(const QuadtreeMesh& mesh,BuildData& out,std::size_t ci,const std::vector<Point>& centers,
    const std::vector<long double>& widths,const std::vector<std::vector<std::size_t>>& actualCellNodes,
    const std::vector<std::vector<int>>& actualNodeCells,const LocalPlicGeometry* geometry,
    const BuildData* source,const GeometricVofConfig& config,GeometricVofReport& report)
{
    auto& plic=out.plicCells[ci];auto& receipt=report.localPlicReuse;
    const bool unchanged=unchangedLocalPlicMixed(mesh,ci,out,source,geometry,config,actualCellNodes,actualNodeCells);
    if(unchanged&&source->localMixedStates&&source->localMixedStates->size()==out.plicCells.size()
        &&(*source->localMixedStates)[ci].matches()){
        const auto& old=source->plicCells[ci];plic.normal=old.normal;plic.offset=old.offset;plic.measured=old.measured;plic.resolved=old.resolved;
        if(out.localMixedStates&&out.localMixedStates->size()==out.plicCells.size())(*out.localMixedStates)[ci]=(*source->localMixedStates)[ci];
        ++receipt.mixedReused;return;}
    if(source&&geometry&&source->plicCells[ci].resolved&&!unchanged)++receipt.mixedDependencyChanged;
    ++receipt.mixedReconstructed;const auto before=LocalPlicQuietState::read();
    literalLocalPlicMixed(mesh,out,ci,centers,widths,actualCellNodes,actualNodeCells,geometry,report);
    if(plic.resolved&&out.localMixedStates&&out.localMixedStates->size()==out.plicCells.size())(*out.localMixedStates)[ci].completed(before);
}
#if defined(_MSC_VER) && !defined(__FAST_MATH__) && (!defined(_M_FP_FAST) || !_M_FP_FAST)
#if !defined(_M_FP_STRICT)
#pragma fenv_access(off)
#endif
#pragma float_control(pop)
#endif
void retainLocalPlicGeometry(const QuadtreeMesh& mesh,BuildData& out,std::vector<Point>& centers,
    const LocalPlicQuietState& vertexState,const LocalPlicQuietState& faceState,const BuildData* predecessor) noexcept
{
    LocalPlicSupportEnvironment environment;if(!environment.captured)return;
    const auto support=out.localNativeSupport;
    const auto recordBytes=localPlicRecordBytes(out.localMixedStates);
    out.report.localPlicReuse.retainedBytes=recordBytes+(support?support->bytes:0);
    if(!support)return;
    try {
        constexpr std::size_t limit=64u*1024u*1024u;
        const auto base=out.report.localPlicReuse.retainedBytes;
        const auto borrowed=predecessor?localPlicRecordBytes(predecessor->localMixedStates):0;
        std::size_t bytes=base;
        const auto account=[&](std::size_t n,std::size_t width){if(bytes>limit||(width&&n>(limit-bytes)/width))return false;bytes+=n*width;return true;};
        if(!account(borrowed,1)){out.localGeometry.reset();return;}
        if(out.localGeometry){if(account(out.localGeometry->bytes,1))out.report.localPlicReuse.retainedBytes=bytes-borrowed;
            else out.localGeometry.reset();return;}
        const auto n=out.cells.size();if(centers.size()!=n)return;
        // Preflight before allocating the new sibling. Ring storage reserves
        // the actual integer incidence count (including duplicates) exactly.
        if(!account(sizeof(LocalPlicGeometry)+512u,1)||!account(out.domain.cellRowByIndex().size(),sizeof(int))
            ||!account(out.domain.cellIndexByRow().size(),sizeof(int))||!account(n,sizeof(Cut::Rectangle))
            ||!account(out.vertices.size(),sizeof(Point))||!account(out.faces.size(),sizeof(GeometricVofFace))
            ||!account(centers.capacity(),sizeof(Point))||!account(n,2*sizeof(std::vector<int>))
            ||(!out.domain.axisymmetric()&&!account(n,sizeof(std::array<int,4>))))return;
        for(std::size_t ci=0;ci<n;++ci)if(out.domain.cellRowByIndex()[ci]>=0){std::size_t count=0;
            for(auto v:support->cellNodes[ci]){const auto more=support->nodeCells[v].size();if(more>limit||count>limit-more)return;count+=more;}
            if(!account(count,2*sizeof(int))||!account(64u,1))return;}
        auto geometry=std::shared_ptr<LocalPlicGeometry>(new LocalPlicGeometry,[](LocalPlicGeometry* p) noexcept {LocalPlicSupportEnvironment released;delete p;});
        geometry->identity=out.domain.meshIdentity();geometry->radial=out.domain.axisymmetric();geometry->origin=out.origin;
        geometry->selectedRows=out.domain.cellRowByIndex();geometry->selectedCells=out.domain.cellIndexByRow();
        geometry->rectangles.reserve(n);for(const auto& p:out.plicCells)geometry->rectangles.push_back(p.rectangle);
        geometry->vertices=out.vertices;geometry->faces=out.faces;geometry->neighbors.resize(n);geometry->reflected.resize(n);
        if(!geometry->radial){geometry->heightSides.resize(n);for(auto& sides:geometry->heightSides)sides.fill(-1);}
        for(std::size_t ci=0;ci<n;++ci){if(geometry->selectedRows[ci]<0)continue;
            std::size_t count=0;for(auto v:support->cellNodes[ci])count+=support->nodeCells[v].size();
            auto& ring=geometry->neighbors[ci];ring.reserve(count);
            for(auto v:support->cellNodes[ci])for(int k:support->nodeCells[v])if(k!=static_cast<int>(ci))ring.push_back(k);
            std::sort(ring.begin(),ring.end());ring.erase(std::unique(ring.begin(),ring.end()),ring.end());
            if(!geometry->radial)for(unsigned side=0;side<4;++side)geometry->heightSides[ci][side]=
                localPlicHeightNeighbor(mesh,out,static_cast<int>(ci),side,nullptr,support->cellNodes,support->nodeCells);
            if(geometry->radial&&mesh.minCoord.x==0&&out.plicCells[ci].rectangle.lower.x==0){auto& reflected=geometry->reflected[ci];reflected.reserve(count);
                for(auto v:support->cellNodes[ci])if(mesh.verts(mesh.vertIndexToLabel[v]).coord().x==0)
                    for(int k:support->nodeCells[v])reflected.push_back(k);
                std::sort(reflected.begin(),reflected.end());reflected.erase(std::unique(reflected.begin(),reflected.end()),reflected.end());}}
        // Recheck actual capacities, not only requested sizes, before publication.
        bytes=base+borrowed;
        if(!account(sizeof(LocalPlicGeometry)+512u,1)||!account(geometry->selectedRows.capacity(),sizeof(int))
            ||!account(geometry->selectedCells.capacity(),sizeof(int))||!account(geometry->rectangles.capacity(),sizeof(Cut::Rectangle))
            ||!account(geometry->vertices.capacity(),sizeof(Point))||!account(geometry->faces.capacity(),sizeof(GeometricVofFace))
            ||!account(centers.capacity(),sizeof(Point))
            ||!account(geometry->neighbors.capacity()+geometry->reflected.capacity(),sizeof(std::vector<int>))
            ||!account(geometry->heightSides.capacity(),sizeof(std::array<int,4>)))return;
        for(const auto* lists:{&geometry->neighbors,&geometry->reflected})for(const auto& list:*lists)
            if(!account(list.capacity(),sizeof(int))||(!list.empty()&&!account(32u,1)))return;
        geometry->centers.swap(centers);
        geometry->vertexState=vertexState;geometry->faceState=faceState;geometry->bytes=bytes-base-borrowed;
        out.localGeometry=std::move(geometry);out.report.localPlicReuse.retainedBytes=bytes-borrowed;
    }catch(...){/* Private ephemeral preparation cannot change a successful publication. */}
}

#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
#pragma float_control(precise,on,push)
#pragma float_control(except,on)
#pragma fenv_access(on)
struct LocalPlicExecutionState
{
    std::fenv_t native{};unsigned csr=0;int error=0;bool captured=false;
    __declspec(noinline) bool capture() noexcept {error=errno;csr=_mm_getcsr();captured=std::fegetenv(&native)==0;errno=error;return captured;}
    __declspec(noinline) void restore() const noexcept {if(captured&&std::fesetenv(&native)!=0)std::terminate();_mm_setcsr(csr);errno=error;}
    bool sameControl(const LocalPlicExecutionState& other) const noexcept {
        auto a=native,b=other.native;a._Fe_stat=0;b._Fe_stat=0;
        return captured&&other.captured&&std::memcmp(&a,&b,sizeof a)==0&&(csr&~0x3fu)==(other.csr&~0x3fu);
    }
};
struct LocalPlicExecutionCleanup
{
    LocalPlicExecutionState state;bool armed=false;
    void arm() noexcept {armed=state.capture();if(!armed)std::terminate();}
    ~LocalPlicExecutionCleanup() noexcept {if(armed)state.restore();}
};
#if !defined(_M_FP_STRICT)
#pragma fenv_access(off)
#endif
#pragma float_control(pop)
#endif
void buildLocalPlic(const QuadtreeMesh& mesh,BuildData& out,
    const std::vector<std::array<long double,2>>& inventory,const GeometricVofConfig& config,
    const std::vector<FiniteVolume::Dyadic>* predictedDefects=nullptr,const BuildData* predecessor=nullptr)
{
    const auto n=static_cast<std::size_t>(mesh.cellIndexNum());const bool rz=out.domain.axisymmetric();
    require(config.localResolution==GeometricVofLocalResolution::VolumeOnly||
        config.localResolution==GeometricVofLocalResolution::NativeInterfaceSupport,"unknown local PLIC resolution law");
    require(config.localResolution==GeometricVofLocalResolution::VolumeOnly||(!rz&&!out.pin),
        "native interface support requires unpinned XY LocalPlic");
    out.localConfig=config;
    require(!predictedDefects||predictedDefects->size()==n,"derived local PLIC closure dimensions differ");
    out.kind=GeometricVofGeometryKind::LocalPlic;out.cells.resize(n);out.plicCells.resize(n);out.alpha.resize(n);
    out.capacityDefects.resize(n);
    const auto& active=out.domain.cellRowByIndex();
    require(active.size()==n,"local PLIC selected map must cover native indices");
    const auto support=out.localNativeSupport;
    const auto* geometry=qualifiedLocalPlicGeometry(mesh,out);
    if(!geometry)out.localGeometry.reset();
    std::vector<Point> coldCenters(geometry?0:n);std::vector<long double> widths(n);
    const auto& centers=geometry?geometry->centers:coldCenters;
    prepareLocalPlicRecords(out,n,predecessor);
    LocalPlicQuietState vertexState,faceState;
    const bool reuseSupport=support&&support->identity==out.domain.meshIdentity()
        &&support->radial==rz&&support->cells==n
        &&support->vertices==static_cast<std::size_t>(mesh.vertIndexNum())
        &&support->faces==static_cast<std::size_t>(out.domain.nativeFaceCount());
    std::vector<std::vector<std::size_t>> cellNodes(reuseSupport?0:n);
    std::vector<std::vector<int>> nodeCells(reuseSupport?0:static_cast<std::size_t>(mesh.vertIndexNum()));
    const auto& actualCellNodes=reuseSupport?support->cellNodes:cellNodes;
    const auto& actualNodeCells=reuseSupport?support->nodeCells:nodeCells;
    out.vertices.reserve(actualNodeCells.size());
    localPlicVertices(mesh,out,geometry,vertexState);
    const auto initializeCell=[&](std::size_t ci,GeometricVofReport& report){
        auto& cell=out.cells[ci];auto& plic=out.plicCells[ci];
        if(active[ci]<0){
            if(!(inventory[ci][0]==0&&inventory[ci][1]==0))require(false,
                "inactive local PLIC cell must have zero independent stocks cell="+std::to_string(ci));
            require(!predictedDefects||(*predictedDefects)[ci].compare(FiniteVolume::Dyadic{})==0,
                "inactive local PLIC cell cannot carry a closure certificate");
            plic.normal={0,0};return;
        }
        const auto& native=mesh.cells(mesh.cellIndexToLabel[ci]);
        if(geometry){plic.rectangle=geometry->rectangles[ci];++report.localPlicReuse.geometryCells;}
        else plic.rectangle={mesh.verts(native.cellToVert[0]).coord(),mesh.verts(native.cellToVert[2]).coord()};
        cell.nativeCell=static_cast<int>(ci);cell.capacity=Cut::capacity(plic.rectangle,rz);
        plic.inventory=inventory[ci];
        if(!predictedDefects){
            // The public raw-stock entry remains strict. No q>C can create
            // its own certificate by presenting its observed residual.
            for(auto value:plic.inventory)if(!(std::isfinite(value)&&value>=0&&value<=cell.capacity))require(false,
                "local PLIC inventory outside actual capacity; no clipping cell="+std::to_string(ci));
        }else for(auto value:plic.inventory)if(!(std::isfinite(value)&&value>=0))require(false,
            "derived local PLIC phase stock is nonfinite/negative cell="+std::to_string(ci));
        auto& exactDefect=out.capacityDefects[ci];exactDefect=FiniteVolume::Dyadic(plic.inventory[0]);
        exactDefect.add(plic.inventory[1]);exactDefect.add(cell.capacity,-1);
        if(predictedDefects)if(!(exactDefect.compare((*predictedDefects)[ci])==0))require(false,
            "derived local PLIC predecessor closure certificate mismatch cell="+std::to_string(ci));
        const auto ceiling=FiniteVolume::Dyadic::product({config.relativeVolumeTolerance,cell.capacity});
        if(!(exactDefect.absolute().compare(ceiling)<=0))require(false,
            "local PLIC exact closure exceeds original geometry budget cell="+std::to_string(ci));
        if(predictedDefects){auto phaseBound=FiniteVolume::Dyadic(cell.capacity);phaseBound.add(exactDefect.absolute());
            for(auto value:plic.inventory)if(!(FiniteVolume::Dyadic(value).compare(phaseBound)<=0))require(false,
                "derived local PLIC phase exceeds certified capacity cell="+std::to_string(ci));}
        const auto stockDefect=std::abs((plic.inventory[0]+plic.inventory[1])-cell.capacity)/cell.capacity;
        if(!(stockDefect<=config.relativeVolumeTolerance))require(false,"local PLIC phase stock closure failed cell="+std::to_string(ci));
        const auto row=active[ci];
        const auto capacityDefect=std::abs(cell.capacity-out.domain.cellVolumes()[row])/cell.capacity;
        require(capacityDefect<=config.relativeVolumeTolerance,"local PLIC Core/native capacity mismatch");
        report.maximumRelativeCapacityResidual=std::max(report.maximumRelativeCapacityResidual,capacityDefect);
        cell.inventory=plic.inventory[0];out.alpha[ci]=static_cast<double>(plic.inventory[0]/cell.capacity);
        if(predictedDefects&&(plic.inventory[0]>cell.capacity||plic.inventory[1]>cell.capacity)){
            // Geometry diagnostic only. A certified pure-liquid stock may be
            // one native ULP above C; never feed q/C>1 to the normal stencil.
            // Both authoritative stocks remain unchanged. Ordinary bounded
            // cells retain the original diagnostic arithmetic bit for bit.
            out.alpha[ci]=plic.inventory[0]==0?0:plic.inventory[1]==0?1:
                static_cast<double>(plic.inventory[0]/(plic.inventory[0]+plic.inventory[1]));
        }
        require(std::isfinite(out.alpha[ci])&&out.alpha[ci]>=0&&out.alpha[ci]<=1,"local PLIC diagnostic alpha outside bounds");
        if(!geometry){const auto center=native.centerCoord();coldCenters[ci]={center.x,center.y};}
        // Two original subtractions/min are cheaper than a per-cell complete
        // FP receipt. Keep this arithmetic at its original place and order.
        widths[ci]=std::min(static_cast<long double>(plic.rectangle.upper.x)-plic.rectangle.lower.x,
            static_cast<long double>(plic.rectangle.upper.y)-plic.rectangle.lower.y);
    };
    const auto prepareFaces=[&]{
    // Only actual native incidence/endpoints, including real hanging nodes.
    // No shared scalar, collar, facet graph, or global interface is built.
    out.faces.reserve(static_cast<std::size_t>(out.domain.nativeFaceCount()));
    localPlicFaces(mesh,out,reuseSupport,cellNodes,nodeCells,geometry,faceState);
    if(!reuseSupport){
        for(auto& list:nodeCells){std::sort(list.begin(),list.end());list.erase(std::unique(list.begin(),list.end()),list.end());}
        for(std::size_t ci=0;ci<n;++ci){auto& list=cellNodes[ci];
            if(active[ci]<0){require(list.empty(),"inactive local PLIC cell has selected incidence");continue;}
            std::sort(list.begin(),list.end());list.erase(std::unique(list.begin(),list.end()),list.end());
            require(list.size()>=4&&list.size()<=8,"local PLIC requires real 2:1 native leaf boundary");}
    }
    };
    PositiveSum stockSum,measureSum;
    // Geometry resolution only: both independent stocks and diagnostic alpha
    // were retained above. Do not inverse-fit a numerical residual interface
    // below this capacity-scaled resolution; transport uses its existing
    // unresolved donor-fraction law for these SAME unmodified stocks.
    const auto geometryResolution=std::min(config.relativeVolumeTolerance,
        64.L*std::numeric_limits<double>::epsilon());
    const auto reconstructCell=[&](std::size_t ci,GeometricVofReport& report){
        if(active[ci]<0)return;
        auto& plic=out.plicCells[ci];auto& cell=out.cells[ci];
        if(out.pin&&out.pin->contact==static_cast<int>(ci)){
            ++report.localPlicReuse.pinnedReconstructed;
            reconstructPinnedContact(mesh,out,inventory,config,false);
            if(plic.inventory[0]>0&&plic.inventory[1]>0)++report.mixedCells;
            else if(plic.inventory[1]==0)++report.pureLiquidCells;else ++report.pureGasCells;
        }else if(plic.inventory[0]==0||plic.inventory[1]==0){
            // Exact pure geometry is the native rectangle or the empty set.
            const auto hx=(static_cast<long double>(plic.rectangle.upper.x)-plic.rectangle.lower.x)/2;
            const bool liquid=plic.inventory[1]==0;plic.offset=liquid?hx:-hx;
            plic.measured=liquid?std::array<long double,2>{cell.capacity,0}:std::array<long double,2>{0,cell.capacity};
            plic.resolved=true;if(liquid)++report.pureLiquidCells;else ++report.pureGasCells;
        }else if(std::min(plic.inventory[0],plic.inventory[1])<=geometryResolution*cell.capacity){
            // A legitimate strict majority plane, not a changed inventory.
            // Its exact rectangle/empty measures carry no resolvable minor
            // interface, so this stock-mixed cell remains explicitly unresolved.
            ++report.mixedCells;
            const auto hx=(static_cast<long double>(plic.rectangle.upper.x)-plic.rectangle.lower.x)/2;
            const bool liquid=plic.inventory[0]>plic.inventory[1];
            plic.normal={1,0};plic.offset=liquid?hx:-hx;
            plic.measured=liquid?std::array<long double,2>{cell.capacity,0}:std::array<long double,2>{0,cell.capacity};
            plic.resolved=false;
        }else{
            ++report.mixedCells;
            localPlicMixed(mesh,out,ci,centers,widths,actualCellNodes,actualNodeCells,geometry,predecessor,config,report);
        }
        for(std::size_t phase=0;phase<2;++phase){
            require(std::isfinite(plic.measured[phase])&&plic.measured[phase]>=0,"invalid actual local PLIC phase measure");
            const auto defect=std::abs(plic.measured[phase]-plic.inventory[phase])/cell.capacity;
            if(!(defect<=config.relativeVolumeTolerance))require(false,"local PLIC measured/stock capacity-scaled gate failed cell="+std::to_string(ci));
            report.maximumRelativeVolumeResidual=std::max(report.maximumRelativeVolumeResidual,defect);}
        const auto closure=std::abs((plic.measured[0]+plic.measured[1])-cell.capacity)/cell.capacity;
        require(closure<=config.relativeVolumeTolerance,"actual local PLIC phase closure failed");
        report.maximumRelativeCapacityResidual=std::max(report.maximumRelativeCapacityResidual,closure);
        cell.measured=plic.measured[0];
    };
    const auto sumCell=[&](std::size_t ci){if(active[ci]>=0){stockSum.add(out.cells[ci].inventory);measureSum.add(out.cells[ci].measured);}};
    bool completed=false;
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
    // Cleanup dies after optional reports and preserves the actual final state.
    LocalPlicExecutionCleanup cleanup;
    LocalPlicExecutionState entry;const auto quiet=LocalPlicQuietState::read();
    const bool eligible=n>=128&&n<=static_cast<std::size_t>(std::numeric_limits<std::ptrdiff_t>::max())
        &&quiet.valid&&!omp_in_parallel()&&omp_get_max_threads()>1
        &&n/static_cast<std::size_t>(omp_get_max_threads())>=16&&entry.capture();
    if(eligible){
        const auto savedReport=out.report;const auto savedPin=out.pin;
        try{
            const int team=omp_get_max_threads();
            if(static_cast<std::size_t>(team)>64u*1024u*1024u/(sizeof(GeometricVofReport)+16384u))
                throw std::length_error("optional local PLIC team scratch cap");
            std::vector<GeometricVofReport> reports(static_cast<std::size_t>(team));
            if(errno!=entry.error)throw std::runtime_error("optional local PLIC setup errno");
            entry.restore();LocalPlicExecutionState phase=entry;std::atomic<int> failed{0};
            unsigned long nativeFlags=0;unsigned csrFlags=0;
#pragma omp parallel num_threads(team)
            {
                LocalPlicExecutionState saved;const bool captured=saved.capture();
                if(!captured)failed.store(1,std::memory_order_relaxed);if(captured)phase.restore();
                auto& report=reports[static_cast<std::size_t>(omp_get_thread_num())];
                const auto rows=[&](std::size_t first,std::size_t last,const auto& kernel){
#pragma omp for schedule(dynamic,8)
                    for(std::ptrdiff_t ci=static_cast<std::ptrdiff_t>(first);ci<static_cast<std::ptrdiff_t>(last);++ci){
                        if(failed.load(std::memory_order_relaxed))continue;
                        try{kernel(static_cast<std::size_t>(ci),report);}catch(...){failed.store(1,std::memory_order_relaxed);}
                    }
                };
                const auto finishPhase=[&](const auto& serial){
                    LocalPlicExecutionState after;const bool got=after.capture();
                    if(!got||after.error!=entry.error||!after.sameControl(entry))failed.store(1,std::memory_order_relaxed);
#pragma omp critical(lamp_local_plic_floating_flags)
                    {if(got){nativeFlags|=after.native._Fe_stat;csrFlags|=after.csr&0x3fu;}}
#pragma omp barrier
#pragma omp master
                    {
                        phase.native._Fe_stat|=nativeFlags;phase.csr|=csrFlags;phase.restore();
                        if(!failed.load(std::memory_order_relaxed))try{serial();}catch(...){failed.store(1,std::memory_order_relaxed);}
                        if(!phase.capture()||phase.error!=entry.error||!phase.sameControl(entry))failed.store(1,std::memory_order_relaxed);
                        nativeFlags=0;csrFlags=0;
                    }
#pragma omp barrier
                    if(captured)phase.restore();
                };
                rows(0,n,initializeCell);finishPhase(prepareFaces);
                const auto pin=out.pin?static_cast<std::size_t>(out.pin->contact):n;
                rows(0,pin,reconstructCell);
                finishPhase([&]{if(pin<n)reconstructCell(pin,out.report);});
                if(pin<n){rows(pin+1,n,reconstructCell);finishPhase([]{});}
                if(captured)saved.restore();
            }
            if(failed.load(std::memory_order_relaxed))throw std::runtime_error("optional local PLIC worker replay");
            phase.restore();
            for(const auto& r:reports){
                out.report.mixedCells+=r.mixedCells;out.report.pureLiquidCells+=r.pureLiquidCells;out.report.pureGasCells+=r.pureGasCells;
                out.report.maximumRelativeVolumeResidual=std::max(out.report.maximumRelativeVolumeResidual,r.maximumRelativeVolumeResidual);
                out.report.maximumRelativeCapacityResidual=std::max(out.report.maximumRelativeCapacityResidual,r.maximumRelativeCapacityResidual);
                auto& a=out.report.localPlicReuse;const auto& b=r.localPlicReuse;
                a.geometryCells+=b.geometryCells;a.nodeRings+=b.nodeRings;a.mixedReused+=b.mixedReused;
                a.mixedReconstructed+=b.mixedReconstructed;a.mixedDependencyChanged+=b.mixedDependencyChanged;
            }
            for(std::size_t ci=0;ci<n;++ci)sumCell(ci);
            completed=true;cleanup.arm();
        }catch(...){
            // A literal mixed row can retain its incoming normal. Restore all
            // unpublished row values, not merely the rows that threw.
            for(std::size_t ci=0;ci<n;++ci){out.cells[ci]=GeometricVofCell{};out.plicCells[ci]=GeometricVofPlicCell{};
                out.alpha[ci]=0;out.capacityDefects[ci]=FiniteVolume::Dyadic{};widths[ci]=0;
                if(!geometry)coldCenters[ci]={};
                if(out.localMixedStates)(*out.localMixedStates)[ci]=LocalPlicQuietState{};}
            out.report=savedReport;out.pin=savedPin;out.faces.clear();faceState={};
            if(!reuseSupport){for(auto& list:cellNodes)list.clear();for(auto& list:nodeCells)list.clear();}
            stockSum={};measureSum={};entry.restore();
        }
        if(completed){cleanup.state.restore();cleanup.armed=false;}
    }
#endif
    if(!completed){
        for(std::size_t ci=0;ci<n;++ci)initializeCell(ci,out.report);
        prepareFaces();
        for(std::size_t ci=0;ci<n;++ci){reconstructCell(ci,out.report);sumCell(ci);}
    }
    // Contact can precede its wet/dry neighbours in native order. Only now
    // are all target planes/measures available. Original ordered sums above
    // remain untouched, and no candidate escapes a failed qualification.
    if(out.pin)qualifyPrescribedPinWalls(out);
    out.report.inventory=stockSum.sum;out.report.measured=measureSum.sum;
    if(!reuseSupport)retainLocalPlicNativeSupport(out,cellNodes,nodeCells);
    retainLocalPlicGeometry(mesh,out,coldCenters,vertexState,faceState,predecessor);
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
    if(completed)cleanup.arm();
#endif
}
std::uint64_t freshIdentity()
{
    static std::atomic<std::uint64_t> next{1};auto value=next.load(std::memory_order_relaxed);
    do {require(value!=std::numeric_limits<std::uint64_t>::max(),"VOF identity space exhausted");}
    while(!next.compare_exchange_weak(value,value+1,std::memory_order_relaxed));return value;
}
}
struct GeometricVofPublication::Data : BuildData {};
GeometricVofPublication::GeometricVofPublication(std::shared_ptr<const Data> data):data_(std::move(data)){}
const GeometricVofPublication::Data& GeometricVofPublication::checked() const
{ require(static_cast<bool>(data_),"access to consumed publication");return *data_; }
namespace
{
using CheckpointWriter=LAMP::IO::Checkpoint::Writer;
using CheckpointReader=LAMP::IO::Checkpoint::Reader;
constexpr const char* checkpointTag="LAMP.LocalPlicPublication.2";
constexpr const char* checkpointLegacyTag="LAMP.LocalPlicPublication.1";
void checkpointPoint(CheckpointWriter& writer,Point p){writer.scalar(p.x);writer.scalar(p.y);}
Point checkpointPoint(CheckpointReader& reader){const auto x=reader.scalar<long double>(),y=reader.scalar<long double>();return {x,y};}
void checkpointVector(CheckpointWriter& writer,Vector2D p){writer.scalar(p.x);writer.scalar(p.y);}
Vector2D checkpointVector(CheckpointReader& reader){const auto x=reader.scalar<double>(),y=reader.scalar<double>();return {x,y};}
void checkpointCompleteDomain(const NativeLeafDomain& domain,std::size_t n)
{
    require(domain.cellIndexByRow().size()==n&&domain.cellRowByIndex().size()==n,"checkpoint requires the complete rectangular domain");
    for(std::size_t ci=0;ci<n;++ci)require(domain.cellRowByIndex()[ci]==static_cast<int>(ci)
        &&domain.cellIndexByRow()[ci]==static_cast<int>(ci),"checkpoint does not support selected-domain masks");
}
std::vector<long double> checkpointCertificate(const FiniteVolume::Dyadic& certificate)
{
    // A native QL+QG-C expansion has at most three native terms. Save the
    // owner's actual certificate, not a public caller-supplied residual.
    auto remainder=certificate;std::vector<long double> terms;
    for(unsigned k=0;k<3&&remainder.compare(FiniteVolume::Dyadic{})!=0;++k){
        const auto term=FiniteVolume::round(remainder).value;require(std::isfinite(term)&&term!=0,"checkpoint capacity certificate is not natively decomposable");
        terms.push_back(term);remainder.add(term,-1);
    }
    require(remainder.compare(FiniteVolume::Dyadic{})==0,"checkpoint capacity certificate exceeds its three-term format");return terms;
}
void checkpointReport(CheckpointWriter& w,const GeometricVofReport& r)
{
    for(auto value:{r.mixedCells,r.pureLiquidCells,r.pureGasCells,r.zeroCollars,r.coincidentFaces})w.scalar<std::uint64_t>(value);
    w.scalar(r.maximumVolumeIterations);w.scalar(r.maximumRelativeVolumeResidual);w.scalar(r.maximumRelativeCapacityResidual);w.scalar(r.inventory);w.scalar(r.measured);
}
GeometricVofReport checkpointReport(CheckpointReader& r)
{
    GeometricVofReport out;out.mixedCells=r.scalar<std::uint64_t>();out.pureLiquidCells=r.scalar<std::uint64_t>();out.pureGasCells=r.scalar<std::uint64_t>();
    out.zeroCollars=r.scalar<std::uint64_t>();out.coincidentFaces=r.scalar<std::uint64_t>();out.maximumVolumeIterations=r.scalar<unsigned>();
    out.maximumRelativeVolumeResidual=r.scalar<long double>();out.maximumRelativeCapacityResidual=r.scalar<long double>();
    out.inventory=r.scalar<long double>();out.measured=r.scalar<long double>();return out;
}
void checkpointFace(CheckpointWriter& w,const GeometricVofFace& f)
{
    w.scalar(f.nativeRow);w.scalar(f.negativeCell);w.scalar(f.positiveCell);
    for(auto value:f.vertices)w.scalar<std::uint64_t>(value);w.boolean(f.coincident);w.scalar(f.rangeCount);
    for(const auto& range:f.ranges){w.scalar(range.begin);w.scalar(range.end);w.boolean(range.minusLiquid);w.boolean(range.plusLiquid);
        for(auto value:range.vertices)w.scalar<std::uint64_t>(value);}
    w.scalar(f.liquidMinus);w.scalar(f.liquidPlus);
}
GeometricVofFace checkpointFace(CheckpointReader& r)
{
    GeometricVofFace f;f.nativeRow=r.scalar<int>();f.negativeCell=r.scalar<int>();f.positiveCell=r.scalar<int>();
    for(auto& value:f.vertices)value=r.scalar<std::uint64_t>();f.coincident=r.boolean();f.rangeCount=r.scalar<unsigned>();
    for(auto& range:f.ranges){range.begin=r.scalar<long double>();range.end=r.scalar<long double>();range.minusLiquid=r.boolean();range.plusLiquid=r.boolean();
        for(auto& value:range.vertices)value=r.scalar<std::uint64_t>();}
    f.liquidMinus=r.scalar<long double>();f.liquidPlus=r.scalar<long double>();return f;
}
void validateCheckpointImage(const BuildData& data,const QuadtreeMesh* mesh)
{
    const auto n=data.cells.size();checkpointCompleteDomain(data.domain,n);data.time.time();
    require(data.kind==GeometricVofGeometryKind::LocalPlic&&data.time.position!=GeometricVofTimeLayer::Position::Midpoint,
        "checkpoint requires actual LocalPlic Initial/End state");
    require(data.localConfig.localResolution==GeometricVofLocalResolution::VolumeOnly||
        (data.localConfig.localResolution==GeometricVofLocalResolution::NativeInterfaceSupport&&!data.domain.axisymmetric()),
        "checkpoint has an unsupported local resolution law");
    require(data.alpha.size()==n&&data.plicCells.size()==n&&data.capacityDefects.size()==n,"checkpoint native cell arrays differ");
    require(data.anchoredVertices.empty()&&data.scalar.empty()&&data.facets.empty()&&data.topology.facetComponent.empty()
        &&data.topology.junctionVertices.empty()&&data.topology.axisEndpoints.empty()&&data.topology.exteriorEndpoints.empty()
        &&data.topology.components==0,"checkpoint LocalPlic image contains foreign fan topology");
    require(std::isfinite(data.origin.x)&&std::isfinite(data.origin.y),"checkpoint frame origin is nonfinite");
    if(mesh){require(n==static_cast<std::size_t>(mesh->cellIndexNum())&&data.vertices.size()==static_cast<std::size_t>(mesh->vertIndexNum()),
            "checkpoint cell/vertex dimensions differ from restored mesh");
        for(std::size_t v=0;v<data.vertices.size();++v){const auto point=mesh->verts(mesh->vertIndexToLabel[v]).coord();
            require(data.vertices[v].x==static_cast<long double>(point.x)-data.origin.x&&data.vertices[v].y==static_cast<long double>(point.y)-data.origin.y,
                "checkpoint native vertex differs from restored mesh");}}
    for(const auto& point:data.vertices)require(std::isfinite(point.x)&&std::isfinite(point.y),"checkpoint contains nonfinite native vertex");
    for(std::size_t ci=0;ci<n;++ci){const auto& p=data.plicCells[ci];const auto& cell=data.cells[ci];
        require(cell.nativeCell==static_cast<int>(ci)&&cell.capacity>0&&std::isfinite(cell.capacity)&&cell.triangles.empty(),"checkpoint invalid LocalPlic cell record");
        const auto normalLength=std::hypot(p.normal.x,p.normal.y);
        require(std::isfinite(p.normal.x)&&std::isfinite(p.normal.y)&&std::isfinite(p.offset)&&std::isfinite(normalLength)&&normalLength>0,
            "checkpoint nonfinite/zero PLIC plane");
        require(std::isfinite(data.alpha[ci])&&data.alpha[ci]>=0&&data.alpha[ci]<=1&&std::isfinite(cell.centerScalar),"checkpoint invalid derived cell scalar");
        FiniteVolume::Dyadic defect(p.inventory[0]);defect.add(p.inventory[1]);defect.add(cell.capacity,-1);
        require(defect.compare(data.capacityDefects[ci])==0&&defect.absolute().compare(FiniteVolume::Dyadic::product({1e-12L,cell.capacity}))<=0,
            "checkpoint exact capacity certificate mismatch");
        auto bound=FiniteVolume::Dyadic(cell.capacity);bound.add(defect.absolute());
        for(unsigned phase=0;phase<2;++phase)require(std::isfinite(p.inventory[phase])&&p.inventory[phase]>=0
            &&FiniteVolume::Dyadic(p.inventory[phase]).compare(bound)<=0&&std::isfinite(p.measured[phase])&&p.measured[phase]>=0
            &&std::abs(p.measured[phase]-p.inventory[phase])<=1e-12L*cell.capacity,"checkpoint invalid independent stock/measure");
        require(cell.inventory==p.inventory[0]&&cell.measured==p.measured[0],"checkpoint public cell/PLIC inventory differs");
        require(std::abs((p.measured[0]+p.measured[1])-cell.capacity)<=1e-12L*cell.capacity,"checkpoint measured phase closure differs");
        for(auto vertex:cell.boundaryVertices)require(vertex<data.vertices.size(),"checkpoint cell vertex outside native layout");
        if(mesh){const auto& native=mesh->cells(mesh->cellIndexToLabel[ci]);const auto lower=mesh->verts(native.cellToVert[0]).coord(),upper=mesh->verts(native.cellToVert[2]).coord();
            require(p.rectangle.lower.x==lower.x&&p.rectangle.lower.y==lower.y&&p.rectangle.upper.x==upper.x&&p.rectangle.upper.y==upper.y,
                "checkpoint PLIC rectangle differs from restored native cell");
            require(std::abs(cell.capacity-data.domain.cellVolumes()[ci])<=1e-12L*cell.capacity,"checkpoint capacity differs from restored native domain");}
    }
    require(data.faces.size()==static_cast<std::size_t>(data.domain.nativeFaceCount()),"checkpoint native faces differ");
    for(std::size_t row=0;row<data.faces.size();++row){const auto& f=data.faces[row];
        require(f.nativeRow==static_cast<int>(row)&&f.vertices[0]<data.vertices.size()&&f.vertices[1]<data.vertices.size()
            &&f.rangeCount==0&&!f.coincident&&f.liquidMinus==0&&f.liquidPlus==0,"checkpoint invalid LocalPlic native face");
        if(mesh){const auto actual=data.domain.face(*mesh,static_cast<int>(row));require(f.negativeCell==actual.negativeCell&&f.positiveCell==actual.positiveCell,
            "checkpoint native face adjacency differs from restored mesh");}}
    if(const auto& support=data.localNativeSupport){require(support->cells==n&&support->vertices==data.vertices.size()&&support->faces==data.faces.size()
            &&support->radial==data.domain.axisymmetric()&&support->cellNodes.size()==n&&support->nodeCells.size()==data.vertices.size()
            &&support->bytes<=32u*1024u*1024u,"checkpoint invalid native support dimensions");
        for(const auto& list:support->cellNodes)for(auto v:list)require(v<data.vertices.size(),"checkpoint support vertex outside native layout");
        for(const auto& list:support->nodeCells)for(int ci:list)require(ci>=0&&static_cast<std::size_t>(ci)<n,"checkpoint support cell outside native layout");}
    require(data.report.mixedCells<=n&&data.report.pureLiquidCells<=n&&data.report.pureGasCells<=n
        &&data.report.mixedCells+data.report.pureLiquidCells+data.report.pureGasCells==n,"checkpoint report cell counts differ");
}
}
std::string GeometricVofPublication::checkpoint() const
{
    const auto& data=checked();require(!data.pin,"pinned selected-domain checkpoint is not supported by the complete-domain image");
    validateCheckpointImage(data,nullptr);
    // Preserve the original default-law image, including its exact version-1
    // bytes. Only the explicit new model needs the version-2 configuration.
    const bool revised=data.localConfig.localResolution!=GeometricVofLocalResolution::VolumeOnly;
    CheckpointWriter w(revised?checkpointTag:checkpointLegacyTag);
    if(revised){w.scalar(data.localConfig.relativeVolumeTolerance);w.scalar(data.localConfig.maximumVolumeIterations);
        w.scalar(data.localConfig.maximumBracketExpansions);w.scalar(data.localConfig.localResolution);}
    w.values(data.domain.meshIdentity()->checkpointWords());w.boolean(data.domain.axisymmetric());
    w.values(data.domain.cellRowByIndex());w.values(data.domain.cellIndexByRow());w.values(data.domain.cellVolumes());
    w.scalar(data.time.step);w.scalar(data.time.begin);w.scalar(data.time.end);w.scalar(data.time.position);checkpointPoint(w,data.origin);w.values(data.alpha);
    w.count(data.vertices.size());for(auto p:data.vertices)checkpointPoint(w,p);
    w.count(data.cells.size());for(std::size_t ci=0;ci<data.cells.size();++ci){const auto& p=data.plicCells[ci];const auto& cell=data.cells[ci];
        checkpointVector(w,p.rectangle.lower);checkpointVector(w,p.rectangle.upper);checkpointVector(w,p.normal);w.scalar(p.offset);
        for(auto value:p.inventory)w.scalar(value);for(auto value:p.measured)w.scalar(value);w.boolean(p.resolved);
        w.values(checkpointCertificate(data.capacityDefects[ci]));w.scalar(cell.nativeCell);w.scalar(cell.capacity);w.scalar(cell.inventory);w.scalar(cell.measured);
        w.scalar(cell.centerScalar);w.values(cell.boundaryVertices);
    }
    w.count(data.faces.size());for(const auto& face:data.faces)checkpointFace(w,face);checkpointReport(w,data.report);
    w.boolean(static_cast<bool>(data.localNativeSupport));if(const auto& support=data.localNativeSupport){w.scalar<std::uint64_t>(support->bytes);
        w.count(support->cellNodes.size());for(const auto& list:support->cellNodes)w.values(list);
        w.count(support->nodeCells.size());for(const auto& list:support->nodeCells)w.values(list);}
    return w.finish();
}
GeometricVofPublication GeometricVofPublication::restoreCheckpoint(const QuadtreeMesh& mesh,NativeLeafDomain domain,const std::string& bytes)
{
    domain.validateMeshIdentity(mesh);checkpointCompleteDomain(domain,static_cast<std::size_t>(mesh.cellIndexNum()));
    const auto tagBytes=std::strlen(checkpointLegacyTag);
    const bool legacy=bytes.size()>=sizeof(std::uint64_t)+tagBytes&&
        bytes.compare(sizeof(std::uint64_t),tagBytes,checkpointLegacyTag)==0;
    CheckpointReader r(bytes,legacy?checkpointLegacyTag:checkpointTag);GeometricVofConfig restoredConfig;
    if(!legacy){restoredConfig.relativeVolumeTolerance=r.scalar<long double>();
        restoredConfig.maximumVolumeIterations=r.scalar<unsigned>();restoredConfig.maximumBracketExpansions=r.scalar<unsigned>();
        restoredConfig.localResolution=r.scalar<GeometricVofLocalResolution>();
        require(restoredConfig.localResolution==GeometricVofLocalResolution::VolumeOnly||
            restoredConfig.localResolution==GeometricVofLocalResolution::NativeInterfaceSupport,"checkpoint resolution law is retired or unsupported");
        require(std::isfinite(restoredConfig.relativeVolumeTolerance)&&restoredConfig.relativeVolumeTolerance>0&&
            restoredConfig.relativeVolumeTolerance<=1e-12L&&restoredConfig.maximumVolumeIterations>0&&restoredConfig.maximumVolumeIterations<=256&&
            restoredConfig.maximumBracketExpansions>0&&restoredConfig.maximumBracketExpansions<=256,"checkpoint has invalid local reconstruction configuration");}
    require(r.values<std::uint64_t>()==domain.meshIdentity()->checkpointWords(),"checkpoint belongs to a different actual native mesh/indexing");
    require(r.boolean()==domain.axisymmetric()&&r.values<int>()==domain.cellRowByIndex()&&r.values<int>()==domain.cellIndexByRow()
        &&r.values<double>()==domain.cellVolumes(),"checkpoint physical domain/measure differs");
    auto data=std::make_shared<Data>();data->kind=GeometricVofGeometryKind::LocalPlic;data->domain=std::move(domain);data->localConfig=restoredConfig;
    data->time.step=r.scalar<std::int64_t>();data->time.begin=r.scalar<double>();data->time.end=r.scalar<double>();data->time.position=r.scalar<GeometricVofTimeLayer::Position>();
    data->origin=checkpointPoint(r);data->alpha=r.values<double>();const auto vertices=r.count(2*sizeof(long double));
    require(vertices==static_cast<std::size_t>(mesh.vertIndexNum()),"checkpoint native vertex count differs");data->vertices.reserve(vertices);
    for(std::size_t v=0;v<vertices;++v)data->vertices.push_back(checkpointPoint(r));const auto cells=r.count();
    require(cells==static_cast<std::size_t>(mesh.cellIndexNum()),"checkpoint native cell count differs");data->cells.resize(cells);data->plicCells.resize(cells);data->capacityDefects.resize(cells);
    for(std::size_t ci=0;ci<cells;++ci){auto& p=data->plicCells[ci];auto& cell=data->cells[ci];
        p.rectangle.lower=checkpointVector(r);p.rectangle.upper=checkpointVector(r);p.normal=checkpointVector(r);p.offset=r.scalar<long double>();
        for(auto& value:p.inventory)value=r.scalar<long double>();for(auto& value:p.measured)value=r.scalar<long double>();p.resolved=r.boolean();
        const auto terms=r.values<long double>();require(terms.size()<=3,"checkpoint capacity certificate exceeds format");
        for(auto value:terms){require(std::isfinite(value),"checkpoint nonfinite certificate term");data->capacityDefects[ci].add(value);}
        cell.nativeCell=r.scalar<int>();cell.capacity=r.scalar<long double>();cell.inventory=r.scalar<long double>();cell.measured=r.scalar<long double>();
        cell.centerScalar=r.scalar<long double>();cell.boundaryVertices=r.values<std::size_t>();
    }
    const auto faces=r.count();require(faces==static_cast<std::size_t>(data->domain.nativeFaceCount()),"checkpoint native face count differs");data->faces.reserve(faces);
    for(std::size_t row=0;row<faces;++row)data->faces.push_back(checkpointFace(r));data->report=checkpointReport(r);
    if(r.boolean()){auto support=std::make_shared<LocalPlicNativeSupport>();support->identity=data->domain.meshIdentity();support->radial=data->domain.axisymmetric();
        support->cells=cells;support->vertices=vertices;support->faces=faces;support->bytes=r.scalar<std::uint64_t>();
        const auto cellLists=r.count();require(cellLists==cells,"checkpoint support cell count differs");support->cellNodes.resize(cells);
        for(auto& list:support->cellNodes)list=r.values<std::size_t>();const auto nodeLists=r.count();require(nodeLists==vertices,"checkpoint support vertex count differs");support->nodeCells.resize(vertices);
        for(auto& list:support->nodeCells)list=r.values<int>();data->localNativeSupport=std::move(support);
    }
    r.finish();validateCheckpointImage(*data,&mesh);data->identity=GeometricVofIdentity(freshIdentity());return GeometricVofPublication(std::move(data));
}
GeometricVofPublication::ConstructionMeshReceipt::ConstructionMeshReceipt(
    const QuadtreeMesh& mesh,const GeometricVofPublication& parent)
    :mesh_(&mesh),parentData_(parent.data_.get()),identity_(parent.domain().meshIdentity()){}
GeometricVofPublication::MeshTransferReceipt::MeshTransferReceipt(
    const QuadtreeMesh& sourceMesh,const GeometricVofPublication& parent,
    const QuadtreeMesh& targetMesh,const Mesh::NativeLeafOverlap& coverage)
    :sourceMesh_(&sourceMesh),targetMesh_(&targetMesh),parentData_(parent.data_.get()),coverage_(&coverage){}
GeometricVofPublication GeometricVofPublication::prepare(const QuadtreeMesh& mesh,NativeLeafDomain domain,
    GeometricVofTimeLayer time,const std::vector<double>& alpha,const GeometricVofConfig& config)
{
    require(config.localResolution==GeometricVofLocalResolution::VolumeOnly,"native interface support requires explicit LocalPlic geometry");
    require(std::isfinite(config.relativeVolumeTolerance)&&config.relativeVolumeTolerance>0&&config.relativeVolumeTolerance<=1e-12L,
        "volume tolerance must be finite, positive and no looser than 1e-12");
    require(config.maximumVolumeIterations>0&&config.maximumVolumeIterations<=256&&config.maximumBracketExpansions>0&&config.maximumBracketExpansions<=256,"invalid finite root budget");
    time.time();domain.validateMeshIdentity(mesh);
    require(alpha.size()==static_cast<std::size_t>(mesh.cellIndexNum()),"alpha must cover all actual native leaf cells");
    require(domain.cellIndexByRow().size()==alpha.size(),"VOF requires the complete actual leaf domain");
    for(double a:alpha)require(std::isfinite(a)&&a>=0&&a<=1,"alpha outside [0,1]; no clipping");
    auto data=std::make_shared<Data>();data->domain=std::move(domain);data->time=time;data->alpha=alpha;
    // One exposed frame for every native site, root, fan and consumer. The
    // radial origin stays on the lower wall; the axial midpoint reduces the
    // maximum axial coordinate magnitude without an alpha-dependent choice.
    data->origin={mesh.minCoord.x,static_cast<long double>(mesh.minCoord.y)
        +(static_cast<long double>(mesh.maxCoord.y)-mesh.minCoord.y)/2};
    Builder builder{mesh,*data,config,{},{},{},{},{},data->domain.axisymmetric()};
    builder.layout();builder.planes();builder.sharedScalars();builder.reconstruct();builder.gamma();builder.faceTraces();builder.connections();
    data->identity=GeometricVofIdentity(freshIdentity());return GeometricVofPublication(std::move(data));
}
GeometricVofPublication GeometricVofPublication::prepareLocalPlic(const QuadtreeMesh& mesh,NativeLeafDomain domain,
    GeometricVofTimeLayer time,const std::vector<double>& alpha,const GeometricVofConfig& config)
{
    time.time();domain.validateMeshIdentity(mesh);
    require(alpha.size()==static_cast<std::size_t>(mesh.cellIndexNum())&&domain.cellRowByIndex().size()==alpha.size(),
        "local PLIC alpha must cover actual native indices, including zero inactive entries");
    std::vector<std::array<long double,2>> inventory;inventory.reserve(alpha.size());
    for(std::size_t ci=0;ci<alpha.size();++ci){const auto value=alpha[ci];
        require(std::isfinite(value)&&value>=0&&value<=1,"local PLIC alpha outside [0,1]; no clipping");
        if(domain.cellRowByIndex()[ci]<0){
            require(value==0,"inactive local PLIC alpha must be exactly zero");inventory.push_back({});continue;
        }
        const auto& cell=mesh.cells(mesh.cellIndexToLabel[ci]);
        const Cut::Rectangle rectangle{mesh.verts(cell.cellToVert[0]).coord(),mesh.verts(cell.cellToVert[2]).coord()};
        const auto capacity=Cut::capacity(rectangle,domain.axisymmetric());const auto a=static_cast<long double>(value);
        inventory.push_back({a*capacity,(1-a)*capacity});}
    return prepareLocalPlicInventory(mesh,std::move(domain),time,inventory,config);
}
GeometricVofPublication GeometricVofPublication::prepareLocalPlicInventory(const QuadtreeMesh& mesh,NativeLeafDomain domain,
    GeometricVofTimeLayer time,const std::vector<std::array<long double,2>>& inventory,const GeometricVofConfig& config)
{
    require(std::isfinite(config.relativeVolumeTolerance)&&config.relativeVolumeTolerance>0&&config.relativeVolumeTolerance<=1e-12L,
        "local PLIC volume tolerance must be finite, positive and no looser than 1e-12");
    time.time();domain.validateMeshIdentity(mesh);
    require(inventory.size()==static_cast<std::size_t>(mesh.cellIndexNum())&&domain.cellRowByIndex().size()==inventory.size(),
        "local PLIC inventory must cover actual native indices, including zero inactive entries");
    auto data=std::make_shared<Data>();data->domain=std::move(domain);data->time=time;
    data->origin={mesh.minCoord.x,static_cast<long double>(mesh.minCoord.y)
        +(static_cast<long double>(mesh.maxCoord.y)-mesh.minCoord.y)/2};
    buildLocalPlic(mesh,*data,inventory,config);data->identity=GeometricVofIdentity(freshIdentity());
    return GeometricVofPublication(std::move(data));
}
GeometricVofPublication GeometricVofPublication::prepareCertifiedLocalPlicInventory(
    const QuadtreeMesh& mesh,NativeLeafDomain domain,GeometricVofTimeLayer time,
    const GeometricVofPublication& parent,CertifiedOperation operation,
    const std::vector<std::array<long double,2>>& inventory,
    const std::vector<FiniteVolume::Dyadic>& predictedDefects,const GeometricVofConfig& config,
    const ConstructionMeshReceipt* construction,const MeshTransferReceipt* transfer)
{
    require(std::isfinite(config.relativeVolumeTolerance)&&config.relativeVolumeTolerance>0&&config.relativeVolumeTolerance<=1e-12L,
        "derived local PLIC volume tolerance must be finite, positive and no looser than 1e-12");
    require(config.maximumVolumeIterations>0&&config.maximumVolumeIterations<=256&&
        config.maximumBracketExpansions>0&&config.maximumBracketExpansions<=256,"invalid finite root budget");
    time.time();
    // Only this same callback-free sweep can offer the receipt. Every other
    // caller/binding keeps the original domain-then-parent live validation.
    const bool meshAlreadyValidated=construction&&construction->mesh_==&mesh
        &&construction->parentData_==parent.data_.get()&&parent.data_&&construction->identity_
        &&domain.meshIdentity()==construction->identity_
        &&parent.data_->domain.meshIdentity()==construction->identity_;
    if(!meshAlreadyValidated)domain.validateMeshIdentity(mesh);
    const auto& predecessor=parent.checked();predecessor.time.time();
    require(predecessor.kind==GeometricVofGeometryKind::LocalPlic&&predecessor.identity.value()!=0&&
        predecessor.capacityDefects.size()==predecessor.plicCells.size(),"invalid typed local PLIC predecessor");
    require(config.localResolution==predecessor.localConfig.localResolution,"derived LocalPlic cannot change its owner's resolution law");
    require(domain.axisymmetric()==predecessor.domain.axisymmetric(),"derived local PLIC geometry kind differs");
    if(operation==CertifiedOperation::MeshTransfer){
        const bool complete=domain.cellIndexByRow().size()==domain.cellRowByIndex().size()&&
            predecessor.domain.cellIndexByRow().size()==predecessor.domain.cellRowByIndex().size();
        require(complete||transfer,"selected-domain mesh transfer requires actual-overlap receipt");
        require(sameClock(time,predecessor.time),"mesh transfer closure requires the same full clock");
        if(!complete){
            require(transfer->sourceMesh_&&transfer->targetMesh_==&mesh&&
                transfer->parentData_==parent.data_.get()&&transfer->coverage_,
                "selected transfer receipt has foreign mesh or predecessor");
            require(time.position==GeometricVofTimeLayer::Position::End,
                "selected transfer requires the same accepted End clock");
            parent.validateFor(*transfer->sourceMesh_,predecessor.identity,time);
            const auto& coverage=*transfer->coverage_;coverage.validateFor(*transfer->sourceMesh_,mesh);
            const auto validateSelection=[&](const auto& active,const auto& groups){
                require(active.size()==groups.size(),"selected transfer coverage dimensions differ");
                for(std::size_t ci=0;ci<active.size();++ci)
                    require((active[ci]>=0)==!groups[ci].empty(),
                        "selected transfer coverage differs from actual fluid map");
            };
            validateSelection(predecessor.domain.cellRowByIndex(),coverage.bySource());
            validateSelection(domain.cellRowByIndex(),coverage.byTarget());
            require(coverage.maximumRelativeCoverageResidual()<=512*std::numeric_limits<long double>::epsilon(),
                "selected transfer lacks complete actual physical coverage");
        }
    }else if(operation==CertifiedOperation::DirectionalScratch){
        if(!meshAlreadyValidated)parent.validateMesh(mesh);require(sameClock(time,predecessor.time),"directional closure requires source full clock");
    }else{
        require(operation==CertifiedOperation::HalfTransport||operation==CertifiedOperation::FullTransport,
            "unknown typed local PLIC closure operation");if(!meshAlreadyValidated)parent.validateMesh(mesh);
        const auto position=operation==CertifiedOperation::HalfTransport?GeometricVofTimeLayer::Position::Midpoint:GeometricVofTimeLayer::Position::End;
        require(predecessor.time.position!=GeometricVofTimeLayer::Position::Midpoint&&
            predecessor.time.step<std::numeric_limits<std::int64_t>::max()&&time.step==predecessor.time.step+1&&
            time.begin==predecessor.time.time()&&time.position==position,
            "transport closure requires a contiguous genuine source clock");
    }
    if(operation!=CertifiedOperation::MeshTransfer)
        require(domain.cellRowByIndex()==predecessor.domain.cellRowByIndex()&&
            domain.cellIndexByRow()==predecessor.domain.cellIndexByRow()&&
            domain.nativeFaceCount()==predecessor.domain.nativeFaceCount(),
            "derived local PLIC selected-domain map differs from predecessor");
    require(inventory.size()==static_cast<std::size_t>(mesh.cellIndexNum())&&domain.cellRowByIndex().size()==inventory.size()&&
        predictedDefects.size()==inventory.size(),"derived local PLIC inventory/certificate must cover native indices");
    auto data=std::make_shared<Data>();data->domain=std::move(domain);data->time=time;
    data->origin={mesh.minCoord.x,static_cast<long double>(mesh.minCoord.y)
        +(static_cast<long double>(mesh.maxCoord.y)-mesh.minCoord.y)/2};
    // Only this actual typed parent can offer native support/geometry. The live
    // source/binding/clock/stock gates above remain unchanged. Exact selected
    // maps are checked before inheritance; node rings contain only active
    // cells. A new domain or MeshTransfer has no inherited preparation.
    if(operation!=CertifiedOperation::MeshTransfer&&predecessor.localNativeSupport
        &&data->domain.meshIdentity()==predecessor.domain.meshIdentity()
        &&data->domain.cellIndexByRow()==predecessor.domain.cellIndexByRow()
        &&data->domain.cellRowByIndex()==predecessor.domain.cellRowByIndex())
    {data->localNativeSupport=predecessor.localNativeSupport;
        data->localGeometry=predecessor.localGeometry;}
    if(predecessor.pin){
        // Same-End AMR uses the retained physical point and signed wall legs;
        // old native rows/cells never identify target support.
        data->pin=operation==CertifiedOperation::MeshTransfer?
            qualifyPinnedCorner(mesh,data->domain,predecessor.pin->point,-1,-1,&*predecessor.pin):predecessor.pin;
    }
    buildLocalPlic(mesh,*data,inventory,config,&predictedDefects,
        operation==CertifiedOperation::MeshTransfer?nullptr:&predecessor);data->identity=GeometricVofIdentity(freshIdentity());
    return GeometricVofPublication(std::move(data));
}
GeometricVofPublication GeometricVofPublication::withPinnedCorner(const QuadtreeMesh& mesh,
    const GeometricVofPublication& source,Vector2D point,int wet,int dry,PinCertificate::WallTraceLaw wallTraceLaw)
{
    source.validateMesh(mesh);const auto& predecessor=source.checked();
    require(predecessor.kind==GeometricVofGeometryKind::LocalPlic,"corner pin requires a LocalPlic source");
    require(predecessor.localConfig.localResolution==GeometricVofLocalResolution::VolumeOnly,
        "native interface support does not support pinned LocalPlic");
    require(!predecessor.pin,"an existing pinned corner cannot be replaced");
    require(wallTraceLaw==PinCertificate::WallTraceLaw::StrictWholeCell||wallTraceLaw==PinCertificate::WallTraceLaw::PrescribedWetDry,
        "invalid pinned wall trace law");
    auto data=std::make_shared<Data>(predecessor);
    data->pin=qualifyPinnedCorner(mesh,data->domain,point,wet,dry);data->pin->wallTraceLaw=wallTraceLaw;
    std::vector<std::array<long double,2>> inventory;inventory.reserve(data->plicCells.size());
    for(const auto& plic:data->plicCells)inventory.push_back(plic.inventory);
    reconstructPinnedContact(mesh,*data,inventory,data->localConfig);
    const auto& pinnedCell=data->plicCells[data->pin->contact];const auto pinnedCapacity=data->cells[data->pin->contact].capacity;
    data->report.maximumRelativeCapacityResidual=std::max(data->report.maximumRelativeCapacityResidual,
        std::abs(pinnedCell.measured[0]+pinnedCell.measured[1]-pinnedCapacity)/pinnedCapacity);
    data->report.maximumRelativeVolumeResidual=0;PositiveSum measured;
    for(int ci:data->domain.cellIndexByRow()){const auto& plic=data->plicCells[ci];const auto& cell=data->cells[ci];
        measured.add(cell.measured);for(unsigned phase=0;phase<2;++phase)
            data->report.maximumRelativeVolumeResidual=std::max(data->report.maximumRelativeVolumeResidual,
                std::abs(plic.measured[phase]-plic.inventory[phase])/cell.capacity);}
    data->report.measured=measured.sum;data->identity=GeometricVofIdentity(freshIdentity());
    return GeometricVofPublication(std::move(data));
}
struct GeometricVofPublication::PinCertificate::Data
{
    std::shared_ptr<const GeometricVofPublication::Data> owner;
};
GeometricVofPublication::PinCertificate::PinCertificate(std::shared_ptr<const Data> data):data_(std::move(data)){}
GeometricVofIdentity GeometricVofPublication::PinCertificate::sourceIdentity() const
{require(bool(data_)&&bool(data_->owner),"consumed pin certificate");return data_->owner->identity;}
GeometricVofTimeLayer GeometricVofPublication::PinCertificate::sourceTime() const
{require(bool(data_)&&bool(data_->owner),"consumed pin certificate");return data_->owner->time;}
Vector2D GeometricVofPublication::PinCertificate::physicalPoint() const
{sourceIdentity();return data_->owner->pin->point;}
int GeometricVofPublication::PinCertificate::nativeVertex() const
{sourceIdentity();return data_->owner->pin->vertex;}
int GeometricVofPublication::PinCertificate::nativeCell() const
{sourceIdentity();return data_->owner->pin->contact;}
Cut::RectangleCorner GeometricVofPublication::PinCertificate::corner() const
{sourceIdentity();return data_->owner->pin->corner;}
int GeometricVofPublication::PinCertificate::wetNativeRow() const
{sourceIdentity();return data_->owner->pin->wet;}
int GeometricVofPublication::PinCertificate::dryNativeRow() const
{sourceIdentity();return data_->owner->pin->dry;}
GeometricVofPublication::PinCertificate::Kind GeometricVofPublication::PinCertificate::kind() const
{sourceIdentity();return data_->owner->pin->kind;}
GeometricVofPublication::PinCertificate::WallTraceLaw GeometricVofPublication::PinCertificate::wallTraceLaw() const
{sourceIdentity();return data_->owner->pin->wallTraceLaw;}
bool GeometricVofPublication::PinCertificate::prescribedWetWallTrace() const
{sourceIdentity();return data_->owner->pin->prescribedWet;}
bool GeometricVofPublication::PinCertificate::prescribedDryWallTrace() const
{sourceIdentity();return data_->owner->pin->prescribedDry;}
void GeometricVofPublication::PinCertificate::validateFor(const QuadtreeMesh& mesh,const GeometricVofPublication& source) const
{
    require(bool(data_)&&bool(data_->owner)&&bool(data_->owner->pin),"consumed pin certificate");
    source.validateFor(mesh,data_->owner->identity,data_->owner->time);
    require(source.data_==data_->owner,"different immutable pin certificate owner");
}
void GeometricVofPublication::PinCertificate::validateForInBatch(const QuadtreeMesh& mesh,
    const GeometricVofPublication& source,const MeshValidationBatch& batch) const
{
    require(bool(data_)&&bool(data_->owner)&&bool(data_->owner->pin),"consumed pin certificate");
    source.validateForInBatch(mesh,data_->owner->identity,data_->owner->time,batch);
    require(source.data_==data_->owner,"different immutable pin certificate owner");
}
std::optional<GeometricVofPublication::PinCertificate> GeometricVofPublication::pinnedCorner() const
{
    const auto& data=checked();if(!data.pin)return std::nullopt;
    auto certificate=std::make_shared<PinCertificate::Data>();certificate->owner=data_;
    return PinCertificate(std::move(certificate));
}
GeometricVofGeometryKind GeometricVofPublication::geometryKind() const{return checked().kind;}
GeometricVofLocalResolution GeometricVofPublication::localPlicResolution() const
{const auto& data=checked();require(data.kind==GeometricVofGeometryKind::LocalPlic,"resolution law requires LocalPlic");return data.localConfig.localResolution;}
const std::vector<GeometricVofPlicCell>& GeometricVofPublication::plicCells() const
{const auto& data=checked();require(data.kind==GeometricVofGeometryKind::LocalPlic,"PLIC cells requested from compatible fan");return data.plicCells;}
std::vector<bool> GeometricVofPublication::localPlicSubgridMask(const QuadtreeMesh& mesh) const
{return localPlicResolutionMasks(mesh).capillarySubgrid;}
GeometricVofLocalResolutionMasks GeometricVofPublication::localPlicResolutionMasks(const QuadtreeMesh& mesh) const
{
    validateMesh(mesh);const auto& data=checked();
    require(data.kind==GeometricVofGeometryKind::LocalPlic,"subgrid resolution mask requires local PLIC");
    const auto n=data.plicCells.size();
    GeometricVofLocalResolutionMasks masks{std::vector<bool>(n,false),std::vector<bool>(n,false)};
    auto& result=masks.capillarySubgrid;
    if(data.localConfig.localResolution==GeometricVofLocalResolution::VolumeOnly)return masks;
    require(data.localConfig.localResolution==GeometricVofLocalResolution::NativeInterfaceSupport&&
        !data.domain.axisymmetric()&&!data.pin,"native interface support requires unpinned XY LocalPlic");
    const auto& active=data.domain.cellRowByIndex();const auto support=data.localNativeSupport;
    const bool reuseSupport=support&&support->identity==data.domain.meshIdentity()&&!support->radial&&support->cells==n&&
        support->vertices==static_cast<std::size_t>(mesh.vertIndexNum())&&support->faces==static_cast<std::size_t>(data.domain.nativeFaceCount())&&
        support->cellNodes.size()==n&&support->nodeCells.size()==support->vertices;
    // A cold/restored owner constructs its actual incidence once, never by a
    // global search per cell. The mask itself is always freshly derived from
    // both raw stocks, so no alpha-only or predecessor classification is used.
    std::vector<std::vector<std::size_t>> cellNodes(reuseSupport?0:n);
    std::vector<std::vector<int>> nodeCells(reuseSupport?0:static_cast<std::size_t>(mesh.vertIndexNum()));
    if(!reuseSupport){
        for(const auto& face:data.faces)for(int ci:{face.negativeCell,face.positiveCell})if(ci>=0)
            for(auto vertex:face.vertices){require(static_cast<std::size_t>(ci)<n&&vertex<nodeCells.size(),
                    "subgrid native incidence outside publication");
                cellNodes[ci].push_back(vertex);nodeCells[vertex].push_back(ci);}
        for(auto& list:cellNodes){std::sort(list.begin(),list.end());list.erase(std::unique(list.begin(),list.end()),list.end());}
        for(auto& list:nodeCells){std::sort(list.begin(),list.end());list.erase(std::unique(list.begin(),list.end()),list.end());}
    }
    const auto& actualCellNodes=reuseSupport?support->cellNodes:cellNodes;
    const auto& actualNodeCells=reuseSupport?support->nodeCells:nodeCells;
    // Eligibility is distinct from a false base mask: incomplete wall/solid
    // support must never become a valid endpoint for the orientation law.
    std::vector<bool> eligible(n,false),base(n,false);
    for(std::size_t ci=0;ci<n;++ci){const auto& source=data.plicCells[ci];
        if(active[ci]<0||!(source.inventory[0]>0&&source.inventory[1]>0))continue;
        const auto& vertices=actualCellNodes[ci];
        if(vertices.size()<4||vertices.size()>8)continue;
        std::array<bool,2> majority{};bool complete=true;
        for(auto vertex:vertices){
            if(vertex>=actualNodeCells.size()){complete=false;break;}
            const auto point=mesh.verts(mesh.vertIndexToLabel[vertex]).coord();unsigned coverage=0;
            const auto& incident=actualNodeCells[vertex];
            if(incident.empty()||incident.size()>8){complete=false;break;}
            for(int native:incident){
                if(native<0||static_cast<std::size_t>(native)>=n||active[native]<0)continue;
                const auto& cell=data.plicCells[native];const auto& box=cell.rectangle;
                if(!(box.lower.x<=point.x&&point.x<=box.upper.x&&box.lower.y<=point.y&&point.y<=box.upper.y)){
                    complete=false;break;}
                for(unsigned y=0;y<2;++y)for(unsigned x=0;x<2;++x)
                    if((x?box.upper.x>point.x:box.lower.x<point.x)&&
                        (y?box.upper.y>point.y:box.lower.y<point.y))coverage|=1u<<(2*y+x);
                majority[0]=majority[0]||cell.inventory[0]>=cell.inventory[1];
                majority[1]=majority[1]||cell.inventory[1]>=cell.inventory[0];
            }
            if(!complete||coverage!=15u){complete=false;break;}
        }
        if(complete){eligible[ci]=true;
            // The declared capillary model also covers an unresolved local
            // orientation; a volume-matching plane alone cannot resolve it.
            base[ci]=!majority[0]||!majority[1]||!source.resolved;
        }
    }
    result=base;
    const auto resolvedMixed=[&](std::size_t native){const auto& cell=data.plicCells[native];
        return cell.resolved&&cell.measured[0]>0&&cell.measured[1]>0&&
            std::isfinite(cell.normal.x)&&std::isfinite(cell.normal.y)&&(cell.normal.x!=0||cell.normal.y!=0);};
    // The complete base decision above is immutable throughout this pass.
    // Only two eligible, actually measured mixed interfaces can expose an
    // opposite normal in their shared native vertex star. Classification is
    // simultaneous and never propagates through a newly marked neighbour.
    for(std::size_t ci=0;ci<n;++ci){
        if(!eligible[ci]||base[ci]||!resolvedMixed(ci))continue;
        const auto normal=data.plicCells[ci].normal;
        for(auto vertex:actualCellNodes[ci]){
            for(int native:actualNodeCells[vertex]){
                if(native<0||static_cast<std::size_t>(native)>=n||native==static_cast<int>(ci)||
                    !eligible[native]||base[native]||!resolvedMixed(static_cast<std::size_t>(native)))continue;
                const auto other=data.plicCells[native].normal;
                const auto dot=static_cast<long double>(normal.x)*other.x+static_cast<long double>(normal.y)*other.y;
                if(std::isfinite(dot)&&dot<=0){result[ci]=true;break;}
            }
            if(result[ci])break;
        }
    }
    // Freeze the original two passes. Resolved witness paths below never
    // depend on marks made by another source in this third pass.
    masks.geometryUnsupported=result;const auto& orientation=masks.geometryUnsupported;
    std::vector<std::array<bool,2>> carrier(n);
    for(std::size_t ci=0;ci<n;++ci)if(active[ci]>=0){const auto& cell=data.plicCells[ci];
        const bool resolved=eligible[ci]&&!orientation[ci]&&resolvedMixed(ci);
        for(unsigned phase=0;phase<2;++phase)
            carrier[ci][phase]=cell.measured[phase]>0&&(cell.measured[1-phase]==0||resolved);
    }
    struct PhaseInterval {long double lower=0,upper=0;};
    const auto traceIntervals=[&](int native,Point a,Point b,bool horizontal){
        std::array<PhaseInterval,2> intervals{};const auto& cell=data.plicCells[native];
        const long double extent=horizontal?b.x-a.x:b.y-a.y;
        require(std::isfinite(extent)&&extent>0,"invalid native witness face interval");
        for(unsigned phase=0;phase<2;++phase)if(carrier[native][phase]&&cell.measured[1-phase]==0){
            intervals[phase]={0,extent};return intervals;}
        if(!carrier[native][0]||!carrier[native][1])return intervals;
        const auto length=std::hypot(static_cast<long double>(cell.normal.x),static_cast<long double>(cell.normal.y));
        require(std::isfinite(length)&&length>0,"invalid native witness plane normal");
        const long double nx=cell.normal.x/length,ny=cell.normal.y/length;
        const long double cx=static_cast<long double>(cell.rectangle.lower.x)+
            (static_cast<long double>(cell.rectangle.upper.x)-cell.rectangle.lower.x)/2;
        const long double cy=static_cast<long double>(cell.rectangle.lower.y)+
            (static_cast<long double>(cell.rectangle.upper.y)-cell.rectangle.lower.y)/2;
        const auto value=[&](Point point){return nx*(point.x-cx)+ny*(point.y-cy)-cell.offset;};
        const auto va=value(a),vb=value(b);
        require(std::isfinite(va)&&std::isfinite(vb),"nonfinite native witness plane trace");
        if(va==0&&vb==0)return intervals; // No strict phase is certified on a coincident mixed face.
        if(va<=0&&vb<=0){intervals[0]={0,extent};return intervals;}
        if(va>=0&&vb>=0){intervals[1]={0,extent};return intervals;}
        // A common face-local coordinate avoids rounding a small occupied
        // interval away when adding a large world-coordinate origin.
        long double crossing=0;
        try{crossing=Cut::linearZeroCrossing({0,0},{extent,0},va,vb).x;}
        catch(const std::runtime_error& error){
            if(std::string(error.what())!="[PlanarCellCut] line/edge intersection is unrepresentable")throw;
            return intervals; // No representable positive witness link, never a phase-measure certificate.
        }
        if(va<0){intervals[0]={0,crossing};intervals[1]={crossing,extent};}
        else{intervals[0]={crossing,extent};intervals[1]={0,crossing};}
        return intervals;
    };
    // Fresh phase adjacency is prepared once per publication query. LocalPlic
    // face trace fields are placeholders, so use the actual planes/subfaces.
    std::vector<std::array<std::vector<int>,2>> phaseNeighbors(n);
    for(const auto& face:data.faces){const int a=face.negativeCell,b=face.positiveCell;
        if(a<0||b<0||active[a]<0||active[b]<0)continue;
        const auto& ca=data.plicCells[a];const auto& cb=data.plicCells[b];
        const bool mixedA=ca.measured[0]>0&&ca.measured[1]>0,mixedB=cb.measured[0]>0&&cb.measured[1]>0;
        if(!mixedA&&!mixedB)continue; // A pure-majority target already terminates its path.
        const auto pa=mesh.verts(mesh.vertIndexToLabel[static_cast<int>(face.vertices[0])]).coord();
        const auto pb=mesh.verts(mesh.vertIndexToLabel[static_cast<int>(face.vertices[1])]).coord();
        const Point first{pa.x,pa.y},last{pb.x,pb.y};const bool horizontal=first.y==last.y;
        const auto minus=traceIntervals(a,first,last,horizontal),plus=traceIntervals(b,first,last,horizontal);
        for(unsigned phase=0;phase<2;++phase)
            if(std::max(minus[phase].lower,plus[phase].lower)<std::min(minus[phase].upper,plus[phase].upper)){
                phaseNeighbors[a][phase].push_back(b);phaseNeighbors[b][phase].push_back(a);}
    }
    for(std::size_t ci=0;ci<n;++ci){
        if(!eligible[ci]||orientation[ci]||!resolvedMixed(ci))continue;
        std::array<int,64> star{};std::size_t count=0;
        for(auto vertex:actualCellNodes[ci])for(int native:actualNodeCells[vertex])
            if(native>=0&&static_cast<std::size_t>(native)<n&&active[native]>=0){
            require(count<star.size(),"native witness star exceeds bounded incidence");star[count++]=native;}
        std::sort(star.begin(),star.begin()+count);
        count=static_cast<std::size_t>(std::unique(star.begin(),star.begin()+count)-star.begin());
        const auto end=star.begin()+count,source=std::lower_bound(star.begin(),end,static_cast<int>(ci));
        require(source!=end&&*source==static_cast<int>(ci),"native witness star omits source");
        for(unsigned phase=0;phase<2;++phase){
            std::array<bool,64> visited{};std::array<std::size_t,64> pending{};std::size_t first=0,last=0;
            const auto sourceIndex=static_cast<std::size_t>(source-star.begin());visited[sourceIndex]=true;pending[last++]=sourceIndex;
            bool witnessed=false;
            while(first<last&&!witnessed){const int native=star[pending[first++]];
                for(int neighbor:phaseNeighbors[native][phase]){const auto found=std::lower_bound(star.begin(),end,neighbor);
                    if(found==end||*found!=neighbor)continue;const auto index=static_cast<std::size_t>(found-star.begin());
                    if(visited[index])continue;visited[index]=true;
                    const auto& cell=data.plicCells[neighbor];
                    if(neighbor!=static_cast<int>(ci)&&cell.inventory[phase]>=cell.inventory[1-phase]){witnessed=true;break;}
                    require(last<pending.size(),"native witness queue exceeds bounded star");pending[last++]=index;
                }
            }
            if(!witnessed){result[ci]=true;break;}
        }
    }
    return masks;
}
const std::vector<FiniteVolume::Dyadic>& GeometricVofPublication::capacityDefects() const
{const auto& data=checked();require(data.kind==GeometricVofGeometryKind::LocalPlic&&data.capacityDefects.size()==data.plicCells.size(),
    "local PLIC capacity certificates unavailable");return data.capacityDefects;}
GeometricVofIdentity GeometricVofPublication::identity() const{return checked().identity;}
GeometricVofTimeLayer GeometricVofPublication::timeLayer() const{return checked().time;}
const NativeLeafDomain& GeometricVofPublication::domain() const{return checked().domain;}
const std::vector<double>& GeometricVofPublication::alpha() const{return checked().alpha;}
Point GeometricVofPublication::frameOrigin() const{return checked().origin;}
const std::vector<AnchoredPoint>& GeometricVofPublication::anchoredVertices() const{const auto& data=checked();require(data.kind==GeometricVofGeometryKind::CompatibleFan,"anchoredVertices unavailable for local PLIC");return data.anchoredVertices;}
const std::vector<Point>& GeometricVofPublication::vertices() const{return checked().vertices;}
const std::vector<long double>& GeometricVofPublication::scalarValues() const{const auto& data=checked();require(data.kind==GeometricVofGeometryKind::CompatibleFan,"scalarValues unavailable for local PLIC");return data.scalar;}
const std::vector<GeometricVofCell>& GeometricVofPublication::cells() const{return checked().cells;}
const std::vector<GeometricVofFace>& GeometricVofPublication::faces() const{return checked().faces;}
const std::vector<GeometricVofFacet>& GeometricVofPublication::facets() const{const auto& data=checked();require(data.kind==GeometricVofGeometryKind::CompatibleFan,"facets unavailable for local PLIC");return data.facets;}
const GeometricVofTopology& GeometricVofPublication::topology() const{const auto& data=checked();require(data.kind==GeometricVofGeometryKind::CompatibleFan,"global topology unavailable for local PLIC");return data.topology;}
const GeometricVofReport& GeometricVofPublication::report() const{return checked().report;}
std::vector<std::vector<AnchoredPoint>> GeometricVofPublication::anchoredPhasePolygons(int nativeCell,bool liquid) const
{
    const auto& data=checked();require(data.kind==GeometricVofGeometryKind::CompatibleFan,"fan phase polygons unavailable for local PLIC");require(nativeCell>=0&&static_cast<std::size_t>(nativeCell)<data.cells.size(),"native phase polygon cell outside publication");
    const auto& cell=data.cells[nativeCell];std::vector<std::vector<AnchoredPoint>> polygons;
    if(cell.triangles.empty()){
        if((data.alpha[nativeCell]==1)==liquid){std::vector<AnchoredPoint> polygon;polygon.reserve(cell.boundaryVertices.size());
            for(auto v:cell.boundaryVertices)polygon.push_back(data.anchoredVertices[v]);polygons.push_back(std::move(polygon));}
    }else for(const auto& triangle:cell.triangles){
        auto cut=triangle.anchoredCut;if(!liquid){auto scalar=triangle.scalar;for(auto& value:scalar)value=-value;cut=Cut::cutAnchoredTriangle(triangle.anchoredVertices,scalar);}
        if(cut.negativeCount>=3){std::vector<AnchoredPoint> polygon(cut.negativeVertices.begin(),cut.negativeVertices.begin()+cut.negativeCount);
            if(Cut::anchoredPolygonMoments(polygon).area>0)polygons.push_back(std::move(polygon));}
    }
    return polygons;
}
std::vector<std::vector<Point>> GeometricVofPublication::phasePolygons(int nativeCell,bool liquid) const
{
    std::vector<std::vector<Point>> polygons;for(const auto& actual:anchoredPhasePolygons(nativeCell,liquid)){
        std::vector<Point> polygon;polygon.reserve(actual.size());for(const auto p:actual)polygon.push_back(Cut::approximateAnchoredPoint(p));
        polygons.push_back(std::move(polygon));}
    return polygons;
}
// D22: read-only geometric sufficient condition. This does not reconstruct
// another interface, select a curvature, or change either phase inventory.
namespace
{
struct FragmentPlane
{
    enum class Kind { Empty,Full,Mixed };
    Kind kind=Kind::Empty;
    long double cx=0,cy=0,hx=0,hy=0,nx=0,ny=0,offset=0;
};
struct FragmentInterval {long double lower=0,upper=0;bool empty=true;};
struct FragmentBounds {long double x0=0,x1=0,y0=0,y1=0;bool set=false;};
struct FragmentProof {unsigned phase=0;long double diameter=0,delta=0;std::size_t faces=0,corners=0;};
std::optional<long double> fragmentExactDifference(long double a,long double b)
{
    FiniteVolume::Dyadic exact(a);exact.add(b,-1);const auto rounded=FiniteVolume::round(exact);
    if(rounded.error.compare(FiniteVolume::Dyadic{})!=0||!std::isfinite(rounded.value))return std::nullopt;
    return rounded.value;
}
std::optional<FragmentPlane> fragmentPlane(const GeometricVofPlicCell& cell,unsigned phase)
{
    if(phase>1||!std::isfinite(cell.measured[0])||!std::isfinite(cell.measured[1])||
        cell.measured[0]<0||cell.measured[1]<0)return std::nullopt;
    FragmentPlane result;
    // Declared empty/full geometry is authoritative even when the retained
    // raw stocks are mixed. Never turn its zero-area plane into a trace.
    if(cell.measured[phase]==0)return result;
    if(cell.measured[1-phase]==0){result.kind=FragmentPlane::Kind::Full;return result;}
    if(!cell.resolved||!std::isfinite(cell.normal.x)||!std::isfinite(cell.normal.y)||
        !std::isfinite(cell.offset))return std::nullopt;
    const auto dx=fragmentExactDifference(cell.rectangle.upper.x,cell.rectangle.lower.x);
    const auto dy=fragmentExactDifference(cell.rectangle.upper.y,cell.rectangle.lower.y);
    if(!dx||!dy||!(*dx>0&&*dy>0))return std::nullopt;
    result.hx=*dx/2;result.hy=*dy/2;
    if(!(result.hx>0&&result.hy>0))return std::nullopt;
    // Align the exact real endpoints with the SAME finite cell-centred frame
    // used by Core PlanarCellCut::shape; unknown translation tails deny proof.
    result.cx=static_cast<long double>(cell.rectangle.lower.x)+result.hx;
    result.cy=static_cast<long double>(cell.rectangle.lower.y)+result.hy;
    for(int axis=0;axis<2;++axis){const auto center=axis?result.cy:result.cx,half=axis?result.hy:result.hx;
        const auto lo=axis?cell.rectangle.lower.y:cell.rectangle.lower.x,hi=axis?cell.rectangle.upper.y:cell.rectangle.upper.x;
        FiniteVolume::Dyadic low(center);low.add(half,-1);FiniteVolume::Dyadic high(center);high.add(half);
        if(low.compare(FiniteVolume::Dyadic(lo))!=0||high.compare(FiniteVolume::Dyadic(hi))!=0)return std::nullopt;}
    const auto length=std::hypot(static_cast<long double>(cell.normal.x),static_cast<long double>(cell.normal.y));
    if(!(std::isfinite(length)&&length>0))return std::nullopt;
    const auto sign=phase==0?1.L:-1.L;
    result.kind=FragmentPlane::Kind::Mixed;result.nx=sign*cell.normal.x/length;
    result.ny=sign*cell.normal.y/length;result.offset=sign*cell.offset;
    if(!std::isfinite(result.nx)||!std::isfinite(result.ny))return std::nullopt;
    return result;
}
std::optional<FiniteVolume::Dyadic> fragmentLevel(const FragmentPlane& plane,Point point)
{
    const auto x=fragmentExactDifference(point.x,plane.cx),y=fragmentExactDifference(point.y,plane.cy);
    if(!x||!y)return std::nullopt;
    auto exact=FiniteVolume::Dyadic::product({plane.nx,*x});
    exact.add(FiniteVolume::Dyadic::product({plane.ny,*y}));exact.add(plane.offset,-1);
    // Native signs and the exact sign must agree. A cancellation rounding to
    // a different side, or onto a node, is not a topological certificate.
    const auto native=plane.nx*(*x)+plane.ny*(*y)-plane.offset;
    if(!std::isfinite(native))return std::nullopt;
    const int sign=exact.compare(FiniteVolume::Dyadic{}),nativeSign=native<0?-1:native>0?1:0;
    if(sign!=nativeSign)return std::nullopt;
    return exact;
}
std::optional<std::array<long double,2>> fragmentPositiveBounds(const FiniteVolume::Dyadic& value)
{
    if(value.compare(FiniteVolume::Dyadic{})<0)return std::nullopt;
    const auto represented=FiniteVolume::round(value);auto lower=represented.value,upper=represented.value;
    const int error=represented.error.compare(FiniteVolume::Dyadic{});
    if(error>0)lower=std::nextafter(lower,-std::numeric_limits<long double>::infinity());
    if(error<0)upper=std::nextafter(upper,std::numeric_limits<long double>::infinity());
    if(!(std::isfinite(lower)&&std::isfinite(upper)&&lower>=0&&upper>=lower))return std::nullopt;
    return std::array<long double,2>{{lower,upper}};
}
std::optional<std::array<long double,2>> fragmentCrossingBounds(
    const FiniteVolume::Dyadic& a,const FiniteVolume::Dyadic& b)
{
    const auto aa=a.absolute(),bb=b.absolute();auto sum=aa;sum.add(bb);
    const auto numerator=fragmentPositiveBounds(aa),denominator=fragmentPositiveBounds(sum);
    if(!numerator||!denominator||!((*denominator)[0]>0))return std::nullopt;
    const auto down=-std::numeric_limits<long double>::infinity(),up=-down;
    auto lower=std::nextafter((*numerator)[0]/(*denominator)[1],down);
    auto upper=std::nextafter((*numerator)[1]/(*denominator)[0],up);
    // Intersect a rigorous parameter enclosure with its exact [0,1] domain;
    // this changes no coordinate, phase shape or material inventory.
    lower=std::max(0.L,lower);upper=std::min(1.L,upper);
    if(!(std::isfinite(lower)&&std::isfinite(upper)&&upper>=lower))return std::nullopt;
    return std::array<long double,2>{{lower,upper}};
}
std::optional<std::array<long double,2>> fragmentCoordinateBounds(
    long double a,long double b,const std::array<long double,2>& parameter)
{
    if(a==b)return std::array<long double,2>{{a,b}};
    const auto span=fragmentExactDifference(b,a);if(!span)return std::nullopt;
    const auto down=-std::numeric_limits<long double>::infinity(),up=-down;
    const auto loT=*span>0?parameter[0]:parameter[1],hiT=*span>0?parameter[1]:parameter[0];
    const auto lower=std::max(std::min(a,b),std::nextafter(a+std::nextafter((*span)*loT,down),down));
    const auto upper=std::min(std::max(a,b),std::nextafter(a+std::nextafter((*span)*hiT,up),up));
    if(!(std::isfinite(lower)&&std::isfinite(upper)&&upper>=lower))return std::nullopt;
    return std::array<long double,2>{{lower,upper}};
}
std::optional<FragmentInterval> fragmentTrace(const FragmentPlane& plane,Point a,Point b,bool horizontal)
{
    const auto lower=horizontal?a.x:a.y,upper=horizontal?b.x:b.y;
    if(!(upper>lower))return std::nullopt;
    if(plane.kind==FragmentPlane::Kind::Empty)return FragmentInterval{};
    if(plane.kind==FragmentPlane::Kind::Full)return FragmentInterval{lower,upper,false};
    const auto va=fragmentLevel(plane,a),vb=fragmentLevel(plane,b);if(!va||!vb)return std::nullopt;
    const auto zero=FiniteVolume::Dyadic{};const int sa=va->compare(zero),sb=vb->compare(zero);
    // A zero endpoint can be a zero-width contact, a corner junction or an
    // exactly coincident branch. This narrow proof never repairs/welds it.
    if(sa==0||sb==0)return std::nullopt;
    if(sa>0&&sb>0)return FragmentInterval{};
    if(sa<0&&sb<0)return FragmentInterval{lower,upper,false};
    const auto parameter=fragmentCrossingBounds(*va,*vb);if(!parameter)return std::nullopt;
    const auto crossing=fragmentCoordinateBounds(lower,upper,*parameter);if(!crossing)return std::nullopt;
    return sa<0?FragmentInterval{lower,(*crossing)[1],false}:FragmentInterval{(*crossing)[0],upper,false};
}
void fragmentInclude(FragmentBounds& bounds,long double x0,long double x1,long double y0,long double y1)
{
    if(!bounds.set){bounds={x0,x1,y0,y1,true};return;}
    bounds.x0=std::min(bounds.x0,x0);bounds.x1=std::max(bounds.x1,x1);
    bounds.y0=std::min(bounds.y0,y0);bounds.y1=std::max(bounds.y1,y1);
}
std::optional<long double> fragmentDiameter(const GeometricVofPlicCell& cell,
    const FragmentPlane& plane,long double delta)
{
    if(plane.kind!=FragmentPlane::Kind::Mixed)return std::nullopt;
    const auto original=Cut::interfaceSegment(cell.rectangle,cell.normal,cell.offset);
    if(!original)return std::nullopt;
    // These local endpoints are the actual Core finite PLIC segment. Include
    // them as well as outward enclosures of the corresponding exact plane.
    FragmentBounds bounds;
    for(const auto point:*original){if(!std::isfinite(point.x)||!std::isfinite(point.y))return std::nullopt;
        fragmentInclude(bounds,point.x,point.x,point.y,point.y);}
    const std::array<Point,4> corner{{{-plane.hx,-plane.hy},{plane.hx,-plane.hy},
        {plane.hx,plane.hy},{-plane.hx,plane.hy}}};
    std::array<FiniteVolume::Dyadic,4> level;std::array<int,4> signs{};
    int inside=0,crossings=0;
    for(unsigned i=0;i<4;++i){const auto q=corner[i];
        level[i]=FiniteVolume::Dyadic::product({plane.nx,q.x});level[i].add(FiniteVolume::Dyadic::product({plane.ny,q.y}));level[i].add(plane.offset,-1);
        signs[i]=level[i].compare(FiniteVolume::Dyadic{});const auto native=plane.nx*q.x+plane.ny*q.y-plane.offset;
        const int nativeSign=native<0?-1:native>0?1:0;
        if(signs[i]==0||signs[i]!=nativeSign||!std::isfinite(native))return std::nullopt;
        if(signs[i]<0){++inside;fragmentInclude(bounds,q.x,q.x,q.y,q.y);}
    }
    if(inside==0||inside==4)return std::nullopt;
    for(unsigned i=0;i<4;++i){const unsigned next=(i+1)%4;if(signs[i]==signs[next])continue;
        ++crossings;const auto parameter=fragmentCrossingBounds(level[i],level[next]);if(!parameter)return std::nullopt;
        const auto x=fragmentCoordinateBounds(corner[i].x,corner[next].x,*parameter);
        const auto y=fragmentCoordinateBounds(corner[i].y,corner[next].y,*parameter);if(!x||!y)return std::nullopt;
        fragmentInclude(bounds,(*x)[0],(*x)[1],(*y)[0],(*y)[1]);
    }
    if(crossings!=2||!bounds.set)return std::nullopt;
    const auto outwardWidth=[](long double a,long double b){
        FiniteVolume::Dyadic exact(b);exact.add(a,-1);const auto rounded=FiniteVolume::round(exact);
        return rounded.error.compare(FiniteVolume::Dyadic{})<0?
            std::nextafter(rounded.value,std::numeric_limits<long double>::infinity()):rounded.value;};
    const auto dx=outwardWidth(bounds.x0,bounds.x1),dy=outwardWidth(bounds.y0,bounds.y1);
    if(!(std::isfinite(dx)&&std::isfinite(dy)&&dx>0&&dy>0))return std::nullopt;
    auto square=FiniteVolume::Dyadic::product({dx,dx});square.add(FiniteVolume::Dyadic::product({dy,dy}));
    if(square.compare(FiniteVolume::Dyadic::product({delta,delta}))>=0)return std::nullopt;
    auto diameter=std::hypot(dx,dy);if(!(std::isfinite(diameter)&&diameter>0))return std::nullopt;
    // Certify the actual returned bound rather than assuming hypot rounds up.
    for(unsigned i=0;i<4&&FiniteVolume::Dyadic::product({diameter,diameter}).compare(square)<0;++i)
        diameter=std::nextafter(diameter,std::numeric_limits<long double>::infinity());
    if(!(diameter<delta)||FiniteVolume::Dyadic::product({diameter,diameter}).compare(square)<0)return std::nullopt;
    return diameter;
}
std::optional<FragmentProof> fragmentProof(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,
    const std::vector<GeometricVofPlicCell>& local,const std::vector<GeometricVofFace>& faces,
    int native,unsigned phase)
{
    const auto& source=local[native];const auto plane=fragmentPlane(source,phase);
    if(!plane||plane->kind!=FragmentPlane::Kind::Mixed)return std::nullopt;
    const auto dx=fragmentExactDifference(source.rectangle.upper.x,source.rectangle.lower.x);
    const auto dy=fragmentExactDifference(source.rectangle.upper.y,source.rectangle.lower.y);
    if(!dx||!dy)return std::nullopt;const auto delta=std::min(*dx,*dy);
    if(!(std::isfinite(delta)&&delta>0))return std::nullopt;
    const auto diameter=fragmentDiameter(source,*plane,delta);if(!diameter)return std::nullopt;
    std::array<std::vector<std::pair<long double,long double>>,4> coverage;
    // LocalPlic cells do not publish boundaryVertices. Recover only these
    // actual source endpoint IDs from the retained native source-face rows.
    std::set<std::size_t> sourceNodes;
    for(std::size_t row=0;row<faces.size();++row){
        if(faces[row].nativeRow!=static_cast<int>(row))return std::nullopt;
        const auto face=domain.face(mesh,static_cast<int>(row));
        if(face.negativeCell==native||face.positiveCell==native)
            sourceNodes.insert(faces[row].vertices.begin(),faces[row].vertices.end());
    }
    if(sourceNodes.size()<4||sourceNodes.size()>8)return std::nullopt;
    std::map<std::size_t,std::set<int>> nodeCells;std::set<int> checkedCorners;
    std::size_t faceCount=0;
    for(std::size_t row=0;row<faces.size();++row){const auto& retained=faces[row];
        if(retained.nativeRow!=static_cast<int>(row))return std::nullopt;
        const auto f=domain.face(mesh,static_cast<int>(row));
        for(const auto v:retained.vertices)if(sourceNodes.count(v)){
            if(f.negativeCell>=0)nodeCells[v].insert(f.negativeCell);
            if(f.positiveCell>=0)nodeCells[v].insert(f.positiveCell);}
        if(f.negativeCell!=native&&f.positiveCell!=native)continue;
        ++faceCount;if(f.negativeCell<0||f.positiveCell<0||!(f.area>0))return std::nullopt;
        const auto& edge=f.horizontal?mesh.hEdges(f.edgeLabel):mesh.vEdges(f.edgeLabel);
        const auto av=mesh.verts(edge.edgeToVert[0]).coord(),bv=mesh.verts(edge.edgeToVert[1]).coord();
        Point a{av.x,av.y},b{bv.x,bv.y};
        if(f.horizontal?a.x>b.x:a.y>b.y)std::swap(a,b);
        const auto constant=f.horizontal?a.y:a.x;
        if((f.horizontal?b.y:b.x)!=constant)return std::nullopt;
        const int side=f.horizontal?(constant==source.rectangle.lower.y?2:constant==source.rectangle.upper.y?3:-1):
            (constant==source.rectangle.lower.x?0:constant==source.rectangle.upper.x?1:-1);
        if(side<0||((side==0||side==2)?f.positiveCell!=native:f.negativeCell!=native))return std::nullopt;
        const auto lo=f.horizontal?a.x:a.y,hi=f.horizontal?b.x:b.y;
        const auto sideLo=f.horizontal?source.rectangle.lower.x:source.rectangle.lower.y;
        const auto sideHi=f.horizontal?source.rectangle.upper.x:source.rectangle.upper.y;
        if(!(hi>lo&&lo>=sideLo&&hi<=sideHi))return std::nullopt;
        coverage[side].push_back({lo,hi});
        const int neighbor=f.negativeCell==native?f.positiveCell:f.negativeCell;
        if(neighbor<0||static_cast<std::size_t>(neighbor)>=local.size())return std::nullopt;
        const auto neighborPlane=fragmentPlane(local[neighbor],phase);if(!neighborPlane)return std::nullopt;
        const auto own=fragmentTrace(*plane,a,b,f.horizontal),other=fragmentTrace(*neighborPlane,a,b,f.horizontal);
        if(!own||!other)return std::nullopt;
        // Closed conservative enclosures must be STRICTLY disjoint. Equality
        // or overlap might be a same-phase connection or unresolved contact.
        if(!own->empty&&!other->empty&&
            std::max(own->lower,other->lower)<=std::min(own->upper,other->upper))return std::nullopt;
    }
    if(faceCount<4||faceCount>8)return std::nullopt;
    for(int side=0;side<4;++side){auto& intervals=coverage[side];if(intervals.empty())return std::nullopt;
        std::sort(intervals.begin(),intervals.end());
        const auto lo=side<2?source.rectangle.lower.y:source.rectangle.lower.x;
        const auto hi=side<2?source.rectangle.upper.y:source.rectangle.upper.x;
        auto cursor=static_cast<long double>(lo);
        for(const auto& interval:intervals){if(interval.first!=cursor)return std::nullopt;cursor=interval.second;}
        if(cursor!=hi)return std::nullopt;
    }
    // Include every actual endpoint, including hanging nodes, and ALL cells
    // incident on it, not only the two cells sharing a source face. Four
    // quadrant coverage rejects a missing diagonal cell or nonmanifold patch.
    for(const auto v:sourceNodes){
        if(v>=static_cast<std::size_t>(mesh.vertIndexNum()))return std::nullopt;
        const auto coord=mesh.verts(mesh.vertIndexToLabel[static_cast<int>(v)]).coord();const Point point{coord.x,coord.y};
        const auto ownLevel=fragmentLevel(*plane,point);if(!ownLevel)return std::nullopt;
        const int ownSign=ownLevel->compare(FiniteVolume::Dyadic{});if(ownSign==0)return std::nullopt;
        const auto found=nodeCells.find(v);if(found==nodeCells.end()||!found->second.count(native))return std::nullopt;
        std::array<unsigned,4> quadrants{};
        for(const int ci:found->second){if(ci<0||static_cast<std::size_t>(ci)>=local.size())return std::nullopt;
            const auto& cell=local[ci];const auto& box=cell.rectangle;
            if(point.x<box.lower.x||point.x>box.upper.x||point.y<box.lower.y||point.y>box.upper.y)return std::nullopt;
            const std::array<bool,2> xs{{box.lower.x<point.x,box.upper.x>point.x}},ys{{box.lower.y<point.y,box.upper.y>point.y}};
            for(unsigned x=0;x<2;++x)for(unsigned y=0;y<2;++y)if(xs[x]&&ys[y])++quadrants[2*x+y];
            if(ci==native)continue;checkedCorners.insert(ci);
            const auto other=fragmentPlane(cell,phase);if(!other)return std::nullopt;
            if(ownSign<0){
                if(other->kind==FragmentPlane::Kind::Full)return std::nullopt;
                if(other->kind==FragmentPlane::Kind::Mixed){const auto value=fragmentLevel(*other,point);
                    if(!value||value->compare(FiniteVolume::Dyadic{})<=0)return std::nullopt;}
            }
        }
        for(const auto count:quadrants)if(count!=1)return std::nullopt;
    }
    return FragmentProof{phase,*diameter,delta,faceCount,checkedCorners.size()};
}
// D23: a distinct two-bin position-deficit proof; singleton diameter law above
// is deliberately unchanged. Everything here is read-only actual geometry.
struct TwoCellProof
{
    std::array<int,2> native{};unsigned phase=0,splitAxis=0;
    FragmentBounds bounds;long double delta=0,bound=0;
    std::size_t faces=0,corners=0;
};
std::optional<FragmentBounds> twoCellLocalBounds(const GeometricVofPlicCell& cell,const FragmentPlane& plane)
{
    if(plane.kind!=FragmentPlane::Kind::Mixed)return std::nullopt;
    const auto original=Cut::interfaceSegment(cell.rectangle,cell.normal,cell.offset);
    if(!original)return std::nullopt;
    FragmentBounds bounds;
    for(const auto point:*original){if(!std::isfinite(point.x)||!std::isfinite(point.y))return std::nullopt;
        fragmentInclude(bounds,point.x,point.x,point.y,point.y);}
    const std::array<Point,4> corner{{{-plane.hx,-plane.hy},{plane.hx,-plane.hy},
        {plane.hx,plane.hy},{-plane.hx,plane.hy}}};
    std::array<FiniteVolume::Dyadic,4> level;std::array<int,4> signs{};int inside=0,crossings=0;
    for(unsigned i=0;i<4;++i){const auto point=corner[i];
        level[i]=FiniteVolume::Dyadic::product({plane.nx,point.x});
        level[i].add(FiniteVolume::Dyadic::product({plane.ny,point.y}));level[i].add(plane.offset,-1);
        signs[i]=level[i].compare(FiniteVolume::Dyadic{});
        const auto native=plane.nx*point.x+plane.ny*point.y-plane.offset;
        const int nativeSign=native<0?-1:native>0?1:0;
        if(signs[i]==0||signs[i]!=nativeSign||!std::isfinite(native))return std::nullopt;
        if(signs[i]<0){++inside;fragmentInclude(bounds,point.x,point.x,point.y,point.y);}}
    if(inside==0||inside==4)return std::nullopt;
    for(unsigned i=0;i<4;++i){const unsigned next=(i+1)%4;if(signs[i]==signs[next])continue;
        ++crossings;const auto parameter=fragmentCrossingBounds(level[i],level[next]);if(!parameter)return std::nullopt;
        const auto x=fragmentCoordinateBounds(corner[i].x,corner[next].x,*parameter);
        const auto y=fragmentCoordinateBounds(corner[i].y,corner[next].y,*parameter);if(!x||!y)return std::nullopt;
        fragmentInclude(bounds,(*x)[0],(*x)[1],(*y)[0],(*y)[1]);}
    if(crossings!=2||!bounds.set)return std::nullopt;
    return bounds;
}
std::optional<long double> twoCellTranslatedBound(long double value,long double shift,bool upper)
{
    FiniteVolume::Dyadic exact(value);exact.add(shift);const auto represented=FiniteVolume::round(exact);
    auto result=represented.value;const int error=represented.error.compare(FiniteVolume::Dyadic{});
    if((upper&&error<0)||(!upper&&error>0))result=std::nextafter(result,
        upper?std::numeric_limits<long double>::infinity():-std::numeric_limits<long double>::infinity());
    if(!std::isfinite(result))return std::nullopt;return result;
}
std::optional<TwoCellProof> twoCellProof(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,
    const std::vector<GeometricVofPlicCell>& local,const std::vector<GeometricVofFace>& faces,
    std::array<int,2> pair,unsigned phase)
{
    if(pair[0]>=pair[1]||pair[0]<0||static_cast<std::size_t>(pair[1])>=local.size())return std::nullopt;
    const std::array<GeometricVofPlicCell,2> source{{local[pair[0]],local[pair[1]]}};
    std::array<FragmentPlane,2> plane;
    for(unsigned i=0;i<2;++i){if(!source[i].resolved||!(source[i].measured[0]>0&&source[i].measured[1]>0))return std::nullopt;
        const auto actual=fragmentPlane(source[i],phase);if(!actual||actual->kind!=FragmentPlane::Kind::Mixed)return std::nullopt;plane[i]=*actual;}
    if(mesh.cells(mesh.cellIndexToLabel[pair[0]]).level()!=mesh.cells(mesh.cellIndexToLabel[pair[1]]).level()||
        plane[0].hx!=plane[1].hx||plane[0].hy!=plane[1].hy)return std::nullopt;
    const auto& a=source[0].rectangle;const auto& b=source[1].rectangle;
    const bool horizontal=a.lower.x==b.lower.x&&a.upper.x==b.upper.x&&(a.lower.y==b.upper.y||a.upper.y==b.lower.y);
    const bool vertical=a.lower.y==b.lower.y&&a.upper.y==b.upper.y&&(a.lower.x==b.upper.x||a.upper.x==b.lower.x);
    if(horizontal==vertical)return std::nullopt;
    const auto dx=fragmentExactDifference(a.upper.x,a.lower.x),dy=fragmentExactDifference(a.upper.y,a.lower.y);
    if(!dx||!dy||!(*dx>0&&*dy>0))return std::nullopt;const auto delta=std::min(*dx,*dy);
    FragmentBounds bounds;
    for(unsigned i=0;i<2;++i){const auto part=twoCellLocalBounds(source[i],plane[i]);if(!part)return std::nullopt;
        const auto sx=fragmentExactDifference(plane[i].cx,plane[0].cx),sy=fragmentExactDifference(plane[i].cy,plane[0].cy);
        if(!sx||!sy)return std::nullopt;
        const auto x0=twoCellTranslatedBound(part->x0,*sx,false),x1=twoCellTranslatedBound(part->x1,*sx,true);
        const auto y0=twoCellTranslatedBound(part->y0,*sy,false),y1=twoCellTranslatedBound(part->y1,*sy,true);
        if(!x0||!x1||!y0||!y1)return std::nullopt;fragmentInclude(bounds,*x0,*x1,*y0,*y1);}
    if(!(bounds.set&&bounds.x1>bounds.x0&&bounds.y1>bounds.y0))return std::nullopt;
    // No rounded width/halfwidth is used in the strict qualification. These
    // are exact dyadic expressions in already outward-enclosing endpoints.
    auto splitX=FiniteVolume::Dyadic::product({.5L,bounds.x1});splitX.add(FiniteVolume::Dyadic::product({.5L,bounds.x0}),-1);
    splitX.add(bounds.y1);splitX.add(bounds.y0,-1);
    auto splitY=FiniteVolume::Dyadic::product({.5L,bounds.y1});splitY.add(FiniteVolume::Dyadic::product({.5L,bounds.y0}),-1);
    splitY.add(bounds.x1);splitY.add(bounds.x0,-1);
    const unsigned splitAxis=splitX.compare(splitY)<=0?0u:1u;const auto& chosen=splitAxis==0?splitX:splitY;
    if(chosen.compare(FiniteVolume::Dyadic(delta))>=0)return std::nullopt;
    const auto represented=FiniteVolume::round(chosen);auto bound=represented.value;
    if(represented.error.compare(FiniteVolume::Dyadic{})<0)bound=std::nextafter(bound,std::numeric_limits<long double>::infinity());
    if(!(std::isfinite(bound)&&bound>0&&bound<delta))return std::nullopt;
    // Bins are geometric halves of the AABB, not two selected sample points.
    // Any three component positions put two in one bin; their Euclidean
    // distance is bounded by this strict L1 bin bound, hence is < Delta.
    std::array<std::array<std::vector<std::pair<long double,long double>>,4>,2> coverage;
    std::set<std::size_t> sourceNodes;
    for(std::size_t row=0;row<faces.size();++row){if(faces[row].nativeRow!=static_cast<int>(row))return std::nullopt;
        const auto f=domain.face(mesh,static_cast<int>(row));
        if(f.negativeCell==pair[0]||f.positiveCell==pair[0]||f.negativeCell==pair[1]||f.positiveCell==pair[1])
            sourceNodes.insert(faces[row].vertices.begin(),faces[row].vertices.end());}
    if(sourceNodes.size()<6||sourceNodes.size()>16)return std::nullopt;
    std::map<std::size_t,std::set<int>> nodeCells;std::set<int> checkedCorners;
    std::size_t faceCount=0,sharedCount=0;std::array<std::size_t,2> ownCounts{};
    for(std::size_t row=0;row<faces.size();++row){const auto& retained=faces[row];const auto f=domain.face(mesh,static_cast<int>(row));
        for(const auto v:retained.vertices)if(sourceNodes.count(v)){
            if(f.negativeCell>=0)nodeCells[v].insert(f.negativeCell);if(f.positiveCell>=0)nodeCells[v].insert(f.positiveCell);}
        const bool touches0=f.negativeCell==pair[0]||f.positiveCell==pair[0],touches1=f.negativeCell==pair[1]||f.positiveCell==pair[1];
        if(!touches0&&!touches1)continue;++faceCount;
        if(f.negativeCell<0||f.positiveCell<0||!(f.area>0))return std::nullopt;
        const auto& edge=f.horizontal?mesh.hEdges(f.edgeLabel):mesh.vEdges(f.edgeLabel);
        const auto av=mesh.verts(edge.edgeToVert[0]).coord(),bv=mesh.verts(edge.edgeToVert[1]).coord();Point p{av.x,av.y},q{bv.x,bv.y};
        if(f.horizontal?p.x>q.x:p.y>q.y)std::swap(p,q);
        const auto constant=f.horizontal?p.y:p.x;if((f.horizontal?q.y:q.x)!=constant)return std::nullopt;
        const auto lo=f.horizontal?p.x:p.y,hi=f.horizontal?q.x:q.y;if(!(hi>lo))return std::nullopt;
        for(unsigned i=0;i<2;++i){if(i==0?!touches0:!touches1)continue;++ownCounts[i];const auto& box=source[i].rectangle;
            const int side=f.horizontal?(constant==box.lower.y?2:constant==box.upper.y?3:-1):(constant==box.lower.x?0:constant==box.upper.x?1:-1);
            if(side<0||((side==0||side==2)?f.positiveCell!=pair[i]:f.negativeCell!=pair[i]))return std::nullopt;
            const auto sideLo=f.horizontal?box.lower.x:box.lower.y,sideHi=f.horizontal?box.upper.x:box.upper.y;
            if(!(lo>=sideLo&&hi<=sideHi))return std::nullopt;coverage[i][side].push_back({lo,hi});}
        if(touches0&&touches1){++sharedCount;
            const auto first=fragmentTrace(plane[0],p,q,f.horizontal),second=fragmentTrace(plane[1],p,q,f.horizontal);
            if(!first||!second||first->empty||second->empty||first->lower!=lo||first->upper!=hi||second->lower!=lo||second->upper!=hi)return std::nullopt;
            // The first scope requires one complete, exactly matching face.
            const auto& box=source[0].rectangle;
            if(lo!=(f.horizontal?box.lower.x:box.lower.y)||hi!=(f.horizontal?box.upper.x:box.upper.y))return std::nullopt;
        }else{const unsigned i=touches0?0u:1u;const int neighbor=f.negativeCell==pair[i]?f.positiveCell:f.negativeCell;
            if(neighbor<0||static_cast<std::size_t>(neighbor)>=local.size())return std::nullopt;
            const auto neighborPlane=fragmentPlane(local[neighbor],phase);if(!neighborPlane)return std::nullopt;
            const auto own=fragmentTrace(plane[i],p,q,f.horizontal),other=fragmentTrace(*neighborPlane,p,q,f.horizontal);
            if(!own||!other)return std::nullopt;
            if(!own->empty&&!other->empty&&std::max(own->lower,other->lower)<=std::min(own->upper,other->upper))return std::nullopt;}
    }
    if(sharedCount!=1)return std::nullopt;
    for(unsigned i=0;i<2;++i){if(ownCounts[i]<4||ownCounts[i]>8)return std::nullopt;
        for(int side=0;side<4;++side){auto& intervals=coverage[i][side];if(intervals.empty())return std::nullopt;std::sort(intervals.begin(),intervals.end());
            const auto& box=source[i].rectangle;const auto lo=side<2?box.lower.y:box.lower.x,hi=side<2?box.upper.y:box.upper.x;
            auto cursor=static_cast<long double>(lo);for(const auto& interval:intervals){if(interval.first!=cursor)return std::nullopt;cursor=interval.second;}
            if(cursor!=hi)return std::nullopt;}}
    for(const auto v:sourceNodes){if(v>=static_cast<std::size_t>(mesh.vertIndexNum()))return std::nullopt;
        const auto coord=mesh.verts(mesh.vertIndexToLabel[static_cast<int>(v)]).coord();const Point point{coord.x,coord.y};
        const auto found=nodeCells.find(v);if(found==nodeCells.end())return std::nullopt;
        std::array<unsigned,4> quadrants{};bool occupied=false;
        for(unsigned i=0;i<2;++i){const auto& box=source[i].rectangle;
            if(point.x<box.lower.x||point.x>box.upper.x||point.y<box.lower.y||point.y>box.upper.y)continue;
            if(!found->second.count(pair[i]))return std::nullopt;
            const auto level=fragmentLevel(plane[i],point);if(!level||level->compare(FiniteVolume::Dyadic{})==0)return std::nullopt;
            occupied=occupied||level->compare(FiniteVolume::Dyadic{})<0;}
        for(const int ci:found->second){if(ci<0||static_cast<std::size_t>(ci)>=local.size())return std::nullopt;
            const auto& cell=local[ci];const auto& box=cell.rectangle;
            if(point.x<box.lower.x||point.x>box.upper.x||point.y<box.lower.y||point.y>box.upper.y)return std::nullopt;
            const std::array<bool,2> xs{{box.lower.x<point.x,box.upper.x>point.x}},ys{{box.lower.y<point.y,box.upper.y>point.y}};
            for(unsigned x=0;x<2;++x)for(unsigned y=0;y<2;++y)if(xs[x]&&ys[y])++quadrants[2*x+y];
            if(ci==pair[0]||ci==pair[1])continue;checkedCorners.insert(ci);
            const auto other=fragmentPlane(cell,phase);if(!other)return std::nullopt;
            if(occupied){if(other->kind==FragmentPlane::Kind::Full)return std::nullopt;
                if(other->kind==FragmentPlane::Kind::Mixed){const auto level=fragmentLevel(*other,point);
                    if(!level||level->compare(FiniteVolume::Dyadic{})<=0)return std::nullopt;}}
        }
        for(const auto count:quadrants)if(count!=1)return std::nullopt;
    }
    return TwoCellProof{pair,phase,splitAxis,bounds,delta,bound,faceCount,checkedCorners.size()};
}
// F25: bounded whole published-geometry component, not a stock topology.
struct VertexStarTrace
{
    FragmentInterval outer,inner; // Outer excludes contacts; inner proves positive overlap.
};
std::optional<VertexStarTrace> vertexStarTrace(const FragmentPlane& plane,Point a,Point b,bool horizontal)
{
    const auto lo=horizontal?a.x:a.y,hi=horizontal?b.x:b.y;if(!(hi>lo))return std::nullopt;
    if(plane.kind==FragmentPlane::Kind::Empty)return VertexStarTrace{};
    if(plane.kind==FragmentPlane::Kind::Full){const FragmentInterval full{lo,hi,false};return VertexStarTrace{full,full};}
    const auto va=fragmentLevel(plane,a),vb=fragmentLevel(plane,b);if(!va||!vb)return std::nullopt;
    const auto zero=FiniteVolume::Dyadic{};const int sa=va->compare(zero),sb=vb->compare(zero);
    if(sa==0||sb==0)return std::nullopt;
    if(sa>0&&sb>0)return VertexStarTrace{};
    if(sa<0&&sb<0){const FragmentInterval full{lo,hi,false};return VertexStarTrace{full,full};}
    const auto parameter=fragmentCrossingBounds(*va,*vb);if(!parameter)return std::nullopt;
    const auto crossing=fragmentCoordinateBounds(lo,hi,*parameter);if(!crossing)return std::nullopt;
    VertexStarTrace result;
    if(sa<0){result.outer={lo,(*crossing)[1],false};result.inner={lo,(*crossing)[0],!(lo<(*crossing)[0])};}
    else{result.outer={(*crossing)[0],hi,false};result.inner={(*crossing)[1],hi,!((*crossing)[1]<hi)};}
    return result;
}
bool vertexStarContact(const FragmentInterval& a,const FragmentInterval& b)
{return !a.empty&&!b.empty&&std::max(a.lower,b.lower)<=std::min(a.upper,b.upper);}
bool vertexStarPositive(const FragmentInterval& a,const FragmentInterval& b)
{return !a.empty&&!b.empty&&std::max(a.lower,b.lower)<std::min(a.upper,b.upper);}
struct VertexStarProof
{
    std::vector<int> native;Point vertex;unsigned phase=0;
    FragmentBounds bounds;long double delta=0,diameter=0;std::size_t faces=0,corners=0;
};
std::optional<VertexStarProof> vertexStarProof(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,
    const std::vector<GeometricVofPlicCell>& local,const std::vector<GeometricVofFace>& faces,
    int native,std::size_t vertex,unsigned phase)
{
    if(vertex>=static_cast<std::size_t>(mesh.vertIndexNum()))return std::nullopt;
    const auto xy=mesh.verts(mesh.vertIndexToLabel[static_cast<int>(vertex)]).coord();const Point point{xy.x,xy.y};
    std::set<int> starSet;
    for(std::size_t row=0;row<faces.size();++row){if(faces[row].nativeRow!=static_cast<int>(row))return std::nullopt;
        if(std::find(faces[row].vertices.begin(),faces[row].vertices.end(),vertex)==faces[row].vertices.end())continue;
        const auto f=domain.face(mesh,static_cast<int>(row));
        if(f.negativeCell>=0)starSet.insert(f.negativeCell);if(f.positiveCell>=0)starSet.insert(f.positiveCell);}
    if(starSet.size()!=4||!starSet.count(native))return std::nullopt;
    const std::vector<int> star(starSet.begin(),starSet.end());std::array<FragmentPlane,4> planes;
    const auto& original=local[native];const auto dx=fragmentExactDifference(original.rectangle.upper.x,original.rectangle.lower.x);
    const auto dy=fragmentExactDifference(original.rectangle.upper.y,original.rectangle.lower.y);
    if(!dx||!dy||!(*dx>0&&*dy>0))return std::nullopt;const auto delta=std::min(*dx,*dy);
    const int level=mesh.cells(mesh.cellIndexToLabel[native]).level();std::array<unsigned,4> quadrants{};
    const auto index=[&](int cell){const auto it=std::lower_bound(star.begin(),star.end(),cell);return it!=star.end()&&*it==cell?static_cast<int>(it-star.begin()):-1;};
    for(unsigned i=0;i<4;++i){const int ci=star[i];if(ci<0||static_cast<std::size_t>(ci)>=local.size())return std::nullopt;
        const auto& cell=local[ci];const auto& b=cell.rectangle;
        if(mesh.cells(mesh.cellIndexToLabel[ci]).level()!=level||
            fragmentExactDifference(b.upper.x,b.lower.x)!=dx||fragmentExactDifference(b.upper.y,b.lower.y)!=dy)return std::nullopt;
        if(!((point.x==b.lower.x||point.x==b.upper.x)&&(point.y==b.lower.y||point.y==b.upper.y)))return std::nullopt;
        const unsigned q=2u*(b.lower.x==point.x?1u:0u)+(b.lower.y==point.y?1u:0u);++quadrants[q];
        const auto plane=fragmentPlane(cell,phase);if(!plane)return std::nullopt;planes[i]=*plane;
    }
    for(unsigned count:quadrants)if(count!=1)return std::nullopt;
    const int source=index(native);if(source<0||planes[source].kind!=FragmentPlane::Kind::Mixed)return std::nullopt;
    const auto sourceLevel=fragmentLevel(planes[source],point);
    if(!sourceLevel||sourceLevel->compare(FiniteVolume::Dyadic{})>=0)return std::nullopt;
    struct Face {std::size_t row;int negative,positive;Point a,b;bool horizontal;};
    std::vector<Face> actual;std::array<std::array<std::vector<std::pair<long double,long double>>,4>,4> coverage;
    std::array<std::vector<int>,4> graph;std::array<std::size_t,4> ownCounts{};
    for(std::size_t row=0;row<faces.size();++row){const auto f=domain.face(mesh,static_cast<int>(row));
        const int minus=index(f.negativeCell),plus=index(f.positiveCell);if(minus<0&&plus<0)continue;
        if(f.negativeCell<0||f.positiveCell<0||!(f.area>0))return std::nullopt;
        const auto& edge=f.horizontal?mesh.hEdges(f.edgeLabel):mesh.vEdges(f.edgeLabel);
        const auto av=mesh.verts(edge.edgeToVert[0]).coord(),bv=mesh.verts(edge.edgeToVert[1]).coord();Point a{av.x,av.y},b{bv.x,bv.y};
        if(f.horizontal?a.x>b.x:a.y>b.y)std::swap(a,b);
        const auto constant=f.horizontal?a.y:a.x,lo=f.horizontal?a.x:a.y,hi=f.horizontal?b.x:b.y;
        if((f.horizontal?b.y:b.x)!=constant||!(hi>lo))return std::nullopt;
        actual.push_back({row,f.negativeCell,f.positiveCell,a,b,f.horizontal});
        for(const int i:{minus,plus})if(i>=0){++ownCounts[i];const auto& box=local[star[i]].rectangle;
            const int side=f.horizontal?(constant==box.lower.y?2:constant==box.upper.y?3:-1):(constant==box.lower.x?0:constant==box.upper.x?1:-1);
            if(side<0||((side==0||side==2)?f.positiveCell!=star[i]:f.negativeCell!=star[i]))return std::nullopt;
            const auto sideLo=f.horizontal?box.lower.x:box.lower.y,sideHi=f.horizontal?box.upper.x:box.upper.y;
            if(!(lo>=sideLo&&hi<=sideHi))return std::nullopt;coverage[i][side].push_back({lo,hi});}
        if(minus>=0&&plus>=0){const auto m=vertexStarTrace(planes[minus],a,b,f.horizontal),q=vertexStarTrace(planes[plus],a,b,f.horizontal);
            if(!m||!q)return std::nullopt;
            if(vertexStarPositive(m->inner,q->inner)){graph[minus].push_back(plus);graph[plus].push_back(minus);}
            else if(vertexStarContact(m->outer,q->outer))return std::nullopt; // Unknown/point contact is not a positive edge.
        }
    }
    for(unsigned i=0;i<4;++i){if(ownCounts[i]<4||ownCounts[i]>8)return std::nullopt;
        for(int side=0;side<4;++side){auto& intervals=coverage[i][side];if(intervals.empty())return std::nullopt;
            std::sort(intervals.begin(),intervals.end());const auto& box=local[star[i]].rectangle;
            const auto lo=side<2?box.lower.y:box.lower.x,hi=side<2?box.upper.y:box.upper.x;
            auto cursor=static_cast<long double>(lo);for(const auto& part:intervals){if(part.first!=cursor)return std::nullopt;cursor=part.second;}
            if(cursor!=hi)return std::nullopt;}}
    std::set<int> member{source};std::vector<int> queue{source};
    for(std::size_t at=0;at<queue.size();++at)for(int other:graph[queue[at]])if(member.insert(other).second)queue.push_back(other);
    std::vector<int> component;for(int i:member){const auto& c=local[star[i]];
        if(planes[i].kind!=FragmentPlane::Kind::Mixed||!c.resolved||!(c.measured[0]>0&&c.measured[1]>0))return std::nullopt;
        const auto value=fragmentLevel(planes[i],point);if(!value||value->compare(FiniteVolume::Dyadic{})>=0)return std::nullopt;
        component.push_back(star[i]);}
    std::sort(component.begin(),component.end());
    const auto inComponent=[&](int ci){return std::binary_search(component.begin(),component.end(),ci);};
    std::set<std::size_t> nodes;std::size_t faceCount=0;
    for(const auto& f:actual){const bool m=inComponent(f.negative),q=inComponent(f.positive);if(!m&&!q)continue;
        ++faceCount;nodes.insert(faces[f.row].vertices.begin(),faces[f.row].vertices.end());
        if(m&&q)continue; // Unequal internal traces add union-boundary facets; never weld/normalize them.
        const int own=m?f.negative:f.positive,neighbor=m?f.positive:f.negative;
        if(neighbor<0||static_cast<std::size_t>(neighbor)>=local.size())return std::nullopt;
        const auto otherPlane=fragmentPlane(local[neighbor],phase);if(!otherPlane)return std::nullopt;
        const auto ownTrace=vertexStarTrace(planes[index(own)],f.a,f.b,f.horizontal),otherTrace=vertexStarTrace(*otherPlane,f.a,f.b,f.horizontal);
        if(!ownTrace||!otherTrace||vertexStarContact(ownTrace->outer,otherTrace->outer))return std::nullopt;
    }
    if(nodes.empty())return std::nullopt;
    std::map<std::size_t,std::set<int>> nodeCells;
    for(std::size_t row=0;row<faces.size();++row){const auto f=domain.face(mesh,static_cast<int>(row));
        for(const auto v:faces[row].vertices)if(nodes.count(v)){
            if(f.negativeCell>=0)nodeCells[v].insert(f.negativeCell);if(f.positiveCell>=0)nodeCells[v].insert(f.positiveCell);}}
    std::set<int> checkedCorners;
    for(const auto v:nodes){if(v>=static_cast<std::size_t>(mesh.vertIndexNum()))return std::nullopt;
        const auto coord=mesh.verts(mesh.vertIndexToLabel[static_cast<int>(v)]).coord();const Point node{coord.x,coord.y};
        const auto found=nodeCells.find(v);if(found==nodeCells.end())return std::nullopt;
        bool occupied=false;
        for(const int ci:component){const auto& box=local[ci].rectangle;
            if(node.x<box.lower.x||node.x>box.upper.x||node.y<box.lower.y||node.y>box.upper.y)continue;
            if(!found->second.count(ci))return std::nullopt;
            const auto value=fragmentLevel(planes[index(ci)],node);if(!value||value->compare(FiniteVolume::Dyadic{})==0)return std::nullopt;
            occupied=occupied||value->compare(FiniteVolume::Dyadic{})<0;}
        std::array<unsigned,4> cover{};
        for(int ci:found->second){if(ci<0||static_cast<std::size_t>(ci)>=local.size())return std::nullopt;
            const auto& box=local[ci].rectangle;
            if(node.x<box.lower.x||node.x>box.upper.x||node.y<box.lower.y||node.y>box.upper.y)return std::nullopt;
            const std::array<bool,2> xs{{box.lower.x<node.x,box.upper.x>node.x}},ys{{box.lower.y<node.y,box.upper.y>node.y}};
            for(unsigned x=0;x<2;++x)for(unsigned y=0;y<2;++y)if(xs[x]&&ys[y])++cover[2*x+y];
            if(inComponent(ci))continue;checkedCorners.insert(ci);
            const auto other=fragmentPlane(local[ci],phase);if(!other)return std::nullopt;
            if(occupied){if(other->kind==FragmentPlane::Kind::Full)return std::nullopt;
                if(other->kind==FragmentPlane::Kind::Mixed){const auto value=fragmentLevel(*other,node);
                    if(!value||value->compare(FiniteVolume::Dyadic{})<=0)return std::nullopt;}}
        }
        for(unsigned n:cover)if(n!=1)return std::nullopt;
    }
    FragmentBounds bounds;const auto& anchor=planes[index(component.front())];
    for(int ci:component){const auto& plane=planes[index(ci)];const auto b=twoCellLocalBounds(local[ci],plane);if(!b)return std::nullopt;
        const auto sx=fragmentExactDifference(plane.cx,anchor.cx),sy=fragmentExactDifference(plane.cy,anchor.cy);if(!sx||!sy)return std::nullopt;
        const auto x0=twoCellTranslatedBound(b->x0,*sx,false),x1=twoCellTranslatedBound(b->x1,*sx,true);
        const auto y0=twoCellTranslatedBound(b->y0,*sy,false),y1=twoCellTranslatedBound(b->y1,*sy,true);
        if(!x0||!x1||!y0||!y1)return std::nullopt;fragmentInclude(bounds,*x0,*x1,*y0,*y1);}
    if(!bounds.set||!(bounds.x1>bounds.x0&&bounds.y1>bounds.y0))return std::nullopt;
    FiniteVolume::Dyadic diameter(bounds.x1);diameter.add(bounds.x0,-1);diameter.add(bounds.y1);diameter.add(bounds.y0,-1);
    if(diameter.compare(FiniteVolume::Dyadic(delta))>=0)return std::nullopt;
    const auto rounded=FiniteVolume::round(diameter);auto upper=rounded.value;
    if(rounded.error.compare(FiniteVolume::Dyadic{})<0)upper=std::nextafter(upper,std::numeric_limits<long double>::infinity());
    if(!(std::isfinite(upper)&&upper>0&&upper<delta))return std::nullopt;
    return VertexStarProof{component,point,phase,bounds,delta,upper,faceCount,checkedCorners.size()};
}
}
struct GeometricVofPublication::IsolatedFragment::Data
{
    std::shared_ptr<const GeometricVofPublication::Data> owner;
    int native=-1;unsigned phase=0;long double diameter=0,delta=0;std::size_t faces=0,corners=0;
};
GeometricVofPublication::IsolatedFragment::IsolatedFragment(std::shared_ptr<const Data> data):data_(std::move(data)){}
GeometricVofIdentity GeometricVofPublication::IsolatedFragment::sourceIdentity() const
{require(bool(data_),"consumed isolated-fragment qualification");return data_->owner->identity;}
GeometricVofTimeLayer GeometricVofPublication::IsolatedFragment::sourceTime() const
{require(bool(data_),"consumed isolated-fragment qualification");return data_->owner->time;}
int GeometricVofPublication::IsolatedFragment::nativeCell() const
{require(bool(data_),"consumed isolated-fragment qualification");return data_->native;}
unsigned GeometricVofPublication::IsolatedFragment::phase() const
{require(bool(data_),"consumed isolated-fragment qualification");return data_->phase;}
long double GeometricVofPublication::IsolatedFragment::diameter() const
{require(bool(data_),"consumed isolated-fragment qualification");return data_->diameter;}
long double GeometricVofPublication::IsolatedFragment::delta() const
{require(bool(data_),"consumed isolated-fragment qualification");return data_->delta;}
std::size_t GeometricVofPublication::IsolatedFragment::coveredNativeFaces() const
{require(bool(data_),"consumed isolated-fragment qualification");return data_->faces;}
std::size_t GeometricVofPublication::IsolatedFragment::checkedCornerCells() const
{require(bool(data_),"consumed isolated-fragment qualification");return data_->corners;}
void GeometricVofPublication::IsolatedFragment::validateFor(const QuadtreeMesh& mesh,const GeometricVofPublication& source) const
{
    require(bool(data_)&&bool(data_->owner),"consumed isolated-fragment qualification");
    source.validateFor(mesh,data_->owner->identity,data_->owner->time);
    require(source.data_==data_->owner,"different immutable isolated-fragment owner");
    require(data_->native>=0&&static_cast<std::size_t>(data_->native)<source.checked().plicCells.size(),"isolated-fragment native index outside owner");
}
std::optional<GeometricVofPublication::IsolatedFragment> GeometricVofPublication::isolatedLocalPlicFragment(
    const QuadtreeMesh& mesh,int native) const
{
    validateMesh(mesh);const auto& data=checked();
    require(native>=0&&static_cast<std::size_t>(native)<data.cells.size(),"isolated-fragment native index outside publication");
    require(data.domain.cellRowByIndex()[native]>=0,"inactive native cell has no local PLIC geometry");
    if(data.kind!=GeometricVofGeometryKind::LocalPlic||data.domain.axisymmetric())return std::nullopt;
    if(data.plicCells.size()!=static_cast<std::size_t>(mesh.cellIndexNum())||
        data.cells.size()!=data.plicCells.size()||data.faces.size()!=static_cast<std::size_t>(data.domain.nativeFaceCount()))return std::nullopt;
    const auto& cell=data.plicCells[native];
    if(!cell.resolved||!(cell.measured[0]>0&&cell.measured[1]>0))return std::nullopt;
    for(unsigned phase=0;phase<2;++phase){
        const auto proof=fragmentProof(mesh,data.domain,data.plicCells,data.faces,native,phase);
        if(!proof)continue;
        auto evidence=std::make_shared<IsolatedFragment::Data>();evidence->owner=data_;
        evidence->native=native;evidence->phase=phase;evidence->diameter=proof->diameter;evidence->delta=proof->delta;
        evidence->faces=proof->faces;evidence->corners=proof->corners;
        return IsolatedFragment(std::move(evidence));
    }
    return std::nullopt;
}
struct GeometricVofPublication::TwoCellPositionDeficit::Data
{
    std::shared_ptr<const GeometricVofPublication::Data> owner;
    int native=-1;std::array<int,2> pair{};unsigned phase=0,axis=0;
    FragmentBounds bounds;long double delta=0,bound=0;std::size_t faces=0,corners=0;
};
GeometricVofPublication::TwoCellPositionDeficit::TwoCellPositionDeficit(std::shared_ptr<const Data> data):data_(std::move(data)){}
GeometricVofIdentity GeometricVofPublication::TwoCellPositionDeficit::sourceIdentity() const
{require(bool(data_),"consumed two-cell position-deficit qualification");return data_->owner->identity;}
GeometricVofTimeLayer GeometricVofPublication::TwoCellPositionDeficit::sourceTime() const
{require(bool(data_),"consumed two-cell position-deficit qualification");return data_->owner->time;}
int GeometricVofPublication::TwoCellPositionDeficit::nativeCell() const
{require(bool(data_),"consumed two-cell position-deficit qualification");return data_->native;}
const std::array<int,2>& GeometricVofPublication::TwoCellPositionDeficit::nativeCells() const
{require(bool(data_),"consumed two-cell position-deficit qualification");return data_->pair;}
unsigned GeometricVofPublication::TwoCellPositionDeficit::phase() const
{require(bool(data_),"consumed two-cell position-deficit qualification");return data_->phase;}
std::array<long double,2> GeometricVofPublication::TwoCellPositionDeficit::bboxLower() const
{require(bool(data_),"consumed two-cell position-deficit qualification");return {{data_->bounds.x0,data_->bounds.y0}};}
std::array<long double,2> GeometricVofPublication::TwoCellPositionDeficit::bboxUpper() const
{require(bool(data_),"consumed two-cell position-deficit qualification");return {{data_->bounds.x1,data_->bounds.y1}};}
long double GeometricVofPublication::TwoCellPositionDeficit::delta() const
{require(bool(data_),"consumed two-cell position-deficit qualification");return data_->delta;}
long double GeometricVofPublication::TwoCellPositionDeficit::twoBinL1Bound() const
{require(bool(data_),"consumed two-cell position-deficit qualification");return data_->bound;}
unsigned GeometricVofPublication::TwoCellPositionDeficit::splitAxis() const
{require(bool(data_),"consumed two-cell position-deficit qualification");return data_->axis;}
unsigned GeometricVofPublication::TwoCellPositionDeficit::maximumIndependentPositions() const
{require(bool(data_),"consumed two-cell position-deficit qualification");return 2;}
std::size_t GeometricVofPublication::TwoCellPositionDeficit::coveredNativeFaces() const
{require(bool(data_),"consumed two-cell position-deficit qualification");return data_->faces;}
std::size_t GeometricVofPublication::TwoCellPositionDeficit::checkedCornerCells() const
{require(bool(data_),"consumed two-cell position-deficit qualification");return data_->corners;}
void GeometricVofPublication::TwoCellPositionDeficit::validateFor(const QuadtreeMesh& mesh,const GeometricVofPublication& source) const
{
    require(bool(data_)&&bool(data_->owner),"consumed two-cell position-deficit qualification");
    source.validateFor(mesh,data_->owner->identity,data_->owner->time);
    require(source.data_==data_->owner,"different immutable two-cell position-deficit owner");
    const auto count=source.checked().plicCells.size();
    require(data_->pair[0]>=0&&data_->pair[0]<data_->pair[1]&&static_cast<std::size_t>(data_->pair[1])<count&&
        (data_->native==data_->pair[0]||data_->native==data_->pair[1]),"two-cell position-deficit native indices outside owner");
}
std::optional<GeometricVofPublication::TwoCellPositionDeficit> GeometricVofPublication::twoCellLocalPlicPositionDeficit(
    const QuadtreeMesh& mesh,int native) const
{
    validateMesh(mesh);const auto& data=checked();
    require(native>=0&&static_cast<std::size_t>(native)<data.cells.size(),"two-cell position-deficit native index outside publication");
    require(data.domain.cellRowByIndex()[native]>=0,"inactive native cell has no local PLIC geometry");
    if(data.kind!=GeometricVofGeometryKind::LocalPlic||data.domain.axisymmetric())return std::nullopt;
    if(data.plicCells.size()!=static_cast<std::size_t>(mesh.cellIndexNum())||data.cells.size()!=data.plicCells.size()||
        data.faces.size()!=static_cast<std::size_t>(data.domain.nativeFaceCount()))return std::nullopt;
    const auto& cell=data.plicCells[native];if(!cell.resolved||!(cell.measured[0]>0&&cell.measured[1]>0))return std::nullopt;
    std::set<int> neighbors;
    for(std::size_t row=0;row<data.faces.size();++row){if(data.faces[row].nativeRow!=static_cast<int>(row))return std::nullopt;
        const auto f=data.domain.face(mesh,static_cast<int>(row));
        if(f.negativeCell==native&&f.positiveCell>=0)neighbors.insert(f.positiveCell);
        if(f.positiveCell==native&&f.negativeCell>=0)neighbors.insert(f.negativeCell);}
    for(unsigned phase=0;phase<2;++phase)for(const int neighbor:neighbors){
        const std::array<int,2> pair{{std::min(native,neighbor),std::max(native,neighbor)}};
        const auto proof=twoCellProof(mesh,data.domain,data.plicCells,data.faces,pair,phase);if(!proof)continue;
        auto evidence=std::make_shared<TwoCellPositionDeficit::Data>();evidence->owner=data_;
        evidence->native=native;evidence->pair=proof->native;evidence->phase=proof->phase;evidence->axis=proof->splitAxis;
        evidence->bounds=proof->bounds;evidence->delta=proof->delta;evidence->bound=proof->bound;
        evidence->faces=proof->faces;evidence->corners=proof->corners;
        return TwoCellPositionDeficit(std::move(evidence));
    }
    return std::nullopt;
}
struct GeometricVofPublication::VertexStarFragment::Data
{
    std::shared_ptr<const GeometricVofPublication::Data> owner;
    int native=-1;std::vector<int> component;Point vertex;unsigned phase=0;
    FragmentBounds bounds;long double diameter=0,delta=0;std::size_t faces=0,corners=0;
};
GeometricVofPublication::VertexStarFragment::VertexStarFragment(std::shared_ptr<const Data> data):data_(std::move(data)){}
GeometricVofIdentity GeometricVofPublication::VertexStarFragment::sourceIdentity() const
{require(bool(data_),"consumed vertex-star qualification");return data_->owner->identity;}
GeometricVofTimeLayer GeometricVofPublication::VertexStarFragment::sourceTime() const
{require(bool(data_),"consumed vertex-star qualification");return data_->owner->time;}
int GeometricVofPublication::VertexStarFragment::nativeCell() const
{require(bool(data_),"consumed vertex-star qualification");return data_->native;}
const std::vector<int>& GeometricVofPublication::VertexStarFragment::nativeCells() const
{require(bool(data_),"consumed vertex-star qualification");return data_->component;}
Point GeometricVofPublication::VertexStarFragment::starVertex() const
{require(bool(data_),"consumed vertex-star qualification");return data_->vertex;}
unsigned GeometricVofPublication::VertexStarFragment::phase() const
{require(bool(data_),"consumed vertex-star qualification");return data_->phase;}
std::array<long double,2> GeometricVofPublication::VertexStarFragment::bboxLower() const
{require(bool(data_),"consumed vertex-star qualification");return {{data_->bounds.x0,data_->bounds.y0}};}
std::array<long double,2> GeometricVofPublication::VertexStarFragment::bboxUpper() const
{require(bool(data_),"consumed vertex-star qualification");return {{data_->bounds.x1,data_->bounds.y1}};}
long double GeometricVofPublication::VertexStarFragment::diameter() const
{require(bool(data_),"consumed vertex-star qualification");return data_->diameter;}
long double GeometricVofPublication::VertexStarFragment::delta() const
{require(bool(data_),"consumed vertex-star qualification");return data_->delta;}
std::size_t GeometricVofPublication::VertexStarFragment::coveredNativeFaces() const
{require(bool(data_),"consumed vertex-star qualification");return data_->faces;}
std::size_t GeometricVofPublication::VertexStarFragment::checkedCornerCells() const
{require(bool(data_),"consumed vertex-star qualification");return data_->corners;}
void GeometricVofPublication::VertexStarFragment::validateFor(const QuadtreeMesh& mesh,const GeometricVofPublication& source) const
{
    require(bool(data_)&&bool(data_->owner),"consumed vertex-star qualification");
    source.validateFor(mesh,data_->owner->identity,data_->owner->time);
    require(source.data_==data_->owner,"different immutable vertex-star owner");
    const auto count=source.checked().plicCells.size();
    require(!data_->component.empty()&&data_->component.size()<=4&&
        std::is_sorted(data_->component.begin(),data_->component.end())&&
        std::adjacent_find(data_->component.begin(),data_->component.end())==data_->component.end()&&
        data_->component.front()>=0&&static_cast<std::size_t>(data_->component.back())<count&&
        std::binary_search(data_->component.begin(),data_->component.end(),data_->native),"vertex-star native indices outside owner");
}
std::optional<GeometricVofPublication::VertexStarFragment> GeometricVofPublication::vertexStarLocalPlicFragment(
    const QuadtreeMesh& mesh,int native) const
{
    validateMesh(mesh);const auto& data=checked();
    require(native>=0&&static_cast<std::size_t>(native)<data.cells.size(),"vertex-star native index outside publication");
    require(data.domain.cellRowByIndex()[native]>=0,"inactive native cell has no local PLIC geometry");
    if(data.kind!=GeometricVofGeometryKind::LocalPlic||data.domain.axisymmetric())return std::nullopt;
    if(data.plicCells.size()!=static_cast<std::size_t>(mesh.cellIndexNum())||data.cells.size()!=data.plicCells.size()||
        data.faces.size()!=static_cast<std::size_t>(data.domain.nativeFaceCount()))return std::nullopt;
    const auto& cell=data.plicCells[native];if(!cell.resolved||!(cell.measured[0]>0&&cell.measured[1]>0))return std::nullopt;
    std::set<std::size_t> vertices;
    for(std::size_t row=0;row<data.faces.size();++row){if(data.faces[row].nativeRow!=static_cast<int>(row))return std::nullopt;
        const auto face=data.domain.face(mesh,static_cast<int>(row));
        if(face.negativeCell==native||face.positiveCell==native)vertices.insert(data.faces[row].vertices.begin(),data.faces[row].vertices.end());}
    if(vertices.size()<4||vertices.size()>8)return std::nullopt;
    for(unsigned phase=0;phase<2;++phase)for(const auto vertex:vertices){
        const auto proof=vertexStarProof(mesh,data.domain,data.plicCells,data.faces,native,vertex,phase);if(!proof)continue;
        auto evidence=std::make_shared<VertexStarFragment::Data>();evidence->owner=data_;
        evidence->native=native;evidence->component=proof->native;evidence->vertex=proof->vertex;evidence->phase=proof->phase;
        evidence->bounds=proof->bounds;evidence->diameter=proof->diameter;evidence->delta=proof->delta;
        evidence->faces=proof->faces;evidence->corners=proof->corners;return VertexStarFragment(std::move(evidence));
    }
    return std::nullopt;
}
void GeometricVofPublication::validateMesh(const QuadtreeMesh& mesh) const{checked().domain.validateMeshIdentity(mesh);}
void GeometricVofPublication::validateFor(const QuadtreeMesh& mesh,GeometricVofIdentity owner,GeometricVofTimeLayer time) const
{ const auto& data=checked();time.time();require(owner==data.identity&&owner.value()!=0,"foreign VOF owner");require(sameClock(time,data.time),"full clock mismatch");data.domain.validateMeshIdentity(mesh); }
void GeometricVofPublication::validateDomainInBatch(const QuadtreeMesh& mesh,
    const NativeLeafDomain& domain,const MeshValidationBatch& batch)
{
    // Equality here authenticates the same immutable construction witness,
    // not merely equal mesh contents. Other bindings keep the original full
    // ordered comparison; equal-content independent meshes remain legal.
    if(batch.mesh_==&mesh&&batch.identity_&&domain.meshIdentity()==batch.identity_)return;
    domain.validateMeshIdentity(mesh);
}
void GeometricVofPublication::validateForInBatch(const QuadtreeMesh& mesh,
    GeometricVofIdentity owner,GeometricVofTimeLayer time,const MeshValidationBatch& batch) const
{
    const auto& data=checked();time.time();
    require(owner==data.identity&&owner.value()!=0,"foreign VOF owner");
    require(sameClock(time,data.time),"full clock mismatch");
    validateDomainInBatch(mesh,data.domain,batch);
}
}





