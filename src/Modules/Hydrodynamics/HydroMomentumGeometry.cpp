#include "Modules/Hydrodynamics/HydroMomentumGeometry.h"
#include "Modules/Hydrodynamics/Detail/HydroFiniteContraction.h"

#include "Modules/Interface/HermiteGalsPhaseRegion.h"
#include "Modules/Hydrodynamics/HydroSpatialAdvector.h"
#include "Modules/Hydrodynamics/Detail/HydroRegularAdvectorGeometry.h"
#include "Core/Geometry/LevelSetMappedSegment.h"
#include "Core/Geometry/Detail/ExactExpansion.h"
#include "Core/Geometry/Detail/FiniteVolumeClosure.h"
#include <Eigen/SparseCholesky>
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/CartesianMeasure.h"
#include <algorithm>
#include <cmath>
#include <chrono>
#include <cfenv>
#include <cerrno>
#if defined(_M_X64) || defined(__SSE2__)
#include <xmmintrin.h>
#endif
#include <exception>
#include <limits>
#include <map>
#include <mutex>
#include "Core/Mesh/NativeLeafOverlap.h"
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <tuple>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace LAMP::Hydrodynamics
{
namespace
{
void require(bool ok, const char* message)
{
    if (!ok) throw std::invalid_argument(std::string("[HydroMomentumGeometry] ")+message);
}
double coord(Vector2D p, int axis) { return axis==0 ? p.x : p.y; }
bool samePoint(Vector2D a, Vector2D b) { return a.x==b.x && a.y==b.y; }
Vector2D onPlane(int axis, double plane, double tangent)
{ return axis==0 ? Vector2D{plane,tangent} : Vector2D{tangent,plane}; }
long double capacity(Vector2D a, Vector2D b, bool rz)
{
    return (static_cast<long double>(b.x)-a.x)*(static_cast<long double>(b.y)-a.y)
        *(rz ? .5L*CartesianMeasure::twoPi*(static_cast<long double>(b.x)+a.x) : 1.L);
}
long double relative(long double a, long double b)
{ return std::abs(a-b)/std::max(std::abs(a),std::abs(b)); }
constexpr long double geometryTolerance=512.L*std::numeric_limits<double>::epsilon();
// Exact physical rectangle primitive; the common RZ pi cancels in a ratio.
// Expand actual published binary coordinates, not rounded widths/volumes.
void addP0Rectangle(Detail::HydroFiniteContraction& value,Vector2D lo,Vector2D hi,bool rz)
{
    if(rz){
        value.addProduct({hi.x,hi.x,hi.y});value.addProduct({-static_cast<long double>(hi.x),hi.x,lo.y});
        value.addProduct({-static_cast<long double>(lo.x),lo.x,hi.y});value.addProduct({lo.x,lo.x,lo.y});
    }else{
        value.addProduct({hi.x,hi.y});value.addProduct({-static_cast<long double>(hi.x),lo.y});
        value.addProduct({-static_cast<long double>(lo.x),hi.y});value.addProduct({lo.x,lo.y});
    }
}
const Cell& cellAt(const QuadtreeMesh& mesh, int index)
{ return mesh.cells(mesh.cellIndexToLabel[index]); }
const Edge& edgeAt(const QuadtreeMesh& mesh, const NativeLeafFace& face)
{ return face.horizontal ? mesh.hEdges(face.edgeLabel) : mesh.vEdges(face.edgeLabel); }
// Local ownership bookkeeping only. At most two actual subfaces per side.
struct Piece { int velocity=-1; Vector2D lower{0,0},upper{0,0}; };
struct HalfSide { std::array<Piece,2> pieces; int count=0; };
using CellPieces=std::array<std::array<HalfSide,2>,2>; // component, normal side
class Builder
{
public:
    Builder(const QuadtreeMesh& mesh, const HydroPhysicalDivergence& metrics,
        std::vector<HydroMomentumVolume>& volumes, std::vector<HydroMomentumSegment>& segments,
        HydroMomentumGeometryReport& report)
        : mesh_(mesh), metrics_(metrics), layout_(metrics.layout()), volumes_(volumes),
          segments_(segments), report_(report), pieces_(layout_.pressureCount()) {}
    void run()
    {
        layout_.validateMesh(mesh_);
        volumes_.resize(layout_.velocityCount());
        makePieces(); insideCells(); acrossPrimalFaces(); validateCoverage();
    }
private:
    const QuadtreeMesh& mesh_;
    const HydroPhysicalDivergence& metrics_;
    const HydroDofLayout& layout_;
    std::vector<HydroMomentumVolume>& volumes_;
    std::vector<HydroMomentumSegment>& segments_;
    HydroMomentumGeometryReport& report_;
    std::vector<CellPieces> pieces_;
    const CellPieces& pieces(int native) const { return pieces_[layout_.pressureDof(native)]; }
    void append(int negative, int positive, int normal, double plane, double low, double high,
        int primal, int owner)
    {
        require(low<high && std::isfinite(low) && std::isfinite(high) && std::isfinite(plane),
            "invalid common segment endpoints");
        require(negative!=positive && (negative>=0 || positive>=0)
            && (primal>=0)!=(owner>=0),"invalid common segment origin");
        if (negative>=0 && positive>=0)
            require(layout_.normalAxis(negative)==layout_.normalAxis(positive),
                "common segment connects different velocity components");
        HydroMomentumSegment s;
        s.negativeVelocity=negative; s.positiveVelocity=positive; s.normalAxis=normal;
        s.begin=onPlane(normal,plane,low); s.end=onPlane(normal,plane,high);
        s.primalVelocity=primal; s.ownerPressure=owner;
        s.area=(static_cast<long double>(high)-low)*(layout_.domain().axisymmetric()
            ? .5L*CartesianMeasure::twoPi*(static_cast<long double>(s.begin.x)+s.end.x) : 1.L);
        require(std::isfinite(s.area) && (s.area>0 || (s.area==0
            && layout_.domain().axisymmetric() && normal==0 && plane==0)),"invalid common segment area");
        const int id=static_cast<int>(segments_.size()); segments_.push_back(s);
        if (negative>=0) volumes_[negative].sides[2*normal+1].push_back(id);
        if (positive>=0) volumes_[positive].sides[2*normal].push_back(id);
        if (negative>=0 && positive>=0) ++report_.sharedSegments; else ++report_.boundarySegments;
    }
    void makePieces()
    {
        for (int p=0;p<layout_.pressureCount();++p) {
            const int nativeCell=layout_.nativeCell(p);
            const auto& cell=cellAt(mesh_,nativeCell);
            const auto low=mesh_.verts(cell.cellToVert[0]).coord();
            const auto high=mesh_.verts(cell.cellToVert[2]).coord();
            const auto center=cell.centerCoord();
            const long double cellVolume=capacity(low,high,layout_.domain().axisymmetric());
            require(std::isfinite(cellVolume) && cellVolume>0 &&
                relative(cellVolume,metrics_.cellVolumes()[p])<=geometryTolerance,"invalid cell volume");
            for (int component=0;component<2;++component) {
                long double total=0;
                for (int side=0;side<2;++side) {
                    auto& half=pieces_[p][component][side];
                    for (int part=0;part<2;++part) {
                        const int label=component==0 ? cell.cellToVEdge[side][part] : cell.cellToHEdge[side][part];
                        if (label<0) continue;
                        const auto& edge=component==0 ? mesh_.vEdges(label) : mesh_.hEdges(label);
                        const int native=edge.index()+(component==0 ? 0 : mesh_.vEdgeIndexNum());
                        const int u=layout_.velocityDof(native);
                        require(u>=0,"missing selected native velocity");
                        const auto face=layout_.domain().face(mesh_,native);
                        require((side==0 ? face.positiveCell : face.negativeCell)==nativeCell,
                            "half-cell owner differs from native face");
                        auto begin=mesh_.verts(edge.edgeToVert[0]).coord(),end=mesh_.verts(edge.edgeToVert[1]).coord();
                        if (coord(begin,1-component)>coord(end,1-component)) std::swap(begin,end);
                        auto& piece=half.pieces[half.count++]; piece.velocity=u;
                        piece.lower=component==0 ? Vector2D{side==0?low.x:center.x,begin.y}
                            : Vector2D{begin.x,side==0?low.y:center.y};
                        piece.upper=component==0 ? Vector2D{side==0?center.x:high.x,end.y}
                            : Vector2D{end.x,side==0?center.y:high.y};
                        const auto v=capacity(piece.lower,piece.upper,layout_.domain().axisymmetric());
                        require(std::isfinite(v) && v>0,"invalid half-cell volume");
                        auto& volume=volumes_[u];
                        if (volume.overlaps.empty()) { volume.lower=piece.lower; volume.upper=piece.upper; }
                        else {
                            require(volume.overlaps.size()==1 &&
                                coord(volume.lower,1-component)==coord(piece.lower,1-component) &&
                                coord(volume.upper,1-component)==coord(piece.upper,1-component),
                                "native owners disagree on tangential coverage");
                            volume.lower={std::min(volume.lower.x,piece.lower.x),std::min(volume.lower.y,piece.lower.y)};
                            volume.upper={std::max(volume.upper.x,piece.upper.x),std::max(volume.upper.y,piece.upper.y)};
                        }
                        volume.overlaps.push_back({p,piece.lower,piece.upper,v,v/cellVolume});
                        volume.volume+=v; total+=v; ++report_.overlaps;
                    }
                    require(half.count>0,"missing native component face");
                    std::sort(half.pieces.begin(),half.pieces.begin()+half.count,[component](const Piece& a,const Piece& b) {
                        return coord(a.lower,1-component)<coord(b.lower,1-component); });
                    double cursor=coord(low,1-component);
                    for (int part=0;part<half.count;++part) {
                        const auto& piece=half.pieces[part];
                        require(coord(piece.lower,1-component)==cursor,"gap or overlap in half-cell pieces");
                        cursor=coord(piece.upper,1-component);
                    }
                    require(cursor==coord(high,1-component),"incomplete half-cell side");
                }
                const auto defect=relative(total,cellVolume);
                report_.maxRelativeCellPartition=std::max(report_.maxRelativeCellPartition,defect);
                require(defect<=geometryTolerance,"momentum pieces do not partition the cell");
            }
        }
    }
    void insideCells()
    {
        for (int p=0;p<layout_.pressureCount();++p) {
            const auto center=cellAt(mesh_,layout_.nativeCell(p)).centerCoord();
            for (int component=0;component<2;++component) {
                const auto& negative=pieces_[p][component][0]; const auto& positive=pieces_[p][component][1];
                for (int i=0;i<negative.count;++i) for (int j=0;j<positive.count;++j) {
                    const auto& a=negative.pieces[i]; const auto& b=positive.pieces[j];
                    const double lo=std::max(coord(a.lower,1-component),coord(b.lower,1-component));
                    const double hi=std::min(coord(a.upper,1-component),coord(b.upper,1-component));
                    if (lo<hi) append(a.velocity,b.velocity,component,coord(center,component),lo,hi,-1,p);
                }
                for (int side=0;side<2;++side) {
                    const auto& half=pieces_[p][component][side]; if (half.count!=2) continue;
                    const auto& a=half.pieces[0]; const auto& b=half.pieces[1];
                    append(a.velocity,b.velocity,1-component,coord(a.upper,1-component),
                        coord(a.lower,component),coord(a.upper,component),-1,p);
                }
            }
        }
    }
    struct Interval { int velocity; double lower,upper; };
    struct Intervals { std::array<Interval,4> values; int count=0; };
    Intervals touching(const NativeLeafFace& cross,int owner,int side) const
    {
        Intervals result; if (owner<0) return result;
        const int normal=cross.horizontal?1:0,component=1-normal;
        const auto& edge=edgeAt(mesh_,cross);
        const double begin=coord(mesh_.verts(edge.edgeToVert[0]).coord(),component);
        const double end=coord(mesh_.verts(edge.edgeToVert[1]).coord(),component);
        const double lo=std::min(begin,end),hi=std::max(begin,end),plane=coord(cross.center,normal);
        for (const auto& half:pieces(owner)[component]) for (int i=0;i<half.count;++i) {
            const auto& piece=half.pieces[i];
            if (coord(side==0 ? piece.lower : piece.upper,normal)!=plane) continue;
            const double a=std::max(lo,coord(piece.lower,component)),b=std::min(hi,coord(piece.upper,component));
            if (a<b) { require(result.count<4,"too many touching intervals"); result.values[result.count++]={piece.velocity,a,b}; }
        }
        require(result.count>0,"no momentum piece touches primal face"); return result;
    }
    void acrossPrimalFaces()
    {
        for (int u=0;u<layout_.velocityCount();++u) {
            const auto face=layout_.domain().face(mesh_,layout_.nativeFace(u));
            const int normal=face.horizontal?1:0; const double plane=coord(face.center,normal);
            const auto negative=touching(face,face.negativeCell,1),positive=touching(face,face.positiveCell,0);
            if (negative.count>0 && positive.count>0) {
                for (int i=0;i<negative.count;++i) for (int j=0;j<positive.count;++j) {
                    const auto& a=negative.values[i]; const auto& b=positive.values[j];
                    const double lo=std::max(a.lower,b.lower),hi=std::min(a.upper,b.upper);
                    if (lo<hi) append(a.velocity,b.velocity,normal,plane,lo,hi,u,-1);
                }
            } else {
                const bool outwardPositive=face.negativeCell>=0;
                const auto& intervals=outwardPositive ? negative : positive;
                for (int i=0;i<intervals.count;++i) {
                    const auto& part=intervals.values[i];
                    append(outwardPositive?part.velocity:-1,outwardPositive?-1:part.velocity,
                        normal,plane,part.lower,part.upper,u,-1);
                }
                const auto& edge=edgeAt(mesh_,face);
                const double a=coord(mesh_.verts(edge.edgeToVert[0]).coord(),1-normal);
                const double b=coord(mesh_.verts(edge.edgeToVert[1]).coord(),1-normal);
                append(outwardPositive?u:-1,outwardPositive?-1:u,normal,plane,std::min(a,b),std::max(a,b),u,-1);
            }
        }
    }
    void validateCoverage()
    {
        std::vector<int> uses(segments_.size());
        for (int u=0;u<layout_.velocityCount();++u) {
            auto& volume=volumes_[u];
            require(!volume.overlaps.empty() && volume.overlaps.size()<=2,"invalid dual owner count");
            const auto box=capacity(volume.lower,volume.upper,layout_.domain().axisymmetric());
            const auto defect=relative(volume.volume,box);
            report_.maxRelativeBoxClosure=std::max(report_.maxRelativeBoxClosure,defect);
            require(std::isfinite(volume.volume) && box>0 && defect<=geometryTolerance,"dual volume has gap or overlap");
            for (int side=0;side<4;++side) {
                const int normal=side/2,tangent=1-normal;
                const double plane=coord(side%2==0?volume.lower:volume.upper,normal);
                auto& ids=volume.sides[side];
                std::sort(ids.begin(),ids.end(),[&](int a,int b){return coord(segments_[a].begin,tangent)<coord(segments_[b].begin,tangent);});
                double cursor=coord(volume.lower,tangent);
                for (int id:ids) {
                    const auto& s=segments_[id];
                    require(s.normalAxis==normal && coord(s.begin,normal)==plane && coord(s.end,normal)==plane
                        && coord(s.begin,tangent)==cursor,"gap, overlap or misplaced common segment");
                    require((side%2==0?s.positiveVelocity:s.negativeVelocity)==u,"wrong outward segment orientation");
                    cursor=coord(s.end,tangent); ++uses[id];
                }
                require(cursor==coord(volume.upper,tangent),"incomplete dual side");
            }
        }
        for (std::size_t id=0;id<segments_.size();++id) {
            const auto& s=segments_[id];
            require(uses[id]==(s.negativeVelocity>=0?1:0)+(s.positiveVelocity>=0?1:0),"missing reciprocal common segment");
        }
    }
};
Interface::HermiteGalsPhaseRegionConfig regionConfig(const HydroPhaseGeometry& phase)
{
    Interface::HermiteGalsPhaseRegionConfig config;
    const auto& snapshot=phase.publication().state().geometry().snapshot();
    for (const auto& curve:phase.publication().topology().curves) {
        require(curve.usesSource(snapshot),"foreign retained curve scalar");
        for (const auto& point:curve.boundaryContacts())
            if (std::none_of(config.boundaryContacts.begin(),config.boundaryContacts.end(),[point](Vector2D other){return samePoint(point,other);}))
                config.boundaryContacts.push_back(point);
        for (const auto& plane:curve.symmetryPlanes())
            if (std::find(config.symmetryPlanes.begin(),config.symmetryPlanes.end(),plane)==config.symmetryPlanes.end())
                config.symmetryPlanes.push_back(plane);
    }
    return config;
}
}
struct HydroMomentumGeometry::Data
{
    Data(HydroPhysicalDivergence m,std::shared_ptr<Detail::HydroRegularOperatorTemplateBank> t)
        :templates(std::move(t)),sparsePlans(makeSparseValidationPlans()),frozenMatrixPlans(makeFrozenMatrixAssemblyPlans()),noSlipRaw(makeStationaryNoSlipStrainCache()),metrics(std::move(m)) {}
    std::shared_ptr<Detail::HydroRegularOperatorTemplateBank> templates;
    std::shared_ptr<Detail::HydroSparseValidationPlanCache> sparsePlans;
    std::shared_ptr<Detail::HydroFrozenMatrixPlanCache> frozenMatrixPlans;
    std::shared_ptr<Detail::HydroStationaryNoSlipStrainCache> noSlipRaw; // Fresh target Data starts empty.
    mutable std::mutex openConvectionMutex;
    mutable bool openConvectionAttempted=false;
    mutable std::shared_ptr<const Detail::HydroOpenConvectionSupport> openConvection;
    HydroPhysicalDivergence metrics;
    std::vector<HydroMomentumVolume> volumes;
    std::vector<HydroMomentumSegment> segments;
    HydroMomentumGeometryReport report;
    mutable std::mutex regularRemapMutex;
    mutable std::shared_ptr<const Detail::HydroRegularRemapPreparation> regularRemap;
    mutable HydroRegularPreparationReuseReport regularReuse;
    mutable bool regularRemapPending=false;
};
HydroMomentumGeometry::HydroMomentumGeometry(std::shared_ptr<const Data> d):data_(std::move(d)) {}
HydroMomentumGeometry HydroMomentumGeometry::build(const QuadtreeMesh& mesh,HydroPhysicalDivergence metrics)
{
    auto d=std::make_shared<Data>(std::move(metrics),makeRegularTemplates());
    Builder(mesh,d->metrics,d->volumes,d->segments,d->report).run(); return HydroMomentumGeometry(std::move(d));
}
HydroMomentumGeometry HydroMomentumGeometry::buildWithRegularTemplates(const QuadtreeMesh& mesh,
    HydroPhysicalDivergence metrics,const HydroMomentumGeometry& sourceTemplateOwner)
{
    auto d=std::make_shared<Data>(std::move(metrics),sourceTemplateOwner.regularTemplates());
    Builder(mesh,d->metrics,d->volumes,d->segments,d->report).run(); return HydroMomentumGeometry(std::move(d));
}
HydroMomentumGeometry HydroMomentumGeometry::buildWithRemappedRegularPreparation(
    const QuadtreeMesh& sourceMesh,const Detail::HydroRegularAdvectorGeometry& source,
    const QuadtreeMesh& targetMesh,HydroPhysicalDivergence metrics,const Mesh::NativeLeafOverlap& overlaps)
{
    auto d=std::make_shared<Data>(std::move(metrics),source.momentumGeometry().regularTemplates());
    Builder(targetMesh,d->metrics,d->volumes,d->segments,d->report).run();
    d->regularReuse.requested=true;d->regularRemapPending=true;
    d->regularRemap=Detail::HydroRegularAdvectorGeometry::makeRemapPreparation(sourceMesh,source,targetMesh,overlaps);
    return HydroMomentumGeometry(std::move(d));
}
HydroRegularPreparationReuseReport HydroMomentumGeometry::regularPreparationReuseReport() const
{std::lock_guard<std::mutex> lock(data_->regularRemapMutex);return data_->regularReuse;}
std::shared_ptr<const Detail::HydroRegularRemapPreparation> HydroMomentumGeometry::takeRegularRemapPreparation(bool& requested) const
{std::lock_guard<std::mutex> lock(data_->regularRemapMutex);requested=data_->regularRemapPending;data_->regularRemapPending=false;auto result=std::move(data_->regularRemap);data_->regularRemap.reset();return result;}
void HydroMomentumGeometry::publishRegularPreparationReuseReport(const HydroRegularPreparationReuseReport& report) const
{std::lock_guard<std::mutex> lock(data_->regularRemapMutex);data_->regularReuse=report;}
const HydroPhysicalDivergence& HydroMomentumGeometry::metrics() const noexcept { return data_->metrics; }
const HydroDofLayout& HydroMomentumGeometry::layout() const noexcept { return metrics().layout(); }
const std::vector<HydroMomentumVolume>& HydroMomentumGeometry::volumes() const noexcept { return data_->volumes; }
const std::vector<HydroMomentumSegment>& HydroMomentumGeometry::segments() const noexcept { return data_->segments; }
const HydroMomentumGeometryReport& HydroMomentumGeometry::report() const noexcept { return data_->report; }
bool HydroMomentumGeometry::sameIdentity(const HydroMomentumGeometry& other) const noexcept { return data_==other.data_; }
std::shared_ptr<Detail::HydroRegularOperatorTemplateBank> HydroMomentumGeometry::regularTemplates() const noexcept
{return data_->templates;}
std::shared_ptr<Detail::HydroSparseValidationPlanCache> HydroMomentumGeometry::sparseValidationPlans() const noexcept
{return data_->sparsePlans;}
std::shared_ptr<Detail::HydroFrozenMatrixPlanCache> HydroMomentumGeometry::frozenMatrixAssemblyPlans() const noexcept
{return data_->frozenMatrixPlans;}
std::shared_ptr<Detail::HydroStationaryNoSlipStrainCache> HydroMomentumGeometry::stationaryNoSlipStrainCache() const noexcept
{return data_->noSlipRaw;}
std::shared_ptr<const Detail::HydroOpenConvectionSupport> HydroMomentumGeometry::openConvectionSupport(const QuadtreeMesh& mesh) const noexcept
{
    // Optional integer preparation is invisible to the caller's numeric state.
    const int error=errno;std::fenv_t environment{};
    const bool saved=std::fegetenv(&environment)==0;
#if defined(_M_X64) || defined(__SSE2__)
    const unsigned csr=_mm_getcsr();
#endif
    const auto restore=[&]() noexcept {if(saved&&std::fesetenv(&environment)!=0)std::terminate();
#if defined(_M_X64) || defined(__SSE2__)
        _mm_setcsr(csr);
#endif
        errno=error;};
    if(!saved){restore();return {};}
    std::shared_ptr<const Detail::HydroOpenConvectionSupport> result;
    try{
        std::lock_guard<std::mutex> lock(data_->openConvectionMutex);
        if(!data_->openConvectionAttempted){
            data_->openConvectionAttempted=true;
            const auto n=static_cast<std::size_t>(layout().velocityCount()),ns=segments().size();
            constexpr std::size_t limit=16u*1024u*1024u;
            const auto cap=std::numeric_limits<std::uint32_t>::max();
            // Retain the conservative array/map-node bound for the smaller flat scratch.
            if(n>=cap||ns>=cap/2||n>limit/64||ns>limit/160)throw std::bad_alloc{};
            if(sizeof(Detail::HydroOpenConvectionSupport)+n*48+ns*128>limit)throw std::bad_alloc{};
            auto out=std::make_shared<Detail::HydroOpenConvectionSupport>();
            out->segmentPairSlot.assign(ns,-1);out->boundaryRowSlot.assign(n,-1);
            out->pairOffsets.resize(n+1);out->boundaryOffsets.resize(n+1);
            {
                // Nonnegative 32-bit rows pack in the original (min,max) lexicographic
                // order. Segment-id tie order is explicit; no numerical data is sorted.
                std::vector<std::pair<std::uint64_t,std::uint32_t>> routed;routed.reserve(ns);
                for(std::size_t id=0;id<ns;++id){const auto& g=segments()[id];
                    if(g.negativeVelocity>=0&&g.positiveVelocity>=0){
                        const auto ends=std::minmax(g.negativeVelocity,g.positiveVelocity);
                        const auto key=(static_cast<std::uint64_t>(ends.first)<<32)|static_cast<std::uint32_t>(ends.second);
                        routed.emplace_back(key,static_cast<std::uint32_t>(id));
                    }else ++out->boundaryOffsets[static_cast<std::size_t>(g.negativeVelocity>=0?g.negativeVelocity:g.positiveVelocity)+1];
                }
                std::sort(routed.begin(),routed.end());
                std::uint64_t previous=0;int slot=-1;
                for(const auto& item:routed){
                    if(slot<0||item.first!=previous){
                        const int first=static_cast<int>(item.first>>32),second=static_cast<int>(static_cast<std::uint32_t>(item.first));
                        slot=static_cast<int>(out->pairs.size());out->pairs.push_back({first,second});previous=item.first;
                        ++out->pairOffsets[static_cast<std::size_t>(first)+1];++out->pairOffsets[static_cast<std::size_t>(second)+1];
                    }
                    out->segmentPairSlot[item.second]=slot;
                }
            } // Release flat temporary before allocating the CSR items/cursors.
            for(std::size_t i=0;i<n;++i){out->pairOffsets[i+1]+=out->pairOffsets[i];out->boundaryOffsets[i+1]+=out->boundaryOffsets[i];
                if(out->boundaryOffsets[i+1]!=out->boundaryOffsets[i]){out->boundaryRowSlot[i]=static_cast<int>(out->boundaryRows.size());out->boundaryRows.push_back(static_cast<int>(i));}
                if(layout().domain().face(mesh,layout().nativeFace(static_cast<int>(i))).outwardSign!=0)out->exteriorVelocityDofs.push_back(static_cast<int>(i));}
            out->pairItems.resize(out->pairOffsets.back());out->boundaryItems.resize(out->boundaryOffsets.back());
            auto cursor=out->pairOffsets;
            for(std::size_t i=0;i<out->pairs.size();++i){const auto& p=out->pairs[i];
                out->pairItems[cursor[p.first]++]={static_cast<std::uint32_t>(i),1};out->pairItems[cursor[p.second]++]={static_cast<std::uint32_t>(i),-1};}
            cursor=out->boundaryOffsets;
            for(std::size_t i=0;i<ns;++i){const auto& g=segments()[i];if(g.negativeVelocity<0||g.positiveVelocity<0)
                out->boundaryItems[cursor[g.negativeVelocity>=0?g.negativeVelocity:g.positiveVelocity]++]=static_cast<std::uint32_t>(i);}
            out->bytes=sizeof(*out)+out->pairs.capacity()*sizeof(Detail::HydroOpenConvectionSupport::Pair)
                +(out->segmentPairSlot.capacity()+out->boundaryRowSlot.capacity()+out->exteriorVelocityDofs.capacity()+out->boundaryRows.capacity())*sizeof(int)
                +(out->pairOffsets.capacity()+out->boundaryOffsets.capacity()+out->boundaryItems.capacity())*sizeof(std::uint32_t)
                +out->pairItems.capacity()*sizeof(Detail::HydroOpenConvectionSupport::Incidence);
            if(out->bytes>limit)throw std::bad_alloc{};data_->openConvection=std::move(out);
        }
        result=data_->openConvection;
    }catch(...){result.reset();}
    restore();return result;
}
void HydroMomentumGeometry::validateMesh(const QuadtreeMesh& mesh) const { layout().validateMesh(mesh); }
struct HydroMomentumMass::Data
{
    Data(HydroMomentumGeometry g,HydroPhaseGeometry p,double l,double v)
        :geometry(std::move(g)),phase(std::move(p)),liquid(l),gas(v) {}
    HydroMomentumGeometry geometry;
    HydroPhaseGeometry phase;
    double liquid,gas;
    DensityRepresentation representation=DensityRepresentation::SharpGeometry;
    std::vector<HydroPhaseMeasure> phases;
    std::vector<std::vector<HydroPhaseMeasure>> overlaps;
    std::vector<double> masses;
    HydroMomentumMassReport report;
    std::array<std::vector<std::array<bool,2>>,2> finiteP0Partitions;
};
HydroMomentumMass::HydroMomentumMass(std::shared_ptr<const Data> d):data_(std::move(d)) {}
HydroMomentumMass HydroMomentumMass::build(const QuadtreeMesh& mesh,HydroMomentumGeometry geometry,
    HydroPhaseGeometry phase,double liquidDensity,double gasDensity)
{
    if(phase.model()==HydroPhaseGeometry::Model::GeometricVof)
        return buildVolumeFractionGeometry(mesh,std::move(geometry),std::move(phase),liquidDensity,gasDensity);
    const auto sourceStart=std::chrono::steady_clock::now();
    geometry.validateMesh(mesh); phase.validateFor(mesh,phase.geometryIdentity(),phase.timeLayer());
    require(geometry.layout().sameIdentity(phase.layout()),"different layout identity");
    require(&geometry.metrics().divergence().matrix()==&phase.metrics().divergence().matrix(),"different retained D_I");
    require(std::isfinite(liquidDensity) && liquidDensity>0 && std::isfinite(gasDensity) && gasDensity>0,"invalid material density");
    LAMP::Geometry::Detail::ExactExpansion::validateEnvironment();
    auto d=std::make_shared<Data>(std::move(geometry),std::move(phase),liquidDensity,gasDensity);
    const Interface::HermiteGalsPhaseRegionSource region(mesh,d->phase.publication(),regionConfig(d->phase));
    require(&region.source().geometry().snapshot()==&d->phase.publication().state().geometry().snapshot(),"different region scalar");
    const auto batchStart=std::chrono::steady_clock::now();
    const auto& layout=d->geometry.layout();
    // Reuse this actual publication's complete topology and one checked
    // immutable-source batch. Each query retains its original leaf clipping,
    // direct moment arithmetic, actual-thread FP checks and error slot.
    std::vector<std::array<Vector2D,2>> boxes;
    boxes.reserve(layout.velocityCount());
    for(int u=0;u<layout.velocityCount();++u){const auto& volume=d->geometry.volumes()[u];
        boxes.push_back({{volume.lower,volume.upper}});}
    const auto measuredVolumes=region.measureBatch(mesh,boxes);
    const auto assemblyStart=std::chrono::steady_clock::now();
    std::array<std::vector<HydroPhaseMeasure>,2> cellSums;
    for (auto& sum:cellSums) sum.resize(layout.pressureCount());
    for (int u=0;u<layout.velocityCount();++u) {
        // Preserve both the original first failure and every addition's order.
        if(measuredVolumes[u].failure)std::rethrow_exception(measuredVolumes[u].failure);
        const auto& volume=d->geometry.volumes()[u];
        const auto& measured=measuredVolumes[u].measure;
        require(measured.sourceIdentity==d->phase.geometryIdentity() && relative(measured.capacity,volume.volume)<=geometryTolerance,
            "dual phase capacity or source mismatch");
        require(measured.leaves.size()==volume.overlaps.size(),"dual phase has different native overlaps");
        std::vector<HydroPhaseMeasure> overlap;
        const int component=layout.normalAxis(u)==HydroDofLayout::NormalAxis::First?0:1;
        for (const auto& piece:volume.overlaps) {
            const auto leaf=std::find_if(measured.leaves.begin(),measured.leaves.end(),[&](const auto& candidate) {
                return candidate.cell==layout.nativeCell(piece.pressure) && samePoint(candidate.lower,piece.lower) && samePoint(candidate.upper,piece.upper); });
            require(leaf!=measured.leaves.end() && relative(leaf->capacity,piece.volume)<=geometryTolerance,"dual phase overlap mismatch");
            overlap.push_back({leaf->liquid,leaf->gas});
            auto& sum=cellSums[component][piece.pressure]; sum.liquid+=leaf->liquid; sum.gas+=leaf->gas;
        }
        d->overlaps.push_back(std::move(overlap)); d->phases.push_back({measured.liquid,measured.gas});
        const long double ml=liquidDensity*measured.liquid,mg=gasDensity*measured.gas;
        require(std::isfinite(ml) && std::isfinite(mg) && (measured.liquid==0 || ml>0) && (measured.gas==0 || mg>0),
            "material product overflow or underflow");
        const long double mass=ml+mg; const double represented=static_cast<double>(mass);
        require(std::isfinite(mass) && std::isfinite(represented) && represented>0,"unrepresentable momentum mass");
        d->masses.push_back(represented);
        auto& total=d->report.totalPhaseVolumes[component]; total.liquid+=measured.liquid; total.gas+=measured.gas;
        d->report.totalMass[component]+=mass;
        d->report.cutOverlaps+=measured.cutLeaves;
        d->report.maxRelativeLedgerDefect=std::max(d->report.maxRelativeLedgerDefect,measured.maxRelativeLedgerDefect);
        d->report.maxRelativeMomentEstimate=std::max(d->report.maxRelativeMomentEstimate,measured.momentErrorEstimate/measured.capacity);
    }
    for (const auto& sum:cellSums) for (int p=0;p<layout.pressureCount();++p) {
        const auto& reference=d->phase.cellVolumes()[p];
        const auto defect=std::max(std::abs(sum[p].liquid-reference.liquid),std::abs(sum[p].gas-reference.gas))/d->geometry.metrics().cellVolumes()[p];
        require(std::isfinite(defect) && defect<=1e-10L,"dual phase cell partition defect");
        d->report.maxRelativeCellPartition=std::max(d->report.maxRelativeCellPartition,defect);
    }
    for (int component=0;component<2;++component)
        require(std::isfinite(d->report.totalMass[component]) && std::isfinite(d->report.totalPhaseVolumes[component].liquid)
            && std::isfinite(d->report.totalPhaseVolumes[component].gas),"material total overflow");
    const auto assemblyEnd=std::chrono::steady_clock::now();
    d->report.sourcePreparationSeconds=std::chrono::duration<double>(batchStart-sourceStart).count();
    d->report.regionBatchSeconds=std::chrono::duration<double>(assemblyStart-batchStart).count();
    d->report.massAssemblySeconds=std::chrono::duration<double>(assemblyEnd-assemblyStart).count();
    const auto work=region.workProfile();
    d->report.cutPairWorkerSeconds=work.pairWorkerSeconds;
    d->report.cutFallbackWorkerSeconds=work.fallbackWorkerSeconds;
    d->report.cutPairCalls=work.pairCalls;d->report.cutFallbackCalls=work.fallbackCalls;
    d->report.cutPairRegionGaugeWorkerSeconds=work.pairRegionGaugeWorkerSeconds;
    d->report.cutPairCutPreparationWorkerSeconds=work.pairCutPreparationWorkerSeconds;
    d->report.cutPairCoreWorkerSeconds=work.pairCoreWorkerSeconds;
    d->report.cutPairRetryWorkerSeconds=work.pairRetryWorkerSeconds;
    d->report.cutPairAssemblyWorkerSeconds=work.pairAssemblyWorkerSeconds;
    d->report.cutPairCutConstructions=work.pairCutConstructions;
    d->report.cutPairCoreAttempts=work.pairCoreAttempts;
    d->report.cutPairFaceCutRefinements=work.pairFaceCutRefinements;
    d->report.cutPairCoreWork=work.pairCoreWork;
    return HydroMomentumMass(std::move(d));
}


HydroMomentumMass HydroMomentumMass::buildVolumeFractionGeometry(const QuadtreeMesh& mesh,
    HydroMomentumGeometry geometry,HydroPhaseGeometry phase,double liquidDensity,double gasDensity)
{
    const auto start=std::chrono::steady_clock::now();
    geometry.validateMesh(mesh);
    require(phase.model()==HydroPhaseGeometry::Model::GeometricVof,
        "requires a retained volume-fraction geometry owner");
    phase.validateVofFor(mesh,phase.vofIdentity(),phase.timeLayer());
    require(geometry.layout().sameIdentity(phase.layout()),"different sharp mass layout identity");
    require(&geometry.metrics().divergence().matrix()==&phase.metrics().divergence().matrix(),
        "different retained sharp mass D_I");
    require(std::isfinite(liquidDensity)&&liquidDensity>0&&std::isfinite(gasDensity)&&gasDensity>0,
        "invalid sharp material density");
    auto d=std::make_shared<Data>(std::move(geometry),std::move(phase),liquidDensity,gasDensity);
    const auto& layout=d->geometry.layout();
    std::vector<HydroVofBox> queries;
    for(const auto& volume:d->geometry.volumes())for(const auto& overlap:volume.overlaps)
        queries.push_back({layout.nativeCell(overlap.pressure),{overlap.lower,overlap.upper}});
    const auto batchStart=std::chrono::steady_clock::now();
    const auto measures=d->phase.measureVofBoxes(mesh,queries);
    const auto assemblyStart=std::chrono::steady_clock::now();
    require(measures.size()==queries.size(),"incomplete sharp dual phase batch");
    struct PositiveSum {
        long double sum=0,error=0;
        void add(long double x) {
            require(std::isfinite(x)&&x>=0,"invalid positive sharp material contribution");
            const auto y=x-error,t=sum+y;error=(t-sum)-y;sum=t;
            require(std::isfinite(sum),"sharp material sum overflow");
        }
    };
    struct PhaseSum {PositiveSum liquid,gas;void add(HydroPhaseMeasure p){liquid.add(p.liquid);gas.add(p.gas);}};
    std::array<std::vector<PhaseSum>,2> cellSums;
    for(auto& component:cellSums)component.resize(layout.pressureCount());
    std::array<PhaseSum,2> totals;std::array<PositiveSum,2> totalMass;
    std::size_t next=0;
    for(int u=0;u<layout.velocityCount();++u) {
        const auto& volume=d->geometry.volumes()[u];
        const int component=layout.normalAxis(u)==HydroDofLayout::NormalAxis::First?0:1;
        PhaseSum sum;std::vector<HydroPhaseMeasure> overlaps;overlaps.reserve(volume.overlaps.size());
        for(const auto& piece:volume.overlaps) {
            const auto p=measures[next++];
            const auto defect=std::abs(p.liquid+p.gas-piece.volume)/piece.volume;
            require(std::isfinite(defect)&&defect<=geometryTolerance,"sharp overlap capacity closure failed");
            d->report.maxRelativeLedgerDefect=std::max(d->report.maxRelativeLedgerDefect,defect);
            overlaps.push_back(p);sum.add(p);cellSums[component][piece.pressure].add(p);
            d->report.cutOverlaps+=p.liquid>0&&p.gas>0;
        }
        const HydroPhaseMeasure p{sum.liquid.sum,sum.gas.sum};
        require(std::abs(p.liquid+p.gas-volume.volume)/volume.volume<=geometryTolerance,
            "sharp dual capacity closure failed");
        const long double ml=liquidDensity*p.liquid,mg=gasDensity*p.gas;
        require(std::isfinite(ml)&&std::isfinite(mg)&&(p.liquid==0||ml>0)&&(p.gas==0||mg>0),
            "sharp material product overflow or underflow");
        const long double mass=ml+mg;const double represented=static_cast<double>(mass);
        require(std::isfinite(mass)&&std::isfinite(represented)&&represented>0,"unrepresentable sharp momentum mass");
        d->overlaps.push_back(std::move(overlaps));d->phases.push_back(p);d->masses.push_back(represented);
        totals[component].add(p);totalMass[component].add(mass);
    }
    require(next==measures.size(),"unused sharp dual phase measurement");
    const bool localPlic=d->phase.vofPublication().geometryKind()==Interface::GeometricVofGeometryKind::LocalPlic;
    // Local PLIC's whole rectangle and bounded partial-box clips take different
    // double-coordinate arithmetic paths. Keep a fixed capacity-scaled rounding
    // allowance, without remapping any measured volume or dropping a small phase.
    constexpr long double localPlicPartitionRelativeTolerance=1e-12L;
    constexpr long double localPlicPartitionCapacityRoundoff=64.L*std::numeric_limits<double>::epsilon();
    for(const auto& component:cellSums)for(int p=0;p<layout.pressureCount();++p) {
        const auto& reference=d->phase.cellVolumes()[p];
        const auto relativePhase=[](long double measured,long double original) {
            return original==0?(measured==0?0.L:std::numeric_limits<long double>::infinity()):std::abs(measured-original)/original;
        };
        const auto liquid=component[p].liquid.sum,gas=component[p].gas.sum;
        const auto defect=std::max(relativePhase(liquid,reference.liquid),relativePhase(gas,reference.gas));
        const auto cellCapacity=d->geometry.metrics().cellVolumes()[p];
        const auto roundoff=localPlicPartitionCapacityRoundoff*cellCapacity;
        const auto admissibleLocalPhase=[&](long double measured,long double original) {
            if(!std::isfinite(measured)||!std::isfinite(original)||measured<0||original<0)return false;
            if(original==0)return measured==0; // Exact pure-phase zero stays exact.
            return measured>0&&std::abs(measured-original)<=localPlicPartitionRelativeTolerance*original+roundoff;
        };
        const bool accepted=localPlic
            ?std::isfinite(cellCapacity)&&cellCapacity>0&&admissibleLocalPhase(liquid,reference.liquid)&&admissibleLocalPhase(gas,reference.gas)
            :std::isfinite(defect)&&defect<=1e-12L;
        const auto capacityDefect=std::max(std::abs(liquid-reference.liquid),std::abs(gas-reference.gas))/cellCapacity;
        if(!accepted) {
            std::ostringstream context;context.precision(std::numeric_limits<long double>::max_digits10);
            context<<"[HydroMomentumGeometry] actual phase cell partition defect model="<<"VOF"
                <<" owner="<<d->phase.vofIdentity().value()
                <<" nativeCell="<<layout.nativeCell(p)<<" component="<<(&component==&cellSums[0]?0:1)
                <<" actualLiquid="<<reference.liquid<<" partitionLiquid="<<liquid
                <<" actualGas="<<reference.gas<<" partitionGas="<<gas<<" relative="<<defect;
            if(localPlic)context<<" localPlicCapacityScaled="<<capacityDefect<<" cellCapacity="<<cellCapacity
                <<" capacityRoundoff="<<roundoff<<" relativeTolerance="<<localPlicPartitionRelativeTolerance;
            throw std::invalid_argument(context.str());
        }
        d->report.maxRelativeCellPartition=std::max(d->report.maxRelativeCellPartition,defect);
        if(localPlic) {
            const auto budgetFraction=std::max(std::abs(liquid-reference.liquid)/(localPlicPartitionRelativeTolerance*reference.liquid+roundoff),
                std::abs(gas-reference.gas)/(localPlicPartitionRelativeTolerance*reference.gas+roundoff));
            d->report.maxLocalPlicCapacityScaledCellPartition=std::max(d->report.maxLocalPlicCapacityScaledCellPartition,capacityDefect);
            d->report.maxLocalPlicPartitionBudgetFraction=std::max(d->report.maxLocalPlicPartitionBudgetFraction,budgetFraction);
        }
    }
    for(int component=0;component<2;++component) {
        d->report.totalPhaseVolumes[component]={totals[component].liquid.sum,totals[component].gas.sum};
        d->report.totalMass[component]=totalMass[component].sum;
    }
    const auto finish=std::chrono::steady_clock::now();
    d->report.sourcePreparationSeconds=std::chrono::duration<double>(batchStart-start).count();
    d->report.regionBatchSeconds=std::chrono::duration<double>(assemblyStart-batchStart).count();
    d->report.massAssemblySeconds=std::chrono::duration<double>(finish-assemblyStart).count();
    return HydroMomentumMass(std::move(d));
}

HydroMomentumMass HydroMomentumMass::buildCellAverageProjection(const QuadtreeMesh& mesh,
    HydroMomentumGeometry geometry,HydroPhaseGeometry phase,double liquidDensity,double gasDensity)
{
    const bool vof=phase.model()==HydroPhaseGeometry::Model::GeometricVof;
    require(vof||phase.model()==HydroPhaseGeometry::Model::FirstJet,
        "cell-average projection requires a supported retained phase owner");
    const auto sourceStart=std::chrono::steady_clock::now();
    geometry.validateMesh(mesh);
    if(vof){
        phase.validateVofFor(mesh,phase.vofIdentity(),phase.timeLayer());
        require(phase.vofPublication().geometryKind()==Interface::GeometricVofGeometryKind::LocalPlic,
            "VOF cell-average projection requires independent LocalPlic phase inventories");
    }else phase.validateFor(mesh,phase.geometryIdentity(),phase.timeLayer());
    require(geometry.layout().sameIdentity(phase.layout()),"different projected mass layout identity");
    require(&geometry.metrics().divergence().matrix()==&phase.metrics().divergence().matrix(),
        "different projected mass D_I");
    require(std::isfinite(liquidDensity)&&liquidDensity>0&&std::isfinite(gasDensity)&&gasDensity>0,
        "invalid projected material density");
    LAMP::Geometry::Detail::ExactExpansion::validateEnvironment();
    auto d=std::make_shared<Data>(std::move(geometry),std::move(phase),liquidDensity,gasDensity);
    d->representation=DensityRepresentation::CellAverageProjection;
    const auto assemblyStart=std::chrono::steady_clock::now();const auto& layout=d->geometry.layout();
    std::array<std::vector<HydroPhaseMeasure>,2> cellSums;
    for(auto& sum:cellSums)sum.resize(layout.pressureCount());
    // Finite-CDF replacement and validation consume the final endpoint owner in place.
    d->overlaps.resize(layout.velocityCount());
    if(vof){
        bool hasFinite=false;
        for(auto& flags:d->finiteP0Partitions)flags.resize(layout.pressureCount());
        for(int u=0;u<layout.velocityCount();++u){const int component=layout.normalAxis(u)==HydroDofLayout::NormalAxis::First?0:1;
            const auto& pieces=d->geometry.volumes()[u].overlaps;d->overlaps[u].resize(pieces.size());
            for(std::size_t k=0;k<pieces.size();++k){const auto& piece=pieces[k];
                const auto& stock=d->phase.vofPublication().plicCells()[layout.nativeCell(piece.pressure)].inventory;
                require(std::isfinite(piece.cellFraction)&&piece.cellFraction>0,"invalid VOF projected mass overlap weight");
                d->overlaps[u][k]={piece.cellFraction*stock[0],piece.cellFraction*stock[1]};
                const std::array<long double,2> raw{{d->overlaps[u][k].liquid,d->overlaps[u][k].gas}};
                for(int phase=0;phase<2;++phase){require(std::isfinite(raw[phase])&&raw[phase]>=0,"invalid ordinary VOF projected phase moment");
                    if(stock[phase]>0&&(raw[phase]==0||std::fpclassify(raw[phase])==FP_SUBNORMAL)){
                        d->finiteP0Partitions[component][piece.pressure][phase]=true;hasFinite=true;}}
            }
        }
        // Ordinary P0 moments never consume the finite-CDF reference table.
        // Build it only after an actual phase has reached the finite range.
        if(hasFinite){
        using Ref=std::pair<int,std::size_t>;
        std::array<std::vector<std::vector<Ref>>,2> groups;
        for(auto& group:groups)group.resize(layout.pressureCount());
        for(int u=0;u<layout.velocityCount();++u){const int component=layout.normalAxis(u)==HydroDofLayout::NormalAxis::First?0:1;
            const auto& pieces=d->geometry.volumes()[u].overlaps;
            for(std::size_t k=0;k<pieces.size();++k)groups[component][pieces[k].pressure].push_back({u,k});
        }
        for(int component=0;component<2;++component)for(int p=0;p<layout.pressureCount();++p){
            const auto& flags=d->finiteP0Partitions[component][p];if(!flags[0]&&!flags[1])continue;
            auto& refs=groups[component][p];require(!refs.empty()&&refs.size()<=4,"invalid finite P0 component partition size");
            const auto key=[&](const Ref& r){const auto& b=d->geometry.volumes()[r.first].overlaps[r.second];
                return component==0?std::make_tuple(b.lower.x,b.upper.x,b.lower.y,b.upper.y):std::make_tuple(b.lower.y,b.upper.y,b.lower.x,b.upper.x);};
            std::sort(refs.begin(),refs.end(),[&](const Ref& a,const Ref& b){return key(a)<key(b);});
            for(std::size_t k=1;k<refs.size();++k)require(key(refs[k-1])!=key(refs[k]),"duplicate finite P0 rectangle");
            const auto& source=d->phase.vofPublication().plicCells()[layout.nativeCell(p)];
            Detail::HydroFiniteContraction zero,whole,prefix;const bool rz=layout.domain().axisymmetric();
            addP0Rectangle(whole,source.rectangle.lower,source.rectangle.upper,rz);require(whole.compare(zero)>0,"nonpositive exact finite P0 parent volume");
            std::vector<Detail::HydroFiniteContraction> endpoints;endpoints.reserve(refs.size());
            for(const auto& ref:refs){const auto& piece=d->geometry.volumes()[ref.first].overlaps[ref.second];
                Detail::HydroFiniteContraction measure;addP0Rectangle(measure,piece.lower,piece.upper,rz);
                require(measure.compare(zero)>0,"nonpositive exact finite P0 piece volume");prefix.add(measure);
                require(prefix.compare(whole)<=0,"finite P0 geometry prefix exceeds actual parent");endpoints.push_back(prefix);}
            // Exact full coverage is proved BEFORE assigning any phase piece.
            // The full endpoint is selected by this geometric identity, never
            // by array position or a posthoc residual/normalization.
            require(prefix.compare(whole)==0,"finite P0 geometry does not exactly partition actual parent");
            for(int phase=0;phase<2;++phase)if(flags[phase]){const auto stock=source.inventory[phase];long double previous=0;
                require(std::isfinite(stock)&&stock>0,"invalid finite P0 source stock");++d->report.finiteP0PhasePartitions;
                for(std::size_t k=0;k<refs.size();++k){const auto& endpoint=endpoints[k];
                    const auto current=endpoint.compare(zero)==0?0.L:endpoint.compare(whole)==0?stock:endpoint.scaledQuotient(stock,whole);
                    require(std::isfinite(current)&&current>=previous&&current<=stock,"invalid finite P0 cumulative phase endpoint");
                    Detail::HydroFiniteContraction difference;difference.add(current);difference.add(-previous);const auto value=difference.rounded();
                    require(std::isfinite(value)&&value>=0&&value<=stock,"invalid finite P0 phase piece");
                    const auto& ref=refs[k];auto& target=d->overlaps[ref.first][ref.second];(phase==0?target.liquid:target.gas)=value;
                    ++d->report.finiteP0PhasePieces;d->report.finiteP0ZeroPhasePieces+=value==0;previous=current;
                }
                require(previous==stock,"finite P0 complete phase endpoint differs from source stock");
            }
        }
        }
    }
    // A derived LocalPLIC publication retains the exact, producer-certified
    // defect kappa = Q_L + Q_G - C. P0 material moments must carry f*kappa;
    // neither independent stock nor physical M is normalized to geometry.
    // Only a failed original fast gate pays for this exact source check. The
    // original 512-epsilon budget still bounds geometry and rounding residuals.
    const auto certifiedCapacity=[&](int component,const HydroMomentumVolume& volume,
        std::size_t first,std::size_t last,const HydroPhaseMeasure& observed,long double full){
        if(!vof||!std::isfinite(observed.liquid)||!std::isfinite(observed.gas)
            ||observed.liquid<0||observed.gas<0||!std::isfinite(full)||full<=0)return false;
        using LAMP::Geometry::Detail::FiniteVolume::Dyadic;
        const auto& publication=d->phase.vofPublication();
        Dyadic residual(observed.liquid);residual.add(observed.gas);residual.add(full,-1);
        for(std::size_t k=first;k<last;++k){const auto& piece=volume.overlaps[k];
            const auto& finite=d->finiteP0Partitions[component][piece.pressure];
            // The exceptional finite-CDF representation has different piece
            // moments. Its existing gates stay unchanged until separately proved.
            if(finite[0]||finite[1])return false;
            const int native=layout.nativeCell(piece.pressure);
            const auto& stock=publication.plicCells()[native].inventory;
            const long double parent=publication.cells()[native].capacity;
            Dyadic actualDefect(stock[0]);actualDefect.add(stock[1]);actualDefect.add(parent,-1);
            require(actualDefect.compare(publication.capacityDefects()[native])==0,
                "projected material differs from retained source capacity certificate");
            // Exact f*kappa expanded from that verified association, with no
            // rounding of kappa or f*kappa and no multiplication of a receipt.
            residual.add(Dyadic::product({piece.cellFraction,stock[0]}),-1);
            residual.add(Dyadic::product({piece.cellFraction,stock[1]}),-1);
            residual.add(Dyadic::product({piece.cellFraction,parent}));
        }
        const long double scale=std::max(std::abs(observed.liquid+observed.gas),std::abs(full));
        return std::isfinite(scale)&&residual.absolute().compare(Dyadic::product({geometryTolerance,scale}))<=0;
    };
    for(int u=0;u<layout.velocityCount();++u){
        const auto& volume=d->geometry.volumes()[u];
        const int component=layout.normalAxis(u)==HydroDofLayout::NormalAxis::First?0:1;
        HydroPhaseMeasure represented;auto& overlaps=d->overlaps[u];
        if(!vof)overlaps.resize(volume.overlaps.size());
        for(std::size_t k=0;k<volume.overlaps.size();++k){const auto& piece=volume.overlaps[k];
            HydroPhaseMeasure value;
            if(vof){
                // Explicit P0 MATERIAL moments from both retained stocks;
                // these are not measurements of the PLIC cut in this overlap.
                const auto& stock=d->phase.vofPublication().plicCells()[layout.nativeCell(piece.pressure)].inventory;
                require(std::isfinite(piece.cellFraction)&&piece.cellFraction>0,
                    "invalid VOF projected mass overlap weight");
                value=d->overlaps[u][k];const auto& finite=d->finiteP0Partitions[component][piece.pressure];
                require((stock[0]==0||value.liquid>0||finite[0])&&(stock[1]==0||value.gas>0||finite[1]),
                    "unrepresentable independent VOF projected phase moment");
            }else{
                const auto& cell=d->phase.cellVolumes()[piece.pressure];
                require(std::isfinite(piece.cellFraction)&&piece.cellFraction>0,
                    "invalid projected mass overlap weight");
                value={piece.cellFraction*cell.liquid,piece.cellFraction*cell.gas};
            }
            if(!(std::isfinite(value.liquid)&&std::isfinite(value.gas)&&value.liquid>=0&&value.gas>=0
                &&(relative(value.liquid+value.gas,piece.volume)<=geometryTolerance
                    ||certifiedCapacity(component,volume,k,k+1,value,piece.volume)))){
                const auto clock=d->phase.timeLayer();
                std::ostringstream why;why.precision(std::numeric_limits<long double>::max_digits10);
                why<<"[HydroMomentumGeometry] invalid projected phase moment capacity; vof="<<vof
                    <<" u="<<u<<" overlap="<<k<<" component="<<component
                    <<" pressure="<<piece.pressure<<" native="<<layout.nativeCell(piece.pressure)
                    <<" clock=("<<clock.step<<','<<clock.begin<<','<<clock.end<<','<<static_cast<int>(clock.position)<<')'
                    <<" represented=("<<value.liquid<<','<<value.gas<<") volume="<<piece.volume
                    <<" relative="<<relative(value.liquid+value.gas,piece.volume)<<" tolerance="<<geometryTolerance
                    <<" fraction="<<piece.cellFraction
                    <<" rectangle=[("<<piece.lower.x<<','<<piece.lower.y<<"),("<<piece.upper.x<<','<<piece.upper.y<<")]";
                if(vof){const auto& source=d->phase.vofPublication().plicCells()[layout.nativeCell(piece.pressure)];
                    const auto& publication=d->phase.vofPublication();const int native=layout.nativeCell(piece.pressure);
                    const auto& finite=d->finiteP0Partitions[component][piece.pressure];
                    why<<" vof_id="<<d->phase.vofIdentity().value()
                        <<" stock=("<<source.inventory[0]<<','<<source.inventory[1]<<") measured=("<<source.measured[0]<<','<<source.measured[1]<<')'
                        <<" native_capacity="<<publication.cells()[native].capacity
                        <<" hydro_capacity="<<d->geometry.metrics().cellVolumes()[piece.pressure]
                        <<" certified_defect="<<publication.capacityDefects()[native].rounded()
                        <<" finite_flags=("<<finite[0]<<','<<finite[1]<<')'
                        <<" parent=[("<<source.rectangle.lower.x<<','<<source.rectangle.lower.y<<"),("<<source.rectangle.upper.x<<','<<source.rectangle.upper.y<<")]";
                }
                throw std::invalid_argument(why.str());
            }
            if(!vof)overlaps[k]=value;
            represented.liquid+=value.liquid;represented.gas+=value.gas;
            auto& sum=cellSums[component][piece.pressure];sum.liquid+=value.liquid;sum.gas+=value.gas;
            if(value.liquid>0&&value.gas>0)++d->report.cutOverlaps;
        }
        require(relative(represented.liquid+represented.gas,volume.volume)<=geometryTolerance
            ||certifiedCapacity(component,volume,0,volume.overlaps.size(),represented,volume.volume),
            "projected mass dual capacity mismatch");
        d->phases.push_back(represented);
        const long double ml=liquidDensity*represented.liquid,mg=gasDensity*represented.gas;
        require(std::isfinite(ml)&&std::isfinite(mg)&&(represented.liquid==0||ml>0)
            &&(represented.gas==0||mg>0),"projected material product overflow or underflow");
        const long double mass=ml+mg;const double stored=static_cast<double>(mass);
        require(std::isfinite(mass)&&std::isfinite(stored)&&stored>0,"unrepresentable projected momentum mass");
        d->masses.push_back(stored);
        auto& total=d->report.totalPhaseVolumes[component];
        total.liquid+=represented.liquid;total.gas+=represented.gas;d->report.totalMass[component]+=mass;
    }
    for(int component=0;component<2;++component)for(int p=0;p<layout.pressureCount();++p){
        const auto& sum=cellSums[component];HydroPhaseMeasure reference;
        if(vof){
            const auto& stock=d->phase.vofPublication().plicCells()[layout.nativeCell(p)].inventory;
            reference={stock[0],stock[1]};
        }else reference=d->phase.cellVolumes()[p];
        const auto defect=std::max(std::abs(sum[p].liquid-reference.liquid),
            std::abs(sum[p].gas-reference.gas))/d->geometry.metrics().cellVolumes()[p];
        require(std::isfinite(defect)&&defect<=1e-10L,"projected dual phase cell partition defect");
        d->report.maxRelativeCellPartition=std::max(d->report.maxRelativeCellPartition,defect);
        if(vof){const std::array<long double,2> actual{{sum[p].liquid,sum[p].gas}},parent{{reference.liquid,reference.gas}};
            for(int phase=0;phase<2;++phase){const auto own=parent[phase]==0?(actual[phase]==0?0.L:std::numeric_limits<long double>::infinity()):std::abs(actual[phase]-parent[phase])/parent[phase];
                require(std::isfinite(own)&&own<=1e-12L,"independent VOF projected own-phase addback failed");
                d->report.maxRelativeP0OwnPhaseCellPartition[component][phase]=std::max(d->report.maxRelativeP0OwnPhaseCellPartition[component][phase],own);}
        }
    }
    for(int component=0;component<2;++component)
        require(std::isfinite(d->report.totalMass[component])
            &&std::isfinite(d->report.totalPhaseVolumes[component].liquid)
            &&std::isfinite(d->report.totalPhaseVolumes[component].gas),"projected material total overflow");
    d->report.sourcePreparationSeconds=std::chrono::duration<double>(assemblyStart-sourceStart).count();
    d->report.massAssemblySeconds=std::chrono::duration<double>(
        std::chrono::steady_clock::now()-assemblyStart).count();
    return HydroMomentumMass(std::move(d));
}

bool HydroMomentumMass::usesFiniteP0Partition(int pressure,int component,int phase) const
{
    require(model()==HydroPhaseGeometry::Model::GeometricVof&&representation()==DensityRepresentation::CellAverageProjection,
        "finite P0 partition flag requires genuine VOF projected material");
    require(pressure>=0&&pressure<layout().pressureCount()&&component>=0&&component<2&&phase>=0&&phase<2,
        "invalid finite P0 partition flag index");return data_->finiteP0Partitions[component][pressure][phase];
}
HydroMomentumMass::DensityRepresentation HydroMomentumMass::representation() const noexcept
{return data_->representation;}
HydroPhaseGeometry::Model HydroMomentumMass::model() const noexcept {return phaseGeometry().model();}
Interface::GeometricVofIdentity HydroMomentumMass::vofIdentity() const {return phaseGeometry().vofIdentity();}
const HydroMomentumGeometry& HydroMomentumMass::geometry() const noexcept { return data_->geometry; }
const HydroPhaseGeometry& HydroMomentumMass::phaseGeometry() const noexcept { return data_->phase; }
const HydroDofLayout& HydroMomentumMass::layout() const noexcept { return geometry().layout(); }
Interface::HermiteGalsIdentity HydroMomentumMass::geometryIdentity() const { return phaseGeometry().geometryIdentity(); }
LAMP::TimeLayer HydroMomentumMass::timeLayer() const
{return phaseGeometry().timeLayer();}
double HydroMomentumMass::liquidDensity() const noexcept { return data_->liquid; }
double HydroMomentumMass::gasDensity() const noexcept { return data_->gas; }
const std::vector<HydroPhaseMeasure>& HydroMomentumMass::phaseVolumes() const
{require(representation()==DensityRepresentation::SharpGeometry,"phaseVolumes requires actual sharp geometry");return data_->phases;}
const std::vector<std::vector<HydroPhaseMeasure>>& HydroMomentumMass::overlapPhaseVolumes() const
{require(representation()==DensityRepresentation::SharpGeometry,"overlapPhaseVolumes requires actual sharp geometry");return data_->overlaps;}
const std::vector<HydroPhaseMeasure>& HydroMomentumMass::representedPhaseVolumes() const noexcept
{return data_->phases;}
const std::vector<std::vector<HydroPhaseMeasure>>& HydroMomentumMass::representedOverlapPhaseVolumes() const noexcept
{return data_->overlaps;}
const std::vector<double>& HydroMomentumMass::masses() const noexcept { return data_->masses; }
const HydroMomentumMassReport& HydroMomentumMass::report() const noexcept { return data_->report; }
void HydroMomentumMass::validateFor(const QuadtreeMesh& mesh,Interface::HermiteGalsIdentity identity,LAMP::TimeLayer time) const
{ geometry().validateMesh(mesh); phaseGeometry().validateFor(mesh,identity,time); }
void HydroMomentumMass::validateVofFor(const QuadtreeMesh& mesh,Interface::GeometricVofIdentity identity,
    Interface::GeometricVofTimeLayer time) const
{
    // Both mass factories certify this same geometry/phase mesh binding.
    // The phase entry checks the live mesh plus the actual owner/full clock.
    phaseGeometry().validateVofFor(mesh,identity,time);
}

namespace
{
constexpr std::array<double,8> kineticNodes{{-.96028985649753623168,-.79666647741362673959,
    -.52553240991632898582,-.18343464249564980494,.18343464249564980494,
    .52553240991632898582,.79666647741362673959,.96028985649753623168}};
constexpr std::array<double,8> kineticWeights{{.10122853629037625915,.22238103445337447054,
    .31370664587788728734,.36268378337836198297,.36268378337836198297,
    .31370664587788728734,.22238103445337447054,.10122853629037625915}};
constexpr std::array<double,16> kineticNodes16{{-.98940093499164993260,-.94457502307323257608,
    -.86563120238783174388,-.75540440835500303390,-.61787624440264374845,-.45801677765722738634,
    -.28160355077925891323,-.09501250983763744019,.09501250983763744019,.28160355077925891323,
    .45801677765722738634,.61787624440264374845,.75540440835500303390,.86563120238783174388,
    .94457502307323257608,.98940093499164993260}};
constexpr std::array<double,16> kineticWeights16{{.02715245941175409485,.06225352393864789286,
    .09515851168249278481,.12462897125553387205,.14959598881657673208,.16915651939500253819,
    .18260341504492358887,.18945061045506849629,.18945061045506849629,.18260341504492358887,
    .16915651939500253819,.14959598881657673208,.12462897125553387205,.09515851168249278481,
    .06225352393864789286,.02715245941175409485}};
template<std::size_t N>
std::pair<std::array<double,N>,std::array<double,N>> kineticGaussian()
{
    std::pair<std::array<double,N>,std::array<double,N>> rule;
    const long double pi=std::acos(-1.L);
    for(std::size_t i=0;i<(N+1)/2;++i){
        long double x=std::cos(pi*(i+.75L)/(N+.5L)),derivative=0;
        bool converged=false;
        for(int iteration=0;iteration<64;++iteration){long double p0=1,p1=x;
            for(std::size_t k=2;k<=N;++k){const long double p=((2*k-1)*x*p1-(k-1)*p0)/k;p0=p1;p1=p;}
            derivative=N*(x*p1-p0)/(x*x-1);const long double delta=p1/derivative;
            x-=delta;if(std::abs(delta)<=4*std::numeric_limits<long double>::epsilon()){converged=true;break;}}
        require(converged,"kinetic Gauss node failed to converge");
        // Re-evaluate the derivative at the retained node.
        long double p0=1,p1=x;
        for(std::size_t k=2;k<=N;++k){const long double p=((2*k-1)*x*p1-(k-1)*p0)/k;p0=p1;p1=p;}
        derivative=N*(x*p1-p0)/(x*x-1);
        const double weight=static_cast<double>(2/((1-x*x)*derivative*derivative));
        rule.first[i]=-static_cast<double>(x);rule.first[N-1-i]=static_cast<double>(x);
        rule.second[i]=rule.second[N-1-i]=weight;
    }
    return rule;
}
const auto kineticRule20=kineticGaussian<20>();
const auto kineticRule32=kineticGaussian<32>();
using KineticMoments=std::vector<long double>; // P8 baseline; up to P34 modal RZ.
struct KineticRule
{
    explicit KineticRule(int degree=8):moments(2*(degree+1)*(degree+1)){}
    std::vector<HydroKineticHodgePoint> points;
    KineticMoments moments{};
    long double error=0;
};
void appendKineticRule(KineticRule& a,KineticRule&& b)
{
    require(a.moments.size()==b.moments.size(),"different kinetic moment spaces");
    for(std::size_t k=0;k<a.moments.size();++k)a.moments[k]+=b.moments[k];
    a.error+=b.error;
    a.points.insert(a.points.end(),b.points.begin(),b.points.end());
}
// Cubature of one actual dual/native-cell intersection. Roots are read-only
// ranges of the original Q4; no scalar fitting or rounded root subtraction.
class KineticBoxRule
{
public:
    KineticBoxRule(std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> scalar,
        int dual,int pressure,int native,Vector2D lower,Vector2D upper,bool rz,int momentDegree=8)
        :scalar_(std::move(scalar)),dual_(dual),pressure_(pressure),native_(native),
         lower_(lower),upper_(upper),radial_(rz),momentDegree_(momentDegree)
    {
        require(momentDegree_>=8&&momentDegree_<=34,"invalid kinetic moment degree");
        int containing=-1;
        for(int cell:scalar_->nativeCellCandidates(lower_,upper_)){
            const auto bounds=scalar_->nativeCellBounds(cell);
            if(bounds[0].x<=lower_.x && upper_.x<=bounds[1].x &&
                bounds[0].y<=lower_.y && upper_.y<=bounds[1].y){
                require(containing<0,"kinetic box has ambiguous Q4 source cell");containing=cell;
            }
        }
        require(containing>=0,"kinetic box crosses Q4 source cells");native_=containing;
        const auto nativeBounds=scalar_->nativeCellBounds(native_);
        const Vector2D nativeWidth=nativeBounds[1]-nativeBounds[0];
        const Vector2D localLower{(lower_.x-nativeBounds[0].x)/nativeWidth.x,
            (lower_.y-nativeBounds[0].y)/nativeWidth.y};
        const Vector2D localUpper{(upper_.x-nativeBounds[0].x)/nativeWidth.x,
            (upper_.y-nativeBounds[0].y)/nativeWidth.y};
        if(scalar_->nativeQuarticBoxSign(native_,lower_,upper_)==0){
            const bool yResolved=scalar_->nativeZeroSetPartialResolved(native_,1,localLower,localUpper);
            const bool xResolved=scalar_->nativeZeroSetPartialResolved(native_,0,localLower,localUpper);
            require(yResolved||xResolved,"kinetic cut has no certified regular slice direction");
            // If phi_y is signed, y is a regular graph over x and x is outer.
            // At a vertical tangent switch the roles, using the same Q4 source.
            // When both are certified, choose the stronger physical derivative
            // at the box center to avoid a badly conditioned graph near a turn.
            if(yResolved&&xResolved){
                const auto jet=scalar_->scalarJetInNativeCell(native_,(lower_+upper_)*.5);
                outerAxis_=std::abs(jet.dy)>=std::abs(jet.dx)?0:1;
            }else outerAxis_=yResolved?0:1;
        }
        const int axis=outerAxis_;
        const auto bounds=scalar_->nativeCellBounds(native_);
        for(double x:{coord(bounds[0],axis),coord(bounds[1],axis)})
            if(x>coord(lower_,axis) && x<coord(upper_,axis))nativeKnots_.push_back(x);
        std::sort(nativeKnots_.begin(),nativeKnots_.end());
        nativeKnots_.erase(std::unique(nativeKnots_.begin(),nativeKnots_.end()),nativeKnots_.end());
    }
    int outerAxis() const noexcept{return outerAxis_;}
    long double outerCapacity(double lo,double hi) const
    {
        return outerAxis_==0?capacity({lo,lower_.y},{hi,upper_.y},radial_)
            :capacity({lower_.x,lo},{upper_.x,hi},radial_);
    }
    KineticRule integrate(double lo,double hi,long double budget)
    {
        auto result=integrateRange(lo,hi,budget,0);
        require(result.error<=budget,"kinetic cubature box exceeded its total error budget");
        return result;
    }
private:
    KineticRule integrateRange(double lo,double hi,long double budget,int depth)
    {
        // Split at certified interface events on the inner-axis box sides.
        // Without these events, phase moments can jump in the outer integral
        // when the interface enters or leaves the integration rectangle.
        const int innerAxis=1-outerAxis_;
        const double innerLo=coord(lower_,innerAxis),innerHi=coord(upper_,innerAxis);
        std::vector<double> cellCuts{lo};
        for(double x:nativeKnots_)if(x>lo && x<hi)cellCuts.push_back(x);
        cellCuts.push_back(hi);std::sort(cellCuts.begin(),cellCuts.end());
        cellCuts.erase(std::unique(cellCuts.begin(),cellCuts.end()),cellCuts.end());
        std::vector<double> cuts=cellCuts;
        std::vector<std::pair<double,double>> eventBands;
        for(double inner:{innerLo,innerHi}){
            for(std::size_t k=1;k<cellCuts.size();++k){
                const double a=cellCuts[k-1],b=cellCuts[k];
                LevelSetMappedSegment edge(scalar_,onPlane(innerAxis,inner,a),onPlane(innerAxis,inner,b));
                require(edge.usesSource(*scalar_),"foreign kinetic interface event source");
                for(std::size_t r=0;r<edge.roots().size();++r){
                    // Endpoint incidence is exact retained geometry. Roots
                    // are interpreted only inside their single native-cell
                    // chart, never as parameters of a multi-cell edge. The
                    // chart orientation is retained by the root but is not
                    // necessarily the input segment orientation; use the
                    // authoritative Q4 box sign to select every candidate
                    // enclosure that can contain the physical edge root.
                    if(edge.rootAtInputEndpoint(r,false)||edge.rootAtInputEndpoint(r,true))continue;
                    const auto& root=edge.roots()[r];
                    const long double width=static_cast<long double>(b)-a;
                    int mapped=0;
                    for(bool reverseChart:{false,true}){
                        const long double tlo=reverseChart?1.L-root.upper:root.lower;
                        const long double thi=reverseChart?1.L-root.lower:root.upper;
                        const long double lower=static_cast<long double>(a)+width*tlo;
                        const long double upper=static_cast<long double>(a)+width*thi;
                        double xlo=static_cast<double>(lower),xhi=static_cast<double>(upper);
                        if(static_cast<long double>(xlo)>lower)xlo=std::nextafter(xlo,-std::numeric_limits<double>::infinity());
                        if(static_cast<long double>(xhi)<upper)xhi=std::nextafter(xhi,std::numeric_limits<double>::infinity());
                        xlo=std::max(a,xlo);xhi=std::min(b,xhi);
                        require(std::isfinite(xlo)&&std::isfinite(xhi)&&xlo<=xhi,"invalid kinetic root enclosure");
                        const auto qlo=onPlane(innerAxis,inner,xlo),qhi=onPlane(innerAxis,inner,xhi);
                        if(scalar_->nativeQuarticBoxSign(native_,qlo,qhi)!=0)continue;
                        ++mapped;
                        if(xlo<xhi){eventBands.emplace_back(xlo,xhi);cuts.push_back(xlo);cuts.push_back(xhi);}
                        else{
                            const double x=static_cast<double>(lower+(upper-lower)/2);
                            if(x>a&&x<b)cuts.push_back(x);
                        }
                    }
                    require(mapped>0,"kinetic root chart cannot be mapped to a certified physical edge enclosure");
                }
            }
        }
        std::sort(eventBands.begin(),eventBands.end());
        std::vector<std::pair<double,double>> mergedBands;
        for(const auto& band:eventBands){
            if(mergedBands.empty() || band.first>mergedBands.back().second)mergedBands.push_back(band);
            else mergedBands.back().second=std::max(mergedBands.back().second,band.second);
        }
        eventBands=std::move(mergedBands);
        std::sort(cuts.begin(),cuts.end());cuts.erase(std::unique(cuts.begin(),cuts.end()),cuts.end());
        if(cuts.size()>2 || !eventBands.empty()){
            KineticRule split(momentDegree_);
            const long double subBudget=budget/static_cast<long double>(cuts.size()-1);
            for(std::size_t i=1;i<cuts.size();++i){
                require(cuts[i]>cuts[i-1],"unresolved coincident kinetic interface events");
                const auto band=std::find_if(eventBands.begin(),eventBands.end(),[&](const auto& e){
                    return e.first<=cuts[i-1] && cuts[i]<=e.second;
                });
                if(band!=eventBands.end()){
                    KineticRule bounded(momentDegree_);bounded.error=outerCapacity(cuts[i-1],cuts[i]);
                    appendKineticRule(split,std::move(bounded));
                }else appendKineticRule(split,integrateRange(cuts[i-1],cuts[i],subBudget,depth));
            }
            return split;
        }
        Vector2D subLower=lower_,subUpper=upper_;
        if(outerAxis_==0){subLower.x=lo;subUpper.x=hi;}
        else{subLower.y=lo;subUpper.y=hi;}
        if(const int sign=scalar_->nativeQuarticBoxSign(native_,subLower,subUpper);sign!=0){
            // A pure-phase tensor rule exactly integrates the configured
            // moments: GL8 for P8, GL20 through P34 plus the RZ Jacobian.
            // Comparing two
            // rules here only measures coordinate roundoff and can drive
            // futile bisection below the geometry's representable scale.
            auto resolved=outer(lo,hi,false);
            resolved.error+=256*std::numeric_limits<double>::epsilon()*outerCapacity(lo,hi);
            return resolved;
        }
        if(depth>=40){
            const auto unresolved=outerCapacity(lo,hi);
            KineticRule bounded(momentDegree_);bounded.error=unresolved;return bounded;
        }
        const double middle=lo+(hi-lo)/2;
        if(!(middle>lo && middle<hi)){
            // At adjacent binary64 abscissae, no quadrature point can be
            // represented inside the interval. Its exact physical capacity
            // bounds every normalized phase moment and is added to the box's
            // final absolute error ledger.
            const auto unresolved=outerCapacity(lo,hi);
            KineticRule bounded(momentDegree_);bounded.error=unresolved;return bounded;
        }
        const auto coarse=outer(lo,hi,false);
        auto fine=outer(lo,hi,true);
        long double error=0;
        for(int k=0;k<2*(momentDegree_+1)*(momentDegree_+1);++k)
            error=std::max(error,std::abs(fine.moments[k]-coarse.moments[k]));
        error+=coarse.error+fine.error;
        require(std::isfinite(error),"nonfinite kinetic cubature estimate");
        if(error<=budget){fine.error=error;return fine;}
        // A Q4 interface may cross the outer integration interval, making
        // phase-separated moments only Lipschitz at a root. Keep the original
        // error budget and allow enough bisections to resolve that root in
        // binary64; the independent query cap still bounds total work.
        if(depth>=40)throw std::runtime_error("[HydroMomentumGeometry] kinetic cubature depth exceeded at x=["+
            std::to_string(lo)+","+std::to_string(hi)+"], box=["+std::to_string(lower_.x)+","+
            std::to_string(upper_.x)+"]x["+std::to_string(lower_.y)+","+std::to_string(upper_.y)+
            "], residual="+std::to_string(static_cast<double>(error))+", budget="+
            std::to_string(static_cast<double>(budget)));
        auto result=integrateRange(lo,middle,budget/2,depth+1);
        appendKineticRule(result,integrateRange(middle,hi,budget/2,depth+1));
        return result;
    }
    KineticRule outer(double lo,double hi,bool highOrder)
    {
        if(momentDegree_>8)return highOrder?outerRule(lo,hi,kineticRule32.first,kineticRule32.second):
            outerRule(lo,hi,kineticRule20.first,kineticRule20.second);
        return highOrder?outerRule(lo,hi,kineticNodes16,kineticWeights16):outerRule(lo,hi,kineticNodes,kineticWeights);
    }
    template<std::size_t N>
    KineticRule outerRule(double lo,double hi,const std::array<double,N>& nodes,const std::array<double,N>& weights)
    {
        KineticRule rule(momentDegree_);
        const auto bounded=[&](){KineticRule result(momentDegree_);result.error=outerCapacity(lo,hi);return result;};
        const int innerAxis=1-outerAxis_;
        const double innerLo=coord(lower_,innerAxis),innerHi=coord(upper_,innerAxis);
        const double* innerNodes=momentDegree_>8?kineticRule20.first.data():kineticNodes.data();
        const double* innerWeights=momentDegree_>8?kineticRule20.second.data():kineticWeights.data();
        const std::size_t innerCount=momentDegree_>8?kineticRule20.first.size():kineticNodes.size();
        const long double width=static_cast<long double>(hi)-lo,height=static_cast<long double>(innerHi)-innerLo;
        for(std::size_t i=0;i<nodes.size();++i){
            const long double distance=width*(1.L-std::abs(nodes[i]))/2;
            const double outerPoint=static_cast<double>(nodes[i]>0?hi-distance:lo+distance);
            if(!(outerPoint>lo && outerPoint<hi))return bounded();
            if(++queries_>65536){
                std::ostringstream message;message.precision(17);
                message<<"[HydroMomentumGeometry] kinetic slice work limit exceeded at outer interval ["
                    <<lo<<','<<hi<<"], box=["<<lower_.x<<','<<upper_.x<<"]x["<<lower_.y<<','<<upper_.y
                    <<"], axis="<<outerAxis_<<", queries="<<queries_<<", native="<<native_;
                throw std::runtime_error(message.str());
            }
            const LevelSetMappedSegment slice(scalar_,onPlane(outerAxis_,outerPoint,innerLo),onPlane(outerAxis_,outerPoint,innerHi));
            require(slice.usesSource(*scalar_),"foreign kinetic phase source");
            for(const auto& range:slice.ranges())for(std::size_t j=0;j<innerCount;++j){
                const long double t=range.offset+range.length*(1.L+innerNodes[j])/2;
                const long double innerPoint=range.fromEnd?innerHi-height*t:innerLo+height*t;
                const Vector2D q=onPlane(outerAxis_,outerPoint,static_cast<double>(innerPoint));
                const long double radialWeight=radial_?CartesianMeasure::twoPi*
                    (outerAxis_==0?static_cast<long double>(outerPoint):innerPoint):1.L;
                const long double physical=width*weights[i]*height*range.length*innerWeights[j]/4*radialWeight;
                if(!(coord(q,innerAxis)>innerLo && coord(q,innerAxis)<innerHi)){
                    // A phase interval narrower than coordinate resolution
                    // has no representable GL8 point. Bound only that positive
                    // subinterval's quadrature mass; do not charge the whole
                    // outer rectangle or synthesize a phase sample.
                    rule.error+=physical;continue;
                }
                const int pointSign=scalar_->nativeQuarticBoxSign(native_,q,q);
                // A binary64 GL point can lie within the outward Q4 sign
                // enclosure of the interface even though its retained slice
                // interval has a certified phase. Reject only a resolved
                // contradiction; a zero point enclosure is not a phase label.
                if(pointSign!=0 && pointSign!=(range.negative?-1:1)){
                    std::ostringstream message;message.precision(17);
                    message<<"[HydroMomentumGeometry] kinetic quadrature point has unresolved/wrong phase at ("
                        <<q.x<<','<<q.y<<"), rangeNegative="<<range.negative<<", pointSign="<<pointSign
                        <<", sliceOuter="<<outerPoint<<", sliceAxis="<<outerAxis_;
                    throw std::runtime_error(message.str());
                }
                const double weight=static_cast<double>(physical);
                require(std::isfinite(weight) && weight>0,"nonpositive kinetic weight");
                rule.points.push_back({dual_,pressure_,q,weight,range.negative});
                require(rule.points.size()<=2000000,"kinetic point work limit exceeded");
                std::array<long double,35> xp{},yp{};xp[0]=yp[0]=1;
                const long double nx=(static_cast<long double>(q.x)-lower_.x)/(upper_.x-lower_.x);
                const long double ny=(static_cast<long double>(q.y)-lower_.y)/(upper_.y-lower_.y);
                const int n=momentDegree_+1;
                for(int k=1;k<n;++k){xp[k]=xp[k-1]*nx;yp[k]=yp[k-1]*ny;}
                for(int a=0;a<n;++a)for(int b=0;b<n;++b)
                    rule.moments[(range.negative?0:n*n)+n*a+b]+=weight*xp[a]*yp[b];
            }
        }
        return rule;
    }
    std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> scalar_;
    int dual_,pressure_,native_;
    Vector2D lower_,upper_;
    bool radial_;
    int outerAxis_=0;
    int momentDegree_=8;
    std::vector<double> nativeKnots_;
    std::size_t queries_=0;
};
long double kineticDot(const Eigen::VectorXd& a,const Eigen::VectorXd& b)
{
    require(a.size()==b.size() && a.allFinite() && b.allFinite(),"invalid kinetic vector product");
    long double sum=0;for(Eigen::Index i=0;i<a.size();++i)sum+=static_cast<long double>(a[i])*b[i];return sum;
}
}
struct HydroKineticHodge::Data
{
    Data(HydroMomentumMass m,HydroAdvectorReconstruction r):mass(std::move(m)),source(std::move(r)) {}
    HydroMomentumMass mass;
    HydroAdvectorReconstruction source;
    Matrix matrix;
    Matrix augmented;
    std::vector<bool> marked;
    std::vector<HydroKineticHodgePoint> points;
    std::vector<HydroKineticHodgePoint> strainPoints;
    HydroKineticHodgeReport report;
};
HydroKineticHodge::HydroKineticHodge(std::shared_ptr<const Data> d):data_(std::move(d)) {}
HydroKineticHodge HydroKineticHodge::build(const QuadtreeMesh& mesh,HydroMomentumMass mass,
    HydroAdvectorReconstruction source)
{
    mass.validateFor(mesh,mass.geometryIdentity(),mass.timeLayer());
    require(mass.representation()==HydroMomentumMass::DensityRepresentation::SharpGeometry,
        "kinetic Hodge requires sharp material cubature");
    require(source.momentumGeometry().sameIdentity(mass.geometry()),"different kinetic source/momentum owner");
    if(const auto* bound=source.phaseGeometry())bound->validateFor(mesh,mass.geometryIdentity(),mass.timeLayer());
    const auto& regular=source.regularGeometry();regular.validateFor(mesh,mass.geometry());
    auto d=std::make_shared<Data>(std::move(mass),std::move(source));
    const auto& layout=d->mass.layout();const bool rz=layout.domain().axisymmetric();
    const auto& phase=d->mass.phaseGeometry();const int nu=layout.velocityCount(),np=layout.pressureCount();
    const int nm=static_cast<int>(d->source.tangentialModes().size());d->report.tangentialModes=nm;
    std::vector<bool> near(np,false);
    for(int p=0;p<np;++p)near[p]=phase.cellVolumes()[p].liquid>0 && phase.cellVolumes()[p].gas>0;
    for(int u=0;u<nu;++u)if(phase.faceTraces()[u].kind!=HydroPhaseFaceTraceKind::Ordinary){
        const auto& f=regular.faces()[u];if(f.negativePressure>=0)near[f.negativePressure]=true;
        if(f.positivePressure>=0)near[f.positivePressure]=true;
    }
    const auto seed=near; // Fixed one-cell reconstruction halo, independent of case/rho.
    for(const auto& f:regular.faces())if(f.negativePressure>=0 && f.positivePressure>=0 &&
        (seed[f.negativePressure] || seed[f.positivePressure]))near[f.negativePressure]=near[f.positivePressure]=true;
    d->marked.assign(nu,false);
    std::vector<std::optional<Detail::HydroRegularCellStencil>> stencils(np);
    const auto& geometry=phase.publication().state().geometry();
    const std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> scalar(
        std::make_shared<const Interface::HermiteGalsGeometryView>(geometry),&geometry.snapshot());
    std::map<std::pair<int,int>,long double> entries;
    for(int u=0;u<nu;++u){
        const auto& volume=d->mass.geometry().volumes()[u];
        for(const auto& overlap:volume.overlaps)d->marked[u]=d->marked[u] || near[overlap.pressure];
        const auto& face=regular.faces()[u];
        // The zero-area axis slot is prescribed zero. Diagonal completion
        // keeps the full algebraic block SPD without assigning physical flow.
        const int component=layout.normalAxis(u)==HydroDofLayout::NormalAxis::First?0:1;
        if(!d->marked[u]){entries[{u,u}]+=d->mass.masses()[u];if(component==0)continue;}
        if(d->marked[u]&&face.area==0)entries[{u,u}]+=d->mass.masses()[u];
        if(d->marked[u])++d->report.markedDuals;
        HydroPhaseMeasure measured;
        for(const auto& overlap:volume.overlaps){
            const int p=overlap.pressure;
            if(!stencils[p])stencils[p]=regular.cellStencil(mesh,p);
            std::vector<double> xs{overlap.lower.x,overlap.upper.x},ys{overlap.lower.y,overlap.upper.y};
            for(double x:regular.cells()[p].supportKnots[0])if(x>overlap.lower.x && x<overlap.upper.x)xs.push_back(x);
            for(double y:regular.cells()[p].supportKnots[1])if(y>overlap.lower.y && y<overlap.upper.y)ys.push_back(y);
            for(int sourceCell:scalar->nativeCellCandidates(overlap.lower,overlap.upper)){
                const auto bounds=scalar->nativeCellBounds(sourceCell);
                for(double x:{bounds[0].x,bounds[1].x})if(x>overlap.lower.x && x<overlap.upper.x)xs.push_back(x);
                for(double y:{bounds[0].y,bounds[1].y})if(y>overlap.lower.y && y<overlap.upper.y)ys.push_back(y);
            }
            std::sort(xs.begin(),xs.end());xs.erase(std::unique(xs.begin(),xs.end()),xs.end());
            std::sort(ys.begin(),ys.end());ys.erase(std::unique(ys.begin(),ys.end()),ys.end());
            for(std::size_t x=1;x<xs.size();++x)for(std::size_t y=1;y<ys.size();++y){
                const Vector2D lo{xs[x-1],ys[y-1]},hi{xs[x],ys[y]};
                const long double cap=capacity(lo,hi,rz);
                const int localMode=d->source.tangentialModeForPressure(p);
                KineticBoxRule box(scalar,u,p,layout.nativeCell(p),lo,hi,rz,localMode>=0?(rz?34:20):8);
                auto rule=box.integrate(coord(lo,box.outerAxis()),coord(hi,box.outerAxis()),cap*1e-11L);
                d->report.maxRelativeMomentEstimate=std::max(d->report.maxRelativeMomentEstimate,rule.error/cap);
                for(const auto& q:rule.points){
                    (q.liquid?measured.liquid:measured.gas)+=q.physicalWeight;
                    if(!d->marked[u])continue;
                    const long double weight=q.physicalWeight*(q.liquid?d->mass.liquidDensity():d->mass.gasDensity());
                    const auto weights=regular.cellJetWeights(*stencils[p],q.point);
                    std::vector<std::pair<int,double>> active;
                    for(const auto& w:weights)if(w.jet[component]!=0)active.emplace_back(w.velocityDof,w.jet[component]);
                    const int mode=d->source.tangentialModeForPressure(p);
                    if(mode>=0){const auto jet=d->source.tangentialModeJet(mode,q.point,q.liquid);
                        const double value=component==0?jet.value.x:jet.value.y;
                        if(value!=0)active.emplace_back(nu+mode,value);}
                    d->report.sourceWeights+=active.size();
                    for(std::size_t a=0;a<active.size();++a)for(std::size_t b=a;b<active.size();++b){
                        const auto [i,wi]=active[a];const auto [j,wj]=active[b];
                        const auto key=std::minmax(i,j);entries[{key.first,key.second}]+=weight*wi*wj;
                    }
                }
                if(component==1)d->strainPoints.insert(d->strainPoints.end(),rule.points.begin(),rule.points.end());
                if(d->marked[u])d->points.insert(d->points.end(),rule.points.begin(),rule.points.end());
                require(d->points.size()<=2000000,"total kinetic point work limit exceeded");
            }
        }
        const auto expected=d->mass.phaseVolumes()[u];
        const auto defect=std::max(std::abs(measured.liquid-expected.liquid),std::abs(measured.gas-expected.gas))/volume.volume;
        require(std::isfinite(defect) && defect<=1e-9L,"kinetic dual quadrature differs from authoritative phase volume");
        d->report.maxRelativePhaseVolumeDifference=std::max(d->report.maxRelativePhaseVolumeDifference,defect);
    }
    std::vector<Eigen::Triplet<double,Eigen::Index>> triplets;
    for(const auto& entry:entries){
        const auto [i,j]=entry.first;const double value=static_cast<double>(entry.second);
        require(std::isfinite(value),"unrepresentable kinetic matrix entry");
        if(value==0)continue;triplets.emplace_back(i,j,value);if(i!=j)triplets.emplace_back(j,i,value);
    }
    d->augmented.resize(nu+nm,nu+nm);d->augmented.setFromTriplets(triplets.begin(),triplets.end());d->augmented.makeCompressed();
    d->matrix=d->augmented.topLeftCorner(nu,nu);d->matrix.makeCompressed();
    const Matrix difference=d->augmented-Matrix(d->augmented.transpose());
    d->report.symmetryDefect=difference.norm()/d->augmented.norm();
    require(std::isfinite(d->report.symmetryDefect) && d->report.symmetryDefect<=1e-12,"nonsymmetric kinetic Hodge");
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> factor;
    factor.compute(Eigen::SparseMatrix<double>(d->augmented));
    require(factor.info()==Eigen::Success,"kinetic Hodge is not positive definite");
    d->report.points=d->points.size();d->report.nonZeros=static_cast<std::size_t>(d->augmented.nonZeros());
    return HydroKineticHodge(std::move(d));
}
const HydroMomentumMass& HydroKineticHodge::materialMass() const noexcept{return data_->mass;}
const HydroAdvectorReconstruction& HydroKineticHodge::reconstruction() const noexcept{return data_->source;}
const HydroDofLayout& HydroKineticHodge::layout() const noexcept{return materialMass().layout();}
Interface::HermiteGalsIdentity HydroKineticHodge::geometryIdentity() const{return materialMass().geometryIdentity();}
LAMP::TimeLayer HydroKineticHodge::timeLayer() const{return materialMass().timeLayer();}
const HydroKineticHodge::Matrix& HydroKineticHodge::matrix() const noexcept{return data_->matrix;}
const HydroKineticHodge::Matrix& HydroKineticHodge::augmentedMatrix() const noexcept{return data_->augmented;}
const std::vector<bool>& HydroKineticHodge::markedDuals() const noexcept{return data_->marked;}
const std::vector<HydroKineticHodgePoint>& HydroKineticHodge::points() const noexcept{return data_->points;}
const std::vector<HydroKineticHodgePoint>& HydroKineticHodge::strainQuadrature() const noexcept{return data_->strainPoints;}
const HydroKineticHodgeReport& HydroKineticHodge::report() const noexcept{return data_->report;}
void HydroKineticHodge::validateFor(const QuadtreeMesh& mesh,Interface::HermiteGalsIdentity id,
    LAMP::TimeLayer time) const
{materialMass().validateFor(mesh,id,time);reconstruction().regularGeometry().validateFor(mesh,materialMass().geometry());}
std::vector<double> HydroKineticHodge::apply(const QuadtreeMesh& mesh,const HydroNativeValues& values,
    Interface::HermiteGalsIdentity id,LAMP::TimeLayer time) const
{
    validateFor(mesh,id,time);values.validateFor(layout(),mesh);
    require(reconstruction().tangentialModes().empty(),"explicit tangential state required for kinetic action");
    for(int u=0;u<layout().velocityCount();++u)if(reconstruction().regularGeometry().faces()[u].area==0)
        require(values.normalVelocity()[u]==0,"nonzero prescribed axis velocity in kinetic field");
    const Eigen::Map<const Eigen::VectorXd> v(values.normalVelocity().data(),layout().velocityCount());
    const Eigen::VectorXd result=matrix()*v;require(result.allFinite(),"nonfinite kinetic action");
    return {result.data(),result.data()+result.size()};
}
long double HydroKineticHodge::kineticEnergy(const QuadtreeMesh& mesh,const HydroNativeValues& values,
    Interface::HermiteGalsIdentity id,LAMP::TimeLayer time) const
{
    const auto action=apply(mesh,values,id,time);long double energy=0;
    for(std::size_t u=0;u<action.size();++u)energy+=static_cast<long double>(values.normalVelocity()[u])*action[u]/2;
    require(std::isfinite(energy) && energy>=0,"negative/nonfinite kinetic energy");return energy;
}
std::vector<double> HydroKineticHodge::apply(const QuadtreeMesh& mesh,const HydroSpatialAdvector& state) const
{
    validateFor(mesh,state.geometryIdentity(),state.timeLayer());
    state.validateFor(mesh,materialMass().geometry(),geometryIdentity(),timeLayer());
    require(state.reconstruction().sameIdentity(reconstruction()),"different tangential kinetic source owner");
    const int nu=layout().velocityCount(),nm=static_cast<int>(reconstruction().tangentialModes().size());
    require(state.tangentialVelocity().size()==static_cast<std::size_t>(nm),"kinetic tangential state size mismatch");
    Eigen::VectorXd x(nu+nm);
    for(int u=0;u<nu;++u){x[u]=state.values().normalVelocity()[u];
        if(reconstruction().regularGeometry().faces()[u].area==0)require(x[u]==0,"nonzero prescribed axis velocity in kinetic field");}
    for(int m=0;m<nm;++m)x[nu+m]=state.tangentialVelocity()[m];
    const Eigen::VectorXd result=augmentedMatrix()*x;require(result.allFinite(),"nonfinite augmented kinetic action");
    return {result.data(),result.data()+result.size()};
}
long double HydroKineticHodge::kineticEnergy(const QuadtreeMesh& mesh,const HydroSpatialAdvector& state) const
{
    const auto action=apply(mesh,state);long double energy=0;const int nu=layout().velocityCount();
    for(int u=0;u<nu;++u)energy+=static_cast<long double>(state.values().normalVelocity()[u])*action[u]/2;
    for(std::size_t m=0;m<state.tangentialVelocity().size();++m)
        energy+=static_cast<long double>(state.tangentialVelocity()[m])*action[nu+m]/2;
    require(std::isfinite(energy)&&energy>=0,"negative/nonfinite augmented kinetic energy");return energy;
}
std::vector<double> HydroKineticHodge::bodyForceLoad(const QuadtreeMesh& mesh,const std::vector<Vector2D>& force) const
{
    validateFor(mesh,geometryIdentity(),timeLayer());require(force.size()==strainQuadrature().size(),"body force sample inventory mismatch");
    const auto& source=reconstruction();const auto& regular=source.regularGeometry();const int nu=layout().velocityCount();
    std::vector<long double> sum(nu+source.tangentialModes().size());
    std::vector<std::optional<Detail::HydroRegularCellStencil>> stencils(layout().pressureCount());
    for(std::size_t i=0;i<force.size();++i){const auto& q=strainQuadrature()[i];const auto f=force[i];
        require(std::isfinite(f.x)&&std::isfinite(f.y),"nonfinite body force sample");
        if(!stencils[q.pressure])stencils[q.pressure]=regular.cellStencil(mesh,q.pressure);
        for(const auto& row:regular.cellJetWeights(*stencils[q.pressure],q.point))
            sum[row.velocityDof]+=q.physicalWeight*(static_cast<long double>(f.x)*row.jet[0]+static_cast<long double>(f.y)*row.jet[1]);
        const int mode=source.tangentialModeForPressure(q.pressure);
        if(mode>=0){const auto jet=source.tangentialModeJet(mode,q.point,q.liquid);
            sum[nu+mode]+=q.physicalWeight*(static_cast<long double>(f.x)*jet.value.x+static_cast<long double>(f.y)*jet.value.y);}
    }
    std::vector<double> load(sum.size());
    for(std::size_t i=0;i<sum.size();++i){load[i]=static_cast<double>(sum[i]);require(std::isfinite(load[i]),"nonfinite augmented body load");}
    return load;
}
HydroKineticHodgeEnergyChange HydroKineticHodge::energyChange(const QuadtreeMesh& mesh,
    const HydroKineticHodge& begin,const HydroKineticHodge& end,const HydroSpatialAdvector& before,
    const HydroSpatialAdvector& after) const
{
    for(const auto* h:{&begin,this,&end}){
        h->validateFor(mesh,h->geometryIdentity(),h->timeLayer());
        require(h->materialMass().geometry().sameIdentity(materialMass().geometry()),"different endpoint kinetic momentum owner");
        require(h->materialMass().liquidDensity()==materialMass().liquidDensity()&&h->materialMass().gasDensity()==materialMass().gasDensity(),
            "different endpoint kinetic materials");
        require(h->augmentedMatrix().rows()==augmentedMatrix().rows(),"modal metric size changed without physical transfer");
        const auto& modes=h->reconstruction().tangentialModes();const auto& actual=reconstruction().tangentialModes();
        require(modes.size()==actual.size(),"modal topology changed without physical transfer");
        for(std::size_t m=0;m<modes.size();++m)require(modes[m].pressure==actual[m].pressure,
            "modal ordering changed without physical transfer");
    }
    const double t0=begin.timeLayer().time(),t1=end.timeLayer().time(),tm=timeLayer().time();
    require(t1>t0&&std::abs(tm-(t0+(t1-t0)/2))<=8*std::numeric_limits<double>::epsilon()*std::max({1.,std::abs(t0),std::abs(t1)}),
        "invalid kinetic endpoint/midpoint clocks");
    HydroKineticHodgeEnergyChange result;
    result.beginEnergy=begin.kineticEnergy(mesh,before);result.endEnergy=end.kineticEnergy(mesh,after);
    const auto pack=[&](const HydroSpatialAdvector& state){Eigen::VectorXd x(augmentedMatrix().rows());const int nu=layout().velocityCount();
        for(int u=0;u<nu;++u)x[u]=state.values().normalVelocity()[u];
        for(std::size_t m=0;m<state.tangentialVelocity().size();++m)x[nu+m]=state.tangentialVelocity()[m];return x;};
    const Eigen::VectorXd x0=pack(before),x1=pack(after),mid=(x0+x1)/2,delta=x1-x0;
    result.energyChange=result.endEnergy-result.beginEnergy;
    result.midpointInertialWork=kineticDot(mid,augmentedMatrix()*delta);
    result.metricVariation=kineticDot(mid,end.augmentedMatrix()*mid-begin.augmentedMatrix()*mid)/2;
    result.endpointQuadratureRemainder=kineticDot(mid,(begin.augmentedMatrix()*delta+end.augmentedMatrix()*delta)/2-augmentedMatrix()*delta)+
        kineticDot(delta,end.augmentedMatrix()*delta-begin.augmentedMatrix()*delta)/8;
    result.identityResidual=result.energyChange-result.midpointInertialWork-result.metricVariation-result.endpointQuadratureRemainder;
    return result;
}
HydroKineticHodgeEnergyChange HydroKineticHodge::energyChange(const QuadtreeMesh& mesh,
    const HydroKineticHodge& begin,const HydroKineticHodge& end,const HydroNativeValues& before,
    const HydroNativeValues& after) const
{
    for(const auto* h:{&begin,this,&end}){
        h->validateFor(mesh,h->geometryIdentity(),h->timeLayer());
        require(h->materialMass().geometry().sameIdentity(materialMass().geometry()),"different endpoint kinetic momentum owner");
        require(h->materialMass().liquidDensity()==materialMass().liquidDensity() &&
            h->materialMass().gasDensity()==materialMass().gasDensity(),"different endpoint kinetic materials");
    }
    const double t0=begin.timeLayer().time(),t1=end.timeLayer().time(),tm=timeLayer().time();
    require(t1>t0 && std::abs(tm-(t0+(t1-t0)/2))<=8*std::numeric_limits<double>::epsilon()*
        std::max({1.,std::abs(t0),std::abs(t1)}),"invalid kinetic endpoint/midpoint clocks");
    before.validateFor(layout(),mesh);after.validateFor(layout(),mesh);
    HydroKineticHodgeEnergyChange result;
    result.beginEnergy=begin.kineticEnergy(mesh,before,begin.geometryIdentity(),begin.timeLayer());
    result.endEnergy=end.kineticEnergy(mesh,after,end.geometryIdentity(),end.timeLayer());
    result.energyChange=result.endEnergy-result.beginEnergy;
    const Eigen::Map<const Eigen::VectorXd> u0(before.normalVelocity().data(),layout().velocityCount()),
        u1(after.normalVelocity().data(),layout().velocityCount());
    const Eigen::VectorXd mid=(u0+u1)/2,delta=u1-u0;
    result.midpointInertialWork=kineticDot(mid,matrix()*delta);
    result.metricVariation=kineticDot(mid,end.matrix()*mid-begin.matrix()*mid)/2;
    result.endpointQuadratureRemainder=kineticDot(mid,(begin.matrix()*delta+end.matrix()*delta)/2-matrix()*delta)+
        kineticDot(delta,end.matrix()*delta-begin.matrix()*delta)/8;
    result.identityResidual=result.energyChange-result.midpointInertialWork-result.metricVariation-result.endpointQuadratureRemainder;
    return result;
}
}
