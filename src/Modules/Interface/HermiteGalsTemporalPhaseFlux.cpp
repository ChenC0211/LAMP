#include "Modules/Interface/HermiteGalsTemporalPhaseFlux.h"
#include "Core/Discrete/BoundedCirculation.h"
#include "Modules/Interface/HermiteGalsPhaseRegion.h"
#include "Modules/Interface/HermiteGalsInventoryRegion.h"
#include "Modules/Interface/Detail/HermiteGalsInventoryAssembly.h"

#include "Core/Mesh/NativeMeshIdentity.h"
#include "Core/Mesh/FaceSegmentGraph.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Geometry/LevelSetHermiteGeometry.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <queue>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>

namespace LAMP::Interface
{
namespace
{
void require(bool ok, const char* message)
{ if (!ok) throw std::invalid_argument(std::string("[Interface/HermiteGalsTemporalPhaseFlux] ") + message); }
void validInterval(double begin, double end)
{
    require(std::isfinite(begin) && std::isfinite(end) && end > begin
        && std::isfinite(end - begin), "invalid finite positive interval");
}
void validVolumes(const HermiteGalsIntegratedFaceVolumes& v)
{
    for (long double x : {v.integratedTotalPositiveVolume, v.integratedTotalNegativeVolume,
            v.integratedLiquidPositiveVolume, v.integratedLiquidNegativeVolume})
        require(std::isfinite(x) && x >= 0, "integrated gross volume must be finite and nonnegative");
    require(v.integratedLiquidPositiveVolume <= v.integratedTotalPositiveVolume
        && v.integratedLiquidNegativeVolume <= v.integratedTotalNegativeVolume,
        "integrated liquid gross volume exceeds total");
}
long double average(long double volume, double dt)
{
    const auto result = volume / static_cast<long double>(dt);
    require(std::isfinite(result), "derived average rate is not representable");
    return result;
}
bool sameInflow(const HermiteGalsConstantInflowSource* a,const HermiteGalsConstantInflowSource* b)
{return a && b?a->sameSource(*b):a==b;}
bool sameBoundary(const HermiteGalsBoundaryPhaseSource* a,const HermiteGalsBoundaryPhaseSource* b)
{return a && b?a->sameSource(*b):a==b;}
bool isInflow(const NativeLeafFace& face,Vector2D velocity)
{
    const double speed=face.horizontal?velocity.y:velocity.x;
    return (face.negativeCell<0 && speed>0) || (face.positiveCell<0 && speed<0);
}
Vector2D constantTranslation(const QuadtreeMesh& mesh,const HermiteGalsPhaseFluxIntervalSource& interval)
{
    interval.validateMesh(mesh);
    const auto velocity=interval.velocitySource().constantValue();
    require(velocity.has_value(),"inflow translation requires an owned constant velocity source");
    require(velocity->x==0 || velocity->y==0,"inflow translation requires single-axis velocity");
    require(interval.measure()!=HermiteGalsMeasure::AxisymmetricVolume || velocity->x==0,"RZ inflow translation requires axial velocity");
    require(!interval.characteristicSource().domain().activeSegment,"inflow translation cannot certify opaque coverage");
    require(interval.previousAccepted().primary().geometryCells.size()==static_cast<std::size_t>(mesh.cellIndexNum()),
        "inflow translation requires full active scalar support");
    return *velocity;
}
}

struct HermiteGalsTemporalVelocitySource::Data
{
    HermiteGalsVelocityView previous, arrival;
    std::optional<Vector2D> constant;
    std::shared_ptr<const LAMP::Geometry::PolynomialVelocitySlab> polynomial;
    Data(const HermiteGalsVelocityView& a, const HermiteGalsVelocityView& b,
        std::optional<Vector2D> value = {},std::shared_ptr<const LAMP::Geometry::PolynomialVelocitySlab> field = {})
        : previous(a), arrival(b), constant(value),polynomial(std::move(field)) {}
};
HermiteGalsTemporalVelocitySource::HermiteGalsTemporalVelocitySource(
    const HermiteGalsState& accepted, HermiteGalsIdentity arrivalGeometry, HermiteGalsTimeLayer arrivalTime,
    HermiteGalsVelocityView::Query previousQuery, HermiteGalsVelocityView::Query arrivalQuery)
{
    const HermiteGalsVelocityView previous(accepted.identity(), accepted.timeLayer(), std::move(previousQuery));
    const HermiteGalsVelocityView arrival(arrivalGeometry, arrivalTime, std::move(arrivalQuery));
    const auto a = previous.timeLayer(), b = arrival.timeLayer();
    validInterval(a.time(), b.time());
    require(a.position != HermiteGalsTimeLayer::Position::Midpoint && b.begin == a.time(),
        "velocity source requires compatible previous/arrival times");
    data_ = std::make_shared<const Data>(previous, arrival);
}
const HermiteGalsTemporalVelocitySource::Data& HermiteGalsTemporalVelocitySource::checked() const
{ require(bool(data_), "moved-from velocity source"); return *data_; }
const HermiteGalsVelocityView& HermiteGalsTemporalVelocitySource::previous() const { return checked().previous; }
const HermiteGalsVelocityView& HermiteGalsTemporalVelocitySource::arrival() const { return checked().arrival; }
bool HermiteGalsTemporalVelocitySource::sameSource(const HermiteGalsTemporalVelocitySource& other) const noexcept
{ return data_ && data_ == other.data_; }
HermiteGalsTemporalVelocitySource HermiteGalsTemporalVelocitySource::constant(
    const HermiteGalsState& accepted, HermiteGalsIdentity arrivalGeometry,
    HermiteGalsTimeLayer arrivalTime, Vector2D value)
{
    require(std::isfinite(value.x) && std::isfinite(value.y), "nonfinite constant velocity");
    const auto query = [value](const Vector2D&) { return LevelSetBoundaryTrace::VelocityJet{value,{0,0,0,0}}; };
    HermiteGalsTemporalVelocitySource result(accepted, arrivalGeometry, arrivalTime, query, query);
    result.data_ = std::make_shared<const Data>(result.previous(), result.arrival(), value);
    return result;
}
std::optional<Vector2D> HermiteGalsTemporalVelocitySource::constantValue() const { return checked().constant; }
HermiteGalsTemporalVelocitySource HermiteGalsTemporalVelocitySource::polynomial(
    const HermiteGalsState& accepted,HermiteGalsIdentity arrivalGeometry,
    HermiteGalsTimeLayer arrivalTime,const LAMP::Geometry::PolynomialVelocityData& data)
{
    const auto field=std::make_shared<const LAMP::Geometry::PolynomialVelocitySlab>(accepted.timeLayer().time(),arrivalTime.time(),data);
    HermiteGalsTemporalVelocitySource result(accepted,arrivalGeometry,arrivalTime,
        [field](Vector2D point){return field->evaluate(point,field->timeBegin());},
        [field](Vector2D point){return field->evaluate(point,field->timeEnd());});
    result.data_=std::make_shared<const Data>(result.previous(),result.arrival(),std::nullopt,field);return result;
}
const LAMP::Geometry::PolynomialVelocitySlab* HermiteGalsTemporalVelocitySource::polynomialVelocity() const
{return checked().polynomial.get();}
HermiteGalsVelocityView HermiteGalsTemporalVelocitySource::polynomialArrivalView(
    HermiteGalsIdentity actualCandidate,HermiteGalsTimeLayer actualArrival) const
{
    const auto& d=checked();require(bool(d.polynomial),"substage view requires an owned polynomial velocity");
    const auto previousTime=d.previous.timeLayer().time(),arrivalTime=actualArrival.time();const auto full=d.arrival.timeLayer();
    require(actualArrival.begin==previousTime && actualArrival.step==full.step && actualArrival.end<=full.end
        && arrivalTime>previousTime && arrivalTime<=full.time(),"polynomial substage outside source time/step extent");
    return HermiteGalsVelocityView(actualCandidate,actualArrival,[owner=data_,arrivalTime](Vector2D point){
        return owner->polynomial->evaluate(point,arrivalTime);
    });
}
HermiteGalsVelocityView HermiteGalsTemporalVelocitySource::constantArrivalView(
    HermiteGalsIdentity actualCandidate,HermiteGalsTimeLayer actualArrival) const
{
    const auto& d=checked();require(d.constant.has_value(),"substage view requires an owned constant velocity");
    const auto previousTime=d.previous.timeLayer().time(),arrivalTime=actualArrival.time();
    const auto full=d.arrival.timeLayer();
    require(actualArrival.begin==previousTime && actualArrival.step==full.step
        && actualArrival.end<=full.end && arrivalTime>previousTime && arrivalTime<=full.time(),
        "constant substage outside source time/step extent");
    return HermiteGalsVelocityView(actualCandidate,actualArrival,[owner=data_](const Vector2D&) {
        return LevelSetBoundaryTrace::VelocityJet{*owner->constant,{0,0,0,0}};
    });
}

struct HermiteGalsTemporalCharacteristicSource::Data
{
    NativeMeshIdentity mesh;
    LevelSetBoundaryTrace::Domain domain;
    Data(const QuadtreeMesh& m, const LevelSetBoundaryTrace::Domain& d) : mesh(m), domain(d) {}
};
HermiteGalsTemporalCharacteristicSource::HermiteGalsTemporalCharacteristicSource(
    const QuadtreeMesh& mesh, const LevelSetBoundaryTrace::Domain& domain)
{
    require(domain.minCoord.x == mesh.minCoord.x && domain.minCoord.y == mesh.minCoord.y
        && domain.maxCoord.x == mesh.maxCoord.x && domain.maxCoord.y == mesh.maxCoord.y,
        "characteristic domain differs from mesh bounds");
    require(std::isfinite(domain.cellSize.x) && std::isfinite(domain.cellSize.y)
        && domain.cellSize.x > 0 && domain.cellSize.y > 0
        && domain.cellSize.x <= mesh.maxCoord.x - mesh.minCoord.x
        && domain.cellSize.y <= mesh.maxCoord.y - mesh.minCoord.y,
        "invalid characteristic domain cell size");
    data_ = std::make_shared<const Data>(mesh, domain);
}
const HermiteGalsTemporalCharacteristicSource::Data& HermiteGalsTemporalCharacteristicSource::checked() const
{ require(bool(data_), "moved-from characteristic source"); return *data_; }
const LevelSetBoundaryTrace::Domain& HermiteGalsTemporalCharacteristicSource::domain() const { return checked().domain; }
const NativeMeshIdentity& HermiteGalsTemporalCharacteristicSource::meshIdentity() const { return checked().mesh; }
HermiteGalsTemporalCharacteristicSource::Method HermiteGalsTemporalCharacteristicSource::method() const
{ checked(); return Method::CoherentExplicitMidpointRK2; }
void HermiteGalsTemporalCharacteristicSource::validateMesh(const QuadtreeMesh& mesh) const { checked().mesh.validate(mesh); }
bool HermiteGalsTemporalCharacteristicSource::sameSource(const HermiteGalsTemporalCharacteristicSource& other) const noexcept
{ return data_ && data_ == other.data_; }

struct HermiteGalsPhaseFluxIntervalSource::Data
{
    HermiteGalsState previous;
    HermiteGalsTemporalVelocitySource velocity;
    HermiteGalsTemporalCharacteristicSource characteristic;
    NativeLeafDomain domain;
    std::optional<HermiteGalsConstantInflowSource> inflow;
    std::optional<HermiteGalsBoundaryPhaseSource> boundary;
    double begin, end;
    Data(const HermiteGalsState& p, const HermiteGalsTemporalVelocitySource& v,
        const HermiteGalsTemporalCharacteristicSource& c, NativeLeafDomain d, double a, double b)
        : previous(p), velocity(v), characteristic(c), domain(std::move(d)), begin(a), end(b) {}
};
HermiteGalsPhaseFluxIntervalSource::HermiteGalsPhaseFluxIntervalSource(const QuadtreeMesh& mesh,
    const HermiteGalsState& previous, const HermiteGalsTemporalVelocitySource& velocity,
    const HermiteGalsTemporalCharacteristicSource& characteristic, double begin, double end)
{
    validInterval(begin, end);
    previous.validateMesh(mesh); characteristic.validateMesh(mesh);
    const auto time = previous.timeLayer();
    require(time.position != HermiteGalsTimeLayer::Position::Midpoint, "previous authority is a midpoint candidate");
    const auto vt = velocity.previous().timeLayer();
    require(velocity.previous().geometryIdentity() == previous.identity()
        && vt.step == time.step && vt.begin == time.begin && vt.end == time.end && vt.position == time.position,
        "velocity source differs from previous geometry/time authority");
    require(begin >= time.time() && end <= velocity.arrival().timeLayer().time(),
        "interval outside retained velocity trajectory");
    const auto& p = previous.primary();
    const bool radial = p.measure == HermiteGalsMeasure::AxisymmetricVolume;
    require(characteristic.domain().radialAxisSymmetry == radial, "characteristic/geometry measure mismatch");
    require(!p.geometryCells.empty(), "interval requires retained geometry cells");
    FaceSegmentGraph graph;
    graph.build(mesh, [&p](const QuadtreeMesh&, const Cell& cell) {
        return std::binary_search(p.geometryCells.begin(), p.geometryCells.end(), cell.index());
    });
    auto domain = NativeLeafDomain::build(mesh, graph, radial);
    data_ = std::make_shared<const Data>(previous, velocity, characteristic, std::move(domain), begin, end);
}
HermiteGalsPhaseFluxIntervalSource::HermiteGalsPhaseFluxIntervalSource(const QuadtreeMesh& mesh,
    const HermiteGalsState& previous,const HermiteGalsTemporalVelocitySource& velocity,
    const HermiteGalsTemporalCharacteristicSource& characteristic,double begin,double end,
    const HermiteGalsConstantInflowSource& inflow)
    :HermiteGalsPhaseFluxIntervalSource(mesh,previous,velocity,characteristic,begin,end)
{
    inflow.validateInterval(mesh,*this);
    auto bound=std::make_shared<Data>(*data_);bound->inflow=inflow;data_=std::move(bound);
}
HermiteGalsPhaseFluxIntervalSource::HermiteGalsPhaseFluxIntervalSource(const QuadtreeMesh& mesh,
    const HermiteGalsState& previous,const HermiteGalsTemporalVelocitySource& velocity,
    const HermiteGalsTemporalCharacteristicSource& characteristic,double begin,double end,
    const HermiteGalsBoundaryPhaseSource& boundary)
    :HermiteGalsPhaseFluxIntervalSource(mesh,previous,velocity,characteristic,begin,end)
{
    boundary.validateInterval(mesh,*this);
    auto bound=std::make_shared<Data>(*data_);bound->boundary=boundary;data_=std::move(bound);
}
const HermiteGalsPhaseFluxIntervalSource::Data& HermiteGalsPhaseFluxIntervalSource::checked() const
{ require(bool(data_), "moved-from interval source"); return *data_; }
const HermiteGalsState& HermiteGalsPhaseFluxIntervalSource::previousAccepted() const { return checked().previous; }
HermiteGalsIdentity HermiteGalsPhaseFluxIntervalSource::previousGeometryIdentity() const { return checked().previous.identity(); }
const HermiteGalsTemporalVelocitySource& HermiteGalsPhaseFluxIntervalSource::velocitySource() const { return checked().velocity; }
const HermiteGalsTemporalCharacteristicSource& HermiteGalsPhaseFluxIntervalSource::characteristicSource() const { return checked().characteristic; }
const HermiteGalsConstantInflowSource* HermiteGalsPhaseFluxIntervalSource::inflowSource() const
{const auto& d=checked();return d.inflow?&*d.inflow:nullptr;}
const HermiteGalsBoundaryPhaseSource* HermiteGalsPhaseFluxIntervalSource::boundaryPhaseSource() const
{const auto& d=checked();return d.boundary?&*d.boundary:nullptr;}
const NativeMeshIdentity& HermiteGalsPhaseFluxIntervalSource::meshIdentity() const { return checked().characteristic.meshIdentity(); }
HermiteGalsMeasure HermiteGalsPhaseFluxIntervalSource::measure() const { return checked().previous.primary().measure; }
double HermiteGalsPhaseFluxIntervalSource::timeBegin() const { return checked().begin; }
double HermiteGalsPhaseFluxIntervalSource::timeEnd() const { return checked().end; }
double HermiteGalsPhaseFluxIntervalSource::dt() const { return checked().end - checked().begin; }
bool HermiteGalsPhaseFluxIntervalSource::sameInterval(const HermiteGalsPhaseFluxIntervalSource& other) const
{
    const auto& a = checked(); const auto& b = other.checked();
    return a.previous.identity() == b.previous.identity() && a.velocity.sameSource(b.velocity)
        && a.characteristic.sameSource(b.characteristic) && a.begin == b.begin && a.end == b.end
        && sameInflow(inflowSource(),other.inflowSource()) && sameBoundary(boundaryPhaseSource(),other.boundaryPhaseSource());
}
void HermiteGalsPhaseFluxIntervalSource::validateMesh(const QuadtreeMesh& mesh) const
{ checked().previous.validateMesh(mesh); checked().characteristic.validateMesh(mesh); checked().domain.validateMeshIdentity(mesh); }
NativeLeafFace HermiteGalsPhaseFluxIntervalSource::face(const QuadtreeMesh& mesh, int nativeFace) const
{
    validateMesh(mesh);
    const auto& d = checked().domain;
    require(nativeFace >= 0 && nativeFace < d.nativeFaceCount(), "invalid native face identity");
    const auto f = d.face(mesh, nativeFace);
    require(f.negativeCell >= 0 || f.positiveCell >= 0, "inactive native face identity");
    return f;
}

struct HermiteGalsCharacteristicEnclosure::Data
{
    HermiteGalsPhaseFluxIntervalSource interval;
    LAMP::Geometry::PolynomialRk2Bounds bounds;
    Data(const HermiteGalsPhaseFluxIntervalSource& source,const LAMP::Geometry::PolynomialRk2Bounds& value):interval(source),bounds(value){}
};
HermiteGalsCharacteristicEnclosure::HermiteGalsCharacteristicEnclosure(std::shared_ptr<const Data> data):data_(std::move(data)){}
const HermiteGalsCharacteristicEnclosure::Data& HermiteGalsCharacteristicEnclosure::checked() const
{require(bool(data_),"moved-from characteristic enclosure");return *data_;}
const HermiteGalsPhaseFluxIntervalSource& HermiteGalsCharacteristicEnclosure::interval() const{return checked().interval;}
const LAMP::Geometry::PolynomialRk2Bounds& HermiteGalsCharacteristicEnclosure::bounds() const{return checked().bounds;}
void HermiteGalsCharacteristicEnclosure::validateMesh(const QuadtreeMesh& mesh) const{checked().interval.validateMesh(mesh);}
HermiteGalsCharacteristicEnclosure encloseHermiteGalsCharacteristic(const QuadtreeMesh& mesh,
    const HermiteGalsPhaseFluxIntervalSource& interval,Vector2D lower,Vector2D upper)
{
    interval.validateMesh(mesh);
    require(interval.previousAccepted().primary().geometryCells.size()==static_cast<std::size_t>(mesh.cellIndexNum()),
        "characteristic enclosure requires full old scalar support");
    const auto* field=interval.velocitySource().polynomialVelocity();
    require(field,"characteristic enclosure requires an owned polynomial source");
    try{
        const auto bounds=field->encloseInteriorRk2(interval.characteristicSource().domain(),lower,upper,interval.timeBegin(),interval.timeEnd());
        return HermiteGalsCharacteristicEnclosure(std::make_shared<const HermiteGalsCharacteristicEnclosure::Data>(interval,bounds));
    }catch(const std::exception& e){throw std::runtime_error("[Interface/HermiteGalsCharacteristicEnclosure previous="
        +std::to_string(interval.previousGeometryIdentity().value())+"] "+e.what());}
}

namespace
{
// Verify the union of CLOSED rectangles, including degenerate faces/points.
// Every x strip and its full y range is covered by endpoint comparisons;
// no rounded midpoint, area sum, or tolerance can conceal a small gap.
void closedPhaseCoverage(Vector2D lo,Vector2D hi,const std::vector<HermiteGalsCharacteristicPhaseCell>& cells)
{
    require(!cells.empty(),"phase coverage has no old native leaves");
    std::vector<double> xs{lo.x,hi.x};std::vector<const HermiteGalsCharacteristicPhaseCell*> ordered;
    for(const auto& c:cells){xs.push_back(c.lower.x);xs.push_back(c.upper.x);ordered.push_back(&c);}
    std::sort(xs.begin(),xs.end());xs.erase(std::unique(xs.begin(),xs.end()),xs.end());
    std::sort(ordered.begin(),ordered.end(),[](const auto* a,const auto* b){
        return a->lower.y!=b->lower.y?a->lower.y<b->lower.y:a->upper.y<b->upper.y;});
    const auto strip=[&](double left,double right){
        double reach=lo.y;bool covered=false;
        for(const auto* c:ordered)if(c->lower.x<=left && c->upper.x>=right){
            require(c->lower.y<=reach,"phase coverage has a gap in old scalar support");
            reach=std::max(reach,c->upper.y);covered=true;
        }
        require(covered && reach==hi.y,"phase coverage leaves old scalar support");
    };
    if(xs.size()==1)strip(xs.front(),xs.front());
    else for(std::size_t i=1;i<xs.size();++i)strip(xs[i-1],xs[i]);
}
struct ClosedPhaseBox
{
    std::vector<HermiteGalsCharacteristicPhaseCell> cells;
    HermiteGalsCertifiedPhase phase=HermiteGalsCertifiedPhase::Unresolved;
    long double minimumLiquidContent=1,maximumLiquidContent=0;
};
// Shared kernel after the caller validates the immutable source and mesh.
// An adaptive integral does not manufacture a separate interval/mesh per tile.
ClosedPhaseBox coverPhaseBox(const LevelSetHermiteGeometry::Snapshot& scalar,
    Vector2D lo,Vector2D hi,const HermiteGalsInventoryRegionSource* inventory)
{
    ClosedPhaseBox result;
    auto candidates=scalar.nativeCellCandidates(lo,hi);std::sort(candidates.begin(),candidates.end());
    require(std::adjacent_find(candidates.begin(),candidates.end())==candidates.end(),"duplicate native phase candidate");
    bool negative=true,positive=true;
    for(int cell:candidates){const auto bounds=scalar.nativeCellBounds(cell);
        const Vector2D a{std::max(lo.x,bounds[0].x),std::max(lo.y,bounds[0].y)},
            b{std::min(hi.x,bounds[1].x),std::min(hi.y,bounds[1].y)};
        if(a.x>b.x || a.y>b.y)continue;
        require(scalar.nativeCellActive(cell),"characteristic phase coverage intersects inactive scalar");
        const int sign=scalar.nativeQuarticBoxSign(cell,a,b);
        require(sign==-1 || sign==0 || sign==1,"invalid native phase certificate");
        long double liquid=1,gas=0;
        if(inventory){const auto& row=inventory->cells().at(cell);require(row.cell==cell,"foreign native inventory row");
            liquid=row.liquidWeight;gas=row.gasWeight;}
        require(std::isfinite(gas) && std::isfinite(liquid) && gas>=0 && gas<=liquid && liquid<=1,
            "characteristic phase content outside strict bounds");
        const long double lower=sign<0?liquid:gas,upper=sign>0?gas:liquid;
        result.cells.push_back({cell,a,b,static_cast<HermiteGalsCertifiedPhase>(sign),lower,upper});
        result.minimumLiquidContent=std::min(result.minimumLiquidContent,lower);result.maximumLiquidContent=std::max(result.maximumLiquidContent,upper);
        negative=negative && sign<0;positive=positive && sign>0;
    }
    closedPhaseCoverage(lo,hi,result.cells);
    result.phase=negative?HermiteGalsCertifiedPhase::Liquid:positive?HermiteGalsCertifiedPhase::Gas:HermiteGalsCertifiedPhase::Unresolved;
    return result;
}
ClosedPhaseBox coverDepartureBox(const LevelSetHermiteGeometry::Snapshot& scalar,
    const LAMP::Geometry::PolynomialRk2Bounds& box,const HermiteGalsInventoryRegionSource* inventory)
{return coverPhaseBox(scalar,box.departureLower,box.departureUpper,inventory);}
}
struct HermiteGalsCharacteristicPhaseCoverage::Data
{
    HermiteGalsCharacteristicEnclosure characteristic;
    std::optional<HermiteGalsInventoryRegionSource> inventory;
    ClosedPhaseBox coverage;
    Data(const QuadtreeMesh& mesh,const HermiteGalsCharacteristicEnclosure& source,
        std::optional<HermiteGalsInventoryRegionSource> weights):characteristic(source),inventory(std::move(weights))
    {
        characteristic.validateMesh(mesh);
        const auto& old=characteristic.interval().previousAccepted();const auto& scalar=old.geometry().snapshot();
        require(scalar.usesNativeLeaves() && scalar.reconstruction()==LevelSetHermiteGeometry::Reconstruction::EnrichedQuartic,
            "characteristic phase coverage requires native Q4");
        if(inventory){inventory->validateMesh(mesh);
            require(inventory->source().identity()==old.identity() && &inventory->source().geometry().snapshot()==&scalar,
                "characteristic phase coverage has foreign inventory/geometry source");}
        coverage=coverDepartureBox(scalar,characteristic.bounds(),inventory?&*inventory:nullptr);
    }
};
HermiteGalsCharacteristicPhaseCoverage::HermiteGalsCharacteristicPhaseCoverage(std::shared_ptr<const Data> data):data_(std::move(data)){}
const HermiteGalsCharacteristicPhaseCoverage::Data& HermiteGalsCharacteristicPhaseCoverage::checked() const
{require(bool(data_),"moved-from characteristic phase coverage");return *data_;}
const HermiteGalsCharacteristicEnclosure& HermiteGalsCharacteristicPhaseCoverage::characteristic() const{return checked().characteristic;}
const std::vector<HermiteGalsCharacteristicPhaseCell>& HermiteGalsCharacteristicPhaseCoverage::cells() const{return checked().coverage.cells;}
HermiteGalsCertifiedPhase HermiteGalsCharacteristicPhaseCoverage::phase() const{return checked().coverage.phase;}
long double HermiteGalsCharacteristicPhaseCoverage::minimumLiquidContent() const{return checked().coverage.minimumLiquidContent;}
long double HermiteGalsCharacteristicPhaseCoverage::maximumLiquidContent() const{return checked().coverage.maximumLiquidContent;}
const HermiteGalsInventoryRegionSource* HermiteGalsCharacteristicPhaseCoverage::inventorySource() const
{const auto& data=checked();return data.inventory?&*data.inventory:nullptr;}
void HermiteGalsCharacteristicPhaseCoverage::validateMesh(const QuadtreeMesh& mesh) const{checked().characteristic.validateMesh(mesh);}
HermiteGalsCharacteristicPhaseCoverage coverHermiteGalsCharacteristicPhase(const QuadtreeMesh& mesh,
    const HermiteGalsCharacteristicEnclosure& source)
{return HermiteGalsCharacteristicPhaseCoverage(std::make_shared<const HermiteGalsCharacteristicPhaseCoverage::Data>(mesh,source,std::nullopt));}
HermiteGalsCharacteristicPhaseCoverage coverHermiteGalsCharacteristicPhase(const QuadtreeMesh& mesh,
    const HermiteGalsCharacteristicEnclosure& source,const HermiteGalsInventoryRegionSource& inventory)
{return HermiteGalsCharacteristicPhaseCoverage(std::make_shared<const HermiteGalsCharacteristicPhaseCoverage::Data>(mesh,source,inventory));}

struct HermiteGalsBoundaryPhaseSource::Data
{
    HermiteGalsPhaseFluxIntervalSource extent;
    std::array<HermiteGalsInflowPhase,4> phases;
    Data(const HermiteGalsPhaseFluxIntervalSource& e,const std::array<HermiteGalsInflowPhase,4>& p):extent(e),phases(p){}
};
HermiteGalsBoundaryPhaseSource::HermiteGalsBoundaryPhaseSource(const QuadtreeMesh& mesh,
    const HermiteGalsPhaseFluxIntervalSource& extent,const std::array<HermiteGalsInflowPhase,4>& phases)
{
    extent.validateMesh(mesh);
    require(!extent.inflowSource() && !extent.boundaryPhaseSource(),"boundary history requires an unbound base interval");
    require(extent.velocitySource().polynomialVelocity(),"boundary history requires an owned polynomial source");
    for(const auto phase:phases)require(phase==HermiteGalsInflowPhase::Gas || phase==HermiteGalsInflowPhase::Liquid,"invalid boundary phase");
    data_=std::make_shared<const Data>(extent,phases);
}
const HermiteGalsBoundaryPhaseSource::Data& HermiteGalsBoundaryPhaseSource::checked() const
{require(bool(data_),"moved-from boundary phase source");return *data_;}
const HermiteGalsPhaseFluxIntervalSource& HermiteGalsBoundaryPhaseSource::baseExtent() const{return checked().extent;}
HermiteGalsInflowPhase HermiteGalsBoundaryPhaseSource::phaseFor(Boundary side) const
{
    const auto& d=checked();
    switch(side){case Boundary::west:return d.phases[0];case Boundary::south:return d.phases[1];
        case Boundary::east:return d.phases[2];case Boundary::north:return d.phases[3];
        default:throw std::invalid_argument("[Interface/HermiteGalsTemporalPhaseFlux] invalid boundary side");}
}
bool HermiteGalsBoundaryPhaseSource::sameSource(const HermiteGalsBoundaryPhaseSource& other) const noexcept
{return data_ && data_==other.data_;}
void HermiteGalsBoundaryPhaseSource::validateInterval(const QuadtreeMesh& mesh,const HermiteGalsPhaseFluxIntervalSource& interval) const
{
    const auto& base=checked().extent;base.validateMesh(mesh);interval.validateMesh(mesh);
    require(base.previousGeometryIdentity()==interval.previousGeometryIdentity(),"boundary history has different old source");
    require(base.velocitySource().sameSource(interval.velocitySource()),"boundary history has different velocity source");
    require(base.characteristicSource().sameSource(interval.characteristicSource()),"boundary history has different characteristic source");
    require(interval.timeBegin()>=base.timeBegin() && interval.timeEnd()<=base.timeEnd(),"boundary history does not cover interval");
}

struct HermiteGalsConstantInflowSource::Data
{
    HermiteGalsPhaseFluxIntervalSource extent;
    std::vector<HermiteGalsInflowFacePhase> entries;
    Data(const HermiteGalsPhaseFluxIntervalSource& e,const std::vector<HermiteGalsInflowFacePhase>& p):extent(e),entries(p){}
};
HermiteGalsConstantInflowSource::HermiteGalsConstantInflowSource(const QuadtreeMesh& mesh,
    const HermiteGalsPhaseFluxIntervalSource& extent,const std::vector<HermiteGalsInflowFacePhase>& entries)
{
    require(!extent.inflowSource() && !extent.boundaryPhaseSource(),"inflow history requires an unbound base interval");
    const auto velocity=constantTranslation(mesh,extent);std::size_t k=0;
    for(int row=0;row<mesh.vEdgeIndexNum()+mesh.hEdgeIndexNum();++row){const auto face=extent.face(mesh,row);
        if(!isInflow(face,velocity))continue;
        require(k<entries.size() && entries[k].nativeFace==row,"missing/duplicate/unordered inflow face entry");
        require(entries[k].phase==HermiteGalsInflowPhase::Gas || entries[k].phase==HermiteGalsInflowPhase::Liquid,"invalid inflow phase");++k;}
    require(k==entries.size(),"foreign/non-inflow face entry");
    data_=std::make_shared<const Data>(extent,entries);
}
const HermiteGalsConstantInflowSource::Data& HermiteGalsConstantInflowSource::checked() const
{require(bool(data_),"moved-from inflow source");return *data_;}
const HermiteGalsPhaseFluxIntervalSource& HermiteGalsConstantInflowSource::baseExtent() const{return checked().extent;}
const std::vector<HermiteGalsInflowFacePhase>& HermiteGalsConstantInflowSource::entries() const{return checked().entries;}
HermiteGalsInflowPhase HermiteGalsConstantInflowSource::phaseFor(int nativeFace) const
{
    const auto& entries=checked().entries;
    const auto found=std::lower_bound(entries.begin(),entries.end(),nativeFace,[](const auto& e,int f){return e.nativeFace<f;});
    require(found!=entries.end() && found->nativeFace==nativeFace,"missing prescribed inflow face");return found->phase;
}
bool HermiteGalsConstantInflowSource::sameSource(const HermiteGalsConstantInflowSource& other) const noexcept
{return data_ && data_==other.data_;}
void HermiteGalsConstantInflowSource::validateInterval(const QuadtreeMesh& mesh,const HermiteGalsPhaseFluxIntervalSource& interval) const
{
    const auto& base=checked().extent;base.validateMesh(mesh);interval.validateMesh(mesh);
    require(base.previousGeometryIdentity()==interval.previousGeometryIdentity(),"inflow history has different old source");
    require(base.velocitySource().sameSource(interval.velocitySource()),"inflow history has different velocity source");
    require(base.characteristicSource().sameSource(interval.characteristicSource()),"inflow history has different characteristic source");
    require(interval.timeBegin()>=base.timeBegin() && interval.timeEnd()<=base.timeEnd(),"inflow history does not cover interval");
}

struct HermiteGalsIntegratedFaceTransfer::Data
{
    HermiteGalsPhaseFluxIntervalSource interval;
    NativeLeafFace face;
    HermiteGalsFaceTransferOrientation orientation;
    HermiteGalsIntegratedFaceVolumes volumes;
    Data(const HermiteGalsPhaseFluxIntervalSource& i, NativeLeafFace f,
        HermiteGalsFaceTransferOrientation o, const HermiteGalsIntegratedFaceVolumes& v)
        : interval(i), face(f), orientation(o), volumes(v) {}
};
HermiteGalsIntegratedFaceTransfer::HermiteGalsIntegratedFaceTransfer(const QuadtreeMesh& mesh,
    const HermiteGalsPhaseFluxIntervalSource& interval, int nativeFace,
    HermiteGalsFaceTransferOrientation orientation, const HermiteGalsIntegratedFaceVolumes& volumes)
{
    const auto face = interval.face(mesh, nativeFace);
    require(orientation == (face.horizontal ? HermiteGalsFaceTransferOrientation::PositiveY
        : HermiteGalsFaceTransferOrientation::PositiveX), "incorrect native face orientation");
    validVolumes(volumes);
    if (face.area == 0)
        require(volumes.integratedTotalPositiveVolume == 0 && volumes.integratedTotalNegativeVolume == 0,
            "nonzero integrated volume on zero-measure face");
    data_ = std::make_shared<const Data>(interval, face, orientation, volumes);
}
const HermiteGalsIntegratedFaceTransfer::Data& HermiteGalsIntegratedFaceTransfer::checked() const
{ require(bool(data_), "moved-from integrated face transfer"); return *data_; }
const HermiteGalsPhaseFluxIntervalSource& HermiteGalsIntegratedFaceTransfer::interval() const { return checked().interval; }
const NativeLeafFace& HermiteGalsIntegratedFaceTransfer::face() const { return checked().face; }
HermiteGalsFaceTransferOrientation HermiteGalsIntegratedFaceTransfer::orientation() const { return checked().orientation; }
const HermiteGalsIntegratedFaceVolumes& HermiteGalsIntegratedFaceTransfer::integratedVolumes() const { return checked().volumes; }
long double HermiteGalsIntegratedFaceTransfer::signedTotalVolume() const
{ const auto& v = checked().volumes; return v.integratedTotalPositiveVolume - v.integratedTotalNegativeVolume; }
long double HermiteGalsIntegratedFaceTransfer::signedLiquidVolume() const
{ const auto& v = checked().volumes; return v.integratedLiquidPositiveVolume - v.integratedLiquidNegativeVolume; }
long double HermiteGalsIntegratedFaceTransfer::averageTotalFlux() const { return average(signedTotalVolume(), interval().dt()); }
long double HermiteGalsIntegratedFaceTransfer::averageLiquidFlux() const { return average(signedLiquidVolume(), interval().dt()); }
HermiteGalsIntegratedFaceTransfer HermiteGalsIntegratedFaceTransfer::appended(
    const QuadtreeMesh& mesh, const HermiteGalsIntegratedFaceTransfer& next) const
{
    const auto& a = checked(); const auto& b = next.checked();
    a.interval.validateMesh(mesh); b.interval.validateMesh(mesh);
    require(a.interval.previousGeometryIdentity() == b.interval.previousGeometryIdentity(), "different previous geometry source");
    require(a.interval.velocitySource().sameSource(b.interval.velocitySource()), "different velocity source");
    require(a.interval.characteristicSource().sameSource(b.interval.characteristicSource()), "different characteristic source");
    require(sameInflow(a.interval.inflowSource(),b.interval.inflowSource()),"different inflow history source");
    require(sameBoundary(a.interval.boundaryPhaseSource(),b.interval.boundaryPhaseSource()),"different boundary history source");
    require(a.face.nativeRow == b.face.nativeRow && a.orientation == b.orientation
        && a.face.negativeCell == b.face.negativeCell && a.face.positiveCell == b.face.positiveCell,
        "different native face identity/orientation");
    require(a.interval.timeEnd() == b.interval.timeBegin(), "intervals are not chronologically contiguous");
    const auto& av = a.volumes; const auto& bv = b.volumes;
    const HermiteGalsIntegratedFaceVolumes sum{
        av.integratedTotalPositiveVolume + bv.integratedTotalPositiveVolume,
        av.integratedTotalNegativeVolume + bv.integratedTotalNegativeVolume,
        av.integratedLiquidPositiveVolume + bv.integratedLiquidPositiveVolume,
        av.integratedLiquidNegativeVolume + bv.integratedLiquidNegativeVolume};
    validVolumes(sum);
    const auto combined=a.interval.inflowSource()
        ?HermiteGalsPhaseFluxIntervalSource(mesh,a.interval.previousAccepted(),a.interval.velocitySource(),a.interval.characteristicSource(),
            a.interval.timeBegin(),b.interval.timeEnd(),*a.interval.inflowSource())
        :a.interval.boundaryPhaseSource()
        ?HermiteGalsPhaseFluxIntervalSource(mesh,a.interval.previousAccepted(),a.interval.velocitySource(),a.interval.characteristicSource(),
            a.interval.timeBegin(),b.interval.timeEnd(),*a.interval.boundaryPhaseSource())
        :HermiteGalsPhaseFluxIntervalSource(mesh,a.interval.previousAccepted(),a.interval.velocitySource(),a.interval.characteristicSource(),
            a.interval.timeBegin(),b.interval.timeEnd());
    return {mesh, combined, a.face.nativeRow, a.orientation, sum};
}

namespace
{
using IntegralBounds=LAMP::Geometry::PolynomialFaceIntegralBounds;
constexpr std::size_t noIntegralNode=std::numeric_limits<std::size_t>::max();
struct IntegralTile
{
    double lower=0,upper=0,begin=0,end=0,phaseSpread=0;
    int depth=0;
    std::size_t parent=noIntegralNode,left=noIntegralNode,right=noIntegralNode,pending=0,unresolvedPhase=0;
    IntegralBounds bounds;
};
struct IntegralComputation {IntegralBounds bounds;HermiteGalsSpaceTimeIntegralReport report;};
IntegralComputation boundedPolynomialFace(const QuadtreeMesh& mesh,const HermiteGalsPhaseFluxIntervalSource& interval,
    int nativeFace,const HermiteGalsInventoryRegionSource* inventory,const HermiteGalsSpaceTimeIntegralConfig& config)
{
    require(std::isfinite(config.absoluteVolumeTolerance) && config.absoluteVolumeTolerance>0
        && config.maximumPanels>0 && config.maximumPanels<(noIntegralNode/2) && config.maximumDepth>0,"invalid space-time integral budget");
    const auto face=interval.face(mesh,nativeFace);const auto& old=interval.previousAccepted();const auto& scalar=old.geometry().snapshot();
    require(old.primary().geometryCells.size()==static_cast<std::size_t>(mesh.cellIndexNum()),"polynomial integral requires full old scalar support");
    require(scalar.usesNativeLeaves() && scalar.reconstruction()==LevelSetHermiteGeometry::Reconstruction::EnrichedQuartic,"polynomial integral requires native Q4");
    const auto* field=interval.velocitySource().polynomialVelocity();require(field,"polynomial integral requires owned polynomial velocity");
    const auto* boundary=interval.boundaryPhaseSource();if(boundary)boundary->validateInterval(mesh,interval);
    if(inventory){inventory->validateMesh(mesh);
        require(inventory->source().identity()==old.identity() && &inventory->source().geometry().snapshot()==&scalar,"polynomial integral has foreign inventory/geometry source");}
    const auto& edge=face.horizontal?mesh.hEdges(face.edgeLabel):mesh.vEdges(face.edgeLabel);
    auto a=mesh.verts(edge.edgeToVert[0]).coord(),b=mesh.verts(edge.edgeToVert[1]).coord();
    if(face.horizontal?a.x>b.x:a.y>b.y)std::swap(a,b);
    const double lo=face.horizontal?a.x:a.y,hi=face.horizontal?b.x:b.y;
    const auto& domain=interval.characteristicSource().domain();const bool radial=interval.measure()==HermiteGalsMeasure::AxisymmetricVolume;
    // Validate even a zero-flux query through the original domain/FP rules.
    field->encloseInteriorRk2(domain,a,a,field->timeBegin(),field->timeBegin());
    if(radial && a.x==0 && b.x==0){const auto axis=field->enclose(a,b,interval.timeBegin(),interval.timeEnd());
        require(axis.value[0].lower==0 && axis.value[0].upper==0 && axis.derivative[1].lower==0 && axis.derivative[1].upper==0,
            "polynomial integral requires exact axis normal velocity/derivative");}
    std::optional<LAMP::Geometry::PolynomialRk2Reachability> reachable;
    if(boundary && boundary->baseExtent().timeBegin()==field->timeBegin()){
        try{reachable=field->encloseRk2Reachability(domain,interval.timeEnd());}
        catch(const LAMP::Geometry::PolynomialRk2EnclosureUnresolved&){/* original branch proofs remain available */}
    }
    IntegralComputation result;std::string lastBranchRefusal;
    const auto position=[&](double s){return face.horizontal?Vector2D{s,a.y}:Vector2D{a.x,s};};
    const auto evaluate=[&](double lower,double upper,double begin,double end,int depth,std::size_t parent){
        IntegralTile tile;tile.lower=lower;tile.upper=upper;tile.begin=begin;tile.end=end;tile.depth=depth;tile.parent=parent;
        ++result.report.evaluatedPanels;result.report.maximumDepthReached=std::max(result.report.maximumDepthReached,depth);
        const auto p=position(lower),q=position(upper);tile.bounds=field->encloseNormalFaceIntegral(p,q,begin,end,radial);
        if(tile.bounds.volumes[0].upper==0 && tile.bounds.volumes[1].upper==0)return tile;
        if(boundary && (face.negativeCell<0 || face.positiveCell<0)){
            try{
                const auto side=field->encloseInitialBoundaryHit(domain,p,q,begin,end);
                const long double content=boundary->phaseFor(side)==HermiteGalsInflowPhase::Liquid?1:0;
                tile.bounds=tile.bounds.withContentRange(content,content);++result.report.initialBoundaryQueries;return tile;
            }catch(const LAMP::Geometry::PolynomialRk2EnclosureUnresolved&){
                // Lack of this certificate is NOT a boundary/old phase
                // fallback. The original complete no-hit proof remains next.
            }
        }
        if(reachable){
            const auto sources=reachable->cover(p,q);
            const auto phase=coverPhaseBox(scalar,sources.lower,sources.upper,inventory);
            auto minimum=phase.minimumLiquidContent,maximum=phase.maximumLiquidContent;
            const std::array<Boundary,4> sides{Boundary::west,Boundary::south,Boundary::east,Boundary::north};
            for(int k=0;k<4;++k)if(sources.boundary[k]){
                const long double content=boundary->phaseFor(sides[k])==HermiteGalsInflowPhase::Liquid?1:0;
                minimum=std::min(minimum,content);maximum=std::max(maximum,content);
            }
            // Exact equality of ALL source contents; no geometric-gas or
            // tolerance fallback, and no assertion that a B/J was computed.
            if(minimum==maximum){tile.bounds=tile.bounds.withContentRange(minimum,maximum);
                ++result.report.uniformReachableQueries;return tile;}
        }
        try{
            const auto characteristic=field->encloseInteriorRk2TightTime(domain,p,q,begin,end);
            const auto phase=coverDepartureBox(scalar,characteristic,inventory);++result.report.phaseQueries;
            tile.bounds=tile.bounds.withContentRange(phase.minimumLiquidContent,phase.maximumLiquidContent);
            tile.unresolvedPhase=phase.phase==HermiteGalsCertifiedPhase::Unresolved?1:0;
            if(tile.unresolvedPhase){
                const auto affine=field->encloseAffineFaceRk2(domain,p,q,begin,end);
                std::vector<LAMP::Geometry::PolynomialAffineCoordinate> forms;forms.reserve(phase.cells.size());
                long double wl0=1,wl1=0,wg0=1,wg1=0;
                for(const auto& c:phase.cells){
                    forms.push_back(scalar.nativeQuarticAffineEnclosure(c.cell,affine.departure));
                    const long double wl=inventory?inventory->cells().at(c.cell).liquidWeight:1,
                        wg=inventory?inventory->cells().at(c.cell).gasWeight:0;
                    wl0=std::min(wl0,wl);wl1=std::max(wl1,wl);wg0=std::min(wg0,wg);wg1=std::max(wg1,wg);
                }
                tile.bounds=tile.bounds.intersected(field->encloseNormalFaceIntegralWithAffinePhase(p,q,begin,end,radial,forms,wl0,wl1,wg0,wg1));
                // A directional work heuristic only. These values never
                // classify phase, bound an integral, or certify an event.
                const auto lo=characteristic.departureLower,hi=characteristic.departureUpper;
                double smallest=std::numeric_limits<double>::infinity(),largest=-smallest;
                for(int i=0;i<3;++i)for(int j=0;j<3;++j){
                    const double value=scalar.value({i==2?hi.x:lo.x+(hi.x-lo.x)*(i*.5),j==2?hi.y:lo.y+(hi.y-lo.y)*(j*.5)});
                    smallest=std::min(smallest,value);largest=std::max(largest,value);}
                tile.phaseSpread=largest-smallest;
            }
        }catch(const LAMP::Geometry::PolynomialRk2EnclosureUnresolved& error){
            // Conditional [0,1] content guides refinement only. A pending
            // branch can never escape into the successful result.
            tile.pending=1;++result.report.branchRefusals;lastBranchRefusal=error.what();
        }
        return tile;
    };
    std::vector<IntegralTile> nodes;nodes.reserve(2*std::min<std::size_t>(config.maximumPanels,1024)-1);
    nodes.push_back(evaluate(lo,hi,interval.timeBegin(),interval.timeEnd(),0,noIntegralNode));
    struct Entry{bool pending;double error;std::size_t node;};
    const auto earlier=[](const Entry& x,const Entry& y){if(x.pending!=y.pending)return x.pending<y.pending;
        if(x.error!=y.error)return x.error<y.error;return x.node>y.node;};
    std::priority_queue<Entry,std::vector<Entry>,decltype(earlier)> queue(earlier);
    const auto enqueue=[&](std::size_t n){queue.push({nodes[n].pending!=0,nodes[n].bounds.errorBound(),n});};
    enqueue(0);std::size_t leaves=1;
    const auto fail=[&](const char* reason){std::ostringstream message;message.precision(17);
        message<<"[Interface/HermiteGalsBoundedFaceIntegral face="<<nativeFace<<"] "<<reason
            <<"; error="<<nodes.front().bounds.errorBound()<<" budget="<<config.absoluteVolumeTolerance
            <<" panels="<<leaves<<" pending="<<nodes.front().pending;
        if(!lastBranchRefusal.empty())message<<"; last branch refusal: "<<lastBranchRefusal;
        throw std::runtime_error(message.str());};
    while(nodes.front().pending || nodes.front().bounds.errorBound()>config.absoluteVolumeTolerance){
        if(leaves>=config.maximumPanels)fail("space-time panel budget exhausted");
        while(!queue.empty() && nodes[queue.top().node].left!=noIntegralNode)queue.pop();
        if(queue.empty())fail("space-time refinement lost its leaf queue");
        const auto index=queue.top().node;queue.pop();const auto parent=nodes[index];
        if(parent.depth>=config.maximumDepth)fail("space-time depth budget exhausted");
        struct Split {IntegralTile left,right;IntegralBounds bounds;std::size_t pending=0;};
        const auto split=[&](bool time,double middle){
            auto left=evaluate(parent.lower,time?parent.upper:middle,parent.begin,time?middle:parent.end,parent.depth+1,index);
            auto right=evaluate(time?parent.lower:middle,parent.upper,time?middle:parent.begin,parent.end,parent.depth+1,index);
            const auto bounds=parent.bounds.intersected(left.bounds.added(right.bounds));const auto pending=left.pending+right.pending;
            return Split{std::move(left),std::move(right),bounds,pending};
        };
        std::optional<Split> spatial,temporal;
        if(const auto middle=LAMP::Geometry::polynomialFaceMidpoint(parent.lower,parent.upper))spatial=split(false,*middle);
        if(const auto middle=LAMP::Geometry::polynomialFaceMidpoint(parent.begin,parent.end))temporal=split(true,*middle);
        if(!spatial && !temporal)fail("space-time subdivision has no represented interior point");
        // Both complete partitions enclose the SAME parent integral. Retain
        // their common certificate even though only one partition becomes
        // the refinement tree. No extra panels/queries or sampled proof.
        const auto common=spatial && temporal?spatial->bounds.intersected(temporal->bounds):
            spatial?spatial->bounds:temporal->bounds;
        bool useTime=!spatial;
        if(spatial && temporal){
            if(spatial->pending!=temporal->pending)useTime=temporal->pending<spatial->pending;
            else{const double sError=spatial->bounds.errorBound(),tError=temporal->bounds.errorBound();
                // A material immediate reduction wins. Otherwise prefer
                // the direction reducing phase variation; tiny gross-speed
                // improvements must not starve the transverse event search.
                // This affects work selection only, never certified bounds.
                const double materialReduction=.125*parent.bounds.errorBound();
                if(std::abs(tError-sError)>materialReduction)useTime=tError<sError;
                else{
                    const auto spread=[](const Split& cut){return cut.left.phaseSpread*cut.left.bounds.errorBound()
                        +cut.right.phaseSpread*cut.right.bounds.errorBound();};
                    const double sSpread=spread(*spatial),tSpread=spread(*temporal);
                    useTime=std::isfinite(sSpread) && std::isfinite(tSpread) && sSpread!=tSpread?tSpread<sSpread:
                        (parent.end-parent.begin)/interval.dt()>=(parent.upper-parent.lower)/(hi-lo);
                }}
        }
        auto selected=useTime?std::move(*temporal):std::move(*spatial);
        const auto left=nodes.size();nodes.push_back(std::move(selected.left));const auto right=nodes.size();nodes.push_back(std::move(selected.right));
        nodes[index].left=left;nodes[index].right=right;nodes[index].bounds=common;nodes[index].pending=selected.pending;
        nodes[index].unresolvedPhase=nodes[left].unresolvedPhase+nodes[right].unresolvedPhase;
        for(auto ancestor=nodes[index].parent;ancestor!=noIntegralNode;ancestor=nodes[ancestor].parent){auto& node=nodes[ancestor];
            node.bounds=node.bounds.intersected(nodes[node.left].bounds.added(nodes[node.right].bounds));
            node.pending=nodes[node.left].pending+nodes[node.right].pending;
            node.unresolvedPhase=nodes[node.left].unresolvedPhase+nodes[node.right].unresolvedPhase;}
        ++leaves;enqueue(left);enqueue(right);
    }
    result.bounds=nodes.front().bounds;result.report.absoluteErrorBound=result.bounds.errorBound();
    result.report.acceptedPanels=leaves;result.report.unresolvedPhasePanels=nodes.front().unresolvedPhase;return result;
}
}
struct HermiteGalsBoundedFaceIntegral::Data
{
    HermiteGalsIntegratedFaceTransfer transfer;
    HermiteGalsSpaceTimeIntegralConfig config;
    IntegralBounds bounds;
    HermiteGalsSpaceTimeIntegralReport report;
    std::optional<HermiteGalsInventoryRegionSource> inventory;
    Data(const QuadtreeMesh& mesh,const HermiteGalsPhaseFluxIntervalSource& interval,int nativeFace,
        const HermiteGalsSpaceTimeIntegralConfig& controls,std::optional<HermiteGalsInventoryRegionSource> weights,IntegralComputation computed)
        :transfer([&]{const auto v=computed.bounds.representative();const auto face=interval.face(mesh,nativeFace);
            return HermiteGalsIntegratedFaceTransfer(mesh,interval,nativeFace,face.horizontal?HermiteGalsFaceTransferOrientation::PositiveY:HermiteGalsFaceTransferOrientation::PositiveX,
                {v[0],v[1],v[2],v[3]});}()),config(controls),bounds(computed.bounds),report(computed.report),inventory(std::move(weights)){}
};
HermiteGalsBoundedFaceIntegral::HermiteGalsBoundedFaceIntegral(std::shared_ptr<const Data> data):data_(std::move(data)){}
const HermiteGalsBoundedFaceIntegral::Data& HermiteGalsBoundedFaceIntegral::checked() const
{require(bool(data_),"moved-from bounded face integral");return *data_;}
const HermiteGalsIntegratedFaceTransfer& HermiteGalsBoundedFaceIntegral::transfer() const{return checked().transfer;}
const HermiteGalsSpaceTimeIntegralConfig& HermiteGalsBoundedFaceIntegral::config() const{return checked().config;}
const LAMP::Geometry::PolynomialFaceIntegralBounds& HermiteGalsBoundedFaceIntegral::bounds() const{return checked().bounds;}
const HermiteGalsSpaceTimeIntegralReport& HermiteGalsBoundedFaceIntegral::report() const{return checked().report;}
const HermiteGalsInventoryRegionSource* HermiteGalsBoundedFaceIntegral::inventorySource() const
{const auto& d=checked();return d.inventory?&*d.inventory:nullptr;}
void HermiteGalsBoundedFaceIntegral::validateMesh(const QuadtreeMesh& mesh) const{checked().transfer.interval().validateMesh(mesh);}
HermiteGalsBoundedFaceIntegral integrateHermiteGalsPolynomialFace(const QuadtreeMesh& mesh,
    const HermiteGalsPhaseFluxIntervalSource& interval,int face,const HermiteGalsSpaceTimeIntegralConfig& config)
{
    auto computed=boundedPolynomialFace(mesh,interval,face,nullptr,config);
    return HermiteGalsBoundedFaceIntegral(std::make_shared<const HermiteGalsBoundedFaceIntegral::Data>(mesh,interval,face,config,std::nullopt,std::move(computed)));
}
HermiteGalsBoundedFaceIntegral integrateHermiteGalsPolynomialFace(const QuadtreeMesh& mesh,
    const HermiteGalsPhaseFluxIntervalSource& interval,int face,const HermiteGalsInventoryRegionSource& inventory,const HermiteGalsSpaceTimeIntegralConfig& config)
{
    auto computed=boundedPolynomialFace(mesh,interval,face,&inventory,config);
    return HermiteGalsBoundedFaceIntegral(std::make_shared<const HermiteGalsBoundedFaceIntegral::Data>(mesh,interval,face,config,inventory,std::move(computed)));
}

HermiteGalsIntegratedFaceTransfer integrateHermiteGalsPlanarFace(const QuadtreeMesh& mesh,
    const HermiteGalsPhaseFluxIntervalSource& interval, int nativeFace)
{
    const auto face = interval.face(mesh, nativeFace);
    const auto velocity = interval.velocitySource().constantValue();
    require(velocity.has_value(), "planar integral requires an owned constant velocity source");
    const auto& domain = interval.characteristicSource().domain();
    require(!domain.activeSegment, "planar integral cannot certify an opaque coverage callback");
    require(!domain.radialAxisSymmetry || velocity->x == 0, "RZ planar integral requires axial translation");
    const auto& edge = face.horizontal ? mesh.hEdges(face.edgeLabel) : mesh.vEdges(face.edgeLabel);
    auto a = mesh.verts(edge.edgeToVert[0]).coord(), b = mesh.verts(edge.edgeToVert[1]).coord();
    if (face.horizontal ? a.x > b.x : a.y > b.y) std::swap(a, b);
    const auto measure = interval.previousAccepted().geometry().snapshot().integrateNativeAffineTranslatedFace(
        a, b, *velocity, interval.previousAccepted().timeLayer().time(), interval.timeBegin(), interval.timeEnd());
    const long double speed = face.horizontal ? velocity->y : velocity->x;
    const long double total = std::abs(speed) * measure.total, liquid = total * measure.liquidFraction;
    HermiteGalsIntegratedFaceVolumes v;
    if (speed >= 0) { v.integratedTotalPositiveVolume = total; v.integratedLiquidPositiveVolume = liquid; }
    else { v.integratedTotalNegativeVolume = total; v.integratedLiquidNegativeVolume = liquid; }
    return {mesh, interval, nativeFace, face.horizontal ? HermiteGalsFaceTransferOrientation::PositiveY
        : HermiteGalsFaceTransferOrientation::PositiveX, v};
}

namespace {
struct SweepQuery {
    HermiteGalsPhaseRegionMeasure geometry;
    long double liquid=0,gas=0,l1Correction=0;
};
SweepQuery sweepQuery(HermiteGalsPhaseRegionMeasure query){
    const auto liquid=query.liquid,gas=query.gas;return {std::move(query),liquid,gas,0};
}
SweepQuery sweepQuery(HermiteGalsInventoryRegionMeasure query){
    return {std::move(query.geometry),query.liquid,query.gas,query.l1Correction};
}
bool inventoryWeighted(const HermiteGalsPhaseRegionSource&){return false;}
bool inventoryWeighted(const HermiteGalsInventoryRegionSource&){return true;}
template<class Regions>
HermiteGalsSweptFaceResult sweptFace(const QuadtreeMesh& mesh,
    const HermiteGalsPhaseFluxIntervalSource& interval,int nativeFace,
    const Regions& regions,const HermiteGalsSweptFaceConfig& config)
{
    const auto face=interval.face(mesh,nativeFace);regions.validateMesh(mesh);
    const auto& old=interval.previousAccepted();
    require(regions.source().identity()==old.identity() && &regions.source().geometry().snapshot()==&old.geometry().snapshot(),
        "swept integral requires the identical old region source");
    for(double value:{config.relativeCoordinateTolerance,config.relativeMomentEstimateTolerance})
        require(std::isfinite(value) && value>0 && value<1,"invalid swept integral tolerance");
    const auto velocity=interval.velocitySource().constantValue();
    require(velocity.has_value(),"swept integral requires an owned constant velocity source");
    const auto& domain=interval.characteristicSource().domain();
    require(!domain.activeSegment,"swept integral cannot certify an opaque coverage callback");
    const bool radial=interval.measure()==HermiteGalsMeasure::AxisymmetricVolume;
    require(!radial || velocity->x==0,"RZ swept integral requires axial translation");
    const double speed=face.horizontal?velocity->y:velocity->x,tangent=face.horizontal?velocity->x:velocity->y;
    const auto orientation=face.horizontal?HermiteGalsFaceTransferOrientation::PositiveY:HermiteGalsFaceTransferOrientation::PositiveX;
    HermiteGalsIntegratedFaceVolumes volumes;HermiteGalsSweptFaceReport report;
    report.inventoryWeighted=inventoryWeighted(regions);
    if(speed==0)return {{mesh,interval,nativeFace,orientation,volumes},report};
    require(tangent==0,"swept integral requires zero tangential velocity");
    const auto& edge=face.horizontal?mesh.hEdges(face.edgeLabel):mesh.vEdges(face.edgeLabel);
    auto a=mesh.verts(edge.edgeToVert[0]).coord(),b=mesh.verts(edge.edgeToVert[1]).coord();
    if(face.horizontal?a.x>b.x:a.y>b.y)std::swap(a,b);
    const auto box=old.geometry().snapshot().nativeNormalTranslationSweep(a,b,speed,
        old.timeLayer().time(),interval.timeBegin(),interval.timeEnd(),radial);
    report.coordinateVolumeBand=box.coordinateVolumeBand;
    require(box.coordinateVolumeBand<=config.relativeCoordinateTolerance*box.totalVolume,"swept coordinate volume band exceeds tolerance");
    const auto inner=sweepQuery(regions.measure(mesh,box.innerLower,box.innerUpper));
    const auto outer=box.coordinateVolumeBand==0?inner:sweepQuery(regions.measure(mesh,box.outerLower,box.outerUpper));
    const std::size_t rectangleQueries=box.coordinateVolumeBand==0?1:2;
    report.rectangleQueries=rectangleQueries;
    report.integratedArcs=inner.geometry.integratedArcs+(rectangleQueries==2?outer.geometry.integratedArcs:0);
    report.rootRefinements=inner.geometry.rootRefinements+(rectangleQueries==2?outer.geometry.rootRefinements:0);
    report.innerLiquid=inner.liquid;report.outerLiquid=outer.liquid;report.innerGas=inner.gas;report.outerGas=outer.gas;
    const auto momentErrorEstimate=.5L*inner.geometry.momentErrorEstimate+.5L*outer.geometry.momentErrorEstimate;
    report.momentErrorEstimate=momentErrorEstimate;
    require(std::isfinite(momentErrorEstimate) && momentErrorEstimate>=0 &&
        momentErrorEstimate<=config.relativeMomentEstimateTolerance*box.totalVolume,"swept moment estimate exceeds tolerance");
    // The geometric inner/outer containment is exact; their floating liquid
    // values are estimates, NOT lower/upper numerical volume bounds. Average
    // their direct smaller phase, preserving pure total/liquid identities.
    const auto negative=.5L*inner.liquid+.5L*outer.liquid,positive=.5L*inner.gas+.5L*outer.gas;
    const auto liquid=negative<=positive?negative:box.totalVolume-positive;
    require(std::isfinite(liquid) && liquid>=0 && liquid<=box.totalVolume,"swept liquid outside strict total volume");
    const auto geometricNegative=.5L*inner.geometry.liquid+.5L*outer.geometry.liquid;
    const auto geometricPositive=.5L*inner.geometry.gas+.5L*outer.geometry.gas;
    report.geometricLiquidVolume=geometricNegative<=geometricPositive?geometricNegative:box.totalVolume-geometricPositive;
    report.signedInventoryCorrection=(speed>0?1:-1)*(liquid-report.geometricLiquidVolume);
    report.l1InventoryCorrection=.5L*inner.l1Correction+.5L*outer.l1Correction;
    if(speed>0){volumes.integratedTotalPositiveVolume=box.totalVolume;volumes.integratedLiquidPositiveVolume=liquid;}
    else{volumes.integratedTotalNegativeVolume=box.totalVolume;volumes.integratedLiquidNegativeVolume=liquid;}
    return {{mesh,interval,nativeFace,orientation,volumes},report};
}

template<class Regions>
HermiteGalsTranslationFaceResult translationFace(const QuadtreeMesh& mesh,
    const HermiteGalsPhaseFluxIntervalSource& interval,int nativeFace,
    const Regions& regions,const HermiteGalsSweptFaceConfig& config)
{
    const auto velocity=constantTranslation(mesh,interval);const auto face=interval.face(mesh,nativeFace);
    if(!isInflow(face,velocity)){
        auto swept=sweptFace(mesh,interval,nativeFace,regions,config);
        return {std::move(swept.transfer),false,swept.report,inventoryWeighted(regions)};
    }
    regions.validateMesh(mesh);const auto& old=interval.previousAccepted();
    require(regions.source().identity()==old.identity() && &regions.source().geometry().snapshot()==&old.geometry().snapshot(),
        "translation requires the identical old region source");
    for(double value:{config.relativeCoordinateTolerance,config.relativeMomentEstimateTolerance})
        require(std::isfinite(value) && value>0 && value<1,"invalid translation integral tolerance");
    const auto* history=interval.inflowSource();require(history,"missing explicit inflow history");history->validateInterval(mesh,interval);
    const auto& edge=face.horizontal?mesh.hEdges(face.edgeLabel):mesh.vEdges(face.edgeLabel);
    auto a=mesh.verts(edge.edgeToVert[0]).coord(),b=mesh.verts(edge.edgeToVert[1]).coord();
    if(face.horizontal?a.x>b.x:a.y>b.y)std::swap(a,b);
    const auto speed=face.horizontal?velocity.y:velocity.x;
    const auto total=LevelSetHermiteGeometry::constantNormalFaceTransferVolume(a,b,speed,interval.timeBegin(),interval.timeEnd(),
        interval.measure()==HermiteGalsMeasure::AxisymmetricVolume);
    const auto liquid=history->phaseFor(nativeFace)==HermiteGalsInflowPhase::Liquid?total:0;
    HermiteGalsIntegratedFaceVolumes volumes;
    if(speed>0){volumes.integratedTotalPositiveVolume=total;volumes.integratedLiquidPositiveVolume=liquid;}
    else{volumes.integratedTotalNegativeVolume=total;volumes.integratedLiquidNegativeVolume=liquid;}
    return {{mesh,interval,nativeFace,face.horizontal?HermiteGalsFaceTransferOrientation::PositiveY:HermiteGalsFaceTransferOrientation::PositiveX,volumes},true,{},inventoryWeighted(regions)};
}
}
HermiteGalsSweptFaceResult integrateHermiteGalsSweptFace(const QuadtreeMesh& mesh,
    const HermiteGalsPhaseFluxIntervalSource& interval,int face,const HermiteGalsPhaseRegionSource& regions,const HermiteGalsSweptFaceConfig& config)
{return sweptFace(mesh,interval,face,regions,config);}
HermiteGalsSweptFaceResult integrateHermiteGalsSweptFace(const QuadtreeMesh& mesh,
    const HermiteGalsPhaseFluxIntervalSource& interval,int face,const HermiteGalsInventoryRegionSource& regions,const HermiteGalsSweptFaceConfig& config)
{return sweptFace(mesh,interval,face,regions,config);}
HermiteGalsTranslationFaceResult integrateHermiteGalsTranslationFace(const QuadtreeMesh& mesh,
    const HermiteGalsPhaseFluxIntervalSource& interval,int face,const HermiteGalsPhaseRegionSource& regions,const HermiteGalsSweptFaceConfig& config)
{return translationFace(mesh,interval,face,regions,config);}
HermiteGalsTranslationFaceResult integrateHermiteGalsTranslationFace(const QuadtreeMesh& mesh,
    const HermiteGalsPhaseFluxIntervalSource& interval,int face,const HermiteGalsInventoryRegionSource& regions,const HermiteGalsSweptFaceConfig& config)
{return translationFace(mesh,interval,face,regions,config);}

struct HermiteGalsIntegratedInventoryResult::Data
{
    explicit Data(const HermiteGalsPhaseFluxIntervalSource& source):interval(source){}
    HermiteGalsPhaseFluxIntervalSource interval;
    std::vector<HermiteGalsIntegratedFaceTransfer> raw,final;
    std::vector<HermiteGalsCellVolume> volumes;
    HermiteGalsInventoryAssemblyReport report;
};
HermiteGalsIntegratedInventoryResult::HermiteGalsIntegratedInventoryResult(std::shared_ptr<const Data> data):data_(std::move(data)){}
const HermiteGalsIntegratedInventoryResult::Data& HermiteGalsIntegratedInventoryResult::checked() const
{require(bool(data_),"moved-from integrated inventory result");return *data_;}
const HermiteGalsPhaseFluxIntervalSource& HermiteGalsIntegratedInventoryResult::interval() const {return checked().interval;}
const std::vector<HermiteGalsIntegratedFaceTransfer>& HermiteGalsIntegratedInventoryResult::rawFaces() const {return checked().raw;}
const std::vector<HermiteGalsIntegratedFaceTransfer>& HermiteGalsIntegratedInventoryResult::faces() const {return checked().final;}
const std::vector<HermiteGalsCellVolume>& HermiteGalsIntegratedInventoryResult::phaseVolumes() const {return checked().volumes;}
const HermiteGalsInventoryAssemblyReport& HermiteGalsIntegratedInventoryResult::report() const {return checked().report;}
void HermiteGalsIntegratedInventoryResult::validateMesh(const QuadtreeMesh& mesh) const {checked().interval.validateMesh(mesh);}

HermiteGalsIntegratedInventoryResult assembleHermiteGalsIntegratedInventory(const QuadtreeMesh& mesh,
    const HermiteGalsPhaseFluxIntervalSource& interval,const std::vector<HermiteGalsIntegratedFaceTransfer>& input,
    const HermiteGalsIntegratedInventoryConfig& config)
{return HermiteGalsIntegratedInventoryResult::assemble(mesh,interval,input,config,nullptr);}

HermiteGalsIntegratedInventoryResult HermiteGalsIntegratedInventoryResult::assemble(const QuadtreeMesh& mesh,
    const HermiteGalsPhaseFluxIntervalSource& interval,const std::vector<HermiteGalsIntegratedFaceTransfer>& input,
    const HermiteGalsIntegratedInventoryConfig& config,const std::vector<HermiteGalsBoundedFaceIntegral>* certificates)
{
    interval.validateMesh(mesh);
    const auto& previous=interval.previousAccepted();const auto& p=previous.primary();
    require(interval.timeBegin()==previous.timeLayer().time(),"inventory interval must begin at the old accepted time");
    for(double tolerance:{config.relativeDivergenceTolerance,config.relativeConservationTolerance})
        require(std::isfinite(tolerance) && tolerance>0 && tolerance<1,"invalid integrated inventory tolerance");
    FaceSegmentGraph graph;graph.build(mesh,[&](const QuadtreeMesh&,const Cell& cell){
        return std::binary_search(p.geometryCells.begin(),p.geometryCells.end(),cell.index());});
    const auto domain=NativeLeafDomain::build(mesh,graph,p.measure==HermiteGalsMeasure::AxisymmetricVolume);
    const Detail::HermiteInventoryAssembly assembly(static_cast<std::size_t>(mesh.cellIndexNum()),
        domain.cellIndexByRow(),domain.cellVolumes(),p.phaseVolumes,previous.identity());
    std::vector<Detail::InventoryFaceBounds> bounds;
    if(certificates){require(certificates->size()==input.size(),"certified assembly face count mismatch");bounds.reserve(input.size());}
    std::vector<Detail::InventoryFace> working;working.reserve(input.size());std::size_t index=0;
    for(const auto face:domain.faces(mesh)){
        if(face.negativeCell<0 && face.positiveCell<0)continue;
        require(index<input.size() && input[index].face().nativeRow==face.nativeRow,
            "missing/duplicate/unordered integrated native face");
        const auto& transfer=input[index++];
        require(transfer.interval().sameInterval(interval),"integrated face has a different interval source");
        require(transfer.face().negativeCell==face.negativeCell && transfer.face().positiveCell==face.positiveCell,
            "integrated native face owner mismatch");
        const auto& v=transfer.integratedVolumes();
        if(certificates){const auto& certificate=(*certificates)[index-1];
            require(certificate.transfer().face().nativeRow==face.nativeRow && certificate.transfer().interval().sameInterval(interval),
                "certified assembly has foreign face/interval");
            const auto& b=certificate.bounds().volumes;
            bounds.push_back({b[2].lower,std::min<long double>(b[2].upper,v.integratedTotalPositiveVolume),
                b[3].lower,std::min<long double>(b[3].upper,v.integratedTotalNegativeVolume)});
        }
        working.push_back({face.negativeCell,face.positiveCell,v.integratedTotalPositiveVolume,v.integratedTotalNegativeVolume,
            v.integratedLiquidPositiveVolume,v.integratedLiquidNegativeVolume,v.integratedLiquidPositiveVolume,v.integratedLiquidNegativeVolume});
    }
    require(index==input.size(),"foreign/inactive integrated native face");
    auto output=std::make_shared<HermiteGalsIntegratedInventoryResult::Data>(interval);
    output->report=assembly.finish(working,1,config.relativeDivergenceTolerance,config.relativeConservationTolerance,output->volumes,
        certificates?&bounds:nullptr);
    output->raw=input;output->final.reserve(input.size());
    for(std::size_t k=0;k<input.size();++k){
        auto transfer=input[k];auto v=transfer.integratedVolumes();
        const auto& final=working[k];
        if(final.liquidPositive!=v.integratedLiquidPositiveVolume || final.liquidNegative!=v.integratedLiquidNegativeVolume){
            v.integratedLiquidPositiveVolume=final.liquidPositive;v.integratedLiquidNegativeVolume=final.liquidNegative;
            validVolumes(v);
            // The interval mesh and ALL face ownership were validated once
            // above. Private construction avoids N complete mesh scans while
            // publishing one batch; no unchecked public factory is exposed.
            transfer.data_=std::make_shared<const HermiteGalsIntegratedFaceTransfer::Data>(
                interval,transfer.face(),transfer.orientation(),v);
        }
        output->final.push_back(std::move(transfer));
    }
    return HermiteGalsIntegratedInventoryResult(std::move(output));
}

struct HermiteGalsCertifiedInventoryResult::Data
{
    std::vector<HermiteGalsBoundedFaceIntegral> certificates;
    HermiteGalsIntegratedInventoryResult inventory;
    HermiteGalsCertifiedInventoryReport report;
    Data(const std::vector<HermiteGalsBoundedFaceIntegral>& c,HermiteGalsIntegratedInventoryResult i,
        HermiteGalsCertifiedInventoryReport r):certificates(c),inventory(std::move(i)),report(r){}
};
HermiteGalsCertifiedInventoryResult::HermiteGalsCertifiedInventoryResult(std::shared_ptr<const Data> data):data_(std::move(data)){}
const HermiteGalsCertifiedInventoryResult::Data& HermiteGalsCertifiedInventoryResult::checked() const
{require(bool(data_),"moved-from certified inventory result");return *data_;}
const HermiteGalsIntegratedInventoryResult& HermiteGalsCertifiedInventoryResult::inventory() const{return checked().inventory;}
const std::vector<HermiteGalsBoundedFaceIntegral>& HermiteGalsCertifiedInventoryResult::certificates() const{return checked().certificates;}
const HermiteGalsCertifiedInventoryReport& HermiteGalsCertifiedInventoryResult::report() const{return checked().report;}

HermiteGalsCertifiedInventoryResult assembleHermiteGalsCertifiedInventory(const QuadtreeMesh& mesh,
    const HermiteGalsPhaseFluxIntervalSource& interval,const std::vector<HermiteGalsBoundedFaceIntegral>& input,
    const HermiteGalsIntegratedInventoryConfig& config)
{
    interval.validateMesh(mesh);const auto& old=interval.previousAccepted();const auto& p=old.primary();
    require(interval.timeBegin()==old.timeLayer().time(),"certified inventory interval must begin at old accepted time");
    require(p.geometryCells.size()==static_cast<std::size_t>(mesh.cellIndexNum())
        && p.phaseVolumes.size()==p.geometryCells.size(),"certified inventory requires complete old support");
    FaceSegmentGraph graph;graph.build(mesh);const auto domain=NativeLeafDomain::build(mesh,graph,p.measure==HermiteGalsMeasure::AxisymmetricVolume);
    std::vector<HermiteGalsIntegratedFaceTransfer> transfers;transfers.reserve(input.size());
    std::vector<long double> exchange(p.phaseVolumes.size(),0);std::size_t index=0;
    const auto nativeFaces=domain.faces(mesh);
    for(const auto face:nativeFaces){
        require(index<input.size() && input[index].transfer().face().nativeRow==face.nativeRow,
            "missing/duplicate/unordered certified native face");
        const auto& certificate=input[index++];const auto& transfer=certificate.transfer();
        require(transfer.interval().sameInterval(interval),"certified face has a different interval source");
        const auto* inventory=certificate.inventorySource();
        require(inventory && inventory->source().identity()==old.identity()
            && &inventory->source().geometry().snapshot()==&old.geometry().snapshot(),"certified face lacks the same old inventory source");
        require(transfer.face().negativeCell==face.negativeCell && transfer.face().positiveCell==face.positiveCell,
            "certified native face owner mismatch");
        transfers.push_back(transfer);const auto volume=transfer.signedLiquidVolume();
        if(face.negativeCell>=0)exchange[face.negativeCell]-=volume;
        if(face.positiveCell>=0)exchange[face.positiveCell]+=volume;
    }
    require(index==input.size(),"foreign/inactive certified native face");
    HermiteGalsCertifiedInventoryReport report;bool needsSelection=false;
    for(std::size_t cell=0;cell<p.phaseVolumes.size();++cell){
        require(p.phaseVolumes[cell].cell==static_cast<int>(cell) && p.geometryCells[cell]==static_cast<int>(cell),"unordered certified old inventory");
        const long double capacity=domain.cellVolumes()[cell],value=p.phaseVolumes[cell].liquid+exchange[cell],
            roundoff=Detail::hermiteInventoryRoundoffBudget*capacity;
        require(std::isfinite(value),"nonfinite raw certified inventory");
        needsSelection=needsSelection || value < -roundoff || value > capacity*(1+Detail::hermiteInventoryRoundoffBudget);
    }
    report.selected=needsSelection;
    if(needsSelection){
        using Arc=LAMP::Discrete::BoundedCirculationArc;
        const auto reservoir=p.phaseVolumes.size();std::vector<Arc> arcs;arcs.reserve(2*input.size()+2*reservoir);
        const auto binary64=[](long double value){const double represented=static_cast<double>(value);
            require(std::isfinite(represented) && static_cast<long double>(represented)==value,"certified circulation requires binary64 inputs");return represented;};
        for(std::size_t k=0;k<transfers.size();++k){const auto& transfer=transfers[k];const auto& f=transfer.face();
            const auto& v=transfer.integratedVolumes();const auto& b=input[k].bounds().volumes;
            const auto negative=f.negativeCell<0?reservoir:static_cast<std::size_t>(f.negativeCell),
                positive=f.positiveCell<0?reservoir:static_cast<std::size_t>(f.positiveCell);
            const auto append=[&](std::size_t from,std::size_t to,long double total,long double liquid,int component){
                const double initial=binary64(liquid),maximum=binary64(total);
                const bool boundary=from==reservoir || to==reservoir;
                const double lo=boundary?initial:b[component].lower,hi=boundary?initial:std::min(maximum,b[component].upper);
                require(lo>=0 && lo<=initial && initial<=hi,"certified representative outside retained subset interval");
                arcs.push_back({from,to,lo,hi,initial});
            };
            append(negative,positive,v.integratedTotalPositiveVolume,v.integratedLiquidPositiveVolume,2);
            append(positive,negative,v.integratedTotalNegativeVolume,v.integratedLiquidNegativeVolume,3);
        }
        for(std::size_t cell=0;cell<reservoir;++cell){const double before=p.phaseVolumes[cell].liquid,capacity=domain.cellVolumes()[cell];
            require(before>=0 && before<=capacity,"certified source inventory outside capacity");
            arcs.push_back({reservoir,cell,before,before,before});arcs.push_back({cell,reservoir,0,capacity,before});
        }
        const auto selected=LAMP::Discrete::solveBoundedCirculation(reservoir+1,arcs);report.augmentations=selected.augmentations;
        for(std::size_t k=0;k<transfers.size();++k){const auto& transfer=transfers[k];const auto& face=transfer.face();
            auto v=transfer.integratedVolumes();const long double a=selected.flows[2*k],b=selected.flows[2*k+1];
            if(a==v.integratedLiquidPositiveVolume && b==v.integratedLiquidNegativeVolume)continue;
            require(face.negativeCell>=0 && face.positiveCell>=0,"certified selection changed fixed boundary");
            const auto change=std::max(std::abs(a-v.integratedLiquidPositiveVolume),std::abs(b-v.integratedLiquidNegativeVolume));
            report.maximumRelativeSelectionChange=std::max(report.maximumRelativeSelectionChange,
                change/std::min(domain.cellVolumes()[face.negativeCell],domain.cellVolumes()[face.positiveCell]));++report.changedFaces;
            v.integratedLiquidPositiveVolume=a;v.integratedLiquidNegativeVolume=b;
            transfers[k]=HermiteGalsIntegratedFaceTransfer(mesh,interval,face.nativeRow,transfer.orientation(),v);
        }
    }
    auto assembled=HermiteGalsIntegratedInventoryResult::assemble(mesh,interval,transfers,config,&input);
    for(std::size_t k=0;k<input.size();++k){const auto& before=input[k].transfer().integratedVolumes();const auto& after=assembled.faces()[k].integratedVolumes();
        const auto& bounds=input[k].bounds().volumes;const auto& face=assembled.faces()[k].face();
        require(after.integratedTotalPositiveVolume==before.integratedTotalPositiveVolume
            && after.integratedTotalNegativeVolume==before.integratedTotalNegativeVolume,"certified assembly changed total flux");
        const auto verifyComponent=[&](long double value,long double selected,int component){
            if(value>=bounds[component].lower && value<=bounds[component].upper)return;
            std::ostringstream why;why.precision(17);
            why<<"[Interface/HermiteGalsTemporalPhaseFlux geometry="<<old.identity().value()<<"] rounded certified assembly escaped original integral bounds"
                <<" face="<<face.nativeRow<<" component="<<component<<" lower="<<bounds[component].lower
                <<" selected="<<selected<<" final="<<value<<" upper="<<bounds[component].upper
                <<" repairs="<<assembled.report().repairedFaces;
            throw std::invalid_argument(why.str());
        };
        verifyComponent(after.integratedLiquidPositiveVolume,transfers[k].integratedVolumes().integratedLiquidPositiveVolume,2);
        verifyComponent(after.integratedLiquidNegativeVolume,transfers[k].integratedVolumes().integratedLiquidNegativeVolume,3);
        require((face.negativeCell>=0 && face.positiveCell>=0) || (after.integratedLiquidPositiveVolume==before.integratedLiquidPositiveVolume
            && after.integratedLiquidNegativeVolume==before.integratedLiquidNegativeVolume),"certified assembly changed boundary history flux");
    }
    return HermiteGalsCertifiedInventoryResult(std::make_shared<const HermiteGalsCertifiedInventoryResult::Data>(input,std::move(assembled),report));
}
}
