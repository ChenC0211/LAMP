#include "Modules/Interface/HermiteGalsFirstJetStep.h"

#include "Modules/Interface/detail/HermiteGalsFirstJet.h"
#include "Core/Geometry/ActiveQuarticHermitePatch.h"
#include "Core/Geometry/ActiveQuarticHermiteSkeleton.h"
#include "Core/Geometry/CompactQuarticDerivative.h"
#include "Core/Geometry/Detail/ExactExpansion.h"
#include "Core/Mesh/UniformLeafGrid.h"
#include <algorithm>
#include <cstring>
#include <cstdint>
#include <cmath>
#include <chrono>
#include <cstddef>
#include <exception>
#include <iomanip>
#include <limits>
#include <map>
#include <mutex>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <type_traits>
#include <string>
#include <utility>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace LAMP::Interface
{
using Detail::requireFirstJet;
HermiteGalsPublication initializeHermiteGalsFirstJet(const QuadtreeMesh& mesh,bool radial,
    const HermiteGalsInitialFirstJetQuery& query,const HermiteGalsPublicationConfig& config)
{
    if(!query)throw std::invalid_argument("FirstJet initial query is empty");
    const auto initial=[&](Vector2D p){
        auto j=query(p);
        for(const auto& plane:config.volume.symmetryPlanes){
            const double coordinate=plane.normalAxis==0?p.x:p.y;
            if(coordinate!=plane.coordinate)continue;
            double& normal=plane.normalAxis==0?j.gradient.x:j.gradient.y;
            const double roundoff=64*std::numeric_limits<double>::epsilon()*
                std::max({1.,std::abs(j.gradient.x),std::abs(j.gradient.y)});
            if(std::abs(normal)>roundoff)throw std::invalid_argument("initial FirstJet violates declared homogeneous Neumann boundary");
            // Canonical exact boundary data; preserve value and tangent.
            normal=0;
        }
        return j;
    };
    const auto point=[&](int v){return mesh.verts(mesh.vertIndexToLabel[v]).coord();};
    LAMP::Geometry::ActiveQuarticHermiteSkeletonInput input;
    input.activeCells.assign(mesh.cellIndexNum(),1);input.nodes.resize(mesh.vertIndexNum());
    std::set<std::array<int,2>> sides;
    for(int c=0;c<mesh.cellIndexNum();++c){const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);
        for(int k=0;k<4;++k){std::array<int,2> key{mesh.verts(cell.cellToVert[k]).index(),mesh.verts(cell.cellToVert[(k+1)%4]).index()};
            const auto a=point(key[0]),b=point(key[1]);if(a.y==b.y?a.x>b.x:a.y>b.y)std::swap(key[0],key[1]);sides.insert(key);}}
    for(const auto& key:sides)input.edgeMidpoints.push_back({key,0,0});
    // Use the existing master/hanging owner engine; the zero data are not a
    // geometry or an initial scalar and never leave this topology query.
    const LAMP::Geometry::ActiveQuarticHermiteSkeleton topology(mesh,input,radial);
    HermiteGalsPrimaryState primary;primary.transport=HermiteGalsTransport::FirstJet;
    primary.measure=radial?HermiteGalsMeasure::AxisymmetricVolume:HermiteGalsMeasure::CartesianArea;
    const std::array<int,2> none{-1,-1};
    for(const auto& v:topology.vertices())if(v.masterSide==none){const auto j=initial(point(v.index));primary.vertices.push_back({v.index,j.phi,j.gradient.x,j.gradient.y});}
    for(const auto& e:topology.edges())if(e.parentSide==none){const auto a=point(e.vertices[0]),b=point(e.vertices[1]);const auto j=initial(a+(b-a)*.5);
        primary.edges.push_back({e.vertices,j.phi,a.y==b.y?j.gradient.y:j.gradient.x});}
    std::sort(primary.vertices.begin(),primary.vertices.end(),[](const auto& a,const auto& b){return a.vertex<b.vertex;});
    std::sort(primary.edges.begin(),primary.edges.end(),[](const auto& a,const auto& b){return a.vertices<b.vertices;});
    std::vector<double> centers;
    for(int c=0;c<mesh.cellIndexNum();++c){primary.geometryCells.push_back(c);centers.push_back(initial(mesh.cells(mesh.cellIndexToLabel[c]).centerCoord()).phi);}
    HermiteGalsCandidate candidate(mesh,{},primary);const auto mixed=closeHermiteGalsMixedAmr(mesh,candidate,config.volume.symmetryPlanes);
    const HermiteGalsCenterPrediction center{candidate.identity(),centers};
    return publishHermiteGalsFirstJetGeometry(mesh,std::move(candidate),mixed,center,config);
}

HermiteGalsResolvedMotion::HermiteGalsResolvedMotion(HermiteGalsIdentity source,HermiteGalsTimeLayer sourceTime,
    HermiteGalsIdentity velocity,HermiteGalsTimeLayer targetTime,Query query)
    :source_(source),velocity_(velocity),sourceTime_(sourceTime),targetTime_(targetTime),query_(std::move(query))
{
    requireFirstJet(source.value()!=0&&velocity.value()!=0&&bool(query_),"unbound resolved motion");
    requireFirstJet(targetTime.time()>sourceTime.time()&&targetTime.begin==sourceTime.time(),"resolved motion interval differs from source");
}
LevelSetBoundaryTrace::JetResult HermiteGalsResolvedMotion::trace(Vector2D point) const
{
    requireFirstJet(LevelSetBoundaryTrace::Detail::finite(point),"nonfinite resolved arrival point");
    auto result=query_(point);const auto& t=result.trajectory;
    const double dt=targetTime_.time()-sourceTime_.time();
    requireFirstJet(LevelSetBoundaryTrace::Detail::finite(t.point)&&std::isfinite(t.elapsed)&&t.elapsed>=0&&t.elapsed<=dt,
        "invalid resolved trajectory point/time");
    requireFirstJet(result.acceptedSubsteps>=0,"invalid resolved substep count");
    if(t.boundaryHit)requireFirstJet(t.side.has_value()&&!result.derivative,"boundary motion lacks side or fabricates interior J");
    else requireFirstJet(t.elapsed==dt&&!t.side&&result.derivative&&LevelSetBoundaryTrace::Detail::finite(*result.derivative),
        "incomplete/nonfinite resolved interior map");
    return result;
}
struct HermiteGalsFirstJetStep::Data
{
    Data(HermiteGalsIdentity old, HermiteGalsIdentity velocity, HermiteGalsPublication next,
        HermiteGalsFirstJetStepReport diagnostic)
        : source(old),arrival(velocity),publication(std::move(next)),report(diagnostic) {}
    HermiteGalsIdentity source,arrival;
    HermiteGalsPublication publication;
    HermiteGalsFirstJetStepReport report;
};
HermiteGalsFirstJetStep::HermiteGalsFirstJetStep(std::shared_ptr<const Data> data):data_(std::move(data)) {}
const HermiteGalsFirstJetStep::Data& HermiteGalsFirstJetStep::checked() const
{requireFirstJet(bool(data_),"[Interface/HermiteGalsFirstJetStep] moved-from result");return *data_;}
HermiteGalsIdentity HermiteGalsFirstJetStep::sourceIdentity() const {return checked().source;}
HermiteGalsIdentity HermiteGalsFirstJetStep::arrivalVelocityIdentity() const {return checked().arrival;}
const HermiteGalsPublication& HermiteGalsFirstJetStep::publication() const {return checked().publication;}
const HermiteGalsFirstJetStepReport& HermiteGalsFirstJetStep::report() const {return checked().report;}

struct HermiteGalsFirstJetSample::Data
{
    Data(HermiteGalsIdentity old,HermiteGalsIdentity velocity,HermiteGalsState scalar,HermiteGalsFirstJetStepReport diagnostic)
        :source(old),arrival(velocity),state(std::move(scalar)),report(diagnostic) {}
    HermiteGalsIdentity source,arrival;
    HermiteGalsState state; // Private owning revision; never exposed for commit.
    HermiteGalsFirstJetStepReport report;
};
HermiteGalsFirstJetSample::HermiteGalsFirstJetSample(std::shared_ptr<const Data> data):data_(std::move(data)) {}
const HermiteGalsFirstJetSample::Data& HermiteGalsFirstJetSample::checked() const
{requireFirstJet(bool(data_),"[Interface/HermiteGalsFirstJetSample] moved-from result");return *data_;}
HermiteGalsIdentity HermiteGalsFirstJetSample::identity() const{return checked().state.identity();}
HermiteGalsTimeLayer HermiteGalsFirstJetSample::timeLayer() const{return checked().state.timeLayer();}
HermiteGalsMeasure HermiteGalsFirstJetSample::measure() const{return checked().state.primary().measure;}
HermiteGalsGeometryView HermiteGalsFirstJetSample::geometry() const{return checked().state.geometry();}
HermiteGalsIdentity HermiteGalsFirstJetSample::sourceIdentity() const{return checked().source;}
HermiteGalsIdentity HermiteGalsFirstJetSample::arrivalVelocityIdentity() const{return checked().arrival;}
const HermiteGalsFirstJetStepReport& HermiteGalsFirstJetSample::report() const{return checked().report;}
void HermiteGalsFirstJetSample::validateMesh(const QuadtreeMesh& mesh) const{checked().state.validateMesh(mesh);}

struct HermiteGalsFirstJetLocalSample::Data
{
    Data(HermiteGalsIdentity old,HermiteGalsIdentity velocity,HermiteGalsState scalar,HermiteGalsFirstJetStepReport diagnostic)
        :source(old),arrival(velocity),state(std::move(scalar)),report(diagnostic) {}
    HermiteGalsIdentity source,arrival;
    HermiteGalsState state; // Only the genuine restricted scalar, never exposed as state.
    HermiteGalsFirstJetStepReport report;
};
HermiteGalsFirstJetLocalSample::HermiteGalsFirstJetLocalSample(std::shared_ptr<const Data> data):data_(std::move(data)) {}
const HermiteGalsFirstJetLocalSample::Data& HermiteGalsFirstJetLocalSample::checked() const
{requireFirstJet(bool(data_),"[Interface/HermiteGalsFirstJetLocalSample] moved-from result");return *data_;}
HermiteGalsIdentity HermiteGalsFirstJetLocalSample::identity() const{return checked().state.identity();}
HermiteGalsTimeLayer HermiteGalsFirstJetLocalSample::timeLayer() const{return checked().state.timeLayer();}
HermiteGalsMeasure HermiteGalsFirstJetLocalSample::measure() const{return checked().state.primary().measure;}
HermiteGalsGeometryView HermiteGalsFirstJetLocalSample::geometry() const{return checked().state.geometry();}
const std::vector<int>& HermiteGalsFirstJetLocalSample::nativeCells() const{return checked().state.primary().geometryCells;}
HermiteGalsIdentity HermiteGalsFirstJetLocalSample::sourceIdentity() const{return checked().source;}
HermiteGalsIdentity HermiteGalsFirstJetLocalSample::arrivalVelocityIdentity() const{return checked().arrival;}
const HermiteGalsFirstJetStepReport& HermiteGalsFirstJetLocalSample::report() const{return checked().report;}
void HermiteGalsFirstJetLocalSample::validateMesh(const QuadtreeMesh& mesh) const{checked().state.validateMesh(mesh);}

bool supportsHermiteGalsFirstJetLocalSample(const QuadtreeMesh& mesh,const HermiteGalsPublication& source)
{
    source.validateMesh(mesh);
    const auto& state=source.state();const auto& scalar=state.geometry().snapshot();
    if(state.primary().measure!=HermiteGalsMeasure::CartesianArea || scalar.radialAxisSymmetry()
        || state.primary().transport!=HermiteGalsTransport::FirstJet || !state.primary().phaseVolumes.empty()
        || state.primary().geometryCells.size()!=static_cast<std::size_t>(mesh.cellIndexNum())
        || !scalar.usesNativeLeaves() || scalar.hasCornerSectors()
        || scalar.reconstruction()!=LevelSetHermiteGeometry::Reconstruction::EnrichedQuartic)return false;
    try {const UniformLeafGrid grid(mesh,CartesianGridScope::UniformLeaves);return grid.nx()>=2 && grid.ny()>=2;}
    catch(const std::invalid_argument&){return false;}
}

HermiteGalsFirstJetStep prepareHermiteGalsFirstJetStep(const QuadtreeMesh& mesh,
    const HermiteGalsPublication& source,const LevelSetBoundaryTrace::Domain& domain,
    const HermiteGalsVelocityView& previous,const HermiteGalsVelocityView& arrival,
    const HermiteGalsBoundaryFirstJetQuery& boundary,const HermiteGalsPublicationConfig& config,
    HermiteGalsPointExecution execution,const HermiteGalsPointSamplingConfig& sampling)
{
    source.validateMesh(mesh);
    using Position=HermiteGalsTimeLayer::Position;
    const auto a=source.state().timeLayer(),b=arrival.timeLayer();
    requireFirstJet(a.position==Position::Initial||a.position==Position::End,"source must be accepted Initial/End");
    requireFirstJet(b.position==Position::Midpoint||b.position==Position::End,"target must be Midpoint/End");
    requireFirstJet(a.step<std::numeric_limits<std::int64_t>::max()&&b.step==a.step+1,"target must be the next step");
    requireFirstJet(source.state().primary().transport==HermiteGalsTransport::FirstJet&&source.state().primary().phaseVolumes.empty(),
        "requires FirstJet source without inventory");
    Detail::validateFirstJetSource(mesh,source.state(),source.state().geometry().snapshot(),domain,previous,arrival);
    const HermiteGalsResolvedMotion motion(source.state().identity(),source.state().timeLayer(),arrival.geometryIdentity(),arrival.timeLayer(),
        [domain,previous,arrival](Vector2D point){
            const auto t=traceHermiteGalsWithJacobian(domain,point,previous,arrival);
            return LevelSetBoundaryTrace::JetResult{{t.departurePoint(),t.elapsed(),t.boundarySide(),t.boundaryHit(),false},t.jacobian(),t.acceptedSubsteps()};
        });
    return prepareHermiteGalsFirstJetStep(mesh,source,motion,boundary,config,HermiteGalsBoundaryTreatment::InflowHistory,
        HermiteGalsTimeProjection::CharacteristicPullback,execution,sampling);
}
namespace
{
// Extra context is constructed only on the existing boundary rejection path.
void requireInteriorFirstJetTrace(const LevelSetBoundaryTrace::JetResult& trace,Vector2D arrival)
{
    if(!trace.trajectory.boundaryHit)return;
    const auto& trajectory=trace.trajectory;
    const char* side="none";
    if(trajectory.side)switch(*trajectory.side) {
        case Boundary::non:side="non";break;
        case Boundary::west:side="west";break;
        case Boundary::south:side="south";break;
        case Boundary::east:side="east";break;
        case Boundary::north:side="north";break;
        case Boundary::southWest:side="southWest";break;
        case Boundary::southEast:side="southEast";break;
        case Boundary::northEast:side="northEast";break;
        case Boundary::northWest:side="northWest";break;
    }
    std::ostringstream message;
    message<<std::setprecision(17)<<"interior Jacobian unavailable: boundary event"
        <<" arrival=("<<arrival.x<<","<<arrival.y<<")"
        <<" hit=("<<trajectory.point.x<<","<<trajectory.point.y<<")"
        <<" side="<<side<<" elapsed="<<trajectory.elapsed
        <<" substeps="<<trace.acceptedSubsteps;
    throw std::invalid_argument(message.str());
}
bool noQueryBoundaryTreatment(HermiteGalsBoundaryTreatment treatment)
{
    return treatment==HermiteGalsBoundaryTreatment::InflowHistory
        ||treatment==HermiteGalsBoundaryTreatment::SymmetryFirstJet;
}
HermiteGalsBoundaryFirstJet applySymmetryFirstJet(Vector2D point,HermiteGalsBoundaryFirstJet jet,
    const HermiteGalsPublicationConfig& config,HermiteGalsBoundaryTreatment treatment)
{
    if(treatment!=HermiteGalsBoundaryTreatment::SymmetryFirstJet)return jet;
    for(const auto& plane:config.volume.symmetryPlanes)
        if((plane.normalAxis==0?point.x:point.y)==plane.coordinate)
            (plane.normalAxis==0?jet.gradient.x:jet.gradient.y)=0;
    return jet;
}
template<class Point>
void validateSymmetryFirstJet(Vector2D lower,Vector2D upper,const HermiteGalsPrimaryState& primary,
    const Point& point,const HermiteGalsBoundaryFirstJetQuery& boundary,
    const HermiteGalsPublicationConfig& config,HermiteGalsBoundaryTreatment treatment)
{
    if(treatment!=HermiteGalsBoundaryTreatment::SymmetryFirstJet)return;
    requireFirstJet(!boundary,"symmetry first-jet does not accept prescribed boundary data");
    requireFirstJet(!config.volume.symmetryPlanes.empty(),"symmetry first-jet requires declared outer planes");
    const auto coordinate=[](Vector2D q,int axis){return axis==0?q.x:q.y;};
    for(std::size_t k=0;k<config.volume.symmetryPlanes.size();++k){const auto& plane=config.volume.symmetryPlanes[k];
        requireFirstJet((plane.normalAxis==0||plane.normalAxis==1)&&std::isfinite(plane.coordinate),"invalid symmetry plane");
        requireFirstJet(plane.coordinate==coordinate(lower,plane.normalAxis)
            ||plane.coordinate==coordinate(upper,plane.normalAxis),"symmetry plane must be an actual outer domain section");
        requireFirstJet(!(primary.measure==HermiteGalsMeasure::AxisymmetricVolume&&plane.normalAxis==0&&plane.coordinate==0),
            "radial axis already has its own symmetry identity");
        for(std::size_t j=0;j<k;++j)requireFirstJet(!(plane==config.volume.symmetryPlanes[j]),"duplicate symmetry plane");
        bool found=false;
        for(const auto& edge:primary.edges){const auto a=point(edge.vertices[0]),b=point(edge.vertices[1]);
            if(coordinate(a,plane.normalAxis)==plane.coordinate&&coordinate(b,plane.normalAxis)==plane.coordinate){
                found=true;requireFirstJet(edge.dn==0,"symmetry source requires exactly zero primary edge normal derivative");}}
        requireFirstJet(found,"symmetry plane has no active boundary trace");
        for(const auto& vertex:primary.vertices)if(coordinate(point(vertex.vertex),plane.normalAxis)==plane.coordinate)
            requireFirstJet((plane.normalAxis==0?vertex.dx:vertex.dy)==0,
                "symmetry source requires exactly zero primary vertex normal derivative");
    }
}
using LocalPointKey=std::pair<double,double>;
struct alignas(64) LocalPointStats
{
    std::size_t pointQueries=0,traceQueries=0,numericalQueries=0,valueFallbacks=0,gradientFallbacks=0;
    double seconds=0,valueError=0,gradientError=0;
    int threads=1;
};
void mergeLocalPointStats(HermiteGalsFirstJetStepReport& report,const LocalPointStats& stats)
{
    report.evaluatedPointQueries+=stats.pointQueries;report.pointTraceQueries+=stats.traceQueries;
    report.numericalPointQueries+=stats.numericalQueries;report.numericalValueFallbacks+=stats.valueFallbacks;
    report.numericalGradientFallbacks+=stats.gradientFallbacks;
    report.pointWorkerSeconds+=stats.seconds;report.pointThreads=std::max(report.pointThreads,stats.threads);
    report.maxAcceptedValueErrorBound=std::max(report.maxAcceptedValueErrorBound,stats.valueError);
    report.maxAcceptedGradientErrorBound=std::max(report.maxAcceptedGradientErrorBound,stats.gradientError);
}
auto localTransportJet(const LevelSetHermiteGeometry::Snapshot& scalar,Vector2D point,
    const HermiteGalsPointSamplingConfig& sampling,LocalPointStats& stats)
{
    ++stats.pointQueries;
    if(sampling.mode==HermiteGalsPointSampling::Exact)return scalar.scalarJetNative(point);
    ++stats.numericalQueries;
    const auto sampled=scalar.numericalScalarJetNative(point,sampling.valueAbsoluteTolerance,sampling.gradientAbsoluteTolerance);
    stats.valueFallbacks+=sampled.exactValueFallback;stats.gradientFallbacks+=sampled.exactGradientFallback;
    requireFirstJet(sampled.boundsResolved,"numerical transport point bounds are unresolved");
    for(int k=0;k<3;++k){const double bound=sampled.absoluteErrorBounds[k];
        requireFirstJet(std::isfinite(bound)&&bound>=0&&bound<=(k==0?sampling.valueAbsoluteTolerance:sampling.gradientAbsoluteTolerance),
            "numerical transport point error exceeds its physical absolute budget");}
    requireFirstJet(std::isfinite(sampled.jet.value)&&std::isfinite(sampled.jet.dx)&&std::isfinite(sampled.jet.dy),
        "nonfinite numerical transport point jet");
    stats.valueError=std::max(stats.valueError,sampled.absoluteErrorBounds[0]);
    stats.gradientError=std::max({stats.gradientError,sampled.absoluteErrorBounds[1],sampled.absoluteErrorBounds[2]});
    return sampled.jet;
}
double localTransportValue(const LevelSetHermiteGeometry::Snapshot& scalar,Vector2D point,
    const HermiteGalsPointSamplingConfig& sampling,LocalPointStats& stats)
{
    ++stats.pointQueries;
    if(sampling.mode==HermiteGalsPointSampling::Exact)return scalar.scalarJetNative(point).value;
    ++stats.numericalQueries;
    const auto sampled=scalar.numericalScalarValueNative(point,sampling.valueAbsoluteTolerance);
    stats.valueFallbacks+=sampled.exactValueFallback;
    requireFirstJet(sampled.boundsResolved,"numerical transport point bounds are unresolved");
    requireFirstJet(std::isfinite(sampled.absoluteErrorBound)&&sampled.absoluteErrorBound>=0
        &&sampled.absoluteErrorBound<=sampling.valueAbsoluteTolerance,
        "numerical transport point error exceeds its physical absolute budget");
    requireFirstJet(std::isfinite(sampled.value),"nonfinite numerical transport point value");
    stats.valueError=std::max(stats.valueError,sampled.absoluteErrorBound);return sampled.value;
}
// One invocation, one authenticated immutable source and target clock.
// Slots belong to existing first-map point indices, never a shared worker map.
struct SourceJetSlot
{
    std::array<std::uint64_t,4> key{}; // physical x/y and the original sampling budgets
    HermiteGalsPointSampling mode=HermiteGalsPointSampling::Exact;
    double value=0;
    std::optional<LevelSetHermiteGeometry::Snapshot::ScalarJet2> jet; // absent for a value-only result
};
std::uint64_t sourceJetDoubleBits(double value) noexcept
{
    static_assert(sizeof(double)==sizeof(std::uint64_t),"binary64 source point key required");
    std::uint64_t bits;std::memcpy(&bits,&value,sizeof(bits));return bits;
}
auto sourceTransportJet(const LevelSetHermiteGeometry::Snapshot& scalar,Vector2D point,
    const HermiteGalsPointSamplingConfig& sampling,LocalPointStats& stats,
    std::optional<SourceJetSlot>& slot,bool enabled,bool retain=true)
{
    if(!enabled)return localTransportJet(scalar,point,sampling,stats);
    LAMP::Geometry::Detail::ExactExpansion::validateActualEnvironment();
    const std::array<std::uint64_t,4> key{sourceJetDoubleBits(point.x),sourceJetDoubleBits(point.y),
        sourceJetDoubleBits(sampling.valueAbsoluteTolerance),sourceJetDoubleBits(sampling.gradientAbsoluteTolerance)};
    if(slot&&slot->key==key&&slot->mode==sampling.mode&&slot->jet){
        LAMP::Geometry::Detail::ExactExpansion::validateActualEnvironment();return *slot->jet;
    }
    const auto jet=localTransportJet(scalar,point,sampling,stats);
    if(retain)slot=SourceJetSlot{key,sampling.mode,jet.value,jet};
    LAMP::Geometry::Detail::ExactExpansion::validateActualEnvironment();return jet;
}
double sourceTransportValue(const LevelSetHermiteGeometry::Snapshot& scalar,Vector2D point,
    const HermiteGalsPointSamplingConfig& sampling,LocalPointStats& stats,
    std::optional<SourceJetSlot>& slot,bool enabled,bool retain=true)
{
    if(!enabled)return localTransportValue(scalar,point,sampling,stats);
    LAMP::Geometry::Detail::ExactExpansion::validateActualEnvironment();
    const std::array<std::uint64_t,4> key{sourceJetDoubleBits(point.x),sourceJetDoubleBits(point.y),
        sourceJetDoubleBits(sampling.valueAbsoluteTolerance),sourceJetDoubleBits(sampling.gradientAbsoluteTolerance)};
    if(slot&&slot->key==key&&slot->mode==sampling.mode){
        LAMP::Geometry::Detail::ExactExpansion::validateActualEnvironment();return slot->value;
    }
    const double value=localTransportValue(scalar,point,sampling,stats);
    if(retain)slot=SourceJetSlot{key,sampling.mode,value,{}};
    LAMP::Geometry::Detail::ExactExpansion::validateActualEnvironment();return value;
}
struct FirstPointRecord
{
    LevelSetBoundaryTrace::JetResult trace;
    std::optional<SourceJetSlot> sourceJet;
};
struct LocalPointTraces
{
    std::vector<LocalPointKey> keys;
    std::vector<std::optional<LevelSetBoundaryTrace::JetResult>> values;
    std::vector<std::optional<SourceJetSlot>> sourceJets;
    std::vector<std::exception_ptr> failures;
    explicit LocalPointTraces(const std::set<LocalPointKey>& points):keys(points.begin(),points.end()),values(keys.size()),sourceJets(keys.size()),failures(keys.size()){}
    std::size_t find(LocalPointKey key) const
    {
        const auto found=std::lower_bound(keys.begin(),keys.end(),key);
        return found!=keys.end()&&*found==key?static_cast<std::size_t>(found-keys.begin()):keys.size();
    }
    const LevelSetBoundaryTrace::JetResult& at(std::size_t index) const
    {
        requireFirstJet(index<keys.size(),"missing prepared local first-map point");
        if(failures[index])std::rethrow_exception(failures[index]);
        requireFirstJet(values[index].has_value(),"missing prepared local first-map result");return *values[index];
    }
};
struct LocalPointJets
{
    std::vector<LocalPointKey> keys;
    std::vector<std::optional<HermiteGalsBoundaryFirstJet>> values;
    std::vector<std::optional<double>> centerValues;
    std::vector<std::exception_ptr> failures;
    explicit LocalPointJets(const std::set<LocalPointKey>& points):keys(points.begin(),points.end()),values(keys.size()),centerValues(keys.size()),failures(keys.size()){}
    std::size_t find(Vector2D point) const
    {
        const auto key=std::make_pair(point.x,point.y);const auto found=std::lower_bound(keys.begin(),keys.end(),key);
        requireFirstJet(found!=keys.end()&&*found==key,"missing prepared local jet point");
        return static_cast<std::size_t>(found-keys.begin());
    }
    HermiteGalsBoundaryFirstJet at(Vector2D point) const
    {
        const auto index=find(point);if(failures[index])std::rethrow_exception(failures[index]);
        requireFirstJet(values[index].has_value(),"missing prepared local point jet");return *values[index];
    }
    double valueAt(Vector2D point) const
    {
        const auto index=find(point);if(failures[index])std::rethrow_exception(failures[index]);
        if(centerValues[index])return *centerValues[index];
        requireFirstJet(values[index].has_value(),"missing prepared local point jet");return values[index]->phi;
    }
};
template<class Job>
void localPointBatch(std::size_t count,std::vector<std::exception_ptr>& failures,
    HermiteGalsFirstJetStepReport& report,const Job& job)
{
    requireFirstJet(count<=static_cast<std::size_t>(std::numeric_limits<std::ptrdiff_t>::max()),"too many local point jobs");
    LAMP::Geometry::Detail::ExactExpansion::validateEnvironment();
    int slots=1;
#ifdef _OPENMP
    slots=std::max(1,omp_get_max_threads());
#endif
    std::vector<LocalPointStats> stats(static_cast<std::size_t>(slots));
#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        int worker=0,threads=1;
#ifdef _OPENMP
        worker=omp_get_thread_num();threads=omp_get_num_threads();
#endif
        auto& own=stats[static_cast<std::size_t>(worker)];
        const auto start=std::chrono::steady_clock::now();bool worked=false;
#ifdef _OPENMP
#pragma omp for schedule(static) nowait
#endif
        for(std::ptrdiff_t q=0;q<static_cast<std::ptrdiff_t>(count);++q){const auto index=static_cast<std::size_t>(q);
            worked=true;
            try{LAMP::Geometry::Detail::ExactExpansion::validateEnvironment();job(index,own);
                LAMP::Geometry::Detail::ExactExpansion::validateEnvironment();}
            catch(...){failures[index]=std::current_exception();}
        }
        // Passive timing is one elapsed interval per active worker, excluding
        // the final team barrier; point failures retain their original slots.
        if(worked){own.threads=threads;
            own.seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();}
    }
    LAMP::Geometry::Detail::ExactExpansion::validateEnvironment();
    report.pointJobs+=count;for(const auto& own:stats)mergeLocalPointStats(report,own);
    // Exceptions are consumed by the ORIGINAL serial payload/closure order,
    // rather than the completion order of the workers.
}
template<class Value> struct FullPointResults
{
    std::vector<std::optional<LevelSetBoundaryTrace::JetResult>> traces;
    std::vector<std::optional<Value>> values;
    std::vector<std::optional<SourceJetSlot>> sourceJets;
    std::vector<std::exception_ptr> failures;
    explicit FullPointResults(std::size_t count):traces(count),values(count),sourceJets(count),failures(count){}
    Value at(std::size_t index) const
    {
        if(failures.at(index))std::rethrow_exception(failures[index]);
        requireFirstJet(traces.at(index).has_value()&&values.at(index).has_value(),"missing prepared full first jet");
        return *values[index];
    }
};
template<class GeometryResult> struct FirstJetPrepared
{
    HermiteGalsIdentity source,arrival;
    GeometryResult result;
    HermiteGalsFirstJetStepReport report;
};

// Immutable physical half-grid support shared by full, local and cardinal
// updates. Only topology is retained between samples, never scalar/time data.
struct IncrementGradientGrid
{
    int nx,ny;
    std::vector<Vector2D> points;
    std::map<LocalPointKey,std::array<int,2>> positions;
    template<class Vertex,class Center>
    IncrementGradientGrid(int x,int y,const Vertex& vertex,const Center& center):nx(x),ny(y)
    {
        requireFirstJet(nx>=2&&ny>=2,"compatible increment requires five half-grid nodes per axis");
        points.resize((2*nx+1)*(2*ny+1));
        for(int j=0;j<=2*ny;++j)for(int i=0;i<=2*nx;++i){Vector2D q;
            if(i%2&&j%2)q=center(i/2,j/2);
            else if(i%2){const auto a=vertex(i/2,j/2),b=vertex(i/2+1,j/2);q=a+(b-a)*.5;}
            else if(j%2){const auto a=vertex(i/2,j/2),b=vertex(i/2,j/2+1);q=a+(b-a)*.5;}
            else q=vertex(i/2,j/2);
            points[j*(2*nx+1)+i]=q;
            requireFirstJet(positions.emplace(LocalPointKey{q.x,q.y},std::array<int,2>{i,j}).second,
                "duplicate physical half-grid node");
        }
    }
    Vector2D at(int i,int j) const{return points[j*(2*nx+1)+i];}
    bool symmetric(Vector2D q,int axis,const HermiteGalsPublicationConfig& config) const
    {for(const auto& p:config.volume.symmetryPlanes)if(p.normalAxis==axis&&(axis==0?q.x:q.y)==p.coordinate)return true;return false;}
    void expand(const std::set<LocalPointKey>& jets,std::set<LocalPointKey>& values,const HermiteGalsPublicationConfig& config) const
    {
        for(const auto& key:jets){const Vector2D q{key.first,key.second};const auto ij=positions.at(key);
            for(int axis=0;axis<2;++axis)if(!symmetric(q,axis,config)){
                const int target=ij[axis],count=2*(axis==0?nx:ny)+1,begin=std::clamp(target-2,0,count-5);
                for(int k=begin;k<begin+5;++k){const auto p=axis==0?at(k,ij[1]):at(ij[0],k);values.emplace(p.x,p.y);}
            }
        }
    }
};
class IncrementGradientProjection
{
    const IncrementGradientGrid& grid;
    const LevelSetHermiteGeometry::Snapshot& old;
    const HermiteGalsPublicationConfig& config;
    std::function<double(Vector2D)> newValue;
    std::map<LocalPointKey,LevelSetHermiteGeometry::Snapshot::ScalarJet2> oldJets;
    std::map<LocalPointKey,double> increments;
    std::map<LocalPointKey,HermiteGalsBoundaryFirstJet> jets;
    const LevelSetHermiteGeometry::Snapshot::ScalarJet2& source(Vector2D q)
    {
        const LocalPointKey key{q.x,q.y};auto found=oldJets.find(key);
        if(found==oldJets.end())found=oldJets.emplace(key,old.scalarJetNative(q)).first;
        return found->second;
    }
    double delta(Vector2D q)
    {
        const LocalPointKey key{q.x,q.y};auto found=increments.find(key);
        if(found==increments.end())found=increments.emplace(key,newValue(q)-source(q).value).first;
        return found->second;
    }
public:
    IncrementGradientProjection(const IncrementGradientGrid& g,const LevelSetHermiteGeometry::Snapshot& s,
        const HermiteGalsPublicationConfig& c,std::function<double(Vector2D)> value):grid(g),old(s),config(c),newValue(std::move(value)){}
    HermiteGalsBoundaryFirstJet at(Vector2D q)
    {
        const LocalPointKey key{q.x,q.y};const auto found=jets.find(key);if(found!=jets.end())return found->second;
        const auto ij=grid.positions.at(key);const auto base=source(q);Vector2D gradient{base.dx,base.dy};
        for(int axis=0;axis<2;++axis){double change=0;
            if(grid.symmetric(q,axis,config)){
                requireFirstJet((axis==0?base.dx:base.dy)==0,"compatible increment requires exact source symmetry");
            }else{
                const auto pointAt=[&](int k){return axis==0?grid.at(k,ij[1]):grid.at(ij[0],k);};
                change=LAMP::Geometry::compactQuarticDerivativeAtNode(ij[axis],2*(axis==0?grid.nx:grid.ny)+1,
                    [&](int k){const auto p=pointAt(k);return axis==0?p.x:p.y;},
                    [&](int k){return delta(pointAt(k));},requireFirstJet);
            }
            if(axis==0)gradient.x+=change;else gradient.y+=change;
        }
        HermiteGalsBoundaryFirstJet result{newValue(q),gradient};
        requireFirstJet(std::isfinite(result.phi)&&LevelSetBoundaryTrace::Detail::finite(gradient),"nonfinite compatible increment");
        jets.emplace(key,result);return result;
    }
    std::size_t pointQueries() const{return oldJets.size();}
};
void validateGradientUpdate(const HermiteGalsPublicationConfig& config,const LevelSetHermiteGeometry::Snapshot& scalar,
    bool compensated,bool boundary,HermiteGalsBoundaryTreatment treatment)
{
    using Update=HermiteGalsGradientUpdate;
    requireFirstJet(config.gradientUpdate==Update::CharacteristicPullback||config.gradientUpdate==Update::CompatibleIncrement,"unknown gradient update");
    if(config.gradientUpdate==Update::CompatibleIncrement){
        requireFirstJet(!compensated,"compatible increment requires characteristic temporal projection");
        requireFirstJet(!boundary&&noQueryBoundaryTreatment(treatment),"compatible increment excludes external boundary callbacks");
        requireFirstJet(config.scalarRepresentation==HermiteGalsScalarRepresentation::SingleQuartic
            &&scalar.usesNativeLeaves()&&!scalar.usesC2Macro()&&!scalar.hasCornerSectors(),"compatible increment requires SingleQuartic source");
    }
}
void projectFullGradient(const QuadtreeMesh& mesh,HermiteGalsPrimaryState& primary,const std::vector<double>& centers,
    const LevelSetHermiteGeometry::Snapshot& old,const HermiteGalsPublicationConfig& config,HermiteGalsFirstJetStepReport& report)
{
    const UniformLeafGrid grid(mesh,CartesianGridScope::UniformLeaves);
    const auto vertex=[&](int i,int j){return mesh.verts(grid.vertices()[j*(grid.nx()+1)+i]).coord();};
    const auto center=[&](int i,int j){return mesh.cells(grid.cells()[j*grid.nx()+i]).centerCoord();};
    const IncrementGradientGrid support(grid.nx(),grid.ny(),vertex,center);
    std::map<LocalPointKey,double> values;
    const auto point=[&](int v){return mesh.verts(mesh.vertIndexToLabel[v]).coord();};
    for(const auto& v:primary.vertices){const auto q=point(v.vertex);values.emplace(LocalPointKey{q.x,q.y},v.phi);}
    for(const auto& e:primary.edges){const auto a=point(e.vertices[0]),b=point(e.vertices[1]),q=a+(b-a)*.5;values.emplace(LocalPointKey{q.x,q.y},e.phi);}
    for(std::size_t k=0;k<primary.geometryCells.size();++k){const auto q=mesh.cells(mesh.cellIndexToLabel[primary.geometryCells[k]]).centerCoord();values.emplace(LocalPointKey{q.x,q.y},centers[k]);}
    IncrementGradientProjection projected(support,old,config,[&](Vector2D q){return values.at({q.x,q.y});});
    for(auto& v:primary.vertices){const auto jet=projected.at(point(v.vertex));v.dx=jet.gradient.x;v.dy=jet.gradient.y;}
    for(auto& e:primary.edges){const auto a=point(e.vertices[0]),b=point(e.vertices[1]);const auto jet=projected.at(a+(b-a)*.5);e.dn=a.y==b.y?jet.gradient.y:jet.gradient.x;}
    report.evaluatedPointQueries+=projected.pointQueries();
}

template<class GeometryResult,class Finalize>
FirstJetPrepared<GeometryResult> prepareResolvedFirstJet(const QuadtreeMesh& mesh,
    const HermiteGalsPublication& source,const HermiteGalsResolvedMotion& motion,
    const HermiteGalsBoundaryFirstJetQuery& boundary,const HermiteGalsPublicationConfig& config,HermiteGalsBoundaryTreatment treatment,
    HermiteGalsTimeProjection projection,Finalize finalize,const char* entryName,
    HermiteGalsPointExecution execution=HermiteGalsPointExecution::Serial,const HermiteGalsPointSamplingConfig& sampling={})
{
    using Clock=std::chrono::steady_clock;const auto start=Clock::now();auto stamp=start;
    const auto lap=[&](){const auto now=Clock::now();const auto seconds=std::chrono::duration<double>(now-stamp).count();stamp=now;return seconds;};
    HermiteGalsIdentity id;const char* stage="input";
    try {
        id=source.state().identity();source.validateMesh(mesh);
        const auto& old=source.state();const auto& scalar=old.geometry().snapshot();
        const auto sourceTime=old.timeLayer(),targetTime=motion.targetTime();
        requireFirstJet(projection==HermiteGalsTimeProjection::CharacteristicPullback||
            projection==HermiteGalsTimeProjection::ClosedProjectionCompensated,"unknown temporal projection rule");
        const bool compensated=projection==HermiteGalsTimeProjection::ClosedProjectionCompensated;
        validateGradientUpdate(config,scalar,compensated,bool(boundary),treatment);
        if(config.gradientUpdate==HermiteGalsGradientUpdate::CompatibleIncrement)
            requireFirstJet(old.primary().measure==HermiteGalsMeasure::CartesianArea,"compatible increment requires uniform XY");
        requireFirstJet(execution==HermiteGalsPointExecution::Serial||execution==HermiteGalsPointExecution::Parallel,"unknown point execution rule");
        requireFirstJet(sampling.mode==HermiteGalsPointSampling::Exact||sampling.mode==HermiteGalsPointSampling::Numerical,"unknown transport point sampling rule");
        if(sampling.mode==HermiteGalsPointSampling::Numerical){
            requireFirstJet(std::isfinite(sampling.valueAbsoluteTolerance)&&sampling.valueAbsoluteTolerance>0
                &&std::isfinite(sampling.gradientAbsoluteTolerance)&&sampling.gradientAbsoluteTolerance>0,"invalid numerical transport point budgets");
            requireFirstJet(!boundary&&noQueryBoundaryTreatment(treatment),"numerical full sample requires no external boundary callback");
        }
        if(execution==HermiteGalsPointExecution::Parallel)
            requireFirstJet(!boundary&&noQueryBoundaryTreatment(treatment),"point-parallel sample requires no external boundary callback");
        bool parallel=false;
#ifdef _OPENMP
        parallel=execution==HermiteGalsPointExecution::Parallel&&!omp_in_parallel();
#endif
        requireFirstJet(!compensated||(!boundary&&noQueryBoundaryTreatment(treatment)),
            "closed projection requires no prescribed boundary query or treatment");
        requireFirstJet(noQueryBoundaryTreatment(treatment)||
            (treatment==HermiteGalsBoundaryTreatment::PrescribedOuterFirstJet&&bool(boundary)),"invalid or missing prescribed outer boundary data");
        using Position=HermiteGalsTimeLayer::Position;
        requireFirstJet(sourceTime.position==Position::Initial || sourceTime.position==Position::End,
            "source must be accepted Initial/End");
        requireFirstJet(targetTime.position==Position::Midpoint || targetTime.position==Position::End,
            "target must be Midpoint/End");
        requireFirstJet(sourceTime.step<std::numeric_limits<std::int64_t>::max() && targetTime.step==sourceTime.step+1,
            "target must be the next step");
        requireFirstJet(old.primary().transport==HermiteGalsTransport::FirstJet && old.primary().phaseVolumes.empty(),
            "requires FirstJet source without inventory");
        requireFirstJet(old.primary().geometryCells.size()==static_cast<std::size_t>(mesh.cellIndexNum()),
            "requires full scalar support");
        requireFirstJet(scalar.usesNativeLeaves() && !scalar.hasCornerSectors()
            && (scalar.reconstruction()==LevelSetHermiteGeometry::Reconstruction::EnrichedQuartic||scalar.usesC2Macro()),
            "requires native quartic source");
        requireFirstJet(scalar.usesC2Macro()==(config.scalarRepresentation==HermiteGalsScalarRepresentation::C2MacroQuartic),
            "source and arrival scalar representations differ");
        const auto bound=motion.sourceTime();
        requireFirstJet(motion.sourceIdentity()==id&&bound.step==sourceTime.step&&bound.begin==sourceTime.begin
            &&bound.end==sourceTime.end&&bound.position==sourceTime.position,"resolved motion source identity/clock differs");
        auto primary=old.primary();HermiteGalsFirstJetStepReport report;
        report.sourceCells=primary.geometryCells.size();report.sourceNativeFaces=source.domain().nativeFaceCount();
        report.sourceCutCells=source.volumes().report().cutCells;report.sourceGraphs=source.report().graphs;
        report.pointExecution=parallel?HermiteGalsPointExecution::Parallel:HermiteGalsPointExecution::Serial;report.pointSampling=sampling.mode;
        report.timings.input=lap();
        const auto point=[&](int v){return mesh.verts(mesh.vertIndexToLabel[v]).coord();};
        validateSymmetryFirstJet(mesh.minCoord,mesh.maxCoord,old.primary(),point,boundary,config,treatment);
        // Only the compensated path reuses first queries. The callable is
        // bound to this immutable original motion and one target clock.
        LocalPointStats serialStats;
        std::map<std::pair<double,double>,FirstPointRecord> firstQueries;
        const auto firstQuery=[&](Vector2D q){
            if(!compensated){++serialStats.traceQueries;return motion.trace(q);}
            const auto key=std::make_pair(q.x,q.y);const auto found=firstQueries.find(key);
            if(found!=firstQueries.end())return found->second.trace;
            ++serialStats.traceQueries;const auto trace=motion.trace(q);firstQueries.emplace(key,FirstPointRecord{trace,{}});return trace;
        };
        const auto preparedPoints=[&](const std::vector<Vector2D>& points,const auto& operation){
            using Value=decltype(operation(std::declval<const LevelSetBoundaryTrace::JetResult&>(),
                std::declval<LocalPointStats&>(),std::declval<std::optional<SourceJetSlot>&>(),std::declval<Vector2D>()));
            FullPointResults<Value> result(points.size());
            // The caller alone mutates this cache, between batches. Workers
            // read the original immutable motion and retained first-map values.
            const auto& retained=firstQueries;
            localPointBatch(points.size(),result.failures,report,[&](std::size_t k,LocalPointStats& stats){
                const auto key=std::make_pair(points[k].x,points[k].y);const auto found=retained.find(key);
                if(compensated&&found!=retained.end()){
                    result.traces[k]=found->second.trace;result.sourceJets[k]=found->second.sourceJet;
                } else {++stats.traceQueries;result.traces[k]=motion.trace(points[k]);}
                result.values[k]=operation(*result.traces[k],stats,result.sourceJets[k],points[k]);
            });
            return result;
        };
        const auto consumePoint=[&](const auto& batch,std::size_t k,Vector2D q){
            const auto jet=batch.at(k);
            if(compensated){
                const auto key=std::make_pair(q.x,q.y);
                auto record=firstQueries.emplace(key,FirstPointRecord{*batch.traces[k],batch.sourceJets[k]});
                if(batch.sourceJets[k])record.first->second.sourceJet=batch.sourceJets[k];
            }
            return jet;
        };
        const auto serialPoint=[&](Vector2D q,const auto& operation){
            const auto begin=Clock::now();++report.pointJobs;
            const auto trace=firstQuery(q);std::optional<SourceJetSlot> temporary;
            auto& slot=compensated?firstQueries.at(std::make_pair(q.x,q.y)).sourceJet:temporary;
            const auto jet=operation(trace,serialStats,slot,q);
            serialStats.seconds+=std::chrono::duration<double>(Clock::now()-begin).count();return jet;
        };
        const auto prescribed=[&](Vector2D q)->std::optional<HermiteGalsBoundaryFirstJet> {
            const bool outer=q.x==mesh.minCoord.x || q.x==mesh.maxCoord.x || q.y==mesh.minCoord.y || q.y==mesh.maxCoord.y;
            if(!outer)return {};
            if(treatment==HermiteGalsBoundaryTreatment::PrescribedOuterFirstJet){
                const bool west=q.x==mesh.minCoord.x,east=q.x==mesh.maxCoord.x,south=q.y==mesh.minCoord.y,north=q.y==mesh.maxCoord.y;
                const auto side=west?(south?Boundary::southWest:north?Boundary::northWest:Boundary::west):
                    east?(south?Boundary::southEast:north?Boundary::northEast:Boundary::east):south?Boundary::south:Boundary::north;
                const auto jet=boundary(q,targetTime,side);
                requireFirstJet(std::isfinite(jet.phi)&&LevelSetBoundaryTrace::Detail::finite(jet.gradient),"nonfinite prescribed outer first jet");
                return jet;
            }
            const auto trace=firstQuery(q);
            if(!trace.trajectory.boundaryHit)return {};
            requireFirstJet(bool(boundary),"outer boundary hit requires prescribed first jet");
            requireFirstJet(trace.trajectory.side.has_value(),"boundary hit lacks side");
            const auto jet=boundary(q,targetTime,*trace.trajectory.side);
            requireFirstJet(std::isfinite(jet.phi) && std::isfinite(jet.gradient.x) && std::isfinite(jet.gradient.y),
                "nonfinite prescribed boundary first jet");
            return jet;
        };
        stage="outer boundary first jet";
        std::vector<int> vertices;std::vector<std::array<int,2>> edges;
        for(auto& v:primary.vertices){const auto jet=prescribed(point(v.vertex));
            if(jet){v.phi=jet->phi;v.dx=jet->gradient.x;v.dy=jet->gradient.y;++report.prescribedVertices;}
            else vertices.push_back(v.vertex);}
        for(auto& e:primary.edges){const auto a=point(e.vertices[0]),b=point(e.vertices[1]);const auto jet=prescribed(a+(b-a)*.5);
            if(jet){e.phi=jet->phi;e.dn=a.y==b.y?jet->gradient.y:jet->gradient.x;++report.prescribedEdges;}
            else edges.push_back(e.vertices);}
        report.timings.boundary=lap();
        const auto pulled=[&](const LevelSetBoundaryTrace::JetResult& trace,LocalPointStats& stats,std::optional<SourceJetSlot>& slot,Vector2D q){
            requireInteriorFirstJetTrace(trace,q);
            const auto jet=sourceTransportJet(scalar,trace.trajectory.point,sampling,stats,slot,compensated);
            const auto gradient=trace.derivative->transposeTimes({jet.dx,jet.dy});
            requireFirstJet(std::isfinite(jet.value)&&LevelSetBoundaryTrace::Detail::finite(gradient),"nonfinite pulled first jet");
            return HermiteGalsBoundaryFirstJet{jet.value,gradient};
        };
        stage="vertex first-jet transport";
        std::optional<FullPointResults<HermiteGalsBoundaryFirstJet>> vertexBatch;
        if(parallel){std::vector<Vector2D> points;points.reserve(vertices.size());for(int v:vertices)points.push_back(point(v));
            vertexBatch.emplace(preparedPoints(points,pulled));}
        for(std::size_t k=0;k<vertices.size();++k){const int v=vertices[k];
            const auto jet=applySymmetryFirstJet(point(v),(parallel?consumePoint(*vertexBatch,k,point(v)):serialPoint(point(v),pulled)),config,treatment);
            auto it=std::lower_bound(primary.vertices.begin(),primary.vertices.end(),v,[](const auto& a,int b){return a.vertex<b;});
            it->phi=jet.phi;it->dx=jet.gradient.x;it->dy=jet.gradient.y;
        }
        vertexBatch.reset();report.timings.vertices=lap();
        stage="edge first-jet transport";
        std::optional<FullPointResults<HermiteGalsBoundaryFirstJet>> edgeBatch;
        if(parallel){std::vector<Vector2D> points;points.reserve(edges.size());
            for(const auto& e:edges){const auto a=point(e[0]),b=point(e[1]);points.push_back(a+(b-a)*.5);}
            edgeBatch.emplace(preparedPoints(points,pulled));}
        for(std::size_t k=0;k<edges.size();++k){const auto& e=edges[k];const auto a=point(e[0]),b=point(e[1]),q=a+(b-a)*.5;
            const auto jet=applySymmetryFirstJet(q,(parallel?consumePoint(*edgeBatch,k,q):serialPoint(q,pulled)),config,treatment);
            auto it=std::lower_bound(primary.edges.begin(),primary.edges.end(),e,[](const auto& a,const auto& b){return a.vertices<b;});
            it->phi=jet.phi;it->dn=a.y==b.y?jet.gradient.y:jet.gradient.x;
        }
        edgeBatch.reset();report.timings.edges=lap();
        stage="center prediction";
        const auto centerValue=[&](const LevelSetBoundaryTrace::JetResult& trace,LocalPointStats& stats,std::optional<SourceJetSlot>& slot,Vector2D){
            requireFirstJet(!trace.trajectory.boundaryHit,"center prediction requires unresolved boundary data");
            return sourceTransportValue(scalar,trace.trajectory.point,sampling,stats,slot,compensated);
        };
        std::optional<FullPointResults<double>> centerBatch;
        if(parallel){std::vector<Vector2D> points;points.reserve(primary.geometryCells.size());
            for(int cell:primary.geometryCells)points.push_back(mesh.cells(mesh.cellIndexToLabel[cell]).centerCoord());
            centerBatch.emplace(preparedPoints(points,centerValue));}
        std::vector<double> centers;centers.reserve(primary.geometryCells.size());
        for(std::size_t k=0;k<primary.geometryCells.size();++k){const int cell=primary.geometryCells[k];
            const auto q=mesh.cells(mesh.cellIndexToLabel[cell]).centerCoord();
            centers.push_back(parallel?consumePoint(*centerBatch,k,q):serialPoint(q,centerValue));
        }
        centerBatch.reset();
        report.transportedVertices=vertices.size();report.transportedEdges=edges.size();report.tracedCenters=centers.size();
        primary.gradientUpdate=config.gradientUpdate;
        if(config.gradientUpdate==HermiteGalsGradientUpdate::CompatibleIncrement){
            stage="compatible gradient increment";projectFullGradient(mesh,primary,centers,scalar,config,report);
        }
        HermiteGalsCandidate candidate(mesh,targetTime,primary);
        report.timings.centers=lap();
        stage="mixed closure";
        const auto mixed=closeHermiteGalsMixedAmr(mesh,candidate,config.volume.symmetryPlanes);
        const HermiteGalsCenterPrediction center{candidate.identity(),centers};
        report.timings.mixed=lap();
        stage="unique publication";
        std::optional<GeometryResult> result;
        if(!compensated){
            result.emplace(finalize(mesh,std::move(candidate),mixed,center,config));
            report.timings.publication=lap();
        }
        if(compensated){
            // Q=P*T*E. This temporary E(Q y) is never an accepted source.
            // Q+.5*(Q*Q-P*T*T*E) cancels the O(dt^2) projection defect;
            // both T factors use the ORIGINAL motion identity and clock.
            // Use the same authenticated scalar builder as full publication.
            // This E(Q y) retains the selected scalar and actual centers;
            // it constructs no topology, volume or surface quadrature.
            const auto projected=prepareHermiteGalsFirstJetScalar(mesh,candidate,mixed,center,config);
            report.timings.publication+=lap();
            const auto corrected=[&](const LevelSetBoundaryTrace::JetResult& first,LocalPointStats& stats,std::optional<SourceJetSlot>& slot,Vector2D){
                requireFirstJet(!first.trajectory.boundaryHit&&first.derivative.has_value(),
                    "closed projection first query requires interior Jacobian");
                ++stats.traceQueries;const auto second=motion.trace(first.trajectory.point);
                requireFirstJet(!second.trajectory.boundaryHit&&second.derivative.has_value(),
                    "closed projection second query requires interior Jacobian");
                const auto one=sourceTransportJet(scalar,first.trajectory.point,sampling,stats,slot,true);
                const auto staged=localTransportJet(projected,first.trajectory.point,sampling,stats);
                const auto twice=sourceTransportJet(scalar,second.trajectory.point,sampling,stats,slot,true,false);
                const auto gtwice=second.derivative->transposeTimes({twice.dx,twice.dy});
                const Vector2D inner{one.dx+.5*(staged.dx-gtwice.x),one.dy+.5*(staged.dy-gtwice.y)};
                const auto gradient=first.derivative->transposeTimes(inner);
                const double value=one.value+.5*(staged.value-twice.value);
                requireFirstJet(std::isfinite(value)&&LevelSetBoundaryTrace::Detail::finite(gradient),"nonfinite compensated first jet");
                return HermiteGalsBoundaryFirstJet{value,gradient};
            };
            const auto correctedValue=[&](const LevelSetBoundaryTrace::JetResult& first,LocalPointStats& stats,std::optional<SourceJetSlot>& slot,Vector2D){
                requireFirstJet(!first.trajectory.boundaryHit&&first.derivative.has_value(),
                    "closed projection first query requires interior Jacobian");
                ++stats.traceQueries;const auto second=motion.trace(first.trajectory.point);
                requireFirstJet(!second.trajectory.boundaryHit&&second.derivative.has_value(),
                    "closed projection second query requires interior Jacobian");
                const double one=sourceTransportValue(scalar,first.trajectory.point,sampling,stats,slot,true),
                    staged=localTransportValue(projected,first.trajectory.point,sampling,stats),
                    twice=sourceTransportValue(scalar,second.trajectory.point,sampling,stats,slot,true,false);
                const double value=one+.5*(staged-twice);
                requireFirstJet(std::isfinite(value),"nonfinite compensated center value");return value;
            };
            auto correctedPrimary=old.primary();stage="closed projection vertex correction";
            std::optional<FullPointResults<HermiteGalsBoundaryFirstJet>> correctedVertexBatch;
            if(parallel){std::vector<Vector2D> points;points.reserve(correctedPrimary.vertices.size());
                for(const auto& v:correctedPrimary.vertices)points.push_back(point(v.vertex));
                correctedVertexBatch.emplace(preparedPoints(points,corrected));}
            for(std::size_t k=0;k<correctedPrimary.vertices.size();++k){auto& v=correctedPrimary.vertices[k];const auto q=point(v.vertex);
                const auto jet=applySymmetryFirstJet(q,(parallel?consumePoint(*correctedVertexBatch,k,q):serialPoint(q,corrected)),config,treatment);
                v.phi=jet.phi;v.dx=jet.gradient.x;v.dy=jet.gradient.y;}
            correctedVertexBatch.reset();report.timings.vertices+=lap();stage="closed projection edge correction";
            std::optional<FullPointResults<HermiteGalsBoundaryFirstJet>> correctedEdgeBatch;
            if(parallel){std::vector<Vector2D> points;points.reserve(correctedPrimary.edges.size());
                for(const auto& e:correctedPrimary.edges){const auto a=point(e.vertices[0]),b=point(e.vertices[1]);points.push_back(a+(b-a)*.5);}
                correctedEdgeBatch.emplace(preparedPoints(points,corrected));}
            for(std::size_t k=0;k<correctedPrimary.edges.size();++k){auto& e=correctedPrimary.edges[k];
                const auto a=point(e.vertices[0]),b=point(e.vertices[1]),q=a+(b-a)*.5;
                const auto jet=applySymmetryFirstJet(q,(parallel?consumePoint(*correctedEdgeBatch,k,q):serialPoint(q,corrected)),config,treatment);
                e.phi=jet.phi;e.dn=a.y==b.y?jet.gradient.y:jet.gradient.x;}
            correctedEdgeBatch.reset();report.timings.edges+=lap();stage="closed projection center correction";
            std::optional<FullPointResults<double>> correctedCenterBatch;
            if(parallel){std::vector<Vector2D> points;points.reserve(correctedPrimary.geometryCells.size());
                for(int cell:correctedPrimary.geometryCells)points.push_back(mesh.cells(mesh.cellIndexToLabel[cell]).centerCoord());
                correctedCenterBatch.emplace(preparedPoints(points,correctedValue));}
            std::vector<double> correctedCenters;correctedCenters.reserve(correctedPrimary.geometryCells.size());
            for(std::size_t k=0;k<correctedPrimary.geometryCells.size();++k){const int cell=correctedPrimary.geometryCells[k];
                const auto q=mesh.cells(mesh.cellIndexToLabel[cell]).centerCoord();
                correctedCenters.push_back(parallel?consumePoint(*correctedCenterBatch,k,q):serialPoint(q,correctedValue));}
            correctedCenterBatch.reset();
            report.timings.centers+=lap();HermiteGalsCandidate finalCandidate(mesh,targetTime,correctedPrimary);
            stage="closed projection mixed closure";
            const auto finalMixed=closeHermiteGalsMixedAmr(mesh,finalCandidate,config.volume.symmetryPlanes);report.timings.mixed+=lap();
            stage="closed projection final publication";
            const HermiteGalsCenterPrediction finalCenters{finalCandidate.identity(),correctedCenters};
            result.emplace(finalize(mesh,std::move(finalCandidate),finalMixed,finalCenters,config));
            report.timings.publication+=lap();
        }
        mergeLocalPointStats(report,serialStats);
        report.timings.total=std::chrono::duration<double>(Clock::now()-start).count();
        return {id,motion.velocityIdentity(),std::move(*result),report};
    } catch(const std::exception& e) {
        throw std::runtime_error(std::string("[Interface/")+entryName+" source="+std::to_string(id.value())
            +" arrival="+std::to_string(motion.velocityIdentity().value())+" stage="+stage+"] "+e.what());
    }
}

// This finalizer preserves the FirstJet publication's input/closure/patch
// arithmetic, then stops at the structural scalar binding. No publication
// certification is fabricated for this temporary rate-only consumer.
HermiteGalsState bindFirstJetSample(const QuadtreeMesh& mesh,HermiteGalsCandidate&& candidate,
    const HermiteGalsAmrClosureResult& mixed,const HermiteGalsCenterPrediction& centers,
    const HermiteGalsPublicationConfig& config)
{
    const auto id=candidate.identity();
    candidate.validateMesh(mesh);
    requireFirstJet(candidate.primary().transport==HermiteGalsTransport::FirstJet && candidate.primary().phaseVolumes.empty(),
        "FirstJet sample requires state without inventory");
    requireFirstJet(!candidate.hasGeometry(),"FirstJet sample requires unpublished candidate");
    requireFirstJet(mixed.primaryIdentity()==id && centers.primaryIdentity==id,"FirstJet primary revision mismatch");
    mixed.skeleton().validateMesh(mesh);
    const auto time=candidate.timeLayer(),mt=mixed.timeLayer();
    requireFirstJet(time.step==mt.step && time.begin==mt.begin && time.end==mt.end && time.position==mt.position,
        "FirstJet mixed closure time mismatch");
    requireFirstJet(mixed.vertexMixedDerivatives().size()==candidate.primary().vertices.size()
        && centers.values.size()==candidate.primary().geometryCells.size(),"FirstJet closure layout mismatch");
    for(double value:centers.values)requireFirstJet(std::isfinite(value),"FirstJet nonfinite center prediction");
    const auto tolerance=[](double value){return std::isfinite(value) && value>0 && value<1;};
    requireFirstJet(tolerance(config.relativeInterfaceResidualTolerance) && tolerance(config.unitNormalTolerance),
        "invalid publication tolerance");
    const HermiteGalsClosure closure{mixed.vertexMixedDerivatives(),centers.values};
    auto scalar=prepareHermiteGalsFirstJetScalar(mesh,candidate,mixed,centers,config);
    candidate.bindGeometry(mesh,id,closure,std::move(scalar));
    return std::move(candidate).finish();
}
} // namespace

HermiteGalsFirstJetStep prepareHermiteGalsFirstJetStep(const QuadtreeMesh& mesh,
    const HermiteGalsPublication& source,const HermiteGalsResolvedMotion& motion,
    const HermiteGalsBoundaryFirstJetQuery& boundary,const HermiteGalsPublicationConfig& config,HermiteGalsBoundaryTreatment treatment,
    HermiteGalsTimeProjection projection,HermiteGalsPointExecution execution,const HermiteGalsPointSamplingConfig& sampling)
{
    auto prepared=prepareResolvedFirstJet<HermiteGalsPublication>(mesh,source,motion,boundary,config,treatment,projection,
        [](const QuadtreeMesh& m,HermiteGalsCandidate&& candidate,const HermiteGalsAmrClosureResult& mixed,
            const HermiteGalsCenterPrediction& centers,const HermiteGalsPublicationConfig& controls) {
            return publishHermiteGalsFirstJetGeometry(m,std::move(candidate),mixed,centers,controls);
        },"HermiteGalsFirstJetStep",execution,sampling);
    return HermiteGalsFirstJetStep(std::make_shared<const HermiteGalsFirstJetStep::Data>(
        prepared.source,prepared.arrival,std::move(prepared.result),prepared.report));
}

HermiteGalsFirstJetSample prepareHermiteGalsFirstJetSample(const QuadtreeMesh& mesh,
    const HermiteGalsPublication& source,const HermiteGalsResolvedMotion& motion,
    const HermiteGalsBoundaryFirstJetQuery& boundary,const HermiteGalsPublicationConfig& config,HermiteGalsBoundaryTreatment treatment,
    HermiteGalsTimeProjection projection)
{
    auto prepared=prepareResolvedFirstJet<HermiteGalsState>(mesh,source,motion,boundary,config,treatment,projection,
        bindFirstJetSample,"HermiteGalsFirstJetSample");
    return HermiteGalsFirstJetSample(std::make_shared<const HermiteGalsFirstJetSample::Data>(
        prepared.source,prepared.arrival,std::move(prepared.result),prepared.report));
}

namespace
{
// Logical supports only. All data below are queried from the actual source/map;
// the output mask never pretends that a dependency halo is a published field.
// A sample's exact native support and point inventory are built once.
// All sets retain the original lexical order; no scalar result is cached here.
struct LocalSampleSupport
{
    std::set<int> vertices;
    std::set<std::array<int,2>> edges;
    std::set<LocalPointKey> jetPoints,points;
};
struct IncrementGridCache {std::once_flag once;std::optional<IncrementGradientGrid> grid;};
struct LocalUniformSupport
{
    std::shared_ptr<IncrementGridCache> incrementCache=std::make_shared<IncrementGridCache>();
    UniformLeafGrid grid;
    std::vector<std::array<int,2>> cellPosition,vertexPosition;
    std::vector<int> gridVertexIndices,gridCellIndices;
    std::vector<Vector2D> nativePoints,cellCenters;
    std::vector<unsigned char> independent;
    explicit LocalUniformSupport(const QuadtreeMesh& m):grid(m,CartesianGridScope::UniformLeaves),
        cellPosition(m.cellIndexNum()),vertexPosition(m.vertIndexNum()),
        gridVertexIndices(grid.vertices().size()),nativePoints(m.vertIndexNum()),
        cellCenters(m.cellIndexNum()),independent(m.vertIndexNum())
    {
        for(int v=0;v<m.vertIndexNum();++v){
            const auto& vertex=m.verts(m.vertIndexToLabel[v]);
            nativePoints[v]=vertex.coord();independent[v]=vertex.type()==TVertType::non;
        }
        for(int c=0;c<m.cellIndexNum();++c)cellCenters[c]=m.cells(m.cellIndexToLabel[c]).centerCoord();
        for(std::size_t k=0;k<grid.vertices().size();++k)gridVertexIndices[k]=m.verts(grid.vertices()[k]).index();
        for(int label:grid.cells())gridCellIndices.push_back(m.cells(label).index());
        for(int j=0;j<grid.ny();++j)for(int i=0;i<grid.nx();++i)
            cellPosition[m.cells(grid.cells()[j*grid.nx()+i]).index()]={{i,j}};
        for(int j=0;j<=grid.ny();++j)for(int i=0;i<=grid.nx();++i)vertexPosition[idAt(i,j)]={{i,j}};
    }
    const IncrementGradientGrid& projectionGrid() const
    {
        std::call_once(incrementCache->once,[&]{incrementCache->grid.emplace(grid.nx(),grid.ny(),
            [&](int i,int j){return point(idAt(i,j));},
            [&](int i,int j){return cellCenters[gridCellIndices[j*grid.nx()+i]];});});
        return *incrementCache->grid;
    }
    int idAt(int i,int j) const{return gridVertexIndices[j*(grid.nx()+1)+i];}
    Vector2D point(int v) const{return nativePoints[v];}
    Vector2D sidePoint(int i,int j,int axis,int sample) const
    {
        const auto a=point(axis==0?idAt(sample/2,j):idAt(i,sample/2));
        if(sample%2==0)return a;
        const auto b=point(axis==0?idAt(sample/2+1,j):idAt(i,sample/2+1));
        return a+(b-a)*.5;
    }
    std::set<int> vertices(const std::vector<int>& cells) const
    {
        std::set<int> result;
        for(int c:cells){const auto ij=cellPosition[c];
            for(int dj=0;dj<2;++dj)for(int di=0;di<2;++di)result.insert(idAt(ij[0]+di,ij[1]+dj));}
        return result;
    }
    std::set<std::array<int,2>> edges(const std::vector<int>& cells) const
    {
        std::set<std::array<int,2>> result;
        for(int c:cells){const auto ij=cellPosition[c];const int i=ij[0],j=ij[1];
            result.insert({{idAt(i,j),idAt(i+1,j)}});result.insert({{idAt(i,j+1),idAt(i+1,j+1)}});
            result.insert({{idAt(i,j),idAt(i,j+1)}});result.insert({{idAt(i+1,j),idAt(i+1,j+1)}});}
        return result;
    }
    LocalSampleSupport sampleSupport(const std::vector<int>& cells,const HermiteGalsPublicationConfig& config={}) const
    {
        LocalSampleSupport result;result.vertices=vertices(cells);result.edges=edges(cells);
        const auto add=[&](Vector2D p){result.jetPoints.emplace(p.x,p.y);};
        for(int v:result.vertices){add(point(v));const auto ij=vertexPosition[v];
            for(int axis=0;axis<2;++axis){const int n=axis==0?grid.nx():grid.ny(),t=2*ij[axis];
                const int begin=std::clamp(t-2,0,2*n-4);
                for(int k=begin;k<begin+5;++k)add(sidePoint(ij[0],ij[1],axis,k));}}
        for(const auto& e:result.edges){const auto a=point(e[0]),b=point(e[1]);add(a+(b-a)*.5);}
        result.points=result.jetPoints;
        for(int c:cells){const auto p=cellCenters[c];result.points.emplace(p.x,p.y);}
        if(config.gradientUpdate==HermiteGalsGradientUpdate::CompatibleIncrement)projectionGrid().expand(result.jetPoints,result.points,config);
        return result;
    }
};
struct LocalFirstJetPayload
{
    HermiteGalsPrimaryState primary;
    HermiteGalsClosure closure;
    LevelSetHermiteGeometry::Snapshot::NativeLeafSubset patch;
};
template<class Jet,class Center>
LocalFirstJetPayload localFirstJetPayload(const LocalUniformSupport& support,const std::vector<int>& cells,
    const LocalSampleSupport& sampleSupport,const Jet& jetAt,const Center& centerAt,const HermiteGalsPublicationConfig& config)
{
    LocalFirstJetPayload result;
    result.primary.measure=HermiteGalsMeasure::CartesianArea;result.primary.transport=HermiteGalsTransport::FirstJet;
    result.primary.geometryCells=cells;result.primary.gradientUpdate=config.gradientUpdate;
    const auto coordinate=[&](Vector2D q,int axis){return axis==0?q.x:q.y;};
    for(std::size_t k=0;k<config.volume.symmetryPlanes.size();++k){const auto& plane=config.volume.symmetryPlanes[k];
        requireFirstJet((plane.normalAxis==0||plane.normalAxis==1)&&std::isfinite(plane.coordinate),"invalid symmetry plane");
        requireFirstJet(plane.coordinate==coordinate(support.grid.minCoord(),plane.normalAxis)
            ||plane.coordinate==coordinate(support.grid.maxCoord(),plane.normalAxis),"symmetry plane must be an actual outer domain section");
        for(std::size_t j=0;j<k;++j)requireFirstJet(!(plane==config.volume.symmetryPlanes[j]),"duplicate symmetry plane");}
    for(int v:sampleSupport.vertices) {
        const auto p=support.point(v);const auto jet=jetAt(p);const auto ij=support.vertexPosition[v];
        result.primary.vertices.push_back({v,jet.phi,jet.gradient.x,jet.gradient.y});
        bool reflected=false;
        for(const auto& plane:config.volume.symmetryPlanes)if(coordinate(p,plane.normalAxis)==plane.coordinate){
            requireFirstJet((plane.normalAxis==0?jet.gradient.x:jet.gradient.y)==0,
                "symmetry requires exactly zero primary vertex normal derivative");reflected=true;}
        const auto derivative=[&](int axis){const int n=axis==0?support.grid.nx():support.grid.ny(),target=2*ij[axis];
            const auto at=[&](int k){return support.sidePoint(ij[0],ij[1],axis,k);};
            const auto coord=[&](int k){return coordinate(at(k),axis);};
            const auto value=[&](int k){const auto sample=jetAt(at(k));return axis==0?sample.gradient.y:sample.gradient.x;};
            bool zero=false;for(const auto& plane:config.volume.symmetryPlanes)
                if(plane.normalAxis==axis&&coordinate(p,axis)==plane.coordinate)zero=true;
            return zero?0:LAMP::Geometry::compactQuarticDerivativeAtNode(target,2*n+1,coord,value,requireFirstJet);};
        const double dx=derivative(0),dy=derivative(1),expected=.5*dx+.5*dy;
        const double mixed=reflected?0:expected;
        requireFirstJet(std::isfinite(expected)&&std::isfinite(mixed),"nonfinite compact mixed derivative");
        const double scale=std::max({1.0,std::abs(dx),std::abs(dy),std::abs(mixed)});
        requireFirstJet(std::abs(mixed/scale-expected/scale)<=1e-12,"direct compact/constraint residual exceeds tolerance");
        result.closure.vertexMixedDerivatives.push_back(mixed);
    }
    for(const auto& e:sampleSupport.edges){const auto a=support.point(e[0]),b=support.point(e[1]);
        const auto jet=jetAt(a+(b-a)*.5);const bool horizontal=a.y==b.y;const double dn=horizontal?jet.gradient.y:jet.gradient.x;
        for(const auto& plane:config.volume.symmetryPlanes)if(coordinate(a,plane.normalAxis)==plane.coordinate
            &&coordinate(b,plane.normalAxis)==plane.coordinate)
            requireFirstJet(dn==0,"symmetry requires exactly zero primary edge normal derivative");
        result.primary.edges.push_back({e,jet.phi,dn});}
    for(int c:cells){const double value=centerAt(support.cellCenters[c]);
        requireFirstJet(std::isfinite(value),"FirstJet nonfinite center prediction");result.closure.cellBubbleCoefficients.push_back(value);}
    // Retain only the actually predicted native support. The immutable source
    // snapshot supplies coordinates/topology/locator; no unused full-mesh
    // payload or zero-filled geometry is synthesized at a local time node.
    result.patch.nativeCells=cells;
    result.patch.cellCenters=result.closure.cellBubbleCoefficients;
    result.patch.nativeVertices.reserve(result.primary.vertices.size());
    result.patch.nodes.reserve(result.primary.vertices.size());
    for(std::size_t k=0;k<result.primary.vertices.size();++k){const auto& v=result.primary.vertices[k];
        result.patch.nativeVertices.push_back(v.vertex);
        result.patch.nodes.push_back({v.phi,v.dx,v.dy,result.closure.vertexMixedDerivatives[k]});}
    for(const auto& e:result.primary.edges)result.patch.edgeMidpoints.push_back({e.vertices,e.phi,e.dn});
    return result;
}


template<class Jet,class Center>
LocalFirstJetPayload updatedLocalPayload(const LocalUniformSupport& support,const std::vector<int>& cells,
    const LocalSampleSupport& selected,const Jet& jet,const Center& value,const HermiteGalsPublicationConfig& config,
    const LevelSetHermiteGeometry::Snapshot& scalar,HermiteGalsFirstJetStepReport& report)
{
    if(config.gradientUpdate==HermiteGalsGradientUpdate::CharacteristicPullback)
        return localFirstJetPayload(support,cells,selected,jet,value,config);
    IncrementGradientProjection projected(support.projectionGrid(),scalar,config,value);
    auto result=localFirstJetPayload(support,cells,selected,[&](Vector2D q){return projected.at(q);},value,config);
    report.evaluatedPointQueries+=projected.pointQueries();return result;
}

LocalPointJets transportLocalPointJets(const LocalUniformSupport&,const LocalSampleSupport&,
    const LevelSetHermiteGeometry::Snapshot&,const HermiteGalsResolvedMotion&,const HermiteGalsPublicationConfig&,
    bool,HermiteGalsPointExecution,const HermiteGalsPointSamplingConfig&,HermiteGalsFirstJetStepReport&,HermiteGalsBoundaryTreatment);
LocalFirstJetPayload pointBatchLocalPayload(const LocalUniformSupport& support,const std::vector<int>& cells,
    const LevelSetHermiteGeometry::Snapshot& scalar,const HermiteGalsResolvedMotion& motion,
    const HermiteGalsPublicationConfig& config,bool compensated,const HermiteGalsPointSamplingConfig& sampling,
    HermiteGalsFirstJetStepReport& report,HermiteGalsBoundaryTreatment treatment)
{
    const auto requestedSupport=support.sampleSupport(cells,config);
    const auto values=transportLocalPointJets(support,requestedSupport,scalar,motion,config,compensated,
        HermiteGalsPointExecution::Parallel,sampling,report,treatment);
    return updatedLocalPayload(support,cells,requestedSupport,[&](Vector2D q){return values.at(q);},
        [&](Vector2D q){return values.valueAt(q);},config,scalar,report);
}
}

namespace
{
struct LocalSampleProducts {HermiteGalsState state;HermiteGalsFirstJetStepReport report;};
template<class Binder>
LocalSampleProducts localFirstJetKernel(const HermiteGalsPublication& source,const LocalUniformSupport& support,
    const HermiteGalsResolvedMotion& motion,const std::vector<int>& requestedNativeCells,
    const HermiteGalsBoundaryFirstJetQuery& boundary,const HermiteGalsPublicationConfig& config,
    HermiteGalsBoundaryTreatment treatment,HermiteGalsTimeProjection projection,HermiteGalsPointExecution execution,
    const HermiteGalsPointSamplingConfig& sampling,const Binder& bind,std::chrono::steady_clock::time_point start)
{
    using Clock=std::chrono::steady_clock;HermiteGalsIdentity id;const char* stage="input";
    try {
        LAMP::Geometry::Detail::ExactExpansion::validateActualEnvironment();id=source.state().identity();
        requireFirstJet(!requestedNativeCells.empty(),"empty requested native cells");
        auto cells=requestedNativeCells;std::sort(cells.begin(),cells.end());
        for(std::size_t k=0;k<cells.size();++k){requireFirstJet(cells[k]>=0&&static_cast<std::size_t>(cells[k])<support.cellPosition.size(),"invalid requested native cell");
            requireFirstJet(k==0||cells[k]!=cells[k-1],"duplicate requested native cell");}
        const auto sourceTime=source.state().timeLayer(),targetTime=motion.targetTime(),bound=motion.sourceTime();
        using Position=HermiteGalsTimeLayer::Position;
        requireFirstJet(sourceTime.position==Position::Initial||sourceTime.position==Position::End,"source must be accepted Initial/End");
        requireFirstJet(targetTime.position==Position::Midpoint||targetTime.position==Position::End,"target must be Midpoint/End");
        requireFirstJet(sourceTime.step<std::numeric_limits<std::int64_t>::max()&&targetTime.step==sourceTime.step+1,"target must be the next step");
        requireFirstJet(motion.sourceIdentity()==id&&bound.step==sourceTime.step&&bound.begin==sourceTime.begin
            &&bound.end==sourceTime.end&&bound.position==sourceTime.position,"resolved motion source identity/clock differs");
        requireFirstJet(projection==HermiteGalsTimeProjection::CharacteristicPullback||projection==HermiteGalsTimeProjection::ClosedProjectionCompensated,
            "unknown temporal projection rule");
        const bool compensated=projection==HermiteGalsTimeProjection::ClosedProjectionCompensated;
        validateSymmetryFirstJet(support.grid.minCoord(),support.grid.maxCoord(),source.state().primary(),
            [&](int vertex){return support.point(vertex);},boundary,config,treatment);
        requireFirstJet(execution==HermiteGalsPointExecution::Serial||execution==HermiteGalsPointExecution::Parallel,"unknown point execution rule");
        requireFirstJet(sampling.mode==HermiteGalsPointSampling::Exact||sampling.mode==HermiteGalsPointSampling::Numerical,"unknown transport point sampling rule");
        if(sampling.mode==HermiteGalsPointSampling::Numerical)
            requireFirstJet(std::isfinite(sampling.valueAbsoluteTolerance)&&sampling.valueAbsoluteTolerance>0
                &&std::isfinite(sampling.gradientAbsoluteTolerance)&&sampling.gradientAbsoluteTolerance>0,"invalid numerical transport point budgets");
        if(execution==HermiteGalsPointExecution::Parallel){
            requireFirstJet(!boundary&&noQueryBoundaryTreatment(treatment),"point-parallel sample requires no external boundary callback");
#ifdef _OPENMP
            requireFirstJet(!omp_in_parallel(),"point-parallel sample cannot run within an outer OpenMP team");
#endif
        }
        requireFirstJet(!compensated||(!boundary&&noQueryBoundaryTreatment(treatment)),
            "closed projection requires no prescribed boundary query or treatment");
        requireFirstJet(noQueryBoundaryTreatment(treatment)
            ||(treatment==HermiteGalsBoundaryTreatment::PrescribedOuterFirstJet&&bool(boundary)),"invalid or missing prescribed outer boundary data");
        const auto tolerance=[](double value){return std::isfinite(value)&&value>0&&value<1;};
        requireFirstJet(tolerance(config.relativeInterfaceResidualTolerance)&&tolerance(config.unitNormalTolerance),"invalid publication tolerance");
        const auto& scalar=source.state().geometry().snapshot();
        validateGradientUpdate(config,scalar,compensated,bool(boundary),treatment);HermiteGalsFirstJetStepReport report;
        report.sourceCells=source.state().primary().geometryCells.size();report.sourceNativeFaces=source.domain().nativeFaceCount();
        report.sourceCutCells=source.volumes().report().cutCells;report.sourceGraphs=source.report().graphs;
        report.requestedNativeCells=cells.size();report.evaluatedNativeCells=cells.size();
        report.pointExecution=execution;report.pointSampling=sampling.mode;
        report.timings.input=std::chrono::duration<double>(Clock::now()-start).count();
        stage="local mixed/departure closure";const auto kernelStart=Clock::now();std::optional<LocalFirstJetPayload> payload;
        if(execution==HermiteGalsPointExecution::Parallel)
            payload.emplace(pointBatchLocalPayload(support,cells,scalar,motion,config,compensated,sampling,report,treatment));
        else {
        const auto requestedSupport=support.sampleSupport(cells,config);
        LocalPointStats serialStats;
        std::map<std::pair<double,double>,FirstPointRecord> firstQueries;
        const auto firstQuery=[&](Vector2D q){
            if(!compensated){++serialStats.traceQueries;return motion.trace(q);}
            const auto key=std::make_pair(q.x,q.y);const auto found=firstQueries.find(key);if(found!=firstQueries.end())return found->second.trace;
            ++serialStats.traceQueries;const auto value=motion.trace(q);firstQueries.emplace(key,FirstPointRecord{value,{}});return value;};
        const auto originalJet=[&](Vector2D q,const LevelSetBoundaryTrace::JetResult& trace){
            std::optional<SourceJetSlot> temporary;
            auto& slot=compensated?firstQueries.at(std::make_pair(q.x,q.y)).sourceJet:temporary;
            return sourceTransportJet(scalar,trace.trajectory.point,sampling,serialStats,slot,compensated);
        };
        const auto originalValue=[&](Vector2D q,const LevelSetBoundaryTrace::JetResult& trace){
            std::optional<SourceJetSlot> temporary;
            auto& slot=compensated?firstQueries.at(std::make_pair(q.x,q.y)).sourceJet:temporary;
            return sourceTransportValue(scalar,trace.trajectory.point,sampling,serialStats,slot,compensated);
        };
        std::map<std::pair<double,double>,HermiteGalsBoundaryFirstJet> rawJets;
        const auto pulled=[&](Vector2D q){const auto key=std::make_pair(q.x,q.y);const auto found=rawJets.find(key);
            if(found!=rawJets.end())return found->second;
            std::optional<HermiteGalsBoundaryFirstJet> prescribed;
            const bool outer=q.x==support.grid.minCoord().x||q.x==support.grid.maxCoord().x||q.y==support.grid.minCoord().y||q.y==support.grid.maxCoord().y;
            if(outer&&treatment==HermiteGalsBoundaryTreatment::PrescribedOuterFirstJet){
                const bool west=q.x==support.grid.minCoord().x,east=q.x==support.grid.maxCoord().x,south=q.y==support.grid.minCoord().y,north=q.y==support.grid.maxCoord().y;
                const auto side=west?(south?Boundary::southWest:north?Boundary::northWest:Boundary::west):
                    east?(south?Boundary::southEast:north?Boundary::northEast:Boundary::east):south?Boundary::south:Boundary::north;
                prescribed=boundary(q,targetTime,side);
            } else if(outer){const auto trace=firstQuery(q);if(trace.trajectory.boundaryHit){
                requireFirstJet(bool(boundary),"outer boundary hit requires prescribed first jet");
                requireFirstJet(trace.trajectory.side.has_value(),"boundary hit lacks side");prescribed=boundary(q,targetTime,*trace.trajectory.side);}}
            HermiteGalsBoundaryFirstJet jet;
            if(prescribed)jet=*prescribed;
            else {const auto trace=firstQuery(q);requireInteriorFirstJetTrace(trace,q);
                const auto old=originalJet(q,trace);
                jet={old.value,trace.derivative->transposeTimes({old.dx,old.dy})};}
            requireFirstJet(std::isfinite(jet.phi)&&LevelSetBoundaryTrace::Detail::finite(jet.gradient),"nonfinite pulled or prescribed first jet");
            jet=applySymmetryFirstJet(q,jet,config,treatment);
            rawJets.emplace(key,jet);return jet;};
        const auto center=[&](Vector2D q){const auto trace=firstQuery(q);
            requireFirstJet(!trace.trajectory.boundaryHit,"center prediction requires unresolved boundary data");
            return originalValue(q,trace);};
        if(!compensated)payload.emplace(updatedLocalPayload(support,cells,requestedSupport,pulled,center,config,scalar,report));
        else {
            std::set<int> stagedCells;
            for(const auto& q:requestedSupport.points){const auto trace=firstQuery({q.first,q.second});
                requireFirstJet(!trace.trajectory.boundaryHit&&trace.derivative.has_value(),"closed projection first query requires interior Jacobian");
                const auto p=trace.trajectory.point;bool found=false;
                for(int c:scalar.nativeCellCandidates(p,p)){const auto box=scalar.nativeCellBounds(c);
                    if(p.x>=box[0].x&&p.x<=box[1].x&&p.y>=box[0].y&&p.y<=box[1].y){
                        requireFirstJet(scalar.nativeCellActive(c),"closed projection missing original departure support");stagedCells.insert(c);found=true;}}
                requireFirstJet(found,"closed projection departure has no native support");}
            const std::vector<int> stageIds(stagedCells.begin(),stagedCells.end());report.stagedNativeCells=stageIds.size();
            const auto stagedSupport=support.sampleSupport(stageIds);
            auto stagedPayload=localFirstJetPayload(support,stageIds,stagedSupport,pulled,center,config);
            const auto projected=LevelSetHermiteGeometry::Snapshot::fromNativeLeafSubset(scalar,stagedPayload.patch);
            std::map<std::pair<double,double>,HermiteGalsBoundaryFirstJet> correctedJets;
            const auto corrected=[&](Vector2D q){const auto key=std::make_pair(q.x,q.y);const auto found=correctedJets.find(key);
                if(found!=correctedJets.end())return found->second;
                const auto first=firstQuery(q);
                requireFirstJet(!first.trajectory.boundaryHit&&first.derivative.has_value(),"closed projection first query requires interior Jacobian");
                ++serialStats.traceQueries;const auto second=motion.trace(first.trajectory.point);
                requireFirstJet(!second.trajectory.boundaryHit&&second.derivative.has_value(),"closed projection second query requires interior Jacobian");
                auto& slot=firstQueries.at(key).sourceJet;
                const auto one=sourceTransportJet(scalar,first.trajectory.point,sampling,serialStats,slot,true),
                    staged=localTransportJet(projected,first.trajectory.point,sampling,serialStats),
                    twice=sourceTransportJet(scalar,second.trajectory.point,sampling,serialStats,slot,true,false);
                const auto gtwice=second.derivative->transposeTimes({twice.dx,twice.dy});
                const Vector2D inner{one.dx+.5*(staged.dx-gtwice.x),one.dy+.5*(staged.dy-gtwice.y)};
                const auto gradient=first.derivative->transposeTimes(inner);const double value=one.value+.5*(staged.value-twice.value);
                requireFirstJet(std::isfinite(value)&&LevelSetBoundaryTrace::Detail::finite(gradient),"nonfinite compensated first jet");
                const auto result=applySymmetryFirstJet(q,{value,gradient},config,treatment);correctedJets.emplace(key,result);return result;};
            std::map<std::pair<double,double>,double> correctedValues;
            const auto correctedValue=[&](Vector2D q){const auto key=std::make_pair(q.x,q.y);
                if(requestedSupport.jetPoints.count(key))return corrected(q).phi;
                const auto found=correctedValues.find(key);if(found!=correctedValues.end())return found->second;
                const auto first=firstQuery(q);
                requireFirstJet(!first.trajectory.boundaryHit&&first.derivative.has_value(),"closed projection first query requires interior Jacobian");
                ++serialStats.traceQueries;const auto second=motion.trace(first.trajectory.point);
                requireFirstJet(!second.trajectory.boundaryHit&&second.derivative.has_value(),"closed projection second query requires interior Jacobian");
                auto& slot=firstQueries.at(key).sourceJet;
                const double one=sourceTransportValue(scalar,first.trajectory.point,sampling,serialStats,slot,true),
                    staged=localTransportValue(projected,first.trajectory.point,sampling,serialStats),
                    twice=sourceTransportValue(scalar,second.trajectory.point,sampling,serialStats,slot,true,false);
                const double value=one+.5*(staged-twice);
                requireFirstJet(std::isfinite(value),"nonfinite compensated center value");
                correctedValues.emplace(key,value);return value;
            };
            payload.emplace(localFirstJetPayload(support,cells,requestedSupport,corrected,correctedValue,config));
        }
        mergeLocalPointStats(report,serialStats);
        }
        report.timings.mixed=std::chrono::duration<double>(Clock::now()-kernelStart).count();
        report.transportedVertices=payload->primary.vertices.size();report.transportedEdges=payload->primary.edges.size();report.tracedCenters=cells.size();
        stage="local scalar binding";const auto bindStart=Clock::now();auto state=bind(*payload,targetTime);
        LAMP::Geometry::Detail::ExactExpansion::validateActualEnvironment();
        report.timings.publication=std::chrono::duration<double>(Clock::now()-bindStart).count();
        report.timings.total=std::chrono::duration<double>(Clock::now()-start).count();
        return LocalSampleProducts{std::move(state),report};
    } catch(const std::exception& e){throw std::runtime_error(std::string("[Interface/HermiteGalsFirstJetLocalSample source=")
        +std::to_string(id.value())+" arrival="+std::to_string(motion.velocityIdentity().value())+" stage="+stage+"] "+e.what());}
}
}

struct HermiteGalsFirstJetLocalPreparation::Data
{
    Data(HermiteGalsPublication actual,LocalUniformSupport topology):source(std::move(actual)),support(std::move(topology)) {}
    HermiteGalsPublication source;
    LocalUniformSupport support;
};
namespace
{
using CardinalPacket=LevelSetMappedSegment::QuarticCardinalData;
using TraceSegment=HermiteGalsFirstJetTraceSample::Segment;
using TraceSegmentKey=std::tuple<double,double,double,double,int>;
TraceSegmentKey traceSegmentKey(Vector2D begin,Vector2D end,int axis)
{return {begin.x,begin.y,end.x,end.y,axis};}
struct CardinalRow
{
    int cell=-1,normalAxis=0,position=0; // 0/1 native side, 2 exact midpoint
    std::array<Vector2D,2> bounds;
    std::array<int,2> vertices{{-1,-1}};
    std::array<std::array<int,2>,2> edges{};
    Vector2D begin,end;
};
double traceCoordinate(Vector2D point,int axis){return axis==0?point.x:point.y;}
bool exactMiddle(double point,double lower,double upper)
{
    namespace E=LAMP::Geometry::Detail::ExactExpansion;
    const auto twice=E::multiply(E::fromDouble(2),E::subtract(E::fromDouble(point),E::fromDouble(lower)));
    return E::sign(E::subtract(twice,E::subtract(E::fromDouble(upper),E::fromDouble(lower))))==0;
}
// Eligibility is a geometric proof, never a tolerance-based snapping rule.
std::optional<CardinalRow> cardinalRow(const LocalUniformSupport& support,
    const LevelSetHermiteGeometry::Snapshot& scalar,const std::vector<int>& cells,const TraceSegment& segment)
{
    const int normal=segment.normalAxis;
    requireFirstJet((normal==0||normal==1)&&LevelSetBoundaryTrace::Detail::finite(segment.begin)
        &&LevelSetBoundaryTrace::Detail::finite(segment.end),"invalid requested trace segment");
    const int tangent=1-normal;const auto lower=support.grid.minCoord(),upper=support.grid.maxCoord();
    for(const auto point:{segment.begin,segment.end})requireFirstJet(point.x>=lower.x&&point.x<=upper.x&&point.y>=lower.y&&point.y<=upper.y,
        "trace segment outside source rectangle");
    const double line=traceCoordinate(segment.begin,normal);
    requireFirstJet(line==traceCoordinate(segment.end,normal)
        &&traceCoordinate(segment.begin,tangent)!=traceCoordinate(segment.end,tangent),"trace must be a nonzero axis-aligned segment");
    const Vector2D lo{std::min(segment.begin.x,segment.end.x),std::min(segment.begin.y,segment.end.y)},
        hi{std::max(segment.begin.x,segment.end.x),std::max(segment.begin.y,segment.end.y)};
    const double first=traceCoordinate(lo,tangent),last=traceCoordinate(hi,tangent);
    std::optional<CardinalRow> result;
    for(int cell:scalar.nativeCellCandidates(lo,hi)){
        const auto box=scalar.nativeCellBounds(cell);
        const double low=traceCoordinate(box[0],normal),high=traceCoordinate(box[1],normal),
            lowT=traceCoordinate(box[0],tangent),highT=traceCoordinate(box[1],tangent);
        if(line<low||line>high||std::max(first,lowT)>=std::min(last,highT))continue;
        requireFirstJet(std::binary_search(cells.begin(),cells.end(),cell),"missing requested trace native support");
        int position=-1;
        if(line==low)position=0;else if(line==high)position=1;else if(exactMiddle(line,low,high))position=2;
        // One complete native row covers the original atom. Every other
        // positive-length owner must also be cardinal; endpoint-only cells
        // do not own a positive phase interval.
        if(position<0||first<lowT||last>highT)return {};
        const auto ij=support.cellPosition[cell];const int i=ij[0],j=ij[1];
        CardinalRow row;row.cell=cell;row.normalAxis=normal;row.position=position;row.bounds=box;
        row.begin=normal==0?Vector2D{line,box[0].y}:Vector2D{box[0].x,line};
        row.end=normal==0?Vector2D{line,box[1].y}:Vector2D{box[1].x,line};
        if(position<2){
            row.vertices=normal==0?std::array<int,2>{{support.idAt(i+position,j),support.idAt(i+position,j+1)}}:
                std::array<int,2>{{support.idAt(i,j+position),support.idAt(i+1,j+position)}};
            row.edges[0]=row.vertices;
            const auto a=support.point(row.vertices[0]),b=support.point(row.vertices[1]);
            if(a.x!=row.begin.x||a.y!=row.begin.y||b.x!=row.end.x||b.y!=row.end.y)return {};
        } else {
            row.edges=normal==0?std::array<std::array<int,2>,2>{{{{support.idAt(i,j),support.idAt(i+1,j)}},
                {{support.idAt(i,j+1),support.idAt(i+1,j+1)}}}}:
                std::array<std::array<int,2>,2>{{{{support.idAt(i,j),support.idAt(i,j+1)}},
                {{support.idAt(i+1,j),support.idAt(i+1,j+1)}}}};
            for(const auto& edge:row.edges){const auto a=support.point(edge[0]),b=support.point(edge[1]),q=a+(b-a)*.5;
                if(traceCoordinate(q,normal)!=line)return {};}
            if(traceCoordinate(support.cellCenters[cell],normal)!=line)return {};
        }
        if(!result||std::tie(box[0].y,box[0].x,cell)<std::tie(result->bounds[0].y,result->bounds[0].x,result->cell))result=row;
    }
    requireFirstJet(result.has_value(),"requested trace has no positive-length native owner");
    return result;
}
struct CardinalSupport
{
    LocalSampleSupport points;
    std::set<int> centers;
};
CardinalSupport cardinalSupport(const LocalUniformSupport& support,const std::vector<CardinalRow>& rows,
    const HermiteGalsPublicationConfig& config)
{
    CardinalSupport result;auto& selected=result.points;
    const auto add=[&](Vector2D q,bool jet){const LocalPointKey key{q.x,q.y};selected.points.insert(key);
        if(jet)selected.jetPoints.insert(key);};
    for(std::size_t k=0;k<config.volume.symmetryPlanes.size();++k){const auto& plane=config.volume.symmetryPlanes[k];
        requireFirstJet((plane.normalAxis==0||plane.normalAxis==1)&&std::isfinite(plane.coordinate),"invalid symmetry plane");
        requireFirstJet(plane.coordinate==traceCoordinate(support.grid.minCoord(),plane.normalAxis)
            ||plane.coordinate==traceCoordinate(support.grid.maxCoord(),plane.normalAxis),"symmetry plane must be an actual outer domain section");
        for(std::size_t j=0;j<k;++j)requireFirstJet(!(plane==config.volume.symmetryPlanes[j]),"duplicate symmetry plane");}
    for(const auto& row:rows){
        if(row.position<2){
            for(int vertex:row.vertices){selected.vertices.insert(vertex);add(support.point(vertex),true);}
            selected.edges.insert(row.edges[0]);
            const auto a=support.point(row.vertices[0]),b=support.point(row.vertices[1]);add(a+(b-a)*.5,false);
        } else {
            for(const auto& edge:row.edges){selected.edges.insert(edge);const auto a=support.point(edge[0]),b=support.point(edge[1]);
                add(a+(b-a)*.5,true);}
            result.centers.insert(row.cell);add(support.cellCenters[row.cell],false);
        }
    }
    // Preserve consumed primary symmetry checks, requesting a genuine
    // normal jet where an otherwise value-only midpoint needs it.
    for(const auto& edge:selected.edges){const auto a=support.point(edge[0]),b=support.point(edge[1]);
        for(const auto& plane:config.volume.symmetryPlanes)
            if(traceCoordinate(a,plane.normalAxis)==plane.coordinate&&traceCoordinate(b,plane.normalAxis)==plane.coordinate)
                add(a+(b-a)*.5,true);}
    return result;
}
template<class Job>
void tracePointBatch(std::size_t count,std::vector<std::exception_ptr>& failures,HermiteGalsFirstJetStepReport& report,
    HermiteGalsPointExecution execution,const Job& job)
{
    if(execution==HermiteGalsPointExecution::Parallel){localPointBatch(count,failures,report,job);return;}
    namespace E=LAMP::Geometry::Detail::ExactExpansion;
    E::validateEnvironment();LocalPointStats stats;const auto start=std::chrono::steady_clock::now();
    for(std::size_t k=0;k<count;++k){
        try{E::validateEnvironment();job(k,stats);E::validateEnvironment();}
        catch(...){failures[k]=std::current_exception();}
    }
    if(count)stats.seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    E::validateEnvironment();report.pointJobs+=count;mergeLocalPointStats(report,stats);
}
// Shared actual two-map transport for full local Q4 and cardinal traces.
// Their requested point inventories differ; the staged source stays full Q4.
LocalPointJets transportLocalPointJets(const LocalUniformSupport& support,const LocalSampleSupport& requestedSupport,
    const LevelSetHermiteGeometry::Snapshot& scalar,const HermiteGalsResolvedMotion& motion,
    const HermiteGalsPublicationConfig& config,bool compensated,HermiteGalsPointExecution execution,
    const HermiteGalsPointSamplingConfig& sampling,HermiteGalsFirstJetStepReport& report,HermiteGalsBoundaryTreatment treatment)
{
    LocalPointTraces first(requestedSupport.points);
    tracePointBatch(first.keys.size(),first.failures,report,execution,[&](std::size_t k,LocalPointStats& stats){
        ++stats.traceQueries;const auto key=first.keys[k];first.values[k].emplace(motion.trace({key.first,key.second}));});
    const auto rawPoints=[&](const LocalSampleSupport& preparedSupport){
        LocalPointJets jets(preparedSupport.points);const auto& jetKeys=preparedSupport.jetPoints;
        tracePointBatch(jets.keys.size(),jets.failures,report,execution,[&](std::size_t k,LocalPointStats& stats){
            const auto key=jets.keys[k];const Vector2D q{key.first,key.second};const auto saved=first.find(key);
            std::optional<LevelSetBoundaryTrace::JetResult> fresh;
            if(saved==first.keys.size()){++stats.traceQueries;fresh.emplace(motion.trace(q));}
            const auto& trace=fresh?*fresh:first.at(saved);std::optional<SourceJetSlot> temporary;
            auto& slot=compensated&&saved<first.keys.size()?first.sourceJets[saved]:temporary;
            if(jetKeys.count(key)==0){
                requireFirstJet(!trace.trajectory.boundaryHit,"center prediction requires unresolved boundary data");
                jets.centerValues[k]=sourceTransportValue(scalar,trace.trajectory.point,sampling,stats,slot,compensated);return;}
            const bool outer=q.x==support.grid.minCoord().x||q.x==support.grid.maxCoord().x||q.y==support.grid.minCoord().y||q.y==support.grid.maxCoord().y;
            if(outer&&trace.trajectory.boundaryHit)requireFirstJet(false,"outer boundary hit requires prescribed first jet");
            requireInteriorFirstJetTrace(trace,q);
            const auto old=sourceTransportJet(scalar,trace.trajectory.point,sampling,stats,slot,compensated);
            const auto jet=applySymmetryFirstJet(q,{old.value,trace.derivative->transposeTimes({old.dx,old.dy})},config,treatment);
            requireFirstJet(std::isfinite(jet.phi)&&LevelSetBoundaryTrace::Detail::finite(jet.gradient),"nonfinite pulled or prescribed first jet");
            jets.values[k].emplace(jet);
        });return jets;
    };
    if(!compensated)return rawPoints(requestedSupport);
    std::set<int> stagedCells;
    for(std::size_t k=0;k<first.keys.size();++k){const auto& trace=first.at(k);
        requireFirstJet(!trace.trajectory.boundaryHit&&trace.derivative.has_value(),"closed projection first query requires interior Jacobian");
        const auto p=trace.trajectory.point;bool found=false;
        for(int c:scalar.nativeCellCandidates(p,p)){const auto box=scalar.nativeCellBounds(c);
            if(p.x>=box[0].x&&p.x<=box[1].x&&p.y>=box[0].y&&p.y<=box[1].y){
                requireFirstJet(scalar.nativeCellActive(c),"closed projection missing original departure support");stagedCells.insert(c);found=true;}}
        requireFirstJet(found,"closed projection departure has no native support");}
    const std::vector<int> stageIds(stagedCells.begin(),stagedCells.end());report.stagedNativeCells=stageIds.size();
    const auto stagedSupport=support.sampleSupport(stageIds);const auto raw=rawPoints(stagedSupport);
    const auto stagedPayload=localFirstJetPayload(support,stageIds,stagedSupport,
        [&](Vector2D q){return raw.at(q);},[&](Vector2D q){return raw.valueAt(q);},config);
    const auto projected=LevelSetHermiteGeometry::Snapshot::fromNativeLeafSubset(scalar,stagedPayload.patch);
    LocalPointJets corrected(requestedSupport.points);
    tracePointBatch(corrected.keys.size(),corrected.failures,report,execution,[&](std::size_t k,LocalPointStats& stats){
        const auto key=corrected.keys[k];const auto index=first.find(key);const auto& trace=first.at(index);
        auto& slot=first.sourceJets[index];
        requireFirstJet(!trace.trajectory.boundaryHit&&trace.derivative.has_value(),"closed projection first query requires interior Jacobian");
        ++stats.traceQueries;const auto second=motion.trace(trace.trajectory.point);
        requireFirstJet(!second.trajectory.boundaryHit&&second.derivative.has_value(),"closed projection second query requires interior Jacobian");
        if(requestedSupport.jetPoints.count(key)==0){
            const double one=sourceTransportValue(scalar,trace.trajectory.point,sampling,stats,slot,true),
                staged=localTransportValue(projected,trace.trajectory.point,sampling,stats),
                twice=sourceTransportValue(scalar,second.trajectory.point,sampling,stats,slot,true,false);
            const double value=one+.5*(staged-twice);
            requireFirstJet(std::isfinite(value),"nonfinite compensated center value");corrected.centerValues[k]=value;return;}
        const auto one=sourceTransportJet(scalar,trace.trajectory.point,sampling,stats,slot,true);
        const auto staged=localTransportJet(projected,trace.trajectory.point,sampling,stats);
        const auto twice=sourceTransportJet(scalar,second.trajectory.point,sampling,stats,slot,true,false);
        const auto gtwice=second.derivative->transposeTimes({twice.dx,twice.dy});
        const Vector2D inner{one.dx+.5*(staged.dx-gtwice.x),one.dy+.5*(staged.dy-gtwice.y)};
        const auto gradient=trace.derivative->transposeTimes(inner);const double value=one.value+.5*(staged.value-twice.value);
        requireFirstJet(std::isfinite(value)&&LevelSetBoundaryTrace::Detail::finite(gradient),"nonfinite compensated first jet");
        corrected.values[k].emplace(applySymmetryFirstJet({key.first,key.second},{value,gradient},config,treatment));
    });return corrected;
}
void mergeTraceAttempt(HermiteGalsFirstJetStepReport& full,const HermiteGalsFirstJetStepReport& attempt)
{
    full.transportedVertices+=attempt.transportedVertices;full.transportedEdges+=attempt.transportedEdges;full.tracedCenters+=attempt.tracedCenters;
    full.stagedNativeCells+=attempt.stagedNativeCells;full.evaluatedNativeCells+=attempt.evaluatedNativeCells;
    full.evaluatedPointQueries+=attempt.evaluatedPointQueries;full.pointJobs+=attempt.pointJobs;full.pointTraceQueries+=attempt.pointTraceQueries;
    full.pointThreads=std::max(full.pointThreads,attempt.pointThreads);full.pointWorkerSeconds+=attempt.pointWorkerSeconds;
    full.numericalPointQueries+=attempt.numericalPointQueries;full.numericalValueFallbacks+=attempt.numericalValueFallbacks;
    full.numericalGradientFallbacks+=attempt.numericalGradientFallbacks;
    full.maxAcceptedValueErrorBound=std::max(full.maxAcceptedValueErrorBound,attempt.maxAcceptedValueErrorBound);
    full.maxAcceptedGradientErrorBound=std::max(full.maxAcceptedGradientErrorBound,attempt.maxAcceptedGradientErrorBound);
    full.timings.input+=attempt.timings.input;full.timings.boundary+=attempt.timings.boundary;
    full.timings.vertices+=attempt.timings.vertices;full.timings.edges+=attempt.timings.edges;full.timings.centers+=attempt.timings.centers;
    full.timings.mixed+=attempt.timings.mixed;full.timings.publication+=attempt.timings.publication;
}
} // namespace

HermiteGalsFirstJetLocalPreparation::HermiteGalsFirstJetLocalPreparation(std::shared_ptr<const Data> data):data_(std::move(data)) {}
const HermiteGalsFirstJetLocalPreparation::Data& HermiteGalsFirstJetLocalPreparation::checked() const
{requireFirstJet(bool(data_),"[Interface/HermiteGalsFirstJetLocalPreparation] moved-from preparation");return *data_;}
HermiteGalsFirstJetLocalPreparation HermiteGalsFirstJetLocalPreparation::prepare(const QuadtreeMesh& mesh,
    const HermiteGalsPublication& source)
{
    namespace Exact=LAMP::Geometry::Detail::ExactExpansion;
    Exact::validateActualEnvironment();source.validateMesh(mesh);
    requireFirstJet(supportsHermiteGalsFirstJetLocalSample(mesh,source),"local sample requires complete uniform XY native Q4 source");
    LocalUniformSupport support(mesh);
    auto data=std::make_shared<const Data>(source,std::move(support));
    Exact::validateActualEnvironment();return HermiteGalsFirstJetLocalPreparation(std::move(data));
}
HermiteGalsIdentity HermiteGalsFirstJetLocalPreparation::sourceIdentity() const{return checked().source.state().identity();}
void HermiteGalsFirstJetLocalPreparation::validateMesh(const QuadtreeMesh& mesh) const{checked().source.validateMesh(mesh);}
HermiteGalsState HermiteGalsFirstJetLocalPreparation::bindLocal(HermiteGalsTimeLayer time,
    const HermiteGalsPrimaryState& primary,const HermiteGalsClosure& closure,
    const LevelSetHermiteGeometry::Snapshot::NativeLeafSubset& patch) const
{
    const auto& prepared=checked();
    requireFirstJet(primary.geometryCells==patch.nativeCells&&primary.vertices.size()==patch.nativeVertices.size()
        &&patch.nodes.size()==primary.vertices.size()&&closure.vertexMixedDerivatives.size()==primary.vertices.size()
        &&patch.cellCenters==closure.cellBubbleCoefficients&&primary.edges.size()==patch.edgeMidpoints.size(),
        "prepared local primary/closure payload mismatch");
    for(std::size_t k=0;k<primary.vertices.size();++k){
        const auto& v=primary.vertices[k];const auto& q=patch.nodes[k];
        requireFirstJet(v.vertex==patch.nativeVertices[k]&&v.phi==q.value&&v.dx==q.dx&&v.dy==q.dy
            &&closure.vertexMixedDerivatives[k]==q.dxy,"prepared local vertex payload mismatch");
    }
    for(std::size_t k=0;k<primary.edges.size();++k){
        const auto& e=primary.edges[k];const auto& q=patch.edgeMidpoints[k];
        requireFirstJet(e.vertices==q.vertices&&e.phi==q.value&&e.dn==q.normalDerivative,"prepared local edge payload mismatch");
    }
    // This is the sole prepared binding path: build from the same private
    // compact payload, never accept a caller-provided scalar or activity mask.
    auto scalar=LevelSetHermiteGeometry::Snapshot::fromNativeLeafSubset(prepared.source.state().geometry().snapshot(),patch);
    return HermiteGalsCandidate::bindPreparedLocal(prepared.source.state(),time,primary,closure,std::move(scalar),
        patch.nativeCells,prepared.support.nativePoints,prepared.support.independent);
}
HermiteGalsFirstJetLocalSample HermiteGalsFirstJetLocalPreparation::sample(const HermiteGalsResolvedMotion& motion,
    const std::vector<int>& cells,const HermiteGalsPublicationConfig& config,HermiteGalsTimeProjection projection,
    HermiteGalsPointExecution execution,const HermiteGalsPointSamplingConfig& sampling,HermiteGalsBoundaryTreatment treatment) const
{
    const auto start=std::chrono::steady_clock::now();const auto& prepared=checked();
    auto result=localFirstJetKernel(prepared.source,prepared.support,motion,cells,{},config,
        treatment,projection,execution,sampling,
        [&](const LocalFirstJetPayload& payload,HermiteGalsTimeLayer time){return bindLocal(time,payload.primary,payload.closure,payload.patch);},start);
    return HermiteGalsFirstJetLocalSample(std::make_shared<const HermiteGalsFirstJetLocalSample::Data>(
        sourceIdentity(),motion.velocityIdentity(),std::move(result.state),result.report));
}
struct HermiteGalsFirstJetTraceSample::Data
{
    Data(HermiteGalsPublication actual,HermiteGalsIdentity velocity,HermiteGalsTimeLayer target,
        std::vector<int> support,std::vector<Segment> inputs)
        :source(std::move(actual)),arrival(velocity),clock(target),cells(std::move(support)),segments(std::move(inputs)){}
    HermiteGalsPublication source;
    HermiteGalsIdentity arrival;
    HermiteGalsTimeLayer clock;
    std::vector<int> cells;
    std::vector<Segment> segments;
    std::vector<std::shared_ptr<const CardinalPacket>> packets;
    // Each actual row shares its own lazy exact coefficients; every atom
    // still creates an independent restriction/root/refinement query.
    std::vector<LevelSetMappedSegment::PreparedCardinalQuarticRow> preparedRows;
    std::map<TraceSegmentKey,std::shared_ptr<const CardinalPacket>> segmentOwners;
    std::optional<HermiteGalsFirstJetLocalSample> fallback;
    std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> fallbackScalar;
    HermiteGalsFirstJetStepReport report;
};
HermiteGalsFirstJetTraceSample::HermiteGalsFirstJetTraceSample(std::shared_ptr<const Data> data):data_(std::move(data)){}
const HermiteGalsFirstJetTraceSample::Data& HermiteGalsFirstJetTraceSample::checked() const
{requireFirstJet(bool(data_),"[Interface/HermiteGalsFirstJetTraceSample] moved-from result");return *data_;}
HermiteGalsIdentity HermiteGalsFirstJetTraceSample::sourceIdentity() const{return checked().source.state().identity();}
HermiteGalsIdentity HermiteGalsFirstJetTraceSample::arrivalVelocityIdentity() const{return checked().arrival;}
HermiteGalsTimeLayer HermiteGalsFirstJetTraceSample::sourceTimeLayer() const{return checked().source.state().timeLayer();}
HermiteGalsTimeLayer HermiteGalsFirstJetTraceSample::timeLayer() const{return checked().clock;}
HermiteGalsMeasure HermiteGalsFirstJetTraceSample::measure() const{return checked().source.state().primary().measure;}
const std::vector<int>& HermiteGalsFirstJetTraceSample::nativeCells() const{return checked().cells;}
std::size_t HermiteGalsFirstJetTraceSample::nativeCellCount() const{return checked().source.state().geometry().snapshot().nativeCellCount();}
std::array<Vector2D,2> HermiteGalsFirstJetTraceSample::nativeCellBounds(int cell) const
{
    const auto& data=checked();requireFirstJet(std::binary_search(data.cells.begin(),data.cells.end(),cell),"missing trace native support");
    return data.source.state().geometry().snapshot().nativeCellBounds(cell);
}
const HermiteGalsFirstJetStepReport& HermiteGalsFirstJetTraceSample::report() const{return checked().report;}
HermiteGalsFirstJetTraceSample::Trace HermiteGalsFirstJetTraceSample::trace(std::size_t index) const
{
    namespace E=LAMP::Geometry::Detail::ExactExpansion;E::validateActualEnvironment();const auto& data=checked();
    requireFirstJet(index<data.segments.size(),"trace segment index out of range");const auto& segment=data.segments[index];
    if(data.fallbackScalar){
        auto value=LevelSetMappedSegment::nativeQuarticTrace(data.fallbackScalar,segment.begin,segment.end,segment.normalAxis);
        E::validateActualEnvironment();return Trace(std::move(value));
    }
    auto value=LevelSetMappedSegment::cardinalQuarticTrace(data.preparedRows.at(index),segment.begin,segment.end);
    requireFirstJet(value.has_value(),"immutable cardinal trace unexpectedly requires complete Taylor data");
    E::validateActualEnvironment();return Trace(std::move(*value));
}
bool HermiteGalsFirstJetTraceSample::ownsTrace(const Trace& trace) const
{
    const auto& data=checked();
    return std::visit([&](const auto& value){
        using Value=std::decay_t<decltype(value)>;
        if constexpr(std::is_same_v<Value,LevelSetMappedSegment::NativeQuarticTrace>)
            return data.fallbackScalar&&value.usesSource(*data.fallbackScalar);
        else {
            const auto found=data.segmentOwners.find(traceSegmentKey(value.begin(),value.end(),value.normalAxis()));
            return found!=data.segmentOwners.end()&&value.usesSource(*found->second);
        }
    },trace);
}
HermiteGalsFirstJetTraceSample HermiteGalsFirstJetLocalPreparation::sampleTraces(const HermiteGalsResolvedMotion& motion,
    const std::vector<int>& requestedNativeCells,const std::vector<HermiteGalsFirstJetTraceSample::Segment>& segments,
    const HermiteGalsPublicationConfig& config,HermiteGalsTimeProjection projection,
    HermiteGalsPointExecution execution,const HermiteGalsPointSamplingConfig& sampling,HermiteGalsBoundaryTreatment treatment) const
{
    using Clock=std::chrono::steady_clock;namespace E=LAMP::Geometry::Detail::ExactExpansion;
    const auto start=Clock::now();const auto& prepared=checked();const char* stage="input";
    try {
        E::validateActualEnvironment();
        requireFirstJet(noQueryBoundaryTreatment(treatment),"trace samples require no-query boundary treatment");
        validateSymmetryFirstJet(prepared.support.grid.minCoord(),prepared.support.grid.maxCoord(),prepared.source.state().primary(),
            [&](int vertex){return prepared.support.point(vertex);},{},config,treatment);
        requireFirstJet(!requestedNativeCells.empty()&&!segments.empty(),"empty requested trace support");
        auto cells=requestedNativeCells;std::sort(cells.begin(),cells.end());
        for(std::size_t k=0;k<cells.size();++k){requireFirstJet(cells[k]>=0&&static_cast<std::size_t>(cells[k])<prepared.support.cellPosition.size(),
            "invalid requested native cell");requireFirstJet(k==0||cells[k]!=cells[k-1],"duplicate requested native cell");}
        const auto sourceTime=prepared.source.state().timeLayer(),targetTime=motion.targetTime(),bound=motion.sourceTime();
        using Position=HermiteGalsTimeLayer::Position;
        requireFirstJet(sourceTime.position==Position::Initial||sourceTime.position==Position::End,"source must be accepted Initial/End");
        requireFirstJet(targetTime.position==Position::Midpoint||targetTime.position==Position::End,"target must be Midpoint/End");
        requireFirstJet(sourceTime.step<std::numeric_limits<std::int64_t>::max()&&targetTime.step==sourceTime.step+1,"target must be the next step");
        requireFirstJet(motion.sourceIdentity()==sourceIdentity()&&bound.step==sourceTime.step&&bound.begin==sourceTime.begin
            &&bound.end==sourceTime.end&&bound.position==sourceTime.position,"resolved motion source identity/clock differs");
        requireFirstJet(projection==HermiteGalsTimeProjection::CharacteristicPullback||projection==HermiteGalsTimeProjection::ClosedProjectionCompensated,
            "unknown temporal projection rule");
        requireFirstJet(execution==HermiteGalsPointExecution::Serial||execution==HermiteGalsPointExecution::Parallel,"unknown point execution rule");
        requireFirstJet(sampling.mode==HermiteGalsPointSampling::Exact||sampling.mode==HermiteGalsPointSampling::Numerical,"unknown transport point sampling rule");
        if(sampling.mode==HermiteGalsPointSampling::Numerical)
            requireFirstJet(std::isfinite(sampling.valueAbsoluteTolerance)&&sampling.valueAbsoluteTolerance>0
                &&std::isfinite(sampling.gradientAbsoluteTolerance)&&sampling.gradientAbsoluteTolerance>0,"invalid numerical transport point budgets");
#ifdef _OPENMP
        if(execution==HermiteGalsPointExecution::Parallel)requireFirstJet(!omp_in_parallel(),"point-parallel sample cannot run within an outer OpenMP team");
#endif
        const auto tolerance=[](double value){return std::isfinite(value)&&value>0&&value<1;};
        requireFirstJet(tolerance(config.relativeInterfaceResidualTolerance)&&tolerance(config.unitNormalTolerance),"invalid publication tolerance");
        auto data=std::make_shared<HermiteGalsFirstJetTraceSample::Data>(prepared.source,motion.velocityIdentity(),targetTime,cells,segments);
        auto& report=data->report;const auto& support=prepared.support;const auto& scalar=prepared.source.state().geometry().snapshot();
        report.sourceCells=prepared.source.state().primary().geometryCells.size();report.sourceNativeFaces=prepared.source.domain().nativeFaceCount();
        report.sourceCutCells=prepared.source.volumes().report().cutCells;report.sourceGraphs=prepared.source.report().graphs;
        report.requestedNativeCells=cells.size();report.pointExecution=execution;report.pointSampling=sampling.mode;
        const auto fullFallback=[&](){
            auto actual=sample(motion,cells,config,projection,execution,sampling,treatment);
            const auto view=std::make_shared<const HermiteGalsGeometryView>(actual.geometry());
            data->fallbackScalar=std::shared_ptr<const LevelSetHermiteGeometry::Snapshot>(view,&view->snapshot());
            auto full=actual.report();mergeTraceAttempt(full,data->report);data->report=full;
            data->fallback.emplace(std::move(actual));data->packets.clear();data->preparedRows.clear();data->segmentOwners.clear();
            data->report.timings.total=std::chrono::duration<double>(Clock::now()-start).count();
            E::validateActualEnvironment();return HermiteGalsFirstJetTraceSample(std::move(data));
        };
        stage="exact cardinal support";std::vector<CardinalRow> rows;rows.reserve(segments.size());bool eligible=true;
        for(const auto& segment:segments){const auto row=cardinalRow(support,scalar,cells,segment);
            if(row)rows.push_back(*row);else eligible=false;}
        report.timings.input=std::chrono::duration<double>(Clock::now()-start).count();
        if(!eligible)return fullFallback();
        auto selected=cardinalSupport(support,rows,config);
        validateGradientUpdate(config,scalar,projection==HermiteGalsTimeProjection::ClosedProjectionCompensated,false,treatment);
        if(config.gradientUpdate==HermiteGalsGradientUpdate::CompatibleIncrement)
            support.projectionGrid().expand(selected.points.jetPoints,selected.points.points,config);
        stage="actual cardinal point transport";const auto kernelStart=Clock::now();
        const auto values=transportLocalPointJets(support,selected.points,scalar,motion,config,
            projection==HermiteGalsTimeProjection::ClosedProjectionCompensated,execution,sampling,report,treatment);
        std::unique_ptr<IncrementGradientProjection> projected;
        if(config.gradientUpdate==HermiteGalsGradientUpdate::CompatibleIncrement)
            projected=std::make_unique<IncrementGradientProjection>(support.projectionGrid(),scalar,config,[&](Vector2D q){return values.valueAt(q);});
        const auto jetAt=[&](Vector2D q){return projected?projected->at(q):values.at(q);};
        for(int vertex:selected.points.vertices){const auto point=support.point(vertex);const auto jet=jetAt(point);
            for(const auto& plane:config.volume.symmetryPlanes)if(traceCoordinate(point,plane.normalAxis)==plane.coordinate)
                requireFirstJet((plane.normalAxis==0?jet.gradient.x:jet.gradient.y)==0,"symmetry requires exactly zero primary vertex normal derivative");}
        for(const auto& edge:selected.points.edges){const auto a=support.point(edge[0]),b=support.point(edge[1]),point=a+(b-a)*.5;
            for(const auto& plane:config.volume.symmetryPlanes)
                if(traceCoordinate(a,plane.normalAxis)==plane.coordinate&&traceCoordinate(b,plane.normalAxis)==plane.coordinate){
                    const auto jet=jetAt(point);const double dn=a.y==b.y?jet.gradient.y:jet.gradient.x;
                    requireFirstJet(dn==0,"symmetry requires exactly zero primary edge normal derivative");}}
        report.transportedVertices=selected.points.vertices.size();report.transportedEdges=selected.points.edges.size();report.tracedCenters=selected.centers.size();
        report.timings.mixed=std::chrono::duration<double>(Clock::now()-kernelStart).count();
        // evaluatedNativeCells stays zero: no final Q4 was prepared. Staged
        // donor cells and all actually executed point work remain reported.
        stage="immutable cardinal packets";const auto packetStart=Clock::now();bool zero=false;
        struct RowOwner {
            std::shared_ptr<const CardinalPacket> packet;
            LevelSetMappedSegment::PreparedCardinalQuarticRow prepared;
        };
        std::map<std::tuple<int,int,int>,RowOwner> rowOwners;
        data->packets.reserve(rows.size());data->preparedRows.reserve(rows.size());
        for(std::size_t k=0;k<rows.size();++k){const auto& row=rows[k];const auto key=std::make_tuple(row.cell,row.normalAxis,row.position);
            auto found=rowOwners.find(key);std::shared_ptr<const CardinalPacket> packet;
            if(found!=rowOwners.end())packet=found->second.packet;
            else {
                CardinalPacket datum;datum.begin=row.begin;datum.end=row.end;datum.normalAxis=row.normalAxis;datum.nativeCellIndex=row.cell;
                const int tangent=1-row.normalAxis;
                if(row.position<2){
                    for(int end=0;end<2;++end){const auto jet=jetAt(support.point(row.vertices[end]));
                        datum.values[end]=jet.phi;datum.tangentDerivatives[end]=tangent==0?jet.gradient.x:jet.gradient.y;}
                    const auto a=support.point(row.vertices[0]),b=support.point(row.vertices[1]);datum.midpointValue=values.valueAt(a+(b-a)*.5);
                } else {
                    for(int end=0;end<2;++end){const auto& edge=row.edges[end];const auto a=support.point(edge[0]),b=support.point(edge[1]);
                        const auto jet=jetAt(a+(b-a)*.5);datum.values[end]=jet.phi;
                        datum.tangentDerivatives[end]=tangent==0?jet.gradient.x:jet.gradient.y;}
                    datum.midpointValue=values.valueAt(support.cellCenters[row.cell]);
                }
                requireFirstJet(std::isfinite(datum.values[0])&&std::isfinite(datum.values[1])&&std::isfinite(datum.tangentDerivatives[0])
                    &&std::isfinite(datum.tangentDerivatives[1])&&std::isfinite(datum.midpointValue),"nonfinite cardinal trace data");
                // Cardinal independence proves exact zero from all five
                // actual data; root isolation occurs only in trace(q).
                zero|=datum.values[0]==0&&datum.values[1]==0&&datum.tangentDerivatives[0]==0&&datum.tangentDerivatives[1]==0&&datum.midpointValue==0;
                packet=std::make_shared<const CardinalPacket>(datum);
                // Preparation retains this authenticated actual packet only;
                // exact row coefficients remain lazy in the shared Core owner.
                found=rowOwners.emplace(key,RowOwner{packet,
                    LevelSetMappedSegment::prepareCardinalQuarticRow(packet)}).first;
            }
            data->packets.push_back(packet);data->preparedRows.push_back(found->second.prepared);
            const auto& segment=segments[k];
            data->segmentOwners.emplace(traceSegmentKey(segment.begin,segment.end,segment.normalAxis),packet);
        }
        report.timings.publication=std::chrono::duration<double>(Clock::now()-packetStart).count();
        if(projected)report.evaluatedPointQueries+=projected->pointQueries();
        if(zero)return fullFallback();
        E::validateActualEnvironment();report.timings.total=std::chrono::duration<double>(Clock::now()-start).count();
        return HermiteGalsFirstJetTraceSample(std::move(data));
    } catch(const std::exception& e){throw std::runtime_error(std::string("[Interface/HermiteGalsFirstJetTraceSample source=")
        +std::to_string(sourceIdentity().value())+" arrival="+std::to_string(motion.velocityIdentity().value())+" stage="+stage+"] "+e.what());}
}

HermiteGalsFirstJetLocalSample prepareHermiteGalsFirstJetLocalSample(const QuadtreeMesh& mesh,
    const HermiteGalsPublication& source,const HermiteGalsResolvedMotion& motion,const std::vector<int>& requestedNativeCells,
    const HermiteGalsBoundaryFirstJetQuery& boundary,const HermiteGalsPublicationConfig& config,
    HermiteGalsBoundaryTreatment treatment,HermiteGalsTimeProjection projection,HermiteGalsPointExecution execution,
    const HermiteGalsPointSamplingConfig& sampling)
{
    const auto start=std::chrono::steady_clock::now();
    source.validateMesh(mesh);
    requireFirstJet(supportsHermiteGalsFirstJetLocalSample(mesh,source),"local sample requires complete uniform XY native Q4 source");
    LocalUniformSupport support(mesh);
    auto result=localFirstJetKernel(source,support,motion,requestedNativeCells,boundary,config,treatment,projection,execution,sampling,
        [&](const LocalFirstJetPayload& payload,HermiteGalsTimeLayer time){
            HermiteGalsCandidate candidate(mesh,time,payload.primary);
            auto scalar=LevelSetHermiteGeometry::Snapshot::fromNativeLeafSubset(source.state().geometry().snapshot(),payload.patch);
            candidate.bindGeometry(mesh,candidate.identity(),payload.closure,std::move(scalar));return std::move(candidate).finish();
        },start);
    return HermiteGalsFirstJetLocalSample(std::make_shared<const HermiteGalsFirstJetLocalSample::Data>(
        source.state().identity(),motion.velocityIdentity(),std::move(result.state),result.report));
}

}
