#include "Modules/Interface/HermiteGalsFirstJetStep.h"

#include "Modules/Interface/detail/HermiteGalsFirstJet.h"
#include <algorithm>
#include <cmath>
#include <chrono>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace Quantum::Interface
{
using Detail::requireFirstJet;
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

HermiteGalsFirstJetStep prepareHermiteGalsFirstJetStep(const QuadtreeMesh& mesh,
    const HermiteGalsPublication& source,const LevelSetBoundaryTrace::Domain& domain,
    const HermiteGalsVelocityView& previous,const HermiteGalsVelocityView& arrival,
    const HermiteGalsBoundaryFirstJetQuery& boundary,const HermiteGalsPublicationConfig& config)
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
    return prepareHermiteGalsFirstJetStep(mesh,source,motion,boundary,config);
}
HermiteGalsFirstJetStep prepareHermiteGalsFirstJetStep(const QuadtreeMesh& mesh,
    const HermiteGalsPublication& source,const HermiteGalsResolvedMotion& motion,
    const HermiteGalsBoundaryFirstJetQuery& boundary,const HermiteGalsPublicationConfig& config,HermiteGalsBoundaryTreatment treatment)
{
    using Clock=std::chrono::steady_clock;const auto start=Clock::now();auto stamp=start;
    const auto lap=[&](){const auto now=Clock::now();const auto seconds=std::chrono::duration<double>(now-stamp).count();stamp=now;return seconds;};
    HermiteGalsIdentity id;const char* stage="input";
    try {
        id=source.state().identity();source.validateMesh(mesh);
        const auto& old=source.state();const auto& scalar=old.geometry().snapshot();
        const auto sourceTime=old.timeLayer(),targetTime=motion.targetTime();
        requireFirstJet(treatment==HermiteGalsBoundaryTreatment::InflowHistory||
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
            && scalar.reconstruction()==LevelSetHermiteGeometry::Reconstruction::EnrichedQuartic,"requires native Q4 source");
        const auto bound=motion.sourceTime();
        requireFirstJet(motion.sourceIdentity()==id&&bound.step==sourceTime.step&&bound.begin==sourceTime.begin
            &&bound.end==sourceTime.end&&bound.position==sourceTime.position,"resolved motion source identity/clock differs");
        auto primary=old.primary();HermiteGalsFirstJetStepReport report;
        report.sourceCells=primary.geometryCells.size();report.sourceNativeFaces=source.domain().nativeFaceCount();
        report.sourceCutCells=source.volumes().report().cutCells;report.sourceGraphs=source.report().graphs;
        report.timings.input=lap();
        const auto point=[&](int v){return mesh.verts(mesh.vertIndexToLabel[v]).coord();};
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
            const auto trace=motion.trace(q);
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
        const auto pulled=[&](Vector2D q){
            const auto trace=motion.trace(q);
            requireFirstJet(!trace.trajectory.boundaryHit,"interior Jacobian unavailable: boundary event");
            const auto jet=scalar.scalarJetNative(trace.trajectory.point);
            const auto gradient=trace.derivative->transposeTimes({jet.dx,jet.dy});
            requireFirstJet(std::isfinite(jet.value)&&LevelSetBoundaryTrace::Detail::finite(gradient),"nonfinite pulled first jet");
            return HermiteGalsBoundaryFirstJet{jet.value,gradient};
        };
        stage="vertex first-jet transport";
        for(int v:vertices){const auto jet=pulled(point(v));
            auto it=std::lower_bound(primary.vertices.begin(),primary.vertices.end(),v,[](const auto& a,int b){return a.vertex<b;});
            it->phi=jet.phi;it->dx=jet.gradient.x;it->dy=jet.gradient.y;
        }
        report.timings.vertices=lap();
        stage="edge first-jet transport";
        for(const auto& e:edges){const auto a=point(e[0]),b=point(e[1]);const auto jet=pulled(a+(b-a)*.5);
            auto it=std::lower_bound(primary.edges.begin(),primary.edges.end(),e,[](const auto& a,const auto& b){return a.vertices<b;});
            it->phi=jet.phi;it->dn=a.y==b.y?jet.gradient.y:jet.gradient.x;
        }
        report.timings.edges=lap();
        stage="center prediction";
        std::vector<double> centers;centers.reserve(primary.geometryCells.size());
        for(int cell:primary.geometryCells){
            const auto trace=motion.trace(mesh.cells(mesh.cellIndexToLabel[cell]).centerCoord());
            requireFirstJet(!trace.trajectory.boundaryHit,"center prediction requires unresolved boundary data");
            centers.push_back(scalar.scalarJetNative(trace.trajectory.point).value);
        }
        report.transportedVertices=vertices.size();report.transportedEdges=edges.size();report.tracedCenters=centers.size();
        HermiteGalsCandidate candidate(mesh,targetTime,primary);
        report.timings.centers=lap();
        stage="mixed closure";
        const auto mixed=closeHermiteGalsMixedAmr(mesh,candidate,config.volume.symmetryPlanes);
        const HermiteGalsCenterPrediction center{candidate.identity(),centers};
        report.timings.mixed=lap();
        stage="unique publication";
        auto publication=publishHermiteGalsFirstJetGeometry(mesh,std::move(candidate),mixed,center,config);
        report.timings.publication=lap();
        report.timings.total=std::chrono::duration<double>(Clock::now()-start).count();
        return HermiteGalsFirstJetStep(std::make_shared<const HermiteGalsFirstJetStep::Data>(id,motion.velocityIdentity(),std::move(publication),report));
    } catch(const std::exception& e) {
        throw std::runtime_error("[Interface/HermiteGalsFirstJetStep source="+std::to_string(id.value())
            +" arrival="+std::to_string(motion.velocityIdentity().value())+" stage="+stage+"] "+e.what());
    }
}
}
