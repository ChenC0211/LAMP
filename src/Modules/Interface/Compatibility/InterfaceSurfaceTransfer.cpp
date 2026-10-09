#include "Modules/Interface/Compatibility/InterfaceSurfaceTransfer.h"
#include "Modules/Interface/Compatibility/InterfaceConservativeGeometry.h"
#include "Modules/Interface/Compatibility/InterfaceTopology.h"
#include "Core/Mesh/CartesianMeasure.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <sstream>
#include <iomanip>

namespace
{
void require(bool ok,const char* why)
{if (!ok) throw std::invalid_argument(std::string("[InterfaceSurfaceTransfer] ")+why);}
double length(Vector2D v) {return std::hypot(v.x,v.y);}
struct Location {std::size_t panel;double fraction;};
bool before(Location a,Location b)
{return a.panel<b.panel || (a.panel==b.panel && a.fraction<b.fraction);}
bool same(Location a,Location b) {return a.panel==b.panel && a.fraction==b.fraction;}
Location canonical(const InterfaceSurfaceTransfer::Boundary& b,std::size_t count,bool closed)
{
    require(b.target.panelIndex>=0 && static_cast<std::size_t>(b.target.panelIndex)<count
        && b.target.nativeGraphFraction>=0 && b.target.nativeGraphFraction<=1,"invalid retained target coordinate");
    Location p{static_cast<std::size_t>(b.target.panelIndex),b.target.nativeGraphFraction};
    if (p.fraction==1) {++p.panel;p.fraction=0;if (closed && p.panel==count) p.panel=0;}
    return p;
}
}
void InterfaceSurfaceTransfer::validateFor(const std::shared_ptr<const LevelSetHermiteGeometry::Snapshot>& source,
    const std::shared_ptr<const LevelSetHermiteGeometry::Snapshot>& target,const InterfaceStepTime& time) const
{
    require(source && target && source==source_ && target==target_ && time.index==time_.index
        && time.begin==time_.begin && time.end==time_.end && time.dt==time_.dt,
        "source, target or physical step does not belong to this material transfer");
}
std::shared_ptr<const InterfaceSurfaceTransfer> InterfaceSurfaceTransfer::create(const InterfaceState& source,
    const InterfaceState& target,const InterfaceStepTime& time,const std::function<Vector2D(Vector2D)>& forward)
{
    require(forward && source.hasConservativePhase() && target.hasConservativeGeometry(),"conservative physical step required");
    require(source.progress().completedSteps+1==time.index && source.progress().endTime==time.begin
        && target.progress().completedSteps==time.index && target.progress().endTime==time.end,"step progress mismatch");
    const auto& old=source.topology().surfaces();const auto& next=target.topology().surfaces();
    require(old.size()==next.size() && source.topology().events().size()==target.topology().events().size(),
        "topology events require a separate material closure");
    auto result=std::shared_ptr<InterfaceSurfaceTransfer>(new InterfaceSurfaceTransfer);
    result->source_=source.geometry();result->predictor_=target.conservativeGeometry().predictor;
    result->target_=target.geometry();result->time_=time;
    for (const auto& previous:old) {
        const auto found=std::find_if(next.begin(),next.end(),[&](const auto& s){return s.id==previous.id;});
        require(found!=next.end() && found->liquidId==previous.liquidId && found->parents.size()==1
            && found->parents.front()==previous.id,"surface has no unique unchanged material parent");
        const auto& oldChain=previous.geometry.nativeChain();const auto& newChain=found->geometry.nativeChain();
        require(oldChain.usesSource(*result->source_) && newChain.usesSource(*result->target_)
            && oldChain.closed()==newChain.closed() && oldChain.axisymmetric()==newChain.axisymmetric(),"retained surface mismatch");
        const bool closed=oldChain.closed();
        using Kind=LevelSetInterfacePanels::EndpointKind;
        if (!closed) for (bool end:{false,true}) {
            const auto kind=end?oldChain.endEndpointKind():oldChain.beginEndpointKind();
            require(kind==(end?newChain.endEndpointKind():newChain.beginEndpointKind())
                && ((kind==Kind::Axis && oldChain.axisymmetric())
                    || (kind==Kind::SymmetryPlane && oldChain.endpointSymmetryPlane(end) && newChain.endpointSymmetryPlane(end)
                        && *oldChain.endpointSymmetryPlane(end)==*newChain.endpointSymmetryPlane(end))
                    || (kind==Kind::DeclaredBoundary && target.conservativeGeometry().fixedWallTraceVertices>0)),
                "material endpoint requires an unchanged axis, symmetry section or explicit fixed wall law");
        }
        Surface surface;surface.sourceId=previous.id;surface.targetId=found->id;surface.liquidId=previous.liquidId;
        surface.closed=closed;surface.axisymmetric=oldChain.axisymmetric();
        for (const auto& panel:previous.geometry.panels()) surface.sourceMeasures.push_back(panel.physicalMeasure);
        for (const auto& panel:found->geometry.panels()) surface.targetMeasures.push_back(panel.physicalMeasure);
        const auto count=surface.sourceMeasures.size(),newCount=surface.targetMeasures.size();
        require(count>0 && newCount>0,"empty material surface");
        const auto spacing=newChain.minimumActiveSpacing();const double h=std::min(spacing.x,spacing.y);
        for (std::size_t i=0;i<count+(closed?0:1);++i) {
            Boundary boundary;const bool end=i==count,endpoint=!closed && (i==0 || end);
            const auto kind=end?oldChain.endEndpointKind():oldChain.beginEndpointKind();
            const bool axis=endpoint && kind==Kind::Axis,fixed=endpoint && kind==Kind::DeclaredBoundary;
            boundary.sourceFrame=oldChain.frameAt(static_cast<int>(end?count-1:i),end?1:0);
            const auto& anchor=boundary.sourceFrame;
            require(anchor.hasCellChartPoint,"material source lost its native chart");
            boundary.physicalPoint=forward(anchor.point);boundary.physicalDisplacement=boundary.physicalPoint-anchor.point;
            const double motion=length(boundary.physicalDisplacement);
            require(std::isfinite(motion) && motion<=.5*h,"material trajectory exceeds half a cell");
            if (axis) require(anchor.point.x==0 && boundary.physicalPoint.x==0,"material axis point left the axis");
            else if (surface.axisymmetric) require(boundary.physicalPoint.x>0,"interior material point reached the axis");
            if (fixed) {
                const auto point=newChain.frameAt(static_cast<int>(end?newCount-1:0),end?1:0).point;
                require(boundary.physicalPoint.x==anchor.point.x && boundary.physicalPoint.y==anchor.point.y
                    && point.x==anchor.point.x && point.y==anchor.point.y,"fixed material endpoint moved");
            }
            if(endpoint && kind==Kind::SymmetryPlane) {
                const auto* plane=oldChain.endpointSymmetryPlane(end);
                require(plane && (plane->normalAxis?boundary.physicalPoint.y:boundary.physicalPoint.x)==plane->coordinate,
                    "material trajectory left its symmetry section");
            }
            const LevelSetInterfacePanels::NativeProjectionQuery query{anchor.nativeCellIndex,anchor.cellChartPoint,
                boundary.physicalDisplacement};
            try {
                boundary.target=endpoint?newChain.projectLocalEndpoint(query,static_cast<int>(end?newCount-1:0),!end):
                    newChain.projectLocal(query);
            } catch(const std::exception& error) {
                std::ostringstream message;message<<std::setprecision(17)<<error.what()
                    <<"; material step="<<time.index<<", surface="<<previous.id<<", boundary="<<i
                    <<", endpoint="<<endpoint<<", source_cell="<<anchor.nativeCellIndex
                    <<", source=("<<anchor.point.x<<','<<anchor.point.y<<"), physical=("
                    <<boundary.physicalPoint.x<<','<<boundary.physicalPoint.y<<')';
                // Failure-only observation of the rounded physical query. It
                // never replaces the exact anchored local query or its gates.
                try {
                    const auto nearest=newChain.projectedFrame(boundary.physicalPoint);
                    message<<"; diagnostic_global_distance="<<nearest.distance<<", global_graph="<<nearest.panelIndex
                        <<", global_point=("<<nearest.frame.point.x<<','<<nearest.frame.point.y<<')';
                } catch(const std::exception& diagnostic) {message<<"; global_diagnostic_failed="<<diagnostic.what();}
                const int cell=anchor.nativeCellIndex;
                if(source.hasConservativePhase() && target.hasConservativePhase())
                    message<<"; source_C="<<source.conservativePhase().liquidVolumes()[cell]/source.conservativePhase().capacities()[cell]
                        <<", target_C="<<target.conservativePhase().liquidVolumes()[cell]/target.conservativePhase().capacities()[cell];
                const auto box=result->target_->nativeCellBounds(cell);
                for(int corner=0;corner<4;++corner) {
                    const Vector2D p{corner&1?box[1].x:box[0].x,corner&2?box[1].y:box[0].y};
                    const auto a=result->source_->scalarJetInNativeCell(cell,p),b=result->target_->scalarJetInNativeCell(cell,p);
                    message<<"; corner=("<<p.x<<','<<p.y<<"), old_jet=("<<a.value<<','<<a.dx<<','<<a.dy
                        <<','<<a.dxy<<"), new_jet=("<<b.value<<','<<b.dx<<','<<b.dy<<','<<b.dxy<<')';
                }
                throw std::runtime_error(message.str());
            }
            const auto& projection=boundary.target;
            require(std::isfinite(projection.distance) && projection.distance<=.25*h
                && std::isfinite(projection.kktResidual) && projection.kktResidual<=1e-9*h
                && std::isfinite(projection.frame.levelSetResidual) && projection.frame.levelSetResidual<=1e-9*h,
                "numerical relocation exceeds retained material projection bounds");
            require(anchor.normal.x*projection.frame.normal.x+anchor.normal.y*projection.frame.normal.y>0,
                "material surface normal reversed");
            // Core already guards nonadjacent branches of this chain. Check
            // competing surfaces too, retaining the sampled-search limitation.
            for (const auto& other:next) if (other.id!=found->id)
                for (const auto& candidate:other.geometry.nativeChain().loadProjectionCandidates(boundary.physicalPoint))
                    require(candidate.distance>projection.distance+.01*h,"material relocation competes with another surface");
            boundary.numericalDisplacement=projection.frame.point-boundary.physicalPoint;
            surface.maxPhysicalDisplacement=std::max(surface.maxPhysicalDisplacement,motion);
            surface.maxNumericalDisplacement=std::max(surface.maxNumericalDisplacement,projection.distance);
            surface.boundaries.push_back(std::move(boundary));
        }
        if (closed) surface.boundaries.push_back(surface.boundaries.front());
        std::vector<Location> locations;
        for (const auto& boundary:surface.boundaries) locations.push_back(canonical(boundary,newCount,closed));
        int winding=0;
        for (std::size_t i=0;i<count;++i) {
            require(!same(locations[i],locations[i+1]),"material boundaries collapsed");
            if (!before(locations[i],locations[i+1])) ++winding;
        }
        if (closed) require(winding==1,"material mapping does not have one positive winding");
        else require(winding==0 && same(locations.front(),{0,0}) && same(locations.back(),{newCount,0}),
            "open material mapping folded or lost its boundary endpoints");
        std::vector<long double> covered(newCount,0);
        for (std::size_t i=0;i<count;++i) {
            Row row;auto current=locations[i];const auto end=locations[i+1];bool finished=false;
            for (std::size_t visited=0;visited<=newCount;++visited) {
                require(current.panel<newCount,"material traversal escaped target panels");
                const bool last=current.panel==end.panel && end.fraction>current.fraction;
                const double upper=last?end.fraction:1;
                long double measure=0;
                if (current.fraction==0 && upper==1) measure=surface.targetMeasures[current.panel];
                else for (const auto& point:newChain.quadrature(static_cast<int>(current.panel),current.fraction,upper))
                    measure+=point.weight*(surface.axisymmetric?CartesianMeasure::twoPi*point.frame.point.x:1.L);
                require(std::isfinite(measure) && measure>0,"unresolved mapped material interval");
                row.overlaps.push_back({current.panel,measure,0});row.mappedMeasure+=measure;covered[current.panel]+=measure;
                if (last) {finished=true;break;}
                ++current.panel;current.fraction=0;if (closed && current.panel==newCount) current.panel=0;
                if (same(current,end)) {finished=true;break;}
            }
            require(finished && std::isfinite(row.mappedMeasure) && row.mappedMeasure>0,"unclosed material interval");
            long double shares=0;
            for (auto& overlap:row.overlaps) {
                overlap.share=overlap.measure/row.mappedMeasure;
                require(overlap.share>0 && overlap.share<=1,"invalid conservative overlap share");shares+=overlap.share;
            }
            require(std::abs(shares-1)<=128*std::numeric_limits<long double>::epsilon()*row.overlaps.size(),
                "material overlap weights do not conserve an integrated quantity");
            surface.rows.push_back(std::move(row));
        }
        for (std::size_t j=0;j<newCount;++j) {
            require(surface.targetMeasures[j]>0 && std::isfinite(surface.targetMeasures[j]),"invalid target surface measure");
            surface.maxRelativeMeasurePartitionError=std::max(surface.maxRelativeMeasurePartitionError,
                static_cast<double>(std::abs(covered[j]/surface.targetMeasures[j]-1)));
        }
        result->surfaces_.push_back(std::move(surface));
    }
    return result;
}
