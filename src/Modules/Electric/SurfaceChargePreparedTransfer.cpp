#include "Modules/Electric/SurfaceChargeField.h"
#include "Core/Geometry/CellChartPoint.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace LAMP::Electric
{
namespace
{
void require(bool ok,const char* why)
{if (!ok) throw std::invalid_argument(std::string("[SurfaceChargePreparedTransfer] ")+why);}
void samePanel(const LevelSetInterfacePanels::Panel& a,const LevelSetInterfacePanels::Panel& b)
{
    require(a.cellIndex==b.cellIndex && a.coordinateIsX==b.coordinateIsX
        && a.measure==b.measure && a.physicalMeasure==b.physicalMeasure
        && a.beginFrame.hasCellChartPoint && b.beginFrame.hasCellChartPoint
        && a.endFrame.hasCellChartPoint && b.endFrame.hasCellChartPoint
        && Geometry::cellChartEqual(a.beginFrame.cellChartPoint,b.beginFrame.cellChartPoint)
        && Geometry::cellChartEqual(a.endFrame.cellChartPoint,b.endFrame.cellChartPoint),
        "source panel partition/orientation differs from charge storage");
}
}
SurfaceChargeMaterialRemap SurfaceChargeField::remapPrepared(const QuadtreeMesh& mesh,
    const SurfaceChargePreparedTransfer& plan,const std::vector<std::vector<double>>& history) const
{
    const auto timeScale=std::max({std::abs(plan.sourceTime),std::abs(plan.targetTime),plan.dt});
    require(plan.step>0 && std::isfinite(plan.sourceTime) && std::isfinite(plan.targetTime)
        && std::isfinite(plan.dt) && plan.dt>0 && plan.targetTime>plan.sourceTime
        && std::abs((plan.targetTime-plan.sourceTime)-plan.dt)<=16*std::numeric_limits<double>::epsilon()*timeScale,
        "prepared transfer has no consistent physical interval");
    require(source_ && panels_.usesNativeLeaves() && plan.source && plan.target
        && plan.source->usesNativeLeaves() && plan.target->usesNativeLeaves()
        && plan.source->usesSource(*source_),"prepared transfer source does not own this charge geometry");
    require(plan.source->panels().size()==panels_.panels().size()
        && plan.source->closed()==panels_.closed() && plan.target->closed()==panels_.closed()
        && plan.source->axisymmetric()==panels_.axisymmetric() && plan.target->axisymmetric()==panels_.axisymmetric(),
        "prepared transfer changed surface topology/metric");
    require(panels_.boundaryContacts().empty() && plan.source->boundaryContacts().empty()
        && plan.target->boundaryContacts().empty(),"prepared contacts require an explicit material boundary law");
    require(panels_.nativeChain().symmetryPlanes()==plan.source->nativeChain().symmetryPlanes()
        && panels_.nativeChain().symmetryPlanes()==plan.target->nativeChain().symmetryPlanes(),
        "prepared transfer changed reflection sections");
    if (!panels_.closed()) for(bool end:{false,true}) {
        using Kind=LevelSetInterfacePanels::EndpointKind;
        const auto kind=end?panels_.endEndpointKind():panels_.beginEndpointKind();
        require((end?plan.source->endEndpointKind():plan.source->beginEndpointKind())==kind
            && (end?plan.target->endEndpointKind():plan.target->beginEndpointKind())==kind,
            "prepared transfer changed endpoint kind");
        const auto* plane=panels_.nativeChain().endpointSymmetryPlane(end);
        const auto* sourcePlane=plan.source->nativeChain().endpointSymmetryPlane(end);
        const auto* targetPlane=plan.target->nativeChain().endpointSymmetryPlane(end);
        require((kind==Kind::Axis && panels_.axisymmetric())
            || (kind==Kind::SymmetryPlane && plane && sourcePlane && targetPlane
                && *plane==*sourcePlane && *plane==*targetPlane
                && (!panels_.axisymmetric() || plane->normalAxis==1)),
            "prepared open transfer requires retained axes or physical reflection sections");
    }
    for (std::size_t i=0;i<panels_.panels().size();++i) samePanel(panels_.panels()[i],plan.source->panels()[i]);
    validateNativeRemapGrid(mesh,*plan.target->nativeChain().source());
    require(std::isfinite(plan.maxMaterialDisplacement) && plan.maxMaterialDisplacement>=0
        && std::isfinite(plan.maxProjectionDistance) && plan.maxProjectionDistance>=0,"invalid material diagnostics");
    const auto spacing=plan.target->nativeChain().minimumActiveSpacing();const double h=std::min(spacing.x,spacing.y);
    require(plan.maxMaterialDisplacement<=.5*h && plan.maxProjectionDistance<=.25*h,
        "prepared transfer exceeds material motion/projection bounds");
    SurfaceChargeMaterialRemap result;
    result.field=std::make_unique<SurfaceChargeField>(mesh,*plan.target,0.0);
    const auto& targetMeasure=result.field->charge_.panelMeasure();
    auto staged=SurfaceChargeTransport::fromIntegratedCharge(charge_.panelMeasure(),charge_.panelCharge(),panels_.closed());
    result.report.transport=staged.remapWeights(plan.weights,targetMeasure);
    require(result.report.transport.success,result.report.transport.reason.c_str());
    result.field->assignIntegratedCharge(staged.panelCharge());
    for (const auto& values:history) {
        auto old=SurfaceChargeTransport::fromIntegratedCharge(charge_.panelMeasure(),values,panels_.closed());
        auto report=old.remapWeights(plan.weights,targetMeasure);
        require(report.success,report.reason.c_str());
        result.history.push_back(old.panelCharge());result.historyReports.push_back(std::move(report));
    }
    result.report.remapped=true;result.report.previousPanelCount=panels_.panels().size();
    result.report.remapConservationResidual=result.report.transport.conservationResidual;
    result.report.maxMaterialDisplacement=plan.maxMaterialDisplacement;
    result.report.maxProjectionDistance=plan.maxProjectionDistance;
    return result;
}
SurfaceChargeMaterialRemap SurfaceChargeField::remapSectionDetachment(const QuadtreeMesh& mesh,
    const SurfaceChargeSectionDetachment& plan,const std::vector<std::vector<double>>& history) const
{
    require(plan.step>=0 && std::isfinite(plan.time) && plan.time>=0
        && std::isfinite(plan.beginZ) && std::isfinite(plan.upperSection) && plan.beginZ<plan.upperSection,
        "invalid section-detachment clock or local interval");
    require(source_ && panels_.usesNativeLeaves() && plan.source && plan.target
        && plan.source->usesNativeLeaves() && plan.target->usesNativeLeaves()
        && plan.source->usesSource(*source_),"section-detachment source does not own this charge geometry");
    require(!panels_.closed() && !plan.source->closed() && !plan.target->closed()
        && panels_.axisymmetric() && plan.source->axisymmetric() && plan.target->axisymmetric()
        && panels_.boundaryContacts().empty() && plan.source->boundaryContacts().empty()
        && plan.target->boundaryContacts().empty(),"section detachment requires open RZ surfaces without contacts");
    require(plan.source->panels().size()==panels_.panels().size()
        && plan.weights.size()==panels_.panels().size(),"section-detachment source partition differs from charge storage");
    for(std::size_t i=0;i<panels_.panels().size();++i)samePanel(panels_.panels()[i],plan.source->panels()[i]);
    validateNativeRemapGrid(mesh,*plan.target->nativeChain().source());
    require(plan.cells[0]>=0 && plan.cells[1]>=0 && plan.cells[0]!=plan.cells[1],"invalid local event cells");
    const auto first=source_->nativeCellBounds(plan.cells[0]),last=source_->nativeCellBounds(plan.cells[1]);
    require(first[0].x==0 && last[0].x==0 && first[1].x==last[1].x
        && first[0].y==plan.beginZ && first[1].y==last[0].y && last[1].y==plan.upperSection,
        "section-detachment cells do not form the declared local axis interval");
    for(int cell:plan.cells) {
        const auto a=source_->nativeCellBounds(cell),b=plan.target->nativeChain().source()->nativeCellBounds(cell);
        require(a[0].x==b[0].x && a[0].y==b[0].y && a[1].x==b[1].x && a[1].y==b[1].y,
            "section-detachment transfer changed mesh coordinates");
    }
    const auto& old=panels_.nativeChain();const auto& next=plan.target->nativeChain();
    int changed=0;
    for(bool end:{false,true}) {
        using Kind=LevelSetInterfacePanels::EndpointKind;
        const auto oldKind=end?old.endEndpointKind():old.beginEndpointKind();
        const auto newKind=end?next.endEndpointKind():next.beginEndpointKind();
        const auto* a=old.endpointSymmetryPlane(end);const auto* b=next.endpointSymmetryPlane(end);
        require(oldKind==Kind::SymmetryPlane && a && a->normalAxis==1,"event source is not an axial reflection column");
        if(a->coordinate==plan.upperSection) {
            ++changed;require(newKind==Kind::Axis && !b,"detached endpoint did not become an axis pole");
            const auto point=next.frameAt(end?static_cast<int>(next.graphs().size())-1:0,end?1:0).point;
            require(point.x==0 && point.y>plan.beginZ && point.y<plan.upperSection,"new axis pole escaped the local event interval");
        } else require(newKind==oldKind && b && *a==*b,"event changed the surviving reflection section");
    }
    require(changed==1,"event must detach exactly one upper reflection section");
    const auto local=[&](const auto& panel){return panel.cellIndex==plan.cells[0] || panel.cellIndex==plan.cells[1];};
    for(std::size_t i=0;i<plan.weights.size();++i) {
        const auto& from=panels_.panels()[i];const auto& row=plan.weights[i];
        require(!row.empty(),"event transfer omitted a source panel");
        for(const auto& weight:row) {
            require(weight.targetPanel<plan.target->panels().size(),"event target panel is outside the retained surface");
            const auto& to=plan.target->panels()[weight.targetPanel];
            if(local(from))require(local(to),"event redistributed local charge outside the cap");
            else {
                require(!local(to) && row.size()==1 && weight.share==1,"event altered a remote material integral");
                samePanel(from,to);
            }
        }
    }
    SurfaceChargeMaterialRemap result;
    result.field=std::make_unique<SurfaceChargeField>(mesh,*plan.target,0.0);
    const auto& measures=result.field->charge_.panelMeasure();
    auto mapped=SurfaceChargeTransport::fromIntegratedCharge(charge_.panelMeasure(),charge_.panelCharge(),false);
    result.report.transport=mapped.remapWeights(plan.weights,measures);
    require(result.report.transport.success,result.report.transport.reason.c_str());
    result.field->assignIntegratedCharge(mapped.panelCharge());
    for(const auto& values:history) {
        auto oldHistory=SurfaceChargeTransport::fromIntegratedCharge(charge_.panelMeasure(),values,false);
        auto report=oldHistory.remapWeights(plan.weights,measures);
        require(report.success,report.reason.c_str());
        result.history.push_back(oldHistory.panelCharge());result.historyReports.push_back(std::move(report));
    }
    result.report.remapped=true;result.report.previousPanelCount=panels_.panels().size();
    result.report.remapConservationResidual=result.report.transport.conservationResidual;
    result.report.sectionDetached=true;
    result.report.sectionRemapConservationResidual=result.report.transport.conservationResidual;
    return result;
}
std::unique_ptr<SurfaceChargeField> SurfaceChargeField::advancePrepared(const QuadtreeMesh& mesh,
    const SurfaceChargePreparedTransfer& plan,const std::vector<double>& cellCurrent,
    double currentScale,double relativeBulkTolerance,double dt,SurfaceChargeFieldStep& report) const
{
    require(panels_.usesNativeLeaves(),"prepared current step requires native panels");
    require(dt==plan.dt,"current update and prepared transfer use different time steps");
    SurfaceChargeField staged(mesh,panels_,0.0);staged.assignIntegratedCharge(charge_.panelCharge());
    const auto current=staged.advanceNormalCurrent(cellCurrent,currentScale,relativeBulkTolerance,dt);
    require(current.transport.success,current.transport.reason.c_str());
    auto mapped=staged.remapPrepared(mesh,plan);
    auto result=mapped.report;
    result.unassignedBulkCurrent=current.unassignedBulkCurrent;result.absoluteBulkCurrent=current.absoluteBulkCurrent;
    result.transport=current.transport;
    result.transport.totalChargeAfter=mapped.report.transport.totalChargeAfter;
    result.transport.conservationResidual+=mapped.report.transport.conservationResidual;
    result.transport.relativeConservationResidual=result.transport.conservationScale==0?0:
        std::abs(result.transport.conservationResidual)/result.transport.conservationScale;
    result.transport.maxAbsDensity=mapped.report.transport.maxAbsDensity;
    report=std::move(result);return std::move(mapped.field);
}
}
