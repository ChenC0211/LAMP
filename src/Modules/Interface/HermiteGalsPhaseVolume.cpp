#include "Modules/Interface/HermiteGalsPhaseVolume.h"

#include "Modules/Interface/InterfaceCellMeasure.h"
#include "Core/Mesh/FaceSegmentGraph.h"
#include "Core/Mesh/NativeLeafDomain.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>

namespace LAMP::Interface
{
namespace
{
void require(bool ok, HermiteGalsIdentity id, const std::string& message)
{
    if(!ok)throw std::invalid_argument("[Interface/HermiteGalsPhaseVolume source="
        +std::to_string(id.value())+"] "+message);
}
struct Sum
{
    long double value=0, correction=0;
    void add(long double x) { const auto a=x-correction,b=value+a; correction=(b-value)-a; value=b; }
};
}
struct HermiteGalsPhaseVolumeResult::Data
{
    explicit Data(const HermiteGalsState& input):source(input){}
    HermiteGalsState source;
    std::vector<HermiteGalsCellVolume> volumes;
    std::vector<HermiteGalsPhaseVolumeCell> cells;
    HermiteGalsPhaseVolumeReport report;
    std::optional<InterfaceCellMeasure::NativeTraceLedger> traceLedger;
};
HermiteGalsPhaseVolumeResult::HermiteGalsPhaseVolumeResult(std::shared_ptr<const Data> data)
    :data_(std::move(data)){}
const HermiteGalsPhaseVolumeResult::Data& HermiteGalsPhaseVolumeResult::checked() const
{ require(static_cast<bool>(data_),{},"access: moved-from result");return *data_; }
HermiteGalsIdentity HermiteGalsPhaseVolumeResult::sourceIdentity() const {return checked().source.identity();}
HermiteGalsTimeLayer HermiteGalsPhaseVolumeResult::timeLayer() const {return checked().source.timeLayer();}
HermiteGalsGeometryView HermiteGalsPhaseVolumeResult::geometry() const {return checked().source.geometry();}
const std::vector<HermiteGalsCellVolume>& HermiteGalsPhaseVolumeResult::phaseVolumes() const {return checked().volumes;}
const std::vector<HermiteGalsPhaseVolumeCell>& HermiteGalsPhaseVolumeResult::cells() const {return checked().cells;}
const HermiteGalsPhaseVolumeReport& HermiteGalsPhaseVolumeResult::report() const {return checked().report;}
const InterfaceCellMeasure::NativeTraceLedger& HermiteGalsPhaseVolumeResult::nativeTraceLedger() const
{
    const auto& d=checked();require(d.traceLedger.has_value(),d.source.identity(),"native trace ledger unavailable on historical port");
    return *d.traceLedger;
}
HermiteGalsCandidate HermiteGalsPhaseVolumeResult::initializeCandidate(const QuadtreeMesh& mesh) const
{
    const auto& d=checked(); const auto t=d.source.timeLayer(); const auto id=d.source.identity();
    require(t.position==HermiteGalsTimeLayer::Position::Initial && t.step==0 && t.begin==0 && t.end==0,
        id,"initialize: requires fresh step-0 time-0 Initial template; restart/advancement forbidden");
    d.source.validateMesh(mesh);
    auto primary=d.source.primary(); primary.phaseVolumes=d.volumes;
    return HermiteGalsCandidate(mesh,t,primary);
}
HermiteGalsPhaseVolumeResult integrateHermiteGalsPhaseVolumes(const QuadtreeMesh& mesh,
    const HermiteGalsState& source, const HermiteGalsPhaseVolumeConfig& config)
{ return HermiteGalsPhaseVolumeResult::integrate(mesh,source,config,nullptr,nullptr,nullptr); }
HermiteGalsPhaseVolumeResult integrateHermiteGalsNativeTracePhaseVolumes(const QuadtreeMesh& mesh,
    const HermiteGalsState& source, const NativeLeafDomain& domain,
    const InterfacePhaseGeometry::MappedFacePartitions& partitions,
    const LevelSetInterfacePanels::NativeComponents& topology, const HermiteGalsPhaseVolumeConfig& config)
{ return HermiteGalsPhaseVolumeResult::integrate(mesh,source,config,&domain,&partitions,&topology); }
HermiteGalsPhaseVolumeResult HermiteGalsPhaseVolumeResult::integratePublicationTrace(const QuadtreeMesh& mesh,
    const HermiteGalsState& source,const NativeLeafDomain& domain,
    const InterfacePhaseGeometry::MappedFacePartitions& partitions,
    const LevelSetInterfacePanels::NativeComponents& topology,const HermiteGalsPhaseVolumeConfig& config,
    const Detail::PublicationTracePermit& permit)
{
    require(permit.mesh_==&mesh && permit.source_==&source && permit.domain_==&domain
        && permit.partitions_==&partitions && permit.topology_==&topology && permit.config_==&config,
        source.identity(),"publication trace permit differs from actual factory inputs");
    return integrate(mesh,source,config,&domain,&partitions,&topology,&permit);
}
HermiteGalsPhaseVolumeResult HermiteGalsPhaseVolumeResult::integrate(const QuadtreeMesh& mesh,
    const HermiteGalsState& source, const HermiteGalsPhaseVolumeConfig& config,
    const NativeLeafDomain* retainedDomain, const InterfacePhaseGeometry::MappedFacePartitions* retainedPartitions,
    const LevelSetInterfacePanels::NativeComponents* retainedTopology, const Detail::PublicationTracePermit* publicationPermit)
{
    const auto id=source.identity();
    try {
        source.validateMesh(mesh);
        const auto& p=source.primary();
        const bool conservative=p.transport==HermiteGalsTransport::ConservativeVolume;
        require((!retainedDomain && !retainedPartitions && !retainedTopology)
            || (!conservative && retainedDomain && retainedPartitions && retainedTopology),id,"integrate: invalid retained trace input/transport");
        const auto validTolerance=[](double v){return std::isfinite(v) && v>0 && v<1;};
        require(validTolerance(config.relativeLedgerTolerance) && validTolerance(config.relativeMomentEstimateTolerance),
            id,"integrate: relative tolerances must be finite and in (0,1)");
        require(!p.geometryCells.empty() && (!conservative || p.phaseVolumes.size()==p.geometryCells.size()),id,
            "integrate: every inventory cell requires scalar support");
        for(std::size_t k=0;conservative && k<p.geometryCells.size();++k)
            require(p.phaseVolumes[k].cell==p.geometryCells[k],id,"integrate: inventory/scalar support mismatch");
        const bool radial=p.measure==HermiteGalsMeasure::AxisymmetricVolume;
        auto owner=std::make_shared<const HermiteGalsGeometryView>(source.geometry());
        // Aliasing ownership retains the immutable view; no scalar copy or
        // externally writable shared_ptr is created.
        const std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> scalar(owner,&owner->snapshot());
        require(scalar->usesNativeLeaves() && (scalar->reconstruction()==LevelSetHermiteGeometry::Reconstruction::EnrichedQuartic||scalar->usesC2Macro()),
            id,"integrate: requires native quartic source");
        require(scalar->radialAxisSymmetry()==radial,id,"integrate: scalar/primary XY/RZ measure mismatch");
        std::optional<NativeLeafDomain> localDomain;
        if(!retainedDomain) {
            FaceSegmentGraph graph;
            graph.build(mesh,[&](const QuadtreeMesh&,const Cell& cell){
                return std::binary_search(p.geometryCells.begin(),p.geometryCells.end(),cell.index());});
            localDomain.emplace(NativeLeafDomain::build(mesh,graph,radial));
        }
        const auto& domain=retainedDomain?*retainedDomain:*localDomain;
        InterfaceCellMeasure::DiscoveredVolumes measured;
        std::optional<InterfaceCellMeasure::NativeTraceLedger> traceLedger;
        if(conservative)measured=InterfaceCellMeasure::discoverVolumes(mesh,domain,scalar,config.boundaryContacts,config.symmetryPlanes);
        else {
            if(retainedPartitions) {
                require(&retainedPartitions->source()==scalar.get(),id,"integrate: retained trace changed scalar owner");
                retainedPartitions->validate(mesh,domain);
            }
            auto trace=publicationPermit
                ? InterfaceCellMeasure::Detail::PublicationTraceVolumeAccess::measure(mesh,domain,*retainedPartitions,*retainedTopology,
                    config.boundaryContacts,config.symmetryPlanes,*publicationPermit)
                : retainedPartitions ? InterfaceCellMeasure::measureNativeTraceVolumes(mesh,domain,*retainedPartitions,*retainedTopology,
                    config.boundaryContacts,config.symmetryPlanes)
                : InterfaceCellMeasure::discoverNativeTraceVolumes(mesh,domain,scalar,config.boundaryContacts,config.symmetryPlanes);
            traceLedger=std::move(trace.ledger);measured=std::move(static_cast<InterfaceCellMeasure::DiscoveredVolumes&>(trace));
        }
        auto output=std::make_shared<HermiteGalsPhaseVolumeResult::Data>(source);
        output->traceLedger=std::move(traceLedger);
        auto& report=output->report; const auto& raw=measured.measure;
        report.components=measured.components;report.graphs=measured.graphs;report.faceRootIncidences=measured.faceRootIncidences;
        Sum capacitySum,rawSum,initialSum,inventorySum;
        for(int cell:p.geometryCells) {
            const auto fail=[&](bool ok,const char* why){require(ok,id,"integrate: cell="+std::to_string(cell)+" "+why);};
            const int row=domain.cellRowByIndex()[cell];fail(row>=0,"missing native volume owner");
            const double capacity=domain.cellVolumes()[row];const auto volume=raw.negativeVolumeByCell[cell];
            fail(std::isfinite(capacity) && capacity>0 && std::isfinite(volume),"nonfinite/invalid volume or capacity");
            const auto capacityDefect=std::abs(raw.fixedBoundary.capacityByCell[cell]-capacity)/capacity;
            fail(std::isfinite(capacityDefect) && capacityDefect<=config.relativeLedgerTolerance,"boundary capacity defect");
            report.maxRelativeCapacityDefect=std::max(report.maxRelativeCapacityDefect,capacityDefect);
            const int phase=measured.cellPhase[cell]; double liquid=0;
            if(phase==0) {
                fail(volume>=0 && volume<=capacity,"cut volume outside strict [0,capacity]; no clipping");
                liquid=static_cast<double>(volume);++report.cutCells;
            } else {
                fail(phase==-1 || phase==1,"invalid pure-phase certificate");
                liquid=phase<0?capacity:0;
                const auto defect=std::abs(volume-liquid)/capacity;
                fail(defect<=config.relativeLedgerTolerance,"pure certificate/raw integral disagreement");
                report.maxRelativePureDefect=std::max(report.maxRelativePureDefect,defect);
                if(phase<0)++report.liquidCells;else ++report.gasCells;
            }
            fail(std::isfinite(liquid) && liquid>=0 && liquid<=capacity,"unrepresentable initial volume");
            output->volumes.push_back({cell,liquid});output->cells.push_back({cell,phase,capacity,volume});
            capacitySum.add(capacity);rawSum.add(volume);initialSum.add(liquid);
            if(conservative){
                const auto inventory=p.phaseVolumes[output->volumes.size()-1].liquid;
                inventorySum.add(inventory);
                report.maxCellRelativeInventoryDifference=std::max(report.maxCellRelativeInventoryDifference,std::abs(volume-inventory)/capacity);
            }
        }
        require(std::isfinite(capacitySum.value) && capacitySum.value>0 && std::isfinite(rawSum.value)
            && std::isfinite(initialSum.value) && std::isfinite(inventorySum.value),id,"integrate: nonfinite totals");
        const auto boundaryClosure=std::abs(raw.fixedBoundary.assemblyClosure)/capacitySum.value;
        const auto momentEstimate=raw.momentErrorEstimate/capacitySum.value;
        require(std::isfinite(boundaryClosure) && boundaryClosure<=config.relativeLedgerTolerance,
            id,"integrate: global boundary ledger does not close");
        if(output->traceLedger) {
            const auto& ledger=*output->traceLedger;
            require(std::isfinite(ledger.fixedRadialClosure) && std::abs(ledger.fixedRadialClosure)/capacitySum.value<=config.relativeLedgerTolerance
                && ledger.maxRelativeGreenDefect<=config.relativeLedgerTolerance
                && ledger.maxRelativeCapacityDefect<=config.relativeLedgerTolerance,id,"integrate: native trace Green/capacity ledger does not close");
        }
        require(std::isfinite(momentEstimate) && momentEstimate>=0
            && momentEstimate<=config.relativeMomentEstimateTolerance,id,"integrate: moment error estimate exceeds tolerance");
        report.totalCapacity=capacitySum.value;report.totalRawLiquid=rawSum.value;report.totalInitialLiquid=initialSum.value;
        report.relativeBoundaryClosure=boundaryClosure;report.relativeMomentEstimate=momentEstimate;
        report.totalSourceInventory=inventorySum.value;
        report.geometricInventoryDifference=rawSum.value-inventorySum.value;
        if(!conservative){
            const auto unavailable=std::numeric_limits<long double>::quiet_NaN();
            report.totalSourceInventory=unavailable;report.geometricInventoryDifference=unavailable;
            report.maxCellRelativeInventoryDifference=unavailable;
        }
        return HermiteGalsPhaseVolumeResult(std::move(output));
    } catch(const std::exception& e) {
        throw std::invalid_argument("[Interface/HermiteGalsPhaseVolume source="+std::to_string(id.value())+"] "+e.what());
    }
}
}
