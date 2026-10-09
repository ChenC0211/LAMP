#include "Modules/Interface/HermiteGalsInventoryRegion.h"
#include "Modules/Interface/Detail/HermiteGalsRegionMeasure.h"
#include "Modules/Interface/Detail/HermiteGalsVolumeRoundoff.h"
#include "Core/Mesh/FaceSegmentGraph.h"
#include "Core/Mesh/NativeLeafDomain.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <utility>

namespace Quantum::Interface {
namespace {
constexpr auto roundoff=Detail::hermiteInventoryRoundoffBudget;
void require(bool ok,const std::string& why){if(!ok)throw std::invalid_argument("[HermiteGalsInventoryRegion] "+why);}
struct Sum {long double value=0,correction=0;void add(long double x){const auto a=x-correction,b=value+a;correction=(b-value)-a;value=b;}};
std::array<long double,2> mix(const HermiteGalsInventoryRegionCell& c,long double cap,long double liquid,long double gas){
    if(c.liquidWeight==1 && c.gasWeight==0)return {{liquid,gas}};
    const auto l=c.liquidWeight*liquid+c.gasWeight*gas;
    const auto g=(1-c.liquidWeight)*liquid+(1-c.gasWeight)*gas;
    require(std::isfinite(l) && std::isfinite(g) && l>=0 && g>=0 && std::abs(l+g-cap)<=roundoff*cap,"weighted phase capacity defect");
    const auto result=Detail::RegionMeasure::complement(cap,l,g);
    require(result[0]>=0 && result[1]>=0 && result[0]<=cap && result[1]<=cap,"weighted phase outside strict capacity");
    return result;
}
}
struct HermiteGalsInventoryRegionSource::Data {
    explicit Data(const HermiteGalsPhaseRegionSource& r):regions(r){}
    HermiteGalsPhaseRegionSource regions;
    std::vector<HermiteGalsInventoryRegionCell> cells;
    HermiteGalsInventoryRegionReport report;
};
HermiteGalsInventoryRegionSource::HermiteGalsInventoryRegionSource(const QuadtreeMesh& mesh,
    const HermiteGalsPhaseRegionSource& regions,const HermiteGalsInventoryRegionConfig& config){
    HermiteGalsIdentity id;int diagnosticCell=-1;
    try{
        regions.validateMesh(mesh);const auto& source=regions.source();id=source.identity();const auto& p=source.primary();
        require(std::isfinite(config.relativeGeometryInventoryTolerance) && config.relativeGeometryInventoryTolerance>0
            && config.relativeGeometryInventoryTolerance<1,"invalid geometric inventory tolerance");
        require(p.geometryCells.size()==static_cast<std::size_t>(mesh.cellIndexNum()) && p.phaseVolumes.size()==p.geometryCells.size(),
            "requires complete native geometry and inventory");
        FaceSegmentGraph graph;graph.build(mesh);const auto domain=NativeLeafDomain::build(mesh,graph,p.measure==HermiteGalsMeasure::AxisymmetricVolume);
        auto data=std::make_shared<Data>(regions);const auto& scalar=source.geometry().snapshot();
        data->report.relativeGeometryInventoryTolerance=config.relativeGeometryInventoryTolerance;
        Sum total,represented,l1;
        for(int cell=0;cell<mesh.cellIndexNum();++cell){diagnosticCell=cell;
            require(p.geometryCells[cell]==cell && p.phaseVolumes[cell].cell==cell,"unordered native inventory");
            // The immutable source and this mesh were checked once above.
            // Reuse the same private query kernel without N full mesh scans.
            const auto bounds=scalar.nativeCellBounds(cell);const auto measured=regions.measureValidatedMesh(mesh,bounds[0],bounds[1]);
            require(measured.leaves.size()==1 && measured.leaves[0].cell==cell,"whole-cell region ownership mismatch");
            const long double cap=domain.cellVolumes()[cell],m=p.phaseVolumes[cell].liquid;
            require(std::isfinite(m) && m>=0 && m<=cap,"source inventory outside strict capacity");
            const auto capacityDefect=std::abs(measured.capacity-cap)/cap;
            require(capacityDefect<=roundoff,"native/geometric capacity exceeds roundoff budget");
            const auto geometryDefect=std::abs(measured.liquid-m)/cap;
            require(geometryDefect<=config.relativeGeometryInventoryTolerance,"source geometry/inventory gap exceeds configured tolerance");
            const auto a=m/cap,v=measured.liquid/measured.capacity;
            require(std::isfinite(a) && std::isfinite(v) && a>=0 && a<=1 && v>=0 && v<=1,"invalid phase fractions");
            HermiteGalsInventoryRegionCell row{cell,cap,measured.capacity,measured.liquid,measured.gas,m};
            if(a<=v){row.liquidWeight=v>0?a/v:1;row.gasWeight=0;}
            else{row.liquidWeight=1;row.gasWeight=(a-v)/(1-v);}
            require(std::isfinite(row.liquidWeight) && std::isfinite(row.gasWeight) && row.gasWeight>=0
                && row.gasWeight<=row.liquidWeight && row.liquidWeight<=1,"phase-content weights outside strict bounds");
            row.representedInventory=mix(row,measured.capacity,measured.liquid,measured.gas)[0];
            const auto representationDefect=std::abs(row.representedInventory-m)/cap;
            require(representationDefect<=roundoff,"whole-cell inventory representation exceeds roundoff budget");
            auto& report=data->report;report.adjustedCells+=row.liquidWeight!=1 || row.gasWeight!=0;
            report.maxRelativeCapacityDefect=std::max(report.maxRelativeCapacityDefect,capacityDefect);
            report.maxRelativeGeometryDefect=std::max(report.maxRelativeGeometryDefect,geometryDefect);
            report.maxRelativeRepresentationDefect=std::max(report.maxRelativeRepresentationDefect,representationDefect);
            total.add(m);represented.add(row.representedInventory);l1.add(std::abs(m-measured.liquid));data->cells.push_back(row);
        }
        data->report.inventory=total.value;data->report.representedInventory=represented.value;data->report.geometryInventoryL1=l1.value;
        data_=std::move(data);
    }catch(const std::exception& e){throw std::invalid_argument("[HermiteGalsInventoryRegion source="+std::to_string(id.value())
        +" cell="+std::to_string(diagnosticCell)+"] "+e.what());}
}
const HermiteGalsInventoryRegionSource::Data& HermiteGalsInventoryRegionSource::checked()const{
    require(bool(data_),"moved-from inventory region source");return *data_;
}
const HermiteGalsState& HermiteGalsInventoryRegionSource::source()const{return checked().regions.source();}
const HermiteGalsPhaseRegionSource& HermiteGalsInventoryRegionSource::regions()const{return checked().regions;}
const std::vector<HermiteGalsInventoryRegionCell>& HermiteGalsInventoryRegionSource::cells()const{return checked().cells;}
const HermiteGalsInventoryRegionReport& HermiteGalsInventoryRegionSource::report()const{return checked().report;}
void HermiteGalsInventoryRegionSource::validateMesh(const QuadtreeMesh& mesh)const{checked().regions.validateMesh(mesh);}
HermiteGalsInventoryRegionMeasure HermiteGalsInventoryRegionSource::measure(const QuadtreeMesh& mesh,Vector2D lo,Vector2D hi)const{
    const auto& data=checked();HermiteGalsInventoryRegionMeasure result;result.geometry=data.regions.measure(mesh,lo,hi);
    Sum negative,positive,l1;bool identical=true;
    for(const auto& leaf:result.geometry.leaves){
        const auto& c=data.cells.at(leaf.cell);const auto values=mix(c,leaf.capacity,leaf.liquid,leaf.gas);
        negative.add(values[0]);positive.add(values[1]);l1.add(std::abs(values[0]-leaf.liquid));
        identical=identical && c.liquidWeight==1 && c.gasWeight==0;
    }
    require(std::isfinite(negative.value) && std::isfinite(positive.value) && negative.value>=0 && positive.value>=0
        && std::abs(negative.value+positive.value-result.geometry.capacity)<=roundoff*result.geometry.capacity,"assembled inventory region capacity defect");
    const auto phases=identical?std::array<long double,2>{{result.geometry.liquid,result.geometry.gas}}:
        Detail::RegionMeasure::complement(result.geometry.capacity,negative.value,positive.value);
    result.liquid=phases[0];result.gas=phases[1];result.signedCorrection=result.liquid-result.geometry.liquid;result.l1Correction=l1.value;
    require(result.liquid>=0 && result.gas>=0 && result.liquid<=result.geometry.capacity && result.gas<=result.geometry.capacity,
        "inventory region outside strict capacity");
    return result;
}
}
