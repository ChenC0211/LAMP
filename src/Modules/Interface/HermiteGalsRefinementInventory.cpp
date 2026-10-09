#include "Modules/Interface/HermiteGalsRefinementInventory.h"

#include "Modules/Interface/Detail/HermiteGalsVolumeRoundoff.h"
#include "Modules/Interface/InterfaceCellMeasure.h"
#include "Core/Geometry/LevelSetHermiteGeometry.h"
#include "Core/Mesh/FaceSegmentGraph.h"
#include "Core/Mesh/NativeLeafDomain.h"
#include "Core/Mesh/NativeMeshIdentity.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <utility>

namespace Quantum::Interface
{
namespace
{
void require(bool ok,const std::string& reason)
{ if(!ok)throw std::invalid_argument(reason); }
struct Sum
{
    long double value=0,correction=0;
    void add(long double x){const auto a=x-correction,b=value+a;correction=(b-value)-a;value=b;}
};
long double canonical(const HermiteGalsPhaseVolumeCell& cell)
{ return cell.certifiedPhase<0?cell.capacity:cell.certifiedPhase>0?0.L:cell.rawLiquid; }
}
struct HermiteGalsRefinementInventory::Data
{
    Data(const QuadtreeMesh& oldMesh,const HermiteGalsState& old,const QuadtreeMesh& target)
        : sourceMesh(oldMesh),targetMesh(target),source(old) {}
    NativeMeshIdentity sourceMesh,targetMesh;
    HermiteGalsState source;
    std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> scalar;
    std::vector<HermiteGalsCellVolume> volumes;
    std::vector<HermiteGalsPhaseVolumeCell> cells;
    std::vector<int> parents;
    std::vector<HermiteGalsRefinementParentBalance> balances;
    HermiteGalsRefinementInventoryReport report;
};
HermiteGalsRefinementInventory::HermiteGalsRefinementInventory(std::shared_ptr<const Data> data):data_(std::move(data)){}
const HermiteGalsRefinementInventory::Data& HermiteGalsRefinementInventory::checked() const
{ require(bool(data_),"[HermiteGalsRefinementInventory] moved-from result");return *data_; }
HermiteGalsIdentity HermiteGalsRefinementInventory::sourceIdentity() const{return checked().source.identity();}
HermiteGalsTimeLayer HermiteGalsRefinementInventory::timeLayer() const{return checked().source.timeLayer();}
const HermiteGalsState& HermiteGalsRefinementInventory::source() const{return checked().source;}
const LevelSetHermiteGeometry::Snapshot& HermiteGalsRefinementInventory::geometry() const{return *checked().scalar;}
const std::vector<HermiteGalsCellVolume>& HermiteGalsRefinementInventory::phaseVolumes() const{return checked().volumes;}
const std::vector<HermiteGalsPhaseVolumeCell>& HermiteGalsRefinementInventory::cells() const{return checked().cells;}
const std::vector<int>& HermiteGalsRefinementInventory::parentCells() const{return checked().parents;}
const std::vector<HermiteGalsRefinementParentBalance>& HermiteGalsRefinementInventory::parentBalances() const{return checked().balances;}
const HermiteGalsRefinementInventoryReport& HermiteGalsRefinementInventory::report() const{return checked().report;}
void HermiteGalsRefinementInventory::validateMeshes(const QuadtreeMesh& sourceMesh,const QuadtreeMesh& targetMesh) const
{ const auto& d=checked();d.sourceMesh.validate(sourceMesh);d.targetMesh.validate(targetMesh); }

HermiteGalsRefinementInventory refineHermiteGalsInventory(const QuadtreeMesh& sourceMesh,
    const HermiteGalsState& source,const QuadtreeMesh& targetMesh,const HermiteGalsPhaseVolumeConfig& config)
{
    HermiteGalsIdentity identity;const char* stage="source";int parent=-1,child=-1;
    try {
        identity=source.identity();source.validateMesh(sourceMesh);
        const auto& primary=source.primary();const auto& oldScalar=source.geometry().snapshot();
        const auto oldCount=static_cast<std::size_t>(sourceMesh.cellIndexNum());
        require(primary.geometryCells.size()==oldCount && primary.phaseVolumes.size()==oldCount,
            "refinement requires full geometry and inventory coverage");
        for(std::size_t p=0;p<oldCount;++p)require(primary.geometryCells[p]==static_cast<int>(p)
            && primary.phaseVolumes[p].cell==static_cast<int>(p),"refinement requires complete native cell ordering");
        const auto oldMeasure=integrateHermiteGalsPhaseVolumes(sourceMesh,source,config);
        const bool radial=primary.measure==HermiteGalsMeasure::AxisymmetricVolume;
        constexpr auto budget=Detail::hermiteInventoryRoundoffBudget;
        auto data=std::make_shared<HermiteGalsRefinementInventory::Data>(sourceMesh,source,targetMesh);
        auto& report=data->report;Sum oldInventory,oldCapacity;
        for(parent=0;parent<sourceMesh.cellIndexNum();++parent){
            const auto& cell=oldMeasure.cells().at(parent);const auto m=primary.phaseVolumes[parent].liquid;
            require(m>=0 && m<=cell.capacity,"source inventory outside strict capacity bounds");
            const auto residual=std::abs(m-canonical(cell))/cell.capacity;
            require(std::isfinite(residual) && residual<=budget,"source inventory exceeds original geometry roundoff budget");
            report.maxRelativeSourceGeometryDefect=std::max(report.maxRelativeSourceGeometryDefect,residual);
            oldInventory.add(m);oldCapacity.add(cell.capacity);
        }
        parent=-1;stage="exact geometry";
        data->scalar=std::make_shared<const LevelSetHermiteGeometry::Snapshot>(
            LevelSetHermiteGeometry::Snapshot::refineNativeQuartic(targetMesh,oldScalar));
        FaceSegmentGraph graph;graph.build(targetMesh);
        const auto domain=NativeLeafDomain::build(targetMesh,graph,radial);
        stage="target measurement";
        const auto measured=InterfaceCellMeasure::discoverVolumes(targetMesh,domain,data->scalar,config.boundaryContacts,config.symmetryPlanes);
        const auto& raw=measured.measure;Sum totalCapacity;
        std::vector<std::vector<int>> children(oldCount);
        for(child=0;child<targetMesh.cellIndexNum();++child){
            const int row=domain.cellRowByIndex().at(child);require(row>=0,"missing target volume owner");
            const auto capacity=domain.cellVolumes().at(row);const auto liquid=raw.negativeVolumeByCell.at(child);
            const int phase=measured.cellPhase.at(child);
            require(std::isfinite(capacity) && capacity>0 && std::isfinite(liquid),"invalid target capacity/volume");
            require(std::abs(raw.fixedBoundary.capacityByCell.at(child)-capacity)/capacity<=config.relativeLedgerTolerance,
                "target boundary capacity defect");
            require(phase>=-1 && phase<=1,"invalid target phase certificate");
            if(!phase)require(liquid>=0 && liquid<=capacity,"target cut geometry outside strict capacity bounds");
            else require(std::abs(liquid-(phase<0?capacity:0))/capacity<=config.relativeLedgerTolerance,
                "target pure certificate/raw integral disagreement");
            data->cells.push_back({child,phase,capacity,liquid});
            data->volumes.push_back({child,static_cast<double>(canonical(data->cells.back()))});totalCapacity.add(capacity);
            const auto box=data->scalar->nativeCellBounds(child);parent=-1;
            for(int p:oldScalar.nativeCellCandidates(box[0],box[1])){const auto b=oldScalar.nativeCellBounds(p);
                if(std::max(box[0].x,b[0].x)>=std::min(box[1].x,b[1].x) || std::max(box[0].y,b[0].y)>=std::min(box[1].y,b[1].y))continue;
                require(parent<0 && box[0].x>=b[0].x && box[0].y>=b[0].y && box[1].x<=b[1].x && box[1].y<=b[1].y,
                    "target leaf has no unique containing parent");parent=p;}
            require(parent>=0,"missing refinement parent");data->parents.push_back(parent);children[parent].push_back(child);
        }
        parent=child=-1;
        require(totalCapacity.value>0 && std::isfinite(totalCapacity.value),"invalid total target capacity");
        require(std::isfinite(raw.fixedBoundary.assemblyClosure)
            && std::abs(raw.fixedBoundary.assemblyClosure)/totalCapacity.value<=config.relativeLedgerTolerance,"target boundary ledger does not close");
        require(std::isfinite(raw.momentErrorEstimate) && raw.momentErrorEstimate>=0
            && raw.momentErrorEstimate/totalCapacity.value<=config.relativeMomentEstimateTolerance,"target moment estimate exceeds original tolerance");
        stage="parent allocation";Sum transferred;
        for(parent=0;parent<sourceMesh.cellIndexNum();++parent){
            child=-1;
            const auto& group=children[parent];require(!group.empty(),"source parent is not covered");
            const auto capacity=oldMeasure.cells()[parent].capacity;const auto target=primary.phaseVolumes[parent].liquid;
            Sum geometry,base;for(int c:group){geometry.add(canonical(data->cells[c]));base.add(data->volumes[c].liquid);}
            const auto delta=static_cast<long double>(target)-base.value;
            report.maxRelativeParentReferenceDifference=std::max(report.maxRelativeParentReferenceDifference,std::abs(target-geometry.value)/capacity);
            if(group.size()==1){
                child=group.front();const auto a=oldScalar.nativeCellBounds(parent),b=data->scalar->nativeCellBounds(child);
                require(a[0].x==b[0].x && a[0].y==b[0].y && a[1].x==b[1].x && a[1].y==b[1].y,"single child does not cover parent");
                data->volumes[child].liquid=target;++report.unchangedParents;
            }else{
                ++report.refinedParents;
                if(delta!=0){
                    std::vector<long double> available;Sum total;
                    for(int c:group){child=c;const auto& cell=data->cells[c];const auto v=canonical(cell);const long double m=data->volumes[c].liquid;
                        // Admissible increments, not clipping an inventory.
                        const auto low=std::max(0.L,v-budget*cell.capacity),high=std::min(static_cast<long double>(cell.capacity),v+budget*cell.capacity);
                        const auto amount=delta>0?high-m:m-low;
                        require(std::isfinite(amount) && amount>=0,"rounded geometry has no admissible allocation interval");
                        available.push_back(amount);total.add(amount);}
                    child=-1; // Feasibility belongs to this parent, not the last visited child.
                    require(total.value>0 && std::abs(delta)<=total.value,"parent inventory has no feasible child allocation within original bounds");
                    const auto fraction=std::abs(delta)/total.value;
                    for(std::size_t k=0;k<group.size();++k){child=group[k];
                        const long double value=data->volumes[child].liquid+(delta>0?1:-1)*fraction*available[k];
                        data->volumes[child].liquid=static_cast<double>(value);}
                }
            }
            Sum sum;
            for(int c:group){child=c;const auto m=data->volumes[c].liquid;const auto& cell=data->cells[c];
                require(std::isfinite(m) && m>=0 && m<=cell.capacity,"allocated inventory violates strict capacity bounds");
                const auto residual=std::abs(m-canonical(cell))/cell.capacity;
                require(std::isfinite(residual) && residual<=budget,"allocated inventory exceeds original geometry roundoff budget");
                report.maxRelativeChildGeometryDefect=std::max(report.maxRelativeChildGeometryDefect,residual);sum.add(m);}
            const auto residual=sum.value-target;
            require(std::isfinite(residual) && std::abs(residual)<=budget*capacity,"parent inventory does not close within original roundoff budget");
            data->balances.push_back({parent,group.size(),capacity,target,geometry.value,sum.value,residual});
            report.maxRelativeParentResidual=std::max(report.maxRelativeParentResidual,std::abs(residual)/capacity);transferred.add(sum.value);
        }
        parent=child=-1;
        require(std::isfinite(transferred.value) && std::abs(transferred.value-oldInventory.value)<=budget*oldCapacity.value,
            "global transferred inventory does not close");
        report.originalInventory=oldInventory.value;report.transferredInventory=transferred.value;
        report.childCells=data->volumes.size();return HermiteGalsRefinementInventory(std::move(data));
    }catch(const std::exception& e){throw std::invalid_argument("[HermiteGalsRefinementInventory source="+std::to_string(identity.value())
        +" stage="+stage+" parent="+std::to_string(parent)+" child="+std::to_string(child)+"] "+e.what());}
}
}
