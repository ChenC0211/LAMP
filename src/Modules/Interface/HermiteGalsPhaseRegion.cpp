#include "Modules/Interface/HermiteGalsPhaseRegion.h"
#include "Modules/Interface/Detail/HermiteGalsRegionMeasure.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>

namespace Quantum::Interface {
namespace {
using I=LevelSetMappedInterface;
using Snapshot=LevelSetHermiteGeometry::Snapshot;
constexpr long double pi=3.1415926535897932384626433832795L;
void require(bool ok,const std::string& message){if(!ok)throw std::invalid_argument("[HermiteGalsPhaseRegion] "+message);}
struct Sum {long double value=0,correction=0;void add(long double x){const auto a=x-correction,b=value+a;correction=(b-value)-a;value=b;}};
long double capacity(Vector2D a,Vector2D b,bool radial){return (static_cast<long double>(b.x)-a.x)*(static_cast<long double>(b.y)-a.y)*
    (radial?pi*(static_cast<long double>(b.x)+a.x):1);}
using Box=Detail::RegionMeasure::Box;
// Exact ordering/coverage of represented rectangle endpoints. No area-sum
// tolerance or rounded interior sample can conceal a small inactive gap.
void coverage(Vector2D lo,Vector2D hi,const std::vector<Box>& boxes){
    std::vector<double> xs{lo.x,hi.x};for(const auto& b:boxes){xs.push_back(b.lower.x);xs.push_back(b.upper.x);}
    std::sort(xs.begin(),xs.end());xs.erase(std::unique(xs.begin(),xs.end()),xs.end());
    for(std::size_t k=1;k<xs.size();++k){std::vector<std::pair<double,double>> ys;
        for(const auto& b:boxes)if(b.lower.x<=xs[k-1] && b.upper.x>=xs[k])ys.emplace_back(b.lower.y,b.upper.y);
        std::sort(ys.begin(),ys.end());double reach=lo.y;
        for(const auto& y:ys){require(y.first==reach,"rectangle has a gap or overlapping active leaves");reach=y.second;}
        require(reach==hi.y,"rectangle leaves active scalar support");}
}
using CellMeasure=Detail::RegionMeasure::CellMeasure;
using Detail::RegionMeasure::integrateCell;
}
struct HermiteGalsPhaseRegionSource::Data {
    explicit Data(const HermiteGalsState& s,const HermiteGalsPhaseRegionConfig& c):source(s),config(c){}
    HermiteGalsState source;
    HermiteGalsPhaseRegionConfig config;
    std::shared_ptr<const Snapshot> scalar;
    LevelSetInterfacePanels::NativeComponents topology;
    std::vector<std::vector<I::Graph>> graphs;
    std::size_t nativeSideGraphs=0;
};
HermiteGalsPhaseRegionSource::HermiteGalsPhaseRegionSource(const QuadtreeMesh& mesh,const HermiteGalsState& source,
    const HermiteGalsPhaseRegionConfig& config){
    source.validateMesh(mesh);
    for(double v:{config.relativeLedgerTolerance,config.relativeMomentEstimateTolerance})
        require(std::isfinite(v) && v>0 && v<1,"invalid region relative tolerance");
    auto data=std::make_shared<Data>(source,config);auto view=std::make_shared<const HermiteGalsGeometryView>(source.geometry());
    data->scalar=std::shared_ptr<const Snapshot>(view,&view->snapshot());const auto& scalar=*data->scalar;
    require(scalar.usesNativeLeaves() && scalar.reconstruction()==LevelSetHermiteGeometry::Reconstruction::EnrichedQuartic,"region requires native Q4");
    const bool radial=source.primary().measure==HermiteGalsMeasure::AxisymmetricVolume;
    require(!scalar.radialAxisSymmetry() || radial,"radial scalar with Cartesian region measure");
    const bool nativeTrace=source.primary().transport==HermiteGalsTransport::FirstJet;
    data->topology=nativeTrace?LevelSetInterfacePanels::discoverNativeTraceComponents(mesh,data->scalar,radial,config.boundaryContacts,config.symmetryPlanes)
        :LevelSetInterfacePanels::discoverNativeComponents(mesh,data->scalar,radial,config.boundaryContacts,config.symmetryPlanes);
    data->graphs.resize(mesh.cellIndexNum());
    for(const auto& chain:data->topology.curves){require(chain.usesSource(scalar) && chain.axisymmetric()==radial,"foreign region component");
        for(const auto& graph:chain.graphs()){require(graph.nativeCellIndex>=0 && graph.nativeCellIndex<mesh.cellIndexNum(),"foreign native region graph");
            if(graph.nativeSide) {
                require(nativeTrace && graph.nativeSide->usesSource(scalar),"foreign/legacy region side support");
                for(const auto* side:{&graph.nativeSide->minus(),&graph.nativeSide->plus()})if(side->available())
                    require((side->phase==-1 || side->phase==1) && data->topology.cellPhase.at(side->nativeCellIndex)==side->phase,"region side lacks complete pure interior certificate");
                ++data->nativeSideGraphs;continue; // Boundary geometry adds no interior box volume.
            }
            data->graphs[graph.nativeCellIndex].push_back(graph);}}
    data_=std::move(data);
}
const HermiteGalsPhaseRegionSource::Data& HermiteGalsPhaseRegionSource::checked()const{require(bool(data_),"moved-from region source");return *data_;}
const HermiteGalsState& HermiteGalsPhaseRegionSource::source()const{return checked().source;}
std::size_t HermiteGalsPhaseRegionSource::components()const{return checked().topology.curves.size();}
std::size_t HermiteGalsPhaseRegionSource::nativeSideGraphs()const{return checked().nativeSideGraphs;}
void HermiteGalsPhaseRegionSource::validateMesh(const QuadtreeMesh& mesh)const{checked().source.validateMesh(mesh);}
HermiteGalsPhaseRegionMeasure HermiteGalsPhaseRegionSource::measure(const QuadtreeMesh& mesh,Vector2D lo,Vector2D hi)const{
    validateMesh(mesh);return measureValidatedMesh(mesh,lo,hi);
}
HermiteGalsPhaseRegionMeasure HermiteGalsPhaseRegionSource::measureValidatedMesh(const QuadtreeMesh& mesh,Vector2D lo,Vector2D hi)const{
    const auto& data=checked();const auto& scalar=*data.scalar;
    const bool radial=data.source.primary().measure==HermiteGalsMeasure::AxisymmetricVolume;
    require(std::isfinite(lo.x) && std::isfinite(lo.y) && std::isfinite(hi.x) && std::isfinite(hi.y) && lo.x<hi.x && lo.y<hi.y && (!radial || lo.x>=0),"invalid rectangle bounds");
    std::vector<Box> boxes;
    for(int cell:scalar.nativeCellCandidates(lo,hi)){const auto b=scalar.nativeCellBounds(cell);
        const Vector2D a{std::max(lo.x,b[0].x),std::max(lo.y,b[0].y)},z{std::min(hi.x,b[1].x),std::min(hi.y,b[1].y)};
        if(a.x>=z.x || a.y>=z.y)continue;require(scalar.nativeCellActive(cell),"rectangle intersects inactive scalar support");
        boxes.push_back({a,z,cell});}
    coverage(lo,hi,boxes);
    HermiteGalsPhaseRegionMeasure result;result.sourceIdentity=data.source.identity();result.lower=lo;result.upper=hi;
    result.capacity=capacity(lo,hi,radial);require(std::isfinite(result.capacity) && result.capacity>0,"unrepresentable rectangle capacity");
    Sum axialLiquid,radialLiquid,axialGas,radialGas,negative,positive,error;bool allNegative=true,allPositive=true;
    for(const auto& box:boxes){++result.intersectedLeaves;const auto cap=capacity(box.lower,box.upper,radial);
        require(std::isfinite(cap) && cap>0,"unrepresentable leaf intersection capacity");
        int phase=data.topology.cellPhase.at(box.cell);if(phase==0)phase=scalar.nativeQuarticBoxSign(box.cell,box.lower,box.upper);
        CellMeasure m;
        if(phase!=0){require(phase==-1 || phase==1,"invalid phase certificate");
            m.axialLiquid=m.radialLiquid=phase<0?cap:0;m.axialGas=m.radialGas=phase>0?cap:0;}
        else{require(!data.graphs[box.cell].empty(),"uncertified graph-free leaf");++result.cutLeaves;
            m=integrateCell(data.scalar,box,data.graphs[box.cell],radial);}
        allNegative=allNegative && phase<0;allPositive=allPositive && phase>0;
        for(long double value:{m.axialLiquid,m.radialLiquid,m.axialGas,m.radialGas,m.error,
            m.fixedAxialAssemblyDefect,m.fixedRadialAssemblyDefect})require(std::isfinite(value),"nonfinite region moment");
        // The new inside-trace clip ledger includes the actual coincident
        // interface term. It is independent of the legacy shared-face ledger
        // and imposes no material-flux or per-cell inventory constraint.
        if(m.assemblyCapacity>0)require(std::max(std::abs(m.fixedAxialAssemblyDefect),std::abs(m.fixedRadialAssemblyDefect))
            <=data.config.relativeLedgerTolerance*m.assemblyCapacity,"region inside-trace fixed/coincident ledger defect");
        const auto defect=std::max({std::abs(m.axialLiquid-m.radialLiquid),std::abs(m.axialGas-m.radialGas),
            std::abs(m.axialLiquid+m.axialGas-cap),std::abs(m.radialLiquid+m.radialGas-cap)})/cap;
        result.maxRelativeLedgerDefect=std::max(result.maxRelativeLedgerDefect,defect);
        require(defect<=data.config.relativeLedgerTolerance,"region Green/capacity ledger defect");
        const bool axial=box.upper.y-box.lower.y<=box.upper.x-box.lower.x;
        const auto liquid=axial?m.axialLiquid:m.radialLiquid,gas=axial?m.axialGas:m.radialGas;
        require(liquid>=0 && liquid<=cap && gas>=0 && gas<=cap,"region phase volume outside strict capacity; no clipping");
        const auto phases=Detail::RegionMeasure::complement(cap,liquid,gas);
        result.leaves.push_back({box.cell,box.lower,box.upper,cap,phases[0],phases[1]});
        negative.add(liquid);positive.add(gas);axialLiquid.add(m.axialLiquid);radialLiquid.add(m.radialLiquid);
        axialGas.add(m.axialGas);radialGas.add(m.radialGas);error.add(m.error);result.integratedArcs+=m.arcs;result.rootRefinements+=m.refinements;
    }
    result.rawAxialLiquid=axialLiquid.value;result.rawRadialLiquid=radialLiquid.value;result.rawAxialGas=axialGas.value;result.rawRadialGas=radialGas.value;
    result.momentErrorEstimate=error.value;
    require(error.value>=0 && error.value<=data.config.relativeMomentEstimateTolerance*result.capacity,"region moment error estimate exceeds tolerance");
    require(std::abs(negative.value+positive.value-result.capacity)<=data.config.relativeLedgerTolerance*result.capacity,"region assembled capacity defect");
    if(allNegative){result.liquid=result.capacity;result.gas=0;}
    else if(allPositive){result.liquid=0;result.gas=result.capacity;}
    else{const auto phases=Detail::RegionMeasure::complement(result.capacity,negative.value,positive.value);result.liquid=phases[0];result.gas=phases[1];}
    require(std::isfinite(result.liquid) && std::isfinite(result.gas) && result.liquid>=0 && result.gas>=0 && result.liquid<=result.capacity && result.gas<=result.capacity,"assembled region outside strict capacity");
    return result;
}
}
