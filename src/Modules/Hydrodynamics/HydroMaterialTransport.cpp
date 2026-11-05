#include "Modules/Hydrodynamics/HydroMaterialTransport.h"
#include "Core/Geometry/Detail/ExactExpansion.h"
#include "Core/Mesh/QuadtreeMesh.h"

#include <algorithm>
#include <cmath>
#include <chrono>
#include <exception>
#include <limits>
#include <map>
#include <optional>
#include <functional>
#include <iomanip>
#include <sstream>
#include <future>
#include <stdexcept>
#include <string>
#include <utility>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace LAMP::Hydrodynamics
{
namespace
{
void require(bool ok,const char* why)
{if(!ok)throw std::invalid_argument(std::string("[HydroMaterialTransport] ")+why);}
void require(bool ok,const std::string& why)
{if(!ok)throw std::invalid_argument("[HydroMaterialTransport] "+why);}
long double finite(long double value,const char* why,const char* suffix="")
{
    if(!std::isfinite(value))require(false,std::string("nonfinite ")+why+suffix);
    return value;
}
long double finite(long double value,const std::string& why)
{require(std::isfinite(value),"nonfinite "+why);return value;}
long double densitySum(HydroPhaseMeasure measure,double liquid,double gas,const char* why)
{
    const auto l=finite(liquid*measure.liquid,why," liquid component");
    const auto g=finite(gas*measure.gas,why," gas component");
    if(!(measure.liquid==0 || l!=0))require(false,std::string("unrepresentable ")+why+" liquid component");
    if(!(measure.gas==0 || g!=0))require(false,std::string("unrepresentable ")+why+" gas component");
    return finite(l+g,why);
}
// Retain dynamic-detail callers without allocating strings for literal paths.
[[maybe_unused]] long double densitySum(HydroPhaseMeasure measure,double liquid,double gas,const std::string& why)
{
    const auto l=finite(liquid*measure.liquid,why+" liquid component");
    const auto g=finite(gas*measure.gas,why+" gas component");
    require(measure.liquid==0 || l!=0,"unrepresentable "+why+" liquid component");
    require(measure.gas==0 || g!=0,"unrepresentable "+why+" gas component");
    return finite(l+g,why);
}
// Private to this call: another owner has just checked the LIVE mesh.
// Reuse only the exact immutable binding plus successful owner/full-clock
// predicates. A foreign binding or failed predicate takes the original public
// validator, preserving its diagnostics. Never memoize across public calls.
void validateVofPhaseAfterMeshCheck(const QuadtreeMesh& mesh,const HydroPhaseGeometry& phase,
    Interface::GeometricVofIdentity owner,LAMP::TimeLayer time,const NativeLeafDomain& checkedDomain)
{
    const auto& source=phase.vofPublication();const auto& domain=source.domain();
    time.time();const auto actual=source.timeLayer();
    if(domain.meshIdentity()==checkedDomain.meshIdentity()&&owner==source.identity()&&owner.value()!=0
        &&time.step==actual.step&&time.begin==actual.begin&&time.end==actual.end&&time.position==actual.position)return;
    phase.validateVofFor(mesh,owner,time);
}
void addFirstJetTimings(Interface::HermiteGalsFirstJetStepTimings& sum,
    const Interface::HermiteGalsFirstJetStepTimings& value)
{
    sum.input+=value.input;sum.boundary+=value.boundary;sum.vertices+=value.vertices;sum.edges+=value.edges;
    sum.centers+=value.centers;sum.mixed+=value.mixed;sum.publication+=value.publication;sum.total+=value.total;
}
void addFullTimings(HydroSpaceTimeTimings& sum,const HydroSpaceTimeTimings& value)
{
    addFirstJetTimings(sum.fullFirstJet,value.fullFirstJet);
    sum.fullGeometrySeconds+=value.fullGeometrySeconds;sum.fullRateSeconds+=value.fullRateSeconds;
    sum.fullSamples+=value.fullSamples;sum.fullCharacteristicSeconds+=value.fullCharacteristicSeconds;
}
void addInternalTimings(HydroSpaceTimeTimings& sum,const HydroSpaceTimeTimings& value)
{
    addFirstJetTimings(sum.internalFirstJet,value.internalFirstJet);
    sum.internalRateSeconds+=value.internalRateSeconds;sum.internalSamples+=value.internalSamples;
    sum.internalCharacteristicSeconds+=value.internalCharacteristicSeconds;
}
struct SamplingProfile : HydroSpaceTimeTimings
{
    HydroSpaceTimeSamplingCost samplingCost;
};
void addSamplingCost(HydroSpaceTimeSamplingCost& sum,const HydroSpaceTimeSamplingCost& value)
{
    sum.localSamples+=value.localSamples;sum.fullSamples+=value.fullSamples;
    sum.localPreparedNativeCells+=value.localPreparedNativeCells;sum.fullPreparedNativeCells+=value.fullPreparedNativeCells;
    sum.localAtomRateEvaluations+=value.localAtomRateEvaluations;sum.fullAtomRateEvaluations+=value.fullAtomRateEvaluations;
    sum.localScalarPointQueries+=value.localScalarPointQueries;
    sum.localPointJobs+=value.localPointJobs;sum.localPointTraceQueries+=value.localPointTraceQueries;
    sum.pointThreads=std::max(sum.pointThreads,value.pointThreads);
    sum.localNumericalPointQueries+=value.localNumericalPointQueries;
    sum.localExactValueFallbacks+=value.localExactValueFallbacks;
    sum.localExactGradientFallbacks+=value.localExactGradientFallbacks;
    sum.maxAcceptedTransportValueError=std::max(sum.maxAcceptedTransportValueError,value.maxAcceptedTransportValueError);
    sum.maxAcceptedTransportGradientError=std::max(sum.maxAcceptedTransportGradientError,value.maxAcceptedTransportGradientError);
}
long double atomSum(const std::vector<long double>& values,const std::vector<int>& indices)
{
    long double sum=0,correction=0;
    for(int i:indices){const auto y=values[i]-correction,t=sum+y;correction=(t-sum)-y;sum=t;}
    return finite(sum,"spatial material atom sum");
}
}
struct HydroMaterialFlux::Data
{
    Data(HydroPhaseRate p,double l,double g):phase(std::move(p)),liquid(l),gas(g){}
    HydroPhaseRate phase;
    double liquid,gas;
    std::vector<long double> atoms,faces,dual;
};
HydroMaterialFlux::HydroMaterialFlux(std::shared_ptr<const Data> d):data_(std::move(d)){}
HydroMaterialFlux HydroMaterialFlux::prepare(const QuadtreeMesh& mesh,HydroPhaseRate phase,double liquid,double gas)
{
    require(phase.isSpatial(),"requires spatial phase owner for material flux");
    phase.geometry().validateFor(mesh,phase.geometry().geometryIdentity(),phase.geometry().timeLayer());
    phase.spatialAdvector().validateFor(mesh,phase.momentumGeometry(),phase.geometry().geometryIdentity(),phase.geometry().timeLayer());
    require(std::isfinite(liquid) && liquid>0 && std::isfinite(gas) && gas>0,"nonpositive/nonfinite density");
    auto data=std::make_shared<Data>(std::move(phase),liquid,gas);
    for(const auto& q:data->phase.atomRates())data->atoms.push_back(densitySum(q,liquid,gas,"atom material rate"));
    const auto& reconstruction=data->phase.spatialAdvector().reconstruction();
    for(const auto& indices:reconstruction.primalAtoms())data->faces.push_back(atomSum(data->atoms,indices));
    for(const auto& indices:reconstruction.dualAtoms())data->dual.push_back(atomSum(data->atoms,indices));
    return HydroMaterialFlux(std::move(data));
}
const HydroPhaseRate& HydroMaterialFlux::phaseRate() const noexcept{return data_->phase;}
const HydroMomentumGeometry& HydroMaterialFlux::momentumGeometry() const{return phaseRate().momentumGeometry();}
const HydroDofLayout& HydroMaterialFlux::layout() const noexcept{return phaseRate().geometry().layout();}
HydroAdvectorIdentity HydroMaterialFlux::advectorIdentity() const{return phaseRate().advectorIdentity();}
Interface::HermiteGalsIdentity HydroMaterialFlux::geometryIdentity() const{return phaseRate().geometry().geometryIdentity();}
LAMP::TimeLayer HydroMaterialFlux::timeLayer() const{return phaseRate().geometry().timeLayer();}
double HydroMaterialFlux::liquidDensity() const noexcept{return data_->liquid;}
double HydroMaterialFlux::gasDensity() const noexcept{return data_->gas;}
const std::vector<long double>& HydroMaterialFlux::atomRates() const noexcept{return data_->atoms;}
const std::vector<long double>& HydroMaterialFlux::faceRates() const noexcept{return data_->faces;}
const std::vector<long double>& HydroMaterialFlux::dualFaceRates() const noexcept{return data_->dual;}
Vector2D HydroMaterialFlux::physicalFlux(int pressure,Vector2D point,bool liquid) const
{
    const auto velocity=phaseRate().spatialAdvector().cellPhysicalVelocityJet(pressure,point,liquid).jet.value;
    const auto flux=velocity*(liquid?liquidDensity():gasDensity());
    require(std::isfinite(flux.x)&&std::isfinite(flux.y),"nonfinite physical material flux");
    return flux;
}
void HydroMaterialFlux::validateFor(const QuadtreeMesh& mesh,const HydroMomentumGeometry& expected,
    Interface::HermiteGalsIdentity geometry,HydroAdvectorIdentity advector,LAMP::TimeLayer clock) const
{
    phaseRate().geometry().validateFor(mesh,geometry,clock);
    phaseRate().spatialAdvector().validateFor(mesh,expected,geometry,clock);
    require(advector==advectorIdentity(),"different actual advector identity");
}
struct HydroMaterialTransport::Data
{
    Data(HydroPhaseGeometry a,HydroPhaseRate m,HydroPhaseGeometry b,double l,double g,double dt,const HydroMaterialFlux* flux)
        :begin(std::move(a)),midpoint(std::move(m)),end(std::move(b)),liquid(l),gas(g),deltaTime(dt),
         spatial(*flux){}
    Data(HydroPhaseGeometry a,HydroPhaseGeometry m,HydroPhaseGeometry b,
        Interface::GeometricVofTransport t,HydroAdvectorReconstruction r,double l,double g,double dt)
        :begin(std::move(a)),end(std::move(b)),liquid(l),gas(g),deltaTime(dt),
         volumeMiddle(std::move(m)),vof(std::move(t)),volumeReconstruction(std::move(r)){}
    HydroPhaseGeometry begin;
    std::optional<HydroPhaseRate> midpoint;
    HydroPhaseGeometry end;
    double liquid,gas,deltaTime;
    std::optional<HydroMaterialFlux> spatial;
    std::optional<HydroPhaseGeometry> volumeMiddle;
    std::optional<Interface::GeometricVofTransport> vof;
    std::optional<HydroAdvectorReconstruction> volumeReconstruction;
    std::vector<long double> integrals,atoms,atomIntegrals,faces,dual,dualIntegrals;
    std::vector<HydroPhaseMeasure> phaseIntegrals;
    std::optional<HydroCharacteristicSource> motion;
    HydroSpaceTimeReport spaceTime;
    std::vector<TimePanel> acceptedTimePanels;
    std::vector<HydroMaterialCell> cells;
    HydroMaterialReport report;
};
HydroMaterialTransport::HydroMaterialTransport(std::shared_ptr<const Data> d):data_(std::move(d)){}
HydroMaterialTransport HydroMaterialTransport::prepare(const QuadtreeMesh& mesh,HydroPhaseGeometry begin,
    HydroPhaseRate midpoint,HydroPhaseGeometry end,double liquid,double gas)
{
    return prepare(mesh,std::move(begin),HydroMaterialFlux::prepare(mesh,midpoint,liquid,gas),std::move(end));
}
HydroMaterialTransport HydroMaterialTransport::prepare(const QuadtreeMesh& mesh,HydroPhaseGeometry begin,
    HydroMaterialFlux midpoint,HydroPhaseGeometry end)
{
    return assemble(mesh,std::move(begin),midpoint.phaseRate(),std::move(end),midpoint.liquidDensity(),midpoint.gasDensity(),&midpoint);
}
HydroMaterialTransport HydroMaterialTransport::prepareGeometricVof(const QuadtreeMesh& mesh,
    HydroMomentumGeometry geometry,HydroPhaseGeometry begin,HydroPhaseGeometry midpoint,
    HydroPhaseGeometry end,Interface::GeometricVofTransport accepted)
{return prepareGeometricVofBound(mesh,std::move(geometry),std::move(begin),std::move(midpoint),std::move(end),std::move(accepted),nullptr);}
HydroMaterialTransport HydroMaterialTransport::prepareGeometricVof(const QuadtreeMesh& mesh,
    HydroMomentumGeometry geometry,HydroPhaseGeometry begin,HydroPhaseGeometry midpoint,
    HydroPhaseGeometry end,Interface::GeometricVofTransport accepted,const HydroAdvectorBoundary& cacheOwner)
{return prepareGeometricVofBound(mesh,std::move(geometry),std::move(begin),std::move(midpoint),std::move(end),std::move(accepted),&cacheOwner);}
HydroMaterialTransport HydroMaterialTransport::prepareGeometricVofBound(const QuadtreeMesh& mesh,
    HydroMomentumGeometry geometry,HydroPhaseGeometry begin,HydroPhaseGeometry midpoint,
    HydroPhaseGeometry end,Interface::GeometricVofTransport accepted,const HydroAdvectorBoundary* cacheOwner)
{
    using Model=HydroPhaseGeometry::Model;
    accepted.validateFor(mesh,accepted.source(),accepted.advector(),accepted.liquidDensity(),accepted.gasDensity());
    const double liquid=static_cast<double>(accepted.liquidDensity()),gas=static_cast<double>(accepted.gasDensity());
    require(std::isfinite(liquid)&&liquid>0&&std::isfinite(gas)&&gas>0
        &&liquid==accepted.liquidDensity()&&gas==accepted.gasDensity(),"VOF densities not exactly representable in Hydro");
    auto midClock=accepted.endTime();midClock.position=LAMP::TimeLayer::Position::Midpoint;
    // The accepted operator has just checked the live mesh. Reuse only that
    // same immutable native binding within this callback-free invocation.
    const auto& checkedDomain=accepted.source().domain();
    validateVofPhaseAfterMeshCheck(mesh,begin,accepted.sourceIdentity(),accepted.sourceTime(),checkedDomain);
    validateVofPhaseAfterMeshCheck(mesh,midpoint,accepted.midpoint().identity(),midClock,checkedDomain);
    validateVofPhaseAfterMeshCheck(mesh,end,accepted.endpoint().identity(),accepted.endTime(),checkedDomain);
    if(geometry.layout().domain().meshIdentity()!=checkedDomain.meshIdentity())geometry.validateMesh(mesh);
    const auto& layout=geometry.layout();
    require(layout.domain().meshIdentity()==accepted.source().domain().meshIdentity(),"different original VOF native domain");
    for(const auto* phase:std::array<const HydroPhaseGeometry*,3>{{&begin,&midpoint,&end}}){
        require(phase->model()==Model::GeometricVof&&phase->layout().sameIdentity(layout),"different VOF native layout");
        require(&phase->metrics().divergence().matrix()==&geometry.metrics().divergence().matrix(),"different retained VOF D_I");
    }
    const double dt=midClock.end-midClock.begin;
    require(std::isfinite(dt)&&dt>0&&midClock.time()>midClock.begin&&midClock.time()<midClock.end,"unrepresentable VOF full interval");
    auto reconstruction=cacheOwner?HydroAdvectorReconstruction::prepareUnboundCached(mesh,std::move(geometry),*cacheOwner):
        HydroAdvectorReconstruction::prepare(mesh,std::move(geometry));
    std::vector<Interface::GeometricVofNativeSubface> nativeQueries;
    std::vector<Interface::GeometricVofInteriorSubface> interiorQueries;
    std::vector<std::pair<bool,std::size_t>> origins;
    for(const auto& atom:reconstruction.atoms()){
        require((atom.primalVelocity>=0)!=(atom.ownerPressure>=0),"canonical VOF atom has no unique origin");
        const Geometry::PlanarCellCut::Point first{atom.begin.x,atom.begin.y},last{atom.end.x,atom.end.y};
        if(atom.primalVelocity>=0){
            origins.emplace_back(true,nativeQueries.size());
            nativeQueries.push_back({layout.nativeFace(atom.primalVelocity),first,last});
        }else{
            origins.emplace_back(false,interiorQueries.size());
            interiorQueries.push_back({layout.nativeCell(atom.ownerPressure),first,last});
        }
    }
    const bool fullOnly=accepted.source().geometryKind()==Interface::GeometricVofGeometryKind::LocalPlic;
    std::vector<Interface::GeometricVofNativeFullFlux> nativeFull;
    std::vector<Interface::GeometricVofInteriorFullFlux> interiorFull;
    std::vector<Interface::GeometricVofTransportFace> native;
    std::vector<Interface::GeometricVofInteriorFlux> interior;
    if(fullOnly){nativeFull=accepted.restrictNativeFacesFull(mesh,nativeQueries);
        interiorFull=accepted.measureInteriorFacesFull(mesh,interiorQueries);}
    else {native=accepted.restrictNativeFaces(mesh,nativeQueries);
        interior=accepted.measureInteriorFaces(mesh,interiorQueries);}
    auto data=std::make_shared<Data>(std::move(begin),std::move(midpoint),std::move(end),
        std::move(accepted),std::move(reconstruction),liquid,gas,dt);
    std::vector<Interface::GeometricVofIntegratedFlux> atomFlux;atomFlux.reserve(origins.size());
    for(const auto& origin:origins){
        if(fullOnly)atomFlux.push_back(origin.first?nativeFull[origin.second].full:interiorFull[origin.second].full);
        else atomFlux.push_back(origin.first?native[origin.second].full:interior[origin.second].full);
    }
    return assembleIntegrated(std::move(data),std::move(atomFlux));
}
HydroMaterialTransport HydroMaterialTransport::assembleIntegrated(std::shared_ptr<Data> data,
    std::vector<Interface::GeometricVofIntegratedFlux> atomFlux)
{
    const auto& layout=data->volumeReconstruction->momentumGeometry().layout();
    const auto liquid=data->liquid,gas=data->gas,dt=data->deltaTime;
    for(const auto& f:atomFlux){
        data->phaseIntegrals.push_back({f.liquid,f.gas});
        data->atomIntegrals.push_back(finite(f.mass,"canonical atom material integral"));
        const auto rate=finite(f.mass/dt,"average atom material rate");
        require(f.mass==0||rate!=0,"unrepresentable atom material rate");data->atoms.push_back(rate);
    }
    const auto& r=*data->volumeReconstruction;
    const std::array<long double Interface::GeometricVofIntegratedFlux::*,4> branches{{
        &Interface::GeometricVofIntegratedFlux::positiveLiquid,&Interface::GeometricVofIntegratedFlux::negativeLiquid,
        &Interface::GeometricVofIntegratedFlux::positiveGas,&Interface::GeometricVofIntegratedFlux::negativeGas}};
    std::array<std::vector<long double>,4> grossAtoms;
    for(std::size_t k=0;k<branches.size();++k){grossAtoms[k].reserve(atomFlux.size());
        for(const auto& f:atomFlux)grossAtoms[k].push_back(f.*branches[k]);}
    const auto& acceptedFaces=data->vof->faces();
    for(int u=0;u<layout.velocityCount();++u){
        const auto& ids=r.primalAtoms()[u];
        // A native accepted payload remains the inventory authority. The one
        // atom table defines BOTH consumer views; never repair a last atom.
        const auto& parent=acceptedFaces[layout.nativeFace(u)].full;
        for(std::size_t k=0;k<branches.size();++k){
            const auto actual=atomSum(grossAtoms[k],ids),expected=parent.*branches[k];
            if(!(expected==0?actual==0:std::abs(actual-expected)/expected<=1e-12L)){
                const std::array<const char*,4> names{{"positive_liquid","negative_liquid","positive_gas","negative_gas"}};
                const auto clock=data->vof->endTime();std::size_t zeros=0;
                for(int id:ids)if(grossAtoms[k][id]==0)++zeros;
                std::ostringstream detail;detail<<std::scientific<<std::setprecision(std::numeric_limits<long double>::max_digits10)
                    <<"VOF native atom own-phase addback failed native_row="<<layout.nativeFace(u)<<" velocity_dof="<<u
                    <<" branch="<<names[k]<<" expected="<<expected<<" actual="<<actual<<" difference="<<(actual-expected)
                    <<" relative="<<(expected==0?std::numeric_limits<long double>::infinity():std::abs(actual-expected)/expected)
                    <<" expected_subnormal="<<(std::fpclassify(expected)==FP_SUBNORMAL)
                    <<" denorm_min="<<std::numeric_limits<long double>::denorm_min()
                    <<" atom_count="<<ids.size()<<" zero_atoms="<<zeros
                    <<" source_identity="<<data->vof->sourceIdentity().value()<<" step="<<clock.step
                    <<" begin="<<clock.begin<<" end="<<clock.end<<" position="<<static_cast<int>(clock.position)
                    <<" expected_hex="<<std::hexfloat<<expected<<" actual_hex="<<actual<<" atoms=";
                for(int id:ids)detail<<id<<':'<<grossAtoms[k][id]<<',';
                require(false,detail.str());
            }
        }
        const auto f=atomSum(data->atomIntegrals,ids);data->integrals.push_back(f);
        const auto rate=finite(f/dt,"VOF primal average material rate");require(f==0||rate!=0,"unrepresentable VOF primal rate");data->faces.push_back(rate);
    }
    for(const auto& ids:r.dualAtoms()){
        const auto f=atomSum(data->atomIntegrals,ids);data->dualIntegrals.push_back(f);
        const auto rate=finite(f/dt,"VOF dual average material rate");require(f==0||rate!=0,"unrepresentable VOF dual rate");data->dual.push_back(rate);
    }
    const bool localPlicInventory=data->vof&&data->begin.vofPublication().geometryKind()==Interface::GeometricVofGeometryKind::LocalPlic;
    data->report.usesPrimaryAlphaInventory=localPlicInventory;
    const auto primalMeasure=[&](const HydroPhaseGeometry& phase,int pressure){
        if(localPlicInventory){
            // The material inventory is independent of the measured PLIC cut.
            // Preserve both authoritative extensive phases, including minors.
            const auto& stock=phase.vofPublication().plicCells()[layout.nativeCell(pressure)].inventory;
            return HydroPhaseMeasure{stock[0],stock[1]};
        }
        return phase.cellVolumes()[pressure];
    };
    data->cells.resize(layout.pressureCount());
    for(int p=0;p<layout.pressureCount();++p){
        auto& cell=data->cells[p];
        const auto va=primalMeasure(data->begin,p),vm=primalMeasure(data->volumeMiddle.value(),p),vb=primalMeasure(data->end,p);
        cell.massBegin=densitySum(va,liquid,gas,"begin geometry mass");
        cell.massMidpoint=densitySum(vm,liquid,gas,"midpoint geometry mass");
        cell.massEnd=densitySum(vb,liquid,gas,"end geometry mass");
        require(cell.massBegin>0 && cell.massMidpoint>0 && cell.massEnd>0,"unrepresentable positive geometry mass");
        // Difference the retained phases first; LocalPlic uses its independent
        // stocks. Never subtract rounded masses or reconstruct a missing phase.
        cell.massChange=densitySum({vb.liquid-va.liquid,vb.gas-va.gas},liquid,gas,"geometry mass change");
        data->report.geometryMass[0]+=cell.massBegin;data->report.geometryMass[1]+=cell.massMidpoint;data->report.geometryMass[2]+=cell.massEnd;
        data->report.massChange+=cell.massChange;
    }
    std::vector<int> incidence(layout.velocityCount()),boundarySign(layout.velocityCount());
    long double balanceScale=0;
    const auto& di=data->volumeMiddle.value().metrics().divergence().matrix();
    for(Eigen::Index p=0;p<di.outerSize();++p){
        auto& cell=data->cells[static_cast<std::size_t>(p)];long double gross=std::abs(cell.massChange);
        for(HydroIntegralDivergence::Matrix::InnerIterator it(di,p);it;++it){
            const int f=static_cast<int>(it.col());
            const double area=data->volumeMiddle.value().metrics().faceMeasures()[f];
            require(area>0 && std::abs(it.value())==area,"different retained native face measure");
            const int sign=it.value()>0?1:-1;
            cell.outwardIntegral+=sign*data->integrals[f];gross+=std::abs(data->integrals[f]);
            ++incidence[f];boundarySign[f]+=sign;
        }
        finite(cell.outwardIntegral,"cell outward material integral");finite(gross,"cell defect scale");
        cell.integratedDefect=finite(cell.massChange+cell.outwardIntegral,"integrated material defect");
        cell.defectRate=finite(cell.integratedDefect/dt,"material continuity defect rate");
        require(gross>0 || cell.integratedDefect==0,"zero defect scale with nonzero defect");
        cell.relativeDefect=gross>0?finite(std::abs(cell.integratedDefect)/gross,"relative material defect"):0;
        data->report.integratedDefect+=cell.integratedDefect;
        data->report.l1IntegratedDefect+=std::abs(cell.integratedDefect);
        data->report.maxAbsDefectRate=std::max(data->report.maxAbsDefectRate,std::abs(cell.defectRate));
        data->report.maxRelativeCellDefect=std::max(data->report.maxRelativeCellDefect,cell.relativeDefect);
        balanceScale+=gross;
    }
    for(int f=0;f<layout.velocityCount();++f){
        if(incidence[f]==0){
            require(data->volumeMiddle.value().metrics().faceMeasures()[f]==0 && data->integrals[f]==0,
                "unrepresented positive-area face");
        }else require((incidence[f]==1 && std::abs(boundarySign[f])==1) ||
            (incidence[f]==2 && boundarySign[f]==0),"invalid shared native incidence");
        if(incidence[f]==1)data->report.boundaryIntegral+=boundarySign[f]*data->integrals[f];
    }
    for(long double mass:data->report.geometryMass)finite(mass,"total geometry mass");
    for(long double value:{data->report.massChange,data->report.boundaryIntegral,data->report.integratedDefect,
        data->report.l1IntegratedDefect,balanceScale})finite(value,"global material balance");
    const auto difference=finite(std::abs(data->report.integratedDefect-data->report.massChange-data->report.boundaryIntegral),"global balance difference");
    require(balanceScale>0 || difference==0,"zero global balance scale with nonzero difference");
    data->report.balanceIdentity=balanceScale>0?difference/balanceScale:0;
    require(data->report.balanceIdentity<=256*std::numeric_limits<double>::epsilon(),"global material balance identity failed");
    return HydroMaterialTransport(std::move(data));
}
HydroMaterialTransport HydroMaterialTransport::prepareSpaceTime(const QuadtreeMesh& mesh,
    Interface::HermiteGalsPublication original,HydroCharacteristicSource motion,double liquid,double gas,
    HydroSpaceTimeConfig config,const Interface::HermiteGalsBoundaryFirstJetQuery& boundary,
    const Interface::HermiteGalsPublicationConfig& publication,const HydroMaterialTransport* previousPlan)
{
    using Clock=std::chrono::steady_clock;const auto profilingStart=Clock::now();
    const auto elapsed=[](const Clock::time_point& start){
        return std::chrono::duration<double>(Clock::now()-start).count();
    };
    using namespace Interface;using Position=LAMP::TimeLayer::Position;
    original.validateMesh(mesh);motion.validateFor(mesh,motion.previous().identity(),motion.arrival().identity());
    const auto& reconstruction=motion.geometry().reconstruction();const auto& geometry=reconstruction.momentumGeometry();
    motion.previous().validateFor(mesh,geometry,original.state().identity(),original.state().timeLayer());
    require(original.state().primary().transport==HermiteGalsTransport::FirstJet&&original.state().primary().phaseVolumes.empty(),"spacetime requires independent Pure FirstJet");
    require(std::isfinite(liquid)&&liquid>0&&std::isfinite(gas)&&gas>0,"nonpositive/nonfinite density");
    require(std::isfinite(config.relativeTolerance)&&config.relativeTolerance>0&&config.relativeTolerance<1
        &&config.initialPanels>0&&config.maxDepth>=0&&config.maxDepth<=50&&config.maxGeometrySamples>2,"invalid spacetime quadrature budget");
    const auto parent=motion.arrival().timeLayer();const double begin=parent.begin,end=parent.end,dt=end-begin,middle=begin+dt*.5;
    using Rule=HydroSpaceTimeIntegrationRule;
    require(config.integrationRule==Rule::EmbeddedGaussLegendre35||config.integrationRule==Rule::EmbeddedSimpsonBoole
        ||config.integrationRule==Rule::EmbeddedSelectiveBoole,"invalid time integration rule");
    const bool selective=config.integrationRule==Rule::EmbeddedSelectiveBoole;
    const bool boole=config.integrationRule==Rule::EmbeddedSimpsonBoole||selective;const int sampleNodes=boole?2:6;
    using BoundaryTreatment=HermiteGalsBoundaryTreatment;
    require(config.firstJetBoundaryTreatment==BoundaryTreatment::InflowHistory
        ||config.firstJetBoundaryTreatment==BoundaryTreatment::PrescribedOuterFirstJet
        ||config.firstJetBoundaryTreatment==BoundaryTreatment::SymmetryFirstJet,
        "invalid FirstJet boundary treatment");
    require(config.firstJetBoundaryTreatment!=BoundaryTreatment::SymmetryFirstJet||!boundary,
        "symmetry FirstJet boundary treatment does not accept prescribed boundary data");
    const auto boundaryTreatment=boundary||config.firstJetBoundaryTreatment==BoundaryTreatment::SymmetryFirstJet
        ?config.firstJetBoundaryTreatment:BoundaryTreatment::InflowHistory;
    HydroSpaceTimeReport report;report.integrationRule=config.integrationRule;
    report.firstJetBoundaryTreatment=boundaryTreatment;
    const std::size_t atoms=reconstruction.atoms().size();
    const bool localMode=config.enableAtomLocalTimeIntegration
        &&!geometry.layout().domain().axisymmetric()&&geometry.layout().domain().uniformLeaves()
        &&supportsHermiteGalsFirstJetLocalSample(mesh,original);
    report.samplingCost.atomLocal=localMode;report.samplingCost.sourceIdentity=original.state().identity();
    report.samplingCost.previousAdvectorIdentity=motion.previous().identity();
    report.samplingCost.arrivalAdvectorIdentity=motion.arrival().identity();
    std::vector<int> allAtomIds;allAtomIds.reserve(atoms);
    for(std::size_t i=0;i<atoms;++i)allAtomIds.push_back(static_cast<int>(i));
    const auto validateLocalPartition=[&](const std::vector<TimePanel>& plan){
        std::vector<std::vector<const TimePanel*>> byAtom(atoms);
        for(const auto& panel:plan){
            require(std::isfinite(panel.lo)&&std::isfinite(panel.hi)&&panel.lo>=begin&&panel.hi<=end&&panel.hi>panel.lo,
                "atom time plan: invalid interval");
            require(panel.depth>=0&&panel.depth<=config.maxDepth,"spacetime refinement depth exhausted by previous time plan");
            require(!panel.atomIds.empty(),"atom time plan: empty atom group");int previous=-1;
            for(int id:panel.atomIds){require(id>previous&&id>=0&&static_cast<std::size_t>(id)<atoms,"atom time plan: duplicate/invalid atom");
                byAtom[id].push_back(&panel);previous=id;}
        }
        for(const auto& intervals:byAtom){auto ordered=intervals;
            std::sort(ordered.begin(),ordered.end(),[](const auto* a,const auto* b){return a->lo<b->lo;});double cursor=begin;
            for(const auto* panel:ordered){require(panel->lo==cursor,"atom time plan: incomplete or overlapping interval partition");cursor=panel->hi;}
            require(cursor==end,"atom time plan: incomplete interval partition");}
    };
    const auto mergeFull=[&](const SamplingProfile& profile){addFullTimings(report.timings,profile);addSamplingCost(report.samplingCost,profile.samplingCost);};
    const auto mergeInternal=[&](const SamplingProfile& profile){addInternalTimings(report.timings,profile);addSamplingCost(report.samplingCost,profile.samplingCost);};
    std::vector<TimePanel> initialPlan;
    if(previousPlan){
        require(previousPlan->isSpaceTime()&&!previousPlan->data_->acceptedTimePanels.empty(),
            "previous time plan: no accepted spacetime partition");
        require(previousPlan->beginGeometry().geometryIdentity()==original.state().identity()
            && &previousPlan->beginGeometry().publication().state().geometry().snapshot()==&original.state().geometry().snapshot(),
            "previous time plan: different original publication");
        require(previousPlan->layout().sameIdentity(geometry.layout())
            &&previousPlan->momentumGeometry().sameIdentity(geometry),
            "previous time plan: different mesh/layout");
        const auto& previousMotion=previousPlan->characteristicSource();
        require(previousMotion.previous().identity()==motion.previous().identity()
            &&previousMotion.geometry().sameIdentity(motion.geometry()),
            "previous time plan: different previous source");
        const auto previousClock=previousMotion.arrival().timeLayer();
        const auto previousMidpoint=previousPlan->timeLayer();
        require(previousClock.step==parent.step&&previousClock.begin==parent.begin&&previousClock.end==parent.end
            &&previousClock.position==parent.position&&previousMidpoint.step==parent.step
            &&previousMidpoint.begin==parent.begin&&previousMidpoint.end==parent.end
            &&previousMidpoint.position==Position::Midpoint,
            "previous time plan: different clock interval");
        previousPlan->validateFor(mesh,geometry,previousPlan->sources(),previousMidpoint);
        previousMotion.previous().validateFor(mesh,geometry,original.state().identity(),original.state().timeLayer());
        require(previousPlan->endGeometry().publication().state().primary().gradientUpdate==publication.gradientUpdate,
            "previous time plan: different first-jet gradient update");
        require(previousPlan->spaceTimeReport().integrationRule==report.integrationRule,
            "previous time plan: different time integration rule");
        require(previousPlan->spaceTimeReport().firstJetBoundaryTreatment==boundaryTreatment,
            "previous time plan: different FirstJet boundary treatment");
        initialPlan=previousPlan->data_->acceptedTimePanels;
        require(previousPlan->spaceTimeReport().samplingCost.atomLocal==localMode,
            "previous time plan: different atom-local route");
        if(localMode){validateLocalPartition(initialPlan);
            for(const auto& p:initialPlan)report.deepestPanel=std::max(report.deepestPanel,p.depth);
        }else{
            auto ordered=initialPlan;
            std::sort(ordered.begin(),ordered.end(),[](const auto& a,const auto& b){return a.lo<b.lo;});
            double cursor=begin;
            for(const auto& p:ordered){
                require(p.atomIds.empty(),"previous time plan: different atom-local route");
                require(std::isfinite(p.lo)&&std::isfinite(p.hi)&&p.lo==cursor&&p.hi>p.lo&&p.hi<=end,
                    "previous time plan: incomplete interval partition");
                require(p.depth>=0&&p.depth<=config.maxDepth,"spacetime refinement depth exhausted by previous time plan");
                report.deepestPanel=std::max(report.deepestPanel,p.depth);cursor=p.hi;
            }
            require(cursor==end,"previous time plan: incomplete interval partition");
        }
        report.reusedTimePlan=true;
    }
    report.initialPlanPanels=previousPlan?initialPlan.size():static_cast<std::size_t>(config.initialPanels);
    std::size_t progressPanels=0;
    const auto notifyProgress=[&](){
        if(!config.progress)return;
        LAMP::Geometry::Detail::ExactExpansion::validateActualEnvironment();
        auto clock=parent;clock.position=Position::Midpoint;
        config.progress({report.geometrySamples,progressPanels,report.deepestPanel,
            report.maxRelativePhaseErrorEstimate,elapsed(profilingStart),clock,report.samplingCost,report.integrationRule});
        LAMP::Geometry::Detail::ExactExpansion::validateActualEnvironment();
        original.validateMesh(mesh);motion.validateFor(mesh,motion.previous().identity(),motion.arrival().identity());
    };
    const auto countSample=[&](){
        require(++report.geometrySamples<=static_cast<std::size_t>(config.maxGeometrySamples),"spacetime geometry sample budget exhausted");
        if(report.geometrySamples%7==0)notifyProgress();
    };
    const auto closedProjectionSupported=[&](){
        const auto& layout=geometry.layout();const auto& domain=layout.domain();
        if(boundary||domain.axisymmetric()||!domain.uniformLeaves()||layout.pressureCount()!=mesh.cellIndexNum())return false;
        std::vector<int> walls;
        for(int u=0;u<layout.velocityCount();++u){const auto f=domain.face(mesh,layout.nativeFace(u));
            if(f.negativeCell<0||f.positiveCell<0)walls.push_back(u);}
        if(walls.empty())return false;
        for(const auto* field:{&motion.previous(),&motion.arrival()}){
            const auto* law=field->boundary();if(!law)return false;
            law->validateFor(mesh,geometry);
            if(law->velocityDofs().size()!=walls.size())return false;
            for(int u:walls)if(std::find(law->velocityDofs().begin(),law->velocityDofs().end(),u)==law->velocityDofs().end()
                ||field->values().normalVelocity()[u]!=0)return false;
            for(double q:field->tangentialVelocity())if(q!=0)return false;
        }
        return true;
    };
    const bool supported=closedProjectionSupported();
    using Scheme=HydroFirstJetTimeScheme;using Projection=HermiteGalsTimeProjection;
    require(config.firstJetTimeScheme==Scheme::BoundaryLaw||config.firstJetTimeScheme==Scheme::CharacteristicPullback
        ||config.firstJetTimeScheme==Scheme::ClosedProjectionCompensated,"invalid FirstJet time scheme");
    if(config.firstJetTimeScheme==Scheme::ClosedProjectionCompensated)
        require(supported,"compensated FirstJet time scheme requires complete uniform XY impermeable base-native sources");
    report.firstJetTimeProjection=(config.firstJetTimeScheme==Scheme::ClosedProjectionCompensated
        ||(config.firstJetTimeScheme==Scheme::BoundaryLaw&&supported))?Projection::ClosedProjectionCompensated:Projection::CharacteristicPullback;
    if(previousPlan)require(previousPlan->spaceTimeReport().firstJetTimeProjection==report.firstJetTimeProjection,
        "previous time plan: different FirstJet time projection");
    const bool pointMode=localMode&&!boundary&&config.enablePointParallel;
    const bool slabMode=localMode&&!boundary;
    std::optional<HydroAdvectorTimeSlab> localSlab;
    std::optional<HermiteGalsFirstJetLocalPreparation> localPreparation;
    if(slabMode){
        localSlab.emplace(HydroAdvectorTimeSlab::prepare(mesh,original,reconstruction,motion.previous(),motion.arrival()));
        localPreparation.emplace(HermiteGalsFirstJetLocalPreparation::prepare(mesh,original));
    }
    const auto nativeStarts=motion.prepareNativeStarts(mesh,original);
    const auto traceSlots=[](){
        int count=1;
#ifdef _OPENMP
        count=std::max(1,omp_get_max_threads());
#endif
        return std::vector<double>(static_cast<std::size_t>(count),0.);
    };
    const auto sumTrace=[](const std::vector<double>& slots){double sum=0;for(double q:slots)sum+=q;return sum;};
    const auto resolvedAt=[&](double time,std::vector<double>& traceSeconds) {
        auto clock=parent;if(time==middle)clock.position=Position::Midpoint;else clock.end=time;
        return HermiteGalsResolvedMotion(original.state().identity(),original.state().timeLayer(),motion.arrival().geometryIdentity(),clock,
            [motion,nativeStarts,clock,config,begin,&traceSeconds](Vector2D point){
                const auto traceStart=Clock::now();
                const auto t=motion.traceBackward(point,clock.time(),begin,config.characteristic,nativeStarts);
                int worker=0;
#ifdef _OPENMP
                if(omp_in_parallel())worker=omp_get_thread_num();
#endif
                traceSeconds[static_cast<std::size_t>(worker)]+=std::chrono::duration<double>(Clock::now()-traceStart).count();
                std::optional<Boundary> side;
                if(t.boundary()){const auto b=*t.boundary();side=b.normalAxis==0?(b.outwardSign<0?Boundary::west:Boundary::east):(b.outwardSign<0?Boundary::south:Boundary::north);}
                return LevelSetBoundaryTrace::JetResult{{t.point(),clock.time()-t.time(),side,t.boundary().has_value(),false},t.jacobian(),t.substeps()};
            });
    };
    const auto sampleGeometryUncounted=[&](double time,SamplingProfile& profile) {
        auto traceTime=traceSlots();
        HermiteGalsPointSamplingConfig sampling;
        if(slabMode&&config.enableNumericalTransportSampling){
            sampling.mode=HermiteGalsPointSampling::Numerical;
            sampling.valueAbsoluteTolerance=config.transportValueAbsoluteTolerance;
            sampling.gradientAbsoluteTolerance=config.transportGradientAbsoluteTolerance;
        }
        const auto step=prepareHermiteGalsFirstJetStep(mesh,original,resolvedAt(time,traceTime),boundary,publication,
            boundaryTreatment,
            report.firstJetTimeProjection,pointMode?HermiteGalsPointExecution::Parallel:HermiteGalsPointExecution::Serial,sampling);
        profile.samplingCost.maxAcceptedTransportValueError=step.report().maxAcceptedValueErrorBound;
        profile.samplingCost.maxAcceptedTransportGradientError=step.report().maxAcceptedGradientErrorBound;
        profile.fullCharacteristicSeconds+=sumTrace(traceTime);
        addFirstJetTimings(profile.fullFirstJet,step.report().timings);
        for(int cell:original.state().primary().geometryCells){(void)cell;++profile.samplingCost.fullPreparedNativeCells;}
        if(report.firstJetTimeProjection==Projection::ClosedProjectionCompensated)
            for(int cell:original.state().primary().geometryCells){(void)cell;++profile.samplingCost.fullPreparedNativeCells;}
        const auto geometryStart=Clock::now();
        auto sampled=HydroPhaseGeometry::build(mesh,geometry.metrics(),step.publication());
        profile.fullGeometrySeconds+=elapsed(geometryStart);++profile.fullSamples;
        return sampled;
    };
    const auto sampleGeometry=[&](double time) {
        countSample();
        SamplingProfile profile;
        auto sampled=sampleGeometryUncounted(time,profile);
        mergeFull(profile);
        return sampled;
    };
    const auto rateAt=[&](const HydroPhaseGeometry& g,SamplingProfile& profile) {
        const auto rateStart=Clock::now();
        auto rate=HydroPhaseRate::prepareSpatial(mesh,g,
            reconstruction.interpolate(mesh,motion.previous(),motion.arrival(),g.geometryIdentity(),g.timeLayer()));
        profile.fullRateSeconds+=elapsed(rateStart);++profile.samplingCost.fullSamples;
        profile.samplingCost.fullAtomRateEvaluations+=rate.atomRates().size();
        return rate;
    };
    const auto sampleRatesUncounted=[&](double time,SamplingProfile& profile,const std::vector<int>* requested=nullptr) {
        if(localMode){
            const auto& ids=requested?*requested:allAtomIds;
            std::optional<HydroAdvectorAtomSupport> support;
            std::vector<int> legacyCells;
            if(slabMode)support.emplace(localSlab->prepareAtomSupport(ids));
            else legacyCells=HydroPhaseRate::requiredNativeCells(mesh,reconstruction,ids);
            const auto& nativeCells=support?support->nativeCells():legacyCells;
            auto traceTime=traceSlots();
            HermiteGalsPointSamplingConfig sampling;
            if(config.enableNumericalTransportSampling){
                sampling.mode=HermiteGalsPointSampling::Numerical;
                sampling.valueAbsoluteTolerance=config.transportValueAbsoluteTolerance;
                sampling.gradientAbsoluteTolerance=config.transportGradientAbsoluteTolerance;
            }
            const auto resolved=resolvedAt(time,traceTime);
            const auto finishLocalSample=[&](const auto& sample,const auto& rateFactory) {
                profile.internalCharacteristicSeconds+=sumTrace(traceTime);
                addFirstJetTimings(profile.internalFirstJet,sample.report().timings);
                profile.samplingCost.localPreparedNativeCells+=sample.report().evaluatedNativeCells+sample.report().stagedNativeCells;
                profile.samplingCost.localScalarPointQueries+=sample.report().evaluatedPointQueries;
                const auto rateStart=Clock::now();
                auto rates=rateFactory();
                profile.samplingCost.localPointJobs+=sample.report().pointJobs;
                profile.samplingCost.localPointTraceQueries+=sample.report().pointTraceQueries;
                profile.samplingCost.pointThreads=std::max(profile.samplingCost.pointThreads,sample.report().pointThreads);
                profile.samplingCost.localNumericalPointQueries+=sample.report().numericalPointQueries;
                profile.samplingCost.localExactValueFallbacks+=sample.report().numericalValueFallbacks;
                profile.samplingCost.localExactGradientFallbacks+=sample.report().numericalGradientFallbacks;
                profile.samplingCost.maxAcceptedTransportValueError=std::max(profile.samplingCost.maxAcceptedTransportValueError,sample.report().maxAcceptedValueErrorBound);
                profile.samplingCost.maxAcceptedTransportGradientError=std::max(profile.samplingCost.maxAcceptedTransportGradientError,sample.report().maxAcceptedGradientErrorBound);
                profile.internalRateSeconds+=elapsed(rateStart);++profile.internalSamples;++profile.samplingCost.localSamples;
                profile.samplingCost.localAtomRateEvaluations+=rates.size();return rates;
            };
            if(localPreparation&&slabMode){
                std::vector<HermiteGalsFirstJetTraceSample::Segment> segments;segments.reserve(ids.size());
                for(int id:ids){const auto& atom=reconstruction.atoms()[id];
                    segments.push_back({atom.begin,atom.end,atom.normalAxis});}
                const auto sample=localPreparation->sampleTraces(resolved,nativeCells,segments,publication,
                    report.firstJetTimeProjection,pointMode?HermiteGalsPointExecution::Parallel:HermiteGalsPointExecution::Serial,sampling,boundaryTreatment);
                return finishLocalSample(sample,[&]{return HydroPhaseRate::sampleAtomRates(sample,*localSlab,*support);});
            }
            const auto sample=localPreparation
                ?localPreparation->sample(resolved,nativeCells,publication,report.firstJetTimeProjection,
                    pointMode?HermiteGalsPointExecution::Parallel:HermiteGalsPointExecution::Serial,sampling,boundaryTreatment)
                :prepareHermiteGalsFirstJetLocalSample(mesh,original,resolved,nativeCells,boundary,publication,
                    boundaryTreatment,
                    report.firstJetTimeProjection,pointMode?HermiteGalsPointExecution::Parallel:HermiteGalsPointExecution::Serial,sampling);
            return finishLocalSample(sample,[&]{return slabMode?HydroPhaseRate::sampleAtomRates(sample,*localSlab,*support):
                HydroPhaseRate::sampleAtomRates(mesh,sample,
                    reconstruction.interpolate(mesh,motion.previous(),motion.arrival(),sample.identity(),sample.timeLayer()),ids);});
        }
        require(requested==nullptr,"full time route cannot consume a local atom group");
        // Internal quadrature nodes consume only the actual Q4 face rates.
        // Midpoint/end retain the complete publication and its original gates.
        auto traceTime=traceSlots();
        const auto sample=prepareHermiteGalsFirstJetSample(mesh,original,resolvedAt(time,traceTime),boundary,publication,
            boundaryTreatment,
            report.firstJetTimeProjection);
        profile.internalCharacteristicSeconds+=sumTrace(traceTime);
        addFirstJetTimings(profile.internalFirstJet,sample.report().timings);
        for(int cell:original.state().primary().geometryCells){(void)cell;++profile.samplingCost.fullPreparedNativeCells;}
        if(report.firstJetTimeProjection==Projection::ClosedProjectionCompensated)
            for(int cell:original.state().primary().geometryCells){(void)cell;++profile.samplingCost.fullPreparedNativeCells;}
        const auto rateStart=Clock::now();
        const auto advector=reconstruction.interpolate(mesh,motion.previous(),motion.arrival(),sample.identity(),sample.timeLayer());
        auto rates=HydroPhaseRate::sampleAtomRates(mesh,sample,advector);
        profile.internalRateSeconds+=elapsed(rateStart);++profile.internalSamples;++profile.samplingCost.fullSamples;
        profile.samplingCost.fullAtomRateEvaluations+=rates.size();return rates;
    };
    const auto panelTimes=[&](double lo,double hi,double center,double width) {
        // Match the future child-center arithmetic exactly, so the two
        // quarter rates become actual cached child centers after a split.
        if(boole)return std::array<double,6>{{lo+(center-lo)*.5,center+(hi-center)*.5}};
        constexpr double nodes[4]={-.90617984593866399280,-.53846931010568309104,.53846931010568309104,.90617984593866399280};
        return std::array<double,6>{{center+width*.5*(-.77459666924148337704),center+width*.5*.77459666924148337704,
            center+width*.5*nodes[0],center+width*.5*nodes[1],center+width*.5*nodes[2],center+width*.5*nodes[3]}};
    };
    const auto eighthTimes=[](double lo,double hi,double center,const std::array<double,6>& quarters) {
        // These are the future children's quarters, using their actual
        // recursive midpoint arithmetic rather than independent width fractions.
        const std::array<double,4> times{{lo+(quarters[0]-lo)*.5,quarters[0]+(center-quarters[0])*.5,
            center+(quarters[1]-center)*.5,quarters[1]+(hi-quarters[1])*.5}};
        require(times[0]>lo&&times[0]<quarters[0]&&times[1]>quarters[0]&&times[1]<center
            &&times[2]>center&&times[2]<quarters[1]&&times[3]>quarters[1]&&times[3]<hi,
            "unresolved spacetime eighth");
        return times;
    };
    const auto booleValue=[](double width,const HydroPhaseMeasure& left,const HydroPhaseMeasure& mid,
        const HydroPhaseMeasure& right,const HydroPhaseMeasure& leftQuarter,const HydroPhaseMeasure& rightQuarter) {
        HydroPhaseMeasure value;
        const auto add=[&](const HydroPhaseMeasure& rate,long double weight) {
            value.liquid+=weight*rate.liquid;value.gas+=weight*rate.gas;
        };
        // Match one actual child panel's original Boole accumulation.
        add(left,width*(7.L/90));add(mid,width*(12.L/90));add(right,width*(7.L/90));
        add(leftQuarter,width*(32.L/90));add(rightQuarter,width*(32.L/90));
        return value;
    };
    const auto initialStart=Clock::now();
    report.timings.workflowWallSeconds[0]=elapsed(profilingStart);
    const auto initial=HydroPhaseGeometry::build(mesh,geometry.metrics(),original);
    report.timings.fullGeometrySeconds+=elapsed(initialStart);
    report.timings.workflowWallSeconds[1]=elapsed(initialStart);
    const auto fullJobsStart=Clock::now();
    std::array<std::optional<HydroPhaseGeometry>,2> fullGeometries;
    std::optional<HydroPhaseRate> midpointRate;
    std::array<std::vector<HydroPhaseMeasure>,6> prefetchedSamples;
    std::array<SamplingProfile,6> prefetchedProfiles{};
    std::array<std::exception_ptr,6> prefetchedFailures{};
    std::array<double,6> firstPanelTimes{};
    double firstPanelLo=begin;
    const bool rootPlan=!previousPlan||(initialPlan.size()==1&&initialPlan.front().lo==begin
        &&initialPlan.front().hi==end&&initialPlan.front().depth==0);
    bool batch=false,prefetchedFirstPanel=false;
#ifdef _OPENMP
    batch=!pointMode&&supported&&!boundary&&!omp_in_parallel()&&(omp_get_max_active_levels()<=1)&&omp_get_max_threads()>1
        &&rootPlan&&config.initialPanels==1&&config.maxGeometrySamples>=2+sampleNodes;
#endif
    if(batch){
        // Use the same first-panel arithmetic as below. Degenerate nodes use
        // the original path, including its midpoint cache and failure order.
        firstPanelLo=begin+dt*(static_cast<double>(0)/config.initialPanels);
        const double width=end-firstPanelLo,center=firstPanelLo+width*.5;
        firstPanelTimes=panelTimes(firstPanelLo,end,center,width);
        batch=center==middle;
        for(int q=0;q<sampleNodes;++q)batch=batch&&firstPanelTimes[q]>begin&&firstPanelTimes[q]<end&&firstPanelTimes[q]!=middle;
    }
    if(batch){
        const int batchTasks=2+sampleNodes;report.timings.fullInternalBatchTasks=batchTasks;
        LAMP::Geometry::Detail::ExactExpansion::validateEnvironment();
        std::array<SamplingProfile,2> fullProfiles{};
        std::array<std::exception_ptr,2> fullFailures{};
        std::exception_ptr midRateFailure;
        // All immutable-source jobs own their profile, result and error
        // slots. The midpoint rate has a distinct error slot: its original
        // serial position follows BOTH full geometry publications.
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif
        for(int q=0;q<batchTasks;++q){
            if(q<2){
                try {
                    LAMP::Geometry::Detail::ExactExpansion::validateEnvironment();
                    fullGeometries[q].emplace(sampleGeometryUncounted(q==0?middle:end,fullProfiles[q]));
                }catch(...){fullFailures[q]=std::current_exception();}
                if(q==0&&!fullFailures[0]){
                    try {midpointRate.emplace(rateAt(*fullGeometries[0],fullProfiles[0]));}
                    catch(...){midRateFailure=std::current_exception();}
                }
            }else{
                const int k=q-2;
                try {
                    LAMP::Geometry::Detail::ExactExpansion::validateEnvironment();
                    prefetchedSamples[k]=sampleRatesUncounted(firstPanelTimes[k],prefetchedProfiles[k]);
                }catch(...){prefetchedFailures[k]=std::current_exception();}
            }
        }
        for(int q=0;q<2;++q){
            countSample();
            if(fullFailures[q])std::rethrow_exception(fullFailures[q]);
            mergeFull(fullProfiles[q]);
        }
        if(midRateFailure)std::rethrow_exception(midRateFailure);
        prefetchedFirstPanel=true;
    }else{
        bool fullPairParallel=false;int fullJobThreads=1;
#ifdef _OPENMP
        fullPairParallel=pointMode&&supported&&!omp_in_parallel()&&(omp_get_max_active_levels()<=1)&&omp_get_max_threads()>1;
        fullJobThreads=std::max(1,omp_get_max_threads()/2);
#endif
        if(fullPairParallel){
            std::array<SamplingProfile,2> profiles{};std::array<std::exception_ptr,2> failures{};
            std::exception_ptr rateFailure;
            std::array<std::future<void>,2> jobs;
            for(int q=0;q<2;++q)jobs[q]=std::async(std::launch::async,[&,q]{
#ifdef _OPENMP
                // Each independent publication gets half the caller's workers.
                // Its existing point/face/cell loops run in a real outer team.
                omp_set_num_threads(fullJobThreads);
#endif
                try{
                    LAMP::Geometry::Detail::ExactExpansion::validateEnvironment();
                    fullGeometries[q].emplace(sampleGeometryUncounted(q==0?middle:end,profiles[q]));
                }catch(...){failures[q]=std::current_exception();}
                // This original independent midpoint job overlaps the end publication.
                // Its failure is still consumed after both geometry failures.
                if(q==0&&!failures[q]){try{midpointRate.emplace(rateAt(*fullGeometries[q],profiles[q]));}
                    catch(...){rateFailure=std::current_exception();}}
            });
            for(auto& job:jobs)job.get();
            for(int q=0;q<2;++q){countSample();if(failures[q])std::rethrow_exception(failures[q]);mergeFull(profiles[q]);}
            if(rateFailure)std::rethrow_exception(rateFailure);
        }else{
            fullGeometries[0].emplace(sampleGeometry(middle));
            fullGeometries[1].emplace(sampleGeometry(end));
            SamplingProfile profile;
            midpointRate.emplace(rateAt(*fullGeometries[0],profile));mergeFull(profile);
        }
    }
    report.timings.workflowWallSeconds[2]=elapsed(fullJobsStart);
    const auto quadratureStart=Clock::now();
    const auto& endGeometry=*fullGeometries[1];
    const auto& midpoint=*midpointRate;
    const auto evaluate=[&](double time) {
        if(time==middle)return midpoint.atomRates();
        countSample();
        SamplingProfile profile;
        auto rates=sampleRatesUncounted(time,profile);
        mergeInternal(profile);
        return rates;
    };
    const auto evaluateUncounted=[&](double time,SamplingProfile& profile) {
        return time==middle?midpoint.atomRates():sampleRatesUncounted(time,profile);
    };
    std::vector<long double> bound(atoms);
    for(std::size_t i=0;i<atoms;++i){long double speed=0;
        for(const auto* a:{&motion.previous(),&motion.arrival()}){const auto& p=a->atomPolynomials()[i];long double work=0,power=1;
            for(auto c:p.beginCoefficients){work+=std::abs(c)*power;power*=p.length;}speed=std::max(speed,work);}
        bound[i]=finite(dt*reconstruction.atoms()[i].area*speed,"spacetime atom work scale");
    }
    // Endpoint/center rates belong only to THIS actual source and prepare.
    // The retained full geometries already cover begin/middle/end. Children
    // reuse their parent's real center; warm-plan seams are evaluated afresh.
    std::map<double,std::vector<HydroPhaseMeasure>> endpointCache;
    endpointCache.emplace(middle,midpoint.atomRates());
    {
        const auto rateStart=Clock::now();
        const auto beginRate=HydroPhaseRate::prepareSpatial(mesh,initial,motion.previous());
        report.timings.fullRateSeconds+=elapsed(rateStart);++report.samplingCost.fullSamples;
        report.samplingCost.fullAtomRateEvaluations+=beginRate.atomRates().size();
        endpointCache.emplace(begin,beginRate.atomRates());
        SamplingProfile profile;
        const auto endRate=rateAt(endGeometry,profile);
        mergeFull(profile);
        endpointCache.emplace(end,endRate.atomRates());
    }
    const auto endpointRates=[&](double time)->const std::vector<HydroPhaseMeasure>& {
        const auto found=endpointCache.find(time);if(found!=endpointCache.end())return found->second;
        auto rates=evaluate(time);
        return endpointCache.emplace(time,std::move(rates)).first->second;
    };
    if(localMode){
        struct RateEntry {std::vector<HydroPhaseMeasure> values;std::vector<unsigned char> present;};
        std::map<double,RateEntry> rateCache;
        for(const auto& entry:endpointCache)rateCache.emplace(entry.first,RateEntry{entry.second,std::vector<unsigned char>(atoms,1)});
        const auto missingAt=[&](double time,const std::vector<int>& ids){
            std::vector<int> missing;const auto found=rateCache.find(time);
            for(int id:ids)if(found==rateCache.end()||!found->second.present[id])missing.push_back(id);
            return missing;
        };
        const auto storeRates=[&](double time,const std::vector<int>& ids,const std::vector<HydroPhaseMeasure>& rates){
            require(ids.size()==rates.size(),"local rate cache layout mismatch");
            auto found=rateCache.find(time);
            if(found==rateCache.end())found=rateCache.emplace(time,RateEntry{std::vector<HydroPhaseMeasure>(atoms),std::vector<unsigned char>(atoms,0)}).first;
            for(std::size_t k=0;k<ids.size();++k){const int id=ids[k];
                require(!found->second.present[id],"local rate cache overwrote actual atom");
                found->second.values[id]=rates[k];found->second.present[id]=1;}
        };
        const auto cachedRates=[&](double time,const std::vector<int>& ids){
            const auto found=rateCache.find(time);require(found!=rateCache.end(),"missing actual local time sample");
            std::vector<HydroPhaseMeasure> rates;rates.reserve(ids.size());
            for(int id:ids){require(found->second.present[id],"missing actual local atom rate");rates.push_back(found->second.values[id]);}
            return rates;
        };
        const auto evaluateRequested=[&](double time,const std::vector<int>& ids){
            const auto missing=missingAt(time,ids);
            if(!missing.empty()){
                countSample();SamplingProfile profile;
                const auto rates=sampleRatesUncounted(time,profile,&missing);
                mergeInternal(profile);storeRates(time,missing,rates);
            }
            return cachedRates(time,ids);
        };
        struct AtomPanel {double lo,hi;int depth;std::vector<int> ids;std::vector<HydroPhaseMeasure> value,error;};
        const auto atomPanel=[&](double lo,double hi,int depth,const std::vector<int>& ids){
            require(hi>lo&&!ids.empty(),"unresolved or empty atom spacetime panel");const double width=hi-lo,center=lo+width*.5;
            require(center>lo&&center<hi,"unresolved spacetime center");
            const auto leftEndpoint=evaluateRequested(lo,ids),rightEndpoint=evaluateRequested(hi,ids),mid=evaluateRequested(center,ids);
            std::vector<HydroPhaseMeasure> coarse(ids.size()),fine(ids.size());
            std::vector<unsigned char> liquidSeen(ids.size()),gasSeen(ids.size());
            const auto add=[&](auto& sum,const auto& value,long double weight){
                require(value.size()==ids.size(),"local panel atom layout mismatch");
                for(std::size_t k=0;k<ids.size();++k){liquidSeen[k]|=value[k].liquid!=0;gasSeen[k]|=value[k].gas!=0;
                    sum[k].liquid+=weight*value[k].liquid;sum[k].gas+=weight*value[k].gas;}
            };
            if(boole){
                // Endpoints contribute to both rules, but not to the interior
                // support inventory used by the original unseen-phase guard.
                const auto addEndpoint=[&](auto& sum,const auto& value,long double weight){
                    require(value.size()==ids.size(),"local panel atom layout mismatch");
                    for(std::size_t k=0;k<ids.size();++k){sum[k].liquid+=weight*value[k].liquid;sum[k].gas+=weight*value[k].gas;}};
                addEndpoint(coarse,leftEndpoint,width/6.L);add(coarse,mid,width*(4.L/6));addEndpoint(coarse,rightEndpoint,width/6.L);
                addEndpoint(fine,leftEndpoint,width*(7.L/90));add(fine,mid,width*(12.L/90));addEndpoint(fine,rightEndpoint,width*(7.L/90));
            }else{add(coarse,mid,width*(4.L/9));add(fine,mid,width*(64.L/225));}
            constexpr long double weights[4]={.11846344252809454376L,.23931433524968323402L,.23931433524968323402L,.11846344252809454376L};
            const auto addNode=[&](int q,const auto& rates){
                if(boole)add(fine,rates,width*(32.L/90));
                else if(q<2)add(coarse,rates,width*(5.L/18));else add(fine,rates,width*weights[q-2]);};
            const auto times=panelTimes(lo,hi,center,width);
            if(boole)require(times[0]>lo&&times[0]<center&&times[1]>center&&times[1]<hi,"unresolved spacetime quarter");
            const bool usePrefetched=prefetchedFirstPanel&&depth==0&&lo==firstPanelLo&&hi==end&&times==firstPanelTimes;
            std::array<std::vector<int>,6> missing;
            std::size_t needed=0;bool distinct=true;
            for(int q=0;q<sampleNodes;++q){missing[q]=missingAt(times[q],ids);if(!missing[q].empty())++needed;
                for(int j=0;j<q;++j)if(times[q]==times[j])distinct=false;}
            bool parallel=false;
#ifdef _OPENMP
            parallel=!pointMode&&supported&&!boundary&&!omp_in_parallel()&&(omp_get_max_active_levels()<=1)&&omp_get_max_threads()>1&&distinct
                &&needed<=static_cast<std::size_t>(config.maxGeometrySamples)-report.geometrySamples;
#endif
            if(usePrefetched){
                // These actual all-atom jobs were evaluated alongside the
                // full midpoint/end jobs. Only the caller commits their cache.
                for(int q=0;q<sampleNodes;++q){countSample();
                    if(prefetchedFailures[q])std::rethrow_exception(prefetchedFailures[q]);
                    mergeInternal(prefetchedProfiles[q]);storeRates(times[q],allAtomIds,prefetchedSamples[q]);
                    const auto rates=cachedRates(times[q],ids);
                    addNode(q,rates);}
                prefetchedFirstPanel=false;
            }else if(parallel){
                LAMP::Geometry::Detail::ExactExpansion::validateEnvironment();
                std::array<std::vector<HydroPhaseMeasure>,6> samples;
                std::array<SamplingProfile,6> profiles{};std::array<std::exception_ptr,6> failures{};
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
                for(int q=0;q<sampleNodes;++q)if(!missing[q].empty()){
                    try {LAMP::Geometry::Detail::ExactExpansion::validateEnvironment();
                        samples[q]=sampleRatesUncounted(times[q],profiles[q],&missing[q]);}
                    catch(...){failures[q]=std::current_exception();}
                }
                for(int q=0;q<sampleNodes;++q){if(!missing[q].empty()){
                        countSample();if(failures[q])std::rethrow_exception(failures[q]);
                        mergeInternal(profiles[q]);storeRates(times[q],missing[q],samples[q]);}
                    const auto rates=cachedRates(times[q],ids);
                    addNode(q,rates);}
            }else{
                if(boole){for(int q=0;q<sampleNodes;++q)addNode(q,evaluateRequested(times[q],ids));}
                else{
                    for(int q=0;q<2;++q)add(coarse,evaluateRequested(times[q],ids),width*(5.L/18));
                    for(int q=0;q<4;++q)add(fine,evaluateRequested(times[q+2],ids),width*weights[q]);}
            }
            std::vector<HydroPhaseMeasure> error(ids.size());bool endpointSupport=false;
            const auto estimateErrors=[&](){
                endpointSupport=false;
                for(std::size_t k=0;k<ids.size();++k){const auto id=static_cast<std::size_t>(ids[k]);
                    error[k]={std::abs(fine[k].liquid-coarse[k].liquid),std::abs(fine[k].gas-coarse[k].gas)};
                    const bool missingLiquid=!liquidSeen[k]&&(leftEndpoint[k].liquid!=0||rightEndpoint[k].liquid!=0);
                    const bool missingGas=!gasSeen[k]&&(leftEndpoint[k].gas!=0||rightEndpoint[k].gas!=0);
                    if(missingLiquid||missingGas){
                        require(bound[id]>0,"zero work with endpoint phase support");
                        const auto guard=finite((static_cast<long double>(width)/dt)*bound[id],"endpoint phase support bound");
                        require(guard>0,"unrepresentable endpoint phase support bound");
                        if(missingLiquid)error[k].liquid=std::max(error[k].liquid,guard);
                        if(missingGas)error[k].gas=std::max(error[k].gas,guard);endpointSupport=true;
                    }
                }
            };
            estimateErrors();
            if(selective){
                std::vector<int> eligibleIds;std::vector<std::size_t> eligiblePositions;
                for(std::size_t k=0;k<ids.size();++k){
                    const auto budget=config.relativeTolerance*(static_cast<long double>(width)/dt)*bound[ids[k]];
                    if(std::max(error[k].liquid,error[k].gas)>budget){
                        eligibleIds.push_back(ids[k]);eligiblePositions.push_back(k);
                    }
                }
                if(!eligibleIds.empty()){
                    // Keep the original five-node Boole value as coarse.
                    // Only these atoms pay for four additional actual samples.
                    const auto extraTimes=eighthTimes(lo,hi,center,times);
                    const auto leftMid=cachedRates(times[0],eligibleIds),rightMid=cachedRates(times[1],eligibleIds);
                    std::array<std::vector<HydroPhaseMeasure>,4> extra;
                    for(int q=0;q<4;++q)extra[q]=evaluateRequested(extraTimes[q],eligibleIds);
                    for(std::size_t j=0;j<eligibleIds.size();++j){const auto k=eligiblePositions[j];
                        coarse[k]=fine[k];
                        const auto left=booleValue(center-lo,leftEndpoint[k],leftMid[j],mid[k],extra[0][j],extra[1][j]);
                        const auto right=booleValue(hi-center,mid[k],rightMid[j],rightEndpoint[k],extra[2][j],extra[3][j]);
                        fine[k]={left.liquid+right.liquid,left.gas+right.gas};
                        for(int q=0;q<4;++q){liquidSeen[k]|=extra[q][j].liquid!=0;gasSeen[k]|=extra[q][j].gas!=0;}
                    }
                    // Do not retain a five-node unseen-phase guard after a
                    // real added interior node has supplied that phase.
                    estimateErrors();
                }
            }
            if(endpointSupport)++report.endpointSupportPanels;
            return AtomPanel{lo,hi,depth,ids,std::move(fine),std::move(error)};
        };
        std::vector<AtomPanel> panels;
        if(previousPlan){for(const auto& p:initialPlan){panels.push_back(atomPanel(p.lo,p.hi,p.depth,p.atomIds));
                progressPanels=panels.size();notifyProgress();}}
        else{for(int i=0;i<config.initialPanels;++i){panels.push_back(atomPanel(begin+dt*(static_cast<double>(i)/config.initialPanels),
                i+1==config.initialPanels?end:begin+dt*(static_cast<double>(i+1)/config.initialPanels),0,allAtomIds));
                progressPanels=panels.size();notifyProgress();}}
        std::vector<HydroPhaseMeasure> integrated(atoms);
        for(;;){
            std::vector<HydroPhaseMeasure> error(atoms);integrated.assign(atoms,{});
            for(const auto& p:panels)for(std::size_t k=0;k<p.ids.size();++k){const auto id=static_cast<std::size_t>(p.ids[k]);
                integrated[id].liquid+=p.value[k].liquid;integrated[id].gas+=p.value[k].gas;
                error[id].liquid+=p.error[k].liquid;error[id].gas+=p.error[k].gas;}
            long double relative=0,absolute=0;std::vector<unsigned char> overdue(atoms,0);
            for(std::size_t i=0;i<atoms;++i){const auto e=std::max(error[i].liquid,error[i].gas);
                require(bound[i]>0||e==0,"zero work with nonzero phase integration error");
                absolute+=error[i].liquid+error[i].gas;
                if(bound[i]>0){const auto r=e/bound[i];relative=std::max(relative,r);overdue[i]=r>config.relativeTolerance;}}
            require(std::isfinite(relative)&&std::isfinite(absolute),"nonfinite spacetime error estimate");
            report.absolutePhaseErrorEstimate=absolute;report.maxRelativePhaseErrorEstimate=relative;notifyProgress();
            if(relative<=config.relativeTolerance)break;
            const auto localEligible=[&](const AtomPanel& panel,std::size_t k){
                const auto id=static_cast<std::size_t>(panel.ids[k]);
                if(!overdue[id])return false;
                const auto budget=config.relativeTolerance*(static_cast<long double>(panel.hi-panel.lo)/dt)*bound[id];
                return std::max(panel.error[k].liquid,panel.error[k].gas)>budget;
            };
            std::size_t worst=0,originalWorst=0;long double worstScore=-1,originalWorstScore=-1;
            for(std::size_t p=0;p<panels.size();++p){long double score=-1,originalScore=-1;
                for(std::size_t k=0;k<panels[p].ids.size();++k){const auto id=static_cast<std::size_t>(panels[p].ids[k]);
                    if(overdue[id]){
                        const auto contribution=std::max(panels[p].error[k].liquid,panels[p].error[k].gas)/bound[id];
                        originalScore=std::max(originalScore,contribution);
                        if(localEligible(panels[p],k))score=std::max(score,contribution);
                    }}
                if(score>worstScore){worst=p;worstScore=score;}
                if(originalScore>originalWorstScore){originalWorst=p;originalWorstScore=originalScore;}}
            // With an exact interval partition, an overdue atom must have an
            // eligible contribution. Keep the original scheduling at a budget
            // rounding boundary or when the filtered choice is depth-limited.
            bool pruneContributions=true;
            if(worstScore<=0||panels[worst].depth>=config.maxDepth){
                worst=originalWorst;worstScore=originalWorstScore;pruneContributions=false;
            }
            require(worstScore>0,"unresolved atom-local cumulative error");
            const auto lo=panels[worst].lo,hi=panels[worst].hi;const int depth=panels[worst].depth+1;
            require(depth<=config.maxDepth,"spacetime refinement depth exhausted");report.deepestPanel=std::max(report.deepestPanel,depth);
            std::vector<int> refineIds;AtomPanel retained{lo,hi,panels[worst].depth,{},{},{}};
            for(std::size_t k=0;k<panels[worst].ids.size();++k){const int id=panels[worst].ids[k];
                if(pruneContributions?localEligible(panels[worst],k):overdue[id]!=0)refineIds.push_back(id);
                else{retained.ids.push_back(id);retained.value.push_back(panels[worst].value[k]);retained.error.push_back(panels[worst].error[k]);}}
            require(!refineIds.empty(),"empty atom-local refinement");
            const double center=lo+(hi-lo)*.5;auto left=atomPanel(lo,center,depth,refineIds);notifyProgress();
            auto right=atomPanel(center,hi,depth,refineIds);
            if(retained.ids.empty()){panels[worst]=std::move(left);panels.push_back(std::move(right));}
            else{panels[worst]=std::move(retained);panels.push_back(std::move(left));panels.push_back(std::move(right));}
            progressPanels=panels.size();notifyProgress();
        }
        report.timings.workflowWallSeconds[3]=elapsed(quadratureStart);
        const auto finalizeStart=Clock::now();
        report.acceptedPanels=panels.size();progressPanels=panels.size();
        std::vector<TimePanel> acceptedTimePanels;acceptedTimePanels.reserve(panels.size());
        for(const auto& p:panels)acceptedTimePanels.push_back({p.lo,p.hi,p.depth,p.ids});
        validateLocalPartition(acceptedTimePanels); // Every atom covers the exact full interval, with no overlap or gap.
        original.validateMesh(mesh);motion.validateFor(mesh,motion.previous().identity(),motion.arrival().identity());
        const auto flux=HydroMaterialFlux::prepare(mesh,midpoint,liquid,gas);notifyProgress();
        report.timings.workflowWallSeconds[4]=elapsed(finalizeStart);
        return assemble(mesh,initial,midpoint,endGeometry,liquid,gas,&flux,&integrated,&motion,&report,&profilingStart,&acceptedTimePanels);
    }
    struct Panel {double lo,hi;int depth;std::vector<HydroPhaseMeasure> value,error;};
    const auto panel=[&](double lo,double hi,int depth){
        require(hi>lo,"unresolved spacetime panel");const double width=hi-lo,center=lo+width*.5;
        require(center>lo&&center<hi,"unresolved spacetime center");
        const auto& leftEndpoint=endpointRates(lo);const auto& rightEndpoint=endpointRates(hi);
        const auto& mid=endpointRates(center);
        std::vector<HydroPhaseMeasure> coarse(atoms),fine(atoms);
        std::vector<unsigned char> liquidSeen(atoms),gasSeen(atoms);
        const auto add=[&](auto& sum,const auto& value,long double weight){for(std::size_t i=0;i<atoms;++i){
            liquidSeen[i]|=value[i].liquid!=0;gasSeen[i]|=value[i].gas!=0;
            sum[i].liquid+=weight*value[i].liquid;sum[i].gas+=weight*value[i].gas;}};
        if(boole){
            const auto addEndpoint=[&](auto& sum,const auto& value,long double weight){for(std::size_t i=0;i<atoms;++i){
                sum[i].liquid+=weight*value[i].liquid;sum[i].gas+=weight*value[i].gas;}};
            addEndpoint(coarse,leftEndpoint,width/6.L);add(coarse,mid,width*(4.L/6));addEndpoint(coarse,rightEndpoint,width/6.L);
            addEndpoint(fine,leftEndpoint,width*(7.L/90));add(fine,mid,width*(12.L/90));addEndpoint(fine,rightEndpoint,width*(7.L/90));
        }else{add(coarse,mid,width*(4.L/9));add(fine,mid,width*(64.L/225));}
        constexpr long double weights[4]={.11846344252809454376L,.23931433524968323402L,.23931433524968323402L,.11846344252809454376L};
        const auto addNode=[&](int q,const auto& rates){
            if(boole)add(fine,rates,width*(32.L/90));
            else if(q<2)add(coarse,rates,width*(5.L/18));else add(fine,rates,width*weights[q-2]);};
        const auto times=panelTimes(lo,hi,center,width);
        if(boole)require(times[0]>lo&&times[0]<center&&times[1]>center&&times[1]<hi,"unresolved spacetime quarter");
        const bool usePrefetched=prefetchedFirstPanel&&depth==0&&lo==firstPanelLo&&hi==end&&times==firstPanelTimes;
        std::array<bool,6> cached{};std::size_t needed=0;
        for(int q=0;q<sampleNodes;++q){cached[q]=boole&&endpointCache.find(times[q])!=endpointCache.end();
            if(!cached[q]&&times[q]!=middle)++needed;}
        bool parallel=false;
#ifdef _OPENMP
        parallel=!pointMode&&supported&&!boundary&&!omp_in_parallel()&&(omp_get_max_active_levels()<=1)&&omp_get_max_threads()>1
            &&needed<=static_cast<std::size_t>(config.maxGeometrySamples)-report.geometrySamples;
#endif
        if(parallel||usePrefetched){
            // Each typed closed-source time has its own FirstJet map/cache and
            // immutable Q4 sample. Reserve only a fully affordable panel, so a work
            // budget failure is still visited at the original serial slot.
            LAMP::Geometry::Detail::ExactExpansion::validateEnvironment();
            std::array<std::vector<HydroPhaseMeasure>,6> samples;
            // Every worker owns one slot; only the caller merges diagnostics
            // after the team, in the original failure and quadrature order.
            std::array<SamplingProfile,6> profiles{};
            std::array<std::exception_ptr,6> failures{};
            if(usePrefetched){
                samples=std::move(prefetchedSamples);profiles=std::move(prefetchedProfiles);
                failures=std::move(prefetchedFailures);prefetchedFirstPanel=false;
            }else{
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
                for(int q=0;q<sampleNodes;++q){
                    try {
                        LAMP::Geometry::Detail::ExactExpansion::validateEnvironment();
                        if(cached[q])samples[q]=endpointCache.at(times[q]);
                        else samples[q]=evaluateUncounted(times[q],profiles[q]);
                    }catch(...){failures[q]=std::current_exception();}
                }
            }
            for(int q=0;q<sampleNodes;++q){
                if(!cached[q]&&times[q]!=middle)countSample();
                if(failures[q])std::rethrow_exception(failures[q]);
                if(!cached[q])mergeInternal(profiles[q]);
                if(boole&&!cached[q])endpointCache.emplace(times[q],samples[q]);
                addNode(q,samples[q]);
            }
        }else{
            // Includes boundary callbacks, active/nested teams, and panels
            // near the remaining sample budget. Preserve their original order.
            if(boole){for(int q=0;q<sampleNodes;++q)addNode(q,endpointRates(times[q]));}
            else{
                for(int q=0;q<2;++q)add(coarse,evaluate(times[q]),width*(5.L/18));
                for(int q=0;q<4;++q)add(fine,evaluate(times[q+2]),width*weights[q]);}
        }
        std::vector<HydroPhaseMeasure> error(atoms);bool endpointSupport=false;
        const auto estimateErrors=[&](){
            endpointSupport=false;
            for(std::size_t i=0;i<atoms;++i){
                error[i]={std::abs(fine[i].liquid-coarse[i].liquid),std::abs(fine[i].gas-coarse[i].gas)};
                const bool missingLiquid=!liquidSeen[i]&&(leftEndpoint[i].liquid!=0||rightEndpoint[i].liquid!=0);
                const bool missingGas=!gasSeen[i]&&(leftEndpoint[i].gas!=0||rightEndpoint[i].gas!=0);
                if(missingLiquid||missingGas){
                    require(bound[i]>0,"zero work with endpoint phase support");
                    // A phase unseen by all interior rule nodes can enter near either
                    // end. Bound its whole panel by the ORIGINAL slab work scale;
                    // endpoint values themselves are not assumed to bound it.
                    const auto guard=finite((static_cast<long double>(width)/dt)*bound[i],"endpoint phase support bound");
                    require(guard>0,"unrepresentable endpoint phase support bound");
                    if(missingLiquid)error[i].liquid=std::max(error[i].liquid,guard);
                    if(missingGas)error[i].gas=std::max(error[i].gas,guard);
                    endpointSupport=true;
                }
            }
        };
        estimateErrors();
        if(selective){
            std::vector<std::size_t> eligible;
            for(std::size_t i=0;i<atoms;++i){
                const auto budget=config.relativeTolerance*(static_cast<long double>(width)/dt)*bound[i];
                if(std::max(error[i].liquid,error[i].gas)>budget)eligible.push_back(i);
            }
            if(!eligible.empty()){
                const auto extraTimes=eighthTimes(lo,hi,center,times);
                const auto& leftMid=endpointRates(times[0]);const auto& rightMid=endpointRates(times[1]);
                // Full geometry routes keep their actual all-atom samples,
                // but upgrade the accepted value only for eligible atoms.
                std::array<const std::vector<HydroPhaseMeasure>*,4> extra{};
                for(int q=0;q<4;++q)extra[q]=&endpointRates(extraTimes[q]);
                for(std::size_t i:eligible){
                    coarse[i]=fine[i];
                    const auto left=booleValue(center-lo,leftEndpoint[i],leftMid[i],mid[i],(*extra[0])[i],(*extra[1])[i]);
                    const auto right=booleValue(hi-center,mid[i],rightMid[i],rightEndpoint[i],(*extra[2])[i],(*extra[3])[i]);
                    fine[i]={left.liquid+right.liquid,left.gas+right.gas};
                    for(int q=0;q<4;++q){liquidSeen[i]|=(*extra[q])[i].liquid!=0;gasSeen[i]|=(*extra[q])[i].gas!=0;}
                }
                estimateErrors();
            }
        }
        if(endpointSupport)++report.endpointSupportPanels;
        return Panel{lo,hi,depth,std::move(fine),std::move(error)};
    };
    std::vector<Panel> panels;
    if(previousPlan){
        for(const auto& p:initialPlan){
            panels.push_back(panel(p.lo,p.hi,p.depth));progressPanels=panels.size();notifyProgress();
        }
    }else{
        for(int i=0;i<config.initialPanels;++i){
            panels.push_back(panel(begin+dt*(static_cast<double>(i)/config.initialPanels),
                i+1==config.initialPanels?end:begin+dt*(static_cast<double>(i+1)/config.initialPanels),0));
            progressPanels=panels.size();notifyProgress();
        }
    }
    std::vector<HydroPhaseMeasure> integrated(atoms);
    for(;;){
        std::vector<HydroPhaseMeasure> error(atoms);integrated.assign(atoms,{});std::size_t worst=0;long double worstScore=-1;
        for(std::size_t p=0;p<panels.size();++p){long double score=0;const auto& item=panels[p];
            for(std::size_t i=0;i<atoms;++i){integrated[i].liquid+=item.value[i].liquid;integrated[i].gas+=item.value[i].gas;
                error[i].liquid+=item.error[i].liquid;error[i].gas+=item.error[i].gas;
                const auto e=std::max(item.error[i].liquid,item.error[i].gas);
                require(bound[i]>0||e==0,"zero work with nonzero phase integration error");if(bound[i]>0)score=std::max(score,e/bound[i]);}
            if(score>worstScore){worst=p;worstScore=score;}
        }
        long double relative=0,absolute=0;
        for(std::size_t i=0;i<atoms;++i){absolute+=error[i].liquid+error[i].gas;if(bound[i]>0)relative=std::max(relative,std::max(error[i].liquid,error[i].gas)/bound[i]);}
        require(std::isfinite(relative)&&std::isfinite(absolute),"nonfinite spacetime error estimate");
        report.absolutePhaseErrorEstimate=absolute;report.maxRelativePhaseErrorEstimate=relative;notifyProgress();
        if(relative<=config.relativeTolerance)break;
        const auto lo=panels[worst].lo,hi=panels[worst].hi;const int depth=panels[worst].depth+1;
        require(depth<=config.maxDepth,"spacetime refinement depth exhausted");report.deepestPanel=std::max(report.deepestPanel,depth);
        const double center=lo+(hi-lo)*.5;auto left=panel(lo,center,depth);notifyProgress();
        auto right=panel(center,hi,depth);
        panels[worst]=std::move(left);panels.push_back(std::move(right));
        progressPanels=panels.size();notifyProgress();
    }
    report.timings.workflowWallSeconds[3]=elapsed(quadratureStart);
    const auto finalizeStart=Clock::now();
    report.acceptedPanels=panels.size();
    std::vector<TimePanel> acceptedTimePanels;acceptedTimePanels.reserve(panels.size());
    for(const auto& p:panels)acceptedTimePanels.push_back({p.lo,p.hi,p.depth});
    // The whole immutable-source batch is checked again before publication;
    // no stale mesh or endpoint identity can escape this prepared consumer.
    original.validateMesh(mesh);motion.validateFor(mesh,motion.previous().identity(),motion.arrival().identity());
    const auto flux=HydroMaterialFlux::prepare(mesh,midpoint,liquid,gas);
    notifyProgress();
    report.timings.workflowWallSeconds[4]=elapsed(finalizeStart);
    return assemble(mesh,initial,midpoint,endGeometry,liquid,gas,&flux,&integrated,&motion,&report,&profilingStart,&acceptedTimePanels);
}
HydroMaterialTransport HydroMaterialTransport::assemble(const QuadtreeMesh& mesh,HydroPhaseGeometry begin,
    HydroPhaseRate midpoint,HydroPhaseGeometry end,double liquid,double gas,const HydroMaterialFlux* flux,
    const std::vector<HydroPhaseMeasure>* integrated,const HydroCharacteristicSource* motion,const HydroSpaceTimeReport* report,
    const std::chrono::steady_clock::time_point* profilingStart,const std::vector<TimePanel>* acceptedTimePanels)
{
    const auto assemblyStart=std::chrono::steady_clock::now();
    using Position=LAMP::TimeLayer::Position;
    const auto& mid=midpoint.geometry();
    for(const auto* geometry:std::array<const HydroPhaseGeometry*,3>{{&begin,&mid,&end}})
        geometry->validateFor(mesh,geometry->geometryIdentity(),geometry->timeLayer());
    require(begin.layout().sameIdentity(mid.layout()) && end.layout().sameIdentity(mid.layout()),
        "different native layout; explicit transfer/rebuild required");
    require(&begin.metrics().divergence().matrix()==&mid.metrics().divergence().matrix() &&
        &end.metrics().divergence().matrix()==&mid.metrics().divergence().matrix(),"different retained D_I");
    require(std::isfinite(liquid) && liquid>0 && std::isfinite(gas) && gas>0,"nonpositive/nonfinite density");
    const auto a=begin.timeLayer(),m=mid.timeLayer(),b=end.timeLayer();
    const double tm=m.time(),ta=a.time(),tb=b.time();
    require(m.position==Position::Midpoint && b.position==Position::End,"requires Midpoint/End candidates");
    require(m.step==b.step && m.begin==b.begin && m.end==b.end,"different full candidate interval");
    require((a.position==Position::Initial || a.position==Position::End) && a.step==m.step-1 && ta==m.begin,
        "begin geometry is not the previous Initial/End layer");
    require(tm>m.begin && tm<m.end && tb==m.end,"unrepresentable midpoint interval");
    const double dt=m.end-m.begin;
    require(flux!=nullptr && midpoint.isSpatial(),"requires prepared spatial material flux");
    flux->validateFor(mesh,midpoint.momentumGeometry(),mid.geometryIdentity(),midpoint.advectorIdentity(),m);
    require(&flux->phaseRate().atomRates()==&midpoint.atomRates(),"different original phase-rate owner");
    auto data=std::make_shared<Data>(std::move(begin),std::move(midpoint),std::move(end),liquid,gas,dt,flux);
    const auto& layout=data->midpoint->geometry().layout();
    if(integrated) {
        require(motion&&report&&acceptedTimePanels&&!acceptedTimePanels->empty()
            &&integrated->size()==data->spatial->atomRates().size(),"invalid spacetime material provenance");
        data->motion=*motion;data->spaceTime=*report;data->phaseIntegrals=*integrated;
        data->acceptedTimePanels=*acceptedTimePanels;
        for(const auto& q:*integrated){const auto f=densitySum(q,liquid,gas,"spacetime atom material integral");
            data->atomIntegrals.push_back(f);data->atoms.push_back(finite(f/dt,"average atom material rate"));}
        const auto& r=data->spatial->phaseRate().spatialAdvector().reconstruction();
        for(const auto& ids:r.primalAtoms()){const auto f=atomSum(data->atomIntegrals,ids);data->integrals.push_back(f);data->faces.push_back(f/dt);}
        for(const auto& ids:r.dualAtoms()){const auto f=atomSum(data->atomIntegrals,ids);data->dualIntegrals.push_back(f);data->dual.push_back(f/dt);}
    } else {
        data->atoms=data->spatial->atomRates();data->faces=data->spatial->faceRates();data->dual=data->spatial->dualFaceRates();
        const auto integrate=[&](const auto& rates,auto& integrals,const std::string& kind){for(auto rate:rates){
            const auto f=finite(dt*rate,kind+" material integral");require(rate==0||f!=0,"unrepresentable "+kind+" material integral");integrals.push_back(f);}};
        // Preserve the original primal failure ordering and diagnostic.
        integrate(data->faces,data->integrals,"face");integrate(data->atoms,data->atomIntegrals,"atom");integrate(data->dual,data->dualIntegrals,"dual");
    }
    data->cells.resize(layout.pressureCount());
    for(int p=0;p<layout.pressureCount();++p){
        auto& cell=data->cells[p];
        const auto va=data->begin.cellVolumes()[p],vm=data->midpoint->geometry().cellVolumes()[p],vb=data->end.cellVolumes()[p];
        cell.massBegin=densitySum(va,liquid,gas,"begin geometry mass");
        cell.massMidpoint=densitySum(vm,liquid,gas,"midpoint geometry mass");
        cell.massEnd=densitySum(vb,liquid,gas,"end geometry mass");
        require(cell.massBegin>0 && cell.massMidpoint>0 && cell.massEnd>0,"unrepresentable positive geometry mass");
        // Difference phase measurements first; never subtract two large
        // independently rounded material masses or manufacture an inventory.
        cell.massChange=densitySum({vb.liquid-va.liquid,vb.gas-va.gas},liquid,gas,"geometry mass change");
        data->report.geometryMass[0]+=cell.massBegin;data->report.geometryMass[1]+=cell.massMidpoint;data->report.geometryMass[2]+=cell.massEnd;
        data->report.massChange+=cell.massChange;
    }
    std::vector<int> incidence(layout.velocityCount()),boundarySign(layout.velocityCount());
    long double balanceScale=0;
    const auto& di=data->midpoint->geometry().metrics().divergence().matrix();
    for(Eigen::Index p=0;p<di.outerSize();++p){
        auto& cell=data->cells[static_cast<std::size_t>(p)];long double gross=std::abs(cell.massChange);
        for(HydroIntegralDivergence::Matrix::InnerIterator it(di,p);it;++it){
            const int f=static_cast<int>(it.col());
            const double area=data->midpoint->geometry().metrics().faceMeasures()[f];
            require(area>0 && std::abs(it.value())==area,"different retained native face measure");
            const int sign=it.value()>0?1:-1;
            cell.outwardIntegral+=sign*data->integrals[f];gross+=std::abs(data->integrals[f]);
            ++incidence[f];boundarySign[f]+=sign;
        }
        finite(cell.outwardIntegral,"cell outward material integral");finite(gross,"cell defect scale");
        cell.integratedDefect=finite(cell.massChange+cell.outwardIntegral,"integrated material defect");
        cell.defectRate=finite(cell.integratedDefect/dt,"material continuity defect rate");
        require(gross>0 || cell.integratedDefect==0,"zero defect scale with nonzero defect");
        cell.relativeDefect=gross>0?finite(std::abs(cell.integratedDefect)/gross,"relative material defect"):0;
        data->report.integratedDefect+=cell.integratedDefect;
        data->report.l1IntegratedDefect+=std::abs(cell.integratedDefect);
        data->report.maxAbsDefectRate=std::max(data->report.maxAbsDefectRate,std::abs(cell.defectRate));
        data->report.maxRelativeCellDefect=std::max(data->report.maxRelativeCellDefect,cell.relativeDefect);
        balanceScale+=gross;
    }
    for(int f=0;f<layout.velocityCount();++f){
        if(incidence[f]==0){
            require(data->midpoint->geometry().metrics().faceMeasures()[f]==0 && data->integrals[f]==0,
                "unrepresented positive-area face");
        }else require((incidence[f]==1 && std::abs(boundarySign[f])==1) ||
            (incidence[f]==2 && boundarySign[f]==0),"invalid shared native incidence");
        if(incidence[f]==1)data->report.boundaryIntegral+=boundarySign[f]*data->integrals[f];
    }
    for(long double mass:data->report.geometryMass)finite(mass,"total geometry mass");
    for(long double value:{data->report.massChange,data->report.boundaryIntegral,data->report.integratedDefect,
        data->report.l1IntegratedDefect,balanceScale})finite(value,"global material balance");
    const auto difference=finite(std::abs(data->report.integratedDefect-data->report.massChange-data->report.boundaryIntegral),"global balance difference");
    require(balanceScale>0 || difference==0,"zero global balance scale with nonzero difference");
    data->report.balanceIdentity=balanceScale>0?difference/balanceScale:0;
    require(data->report.balanceIdentity<=256*std::numeric_limits<double>::epsilon(),"global material balance identity failed");
    if(profilingStart)data->spaceTime.timings.workflowWallSeconds[5]=
        std::chrono::duration<double>(std::chrono::steady_clock::now()-assemblyStart).count();
    if(profilingStart)data->spaceTime.timings.totalWallSeconds=
        std::chrono::duration<double>(std::chrono::steady_clock::now()-*profilingStart).count();
    return HydroMaterialTransport(std::move(data));
}
HydroPhaseGeometry::Model HydroMaterialTransport::model() const noexcept{return data_->begin.model();}
const HydroPhaseGeometry& HydroMaterialTransport::midpointGeometry() const noexcept
{return data_->volumeMiddle?*data_->volumeMiddle:data_->midpoint->geometry();}
const Interface::GeometricVofTransport& HydroMaterialTransport::vofTransport() const
{require(data_->vof.has_value(),"not a geometric VOF interval");return *data_->vof;}
const HydroAdvectorReconstruction& HydroMaterialTransport::fluxReconstruction() const
{return data_->volumeReconstruction?*data_->volumeReconstruction:midpointRate().spatialAdvector().reconstruction();}
bool HydroMaterialTransport::isSpaceTime() const noexcept{return data_->motion.has_value();}
const HydroCharacteristicSource& HydroMaterialTransport::characteristicSource() const
{require(isSpaceTime(),"not a FirstJet spacetime interval");return *data_->motion;}
const HydroSpaceTimeReport& HydroMaterialTransport::spaceTimeReport() const
{require(isSpaceTime(),"not a FirstJet spacetime interval");return data_->spaceTime;}
const std::vector<HydroPhaseMeasure>& HydroMaterialTransport::phaseAtomIntegrals() const
{require(isSpaceTime()||data_->vof.has_value(),"not an integrated phase interval");return data_->phaseIntegrals;}
const std::vector<long double>& HydroMaterialTransport::atomIntegrals() const noexcept{return data_->atomIntegrals;}
const std::vector<long double>& HydroMaterialTransport::dualFaceIntegrals() const noexcept{return data_->dualIntegrals;}
const HydroPhaseGeometry& HydroMaterialTransport::beginGeometry() const noexcept{return data_->begin;}
const HydroPhaseRate& HydroMaterialTransport::midpointRate() const
{require(data_->midpoint.has_value(),"geometric VOF has no FirstJet midpoint rate");return *data_->midpoint;}
const HydroPhaseGeometry& HydroMaterialTransport::endGeometry() const noexcept{return data_->end;}
const HydroDofLayout& HydroMaterialTransport::layout() const noexcept{return midpointGeometry().layout();}
HydroMaterialSources HydroMaterialTransport::sources() const
{
    if(data_->vof){HydroMaterialSources s;s.model=model();s.vofIds={beginGeometry().vofIdentity(),midpointGeometry().vofIdentity(),endGeometry().vofIdentity()};
        s.vofStamp=data_->vof;return s;}
    return {beginGeometry().geometryIdentity(),midpointRate().geometry().geometryIdentity(),endGeometry().geometryIdentity(),midpointRate().advectorIdentity()};
}
LAMP::TimeLayer HydroMaterialTransport::timeLayer() const{return midpointGeometry().timeLayer();}
double HydroMaterialTransport::deltaTime() const noexcept{return data_->deltaTime;}
double HydroMaterialTransport::liquidDensity() const noexcept{return data_->liquid;}
double HydroMaterialTransport::gasDensity() const noexcept{return data_->gas;}
bool HydroMaterialTransport::isSpatial() const noexcept{return data_->vof.has_value()||data_->midpoint->isSpatial();}
const HydroMomentumGeometry& HydroMaterialTransport::momentumGeometry() const
{return data_->volumeReconstruction?data_->volumeReconstruction->momentumGeometry():midpointRate().momentumGeometry();}
const HydroMaterialFlux& HydroMaterialTransport::midpointFlux() const
{require(data_->spatial.has_value(),"geometric VOF has no FirstJet midpoint flux");require(isSpatial(),"requires spatial material owner");return *data_->spatial;}
const std::vector<long double>& HydroMaterialTransport::atomRates() const
{return isSpaceTime()||data_->vof?data_->atoms:data_->spatial->atomRates();}
const std::vector<long double>& HydroMaterialTransport::dualFaceRates() const
{return isSpaceTime()||data_->vof?data_->dual:data_->spatial->dualFaceRates();}
const std::vector<long double>& HydroMaterialTransport::faceRates() const noexcept
{return isSpaceTime()||data_->vof?data_->faces:data_->spatial->faceRates();}
const std::vector<long double>& HydroMaterialTransport::faceIntegrals() const noexcept{return data_->integrals;}
const std::vector<HydroMaterialCell>& HydroMaterialTransport::cells() const noexcept{return data_->cells;}
const HydroMaterialReport& HydroMaterialTransport::report() const noexcept{return data_->report;}
void HydroMaterialTransport::validateFor(const QuadtreeMesh& mesh,HydroMaterialSources source,LAMP::TimeLayer clock) const
{
    require(source.model==model(),"different material geometry model");
    if(data_->vof){
        require(source.begin.value()==0&&source.midpoint.value()==0&&source.end.value()==0&&source.advector.value()==0,"VOF source carries a FirstJet identity");
        require(source.vofStamp&&data_->vof->sameInstance(*source.vofStamp),"different accepted VOF transport instance");
        beginGeometry().validateVofFor(mesh,source.vofIds[0],beginGeometry().timeLayer());
        const auto& checkedDomain=beginGeometry().vofPublication().domain();
        validateVofPhaseAfterMeshCheck(mesh,midpointGeometry(),source.vofIds[1],clock,checkedDomain);
        auto endClock=clock;endClock.position=LAMP::TimeLayer::Position::End;
        validateVofPhaseAfterMeshCheck(mesh,endGeometry(),source.vofIds[2],endClock,checkedDomain);
        data_->vof->validateFor(mesh,beginGeometry().vofPublication(),data_->vof->advector(),liquidDensity(),gasDensity());
        // The factories certify this immutable phase/momentum mesh binding.
        // Keep the original live check for any different binding.
        if(momentumGeometry().layout().domain().meshIdentity()!=checkedDomain.meshIdentity())
            momentumGeometry().validateMesh(mesh);
        return;
    }
    require(!source.vofStamp&&source.vofIds[0].value()==0&&source.vofIds[1].value()==0&&source.vofIds[2].value()==0,"FirstJet source carries a VOF identity");
    if(data_->motion)data_->motion->validateFor(mesh,data_->motion->previous().identity(),data_->motion->arrival().identity());
    beginGeometry().validateFor(mesh,source.begin,beginGeometry().timeLayer());
    midpointRate().geometry().validateFor(mesh,source.midpoint,clock);
    require(source.advector==sources().advector,"different actual advector identity");
    if(isSpatial())midpointRate().spatialAdvector().validateFor(mesh,momentumGeometry(),source.midpoint,clock);
    clock.position=LAMP::TimeLayer::Position::End;
    endGeometry().validateFor(mesh,source.end,clock);
}
void HydroMaterialTransport::validateFor(const QuadtreeMesh& mesh,const HydroMomentumGeometry& expected,
    HydroMaterialSources source,LAMP::TimeLayer clock) const
{
    require(isSpatial(),"requires spatial material owner");
    if(data_->vof){expected.validateMesh(mesh);require(momentumGeometry().sameIdentity(expected),"different volume-method momentum geometry instance");}
    else midpointRate().spatialAdvector().validateFor(mesh,expected,source.midpoint,clock);
    validateFor(mesh,source,clock);
}
}

