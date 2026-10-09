#include "Modules/Hydrodynamics/CchMaterialTransport.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>
#include <functional>
#include <stdexcept>
#include <string>
#include <utility>

namespace Quantum::Hydrodynamics
{
namespace
{
void require(bool ok,const std::string& why)
{if(!ok)throw std::invalid_argument("[CchMaterialTransport] "+why);}
long double finite(long double value,const std::string& why)
{require(std::isfinite(value),"nonfinite "+why);return value;}
long double densitySum(CchPhaseMeasure measure,double liquid,double gas,const std::string& why)
{
    const auto l=finite(liquid*measure.liquid,why+" liquid component");
    const auto g=finite(gas*measure.gas,why+" gas component");
    require(measure.liquid==0 || l!=0,"unrepresentable "+why+" liquid component");
    require(measure.gas==0 || g!=0,"unrepresentable "+why+" gas component");
    return finite(l+g,why);
}
long double atomSum(const std::vector<long double>& values,const std::vector<int>& indices)
{
    long double sum=0,correction=0;
    for(int i:indices){const auto y=values[i]-correction,t=sum+y;correction=(t-sum)-y;sum=t;}
    return finite(sum,"spatial material atom sum");
}
}
struct CchMaterialFlux::Data
{
    Data(CchPhaseRate p,double l,double g):phase(std::move(p)),liquid(l),gas(g){}
    CchPhaseRate phase;
    double liquid,gas;
    std::vector<long double> atoms,faces,dual;
};
CchMaterialFlux::CchMaterialFlux(std::shared_ptr<const Data> d):data_(std::move(d)){}
CchMaterialFlux CchMaterialFlux::prepare(const QuadtreeMesh& mesh,CchPhaseRate phase,double liquid,double gas)
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
    return CchMaterialFlux(std::move(data));
}
const CchPhaseRate& CchMaterialFlux::phaseRate() const noexcept{return data_->phase;}
const CchMomentumGeometry& CchMaterialFlux::momentumGeometry() const{return phaseRate().momentumGeometry();}
const CchDofLayout& CchMaterialFlux::layout() const noexcept{return phaseRate().geometry().layout();}
CchAdvectorIdentity CchMaterialFlux::advectorIdentity() const{return phaseRate().advectorIdentity();}
Interface::HermiteGalsIdentity CchMaterialFlux::geometryIdentity() const{return phaseRate().geometry().geometryIdentity();}
Interface::HermiteGalsTimeLayer CchMaterialFlux::timeLayer() const{return phaseRate().geometry().timeLayer();}
double CchMaterialFlux::liquidDensity() const noexcept{return data_->liquid;}
double CchMaterialFlux::gasDensity() const noexcept{return data_->gas;}
const std::vector<long double>& CchMaterialFlux::atomRates() const noexcept{return data_->atoms;}
const std::vector<long double>& CchMaterialFlux::faceRates() const noexcept{return data_->faces;}
const std::vector<long double>& CchMaterialFlux::dualFaceRates() const noexcept{return data_->dual;}
Vector2D CchMaterialFlux::physicalFlux(int pressure,Vector2D point,bool liquid) const
{
    const auto velocity=phaseRate().spatialAdvector().cellPhysicalVelocityJet(pressure,point,liquid).jet.value;
    const auto flux=velocity*(liquid?liquidDensity():gasDensity());
    require(std::isfinite(flux.x)&&std::isfinite(flux.y),"nonfinite physical material flux");
    return flux;
}
void CchMaterialFlux::validateFor(const QuadtreeMesh& mesh,const CchMomentumGeometry& expected,
    Interface::HermiteGalsIdentity geometry,CchAdvectorIdentity advector,Interface::HermiteGalsTimeLayer clock) const
{
    phaseRate().geometry().validateFor(mesh,geometry,clock);
    phaseRate().spatialAdvector().validateFor(mesh,expected,geometry,clock);
    require(advector==advectorIdentity(),"different actual advector identity");
}
struct CchMaterialTransport::Data
{
    Data(CchPhaseGeometry a,CchPhaseRate m,CchPhaseGeometry b,double l,double g,double dt,const CchMaterialFlux* flux)
        :begin(std::move(a)),midpoint(std::move(m)),end(std::move(b)),liquid(l),gas(g),deltaTime(dt),
         spatial(*flux){}
    CchPhaseGeometry begin;
    CchPhaseRate midpoint;
    CchPhaseGeometry end;
    double liquid,gas,deltaTime;
    CchMaterialFlux spatial;
    std::vector<long double> integrals,atoms,atomIntegrals,faces,dual,dualIntegrals;
    std::vector<CchPhaseMeasure> phaseIntegrals;
    std::optional<CchCharacteristicSource> motion;
    CchSpaceTimeReport spaceTime;
    std::vector<CchMaterialCell> cells;
    CchMaterialReport report;
};
CchMaterialTransport::CchMaterialTransport(std::shared_ptr<const Data> d):data_(std::move(d)){}
CchMaterialTransport CchMaterialTransport::prepare(const QuadtreeMesh& mesh,CchPhaseGeometry begin,
    CchPhaseRate midpoint,CchPhaseGeometry end,double liquid,double gas)
{
    return prepare(mesh,std::move(begin),CchMaterialFlux::prepare(mesh,midpoint,liquid,gas),std::move(end));
}
CchMaterialTransport CchMaterialTransport::prepare(const QuadtreeMesh& mesh,CchPhaseGeometry begin,
    CchMaterialFlux midpoint,CchPhaseGeometry end)
{
    return assemble(mesh,std::move(begin),midpoint.phaseRate(),std::move(end),midpoint.liquidDensity(),midpoint.gasDensity(),&midpoint);
}
CchMaterialTransport CchMaterialTransport::prepareSpaceTime(const QuadtreeMesh& mesh,
    Interface::HermiteGalsPublication original,CchCharacteristicSource motion,double liquid,double gas,
    CchSpaceTimeConfig config,const Interface::HermiteGalsBoundaryFirstJetQuery& boundary,
    const Interface::HermiteGalsPublicationConfig& publication)
{
    using namespace Interface;using Position=HermiteGalsTimeLayer::Position;
    original.validateMesh(mesh);motion.validateFor(mesh,motion.previous().identity(),motion.arrival().identity());
    const auto& reconstruction=motion.geometry().reconstruction();const auto& geometry=reconstruction.momentumGeometry();
    motion.previous().validateFor(mesh,geometry,original.state().identity(),original.state().timeLayer());
    require(original.state().primary().transport==HermiteGalsTransport::FirstJet&&original.state().primary().phaseVolumes.empty(),"spacetime requires independent Pure FirstJet");
    require(std::isfinite(liquid)&&liquid>0&&std::isfinite(gas)&&gas>0,"nonpositive/nonfinite density");
    require(std::isfinite(config.relativeTolerance)&&config.relativeTolerance>0&&config.relativeTolerance<1
        &&config.initialPanels>0&&config.maxDepth>=0&&config.maxDepth<=50&&config.maxGeometrySamples>2,"invalid spacetime quadrature budget");
    const auto parent=motion.arrival().timeLayer();const double begin=parent.begin,end=parent.end,dt=end-begin,middle=begin+dt*.5;
    CchSpaceTimeReport report;
    const auto sampleGeometry=[&](double time) {
        require(++report.geometrySamples<=static_cast<std::size_t>(config.maxGeometrySamples),"spacetime geometry sample budget exhausted");
        auto clock=parent;if(time==middle)clock.position=Position::Midpoint;else clock.end=time;
        const HermiteGalsResolvedMotion resolved(original.state().identity(),original.state().timeLayer(),motion.arrival().geometryIdentity(),clock,
            [motion,clock,config,begin](Vector2D point){
                const auto t=motion.traceBackward(point,clock.time(),begin,config.characteristic);
                std::optional<Boundary> side;
                if(t.boundary()){const auto b=*t.boundary();side=b.normalAxis==0?(b.outwardSign<0?Boundary::west:Boundary::east):(b.outwardSign<0?Boundary::south:Boundary::north);}
                return LevelSetBoundaryTrace::JetResult{{t.point(),clock.time()-t.time(),side,t.boundary().has_value(),false},t.jacobian(),t.substeps()};
            });
        const auto step=prepareHermiteGalsFirstJetStep(mesh,original,resolved,boundary,publication,
            boundary?HermiteGalsBoundaryTreatment::PrescribedOuterFirstJet:HermiteGalsBoundaryTreatment::InflowHistory);
        return CchPhaseGeometry::build(mesh,geometry.metrics(),step.publication());
    };
    const auto initial=CchPhaseGeometry::build(mesh,geometry.metrics(),original);
    const auto midGeometry=sampleGeometry(middle),endGeometry=sampleGeometry(end);
    const auto rateAt=[&](const CchPhaseGeometry& g){return CchPhaseRate::prepareSpatial(mesh,g,
        reconstruction.interpolate(mesh,motion.previous(),motion.arrival(),g.geometryIdentity(),g.timeLayer()));};
    const auto midpoint=rateAt(midGeometry);
    const std::size_t atoms=reconstruction.atoms().size();
    const auto evaluate=[&](double time){return time==middle?midpoint.atomRates():rateAt(sampleGeometry(time)).atomRates();};
    std::vector<long double> bound(atoms);
    for(std::size_t i=0;i<atoms;++i){long double speed=0;
        for(const auto* a:{&motion.previous(),&motion.arrival()}){const auto& p=a->atomPolynomials()[i];long double work=0,power=1;
            for(auto c:p.beginCoefficients){work+=std::abs(c)*power;power*=p.length;}speed=std::max(speed,work);}
        bound[i]=finite(dt*reconstruction.atoms()[i].area*speed,"spacetime atom work scale");
    }
    struct Panel {double lo,hi;int depth;std::vector<CchPhaseMeasure> value,error;};
    const auto panel=[&](double lo,double hi,int depth){
        require(hi>lo,"unresolved spacetime panel");const double width=hi-lo,center=lo+width*.5;
        require(center>lo&&center<hi,"unresolved spacetime center");const auto mid=evaluate(center);
        std::vector<CchPhaseMeasure> coarse(atoms),fine(atoms);
        const auto add=[&](auto& sum,const auto& value,long double weight){for(std::size_t i=0;i<atoms;++i){sum[i].liquid+=weight*value[i].liquid;sum[i].gas+=weight*value[i].gas;}};
        add(coarse,mid,width*(4.L/9));add(fine,mid,width*(64.L/225));
        for(double z:{-.77459666924148337704,.77459666924148337704})add(coarse,evaluate(center+width*.5*z),width*(5.L/18));
        constexpr double nodes[4]={-.90617984593866399280,-.53846931010568309104,.53846931010568309104,.90617984593866399280};
        constexpr long double weights[4]={.11846344252809454376L,.23931433524968323402L,.23931433524968323402L,.11846344252809454376L};
        for(int q=0;q<4;++q)add(fine,evaluate(center+width*.5*nodes[q]),width*weights[q]);
        std::vector<CchPhaseMeasure> error(atoms);for(std::size_t i=0;i<atoms;++i)error[i]={std::abs(fine[i].liquid-coarse[i].liquid),std::abs(fine[i].gas-coarse[i].gas)};
        return Panel{lo,hi,depth,std::move(fine),std::move(error)};
    };
    std::vector<Panel> panels;
    for(int i=0;i<config.initialPanels;++i)panels.push_back(panel(begin+dt*(static_cast<double>(i)/config.initialPanels),i+1==config.initialPanels?end:begin+dt*(static_cast<double>(i+1)/config.initialPanels),0));
    std::vector<CchPhaseMeasure> integrated(atoms);
    for(;;){
        std::vector<CchPhaseMeasure> error(atoms);integrated.assign(atoms,{});std::size_t worst=0;long double worstScore=-1;
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
        if(relative<=config.relativeTolerance){report.absolutePhaseErrorEstimate=absolute;report.maxRelativePhaseErrorEstimate=relative;break;}
        const auto lo=panels[worst].lo,hi=panels[worst].hi;const int depth=panels[worst].depth+1;
        require(depth<=config.maxDepth,"spacetime refinement depth exhausted");report.deepestPanel=std::max(report.deepestPanel,depth);
        const double center=lo+(hi-lo)*.5;auto left=panel(lo,center,depth),right=panel(center,hi,depth);
        panels[worst]=std::move(left);panels.push_back(std::move(right));
    }
    report.acceptedPanels=panels.size();
    const auto flux=CchMaterialFlux::prepare(mesh,midpoint,liquid,gas);
    return assemble(mesh,initial,midpoint,endGeometry,liquid,gas,&flux,&integrated,&motion,&report);
}
CchMaterialTransport CchMaterialTransport::assemble(const QuadtreeMesh& mesh,CchPhaseGeometry begin,
    CchPhaseRate midpoint,CchPhaseGeometry end,double liquid,double gas,const CchMaterialFlux* flux,
    const std::vector<CchPhaseMeasure>* integrated,const CchCharacteristicSource* motion,const CchSpaceTimeReport* report)
{
    using Position=Interface::HermiteGalsTimeLayer::Position;
    const auto& mid=midpoint.geometry();
    for(const auto* geometry:std::array<const CchPhaseGeometry*,3>{{&begin,&mid,&end}})
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
    const auto& layout=data->midpoint.geometry().layout();
    if(integrated) {
        require(motion&&report&&integrated->size()==data->spatial.atomRates().size(),"invalid spacetime material provenance");
        data->motion=*motion;data->spaceTime=*report;data->phaseIntegrals=*integrated;
        for(const auto& q:*integrated){const auto f=densitySum(q,liquid,gas,"spacetime atom material integral");
            data->atomIntegrals.push_back(f);data->atoms.push_back(finite(f/dt,"average atom material rate"));}
        const auto& r=data->spatial.phaseRate().spatialAdvector().reconstruction();
        for(const auto& ids:r.primalAtoms()){const auto f=atomSum(data->atomIntegrals,ids);data->integrals.push_back(f);data->faces.push_back(f/dt);}
        for(const auto& ids:r.dualAtoms()){const auto f=atomSum(data->atomIntegrals,ids);data->dualIntegrals.push_back(f);data->dual.push_back(f/dt);}
    } else {
        data->atoms=data->spatial.atomRates();data->faces=data->spatial.faceRates();data->dual=data->spatial.dualFaceRates();
        const auto integrate=[&](const auto& rates,auto& integrals,const std::string& kind){for(auto rate:rates){
            const auto f=finite(dt*rate,kind+" material integral");require(rate==0||f!=0,"unrepresentable "+kind+" material integral");integrals.push_back(f);}};
        // Preserve the original primal failure ordering and diagnostic.
        integrate(data->faces,data->integrals,"face");integrate(data->atoms,data->atomIntegrals,"atom");integrate(data->dual,data->dualIntegrals,"dual");
    }
    data->cells.resize(layout.pressureCount());
    for(int p=0;p<layout.pressureCount();++p){
        auto& cell=data->cells[p];
        const auto va=data->begin.cellVolumes()[p],vm=data->midpoint.geometry().cellVolumes()[p],vb=data->end.cellVolumes()[p];
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
    const auto& di=data->midpoint.geometry().metrics().divergence().matrix();
    for(Eigen::Index p=0;p<di.outerSize();++p){
        auto& cell=data->cells[static_cast<std::size_t>(p)];long double gross=std::abs(cell.massChange);
        for(CchIntegralDivergence::Matrix::InnerIterator it(di,p);it;++it){
            const int f=static_cast<int>(it.col());
            const double area=data->midpoint.geometry().metrics().faceMeasures()[f];
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
            require(data->midpoint.geometry().metrics().faceMeasures()[f]==0 && data->integrals[f]==0,
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
    return CchMaterialTransport(std::move(data));
}
bool CchMaterialTransport::isSpaceTime() const noexcept{return data_->motion.has_value();}
const CchCharacteristicSource& CchMaterialTransport::characteristicSource() const
{require(isSpaceTime(),"not a spacetime interval");return *data_->motion;}
const CchSpaceTimeReport& CchMaterialTransport::spaceTimeReport() const
{require(isSpaceTime(),"not a spacetime interval");return data_->spaceTime;}
const std::vector<CchPhaseMeasure>& CchMaterialTransport::phaseAtomIntegrals() const
{require(isSpaceTime(),"not a spacetime interval");return data_->phaseIntegrals;}
const std::vector<long double>& CchMaterialTransport::atomIntegrals() const noexcept{return data_->atomIntegrals;}
const std::vector<long double>& CchMaterialTransport::dualFaceIntegrals() const noexcept{return data_->dualIntegrals;}
const CchPhaseGeometry& CchMaterialTransport::beginGeometry() const noexcept{return data_->begin;}
const CchPhaseRate& CchMaterialTransport::midpointRate() const noexcept{return data_->midpoint;}
const CchPhaseGeometry& CchMaterialTransport::endGeometry() const noexcept{return data_->end;}
const CchDofLayout& CchMaterialTransport::layout() const noexcept{return midpointRate().geometry().layout();}
CchMaterialSources CchMaterialTransport::sources() const{return {beginGeometry().geometryIdentity(),midpointRate().geometry().geometryIdentity(),endGeometry().geometryIdentity(),
    midpointRate().advectorIdentity()};}
Interface::HermiteGalsTimeLayer CchMaterialTransport::timeLayer() const{return midpointRate().geometry().timeLayer();}
double CchMaterialTransport::deltaTime() const noexcept{return data_->deltaTime;}
double CchMaterialTransport::liquidDensity() const noexcept{return data_->liquid;}
double CchMaterialTransport::gasDensity() const noexcept{return data_->gas;}
bool CchMaterialTransport::isSpatial() const noexcept{return midpointRate().isSpatial();}
const CchMomentumGeometry& CchMaterialTransport::momentumGeometry() const{return midpointRate().momentumGeometry();}
const CchMaterialFlux& CchMaterialTransport::midpointFlux() const
{require(isSpatial(),"requires spatial material owner");return data_->spatial;}
const std::vector<long double>& CchMaterialTransport::atomRates() const
{return isSpaceTime()?data_->atoms:data_->spatial.atomRates();}
const std::vector<long double>& CchMaterialTransport::dualFaceRates() const
{return isSpaceTime()?data_->dual:data_->spatial.dualFaceRates();}
const std::vector<long double>& CchMaterialTransport::faceRates() const noexcept{return isSpaceTime()?data_->faces:data_->spatial.faceRates();}
const std::vector<long double>& CchMaterialTransport::faceIntegrals() const noexcept{return data_->integrals;}
const std::vector<CchMaterialCell>& CchMaterialTransport::cells() const noexcept{return data_->cells;}
const CchMaterialReport& CchMaterialTransport::report() const noexcept{return data_->report;}
void CchMaterialTransport::validateFor(const QuadtreeMesh& mesh,CchMaterialSources source,Interface::HermiteGalsTimeLayer clock) const
{
    if(data_->motion)data_->motion->validateFor(mesh,data_->motion->previous().identity(),data_->motion->arrival().identity());
    beginGeometry().validateFor(mesh,source.begin,beginGeometry().timeLayer());
    midpointRate().geometry().validateFor(mesh,source.midpoint,clock);
    require(source.advector==sources().advector,"different actual advector identity");
    if(isSpatial())midpointRate().spatialAdvector().validateFor(mesh,momentumGeometry(),source.midpoint,clock);
    clock.position=Interface::HermiteGalsTimeLayer::Position::End;
    endGeometry().validateFor(mesh,source.end,clock);
}
void CchMaterialTransport::validateFor(const QuadtreeMesh& mesh,const CchMomentumGeometry& expected,
    CchMaterialSources source,Interface::HermiteGalsTimeLayer clock) const
{
    require(isSpatial(),"requires spatial material owner");
    midpointRate().spatialAdvector().validateFor(mesh,expected,source.midpoint,clock);
    validateFor(mesh,source,clock);
}
}
