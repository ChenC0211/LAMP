#include "Modules/Hydrodynamics/HydroDualMaterial.h"
#include "Modules/Hydrodynamics/Detail/HydroFiniteContraction.h"
#include "Modules/Hydrodynamics/HydroSpatialAdvector.h"

#include "Core/Mesh/QuadtreeMesh.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>

namespace LAMP::Hydrodynamics
{
namespace
{
constexpr long double tolerance=256.L*std::numeric_limits<double>::epsilon();
void require(bool ok,const char* why)
{ if(!ok)throw std::invalid_argument(std::string("[HydroDualMaterial] ")+why); }
long double finite(long double v,const char* why)
{ require(std::isfinite(v),why);return v; }
long double defect(long double difference,long double scale,const char* why)
{
    finite(difference,why);finite(scale,why);require(scale>=0,why);
    require(scale>0 || difference==0,why);
    return scale>0?std::abs(difference)/scale:0;
}
// Backward error of one integral balance. Positive stored endpoint masses
// and the actual exchanged mass define its scale; no absolute floor or dt.
long double continuityRatio(long double residual,long double m0,long double m1,long double gross)
{
    const auto scale=std::max({m0,m1,gross});
    require(scale>0,"zero conservative endpoint mass scale");
    return std::abs(residual/scale)/(m0/scale+m1/scale+gross/scale);
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
long double change(HydroPhaseMeasure a,HydroPhaseMeasure b,double l,double g)
{
    const auto dl=b.liquid-a.liquid,dg=b.gas-a.gas;
    const auto ml=finite(l*dl,"nonfinite dual liquid mass change"),mg=finite(g*dg,"nonfinite dual gas mass change");
    require((dl==0 || ml!=0) && (dg==0 || mg!=0),"unrepresentable dual mass change");
    return finite(ml+mg,"nonfinite dual mass change");
}
}
struct HydroDualMaterial::Data
{
    Data(HydroMaterialTransport t,HydroMomentumMass a,HydroMomentumMass m,HydroMomentumMass b)
        :transport(std::move(t)),begin(std::move(a)),midpoint(std::move(m)),end(std::move(b)){}
    HydroMaterialTransport transport;
    HydroMomentumMass begin,midpoint,end;
    std::vector<HydroDualMaterialCell> cells;
    HydroDualMaterialReport report;
    bool conservativeVofMomentum=false;
};
HydroDualMaterial::HydroDualMaterial(std::shared_ptr<const Data> d):data_(std::move(d)){}
HydroDualMaterial HydroDualMaterial::prepare(const QuadtreeMesh& mesh,HydroMaterialTransport material,
    HydroMomentumMass begin,HydroMomentumMass midpoint,HydroMomentumMass end)
{
    material.validateFor(mesh,material.sources(),material.timeLayer());
    const auto source=material.sources();
    if(material.model()==HydroPhaseGeometry::Model::GeometricVof){
        const auto& checkedDomain=material.beginGeometry().vofPublication().domain();
        validateVofPhaseAfterMeshCheck(mesh,begin.phaseGeometry(),source.vofIds[0],material.beginGeometry().timeLayer(),checkedDomain);
        validateVofPhaseAfterMeshCheck(mesh,midpoint.phaseGeometry(),source.vofIds[1],material.timeLayer(),checkedDomain);
        validateVofPhaseAfterMeshCheck(mesh,end.phaseGeometry(),source.vofIds[2],material.endGeometry().timeLayer(),checkedDomain);
    }else{
        begin.validateFor(mesh,source.begin,material.beginGeometry().timeLayer());
        midpoint.validateFor(mesh,source.midpoint,material.timeLayer());
        end.validateFor(mesh,source.end,material.endGeometry().timeLayer());
    }
    require(begin.geometry().sameIdentity(midpoint.geometry()) && end.geometry().sameIdentity(midpoint.geometry()),"different dual geometry identity");
    for(const auto* mass:std::array<const HydroMomentumMass*,3>{{&begin,&midpoint,&end}}){
        require(mass->model()==material.model(),"different material mass geometry model");
        require(mass->representation()==midpoint.representation(),"different dual density representation");
        require(mass->layout().sameIdentity(material.layout()),"different native layout");
        require(&mass->geometry().metrics().divergence().matrix()==&material.midpointGeometry().metrics().divergence().matrix(),"different retained D_I");
        require(mass->liquidDensity()==material.liquidDensity() && mass->gasDensity()==material.gasDensity(),"different material densities");
    }
    // The entry validation already checked the same material/source/clock.
    // Exact geometry-owner equality also retains its checked mesh binding;
    // any foreign owner takes the original complete validator and diagnostics.
    if(material.model()!=HydroPhaseGeometry::Model::GeometricVof
        ||!material.momentumGeometry().sameIdentity(midpoint.geometry()))
        material.validateFor(mesh,midpoint.geometry(),source,material.timeLayer());
    auto d=std::make_shared<Data>(std::move(material),std::move(begin),std::move(midpoint),std::move(end));
    const auto& geometry=d->midpoint.geometry();const auto& layout=geometry.layout();const double dt=d->transport.deltaTime();
    d->conservativeVofMomentum=d->transport.model()==HydroPhaseGeometry::Model::GeometricVof
        &&d->transport.beginGeometry().vofPublication().geometryKind()==Interface::GeometricVofGeometryKind::LocalPlic
        &&d->midpoint.representation()==HydroMomentumMass::DensityRepresentation::CellAverageProjection;
    const auto& rates=d->transport.dualFaceRates();
    require(rates.size()==geometry.segments().size(),"different original dual flux view");
    const auto& integrals=d->transport.dualFaceIntegrals();
    require(integrals.size()==rates.size(),"different original dual integral view");
    for(auto integral:integrals)finite(integral,"nonfinite dual material integral");
    std::vector<long double> primalRate(layout.pressureCount()),primalGross(layout.pressureCount()),primalGrossIntegral(layout.pressureCount());
    const auto& di=d->transport.midpointGeometry().metrics().divergence().matrix();
    for(Eigen::Index p=0;p<di.outerSize();++p)for(HydroIntegralDivergence::Matrix::InnerIterator it(di,p);it;++it){
        const auto rate=d->transport.faceRates()[static_cast<std::size_t>(it.col())];
        primalRate[p]+=(it.value()>0?1:-1)*rate;primalGross[p]+=std::abs(rate);
        primalGrossIntegral[p]+=std::abs(d->transport.faceIntegrals()[static_cast<std::size_t>(it.col())]);
    }
    if(d->conservativeVofMomentum)for(int p=0;p<layout.pressureCount();++p){
        const auto& cell=d->transport.cells()[p];
        const auto ratio=continuityRatio(cell.integratedDefect,cell.massBegin,cell.massEnd,primalGrossIntegral[p]);
        d->report.maxRelativePrimalContinuityDefect=std::max(d->report.maxRelativePrimalContinuityDefect,ratio);
        if(ratio>tolerance){std::ostringstream why;why.precision(std::numeric_limits<long double>::max_digits10);
            why<<"[HydroDualMaterial] conservative VOF primal continuity exceeds roundoff native_cell="<<layout.nativeCell(p)
                <<" residual="<<ratio<<" limit="<<tolerance;
            throw std::invalid_argument(why.str());}
    }
    d->cells.resize(layout.velocityCount());std::array<long double,2> balanceScale{};
    for(int u=0;u<layout.velocityCount();++u){
        const int component=layout.normalAxis(u)==HydroDofLayout::NormalAxis::First?0:1;
        auto& cell=d->cells[u];auto& report=d->report.components[component];
        if(d->transport.model()==HydroPhaseGeometry::Model::GeometricVof&&
            d->midpoint.representation()==HydroMomentumMass::DensityRepresentation::CellAverageProjection){
            const auto& before=d->begin.phaseGeometry().vofPublication().plicCells();
            const auto& after=d->end.phaseGeometry().vofPublication().plicCells();
            long double correction=0;
            const auto& pieces=geometry.volumes()[u].overlaps;
            for(std::size_t k=0;k<pieces.size();++k){const auto& overlap=pieces[k];
                const auto native=layout.nativeCell(overlap.pressure);
                const auto& a=before[native].inventory;const auto& b=after[native].inventory;
                const auto affected=[&](int phase){return d->begin.usesFiniteP0Partition(overlap.pressure,component,phase)
                    ||d->end.usesFiniteP0Partition(overlap.pressure,component,phase);};
                const bool finiteLiquid=affected(0),finiteGas=affected(1);long double term=0;
                if(!finiteLiquid&&!finiteGas){
                    // Preserve the original delta-first arithmetic/order for
                    // every ordinary source/component, without subtracting M.
                    const auto delta=change({a[0],a[1]},{b[0],b[1]},d->transport.liquidDensity(),d->transport.gasDensity());
                    term=finite(overlap.cellFraction*delta,"nonfinite canonical VOF dual mass change");
                    require(delta==0||term!=0,"unrepresentable canonical VOF dual mass change");
                }else{
                    const auto& qa=d->begin.representedOverlapPhaseVolumes()[u][k];
                    const auto& qb=d->end.representedOverlapPhaseVolumes()[u][k];
                    const std::array<long double,2> projectedA{{qa.liquid,qa.gas}},projectedB{{qb.liquid,qb.gas}};
                    // RN(S*F) is nonlinear in S: a flagged phase MUST use the
                    // difference of the actual endpoint table. Never apply a
                    // fresh CDF to delta stock or difference rounded total M.
                    for(int phase=0;phase<2;++phase){const bool useEndpoint=phase==0?finiteLiquid:finiteGas;
                        const auto delta=useEndpoint?projectedB[phase]-projectedA[phase]:b[phase]-a[phase];
                        const auto rho=phase==0?d->transport.liquidDensity():d->transport.gasDensity();
                        const auto product=finite(rho*delta,"nonfinite finite-P0 canonical phase mass change");
                        require(delta==0||product!=0,"unrepresentable finite-P0 canonical phase mass change");
                        const auto mapped=useEndpoint?product:finite(overlap.cellFraction*product,"nonfinite ordinary canonical phase mapping");
                        require(product==0||mapped!=0,"unrepresentable ordinary canonical phase mapping");
                        term=finite(term+mapped,"nonfinite finite-P0 canonical mass sum");
                    }
                }
                const auto adjusted=term-correction,next=cell.massChange+adjusted;
                correction=(next-cell.massChange)-adjusted;cell.massChange=next;
            }
            finite(cell.massChange,"nonfinite canonical VOF dual mass sum");
        }else cell.massChange=change(d->begin.representedPhaseVolumes()[u],d->end.representedPhaseVolumes()[u],
            d->transport.liquidDensity(),d->transport.gasDensity());
        long double outwardRate=0,grossRate=0,grossIntegral=0,mappedRate=0,mappedGross=0,mappedIntegral=0;
        for(int side=0;side<4;++side)for(int id:geometry.volumes()[u].sides[side]){
            const int sign=side%2==0?-1:1;
            outwardRate+=sign*rates[id];cell.outwardIntegral+=sign*integrals[id];
            grossRate+=std::abs(rates[id]);grossIntegral+=std::abs(integrals[id]);
        }
        for(const auto& overlap:geometry.volumes()[u].overlaps){
            const auto& primal=d->transport.cells()[overlap.pressure];const auto w=overlap.cellFraction;
            mappedRate+=w*primalRate[overlap.pressure];mappedGross+=w*primalGross[overlap.pressure];mappedIntegral+=w*primal.outwardIntegral;
            cell.mappedPrimalMassChange+=w*primal.massChange;cell.mappedPrimalIntegratedDefect+=w*primal.integratedDefect;
        }
        // The retained time integrals define net interval rates. Keep the
        // original sum of separately divided face rates as raw diagnostics;
        // strong cancellation must not choose a different net R/KEP clock.
        cell.rawOutwardRate=finite(outwardRate,"nonfinite actual dual outward rate");
        cell.rawMappedPrimalOutwardRate=finite(mappedRate,"nonfinite mapped primal outward rate");
        cell.rawSpatialExchangeRate=finite(outwardRate-mappedRate,"nonfinite spatial exchange rate");
        cell.outwardRate=finite(cell.outwardIntegral/dt,"nonfinite canonical dual outward rate");
        cell.mappedPrimalOutwardRate=finite(mappedIntegral/dt,"nonfinite canonical mapped primal outward rate");
        cell.spatialExchangeRate=finite((cell.outwardIntegral-mappedIntegral)/dt,"nonfinite canonical spatial exchange rate");
        const auto exchange=defect(cell.spatialExchangeRate,grossRate+mappedGross,"invalid spatial exchange diagnostic");
        d->report.maxRelativeSpatialExchange=std::max(d->report.maxRelativeSpatialExchange,exchange);
        d->report.maxAbsSpatialExchange=std::max(d->report.maxAbsSpatialExchange,std::abs(cell.spatialExchangeRate));
        cell.geometryTransferChange=finite(cell.massChange-cell.mappedPrimalMassChange,"nonfinite geometry transfer change");
        cell.integratedDefect=finite(cell.massChange+cell.outwardIntegral,"nonfinite dual integrated defect");
        cell.defectRate=finite(cell.integratedDefect/dt,"nonfinite dual defect rate");
        const auto gross=finite(std::abs(cell.massChange)+grossIntegral,"nonfinite dual defect scale");
        cell.relativeDefect=defect(cell.integratedDefect,gross,"invalid relative dual defect");
        if(d->conservativeVofMomentum){
            // A backward-error scale uses the actual positive endpoint mass
            // and exchanged mass, not a nearly cancelling DeltaM or R/dt.
            // Normalize before adding so a finite endpoint near DBL_MAX does
            // not manufacture an overflowing scale. No absolute mass floor.
            const long double m0=d->begin.masses()[u],m1=d->end.masses()[u];
            const auto canonical=continuityRatio(cell.integratedDefect,m0,m1,grossIntegral);
            const auto storedDefect=finite((m1-m0)+cell.outwardIntegral,"nonfinite stored mass continuity defect");
            const auto stored=continuityRatio(storedDefect,m0,m1,grossIntegral);
            const auto massTransfer=continuityRatio(cell.geometryTransferChange,m0,m1,grossIntegral);
            const auto spatial=continuityRatio(cell.outwardIntegral-mappedIntegral,m0,m1,grossIntegral);
            d->report.maxRelativeMassTransferDefect=std::max(d->report.maxRelativeMassTransferDefect,massTransfer);
            d->report.maxRelativeFluxCommutatorDefect=std::max(d->report.maxRelativeFluxCommutatorDefect,spatial);
            d->report.maxRelativeContinuityDefect=std::max(d->report.maxRelativeContinuityDefect,canonical);
            d->report.maxRelativeStoredMassContinuityDefect=std::max(d->report.maxRelativeStoredMassContinuityDefect,stored);
            if(canonical>tolerance||stored>tolerance||massTransfer>tolerance||spatial>tolerance){
                std::ostringstream why;why.precision(std::numeric_limits<long double>::max_digits10);
                why<<"[HydroDualMaterial] conservative VOF continuity exceeds roundoff native_face="<<layout.nativeFace(u)
                    <<" canonical="<<canonical<<" stored="<<stored<<" limit="<<tolerance
                    <<" mapped_primal="<<cell.mappedPrimalIntegratedDefect<<" geometry_transfer="<<cell.geometryTransferChange
                    <<" spatial_exchange_integral="<<(cell.outwardIntegral-mappedIntegral);
                throw std::invalid_argument(why.str());
            }
        }
        const auto decomposition=defect(cell.integratedDefect-cell.mappedPrimalIntegratedDefect-cell.geometryTransferChange-dt*cell.spatialExchangeRate,
            gross+std::abs(cell.mappedPrimalMassChange)+dt*mappedGross,"invalid defect decomposition");
        d->report.maxDefectDecomposition=std::max(d->report.maxDefectDecomposition,decomposition);
        require(decomposition<=tolerance,"dual defect decomposition failed");
        report.geometryMass[0]+=d->begin.masses()[u];report.geometryMass[1]+=d->midpoint.masses()[u];report.geometryMass[2]+=d->end.masses()[u];
        report.massChange+=cell.massChange;report.integratedDefect+=cell.integratedDefect;report.l1IntegratedDefect+=std::abs(cell.integratedDefect);
        report.maxAbsDefectRate=std::max(report.maxAbsDefectRate,std::abs(cell.defectRate));report.maxRelativeCellDefect=std::max(report.maxRelativeCellDefect,cell.relativeDefect);
        balanceScale[component]+=gross;
    }
    for(std::size_t id=0;id<geometry.segments().size();++id){
        const auto& s=geometry.segments()[id];if(s.negativeVelocity>=0 && s.positiveVelocity>=0)continue;
        const int u=s.negativeVelocity>=0?s.negativeVelocity:s.positiveVelocity;
        const int component=layout.normalAxis(u)==HydroDofLayout::NormalAxis::First?0:1;
        d->report.components[component].boundaryIntegral+=(s.negativeVelocity>=0?1:-1)*integrals[id];
    }
    for(int component=0;component<2;++component){auto& report=d->report.components[component];
        for(auto v:report.geometryMass)finite(v,"nonfinite total dual mass");
        for(auto v:{report.massChange,report.boundaryIntegral,report.integratedDefect,report.l1IntegratedDefect})finite(v,"nonfinite total dual balance");
        report.balanceIdentity=defect(report.integratedDefect-report.massChange-report.boundaryIntegral,balanceScale[component],"dual global balance defect");
        require(report.balanceIdentity<=tolerance,"dual global balance identity failed");
    }
    return HydroDualMaterial(std::move(d));
}
const HydroMaterialTransport& HydroDualMaterial::material() const noexcept{return data_->transport;}
const HydroMomentumMass& HydroDualMaterial::beginMass() const noexcept{return data_->begin;}
const HydroMomentumMass& HydroDualMaterial::midpointMass() const noexcept{return data_->midpoint;}
const HydroMomentumMass& HydroDualMaterial::endMass() const noexcept{return data_->end;}
const HydroMomentumGeometry& HydroDualMaterial::geometry() const noexcept{return midpointMass().geometry();}
const HydroDofLayout& HydroDualMaterial::layout() const noexcept{return material().layout();}
HydroMaterialSources HydroDualMaterial::sources() const{return material().sources();}
LAMP::TimeLayer HydroDualMaterial::timeLayer() const{return material().timeLayer();}
double HydroDualMaterial::deltaTime() const noexcept{return material().deltaTime();}
const std::vector<long double>& HydroDualMaterial::faceRates() const noexcept{return material().dualFaceRates();}
const std::vector<long double>& HydroDualMaterial::faceIntegrals() const noexcept{return material().dualFaceIntegrals();}
const std::vector<HydroDualMaterialCell>& HydroDualMaterial::cells() const noexcept{return data_->cells;}
const HydroDualMaterialReport& HydroDualMaterial::report() const noexcept{return data_->report;}
bool HydroDualMaterial::usesConservativeVofMomentum() const noexcept{return data_->conservativeVofMomentum;}
void HydroDualMaterial::validateFor(const QuadtreeMesh& mesh,HydroMaterialSources source,LAMP::TimeLayer clock) const
{
    material().validateFor(mesh,source,clock);
    require(material().momentumGeometry().sameIdentity(geometry()),"different momentum geometry instance");
    // VOF material validation already checks this exact momentum geometry
    // against the current live mesh. Other methods keep their original scan.
    if(material().model()!=HydroPhaseGeometry::Model::GeometricVof)geometry().validateMesh(mesh);
}
std::vector<long double> HydroDualMaterial::consistencyImpulse(const QuadtreeMesh& mesh,const HydroNativeValues& u,
    HydroMaterialSources source,LAMP::TimeLayer clock) const
{
    validateFor(mesh,source,clock);u.validateFor(layout(),mesh);
    return consistencyImpulseAfterValidation(u);
}
std::vector<long double> HydroDualMaterial::consistencyImpulseAfterValidation(const HydroNativeValues& u) const
{
    for(int i:geometry().metrics().axisVelocityDofs())require(u.normalVelocity()[i]==0,"nonzero radial axis velocity");
    std::vector<long double> result;result.reserve(cells().size());
    for(std::size_t i=0;i<cells().size();++i){
        const auto velocity=static_cast<long double>(u.normalVelocity()[i]),defect=cells()[i].integratedDefect;
        const auto raw=velocity*defect;
        if(material().model()==HydroPhaseGeometry::Model::GeometricVof){
            // An enclosing typed closed consumer must retain u*R/dt and
            // .5*u*u*R whole, rather than consume this rounded impulse.
            Detail::hydroFiniteProductRange(velocity,defect,raw);
            if(!std::isfinite(raw)||(velocity!=0&&defect!=0&&(raw==0||std::fpclassify(raw)==FP_SUBNORMAL))){
                Detail::HydroFiniteContraction exact;exact.addProduct({velocity,defect});result.push_back(exact.rounded());continue;}
        }
        const auto value=finite(raw,"nonfinite consistency impulse");
        require(velocity==0 || defect==0 || value!=0,"unrepresentable consistency impulse");result.push_back(value);}
    return result;
}
std::array<long double,2> HydroDualMaterial::materialKineticChange(const QuadtreeMesh& mesh,const HydroNativeValues& u,
    HydroMaterialSources source,LAMP::TimeLayer clock) const
{
    const auto ordinary=[&](){
        const auto impulse=consistencyImpulse(mesh,u,source,clock);std::array<long double,2> result{};
        for(int i=0;i<layout().velocityCount();++i){const int component=layout().normalAxis(i)==HydroDofLayout::NormalAxis::First?0:1;
            const auto half=.5L*u.normalVelocity()[i],raw=half*impulse[i];
            Detail::hydroFiniteProductRange(.5L,u.normalVelocity()[i],half);
            Detail::hydroFiniteProductRange(half,impulse[i],raw);
            const auto value=finite(raw,"nonfinite material kinetic term");
            require(u.normalVelocity()[i]==0 || impulse[i]==0 || value!=0,"unrepresentable material kinetic term");result[component]+=value;
            if(Detail::hydroFiniteRangeEnabled()&&!std::isfinite(result[component]))throw Detail::HydroFiniteRange{};}
        for(auto value:result)finite(value,"nonfinite material kinetic total");return result;
    };
    if(material().model()!=HydroPhaseGeometry::Model::GeometricVof)return ordinary();
    try{Detail::HydroFiniteRangeScope range(true);return ordinary();}
    catch(const Detail::HydroFiniteRange&){
        // Validation/axis checks already ran in the ordinary impulse path.
        std::array<Detail::HydroFiniteContraction,2> exact;
        for(int i=0;i<layout().velocityCount();++i){const int a=layout().normalAxis(i)==HydroDofLayout::NormalAxis::First?0:1;
            exact[a].addProduct({.5L,u.normalVelocity()[i],u.normalVelocity()[i],cells()[i].integratedDefect});}
        return {exact[0].rounded(),exact[1].rounded()};
    }
}
HydroModalMaterialEnergyReport HydroDualMaterial::modalEnergyChange(const QuadtreeMesh& mesh,
    const HydroKineticHodge& begin,const HydroKineticHodge& midpoint,const HydroKineticHodge& end,
    const HydroSpatialAdvector& before,const HydroSpatialAdvector& after) const
{
    validateFor(mesh,sources(),timeLayer());
    require(&begin.materialMass().masses()==&beginMass().masses() && &midpoint.materialMass().masses()==&midpointMass().masses()
        && &end.materialMass().masses()==&endMass().masses(),"different exact modal material mass owners");
    require(material().model()==HydroPhaseGeometry::Model::FirstJet,"modal material energy requires FirstJet spatial modes");
    const auto& actual=material().midpointRate().spatialAdvector();
    require(actual.reconstruction().sameIdentity(midpoint.reconstruction()),"different midpoint modal material source");
    HydroModalMaterialEnergyReport report;
    report.kinetic=midpoint.energyChange(mesh,begin,end,before,after);
    report.legacyDiagonalMaterialChange=materialKineticChange(mesh,actual.values(),sources(),timeLayer());
    return report;
}
}
