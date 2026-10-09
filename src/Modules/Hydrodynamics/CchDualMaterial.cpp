#include "Modules/Hydrodynamics/CchDualMaterial.h"
#include "Modules/Hydrodynamics/CchSpatialAdvector.h"

#include "Core/Mesh/QuadtreeMesh.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

namespace Quantum::Hydrodynamics
{
namespace
{
constexpr long double tolerance=256.L*std::numeric_limits<double>::epsilon();
void require(bool ok,const char* why)
{ if(!ok)throw std::invalid_argument(std::string("[CchDualMaterial] ")+why); }
long double finite(long double v,const char* why)
{ require(std::isfinite(v),why);return v; }
long double defect(long double difference,long double scale,const char* why)
{
    finite(difference,why);finite(scale,why);require(scale>=0,why);
    require(scale>0 || difference==0,why);
    return scale>0?std::abs(difference)/scale:0;
}
long double change(CchPhaseMeasure a,CchPhaseMeasure b,double l,double g)
{
    const auto dl=b.liquid-a.liquid,dg=b.gas-a.gas;
    const auto ml=finite(l*dl,"nonfinite dual liquid mass change"),mg=finite(g*dg,"nonfinite dual gas mass change");
    require((dl==0 || ml!=0) && (dg==0 || mg!=0),"unrepresentable dual mass change");
    return finite(ml+mg,"nonfinite dual mass change");
}
}
struct CchDualMaterial::Data
{
    Data(CchMaterialTransport t,CchMomentumMass a,CchMomentumMass m,CchMomentumMass b)
        :transport(std::move(t)),begin(std::move(a)),midpoint(std::move(m)),end(std::move(b)){}
    CchMaterialTransport transport;
    CchMomentumMass begin,midpoint,end;
    std::vector<long double> integrals;
    std::vector<CchDualMaterialCell> cells;
    CchDualMaterialReport report;
};
CchDualMaterial::CchDualMaterial(std::shared_ptr<const Data> d):data_(std::move(d)){}
CchDualMaterial CchDualMaterial::prepare(const QuadtreeMesh& mesh,CchMaterialTransport material,
    CchMomentumMass begin,CchMomentumMass midpoint,CchMomentumMass end)
{
    material.validateFor(mesh,material.sources(),material.timeLayer());
    const auto source=material.sources();
    begin.validateFor(mesh,source.begin,material.beginGeometry().timeLayer());
    midpoint.validateFor(mesh,source.midpoint,material.timeLayer());
    end.validateFor(mesh,source.end,material.endGeometry().timeLayer());
    require(begin.geometry().sameIdentity(midpoint.geometry()) && end.geometry().sameIdentity(midpoint.geometry()),"different dual geometry identity");
    for(const auto* mass:std::array<const CchMomentumMass*,3>{{&begin,&midpoint,&end}}){
        require(mass->layout().sameIdentity(material.layout()),"different native layout");
        require(&mass->geometry().metrics().divergence().matrix()==&material.midpointRate().geometry().metrics().divergence().matrix(),"different retained D_I");
        require(mass->liquidDensity()==material.liquidDensity() && mass->gasDensity()==material.gasDensity(),"different material densities");
    }
    material.validateFor(mesh,midpoint.geometry(),source,material.timeLayer());
    auto d=std::make_shared<Data>(std::move(material),std::move(begin),std::move(midpoint),std::move(end));
    const auto& geometry=d->midpoint.geometry();const auto& layout=geometry.layout();const double dt=d->transport.deltaTime();
    const auto& rates=d->transport.dualFaceRates();
    require(rates.size()==geometry.segments().size(),"different original dual flux view");
    d->integrals=d->transport.dualFaceIntegrals();
    require(d->integrals.size()==rates.size(),"different original dual integral view");
    for(auto integral:d->integrals)finite(integral,"nonfinite dual material integral");
    std::vector<long double> primalRate(layout.pressureCount()),primalGross(layout.pressureCount());
    const auto& di=d->transport.midpointRate().geometry().metrics().divergence().matrix();
    for(Eigen::Index p=0;p<di.outerSize();++p)for(CchIntegralDivergence::Matrix::InnerIterator it(di,p);it;++it){
        const auto rate=d->transport.faceRates()[static_cast<std::size_t>(it.col())];
        primalRate[p]+=(it.value()>0?1:-1)*rate;primalGross[p]+=std::abs(rate);
    }
    d->cells.resize(layout.velocityCount());std::array<long double,2> balanceScale{};
    for(int u=0;u<layout.velocityCount();++u){
        const int component=layout.normalAxis(u)==CchDofLayout::NormalAxis::First?0:1;
        auto& cell=d->cells[u];auto& report=d->report.components[component];
        cell.massChange=change(d->begin.phaseVolumes()[u],d->end.phaseVolumes()[u],d->transport.liquidDensity(),d->transport.gasDensity());
        long double outwardRate=0,grossRate=0,grossIntegral=0,mappedRate=0,mappedGross=0;
        for(int side=0;side<4;++side)for(int id:geometry.volumes()[u].sides[side]){
            const int sign=side%2==0?-1:1;
            outwardRate+=sign*rates[id];cell.outwardIntegral+=sign*d->integrals[id];
            grossRate+=std::abs(rates[id]);grossIntegral+=std::abs(d->integrals[id]);
        }
        for(const auto& overlap:geometry.volumes()[u].overlaps){
            const auto& primal=d->transport.cells()[overlap.pressure];const auto w=overlap.cellFraction;
            mappedRate+=w*primalRate[overlap.pressure];mappedGross+=w*primalGross[overlap.pressure];
            cell.mappedPrimalMassChange+=w*primal.massChange;cell.mappedPrimalIntegratedDefect+=w*primal.integratedDefect;
        }
        cell.outwardRate=finite(outwardRate,"nonfinite actual dual outward rate");
        cell.mappedPrimalOutwardRate=finite(mappedRate,"nonfinite mapped primal outward rate");
        cell.spatialExchangeRate=finite(outwardRate-mappedRate,"nonfinite spatial exchange rate");
        const auto exchange=defect(cell.spatialExchangeRate,grossRate+mappedGross,"invalid spatial exchange diagnostic");
        d->report.maxRelativeSpatialExchange=std::max(d->report.maxRelativeSpatialExchange,exchange);
        d->report.maxAbsSpatialExchange=std::max(d->report.maxAbsSpatialExchange,std::abs(cell.spatialExchangeRate));
        cell.geometryTransferChange=finite(cell.massChange-cell.mappedPrimalMassChange,"nonfinite geometry transfer change");
        cell.integratedDefect=finite(cell.massChange+cell.outwardIntegral,"nonfinite dual integrated defect");
        cell.defectRate=finite(cell.integratedDefect/dt,"nonfinite dual defect rate");
        const auto gross=finite(std::abs(cell.massChange)+grossIntegral,"nonfinite dual defect scale");
        cell.relativeDefect=defect(cell.integratedDefect,gross,"invalid relative dual defect");
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
        const int component=layout.normalAxis(u)==CchDofLayout::NormalAxis::First?0:1;
        d->report.components[component].boundaryIntegral+=(s.negativeVelocity>=0?1:-1)*d->integrals[id];
    }
    for(int component=0;component<2;++component){auto& report=d->report.components[component];
        for(auto v:report.geometryMass)finite(v,"nonfinite total dual mass");
        for(auto v:{report.massChange,report.boundaryIntegral,report.integratedDefect,report.l1IntegratedDefect})finite(v,"nonfinite total dual balance");
        report.balanceIdentity=defect(report.integratedDefect-report.massChange-report.boundaryIntegral,balanceScale[component],"dual global balance defect");
        require(report.balanceIdentity<=tolerance,"dual global balance identity failed");
    }
    return CchDualMaterial(std::move(d));
}
const CchMaterialTransport& CchDualMaterial::material() const noexcept{return data_->transport;}
const CchMomentumMass& CchDualMaterial::beginMass() const noexcept{return data_->begin;}
const CchMomentumMass& CchDualMaterial::midpointMass() const noexcept{return data_->midpoint;}
const CchMomentumMass& CchDualMaterial::endMass() const noexcept{return data_->end;}
const CchMomentumGeometry& CchDualMaterial::geometry() const noexcept{return midpointMass().geometry();}
const CchDofLayout& CchDualMaterial::layout() const noexcept{return material().layout();}
CchMaterialSources CchDualMaterial::sources() const{return material().sources();}
Interface::HermiteGalsTimeLayer CchDualMaterial::timeLayer() const{return material().timeLayer();}
double CchDualMaterial::deltaTime() const noexcept{return material().deltaTime();}
const std::vector<long double>& CchDualMaterial::faceRates() const noexcept{return material().dualFaceRates();}
const std::vector<long double>& CchDualMaterial::faceIntegrals() const noexcept{return data_->integrals;}
const std::vector<CchDualMaterialCell>& CchDualMaterial::cells() const noexcept{return data_->cells;}
const CchDualMaterialReport& CchDualMaterial::report() const noexcept{return data_->report;}
void CchDualMaterial::validateFor(const QuadtreeMesh& mesh,CchMaterialSources source,Interface::HermiteGalsTimeLayer clock) const
{material().validateFor(mesh,source,clock);require(material().momentumGeometry().sameIdentity(geometry()),"different momentum geometry instance");geometry().validateMesh(mesh);}
std::vector<long double> CchDualMaterial::consistencyImpulse(const QuadtreeMesh& mesh,const CchNativeValues& u,
    CchMaterialSources source,Interface::HermiteGalsTimeLayer clock) const
{
    validateFor(mesh,source,clock);u.validateFor(layout(),mesh);
    for(int i:geometry().metrics().axisVelocityDofs())require(u.normalVelocity()[i]==0,"nonzero radial axis velocity");
    std::vector<long double> result;result.reserve(cells().size());
    for(std::size_t i=0;i<cells().size();++i){const auto value=finite(u.normalVelocity()[i]*cells()[i].integratedDefect,"nonfinite consistency impulse");
        require(u.normalVelocity()[i]==0 || cells()[i].integratedDefect==0 || value!=0,"unrepresentable consistency impulse");result.push_back(value);}
    return result;
}
std::array<long double,2> CchDualMaterial::materialKineticChange(const QuadtreeMesh& mesh,const CchNativeValues& u,
    CchMaterialSources source,Interface::HermiteGalsTimeLayer clock) const
{
    const auto impulse=consistencyImpulse(mesh,u,source,clock);std::array<long double,2> result{};
    for(int i=0;i<layout().velocityCount();++i){const int component=layout().normalAxis(i)==CchDofLayout::NormalAxis::First?0:1;
        const auto value=finite((.5L*u.normalVelocity()[i])*impulse[i],"nonfinite material kinetic term");
        require(u.normalVelocity()[i]==0 || impulse[i]==0 || value!=0,"unrepresentable material kinetic term");result[component]+=value;}
    for(auto value:result)finite(value,"nonfinite material kinetic total");return result;
}
CchModalMaterialEnergyReport CchDualMaterial::modalEnergyChange(const QuadtreeMesh& mesh,
    const CchKineticHodge& begin,const CchKineticHodge& midpoint,const CchKineticHodge& end,
    const CchSpatialAdvector& before,const CchSpatialAdvector& after) const
{
    validateFor(mesh,sources(),timeLayer());
    require(&begin.materialMass().masses()==&beginMass().masses() && &midpoint.materialMass().masses()==&midpointMass().masses()
        && &end.materialMass().masses()==&endMass().masses(),"different exact modal material mass owners");
    const auto& actual=material().midpointRate().spatialAdvector();
    require(actual.reconstruction().sameIdentity(midpoint.reconstruction()),"different midpoint modal material source");
    CchModalMaterialEnergyReport report;
    report.kinetic=midpoint.energyChange(mesh,begin,end,before,after);
    report.legacyDiagonalMaterialChange=materialKineticChange(mesh,actual.values(),sources(),timeLayer());
    return report;
}
}
