#include "Modules/Hydrodynamics/CchPhaseViscosity.h"
#include "Modules/Hydrodynamics/CchMomentumGeometry.h"
#include "Modules/Hydrodynamics/CchSpatialAdvector.h"
#include "Modules/Hydrodynamics/Detail/CchRegularAdvectorGeometry.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <optional>
#include <stdexcept>
#include <utility>

namespace Quantum::Hydrodynamics
{
CchEnergyViscosity::Matrix CchPhaseViscosity::modalBulkMatrix(const QuadtreeMesh& mesh,const CchKineticHodge& hodge,
    double liquid,double gas)
{
    if(!std::isfinite(liquid)||liquid<0||!std::isfinite(gas)||gas<0)
        throw std::invalid_argument("[CchPhaseViscosity] negative/nonfinite phase viscosity");
    hodge.validateFor(mesh,hodge.geometryIdentity(),hodge.timeLayer());
    const auto& source=hodge.reconstruction();const auto& regular=source.regularGeometry();
    const auto& layout=hodge.layout();const int nu=layout.velocityCount(),nm=static_cast<int>(source.tangentialModes().size());
    CchEnergyViscosity::Matrix matrix(nu+nm,nu+nm);
    if(liquid==0&&gas==0){matrix.makeCompressed();return matrix;}
    const bool rz=layout.domain().axisymmetric();const double sqrt2=std::sqrt(2.);
    std::vector<std::optional<Detail::CchRegularCellStencil>> stencils(layout.pressureCount());
    std::map<std::pair<int,int>,long double> entries;
    for(const auto& q:hodge.strainQuadrature()){
        const double mu=q.liquid?liquid:gas;if(mu==0)continue;
        if(!stencils[q.pressure])stencils[q.pressure]=regular.cellStencil(mesh,q.pressure);
        struct Row{int column;std::array<double,4> strain;};std::vector<Row> rows;
        const auto strain=[&](double ux,double xx,double xy,double yx,double yy){
            const double hoop=rz?(q.point.x==0?xx:ux/q.point.x):0.;
            return std::array<double,4>{xx,yy,(xy+yx)/sqrt2,hoop};};
        for(const auto& weight:regular.cellJetWeights(*stencils[q.pressure],q.point))
            rows.push_back({weight.velocityDof,strain(weight.jet[0],weight.jet[2],weight.jet[3],weight.jet[4],weight.jet[5])});
        const int mode=source.tangentialModeForPressure(q.pressure);
        if(mode>=0){const auto jet=source.tangentialModeJet(mode,q.point,q.liquid);
            rows.push_back({nu+mode,strain(jet.value.x,jet.derivative.xx,jet.derivative.xy,jet.derivative.yx,jet.derivative.yy)});}
        for(std::size_t a=0;a<rows.size();++a)for(std::size_t b=a;b<rows.size();++b){long double product=0;
            for(int k=0;k<(rz?4:3);++k)product+=static_cast<long double>(rows[a].strain[k])*rows[b].strain[k];
            const auto key=std::minmax(rows[a].column,rows[b].column);
            entries[{key.first,key.second}]+=2*mu*q.physicalWeight*product;}
    }
    std::vector<Eigen::Triplet<double,Eigen::Index>> triplets;
    for(const auto& entry:entries){const auto [i,j]=entry.first;const double value=static_cast<double>(entry.second);
        if(!std::isfinite(value))throw std::runtime_error("[CchPhaseViscosity] nonfinite modal bulk entry");
        if(value==0)continue;triplets.emplace_back(i,j,value);if(i!=j)triplets.emplace_back(j,i,value);}
    matrix.setFromTriplets(triplets.begin(),triplets.end());matrix.makeCompressed();return matrix;
}
struct CchPhaseViscosity::Data
{
    Data(CchPhaseGeometry g,Interface::HermiteGalsPhaseQuadrature q,CchEnergyViscosity v,double l,double a)
        :geometry(std::move(g)),quadrature(std::move(q)),viscosity(std::move(v)),liquid(l),gas(a){}
    CchPhaseGeometry geometry;
    Interface::HermiteGalsPhaseQuadrature quadrature;
    CchEnergyViscosity viscosity;
    double liquid,gas;
};
CchPhaseViscosity::CchPhaseViscosity(std::shared_ptr<const Data> d):data_(std::move(d)){}
CchPhaseViscosity CchPhaseViscosity::build(const QuadtreeMesh& mesh,CchPhaseGeometry geometry,
    Interface::HermiteGalsPhaseQuadrature quadrature,double liquid,double gas)
{
    geometry.validateFor(mesh,quadrature.geometryIdentity(),quadrature.timeLayer());quadrature.validateMesh(mesh);
    if(!std::isfinite(liquid) || liquid<0 || !std::isfinite(gas) || gas<0)
        throw std::invalid_argument("[CchPhaseViscosity] negative/nonfinite phase viscosity");
    std::vector<CchStrainPoint> points;std::vector<double> mu;
    std::vector<bool> covered(geometry.layout().pressureCount(),false);
    for(const auto& q:quadrature.points()) {
        const int p=geometry.layout().pressureDof(q.cell);if(p<0)continue; // Explicit selected CCH subset.
        points.push_back({p,q.point,q.physicalWeight});mu.push_back(q.liquid?liquid:gas);covered[p]=true;
    }
    for(bool cell:covered)if(!cell)throw std::invalid_argument("[CchPhaseViscosity] missing selected-cell quadrature");
    auto viscosity=CchEnergyViscosity::build(CchStrainOperator::build(mesh,geometry.metrics(),points),mu);
    return CchPhaseViscosity(std::make_shared<Data>(std::move(geometry),std::move(quadrature),std::move(viscosity),liquid,gas));
}
const CchPhaseGeometry& CchPhaseViscosity::geometry() const noexcept{return data_->geometry;}
const Interface::HermiteGalsPhaseQuadrature& CchPhaseViscosity::quadrature() const noexcept{return data_->quadrature;}
const CchEnergyViscosity& CchPhaseViscosity::viscosity() const noexcept{return data_->viscosity;}
double CchPhaseViscosity::liquidViscosity() const noexcept{return data_->liquid;}
double CchPhaseViscosity::gasViscosity() const noexcept{return data_->gas;}
void CchPhaseViscosity::validateFor(const QuadtreeMesh& mesh,Interface::HermiteGalsIdentity id,Interface::HermiteGalsTimeLayer t) const
{geometry().validateFor(mesh,id,t);quadrature().validateMesh(mesh);}
std::vector<double> CchPhaseViscosity::apply(const QuadtreeMesh& mesh,const CchNativeValues& u,
    Interface::HermiteGalsIdentity id,Interface::HermiteGalsTimeLayer t) const
{validateFor(mesh,id,t);return viscosity().apply(mesh,u);}
double CchPhaseViscosity::dissipation(const QuadtreeMesh& mesh,const CchNativeValues& u,
    Interface::HermiteGalsIdentity id,Interface::HermiteGalsTimeLayer t) const
{validateFor(mesh,id,t);return viscosity().dissipation(mesh,u);}
}
