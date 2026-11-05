#include "Modules/Hydrodynamics/HydroPhaseViscosity.h"
#include "Modules/Hydrodynamics/HydroMomentumGeometry.h"
#include "Modules/Hydrodynamics/HydroSpatialAdvector.h"
#include "Modules/Hydrodynamics/Detail/HydroRegularAdvectorGeometry.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <optional>
#include <stdexcept>
#include <utility>

namespace LAMP::Hydrodynamics
{
HydroEnergyViscosity::Matrix HydroPhaseViscosity::modalBulkMatrix(const QuadtreeMesh& mesh,const HydroKineticHodge& hodge,
    double liquid,double gas)
{
    if(!std::isfinite(liquid)||liquid<0||!std::isfinite(gas)||gas<0)
        throw std::invalid_argument("[HydroPhaseViscosity] negative/nonfinite phase viscosity");
    hodge.validateFor(mesh,hodge.geometryIdentity(),hodge.timeLayer());
    const auto& source=hodge.reconstruction();const auto& regular=source.regularGeometry();
    const auto& layout=hodge.layout();const int nu=layout.velocityCount(),nm=static_cast<int>(source.tangentialModes().size());
    HydroEnergyViscosity::Matrix matrix(nu+nm,nu+nm);
    if(liquid==0&&gas==0){matrix.makeCompressed();return matrix;}
    const bool rz=layout.domain().axisymmetric();const double sqrt2=std::sqrt(2.);
    std::vector<std::optional<Detail::HydroRegularCellStencil>> stencils(layout.pressureCount());
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
        if(!std::isfinite(value))throw std::runtime_error("[HydroPhaseViscosity] nonfinite modal bulk entry");
        if(value==0)continue;triplets.emplace_back(i,j,value);if(i!=j)triplets.emplace_back(j,i,value);}
    matrix.setFromTriplets(triplets.begin(),triplets.end());matrix.makeCompressed();return matrix;
}
struct HydroPhaseViscosity::Data
{
    Data(HydroPhaseGeometry g,std::optional<Interface::HermiteGalsPhaseQuadrature> q,HydroEnergyViscosity v,double l,double a,
        MaterialModel m=MaterialModel::HermitePhaseQuadrature)
        :geometry(std::move(g)),quadrature(std::move(q)),viscosity(std::move(v)),liquid(l),gas(a),materialModel(m){}
    HydroPhaseGeometry geometry;
    std::optional<Interface::HermiteGalsPhaseQuadrature> quadrature;
    HydroEnergyViscosity viscosity;
    double liquid,gas;
    MaterialModel materialModel;
};
HydroPhaseViscosity::HydroPhaseViscosity(std::shared_ptr<const Data> d):data_(std::move(d)){}
HydroPhaseViscosity HydroPhaseViscosity::zero(const QuadtreeMesh& mesh,HydroPhaseGeometry geometry)
{
    if(geometry.model()==HydroPhaseGeometry::Model::GeometricVof)
        geometry.validateVofFor(mesh,geometry.vofIdentity(),geometry.timeLayer());
    else geometry.validateFor(mesh,geometry.geometryIdentity(),geometry.timeLayer());
    auto viscosity=HydroEnergyViscosity::zero(mesh,geometry.metrics());
    return HydroPhaseViscosity(std::make_shared<Data>(std::move(geometry),std::nullopt,std::move(viscosity),0,0,MaterialModel::ExactZero));
}
HydroPhaseViscosity HydroPhaseViscosity::build(const QuadtreeMesh& mesh,HydroPhaseGeometry geometry,
    Interface::HermiteGalsPhaseQuadrature quadrature,double liquid,double gas)
{
    geometry.validateFor(mesh,quadrature.geometryIdentity(),quadrature.timeLayer());quadrature.validateMesh(mesh);
    if(!std::isfinite(liquid) || liquid<0 || !std::isfinite(gas) || gas<0)
        throw std::invalid_argument("[HydroPhaseViscosity] negative/nonfinite phase viscosity");
    std::vector<HydroStrainPoint> points;std::vector<double> mu;
    std::vector<bool> covered(geometry.layout().pressureCount(),false);
    for(const auto& q:quadrature.points()) {
        const int p=geometry.layout().pressureDof(q.cell);if(p<0)continue; // Explicit selected HYDRO subset.
        points.push_back({p,q.point,q.physicalWeight});mu.push_back(q.liquid?liquid:gas);covered[p]=true;
    }
    for(bool cell:covered)if(!cell)throw std::invalid_argument("[HydroPhaseViscosity] missing selected-cell quadrature");
    auto viscosity=HydroEnergyViscosity::build(HydroStrainOperator::build(mesh,geometry.metrics(),points),mu);
    return HydroPhaseViscosity(std::make_shared<Data>(std::move(geometry),std::move(quadrature),std::move(viscosity),liquid,gas));
}
HydroPhaseViscosity HydroPhaseViscosity::buildVofCellAverage(const QuadtreeMesh& mesh,
    HydroPhaseGeometry geometry,double liquid,double gas)
{return buildVofCellAverageImpl(mesh,std::move(geometry),liquid,gas,nullptr);}
HydroPhaseViscosity HydroPhaseViscosity::buildVofCellAverage(const QuadtreeMesh& mesh,
    HydroPhaseGeometry geometry,double liquid,double gas,const HydroPhaseViscosity& previous)
{return buildVofCellAverageImpl(mesh,std::move(geometry),liquid,gas,&previous);}
HydroPhaseViscosity HydroPhaseViscosity::buildVofCellAverage(const QuadtreeMesh& mesh,
    HydroPhaseGeometry geometry,double liquid,double gas,const HydroStrainOperator& canonicalStressGeometry)
{return buildVofCellAverageImpl(mesh,std::move(geometry),liquid,gas,nullptr,&canonicalStressGeometry);}
HydroPhaseViscosity HydroPhaseViscosity::buildVofCellAverage(const QuadtreeMesh& mesh,
    HydroPhaseGeometry geometry,double liquid,double gas,const std::vector<int>& walls)
{return buildVofCellAverageImpl(mesh,std::move(geometry),liquid,gas,nullptr,nullptr,&walls);}
HydroPhaseViscosity HydroPhaseViscosity::buildVofCellAverage(const QuadtreeMesh& mesh,
    HydroPhaseGeometry geometry,double liquid,double gas,const std::vector<int>& walls,const HydroPhaseViscosity& previous)
{return buildVofCellAverageImpl(mesh,std::move(geometry),liquid,gas,&previous,nullptr,&walls);}
HydroPhaseViscosity HydroPhaseViscosity::buildVofCellAverage(const QuadtreeMesh& mesh,
    HydroPhaseGeometry geometry,double liquid,double gas,const std::vector<int>& walls,
    const HydroMomentumGeometry& owner)
{return buildVofCellAverageImpl(mesh,std::move(geometry),liquid,gas,nullptr,nullptr,&walls,&owner);}
HydroPhaseViscosity HydroPhaseViscosity::buildVofCellAverageImpl(const QuadtreeMesh& mesh,
    HydroPhaseGeometry geometry,double liquid,double gas,const HydroPhaseViscosity* previous,
    const HydroStrainOperator* canonicalStressGeometry,const std::vector<int>* stationaryWalls,
    const HydroMomentumGeometry* rawGeometryOwner)
{
    if(geometry.model()!=HydroPhaseGeometry::Model::GeometricVof)
        throw std::invalid_argument("[HydroPhaseViscosity] cell-average VOF material requires its genuine phase owner");
    geometry.validateVofFor(mesh,geometry.vofIdentity(),geometry.timeLayer());
    const auto& source=geometry.vofPublication();
    if(source.geometryKind()!=Interface::GeometricVofGeometryKind::LocalPlic)
        throw std::invalid_argument("[HydroPhaseViscosity] cell-average VOF material requires actual LocalPlic inventories");
    if(!std::isfinite(liquid)||liquid<0||!std::isfinite(gas)||gas<0)
        throw std::invalid_argument("[HydroPhaseViscosity] negative/nonfinite phase viscosity");
    const std::vector<int> emptyWalls;
    const auto& walls=stationaryWalls?*stationaryWalls:emptyWalls;
    if(!std::is_sorted(walls.begin(),walls.end())||std::adjacent_find(walls.begin(),walls.end())!=walls.end())
        throw std::invalid_argument("[HydroPhaseViscosity] stationary no-slip wall rows must be sorted and unique");
    for(const int u:walls){
        if(u<0||u>=geometry.layout().velocityCount()||geometry.metrics().faceMeasures()[u]<=0||
            geometry.layout().domain().face(mesh,geometry.layout().nativeFace(u)).outwardSign==0)
            throw std::invalid_argument("[HydroPhaseViscosity] stationary no-slip requires positive-area physical wall rows");
    }
    if(previous){
        if(previous->materialModel()!=MaterialModel::VofCellAverage)
            throw std::invalid_argument("[HydroPhaseViscosity] prepared VOF viscosity requires a previous VofCellAverage source");
        const auto& prior=previous->geometry();
        prior.validateVofFor(mesh,prior.vofIdentity(),prior.timeLayer());
        if(!prior.layout().sameIdentity(geometry.layout()))
            throw std::invalid_argument("[HydroPhaseViscosity] prepared VOF viscosity has a different actual native layout");
        if(&prior.metrics().divergence().matrix()!=&geometry.metrics().divergence().matrix())
            throw std::invalid_argument("[HydroPhaseViscosity] prepared VOF viscosity must retain the actual phase divergence owner");
        const auto& strain=previous->viscosity().strain();
        if(strain.stationaryNoSlipVelocityDofs()!=walls)
            throw std::invalid_argument("[HydroPhaseViscosity] previous raw strain has a different stationary no-slip boundary");
        if(strain.isStressCompatible()||!strain.layout().sameIdentity(geometry.layout())||
            &strain.metrics().divergence().matrix()!=&geometry.metrics().divergence().matrix())
            throw std::invalid_argument("[HydroPhaseViscosity] previous VofCellAverage must retain its canonical raw strain owner");
    }
    if(liquid==0&&gas==0)return zero(mesh,std::move(geometry));
    std::optional<HydroStrainOperator> preparedRaw;
    if(canonicalStressGeometry){
        if(!canonicalStressGeometry->isStressCompatible()||!canonicalStressGeometry->hasCanonicalRawStrain())
            throw std::invalid_argument("[HydroPhaseViscosity] prepared stress source lacks the builtin canonical raw strain");
        canonicalStressGeometry->layout().validateMesh(mesh);
        if(!canonicalStressGeometry->layout().sameIdentity(geometry.layout()))
            throw std::invalid_argument("[HydroPhaseViscosity] prepared stress raw strain has a different actual native layout");
        if(&canonicalStressGeometry->metrics().divergence().matrix()!=&geometry.metrics().divergence().matrix())
            throw std::invalid_argument("[HydroPhaseViscosity] prepared stress raw strain must retain the actual phase divergence owner");
        preparedRaw.emplace(canonicalStressGeometry->canonicalRawStrain());
        if(preparedRaw->isStressCompatible()||!preparedRaw->layout().sameIdentity(geometry.layout())||
            &preparedRaw->metrics().divergence().matrix()!=&geometry.metrics().divergence().matrix())
            throw std::invalid_argument("[HydroPhaseViscosity] builtin stress source lost its actual canonical raw strain owner");
    }
    const auto points=previous?previous->viscosity().strain().points():preparedRaw?preparedRaw->points():
        HydroStrainOperator::volumeQuadrature(mesh,geometry.metrics());
    std::vector<double> mu;mu.reserve(points.size());
    for(const auto& q:points)
    {
        const auto native=geometry.layout().nativeCell(q.cellDof);
        const auto& stock=source.plicCells().at(native).inventory;
        const auto capacity=source.cells().at(native).capacity;
        // Both phase stocks are authoritative. Do not recover either from
        // alpha, the other stock, or the actual PLIC plane's measured phase.
        const auto value=(stock[0]/capacity)*liquid+(stock[1]/capacity)*gas;
        const auto stored=static_cast<double>(value);
        if(!std::isfinite(value)||value<0||!std::isfinite(stored)||
            (value>0&&stored==0)||(value==0&&((stock[0]>0&&liquid>0)||(stock[1]>0&&gas>0))))
            throw std::invalid_argument("[HydroPhaseViscosity] unrepresentable nonnegative inventory P0 viscosity");
        mu.push_back(stored);
    }
    auto viscosity=[&]{
        if(previous)return HydroEnergyViscosity::build(previous->viscosity().strain(),mu);
        if(preparedRaw)return HydroEnergyViscosity::build(*preparedRaw,mu);
        if(rawGeometryOwner&&!walls.empty())return HydroEnergyViscosity::buildStationaryNoSlipOwned(
            mesh,*rawGeometryOwner,geometry.metrics(),points,walls,mu);
        return HydroEnergyViscosity::build(walls.empty()?HydroStrainOperator::build(mesh,geometry.metrics(),points)
            :HydroStrainOperator::buildStationaryNoSlip(mesh,geometry.metrics(),points,walls),mu);
    }();
    return HydroPhaseViscosity(std::make_shared<Data>(std::move(geometry),std::nullopt,std::move(viscosity),
        liquid,gas,MaterialModel::VofCellAverage));
}
HydroPhaseViscosity::MaterialModel HydroPhaseViscosity::materialModel() const noexcept{return data_->materialModel;}
const HydroPhaseGeometry& HydroPhaseViscosity::geometry() const noexcept{return data_->geometry;}
bool HydroPhaseViscosity::hasQuadrature() const noexcept{return data_->quadrature.has_value();}
const Interface::HermiteGalsPhaseQuadrature& HydroPhaseViscosity::quadrature() const
{
    if(!hasQuadrature())throw std::invalid_argument("[HydroPhaseViscosity] this material model has no Hermite phase quadrature");
    return *data_->quadrature;
}
const HydroEnergyViscosity& HydroPhaseViscosity::viscosity() const noexcept{return data_->viscosity;}
double HydroPhaseViscosity::liquidViscosity() const noexcept{return data_->liquid;}
double HydroPhaseViscosity::gasViscosity() const noexcept{return data_->gas;}
void HydroPhaseViscosity::validateFor(const QuadtreeMesh& mesh,Interface::HermiteGalsIdentity id,LAMP::TimeLayer t) const
{geometry().validateFor(mesh,id,t);if(hasQuadrature())quadrature().validateMesh(mesh);}
std::vector<double> HydroPhaseViscosity::apply(const QuadtreeMesh& mesh,const HydroNativeValues& u,
    Interface::HermiteGalsIdentity id,LAMP::TimeLayer t) const
{validateFor(mesh,id,t);return viscosity().apply(mesh,u);}
double HydroPhaseViscosity::dissipation(const QuadtreeMesh& mesh,const HydroNativeValues& u,
    Interface::HermiteGalsIdentity id,LAMP::TimeLayer t) const
{validateFor(mesh,id,t);return viscosity().dissipation(mesh,u);}
}
