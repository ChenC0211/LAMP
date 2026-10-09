#include "Modules/Hydrodynamics/CchInterfaceNormalLoad.h"
#include "Modules/Hydrodynamics/CchSpatialAdvector.h"
#include "Modules/Hydrodynamics/Detail/CchRegularAdvectorGeometry.h"
#include <algorithm>
#include <cmath>
#include <map>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>

namespace Quantum::Hydrodynamics
{
namespace
{
void require(bool value,const char* reason)
{
    if(!value)throw std::invalid_argument(std::string("[CchInterfaceNormalLoad] ")+reason);
}
}
namespace Detail
{
CchInterfaceStressAssembly assembleCchInterfaceStress(const QuadtreeMesh& mesh,
    const CchPressureLoad& pressure,const CchPhaseGeometry& geometry,
    const CchAdvectorReconstruction& source,const CchAdvectorBoundary* boundary,
    const std::vector<double>& stress)
{
    require(pressure.layout().sameIdentity(geometry.layout()),"different pressure/geometry layout");
    require(&pressure.divergence().matrix()==&geometry.metrics().divergence().matrix(),"different authoritative divergence owner");
    require(source.momentumGeometry().layout().sameIdentity(pressure.layout()),"different source layout");
    require(&source.momentumGeometry().metrics().divergence().matrix()==&pressure.divergence().matrix(),"different source divergence owner");
    const auto& regular=source.regularGeometry();regular.validateFor(mesh,source.momentumGeometry());
    if(boundary)boundary->validateFor(mesh,source.momentumGeometry());geometry.validateFor(mesh,geometry.geometryIdentity(),geometry.timeLayer());
    for(double j:stress)require(std::isfinite(j),"nonfinite surface normal stress");
    const auto& layout=pressure.layout();const int count=layout.velocityCount();const auto& publication=geometry.publication();
    const auto incidence=pressure.apply(mesh,CchNativeValues(layout,std::vector<double>(count),geometry.liquidFractions()));
    struct Point { long double weight=0,delta=0,incidenceScale=0;
        std::vector<std::pair<int,double>> normal;
        std::vector<std::pair<int,long double>> functional; };
    std::vector<Point> points;std::map<int,CchRegularCellStencil> stencils;
    const long double anchor=stress.empty()?0:stress.front();std::vector<long double> basisScale(count),weightScale(count),surfaceMeasure(count);
    std::size_t index=0;
    for(std::size_t c=0;c<publication.topology().curves.size();++c)for(std::size_t g=0;g<publication.topology().curves[c].graphs().size();++g)for(const auto& q:publication.quadrature(c,g)){
        require(index<stress.size(),"surface stress sample inventory");const long double delta=static_cast<long double>(stress[index++])-anchor;
        const auto& frame=q.sample.frame;const int p=layout.pressureDof(frame.nativeCellIndex);if(p<0)continue;
        require(std::isfinite(q.physicalWeight)&&q.physicalWeight>=0,"invalid physical surface weight");if(q.physicalWeight==0)continue;
        auto prepared=stencils.find(p);if(prepared==stencils.end())prepared=stencils.emplace(p,regular.cellStencil(mesh,p,boundary)).first;
        Point point;point.weight=q.physicalWeight;point.delta=delta;
        for(const auto& row:regular.cellJetWeights(prepared->second,frame.point)){
            const double b=frame.normal.x*row.jet[0]+frame.normal.y*row.jet[1];require(std::isfinite(b),"nonfinite normal source weight");if(b==0)continue;
            point.normal.emplace_back(row.velocityDof,b);basisScale[row.velocityDof]=std::max(basisScale[row.velocityDof],std::abs(static_cast<long double>(b)));
            weightScale[row.velocityDof]=std::max(weightScale[row.velocityDof],point.weight);surfaceMeasure[row.velocityDof]+=point.weight;
            point.incidenceScale=std::max(point.incidenceScale,std::abs(static_cast<long double>(incidence[row.velocityDof])));}
        if(!point.normal.empty()){
            require(point.incidenceScale>0,"normal stress footprint has zero local pressure incidence");long double norm=0;
            for(const auto& row:point.normal){const long double v=incidence[row.first]/point.incidenceScale;norm+=v*v;}
            // Dimensionless functional; weight/scale is combined only when
            // forming a matrix entry. Never square a tiny physical incidence.
            for(const auto& row:point.normal)if(incidence[row.first]!=0)point.functional.emplace_back(row.first,(incidence[row.first]/point.incidenceScale)/norm);
        }
        points.push_back(std::move(point));
    }
    require(index==stress.size(),"surface stress sample inventory");
    std::vector<long double> positive(count),offset(count),averageDelta(count),correction(count);
    for(const auto& point:points)for(const auto& row:point.normal){const int u=row.first;const long double b=row.second/basisScale[u];
        const long double weight=(point.weight/weightScale[u])*b*b;positive[u]+=weight;offset[u]+=weight*point.delta;}
    for(int u=0;u<count;++u)if(surfaceMeasure[u]>0){require(positive[u]>0,"lost positive surface anchor");averageDelta[u]=offset[u]/positive[u];}
    CchInterfaceStressAssembly out;out.constantChannel=static_cast<double>(anchor);out.matrix.resize(count,count);
    out.load.resize(count);out.averages.resize(count);out.surfaceMeasures.resize(count);out.anchorMeasures.resize(count);out.constantReactions.assign(count,0);
    std::vector<Eigen::Triplet<double>> triplets;
    for(int u=0;u<count;++u){const long double mean=anchor+averageDelta[u];const double value=static_cast<double>(mean);require(std::isfinite(value),"nonfinite stress anchor");
        if(value!=0)triplets.emplace_back(u,u,value);out.averages[u]=surfaceMeasure[u]>0?value:0;
        out.surfaceMeasures[u]=static_cast<double>(surfaceMeasure[u]);out.anchorMeasures[u]=static_cast<double>(weightScale[u]*basisScale[u]*basisScale[u]*positive[u]);}
    for(const auto& point:points)for(const auto& row:point.normal){const int u=row.first;
        const long double delta=point.delta-averageDelta[u],force=-point.weight*row.second*delta;correction[u]+=force;
        if(delta==0)continue;const long double factor=-(point.weight/point.incidenceScale)*row.second*delta;
        for(const auto& column:point.functional){const double entry=static_cast<double>(factor*column.second);require(std::isfinite(entry),"nonfinite normal stress matrix entry");if(entry!=0)triplets.emplace_back(u,column.first,entry);}
    }
    out.matrix.setFromTriplets(triplets.begin(),triplets.end());out.matrix.makeCompressed();
    for(int k=0;k<out.matrix.nonZeros();++k)require(std::isfinite(out.matrix.valuePtr()[k]),"nonfinite assembled normal stress matrix");
    for(int u=0;u<count;++u){out.load[u]=static_cast<double>((anchor+averageDelta[u])*incidence[u]+correction[u]);
        require(std::isfinite(out.load[u])&&std::isfinite(out.surfaceMeasures[u])&&std::isfinite(out.anchorMeasures[u]),"nonfinite anchored weak normal load");}
    if(boundary)for(int u:boundary->velocityDofs())out.constantReactions[u]=static_cast<double>(anchor*incidence[u]);
    return out;
}
}
struct CchInterfaceNormalLoad::Data
{
    Data(CchPressureLoad p,CchPhaseGeometry g,double j)
        :pressure(std::move(p)),geometry(std::move(g)),jump(j){}
    CchPressureLoad pressure;
    CchPhaseGeometry geometry;
    double jump;
    bool capillary=false;
    double sigma=0;
    std::vector<double> curvature,measure;
    std::vector<double> load;
    CchPressureLoad::Matrix stress;
    std::vector<double> reactions;
    std::vector<double> tangential;
    std::optional<CchAdvectorReconstruction> source;
    std::optional<CchAdvectorBoundary> boundary;
};
CchInterfaceNormalLoad::CchInterfaceNormalLoad(std::shared_ptr<const Data> data):data_(std::move(data)){}
CchInterfaceNormalLoad CchInterfaceNormalLoad::constant(const QuadtreeMesh& mesh,
    CchPressureLoad pressure,CchPhaseGeometry geometry,double jump)
{
    require(std::isfinite(jump),"nonfinite constant normal stress");
    require(pressure.layout().sameIdentity(geometry.layout()),"different pressure/geometry layout");
    require(&pressure.divergence().matrix()==&geometry.metrics().divergence().matrix(),
        "different authoritative divergence owner");
    geometry.validateFor(mesh,geometry.geometryIdentity(),geometry.timeLayer());
    auto data=std::make_shared<Data>(std::move(pressure),std::move(geometry),jump);
    std::vector<double> potential=data->geometry.liquidFractions();
    for(double& value:potential){value*=jump;require(std::isfinite(value),"nonfinite normal stress potential");}
    data->load=data->pressure.apply(mesh,CchNativeValues(data->pressure.layout(),
        std::vector<double>(data->pressure.layout().velocityCount(),0),potential));
    const int n=data->pressure.layout().velocityCount();data->stress.resize(n,n);
    for(int u=0;u<n;++u)if(jump!=0)data->stress.insert(u,u)=jump;data->stress.makeCompressed();
    return CchInterfaceNormalLoad(std::move(data));
}
CchInterfaceNormalLoad CchInterfaceNormalLoad::capillary(const QuadtreeMesh& mesh,
    CchPressureLoad pressure,CchPhaseGeometry geometry,double sigma,
    const CchAdvectorReconstruction& source,const CchAdvectorBoundary* boundary)
{
    require(std::isfinite(sigma) && sigma>=0,"invalid surface tension");
    require(pressure.layout().sameIdentity(geometry.layout()),"different pressure/geometry layout");
    require(&pressure.divergence().matrix()==&geometry.metrics().divergence().matrix(),
        "different authoritative divergence owner");
    geometry.validateFor(mesh,geometry.geometryIdentity(),geometry.timeLayer());
    auto data=std::make_shared<Data>(std::move(pressure),std::move(geometry),0);
    data->capillary=true;data->sigma=sigma;
    data->source=source;if(boundary)data->boundary=*boundary;
    if(const auto* bound=source.phaseGeometry())bound->validateFor(mesh,data->geometry.geometryIdentity(),data->geometry.timeLayer());
    data->tangential.assign(source.tangentialModes().size(),0.);
    const auto& publication=data->geometry.publication();std::vector<double> stress;
    for(std::size_t c=0;c<publication.topology().curves.size();++c)
        for(std::size_t g=0;g<publication.topology().curves[c].graphs().size();++g)
            for(const auto& q:publication.quadrature(c,g)){
                stress.push_back(-q.sample.curvature);
            }
    auto assembled=Detail::assembleCchInterfaceStress(mesh,data->pressure,data->geometry,source,boundary,stress);
    data->stress=sigma*assembled.matrix;data->load=std::move(assembled.load);data->measure=std::move(assembled.surfaceMeasures);
    data->curvature=std::move(assembled.averages);data->reactions=std::move(assembled.constantReactions);
    for(std::size_t u=0;u<data->load.size();++u){data->load[u]*=sigma;data->curvature[u]=-data->curvature[u];data->reactions[u]*=sigma;
        require(std::isfinite(data->load[u])&&std::isfinite(data->reactions[u]),"nonfinite capillary load");
    }
    for(int k=0;k<data->stress.nonZeros();++k)require(std::isfinite(data->stress.valuePtr()[k]),"nonfinite capillary stress matrix");
    return CchInterfaceNormalLoad(std::move(data));
}
const CchPressureLoad& CchInterfaceNormalLoad::pressureLoad() const noexcept{return data_->pressure;}
const CchPhaseGeometry& CchInterfaceNormalLoad::geometry() const noexcept{return data_->geometry;}
const CchDofLayout& CchInterfaceNormalLoad::layout() const noexcept{return data_->pressure.layout();}
bool CchInterfaceNormalLoad::isCapillary() const noexcept{return data_->capillary;}
double CchInterfaceNormalLoad::constantJump() const{require(!data_->capillary,"capillary load has no prescribed constant jump");return data_->jump;}
double CchInterfaceNormalLoad::surfaceTension() const{require(data_->capillary,"constant load has no surface tension");return data_->sigma;}
const std::vector<double>& CchInterfaceNormalLoad::curvatureAverages() const{require(data_->capillary,"constant load has no curvature averages");return data_->curvature;}
const std::vector<double>& CchInterfaceNormalLoad::supportMeasures() const{require(data_->capillary,"constant load has no surface supports");return data_->measure;}
const std::vector<double>& CchInterfaceNormalLoad::balanceLoads() const noexcept{return data_->load;}
const std::vector<double>& CchInterfaceNormalLoad::tangentialBalanceLoads() const
{require(data_->capillary,"constant load has no tangential source binding");return data_->tangential;}
const CchPressureLoad::Matrix& CchInterfaceNormalLoad::stressMatrix() const noexcept{return data_->stress;}
const CchAdvectorReconstruction& CchInterfaceNormalLoad::reconstruction() const{require(data_->capillary,"constant load has no source reconstruction");return *data_->source;}
const CchAdvectorBoundary* CchInterfaceNormalLoad::boundary() const{require(data_->capillary,"constant load has no source boundary");return data_->boundary?&*data_->boundary:nullptr;}
const std::vector<double>& CchInterfaceNormalLoad::constantReactions() const{require(data_->capillary,"constant load has no bound wall reactions");return data_->reactions;}
void CchInterfaceNormalLoad::validateFor(const QuadtreeMesh& mesh,const CchPressureLoad& pressure,
    Interface::HermiteGalsIdentity identity,Interface::HermiteGalsTimeLayer clock) const
{
    require(&pressure.matrix()==&data_->pressure.matrix(),"different authoritative pressure owner");
    data_->geometry.validateFor(mesh,identity,clock);
}
void CchInterfaceNormalLoad::validateFor(const QuadtreeMesh& mesh,const CchPressureLoad& pressure,
    Interface::HermiteGalsIdentity identity,Interface::HermiteGalsTimeLayer clock,
    const CchAdvectorReconstruction& source,const CchAdvectorBoundary* boundary) const
{
    validateFor(mesh,pressure,identity,clock);require(data_->capillary,"constant load has no source reconstruction");
    require(reconstruction().sameIdentity(source),"different source reconstruction owner");
    require((this->boundary()==nullptr)==(boundary==nullptr),"different source boundary owner");
    if(boundary)require(this->boundary()->sameIdentity(*boundary),"different source boundary owner");
}
}
