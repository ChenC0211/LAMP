#include "Modules/Hydrodynamics/HydroInterfaceNormalLoad.h"
#include "Modules/Hydrodynamics/HydroSpatialAdvector.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <algorithm>
#include <cmath>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>

namespace LAMP::Hydrodynamics
{
namespace
{
void require(bool value,const char* reason)
{
    if(!value)throw std::invalid_argument(std::string("[HydroInterfaceNormalLoad] ")+reason);
}
}
namespace
{
using ScalarJet=LevelSetHermiteGeometry::Snapshot::ScalarJet2;
bool contains(const std::array<Vector2D,2>& box,Vector2D point)
{return point.x>=box[0].x&&point.x<=box[1].x&&point.y>=box[0].y&&point.y<=box[1].y;}
struct LinkSample {Vector2D point;ScalarJet jet;};
LinkSample interfaceOnLink(const LevelSetHermiteGeometry::Snapshot& scalar,
    int negative,int positive,Vector2D begin,Vector2D end,ScalarJet first,ScalarJet last)
{
    const auto negativeBounds=scalar.nativeCellBounds(negative),positiveBounds=scalar.nativeCellBounds(positive);
    const auto sample=[&](Vector2D q){
        // The two actual adjacent leaves cover the center link, including AMR.
        const int cell=contains(positiveBounds,q)?positive:negative;
        require(contains(cell==positive?positiveBounds:negativeBounds,q),"interface link leaves native support");
        const auto jet=scalar.scalarJetInNativeCell(cell,q);
        require(std::isfinite(jet.value),"nonfinite interface link scalar");return LinkSample{q,jet};
    };
    if(first.value==0)return {begin,first};
    if(last.value==0)return {end,last};
    require((first.value<0)!=(last.value<0),"interface link lacks phase bracket");
    LinkSample best=std::abs(first.value)<std::abs(last.value)?LinkSample{begin,first}:LinkSample{end,last};
    for(int iteration=0;iteration<52;++iteration){
        const auto middle=begin+(end-begin)*.5;
        if((middle.x==begin.x&&middle.y==begin.y)||(middle.x==end.x&&middle.y==end.y))break;
        const auto current=sample(middle);
        if(std::abs(current.jet.value)<std::abs(best.jet.value))best=current;
        if(current.jet.value==0)return current;
        if((current.jet.value<0)==(first.value<0)){begin=middle;first=current.jet;}
        else{end=middle;last=current.jet;}
    }
    const long double gradient=std::hypot(static_cast<long double>(best.jet.dx),static_cast<long double>(best.jet.dy));
    const double h=std::min({negativeBounds[1].x-negativeBounds[0].x,negativeBounds[1].y-negativeBounds[0].y,
        positiveBounds[1].x-positiveBounds[0].x,positiveBounds[1].y-positiveBounds[0].y});
    require(std::isfinite(gradient)&&gradient>0&&std::abs(best.jet.value)/gradient<=1e-10*h,
        "interface link residual exceeds tolerance");
    return best;
}
double physicalCurvature(const LinkSample& sample,bool radial)
{
    const auto& jet=sample.jet;
    const long double gradient=std::hypot(static_cast<long double>(jet.dx),static_cast<long double>(jet.dy));
    require(std::isfinite(gradient)&&gradient>0,"unresolved interface link gradient");
    const long double nx=jet.dx/gradient,ny=jet.dy/gradient;
    const long double meridional=(ny*ny*jet.dxx-2*nx*ny*jet.dxy+nx*nx*jet.dyy)/gradient;
    const long double azimuthal=radial?(sample.point.x==0?meridional:nx/sample.point.x):0;
    require(!radial||sample.point.x>0||(sample.point.x==0&&jet.dx==0),"invalid radial interface link");
    const double curvature=static_cast<double>(meridional+azimuthal);
    require(std::isfinite(curvature),"nonfinite physical curvature");return curvature;
}
}
struct HydroInterfaceNormalLoad::Data
{
    Data(HydroPressureLoad p,HydroPhaseGeometry g,double j)
        :pressure(std::move(p)),geometry(std::move(g)),jump(j){}
    HydroPressureLoad pressure;
    HydroPhaseGeometry geometry;
    double jump;
    bool capillary=false;
    double sigma=0;
    std::vector<double> curvature,measure,phaseIndicator;
    std::vector<double> load;
    HydroPressureLoad::Matrix stress;
    std::vector<double> reactions;
    std::vector<double> tangential;
    std::optional<HydroAdvectorReconstruction> source;
    std::optional<HydroAdvectorBoundary> boundary;
};
HydroInterfaceNormalLoad::HydroInterfaceNormalLoad(std::shared_ptr<const Data> data):data_(std::move(data)){}
HydroInterfaceNormalLoad HydroInterfaceNormalLoad::constant(const QuadtreeMesh& mesh,
    HydroPressureLoad pressure,HydroPhaseGeometry geometry,double jump)
{
    require(std::isfinite(jump),"nonfinite constant normal stress");
    require(pressure.layout().sameIdentity(geometry.layout()),"different pressure/geometry layout");
    require(&pressure.divergence().matrix()==&geometry.metrics().divergence().matrix(),
        "different authoritative divergence owner");
    geometry.validateFor(mesh,geometry.geometryIdentity(),geometry.timeLayer());
    auto data=std::make_shared<Data>(std::move(pressure),std::move(geometry),jump);
    std::vector<double> potential=data->geometry.liquidFractions();
    for(double& value:potential){value*=jump;require(std::isfinite(value),"nonfinite normal stress potential");}
    data->load=data->pressure.apply(mesh,HydroNativeValues(data->pressure.layout(),
        std::vector<double>(data->pressure.layout().velocityCount(),0),potential));
    const int n=data->pressure.layout().velocityCount();data->stress.resize(n,n);
    for(int u=0;u<n;++u)if(jump!=0)data->stress.insert(u,u)=jump;data->stress.makeCompressed();
    return HydroInterfaceNormalLoad(std::move(data));
}
HydroInterfaceNormalLoad HydroInterfaceNormalLoad::capillary(const QuadtreeMesh& mesh,
    HydroPressureLoad pressure,HydroPhaseGeometry geometry,double sigma,
    const HydroAdvectorReconstruction& source,const HydroAdvectorBoundary* boundary)
{
    require(std::isfinite(sigma) && sigma>=0,"invalid surface tension");
    require(pressure.layout().sameIdentity(geometry.layout()),"different pressure/geometry layout");
    require(&pressure.divergence().matrix()==&geometry.metrics().divergence().matrix(),
        "different authoritative divergence owner");
    geometry.validateFor(mesh,geometry.geometryIdentity(),geometry.timeLayer());
    require(source.momentumGeometry().layout().sameIdentity(pressure.layout()),"different source layout");
    require(&source.momentumGeometry().metrics().divergence().matrix()==&pressure.divergence().matrix(),"different source divergence owner");
    if(boundary)boundary->validateFor(mesh,source.momentumGeometry());
    auto data=std::make_shared<Data>(std::move(pressure),std::move(geometry),0);
    data->capillary=true;data->sigma=sigma;
    data->source=source;if(boundary)data->boundary=*boundary;
    if(const auto* bound=source.phaseGeometry())bound->validateFor(mesh,data->geometry.geometryIdentity(),data->geometry.timeLayer());
    data->tangential.assign(source.tangentialModes().size(),0.);
    const auto& scalar=data->geometry.publication().state().geometry().snapshot();
    const auto& layout=data->pressure.layout();const int count=layout.velocityCount();
    std::vector<Vector2D> centers(layout.pressureCount());std::vector<ScalarJet> jets(layout.pressureCount());
    data->phaseIndicator.resize(layout.pressureCount());
    for(int p=0;p<layout.pressureCount();++p){
        const int native=layout.nativeCell(p);centers[p]=mesh.cells(mesh.cellIndexToLabel[native]).centerCoord();
        jets[p]=scalar.scalarJetInNativeCell(native,centers[p]);const double value=jets[p].value;
        require(std::isfinite(value),"nonfinite pressure-center scalar");
        data->phaseIndicator[p]=value<0?1:value>0?0:.5;
    }
    const auto incidence=data->pressure.apply(mesh,HydroNativeValues(layout,std::vector<double>(count),data->phaseIndicator));
    data->load.assign(count,0);data->curvature.assign(count,0);data->measure.assign(count,0);
    data->reactions.assign(count,0);data->stress.resize(count,count);
    for(int u=0;u<count;++u){
        const auto face=layout.domain().face(mesh,layout.nativeFace(u));
        // Prescribed outer faces and the zero-measure axis are reactions, not
        // internal capillary links. Same-phase links carry exactly no force.
        if(face.negativeCell<0||face.positiveCell<0||incidence[u]==0)continue;
        const int a=layout.pressureDof(face.negativeCell),b=layout.pressureDof(face.positiveCell);
        require(a>=0&&b>=0,"interface link has unselected pressure center");
        const auto sample=interfaceOnLink(scalar,face.negativeCell,face.positiveCell,centers[a],centers[b],jets[a],jets[b]);
        const double curvature=physicalCurvature(sample,layout.domain().axisymmetric()),coefficient=-sigma*curvature;
        data->curvature[u]=curvature;data->measure[u]=std::hypot(centers[b].x-centers[a].x,centers[b].y-centers[a].y);
        data->load[u]=coefficient*incidence[u];
        require(std::isfinite(coefficient)&&std::isfinite(data->load[u]),"nonfinite capillary load");
        if(coefficient!=0)data->stress.insert(u,u)=coefficient;
    }
    data->stress.makeCompressed();
    return HydroInterfaceNormalLoad(std::move(data));
}
const HydroPressureLoad& HydroInterfaceNormalLoad::pressureLoad() const noexcept{return data_->pressure;}
const HydroPhaseGeometry& HydroInterfaceNormalLoad::geometry() const noexcept{return data_->geometry;}
const HydroDofLayout& HydroInterfaceNormalLoad::layout() const noexcept{return data_->pressure.layout();}
bool HydroInterfaceNormalLoad::isCapillary() const noexcept{return data_->capillary;}
double HydroInterfaceNormalLoad::constantJump() const{require(!data_->capillary,"capillary load has no prescribed constant jump");return data_->jump;}
double HydroInterfaceNormalLoad::surfaceTension() const{require(data_->capillary,"constant load has no surface tension");return data_->sigma;}
const std::vector<double>& HydroInterfaceNormalLoad::phaseIndicators() const{require(data_->capillary,"constant load has no pressure-center phase");return data_->phaseIndicator;}
const std::vector<double>& HydroInterfaceNormalLoad::linkCurvatures() const{require(data_->capillary,"constant load has no curvature links");return data_->curvature;}
const std::vector<double>& HydroInterfaceNormalLoad::linkLengths() const{require(data_->capillary,"constant load has no surface supports");return data_->measure;}
const std::vector<double>& HydroInterfaceNormalLoad::balanceLoads() const noexcept{return data_->load;}
const std::vector<double>& HydroInterfaceNormalLoad::tangentialBalanceLoads() const
{require(data_->capillary,"constant load has no tangential source binding");return data_->tangential;}
const HydroPressureLoad::Matrix& HydroInterfaceNormalLoad::stressMatrix() const noexcept{return data_->stress;}
const HydroAdvectorReconstruction& HydroInterfaceNormalLoad::reconstruction() const{require(data_->capillary,"constant load has no source reconstruction");return *data_->source;}
const HydroAdvectorBoundary* HydroInterfaceNormalLoad::boundary() const{require(data_->capillary,"constant load has no source boundary");return data_->boundary?&*data_->boundary:nullptr;}
const std::vector<double>& HydroInterfaceNormalLoad::constantReactions() const{require(data_->capillary,"constant load has no bound wall reactions");return data_->reactions;}
void HydroInterfaceNormalLoad::validateFor(const QuadtreeMesh& mesh,const HydroPressureLoad& pressure,
    Interface::HermiteGalsIdentity identity,Interface::HermiteGalsTimeLayer clock) const
{
    require(&pressure.matrix()==&data_->pressure.matrix(),"different authoritative pressure owner");
    data_->geometry.validateFor(mesh,identity,clock);
}
void HydroInterfaceNormalLoad::validateFor(const QuadtreeMesh& mesh,const HydroPressureLoad& pressure,
    Interface::HermiteGalsIdentity identity,Interface::HermiteGalsTimeLayer clock,
    const HydroAdvectorReconstruction& source,const HydroAdvectorBoundary* boundary) const
{
    validateFor(mesh,pressure,identity,clock);require(data_->capillary,"constant load has no source reconstruction");
    require(reconstruction().sameIdentity(source),"different source reconstruction owner");
    require((this->boundary()==nullptr)==(boundary==nullptr),"different source boundary owner");
    if(boundary)require(this->boundary()->sameIdentity(*boundary),"different source boundary owner");
}
}
