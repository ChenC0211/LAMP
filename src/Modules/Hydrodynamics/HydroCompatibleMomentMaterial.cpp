#include "Modules/Hydrodynamics/HydroCompatibleMomentMaterial.h"
#include "Core/Mesh/NativeLeafDomain.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace Quantum::Hydrodynamics {
namespace {
void require(bool ok,const char* message){if(!ok)throw std::invalid_argument(message);}
Eigen::Matrix3d gram(const std::array<double,6>& m){
    Eigen::Matrix3d result;result<<m[0],m[1],m[2],m[1],m[3],m[4],m[2],m[4],m[5];return result;
}
Eigen::Matrix<double,6,6> vectorGram(const Eigen::Matrix3d& measure,Vector2D h){
    const auto b0=HydroCompatibleMomentSpace::affineVelocityBasis({0,0},h);
    const auto bx=(HydroCompatibleMomentSpace::affineVelocityBasis({1,0},h)-b0).eval();
    const auto by=(HydroCompatibleMomentSpace::affineVelocityBasis({0,1},h)-b0).eval();
    Eigen::Matrix<double,6,6> result=Eigen::Matrix<double,6,6>::Zero();
    for(int component=0;component<2;component++){
        Eigen::Matrix<double,3,6> b;b.row(0)=b0.row(component);b.row(1)=bx.row(component);b.row(2)=by.row(component);
        result+=b.transpose()*measure*b;
    }
    return result;
}
}
HydroCompatibleMomentMaterial HydroCompatibleMomentMaterial::buildFrozen(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,
    const HydroCompatibleMomentSpace& space,std::shared_ptr<const Curve> geometry,double time,Properties properties){
    using Boundary=Geometry::Parametric::Boundary;
    return buildFrozen(mesh,domain,space,std::make_shared<const Boundary>(std::vector<Boundary::Component>{{1,std::move(geometry)}}),time,properties);
}
HydroCompatibleMomentMaterial HydroCompatibleMomentMaterial::buildFrozen(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,
    const HydroCompatibleMomentSpace& space,std::shared_ptr<const Geometry::Parametric::Boundary> geometry,double time,Properties properties){
    space.validate(mesh,domain);
    require(geometry&&std::isfinite(time),"invalid frozen material geometry/time");
    for(double rho:{properties.liquidDensity,properties.gasDensity})require(std::isfinite(rho)&&rho>0,"phase density must be positive");
    for(double mu:{properties.liquidViscosity,properties.gasViscosity})require(std::isfinite(mu)&&mu>=0,"phase viscosity must be nonnegative");
    HydroCompatibleMomentMaterial result;result.identity_=domain.meshIdentity();result.boundary_=std::move(geometry);result.time_=time;result.properties_=properties;
    if(properties.liquidDensity==properties.gasDensity)result.uniformDensity_=properties.liquidDensity;
    if(result.boundary_->components().empty())result.uniformDensity_=properties.gasDensity;
    std::vector<Geometry::Parametric::PhaseMeasure> measures;
    for(const auto& component:result.boundary_->components())measures.emplace_back(*component.curve);
    result.cells_.resize(mesh.cellIndexNum());std::vector<Eigen::Triplet<double>> entries;
    for(int c=0;c<mesh.cellIndexNum();c++){
        const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);const auto h=cell.size(),x=cell.centerCoord();
        const double volume=domain.cellVolumes()[domain.cellRowByIndex()[c]];
        auto& saved=result.cells_[c];saved.volume=Eigen::Matrix3d::Zero();saved.volume.diagonal()<<volume,volume/12,volume/12;
        if(domain.axisymmetric())saved.volume(0,1)=saved.volume(1,0)=volume*h.x/(12*x.x);
        saved.liquid=Eigen::Matrix3d::Zero();
        for(auto& measure:measures){const auto phase=measure.rectanglePhaseMoments(x.x-.5*h.x,x.x+.5*h.x,x.y-.5*h.y,x.y+.5*h.y,domain.axisymmetric());
            saved.componentLiquid.push_back(phase.phase==Geometry::Parametric::RectanglePhase::Full?saved.volume:gram(phase.values));saved.liquid+=saved.componentLiquid.back();}
        // Validate the phase and its complement; never clip eigenvalues or
        // turn a geometric defect into a change of conserved phase volume.
        for(const Eigen::Matrix3d part:{saved.liquid,Eigen::Matrix3d(saved.volume-saved.liquid)}){
            Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigen(part);
            require(part.allFinite()&&eigen.info()==Eigen::Success&&eigen.eigenvalues()[0]>=-4096*std::numeric_limits<double>::epsilon()*volume,
                "phase material integral is not a positive submeasure");
        }
        const Eigen::Matrix3d gas=saved.volume-saved.liquid;
        saved.density=properties.gasDensity*gas+properties.liquidDensity*saved.liquid;
        saved.viscosity=properties.gasViscosity*gas+properties.liquidViscosity*saved.liquid;
        const auto& block=space.cells()[c];const auto& projection=block.affineProjection;
        // R kills every affine velocity. The existing positive H0 supplies
        // its complement norm, weighted by the actual cell-average density.
        // Constant material is exactly the existing rho*H0 construction.
        Eigen::MatrixXd mass;
        if(result.uniformDensity_)mass=*result.uniformDensity_*block.mass;
        else {
            const Eigen::MatrixXd residual=Eigen::MatrixXd::Identity(block.rows.size(),block.rows.size())-block.affineSamples*projection;
            mass=projection.transpose()*vectorGram(saved.density,h)*projection
                +(saved.density(0,0)/volume)*residual.transpose()*block.mass*residual;
            mass=.5*(mass+mass.transpose()).eval();
        }
        require(mass.allFinite(),"nonfinite material kinetic mass");
        for(int i=0;i<mass.rows();i++)for(int j=0;j<mass.cols();j++)if(mass(i,j)!=0)entries.emplace_back(block.rows[i],block.rows[j],mass(i,j));
    }
    result.mass_.resize(space.mass().rows(),space.mass().cols());result.mass_.setFromTriplets(entries.begin(),entries.end());return result;
}
void HydroCompatibleMomentMaterial::validate(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,const HydroCompatibleMomentSpace& space) const{
    space.validate(mesh,domain);require(identity_&&identity_==domain.meshIdentity(),"material geometry/mesh source differs");
    require(mass_.rows()==space.mass().rows(),"material velocity space differs");
}
HydroCompatibleMomentMaterial HydroCompatibleMomentMaterial::buildSinglePhase(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,
    const HydroCompatibleMomentSpace& space,double density,const std::vector<double>& viscosity,double time){
    ViscosityField field;field.reserve(viscosity.size());
    for(double mu:viscosity)field.push_back(Eigen::Matrix3d::Constant(mu));
    return buildSinglePhase(mesh,domain,space,density,field,time);
}
HydroCompatibleMomentMaterial HydroCompatibleMomentMaterial::buildSinglePhase(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,
    const HydroCompatibleMomentSpace& space,double density,const ViscosityField& viscosity,double time){
    require(viscosity.size()==std::size_t(mesh.cellIndexNum()),"prescribed viscosity cell count differs");
    double bound=0;for(const auto& coefficients:viscosity){
        require(coefficients.allFinite()&&coefficients.minCoeff()>=0,"prescribed viscosity needs finite nonnegative Bernstein coefficients");
        bound=std::max(bound,coefficients.maxCoeff());
    }
    using Boundary=Geometry::Parametric::Boundary;
    auto result=buildFrozen(mesh,domain,space,std::make_shared<const Boundary>(std::vector<Boundary::Component>{}),time,{density,density,0,0});
    result.viscosity_=std::make_shared<const ViscosityField>(viscosity);result.viscosityBound_=bound;
    // mu is degree (2,2). Multiplication by the quadratic moment products
    // and the RZ radial measure gives degree at most (5,4), exactly covered
    // by this positive 3x3 Gauss rule. No recovery from the old Gram is used.
    const double root=std::sqrt(3./5),nodes[3]={.5*(1-root),.5,.5*(1+root)},weights[3]={5./18,4./9,5./18};
    for(int c=0;c<mesh.cellIndexNum();c++){
        const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);const auto center=cell.centerCoord(),h=cell.size();
        const double volume=domain.cellVolumes()[domain.cellRowByIndex()[c]];
        auto& moment=result.cells_[c].viscosity;moment.setZero();
        for(int i=0;i<3;i++)for(int j=0;j<3;j++){
            const Vector2D unit{nodes[i],nodes[j]};const Eigen::Vector3d shape(1,unit.x-.5,unit.y-.5);
            const double radial=domain.axisymmetric()?1+h.x*(unit.x-.5)/center.x:1;
            moment.noalias()+=volume*weights[i]*weights[j]*radial*result.viscosityAt(c,unit,false)*shape*shape.transpose();
        }
        require(moment.allFinite(),"nonfinite prescribed viscosity integral");
    }
    return result;
}
double HydroCompatibleMomentMaterial::viscosityUpperBound() const noexcept{
    return viscosity_?viscosityBound_:std::max(properties_.liquidViscosity,properties_.gasViscosity);
}
HydroCompatibleMomentMaterial HydroCompatibleMomentMaterial::atTime(double time) const{
    require(identity_&&boundary_&&std::isfinite(time),"invalid fixed-material evaluation time");
    auto result=*this;result.time_=time;return result;
}
double HydroCompatibleMomentMaterial::viscosityAt(int cell,Vector2D unit,bool liquid) const{
    if(!viscosity_)return liquid?properties_.liquidViscosity:properties_.gasViscosity;
    require(std::isfinite(unit.x)&&std::isfinite(unit.y)&&unit.x>=0&&unit.x<=1&&unit.y>=0&&unit.y<=1,
        "prescribed viscosity query outside its cell chart");
    const Eigen::Vector3d bx((1-unit.x)*(1-unit.x),2*unit.x*(1-unit.x),unit.x*unit.x);
    const Eigen::Vector3d by((1-unit.y)*(1-unit.y),2*unit.y*(1-unit.y),unit.y*unit.y);
    return bx.dot(viscosity_->at(cell)*by);
}
std::vector<Eigen::Matrix3d> HydroCompatibleMomentMaterial::viscosityMoments() const{
    std::vector<Eigen::Matrix3d> result;result.reserve(cells_.size());for(const auto& cell:cells_)result.push_back(cell.viscosity);return result;
}
}
