#include "Modules/Hydrodynamics/HydroCompatibleMomentSpace.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/NativeLeafDomain.h"
#include "Core/Mesh/CartesianMeasure.h"
#include "Core/Geometry/ParametricHermiteCurve.h"
#include <stdexcept>

namespace Quantum::Hydrodynamics {
namespace {
double weight(Vector2D p,bool rz){return rz?CartesianMeasure::twoPi*p.x:1.;}
Eigen::Vector2d value(const HydroCompatibleMomentSpace::PhysicalVector& input,Vector2D p){
    if(!input)throw std::invalid_argument("missing physical moment field");const Eigen::Vector2d result=input(p);
    if(!result.allFinite())throw std::invalid_argument("nonfinite physical moment field");return result;
}
template<class F>void integrateCell(Vector2D center,Vector2D h,F&& f){
    const auto& g=Quantum::Geometry::Parametric::gaussEight();
    for(int j=0;j<8;j++)for(int i=0;i<8;i++)f(Vector2D{center.x+.5*h.x*g.x[i],center.y+.5*h.y*g.x[j]},.25*h.x*h.y*g.w[i]*g.w[j]);
}
std::vector<double> splitMomentInterval(double a,double b,const std::vector<double>& planes){
    std::vector<double> out{a};for(double p:planes)if(p>a&&p<b)out.push_back(p);out.push_back(b);
    std::sort(out.begin(),out.end());out.erase(std::unique(out.begin(),out.end()),out.end());return out;
}
}
Eigen::VectorXd HydroCompatibleMomentSpace::interpolateVelocity(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,const PhysicalVector& velocity) const{
    return interpolateVelocity(mesh,domain,velocity,{});
}
Eigen::VectorXd HydroCompatibleMomentSpace::interpolateVelocity(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,const PhysicalVector& velocity,const Partition& planes) const{
    validate(mesh,domain);Eigen::VectorXd result=Eigen::VectorXd::Zero(mass_.rows());const auto& g=Quantum::Geometry::Parametric::gaussEight();
    for(const auto face:domain.faces(mesh)){
        const auto& info=faces_[face.nativeRow];if(info.first<0)continue;
        const auto& edge=face.horizontal?mesh.hEdges(face.edgeLabel):mesh.vEdges(face.edgeLabel);
        const auto a=mesh.verts(edge.edgeToVert[0]).coord(),b=mesh.verts(edge.edgeToVert[1]).coord();const double length=face.horizontal?b.x-a.x:b.y-a.y;
        const auto cuts=splitMomentInterval(0,1,[&]{std::vector<double> cuts;for(double p:planes[face.horizontal?0:1])cuts.push_back((p-(face.horizontal?a.x:a.y))/length);return cuts;}());
        for(std::size_t part=1;part<cuts.size();part++)for(int i=0;i<8;i++){
            const double t=cuts[part-1]+(cuts[part]-cuts[part-1])*.5*(1+g.x[i]);const Vector2D p{a.x+t*(b.x-a.x),a.y+t*(b.y-a.y)};
            const double u=value(velocity,p)[face.horizontal?1:0],dw=.5*length*(cuts[part]-cuts[part-1])*g.w[i]*weight(p,axisymmetric_)/info.area;
            result[info.first]+=dw*u;result[info.first+1]+=dw*u*((face.horizontal?p.x:p.y)-info.tangentCenter)/info.tangentScale;
        }
    }
    for(int c=0;c<mesh.cellIndexNum();c++){
        const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);const auto center=cell.centerCoord(),h=cell.size();const int first=firstCellVelocity(c);
        const double volume=domain.cellVolumes()[domain.cellRowByIndex()[c]];
        const auto add=[&](Vector2D p,double da){const auto u=value(velocity,p);const double dv=da*weight(p,axisymmetric_)/volume;
            const double x=(p.x-center.x)/h.x,y=(p.y-center.y)/h.y;
            result[first]+=dv*u[0];result[first+1]+=dv*u[1];result[first+2]+=dv*(-y*u[0]+x*u[1]);
            };
        if(planes[0].empty()&&planes[1].empty())integrateCell(center,h,add);
        else {const auto xs=splitMomentInterval(center.x-.5*h.x,center.x+.5*h.x,planes[0]),ys=splitMomentInterval(center.y-.5*h.y,center.y+.5*h.y,planes[1]);
            for(std::size_t j=1;j<ys.size();j++)for(std::size_t i=1;i<xs.size();i++)
                integrateCell({.5*(xs[i-1]+xs[i]),.5*(ys[j-1]+ys[j])},{xs[i]-xs[i-1],ys[j]-ys[j-1]},add);}
    }
    if(axisymmetric_){
        for(int c=0;c<mesh.cellIndexNum();c++){
            const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);const auto center=cell.centerCoord(),h=cell.size();
            result.segment<3>(firstCellVelocity(c)+3)=rzExtraMoments(center,h,[&](Vector2D p)->Eigen::MatrixXd{return value(velocity,p);},planes);
        }
    }
    return result;
}
Eigen::VectorXd HydroCompatibleMomentSpace::integrateBodyForce(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,const PhysicalVector& force) const{
    validate(mesh,domain);Eigen::VectorXd result=Eigen::VectorXd::Zero(mass_.rows());
    for(int c=0;c<mesh.cellIndexNum();c++){
        const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);const auto center=cell.centerCoord(),h=cell.size();const auto& block=cells_[c];
        Eigen::Matrix<double,6,1> integral=Eigen::Matrix<double,6,1>::Zero();
        integrateCell(center,h,[&](Vector2D p,double da){integral+=da*weight(p,axisymmetric_)*affineVelocityBasis({(p.x-center.x)/h.x,(p.y-center.y)/h.y},h).transpose()*value(force,p);});
        const Eigen::VectorXd local=block.affineProjection.transpose()*integral;for(int j=0;j<local.size();j++)result[block.rows[j]]+=local[j];
    }
    return result;
}
}
