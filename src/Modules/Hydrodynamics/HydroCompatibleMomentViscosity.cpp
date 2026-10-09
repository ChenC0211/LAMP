#include "Modules/Hydrodynamics/HydroCompatibleMomentViscosity.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/NativeLeafDomain.h"
#include <cmath>
#include <stdexcept>
#include <limits>

namespace Quantum::Hydrodynamics {
namespace {
void require(bool ok,const char* message){if(!ok)throw std::invalid_argument(message);}
constexpr double gx[3]={-.774596669241483377,0,.774596669241483377};
constexpr double gw[3]={5./9,8./9,5./9};
}
HydroCompatibleMomentViscosity HydroCompatibleMomentViscosity::build(const QuadtreeMesh& mesh,
    const NativeLeafDomain& domain,const HydroCompatibleMomentSpace& space,const std::vector<Eigen::Matrix3d>& mu){
    space.validate(mesh,domain);
    const int nc=mesh.cellIndexNum(),nf=domain.nativeFaceCount();
    require(mu.size()==std::size_t(nc),"moment viscosity needs one owned Gram per native cell");
    for(const auto& value:mu){
        require(value.allFinite()&&(value-value.transpose()).norm()<=256*std::numeric_limits<double>::epsilon()*value.norm(),"invalid viscosity moment symmetry");
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eig(value);
        require(eig.info()==Eigen::Success&&eig.eigenvalues()[0]>=0,"viscosity moments are not a nonnegative measure");
    }
    if(domain.axisymmetric())return buildAxisymmetric(mesh,domain,space,mu);
    HydroCompatibleMomentViscosity result;result.identity_=domain.meshIdentity();result.bulkSize_=int(space.mass().rows());
    result.tangentRows_.resize(nf);for(int f=0;f<nf;f++)result.tangentRows_[f]=result.bulkSize_+2*f;
    result.cells_.resize(nc);std::vector<std::vector<int>> incident(nc);
    for(const auto face:domain.faces(mesh)){
        if(face.negativeCell>=0)incident[face.negativeCell].push_back(face.nativeRow);
        if(face.positiveCell>=0)incident[face.positiveCell].push_back(face.nativeRow);
    }
    std::vector<Eigen::Triplet<double>> entries;
    const double rootTwo=std::sqrt(2.);
    for(int c=0;c<nc;c++){
        const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);const auto center=cell.centerCoord(),h=cell.size();
        const double volume=h.x*h.y,diameter=std::hypot(h.x,h.y);
        const auto& bulk=space.cells()[c];auto& block=result.cells_[c];block.rows=bulk.rows;
        const int count=int(incident[c].size()),nb=int(bulk.rows.size()),n=nb+2*count;
        for(int f:incident[c]){block.rows.push_back(result.firstTangent(f));block.rows.push_back(result.firstTangent(f)+1);}
        // All strain functionals use the same integration-by-parts identity.
        // For linear stress tests, div(stress) is constant, so cell means suffice.
        Eigen::MatrixXd work=Eigen::MatrixXd::Zero(9,n);
        Eigen::MatrixXd faceSamples=Eigen::MatrixXd::Zero(4*count,12),trace=Eigen::MatrixXd::Zero(4*count,n);
        Eigen::VectorXd faceMetric(4*count);
        work(1,nb-3)=-volume/h.x;work(5,nb-2)=-volume/h.y;
        work(7,nb-2)=-volume/(rootTwo*h.x);work(8,nb-3)=-volume/(rootTwo*h.y);
        for(int i=0;i<count;i++){
            const auto face=domain.face(mesh,incident[c][i]);const auto& info=space.faces()[face.nativeRow];
            const auto& edge=face.horizontal?mesh.hEdges(face.edgeLabel):mesh.vEdges(face.edgeLabel);
            const auto a=mesh.verts(edge.edgeToVert[0]).coord(),b=mesh.verts(edge.edgeToVert[1]).coord();
            const double sign=face.negativeCell==c?1.:-1.;
            for(int m=0;m<2;m++){
                trace(4*i+m,2*i+m)=1;trace(4*i+2+m,nb+2*i+m)=1;
                faceMetric[4*i+m]=faceMetric[4*i+2+m]=face.area/diameter;
            }
            for(int g=0;g<3;g++){
                const double t=.5*(1+gx[g]),ds=.5*face.area*gw[g];
                const Vector2D p{a.x+t*(b.x-a.x),a.y+t*(b.y-a.y)},local{(p.x-center.x)/h.x,(p.y-center.y)/h.y};
                const double theta=((face.horizontal?p.x:p.y)-info.tangentCenter)/info.tangentScale;
                Eigen::VectorXd normal=Eigen::VectorXd::Zero(n),tangent=normal;
                normal[2*i]=1;normal[2*i+1]=theta;tangent[nb+2*i]=1;tangent[nb+2*i+1]=theta;
                Eigen::Vector3d q;q<<1,local.x,local.y;
                const auto& ux=face.horizontal?tangent:normal;const auto& uy=face.horizontal?normal:tangent;
                if(face.horizontal){work.block(3,0,3,n)+=sign*ds*q*uy.transpose();work.block(6,0,3,n)+=sign*ds/rootTwo*q*ux.transpose();}
                else {work.block(0,0,3,n)+=sign*ds*q*ux.transpose();work.block(6,0,3,n)+=sign*ds/rootTwo*q*uy.transpose();}
                Eigen::Matrix<double,6,1> polynomial;polynomial<<1,local.x,local.y,local.x*local.x,local.x*local.y,local.y*local.y;
                for(int m=0;m<2;m++){
                    const double dw=ds/face.area*(m?theta:1.);
                    faceSamples.block(4*i+m,face.horizontal?6:0,1,6)+=dw*polynomial.transpose();
                    faceSamples.block(4*i+2+m,face.horizontal?0:6,1,6)+=dw*polynomial.transpose();
                }
            }
        }
        Eigen::VectorXd gram(9);for(int k=0;k<3;k++)gram.segment<3>(3*k)<<volume,volume/12,volume/12;
        block.strain=gram.cwiseInverse().asDiagonal()*work;
        // Every linear symmetric strain in 2D is the strain of a quadratic
        // velocity, unique up to a rigid motion. Fix that motion with the same
        // two means and rotation moment already owned by the bulk space.
        Eigen::Matrix<double,12,12> constraints=Eigen::Matrix<double,12,12>::Zero();
        constraints(0,1)=1/h.x;constraints(1,3)=2/h.x;constraints(2,4)=1/h.x;
        constraints(3,8)=1/h.y;constraints(4,10)=1/h.y;constraints(5,11)=2/h.y;
        constraints(6,2)=1/(rootTwo*h.y);constraints(6,7)=1/(rootTwo*h.x);
        constraints(7,4)=1/(rootTwo*h.y);constraints(7,9)=rootTwo/h.x;
        constraints(8,5)=rootTwo/h.y;constraints(8,10)=1/(rootTwo*h.x);
        constraints(9,0)=1;constraints(9,3)=constraints(9,5)=1./12;
        constraints(10,6)=1;constraints(10,9)=constraints(10,11)=1./12;
        constraints(11,2)=-1./12;constraints(11,7)=1./12;
        Eigen::MatrixXd rhs=Eigen::MatrixXd::Zero(12,n);rhs.topRows(9)=block.strain;
        for(int j=0;j<3;j++)rhs(9+j,nb-3+j)=1;
        block.reconstruction=constraints.partialPivLu().solve(rhs);
        const Eigen::MatrixXd difference=trace-faceSamples*block.reconstruction;
        Eigen::Matrix<double,9,9> material=Eigen::Matrix<double,9,9>::Zero();
        for(int j=0;j<3;j++)material.block<3,3>(3*j,3*j)=mu[c];
        block.dissipation=2*(block.strain.transpose()*material*block.strain
            +mu[c](0,0)/volume*difference.transpose()*faceMetric.asDiagonal()*difference);
        for(int i=0;i<n;i++)for(int j=0;j<n;j++)if(block.dissipation(i,j)!=0)
            entries.emplace_back(block.rows[i],block.rows[j],block.dissipation(i,j));
    }
    result.matrix_.resize(result.bulkSize_+2*nf,result.bulkSize_+2*nf);result.matrix_.setFromTriplets(entries.begin(),entries.end());
    return result;
}
HydroCompatibleMomentViscosity HydroCompatibleMomentViscosity::buildCellwise(const QuadtreeMesh& mesh,
    const NativeLeafDomain& domain,const HydroCompatibleMomentSpace& space,const std::vector<double>& viscosity){
    space.validate(mesh,domain);
    require(viscosity.size()==std::size_t(mesh.cellIndexNum()),"moment viscosity cell count differs");
    std::vector<Eigen::Matrix3d> moments(viscosity.size(),Eigen::Matrix3d::Zero());
    for(int c=0;c<int(viscosity.size());c++){
        require(std::isfinite(viscosity[c])&&viscosity[c]>=0,"invalid piecewise constant viscosity");
        const double integral=viscosity[c]*domain.cellVolumes()[domain.cellRowByIndex()[c]];
        moments[c].diagonal()<<integral,integral/12,integral/12;
        if(domain.axisymmetric()){
            const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);
            moments[c](0,1)=moments[c](1,0)=integral*cell.size().x/(12*cell.centerCoord().x);
        }
    }
    auto result=build(mesh,domain,space,moments);
    result.cellViscosity_=std::make_shared<const std::vector<double>>(viscosity);return result;
}
void HydroCompatibleMomentViscosity::validate(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,const HydroCompatibleMomentSpace& space) const{
    space.validate(mesh,domain);require(identity_&&identity_==domain.meshIdentity()
        &&bulkSize_==space.mass().rows(),"moment viscosity/domain source differs");
}
}
