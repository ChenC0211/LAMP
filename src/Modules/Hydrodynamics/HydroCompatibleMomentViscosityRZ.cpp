#include "Modules/Hydrodynamics/HydroCompatibleMomentViscosity.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/NativeLeafDomain.h"
#include "Core/Mesh/CartesianMeasure.h"
#include "Core/Geometry/ParametricHermiteCurve.h"
#include <cmath>
#include <stdexcept>

namespace Quantum::Hydrodynamics {
namespace {
void require(bool ok,const char* message){if(!ok)throw std::invalid_argument(message);}
Eigen::Matrix<double,6,1> polynomial(double x,double y){Eigen::Matrix<double,6,1> p;p<<1,x,y,x*x,x*y,y*y;return p;}
}
HydroCompatibleMomentViscosity HydroCompatibleMomentViscosity::buildAxisymmetric(const QuadtreeMesh& mesh,
    const NativeLeafDomain& domain,const HydroCompatibleMomentSpace& space,const std::vector<Eigen::Matrix3d>& mu){
    require(space.cellVelocityCount()==6,"RZ weak strain requires independent stable moments");
    HydroCompatibleMomentViscosity result;result.identity_=domain.meshIdentity();result.bulkSize_=int(space.mass().rows());
    const int nc=mesh.cellIndexNum(),nf=domain.nativeFaceCount();result.cells_.resize(nc);result.tangentRows_.assign(nf,-1);
    int size=result.bulkSize_;std::vector<std::vector<int>> incident(nc);std::vector<bool> axisCell(nc,false);
    for(const auto face:domain.faces(mesh))if(face.area>0){
        result.tangentRows_[face.nativeRow]=size;size+=2;
        if(face.negativeCell>=0)incident[face.negativeCell].push_back(face.nativeRow);
        if(face.positiveCell>=0)incident[face.positiveCell].push_back(face.nativeRow);
    }else {if(face.negativeCell>=0)axisCell[face.negativeCell]=true;if(face.positiveCell>=0)axisCell[face.positiveCell]=true;}
    const auto& gauss=Quantum::Geometry::Parametric::gaussEight();const double twoPi=CartesianMeasure::twoPi,rootTwo=std::sqrt(2.);
    std::vector<Eigen::Triplet<double>> entries;
    for(int c=0;c<nc;c++){
        const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);const auto center=cell.centerCoord(),h=cell.size();
        const double flat=h.x*h.y,volume=domain.cellVolumes()[domain.cellRowByIndex()[c]],diameter=std::hypot(h.x,h.y);
        const auto& bulk=space.cells()[c];auto& block=result.cells_[c];block.rows=bulk.rows;
        const int count=int(incident[c].size()),nb=int(bulk.rows.size()),first=nb-6,n=nb+2*count;
        for(int face:incident[c]){block.rows.push_back(result.firstTangent(face));block.rows.push_back(result.firstTangent(face)+1);}
        // Flat velocity integrals against (1,X,Y). Use first moments directly:
        // mean_flat(u)=mean_r(u)-(hx/rc)*mean_flat(X*u), avoiding r/h subtraction.
        Eigen::MatrixXd radial=Eigen::MatrixXd::Zero(3,n),axial=radial;
        radial.leftCols(nb)=flat*bulk.flatVelocityMoments.topRows(3);
        axial.leftCols(nb)=flat*bulk.flatVelocityMoments.bottomRows(3);
        Eigen::MatrixXd work=Eigen::MatrixXd::Zero(12,n);
        work.topRows(3)=-twoPi*radial;work(1,first)-=volume/h.x;
        work(5,first+1)-=volume/h.y;
        work.middleRows(6,3)=twoPi*radial;
        work.bottomRows(3)=-twoPi/rootTwo*axial;
        work(10,first+1)-=volume/(rootTwo*h.x);work(11,first)-=volume/(rootTwo*h.y);
        Eigen::MatrixXd samples=Eigen::MatrixXd::Zero(n,12);
        Eigen::VectorXd metric=Eigen::VectorXd::Zero(n);
        for(int i=0;i<count;i++){
            const auto face=domain.face(mesh,incident[c][i]);const auto& info=space.faces()[face.nativeRow];
            const auto& edge=face.horizontal?mesh.hEdges(face.edgeLabel):mesh.vEdges(face.edgeLabel);
            const auto a=mesh.verts(edge.edgeToVert[0]).coord(),b=mesh.verts(edge.edgeToVert[1]).coord();
            const double sign=face.negativeCell==c?1.:-1.,length=face.horizontal?b.x-a.x:b.y-a.y;
            for(int m=0;m<2;m++)metric[2*i+m]=metric[nb+2*i+m]=face.area/diameter;
            for(int g=0;g<8;g++){
                const double t=.5*(1+gauss.x[g]);const Vector2D p{a.x+t*(b.x-a.x),a.y+t*(b.y-a.y)};
                const double x=(p.x-center.x)/h.x,y=(p.y-center.y)/h.y,dw=.5*length*gauss.w[g]*twoPi*p.x;
                const double theta=((face.horizontal?p.x:p.y)-info.tangentCenter)/info.tangentScale;
                const Eigen::Vector3d q(1,x,y);const auto poly=polynomial(x,y);
                for(int m=0;m<2;m++){
                    const double shape=m?theta:1.,weight=sign*dw*shape;
                    const int nr=2*i+m,tr=nb+2*i+m;
                    if(face.horizontal){work.block(3,nr,3,1)+=weight*q;work.block(9,tr,3,1)+=weight/rootTwo*q;}
                    else {work.block(0,nr,3,1)+=weight*q;work.block(9,tr,3,1)+=weight/rootTwo*q;}
                    samples.block(nr,face.horizontal?6:0,1,6)+=dw/face.area*shape*poly.transpose();
                    samples.block(tr,face.horizontal?0:6,1,6)+=dw/face.area*shape*poly.transpose();
                }
            }
        }
        Eigen::Matrix3d gram=Eigen::Matrix3d::Zero();gram.diagonal()<<volume,volume/12,volume/12;
        gram(0,1)=gram(1,0)=volume*h.x/(12*center.x);
        block.strain.resize(12,n);
        for(int component=0;component<4;component++)block.strain.middleRows(3*component,3)=gram.ldlt().solve(work.middleRows(3*component,3));
        for(int j=0;j<8;j++)for(int i=0;i<8;i++){
            const double x=.5*gauss.x[i],y=.5*gauss.x[j],r=center.x+h.x*x,da=.25*flat*gauss.w[i]*gauss.w[j],dv=twoPi*r*da;
            const auto p=polynomial(x,y);
            samples.block(first,0,1,6)+=dv/volume*p.transpose();samples.block(first+1,6,1,6)+=dv/volume*p.transpose();
            samples.block(first+2,0,1,6)-=dv/volume*y*p.transpose();samples.block(first+2,6,1,6)+=dv/volume*x*p.transpose();
        }
        samples.middleRows(first+3,3)=HydroCompatibleMomentSpace::rzExtraMoments(center,h,[&](Vector2D point)->Eigen::MatrixXd{
            const auto p=polynomial((point.x-center.x)/h.x,(point.y-center.y)/h.y);
            Eigen::MatrixXd v=Eigen::MatrixXd::Zero(2,12);v.block(0,0,1,6)=p.transpose();v.block(1,6,1,6)=p.transpose();return v;});
        metric[first]=metric[first+1]=volume/(diameter*diameter);metric[first+2]=6*metric[first];
        metric.segment(first+3,3).setConstant(12*metric[first]);
        // One energy projection onto quadratic meridional velocities. On-axis
        // ur=r*(a+b*X+c*Y), so its hoop energy is finite for EVERY state.
        // Off-axis use all P2 velocities. Preserve the SAME two weighted means
        // and rotation moment, then fit strain on their nullspace. Constraining
        // these moments avoids dividing by the small hoop stiffness of nearly
        // planar rigid modes far from the axis.
        Eigen::MatrixXd allowed=Eigen::MatrixXd::Identity(12,12);
        if(axisCell[c]){allowed=Eigen::MatrixXd::Zero(12,9);const double rc=center.x/h.x;
            allowed(0,0)=rc;allowed(1,0)=1;allowed(1,1)=rc;allowed(3,1)=1;allowed(2,2)=rc;allowed(4,2)=1;
            allowed.bottomRightCorner(6,6).setIdentity();}
        const Eigen::MatrixXd target=block.strain*samples*allowed;
        Eigen::MatrixXd weighted=Eigen::MatrixXd::Zero(12,allowed.cols()),rhs=Eigen::MatrixXd::Zero(12,n);
        const Eigen::Matrix3d root=gram.llt().matrixL().toDenseMatrix().transpose();
        for(int component=0;component<4;component++){
            weighted.middleRows(3*component,3)=diameter/std::sqrt(volume)*root*target.middleRows(3*component,3);
            rhs.middleRows(3*component,3)=diameter/std::sqrt(volume)*root*block.strain.middleRows(3*component,3);
        }
        const Eigen::MatrixXd constraints=samples.middleRows(first,3)*allowed;
        Eigen::MatrixXd moments=Eigen::MatrixXd::Zero(3,n);moments.block(0,first,3,3).setIdentity();
        Eigen::JacobiSVD<Eigen::MatrixXd> momentSvd(constraints,Eigen::ComputeThinU|Eigen::ComputeFullV);
        require(momentSvd.singularValues()[2]>1e-12*momentSvd.singularValues()[0],"RZ reconstruction cell moments lost rank");
        const Eigen::MatrixXd lift=constraints.transpose()*(constraints*constraints.transpose()).ldlt().solve(moments);
        const Eigen::MatrixXd nullspace=momentSvd.matrixV().rightCols(allowed.cols()-3);
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(weighted*nullspace,Eigen::ComputeThinU|Eigen::ComputeThinV);
        require(svd.singularValues().tail(1)[0]>1e-12*svd.singularValues()[0],"RZ viscous polynomial projection lost rank");
        block.reconstruction=allowed*(lift+nullspace*svd.solve(rhs-weighted*lift));
        const Eigen::MatrixXd difference=Eigen::MatrixXd::Identity(n,n)-samples*block.reconstruction;
        Eigen::MatrixXd material=Eigen::MatrixXd::Zero(12,12);for(int component=0;component<4;component++)material.block<3,3>(3*component,3*component)=mu[c];
        block.dissipation=2*(block.strain.transpose()*material*block.strain
            +mu[c](0,0)/volume*difference.transpose()*metric.asDiagonal()*difference);
        for(int i=0;i<n;i++)for(int j=0;j<n;j++)if(block.dissipation(i,j)!=0)entries.emplace_back(block.rows[i],block.rows[j],block.dissipation(i,j));
    }
    result.matrix_.resize(size,size);result.matrix_.setFromTriplets(entries.begin(),entries.end());return result;
}
}
