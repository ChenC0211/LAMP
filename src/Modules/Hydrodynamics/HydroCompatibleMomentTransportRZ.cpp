#include "Modules/Hydrodynamics/HydroCompatibleMomentTransport.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/NativeLeafDomain.h"
#include "Core/Mesh/CartesianMeasure.h"
#include "Core/Geometry/ParametricHermiteCurve.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <limits>
#include <sstream>

namespace Quantum::Hydrodynamics {
namespace {
using Transport=HydroCompatibleMomentTransport;
void require(bool ok,const char* message){if(!ok)throw std::invalid_argument(message);}
struct Legendre {Eigen::VectorXd value,derivative;};
Legendre legendre(double x,int degree){
    Legendre p;p.value.resize(degree+1);p.derivative.resize(degree+1);p.value[0]=1;p.derivative[0]=0;
    if(degree){p.value[1]=2*x;p.derivative[1]=2;}
    for(int n=2;n<=degree;n++){p.value[n]=((2*n-1)*2*x*p.value[n-1]-(n-1)*p.value[n-2])/n;
        p.derivative[n]=((2*n-1)*(2*p.value[n-1]+2*x*p.derivative[n-1])-(n-1)*p.derivative[n-2])/n;}
    for(int n=1;n<=degree;n++){p.value[n]*=std::sqrt(2.*n+1);p.derivative[n]*=std::sqrt(2.*n+1);}return p;
}
int size(const Transport::Piece& p){return p.axisRegular?p.degree*(p.degree+2):(p.degree+2)*(p.degree+2)-1;}
Eigen::RowVectorXd streamBasis(const Transport::Piece& p,Vector2D local){
    const int k=p.degree;const auto x=legendre(local.x,k+1),y=legendre(local.y,k+1);
    const double length=std::sqrt(p.size.x*p.size.y),radial=p.center.x/p.size.x+local.x;
    Eigen::RowVectorXd psi(size(p));int col=0;
    if(p.axisRegular){for(int j=0;j<=k+1;j++)for(int i=0;i<k;i++)psi[col++]=length*radial*radial*x.value[i]*y.value[j];}
    else for(int j=0;j<=k+1;j++)for(int i=0;i<=k+1;i++)if(i||j)psi[col++]=length*x.value[i]*y.value[j];
    return psi;
}
Eigen::MatrixXd fluxBasis(const Transport::Piece& p,Vector2D local){
    const int k=p.degree;const auto x=legendre(local.x,k+1),y=legendre(local.y,k+1);const double r=p.center.x+p.size.x*local.x;
    const double length=std::sqrt(p.size.x*p.size.y);
    Eigen::MatrixXd b=Eigen::MatrixXd::Zero(2,size(p));int col=0;
    // q=(psi_z,-psi_r) is solenoidal by construction, rather than by a
    // nearly cancelling numerical derivative constraint. On axis psi=r^2*a.
    if(p.axisRegular){
        const double radial=r/p.size.x;
        for(int j=0;j<=k+1;j++)for(int i=0;i<k;i++){
            b(0,col)=length/p.size.y*radial*radial*x.value[i]*y.derivative[j];
            b(1,col++)=-length/p.size.x*(2*radial*x.value[i]+radial*radial*x.derivative[i])*y.value[j];}
    }else{
        for(int j=0;j<=k+1;j++)for(int i=0;i<=k+1;i++)if(i||j){
            b(0,col)=length/p.size.y*x.value[i]*y.derivative[j];b(1,col++)=-length/p.size.x*x.derivative[i]*y.value[j];}
    }return b;
}
template<class F>void volumeQuadrature(const Transport::Piece& p,F&& f){
    const auto& g=Quantum::Geometry::Parametric::gaussEight();const int panels=p.axisRegular?1:4;
    for(int panel=0;panel<panels;panel++)for(int j=0;j<8;j++)for(int i=0;i<8;i++){
        const Vector2D local{-.5+(panel+.5+.5*g.x[i])/panels,.5*g.x[j]};
        f(local,Vector2D{p.center.x+p.size.x*local.x,p.center.y+p.size.y*local.y},.25*p.size.x*p.size.y*g.w[i]*g.w[j]/panels);
    }
}
Eigen::Matrix<double,2,3> particularFlux(Vector2D point,Vector2D center,Vector2D h){
    const double X=(point.x-center.x)/h.x,Y=(point.y-center.y)/h.y;
    Eigen::Matrix<double,2,3> q=Eigen::Matrix<double,2,3>::Zero();
    q.row(0)<<h.x*(X+.5),.5*h.x*(X*X-.25),h.x*(X+.5)*Y;return q;
}
Eigen::Matrix<double,2,3> particularVelocity(Vector2D point,Vector2D center,Vector2D h){
    // If the parent touches the axis, cancel r analytically. An arbitrary
    // divergent algebraic vector need not satisfy physical axis parity;
    // accepted incompressible states have exactly zero particular part.
    if(center.x==.5*h.x){
        const double X=(point.x-center.x)/h.x,Y=(point.y-center.y)/h.y;
        Eigen::Matrix<double,2,3> v=Eigen::Matrix<double,2,3>::Zero();
        v.row(0)<<center.x,.5*center.x*(X-.5),center.x*Y;return v;
    }
    return (center.x/point.x)*particularFlux(point,center,h);
}
}
Eigen::MatrixXd HydroCompatibleMomentTransport::physicalBasis(const Piece& p,Vector2D local){
    if(!p.radialWeighted)return rtBasis(local);
    const double r=p.center.x+p.size.x*local.x;
    Eigen::MatrixXd b;
    if(!p.axisRegular)b=p.radialScale/r*fluxBasis(p,local);
    else {
    const int k=p.degree;const auto x=legendre(local.x,k+1),y=legendre(local.y,k+1);
    const double length=std::sqrt(p.size.x*p.size.y);b=Eigen::MatrixXd::Zero(2,size(p));int col=0;
    for(int j=0;j<=k+1;j++)for(int i=0;i<k;i++){
        b(0,col)=p.radialScale*length/(p.size.x*p.size.x*p.size.y)*r*x.value[i]*y.derivative[j];
        b(1,col++)=-p.radialScale*length/(p.size.x*p.size.x)*(2*x.value[i]+r/p.size.x*x.derivative[i])*y.value[j];}}
    return b;
}
HydroCompatibleMomentTransport HydroCompatibleMomentTransport::buildAxisymmetric(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,const HydroCompatibleMomentSpace& space){
    require(space.cellVelocityCount()==6,"RZ transport requires six independent cell moments");
    HydroCompatibleMomentTransport result;result.identity_=domain.meshIdentity();result.axisymmetric_=true;result.velocityCount_=int(space.mass().rows());result.cells_.resize(mesh.cellIndexNum());
    std::vector<std::vector<int>> incident(mesh.cellIndexNum());
    for(const auto f:domain.faces(mesh)){if(f.negativeCell>=0)incident[f.negativeCell].push_back(f.nativeRow);if(f.positiveCell>=0)incident[f.positiveCell].push_back(f.nativeRow);}
    const auto& g=Quantum::Geometry::Parametric::gaussEight();
    for(int c=0;c<mesh.cellIndexNum();c++){
        const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);const auto center=cell.centerCoord(),h=cell.size();const auto& bulk=space.cells()[c];
        // RT2 bubbles have only two even-in-z radial modes, insufficient for
        // the distinct streamfunction moments weighted by 1,1/r,1/r^2.
        // RT3 supplies three radial bubble modes away from axis. Regularity
        // imposes a second radial zero on axis, requiring RT4 there for the
        // same three modes. This is a fixed space definition, not a retry.
        const int side=incident[c].size()>4?2:1,pieces=side*side,nlocal=int(bulk.rows.size()),k=center.x==.5*h.x?4:3;
        const Vector2D hs{h.x/side,h.y/side};auto& saved=result.cells_[c];saved.pieces.resize(pieces);
        std::vector<int> begin(pieces);int columns=0;
        for(int j=0;j<side;j++)for(int i=0;i<side;i++){
            const int index=j*side+i;auto& p=saved.pieces[index];p.center={center.x+h.x*((i+.5)/side-.5),center.y+h.y*((j+.5)/side-.5)};p.size=hs;
            p.degree=k;p.radialWeighted=true;p.axisRegular=p.center.x==.5*hs.x;p.radialScale=center.x;begin[index]=columns;columns+=size(p);
        }
        const int constraints=4*side*(k+1)+2*side*(side-1)*(k+1)+6;
        Eigen::MatrixXd localD(3,nlocal);for(int i=0;i<3;i++)for(int j=0;j<nlocal;j++)localD(i,j)=space.divergence().coeff(3*c+i,bulk.rows[j]);
        Eigen::Vector3d inverseD;inverseD<<1,12,12;
        const Eigen::MatrixXd divergence=inverseD.asDiagonal()*localD/(CartesianMeasure::twoPi*center.x*h.x*h.y);
        Eigen::JacobiSVD<Eigen::MatrixXd> divergenceSvd(localD,Eigen::ComputeFullV);
        require(divergenceSvd.singularValues()[2]>1e-12*divergenceSvd.singularValues()[0],"RZ divergence constraints lost rank");
        const Eigen::MatrixXd nullspace=divergenceSvd.matrixV().rightCols(nlocal-3);
        const Eigen::MatrixXd solenoidal=nullspace*nullspace.transpose();
        Eigen::MatrixXd a=Eigen::MatrixXd::Zero(constraints,columns),b=Eigen::MatrixXd::Zero(constraints,nlocal),gram=Eigen::MatrixXd::Zero(columns,columns);
        Eigen::MatrixXd cross=Eigen::MatrixXd::Zero(columns,3);int row=0;
        for(int axis=0;axis<2;axis++)for(int sign:{-1,1})for(int s=0;s<side;s++){
            const int i=axis==1?s:(sign<0?0:side-1),j=axis==1?(sign<0?0:side-1):s,piece=j*side+i;
            const auto& p=saved.pieces[piece];const double tangentCenter=axis==1?p.center.x:p.center.y,length=axis==1?hs.x:hs.y;int native=-1;
            for(int f:incident[c]){const auto face=domain.face(mesh,f);if(face.horizontal!=(axis==1)||(face.negativeCell==c?1:-1)!=sign)continue;
                const auto& edge=face.horizontal?mesh.hEdges(face.edgeLabel):mesh.vEdges(face.edgeLabel);const auto va=mesh.verts(edge.edgeToVert[0]).coord(),vb=mesh.verts(edge.edgeToVert[1]).coord();
                const double low=axis==1?va.x:va.y,high=axis==1?vb.x:vb.y;
                if(tangentCenter-.5*length>=low-1e-12*length&&tangentCenter+.5*length<=high+1e-12*length){native=f;break;}}
            require(native>=0,"RZ transport boundary piece crosses native face");const auto& info=space.faces()[native];
            int index=-1;if(info.first>=0){const auto found=std::find(bulk.rows.begin(),bulk.rows.end(),info.first);require(found!=bulk.rows.end(),"RZ face ownership differs");index=int(found-bulk.rows.begin());}
            for(int q=0;q<8;q++){
                const double t=.5*g.x[q],w=.5*g.w[q];const Vector2D local=axis==0?Vector2D{.5*sign,t}:Vector2D{t,.5*sign};
                const auto tests=legendre(t,k);const auto flux=fluxBasis(p,local);const double r=p.center.x+hs.x*local.x,scoord=tangentCenter+length*t;
                for(int m=0;m<=k;m++){
                    a.block(row+m,begin[piece],1,size(p))+=w*tests.value[m]*flux.row(axis);
                    if(index>=0){b(row+m,index)+=w*tests.value[m]*r/center.x;b(row+m,index+1)+=w*tests.value[m]*r/center.x*(scoord-info.tangentCenter)/info.tangentScale;}
                    const Vector2D point{p.center.x+hs.x*local.x,p.center.y+hs.y*local.y};
                    b.row(row+m)-=w*tests.value[m]*particularFlux(point,center,h).row(axis)*divergence;
                }
            }row+=k+1;
        }
        for(int axis=0;axis<2;axis++)for(int j=0;j<side;j++)for(int i=0;i+1<side;i++){
            const int firstPiece=axis==0?j*side+i:i*side+j,secondPiece=axis==0?firstPiece+1:firstPiece+side;
            for(int q=0;q<8;q++){const double t=.5*g.x[q];const auto tests=legendre(t,k);
                for(int which=0;which<2;which++){const int piece=which?secondPiece:firstPiece;const auto& p=saved.pieces[piece];const double sign=which?-1.:1.;
                    const Vector2D local=axis==0?Vector2D{.5*sign,t}:Vector2D{t,.5*sign};const auto flux=fluxBasis(p,local);
                    for(int m=0;m<=k;m++)a.block(row+m,begin[piece],1,size(p))+=sign*.5*g.w[q]*tests.value[m]*flux.row(axis);}
            }row+=k+1;
        }
        const int momentRow=row;
        // Integrate the SAME streamfunction with one cell-wide gauge. Normal
        // continuity fixes differences between neighbouring piece gauges.
        // The spanning path runs along the bottom, then up each column.
        std::vector<Eigen::RowVectorXd> corner(pieces,Eigen::RowVectorXd::Zero(columns)),offsetPsi=corner;
        for(int j=0;j<side;j++)for(int i=0;i<side;i++){
            const int index=j*side+i;
            if(j){const int below=index-side;corner[index]=corner[below];corner[index].segment(begin[below],size(saved.pieces[below]))
                +=streamBasis(saved.pieces[below],{-.5,.5})-streamBasis(saved.pieces[below],{-.5,-.5});}
            else if(i){const int left=index-1;corner[index]=corner[left];corner[index].segment(begin[left],size(saved.pieces[left]))
                +=streamBasis(saved.pieces[left],{.5,-.5})-streamBasis(saved.pieces[left],{-.5,-.5});}
            offsetPsi[index]=corner[index];offsetPsi[index].segment(begin[index],size(saved.pieces[index]))-=streamBasis(saved.pieces[index],{-.5,-.5});
        }
        for(int piece=0;piece<pieces;piece++){
            const auto& p=saved.pieces[piece];const int offset=begin[piece],n=size(p);
            volumeQuadrature(p,[&](Vector2D local,Vector2D point,double da){
                const auto physical=physicalBasis(p,local);
                const auto particular=particularVelocity(point,center,h);
                const double X=(point.x-center.x)/h.x,Y=(point.y-center.y)/h.y,flat=da/(h.x*h.y),weighted=flat*point.x/center.x;
                a.block(momentRow,offset,1,n)+=weighted*physical.row(0);a.block(momentRow+1,offset,1,n)+=weighted*physical.row(1);
                a.block(momentRow+2,offset,1,n)+=weighted*(-Y*physical.row(0)+X*physical.row(1));
                b.row(momentRow)-=weighted*particular.row(0)*divergence;
                b.row(momentRow+1)-=weighted*particular.row(1)*divergence;
                b.row(momentRow+2)-=weighted*(-Y*particular.row(0)+X*particular.row(1))*divergence;
                const double s=center.x/point.x;
                const Eigen::Vector3d tests(X*s,X*X*s*s,Y*s*s);
                Eigen::RowVectorXd psi=offsetPsi[piece];psi.segment(offset,n)+=streamBasis(p,local);
                a.middleRows(momentRow+3,3)+=flat/std::sqrt(h.x*h.y)*tests*psi;
                gram.block(offset,offset,n,n)+=weighted*physical.transpose()*physical;
                cross.middleRows(offset,n)+=weighted*physical.transpose()*particular;
            });
        }
        require(row==momentRow,"RZ transport constraint count differs");for(int i=0;i<6;i++)b(row+i,nlocal-6+i)+=1;row+=6;
        Eigen::LLT<Eigen::MatrixXd> metric(gram);require(metric.info()==Eigen::Success,"RZ reconstruction L2 metric not positive");
        const Eigen::MatrixXd inverseRoot=metric.matrixU().solve(Eigen::MatrixXd::Identity(columns,columns));
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(a*inverseRoot,Eigen::ComputeThinU|Eigen::ComputeThinV);
        const double tolerance=4096*std::numeric_limits<double>::epsilon();int rank=0;
        for(int i=0;i<svd.singularValues().size();i++)if(svd.singularValues()[i]>tolerance*svd.singularValues()[0])rank++;
        const Eigen::MatrixXd inverse=inverseRoot*svd.matrixV().leftCols(rank)*svd.singularValues().head(rank).cwiseInverse().asDiagonal()*svd.matrixU().leftCols(rank).transpose();
        // Minimize the energy of the WHOLE physical velocity, not just its
        // solenoidal remainder. The unconstrained term accounts for the
        // cross work with the fixed divergence lift. All owned DOFs remain.
        const Eigen::MatrixXd unconstrained=-metric.solve(cross)*divergence;
        Eigen::MatrixXd fitted=unconstrained+inverse*(b-a*unconstrained);
        // Two linear refinement corrections, in the same minimum-norm range.
        for(int refinement=0;refinement<2;refinement++)fitted+=inverse*(b-a*fitted);
        saved.constraintDefect=(a*fitted-b).norm()/b.norm();
        saved.constraintCondition=svd.singularValues()[0]/svd.singularValues()[rank-1];saved.coefficientAmplification=fitted.norm()/b.norm();
        if(saved.constraintDefect>=8*tolerance){std::ostringstream message;message<<"RZ solenoidal reconstruction cannot preserve moments: cell="<<c<<" r/h="<<center.x/h.x<<" rank="<<rank<<" columns="<<columns<<" defect="<<saved.constraintDefect
            <<" condition="<<svd.singularValues()[0]/svd.singularValues()[rank-1]<<" coefficient_amplification="<<saved.coefficientAmplification;throw std::invalid_argument(message.str());}
        saved.kineticMass=Eigen::MatrixXd::Zero(nlocal,nlocal);
        for(int piece=0;piece<pieces;piece++){
            auto& p=saved.pieces[piece];
            const Eigen::MatrixXd complete=fitted.middleRows(begin[piece],size(p));
            // The motion/transport query is defined on null(D) and retains
            // its exact streamfunction axis parity. The complete divergent
            // algebraic extension below defines mass, not a second velocity.
            p.coefficients=complete*solenoidal;
            volumeQuadrature(p,[&](Vector2D local,Vector2D point,double da){
                const Eigen::MatrixXd physical=physicalBasis(p,local)*complete+particularVelocity(point,center,h)*divergence;
                saved.kineticMass+=CartesianMeasure::twoPi*point.x*da*physical.transpose()*physical;
            });
        }
        saved.kineticMass=.5*(saved.kineticMass+saved.kineticMass.transpose()).eval();
    }return result;
}
}
