#include "Modules/Hydrodynamics/HydroCompatibleMomentSpace.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/NativeLeafDomain.h"
#include "Core/Mesh/CartesianMeasure.h"
#include "Core/Geometry/ParametricHermiteCurve.h"
#include <cmath>
#include <algorithm>
#include <stdexcept>

namespace Quantum::Hydrodynamics {
namespace {
// phi=X-(rc/h)*log(r/rc); r*grad(phi)=(X,0). Its smooth small-cell
// evaluation avoids subtracting two nearly identical O(X) terms far from axis.
double potential(double r,double center,double h){
    const double t=(r-center)/center;
    if(std::abs(t)<.125){double power=t*t,sum=0;for(int k=2;k<24;k++){sum+=(k%2?-1.:1.)*power/k;power*=t;}return center/h*sum;}
    return center/h*(t-std::log1p(t));
}
// Integral_0^y s^k [s-1/2-log(2s)/2] ds, including the true axis limit.
double axisPrimitive(int k,double y){
    if(y==0)return 0;
    return std::pow(y,k+2)/(k+2)-.5*std::pow(y,k+1)/(k+1)
        -.5*std::pow(y,k+1)/(k+1)*(std::log(2*y)-1./(k+1));
}
template<class F>void integrate(double a,double b,F&& f){
    const auto& g=Quantum::Geometry::Parametric::gaussEight();
    for(int panel=0;panel<4;panel++){const double left=a+(b-a)*panel/4,right=a+(b-a)*(panel+1)/4;
        for(int i=0;i<8;i++)f(.5*(left+right)+.5*(right-left)*g.x[i],.5*(right-left)*g.w[i]);}
}
}
Eigen::MatrixXd HydroCompatibleMomentSpace::rzExtraMoments(Vector2D center,Vector2D h,
    const PhysicalFields& velocity){
    if(!velocity)throw std::invalid_argument("invalid RZ streamfunction input");
    const int n=int(velocity(center).cols());const double r0=center.x-.5*h.x,z0=center.y-.5*h.y,length=std::sqrt(h.x*h.y);
    Eigen::MatrixXd moments=Eigen::MatrixXd::Zero(3,n);const auto& g=Quantum::Geometry::Parametric::gaussEight();
    for(int j=0;j<8;j++){
        const double Y=.5*g.x[j],z=center.y+h.y*Y;
        Eigen::RowVectorXd left=Eigen::RowVectorXd::Zero(n);
        if(r0>0)for(int t=0;t<8;t++){
            const double zz=z0+.5*(z-z0)*(1+g.x[t]);const auto v=velocity({r0,zz});
            if(v.rows()!=2||v.cols()!=n||!v.allFinite())throw std::invalid_argument("invalid RZ physical fields");
            left+=.5*(z-z0)*g.w[t]*r0/center.x*v.row(0);
        }
        for(int panel=0;panel<4;panel++)for(int i=0;i<8;i++){
            const double X=-.5+(panel+.5+.5*g.x[i])/4,r=center.x+h.x*X,s=center.x/r;
            Eigen::RowVectorXd psi=left;
            for(int t=0;t<8;t++){
                const double rr=r0+.5*(r-r0)*(1+g.x[t]);const auto v=velocity({rr,z});
                if(v.rows()!=2||v.cols()!=n||!v.allFinite())throw std::invalid_argument("invalid RZ physical fields");
                psi-=.5*(r-r0)*g.w[t]*rr/center.x*v.row(1);
            }
            const Eigen::Vector3d tests(X*s,X*X*s*s,Y*s*s);
            moments+=.25*g.w[i]*g.w[j]/(4*length)*tests*psi;
        }
    }return moments;
}
Eigen::MatrixXd HydroCompatibleMomentSpace::rzExtraMoments(Vector2D center,Vector2D h,
    const PhysicalFields& velocity,const Partition& planes){
    if(planes[0].empty()&&planes[1].empty())return rzExtraMoments(center,h,velocity);
    if(!velocity)throw std::invalid_argument("invalid partitioned RZ streamfunction input");
    const auto rule=[](double a,double b,const std::vector<double>& cuts,int panels){
        std::vector<double> points;for(int k=0;k<=panels;k++)points.push_back(a+(b-a)*k/panels);
        for(double p:cuts)if(p>a&&p<b)points.push_back(p);std::sort(points.begin(),points.end());points.erase(std::unique(points.begin(),points.end()),points.end());
        std::vector<std::pair<double,double>> out;const auto& g=Geometry::Parametric::gaussEight();
        for(std::size_t k=1;k<points.size();k++)for(int q=0;q<8;q++)out.emplace_back(.5*(points[k]+points[k-1])+.5*(points[k]-points[k-1])*g.x[q],.5*(points[k]-points[k-1])*g.w[q]);return out;
    };
    const int n=int(velocity(center).cols());const double r0=center.x-.5*h.x,z0=center.y-.5*h.y,length=std::sqrt(h.x*h.y);
    Eigen::MatrixXd moments=Eigen::MatrixXd::Zero(3,n);
    const auto sample=[&](Vector2D p){Eigen::MatrixXd v=velocity(p);if(v.rows()!=2||v.cols()!=n||!v.allFinite())throw std::invalid_argument("invalid partitioned RZ physical fields");return v;};
    for(const auto z:rule(z0,z0+h.y,planes[1],1)){
        Eigen::RowVectorXd left=Eigen::RowVectorXd::Zero(n);
        if(r0>0)for(const auto zz:rule(z0,z.first,planes[1],1))left+=zz.second*r0/center.x*sample({r0,zz.first}).row(0);
        for(const auto r:rule(r0,r0+h.x,planes[0],4)){
            Eigen::RowVectorXd psi=left;
            for(const auto rr:rule(r0,r.first,planes[0],1))psi-=rr.second*rr.first/center.x*sample({rr.first,z.first}).row(1);
            const double X=(r.first-center.x)/h.x,Y=(z.first-center.y)/h.y,s=center.x/r.first;
            moments+=(r.second*z.second/(h.x*h.y*length))*Eigen::Vector3d(X*s,X*X*s*s,Y*s*s)*psi;
        }
    }return moments;
}
Eigen::MatrixXd HydroCompatibleMomentSpace::rzFlatMoments(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,
    int cellIndex,const std::vector<int>& faces,const Eigen::MatrixXd& d) const{
    const auto& cell=mesh.cells(mesh.cellIndexToLabel[cellIndex]);const auto center=cell.centerCoord(),h=cell.size();
    const double flat=h.x*h.y,twoPi=CartesianMeasure::twoPi;const int first=2*int(faces.size()),n=int(d.cols());
    const bool axis=center.x==.5*h.x;
    Eigen::MatrixXd result=Eigen::MatrixXd::Zero(6,n);Eigen::Vector3d cellPotential=Eigen::Vector3d::Zero();
    if(axis){cellPotential[0]=flat*axisPrimitive(0,1);cellPotential[1]=flat*(axisPrimitive(1,1)-.5*axisPrimitive(0,1));}
    else integrate(center.x-.5*h.x,center.x+.5*h.x,[&](double r,double dr){const double phi=potential(r,center.x,h.x);
        cellPotential[0]+=h.y*dr*phi;cellPotential[1]+=h.y*dr*phi*(r-center.x)/h.x;});
    Eigen::Vector3d inverseGram;inverseGram<<1/flat,12/flat,12/flat;
    result.row(1)=-cellPotential.transpose()*inverseGram.asDiagonal()*d/(twoPi*flat);
    for(int i=0;i<int(faces.size());i++){
        const auto f=domain.face(mesh,faces[i]);const auto& info=faces_[f.nativeRow];const double sign=f.negativeCell==cellIndex?1.:-1.;
        if(!f.horizontal){result(1,2*i)+=sign*f.area*potential(f.center.x,center.x,h.x)/(twoPi*flat);continue;}
        const auto& edge=mesh.hEdges(f.edgeLabel);const double a=mesh.verts(edge.edgeToVert[0]).coord().x,b=mesh.verts(edge.edgeToVert[1]).coord().x;
        if(axis){const double j1=axisPrimitive(1,b/h.x)-axisPrimitive(1,a/h.x),j2=axisPrimitive(2,b/h.x)-axisPrimitive(2,a/h.x);
            result(1,2*i)+=sign*h.x*h.x*j1/flat;
            result(1,2*i+1)+=sign*h.x*h.x*(h.x*j2-info.tangentCenter*j1)/(flat*info.tangentScale);
        }else integrate(a,b,[&](double r,double dr){const double w=sign*r*dr*potential(r,center.x,h.x)/flat;
            result(1,2*i)+=w;result(1,2*i+1)+=w*(r-info.tangentCenter)/info.tangentScale;});
    }
    // The same retained divergence defines a polynomial particular flux.
    // Integrate d from the CELL LEFT boundary, not the global axis. This
    // bounded lift has no r/h amplification and is zero on that boundary.
    const Eigen::MatrixXd divergence=inverseGram.asDiagonal()*d/(twoPi*center.x);
    const auto particular=[&](Vector2D p)->Eigen::RowVectorXd{
        const double X=(p.x-center.x)/h.x,Y=(p.y-center.y)/h.y;
        return h.x*((X+.5)*divergence.row(0)+.5*(X*X-.25)*divergence.row(1)+(X+.5)*Y*divergence.row(2));};
    struct Segment {int index,side;Vector2D a,b;double length;};std::vector<Segment> boundary;
    for(int i=0;i<int(faces.size());i++){
        const auto f=domain.face(mesh,faces[i]);const int sign=f.negativeCell==cellIndex?1:-1;
        const auto& edge=f.horizontal?mesh.hEdges(f.edgeLabel):mesh.vEdges(f.edgeLabel);
        auto a=mesh.verts(edge.edgeToVert[0]).coord(),b=mesh.verts(edge.edgeToVert[1]).coord();
        const int side=f.horizontal?(sign<0?0:2):(sign>0?1:3);
        const double length=f.horizontal?b.x-a.x:b.y-a.y;if(side>=2)std::swap(a,b);
        boundary.push_back({i,side,a,b,length});
    }
    std::sort(boundary.begin(),boundary.end(),[](const Segment& a,const Segment& b){
        if(a.side!=b.side)return a.side<b.side;
        return a.side==0?a.a.x<b.a.x:a.side==1?a.a.y<b.a.y:a.side==2?a.a.x>b.a.x:a.a.y>b.a.y;});
    Eigen::RowVectorXd start=Eigen::RowVectorXd::Zero(n),rotationBoundary=start,yurBoundary=start,xuzBoundary=start,yuzBoundary=start;
    const auto& g=Quantum::Geometry::Parametric::gaussEight();
    for(const auto& segment:boundary){
        const auto f=domain.face(mesh,faces[segment.index]);const auto& info=faces_[f.nativeRow];const double sign=f.negativeCell==cellIndex?1.:-1.;
        const auto point=[&](double t){return Vector2D{segment.a.x+t*(segment.b.x-segment.a.x),segment.a.y+t*(segment.b.y-segment.a.y)};};
        const auto increment=[&](double end)->Eigen::RowVectorXd{
            Eigen::RowVectorXd value=Eigen::RowVectorXd::Zero(n);
            for(int k=0;k<8;k++){
                const auto p=point(.5*end*(1+g.x[k]));const double dw=sign*segment.length*.5*end*g.w[k];
                value[2*segment.index]+=dw*p.x/center.x;
                value[2*segment.index+1]+=dw*p.x/center.x*((f.horizontal?p.x:p.y)-info.tangentCenter)/info.tangentScale;
                if(!f.horizontal)value-=dw*particular(p);
            }return value;};
        integrate(0,1,[&](double t,double dt){const auto p=point(t);const Eigen::RowVectorXd psi=start+increment(t);
            const double X=(p.x-center.x)/h.x,Y=(p.y-center.y)/h.y,dw=sign*segment.length*dt/flat,s=center.x/p.x;
            rotationBoundary+=dw*(f.horizontal?Y:X)*psi;
            if(f.horizontal)yurBoundary+=dw*s*Y*psi;
            else {xuzBoundary-=dw*s*X*psi;yuzBoundary-=dw*s*Y*psi;}
        });start+=increment(1);
    }
    Eigen::RowVectorXd rotationPart=Eigen::RowVectorXd::Zero(n),yurPart=rotationPart;
    for(int j=0;j<8;j++)for(int i=0;i<8;i++){
        const Vector2D p{center.x+.5*h.x*g.x[i],center.y+.5*h.y*g.x[j]};const double Y=.5*g.x[j],dw=.25*g.w[i]*g.w[j];
        const auto qp=particular(p);rotationPart-=dw*Y*qp;yurPart+=dw*Y*center.x/p.x*qp;
    }
    const double length=std::sqrt(flat),R=center.x/h.x;
    Eigen::RowVectorXd i0=rotationBoundary-rotationPart;i0[first+2]+=1;i0/=length*(1/h.x+1/h.y);
    Eigen::RowVectorXd i1=i0,i2=i0;i1[first+3]-=1/R;i2[first+3]-=2/R;i2[first+4]+=1/(R*R);
    result.row(0)=-h.x/center.x*result.row(1);result(0,first)+=1;
    result.row(2)=yurBoundary+yurPart-length/h.y*i1;
    result.row(4)=xuzBoundary+length/h.x*i2;
    result.row(5)=yuzBoundary;result(5,first+5)-=length/center.x;
    result.row(3)=-h.x/center.x*result.row(4);result(3,first+1)+=1;return result;
}
}
