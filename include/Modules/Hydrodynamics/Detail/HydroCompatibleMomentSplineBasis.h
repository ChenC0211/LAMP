#pragma once
#include <array>
#include <cmath>
#include <vector>
#include <stdexcept>
#include <Eigen/Core>
namespace Quantum::Hydrodynamics::MomentSplineBasis {
namespace Bernstein {
using Polynomial=std::vector<long double>;
inline Polynomial elevate(Polynomial p,int order){
    if(int(p.size())>order+1)throw std::invalid_argument("Bernstein degree exceeds its declared space");
    while(int(p.size())<=order){const int n=int(p.size());Polynomial q(n+1);q[0]=p.front();q[n]=p.back();
        for(int k=1;k<n;k++)q[k]=(k*p[k-1]+(n-k)*p[k])/n;p=std::move(q);}return p;
}
inline Polynomial linear(const Polynomial& p,long double a,long double b){
    const int n=int(p.size());Polynomial q(n+1);
    for(int k=0;k<=n;k++){if(k<n)q[k]+=(n-k)*a*p[k]/n;if(k)q[k]+=k*b*p[k-1]/n;}return q;
}
inline long double choose(int n,int k){long double out=1;for(int j=1;j<=k;j++)out=out*(n-j+1)/j;return out;}
}
inline double basis(int k,double s,int derivative){if(k<derivative)return 0;double factor=1;for(int j=0;j<derivative;j++)factor*=k-j;return factor*std::pow(s,k-derivative);}
// The complete cell moment space belongs to the SAME variational problem.
// Keeping all four directions contains the previous admissible space, so
// energy minimization cannot discard its feasible low-energy competitors.
// All C2 edge jets vanish. At the axis its third radial jet vanishes too.
inline std::array<double,6> interiorJet(double x,double y,bool axis){
    const auto factor=[](double t,int p){const double v=1-t,a=std::pow(t,p),b=std::pow(t,p-1),c=std::pow(t,p-2);
        return std::array<double,3>{a*v*v*v,p*b*v*v*v-3*a*v*v,p*(p-1)*c*v*v*v-6*p*b*v*v+6*a*v};};
    const auto a=factor(x,axis?4:3),b=factor(y,3);
    return {a[0]*b[0],a[1]*b[0],a[0]*b[1],a[2]*b[0],a[1]*b[1],a[0]*b[2]};
}
inline std::array<std::array<double,6>,4> interiorJets(double x,double y,double rc,double hx,bool axis){
    const auto b=interiorJet(x,y,axis);const double X=x-.5,Y=y-.5,d=hx/rc,s=1+d*X;
    const double q[4]={s*s,X*s,X*X,Y},qx[4]={2*d*s,s+d*X,2*X,0},qxx[4]={2*d*d,2*d,2,0},qy[4]={0,0,0,1};
    std::array<std::array<double,6>,4> out;
    for(int k=0;k<4;k++)out[k]={b[0]*q[k],b[1]*q[k]+b[0]*qx[k],b[2]*q[k]+b[0]*qy[k],
        b[3]*q[k]+2*b[1]*qx[k]+b[0]*qxx[k],b[4]*q[k]+b[1]*qy[k]+b[2]*qx[k],b[5]*q[k]+2*b[2]*qy[k]};return out;
}
struct GaussTwelve {double x[12],w[12];GaussTwelve(){
    for(int k=0;k<12;k++){double t=std::cos(std::acos(-1.)*(k+.75)/12.5),d=0;
        for(int iteration=0;iteration<30;iteration++){double a=1,b=t;for(int n=2;n<=12;n++){const double c=((2*n-1)*t*b-(n-1)*a)/n;a=b;b=c;}
            d=12*(t*b-a)/(t*t-1);const double next=t-b/d;if(next==t)break;t=next;}
        x[k]=t;w[k]=2/((1-t*t)*d*d);}
}};
// Same positive Cartesian quadrature, contracted first in x. For a tensor
// streamfunction A_a(x) B_b(y), u=(A B_y,-A_x B). The native 81 by 81
// Gram blocks then need nine by nine radial moments, not one dense outer
// product per point. No polynomial fitting, signed weights or mass lumping.
struct CartesianGram {
    using Axes=Eigen::Matrix<double,3,9>;
    using Block=Eigen::Matrix<double,9,9>;
    Axes axial;
    Block rho00=Block::Zero(),rho11=Block::Zero();
    Block mu00=Block::Zero(),mu11=Block::Zero(),mu22=Block::Zero(),mu02=Block::Zero();
    bool viscous=false;
    explicit CartesianGram(const Axes& y):axial(y){}
    void add(const Axes& x,double massWeight,double viscosityWeight){
        rho00.noalias()+=massWeight*x.row(0).transpose()*x.row(0);
        rho11.noalias()+=massWeight*x.row(1).transpose()*x.row(1);
        if(viscosityWeight>0){
            viscous=true;
            mu00.noalias()+=viscosityWeight*x.row(0).transpose()*x.row(0);
            mu11.noalias()+=viscosityWeight*x.row(1).transpose()*x.row(1);
            mu22.noalias()+=viscosityWeight*x.row(2).transpose()*x.row(2);
            mu02.noalias()+=viscosityWeight*x.row(0).transpose()*x.row(2);
        }
    }
    void addTo(Eigen::MatrixXd& mass,Eigen::MatrixXd& viscosity)const{
        const Block yy00=axial.row(0).transpose()*axial.row(0),yy11=axial.row(1).transpose()*axial.row(1);
        const Block yy22=axial.row(2).transpose()*axial.row(2),yy20=axial.row(2).transpose()*axial.row(0);
        for(int a=0;a<9;a++)for(int c=0;c<9;c++){
            mass.block<9,9>(9*a,9*c)+=rho00(a,c)*yy11+rho11(a,c)*yy00;
            if(viscous)viscosity.block<9,9>(9*a,9*c)+=4*mu11(a,c)*yy11+mu00(a,c)*yy22+mu22(a,c)*yy00
                -(mu02(a,c)*yy20+mu02(c,a)*yy20.transpose());
        }
    }
};
}
