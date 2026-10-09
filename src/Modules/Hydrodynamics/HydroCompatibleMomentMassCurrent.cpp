#include "Modules/Hydrodynamics/HydroCompatibleMomentMassCurrent.h"
#include "Modules/Hydrodynamics/HydroCompatibleMomentTransport.h"
#include "Modules/Hydrodynamics/HydroCompatibleMomentComplement.h"
#include "Modules/Interface/Compatibility/InterfaceParametricMotion.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/NativeLeafDomain.h"
#include "Core/Mesh/CartesianMeasure.h"
#include <Eigen/Cholesky>
#include <algorithm>
#include <cmath>
#include <functional>
#include <exception>
#include <limits>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <tuple>

namespace Quantum::Hydrodynamics {
namespace {
using Current=HydroCompatibleMomentMassCurrent;
using Cell=Current::Cell;
namespace P=Geometry::Parametric;
void require(bool ok,const char* why){if(!ok)throw std::invalid_argument(why);}
struct Gauss {
    std::vector<double> x,w;
    explicit Gauss(int n):x(n),w(n){
        for(int i=0;i<n;i++){double z=std::cos(std::acos(-1.)*(i+.75)/(n+.5)),d=0;
            for(int k=0;k<32;k++){double a=1,b=z;for(int m=2;m<=n;m++){double c=((2*m-1)*z*b-(m-1)*a)/m;a=b;b=c;}
                d=n*(z*b-a)/(z*z-1);double next=z-b/d;if(next==z)break;z=next;}
            x[i]=z;w[i]=2/((1-z*z)*d*d);
        }
    }
};
// Rows: value and derivative in physical coordinates. The zero mode is
// exactly one, including the rational chart, so its gradient is exactly zero.
Eigen::MatrixXd axis(double t,double h,int degree,double beta=0){
    Eigen::MatrixXd out=Eigen::MatrixXd::Zero(2,degree+1);out(0,0)=1;
    if(!degree)return out;
    const double z=2*t-1;out(0,1)=z;out(1,1)=2/h;
    for(int k=1;k<degree;k++){out(0,k+1)=((2*k+1)*z*out(0,k)-k*out(0,k-1))/(k+1);
        out(1,k+1)=((2*k+1)*(2/h*out(0,k)+z*out(1,k))-k*out(1,k-1))/(k+1);}
    if(beta){const double d=1+beta*t;for(int k=1;k<=degree;k++){
        out(1,k)=out(1,k)/(d*d)-2*beta*out(0,k)/(h*d*d*d);out(0,k)/=d*d;}}
    return out;
}
Eigen::MatrixXd scalarJets(const Cell& c,Vector2D p){
    const auto x=axis((p.x-c.lower.x)/c.size.x,c.size.x,c.radialDegree,c.rational?c.size.x/c.lower.x:0);
    const auto y=axis((p.y-c.lower.y)/c.size.y,c.size.y,c.axialDegree);
    Eigen::MatrixXd out(3,(c.radialDegree+1)*(c.axialDegree+1));
    for(int i=0;i<=c.radialDegree;i++)for(int j=0;j<=c.axialDegree;j++){
        const int k=i*(c.axialDegree+1)+j;out(0,k)=x(0,i)*y(0,j);out(1,k)=x(1,i)*y(0,j);out(2,k)=x(0,i)*y(1,j);
    }return out;
}
using Poly=std::vector<long double>;
Poly elevate(Poly p,int n){
    while(int(p.size())<=n){const int d=int(p.size());Poly q(d+1);q[0]=p[0];q[d]=p[d-1];
        for(int k=1;k<d;k++)q[k]=(static_cast<long double>(k)*p[k-1]+(d-k)*p[k])/d;p=std::move(q);}
    require(int(p.size())==n+1,"mass-current face exceeds the shared functional degree");return p;
}
Poly linearProduct(const Poly& p,long double a,long double b){
    const int n=int(p.size());Poly q(n+1);
    for(int k=0;k<=n;k++){if(k<n)q[k]+=(n-k)*a*p[k]/n;if(k)q[k]+=k*b*p[k-1]/n;}return q;
}
std::vector<Poly> legendreTrace(int degree,double a,double b){
    std::vector<Poly> p{{1}};if(!degree)return p;p.push_back({2*a-1,2*(a+b)-1});
    for(int k=1;k<degree;k++){auto q=linearProduct(p[k],2*a-1,2*(a+b)-1);const auto prev=elevate(p[k-1],k+1);
        for(int j=0;j<=k+1;j++)q[j]=((2*k+1)*q[j]-k*prev[j])/(k+1);p.push_back(std::move(q));}
    return p;
}
Eigen::VectorXd faceMoments(const Cell& c,bool horizontal,bool rz,double normal,double lo,double hi,const Eigen::VectorXd& mass){
    constexpr int degree=InterfaceParametricTransport::FaceVolumes::productDegree;
    require(mass.size()==degree+1,"mass current requires the complete shared phase functionals");
    const double a=(lo-(horizontal?c.lower.x:c.lower.y))/(horizontal?c.size.x:c.size.y);
    const double b=(hi-lo)/(horizontal?c.size.x:c.size.y);
    auto polynomials=legendreTrace(horizontal?c.radialDegree:c.axialDegree,a,b);
    const auto fixed=horizontal?axis((normal-c.lower.y)/c.size.y,c.size.y,c.axialDegree)
        :axis((normal-c.lower.x)/c.size.x,c.size.x,c.radialDegree,c.rational?c.size.x/c.lower.x:0);
    Eigen::VectorXd contracted(polynomials.size());
    for(int k=0;k<int(polynomials.size());k++){
        auto p=polynomials[k];
        if(horizontal&&rz&&lo>0){
            if(c.rational&&k>0){const long double scale=c.lower.x/lo;for(auto& v:p)v*=scale*scale;}
            else {p=linearProduct(p,1,hi/lo);p=linearProduct(p,1,hi/lo);}
        }
        // A rational cell never touches r=0; a horizontal subface of an
        // axis cell may begin away from the axis, hence radial degree 18.
        p=elevate(std::move(p),degree);long double sum=0;for(int j=0;j<=degree;j++)sum+=p[j]*mass[j];contracted[k]=double(sum);
    }
    Eigen::VectorXd out((c.radialDegree+1)*(c.axialDegree+1));
    for(int i=0;i<=c.radialDegree;i++)for(int j=0;j<=c.axialDegree;j++)
        out[i*(c.axialDegree+1)+j]=horizontal?contracted[i]*fixed(0,j):fixed(0,i)*contracted[j];
    return out;
}
struct Integral {
    Eigen::VectorXd moments;
    std::vector<Current::MaterialPoint> nodes;
    double estimate=0;
    explicit Integral(int n):moments(Eigen::VectorXd::Zero(n)){}
    void add(Integral other,double scale=1){
        moments+=scale*other.moments;estimate+=std::abs(scale)*other.estimate;
        for(auto& p:other.nodes){p.massWeight*=scale;nodes.push_back(p);}
    }
};
// Positive cut-phase quadrature of the SAME immutable Hermite boundary.
// Monitor the complete scalar product space and (when requested) reference
// current moments. Retained nodes can subsequently assemble velocity work.
Integral integrate(const Cell& cell,const P::Boundary& boundary,bool rz,
    HydroCompatibleMomentMaterial::Properties properties,double tolerance,
    const HydroCompatibleMomentSpline::NativeVelocity* wind=nullptr,double dt=0){
    const int n=(cell.radialDegree+1)*(cell.axialDegree+1),count=wind?2*n:n;
    const auto lo=cell.lower,h=cell.size,hi=lo+h;
    const double volume=h.x*h.y*(rz?CartesianMeasure::twoPi*(lo.x+.5*h.x):1.);
    const double scale=volume*std::max(properties.liquidDensity,properties.gasDensity);
    const double currentScale=wind?std::max(1.,dt*400/std::min(h.x,h.y)):1.;
    const auto difference=[&](const Integral& a,const Integral& b){
        double error=(a.moments.head(n)-b.moments.head(n)).cwiseAbs().maxCoeff()/scale;
        if(wind)error=std::max(error,(a.moments.tail(n)-b.moments.tail(n)).cwiseAbs().maxCoeff()/(scale*currentScale));return error;
    };
    std::vector<P::PhaseMeasure> measures;for(const auto& part:boundary.components())measures.emplace_back(*part.curve);
    const auto& xs=cell.partition[0];
    bool cut=false;int full=0;std::vector<double> cuts=cell.partition[1];
    for(auto& m:measures){
        const auto phase=m.rectanglePhaseMoments(lo.x,hi.x,lo.y,hi.y,rz);cut=cut||phase.phase==P::RectanglePhase::Cut;full+=phase.phase==P::RectanglePhase::Full;
        if(phase.phase!=P::RectanglePhase::Cut)continue;
        const auto add=[&](const auto& arc){const auto y=arc.coordinate(1);std::array<double,5> derivative{};
            cuts.push_back(arc.sampleFraction(0).y);cuts.push_back(arc.sampleFraction(1).y);
            for(int k=1;k<5;k++)derivative[k-1]=k*y[k];for(double t:P::roots(derivative,0,1))cuts.push_back(arc.sampleFraction(t).y);
            for(double x:xs){auto equation=arc.coordinate(0);equation[0]-=x;for(double t:P::roots(equation,0,1))cuts.push_back(arc.sampleFraction(t).y);}
        };
        for(std::size_t k=0;k<m.curve.arcs.size();k++){const auto& b=m.bounds[k];if(b.x[1]<lo.x||b.x[0]>hi.x||b.y[1]<lo.y||b.y[0]>hi.y)continue;add(m.curve.arcs[k]);cuts.push_back(m.sharedEnds[k].y);}
        for(const auto& arc:m.closureEdges)add(arc);
    }
    require(full<=1&&!(full&&cut),"mass-current material components overlap");
    static const Gauss rule(16);
    const auto radialRule=[&](double a,double b,double y,bool liquid){Integral out(count);
        const double rho=liquid?properties.liquidDensity:properties.gasDensity;
        for(int k=0;k<16;k++){Vector2D p{a+.5*(b-a)*(1+rule.x[k]),y};
            const double weight=.5*(b-a)*rule.w[k]*rho*(rz?CartesianMeasure::twoPi*p.x:1.);
            const auto v=scalarJets(cell,p);out.moments.head(n).noalias()+=weight*v.row(0).transpose();
            if(wind){const auto u=wind->firstJet(p).value;out.moments.tail(n).noalias()+=weight*dt*(u.x*v.row(1)+u.y*v.row(2)).transpose();}
            out.nodes.push_back({p,weight,rho});
        }return out;
    };
    std::function<Integral(double,double,double,bool,double,int)> radial;
    radial=[&](double a,double b,double y,bool liquid,double budget,int depth){
        // On each XY phase/support panel, scalar moments have degree 16
        // and rho*u.grad(phi) has degree at most 24 in either coordinate.
        // Gauss-16 is exact through degree 31. RZ quotients still adapt.
        if(!rz)return radialRule(a,b,y,liquid);
        auto coarse=radialRule(a,b,y,liquid);double mid=.5*(a+b);auto fine=radialRule(a,mid,y,liquid);fine.add(radialRule(mid,b,y,liquid));
        const double error=difference(coarse,fine);if(error<=budget){fine.estimate+=error;return fine;}
        require(depth<14&&mid>a&&mid<b,"mass-current radial quadrature unresolved");
        auto out=radial(a,mid,y,liquid,.5*budget,depth+1);out.add(radial(mid,b,y,liquid,.5*budget,depth+1));return out;
    };
    const auto line=[&](double y){Integral out(count);
        if(!cut){for(std::size_t k=1;k<xs.size();k++)out.add(radial(xs[k-1],xs[k],y,full==1,.25*tolerance/(h.y*(xs.size()-1)),0));return out;}
        std::vector<std::array<double,2>> intervals;
        for(auto& m:measures)for(const auto& ab:m.intervals(y)){double a=std::max(lo.x,ab[0]),b=std::min(hi.x,ab[1]);if(b>a)intervals.push_back({a,b});}
        std::sort(intervals.begin(),intervals.end());double begin=lo.x;const double budget=.25*tolerance/(h.y*(2*intervals.size()+1)*(xs.size()-1));
        const auto part=[&](double a,double b,bool liquid){for(std::size_t k=1;k<xs.size();k++){double first=std::max(a,xs[k-1]),last=std::min(b,xs[k]);if(last>first)out.add(radial(first,last,y,liquid,budget,0));}};
        for(const auto& ab:intervals){require(ab[0]>=begin,"mass-current phase slices overlap");part(begin,ab[0],false);part(ab[0],ab[1],true);begin=ab[1];}
        part(begin,hi.x,false);return out;
    };
    cuts.erase(std::remove_if(cuts.begin(),cuts.end(),[&](double y){return y<lo.y||y>hi.y;}),cuts.end());
    std::sort(cuts.begin(),cuts.end());cuts.erase(std::unique(cuts.begin(),cuts.end()),cuts.end());
    Integral total(count);const double pi=std::acos(-1.);
    for(std::size_t k=1;k<cuts.size();k++){double a=cuts[k-1],b=cuts[k];if(!(b>a))continue;
        const auto axialRule=[&](double first,double last){Integral out(count);
            for(int q=0;q<16;q++){double t=first+.5*(last-first)*(1+rule.x[q]);
                double y=cut?a+(b-a)*std::pow(std::sin(.5*pi*t),2):a+(b-a)*t;
                double jacobian=cut?.5*pi*(b-a)*std::sin(pi*t):b-a;
                out.add(line(y),.5*(last-first)*rule.w[q]*jacobian);
            }return out;
        };
        std::function<Integral(double,double,double,int)> axial;
        axial=[&](double first,double last,double budget,int depth){
            auto coarse=axialRule(first,last);double mid=.5*(first+last);auto fine=axialRule(first,mid);fine.add(axialRule(mid,last));
            const double error=difference(coarse,fine);if(error<=budget){fine.estimate+=error;return fine;}
            if(!(depth<14&&mid>first&&mid<last)){std::ostringstream message;message.precision(17);message<<"mass-current axial quadrature unresolved: cell="<<cell.index<<" error="<<error<<" budget="<<budget;throw std::invalid_argument(message.str());}
            auto out=axial(first,mid,.5*budget,depth+1);out.add(axial(mid,last,.5*budget,depth+1));return out;
        };
        // Only an uncut XY panel has polynomial outer dependence; cut
        // endpoints and the sine chart retain the geometric error check.
        total.add(!rz&&!cut?axialRule(0,1):axial(0,1,.75*tolerance/(cuts.size()-1),0));
    }return total;
}
// Geometry-only H1 Riesz map. Constant is the sole null mode. This is a
// definition of the mass-current moments, not a momentum-matrix repair.
struct RieszMap {
    Eigen::MatrixXd gram,gradient;
    Eigen::VectorXd diagonal;
    Eigen::LLT<Eigen::MatrixXd> solver;
    RieszMap(const Cell& cell,bool rz){
        const int nx=cell.radialDegree+1,ny=cell.axialDegree+1,n=nx*ny;
        Eigen::MatrixXd mx=Eigen::MatrixXd::Zero(nx,nx),kx=mx,my=Eigen::MatrixXd::Zero(ny,ny),ky=my;
        static const Gauss rule(24);
        for(int panel=0;panel<2;panel++)for(int k=0;k<24;k++){
            const double t=.5*(panel+.5+.5*rule.x[k]),r=cell.lower.x+cell.size.x*t;
            const double w=.25*cell.size.x*rule.w[k]*(rz?CartesianMeasure::twoPi*r:1.);
            const auto v=axis(t,cell.size.x,cell.radialDegree,cell.rational?cell.size.x/cell.lower.x:0);
            mx.noalias()+=w*v.row(0).transpose()*v.row(0);kx.noalias()+=w*v.row(1).transpose()*v.row(1);
        }
        for(int k=0;k<24;k++){const auto v=axis(.5*(1+rule.x[k]),cell.size.y,cell.axialDegree);const double w=.5*cell.size.y*rule.w[k];
            my.noalias()+=w*v.row(0).transpose()*v.row(0);ky.noalias()+=w*v.row(1).transpose()*v.row(1);}
        gram.resize(n,n);
        for(int i=0;i<nx;i++)for(int j=0;j<nx;j++)gram.block(i*ny,j*ny,ny,ny)=kx(i,j)*my+mx(i,j)*ky;
        gradient=gram.bottomRightCorner(n-1,n-1);
        diagonal=gradient.diagonal().cwiseSqrt().cwiseInverse();
        const Eigen::MatrixXd scaled=diagonal.asDiagonal()*gradient*diagonal.asDiagonal();
        solver.compute(scaled);require(solver.info()==Eigen::Success,"mass-current gradient Riesz map is not positive");
    }
};
void closeCurrent(Cell& cell,const RieszMap& riesz,bool rz,const HydroCompatibleMomentSpline::NativeVelocity& wind,double dt,double scale,double tolerance){
    const int n=int(riesz.gram.rows());
    const auto& gram=riesz.gram;const auto& g=riesz.gradient;const auto& diagonal=riesz.diagonal;const auto& solver=riesz.solver;
    const Eigen::VectorXd rhs=cell.afterMoments-cell.beforeMoments+cell.boundaryMoments-cell.referenceMoments;
    cell.constantDefect=std::abs(rhs[0])/scale;cell.referenceDefect=rhs.cwiseAbs().maxCoeff()/scale;
    // A constant has zero gradient. Do not hide a geometric inventory
    // mismatch by projecting it onto another mode or modifying face mass.
    require(cell.constantDefect<=8*tolerance,"mass-current constant compatibility failed");
    cell.potential=Eigen::VectorXd::Zero(n);
    cell.potential.tail(n-1)=diagonal.asDiagonal()*solver.solve(diagonal.asDiagonal()*rhs.tail(n-1));
    cell.potential.tail(n-1)+=diagonal.asDiagonal()*solver.solve(diagonal.asDiagonal()*(rhs.tail(n-1)-g*cell.potential.tail(n-1)));
    cell.gradientDefect=(gram*cell.potential-rhs).tail(n-1).cwiseAbs().maxCoeff()/scale;
    require(cell.potential.allFinite()&&cell.gradientDefect<512*std::numeric_limits<double>::epsilon(),"mass-current weak continuity solve failed");
    cell.correctionNorm=std::sqrt(std::max(0.,cell.potential.dot(gram*cell.potential)));
    for(const auto& p:cell.midpoint){const auto u=wind.firstJet(p.point).value;
        cell.current.push_back({p.point,{dt*p.massWeight*u.x,dt*p.massWeight*u.y}});
        cell.referenceNorm+=dt*dt*p.density*p.massWeight*(u.x*u.x+u.y*u.y);
    }
    cell.referenceNorm=std::sqrt(cell.referenceNorm);
    const auto& xs=cell.partition[0];const auto& ys=cell.partition[1];
    static const Gauss rule(24);
    for(std::size_t x=1;x<xs.size();x++)for(std::size_t y=1;y<ys.size();y++)
    for(int panel=0;panel<2;panel++)for(int i=0;i<24;i++)for(int j=0;j<24;j++){
        const double hx=xs[x]-xs[x-1],hy=ys[y]-ys[y-1];
        Vector2D p{xs[x-1]+.5*hx*(panel+.5+.5*rule.x[i]),ys[y-1]+.5*hy*(1+rule.x[j])};
        const double w=.125*hx*hy*rule.w[i]*rule.w[j]*(rz?CartesianMeasure::twoPi*p.x:1.);
        const auto v=scalarJets(cell,p);Current::CurrentPoint q{p,{w*v.row(1).dot(cell.potential),w*v.row(2).dot(cell.potential)}};
        cell.correction.push_back(q);cell.current.push_back(q);
    }
}
}

HydroCompatibleMomentMassCurrent HydroCompatibleMomentMassCurrent::build(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,
    const HydroCompatibleMomentSpace& space,const HydroCompatibleMomentSpline::NativeVelocity& wind,
    const InterfaceState& previous,const InterfaceState& next,const InterfaceStepTime& time,
    const HydroCompatibleMomentMaterial& midpoint,double tolerance,const HydroCompatibleMomentComplement* complement){
    midpoint.validate(mesh,domain,space);wind.validate(mesh,domain,space);
    const auto properties=midpoint.properties();
    auto mass=HydroCompatibleMomentTransport::phaseMassFlux(mesh,domain,space,previous,next,time,properties.liquidDensity,properties.gasDensity);
    std::vector<FacePhase> phases(domain.nativeFaceCount());
    for(const auto f:domain.faces(mesh)){const auto& v=InterfaceParametricMotion::faceVolumes(next,f.nativeRow);
        phases[f.nativeRow]={v.total,v.totalPositive,v.totalNegative,v.liquidPositive,v.liquidNegative};}
    return assemble(mesh,domain,space,wind,*previous.parametricGeometry().boundaryHandle(),*next.parametricGeometry().boundaryHandle(),
        time,midpoint,tolerance,complement,std::move(mass),std::move(phases),true);
}
HydroCompatibleMomentMassCurrent HydroCompatibleMomentMassCurrent::buildConstantDensity(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,
    const HydroCompatibleMomentSpace& space,const HydroCompatibleMomentSpline::NativeVelocity& wind,const InterfaceStepTime& time,
    const HydroCompatibleMomentMaterial& midpoint,double tolerance,const HydroCompatibleMomentComplement* complement){
    midpoint.validate(mesh,domain,space);wind.validate(mesh,domain,space);
    require(bool(midpoint.uniformDensity()),"fixed-material convection needs uniform density");
    const auto trace=wind.normalFlux();require(trace.retainedMomentDegree()==2&&trace.retainedVelocityProducts(),"single-phase current needs complete face products");
    const auto volumes=InterfaceParametricMotion::backgroundFaceVolumes(mesh,domain,trace,time.dt);
    auto mass=HydroCompatibleMomentTransport::backgroundMassFlux(mesh,domain,space,volumes,*midpoint.uniformDensity());
    std::vector<FacePhase> phases(volumes.size());for(std::size_t i=0;i<volumes.size();i++){
        const auto& v=volumes[i];phases[i]={v.total,v.totalPositive,v.totalNegative,0,0};}
    return assemble(mesh,domain,space,wind,*midpoint.boundaryHandle(),*midpoint.boundaryHandle(),time,midpoint,tolerance,complement,
        std::move(mass),std::move(phases),false);
}
HydroCompatibleMomentMassCurrent HydroCompatibleMomentMassCurrent::assemble(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,
    const HydroCompatibleMomentSpace& space,const HydroCompatibleMomentSpline::NativeVelocity& wind,
    const Geometry::Parametric::Boundary& previous,const Geometry::Parametric::Boundary& next,const InterfaceStepTime& time,
    const HydroCompatibleMomentMaterial& midpoint,double tolerance,const HydroCompatibleMomentComplement* complement,
    HydroCompatibleMomentTransport::PhaseMassFlux mass,std::vector<FacePhase> phases,bool tracksInterface){
    midpoint.validate(mesh,domain,space);wind.validate(mesh,domain,space);
    if(complement)complement->validate(mesh,domain,space);
    require(std::isfinite(tolerance)&&tolerance>=64*std::numeric_limits<double>::epsilon()&&tolerance<1e-3,"invalid mass-current quadrature tolerance");
    require(std::isfinite(time.begin)&&std::isfinite(time.end)&&std::isfinite(time.dt)&&time.dt>0&&time.end>time.begin
        &&std::abs(time.end-time.begin-time.dt)<=32*std::numeric_limits<double>::epsilon()*std::max({1.,std::abs(time.begin),std::abs(time.end)}),
        "invalid mass-current time interval");
    require(std::abs(wind.time()-(time.begin+.5*time.dt))<=16*std::numeric_limits<double>::epsilon()*std::max(1.,std::abs(wind.time()))
        &&wind.time()==midpoint.geometryTime(),"mass-current midpoint time differs");
    const auto properties=midpoint.properties();
    require(mass.products.size()==std::size_t(domain.nativeFaceCount()),"mass current requires actual high-order shared mass flux");
    HydroCompatibleMomentMassCurrent result;result.inventoryClosure_=mass.maximumRelativeClosure;result.cells_.resize(mesh.cellIndexNum());
    result.identity_=domain.meshIdentity();result.time_=wind.time();result.tracksInterface_=tracksInterface;result.facePhase_=std::move(phases);
    for(int c=0;c<mesh.cellIndexNum();c++){const auto& native=mesh.cells(mesh.cellIndexToLabel[c]);auto& cell=result.cells_[c];
        cell.index=c;cell.size=native.size();cell.lower=native.centerCoord()-cell.size*.5;
        cell.rational=domain.axisymmetric()&&cell.lower.x>0;
        cell.radialDegree=domain.axisymmetric()?(cell.rational?20:18):16;
        cell.partition[0]={cell.lower.x,cell.lower.x+cell.size.x};cell.partition[1]={cell.lower.y,cell.lower.y+cell.size.y};
        if(complement)cell.partition=complement->cells().at(c).partition;
        const int n=(cell.radialDegree+1)*(cell.axialDegree+1);cell.boundaryMoments=Eigen::VectorXd::Zero(n);
    }
    for(const auto f:domain.faces(mesh)){const double lo=(f.horizontal?f.center.x:f.center.y)-.5*f.tangentSpacing;
        for(int c:{f.negativeCell,f.positiveCell})if(c>=0){auto& cell=result.cells_[c];
            cell.boundaryMoments+=(f.negativeCell==c?1:-1)*faceMoments(cell,f.horizontal,domain.axisymmetric(),
                f.horizontal?f.center.y:f.center.x,lo,lo+f.tangentSpacing,mass.products[f.nativeRow]);}}
    // The Riesz matrix depends only on cell size, metric and scalar chart.
    // XY translations leave it unchanged; RZ also needs the exact radius.
    // Build each geometry once, then share immutable factors with workers.
    // This pool belongs to this assembly, so it cannot outlive a remesh.
    using RieszKey=std::tuple<double,double,double,int,int,bool,bool>;
    std::map<RieszKey,std::unique_ptr<const RieszMap>> rieszMaps;
    std::vector<const RieszMap*> cellRiesz(result.cells_.size());
    for(const auto& cell:result.cells_){
        const RieszKey key{cell.size.x,cell.size.y,domain.axisymmetric()?cell.lower.x:0.,cell.radialDegree,cell.axialDegree,cell.rational,domain.axisymmetric()};
        auto entry=rieszMaps.try_emplace(key);
        if(entry.second)entry.first->second=std::make_unique<const RieszMap>(cell,domain.axisymmetric());
        cellRiesz[cell.index]=entry.first->second.get();
    }
    // Face contributions have already been gathered in native order. Each
    // worker now owns one cell's phase quadrature, Riesz solve and nodes;
    // the boundary and complete wind are immutable. No shared reduction.
    std::vector<std::exception_ptr> failures(result.cells_.size());
#pragma omp parallel for schedule(dynamic,1) if(result.cells_.size()>1)
    for(int index=0;index<int(result.cells_.size());index++){try{
        auto& cell=result.cells_[index];
        const int n=int(cell.boundaryMoments.size());
        auto first=integrate(cell,previous,domain.axisymmetric(),properties,tolerance);
        auto last=integrate(cell,next,domain.axisymmetric(),properties,tolerance);
        auto middle=integrate(cell,*midpoint.boundaryHandle(),domain.axisymmetric(),properties,tolerance,&wind,time.dt);
        cell.beforeMoments=std::move(first.moments);cell.afterMoments=std::move(last.moments);cell.referenceMoments=middle.moments.tail(n);
        cell.before=std::move(first.nodes);cell.after=std::move(last.nodes);cell.midpoint=std::move(middle.nodes);
        cell.quadratureEstimate=std::max({first.estimate,last.estimate,middle.estimate});
        const double volume=domain.cellVolumes()[domain.cellRowByIndex()[cell.index]];
        closeCurrent(cell,*cellRiesz[index],domain.axisymmetric(),wind,time.dt,volume*std::max(properties.liquidDensity,properties.gasDensity),tolerance);
    }catch(...){failures[index]=std::current_exception();}}
    // Preserve serial failure selection independently of worker scheduling.
    for(const auto& failure:failures)if(failure)std::rethrow_exception(failure);
    result.mass_=std::move(mass);return result;
}
void HydroCompatibleMomentMassCurrent::validate(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,const HydroCompatibleMomentSpace& space) const{
    space.validate(mesh,domain);require(identity_&&identity_==domain.meshIdentity(),"volume mass current changed its source layout");
}
}
