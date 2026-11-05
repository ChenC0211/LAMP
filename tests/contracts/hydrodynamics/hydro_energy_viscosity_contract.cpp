#include "Modules/Hydrodynamics/HydroEnergyViscosity.h" // Public header is self-contained.
#include "Modules/Hydrodynamics/HydroMaxwellLoad.h"
#include "Modules/Interface/GeometricVofInitialization.h"
#include "Core/Geometry/Detail/FiniteVolumeClosure.h"
#include <chrono>
#include <exception>
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
#include <omp.h>
#endif
#include "Modules/Hydrodynamics/HydroSinglePhaseSaddle.h"
#include "Modules/Hydrodynamics/HydroMomentumGeometry.h"
#include <cerrno>
#include "Core/IO/DictionaryMesh.h"
#include "Core/Mesh/FaceSegmentGraph.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/CartesianMeasure.h"
#include "Core/Mesh/NativeLeafOverlap.h"
#include "Core/Field/NativeTargetLevelField.h"
#include "Core/Field/VertField.h"
#include "Core/Field/HEdgeField.h"
#include "Core/Field/VEdgeField.h"
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <map> // Task41 independent physical-face correspondence.
#include <type_traits> // Task44 independent native-product oracle.
#include <array>
#include <cfenv>
#include <cstdint>
#include <cstring>
#if defined(_MSC_VER) && defined(_M_X64)
#include <float.h>
#include <xmmintrin.h>
#endif
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <string>

namespace
{
using namespace LAMP::Hydrodynamics;
const double pi=std::acos(-1.0);
int algebraCases=0,polynomialCases=0,spectra=0,scalings=0,mmsCases=0,axisLimits=0,saddles=0,rejected=0;
double peakPolynomial=0,peakSymmetry=0,peakEnergy=0,peakRigid=0,peakScaling=0,peakSaddle=0,peakQuadrature=0,peakCentral=0;
void check(bool ok,const char* message) {if (!ok) throw std::runtime_error(message);}
void close(double error,double scale,double& peak,double tolerance=1e-11)
{
    const double d=std::abs(error)/std::max(1.0,scale);peak=std::max(peak,d);
    check(std::isfinite(d) && d<=tolerance,"independent operator gate failed");
}
template<class F> void rejects(F action)
{
    try {action();} catch (const std::invalid_argument&) {++rejected;return;}
    throw std::runtime_error("invalid viscosity input accepted after guard "+std::to_string(rejected));
}
DictionaryMesh input(int n,int metric,double radius=1)
{
    DictionaryMesh d;d.Nx=n;d.Ny=n;d.minCoord={metric==2 ? .25 : 0,0};
    d.maxCoord={metric==2 ? 1.25 : radius,1};d.maxInitialIterations=0;d.maxCellNum=50000;return d;
}
void refine(QuadtreeMesh& mesh,int depth)
{
    std::vector<DictionaryMesh::RefineBox> boxes;
    for (int level=1;level<=depth;++level) {
        DictionaryMesh::RefineBox b;const double lo=level==1 ? .25 : .375;
        const auto w=mesh.maxCoord-mesh.minCoord;b.min=mesh.minCoord+w*lo;b.max=mesh.minCoord+w*(1-lo);b.level=level;boxes.push_back(b);
    }
    if (!boxes.empty()) mesh.applyRefineBoxes(boxes);
}
struct NativeReference {std::vector<Vector2D> centers;std::vector<int> incidence;};
NativeReference reference(const QuadtreeMesh& mesh,const FaceSegmentGraph& graph,const HydroDofLayout& layout)
{
    NativeReference r;r.centers.resize(layout.velocityCount());r.incidence.resize(layout.velocityCount());
    for (const auto& segment : graph.segments()) {
        const auto& edge=segment.isHorizontal ? mesh.hEdges(segment.edgeLabel) : mesh.vEdges(segment.edgeLabel);
        const auto a=mesh.verts(edge.edgeToVert[0]).coord(),b=mesh.verts(edge.edgeToVert[1]).coord();
        const int native=segment.isHorizontal ? mesh.vEdgeIndexNum()+edge.index() : edge.index();const int u=layout.velocityDof(native);
        r.centers[u]=(a+b)*.5;++r.incidence[u];
    }
    return r;
}
struct AnalyticField {Vector2D velocity;std::array<double,4> strain;};
AnalyticField polynomial(Vector2D p,bool rz,int type)
{
    const double x=p.x,y=p.y;
    if (type==1) return rz ? AnalyticField{{0,2},{0,0,0,0}} : AnalyticField{{2,-3},{0,0,0,0}};
    if (type==2 && !rz) return {{-y,x},{0,0,0,0}};
    if (type==2 && rz) {
        auto f=polynomial(p,false,0);f.strain[3]=f.strain[2];f.strain[2]=f.velocity.x/x;return f;
    }
    if (rz) return {{x*(1+2*y),3+x*x+y+2*y*y},{1+2*y,1+4*y,1+2*y,4*x/std::sqrt(2.0)}};
    return {{1+x+2*y+3*x*x+4*x*y+5*y*y,-1+2*x-y+4*x*x+3*x*y+2*y*y},
        {1+6*x+4*y,-1+3*x+4*y,(4+12*x+13*y)/std::sqrt(2.0),0}};
}
AnalyticField smooth(Vector2D p,bool rz)
{
    const double x=p.x,y=p.y,sx=std::sin(pi*x),cx=std::cos(pi*x),sy=std::sin(pi*y),cy=std::cos(pi*y);
    if (rz) return {{x*sy,cx*sy},{sy,pi*cx*cy,sy,(pi*x*cy-pi*sx*sy)/std::sqrt(2.0)}};
    return {{sx*cy,cx*sy},{pi*cx*cy,pi*cx*cy,-std::sqrt(2.0)*pi*sx*sy,0}};
}
template<class F> HydroNativeValues values(const HydroDofLayout& layout,const NativeReference& ref,F field)
{
    std::vector<double> u(layout.velocityCount());
    for (int f=0;f<layout.velocityCount();++f) {const auto v=field(ref.centers[f]).velocity;u[f]=layout.normalAxis(f)==HydroDofLayout::NormalAxis::First ? v.x : v.y;}
    return HydroNativeValues(layout,u,std::vector<double>(layout.pressureCount()));
}
double dot(const std::vector<double>& a,const std::vector<double>& b)
{
    long double result=0;for (std::size_t i=0;i<a.size();++i) result+=static_cast<long double>(a[i])*b[i];return static_cast<double>(result);
}
void saddle(const QuadtreeMesh& mesh,const HydroPhysicalDivergence& metrics,const HydroEnergyViscosity& k,const NativeReference& ref)
{
    const auto& l=metrics.layout();const bool rz=l.domain().axisymmetric();
    std::vector<double> u(l.velocityCount()),p(l.pressureCount()),prescribed(l.velocityCount());
    for (int f=0;f<l.velocityCount();++f) {
        const auto x=ref.centers[f];const auto v=rz ? Vector2D{x.x*x.y,-x.y*x.y} : Vector2D{x.x*x.x,-2*x.x*x.y};
        u[f]=l.normalAxis(f)==HydroDofLayout::NormalAxis::First ? v.x : v.y;
        if (ref.incidence[f]==1 || metrics.faceMeasures()[f]==0) prescribed[f]=u[f];
    }
    for (int c=0;c<l.pressureCount();++c) p[c]=mesh.cells(mesh.cellIndexToLabel[l.nativeCell(c)]).centerCoord().y-.5;
    HydroSaddleVelocityBlock::Matrix a=k.matrix();for (int f=0;f<l.velocityCount();++f) a.coeffRef(f,f)+=1;a.makeCompressed();
    const auto pressure=HydroPressureLoad::build(metrics.divergence());
    const Eigen::Map<const Eigen::VectorXd> eu(u.data(),u.size()),ep(p.data(),p.size());const Eigen::VectorXd rhs=a*eu+pressure.matrix()*ep;
    const auto solve=HydroSinglePhaseSaddle::buildClosed(pressure,HydroSaddleVelocityBlock::build(l,a));
    const auto result=solve.solve(mesh,HydroSaddleRhs(l,{rhs.data(),rhs.data()+rhs.size()},prescribed));
    for (int f=0;f<l.velocityCount();++f) close(result.values().normalVelocity()[f]-u[f],1,peakSaddle);
    for (int c=0;c<l.pressureCount();++c) close(result.values().pressure()[c]-p[c],1,peakSaddle);
    close(result.report().maxCellContinuityDefect,1,peakSaddle);close(result.report().maxPressureMeanDefect,1,peakSaddle);++saddles;
}
void algebra(int n,int metric,int depth,int selection=0)
{
    QuadtreeMesh mesh(input(n,metric));refine(mesh,depth);FaceSegmentGraph graph;
    graph.build(mesh,[selection](const QuadtreeMesh& m,const Cell& c) {
        const auto p=c.centerCoord();const double x=(p.x-m.minCoord.x)/(m.maxCoord.x-m.minCoord.x);
        if (selection==1) return x<.5;
        if (selection==2) return x<.25 || x>.75;
        if (selection==3) return !(x>.375 && x<.625 && p.y>.375 && p.y<.625);
        return true;
    });
    const auto l=HydroDofLayout::build(mesh,NativeLeafDomain::build(mesh,graph,metric!=0));
    const auto metrics=HydroPhysicalDivergence::build(mesh,HydroIntegralDivergence::build(mesh,l));const auto ref=reference(mesh,graph,l);
    const auto points=HydroStrainOperator::volumeQuadrature(mesh,metrics);const auto e=HydroStrainOperator::build(mesh,metrics,points);
    std::vector<std::array<long double,3>> moments(l.pressureCount(),{0,0,0});
    for (const auto& q : points) {moments[q.cellDof][0]+=q.physicalWeight;moments[q.cellDof][1]+=q.physicalWeight*q.point.x;moments[q.cellDof][2]+=q.physicalWeight*q.point.y;}
    for (int p=0;p<l.pressureCount();++p) {
        const auto& cell=mesh.cells(mesh.cellIndexToLabel[l.nativeCell(p)]);const auto c=cell.centerCoord(),w=cell.size();
        const long double a=c.x-w.x/2,b=c.x+w.x/2,volume=metric==0 ? static_cast<long double>(w.x)*w.y : pi*(b*b-a*a)*w.y;
        const long double radialMoment=metric==0 ? volume*c.x : 2*pi*(b*b*b-a*a*a)*w.y/3;
        close(static_cast<double>(moments[p][0]-volume),static_cast<double>(volume),peakQuadrature,1e-13);
        close(static_cast<double>(moments[p][1]-radialMoment),static_cast<double>(std::abs(radialMoment)),peakQuadrature,1e-13);
        close(static_cast<double>(moments[p][2]-volume*c.y),static_cast<double>(std::abs(volume*c.y)),peakQuadrature,1e-13);
    }
    std::vector<double> mu(points.size());for (std::size_t q=0;q<mu.size();++q) mu[q]=q%5==0 ? 0 : .25+static_cast<double>(q%7)/4;
    const auto k=HydroEnergyViscosity::build(e,mu);
    for (int type=0;type<(metric==1 ? 2 : 3);++type) {
        const auto v=values(l,ref,[&](Vector2D p) {return polynomial(p,metric!=0,type);});const auto strain=e.apply(mesh,v);
        long double energy=0;
        for (std::size_t q=0;q<points.size();++q) {
            const auto exact=polynomial(points[q].point,metric!=0,type).strain;
            for (int c=0;c<e.components();++c) {
                close(strain[q*e.components()+c]-exact[c],std::abs(exact[c]),peakPolynomial);
                energy+=2*static_cast<long double>(mu[q])*points[q].physicalWeight*exact[c]*exact[c];
            }
        }
        close(k.dissipation(mesh,v)-static_cast<double>(energy),std::abs(static_cast<double>(energy)),peakEnergy);
        close(dot(v.normalVelocity(),k.apply(mesh,v))-static_cast<double>(energy),std::abs(static_cast<double>(energy)),peakEnergy);
        if (type==1 || (type==2 && metric==0)) close(k.dissipation(mesh,v),1,peakRigid);
        ++polynomialCases;
    }
    std::vector<double> random(l.velocityCount());for (int f=0;f<l.velocityCount();++f) random[f]=std::sin(3.7*f)+std::cos(.31*f);
    for (int f : metrics.axisVelocityDofs()) random[f]=0;
    const HydroNativeValues v(l,random,std::vector<double>(l.pressureCount()));const auto strain=e.apply(mesh,v);
    if (depth==0 && metric!=1) for (int p=0;p<l.pressureCount();++p) {
        const auto& cell=mesh.cells(mesh.cellIndexToLabel[l.nativeCell(p)]);const auto c=cell.centerCoord(),w=cell.size();
        std::array<int,4> face{-1,-1,-1,-1};
        for (int f=0;f<l.velocityCount();++f) {
            const auto x=ref.centers[f];
            if (l.normalAxis(f)==HydroDofLayout::NormalAxis::First && x.y==c.y) {if (x.x==c.x-w.x/2) face[0]=f;if (x.x==c.x+w.x/2) face[1]=f;}
            if (l.normalAxis(f)==HydroDofLayout::NormalAxis::Second && x.x==c.x) {if (x.y==c.y-w.y/2) face[2]=f;if (x.y==c.y+w.y/2) face[3]=f;}
        }
        check(std::all_of(face.begin(),face.end(),[](int f) {return f>=0;}),"independent regular face stencil");
        double dx=0,dy=0;for (std::size_t q=0;q<points.size();++q) if (points[q].cellDof==p) {dx+=strain[q*e.components()]/4;dy+=strain[q*e.components()+1]/4;}
        close(dx-(random[face[1]]-random[face[0]])/w.x,1,peakCentral);
        close(dy-(random[face[3]]-random[face[2]])/w.y,1,peakCentral);
    }
    long double independent=0;for (std::size_t q=0;q<mu.size();++q) for (int c=0;c<e.components();++c) {
        const double s=strain[q*e.components()+c];independent+=2*static_cast<long double>(mu[q])*points[q].physicalWeight*s*s;
    }
    close(dot(random,k.apply(mesh,v))-static_cast<double>(independent),std::abs(static_cast<double>(independent)),peakEnergy);
    check(k.dissipation(mesh,v)>=0,"negative energy");
    const Eigen::MatrixXd dense=k.matrix();close((dense-dense.transpose()).norm(),dense.norm(),peakSymmetry,1e-12);
    const auto full=HydroEnergyViscosity::build(e,std::vector<double>(mu.size(),1));
    if (n==4 && depth<=1) {
        std::vector<int> free;for (int f=0;f<l.velocityCount();++f) if (std::find(metrics.axisVelocityDofs().begin(),metrics.axisVelocityDofs().end(),f)==metrics.axisVelocityDofs().end()) free.push_back(f);
        const Eigen::MatrixXd d=full.matrix();Eigen::MatrixXd block(free.size(),free.size());
        for (Eigen::Index i=0;i<block.rows();++i) for (Eigen::Index j=0;j<block.cols();++j) block(i,j)=d(free[static_cast<std::size_t>(i)],free[static_cast<std::size_t>(j)]);
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigen(block);check(eigen.info()==Eigen::Success,"spectral verification failed");
        const double threshold=1e-11*std::max(1.0,eigen.eigenvalues().maxCoeff());int zeros=0;
        for (Eigen::Index i=0;i<eigen.eigenvalues().size();++i) {const double lambda=eigen.eigenvalues()[i];check(lambda>=-threshold,"not PSD");if (std::abs(lambda)<=threshold) ++zeros;}
        std::cout<<"SPECTRUM,"<<metric<<','<<depth<<','<<block.rows()<<','<<eigen.eigenvalues().minCoeff()<<','<<zeros<<'\n';
        check(zeros==(metric==0 ? 3 : 1),"extra non-rigid zero energy mode");++spectra;
    }
    const auto baseForce=full.apply(mesh,v);const double baseEnergy=full.dissipation(mesh,v);
    for (double scale : {0.,1e-3,1e-6,1e-12,1e-300}) {
        const auto scaled=HydroEnergyViscosity::build(e,std::vector<double>(mu.size(),scale));const auto force=scaled.apply(mesh,v);
        if (scale==0) {check(scaled.matrix().nonZeros()==0 && scaled.dissipation(mesh,v)==0,"zero viscosity residual");for (double f : force) check(f==0,"zero force residual");}
        else {
            double defect=0,referenceNorm=0;for (std::size_t i=0;i<force.size();++i) {defect=std::max(defect,std::abs(force[i]/scale-baseForce[i]));referenceNorm=std::max(referenceNorm,std::abs(baseForce[i]));}
            close(defect,referenceNorm,peakScaling);close(scaled.dissipation(mesh,v)/scale-baseEnergy,std::abs(baseEnergy),peakScaling);
        }
        ++scalings;
    }
    if (n==4 && depth<=1) saddle(mesh,metrics,full,ref);
    std::cout<<"ALGEBRA,"<<metric<<','<<n<<','<<depth<<','<<selection<<','<<l.velocityCount()<<','<<points.size()<<','<<e.matrix().nonZeros()<<'\n';++algebraCases;
}
// Independent GL8 volume reference; never reads the production quadrature or E.
double referenceEnergy(bool rz)
{
    const std::array<double,4> node{.1834346424956498,.5255324099163290,.7966664774136267,.9602898564975363};
    const std::array<double,4> weight{.3626837833783620,.3137066458778873,.2223810344533745,.1012285362903763};
    long double result=0;
    for (int a=0;a<8;++a) for (int b=0;b<8;++b) {
        const double x=(1+(a<4 ? -1 : 1)*node[a%4])/2,y=(1+(b<4 ? -1 : 1)*node[b%4])/2;
        const auto s=smooth({x,y},rz).strain;long double norm=0;for (int c=0;c<(rz ? 4 : 3);++c) norm+=s[c]*s[c];
        result+=2*norm*weight[a%4]*weight[b%4]/4*(rz ? 2*pi*x : 1);
    }
    return static_cast<double>(result);
}
void mms(bool rz,bool amr)
{
    double previousStrain=0,previousEnergy=0,strainOrder=0,energyOrder=0;
    const double exactEnergy=referenceEnergy(rz);
    if (!rz) check(std::abs(exactEnergy-2*pi*pi)<1e-7,"independent GL8 reference");
    for (int n : {4,8,16}) {
        QuadtreeMesh mesh(input(n,rz ? 1 : 0));if (amr) refine(mesh,1);FaceSegmentGraph graph;graph.build(mesh);
        const auto l=HydroDofLayout::build(mesh,NativeLeafDomain::build(mesh,graph,rz));const auto metrics=HydroPhysicalDivergence::build(mesh,HydroIntegralDivergence::build(mesh,l));
        const auto points=HydroStrainOperator::volumeQuadrature(mesh,metrics);const auto e=HydroStrainOperator::build(mesh,metrics,points);
        const auto v=values(l,reference(mesh,graph,l),[rz](Vector2D p) {return smooth(p,rz);});const auto s=e.apply(mesh,v);
        long double error=0,norm=0;for (std::size_t q=0;q<points.size();++q) {
            const auto exact=smooth(points[q].point,rz).strain;
            for (int c=0;c<e.components();++c) {const double diff=s[q*e.components()+c]-exact[c];error+=points[q].physicalWeight*diff*diff;norm+=points[q].physicalWeight*exact[c]*exact[c];}
        }
        const double strainError=std::sqrt(static_cast<double>(error/norm));
        const auto k=HydroEnergyViscosity::build(e,std::vector<double>(points.size(),1));const double energyError=std::abs(k.dissipation(mesh,v)-exactEnergy)/exactEnergy;
        std::cout<<"MMS,"<<rz<<','<<amr<<','<<n<<','<<strainError<<','<<energyError<<'\n';
        if (n>4) {check(strainError<previousStrain && energyError<previousEnergy,"MMS error did not decrease");strainOrder=std::log2(previousStrain/strainError);energyOrder=std::log2(previousEnergy/energyError);}
        previousStrain=strainError;previousEnergy=energyError;++mmsCases;
    }
    std::cout<<"ORDER,"<<rz<<','<<amr<<','<<strainOrder<<','<<energyOrder<<'\n';check(strainOrder>=1.8 && energyOrder>=1.8,"strain/energy second-order MMS gate");
}
void axis()
{
    for (double radius : {1.,1e-3,1e-6}) {
        QuadtreeMesh mesh(input(4,1,radius));FaceSegmentGraph graph;graph.build(mesh);
        const auto l=HydroDofLayout::build(mesh,NativeLeafDomain::build(mesh,graph,true));const auto metrics=HydroPhysicalDivergence::build(mesh,HydroIntegralDivergence::build(mesh,l));
        const auto ref=reference(mesh,graph,l);const auto v=values(l,ref,[](Vector2D p) {return polynomial(p,true,0);});
        std::vector<HydroStrainPoint> points;
        for (int p=0;p<l.pressureCount();++p) {const auto& c=mesh.cells(mesh.cellIndexToLabel[l.nativeCell(p)]);if (c.centerCoord().x-c.size().x/2==0) points.push_back({p,{0,c.centerCoord().y},0});}
        const auto e=HydroStrainOperator::build(mesh,metrics,points);const auto s=e.apply(mesh,v);
        for (std::size_t q=0;q<points.size();++q) {close(s[4*q+2]-(1+2*points[q].point.y),4,peakPolynomial);close(s[4*q+3],1,peakPolynomial);check(s[4*q]==s[4*q+2],"hoop limit differs from radial strain");}
        ++axisLimits;
    }
}
void guards()
{
    QuadtreeMesh mesh(input(4,1));FaceSegmentGraph graph;graph.build(mesh);
    const auto l=HydroDofLayout::build(mesh,NativeLeafDomain::build(mesh,graph,true));const auto metrics=HydroPhysicalDivergence::build(mesh,HydroIntegralDivergence::build(mesh,l));
    const auto points=HydroStrainOperator::volumeQuadrature(mesh,metrics);const auto e=HydroStrainOperator::build(mesh,metrics,points);const auto k=HydroEnergyViscosity::build(e,std::vector<double>(points.size(),1));
    const HydroNativeValues zero(l,std::vector<double>(l.velocityCount()),std::vector<double>(l.pressureCount()));const auto accepted=k.apply(mesh,zero);
    rejects([&]{HydroStrainOperator::build(mesh,metrics,{});});
    for (int kind=0;kind<5;++kind) {auto q=points;
        if (kind==0) q[0].cellDof=-1;
        if (kind==1) q[0].point.x=-1;
        if (kind==2) q[0].physicalWeight=-1;
        if (kind==3) q[0].point.y=std::numeric_limits<double>::quiet_NaN();
        if (kind==4) q[0].physicalWeight=std::numeric_limits<double>::infinity();
        rejects([&]{HydroStrainOperator::build(mesh,metrics,q);});
    }
    rejects([&]{HydroEnergyViscosity::build(e,{});});
    for (double bad : {-1.,std::numeric_limits<double>::quiet_NaN(),std::numeric_limits<double>::infinity()}) {std::vector<double> mu(points.size(),1);mu[0]=bad;rejects([&]{HydroEnergyViscosity::build(e,mu);});}
    auto overflow=points;overflow[0].physicalWeight=std::numeric_limits<double>::max();const auto largeE=HydroStrainOperator::build(mesh,metrics,overflow);
    rejects([&]{HydroEnergyViscosity::build(largeE,std::vector<double>(points.size(),1));});
    auto axisU=zero.normalVelocity();axisU[metrics.axisVelocityDofs()[0]]=std::numeric_limits<double>::denorm_min();const HydroNativeValues badAxis(l,axisU,zero.pressure());
    rejects([&]{e.apply(mesh,badAxis);});rejects([&]{k.apply(mesh,badAxis);});rejects([&]{k.dissipation(mesh,badAxis);});
    const auto other=HydroDofLayout::build(mesh,NativeLeafDomain::build(mesh,graph,true));const HydroNativeValues wrong(other,zero.normalVelocity(),zero.pressure());
    rejects([&]{e.apply(mesh,wrong);});rejects([&]{k.apply(mesh,wrong);});
    int largestColumn=0;double largest=0;
    for (Eigen::Index row=0;row<e.matrix().outerSize();++row) for (HydroStrainOperator::Matrix::InnerIterator entry(e.matrix(),row);entry;++entry)
        if (std::abs(entry.value())>largest) {largest=std::abs(entry.value());largestColumn=static_cast<int>(entry.col());}
    check(largest>1,"overflow fixture has no amplifying derivative");
    auto huge=zero.normalVelocity();huge[largestColumn]=std::numeric_limits<double>::max();const HydroNativeValues hugeValues(l,huge,zero.pressure());
    rejects([&]{e.apply(mesh,hugeValues);});rejects([&]{k.apply(mesh,hugeValues);});rejects([&]{k.dissipation(mesh,hugeValues);});
    const auto inviscid=HydroEnergyViscosity::build(e,std::vector<double>(points.size()));
    check(inviscid.dissipation(mesh,hugeValues)==0,"zero viscosity evaluated an overflowing strain");
    for (double f : inviscid.apply(mesh,hugeValues)) check(f==0,"inviscid extreme velocity force");
    QuadtreeMesh stretched(input(4,1,2));rejects([&]{e.apply(stretched,zero);});rejects([&]{k.apply(stretched,zero);});rejects([&]{HydroStrainOperator::build(stretched,metrics,points);});
    QuadtreeMesh tiny(input(1,0));FaceSegmentGraph tg;tg.build(tiny);const auto tl=HydroDofLayout::build(tiny,NativeLeafDomain::build(tiny,tg,false));
    const auto tm=HydroPhysicalDivergence::build(tiny,HydroIntegralDivergence::build(tiny,tl));rejects([&]{HydroStrainOperator::build(tiny,tm,HydroStrainOperator::volumeQuadrature(tiny,tm));});
    refine(mesh,1);rejects([&]{k.apply(mesh,zero);});rejects([&]{HydroStrainOperator::volumeQuadrature(mesh,metrics);});
    QuadtreeMesh equivalent(input(4,1));check(k.apply(equivalent,zero)==accepted,"failed operation changed source");
    std::optional<HydroEnergyViscosity> retained;std::optional<HydroNativeValues> rv;
    {QuadtreeMesh temporary(input(4,1));FaceSegmentGraph g;g.build(temporary);const auto dl=HydroDofLayout::build(temporary,NativeLeafDomain::build(temporary,g,true));
     const auto m=HydroPhysicalDivergence::build(temporary,HydroIntegralDivergence::build(temporary,dl));const auto s=HydroStrainOperator::build(temporary,m,HydroStrainOperator::volumeQuadrature(temporary,m));
     retained=HydroEnergyViscosity::build(s,std::vector<double>(s.points().size(),1));rv.emplace(dl,std::vector<double>(dl.velocityCount()),std::vector<double>(dl.pressureCount()));}
    check(retained->apply(equivalent,*rv)==accepted,"borrowed mesh lifetime");
}

// Requested RT0 capability: use the existing native/P2 fixtures and full
// factory as a value/CSR reference, never as a cross-request FP reference.
namespace Rt0GreenContract {
using Matrix=HydroStrainOperator::Matrix;
std::size_t geometries=0,missingCapabilities=0,loadCases=0,ownFpPairs=0;
std::uint64_t bits(double x){std::uint64_t b=0;std::memcpy(&b,&x,sizeof b);return b;}
void scalar(double a,double b,const char* why){check(bits(a)==bits(b),why);}
void matrix(const Matrix& a,const Matrix& b){
    check(a.rows()==b.rows()&&a.cols()==b.cols()&&a.nonZeros()==b.nonZeros()&&a.isCompressed()&&b.isCompressed(),"Green CSR shape differs from original Full");
    for(Eigen::Index i=0;i<=a.outerSize();++i)check(a.outerIndexPtr()[i]==b.outerIndexPtr()[i],"Green CSR row offsets differ from original Full");
    for(Eigen::Index i=0;i<a.nonZeros();++i){check(a.innerIndexPtr()[i]==b.innerIndexPtr()[i],"Green CSR columns differ from original Full");scalar(a.valuePtr()[i],b.valuePtr()[i],"Green CSR coefficient differs from original Full bits");}
}
void vector(const std::vector<double>& a,const std::vector<double>& b){
    check(a.size()==b.size(),"Green native load shape changed");for(std::size_t i=0;i<a.size();++i)scalar(a[i],b[i],"Green native load value/sign differs from original Full");
}
void same(const HydroStrainOperator& a,const HydroStrainOperator& b){
    check(a.isRt0GreenCompatible()&&b.isRt0GreenCompatible()&&a.layout().sameIdentity(b.layout())
        &&&a.metrics().divergence().matrix()==&b.metrics().divergence().matrix(),"Green geometry lost actual layout/D owner");
    matrix(a.stressVolumeVelocityMatrix(),b.stressVolumeVelocityMatrix());matrix(a.stressFaceVelocityMatrix(),b.stressFaceVelocityMatrix());
    const auto& ap=a.points();const auto& bp=b.points();check(ap.size()==bp.size(),"Green volume rule changed");
    for(std::size_t i=0;i<ap.size();++i){check(ap[i].cellDof==bp[i].cellDof,"Green volume owner changed");scalar(ap[i].point.x,bp[i].point.x,"Green volume x changed");scalar(ap[i].point.y,bp[i].point.y,"Green volume y changed");scalar(ap[i].physicalWeight,bp[i].physicalWeight,"Green volume measure changed");}
    const auto& af=a.stressFacePoints();const auto& bf=b.stressFacePoints();check(af.size()==bf.size(),"Green face rule changed");
    for(std::size_t i=0;i<af.size();++i){check(af[i].velocityDof==bf[i].velocityDof&&af[i].nativeFace==bf[i].nativeFace&&af[i].negativeCellDof==bf[i].negativeCellDof&&af[i].positiveCellDof==bf[i].positiveCellDof,"Green face owner/incidence changed");scalar(af[i].point.x,bf[i].point.x,"Green face x changed");scalar(af[i].point.y,bf[i].point.y,"Green face y changed");scalar(af[i].physicalWeight,bf[i].physicalWeight,"Green face measure changed");}
    ++geometries;
}
template<class F>void missing(F f){bool rejected=false;try{f();}catch(const std::invalid_argument&){rejected=true;}check(rejected,"Green-only exposed an unrequested raw/Ehat capability");++missingCapabilities;}
void capabilities(const QuadtreeMesh& mesh,const HydroStrainOperator& g){
    check(g.isRt0GreenCompatible()&&!g.isStressCompatible()&&!g.hasCanonicalRawStrain(),"Green-only capability was relabelled Full/raw");
    const auto& l=g.layout();const HydroNativeValues zero(l,std::vector<double>(l.velocityCount()),std::vector<double>(l.pressureCount()));
    missing([&]{(void)g.matrix();});missing([&]{(void)g.apply(mesh,zero);});missing([&]{(void)g.canonicalRawStrain();});
    missing([&]{(void)g.stressComponentIntegralMatrix();});missing([&]{(void)g.stressCompatibilityReport();});
    missing([&]{(void)HydroEnergyViscosity::build(g,std::vector<double>(g.points().size(),1.));});
}
void reused(const HydroStrainOperator& g){
    const auto& r=g.geometryReuseReport();check(r.requested&&r.enabled&&r.targetCells==static_cast<std::size_t>(g.layout().pressureCount())
        &&r.reusedFaces>0&&r.rebuiltFaces>0&&r.recomputedFitCells>0,"Green remap skipped original clean/changed/support work");
    check(r.reusedRawCells==0&&r.reusedStressCells==0&&r.rawDirectReboundRows==0,
        "Green remap reports omitted raw/Ehat work");
}
void physical(const QuadtreeMesh& mesh,const FaceSegmentGraph& graph,const HydroStrainOperator& g){
    const auto& l=g.layout();const bool rz=l.domain().axisymmetric();const auto ref=reference(mesh,graph,l);
    const auto poly=values(l,ref,[rz](Vector2D p){return polynomial(p,rz,0);});
    const Eigen::Map<const Eigen::VectorXd> u(poly.normalVelocity().data(),l.velocityCount());
    const Eigen::VectorXd volume=g.stressVolumeVelocityMatrix()*u,face=g.stressFaceVelocityMatrix()*u;
    const auto bound=[](long double e,long double scale){check(std::isfinite(e)&&std::isfinite(scale)&&scale>=0&&(scale==0?e==0:std::abs(e)<=1e-12L*scale),"Green independent original P2/trace gate failed");};
    for(std::size_t q=0;q<g.points().size();++q){const auto exact=polynomial(g.points()[q].point,rz,0).velocity;const long double scale=std::abs(exact.x)+std::abs(exact.y);bound(volume[2*q]-exact.x,scale);bound(volume[2*q+1]-exact.y,scale);}
    for(std::size_t q=0;q<g.stressFacePoints().size();++q){const auto& p=g.stressFacePoints()[q];const int normal=l.normalAxis(p.velocityDof)==HydroDofLayout::NormalAxis::First?0:1,tangent=1-normal;const auto exact=polynomial(p.point,rz,0).velocity;
        check(face[2*q+normal]==poly.normalVelocity()[p.velocityDof],"Green normal trace changed actual native mean");bound(face[2*q+tangent]-(tangent==0?exact.x:exact.y),std::abs(exact.x)+std::abs(exact.y));}
}
struct Floating {
    std::fenv_t native{};int error=errno;
#if defined(_MSC_VER) && defined(_M_X64)
    unsigned csr=_mm_getcsr();
#endif
    Floating(){check(std::fegetenv(&native)==0,"Green FP capture failed");errno=error;}
    void restore()const{check(std::fesetenv(&native)==0,"Green FP restore failed");
#if defined(_MSC_VER) && defined(_M_X64)
        _mm_setcsr(csr);
#endif
        errno=error;}
    void same(const Floating& b)const{check(std::memcmp(&native,&b.native,sizeof native)==0&&error==b.error,"Green own-request full FP/errno differs across threads");
#if defined(_MSC_VER) && defined(_M_X64)
        check(csr==b.csr,"Green own-request raw MXCSR differs across threads");
#endif
    }
};
} // namespace Rt0GreenContract

namespace StressCompatibleContract {
using Matrix=HydroStrainOperator::Matrix;
using DenseRows=std::vector<std::vector<long double>>;
void bound(long double error,long double scale,double& peak,const char* why) {
    check(std::isfinite(error)&&std::isfinite(scale)&&scale>=0,why);
    if(scale==0){check(error==0,why);return;}
    const auto relative=std::abs(error)/scale;peak=std::max(peak,static_cast<double>(relative));check(relative<=1e-12L,why);
}
void add(std::vector<long double>& to,const Matrix& from,Eigen::Index row,long double factor) {
    for(Matrix::InnerIterator e(from,row);e;++e)to[static_cast<std::size_t>(e.col())]+=factor*e.value();
}
void scene(bool rz) {
    QuadtreeMesh mesh(input(8,rz?1:0));refine(mesh,1);FaceSegmentGraph graph;graph.build(mesh);
    const auto l=HydroDofLayout::build(mesh,NativeLeafDomain::build(mesh,graph,rz));
    check(!l.domain().uniformLeaves(),"stress geometry fixture is not actual 2:1 AMR");
    const auto metrics=HydroPhysicalDivergence::build(mesh,HydroIntegralDivergence::build(mesh,l));
    const auto points=HydroStrainOperator::volumeQuadrature(mesh,metrics);
    const auto rawStrain=HydroStrainOperator::build(mesh,metrics,points);const auto before=rawStrain.matrix();
    const auto strain=HydroStrainOperator::buildStressCompatible(mesh,metrics);
    check(strain.isStressCompatible()&&!rawStrain.isStressCompatible(),"stress-compatible factory/tag lost");
    check((rawStrain.matrix()-before).norm()==0,"stress geometry changed existing raw strain owner");
    check(strain.isRt0GreenCompatible()&&!rawStrain.isRt0GreenCompatible(),"Full/raw Green capability differs from actual factory");
    const auto green=HydroStrainOperator::buildRt0Green(mesh,metrics);
    Rt0GreenContract::same(green,strain);Rt0GreenContract::capabilities(mesh,green);Rt0GreenContract::physical(mesh,graph,green);
    const auto& facePoints=strain.stressFacePoints();const auto& fv=strain.stressFaceVelocityMatrix();
    const auto& vv=strain.stressVolumeVelocityMatrix();const auto& integral=strain.stressComponentIntegralMatrix();
    const int pc=l.pressureCount(),nu=l.velocityCount(),nc=strain.components();const long double root2=std::sqrt(2.L);
    check(strain.points().size()==points.size()&&facePoints.size()==2*static_cast<std::size_t>(nu)&&
        fv.rows()==static_cast<Eigen::Index>(facePoints.size())*2&&vv.rows()==static_cast<Eigen::Index>(points.size())*2&&
        integral.rows()==pc*nc&&fv.cols()==nu&&vv.cols()==nu&&integral.cols()==nu,"stress geometry native shape mismatch");
    DenseRows target(static_cast<std::size_t>(pc)*nc,std::vector<long double>(nu)),actual=target,rawScale=target;
    DenseRows dx(pc,std::vector<long double>(nu)),dy=dx,hr=dx,hz=dx;
    std::vector<long double> weights(pc);std::vector<int> pointCounts(nu);std::size_t hanging=0,axisPoints=0;
    double peakMoment=0,peakTrace=0,peakConstant=0,peakPressure=0,peakWork=0;
    for(int u=0;u<nu;++u){const auto f=l.domain().face(mesh,l.nativeFace(u));
        const auto& edge=f.horizontal?mesh.hEdges(f.edgeLabel):mesh.vEdges(f.edgeLabel);
        const auto a=mesh.verts(edge.edgeToVert[0]).coord(),b=mesh.verts(edge.edgeToVert[1]).coord();
        const long double length=f.horizontal?std::abs(b.x-a.x):std::abs(b.y-a.y);
        const long double area=length*(rz?CartesianMeasure::twoPi*(f.horizontal?(a.x+b.x)/2:a.x):1.L);
        bound(f.area-area,area,peakTrace,"independent actual native stress face measure");
        for(int side=0;side<2;++side){const int ci=side==0?f.negativeCell:f.positiveCell;if(ci<0)continue;
            const int p=l.pressureDof(ci);(f.horizontal?dy[p]:dx[p])[u]+=(side==0?1.L:-1.L)*area;}
        if(f.negativeCell>=0&&f.positiveCell>=0){const auto aSize=mesh.cells(mesh.cellIndexToLabel[f.negativeCell]).size();
            const auto bSize=mesh.cells(mesh.cellIndexToLabel[f.positiveCell]).size();if(aSize.x!=bSize.x||aSize.y!=bSize.y)++hanging;}
    }
    check(hanging>0,"stress geometry used no actual hanging face");
    for(std::size_t q=0;q<points.size();++q){const int p=points[q].cellDof;const long double w=points[q].physicalWeight;weights[p]+=w;
        for(int k=0;k<nc;++k){const auto row=static_cast<std::size_t>(p)*nc+k;
            add(actual[row],strain.matrix(),static_cast<Eigen::Index>(q)*nc+k,w);
            for(Matrix::InnerIterator e(rawStrain.matrix(),static_cast<Eigen::Index>(q)*nc+k);e;++e)
                rawScale[row][static_cast<std::size_t>(e.col())]+=w*std::abs(e.value());}
        if(rz){check(points[q].point.x>0,"actual volume Gauss point on/below RZ axis");
            add(hr[p],vv,static_cast<Eigen::Index>(q)*2,w/points[q].point.x);
            add(hz[p],vv,static_cast<Eigen::Index>(q)*2+1,w/points[q].point.x);}
    }
    for(std::size_t q=0;q<facePoints.size();++q){const auto& s=facePoints[q];check(s.velocityDof>=0&&s.velocityDof<nu,"stress trace foreign native velocity");
        ++pointCounts[s.velocityDof];const auto f=l.domain().face(mesh,s.nativeFace);check(s.nativeFace==l.nativeFace(s.velocityDof),"stress trace native row identity");
        const auto& edge=f.horizontal?mesh.hEdges(f.edgeLabel):mesh.vEdges(f.edgeLabel);
        const auto a=mesh.verts(edge.edgeToVert[0]).coord(),b=mesh.verts(edge.edgeToVert[1]).coord();
        const long double length=f.horizontal?std::abs(b.x-a.x):std::abs(b.y-a.y);
        const long double weight=length/2*(rz?CartesianMeasure::twoPi*s.point.x:1.L);
        bound(s.physicalWeight-weight,weight,peakTrace,"independent two-point physical stress face quadrature");
        if(f.area==0){check(rz&&!f.horizontal&&s.point.x==0&&s.physicalWeight==0,"stress axis has fake radius/weight");++axisPoints;}
        const int tangent=f.horizontal?0:1;
        check(s.negativeCellDof==(f.negativeCell<0?-1:l.pressureDof(f.negativeCell))&&
            s.positiveCellDof==(f.positiveCell<0?-1:l.pressureDof(f.positiveCell)),"stress shared trace incidence changed");
        for(int side=0;side<2;++side){const int p=side==0?s.negativeCellDof:s.positiveCellDof;if(p<0)continue;
            add(target[static_cast<std::size_t>(p)*nc+nc-1],fv,static_cast<Eigen::Index>(q)*2+tangent,(side==0?1.L:-1.L)*weight/root2);}
    }
    check(std::all_of(pointCounts.begin(),pointCounts.end(),[](int n){return n==2;}),"native stress face has no unique shared two-point rule");
    if(rz)check(axisPoints>0,"stress RZ witness retained no actual axis");
    const auto& report=strain.stressCompatibilityReport();check(report.cellQuadratureWeightSums.size()==static_cast<std::size_t>(pc)&&
        report.cellQuadratureWeightDefects.size()==static_cast<std::size_t>(pc),"stress true weight report shape");
    for(int p=0;p<pc;++p){const auto base=static_cast<std::size_t>(p)*nc;target[base]=dx[p];target[base+1]=dy[p];
        if(rz)for(int u=0;u<nu;++u){target[base][u]-=hr[p][u];target[base+2][u]=hr[p][u];target[base+3][u]-=hz[p][u]/root2;}
        check(report.cellQuadratureWeightSums[p]==static_cast<double>(weights[p])&&
            report.cellQuadratureWeightDefects[p]==static_cast<double>(weights[p]-metrics.cellVolumes()[p]),"stress W-nativeC was not reported from actual quadrature");
    }
    const Eigen::MatrixXd declared=integral;
    for(std::size_t row=0;row<target.size();++row){long double scale=0;for(int u=0;u<nu;++u)scale+=std::abs(target[row][u])+rawScale[row][u];
        for(int u=0;u<nu;++u){bound(actual[row][u]-target[row][u],scale,peakMoment,"stress per-component conservative mean differs from independent native moment");
            bound(declared(static_cast<Eigen::Index>(row),u)-target[row][u],scale,peakMoment,"stress declared component moment differs from actual geometry");}}
    // Independent P2 oracle verifies that exposed value/trace maps use the
    // original polynomial fit, with no unit-vector probing or new fit.
    const auto ref=reference(mesh,graph,l);const auto poly=values(l,ref,[&](Vector2D p){return polynomial(p,rz,0);});
    const Eigen::Map<const Eigen::VectorXd> pu(poly.normalVelocity().data(),nu);const Eigen::VectorXd pv=vv*pu,pf=fv*pu;
    for(std::size_t q=0;q<points.size();++q){const auto exact=polynomial(points[q].point,rz,0).velocity;
        const long double scale=std::abs(exact.x)+std::abs(exact.y);
        bound(pv[2*q]-exact.x,scale,peakTrace,"stress volume radial/x value is not original P2");
        bound(pv[2*q+1]-exact.y,scale,peakTrace,"stress volume axial/y value is not original P2");}
    for(std::size_t q=0;q<facePoints.size();++q){const auto& s=facePoints[q];const bool vertical=l.normalAxis(s.velocityDof)==HydroDofLayout::NormalAxis::First;
        const int normal=vertical?0:1,tangent=1-normal;const auto exact=polynomial(s.point,rz,0).velocity;
        check(pf[2*q+normal]==poly.normalVelocity()[s.velocityDof],"stress normal trace is not exact native mean");
        bound(pf[2*q+tangent]-(tangent==0?exact.x:exact.y),std::abs(exact.x)+std::abs(exact.y),peakTrace,"stress shared tangent is not original P2 value");}
    if(rz)for(int p=0;p<pc;++p){const auto& c=mesh.cells(mesh.cellIndexToLabel[l.nativeCell(p)]);const auto center=c.centerCoord(),size=c.size();
        const long double a=center.x-size.x/2,b=center.x+size.x/2,z0=center.y-size.y/2,z1=center.y+size.y/2;
        const long double volume=CartesianMeasure::twoPi*(b*b-a*a)*(z1-z0)/2;
        const long double expectedR=volume*(1+2*center.y),expectedZ=CartesianMeasure::twoPi*(b-a)*(z1-z0)*
            (3+(a*a+a*b+b*b)/3+center.y+2*(z0*z0+z0*z1+z1*z1)/3);
        long double r=0,z=0;for(int u=0;u<nu;++u){r+=hr[p][u]*poly.normalVelocity()[u];z+=hz[p][u]*poly.normalVelocity()[u];}
        bound(r-expectedR,std::abs(expectedR),peakTrace,"independent RZ raw u_r/r moment");
        bound(z-expectedZ,std::abs(expectedZ),peakTrace,"independent RZ raw u_z/r shear moment");}
    // Constant admissible physical tensors: arbitrary XY tensor, or a true
    // axis-regular constant axial electric-field tensor (Trr=Ttheta) in RZ.
    const double txx=rz?-2.:2.,tyy=rz?2.:-1.,txy=rz?0.:.37,hoop=-2.;
    Eigen::VectorXd wt=Eigen::VectorXd::Zero(strain.matrix().rows()),bt=Eigen::VectorXd::Zero(fv.rows());std::vector<long double> forceScale(nu);
    for(std::size_t q=0;q<points.size();++q){const std::array<double,4> tensor{txx,tyy,rz?hoop:std::sqrt(2.)*txy,std::sqrt(2.)*txy};
        for(int k=0;k<nc;++k){const auto row=static_cast<Eigen::Index>(q)*nc+k;wt[row]=points[q].physicalWeight*tensor[k];
            for(Matrix::InnerIterator e(strain.matrix(),row);e;++e)forceScale[e.col()]+=std::abs(static_cast<long double>(e.value())*wt[row]);}}
    for(std::size_t q=0;q<facePoints.size();++q){const auto& s=facePoints[q];if(s.negativeCellDof>=0&&s.positiveCellDof>=0)continue;
        const auto f=l.domain().face(mesh,s.nativeFace);const double sign=s.negativeCellDof>=0?1.:-1.;
        bt[2*q]=s.physicalWeight*sign*(f.horizontal?txy:txx);bt[2*q+1]=s.physicalWeight*sign*(f.horizontal?tyy:txy);
        for(int axis=0;axis<2;++axis)for(Matrix::InnerIterator e(fv,static_cast<Eigen::Index>(q)*2+axis);e;++e)
            forceScale[e.col()]+=std::abs(static_cast<long double>(e.value())*bt[2*q+axis]);}
    const Eigen::VectorXd constantForce=-strain.matrix().transpose()*wt+fv.transpose()*bt;
    for(int u=0;u<nu;++u)bound(constantForce[u],forceScale[u],peakConstant,"constant physical tensor plus same boundary traction has native force");
    // Arbitrary isotropic cell stress must reproduce the ORIGINAL C=-DI^T
    // before boundary elimination; no separately assembled pressure stencil.
    std::vector<double> pressure(pc);for(int p=0;p<pc;++p)pressure[p]=.13+std::sin(.2*p);
    const HydroNativeValues pressureValues(l,std::vector<double>(nu),pressure);const auto cp=HydroPressureLoad::build(metrics.divergence()).apply(mesh,pressureValues);
    std::fill(forceScale.begin(),forceScale.end(),0);wt.setZero();
    for(std::size_t q=0;q<points.size();++q)for(int k=0;k<(rz?3:2);++k){const auto row=static_cast<Eigen::Index>(q)*nc+k;
        wt[row]=points[q].physicalWeight*pressure[points[q].cellDof];for(Matrix::InnerIterator e(strain.matrix(),row);e;++e)
            forceScale[e.col()]+=std::abs(static_cast<long double>(e.value())*wt[row]);}
    const Eigen::VectorXd isotropic=-strain.matrix().transpose()*wt;
    for(int u=0;u<nu;++u)bound(isotropic[u]-cp[u],forceScale[u]+std::abs(cp[u]),peakPressure,"stress isotropic load differs from original C");
    // Actual native virtual work, with variable tensor and the SAME boundary
    // trace; no boundary row is silently dropped or called an interior force.
    std::vector<double> random(nu);for(int u=0;u<nu;++u)random[u]=std::sin(.37*u)+std::cos(.91*u);
    for(int u:metrics.axisVelocityDofs())random[u]=0;
    const HydroNativeValues velocity(l,random,std::vector<double>(pc));const auto strainValue=strain.apply(mesh,velocity);
    const Eigen::Map<const Eigen::VectorXd> eu(random.data(),nu);const Eigen::VectorXd faceValue=fv*eu;
    for(Eigen::Index row=0;row<wt.size();++row)wt[row]=points[static_cast<std::size_t>(row)/nc].physicalWeight*std::cos(.17*row);
    for(Eigen::Index row=0;row<bt.size();++row)bt[row]*=std::sin(.23*row); // Only the actual physical-boundary rows are active.
    const Eigen::VectorXd force=-strain.matrix().transpose()*wt+fv.transpose()*bt;long double direct=0,independent=0,workScale=0;
    for(int u=0;u<nu;++u){const long double term=static_cast<long double>(random[u])*force[u];direct+=term;workScale+=std::abs(term);}
    for(std::size_t row=0;row<strainValue.size();++row){const long double term=static_cast<long double>(wt[row])*strainValue[row];independent-=term;workScale+=std::abs(term);}
    for(Eigen::Index row=0;row<bt.size();++row){const long double term=static_cast<long double>(bt[row])*faceValue[row];independent+=term;workScale+=std::abs(term);}
    bound(direct-independent,workScale,peakWork,"stress native virtual work with boundary trace failed");
    const auto zero=HydroEnergyViscosity::zero(mesh,metrics);check(zero.matrix().nonZeros()==0&&zero.strain().matrix().rows()==0&&
        zero.strain().points().empty()&&zero.dissipation(mesh,velocity)==0,"stress geometry inserted fake zero-viscosity strain/work");
    for(double f:zero.apply(mesh,velocity))check(f==0,"stress geometry inserted fake zero-viscosity load");
    std::cout<<"STRESS_GEOMETRY,"<<(rz?"RZ":"XY")<<",cells,"<<pc<<",hanging_faces,"<<hanging<<",axis_points,"<<axisPoints
        <<",face_points,"<<facePoints.size()<<",moment,"<<peakMoment<<",trace,"<<peakTrace<<",constant_tensor,"<<peakConstant
        <<",same_C,"<<peakPressure<<",virtual_work,"<<peakWork<<",W_nativeC,"<<report.maximumRelativeQuadratureWeightDefect<<'\n';
}
void run(){scene(false);scene(true);std::cout<<"HYDRO_STRESS_GEOMETRY_PASS,XY_RZ_true2to1,axis,component_moments,same_C,boundary,virtualwork,mu0_exact0,no_Maxwell\n";}
} // namespace StressCompatibleContract



namespace StressAmrLocalReuseContract {
using Matrix=HydroStrainOperator::Matrix;
using LAMP::Mesh::NativeLeafOverlap;
using LAMP::Mesh::NativeLeafOverlapPiece;
using LAMP::Mesh::NativeTargetLevelField;
using StressCompatibleContract::DenseRows;
using StressCompatibleContract::bound;
using StressCompatibleContract::add;
std::size_t scalarChecks=0,matrixEntries=0,reusedRaw=0,reusedStress=0,reusedFaces=0,fitWork=0,changedPairs=0;
std::uint64_t bits(double value) {std::uint64_t result=0;std::memcpy(&result,&value,sizeof result);return result;}
void exact(double actual,double expected,const char* why) {check(bits(actual)==bits(expected),why);++scalarChecks;}
void matrix(const Matrix& actual,const Matrix& expected,const char* why) {
    check(actual.rows()==expected.rows()&&actual.cols()==expected.cols()&&actual.isCompressed()&&expected.isCompressed()
        &&actual.nonZeros()==expected.nonZeros(),why);
    for(Eigen::Index row=0;row<=actual.outerSize();++row)check(actual.outerIndexPtr()[row]==expected.outerIndexPtr()[row],why);
    for(Eigen::Index i=0;i<actual.nonZeros();++i){check(actual.innerIndexPtr()[i]==expected.innerIndexPtr()[i],why);
        exact(actual.valuePtr()[i],expected.valuePtr()[i],why);++matrixEntries;}
}
void same(const HydroStrainOperator& actual,const HydroStrainOperator& expected) {
    check(actual.isStressCompatible()&&actual.hasCanonicalRawStrain()&&expected.hasCanonicalRawStrain()
        &&actual.layout().sameIdentity(expected.layout())&&actual.components()==expected.components(),"AMR stress target qualification differs from cold target");
    check(&actual.metrics().divergence().matrix()==&expected.metrics().divergence().matrix(),"AMR stress target lost actual fresh D owner");
    matrix(actual.matrix(),expected.matrix(),"AMR Ehat CSR differs from independent cold target bits");
    matrix(actual.canonicalRawStrain().matrix(),expected.canonicalRawStrain().matrix(),"AMR raw CSR differs from independent cold target bits");
    matrix(actual.stressVolumeVelocityMatrix(),expected.stressVolumeVelocityMatrix(),"AMR volume velocity CSR differs from cold target bits");
    matrix(actual.stressFaceVelocityMatrix(),expected.stressFaceVelocityMatrix(),"AMR face trace CSR differs from cold target bits");
    matrix(actual.stressComponentIntegralMatrix(),expected.stressComponentIntegralMatrix(),"AMR component moment CSR differs from cold target bits");
    const auto& a=actual.points();const auto& b=expected.points();check(a.size()==b.size(),"AMR actual quadrature count changed");
    for(std::size_t q=0;q<a.size();++q){check(a[q].cellDof==b[q].cellDof,"AMR actual volume quadrature owner changed");
        exact(a[q].point.x,b[q].point.x,"AMR actual volume x changed");exact(a[q].point.y,b[q].point.y,"AMR actual volume y changed");
        exact(a[q].physicalWeight,b[q].physicalWeight,"AMR actual volume measure changed");}
    const auto& af=actual.stressFacePoints();const auto& bf=expected.stressFacePoints();check(af.size()==bf.size(),"AMR actual face quadrature count changed");
    for(std::size_t q=0;q<af.size();++q){check(af[q].velocityDof==bf[q].velocityDof&&af[q].nativeFace==bf[q].nativeFace
        &&af[q].negativeCellDof==bf[q].negativeCellDof&&af[q].positiveCellDof==bf[q].positiveCellDof,"AMR face quadrature actual incidence changed");
        exact(af[q].point.x,bf[q].point.x,"AMR actual face x changed");exact(af[q].point.y,bf[q].point.y,"AMR actual face y changed");
        exact(af[q].physicalWeight,bf[q].physicalWeight,"AMR actual face measure changed");}
    const auto& ar=actual.stressCompatibilityReport();const auto& br=expected.stressCompatibilityReport();
    check(ar.cellQuadratureWeightSums.size()==br.cellQuadratureWeightSums.size()&&ar.cellQuadratureWeightDefects.size()==br.cellQuadratureWeightDefects.size(),"AMR weight report dimensions changed");
    for(std::size_t p=0;p<ar.cellQuadratureWeightSums.size();++p){exact(ar.cellQuadratureWeightSums[p],br.cellQuadratureWeightSums[p],"AMR weight report sum changed");
        exact(ar.cellQuadratureWeightDefects[p],br.cellQuadratureWeightDefects[p],"AMR weight report defect changed");}
    exact(ar.maximumRelativeQuadratureWeightDefect,br.maximumRelativeQuadratureWeightDefect,"AMR max volume weight defect changed");
    exact(ar.maximumRelativeFaceQuadratureWeightDefect,br.maximumRelativeFaceQuadratureWeightDefect,"AMR max face weight defect changed");
}
bool insideBaseBox(const QuadtreeMesh& mesh,const Cell& cell) {
    const auto x=cell.centerCoord(),h=cell.size(),w=mesh.maxCoord-mesh.minCoord;
    const auto a=mesh.minCoord+w*.25,b=mesh.minCoord+w*.375;
    return std::max(x.x-h.x/2,a.x)<std::min(x.x+h.x/2,b.x)&&std::max(x.y-h.y/2,a.y)<std::min(x.y+h.y/2,b.y);
}
void coarsenAddedRegion(QuadtreeMesh& mesh) {
    int sweeps=0;
    for(;;){NativeTargetLevelField targets(mesh,0.),gradient(mesh,0.);
        for(int p=0;p<mesh.cellIndexNum();++p)targets.value[p]=insideBaseBox(mesh,mesh.cells(mesh.cellIndexToLabel[p]))?1.:0.;
        mesh.completeRefinementPrerequisites(targets);if(!mesh.mayAdaptToTargets(targets))break;
        check(sweeps++<mesh.maxLevel()+2,"AMR stress real coarsening did not settle");
        std::vector<VertField<double>> levelSets;std::vector<VertField<double>*> vertices;
        std::vector<HEdgeField<double>*> horizontal;std::vector<VEdgeField<double>*> vertical;
        std::vector<NativeTargetLevelField*> fields{&targets};mesh.adapt(levelSets,gradient,vertices,horizontal,vertical,fields,&targets);}
    check(sweeps>0,"AMR stress fixture skipped actual coarsening operation");
}
HydroPhysicalDivergence metrics(const QuadtreeMesh& mesh,bool rz) {
    FaceSegmentGraph graph;graph.build(mesh);const auto layout=HydroDofLayout::build(mesh,NativeLeafDomain::build(mesh,graph,rz));
    return HydroPhysicalDivergence::build(mesh,HydroIntegralDivergence::build(mesh,layout));
}
void physical(const QuadtreeMesh& mesh,const FaceSegmentGraph& graph,const HydroStrainOperator& strain) {
    const auto& metrics=strain.metrics();const auto& l=strain.layout();const bool rz=l.domain().axisymmetric();
    const auto points=HydroStrainOperator::volumeQuadrature(mesh,metrics);
    const auto rawStrain=HydroStrainOperator::build(mesh,metrics,points);
    const auto rawRef=reference(mesh,graph,l);const auto rawPoly=values(l,rawRef,[rz](Vector2D p){return polynomial(p,rz,0);});
    const auto rawValues=strain.canonicalRawStrain().apply(mesh,rawPoly);
    for(std::size_t q=0;q<points.size();++q){const auto expected=polynomial(points[q].point,rz,0).strain;
        for(int k=0;k<strain.components();++k)close(rawValues[q*strain.components()+k]-expected[k],std::abs(expected[k]),peakPolynomial);}
    const auto& facePoints=strain.stressFacePoints();const auto& fv=strain.stressFaceVelocityMatrix();
    const auto& vv=strain.stressVolumeVelocityMatrix();const auto& integral=strain.stressComponentIntegralMatrix();
    const int pc=l.pressureCount(),nu=l.velocityCount(),nc=strain.components();const long double root2=std::sqrt(2.L);
    check(strain.points().size()==points.size()&&facePoints.size()==2*static_cast<std::size_t>(nu)&&
        fv.rows()==static_cast<Eigen::Index>(facePoints.size())*2&&vv.rows()==static_cast<Eigen::Index>(points.size())*2&&
        integral.rows()==pc*nc&&fv.cols()==nu&&vv.cols()==nu&&integral.cols()==nu,"stress geometry native shape mismatch");
    DenseRows target(static_cast<std::size_t>(pc)*nc,std::vector<long double>(nu)),actual=target,rawScale=target;
    DenseRows dx(pc,std::vector<long double>(nu)),dy=dx,hr=dx,hz=dx;
    std::vector<long double> weights(pc);std::vector<int> pointCounts(nu);std::size_t hanging=0,axisPoints=0;
    double peakMoment=0,peakTrace=0,peakConstant=0,peakPressure=0,peakWork=0;
    for(int u=0;u<nu;++u){const auto f=l.domain().face(mesh,l.nativeFace(u));
        const auto& edge=f.horizontal?mesh.hEdges(f.edgeLabel):mesh.vEdges(f.edgeLabel);
        const auto a=mesh.verts(edge.edgeToVert[0]).coord(),b=mesh.verts(edge.edgeToVert[1]).coord();
        const long double length=f.horizontal?std::abs(b.x-a.x):std::abs(b.y-a.y);
        const long double area=length*(rz?CartesianMeasure::twoPi*(f.horizontal?(a.x+b.x)/2:a.x):1.L);
        bound(f.area-area,area,peakTrace,"independent actual native stress face measure");
        for(int side=0;side<2;++side){const int ci=side==0?f.negativeCell:f.positiveCell;if(ci<0)continue;
            const int p=l.pressureDof(ci);(f.horizontal?dy[p]:dx[p])[u]+=(side==0?1.L:-1.L)*area;}
        if(f.negativeCell>=0&&f.positiveCell>=0){const auto aSize=mesh.cells(mesh.cellIndexToLabel[f.negativeCell]).size();
            const auto bSize=mesh.cells(mesh.cellIndexToLabel[f.positiveCell]).size();if(aSize.x!=bSize.x||aSize.y!=bSize.y)++hanging;}
    }
    check(hanging>0,"stress geometry used no actual hanging face");
    for(std::size_t q=0;q<points.size();++q){const int p=points[q].cellDof;const long double w=points[q].physicalWeight;weights[p]+=w;
        for(int k=0;k<nc;++k){const auto row=static_cast<std::size_t>(p)*nc+k;
            add(actual[row],strain.matrix(),static_cast<Eigen::Index>(q)*nc+k,w);
            for(Matrix::InnerIterator e(rawStrain.matrix(),static_cast<Eigen::Index>(q)*nc+k);e;++e)
                rawScale[row][static_cast<std::size_t>(e.col())]+=w*std::abs(e.value());}
        if(rz){check(points[q].point.x>0,"actual volume Gauss point on/below RZ axis");
            add(hr[p],vv,static_cast<Eigen::Index>(q)*2,w/points[q].point.x);
            add(hz[p],vv,static_cast<Eigen::Index>(q)*2+1,w/points[q].point.x);}
    }
    for(std::size_t q=0;q<facePoints.size();++q){const auto& s=facePoints[q];check(s.velocityDof>=0&&s.velocityDof<nu,"stress trace foreign native velocity");
        ++pointCounts[s.velocityDof];const auto f=l.domain().face(mesh,s.nativeFace);check(s.nativeFace==l.nativeFace(s.velocityDof),"stress trace native row identity");
        const auto& edge=f.horizontal?mesh.hEdges(f.edgeLabel):mesh.vEdges(f.edgeLabel);
        const auto a=mesh.verts(edge.edgeToVert[0]).coord(),b=mesh.verts(edge.edgeToVert[1]).coord();
        const long double length=f.horizontal?std::abs(b.x-a.x):std::abs(b.y-a.y);
        const long double weight=length/2*(rz?CartesianMeasure::twoPi*s.point.x:1.L);
        bound(s.physicalWeight-weight,weight,peakTrace,"independent two-point physical stress face quadrature");
        if(f.area==0){check(rz&&!f.horizontal&&s.point.x==0&&s.physicalWeight==0,"stress axis has fake radius/weight");++axisPoints;}
        const int tangent=f.horizontal?0:1;
        check(s.negativeCellDof==(f.negativeCell<0?-1:l.pressureDof(f.negativeCell))&&
            s.positiveCellDof==(f.positiveCell<0?-1:l.pressureDof(f.positiveCell)),"stress shared trace incidence changed");
        for(int side=0;side<2;++side){const int p=side==0?s.negativeCellDof:s.positiveCellDof;if(p<0)continue;
            add(target[static_cast<std::size_t>(p)*nc+nc-1],fv,static_cast<Eigen::Index>(q)*2+tangent,(side==0?1.L:-1.L)*weight/root2);}
    }
    check(std::all_of(pointCounts.begin(),pointCounts.end(),[](int n){return n==2;}),"native stress face has no unique shared two-point rule");
    if(rz&&mesh.minCoord.x==0)check(axisPoints>0,"stress RZ witness retained no actual axis");
    if(rz&&mesh.minCoord.x>0)check(axisPoints==0,"annular stress geometry invented an axis");
    const auto& report=strain.stressCompatibilityReport();check(report.cellQuadratureWeightSums.size()==static_cast<std::size_t>(pc)&&
        report.cellQuadratureWeightDefects.size()==static_cast<std::size_t>(pc),"stress true weight report shape");
    for(int p=0;p<pc;++p){const auto base=static_cast<std::size_t>(p)*nc;target[base]=dx[p];target[base+1]=dy[p];
        if(rz)for(int u=0;u<nu;++u){target[base][u]-=hr[p][u];target[base+2][u]=hr[p][u];target[base+3][u]-=hz[p][u]/root2;}
        check(report.cellQuadratureWeightSums[p]==static_cast<double>(weights[p])&&
            report.cellQuadratureWeightDefects[p]==static_cast<double>(weights[p]-metrics.cellVolumes()[p]),"stress W-nativeC was not reported from actual quadrature");
    }
    const Eigen::MatrixXd declared=integral;
    for(std::size_t row=0;row<target.size();++row){long double scale=0;for(int u=0;u<nu;++u)scale+=std::abs(target[row][u])+rawScale[row][u];
        for(int u=0;u<nu;++u){bound(actual[row][u]-target[row][u],scale,peakMoment,"stress per-component conservative mean differs from independent native moment");
            bound(declared(static_cast<Eigen::Index>(row),u)-target[row][u],scale,peakMoment,"stress declared component moment differs from actual geometry");}}
    // Independent P2 oracle verifies that exposed value/trace maps use the
    // original polynomial fit, with no unit-vector probing or new fit.
    const auto ref=reference(mesh,graph,l);const auto poly=values(l,ref,[&](Vector2D p){return polynomial(p,rz,0);});
    const Eigen::Map<const Eigen::VectorXd> pu(poly.normalVelocity().data(),nu);const Eigen::VectorXd pv=vv*pu,pf=fv*pu;
    for(std::size_t q=0;q<points.size();++q){const auto exact=polynomial(points[q].point,rz,0).velocity;
        const long double scale=std::abs(exact.x)+std::abs(exact.y);
        bound(pv[2*q]-exact.x,scale,peakTrace,"stress volume radial/x value is not original P2");
        bound(pv[2*q+1]-exact.y,scale,peakTrace,"stress volume axial/y value is not original P2");}
    for(std::size_t q=0;q<facePoints.size();++q){const auto& s=facePoints[q];const bool vertical=l.normalAxis(s.velocityDof)==HydroDofLayout::NormalAxis::First;
        const int normal=vertical?0:1,tangent=1-normal;const auto exact=polynomial(s.point,rz,0).velocity;
        check(pf[2*q+normal]==poly.normalVelocity()[s.velocityDof],"stress normal trace is not exact native mean");
        bound(pf[2*q+tangent]-(tangent==0?exact.x:exact.y),std::abs(exact.x)+std::abs(exact.y),peakTrace,"stress shared tangent is not original P2 value");}
    if(rz)for(int p=0;p<pc;++p){const auto& c=mesh.cells(mesh.cellIndexToLabel[l.nativeCell(p)]);const auto center=c.centerCoord(),size=c.size();
        const long double a=center.x-size.x/2,b=center.x+size.x/2,z0=center.y-size.y/2,z1=center.y+size.y/2;
        const long double volume=CartesianMeasure::twoPi*(b*b-a*a)*(z1-z0)/2;
        const long double expectedR=volume*(1+2*center.y),expectedZ=CartesianMeasure::twoPi*(b-a)*(z1-z0)*
            (3+(a*a+a*b+b*b)/3+center.y+2*(z0*z0+z0*z1+z1*z1)/3);
        long double r=0,z=0;for(int u=0;u<nu;++u){r+=hr[p][u]*poly.normalVelocity()[u];z+=hz[p][u]*poly.normalVelocity()[u];}
        bound(r-expectedR,std::abs(expectedR),peakTrace,"independent RZ raw u_r/r moment");
        bound(z-expectedZ,std::abs(expectedZ),peakTrace,"independent RZ raw u_z/r shear moment");}
    // Constant admissible physical tensors: arbitrary XY tensor, or a true
    // axis-regular constant axial electric-field tensor (Trr=Ttheta) in RZ.
    const double txx=rz?-2.:2.,tyy=rz?2.:-1.,txy=rz?0.:.37,hoop=-2.;
    Eigen::VectorXd wt=Eigen::VectorXd::Zero(strain.matrix().rows()),bt=Eigen::VectorXd::Zero(fv.rows());std::vector<long double> forceScale(nu);
    for(std::size_t q=0;q<points.size();++q){const std::array<double,4> tensor{txx,tyy,rz?hoop:std::sqrt(2.)*txy,std::sqrt(2.)*txy};
        for(int k=0;k<nc;++k){const auto row=static_cast<Eigen::Index>(q)*nc+k;wt[row]=points[q].physicalWeight*tensor[k];
            for(Matrix::InnerIterator e(strain.matrix(),row);e;++e)forceScale[e.col()]+=std::abs(static_cast<long double>(e.value())*wt[row]);}}
    for(std::size_t q=0;q<facePoints.size();++q){const auto& s=facePoints[q];if(s.negativeCellDof>=0&&s.positiveCellDof>=0)continue;
        const auto f=l.domain().face(mesh,s.nativeFace);const double sign=s.negativeCellDof>=0?1.:-1.;
        bt[2*q]=s.physicalWeight*sign*(f.horizontal?txy:txx);bt[2*q+1]=s.physicalWeight*sign*(f.horizontal?tyy:txy);
        for(int axis=0;axis<2;++axis)for(Matrix::InnerIterator e(fv,static_cast<Eigen::Index>(q)*2+axis);e;++e)
            forceScale[e.col()]+=std::abs(static_cast<long double>(e.value())*bt[2*q+axis]);}
    const Eigen::VectorXd constantForce=-strain.matrix().transpose()*wt+fv.transpose()*bt;
    for(int u=0;u<nu;++u)bound(constantForce[u],forceScale[u],peakConstant,"constant physical tensor plus same boundary traction has native force");
    // Arbitrary isotropic cell stress must reproduce the ORIGINAL C=-DI^T
    // before boundary elimination; no separately assembled pressure stencil.
    std::vector<double> pressure(pc);for(int p=0;p<pc;++p)pressure[p]=.13+std::sin(.2*p);
    const HydroNativeValues pressureValues(l,std::vector<double>(nu),pressure);const auto cp=HydroPressureLoad::build(metrics.divergence()).apply(mesh,pressureValues);
    std::fill(forceScale.begin(),forceScale.end(),0);wt.setZero();
    for(std::size_t q=0;q<points.size();++q)for(int k=0;k<(rz?3:2);++k){const auto row=static_cast<Eigen::Index>(q)*nc+k;
        wt[row]=points[q].physicalWeight*pressure[points[q].cellDof];for(Matrix::InnerIterator e(strain.matrix(),row);e;++e)
            forceScale[e.col()]+=std::abs(static_cast<long double>(e.value())*wt[row]);}
    const Eigen::VectorXd isotropic=-strain.matrix().transpose()*wt;
    for(int u=0;u<nu;++u)bound(isotropic[u]-cp[u],forceScale[u]+std::abs(cp[u]),peakPressure,"stress isotropic load differs from original C");
    // Actual native virtual work, with variable tensor and the SAME boundary
    // trace; no boundary row is silently dropped or called an interior force.
    std::vector<double> random(nu);for(int u=0;u<nu;++u)random[u]=std::sin(.37*u)+std::cos(.91*u);
    for(int u:metrics.axisVelocityDofs())random[u]=0;
    const HydroNativeValues velocity(l,random,std::vector<double>(pc));const auto strainValue=strain.apply(mesh,velocity);
    const Eigen::Map<const Eigen::VectorXd> eu(random.data(),nu);const Eigen::VectorXd faceValue=fv*eu;
    for(Eigen::Index row=0;row<wt.size();++row)wt[row]=points[static_cast<std::size_t>(row)/nc].physicalWeight*std::cos(.17*row);
    for(Eigen::Index row=0;row<bt.size();++row)bt[row]*=std::sin(.23*row); // Only the actual physical-boundary rows are active.
    const Eigen::VectorXd force=-strain.matrix().transpose()*wt+fv.transpose()*bt;long double direct=0,independent=0,workScale=0;
    for(int u=0;u<nu;++u){const long double term=static_cast<long double>(random[u])*force[u];direct+=term;workScale+=std::abs(term);}
    for(std::size_t row=0;row<strainValue.size();++row){const long double term=static_cast<long double>(wt[row])*strainValue[row];independent-=term;workScale+=std::abs(term);}
    for(Eigen::Index row=0;row<bt.size();++row){const long double term=static_cast<long double>(bt[row])*faceValue[row];independent+=term;workScale+=std::abs(term);}
    bound(direct-independent,workScale,peakWork,"stress native virtual work with boundary trace failed");
    const auto zero=HydroEnergyViscosity::zero(mesh,metrics);check(zero.matrix().nonZeros()==0&&zero.strain().matrix().rows()==0&&
        zero.strain().points().empty()&&zero.dissipation(mesh,velocity)==0,"stress geometry inserted fake zero-viscosity strain/work");
    for(double f:zero.apply(mesh,velocity))check(f==0,"stress geometry inserted fake zero-viscosity load");
    std::cout<<"STRESS_GEOMETRY,"<<(rz?"RZ":"XY")<<",cells,"<<pc<<",hanging_faces,"<<hanging<<",axis_points,"<<axisPoints
        <<",face_points,"<<facePoints.size()<<",moment,"<<peakMoment<<",trace,"<<peakTrace<<",constant_tensor,"<<peakConstant
        <<",same_C,"<<peakPressure<<",virtual_work,"<<peakWork<<",W_nativeC,"<<report.maximumRelativeQuadratureWeightDefect<<'\n';

}

// Task41 is opt-in to the SAME original six actual AMR scenes/oracles below.
// No new mesh, source certificate or test-only direct/fallback switch exists.
namespace RawVolumeCsrRebind
{
void absent(const HydroStrainOperator&);
void snapshots(const HydroStressGeometryReuseReport&,const HydroStressGeometryReuseReport&,const HydroStressGeometryReuseReport&);
}
namespace CertifiedCsrRebind
{
bool enabled=false;
std::size_t scenes=0,directRows=0,literalRows=0,publicOrderedRows=0,publicMissingRows=0,publicUnorderedRows=0,fpPairs=0,coldFallbacks=0;
using FaceKey=std::array<std::uint64_t,5>;
FaceKey faceKey(const QuadtreeMesh& mesh,const HydroDofLayout& layout,int u)
{
    const auto f=layout.domain().face(mesh,layout.nativeFace(u));const auto& e=f.horizontal?mesh.hEdges(f.edgeLabel):mesh.vEdges(f.edgeLabel);
    auto a=mesh.verts(e.edgeToVert[0]).coord(),b=mesh.verts(e.edgeToVert[1]).coord();
    if(a.x>b.x||(a.x==b.x&&a.y>b.y))std::swap(a,b);
    return {{f.horizontal?1u:0u,bits(a.x),bits(a.y),bits(b.x),bits(b.y)}};
}
struct Floating
{
    std::fenv_t native{};int error=errno;
#if defined(_MSC_VER) && defined(_M_X64)
    unsigned csr=_mm_getcsr();
#endif
    Floating(){check(std::fegetenv(&native)==0,"CSR rebind oracle cannot capture actual FP environment");}
    void restore()const{check(std::fesetenv(&native)==0,"CSR rebind oracle cannot restore actual FP environment");
#if defined(_MSC_VER) && defined(_M_X64)
        _mm_setcsr(csr);
#endif
        errno=error;
    }
};
void sameFp(const Floating& a,const Floating& b)
{
    check(std::memcmp(&a.native,&b.native,sizeof a.native)==0&&a.error==b.error,"CSR repeated actual factory changed full FP environment/errno");
#if defined(_MSC_VER) && defined(_M_X64)
    check(a.csr==b.csr,"CSR repeated actual factory changed raw MXCSR status/controls");
#endif
    ++fpPairs;
}
void absent(const HydroStrainOperator& geometry)
{
    RawVolumeCsrRebind::absent(geometry);
    if(!enabled)return;const auto r=geometry.geometryReuseReport();
    check(r.directReboundRows==0&&r.literalReboundRows==0&&r.csrRebindScratchBytes==0,
        "CSR original cold/unsupported source invented copied rows/scratch");++coldFallbacks;
}
void inspect(const QuadtreeMesh& sourceMesh,const HydroStrainOperator& source,
    const QuadtreeMesh& targetMesh,const HydroPhysicalDivergence& targetMetrics,const NativeLeafOverlap& overlap,
    const HydroStrainOperator& cold,const HydroStrainOperator& reused)
{
    if(!enabled)return;absent(cold);const auto r=reused.geometryReuseReport();
    const auto& layout=reused.layout();const auto nc=static_cast<std::size_t>(reused.components());
    std::vector<std::size_t> volume(layout.pressureCount()),face(layout.velocityCount());
    for(const auto& q:reused.points()){check(q.cellDof>=0&&q.cellDof<layout.pressureCount(),"CSR target volume point foreign owner");++volume[q.cellDof];}
    for(const auto& q:reused.stressFacePoints()){check(q.velocityDof>=0&&q.velocityDof<layout.velocityCount(),"CSR target face point foreign owner");++face[q.velocityDof];}
    check(!volume.empty()&&!face.empty()&&std::all_of(volume.begin(),volume.end(),[&](std::size_t n){return n==volume.front();})
        &&std::all_of(face.begin(),face.end(),[&](std::size_t n){return n==face.front();}),"CSR actual canonical scalar-row multiplicities vary");
    // Count actual public quadrature multiplicities first; do not hard-code a
    // planar/RZ component number or call a net leaf difference copied work.
    const auto expected=2*face.front()*r.reusedFaces+nc*(volume.front()+1)*r.reusedStressCells;
    check(r.directReboundRows+r.literalReboundRows==expected&&expected>0,
        "CSR successful copy counters differ from actual face/integral/Ehat scalar-row multiplicities");
    check(r.csrRebindScratchBytes<=64u*1024u*1024u,"CSR representation scratch exceeded actual 64MiB stock limit");
    // Independent public geometry correspondence uses exact physical endpoints,
    // not producer toTarget arrays, source native labels, hashes, or its key.
    std::map<FaceKey,int> byPhysicalFace;
    for(int u=0;u<layout.velocityCount();++u)check(byPhysicalFace.emplace(faceKey(targetMesh,layout,u),u).second,
        "CSR public target faces have duplicate exact physical geometry");
    std::vector<int> mapping(source.layout().velocityCount(),-1),mapped;
    for(int u=0;u<source.layout().velocityCount();++u){const auto found=byPhysicalFace.find(faceKey(sourceMesh,source.layout(),u));
        if(found!=byPhysicalFace.end()){mapping[u]=found->second;mapped.push_back(found->second);}}
    auto unique=mapped;std::sort(unique.begin(),unique.end());check(std::adjacent_find(unique.begin(),unique.end())==unique.end(),
        "CSR true public face correspondence is not injective");
    std::size_t ordered=0,missing=0,unordered=0;
    for(const Matrix* matrix:{&source.stressFaceVelocityMatrix(),&source.stressComponentIntegralMatrix(),&source.matrix()})
        for(Eigen::Index row=0;row<matrix->rows();++row){std::vector<int> columns;bool bound=true,normal=true;
            for(Matrix::InnerIterator e(*matrix,row);e;++e){check(e.col()>=0&&static_cast<std::size_t>(e.col())<mapping.size(),"CSR public source column out of range");
                const int target=mapping[static_cast<std::size_t>(e.col())];if(target<0)bound=false;columns.push_back(target);
                const auto magnitude=bits(e.value())&UINT64_C(0x7fffffffffffffff),exponent=magnitude&UINT64_C(0x7ff0000000000000);
                normal=normal&&!((exponent==0&&magnitude!=0)||exponent==UINT64_C(0x7ff0000000000000));}
            if(!bound){++missing;continue;}
            const bool increasing=std::is_sorted(columns.begin(),columns.end())&&std::adjacent_find(columns.begin(),columns.end())==columns.end();
            if(!increasing){++unordered;continue;}if(normal)++ordered;
        }
    // This is an independent upper bound on possible source-copy rows, not a
    // claim that every geometrically unchanged row passed the private stencil
    // certificate or that the tests can manufacture an invalid private map.
    check(r.directReboundRows<=ordered,"CSR direct counter exceeds all publicly possible strictly ordered normal source rows");
#if defined(_MSC_VER) && defined(_M_X64) && !defined(__FAST_MATH__) && (!defined(_M_FP_FAST) || !_M_FP_FAST) && EIGEN_WORLD_VERSION==3 && EIGEN_MAJOR_VERSION==4 && EIGEN_MINOR_VERSION==1 && !defined(EIGEN_SPARSEMATRIX_PLUGIN) && !defined(EIGEN_SPARSE_COMPRESSEDBASE_PLUGIN)
    check(sizeof(long double)==sizeof(double)&&r.directReboundRows>0&&r.csrRebindScratchBytes>0,
        "CSR actual supported small AMR scene did not execute certified direct rows");
#else
    check(r.directReboundRows==0,"CSR unsupported implementation silently used direct representation");
#endif
    // Repeated real reuse from the SAME complete entry verifies representation
    // stock/flags determinism. It is not described as a disabled-literal API.
    const Floating entry;
    const auto first=HydroStrainOperator::buildStressCompatibleRemapped(sourceMesh,source,targetMesh,targetMetrics,overlap);const Floating firstAfter;
    entry.restore();
    const auto second=HydroStrainOperator::buildStressCompatibleRemapped(sourceMesh,source,targetMesh,targetMetrics,overlap);const Floating secondAfter;
    sameFp(firstAfter,secondAfter);same(first,cold);same(second,cold);
    const auto a=first.geometryReuseReport(),b=second.geometryReuseReport();
    check(a.directReboundRows==r.directReboundRows&&a.literalReboundRows==r.literalReboundRows&&a.csrRebindScratchBytes==r.csrRebindScratchBytes
        &&a.directReboundRows==b.directReboundRows&&a.literalReboundRows==b.literalReboundRows&&a.csrRebindScratchBytes==b.csrRebindScratchBytes,
        "CSR construction snapshot was accumulated or changed across identical actual factory calls");
    RawVolumeCsrRebind::snapshots(r,a,b);
    directRows+=r.directReboundRows;literalRows+=r.literalReboundRows;publicOrderedRows+=ordered;publicMissingRows+=missing;publicUnorderedRows+=unordered;++scenes;
}
void finish()
{
    check(scenes==6&&fpPairs==6&&coldFallbacks>=9&&publicOrderedRows>0&&publicMissingRows>0,
        "CSR six real directional scenes/public source gaps/FP coverage incomplete");
    std::cout<<"STRESS_CERTIFIED_CSR_REBIND_PASS,actual_AMR_scenes,"<<scenes<<",direct_copy_scalar_rows,"<<directRows
        <<",literal_copy_scalar_rows,"<<literalRows<<",public_ordered_normal_source_rows,"<<publicOrderedRows
        <<",public_unbound_source_rows,"<<publicMissingRows<<",public_nonascending_source_rows,"<<publicUnorderedRows
        <<",whole_actual_repeat_FP_errno_pairs,"<<fpPairs<<",cold_or_unsupported_snapshots,"<<coldFallbacks
        <<",private_duplicate_missing_nonascending_copy_map_dynamic_unavailable,optional_resource_cap_not_injected,no_speed_claim\n";
}
}

namespace RawVolumeCsrRebind
{
bool enabled=false;
std::size_t scenes=0,fpPairs=0,snapshotPairs=0,coldFallbacks=0,rawRows=0,volumeRows=0;
std::size_t potentialRawRows=0,potentialVolumeRows=0,orderedRawRows=0,orderedVolumeRows=0;
std::array<std::size_t,2> publicMissingRows{},publicUnorderedRows{},emptySourceRows{};
using Floating=CertifiedCsrRebind::Floating;
// Reuse the two existing actual factories. Capture each whole post-environment
// BEFORE its independent matrix/physics assertions; no synthetic fenv bytes.
struct ColdPair
{
    Matrix sourceVolume;
    Floating entry;
    std::optional<Floating> coldAfter;
    explicit ColdPair(const HydroStrainOperator& source):sourceVolume(source.stressVolumeVelocityMatrix()){}
    void beforeRemap(){coldAfter.emplace();entry.restore();}
    void afterRemap(const HydroStrainOperator& source)
    {
        const Floating actual;
        check(coldAfter.has_value(),"raw/volume CSR pair lost actual cold factory receipt");
        check(std::memcmp(&coldAfter->native,&actual.native,sizeof actual.native)==0&&coldAfter->error==actual.error,
            "raw/volume CSR same-entry cold/remap complete native FP/errno differs");
#if defined(_MSC_VER) && defined(_M_X64)
        check(coldAfter->csr==actual.csr,"raw/volume CSR same-entry cold/remap raw MXCSR differs");
#endif
        ++fpPairs;
        matrix(source.stressVolumeVelocityMatrix(),sourceVolume,"raw/volume CSR remap modified source volume matrix");
    }
};
void absent(const HydroStrainOperator& geometry)
{
    if(!enabled)return;const auto& r=geometry.geometryReuseReport();
    check(r.rawDirectReboundRows==0&&r.volumeDirectReboundRows==0,
        "raw/volume CSR cold or unsupported source invented direct-copy rows");++coldFallbacks;
}
void snapshots(const HydroStressGeometryReuseReport& initial,const HydroStressGeometryReuseReport& first,
    const HydroStressGeometryReuseReport& second)
{
    if(!enabled)return;
    check(initial.rawDirectReboundRows==first.rawDirectReboundRows&&initial.rawDirectReboundRows==second.rawDirectReboundRows
        &&initial.volumeDirectReboundRows==first.volumeDirectReboundRows&&initial.volumeDirectReboundRows==second.volumeDirectReboundRows,
        "raw/volume CSR construction snapshot accumulated across repeated real remaps");++snapshotPairs;
}
void inspect(const QuadtreeMesh& sourceMesh,const HydroStrainOperator& source,const QuadtreeMesh& targetMesh,
    const HydroStrainOperator& actual)
{
    if(!enabled)return;const auto& r=actual.geometryReuseReport();const auto& layout=actual.layout();
    std::vector<std::size_t> pointRows(layout.pressureCount());
    for(const auto& q:actual.points()){
        check(q.cellDof>=0&&q.cellDof<layout.pressureCount(),"raw/volume CSR actual quadrature has foreign cell");++pointRows[q.cellDof];}
    check(!pointRows.empty()&&pointRows.front()>0&&std::all_of(pointRows.begin(),pointRows.end(),
        [&](std::size_t n){return n==pointRows.front();}),"raw/volume CSR canonical point multiplicities differ");
    const auto components=static_cast<std::size_t>(actual.components());
    check(r.reusedRawCells>0&&r.reusedRawCells<=pointRows.size(),"raw/volume CSR scene has no certified copied raw cells");
    const auto copiedPoints=pointRows.front()*r.reusedRawCells;
    const auto potentialRaw=copiedPoints*components,potentialVolume=copiedPoints*2;
    // Include empty scalar rows in these shape-derived upper bounds; do not
    // equate net leaf difference, CSR nnz or unreported private eligibility.
    check(actual.canonicalRawStrain().matrix().rows()==static_cast<Eigen::Index>(actual.points().size()*components)
        &&actual.stressVolumeVelocityMatrix().rows()==static_cast<Eigen::Index>(actual.points().size()*2),
        "raw/volume CSR observed scalar-row shape is not actual quadrature based");
    check(r.rawDirectReboundRows<=potentialRaw&&r.volumeDirectReboundRows<=potentialVolume,
        "raw/volume direct rows exceed actual copied-cell scalar-row potential");
    check(r.csrRebindScratchBytes<=64u*1024u*1024u,"raw/volume shared representation scratch exceeds 64MiB");
    // Independently match exact physical face endpoints to fresh target DOFs;
    // the producer's column map, normalized cache key and row selector are not read.
    std::map<CertifiedCsrRebind::FaceKey,int> targetByFace;
    for(int u=0;u<layout.velocityCount();++u)
        check(targetByFace.emplace(CertifiedCsrRebind::faceKey(targetMesh,layout,u),u).second,
            "raw/volume CSR public target face geometry is not unique");
    std::vector<int> toTarget(source.layout().velocityCount(),-1),mapped;
    for(int u=0;u<source.layout().velocityCount();++u){
        const auto found=targetByFace.find(CertifiedCsrRebind::faceKey(sourceMesh,source.layout(),u));
        if(found!=targetByFace.end()){toTarget[u]=found->second;mapped.push_back(found->second);}}
    std::sort(mapped.begin(),mapped.end());
    check(std::adjacent_find(mapped.begin(),mapped.end())==mapped.end(),"raw/volume public correspondence is not injective");
    const auto raw=source.canonicalRawStrain();
    const Matrix* matrices[]{&raw.matrix(),&source.stressVolumeVelocityMatrix()};
    std::array<std::size_t,2> possible{},missing{},unordered{},empty{};
    for(unsigned kind=0;kind<2;++kind)for(Eigen::Index row=0;row<matrices[kind]->rows();++row){
        bool bound=true,increasing=true,normal=true,any=false;int previous=-1;
        for(Matrix::InnerIterator e(*matrices[kind],row);e;++e){
            any=true;check(e.col()>=0&&static_cast<std::size_t>(e.col())<toTarget.size(),"raw/volume source CSR column is out of range");
            const int column=toTarget[static_cast<std::size_t>(e.col())];
            if(column<0)bound=false;if(column<=previous)increasing=false;previous=column;
            const auto magnitude=bits(e.value())&UINT64_C(0x7fffffffffffffff),exponent=magnitude&UINT64_C(0x7ff0000000000000);
            if((exponent==0&&magnitude!=0)||exponent==UINT64_C(0x7ff0000000000000))normal=false;
        }
        if(!any)++empty[kind];
        if(!bound){++missing[kind];continue;}
        if(!increasing){++unordered[kind];continue;}
        if(normal)++possible[kind];
    }
    // This bounds all public source rows, not just the unknown private subset.
    // Existing complete cold CSR + polynomial/moment/work assertions below
    // establish numerical correctness independently of these metadata checks.
    check(r.rawDirectReboundRows<=possible[0]&&r.volumeDirectReboundRows<=possible[1],
        "raw/volume direct report exceeds every publicly eligible ordered normal source row");
#if defined(_MSC_VER) && defined(_M_X64) && !defined(__FAST_MATH__) && (!defined(_M_FP_FAST) || !_M_FP_FAST) && EIGEN_WORLD_VERSION==3 && EIGEN_MAJOR_VERSION==4 && EIGEN_MINOR_VERSION==1 && !defined(EIGEN_SPARSEMATRIX_PLUGIN) && !defined(EIGEN_SPARSE_COMPRESSEDBASE_PLUGIN)
    check(sizeof(long double)==sizeof(double)&&r.rawDirectReboundRows>0&&r.volumeDirectReboundRows>0,
        "raw/volume supported real AMR scene exercised no direct raw/value rows");
#else
    check(r.rawDirectReboundRows==0&&r.volumeDirectReboundRows==0,
        "raw/volume unsupported implementation used a certified CSR representation");
#endif
    rawRows+=r.rawDirectReboundRows;volumeRows+=r.volumeDirectReboundRows;
    potentialRawRows+=potentialRaw;potentialVolumeRows+=potentialVolume;orderedRawRows+=possible[0];orderedVolumeRows+=possible[1];
    for(unsigned k=0;k<2;++k){publicMissingRows[k]+=missing[k];publicUnorderedRows[k]+=unordered[k];emptySourceRows[k]+=empty[k];}
    ++scenes;
}
void finish()
{
    check(scenes==6&&fpPairs==6&&snapshotPairs==6&&coldFallbacks>=9&&potentialRawRows>0&&potentialVolumeRows>0,
        "raw/volume CSR six actual directional scene/full-FP/snapshot/fallback coverage incomplete");
    std::cout<<"STRESS_RAW_VOLUME_CSR_REBIND_PASS,actual_AMR_scenes,"<<scenes<<",actual_raw_direct_rows,"<<rawRows
        <<",actual_volume_direct_rows,"<<volumeRows<<",point_component_potential_raw_rows,"<<potentialRawRows
        <<",point_component_potential_volume_rows,"<<potentialVolumeRows<<",public_ordered_raw_source_rows,"<<orderedRawRows
        <<",public_ordered_volume_source_rows,"<<orderedVolumeRows<<",public_unbound_raw_rows,"<<publicMissingRows[0]
        <<",public_unbound_volume_rows,"<<publicMissingRows[1]<<",public_nonascending_raw_rows,"<<publicUnorderedRows[0]
        <<",public_nonascending_volume_rows,"<<publicUnorderedRows[1]<<",observed_empty_raw_source_rows,"<<emptySourceRows[0]
        <<",observed_empty_volume_source_rows,"<<emptySourceRows[1]<<",same_entry_cold_reuse_full_FP_errno_pairs,"<<fpPairs
        <<",repeat_construction_snapshot_pairs,"<<snapshotPairs<<",cold_unsupported_snapshots,"<<coldFallbacks
        <<",original_41_three_product_counters_and_all_cold_CSR_P2_moment_axis_virtualwork_source_release_gates_retained"
        <<",private_nonascending_missing_duplicate_cap_fault_injection_unavailable,no_speed_claim\n";
}
}

// A failed original directed-rounding factory remains a real failure. The
// guard only restores the caller's complete floating environment on exit.
struct FloatingEnvironmentRestore {
    std::fenv_t saved{};
#if defined(_MSC_VER) && defined(_M_X64)
    unsigned csr=0;
#endif
    bool active=true;
    FloatingEnvironmentRestore() {
        check(std::fegetenv(&saved)==0,"AMR source FP fixture cannot capture caller environment");
#if defined(_MSC_VER) && defined(_M_X64)
        csr=_mm_getcsr();
#endif
    }
    ~FloatingEnvironmentRestore() noexcept {if(active){std::fesetenv(&saved);
#if defined(_MSC_VER) && defined(_M_X64)
        _mm_setcsr(csr);
#endif
    }}
    void restore() {
        check(std::fesetenv(&saved)==0,"AMR source FP fixture cannot restore caller environment");
#if defined(_MSC_VER) && defined(_M_X64)
        _mm_setcsr(csr);
#endif
        active=false;
    }
};
void sourceFpFallback(const QuadtreeMesh& sourceMesh,const HydroPhysicalDivergence& sourceMetrics,
    const QuadtreeMesh& targetMesh,const HydroPhysicalDivergence& targetMetrics,const NativeLeafOverlap& overlap) {
    FloatingEnvironmentRestore restore;
    check(std::fesetround(FE_UPWARD)==0&&std::fegetround()==FE_UPWARD,"AMR source FP fixture cannot establish directed rounding");
#if defined(_MSC_VER) && defined(_M_X64)
    check((_mm_getcsr()&0x6000u)==0x4000u,"AMR source FP fixture did not establish actual SSE upward rounding");
#endif
    // Deliberately no exception-catching/skip: original factory failure must
    // escape the selector, rather than masquerade as a reuse-guard pass.
    const auto upwardSource=HydroStrainOperator::buildStressCompatible(sourceMesh,sourceMetrics);
    check(std::fesetround(FE_TONEAREST)==0&&std::fegetround()==FE_TONEAREST,"AMR source FP fixture cannot restore target nearest rounding");
#if defined(_MSC_VER) && defined(_M_X64)
    check((_mm_getcsr()&0x6000u)==0,"AMR source FP fixture did not establish actual SSE nearest rounding");
#endif
    const auto nearestCold=HydroStrainOperator::buildStressCompatible(targetMesh,targetMetrics);
    const auto fallback=HydroStrainOperator::buildStressCompatibleRemapped(sourceMesh,upwardSource,targetMesh,targetMetrics,overlap);
    check(fallback.geometryReuseReport().requested&&!fallback.geometryReuseReport().enabled,
        "nearest target reused coefficients built under an ineligible source FP environment");
    CertifiedCsrRebind::absent(fallback);
    same(fallback,nearestCold);restore.restore();
    std::cout<<"STRESS_AMR_SOURCE_FP_FALLBACK,actual_upward_source,actual_nearest_target,literal_cold_bits\n";
}
void scene(int metric,bool coarsen) {
    const bool rz=metric!=0;auto sourceMesh=std::make_unique<QuadtreeMesh>(input(16,metric));
    DictionaryMesh::RefineBox base;base.min=sourceMesh->minCoord+(sourceMesh->maxCoord-sourceMesh->minCoord)*.25;
    base.max=sourceMesh->minCoord+(sourceMesh->maxCoord-sourceMesh->minCoord)*.375;base.level=1;sourceMesh->applyRefineBoxes({base});
    DictionaryMesh::RefineBox added;const auto origin=sourceMesh->minCoord,w=sourceMesh->maxCoord-sourceMesh->minCoord;
    added.min={origin.x+(metric==1?0.:.625)*w.x,origin.y+.625*w.y};
    added.max={origin.x+(metric==1?.125:.75)*w.x,origin.y+.75*w.y};added.level=1;
    if(coarsen)sourceMesh->applyRefineBoxes({added});
    auto targetMesh=sourceMesh->clone();if(coarsen)coarsenAddedRegion(*targetMesh);else targetMesh->applyRefineBoxes({added});
    check(sourceMesh->cellIndexNum()!=targetMesh->cellIndexNum(),"AMR stress fixture has no actual topology change");
    const auto sourceMetrics=metrics(*sourceMesh,rz),targetMetrics=metrics(*targetMesh,rz);
    check(!sourceMetrics.layout().domain().uniformLeaves()&&!targetMetrics.layout().domain().uniformLeaves(),"AMR stress source/target lost true nonuniform 2:1 leaves");
    const auto overlap=NativeLeafOverlap::build(*sourceMesh,sourceMetrics.layout().domain(),*targetMesh,targetMetrics.layout().domain());
    std::size_t changed=0,unchanged=0;for(const auto& p:overlap.pieces()){
        if(p.kind==NativeLeafOverlapPiece::Kind::Unchanged)++unchanged;else {++changed;
            check(p.kind==(coarsen?NativeLeafOverlapPiece::Kind::Coarsened:NativeLeafOverlapPiece::Kind::Refined),"AMR stress fixture changed in wrong actual direction");}}
    check(changed>0&&unchanged>0,"AMR stress overlap lacks changed and unchanged witnesses");++changedPairs;
    std::optional<HydroStrainOperator> sourceGeometry(HydroStrainOperator::buildStressCompatible(*sourceMesh,sourceMetrics));
    const Matrix sourceE=sourceGeometry->matrix(),sourceRaw=sourceGeometry->canonicalRawStrain().matrix();
    std::optional<RawVolumeCsrRebind::ColdPair> rawVolumePair;
    if(RawVolumeCsrRebind::enabled)rawVolumePair.emplace(*sourceGeometry);
    const auto cold=HydroStrainOperator::buildStressCompatible(*targetMesh,targetMetrics);
    if(rawVolumePair)rawVolumePair->beforeRemap();
    const auto reused=HydroStrainOperator::buildStressCompatibleRemapped(*sourceMesh,*sourceGeometry,*targetMesh,targetMetrics,overlap);
    if(rawVolumePair)rawVolumePair->afterRemap(*sourceGeometry);
    same(reused,cold);const auto& rr=reused.geometryReuseReport();
    CertifiedCsrRebind::inspect(*sourceMesh,*sourceGeometry,*targetMesh,targetMetrics,overlap,cold,reused);
    RawVolumeCsrRebind::inspect(*sourceMesh,*sourceGeometry,*targetMesh,reused);
    check(rr.requested&&rr.enabled&&rr.targetCells==static_cast<std::size_t>(targetMetrics.layout().pressureCount()),"AMR reuse actual target report was not enabled");
    check(rr.reusedRawCells>0&&rr.reusedRawCells<=rr.targetCells&&rr.reusedStressCells>0&&rr.reusedStressCells<=rr.targetCells
        &&rr.reusedFaces>0&&rr.rebuiltFaces>0&&rr.recomputedFitCells>0,"AMR reuse did not exercise clean blocks and changed support work");
    check(!reused.layout().sameIdentity(sourceGeometry->layout())&& &reused.metrics().divergence().matrix()!=&sourceGeometry->metrics().divergence().matrix(),"AMR reuse retained source layout/D owner");
    matrix(sourceGeometry->matrix(),sourceE,"AMR reuse mutated source Ehat");matrix(sourceGeometry->canonicalRawStrain().matrix(),sourceRaw,"AMR reuse mutated source raw CSR");
    const auto sourceGreen=HydroStrainOperator::buildRt0Green(*sourceMesh,sourceMetrics);
    const auto coldGreen=HydroStrainOperator::buildRt0Green(*targetMesh,targetMetrics);
    const auto reusedGreen=HydroStrainOperator::buildRt0GreenRemapped(*sourceMesh,sourceGreen,*targetMesh,targetMetrics,overlap);
    Rt0GreenContract::same(reusedGreen,cold);Rt0GreenContract::same(reusedGreen,coldGreen);Rt0GreenContract::reused(reusedGreen);
    if(metric==0&&!coarsen){
        const auto fullFromGreen=HydroStrainOperator::buildStressCompatibleRemapped(*sourceMesh,sourceGreen,*targetMesh,targetMetrics,overlap);
        check(fullFromGreen.isStressCompatible()&&fullFromGreen.hasCanonicalRawStrain()&&fullFromGreen.geometryReuseReport().requested&&!fullFromGreen.geometryReuseReport().enabled,"Full request relabelled Green source or invented Full reuse");
        Rt0GreenContract::same(fullFromGreen,cold);Rt0GreenContract::matrix(fullFromGreen.matrix(),cold.matrix());
        Rt0GreenContract::matrix(fullFromGreen.canonicalRawStrain().matrix(),cold.canonicalRawStrain().matrix());
        const auto reverseGreen=NativeLeafOverlap::build(*targetMesh,targetMetrics.layout().domain(),*sourceMesh,sourceMetrics.layout().domain());
        Rt0GreenContract::missing([&]{(void)HydroStrainOperator::buildRt0GreenRemapped(*sourceMesh,sourceGreen,*targetMesh,targetMetrics,reverseGreen);});
        Rt0GreenContract::missing([&]{(void)HydroStrainOperator::buildRt0GreenRemapped(*sourceMesh,sourceGreen,*targetMesh,sourceMetrics,overlap);});
    }
    if(metric==0&&!coarsen){
        auto unrelated=std::make_unique<QuadtreeMesh>(input(8,0));const auto unrelatedMetrics=metrics(*unrelated,false);
        const auto unrelatedGeometry=HydroStrainOperator::buildStressCompatible(*unrelated,unrelatedMetrics);
        Rt0GreenContract::missing([&]{(void)HydroStrainOperator::buildRt0GreenRemapped(*unrelated,sourceGreen,*targetMesh,targetMetrics,overlap);});
        Rt0GreenContract::missing([&]{(void)HydroStrainOperator::buildRt0GreenRemapped(*sourceMesh,unrelatedGeometry,*targetMesh,targetMetrics,overlap);});
        Rt0GreenContract::missing([&]{(void)HydroStrainOperator::buildRt0GreenRemapped(*sourceMesh,sourceGreen,*unrelated,unrelatedMetrics,overlap);});
        rejects([&]{HydroStrainOperator::buildStressCompatibleRemapped(*unrelated,*sourceGeometry,*targetMesh,targetMetrics,overlap);});
        rejects([&]{HydroStrainOperator::buildStressCompatibleRemapped(*sourceMesh,unrelatedGeometry,*targetMesh,targetMetrics,overlap);});
        rejects([&]{HydroStrainOperator::buildStressCompatibleRemapped(*sourceMesh,*sourceGeometry,*unrelated,unrelatedMetrics,overlap);});
        const auto reverse=NativeLeafOverlap::build(*targetMesh,targetMetrics.layout().domain(),*sourceMesh,sourceMetrics.layout().domain());
        rejects([&]{HydroStrainOperator::buildStressCompatibleRemapped(*sourceMesh,*sourceGeometry,*targetMesh,targetMetrics,reverse);});
        rejects([&]{HydroStrainOperator::buildStressCompatibleRemapped(*sourceMesh,*sourceGeometry,*targetMesh,sourceMetrics,overlap);});
        sourceFpFallback(*sourceMesh,sourceMetrics,*targetMesh,targetMetrics,overlap);
        const auto raw=sourceGeometry->canonicalRawStrain();
        const auto rawFallback=HydroStrainOperator::buildStressCompatibleRemapped(*sourceMesh,raw,*targetMesh,targetMetrics,overlap);
        check(rawFallback.geometryReuseReport().requested&&!rawFallback.geometryReuseReport().enabled,"noncanonical raw source pretended to reuse stress blocks");same(rawFallback,cold);CertifiedCsrRebind::absent(rawFallback);
        const auto partitioned=HydroStrainOperator::buildStressCompatible(*sourceMesh,sourceMetrics,sourceGeometry->points(),sourceGeometry->stressFacePoints());
        check(!partitioned.hasCanonicalRawStrain(),"partitioned source acquired canonical provenance");
        const auto partitionFallback=HydroStrainOperator::buildStressCompatibleRemapped(*sourceMesh,partitioned,*targetMesh,targetMetrics,overlap);
        check(partitionFallback.geometryReuseReport().requested&&!partitionFallback.geometryReuseReport().enabled,"partitioned source pretended to reuse canonical blocks");same(partitionFallback,cold);CertifiedCsrRebind::absent(partitionFallback);
    }
    sourceGeometry.reset();sourceMesh.reset(); // No caller-owned source mesh/geometry remains for the numerical oracle below.
    FaceSegmentGraph targetGraph;targetGraph.build(*targetMesh);physical(*targetMesh,targetGraph,reused);Rt0GreenContract::physical(*targetMesh,targetGraph,reusedGreen);
    reusedRaw+=rr.reusedRawCells;reusedStress+=rr.reusedStressCells;reusedFaces+=rr.reusedFaces;fitWork+=rr.recomputedFitCells;
    std::cout<<"STRESS_AMR_REUSE,"<<metric<<','<<(coarsen?"coarsen":"refine")<<",target_cells,"<<rr.targetCells
        <<",changed_overlap_pieces,"<<changed<<",reused_raw,"<<rr.reusedRawCells<<",reused_stress,"<<rr.reusedStressCells
        <<",reused_faces,"<<rr.reusedFaces<<",rebuilt_faces,"<<rr.rebuiltFaces<<",actual_fit_cells,"<<rr.recomputedFitCells<<'\n';
}
void run() {
    for(int metric:{0,1,2})for(bool coarsen:{false,true})scene(metric,coarsen);
    check(changedPairs==6,"AMR stress directional scene coverage incomplete");
    std::cout<<"HYDRO_STRESS_AMR_LOCAL_REUSE_PASS,scenes,6,scalar_bits,"<<scalarChecks<<",matrix_entries,"<<matrixEntries
        <<",reused_raw,"<<reusedRaw<<",reused_stress,"<<reusedStress<<",reused_faces,"<<reusedFaces<<",fit_work,"<<fitWork
        <<",guards,"<<rejected<<",independent_P2_moments_constant_tensor_same_C_virtualwork_original_1e-12,no_speed_claim\n";
}
} // namespace StressAmrLocalReuseContract

namespace SelectedStressAmrReuseContract
{
namespace Original=StressAmrLocalReuseContract;
using Matrix=HydroStrainOperator::Matrix;
using LAMP::Mesh::NativeLeafOverlap;
using LAMP::Mesh::NativeLeafOverlapPiece;
using Floating=Original::CertifiedCsrRebind::Floating;
std::size_t scenes=0,guardPairs=0,fpPairs=0,internalFaces=0,excludedFaces=0,reusedRaw=0,reusedStress=0,fitWork=0;
bool diagnosticCoarsen=false;
bool selected(Vector2D p,bool alternative)
{
    // Keep the full RZ axis in both equal-count, different selections: only
    // reflect z. A central solid hole has six coarse cells to each outer wall;
    // the nozzle notch is separate. Both actual refine boxes only touch the
    // hole at a corner, so every refined leaf remains selected in A and B.
    if(alternative)p.y=1-p.y;
    return !(p.x>.75&&p.y<.125)
        &&!(p.x>.375&&p.x<.625&&p.y>.375&&p.y<.625);
}
// Test diagnostics only. Every call returns BEFORE the corresponding complete
// entry.restore(); stream formatting/allocation/errno cannot enter the oracle.
void diagnosticStage(const char* stage,const QuadtreeMesh& mesh,const HydroPhysicalDivergence& metrics)
{
    const auto& layout=metrics.layout();
    std::cerr<<std::setprecision(17)<<"F16_STAGE,"<<(layout.domain().axisymmetric()?"RZ":"XY")
        <<','<<(diagnosticCoarsen?"coarsen":"refine")<<','<<stage
        <<",selected_cells,"<<layout.pressureCount()<<",native_cells,"<<mesh.cellIndexNum();
    if(layout.pressureCount()>37){const int native=layout.nativeCell(37);
        const auto& cell=mesh.cells(mesh.cellIndexToLabel[native]);const auto center=cell.centerCoord(),width=cell.size();
        std::cerr<<",selected_p,37,native,"<<native<<",center,"<<center.x<<','<<center.y
            <<",width,"<<width.x<<','<<width.y<<",A_selected,"<<selected(center,false)<<",B_selected,"<<selected(center,true);}
    std::cerr<<'\n'<<std::flush;
}
FaceSegmentGraph graph(const QuadtreeMesh& mesh,bool alternative)
{
    FaceSegmentGraph result;
    result.build(mesh,[alternative](const QuadtreeMesh&,const Cell& cell){return selected(cell.centerCoord(),alternative);});
    return result;
}
HydroPhysicalDivergence metrics(const QuadtreeMesh& mesh,const FaceSegmentGraph& graph,bool rz)
{
    const auto layout=HydroDofLayout::build(mesh,NativeLeafDomain::build(mesh,graph,rz));
    check(layout.pressureCount()<mesh.cellIndexNum(),"selected stress fixture retained the full native domain");
    return HydroPhysicalDivergence::build(mesh,HydroIntegralDivergence::build(mesh,layout));
}
void sameSelectedFp(const Floating& a,const Floating& b)
{
    check(std::memcmp(&a.native,&b.native,sizeof a.native)==0&&a.error==b.error,
        "selected stress cold/reuse changed complete native FP or errno");
#if defined(_MSC_VER) && defined(_M_X64)
    check(a.csr==b.csr,"selected stress cold/reuse changed full MXCSR controls/status");
#endif
    ++fpPairs;
}
std::array<Matrix,5> matrices(const HydroStrainOperator& source)
{
    return {source.matrix(),source.canonicalRawStrain().matrix(),source.stressVolumeVelocityMatrix(),
        source.stressFaceVelocityMatrix(),source.stressComponentIntegralMatrix()};
}
void selectionWitness(const HydroPhysicalDivergence& a,const HydroPhysicalDivergence& b)
{
    check(a.layout().pressureCount()==b.layout().pressureCount(),"wrong-selection stress witness lacks equal selected counts");
    check(a.layout().domain().cellRowByIndex()!=b.layout().domain().cellRowByIndex(),
        "wrong-selection stress witness has identical native membership");
}
void shapeWitness(const QuadtreeMesh& mesh,const HydroPhysicalDivergence& metrics)
{
    const auto& layout=metrics.layout();std::size_t internal=0,excluded=0,axis=0;
    for(int native=0;native<layout.domain().nativeFaceCount();++native){
        const auto f=layout.domain().face(mesh,native);const int u=layout.velocityDof(native);
        if(f.negativeCell<0&&f.positiveCell<0){check(u<0,"selected stress solid-only face acquired a velocity column");++excluded;continue;}
        check(u>=0&&layout.nativeFace(u)==native,"selected stress native velocity round trip failed");
        if(f.outwardSign!=0&&f.center.x>mesh.minCoord.x&&f.center.x<mesh.maxCoord.x
            &&f.center.y>mesh.minCoord.y&&f.center.y<mesh.maxCoord.y){
            check((f.negativeCell>=0)!=(f.positiveCell>=0),"selected stress CSG wall lost one-sided incidence");++internal;}
        if(f.area==0){check(layout.domain().axisymmetric()&&!f.horizontal&&f.center.x==0,
            "selected stress invented a zero-area nonaxis face");++axis;}
    }
    check(internal>0&&excluded>0,"selected stress lacks actual internal solid walls or excluded faces");
    if(layout.domain().axisymmetric())check(axis>0,"selected stress RZ scene lost the real axis");
    internalFaces+=internal;excludedFaces+=excluded;
}
void wrongSelection(const QuadtreeMesh& sourceMesh,const HydroStrainOperator& source,
    const QuadtreeMesh& targetMesh,const HydroPhysicalDivergence& targetB,
    const NativeLeafOverlap& overlapB,const NativeLeafOverlap& overlapA,const Floating& entry)
{
    // A(source operator)/B(overlap source), B(target)/B(overlap target): only
    // the SOURCE consumer membership is wrong, although all meshes are live.
    diagnosticStage("cold_B",targetMesh,targetB);entry.restore();const auto coldB=HydroStrainOperator::buildStressCompatible(targetMesh,targetB);const Floating coldExit;
    diagnosticStage("source_only_wrong_selection",targetMesh,targetB);entry.restore();const auto sourceOnly=HydroStrainOperator::buildStressCompatibleRemapped(sourceMesh,source,targetMesh,targetB,overlapB);
    const Floating sourceExit;sameSelectedFp(sourceExit,coldExit);Original::same(sourceOnly,coldB);
    check(sourceOnly.geometryReuseReport().requested&&!sourceOnly.geometryReuseReport().enabled,
        "selected stress accepted a same-mesh source-only wrong selection");
    Original::CertifiedCsrRebind::absent(sourceOnly);Original::RawVolumeCsrRebind::absent(sourceOnly);++guardPairs;
    // A(source)/A(overlap source), B(target)/A(overlap target): only TARGET is
    // wrong. The reverse overlap and stale-mesh guards remain in old fixtures.
    diagnosticStage("target_only_wrong_selection",targetMesh,targetB);entry.restore();const auto targetOnly=HydroStrainOperator::buildStressCompatibleRemapped(sourceMesh,source,targetMesh,targetB,overlapA);
    const Floating targetExit;sameSelectedFp(targetExit,coldExit);Original::same(targetOnly,coldB);
    check(targetOnly.geometryReuseReport().requested&&!targetOnly.geometryReuseReport().enabled,
        "selected stress accepted a same-mesh target-only wrong selection");
    Original::CertifiedCsrRebind::absent(targetOnly);Original::RawVolumeCsrRebind::absent(targetOnly);++guardPairs;
}
void scene(bool rz,bool coarsen)
{
    diagnosticCoarsen=coarsen;
    auto sourceMesh=std::make_unique<QuadtreeMesh>(input(16,rz?1:0));
    DictionaryMesh::RefineBox base;base.min={.25,.25};base.max={.375,.375};base.level=1;
    sourceMesh->applyRefineBoxes({base});
    DictionaryMesh::RefineBox added;added.min={.625,.625};added.max={.75,.75};added.level=1;
    if(coarsen)sourceMesh->applyRefineBoxes({added});
    auto targetMesh=sourceMesh->clone();if(coarsen)Original::coarsenAddedRegion(*targetMesh);else targetMesh->applyRefineBoxes({added});
    check(sourceMesh->cellIndexNum()!=targetMesh->cellIndexNum(),"selected stress skipped actual tree change");
    const auto sourceGraphA=graph(*sourceMesh,false),sourceGraphB=graph(*sourceMesh,true);
    const auto targetGraphA=graph(*targetMesh,false),targetGraphB=graph(*targetMesh,true);
    const auto sourceA=metrics(*sourceMesh,sourceGraphA,rz),sourceB=metrics(*sourceMesh,sourceGraphB,rz);
    const auto targetA=metrics(*targetMesh,targetGraphA,rz),targetB=metrics(*targetMesh,targetGraphB,rz);
    selectionWitness(sourceA,sourceB);selectionWitness(targetA,targetB);
    check(!sourceA.layout().domain().uniformLeaves()&&!targetA.layout().domain().uniformLeaves(),
        "selected stress source/target lack actual 2:1 AMR");
    shapeWitness(*sourceMesh,sourceA);shapeWitness(*targetMesh,targetA);
    const auto overlapA=NativeLeafOverlap::build(*sourceMesh,sourceA.layout().domain(),*targetMesh,targetA.layout().domain());
    const auto overlapB=NativeLeafOverlap::build(*sourceMesh,sourceB.layout().domain(),*targetMesh,targetB.layout().domain());
    std::size_t changed=0,unchanged=0;
    for(const auto& piece:overlapA.pieces())if(piece.kind==NativeLeafOverlapPiece::Kind::Unchanged)++unchanged;
        else{++changed;check(piece.kind==(coarsen?NativeLeafOverlapPiece::Kind::Coarsened:NativeLeafOverlapPiece::Kind::Refined),
            "selected stress tree changed in the wrong direction");}
    check(changed>0&&unchanged>0,"selected stress has no changed and unchanged physical overlap");
    const Floating sourceCallEntry;diagnosticStage("source_A",*sourceMesh,sourceA);sourceCallEntry.restore();
    std::optional<HydroStrainOperator> source(HydroStrainOperator::buildStressCompatible(*sourceMesh,sourceA));const Floating entry;
    const auto sourceBefore=matrices(*source);
    diagnosticStage("cold_A",*targetMesh,targetA);entry.restore();const auto cold=HydroStrainOperator::buildStressCompatible(*targetMesh,targetA);const Floating coldExit;
    diagnosticStage("reuse_A",*targetMesh,targetA);entry.restore();const auto reused=HydroStrainOperator::buildStressCompatibleRemapped(*sourceMesh,*source,*targetMesh,targetA,overlapA);
    const Floating reuseExit;sameSelectedFp(reuseExit,coldExit);Original::same(reused,cold);
    const auto report=reused.geometryReuseReport();
    check(report.requested&&report.enabled&&report.targetCells==static_cast<std::size_t>(targetA.layout().pressureCount()),
        "selected stress actual remap was not enabled on the new selected layout");
    check(report.reusedRawCells>0&&report.reusedStressCells>0&&report.reusedFaces>0&&report.rebuiltFaces>0&&report.recomputedFitCells>0,
        "selected stress did not exercise certified copies and changed/halo recomputation");
    check(!reused.layout().sameIdentity(source->layout())&& &reused.metrics().divergence().matrix()!=&source->metrics().divergence().matrix(),
        "selected stress relabelled the old native layout or D");
    const Floating inspectCallEntry;diagnosticStage("repeat_reuse_CSR_oracle",*targetMesh,targetA);inspectCallEntry.restore();
    Original::CertifiedCsrRebind::inspect(*sourceMesh,*source,*targetMesh,targetA,overlapA,cold,reused);
    Original::RawVolumeCsrRebind::inspect(*sourceMesh,*source,*targetMesh,reused);
    wrongSelection(*sourceMesh,*source,*targetMesh,targetB,overlapB,overlapA,entry);
    const auto sourceGreen=HydroStrainOperator::buildRt0Green(*sourceMesh,sourceA);
    const auto greenCold=HydroStrainOperator::buildRt0Green(*targetMesh,targetA);
    const auto greenReused=HydroStrainOperator::buildRt0GreenRemapped(*sourceMesh,sourceGreen,*targetMesh,targetA,overlapA);
    Rt0GreenContract::same(greenReused,cold);Rt0GreenContract::same(greenReused,greenCold);Rt0GreenContract::reused(greenReused);
    const auto greenColdB=HydroStrainOperator::buildRt0Green(*targetMesh,targetB);
    const auto wrongSource=HydroStrainOperator::buildRt0GreenRemapped(*sourceMesh,sourceGreen,*targetMesh,targetB,overlapB);
    const auto wrongTarget=HydroStrainOperator::buildRt0GreenRemapped(*sourceMesh,sourceGreen,*targetMesh,targetB,overlapA);
    for(const auto* g:{&wrongSource,&wrongTarget}){Rt0GreenContract::same(*g,greenColdB);check(g->geometryReuseReport().requested&&!g->geometryReuseReport().enabled,"Green remap accepted one-sided wrong selected membership");}
    const auto sourceAfter=matrices(*source);
    for(std::size_t k=0;k<sourceBefore.size();++k)Original::matrix(sourceAfter[k],sourceBefore[k],"selected stress remap mutated an actual source operator");
    source.reset();sourceMesh.reset();
    // The independent oracle must receive the SELECTED graph. A full graph
    // would incorrectly feed excluded native faces (velocityDof==-1) to it.
    const Floating physicalCallEntry;diagnosticStage("physical_oracle_including_cold_raw",*targetMesh,targetA);physicalCallEntry.restore();
    Original::physical(*targetMesh,targetGraphA,reused);Rt0GreenContract::physical(*targetMesh,targetGraphA,greenReused);
    reusedRaw+=report.reusedRawCells;reusedStress+=report.reusedStressCells;fitWork+=report.recomputedFitCells;++scenes;
    std::cout<<"SELECTED_STRESS_AMR_REUSE,"<<(rz?"RZ":"XY")<<','<<(coarsen?"coarsen":"refine")
        <<",selected_cells,"<<report.targetCells<<",native_cells,"<<targetMesh->cellIndexNum()<<",changed_pieces,"<<changed
        <<",reused_raw,"<<report.reusedRawCells<<",reused_stress,"<<report.reusedStressCells<<",actual_fit_cells,"<<report.recomputedFitCells<<'\n';
}
void run()
{
    // Reuse the existing independent public physical-face/CSR oracles. Their
    // six-scene finish() belongs to the unchanged original selectors only.
    Original::CertifiedCsrRebind::enabled=true;Original::RawVolumeCsrRebind::enabled=true;
    for(bool rz:{false,true})for(bool coarsen:{false,true})scene(rz,coarsen);
    check(scenes==4&&guardPairs==8&&fpPairs==12&&internalFaces>0&&excludedFaces>0,
        "selected stress actual directional/one-sided-selection/full-FP coverage incomplete");
    std::cout<<"HYDRO_SELECTED_STRESS_AMR_REUSE_PASS,scenes,"<<scenes<<",one_sided_equal_count_wrong_selection_pairs,"<<guardPairs
        <<",same_entry_cold_full_FP_errno_pairs,"<<fpPairs<<",internal_wall_faces,"<<internalFaces<<",excluded_faces,"<<excludedFaces
        <<",reused_raw,"<<reusedRaw<<",reused_stress,"<<reusedStress<<",actual_fit_cells,"<<fitWork
        <<",original_P2_native_moments_same_C_boundary_virtualwork_source_release_gates,no_speed_claim\n";
}
} // namespace SelectedStressAmrReuseContract


namespace GramCoercivityContract {
namespace FV=LAMP::Geometry::Detail::FiniteVolume;
using Matrix=HydroEnergyViscosity::Matrix;
std::size_t exactRows=0,exactTerms=0,mathCases=0,positiveBlocks=0,guardCases=0;
std::uint64_t word(double value){std::uint64_t bits=0;std::memcpy(&bits,&value,sizeof bits);return bits;}
std::optional<std::vector<double>> query(const HydroEnergyViscosity& energy,const std::vector<double>& scale) {
    std::fenv_t before{},after{};check(std::fegetenv(&before)==0,"Gram oracle cannot capture FP environment");
#if defined(_MSC_VER) && defined(_M_X64)
    const unsigned csr=_mm_getcsr();
#endif
    const int error=errno;auto result=energy.scaledGramLowerErrorBounds(scale);const int afterError=errno;
    check(std::fegetenv(&after)==0&&std::memcmp(&before,&after,sizeof before)==0&&afterError==error,"Gram bound query changed caller FENV/errno");
#if defined(_MSC_VER) && defined(_M_X64)
    check(csr==_mm_getcsr(),"Gram bound query changed caller raw MXCSR");
#endif
    return result;
}
std::size_t coarseRows=0;
std::optional<double> coarseQuery(const HydroEnergyViscosity& energy,double maximumScale) {
    std::fenv_t before{},after{};check(std::fegetenv(&before)==0,"coarse Gram oracle cannot capture FP environment");
#if defined(_MSC_VER) && defined(_M_X64)
    const unsigned csr=_mm_getcsr();
#endif
    const int error=errno;auto result=energy.coarseScaledGramLowerErrorBound(maximumScale);const int afterError=errno;
    check(std::fegetenv(&after)==0&&std::memcmp(&before,&after,sizeof before)==0&&afterError==error,"coarse Gram query changed caller FENV/errno");
#if defined(_MSC_VER) && defined(_M_X64)
    check(csr==_mm_getcsr(),"coarse Gram query changed caller raw MXCSR");
#endif
    return result;
}
void exactLowerOracle(const HydroEnergyViscosity& energy,const std::vector<double>& scale,bool tiny) {
    const auto bound=query(energy,scale);const auto& e=energy.strain().matrix();const auto& k=energy.matrix();const auto& w=energy.rowWeights();
    const auto n=static_cast<std::size_t>(e.cols());check(bound&&bound->size()==n&&n<=100,"small actual Gram oracle lacks certified bounds");
    const auto coarse=coarseQuery(energy,*std::max_element(scale.begin(),scale.end()));
    check(coarse&&std::isfinite(*coarse)&&*coarse>=0,"small actual Gram oracle lacks finite coarse scalar bound");
    std::vector<std::vector<FV::Dyadic>> gram(n,std::vector<FV::Dyadic>(n));
    bool subnormalWeight=false,subnormalK=false;
    for(double x:w)subnormalWeight=subnormalWeight||std::fpclassify(x)==FP_SUBNORMAL;
    for(Eigen::Index z=0;z<k.nonZeros();++z)subnormalK=subnormalK||std::fpclassify(k.valuePtr()[z])==FP_SUBNORMAL;
    if(tiny)check(subnormalWeight&&subnormalK,"underflow oracle did not retain real subnormal W and K");
    // Independent exact algebra, not a transcription of gamma/eta or interval
    // majorants. Every actual represented W/E factor enters the complete Gram.
    for(Eigen::Index row=0;row<e.outerSize();++row)
        for(Matrix::InnerIterator a(e,row);a;++a)for(Matrix::InnerIterator b(e,row);b;++b){
            gram[static_cast<std::size_t>(a.col())][static_cast<std::size_t>(b.col())].add(
                FV::Dyadic::product({w[static_cast<std::size_t>(row)],a.value(),b.value()}));++exactTerms;}
    for(std::size_t i=0;i<n;++i){FV::Dyadic rowError;
        for(std::size_t j=0;j<n;++j){
            // Exactly the symmetric Lower completion used by the original LLT.
            FV::Dyadic difference(k.coeff(static_cast<Eigen::Index>(std::max(i,j)),static_cast<Eigen::Index>(std::min(i,j))));
            difference.add(gram[i][j],-1);rowError.add(difference.absolute().multiplied(scale[i]).multiplied(scale[j]));}
        check(rowError.compare(FV::Dyadic(*coarse))<=0,"coarse scalar is below independent exact scaled Lower Gram row error");++coarseRows;
        check(std::isfinite((*bound)[i])&&(*bound)[i]>=0&&rowError.compare(FV::Dyadic((*bound)[i]))<=0,
            "certified Lower Gram row bound is below independent exact Dyadic error");++exactRows;}
    ++mathCases;
}
template<class Sparse> void matrixBits(const Sparse& a,const Sparse& b) {
    check(a.rows()==b.rows()&&a.cols()==b.cols()&&a.nonZeros()==b.nonZeros()&&a.isCompressed()&&b.isCompressed(),"energy SPD path changed actual matrix shape");
    for(Eigen::Index i=0;i<=a.outerSize();++i)check(a.outerIndexPtr()[i]==b.outerIndexPtr()[i],"energy SPD path changed actual matrix outer pattern");
    for(Eigen::Index i=0;i<a.nonZeros();++i)check(a.innerIndexPtr()[i]==b.innerIndexPtr()[i]&&word(a.valuePtr()[i])==word(b.valuePtr()[i]),"energy SPD path changed actual coefficient bits");
}
HydroSaddleVelocityBlock::Matrix actualA(const HydroMomentumMass& mass,const HydroEnergyViscosity& energy,double dt) {
    HydroSaddleVelocityBlock::Matrix a=.5*energy.matrix();
    for(Eigen::Index u=0;u<a.rows();++u)a.coeffRef(u,u)+=mass.masses()[static_cast<std::size_t>(u)]/dt;
    a.makeCompressed();return a;
}
void positive(const QuadtreeMesh& mesh,const HydroPhysicalDivergence& metrics,const HydroMomentumMass& mass,
    const HydroEnergyViscosity& energy,bool guards,const std::vector<double>& mu) {
    const auto& layout=metrics.layout();const double dt=1.;const auto a=actualA(mass,energy,dt);
    const HydroSaddleVelocityBlock::Matrix transpose(a.transpose()),symmetric=.5*(a+transpose);
    Eigen::MatrixXd lower=symmetric;for(Eigen::Index i=0;i<lower.rows();++i)for(Eigen::Index j=i+1;j<lower.cols();++j)lower(i,j)=lower(j,i);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> spectrum(lower);
    check(spectrum.info()==Eigen::Success&&spectrum.eigenvalues().allFinite()&&spectrum.eigenvalues().minCoeff()>0,"independent actual Lower eigenvalue oracle is not SPD");
    bool notDominant=false;for(Eigen::Index i=0;i<lower.rows();++i){long double radius=0;
        for(Eigen::Index j=0;j<lower.cols();++j)if(i!=j)radius+=std::abs(static_cast<long double>(lower(i,j)));
        notDominant=notDominant||radius>static_cast<long double>(lower(i,i))*1.001L;}
    check(notDominant,"typed Gram fixture does not exercise beyond strict diagonal dominance");
    const auto old=HydroSaddleVelocityBlock::buildWithTransportSplit(layout,a);
    const auto certified=HydroSaddleVelocityBlock::buildWithEnergyTransportSplit(mesh,layout,a,mass,energy,dt);
    const auto& report=certified.spdReport();
    check(report.tier==HydroVelocitySpdTier::EnergyGram&&report.gramAttempted&&report.gramBoundsAvailable
        &&std::isfinite(report.gramBoundUpper)&&report.gramBoundUpper>=0&&report.scaledMassLower>report.gramBoundUpper&&report.lltSeconds==0,
        "actual non-DD positive block did not obtain a strict energy Gram certificate");
    matrixBits(certified.matrix(),old.matrix());++positiveBlocks;
    std::cout<<"GRAM_SPD_POSITIVE,"<<(layout.domain().axisymmetric()?"RZ":"XY")<<",lambda_min,"<<spectrum.eigenvalues().minCoeff()
        <<",scaled_mass_lower,"<<report.scaledMassLower<<",error_upper,"<<report.gramBoundUpper<<'\n';
    if(!guards)return;
    auto indefinite=a;indefinite.coeffRef(0,0)=-std::abs(indefinite.coeff(0,0))-1.;indefinite.makeCompressed();
    check(indefinite.coeff(0,0)<0,"indefinite fixture has no negative unit-vector quadratic form");
    rejects([&]{HydroSaddleVelocityBlock::buildWithEnergyTransportSplit(mesh,layout,indefinite,mass,energy,dt);});++guardCases;
    const auto otherLayout=HydroDofLayout::build(mesh,layout.domain());
    rejects([&]{HydroSaddleVelocityBlock::buildWithEnergyTransportSplit(mesh,otherLayout,a,mass,energy,dt);});++guardCases;
    const auto otherMetrics=HydroPhysicalDivergence::build(mesh,HydroIntegralDivergence::build(mesh,layout));
    const auto otherEnergy=HydroEnergyViscosity::build(HydroStrainOperator::build(mesh,otherMetrics,
        HydroStrainOperator::volumeQuadrature(mesh,otherMetrics)),mu);
    rejects([&]{HydroSaddleVelocityBlock::buildWithEnergyTransportSplit(mesh,layout,a,mass,otherEnergy,dt);});++guardCases;
    QuadtreeMesh foreign(input(2,0));rejects([&]{HydroSaddleVelocityBlock::buildWithEnergyTransportSplit(foreign,layout,a,mass,energy,dt);});++guardCases;
    rejects([&]{HydroSaddleVelocityBlock::buildWithEnergyTransportSplit(mesh,layout,a,mass,energy,0.);});++guardCases;
    // Original directed-rounding builds/LLT must genuinely succeed: no catch
    // or skip can turn an unsupported original operation into a guard pass.
    StressAmrLocalReuseContract::FloatingEnvironmentRestore restore;
    check(std::fesetround(FE_UPWARD)==0&&std::fegetround()==FE_UPWARD,"Gram source FP fixture cannot set actual upward mode");
    const auto upwardEnergy=HydroEnergyViscosity::build(energy.strain(),mu);
    check(std::fesetround(FE_TONEAREST)==0&&std::fegetround()==FE_TONEAREST,"Gram source FP fixture cannot restore actual nearest mode");
    check(!query(upwardEnergy,std::vector<double>(layout.velocityCount(),1)),"nearest query certified a K built under unsupported source rounding");++guardCases;
    check(!coarseQuery(upwardEnergy,1.),"coarse query certified K built under unsupported source rounding");++guardCases;
    const auto upwardA=actualA(mass,upwardEnergy,dt);
    const auto literalFallback=HydroSaddleVelocityBlock::buildWithEnergyTransportSplit(mesh,layout,upwardA,mass,upwardEnergy,dt);
    const auto literalOld=HydroSaddleVelocityBlock::buildWithTransportSplit(layout,upwardA);
    check(literalFallback.spdReport().tier==HydroVelocitySpdTier::LLT&&literalFallback.spdReport().gramAttempted
        &&!literalFallback.spdReport().gramBoundsAvailable,"unsupported source FP did not preserve original LLT fallback");
    matrixBits(literalFallback.matrix(),literalOld.matrix());restore.restore();++guardCases;
}
void scaleAndFpGuards(const HydroEnergyViscosity& energy,const std::vector<double>& scale) {
    auto invalid=scale;invalid.pop_back();check(!query(energy,invalid),"Gram invalid scale dimension accepted");++guardCases;
    for(double value:{0.,-1.,std::numeric_limits<double>::quiet_NaN(),std::numeric_limits<double>::infinity(),std::numeric_limits<double>::max()}){
        invalid=scale;invalid[0]=value;check(!query(energy,invalid),"Gram invalid/overflow scale certified");++guardCases;}
    for(double maximum:{0.,-1.,std::numeric_limits<double>::quiet_NaN(),std::numeric_limits<double>::infinity(),std::numeric_limits<double>::max()}){
        check(!coarseQuery(energy,maximum),"coarse Gram invalid/overflow maximum scale certified");++guardCases;}
    StressAmrLocalReuseContract::FloatingEnvironmentRestore restore;
    check(std::fesetround(FE_UPWARD)==0&&std::fegetround()==FE_UPWARD,"Gram query FP fixture cannot set upward mode");
    check(!query(energy,scale),"Gram query certified unsupported current rounding");++guardCases;
    check(!coarseQuery(energy,*std::max_element(scale.begin(),scale.end())),"coarse Gram query certified unsupported current rounding");++guardCases;
    check(std::fesetround(FE_TONEAREST)==0,"Gram query FP fixture cannot restore nearest mode");
#if defined(_MSC_VER) && defined(_M_X64)
    const unsigned current=_mm_getcsr();_mm_setcsr(current|0x8040u);
    check(!coarseQuery(energy,*std::max_element(scale.begin(),scale.end())),"coarse Gram query certified FTZ/DAZ controls");++guardCases;
    check(!query(energy,scale),"Gram query certified FTZ/DAZ current controls");++guardCases;_mm_setcsr(current);
#endif
    restore.restore();
}
void scene(bool rz) {
    QuadtreeMesh mesh(input(4,rz?1:0));FaceSegmentGraph graph;graph.build(mesh);
    const auto layout=HydroDofLayout::build(mesh,NativeLeafDomain::build(mesh,graph,rz));
    const auto metrics=HydroPhysicalDivergence::build(mesh,HydroIntegralDivergence::build(mesh,layout));
    const auto raw=HydroStrainOperator::build(mesh,metrics,HydroStrainOperator::volumeQuadrature(mesh,metrics));
    std::vector<double> scale(layout.velocityCount()),mu(raw.points().size());
    for(std::size_t i=0;i<scale.size();++i)scale[i]=std::ldexp(1.,static_cast<int>(i%7)-3);
    for(std::size_t q=0;q<mu.size();++q)mu[q]=q%7==0?0.:.25+.125*(q%4);
    const auto energy=HydroEnergyViscosity::build(raw,mu);exactLowerOracle(energy,scale,false);
    const auto tiny=HydroEnergyViscosity::build(raw,std::vector<double>(mu.size(),1e-319));exactLowerOracle(tiny,scale,true);
    const auto momentum=HydroMomentumGeometry::build(mesh,metrics);
    const LAMP::TimeLayer clock{0,0,0,LAMP::TimeLayer::Position::Initial};
    const auto pub=LAMP::Interface::GeometricVofPublication::prepareLocalPlic(mesh,layout.domain(),clock,std::vector<double>(mesh.cellIndexNum(),1.));
    const auto mass=HydroMomentumMass::buildCellAverageProjection(mesh,momentum,HydroPhaseGeometry::build(mesh,metrics,pub),1.,1.);
    positive(mesh,metrics,mass,energy,!rz,mu);
    if(!rz){scaleAndFpGuards(energy,scale);
        const auto stress=HydroStrainOperator::buildStressCompatible(mesh,metrics);
        exactLowerOracle(HydroEnergyViscosity::build(stress,std::vector<double>(stress.points().size(),.5)),scale,false);
        const auto zero=HydroEnergyViscosity::zero(mesh,metrics);const auto bound=query(zero,scale);
        check(bound&&bound->size()==scale.size()&&std::all_of(bound->begin(),bound->end(),[](double x){return x==0;}),"exact zero Gram form has a nonzero or unavailable bound");
        const auto coarse=coarseQuery(zero,*std::max_element(scale.begin(),scale.end()));
        check(coarse&&*coarse==0,"exact zero Gram form has nonzero/unavailable coarse scalar");
    }
}
void run(){scene(false);scene(true);
    check(mathCases==5&&positiveBlocks==2,"Gram small algebra/physical geometry coverage incomplete");
    check(coarseRows==exactRows,"coarse Gram oracle did not visit all existing exact rows");
    std::cout<<"GRAM_COARSE_BOUND_PASS,existing_math_cases,"<<mathCases<<",independent_exact_scaled_rows,"<<coarseRows<<",source_current_FP_and_invalid_scale,1\n";
    std::cout<<"GRAM_COERCIVITY_PASS,math_cases,"<<mathCases<<",exact_scaled_rows,"<<exactRows<<",exact_terms,"<<exactTerms
        <<",positive_XY_RZ,"<<positiveBlocks<<",guards,"<<guardCases<<",actual_subnormal_W_K,source_current_FP,strict_original_SPD_gates,no_speed_claim\n";
}
} // namespace GramCoercivityContract

namespace ViscosityProductPlanContract {
using Matrix=HydroEnergyViscosity::Matrix;
struct Environment {
    int flags=0,rounding=0;
    unsigned control=0,csr=0;
    bool eligible=false;
};
Environment environment()
{
    Environment r;r.flags=std::fetestexcept(FE_ALL_EXCEPT);r.rounding=std::fegetround();
#if defined(_MSC_VER) && defined(_M_X64)
    const bool valid=_controlfp_s(&r.control,0,0)==0;r.csr=_mm_getcsr();
#if !defined(_M_FP_FAST) && !defined(__FAST_MATH__) && EIGEN_WORLD_VERSION==3 && ((EIGEN_MAJOR_VERSION==3 && EIGEN_MINOR_VERSION==9) || (EIGEN_MAJOR_VERSION==4 && EIGEN_MINOR_VERSION==1))
    r.eligible=valid&&(r.control&_MCW_EM)==_MCW_EM&&(r.csr&0x1f80u)==0x1f80u;
#if defined(_EM_AMBIGUOUS)
    r.eligible=r.eligible&&!(r.control&_EM_AMBIGUOUS);
#endif
#endif
#endif
    static bool announced=false;
    if(!announced){
#if defined(_M_FP_FAST)
        const int fast=_M_FP_FAST;
#else
        const int fast=-1;
#endif
#if defined(_M_FP_PRECISE)
        const int precise=_M_FP_PRECISE;
#else
        const int precise=-1;
#endif
#if defined(_M_FP_STRICT)
        const int strict=_M_FP_STRICT;
#else
        const int strict=-1;
#endif
#if defined(__FAST_MATH__)
        const int fastMath=__FAST_MATH__;
#else
        const int fastMath=-1;
#endif
        std::cout<<"VISCOSITY_PRODUCT_COMPILE,eigen,"<<EIGEN_WORLD_VERSION<<'/'<<EIGEN_MAJOR_VERSION<<'/'<<EIGEN_MINOR_VERSION
            <<",M_FP_FAST,"<<fast<<",M_FP_PRECISE,"<<precise<<",M_FP_STRICT,"<<strict
            <<",FAST_MATH,"<<fastMath<<",undefined_macro_value,-1";
#if defined(_MSC_VER) && defined(_M_X64)
        std::cout<<",MCW_EM,"<<std::hex<<_MCW_EM<<",EM_AMBIGUOUS,"<<_EM_AMBIGUOUS<<std::dec
            <<",control_masks_ok,"<<((r.control&_MCW_EM)==_MCW_EM)
            <<",sse_masks_ok,"<<((r.csr&0x1f80u)==0x1f80u)<<",ambiguous,"<<bool(r.control&_EM_AMBIGUOUS);
#endif
        std::cout<<'\n';announced=true;
    }
    return r;
}
void sameControl(const Environment& before,const Environment& after)
{
    check(before.rounding==after.rounding&&before.control==after.control
        &&(before.csr&~0x3fu)==(after.csr&~0x3fu),"viscosity product changed FP control/rounding");
}
struct OriginalProduct {std::vector<double> weights;Matrix matrix;std::size_t symbolicStoredValues=0,symbolicZeroValues=0;};
OriginalProduct original(const HydroStrainOperator& raw,const std::vector<double>& mu)
{
    OriginalProduct r;check(mu.size()==raw.points().size(),"original product fixture dimensions");
    // Literal original public bilinear form: no symbolic-plan traversal.
    for(std::size_t q=0;q<mu.size();++q){
        const long double product=2*static_cast<long double>(mu[q])*raw.points()[q].physicalWeight;
        const double weight=static_cast<double>(product);
        for(int k=0;k<raw.components();++k)r.weights.push_back(weight);
    }
    Matrix weighted=raw.matrix();
    for(Eigen::Index row=0;row<weighted.outerSize();++row)
        for(Matrix::InnerIterator e(weighted,row);e;++e)e.valueRef()*=r.weights[static_cast<std::size_t>(row)];
    r.matrix=raw.matrix().transpose()*weighted;
    r.symbolicStoredValues=static_cast<std::size_t>(r.matrix.nonZeros());
    for(Eigen::Index k=0;k<r.matrix.nonZeros();++k){std::uint64_t bits=0;
        std::memcpy(&bits,r.matrix.valuePtr()+k,sizeof bits);r.symbolicZeroValues+=(bits<<1)==0;}
    r.matrix.prune(0.0);r.matrix.makeCompressed();return r;
}
void exactMatrix(const Matrix& actual,const Matrix& expected,std::size_t& coefficients)
{
    check(actual.isCompressed()&&expected.isCompressed()&&actual.rows()==expected.rows()
        &&actual.cols()==expected.cols()&&actual.nonZeros()==expected.nonZeros(),"viscosity product dimensions/compression/nnz differ");
    const auto outer=static_cast<std::size_t>(actual.outerSize()+1),nonzeros=static_cast<std::size_t>(actual.nonZeros());
    check(std::memcmp(actual.outerIndexPtr(),expected.outerIndexPtr(),outer*sizeof(Matrix::StorageIndex))==0,
        "viscosity product compressed outer indices differ");
    if(nonzeros)check(std::memcmp(actual.innerIndexPtr(),expected.innerIndexPtr(),nonzeros*sizeof(Matrix::StorageIndex))==0,
        "viscosity product compressed inner indices differ");
    for(std::size_t k=0;k<nonzeros;++k){
        std::uint64_t a=0,b=0;static_assert(sizeof(a)==sizeof(double),"coefficient oracle requires binary64");
        std::memcpy(&a,actual.valuePtr()+k,sizeof(a));std::memcpy(&b,expected.valuePtr()+k,sizeof(b));
        if(a!=b)std::cerr<<"VISCOSITY_PRODUCT_BITS,index,"<<k<<",actual,"<<std::hex<<a<<",expected,"<<b<<std::dec<<'\n';
        check(a==b,"viscosity product coefficient bits differ");++coefficients;
    }
}
void exactWeights(const std::vector<double>& actual,const std::vector<double>& expected)
{
    check(actual.size()==expected.size(),"viscosity original row-weight dimensions differ");
    if(!actual.empty())check(std::memcmp(actual.data(),expected.data(),actual.size()*sizeof(double))==0,
        "viscosity original row-weight bits differ");
}
void scene(bool rz,std::size_t& coefficients)
{
    QuadtreeMesh mesh(input(8,rz?1:0));refine(mesh,1);FaceSegmentGraph graph;graph.build(mesh);
    const auto domain=NativeLeafDomain::build(mesh,graph,rz);const auto layout=HydroDofLayout::build(mesh,domain);
    const auto metrics=HydroPhysicalDivergence::build(mesh,HydroIntegralDivergence::build(mesh,layout));
    std::size_t hanging=0;
    for(int row=0;row<domain.nativeFaceCount();++row){const auto f=domain.face(mesh,row);
        if(f.negativeCell>=0&&f.positiveCell>=0
            &&mesh.cells(mesh.cellIndexToLabel[f.negativeCell]).level()!=mesh.cells(mesh.cellIndexToLabel[f.positiveCell]).level())++hanging;
    }
    check(hanging>0,"product-plan fixture has no actual 2:1 native faces");
    const auto points=HydroStrainOperator::volumeQuadrature(mesh,metrics);
    const auto raw=HydroStrainOperator::build(mesh,metrics,points);const Matrix savedRaw=raw.matrix();
    check(!raw.isStressCompatible(),"product-plan fixture is not actual raw strain");
    std::vector<double> initial(points.size()),changed(points.size());std::size_t zeros=0;
    for(std::size_t q=0;q<points.size();++q){initial[q]=.5+static_cast<double>(q%7)/8;
        changed[q]=q%3==0?0.:.125+static_cast<double>(q%11)/16;zeros+=changed[q]==0;}
    check(zeros>0&&zeros<changed.size()&&initial!=changed,"product-plan changed/zero-mu coverage missing");
    const auto before=environment();
    std::cout<<"VISCOSITY_PRODUCT_ENV,"<<(rz?"RZ":"XY")<<",eligible,"<<before.eligible
        <<",control,"<<std::hex<<before.control<<",mxcsr,"<<before.csr<<std::dec
        <<",round,"<<before.rounding<<",flags,"<<before.flags<<'\n';
    check(before.eligible,"product-plan environment ineligible; this selector cannot certify a cache hit");
    // First nonzero build must match the independent original Eigen expression.
    const auto first=HydroEnergyViscosity::build(raw,initial);const auto afterFirst=environment();sameControl(before,afterFirst);
    const auto firstOracle=original(raw,initial);exactWeights(first.rowWeights(),firstOracle.weights);exactMatrix(first.matrix(),firstOracle.matrix,coefficients);
    // Evaluate the independent original expression first without clearing any
    // FP flags. The cache call must not introduce flags beyond that expression.
    const auto changedOracle=original(raw,changed);const auto beforeReuse=environment();sameControl(before,beforeReuse);
    const auto second=HydroEnergyViscosity::build(raw,changed);const auto afterReuse=environment();sameControl(beforeReuse,afterReuse);
    check(beforeReuse.flags==afterReuse.flags&&(beforeReuse.csr&0x3fu)==(afterReuse.csr&0x3fu),
        "cached product raised flags absent from original Eigen work");
    check(&first.strain().matrix()==&raw.matrix()&&&second.strain().matrix()==&raw.matrix(),
        "product-plan builds lost the actual shared raw owner");
    exactWeights(second.rowWeights(),changedOracle.weights);exactMatrix(second.matrix(),changedOracle.matrix,coefficients);
    exactMatrix(raw.matrix(),savedRaw,coefficients);
    std::cout<<"VISCOSITY_PRODUCT_PLAN,"<<(rz?"RZ":"XY")<<",cells,"<<layout.pressureCount()<<",hanging,"<<hanging
        <<",first_mu_nonzero,"<<initial.size()<<",changed_mu_zero,"<<zeros<<",nnz,"<<second.matrix().nonZeros()
        <<",flags_initial,"<<before.flags<<",flags_before_reuse,"<<beforeReuse.flags<<",flags_after_reuse,"<<afterReuse.flags
        <<",same_raw_owner,1,index_and_coefficient_bits_exact,1\n";
}
void run(){std::size_t coefficients=0;scene(false,coefficients);scene(true,coefficients);
    std::cout<<"VISCOSITY_PRODUCT_PLAN_PASS,scenes,2,coefficient_checks,"<<coefficients
        <<",original_Eigen_oracle,1,changed_and_zero_mu,1,no_environment_mutation,1\n";}
} // namespace ViscosityProductPlanContract

namespace MaxwellLoadContract {
namespace E=LAMP::Electric;namespace I=LAMP::Interface;using Clock=LAMP::TimeLayer;
using Stress=E::ElectricVofMaxwellStress;
void bounded(long double error,long double scale,double& peak,const char* why,long double gate=1e-12L) {
    check(std::isfinite(error)&&std::isfinite(scale)&&scale>=0,why);if(scale==0){check(error==0,why);return;}
    const auto relative=std::abs(error)/scale;peak=std::max(peak,static_cast<double>(relative));check(relative<=gate,why);
}
template<class F>void rejected(F f,const char* why){bool caught=false;try{f();}catch(const std::exception&){caught=true;}check(caught,why);}
double potential(Vector2D p,bool rz,bool jump) {
    if(jump){const double d=1/(.5/2.+.5/5.);return 1-d*(p.x<=.5?p.x/2.:.5/2.+(p.x-.5)/5.);}
    return .17+(rz?p.y:p.x+.37*p.y);
}
namespace FV=LAMP::Geometry::Detail::FiniteVolume;
void exactScalingChecks()
{
    const long double tiny=std::numeric_limits<long double>::denorm_min(),maximum=std::numeric_limits<long double>::max();
    const FV::Dyadic one(1),zero;FV::Dyadic low(1);low.add(tiny);const auto saved=low;
    check(low.multiplied(.5L).multiplied(2).compare(low)==0&&low.compare(saved)==0,
        "F12 dyadic scaling lost low bits or mutated the receiver");
    check(FV::Dyadic(tiny).multiplied(.5L).multiplied(2).rounded()==tiny,
        "F12 scaling rounded an underflowing intermediate");
    check(FV::Dyadic(maximum).multiplied(2).multiplied(.25L).rounded()==maximum*.5L,
        "F12 scaling rounded/rejected an overflowing intermediate with finite final result");
    check(one.multiplied(-3).compare(FV::Dyadic(-3))==0&&one.multiplied(0).compare(zero)==0
        &&zero.multiplied(-2).compare(zero)==0,"F12 negative/zero exact scaling changed");
    bool caught=false;try{zero.multiplied(std::numeric_limits<long double>::infinity());}catch(const std::exception&){caught=true;}
    check(caught,"F12 zero receiver silently accepted nonfinite multiplier");
}
// Independent literal Green quadrature reference. The old Ehat form remains
// only in the fixed pure-gas probe comparison below, never a production switch.
struct ReferenceMaxwellWork
{
    long double momentumWork=0,volumeWork=0,boundaryWork=0,residual=0,relativeResidual=0,absoluteScale=0;
};
long double nativePower(const HydroMaxwellLoad& load,const HydroNativeValues& values)
{
    const auto& u=values.normalVelocity();const auto& f=load.momentumRhsLoad();
    check(u.size()==f.size(),"Maxwell actual native power size mismatch");FV::Dyadic sum;
    for(std::size_t i=0;i<u.size();++i)sum.add(FV::Dyadic::product({u[i],f[i]}));
    return sum.rounded();
}
struct LiteralGreenWork
{
    struct Term {int column;double weight,sample,coefficient;};
    const HydroMaxwellLoad& load;
    std::vector<std::array<double,4>> tensor;
    std::vector<std::array<double,2>> traction;
    std::vector<Term> volumeTerms,boundaryTerms;
    explicit LiteralGreenWork(const HydroMaxwellLoad& actual):load(actual)
    {
        check(!load.isSharp()&&load.report().pairing==HydroMaxwellPairing::Rt0Green,"A load did not declare RT0 Green pairing");
        const auto& g=load.geometry();const auto& layout=g.layout();const bool rz=layout.domain().axisymmetric();
        const auto& points=g.points();const auto& faces=g.stressFacePoints();tensor.resize(points.size());traction.resize(faces.size());
        using Matrix=HydroStrainOperator::Matrix;
        const auto append=[&](std::vector<Term>& terms,const Matrix& matrix,Eigen::Index row,double weight,double sample){
            for(Matrix::InnerIterator entry(matrix,row);entry;++entry)terms.push_back({static_cast<int>(entry.col()),weight,sample,entry.value()});};
        for(std::size_t q=0;q<points.size();++q){const auto sample=load.stress().sample(layout.nativeCell(points[q].cellDof),points[q].point);const auto t=sample.tensor;
            const auto shear=FV::roundDouble(FV::Dyadic::product({std::sqrt(2.),t.xy})).value;
            tensor[q]={t.xx,t.yy,rz?t.hoop:shear,shear};
            append(volumeTerms,g.stressVolumeVelocityMatrix(),static_cast<Eigen::Index>(q)*2,points[q].physicalWeight,sample.tensorDivergence.x);
            append(volumeTerms,g.stressVolumeVelocityMatrix(),static_cast<Eigen::Index>(q)*2+1,points[q].physicalWeight,sample.tensorDivergence.y);}
        std::size_t sidePoints=0;
        for(std::size_t q=0;q<faces.size();++q){const auto& f=faces[q];if(f.physicalWeight==0)continue;
            const bool vertical=layout.normalAxis(f.velocityDof)==HydroDofLayout::NormalAxis::First;
            for(int side=0;side<2;++side){const int p=side==0?f.negativeCellDof:f.positiveCellDof;if(p<0)continue;++sidePoints;
                const auto t=load.stress().sample(layout.nativeCell(p),f.point).tensor;
                const std::array<double,2> canonical{vertical?t.xx:t.xy,vertical?t.xy:t.yy};
                for(int axis=0;axis<2;++axis)append(volumeTerms,g.stressFaceVelocityMatrix(),static_cast<Eigen::Index>(q)*2+axis,
                    (side==0?-1.:1.)*f.physicalWeight,canonical[axis]);
                if(f.negativeCellDof<0||f.positiveCellDof<0){const double sign=side==0?1.:-1.;
                    for(int axis=0;axis<2;++axis){traction[q][axis]=sign*canonical[axis];
                        append(boundaryTerms,g.stressFaceVelocityMatrix(),static_cast<Eigen::Index>(q)*2+axis,f.physicalWeight,traction[q][axis]);}}}}
        check(sidePoints==load.report().cellSidePoints,"Green pairing omitted/duplicated an actual positive-area cell side");
    }
    ReferenceMaxwellWork evaluate(const QuadtreeMesh& mesh,const HydroNativeValues& values) const
    {
        load.validateFor(mesh,load.pressureLoad(),load.phaseGeometry(),load.stress(),load.phaseGeometry().timeLayer());
        const auto& g=load.geometry();values.validateFor(g.layout(),mesh);const auto& u=values.normalVelocity();
        for(int axis:g.metrics().axisVelocityDofs())check(u[axis]==0,"F12 literal reference nonzero radial axis");
        FV::Dyadic momentum,volume,boundary,scale,budget;
        const auto add=[&](FV::Dyadic& sum,std::initializer_list<long double> factors){const auto value=FV::Dyadic::product(factors);sum.add(value);scale.add(value.absolute());};
        for(std::size_t i=0;i<u.size();++i){add(momentum,{u[i],load.momentumRhsLoad()[i]});budget.add(FV::Dyadic::product({1e-12L,std::abs(u[i]),std::abs(load.momentumRhsLoad()[i])}));}
        const auto contract=[&](FV::Dyadic& sum,const std::vector<Term>& terms){for(const auto& t:terms){
            add(sum,{t.weight,t.sample,t.coefficient,u[t.column]});
            budget.add(FV::Dyadic::product({1e-12L,std::abs(t.weight),std::abs(t.sample),std::abs(t.coefficient),std::abs(u[t.column])}));}};
        contract(volume,volumeTerms);contract(boundary,boundaryTerms);
        auto residual=momentum;residual.add(volume,-1);residual.add(boundary,-1);
        check(budget.compare(scale.multiplied(1e-12L))==0,"F12 literal original tau budget differs from tau times full absolute scale");
        check(residual.absolute().compare(budget)<=0,"F12 literal per-term 1e-12 Maxwell work gate failed");
        ReferenceMaxwellWork result;result.momentumWork=momentum.rounded();result.volumeWork=volume.rounded();result.boundaryWork=boundary.rounded();
        result.residual=residual.rounded();result.absoluteScale=scale.rounded();result.relativeResidual=result.absoluteScale==0?0:residual.absolute().rounded()/result.absoluteScale;
        return result;
    }
};
void compareOriginalWork(const QuadtreeMesh& mesh,const HydroMaxwellLoad& load,const HydroNativeValues& velocity,int kind)
{
    const LiteralGreenWork original(load);const auto reference=original.evaluate(mesh,velocity);
    const auto actual=nativePower(load,velocity);
    check(actual==reference.momentumWork&&std::signbit(actual)==std::signbit(reference.momentumWork),
        "F12 canonical RHS exact-dot power differs from independent Green quadrature reference");
    std::cout<<"MAXWELL_RHS_POWER_REFERENCE,kind,"<<kind<<",literal_Green_terms,"<<original.volumeTerms.size()+original.boundaryTerms.size()
        <<",independent_virtual_work_gate,1e-12,retired_public_work_port,1,no_microbenchmark\n";
}
void greenLoad(const QuadtreeMesh& mesh,const HydroMaxwellLoad& actual,const HydroMaxwellLoad& full,const HydroNativeValues& velocity){
    using namespace Rt0GreenContract;actual.validateFor(mesh,full.pressureLoad(),full.phaseGeometry(),full.stress(),full.phaseGeometry().timeLayer());
    same(actual.geometry(),full.geometry());vector(actual.momentumRhsLoad(),full.momentumRhsLoad());
    const LiteralGreenWork original(actual);const int n=actual.geometry().layout().velocityCount();std::vector<FV::Dyadic> volume(n),boundary(n);
    for(const auto& t:original.volumeTerms)volume[t.column].add(FV::Dyadic::product({t.weight,t.sample,t.coefficient}));
    for(const auto& t:original.boundaryTerms)boundary[t.column].add(FV::Dyadic::product({t.weight,t.sample,t.coefficient}));
    for(int i=0;i<n;++i){auto complete=volume[i];complete.add(boundary[i]);scalar(actual.momentumRhsLoad()[i],complete.roundedDouble(),"Green force independent complete RNE");}
    const auto reference=original.evaluate(mesh,velocity);const auto power=nativePower(actual,velocity);
    check(power==reference.momentumWork&&std::signbit(power)==std::signbit(reference.momentumWork),"Green actual RHS exact-dot power/sign");
    check(actual.report().pairing==HydroMaxwellPairing::Rt0Green&&reference.relativeResidual<=1e-12L,"Green pairing or independent original virtual-work proof");++loadCases;
}
void greenRequests(const QuadtreeMesh& mesh,const HydroMaxwellLoad& full,const HydroNativeValues& velocity,bool ownThreads){
    const auto cold=HydroMaxwellLoad::prepareRt0Green(mesh,full.pressureLoad(),full.phaseGeometry(),full.stress());
    Rt0GreenContract::capabilities(mesh,cold.geometry());greenLoad(mesh,cold,full,velocity);
    auto wrong=full.phaseGeometry().timeLayer();++wrong.step;
    Rt0GreenContract::missing([&]{cold.validateFor(mesh,full.pressureLoad(),full.phaseGeometry(),full.stress(),wrong);});
    Rt0GreenContract::missing([&]{(void)HydroMaxwellLoad::prepare(mesh,full.pressureLoad(),full.phaseGeometry(),full.stress(),full.geometry().canonicalRawStrain());});
    const auto prepared=HydroMaxwellLoad::prepare(mesh,full.pressureLoad(),full.phaseGeometry(),full.stress(),cold.geometry());
    check(&prepared.geometry().stressVolumeVelocityMatrix()==&cold.geometry().stressVolumeVelocityMatrix()&&&prepared.geometry().stressFaceVelocityMatrix()==&cold.geometry().stressFaceVelocityMatrix(),"prepared Green load lost actual Rv/B owner");greenLoad(mesh,prepared,full,velocity);
    if(!ownThreads)return;
    using Floating=Rt0GreenContract::Floating;const Floating caller;
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
    const int oldThreads=omp_get_max_threads(),oldDynamic=omp_get_dynamic();
#endif
    std::optional<HydroMaxwellLoad> reference;std::optional<Floating> exit;
    for(int threads:{1,2,8}){
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
        omp_set_dynamic(0);omp_set_num_threads(threads);
#endif
        caller.restore();const auto load=HydroMaxwellLoad::prepareRt0Green(mesh,full.pressureLoad(),full.phaseGeometry(),full.stress());
        const Floating after;
        if(!reference){reference.emplace(load);exit.emplace(after);}else{exit->same(after);Rt0GreenContract::same(load.geometry(),reference->geometry());Rt0GreenContract::vector(load.momentumRhsLoad(),reference->momentumRhsLoad());++Rt0GreenContract::ownFpPairs;}
        greenLoad(mesh,load,full,velocity);
    }
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
    omp_set_num_threads(oldThreads);omp_set_dynamic(oldDynamic);
#endif
    caller.restore();
    std::cout<<"RT0_GREEN_OWN_FP_PASS metric="<<(full.geometry().layout().domain().axisymmetric()?"RZ":"XY")<<" requests=3 thread_requests=1/2/8 complete_FP_errno_pairs=2 independent_literal_force_and_work=1\n";
}
// Same real RT0 quadrature and native velocity as the existing Maxwell scene;
// no synthetic force, phase owner, or geometry fixture is introduced.
void requestedGradientParity(const QuadtreeMesh& mesh,const HydroMaxwellLoad& recovered,
    const HydroNativeValues& velocity,const E::ElectricVofField& omittedField)
{
    bool missing=false;try{(void)omittedField.cellPotentialGradients();}
    catch(const std::invalid_argument& e){missing=std::string(e.what()).find("cell potential gradients were not requested")!=std::string::npos;}
    check(missing,"Maxwell Omit exposed fabricated cell-gradient output");
    const auto omittedStress=Stress::prepare(mesh,omittedField);
    const auto omitted=HydroMaxwellLoad::prepare(mesh,recovered.pressureLoad(),recovered.phaseGeometry(),omittedStress,recovered.geometry());
    omitted.validateFor(mesh,recovered.pressureLoad(),recovered.phaseGeometry(),omittedStress,recovered.phaseGeometry().timeLayer());
    const auto equalVector=[](const std::vector<double>& a,const std::vector<double>& b){
        check(a.size()==b.size()&&(a.empty()||std::memcmp(a.data(),b.data(),a.size()*sizeof(double))==0),
            "Omit changed actual canonical Maxwell native load bits");};
    equalVector(recovered.momentumRhsLoad(),omitted.momentumRhsLoad());
    const auto equalSample=[&](int native,Vector2D point){
        const auto a=recovered.stress().sample(native,point),b=omittedStress.sample(native,point);
        const std::array<double,12> x{a.point.x,a.point.y,a.displacement.x,a.displacement.y,a.electricField.x,a.electricField.y,
            a.tensor.xx,a.tensor.yy,a.tensor.xy,a.tensor.hoop,a.tensorDivergence.x,a.tensorDivergence.y},
            y{b.point.x,b.point.y,b.displacement.x,b.displacement.y,b.electricField.x,b.electricField.y,
            b.tensor.xx,b.tensor.yy,b.tensor.xy,b.tensor.hoop,b.tensorDivergence.x,b.tensorDivergence.y};
        check(a.nativeCell==b.nativeCell&&std::memcmp(x.data(),y.data(),sizeof x)==0,"Omit changed real volume/side RT0 Maxwell sample bits");};
    const auto& g=recovered.geometry();
    for(const auto& p:g.points())equalSample(g.layout().nativeCell(p.cellDof),p.point);
    for(const auto& f:g.stressFacePoints())if(f.physicalWeight>0){
        if(f.negativeCellDof>=0)equalSample(g.layout().nativeCell(f.negativeCellDof),f.point);
        if(f.positiveCellDof>=0)equalSample(g.layout().nativeCell(f.positiveCellDof),f.point);}
    const auto a=LiteralGreenWork(recovered).evaluate(mesh,velocity),b=LiteralGreenWork(omitted).evaluate(mesh,velocity);
    const std::array<long double,6> x{a.momentumWork,a.volumeWork,a.boundaryWork,a.residual,a.relativeResidual,a.absoluteScale},
        y{b.momentumWork,b.volumeWork,b.boundaryWork,b.residual,b.relativeResidual,b.absoluteScale};
    for(std::size_t k=0;k<x.size();++k)check(x[k]==y[k]&&std::signbit(x[k])==std::signbit(y[k]),"Omit changed independent Green virtual-work reference field/sign");
    check(b.relativeResidual<=1e-12L&&nativePower(recovered,velocity)==a.momentumWork&&nativePower(omitted,velocity)==b.momentumWork,"Omit independent virtual work or actual RHS power changed");
}
void scene(int kind,std::size_t& guards,bool exactProbe=false) {
    const bool rz=kind==1,jump=kind==2;QuadtreeMesh mesh(input(8,rz?1:0));
    if(jump){DictionaryMesh::RefineBox box;box.min={.75,0};box.max={1,1};box.level=1;mesh.applyRefineBoxes({box});}else refine(mesh,1);
    FaceSegmentGraph graph;graph.build(mesh);const auto domain=NativeLeafDomain::build(mesh,graph,rz);
    const auto layout=HydroDofLayout::build(mesh,domain);check(!domain.uniformLeaves(),"Hydro Maxwell witness needs actual 2:1 AMR");
    const auto divergence=HydroIntegralDivergence::build(mesh,layout);const auto metrics=HydroPhysicalDivergence::build(mesh,divergence);
    const auto pressureLoad=HydroPressureLoad::build(divergence);
    using Stocks=std::array<long double,2>;std::vector<Stocks> stocks(mesh.cellIndexNum());
    for(int ci=0;ci<mesh.cellIndexNum();++ci){const auto c=mesh.cells(mesh.cellIndexToLabel[ci]).centerCoord();const auto volume=domain.cellVolumes()[domain.cellRowByIndex()[ci]];
        stocks[ci]=c.x<.5?Stocks{volume,0}:Stocks{0,volume};}
    const Clock acceptedClock{0,.2,.2,Clock::Position::Initial},clock{1,.2,.22,Clock::Position::Midpoint};
    const auto accepted=I::GeometricVofPublication::prepareLocalPlicInventory(mesh,domain,acceptedClock,stocks,{});
    const auto pub=I::GeometricVofPublication::prepareLocalPlicInventory(mesh,domain,clock,stocks,{});
    const auto phase=HydroPhaseGeometry::build(mesh,metrics,pub);
    const auto charge=E::ElectricVofCharge::fromIntegrated(mesh,accepted,std::vector<double>(mesh.cellIndexNum()));
    const E::ElectricVofMaterials materials=jump?E::ElectricVofMaterials{2,5,0,0}:E::ElectricVofMaterials{3,3,0,0};
    std::vector<E::ElectricVofBoundary> boundary(domain.nativeFaceCount());std::size_t hanging=0;
    for(int row=0;row<domain.nativeFaceCount();++row){const auto f=domain.face(mesh,row);
        if(f.negativeCell>=0&&f.positiveCell>=0){boundary[row]={E::ElectricVofBoundary::Type::Interior,0};
            const auto a=mesh.cells(mesh.cellIndexToLabel[f.negativeCell]).size(),b=mesh.cells(mesh.cellIndexToLabel[f.positiveCell]).size();
            if(a.x!=b.x||a.y!=b.y){++hanging;if(jump)check(f.center.x>.5,"Hydro Maxwell jump's AMR transition must be remote in same gas phase");}
        }else boundary[row]=f.area==0?E::ElectricVofBoundary{E::ElectricVofBoundary::Type::HomogeneousNeumann,0}:
            E::ElectricVofBoundary{E::ElectricVofBoundary::Type::Dirichlet,potential(f.center,rz,jump)};
    }
    check(hanging>0,"Hydro Maxwell has no genuine hanging faces");
    const auto field=E::ElectricVofField::prepareOhmic(mesh,pub,charge,clock,materials,boundary);
    check((field.report().usedSuperLU!=field.report().usedGMRES)&&field.report().relativeLinearResidual<=1e-11&&field.report().maximumRelativeGaussResidual<=1e-10,"Hydro Maxwell changed retained Electric linear/Gauss gates");
    check(field.chargeCandidate().integratedCharge()==charge.integratedCharge(),"sigma0 Maxwell candidate changed independent Q");
    const auto stress=Stress::prepare(mesh,field);const auto savedPhi=field.potential();const auto savedPhase=phase.cellVolumes();
    std::vector<std::array<double,3>> savedD;for(const auto& f:field.faces())savedD.push_back({f.integratedField,f.displacement,f.current});
    const auto load=HydroMaxwellLoad::prepare(mesh,pressureLoad,phase,stress);load.validateFor(mesh,pressureLoad,phase,stress,clock);
    const auto& g=load.geometry();const auto& points=g.points();const auto& faces=g.stressFacePoints();
    const int nu=layout.velocityCount(),pc=layout.pressureCount();
    check(g.isStressCompatible()&&load.momentumRhsLoad().size()==static_cast<std::size_t>(nu)&&&load.pressureLoad().divergence().matrix()==&phase.metrics().divergence().matrix()&&
        &load.phaseGeometry().cellVolumes()==&phase.cellVolumes()&&&load.stress().cells()==&stress.cells(),"Hydro Maxwell genuine owners/native dimensions lost");
    const LiteralGreenWork reference(load);
    std::vector<long double> scale(nu);
    for(const auto& t:reference.volumeTerms)scale[t.column]+=std::abs(static_cast<long double>(t.weight)*t.sample*t.coefficient);
    for(const auto& t:reference.boundaryTerms)scale[t.column]+=std::abs(static_cast<long double>(t.weight)*t.sample*t.coefficient);
    std::size_t exterior=0;for(const auto& face:faces)if((face.negativeCellDof<0||face.positiveCellDof<0)&&face.physicalWeight>0)++exterior;
    check(exterior==load.report().boundaryPoints&&load.report().volumePoints==points.size(),"Hydro Maxwell physical boundary ownership/count wrong");
    if(rz)check(load.report().axisPoints>0,"Hydro Maxwell RZ omitted true zero-area axis");
    std::vector<double> p(pc);const long double d=1/(.5L/2+.5L/5);
    for(int i=0;i<pc;++i){const auto c=mesh.cells(mesh.cellIndexToLabel[layout.nativeCell(i)]).centerCoord();
        p[i]=jump?static_cast<double>(.5L*d*d/(c.x<.5?2:5)):0.;}
    const HydroNativeValues pv(layout,std::vector<double>(nu),p);const auto cp=pressureLoad.apply(mesh,pv);
    double peakForce=0,peakWall=0,peakWork=0;std::size_t forced=0,reactions=0;
    for(int u=0;u<nu;++u){const auto f=domain.face(mesh,layout.nativeFace(u));const auto actual=load.momentumRhsLoad()[u];
        if(f.negativeCell>=0&&f.positiveCell>=0){bounded(actual-(jump?cp[u]:0),scale[u]+std::abs(cp[u]),peakForce,
                "Hydro Maxwell independent 1D pressure balance/constant-field interior force");if(jump&&actual!=0)++forced;
        }else if(f.area>0){const long double reaction=static_cast<long double>(cp[u])-actual;check(std::isfinite(reaction),"Hydro Maxwell nonfinite closed-wall reaction");
            if(!jump)bounded(actual,scale[u],peakForce,"Hydro Maxwell constant-field force with actual wall traction");
            if(jump){const int cell=f.negativeCell>=0?f.negativeCell:f.positiveCell;const int dof=layout.pressureDof(cell);
                const long double exact=-(f.negativeCell>=0?1.L:-1.L)*f.area*p[dof];
                bounded(reaction-exact,scale[u]+std::abs(cp[u])+std::abs(exact),peakWall,"Hydro Maxwell independent wall constraint reaction");if(reaction!=0)++reactions;}
        }
    }
    if(jump)check(forced>0&&reactions>0,"Hydro Maxwell jump was not a real force plus wall reaction witness");
    std::vector<double> u(nu);for(int i=0;i<nu;++i)u[i]=std::sin(.37*i)+std::cos(.91*i);for(int i:metrics.axisVelocityDofs())u[i]=0;
    const HydroNativeValues velocity(layout,u,std::vector<double>(pc));const auto work=LiteralGreenWork(load).evaluate(mesh,velocity);
    check(nativePower(load,velocity)==work.momentumWork,"Hydro Maxwell actual RHS exact-dot power");
    check(field.cellPotentialGradients().size()==static_cast<std::size_t>(mesh.cellIndexNum()),"default Recover lost original native cell-gradient output");
    const auto omittedField=E::ElectricVofField::prepareOhmic(mesh,pub,charge,clock,materials,boundary,E::ElectricVofControls{},
        E::ElectricVofField::CellGradientRecovery::Omit);
    requestedGradientParity(mesh,load,velocity,omittedField);greenRequests(mesh,load,velocity,kind<2);
    if(exactProbe)compareOriginalWork(mesh,load,velocity,kind);
    long double dotUF=0,volume=0,wall=0,workScale=0;
    for(int i=0;i<nu;++i){const long double term=static_cast<long double>(u[i])*load.momentumRhsLoad()[i];dotUF+=term;workScale+=std::abs(term);}
    for(const auto& t:reference.volumeTerms){const long double term=static_cast<long double>(t.weight)*t.sample*t.coefficient*u[t.column];volume+=term;workScale+=std::abs(term);}
    for(const auto& t:reference.boundaryTerms){const long double term=static_cast<long double>(t.weight)*t.sample*t.coefficient*u[t.column];wall+=term;workScale+=std::abs(term);}
    bounded(dotUF-(volume+wall),workScale,peakWork,"Hydro Maxwell independent RT0 Green/boundary virtual work");
    bounded(work.momentumWork-dotUF,workScale,peakWork,"Hydro Maxwell work receipt did not use actual native F");
    bounded(work.volumeWork-volume,workScale,peakWork,"Hydro Maxwell independent Green volume work differs from actual tensor");
    bounded(work.boundaryWork-wall,workScale,peakWork,"Hydro Maxwell independent boundary work differs from same RT0 tensor");
    check(work.relativeResidual<=1e-12L,"Hydro Maxwell independent virtual-work gate failed");
    const auto zero=HydroEnergyViscosity::zero(mesh,metrics);check(zero.matrix().nonZeros()==0&&zero.strain().matrix().rows()==0&&zero.dissipation(mesh,velocity)==0,
        "Hydro Maxwell inserted a fake viscosity to obtain strain");for(double value:zero.apply(mesh,velocity))check(value==0,"Hydro Maxwell changed true mu0 load");
    if(kind==0){auto wrong=clock;++wrong.step;rejected([&]{load.validateFor(mesh,pressureLoad,phase,stress,wrong);},"Hydro Maxwell same-time wrong fullclock accepted");++guards;
        const auto freshDiv=HydroIntegralDivergence::build(mesh,layout);const auto freshP=HydroPressureLoad::build(freshDiv);
        rejected([&]{HydroMaxwellLoad::prepare(mesh,freshP,phase,stress);},"Hydro Maxwell fresh DI with same coefficients accepted");++guards;
        Rt0GreenContract::missing([&]{(void)HydroMaxwellLoad::prepareRt0Green(mesh,freshP,phase,stress);});
        const auto freshPhase=HydroPhaseGeometry::build(mesh,metrics,pub);
        rejected([&]{load.validateFor(mesh,pressureLoad,freshPhase,stress,clock);},"Hydro Maxwell fresh phase measurement owner accepted");++guards;
        const auto freshPub=I::GeometricVofPublication::prepareLocalPlicInventory(mesh,domain,clock,stocks,{});
        const auto foreignPhase=HydroPhaseGeometry::build(mesh,metrics,freshPub);
        rejected([&]{HydroMaxwellLoad::prepare(mesh,pressureLoad,foreignPhase,stress);},"Hydro Maxwell foreign VOF owner accepted");++guards;
        const auto freshStress=Stress::prepare(mesh,field);
        rejected([&]{load.validateFor(mesh,pressureLoad,phase,freshStress,clock);},"Hydro Maxwell fresh stress owner accepted");++guards;
        QuadtreeMesh foreignMesh(input(4,0));rejected([&]{load.validateFor(foreignMesh,pressureLoad,phase,stress,clock);},"Hydro Maxwell foreign mesh accepted");++guards;
        auto end=clock;end.position=Clock::Position::End;const auto endPub=I::GeometricVofPublication::prepareLocalPlicInventory(mesh,domain,end,stocks,{});
        const auto endPhase=HydroPhaseGeometry::build(mesh,metrics,endPub);const auto endField=E::ElectricVofField::prepareOhmic(mesh,endPub,charge,end,materials,boundary);
        const auto endStress=Stress::prepare(mesh,endField);
        rejected([&]{HydroMaxwellLoad::prepare(mesh,pressureLoad,endPhase,endStress);},"Hydro Maxwell genuine coeval End force accepted as Midpoint");++guards;
        Rt0GreenContract::missing([&]{(void)HydroMaxwellLoad::prepareRt0Green(mesh,pressureLoad,endPhase,endStress);});
    }
    check(savedPhi==field.potential()&&charge.integratedCharge()==field.chargeCandidate().integratedCharge(),"Hydro Maxwell changed accepted electric state");
    for(std::size_t row=0;row<field.faces().size();++row)check(savedD[row]==std::array<double,3>{field.faces()[row].integratedField,field.faces()[row].displacement,field.faces()[row].current},"Hydro Maxwell changed H/D/J");
    for(std::size_t i=0;i<savedPhase.size();++i)check(savedPhase[i].liquid==phase.cellVolumes()[i].liquid&&savedPhase[i].gas==phase.cellVolumes()[i].gas,"Hydro Maxwell changed independent phase geometry");
    std::cout<<"HYDRO_MAXWELL_LOAD,"<<(rz?"RZ":"XY")<<",jump,"<<jump<<",cells,"<<pc<<",hanging,"<<hanging<<",volume_points,"<<load.report().volumePoints
        <<",exterior_points,"<<load.report().boundaryPoints<<",axis_points,"<<load.report().axisPoints<<",Cp_F,"<<peakForce<<",wall_reaction,"<<peakWall
        <<",virtual_work,"<<peakWork<<",native_work,"<<work.relativeResidual<<",mu0,exact0,no_NS_or_App_commit,requested_gradient_RT0_load_work_parity,1\n";
}
void run(){std::size_t guards=0;scene(0,guards);scene(1,guards);scene(2,guards);check(guards==7,"Hydro Maxwell remaining prepare/owner/fullclock guards missing");
    std::cout<<"HYDRO_MAXWELL_LOAD_PASS,scenes,3,guards,"<<guards<<",same_C,actual_boundary,virtualwork_1e-12,mu0_exact0,retired_work_only_guards,2,no_NS_or_App_commit\n";}
void runPreaggregated(){exactScalingChecks();std::size_t guards=0;scene(0,guards,true);scene(1,guards,true);scene(2,guards,true);
    check(guards==7,"F12 remaining typed Maxwell guards missing");
    std::cout<<"MAXWELL_WORK_PREAGGREGATE_PASS,scenes,3,remaining_source_guards,7,retired_work_only_guards,2,independent_Green_and_rhs_power,1,no_step_speed_claim\n";}
} // namespace MaxwellLoadContract

namespace MaxwellPureGasSolenoidalContract {
namespace E=LAMP::Electric;namespace I=LAMP::Interface;using Clock=LAMP::TimeLayer;
constexpr double inclusionRadius=.23,centerZ=.513,contrast=2.5;
using MaxwellLoadContract::bounded;
bool sameClock(Clock a,Clock b){return a.step==b.step&&a.begin==b.begin&&a.end==b.end&&a.position==b.position;}
double analyticPotential(Vector2D p,bool rz)
{
    const double x=p.x-(rz?0.:.47),z=p.y-centerZ,s=std::hypot(x,z);const int m=rz?3:2;
    if(s<inclusionRadius)return -m*z/(contrast+m-1);
    return -z+(contrast-1)/(contrast+m-1)*std::pow(inclusionRadius/s,m)*z;
}
struct Probe {
    std::vector<double> velocity;
    std::vector<long double> analyticBalance,physicalBalance,gross;
    std::size_t nonzeroFaces=0,hangingNonzeroFaces=0;
    double areaDefect=0,divergenceDefect=0;
};
Probe cycle(const QuadtreeMesh& mesh,const HydroDofLayout& layout,Vector2D site,bool rz)
{
    Probe p;p.velocity.resize(layout.velocityCount());p.analyticBalance.resize(layout.pressureCount());
    p.physicalBalance.resize(layout.pressureCount());p.gross.resize(layout.pressureCount());
    const auto psi=[site](Vector2D x){return x.x==site.x&&x.y==site.y?1.L:0.L;};
    for(int u=0;u<layout.velocityCount();++u){const auto f=layout.domain().face(mesh,layout.nativeFace(u));
        const auto& edge=f.horizontal?mesh.hEdges(f.edgeLabel):mesh.vEdges(f.edgeLabel);
        auto a=mesh.verts(edge.edgeToVert[0]).coord(),b=mesh.verts(edge.edgeToVert[1]).coord();
        if(f.horizontal?b.x<a.x:b.y<a.y)std::swap(a,b);
        const long double length=f.horizontal?static_cast<long double>(b.x)-a.x:static_cast<long double>(b.y)-a.y;
        const long double area=length*(rz?static_cast<long double>(CartesianMeasure::twoPi)*(f.horizontal?(static_cast<long double>(a.x)+b.x)/2:a.x):1.L);
        bounded(area-f.area,area,p.areaDefect,"pure-gas independent native physical measure");
        const long double q=(f.horizontal?-1.L:1.L)*(rz?static_cast<long double>(CartesianMeasure::twoPi):1.L)*(psi(b)-psi(a));
        if(area==0){check(rz&&q==0,"compact cycle placed flux on the zero-area true axis");continue;}
        p.velocity[u]=static_cast<double>(q/area);
        if(q!=0){check(f.negativeCell>=0&&f.positiveCell>=0,"compact cycle touches a physical boundary");++p.nonzeroFaces;
            if(mesh.cells(mesh.cellIndexToLabel[f.negativeCell]).level()!=mesh.cells(mesh.cellIndexToLabel[f.positiveCell]).level())++p.hangingNonzeroFaces;}
        const long double actual=area*p.velocity[u];
        for(int side=0;side<2;++side){const int native=side==0?f.negativeCell:f.positiveCell;if(native<0)continue;
            const int cell=layout.pressureDof(native);const int sign=side==0?1:-1;
            p.analyticBalance[cell]+=sign*q;p.physicalBalance[cell]+=sign*actual;p.gross[cell]+=std::abs(actual);}
    }
    check(p.nonzeroFaces>=3&&p.nonzeroFaces<=8,"compact streamfunction vertex has no bounded actual internal native cycle");
    if(layout.domain().uniformLeaves())check(p.nonzeroFaces==4,"uniform compact cycle does not have four internal native faces");
    for(std::size_t c=0;c<p.gross.size();++c){check(p.analyticBalance[c]==0,"endpoint streamfunction flux does not telescope independently");
        bounded(p.physicalBalance[c],p.gross[c],p.divergenceDefect,"represented compact cycle is not physically divergence free");}
    return p;
}
void probe(const QuadtreeMesh& mesh,const HydroMaxwellLoad& load,Vector2D site,bool rz,int number)
{
    const auto& phase=load.phaseGeometry();const auto& layout=phase.layout();const auto& geometry=load.geometry();
    const auto& field=load.stress().field();const auto& pub=phase.vofPublication();
    auto compact=cycle(mesh,layout,site,rz);const auto& u=compact.velocity;const int pc=layout.pressureCount();
    const HydroNativeValues values(layout,u,std::vector<double>(pc));
    const auto actualDiv=load.pressureLoad().divergence().apply(mesh,values);
    for(int p=0;p<pc;++p)bounded(actualDiv[p],compact.gross[p],compact.divergenceDefect,"compact cycle fails the authoritative DI");
    std::vector<unsigned char> affected(mesh.cellIndexNum());
    const auto touches=[&](const HydroStrainOperator::Matrix& m,Eigen::Index row){
        for(HydroStrainOperator::Matrix::InnerIterator e(m,row);e;++e)if(e.value()!=0&&u[e.col()]!=0)return true;return false;};
    for(std::size_t q=0;q<geometry.points().size();++q){bool active=false;
        for(int k=0;k<geometry.components();++k)active=active||touches(geometry.matrix(),static_cast<Eigen::Index>(q)*geometry.components()+k);
        for(int k=0;k<2;++k)active=active||touches(geometry.stressVolumeVelocityMatrix(),static_cast<Eigen::Index>(q)*2+k);
        if(active)affected[layout.nativeCell(geometry.points()[q].cellDof)]=1;}
    for(std::size_t q=0;q<geometry.stressFacePoints().size();++q){
        if(!touches(geometry.stressFaceVelocityMatrix(),static_cast<Eigen::Index>(q)*2)&&
            !touches(geometry.stressFaceVelocityMatrix(),static_cast<Eigen::Index>(q)*2+1))continue;
        const auto& f=geometry.stressFacePoints()[q];
        check((f.negativeCellDof>=0&&f.positiveCellDof>=0)||(rz&&f.point.x==0&&f.physicalWeight==0),
            "compact raw/shared trace stencil reaches positive-area physical boundary traction");
        if(f.negativeCellDof>=0)affected[layout.nativeCell(f.negativeCellDof)]=1;
        if(f.positiveCellDof>=0)affected[layout.nativeCell(f.positiveCellDof)]=1;
    }
    std::size_t supportCells=0,axisSupportCells=0;double minimumGap=std::numeric_limits<double>::infinity(),rt0Gauss=0;
    for(int native=0;native<mesh.cellIndexNum();++native)if(affected[native]){++supportCells;const auto& cell=pub.plicCells()[native];
        check(cell.inventory[0]==0&&cell.inventory[1]>0&&cell.resolved,"compact reconstruction support includes liquid/mixed/unresolved geometry");
        check(field.epsilonByCell()[native]==1&&field.chargeCandidate().integratedCharge()[native]==0&&field.sourceCharge().integratedCharge()[native]==0,
            "compact support is not exact constant-epsilon zero-free-Q gas");
        const auto& b=cell.rectangle;const double dx=b.lower.x-(rz?0.:.47),dy=b.lower.y-centerZ;
        const double dz=b.upper.y-centerZ;const double closestZ=dy>0?dy:(dz<0?-dz:0.);
        const double nearX=dx>0?dx:(b.upper.x<(rz?0.:.47)?(rz?0.:.47)-b.upper.x:0.);
        minimumGap=std::min(minimumGap,std::hypot(nearX,closestZ)-inclusionRadius);if(rz&&b.lower.x==0)++axisSupportCells;
        const auto& c=load.stress().cells()[native];const auto& d=c.sideDisplacement;
        const long double balance=d[1]-d[0]+d[3]-d[2],gross=std::abs(d[0])+std::abs(d[1])+std::abs(d[2])+std::abs(d[3]);
        bounded(balance,gross,rt0Gauss,"pure-gas whole-side RT0 Gauss disagrees with zero free charge",1e-10L);
    }
    check(supportCells>0&&minimumGap>0,"compact support has no independent positive separation from analytic inclusion");
    if(rz&&number==2)check(compact.hangingNonzeroFaces>0,"RZ compact probe misses its actual2:1 internal child faces");
    if(rz&&number==3)check(axisSupportCells>0,"RZ near-axis probe does not consume true axis-adjacent reconstruction support");
    std::vector<double> pressure(pc);for(int p=0;p<pc;++p)pressure[p]=.13+std::sin(.19*p);
    const auto cp=load.pressureLoad().apply(mesh,HydroNativeValues(layout,std::vector<double>(layout.velocityCount()),pressure));
    long double cpWork=0,cpScale=0,adjoint=0;
    for(std::size_t j=0;j<u.size();++j){cpWork+=static_cast<long double>(u[j])*cp[j];cpScale+=std::abs(static_cast<long double>(u[j])*cp[j]);}
    for(int p=0;p<pc;++p)adjoint+=static_cast<long double>(pressure[p])*actualDiv[p];
    double nullDefect=0;bounded(cpWork,cpScale,nullDefect,"analytic Cp-force compact divfree work baseline is not zero");
    bounded(cpWork+adjoint,cpScale,nullDefect,"compact probe lost same-C pressure adjoint identity");
    const auto work=MaxwellLoadContract::LiteralGreenWork(load).evaluate(mesh,values);check(work.boundaryWork==0,"certified internal probe has nonzero independent physical boundary work");
    check(geometry.hasCanonicalRawStrain(),"pure-gas work split needs the actual retained raw P2 geometry");
    const auto strain=geometry.apply(mesh,values),raw=geometry.canonicalRawStrain().apply(mesh,values);
    long double rawPower=0,meanCorrectionPower=0;
    for(std::size_t q=0;q<geometry.points().size();++q){const auto& point=geometry.points()[q];const auto t=load.stress().sample(layout.nativeCell(point.cellDof),point.point).tensor;
        const std::array<double,4> packed{t.xx,t.yy,rz?t.hoop:std::sqrt(2.)*t.xy,std::sqrt(2.)*t.xy};
        for(int k=0;k<geometry.components();++k){const auto row=q*geometry.components()+k;
            const long double factor=-static_cast<long double>(point.physicalWeight)*packed[k];
            rawPower+=factor*raw[row];meanCorrectionPower+=factor*(static_cast<long double>(strain[row])-raw[row]);
        }}
    const long double legacyEhatPower=rawPower+meanCorrectionPower;
    const Eigen::Map<const Eigen::VectorXd> native(u.data(),static_cast<Eigen::Index>(u.size()));
    const Eigen::VectorXd rawVelocity=geometry.stressVolumeVelocityMatrix()*native,sharedVelocity=geometry.stressFaceVelocityMatrix()*native;
    long double divergencePower=0,sidePower=0,greenScale=0;double derivativeDefect=0;
    for(std::size_t q=0;q<geometry.points().size();++q){const auto& point=geometry.points()[q];const int cell=layout.nativeCell(point.cellDof);
        const auto sample=load.stress().sample(cell,point.point);const auto& c=load.stress().cells()[cell];
        const auto density=[&](int side){return c.sidePhysicalMeasure[side]==0?0.L:c.sideDisplacement[side]/c.sidePhysicalMeasure[side];};
        const long double dx=static_cast<long double>(c.upper.x)-c.lower.x,dz=static_cast<long double>(c.upper.y)-c.lower.y;
        const long double independentDivD=(rz?2*(c.upper.x*density(1)-c.lower.x*density(0))/(dx*(c.upper.x+c.lower.x)):
            (density(1)-density(0))/dx)+(density(3)-density(2))/dz;
        for(int axis=0;axis<2;++axis){const double fieldComponent=axis==0?sample.electricField.x:sample.electricField.y,
                actual=axis==0?sample.tensorDivergence.x:sample.tensorDivergence.y;
            const long double expected=independentDivD*fieldComponent;
            // Native whole-side arithmetic differs from the exact represented
            // derivative factors by roundoff. Scale by uncancelled side data.
            const long double derivativeScale=(std::abs(density(0))+std::abs(density(1)))/dx+
                (std::abs(density(2))+std::abs(density(3)))/dz;
            bounded(actual-expected,derivativeScale*std::abs(fieldComponent),derivativeDefect,"same RT0 analytic tensor divergence/hoop derivative changed");
            const long double term=static_cast<long double>(point.physicalWeight)*actual*rawVelocity[2*q+axis];
            divergencePower+=term;greenScale+=std::abs(term);}}
    std::size_t actualSidePoints=0;
    for(std::size_t q=0;q<geometry.stressFacePoints().size();++q){const auto& f=geometry.stressFacePoints()[q];if(f.physicalWeight==0)continue;
        const bool vertical=layout.normalAxis(f.velocityDof)==HydroDofLayout::NormalAxis::First;
        for(int side=0;side<2;++side){const int cell=side==0?f.negativeCellDof:f.positiveCellDof;if(cell<0)continue;++actualSidePoints;
            const auto t=load.stress().sample(layout.nativeCell(cell),f.point).tensor;
            const std::array<double,2> canonical{vertical?t.xx:t.xy,vertical?t.xy:t.yy};
            for(int axis=0;axis<2;++axis){const long double term=(side==0?-1.L:1.L)*f.physicalWeight*canonical[axis]*sharedVelocity[2*q+axis];
                sidePower+=term;greenScale+=std::abs(term);}}}
    if(rz)for(int cell=0;cell<mesh.cellIndexNum();++cell){const auto& c=load.stress().cells()[cell];if(c.lower.x!=0)continue;
        const auto axis=load.stress().sample(cell,{0,(c.lower.y+c.upper.y)/2});
        check(axis.displacement.x==0&&axis.tensorDivergence.x==0&&std::isfinite(axis.tensorDivergence.y),"same RT0 Green axis limit lost hoop cancellation");}
    check(load.report().pairing==HydroMaxwellPairing::Rt0Green&&load.report().cellSidePoints==actualSidePoints,
        "pure-gas consumer did not retain all actual-side Green pairing");
    double splitDefect=0;bounded(work.volumeWork-divergencePower-sidePower,greenScale,splitDefect,"actual Green work lost raw-divergence/independent cell-side ledger");
    check(MaxwellLoadContract::nativePower(load,values)==work.momentumWork,
        "fixed-probe actual canonical RHS power differs from independent Green reference");
    double tangentJump=0,tractionJump=0,normalDJump=0;std::size_t jumpFaces=0;
    for(int row=0;row<layout.domain().nativeFaceCount();++row){const auto f=layout.domain().face(mesh,row);
        if(f.negativeCell<0||f.positiveCell<0||!affected[f.negativeCell]||!affected[f.positiveCell])continue;
        const auto a=load.stress().sample(f.negativeCell,f.center),b=load.stress().sample(f.positiveCell,f.center);
        tangentJump=std::max(tangentJump,std::abs(f.horizontal?a.electricField.x-b.electricField.x:a.electricField.y-b.electricField.y));
        normalDJump=std::max(normalDJump,std::abs(f.horizontal?a.displacement.y-b.displacement.y:a.displacement.x-b.displacement.x));
        tractionJump=std::max(tractionJump,f.horizontal?std::hypot(a.tensor.xy-b.tensor.xy,a.tensor.yy-b.tensor.yy):
            std::hypot(a.tensor.xx-b.tensor.xx,a.tensor.xy-b.tensor.xy));++jumpFaces;}
    check(jumpFaces>0,"pure-gas compact patch has no qualified internal traction samples");
    const long double scaled=work.absoluteScale==0?0:std::abs(work.momentumWork)/work.absoluteScale;
    std::cout<<"MAXWELL_PURE_GAS_SOLENOIDAL,"<<(rz?"RZ":"XY")<<",probe,"<<number<<",vertex_x,"<<site.x<<",vertex_y,"<<site.y
        <<",native_faces,"<<compact.nonzeroFaces<<",hanging_nonzero_faces,"<<compact.hangingNonzeroFaces<<",support_cells,"<<supportCells
        <<",axis_support_cells,"<<axisSupportCells<<",analytic_interface_gap,"<<minimumGap
        <<",physical_area_relative,"<<compact.areaDefect<<",DI_relative,"<<compact.divergenceDefect<<",Cp_work,"<<cpWork<<",Cp_null_relative,"<<nullDefect
        <<",actual_Mid_gid,"<<pub.identity().value()<<",Q0,1,epsilon_gas,1,hoop_retained,"<<rz<<",boundary_work,"<<work.boundaryWork
        <<",native_work,"<<work.momentumWork<<",positive_work_scale,"<<work.absoluteScale<<",scaled_work,"<<scaled
        <<",virtual_work_relative,"<<work.relativeResidual<<",raw_P2_power,"<<rawPower<<",Ehat_mean_correction_power,"<<meanCorrectionPower
        <<",legacy_Ehat_power,"<<legacyEhatPower<<",Green_divergence_power,"<<divergencePower<<",Green_cell_side_power,"<<sidePower
        <<",Green_derivative_relative,"<<derivativeDefect<<",Green_work_split_relative,"<<splitDefect<<",same_field_sites_gaps,1,whole_side_Gauss,"<<rt0Gauss<<",internal_midpoint_faces,"<<jumpFaces
        <<",max_midpoint_tangential_E_jump,"<<tangentJump<<",max_midpoint_normal_D_jump,"<<normalDJump
        <<",max_midpoint_traction_jump,"<<tractionJump<<",new_physical_rejection_gate,0\n";
}
void scene(bool rz)
{
    QuadtreeMesh mesh(input(32,rz?1:0));
    if(rz){DictionaryMesh::RefineBox box;box.min={.5,0};box.max={1,1};box.level=1;mesh.applyRefineBoxes({box});}
    FaceSegmentGraph graph;graph.build(mesh);const auto domain=NativeLeafDomain::build(mesh,graph,rz);
    check(domain.uniformLeaves()!=rz,"inclusion compact probe lost its fixed XY or real RZ2:1 mesh");
    const auto initial=I::GeometricVofInitialization::prepare(mesh,domain,[rz](Vector2D p){return std::hypot(p.x-(rz?0.:.47),p.y-centerZ)-inclusionRadius;});
    using Stock=std::array<long double,2>;std::vector<Stock> stocks;for(const auto& p:initial.plicCells())stocks.push_back(p.inventory);
    const Clock clock{1,0,.01,Clock::Position::Midpoint};
    const auto pub=I::GeometricVofPublication::prepareLocalPlicInventory(mesh,domain,clock,stocks,{});
    const auto savedInitial=initial.plicCells(),savedMid=pub.plicCells();
    const auto charge=E::ElectricVofCharge::fromIntegrated(mesh,initial,std::vector<double>(mesh.cellIndexNum()));const auto savedQ=charge.integratedCharge();
    const auto layout=HydroDofLayout::build(mesh,domain);const auto divergence=HydroIntegralDivergence::build(mesh,layout);
    const auto metrics=HydroPhysicalDivergence::build(mesh,divergence);const auto pressure=HydroPressureLoad::build(divergence);
    const auto phase=HydroPhaseGeometry::build(mesh,metrics,pub);const E::ElectricVofMaterials materials{contrast,1,0,0};
    std::vector<E::ElectricVofBoundary> boundaries(domain.nativeFaceCount());
    for(int row=0;row<domain.nativeFaceCount();++row){const auto f=domain.face(mesh,row);
        boundaries[row]=f.negativeCell>=0&&f.positiveCell>=0?E::ElectricVofBoundary{E::ElectricVofBoundary::Type::Interior,0}:
            (f.area==0?E::ElectricVofBoundary{E::ElectricVofBoundary::Type::HomogeneousNeumann,0}:
                E::ElectricVofBoundary{E::ElectricVofBoundary::Type::Dirichlet,analyticPotential(f.center,rz)});}
    const auto field=E::ElectricVofField::prepareOhmic(mesh,pub,charge,clock,materials,boundaries);field.validateFor(mesh,pub,charge,clock);
    check(field.report().relativeLinearResidual<=1e-11&&field.report().maximumRelativeGaussResidual<=1e-10,
        "pure-gas probe changed retained Electric linear/Gauss gates");
    check(sameClock(field.geometryTime(),clock)&&sameClock(field.candidateTime(),clock)&&field.geometry().identity()==pub.identity()&&
        field.sourceCharge().sameInstance(charge)&&field.chargeCandidate().integratedCharge()==savedQ,"pure-gas probe lost actual Mid/source Q0 owners");
    const auto stress=E::ElectricVofMaxwellStress::prepare(mesh,field);const auto load=HydroMaxwellLoad::prepare(mesh,pressure,phase,stress);
    load.validateFor(mesh,pressure,phase,stress,clock);
    std::cout<<"MAXWELL_PURE_GAS_SCENE,"<<(rz?"RZ":"XY")<<",N,32,epsilon_ratio,2.5,K,0,free_Q,0,field_SuperLU,"<<field.report().usedSuperLU
        <<",field_GMRES,"<<field.report().usedGMRES<<",actual_Mid_time,"<<clock.time()<<",actual_Mid_gid,"<<pub.identity().value()
        <<",linear_residual,"<<field.report().relativeLinearResidual<<",Gauss,"<<field.report().maximumRelativeGaussResidual<<"\n";
    if(rz){probe(mesh,load,{.65625,.75},true,0);probe(mesh,load,{.71875,.6875},true,1);
        probe(mesh,load,{.5,.703125},true,2);probe(mesh,load,{.03125,.875},true,3);}
    else{probe(mesh,load,{.8125,.78125},false,0);probe(mesh,load,{.78125,.8125},false,1);}
    for(std::size_t k=0;k<stocks.size();++k){const auto unchanged=[](const auto& a,const auto& b){return a.inventory==b.inventory&&
            a.normal.x==b.normal.x&&a.normal.y==b.normal.y&&a.offset==b.offset&&a.resolved==b.resolved;};
        check(unchanged(initial.plicCells()[k],savedInitial[k])&&unchanged(pub.plicCells()[k],savedMid[k]),"readonly compact probes mutated source stocks or PLIC planes");}
    check(charge.integratedCharge()==savedQ&&field.chargeCandidate().integratedCharge()==savedQ,"readonly compact probes mutated signed Q0");
}
void run(){scene(false);scene(true);std::cout<<"MAXWELL_PURE_GAS_PROBE_CERTIFICATES_PASS,scenes,2,probes,6,RZ_actual2to1_nearaxis,actual_Mid_owners,exact_Q0_constant_epsilon_gas_support,"
    <<"independent_divfree_same_Cp_null_baseline,original_virtualwork_gates,new_physical_threshold,none,no_NS_no_time_step_no_case_change\n";}
} // namespace MaxwellPureGasSolenoidalContract

namespace MaxwellCompiledContract {
namespace FV=LAMP::Geometry::Detail::FiniteVolume;
using Port=FV::Binary64ExactProducts;
namespace E=LAMP::Electric;namespace I=LAMP::Interface;using Clock=LAMP::TimeLayer;
using Matrix=HydroStrainOperator::Matrix;
std::size_t tapes=0,terms=0,guards=0,loads=0,forceRows=0,powerChecks=0,subnormalWitnesses=0;
bool same(double a,double b){return GramCoercivityContract::word(a)==GramCoercivityContract::word(b);}
template<class F>void rejection(F action,const char* literal) {
    bool caught=false;try{action();}catch(const std::invalid_argument& e){caught=std::string(e.what()).find(literal)!=std::string::npos;}
    check(caught,"compiled exact port lost its original explicit rejection");++guards;
}
void tape(const std::vector<std::array<double,3>>& source,int finalZeroSign=-1) {
    std::vector<Port::Factor> factors;std::vector<Port::Value> values(source.size());std::vector<std::uint32_t> rows;
    FV::Dyadic literal,absolute;
    for(std::size_t i=0;i<source.size();++i){const auto& t=source[i];const auto row=source.size()-1-i;
        factors.push_back(Port::compile(t[0],t[1]));values[row]=Port::decode(t[2]);rows.push_back(static_cast<std::uint32_t>(row));
        const auto exact=FV::Dyadic::product({t[0],t[1],t[2]});literal.add(exact);absolute.add(exact.absolute());++terms;}
    const auto actual=Port::sum(factors.data(),rows.data(),rows.size(),values.data(),values.size());
    check(actual.sum.compare(literal)==0&&actual.absolute.compare(absolute)==0,
        "compiled three-factor integer sum or original-term absolute scale changed");
    check(same(actual.sum.roundedDouble(),literal.roundedDouble()),"compiled terminal binary64 RNE/sign differs from original Dyadic");
    if(finalZeroSign>=0)check(actual.sum.roundedDouble()==0&&std::signbit(actual.sum.roundedDouble())==(finalZeroSign!=0),
        "compiled exact cancellation/negative underflow lost declared zero sign");++tapes;
}
void arithmetic() {
    const double tiny=std::numeric_limits<double>::denorm_min(),maximum=std::numeric_limits<double>::max();
    tape({{.5,-3.,.75},{-.25,7.,-.5},{-0.,maximum,-1.},{.5,-3.,-.75},{1.,1.,.125}});
    tape({{1.,1.,1.},{1.,1.,std::ldexp(1.,-53)}}); // even lower endpoint
    tape({{1.,1.,1.},{3.,1.,std::ldexp(1.,-53)}}); // even upper endpoint
    tape({{std::ldexp(1.,900),std::ldexp(1.,900),1.},
        {-std::ldexp(1.,900),std::ldexp(1.,900),1.},{.5,.25,.125}});
    // The individual products overflow binary64, then cancel exactly. The
    // remaining tiny cubed term spans over 6000 binary exponents from max^3.
    tape({{maximum,maximum,maximum},{-maximum,maximum,maximum},{tiny,tiny,tiny}},0);
    tape({{maximum,maximum,maximum},{-maximum,maximum,maximum},{-tiny,tiny,tiny}},1);
    tape({{-std::ldexp(1.,-538),std::ldexp(1.,-538),.25}},1);
    tape({{maximum,maximum,maximum},{-maximum,maximum,maximum}},0);
    tape({{0.,-1.,maximum},{-0.,2.,-tiny}},0);tape({},0);
    // Duplicate routing really consumes the same decoded current row twice.
    const std::array<Port::Factor,3> factors{{Port::compile(2.,3.),Port::compile(-7.,.5),Port::compile(.25,-4.)}};
    const std::array<Port::Value,2> values{{Port::decode(.375),Port::decode(-.125)}};
    const std::array<std::uint32_t,3> rows{{1,0,1}};
    const auto routed=Port::sum(factors.data(),rows.data(),3,values.data(),2);
    FV::Dyadic expected;expected.add(FV::Dyadic::product({2.,3.,-.125}));expected.add(FV::Dyadic::product({-7.,.5,.375}));expected.add(FV::Dyadic::product({.25,-4.,-.125}));
    check(routed.sum.compare(expected)==0,"compiled row routing substituted a term index for actual fresh row");
    rejection([&]{Port::decode(std::numeric_limits<double>::infinity());},"nonfinite compiled product input");
    rejection([&]{Port::compile(1.,std::numeric_limits<double>::quiet_NaN());},"nonfinite compiled product input");
    const std::uint32_t invalid=2;rejection([&]{Port::sum(factors.data(),&invalid,1,values.data(),2);},"compiled product row out of range");
    rejection([&]{Port::sum(nullptr,rows.data(),1,values.data(),2);},"missing compiled product terms");
    StressAmrLocalReuseContract::FloatingEnvironmentRestore restore;
    check(std::fesetround(FE_UPWARD)==0&&std::fegetround()==FE_UPWARD,"compiled arithmetic FP fixture cannot set actual upward mode");
    std::fenv_t before{},after{};check(std::fegetenv(&before)==0,"compiled arithmetic cannot capture caller FP");
#if defined(_MSC_VER) && defined(_M_X64)
    const unsigned csr=_mm_getcsr();
#endif
    const int savedError=errno;const auto directed=Port::sum(factors.data(),rows.data(),3,values.data(),2);
    check(errno==savedError&&std::fegetenv(&after)==0&&std::memcmp(&before,&after,sizeof before)==0,
        "integer-only compiled accumulation changed native FP state/errno");
#if defined(_MSC_VER) && defined(_M_X64)
    check(csr==_mm_getcsr(),"integer-only compiled accumulation changed raw MXCSR");
#endif
    check(directed.sum.compare(expected)==0,"directed FP controls changed exact integer accumulation");
    rejection([&]{directed.sum.roundedDouble();},"binary64 round-to-nearest-even required");
    rejection([&]{expected.roundedDouble();},"binary64 round-to-nearest-even required");restore.restore();
    std::cout<<"MAXWELL_COMPILED_ARITHMETIC,tapes,"<<tapes<<",terms,"<<terms
        <<",routing,actual_rows,exact_original_absolute_scale,1,extreme_span_over_6000,1,tie_even_and_signed_zero,1\n";
}
void literal(const QuadtreeMesh& mesh,const HydroMaxwellLoad& load,const HydroNativeValues& velocity,bool needSubnormal) {
    const MaxwellLoadContract::LiteralGreenWork sample(load);const auto& g=load.geometry();
    const int nu=g.layout().velocityCount();
    std::vector<FV::Dyadic> volume(nu),boundary(nu);std::size_t actualBoundaryTerms=0,witness=0;
    for(const auto& t:sample.volumeTerms){witness+=std::fpclassify(t.sample)==FP_SUBNORMAL;
        volume[static_cast<std::size_t>(t.column)].add(FV::Dyadic::product({t.weight,t.sample,t.coefficient}));}
    for(const auto& t:sample.boundaryTerms){witness+=std::fpclassify(t.sample)==FP_SUBNORMAL;
        boundary[static_cast<std::size_t>(t.column)].add(FV::Dyadic::product({t.weight,t.sample,t.coefficient}));++actualBoundaryTerms;}
    if(needSubnormal){check(witness>0,"fallback fixture has no actual represented subnormal tensor/traction");subnormalWitnesses+=witness;}
    for(int u=0;u<nu;++u){auto combined=volume[u];combined.add(boundary[u]);
        check(same(load.momentumRhsLoad()[u],combined.roundedDouble()),"compiled native force differs from independent complete original-row RNE/sign");++forceRows;}
    const auto reference=sample.evaluate(mesh,velocity);const auto actual=MaxwellLoadContract::nativePower(load,velocity);
    check(actual==reference.momentumWork&&std::signbit(actual)==std::signbit(reference.momentumWork),"compiled actual RHS exact-dot power/sign");++powerChecks;
    check(reference.relativeResidual<=1e-12L,"compiled independent original physical virtual-work gate");
    if(load.report().compiledContraction){check(load.report().compiledVolumeTerms==sample.volumeTerms.size()
        &&load.report().compiledBoundaryTerms==actualBoundaryTerms&&load.report().compiledRetainedBytes>0
        &&load.report().compiledRetainedBytes<=256u*1024u*1024u,"compiled report does not describe actual retained geometry terms");}
    else check(!load.report().compiledGeometryHit&&!load.report().compiledTimingAvailable&&load.report().compiledVolumeTerms==0
        &&load.report().compiledBoundaryTerms==0&&load.report().compiledRetainedBytes==0,"fallback falsely reports executed compiled contraction");
    ++loads;
}
void scene(bool rz) {
    QuadtreeMesh mesh(input(8,rz?1:0));refine(mesh,1);FaceSegmentGraph graph;graph.build(mesh);
    const auto domain=NativeLeafDomain::build(mesh,graph,rz);const auto layout=HydroDofLayout::build(mesh,domain);
    check(!domain.uniformLeaves(),"compiled Maxwell fixture lost actual 2:1 AMR");
    const auto divergence=HydroIntegralDivergence::build(mesh,layout);const auto metrics=HydroPhysicalDivergence::build(mesh,divergence);
    const auto pressure=HydroPressureLoad::build(divergence);using Stocks=std::array<long double,2>;
    std::vector<Stocks> stocks(mesh.cellIndexNum());for(int ci=0;ci<mesh.cellIndexNum();++ci){
        const auto c=mesh.cells(mesh.cellIndexToLabel[ci]).centerCoord();const auto v=domain.cellVolumes()[domain.cellRowByIndex()[ci]];stocks[ci]=c.x<.5?Stocks{v,0}:Stocks{0,v};}
    const Clock begin{0,.2,.2,Clock::Position::Initial},mid{1,.2,.22,Clock::Position::Midpoint};
    const auto accepted=I::GeometricVofPublication::prepareLocalPlicInventory(mesh,domain,begin,stocks,{});
    const auto pub=I::GeometricVofPublication::prepareLocalPlicInventory(mesh,domain,mid,stocks,{});
    const auto phase=HydroPhaseGeometry::build(mesh,metrics,pub);const auto charge=E::ElectricVofCharge::fromIntegrated(mesh,accepted,std::vector<double>(mesh.cellIndexNum()));
    const E::ElectricVofMaterials material{3,3,0,0};const auto geometry=HydroStrainOperator::buildStressCompatible(mesh,metrics);
    std::vector<double> u(layout.velocityCount());for(std::size_t i=0;i<u.size();++i)u[i]=std::sin(.37*i)+std::cos(.91*i);
    for(int axis:metrics.axisVelocityDofs())u[axis]=0;const HydroNativeValues velocity(layout,u,std::vector<double>(layout.pressureCount()));
    const auto field=[&](double amplitude){std::vector<E::ElectricVofBoundary> boundary(domain.nativeFaceCount());
        for(int row=0;row<domain.nativeFaceCount();++row){const auto f=domain.face(mesh,row);
            if(f.negativeCell>=0&&f.positiveCell>=0)boundary[row]={E::ElectricVofBoundary::Type::Interior,0};
            else boundary[row]=f.area==0?E::ElectricVofBoundary{E::ElectricVofBoundary::Type::HomogeneousNeumann,0}:
                E::ElectricVofBoundary{E::ElectricVofBoundary::Type::Dirichlet,amplitude*MaxwellLoadContract::potential(f.center,rz,false)};}
        const auto current=E::ElectricVofField::prepareOhmic(mesh,pub,charge,mid,material,boundary);
        check((current.report().usedSuperLU!=current.report().usedGMRES)&&current.report().relativeLinearResidual<=1e-11&&current.report().maximumRelativeGaussResidual<=1e-10,
            "compiled Maxwell fixture changed retained electric solver/Gauss gates");return current;};
    const auto firstField=field(1.);const auto firstStress=E::ElectricVofMaxwellStress::prepare(mesh,firstField);
    const auto first=HydroMaxwellLoad::prepare(mesh,pressure,phase,firstStress,geometry);literal(mesh,first,velocity,false);
    const auto changedField=field(.5);const auto changedStress=E::ElectricVofMaxwellStress::prepare(mesh,changedField);
    const auto changed=HydroMaxwellLoad::prepare(mesh,pressure,phase,changedStress,geometry);literal(mesh,changed,velocity,false);
    const auto greenGeometry=HydroStrainOperator::buildRt0Green(mesh,metrics);
    const auto firstGreen=HydroMaxwellLoad::prepare(mesh,pressure,phase,firstStress,greenGeometry);
    const auto changedGreen=HydroMaxwellLoad::prepare(mesh,pressure,phase,changedStress,greenGeometry);
    MaxwellLoadContract::greenLoad(mesh,firstGreen,first,velocity);MaxwellLoadContract::greenLoad(mesh,changedGreen,changed,velocity);
    const auto actualTensor=[](const HydroMaxwellLoad& load){const auto& p=load.geometry().points().front();
        const auto t=load.stress().sample(load.geometry().layout().nativeCell(p.cellDof),p.point).tensor;
        return std::array<double,4>{t.xx,t.yy,t.xy,t.hoop};};
    check(&firstGreen.geometry().stressVolumeVelocityMatrix()==&changedGreen.geometry().stressVolumeVelocityMatrix()
        &&&firstGreen.stress().cells()!=&changedGreen.stress().cells()&&actualTensor(firstGreen)!=actualTensor(changedGreen),"Green reused geometry froze physical stress");
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
    check(firstGreen.report().compiledContraction&&changedGreen.report().compiledGeometryHit,"Green actual Rv/B compiled identity did not reuse");
#endif
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
    check(first.report().compiledContraction&&changed.report().compiledContraction&&changed.report().compiledGeometryHit,
        "qualified normal actual-load scene did not execute/reuse compiled geometry plan");
#endif
    check(&first.geometry().matrix()==&geometry.matrix()&&&changed.geometry().matrix()==&geometry.matrix()
        &&&first.stress().cells()!=&changed.stress().cells()&&firstField.potential()!=changedField.potential()
        &&actualTensor(first)!=actualTensor(changed),"compiled geometry reuse froze old stress values or changed actual geometry owner");
    if(!rz){const auto tinyField=field(1e-154);const auto tinyStress=E::ElectricVofMaxwellStress::prepare(mesh,tinyField);
        const auto fallback=HydroMaxwellLoad::prepare(mesh,pressure,phase,tinyStress,geometry);
        const auto tinyGreen=HydroMaxwellLoad::prepare(mesh,pressure,phase,tinyStress,greenGeometry);
        MaxwellLoadContract::greenLoad(mesh,tinyGreen,fallback,velocity);check(!tinyGreen.report().compiledContraction,"Green subnormal bypassed original literal fallback");
        check(!fallback.report().compiledContraction,"actual subnormal tensor load bypassed original serial fallback");literal(mesh,fallback,velocity,true);
        const auto zeroField=field(0.);const auto zeroStress=E::ElectricVofMaxwellStress::prepare(mesh,zeroField);
        const auto zero=HydroMaxwellLoad::prepare(mesh,pressure,phase,zeroStress,geometry);literal(mesh,zero,velocity,false);
        for(double v:zero.momentumRhsLoad())check(same(v,0.),"exact zero stress produced nonpositive-zero force");
        const auto zeroGreen=HydroMaxwellLoad::prepare(mesh,pressure,phase,zeroStress,greenGeometry);
        MaxwellLoadContract::greenLoad(mesh,zeroGreen,zero,velocity);
        for(double v:zeroGreen.momentumRhsLoad())check(same(v,0.),"Green exact zero stress changed positive-zero force");
        auto wrong=mid;++wrong.step;rejection([&]{changed.validateFor(mesh,pressure,phase,changedStress,wrong);},"complete clock");
        rejection([&]{changed.validateFor(mesh,pressure,phase,firstStress,mid);},"different actual Maxwell stress owner");
    }
    std::cout<<"MAXWELL_COMPILED_LOAD,"<<(rz?"RZ":"XY")<<",cells,"<<layout.pressureCount()<<",volume_terms,"<<changed.report().compiledVolumeTerms
        <<",boundary_terms,"<<changed.report().compiledBoundaryTerms<<",compiled,"<<changed.report().compiledContraction
        <<",same_geometry_fresh_tau,"<<changed.report().compiledGeometryHit<<",literal_Green_row_RNE_and_work,1\n";
}
void run(){arithmetic();scene(false);scene(true);
    check(tapes==10&&loads==6&&powerChecks==6&&subnormalWitnesses>0,"compiled Maxwell finite scene coverage incomplete");
    std::cout<<"MAXWELL_COMPILED_CONTRACTION_PASS,tapes,"<<tapes<<",original_terms,"<<terms<<",actual_loads,"<<loads
        <<",native_force_rows,"<<forceRows<<",canonical_RHS_exact_power_checks,"<<powerChecks<<",guards,"<<guards
        <<",actual_subnormal_samples,"<<subnormalWitnesses<<",cap_dynamic_not_exercised,no_step_speed_claim\n";
}
} // namespace MaxwellCompiledContract


namespace GramGeometryStatisticsContract {
using Matrix=HydroEnergyViscosity::Matrix;
std::size_t scenes=0,energies=0,repeatedBounds=0,coefficientChecks=0,guardChecks=0;
HydroEnergyViscosity buildAndCheck(const HydroStrainOperator& strain,const std::vector<double>& mu,const std::vector<double>& scale) {
    // Original numerical product, independent of either geometry statistics or
    // the actual K assembly plan. Original entry order and explicit weights.
    const auto literal=ViscosityProductPlanContract::original(strain,mu);
    const auto energy=HydroEnergyViscosity::build(strain,mu);
    check(&energy.strain().matrix()==&strain.matrix(),"coarse statistics fixture lost the same actual immutable strain owner");
    ViscosityProductPlanContract::exactWeights(energy.rowWeights(),literal.weights);
    ViscosityProductPlanContract::exactMatrix(energy.matrix(),literal.matrix,coefficientChecks);
    // Complete exact Dyadic Lower-minus-Gram error, not a restatement of the
    // coarse counting/gamma formula and not a cache-pointer/hit probe.
    GramCoercivityContract::exactLowerOracle(energy,scale,false);
    const double maximum=*std::max_element(scale.begin(),scale.end());
    const auto first=GramCoercivityContract::coarseQuery(energy,maximum);
    const auto second=GramCoercivityContract::coarseQuery(energy,maximum);
    check(first&&second&&GramCoercivityContract::word(*first)==GramCoercivityContract::word(*second),
        "repeated same-source coarse bound changed its original result bits");
    ++repeatedBounds;++energies;return energy;
}
void scene(bool stressRz) {
    StressAmrLocalReuseContract::FloatingEnvironmentRestore restore;
    check(std::fesetround(FE_TONEAREST)==0,"coarse statistics fixture cannot set nearest mode");
    QuadtreeMesh mesh(input(4,stressRz?1:0));FaceSegmentGraph graph;graph.build(mesh);
    const auto layout=HydroDofLayout::build(mesh,NativeLeafDomain::build(mesh,graph,stressRz));
    const auto metrics=HydroPhysicalDivergence::build(mesh,HydroIntegralDivergence::build(mesh,layout));
    const auto strain=stressRz?HydroStrainOperator::buildStressCompatible(mesh,metrics)
        :HydroStrainOperator::build(mesh,metrics,HydroStrainOperator::volumeQuadrature(mesh,metrics));
    check(strain.isStressCompatible()==stressRz,"coarse statistics fixture lost raw/stress distinction");
    const Matrix originalStrain=strain.matrix();
    std::vector<double> scale(layout.velocityCount()),zero(strain.points().size()),mu(zero.size()),high(zero.size()),low(zero.size());
    for(std::size_t i=0;i<scale.size();++i)scale[i]=std::ldexp(1.,static_cast<int>(i%5)-2);
    for(std::size_t q=0;q<mu.size();++q){zero[q]=q%2?0.:-0.;mu[q]=q%7==0?0.:.25+.125*(q%4);high[q]=std::ldexp(mu[q],40);low[q]=std::ldexp(mu[q],-40);}
    // The first query has a genuine nonempty strain but an all-zero K. The
    // original build cannot have constructed its nonzero viscosity K plan.
    // Stress-compatible K has no such plan at all, including later weights.
    const auto inactive=buildAndCheck(strain,zero,scale);
    const double maximum=*std::max_element(scale.begin(),scale.end());
    const auto inactiveBound=GramCoercivityContract::coarseQuery(inactive,maximum);
    check(inactiveBound&&GramCoercivityContract::word(*inactiveBound)==GramCoercivityContract::word(0.),
        "actual plus/minus-zero row weights lost the exact zero coarse bound");
    const auto regular=buildAndCheck(strain,mu,scale),large=buildAndCheck(strain,high,scale);
    auto changedScale=scale;for(double& s:changedScale)s*=4.;
    const auto small=buildAndCheck(strain,low,changedScale);
    const auto regularBound=GramCoercivityContract::coarseQuery(regular,maximum),largeBound=GramCoercivityContract::coarseQuery(large,maximum);
    const auto smallBound=GramCoercivityContract::coarseQuery(small,*std::max_element(changedScale.begin(),changedScale.end()));
    check(regularBound&&largeBound&&smallBound&&*regularBound>0&&*largeBound>*regularBound&&*smallBound<*regularBound,
        "same-geometry coarse statistics froze prior zero/nonzero row weights or maximum scale");
    const auto changedBound=GramCoercivityContract::coarseQuery(regular,*std::max_element(changedScale.begin(),changedScale.end()));
    check(changedBound&&*changedBound>*regularBound,"coarse query reused an earlier maximum scale");
    GramCoercivityContract::scaleAndFpGuards(regular,scale); // Original warm current-FP/invalid-scale assertions.
    {
        StressAmrLocalReuseContract::FloatingEnvironmentRestore directed;
        check(std::fesetround(FE_UPWARD)==0&&std::fegetround()==FE_UPWARD,"warm source FP fixture cannot set actual upward mode");
        // The original factory must really succeed; do not swallow an exception
        // and call it an unsupported-source certificate test.
        const auto upward=HydroEnergyViscosity::build(strain,mu);
        check(std::fesetround(FE_TONEAREST)==0,"warm source FP fixture cannot restore nearest mode");
        check(!GramCoercivityContract::coarseQuery(upward,maximum)&&!GramCoercivityContract::query(upward,scale),
            "geometry metadata laundered an unsupported actual K-build environment");++guardChecks;
    }
    const auto afterFallback=GramCoercivityContract::coarseQuery(regular,maximum);
    check(afterFallback&&GramCoercivityContract::word(*afterFallback)==GramCoercivityContract::word(*regularBound),
        "unsupported query/source attempt poisoned the actual qualified geometry metadata");
    ViscosityProductPlanContract::exactMatrix(strain.matrix(),originalStrain,coefficientChecks);
    bool caught=false;auto invalid=mu;invalid.front()=std::numeric_limits<double>::quiet_NaN();
    try{HydroEnergyViscosity::build(strain,invalid);}catch(const std::invalid_argument& e){
        check(std::string(e.what())=="[HydroEnergyViscosity] negative/nonfinite viscosity","literal diagnostic overload changed the original failure text/type");caught=true;}
    check(caught,"literal diagnostic success optimization bypassed original false guard");++guardChecks;
    ++scenes;std::cout<<"GRAM_GEOMETRY_STATISTICS,"<<(stressRz?"RZ_stress":"XY_raw")<<",same_actual_E,1,zero_first,1,fresh_large_small_weights,1"
        <<",regular_bound,"<<*regularBound<<",large_bound,"<<*largeBound<<",small_new_scale_bound,"<<*smallBound<<",original_K_and_exact_Lower_oracles,PASS\n";
}
void run(){scene(false);scene(true);check(scenes==2&&energies==8&&repeatedBounds==8&&guardChecks==4
        &&GramCoercivityContract::mathCases==8&&GramCoercivityContract::exactRows==GramCoercivityContract::coarseRows,
        "coarse geometry statistics finite coverage incomplete");
    std::cout<<"GRAM_GEOMETRY_STATISTICS_PASS,scenes,"<<scenes<<",fresh_energy_forms,"<<energies<<",repeated_bound_bits,"<<repeatedBounds
        <<",original_K_coefficient_checks,"<<coefficientChecks<<",independent_exact_rows,"<<GramCoercivityContract::exactRows
        <<",source_FP_and_original_literal_guards,"<<guardChecks<<",current_FP_and_invalid_scale_original_guards,1,no_speed_claim\n";
}
} // namespace GramGeometryStatisticsContract
}
namespace Binary64LimbAccumulatorContract {
namespace FV=LAMP::Geometry::Detail::FiniteVolume;
using Port=FV::Binary64ExactProducts;
using Term=std::array<double,3>;
std::size_t tapes=0,terms=0,exactResults=0,rneResults=0,overflowResults=0,fpChecks=0,guards=0;
std::uint64_t bits(double value){std::uint64_t result=0;std::memcpy(&result,&value,sizeof result);return result;}
double fromBits(std::uint64_t value){double result=0;std::memcpy(&result,&value,sizeof result);return result;}
template<class F>void reject(F action,const char* literal) {
    bool caught=false;try{action();}catch(const std::invalid_argument& e){caught=std::string(e.what()).find(literal)!=std::string::npos;}
    check(caught,"fixed-limb exact products changed an original invalid-input guard");++guards;
}
void rounded(const FV::Dyadic& actual,const FV::Dyadic& expected) {
    std::optional<double> a,b;std::string ae,be;
    try{a=actual.roundedDouble();}catch(const std::invalid_argument& e){ae=e.what();}
    try{b=expected.roundedDouble();}catch(const std::invalid_argument& e){be=e.what();}
    check(a.has_value()==b.has_value(),"fixed-limb sum changed terminal overflow availability");
    if(a){check(bits(*a)==bits(*b),"fixed-limb sum changed terminal RNE/sign bits");++rneResults;}
    else{check(!ae.empty()&&ae==be,"fixed-limb sum changed original terminal overflow rejection");++overflowResults;}
}
void tape(const std::vector<Term>& source,int finalZeroSign=-1) {
    std::vector<Port::Factor> factors;factors.reserve(source.size());
    std::vector<Port::Value> values(source.size());std::vector<std::uint32_t> rows;rows.reserve(source.size());
    std::vector<Port::Factor> negatedFactors(source.size());
    FV::Dyadic original,absolute,negatedOriginal;
    for(std::size_t i=0;i<source.size();++i){const auto& t=source[i];const auto row=source.size()-1-i;
        factors.push_back(Port::compile(t[0],t[1]));values[row]=Port::decode(t[2]);rows.push_back(static_cast<std::uint32_t>(row));
        const auto term=FV::Dyadic::product({t[0],t[1],t[2]});original.add(term);absolute.add(term.absolute());
        negatedOriginal.add(FV::Dyadic::product({-t[0],t[1],t[2]}));++terms;}
    std::fenv_t before{},after{};check(std::fegetenv(&before)==0,"fixed-limb test cannot capture FP state");
#if defined(_MSC_VER) && defined(_M_X64)
    const auto csr=_mm_getcsr();
#endif
    const int savedError=errno;
    const auto actual=Port::sum(factors.data(),rows.data(),rows.size(),values.data(),values.size());
    for(std::size_t i=0;i<factors.size();++i)negatedFactors[i]=factors[i].negated();
    const auto negated=Port::sum(negatedFactors.data(),rows.data(),rows.size(),values.data(),values.size());
    check(errno==savedError&&std::fegetenv(&after)==0&&std::memcmp(&before,&after,sizeof before)==0,
        "fixed-limb integer accumulation changed complete FP environment/errno");
#if defined(_MSC_VER) && defined(_M_X64)
    check(csr==_mm_getcsr(),"fixed-limb integer accumulation changed raw MXCSR");
#endif
    ++fpChecks;
    check(actual.sum.compare(original)==0&&actual.absolute.compare(absolute)==0,
        "fixed-limb signed sum or original absolute-term sum differs from independent Dyadic products");exactResults+=2;
    rounded(actual.sum,original);rounded(actual.absolute,absolute);
    check(negated.sum.compare(negatedOriginal)==0&&negated.absolute.compare(absolute)==0,
        "negated compiled factors differ from independent signed Dyadic products or absolute terms");exactResults+=2;
    rounded(negated.sum,negatedOriginal);rounded(negated.absolute,absolute);
    if(finalZeroSign>=0){const double value=actual.sum.roundedDouble();
        check(value==0&&std::signbit(value)==(finalZeroSign!=0),"fixed-limb exact cancellation/tiny negative sign changed");}
    ++tapes;
}
Term power(int exponent) {
    const int a=std::max(-1074,std::min(1023,exponent));
    const int b=std::max(-1074,std::min(1023,exponent-a)),c=exponent-a-b;
    check(c>=-1074&&c<=1023,"power witness is not three actual finite binary64 factors");
    return {{std::ldexp(1.,a),std::ldexp(1.,b),std::ldexp(1.,c)}};
}
std::uint64_t randomWord() {
    static std::uint64_t state=UINT64_C(0x6a09e667f3bcc909);
    state^=state>>12;state^=state<<25;state^=state>>27;
    return state*UINT64_C(2685821657736338717);
}
double randomFinite() {
    const auto fraction=randomWord()&UINT64_C(0x000fffffffffffff);
    const auto exponent=(randomWord()%2047)<<52;
    return fromBits((randomWord()&UINT64_C(0x8000000000000000))|exponent|fraction);
}
void routingAndGuards() {
    const std::array<Port::Factor,4> factor{{Port::compile(2.,3.),Port::compile(-7.,.5),
        Port::compile(0.,std::numeric_limits<double>::max()),Port::compile(.25,-4.)}};
    const std::array<Port::Value,2> value{{Port::decode(.375),Port::decode(-.125)}};
    const std::array<std::uint32_t,4> rows{{1,0,1,1}};
    FV::Dyadic expected;expected.add(FV::Dyadic::product({2.,3.,-.125}));
    expected.add(FV::Dyadic::product({-7.,.5,.375}));expected.add(FV::Dyadic::product({.25,-4.,-.125}));
    const auto actual=Port::sum(factor.data(),rows.data(),rows.size(),value.data(),value.size());
    check(actual.sum.compare(expected)==0,"fixed-limb accumulator lost duplicate/permuted actual row routing");++exactResults;
    const std::uint32_t invalid=2,zeroInvalid=1;const auto zeroFactor=Port::compile(0.,1.);
    const auto zeroValue=Port::decode(-0.);
    reject([&]{Port::sum(factor.data(),&invalid,1,value.data(),value.size());},"compiled product row out of range");
    reject([&]{Port::sum(&zeroFactor,&zeroInvalid,1,&zeroValue,1);},"compiled product row out of range");
    reject([&]{Port::sum(factor.data(),&zeroInvalid,1,&zeroValue,1);},"compiled product row out of range");
    reject([&]{Port::sum(nullptr,nullptr,1,nullptr,0);},"missing compiled product terms");
    reject([&]{Port::sum(nullptr,nullptr,0,nullptr,1);},"missing compiled product values");
    reject([&]{Port::sum(factor.data(),rows.data(),1,nullptr,1);},"missing compiled product values");
    reject([&]{Port::decode(std::numeric_limits<double>::infinity());},"nonfinite compiled product input");
    reject([&]{Port::compile(1.,std::numeric_limits<double>::quiet_NaN());},"nonfinite compiled product input");
    const auto empty=Port::sum(nullptr,nullptr,0,nullptr,0);
    check(empty.sum.compare(FV::Dyadic{})==0&&empty.absolute.compare(FV::Dyadic{})==0,"count0 changed original exact zeros");
    // Pure integer summation supports directed controls; only the unchanged
    // terminal binary64 conversion requires nearest mode.
    StressAmrLocalReuseContract::FloatingEnvironmentRestore restore;
    const int originalError=errno;check(std::fesetround(FE_UPWARD)==0,"fixed-limb test cannot set upward mode");
    std::fenv_t before{},after{};check(std::fegetenv(&before)==0,"fixed-limb directed capture failed");
#if defined(_MSC_VER) && defined(_M_X64)
    const auto csr=_mm_getcsr();
#endif
    errno=EDOM;const auto directed=Port::sum(factor.data(),rows.data(),rows.size(),value.data(),value.size());
    check(errno==EDOM&&std::fegetenv(&after)==0&&std::memcmp(&before,&after,sizeof before)==0,
        "directed integer-only accumulator changed native control/status/errno");
#if defined(_MSC_VER) && defined(_M_X64)
    check(csr==_mm_getcsr(),"directed integer-only accumulator changed raw CSR");
#endif
    check(directed.sum.compare(expected)==0,"directed controls changed exact limb sum");++fpChecks;
    reject([&]{directed.sum.roundedDouble();},"binary64 round-to-nearest-even required");
    errno=originalError;restore.restore();
}
void run() {
    StressAmrLocalReuseContract::FloatingEnvironmentRestore restore;
    check(std::fesetround(FE_TONEAREST)==0,"fixed-limb test requires nearest for oracle/terminal conversion");
    const double tiny=std::numeric_limits<double>::denorm_min(),maximum=std::numeric_limits<double>::max();
    const double fullSignificand=fromBits(UINT64_C(0x3fffffffffffffff));
    tape({});tape({{0.,maximum,-1.},{-0.,maximum,tiny}},0);
    tape({{1.,1.,1.},{1.,1.,std::ldexp(1.,-53)}});
    tape({{1.,1.,1.},{3.,1.,std::ldexp(1.,-53)}});
    tape({{fullSignificand,fullSignificand,fullSignificand}});
    tape({{maximum,maximum,maximum},{-maximum,maximum,maximum},{tiny,tiny,tiny}},0);
    tape({{maximum,maximum,maximum},{-maximum,maximum,maximum},{-tiny,tiny,tiny}},1);
    tape({{maximum,maximum,maximum},{-maximum,maximum,maximum}},0);
    tape({{-std::ldexp(1.,-538),std::ldexp(1.,-538),.25}},1);
    // Every aligned bit is one, then the least bit is added again: carry
    // traverses the entire >6000-bit span, not just one 64-bit boundary.
    std::vector<Term> carry;carry.reserve(6293);
    for(int e=-3222;e<=3069;++e)carry.push_back(power(e));
    carry.push_back(power(-3222));tape(carry);
    for(auto& t:carry)t[0]=-t[0];tape(carry);
    std::vector<Term> many(65537,Term{{fullSignificand,fullSignificand,fullSignificand}});
    tape(many);for(auto& t:many)t[0]=-t[0];tape(many);
    // Maximum finite same-sign terms are valid exact sums even though their
    // terminal binary64 rounding must preserve the original overflow failure.
    tape(std::vector<Term>(4097,Term{{maximum,maximum,maximum}}));
    for(unsigned n=0;n<96;++n){std::vector<Term> random(1+n%64);
        for(auto& t:random)for(double& x:t)x=randomFinite();tape(random);}
    routingAndGuards();
    check(tapes==110&&guards==9&&fpChecks==111,"fixed-limb finite test coverage incomplete");
    std::cout<<"BINARY64_LIMB_ACCUMULATOR_PASS,tapes,"<<tapes<<",terms,"<<terms<<",exact_results,"<<exactResults
        <<",RNE_fields,"<<rneResults<<",original_overflow_results,"<<overflowResults<<",FP_and_errno,"<<fpChecks
        <<",original_guards,"<<guards<<",generic_count_no_new_limit,1,full_span_carry,1,absolute_original_terms,1";
#if defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
    std::cout<<",implementation,MSVC_x64_fixed_limbs\n";
#else
    std::cout<<",implementation,original_generic_fallback\n";
#endif
}
} // namespace Binary64LimbAccumulatorContract


namespace ViscosityProductParallelContract {
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
using Matrix=HydroEnergyViscosity::Matrix;
std::size_t scenes=0,cases=0,coefficients=0,fpComparisons=0,nestedFallbacks=0,subnormalScenes=0,
    coldOwners=0,coldLocalHits=0,coldSymbolicZeroSlots=0;
struct Floating {
    std::fenv_t native{};unsigned csr=_mm_getcsr();int error=errno;
    Floating(){check(std::fegetenv(&native)==0,"parallel K oracle cannot capture caller FP");}
    void restore()const{check(std::fesetenv(&native)==0,"parallel K oracle cannot restore caller FP");_mm_setcsr(csr);errno=error;}
};
struct RestoreFloating {Floating saved;~RestoreFloating()noexcept{if(std::fesetenv(&saved.native)!=0)std::terminate();_mm_setcsr(saved.csr);errno=saved.error;}};
struct RestoreRuntime {Floating saved;int threads=omp_get_max_threads(),dynamic=omp_get_dynamic();
    ~RestoreRuntime()noexcept{omp_set_num_threads(threads);omp_set_dynamic(dynamic);if(std::fesetenv(&saved.native)!=0)std::terminate();_mm_setcsr(saved.csr);errno=saved.error;}};
Floating cleanInput() {
    Floating current;auto native=current.native;native._Fe_stat=0;
    check(std::fesetenv(&native)==0,"parallel K oracle cannot clear native sticky flags");
    _mm_setcsr(current.csr&~0x3fu);errno=E2BIG;return Floating{};
}
void sameFloating(const Floating& a,const Floating& b) {
    check(std::memcmp(&a.native,&b.native,sizeof a.native)==0&&a.csr==b.csr&&a.error==b.error,
        "parallel K differs from original serial/literal full native FP, raw MXCSR or errno");++fpComparisons;
}
struct Pair {HydroEnergyViscosity serial;Floating after,input;};
Pair compare(const HydroStrainOperator& raw,const std::vector<double>& mu,int team,const char* name,bool needSubnormal) {
    const Floating input=cleanInput();input.restore();
    const auto literal=ViscosityProductPlanContract::original(raw,mu);const Floating literalAfter;
    omp_set_num_threads(1);input.restore();const auto serial=HydroEnergyViscosity::build(raw,mu);const Floating serialAfter;
    omp_set_num_threads(team);input.restore();const auto parallel=HydroEnergyViscosity::build(raw,mu);const Floating parallelAfter;
    input.restore();const auto repeated=HydroEnergyViscosity::build(raw,mu);const Floating repeatedAfter;
    sameFloating(literalAfter,serialAfter);sameFloating(serialAfter,parallelAfter);sameFloating(parallelAfter,repeatedAfter);
    for(const auto* energy:{&serial,&parallel,&repeated}){
        check(&energy->strain().matrix()==&raw.matrix(),"parallel product lost the actual shared raw strain owner");
        ViscosityProductPlanContract::exactWeights(energy->rowWeights(),literal.weights);
        ViscosityProductPlanContract::exactMatrix(energy->matrix(),literal.matrix,coefficients);
    }
    if(needSubnormal){std::size_t w=0,k=0;for(double value:parallel.rowWeights())w+=std::fpclassify(value)==FP_SUBNORMAL;
        for(Eigen::Index i=0;i<parallel.matrix().nonZeros();++i)k+=std::fpclassify(parallel.matrix().valuePtr()[i])==FP_SUBNORMAL;
        check(w>0&&k>0,"parallel K tiny scene lacks actual represented subnormal rowWeights and K");++subnormalScenes;}
    ++cases;std::cout<<"VISCOSITY_PRODUCT_PARALLEL,"<<name<<",team,"<<team<<",nnz,"<<parallel.matrix().nonZeros()
        <<",original_Eigen_and_OMP1_matrix_bits,1,clean_to_final_native_CSR_errno_exact,1\n";
    return {serial,serialAfter,input};
}
void cold(const HydroStrainOperator& raw,const std::vector<double>& mu,int team,std::size_t zeroStructure) {
    // An actual new E owner reaches the FIRST numeric call at this request.
    // Both calls use the independent old weighted-Eigen expression, never a
    // newly prepared production object as a literal arithmetic reference.
    omp_set_num_threads(team);const Floating entry=cleanInput();entry.restore();
    const auto literal=ViscosityProductPlanContract::original(raw,mu);const Floating literalAfter;
    entry.restore();const auto first=HydroEnergyViscosity::build(raw,mu);const Floating firstAfter;
    sameFloating(literalAfter,firstAfter);
    ViscosityProductPlanContract::exactWeights(first.rowWeights(),literal.weights);
    ViscosityProductPlanContract::exactMatrix(first.matrix(),literal.matrix,coefficients);
    const auto& f=first.productReport();
    check(f.requested&&f.eligible&&f.currentFpQualified&&!f.sourceAvailable&&!f.sourceFpReceipt
        &&f.fullAssembly&&!f.localUsed&&f.memoPublished&&f.fallbackReason==0
        &&f.recomputedRows==static_cast<std::size_t>(raw.matrix().cols())&&f.reusedRows==0,
        "first qualified numeric call did not compute and publish real complete row receipts");
    check(f.symbolicStoredValues==literal.symbolicStoredValues&&f.symbolicStoredValues==zeroStructure,
        "cold symbolic CSR omitted independent original Eigen structural zero slots");
    check(f.peakOptionalBytes<=f.byteLimit&&f.retainedBytes<=f.byteLimit,"cold product exceeded its actual optional byte limit");
    auto changed=mu;changed[changed.size()/2]+=.125;entry.restore();
    const auto changedLiteral=ViscosityProductPlanContract::original(raw,changed);const Floating changedLiteralAfter;
    entry.restore();const auto second=HydroEnergyViscosity::build(raw,changed);const Floating secondAfter;
    sameFloating(changedLiteralAfter,secondAfter);
    ViscosityProductPlanContract::exactWeights(second.rowWeights(),changedLiteral.weights);
    ViscosityProductPlanContract::exactMatrix(second.matrix(),changedLiteral.matrix,coefficients);
    const auto& c=second.productReport();
    check(c.eligible&&c.sourceAvailable&&c.sourceFpReceipt&&c.memoPublished&&c.localUsed&&!c.fullAssembly
        &&c.changedWeightRows==static_cast<std::size_t>(raw.components())&&c.recomputedRows>0&&c.reusedRows>0,
        "second cold-owner call did not use its actual first source and changed local weights");
    check(c.symbolicStoredValues==zeroStructure,"changed local weights changed the retained complete symbolic CSR");
    ++coldOwners;++coldLocalHits;coldSymbolicZeroSlots+=zeroStructure;
    std::cout<<"VISCOSITY_PRODUCT_COLD,requested_team,"<<team<<",source_first,0,first_complete_row_memo,1,second_actual_local_hit,1"
        <<",unpruned_zero_weight_Eigen_slots,"<<zeroStructure<<",matrix_CSR_bits_and_full_FP_errno,1\n";
}
void nested(const HydroStrainOperator& raw,const std::vector<double>& mu,const Pair& reference) {
    std::optional<HydroEnergyViscosity> result;std::optional<Floating> after;std::exception_ptr failure;bool active=false;
#pragma omp parallel num_threads(2)
    {
#pragma omp single
        {
            try{RestoreFloating restore;reference.input.restore();active=omp_in_parallel()!=0;
                check(active,"parallel K nested serial-fallback fixture did not enter an actual active region");
                result.emplace(HydroEnergyViscosity::build(raw,mu));after.emplace();
            }catch(...){failure=std::current_exception();}
        }
    }
    if(failure)std::rethrow_exception(failure);
    check(active&&result&&after,"parallel K nested fallback did not actually produce an original result");
    sameFloating(reference.after,*after);
    ViscosityProductPlanContract::exactWeights(result->rowWeights(),reference.serial.rowWeights());
    ViscosityProductPlanContract::exactMatrix(result->matrix(),reference.serial.matrix(),coefficients);++nestedFallbacks;
}
void scene(bool rz,int team) {
    omp_set_num_threads(1);QuadtreeMesh mesh(input(8,rz?1:0));refine(mesh,1);FaceSegmentGraph graph;graph.build(mesh);
    const auto domain=NativeLeafDomain::build(mesh,graph,rz);const auto layout=HydroDofLayout::build(mesh,domain);
    check(!domain.uniformLeaves(),"parallel product fixture lost actual XY/RZ 2:1 AMR");
    const auto metrics=HydroPhysicalDivergence::build(mesh,HydroIntegralDivergence::build(mesh,layout));
    const auto raw=HydroStrainOperator::build(mesh,metrics,HydroStrainOperator::volumeQuadrature(mesh,metrics));const Matrix saved=raw.matrix();
    std::vector<double> initial(raw.points().size()),changed(initial.size()),zero(initial.size()),tiny(initial.size(),1e-319);
    for(std::size_t q=0;q<initial.size();++q){initial[q]=.5+static_cast<double>(q%7)/8;
        changed[q]=q%3==0?(q%2?0.:-0.):.125+static_cast<double>(q%11)/16;zero[q]=q%2?0.:-0.;}
    // The same existing fixture independently witnesses the complete old
    // Eigen zero-weight pattern; zero slots must survive symbolic preparation.
    const auto zeroStructure=ViscosityProductPlanContract::original(raw,zero);
    check(zeroStructure.symbolicStoredValues>0&&zeroStructure.symbolicZeroValues==zeroStructure.symbolicStoredValues
        &&zeroStructure.matrix.nonZeros()==0,"cold fixture lacks an independent unpruned structural-zero witness");
    for(int coldTeam:{1,2,8}){
        const auto coldRaw=coldTeam==1?raw:HydroStrainOperator::build(mesh,metrics,raw.points());
        check(coldTeam==1||&coldRaw.matrix()!=&raw.matrix(),"cold parallel fixture reused an old E owner");
        cold(coldRaw,initial,coldTeam,zeroStructure.symbolicStoredValues);
    }
    // Restore the original initial coefficients; retain every prior mixed,
    // inactive, subnormal and nested-fallback comparison below.
    omp_set_num_threads(1);const auto prime=HydroEnergyViscosity::build(raw,initial);
    const auto literalPrime=ViscosityProductPlanContract::original(raw,initial);
    ViscosityProductPlanContract::exactWeights(prime.rowWeights(),literalPrime.weights);
    ViscosityProductPlanContract::exactMatrix(prime.matrix(),literalPrime.matrix,coefficients);
    const auto normal=compare(raw,changed,team,rz?"RZ fresh mixed weights":"XY fresh mixed weights",false);
    nested(raw,changed,normal);
    compare(raw,zero,team,"fresh plus/minus-zero inactive form",false);
    compare(raw,tiny,team,"fresh actual subnormal W and K",true);
    ViscosityProductPlanContract::exactMatrix(raw.matrix(),saved,coefficients);++scenes;
}
void run() {
    RestoreRuntime restore;omp_set_dynamic(0);omp_set_num_threads(3);
    check(std::fesetround(FE_TONEAREST)==0,"parallel product qualification cannot set nearest mode");
    const auto environment=ViscosityProductPlanContract::environment();
    check(environment.eligible&&!omp_in_parallel(),"parallel product selector cannot certify the actual qualified native path");
    int actual=0;
#pragma omp parallel num_threads(3)
    {
#pragma omp atomic update
        ++actual;
    }
    check(actual>=2,"parallel product selector has no actual worker team; it cannot qualify parallel execution");
    scene(false,3);scene(true,3);
    check(scenes==2&&cases==6&&fpComparisons==32&&nestedFallbacks==2&&subnormalScenes==2
        &&coldOwners==6&&coldLocalHits==6&&coldSymbolicZeroSlots>0,"parallel product finite coverage incomplete");
    std::cout<<"VISCOSITY_PRODUCT_PARALLEL_PASS,actual_probe_team,"<<actual<<",XY_RZ_AMR_scenes,"<<scenes<<",fresh_weight_cases,"<<cases
        <<",original_coefficient_bits,"<<coefficients<<",full_FP_errno_comparisons,"<<fpComparisons<<",active_nested_serial_fallbacks,"<<nestedFallbacks
        <<",actual_subnormal_W_K_scenes,"<<subnormalScenes<<",cold_1_2_8_new_E_owners,"<<coldOwners
        <<",second_call_local_hits,"<<coldLocalHits<<",independent_symbolic_zero_slots,"<<coldSymbolicZeroSlots
        <<",64MiB_cap_and_worker_fault_dynamic_not_exercised,no_speed_claim\n";
}
#else
void run(){throw std::runtime_error("parallel product selector requires its actual qualified native MSVC x64 OpenMP execution contract");}
#endif
} // namespace ViscosityProductParallelContract

namespace ViscosityProductLocalMemoContract {
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
using Matrix=HydroEnergyViscosity::Matrix;
using Floating=ViscosityProductParallelContract::Floating;
std::size_t scenes=0,cases=0,coefficients=0,fpChecks=0,localCases=0,controlFallbacks=0,coldOwners=0,guards=0;
Floating inputState(bool sticky) {
    Floating old;auto native=old.native;native._Fe_stat=0;
    check(std::fesetenv(&native)==0,"local K oracle cannot set entry status");
    _mm_setcsr(old.csr&~0x3fu);
    if(sticky){
        // Public FE_* masks are NOT a serializable UCRT _Fe_stat encoding.
        // Raise real masked exceptions through its public API, then capture
        // the native representation; raw CSR additionally retains denormal.
        check(std::feraiseexcept(FE_ALL_EXCEPT)==0,"local K oracle cannot raise actual entry exceptions");
        check(std::fetestexcept(FE_ALL_EXCEPT)==FE_ALL_EXCEPT,"local K oracle did not actually raise entry exceptions");
        _mm_setcsr(_mm_getcsr()|0x02u);
    }
    errno=sticky?EDOM:E2BIG;return Floating{};
}
void sameFp(const Floating& a,const Floating& b) {
    check(std::memcmp(&a.native,&b.native,sizeof a.native)==0&&a.csr==b.csr&&a.error==b.error,
        "local K did not reproduce literal full native/raw-CSR sticky and errno");++fpChecks;
}
HydroEnergyViscosity compare(const HydroStrainOperator& raw,const std::vector<double>& mu,bool sticky,const char* label) {
    const Floating entry=inputState(sticky);entry.restore();
    const auto literal=ViscosityProductPlanContract::original(raw,mu);const Floating literalAfter;
    entry.restore();const auto actual=HydroEnergyViscosity::build(raw,mu);const Floating actualAfter;
    sameFp(literalAfter,actualAfter);
    check(&actual.strain().matrix()==&raw.matrix(),"local K factory lost actual immutable E owner");
    ViscosityProductPlanContract::exactWeights(actual.rowWeights(),literal.weights);
    ViscosityProductPlanContract::exactMatrix(actual.matrix(),literal.matrix,coefficients);
    const auto& r=actual.productReport();
    check(r.peakOptionalBytes<=r.byteLimit&&r.retainedBytes<=r.byteLimit,"local K optional stocks exceed owned byte limit");
    if(r.eligible){check(r.recomputedRows+r.reusedRows==r.outputRows&&r.outputRows==static_cast<std::size_t>(raw.matrix().cols()),
        "local K real output-row accounting is incomplete");}
    ++cases;std::cout<<"VISCOSITY_PRODUCT_LOCAL,"<<label<<",eligible,"<<r.eligible<<",source,"<<r.sourceAvailable
        <<",weights_changed,"<<r.changedWeightRows<<",recomputed,"<<r.recomputedRows<<",reused,"<<r.reusedRows
        <<",matrix_bits_and_full_FP_errno,1\n";
    return actual;
}
void allClean(const HydroEnergyViscosity& a,const char* why) {
    const auto& r=a.productReport();check(r.requested&&r.eligible&&r.sourceAvailable&&r.sourceFpReceipt&&r.currentFpQualified
        &&r.localUsed&&!r.fullAssembly&&r.memoPublished&&r.changedWeightRows==0&&r.recomputedRows==0&&r.reusedRows==r.outputRows,why);
    ++localCases;
}
void unsupported(const HydroStrainOperator& raw,const std::vector<double>& mu) {
    ViscosityProductParallelContract::RestoreFloating restore;
    check(std::fesetround(FE_UPWARD)==0,"local K directed-control negative could not set UPWARD");
    const auto directed=compare(raw,mu,true,"actual current UPWARD complete old path");const auto& r=directed.productReport();
    check(r.requested&&!r.currentFpQualified&&!r.eligible&&!r.localUsed&&!r.memoPublished&&r.fallbackReason==2,
        "local K reused source under directed current controls");++controlFallbacks;
    restore.saved.restore();_mm_setcsr(_mm_getcsr()|(1u<<15)|(1u<<6));
    const auto nongradual=compare(raw,mu,false,"actual current FTZ/DAZ complete old path");const auto& f=nongradual.productReport();
    check(f.requested&&!f.currentFpQualified&&!f.eligible&&!f.localUsed&&!f.memoPublished&&f.fallbackReason==2,
        "local K reused source under FTZ/DAZ current controls");++controlFallbacks;
}
void sourceControl(const QuadtreeMesh& mesh,const HydroPhysicalDivergence& metrics,const std::vector<double>& mu) {
    // This is an actually new E owner. Its first two original Energy factories
    // must succeed under UPWARD; an exception is a failed test, never a skip.
    const auto foreign=HydroStrainOperator::build(mesh,metrics,HydroStrainOperator::volumeQuadrature(mesh,metrics));
    {ViscosityProductParallelContract::RestoreFloating restore;
        check(std::fesetround(FE_UPWARD)==0,"local K source-control negative could not set UPWARD");
        const auto first=compare(foreign,mu,false,"new source UPWARD Eigen");
        check(!first.productReport().memoPublished,"directed first Eigen invented row receipt");
        const auto second=compare(foreign,mu,true,"new source UPWARD Plan fallback");
        check(!second.productReport().eligible&&!second.productReport().memoPublished,"UPWARD source seeded a qualified local memo");
        ++controlFallbacks;
    }
    const auto nearest=compare(foreign,mu,false,"restored RNE needs full actual row receipt");const auto& r=nearest.productReport();
    check(r.eligible&&!r.sourceAvailable&&!r.sourceFpReceipt&&r.fullAssembly&&r.reusedRows==0&&r.memoPublished,
        "current nearest controls relabelled an UPWARD source as qualified");
    allClean(compare(foreign,mu,true,"RNE new source next actual hit"),"new RNE source was not independently published");++coldOwners;
}
void scene(bool rz) {
    omp_set_num_threads(1);QuadtreeMesh mesh(input(8,rz?1:0));refine(mesh,1);FaceSegmentGraph graph;graph.build(mesh);
    const auto domain=NativeLeafDomain::build(mesh,graph,rz);const auto layout=HydroDofLayout::build(mesh,domain);
    check(!domain.uniformLeaves(),"local K fixture lacks actual XY/RZ 2:1 AMR");
    const auto metrics=HydroPhysicalDivergence::build(mesh,HydroIntegralDivergence::build(mesh,layout));
    const auto points=HydroStrainOperator::volumeQuadrature(mesh,metrics);
    const auto raw=HydroStrainOperator::build(mesh,metrics,points);const Matrix saved=raw.matrix();
    std::vector<double> mu(points.size());for(std::size_t q=0;q<mu.size();++q)mu[q]=.375+static_cast<double>(q%7)/16;
    check(mu.size()>8,"local K fixture too small for distinct local weights");mu[1]=0.;mu[2]=-0.;
    omp_set_num_threads(3);
    const auto first=compare(raw,mu,false,"first numeric call records actual complete row flags");const auto& f=first.productReport();
    check(f.requested&&f.eligible&&f.currentFpQualified&&!f.sourceAvailable&&!f.sourceFpReceipt&&f.memoPublished
        &&f.fullAssembly&&!f.localUsed&&f.reusedRows==0&&f.recomputedRows==f.outputRows&&f.fallbackReason==0,
        "first numeric row kernel failed to publish its real complete source receipt");
    auto seeded=mu;seeded[seeded.size()/2]+=.125;
    const auto seed=compare(raw,seeded,false,"second call uses first numeric source with local changed weight");const auto& s=seed.productReport();
    check(s.eligible&&s.sourceAvailable&&s.sourceFpReceipt&&s.localUsed&&!s.fullAssembly&&s.recomputedRows>0&&s.reusedRows>0
        &&s.changedWeightRows==static_cast<std::size_t>(raw.components())&&s.memoPublished,
        "second call did not reuse the real cold source and recompute actual local support");++localCases;
    const auto clean=compare(raw,seeded,true,"same fresh W different entry sticky");
    allClean(clean,"same weights with different sticky failed complete clean-row reuse");
    check(&clean.matrix()==&seed.matrix(),"all-clean K did not share actual final pruned matrix owner");
    auto local=seeded;local[local.size()/3]+=.125;
    const auto changed=compare(raw,local,false,"one actual local quadrature viscosity change");const auto& c=changed.productReport();
    check(c.eligible&&c.sourceFpReceipt&&c.changedWeightRows==static_cast<std::size_t>(raw.components())
        &&c.recomputedRows>0&&c.reusedRows>0&&c.localUsed&&!c.fullAssembly,
        "one local fresh weight did not actually recompute its bounded original row support and reuse other rows");++localCases;
    auto zeros=local;zeros[1]=-0.;zeros[2]=0.;
    const auto signedZero=compare(raw,zeros,true,"fresh plus/minus zero weight sign changes");const auto& z=signedZero.productReport();
    check(z.eligible&&z.changedWeightRows==2*static_cast<std::size_t>(raw.components())&&z.recomputedRows>0&&z.reusedRows>0,
        "signed-zero fresh weight bits were incorrectly treated as identical");++localCases;
    const std::vector<double> inactive(mu.size(),-0.);
    const auto zero=compare(raw,inactive,false,"original exact inactive form");
    check(zero.matrix().nonZeros()==0&&!zero.productReport().memoPublished,"inactive form rewrote a qualified numeric slot");
    const auto afterInactive=compare(raw,zeros,false,"active source survives inactive form");
    allClean(afterInactive,"inactive form contaminated previous actual active source");
    check(&afterInactive.matrix()==&signedZero.matrix(),"inactive zero K replaced the real active immutable stock");
    unsupported(raw,zeros);
    allClean(compare(raw,zeros,true,"qualified source survives unsupported current controls"),"unsupported current controls poisoned the retained qualified source");
    std::vector<double> tiny(mu.size(),1e-319);
    const auto small=compare(raw,tiny,false,"fresh actual subnormal weighted product");
    std::size_t sub=0,subK=0;for(double w:small.rowWeights())sub+=std::fpclassify(w)==FP_SUBNORMAL;
    for(Eigen::Index k=0;k<small.matrix().nonZeros();++k)subK+=std::fpclassify(small.matrix().valuePtr()[k])==FP_SUBNORMAL;
    check(sub>0&&subK>0&&small.productReport().eligible,"local K small scene lacks actual qualified subnormal weights and K");
    const auto smallClean=compare(raw,tiny,true,"subnormal clean rows replay all native and CSR raises");
    allClean(smallClean,"subnormal source did not replay actual per-row underflow/denormal/inexact flags");
    check(&smallClean.matrix()==&small.matrix(),"subnormal all-clean K copied instead of sharing final immutable stock");
    const std::vector<double> large(mu.size(),std::ldexp(1.,400));
    compare(raw,large,false,"finite large-exponent fresh inputs");
    allClean(compare(raw,large,true,"finite large-exponent repeated inputs"),"large finite source failed exact clean-row reuse");
    // Same mesh/points/layout with a new actual strain owner is a legal cold
    // factory, never a source borrowed through content equality or mesh hash.
    const auto fresh=HydroStrainOperator::build(mesh,metrics,points);
    const auto independent=compare(fresh,mu,false,"equal geometry new actual E owner");
    const auto& freshReport=independent.productReport();
    check(freshReport.eligible&&!freshReport.sourceAvailable&&!freshReport.sourceFpReceipt&&!freshReport.localUsed
        &&freshReport.fullAssembly&&freshReport.memoPublished&&freshReport.fallbackReason==0,
        "numeric K crossed actual immutable E owners or failed independent cold preparation");++coldOwners;
    sourceControl(mesh,metrics,mu);
    for(int kind=0;kind<3;++kind){auto invalid=mu;if(kind==0)invalid.pop_back();if(kind==1)invalid[0]=-1.;if(kind==2)invalid[0]=std::numeric_limits<double>::infinity();
        bool caught=false;try{HydroEnergyViscosity::build(raw,invalid);}catch(const std::invalid_argument& e){
            const std::string text=e.what();caught=text.find(kind==0?"viscosity dimensions do not match quadrature":"negative/nonfinite viscosity")!=std::string::npos;}
        check(caught,"local K changed original viscosity input rejection");++guards;}
    ViscosityProductPlanContract::exactMatrix(raw.matrix(),saved,coefficients);++scenes;
}
void run() {
    ViscosityProductParallelContract::RestoreRuntime restore;omp_set_dynamic(0);omp_set_num_threads(3);
    check(std::fesetround(FE_TONEAREST)==0,"local K selector cannot set nearest controls");
    _mm_setcsr(_mm_getcsr()&~((1u<<15)|(1u<<6)));
    check(ViscosityProductPlanContract::environment().eligible,"local K selector lacks actual known masked binary64 implementation");
    scene(false);scene(true);
    check(scenes==2&&cases==38&&fpChecks==cases&&localCases==18&&controlFallbacks==6&&coldOwners==4&&guards==6,
        "local K finite actual-work coverage incomplete");
    std::cout<<"VISCOSITY_PRODUCT_LOCAL_MEMO_PASS,XY_RZ_AMR_scenes,"<<scenes<<",literal_matrix_and_full_FP_errno_pairs,"<<cases
        <<",coefficient_bits,"<<coefficients<<",actual_local_or_clean_hits,"<<localCases<<",actual_control_fallbacks,"<<controlFallbacks
        <<",independent_cold_owners,"<<coldOwners<<",original_input_guards,"<<guards
        <<",128MiB_cap_concurrent_busy_and_worker_fault_dynamic_not_exercised,no_speed_claim\n";
}
#else
void run(){throw std::runtime_error("local K memo selector requires qualified native MSVC x64 OpenMP");}
#endif
} // namespace ViscosityProductLocalMemoContract


namespace NativeSmallProductContract
{
namespace FV=LAMP::Geometry::Detail::FiniteVolume;
using Integer=boost::multiprecision::cpp_int;
std::size_t cases=0,exactChains=0,wideChecks=0,fpPairs=0,terminalPairs=0,rneChecks=0,terminalOverflows=0,
    originalGuards=0,zeroShortCircuits=0,controlCases=0,fullWidthProducts=0;
std::uint64_t bits(double x){std::uint64_t n=0;std::memcpy(&n,&x,sizeof n);return n;}
double fromBits(std::uint64_t n){double x=0;std::memcpy(&x,&n,sizeof x);return x;}
struct Floating
{
    std::fenv_t native{};int error=errno;
#if defined(_MSC_VER) && defined(_M_X64)
    unsigned csr=_mm_getcsr();
#endif
    Floating(){check(std::fegetenv(&native)==0,"small product oracle cannot capture actual FP environment");}
    void restore()const{check(std::fesetenv(&native)==0,"small product oracle cannot restore actual FP environment");
#if defined(_MSC_VER) && defined(_M_X64)
        _mm_setcsr(csr);
#endif
        errno=error;
    }
};
struct Restore
{
    Floating saved;~Restore()noexcept{std::fesetenv(&saved.native);
#if defined(_MSC_VER) && defined(_M_X64)
        _mm_setcsr(saved.csr);
#endif
        errno=saved.error;}
};
void sameFp(const Floating& a,const Floating& b)
{
    check(std::memcmp(&a.native,&b.native,sizeof a.native)==0&&a.error==b.error,
        "native small product changed literal full FP controls/status/errno");
#if defined(_MSC_VER) && defined(_M_X64)
    check(a.csr==b.csr,"native small product changed literal complete raw MXCSR");
#endif
    ++fpPairs;
}
// The unmodified constructor is invoked ONCE per consumed factor by multiplied.
// Stop immediately after first exact zero; do not let multiplied inspect a
// later NaN merely because the running value was already zero.
FV::Dyadic literal(std::initializer_list<long double> factors)
{
    FV::Dyadic result(1);
    for(const auto f:factors){result=result.multiplied(f);if(result.compare(FV::Dyadic{})==0)return {};}
    return result;
}
template<class T>struct Outcome {std::optional<T> value;std::string error;};
template<class T,class F>Outcome<T> invoke(F action)
{Outcome<T> result;try{result.value.emplace(action());}catch(const std::invalid_argument& e){result.error=e.what();}return result;}
template<class T,class A,class B>Outcome<T> paired(A actual,B reference)
{
    const Floating entry;auto a=invoke<T>(actual);const Floating afterA;entry.restore();auto b=invoke<T>(reference);const Floating afterB;
    sameFp(afterA,afterB);check(a.value.has_value()==b.value.has_value()&&a.error==b.error,
        "native small product changed original availability/first error string");
    if constexpr(std::is_same_v<T,FV::Dyadic>){if(a.value){check(a.value->compare(*b.value)==0,"native product differs from original constructor.multiplied exact chain");++exactChains;}}
    else if(a.value){if constexpr(std::is_same_v<T,double>)check(bits(*a.value)==bits(*b.value),"native product changed original RNE binary64/sign bits");
        else check(*a.value==*b.value&&std::signbit(*a.value)==std::signbit(*b.value),"native product changed original native rounding/sign");}
    return a;
}
struct Wide {Integer magnitude=1;int exponent=0;bool negative=false,valid=true;};
Wide independent(std::initializer_list<long double> factors)
{
    Wide r;
    for(const auto f:factors){
        // Every test finite scalar originates in a representable binary64.
        // IEEE decoding and multiplication are independent of all FV integers,
        // factors, private limbs, product(), or compiled column-sum ports.
        const auto word=bits(static_cast<double>(f));const auto field=static_cast<unsigned>((word>>52)&0x7ffu);
        const auto fraction=word&UINT64_C(0x000fffffffffffff);
        if(field==0x7ffu){r.valid=false;return r;}
        const std::uint64_t integer=field?(fraction|UINT64_C(0x0010000000000000)):fraction;
        if(integer==0){r.magnitude=0;r.exponent=0;r.negative=false;return r;}
        r.magnitude*=integer;r.exponent+=field?static_cast<int>(field)-1075:-1074;r.negative=r.negative!=((word>>63)!=0);
    }
    return r;
}
FV::Dyadic power(int exponent)
{
    FV::Dyadic result(1);
    while(exponent){const int step=std::max(-512,std::min(512,exponent));result=result.multiplied(std::ldexp(1.L,step));exponent-=step;}
    return result;
}
FV::Dyadic expanded(const Wide& w)
{
    FV::Dyadic result;if(w.magnitude==0)return result;
    const auto count=boost::multiprecision::msb(w.magnitude)+1;
    for(unsigned bit=0;bit<count;++bit)if(boost::multiprecision::bit_test(w.magnitude,bit))result.add(power(w.exponent+static_cast<int>(bit)),w.negative?-1:1);
    return result;
}
std::optional<std::uint64_t> ieeeRne(const Wide& w)
{
    const auto sign=w.negative?UINT64_C(0x8000000000000000):UINT64_C(0);
    if(w.magnitude==0)return UINT64_C(0);
    const int top=static_cast<int>(boost::multiprecision::msb(w.magnitude))+w.exponent;
    if(top>=1024)return {};
    const int quantum=std::max(top-52,-1074),shift=quantum-w.exponent;Integer rounded;
    if(shift>0){rounded=w.magnitude>>shift;const Integer rest=w.magnitude-(rounded<<shift),half=Integer(1)<<(shift-1);
        if(rest>half||(rest==half&&boost::multiprecision::bit_test(rounded,0)))++rounded;}
    else rounded=w.magnitude<<(-shift);
    if(rounded==0)return sign;
    const int last=static_cast<int>(boost::multiprecision::msb(rounded)),roundedTop=last+quantum;
    if(roundedTop>=1024)return {};
    if(roundedTop>=-1022){Integer significand;if(last>52)significand=rounded>>(last-52);else significand=rounded<<(52-last);
        const auto fraction=(significand-(Integer(1)<<52)).convert_to<std::uint64_t>();
        return sign|(static_cast<std::uint64_t>(roundedTop+1023)<<52)|fraction;}
    const Integer fraction=rounded<<(quantum+1074);return sign|fraction.convert_to<std::uint64_t>();
}
bool gradualControls()
{
#if defined(_MSC_VER) && defined(_M_X64)
    return (_mm_getcsr()&((1u<<15)|(1u<<6)))==0;
#else
    return true;
#endif
}
void tape(std::initializer_list<long double> factors,int expectedFailure=-1,bool expectPositiveZero=false)
{
    // Independent oracle work is outside both measured entries. Invalid inputs
    // are still tested by the real constructor prefix; no exception is skipped.
    const bool wideEligible=gradualControls();const auto wide=wideEligible?independent(factors):Wide{};
    std::optional<FV::Dyadic> expected;if(wideEligible&&wide.valid)expected=expanded(wide);
    check(std::feclearexcept(FE_ALL_EXCEPT)==0,"small product oracle cannot clear real entry flags");
#if defined(_MSC_VER) && defined(_M_X64)
    _mm_setcsr(_mm_getcsr()&~0x3fu);
#endif
    if(cases%2)check(std::feraiseexcept(FE_INEXACT)==0,"small product oracle cannot prime real sticky inexact");errno=E2BIG;
    const auto actual=paired<FV::Dyadic>([&]{return FV::Dyadic::product(factors);},[&]{return literal(factors);});
    if(expectedFailure>=0)check(actual.value.has_value()==(expectedFailure==0),"small product first-zero/first-invalid order changed");
    if(!actual.value){check(actual.error=="[FiniteVolume] nonfinite dyadic input","small product original nonfinite input message changed");++originalGuards;++cases;return;}
    if(expected){check(actual.value->compare(*expected)==0,"small native product differs from independent unbounded IEEE integer bit expansion");++wideChecks;
        if(wide.magnitude!=0&&boost::multiprecision::msb(wide.magnitude)-boost::multiprecision::lsb(wide.magnitude)>=192)++fullWidthProducts;}
    if(expectPositiveZero){check(actual.value->compare(FV::Dyadic{})==0,"first zero did not short-circuit exact product");++zeroShortCircuits;}
    const auto old=literal(factors);
    const auto d=paired<double>([&]{return actual.value->roundedDouble();},[&]{return old.roundedDouble();});++terminalPairs;
    const auto n=paired<long double>([&]{return actual.value->rounded();},[&]{return old.rounded();});++terminalPairs;
    if(std::fegetround()!=FE_TONEAREST){check(!d.value&&d.error=="[FiniteVolume] binary64 round-to-nearest-even required"
        &&!n.value&&n.error=="[FiniteVolume] native round-to-nearest-even required","original directed terminal RNE gates were bypassed");}
    else if(!wideEligible){check(!d.value&&d.error=="[FiniteVolume] binary64 gradual underflow required"
        &&!n.value&&n.error=="[FiniteVolume] gradual underflow required","original actual FTZ/DAZ terminal gates were bypassed");}
    if(wideEligible&&wide.valid&&std::fegetround()==FE_TONEAREST){const auto wanted=ieeeRne(wide);
        if(wanted){check(d.value&&bits(*d.value)==*wanted,"terminal binary64 RNE differs from independent exact integer tie/quantum oracle");++rneChecks;
            if constexpr(sizeof(long double)==8&&std::numeric_limits<long double>::digits==53){check(n.value&&bits(static_cast<double>(*n.value))==*wanted,"terminal native binary64 differs from independent integer oracle");}}
        else{check(!d.value&&d.error=="[FiniteVolume] final binary64 dyadic overflow","original terminal binary64 overflow was changed");++terminalOverflows;}}
    if(expectPositiveZero&&d.value)check(bits(*d.value)==0,"first signed-zero product lost original canonical positive zero");
    ++cases;
}
std::uint64_t randomWord()
{static std::uint64_t state=UINT64_C(0xd1b54a32d192ed03);state^=state>>12;state^=state<<25;state^=state>>27;return state*UINT64_C(2685821657736338717);}
double randomFinite()
{return fromBits((randomWord()&UINT64_C(0x800fffffffffffff))|((randomWord()%2047)<<52));}
void normalAndExtreme()
{
    const double tiny=std::numeric_limits<double>::denorm_min(),maximum=std::numeric_limits<double>::max(),minimum=std::numeric_limits<double>::min(),full=fromBits(UINT64_C(0x3fffffffffffffff));
    tape({});tape({-1.});tape({1.,2.,3.,4.,5.});tape({-1.,2.,3.,4.,5.,6.});
    for(int count:{2,3,4})for(int signs=0;signs<(1<<count);++signs){const long double a=(signs&1)?-full:full,b=(signs&2)?-1.5L:1.5L,c=(signs&4)?-1.25L:1.25L,d=(signs&8)?-1.75L:1.75L;
        if(count==2)tape({a,b});else if(count==3)tape({a,b,c});else tape({a,b,c,d});}
    tape({full,full});tape({full,full,full});tape({full,full,full,full});
    tape({maximum,maximum});tape({maximum,maximum,maximum});tape({maximum,maximum,maximum,maximum});
    tape({-maximum,maximum,maximum,maximum});tape({tiny,tiny});tape({tiny,tiny,tiny});tape({-tiny,tiny,tiny,tiny});
    tape({maximum,tiny,maximum,tiny});tape({minimum,tiny,2.,.5L});tape({-maximum,tiny});
    tape({std::ldexp(1.,-537),std::ldexp(1.,-538)});tape({-std::ldexp(1.,-537),std::ldexp(1.,-538)});
    tape({3*std::ldexp(1.,-537),std::ldexp(1.,-538)});tape({1+std::ldexp(1.,-52),1+std::ldexp(1.,-52)});
    tape({0.,-maximum,tiny},0,true);tape({-0.,maximum,maximum,maximum},0,true);
    for(int j=0;j<18;++j){const long double a=randomFinite(),b=randomFinite(),c=randomFinite(),d=randomFinite();tape({a,b});tape({a,b,c});tape({a,b,c,d});}
}
void invalidAndShortCircuit()
{
    const long double q=std::numeric_limits<double>::quiet_NaN(),inf=std::numeric_limits<double>::infinity(),s=fromBits(UINT64_C(0x7ff0000000000001));
    tape({q,0.},1);tape({inf,-0.},1);tape({s,0.},1);tape({0.,q},0,true);tape({-0.,inf},0,true);tape({-0.,s},0,true);
    tape({1.,q,0.},1);tape({1.,0.,q},0,true);tape({1.,s,0.},1);tape({1.,-0.,s},0,true);
    tape({1.,2.,q,0.},1);tape({1.,2.,0.,q},0,true);tape({1.,2.,s,0.},1);tape({1.,2.,-0.,s},0,true);
    tape({1.,q,inf,0.},1);tape({1.,inf,q,0.},1);
}
void controls()
{
    const double tiny=std::numeric_limits<double>::denorm_min();const Floating nearest;
    for(int mode:{FE_UPWARD,FE_DOWNWARD,FE_TOWARDZERO}){check(std::fesetround(mode)==0&&std::fegetround()==mode,"small product cannot establish actual directed controls");
        tape({1.25L,-1.5L,1.75L});tape({tiny,1.5L});tape({-0.,std::numeric_limits<double>::quiet_NaN()},0,true);++controlCases;nearest.restore();}
#if defined(_MSC_VER) && defined(_M_X64)
    const unsigned base=_mm_getcsr();for(unsigned mode:{1u<<15,1u<<6,(1u<<15)|(1u<<6)}){_mm_setcsr(base|mode);
        tape({1.25L,-1.5L,1.75L});tape({tiny,1.5L});tape({-tiny,tiny,2.,.5L});tape({0.,std::numeric_limits<double>::infinity()},0,true);++controlCases;nearest.restore();}
#endif
}
void run()
{
    Restore restore;check(std::fesetround(FE_TONEAREST)==0,"small product requires actual nearest baseline");
#if defined(_MSC_VER) && defined(_M_X64)
    _mm_setcsr(_mm_getcsr()&~((1u<<15)|(1u<<6)));
#endif
    normalAndExtreme();invalidAndShortCircuit();controls();
    check(cases>=120&&exactChains>100&&wideChecks>90&&fpPairs>cases&&terminalPairs>150&&rneChecks>=50&&terminalOverflows>0
        &&originalGuards>=8&&zeroShortCircuits>=10&&controlCases>=3&&fullWidthProducts>0,"native small-product finite exact/FP/error coverage incomplete");
    std::cout<<"NATIVE_SMALL_PRODUCT_PASS,cases,"<<cases<<",constructor_multiplied_exact_chains,"<<exactChains<<",independent_wide_integer_products,"<<wideChecks
        <<",full_native_rawCSR_errno_pairs,"<<fpPairs<<",terminal_pairs,"<<terminalPairs<<",independent_RNE_bits,"<<rneChecks<<",original_terminal_overflows,"<<terminalOverflows
        <<",original_nonfinite_guards,"<<originalGuards<<",first_zero_short_circuits,"<<zeroShortCircuits<<",actual_control_modes,"<<controlCases
        <<",at_least_193bit_products,"<<fullWidthProducts<<",no_32a_sum_oracle,no_trap_or_resource_injection,no_speed_claim\n";
}
}

namespace DyadicTerminalRneContract
{
namespace FV=LAMP::Geometry::Detail::FiniteVolume;
using Integer=boost::multiprecision::cpp_int;
using UInt=FV::Internal::UInt;
using Wide=NativeSmallProductContract::Wide;
using Floating=NativeSmallProductContract::Floating;
std::size_t cases=0,literalPairs=0,independentBits=0,overflows=0,topOverflows=0,carryOverflows=0,
    negativeZeros=0,canonicalZeros=0,normalBoundary=0,largestBits=0,guardCases=0,
    inputCancellationPairs=0,exactCancellationZeros=0,tieEven=0,tieOdd=0,belowHalf=0,aboveHalf=0;
std::array<std::size_t,6> entryStates{};
std::vector<int> actualNormalizedShifts;
struct Fixture {Wide wide;FV::Dyadic value;UInt magnitude;};
// Only public exact operations construct the subject. No new terminal helper,
// private Dyadic representation, or production limb decoder is used by an oracle.
FV::Dyadic publicDyadic(const Integer& magnitude,int exponent,bool negative)
{
    check(magnitude>=0,"terminal oracle requires nonnegative independent integer");
    if(magnitude==0)return FV::Dyadic(negative?-0.L:0.L);
    FV::Dyadic value;const int width=static_cast<int>(boost::multiprecision::msb(magnitude))+1;
    for(int chunk=(width+31)/32;chunk-->0;){value=value.multiplied(4294967296.L);
        const Integer digit=(magnitude>>(32*chunk))&Integer(UINT64_C(0xffffffff));
        value.add(static_cast<long double>(digit.convert_to<std::uint32_t>()));}
    while(exponent){const int step=std::max(-512,std::min(512,exponent));
        value=value.multiplied(std::ldexp(1.L,step));exponent-=step;}
    return negative?value.multiplied(-1.L):value;
}
Fixture fixture(Integer magnitude,int exponent,bool negative=false,bool cancelNoise=false)
{
    Fixture result;result.wide.magnitude=magnitude;result.wide.exponent=exponent;result.wide.negative=negative;
    if(cancelNoise&&magnitude!=0){const Integer noise=Integer(1)<<(static_cast<int>(boost::multiprecision::msb(magnitude))+130);
        result.value=publicDyadic(magnitude+noise,exponent,false);result.value.add(publicDyadic(noise,exponent,false),-1);
        if(negative)result.value=result.value.multiplied(-1.L);++inputCancellationPairs;}
    else result.value=publicDyadic(magnitude,exponent,negative);
    // Public Dyadic normalization removes powers of two. A true half tape with
    // a large raw shift consequently has actual normalized shift one; do not
    // claim it tests an unreachable private wide-even representation.
    if(magnitude==0){result.wide.exponent=0;result.wide.negative=false;}
    else{const auto trailing=boost::multiprecision::lsb(magnitude);result.wide.magnitude>>=trailing;result.wide.exponent+=static_cast<int>(trailing);}
    result.magnitude=UInt(UInt::Integer(result.wide.magnitude));return result;
}
// Original terminal UInt algorithm, including its original environment prefix,
// right/remainder/compare/odd/convert/scalbn sequence and first error strings.
// Inputs come from the independent cpp_int tape, never the subject's internals.
template<class T>T originalUIntRne(const Fixture& source)
{
    if constexpr(std::is_same_v<T,double>)FV::Internal::binary64Environment();else FV::Internal::environment();
    if(source.magnitude.zero())return 0;
    constexpr int precision=std::numeric_limits<T>::digits;
    constexpr int least=std::numeric_limits<T>::min_exponent-precision;
    const char* overflow=std::is_same_v<T,double>?"[FiniteVolume] final binary64 dyadic overflow":"[FiniteVolume] final native dyadic overflow";
    const char* lost=std::is_same_v<T,double>?"[FiniteVolume] binary64 dyadic conversion lost positive quantum":"[FiniteVolume] native dyadic conversion lost positive quantum";
    const int top=source.magnitude.bits()-1+source.wide.exponent;
    FV::Internal::require(top<std::numeric_limits<T>::max_exponent,overflow);
    const int quantum=std::max(top-(precision-1),least),shift=quantum-source.wide.exponent;UInt integer;
    if(shift>0){integer=source.magnitude.right(shift);auto remainder=source.magnitude;remainder.subtract(integer.shifted(shift));
        const int half=remainder.compare(UInt(1).shifted(shift-1));
        if(half>0||(half==0&&integer.odd()))integer.add(UInt(1));}
    else integer=source.magnitude.shifted(-shift);
    if(integer.zero())return std::copysign(T(0),source.wide.negative?T(-1):T(1));
    const T significand=integer.n.template convert_to<T>();const T value=std::scalbn(significand,quantum);
    FV::Internal::require(std::isfinite(value),overflow);FV::Internal::require(value>0,lost);
    return source.wide.negative?-value:value;
}
void prime(unsigned variant)
{
    const std::array<int,6> flags{{0,FE_INEXACT,FE_UNDERFLOW|FE_INEXACT,FE_INVALID,FE_DIVBYZERO,FE_OVERFLOW|FE_INEXACT}};
    variant%=static_cast<unsigned>(flags.size());
    check(std::feclearexcept(FE_ALL_EXCEPT)==0,"terminal oracle cannot clear actual entry status");
#if defined(_MSC_VER) && defined(_M_X64)
    _mm_setcsr(_mm_getcsr()&~0x3fu);
#endif
    if(flags[variant])check(std::feraiseexcept(flags[variant])==0,"terminal oracle cannot raise actual entry status");
#if defined(_MSC_VER) && defined(_M_X64)
    if(variant&1u)_mm_setcsr(_mm_getcsr()|0x2u); // Real denormal sticky, not synthetic fenv bytes.
#endif
    errno=(variant&1u)?EDOM:E2BIG;++entryStates[variant];
}
void observe(const Fixture& source)
{
    if(source.magnitude.zero())return;const int width=source.magnitude.bits();largestBits=std::max(largestBits,static_cast<std::size_t>(width));
    const int top=width-1+source.wide.exponent,shift=std::max(top-52,-1074)-source.wide.exponent;
    actualNormalizedShifts.push_back(shift);
    if(shift>0){const Integer q=source.wide.magnitude>>shift,rest=source.wide.magnitude-(q<<shift),half=Integer(1)<<(shift-1);
        if(rest<half)++belowHalf;else if(rest>half)++aboveHalf;else if(boost::multiprecision::bit_test(q,0))++tieOdd;else ++tieEven;}
}
void tape(const Fixture& source,unsigned variant=0)
{
    const bool nearest=std::fegetround()==FE_TONEAREST,gradual=NativeSmallProductContract::gradualControls();
    const auto wanted=NativeSmallProductContract::ieeeRne(source.wide);
    if(nearest&&gradual)observe(source);
    prime(variant);const auto d=NativeSmallProductContract::paired<double>([&]{return source.value.roundedDouble();},[&]{return originalUIntRne<double>(source);});++literalPairs;
    prime(variant);const auto n=NativeSmallProductContract::paired<long double>([&]{return source.value.rounded();},[&]{return originalUIntRne<long double>(source);});++literalPairs;
    if(!nearest){check(!d.value&&d.error=="[FiniteVolume] binary64 round-to-nearest-even required"
        &&!n.value&&n.error=="[FiniteVolume] native round-to-nearest-even required","terminal original directed-mode first gate changed");++guardCases;}
    else if(!gradual){check(!d.value&&d.error=="[FiniteVolume] binary64 gradual underflow required"
        &&!n.value&&n.error=="[FiniteVolume] gradual underflow required","terminal original FTZ/DAZ first gate changed");++guardCases;}
    else if(wanted){check(d.value&&NativeSmallProductContract::bits(*d.value)==*wanted,"terminal result differs from independent integer IEEE encoding");++independentBits;
        if constexpr(sizeof(long double)==8&&std::numeric_limits<long double>::digits==53)
            check(n.value&&NativeSmallProductContract::bits(static_cast<double>(*n.value))==*wanted,"native terminal differs from independent binary64 integer encoding");
        if(*wanted==UINT64_C(0x8000000000000000))++negativeZeros;
        if(source.magnitude.zero()){check(*wanted==0,"exact zero must retain original canonical positive zero");++canonicalZeros;}
        if((*wanted&UINT64_C(0x7fffffffffffffff))==UINT64_C(0x0010000000000000))++normalBoundary;
    }else{check(!d.value&&d.error=="[FiniteVolume] final binary64 dyadic overflow","terminal original overflow changed");
        if constexpr(sizeof(long double)==8&&std::numeric_limits<long double>::digits==53)
            check(!n.value&&n.error=="[FiniteVolume] final native dyadic overflow","native original overflow changed");
        const int top=source.magnitude.bits()-1+source.wide.exponent;if(top>=1024)++topOverflows;else ++carryOverflows;++overflows;}
    ++cases;
}
void finiteAndExtreme()
{
    unsigned entry=0;auto emit=[&](const Integer& m,int exponent,bool both=false,bool noise=false){tape(fixture(m,exponent,false,noise),entry++);
        if(both)tape(fixture(m,exponent,true,noise),entry++);};
    emit(0,0,true);emit(1,0);emit(3,-1);emit((Integer(1)<<52)+1,-52);
    emit(1,-1074,true);emit(1,-1075,true);emit(1,-1076,true);
    for(int m:{3,5,7})emit(m,-1075,true);
    emit((Integer(1)<<52)-1,-1074,true);emit((Integer(1)<<53)-1,-1075,true);
    emit((Integer(1)<<54)-3,-1076);emit((Integer(1)<<54)-1,-1076);
    emit((Integer(1)<<53)-1,971,true);emit((Integer(1)<<54)-1,970,true);
    emit((Integer(1)<<55)-3,969);emit((Integer(1)<<55)-1,969);emit(1,1024,true);
    emit((Integer(1)<<12000)-1,-13074,true);emit((Integer(1)<<11999)+1,-13074);
    emit((Integer(1)<<12000)-1,-25000,true);emit((Integer(1)<<12000)-1,-12052);
    // Odd near-half remainders retain each actual limb-boundary shift. True
    // ties also exercise both even/odd q, but normalize to shift one.
    const std::array<int,13> shifts{{63,64,65,127,128,129,4095,4096,4097,8191,8192,8193,11947}};
    for(const int shift:shifts)for(int delta:{-1,0,1}){const Integer q=(Integer(1)<<52)+((shift+delta)&1);
        const Integer m=(q<<shift)+(Integer(1)<<(shift-1))+delta;
        const bool negative=(entry&1u)!=0;tape(fixture(m,-52-shift,negative),entry++);}
    std::uint64_t state=UINT64_C(0xd1b54a32d192ed03);
    for(int width:{64,65,4096,4097,12000}){Integer m=0;
        for(int bit=0;bit<width;bit+=64){state^=state>>12;state^=state<<25;state^=state>>27;m|=Integer(state*UINT64_C(2685821657736338717))<<bit;}
        m&=(Integer(1)<<width)-1;m|=(Integer(1)<<(width-1))|1;emit(m,-width,true);}
    emit((Integer(1)<<11999)+12345,-11999,true,true);emit((Integer(1)<<127)+31,-1127,true,true);
    const auto repeated=fixture((Integer(1)<<11999)+(Integer(1)<<6401)+1,-11999,true);
    for(unsigned flags=0;flags<6;++flags)tape(repeated,flags);
    for(bool sign:{false,true}){auto zero=fixture(0,0,sign);auto large=publicDyadic((Integer(1)<<11999)+7,-11999,sign);const auto copy=large;large.add(copy,-1);
        check(large.compare(FV::Dyadic{})==0,"public wide exact cancellation failed");zero.value=std::move(large);tape(zero,entry++);++exactCancellationZeros;}
}
void controls()
{
    const std::array<Fixture,3> probes{{fixture(0,0,true),fixture((Integer(1)<<11999)+1,-11999),fixture(1,1024)}};
    const Floating baseline;
    for(int mode:{FE_UPWARD,FE_DOWNWARD,FE_TOWARDZERO}){check(std::fesetround(mode)==0&&std::fegetround()==mode,"terminal oracle cannot establish actual directed controls");
        for(const auto& source:probes)tape(source,static_cast<unsigned>(guardCases));baseline.restore();}
#if defined(_MSC_VER) && defined(_M_X64)
    for(unsigned mode:{1u<<15,1u<<6,(1u<<15)|(1u<<6)}){_mm_setcsr(_mm_getcsr()|mode);
        for(const auto& source:probes)tape(source,static_cast<unsigned>(guardCases));baseline.restore();}
    check(std::fesetround(FE_UPWARD)==0,"terminal oracle cannot establish combined unsupported controls");_mm_setcsr(_mm_getcsr()|(1u<<15));
    for(const auto& source:probes)tape(source,static_cast<unsigned>(guardCases));baseline.restore();
#endif
}
void run()
{
    NativeSmallProductContract::Restore restore;check(std::fesetround(FE_TONEAREST)==0,"terminal oracle requires actual nearest baseline");
#if defined(_MSC_VER) && defined(_M_X64)
    _mm_setcsr(_mm_getcsr()&~((1u<<15)|(1u<<6)));unsigned native=0;
    check(_controlfp_s(&native,0,0)==0&&(native&_MCW_EM)==_MCW_EM&&(_mm_getcsr()&0x1f80u)==0x1f80u,
        "terminal oracle requires masked original exception controls; traps are not injected");
#endif
    const auto fpBefore=NativeSmallProductContract::fpPairs;finiteAndExtreme();controls();
    for(int shift:{-52,0,1,63,64,65,127,128,129,4095,4096,4097,8191,8192,8193,11947})
        check(std::find(actualNormalizedShifts.begin(),actualNormalizedShifts.end(),shift)!=actualNormalizedShifts.end(),"terminal normalized limb/quantum boundary tape missing");
    check(cases>80&&literalPairs==2*cases&&NativeSmallProductContract::fpPairs-fpBefore==literalPairs&&independentBits>50
        &&overflows>0&&topOverflows>0&&carryOverflows>0&&negativeZeros>0&&canonicalZeros>0&&normalBoundary>0&&largestBits>=12000
        &&tieEven>0&&tieOdd>0&&belowHalf>0&&aboveHalf>0&&guardCases>=9&&inputCancellationPairs==4&&exactCancellationZeros==2,
        "terminal independent literal/IEEE/FP/edge coverage incomplete");
    for(const auto count:entryStates)check(count>0,"terminal real entry-sticky variant missing");
    std::cout<<"DYADIC_TERMINAL_RNE_PASS,tapes,"<<cases<<",original_UInt_terminal_pairs,"<<literalPairs
        <<",full_native_rawCSR_errno_pairs,"<<(NativeSmallProductContract::fpPairs-fpBefore)<<",independent_IEEE_bits,"<<independentBits
        <<",largest_normalized_bits,"<<largestBits<<",top_overflows,"<<topOverflows<<",rounding_carry_overflows,"<<carryOverflows
        <<",negative_rounded_zeros,"<<negativeZeros<<",canonical_exact_zeros,"<<canonicalZeros<<",normal_boundary_results,"<<normalBoundary
        <<",tie_even,"<<tieEven<<",tie_odd,"<<tieOdd<<",below_half,"<<belowHalf<<",above_half,"<<aboveHalf
        <<",original_control_guards,"<<guardCases<<",public_noisy_cancellations,"<<inputCancellationPairs<<",public_exact_zero_cancellations,"<<exactCancellationZeros
        <<",no_private_Dyadic_or_new_helper_probe,no_trap_or_allocation_fault_injection,no_speed_claim\n";
}
}

namespace DyadicIntegerOperationsContract
{
namespace FV=LAMP::Geometry::Detail::FiniteVolume;
using Integer=boost::multiprecision::cpp_int;
using Rational=boost::multiprecision::cpp_rational;
using Wide=NativeSmallProductContract::Wide;
using Floating=NativeSmallProductContract::Floating;
using Subject=std::pair<FV::Dyadic,Rational>;
std::size_t additions=0,aliases=0,comparisons=0,observations=0,independentHeads=0,
    zeroResults=0,guards=0,scalarPrefixes=0,controlModes=0,entryVariant=0;
std::array<std::size_t,6> entries{};
Rational exact(const Wide& w)
{
    check(w.valid,"integer-flow oracle requires finite original IEEE input");
    Rational value(w.magnitude);
    if(w.exponent>=0)value*=Integer(1)<<w.exponent;else value/=Integer(1)<<(-w.exponent);
    if(w.negative)value=-value;return value;
}
Rational ieee(double value)
{return exact(NativeSmallProductContract::independent({static_cast<long double>(value)}));}
Rational power(int exponent)
{
    Rational value(1);
    if(exponent>=0)value*=Integer(1)<<exponent;else value/=Integer(1)<<(-exponent);
    return value;
}
Wide independent(const Rational& value)
{
    Wide result;Integer numerator=boost::multiprecision::numerator(value);
    const Integer denominator=boost::multiprecision::denominator(value);
    check(denominator>0&&(denominator&(denominator-1))==0,"integer-flow rational denominator is not a binary power");
    result.negative=numerator<0;if(result.negative)numerator=-numerator;
    result.magnitude=std::move(numerator);result.exponent=-static_cast<int>(boost::multiprecision::msb(denominator));
    if(result.magnitude==0){result.exponent=0;result.negative=false;}
    return result;
}
Subject scalar(double value)
{return {FV::Dyadic(static_cast<long double>(value)),ieee(value)};}
Subject wide(const Integer& magnitude,int exponent,bool negative=false)
{
    // Public construction is separately checked against the rational tape;
    // neither private Dyadic storage nor its compare() supplies an oracle.
    Wide reference;reference.magnitude=magnitude;reference.exponent=exponent;reference.negative=negative;
    return {DyadicTerminalRneContract::publicDyadic(magnitude,exponent,negative),exact(reference)};
}
void prime()
{
    const auto variant=static_cast<unsigned>(entryVariant++%entries.size());++entries[variant];
    DyadicTerminalRneContract::prime(variant);
}
void observe(const Subject& source)
{
    auto actual=source.first;auto remainder=source.second;++observations;
    // Each fixture has at most a few nonzero binary64 heads. Normalize its
    // independent rational scale, check IEEE RNE, then remove that known head.
    // This exposes tails 12000 bits below the original top without rounding
    // them away or constructing an expected value with the subject's add().
    for(unsigned head=0;head<8;++head){
        if(remainder==0){check(NativeSmallProductContract::bits(actual.roundedDouble())==0,
                "integer-flow exact zero lost its canonical positive sign");++zeroResults;return;}
        const auto reference=independent(remainder);
        const int top=static_cast<int>(boost::multiprecision::msb(reference.magnitude))+reference.exponent;
        int shift=-top;remainder*=power(shift);
        while(shift){const int part=std::max(-512,std::min(512,shift));
            actual=actual.multiplied(std::ldexp(1.L,part));shift-=part;}
        const auto wanted=NativeSmallProductContract::ieeeRne(independent(remainder));
        check(wanted.has_value(),"integer-flow normalized rational head unexpectedly overflowed");
        const auto got=actual.roundedDouble();
        check(NativeSmallProductContract::bits(got)==*wanted,"integer-flow sum differs from independent rational IEEE head");
        if constexpr(sizeof(long double)==8&&std::numeric_limits<long double>::digits==53)
            check(NativeSmallProductContract::bits(static_cast<double>(actual.rounded()))==*wanted,
                "integer-flow native head differs from independent rational IEEE bits");
        const auto known=NativeSmallProductContract::fromBits(*wanted);
        actual.add(static_cast<long double>(known),-1);remainder-=ieee(known);++independentHeads;
    }
    check(false,"integer-flow fixture exceeded its bounded independent head tape");
}
void compare(const Subject& a,const Subject& b)
{
    const int expected=a.second<b.second?-1:a.second>b.second?1:0;
    prime();const Floating before;const int forward=a.first.compare(b.first),reverse=b.first.compare(a.first);
    const Floating after;NativeSmallProductContract::sameFp(before,after);before.restore();
    check(forward==expected&&reverse==-expected,"integer-flow compare differs from independent signed rational order");
    ++comparisons;observe(a);observe(b); // compare must preserve both operands.
}
Subject add(const Subject& a,const Subject& b,int sign,bool alias=false)
{
    Subject result=a;
    Rational expected=a.second;if(alias)expected*=1+sign;else expected+=sign*b.second;
    prime();const Floating before;
    if(alias)result.first.add(result.first,sign);else result.first.add(b.first,sign);
    const Floating after;NativeSmallProductContract::sameFp(before,after);before.restore();
    result.second=expected;++additions;if(alias)++aliases;
    observe(result);observe(a);observe(b); // result copy and actual RHS have distinct ownership.
    compare(result,scalar(0));return result;
}
void scalarAdd(const Subject& a,double value,int sign)
{
    Subject result=a;const Rational expected=a.second+sign*ieee(value);
    prime();const Floating entry;result.first.add(static_cast<long double>(value),sign);const Floating actualAfter;
    entry.restore();const FV::Dyadic unchangedPrefix(static_cast<long double>(value));(void)unchangedPrefix;
    const Floating prefixAfter;NativeSmallProductContract::sameFp(actualAfter,prefixAfter);entry.restore();
    result.second=expected;observe(result);observe(a);++scalarPrefixes;
}
void guard(const Subject& a,long double input,int sign,bool typed)
{
    auto result=a.first;const FV::Dyadic zero;
    const char* expected=typed||std::isfinite(input)?"[FiniteVolume] invalid dyadic addition sign":"[FiniteVolume] nonfinite dyadic input";
    prime();const Floating entry;
    const auto got=NativeSmallProductContract::invoke<int>([&]{if(typed)result.add(zero,sign);else result.add(input,sign);return 1;});
    const Floating actualAfter;entry.restore();
    // Only unchanged public constructor/gate prefixes form the FP/first-error
    // reference. No old add/compare/normalization algorithm is copied here.
    const auto prefix=NativeSmallProductContract::invoke<int>([&]{
        if(!typed){const FV::Dyadic originalPrefix(input);(void)originalPrefix;}
        FV::Internal::require(false,"[FiniteVolume] invalid dyadic addition sign");return 1;});
    const Floating prefixAfter;NativeSmallProductContract::sameFp(actualAfter,prefixAfter);entry.restore();
    check(!got.value&&!prefix.value&&got.error==expected&&got.error==prefix.error,
        "integer-flow invalid sign/nonfinite original first error changed");
    observe({std::move(result),a.second});observe(a);++guards;
}
void finiteTapes()
{
    const auto one=scalar(1),zero=scalar(-0.);
    for(int boundary:{0,1,31,32,63,64,127,128,255,256,257,12000}){
        const Integer high=Integer(1)<<boundary;
        const auto a=wide(high+1,-boundary),tail=wide(Integer(1),-boundary);
        observe(a);observe(tail);compare(a,one);compare(tail,zero);
        const auto negative=wide(high+1,-boundary,true);compare(negative,a);
        for(int sign:{1,-1}){
            (void)add(a,tail,sign);(void)add(tail,a,sign);(void)add(a,one,sign);
            (void)add(negative,tail,sign);(void)add(a,a,sign,true);
        }
        if(boundary){
            const auto borrow=wide(high-1,-boundary);(void)add(borrow,tail,1);(void)add(tail,borrow,-1);
            compare(wide(Integer(3),-boundary),wide(Integer(5),-boundary)); // same exponent, differing magnitudes.
            compare(a,wide((Integer(1)<<(boundary+1))+1,-boundary)); // different top bit.
            if(boundary>1)compare(a,wide((Integer(1)<<(boundary-1))+1,1-boundary)); // same top, differing exponents/tails.
        }
    }
    const auto huge=wide(Integer(1),1000),micro=wide(Integer(1),-12000);
    const auto combined=add(huge,micro,1);const auto exposed=add(combined,huge,-1);compare(exposed,micro);
    for(int sign:{1,-1}){(void)add(zero,zero,sign);(void)add(zero,one,sign);(void)add(one,zero,sign);}
    const auto tiny=scalar(std::numeric_limits<double>::denorm_min());
    for(double value:{-0.,1.5,-1.5,std::numeric_limits<double>::denorm_min(),std::numeric_limits<double>::max()})
        for(int sign:{1,-1})scalarAdd(tiny,value,sign);
    for(std::uint64_t word:{UINT64_C(0x0010000000000000),UINT64_C(0x7fefffffffffffff),
            UINT64_C(0x8000000000000001),UINT64_C(0xbff0000000000001),UINT64_C(0x3fffffffffffffff)}){
        const auto a=scalar(NativeSmallProductContract::fromBits(word));
        for(int sign:{1,-1}){(void)add(a,tiny,sign);(void)add(a,a,sign,true);}compare(a,one);
    }
}
void invalidTapes()
{
    const auto source=wide((Integer(1)<<257)+1,-257);
    for(int sign:{0,2,-2}){guard(source,0.L,sign,true);guard(source,0.L,sign,false);guard(source,1.25L,sign,false);}
    for(const auto word:{UINT64_C(0x7ff0000000000000),UINT64_C(0x7ff8000000000001),UINT64_C(0x7ff0000000000001)})
        for(int sign:{0,1})guard(source,static_cast<long double>(NativeSmallProductContract::fromBits(word)),sign,false);
}
void controls()
{
    const auto left=wide((Integer(1)<<257)+1,-257,true),right=wide(Integer(1),-12000);
    const Floating baseline;
    const auto one=[&]{
        auto result=left;prime();const Floating entry;result.first.add(right.first);
        const int ordering=left.first.compare(right.first);const Floating after;
        NativeSmallProductContract::sameFp(entry,after);entry.restore();
        check(ordering==-1,"integer-flow changed exact order under nondefault controls");
        baseline.restore();result.second+=right.second;observe(result);observe(left);observe(right);++controlModes;
    };
    for(int mode:{FE_UPWARD,FE_DOWNWARD,FE_TOWARDZERO}){
        check(std::fesetround(mode)==0,"integer-flow cannot establish directed FP controls");one();}
#if defined(_MSC_VER) && defined(_M_X64)
    for(unsigned mode:{1u<<15,1u<<6,(1u<<15)|(1u<<6)}){_mm_setcsr(_mm_getcsr()|mode);one();}
#endif
    baseline.restore();
}
void run()
{
    NativeSmallProductContract::Restore restore;
    check(sizeof(long double)==8&&std::numeric_limits<long double>::digits==53,
        "integer-flow IEEE fixture requires actual binary64 native operands");
    check(std::fesetround(FE_TONEAREST)==0,"integer-flow requires nearest baseline");
#if defined(_MSC_VER) && defined(_M_X64)
    _mm_setcsr(_mm_getcsr()&~((1u<<15)|(1u<<6)));unsigned control=0;
    check(_controlfp_s(&control,0,0)==0&&(control&_MCW_EM)==_MCW_EM&&(_mm_getcsr()&0x1f80u)==0x1f80u,
        "integer-flow requires masked exception controls; traps are not injected");
#endif
    const auto fpBefore=NativeSmallProductContract::fpPairs;finiteTapes();invalidTapes();controls();
    for(auto count:entries)check(count>0,"integer-flow real sticky/errno entry variant missing");
    check(additions>130&&aliases>=30&&comparisons>170&&independentHeads>200&&guards==15&&scalarPrefixes==10&&controlModes>=3,
        "integer-flow independent rational/alias/boundary/FP coverage incomplete");
    std::cout<<"DYADIC_INTEGER_OPERATIONS_PASS,additions,"<<additions<<",self_alias_additions,"<<aliases
        <<",independent_rational_order_pairs,"<<comparisons<<",independent_IEEE_heads,"<<independentHeads
        <<",operand_preservation_observations,"<<observations<<",exact_zero_results,"<<zeroResults
        <<",original_first_error_guards,"<<guards<<",original_scalar_FP_prefix_pairs,"<<scalarPrefixes
        <<",nondefault_integer_control_modes,"<<controlModes<<",full_native_rawCSR_errno_pairs,"<<(NativeSmallProductContract::fpPairs-fpBefore)
        <<",maximum_alignment_boundary,12000,no_private_Dyadic_or_optimized_compare_oracle,no_old_integer_algorithm_copy,no_fault_injection,no_speed_claim\n";
}
}
namespace NativeSmallExactContract
{
namespace FV=LAMP::Geometry::Detail::FiniteVolume;
#if defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE) && !defined(__FAST_MATH__) && (!defined(_M_FP_FAST) || !_M_FP_FAST)
using Small=FV::NativeSmallExact;
using Integer=boost::multiprecision::cpp_int;
using Wide=NativeSmallProductContract::Wide;
using Floating=NativeSmallProductContract::Floating;
std::size_t tapes=0,productPairs=0,combinePairs=0,terminalPairs=0,literalTerminalPairs=0,independentBits=0,
    exactReferences=0,firstErrors=0,firstZeros=0,additionGuards=0,countGuards=0,controlGuards=0,
    negativeZeros=0,canonicalZeros=0,terminalOverflows=0,hugeCancellations=0,copyMovePairs=0,
    maximumProductBits=0,maximumAlignedSumBits=0;
unsigned entryVariant=0;
struct Built
{
    NativeSmallProductContract::Outcome<Small::Product> actual;
    NativeSmallProductContract::Outcome<FV::Dyadic> literal;
    Wide wide;
};
bool mathQualified(){return std::fegetround()==FE_TONEAREST&&NativeSmallProductContract::gradualControls();}
Built build(std::initializer_list<long double> factors)
{
    Built r;r.wide.valid=false;
    if(mathQualified())r.wide=NativeSmallProductContract::independent(factors);
    DyadicTerminalRneContract::prime(entryVariant++);const Floating entry;
    r.actual=NativeSmallProductContract::invoke<Small::Product>([&]{return Small::product(factors);});const Floating actualAfter;
    entry.restore();r.literal=NativeSmallProductContract::invoke<FV::Dyadic>([&]{return FV::Dyadic::product(factors);});const Floating literalAfter;
    NativeSmallProductContract::sameFp(actualAfter,literalAfter);++productPairs;
    check(r.actual.value.has_value()==r.literal.value.has_value()&&r.actual.error==r.literal.error,
        "small exact product changed original scalar prefix or first error");
    if(!r.actual.value){check(r.actual.error=="[FiniteVolume] nonfinite dyadic input","small exact product changed original nonfinite message");++firstErrors;}
    else if(r.wide.valid){const auto expected=NativeSmallProductContract::expanded(r.wide);
        check(r.literal.value->compare(expected)==0,"literal product differs from independent IEEE-decoded integer product");++exactReferences;
        if(r.wide.magnitude==0)++firstZeros;
        else{const auto width=boost::multiprecision::msb(r.wide.magnitude)-boost::multiprecision::lsb(r.wide.magnitude)+1;
            maximumProductBits=std::max(maximumProductBits,static_cast<std::size_t>(width));}}
    return r;
}
Wide sumWide(const Wide& a,int signA,const Wide& b,int signB)
{
    Wide r;r.valid=a.valid&&b.valid;
    if(!r.valid)return r;
    const int common=std::min(a.exponent,b.exponent);
    Integer x=a.magnitude<<(a.exponent-common),y=b.magnitude<<(b.exponent-common);
    if(a.negative!=(signA<0))x=-x;if(b.negative!=(signB<0))y=-y;
    const Integer total=x+y;r.negative=total<0;r.magnitude=total;if(r.negative)r.magnitude=-r.magnitude;r.exponent=common;
    if(r.magnitude==0){r.negative=false;r.exponent=0;}
    else maximumAlignedSumBits=std::max(maximumAlignedSumBits,static_cast<std::size_t>(boost::multiprecision::msb(r.magnitude)+1));
    return r;
}
DyadicTerminalRneContract::Fixture originalIntegerInput(const Wide& w)
{
    DyadicTerminalRneContract::Fixture input;input.wide=w;
    // This reference needs no subject value: reconstruct its original UInt from
    // independent exact inputs, never Product/Sum's private word representation.
    input.magnitude=FV::Internal::UInt(FV::Internal::UInt::Integer(w.magnitude));return input;
}
template<class F>void terminal(F actual,const FV::Dyadic& old,const Wide& wide)
{
    const bool nearest=std::fegetround()==FE_TONEAREST,gradual=NativeSmallProductContract::gradualControls();
    DyadicTerminalRneContract::prime(entryVariant++);
    const auto observed=NativeSmallProductContract::paired<double>(actual,[&]{return old.roundedDouble();});++terminalPairs;
    if(!nearest){check(!observed.value&&observed.error=="[FiniteVolume] binary64 round-to-nearest-even required",
        "small exact terminal changed original directed-control first gate");++controlGuards;return;}
    if(!gradual){check(!observed.value&&observed.error=="[FiniteVolume] binary64 gradual underflow required",
        "small exact terminal changed original FTZ/DAZ first gate");++controlGuards;return;}
    check(wide.valid,"small exact nearest/gradual tape lacks independent integer input");
    const auto input=originalIntegerInput(wide);const auto wanted=NativeSmallProductContract::ieeeRne(wide);
    DyadicTerminalRneContract::prime(entryVariant++);
    const auto literal=NativeSmallProductContract::paired<double>(actual,[&]{return DyadicTerminalRneContract::originalUIntRne<double>(input);});++literalTerminalPairs;
    if(wanted){check(observed.value&&literal.value&&NativeSmallProductContract::bits(*observed.value)==*wanted,
        "small exact terminal differs from independent unbounded IEEE/RNE result");++independentBits;
        if(*wanted==UINT64_C(0x8000000000000000))++negativeZeros;
        if(wide.magnitude==0){check(*wanted==0,"small exact cancellation/first zero lost canonical positive zero");++canonicalZeros;}}
    else{check(!observed.value&&observed.error=="[FiniteVolume] final binary64 dyadic overflow",
        "small exact terminal changed original overflow result");++terminalOverflows;}
}
void productTape(std::initializer_list<long double> factors)
{
    const auto r=build(factors);if(r.actual.value)terminal([&]{return r.actual.value->roundedDouble();},*r.literal.value,r.wide);++tapes;
}
void sumTape(std::initializer_list<long double> factorsA,int signA,std::initializer_list<long double> factorsB,int signB,bool exerciseCopies=false)
{
    const auto a=build(factorsA),b=build(factorsB);check(a.actual.value&&b.actual.value,"sum fixture requires two successfully constructed products");
    DyadicTerminalRneContract::prime(entryVariant++);const Floating entry;
    auto actual=NativeSmallProductContract::invoke<Small::Sum>([&]{return Small::combine(*a.actual.value,signA,*b.actual.value,signB);});const Floating actualAfter;
    entry.restore();auto literal=NativeSmallProductContract::invoke<FV::Dyadic>([&]{FV::Dyadic sum;sum.add(*a.literal.value,signA);sum.add(*b.literal.value,signB);return sum;});const Floating literalAfter;
    NativeSmallProductContract::sameFp(actualAfter,literalAfter);++combinePairs;
    check(actual.value.has_value()==literal.value.has_value()&&actual.error==literal.error,"small exact combine changed original two-add availability/first error");
    if(!actual.value){check(actual.error=="[FiniteVolume] invalid dyadic addition sign","small exact combine changed original sign guard");++additionGuards;++tapes;return;}
    const auto wide=sumWide(a.wide,signA,b.wide,signB);
    if(wide.valid){const auto expected=NativeSmallProductContract::expanded(wide);
        check(literal.value->compare(expected)==0,"literal two-add differs from independent signed integer aligned sum");++exactReferences;
        const auto over=[](const Wide& w){return w.magnitude!=0&&static_cast<int>(boost::multiprecision::msb(w.magnitude))+w.exponent>=1024;};
        if(wide.magnitude==0&&over(a.wide)&&over(b.wide))++hugeCancellations;}
    // Do not round either product before combining: huge exact cancellation
    // must reach its original final terminal without an artificial early overflow.
    terminal([&]{return actual.value->roundedDouble();},*literal.value,wide);
    if(exerciseCopies){const auto copy=*actual.value;auto moved=Small::Sum(copy);Small::Sum assigned;assigned=copy;Small::Sum movedAssigned;movedAssigned=std::move(moved);
        terminal([&]{return copy.roundedDouble();},*literal.value,wide);terminal([&]{return assigned.roundedDouble();},*literal.value,wide);
        terminal([&]{return movedAssigned.roundedDouble();},*literal.value,wide);copyMovePairs+=3;}
    ++tapes;
}
void productsAndFirstErrors()
{
    const long double tiny=std::numeric_limits<double>::denorm_min(),minimum=std::numeric_limits<double>::min(),maximum=std::numeric_limits<double>::max();
    const long double full=NativeSmallProductContract::fromBits(UINT64_C(0x3fffffffffffffff));
    for(int n:{2,3,4})for(int signs=0;signs<4;++signs){const long double a=(signs&1)?-full:full,b=(signs&2)?-1.5L:1.5L;
        if(n==2)productTape({a,b});else if(n==3)productTape({a,b,full});else productTape({a,b,full,full});}
    productTape({full,full,full,full});productTape({maximum,maximum,maximum,maximum});productTape({-maximum,maximum,maximum,maximum});
    productTape({tiny,tiny,tiny,tiny});productTape({-tiny,tiny,tiny,tiny});productTape({maximum,tiny,maximum,tiny});
    productTape({minimum,tiny,2.,.5L});productTape({tiny,.5L});productTape({-tiny,.5L});productTape({tiny,1.5L});
    productTape({minimum,1.});productTape({maximum,1.});productTape({maximum,2.});
    const long double q=std::numeric_limits<double>::quiet_NaN(),inf=std::numeric_limits<double>::infinity();
    const long double s=NativeSmallProductContract::fromBits(UINT64_C(0x7ff0000000000001));
    productTape({q,0.});productTape({inf,-0.});productTape({s,0.});productTape({0.,q});productTape({-0.,inf});productTape({-0.,s});
    productTape({1.,q,0.});productTape({1.,0.,q});productTape({1.,s,0.});productTape({1.,-0.,s});
    productTape({1.,2.,q,0.});productTape({1.,2.,0.,q});productTape({1.,2.,s,0.});productTape({1.,2.,-0.,s});
}
void sumsAndBoundaries()
{
    const long double tiny=std::numeric_limits<double>::denorm_min(),minimum=std::numeric_limits<double>::min(),maximum=std::numeric_limits<double>::max();
    const long double full=NativeSmallProductContract::fromBits(UINT64_C(0x3fffffffffffffff)),next=NativeSmallProductContract::fromBits(UINT64_C(0x3ff0000000000001));
    for(int signA:{-1,1})for(int signB:{-1,1}){sumTape({1.25L,1.5L},signA,{1.75L,.5L},signB);
        sumTape({maximum,maximum,maximum,maximum},signA,{tiny,tiny,tiny,tiny},signB,true);}
    sumTape({maximum,maximum,maximum,maximum},1,{maximum,maximum,maximum,maximum},-1,true);
    sumTape({maximum,maximum,maximum,maximum},-1,{maximum,maximum,maximum,maximum},1);
    sumTape({maximum,maximum,maximum,maximum},1,{maximum,maximum,maximum,maximum},1);
    sumTape({std::ldexp(1.L,1023),2.L},1,{std::ldexp(1.L,1023),2.L},-1);
    for(int sign:{-1,1}){sumTape({1.L,1.L},sign,{std::ldexp(1.L,-53),1.L},sign);
        sumTape({next,1.L},sign,{std::ldexp(1.L,-53),1.L},sign);
        sumTape({minimum-tiny,1.L},sign,{tiny,.5L},sign);
        sumTape({minimum-tiny,1.L},sign,{tiny,.25L},sign);sumTape({minimum-tiny,1.L},sign,{tiny,.75L},sign);
        sumTape({maximum,1.L},sign,{std::ldexp(1.L,970),1.L},sign);
        sumTape({maximum,1.L},sign,{std::ldexp(1.L,969),1.L},sign);
        sumTape({maximum,1.L},sign,{std::ldexp(1.L,970),1.5L},sign);}
    sumTape({tiny,.5L},-1,{0.,1.},1);sumTape({tiny,1.},1,{tiny,.5L},1);sumTape({tiny,1.},1,{tiny,1.},-1,true);
    for(int shift:{63,64,65,127,128,129})sumTape({1.L,1.L},1,{std::ldexp(1.L,-shift),1.L},(shift&1)?1:-1);
    sumTape({full,full,full,full},1,{full,full,full,full},-1);sumTape({0.,1.},-1,{-0.,1.},-1,true);
    sumTape({1.L,2.L},0,{3.L,4.L},1);sumTape({1.L,2.L},1,{3.L,4.L},0);
    sumTape({0.,1.},0,{0.,1.},0);
}
void controls()
{
    const Floating nearest;const long double tiny=std::numeric_limits<double>::denorm_min(),maximum=std::numeric_limits<double>::max();
    const auto exercise=[&]{productTape({1.25L,-1.5L,1.75L});productTape({tiny,1.5L});productTape({-0.,std::numeric_limits<double>::quiet_NaN()});
        sumTape({maximum,maximum,maximum,maximum},1,{maximum,maximum,maximum,maximum},-1);
        sumTape({tiny,1.5L},1,{1.L,1.L},-1);};
    for(int mode:{FE_UPWARD,FE_DOWNWARD,FE_TOWARDZERO}){check(std::fesetround(mode)==0&&std::fegetround()==mode,"small exact cannot establish real directed controls");exercise();nearest.restore();}
    for(unsigned mode:{1u<<15,1u<<6,(1u<<15)|(1u<<6)}){_mm_setcsr(_mm_getcsr()|mode);exercise();nearest.restore();}
    check(std::fesetround(FE_UPWARD)==0,"small exact cannot establish combined unsupported controls");_mm_setcsr(_mm_getcsr()|(1u<<15));exercise();nearest.restore();
}
void invalidCounts()
{
    const auto q=std::numeric_limits<double>::quiet_NaN();
    const auto checkCount=[&](std::initializer_list<long double> factors){
        DyadicTerminalRneContract::prime(entryVariant++);const Floating before;
        const auto actual=NativeSmallProductContract::invoke<Small::Product>([&]{return Small::product(factors);});const Floating after;
        check(!actual.value&&actual.error=="[FiniteVolume] native small exact product requires supported two-to-four-factor format",
            "small exact own count guard failed before scalar consumption");NativeSmallProductContract::sameFp(before,after);++countGuards;};
    checkCount({});checkCount({q});checkCount({q,1.,2.,3.,4.});
}
void run()
{
    if constexpr(!Small::supported){std::cout<<"NATIVE_SMALL_EXACT_UNAVAILABLE,native_binary64_format\n";return;}
    NativeSmallProductContract::Restore restore;check(std::fesetround(FE_TONEAREST)==0,"small exact requires nearest baseline");
    _mm_setcsr(_mm_getcsr()&~((1u<<15)|(1u<<6)));unsigned native=0;
    check(_controlfp_s(&native,0,0)==0&&(native&_MCW_EM)==_MCW_EM&&(_mm_getcsr()&0x1f80u)==0x1f80u,"small exact test requires masked controls; traps are not injected");
    const auto fpBefore=NativeSmallProductContract::fpPairs;
    productsAndFirstErrors();sumsAndBoundaries();controls();invalidCounts();
    check(tapes>=90&&productPairs>tapes&&combinePairs>30&&terminalPairs>75&&literalTerminalPairs>50&&independentBits>50
        &&exactReferences>50&&firstErrors>=6&&firstZeros>=6&&additionGuards==3&&countGuards==3&&controlGuards>=30
        &&negativeZeros>0&&canonicalZeros>0&&terminalOverflows>0&&hugeCancellations>=2&&copyMovePairs>=12
        &&maximumProductBits>=212&&maximumAlignedSumBits>=8392
        &&NativeSmallProductContract::fpPairs-fpBefore==productPairs+combinePairs+terminalPairs+literalTerminalPairs+countGuards,
        "small exact original-prefix/exact-integer/RNE/FP/boundary coverage incomplete");
    std::cout<<"NATIVE_SMALL_EXACT_PASS,tapes,"<<tapes<<",original_product_prefix_pairs,"<<productPairs<<",original_two_add_pairs,"<<combinePairs
        <<",original_Dyadic_terminal_pairs,"<<terminalPairs<<",original_UInt_terminal_pairs,"<<literalTerminalPairs
        <<",full_native_rawCSR_errno_pairs,"<<(NativeSmallProductContract::fpPairs-fpBefore)<<",independent_IEEE_bits,"<<independentBits
        <<",independent_exact_reference_values,"<<exactReferences<<",maximum_product_significand_bits,"<<maximumProductBits
        <<",maximum_actual_aligned_sum_bits,"<<maximumAlignedSumBits<<",huge_exact_cancellations_without_product_rounding,"<<hugeCancellations
        <<",copy_move_terminal_pairs,"<<copyMovePairs<<",original_nonfinite_first_errors,"<<firstErrors<<",original_first_zero_tapes,"<<firstZeros
        <<",original_sign_guards,"<<additionGuards<<",new_count_guards,"<<countGuards<<",original_terminal_control_guards,"<<controlGuards
        <<",negative_rounded_zeros,"<<negativeZeros<<",canonical_exact_zeros,"<<canonicalZeros<<",original_terminal_overflows,"<<terminalOverflows
        <<",no_private_limb_oracle,no_trap_or_resource_injection,no_speed_claim\n";
}
#else
void run(){std::cout<<"NATIVE_SMALL_EXACT_UNAVAILABLE,unsupported_compile_platform\n";}
#endif
}

namespace NativeScalarPrefixContract
{
namespace FV=LAMP::Geometry::Detail::FiniteVolume;
#if defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE) && !defined(__FAST_MATH__) && (!defined(_M_FP_FAST) || !_M_FP_FAST)
using Wide=NativeSmallProductContract::Wide;
using Floating=NativeSmallProductContract::Floating;
struct Part {std::uint64_t magnitude=0;int exponent=0;bool negative=false;};
struct Tape {std::array<Part,4> parts{};unsigned count=0;bool zero=false;};
std::size_t pairs=0,firstErrors=0,zeros=0,outputs=0,controlModes=0;
unsigned entryVariant=0;
Part literalScalar(long double value)
{
    // The unchanged scalar FP prefix, independently retained for this test.
    // No Dyadic integer representation, normalization or bit decoder is copied.
    FV::Internal::require(std::isfinite(value),"[FiniteVolume] nonfinite dyadic input");
    FV::Internal::require(std::numeric_limits<long double>::radix==2,"[FiniteVolume] binary native arithmetic required");
    Part result;if(value==0)return result;result.negative=std::signbit(value);int exponent=0;
    const auto fraction=std::frexp(std::abs(value),&exponent);
    const auto integer=std::scalbn(fraction,std::numeric_limits<long double>::digits);
    FV::Internal::require(integer>=1&&integer<std::scalbn(1.L,64),"[FiniteVolume] native significand outside uint64");
    result.magnitude=static_cast<std::uint64_t>(integer);
    result.exponent=exponent-std::numeric_limits<long double>::digits;return result;
}
Tape literal(std::initializer_list<long double> factors,bool initialOne)
{
    if(initialOne)(void)literalScalar(1.L); // Both product entries begin with the original scalar one.
    Tape result;for(const auto factor:factors){const auto part=literalScalar(factor);
        if(part.magnitude==0){result.zero=true;return result;}
        result.parts[result.count++]=part;}return result;
}
Wide exact(const Tape& tape)
{
    Wide value;if(tape.zero){value.magnitude=0;return value;}
    for(unsigned i=0;i<tape.count;++i){const auto& part=tape.parts[i];value.magnitude*=part.magnitude;
        value.exponent+=part.exponent;value.negative=value.negative!=part.negative;}return value;
}
template<class T,class Make>void paired(std::initializer_list<long double> factors,bool initialOne,Make make)
{
    DyadicTerminalRneContract::prime(entryVariant++);const Floating entry;
    const auto actual=NativeSmallProductContract::invoke<T>(make);const Floating actualAfter;
    entry.restore();const auto original=NativeSmallProductContract::invoke<Tape>([&]{return literal(factors,initialOne);});
    const Floating originalAfter;NativeSmallProductContract::sameFp(actualAfter,originalAfter);entry.restore();++pairs;
    check(actual.value.has_value()==original.value.has_value()&&actual.error==original.error,
        "normal scalar decoding changed original literal prefix/first error");
    if(!actual.value){check(actual.error=="[FiniteVolume] nonfinite dyadic input","normal scalar nonfinite first error changed");++firstErrors;return;}
    const auto reference=exact(*original.value);if(reference.magnitude==0)++zeros;
    if(std::fegetround()==FE_TONEAREST&&NativeSmallProductContract::gradualControls()){
        const auto independent=NativeSmallProductContract::independent(factors);
        check(independent.valid&&DyadicIntegerOperationsContract::exact(independent)==DyadicIntegerOperationsContract::exact(reference),
            "literal scalar prefix disagrees with independent original IEEE inputs");}
    // Construction is already checked under its actual controls. Observe the
    // resulting exact value under a separate nearest/gradual terminal entry.
    const Floating construction;check(std::fesetround(FE_TONEAREST)==0,"scalar oracle cannot establish observation rounding");
    _mm_setcsr(_mm_getcsr()&~((1u<<15)|(1u<<6)));
    if constexpr(std::is_same_v<T,FV::Dyadic>)
        DyadicIntegerOperationsContract::observe({*actual.value,DyadicIntegerOperationsContract::exact(reference)});
    else{
        const auto wanted=NativeSmallProductContract::ieeeRne(reference);
        const auto got=NativeSmallProductContract::invoke<double>([&]{return actual.value->roundedDouble();});
        if(wanted)check(got.value&&NativeSmallProductContract::bits(*got.value)==*wanted,"normal scalar exact product changed independent IEEE bits");
        else check(!got.value&&got.error=="[FiniteVolume] final binary64 dyadic overflow","normal scalar product changed original terminal overflow");}
    construction.restore();++outputs;
}
void one(long double value)
{
    paired<FV::Dyadic>({value},false,[&]{return FV::Dyadic(value);});
    paired<FV::Dyadic>({1.L,value,1.L},true,[&]{return FV::Dyadic::product({1.L,value,1.L});});
    paired<FV::NativeSmallExact::Product>({1.L,1.L,value,1.L},true,[&]{return FV::NativeSmallExact::product({1.L,1.L,value,1.L});});
}
void nonfiniteOrder()
{
    const auto q=std::numeric_limits<long double>::quiet_NaN(),inf=std::numeric_limits<long double>::infinity();
    const long double signal=NativeSmallProductContract::fromBits(UINT64_C(0x7ff0000000000001));
    const auto product=[&](std::initializer_list<long double> factors){
        paired<FV::Dyadic>(factors,true,[&]{return FV::Dyadic::product(factors);});
        paired<FV::NativeSmallExact::Product>(factors,true,[&]{return FV::NativeSmallExact::product(factors);});};
    product({q,0.L});product({0.L,q});product({signal,0.L});product({-0.L,signal});
    product({1.L,inf,0.L});product({1.L,0.L,inf});
    product({1.L+std::ldexp(1.L,-52),1.L+std::ldexp(1.L,-52)});
}
void run()
{
    NativeSmallProductContract::Restore restore;
    check(FV::NativeSmallExact::supported,"scalar prefix oracle requires actual native binary64 format");
    check(std::fesetround(FE_TONEAREST)==0,"scalar prefix oracle cannot establish nearest baseline");
    _mm_setcsr(_mm_getcsr()&~((1u<<15)|(1u<<6)));unsigned control=0;
    check(_controlfp_s(&control,0,0)==0&&(control&_MCW_EM)==_MCW_EM&&(_mm_getcsr()&0x1f80u)==0x1f80u,
        "scalar prefix oracle requires masked controls; traps are not injected");
    const auto fpBefore=NativeSmallProductContract::fpPairs;
    for(auto word:{UINT64_C(0x0010000000000000),UINT64_C(0x0010000000000001),UINT64_C(0x3fefffffffffffff),
            UINT64_C(0x3ff0000000000000),UINT64_C(0x3ff0000000000001),UINT64_C(0x3ff4000000000000),
            UINT64_C(0x3ff8000000000000),UINT64_C(0x4020000000000000),UINT64_C(0x7fefffffffffffff)})
        for(auto sign:{UINT64_C(0),UINT64_C(0x8000000000000000)})one(NativeSmallProductContract::fromBits(word|sign));
    for(auto word:{UINT64_C(0),UINT64_C(0x8000000000000000),UINT64_C(1),UINT64_C(0x8000000000000001),
            UINT64_C(0x000fffffffffffff),UINT64_C(0x7ff0000000000000),UINT64_C(0x7ff8000000000001),UINT64_C(0x7ff0000000000001)})
        one(NativeSmallProductContract::fromBits(word));
    nonfiniteOrder();const Floating baseline;
    const auto exercise=[&]{one(std::numeric_limits<long double>::min());one(std::numeric_limits<long double>::max());
        one(std::numeric_limits<long double>::denorm_min());one(-0.L);
        one(NativeSmallProductContract::fromBits(UINT64_C(0x7ff0000000000001)));++controlModes;};
    for(int mode:{FE_UPWARD,FE_DOWNWARD,FE_TOWARDZERO}){check(std::fesetround(mode)==0,"scalar prefix directed control setup");exercise();baseline.restore();}
    for(unsigned mode:{1u<<15,1u<<6,(1u<<15)|(1u<<6)}){_mm_setcsr(_mm_getcsr()|mode);exercise();baseline.restore();}
    check(std::fesetround(FE_UPWARD)==0,"scalar prefix combined control setup");_mm_setcsr(_mm_getcsr()|(1u<<15));exercise();baseline.restore();
    check(pairs>180&&firstErrors>=30&&zeros>20&&outputs>100&&controlModes==7
        &&NativeSmallProductContract::fpPairs-fpBefore==pairs,"scalar prefix normal/non-normal/first-error/FP coverage incomplete");
    std::cout<<"NATIVE_SCALAR_PREFIX_PASS,literal_prefix_pairs,"<<pairs<<",actual_entry_points,3,independent_exact_output_observations,"<<outputs
        <<",original_nonfinite_first_errors,"<<firstErrors<<",first_zero_results,"<<zeros<<",actual_control_modes,"<<controlModes
        <<",full_native_rawCSR_errno_pairs,"<<(NativeSmallProductContract::fpPairs-fpBefore)
        <<",original_frexp_scalbn_prefix_reference,1,no_Dyadic_algorithm_copy,no_trap_or_resource_injection,no_speed_claim\n";
}
#else
void run(){std::cout<<"NATIVE_SCALAR_PREFIX_UNAVAILABLE,unsupported_compile_platform\n";}
#endif
}
namespace StationaryNoSlipStrainContract
{
using namespace LAMP::Hydrodynamics;
std::size_t parallelComparisons=0,literalComparisons=0,fpComparisons=0,multiwallGuards=0,coefficients=0;
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
using Floating=ViscosityProductParallelContract::Floating;
void sameFp(const Floating& a,const Floating& b){
    check(std::memcmp(&a.native,&b.native,sizeof a.native)==0&&a.csr==b.csr&&a.error==b.error,
        "no-slip fit changed original complete native FP, raw MXCSR or errno");++fpComparisons;
}
void sameRaw(const HydroStrainOperator& actual,const HydroStrainOperator& expected){
    ViscosityProductPlanContract::exactMatrix(actual.matrix(),expected.matrix(),coefficients);
    check(actual.components()==expected.components()&&actual.stationaryNoSlipVelocityDofs()==expected.stationaryNoSlipVelocityDofs()
        &&actual.points().size()==expected.points().size()&&!actual.isStressCompatible()&&!expected.isStressCompatible(),
        "no-slip thread request changed raw component/point/wall capability");
    for(std::size_t q=0;q<actual.points().size();++q){const auto& a=actual.points()[q];const auto& b=expected.points()[q];
        check(a.cellDof==b.cellDof&&std::memcmp(&a.point.x,&b.point.x,sizeof(double))==0
            &&std::memcmp(&a.point.y,&b.point.y,sizeof(double))==0&&std::memcmp(&a.physicalWeight,&b.physicalWeight,sizeof(double))==0,
            "no-slip fit changed original quadrature cell, point or physical weight bits");}
}
std::pair<HydroStrainOperator,HydroStrainOperator> compare(const QuadtreeMesh& mesh,const HydroPhysicalDivergence& metrics,
    const std::vector<HydroStrainPoint>& points,const std::vector<int>& walls){
    ViscosityProductParallelContract::RestoreRuntime runtime;omp_set_dynamic(0);
    check(std::fesetround(FE_TONEAREST)==0,"no-slip fit comparison cannot set RNE");_mm_setcsr(_mm_getcsr()&~0x8040u);
    const Floating entry=ViscosityProductParallelContract::cleanInput();
    std::optional<HydroStrainOperator> literal;std::optional<Floating> literalAfter;std::exception_ptr failure;bool active=false;
#pragma omp parallel num_threads(2)
    {
#pragma omp single
        {try{ViscosityProductParallelContract::RestoreFloating worker;entry.restore();active=omp_in_parallel()!=0;
            literal.emplace(HydroStrainOperator::buildStationaryNoSlip(mesh,metrics,points,walls));literalAfter.emplace();}
        catch(...){failure=std::current_exception();}}
    }
    if(failure)std::rethrow_exception(failure);check(active&&literal&&literalAfter,"no-slip original nested literal route was not actually active");
    omp_set_num_threads(1);entry.restore();const auto one=HydroStrainOperator::buildStationaryNoSlip(mesh,metrics,points,walls);const Floating afterOne;
    omp_set_num_threads(8);entry.restore();auto eight=HydroStrainOperator::buildStationaryNoSlip(mesh,metrics,points,walls);const Floating afterEight;
    sameFp(*literalAfter,afterOne);sameFp(afterOne,afterEight);sameRaw(one,*literal);sameRaw(eight,*literal);
    ++literalComparisons;++parallelComparisons;return {std::move(eight),std::move(*literal)};
}
void multiwall(const QuadtreeMesh& mesh,const HydroPhysicalDivergence& metrics,const std::vector<HydroStrainPoint>& points,
    const std::vector<int>& walls){
    const auto& layout=metrics.layout();const auto& domain=layout.domain();std::vector<int> corner;
    for(int w:walls){const auto a=domain.face(mesh,layout.nativeFace(w));const int owner=a.negativeCell>=0?a.negativeCell:a.positiveCell;
        for(int u=0;u<layout.velocityCount();++u){const auto b=domain.face(mesh,layout.nativeFace(u));
            if(b.outwardSign!=0&&b.area>0&&b.horizontal!=a.horizontal&&(b.negativeCell>=0?b.negativeCell:b.positiveCell)==owner){corner={w,u};break;}}
        if(!corner.empty())break;}
    check(corner.size()==2,"existing no-slip fixture lacks two real perpendicular positive-area walls on one active owner");
    ViscosityProductParallelContract::RestoreRuntime runtime;omp_set_dynamic(0);
    const Floating entry=ViscosityProductParallelContract::cleanInput();
    std::string literalError;std::optional<Floating> literalAfter;std::exception_ptr failure;bool active=false;
#pragma omp parallel num_threads(2)
    {
#pragma omp single
        {try{ViscosityProductParallelContract::RestoreFloating worker;entry.restore();active=omp_in_parallel()!=0;
            try{(void)HydroStrainOperator::buildStationaryNoSlip(mesh,metrics,points,corner);}
            catch(const std::invalid_argument& e){literalAfter.emplace();literalError=e.what();}}
        catch(...){failure=std::current_exception();}}
    }
    if(failure)std::rethrow_exception(failure);
    check(active&&literalAfter&&literalError.find("genuine multiwall corner")!=std::string::npos,
        "original no-slip literal route failed to reject the actual two-wall corner");
    for(int team:{1,8}){omp_set_num_threads(team);entry.restore();std::string message;std::optional<Floating> after;
        try{(void)HydroStrainOperator::buildStationaryNoSlip(mesh,metrics,points,corner);}
        catch(const std::invalid_argument& e){after.emplace();message=e.what();}
        check(after&&message==literalError,"no-slip serial/parallel fit changed earliest original physical corner exception");sameFp(*literalAfter,*after);}
    ++multiwallGuards;
}
#else
std::pair<HydroStrainOperator,HydroStrainOperator> compare(const QuadtreeMesh& mesh,const HydroPhysicalDivergence& metrics,
    const std::vector<HydroStrainPoint>& points,const std::vector<int>& walls){
    auto original=HydroStrainOperator::buildStationaryNoSlip(mesh,metrics,points,walls);return {original,original};
}
void multiwall(const QuadtreeMesh&,const HydroPhysicalDivergence&,const std::vector<HydroStrainPoint>&,const std::vector<int>&){}
#endif
void run()
{
    int scenes=0,guards=0,wallChecks=0;double peakStrain=0,peakDissipation=0,smallestPlug=std::numeric_limits<double>::infinity();
    for(bool rz:{false,true})for(bool amr:{false,true})for(bool horizontalWall:{false,true}){
        DictionaryMesh d;d.Nx=12;d.Ny=16;d.minCoord={0,0};d.maxCoord={1.5,2};d.maxInitialIterations=0;d.maxCellNum=10000;
        QuadtreeMesh mesh(d);
        if(amr){DictionaryMesh::RefineBox box;box.min={.75,.5};box.max={1.25,1.5};box.level=1;mesh.applyRefineBoxes({box});}
        FaceSegmentGraph graph;graph.build(mesh,[](const QuadtreeMesh&,const Cell& c){return c.centerCoord().x<1;});
        const auto domain=NativeLeafDomain::build(mesh,graph,rz);const auto layout=HydroDofLayout::build(mesh,domain);
        const auto metrics=HydroPhysicalDivergence::build(mesh,HydroIntegralDivergence::build(mesh,layout));
        const auto points=HydroStrainOperator::volumeQuadrature(mesh,metrics);std::vector<int> walls;
        int interior=-1,axis=-1;
        for(int u=0;u<layout.velocityCount();++u){const auto f=domain.face(mesh,layout.nativeFace(u));
            if(f.outwardSign!=0&&(horizontalWall?(f.horizontal&&f.center.y==2):(!f.horizontal&&f.center.x==1)))walls.push_back(u);
            if(f.outwardSign==0)interior=u;if(f.area==0)axis=u;}
        check(!walls.empty(),"no actual internal stationary wall");
        std::cerr<<"NO_SLIP_SCENE rz="<<rz<<" amr="<<amr<<" horizontal_wall="<<horizontalWall<<" building\n"<<std::flush;
        const auto paired=compare(mesh,metrics,points,walls);const auto& e=paired.first;const auto& literal=paired.second;
        std::cerr<<"NO_SLIP_SCENE built\n"<<std::flush;
        check(e.stationaryNoSlipVelocityDofs()==walls&&!e.isStressCompatible(),"stationary wall lost raw strain owner");
        std::vector<double> u(layout.velocityCount()),plug(u.size());
        for(int j=0;j<layout.velocityCount();++j){const auto f=domain.face(mesh,layout.nativeFace(j));
            if(horizontalWall){u[j]=f.horizontal?0:(rz?f.center.x:1)*(2-f.center.y);continue;}
            if(!f.horizontal)continue;
            const auto& edge=mesh.hEdges(f.edgeLabel);const auto a=mesh.verts(edge.edgeToVert[0]).coord(),b=mesh.verts(edge.edgeToVert[1]).coord();
            const double length=std::abs(b.x-a.x);
            // Independent physical face integral of 1-r^2, with r dr in RZ.
            u[j]=1-f.center.x*f.center.x-length*length/(rz?4.:12.);plug[j]=1;}
        const HydroNativeValues values(layout,u,std::vector<double>(layout.pressureCount()));const auto strain=e.apply(mesh,values);
        ViscosityProductPlanContract::exactWeights(strain,literal.apply(mesh,values));
        for(std::size_t q=0;q<points.size();++q)for(int k=0;k<e.components();++k){
            double exact=k==e.components()-1?-std::sqrt(2.)*points[q].point.x:0;
            if(horizontalWall){exact=0;if(rz&&(k==0||k==2))exact=2-points[q].point.y;
                if(k==e.components()-1)exact=-(rz?points[q].point.x:1)/std::sqrt(2.);}

            const double error=std::abs(strain[q*e.components()+k]-exact);peakStrain=std::max(peakStrain,error);
            check(error<=1e-11,"stationary no-slip Poiseuille physical-mean P2 strain mismatch");}
        const double mu=.3;const std::vector<double> viscositySamples(points.size(),mu);const auto viscosity=HydroEnergyViscosity::build(e,viscositySamples);
        const auto oldProduct=ViscosityProductPlanContract::original(literal,viscositySamples);
        ViscosityProductPlanContract::exactWeights(viscosity.rowWeights(),oldProduct.weights);
        ViscosityProductPlanContract::exactMatrix(viscosity.matrix(),oldProduct.matrix,coefficients);
        check(!viscosity.productReport().noSlip.available,"public Energy build repeated a historical raw no-slip preparation receipt");
        const double energy=viscosity.dissipation(mesh,values),expected=horizontalWall?(rz?35*std::acos(-1.)*mu/3:2*mu):(rz?4*std::acos(-1.)*mu:8*mu/3);
        const double relative=std::abs(energy-expected)/expected;peakDissipation=std::max(peakDissipation,relative);
        check(relative<=1e-11,"stationary no-slip physical Poiseuille dissipation mismatch");
        if(!horizontalWall){const HydroNativeValues plugValues(layout,plug,std::vector<double>(layout.pressureCount()));
        const double resistance=viscosity.dissipation(mesh,plugValues);smallestPlug=std::min(smallestPlug,resistance);
        check(std::isfinite(resistance)&&resistance>mu,"stationary no-slip accepted axial plug as zero strain");
        const auto freeE=HydroStrainOperator::build(mesh,metrics,points);const auto freeK=HydroEnergyViscosity::build(freeE,std::vector<double>(points.size(),mu));
        check(freeK.dissipation(mesh,plugValues)<1e-20,"original free-slip plug comparison is not a null mode");}
        const auto applied=viscosity.apply(mesh,values);long double actual=0,gross=0;
        for(std::size_t j=0;j<u.size();++j){actual+=static_cast<long double>(u[j])*applied[j];gross+=std::abs(static_cast<long double>(u[j])*applied[j]);}
        check(std::abs(actual-energy)<=1e-11L*gross,"no-slip stored Gram work mismatch");
        for(int j:walls){check(u[j]==0&&plug[j]==0,"stationary wall normal is nonzero");++wallChecks;}
        const auto bad=[&](const auto& rows){bool rejected=false;try{HydroStrainOperator::buildStationaryNoSlip(mesh,metrics,points,rows);}catch(const std::invalid_argument&){rejected=true;}
            check(rejected,"invalid no-slip rows accepted");++guards;};
        bad(std::vector<int>{});bad(std::vector<int>{walls[0],walls[0]});bad(std::vector<int>{interior});bad(std::vector<int>{layout.velocityCount()});if(rz)bad(std::vector<int>{axis});
        multiwall(mesh,metrics,points,walls);++scenes;
    }
    check(scenes==8&&guards==36&&wallChecks>0,"no-slip XY/RZ/2:1/solid coverage incomplete");
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
    check(parallelComparisons==8&&literalComparisons==8&&multiwallGuards==8&&fpComparisons==32,
        "no-slip actual nested literal/request1/request8 or earliest-corner FP comparisons incomplete");
#endif
    std::cout<<"STATIONARY_NO_SLIP_STRAIN_PASS scenes="<<scenes<<" guards="<<guards<<" actual_internal_wall_rows="<<wallChecks
        <<" maximum_physical_mean_strain_error="<<peakStrain<<" maximum_dissipation_relative_error="<<peakDissipation
        <<" minimum_plug_dissipation="<<smallestPlug<<" requested1_8_pairs="<<parallelComparisons<<" active_nested_literal_pairs="<<literalComparisons
        <<" exact_E_K_coefficients="<<coefficients<<" full_native_FP_rawCSR_errno_pairs="<<fpComparisons<<" real_multiwall_guards="<<multiwallGuards
        <<" original36guards_and_physical_gates_retained_no_internal_worker_count_claim\n";
}
}

// BEGIN HELMHOLTZ_SOURCE_RESEARCH_PHASE1: test-only law, no production producer.
namespace HelmholtzSourceResearch {
namespace E=LAMP::Electric;namespace I=LAMP::Interface;
namespace FV=LAMP::Geometry::Detail::FiniteVolume;
namespace Gas=MaxwellPureGasSolenoidalContract;
using Clock=LAMP::TimeLayer;using Matrix=HydroStrainOperator::Matrix;
constexpr long double originalWorkTolerance=1e-12L;
std::size_t slabCases=0,counterexamples=0,gasProbes=0;
bool physicalReferences=true,counterexampleIdentities=true;
double peakReference=0,peakCounterexampleIdentity=0;
FV::Dyadic phi(const E::ElectricVofField& f,int c) {
    FV::Dyadic value;for(double term:f.potentialTerms()[c])value.add(term);return value;
}
bool reference(long double actual,long double wanted,long double scale,double& peak) {
    const long double error=std::abs(actual-wanted),relative=scale==0?(error==0?0:std::numeric_limits<long double>::infinity()):error/scale;
    peak=std::max(peak,static_cast<double>(relative));return std::isfinite(relative)&&error<=originalWorkTolerance*scale;
}
struct Work {
    long double charge=0,volume=0,material=0,adjointNative=0,alternativeNative=0;
    long double adjointRelative=0,alternativeRelative=0,adjointScale=0,alternativeScale=0;
};
struct ResearchLaw {
    struct ChargeFace {int column=-1,negative=-1,positive=-1;double area=0,outward=0,boundaryPhi=0;
        long double negativeWeightByVolume=0,positiveWeightByVolume=0;};
    struct MaterialFace {int column=-1;FV::Dyadic force;};
    const E::ElectricVofField& field;const HydroStrainOperator& geometry;const E::ElectricVofMaxwellStress& stress;
    std::vector<ChargeFace> chargeFaces;std::vector<MaterialFace> materialFaces;
    std::vector<FV::Dyadic> adjointForce,volumeForce,materialForce;
    std::vector<std::array<FV::Dyadic,2>> volumeFactors;
    std::vector<double> adjointRhs,alternativeRhs;
    ResearchLaw(const QuadtreeMesh& mesh,const HydroMaxwellLoad& load,const std::vector<E::ElectricVofBoundary>& actualFixtureBoundary)
        :field(load.stress().field()),geometry(load.geometry()),stress(load.stress()) {
        const auto& layout=geometry.layout();const auto& domain=layout.domain();const int nu=layout.velocityCount();
        load.validateFor(mesh,load.pressureLoad(),load.phaseGeometry(),stress,field.candidateTime());
        check(geometry.hasCanonicalRawStrain()&&actualFixtureBoundary.size()==field.faces().size(),"research lost actual raw geometry/explicit fixture boundary");
        adjointForce.resize(nu);volumeForce.resize(nu);materialForce.resize(nu);
        const auto& sourceQ=field.chargeCandidate().integratedCharge();
        const auto volume=[&](int c){const int row=domain.cellRowByIndex()[c];check(row>=0,"research Q source outside actual domain");return domain.cellVolumes()[row];};
        std::map<int,std::vector<int>> owners;
        for(int row=0;row<domain.nativeFaceCount();++row){const auto f=domain.face(mesh,row);if(f.negativeCell<0&&f.positiveCell<0)continue;
            const auto& edge=f.horizontal?mesh.hEdges(f.edgeLabel):mesh.vEdges(f.edgeLabel);
            for(int v:edge.edgeToVert)for(int c:{f.negativeCell,f.positiveCell})if(c>=0)owners[v].push_back(c);}
        std::map<int,FV::Dyadic> vertexPhi;
        const auto projected=[&](int c,Vector2D v){auto value=phi(field,c);const auto center=mesh.cells(mesh.cellIndexToLabel[c]).centerCoord();
            const auto grad=field.cellPotentialGradients()[c];value.add(FV::Dyadic::product({grad.x,static_cast<long double>(v.x)-center.x}));
            value.add(FV::Dyadic::product({grad.y,static_cast<long double>(v.y)-center.y}));return value;};
        for(auto& item:owners){auto& cells=item.second;std::sort(cells.begin(),cells.end());cells.erase(std::unique(cells.begin(),cells.end()),cells.end());
            const auto point=mesh.verts(item.first).coord();const auto base=projected(cells.front(),point);auto value=base;
            const long double weight=1.L/cells.size();for(std::size_t k=1;k<cells.size();++k){auto delta=projected(cells[k],point);delta.add(base,-1);value.add(delta.multiplied(weight));}
            vertexPhi.emplace(item.first,std::move(value));}
        for(int row=0;row<domain.nativeFaceCount();++row){const auto f=domain.face(mesh,row);if(f.negativeCell<0&&f.positiveCell<0)continue;
            const int u=layout.velocityDof(row);check(u>=0,"research native face lacks same Hydro unknown");if(f.area==0)continue;
            ChargeFace entry;entry.column=u;entry.negative=f.negativeCell;entry.positive=f.positiveCell;entry.area=f.area;
            const bool internal=f.negativeCell>=0&&f.positiveCell>=0;
            if(internal){const auto a=mesh.cells(mesh.cellIndexToLabel[f.negativeCell]).centerCoord(),b=mesh.cells(mesh.cellIndexToLabel[f.positiveCell]).centerCoord();
                const long double da=f.horizontal?static_cast<long double>(f.center.y)-a.y:static_cast<long double>(f.center.x)-a.x;
                const long double db=f.horizontal?static_cast<long double>(b.y)-f.center.y:static_cast<long double>(b.x)-f.center.x;
                check(da>0&&db>0,"research virtual Q weights have foreign face distances");
                entry.negativeWeightByVolume=(db/(da+db))/volume(f.negativeCell);entry.positiveWeightByVolume=(da/(da+db))/volume(f.positiveCell);
            }else {const int c=f.negativeCell>=0?f.negativeCell:f.positiveCell;const auto& boundary=actualFixtureBoundary[row];
                check(boundary.type==E::ElectricVofBoundary::Type::Dirichlet||sourceQ[c]==0,"research lacks explicit nonzero-Q boundary potential");
                entry.outward=f.negativeCell>=0?1.:-1.;entry.boundaryPhi=boundary.value;
                if(f.negativeCell>=0)entry.negativeWeightByVolume=1.L/volume(c);else entry.positiveWeightByVolume=1.L/volume(c);}
            chargeFaces.push_back(entry);
            const auto addCharge=[&](int c,long double weight){if(c<0)return;auto base=FV::Dyadic::product({entry.area,weight,sourceQ[c]});
                if(internal){for(double term:field.potentialTerms()[f.negativeCell])adjointForce[u].add(base.multiplied(term));
                    for(double term:field.potentialTerms()[f.positiveCell])adjointForce[u].add(base.multiplied(term),-1);}
                else {base=base.multiplied(entry.outward);for(double term:field.potentialTerms()[c])adjointForce[u].add(base.multiplied(term));
                    adjointForce[u].add(base.multiplied(entry.boundaryPhi),-1);}};
            addCharge(entry.negative,entry.negativeWeightByVolume);addCharge(entry.positive,entry.positiveWeightByVolume);
            if(!internal)continue;
            const auto& edge=f.horizontal?mesh.hEdges(f.edgeLabel):mesh.vEdges(f.edgeLabel);int va=edge.edgeToVert[0],vb=edge.edgeToVert[1];
            auto a=mesh.verts(va).coord(),b=mesh.verts(vb).coord();if(f.horizontal?b.x<a.x:b.y<a.y){std::swap(va,vb);std::swap(a,b);}
            const long double length=f.horizontal?static_cast<long double>(b.x)-a.x:static_cast<long double>(b.y)-a.y;
            auto integral=vertexPhi.at(va);integral.add(vertexPhi.at(vb),-1);
            const long double et=integral.rounded()/length,dn=static_cast<long double>(field.faces()[row].displacement)/f.area;
            const long double em=field.epsilonByCell()[f.negativeCell],ep=field.epsilonByCell()[f.positiveCell],enMinus=dn/em,enPlus=dn/ep;
            check(std::isfinite(et)&&std::isfinite(enMinus)&&std::isfinite(enPlus),"research shared field trace is unrepresentable");
            auto secant=FV::Dyadic::product({enMinus,enPlus});secant.add(FV::Dyadic::product({et,et}));
            auto force=secant.multiplied(-.5L).multiplied(ep-em).multiplied(f.area);
            materialForce[u].add(force);materialFaces.push_back({u,std::move(force)});
        }
        const auto& raw=geometry.stressVolumeVelocityMatrix();volumeFactors.resize(geometry.points().size());
        for(std::size_t q=0;q<geometry.points().size();++q){const auto& point=geometry.points()[q];const int c=layout.nativeCell(point.cellDof);
            const auto sample=stress.sample(c,point.point);const long double fraction=static_cast<long double>(point.physicalWeight)/volume(c);
            for(int axis=0;axis<2;++axis){auto factor=FV::Dyadic::product({sourceQ[c],fraction,axis==0?sample.electricField.x:sample.electricField.y});
                volumeFactors[q][axis]=factor;for(Matrix::InnerIterator e(raw,static_cast<Eigen::Index>(q)*2+axis);e;++e)volumeForce[e.col()].add(factor.multiplied(e.value()));}}
        adjointRhs.resize(nu);alternativeRhs.resize(nu);for(int u=0;u<nu;++u){auto a=adjointForce[u];a.add(materialForce[u]);adjointRhs[u]=a.roundedDouble();
            auto b=volumeForce[u];b.add(materialForce[u]);alternativeRhs[u]=b.roundedDouble();}
    }
    Work work(const std::vector<double>& u) const {
        check(u.size()==adjointRhs.size(),"research virtual velocity has foreign native layout");
        const auto& q=field.chargeCandidate().integratedCharge();std::vector<FV::Dyadic> dotQ(q.size());
        FV::Dyadic boundaryWork,chargeWork,chargeScale,volumeWork,volumeScale,materialWork,materialScale,totalOutward;
        for(const auto& f:chargeFaces){FV::Dyadic flux;
            if(f.negative>=0)flux.add(FV::Dyadic::product({f.area,f.negativeWeightByVolume,q[f.negative],u[f.column]}));
            if(f.positive>=0)flux.add(FV::Dyadic::product({f.area,f.positiveWeightByVolume,q[f.positive],u[f.column]}));
            if(f.outward==0){dotQ[f.negative].add(flux,-1);dotQ[f.positive].add(flux);}
            else {const auto outward=flux.multiplied(f.outward);totalOutward.add(outward);const int c=f.negative>=0?f.negative:f.positive;dotQ[c].add(outward,-1);
                const auto term=outward.multiplied(-f.boundaryPhi);boundaryWork.add(term);chargeScale.add(term.absolute());}}
        auto virtualConservation=totalOutward;for(const auto& value:dotQ)virtualConservation.add(value);
        check(virtualConservation.compare(FV::Dyadic{})==0,"research shared virtual-Q flux lost exact conservation/gauge identity");
        for(std::size_t c=0;c<dotQ.size();++c)for(double term:field.potentialTerms()[c]){const auto value=dotQ[c].multiplied(-term);chargeWork.add(value);chargeScale.add(value.absolute());}
        chargeWork.add(boundaryWork);
        const auto& raw=geometry.stressVolumeVelocityMatrix();
        for(std::size_t qpoint=0;qpoint<volumeFactors.size();++qpoint)for(int axis=0;axis<2;++axis)
            for(Matrix::InnerIterator e(raw,static_cast<Eigen::Index>(qpoint)*2+axis);e;++e){const auto term=volumeFactors[qpoint][axis].multiplied(e.value()).multiplied(u[e.col()]);volumeWork.add(term);volumeScale.add(term.absolute());}
        for(const auto& f:materialFaces){const auto term=f.force.multiplied(u[f.column]);materialWork.add(term);materialScale.add(term.absolute());}
        FV::Dyadic actualAdjoint,actualAlternative,adjointScale=chargeScale,alternativeScale=volumeScale;
        adjointScale.add(materialScale);alternativeScale.add(materialScale);
        for(std::size_t j=0;j<u.size();++j){const auto a=FV::Dyadic::product({u[j],adjointRhs[j]}),b=FV::Dyadic::product({u[j],alternativeRhs[j]});
            actualAdjoint.add(a);actualAlternative.add(b);adjointScale.add(a.absolute());alternativeScale.add(b.absolute());}
        auto da=actualAdjoint;da.add(chargeWork,-1);da.add(materialWork,-1);auto db=actualAlternative;db.add(volumeWork,-1);db.add(materialWork,-1);
        check(da.absolute().compare(adjointScale.multiplied(originalWorkTolerance))<=0,"research Q-adjoint literal original 1e-12 work gate failed");
        check(db.absolute().compare(alternativeScale.multiplied(originalWorkTolerance))<=0,"research RT0-volumeQ literal original 1e-12 work gate failed");
        Work result;result.charge=chargeWork.rounded();result.volume=volumeWork.rounded();result.material=materialWork.rounded();
        result.adjointNative=actualAdjoint.rounded();result.alternativeNative=actualAlternative.rounded();result.adjointScale=adjointScale.rounded();result.alternativeScale=alternativeScale.rounded();
        result.adjointRelative=result.adjointScale==0?0:da.absolute().rounded()/result.adjointScale;
        result.alternativeRelative=result.alternativeScale==0?0:db.absolute().rounded()/result.alternativeScale;return result;
    }
};
std::array<long double,2> physicalBoundaryTensor(const QuadtreeMesh& mesh,const HydroMaxwellLoad& load) {
    std::array<FV::Dyadic,2> force;const auto& g=load.geometry();const auto& layout=g.layout();
    for(const auto& f:g.stressFacePoints()){if((f.negativeCellDof>=0&&f.positiveCellDof>=0)||f.physicalWeight==0)continue;
        const int c=layout.nativeCell(f.negativeCellDof>=0?f.negativeCellDof:f.positiveCellDof);const double sign=f.negativeCellDof>=0?1.:-1.;
        const auto t=load.stress().sample(c,f.point).tensor;const bool vertical=layout.normalAxis(f.velocityDof)==HydroDofLayout::NormalAxis::First;
        force[0].add(FV::Dyadic::product({sign,f.physicalWeight,vertical?t.xx:t.xy}));force[1].add(FV::Dyadic::product({sign,f.physicalWeight,vertical?t.xy:t.yy}));}
    return {force[0].rounded(),force[1].rounded()};
}
void slab(int n,int signedSigma) {
    constexpr double epsilonLiquid=4,epsilonGas=1,da=2,et=.5;const double sigma=signedSigma,db=da+sigma,h=1./n;
    auto dictionary=input(n,0);dictionary.Ny=8;QuadtreeMesh mesh(dictionary);FaceSegmentGraph graph;graph.build(mesh);
    const auto domain=NativeLeafDomain::build(mesh,graph,false);const auto layout=HydroDofLayout::build(mesh,domain);
    const auto divergence=HydroIntegralDivergence::build(mesh,layout);const auto metrics=HydroPhysicalDivergence::build(mesh,divergence);
    const auto pressure=HydroPressureLoad::build(divergence);using Stock=std::array<long double,2>;std::vector<Stock> stocks(mesh.cellIndexNum());
    std::vector<double> q(mesh.cellIndexNum());std::size_t charged=0;
    for(int c=0;c<mesh.cellIndexNum();++c){const auto center=mesh.cells(mesh.cellIndexToLabel[c]).centerCoord();const auto volume=domain.cellVolumes()[domain.cellRowByIndex()[c]];
        stocks[c]=center.x<.5?Stock{volume,0}:Stock{0,volume};if(center.x>.5-h&&center.x<.5){q[c]=sigma*volume/h;if(sigma!=0)++charged;}}
    const Clock initialClock{0,0,0,Clock::Position::Initial},clock{1,0,.01,Clock::Position::Midpoint};
    const auto initial=I::GeometricVofPublication::prepareLocalPlicInventory(mesh,domain,initialClock,stocks,{});
    const auto pub=I::GeometricVofPublication::prepareLocalPlicInventory(mesh,domain,clock,stocks,{});const auto phase=HydroPhaseGeometry::build(mesh,metrics,pub);
    const auto charge=E::ElectricVofCharge::fromIntegrated(mesh,initial,q);
    // Exact discrete one-dimensional chain: explicit top/bottom cell-center
    // potential traces include the known O(h) concentrated-cell offset. This
    // is not smooth MMS and cannot qualify a general physical-boundary port.
    const auto boundaryPhi=[&](Vector2D p){const double xPart=p.x<.5?-da*p.x/epsilonLiquid:
        -da*.5/epsilonLiquid-sigma*h/(2*epsilonLiquid)-db*(p.x-.5)/epsilonGas;return xPart-et*p.y;};
    std::vector<E::ElectricVofBoundary> boundaries(domain.nativeFaceCount());
    for(int row=0;row<domain.nativeFaceCount();++row){const auto f=domain.face(mesh,row);boundaries[row]=f.negativeCell>=0&&f.positiveCell>=0?
        E::ElectricVofBoundary{E::ElectricVofBoundary::Type::Interior,0}:E::ElectricVofBoundary{E::ElectricVofBoundary::Type::Dirichlet,boundaryPhi(f.center)};}
    const auto field=E::ElectricVofField::prepareOhmic(mesh,pub,charge,clock,{epsilonLiquid,epsilonGas,0,0},boundaries);
    check(Gas::sameClock(field.geometryTime(),clock)&&Gas::sameClock(field.candidateTime(),clock)&&field.sourceCharge().sameInstance(charge)&&field.chargeCandidate().integratedCharge()==q,
        "research slab changed actual source/candidate Q or full Mid clock");
    check(field.report().relativeLinearResidual<=1e-11&&field.report().maximumRelativeGaussResidual<=1e-10,
        "research slab changed retained electric solver/Gauss gates");
    const auto stress=E::ElectricVofMaxwellStress::prepare(mesh,field);const auto load=HydroMaxwellLoad::prepare(mesh,pressure,phase,stress);const ResearchLaw law(mesh,load,boundaries);
    std::vector<double> normal(layout.velocityCount()),tangent(normal.size());
    for(int u=0;u<layout.velocityCount();++u)(layout.normalAxis(u)==HydroDofLayout::NormalAxis::First?normal[u]:tangent[u])=1;
    const auto wn=law.work(normal),wt=law.work(tangent);const auto boundary=physicalBoundaryTensor(mesh,load);
    // Independent field average uses physical quadrature and no native velocity map.
    std::array<FV::Dyadic,2> columnField;FV::Dyadic columnVolume;
    for(const auto& point:load.geometry().points()){const int c=layout.nativeCell(point.cellDof);const auto center=mesh.cells(mesh.cellIndexToLabel[c]).centerCoord();
        if(!(center.x>.5-h&&center.x<.5))continue;const auto sample=stress.sample(c,point.point);columnVolume.add(point.physicalWeight);
        columnField[0].add(FV::Dyadic::product({point.physicalWeight,sample.electricField.x}));columnField[1].add(FV::Dyadic::product({point.physicalWeight,sample.electricField.y}));}
    const long double averageNormal=columnField[0].rounded()/columnVolume.rounded(),averageTangent=columnField[1].rounded()/columnVolume.rounded();
    double maximumActualDerror=0;bool actualD=true;
    for(int row=0;row<domain.nativeFaceCount();++row){const auto f=domain.face(mesh,row);if(f.horizontal||(f.center.x!=.5-h&&f.center.x!=.5))continue;
        const long double expected=f.center.x==.5?db:da,actual=static_cast<long double>(field.faces()[row].displacement)/f.area;
        maximumActualDerror=std::max(maximumActualDerror,static_cast<double>(std::abs(actual-expected)));
        actualD=reference(actual,expected,std::abs(actual)+std::abs(expected),peakReference)&&actualD;}
    const long double bodyNormal=sigma*(da+db)/(2*epsilonLiquid),bodyTangent=sigma*et;
    const long double materialNormal=.5L*db*db*(1.L/epsilonGas-1.L/epsilonLiquid)-.5L*et*et*(epsilonGas-epsilonLiquid);
    const long double totalNormal=.5L*(db*db/epsilonGas-da*da/epsilonLiquid)-.5L*et*et*(epsilonGas-epsilonLiquid),totalTangent=sigma*et;
    const long double twoPointNormal=.5L*sigma*da/epsilonLiquid+.25L*sigma*db*(1.L/epsilonLiquid+1.L/epsilonGas);
    const long double predictedError=.25L*sigma*db*(1.L/epsilonGas-1.L/epsilonLiquid);
    const long double scale=std::abs(bodyNormal)+std::abs(bodyTangent)+std::abs(materialNormal)+std::abs(totalNormal)+std::abs(totalTangent)+std::abs(twoPointNormal)+std::abs(predictedError);
    bool exactReferences=true;const auto observe=[&](long double actual,long double wanted){const bool ok=reference(actual,wanted,scale,peakReference);exactReferences=ok&&exactReferences;};
    observe(wn.volume,bodyNormal);observe(wt.volume,bodyTangent);observe(wn.material,materialNormal);observe(wt.material,0);
    observe(wn.alternativeNative,totalNormal);observe(wt.alternativeNative,totalTangent);observe(boundary[0],totalNormal);observe(boundary[1],totalTangent);
    exactReferences=reference(averageNormal,(da+db)/(2*epsilonLiquid),std::abs(averageNormal)+(da+db)/(2*epsilonLiquid),peakReference)&&exactReferences;
    exactReferences=reference(averageTangent,et,std::abs(averageTangent)+std::abs(et),peakReference)&&actualD&&exactReferences;
    const bool identity=reference(wn.charge,twoPointNormal,scale,peakCounterexampleIdentity)&&reference(wn.charge-bodyNormal,predictedError,scale,peakCounterexampleIdentity)
        &&reference(wt.charge,bodyTangent,scale,peakCounterexampleIdentity);
    const bool rejected=sigma!=0&&std::abs(wn.charge-bodyNormal)>originalWorkTolerance*scale;
    physicalReferences=physicalReferences&&exactReferences;counterexampleIdentities=counterexampleIdentities&&identity&&(sigma==0||rejected);
    if(rejected)++counterexamples;++slabCases;
    std::cout<<"HELMHOLTZ_RESEARCH_1D,n,"<<n<<",h,"<<h<<",sigma,"<<sigma<<",charged_cells,"<<charged<<",Da,"<<da<<",Db,"<<db<<",Et,"<<et
        <<",actual_Mid_gid,"<<pub.identity().value()<<",actual_field_Gauss,"<<field.report().maximumRelativeGaussResidual
        <<",Q_adjoint_normal,"<<wn.charge<<",correct_bulk_normal_force,"<<bodyNormal<<",Q_adjoint_tangent,"<<wt.charge<<",correct_bulk_tangent,"<<bodyTangent
        <<",RT0_charged_cell_En_average,"<<averageNormal<<",independent_En_average,"<<(da+db)/(2*epsilonLiquid)<<",RT0_charged_cell_Et_average,"<<averageTangent
        <<",independent_Et_average,"<<et<<",maximum_actual_D_error,"<<maximumActualDerror
        <<",counterexample_error,"<<wn.charge-bodyNormal<<",h_independent_error_formula,"<<predictedError<<",two_point_physical_assertion,"<<(sigma==0?"Q0_CONTROL":(rejected?"FAIL":"UNEXPECTED_NO_FAILURE"))
        <<",counterexample_identity,"<<(identity?"PASS":"FAIL")<<",RT0_volumeQ_normal,"<<wn.volume<<",RT0_volumeQ_tangent,"<<wt.volume
        <<",secant_normal,"<<wn.material<<",secant_tangent,"<<wt.material<<",alternative_total_normal,"<<wn.alternativeNative<<",alternative_total_tangent,"<<wt.alternativeNative
        <<",independent_Maxwell_farfield_normal,"<<totalNormal<<",independent_Maxwell_farfield_tangent,"<<totalTangent
        <<",same_field_full_boundary_tensor_normal,"<<boundary[0]<<",same_field_full_boundary_tensor_tangent,"<<boundary[1]<<",alternative_1D_reference,"<<(exactReferences?"PASS":"FAIL")
        <<",original_literal_work_tolerance,"<<originalWorkTolerance<<",Q_literal_relative,"<<std::max(wn.adjointRelative,wt.adjointRelative)
        <<",volumeQ_literal_relative,"<<std::max(wn.alternativeRelative,wt.alternativeRelative)<<",explicit_fixture_boundary_half_cells,1,discrete_chain_not_smooth_MMS,1\n";
}
void gasProbe(const QuadtreeMesh& mesh,const HydroMaxwellLoad& load,const ResearchLaw& law,Vector2D site,bool rz,int number) {
    // The original six certificates/gap run on this SAME const Field/load.
    Gas::probe(mesh,load,site,rz,number);const auto compact=Gas::cycle(mesh,load.geometry().layout(),site,rz);const auto work=law.work(compact.velocity);
    const auto& layout=load.geometry().layout();const auto& g=load.geometry();const auto& pub=load.phaseGeometry().vofPublication();
    std::vector<unsigned char> affected(mesh.cellIndexNum());
    const auto touches=[&](const Matrix& matrix,Eigen::Index row){for(Matrix::InnerIterator e(matrix,row);e;++e)if(e.value()!=0&&compact.velocity[e.col()]!=0)return true;return false;};
    for(std::size_t q=0;q<g.points().size();++q){bool active=false;for(int k=0;k<g.components();++k)active=active||touches(g.matrix(),static_cast<Eigen::Index>(q)*g.components()+k);
        for(int k=0;k<2;++k)active=active||touches(g.stressVolumeVelocityMatrix(),static_cast<Eigen::Index>(q)*2+k);if(active)affected[layout.nativeCell(g.points()[q].cellDof)]=1;}
    for(std::size_t q=0;q<g.stressFacePoints().size();++q){if(!touches(g.stressFaceVelocityMatrix(),static_cast<Eigen::Index>(q)*2)&&!touches(g.stressFaceVelocityMatrix(),static_cast<Eigen::Index>(q)*2+1))continue;
        const auto& f=g.stressFacePoints()[q];if(f.negativeCellDof>=0)affected[layout.nativeCell(f.negativeCellDof)]=1;if(f.positiveCellDof>=0)affected[layout.nativeCell(f.positiveCellDof)]=1;}
    double gap=std::numeric_limits<double>::infinity();std::size_t support=0;
    for(int c=0;c<mesh.cellIndexNum();++c)if(affected[c]){++support;const auto& cell=pub.plicCells()[c];check(cell.inventory[0]==0&&cell.inventory[1]>0&&law.field.epsilonByCell()[c]==1&&law.field.chargeCandidate().integratedCharge()[c]==0,
        "research fixed probe lost original exact zero-source support");const auto& b=cell.rectangle;const double cx=rz?0.:.47;
        const double dx=b.lower.x-cx,dy=b.lower.y-Gas::centerZ,dz=b.upper.y-Gas::centerZ,nearZ=dy>0?dy:(dz<0?-dz:0.);
        const double nearX=dx>0?dx:(b.upper.x<cx?cx-b.upper.x:0.);gap=std::min(gap,std::hypot(nearX,nearZ)-Gas::inclusionRadius);}
    check(support>0&&gap>0&&work.charge==0&&work.volume==0&&work.material==0&&work.adjointNative==0&&work.alternativeNative==0,
        "research zero-source law is not exact on fixed same-field gas probe");++gasProbes;
    std::cout<<"HELMHOLTZ_RESEARCH_FIXED_GAS,"<<(rz?"RZ":"XY")<<",probe,"<<number<<",vertex_x,"<<site.x<<",vertex_y,"<<site.y
        <<",actual_Mid_gid,"<<pub.identity().value()<<",support_cells,"<<support<<",original_analytic_interface_gap,"<<gap
        <<",Q_adjoint_work,"<<work.charge<<",RT0_volumeQ_work,"<<work.volume<<",secant_material_work,"<<work.material
        <<",Q_adjoint_native_work,"<<work.adjointNative<<",alternative_native_work,"<<work.alternativeNative
        <<",same_actual_Field_geometry_fullclock,1,zero_source_algebra_not_gas_mask,1,physical_boundary_compact_support,1,production_2D_RZ_qualification,0\n";
}
void gasScene(bool rz) {
    QuadtreeMesh mesh(input(32,rz?1:0));if(rz){DictionaryMesh::RefineBox box;box.min={.5,0};box.max={1,1};box.level=1;mesh.applyRefineBoxes({box});}
    FaceSegmentGraph graph;graph.build(mesh);const auto domain=NativeLeafDomain::build(mesh,graph,rz);
    check(domain.uniformLeaves()!=rz,"research changed fixed XY/RZ2:1 gas mesh");
    const auto initial=I::GeometricVofInitialization::prepare(mesh,domain,[rz](Vector2D p){return std::hypot(p.x-(rz?0.:.47),p.y-Gas::centerZ)-Gas::inclusionRadius;});
    using Stock=std::array<long double,2>;std::vector<Stock> stocks;for(const auto& p:initial.plicCells())stocks.push_back(p.inventory);
    const Clock clock{1,0,.01,Clock::Position::Midpoint};const auto pub=I::GeometricVofPublication::prepareLocalPlicInventory(mesh,domain,clock,stocks,{});
    const auto charge=E::ElectricVofCharge::fromIntegrated(mesh,initial,std::vector<double>(mesh.cellIndexNum()));
    const auto layout=HydroDofLayout::build(mesh,domain);const auto divergence=HydroIntegralDivergence::build(mesh,layout);
    const auto metrics=HydroPhysicalDivergence::build(mesh,divergence);const auto pressure=HydroPressureLoad::build(divergence);const auto phase=HydroPhaseGeometry::build(mesh,metrics,pub);
    std::vector<E::ElectricVofBoundary> boundaries(domain.nativeFaceCount());
    for(int row=0;row<domain.nativeFaceCount();++row){const auto f=domain.face(mesh,row);boundaries[row]=f.negativeCell>=0&&f.positiveCell>=0?
        E::ElectricVofBoundary{E::ElectricVofBoundary::Type::Interior,0}:(f.area==0?E::ElectricVofBoundary{E::ElectricVofBoundary::Type::HomogeneousNeumann,0}:
        E::ElectricVofBoundary{E::ElectricVofBoundary::Type::Dirichlet,Gas::analyticPotential(f.center,rz)});}
    const auto field=E::ElectricVofField::prepareOhmic(mesh,pub,charge,clock,{Gas::contrast,1,0,0},boundaries);field.validateFor(mesh,pub,charge,clock);
    check(field.report().relativeLinearResidual<=1e-11&&field.report().maximumRelativeGaussResidual<=1e-10,
        "research fixed gas fixture changed retained Electric solver/Gauss gates");
    std::cout<<"HELMHOLTZ_RESEARCH_FIXED_SCENE,"<<(rz?"RZ":"XY")<<",N,32,epsilon_ratio,2.5,K,0,Q,0,actual_Mid_gid,"<<pub.identity().value()
        <<",actual_Mid_time,"<<clock.time()<<",field_SuperLU,"<<field.report().usedSuperLU<<",field_GMRES,"<<field.report().usedGMRES
        <<",linear_residual,"<<field.report().relativeLinearResidual<<",Gauss,"<<field.report().maximumRelativeGaussResidual<<",same_original_field_fixture,1\n";
    check(Gas::sameClock(field.geometryTime(),clock)&&Gas::sameClock(field.candidateTime(),clock)&&field.sourceCharge().sameInstance(charge)&&field.chargeCandidate().integratedCharge()==charge.integratedCharge(),
        "research fixed gas field changed source Q/full clock");
    const auto stress=E::ElectricVofMaxwellStress::prepare(mesh,field);const auto load=HydroMaxwellLoad::prepare(mesh,pressure,phase,stress);const ResearchLaw law(mesh,load,boundaries);
    if(rz){gasProbe(mesh,load,law,{.65625,.75},true,0);gasProbe(mesh,load,law,{.71875,.6875},true,1);gasProbe(mesh,load,law,{.5,.703125},true,2);gasProbe(mesh,load,law,{.03125,.875},true,3);}
    else {gasProbe(mesh,load,law,{.8125,.78125},false,0);gasProbe(mesh,load,law,{.78125,.8125},false,1);}
}
void run() {
    for(int n:{8,16,32})for(int sigma:{0,-1,1})slab(n,sigma);gasScene(false);gasScene(true);
    std::cout<<"HELMHOLTZ_SOURCE_RESEARCH_COMPLETE,slab_cases,"<<slabCases<<",expected_counterexamples,"<<counterexamples<<",fixed_gas_probes,"<<gasProbes
        <<",two_point_physical_status,FAIL,rejected_two_point,1,counterexample_identities,"<<(counterexampleIdentities?"PASS":"FAIL")
        <<",alternative_1D_reference,"<<(physicalReferences?"PASS":"FAIL")<<",peak_1D_reference_relative,"<<peakReference<<",peak_counterexample_identity_relative,"<<peakCounterexampleIdentity
        <<",original_virtualwork_tolerance,"<<originalWorkTolerance<<",virtual_Q_shared_face_conservation_exact,1,alternative_1D_evidence_only,1,MMS_not_run,1,production_2D_RZ_qualified,0,no_NS_no_App_no_inventory_update\n";
    check(slabCases==9&&counterexamples==6&&gasProbes==6&&counterexampleIdentities,"research did not preserve the rejected two-point concentrated-Q counterexample");
    check(physicalReferences,"research RT0-volumeQ/secant 1D independent physical reference FAIL; do not change weights/tolerance");
}
} // namespace HelmholtzSourceResearch
// END HELMHOLTZ_SOURCE_RESEARCH_PHASE1

// BEGIN HELMHOLTZ_SMOOTH_SIGNEDQ_N32_DIAGNOSTIC: test only; no accuracy PASS gate.
namespace HelmholtzSmoothSignedQDiagnostic {
namespace R=HelmholtzSourceResearch;namespace E=LAMP::Electric;namespace I=LAMP::Interface;
namespace FV=LAMP::Geometry::Detail::FiniteVolume;
using Clock=LAMP::TimeLayer;using Matrix=HydroStrainOperator::Matrix;
constexpr int baseN=32;
struct Exact {
    long double epsilon=0,potential=0,charge=0,psi=0;
    std::array<long double,2> electric{},velocity{};
    // (d_x v_x,d_z v_x,d_x v_z,d_z v_z,v_r/r); last component XY=0.
    std::array<long double,5> gradient{};
};
std::array<long double,3> bump(long double x,long double a) {
    if(std::abs(x)>=a)return {};
    const long double aa=a*a,t=1-x*x/aa;
    return {t*t*t,-6*x*t*t/aa,-6*t*t/aa+24*x*x*t/(aa*aa)};
}
Exact exact(Vector2D point,bool rz) {
    Exact f;const long double x=point.x,z=static_cast<long double>(point.y)-.5L;
    const auto bz=bump(z,.25L);const long double c=z*bz[0],dc=bz[0]+z*bz[1],ddc=2*bz[1]+z*bz[2];
    if(rz){const auto b=bump(x,.5L);f.epsilon=1+4*x*x;f.potential=x*x*z;f.charge=-4*z*(1+8*x*x);
        f.electric={-2*x*z,-x*x};f.psi=x*x*b[0]*c;
        f.velocity={x*b[0]*dc,-(2*b[0]+x*b[1])*c};
        f.gradient={(b[0]+x*b[1])*dc,x*b[0]*ddc,-(3*b[1]+x*b[2])*c,-(2*b[0]+x*b[1])*dc,b[0]*dc};
    }else {const long double xi=x-.5L;const auto b=bump(xi,.25L);
        const long double a=(xi+.25L)*b[0],da=b[0]+(xi+.25L)*b[1],dda=2*b[1]+(xi+.25L)*b[2];
        f.epsilon=1+xi;f.potential=xi*xi*z;f.charge=-2*z*(1+2*xi);f.electric={-2*xi*z,-xi*xi};
        f.psi=a*c;f.velocity={a*dc,-da*c};f.gradient={da*dc,a*ddc,-dda*c,-da*dc,0};}
    return f;
}
struct Probe {
    std::vector<double> velocity;
    std::vector<long double> gross;
    std::size_t nonzero=0,hanging=0,axisFaces=0;
    double measureDefect=0,divergenceDefect=0;
};
Probe probe(const QuadtreeMesh& mesh,const HydroDofLayout& layout,bool rz) {
    Probe result;result.velocity.resize(layout.velocityCount());result.gross.resize(layout.pressureCount());
    std::vector<FV::Dyadic> endpointBalance(layout.pressureCount()),representedBalance(layout.pressureCount());
    std::map<int,long double> vertexPsi;
    const auto psi=[&](int v){const auto found=vertexPsi.find(v);if(found!=vertexPsi.end())return found->second;
        const auto value=exact(mesh.verts(v).coord(),rz).psi;vertexPsi.emplace(v,value);return value;};
    for(int u=0;u<layout.velocityCount();++u){const auto f=layout.domain().face(mesh,layout.nativeFace(u));
        const auto& edge=f.horizontal?mesh.hEdges(f.edgeLabel):mesh.vEdges(f.edgeLabel);int va=edge.edgeToVert[0],vb=edge.edgeToVert[1];
        auto a=mesh.verts(va).coord(),b=mesh.verts(vb).coord();if(f.horizontal?b.x<a.x:b.y<a.y){std::swap(va,vb);std::swap(a,b);}
        const long double length=f.horizontal?static_cast<long double>(b.x)-a.x:static_cast<long double>(b.y)-a.y;
        const long double metric=rz?static_cast<long double>(CartesianMeasure::twoPi)*(f.horizontal?(static_cast<long double>(a.x)+b.x)/2:a.x):1.L;
        const long double area=length*metric;
        MaxwellLoadContract::bounded(area-f.area,area,result.measureDefect,"smooth probe independent native face physical measure");
        FV::Dyadic delta(psi(vb));delta.add(psi(va),-1);auto flux=delta.multiplied(f.horizontal?-1.L:1.L);
        if(rz)flux=flux.multiplied(CartesianMeasure::twoPi);
        if(area==0){check(rz&&a.x==0&&flux.compare(FV::Dyadic{})==0,"smooth probe nonzero flux at the true zero-area axis");++result.axisFaces;continue;}
        // Independent physical face mean, with common 2*pi cancelled in RZ.
        const long double denominator=rz?(f.horizontal?(static_cast<long double>(b.x)*b.x-static_cast<long double>(a.x)*a.x)/2:a.x*length):length;
        result.velocity[u]=static_cast<double>((f.horizontal?-1.L:1.L)*delta.rounded()/denominator);
        if(flux.compare(FV::Dyadic{})!=0){check(f.negativeCell>=0&&f.positiveCell>=0,"smooth analytic probe reaches a physical exterior face");++result.nonzero;
            const int lm=mesh.cells(mesh.cellIndexToLabel[f.negativeCell]).level(),lp=mesh.cells(mesh.cellIndexToLabel[f.positiveCell]).level();
            check(std::abs(lm-lp)<=1,"smooth diagnostic native faces violate actual 2:1");if(lm!=lp)++result.hanging;}
        const auto physical=FV::Dyadic::product({f.area,result.velocity[u]});
        for(int side=0;side<2;++side){const int c=side==0?f.negativeCell:f.positiveCell;if(c<0)continue;
            const int row=layout.pressureDof(c),sign=side==0?1:-1;endpointBalance[row].add(flux,sign);representedBalance[row].add(physical,sign);
            result.gross[row]+=physical.absolute().rounded();}
    }
    check(result.nonzero>0&&result.hanging>0&&!layout.domain().uniformLeaves(),"smooth N32 family misses nonzero native hanging-face support");
    for(int p=0;p<layout.pressureCount();++p){check(endpointBalance[p].compare(FV::Dyadic{})==0,"shared smooth endpoint flux ledger fails exact telescoping");
        MaxwellLoadContract::bounded(representedBalance[p].rounded(),result.gross[p],result.divergenceDefect,"represented smooth face-mean physical divergence");}
    check(!rz||result.axisFaces>0,"smooth RZ family misses the true axis");return result;
}
void scene(bool rz,int n=baseN,bool convergence=false) {
    check(n==32||n==64||n==128,"smooth diagnostic requires the frozen N32/N64/N128 resolutions");
    auto dictionary=input(n,rz?1:0);dictionary.maxCellNum=std::max(dictionary.maxCellNum,8*n*n); // Storage only; N32/N64 unchanged.
    QuadtreeMesh mesh(dictionary);DictionaryMesh::RefineBox box;box.min={rz?.25:.625,0};box.max={1,1};box.level=1;mesh.applyRefineBoxes({box});
    FaceSegmentGraph graph;graph.build(mesh);const auto domain=NativeLeafDomain::build(mesh,graph,rz);const auto layout=HydroDofLayout::build(mesh,domain);
    const int expectedLeaves=rz?13*n*n/4:17*n*n/8;
    check(mesh.cellIndexNum()==expectedLeaves&&layout.pressureCount()==expectedLeaves,"smooth refinement was truncated or changed the frozen actual leaf family");
    const auto divergence=HydroIntegralDivergence::build(mesh,layout);const auto metrics=HydroPhysicalDivergence::build(mesh,divergence);const auto pressure=HydroPressureLoad::build(divergence);
    using Stock=std::array<long double,2>;std::vector<Stock> stocks(mesh.cellIndexNum());std::vector<double> q(mesh.cellIndexNum());
    std::size_t positiveQ=0,negativeQ=0;long double chargeStorageError=0;
    for(int c:domain.cellIndexByRow()){const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);const auto center=cell.centerCoord(),size=cell.size();
        const long double r0=static_cast<long double>(center.x)-size.x/2,r1=static_cast<long double>(center.x)+size.x/2;
        const long double mean=rz?(r0*r0+r1*r1)/2:static_cast<long double>(center.x),volume=domain.cellVolumes()[domain.cellRowByIndex()[c]];
        stocks[c]={mean*volume,(1-mean)*volume};const long double eta=static_cast<long double>(center.y)-.5L;
        const long double density=rz?-4*eta*(1+8*mean):-2*eta*(1+2*(mean-.5L));q[c]=static_cast<double>(density*volume);
        chargeStorageError=std::max(chargeStorageError,std::abs(static_cast<long double>(q[c])/volume-density));if(q[c]>0)++positiveQ;if(q[c]<0)++negativeQ;}
    check(positiveQ>0&&negativeQ>0,"smooth fixture lost nonuniform signed bulk Q");
    const Clock initialClock{0,0,0,Clock::Position::Initial},clock{1,0,.01,Clock::Position::Midpoint};
    const auto initial=I::GeometricVofPublication::prepareLocalPlicInventory(mesh,domain,initialClock,stocks,{});
    const auto pub=I::GeometricVofPublication::prepareLocalPlicInventory(mesh,domain,clock,stocks,{});const auto phase=HydroPhaseGeometry::build(mesh,metrics,pub);
    const auto charge=E::ElectricVofCharge::fromIntegrated(mesh,initial,q);std::vector<E::ElectricVofBoundary> boundaries(domain.nativeFaceCount());
    for(int row=0;row<domain.nativeFaceCount();++row){const auto f=domain.face(mesh,row);
        boundaries[row]=f.negativeCell>=0&&f.positiveCell>=0?E::ElectricVofBoundary{E::ElectricVofBoundary::Type::Interior,0}:
            (f.area==0?E::ElectricVofBoundary{E::ElectricVofBoundary::Type::HomogeneousNeumann,0}:
                E::ElectricVofBoundary{E::ElectricVofBoundary::Type::Dirichlet,static_cast<double>(exact(f.center,rz).potential)});}
    // Actual public Dirichlet data are continuous Phi at actual face centres.
    // This is not the h-dependent discrete-chain boundary used in phase 1.
    const E::ElectricVofMaterials material=rz?E::ElectricVofMaterials{5,1,0,0}:E::ElectricVofMaterials{1.5,.5,0,0};
    const auto field=E::ElectricVofField::prepareOhmic(mesh,pub,charge,clock,material,boundaries);field.validateFor(mesh,pub,charge,clock);
    check(R::Gas::sameClock(field.geometryTime(),clock)&&R::Gas::sameClock(field.candidateTime(),clock)&&field.geometry().identity()==pub.identity()&&
        field.sourceCharge().sameInstance(charge)&&field.chargeCandidate().integratedCharge()==q,"smooth actual Field/source/candidate Q/full clock owner mismatch");
    check(field.report().relativeLinearResidual<=1e-11&&field.report().maximumRelativeGaussResidual<=1e-10,"smooth fixture changed retained Electric linear/Gauss gates");
    const auto stress=E::ElectricVofMaxwellStress::prepare(mesh,field);const auto load=HydroMaxwellLoad::prepare(mesh,pressure,phase,stress);const R::ResearchLaw law(mesh,load,boundaries);
    auto compact=probe(mesh,layout,rz);const HydroNativeValues values(layout,compact.velocity,std::vector<double>(layout.pressureCount()));
    const auto actualDI=divergence.apply(mesh,values);
    for(int p=0;p<layout.pressureCount();++p)MaxwellLoadContract::bounded(actualDI[p],compact.gross[p],compact.divergenceDefect,"smooth probe authoritative DI of independent face means");
    const auto& geometry=load.geometry();check(geometry.hasCanonicalRawStrain(),"smooth diagnostic lost public canonical raw CSR");
    const auto& rawV=geometry.stressVolumeVelocityMatrix();const auto rawS=geometry.canonicalRawStrain();
    const auto touches=[&](const Matrix& matrix,Eigen::Index row){for(Matrix::InnerIterator e(matrix,row);e;++e)if(e.value()!=0&&compact.velocity[e.col()]!=0)return true;return false;};
    std::vector<unsigned char> affected(mesh.cellIndexNum());
    for(std::size_t i=0;i<geometry.points().size();++i){bool active=touches(rawV,i*2)||touches(rawV,i*2+1);
        for(int k=0;k<geometry.components();++k)active=active||touches(rawS.matrix(),i*geometry.components()+k);
        if(active)affected[layout.nativeCell(geometry.points()[i].cellDof)]=1;}
    for(std::size_t i=0;i<geometry.stressFacePoints().size();++i){if(!touches(geometry.stressFaceVelocityMatrix(),i*2)&&!touches(geometry.stressFaceVelocityMatrix(),i*2+1))continue;
        const auto& f=geometry.stressFacePoints()[i];check((f.negativeCellDof>=0&&f.positiveCellDof>=0)||(rz&&f.point.x==0&&f.physicalWeight==0),
            "smooth reconstructed/shared-face compact closure reaches positive-area physical boundary");
        for(int p:{f.negativeCellDof,f.positiveCellDof})if(p>=0)affected[layout.nativeCell(p)]=1;}
    std::size_t support=0,axisSupport=0;long double wallGap=1;
    for(int c:domain.cellIndexByRow())if(affected[c]){++support;const auto& b=pub.plicCells()[c].rectangle;
        long double gap=std::min({1.L-b.upper.x,static_cast<long double>(b.lower.y),1.L-b.upper.y});if(!rz)gap=std::min(gap,static_cast<long double>(b.lower.x));
        wallGap=std::min(wallGap,gap);check(gap>0,"smooth raw-volume/raw-strain compact cell closure touches a positive-area exterior wall");
        if(rz&&b.lower.x==0){++axisSupport;const auto sample=stress.sample(c,{0,(b.lower.y+b.upper.y)/2});
            check(sample.electricField.x==0&&sample.displacement.x==0&&sample.tensorDivergence.x==0&&std::isfinite(sample.electricField.y)&&std::isfinite(sample.tensorDivergence.y),
                "smooth actual RT0 axis/hoop limit is not finite and radially regular");}}
    check(support>0&&(!rz||axisSupport>0),"smooth compact source/axis closure has no actual support witness");
    const auto work=law.work(compact.velocity);const Eigen::Map<const Eigen::VectorXd> native(compact.velocity.data(),compact.velocity.size());
    const Eigen::VectorXd reconstructed=rawV*native;const auto rawStrain=rawS.apply(mesh,values);
    FV::Dyadic oracleQbarRV,oracleQexactRV,oracleQexactV,oracleMaterialV;
    FV::Dyadic analyticMeridional,analyticHoop,rawMeridional,rawHoop,analyticPressureMeridional,analyticPressureHoop,rawPressureMeridional,rawPressureHoop;
    long double fieldError2=0,fieldNorm2=0,velocityError2=0,velocityNorm2=0,maximumFieldError=0,epsilonMeanError=0;
    for(std::size_t i=0;i<geometry.points().size();++i){const auto& point=geometry.points()[i];const int c=layout.nativeCell(point.cellDof);const auto f=exact(point.point,rz);
        const auto sample=stress.sample(c,point.point);const long double w=point.physicalWeight,volume=domain.cellVolumes()[domain.cellRowByIndex()[c]],qbar=field.chargeCandidate().integratedCharge()[c]/volume;
        const long double ex=f.electric[0],ez=f.electric[1],e2=ex*ex+ez*ez,pr=.5L*f.epsilon*e2;
        const std::array<long double,4> tensor{.5L*f.epsilon*(ex*ex-ez*ez),.5L*f.epsilon*(ez*ez-ex*ex),f.epsilon*ex*ez,-pr};
        for(int k=0;k<2;++k){const long double v=reconstructed[static_cast<Eigen::Index>(i)*2+k],e=k==0?sample.electricField.x:sample.electricField.y;
            const long double de=e-f.electric[k],dv=v-f.velocity[k];fieldError2+=w*de*de;fieldNorm2+=w*f.electric[k]*f.electric[k];
            velocityError2+=w*dv*dv;velocityNorm2+=w*f.velocity[k]*f.velocity[k];maximumFieldError=std::max(maximumFieldError,std::abs(de));
            oracleQbarRV.add(FV::Dyadic::product({w,qbar,f.electric[k],v}));oracleQexactRV.add(FV::Dyadic::product({w,f.charge,f.electric[k],v}));
            oracleQexactV.add(FV::Dyadic::product({w,f.charge,f.electric[k],f.velocity[k]}));}
        oracleMaterialV.add(FV::Dyadic::product({w,-.5L*e2,rz?8.L*point.point.x:1.L,f.velocity[0]}));
        analyticMeridional.add(FV::Dyadic::product({-w,tensor[0],f.gradient[0]}));analyticMeridional.add(FV::Dyadic::product({-w,tensor[1],f.gradient[3]}));
        analyticMeridional.add(FV::Dyadic::product({-w,tensor[2],f.gradient[1]+f.gradient[2]}));
        analyticPressureMeridional.add(FV::Dyadic::product({-w,pr,f.gradient[0]+f.gradient[3]}));
        const auto row=static_cast<Eigen::Index>(i)*geometry.components();const int shear=rz?3:2;
        rawMeridional.add(FV::Dyadic::product({-w,tensor[0],rawStrain[row]}));rawMeridional.add(FV::Dyadic::product({-w,tensor[1],rawStrain[row+1]}));
        rawMeridional.add(FV::Dyadic::product({-w,std::sqrt(2.L)*tensor[2],rawStrain[row+shear]}));
        rawPressureMeridional.add(FV::Dyadic::product({-w,pr,static_cast<long double>(rawStrain[row])+rawStrain[row+1]}));
        if(rz){analyticHoop.add(FV::Dyadic::product({-w,tensor[3],f.gradient[4]}));rawHoop.add(FV::Dyadic::product({-w,tensor[3],rawStrain[row+2]}));
            analyticPressureHoop.add(FV::Dyadic::product({-w,pr,f.gradient[4]}));rawPressureHoop.add(FV::Dyadic::product({-w,pr,rawStrain[row+2]}));}
    }
    for(int c:domain.cellIndexByRow()){const auto& b=pub.plicCells()[c].rectangle;const long double mean=rz?(static_cast<long double>(b.lower.x)*b.lower.x+static_cast<long double>(b.upper.x)*b.upper.x)/2:
            (static_cast<long double>(b.lower.x)+b.upper.x)/2;
        const long double wanted=rz?1+4*mean:.5L+mean;epsilonMeanError=std::max(epsilonMeanError,std::abs(field.epsilonByCell()[c]-wanted));}
    std::vector<double> p(layout.pressureCount());for(int row=0;row<layout.pressureCount();++row){const auto f=exact(mesh.cells(mesh.cellIndexToLabel[layout.nativeCell(row)]).centerCoord(),rz);
        p[row]=static_cast<double>(.5L*f.epsilon*(f.electric[0]*f.electric[0]+f.electric[1]*f.electric[1]));}
    const auto cp=pressure.apply(mesh,HydroNativeValues(layout,std::vector<double>(layout.velocityCount()),p));FV::Dyadic cpWork,cpScale,pressureAdjoint;
    for(std::size_t u=0;u<compact.velocity.size();++u){const auto term=FV::Dyadic::product({compact.velocity[u],cp[u]});cpWork.add(term);cpScale.add(term.absolute());}
    for(std::size_t row=0;row<p.size();++row)pressureAdjoint.add(FV::Dyadic::product({p[row],actualDI[row]}));
    double pressureNull=0;MaxwellLoadContract::bounded(cpWork.rounded(),cpScale.rounded(),pressureNull,"smooth native C pressure-null work");
    auto cpIdentity=cpWork;cpIdentity.add(pressureAdjoint);MaxwellLoadContract::bounded(cpIdentity.rounded(),cpScale.rounded(),pressureNull,"smooth native C/DI pressure adjoint identity");
    const long double piL=std::acos(-1.L),refQ=rz?-piL/25200:-1.L/79380,refMaterial=rz?piL/75600:1.L/396900,refTotal=rz?-piL/37800:-1.L/99225;
    const long double refMeridional=rz?0:refTotal,refHoop=rz?refTotal:0,refPressureMeridional=rz?refTotal:0,refPressureHoop=rz?-refTotal:0;
    const long double b=oracleQbarRV.rounded(),c=oracleQexactRV.rounded(),d=oracleQexactV.rounded();
    const auto metric=rz?"RZ":"XY";
    const char* prefix=convergence?"HELMHOLTZ_SMOOTH_CONVERGENCE":"HELMHOLTZ_SMOOTH_N32";
    const auto resolution=convergence?std::string(",N,")+std::to_string(n):std::string{};
    std::cout<<prefix<<"_FIXTURE,"<<metric<<",N,"<<n<<",refine_level,1,transition,"<<(rz?.25:.625)<<",actual_cells,"<<layout.pressureCount()
        <<",nonzero_native_faces,"<<compact.nonzero<<",nonzero_hanging_faces,"<<compact.hanging<<",support_cells,"<<support<<",axis_support_cells,"<<axisSupport
        <<",minimum_positive_wall_gap,"<<wallGap<<",cached_finite_endpoint_telescoping_exact,1,native_face_mean_roundoff_checked,1,native_DI_relative,"<<compact.divergenceDefect<<",native_C_pressure_null_relative,"<<pressureNull
        <<",positive_Q_cells,"<<positiveQ<<",negative_Q_cells,"<<negativeQ<<",Q_density_storage_abs_error,"<<chargeStorageError<<",epsilon_cell_mean_abs_error,"<<epsilonMeanError
        <<",actual_Mid_gid,"<<pub.identity().value()<<",source_Q_id,"<<charge.identity()<<",clock,("<<clock.step<<','<<clock.begin<<','<<clock.end<<','<<static_cast<int>(clock.position)
        <<"),same_actual_Field_chargeCandidate_RawCSR,1,continuous_face_Dirichlet_fixture,1,axis_actual_zero_area_Neumann,"<<rz<<",smooth_material_inventory_not_sharp_interface,1\n";
    std::cout<<prefix<<"_PHYSICAL_DIAGNOSTIC,"<<metric<<resolution<<",actual_RT0_E_L2_relative,"<<std::sqrt(fieldError2/fieldNorm2)<<",actual_RT0_E_max_abs_error,"<<maximumFieldError
        <<",raw_RV_velocity_L2_relative,"<<std::sqrt(velocityError2/velocityNorm2)<<",actual_volumeQ,"<<work.volume<<",analytic_WQ,"<<refQ<<",WQ_signed_error,"<<work.volume-refQ<<",WQ_relative_error,"<<std::abs((work.volume-refQ)/refQ)
        <<",actual_secant_material,"<<work.material<<",analytic_Wepsilon,"<<refMaterial<<",Wepsilon_signed_error,"<<work.material-refMaterial<<",Wepsilon_relative_error,"<<std::abs((work.material-refMaterial)/refMaterial)
        <<",actual_alternative_native,"<<work.alternativeNative<<",analytic_Wtotal,"<<refTotal<<",Wtotal_signed_error,"<<work.alternativeNative-refTotal<<",Wtotal_relative_error,"<<std::abs((work.alternativeNative-refTotal)/refTotal)
        <<",original_literal_tolerance,"<<R::originalWorkTolerance<<",literal_relative,"<<work.alternativeRelative<<",physical_status,DIAGNOSTIC_NO_ACCURACY_GATE,physical_PASS_claim,0\n";
    std::cout<<prefix<<"_REFERENCE_ONLY,"<<metric<<resolution<<",analyticE_Qbar_RV,"<<b<<",analyticE_qexact_RV,"<<c<<",analyticE_qexact_vexact_public_quadrature,"<<d
        <<",actual_field_and_native_RNE_work_component,"<<work.volume-b<<",Q_P0_work_component,"<<b-c<<",RV_work_component,"<<c-d<<",quadrature_work_component,"<<d-refQ
        <<",analytic_material_vexact_public_quadrature,"<<oracleMaterialV.rounded()<<",material_quadrature_error,"<<oracleMaterialV.rounded()-refMaterial<<",reference_only_no_Field_replacement,1\n";
    std::cout<<prefix<<"_HOOP_PRESSURE,"<<metric<<resolution<<",analyticT_vexact_public_meridional,"<<analyticMeridional.rounded()<<",meridional_reference,"<<refMeridional
        <<",analyticT_vexact_public_hoop,"<<analyticHoop.rounded()<<",hoop_reference,"<<refHoop<<",analyticT_vexact_public_full,"<<analyticMeridional.rounded()+analyticHoop.rounded()<<",full_reference,"<<refTotal
        <<",analyticT_rawS_meridional,"<<rawMeridional.rounded()<<",analyticT_rawS_hoop,"<<rawHoop.rounded()<<",analyticT_rawS_full,"<<rawMeridional.rounded()+rawHoop.rounded()
        <<",pressure_vexact_meridional,"<<analyticPressureMeridional.rounded()<<",pressure_meridional_reference,"<<refPressureMeridional<<",pressure_vexact_hoop,"<<analyticPressureHoop.rounded()
        <<",pressure_hoop_reference,"<<refPressureHoop<<",pressure_vexact_full,"<<analyticPressureMeridional.rounded()+analyticPressureHoop.rounded()<<",pressure_full_reference,0"
        <<",pressure_rawS_meridional,"<<rawPressureMeridional.rounded()<<",pressure_rawS_hoop,"<<rawPressureHoop.rounded()<<",pressure_rawS_full,"<<rawPressureMeridional.rounded()+rawPressureHoop.rounded()
        <<",native_C_pressure_work,"<<cpWork.rounded()<<",native_C_pressure_reference,0,reference_only_rawS_error_not_suppressed,1,no_extra_hoop_body_force,1\n";
    check(charge.integratedCharge()==q&&field.chargeCandidate().integratedCharge()==q,"readonly smooth diagnostic changed signed Q");
}
void run() {
    scene(false);scene(true);
    std::cout<<"HELMHOLTZ_SMOOTH_SIGNEDQ_DIAGNOSTIC_COMPLETE,scene_count,2,base_N,32,layer_count,1,XY_RZ_nonzero_hanging_support,1"
        <<",rc0_means_structural_and_literal_controls_only,1,physical_accuracy_threshold,NONE,MMS_convergence_qualified,0,production_qualified,0,old_Green_equivalence,0,no_NS_no_App_no_phase_or_charge_update\n";
}
void runConvergence() {
    for(int n:{32,64,128}){scene(false,n,true);scene(true,n,true);}
    std::cout<<"HELMHOLTZ_SMOOTH_SIGNEDQ_CONVERGENCE_DIAGNOSTIC_COMPLETE,scene_count,6,resolution_count,3,base_N_values,32|64|128,local_refine_level,1"
        <<",rc0_means_structural_and_literal_controls_only,1,physical_accuracy_threshold,NONE,observed_orders,POSTPROCESS_ONLY,smooth_only,1,sharp_material_qualified,0,physical_boundary_qualified,0,production_qualified,0,old_Green_equivalence,0,no_NS_no_App_no_phase_or_charge_update\n";
}
} // namespace HelmholtzSmoothSignedQDiagnostic
// END HELMHOLTZ_SMOOTH_SIGNEDQ_N32_DIAGNOSTIC

// BEGIN HELMHOLTZ_OBLIQUE_SHEETQ_N32_DIAGNOSTIC: test-only, no physical PASS gate.
namespace HelmholtzObliqueSheetQDiagnostic {
namespace R=HelmholtzSourceResearch;namespace E=LAMP::Electric;namespace I=LAMP::Interface;
namespace FV=LAMP::Geometry::Detail::FiniteVolume;namespace Cut=LAMP::Geometry::PlanarCellCut;
using Clock=LAMP::TimeLayer;using Matrix=HydroStrainOperator::Matrix;
constexpr int baseN=32;
struct RationalCell {long long liquid=0,gas=0,length=0;};
long long scaled(double x,int coordinateScale) {
    const long double value=static_cast<long double>(x)*coordinateScale;const auto integer=static_cast<long long>(value);
    check(value==integer,"oblique fixture coordinates left the frozen exact 1/(2*N) lattice");return integer;
}
long long positiveSquare(long long x) {return x>0?x*x:0;}
RationalCell rationalCell(const Cut::Rectangle& b,int coordinateScale) {
    const long long x0=scaled(b.lower.x,coordinateScale),x1=scaled(b.upper.x,coordinateScale),z0=scaled(b.lower.y,coordinateScale),z1=scaled(b.upper.y,coordinateScale);
    const long long a=3*(x1-x0),d=4*(z1-z0),w=7LL*coordinateScale/2-3*x0-4*z0;
    const auto cdf=[&](long long v){return positiveSquare(v)-positiveSquare(v-a)-positiveSquare(v-d)+positiveSquare(v-a-d);};
    // Independent closed-form CDF of 3*x+4*z, not a copy of polygon clipping.
    const long long liquid=cdf(w),gas=cdf(a+d-w);
    check(liquid>=0&&gas>=0&&liquid+gas==24*(x1-x0)*(z1-z0),"oblique rational two-phase area partition failed");
    const long long lo=std::max(3*x0,7LL*coordinateScale/2-4*z1),hi=std::min(3*x1,7LL*coordinateScale/2-4*z0);
    const long long width=std::max(0LL,hi-lo);
    // delta-x=width/(3*coordinateScale), ds=(5/4)*delta-x: exact rational length, no sqrt.
    return {liquid,gas,5*width};
}
struct Exact {long double potential=0;std::array<long double,2> electric{};};
Exact exact(Vector2D p) {
    const long double x=static_cast<long double>(p.x)-.5L,z=static_cast<long double>(p.y)-.5L;
    if(3*x+4*z<=0)return {(23*x-36*z)/80,{-23.L/80,9.L/20}};
    return {(-13*x-34*z)/20,{13.L/20,17.L/10}};
}
std::array<long double,2> analyticExterior(const QuadtreeMesh& mesh,const HydroDofLayout& layout) {
    // Independently integrated physical tensor, ordered xx,zz,xz. The line
    // endpoints are (0,7/8),(1,1/8); no actual Field data supplies this oracle.
    const std::array<long double,3> tl{-767.L/3200,767.L/3200,-207.L/400},tg{-987.L/800,987.L/800,221.L/200};
    std::array<FV::Dyadic,2> power;
    for(int row=0;row<layout.domain().nativeFaceCount();++row){const auto f=layout.domain().face(mesh,row);
        if(f.negativeCell>=0&&f.positiveCell>=0)continue;check(f.area>0,"XY oblique fixture has an unexpected zero-area exterior");
        const auto& edge=f.horizontal?mesh.hEdges(f.edgeLabel):mesh.vEdges(f.edgeLabel);const auto a=mesh.verts(edge.edgeToVert[0]).coord(),b=mesh.verts(edge.edgeToVert[1]).coord();
        const long double lo=f.horizontal?std::min(a.x,b.x):std::min(a.y,b.y),hi=f.horizontal?std::max(a.x,b.x):std::max(a.y,b.y);
        long double liquid=0;if(f.horizontal){check(a.y==0||a.y==1,"oblique exterior oracle encountered a foreign horizontal boundary");liquid=a.y==0?hi-lo:0;}
        else {check(a.x==0||a.x==1,"oblique exterior oracle encountered a foreign vertical boundary");const long double cut=a.x==0?7.L/8:1.L/8;liquid=std::max(0.L,std::min(hi,cut)-lo);}
        const long double gas=hi-lo-liquid,sign=f.negativeCell>=0?1.L:-1.L;
        check(hi-lo==f.area,"oblique boundary oracle does not cover actual native face measure");
        for(int k=0;k<2;++k){const int entry=f.horizontal?(k==0?2:1):(k==0?0:2);
            power[k].add(FV::Dyadic::product({sign,liquid,tl[entry]}));power[k].add(FV::Dyadic::product({sign,gas,tg[entry]}));}
    }
    return {power[0].rounded(),power[1].rounded()};
}
void scene(int n,bool convergence,bool compareGreen=false) {
    check(n==32||n==64||n==128,"oblique diagnostic resolution is outside the frozen same-family scope");
    const int coordinateScale=2*n,expectedLeaves=17*n*n/8;
    const long long areaDenominator=24LL*coordinateScale*coordinateScale,lengthDenominator=12LL*coordinateScale;
    auto dictionary=input(n,0);dictionary.maxCellNum=std::max(50000,8*n*n);
    QuadtreeMesh mesh(dictionary);DictionaryMesh::RefineBox box;box.min={.625,0};box.max={1,1};box.level=1;mesh.applyRefineBoxes({box});
    FaceSegmentGraph graph;graph.build(mesh);const auto domain=NativeLeafDomain::build(mesh,graph,false);const auto layout=HydroDofLayout::build(mesh,domain);
    check(mesh.cellIndexNum()==expectedLeaves&&layout.pressureCount()==expectedLeaves&&!domain.uniformLeaves(),"oblique frozen actual 17*N*N/8 2:1 leaf family changed");
    const auto divergence=HydroIntegralDivergence::build(mesh,layout);const auto metrics=HydroPhysicalDivergence::build(mesh,divergence);const auto pressure=HydroPressureLoad::build(divergence);
    using Stock=std::array<long double,2>;std::vector<Stock> stocks(mesh.cellIndexNum());std::vector<double> q(mesh.cellIndexNum());std::vector<unsigned char> endpoint(mesh.cellIndexNum());
    FV::Dyadic inputQ,inputL,inputG,qStorageAbsolute;long long sumAreaL=0,sumAreaG=0,sumLength=0;
    std::size_t expectedMixed=0,endpointCells=0;double areaInputDefect=0,lengthInputDefect=0,capacityDefect=0,qTotalDefect=0;
    for(int c:domain.cellIndexByRow()){const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);
        const Cut::Rectangle b{mesh.verts(cell.cellToVert[0]).coord(),mesh.verts(cell.cellToVert[2]).coord()};const auto rational=rationalCell(b,coordinateScale);
        const std::vector<Cut::Point> polygon{{b.lower.x-.5L,b.lower.y-.5L},{b.upper.x-.5L,b.lower.y-.5L},{b.upper.x-.5L,b.upper.y-.5L},{b.lower.x-.5L,b.upper.y-.5L}};
        stocks[c]={Cut::polygonMoments(Cut::clipConvexPolygon(polygon,{3,4},0)).area,Cut::polygonMoments(Cut::clipConvexPolygon(polygon,{-3,-4},0)).area};
        const long double volume=domain.cellVolumes()[domain.cellRowByIndex()[c]];
        for(int s=0;s<2;++s){const long long numerator=s==0?rational.liquid:rational.gas;auto error=FV::Dyadic::product({stocks[c][s],static_cast<long double>(areaDenominator)});error.add(-static_cast<long double>(numerator));
            MaxwellLoadContract::bounded(error.rounded(),volume*static_cast<long double>(areaDenominator),areaInputDefect,"public oblique stock clip disagrees with independent rational CDF");}
        auto closure=FV::Dyadic(stocks[c][0]);closure.add(stocks[c][1]);closure.add(-volume);
        MaxwellLoadContract::bounded(closure.rounded(),volume,capacityDefect,"oblique public two-phase stock capacity gate");
        q[c]=static_cast<double>(static_cast<long double>(rational.length)/static_cast<long double>(lengthDenominator));
        auto error=FV::Dyadic::product({q[c],static_cast<long double>(lengthDenominator)});error.add(-static_cast<long double>(rational.length));qStorageAbsolute.add(error.absolute());
        MaxwellLoadContract::bounded(error.rounded(),static_cast<long double>(rational.length),lengthInputDefect,"oblique independent per-cell rational Q RNE input gate");
        const bool mixed=rational.liquid>0&&rational.gas>0;check((rational.length>0)==mixed,"oblique rational interface length does not agree with actual mixed area");
        check((stocks[c][0]>0&&stocks[c][1]>0)==mixed,"oblique public clip lost a real mixed cell or manufactured one");if(mixed)++expectedMixed;
        endpoint[c]=mixed&&(b.lower.x==0||b.upper.x==1);if(endpoint[c])++endpointCells;
        inputQ.add(q[c]);inputL.add(stocks[c][0]);inputG.add(stocks[c][1]);sumAreaL+=rational.liquid;sumAreaG+=rational.gas;sumLength+=rational.length;
    }
    check(expectedMixed>0&&endpointCells==2&&sumAreaL==areaDenominator/2&&sumAreaG==areaDenominator/2&&sumLength==5*lengthDenominator/4,
        "oblique independent rational total area/length or boundary endpoint cover failed");
    auto totalQError=inputQ;totalQError.add(-1.25L);MaxwellLoadContract::bounded(totalQError.rounded(),2.5L,qTotalDefect,"oblique total native Q differs from independent 5/4 reference");
    const Clock initialClock{0,0,0,Clock::Position::Initial},clock{1,0,.01,Clock::Position::Midpoint};
    const auto initial=I::GeometricVofPublication::prepareLocalPlicInventory(mesh,domain,initialClock,stocks,{});
    const auto pub=I::GeometricVofPublication::prepareLocalPlicInventory(mesh,domain,clock,stocks,{});const auto phase=HydroPhaseGeometry::build(mesh,metrics,pub);
    const auto charge=E::ElectricVofCharge::fromIntegrated(mesh,initial,q);std::vector<E::ElectricVofBoundary> boundaries(domain.nativeFaceCount());
    std::size_t publishedMixed=0,resolvedMixed=0,exteriorFaces=0;
    for(int c:domain.cellIndexByRow()){check(pub.plicCells()[c].inventory==stocks[c],"oblique publication changed independent public stock inputs");
        if(stocks[c][0]>0&&stocks[c][1]>0){++publishedMixed;if(pub.plicCells()[c].resolved)++resolvedMixed;}}
    check(publishedMixed==expectedMixed,"oblique real published mixed-cell count differs from rational geometry");
    for(int row=0;row<domain.nativeFaceCount();++row){const auto f=domain.face(mesh,row);const bool internal=f.negativeCell>=0&&f.positiveCell>=0;
        boundaries[row]=internal?E::ElectricVofBoundary{E::ElectricVofBoundary::Type::Interior,0}:
            E::ElectricVofBoundary{E::ElectricVofBoundary::Type::Dirichlet,static_cast<double>(exact(f.center).potential)};if(!internal)++exteriorFaces;}
    const E::ElectricVofMaterials materials{4,1,0,0};const auto field=E::ElectricVofField::prepareOhmic(mesh,pub,charge,clock,materials,boundaries);field.validateFor(mesh,pub,charge,clock);
    check(R::Gas::sameClock(field.geometryTime(),clock)&&R::Gas::sameClock(field.candidateTime(),clock)&&field.geometry().identity()==pub.identity()&&field.sourceCharge().sameInstance(charge)
        &&field.chargeCandidate().integratedCharge()==q,"oblique actual Mid Field/source/chargeCandidate owner or full clock mismatch");
    check(field.report().relativeLinearResidual<=1e-11&&field.report().maximumRelativeGaussResidual<=1e-10,"oblique actual Electric linear/Gauss gate changed");
    const auto stress=E::ElectricVofMaxwellStress::prepare(mesh,field);const auto load=HydroMaxwellLoad::prepare(mesh,pressure,phase,stress);const R::ResearchLaw law(mesh,load,boundaries);
    const auto& geometry=load.geometry();check(geometry.hasCanonicalRawStrain(),"oblique actual canonical raw CSR is missing");const auto& rawV=geometry.stressVolumeVelocityMatrix();const auto rawS=geometry.canonicalRawStrain();
    std::size_t hanging=0,sharpHanging=0,materialWitness=0;std::array<std::size_t,2> witnessByComponent{};
    for(const auto& mf:law.materialFaces){const auto f=domain.face(mesh,layout.nativeFace(mf.column));check(f.negativeCell>=0&&f.positiveCell>=0,"oblique material source includes a physical exterior");
        const int lm=mesh.cells(mesh.cellIndexToLabel[f.negativeCell]).level(),lp=mesh.cells(mesh.cellIndexToLabel[f.positiveCell]).level();check(std::abs(lm-lp)<=1,"oblique actual face violates 2:1");
        if(lm==lp)continue;++hanging;if(q[f.negativeCell]==0&&q[f.positiveCell]==0)continue;++sharpHanging;
        if(field.epsilonByCell()[f.negativeCell]!=field.epsilonByCell()[f.positiveCell]&&mf.force.compare(FV::Dyadic{})!=0){++materialWitness;++witnessByComponent[f.horizontal?1:0];}}
    check(hanging>0&&sharpHanging>0&&materialWitness>0,"oblique real hanging child lacks nonzero Q-adjacent delta-epsilon/material source witness");
    const std::array<long double,2> reference{449.L/512,283.L/128};const auto exteriorReference=analyticExterior(mesh,layout),actualExterior=R::physicalBoundaryTensor(mesh,load);
    double exteriorOracleDefect=0;for(int k=0;k<2;++k)MaxwellLoadContract::bounded(exteriorReference[k]-reference[k],2*std::abs(reference[k]),exteriorOracleDefect,"oblique independent exterior tensor integral differs from exact total source reference");
    long double fieldError2=0,fieldNorm2=0,pureFieldError2=0,pureFieldNorm2=0;
    for(const auto& point:geometry.points()){const int c=layout.nativeCell(point.cellDof);const auto sample=stress.sample(c,point.point);const auto wanted=exact(point.point);
        for(int k=0;k<2;++k){const long double e=k==0?sample.electricField.x:sample.electricField.y,de=e-wanted.electric[k],w=point.physicalWeight;
            fieldError2+=w*de*de;fieldNorm2+=w*wanted.electric[k]*wanted.electric[k];if(q[c]==0){pureFieldError2+=w*de*de;pureFieldNorm2+=w*wanted.electric[k]*wanted.electric[k];}}}
    const char* prefix=compareGreen?"HELMHOLTZ_OBLIQUE_GREEN_COMPARISON":(convergence?"HELMHOLTZ_OBLIQUE_CONVERGENCE":"HELMHOLTZ_OBLIQUE_N32");
    const auto resolution=convergence?std::string(",N,")+std::to_string(n):std::string{};
    std::cout<<prefix<<"_FIXTURE,XY,N,"<<n<<",transition,0.625,actual_cells,"<<layout.pressureCount()<<",published_mixed,"<<publishedMixed<<",resolved_mixed,"<<resolvedMixed
        <<",true_hanging_faces,"<<hanging<<",Q_adjacent_hanging_faces,"<<sharpHanging<<",nonzero_hanging_material_faces,"<<materialWitness<<",hanging_material_x_probe_witness,"<<witnessByComponent[0]<<",hanging_material_y_probe_witness,"<<witnessByComponent[1]
        <<",endpoint_cells,"<<endpointCells<<",all_exterior_faces,"<<exteriorFaces<<",actual_Mid_gid,"<<pub.identity().value()<<",source_Q_id,"<<charge.identity()<<",clock,("<<clock.step<<','<<clock.begin<<','<<clock.end<<','<<static_cast<int>(clock.position)<<")"
        <<",public_stock_clip_rational_CDF_relative,"<<areaInputDefect<<",native_stock_capacity_relative,"<<capacityDefect<<",per_cell_Q_rational_RNE_relative,"<<lengthInputDefect
        <<",sum_native_Q,"<<inputQ.rounded()<<",sum_Q_reference,1.25,sum_Q_error,"<<totalQError.rounded()<<",sum_Q_scaled_absolute_storage_error,"<<qStorageAbsolute.rounded()<<",Q_storage_error_denominator,"<<lengthDenominator
        <<",sum_liquid,"<<inputL.rounded()<<",sum_gas,"<<inputG.rounded()<<",independent_area_reference_each,0.5,independent_interface_length,1.25,derived_PLIC_normal_not_injected,1\n";
    std::cout<<prefix<<"_FIELD_DIAGNOSTIC"<<resolution<<",actual_RT0_E_all_quadrature_L2_relative,"<<std::sqrt(fieldError2/fieldNorm2)<<",actual_RT0_E_pure_cells_L2_relative,"<<std::sqrt(pureFieldError2/pureFieldNorm2)
        <<",field_error_is_diagnostic_not_gate,1,continuous_sharp_face_center_Dirichlet,1,boundary_endpoint_vertex_Et_qualified,0,common_tensor_qualified,0\n";
    for(int k=0;k<2;++k){std::vector<double> velocity(layout.velocityCount());std::vector<long double> gross(layout.pressureCount());
        for(int u=0;u<layout.velocityCount();++u){const auto f=domain.face(mesh,layout.nativeFace(u));velocity[u]=(f.horizontal?1:0)==k?1.:0.;
            for(int c:{f.negativeCell,f.positiveCell})if(c>=0)gross[layout.pressureDof(c)]+=std::abs(static_cast<long double>(f.area)*velocity[u]);}
        const HydroNativeValues values(layout,velocity,std::vector<double>(layout.pressureCount()));const auto actualDI=divergence.apply(mesh,values);
        double divergenceDefect=0,valueDefect=0,strainDefect=0;for(int p=0;p<layout.pressureCount();++p)MaxwellLoadContract::bounded(actualDI[p],gross[p],divergenceDefect,"oblique constant native virtual velocity fails actual DI=0");
        const Eigen::Map<const Eigen::VectorXd> native(velocity.data(),velocity.size());const Eigen::VectorXd reconstructed=rawV*native;const auto strain=rawS.apply(mesh,values);
        for(Eigen::Index row=0;row<rawV.rows();++row){const long double wanted=row%2==k?1.L:0.L;FV::Dyadic scale(std::abs(wanted));
            for(Matrix::InnerIterator e(rawV,row);e;++e)scale.add(FV::Dyadic::product({std::abs(e.value()),std::abs(velocity[e.col()])}));
            MaxwellLoadContract::bounded(static_cast<long double>(reconstructed[row])-wanted,scale.rounded(),valueDefect,"oblique actual R_V lost a constant native vector");}
        for(Eigen::Index row=0;row<rawS.matrix().rows();++row){FV::Dyadic scale;for(Matrix::InnerIterator e(rawS.matrix(),row);e;++e)scale.add(FV::Dyadic::product({std::abs(e.value()),std::abs(velocity[e.col()])}));
            MaxwellLoadContract::bounded(strain[row],scale.rounded(),strainDefect,"oblique raw strain of a constant native vector is not zero");}
        const auto work=law.work(velocity);FV::Dyadic endpointVolume,endpointMaterial;
        for(std::size_t qi=0;qi<geometry.points().size();++qi){const int c=layout.nativeCell(geometry.points()[qi].cellDof);if(!endpoint[c])continue;
            for(int axis=0;axis<2;++axis)for(Matrix::InnerIterator e(rawV,static_cast<Eigen::Index>(qi)*2+axis);e;++e)endpointVolume.add(law.volumeFactors[qi][axis].multiplied(e.value()).multiplied(velocity[e.col()]));}
        for(const auto& mf:law.materialFaces){const auto f=domain.face(mesh,layout.nativeFace(mf.column));if(endpoint[f.negativeCell]||endpoint[f.positiveCell])endpointMaterial.add(mf.force.multiplied(velocity[mf.column]));}
        std::cout<<prefix<<"_PHYSICAL_DIAGNOSTIC"<<resolution<<",component,"<<(k==0?"x":"y")<<",native_DI_relative,"<<divergenceDefect<<",R_V_constant_relative,"<<valueDefect<<",raw_strain_zero_relative,"<<strainDefect
            <<",actual_volumeQ,"<<work.volume<<",actual_material_faces,"<<work.material<<",individual_terms_have_no_unique_sharp_reference,1,actual_alternative_native,"<<work.alternativeNative<<",analytic_total_source_reference,"<<reference[k]
            <<",total_source_signed_error,"<<work.alternativeNative-reference[k]<<",total_source_relative_error,"<<std::abs((work.alternativeNative-reference[k])/reference[k])<<",original_literal_tolerance,"<<R::originalWorkTolerance<<",source_literal_relative,"<<work.alternativeRelative
            <<",analytic_exterior_tensor_reference,"<<exteriorReference[k]<<",analytic_exterior_oracle_relative,"<<exteriorOracleDefect<<",actual_RT0_exterior_tensor_diagnostic,"<<actualExterior[k]<<",source_exterior_not_added,1"
            <<",endpoint_cell_volumeQ_diagnostic,"<<endpointVolume.rounded()<<",endpoint_touching_material_diagnostic,"<<endpointMaterial.rounded()<<",endpoint_terms_not_subtracted,1,noncompact_constant_probe,1,physical_status,DIAGNOSTIC_NO_ACCURACY_GATE,physical_PASS_claim,0\n";
        if(compareGreen){const auto& shared=geometry.stressFaceVelocityMatrix();
            check(load.report().pairing==HydroMaxwellPairing::Rt0Green&&shared.rows()==static_cast<Eigen::Index>(geometry.stressFacePoints().size())*2,
                "oblique Green comparison lost actual RT0 pairing/shared trace owner");
            const Eigen::VectorXd trace=shared*native;double traceDefect=0;
            // Check the actual shared trace independently. Raw strain=0 does
            // not certify this matrix's full boundary or hanging-face value.
            for(Eigen::Index row=0;row<shared.rows();++row){const long double wanted=row%2==k?1.L:0.L;FV::Dyadic scale(std::abs(wanted));
                for(Matrix::InnerIterator e(shared,row);e;++e)scale.add(FV::Dyadic::product({std::abs(e.value()),std::abs(velocity[e.col()])}));
                MaxwellLoadContract::bounded(static_cast<long double>(trace[row])-wanted,scale.rounded(),traceDefect,"oblique actual shared face trace lost a constant native vector");}
            const auto green=MaxwellLoadContract::LiteralGreenWork(load).evaluate(mesh,values);
            check(MaxwellLoadContract::nativePower(load,values)==green.momentumWork,"oblique actual RHS exact-dot power differs from independent Green reference");
            std::cout<<"HELMHOLTZ_OBLIQUE_GREEN_COMPARISON_WORK,N,"<<n<<",component,"<<(k==0?"x":"y")<<",actual_Mid_gid,"<<pub.identity().value()<<",actual_shared_trace_constant_relative,"<<traceDefect
                <<",actual_Green_momentum_work,"<<green.momentumWork<<",actual_Green_volume_work,"<<green.volumeWork<<",actual_Green_boundary_work,"<<green.boundaryWork
                <<",actual_Green_literal_residual,"<<green.residual<<",actual_Green_literal_relative,"<<green.relativeResidual<<",actual_Green_literal_absolute_scale,"<<green.absoluteScale
                <<",independent_physical_reference,"<<reference[k]<<",Green_signed_error,"<<green.momentumWork-reference[k]<<",Green_relative_error,"<<std::abs((green.momentumWork-reference[k])/reference[k])
                <<",actual_source_alternative_work,"<<work.alternativeNative<<",source_signed_error,"<<work.alternativeNative-reference[k]<<",source_relative_error,"<<std::abs((work.alternativeNative-reference[k])/reference[k])
                <<",actual_RT0_exterior_tensor_diagnostic,"<<actualExterior[k]<<",Green_minus_actual_RT0_exterior_tensor,"<<green.momentumWork-actualExterior[k]
                <<",source_and_exterior_not_added,1,same_actual_Field_no_extra_Field_or_NS,1,Green_work_is_own_literal_ledger,1,source_common_tensor_claim,0,physical_status,DIAGNOSTIC_NO_ACCURACY_GATE,physical_PASS_claim,0\n";}
    }
    check(charge.integratedCharge()==q&&field.chargeCandidate().integratedCharge()==q,"readonly oblique sheet-Q diagnostic changed any native charge");
    if(!convergence)std::cout<<"HELMHOLTZ_OBLIQUE_SHEETQ_DIAGNOSTIC_COMPLETE,family_count,1,layer_count,1,N,32,actual_cells,2176,constant_probe_count,2,rc0_means_structural_and_literal_controls_only,1"
        <<",physical_accuracy_threshold,NONE,sharp_convergence_qualified,0,Neumann_qualified,0,noSlip_qualified,0,electric_energy_qualified,0,production_qualified,0,old_Green_equivalence,0,no_NS_no_App_no_phase_or_charge_update\n";
}
void run() {scene(baseN,false);}
void runConvergence() {
    for(int n:{32,64,128})scene(n,true);
    std::cout<<"HELMHOLTZ_OBLIQUE_SHEETQ_CONVERGENCE_DIAGNOSTIC_COMPLETE,family_count,1,layer_count,3,N_values,32|64|128,local_refine_level,1,transition,0.625,actual_leaf_rule,17*N*N/8,constant_probe_count_per_layer,2"
        <<",rc0_means_structural_and_literal_controls_only,1,physical_accuracy_threshold,NONE,observed_orders,POSTPROCESS_ONLY,sharp_convergence_qualified,0,Neumann_qualified,0,noSlip_qualified,0,electric_energy_qualified,0,production_qualified,0,old_Green_equivalence,0,no_NS_no_App_no_phase_or_charge_update\n";
}
void runGreenComparison() {
    for(int n:{32,64,128})scene(n,true,true);
    std::cout<<"HELMHOLTZ_OBLIQUE_GREEN_COMPARISON_COMPLETE,family_count,1,layer_count,3,N_values,32|64|128,actual_Field_prepares_per_layer,1,extra_Field_or_NS,0,constant_probe_count_per_layer,2"
        <<",rc0_means_structural_and_literal_controls_only,1,physical_accuracy_threshold,NONE,observed_orders,POSTPROCESS_ONLY,general_physical_PASS,0,Neumann_qualified,0,noSlip_qualified,0,electric_energy_qualified,0,production_force_changed,0,source_common_tensor_claim,0\n";
}
} // namespace HelmholtzObliqueSheetQDiagnostic
// END HELMHOLTZ_OBLIQUE_SHEETQ_N32_DIAGNOSTIC

int main(int argc,char** argv)
{
    try {
        std::cout<<std::setprecision(17);
        if(argc==2&&std::string(argv[1])=="--helmholtz-oblique-green-comparison-only"){HelmholtzObliqueSheetQDiagnostic::runGreenComparison();return 0;}
        if(argc==2&&std::string(argv[1])=="--helmholtz-oblique-sheetq-convergence-only"){HelmholtzObliqueSheetQDiagnostic::runConvergence();return 0;}
        if(argc==2&&std::string(argv[1])=="--helmholtz-oblique-sheetq-diagnostic-only"){HelmholtzObliqueSheetQDiagnostic::run();return 0;}
        if(argc==2&&std::string(argv[1])=="--helmholtz-smooth-signedq-convergence-only"){HelmholtzSmoothSignedQDiagnostic::runConvergence();return 0;}
        if(argc==2&&std::string(argv[1])=="--helmholtz-smooth-signedq-diagnostic-only"){HelmholtzSmoothSignedQDiagnostic::run();return 0;}
        if(argc==2&&std::string(argv[1])=="--helmholtz-source-research-only"){HelmholtzSourceResearch::run();return 0;}
        if(argc==2&&std::string(argv[1])=="--no-slip-strain-only"){StationaryNoSlipStrainContract::run();return 0;}
        if(argc==2&&std::string(argv[1])=="--native-scalar-prefix-only"){NativeScalarPrefixContract::run();return 0;}
        if(argc==2&&std::string(argv[1])=="--dyadic-integer-operations-only"){DyadicIntegerOperationsContract::run();return 0;}
        if(argc==2&&std::string(argv[1])=="--native-small-exact-only"){NativeSmallExactContract::run();return 0;}
        if(argc==2&&std::string(argv[1])=="--dyadic-terminal-rne-only"){DyadicTerminalRneContract::run();return 0;}
        if(argc==2&&std::string(argv[1])=="--viscosity-product-local-memo-only"){ViscosityProductLocalMemoContract::run();return 0;}
        if(argc==2&&std::string(argv[1])=="--viscosity-product-parallel-only"){ViscosityProductParallelContract::run();return 0;}
        if(argc==2&&std::string(argv[1])=="--native-small-product-only"){NativeSmallProductContract::run();return 0;}
        if(argc==2&&std::string(argv[1])=="--binary64-limb-accumulator-only"){Binary64LimbAccumulatorContract::run();return 0;}
        if(argc==2&&std::string(argv[1])=="--gram-geometry-cache-only"){GramGeometryStatisticsContract::run();return 0;}
        if(argc==2&&std::string(argv[1])=="--viscosity-product-plan-only"){ViscosityProductPlanContract::run();return 0;}
        if(argc==2&&std::string(argv[1])=="--maxwell-work-preaggregate-only"){MaxwellLoadContract::runPreaggregated();return 0;}
        if(argc==2&&std::string(argv[1])=="--maxwell-pure-gas-solenoidal-only"){MaxwellPureGasSolenoidalContract::run();return 0;}
        if(argc==2&&std::string(argv[1])=="--maxwell-load-only"){MaxwellLoadContract::run();return 0;}
        if(argc==2&&std::string(argv[1])=="--stress-compatible-only"){StressCompatibleContract::run();return 0;}
        if(argc==2&&std::string(argv[1])=="--stress-raw-volume-csr-rebind-only"){StressAmrLocalReuseContract::RawVolumeCsrRebind::enabled=true;StressAmrLocalReuseContract::CertifiedCsrRebind::enabled=true;StressAmrLocalReuseContract::run();StressAmrLocalReuseContract::CertifiedCsrRebind::finish();StressAmrLocalReuseContract::RawVolumeCsrRebind::finish();return 0;}
        if(argc==2&&std::string(argv[1])=="--stress-certified-csr-rebind-only"){StressAmrLocalReuseContract::CertifiedCsrRebind::enabled=true;StressAmrLocalReuseContract::run();StressAmrLocalReuseContract::CertifiedCsrRebind::finish();return 0;}
        if(argc==2&&std::string(argv[1])=="--stress-selected-amr-reuse-only"){SelectedStressAmrReuseContract::run();return 0;}
        if(argc==2&&std::string(argv[1])=="--stress-amr-local-reuse-only"){StressAmrLocalReuseContract::run();return 0;}
        if(argc==2&&std::string(argv[1])=="--gram-coercivity-only"){GramCoercivityContract::run();return 0;}
        if(argc==2&&std::string(argv[1])=="--maxwell-compiled-contraction-only"){MaxwellCompiledContract::run();return 0;}
        for (int metric : {0,1,2}) for (int n : {4,8}) for (int depth : {0,1}) algebra(n,metric,depth);
        for (int metric : {0,1,2}) algebra(4,metric,2);
        for (int metric : {0,1}) for (int selection : {1,2,3}) algebra(12,metric,1,selection);
        for (bool rz : {false,true}) for (bool amr : {false,true}) mms(rz,amr);
        axis();guards();
        std::cout<<"SUMMARY,"<<algebraCases<<','<<polynomialCases<<','<<spectra<<','<<scalings<<','<<mmsCases<<','<<axisLimits<<','<<saddles<<','<<rejected
            <<','<<peakPolynomial<<','<<peakSymmetry<<','<<peakEnergy<<','<<peakRigid<<','<<peakScaling<<','<<peakSaddle<<','<<peakQuadrature<<','<<peakCentral<<'\n';return 0;
    } catch (const std::exception& e) {std::cerr<<"HYDRO energy viscosity FAIL: "<<e.what()<<'\n';return 1;}
}
