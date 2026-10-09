#include "Modules/Hydrodynamics/CchEnergyViscosity.h" // Public header is self-contained.
#include "Modules/Hydrodynamics/CchSinglePhaseSaddle.h"
#include "Core/IO/DictionaryMesh.h"
#include "Core/Mesh/FaceSegmentGraph.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <array>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <string>

namespace
{
using namespace Quantum::Hydrodynamics;
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
NativeReference reference(const QuadtreeMesh& mesh,const FaceSegmentGraph& graph,const CchDofLayout& layout)
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
struct Field {Vector2D velocity;std::array<double,4> strain;};
Field polynomial(Vector2D p,bool rz,int type)
{
    const double x=p.x,y=p.y;
    if (type==1) return rz ? Field{{0,2},{0,0,0,0}} : Field{{2,-3},{0,0,0,0}};
    if (type==2 && !rz) return {{-y,x},{0,0,0,0}};
    if (type==2 && rz) {
        auto f=polynomial(p,false,0);f.strain[3]=f.strain[2];f.strain[2]=f.velocity.x/x;return f;
    }
    if (rz) return {{x*(1+2*y),3+x*x+y+2*y*y},{1+2*y,1+4*y,1+2*y,4*x/std::sqrt(2.0)}};
    return {{1+x+2*y+3*x*x+4*x*y+5*y*y,-1+2*x-y+4*x*x+3*x*y+2*y*y},
        {1+6*x+4*y,-1+3*x+4*y,(4+12*x+13*y)/std::sqrt(2.0),0}};
}
Field smooth(Vector2D p,bool rz)
{
    const double x=p.x,y=p.y,sx=std::sin(pi*x),cx=std::cos(pi*x),sy=std::sin(pi*y),cy=std::cos(pi*y);
    if (rz) return {{x*sy,cx*sy},{sy,pi*cx*cy,sy,(pi*x*cy-pi*sx*sy)/std::sqrt(2.0)}};
    return {{sx*cy,cx*sy},{pi*cx*cy,pi*cx*cy,-std::sqrt(2.0)*pi*sx*sy,0}};
}
template<class F> CchNativeValues values(const CchDofLayout& layout,const NativeReference& ref,F field)
{
    std::vector<double> u(layout.velocityCount());
    for (int f=0;f<layout.velocityCount();++f) {const auto v=field(ref.centers[f]).velocity;u[f]=layout.normalAxis(f)==CchDofLayout::NormalAxis::First ? v.x : v.y;}
    return CchNativeValues(layout,u,std::vector<double>(layout.pressureCount()));
}
double dot(const std::vector<double>& a,const std::vector<double>& b)
{
    long double result=0;for (std::size_t i=0;i<a.size();++i) result+=static_cast<long double>(a[i])*b[i];return static_cast<double>(result);
}
void saddle(const QuadtreeMesh& mesh,const CchPhysicalDivergence& metrics,const CchEnergyViscosity& k,const NativeReference& ref)
{
    const auto& l=metrics.layout();const bool rz=l.domain().axisymmetric();
    std::vector<double> u(l.velocityCount()),p(l.pressureCount()),prescribed(l.velocityCount());
    for (int f=0;f<l.velocityCount();++f) {
        const auto x=ref.centers[f];const auto v=rz ? Vector2D{x.x*x.y,-x.y*x.y} : Vector2D{x.x*x.x,-2*x.x*x.y};
        u[f]=l.normalAxis(f)==CchDofLayout::NormalAxis::First ? v.x : v.y;
        if (ref.incidence[f]==1 || metrics.faceMeasures()[f]==0) prescribed[f]=u[f];
    }
    for (int c=0;c<l.pressureCount();++c) p[c]=mesh.cells(mesh.cellIndexToLabel[l.nativeCell(c)]).centerCoord().y-.5;
    CchSaddleVelocityBlock::Matrix a=k.matrix();for (int f=0;f<l.velocityCount();++f) a.coeffRef(f,f)+=1;a.makeCompressed();
    const auto pressure=CchPressureLoad::build(metrics.divergence());
    const Eigen::Map<const Eigen::VectorXd> eu(u.data(),u.size()),ep(p.data(),p.size());const Eigen::VectorXd rhs=a*eu+pressure.matrix()*ep;
    const auto solve=CchSinglePhaseSaddle::buildClosed(pressure,CchSaddleVelocityBlock::build(l,a));
    const auto result=solve.solve(mesh,CchSaddleRhs(l,{rhs.data(),rhs.data()+rhs.size()},prescribed));
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
    const auto l=CchDofLayout::build(mesh,NativeLeafDomain::build(mesh,graph,metric!=0));
    const auto metrics=CchPhysicalDivergence::build(mesh,CchIntegralDivergence::build(mesh,l));const auto ref=reference(mesh,graph,l);
    const auto points=CchStrainOperator::volumeQuadrature(mesh,metrics);const auto e=CchStrainOperator::build(mesh,metrics,points);
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
    const auto k=CchEnergyViscosity::build(e,mu);
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
    const CchNativeValues v(l,random,std::vector<double>(l.pressureCount()));const auto strain=e.apply(mesh,v);
    if (depth==0 && metric!=1) for (int p=0;p<l.pressureCount();++p) {
        const auto& cell=mesh.cells(mesh.cellIndexToLabel[l.nativeCell(p)]);const auto c=cell.centerCoord(),w=cell.size();
        std::array<int,4> face{-1,-1,-1,-1};
        for (int f=0;f<l.velocityCount();++f) {
            const auto x=ref.centers[f];
            if (l.normalAxis(f)==CchDofLayout::NormalAxis::First && x.y==c.y) {if (x.x==c.x-w.x/2) face[0]=f;if (x.x==c.x+w.x/2) face[1]=f;}
            if (l.normalAxis(f)==CchDofLayout::NormalAxis::Second && x.x==c.x) {if (x.y==c.y-w.y/2) face[2]=f;if (x.y==c.y+w.y/2) face[3]=f;}
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
    const auto full=CchEnergyViscosity::build(e,std::vector<double>(mu.size(),1));
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
        const auto scaled=CchEnergyViscosity::build(e,std::vector<double>(mu.size(),scale));const auto force=scaled.apply(mesh,v);
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
        const auto l=CchDofLayout::build(mesh,NativeLeafDomain::build(mesh,graph,rz));const auto metrics=CchPhysicalDivergence::build(mesh,CchIntegralDivergence::build(mesh,l));
        const auto points=CchStrainOperator::volumeQuadrature(mesh,metrics);const auto e=CchStrainOperator::build(mesh,metrics,points);
        const auto v=values(l,reference(mesh,graph,l),[rz](Vector2D p) {return smooth(p,rz);});const auto s=e.apply(mesh,v);
        long double error=0,norm=0;for (std::size_t q=0;q<points.size();++q) {
            const auto exact=smooth(points[q].point,rz).strain;
            for (int c=0;c<e.components();++c) {const double diff=s[q*e.components()+c]-exact[c];error+=points[q].physicalWeight*diff*diff;norm+=points[q].physicalWeight*exact[c]*exact[c];}
        }
        const double strainError=std::sqrt(static_cast<double>(error/norm));
        const auto k=CchEnergyViscosity::build(e,std::vector<double>(points.size(),1));const double energyError=std::abs(k.dissipation(mesh,v)-exactEnergy)/exactEnergy;
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
        const auto l=CchDofLayout::build(mesh,NativeLeafDomain::build(mesh,graph,true));const auto metrics=CchPhysicalDivergence::build(mesh,CchIntegralDivergence::build(mesh,l));
        const auto ref=reference(mesh,graph,l);const auto v=values(l,ref,[](Vector2D p) {return polynomial(p,true,0);});
        std::vector<CchStrainPoint> points;
        for (int p=0;p<l.pressureCount();++p) {const auto& c=mesh.cells(mesh.cellIndexToLabel[l.nativeCell(p)]);if (c.centerCoord().x-c.size().x/2==0) points.push_back({p,{0,c.centerCoord().y},0});}
        const auto e=CchStrainOperator::build(mesh,metrics,points);const auto s=e.apply(mesh,v);
        for (std::size_t q=0;q<points.size();++q) {close(s[4*q+2]-(1+2*points[q].point.y),4,peakPolynomial);close(s[4*q+3],1,peakPolynomial);check(s[4*q]==s[4*q+2],"hoop limit differs from radial strain");}
        ++axisLimits;
    }
}
void guards()
{
    QuadtreeMesh mesh(input(4,1));FaceSegmentGraph graph;graph.build(mesh);
    const auto l=CchDofLayout::build(mesh,NativeLeafDomain::build(mesh,graph,true));const auto metrics=CchPhysicalDivergence::build(mesh,CchIntegralDivergence::build(mesh,l));
    const auto points=CchStrainOperator::volumeQuadrature(mesh,metrics);const auto e=CchStrainOperator::build(mesh,metrics,points);const auto k=CchEnergyViscosity::build(e,std::vector<double>(points.size(),1));
    const CchNativeValues zero(l,std::vector<double>(l.velocityCount()),std::vector<double>(l.pressureCount()));const auto accepted=k.apply(mesh,zero);
    rejects([&]{CchStrainOperator::build(mesh,metrics,{});});
    for (int kind=0;kind<5;++kind) {auto q=points;
        if (kind==0) q[0].cellDof=-1;
        if (kind==1) q[0].point.x=-1;
        if (kind==2) q[0].physicalWeight=-1;
        if (kind==3) q[0].point.y=std::numeric_limits<double>::quiet_NaN();
        if (kind==4) q[0].physicalWeight=std::numeric_limits<double>::infinity();
        rejects([&]{CchStrainOperator::build(mesh,metrics,q);});
    }
    rejects([&]{CchEnergyViscosity::build(e,{});});
    for (double bad : {-1.,std::numeric_limits<double>::quiet_NaN(),std::numeric_limits<double>::infinity()}) {std::vector<double> mu(points.size(),1);mu[0]=bad;rejects([&]{CchEnergyViscosity::build(e,mu);});}
    auto overflow=points;overflow[0].physicalWeight=std::numeric_limits<double>::max();const auto largeE=CchStrainOperator::build(mesh,metrics,overflow);
    rejects([&]{CchEnergyViscosity::build(largeE,std::vector<double>(points.size(),1));});
    auto axisU=zero.normalVelocity();axisU[metrics.axisVelocityDofs()[0]]=std::numeric_limits<double>::denorm_min();const CchNativeValues badAxis(l,axisU,zero.pressure());
    rejects([&]{e.apply(mesh,badAxis);});rejects([&]{k.apply(mesh,badAxis);});rejects([&]{k.dissipation(mesh,badAxis);});
    const auto other=CchDofLayout::build(mesh,NativeLeafDomain::build(mesh,graph,true));const CchNativeValues wrong(other,zero.normalVelocity(),zero.pressure());
    rejects([&]{e.apply(mesh,wrong);});rejects([&]{k.apply(mesh,wrong);});
    int largestColumn=0;double largest=0;
    for (Eigen::Index row=0;row<e.matrix().outerSize();++row) for (CchStrainOperator::Matrix::InnerIterator entry(e.matrix(),row);entry;++entry)
        if (std::abs(entry.value())>largest) {largest=std::abs(entry.value());largestColumn=static_cast<int>(entry.col());}
    check(largest>1,"overflow fixture has no amplifying derivative");
    auto huge=zero.normalVelocity();huge[largestColumn]=std::numeric_limits<double>::max();const CchNativeValues hugeValues(l,huge,zero.pressure());
    rejects([&]{e.apply(mesh,hugeValues);});rejects([&]{k.apply(mesh,hugeValues);});rejects([&]{k.dissipation(mesh,hugeValues);});
    const auto inviscid=CchEnergyViscosity::build(e,std::vector<double>(points.size()));
    check(inviscid.dissipation(mesh,hugeValues)==0,"zero viscosity evaluated an overflowing strain");
    for (double f : inviscid.apply(mesh,hugeValues)) check(f==0,"inviscid extreme velocity force");
    QuadtreeMesh stretched(input(4,1,2));rejects([&]{e.apply(stretched,zero);});rejects([&]{k.apply(stretched,zero);});rejects([&]{CchStrainOperator::build(stretched,metrics,points);});
    QuadtreeMesh tiny(input(1,0));FaceSegmentGraph tg;tg.build(tiny);const auto tl=CchDofLayout::build(tiny,NativeLeafDomain::build(tiny,tg,false));
    const auto tm=CchPhysicalDivergence::build(tiny,CchIntegralDivergence::build(tiny,tl));rejects([&]{CchStrainOperator::build(tiny,tm,CchStrainOperator::volumeQuadrature(tiny,tm));});
    refine(mesh,1);rejects([&]{k.apply(mesh,zero);});rejects([&]{CchStrainOperator::volumeQuadrature(mesh,metrics);});
    QuadtreeMesh equivalent(input(4,1));check(k.apply(equivalent,zero)==accepted,"failed operation changed source");
    std::optional<CchEnergyViscosity> retained;std::optional<CchNativeValues> rv;
    {QuadtreeMesh temporary(input(4,1));FaceSegmentGraph g;g.build(temporary);const auto dl=CchDofLayout::build(temporary,NativeLeafDomain::build(temporary,g,true));
     const auto m=CchPhysicalDivergence::build(temporary,CchIntegralDivergence::build(temporary,dl));const auto s=CchStrainOperator::build(temporary,m,CchStrainOperator::volumeQuadrature(temporary,m));
     retained=CchEnergyViscosity::build(s,std::vector<double>(s.points().size(),1));rv.emplace(dl,std::vector<double>(dl.velocityCount()),std::vector<double>(dl.pressureCount()));}
    check(retained->apply(equivalent,*rv)==accepted,"borrowed mesh lifetime");
}
}
int main()
{
    try {
        std::cout<<std::setprecision(17);
        for (int metric : {0,1,2}) for (int n : {4,8}) for (int depth : {0,1}) algebra(n,metric,depth);
        for (int metric : {0,1,2}) algebra(4,metric,2);
        for (int metric : {0,1}) for (int selection : {1,2,3}) algebra(12,metric,1,selection);
        for (bool rz : {false,true}) for (bool amr : {false,true}) mms(rz,amr);
        axis();guards();
        std::cout<<"SUMMARY,"<<algebraCases<<','<<polynomialCases<<','<<spectra<<','<<scalings<<','<<mmsCases<<','<<axisLimits<<','<<saddles<<','<<rejected
            <<','<<peakPolynomial<<','<<peakSymmetry<<','<<peakEnergy<<','<<peakRigid<<','<<peakScaling<<','<<peakSaddle<<','<<peakQuadrature<<','<<peakCentral<<'\n';return 0;
    } catch (const std::exception& e) {std::cerr<<"CCH energy viscosity FAIL: "<<e.what()<<'\n';return 1;}
}
