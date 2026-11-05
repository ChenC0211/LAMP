#include "Modules/Hydrodynamics/HydroSinglePhaseSaddle.h"
#include "Modules/Hydrodynamics/HydroEnergyViscosity.h"
#include "Modules/Hydrodynamics/HydroMomentumGeometry.h"
#include "Modules/Hydrodynamics/Detail/HydroPressureCoarseCorrection.h"
#include "Solvers/Linear/ParallelBlockJacobi.h"
#include <cfenv>
#include <cstdint>
#include <cstring>
#include <cerrno>
#if defined(_MSC_VER) && defined(_M_X64)
#include <float.h>
#include <xmmintrin.h>
#endif
#include <Eigen/SparseLU>
#include <Eigen/LU>
#include "Core/IO/DictionaryMesh.h"
#include "Core/Mesh/FaceSegmentGraph.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <functional>
#include <iostream>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>
#include <exception>
#ifdef _OPENMP
#include <omp.h>
#endif


namespace
{
using namespace LAMP::Hydrodynamics;
using Matrix=HydroSaddleVelocityBlock::Matrix;
const double pi=std::acos(-1.0);
int algebraCases=0,mmsCases=0,rejected=0;
double peakVelocity=0,peakPressure=0,peakEnergy=0,peakContinuity=0,peakAlgebraic=0,peakMean=0;
double velocityOrder=0,pressureOrder=0;
void check(bool ok,const char* reason) {if (!ok) throw std::runtime_error(reason);}
template<class F> void rejects(F action,const char* stage=nullptr) {
    try {action();} catch (const std::exception& error) {
        if (stage) check(std::string(error.what()).find(stage)!=std::string::npos,"wrong rejection stage");
        ++rejected;return;
    }
    throw std::runtime_error("invalid saddle input accepted");
}
DictionaryMesh input(int n,int metric) {
    DictionaryMesh data;data.Nx=n;data.Ny=n;data.minCoord={metric==2 ? .25 : 0,0};data.maxCoord={metric==2 ? 1.25 : 1,1};
    data.maxInitialIterations=0;data.maxCellNum=10000;return data;
}
void refine(QuadtreeMesh& mesh,int depth) {
    const auto width=mesh.maxCoord-mesh.minCoord;std::vector<DictionaryMesh::RefineBox> boxes;
    for (int level=1;level<=depth;++level) {
        DictionaryMesh::RefineBox box;const double lo=level==1 ? .25 : .375;
        box.min=mesh.minCoord+width*lo;box.max=mesh.minCoord+width*(1-lo);box.level=level;boxes.push_back(box);
    }
    if (!boxes.empty()) mesh.applyRefineBoxes(boxes);
}
FaceSegmentGraph selected(const QuadtreeMesh& mesh,int selection) {
    FaceSegmentGraph graph;graph.build(mesh,[selection](const QuadtreeMesh& m,const Cell& cell) {
        const auto p=cell.centerCoord();const double x=(p.x-m.minCoord.x)/(m.maxCoord.x-m.minCoord.x);
        if (selection==1) return x<.25 || x>.75;
        if (selection==2) return !(x>.25 && x<.75 && p.y>.25 && p.y<.75);
        return true;
    });return graph;
}
struct Face {double x0=0,x1=0,y0=0,y1=0,area=0;bool horizontal=false;int incidence=0;};
struct Term {int p,u;double outwardArea;};
struct Reference {std::vector<Face> faces;std::vector<double> volume;std::vector<Vector2D> center;std::vector<Term> terms;};
Reference reference(const QuadtreeMesh& mesh,const FaceSegmentGraph& graph,const HydroDofLayout& layout,bool rz) {
    Reference ref;ref.faces.resize(layout.velocityCount());ref.volume.resize(layout.pressureCount());ref.center.resize(layout.pressureCount());
    for (int p=0;p<layout.pressureCount();++p) {
        const auto& cell=mesh.cells(mesh.cellIndexToLabel[layout.nativeCell(p)]);const auto c=cell.centerCoord(),w=cell.size();
        ref.center[p]=c;ref.volume[p]=rz ? pi*((c.x+w.x/2)*(c.x+w.x/2)-(c.x-w.x/2)*(c.x-w.x/2))*w.y : w.x*w.y;
    }
    for (const auto& segment : graph.segments()) {
        const auto& edge=segment.isHorizontal ? mesh.hEdges(segment.edgeLabel) : mesh.vEdges(segment.edgeLabel);
        const auto a=mesh.verts(edge.edgeToVert[0]).coord(),b=mesh.verts(edge.edgeToVert[1]).coord();
        const int native=segment.isHorizontal ? mesh.vEdgeIndexNum()+edge.index() : edge.index();
        const int u=layout.velocityDof(native),p=layout.pressureDof(mesh.cells(segment.K).index());
        auto& f=ref.faces[u];f.x0=std::min(a.x,b.x);f.x1=std::max(a.x,b.x);f.y0=std::min(a.y,b.y);f.y1=std::max(a.y,b.y);
        f.horizontal=segment.isHorizontal;
        f.area=rz ? (f.horizontal ? pi*(f.x1*f.x1-f.x0*f.x0) : 2*pi*f.x0*(f.y1-f.y0)) : (f.horizontal ? f.x1-f.x0 : f.y1-f.y0);
        ++f.incidence;ref.terms.push_back({p,u,f.area*(f.horizontal ? segment.n.y : segment.n.x)});
    }
    return ref;
}
Matrix diagonal(const std::vector<double>& values) {
    std::vector<Eigen::Triplet<double,Eigen::Index>> entries;
    for (std::size_t i=0;i<values.size();++i) entries.emplace_back(static_cast<Eigen::Index>(i),static_cast<Eigen::Index>(i),values[i]);
    Matrix result(values.size(),values.size());result.setFromTriplets(entries.begin(),entries.end());return result;
}
double dot(const std::vector<double>& a,const std::vector<double>& b) {
    double result=0;for (std::size_t i=0;i<a.size();++i) result+=a[i]*b[i];return result;
}
std::vector<double> multiply(const Matrix& a,const std::vector<double>& u) {
    Eigen::Map<const Eigen::VectorXd> x(u.data(),u.size());const Eigen::VectorXd y=a*x;return {y.data(),y.data()+y.size()};
}
double relativeError(const std::vector<double>& actual,const std::vector<double>& exact,const std::vector<double>& weights) {
    double error=0,scale=0;
    for (std::size_t i=0;i<exact.size();++i) {error=std::hypot(error,std::sqrt(weights[i])*(actual[i]-exact[i]));scale=std::hypot(scale,std::sqrt(weights[i])*exact[i]);}
    check(std::isfinite(error) && std::isfinite(scale) && scale>0,"invalid independent solution error");return error/scale;
}
void report(const HydroSaddleResult& result) {
    const auto& r=result.report();peakContinuity=std::max(peakContinuity,r.maxCellContinuityDefect);
    peakAlgebraic=std::max(peakAlgebraic,r.algebraicRelativeResidual);peakMean=std::max(peakMean,r.maxPressureMeanDefect);
    check(r.algebraicRelativeResidual<=1e-11 && r.freeMomentumDefect<=1e-11 && r.maxCellContinuityDefect<=1e-11 && r.maxPressureMeanDefect<=1e-11,
        "unaccepted physical/algebraic result returned");
    check((r.usedSuperLU||r.usedBlockGmres)==(r.factorDimension>0),"reported backend does not match actual nonempty quotient");
    check(!r.usedDirectFallback||r.usedSuperLU,"reported fallback did not execute direct backend");
    check(!r.usedBlockGmres||r.krylovSolves>0,"reported Krylov backend did not execute a solve");
    for(double seconds:{r.velocityBlockSeconds,r.factorSeconds,r.solveSeconds})
        check(std::isfinite(seconds)&&seconds>=0,"invalid direct-backend diagnostic timing");
    std::cout<<"H04_BACKEND,"<<(r.usedSuperLU?"SuperLU":"BlockGMRES")<<",gmres,"<<r.usedBlockGmres<<",fallback,"<<r.usedDirectFallback<<",iterations,"<<r.krylovIterations<<",dimension,"<<r.factorDimension<<','<<r.factorNonZeros<<','
        <<r.velocityBlockSeconds<<','<<r.factorSeconds<<','<<r.solveSeconds<<'\n';
}
double polynomialPsi(double x,double y,const QuadtreeMesh& mesh) {
    const double a=(x-mesh.minCoord.x)*(mesh.maxCoord.x-x),b=y*(1-y);return a*a*b*b;
}
void algebra(int n,int metric,int depth,int selection,bool nonDiagonal=false,bool axial=false,bool gmresOnlyParallel=false) {
    QuadtreeMesh mesh(input(n,metric));refine(mesh,depth);const auto graph=selected(mesh,selection);
    const auto layout=HydroDofLayout::build(mesh,NativeLeafDomain::build(mesh,graph,metric!=0));
    const auto load=HydroPressureLoad::build(HydroIntegralDivergence::build(mesh,layout));const auto ref=reference(mesh,graph,layout,metric!=0);
    std::vector<double> u(layout.velocityCount()),p(layout.pressureCount()),known(u.size()),weights(u.size(),1);
    for (std::size_t i=0;i<u.size();++i) {
        const auto& f=ref.faces[i];
        const double flow=f.horizontal ? polynomialPsi(f.x0,f.y0,mesh)-polynomialPsi(f.x1,f.y0,mesh)
                                       : polynomialPsi(f.x0,f.y1,mesh)-polynomialPsi(f.x0,f.y0,mesh);
        u[i]=f.area==0 ? 0 : (metric!=0 ? 2*pi : 1)*flow/f.area;
        if (axial && f.horizontal) u[i]+=1;
        if (f.incidence==1 || f.area==0) known[i]=u[i];
    }
    std::vector<double> mean(selection==1 ? 2 : 1,0),volume(mean.size(),0);
    for (std::size_t i=0;i<p.size();++i) {
        const auto c=ref.center[i];const int component=selection==1 && c.x>mesh.minCoord.x+.5 ? 1 : 0;
        p[i]=1+.3*c.x+.7*c.y+.4*c.x*c.y;mean[component]+=ref.volume[i]*p[i];volume[component]+=ref.volume[i];
    }
    for (std::size_t i=0;i<p.size();++i) p[i]-=mean[selection==1 && ref.center[i].x>mesh.minCoord.x+.5 ? 1 : 0]/volume[selection==1 && ref.center[i].x>mesh.minCoord.x+.5 ? 1 : 0];
    Matrix a=diagonal(std::vector<double>(u.size(),3));
    if (nonDiagonal) for (int i=1;i<layout.velocityCount();++i) {a.coeffRef(i-1,i)=-.125;a.coeffRef(i,i-1)=-.125;}
    HydroSaddleConfig config;config.gmresOnlyParallel=gmresOnlyParallel;
    const auto velocity=HydroSaddleVelocityBlock::build(layout,a);const auto solver=HydroSinglePhaseSaddle::buildClosed(load,velocity,config);
    check(solver.pressureGaugeCount()==static_cast<int>(mean.size()),"wrong component pressure gauge count");
    check(&solver.pressureLoad().divergence().matrix()==&load.divergence().matrix(),"solver lost unique divergence");
    const HydroSinglePhaseSaddle::Matrix asym=solver.matrix()-HydroSinglePhaseSaddle::Matrix(solver.matrix().transpose());
    check(asym.norm()==0,"closed augmented saddle is not exactly symmetric");
    auto momentum=multiply(a,u);for (const auto& term : ref.terms) momentum[term.u]-=term.outwardArea*p[term.p];
    const HydroSaddleRhs rhs(layout,momentum,known);const auto result=solver.solve(mesh,rhs);report(result);
    if(gmresOnlyParallel){const auto& r=result.report();
        check(r.usedBlockGmres&&!r.usedSuperLU&&!r.usedDirectFallback&&r.krylovIterations>0,
            "GMRES-only saddle did not execute a nonzero solve without direct fallback");
        check(r.parallelPressureBlockCount>0&&r.parallelPressureMaximumWorkers>0,
            "GMRES-only saddle did not use its parallel-capable pressure block preconditioner");}
    const double eu=relativeError(result.values().normalVelocity(),u,weights),ep=relativeError(result.values().pressure(),p,ref.volume);
    peakVelocity=std::max(peakVelocity,eu);peakPressure=std::max(peakPressure,ep);check(eu<=1e-9 && ep<=1e-9,"independent algebra recovery failed");
    const auto repeated=solver.solve(mesh,rhs);
    check(repeated.values().normalVelocity()==result.values().normalVelocity() && repeated.values().pressure()==result.values().pressure(),"reused factor changed deterministic RHS");
    if (selection==0 && !axial) {
        const auto actual=result.values().normalVelocity();const double energy=dot(actual,multiply(a,actual)),work=dot(actual,momentum);
        const double defect=std::abs(energy-work)/(std::abs(energy)+std::abs(work));peakEnergy=std::max(peakEnergy,defect);check(defect<=1e-10,"zero-boundary saddle energy relation failed");
    }
    ++algebraCases;
    std::cout << "HYDRO_SADDLE_ALG," << metric << ',' << n << ',' << depth << ',' << selection << ',' << eu << ',' << ep << '\n';
}
double sinc(double x) {return x==0 ? 1 : std::sin(x)/x;}
double analyticPsi(double x,double y) {return x==0 || x==1 || y==0 || y==1 ? 0 : std::sin(pi*x)*std::sin(pi*y);}
std::pair<double,double> mms(int n) {
    QuadtreeMesh mesh(input(n,0));const auto graph=selected(mesh,0);
    const auto layout=HydroDofLayout::build(mesh,NativeLeafDomain::build(mesh,graph,false));const auto ref=reference(mesh,graph,layout,false);
    const auto load=HydroPressureLoad::build(HydroIntegralDivergence::build(mesh,layout));
    constexpr double reaction=2.5;const double h=1.0/n;
    std::vector<double> u(layout.velocityCount()),p(layout.pressureCount()),mass(u.size()),dual(u.size()),momentum(u.size()),known(u.size(),0);
    for (std::size_t i=0;i<u.size();++i) {
        const auto& f=ref.faces[i];const double x=(f.x0+f.x1)/2,y=(f.y0+f.y1)/2;
        u[i]=(f.horizontal ? analyticPsi(f.x0,f.y0)-analyticPsi(f.x1,f.y0) : analyticPsi(f.x0,f.y1)-analyticPsi(f.x0,f.y0))/f.area;
        dual[i]=f.area*h*(f.incidence==1 ? .5 : 1);mass[i]=reaction*dual[i];
        const double gradient=f.horizontal ? -2*pi*std::cos(pi*x)*sinc(pi*h/2)*std::sin(2*pi*y)
                                           : -pi*std::sin(pi*x)*std::cos(2*pi*y)*sinc(pi*h);
        // Exact continuum FACE-AVERAGE forcing a*u + grad(p), projected by
        // the uniform dual-volume Hodge. No numerical C/D feeds this RHS.
        momentum[i]=dual[i]*(reaction*u[i]+gradient);
    }
    for (std::size_t i=0;i<p.size();++i) p[i]=std::cos(pi*ref.center[i].x)*std::cos(2*pi*ref.center[i].y)*sinc(pi*h/2)*sinc(pi*h);
    const auto solver=HydroSinglePhaseSaddle::buildClosed(load,HydroSaddleVelocityBlock::build(layout,diagonal(mass)));
    const auto result=solver.solve(mesh,HydroSaddleRhs(layout,momentum,known));report(result);
    const double eu=relativeError(result.values().normalVelocity(),u,dual),ep=relativeError(result.values().pressure(),p,ref.volume);
    ++mmsCases;std::cout << "HYDRO_SADDLE_MMS," << n << ',' << std::setprecision(16) << eu << ',' << ep << '\n';return {eu,ep};
}
// Independent original augmented solve: geometry incidence, physical RZ
// volumes and boundary elimination are assembled here, without the quotient
// basis, reduced factor or its correction/lifting implementation.
void quotientReferenceWitness() {
    QuadtreeMesh mesh(input(8,2));const auto graph=selected(mesh,1);
    const auto layout=HydroDofLayout::build(mesh,NativeLeafDomain::build(mesh,graph,true));
    const auto ref=reference(mesh,graph,layout,true);
    const auto load=HydroPressureLoad::build(HydroIntegralDivergence::build(mesh,layout));
    const int nu=layout.velocityCount(),np=layout.pressureCount();
    std::vector<double> u(nu),p(np),known(nu),dual(nu),mass(nu);
    std::vector<int> free,byVelocity(nu,-1),component(np);int numbers[2]={-1,-1},ng=0;
    for(int i=0;i<nu;++i){const auto& f=ref.faces[i];
        const double flow=f.horizontal?polynomialPsi(f.x0,f.y0,mesh)-polynomialPsi(f.x1,f.y0,mesh)
            :polynomialPsi(f.x0,f.y1,mesh)-polynomialPsi(f.x0,f.y0,mesh);
        u[i]=2*pi*flow/f.area+(f.horizontal?1:0);
        if(f.incidence==1)known[i]=u[i];
        else {check(f.incidence==2,"quotient witness has invalid face incidence");byVelocity[i]=static_cast<int>(free.size());free.push_back(i);}
        dual[i]=f.area/8*(f.incidence==1?.5:1);
        // Actual unequal density coefficients multiply the independently
        // measured dual volume; no equal-density substitute enters this A.
        mass[i]=dual[i]*((f.x0+f.x1)*.5<.75?1000:1);
    }
    double mean[2]={0,0},volume[2]={0,0};
    for(int i=0;i<np;++i){const auto x=ref.center[i];const int side=x.x<.75?0:1;
        if(numbers[side]<0)numbers[side]=ng++;component[i]=numbers[side];
        p[i]=1+.3*x.x+.7*x.y+.4*x.x*x.y;mean[component[i]]+=ref.volume[i]*p[i];volume[component[i]]+=ref.volume[i];}
    check(ng==2&&*std::min_element(ref.volume.begin(),ref.volume.end())<*std::max_element(ref.volume.begin(),ref.volume.end()),
        "quotient witness lacks disconnected nonuniform physical RZ volumes");
    for(int i=0;i<np;++i)p[i]-=mean[component[i]]/volume[component[i]];
    Matrix spd=diagonal(mass),skew(nu,nu);
    for(int i=1;i<nu;++i){const double scale=std::sqrt(mass[i-1]*mass[i]);
        spd.coeffRef(i-1,i)=spd.coeffRef(i,i-1)=-.03125*scale;
        skew.coeffRef(i-1,i)=.125*scale;skew.coeffRef(i,i-1)=-.125*scale;}
    const auto velocity=HydroSaddleVelocityBlock::build(layout,spd).withSkew(skew);
    const auto solver=HydroSinglePhaseSaddle::buildClosed(load,velocity);const auto& a=velocity.matrix();
    check(solver.freeVelocityDofs()==free&&solver.pressureComponents()==component&&solver.pressureGaugeCount()==2,
        "quotient witness changed independently classified native components");
    check((a-Matrix(a.transpose())).norm()>0,"quotient witness did not exercise SPD plus nonzero skew");
    const Eigen::VectorXd constant=load.matrix()*Eigen::VectorXd::Ones(np);
    for(int face:free)check(constant[face]==0,"actual free C does not annihilate component constants");
    const int nf=static_cast<int>(free.size());using Sparse=HydroSinglePhaseSaddle::Matrix;
    std::vector<Eigen::Triplet<double,Eigen::Index>> entries;
    for(int i=0;i<nf;++i)for(int j=0;j<nf;++j)if(a.coeff(free[i],free[j])!=0)
        entries.emplace_back(i,j,a.coeff(free[i],free[j]));
    for(const auto& term:ref.terms)if(byVelocity[term.u]>=0){const int i=byVelocity[term.u];
        entries.emplace_back(i,nf+term.p,-term.outwardArea);entries.emplace_back(nf+term.p,i,-term.outwardArea);}
    for(int i=0;i<np;++i){entries.emplace_back(nf+i,nf+np+component[i],ref.volume[i]);
        entries.emplace_back(nf+np+component[i],nf+i,ref.volume[i]);}
    Sparse augmented(nf+np+ng,nf+np+ng);augmented.setFromTriplets(entries.begin(),entries.end());augmented.makeCompressed();
    check((augmented-solver.matrix()).norm()<=1e-11*augmented.norm(),"public augmented matrix changed original physical equations");
    Eigen::SparseLU<Sparse> oracle;oracle.compute(augmented);check(oracle.info()==Eigen::Success,"independent augmented SparseLU factor failed");
    auto momentum=multiply(a,u);for(const auto& term:ref.terms)momentum[term.u]-=term.outwardArea*p[term.p];
    double peakU=0,peakP=0,peakGauge=0;int nonzeroBoundary=0;for(double v:known)if(v!=0)++nonzeroBoundary;
    check(nonzeroBoundary>0,"quotient witness has zero prescribed boundary");
    const auto compare=[&](const std::vector<double>& prescribed,bool tinyFlux){
        const auto knownMomentum=multiply(a,prescribed);Eigen::VectorXd b=Eigen::VectorXd::Zero(augmented.rows());
        for(int i=0;i<nf;++i)b[i]=momentum[free[i]]-knownMomentum[free[i]];
        for(const auto& term:ref.terms)b[nf+term.p]+=term.outwardArea*prescribed[term.u];
        const Eigen::VectorXd exact=oracle.solve(b);check(oracle.info()==Eigen::Success&&exact.allFinite(),"independent augmented solve failed");
        check((augmented*exact-b).stableNorm()<=1e-11*b.stableNorm(),"independent original augmented residual exceeds original gate");
        auto exactU=prescribed;std::vector<double> exactP(np);
        for(int i=0;i<nf;++i)exactU[free[i]]=exact[i];for(int i=0;i<np;++i)exactP[i]=exact[nf+i];
        const auto result=solver.solve(mesh,HydroSaddleRhs(layout,momentum,prescribed));report(result);
        peakU=std::max(peakU,relativeError(result.values().normalVelocity(),exactU,dual));
        peakP=std::max(peakP,relativeError(result.values().pressure(),exactP,ref.volume));
        check(peakU<=1e-9&&peakP<=1e-9,"quotient differs from independent original augmented solve");
        check(result.report().factorDimension==nf+np-ng&&result.report().factorNonZeros>0,
            "reported factor is not the actual component pressure quotient");
        double actualMean[2]={0,0},meanWork[2]={0,0};bool first[2]={false,false};
        for(int i=0;i<np;++i){const double v=result.values().pressure()[i];actualMean[component[i]]+=ref.volume[i]*v;
            meanWork[component[i]]+=ref.volume[i]*std::abs(v);
            if(!first[component[i]]){check(v!=0,"output pressure was pinned instead of volume gauged");first[component[i]]=true;}}
        for(int g=0;g<ng;++g){check(meanWork[g]>0,"quotient witness pressure gauge has no positive work");
            peakGauge=std::max(peakGauge,std::abs(actualMean[g])/meanWork[g]);}
        check(peakGauge<=1e-11,"independent physical volume pressure means exceed original gate");
        if(!tinyFlux){check(relativeError(result.values().normalVelocity(),u,dual)<=1e-9
            &&relativeError(result.values().pressure(),p,ref.volume)<=1e-9,"quotient missed manufactured physical state");}
        else check(result.report().maxAbsGaugeMultiplier>0&&exact.tail(ng).cwiseAbs().maxCoeff()>0,
            "tolerance-compatible tiny actual boundary flux did not exercise nonzero gauge multiplier");
        std::cout<<"QUOTIENT_AUGMENTED_ORACLE,RZ_DISCONNECTED,RHO_1000_1,SPD_SKEW,"<<tinyFlux<<','
            <<result.report().factorDimension<<','<<result.report().factorNonZeros<<','<<result.report().maxBoundaryFluxDefect<<','
            <<result.report().maxAbsGaugeMultiplier<<','<<result.report().linearRefinements<<','<<peakU<<','<<peakP<<','<<peakGauge<<"\n";
    };
    compare(known,false);
    auto shifted=known;bool changed=false;
    for(int i=0;i<nu;++i)if(ref.faces[i].incidence==1&&ref.faces[i].horizontal&&ref.faces[i].y0==1){
        // Within the ORIGINAL 1e-12 component-compatibility budget. This
        // changes real prescribed data and keeps every physical gate intact.
        shifted[i]+=std::ldexp(1.,-42);changed=true;break;}
    check(changed,"quotient witness missed compatible top boundary");compare(shifted,true);
    QuadtreeMesh single(input(1,2));const auto singleGraph=selected(single,0);
    const auto singleLayout=HydroDofLayout::build(single,NativeLeafDomain::build(single,singleGraph,true));
    const auto singleSolver=HydroSinglePhaseSaddle::buildClosed(HydroPressureLoad::build(HydroIntegralDivergence::build(single,singleLayout)),
        HydroSaddleVelocityBlock::build(singleLayout,diagonal(std::vector<double>(singleLayout.velocityCount(),1000))));
    const std::vector<double> zero(singleLayout.velocityCount());const auto singleton=singleSolver.solve(single,HydroSaddleRhs(singleLayout,zero,zero));report(singleton);
    check(singleton.report().factorDimension==0&&singleton.report().factorNonZeros==0
        &&singleton.values().normalVelocity()==zero&&singleton.values().pressure()==std::vector<double>(1,0),
        "singleton zero-dimensional quotient did not retain exact zero RHS");
    std::cout<<"QUOTIENT_SINGLETON,ZERO_DIMENSION,EXACT_ZERO,PASS\n";
}
// A single mixed RZ system: density contrast occurs INSIDE the flowing
// component, one pressure outlet anchors it, and the disconnected component
// retains its true volume gauge. The independent oracle uses dense LU of the
// original physical equations, without the production pressure quotient.
void pressureOutletWitness() {
    QuadtreeMesh mesh(input(8,1));const auto graph=selected(mesh,1);
    const auto layout=HydroDofLayout::build(mesh,NativeLeafDomain::build(mesh,graph,true));
    const auto ref=reference(mesh,graph,layout,true);
    const auto load=HydroPressureLoad::build(HydroIntegralDivergence::build(mesh,layout));
    const int nu=layout.velocityCount(),np=layout.pressureCount();constexpr double h=.125,pout=2.5;
    const auto closeWork=[](long double actual,long double expected,long double work,const char* why,double tolerance=1e-12){
        const auto gap=std::abs(actual-expected),relative=work==0?(gap==0?0:std::numeric_limits<long double>::infinity()):gap/work;
        check(std::isfinite(relative)&&relative<=tolerance,why);return relative;};
    std::vector<double> u(nu),p(np),known(nu),dual(nu),mass(nu),traction(nu);
    std::vector<int> free,byVelocity(nu,-1),component(np),axis;
    std::vector<HydroSaddlePressureOutlet> outlets;int numbers[2]={-1,-1},components=0;
    bool flowingLiquid=false,flowingGas=false;
    for(int i=0;i<nu;++i){const auto& f=ref.faces[i];const bool open=(f.x0+f.x1)*.5<.5;
        const double rlo=f.horizontal?f.x0:std::max(open?0.:.75,f.x0-h/2);
        const double rhi=f.horizontal?f.x1:std::min(open?.25:1.,f.x0+h/2);
        const double zlo=f.horizontal?std::max(0.,f.y0-h/2):f.y0;
        const double zhi=f.horizontal?std::min(1.,f.y0+h/2):f.y1;
        const double crossSection=pi*(rhi*rhi-rlo*rlo),liquid=crossSection*std::max(0.,std::min(zhi,.5)-zlo);
        dual[i]=crossSection*(zhi-zlo);check(dual[i]>0,"pressure outlet lost the actual positive RZ half-dual volume");
        mass[i]=1000*liquid+(dual[i]-liquid);
        if(open&&f.horizontal){u[i]=.2;flowingLiquid=flowingLiquid||liquid==dual[i];flowingGas=flowingGas||liquid==0;}
        const bool outlet=open&&f.horizontal&&f.incidence==1&&f.y0==mesh.maxCoord.y;
        if(outlet){outlets.push_back({i,pout});traction[i]=-f.area*pout;}
        if(f.incidence==2||outlet){byVelocity[i]=static_cast<int>(free.size());free.push_back(i);}
        else known[i]=u[i];
        if(f.area==0){check(!f.horizontal&&f.x0==0&&u[i]==0,"pressure outlet fixture has a false radial axis");axis.push_back(i);}
    }
    check(flowingLiquid&&flowingGas&&!axis.empty()&&!outlets.empty(),"pressure outlet lacks true flowing 1000:1 contrast or axis");
    double closedMean=0,closedVolume=0;
    for(int i=0;i<np;++i){const int side=ref.center[i].x<.5?0:1;
        if(numbers[side]<0)numbers[side]=components++;component[i]=numbers[side];
        p[i]=side==0?pout:1+.3*ref.center[i].x+.7*ref.center[i].y+.4*ref.center[i].x*ref.center[i].y;
        if(side==1){closedMean+=ref.volume[i]*p[i];closedVolume+=ref.volume[i];}}
    check(components==2&&closedVolume>0,"pressure outlet fixture lost disconnected components");
    for(int i=0;i<np;++i)if(component[i]==numbers[1])p[i]-=closedMean/closedVolume;
    Matrix spd=diagonal(mass),skew(nu,nu);
    for(int i=1;i<nu;++i){const double scale=std::sqrt(mass[i-1]*mass[i]);
        spd.coeffRef(i-1,i)=spd.coeffRef(i,i-1)=-.03125*scale;
        skew.coeffRef(i-1,i)=.125*scale;skew.coeffRef(i,i-1)=-.125*scale;}
    const auto velocity=HydroSaddleVelocityBlock::build(layout,spd).withSkew(skew);const auto& a=velocity.matrix();
    const auto solver=HydroSinglePhaseSaddle::buildWithPressureOutlets(load,velocity,outlets);
    check(solver.freeVelocityDofs()==free&&solver.pressureComponents()==component&&solver.pressureGaugeCount()==1,
        "pressure outlet changed the actual free faces or mixed component gauges");
    const auto& gauge=solver.pressureGaugeByComponent();check(gauge.size()==2&&gauge[numbers[0]]==-1&&gauge[numbers[1]]==0,
        "pressure outlet either volume gauged anchored pressure or dropped the closed gauge");
    check((a-Matrix(a.transpose())).norm()>0&&solver.pressureOutlets().size()==outlets.size()
        &&solver.pressureOutletLoads().size()==static_cast<std::size_t>(nu),"pressure outlet lost SPD/skew or actual traction inventory");
    for(std::size_t k=0;k<outlets.size();++k)check(solver.pressureOutlets()[k].velocityDof==outlets[k].velocityDof
        &&solver.pressureOutlets()[k].pressure==pout,"pressure outlet replaced its actual prescribed pressure");
    for(int i=0;i<nu;++i)closeWork(solver.pressureOutletLoads()[i],traction[i],std::abs(traction[i]),"pressure outlet wrong signed physical traction");
    const int nf=static_cast<int>(free.size()),size=nf+np+1;
    Eigen::MatrixXd dense=Eigen::MatrixXd::Zero(size,size);
    for(int i=0;i<nf;++i)for(int j=0;j<nf;++j)dense(i,j)=a.coeff(free[i],free[j]);
    for(const auto& t:ref.terms)if(byVelocity[t.u]>=0){const int i=byVelocity[t.u];
        dense(i,nf+t.p)-=t.outwardArea;dense(nf+t.p,i)-=t.outwardArea;}
    for(int i=0;i<np;++i)if(component[i]==numbers[1])dense(nf+i,nf+np)=dense(nf+np,nf+i)=ref.volume[i];
    const Eigen::MatrixXd publicMatrix(solver.matrix());
    closeWork((dense-publicMatrix).norm(),0,dense.norm(),"pressure outlet public mixed K differs from original physical equations",1e-11);
    Eigen::FullPivLU<Eigen::MatrixXd> oracle(dense);check(oracle.isInvertible(),"independent pressure outlet dense mixed K is singular");
    // Physical bulk RHS is formed independently; T is inserted separately.
    auto momentum=multiply(a,u);
    for(const auto& t:ref.terms)if(byVelocity[t.u]>=0)momentum[t.u]-=t.outwardArea*p[t.p];
    // On imposed normals the supplied bulk contains Au only: the actual
    // unprescribed pressure force produces a nonzero inlet reaction work.
    for(int i=0;i<nu;++i)momentum[i]-=traction[i];
    const auto eliminated=multiply(a,known);Eigen::VectorXd rhs=Eigen::VectorXd::Zero(size);
    for(int i=0;i<nf;++i)rhs[i]=momentum[free[i]]+traction[free[i]]-eliminated[free[i]];
    for(const auto& t:ref.terms)rhs[nf+t.p]+=t.outwardArea*known[t.u];
    const Eigen::VectorXd exact=oracle.solve(rhs);check(exact.allFinite(),"independent pressure outlet dense LU produced nonfinite values");
    closeWork((dense*exact-rhs).stableNorm(),0,rhs.stableNorm(),"independent pressure outlet dense LU residual",1e-11);
    auto exactU=known;std::vector<double> exactP(np);
    for(int i=0;i<nf;++i)exactU[free[i]]=exact[i];for(int i=0;i<np;++i)exactP[i]=exact[nf+i];
    const auto result=solver.solve(mesh,HydroSaddleRhs(layout,momentum,known));report(result);
    const auto& actual=result.values().normalVelocity();const auto& pressure=result.values().pressure();
    const auto eu=relativeError(actual,exactU,dual),ep=relativeError(pressure,exactP,ref.volume);
    check(eu<=1e-9&&ep<=1e-9&&relativeError(actual,u,dual)<=1e-9&&relativeError(pressure,p,ref.volume)<=1e-9,
        "pressure outlet differs from independent dense or manufactured physical state");
    check(result.report().pressureGaugeCount==1&&result.report().factorDimension==nf+np-1
        &&result.report().maxBoundaryFluxDefect<=1e-12,"pressure outlet altered the original quotient/closed flux gate");
    const auto repeated=solver.solve(mesh,HydroSaddleRhs(layout,momentum,known));
    check(repeated.values().normalVelocity()==actual&&repeated.values().pressure()==pressure,"pressure outlet reuse changed deterministic solution");
    std::vector<long double> au(nu),aWork(nu),cp(nu),cpWork(nu),div(np),divWork(np);
    for(int outer=0;outer<a.outerSize();++outer)for(Matrix::InnerIterator it(a,outer);it;++it){const auto term=static_cast<long double>(it.value())*actual[it.col()];
        au[it.row()]+=term;aWork[it.row()]+=std::abs(term);}
    for(const auto& t:ref.terms){const auto force=-static_cast<long double>(t.outwardArea)*pressure[t.p];cp[t.u]+=force;cpWork[t.u]+=std::abs(force);
        const auto flux=static_cast<long double>(t.outwardArea)*actual[t.u];div[t.p]+=flux;divWork[t.p]+=std::abs(flux);}
    long double freeGap=0,freeWork=0,energy=0,bulkWork=0,tractionWork=0,reactionWork=0,pressureWork=0,totalWork=0;
    long double netOpen=0,grossOpen=0,inflow=0,outflow=0,openMean=0,openMeanWork=0,meanClosed=0,meanClosedWork=0;
    for(int i=0;i<nu;++i){const auto residual=au[i]+cp[i]-momentum[i]-traction[i];
        if(byVelocity[i]>=0){freeGap=std::hypot(freeGap,residual);freeWork=std::hypot(freeWork,aWork[i]+cpWork[i]+std::abs(momentum[i])+std::abs(traction[i]));}
        else {check(actual[i]==known[i],"pressure outlet solve changed an imposed normal or axis");reactionWork+=actual[i]*residual;}
        energy+=actual[i]*au[i];bulkWork+=actual[i]*momentum[i];tractionWork+=actual[i]*traction[i];pressureWork+=actual[i]*cp[i];
        totalWork+=std::abs(actual[i])*(aWork[i]+cpWork[i]+std::abs(momentum[i])+std::abs(traction[i]));}
    for(int i=0;i<np;++i){closeWork(div[i],0,divWork[i],"pressure outlet actual native continuity",1e-11);
        const auto term=static_cast<long double>(ref.volume[i])*pressure[i];
        if(component[i]==numbers[0]){openMean+=term;openMeanWork+=std::abs(term);}else{meanClosed+=term;meanClosedWork+=std::abs(term);}}
    for(const auto& t:ref.terms)if(ref.faces[t.u].incidence==1&&component[t.p]==numbers[0]){
        const auto flux=static_cast<long double>(t.outwardArea)*actual[t.u];netOpen+=flux;grossOpen+=std::abs(flux);
        if(std::find_if(outlets.begin(),outlets.end(),[&](const auto& o){return o.velocityDof==t.u;})!=outlets.end())outflow+=flux;
        else if(ref.faces[t.u].horizontal&&ref.faces[t.u].y0==0)inflow-=flux;}
    check(inflow>0&&outflow>0&&openMean!=0&&openMeanWork>0&&tractionWork<0&&reactionWork!=0,
        "pressure outlet lost actual supply/free outflow, absolute pressure, traction or fixed reaction work");
    closeWork(netOpen,0,grossOpen,"pressure outlet actual free outflow balances prescribed supply",1e-11);
    closeWork(outflow,inflow,std::abs(inflow)+std::abs(outflow),"pressure outlet free physical flux balance",1e-11);
    closeWork(meanClosed,0,meanClosedWork,"pressure outlet untouched true-volume closed pressure gauge",1e-11);
    closeWork(freeGap,0,freeWork,"pressure outlet original free physical momentum",1e-11);
    closeWork(tractionWork,-pout*outflow,std::abs(tractionWork)+std::abs(pout*outflow),"pressure outlet independent nonzero physical traction power");
    closeWork(energy+pressureWork-bulkWork-tractionWork-reactionWork,0,totalWork+std::abs(reactionWork),
        "pressure outlet physical power with actual fixed-normal reaction",1e-11);
    // SAME operator, zero supply and bulk: constant absolute open pressure
    // balances T, while the disconnected closed pressure remains volume zero.
    const std::vector<double> zero(nu);const auto staticResult=solver.solve(mesh,HydroSaddleRhs(layout,zero,zero));report(staticResult);
    std::vector<double> staticPressure(np);for(int i=0;i<np;++i)if(component[i]==numbers[0])staticPressure[i]=pout;
    check(relativeError(staticResult.values().pressure(),staticPressure,ref.volume)<=1e-9,"static pressure outlet lost Cp=T absolute pressure sign");
    const auto staticAu=multiply(a,staticResult.values().normalVelocity());long double staticGap=0,staticWork=0;
    for(int i:free){long double force=0,work=0;for(const auto& t:ref.terms)if(t.u==i){const auto term=-static_cast<long double>(t.outwardArea)*staticResult.values().pressure()[t.p];force+=term;work+=std::abs(term);}
        staticGap=std::hypot(staticGap,static_cast<long double>(staticAu[i])+force-traction[i]);staticWork=std::hypot(staticWork,std::abs(staticAu[i])+work+std::abs(traction[i]));}
    closeWork(staticGap,0,staticWork,"static pressure outlet actual Cp=T free equation",1e-11);
    int guards=0;const auto refuse=[&](const std::function<void()>& action){bool caught=false;try{action();}catch(const std::exception&){caught=true;}
        check(caught,"invalid pressure outlet data accepted");++guards;};
    auto bad=outlets;bad.push_back({axis.front(),pout});refuse([&]{HydroSinglePhaseSaddle::buildWithPressureOutlets(load,velocity,bad);});
    bad=outlets;int internal=-1;for(int i=0;i<nu;++i)if(ref.faces[i].incidence==2){internal=i;break;}
    check(internal>=0,"pressure outlet guard has no actual internal face");bad.push_back({internal,pout});refuse([&]{HydroSinglePhaseSaddle::buildWithPressureOutlets(load,velocity,bad);});
    bad=outlets;bad.push_back(outlets.front());refuse([&]{HydroSinglePhaseSaddle::buildWithPressureOutlets(load,velocity,bad);});
    bad=outlets;bad.front().pressure=std::numeric_limits<double>::quiet_NaN();refuse([&]{HydroSinglePhaseSaddle::buildWithPressureOutlets(load,velocity,bad);});
    auto wrongKnown=known;wrongKnown[outlets.front().velocityDof]=.2;refuse([&]{solver.solve(mesh,HydroSaddleRhs(layout,momentum,wrongKnown));});
    const auto closed=HydroSinglePhaseSaddle::buildClosed(load,velocity);refuse([&]{closed.solve(mesh,HydroSaddleRhs(layout,momentum,known));});
    wrongKnown=known;bool changed=false;for(const auto& t:ref.terms)if(component[t.p]==numbers[1]&&ref.faces[t.u].incidence==1
        &&ref.faces[t.u].horizontal&&ref.faces[t.u].y0==0){wrongKnown[t.u]=.2;changed=true;break;}
    check(changed,"pressure outlet guard missed the disconnected closed boundary");refuse([&]{solver.solve(mesh,HydroSaddleRhs(layout,momentum,wrongKnown));});
    const auto after=solver.solve(mesh,HydroSaddleRhs(layout,momentum,known));check(after.values().normalVelocity()==actual&&after.values().pressure()==pressure,
        "failed pressure outlet solve poisoned retained system/result");check(guards==7,"pressure outlet guard coverage incomplete");
    std::cout<<"PRESSURE_OUTLET_ONLY,RZ_MIXED,RHO_1000_1_INTERNAL,SPD_SKEW,DENSE_LU,"<<eu<<','<<ep<<','<<inflow<<','<<outflow
        <<','<<tractionWork<<','<<reactionWork<<",rejected="<<guards<<",PASS\n";
}

void guardsAndLifetime() {
    static_assert(!std::is_default_constructible_v<HydroSinglePhaseSaddle>);
    static_assert(std::is_same_v<decltype(std::declval<const HydroSaddleResult&>().values()),const HydroNativeValues&>);
    QuadtreeMesh mesh(input(4,0));const auto graph=selected(mesh,0);const auto layout=HydroDofLayout::build(mesh,NativeLeafDomain::build(mesh,graph,false));
    Matrix a=diagonal(std::vector<double>(layout.velocityCount(),1));const auto load=HydroPressureLoad::build(HydroIntegralDivergence::build(mesh,layout));
    auto withZeros=a;withZeros.coeffRef(0,1)=0.;withZeros.coeffRef(1,0)=-0.;withZeros.makeCompressed();
    check(withZeros.nonZeros()==layout.velocityCount()+2,"diagonal guard lacks stored off-diagonal zeros");
    const auto retainedZeros=HydroSaddleVelocityBlock::build(layout,withZeros);
    check(retainedZeros.matrix().nonZeros()==withZeros.nonZeros()&&retainedZeros.matrix().coeff(0,1)==0&&retainedZeros.matrix().coeff(1,0)==0
        &&std::signbit(retainedZeros.matrix().coeff(1,0)),"diagonal SPD certification changed explicit zero entries");
    const auto velocity=HydroSaddleVelocityBlock::build(layout,a);const auto solver=HydroSinglePhaseSaddle::buildClosed(load,velocity);const auto copy=solver;
    check(&copy.matrix()==&solver.matrix(),"saddle copy is not shared read-only");a.coeffRef(0,0)=7;check(velocity.matrix().coeff(0,0)==1,"velocity block aliases caller matrix");
    std::vector<double> zero(layout.velocityCount(),0);const HydroSaddleRhs rhs(layout,zero,zero);
    const auto accepted=solver.solve(mesh,rhs);check(accepted.values().normalVelocity()==zero && accepted.values().pressure()==std::vector<double>(layout.pressureCount(),0),"zero RHS is not exact zero");
    auto bad=a;bad.coeffRef(0,0)=-1;rejects([&] {HydroSaddleVelocityBlock::build(layout,bad);});
    bad=a;bad.coeffRef(0,0)=0;rejects([&] {HydroSaddleVelocityBlock::build(layout,bad);});
    bad.prune(0.);rejects([&] {HydroSaddleVelocityBlock::build(layout,bad);},"not positive definite");
    bad=a;bad.coeffRef(0,1)=.1;rejects([&] {HydroSaddleVelocityBlock::build(layout,bad);});
    bad=a;bad.coeffRef(0,0)=std::numeric_limits<double>::infinity();rejects([&] {HydroSaddleVelocityBlock::build(layout,bad);});
    rejects([&] {HydroSaddleVelocityBlock::build(layout,Matrix(2,2));});
    auto config=HydroSaddleConfig{};config.algebraicRelativeTolerance=0;rejects([&] {HydroSinglePhaseSaddle::buildClosed(load,velocity,config);});
    config=HydroSaddleConfig{};config.continuityRelativeTolerance=std::numeric_limits<double>::quiet_NaN();rejects([&] {HydroSinglePhaseSaddle::buildClosed(load,velocity,config);});
    const auto fresh=HydroDofLayout::build(mesh,NativeLeafDomain::build(mesh,graph,false));
    rejects([&] {HydroSinglePhaseSaddle::buildClosed(load,HydroSaddleVelocityBlock::build(fresh,velocity.matrix()));});
    rejects([&] {solver.solve(mesh,HydroSaddleRhs(fresh,zero,zero));});
    rejects([&] {HydroSaddleRhs(layout,std::vector<double>(2),zero);});
    auto payload=zero;payload[0]=std::numeric_limits<double>::quiet_NaN();rejects([&] {HydroSaddleRhs(layout,payload,zero);});
    payload=zero;payload[solver.freeVelocityDofs().front()]=1;rejects([&] {solver.solve(mesh,HydroSaddleRhs(layout,zero,payload));});
    const auto ref=reference(mesh,graph,layout,false);payload=zero;
    for (std::size_t i=0;i<payload.size();++i) if (ref.faces[i].incidence==1) {payload[i]=1;break;}
    rejects([&] {solver.solve(mesh,HydroSaddleRhs(layout,zero,payload));});
    auto stretchedInput=input(4,0);stretchedInput.maxCoord.x=1.25;QuadtreeMesh stretched(stretchedInput);rejects([&] {solver.solve(stretched,rhs);});
    auto& vertex=mesh.verts(mesh.vertIndexToLabel[0]);const auto saved=vertex.coord();vertex.setCoord(saved+Vector2D{.01,0});rejects([&] {solver.solve(mesh,rhs);});vertex.setCoord(saved);
    payload=zero;payload[solver.freeVelocityDofs().front()]=std::numeric_limits<double>::max();rejects([&] {solver.solve(mesh,HydroSaddleRhs(layout,payload,zero));});
    payload=zero;payload[solver.freeVelocityDofs().front()]=1;
    config=HydroSaddleConfig{};config.algebraicRelativeTolerance=1e-30;
    const auto strictAlgebra=HydroSinglePhaseSaddle::buildClosed(load,velocity,config);
    rejects([&] {strictAlgebra.solve(mesh,HydroSaddleRhs(layout,payload,zero));},"complete algebraic");
    config=HydroSaddleConfig{};config.continuityRelativeTolerance=1e-30;
    const auto strictContinuity=HydroSinglePhaseSaddle::buildClosed(load,velocity,config);
    rejects([&] {strictContinuity.solve(mesh,HydroSaddleRhs(layout,payload,zero));},"native cell continuity");
    check(solver.solve(mesh,rhs).values().normalVelocity()==accepted.values().normalVelocity(),"failed solve poisoned retained system/result");
    refine(mesh,1);rejects([&] {solver.solve(mesh,rhs);});
    QuadtreeMesh axis(input(4,1));const auto axisLayout=HydroDofLayout::build(axis,NativeLeafDomain::build(axis,selected(axis,0),true));
    const auto axisLoad=HydroPressureLoad::build(HydroIntegralDivergence::build(axis,axisLayout));
    const auto axisSolver=HydroSinglePhaseSaddle::buildClosed(axisLoad,HydroSaddleVelocityBlock::build(axisLayout,diagonal(std::vector<double>(axisLayout.velocityCount(),1))));
    std::vector<double> axisKnown(axisLayout.velocityCount(),0);
    for (Eigen::Index u=0;u<axisLoad.matrix().outerSize();++u) if (!HydroPressureLoad::Matrix::InnerIterator(axisLoad.matrix(),u)) {axisKnown[u]=1;break;}
    rejects([&] {axisSolver.solve(axis,HydroSaddleRhs(axisLayout,std::vector<double>(axisKnown.size()),axisKnown));});
    QuadtreeMesh disconnected(input(4,0));const auto parts=selected(disconnected,1);
    const auto partLayout=HydroDofLayout::build(disconnected,NativeLeafDomain::build(disconnected,parts,false));
    const auto partLoad=HydroPressureLoad::build(HydroIntegralDivergence::build(disconnected,partLayout));
    const auto partSolver=HydroSinglePhaseSaddle::buildClosed(partLoad,HydroSaddleVelocityBlock::build(partLayout,diagonal(std::vector<double>(partLayout.velocityCount(),1))));
    std::vector<double> incompatible(partLayout.velocityCount(),0);std::vector<bool> marked(2,false);
    for (int i=0;i<partLayout.velocityCount();++i) {
        HydroPressureLoad::Matrix::InnerIterator it(partLoad.matrix(),i);if (!it) continue;
        const int p=static_cast<int>(it.col());const double coefficient=it.value();++it;
        if (!it && !marked[partSolver.pressureComponents()[p]]) {
            const int component=partSolver.pressureComponents()[p];incompatible[i]=(component==0 ? 1 : -1)/(-coefficient);marked[component]=true;
        }
    }
    check(marked[0] && marked[1],"fixture missed both component boundaries");
    const auto globalFlow=partLoad.divergence().apply(disconnected,HydroNativeValues(partLayout,incompatible,std::vector<double>(partLayout.pressureCount())));
    double total=0;for (double flow : globalFlow) total+=flow;check(total==0,"component incompatibility fixture is not globally balanced");
    rejects([&] {partSolver.solve(disconnected,HydroSaddleRhs(partLayout,std::vector<double>(incompatible.size()),incompatible));},"component boundary flux compatibility");
    std::optional<HydroSinglePhaseSaddle> retained;std::optional<HydroSaddleRhs> retainedRhs;
    {QuadtreeMesh temporary(input(4,0));const auto l=HydroDofLayout::build(temporary,NativeLeafDomain::build(temporary,selected(temporary,0),false));
     retained=HydroSinglePhaseSaddle::buildClosed(HydroPressureLoad::build(HydroIntegralDivergence::build(temporary,l)),HydroSaddleVelocityBlock::build(l,diagonal(std::vector<double>(l.velocityCount(),1))));
     retainedRhs.emplace(l,std::vector<double>(l.velocityCount()),std::vector<double>(l.velocityCount()));}
    QuadtreeMesh equivalent(input(4,0));check(retained->solve(equivalent,*retainedRhs).values().normalVelocity()==zero,"saddle borrowed source lifetime");
}

namespace SparseValidationPlanContract {
std::size_t cases=0,coefficientChecks=0,negativeCases=0,flagComparisons=0,actualSolveChecks=0,appendChecks=0;
std::uint64_t bits(double x){std::uint64_t r=0;std::memcpy(&r,&x,sizeof r);return r;}
struct Floating {
    std::fenv_t native{};int error=errno;
#if defined(_MSC_VER) && defined(_M_X64)
    unsigned csr=_mm_getcsr();
#endif
    Floating(){check(std::fegetenv(&native)==0,"sparse validation cannot capture native FP");}
    void restore()const {check(std::fesetenv(&native)==0,"sparse validation cannot restore native FP");
#if defined(_MSC_VER) && defined(_M_X64)
        _mm_setcsr(csr);
#endif
        errno=error;}
};
struct RestoreFloating {Floating saved;~RestoreFloating()noexcept{std::fesetenv(&saved.native);
#if defined(_MSC_VER) && defined(_M_X64)
    _mm_setcsr(saved.csr);
#endif
    errno=saved.error;}};
void sameFlags(const Floating& a,const Floating& b){
    check(std::memcmp(&a.native,&b.native,sizeof a.native)==0&&a.error==b.error,"sparse validation changed original caller FP flags/controls/errno");
#if defined(_MSC_VER) && defined(_M_X64)
    check(a.csr==b.csr,"sparse validation changed original raw MXCSR flags/controls");
#endif
    ++flagComparisons;
}
void sameMatrix(const Matrix& a,const Matrix& b) {
    check(a.rows()==b.rows()&&a.cols()==b.cols()&&a.nonZeros()==b.nonZeros()&&a.isCompressed()&&b.isCompressed(),"sparse validation changed matrix shape/storage");
    for(Eigen::Index i=0;i<=a.outerSize();++i)check(a.outerIndexPtr()[i]==b.outerIndexPtr()[i],"sparse validation changed exact outer pattern");
    for(Eigen::Index i=0;i<a.nonZeros();++i){check(a.innerIndexPtr()[i]==b.innerIndexPtr()[i]&&bits(a.valuePtr()[i])==bits(b.valuePtr()[i]),
        "sparse validation changed original coefficient/index bits, including signed zero");++coefficientChecks;}
}
Matrix specimen(int n,int firstPair=1) {
    Matrix a=diagonal(std::vector<double>(n,4.));a.coeffRef(0,firstPair)=.375;a.coeffRef(firstPair,0)=.125;
    a.coeffRef(2,3)=-.125;a.coeffRef(3,2)=.125;a.coeffRef(4,5)=0.;a.coeffRef(5,4)=-0.;a.makeCompressed();return a;
}
HydroSaddleVelocityBlock compare(const QuadtreeMesh& mesh,const HydroDofLayout& layout,const HydroMomentumMass& mass,
    const HydroEnergyViscosity& energy,const HydroMomentumGeometry& owner,const Matrix& a,bool expectUsed,bool expectHit,const char* scene) {
    RestoreFloating restore;check(std::feraiseexcept(FE_INEXACT)==0,"cannot prime sticky inexact for passive timer comparison");
    const Matrix original=a;
    // Complete original full-norm expression: no sampled coefficients or
    // tolerance-based proof substitutes for the old public diagnostic.
    const Matrix at(a.transpose()),symmetric=.5*(a+at),skew=.5*(a-at),combined=symmetric+skew;
    const Matrix transpose(combined.transpose()),difference=combined-transpose;
    const double scale=combined.norm(),expected=difference.norm()/scale;
    check(std::isfinite(scale)&&scale>0&&std::isfinite(expected),"independent sparse full norm fixture is unresolved");
    const Floating before;
    const auto old=HydroSaddleVelocityBlock::buildWithEnergyTransportSplit(mesh,layout,a,mass,energy,.125);const Floating oldAfter;
    before.restore();
    const auto prepared=HydroSaddleVelocityBlock::buildWithEnergyTransportSplit(mesh,layout,a,mass,energy,.125,owner);const Floating newAfter;
    sameFlags(oldAfter,newAfter);sameMatrix(prepared.matrix(),old.matrix());sameMatrix(prepared.matrix(),combined);sameMatrix(a,original);
    check(bits(prepared.symmetryDefect())==bits(expected)&&bits(old.symmetryDefect())==bits(expected),"sparse plan skipped/changed original final complete norm diagnostic");
    const auto& r=prepared.sparseValidationReport();check(r.requested,"typed sparse owner request was not reported");
    check(r.eligible==expectUsed&&r.planHit==expectHit&&r.exactSymmetryUsed==expectUsed&&r.exactSkewUsed==expectUsed,
        "sparse plan report claims a shortcut/hit inconsistent with actual qualified fixture");
    if(expectUsed)check(r.eligible&&r.nonZeros>0&&r.retainedBytes>0&&r.retainedBytes<=64u*1024u*1024u,
        "sparse exact plan does not report actual bounded retained structure");
    else check(!r.exactSymmetryUsed&&!r.exactSkewUsed,"unsupported sparse validation bypassed original full expressions");
    const auto& noRequest=old.sparseValidationReport();check(!noRequest.requested&&!noRequest.planHit&&!noRequest.exactSymmetryUsed&&!noRequest.exactSkewUsed,
        "old six-argument factory silently entered new sparse plan");
    bool appendExpected=expectUsed;
#if EIGEN_WORLD_VERSION!=3 || EIGEN_MAJOR_VERSION!=4 || EIGEN_MINOR_VERSION!=1 || defined(EIGEN_SCALAR_BINARY_OP_PLUGIN)
    appendExpected=false;
#endif
    check(r.appendPlanUsed==appendExpected,"Task39 actual append use differs from supported complete pattern/FP");
    if(appendExpected){
        check(r.appendScratchBytes==static_cast<std::size_t>(prepared.matrix().nonZeros())*sizeof(double)
            &&r.appendScratchBytes>0&&r.appendScratchBytes<=64u*1024u*1024u,
            "Task39 actual fresh difference scratch count/cap");++appendChecks;
    }else check(r.appendScratchBytes==0,"Task39 literal fallback fabricated used scratch");
    check(!noRequest.appendPlanUsed&&noRequest.appendScratchBytes==0,"Task39 unrequested literal silently routed append");
    ++cases;std::cout<<"SPARSE_VALIDATION_PLAN,"<<scene<<",dofs,"<<layout.velocityCount()<<",nnz,"<<a.nonZeros()
        <<",used,"<<r.exactSymmetryUsed<<",hit,"<<r.planHit<<",full_norm,"<<expected<<",matrix_bits_and_FP,original\n";
    return prepared;
}
template<class F>void invalid(F action){bool caught=false;try{action();}catch(const std::invalid_argument&){caught=true;}
    check(caught,"typed sparse validation bypassed an original input/owner guard");++negativeCases;}
void scene(bool rz) {
    QuadtreeMesh mesh(input(4,rz?1:0));refine(mesh,1);FaceSegmentGraph graph;graph.build(mesh);
    const auto layout=HydroDofLayout::build(mesh,NativeLeafDomain::build(mesh,graph,rz));
    check(!layout.domain().uniformLeaves(),"sparse validation typed owner fixture lost actual AMR");
    const auto metrics=HydroPhysicalDivergence::build(mesh,HydroIntegralDivergence::build(mesh,layout));
    const auto geometry=HydroMomentumGeometry::build(mesh,metrics),alias=geometry;
    const LAMP::TimeLayer clock{0,0,0,LAMP::TimeLayer::Position::Initial};
    const auto source=LAMP::Interface::GeometricVofPublication::prepareLocalPlic(mesh,layout.domain(),clock,std::vector<double>(mesh.cellIndexNum(),1.));
    const auto mass=HydroMomentumMass::buildCellAverageProjection(mesh,geometry,HydroPhaseGeometry::build(mesh,metrics,source),1.,1.);
    const auto energy=HydroEnergyViscosity::zero(mesh,metrics);const auto a=specimen(layout.velocityCount());
#if defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
    constexpr bool nativePlan=true;
#else
    constexpr bool nativePlan=false;
#endif
    const auto initial=compare(mesh,layout,mass,energy,alias,a,nativePlan,false,rz?"RZ first":"XY first");
    auto changed=a;changed.coeffRef(0,1)=.625;changed.coeffRef(1,0)=.375;changed.makeCompressed();
    const auto reused=compare(mesh,layout,mass,energy,geometry,changed,nativePlan,nativePlan,"same-pattern fresh values");
    check(initial.matrix().coeff(0,1)!=reused.matrix().coeff(0,1),"sparse index plan froze an earlier numeric coefficient");
    auto zeros=a;for(Eigen::Index row=0;row<zeros.outerSize();++row)for(Matrix::InnerIterator e(zeros,row);e;++e)
        if(e.row()!=e.col())e.valueRef()=e.row()<e.col()?0.:-0.;
    check(zeros.nonZeros()==a.nonZeros()&&std::signbit(zeros.coeff(5,4)),"sparse zero fixture discarded stored negative zero");
    compare(mesh,layout,mass,energy,geometry,zeros,nativePlan,nativePlan,"stored signed zero and zero skew");
    auto large=a;for(Eigen::Index i=0;i<large.nonZeros();++i)large.valuePtr()[i]=std::ldexp(large.valuePtr()[i],480);
    compare(mesh,layout,mass,energy,geometry,large,nativePlan,nativePlan,"large finite normal");
    auto small=a;for(Eigen::Index i=0;i<small.nonZeros();++i)small.valuePtr()[i]=std::ldexp(small.valuePtr()[i],-480);
    compare(mesh,layout,mass,energy,geometry,small,nativePlan,nativePlan,"small finite normal");
    const auto otherPattern=specimen(layout.velocityCount(),6);check(otherPattern.nonZeros()==a.nonZeros(),"foreign-pattern witness changed stored count");
    const Matrix sourceTranspose(a.transpose()),targetTranspose(otherPattern.transpose());
    const Matrix sourceSymmetric=.5*(a+sourceTranspose),targetSymmetric=.5*(otherPattern+targetTranspose);
    check(sourceSymmetric.nonZeros()==targetSymmetric.nonZeros()&&
        (!std::equal(sourceSymmetric.outerIndexPtr(),sourceSymmetric.outerIndexPtr()+sourceSymmetric.outerSize()+1,targetSymmetric.outerIndexPtr())
        ||!std::equal(sourceSymmetric.innerIndexPtr(),sourceSymmetric.innerIndexPtr()+sourceSymmetric.nonZeros(),targetSymmetric.innerIndexPtr())),
        "foreign-pattern witness did not change actual S union pattern at the same stored count");
    compare(mesh,layout,mass,energy,geometry,otherPattern,false,false,"same-count different exact pattern");
    auto subnormal=a;subnormal.coeffRef(0,1)=std::numeric_limits<double>::denorm_min();subnormal.coeffRef(1,0)=std::numeric_limits<double>::denorm_min();subnormal.makeCompressed();
    compare(mesh,layout,mass,energy,geometry,subnormal,false,nativePlan,"actual subnormal full-expression fallback");
    {
        RestoreFloating restore;check(std::fesetround(FE_UPWARD)==0&&std::fegetround()==FE_UPWARD,"sparse fallback cannot set directed FP mode");
        // Original directed factory must really succeed; no exception is
        // swallowed or counted as unsupported-path qualification.
        compare(mesh,layout,mass,energy,geometry,a,false,false,"current upward full-expression fallback");
    }
    const auto foreign=HydroMomentumGeometry::build(mesh,metrics);check(!foreign.sameIdentity(mass.geometry()),"foreign owner witness is the actual mass owner");
    invalid([&]{HydroSaddleVelocityBlock::buildWithEnergyTransportSplit(mesh,layout,a,mass,energy,.125,foreign);});
    auto bad=a;bad.coeffRef(0,0)=std::numeric_limits<double>::infinity();
    invalid([&]{HydroSaddleVelocityBlock::buildWithEnergyTransportSplit(mesh,layout,bad,mass,energy,.125);});
    invalid([&]{HydroSaddleVelocityBlock::buildWithEnergyTransportSplit(mesh,layout,bad,mass,energy,.125,geometry);});
    bad=a;bad.coeffRef(0,0)=std::numeric_limits<double>::max();
    invalid([&]{HydroSaddleVelocityBlock::buildWithEnergyTransportSplit(mesh,layout,bad,mass,energy,.125);});
    invalid([&]{HydroSaddleVelocityBlock::buildWithEnergyTransportSplit(mesh,layout,bad,mass,energy,.125,geometry);});
    const auto pressure=HydroPressureLoad::build(metrics.divergence());
    const auto solver=HydroSinglePhaseSaddle::buildClosed(pressure,reused);
    const std::vector<double> zero(layout.velocityCount());const auto result=solver.solve(mesh,HydroSaddleRhs(layout,zero,zero));
    check(result.values().normalVelocity()==zero&&result.values().pressure()==std::vector<double>(layout.pressureCount()),"sparse validation changed exact zero physical solve");
    const auto& rr=result.report().velocitySparseValidation;const auto& vr=reused.sparseValidationReport();
    check(rr.requested==vr.requested&&rr.eligible==vr.eligible&&rr.planHit==vr.planHit&&rr.exactSymmetryUsed==vr.exactSymmetryUsed
        &&rr.exactSkewUsed==vr.exactSkewUsed&&rr.nonZeros==vr.nonZeros&&rr.retainedBytes==vr.retainedBytes
        &&rr.appendPlanUsed==vr.appendPlanUsed&&rr.appendScratchBytes==vr.appendScratchBytes,
        "actual Saddle report lost the consumed sparse validation receipt");++actualSolveChecks;
}
void run(){scene(false);scene(true);check(cases==16&&negativeCases==10&&flagComparisons==16&&actualSolveChecks==2,"sparse validation finite scene coverage incomplete");
#if defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE) && EIGEN_WORLD_VERSION==3 && EIGEN_MAJOR_VERSION==4 && EIGEN_MINOR_VERSION==1 && !defined(EIGEN_SCALAR_BINARY_OP_PLUGIN)
    check(appendChecks==10,"Task39 required actual append numerical/FP cases missing");
#else
    check(appendChecks==0,"Task39 unsupported compile environment routed append");
#endif
    std::cout<<"VELOCITY_APPEND_PLAN_PASS,actual_append_original_CSR_norm_FP_pairs,"<<appendChecks
        <<",literal_and_unsupported_pairs,"<<cases-appendChecks<<",accepted_saddle_receipts,"<<actualSolveChecks
        <<",resource_cap_not_injected,no_speed_claim\n";
    std::cout<<"SPARSE_VALIDATION_PLAN_PASS,cases,"<<cases<<",coefficient_bits,"<<coefficientChecks<<",original_FP_comparisons,"<<flagComparisons
        <<",input_owner_rejections,"<<negativeCases<<",actual_zero_solves,"<<actualSolveChecks<<",64MiB_cap_dynamic_not_exercised,no_speed_claim\n";
}
}

namespace SaddleScatterPlanContract {
using System=HydroSinglePhaseSaddle::Matrix;
std::size_t scenes=0,cases=0,coefficientChecks=0,fieldChecks=0,independentSolves=0,negativeCases=0,flagComparisons=0;
std::size_t coldThreadGroups=0,coldThreadPairs=0;
std::uint64_t bits(double x){std::uint64_t result=0;std::memcpy(&result,&x,sizeof result);return result;}
#if defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE) && EIGEN_WORLD_VERSION==3 && EIGEN_MAJOR_VERSION==4 && EIGEN_MINOR_VERSION==1
constexpr bool nativePlan=true;
#else
constexpr bool nativePlan=false;
#endif
void sameMatrix(const System& a,const System& b) {
    check(a.isCompressed()&&b.isCompressed()&&a.rows()==b.rows()&&a.cols()==b.cols()&&a.nonZeros()==b.nonZeros(),
        "scatter plan changed complete augmented matrix shape/storage");
    for(Eigen::Index i=0;i<=a.outerSize();++i)check(a.outerIndexPtr()[i]==b.outerIndexPtr()[i],"scatter plan changed exact augmented outer indices");
    for(Eigen::Index i=0;i<a.nonZeros();++i){check(a.innerIndexPtr()[i]==b.innerIndexPtr()[i]&&bits(a.valuePtr()[i])==bits(b.valuePtr()[i]),
        "scatter plan changed augmented coefficient bits or explicit signed zero");++coefficientChecks;}
}
void sameValues(const std::vector<double>& a,const std::vector<double>& b) {
    check(a.size()==b.size(),"scatter plan changed original physical field extent");
    for(std::size_t i=0;i<a.size();++i){check(bits(a[i])==bits(b[i]),"scatter plan changed old-factory deterministic physical field bits");++fieldChecks;}
}
void originalGates(const HydroSaddleResult& result) {
    const auto& r=result.report();check(r.algebraicRelativeResidual<=1e-11&&r.freeMomentumDefect<=1e-11&&r.maxCellContinuityDefect<=1e-11
        &&r.maxPressureMeanDefect<=1e-11&&r.maxBoundaryFluxDefect<=1e-12,"scatter solve escaped an original physical/algebraic gate");
}
struct Trial {HydroSinglePhaseSaddle solver;HydroSaddleResult result;};
Trial compare(const QuadtreeMesh& mesh,const HydroDofLayout& layout,const HydroPressureLoad& load,const HydroMomentumGeometry& owner,
    const Reference& ref,const Matrix& actualA,const std::vector<HydroSaddlePressureOutlet>& outlets,
    bool expectedEligible,bool expectedHit,bool expectedScatter,const char* name,bool restoreNearestForSolve=false,
    SparseValidationPlanContract::Floating* factoryCompleted=nullptr) {
    SparseValidationPlanContract::RestoreFloating restore;
    const auto velocity=HydroSaddleVelocityBlock::build(layout,actualA);
    check(std::feraiseexcept(FE_INEXACT)==0,"scatter plan cannot prime original sticky inexact");
    const SparseValidationPlanContract::Floating before;
    const auto old=outlets.empty()?HydroSinglePhaseSaddle::buildClosed(load,velocity):HydroSinglePhaseSaddle::buildWithPressureOutlets(load,velocity,outlets);
    const SparseValidationPlanContract::Floating oldAfter;before.restore();
    const auto prepared=outlets.empty()?HydroSinglePhaseSaddle::buildClosed(load,velocity,HydroSaddleConfig{},HydroSaddlePreparationTimings{},owner)
        :HydroSinglePhaseSaddle::buildWithPressureOutlets(load,velocity,outlets,HydroSaddleConfig{},HydroSaddlePreparationTimings{},owner);
    const SparseValidationPlanContract::Floating newAfter;
    if(factoryCompleted)*factoryCompleted=newAfter; // Capture before any oracle/comparison arithmetic.
    SparseValidationPlanContract::sameFlags(oldAfter,newAfter);++flagComparisons;
    sameMatrix(old.matrix(),prepared.matrix());
    check(old.freeVelocityDofs()==prepared.freeVelocityDofs()&&old.pressureComponents()==prepared.pressureComponents()
        &&old.pressureGaugeByComponent()==prepared.pressureGaugeByComponent()&&old.pressureGaugeCount()==prepared.pressureGaugeCount(),
        "scatter plan changed actual free partition, component order or pressure gauges");
    sameValues(old.pressureOutletLoads(),prepared.pressureOutletLoads());
    check(old.pressureOutlets().size()==outlets.size()&&prepared.pressureOutlets().size()==outlets.size(),"scatter plan changed prescribed outlet inventory");
    for(std::size_t i=0;i<outlets.size();++i)check(prepared.pressureOutlets()[i].velocityDof==outlets[i].velocityDof
        &&bits(prepared.pressureOutlets()[i].pressure)==bits(outlets[i].pressure),"scatter plan froze old outlet roles/pressure values");
    if(restoreNearestForSolve)check(std::fesetround(FE_TONEAREST)==0,"scatter unsupported-FP fixture cannot restore nearest for the physical solve");
    // Independent geometry incidence and physical cell volumes assemble the
    // ORIGINAL augmented equations. No private quotient map or plan enters it.
    const int nu=layout.velocityCount(),np=layout.pressureCount();
    std::vector<int> free,byVelocity(nu,-1);std::vector<bool> isOutlet(nu,false);
    std::vector<double> traction(nu),u(nu),p(np),known(nu),weights(nu,1.);
    for(const auto& outlet:outlets){check(outlet.velocityDof>=0&&outlet.velocityDof<nu,"scatter fixture has a bad outlet role");
        const auto& f=ref.faces[outlet.velocityDof];check(f.incidence==1&&f.horizontal&&f.y0==mesh.maxCoord.y&&f.area>0,
            "scatter outlet fixture is not an actual top physical boundary");
        isOutlet[outlet.velocityDof]=true;traction[outlet.velocityDof]=-f.area*outlet.pressure;}
    const bool rz=layout.domain().axisymmetric();
    for(int i=0;i<nu;++i){const auto& f=ref.faces[i];
        const double flow=f.horizontal?polynomialPsi(f.x0,f.y0,mesh)-polynomialPsi(f.x1,f.y0,mesh)
            :polynomialPsi(f.x0,f.y1,mesh)-polynomialPsi(f.x0,f.y0,mesh);
        u[i]=f.area==0?0:(rz?2*pi:1)*flow/f.area;
        if(!outlets.empty()&&f.horizontal)u[i]+=.125;
        if(f.incidence==2||isOutlet[i]){byVelocity[i]=static_cast<int>(free.size());free.push_back(i);}else known[i]=u[i];
    }
    check(free==prepared.freeVelocityDofs()&&free.size()>=5,"scatter fixture changed independent free normal classification");
    double mean=0,volume=0;
    for(int i=0;i<np;++i){const auto c=ref.center[i];p[i]=1+.3*c.x+.7*c.y+.4*c.x*c.y;mean+=ref.volume[i]*p[i];volume+=ref.volume[i];}
    if(outlets.empty())for(double& value:p)value-=mean/volume;
    const int nf=static_cast<int>(free.size()),ng=outlets.empty()?1:0,size=nf+np+ng;
    check(prepared.pressureGaugeCount()==ng&&prepared.pressureComponents()==std::vector<int>(np,0),"scatter fixture lost single connected physical component/gauge role");
    Eigen::MatrixXd dense=Eigen::MatrixXd::Zero(size,size);
    for(int i=0;i<nf;++i)for(int j=0;j<nf;++j)dense(i,j)=actualA.coeff(free[i],free[j]);
    for(const auto& t:ref.terms)if(byVelocity[t.u]>=0){const int i=byVelocity[t.u];dense(i,nf+t.p)-=t.outwardArea;dense(nf+t.p,i)-=t.outwardArea;}
    if(ng)for(int i=0;i<np;++i)dense(nf+i,nf+np)=dense(nf+np,nf+i)=ref.volume[i];
    const Eigen::MatrixXd actual(prepared.matrix());check((actual-dense).norm()<=1e-11*dense.norm(),"scatter matrix changed independent original physical equations");
    auto momentum=multiply(actualA,u);for(const auto& t:ref.terms)momentum[t.u]-=t.outwardArea*p[t.p];
    for(int i=0;i<nu;++i)momentum[i]-=traction[i];
    const auto eliminated=multiply(actualA,known);Eigen::VectorXd b=Eigen::VectorXd::Zero(size);
    for(int i=0;i<nf;++i)b[i]=momentum[free[i]]+traction[free[i]]-eliminated[free[i]];
    for(const auto& t:ref.terms)b[nf+t.p]+=t.outwardArea*known[t.u];
    check(b.stableNorm()>0,"scatter contract accidentally used a zero RHS");
    const Eigen::FullPivLU<Eigen::MatrixXd> oracle(dense);check(oracle.isInvertible(),"scatter independent original augmented matrix is singular");
    const Eigen::VectorXd exact=oracle.solve(b);check(exact.allFinite()&&(dense*exact-b).stableNorm()<=1e-11*b.stableNorm(),
        "scatter independent original augmented solve failed its residual gate");
    auto exactU=known;std::vector<double> exactP(np);for(int i=0;i<nf;++i)exactU[free[i]]=exact[i];for(int i=0;i<np;++i)exactP[i]=exact[nf+i];
    const HydroSaddleRhs rhs(layout,momentum,known);
    const auto oldResult=old.solve(mesh,rhs),result=prepared.solve(mesh,rhs);originalGates(oldResult);originalGates(result);
    sameValues(oldResult.values().normalVelocity(),result.values().normalVelocity());sameValues(oldResult.values().pressure(),result.values().pressure());
    check(relativeError(result.values().normalVelocity(),exactU,weights)<=1e-9&&relativeError(result.values().pressure(),exactP,ref.volume)<=1e-9
        &&relativeError(result.values().normalVelocity(),u,weights)<=1e-9&&relativeError(result.values().pressure(),p,ref.volume)<=1e-9,
        "scatter solution differs from independent dense or manufactured physical state");
    const auto& r=result.report().sparseAssemblyPlan;const auto& noRequest=oldResult.report().sparseAssemblyPlan;
    check(!noRequest.requested&&!noRequest.eligible&&!noRequest.planHit&&!noRequest.systemScatterUsed&&!noRequest.quotientScatterUsed,
        "old factory silently entered typed scatter planning");
    const bool aliases=expectedEligible&&!outlets.empty();
    check(r.requested&&r.eligible==expectedEligible&&r.planHit==expectedHit&&r.systemScatterUsed==expectedScatter&&
        r.quotientScatterUsed==(expectedScatter&&!aliases)&&r.quotientAliasesSystem==aliases,
        "scatter report claims a shortcut/hit inconsistent with actual pattern/role/FP qualification");
    check(r.systemNonZeros==static_cast<std::size_t>(prepared.matrix().nonZeros())&&r.quotientNonZeros==static_cast<std::size_t>(result.report().factorNonZeros),
        "scatter report lost actual original system/quotient stored counts");
    if(expectedEligible)check(r.retainedBytes>0&&r.retainedBytes<=128u*1024u*1024u,"scatter plan storage is not actually bounded");
    check(!r.symbolicBuilt||(!r.planHit&&r.eligible),"cold scatter compilation was mislabeled as a hit");
    check(r.symbolicWorkers>=0&&r.scatterMaximumWorkers>=0&&
        r.scatterParallelUsed==(r.scatterMaximumWorkers>1),"scatter receipt fabricates an actual numeric team");
    check(!noRequest.symbolicBuilt&&!noRequest.quotientAliasesSystem&&!noRequest.scatterParallelUsed&&
        !noRequest.scatterParallelFallback&&noRequest.scatterMaximumWorkers==0,"old public literal factory entered typed pipeline");
    check(result.report().factorDimension==oldResult.report().factorDimension&&result.report().factorNonZeros==oldResult.report().factorNonZeros,
        "scatter path changed original private quotient dimensions/stored count");
    ++cases;++independentSolves;
    std::cout<<"SADDLE_SCATTER_PLAN,"<<name<<",rz,"<<rz<<",outlets,"<<outlets.size()<<",hit,"<<r.planHit
        <<",system_scatter,"<<r.systemScatterUsed<<",quotient_scatter,"<<r.quotientScatterUsed<<",nonzero_RHS_independent_original_oracle,PASS\n";
    return {prepared,result};
}
Matrix specimen(int n,const std::vector<int>& free,int partner=1) {
    check(free.size()>4,"scatter specimen lacks independent actual free rows");Matrix a=diagonal(std::vector<double>(n,4.));
    a.coeffRef(free[0],free[partner])=-.125;a.coeffRef(free[partner],free[0])=-.125;
    a.coeffRef(free[2],free[3])=0.;a.coeffRef(free[3],free[2])=-0.;a.makeCompressed();return a;
}
template<class F>void invalid(F action){bool caught=false;try{action();}catch(const std::invalid_argument&){caught=true;}
    check(caught,"typed scatter accepted a foreign actual layout or divergence owner");++negativeCases;}
void scene(bool rz,bool outletMode) {
    SparseValidationPlanContract::RestoreFloating restore;check(std::fesetround(FE_TONEAREST)==0,"scatter fixture requires nearest mode");
    QuadtreeMesh mesh(input(4,rz?1:0));refine(mesh,1);FaceSegmentGraph graph;graph.build(mesh);
    const auto layout=HydroDofLayout::build(mesh,NativeLeafDomain::build(mesh,graph,rz));
    check(!layout.domain().uniformLeaves(),"scatter fixture lost actual nonuniform leaves");
    const auto integral=HydroIntegralDivergence::build(mesh,layout);const auto metrics=HydroPhysicalDivergence::build(mesh,integral);
    const auto load=HydroPressureLoad::build(integral);const auto geometry=HydroMomentumGeometry::build(mesh,metrics);
    const auto ref=reference(mesh,graph,layout,rz);std::vector<int> free;std::vector<HydroSaddlePressureOutlet> top,outlets;
    for(int i=0;i<layout.velocityCount();++i){const auto& f=ref.faces[i];if(f.incidence==2)free.push_back(i);
        if(f.incidence==1&&f.horizontal&&f.y0==mesh.maxCoord.y&&f.area>0)top.push_back({i,2.5});}
    check(top.size()>1,"scatter fixture lacks ordered multiple physical outlet roles");if(outletMode)outlets=top;
    const auto a=specimen(layout.velocityCount(),free);
    check(std::signbit(a.coeff(free[3],free[2]))&&a.nonZeros()==layout.velocityCount()+4,"scatter fixture discarded explicit signed zero");
    const auto cold=compare(mesh,layout,load,geometry,ref,a,outlets,nativePlan,false,nativePlan,"cold integer route used in current call");
    auto fresh=a;for(int i=0;i<layout.velocityCount();++i)fresh.coeffRef(i,i)=5.;
    fresh.coeffRef(free[0],free[1])=-.25;fresh.coeffRef(free[1],free[0])=-.25;fresh.makeCompressed();
    auto freshOutlets=outlets;for(auto& outlet:freshOutlets)outlet.pressure+=.25;
    const auto warm=compare(mesh,layout,load,geometry,ref,fresh,freshOutlets,nativePlan,nativePlan,nativePlan,"same stored pattern, fresh A and outlet pressure");
    check(bits(cold.solver.matrix().coeff(0,0))!=bits(warm.solver.matrix().coeff(0,0)),"scatter plan froze previous velocity values");
    if(outletMode)check(cold.solver.pressureOutletLoads()!=warm.solver.pressureOutletLoads(),"scatter plan froze previous physical traction values");
    auto zeros=a;for(Eigen::Index row=0;row<zeros.outerSize();++row)for(Matrix::InnerIterator e(zeros,row);e;++e)
        if(e.row()!=e.col())e.valueRef()=e.row()<e.col()?0.:-0.;
    check(zeros.nonZeros()==a.nonZeros(),"scatter zero scene changed stored pattern");
    compare(mesh,layout,load,geometry,ref,zeros,outlets,nativePlan,nativePlan,nativePlan,"fresh explicit plus/minus zero values");
    const auto other=specimen(layout.velocityCount(),free,4);
    check(other.nonZeros()==a.nonZeros()&&(!std::equal(a.outerIndexPtr(),a.outerIndexPtr()+a.outerSize()+1,other.outerIndexPtr())
        ||!std::equal(a.innerIndexPtr(),a.innerIndexPtr()+a.nonZeros(),other.innerIndexPtr())),"scatter foreign pattern witness did not change complete pattern at equal NNZ");
    const auto mismatch=compare(mesh,layout,load,geometry,ref,other,outlets,nativePlan,false,nativePlan,"same count, changed actual free-free pattern compiles current replacement");
    check((mismatch.solver.matrix()-cold.solver.matrix()).norm()>0,"scatter foreign pattern changed only discarded prescribed rows");
    compare(mesh,layout,load,geometry,ref,other,outlets,nativePlan,nativePlan,nativePlan,"replacement pattern now reuses both scatters");
    compare(mesh,layout,load,geometry,ref,a,outlets,nativePlan,false,nativePlan,"returning old pattern recompiles and uses the single slot");
    // Another real geometry with identical layout and authoritative D is legal;
    // its independent empty slot compiles cold, rather than failing by address.
    const auto independentOwner=HydroMomentumGeometry::build(mesh,metrics);
    check(!independentOwner.sameIdentity(geometry),"scatter independent owner scene accidentally used the same geometry Data");
    compare(mesh,layout,load,independentOwner,ref,a,outlets,nativePlan,false,nativePlan,"same actual layout/D, independent legitimate cold owner");
    auto changedRoles=top;if(outletMode)std::reverse(changedRoles.begin(),changedRoles.end());
    compare(mesh,layout,load,geometry,ref,a,changedRoles,nativePlan,false,nativePlan,outletMode?"changed ordered outlet roles":"closed to physical outlet and gauge role change");
    if(!outletMode){SparseValidationPlanContract::RestoreFloating directed;check(std::fesetround(FE_UPWARD)==0&&std::fegetround()==FE_UPWARD,
            "scatter fallback cannot set actual directed mode");
        // Both original factories must truly succeed in this environment.
        // Restore RNE only for the subsequent independent physical solve.
        compare(mesh,layout,load,geometry,ref,a,outlets,false,false,false,"unsupported current FP retains original packing",true);}
    if(rz&&outletMode){
        SparseValidationPlanContract::RestoreFloating threadRestore;
#ifdef _OPENMP
        const int savedThreads=omp_get_max_threads(),savedDynamic=omp_get_dynamic();
        struct ThreadRestore {int count,dynamic;~ThreadRestore(){omp_set_num_threads(count);omp_set_dynamic(dynamic);}} keep{savedThreads,savedDynamic};
        omp_set_dynamic(0);
#endif
        const SparseValidationPlanContract::Floating entry;
        std::optional<Trial> first;std::optional<SparseValidationPlanContract::Floating> firstFlags;
        for(int requested:{1,2,8}){
#ifdef _OPENMP
            omp_set_num_threads(requested);
#endif
            const auto freshOwner=HydroMomentumGeometry::build(mesh,metrics);
            check(!freshOwner.sameIdentity(geometry),"cold scatter thread fixture reused the original owner");
            entry.restore();SparseValidationPlanContract::Floating completed;
            auto result=compare(mesh,layout,load,freshOwner,ref,a,outlets,nativePlan,false,nativePlan,
                "actual new-owner cold route own 1/2/8",false,&completed);
            const auto& actual=result.result.report().sparseAssemblyPlan;
            check(!actual.planHit&&actual.symbolicBuilt==nativePlan,"fresh actual owner did not report current cold compilation");
            if(first){sameMatrix(first->solver.matrix(),result.solver.matrix());
                sameValues(first->result.values().normalVelocity(),result.result.values().normalVelocity());
                sameValues(first->result.values().pressure(),result.result.values().pressure());
                SparseValidationPlanContract::sameFlags(*firstFlags,completed);++coldThreadPairs;}
            else{first=std::move(result);firstFlags=completed;}
        }
        ++coldThreadGroups;
    }
    const auto foreignLayout=HydroDofLayout::build(mesh,NativeLeafDomain::build(mesh,graph,rz));
    const auto foreignMetrics=HydroPhysicalDivergence::build(mesh,HydroIntegralDivergence::build(mesh,foreignLayout));
    const auto foreignGeometry=HydroMomentumGeometry::build(mesh,foreignMetrics);
    const auto velocity=HydroSaddleVelocityBlock::build(layout,a);
    invalid([&]{HydroSinglePhaseSaddle::buildClosed(load,velocity,HydroSaddleConfig{},HydroSaddlePreparationTimings{},foreignGeometry);});
    const auto anotherIntegral=HydroIntegralDivergence::build(mesh,layout);
    check(&anotherIntegral.matrix()!=&integral.matrix(),"scatter foreign divergence witness reused the actual D owner");
    const auto anotherGeometry=HydroMomentumGeometry::build(mesh,HydroPhysicalDivergence::build(mesh,anotherIntegral));
    invalid([&]{HydroSinglePhaseSaddle::buildClosed(load,velocity,HydroSaddleConfig{},HydroSaddlePreparationTimings{},anotherGeometry);});
    ++scenes;
}
void run(){scene(false,false);scene(false,true);scene(true,false);scene(true,true);
    check(scenes==4&&cases==37&&independentSolves==37&&negativeCases==8&&flagComparisons==37&&
        coldThreadGroups==1&&coldThreadPairs==2,"scatter finite cold/own-thread coverage is incomplete");
    std::cout<<"SADDLE_SCATTER_PLAN_PASS,scenes,"<<scenes<<",old_new_nonzero_RHS_pairs,"<<cases<<",augmented_coefficient_bits,"<<coefficientChecks
        <<",physical_field_bits,"<<fieldChecks<<",independent_original_solves,"<<independentSolves<<",FP_comparisons,"<<flagComparisons
        <<",foreign_layout_D_rejections,"<<negativeCases<<",actual_new_owner_cold_thread_groups,"<<coldThreadGroups
        <<",own_1_2_8_FP_pairs,"<<coldThreadPairs<<",128MiB_cap_dynamic_not_exercised,no_speed_claim\n";
}
}
namespace VelocitySplitPlanContract {
using SparseValidationPlanContract::Floating;
using SparseValidationPlanContract::RestoreFloating;
using SparseValidationPlanContract::bits;
std::size_t scenes=0,pairs=0,rejections=0,physicalSolves=0;
#if defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE) && !defined(__FAST_MATH__) && !defined(_M_FP_FAST) && EIGEN_WORLD_VERSION==3 && EIGEN_MAJOR_VERSION==4 && EIGEN_MINOR_VERSION==1 && !defined(EIGEN_SCALAR_BINARY_OP_PLUGIN)
constexpr bool nativePlan=true;
#else
constexpr bool nativePlan=false;
#endif
Matrix specimen(int n,bool reversedMissingSide=false) {
    check(n>8,"split fixture lacks actual velocity rows");Matrix a=diagonal(std::vector<double>(n,4.));
    a.coeffRef(0,1)=.375;a.coeffRef(1,0)=.125;
    if(reversedMissingSide)a.coeffRef(3,2)=-.125;else a.coeffRef(2,3)=-.125;
    a.coeffRef(5,4)=-0.;a.makeCompressed();return a;
}
void sameReceipt(const HydroVelocitySplitPlanReport& a,const HydroVelocitySplitPlanReport& b) {
    check(a.requested==b.requested&&a.eligible==b.eligible&&a.planHit==b.planHit
        &&a.symmetricValuesUsed==b.symmetricValuesUsed&&a.skewValuesUsed==b.skewValuesUsed
        &&a.inputNonZeros==b.inputNonZeros&&a.unionNonZeros==b.unionNonZeros&&a.retainedBytes==b.retainedBytes
        &&a.parallelUsed==b.parallelUsed&&a.parallelFallback==b.parallelFallback&&a.maximumWorkers==b.maximumWorkers,
        "actual Saddle result lost the actual split receipt");
}
HydroSaddleVelocityBlock compare(const QuadtreeMesh& mesh,const HydroDofLayout& layout,const HydroMomentumMass& mass,
    const HydroEnergyViscosity& energy,const HydroMomentumGeometry& owner,const Matrix& a,
    bool eligible,bool hit,bool used,const char* name) {
    RestoreFloating restore;
    const Matrix saved=a,at(a.transpose()),s=.5*(a+at),w=.5*(a-at),combined=s+w;
    check(s.nonZeros()==w.nonZeros(),"literal split union patterns differ");
    const Matrix transpose(combined.transpose()),difference=combined-transpose;
    const double expected=difference.norm()/combined.norm();
    check(std::feraiseexcept(FE_INEXACT)==0,"split comparison cannot prime sticky flags");errno=ERANGE;
    const Floating before;
    const auto old=HydroSaddleVelocityBlock::buildWithEnergyTransportSplit(mesh,layout,a,mass,energy,.125);const Floating oldAfter;
    before.restore();
    const auto fresh=HydroSaddleVelocityBlock::buildWithEnergyTransportSplit(mesh,layout,a,mass,energy,.125,owner);const Floating freshAfter;
    SparseValidationPlanContract::sameFlags(oldAfter,freshAfter);
    SparseValidationPlanContract::sameMatrix(old.matrix(),fresh.matrix());
    SparseValidationPlanContract::sameMatrix(fresh.matrix(),combined);
    SparseValidationPlanContract::sameMatrix(a,saved);
    check(bits(old.symmetryDefect())==bits(expected)&&bits(fresh.symmetryDefect())==bits(expected),
        "split changed the complete original final symmetry norm");
    check(old.spdReport().tier!=HydroVelocitySpdTier::None&&fresh.spdReport().tier!=HydroVelocitySpdTier::None,
        "explicit public S/W factories lost their requested original SPD audit after production raw-A separation");
    const auto& r=fresh.splitPlanReport();
    check(r.requested&&r.eligible==eligible&&r.planHit==hit&&r.symmetricValuesUsed==used&&r.skewValuesUsed==used,
        "split report does not describe the actual source/pattern/FP path");
    check(r.inputNonZeros==static_cast<std::size_t>(a.nonZeros())&&r.unionNonZeros==static_cast<std::size_t>(s.nonZeros()),
        "split report confused original A and full S/W union stored counts");
    if(eligible)check(r.retainedBytes>0&&r.retainedBytes<=128u*1024u*1024u,"split retained stock exceeds its actual bound");
    const auto& legacy=old.splitPlanReport();check(!legacy.requested&&!legacy.eligible&&!legacy.planHit
        &&!legacy.symmetricValuesUsed&&!legacy.skewValuesUsed&&legacy.retainedBytes==0,
        "six-argument old factory silently requested split planning");
    ++pairs;std::cout<<"VELOCITY_SPLIT_PLAN,"<<name<<",nnzA,"<<a.nonZeros()<<",nnzUnion,"<<s.nonZeros()
        <<",eligible,"<<r.eligible<<",hit,"<<r.planHit<<",used,"<<r.symmetricValuesUsed<<",literal_bits_and_FP,PASS\n";
    return fresh;
}
template<class F> void invalid(F action) {
    bool rejected=false;try{action();}catch(const std::invalid_argument&){rejected=true;}
    check(rejected,"split bypassed an original source/finite/dimension/mesh guard");++rejections;
}
void scene(bool rz) {
    RestoreFloating restore;check(std::fesetround(FE_TONEAREST)==0,"split fixture cannot select original RNE");
    QuadtreeMesh mesh(input(4,rz?1:0));refine(mesh,1);FaceSegmentGraph graph;graph.build(mesh);
    const auto layout=HydroDofLayout::build(mesh,NativeLeafDomain::build(mesh,graph,rz));
    check(!layout.domain().uniformLeaves(),"split fixture lost real 2:1 AMR");
    const auto metrics=HydroPhysicalDivergence::build(mesh,HydroIntegralDivergence::build(mesh,layout));
    const auto owner=HydroMomentumGeometry::build(mesh,metrics);
    const LAMP::TimeLayer initial{0,0,0,LAMP::TimeLayer::Position::Initial};
    const auto pub=LAMP::Interface::GeometricVofPublication::prepareLocalPlic(mesh,layout.domain(),initial,
        std::vector<double>(mesh.cellIndexNum(),1.));
    const auto mass=HydroMomentumMass::buildCellAverageProjection(mesh,owner,HydroPhaseGeometry::build(mesh,metrics,pub),1.,1.);
    const auto energy=HydroEnergyViscosity::zero(mesh,metrics);const auto a=specimen(layout.velocityCount());
    check(a.nonZeros()==layout.velocityCount()+4&&std::signbit(a.coeff(5,4)),"split missing-side/signedzero witness was pruned");
    const auto cold=compare(mesh,layout,mass,energy,owner,a,nativePlan,false,false,"cold literal then compile");
    auto changed=a;changed.coeffRef(0,1)=.625;changed.coeffRef(1,0)=.375;changed.makeCompressed();
    const auto warm=compare(mesh,layout,mass,energy,owner,changed,nativePlan,nativePlan,nativePlan,"same pattern, fresh values");
    check(bits(cold.matrix().coeff(0,1))!=bits(warm.matrix().coeff(0,1)),"split reused previous physical values");
    auto zeros=a;for(Eigen::Index p=0;p<zeros.outerSize();++p)for(Matrix::InnerIterator e(zeros,p);e;++e)
        if(e.row()!=e.col())e.valueRef()=e.row()<e.col()?0.:-0.;
    compare(mesh,layout,mass,energy,owner,zeros,nativePlan,nativePlan,nativePlan,"stored plus/minuszero and missing positivezero");
    auto halfUnderflow=a;halfUnderflow.coeffRef(0,1)=std::nextafter(std::numeric_limits<double>::min(),std::numeric_limits<double>::infinity());
    halfUnderflow.coeffRef(1,0)=-std::numeric_limits<double>::min();halfUnderflow.makeCompressed();
    compare(mesh,layout,mass,energy,owner,halfUnderflow,nativePlan,nativePlan,nativePlan,"normal inputs with sum-then-half subnormal rounding");
    auto large=a;for(Eigen::Index k=0;k<large.nonZeros();++k)large.valuePtr()[k]=std::ldexp(large.valuePtr()[k],480);
    compare(mesh,layout,mass,energy,owner,large,nativePlan,nativePlan,nativePlan,"finite high exponent");
    auto small=a;for(Eigen::Index k=0;k<small.nonZeros();++k)small.valuePtr()[k]=std::ldexp(small.valuePtr()[k],-480);
    compare(mesh,layout,mass,energy,owner,small,nativePlan,nativePlan,nativePlan,"finite low exponent");
    const auto other=specimen(layout.velocityCount(),true);
    check(other.nonZeros()==a.nonZeros()&&!std::equal(a.innerIndexPtr(),a.innerIndexPtr()+a.nonZeros(),other.innerIndexPtr()),
        "split foreign original pattern witness did not change original input A");
    const Matrix at(a.transpose()),bt(other.transpose()),sa=.5*(a+at),sb=.5*(other+bt);
    check(sa.nonZeros()==sb.nonZeros()&&std::equal(sa.outerIndexPtr(),sa.outerIndexPtr()+sa.outerSize()+1,sb.outerIndexPtr())
        &&std::equal(sa.innerIndexPtr(),sa.innerIndexPtr()+sa.nonZeros(),sb.innerIndexPtr()),
        "split witness did not preserve union while changing actual input A");
    compare(mesh,layout,mass,energy,owner,other,nativePlan,false,false,"same union but different actual A compiles replacement");
    compare(mesh,layout,mass,energy,owner,other,nativePlan,nativePlan,nativePlan,"replacement actual A now warm");
    compare(mesh,layout,mass,energy,owner,a,nativePlan,false,false,"return original A recompiles single slot");
    auto subnormal=a;subnormal.coeffRef(0,1)=std::numeric_limits<double>::denorm_min();subnormal.makeCompressed();
    compare(mesh,layout,mass,energy,owner,subnormal,false,nativePlan,false,"same pattern subnormal literal fallback");
    {RestoreFloating directed;check(std::fesetround(FE_UPWARD)==0,"split cannot select directed mode");
        compare(mesh,layout,mass,energy,owner,a,false,false,false,"unsupported FP literal fallback");}
    // Boundary roles are not input to A -> S/W. Actual Saddle scatter still
    // independently binds each complete ordered outlet/free/gauge inventory.
    const auto pressure=HydroPressureLoad::build(metrics.divergence());const auto ref=reference(mesh,graph,layout,rz);
    std::vector<HydroSaddlePressureOutlet> outlets;
    for(int u=0;u<layout.velocityCount();++u){const auto& f=ref.faces[u];
        if(f.incidence==1&&f.horizontal&&f.y0==mesh.maxCoord.y&&f.area>0)outlets.push_back({u,.25});}
    check(outlets.size()>1,"split role fixture lacks actual multiple pressure outlets");
    HydroSinglePhaseSaddle::buildWithPressureOutlets(pressure,warm,outlets,HydroSaddleConfig{},HydroSaddlePreparationTimings{},owner);
    std::reverse(outlets.begin(),outlets.end());
    const auto roleBlock=compare(mesh,layout,mass,energy,owner,a,nativePlan,nativePlan,nativePlan,"actual outlet roles change keeps pure A plan");
    const auto reversed=HydroSinglePhaseSaddle::buildWithPressureOutlets(pressure,roleBlock,outlets,HydroSaddleConfig{},HydroSaddlePreparationTimings{},owner);
    check(reversed.pressureOutlets().front().velocityDof==outlets.front().velocityDof,"split bypassed actual ordered outlet binding");
    const auto oldBlock=HydroSaddleVelocityBlock::buildWithEnergyTransportSplit(mesh,layout,a,mass,energy,.125);
    const auto oldSolver=HydroSinglePhaseSaddle::buildClosed(pressure,oldBlock);
    const auto newSolver=HydroSinglePhaseSaddle::buildClosed(pressure,roleBlock);
    SaddleScatterPlanContract::sameMatrix(oldSolver.matrix(),newSolver.matrix());
    std::vector<double> exactP(layout.pressureCount()),momentum(layout.velocityCount()),zero(layout.velocityCount());
    double sum=0,capacity=0;for(std::size_t p=0;p<exactP.size();++p){exactP[p]=ref.center[p].x+.37*ref.center[p].y;
        sum+=ref.volume[p]*exactP[p];capacity+=ref.volume[p];}
    for(double& value:exactP)value-=sum/capacity;
    for(const auto& t:ref.terms)momentum[t.u]-=t.outwardArea*exactP[t.p];
    const HydroSaddleRhs rhs(layout,momentum,zero);const auto oldResult=oldSolver.solve(mesh,rhs),newResult=newSolver.solve(mesh,rhs);
    SaddleScatterPlanContract::sameValues(oldResult.values().normalVelocity(),newResult.values().normalVelocity());
    SaddleScatterPlanContract::sameValues(oldResult.values().pressure(),newResult.values().pressure());
    check(relativeError(newResult.values().pressure(),exactP,ref.volume)<=1e-9,"split changed independent pressure-gradient manufactured solution");
    for(double value:newResult.values().normalVelocity())check(std::abs(value)<=1e-11,"split changed exact zero-velocity manufactured state");
    sameReceipt(newResult.report().velocitySplitPlan,roleBlock.splitPlanReport());++physicalSolves;
    Matrix zeroSkew(layout.velocityCount(),layout.velocityCount());
    check(!roleBlock.withSkew(zeroSkew).splitPlanReport().requested,"withSkew falsely reported a fresh input-A split");
    const auto foreign=HydroMomentumGeometry::build(mesh,metrics);
    invalid([&]{HydroSaddleVelocityBlock::buildWithEnergyTransportSplit(mesh,layout,a,mass,energy,.125,foreign);});
    const auto foreignLayout=HydroDofLayout::build(mesh,NativeLeafDomain::build(mesh,graph,rz));
    invalid([&]{HydroSaddleVelocityBlock::buildWithEnergyTransportSplit(mesh,foreignLayout,a,mass,energy,.125,owner);});
    QuadtreeMesh changedMesh(input(4,rz?1:0));refine(changedMesh,2);
    invalid([&]{HydroSaddleVelocityBlock::buildWithEnergyTransportSplit(changedMesh,layout,a,mass,energy,.125,owner);});
    for(double badValue:{std::numeric_limits<double>::infinity(),std::numeric_limits<double>::max()}){
        auto bad=a;bad.coeffRef(0,0)=badValue;
        invalid([&]{HydroSaddleVelocityBlock::buildWithEnergyTransportSplit(mesh,layout,bad,mass,energy,.125);});
        invalid([&]{HydroSaddleVelocityBlock::buildWithEnergyTransportSplit(mesh,layout,bad,mass,energy,.125,owner);});}
    ++scenes;
}
void run(){scene(false);scene(true);check(scenes==2&&pairs==24&&rejections==14&&physicalSolves==2,"split finite coverage is incomplete");
    std::cout<<"VELOCITY_SPLIT_PLAN_PASS,actual_AMR_XY_RZ,"<<scenes<<",literal_matrix_and_FP_pairs,"<<pairs
        <<",original_rejections,"<<rejections<<",nonzero_RHS_manufactured_solves,"<<physicalSolves
        <<",128MiB_cap_dynamic_not_exercised,no_speed_claim\n";
}
}
}
namespace OrderedRowSpmvContract
{
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE) \
    && !defined(__FAST_MATH__) && (!defined(_M_FP_FAST) || !_M_FP_FAST) \
    && !defined(__AVX__) && !defined(__AVX2__) && !defined(__AVX512F__) && !defined(__FMA__) \
    && !defined(EIGEN_SCALAR_BINARY_OP_PLUGIN) \
    && EIGEN_WORLD_VERSION==3 && EIGEN_MAJOR_VERSION==4 && EIGEN_MINOR_VERSION==1
using Floating=SparseValidationPlanContract::Floating;
using Result=HydroSaddleResult;
std::size_t scenes=0,nonzeroPairs=0,parallelPairs=0,fpPairs=0,matrixValues=0,fieldValues=0,
    denseOracles=0,zeroPairs=0,nestedPairs=0,lateThreadPairs=0,unsupportedAcceptedPairs=0,unsupportedRejectedPairs=0,guardPairs=0,orthThreadGroups=0,orthThreadPairs=0;
unsigned entryVariant=0;
struct OmpRestore{int count=omp_get_max_threads(),dynamic=omp_get_dynamic();~OmpRestore()noexcept{omp_set_dynamic(dynamic);omp_set_num_threads(count);}};
void threads(int n){omp_set_dynamic(0);omp_set_num_threads(n);}
void prime(){check(std::feclearexcept(FE_ALL_EXCEPT)==0,"ordered SpMV cannot clear real status");_mm_setcsr(_mm_getcsr()&~0x3fu);
    check(std::feraiseexcept((entryVariant&1u)?FE_UNDERFLOW|FE_INEXACT:FE_INEXACT)==0,"ordered SpMV cannot raise real sticky flags");
    if(entryVariant&1u)_mm_setcsr(_mm_getcsr()|0x2u);errno=(entryVariant++&1u)?EDOM:E2BIG;}
void sameFp(const Floating& a,const Floating& b){SparseValidationPlanContract::sameFlags(a,b);++fpPairs;}
template<class T>struct Outcome{std::optional<T> value;std::string error;int kind=0;};
template<class T,class F>Outcome<T> invoke(F action){Outcome<T> r;try{r.value.emplace(action());}
    catch(const std::invalid_argument& e){r.error=e.what();r.kind=1;}catch(const std::runtime_error& e){r.error=e.what();r.kind=2;}return r;}
struct Pair{HydroSinglePhaseSaddle literal,parallel;};
Pair factories(const HydroPressureLoad& load,const HydroSaddleVelocityBlock& velocity,const std::vector<HydroSaddlePressureOutlet>& outlets){
    threads(1);prime();const Floating before;
    const auto a=outlets.empty()?HydroSinglePhaseSaddle::buildClosed(load,velocity):HydroSinglePhaseSaddle::buildWithPressureOutlets(load,velocity,outlets);const Floating afterA;
    threads(2);before.restore();const auto b=outlets.empty()?HydroSinglePhaseSaddle::buildClosed(load,velocity):HydroSinglePhaseSaddle::buildWithPressureOutlets(load,velocity,outlets);const Floating afterB;
    sameFp(afterA,afterB);SaddleScatterPlanContract::sameMatrix(a.matrix(),b.matrix());matrixValues+=static_cast<std::size_t>(a.matrix().nonZeros());
    check(a.freeVelocityDofs()==b.freeVelocityDofs()&&a.pressureComponents()==b.pressureComponents()&&a.pressureGaugeByComponent()==b.pressureGaugeByComponent()
        &&a.pressureGaugeCount()==b.pressureGaugeCount(),"ordered SpMV changed original partition/component/gauge equations");
    SaddleScatterPlanContract::sameValues(a.pressureOutletLoads(),b.pressureOutletLoads());return {a,b};}
struct Manufactured{std::vector<double> momentum,known,u,p,weights;};
Manufactured manufactured(const QuadtreeMesh& mesh,const HydroDofLayout& layout,const Reference& ref,const Matrix& a,
    const HydroSinglePhaseSaddle& solver,const std::vector<HydroSaddlePressureOutlet>& outlets,double scale){
    const int nu=layout.velocityCount(),np=layout.pressureCount(),ng=outlets.empty()?1:0;Manufactured m;
    m.known.resize(nu);m.u.resize(nu);m.p.resize(np);m.weights.assign(nu,1.);std::vector<bool> isOutlet(nu);std::vector<double> traction(nu);std::vector<int> free,byVelocity(nu,-1);
    for(const auto& o:outlets){const auto& f=ref.faces[o.velocityDof];check(f.incidence==1&&f.area>0&&f.horizontal&&f.y0==mesh.maxCoord.y,"ordered SpMV fixture outlet is not physical top normal");
        isOutlet[o.velocityDof]=true;traction[o.velocityDof]=-f.area*o.pressure;}
    for(int i=0;i<nu;++i){const auto& face=ref.faces[i];const double flow=face.horizontal?polynomialPsi(face.x0,face.y0,mesh)-polynomialPsi(face.x1,face.y0,mesh)
        :polynomialPsi(face.x0,face.y1,mesh)-polynomialPsi(face.x0,face.y0,mesh);
        m.u[i]=scale*(face.area==0?0:(layout.domain().axisymmetric()?2*pi:1)*flow/face.area);if(!outlets.empty()&&face.horizontal)m.u[i]+=scale*.125;
        if(face.incidence==2||isOutlet[i]){byVelocity[i]=static_cast<int>(free.size());free.push_back(i);}else m.known[i]=m.u[i];}
    check(free==solver.freeVelocityDofs()&&free.size()>5,"ordered SpMV independent free normal classification changed");double mean=0,volume=0;
    for(int p=0;p<np;++p){const auto c=ref.center[p];m.p[p]=scale*(1+.3*c.x+.7*c.y+.4*c.x*c.y);mean+=ref.volume[p]*m.p[p];volume+=ref.volume[p];}
    if(ng)for(double& value:m.p)value-=mean/volume;
    const int nf=static_cast<int>(free.size()),size=nf+np+ng;Eigen::MatrixXd dense=Eigen::MatrixXd::Zero(size,size);
    for(int i=0;i<nf;++i)for(int j=0;j<nf;++j)dense(i,j)=a.coeff(free[i],free[j]);
    for(const auto& t:ref.terms)if(byVelocity[t.u]>=0){dense(byVelocity[t.u],nf+t.p)-=t.outwardArea;dense(nf+t.p,byVelocity[t.u])-=t.outwardArea;}
    if(ng)for(int p=0;p<np;++p)dense(nf+p,nf+np)=dense(nf+np,nf+p)=ref.volume[p];
    const Eigen::MatrixXd actual(solver.matrix());check((actual-dense).norm()<=1e-11*dense.norm(),"ordered SpMV changed independent original augmented equations");
    m.momentum=multiply(a,m.u);for(const auto& t:ref.terms)m.momentum[t.u]-=t.outwardArea*m.p[t.p];for(int i=0;i<nu;++i)m.momentum[i]-=traction[i];
    const auto eliminated=multiply(a,m.known);Eigen::VectorXd rhs=Eigen::VectorXd::Zero(size);
    for(int i=0;i<nf;++i)rhs[i]=m.momentum[free[i]]+traction[free[i]]-eliminated[free[i]];for(const auto& t:ref.terms)rhs[nf+t.p]+=t.outwardArea*m.known[t.u];
    check(rhs.stableNorm()>0,"ordered SpMV manufactured RHS is zero");const Eigen::FullPivLU<Eigen::MatrixXd> oracle(dense);check(oracle.isInvertible(),"ordered SpMV independent augmented matrix is singular");
    const Eigen::VectorXd exact=oracle.solve(rhs);check(exact.allFinite()&&(dense*exact-rhs).stableNorm()<=1e-11*rhs.stableNorm(),"ordered SpMV independent dense residual gate failed");
    auto denseU=m.known;std::vector<double> denseP(np);for(int i=0;i<nf;++i)denseU[free[i]]=exact[i];for(int p=0;p<np;++p)denseP[p]=exact[nf+p];
    check(relativeError(denseU,m.u,m.weights)<=1e-9&&relativeError(denseP,m.p,ref.volume)<=1e-9,"ordered SpMV dense/manufactured solution disagreement");++denseOracles;return m;}
void cgs2Receipt(const HydroSaddleReport& r){
    check(r.orthogonalizationParallelCalls>=0&&r.orthogonalizationParallelCalls<=r.krylovIterations
        &&r.orthogonalizationMaximumWorkers>=0
        &&(r.orthogonalizationParallelCalls>0)==(r.orthogonalizationMaximumWorkers>1),"CGS2 receipt fabricates successful parallel work");
    if(r.usedBlockGmres&&r.krylovIterations>0)check(r.orthogonalizationMaximumWorkers>=1,"CGS2 lost actual serial execution");
    if(r.solveBreakdownTimingsAvailable){const double times[]{r.krylovPreconditionApplySeconds,r.krylovSpmvSeconds,
        r.krylovOrthogonalizationSeconds,r.krylovTrueResidualSeconds,r.physicalAssessmentSeconds};
        for(double t:times)check(std::isfinite(t)&&t>=0,"CGS2 available solve breakdown is nonfinite or negative");}
}
void oldReport(const Result& result){SaddleScatterPlanContract::originalGates(result);const auto& r=result.report();
    check(r.usedBlockGmres&&!r.usedDirectFallback&&r.krylovIterations>0&&r.krylovSolves>0,"ordered SpMV nonzero fixture did not execute genuine accepted GMRES iterations");
    check(std::isfinite(r.krylovRelativeResidual)&&std::isfinite(r.krylovScaledRelativeResidual)&&r.krylovScaledRelativeResidual<=1e-12,
        "CGS2 accepted only an Arnoldi estimate without the original true residual target");cgs2Receipt(r);}
void sameResult(const Result& a,const Result& b){
    SaddleScatterPlanContract::sameValues(a.values().normalVelocity(),b.values().normalVelocity());SaddleScatterPlanContract::sameValues(a.values().pressure(),b.values().pressure());
    SaddleScatterPlanContract::sameValues(a.tangentialVelocity(),b.tangentialVelocity());fieldValues+=a.values().normalVelocity().size()+a.values().pressure().size()+a.tangentialVelocity().size();
    const auto& x=a.report();const auto& y=b.report();
    check(x.usedBlockGmres==y.usedBlockGmres&&x.usedDirectFallback==y.usedDirectFallback&&x.directFallbackReason==y.directFallbackReason
        &&x.krylovIterations==y.krylovIterations&&x.krylovRestarts==y.krylovRestarts&&x.krylovSolves==y.krylovSolves&&x.linearRefinements==y.linearRefinements
        &&x.factorDimension==y.factorDimension&&x.factorNonZeros==y.factorNonZeros,"ordered SpMV changed old backend, iteration/refinement path or quotient size");
    const double before[]{x.algebraicRelativeResidual,x.freeMomentumDefect,x.maxCellContinuityDefect,x.maxPressureMeanDefect,x.maxBoundaryFluxDefect,x.maxAbsGaugeMultiplier,x.krylovRelativeResidual,x.krylovScaledRelativeResidual};
    const double after[]{y.algebraicRelativeResidual,y.freeMomentumDefect,y.maxCellContinuityDefect,y.maxPressureMeanDefect,y.maxBoundaryFluxDefect,y.maxAbsGaugeMultiplier,y.krylovRelativeResidual,y.krylovScaledRelativeResidual};
    for(std::size_t i=0;i<8;++i)check(SparseValidationPlanContract::bits(before[i])==SparseValidationPlanContract::bits(after[i]),"ordered SpMV changed an original physical/algebraic diagnostic bit");}
enum class Mode{Parallel,Nested,NoView,Unsupported,Zero};
void compareSolve(const QuadtreeMesh& mesh,const HydroDofLayout&,const Reference& ref,const Pair& solvers,const HydroSaddleRhs& rhs,const Manufactured* exact,Mode mode){
    const bool nested=mode==Mode::Nested;if(!nested)threads(1);prime();const Floating entry;
    const auto a=invoke<Result>([&]{return solvers.literal.solve(mesh,rhs);});const Floating afterA;
    if(!nested)threads(2);entry.restore();const auto b=invoke<Result>([&]{return solvers.parallel.solve(mesh,rhs);});const Floating afterB;sameFp(afterA,afterB);
    check(a.value.has_value()==b.value.has_value()&&a.kind==b.kind&&a.error==b.error,"ordered SpMV changed original public solve availability/type/first error");
    if(!a.value){check(mode==Mode::Unsupported,"ordered SpMV ordinary/nested fixture failed original public gate");++unsupportedRejectedPairs;return;}
    sameResult(*a.value,*b.value);SaddleScatterPlanContract::originalGates(*a.value);SaddleScatterPlanContract::originalGates(*b.value);
    const auto& old=a.value->report();const auto& fresh=b.value->report();
    check(!old.orderedSpmvEligible&&!old.orderedSpmvUsed&&old.orderedSpmvParallelCalls==0&&old.orderedSpmvMaximumWorkers==0,"OMP1 source factor silently acquired a parallel view");
    check(fresh.orderedSpmvParallelCalls<=fresh.orderedSpmvCalls&&fresh.orderedSpmvUsed==(fresh.orderedSpmvParallelCalls>0)
        &&fresh.orderedSpmvFallback==(fresh.orderedSpmvCalls>fresh.orderedSpmvParallelCalls),"ordered SpMV actual per-solve flags/calls disagree");
    check(old.orderedSpmvCalls==fresh.orderedSpmvCalls,"ordered SpMV changed actual main-product call count at identical original iterations");
    if(mode==Mode::Zero){check(old.orderedSpmvCalls==0&&fresh.orderedSpmvCalls==0&&!fresh.orderedSpmvUsed&&!fresh.orderedSpmvFallback,"zero RHS falsely counted as actual SpMV use");
        for(double v:b.value->values().normalVelocity())check(v==0,"zero RHS velocity not exact zero");for(double v:b.value->values().pressure())check(v==0,"zero RHS pressure not exact zero");++zeroPairs;return;}
    if(mode==Mode::Unsupported){check(fresh.orderedSpmvEligible&&fresh.orderedSpmvViewBytes>0&&!fresh.orderedSpmvUsed
        &&fresh.orderedSpmvParallelCalls==0&&fresh.orderedSpmvMaximumWorkers==0,"unsupported-current fixture changed an immutable prepared view or used parallel arithmetic");
        if(fresh.orderedSpmvCalls>0)check(fresh.orderedSpmvFallback,"unsupported actual main product did not retain Eigen fallback");
        ++unsupportedAcceptedPairs;return;}
    oldReport(*a.value);oldReport(*b.value);check(old.orderedSpmvCalls>0&&old.orderedSpmvFallback,"OMP1 literal nonzero GMRES missing actual main products");++nonzeroPairs;
    if(exact)check(relativeError(b.value->values().normalVelocity(),exact->u,exact->weights)<=1e-9&&relativeError(b.value->values().pressure(),exact->p,ref.volume)<=1e-9,"ordered SpMV lost independent manufactured physical state");
    if(mode==Mode::Parallel){check(fresh.orderedSpmvRequested&&fresh.orderedSpmvEligible&&fresh.orderedSpmvUsed&&fresh.orderedSpmvParallelCalls>0
        &&fresh.orderedSpmvMaximumWorkers==2&&fresh.orderedSpmvViewBytes>0&&fresh.orderedSpmvViewBytes<=64u*1024u*1024u,"actual two-worker nonzero GMRES did not execute bounded row SpMV");++parallelPairs;}
    else{check(!fresh.orderedSpmvUsed&&fresh.orderedSpmvParallelCalls==0&&fresh.orderedSpmvMaximumWorkers==0&&fresh.orderedSpmvFallback,"unsupported/nested/unprepared SpMV bypassed original Eigen fallback");
        if(mode==Mode::Nested){check(fresh.orderedSpmvEligible&&fresh.orderedSpmvViewBytes>0,"nested fixture lacked genuine prepared source view");++nestedPairs;}
        else if(mode==Mode::NoView){check(!fresh.orderedSpmvEligible,"late thread change fabricated a factor view");++lateThreadPairs;}else ++unsupportedAcceptedPairs;}}
// One existing XY/closed manufactured specimen, no new matrix or old-MGS oracle.
void ownCgs2(const QuadtreeMesh& mesh,const HydroPressureLoad& load,const HydroSaddleVelocityBlock& velocity,
    const HydroDofLayout& layout,const Reference& ref,const HydroSaddleRhs& rhs,const Manufactured& exact){
    threads(1);prime();const Floating entry;std::optional<Result> first;std::optional<Floating> firstFp;
    for(int count:{1,2,8}){threads(count);entry.restore();const auto solver=HydroSinglePhaseSaddle::buildClosed(load,velocity);
        const auto result=solver.solve(mesh,rhs);const Floating completed; // Before receipt/reference arithmetic.
        if(first){sameFp(*firstFp,completed);sameResult(*first,result);++orthThreadPairs;}else{first=result;firstFp=completed;}
        oldReport(result);const auto& r=result.report();
        check(r.orthogonalizationMaximumWorkers==count&&(count==1?r.orthogonalizationParallelCalls==0:r.orthogonalizationParallelCalls>0),
            "CGS2 own-request did not execute its actual requested1/2/8 team");
        check(relativeError(result.values().normalVelocity(),exact.u,exact.weights)<=1e-9
            &&relativeError(result.values().pressure(),exact.p,ref.volume)<=1e-9,"CGS2 own-request lost independent original manufactured U/p");
    }
    ++orthThreadGroups;threads(2);
}
template<class F>void rejectedPair(F action){threads(1);prime();const Floating entry;const auto a=invoke<int>([&]{action();return 1;});const Floating afterA;
    threads(2);entry.restore();const auto b=invoke<int>([&]{action();return 1;});const Floating afterB;sameFp(afterA,afterB);
    check(!a.value&&!b.value&&a.kind==b.kind&&a.error==b.error&&!a.error.empty(),"ordered SpMV changed/bypassed original public rejection");++guardPairs;}
void scene(bool rz,bool outlet){
    SparseValidationPlanContract::RestoreFloating floating;OmpRestore omp;check(std::fesetround(FE_TONEAREST)==0,"ordered SpMV cannot select RNE");_mm_setcsr(_mm_getcsr()&~((1u<<15)|(1u<<6)));unsigned native=0;
    check(_controlfp_s(&native,0,0)==0&&(native&_MCW_EM)==_MCW_EM&&(_mm_getcsr()&0x1f80u)==0x1f80u,"ordered SpMV requires masked native/SSE controls");
    threads(1);QuadtreeMesh mesh(input(4,rz?1:0));refine(mesh,1);FaceSegmentGraph graph;graph.build(mesh);
    const auto layout=HydroDofLayout::build(mesh,NativeLeafDomain::build(mesh,graph,rz));const auto integral=HydroIntegralDivergence::build(mesh,layout);const auto load=HydroPressureLoad::build(integral);
    const auto ref=reference(mesh,graph,layout,rz);check(!layout.domain().uniformLeaves(),"ordered SpMV scene lost real 2:1 leaves");std::vector<int> interior;std::vector<HydroSaddlePressureOutlet> outlets;
    for(int i=0;i<layout.velocityCount();++i){const auto& f=ref.faces[i];if(f.incidence==2)interior.push_back(i);if(outlet&&f.incidence==1&&f.horizontal&&f.y0==mesh.maxCoord.y&&f.area>0)outlets.push_back({i,2.5});}
    check(interior.size()>5&&(!outlet||outlets.size()>1),"ordered SpMV scene lacks actual interior/outlet roles");
    Matrix a=diagonal(std::vector<double>(layout.velocityCount(),4.));a.coeffRef(interior[0],interior[1])=-.125;a.coeffRef(interior[1],interior[0])=-.125;
    a.coeffRef(interior[2],interior[3])=0.;a.coeffRef(interior[3],interior[2])=-0.;a.makeCompressed();const Matrix saved=a;
    check(std::signbit(a.coeff(interior[3],interior[2]))&&a.nonZeros()==layout.velocityCount()+4,"ordered SpMV fixture discarded explicit signed zero");
    const auto velocity=HydroSaddleVelocityBlock::build(layout,a);SparseValidationPlanContract::sameMatrix(velocity.matrix(),a);
    const auto solvers=factories(load,velocity,outlets);const auto first=manufactured(mesh,layout,ref,a,solvers.literal,outlets,1.);const HydroSaddleRhs rhs(layout,first.momentum,first.known);
    compareSolve(mesh,layout,ref,solvers,rhs,&first,Mode::Parallel);
    const auto second=manufactured(mesh,layout,ref,a,solvers.literal,outlets,.5);compareSolve(mesh,layout,ref,solvers,HydroSaddleRhs(layout,second.momentum,second.known),&second,Mode::Parallel);
    if(!rz&&!outlet){
        ownCgs2(mesh,load,velocity,layout,ref,rhs,first);
        compareSolve(mesh,layout,ref,Pair{solvers.literal,solvers.literal},rhs,&first,Mode::NoView);
        threads(2);std::exception_ptr failure;
#pragma omp parallel num_threads(2)
        {if(omp_get_thread_num()==0){try{compareSolve(mesh,layout,ref,solvers,rhs,&first,Mode::Nested);}catch(...){failure=std::current_exception();}}}
        if(failure)std::rethrow_exception(failure);
        const Floating nearest;check(std::fesetround(FE_UPWARD)==0,"ordered SpMV cannot select directed mode");compareSolve(mesh,layout,ref,solvers,rhs,nullptr,Mode::Unsupported);nearest.restore();
        for(unsigned mode:{1u<<15,1u<<6}){_mm_setcsr(_mm_getcsr()|mode);compareSolve(mesh,layout,ref,solvers,rhs,nullptr,Mode::Unsupported);nearest.restore();}
        const std::vector<double> zero(layout.velocityCount());compareSolve(mesh,layout,ref,solvers,HydroSaddleRhs(layout,zero,zero),nullptr,Mode::Zero);
        const auto fresh=HydroDofLayout::build(mesh,NativeLeafDomain::build(mesh,graph,false));rejectedPair([&]{solvers.parallel.solve(mesh,HydroSaddleRhs(fresh,first.momentum,first.known));});
        for(double bad:{std::numeric_limits<double>::infinity(),std::numeric_limits<double>::quiet_NaN()})rejectedPair([&]{auto invalid=a;invalid.coeffRef(0,0)=bad;HydroSaddleVelocityBlock::build(layout,invalid);});
        rejectedPair([&]{auto invalid=first.momentum;invalid[interior[0]]=std::numeric_limits<double>::quiet_NaN();HydroSaddleRhs(layout,invalid,first.known);});
        auto& v=mesh.verts(mesh.vertIndexToLabel[0]);const auto coordinate=v.coord();v.setCoord(coordinate+Vector2D{.01,0});
        try{rejectedPair([&]{solvers.parallel.solve(mesh,rhs);});}catch(...){v.setCoord(coordinate);throw;}v.setCoord(coordinate);
        compareSolve(mesh,layout,ref,solvers,rhs,&first,Mode::Parallel);
    }
    SparseValidationPlanContract::sameMatrix(a,saved);++scenes;
    std::cout<<"ORDERED_ROW_SPMV_SCENE,rz,"<<rz<<",outlet,"<<outlet<<",leaves,"<<layout.pressureCount()<<",nonzero_original_GMRES_and_dense_oracle,PASS\n";
}
void run(){scene(false,false);scene(false,true);scene(true,false);scene(true,true);
    check(scenes==4&&parallelPairs==9&&denseOracles==8&&zeroPairs==1&&nestedPairs==1&&lateThreadPairs==1&&unsupportedAcceptedPairs+unsupportedRejectedPairs==3
        &&guardPairs==5&&nonzeroPairs>=11&&fpPairs>nonzeroPairs&&matrixValues>0&&fieldValues>0&&orthThreadGroups==1&&orthThreadPairs==2,"ordered SpMV actual parallel/literal/FP/oracle/negative coverage incomplete");
    std::cout<<"ORDERED_ROW_SPMV_PASS,scenes,"<<scenes<<",nonzero_GMRES_pairs,"<<nonzeroPairs<<",actual_two_worker_pairs,"<<parallelPairs
        <<",independent_dense_manufactured_oracles,"<<denseOracles<<",matrix_coefficient_bits,"<<matrixValues<<",physical_field_bits,"<<fieldValues<<",full_native_rawCSR_errno_pairs,"<<fpPairs
        <<",nested_original_Eigen_fallback_pairs,"<<nestedPairs<<",late_threads_no_new_view_pairs,"<<lateThreadPairs<<",unsupported_original_accept_pairs,"<<unsupportedAcceptedPairs
        <<",unsupported_original_rejection_pairs,"<<unsupportedRejectedPairs<<",zero_RHS_not_positive_coverage_pairs,"<<zeroPairs<<",original_public_rejection_pairs,"<<guardPairs
        <<",CGS2_own_1_2_8_groups,"<<orthThreadGroups<<",CGS2_own_FP_pairs,"<<orthThreadPairs<<",private_quotient_capacity_worker_fault_not_injected,no_speed_claim\n";}
#else
void run(){std::cout<<"ORDERED_ROW_SPMV_UNAVAILABLE,unsupported_compile_Eigen_or_platform\n";}
#endif
}

namespace LazyDirectInputContract {
#if defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE) \
    && !defined(__FAST_MATH__) && (!defined(_M_FP_FAST) || !_M_FP_FAST) \
    && EIGEN_WORLD_VERSION==3 && EIGEN_MAJOR_VERSION==4 && EIGEN_MINOR_VERSION==1 \
    && !defined(EIGEN_SPARSEMATRIX_PLUGIN) && !defined(EIGEN_SPARSEMATRIXBASE_PLUGIN) \
    && !defined(EIGEN_SPARSE_CREATE_TEMPORARY_PLUGIN) \
    && !defined(EIGEN_SPARSE_TRANSPOSED_COPY_PLUGIN) \
    && !defined(EIGEN_SPARSE_COMPRESSED_STORAGE_REALLOCATE_PLUGIN) \
    && !defined(EIGEN_SCALAR_BINARY_OP_PLUGIN)
using Floating=SparseValidationPlanContract::Floating;
std::size_t scenes=0,gmres=0,direct=0,copyFpPairs=0,factoryFpPairs=0,solveFpPairs=0,typedBorrowCases=0;
void prime(){check(std::feclearexcept(FE_ALL_EXCEPT)==0,"lazy Direct cannot clear real flags");
    _mm_setcsr(_mm_getcsr()&~0x3fu);check(std::feraiseexcept(FE_INEXACT)==0,"lazy Direct cannot prime passive timer flag");errno=E2BIG;}
void copiedInputIsNeutral(const HydroSinglePhaseSaddle::Matrix& matrix){
    // This actual one-cell outlet has no quotient gauge elimination: its full
    // public matrix IS the Direct input. Independently execute the old Eigen
    // index conversion and stable IEEE-zero compaction, never the new helper.
    const Floating before;
    {
        using D=Eigen::SparseMatrix<double,Eigen::ColMajor,int>;
        D converted=matrix;converted.makeCompressed();
        check(converted.rows()==matrix.rows()&&converted.cols()==matrix.cols()&&converted.nonZeros()==matrix.nonZeros(),"lazy Direct old conversion dimensions changed");
        for(Eigen::Index j=0;j<=matrix.outerSize();++j)check(converted.outerIndexPtr()[j]==matrix.outerIndexPtr()[j],"lazy Direct old conversion outer changed");
        for(Eigen::Index k=0;k<matrix.nonZeros();++k)check(converted.innerIndexPtr()[k]==matrix.innerIndexPtr()[k]
            &&std::memcmp(converted.valuePtr()+k,matrix.valuePtr()+k,sizeof(double))==0,"lazy Direct old conversion entry changed");
        const auto nonzero=[](const double& x){std::uint64_t word=0;std::memcpy(&word,&x,sizeof(word));return (word&UINT64_C(0x7fffffffffffffff))!=0;};
        converted.prune([&](Eigen::Index,Eigen::Index,const double& x){return nonzero(x);});
        for(Eigen::Index col=0;col<matrix.outerSize();++col){int p=converted.outerIndexPtr()[col];
            for(Eigen::Index k=matrix.outerIndexPtr()[col];k<matrix.outerIndexPtr()[col+1];++k)if(nonzero(matrix.valuePtr()[k])){
                check(p<converted.outerIndexPtr()[col+1]&&converted.innerIndexPtr()[p]==matrix.innerIndexPtr()[k]
                    &&std::memcmp(converted.valuePtr()+p,matrix.valuePtr()+k,sizeof(double))==0,"lazy Direct old stable prune changed current nonzero bytes");++p;}
            check(p==converted.outerIndexPtr()[col+1],"lazy Direct old stable prune inserted entries");}
    }
    const Floating after;SparseValidationPlanContract::sameFlags(before,after);++copyFpPairs;
}
void scene(bool rz,bool forceDirect){
    SparseValidationPlanContract::RestoreFloating restore;
#ifdef _OPENMP
    const int oldThreads=omp_get_max_threads();struct Threads{int n;~Threads(){omp_set_num_threads(n);}} threads{oldThreads};omp_set_num_threads(1);
#endif
    check(std::fegetround()==FE_TONEAREST&&(_mm_getcsr()&0x8040u)==0,"lazy Direct fixture needs original nearest/gradual execution");
    QuadtreeMesh mesh(input(1,rz?2:0));const auto graph=selected(mesh,0);
    const auto layout=HydroDofLayout::build(mesh,NativeLeafDomain::build(mesh,graph,rz));
    const auto load=HydroPressureLoad::build(HydroIntegralDivergence::build(mesh,layout));const auto ref=reference(mesh,graph,layout,rz);
    const int nu=layout.velocityCount();int top=-1;std::vector<double> u(nu),known(nu),p(1,.5),diagonalValues(nu,1.);
    check(layout.pressureCount()==1,"lazy Direct fixture lost actual one-cell geometry");
    for(int i=0;i<nu;++i){const auto& face=ref.faces[i];if(face.horizontal)u[i]=.125;
        if(face.horizontal&&face.y0==mesh.maxCoord.y){check(top<0,"lazy Direct fixture has multiple top rows");top=i;}else known[i]=u[i];}
    check(top>=0&&ref.faces[top].incidence==1&&ref.faces[top].area>0,"lazy Direct fixture lacks actual pressure outlet");
    // Other known-face diagonal entries keep the old full-block norm resolved;
    // only the real free outlet coefficient is subnormal. Its reciprocal
    // overflows, invoking the existing preconditioner-ineligible Direct route.
    const double d=forceDirect?std::ldexp(1.,-1024):1.;diagonalValues[top]=d;
    const Matrix a=diagonal(diagonalValues);const auto velocity=HydroSaddleVelocityBlock::build(layout,a);
    const std::vector<HydroSaddlePressureOutlet> outlets{{top,0.}};
    prime();const Floating entry;
    const auto first=HydroSinglePhaseSaddle::buildWithPressureOutlets(load,velocity,outlets);const Floating afterFirst;
    entry.restore();const auto second=HydroSinglePhaseSaddle::buildWithPressureOutlets(load,velocity,outlets);const Floating afterSecond;
    SparseValidationPlanContract::sameFlags(afterFirst,afterSecond);++factoryFpPairs;
    const auto owner=HydroMomentumGeometry::build(mesh,HydroPhysicalDivergence::build(mesh,load.divergence()));
    entry.restore();const auto typed=HydroSinglePhaseSaddle::buildWithPressureOutlets(load,velocity,outlets,
        HydroSaddleConfig{},HydroSaddlePreparationTimings{},owner);const Floating typedFactoryFp;
    SparseValidationPlanContract::sameFlags(afterFirst,typedFactoryFp);
    SaddleScatterPlanContract::sameMatrix(first.matrix(),typed.matrix());const auto savedTyped=typed.matrix();
    SaddleScatterPlanContract::sameMatrix(first.matrix(),second.matrix());
    check(first.pressureGaugeCount()==0&&first.freeVelocityDofs()==std::vector<int>{top}&&first.matrix().rows()==2,"lazy Direct fixture is not the true unbordered two-equation quotient");
    prime();copiedInputIsNeutral(first.matrix());
    auto momentum=multiply(a,u);for(const auto& term:ref.terms)momentum[term.u]-=term.outwardArea*p[term.p];
    const HydroSaddleRhs rhs(layout,momentum,known);const auto knownMomentum=multiply(a,known);
    Eigen::Vector2d b;b[0]=momentum[top]-knownMomentum[top];b[1]=0.;for(const auto& term:ref.terms)b[1]+=term.outwardArea*known[term.u];
    Eigen::Matrix2d independent;independent<<d,-ref.faces[top].area,-ref.faces[top].area,0.;
    check((Eigen::MatrixXd(first.matrix())-independent).norm()<=1e-12*independent.norm()&&b.stableNorm()>0,"lazy Direct independent physical matrix/RHS mismatch");
    const Eigen::FullPivLU<Eigen::Matrix2d> oracle(independent);check(oracle.isInvertible(),"lazy Direct independent physical system singular");
    const Eigen::Vector2d exact=oracle.solve(b);check(exact.allFinite()&&(independent*exact-b).stableNorm()<=1e-12*b.stableNorm(),"lazy Direct independent original residual failed");
    prime();const Floating solveEntry;const auto one=first.solve(mesh,rhs);const Floating oneFp;
    solveEntry.restore();const auto two=second.solve(mesh,rhs);const Floating twoFp;
    SparseValidationPlanContract::sameFlags(oneFp,twoFp);++solveFpPairs;
    SaddleScatterPlanContract::sameValues(one.values().normalVelocity(),two.values().normalVelocity());
    SaddleScatterPlanContract::sameValues(one.values().pressure(),two.values().pressure());
    SaddleScatterPlanContract::originalGates(one);SaddleScatterPlanContract::originalGates(two);
    prime();const auto typedResult=typed.solve(mesh,rhs);
    SaddleScatterPlanContract::sameMatrix(savedTyped,typed.matrix());
    SaddleScatterPlanContract::sameValues(one.values().normalVelocity(),typedResult.values().normalVelocity());
    SaddleScatterPlanContract::sameValues(one.values().pressure(),typedResult.values().pressure());
    SaddleScatterPlanContract::originalGates(typedResult);
    const auto& borrowed=typedResult.report().sparseAssemblyPlan;
    check(borrowed.quotientAliasesSystem&&!borrowed.quotientScatterUsed&&typedResult.report().factorDimension==typed.matrix().rows()
        &&typedResult.report().factorNonZeros==typed.matrix().nonZeros(),"typed no-gauge factor did not borrow the unchanged current system");++typedBorrowCases;
    for(const auto* result:{&one,&two,&typedResult}){const auto& r=result->report();
        check(r.directInputDeferred,"lazy Direct actual pre-certificate did not omit eager input");
        check(r.directInputMaterialized==forceDirect,"lazy Direct actual materialization branch differs from requested physical witness");
        check(r.factorDimension==2&&std::abs(result->values().normalVelocity()[top]-exact[0])<=1e-12
            &&std::abs(result->values().pressure()[0]-exact[1])<=1e-12,"lazy Direct changed original two-equation solution");
        if(forceDirect)check(r.usedSuperLU&&r.usedDirectFallback&&r.directFallbackReason==1&&!r.usedBlockGmres,"lazy Direct witness did not actually execute original Direct route");
        else check(r.usedBlockGmres&&!r.usedSuperLU&&!r.usedDirectFallback&&r.krylovIterations>0,"lazy Direct ordinary witness did not actually execute nonzero GMRES");}
    if(forceDirect)++direct;else ++gmres;++scenes;
}
void run(){for(bool rz:{false,true})for(bool force:{false,true})scene(rz,force);
    check(scenes==4&&gmres==2&&direct==2&&copyFpPairs==4&&factoryFpPairs==4&&solveFpPairs==4&&typedBorrowCases==4,"lazy Direct finite branch coverage incomplete");
    std::cout<<"LAZY_DIRECT_INPUT_PASS,scenes,"<<scenes<<",actual_nonzero_gmres,"<<gmres<<",actual_nonzero_direct,"<<direct
        <<",old_Eigen_copy_prune_FP_errno_pairs,"<<copyFpPairs<<",factory_FP_pairs,"<<factoryFpPairs<<",solve_FP_pairs,"<<solveFpPairs
        <<",actual_typed_system_borrow_cases,"<<typedBorrowCases<<",allocator_failure_order_not_promised,no_speed_claim\n";}
#else
void run(){std::cout<<"LAZY_DIRECT_INPUT_UNAVAILABLE,unsupported_compile_Eigen_or_platform\n";}
#endif
}

namespace ParallelBlockGmresContract {
void run(){
#ifdef _OPENMP
    struct ThreadsRestore {
        int maximum=omp_get_max_threads(),dynamic=omp_get_dynamic();
        ~ThreadsRestore(){omp_set_num_threads(maximum);omp_set_dynamic(dynamic);}
    } threadsRestore;
    omp_set_dynamic(0);omp_set_num_threads(2);
#endif
    constexpr int size=128,blockRows=16;
    std::vector<Eigen::Triplet<double>> entries;
    Eigen::VectorXd exact(size);
    for(int i=0;i<size;++i){
        exact[i]=1+.01*i;
        entries.emplace_back(i,i,4.);
        if(i>0)entries.emplace_back(i,i-1,-1.);
        if(i+1<size)entries.emplace_back(i,i+1,-.25);
    }
    LampLinear::Sparse matrix(size,size);matrix.setFromTriplets(entries.begin(),entries.end());matrix.makeCompressed();
    const Eigen::VectorXd rhs=matrix*exact;
    LampLinear::ParallelBlockJacobi preconditioner;
    check(preconditioner.compute(matrix,blockRows)&&preconditioner.blockCount()==size/blockRows,
        "parallel block Jacobi did not prepare all independent blocks");
    LampLinear::GmresReport gmres;
    Eigen::VectorXd x;
    check(LampLinear::solveGmres(matrix,preconditioner,rhs,x,1e-12,100,32,gmres),
        "parallel block GMRES failed the nonsymmetric cross-block system");
    const double trueResidual=(rhs-matrix*x).stableNorm()/rhs.stableNorm();
    const double error=(x-exact).stableNorm()/exact.stableNorm();
    check(gmres.iterations>0&&std::isfinite(trueResidual)&&trueResidual<=1e-12&&error<=1e-11,
        "parallel block GMRES did not recover the independent known solution/original residual");
#ifdef _OPENMP
    check(preconditioner.maximumWorkers()>=2&&gmres.maximumWorkers>=2,
        "parallel block GMRES did not execute multiple workers");
#endif
    const int iterations=gmres.iterations,workers=gmres.maximumWorkers;
    const Eigen::VectorXd zero=Eigen::VectorXd::Zero(size);
    check(LampLinear::solveGmres(matrix,preconditioner,zero,x,1e-12,100,32,gmres)
        &&x.isZero(0)&&gmres.iterations==0&&gmres.relativeResidual==0,
        "parallel block GMRES rejected or fabricated work for a zero RHS");
    LampLinear::Sparse singular(size,size);
    for(int i=0;i<size;++i)if(i!=17)singular.insert(i,i)=1.;
    singular.makeCompressed();
    LampLinear::ParallelBlockJacobi rejectedFactor;
    check(!rejectedFactor.compute(singular,blockRows),"parallel block Jacobi accepted a singular local block");
    check(!rejectedFactor.solve(rhs,x),"failed block Jacobi preparation remained usable");
    // Reuse the existing independently manufactured U/p fixture to exercise
    // the production selection, original physical gates and no-direct policy.
    algebra(4,0,0,0,false,false,true);
    std::cout<<"PARALLEL_BLOCK_GMRES_PASS,rows,"<<size<<",blocks,"<<size/blockRows
        <<",iterations,"<<iterations<<",workers,"<<workers<<",true_residual,"<<trueResidual
        <<",solution_error,"<<error<<",zero_rhs,1,singular_block_rejected,1,production_saddle,1\n";
}
}


namespace PressureCoarseContract
{
using Coarse=LAMP::Hydrodynamics::Detail::HydroPressureCoarseCorrection;
using Sparse=Coarse::Sparse;
void run(){
    struct EnvironmentRestore{
        std::fenv_t state{};int error=errno;
#if defined(_M_X64) || defined(__SSE2__)
        unsigned csr=_mm_getcsr();
#endif
        EnvironmentRestore(){check(std::fegetenv(&state)==0,"coarse fixture cannot retain environment");}
        ~EnvironmentRestore(){if(std::fesetenv(&state)!=0)std::terminate();
#if defined(_M_X64) || defined(__SSE2__)
            _mm_setcsr(csr);
#endif
            errno=error;}
    } environmentRestore;
    check(std::fesetround(FE_TONEAREST)==0,"coarse fixture requires RNE");
#if defined(_M_X64) || defined(__SSE2__)
    _mm_setcsr(_mm_getcsr()&~0x8040u);
#endif
#ifdef _OPENMP
    struct RestoreThreads{int maximum=omp_get_max_threads(),dynamic=omp_get_dynamic();
        ~RestoreThreads(){omp_set_dynamic(dynamic);omp_set_num_threads(maximum);}} restoreThreads;
    omp_set_dynamic(0);omp_set_num_threads(2);
#endif
    // The same central refine boxes and open-top RZ MMS as the original ordered
    // SpMV witness, enlarged just enough to contain multiple actual core256s.
    QuadtreeMesh mesh(input(16,1));refine(mesh,2);FaceSegmentGraph graph;graph.build(mesh);
    const auto layout=HydroDofLayout::build(mesh,NativeLeafDomain::build(mesh,graph,true));
    const auto load=HydroPressureLoad::build(HydroIntegralDivergence::build(mesh,layout));
    const auto ref=reference(mesh,graph,layout,true);
    check(!layout.domain().uniformLeaves()&&layout.pressureCount()>512,
        "pressure coarse fixture lacks actual RZ2:1/multiple256 cores");
    std::vector<HydroSaddlePressureOutlet> outlets;
    for(int i=0;i<layout.velocityCount();++i){const auto& f=ref.faces[i];
        if(f.incidence==1&&f.horizontal&&f.y0==mesh.maxCoord.y&&f.area>0)outlets.push_back({i,2.5});}
    check(outlets.size()>1,"pressure coarse fixture lost its physical outlet");
    Matrix a=diagonal(std::vector<double>(layout.velocityCount(),4.));
    const auto velocity=HydroSaddleVelocityBlock::build(layout,a);
    const auto original=HydroSinglePhaseSaddle::buildWithPressureOutlets(load,velocity,outlets);
    HydroSaddleConfig config;check(!config.gmresGlobalPressure,"global pressure factor ceased to be an explicit opt-in");
    config.gmresOnlyParallel=true;
    const auto solver=HydroSinglePhaseSaddle::buildWithPressureOutlets(load,velocity,outlets,config);
    SaddleScatterPlanContract::sameMatrix(original.matrix(),solver.matrix());
    check(original.freeVelocityDofs()==solver.freeVelocityDofs()&&solver.pressureGaugeCount()==0,
        "pressure coarse changed original open quotient/partition/gauge");
    const Matrix savedK=solver.matrix();const auto savedA=a;
    const int nu=layout.velocityCount(),np=layout.pressureCount(),nf=static_cast<int>(solver.freeVelocityDofs().size());
    std::vector<double> u(nu),known(nu),p(np),weights(nu,1.),traction(nu);
    std::vector<bool> isOutlet(nu);for(const auto& o:outlets){isOutlet[o.velocityDof]=true;traction[o.velocityDof]=-ref.faces[o.velocityDof].area*o.pressure;}
    for(int i=0;i<nu;++i){const auto& f=ref.faces[i];
        const double flow=f.horizontal?polynomialPsi(f.x0,f.y0,mesh)-polynomialPsi(f.x1,f.y0,mesh)
            :polynomialPsi(f.x0,f.y1,mesh)-polynomialPsi(f.x0,f.y0,mesh);
        u[i]=f.area==0?0:2*pi*flow/f.area;if(f.horizontal)u[i]+=.125;
        if(f.incidence==1&&!isOutlet[i])known[i]=u[i];}
    for(int i=0;i<np;++i){const auto c=ref.center[i];p[i]=1+.3*c.x+.7*c.y+.4*c.x*c.y;}
    auto momentum=multiply(a,u);for(const auto& t:ref.terms)momentum[t.u]-=t.outwardArea*p[t.p];
    for(int i=0;i<nu;++i)momentum[i]-=traction[i];
    const HydroSaddleRhs rhs(layout,momentum,known);
    const auto result=solver.solve(mesh,rhs);report(result);SaddleScatterPlanContract::originalGates(result);
    const auto& diagnostics=result.report();
    check(diagnostics.usedBlockGmres&&!diagnostics.usedDirectFallback&&!diagnostics.usedSuperLU
        &&diagnostics.pressureCoarsePrepared&&diagnostics.pressureCoarseApplications>0
        &&diagnostics.pressureCoarseDimension==(np+255)/256,
        "production pressure coarse did not execute its qualified same-core policy");
    check(!diagnostics.pressureGlobalRequested&&!diagnostics.pressureGlobalPrepared
        &&diagnostics.pressureGlobalSolveCalls==0,
        "default coarse policy performed an unrequested global pressure action");
    const auto eu=relativeError(result.values().normalVelocity(),u,weights);
    const auto ep=relativeError(result.values().pressure(),p,ref.volume);
    check(eu<=1e-9&&ep<=1e-9,"pressure coarse lost original independent manufactured U/p gates");
    // In this anchored open fixture the original pressure quotient retains
    // every native pressure in the same order. Rebuild S with the production
    // operation order from the SAME public K; never use C^T A^-1 C instead.
    Sparse coupling=solver.matrix().block(0,nf,nf,np),weighted=coupling;
    for(Eigen::Index col=0;col<weighted.outerSize();++col)
        for(Sparse::InnerIterator e(weighted,col);e;++e)e.valueRef()*=1/solver.matrix().coeff(e.row(),e.row());
    Sparse s=coupling.transpose()*weighted;s.makeCompressed();
    // The explicit global policy consumes the SAME public K, source RHS and
    // current S. Its success still requires the original independent U/p gates.
    HydroSaddleConfig invalidGlobal;invalidGlobal.gmresGlobalPressure=true;
    bool invalidGlobalRejected=false;
    try{HydroSinglePhaseSaddle::buildWithPressureOutlets(load,velocity,outlets,invalidGlobal);}
    catch(const std::invalid_argument&){invalidGlobalRejected=true;}
    check(invalidGlobalRejected,"global pressure policy accepted a non-GMRES-only configuration");
    auto globalConfig=config;globalConfig.gmresGlobalPressure=true;
    const auto globalSolver=HydroSinglePhaseSaddle::buildWithPressureOutlets(load,velocity,outlets,globalConfig);
    SaddleScatterPlanContract::sameMatrix(globalSolver.matrix(),savedK);
    check(globalSolver.freeVelocityDofs()==solver.freeVelocityDofs()
        &&globalSolver.pressureGaugeCount()==solver.pressureGaugeCount()
        &&globalSolver.pressureOutletLoads()==solver.pressureOutletLoads(),
        "global pressure policy changed original open quotient/partition/traction");
    const auto globalCheck=[&](const HydroSaddleResult& value){
        report(value);SaddleScatterPlanContract::originalGates(value);const auto& r=value.report();
        check(r.usedBlockGmres&&!r.usedSuperLU&&!r.usedDirectFallback&&r.krylovIterations>0
            &&r.pressureGlobalRequested&&r.pressureGlobalAvailable&&r.pressureGlobalPrepared
            &&r.pressureGlobalTimingsAvailable&&r.pressureGlobalSolveCalls>0,
            "global pressure policy did not consume its qualified actual factor without direct fallback");
        check(r.parallelPressureBlockCount==0&&r.parallelPressureMaximumWorkers==0
            &&!r.pressureCoarsePrepared&&r.pressureCoarseApplications==0,
            "global pressure policy also executed the original RAS/coarse action");
        constexpr std::size_t limit=128u*1024u*1024u;
        check(r.pressureGlobalRows==static_cast<std::size_t>(np)
            &&r.pressureGlobalNonZeros==static_cast<std::size_t>(s.nonZeros())
            &&r.pressureGlobalFactorNonZeros>0&&r.pressureGlobalStorageBytes>0
            &&r.pressureGlobalStorageBytes<=r.pressureGlobalAdmissionBytes
            &&r.pressureGlobalAdmissionBytes<=r.pressureGlobalByteLimit&&r.pressureGlobalByteLimit==limit,
            "global pressure factor stock/admission did not describe the bounded actual S factor");
        for(double seconds:{r.pressureGlobalSymbolicSeconds,r.pressureGlobalNumericSeconds,r.pressureGlobalSolveSeconds})
            check(std::isfinite(seconds)&&seconds>=0,"global pressure factor supplied invalid available timings");
        check(relativeError(value.values().normalVelocity(),u,weights)<=1e-9
            &&relativeError(value.values().pressure(),p,ref.volume)<=1e-9,
            "global pressure policy lost original independent manufactured U/p gates");
    };
    const auto globalResult=globalSolver.solve(mesh,rhs);globalCheck(globalResult);
    const auto globalAgain=globalSolver.solve(mesh,rhs);globalCheck(globalAgain);
    const auto& globalDiagnostics=globalResult.report();const auto& repeatedGlobalDiagnostics=globalAgain.report();
    check(globalDiagnostics.pressureGlobalRows==repeatedGlobalDiagnostics.pressureGlobalRows
        &&globalDiagnostics.pressureGlobalNonZeros==repeatedGlobalDiagnostics.pressureGlobalNonZeros
        &&globalDiagnostics.pressureGlobalFactorNonZeros==repeatedGlobalDiagnostics.pressureGlobalFactorNonZeros
        &&globalDiagnostics.pressureGlobalStorageBytes==repeatedGlobalDiagnostics.pressureGlobalStorageBytes
        &&globalDiagnostics.pressureGlobalAdmissionBytes==repeatedGlobalDiagnostics.pressureGlobalAdmissionBytes,
        "repeated same-owner RHS changed retained global factor stock");
    check(rhs.integralMomentum()==momentum&&rhs.prescribedNormalVelocity()==known,
        "global pressure action mutated the original shared source RHS");
    SaddleScatterPlanContract::sameMatrix(globalSolver.matrix(),savedK);
    LampLinear::ParallelBlockJacobi ras;check(ras.compute(s,256,true),"original core256 RAS fixture failed preparation");
    Coarse coarse;check(coarse.compute(s)&&coarse.ready(),"same actual S coarse qualification failed");
    check(coarse.report().dimension==ras.blockCount()&&coarse.coreRows()==256,
        "coarse Z did not use original non-overlap core256 mapping");
    const auto savedS=s;const auto savedG=coarse.coarseMatrix();
    const auto local=[&](const Eigen::Ref<const Eigen::VectorXd>& in,Eigen::Ref<Eigen::VectorXd> out){return ras.solve(in,out);};
    const int k=coarse.report().dimension;
    Eigen::MatrixXd z=Eigen::MatrixXd::Zero(np,k);for(int i=0;i<np;++i)z(i,i/256)=1;
    const Eigen::MatrixXd independentG=z.transpose()*Eigen::MatrixXd(s)*z;
    const double eps=std::numeric_limits<double>::epsilon(),gate=std::sqrt(eps);
    check((savedG-independentG).lpNorm<Eigen::Infinity>()<=gate*std::max(1.,independentG.lpNorm<Eigen::Infinity>()),
        "retained Galerkin G differs from independent actual S/core aggregation");
    Eigen::VectorXd fine(np),out(np),sr(np);double coarsePeak=0;
    for(int j=0;j<k;++j){fine=z.col(j);sr=s*fine;check(coarse.apply(sr,out,local),"actual coarse mode application failed");
        const double defect=(out-fine).lpNorm<Eigen::Infinity>()/std::max(1.,out.lpNorm<Eigen::Infinity>());
        check(defect<=gate,"actual PSZ did not preserve core mode");coarsePeak=std::max(coarsePeak,defect);}
    for(int i=0;i<np;++i)fine[i]=(static_cast<double>(i%17)-8)/17;
    check(coarse.apply(fine,out,local),"actual general RHS coarse application failed");
    const Eigen::VectorXd residual=fine-s*out;
    const double scale=std::max({1.,(z.transpose()*fine).lpNorm<Eigen::Infinity>(),(z.transpose()*(s*out)).lpNorm<Eigen::Infinity>()});
    const double leftDefect=(z.transpose()*residual).lpNorm<Eigen::Infinity>()/scale;
    check(leftDefect<=gate,"actual ZtSP identity failed on general RHS");
    // Independent explicit small-G solves evaluate the two distinct projectors.
    // No full Q/E/F is built in production.
    const Eigen::FullPivLU<Eigen::MatrixXd> independentFactor(independentG);
    check(independentFactor.rank()==k,"independent actual Galerkin matrix is singular");
    const auto q=[&](const Eigen::VectorXd& v)->Eigen::VectorXd{
        const Eigen::VectorXd c=z.transpose()*v;return z*independentFactor.solve(c);};
    const auto e=[&](const Eigen::VectorXd& v)->Eigen::VectorXd{const Eigen::VectorXd sv=s*v;return v-q(sv);};
    const auto f=[&](const Eigen::VectorXd& v)->Eigen::VectorXd{return v-s*q(v);};
    Eigen::VectorXd ps(np),sp(np),b(np);
    sr=s*fine;check(coarse.apply(sr,ps,local),"actual PS general application failed");
    const auto er=e(fine);sr=s*er;check(local(sr,b),"original RAS E identity RHS failed");
    const auto rightE=e(Eigen::VectorXd(er-b));
    const double rightDefect=((fine-ps)-rightE).lpNorm<Eigen::Infinity>()/
        std::max({1.,fine.lpNorm<Eigen::Infinity>(),ps.lpNorm<Eigen::Infinity>(),rightE.lpNorm<Eigen::Infinity>()});
    const auto fr=f(fine);check(local(fr,b),"original RAS F identity RHS failed");
    const auto rightF=f(Eigen::VectorXd(fr-s*b));sp=s*out;
    const double residualDefect=((fine-sp)-rightF).lpNorm<Eigen::Infinity>()/
        std::max({1.,fine.lpNorm<Eigen::Infinity>(),sp.lpNorm<Eigen::Infinity>(),rightF.lpNorm<Eigen::Infinity>()});
    check(rightDefect<=gate&&residualDefect<=gate,"actual same-S left/right balanced factorizations failed");
    // Retained factors/actions are immutable; failed applies never change the
    // policy to original RAS. Distinct callbacks reveal accidental fallback.
    int unexpected=0;auto invalid=fine;invalid[0]=std::numeric_limits<double>::quiet_NaN();
    check(!coarse.apply(invalid,ps,[&](const auto&,auto){++unexpected;return true;})&&unexpected==0&&coarse.ready(),
        "nonfinite coarse RHS switched policy or executed a local solve");
    check(!coarse.apply(fine,ps,[&](const auto&,auto){++unexpected;return false;})&&unexpected==1&&coarse.ready(),
        "failed active coarse apply silently switched to original RAS");
    check(coarse.apply(fine,ps,local)&&(ps.array()==out.array()).all(),
        "failed apply corrupted retained factor/scratch or deterministic RHS");
    {
        EnvironmentRestore stickyRestore;
        check(std::feclearexcept(FE_ALL_EXCEPT)==0&&std::feraiseexcept(FE_INEXACT)==0,
            "coarse fixture cannot prime native flags");errno=E2BIG;
        const int flags=std::fetestexcept(FE_ALL_EXCEPT),rounding=std::fegetround();
#if defined(_M_X64) || defined(__SSE2__)
        const unsigned csr=_mm_getcsr();
#endif
        check(coarse.apply(fine,ps,local),"primed environment coarse action failed");
        check(errno==E2BIG&&std::fetestexcept(FE_ALL_EXCEPT)==flags&&std::fegetround()==rounding,
            "coarse action leaked FP flags/rounding/errno");
#if defined(_M_X64) || defined(__SSE2__)
        check(_mm_getcsr()==csr,"coarse action leaked MXCSR");
#endif
        check(std::fesetround(FE_UPWARD)==0,"coarse fixture cannot select unsupported rounding");
        check(!coarse.apply(fine,ps,local)&&coarse.ready(),
            "active coarse policy changed under unsupported rounding");
        Coarse unsupported;
        check(!unsupported.compute(s)&&unsupported.report().reason==Coarse::Reason::UnsupportedEnvironment,
            "coarse setup accepted unsupported rounding");
    }
    SaddleScatterPlanContract::sameMatrix(s,savedS);
    check((coarse.coarseMatrix().array()==savedG.array()).all(),"coarse applications mutated original G");
    auto nonfinite=s;nonfinite.coeffRef(0,0)=std::numeric_limits<double>::quiet_NaN();
    Coarse rejected;check(!rejected.compute(nonfinite)&&!rejected.ready()
        &&rejected.report().reason==Coarse::Reason::InvalidMatrix,"nonfinite S coarse setup was accepted");
    Sparse overflowG(np,np);
    for(int i=0;i<np;++i)overflowG.insert(i,i)=std::numeric_limits<double>::max();
    overflowG.makeCompressed();
    check(!rejected.compute(overflowG)&&!rejected.ready()&&rejected.report().reason==Coarse::Reason::InvalidMatrix,
        "finite S producing nonfinite Galerkin G was accepted");
    Sparse singular(np,np);singular.makeCompressed();
    check(!rejected.compute(singular)&&!rejected.ready()&&rejected.report().reason==Coarse::Reason::Rank,
        "singular Galerkin G was accepted or repaired");
    auto mutated=s;Coarse retained;check(retained.compute(mutated),"retained coarse input qualification failed");
    mutated.coeffRef(0,0)+=1;
    check(retained.apply(fine,ps,local)&&(ps.array()==out.array()).all(),
        "coarse helper borrowed a mutable source matrix");
    SaddleScatterPlanContract::sameMatrix(a,savedA);SaddleScatterPlanContract::sameMatrix(solver.matrix(),savedK);
    std::cout<<"PRESSURE_COARSE_PASS,rz_real2to1_open,1,pressure_rows,"<<np<<",actual256cores,"<<k
        <<",coarse_mode_defect,"<<coarsePeak<<",ZtSP_general_defect,"<<leftDefect
        <<",I_PS_factor_defect,"<<rightDefect<<",I_SP_factor_defect,"<<residualDefect
        <<",U_error,"<<eu<<",p_error,"<<ep<<",iterations,"<<diagnostics.krylovIterations
        <<",coarse_applications,"<<diagnostics.pressureCoarseApplications
        <<",global_rows,"<<globalDiagnostics.pressureGlobalRows
        <<",global_S_nnz,"<<globalDiagnostics.pressureGlobalNonZeros
        <<",global_strict_L_nnz,"<<globalDiagnostics.pressureGlobalFactorNonZeros
        <<",global_storage_bytes,"<<globalDiagnostics.pressureGlobalStorageBytes
        <<",global_admission_bytes,"<<globalDiagnostics.pressureGlobalAdmissionBytes
        <<",global_solve_calls,"<<globalDiagnostics.pressureGlobalSolveCalls
        <<",global_iterations,"<<globalDiagnostics.krylovIterations
        <<",same_owner_RHS_reuse,1,invalid_global_config_rejected,1"
        <<",singular_nonfinite_rejected,1,original_physical_gates,no_speed_or_1s_claim\n";
}
}


// F09 additive fragment: insert this namespace before the existing main.
// Add one main selector branch shown after the namespace; default/old tests
// stay untouched. No new production hook, header, timing gate or counter.
namespace SpdOrderedRowsContract {
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE) \
    && !defined(__FAST_MATH__) && (!defined(_M_FP_FAST) || !_M_FP_FAST)
using Floating=SparseValidationPlanContract::Floating;
using RestoreFloating=SparseValidationPlanContract::RestoreFloating;
using Block=HydroSaddleVelocityBlock;
using Tier=HydroVelocitySpdTier;
std::size_t pairs=0,negativePairs=0,physicalSolves=0,coefficientWords=0;
struct OmpRestore {int threads=omp_get_max_threads(),dynamic=omp_get_dynamic();
    ~OmpRestore()noexcept{omp_set_dynamic(dynamic);omp_set_num_threads(threads);}};
void threads(int count){omp_set_dynamic(0);omp_set_num_threads(count);}
void prime(bool extra){check(std::feclearexcept(FE_ALL_EXCEPT)==0,"SPD row test cannot clear flags");
    _mm_setcsr(_mm_getcsr()&~0x3fu);
    check(std::feraiseexcept(extra?FE_INEXACT|FE_UNDERFLOW:FE_INEXACT)==0,"SPD row test cannot prime sticky status");
    _mm_setcsr(_mm_getcsr()|0x20u);errno=extra?EDOM:ERANGE;
    Floating status;check((status.native._Fe_stat&FE_INEXACT)!=0&&(status.csr&0x20u)!=0,"SPD row fixture lacks actual native and raw inexact");}
void sameSpd(const HydroVelocitySpdReport& a,const HydroVelocitySpdReport& b){
    check(a.tier==b.tier&&a.gramAttempted==b.gramAttempted&&a.gramBoundsAvailable==b.gramBoundsAvailable,
        "SPD rows changed tier or Gram availability");
    check(SparseValidationPlanContract::bits(a.gramBoundUpper)==SparseValidationPlanContract::bits(b.gramBoundUpper)
        &&SparseValidationPlanContract::bits(a.scaledMassLower)==SparseValidationPlanContract::bits(b.scaledMassLower),
        "SPD rows changed exact Gram bound or scaled mass bits");}
struct Pair {Block one,eight;};
template<class Factory>Pair factoryPair(const Matrix& source,const Matrix& expected,Factory factory,Tier tier,const char* name){
    RestoreFloating restore;const Matrix saved=source;prime((pairs&1u)!=0);const Floating entry;
    threads(1);entry.restore();const auto one=factory();const Floating afterOne;
    threads(8);entry.restore();const auto eight=factory();const Floating afterEight;
    SparseValidationPlanContract::sameFlags(afterOne,afterEight);
    SparseValidationPlanContract::sameMatrix(one.matrix(),eight.matrix());
    SparseValidationPlanContract::sameMatrix(one.matrix(),expected);SparseValidationPlanContract::sameMatrix(source,saved);
    sameSpd(one.spdReport(),eight.spdReport());
    check(one.spdReport().tier==tier,"SPD row independent fixture obtained wrong certificate tier");
    check(SparseValidationPlanContract::bits(one.symmetryDefect())==SparseValidationPlanContract::bits(eight.symmetryDefect()),
        "SPD rows changed the original symmetry diagnostic bits");
    coefficientWords+=static_cast<std::size_t>(one.matrix().nonZeros());++pairs;
    std::cout<<"SPD_ORDERED_ROWS_PAIR,"<<name<<",rows,"<<source.rows()<<",requested_workers,8,tier,"<<static_cast<int>(tier)
        <<",CSR_values_and_full_FENV_CSR_errno,original\n";
    return {one,eight};}
template<class Factory>void rejectionPair(Factory factory,const char* expected=nullptr){
    RestoreFloating restore;prime(true);const Floating entry;
    const auto invoke=[&](){std::string message;int kind=0;try{factory();}
        catch(const std::invalid_argument& e){kind=1;message=e.what();}
        catch(const std::runtime_error& e){kind=2;message=e.what();}return std::make_pair(kind,message);};
    threads(1);entry.restore();const auto one=invoke();const Floating afterOne;
    threads(8);entry.restore();const auto eight=invoke();const Floating afterEight;
    SparseValidationPlanContract::sameFlags(afterOne,afterEight);
    check(one.first!=0&&one==eight,"SPD row public guard changed rejection type/first message");
    if(expected)check(one.second.find(expected)!=std::string::npos,"SPD energy negative failed another original guard");
    ++negativePairs;}
Matrix tridiagonal(int n,double d){Matrix a=diagonal(std::vector<double>(n,d));
    for(int i=1;i<n;++i){a.coeffRef(i-1,i)=-.125;a.coeffRef(i,i-1)=-.125;}a.makeCompressed();return a;}
// The original RZ reference supplies real face areas/incidence and cell volume.
// The manufactured pressure load is assembled from those independent terms,
// without the production C/D or quotient matrix. Retain the original u/p gates.
void manufactured(const QuadtreeMesh& mesh,const HydroDofLayout& layout,const Reference& ref,const Block& velocity){
    threads(8);const int nu=layout.velocityCount(),np=layout.pressureCount();
    std::vector<double> u(nu),p(np),known(nu),weights(nu,1.);double mean=0,volume=0;
    for(int i=0;i<nu;++i){const auto& f=ref.faces[i];
        const double flow=f.horizontal?polynomialPsi(f.x0,f.y0,mesh)-polynomialPsi(f.x1,f.y0,mesh)
            :polynomialPsi(f.x0,f.y1,mesh)-polynomialPsi(f.x0,f.y0,mesh);
        u[i]=f.area==0?0:2*pi*flow/f.area;if(f.incidence==1||f.area==0)known[i]=u[i];}
    for(int i=0;i<np;++i){const auto c=ref.center[i];p[i]=1+.3*c.x+.7*c.y+.4*c.x*c.y;
        mean+=ref.volume[i]*p[i];volume+=ref.volume[i];}
    check(volume>0,"SPD row MMS lacks positive physical RZ volume");for(double& value:p)value-=mean/volume;
    auto momentum=multiply(velocity.matrix(),u);for(const auto& term:ref.terms)momentum[term.u]-=term.outwardArea*p[term.p];
    const auto load=HydroPressureLoad::build(HydroIntegralDivergence::build(mesh,layout));
    const auto solver=HydroSinglePhaseSaddle::buildClosed(load,velocity);
    const auto result=solver.solve(mesh,HydroSaddleRhs(layout,momentum,known));SaddleScatterPlanContract::originalGates(result);
    check(relativeError(result.values().normalVelocity(),u,weights)<=1e-9
        &&relativeError(result.values().pressure(),p,ref.volume)<=1e-9,"SPD ordered rows failed original manufactured u/p gates");++physicalSolves;}
void run(){
    RestoreFloating floating;OmpRestore omp;check(std::fesetround(FE_TONEAREST)==0,"SPD row fixture requires RNE");
    _mm_setcsr(_mm_getcsr()&~0x8040u);unsigned controls=0;
    check(_controlfp_s(&controls,0,0)==0&&(controls&_MCW_EM)==_MCW_EM&&(_mm_getcsr()&0x1f80u)==0x1f80u,
        "SPD row fixture requires masked native/raw controls");
    threads(8);int formed=0;
#pragma omp parallel num_threads(8)
    {
        if(omp_get_thread_num()==0)formed=omp_get_num_threads();
    }
    check(formed==8,"SPD row positive fixture cannot form requested eight-worker team");
    // Reuse the existing PressureCoarse input/refine/reference geometry scale.
    QuadtreeMesh mesh(input(16,1));refine(mesh,2);FaceSegmentGraph graph;graph.build(mesh);
    const auto layout=HydroDofLayout::build(mesh,NativeLeafDomain::build(mesh,graph,true));
    const auto ref=reference(mesh,graph,layout,true);const int n=layout.velocityCount();
    check(!layout.domain().uniformLeaves()&&layout.pressureCount()>512&&n>=512,"SPD row fixture lacks actual RZ2:1/multiworker rows");
    const auto build=[&](const Matrix& a,Tier tier,const char* name){return factoryPair(a,a,[&]{return Block::build(layout,a);},tier,name);};
    const auto dominant=tridiagonal(n,4.);const auto first=build(dominant,Tier::StrictDiagonalDominance,"normal dominant");
    manufactured(mesh,layout,ref,first.eight);
    // Direct sum of [[1,2],[2,5]] and positive scalar blocks: leading minors
    // 1 and 1 prove SPD; its first row is genuinely not diagonally dominant.
    auto nonDd=diagonal(std::vector<double>(n,4.));nonDd.coeffRef(0,0)=1.;nonDd.coeffRef(1,1)=5.;
    nonDd.coeffRef(0,1)=nonDd.coeffRef(1,0)=2.;nonDd.makeCompressed();build(nonDd,Tier::LLT,"analytic SPD non-DD normal false");
    const double d1=std::nextafter(.25,std::numeric_limits<double>::infinity());
    const double d2=std::nextafter(d1,std::numeric_limits<double>::infinity());
    // For each interior row, literal .125/.125 additions and successors give
    // radius .25+one-ulp. Both Toeplitz matrices are independently SPD since
    // d>.25=2*.125, but only the two-ulp diagonal passes this strict bound.
    build(tridiagonal(n,d1),Tier::LLT,"near-bound one-ulp diagonal");
    build(tridiagonal(n,d2),Tier::StrictDiagonalDominance,"near-bound two-ulp diagonal");
    // At the centre, original two-segment orders are [1,eps/2,eps/2]
    // and [eps/2,eps/2,1], yielding literal radii 1+5eps and 1+2eps.
    // In both cases diag=1+4eps exceeds the exact sum 1+eps: strict SPD.
    const int middle=n/2;const double eps=std::numeric_limits<double>::epsilon();
    auto largeFirst=diagonal(std::vector<double>(n,4.));largeFirst.coeffRef(middle,middle)=1+4*eps;
    largeFirst.coeffRef(middle,0)=largeFirst.coeffRef(0,middle)=1.;
    for(int j:{middle+1,middle+2})largeFirst.coeffRef(j,middle)=largeFirst.coeffRef(middle,j)=.5*eps;
    largeFirst.makeCompressed();build(largeFirst,Tier::LLT,"two-segment large-first exact order");
    auto smallFirst=diagonal(std::vector<double>(n,4.));smallFirst.coeffRef(middle,middle)=1+4*eps;
    for(int j:{0,1})smallFirst.coeffRef(middle,j)=smallFirst.coeffRef(j,middle)=.5*eps;
    smallFirst.coeffRef(middle+1,middle)=smallFirst.coeffRef(middle,middle+1)=1.;
    smallFirst.makeCompressed();build(smallFirst,Tier::StrictDiagonalDominance,"two-segment small-first exact order");
    auto zeros=dominant;zeros.coeffRef(0,2)=0.;zeros.coeffRef(2,0)=-0.;zeros.makeCompressed();
    check(std::signbit(zeros.coeff(2,0)),"SPD row fixture lost stored negative zero");build(zeros,Tier::StrictDiagonalDominance,"explicit signed zero");
    auto normalTiny=dominant;for(Eigen::Index i=0;i<normalTiny.nonZeros();++i)normalTiny.valuePtr()[i]=std::ldexp(normalTiny.valuePtr()[i],-480);
    build(normalTiny,Tier::StrictDiagonalDominance,"normal tiny no-underflow certificate");
    const double q=std::numeric_limits<double>::denorm_min();check(q>0,"SPD row fallback fixture lacks real subnormal");
    auto tiny=dominant;tiny.coeffRef(0,2)=tiny.coeffRef(2,0)=q;tiny.makeCompressed();
    build(tiny,Tier::LLT,"subnormal raw-entry literal fallback");
    // Upper is not used by Lower LLT, but the original certificate validates
    // every raw entry. Its q-sized skew is beneath the original norm gate.
    auto upperTiny=dominant;upperTiny.coeffRef(0,2)=q;upperTiny.coeffRef(2,0)=0.;upperTiny.makeCompressed();
    build(upperTiny,Tier::LLT,"upper-only subnormal still forces original guard");
    for(double bad:{std::numeric_limits<double>::infinity(),std::numeric_limits<double>::quiet_NaN()}){
        auto invalid=dominant;invalid.coeffRef(0,0)=bad;rejectionPair([&]{return Block::build(layout,invalid);});}
    // Same actual geometry/mass/raw-volume strain construction as the existing
    // sparse-validation and Gram fixtures. No large dense eigenproblem.
    const auto metrics=HydroPhysicalDivergence::build(mesh,HydroIntegralDivergence::build(mesh,layout));
    const auto owner=HydroMomentumGeometry::build(mesh,metrics);
    const LAMP::TimeLayer clock{0,0,0,LAMP::TimeLayer::Position::Initial};
    const auto pub=LAMP::Interface::GeometricVofPublication::prepareLocalPlic(mesh,layout.domain(),clock,std::vector<double>(mesh.cellIndexNum(),1.));
    const auto mass=HydroMomentumMass::buildCellAverageProjection(mesh,owner,HydroPhaseGeometry::build(mesh,metrics,pub),1.,1.);
    const auto strain=HydroStrainOperator::build(mesh,metrics,HydroStrainOperator::volumeQuadrature(mesh,metrics));
    std::vector<double> mu(strain.points().size());for(std::size_t i=0;i<mu.size();++i)mu[i]=i%7==0?0.:.25+.125*(i%4);
    const auto energy=HydroEnergyViscosity::build(strain,mu);constexpr double dt=1.;Matrix a=.5*energy.matrix();
    for(int i=0;i<n;++i){check(mass.masses()[i]>0,"actual Gram mass is not strictly positive");a.coeffRef(i,i)+=mass.masses()[i]/dt;}a.makeCompressed();
    const Matrix transpose(a.transpose()),s=.5*(a+transpose),w=.5*(a-transpose),combined=s+w;
    // A wide-margin independent row-sum witness ensures this
    // fixture reaches beyond the DD sufficient gate before the Gram path.
    std::vector<long double> radii(n,0.);for(Eigen::Index i=0;i<s.outerSize();++i)for(Matrix::InnerIterator e(s,i);e;++e)
        if(e.row()>e.col()){const long double magnitude=std::abs(static_cast<long double>(e.value()));radii[e.row()]+=magnitude;radii[e.col()]+=magnitude;}
    bool nonDominant=false;for(int i=0;i<n;++i)nonDominant=nonDominant||radii[i]>1.001L*static_cast<long double>(s.coeff(i,i));
    check(nonDominant,"actual mass+Gram fixture did not exercise non-DD");
    const auto gram=factoryPair(a,combined,[&]{return Block::buildWithEnergyTransportSplit(mesh,layout,a,mass,energy,dt,owner);},Tier::EnergyGram,"actual mass and energy Gram maxima");
    const auto& gr=gram.eight.spdReport();check(gr.gramAttempted&&gr.gramBoundsAvailable&&std::isfinite(gr.gramBoundUpper)
        &&gr.gramBoundUpper>=0&&gr.scaledMassLower>gr.gramBoundUpper&&gr.lltSeconds==0,"actual eight-thread request lacks retained strict Gram certificate");
    // Strictly positive extra diagonal is a valid coercive contribution, not
    // an error requiring an absolute-value penalty. The SAME stored Lower S
    // has an independent original LLT reference, without an energy hint.
    Matrix dissipative=a;
    for(int i=0;i<n;++i)dissipative.coeffRef(i,i)+=4.*(mass.masses()[i]/dt);
    dissipative.makeCompressed();
    const Matrix dTranspose(dissipative.transpose()),dSym=.5*(dissipative+dTranspose),
        dSkew=.5*(dissipative-dTranspose),dCombined=dSym+dSkew;
    bool extraNonDominant=false;for(int i=0;i<n;++i)
        extraNonDominant=extraNonDominant||radii[i]>1.001L*static_cast<long double>(dSym.coeff(i,i));
    check(extraNonDominant,"positive diagonal energy witness became DD and missed the required path");
    const auto independent=factoryPair(dSym,dSym,[&]{return Block::build(layout,dSym);},Tier::LLT,
        "actual 4mass positive diagonal independent original LLT");
    (void)independent;
    // One independently enclosed exact stored diagonal Rayleigh quotient
    // excludes BOTH old coarse |R| and old tight |R| certificates. It does not
    // reproduce their support scan/roundoff algorithm or inspect a private API.
    // For row zero: D*R*D > 2.5 while D*mass*D < 1.01; the minimum scaled mass
    // can only be smaller. Both old error norms include this positive entry.
    const double infinity=std::numeric_limits<double>::infinity(),d0=mass.masses()[0]/dt;
    volatile double halfK0=.5*energy.matrix().coeff(0,0);
    volatile double firstDifference=dSym.coeff(0,0)-halfK0;
    const double firstLower=std::nextafter(firstDifference,-infinity);
    volatile double secondDifference=firstLower-d0;
    const double residualLower=std::nextafter(secondDifference,-infinity);
    volatile double root0=std::sqrt(d0),scale0=1./root0;
    volatile double residualFirst=scale0*residualLower;
    const double residualFirstLower=std::nextafter(residualFirst,0.);
    volatile double residualSecond=residualFirstLower*scale0;
    const double scaledResidualLower=std::nextafter(residualSecond,0.);
    volatile double massFirst=scale0*d0;
    const double massFirstUpper=std::nextafter(massFirst,infinity);
    volatile double massSecond=massFirstUpper*scale0;
    const double scaledMassUpper=std::nextafter(massSecond,infinity);
    check(std::isfinite(scaledResidualLower)&&std::isfinite(scaledMassUpper)&&residualLower>0&&
        scaledResidualLower>2.5&&scaledMassUpper<1.01,"4mass stored Rayleigh witness did not exclude old absolute-residual paths");
    const auto positive=factoryPair(dissipative,dCombined,
        [&]{return Block::buildWithEnergyTransportSplit(mesh,layout,dissipative,mass,energy,dt,owner);},
        Tier::EnergyGram,"actual positive 4mass new coarse spectral certificate");
    const auto& positiveReport=positive.eight.spdReport();
    check(positiveReport.gramAttempted&&positiveReport.gramBoundsAvailable&&positiveReport.lltSeconds==0&&
        positiveReport.gramBoundUpper<positiveReport.scaledMassLower&&positiveReport.scaledMassLower<1.01,
        "positive diagonal did not consume the new sufficient coarse certificate");
    // Genuine indefinite inputs, not unsupported FP or an unrelated guard.
    // A negative diagonal gives a negative coordinate Rayleigh quotient.
    auto negativeDiagonal=dissipative;negativeDiagonal.coeffRef(0,0)=-std::abs(dSym.coeff(0,0));negativeDiagonal.makeCompressed();
    check(negativeDiagonal.coeff(0,0)<0,"negative diagonal Rayleigh witness");
    rejectionPair([&]{return Block::buildWithEnergyTransportSplit(mesh,layout,negativeDiagonal,mass,energy,dt,owner);},
        "velocity block is not positive definite");
    // Positive diagonals but a negative 2x2 principal minor: a>|d0|,|d1|
    // implies a*a>d0*d1 without any potentially overflowing product oracle.
    auto indefinitePair=dissipative;const double largeOff=2.*std::max(dSym.coeff(0,0),dSym.coeff(1,1));
    check(std::isfinite(largeOff)&&dSym.coeff(0,0)>0&&dSym.coeff(1,1)>0&&
        largeOff>dSym.coeff(0,0)&&largeOff>dSym.coeff(1,1),"offdiagonal principal-minor witness");
    indefinitePair.coeffRef(0,1)=indefinitePair.coeffRef(1,0)=largeOff;indefinitePair.makeCompressed();
    rejectionPair([&]{return Block::buildWithEnergyTransportSplit(mesh,layout,indefinitePair,mass,energy,dt,owner);},
        "velocity block is not positive definite");
    check(pairs==13&&negativePairs==4&&physicalSolves==1&&coefficientWords>0,"SPD row fixture coverage incomplete");
    std::cout<<"SPD_ORDERED_ROWS_PASS,RZ_real2to1,1,pressure_rows,"<<layout.pressureCount()<<",velocity_rows,"<<n
        <<",normal_team_formation,"<<formed<<",factory_1_vs_8_pairs,"<<pairs<<",original_public_rejection_pairs,"<<negativePairs
        <<",nonzero_original_u_p_and_physical_gate_solves,"<<physicalSolves<<",matrix_words,"<<coefficientWords
        <<",actual_mass_energy_Gram_pair,1,positive_4mass_new_coarse_pair,1,independent_LLT_pair,1,indefinite_energy_guards,2"
        <<",no_public_used_counter_or_fault_injection,no_speed_claim\n";
}
#else
void run(){std::cout<<"SPD_ORDERED_ROWS_UNAVAILABLE,requires_original_strict_MSVC_x64_OpenMP; old selectors unchanged\n";}
#endif
}

namespace NativeInitialGuessContract {
void runCase(bool anchored) {
    QuadtreeMesh mesh(input(8,0));const auto graph=selected(mesh,anchored?0:1);
    const auto layout=HydroDofLayout::build(mesh,NativeLeafDomain::build(mesh,graph,false));
    const auto load=HydroPressureLoad::build(HydroIntegralDivergence::build(mesh,layout));
    const auto ref=reference(mesh,graph,layout,false);const int nu=layout.velocityCount(),np=layout.pressureCount();
    Matrix a=diagonal(std::vector<double>(nu,3.));
    for(int i=1;i<nu;++i){a.coeffRef(i-1,i)=-.25;a.coeffRef(i,i-1)=-.25;}a.makeCompressed();
    const auto velocity=HydroSaddleVelocityBlock::build(layout,a);
    HydroSaddleConfig config;config.gmresOnlyParallel=true;config.gmresGlobalPressure=true;
    std::vector<HydroSaddlePressureOutlet> outlets;
    if(anchored)for(int u=0;u<nu;++u){const auto& f=ref.faces[u];
        if(f.incidence==1&&f.horizontal&&f.y0==mesh.maxCoord.y)outlets.push_back({u,2.5});}
    const auto solver=anchored?HydroSinglePhaseSaddle::buildWithPressureOutlets(load,velocity,outlets,config):
        HydroSinglePhaseSaddle::buildClosed(load,velocity,config);
    check(solver.pressureGaugeCount()==(anchored?0:2),"initial-guess fixture lost anchored or disconnected closed components");
    std::vector<bool> free(nu);for(int u:solver.freeVelocityDofs())free[u]=true;
    std::vector<double> u(nu),p(np),known(nu),weights(nu,1.);
    for(int i=0;i<nu;++i){u[i]=ref.faces[i].horizontal?1.:0.;if(!free[i])known[i]=u[i];}
    for(int i=0;i<np;++i)p[i]=3.+.25*ref.center[i].x+.5*ref.center[i].y;
    if(!anchored){std::vector<double> mean(2),volume(2);
        for(int i=0;i<np;++i){const int group=solver.pressureComponents()[i];mean[group]+=ref.volume[i]*p[i];volume[group]+=ref.volume[i];}
        for(int i=0;i<np;++i){const int group=solver.pressureComponents()[i];p[i]-=mean[group]/volume[group];}}
    const auto momentumFor=[&](const std::vector<double>& actualU,const std::vector<double>& actualP){
        auto value=multiply(a,actualU);for(const auto& term:ref.terms)value[term.u]-=term.outwardArea*actualP[term.p];
        for(int i=0;i<nu;++i)value[i]-=solver.pressureOutletLoads()[i];return value;};
    const auto verify=[&](const HydroSaddleResult& result,const std::vector<double>& expectedU,
        const std::vector<double>& expectedP,const std::vector<double>& momentum){
        report(result);check(result.report().usedBlockGmres&&!result.report().usedSuperLU&&!result.report().usedDirectFallback,
            "native initial iterate changed the explicit original backend");
        check(relativeError(result.values().normalVelocity(),expectedU,weights)<=1e-9
            &&relativeError(result.values().pressure(),expectedP,ref.volume)<=1e-9,"native initial iterate lost independent U/p");
        auto physical=momentum;for(int i=0;i<nu;++i)physical[i]+=solver.pressureOutletLoads()[i];
        const auto& actual=result.values().normalVelocity();const auto energy=dot(actual,multiply(a,actual)),work=dot(actual,physical);
        check(std::abs(energy-work)/(std::abs(energy)+std::abs(work))<=1e-10,"native initial iterate lost original global power identity");};
    const auto momentum=momentumFor(u,p);const HydroSaddleRhs rhs(layout,momentum,known);
    const auto cold=solver.solve(mesh,rhs);verify(cold,u,p,momentum);
    auto guessedU=u,guessedP=p;
    for(int i=0;i<nu;++i)if(!free[i])guessedU[i]=13.; // Current fixed rows come only from RHS.
    if(!anchored)for(int i=0;i<np;++i)guessedP[i]+=solver.pressureComponents()[i]==0?3.25:-7.5;
    const HydroNativeValues exactGuess(layout,guessedU,guessedP);
    const auto exact=solver.solve(mesh,rhs,exactGuess);verify(exact,u,p,momentum);
    check(exact.report().krylovIterations==0&&exact.report().linearRefinements==0,
        "exact native free/gauge initial iterate was discarded instead of executing zero Arnoldi iterations");
    auto nextU=u,nextP=p,nextKnown=known;const double scale=1.+std::ldexp(1.,-20);
    for(auto* values:{&nextU,&nextP,&nextKnown})for(double& value:*values)value*=scale;
    const auto nextMomentum=momentumFor(nextU,nextP);const HydroSaddleRhs nextRhs(layout,nextMomentum,nextKnown);
    const auto nextCold=solver.solve(mesh,nextRhs),nextWarm=solver.solve(mesh,nextRhs,HydroNativeValues(layout,u,p));
    verify(nextCold,nextU,nextP,nextMomentum);verify(nextWarm,nextU,nextP,nextMomentum);
    check(nextCold.report().krylovIterations>0&&nextWarm.report().krylovIterations<nextCold.report().krylovIterations,
        "nearby original-RHS initial iterate did not reduce actual Arnoldi work");
    std::vector<double> zeroU(nu),zeroP(np),zeroMomentum(nu);
    for(int i=0;i<nu;++i)zeroMomentum[i]=-solver.pressureOutletLoads()[i];
    const auto zero=solver.solve(mesh,HydroSaddleRhs(layout,zeroMomentum,zeroU),exactGuess);
    check(zero.values().normalVelocity()==zeroU&&zero.values().pressure()==zeroP&&zero.report().krylovIterations==0,
        "zero original RHS retained a nonzero native initial iterate");
    const auto fresh=HydroDofLayout::build(mesh,NativeLeafDomain::build(mesh,graph,false));
    rejects([&]{solver.solve(mesh,rhs,HydroNativeValues(fresh,u,p));},"different layout identity");
    auto poisoned=exactGuess;const_cast<std::vector<double>&>(poisoned.normalVelocity())[0]=std::numeric_limits<double>::quiet_NaN();
    rejects([&]{solver.solve(mesh,rhs,poisoned);},"nonfinite native initial guess");
    auto poisonedP=exactGuess;const_cast<std::vector<double>&>(poisonedP.pressure())[0]=std::numeric_limits<double>::infinity();
    rejects([&]{solver.solve(mesh,rhs,poisonedP);},"nonfinite native initial guess");
    const auto ordinary=anchored?HydroSinglePhaseSaddle::buildWithPressureOutlets(load,velocity,outlets):
        HydroSinglePhaseSaddle::buildClosed(load,velocity);
    rejects([&]{ordinary.solve(mesh,rhs,exactGuess);},"explicit GMRES-only mode");
    const auto mixedVelocity=HydroSaddleVelocityBlock::build(layout,diagonal(std::vector<double>(nu+1,3.)),1);
    const auto mixed=HydroSinglePhaseSaddle::buildClosed(load,mixedVelocity,config);
    rejects([&]{mixed.solve(mesh,HydroSaddleRhs(layout,momentum,known,{0.}),exactGuess);},"tangential unknowns");
    auto& vertex=mesh.verts(mesh.vertIndexToLabel[0]);const auto saved=vertex.coord();vertex.setCoord(saved+Vector2D{.01,0});
    rejects([&]{solver.solve(mesh,rhs,exactGuess);});vertex.setCoord(saved);
    std::cout<<"NATIVE_INITIAL_GUESS_CASE,"<<(anchored?"anchored":"two_closed_gauges")
        <<",exact_iterations,"<<exact.report().krylovIterations<<",cold_iterations,"<<nextCold.report().krylovIterations
        <<",nearby_iterations,"<<nextWarm.report().krylovIterations<<",original_matrix_rhs_and_power_gates,1\n";
}
void run(){runCase(false);runCase(true);std::cout<<"NATIVE_INITIAL_GUESS_PASS\n";}
}

int main(int argc,char** argv) try {
    if(argc==2&&std::string(argv[1])=="--native-initial-guess-only"){std::cout<<std::setprecision(17);NativeInitialGuessContract::run();return 0;}
    if(argc==2&&std::string(argv[1])=="--spd-ordered-rows-only"){std::cout<<std::setprecision(17);SpdOrderedRowsContract::run();return 0;}
    if(argc==2&&std::string(argv[1])=="--pressure-coarse-only"){std::cout<<std::setprecision(17);PressureCoarseContract::run();return 0;}
    if(argc==2&&std::string(argv[1])=="--parallel-block-gmres-only"){std::cout<<std::setprecision(17);ParallelBlockGmresContract::run();return 0;}
    if(argc==2&&std::string(argv[1])=="--lazy-direct-input-only"){std::cout<<std::setprecision(17);LazyDirectInputContract::run();return 0;}
    if(argc==2&&std::string(argv[1])=="--ordered-row-spmv-only"){std::cout<<std::setprecision(17);OrderedRowSpmvContract::run();return 0;}
    if(argc==2&&std::string(argv[1])=="--velocity-split-plan-only"){std::cout<<std::setprecision(17);VelocitySplitPlanContract::run();return 0;}
    if(argc==2&&std::string(argv[1])=="--saddle-scatter-plan-only"){std::cout<<std::setprecision(17);SaddleScatterPlanContract::run();return 0;}
    if(argc==2&&std::string(argv[1])=="--sparse-validation-plan-only"){std::cout<<std::setprecision(17);SparseValidationPlanContract::run();return 0;}
    const bool quotientOnly=argc==2&&std::string(argv[1])=="--quotient-only";
    const bool outletOnly=argc==2&&std::string(argv[1])=="--pressure-outlet-only";
    check(argc==1||quotientOnly||outletOnly,"unknown saddle contract argument");
    if(quotientOnly){std::cout<<std::setprecision(17);quotientReferenceWitness();return 0;}
    if(outletOnly){std::cout<<std::setprecision(17);pressureOutletWitness();return 0;}
    for (int metric : {0,1,2}) for (int n : {4,8}) for (int depth : {0,2}) for (int selection : {0,1,2}) algebra(n,metric,depth,selection);
    for (int metric : {0,1,2}) algebra(4,metric,0,0,true,true);
    const auto coarse=mms(4),middle=mms(8),fine=mms(16);
    check(coarse.first>middle.first && middle.first>fine.first && coarse.second>middle.second && middle.second>fine.second,"MMS errors do not decrease");
    velocityOrder=std::log(middle.first/fine.first)/std::log(2.0);pressureOrder=std::log(middle.second/fine.second)/std::log(2.0);
    check(velocityOrder>=1.8 && pressureOrder>=1.8,"pressure/velocity MMS misses frozen order gate");
    guardsAndLifetime();check(algebraCases==39 && mmsCases==3 && rejected==22,"incomplete saddle coverage");
    quotientReferenceWitness();pressureOutletWitness();
    std::cout << std::scientific << std::setprecision(16) << "H04 HYDRO saddle PASS algebra=" << algebraCases << " mms=" << mmsCases << " rejected=" << rejected
              << " velocity=" << peakVelocity << " pressure=" << peakPressure << " energy=" << peakEnergy << " continuity=" << peakContinuity
              << " residual=" << peakAlgebraic << " mean=" << peakMean << " velocity_order=" << velocityOrder << " pressure_order=" << pressureOrder << '\n';return 0;
}
catch (const std::exception& error) {std::cerr << "H04 FAIL: " << error.what() << '\n';return 1;}
