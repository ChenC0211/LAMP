#include "Modules/Hydrodynamics/CchSinglePhaseSaddle.h"
#include "Core/IO/DictionaryMesh.h"
#include "Core/Mesh/FaceSegmentGraph.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace
{
using namespace Quantum::Hydrodynamics;
using Matrix=CchSaddleVelocityBlock::Matrix;
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
Reference reference(const QuadtreeMesh& mesh,const FaceSegmentGraph& graph,const CchDofLayout& layout,bool rz) {
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
void report(const CchSaddleResult& result) {
    const auto& r=result.report();peakContinuity=std::max(peakContinuity,r.maxCellContinuityDefect);
    peakAlgebraic=std::max(peakAlgebraic,r.algebraicRelativeResidual);peakMean=std::max(peakMean,r.maxPressureMeanDefect);
    check(r.algebraicRelativeResidual<=1e-11 && r.freeMomentumDefect<=1e-11 && r.maxCellContinuityDefect<=1e-11 && r.maxPressureMeanDefect<=1e-11,
        "unaccepted physical/algebraic result returned");
}
double polynomialPsi(double x,double y,const QuadtreeMesh& mesh) {
    const double a=(x-mesh.minCoord.x)*(mesh.maxCoord.x-x),b=y*(1-y);return a*a*b*b;
}
void algebra(int n,int metric,int depth,int selection,bool nonDiagonal=false,bool axial=false) {
    QuadtreeMesh mesh(input(n,metric));refine(mesh,depth);const auto graph=selected(mesh,selection);
    const auto layout=CchDofLayout::build(mesh,NativeLeafDomain::build(mesh,graph,metric!=0));
    const auto load=CchPressureLoad::build(CchIntegralDivergence::build(mesh,layout));const auto ref=reference(mesh,graph,layout,metric!=0);
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
    const auto velocity=CchSaddleVelocityBlock::build(layout,a);const auto solver=CchSinglePhaseSaddle::buildClosed(load,velocity);
    check(solver.pressureGaugeCount()==static_cast<int>(mean.size()),"wrong component pressure gauge count");
    check(&solver.pressureLoad().divergence().matrix()==&load.divergence().matrix(),"solver lost unique divergence");
    const CchSinglePhaseSaddle::Matrix asym=solver.matrix()-CchSinglePhaseSaddle::Matrix(solver.matrix().transpose());
    check(asym.norm()==0,"closed augmented saddle is not exactly symmetric");
    auto momentum=multiply(a,u);for (const auto& term : ref.terms) momentum[term.u]-=term.outwardArea*p[term.p];
    const CchSaddleRhs rhs(layout,momentum,known);const auto result=solver.solve(mesh,rhs);report(result);
    const double eu=relativeError(result.values().normalVelocity(),u,weights),ep=relativeError(result.values().pressure(),p,ref.volume);
    peakVelocity=std::max(peakVelocity,eu);peakPressure=std::max(peakPressure,ep);check(eu<=1e-9 && ep<=1e-9,"independent algebra recovery failed");
    const auto repeated=solver.solve(mesh,rhs);
    check(repeated.values().normalVelocity()==result.values().normalVelocity() && repeated.values().pressure()==result.values().pressure(),"reused factor changed deterministic RHS");
    if (selection==0 && !axial) {
        const auto actual=result.values().normalVelocity();const double energy=dot(actual,multiply(a,actual)),work=dot(actual,momentum);
        const double defect=std::abs(energy-work)/(std::abs(energy)+std::abs(work));peakEnergy=std::max(peakEnergy,defect);check(defect<=1e-10,"zero-boundary saddle energy relation failed");
    }
    ++algebraCases;
    std::cout << "CCH_SADDLE_ALG," << metric << ',' << n << ',' << depth << ',' << selection << ',' << eu << ',' << ep << '\n';
}
double sinc(double x) {return x==0 ? 1 : std::sin(x)/x;}
double analyticPsi(double x,double y) {return x==0 || x==1 || y==0 || y==1 ? 0 : std::sin(pi*x)*std::sin(pi*y);}
std::pair<double,double> mms(int n) {
    QuadtreeMesh mesh(input(n,0));const auto graph=selected(mesh,0);
    const auto layout=CchDofLayout::build(mesh,NativeLeafDomain::build(mesh,graph,false));const auto ref=reference(mesh,graph,layout,false);
    const auto load=CchPressureLoad::build(CchIntegralDivergence::build(mesh,layout));
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
    const auto solver=CchSinglePhaseSaddle::buildClosed(load,CchSaddleVelocityBlock::build(layout,diagonal(mass)));
    const auto result=solver.solve(mesh,CchSaddleRhs(layout,momentum,known));report(result);
    const double eu=relativeError(result.values().normalVelocity(),u,dual),ep=relativeError(result.values().pressure(),p,ref.volume);
    ++mmsCases;std::cout << "CCH_SADDLE_MMS," << n << ',' << std::setprecision(16) << eu << ',' << ep << '\n';return {eu,ep};
}
void guardsAndLifetime() {
    static_assert(!std::is_default_constructible_v<CchSinglePhaseSaddle>);
    static_assert(std::is_same_v<decltype(std::declval<const CchSaddleResult&>().values()),const CchNativeValues&>);
    QuadtreeMesh mesh(input(4,0));const auto graph=selected(mesh,0);const auto layout=CchDofLayout::build(mesh,NativeLeafDomain::build(mesh,graph,false));
    Matrix a=diagonal(std::vector<double>(layout.velocityCount(),1));const auto load=CchPressureLoad::build(CchIntegralDivergence::build(mesh,layout));
    const auto velocity=CchSaddleVelocityBlock::build(layout,a);const auto solver=CchSinglePhaseSaddle::buildClosed(load,velocity);const auto copy=solver;
    check(&copy.matrix()==&solver.matrix(),"saddle copy is not shared read-only");a.coeffRef(0,0)=7;check(velocity.matrix().coeff(0,0)==1,"velocity block aliases caller matrix");
    std::vector<double> zero(layout.velocityCount(),0);const CchSaddleRhs rhs(layout,zero,zero);
    const auto accepted=solver.solve(mesh,rhs);check(accepted.values().normalVelocity()==zero && accepted.values().pressure()==std::vector<double>(layout.pressureCount(),0),"zero RHS is not exact zero");
    auto bad=a;bad.coeffRef(0,0)=-1;rejects([&] {CchSaddleVelocityBlock::build(layout,bad);});
    bad=a;bad.coeffRef(0,0)=0;rejects([&] {CchSaddleVelocityBlock::build(layout,bad);});
    bad=a;bad.coeffRef(0,1)=.1;rejects([&] {CchSaddleVelocityBlock::build(layout,bad);});
    bad=a;bad.coeffRef(0,0)=std::numeric_limits<double>::infinity();rejects([&] {CchSaddleVelocityBlock::build(layout,bad);});
    rejects([&] {CchSaddleVelocityBlock::build(layout,Matrix(2,2));});
    auto config=CchSaddleConfig{};config.algebraicRelativeTolerance=0;rejects([&] {CchSinglePhaseSaddle::buildClosed(load,velocity,config);});
    config=CchSaddleConfig{};config.continuityRelativeTolerance=std::numeric_limits<double>::quiet_NaN();rejects([&] {CchSinglePhaseSaddle::buildClosed(load,velocity,config);});
    const auto fresh=CchDofLayout::build(mesh,NativeLeafDomain::build(mesh,graph,false));
    rejects([&] {CchSinglePhaseSaddle::buildClosed(load,CchSaddleVelocityBlock::build(fresh,velocity.matrix()));});
    rejects([&] {solver.solve(mesh,CchSaddleRhs(fresh,zero,zero));});
    rejects([&] {CchSaddleRhs(layout,std::vector<double>(2),zero);});
    auto payload=zero;payload[0]=std::numeric_limits<double>::quiet_NaN();rejects([&] {CchSaddleRhs(layout,payload,zero);});
    payload=zero;payload[solver.freeVelocityDofs().front()]=1;rejects([&] {solver.solve(mesh,CchSaddleRhs(layout,zero,payload));});
    const auto ref=reference(mesh,graph,layout,false);payload=zero;
    for (std::size_t i=0;i<payload.size();++i) if (ref.faces[i].incidence==1) {payload[i]=1;break;}
    rejects([&] {solver.solve(mesh,CchSaddleRhs(layout,zero,payload));});
    auto stretchedInput=input(4,0);stretchedInput.maxCoord.x=1.25;QuadtreeMesh stretched(stretchedInput);rejects([&] {solver.solve(stretched,rhs);});
    auto& vertex=mesh.verts(mesh.vertIndexToLabel[0]);const auto saved=vertex.coord();vertex.setCoord(saved+Vector2D{.01,0});rejects([&] {solver.solve(mesh,rhs);});vertex.setCoord(saved);
    payload=zero;payload[solver.freeVelocityDofs().front()]=std::numeric_limits<double>::max();rejects([&] {solver.solve(mesh,CchSaddleRhs(layout,payload,zero));});
    payload=zero;payload[solver.freeVelocityDofs().front()]=1;
    config=CchSaddleConfig{};config.algebraicRelativeTolerance=1e-30;
    const auto strictAlgebra=CchSinglePhaseSaddle::buildClosed(load,velocity,config);
    rejects([&] {strictAlgebra.solve(mesh,CchSaddleRhs(layout,payload,zero));},"complete algebraic");
    config=CchSaddleConfig{};config.continuityRelativeTolerance=1e-30;
    const auto strictContinuity=CchSinglePhaseSaddle::buildClosed(load,velocity,config);
    rejects([&] {strictContinuity.solve(mesh,CchSaddleRhs(layout,payload,zero));},"native cell continuity");
    check(solver.solve(mesh,rhs).values().normalVelocity()==accepted.values().normalVelocity(),"failed solve poisoned retained system/result");
    refine(mesh,1);rejects([&] {solver.solve(mesh,rhs);});
    QuadtreeMesh axis(input(4,1));const auto axisLayout=CchDofLayout::build(axis,NativeLeafDomain::build(axis,selected(axis,0),true));
    const auto axisLoad=CchPressureLoad::build(CchIntegralDivergence::build(axis,axisLayout));
    const auto axisSolver=CchSinglePhaseSaddle::buildClosed(axisLoad,CchSaddleVelocityBlock::build(axisLayout,diagonal(std::vector<double>(axisLayout.velocityCount(),1))));
    std::vector<double> axisKnown(axisLayout.velocityCount(),0);
    for (Eigen::Index u=0;u<axisLoad.matrix().outerSize();++u) if (!CchPressureLoad::Matrix::InnerIterator(axisLoad.matrix(),u)) {axisKnown[u]=1;break;}
    rejects([&] {axisSolver.solve(axis,CchSaddleRhs(axisLayout,std::vector<double>(axisKnown.size()),axisKnown));});
    QuadtreeMesh disconnected(input(4,0));const auto parts=selected(disconnected,1);
    const auto partLayout=CchDofLayout::build(disconnected,NativeLeafDomain::build(disconnected,parts,false));
    const auto partLoad=CchPressureLoad::build(CchIntegralDivergence::build(disconnected,partLayout));
    const auto partSolver=CchSinglePhaseSaddle::buildClosed(partLoad,CchSaddleVelocityBlock::build(partLayout,diagonal(std::vector<double>(partLayout.velocityCount(),1))));
    std::vector<double> incompatible(partLayout.velocityCount(),0);std::vector<bool> marked(2,false);
    for (int i=0;i<partLayout.velocityCount();++i) {
        CchPressureLoad::Matrix::InnerIterator it(partLoad.matrix(),i);if (!it) continue;
        const int p=static_cast<int>(it.col());const double coefficient=it.value();++it;
        if (!it && !marked[partSolver.pressureComponents()[p]]) {
            const int component=partSolver.pressureComponents()[p];incompatible[i]=(component==0 ? 1 : -1)/(-coefficient);marked[component]=true;
        }
    }
    check(marked[0] && marked[1],"fixture missed both component boundaries");
    const auto globalFlow=partLoad.divergence().apply(disconnected,CchNativeValues(partLayout,incompatible,std::vector<double>(partLayout.pressureCount())));
    double total=0;for (double flow : globalFlow) total+=flow;check(total==0,"component incompatibility fixture is not globally balanced");
    rejects([&] {partSolver.solve(disconnected,CchSaddleRhs(partLayout,std::vector<double>(incompatible.size()),incompatible));},"component boundary flux compatibility");
    std::optional<CchSinglePhaseSaddle> retained;std::optional<CchSaddleRhs> retainedRhs;
    {QuadtreeMesh temporary(input(4,0));const auto l=CchDofLayout::build(temporary,NativeLeafDomain::build(temporary,selected(temporary,0),false));
     retained=CchSinglePhaseSaddle::buildClosed(CchPressureLoad::build(CchIntegralDivergence::build(temporary,l)),CchSaddleVelocityBlock::build(l,diagonal(std::vector<double>(l.velocityCount(),1))));
     retainedRhs.emplace(l,std::vector<double>(l.velocityCount()),std::vector<double>(l.velocityCount()));}
    QuadtreeMesh equivalent(input(4,0));check(retained->solve(equivalent,*retainedRhs).values().normalVelocity()==zero,"saddle borrowed source lifetime");
}
}
int main() try {
    for (int metric : {0,1,2}) for (int n : {4,8}) for (int depth : {0,2}) for (int selection : {0,1,2}) algebra(n,metric,depth,selection);
    for (int metric : {0,1,2}) algebra(4,metric,0,0,true,true);
    const auto coarse=mms(4),middle=mms(8),fine=mms(16);
    check(coarse.first>middle.first && middle.first>fine.first && coarse.second>middle.second && middle.second>fine.second,"MMS errors do not decrease");
    velocityOrder=std::log(middle.first/fine.first)/std::log(2.0);pressureOrder=std::log(middle.second/fine.second)/std::log(2.0);
    check(velocityOrder>=1.8 && pressureOrder>=1.8,"pressure/velocity MMS misses frozen order gate");
    guardsAndLifetime();check(algebraCases==39 && mmsCases==3 && rejected==21,"incomplete saddle coverage");
    std::cout << std::scientific << std::setprecision(16) << "H04 CCH saddle PASS algebra=" << algebraCases << " mms=" << mmsCases << " rejected=" << rejected
              << " velocity=" << peakVelocity << " pressure=" << peakPressure << " energy=" << peakEnergy << " continuity=" << peakContinuity
              << " residual=" << peakAlgebraic << " mean=" << peakMean << " velocity_order=" << velocityOrder << " pressure_order=" << pressureOrder << '\n';return 0;
}
catch (const std::exception& error) {std::cerr << "H04 FAIL: " << error.what() << '\n';return 1;}
