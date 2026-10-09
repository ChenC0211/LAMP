#include "Modules/Hydrodynamics/CchPressureLoad.h" // Public header stands alone.
#include "Core/IO/DictionaryMesh.h"
#include "Core/Mesh/FaceSegmentGraph.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <random>
#include <sstream>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace
{
using Layout = Quantum::Hydrodynamics::CchDofLayout;
using Values = Quantum::Hydrodynamics::CchNativeValues;
using Divergence = Quantum::Hydrodynamics::CchIntegralDivergence;
using Load = Quantum::Hydrodynamics::CchPressureLoad;
constexpr long double verificationGate = 1e-13L, hardGate = 1e-12L;
const long double pi = std::acos(-1.0L);
int scenarios=0, pairs=0, rejected=0, coarseSegments=0, axisModes=0;
long double adjointDefect=0, loadDefect=0;
void check(bool ok, const char* reason) { if (!ok) throw std::runtime_error(reason); }
template<class F> void rejects(F action)
{
    try { action(); } catch (const std::invalid_argument&) { ++rejected; return; }
    throw std::runtime_error("invalid pressure load source/output accepted");
}
long double norm(const std::vector<double>& values)
{
    long double result=0;
    for (double value : values) result=std::hypot(result,static_cast<long double>(value));
    return result;
}
long double dot(const std::vector<double>& a, const std::vector<double>& b)
{
    check(a.size()==b.size(), "invalid pressure work dimensions");
    long double result=0;
    for (std::size_t i=0; i<a.size(); ++i) result+=static_cast<long double>(a[i])*b[i];
    return result;
}
void observe(long double error, long double scale, long double& peak, bool adjoint=false)
{
    check(std::isfinite(error) && std::isfinite(scale) && scale>=0, "invalid independent work/load scale");
    if (scale==0) { check(error==0, "zero-scale work/load is not exactly zero"); return; }
    const auto defect=std::abs(error)/scale; peak=std::max(peak,defect);
    if ((adjoint && defect>=hardGate) || defect>=verificationGate)
    {
        std::ostringstream message;
        message << std::setprecision(17) << "pressure load exceeds frozen gate: adjoint=" << adjoint
                << " defect=" << defect << " error=" << error << " scale=" << scale << " scenarios=" << scenarios;
        throw std::runtime_error(message.str());
    }
}
DictionaryMesh input(int nx, int ny, int metric)
{
    DictionaryMesh data; data.Nx=nx; data.Ny=ny;
    data.minCoord={metric==0 ? -.5 : metric==1 ? 0.0 : .25,-.25};
    data.maxCoord={metric==0 ? 1.5 : metric==1 ? 1.0 : 1.25,1.25};
    data.maxInitialIterations=0; data.maxCellNum=20000; return data;
}
void refine(QuadtreeMesh& mesh, int depth)
{
    std::vector<DictionaryMesh::RefineBox> boxes;
    const auto width=mesh.maxCoord-mesh.minCoord;
    for (int level=1; level<=depth; ++level)
    {
        const double lo=level==1 ? .25 : .375;
        DictionaryMesh::RefineBox box;
        box.min=mesh.minCoord+width*lo; box.max=mesh.minCoord+width*(1-lo); box.level=level;
        boxes.push_back(box);
    }
    if (!boxes.empty()) mesh.applyRefineBoxes(boxes);
}
FaceSegmentGraph selected(const QuadtreeMesh& mesh, int selection)
{
    FaceSegmentGraph graph;
    graph.build(mesh,[selection](const QuadtreeMesh& m,const Cell& c) {
        const auto p=c.centerCoord(), w=m.maxCoord-m.minCoord;
        const double x=(p.x-m.minCoord.x)/w.x,y=(p.y-m.minCoord.y)/w.y;
        if (selection==1) return x<.5;
        if (selection==2) return x<.25 || x>.75;
        if (selection==3) return !(x>.375 && x<.625 && y>.375 && y<.625);
        return true;
    });
    return graph;
}
struct Contribution { int p,u; long double signedArea; };
struct Reference
{
    std::vector<Contribution> contributions;
    std::vector<int> incidence;
    std::vector<long double> area;
    int mixed=0;
};
Reference reference(const QuadtreeMesh& mesh, const FaceSegmentGraph& graph, const Layout& layout, bool rz)
{
    Reference result; result.incidence.assign(layout.velocityCount(),0); result.area.assign(layout.velocityCount(),0);
    for (const auto& segment : graph.segments())
    {
        const auto& edge=segment.isHorizontal ? mesh.hEdges(segment.edgeLabel) : mesh.vEdges(segment.edgeLabel);
        const auto a=mesh.verts(edge.edgeToVert[0]).coord(), b=mesh.verts(edge.edgeToVert[1]).coord();
        const long double x0=std::min(a.x,b.x),x1=std::max(a.x,b.x),y0=std::min(a.y,b.y),y1=std::max(a.y,b.y);
        // Endpoint geometry and the independent per-cell outward normal;
        // do not call domain.face(), CartesianMeasure or the D_I matrix.
        const long double area=rz ? (segment.isHorizontal ? pi*(x1*x1-x0*x0) : 2*pi*x0*(y1-y0))
                                  : (segment.isHorizontal ? x1-x0 : y1-y0);
        const int native=segment.isHorizontal ? mesh.vEdgeIndexNum()+edge.index() : edge.index();
        const int u=layout.velocityDof(native), p=layout.pressureDof(mesh.cells(segment.K).index());
        check(u>=0 && p>=0 && std::isfinite(area) && area>=0, "invalid independent native contribution");
        const long double sign=segment.isHorizontal ? segment.n.y : segment.n.x;
        result.contributions.push_back({p,u,sign*area}); ++result.incidence[u]; result.area[u]=area;
        if (segment.L>=0 && mesh.cells(segment.L).level()>mesh.cells(segment.K).level()) ++result.mixed;
    }
    coarseSegments+=result.mixed; return result;
}
void transpose(const Load& op, const Divergence& source, const Reference& ref)
{
    const auto& d=source.matrix(); const auto& c=op.matrix();
    check(&op.divergence().matrix()==&d && op.layout().sameIdentity(source.layout()), "lost authoritative divergence handle");
    check(c.rows()==d.cols() && c.cols()==d.rows() && c.nonZeros()==d.nonZeros() && c.isCompressed(),
        "invalid compressed negative transpose dimensions/storage");
    for (Eigen::Index p=0; p<d.outerSize(); ++p)
        for (Divergence::Matrix::InnerIterator it(d,p); it; ++it)
            check(c.coeff(it.col(),it.row())==-it.value(), "C is not the exact negative transpose of D_I");
    for (Eigen::Index u=0; u<c.outerSize(); ++u)
    {
        int entries=0;
        for (Load::Matrix::InnerIterator it(c,u); it; ++it)
        { check(std::isfinite(it.value()) && it.value()!=0, "nonfinite/explicit-zero pressure coefficient"); ++entries; }
        check(ref.incidence[u]==1 || ref.incidence[u]==2, "invalid independent face incidence");
        check(entries==(ref.area[u]==0 ? 0 : ref.incidence[u]), "missing/extra boundary, subface or axis row");
    }
}
void independentLoad(const Reference& ref, const std::vector<double>& p, const std::vector<double>& actual, bool constant)
{
    std::vector<long double> exact(actual.size(),0),scale(actual.size(),0);
    for (const auto& term : ref.contributions)
    {
        const auto load=-term.signedArea*p[term.p];
        exact[term.u]+=load; scale[term.u]+=std::abs(load);
    }
    for (std::size_t u=0; u<actual.size(); ++u)
    {
        observe(actual[u]-exact[u],scale[u],loadDefect);
        if (constant && ref.incidence[u]==2) check(actual[u]==0, "constant pressure does not cancel on an internal face");
        if (ref.area[u]==0) check(actual[u]==0, "RZ axis pressure load is not exactly zero");
        if (constant && ref.incidence[u]==1 && ref.area[u]>0) check(actual[u]!=0, "uneliminated boundary load was silently removed");
    }
}
struct Retained { Load op; Values values; };
Retained run(const QuadtreeMesh& mesh, const FaceSegmentGraph& graph, bool rz, bool requireMixed)
{
    const auto layout=Layout::build(mesh,NativeLeafDomain::build(mesh,graph,rz));
    const auto divergence=Divergence::build(mesh,layout); const auto load=Load::build(divergence); const auto copy=load;
    check(&copy.matrix()==&load.matrix(), "pressure load copy lost immutable sharing");
    const auto ref=reference(mesh,graph,layout,rz); transpose(load,divergence,ref);
    if (requireMixed) check(ref.mixed>0, "fixture missed real AMR subfaces");
    std::mt19937_64 rng(810237); std::uniform_real_distribution<double> random(-1,1);
    std::vector<double> u(layout.velocityCount()),p(layout.pressureCount());
    for (int trial=0; trial<32; ++trial)
    {
        const double us=trial%3==0 ? 1e-100 : trial%3==1 ? 1 : 1e100;
        const double ps=(trial/3)%3==0 ? 1e-100 : (trial/3)%3==1 ? 1 : 1e100;
        for (auto& value : u) value=trial==0 ? 0 : trial==2 ? 1 : us*random(rng);
        for (auto& value : p) value=trial==1 ? 0 : trial==2 ? 1 : ps*random(rng);
        const Values values(layout,u,p); const auto cp=load.apply(mesh,values),du=divergence.apply(mesh,values);
        independentLoad(ref,p,cp,trial==2);
        observe(dot(u,cp)+dot(du,p),norm(u)*norm(p),adjointDefect,true);
        if (trial==3)
            check(cp==load.apply(mesh,Values(layout,std::vector<double>(u.size(),-23),p)), "pressure load depends on velocity");
        ++pairs;
    }
    if (rz && std::find(ref.area.begin(),ref.area.end(),0)!=ref.area.end())
    {
        std::vector<double> axis(u.size(),0);
        for (std::size_t i=0; i<axis.size(); ++i) if (ref.area[i]==0) axis[i]=1;
        const Values values(layout,axis,p);
        check(divergence.apply(mesh,values)==std::vector<double>(p.size(),0) && dot(axis,load.apply(mesh,values))==0,
            "zero-axis velocity mode changes divergence or pressure work");
        ++axisModes;
    }
    ++scenarios;
    std::cout << "CCH_PRESSURE," << rz << ',' << mesh.cellIndexNum() << ',' << layout.velocityCount() << ','
              << layout.pressureCount() << ',' << load.matrix().nonZeros() << ',' << ref.mixed << '\n';
    return {load,Values(layout,u,p)};
}
void grids()
{
    for (int metric : {0,1,2})
    {
        for (int n : {2,4,8,16,32})
        { QuadtreeMesh mesh(input(n,n,metric)); run(mesh,selected(mesh,0),metric!=0,false); }
        QuadtreeMesh rectangular(input(4,8,metric)); run(rectangular,selected(rectangular,0),metric!=0,false);
        for (int n : {4,8,16})
        {
            QuadtreeMesh mesh(input(n,n,metric)); refine(mesh,2);
            for (int selection : {0,1,2,3}) run(mesh,selected(mesh,selection),metric!=0,selection==0);
        }
        auto nonDyadic=input(8,6,metric);
        nonDyadic.minCoord={metric==0 ? -.37 : metric==1 ? 0.0 : .17,-.23};
        nonDyadic.maxCoord={metric==0 ? 1.19 : 1.31,1.07};
        QuadtreeMesh irregular(nonDyadic); refine(irregular,2);
        for (int selection : {0,3}) run(irregular,selected(irregular,selection),metric!=0,selection==0);
    }
}
void rebuilds()
{
    for (int metric : {0,1,2})
    {
        std::optional<Retained> previous;
        for (int depth : {0,1,2,0,1})
        {
            QuadtreeMesh mesh(input(8,8,metric)); refine(mesh,depth);
            auto current=run(mesh,selected(mesh,0),metric!=0,depth>0);
            if (previous)
            {
                rejects([&] {previous->op.apply(mesh,previous->values);});
                rejects([&] {current.op.apply(mesh,previous->values);});
            }
            previous=std::move(current);
        }
    }
}
void guardsAndLifetime()
{
    static_assert(!std::is_default_constructible_v<Load>);
    static_assert(std::is_same_v<decltype(std::declval<const Load&>().matrix()),const Load::Matrix&>);
    static_assert(std::is_same_v<decltype(std::declval<const Load&>().divergence()),const Divergence&>);
    QuadtreeMesh mesh(input(4,4,1)); const auto graph=selected(mesh,0);
    const auto layout=Layout::build(mesh,NativeLeafDomain::build(mesh,graph,false));
    const auto load=Load::build(Divergence::build(mesh,layout));
    const std::vector<double> u(layout.velocityCount(),1),p(layout.pressureCount(),1); const Values values(layout,u,p);
    const auto sameSize=Layout::build(mesh,NativeLeafDomain::build(mesh,graph,false));
    rejects([&] {load.apply(mesh,Values(sameSize,u,p));});
    auto stretchedInput=input(4,4,1); stretchedInput.maxCoord.x+=.25; QuadtreeMesh stretched(stretchedInput);
    check(stretched.cellIndexNum()==mesh.cellIndexNum(), "same-count geometry fixture differs in size");
    rejects([&] {load.apply(stretched,values);});
    auto& vertex=mesh.verts(mesh.vertIndexToLabel[0]); const auto saved=vertex.coord();
    vertex.setCoord(saved+Vector2D{.01,0}); rejects([&] {load.apply(mesh,values);});
    vertex.setCoord(saved); load.apply(mesh,values);
    const auto rz=Layout::build(mesh,NativeLeafDomain::build(mesh,graph,true));
    rejects([&] {load.apply(mesh,Values(rz,u,p));});
    const auto left=Layout::build(mesh,NativeLeafDomain::build(mesh,selected(mesh,1),false));
    FaceSegmentGraph rightGraph; rightGraph.build(mesh,[](const QuadtreeMesh&,const Cell& c) {return c.centerCoord().x>.5;});
    const auto right=Layout::build(mesh,NativeLeafDomain::build(mesh,rightGraph,false));
    check(left.velocityCount()==right.velocityCount() && left.pressureCount()==right.pressureCount(), "selection dimensions differ");
    rejects([&] {Load::build(Divergence::build(mesh,left)).apply(mesh,Values(right,std::vector<double>(right.velocityCount()),std::vector<double>(right.pressureCount())));});
    auto invalidPressure=p; invalidPressure[0]=std::numeric_limits<double>::infinity();
    rejects([&] {load.apply(mesh,Values(layout,u,invalidPressure));});
    refine(mesh,1); rejects([&] {load.apply(mesh,values);});
    auto wideInput=input(2,2,0); wideInput.minCoord={0,0}; wideInput.maxCoord={4,4}; QuadtreeMesh wide(wideInput);
    const auto wideLayout=Layout::build(wide,NativeLeafDomain::build(wide,selected(wide,0),false));
    const auto wideLoad=Load::build(Divergence::build(wide,wideLayout));
    std::vector<double> huge(wideLayout.pressureCount(),0); huge[0]=std::numeric_limits<double>::max();
    rejects([&] {wideLoad.apply(wide,Values(wideLayout,std::vector<double>(wideLayout.velocityCount()),huge));});
    check(wideLoad.apply(wide,Values(wideLayout,std::vector<double>(wideLayout.velocityCount()),std::vector<double>(huge.size())))==std::vector<double>(wideLayout.velocityCount(),0),
        "overflow poisoned immutable pressure load");
    std::optional<Retained> retained;
    {
        QuadtreeMesh temporary(input(4,4,1)); const auto domain=NativeLeafDomain::build(temporary,selected(temporary,0),true);
        const auto l=Layout::build(temporary,domain); const auto divergence=Divergence::build(temporary,l);
        retained.emplace(Retained{Load::build(divergence),Values(l,std::vector<double>(l.velocityCount()),std::vector<double>(l.pressureCount()))});
    }
    QuadtreeMesh replacement(input(4,4,1));
    check(retained->op.matrix().cols()==16 && retained->op.divergence().matrix().rows()==16 &&
        retained->op.apply(replacement,retained->values)==std::vector<double>(retained->op.layout().velocityCount(),0),
        "pressure load borrowed mesh/domain/divergence lifetime");
}
}
int main() try
{
    grids(); rebuilds(); guardsAndLifetime();
    check(scenarios==75 && pairs==2400 && rejected==32 && coarseSegments>0 && axisModes==25, "incomplete pressure load coverage");
    std::cout << std::scientific << std::setprecision(16)
              << "H03 CCH pressure load PASS scenarios=" << scenarios << " pairs=" << pairs << " rejected=" << rejected
              << " coarse_segments=" << coarseSegments << " axis_modes=" << axisModes
              << " adjoint=" << adjointDefect << " load=" << loadDefect << '\n';
    return 0;
}
catch (const std::exception& error) {std::cerr << "H03 FAIL: " << error.what() << '\n'; return 1;}
