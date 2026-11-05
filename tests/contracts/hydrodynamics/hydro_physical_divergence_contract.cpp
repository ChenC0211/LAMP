#include "Modules/Hydrodynamics/HydroPhysicalDivergence.h" // Public header stands alone.
#include "Core/IO/DictionaryMesh.h"
#include "Core/Mesh/FaceSegmentGraph.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <string>
#include <type_traits>

namespace
{
using namespace LAMP::Hydrodynamics;
using Axis=HydroDofLayout::NormalAxis;
const long double pi=std::acos(-1.0L);
constexpr long double tolerance=1e-13L;
int metricCases=0,fieldCases=0,nearAxis=0,mmsCases=0,reflections=0,rejected=0;
long double metricDefect=0,analyticDefect=0,globalDefect=0,axisDefect=0;
void check(bool ok,const char* reason) {if (!ok) throw std::runtime_error(reason);}
void observe(long double error,long double scale,long double& peak) {
    check(std::isfinite(error) && std::isfinite(scale) && scale>0,"invalid independent scale");
    const long double defect=std::abs(error)/scale;peak=std::max(peak,defect);
    check(defect<=tolerance,"independent true-measure/regularity gate failed");
}
template<class F> void rejects(F action,const char* stage=nullptr) {
    try {action();} catch (const std::invalid_argument& error) {
        if (stage) check(std::string(error.what()).find(stage)!=std::string::npos,"wrong rejection stage");
        ++rejected;return;
    }
    throw std::runtime_error("invalid physical-divergence/axis input accepted");
}
DictionaryMesh input(int n,int metric) {
    DictionaryMesh data;data.Nx=n;data.Ny=n;data.minCoord={metric==0 ? -.5 : metric==1 ? 0 : .25,-.25};
    data.maxCoord={metric==0 ? 1.5 : metric==1 ? 1 : 1.25,1.25};data.maxInitialIterations=0;data.maxCellNum=20000;return data;
}
void refine(QuadtreeMesh& mesh,int depth) {
    std::vector<DictionaryMesh::RefineBox> boxes;const auto w=mesh.maxCoord-mesh.minCoord;
    for (int level=1;level<=depth;++level) {
        DictionaryMesh::RefineBox b;const double upper=level==1 ? .5 : .25,lowerZ=level==1 ? .25 : .375;
        b.min={mesh.minCoord.x,mesh.minCoord.y+lowerZ*w.y};b.max={mesh.minCoord.x+upper*w.x,mesh.minCoord.y+(1-lowerZ)*w.y};b.level=level;boxes.push_back(b);
    }
    if (!boxes.empty()) mesh.applyRefineBoxes(boxes);
}
FaceSegmentGraph selected(const QuadtreeMesh& mesh,int selection) {
    FaceSegmentGraph graph;graph.build(mesh,[selection](const QuadtreeMesh& m,const Cell& c) {
        const auto p=c.centerCoord(),w=m.maxCoord-m.minCoord;const double x=(p.x-m.minCoord.x)/w.x,y=(p.y-m.minCoord.y)/w.y;
        if (selection==1) return x<.5;
        if (selection==2) return x<.25 || x>.75;
        if (selection==3) return !(x>.375 && x<.625 && y>.375 && y<.625);
        return true;
    });return graph;
}
struct Face {long double a=0,b=0,c=0,d=0,area=0;bool horizontal=false;int incidence=0;};
struct CellBounds {long double a,b,c,d,volume;};
struct Term {int p,u;long double signedArea;};
struct Reference {std::vector<Face> faces;std::vector<CellBounds> cells;std::vector<Term> terms;};
Reference reference(const QuadtreeMesh& mesh,const FaceSegmentGraph& graph,const HydroDofLayout& layout,bool rz) {
    Reference ref;ref.faces.resize(layout.velocityCount());
    for (int p=0;p<layout.pressureCount();++p) {
        const auto& cell=mesh.cells(mesh.cellIndexToLabel[layout.nativeCell(p)]);const auto x=cell.centerCoord(),h=cell.size();
        CellBounds b{x.x-h.x/2,x.x+h.x/2,x.y-h.y/2,x.y+h.y/2,0};
        b.volume=(rz ? pi*(b.b*b.b-b.a*b.a) : b.b-b.a)*(b.d-b.c);ref.cells.push_back(b);
    }
    for (const auto& s : graph.segments()) {
        const auto& e=s.isHorizontal ? mesh.hEdges(s.edgeLabel) : mesh.vEdges(s.edgeLabel);
        const auto a=mesh.verts(e.edgeToVert[0]).coord(),b=mesh.verts(e.edgeToVert[1]).coord();
        const int u=layout.velocityDof(e.index()+(s.isHorizontal ? mesh.vEdgeIndexNum() : 0));auto& f=ref.faces[u];
        f.a=std::min(a.x,b.x);f.b=std::max(a.x,b.x);f.c=std::min(a.y,b.y);f.d=std::max(a.y,b.y);f.horizontal=s.isHorizontal;++f.incidence;
        f.area=rz ? (f.horizontal ? pi*(f.b*f.b-f.a*f.a) : 2*pi*f.a*(f.d-f.c)) : (f.horizontal ? f.b-f.a : f.d-f.c);
        ref.terms.push_back({layout.pressureDof(mesh.cells(s.K).index()),u,f.area*(f.horizontal ? s.n.y : s.n.x)});
    }
    return ref;
}
long double squareMean(long double a,long double b,bool weighted=false) {
    return weighted ? (a*a+b*b)/2 : (a*a+a*b+b*b)/3;
}
// Regular polynomial: radial odd, axial even. Integrals come from endpoints.
long double radial(long double r,long double z) {return .7L*r+.2L*r*r*r-.15L*r*z*z;}
long double axial(long double r,long double z) {return .45L+.3L*r*r*z-.25L*z*z*z;}
long double pressure(long double r,long double z) {return r*r+z*z+std::cos(pi*r);}
std::vector<double> polynomial(const Reference& ref,bool rz,int field) {
    std::vector<double> u;
    for (const auto& f : ref.faces) {
        const long double r=(f.a+f.b)/2,z=(f.c+f.d)/2;
        const long double r2=squareMean(f.a,f.b,rz),z2=squareMean(f.c,f.d);
        long double v=0;
        if (field==0) v=f.horizontal ? .45L+.3L*r2*z-.25L*z*z*z : .7L*r+.2L*r*r*r-.15L*r*z2;
        if (field==1) v=f.horizontal ? (rz ? 1 : .5L)*z*z : -r*z;
        if (field==2) v=f.horizontal ? .45L : .7L*r;
        u.push_back(static_cast<double>(v));
    }
    return u;
}
void metric(int n,int measure,int depth,int selection) {
    QuadtreeMesh mesh(input(n,measure));refine(mesh,depth);const auto graph=selected(mesh,selection);
    const bool rz=measure!=0;const auto layout=HydroDofLayout::build(mesh,NativeLeafDomain::build(mesh,graph,rz));
    const auto di=HydroIntegralDivergence::build(mesh,layout);const auto op=HydroPhysicalDivergence::build(mesh,di),copy=op;
    const auto ref=reference(mesh,graph,layout,rz);
    check(&op.divergence().matrix()==&di.matrix(),"physical divergence lost authoritative integral source");
    check(&op.cellVolumes()==&di.layout().domain().cellVolumes(),"physical divergence duplicated cell measures");
    check(&copy.faceMeasures()==&op.faceMeasures() && &copy.axisVelocityDofs()==&op.axisVelocityDofs(),"copy lost shared immutable metrics");
    int axes=0;
    for (std::size_t u=0;u<ref.faces.size();++u) {
        const auto area=ref.faces[u].area;
        if (area==0) {check(op.faceMeasures()[u]==0,"axis has artificial area");++axes;}
        else observe(op.faceMeasures()[u]-area,area,metricDefect);
    }
    check(axes==static_cast<int>(op.axisVelocityDofs().size()),"wrong axis mask");
    for (std::size_t p=0;p<ref.cells.size();++p) observe(op.cellVolumes()[p]-ref.cells[p].volume,ref.cells[p].volume,metricDefect);
    for (int field : {0,1,2}) {
        const auto u=polynomial(ref,rz,field);const auto result=op.apply(mesh,HydroNativeValues(layout,u,std::vector<double>(layout.pressureCount())));
        long double integral=0,boundary=0,totalScale=0;
        for (std::size_t p=0;p<ref.cells.size();++p) {
            const auto b=ref.cells[p];const long double exact=field==0 ? (rz ? 1.4L : .7L)+(rz ? 1.1L : .9L)*squareMean(b.a,b.b,rz)-(rz ? 1.05L : .9L)*squareMean(b.c,b.d)
                : field==1 ? 0 : rz ? 1.4L : .7L;
            observe(result[p]-exact,std::max(1.0L,std::abs(exact)),analyticDefect);integral+=ref.cells[p].volume*result[p];
        }
        for (const auto& term : ref.terms) {
            const auto flow=term.signedArea*u[term.u];totalScale+=std::abs(flow);
            if (ref.faces[term.u].incidence==1) boundary+=flow;
        }
        observe(integral-boundary,std::max(1.0L,totalScale),globalDefect);++fieldCases;
    }
    if (axes) for (long double r : {.05L,.2L,.7L}) for (long double z : {-.1L,.3L,.8L}) {
        observe(op.reflectAxisNormalVelocity(Axis::First,static_cast<double>(radial(r,z)))-radial(-r,z),1,axisDefect);
        const double uz=static_cast<double>(axial(r,z)),ghost=op.reflectAxisNormalVelocity(Axis::Second,uz);
        observe(ghost-axial(-r,z),1,axisDefect);check((uz-ghost)/(2*r)==0,"axial even reflection does not enforce axis Neumann");
        observe(op.reflectAxisPressure(static_cast<double>(pressure(r,z)))-pressure(-r,z),1,axisDefect);++reflections;
    }
    ++metricCases;
}
void nearAxisLimit(int n,double radius) {
    auto data=input(n,1);data.minCoord={0,0};data.maxCoord={radius,1};QuadtreeMesh mesh(data);const auto graph=selected(mesh,0);
    const auto layout=HydroDofLayout::build(mesh,NativeLeafDomain::build(mesh,graph,true));
    const auto op=HydroPhysicalDivergence::build(mesh,HydroIntegralDivergence::build(mesh,layout));
    const auto ref=reference(mesh,graph,layout,true);const auto u=polynomial(ref,true,2);
    const auto result=op.apply(mesh,HydroNativeValues(layout,u,std::vector<double>(layout.pressureCount())));
    for (double d : result) observe(d-1.4L,1.4L,axisDefect);++nearAxis;
}
long double cosinePrimitive(long double r) {return r*std::sin(pi*r)/pi+std::cos(pi*r)/(pi*pi);}
double mms(int n) {
    auto data=input(n,1);data.minCoord={0,0};data.maxCoord={1,1};QuadtreeMesh mesh(data);const auto graph=selected(mesh,0);
    const auto layout=HydroDofLayout::build(mesh,NativeLeafDomain::build(mesh,graph,true));
    const auto op=HydroPhysicalDivergence::build(mesh,HydroIntegralDivergence::build(mesh,layout));const auto ref=reference(mesh,graph,layout,true);
    std::vector<double> u;
    for (const auto& f : ref.faces) {const auto r=(f.a+f.b)/2,z=(f.c+f.d)/2;u.push_back(static_cast<double>((f.horizontal ? std::cos(pi*r) : r)*std::sin(pi*z)));}
    const auto d=op.apply(mesh,HydroNativeValues(layout,u,std::vector<double>(layout.pressureCount())));long double error=0,scale=0;
    for (std::size_t p=0;p<d.size();++p) {
        const auto b=ref.cells[p];const auto z=(b.c+b.d)/2,hz=b.d-b.c;
        const auto sinc=std::sin(pi*hz/2)/(pi*hz/2);
        const auto weightedCos=2*(cosinePrimitive(b.b)-cosinePrimitive(b.a))/(b.b*b.b-b.a*b.a);
        const auto exact=sinc*(2*std::sin(pi*z)+pi*weightedCos*std::cos(pi*z));
        error+=b.volume*(d[p]-exact)*(d[p]-exact);scale+=b.volume*exact*exact;
    }
    const double result=static_cast<double>(std::sqrt(error/scale));check(std::isfinite(result) && result>0,"invalid continuous RZ error");
    ++mmsCases;std::cout << "HYDRO_RZ_DIVERGENCE_MMS," << n << ',' << std::setprecision(17) << result << '\n';return result;
}
void guards() {
    static_assert(!std::is_default_constructible_v<HydroPhysicalDivergence>);
    static_assert(std::is_same_v<decltype(std::declval<const HydroPhysicalDivergence&>().faceMeasures()),const std::vector<double>&>);
    QuadtreeMesh mesh(input(4,1));const auto graph=selected(mesh,0);const auto layout=HydroDofLayout::build(mesh,NativeLeafDomain::build(mesh,graph,true));
    const auto di=HydroIntegralDivergence::build(mesh,layout);const auto op=HydroPhysicalDivergence::build(mesh,di);
    std::vector<double> zero(layout.velocityCount()),p(layout.pressureCount());const HydroNativeValues values(layout,zero,p);
    const auto accepted=op.apply(mesh,values);check(accepted==p,"zero field is not exactly zero divergence");
    const auto fresh=HydroDofLayout::build(mesh,NativeLeafDomain::build(mesh,graph,true));rejects([&]{op.apply(mesh,HydroNativeValues(fresh,zero,p));});
    auto bad=zero;bad[op.axisVelocityDofs().front()]=1;const HydroNativeValues blind(layout,bad,p);
    check(di.apply(mesh,blind)==p,"fixture is not a zero-area blind mode");rejects([&]{op.apply(mesh,blind);},"nonzero radial axis velocity");
    bad[op.axisVelocityDofs().front()]=std::numeric_limits<double>::denorm_min();rejects([&]{op.apply(mesh,HydroNativeValues(layout,bad,p));},"nonzero radial axis velocity");
    auto stretch=input(4,1);stretch.maxCoord.x=1.25;QuadtreeMesh stretched(stretch);rejects([&]{HydroPhysicalDivergence::build(stretched,di);});rejects([&]{op.apply(stretched,values);});
    auto& vertex=mesh.verts(mesh.vertIndexToLabel[0]);const auto saved=vertex.coord();vertex.setCoord(saved+Vector2D{.01,0});rejects([&]{op.apply(mesh,values);});vertex.setCoord(saved);
    for (int metric : {0,2}) {
        QuadtreeMesh other(input(4,metric));const auto l=HydroDofLayout::build(other,NativeLeafDomain::build(other,selected(other,0),metric!=0));
        const auto noAxis=HydroPhysicalDivergence::build(other,HydroIntegralDivergence::build(other,l));
        rejects([&]{noAxis.reflectAxisNormalVelocity(Axis::First,1);});rejects([&]{noAxis.reflectAxisPressure(1);});
    }
    rejects([&]{op.reflectAxisNormalVelocity(Axis::Second,std::numeric_limits<double>::infinity());});
    rejects([&]{op.reflectAxisPressure(std::numeric_limits<double>::quiet_NaN());});
    rejects([&]{op.reflectAxisNormalVelocity(static_cast<Axis>(2),1);});
    auto tiny=input(4,1);tiny.minCoord={0,0};tiny.maxCoord={.01,.01};QuadtreeMesh small(tiny);
    const auto sl=HydroDofLayout::build(small,NativeLeafDomain::build(small,selected(small,0),true));const auto sd=HydroIntegralDivergence::build(small,sl);
    const auto so=HydroPhysicalDivergence::build(small,sd);std::vector<double> huge(sl.velocityCount());
    for (int u=0;u<sl.velocityCount();++u) if (sl.normalAxis(u)==Axis::First && so.faceMeasures()[u]>0) {huge[u]=std::numeric_limits<double>::max();break;}
    const HydroNativeValues large(sl,huge,std::vector<double>(sl.pressureCount()));const auto finite=sd.apply(small,large);
    for (double x : finite) check(std::isfinite(x),"overflow fixture already failed integral divergence");
    rejects([&]{so.apply(small,large);},"nonfinite volume-normalized divergence");
    refine(mesh,1);rejects([&]{HydroPhysicalDivergence::build(mesh,di);});rejects([&]{op.apply(mesh,values);});
    QuadtreeMesh equivalent(input(4,1));check(op.apply(equivalent,values)==accepted,"failed operation changed source/result");
    std::optional<HydroPhysicalDivergence> retained;std::optional<HydroNativeValues> retainedValues;
    {QuadtreeMesh temporary(input(4,1));const auto l=HydroDofLayout::build(temporary,NativeLeafDomain::build(temporary,selected(temporary,0),true));
     retained=HydroPhysicalDivergence::build(temporary,HydroIntegralDivergence::build(temporary,l));retainedValues.emplace(l,std::vector<double>(l.velocityCount()),std::vector<double>(l.pressureCount()));}
    check(retained->apply(equivalent,*retainedValues)==accepted,"physical divergence borrowed mesh lifetime");
}
}
int main() try {
    for (int measure : {0,1,2}) for (int n : {4,8}) for (int depth : {0,1,2}) for (int selection : {0,1,2,3}) metric(n,measure,depth,selection);
    for (int n : {4,8,16}) for (double radius : {1.0,1e-3,1e-6}) nearAxisLimit(n,radius);
    const double coarse=mms(4),middle=mms(8),fine=mms(16),order=std::log(middle/fine)/std::log(2.0);
    check(coarse>middle && middle>fine && order>=1.8,"RZ divergence MMS misses frozen convergence gate");guards();
    check(metricCases==72 && fieldCases==216 && nearAxis==9 && mmsCases==3 && reflections==216 && rejected==16,"incomplete H05 coverage");
    std::cout << std::scientific << std::setprecision(16) << "H05 HYDRO physical divergence PASS metrics=" << metricCases << " fields=" << fieldCases << " near_axis=" << nearAxis
        << " mms=" << mmsCases << " reflections=" << reflections << " rejected=" << rejected << " measure=" << metricDefect << " analytic=" << analyticDefect
        << " global=" << globalDefect << " axis=" << axisDefect << " order=" << order << '\n';return 0;
} catch (const std::exception& error) {std::cerr << "H05 FAIL: " << error.what() << '\n';return 1;}
