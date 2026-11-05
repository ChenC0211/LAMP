#include "Core/Geometry/QuadraticSurfaceMeasure.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>
namespace Q=LAMP::Geometry::QuadraticCellCut;
namespace S=LAMP::Geometry::QuadraticSurfaceMeasure;
void check(bool ok,const char* message){if(!ok)throw std::runtime_error(message);}
void near(long double a,long double b,long double scale,const char*message){check(std::isfinite(a)&&std::isfinite(b)&&std::abs(a-b)<=1e-12L*scale,message);}
template<class F>void rejects(F&&f,const char*message){bool bad=false;try{f();}catch(const std::exception&){bad=true;}check(bad,message);}
struct Node{long double s,w;};
std::vector<Node> rule(){std::vector<Node> nodes;constexpr int n=96;
    for(int i=0;i<n/2;++i){auto x=std::cos(std::acos(-1.L)*(i+.75L)/(n+.5L));long double d=0;
        for(int iteration=0;iteration<40;++iteration){long double p0=1,p1=x;for(int k=2;k<=n;++k){const auto p=((2*k-1)*x*p1-(k-1)*p0)/k;p0=p1;p1=p;}
            d=n*(x*p1-p0)/(x*x-1);const auto dx=p1/d;x-=dx;if(std::abs(dx)<=4*std::numeric_limits<long double>::epsilon())break;check(iteration<39,"independent moment reference convergence");}
        long double p0=1,p1=x;for(int k=2;k<=n;++k){const auto p=((2*k-1)*x*p1-(k-1)*p0)/k;p0=p1;p1=p;}d=n*(x*p1-p0)/(x*x-1);
        const auto w=1/((1-x*x)*d*d);nodes.push_back({(1-x)/2,w});nodes.push_back({(1+x)/2,w});}return nodes;}
S::MomentValues moments(const S::MomentPatch&p,const std::vector<Node>&nodes){S::MomentValues value;
    for(const auto&n:nodes){const auto m=(1-n.s)*p.m0+n.s*p.m1,g=std::hypot(1.L,m);long double power=1;
        for(int j=0;j<3;++j){value.moments[j]+=n.w*power*g;value.firstDerivative[j]+=n.w*power*m/g*(1-n.s);value.lastDerivative[j]+=n.w*power*m/g*n.s;power*=n.s;}}return value;}
S::Sensitivity inspect(const Q::Rectangle&source,const Q::Reconstruction&graph,const Q::Rectangle&region,bool rz,const std::vector<Node>&nodes){
    const auto plan=S::MomentPlan::prepare(source,graph,region,rz);const auto original=S::differentiate(source,graph,region,rz);
    check(plan.axisymmetric()==rz&&plan.patches().size()==original.surface.patches.size(),"moment plan preserves original patch inventory");
    std::vector<S::MomentValues> values;for(std::size_t i=0;i<plan.patches().size();++i){const auto&p=plan.patches()[i];const auto&q=original.surface.patches[i];
        check(p.lower==q.lower&&p.upper==q.upper&&p.begin.x==q.begin.x&&p.begin.y==q.begin.y&&p.end.x==q.end.x&&p.end.y==q.end.y&&p.beginBoundary==q.beginBoundary&&p.endBoundary==q.endBoundary,"same actual clip roots/ports/ownership metadata");values.push_back(moments(p,nodes));}
    const auto applied=plan.apply(values);const auto scale=std::max(1.L,original.surface.integral.physicalMeasure);
    near(applied.surface.integral.length,original.surface.integral.length,scale,"moment length versus original S");
    near(applied.surface.integral.physicalMeasure,original.surface.integral.physicalMeasure,scale,"moment physical area versus original S");
    for(int k=0;k<3;++k)near(applied.derivative[k],original.derivative[k],std::max(scale,std::abs(original.derivative[k])),"moment physical moving-endpoint sensitivity versus original S");
    check(applied.surface.integral.quadratureErrorEstimate==0,"plan apply performs no quadrature; not an interpolation certificate");
    auto wrong=values;wrong.push_back({});rejects([&]{(void)plan.apply(wrong);},"moment-count mismatch accepted");
    if(!values.empty()){auto bad=values;bad[0].firstDerivative[0]=std::numeric_limits<long double>::infinity();rejects([&]{(void)plan.apply(bad);},"nonfinite moment partial accepted");
        bad=values;bad[0].moments[0]=-1;rejects([&]{(void)plan.apply(bad);},"negative model length accepted");
        if(rz)for(int j=0;j<3;++j)if(plan.patches()[0].radius[j]<0){bad=values;bad[0].moments={1,1,1};bad[0].moments[j]=1e12L;
            rejects([&]{(void)plan.apply(bad);},"negative model radial area was clipped/accepted");break;}}
    return applied;
}
Q::Rectangle rectangle(const QuadtreeMesh&mesh,int index){const auto&c=mesh.cells(mesh.cellIndexToLabel[index]);const auto a=mesh.verts(c.cellToVert[0]).coord(),b=mesh.verts(c.cellToVert[2]).coord();return {{a.x,a.y},{b.x,b.y}};}
int main(){try{
    const auto nodes=rule();DictionaryMesh input;input.Nx=input.Ny=4;input.maxCellNum=4096;input.maxInitialIterations=0;input.maxLevel_local={1};QuadtreeMesh mesh(input);
    Q::Rectangle source{{.25,.25},{.5,.5}};unsigned coarse=0;for(int i=0;i<mesh.cellIndexNum();++i){const auto b=rectangle(mesh,i);coarse+=b.lower.x==source.lower.x&&b.lower.y==source.lower.y&&b.upper.x==source.upper.x&&b.upper.y==source.upper.y;}check(coarse==1,"actual coarse donor rectangle");
    Q::Shape shape;shape.frameNormal=Vector2D{.27488060480878,.9614783685034};const auto graph=Q::fromFixedPorts(source,shape,{.25,.34},{.5,.4},.012345L);
    DictionaryMesh::RefineBox refined;refined.level=1;refined.min={.25,.25};refined.max={.5,.5};mesh.applyRefineBoxes({refined});
    std::vector<Q::Rectangle> clips;for(int i=0;i<mesh.cellIndexNum();++i){const auto b=rectangle(mesh,i);if(b.lower.x>=source.lower.x&&b.lower.y>=source.lower.y&&b.upper.x<=source.upper.x&&b.upper.y<=source.upper.y)clips.push_back(b);}
    check(clips.size()==4,"actual 2:1 AMR four-child clip inventory");
    for(bool rz:{false,true}){const auto whole=inspect(source,graph,source,rz,nodes);long double area=0;std::array<long double,3> derivative{};
        for(const auto&clip:clips){const auto part=inspect(source,graph,clip,rz,nodes);area+=part.surface.integral.physicalMeasure;for(int k=0;k<3;++k)derivative[k]+=part.derivative[k];}
        near(area,whole.surface.integral.physicalMeasure,1,"actual AMR clipped moment additivity");for(int k=0;k<3;++k)near(derivative[k],whole.derivative[k],std::max(1.L,std::abs(whole.derivative[k])),"actual AMR moving-boundary derivative additivity");}
    const Q::Rectangle axisBox{{0,.25},{.25,.5}};const auto axis=Q::fromAxisPorts(axisBox,1,{0,.34},{.25,.4});(void)inspect(axisBox,axis,axisBox,true,nodes);
    auto stale=graph;stale.offset+=.1L;rejects([&]{(void)S::MomentPlan::prepare(source,stale,source,false);},"stale encoded view accepted");
    auto legacy=graph;legacy.exact.reset();rejects([&]{(void)S::MomentPlan::prepare(source,legacy,source,false);},"legacy moment chart silently accepted");
    const auto corner=Q::fromFixedPorts(source,Q::Shape{}, {.25,.25},{.5,.5},0);rejects([&]{(void)S::MomentPlan::prepare(source,corner,source,false);},"corner/coincident physical derivative silently accepted");
    rejects([&]{(void)S::MomentPlan::prepare(source,graph,source,false,32);},"invalid clip ownership accepted");
    std::cout<<"quadratic moment plan contract PASS; same-polynomial moments, moving endpoints and actual AMR clips\n";return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
