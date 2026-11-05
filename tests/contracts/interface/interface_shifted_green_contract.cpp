#include "Core/Geometry/LevelSetMappedInterface.h"
#include "Core/Geometry/LevelSetMappedSegment.h"
#include "Core/Geometry/LevelSetHermiteGeometry.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <algorithm>
#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
using I=LevelSetMappedInterface;
using Snapshot=LevelSetHermiteGeometry::Snapshot;
using Patch=LAMP::Geometry::ActiveQuarticHermitePatch;
using Poly=std::vector<long double>;
constexpr long double pi=3.1415926535897932384626433832795L;
int cases=0,rejections=0;
void check(bool value,const std::string& why){if(!value)throw std::runtime_error(why);}
void near(long double a,long double b,long double scale,const std::string& why)
{check(std::abs(a-b)<=2e-12L*scale,why+" actual="+std::to_string(a)+" expected="+std::to_string(b));}
template<class F> void rejects(F f,const std::string& why){try{f();}catch(const std::exception& e){
    check(std::string(e.what()).find(why)!=std::string::npos,std::string("unexpected rejection: ")+e.what());++rejections;return;}
    throw std::runtime_error("missing rejection: "+why);}
Poly mul(const Poly& a,const Poly& b){Poly c(a.size()+b.size()-1,0);for(std::size_t i=0;i<a.size();++i)
    for(std::size_t j=0;j<b.size();++j)c[i+j]+=a[i]*b[j];return c;}
Poly offset(Poly a,long double v){a[0]+=v;return a;}
long double integral(const Poly& a,long double l,long double h){long double result=0;
    for(std::size_t k=0;k<a.size();++k)result+=a[k]*(std::pow(h,static_cast<int>(k+1))-std::pow(l,static_cast<int>(k+1)))/(k+1);return result;}

struct Fixture {
    std::unique_ptr<QuadtreeMesh> mesh;
    std::shared_ptr<const Snapshot> source;
    std::vector<LevelSetMappedSegment> retained;
    I::Graph graph;
    Vector2D lo,hi;
    Poly g;
    bool xGraph,radial;
    int sign;
    Fixture(bool x,bool rz,int s,bool controlled,bool thin=false):lo(thin?Vector2D{0x1p30,0x1p30}:Vector2D{1,2}),
        hi(lo+Vector2D{1,1}),g(thin?Poly{0x1p-25L,0x1p-26L,0}:Poly{.25L,.125L,.125L}),xGraph(x),radial(rz),sign(s)
    {
        DictionaryMesh d;d.Nx=d.Ny=1;d.minCoord=lo;d.maxCoord=hi;d.maxCellNum=64;d.maxInitialIterations=0;
        mesh=std::make_unique<QuadtreeMesh>(d);check(mesh->cellIndexNum()==1,"one-cell fixture");
        const auto jet=[&](Vector2D p){const double q=x?p.x-lo.x:p.y-lo.y,v=x?p.y-lo.y:p.x-lo.x;
            const double f=static_cast<double>(g[0]+g[1]*q+g[2]*q*q),df=static_cast<double>(g[1]+2*g[2]*q);
            return Patch::NodalDof{sign*(v-f),sign*(x?-df:1),sign*(x?1:-df),0};};
        Patch p;p.activeCells={1};p.nodes.resize(mesh->vertIndexNum());p.cellCenters={jet(lo+Vector2D{.5,.5}).value};
        for(int v=0;v<mesh->vertIndexNum();++v)p.nodes[v]=jet(mesh->verts(mesh->vertIndexToLabel[v]).coord());
        const auto& cell=mesh->cells(mesh->cellIndexToLabel[0]);
        for(int k=0;k<4;++k){int a=mesh->verts(cell.cellToVert[k]).index(),b=mesh->verts(cell.cellToVert[(k+1)%4]).index();
            auto pa=mesh->verts(mesh->vertIndexToLabel[a]).coord(),pb=mesh->verts(mesh->vertIndexToLabel[b]).coord();
            const bool horizontal=pa.y==pb.y;if(horizontal?pa.x>pb.x:pa.y>pb.y){std::swap(a,b);std::swap(pa,pb);}
            const auto j=jet(pa+(pb-pa)*.5);p.edgeMidpoints.push_back({{{a,b}},j.value,horizontal?j.dy:j.dx});}
        // This annular domain has no r=0 boundary requiring an even extension.
        source=std::make_shared<const Snapshot>(Snapshot::fromNativeLeaves(*mesh,p,false));
        const std::array<Vector2D,4> corners{{lo,{hi.x,lo.y},hi,{lo.x,hi.y}}};
        for(int k=0;k<4;++k)retained.emplace_back(source,corners[k],corners[(k+1)%4]);
        graph.nativeCellIndex=0;graph.coordinateIsX=x;graph.begin={x?0:1,0};graph.end={x?2:3,0};graph.controlledIntegration=controlled;
    }
    std::vector<const LevelSetMappedSegment*> cuts()const{std::vector<const LevelSetMappedSegment*> c;
        for(const auto& r:retained)c.push_back(&r);return c;}
    std::vector<I::Region> regions(bool split)const{if(!split)return {{0,{{lo,{hi.x,lo.y},hi,{lo.x,hi.y}}}}};
        const Vector2D a=xGraph?Vector2D{lo.x+.5,lo.y}:Vector2D{lo.x,lo.y+.5};
        const Vector2D b=xGraph?Vector2D{lo.x+.5,hi.y}:Vector2D{hi.x,lo.y+.5};
        if(xGraph)return {I::Region{0,{{lo,a,b,{lo.x,hi.y}}}},I::Region{1,{{a,{hi.x,lo.y},hi,b}}}};
        return {I::Region{0,{{lo,{hi.x,lo.y},b,a}}},I::Region{1,{{a,b,hi,{lo.x,hi.y}}}}};
    }
    void split(){retained.clear();for(const auto& r:regions(true))for(int k=0;k<4;++k){
        const auto a=r.vertices[k],b=r.vertices[(k+1)%4];bool present=false;
        for(const auto& cut:retained){const auto c=cut.begin(),d=cut.end();
            present=present || (a.x==c.x && a.y==c.y && b.x==d.x && b.y==d.y) || (a.x==d.x && a.y==d.y && b.x==c.x && b.y==c.y);}
        if(!present)retained.emplace_back(source,a,b);}}
    std::array<long double,2> exact(Vector2D origin,long double a,long double b)const{
        // Independent polynomial antiderivatives in q. No production roots,
        // graph quadrature or rounded physical interface positions are used.
        const Poly q{0,1},dg{g[1],2*g[2]},one{1};
        const auto r=xGraph?q:g,z=xGraph?g:q;
        const auto physicalR=offset(r,lo.x),shiftedZ=offset(z,static_cast<long double>(lo.y)-origin.y);
        const auto shiftedR=offset(r,static_cast<long double>(lo.x)-origin.x);
        const auto summedR=offset(r,static_cast<long double>(lo.x)+origin.x);
        const auto axial=mul(mul(radial?physicalR:one,shiftedZ),xGraph?one:dg);
        const auto radialP=mul(mul(shiftedR,radial?summedR:one),xGraph?dg:one);
        return {{sign*(xGraph?1:-1)*(radial?2*pi:1)*integral(axial,a,b),
            sign*(xGraph?-1:1)*(radial?pi:1)*integral(radialP,a,b)}};
    }
};
void same(const I::QuarticResult& a,const I::QuarticResult& b){
    check(a.arcs.size()==b.arcs.size() && a.sourceRootIncidences==b.sourceRootIncidences && a.retainedCuts==b.retainedCuts,"zero origin changed structure");
    for(std::size_t k=0;k<a.arcs.size();++k){const auto& x=a.arcs[k];const auto& y=b.arcs[k];
        check(x.regionId==y.regionId && x.fromEnd==y.fromEnd && x.offset==y.offset && x.width==y.width &&
            x.axialMoment==y.axialMoment && x.radialMoment==y.radialMoment && x.momentErrorEstimate==y.momentErrorEstimate &&
            x.sensitivity==y.sensitivity && x.minTransverseDerivative==y.minTransverseDerivative &&
            x.maxScalarResidual==y.maxScalarResidual && x.quadraturePoints==y.quadraturePoints &&
            x.nativeCellIndex==y.nativeCellIndex && x.requiresNativeDofConstraints==y.requiresNativeDofConstraints &&
            x.independentSensitivity.size()==y.independentSensitivity.size(),"zero origin changed arc fields");
        for(std::size_t j=0;j<x.independentSensitivity.size();++j){const auto& u=x.independentSensitivity[j];const auto& v=y.independentSensitivity[j];
            check(u.kind==v.kind && u.index==v.index && u.component==v.component && u.coefficient==v.coefficient,"zero origin changed constrained sensitivity");}}
}
void verify(Fixture& f,Vector2D origin,bool split){
    const auto cuts=f.cuts();const auto regions=f.regions(split);
    const auto old=I::integrateQuartic(*f.source,f.graph,cuts,regions,f.radial);
    same(old,I::integrateQuartic(*f.source,f.graph,cuts,regions,f.radial,{0,0}));
    const auto shifted=I::integrateQuartic(*f.source,f.graph,cuts,regions,f.radial,origin);
    check(shifted.arcs.size()==old.arcs.size(),"translated potential changed cut count");
    long double totalA=0,totalR=0;
    for(std::size_t k=0;k<shifted.arcs.size();++k){const auto& a=shifted.arcs[k];const auto& b=old.arcs[k];
        check(a.regionId==b.regionId && a.fromEnd==b.fromEnd && a.offset==b.offset && a.width==b.width && a.nativeCellIndex==b.nativeCellIndex,"translated potential changed retained arc");
        const auto expected=f.exact(origin,split?(a.regionId==0?0:.5L):0,split?(a.regionId==0?.5L:1):1);
        near(a.axialMoment,expected[0],std::max(std::abs(expected[0]),1e-30L),"axial polynomial integral");
        near(a.radialMoment,expected[1],std::max(std::abs(expected[1]),1e-30L),"radial polynomial integral");
        for(std::size_t j=0;j<25;++j)near(a.sensitivity[j],b.sensitivity[j],std::max(1.L,std::abs(b.sensitivity[j])),"volume sensitivity");
        check(a.independentSensitivity.size()==b.independentSensitivity.size(),"sensitivity identity count");
        for(std::size_t j=0;j<a.independentSensitivity.size();++j){const auto& u=a.independentSensitivity[j];const auto& v=b.independentSensitivity[j];
            check(u.kind==v.kind && u.index==v.index && u.component==v.component,"sensitivity identity changed");
            near(u.coefficient,v.coefficient,std::max(1.L,std::abs(v.coefficient)),"constrained sensitivity");}
        totalA+=a.axialMoment;totalR+=a.radialMoment;
    }
    const auto all=f.exact(origin,0,1);near(totalA,all[0],std::max(std::abs(all[0]),1e-30L),"partition axial sum");
    near(totalR,all[1],std::max(std::abs(all[1]),1e-30L),"partition radial sum");
    // Close the free curve with the actual rectangle boundary, using the
    // same potential origin. Both independent Green directions give volume.
    const long double g0=f.g[0],g1=f.g[0]+f.g[1]+f.g[2],r0=f.lo.x;
    const long double capacity=f.radial?pi*(2*r0+1):1;
    const auto radialPotential=[&](long double x){return (static_cast<long double>(f.lo.x)-origin.x+x)*
        (f.radial?pi*(static_cast<long double>(f.lo.x)+origin.x+x):1);};
    const auto slice=[&](long double g){return f.radial?pi*g*(2*r0+g):g;};
    const long double bottom=static_cast<long double>(f.lo.y)-origin.y,top=bottom+1;
    const long double belowA=f.xGraph?-bottom*capacity:top*slice(g1)-bottom*slice(g0);
    const long double belowR=f.xGraph?radialPotential(1)*g1-radialPotential(0)*g0:-radialPotential(0);
    const long double below=f.xGraph?(f.radial?2*pi*integral(mul(f.g,Poly{r0,1}),0,1):integral(f.g,0,1)):
        (f.radial?pi*integral(mul(f.g,offset(f.g,2*r0)),0,1):integral(f.g,0,1));
    const long double volume=f.sign>0?below:capacity-below;
    near(totalA+(f.sign>0?belowA:capacity-belowA),volume,capacity,"closed axial volume");
    near(totalR+(f.sign>0?belowR:capacity-belowR),volume,capacity,"closed radial volume");++cases;
}
}
int main(){try{
    for(bool radial:{false,true})for(bool x:{false,true})for(int sign:{-1,1})for(bool controlled:{false,true}){
        Fixture f(x,radial,sign,controlled);
        for(bool split:{false,true}){if(split)f.split();for(Vector2D o:{Vector2D{0,0},f.lo,f.lo+Vector2D{.375,.625},Vector2D{-2,3}})verify(f,o,split);}
    }
    for(bool radial:{false,true})for(bool x:{false,true})for(int sign:{-1,1}){
        Fixture thin(x,radial,sign,true,true);verify(thin,thin.lo,false);thin.split();verify(thin,thin.lo,true);
        check(thin.lo.x+static_cast<double>(thin.g[0])==thin.lo.x,"thin fixture must lose displacement in global double");
    }
    Fixture f(true,false,1,true);const auto cuts=f.cuts();const auto regions=f.regions(false);
    for(double bad:{std::numeric_limits<double>::quiet_NaN(),std::numeric_limits<double>::infinity(),-std::numeric_limits<double>::infinity()})
        for(Vector2D o:{Vector2D{bad,0},Vector2D{0,bad}})rejects([&]{I::integrateQuartic(*f.source,f.graph,cuts,regions,false,o);},"nonfinite Green-potential origin");
    Fixture other(true,false,1,true);auto wrong=cuts;wrong[0]=&other.retained[0];
    rejects([&]{I::integrateQuartic(*f.source,f.graph,wrong,regions,false,f.lo);},"different snapshot");
    wrong=cuts;wrong[0]=nullptr;rejects([&]{I::integrateQuartic(*f.source,f.graph,wrong,regions,false,f.lo);},"null retained face cut");
    verify(f,f.lo,false);
    std::cout<<"I12a5 shifted Green PASS cases="<<cases<<" rejections="<<rejections<<'\n';return 0;
}catch(const std::exception& e){std::cerr<<"I12a5 shifted Green FAIL: "<<e.what()<<'\n';return 1;}}
