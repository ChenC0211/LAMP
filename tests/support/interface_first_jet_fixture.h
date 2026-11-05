#pragma once
#include "Modules/Interface/HermiteGalsPublication.h"
#include "Modules/Interface/HermiteGalsMeshTransfer.h"
#include "Modules/Interface/HermiteGalsFirstJetStep.h"
#include "Modules/Interface/HermiteGalsPhaseFlux.h"
#include "Modules/Interface/HermiteGalsTemporalPhaseFlux.h"
#include "Modules/Interface/HermiteGalsPhaseRegion.h"
#include "Modules/Interface/HermiteGalsInventoryRegion.h"
#include "Modules/Interface/Detail/HermiteGalsVolumeRoundoff.h"
#include "Modules/Interface/HermiteGalsVertexTransport.h"
#include "Modules/Interface/HermiteGalsEdgeTransport.h"
#include "Modules/Interface/HermiteGalsReflux.h"
#include "Core/Geometry/LevelSetMappedSegment.h"
#include "Core/Geometry/UnitPolynomialRoots.h"

#include "Modules/Interface/HermiteGalsPhaseVolume.h"
#include "Modules/Interface/InterfaceCellMeasure.h"
#include "Core/Geometry/ActiveQuarticHermitePatch.h"
#include "Core/Geometry/LevelSetHermiteGeometry.h"
#include "Core/Mesh/FaceSegmentGraph.h"
#include "Core/Mesh/NativeLeafDomain.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <algorithm>
#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <type_traits>

namespace
{
using namespace LAMP::Interface;
using Patch=LAMP::Geometry::ActiveQuarticHermitePatch;
using Jet=Patch::NodalDof;
using Profile=std::function<Jet(Vector2D)>;
using Key=std::array<int,2>;
using Position=HermiteGalsTimeLayer::Position;
constexpr long double pi=3.1415926535897932384626433832795L;
void check(bool ok,const std::string& why){if(!ok)throw std::runtime_error(why);}
long double close(long double a,long double b,long double tolerance,const std::string& why)
{const auto error=std::abs(a-b);check(std::isfinite(a) && error<=tolerance,why+" error="+std::to_string(static_cast<double>(error)));return error;}
template<class F> void rejects(F f,const std::string& expected)
{try{f();}catch(const std::exception& e){check(std::string(e.what()).find(expected)!=std::string::npos,std::string("unexpected rejection: ")+e.what());std::cout<<"REJECT,"<<expected<<'\n';return;}throw std::runtime_error("missing rejection: "+expected);}
std::unique_ptr<QuadtreeMesh> makeMesh(int n,int amr=0)
{DictionaryMesh d;d.Nx=d.Ny=n;d.maxCellNum=16384;d.maxInitialIterations=0;auto mesh=std::make_unique<QuadtreeMesh>(d);
    if(amr>1){DictionaryMesh::RefineBox b;b.level=1;b.min={.125,.125};b.max={.875,.875};mesh->applyRefineBoxes({b});}
    if(amr){DictionaryMesh::RefineBox b;b.level=amr;b.min={.25,.25};b.max={.75,.75};mesh->applyRefineBoxes({b});}return mesh;}
Vector2D point(const QuadtreeMesh& mesh,int vertex){return mesh.verts(mesh.vertIndexToLabel[vertex]).coord();}
NativeLeafDomain domainFor(const QuadtreeMesh& mesh,bool radial)
{FaceSegmentGraph graph;graph.build(mesh);return NativeLeafDomain::build(mesh,graph,radial);}
bool fineSide(const QuadtreeMesh& m,Key key)
{const auto a=point(m,key[0]),b=point(m,key[1]);const bool h=a.y==b.y;
    for(int v:key){const auto owner=m.verts(m.vertIndexToLabel[v]).supCell;if(owner<0)continue;
        const auto& c=m.cells(owner);const auto lo=m.verts(c.cellToVert[0]).coord(),hi=m.verts(c.cellToVert[2]).coord();
        if(h && (a.y==lo.y || a.y==hi.y) && a.x>=lo.x && b.x<=hi.x && b.x-a.x<hi.x-lo.x)return true;
        if(!h && (a.x==lo.x || a.x==hi.x) && a.y>=lo.y && b.y<=hi.y && b.y-a.y<hi.y-lo.y)return true;}return false;}
using Inventory=std::function<double(int,Vector2D,Vector2D)>;
// Operator fixture with explicit primary M. Q4 geometry uses I06 mixed
// closure and analytic center samples; this does not implement I11.
HermiteGalsState stateFor(const QuadtreeMesh& mesh,const Profile& f,bool radial,HermiteGalsTimeLayer time,const Inventory& inventory,
    const std::vector<LevelSetInterfacePanels::SymmetryPlane>& symmetryPlanes={},
    HermiteGalsTransport transport=HermiteGalsTransport::ConservativeVolume)
{
    HermiteGalsPrimaryState p;p.measure=radial?HermiteGalsMeasure::AxisymmetricVolume:HermiteGalsMeasure::CartesianArea;p.transport=transport;
    Patch patch;patch.activeCells.assign(mesh.cellIndexNum(),1);patch.nodes.resize(mesh.vertIndexNum());patch.cellCenters.resize(mesh.cellIndexNum());
    for(int v=0;v<mesh.vertIndexNum();++v)if(mesh.verts(mesh.vertIndexToLabel[v]).type()==TVertType::non){const auto j=f(point(mesh,v));patch.nodes[v]=j;p.vertices.push_back({v,j.value,j.dx,j.dy});}
    std::set<Key> keys;
    for(int c=0;c<mesh.cellIndexNum();++c){const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);const auto a=mesh.verts(cell.cellToVert[0]).coord(),b=mesh.verts(cell.cellToVert[2]).coord();
        p.geometryCells.push_back(c);if(transport==HermiteGalsTransport::ConservativeVolume)p.phaseVolumes.push_back({c,inventory(c,a,b)});patch.cellCenters[c]=f(a+(b-a)*.5).value;
        for(int k=0;k<4;++k){Key key{{mesh.verts(cell.cellToVert[k]).index(),mesh.verts(cell.cellToVert[(k+1)%4]).index()}};
            const auto x=point(mesh,key[0]),y=point(mesh,key[1]);if(x.y==y.y?x.x>y.x:x.y>y.y)std::swap(key[0],key[1]);keys.insert(key);}}
    for(const auto& key:keys)if(!fineSide(mesh,key)){const auto a=point(mesh,key[0]),b=point(mesh,key[1]);const auto j=f(a+(b-a)*.5);const auto dn=a.y==b.y?j.dy:j.dx;
        p.edges.push_back({key,j.value,dn});patch.edgeMidpoints.push_back({key,j.value,dn});}
    HermiteGalsCandidate candidate(mesh,time,p);const auto mixed=closeHermiteGalsMixedAmr(mesh,candidate,symmetryPlanes);
    HermiteGalsClosure closure;closure.vertexMixedDerivatives=mixed.vertexMixedDerivatives();closure.cellBubbleCoefficients=patch.cellCenters;
    for(std::size_t k=0;k<p.vertices.size();++k)patch.nodes[p.vertices[k].vertex].dxy=closure.vertexMixedDerivatives[k];
    candidate.bindGeometry(mesh,candidate.identity(),closure,LevelSetHermiteGeometry::Snapshot::fromNativeLeaves(mesh,patch,radial));return std::move(candidate).finish();
}
using Snapshot=LevelSetHermiteGeometry::Snapshot;
using Interval=Snapshot::QuarticCenterInterval;
static_assert(!std::is_default_constructible_v<HermiteGalsBubbleResult>);
static_assert(std::is_same_v<decltype(std::declval<HermiteGalsBubbleResult&>().closure()),const HermiteGalsClosure&>);
Patch patchFor(const QuadtreeMesh& mesh,const HermiteGalsCandidate& candidate,const HermiteGalsAmrClosureResult& mixed,const std::vector<double>& c)
{
    Patch patch;const auto& p=candidate.primary();patch.activeCells.assign(mesh.cellIndexNum(),0);
    patch.nodes.resize(mesh.vertIndexNum());patch.cellCenters.resize(mesh.cellIndexNum());
    for(std::size_t k=0;k<p.vertices.size();++k){const auto& v=p.vertices[k];patch.nodes[v.vertex]={v.phi,v.dx,v.dy,mixed.vertexMixedDerivatives()[k]};}
    for(const auto& e:p.edges)patch.edgeMidpoints.push_back({e.vertices,e.phi,e.dn});
    for(std::size_t k=0;k<p.geometryCells.size();++k){patch.activeCells[p.geometryCells[k]]=1;patch.cellCenters[p.geometryCells[k]]=c[k];}
    return patch;
}
std::shared_ptr<const Snapshot> scalarFor(const QuadtreeMesh& mesh,const Patch& p,bool radial)
{return std::make_shared<const Snapshot>(Snapshot::fromNativeLeaves(mesh,p,radial));}
std::vector<double> centers(const QuadtreeMesh& mesh,const Profile& f)
{std::vector<double> c;for(int k=0;k<mesh.cellIndexNum();++k)c.push_back(f(mesh.cells(mesh.cellIndexToLabel[k]).centerCoord()).value);return c;}
// Test-only second-order analytic differentiation for the exact continuum
// reference, independent of Snapshot/closure/tracer and its derivatives.
struct D2 {double v=0,x=0,y=0,xx=0,xy=0,yy=0;D2()=default;D2(double a):v(a){};D2(double a,double b,double c,double d,double e,double f):v(a),x(b),y(c),xx(d),xy(e),yy(f){}};
D2 operator+(D2 a,D2 b){return {a.v+b.v,a.x+b.x,a.y+b.y,a.xx+b.xx,a.xy+b.xy,a.yy+b.yy};}
D2 operator-(D2 a){return {-a.v,-a.x,-a.y,-a.xx,-a.xy,-a.yy};}
D2 operator-(D2 a,D2 b){return a+-b;}
D2 operator*(D2 a,D2 b){return {a.v*b.v,a.x*b.v+a.v*b.x,a.y*b.v+a.v*b.y,
    a.xx*b.v+2*a.x*b.x+a.v*b.xx,a.xy*b.v+a.x*b.y+a.y*b.x+a.v*b.xy,a.yy*b.v+2*a.y*b.y+a.v*b.yy};}
D2 compose(D2 a,double f,double fp,double fpp){return {f,fp*a.x,fp*a.y,fp*a.xx+fpp*a.x*a.x,fp*a.xy+fpp*a.x*a.y,fp*a.yy+fpp*a.y*a.y};}
D2 sine(D2 a){return compose(a,std::sin(a.v),std::cos(a.v),-std::sin(a.v));}
D2 cosine(D2 a){return compose(a,std::cos(a.v),-std::sin(a.v),-std::cos(a.v));}
D2 exponential(D2 a){const auto e=std::exp(a.v);return compose(a,e,e,e);}
struct Case
{
    std::string name;bool radial=false,amr=false;int motion=0;bool nonpolynomial=false;
    double speed=.2,period=.2,radius=.19,cx=.43,cy=.51;
    bool planar=false,axialTranslation=false;int planeSign=1;double planeOffset=.5;
    bool boundedInventory=false; // explicit benchmark contract, no failure fallback
    bool meshScaledGeometryBudget=false; // explicit fixed verification policy
    bool polynomialIntervalInventory=false; // explicit owned-source interval driver
    bool pureFirstJet=false; // explicitly selected nonconservative method
    bool axial()const{return radial || axialTranslation;}
    std::vector<LevelSetInterfacePanels::SymmetryPlane> symmetry(const QuadtreeMesh& mesh)const{
        if(!planar)return {};
        const int transverse=axial()?0:1;std::vector<LevelSetInterfacePanels::SymmetryPlane> result;
        if(!radial)result.push_back({transverse,transverse==0?mesh.minCoord.x:mesh.minCoord.y});
        result.push_back({transverse,transverse==0?mesh.maxCoord.x:mesh.maxCoord.y});return result;
    }
    D2 exact(Vector2D p,double t) const {
        D2 x{p.x,1,0,0,0,0},y{p.y,0,1,0,0,0};
        if(planar)return double(planeSign)*((axial()?y:x)-(planeOffset+(motion==1?speed*t:0)));
        if(motion==1){if(axial())y=y-speed*t;else x=x-speed*t;}
        if(motion>=2){const auto a=x-.5,b=y-.5;
            const D2 angle=motion==2?D2(-speed*t):(-speed*period/double(pi)*std::sin(double(pi)*t/period))*(a*a+b*b);
            x=.5+cosine(angle)*a-sine(angle)*b;y=.5+sine(angle)*a+cosine(angle)*b;}
        const auto a=x-(radial?0:cx),b=y-cy,q=a*a+b*b-radius*radius;
        return nonpolynomial?exponential(q)-1:q;
    }
    LevelSetBoundaryTrace::VelocityJet velocity(Vector2D p,double t) const {
        if(motion==0)return {};
        if(motion==1)return {{axial()?0:speed,axial()?speed:0},{0,0,0,0}};
        const auto x=p.x-.5,y=p.y-.5;
        if(motion==2)return {{-speed*y,speed*x},{0,-speed,speed,0}};
        const double k=speed*std::cos(double(pi)*t/period),r2=x*x+y*y;
        return {{-k*r2*y,k*r2*x},{-2*k*x*y,-k*(x*x+3*y*y),k*(3*x*x+y*y),2*k*x*y}};
    }
    LAMP::Geometry::PolynomialVelocityData polynomialData(double previous,double arrival) const {
        check(!radial && (motion==2 || motion==3),"polynomial driver requires original XY rotation/vortex");
        LAMP::Geometry::PolynomialVelocityData data;data.origin={.5,.5};
        const auto coefficients=[&](double time){std::array<Vector2D,16> c{};
            const double k=motion==2?speed:speed*std::cos(double(pi)*time/period);
            if(motion==2){c[1].x=-k;c[4].y=k;}
            else {c[3].x=-k;c[9].x=-k;c[12].y=k;c[6].y=k;}return c;};
        data.previous=coefficients(previous);data.arrival=coefficients(arrival);return data;
    }
    Profile profile(double t) const {return [s=*this,t](Vector2D p){const auto a=s.exact(p,t);return Jet{a.v,a.x,a.y,a.xy};};}
};
LevelSetBoundaryTrace::Domain traceDomain(const QuadtreeMesh& mesh,bool radial)
{
    Vector2D h{1,1};for(int c=0;c<mesh.cellIndexNum();++c){const auto size=mesh.cells(mesh.cellIndexToLabel[c]).size();h.x=std::min(h.x,size.x);h.y=std::min(h.y,size.y);}
    return {mesh.minCoord,mesh.maxCoord,h,radial};
}
struct GeometryBudget {double minimumWidth=0,domainLength=0,tolerance=1e-10;};
GeometryBudget geometryBudget(const QuadtreeMesh& mesh,const Case& s)
{
    GeometryBudget budget;if(!s.meshScaledGeometryBudget)return budget;
    check(s.boundedInventory,"mesh geometry budget requires bounded inventory transport");
    budget.domainLength=std::max(mesh.maxCoord.x-mesh.minCoord.x,mesh.maxCoord.y-mesh.minCoord.y);
    budget.minimumWidth=std::numeric_limits<double>::infinity();
    for(int c=0;c<mesh.cellIndexNum();++c){const auto size=mesh.cells(mesh.cellIndexToLabel[c]).size();
        budget.minimumWidth=std::min({budget.minimumWidth,size.x,size.y});}
    check(std::isfinite(budget.domainLength) && budget.domainLength>0 && std::isfinite(budget.minimumWidth)
        && budget.minimumWidth>0,"invalid mesh geometry budget scale");
    const double ratio=budget.minimumWidth/budget.domainLength;
    budget.tolerance=std::max(1e-10,ratio*ratio*ratio);
    check(std::isfinite(budget.tolerance) && budget.tolerance<1,"invalid mesh geometry tolerance");
    return budget;
}
HermiteGalsBubbleConfig bubbleControls(const QuadtreeMesh& mesh,const Case& s,const HermiteGalsPrimaryState& p)
{
    HermiteGalsBubbleConfig config;config.symmetryPlanes=s.symmetry(mesh);
    config.relativeVolumeTolerance=geometryBudget(mesh,s).tolerance;
    if(s.boundedInventory){config.volumeConstraint=HermiteGalsBubbleVolumeConstraint::GeometryAccuracy;
        for(int c:p.geometryCells){const auto size=mesh.cells(mesh.cellIndexToLabel[c]).size();const double h=std::min(size.x,size.y);
            config.accuracyBudgets.push_back({c,h*h});}}
    return config;
}
void recordBoundedClosure(const Case& s,HermiteGalsTimeLayer time,const HermiteGalsBubbleConfig& config,const HermiteGalsBubbleResult& bubble)
{
    if(!s.boundedInventory)return;double maximumBudget=0;
    for(const auto& c:config.accuracyBudgets)maximumBudget=std::max(maximumBudget,c.maximumHessianComponentChange);
    const auto& r=bubble.report();std::cout<<"BOUNDED_CLOSURE,"<<s.name<<','<<s.radial<<','<<s.amr<<','<<time.step<<','<<time.time()<<','
        <<(time.step==0?"initial":time.position==Position::Midpoint?"midpoint":"end")<<','<<config.relativeVolumeTolerance<<','<<maximumBudget<<','
        <<r.initialMaxRelativeResidual<<','<<r.finalMaxRelativeResidual<<','<<r.maxPureRelativeResidual<<','<<r.maxHessianComponentChangeBound<<','
        <<r.correctedCells<<','<<r.signedResidual<<','<<r.l1Residual<<std::endl;
}
HermiteGalsPublication initial(const QuadtreeMesh& mesh,const Case& s)
{
    if(s.pureFirstJet){
        check(!s.boundedInventory && !s.meshScaledGeometryBudget && !s.polynomialIntervalInventory,"FirstJet cannot enable conserved inventory route");
        const auto f=s.profile(0);const auto seed=stateFor(mesh,f,s.radial,{},[](int,Vector2D,Vector2D){return 0.;},
            s.symmetry(mesh),HermiteGalsTransport::FirstJet);
        HermiteGalsCandidate candidate(mesh,{},seed.primary());const auto mixed=closeHermiteGalsMixedAmr(mesh,candidate,s.symmetry(mesh));
        HermiteGalsPublicationConfig config;config.volume.symmetryPlanes=s.symmetry(mesh);
        return publishHermiteGalsFirstJetGeometry(mesh,std::move(candidate),mixed,{candidate.identity(),centers(mesh,f)},config);
    }
    if(s.meshScaledGeometryBudget){const auto b=geometryBudget(mesh,s);
        std::cout<<"GEOMETRY_BUDGET,"<<s.name<<','<<s.radial<<','<<s.amr<<','<<b.minimumWidth<<','<<b.domainLength<<','<<b.tolerance<<std::endl;}
    const auto f=s.profile(0);const auto seed=stateFor(mesh,f,s.radial,{},[](int,Vector2D,Vector2D){return 0.;},s.symmetry(mesh));
    HermiteGalsPhaseVolumeConfig volumeConfig;volumeConfig.symmetryPlanes=s.symmetry(mesh);
    const auto measure=integrateHermiteGalsPhaseVolumes(mesh,seed,volumeConfig);auto p=seed.primary();p.phaseVolumes=measure.phaseVolumes();
    HermiteGalsCandidate candidate(mesh,{},p);const auto mixed=closeHermiteGalsMixedAmr(mesh,candidate,s.symmetry(mesh));
    const auto bubbleConfig=bubbleControls(mesh,s,p);
    const auto bubble=closeHermiteGalsBubble(mesh,candidate,mixed,{candidate.identity(),centers(mesh,f)},bubbleConfig);
    recordBoundedClosure(s,candidate.timeLayer(),bubbleConfig,bubble);
    HermiteGalsPublicationConfig publicationConfig;publicationConfig.volume=volumeConfig;
    publicationConfig.relativeInventoryTolerance=bubbleConfig.relativeVolumeTolerance;
    return publishHermiteGalsGeometry(mesh,std::move(candidate),bubble,publicationConfig);
}
long double circleCell(Vector2D a,Vector2D b,Vector2D center,double radius,bool radial)
{
    const long double cx=center.x,cy=center.y,r=radius,r2=r*r;
    const long double lo=std::max((long double)a.x,cx-r),hi=std::min((long double)b.x,cx+r);if(hi<=lo)return 0;
    std::vector<long double> cuts{lo,hi};
    for(long double y:{(long double)a.y,(long double)b.y}){const auto d=y-cy;if(std::abs(d)<r){const auto s=std::sqrt(r2-d*d);
        for(const auto x:{cx-s,cx+s})if(x>lo && x<hi)cuts.push_back(x);}}
    std::sort(cuts.begin(),cuts.end());long double result=0;
    // Near |d|=r, asin(d/r) amplifies the independent ratio rounding; its
    // leading endpoint term then fails to cancel d*s. atan2 uses this SAME
    // chord height, preserving the analytic primitive in binary64.
    const auto arc=[&](long double x){const auto d=x-cx,s=std::sqrt(std::max(0.L,r2-d*d));return (d*s+r2*std::atan2(d,s))/2;};
    for(std::size_t k=1;k<cuts.size();++k){const auto l=cuts[k-1],h=cuts[k],mid=(l+h)/2,s=std::sqrt(std::max(0.L,r2-(mid-cx)*(mid-cx)));
        if(std::min((long double)b.y,cy+s)<=std::max((long double)a.y,cy-s))continue;
        const bool topArc=cy+s<b.y,bottomArc=cy-s>a.y;const auto constant=(topArc?cy:b.y)-(bottomArc?cy:a.y);const int factor=int(topArc)+int(bottomArc);
        if(radial){check(center.x==0,"RZ reference requires on-axis sphere");result+=2*pi*(constant*(h*h-l*l)/2+factor*(std::pow(std::max(0.L,r2-l*l),1.5L)-std::pow(std::max(0.L,r2-h*h),1.5L))/3);}
        else result+=constant*(h-l)+factor*(arc(h)-arc(l));}
    return result;
}
struct Errors {double position=0,gradient=0,curvature=0,volume=0,drift=0,cellVolume=0;};
Errors errors(const HermiteGalsPublication& pub,const Case& s)
{
    Errors e;double sum=0;const auto time=pub.state().timeLayer().time();
    const auto& phaseVolumes=s.pureFirstJet?pub.volumes().phaseVolumes():pub.state().primary().phaseVolumes;
    for(std::size_t c=0;c<pub.topology().curves.size();++c)for(std::size_t g=0;g<pub.topology().curves[c].graphs().size();++g)for(const auto& q:pub.quadrature(c,g)){
        const auto& f=q.sample.frame;const auto j=s.exact(f.point,time);const auto n=std::hypot(j.x,j.y);
        const auto actual=pub.state().geometry().snapshot().scalarJetAtNativeCellChart(f.nativeCellIndex,f.cellChartPoint);
        e.position=std::max(e.position,std::abs(j.v)/n);e.gradient=std::max(e.gradient,std::hypot(actual.dx-j.x,actual.dy-j.y)/n);
        double curvature=(j.y*j.y*j.xx-2*j.x*j.y*j.xy+j.x*j.x*j.yy)/(n*n*n);if(s.radial)curvature+=j.x/(n*f.point.x);
        e.curvature=std::max(e.curvature,std::abs(q.sample.curvature-curvature));}
    for(const auto& m:phaseVolumes)sum+=m.liquid;
    const double expected=s.radial?4*double(pi)*s.radius*s.radius*s.radius/3:double(pi)*s.radius*s.radius;e.volume=std::abs(sum-expected);
    if(s.motion<3){Vector2D center{s.radial?0:s.cx,s.cy};if(s.motion==1){if(s.radial)center.y+=s.speed*time;else center.x+=s.speed*time;}
        if(s.motion==2){const auto d=center-Vector2D(.5,.5);const auto angle=s.speed*time;center={.5+std::cos(angle)*d.x-std::sin(angle)*d.y,.5+std::sin(angle)*d.x+std::cos(angle)*d.y};}
        long double exactSum=0;for(const auto& m:phaseVolumes){const auto bounds=pub.state().geometry().snapshot().nativeCellBounds(m.cell);
            const auto exact=circleCell(bounds[0],bounds[1],center,s.radius,s.radial);exactSum+=exact;
            e.cellVolume=std::max(e.cellVolume,double(std::abs(m.liquid-exact)/pub.domain().cellVolumes()[m.cell]));}
        close(exactSum,expected,1e-12,"independent analytic cell volumes do not sum to circle/sphere");
    }else e.cellVolume=std::numeric_limits<double>::quiet_NaN(); // no per-cell vortex reference claimed
    return e;
}
// Explicit I12r4 diagnostic. Analytic inventories never enter advance().
double curvatureOf(const Snapshot::ScalarJet2& j)
{const double n=std::hypot(j.dx,j.dy);check(n>0,"audit zero gradient");
    return (j.dy*j.dy*j.dxx-2*j.dx*j.dy*j.dxy+j.dx*j.dx*j.dyy)/(n*n*n);}
// Same-point diagnostic only. The mixed-pullback patch is never published
HermiteGalsPublication advanceFirstJet(const QuadtreeMesh& mesh,const HermiteGalsPublication& source,
    const Case& s,HermiteGalsTimeLayer time,std::string& stage,
    const HermiteGalsTemporalVelocitySource* sharedVelocity=nullptr,
    const HermiteGalsTemporalCharacteristicSource* sharedCharacteristic=nullptr)
{
    check(s.pureFirstJet && bool(sharedVelocity)==bool(sharedCharacteristic),"explicit pure velocity inputs required");
    const auto& old=source.state();HermiteGalsCandidate velocityGeometry(mesh,time,old.primary());
    const auto velocity=sharedVelocity?*sharedVelocity:s.motion>=2
        ?HermiteGalsTemporalVelocitySource::polynomial(old,velocityGeometry.identity(),time,s.polynomialData(old.timeLayer().time(),time.time()))
        :HermiteGalsTemporalVelocitySource::constant(old,velocityGeometry.identity(),time,s.velocity({},old.timeLayer().time()).value);
    const auto characteristic=sharedCharacteristic?*sharedCharacteristic:HermiteGalsTemporalCharacteristicSource(mesh,traceDomain(mesh,s.radial));
    characteristic.validateMesh(mesh);
    const auto arrival=sharedVelocity?(velocity.polynomialVelocity()?velocity.polynomialArrivalView(velocityGeometry.identity(),time)
        :velocity.constantArrivalView(velocityGeometry.identity(),time)):velocity.arrival();
    const HermiteGalsBoundaryFirstJetQuery boundary=[&](Vector2D q,HermiteGalsTimeLayer t,Boundary){
        const auto j=s.exact(q,t.time());return HermiteGalsBoundaryFirstJet{j.v,{j.x,j.y}};};
    HermiteGalsPublicationConfig config;config.volume.symmetryPlanes=s.symmetry(mesh);
    stage="first_jet_step_producer";
    return prepareHermiteGalsFirstJetStep(mesh,source,characteristic.domain(),velocity.previous(),arrival,boundary,config).publication();
}
int failures=0;
} // namespace
