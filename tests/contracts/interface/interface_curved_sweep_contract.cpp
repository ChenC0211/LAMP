#include "Modules/Interface/HermiteGalsTemporalPhaseFlux.h"
#include "Modules/Interface/HermiteGalsPhaseRegion.h"
#include "Modules/Interface/HermiteGalsAmrClosure.h"
#include "Core/Geometry/ActiveQuarticHermitePatch.h"
#include "Core/Geometry/LevelSetHermiteGeometry.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <algorithm>
#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
namespace {
using namespace LAMP::Interface;
using Patch=LAMP::Geometry::ActiveQuarticHermitePatch;
using Snapshot=LevelSetHermiteGeometry::Snapshot;
using Jet=Patch::NodalDof;
using Profile=std::function<Jet(Vector2D)>;
using Key=std::array<int,2>;
using Velocity=HermiteGalsTemporalVelocitySource;
using Characteristic=HermiteGalsTemporalCharacteristicSource;
using Interval=HermiteGalsPhaseFluxIntervalSource;
using Position=HermiteGalsTimeLayer::Position;
constexpr long double pi=3.1415926535897932384626433832795L;
int cases=0,rejections=0,decimalBands=0;
void check(bool ok,const std::string& why){if(!ok)throw std::runtime_error(why);}
template<class F>void rejects(F f,const char* reason){try{f();}catch(const std::exception& e){
    check(std::string(e.what()).find(reason)!=std::string::npos,std::string("unexpected rejection: ")+e.what());++rejections;return;}
    throw std::runtime_error(std::string("missing rejection: ")+reason);}
void near(long double a,long double b,long double scale,const char* why){if(std::abs(a-b)>1e-9L*scale){
    std::cerr<<std::setprecision(18)<<why<<" actual="<<a<<" expected="<<b<<" scale="<<scale<<'\n';throw std::runtime_error(why);}}
std::unique_ptr<QuadtreeMesh> meshFor(int n,bool amr=false,double radialMin=0)
{
    DictionaryMesh d;d.Nx=d.Ny=n;d.minCoord={radialMin,0};d.maxCoord={1,1};d.maxCellNum=32768;d.maxInitialIterations=0;
    auto m=std::make_unique<QuadtreeMesh>(d);
    if(amr){DictionaryMesh::RefineBox b;b.level=1;b.min={radialMin,0.25};b.max={0.75,0.75};m->applyRefineBoxes({b});}
    return m;
}
Vector2D point(const QuadtreeMesh& m,int v){return m.verts(m.vertIndexToLabel[v]).coord();}
std::set<Key> sides(const QuadtreeMesh& m)
{
    std::set<Key> keys;
    for(int c=0;c<m.cellIndexNum();++c){const auto& cell=m.cells(m.cellIndexToLabel[c]);
        for(int k=0;k<4;++k){Key key{{m.verts(cell.cellToVert[k]).index(),m.verts(cell.cellToVert[(k+1)%4]).index()}};
            auto a=point(m,key[0]),b=point(m,key[1]);if(a.y==b.y?a.x>b.x:a.y>b.y)std::swap(key[0],key[1]);keys.insert(key);}}
    return keys;
}
bool fineSide(const QuadtreeMesh& m,Key key)
{
    const auto a=point(m,key[0]),b=point(m,key[1]);const bool h=a.y==b.y;
    for(int v:key){const int owner=m.verts(m.vertIndexToLabel[v]).supCell;if(owner<0)continue;
        const auto& c=m.cells(owner);const auto lo=m.verts(c.cellToVert[0]).coord(),hi=m.verts(c.cellToVert[2]).coord();
        if(h && (a.y==lo.y || a.y==hi.y) && a.x>=lo.x && b.x<=hi.x && b.x-a.x<hi.x-lo.x)return true;
        if(!h && (a.x==lo.x || a.x==hi.x) && a.y>=lo.y && b.y<=hi.y && b.y-a.y<hi.y-lo.y)return true;}
    return false;
}
// Initial GEOMETRY TEMPLATE only, never an accepted solution. Center samples
// supply the analytic Q4 interpolant. I06 supplies the actual compact mixed
// closure; I07 initializes the zero placeholder inventory. I11 remains absent.
HermiteGalsState initialTemplate(const QuadtreeMesh& mesh,const Profile& f,bool radial,
    HermiteGalsTimeLayer time={},bool mismatchedMeasure=false,bool extraInventory=false)
{
    HermiteGalsPrimaryState p;p.measure=radial?HermiteGalsMeasure::AxisymmetricVolume:HermiteGalsMeasure::CartesianArea;
    Patch patch;patch.nodes.resize(mesh.vertIndexNum());patch.activeCells.assign(mesh.cellIndexNum(),1);
    patch.cellCenters.resize(mesh.cellIndexNum());
    for(int v=0;v<mesh.vertIndexNum();++v)if(mesh.verts(mesh.vertIndexToLabel[v]).type()==TVertType::non){
        const auto j=f(point(mesh,v));p.vertices.push_back({v,j.value,j.dx,j.dy});patch.nodes[v]=j;}
    for(const auto& key:sides(mesh))if(!fineSide(mesh,key)){
        const auto a=point(mesh,key[0]),b=point(mesh,key[1]);const auto j=f(a+(b-a)*.5);const auto dn=a.y==b.y?j.dy:j.dx;
        p.edges.push_back({key,j.value,dn});patch.edgeMidpoints.push_back({key,j.value,dn});}
    for(int c=0;c<mesh.cellIndexNum();++c){p.geometryCells.push_back(c);p.phaseVolumes.push_back({c,0});
        const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);const auto a=mesh.verts(cell.cellToVert[0]).coord(),b=mesh.verts(cell.cellToVert[2]).coord();
        patch.cellCenters[c]=f(a+(b-a)*.5).value;}
    HermiteGalsCandidate seed(mesh,time,p);const auto mixed=closeHermiteGalsMixedAmr(mesh,seed);
    HermiteGalsClosure closure;closure.vertexMixedDerivatives=mixed.vertexMixedDerivatives();closure.cellBubbleCoefficients=patch.cellCenters;
    for(std::size_t i=0;i<p.vertices.size();++i)patch.nodes[p.vertices[i].vertex].dxy=closure.vertexMixedDerivatives[i];
    if(extraInventory){p.geometryCells.pop_back();patch.activeCells.back()=0;closure.cellBubbleCoefficients.pop_back();
        std::set<Key> activeSides;
        for(int c:p.geometryCells){const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);
            for(int k=0;k<4;++k){Key key{{mesh.verts(cell.cellToVert[k]).index(),mesh.verts(cell.cellToVert[(k+1)%4]).index()}};
                const auto a=point(mesh,key[0]),b=point(mesh,key[1]);if(a.y==b.y?a.x>b.x:a.y>b.y)std::swap(key[0],key[1]);activeSides.insert(key);}}
        patch.edgeMidpoints.erase(std::remove_if(patch.edgeMidpoints.begin(),patch.edgeMidpoints.end(),
            [&](const auto& edge){return activeSides.count(edge.vertices)==0;}),patch.edgeMidpoints.end());}
    if(mismatchedMeasure)p.measure=radial?HermiteGalsMeasure::CartesianArea:HermiteGalsMeasure::AxisymmetricVolume;
    HermiteGalsCandidate candidate(mesh,time,p);
    candidate.bindGeometry(mesh,candidate.identity(),closure,Snapshot::fromNativeLeaves(mesh,patch,radial));
    return std::move(candidate).finish();
}
struct Circle{double x,y,r;};
Jet circleJet(Vector2D p,Circle c)
{const auto x=p.x-c.x,y=p.y-c.y;return {x*x+y*y-c.r*c.r,2*x,2*y,0};}
Jet product(Jet a,Jet b)
{return {a.value*b.value,a.dx*b.value+a.value*b.dx,a.dy*b.value+a.value*b.dy,
    a.dxy*b.value+a.dx*b.dy+a.dy*b.dx+a.value*b.dxy};}
long double diskPrimitive(long double u,long double r)
{u=std::clamp(u,-r,r);const auto s=std::sqrt(std::max(0.L,r*r-u*u));return .5L*(u*s+r*r*std::atan2(u,s));}
// Independent analytic rectangle/disk intersection: split only at exact
// circle crossings of the rectangle's horizontal bounds, then integrate
// constant + k*sqrt(R^2-(x-cx)^2). No production graph, roots or quadrature.
long double diskVolume(Vector2D lo,Vector2D hi,Circle c,bool radial)
{
    const long double cx=c.x,cy=c.y,r=c.r;
    const long double a=std::max<long double>(lo.x,cx-r),b=std::min<long double>(hi.x,cx+r);
    if(b<=a)return 0;
    std::vector<long double> points{a,b};
    for(long double y:{static_cast<long double>(lo.y),static_cast<long double>(hi.y)})
        if(std::abs(y-cy)<r){const auto u=std::sqrt(r*r-(y-cy)*(y-cy));for(auto x:{cx-u,cx+u})if(x>a && x<b)points.push_back(x);}
    std::sort(points.begin(),points.end());long double sum=0;
    const auto primitive=[&](long double x){const auto u=x-cx;
        return radial?cx*diskPrimitive(u,r)-std::pow(std::max(0.L,r*r-u*u),1.5L)/3:diskPrimitive(u,r);};
    for(std::size_t i=1;i<points.size();++i){const auto x0=points[i-1],x1=points[i],xm=(x0+x1)/2;
        const auto s=std::sqrt(std::max(0.L,r*r-(xm-cx)*(xm-cx)));
        if(std::min<long double>(hi.y,cy+s)<=std::max<long double>(lo.y,cy-s))continue;
        long double constant=0,k=0;
        if(cy+s<hi.y){constant+=cy;k+=1;}else constant+=hi.y;
        if(cy-s>lo.y){constant-=cy;k+=1;}else constant-=lo.y;
        sum+=constant*(radial?(x1*x1-x0*x0)/2:x1-x0)+k*(primitive(x1)-primitive(x0));}
    return radial?2*pi*sum:sum;
}

LevelSetBoundaryTrace::Domain traceDomain(const QuadtreeMesh& mesh,bool radial){
    Vector2D h{1,1};for(int c=0;c<mesh.cellIndexNum();++c){const auto d=mesh.cells(mesh.cellIndexToLabel[c]).size();h.x=std::min(h.x,d.x);h.y=std::min(h.y,d.y);}
    return {mesh.minCoord,mesh.maxCoord,h,radial};
}
int faceFor(const QuadtreeMesh& mesh,const Interval& interval,bool horizontal,double normal,double tangent){
    int best=-1;double distance=2;
    for(int row=0;row<mesh.vEdgeIndexNum()+mesh.hEdgeIndexNum();++row){const auto f=interval.face(mesh,row);
        if(f.horizontal!=horizontal || (horizontal?f.center.y:f.center.x)!=normal)continue;
        const auto d=std::abs((horizontal?f.center.x:f.center.y)-tangent);
        if(d<=distance){distance=d;best=row;}}
    check(best>=0,"missing native face");return best;
}
std::array<Vector2D,2> ends(const QuadtreeMesh& mesh,const NativeLeafFace& f){
    const auto& edge=f.horizontal?mesh.hEdges(f.edgeLabel):mesh.vEdges(f.edgeLabel);
    auto a=mesh.verts(edge.edgeToVert[0]).coord(),b=mesh.verts(edge.edgeToVert[1]).coord();
    if(f.horizontal?a.x>b.x:a.y>b.y)std::swap(a,b);return {{a,b}};
}
void matrix(){
    for(bool radial:{false,true})for(bool horizontal:{false,true}){
        if(radial && !horizontal)continue;
        for(bool amr:{false,true}){auto mesh=meshFor(4,amr);const Circle c{radial?0:.4375,.484375,.203125};
            for(int phase:{-1,1}){const auto state=initialTemplate(*mesh,[&](Vector2D p){auto j=circleJet(p,c);j.value*=phase;j.dx*=phase;j.dy*=phase;return j;},radial);
                const HermiteGalsPhaseRegionSource regions(*mesh,state);const Characteristic map(*mesh,traceDomain(*mesh,radial));
                for(int direction:{-1,1})for(double dt:{1./64,.02}){
                    const double speed=.25*direction;const auto v=Velocity::constant(state,state.identity(),{1,0,dt,Position::End},horizontal?Vector2D{0,speed}:Vector2D{speed,0});
                    const Interval full(*mesh,state,v,map,0,dt),first(*mesh,state,v,map,0,dt*.375),last(*mesh,state,v,map,dt*.375,dt);
                    const int row=faceFor(*mesh,full,horizontal,.5,radial?.1875:horizontal?.625:.375);
                    const auto face=full.face(*mesh,row);const auto ab=ends(*mesh,face);
                    auto lo=ab[0],hi=ab[1];if(horizontal){lo.y=std::min(ab[0].y,ab[0].y-speed*dt);hi.y=std::max(ab[0].y,ab[0].y-speed*dt);}
                    else{lo.x=std::min(ab[0].x,ab[0].x-speed*dt);hi.x=std::max(ab[0].x,ab[0].x-speed*dt);}
                    const long double area=horizontal?(static_cast<long double>(ab[1].x)-ab[0].x)*(radial?pi*(static_cast<long double>(ab[1].x)+ab[0].x):1):ab[1].y-ab[0].y;
                    const long double total=area*std::abs(speed)*dt,inside=diskVolume(lo,hi,c,radial);
                    const auto measured=integrateHermiteGalsSweptFace(*mesh,full,row,regions);
                    near(measured.transfer.signedTotalVolume(),direction*total,total,"swept total");
                    near(measured.transfer.signedLiquidVolume(),direction*(phase>0?inside:total-inside),total,"swept analytic phase");
                    const auto& gross=measured.transfer.integratedVolumes();
                    check(gross.integratedLiquidPositiveVolume>=0 && gross.integratedLiquidPositiveVolume<=gross.integratedTotalPositiveVolume &&
                        gross.integratedLiquidNegativeVolume>=0 && gross.integratedLiquidNegativeVolume<=gross.integratedTotalNegativeVolume,"gross bounds");
                    const auto joined=integrateHermiteGalsSweptFace(*mesh,first,row,regions).transfer.appended(*mesh,
                        integrateHermiteGalsSweptFace(*mesh,last,row,regions).transfer);
                    near(joined.signedLiquidVolume(),measured.transfer.signedLiquidVolume(),total,"swept interval composition");
                    near(joined.signedTotalVolume(),measured.transfer.signedTotalVolume(),total,"swept total composition");
                    check(joined.interval().sameInterval(full) && measured.transfer.interval().sameInterval(full),"interval provenance");
                    check(measured.report.rectangleQueries>=1 && measured.report.rectangleQueries<=2 && measured.report.coordinateVolumeBand>=0 &&
                        measured.report.coordinateVolumeBand<=1e-10L*total && measured.report.momentErrorEstimate<=1e-10L*total,"separate sweep error diagnostics");
                    if(dt==.02 && measured.report.coordinateVolumeBand>0)++decimalBands;
                    for(const auto& m:state.primary().phaseVolumes)check(m.liquid==0,"sweep changed independent M");
                    check(&regions.source().geometry().snapshot()==&state.geometry().snapshot(),"sweep replaced old scalar");
                    ++cases;
                }
            }
        }
    }
    check(decimalBands>0,"decimal intervals did not exercise coordinate envelopes");
}
void missedEvents(){
    auto mesh=meshFor(4);const double dt=1./32,r=dt/32;
    for(bool radial:{false,true})for(bool horizontal:{false,true}){
        if(radial && !horizontal)continue;
        for(int direction:{-1,1}){
            const double normalCenter=.5-direction*3*dt/32;
            const Circle c{radial?0:horizontal?.375:normalCenter,horizontal?normalCenter:.375,r};
            // XY: a thin ellipse crosses existing tangent-direction leaf
            // boundaries, so complete negative-fragment ownership is resolved.
            // A tiny circle wholly inside one old leaf still explicitly rejects
            // in Core; that distinct spatial limitation is tested below.
            const double tangentRadius=.25,k=(r/tangentRadius)*(r/tangentRadius);
            const auto state=initialTemplate(*mesh,[&](Vector2D p){if(radial)return circleJet(p,c);
                const double n=(horizontal?p.y:p.x)-normalCenter,t=(horizontal?p.x:p.y)-.375;
                return Jet{n*n+k*t*t-r*r,horizontal?2*k*t:2*n,horizontal?2*n:2*k*t,0};},radial);
            const HermiteGalsPhaseRegionSource regions(*mesh,state);const Characteristic map(*mesh,traceDomain(*mesh,radial));
            const auto v=Velocity::constant(state,state.identity(),{1,0,dt,Position::End},horizontal?Vector2D{0,double(direction)}:Vector2D{double(direction),0});
            const Interval interval(*mesh,state,v,map,0,dt);const int row=faceFor(*mesh,interval,horizontal,.5,radial?.125:.375);
            // The entire face is gas at ALL coarse and split midpoint nodes:
            // its normal distance alone exceeds the radius. This is not a
            // point-sampled claim about the rest of a face.
            for(double fraction:{.25,.5,.75})check(std::abs(.5-direction*dt*fraction-normalCenter)>r,"manufactured hidden event not hidden");
            const auto measured=integrateHermiteGalsSweptFace(*mesh,interval,row,regions);
            const auto ab=ends(*mesh,interval.face(*mesh,row));
            const double tangentLow=horizontal?ab[0].x:ab[0].y,tangentHigh=horizontal?ab[1].x:ab[1].y;
            const long double expected=radial?4*pi*r*r*r/3:(r/tangentRadius)*
                diskVolume({-tangentRadius,tangentLow},{tangentRadius,tangentHigh},{0,.375,tangentRadius},false);
            check(std::abs(measured.transfer.signedLiquidVolume())>0,"missed positive-measure time event");
            near(measured.transfer.signedLiquidVolume(),direction*expected,expected,"complete hidden component");
            std::cout<<std::setprecision(17)<<"HIDDEN_EVENT,"<<radial<<','<<horizontal<<','<<direction<<','<<measured.transfer.signedLiquidVolume()<<','<<expected<<",midpoint_liquid=0\n";
            ++cases;
        }
    }
}
void guards(){
    auto mesh=meshFor(4);const Circle c{.4375,.484375,.203125};
    const auto state=initialTemplate(*mesh,[&](Vector2D p){return circleJet(p,c);},false);
    const HermiteGalsPhaseRegionSource regions(*mesh,state);const auto domain=traceDomain(*mesh,false);const Characteristic map(*mesh,domain);
    const auto v=Velocity::constant(state,state.identity(),{1,0,.02,Position::End},{1,0});const Interval interval(*mesh,state,v,map,0,.02);
    const int row=faceFor(*mesh,interval,false,.5,.375);
    // Independent known decimal rounding: .5 - .25*.02 lies strictly
    // between .495 and its successor; the negative-speed endpoint lies
    // between the predecessor of .505 and .505. Tangent endpoints are exact.
    for(int sign:{-1,1}){const auto sweep=state.geometry().snapshot().nativeNormalTranslationSweep(
        {.5,.25},{.5,.5},sign*.25,0,0,.02,false);
        check(sweep.outerLower.x==(sign>0?.495:.5) && sweep.outerUpper.x==(sign>0?.5:.505) &&
            sweep.innerLower.x==(sign>0?std::nextafter(.495,1.):.5) &&
            sweep.innerUpper.x==(sign>0?.5:std::nextafter(.505,0.)),"independent decimal endpoint enclosures");
        check(sweep.innerLower.y==.25 && sweep.outerLower.y==.25 && sweep.innerUpper.y==.5 && sweep.outerUpper.y==.5,
            "sweep changed exact tangent bounds");
        check(sweep.coordinateVolumeBand==(sign>0?0x1p-56L:0x1p-55L),"independent coordinate volume band");
        near(sweep.totalVolume,.00125,.00125,"independent exact total volume");++cases;
    }
    const auto different=initialTemplate(*mesh,[&](Vector2D p){return circleJet(p,c);},false);const HermiteGalsPhaseRegionSource foreign(*mesh,different);
    rejects([&]{integrateHermiteGalsSweptFace(*mesh,interval,row,foreign);},"identical old region source");
    auto badMesh=meshFor(8);rejects([&]{integrateHermiteGalsSweptFace(*badMesh,interval,row,regions);},"mesh");
    const auto query=[](Vector2D){return LevelSetBoundaryTrace::VelocityJet{{1,0},{0,0,0,0}};};
    const Velocity opaque(state,state.identity(),{1,0,.02,Position::End},query,query);const Interval oi(*mesh,state,opaque,map,0,.02);
    rejects([&]{integrateHermiteGalsSweptFace(*mesh,oi,row,regions);},"owned constant velocity");
    auto unsupported=domain;unsupported.activeSegment=[](Vector2D,Vector2D){return true;};const Characteristic unknown(*mesh,unsupported);
    const Interval ui(*mesh,state,v,unknown,0,.02);rejects([&]{integrateHermiteGalsSweptFace(*mesh,ui,row,regions);},"opaque coverage");
    const auto tangent=Velocity::constant(state,state.identity(),{1,0,.02,Position::End},{1,.5});const Interval ti(*mesh,state,tangent,map,0,.02);
    rejects([&]{integrateHermiteGalsSweptFace(*mesh,ti,row,regions);},"zero tangential");
    rejects([&]{integrateHermiteGalsSweptFace(*mesh,interval,faceFor(*mesh,interval,false,0,.375),regions);},"exterior history required");
    const double tiny=0x1p-60;const auto small=Velocity::constant(state,state.identity(),{1,0,tiny,Position::End},{1,0});
    const Interval si(*mesh,state,small,map,0,tiny);rejects([&]{integrateHermiteGalsSweptFace(*mesh,si,row,regions);},"not distinguishable");
    HermiteGalsSweptFaceConfig invalid;invalid.relativeCoordinateTolerance=0;
    rejects([&]{integrateHermiteGalsSweptFace(*mesh,interval,row,regions,invalid);},"invalid swept integral tolerance");
    invalid.relativeCoordinateTolerance=1e-30;
    rejects([&]{integrateHermiteGalsSweptFace(*mesh,interval,row,regions,invalid);},"coordinate volume band");
    const auto stationary=Velocity::constant(state,state.identity(),{1,0,.02,Position::End},{0,0});const Interval rest(*mesh,state,stationary,map,0,.02);
    const auto zero=integrateHermiteGalsSweptFace(*mesh,rest,row,regions);
    check(zero.transfer.signedTotalVolume()==0 && zero.transfer.signedLiquidVolume()==0 && zero.report.rectangleQueries==0,"stationary volume");++cases;
    const auto cross=integrateHermiteGalsSweptFace(*mesh,interval,faceFor(*mesh,interval,true,.5,.375),regions);
    check(cross.transfer.signedLiquidVolume()==0 && cross.transfer.signedTotalVolume()==0,"zero normal/tangential-only volume");++cases;
    const Circle unresolved{.5-3./1024,.375,1./1024};
    const auto smallCircle=initialTemplate(*mesh,[&](Vector2D p){return circleJet(p,unresolved);},false);
    rejects([&]{HermiteGalsPhaseRegionSource unsupportedCycle(*mesh,smallCircle);},"closed subcell cycle");
    const Circle sphere{0,.484375,.203125};const auto rz=initialTemplate(*mesh,[&](Vector2D p){return circleJet(p,sphere);},true);
    const HermiteGalsPhaseRegionSource rzRegions(*mesh,rz);const Characteristic rzMap(*mesh,traceDomain(*mesh,true));
    const auto rv=Velocity::constant(rz,rz.identity(),{1,0,.02,Position::End},{1,0});const Interval ri(*mesh,rz,rv,rzMap,0,.02);
    rejects([&]{integrateHermiteGalsSweptFace(*mesh,ri,row,rzRegions);},"axial translation");
    for(bool radial:{false,true})for(int phase:{-1,1}){
        const auto pure=initialTemplate(*mesh,[&](Vector2D){return Jet{static_cast<double>(phase),0,0,0};},radial);
        const HermiteGalsPhaseRegionSource p(*mesh,pure);const Characteristic m(*mesh,traceDomain(*mesh,radial));
        for(int sign:{-1,1}){const auto vel=Velocity::constant(pure,pure.identity(),{1,0,.02,Position::End},
            radial?Vector2D{0,sign*.25}:Vector2D{sign*.25,0});const Interval extent(*mesh,pure,vel,m,0,.02);
            const auto transfer=integrateHermiteGalsSweptFace(*mesh,extent,faceFor(*mesh,extent,radial,.5,.375),p).transfer;
            check(transfer.signedLiquidVolume()==(phase<0?transfer.signedTotalVolume():0),"exact pure swept identity");++cases;}
    }
    mesh->verts(mesh->vertIndexToLabel[0]).setCoord({-.125,0});
    rejects([&]{integrateHermiteGalsSweptFace(*mesh,interval,row,regions);},"mesh");
}
}
int main(int argc,char** argv){try{if(argc==2 && std::string(argv[1])=="--guards")guards();else{matrix();missedEvents();guards();}
    std::cout<<"I12a7 curved sweep PASS cases="<<cases<<" rejections="<<rejections<<" decimal_bands="<<decimalBands<<'\n';return 0;
}catch(const std::exception& e){std::cerr<<"I12a7 curved sweep FAIL: "<<e.what()<<'\n';return 1;}}
