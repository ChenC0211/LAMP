#include "Modules/Interface/HermiteGalsTemporalPhaseFlux.h"
#include "Core/Geometry/ActiveQuarticHermitePatch.h"
#include "Core/Geometry/LevelSetHermiteGeometry.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <limits>
#include <set>
#include <stdexcept>

namespace
{
using namespace Quantum::Interface;
using Velocity=HermiteGalsTemporalVelocitySource;
using Characteristic=HermiteGalsTemporalCharacteristicSource;
using Interval=HermiteGalsPhaseFluxIntervalSource;
using Position=HermiteGalsTimeLayer::Position;
int cases=0,rejections=0;
void check(bool ok,const char* why){if(!ok)throw std::runtime_error(why);}
void near(long double actual,long double expected,long double scale,const char* why)
{
    if(std::abs(actual-expected)>64*std::numeric_limits<double>::epsilon()*scale) {
        std::cerr<<std::setprecision(18)<<why<<" actual="<<actual<<" expected="<<expected<<'\n';
        throw std::runtime_error(why);
    }
}
template<class F>void rejects(F f,const char* message)
{
    try{f();}catch(const std::exception& e){check(std::string(e.what()).find(message)!=std::string::npos,e.what());++rejections;return;}
    throw std::runtime_error(std::string("missing rejection: ")+message);
}
std::unique_ptr<QuadtreeMesh> meshFor()
{
    DictionaryMesh d;d.Nx=d.Ny=4;d.maxCellNum=128;d.maxInitialIterations=0;
    return std::make_unique<QuadtreeMesh>(d);
}
// Analytic data enter ONLY this old-state fixture, never a time-node rebuild.
HermiteGalsState stateFor(const QuadtreeMesh& mesh,bool radial,double gx,double gy,double c,
    bool perturb=false,HermiteGalsTimeLayer time={})
{
    Quantum::Geometry::ActiveQuarticHermitePatch patch;
    HermiteGalsPrimaryState p;HermiteGalsClosure closure;
    p.measure=radial?HermiteGalsMeasure::AxisymmetricVolume:HermiteGalsMeasure::CartesianArea;
    patch.activeCells.assign(mesh.cellIndexNum(),1);patch.cellCenters.resize(mesh.cellIndexNum());
    const auto value=[&](Vector2D x){return gx*x.x+gy*x.y-c;};
    for(int v=0;v<mesh.vertIndexNum();++v){
        const auto x=mesh.verts(mesh.vertIndexToLabel[v]).coord();
        const double dx=radial && x.x==0?0:gx;
        p.vertices.push_back({v,value(x),dx,gy});patch.nodes.push_back({value(x),dx,gy,0});
    }
    std::set<std::array<int,2>> sides;
    for(int i=0;i<mesh.cellIndexNum();++i){
        const auto& cell=mesh.cells(mesh.cellIndexToLabel[i]);
        p.geometryCells.push_back(i);p.phaseVolumes.push_back({i,0});
        patch.cellCenters[i]=value(cell.centerCoord())+(perturb?std::ldexp(1.,-30):0);
        for(int k=0;k<4;++k){const auto& a=mesh.verts(cell.cellToVert[k]);const auto& b=mesh.verts(cell.cellToVert[(k+1)%4]);
            std::array<int,2> key{{a.index(),b.index()}};
            if(a.coord().y==b.coord().y?a.coord().x>b.coord().x:a.coord().y>b.coord().y)std::swap(key[0],key[1]);
            sides.insert(key);}
    }
    for(const auto& key:sides){
        const auto a=mesh.verts(mesh.vertIndexToLabel[key[0]]).coord(),b=mesh.verts(mesh.vertIndexToLabel[key[1]]).coord();
        const auto x=Vector2D{.5*(a.x+b.x),.5*(a.y+b.y)};
        const double dn=a.y==b.y?gy:(radial && x.x==0?0:gx);
        p.edges.push_back({key,value(x),dn});patch.edgeMidpoints.push_back({key,value(x),dn});
    }
    closure.vertexMixedDerivatives.assign(p.vertices.size(),0);closure.cellBubbleCoefficients=patch.cellCenters;
    HermiteGalsCandidate candidate(mesh,time,p);
    candidate.bindGeometry(mesh,candidate.identity(),closure,LevelSetHermiteGeometry::Snapshot::fromNativeLeaves(mesh,patch,radial));
    return std::move(candidate).finish();
}
LevelSetBoundaryTrace::Domain domainFor(const QuadtreeMesh& mesh,bool radial)
{return {mesh.minCoord,mesh.maxCoord,{.25,.25},radial};}
int faceFor(const QuadtreeMesh& mesh,const Interval& interval,bool horizontal,double normal=.5,double tangent=.375)
{
    for(int row=0;row<mesh.vEdgeIndexNum()+mesh.hEdgeIndexNum();++row){const auto f=interval.face(mesh,row);
        if(f.horizontal==horizontal && (horizontal?f.center.y:f.center.x)==normal &&
            (horizontal?f.center.x:f.center.y)==tangent)return row;}
    throw std::runtime_error("missing native test face");
}
void crossingMatrix()
{
    auto mesh=meshFor();
    for(bool radial:{false,true})for(int band=0;band<(radial?4:1);++band)for(int sign:{-1,1})for(int phase=-4;phase<=68;++phase){
        const double dt=1./64,q=phase/64.;const bool horizontal=radial;
        const double c=sign*.5-q*dt;
        const auto old=stateFor(*mesh,radial,horizontal?0:sign,horizontal?sign:0,c);
        const auto* snapshot=&old.geometry().snapshot();const auto oldId=old.identity();
        const auto velocity=Velocity::constant(old,oldId,{1,0,dt,Position::End},horizontal?Vector2D{0,double(sign)}:Vector2D{double(sign),0});
        const Characteristic characteristic(*mesh,domainFor(*mesh,radial));
        const Interval full(*mesh,old,velocity,characteristic,0,dt),first(*mesh,old,velocity,characteristic,0,dt*.375),last(*mesh,old,velocity,characteristic,dt*.375,dt);
        const int face=faceFor(*mesh,full,horizontal,.5,radial?(band+.5)*.25:.375);
        const auto integral=integrateHermiteGalsPlanarFace(*mesh,full,face);
        const long double r0=band*.25L,r1=(band+1)*.25L;
        const long double area=radial?std::acos(-1.L)*(r1*r1-r0*r0):.25L;
        const long double expected=area*dt*std::clamp(1-q,0.,1.);
        near(integral.signedLiquidVolume(),sign*expected,area*dt,"crossing phase integral");
        near(integral.signedTotalVolume(),sign*area*dt,area*dt,"total integral");
        const auto joined=integrateHermiteGalsPlanarFace(*mesh,first,face).appended(*mesh,integrateHermiteGalsPlanarFace(*mesh,last,face));
        near(joined.signedLiquidVolume(),integral.signedLiquidVolume(),area*dt,"interval composition");
        check(joined.interval().sameInterval(full) && integral.interval().sameInterval(full),"lost interval provenance");
        check(old.identity()==oldId && &old.geometry().snapshot()==snapshot,"changed accepted geometry");
        for(const auto& volume:old.primary().phaseVolumes)check(volume.liquid==0,"changed accepted M");
        ++cases;
    }
    for(double q:{.25,.75}){
        const double dt=.02;const auto old=stateFor(*mesh,false,1,0,.5-q*dt);
        const auto velocity=Velocity::constant(old,old.identity(),{1,0,dt,Position::End},{1,0});
        const Characteristic characteristic(*mesh,domainFor(*mesh,false));const Interval interval(*mesh,old,velocity,characteristic,0,dt);
        const auto t=integrateHermiteGalsPlanarFace(*mesh,interval,faceFor(*mesh,interval,false));
        near(t.signedLiquidVolume(),.25L*(1-q)*dt,.25L*dt,"original decimal crossing");
        std::cout<<std::setprecision(17)<<"ORIGINAL_CROSSING,"<<q<<','<<t.signedLiquidVolume()<<",PASS\n";++cases;
    }
}
void obliqueAndDegenerate()
{
    auto mesh=meshFor();
    for(bool radial:{false,true})for(int sign:{-1,1}){
        // On z=.5, r in [.5,.75], phi=r+z-1 and u_z=1/4:
        // liquid r-.5 < t/4, t in [0,1/2]. Independent triangular moments.
        const auto old=stateFor(*mesh,radial,sign,sign,sign);
        const auto velocity=Velocity::constant(old,old.identity(),{1,0,.5,Position::End},{0,.25});
        const Characteristic characteristic(*mesh,domainFor(*mesh,radial));const Interval interval(*mesh,old,velocity,characteristic,0,.5);
        const int face=faceFor(*mesh,interval,true,.5,.625);
        const auto t=integrateHermiteGalsPlanarFace(*mesh,interval,face);
        const long double triangle=radial?2*std::acos(-1.L)*(.5L*.25L*.5L*.5L/2+.25L*.25L*.5L*.5L*.5L/6):.25L*.5L*.5L/2;
        const long double total=(radial?std::acos(-1.L)*(.75L*.75L-.5L*.5L):.25L)*.5L*.25L;
        near(t.signedLiquidVolume(),sign>0?.25L*triangle:total-.25L*triangle,total,"oblique exact first moment");++cases;
    }
    // A stationary transverse cut is supported; a zero face trace over the
    // entire time interval is not silently assigned to either phase.
    const auto old=stateFor(*mesh,false,0,1,.375);
    auto velocity=Velocity::constant(old,old.identity(),{1,0,.125,Position::End},{1,0});
    const Characteristic characteristic(*mesh,domainFor(*mesh,false));const Interval interval(*mesh,old,velocity,characteristic,0,.125);
    const auto t=integrateHermiteGalsPlanarFace(*mesh,interval,faceFor(*mesh,interval,false));
    near(t.signedLiquidVolume(),.125L*.125L,.25L*.125L,"stationary transverse cut");++cases;
    const auto zero=stateFor(*mesh,false,1,0,.5);
    velocity=Velocity::constant(zero,zero.identity(),{1,0,.125,Position::End},{0,0});
    const Interval zi(*mesh,zero,velocity,characteristic,0,.125);
    rejects([&]{integrateHermiteGalsPlanarFace(*mesh,zi,faceFor(*mesh,zi,false));},"identically zero face-time trace");
    const auto moving=Velocity::constant(zero,zero.identity(),{1,0,.125,Position::End},{1,0});
    const Interval event(*mesh,zero,moving,characteristic,0,.125);
    near(integrateHermiteGalsPlanarFace(*mesh,event,faceFor(*mesh,event,false)).signedLiquidVolume(),.03125,.03125,"endpoint zero event");++cases;
    const auto restarted=stateFor(*mesh,false,1,0,.5,{},{7,1.75,2,Position::End});
    const auto restartV=Velocity::constant(restarted,restarted.identity(),{8,2,2.125,Position::End},{1,0});
    const Interval ri(*mesh,restarted,restartV,characteristic,2.03125,2.125);
    near(integrateHermiteGalsPlanarFace(*mesh,ri,faceFor(*mesh,ri,false)).signedLiquidVolume(),.25L*.09375L,.03125,"restart interval");++cases;
}
void guards()
{
    auto mesh=meshFor();const auto old=stateFor(*mesh,false,1,0,.5),curved=stateFor(*mesh,false,1,0,.5,true);
    const auto domain=domainFor(*mesh,false);const Characteristic characteristic(*mesh,domain);
    const auto query=[](const Vector2D&){return LevelSetBoundaryTrace::VelocityJet{{1,0},{0,0,0,0}};};
    const Velocity opaque(old,old.identity(),{1,0,.125,Position::End},query,query);
    const Interval oi(*mesh,old,opaque,characteristic,0,.125);const int face=faceFor(*mesh,oi,false);
    rejects([&]{integrateHermiteGalsPlanarFace(*mesh,oi,face);},"owned constant velocity");
    const auto cv=Velocity::constant(curved,curved.identity(),{1,0,.125,Position::End},{1,0});
    const Interval ci(*mesh,curved,cv,characteristic,0,.125);
    rejects([&]{integrateHermiteGalsPlanarFace(*mesh,ci,face);},"non-affine native Q4");
    const auto v=Velocity::constant(old,old.identity(),{1,0,.125,Position::End},{1,0});
    const auto traced=traceHermiteGalsWithJacobian(domain,{.5,.375},v.previous(),v.arrival());
    check(traced.departurePoint().x==.375 && traced.departurePoint().y==.375,
        "constant capability differs from the actual retained RK2 views");
    const auto copy=v;check(copy.sameSource(v) && copy.constantValue()->x==1,"constant capability copy lost owner");
    auto unknown=domain;unknown.activeSegment=[](Vector2D,Vector2D){return true;};
    const Characteristic callback(*mesh,unknown);const Interval cb(*mesh,old,v,callback,0,.125);
    rejects([&]{integrateHermiteGalsPlanarFace(*mesh,cb,face);},"opaque coverage");
    const Interval normal(*mesh,old,v,characteristic,0,.125);
    rejects([&]{integrateHermiteGalsPlanarFace(*mesh,normal,faceFor(*mesh,normal,false,0));},"leaves the retained domain");
    const auto radial=stateFor(*mesh,true,0,1,.5);const auto rv=Velocity::constant(radial,radial.identity(),{1,0,.125,Position::End},{1,0});
    const Characteristic rc(*mesh,domainFor(*mesh,true));const Interval r(*mesh,radial,rv,rc,0,.125);
    rejects([&]{integrateHermiteGalsPlanarFace(*mesh,r,faceFor(*mesh,r,true));},"axial translation");
    const auto axial=Velocity::constant(radial,radial.identity(),{1,0,.125,Position::End},{0,1});
    const Interval axis(*mesh,radial,axial,rc,0,.125);
    const auto axisTransfer=integrateHermiteGalsPlanarFace(*mesh,axis,faceFor(*mesh,axis,false,0));
    check(axisTransfer.signedLiquidVolume()==0 && axisTransfer.signedTotalVolume()==0,"nonzero RZ axis transfer");++cases;
    const auto resting=Velocity::constant(old,old.identity(),{1,0,.125,Position::End},{0,0});
    const Interval rest(*mesh,old,resting,characteristic,0,.125);
    check(integrateHermiteGalsPlanarFace(*mesh,rest,faceFor(*mesh,rest,false,.25)).signedTotalVolume()==0,
        "nonzero stationary transfer");++cases;
    rejects([&]{Velocity::constant(old,old.identity(),{1,0,.125,Position::End},{std::numeric_limits<double>::infinity(),0});},"nonfinite constant");
    mesh->verts(mesh->vertIndexToLabel[0]).setCoord({-.125,0});
    rejects([&]{integrateHermiteGalsPlanarFace(*mesh,normal,face);},"mesh");
    check(old.primary().phaseVolumes[0].liquid==0 && old.geometry().identity()==old.identity(),"failure changed old state");
}
}
int main()
{
    try{crossingMatrix();obliqueAndDegenerate();guards();
        std::cout<<"I12a2 planar temporal flux PASS cases="<<cases<<" rejections="<<rejections<<'\n';return 0;}
    catch(const std::exception& e){std::cerr<<"I12a2 FAIL: "<<e.what()<<'\n';return 1;}
}
