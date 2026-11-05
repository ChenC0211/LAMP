#include "Modules/Interface/HermiteGalsTemporalPhaseFlux.h"
#include "Modules/Interface/Detail/HermiteGalsInventoryAssembly.h"
#include "Modules/Interface/HermiteGalsAmrClosure.h"
#include "Core/Geometry/ActiveQuarticHermitePatch.h"
#include "Core/Geometry/LevelSetHermiteGeometry.h"
#include "Core/Mesh/FaceSegmentGraph.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <iomanip>
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
using Velocity=HermiteGalsTemporalVelocitySource;
using Characteristic=HermiteGalsTemporalCharacteristicSource;
using Interval=HermiteGalsPhaseFluxIntervalSource;
using Transfer=HermiteGalsIntegratedFaceTransfer;
using Volumes=HermiteGalsIntegratedFaceVolumes;
using Result=HermiteGalsIntegratedInventoryResult;
using Orientation=HermiteGalsFaceTransferOrientation;
static_assert(!std::is_default_constructible_v<Result>);
static_assert(std::is_same_v<decltype(std::declval<Result&>().faces()),const std::vector<Transfer>&>);
int cases=0,rejections=0;
int increasedPositive=0,increasedNegative=0;
std::string currentCase;
void check(bool ok,const std::string& why){if(!ok)throw std::runtime_error(why);}
void near(long double a,long double b,long double scale,const char* why)
{check(std::abs(a-b)<=64*std::numeric_limits<double>::epsilon()*scale,why);}
template<class F>void rejects(F f,const char* why)
{
    try{f();}catch(const std::exception& e){check(std::string(e.what()).find(why)!=std::string::npos,e.what());++rejections;return;}
    throw std::runtime_error(std::string("missing rejection: ")+why);
}
std::unique_ptr<QuadtreeMesh> makeMesh(int n,bool amr=false)
{DictionaryMesh d;d.Nx=d.Ny=n;d.maxCellNum=16384;d.maxInitialIterations=0;auto mesh=std::make_unique<QuadtreeMesh>(d);
    if(amr){DictionaryMesh::RefineBox b;b.level=1;b.min={.25,.25};b.max={.75,.75};mesh->applyRefineBoxes({b});}return mesh;}
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
HermiteGalsState stateFor(const QuadtreeMesh& mesh,const Profile& f,bool radial,HermiteGalsTimeLayer time,const Inventory& inventory)
{
    HermiteGalsPrimaryState p;p.measure=radial?HermiteGalsMeasure::AxisymmetricVolume:HermiteGalsMeasure::CartesianArea;
    Patch patch;patch.activeCells.assign(mesh.cellIndexNum(),1);patch.nodes.resize(mesh.vertIndexNum());patch.cellCenters.resize(mesh.cellIndexNum());
    for(int v=0;v<mesh.vertIndexNum();++v)if(mesh.verts(mesh.vertIndexToLabel[v]).type()==TVertType::non){const auto j=f(point(mesh,v));patch.nodes[v]=j;p.vertices.push_back({v,j.value,j.dx,j.dy});}
    std::set<Key> keys;
    for(int c=0;c<mesh.cellIndexNum();++c){const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);const auto a=mesh.verts(cell.cellToVert[0]).coord(),b=mesh.verts(cell.cellToVert[2]).coord();
        p.geometryCells.push_back(c);p.phaseVolumes.push_back({c,inventory(c,a,b)});patch.cellCenters[c]=f(a+(b-a)*.5).value;
        for(int k=0;k<4;++k){Key key{{mesh.verts(cell.cellToVert[k]).index(),mesh.verts(cell.cellToVert[(k+1)%4]).index()}};
            const auto x=point(mesh,key[0]),y=point(mesh,key[1]);if(x.y==y.y?x.x>y.x:x.y>y.y)std::swap(key[0],key[1]);keys.insert(key);}}
    for(const auto& key:keys)if(!fineSide(mesh,key)){const auto a=point(mesh,key[0]),b=point(mesh,key[1]);const auto j=f(a+(b-a)*.5);const auto dn=a.y==b.y?j.dy:j.dx;
        p.edges.push_back({key,j.value,dn});patch.edgeMidpoints.push_back({key,j.value,dn});}
    HermiteGalsCandidate candidate(mesh,time,p);const auto mixed=closeHermiteGalsMixedAmr(mesh,candidate);
    HermiteGalsClosure closure;closure.vertexMixedDerivatives=mixed.vertexMixedDerivatives();closure.cellBubbleCoefficients=patch.cellCenters;
    for(std::size_t k=0;k<p.vertices.size();++k)patch.nodes[p.vertices[k].vertex].dxy=closure.vertexMixedDerivatives[k];
    candidate.bindGeometry(mesh,candidate.identity(),closure,LevelSetHermiteGeometry::Snapshot::fromNativeLeaves(mesh,patch,radial));return std::move(candidate).finish();
}

LevelSetBoundaryTrace::Domain traceDomain(const QuadtreeMesh& mesh,bool radial)
{return {mesh.minCoord,mesh.maxCoord,{.125,.125},radial};}
double fraction(Vector2D lo,Vector2D hi,bool horizontal,int sign,double position)
{
    const double a=horizontal?lo.y:lo.x,b=horizontal?hi.y:hi.x;
    return sign>0?std::clamp((position-a)/(b-a),0.,1.):std::clamp((b-position)/(b-a),0.,1.);
}
HermiteGalsState plane(const QuadtreeMesh& mesh,bool radial,int sign,double position,bool empty=false)
{
    const auto domain=domainFor(mesh,radial);
    const Profile profile=[=](Vector2D x){return Jet{sign*((radial?x.y:x.x)-position),
        radial?0.:double(sign),radial?double(sign):0.,0};};
    return stateFor(mesh,profile,radial,{},[&](int cell,Vector2D lo,Vector2D hi){
        return empty?0.:domain.cellVolumes()[domain.cellRowByIndex()[cell]]*fraction(lo,hi,radial,sign,position);});
}
Transfer with(const QuadtreeMesh& mesh,const Transfer& t,const Volumes& v)
{return {mesh,t.interval(),t.face().nativeRow,t.orientation(),v};}
std::vector<Transfer> zeros(const QuadtreeMesh& mesh,const Interval& interval)
{
    const auto domain=domainFor(mesh,interval.measure()==HermiteGalsMeasure::AxisymmetricVolume);
    std::vector<Transfer> result;
    for(const auto f:domain.faces(mesh))result.emplace_back(mesh,interval,f.nativeRow,f.horizontal?Orientation::PositiveY:Orientation::PositiveX,Volumes{});
    return result;
}
std::vector<Transfer> crossingFaces(const QuadtreeMesh& mesh,const Interval& interval,bool radial,int sign,double position)
{
    auto result=zeros(mesh,interval);
    for(auto& t:result){
        const auto& f=t.face();if(f.horizontal!=radial)continue;
        if(f.negativeCell>=0 && f.positiveCell>=0){t=integrateHermiteGalsPlanarFace(mesh,interval,f.nativeRow);continue;}
        // Explicit analytic boundary HISTORY for this manufactured test.
        // I12a2 deliberately supplies no exterior extrapolation/history.
        const double coordinate=radial?f.center.y:f.center.x;
        const long double crossing=sign*(static_cast<long double>(coordinate)-position);
        const long double duration=interval.dt();
        const long double wet=std::max(0.L,static_cast<long double>(interval.timeEnd())-
            std::max(static_cast<long double>(interval.timeBegin()),crossing));
        const long double total=f.area*duration,liquid=f.area*wet;
        t=with(mesh,t,sign>0?Volumes{total,0,liquid,0}:Volumes{0,total,0,liquid});
    }
    return result;
}
void verify(const HermiteGalsState& old,const Result& result,const NativeLeafDomain& domain)
{
    check(result.interval().previousGeometryIdentity()==old.identity(),"lost previous authority");
    std::vector<long double> exchange(old.primary().phaseVolumes.size(),0);
    for(std::size_t row=0;row<result.faces().size();++row){const auto& t=result.faces()[row];
        const auto f=t.face();const auto volume=t.signedLiquidVolume();
        const auto& raw=result.rawFaces()[row].integratedVolumes();const auto& final=t.integratedVolumes();
        check(raw.integratedTotalPositiveVolume==final.integratedTotalPositiveVolume && raw.integratedTotalNegativeVolume==final.integratedTotalNegativeVolume,
            "roundoff changed total volume");
        check(final.integratedLiquidPositiveVolume>=0 && final.integratedLiquidPositiveVolume<=final.integratedTotalPositiveVolume &&
            final.integratedLiquidNegativeVolume>=0 && final.integratedLiquidNegativeVolume<=final.integratedTotalNegativeVolume,"final gross bounds");
        const auto correction=std::max(std::abs(final.integratedLiquidPositiveVolume-raw.integratedLiquidPositiveVolume),
            std::abs(final.integratedLiquidNegativeVolume-raw.integratedLiquidNegativeVolume));
        if(f.negativeCell<0 || f.positiveCell<0)check(correction==0,"roundoff changed prescribed boundary volume");
        else{const auto scale=std::min(domain.cellVolumes()[domain.cellRowByIndex()[f.negativeCell]],domain.cellVolumes()[domain.cellRowByIndex()[f.positiveCell]]);
            check(correction<=64*std::numeric_limits<double>::epsilon()*scale,"face correction exceeds original budget");}
        increasedPositive+=final.integratedLiquidPositiveVolume>raw.integratedLiquidPositiveVolume;
        increasedNegative+=final.integratedLiquidNegativeVolume>raw.integratedLiquidNegativeVolume;
        check(t.interval().sameInterval(result.interval()),"lost final face interval");
        if(f.negativeCell>=0)exchange[f.negativeCell]-=volume;
        if(f.positiveCell>=0)exchange[f.positiveCell]+=volume;
    }
    for(std::size_t k=0;k<result.phaseVolumes().size();++k){
        const auto& v=result.phaseVolumes()[k];const auto cap=domain.cellVolumes()[domain.cellRowByIndex()[v.cell]];
        check(v.liquid>=0 && v.liquid<=cap,"final strict capacity violation");
        const double raw=static_cast<double>(old.primary().phaseVolumes[k].liquid+exchange[v.cell]);
        check(v.liquid==raw,"final inventory did not use shared final volume");
    }
    near(result.report().closure,0,result.report().liquidBefore+1,"global ledger");
}
void crossings()
{
    for(bool radial:{false,true})for(bool amr:{false,true}){
        auto mesh=makeMesh(4,amr);const auto domain=domainFor(*mesh,radial);
        for(int sign:{-1,1})for(int phase=0;phase<=16;++phase){
            currentCase=std::string(radial?"RZ":"XY")+",amr="+std::to_string(amr)+",sign="+std::to_string(sign)+",phase="+std::to_string(phase);
            const double dt=1./64,q=phase/16.,position=.5-sign*q*dt;
            const auto old=plane(*mesh,radial,sign,position);const auto oldM=old.primary().phaseVolumes;
            const auto* geometry=&old.geometry().snapshot();
            const auto velocity=Velocity::constant(old,old.identity(),{1,0,dt,Position::End},radial?Vector2D{0,double(sign)}:Vector2D{double(sign),0});
            const Characteristic characteristic(*mesh,traceDomain(*mesh,radial));
            const Interval full(*mesh,old,velocity,characteristic,0,dt),half(*mesh,old,velocity,characteristic,0,dt*.5),second(*mesh,old,velocity,characteristic,dt*.5,dt);
            const auto inputs=crossingFaces(*mesh,full,radial,sign,position);
            const auto result=[&]{try{return assembleHermiteGalsIntegratedInventory(*mesh,full,inputs);}
                catch(const std::exception&){
                    std::vector<long double> raw(oldM.size());
                    for(const auto& t:inputs){const auto f=t.face();if(f.negativeCell>=0)raw[f.negativeCell]-=t.signedLiquidVolume();if(f.positiveCell>=0)raw[f.positiveCell]+=t.signedLiquidVolume();}
                    std::cerr<<std::setprecision(18);
                    for(std::size_t k=0;k<raw.size();++k)if(oldM[k].liquid+raw[k]<0 || oldM[k].liquid+raw[k]>domain.cellVolumes()[k]){
                        std::cerr<<"RAW_INVENTORY_DEFECT,"<<k<<','<<oldM[k].liquid<<','<<raw[k]<<','<<domain.cellVolumes()[k]
                            <<",defect="<<oldM[k].liquid+raw[k]-domain.cellVolumes()[k]<<'\n';
                        for(const auto& t:inputs)if(t.face().negativeCell==int(k) || t.face().positiveCell==int(k))
                            std::cerr<<"INCIDENT_TRANSFER,"<<t.face().nativeRow<<','<<t.signedTotalVolume()<<','<<t.signedLiquidVolume()<<",native_total="
                                <<(t.face().horizontal==radial?sign*t.face().area*dt:0)<<'\n';
                    }
                    throw;}}();verify(old,result,domain);
            for(const auto& v:result.phaseVolumes()){
                const auto& cell=mesh->cells(mesh->cellIndexToLabel[v.cell]);
                const auto lo=mesh->verts(cell.cellToVert[0]).coord(),hi=mesh->verts(cell.cellToVert[2]).coord();
                const auto cap=domain.cellVolumes()[domain.cellRowByIndex()[v.cell]];
                near(v.liquid,cap*fraction(lo,hi,radial,sign,position+sign*dt),cap,"crossing inventory versus independent end reference");
            }
            auto composed=crossingFaces(*mesh,half,radial,sign,position);
            const auto tail=crossingFaces(*mesh,second,radial,sign,position);
            for(std::size_t k=0;k<composed.size();++k)composed[k]=composed[k].appended(*mesh,tail[k]);
            const auto joined=assembleHermiteGalsIntegratedInventory(*mesh,full,composed);verify(old,joined,domain);
            for(std::size_t k=0;k<result.phaseVolumes().size();++k)
                near(result.phaseVolumes()[k].liquid,joined.phaseVolumes()[k].liquid,domain.cellVolumes()[k],"composed inventory");
            check(&old.geometry().snapshot()==geometry,"changed old snapshot");
            for(std::size_t k=0;k<oldM.size();++k)check(oldM[k].liquid==old.primary().phaseVolumes[k].liquid,"changed old M");
            ++cases;
        }
        std::cout<<"INTEGRATED_CROSSING,"<<(radial?"RZ":"XY")<<','<<amr<<",PASS\n";
    }
}
void guardsAndRepair()
{
    auto mesh=makeMesh(4);const auto domain=domainFor(*mesh,false);
    const double dt=1./64,position=.5-dt*.25;const auto old=plane(*mesh,false,1,position);
    const Characteristic characteristic(*mesh,traceDomain(*mesh,false));
    const auto velocity=Velocity::constant(old,old.identity(),{1,0,dt,Position::End},{1,0});
    const Interval full(*mesh,old,velocity,characteristic,0,dt);
    const auto input=crossingFaces(*mesh,full,false,1,position);
    auto missing=input;missing.pop_back();
    rejects([&]{assembleHermiteGalsIntegratedInventory(*mesh,full,missing);},"missing/duplicate/unordered");
    auto duplicate=input;duplicate[1]=duplicate[0];
    rejects([&]{assembleHermiteGalsIntegratedInventory(*mesh,full,duplicate);},"missing/duplicate/unordered");
    auto extra=input;extra.push_back(input.back());
    rejects([&]{assembleHermiteGalsIntegratedInventory(*mesh,full,extra);},"foreign/inactive");
    const auto foreignV=Velocity::constant(old,old.identity(),{1,0,dt,Position::End},{1,0});
    const Interval foreign(*mesh,old,foreignV,characteristic,0,dt);
    auto mixed=input;mixed[0]=Transfer(*mesh,foreign,input[0].face().nativeRow,input[0].orientation(),input[0].integratedVolumes());
    rejects([&]{assembleHermiteGalsIntegratedInventory(*mesh,full,mixed);},"different interval source");
    const Interval late(*mesh,old,velocity,characteristic,dt*.5,dt);
    rejects([&]{assembleHermiteGalsIntegratedInventory(*mesh,late,zeros(*mesh,late));},"begin at the old accepted time");
    auto large=input;
    for(auto& t:large){auto v=t.integratedVolumes();v.integratedTotalPositiveVolume*=32;v.integratedLiquidPositiveVolume*=32;t=with(*mesh,t,v);}
    rejects([&]{assembleHermiteGalsIntegratedInventory(*mesh,full,large);},"gross outflow CFL");
    std::size_t interior=0,downstream=0;
    for(std::size_t k=0;k<input.size();++k){const auto f=input[k].face();
        if(!f.horizontal && f.negativeCell>=0 && f.positiveCell>=0){interior=k;if(f.center.x==.75)downstream=k;}}
    auto divergent=input;auto v=divergent[interior].integratedVolumes();v.integratedTotalPositiveVolume+=.001;
    divergent[interior]=with(*mesh,divergent[interior],v);
    rejects([&]{assembleHermiteGalsIntegratedInventory(*mesh,full,divergent);},"not divergence free");
    auto invalidM=input;v=invalidM[downstream].integratedVolumes();v.integratedLiquidPositiveVolume=v.integratedTotalPositiveVolume;
    invalidM[downstream]=with(*mesh,invalidM[downstream],v);
    rejects([&]{assembleHermiteGalsIntegratedInventory(*mesh,full,invalidM);},"inventory violation exceeds roundoff budget");
    for(double tolerance:{0.,-1.,1.,std::numeric_limits<double>::quiet_NaN()}){
        rejects([&]{assembleHermiteGalsIntegratedInventory(*mesh,full,input,{tolerance,1e-12});},"invalid integrated inventory tolerance");
        rejects([&]{assembleHermiteGalsIntegratedInventory(*mesh,full,input,{1e-12,tolerance});},"invalid integrated inventory tolerance");
    }
    const auto empty=plane(*mesh,false,1,-1,true);
    const auto ev=Velocity::constant(empty,empty.identity(),{1,0,dt,Position::End},{0,0});
    const Interval ei(*mesh,empty,ev,characteristic,0,dt);auto roundoff=zeros(*mesh,ei);
    const auto f=roundoff[interior].face();const auto cap=domain.cellVolumes()[domain.cellRowByIndex()[f.negativeCell]];
    const long double tiny=32*std::numeric_limits<double>::epsilon()*cap;
    roundoff[interior]=with(*mesh,roundoff[interior],{.001,.001,tiny,0});
    const auto repaired=assembleHermiteGalsIntegratedInventory(*mesh,ei,roundoff);verify(empty,repaired,domain);
    check(repaired.report().repairedFaces==1 && repaired.report().roundoffRepairSweeps==1,"missing shared roundoff repair");
    check(repaired.rawFaces()[interior].signedLiquidVolume()==tiny && repaired.faces()[interior].signedLiquidVolume()==0,"raw/final repair provenance");
    check(roundoff[interior].signedLiquidVolume()==tiny,"assembly changed caller input");
    auto excessive=roundoff;excessive[interior]=with(*mesh,excessive[interior],{.001,.001,4*tiny,0});
    rejects([&]{assembleHermiteGalsIntegratedInventory(*mesh,ei,excessive);},"inventory violation exceeds roundoff budget");
    auto counter=zeros(*mesh,ei);
    for(auto& t:counter)if(t.face().negativeCell>=0 && t.face().positiveCell>=0)t=with(*mesh,t,{.001,.001,.0005,.0005});
    const auto cf=assembleHermiteGalsIntegratedInventory(*mesh,ei,counter);verify(empty,cf,domain);
    check(cf.report().counterflowFaces>0 && cf.report().maxOutflowCfl>0 && cf.report().maxRelativeDivergence==0,"counterflow hidden by signed cancellation");
    auto moved=cf;auto retained=std::move(moved);
    rejects([&]{moved.faces();},"moved-from integrated inventory");
    check(retained.interval().sameInterval(ei),"lost retained interval");
    mesh->verts(mesh->vertIndexToLabel[0]).setCoord({-.125,0});
    rejects([&]{assembleHermiteGalsIntegratedInventory(*mesh,full,input);},"mesh");
    rejects([&]{retained.validateMesh(*mesh);},"mesh");
    cases+=2;
    std::cout<<"INVENTORY_GUARDS,COUNTERFLOW,SHARED_ROUNDOFF,PASS\n";
}
void availableGrossRepair()
{
    auto mesh=makeMesh(4);const auto domain=domainFor(*mesh,false);
    for(int sign:{-1,1}){
        // Fixed boundary liquid history, partial internal liquid gross.
        // sign+ fills a full cell by one ulp; sign- empties an empty cell.
        const auto old=plane(*mesh,false,1,sign>0?2:-1,sign<0);
        const auto velocity=Velocity::constant(old,old.identity(),{1,0,.125,Position::End},{0,0});
        const Characteristic characteristic(*mesh,traceDomain(*mesh,false));const Interval interval(*mesh,old,velocity,characteristic,0,.125);
        auto input=zeros(*mesh,interval);const long double tiny=std::ldexp(1.L,-54);
        for(auto& t:input)if(!t.face().horizontal && t.face().center.y==.375){
            const auto liquid=.0005L-(t.face().center.x==.5?tiny:0);
            t=with(*mesh,t,sign>0?Volumes{.001,0,liquid,0}:Volumes{0,.001,0,liquid});}
        const auto result=assembleHermiteGalsIntegratedInventory(*mesh,interval,input);verify(old,result,domain);
        check(result.report().repairedFaces>0,"available gross repair was not exercised");
        bool increased=false;for(std::size_t k=0;k<input.size();++k)
            increased=increased || (sign>0?result.faces()[k].integratedVolumes().integratedLiquidPositiveVolume>input[k].integratedVolumes().integratedLiquidPositiveVolume
                :result.faces()[k].integratedVolumes().integratedLiquidNegativeVolume>input[k].integratedVolumes().integratedLiquidNegativeVolume);
        check(increased,"available internal gross component was not used");++cases;
    }
    std::cout<<"AVAILABLE_GROSS_REPAIR,BOTH_SIGNS,FIXED_BOUNDARY,PASS\n";
}
void decimalCrossings()
{
    for(bool radial:{false,true}){
        // Retain the original h=.25 crossing cells around coordinate .5.
        // The normal domain [.25,.75] keeps the old affine scalar exactly
        // representable throughout its entire retained support; boundary
        // histories remain explicit. No new-time scalar is constructed.
        DictionaryMesh d;d.Nx=radial?4:2;d.Ny=radial?2:4;d.maxCellNum=128;d.maxInitialIterations=0;
        d.minCoord=radial?Vector2D{0,.25}:Vector2D{.25,0};d.maxCoord=radial?Vector2D{1,.75}:Vector2D{.75,1};
        QuadtreeMesh mesh(d);const auto domain=domainFor(mesh,radial);
        for(int sign:{-1,1})for(double q:{.25,.75}){
            currentCase=std::string("decimal ")+(radial?"RZ":"XY")+",sign="+std::to_string(sign)+",q="+std::to_string(q);
            const double dt=.02,position=.5-sign*q*dt;const auto old=plane(mesh,radial,sign,position);
            const auto velocity=Velocity::constant(old,old.identity(),{1,0,dt,Position::End},radial?Vector2D{0,double(sign)}:Vector2D{double(sign),0});
            const Characteristic characteristic(mesh,traceDomain(mesh,radial));const Interval interval(mesh,old,velocity,characteristic,0,dt);
            const auto input=crossingFaces(mesh,interval,radial,sign,position);
            const auto result=assembleHermiteGalsIntegratedInventory(mesh,interval,input);verify(old,result,domain);
            for(const auto& v:result.phaseVolumes()){
                const auto& cell=mesh.cells(mesh.cellIndexToLabel[v.cell]);
                const auto lo=mesh.verts(cell.cellToVert[0]).coord(),hi=mesh.verts(cell.cellToVert[2]).coord();
                const auto cap=domain.cellVolumes()[domain.cellRowByIndex()[v.cell]];
                near(v.liquid,cap*fraction(lo,hi,radial,sign,position+sign*dt),cap,"decimal crossing inventory");}
            if(!radial)for(const auto& t:input)if(!t.face().horizontal && t.face().center.x==.5)
                near(t.signedLiquidVolume(),sign*.25L*(1-q)*dt,.25L*dt,"original decimal face integral");
            ++cases;
        }
    }
    std::cout<<"DECIMAL_CROSSING,XY_RZ,BOTH_SIGNS,Q_1_4_3_4,PASS\n";
}
void certifiedRoundoffRepair()
{
    using Face=Detail::InventoryFace;using Bounds=Detail::InventoryFaceBounds;using Assembly=Detail::HermiteInventoryAssembly;
    const std::vector<int> cells{0,1};int repaired=0,oldEscapes=0,growth=0;const int oldRejections=rejections;
    for(double capacity:{1.,.00390625})for(bool full:{false,true})for(int sign:{-1,1})for(bool split:{false,true})
    for(double baseFraction:{0.,.03125})for(int exponent:{-49,-90}){
        if(full && exponent==-90)continue;
        const double total=.125*capacity,base=baseFraction*capacity,tiny=std::ldexp(capacity,exponent);
        if(base+tiny==base)continue;
        const std::vector<double> capacities{capacity,capacity};
        const std::vector<HermiteGalsCellVolume> old{{0,full?capacity:0.},{1,full?capacity:0.}};
        const Assembly assembly(2,cells,capacities,old,{});
        std::vector<Face> input{{-1,0,total,total,0,0,0,0},{1,-1,total,total,0,0,0,0}};
        std::vector<Bounds> limits(2);
        for(int piece=0;piece<(split?2:1);++piece){
            const auto positive=base+(sign>0?tiny:0),negative=base+(sign<0?tiny:0);
            input.push_back({0,1,total,total,positive,negative,positive,negative});
            limits.push_back({base,positive,base,negative});
        }
        auto prior=input;std::vector<HermiteGalsCellVolume> priorVolumes;
        const auto priorReport=assembly.finish(prior,1,1e-12,1e-12,priorVolumes);
        bool escaped=false;for(std::size_t k=0;k<prior.size();++k)escaped=escaped
            || prior[k].liquidPositive>limits[k].upperPositive || prior[k].liquidNegative>limits[k].upperNegative;
        check(escaped && priorReport.repairedFaces>0,"fixture did not reproduce original certificate escape");++oldEscapes;
        auto faces=input;std::vector<HermiteGalsCellVolume> volumes;
        const auto report=assembly.finish(faces,1,1e-12,1e-12,volumes,&limits);
        check(report.roundoffRepairSweeps>0 && report.repairedFaces>0 && report.closure==0
            && report.maxRelativeFluxCorrection<=Detail::hermiteInventoryRoundoffBudget,"certified repair lost old budget/ledger");
        for(std::size_t k=0;k<faces.size();++k){const auto& f=faces[k];const auto& b=limits[k];
            check(f.liquidPositive>=b.lowerPositive && f.liquidPositive<=b.upperPositive
                && f.liquidNegative>=b.lowerNegative && f.liquidNegative<=b.upperNegative,"repair escaped certificate");
            check(f.totalPositive==input[k].totalPositive && f.totalNegative==input[k].totalNegative,"repair changed total");
            if(k<2)check(f.liquidPositive==input[k].liquidPositive && f.liquidNegative==input[k].liquidNegative,"repair changed boundary");
            else check(f.liquidPositive==base && f.liquidNegative==base,"repair failed to preserve nonzero lower bound");
        }
        check(volumes.size()==2 && volumes[0].liquid==old[0].liquid && volumes[1].liquid==old[1].liquid,"shared repaired flux did not restore inventory");
        auto again=input;std::vector<HermiteGalsCellVolume> repeated;
        const auto repeat=assembly.finish(again,1,1e-12,1e-12,repeated,&limits);
        check(repeat.roundoffRepairSweeps==report.roundoffRepairSweeps && repeat.repairedFaces==report.repairedFaces,"repair nondeterministic");
        for(std::size_t k=0;k<faces.size();++k)check(again[k].liquidPositive==faces[k].liquidPositive
            && again[k].liquidNegative==faces[k].liquidNegative,"repeat selected different flux");++repaired;
    }
    for(bool full:{false,true})for(int sign:{-1,1}){
        const std::vector<double> capacity{1,1};const std::vector<HermiteGalsCellVolume> old{{0,full?1.:0.},{1,full?1.:0.}};
        const Assembly assembly(2,cells,capacity,old,{});const double base=.03125,tiny=std::ldexp(1.,-49);
        const double positive=base+(sign>0?tiny:0),negative=base+(sign<0?tiny:0);
        std::vector<Face> faces{{0,1,.125,.125,positive,negative,positive,negative}};
        const std::vector<Bounds> limits{{base,base+tiny,base,base+tiny}};
        std::vector<HermiteGalsCellVolume> volumes;const auto report=assembly.finish(faces,1,1e-12,1e-12,volumes,&limits);
        check(faces[0].liquidPositive==base+tiny && faces[0].liquidNegative==base+tiny
            && report.repairedFaces==1 && report.closure==0,"legal certified growth unavailable");++growth;
    }
    const std::vector<double> capacities{1,1};const std::vector<HermiteGalsCellVolume> old{{0,0},{1,0}};
    const Assembly assembly(2,cells,capacities,old,{});const double tiny=std::ldexp(1.,-49);
    const std::vector<Face> initial{{0,1,.125,.125,tiny,0,tiny,0}};
    const std::vector<Bounds> fixed{{tiny,tiny,0,0}},valid{{0,tiny,0,0}};
    const auto attempt=[&](const std::vector<Bounds>& limits,long double unit=1){auto faces=initial;std::vector<HermiteGalsCellVolume> volumes;
        return assembly.finish(faces,unit,1e-12,1e-12,volumes,&limits);};
    rejects([&]{attempt(fixed);},"no certified flux capacity");
    rejects([&]{attempt({});},"layout/units");
    rejects([&]{attempt(valid,2);},"layout/units");
    rejects([&]{attempt({{0,0,0,0}});},"certified phase bounds");
    rejects([&]{attempt({{-tiny,tiny,0,0}});},"certified phase bounds");
    rejects([&]{attempt({{0,.25,0,0}});},"certified phase bounds");
    rejects([&]{attempt({{0,std::numeric_limits<long double>::infinity(),0,0}});},"certified phase bounds");
    rejects([&]{auto faces=initial;faces[0].liquidPositive=faces[0].rawPositive=4*Detail::hermiteInventoryRoundoffBudget;
        const std::vector<Bounds> limits{{0,faces[0].liquidPositive,0,0}};std::vector<HermiteGalsCellVolume> volumes;
        assembly.finish(faces,1,1e-12,1e-12,volumes,&limits);},"inventory violation exceeds roundoff budget");
    check(repaired==40 && oldEscapes==40 && growth==4,"incomplete certified repair matrix");
    std::cout<<"I12a16r7 certified roundoff PASS reductions="<<repaired<<" original_escapes="<<oldEscapes
        <<" growth="<<growth<<" guards="<<rejections-oldRejections<<'\n';
}
}
int main(int argc,char** argv)
{
    try{
        if(argc==2 && std::string(argv[1])=="--algebra-contract")guardsAndRepair();
        else{check(argc==1,"unknown integrated inventory argument");crossings();guardsAndRepair();availableGrossRepair();decimalCrossings();}
        std::cout<<"INCREASED_COMPONENTS,"<<increasedPositive<<','<<increasedNegative<<'\n';
        std::cout<<"I12a3 integrated inventory PASS cases="<<cases<<" rejections="<<rejections<<'\n';
        if(argc==1)certifiedRoundoffRepair();return 0;}
    catch(const std::exception& e){std::cerr<<"I12a3 FAIL: "<<currentCase<<' '<<e.what()<<'\n';return 1;}
}
