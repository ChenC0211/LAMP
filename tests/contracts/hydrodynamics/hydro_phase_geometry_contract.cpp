#include "Modules/Hydrodynamics/HydroPhaseGeometry.h" // Public header stands alone.
#include "Modules/Hydrodynamics/HydroSpatialAdvector.h"
#include "Modules/Hydrodynamics/HydroMomentumGeometry.h"
#include "Modules/Hydrodynamics/HydroVofAdvector.h"
#include <sstream>
#include "Modules/Hydrodynamics/Detail/HydroRegularAdvectorGeometry.h"
#include "Modules/Hydrodynamics/Detail/HydroNormalRatePolynomial.h"
#include <chrono>
#include "interface_first_jet_fixture.h"
#include "Core/Mesh/CartesianMeasure.h"
#include <map>
#include <cfenv>
#include <optional>
#include <type_traits>

namespace
{
using namespace LAMP::Hydrodynamics;
int scenarios=0, faceChecks=0, cellChecks=0, coarseChecks=0, rejectionChecks=0, convergenceChecks=0;
long double maximumCell=0, maximumFace=0, maximumRate=0, maximumClosure=0;
static_assert(!std::is_default_constructible_v<HydroPhaseGeometry>);
static_assert(!std::is_default_constructible_v<HydroPhaseRate>);
static_assert(std::is_same_v<decltype(std::declval<HydroPhaseGeometry&>().cellVolumes()),
    const std::vector<HydroPhaseMeasure>&>);

HydroPhysicalDivergence metricsFor(const QuadtreeMesh& mesh,bool radial,bool subset=false)
{
    FaceSegmentGraph graph;
    graph.build(mesh,[subset](const QuadtreeMesh&,const Cell& cell){return !subset || cell.centerCoord().x<.5;});
    return HydroPhysicalDivergence::build(mesh,HydroIntegralDivergence::build(mesh,
        HydroDofLayout::build(mesh,NativeLeafDomain::build(mesh,graph,radial))));
}
HermiteGalsPublication publishAt(const QuadtreeMesh& mesh,const Case& s,HermiteGalsTimeLayer time,bool subset=false)
{
    const auto f=s.profile(time.time());
    const auto seed=stateFor(mesh,f,s.radial,time,[](int,Vector2D,Vector2D){return 0.;},
        s.symmetry(mesh),HermiteGalsTransport::FirstJet);
    auto primary=seed.primary();
    if(subset) {
        primary.geometryCells.erase(std::remove_if(primary.geometryCells.begin(),primary.geometryCells.end(),
            [&](int c){return mesh.cells(mesh.cellIndexToLabel[c]).centerCoord().x>=.75;}),primary.geometryCells.end());
        std::set<Key> activeSides;
        for(int c:primary.geometryCells) {
            const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);
            for(int k=0;k<4;++k) {
                Key key{mesh.verts(cell.cellToVert[k]).index(),mesh.verts(cell.cellToVert[(k+1)%4]).index()};
                const auto a=point(mesh,key[0]),b=point(mesh,key[1]);
                if(a.y==b.y?a.x>b.x:a.y>b.y)std::swap(key[0],key[1]);activeSides.insert(key);
            }
        }
        primary.edges.erase(std::remove_if(primary.edges.begin(),primary.edges.end(),
            [&](const auto& e){return activeSides.count(e.vertices)==0;}),primary.edges.end());
        std::set<int> activeVertices;for(const auto& e:activeSides){activeVertices.insert(e[0]);activeVertices.insert(e[1]);}
        primary.vertices.erase(std::remove_if(primary.vertices.begin(),primary.vertices.end(),
            [&](const auto& v){return activeVertices.count(v.vertex)==0;}),primary.vertices.end());
    }
    HermiteGalsCandidate candidate(mesh,time,primary);
    const auto mixed=closeHermiteGalsMixedAmr(mesh,candidate,subset ? std::vector<LevelSetInterfacePanels::SymmetryPlane>{} : s.symmetry(mesh));
    std::vector<double> centerValues;for(int c:primary.geometryCells)centerValues.push_back(f(mesh.cells(mesh.cellIndexToLabel[c]).centerCoord()).value);
    const HermiteGalsCenterPrediction center{candidate.identity(),centerValues};
    HermiteGalsPublicationConfig config;config.volume.symmetryPlanes=s.symmetry(mesh);
    if(subset) {
        check(!s.radial && s.planar && s.axial(),"partial fixture is an explicitly open horizontal plane");
        config.volume.symmetryPlanes.clear();
        config.volume.boundaryContacts={{0,s.planeOffset},{.75,s.planeOffset}};
    }
    return publishHermiteGalsFirstJetGeometry(mesh,std::move(candidate),mixed,center,config);
}
// Independent exact interval/disk intersection, not a query of the production
// scalar or mapped phase ranges. Every reference uses the continuum Case.
long double intervalArea(long double a,long double b,long double normal,bool horizontal,const Case& s)
{
    const long double capacity=(b-a)*(s.radial ? CartesianMeasure::twoPi*(horizontal ? (a+b)/2 : normal) : 1.L);
    long double lo=a,hi=b;
    if(s.planar) {
        const bool planeAlongTangent=horizontal!=s.axial();
        if(planeAlongTangent) {
            if(s.planeSign>0)hi=std::min(hi,static_cast<long double>(s.planeOffset));
            else lo=std::max(lo,static_cast<long double>(s.planeOffset));
        }else if(s.planeSign*(normal-s.planeOffset)>=0)return 0;
    }else {
        const long double normalCenter=horizontal ? s.cy : (s.radial ? 0 : s.cx);
        const long double tangentCenter=horizontal ? (s.radial ? 0 : s.cx) : s.cy;
        const long double delta=normal-normalCenter, q=static_cast<long double>(s.radius)*s.radius-delta*delta;
        if(q<=0)return 0;
        const long double width=std::sqrt(q);
        lo=std::max(lo,tangentCenter-width);hi=std::min(hi,tangentCenter+width);
    }
    if(hi<=lo)return 0;
    const long double result=(hi-lo)*(s.radial ? CartesianMeasure::twoPi*(horizontal ? (lo+hi)/2 : normal) : 1.L);
    check(result>=0 && result<=capacity*(1+1e-14L),"invalid independent interval reference");
    return result;
}
long double exactFace(const QuadtreeMesh& mesh,const NativeLeafFace& face,const Case& s)
{
    const auto& edge=face.horizontal ? mesh.hEdges(face.edgeLabel) : mesh.vEdges(face.edgeLabel);
    auto a=mesh.verts(edge.edgeToVert[0]).coord(), b=mesh.verts(edge.edgeToVert[1]).coord();
    if(face.horizontal ? a.x>b.x : a.y>b.y)std::swap(a,b);
    return intervalArea(face.horizontal ? a.x : a.y,face.horizontal ? b.x : b.y,
        face.horizontal ? a.y : a.x,face.horizontal,s);
}
long double exactCell(Vector2D a,Vector2D b,const Case& s)
{
    if(!s.planar)return circleCell(a,b,{s.radial?0:s.cx,s.cy},s.radius,s.radial);
    if(s.axial()) {
        if(s.planeSign>0)b.y=std::min(b.y,s.planeOffset);else a.y=std::max(a.y,s.planeOffset);
    }else {
        if(s.planeSign>0)b.x=std::min(b.x,s.planeOffset);else a.x=std::max(a.x,s.planeOffset);
    }
    if(b.x<=a.x || b.y<=a.y)return 0;
    return (static_cast<long double>(b.x)-a.x)*(static_cast<long double>(b.y)-a.y)*
        (s.radial ? CartesianMeasure::twoPi*(static_cast<long double>(b.x)+a.x)/2 : 1.L);
}
std::vector<double> velocityFor(const QuadtreeMesh& mesh,const HydroDofLayout& layout)
{
    // Declared physical-mean affine/even-r^2 family, independently tested below.
    std::vector<double> u;
    for(int i=0;i<layout.velocityCount();++i){const auto f=layout.domain().face(mesh,layout.nativeFace(i));
        if(!layout.domain().axisymmetric())u.push_back(f.horizontal?-.2+.31*f.center.x+.19*f.center.y:.3+.17*f.center.x-.23*f.center.y);
        else if(!f.horizontal)u.push_back(.13*f.center.x);
        else{const auto& e=mesh.hEdges(f.edgeLabel);const auto a=mesh.verts(e.edgeToVert[0]).coord(),b=mesh.verts(e.edgeToVert[1]).coord();u.push_back(.2+.155*.5*(a.x*a.x+b.x*b.x)+.19*f.center.y);}}
    return u;
}
HydroPhaseRate phaseRate(const QuadtreeMesh& mesh,HydroPhaseGeometry g,const HydroNativeValues& values,
    HermiteGalsIdentity id,HermiteGalsTimeLayer time)
{
    const auto geometry=HydroMomentumGeometry::build(mesh,g.metrics());
    const auto advector=HydroAdvectorReconstruction::prepare(mesh,geometry).evaluate(mesh,values,id,time);
    return HydroPhaseRate::prepareSpatial(mesh,g,advector);
}
HydroPhaseMeasure spatialReference(const HydroAdvectorAtom&,const Case&,int);
struct Error {long double cell=0,face=0;};
Error inspect(const QuadtreeMesh& mesh,const Case& s,const HydroPhaseGeometry& geometry,bool exact)
{
    const auto& layout=geometry.layout();
    const auto u=velocityFor(mesh,layout);
    const HydroNativeValues values(layout,u,std::vector<double>(layout.pressureCount(),0));
    const auto rate=phaseRate(mesh,geometry,values,geometry.geometryIdentity(),geometry.timeLayer());
    std::vector<HydroPhaseMeasure> cellReference(layout.pressureCount());
    HydroPhaseMeasure exterior;
    long double totalArea=0,totalVolume=0;
    Error error;
    std::map<std::pair<int,int>,std::vector<int>> sides;
    for(int i=0;i<layout.velocityCount();++i) {
        const auto face=layout.domain().face(mesh,layout.nativeFace(i));
        const auto measured=geometry.faceMeasures()[i], q=rate.faceRates()[i];
        const long double ref=exactFace(mesh,face,s);
        const long double difference=std::abs(measured.liquid-ref);
        if(exact) {
            const auto e=close(measured.liquid,ref,1e-10L*face.area,"analytic face area");
            maximumFace=std::max(maximumFace,face.area>0?e/face.area:0);
        }
        if(exact){const auto& edge=face.horizontal?mesh.hEdges(face.edgeLabel):mesh.vEdges(face.edgeLabel);auto a=mesh.verts(edge.edgeToVert[0]).coord(),b=mesh.verts(edge.edgeToVert[1]).coord();
            if(face.horizontal?a.x>b.x:a.y>b.y)std::swap(a,b);
            const auto expected=spatialReference(HydroAdvectorAtom{face.horizontal?1:0,-1,-1,a,b,face.area},s,2);
            maximumRate=std::max({maximumRate,close(q.liquid,expected.liquid,1e-10*face.area,"independent spatial liquid rate"),close(q.gas,expected.gas,1e-10*face.area,"independent spatial gas rate")});}
        error.face+=difference;totalArea+=face.area;++faceChecks;
        check(rate.isSpatial(),"old face-constant owner used by migrated contract");
        maximumClosure=std::max(maximumClosure,close(q.liquid+q.gas,rate.totalFaceRates()[i],
            1e-12*std::max(1.,std::abs(rate.totalFaceRates()[i])),"phase total rate"));
        for(int sign : {1,-1}) {
            const int native=sign>0 ? face.negativeCell : face.positiveCell;
            if(native>=0) {
                const int p=layout.pressureDof(native);
                cellReference[p].liquid+=sign*q.liquid;cellReference[p].gas+=sign*q.gas;
                sides[{p,(face.horizontal?2:0)+(sign<0?1:0)}].push_back(i);
            }
        }
        if(face.outwardSign!=0){exterior.liquid+=face.outwardSign*q.liquid;exterior.gas+=face.outwardSign*q.gas;}
    }
    HydroPhaseMeasure sum;
    for(int p=0;p<layout.pressureCount();++p) {
        const auto bounds=geometry.publication().state().geometry().snapshot().nativeCellBounds(layout.nativeCell(p));
        const auto m=geometry.cellVolumes()[p];const long double ref=exactCell(bounds[0],bounds[1],s);
        const double volume=geometry.metrics().cellVolumes()[p];
        if(exact)maximumCell=std::max(maximumCell,close(m.liquid,ref,1e-10*volume,"analytic cell volume")/volume);
        error.cell+=std::abs(m.liquid-ref);totalVolume+=volume;++cellChecks;
        close(m.liquid+m.gas,volume,1e-13*volume,"cell capacity");
        close(geometry.liquidFractions()[p],m.liquid/volume,1e-15,"alpha is geometric volume fraction");
        const auto q=rate.cellOutwardRates()[p];
        maximumClosure=std::max(maximumClosure,close(q.liquid,cellReference[p].liquid,1e-12,"native cell liquid outflow"));
        maximumClosure=std::max(maximumClosure,close(q.gas,cellReference[p].gas,1e-12,"native cell gas outflow"));
        sum.liquid+=q.liquid;sum.gas+=q.gas;
    }
    maximumClosure=std::max(maximumClosure,close(sum.liquid,exterior.liquid,1e-12,"internal liquid cancellation"));
    maximumClosure=std::max(maximumClosure,close(sum.gas,exterior.gas,1e-12,"internal gas cancellation"));
    for(const auto& [key,faces] : sides)if(faces.size()>1) {
        const int p=key.first;const auto bounds=geometry.publication().state().geometry().snapshot().nativeCellBounds(layout.nativeCell(p));
        long double measured=0,ref=0;
        for(int i:faces){measured+=geometry.faceMeasures()[i].liquid;ref+=exactFace(mesh,layout.domain().face(mesh,layout.nativeFace(i)),s);}
        if(exact)close(measured,ref,1e-10,"actual fine subface aggregation");
        // This comparison uses the whole physical side, independently of its
        // native subface split, including RZ radial variation on horizontal sides.
        const bool horizontal=key.second>=2;
        const double normal=horizontal ? (key.second==2?bounds[1].y:bounds[0].y) :
            (key.second==0?bounds[1].x:bounds[0].x);
        const auto whole=intervalArea(horizontal?bounds[0].x:bounds[0].y,
            horizontal?bounds[1].x:bounds[1].y,normal,horizontal,s);
        close(ref,whole,1e-12,"independent coarse-side measure");++coarseChecks;
    }
    auto reverse=u;for(auto& value:reverse)value=-value;
    const auto negative=phaseRate(mesh,geometry,HydroNativeValues(layout,reverse,
        std::vector<double>(layout.pressureCount(),0)),geometry.geometryIdentity(),geometry.timeLayer());
    for(int i=0;i<layout.velocityCount();++i) {
        close(negative.faceRates()[i].liquid,-rate.faceRates()[i].liquid,0,"velocity reversal liquid");
        close(negative.faceRates()[i].gas,-rate.faceRates()[i].gas,0,"velocity reversal gas");
    }
    check(geometry.publication().state().primary().phaseVolumes.empty(),"material consumption restored Interface M");
    ++scenarios;
    std::cout<<"GEOMETRY,"<<s.radial<<','<<s.amr<<','<<s.planar<<','<<s.nonpolynomial<<','<<s.planeOffset<<','
        <<layout.pressureCount()<<','<<layout.velocityCount()<<','<<geometry.report().cutCells<<','<<geometry.report().splitFaces<<','
        <<geometry.report().maxRelativeFaceClosure<<'\n';
    error.cell/=totalVolume;error.face/=totalArea;return error;
}
template<class F>void reject(F f,const std::string& message)
{rejects(f,message);++rejectionChecks;}
void rejections()
{
    auto mesh=makeMesh(4);Case s;s.planar=true;s.planeOffset=.43;s.pureFirstJet=true;
    const auto pub=publishAt(*mesh,s,{});
    const auto geometry=HydroPhaseGeometry::build(*mesh,metricsFor(*mesh,false),pub);
    const auto& layout=geometry.layout();
    const auto u=velocityFor(*mesh,layout);
    const HydroNativeValues values(layout,u,std::vector<double>(layout.pressureCount()));
    const auto evaluate=[&](auto id,auto time){return phaseRate(*mesh,geometry,values,id,time);};
    reject([&]{evaluate(HermiteGalsIdentity{},geometry.timeLayer());},"unbound provenance");
    const auto other=publishAt(*mesh,s,{});
    reject([&]{evaluate(other.state().identity(),geometry.timeLayer());},"different geometry identity");
    // Same physical time, different clock fields: cannot alias a candidate.
    for(auto t : {HermiteGalsTimeLayer{1,0,0,Position::Initial},HermiteGalsTimeLayer{0,0,0,Position::End},
        HermiteGalsTimeLayer{1,-1,1,Position::Midpoint},HermiteGalsTimeLayer{0,0,1,Position::Initial}})
        reject([&]{evaluate(geometry.geometryIdentity(),t);},"time");
    auto nanTime=geometry.timeLayer();nanTime.end=std::numeric_limits<double>::quiet_NaN();
    reject([&]{evaluate(geometry.geometryIdentity(),nanTime);},"time");
    const auto newMetrics=metricsFor(*mesh,false);
    reject([&]{phaseRate(*mesh,geometry,HydroNativeValues(newMetrics.layout(),u,
        std::vector<double>(layout.pressureCount())),geometry.geometryIdentity(),geometry.timeLayer());},"layout identity");
    reject([&]{HydroPhaseGeometry::build(*mesh,metricsFor(*mesh,true),pub);},"different XY/RZ measure");
    Case partialShape=s;partialShape.axialTranslation=true;partialShape.planeOffset=.4375;
    const auto partial=publishAt(*mesh,partialShape,{},true);
    reject([&]{HydroPhaseGeometry::build(*mesh,metricsFor(*mesh,false),partial);},"missing selected-cell phase geometry");
    Case historical=s;historical.pureFirstJet=false;
    reject([&]{HydroPhaseGeometry::build(*mesh,metricsFor(*mesh,false),initial(*mesh,historical));},"requires FirstJet");
    auto otherMesh=makeMesh(4);
    // Actual identity is geometric, so an equivalent mesh is permitted; a
    // changed coordinate or topology of the same dimensions must reject.
    auto& altered=otherMesh->verts(otherMesh->vertIndexToLabel[0]);altered.setCoord(altered.coord()+Vector2D{.125,0});
    reject([&]{phaseRate(*otherMesh,geometry,values,geometry.geometryIdentity(),geometry.timeLayer());},"NativeMeshIdentity");
    DictionaryMesh::RefineBox box;box.min={.25,.25};box.max={.75,.75};box.level=1;
    mesh->applyRefineBoxes({box});
    reject([&]{evaluate(geometry.geometryIdentity(),geometry.timeLayer());},"NativeMeshIdentity");
    const auto rebuilt=HydroPhaseGeometry::build(*mesh,metricsFor(*mesh,false),publishAt(*mesh,s,{}));
    inspect(*mesh,s,rebuilt,true);
    reject([&]{HydroNativeValues(rebuilt.layout(),{0},{});},"dimensions");
    auto invalid=velocityFor(*mesh,rebuilt.layout());invalid[0]=std::numeric_limits<double>::infinity();
    reject([&]{HydroNativeValues(rebuilt.layout(),invalid,std::vector<double>(rebuilt.layout().pressureCount()));},"nonfinite");
    Case rz=s;rz.radial=true;
    const auto radial=HydroPhaseGeometry::build(*mesh,metricsFor(*mesh,true),publishAt(*mesh,rz,{}));
    auto radialU=velocityFor(*mesh,radial.layout());radialU[radial.metrics().axisVelocityDofs().front()]=1e-300;
    reject([&]{phaseRate(*mesh,radial,HydroNativeValues(radial.layout(),radialU,
        std::vector<double>(radial.layout().pressureCount())),radial.geometryIdentity(),radial.timeLayer());},"nonzero radial axis");
    radialU=velocityFor(*mesh,radial.layout());
    for(int i=0;i<radial.layout().velocityCount();++i)if(radial.metrics().faceMeasures()[i]>1)radialU[i]=std::numeric_limits<double>::max();
    reject([&]{phaseRate(*mesh,radial,HydroNativeValues(radial.layout(),radialU,
        std::vector<double>(radial.layout().pressureCount())),radial.geometryIdentity(),radial.timeLayer());},"nonfinite");
    // Results retain the publication/metrics after caller handles leave scope.
    const auto retained=[&]{const auto p=publishAt(*mesh,s,{});return HydroPhaseGeometry::build(*mesh,metricsFor(*mesh,false),p);}();
    inspect(*mesh,s,retained,true);
    for(auto t : {HermiteGalsTimeLayer{7,.2,.2,Position::Initial},HermiteGalsTimeLayer{8,.2,.4,Position::Midpoint},
        HermiteGalsTimeLayer{8,.2,.4,Position::End}}) {
        const auto timed=HydroPhaseGeometry::build(*mesh,metricsFor(*mesh,false),publishAt(*mesh,s,t));
        inspect(*mesh,s,timed,true);
        timed.validateFor(*mesh,timed.geometryIdentity(),t);
        // This same instant has a different valid midpoint interval and step.
        auto wrong=t;wrong.step+=1;
        reject([&]{timed.validateFor(*mesh,timed.geometryIdentity(),wrong);},"different full time layer");
    }
    const auto copied=retained;
    check(&copied.publication().state().geometry().snapshot()==&retained.publication().state().geometry().snapshot(),
        "copy duplicated unique scalar");
    auto callerU=velocityFor(*mesh,copied.layout());
    const HydroNativeValues copiedValues(copied.layout(),callerU,std::vector<double>(copied.layout().pressureCount()));
    const auto savedRate=phaseRate(*mesh,copied,copiedValues,copied.geometryIdentity(),copied.timeLayer());
    const auto saved=savedRate.faceRates();callerU.assign(callerU.size(),99);
    const auto repeat=phaseRate(*mesh,copied,copiedValues,copied.geometryIdentity(),copied.timeLayer());
    for(std::size_t i=0;i<saved.size();++i) {
        close(savedRate.faceRates()[i].liquid,saved[i].liquid,0,"immutable rate");
        close(repeat.faceRates()[i].liquid,saved[i].liquid,0,"velocity caller alias");
    }
}
}
HermiteGalsPublication tracePublication(const QuadtreeMesh& mesh,bool radial,int axis,double line,int sign,
    HermiteGalsTimeLayer time={},bool partial=false,int shape=0)
{
    const Profile f=[=](Vector2D p){
        if(shape==1){const double t=std::max(0.,p.y-.5);return Jet{sign*(p.x-.5+t*t),static_cast<double>(sign),sign*2*t,0};}
        if(shape==2)return Jet{sign*(p.x-.25)*(p.x-.75),sign*(2*p.x-1),0,0};
        const double n=axis==0?(radial?p.x*p.x-line*line:p.x-line):p.y-line;
        const double t=axis==0?p.y:p.x,a=radial && axis==1?1+t*t:1+t+t*t,da=radial && axis==1?2*t:1+2*t;
        const double dn=axis==0 && radial?2*p.x:1;
        return Jet{sign*n*a,sign*(axis==0?dn*a:n*da),sign*(axis==1?a:n*da),sign*dn*da};
    };
    const auto seed=stateFor(mesh,f,radial,time,[](int,Vector2D,Vector2D){return 0.;},{},HermiteGalsTransport::FirstJet);
    auto p=seed.primary();
    if(partial){p.geometryCells.erase(std::remove_if(p.geometryCells.begin(),p.geometryCells.end(),[&](int c){return mesh.cells(mesh.cellIndexToLabel[c]).centerCoord().x>=.5;}),p.geometryCells.end());
        std::set<Key> edges;std::set<int> vertices;
        for(int c:p.geometryCells){const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);for(int k=0;k<4;++k){
            Key key{mesh.verts(cell.cellToVert[k]).index(),mesh.verts(cell.cellToVert[(k+1)%4]).index()};const auto a=point(mesh,key[0]),b=point(mesh,key[1]);
            if(a.y==b.y?a.x>b.x:a.y>b.y)std::swap(key[0],key[1]);edges.insert(key);vertices.insert(key[0]);vertices.insert(key[1]);}}
        p.edges.erase(std::remove_if(p.edges.begin(),p.edges.end(),[&](const auto& e){return !edges.count(e.vertices);}),p.edges.end());
        p.vertices.erase(std::remove_if(p.vertices.begin(),p.vertices.end(),[&](const auto& v){return !vertices.count(v.vertex);}),p.vertices.end());}
    HermiteGalsCandidate candidate(mesh,time,p);const auto id=candidate.identity();const auto mixed=closeHermiteGalsMixedAmr(mesh,candidate);
    std::vector<double> center;for(int c:p.geometryCells)center.push_back(f(mesh.cells(mesh.cellIndexToLabel[c]).centerCoord()).value);
    HermiteGalsPublicationConfig config;
    if(shape==1)config.volume.boundaryContacts={{.5,0},{.25,1}};
    else if(shape==2)config.volume.boundaryContacts={{.25,0},{.25,1},{.75,0},{.75,1}};
    else config.volume.boundaryContacts=axis==0?std::vector<Vector2D>{{line,0},{line,1}}:(radial?std::vector<Vector2D>{{1,line}}:std::vector<Vector2D>{{0,line},{1,line}});
    return publishHermiteGalsFirstJetGeometry(mesh,std::move(candidate),mixed,{id,center},config);
}
void nativePhaseTraces()
{
    int cases=0,faces=0,zeroFaces=0,rateGuards=0,scaleGuards=0,exterior=0,selected=0;long double peak=0;
    const std::vector<HermiteGalsTimeLayer> times{{},{7,.2,.3,Position::Midpoint},{7,.2,.3,Position::End}};
    const auto inspectTrace=[&](const QuadtreeMesh& mesh,const HydroPhaseGeometry& g,int axis,double line,int sign,bool radial){
        check(g.faceTraces().size()==g.faceMeasures().size(),"native trace typed inventory size");
        const auto& scalar=g.publication().state().geometry().snapshot();
        check(&g.publication().facePartitions().source()==&scalar && g.publication().state().primary().phaseVolumes.empty(),"native H07 source/inventory");
        std::vector<double> velocities(g.layout().velocityCount());const auto pressures=std::vector<double>(g.layout().pressureCount());
        const auto zero=phaseRate(mesh,g,HydroNativeValues(g.layout(),velocities,pressures),g.geometryIdentity(),g.timeLayer());
        for(const auto& rate:zero.faceRates())check(rate.liquid==0 && rate.gas==0,"zero normal velocity has nonzero phase rate");
        for(int u=0;u<g.layout().velocityCount();++u){const auto face=g.layout().domain().face(mesh,g.layout().nativeFace(u));const auto& t=g.faceTraces()[u];const auto m=g.faceMeasures()[u];
            const auto& edge=face.horizontal?mesh.hEdges(face.edgeLabel):mesh.vEdges(face.edgeLabel);auto a=mesh.verts(edge.edgeToVert[0]).coord(),b=mesh.verts(edge.edgeToVert[1]).coord();
            const double coordinate=axis==0?a.x:a.y,delta=axis==0?b.x-a.x:b.y-a.y;
            const bool coincident=delta==0 && coordinate==line;long double expected=0;
            check(t.negativeCell==face.negativeCell && t.positiveCell==face.positiveCell,"selected trace availability changed native owner");
            if(coincident){++zeroFaces;check(t.coincidentHalfFaces==2,"zero face lost actual two halves");
                const bool both=t.minusAvailable() && t.plusAvailable();
                check(t.kind==(both?HydroPhaseFaceTraceKind::InteriorCoincident:HydroPhaseFaceTraceKind::BoundaryCoincident),"wrong coincident face kind");
                if(t.minusAvailable())close(t.minus.liquid,sign>0?face.area:0,1e-12L*face.area,"independent minus spatial trace");
                if(t.plusAvailable())close(t.plus.liquid,sign<0?face.area:0,1e-12L*face.area,"independent plus spatial trace");
                expected=both?face.area*.5L:(t.minusAvailable()?(sign>0?face.area:0):(sign<0?face.area:0));
                for(double speed:{-1e-300,1e-300,-.4,.4}){
                    for(int v=0;v<g.layout().velocityCount();++v){const auto f=g.layout().domain().face(mesh,g.layout().nativeFace(v));velocities[v]=(f.horizontal?1:0)==axis?(radial && axis==0?speed*f.center.x:speed):0;}
                    const bool available=speed>0?t.minusAvailable():t.plusAvailable();
                    std::optional<HydroPhaseRate> causal;
                    try{causal=phaseRate(mesh,g,HydroNativeValues(g.layout(),velocities,pressures),g.geometryIdentity(),g.timeLayer());}
                    catch(const std::exception& e){const std::string reason=e.what();
                        if(std::abs(speed)<1e-299 && reason.find("unsupported arithmetic scale")!=std::string::npos)++scaleGuards;
                        else check(!available && reason.find("without a formal selected/exterior phase port")!=std::string::npos,"unexpected native spatial rejection: "+reason);}
                    if(causal){check(available,"accepted unavailable causal trace");
                        const auto total=face.area*(radial && axis==0?speed*face.center.x:speed);const bool liquid=sign*(speed>0?-1:1)<0;
                        close(causal->faceRates()[u].liquid,liquid?total:0,1e-12L*std::abs(total),"independent causal liquid");
                        close(causal->faceRates()[u].gas,liquid?0:total,1e-12L*std::abs(total),"independent causal gas");}++rateGuards;}
                velocities.assign(velocities.size(),0);
            }else{check(t.kind==HydroPhaseFaceTraceKind::Ordinary && t.coincidentHalfFaces==0,"ordinary face changed kind");
                if(delta==0)expected=sign*(coordinate-line)<0?face.area:0;
                else {auto lo=a,hi=b;if(sign>0){if(axis==0)hi.x=std::max(a.x,std::min(b.x,line));else hi.y=std::max(a.y,std::min(b.y,line));}
                    else{if(axis==0)lo.x=std::min(b.x,std::max(a.x,line));else lo.y=std::min(b.y,std::max(a.y,line));}
                    expected=(face.horizontal?static_cast<long double>(hi.x)-lo.x:static_cast<long double>(hi.y)-lo.y)*(radial?pi*(static_cast<long double>(hi.x)+lo.x):1.L);}
            }
            const auto error=face.area>0?std::abs(m.liquid-expected)/face.area:std::abs(m.liquid-expected);
            check(error<=1e-12L,"actual native H07 independent face reference");peak=std::max(peak,error);++faces;
        }
        auto wrong=g.timeLayer();++wrong.step;
        rejects([&]{phaseRate(mesh,g,HydroNativeValues(g.layout(),velocities,pressures),g.geometryIdentity(),wrong);},"different full time layer");
        check(g.publication().state().identity()==g.geometryIdentity(),"failed native rate changed publication");
    };
    for(bool radial:{false,true})for(bool amr:{false,true})for(int axis:{0,1})for(int sign:{-1,1})for(const auto& time:times){
        auto mesh=makeMesh(4,amr?1:0);const auto pub=tracePublication(*mesh,radial,axis,.5,sign,time);
        const auto g=HydroPhaseGeometry::build(*mesh,metricsFor(*mesh,radial),pub);inspectTrace(*mesh,g,axis,.5,sign,radial);++cases;
    }
    for(bool radial:{false,true})for(int sign:{-1,1})for(bool partial:{false,true}){
        auto mesh=makeMesh(4,1);const auto pub=tracePublication(*mesh,radial,0,.5,sign,{},partial);
        const auto g=HydroPhaseGeometry::build(*mesh,metricsFor(*mesh,radial,true),pub);inspectTrace(*mesh,g,0,.5,sign,radial);
        check(g.report().boundaryCoincidentFaces>0,"selected boundary trace not exercised");++selected;
    }
    for(bool radial:{false,true})for(double line:{0.,1.})for(int sign:{-1,1}){
        if(radial && line==0)continue;auto mesh=makeMesh(4);const auto pub=tracePublication(*mesh,radial,0,line,sign);
        const auto g=HydroPhaseGeometry::build(*mesh,metricsFor(*mesh,radial),pub);inspectTrace(*mesh,g,0,line,sign,radial);
        check(g.report().boundaryCoincidentFaces>0,"physical boundary trace not exercised");++exterior;
    }
    auto mesh=makeMesh(4);for(int shape:{1,2}){const auto pub=tracePublication(*mesh,false,0,.5,1,{},false,shape);
        const auto g=HydroPhaseGeometry::build(*mesh,metricsFor(*mesh,false),pub);check(g.report().coincidentFaces>0,"joined/multiple H07 trace missing");
        check(pub.topology().curves.size()==static_cast<std::size_t>(shape),"joined/multiple H07 components changed");}
    const auto pub=tracePublication(*mesh,false,0,.5,1);const auto id=pub.state().identity();const auto* scalar=&pub.state().geometry().snapshot();
    const int rounding=std::fegetround();std::fesetround(FE_DOWNWARD);bool fp=false;
    try{HydroPhaseGeometry::build(*mesh,metricsFor(*mesh,false),pub);}catch(const std::exception&){fp=true;}std::fesetround(rounding);
    check(fp && pub.state().identity()==id && &pub.state().geometry().snapshot()==scalar,"native H07 FP failure changed original source");
    std::optional<HydroPhaseGeometry> retained;const Snapshot* retainedScalar=nullptr;
    {const auto temporary=tracePublication(*mesh,false,0,.5,1);retainedScalar=&temporary.state().geometry().snapshot();retained=HydroPhaseGeometry::build(*mesh,metricsFor(*mesh,false),temporary);}
    check(&retained->publication().state().geometry().snapshot()==retainedScalar,"native H07 lifetime/source identity");
    inspectTrace(*mesh,*retained,0,.5,1,false);
    std::cout<<"NATIVE_PHASE_TRACE,"<<cases<<','<<selected<<','<<exterior<<','<<faces<<','<<zeroFaces<<','<<rateGuards<<','<<peak<<",PASS\n";
    std::cout<<"NATIVE_SPATIAL_SCALE_REJECTIONS,"<<scaleGuards<<"\n";
}
namespace
{
static_assert(std::is_same_v<decltype(std::declval<HydroPhaseRate&>().atomRates()),const std::vector<HydroPhaseMeasure>&>);
// Independent physical polynomial and primitive for the declared affine MMS.
// No mapped scalar ranges, advector coefficients or new rate enter references.
Vector2D spatialField(Vector2D p,bool rz,int mode=2)
{
    if(mode==0)return {};
    if(mode==1)return rz?Vector2D{-.13*p.x,-.35}:Vector2D{-.3,-.35};
    if(mode==3)return rz?Vector2D{0,.5*p.x*p.x-.015625}:Vector2D{p.y-.375,0};
    return rz?Vector2D{.13*p.x,.2+.155*p.x*p.x+.19*p.y}:Vector2D{.3+.17*p.x-.23*p.y,-.2+.31*p.x+.19*p.y};
}
std::vector<double> spatialValues(const QuadtreeMesh& mesh,const HydroDofLayout& l,int mode=2)
{
    const bool rz=l.domain().axisymmetric();std::vector<double> out;
    for(int u=0;u<l.velocityCount();++u){const auto f=l.domain().face(mesh,l.nativeFace(u));const auto v=spatialField(f.center,rz,mode);
        if(rz && f.horizontal && (mode==2 || mode==3)){
            const auto& e=mesh.hEdges(f.edgeLabel);const auto a=mesh.verts(e.edgeToVert[0]).coord(),b=mesh.verts(e.edgeToVert[1]).coord();
            out.push_back(mode==2?.2+.155*.5*(a.x*a.x+b.x*b.x)+.19*f.center.y:.25*(a.x*a.x+b.x*b.x)-.015625);
        }else out.push_back(f.horizontal?v.y:v.x);
    }return out;
}
long double physicalIntegral(Vector2D a,Vector2D b,int normal,bool rz,int mode)
{
    const long double lo=normal==0?a.y:a.x,hi=normal==0?b.y:b.x,w=hi-lo;
    const auto v=spatialField(a,rz,mode);long double c0=normal==0?v.x:v.y,c1=0,c2=0;
    if(mode==2){if(rz && normal==1){c1=.31L*a.x;c2=.155L;}else if(!rz)c1=normal==0?-.23L:.31L;}
    if(mode==3){if(rz && normal==1){c1=a.x;c2=.5L;}else if(!rz && normal==0)c1=1;}
    const auto primitive=c0*w+c1*w*w/2+c2*w*w*w/3;
    return rz?2*pi*(a.x*primitive+(normal==1?c0*w*w/2+c1*w*w*w/3+c2*w*w*w*w/4:0)):primitive;
}
HydroPhaseMeasure spatialReference(const HydroAdvectorAtom& atom,const Case& shape,int mode)
{
    const auto a=atom.begin,b=atom.end;const int n=atom.normalAxis,t=1-n;const bool rz=shape.radial;
    const auto pointAt=[&](long double x){return n==0?Vector2D{a.x,static_cast<double>(x)}:Vector2D{static_cast<double>(x),a.y};};
    const long double start=t==0?a.x:a.y,end=t==0?b.x:b.y;
    std::vector<long double> cuts{start,end};
    const bool zero=shape.planar && (shape.axial()?n==1 && a.y==shape.planeOffset:n==0 && a.x==shape.planeOffset);
    if(shape.planar){if((shape.axial()?t==1:t==0) && start<shape.planeOffset && shape.planeOffset<end)cuts.push_back(shape.planeOffset);}
    else {const long double centerN=n==0?(rz?0:shape.cx):shape.cy,centerT=t==0?(rz?0:shape.cx):shape.cy;
        const long double delta=(n==0?a.x:a.y)-centerN,q=static_cast<long double>(shape.radius)*shape.radius-delta*delta;
        if(q>0)for(auto root:{centerT-std::sqrt(q),centerT+std::sqrt(q)})if(start<root && root<end)cuts.push_back(root);}
    // Only a coincident phase needs the independent positive/negative split.
    if(zero && mode==3){const long double root=rz?std::sqrt(.03125L):.375L;if(start<root && root<end)cuts.push_back(root);}
    std::sort(cuts.begin(),cuts.end());HydroPhaseMeasure out;
    for(std::size_t k=1;k<cuts.size();++k){const auto mid=pointAt((cuts[k-1]+cuts[k])/2);bool liquid=false;
        if(zero){const auto v=spatialField(mid,rz,mode);liquid=shape.exact(mid-v*1e-7,0).v<0;}
        else liquid=shape.exact(mid,0).v<0;
        const auto value=physicalIntegral(pointAt(cuts[k-1]),pointAt(cuts[k]),n,rz,mode);(liquid?out.liquid:out.gas)+=value;
    }return out;
}
long double spatialError(long double a,long double b,long double scale,const std::string& why)
{
    const auto e=std::abs(a-b)/(scale>0?scale:1.L);check(std::isfinite(e) && e<=1e-12L,why+" error="+std::to_string(static_cast<double>(e)));return e;
}
void inspectSpatial(const QuadtreeMesh& mesh,const Case& shape,const HydroPhaseRate& rate,int mode,
    std::size_t& rows,long double& peak)
{
    check(rate.isSpatial() && rate.advectorIdentity()==rate.spatialAdvector().identity(),"spatial rate actual provenance");
    check(rate.momentumGeometry().sameIdentity(rate.spatialAdvector().momentumGeometry()) && rate.geometry().publication().state().primary().phaseVolumes.empty(),"spatial rate geometry/inventory");
    const auto& r=rate.spatialAdvector().reconstruction();check(rate.spatialReport().sourceQueries==r.atoms().size(),"canonical scalar queried more than once");
    for(std::size_t i=0;i<r.atoms().size();++i){const auto expected=spatialReference(r.atoms()[i],shape,mode),actual=rate.atomRates()[i];const auto scale=r.atoms()[i].area;
        peak=std::max({peak,spatialError(actual.liquid,expected.liquid,scale,"spatial independent liquid"),spatialError(actual.gas,expected.gas,scale,"spatial independent gas")});++rows;}
    for(int u=0;u<rate.geometry().layout().velocityCount();++u){HydroPhaseMeasure expected;for(int i:r.primalAtoms()[u]){expected.liquid+=rate.atomRates()[i].liquid;expected.gas+=rate.atomRates()[i].gas;}
        spatialError(rate.faceRates()[u].liquid,expected.liquid,rate.geometry().metrics().faceMeasures()[u],"actual primal addback");
        spatialError(rate.faceRates()[u].gas,expected.gas,rate.geometry().metrics().faceMeasures()[u],"actual primal gas addback");
        spatialError(expected.liquid+expected.gas,rate.totalFaceRates()[u],rate.geometry().metrics().faceMeasures()[u],"original A*u total");}
    for(std::size_t g=0;g<rate.momentumGeometry().segments().size();++g){const auto& segment=rate.momentumGeometry().segments()[g];HydroPhaseMeasure expected;
        for(int i:r.dualAtoms()[g]){expected.liquid+=rate.atomRates()[i].liquid;expected.gas+=rate.atomRates()[i].gas;}
        spatialError(rate.dualFaceRates()[g].liquid,expected.liquid,segment.area,"actual dual liquid addback");spatialError(rate.dualFaceRates()[g].gas,expected.gas,segment.area,"actual dual gas addback");
    }
    std::vector<int> incidents(rate.dualFaceRates().size()),balance(incidents.size());
    for(std::size_t u=0;u<rate.momentumGeometry().volumes().size();++u)for(int side=0;side<4;++side)
        for(int g:rate.momentumGeometry().volumes()[u].sides[side]){
            const auto& segment=rate.momentumGeometry().segments()[g];const int sign=side%2?1:-1;
            check((sign>0?segment.negativeVelocity:segment.positiveVelocity)==static_cast<int>(u),"actual dual incidence owner");
            ++incidents[g];balance[g]+=sign;
        }
    for(std::size_t g=0;g<incidents.size();++g){const auto& segment=rate.momentumGeometry().segments()[g];
        const bool internal=segment.negativeVelocity>=0 && segment.positiveVelocity>=0;
        check(incidents[g]==(internal?2:1) && (internal?balance[g]==0:std::abs(balance[g])==1),"actual shared dual incidence");}
}
void spatialPhaseContract()
{
    int scenes=0,zeroNet=0,tiny=0;std::size_t rows=0;long double peak=0;
    for(bool rz:{false,true})for(int amr:{0,1}){auto mesh=makeMesh(4,amr);const auto momentum=HydroMomentumGeometry::build(*mesh,metricsFor(*mesh,rz));const auto reconstruction=HydroAdvectorReconstruction::prepare(*mesh,momentum);
        for(int kind=0;kind<4;++kind){Case shape;shape.pureFirstJet=true;shape.radial=rz;shape.planar=kind!=0;shape.radius=rz?.34:.19;shape.planeOffset=kind==1?.43:kind==2?2:-1;
            const auto phase=HydroPhaseGeometry::build(*mesh,momentum.metrics(),publishAt(*mesh,shape,{}));const auto advector=reconstruction.evaluate(*mesh,HydroNativeValues(momentum.layout(),spatialValues(*mesh,momentum.layout()),std::vector<double>(momentum.layout().pressureCount())),phase.geometryIdentity(),phase.timeLayer());
            const auto rate=HydroPhaseRate::prepareSpatial(*mesh,phase,advector);inspectSpatial(*mesh,shape,rate,2,rows,peak);++scenes;}
        for(bool native:{false,true})for(int sign:{-1,1}){const double line=native?.5:(amr?.125:.375);Case shape;shape.pureFirstJet=true;shape.radial=rz;shape.planar=true;shape.planeOffset=line;shape.planeSign=sign;
            const auto pub=tracePublication(*mesh,rz,rz?1:0,line,sign);const auto phase=HydroPhaseGeometry::build(*mesh,momentum.metrics(),pub);
            for(int mode=0;mode<4;++mode){const auto advector=reconstruction.evaluate(*mesh,HydroNativeValues(momentum.layout(),spatialValues(*mesh,momentum.layout(),mode),std::vector<double>(momentum.layout().pressureCount())),phase.geometryIdentity(),phase.timeLayer());const auto rate=HydroPhaseRate::prepareSpatial(*mesh,phase,advector);
                inspectSpatial(*mesh,shape,rate,mode,rows,peak);if(mode==3 && native){for(int u=0;u<momentum.layout().velocityCount();++u)if(rate.totalFaceRates()[u]==0 && rate.faceRates()[u].liquid!=0 && rate.faceRates()[u].gas!=0)++zeroNet;}++scenes;}
        }
        for(bool gas:{false,true}){Case shape;shape.pureFirstJet=true;shape.radial=rz;shape.planar=true;shape.planeOffset=gas?1-std::ldexp(1.,-40):std::ldexp(1.,-40);
            const auto phase=HydroPhaseGeometry::build(*mesh,momentum.metrics(),publishAt(*mesh,shape,{}));const auto advector=reconstruction.evaluate(*mesh,HydroNativeValues(momentum.layout(),spatialValues(*mesh,momentum.layout()),std::vector<double>(momentum.layout().pressureCount())),phase.geometryIdentity(),phase.timeLayer());const auto rate=HydroPhaseRate::prepareSpatial(*mesh,phase,advector);
            for(std::size_t i=0;i<rate.atomRates().size();++i){const auto expected=spatialReference(reconstruction.atoms()[i],shape,2);const auto actual=rate.atomRates()[i];const auto small=gas?expected.gas:expected.liquid;
                if(std::abs(small)>0 && std::abs(small)<1e-10L*reconstruction.atoms()[i].area){spatialError(gas?actual.gas:actual.liquid,small,std::abs(small),"tiny phase actual scale");++tiny;}}
            inspectSpatial(*mesh,shape,rate,2,rows,peak);++scenes;}
    }
    check(scenes==88 && rows>0 && zeroNet>0 && tiny>0,"spatial phase source/counterflow/tiny coverage");std::cout<<"SPATIAL_PHASE_SUMMARY,"<<scenes<<','<<rows<<','<<zeroNet<<','<<tiny<<','<<peak<<",PASS\n";
}
void spatialPhaseShapes()
{
    int scenes=0;std::size_t rows=0;long double peak=0;
    for(int amr:{0,1}){auto mesh=makeMesh(4,amr);const auto momentum=HydroMomentumGeometry::build(*mesh,metricsFor(*mesh,false));const auto reconstruction=HydroAdvectorReconstruction::prepare(*mesh,momentum);
        for(int kind=0;kind<4;++kind){const int sign=kind==1?-1:1;
            const Profile field=[=](Vector2D p){if(kind<2)return Jet{sign*(p.x+.25*p.y-.375),static_cast<double>(sign),sign*.25,0};
                if(kind==2){const double t=std::max(0.,p.y-.5);return Jet{p.x-.5+t*t,1,2*t,0};}
                return Jet{(p.x-.25)*(p.x-.75),2*p.x-1,0,0};};
            const auto pub=[&]{if(kind>=2)return tracePublication(*mesh,false,0,.5,1,{},false,kind-1);
                const auto seed=stateFor(*mesh,field,false,{},[](int,Vector2D,Vector2D){return 0.;},{},HermiteGalsTransport::FirstJet);
                HermiteGalsCandidate candidate(*mesh,{},seed.primary());const auto id=candidate.identity();const auto mixed=closeHermiteGalsMixedAmr(*mesh,candidate);
                std::vector<double> centers;for(int c:candidate.primary().geometryCells)centers.push_back(field(mesh->cells(mesh->cellIndexToLabel[c]).centerCoord()).value);
                HermiteGalsPublicationConfig config;config.volume.boundaryContacts={{.375,0},{.125,1}};
                return publishHermiteGalsFirstJetGeometry(*mesh,std::move(candidate),mixed,{id,centers},config);
            }();
            const auto phase=HydroPhaseGeometry::build(*mesh,momentum.metrics(),pub);const auto advector=reconstruction.evaluate(*mesh,HydroNativeValues(momentum.layout(),spatialValues(*mesh,momentum.layout()),std::vector<double>(momentum.layout().pressureCount())),phase.geometryIdentity(),phase.timeLayer());
            const auto rate=HydroPhaseRate::prepareSpatial(*mesh,phase,advector);
            for(std::size_t i=0;i<reconstruction.atoms().size();++i){const auto& atom=reconstruction.atoms()[i];const bool vertical=atom.normalAxis==0;
                const double lo=vertical?atom.begin.y:atom.begin.x,hi=vertical?atom.end.y:atom.end.x;
                const auto at=[&](double t){return vertical?Vector2D{atom.begin.x,t}:Vector2D{t,atom.begin.y};};
                std::vector<double> cuts{lo,hi};const auto cut=[&](double t){if(lo<t && t<hi)cuts.push_back(t);};
                if(kind<2)cut(vertical?(.375-atom.begin.x)/.25:.375-.25*atom.begin.y);
                else if(kind==2){if(vertical){cut(.5);if(atom.begin.x<.5)cut(.5+std::sqrt(.5-atom.begin.x));}else{const double t=std::max(0.,atom.begin.y-.5);cut(.5-t*t);}}
                else if(!vertical){cut(.25);cut(.75);}
                std::sort(cuts.begin(),cuts.end());HydroPhaseMeasure expected;
                for(std::size_t k=1;k<cuts.size();++k){const auto mid=at((cuts[k-1]+cuts[k])/2);const auto value=field(mid).value;
                    const bool liquid=(value==0?field(mid-spatialField(mid,false)*1e-7).value:value)<0;
                    (liquid?expected.liquid:expected.gas)+=physicalIntegral(at(cuts[k-1]),at(cuts[k]),atom.normalAxis,false,2);}
                const auto actual=rate.atomRates()[i];peak=std::max({peak,spatialError(actual.liquid,expected.liquid,atom.area,"oblique/joined/multiple liquid"),spatialError(actual.gas,expected.gas,atom.area,"oblique/joined/multiple gas")});++rows;
            }++scenes;
        }
    }
    check(scenes==8 && rows>0,"oblique/joined/multiple coverage");std::cout<<"SPATIAL_PHASE_SHAPES,"<<scenes<<','<<rows<<','<<peak<<",PASS\n";
}
void spatialPhaseGuards()
{
    for(bool rz:{false,true})for(bool partial:{false,true}){auto mesh=makeMesh(4,1);const auto momentum=HydroMomentumGeometry::build(*mesh,metricsFor(*mesh,rz,partial));const auto reconstruction=HydroAdvectorReconstruction::prepare(*mesh,momentum);
        const auto pub=tracePublication(*mesh,rz,0,partial?.5:1,1);const auto phase=HydroPhaseGeometry::build(*mesh,momentum.metrics(),pub);
        const auto evaluate=[&](double speed){std::vector<double> values(momentum.layout().velocityCount());for(int u=0;u<momentum.layout().velocityCount();++u){const auto f=momentum.layout().domain().face(*mesh,momentum.layout().nativeFace(u));values[u]=f.horizontal?0:(rz?speed*f.center.x:speed);}return reconstruction.evaluate(*mesh,HydroNativeValues(momentum.layout(),values,std::vector<double>(momentum.layout().pressureCount())),phase.geometryIdentity(),phase.timeLayer());};
        const auto outgoing=HydroPhaseRate::prepareSpatial(*mesh,phase,evaluate(1));const auto saved=outgoing.atomRates();rejects([&]{HydroPhaseRate::prepareSpatial(*mesh,phase,evaluate(-1));},"without a formal selected/exterior phase port");
        check(outgoing.atomRates().size()==saved.size() && outgoing.atomRates()[0].liquid==saved[0].liquid,"failed factory changed original rate");
        const auto stationary=HydroPhaseRate::prepareSpatial(*mesh,phase,evaluate(0));for(auto v:stationary.atomRates())check(v.liquid==0 && v.gas==0,"stationary spatial zero trace rate");
    }
    auto mesh=makeMesh(4,1);Case shape;shape.pureFirstJet=true;shape.planar=true;shape.planeOffset=.43;
    const auto geometry=HydroMomentumGeometry::build(*mesh,metricsFor(*mesh,false));const auto reconstruction=HydroAdvectorReconstruction::prepare(*mesh,geometry);const auto phase=HydroPhaseGeometry::build(*mesh,geometry.metrics(),publishAt(*mesh,shape,{}));
    const HydroNativeValues values(geometry.layout(),spatialValues(*mesh,geometry.layout()),std::vector<double>(geometry.layout().pressureCount()));const auto advector=reconstruction.evaluate(*mesh,values,phase.geometryIdentity(),phase.timeLayer());
    const auto rate=HydroPhaseRate::prepareSpatial(*mesh,phase,advector);const auto copy=rate;
    const auto foreign=HydroPhaseGeometry::build(*mesh,geometry.metrics(),publishAt(*mesh,shape,{}));rejects([&]{HydroPhaseRate::prepareSpatial(*mesh,foreign,advector);},"different geometry identity");
    const auto wrongTime=reconstruction.evaluate(*mesh,values,phase.geometryIdentity(),{7,.2,.3,Position::Midpoint});rejects([&]{HydroPhaseRate::prepareSpatial(*mesh,phase,wrongTime);},"different full time layer");
    const auto fresh=HydroMomentumGeometry::build(*mesh,metricsFor(*mesh,false));const auto freshAdvector=HydroAdvectorReconstruction::prepare(*mesh,fresh).evaluate(*mesh,HydroNativeValues(fresh.layout(),values.normalVelocity(),std::vector<double>(fresh.layout().pressureCount())),phase.geometryIdentity(),phase.timeLayer());rejects([&]{HydroPhaseRate::prepareSpatial(*mesh,phase,freshAdvector);},"different spatial native layout");
    const auto otherDi=HydroPhysicalDivergence::build(*mesh,HydroIntegralDivergence::build(*mesh,geometry.layout()));const auto otherPhase=HydroPhaseGeometry::build(*mesh,otherDi,phase.publication());rejects([&]{HydroPhaseRate::prepareSpatial(*mesh,otherPhase,advector);},"different original D_I");
    check(rate.isSpatial() && !rate.atomRates().empty(),"migrated phase owner lost atom port");
    const auto rounding=std::fegetround();std::fesetround(FE_DOWNWARD);bool rejected=false;try{HydroPhaseRate::prepareSpatial(*mesh,phase,advector);}catch(const std::exception&){rejected=true;}std::fesetround(rounding);check(rejected,"spatial factory accepted wrong FP");
    const auto retained=[&]{const auto g=HydroMomentumGeometry::build(*mesh,metricsFor(*mesh,false));const auto p=HydroPhaseGeometry::build(*mesh,g.metrics(),publishAt(*mesh,shape,{}));const auto r=HydroAdvectorReconstruction::prepare(*mesh,g);return HydroPhaseRate::prepareSpatial(*mesh,p,r.evaluate(*mesh,HydroNativeValues(g.layout(),spatialValues(*mesh,g.layout()),std::vector<double>(g.layout().pressureCount())),p.geometryIdentity(),p.timeLayer()));}();
    check(retained.advectorIdentity().value()!=0 && copy.advectorIdentity()==rate.advectorIdentity() && &copy.atomRates()==&rate.atomRates(),"spatial copy/source lifetime");
    auto equivalent=makeMesh(4,1);mesh.reset();retained.geometry().validateFor(*equivalent,retained.geometry().geometryIdentity(),retained.geometry().timeLayer());
    auto& vertex=equivalent->verts(equivalent->vertIndexToLabel[0]);const auto old=vertex.coord();vertex.setCoord(old+Vector2D{.01,0});rejects([&]{HydroPhaseRate::prepareSpatial(*equivalent,phase,advector);},"NativeMeshIdentity");vertex.setCoord(old);
    DictionaryMesh::RefineBox box;box.level=2;box.min={.25,.25};box.max={.75,.75};equivalent->applyRefineBoxes({box});rejects([&]{HydroPhaseRate::prepareSpatial(*equivalent,phase,advector);},"NativeMeshIdentity");
    std::cout<<"SPATIAL_PHASE_GUARDS,PASS\n";
}
}
namespace {
// Independent GL3 integration of the public physical polynomial. The work
// bound uses absolute coefficients, never the possibly-zero signed mean.
std::array<long double,2> advectorIntegral(const HydroAdvectorAtom& a,const HydroNormalPolynomial& p)
{
    const auto w=p.length;const long double g=std::sqrt(3.L/5),nodes[3]={-g,0,g},weights[3]={5.L/9,8.L/9,5.L/9};long double q=0;
    for(int j=0;j<3;++j){const auto t=w*(1+nodes[j])/2;const auto& c=p.beginCoefficients;long double v=0;for(std::size_t k=0;k<c.size();++k){auto term=c[k];for(std::size_t l=0;l<k;++l)term*=t;v+=term;}
        q+=weights[j]*v*(p.radialMeasure?2*pi*(a.begin.x+(a.normalAxis==1?t:0)):1)*w/2;}
    const auto& c=p.beginCoefficients;long double work=0;for(std::size_t k=0;k<c.size();++k){auto term=std::abs(c[k]);for(std::size_t l=0;l<k;++l)term*=w;work+=term;}return {q,a.area*work};
}
// Physical divergence density: polynomial through degree five on each actual
// piece for both supported sources. GL3 is exact; the work scale retains each
// signed derivative contribution before cancellation. No volume-coordinate fit.
using PhysicalCuts=std::array<std::vector<double>,2>;
std::vector<PhysicalCuts> physicalPartitions(const HydroAdvectorReconstruction& r)
{
    std::vector<PhysicalCuts> out(r.momentumGeometry().layout().pressureCount());
    for(const auto& p:r.pieces())for(int d=0;d<2;++d){out[p.pressure][d].push_back(d?p.lower.y:p.lower.x);out[p.pressure][d].push_back(d?p.upper.y:p.upper.x);}
    for(auto& cell:out)for(auto& line:cell){std::sort(line.begin(),line.end());line.erase(std::unique(line.begin(),line.end()),line.end());check(line.size()>=2,"missing physical partition");}
    return out;
}
template<class Query> std::array<long double,2> physicalVolumeIntegral(bool rz,int p,Vector2D lo,Vector2D hi,const PhysicalCuts& partition,Query query)
{
    PhysicalCuts cuts;
    for(int d=0;d<2;++d){const auto a=d?lo.y:lo.x,b=d?hi.y:hi.x;check(b>a,"empty physical integration region");cuts[d].push_back(a);for(auto k:partition[d])if(a<k&&k<b)cuts[d].push_back(k);cuts[d].push_back(b);}
    const long double g=std::sqrt(.6L),nodes[3]={-g,0,g},weights[3]={5.L/9,8.L/9,5.L/9};long double value=0,work=0;
    for(std::size_t x=1;x<cuts[0].size();++x)for(std::size_t y=1;y<cuts[1].size();++y){const long double hx=cuts[0][x]-cuts[0][x-1],hy=cuts[1][y]-cuts[1][y-1];
        for(int i=0;i<3;++i)for(int j=0;j<3;++j){const Vector2D point{static_cast<double>(cuts[0][x-1]+hx*(1+nodes[i])/2),static_cast<double>(cuts[1][y-1]+hy*(1+nodes[j])/2)};
            const auto a=query(p,point).jet;const long double xx=a.derivative.xx,yy=a.derivative.yy;
            const long double terms[3]={rz?point.x*xx:xx,rz?point.x*yy:yy,rz?a.value.x:0.L};const auto w=weights[i]*weights[j]*hx*hy/4*(rz?2*pi:1);
            for(auto term:terms){value+=w*term;work+=w*std::abs(term);}}
    }return {value,work};
}

// All terms of a physical-distance chart are compared in velocity units.
using Coefficients=std::array<long double,5>;
long double chartNorm(const Coefficients& a,long double length,int derivative=0)
{
    long double result=0,power=1;
    for(int k=0;k<5;++k){result+=(derivative?k:1)*std::abs(a[k])*power;power*=length;}
    return result;
}
long double budgetError(long double residual,long double work,const char* why)
{
    check(std::isfinite(residual)&&std::isfinite(work)&&work>=0,"invalid physical error budget");
    const auto error=work>0?std::abs(residual)/work:(residual==0?0:std::numeric_limits<long double>::infinity());
    check(error<=1e-12L,std::string(why)+" normalized="+std::to_string(static_cast<double>(error)));return error;
}
Coefficients chartDifference(const Coefficients& sum,const Coefficients& a,const Coefficients& b)
{Coefficients d{};for(int k=0;k<5;++k)d[k]=sum[k]-a[k]-b[k];return d;}

void advectorClosureContract()
{
    int bases=0,updates=0,zeroCells=0,faceRows=0,pointRows=0;long double peakMean=0,peakDiv=0,peakPoint=0,peakLinear=0;
    for(bool rz:{false,true})for(int amr:{0,1}){auto mesh=makeMesh(4,amr);Case shape;shape.pureFirstJet=true;shape.radial=rz;shape.planar=true;shape.planeOffset=2;
        const auto g=HydroMomentumGeometry::build(*mesh,metricsFor(*mesh,rz));const auto r=HydroAdvectorReconstruction::prepare(*mesh,g);const auto& l=g.layout();const auto pub=publishAt(*mesh,shape,{});const auto partitions=physicalPartitions(r);
        std::vector<std::optional<HydroSpatialAdvector>> unitFields(l.velocityCount());
        std::vector<std::vector<HydroNormalPolynomial>> unitResponses(l.velocityCount());std::vector<long double> columnNorm(l.velocityCount());
        for(int u=0;u<l.velocityCount();++u){const auto area=g.metrics().faceMeasures()[u];if(area==0)continue;std::vector<double> v(l.velocityCount());v[u]=1/area;
            unitFields[u]=r.evaluate(*mesh,HydroNativeValues(l,v,std::vector<double>(l.pressureCount())),pub.state().identity(),pub.state().timeLayer());unitResponses[u]=unitFields[u]->atomPolynomials();
            for(const auto& p:unitResponses[u])for(int end=0;end<2;++end)columnNorm[u]=std::max(columnNorm[u],area*chartNorm(end?p.endCoefficients:p.beginCoefficients,p.length));}
        const auto inspectField=[&](const HydroSpatialAdvector& b){const auto& values=b.values().normalVelocity();long double fieldWork=0;
            for(int u=0;u<l.velocityCount();++u)fieldWork+=std::abs(values[u])*columnNorm[u];
            check(b.boundary()==nullptr,"ordinary evaluation inferred boundary law");
            std::vector<long double> native(l.pressureCount()),integrated(l.pressureCount()),work(l.pressureCount()),gross(l.pressureCount());
            for(int u=0;u<l.velocityCount();++u){long double q=0,w=0;for(int i:r.primalAtoms()[u]){const auto ref=advectorIntegral(r.atoms()[i],b.atomPolynomials()[i]);q+=ref[0];w+=ref[1];const auto& p=b.atomPolynomials()[i];
                    spatialError(p.integral(false,0,p.length),ref[0],ref[1],"independent shifted normal integral");spatialError(p.integral(true,0,p.length),ref[0],ref[1],"independent reverse normal integral");}
                const auto expected=static_cast<long double>(g.metrics().faceMeasures()[u])*values[u];peakMean=std::max(peakMean,budgetError(q-expected,w+std::abs(expected)+g.metrics().faceMeasures()[u]*fieldWork,"all-basis physical native mean"));++faceRows;
                const auto f=l.domain().face(*mesh,l.nativeFace(u));for(auto side:{std::pair<int,int>{f.negativeCell,1},{f.positiveCell,-1}})if(side.first>=0){const int p=l.pressureDof(side.first);native[p]+=side.second*expected;integrated[p]+=side.second*q;gross[p]+=std::abs(expected);work[p]+=w+g.metrics().faceMeasures()[u]*fieldWork;}}
            const auto original=g.metrics().apply(*mesh,b.values());for(int p=0;p<l.pressureCount();++p){const auto v=g.metrics().cellVolumes()[p];peakDiv=std::max(peakDiv,spatialError(b.cellDivergence()[p]*v,native[p],gross[p],"actual B divergence independent incidence"));
                spatialError(original[p]*v,native[p],gross[p],"original DI independent incidence");budgetError(integrated[p]-native[p],work[p],"actual trace divergence integral");if(gross[p]==0){check(b.cellDivergence()[p]==0,"zero native-gross divergence not exact");++zeroCells;}
                const auto& cell=mesh->cells(mesh->cellIndexToLabel[l.nativeCell(p)]);const auto lo=mesh->verts(cell.cellToVert[0]).coord(),hi=mesh->verts(cell.cellToVert[2]).coord();
                const auto integral=physicalVolumeIntegral(rz,p,lo,hi,partitions[p],[&](int cell,Vector2D point){return b.cellVelocityJet(cell,point);});
                peakPoint=std::max(peakPoint,budgetError(integral[0]-native[p],integral[1]+work[p]+gross[p],"physical piecewise divergence integral"));++pointRows;
            }++updates;return b;
        };
        const auto inspect=[&](const std::vector<double>& v){return inspectField(r.evaluate(*mesh,HydroNativeValues(l,v,std::vector<double>(l.pressureCount())),pub.state().identity(),pub.state().timeLayer()));};
        for(const auto& field:unitFields)if(field){inspectField(*field);++bases;}
        std::vector<double> zero(l.velocityCount());const auto allZero=inspect(zero);for(const auto& p:allZero.atomPolynomials())for(auto c:p.beginCoefficients)check(c==0,"all-zero actual field not exact");
        auto a=spatialValues(*mesh,l),b=spatialValues(*mesh,l,3);std::vector<double> sum(a.size());for(std::size_t i=0;i<a.size();++i)sum[i]=a[i]+b[i];
        const auto fa=inspect(a),fb=inspect(b),fs=inspect(sum);
        for(std::size_t i=0;i<r.atoms().size();++i)for(int end=0;end<2;++end){const auto& pa=fa.atomPolynomials()[i];const auto& pb=fb.atomPolynomials()[i];const auto& ps=fs.atomPolynomials()[i];long double work=0;
            for(int u=0;u<l.velocityCount();++u)if(!unitResponses[u].empty()){const auto& q=unitResponses[u][i];work+=(std::abs(a[u])+std::abs(b[u]))*g.metrics().faceMeasures()[u]*chartNorm(end?q.endCoefficients:q.beginCoefficients,pa.length);}
            const auto delta=chartDifference(end?ps.endCoefficients:ps.beginCoefficients,end?pa.endCoefficients:pa.beginCoefficients,end?pb.endCoefficients:pb.beginCoefficients);
            peakLinear=std::max(peakLinear,budgetError(chartNorm(delta,pa.length),work,"actual polynomial linearity"));}
        for(double factor:{1e-200,1e200}){auto scaled=sum;for(auto& v:scaled)v*=factor;inspect(scaled);}
        for(std::size_t i=0;i<sum.size();++i)sum[i]=a[i]+(-a[i]+std::ldexp(b[i],-40));inspect(sum);
    }
    check(bases==208 && zeroCells>0,"all basis and zero-gross coverage");std::cout<<"ADVECTOR_CLOSURE,"<<bases<<','<<updates<<','<<faceRows<<','<<pointRows<<','<<zeroCells<<','<<peakMean<<','<<peakDiv<<','<<peakPoint<<','<<peakLinear<<",PASS\n";
}
void advectorBoundaryGuards()
{
    static_assert(!std::is_default_constructible_v<HydroAdvectorBoundary>);
    int guards=0;
    for(bool rz:{false,true})for(bool selected:{false,true}){auto mesh=makeMesh(4,1);const auto g=HydroMomentumGeometry::build(*mesh,metricsFor(*mesh,rz,selected));const auto& l=g.layout();const auto r=HydroAdvectorReconstruction::prepare(*mesh,g);Case shape;shape.pureFirstJet=true;shape.radial=rz;shape.planar=true;shape.planeOffset=2;const auto pub=publishAt(*mesh,shape,{});
        std::vector<int> walls;int inner=-1;for(int u=0;u<l.velocityCount();++u){const auto f=l.domain().face(*mesh,l.nativeFace(u));if(f.outwardSign!=0)walls.push_back(u);else inner=u;}
        const auto law=HydroAdvectorBoundary::impermeable(*mesh,g,walls),copy=law;check(copy.sameIdentity(law) && law.momentumGeometry().sameIdentity(g),"boundary copy/source identity");
        auto v=spatialValues(*mesh,l);for(int u:walls)v[u]=0;const HydroNativeValues values(l,v,std::vector<double>(l.pressureCount()));const auto b=r.evaluate(*mesh,values,pub.state().identity(),pub.state().timeLayer(),copy);
        check(b.boundary() && b.boundary()->sameIdentity(law) && b.values().normalVelocity()==v,"actual boundary binding changed values");
        const auto reject=[&](auto f,const char* why){rejects(f,why);++guards;};
        reject([&]{HydroAdvectorBoundary::impermeable(*mesh,g,{walls[0],walls[0]});},"duplicate impermeable face");reject([&]{HydroAdvectorBoundary::impermeable(*mesh,g,{-1});},"impermeable velocity index");reject([&]{HydroAdvectorBoundary::impermeable(*mesh,g,{l.velocityCount()});},"impermeable velocity index");reject([&]{HydroAdvectorBoundary::impermeable(*mesh,g,{inner});},"selected-domain boundary face");
        auto wrong=v;const int wall=*std::find_if(walls.begin(),walls.end(),[&](int u){return g.metrics().faceMeasures()[u]>0;});wrong[wall]=std::numeric_limits<double>::denorm_min();
        reject([&]{r.evaluate(*mesh,HydroNativeValues(l,wrong,std::vector<double>(l.pressureCount())),pub.state().identity(),pub.state().timeLayer(),law);},"nonzero impermeable native mean");
        const auto other=HydroMomentumGeometry::build(*mesh,g.metrics());const auto wrongLaw=HydroAdvectorBoundary::impermeable(*mesh,other,walls);reject([&]{r.evaluate(*mesh,values,pub.state().identity(),pub.state().timeLayer(),wrongLaw);},"different impermeable momentum geometry instance");
        reject([&]{b.validateFor(*mesh,g,{},pub.state().timeLayer());},"different velocity geometry identity");reject([&]{b.validateFor(*mesh,g,pub.state().identity(),{1,0,.1,Position::Midpoint});},"different full time layer");
        const auto retained=[&]{const auto temporary=HydroAdvectorBoundary::impermeable(*mesh,g,walls);return r.evaluate(*mesh,values,pub.state().identity(),pub.state().timeLayer(),temporary);}();
        check(retained.boundary()->velocityDofs()==walls,"boundary lifetime");auto equivalent=makeMesh(4,1);mesh.reset();retained.validateFor(*equivalent,g,pub.state().identity(),pub.state().timeLayer());
        DictionaryMesh::RefineBox box;box.level=2;box.min={.25,.25};box.max={.75,.75};equivalent->applyRefineBoxes({box});reject([&]{law.validateFor(*equivalent,g);},"NativeMeshIdentity");reject([&]{retained.validateFor(*equivalent,g,pub.state().identity(),pub.state().timeLayer());},"NativeMeshIdentity");
    }
    std::cout<<"ADVECTOR_BOUNDARY_GUARDS,"<<guards<<",PASS\n";
}
// The actual production conversion is tested through the existing exact
// range and positive-range consumers. GL3 and known factor roots are reference
// data only; no independently reconstructed source is published.
void completeNormalRateContract()
{
    using namespace LAMP::Geometry;
    const auto began=std::chrono::steady_clock::now();
    int ranges=0,positive=0,tiny=0;long double peak=0;
    const auto reference=[](const HydroNormalPolynomial& p,long double lo,long double hi) {
        const long double g=std::sqrt(3.L/5),nodes[3]={-g,0,g},weights[3]={5.L/9,8.L/9,5.L/9};long double out=0,work=0;
        for(int i=0;i<3;++i){const auto x=(lo+hi)/2+(hi-lo)*nodes[i]/2;long double value=0,scale=0;
            for(std::size_t k=0;k<p.beginCoefficients.size();++k){long double term=p.beginCoefficients[k];for(std::size_t j=0;j<k;++j)term*=x;value+=term;scale+=std::abs(term);}
            const auto measure=p.radialMeasure?CartesianMeasure::twoPi*(p.radiusBegin+(p.radialTangent?x:0)):1;
            out+=weights[i]*value*measure*(hi-lo)/2;work+=weights[i]*scale*measure*(hi-lo)/2;}
        return std::array<long double,2>{out,work};
    };
    std::cout<<"RATE_CONSTANT,"<<CartesianMeasure::twoPi<<'\n';
    for(int family=0;family<6;++family)for(long double length:{.5L,2.L})for(int metric=0;metric<4;++metric) {
        HydroNormalPolynomial p;p.length=length;p.radialMeasure=metric>0;p.radialTangent=metric>=2;p.radiusBegin=metric==3?0:.25L;
        if(family<5)p.beginCoefficients[family]=1;else p.beginCoefficients={1,-2,3,-4,5};
        const auto terms=LAMP::Hydrodynamics::Detail::normalRatePolynomial(p);
        if(family>=3){bool high=false;for(const auto& t:terms)for(int k=family==3?3:4;k<=5;++k)high=high || t.coefficients[k]*t.factor!=0;check(high,"actual rate conversion lost high degree");}
        if(metric>=2 && family==4){bool fifth=false;for(const auto& t:terms)fifth=fifth || t.coefficients[5]*t.factor!=0;check(fifth,"missing RZ fifth-degree density");}
        for(double s:{0.,.125,.5,.875,1.}){const long double x=length*s;long double value=0,work=0;
            for(std::size_t k=0;k<p.beginCoefficients.size();++k){const auto term=p.beginCoefficients[k]*std::pow(x,static_cast<int>(k));value+=term;work+=std::abs(term);}
            const auto measure=length*(metric?CartesianMeasure::twoPi*(p.radiusBegin+(p.radialTangent?x:0)):1);
            peak=std::max(peak,spatialError(evaluateUnitPolynomial(terms,s),value*measure,work*measure,"complete converted density value"));}
        for(bool end:{false,true})for(auto range:std::array<std::array<long double,2>,3>{{{0,1},{.125L,.375L},{.75L,.125L}}}) {
            const auto off=range[0],width=range[1],lo=length*(end?1-off-width:off),hi=lo+length*width;
            const auto ref=reference(p,lo,hi);const auto actual=integrateUnitPolynomialRange(terms,0,1,end,off,width);
            peak=std::max(peak,spatialError(actual,ref[0],ref[1],"complete rate physical GL3"));
            std::cout<<"RATE_RANGE,"<<family<<','<<length<<','<<metric<<','<<end<<','<<off<<','<<width<<','<<actual<<','<<ref[0]<<','<<ref[1]<<'\n';++ranges;
        }
    }
    const std::array<std::array<long double,4>,3> roots{{{.125L,.375L,.625L,.875L},{.25L,.25L,.75L,.75L},{0,.25L,.75L,1}}};
    for(int fixture=0;fixture<3;++fixture)for(int metric=0;metric<4;++metric)for(int sign:{-1,1}) {
        HydroNormalPolynomial p;p.length=1;p.radialMeasure=metric>0;p.radialTangent=metric>=2;p.radiusBegin=metric==3?0:.25L;
        p.beginCoefficients[0]=sign;int degree=0;
        for(auto root:roots[fixture]){auto c=p.beginCoefficients;p.beginCoefficients.fill(0);for(int k=0;k<=degree;++k){p.beginCoefficients[k]-=root*c[k];p.beginCoefficients[k+1]+=c[k];}++degree;}
        const auto terms=LAMP::Hydrodynamics::Detail::normalRatePolynomial(p);
        for(bool end:{false,true})for(auto range:std::array<std::array<long double,2>,2>{{{0,1},{.125L,.625L}}}) {
            const auto off=range[0],width=range[1],lo=end?1-off-width:off,hi=lo+width;
            std::vector<long double> cuts{lo,hi};for(auto root:roots[fixture])if(lo<root && root<hi)cuts.push_back(root);
            std::sort(cuts.begin(),cuts.end());cuts.erase(std::unique(cuts.begin(),cuts.end()),cuts.end());
            long double expected=0;for(std::size_t k=1;k<cuts.size();++k){long double value=sign;for(auto root:roots[fixture])value*=(cuts[k-1]+cuts[k])/2-root;
                if(value>0)expected+=reference(p,cuts[k-1],cuts[k])[0];}
            const auto work=reference(p,lo,hi)[1];const auto actual=integratePositiveUnitPolynomialRange(terms,0,1,end,off,width,1e-13L*work);
            check(actual.value>=0 && actual.absoluteError>=0 && std::isfinite(actual.absoluteError),"invalid complete positive rate");
            peak=std::max(peak,spatialError(actual.value,expected,work,"quartic causal rate from known roots"));
            std::cout<<"RATE_POSITIVE,"<<fixture<<','<<metric<<','<<sign<<','<<end<<','<<off<<','<<width<<','<<actual.value<<','<<expected<<','<<work<<','<<actual.absoluteError<<'\n';++positive;
        }
    }
    for(long double root:{0.L,.5L})for(int metric=0;metric<4;++metric)for(bool end:{false,true}) {
        if(root==0 && end)continue;
        HydroNormalPolynomial p;p.length=1;p.radialMeasure=metric>0;p.radialTangent=metric>=2;p.radiusBegin=metric==3?0:.25L;
        p.beginCoefficients={std::pow(root,4),-4*std::pow(root,3),6*root*root,-4*root,1};
        const auto terms=LAMP::Hydrodynamics::Detail::normalRatePolynomial(p);const auto offset=end?1-root:root;
        for(int exponent:{60,80}) {
            const auto width=std::ldexp(1.L,-exponent),base=std::pow(width,5)/5;
            const auto radius=p.radiusBegin+(p.radialTangent?root:0);
            const auto expected=metric?CartesianMeasure::twoPi*(radius*base+(p.radialTangent?(end?-1:1)*std::pow(width,6)/6:0)):base;
            const auto actual=integrateUnitPolynomialRange(terms,0,1,end,offset,width);
            const auto pos=integratePositiveUnitPolynomialRange(terms,0,1,end,offset,width,1e-13L*expected);
            check(actual>0 && pos.value>0,"tiny complete converted rate lost");
            peak=std::max({peak,spatialError(actual,expected,expected,"tiny complete rate"),spatialError(pos.value,expected,expected,"tiny complete positive rate")});
            std::cout<<"RATE_TINY,"<<root<<','<<metric<<','<<end<<','<<exponent<<','<<actual<<','<<pos.value<<','<<expected<<'\n';++tiny;
        }
    }
    check(ranges==288 && positive==96 && tiny==24,"complete rate coverage");
    std::cout<<"RATE_COMPLETE_SUMMARY,"<<ranges<<','<<positive<<','<<tiny<<','<<peak<<",PASS\n";
    std::cout<<"RATE_COMPLETE_SECONDS,"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-began).count()<<'\n';
}

}
namespace {
struct ArithmeticRegion {int atom,pressure,normal;Vector2D begin,end;};
Coefficients shiftedChart(const Coefficients& a,long double offset)
{
    constexpr int choose[5][5]={{1,0,0,0,0},{1,1,0,0,0},{1,2,1,0,0},{1,3,3,1,0},{1,4,6,4,1}};
    Coefficients out{};
    for(int k=0;k<5;++k)for(int j=k;j<5;++j)out[k]+=choose[j][k]*a[j]*std::pow(offset,j-k);
    return out;
}
void regularArithmeticContract()
{
    int layouts=0,fields=0,originalAtoms=0,rawFailed=0;std::size_t comparisons=0,unitChecks=0;
    // V4a1c current canonical atoms and physical chart coverage, in XY/RZ x uniform/AMR order.
    std::array<int,4> originalContextAtoms{},originalContextCharts{};
    std::array<long double,6> peak{};std::array<long double,4> mutation{};
    const long double node=std::sqrt(.6L),nodes[3]={-node,0,node},weights[3]={5.L/9,8.L/9,5.L/9};
    for(bool rz:{false,true})for(int amr:{0,1}){auto mesh=makeMesh(4,amr);Case shape;shape.pureFirstJet=true;shape.radial=rz;shape.planar=true;shape.planeOffset=2;const auto pub=publishAt(*mesh,shape,{});
        for(int selection=0;selection<4;++selection){FaceSegmentGraph graph;graph.build(*mesh,[&](const QuadtreeMesh&,const Cell& cell){const auto p=cell.centerCoord();return selection==0?true:selection==1?p.x<.5:selection==2?!(p.x>.25&&p.x<.75&&p.y>.25&&p.y<.75):p.x<.25||p.x>.75;});
            const auto layout=HydroDofLayout::build(*mesh,NativeLeafDomain::build(*mesh,graph,rz));
            const auto g=HydroMomentumGeometry::build(*mesh,HydroPhysicalDivergence::build(*mesh,HydroIntegralDivergence::build(*mesh,layout)));
            const auto unbound=HydroAdvectorReconstruction::prepare(*mesh,g);const auto& unboundFaces=unbound.regularGeometry().faces();
            std::vector<int> walls;for(std::size_t u=0;u<unboundFaces.size();++u)if(unboundFaces[u].external&&unboundFaces[u].area>0)walls.push_back(static_cast<int>(u));
            for(int boundary=0;boundary<3;++boundary){const auto chosen=boundary==2?std::vector<int>{walls.front()}:walls;const auto law=HydroAdvectorBoundary::impermeable(*mesh,g,chosen);
                const auto r=boundary?HydroAdvectorReconstruction::prepare(*mesh,g,law):unbound;const auto& regular=r.regularGeometry();const auto& faces=regular.faces();
                const auto coord=[](Vector2D p,int d){return d?p.y:p.x;};std::vector<ArithmeticRegion> regions;
                for(std::size_t i=0;i<r.atoms().size();++i){const auto& atom=r.atoms()[i];int p=atom.ownerPressure;
                    if(atom.primalVelocity>=0){const auto& f=faces[atom.primalVelocity];p=f.negativePressure>=0?f.negativePressure:f.positivePressure;}
                    const int t=1-atom.normalAxis;std::vector<double> cuts{coord(atom.begin,t)};for(double k:regular.cells()[p].supportKnots[t])if(k>cuts.front()&&k<coord(atom.end,t))cuts.push_back(k);cuts.push_back(coord(atom.end,t));
                    for(std::size_t k=1;k<cuts.size();++k){auto a=atom.begin,b=atom.end;(t?a.y:a.x)=cuts[k-1];(t?b.y:b.x)=cuts[k];regions.push_back({static_cast<int>(i),p,atom.normalAxis,a,b});}}
                std::vector<bool> excluded(faces.size());for(std::size_t u=0;u<faces.size();++u)excluded[u]=faces[u].area==0;
                if(boundary)for(int u:chosen)excluded[u]=true;
                const auto curves=[&](int source,const std::vector<double>& v){
                    const HydroNativeValues values(layout,v,std::vector<double>(layout.pressureCount()));std::vector<HydroNormalPolynomial> out;
                    if(source==0){const auto field=boundary?r.evaluate(*mesh,values,pub.state().identity(),pub.state().timeLayer(),law):r.evaluate(*mesh,values,pub.state().identity(),pub.state().timeLayer());
                        for(const auto& region:regions){auto poly=field.atomPolynomials()[region.atom];const auto& atom=r.atoms()[region.atom];const int t=1-region.normal;
                            poly.beginCoefficients=shiftedChart(poly.beginCoefficients,coord(region.begin,t)-coord(atom.begin,t));poly.endCoefficients=shiftedChart(poly.endCoefficients,coord(atom.end,t)-coord(region.end,t));
                            poly.length=coord(region.end,t)-coord(region.begin,t);poly.radiusBegin=region.begin.x;poly.radiusEnd=region.end.x;out.push_back(poly);}}
                    else{const auto data=regular.apply(*mesh,values,boundary?&law:nullptr);for(const auto& q:regions)out.push_back(regular.normalPolynomial(data,q.pressure,q.normal,q.begin,q.end));}
                    return out;
                };
                for(int source=0;source<2;++source){
                    std::vector<std::vector<HydroNormalPolynomial>> basis(faces.size(),std::vector<HydroNormalPolynomial>(regions.size()));
                    for(std::size_t u=0;u<faces.size();++u)if(!excluded[u]){std::vector<double> v(faces.size());v[u]=1;basis[u]=curves(source,v);
                    }
                    // Column norms of the WHOLE vector of physical charts.
                    // A zero functional can have cancelling intermediates;
                    // its final unit responses alone are not a roundoff scale.
                    std::vector<long double> columnNorm(faces.size());
                    for(std::size_t u=0;u<faces.size();++u)for(std::size_t i=0;i<regions.size();++i)for(int end=0;end<2;++end){const auto length=coord(regions[i].end,1-regions[i].normal)-coord(regions[i].begin,1-regions[i].normal);columnNorm[u]=std::max(columnNorm[u],chartNorm(end?basis[u][i].endCoefficients:basis[u][i].beginCoefficients,length));}
                    for(int family=0;family<3;++family)for(double magnitude:{0.,1e-200,1.,1e200}){
                        auto a=spatialValues(*mesh,layout),b=spatialValues(*mesh,layout,3);
                        for(std::size_t u=0;u<a.size();++u){if(family){a[u]=.17*std::sin(.71*u)+.11*std::cos(1.21*u);b[u]=family==1?.13*std::cos(.31*u):(-a[u]+std::ldexp(.07*std::sin(.93*u),-40));}
                            a[u]=excluded[u]?0:a[u]*magnitude;b[u]=excluded[u]?0:b[u]*magnitude;}
                        std::vector<double> sum(a.size()),low(a.size());
                        for(std::size_t u=0;u<a.size();++u){sum[u]=a[u]+b[u];const double virtualB=sum[u]-a[u];low[u]=(a[u]-(sum[u]-virtualB))+(b[u]-virtualB);}
                        const std::array<std::vector<HydroNormalPolynomial>,4> polys{curves(source,a),curves(source,b),curves(source,sum),curves(source,low)};
                        const HydroNativeValues sumValues(layout,sum,std::vector<double>(layout.pressureCount()));
                        const auto accepted=boundary?r.evaluate(*mesh,sumValues,pub.state().identity(),pub.state().timeLayer(),law):r.evaluate(*mesh,sumValues,pub.state().identity(),pub.state().timeLayer());const auto data=regular.apply(*mesh,sumValues,boundary?&law:nullptr);
                        const auto jet=[&](int p,Vector2D x){return source==0?accepted.cellVelocityJet(p,x).jet:regular.cellJet(data,p,x).jet;};
                        const bool original=family==0&&magnitude==1&&selection==0&&boundary==0;
                        if(original&&source==1)originalContextAtoms[2*rz+amr]+=static_cast<int>(r.atoms().size());
                        for(std::size_t i=0;i<regions.size();++i){const auto& q=regions[i];const auto& p=polys[2][i];const auto length=p.length;
                            if(original&&source==1){long double rawWork=0;for(int k=0;k<5;++k)rawWork+=std::abs(polys[0][i].beginCoefficients[k])+std::abs(polys[1][i].beginCoefficients[k]);
                                for(int k=0;k<5;++k)if(std::abs(p.beginCoefficients[k]-(polys[0][i].beginCoefficients[k]+polys[1][i].beginCoefficients[k]))>1e-12L*(rawWork>0?rawWork:1))++rawFailed;++originalAtoms;++originalContextCharts[2*rz+amr];}
                            for(int end=0;end<2;++end){const auto coeff=[&](int field)->const Coefficients& {return end?polys[field][i].endCoefficients:polys[field][i].beginCoefficients;};
                                Coefficients reference{};long double work=0,sumWork=0,nativeRounding=0;
                                for(std::size_t u=0;u<a.size();++u){const auto& response=end?basis[u][i].endCoefficients:basis[u][i].beginCoefficients;
                                    const auto w=columnNorm[u];work+=(std::abs(a[u])+std::abs(b[u]))*w;sumWork+=std::abs(sum[u])*w;nativeRounding+=std::abs(low[u])*w;
                                    for(int k=0;k<5;++k)reference[k]+=sum[u]*response[k];}
                                auto difference=chartDifference(coeff(2),coeff(0),coeff(1));for(int k=0;k<5;++k)difference[k]+=coeff(3)[k];
                                peak[0]=std::max(peak[0],budgetError(chartNorm(difference,length),work,"TwoSum field linearity"));
                                Coefficients residual{};for(int k=0;k<5;++k)residual[k]=coeff(2)[k]-reference[k];
                                peak[1]=std::max(peak[1],budgetError(chartNorm(residual,length),sumWork,"unit-response field reconstruction"));
                                if(work>0)peak[2]=std::max(peak[2],nativeRounding/work);
                                for(int exponent:{-8,8}){const long double unit=std::ldexp(1.L,exponent),speed=32;auto changed=p;changed.length*=unit;changed.radiusBegin*=unit;changed.radiusEnd*=unit;
                                    for(int k=0;k<5;++k){const auto factor=speed/std::pow(unit,k);changed.beginCoefficients[k]*=factor;changed.endCoefficients[k]*=factor;}
                                    const auto& c=end?changed.endCoefficients:changed.beginCoefficients;
                                    check(chartNorm(c,changed.length)==speed*chartNorm(coeff(2),length),"chart norm depends on length units");
                                    const auto metric=std::pow(unit,rz?2:1)*speed;
                                    peak[3]=std::max(peak[3],budgetError(changed.integral(end,changed.length*.21L,changed.length*.31L)/metric-p.integral(end,length*.21L,length*.31L),length*(rz?2*pi*std::max(q.begin.x,q.end.x):1)*sumWork,"integral units"));++unitChecks;}
                                ++comparisons;
                            }
                            long double physical=0,physicalWork=0;long double work=0;for(std::size_t u=0;u<sum.size();++u)work+=std::abs(sum[u])*columnNorm[u];
                            for(int j=0;j<3;++j){const long double z=(1+nodes[j])/2;const auto x=q.begin+(q.end-q.begin)*static_cast<double>(z);const auto actual=jet(q.pressure,x);const auto v=q.normal?actual.value.y:actual.value.x,d=q.normal?actual.derivative.yx:actual.derivative.xy;
                                long double derivative=0;for(int k=1;k<5;++k)derivative+=k*p.beginCoefficients[k]*std::pow(length*z,k-1);
                                peak[4]=std::max({peak[4],budgetError(p.value(false,length*z)-v,work,"physical value"),budgetError(p.value(true,length*(1-z))-v,work,"reverse physical value"),budgetError(length*(derivative-d),4*work,"physical tangent derivative")});
                                const auto weight=weights[j]*length*.5L*(rz?2*pi*x.x:1);physical+=weight*v;physicalWork+=weight*work;}
                            for(int end=0;end<2;++end)peak[5]=std::max(peak[5],budgetError(p.integral(end,0,length)-physical,physicalWork,"independent GL3 integral"));
                            if(family==1&&magnitude==1&&source==1){const auto denom=physicalWork;
                                if(denom>0){auto cut=p;cut.beginCoefficients[3]=cut.beginCoefficients[4]=0;
                                    mutation[0]=std::max(mutation[0],std::abs(-p.integral(false,0,length)-physical)/denom);
                                    mutation[1]=std::max(mutation[1],std::abs(cut.integral(false,0,length)-physical)/denom);
                                    mutation[2]=std::max(mutation[2],std::abs(p.integral(false,0,length*.5L)-physical)/denom);}
                            }
                        }
                        if(family==1&&magnitude==1&&source==1)for(std::size_t cell=0;cell<regular.cells().size();++cell){const auto& c=regular.cells()[cell];const auto query=[&](int p,Vector2D x){return regular.cellJet(data,p,x);};
                            const auto correct=physicalVolumeIntegral(rz,static_cast<int>(cell),c.lower,c.upper,c.supportKnots,query),uncut=physicalVolumeIntegral(rz,static_cast<int>(cell),c.lower,c.upper,PhysicalCuts{},query);
                            if(correct[1]>0)mutation[3]=std::max(mutation[3],std::abs(uncut[0]-correct[0])/correct[1]);}
                        ++fields;
                    }
                }++layouts;
            }
        }
    }
    constexpr std::array<int,4> expectedOriginalContexts{{288,248,144,248}};
    int originalRawAtoms=0;
    for(std::size_t context=0;context<expectedOriginalContexts.size();++context){
        check(originalContextAtoms[context]==expectedOriginalContexts[context]&&
              originalContextCharts[context]==expectedOriginalContexts[context]&&
              originalContextAtoms[context]==originalContextCharts[context],"arithmetic original geometry context coverage changed");
        originalRawAtoms+=originalContextAtoms[context];
    }
    check(layouts==48&&fields==1152&&originalAtoms==928&&originalRawAtoms==928,"arithmetic coverage or original witness changed");
    for(auto gap:mutation)check(gap>1e-8,"arithmetic negative control missed");
    for(std::size_t context=0;context<expectedOriginalContexts.size();++context)
        std::cout<<"ARITHMETIC_CONTEXT,"<<context<<','<<originalContextAtoms[context]<<','<<originalContextCharts[context]<<",PASS\n";
    std::cout<<"ARITHMETIC_BUDGET,"<<layouts<<','<<fields<<','<<comparisons<<','<<unitChecks<<','<<originalAtoms<<','<<rawFailed;for(auto p:peak)std::cout<<','<<p;std::cout<<",PASS\n";
    std::cout<<"ARITHMETIC_NEGATIVE";for(auto p:mutation)std::cout<<','<<p;std::cout<<",PASS\n";
}


// V4a: same retained VOF fan -> actual native dual material mass. The
// independent reference keeps twofold clipping intersections and integrates
// both phases, without calling the production Core clip/measure routines.
#if defined(_MSC_VER)
#pragma float_control(precise,on,push)
#endif
struct VofMaterialWide { long double high=0,low=0; };
VofMaterialWide vofMaterialSum(long double a,long double b)
{
    const auto s=a+b,p=s-a;return {s,(a-(s-p))+(b-p)};
}
VofMaterialWide vofMaterialAdd(VofMaterialWide a,VofMaterialWide b)
{
    const auto h=vofMaterialSum(a.high,b.high),l=vofMaterialSum(a.low,b.low);
    const auto m=vofMaterialSum(h.low,l.high),n=vofMaterialSum(h.high,m.high);
    return vofMaterialSum(n.high,n.low+m.low+l.low);
}
VofMaterialWide vofMaterialNeg(VofMaterialWide a){return {-a.high,-a.low};}
VofMaterialWide vofMaterialSub(VofMaterialWide a,VofMaterialWide b)
{return vofMaterialAdd(a,vofMaterialNeg(b));}
VofMaterialWide vofMaterialMul(VofMaterialWide a,VofMaterialWide b)
{
    const auto p=a.high*b.high;auto r=VofMaterialWide{p,std::fma(a.high,b.high,-p)};
    r=vofMaterialAdd(r,{a.high*b.low,0});r=vofMaterialAdd(r,{a.low*b.high,0});
    return vofMaterialAdd(r,{a.low*b.low,0});
}
VofMaterialWide vofMaterialDiv(VofMaterialWide a,VofMaterialWide b)
{
    check(b.high!=0,"VOF reference division by zero");const auto q=a.high/b.high;
    const auto r=vofMaterialSub(a,vofMaterialMul({q,0},b));
    const auto q1=(r.high+r.low)/b.high;
    const auto r1=vofMaterialSub(r,vofMaterialMul({q1,0},b));
    return vofMaterialAdd(vofMaterialAdd({q,0},{q1,0}),{(r1.high+r1.low)/b.high,0});
}
long double vofMaterialValue(VofMaterialWide a){return a.high+a.low;}
int vofMaterialSign(VofMaterialWide a)
{return a.high>0?1:a.high<0?-1:a.low>0?1:a.low<0?-1:0;}
using VofMaterialPoint=std::array<VofMaterialWide,2>;
using VofMaterialPolygon=std::vector<VofMaterialPoint>;
using VofMaterialBounds=std::array<VofMaterialWide,4>; // x0,x1,y0,y1
using VofPublication=LAMP::Interface::GeometricVofPublication;
using VofClock=LAMP::Interface::GeometricVofTimeLayer;
using VofTransport=LAMP::Interface::GeometricVofTransport;
unsigned vofMaterialLayers=0,vofMaterialComparisons=0,vofMaterialGuards=0,vofMaterialSlabs=0;
std::array<unsigned,2> vofMaterialAmrCuts{};
long double vofMaterialPeakPhase=0,vofMaterialPeakMass=0,vofMaterialPeakPartition=0,
    vofMaterialPeakInventory=0,vofMaterialPeakGasInventory=0,vofMaterialProjectionGap=0;
bool vofMaterialDifferentCandidates=false;
unsigned vofMaterialFaceChecks=0,vofMaterialCoincidentFaces=0,vofMaterialAxisFaces=0;
void vofMaterialNear(long double actual,VofMaterialWide expected,const std::string& what)
{
    ++vofMaterialComparisons;const auto ref=vofMaterialValue(expected);
    check(std::isfinite(actual)&&std::isfinite(ref)&&vofMaterialSign(expected)>=0,"VOF invalid "+what);
    if(vofMaterialSign(expected)==0){check(actual==0,"VOF true zero changed: "+what);return;}
    const auto e=std::abs(actual-ref)/ref;vofMaterialPeakPhase=std::max(vofMaterialPeakPhase,e);
    std::ostringstream message;message<<std::setprecision(17)<<"VOF "<<what<<" actual="<<actual<<" expected="<<ref<<" relative="<<e;
    check(actual>0&&e<=1e-12L,message.str());
}
VofMaterialBounds vofMaterialLocalBounds(Vector2D lo,Vector2D hi,LAMP::Geometry::PlanarCellCut::Point origin)
{return {vofMaterialSum(lo.x,-origin.x),vofMaterialSum(hi.x,-origin.x),vofMaterialSum(lo.y,-origin.y),vofMaterialSum(hi.y,-origin.y)};}
VofMaterialPolygon vofMaterialClip(VofMaterialPolygon polygon,const VofMaterialBounds& b)
{
    for(int side=0;side<4&&!polygon.empty();++side){
        const int axis=side/2;const auto bound=b[side];const int direction=side%2?-1:1;
        const auto inside=[&](const VofMaterialPoint& p){return direction*vofMaterialSign(vofMaterialSub(p[axis],bound))>=0;};
        VofMaterialPolygon result;auto a=polygon.back();bool ain=inside(a);
        for(const auto& c:polygon){const bool cin=inside(c);
            if(ain!=cin){const auto t=vofMaterialDiv(vofMaterialSub(bound,a[axis]),vofMaterialSub(c[axis],a[axis]));
                VofMaterialPoint cut;cut[axis]=bound;cut[1-axis]=vofMaterialAdd(a[1-axis],vofMaterialMul(t,vofMaterialSub(c[1-axis],a[1-axis])));result.push_back(cut);}
            if(cin)result.push_back(c);a=c;ain=cin;
        }polygon=std::move(result);
    }return polygon;
}
VofMaterialWide vofMaterialGreen(const VofMaterialPolygon& p,long double originX,bool rz)
{
    if(p.size()<3)return {};VofMaterialWide result;const auto& a=p.front();
    for(std::size_t i=1;i+1<p.size();++i){
        const auto bx=vofMaterialSub(p[i][0],a[0]),by=vofMaterialSub(p[i][1],a[1]);
        const auto cx=vofMaterialSub(p[i+1][0],a[0]),cy=vofMaterialSub(p[i+1][1],a[1]);
        const auto det=vofMaterialSub(vofMaterialMul(bx,cy),vofMaterialMul(by,cx));
        check(vofMaterialSign(det)>=0,"VOF reference polygon not convex CCW");
        auto measure=vofMaterialDiv(det,{2,0});
        if(rz){const auto centroidRadius=vofMaterialAdd(vofMaterialAdd({originX,0},a[0]),vofMaterialDiv(vofMaterialAdd(bx,cx),{3,0}));
            measure=vofMaterialMul(vofMaterialMul(measure,centroidRadius),{CartesianMeasure::twoPi,0});}
        result=vofMaterialAdd(result,measure);
    }return result;
}
VofMaterialWide vofMaterialReference(const VofPublication& p,int native,bool liquid,const VofMaterialBounds* b=nullptr)
{
    VofMaterialWide sum;for(const auto& polygon:p.anchoredPhasePolygons(native,liquid)){
        VofMaterialPolygon wide;for(const auto& q:polygon)wide.push_back({{vofMaterialSum(q.anchor.x,q.offset.x),vofMaterialSum(q.anchor.y,q.offset.y)}});
        if(b)wide=vofMaterialClip(std::move(wide),*b);
        sum=vofMaterialAdd(sum,vofMaterialGreen(wide,p.frameOrigin().x,p.domain().axisymmetric()));
    }return sum;
}
long double vofMaterialCapacity(Vector2D lo,Vector2D hi,bool rz)
{return (static_cast<long double>(hi.x)-lo.x)*(static_cast<long double>(hi.y)-lo.y)*(rz?CartesianMeasure::twoPi*(static_cast<long double>(lo.x)+hi.x)/2:1.L);}
std::unique_ptr<QuadtreeMesh> vofMaterialMesh(int n)
{
    DictionaryMesh d;d.Nx=d.Ny=n;d.maxCellNum=4096;d.maxInitialIterations=0;d.maxLevel_local={2};
    auto mesh=std::make_unique<QuadtreeMesh>(d);DictionaryMesh::RefineBox b;b.level=1;b.min={.5,0};b.max={1,1};mesh->applyRefineBoxes({b});
    check(mesh->cellIndexNum()==5*n*n/2,"VOF material fixture leaf count changed");return mesh;
}
long double vofMaterialLine(Vector2D lo,Vector2D hi,bool rz,long double c)
{
    if(c>=static_cast<long double>(hi.x)+hi.y)return vofMaterialCapacity(lo,hi,rz);
    if(c<=static_cast<long double>(lo.x)+lo.y)return 0;
    const auto full=std::clamp(c-hi.y,static_cast<long double>(lo.x),static_cast<long double>(hi.x));
    const auto end=std::clamp(c-lo.y,static_cast<long double>(lo.x),static_cast<long double>(hi.x));
    auto a=(full-lo.x)*(hi.y-lo.y),r=(full*full-lo.x*lo.x)*(hi.y-lo.y)/2;
    if(end>full){const auto h=c-lo.y,s=(end*end-full*full)/2,t=(end*end*end-full*full*full)/3;
        a+=h*(end-full)-s;r+=h*s-t;}
    return rz?CartesianMeasure::twoPi*r:a;
}
VofPublication vofMaterialPublish(const QuadtreeMesh& mesh,bool rz,long double c,int kind=-1)
{
    std::vector<double> alpha(mesh.cellIndexNum());
    for(int ci=0;ci<mesh.cellIndexNum();++ci){const auto& cell=mesh.cells(mesh.cellIndexToLabel[ci]);
        const auto lo=mesh.verts(cell.cellToVert[0]).coord(),hi=mesh.verts(cell.cellToVert[2]).coord();
        const auto cap=vofMaterialCapacity(lo,hi,rz);long double liquid=0;
        if(kind==0)liquid=0;else if(kind==1)liquid=cap;else if(kind==2){auto cut=hi;cut.x=std::min(cut.x,.5);liquid=cut.x<=lo.x?0:vofMaterialCapacity(lo,cut,rz);}
        else liquid=vofMaterialLine(lo,hi,rz,c);
        check(liquid>=0&&liquid<=cap,"VOF material analytic alpha invalid");alpha[ci]=liquid==0?0:liquid==cap?1:static_cast<double>(liquid/cap);
    }
    return VofPublication::prepare(mesh,domainFor(mesh,rz),{7,.125,.25,VofClock::Position::End},alpha);
}
HydroPhysicalDivergence vofMaterialMetrics(const QuadtreeMesh& mesh,const VofPublication& p)
{return HydroPhysicalDivergence::build(mesh,HydroIntegralDivergence::build(mesh,HydroDofLayout::build(mesh,p.domain())));}
bool vofMaterialSameClock(VofClock a,VofClock b)
{return a.step==b.step&&a.begin==b.begin&&a.end==b.end&&a.position==b.position;}

std::array<VofMaterialWide,2> vofMaterialTraceReference(const VofPublication& p,const LAMP::Interface::GeometricVofFace& face,bool minus)
{
    std::array<VofMaterialWide,2> result;
    for(unsigned i=0;i<face.rangeCount;++i){const auto& range=face.ranges[i];
        const auto a=p.anchoredVertices()[range.vertices[0]],b=p.anchoredVertices()[range.vertices[1]];
        const auto ax=vofMaterialSum(a.anchor.x,a.offset.x),ay=vofMaterialSum(a.anchor.y,a.offset.y);
        const auto dx=vofMaterialSub(vofMaterialSum(b.anchor.x,b.offset.x),ax),dy=vofMaterialSub(vofMaterialSum(b.anchor.y,b.offset.y),ay);
        const auto length=vofMaterialSign(dx)==0?dy:dx;
        check(vofMaterialSign(length)>0,"VOF retained native trace orientation invalid");auto physical=length;
        if(p.domain().axisymmetric()){
            const auto radius=vofMaterialAdd(vofMaterialAdd({p.frameOrigin().x,0},ax),vofMaterialDiv(dx,{2,0}));
            physical=vofMaterialMul(vofMaterialMul(physical,radius),{CartesianMeasure::twoPi,0});
        }
        const bool liquid=minus?range.minusLiquid:range.plusLiquid;
        result[liquid?0:1]=vofMaterialAdd(result[liquid?0:1],physical);
    }return result;
}
std::array<VofMaterialWide,2> vofMaterialAnalyticTrace(const QuadtreeMesh& mesh,const NativeLeafFace& f,int native,int kind,bool rz)
{
    std::array<VofMaterialWide,2> result;if(native<0)return result;
    const auto& edge=f.horizontal?mesh.hEdges(f.edgeLabel):mesh.vEdges(f.edgeLabel);
    auto a=mesh.verts(edge.edgeToVert[0]).coord(),b=mesh.verts(edge.edgeToVert[1]).coord();if(f.horizontal?a.x>b.x:a.y>b.y)std::swap(a,b);
    if(kind<2||a.x==b.x){bool liquid=kind==1;
        if(kind==2){const auto& cell=mesh.cells(mesh.cellIndexToLabel[native]);liquid=cell.centerCoord().x<.5;}
        result[liquid?0:1]={f.area,0};return result;
    }
    const auto measure=[&](long double lo,long double hi){
        if(!(hi>lo))return VofMaterialWide{};const auto width=vofMaterialSum(hi,-lo);
        return rz?vofMaterialMul(vofMaterialMul(width,vofMaterialDiv(vofMaterialSum(lo,hi),{2,0})),{CartesianMeasure::twoPi,0}):width;
    };
    result[0]=measure(a.x,std::min(static_cast<long double>(b.x),.5L));
    result[1]=measure(std::max(static_cast<long double>(a.x),.5L),b.x);return result;
}
void vofMaterialFaceContract(const QuadtreeMesh& mesh,const HydroPhaseGeometry& phase,const VofPublication& p,int analytic)
{
    const auto& layout=phase.layout();for(int u=0;u<layout.velocityCount();++u){
        const int row=layout.nativeFace(u);const auto& source=p.faces()[row];const auto& trace=phase.faceTraces()[u];const auto native=layout.domain().face(mesh,row);
        check(trace.negativeCell==source.negativeCell&&trace.positiveCell==source.positiveCell&&trace.coincidentHalfFaces==0,"VOF whole-face metadata fabricated/changed");
        const auto minus=vofMaterialTraceReference(p,source,true),plus=vofMaterialTraceReference(p,source,false);
        std::array<VofMaterialWide,2> spatial;
        if(trace.minusAvailable()){
            vofMaterialNear(trace.minus.liquid,minus[0],"minus face liquid row="+std::to_string(row));
            vofMaterialNear(trace.minus.gas,minus[1],"minus face gas row="+std::to_string(row));
            vofMaterialNear(trace.minus.liquid,{source.liquidMinus,0},"actual source minus liquid");}
        if(trace.plusAvailable()){
            vofMaterialNear(trace.plus.liquid,plus[0],"plus face liquid row="+std::to_string(row));
            vofMaterialNear(trace.plus.gas,plus[1],"plus face gas row="+std::to_string(row));
            vofMaterialNear(trace.plus.liquid,{source.liquidPlus,0},"actual source plus liquid");}
        const auto expectedKind=!source.coincident?HydroPhaseFaceTraceKind::Ordinary:
            trace.minusAvailable()&&trace.plusAvailable()?HydroPhaseFaceTraceKind::InteriorCoincident:HydroPhaseFaceTraceKind::BoundaryCoincident;
        check(trace.kind==expectedKind,"VOF actual whole-face kind changed");vofMaterialCoincidentFaces+=source.coincident;
        for(int k=0;k<2;++k){
            spatial[k]=trace.minusAvailable()&&trace.plusAvailable()?vofMaterialDiv(vofMaterialAdd(minus[k],plus[k]),{2,0}):trace.minusAvailable()?minus[k]:plus[k];
            vofMaterialNear(k==0?phase.faceMeasures()[u].liquid:phase.faceMeasures()[u].gas,spatial[k],"spatial same-face phase");}
        if(analytic>=0){const auto am=vofMaterialAnalyticTrace(mesh,native,source.negativeCell,analytic,p.domain().axisymmetric());
            const auto ap=vofMaterialAnalyticTrace(mesh,native,source.positiveCell,analytic,p.domain().axisymmetric());
            if(trace.minusAvailable()){vofMaterialNear(trace.minus.liquid,am[0],"analytic minus face liquid");vofMaterialNear(trace.minus.gas,am[1],"analytic minus face gas");}
            if(trace.plusAvailable()){vofMaterialNear(trace.plus.liquid,ap[0],"analytic plus face liquid");vofMaterialNear(trace.plus.gas,ap[1],"analytic plus face gas");}
        }
        if(p.domain().axisymmetric()&&native.area==0){++vofMaterialAxisFaces;check(phase.faceMeasures()[u].liquid==0&&phase.faceMeasures()[u].gas==0
            &&trace.minus.liquid==0&&trace.minus.gas==0&&trace.plus.liquid==0&&trace.plus.gas==0,"VOF radial axis physical phase not zero");}
        ++vofMaterialFaceChecks;
    }
}

HydroMomentumMass vofMaterialInspect(const QuadtreeMesh& mesh,const HydroMomentumGeometry& geometry,
    const VofPublication& publication,int analytic=-1)
{
    const auto phase=HydroPhaseGeometry::build(mesh,geometry.metrics(),publication);
    const auto mass=HydroMomentumMass::build(mesh,geometry,phase,100,1);const auto& layout=geometry.layout();const bool rz=publication.domain().axisymmetric();
    check(phase.model()==HydroPhaseGeometry::Model::GeometricVof&&mass.model()==phase.model(),"VOF material model lost");
    check(phase.vofIdentity()==publication.identity()&&mass.vofIdentity()==publication.identity(),"VOF material typed identity changed");
    check(&phase.vofPublication().alpha()==&publication.alpha()&&&mass.phaseGeometry().vofPublication().alpha()==&publication.alpha(),"VOF material detached actual alpha owner");
    check(vofMaterialSameClock(phase.timeLayer(),publication.timeLayer())&&vofMaterialSameClock(mass.timeLayer(),publication.timeLayer()),"VOF material full clock changed");
    phase.validateVofFor(mesh,publication.identity(),publication.timeLayer());mass.validateVofFor(mesh,publication.identity(),publication.timeLayer());
    check(mass.geometry().sameIdentity(geometry)&&mass.layout().sameIdentity(layout)&&phase.layout().sameIdentity(layout)
        &&&mass.geometry().metrics().divergence().matrix()==&phase.metrics().divergence().matrix(),"VOF material native geometry/D changed");
    check(mass.representation()==HydroMomentumMass::DensityRepresentation::SharpGeometry,"VOF actual fan mass projected");
    check(mass.liquidDensity()==100&&mass.gasDensity()==1,"VOF actual density ratio changed");
    vofMaterialFaceContract(mesh,phase,publication,analytic);
    std::array<std::vector<std::array<VofMaterialWide,2>>,2> partition;for(auto& row:partition)row.resize(layout.pressureCount());
    std::array<VofMaterialWide,2> totalMass;std::array<std::array<VofMaterialWide,2>,2> totalPhase{};
    for(int p=0;p<layout.pressureCount();++p){const auto native=layout.nativeCell(p);const auto& source=publication.cells()[native];
        check(phase.liquidFractions()[p]==publication.alpha()[native],"VOF material copied/rederived alpha changed bits");
        vofMaterialNear(phase.cellVolumes()[p].liquid,vofMaterialReference(publication,native,true),"native cell liquid");
        vofMaterialNear(phase.cellVolumes()[p].gas,vofMaterialReference(publication,native,false),"native cell gas");
        const auto minor=std::min(source.inventory,source.capacity-source.inventory);
        if(minor>0){const auto e=std::abs(phase.cellVolumes()[p].liquid-source.inventory)/minor;
            check(e<=1e-12L,"VOF measured liquid/inventory original minor gate");vofMaterialPeakInventory=std::max(vofMaterialPeakInventory,e);
            vofMaterialPeakGasInventory=std::max(vofMaterialPeakGasInventory,std::abs(phase.cellVolumes()[p].gas-(source.capacity-source.inventory))/minor);}
    }
    for(int u=0;u<layout.velocityCount();++u){const auto& volume=geometry.volumes()[u];const int axis=layout.normalAxis(u)==HydroDofLayout::NormalAxis::First?0:1;
        std::array<VofMaterialWide,2> dual;const auto wholeBounds=vofMaterialLocalBounds(volume.lower,volume.upper,publication.frameOrigin());
        long double projected=0;
        check(mass.overlapPhaseVolumes()[u].size()==volume.overlaps.size(),"VOF overlap ordering changed");
        for(std::size_t k=0;k<volume.overlaps.size();++k){const auto& overlap=volume.overlaps[k];const auto native=layout.nativeCell(overlap.pressure);
            const auto b=vofMaterialLocalBounds(overlap.lower,overlap.upper,publication.frameOrigin());
            const auto l=vofMaterialReference(publication,native,true,&b),g=vofMaterialReference(publication,native,false,&b);
            const auto actual=mass.overlapPhaseVolumes()[u][k];vofMaterialNear(actual.liquid,l,"dual overlap liquid");vofMaterialNear(actual.gas,g,"dual overlap gas");
            check(std::abs(actual.liquid+actual.gas-overlap.volume)<=1e-12L*overlap.volume,"VOF actual overlap closure");
            partition[axis][overlap.pressure][0]=vofMaterialAdd(partition[axis][overlap.pressure][0],l);
            partition[axis][overlap.pressure][1]=vofMaterialAdd(partition[axis][overlap.pressure][1],g);
            dual[0]=vofMaterialAdd(dual[0],vofMaterialReference(publication,native,true,&wholeBounds));
            dual[1]=vofMaterialAdd(dual[1],vofMaterialReference(publication,native,false,&wholeBounds));
            projected+=overlap.cellFraction*(100*phase.cellVolumes()[overlap.pressure].liquid+phase.cellVolumes()[overlap.pressure].gas);
        }
        vofMaterialNear(mass.phaseVolumes()[u].liquid,dual[0],"whole dual liquid");vofMaterialNear(mass.phaseVolumes()[u].gas,dual[1],"whole dual gas");
        const auto expected=vofMaterialAdd(vofMaterialMul({100,0},dual[0]),dual[1]);const auto expectedValue=vofMaterialValue(expected);
        const auto e=std::abs(mass.masses()[u]-expectedValue)/expectedValue;vofMaterialPeakMass=std::max(vofMaterialPeakMass,e);
        check(mass.masses()[u]>0&&std::isfinite(mass.masses()[u])&&e<=1e-12L,"VOF actual rho100:1 native mass");
        vofMaterialProjectionGap=std::max(vofMaterialProjectionGap,std::abs(projected-expectedValue)/(101*volume.volume));
        check(std::abs(mass.phaseVolumes()[u].liquid+mass.phaseVolumes()[u].gas-volume.volume)<=1e-12L*volume.volume,"VOF actual dual closure");
        const auto nativeFace=layout.domain().face(mesh,layout.nativeFace(u));
        if(nativeFace.negativeCell>=0&&nativeFace.positiveCell>=0&&vofMaterialSign(dual[0])>0&&vofMaterialSign(dual[1])>0){
            const auto& a=mesh.cells(mesh.cellIndexToLabel[nativeFace.negativeCell]);const auto& b=mesh.cells(mesh.cellIndexToLabel[nativeFace.positiveCell]);
            if(a.size().x!=b.size().x)++vofMaterialAmrCuts[rz?1:0];}
        totalMass[axis]=vofMaterialAdd(totalMass[axis],expected);for(int k=0;k<2;++k)totalPhase[axis][k]=vofMaterialAdd(totalPhase[axis][k],dual[k]);
        if(analytic>=0){auto leftHi=volume.upper,rightLo=volume.lower;long double l=0,g=0;
            if(analytic==0)g=vofMaterialCapacity(volume.lower,volume.upper,rz);
            else if(analytic==1)l=vofMaterialCapacity(volume.lower,volume.upper,rz);
            else {leftHi.x=std::min(leftHi.x,.5);rightLo.x=std::max(rightLo.x,.5);
                if(leftHi.x>volume.lower.x)l=vofMaterialCapacity(volume.lower,leftHi,rz);
                if(rightLo.x<volume.upper.x)g=vofMaterialCapacity(rightLo,volume.upper,rz);}
            vofMaterialNear(mass.phaseVolumes()[u].liquid,{l,0},"analytic pure/vertical dual liquid");
            vofMaterialNear(mass.phaseVolumes()[u].gas,{g,0},"analytic pure/vertical dual gas");}
    }
    for(int axis=0;axis<2;++axis){
        for(int p=0;p<layout.pressureCount();++p)for(int k=0;k<2;++k){const auto ref=vofMaterialReference(publication,layout.nativeCell(p),k==0),sum=partition[axis][p][k];
            vofMaterialNear(vofMaterialValue(sum),ref,"component dual/cell phase partition");
            if(vofMaterialSign(ref)>0)vofMaterialPeakPartition=std::max(vofMaterialPeakPartition,std::abs(vofMaterialValue(vofMaterialSub(sum,ref)))/vofMaterialValue(ref));}
        vofMaterialNear(mass.report().totalPhaseVolumes[axis].liquid,totalPhase[axis][0],"reported liquid total");
        vofMaterialNear(mass.report().totalPhaseVolumes[axis].gas,totalPhase[axis][1],"reported gas total");
        vofMaterialNear(mass.report().totalMass[axis],totalMass[axis],"reported native mass total");
    }
    ++vofMaterialLayers;std::cout<<"V4A_MASS_LAYER,"<<vofMaterialLayers<<','<<rz<<','<<layout.pressureCount()<<','<<layout.velocityCount()<<','
        <<publication.timeLayer().step<<','<<static_cast<int>(publication.timeLayer().position)<<','<<mass.report().cutOverlaps<<",PASS\n";return mass;
}
long double vofMaterialPsi(Vector2D p,int variant)
{
    const long double x=p.x,y=p.y;const auto base=x*x*(1-x)*(1-x)*y*y*(1-y)*(1-y);
    return base*(variant==0?1:variant==1?x+.4L*y+.2L:-2.L);
}
std::vector<double> vofMaterialNormals(const QuadtreeMesh& mesh,const VofPublication& p,int variant)
{
    std::vector<double> result(p.domain().nativeFaceCount());
    for(int row=0;row<p.domain().nativeFaceCount();++row){const auto f=p.domain().face(mesh,row);const auto& edge=f.horizontal?mesh.hEdges(f.edgeLabel):mesh.vEdges(f.edgeLabel);
        auto a=mesh.verts(edge.edgeToVert[0]).coord(),b=mesh.verts(edge.edgeToVert[1]).coord();if(f.horizontal?a.x>b.x:a.y>b.y)std::swap(a,b);
        const auto q=(f.horizontal?-1:1)*(vofMaterialPsi(b,variant)-vofMaterialPsi(a,variant));
        if(f.outwardSign!=0||f.area==0){check(q==0,"VOF material streamfunction boundary leakage");result[row]=0;}else result[row]=static_cast<double>(q/f.area);
    }return result;
}
long double vofMaterialOutflow(const QuadtreeMesh& mesh,const VofPublication& p,const std::vector<double>& a,const std::vector<double>& b)
{
    std::vector<long double> sums(p.cells().size());const auto positive=[](long double x,long double y){
        if(x<=0&&y<=0)return 0.L;if(x>=0&&y>=0)return (x+y)/2;const auto z=std::max(x,y);return z*z/(2*std::abs(y-x));};
    for(int row=0;row<p.domain().nativeFaceCount();++row){const auto f=p.domain().face(mesh,row);
        if(f.negativeCell>=0)sums[f.negativeCell]+=f.area*positive(a[row],b[row]);
        if(f.positiveCell>=0)sums[f.positiveCell]+=f.area*positive(-a[row],-b[row]);}
    long double maximum=0;for(const auto& c:p.cells())maximum=std::max(maximum,sums[c.nativeCell]/c.capacity);check(maximum>0,"VOF material slab zero motion");return maximum;
}
LAMP::Interface::GeometricVofAdvectorSlab vofMaterialSlab(const QuadtreeMesh& mesh,const VofPublication& p,const HydroDofLayout& layout,
    VofClock next,const std::vector<double>& a,const std::vector<double>& b)
{
    std::vector<double> av(layout.velocityCount()),bv(av.size());for(int u=0;u<layout.velocityCount();++u){av[u]=a[layout.nativeFace(u)];bv[u]=b[layout.nativeFace(u)];}
    const HydroNativeValues begin(layout,av,std::vector<double>(layout.pressureCount(),.125)),end(layout,bv,std::vector<double>(layout.pressureCount(),-.25));
    const auto result=prepareHydroVofAdvectorSlab(mesh,p,next,begin,end);
    check(result.beginNormal()==a&&result.endNormal()==b&&begin.normalVelocity()==av&&end.normalVelocity()==bv,"VOF material native slab input changed");return result;
}
template<class F>void vofMaterialReject(F action,const std::string& what)
{
    bool caught=false;try{action();}catch(const std::invalid_argument&){caught=true;}
    check(caught,"missing explicit VOF material rejection: "+what);++vofMaterialGuards;
}
void vofMaterialGuardContract(QuadtreeMesh& mesh,const HydroMomentumGeometry& geometry,const VofPublication& p,const HydroMomentumMass& mass,bool mutate=false)
{
    const auto saved=p.alpha(),savedMass=mass.masses();const auto identity=p.identity();const auto clock=p.timeLayer();const auto& phase=mass.phaseGeometry();
    const auto foreign=VofPublication::prepare(mesh,p.domain(),clock,p.alpha());
    vofMaterialReject([&]{phase.validateVofFor(mesh,foreign.identity(),clock);},"same clock foreign owner");
    vofMaterialReject([&]{mass.validateVofFor(mesh,foreign.identity(),clock);},"same clock foreign mass owner");
    auto wrong=clock;wrong.begin-=.0625;vofMaterialReject([&]{mass.validateVofFor(mesh,identity,wrong);},"full clock begin");
    wrong=clock;wrong.position=VofClock::Position::Midpoint;vofMaterialReject([&]{phase.validateVofFor(mesh,identity,wrong);},"full clock position");
    vofMaterialReject([&]{(void)phase.publication();},"VOF cannot pretend FirstJet publication");
    vofMaterialReject([&]{(void)phase.geometryIdentity();},"VOF cannot pretend FirstJet identity");
    vofMaterialReject([&]{(void)mass.geometryIdentity();},"VOF mass cannot pretend FirstJet identity");
    vofMaterialReject([&]{HydroMomentumMass::buildCellAverageProjection(mesh,geometry,phase,100,1);},"VOF cannot use cell average projection");
    for(double rho:{0.,-1.,std::numeric_limits<double>::infinity(),std::numeric_limits<double>::quiet_NaN()})
        vofMaterialReject([&]{HydroMomentumMass::build(mesh,geometry,phase,rho,1);},"invalid density");
    vofMaterialReject([&]{HydroPhaseGeometry::build(mesh,metricsFor(mesh,!p.domain().axisymmetric()),p);},"XY/RZ measure");
    vofMaterialReject([&]{HydroPhaseGeometry::build(mesh,metricsFor(mesh,p.domain().axisymmetric()),p);},"foreign equivalent native mesh domain");
    const auto native=geometry.layout().nativeCell(0);const auto& cell=mesh.cells(mesh.cellIndexToLabel[native]);
    const auto lo=mesh.verts(cell.cellToVert[0]).coord(),hi=mesh.verts(cell.cellToVert[2]).coord();
    HydroVofBox bad{native,{lo,hi}};bad.rectangle.lower.x-=.125;
    vofMaterialReject([&]{phase.measureVofBoxes(mesh,{bad});},"outside selected native leaf box");
    bad={-1,{lo,hi}};vofMaterialReject([&]{phase.measureVofBoxes(mesh,{bad});},"invalid native leaf box");
    bad={native,{lo,hi}};bad.rectangle.lower.x=hi.x;vofMaterialReject([&]{phase.measureVofBoxes(mesh,{bad});},"zero width native leaf box");
    const auto otherD=HydroPhysicalDivergence::build(mesh,HydroIntegralDivergence::build(mesh,geometry.layout()));
    const auto otherPhase=HydroPhaseGeometry::build(mesh,otherD,p);
    vofMaterialReject([&]{HydroMomentumMass::build(mesh,geometry,otherPhase,100,1);},"different retained D_I");
    if(mutate){DictionaryMesh::RefineBox b;b.level=1;b.min={0,.5};b.max={.5,1};mesh.applyRefineBoxes({b});
        vofMaterialReject([&]{mass.validateVofFor(mesh,identity,clock);},"mesh revision");}
    check(p.alpha()==saved&&p.identity()==identity&&vofMaterialSameClock(p.timeLayer(),clock)&&mass.masses()==savedMass,"VOF material rejection polluted old source/mass");
}
void vofMaterialMassContract()
{
    for(int n:{4,8})for(bool rz:{false,true})for(long double c:{.61L,.913L,1.19L}){
        auto mesh=vofMaterialMesh(n);const auto old=vofMaterialPublish(*mesh,rz,c);const auto saved=old.alpha();const auto identity=old.identity();const auto clock=old.timeLayer();
        const auto geometry=HydroMomentumGeometry::build(*mesh,vofMaterialMetrics(*mesh,old));const auto oldMass=vofMaterialInspect(*mesh,geometry,old);const auto savedMass=oldMass.masses();
        const auto a=vofMaterialNormals(*mesh,old,0),b=vofMaterialNormals(*mesh,old,1),corrected=vofMaterialNormals(*mesh,old,2);
        const auto out=std::max(vofMaterialOutflow(*mesh,old,a,b),vofMaterialOutflow(*mesh,old,a,corrected));
        const auto dt=static_cast<double>(.2L/out);const VofClock next{8,.25,.25+dt,VofClock::Position::End};
        const auto slabA=vofMaterialSlab(*mesh,old,geometry.layout(),next,a,b),slabC=vofMaterialSlab(*mesh,old,geometry.layout(),next,a,corrected);
        const auto A=VofTransport::prepare(*mesh,old,slabA,100,1),C=VofTransport::prepare(*mesh,old,slabC,100,1);vofMaterialSlabs+=2;
        for(const auto* transport:{&A,&C}){const auto mid=next;auto half=mid;half.position=VofClock::Position::Midpoint;
            check(transport->sourceIdentity()==identity&&vofMaterialSameClock(transport->sourceTime(),clock),"VOF material slab did not start from same old owner");
            check(vofMaterialSameClock(transport->midpoint().timeLayer(),half)&&vofMaterialSameClock(transport->endpoint().timeLayer(),next),"VOF material candidates lost complete interval");
            check(transport->midpoint().identity()!=identity&&transport->endpoint().identity()!=identity&&transport->midpoint().identity()!=transport->endpoint().identity(),"VOF material candidate reused source owner");
            vofMaterialInspect(*mesh,geometry,transport->midpoint());vofMaterialInspect(*mesh,geometry,transport->endpoint());}
        vofMaterialDifferentCandidates=vofMaterialDifferentCandidates||A.endpoint().alpha()!=C.endpoint().alpha();
        check(old.alpha()==saved&&old.identity()==identity&&vofMaterialSameClock(old.timeLayer(),clock)&&oldMass.masses()==savedMass,"VOF material candidate mutated accepted source");
        if(n==4&&c==.913L)vofMaterialGuardContract(*mesh,geometry,old,oldMass);
    }
    for(int n:{4,8})for(bool rz:{false,true})for(int kind:{0,1,2}){
        auto mesh=vofMaterialMesh(n);const auto old=vofMaterialPublish(*mesh,rz,0,kind);const auto geometry=HydroMomentumGeometry::build(*mesh,vofMaterialMetrics(*mesh,old));
        const auto mass=vofMaterialInspect(*mesh,geometry,old,kind);if(n==8&&rz&&kind==2)vofMaterialGuardContract(*mesh,geometry,old,mass,true);
    }
    check(vofMaterialLayers==72&&vofMaterialSlabs==24,"VOF material frozen matrix coverage changed");
    check(vofMaterialAmrCuts[0]>0&&vofMaterialAmrCuts[1]>0,"VOF material missing real XY/RZ mixed coarse-fine duals");
    check(vofMaterialProjectionGap>1e-4L,"VOF material fixture missed cell-average negative control");
    check(vofMaterialDifferentCandidates,"VOF material A/C candidate velocities did not differ");
    check(vofMaterialCoincidentFaces>0&&vofMaterialAxisFaces>0,"VOF face coincidence/axis coverage missing");
    std::cout<<"V4A_VOF_NATIVE_FACES,"<<vofMaterialFaceChecks<<','<<vofMaterialCoincidentFaces<<','<<vofMaterialAxisFaces<<",PASS\n";
    std::cout<<"V4A_VOF_NATIVE_MASS,"<<vofMaterialLayers<<','<<vofMaterialSlabs<<','<<vofMaterialComparisons<<','<<vofMaterialGuards<<','
        <<vofMaterialAmrCuts[0]<<','<<vofMaterialAmrCuts[1]<<','<<vofMaterialPeakPhase<<','<<vofMaterialPeakMass<<','<<vofMaterialPeakPartition<<','
        <<vofMaterialPeakInventory<<','<<vofMaterialPeakGasInventory<<','<<vofMaterialProjectionGap<<",PASS\n";
}
#if defined(_MSC_VER)
#pragma float_control(pop)
#endif

}
int main()
{
    try {
        std::cout<<std::setprecision(17);advectorClosureContract();regularArithmeticContract();advectorBoundaryGuards();
        for(bool radial : {false,true})for(int amr : {0,1}) {
            auto mesh=makeMesh(4,amr);
            for(int kind=0;kind<4;++kind) {
                Case s;s.radial=radial;s.amr=amr!=0;s.pureFirstJet=true;
                s.planar=kind!=0;s.planeOffset=kind==1?.43:kind==2?-1:2;
                const auto g=HydroPhaseGeometry::build(*mesh,metricsFor(*mesh,radial),publishAt(*mesh,s,{}));
                inspect(*mesh,s,g,true);
                if(kind>=2)for(const auto& m:g.cellVolumes())check(kind==2?m.liquid==0:m.gas==0,"pure cell not exact");
            }
            Case subset;subset.radial=radial;subset.amr=amr!=0;subset.pureFirstJet=true;
            inspect(*mesh,subset,HydroPhaseGeometry::build(*mesh,metricsFor(*mesh,radial,true),publishAt(*mesh,subset,{})),true);
        }
        for(bool gas : {false,true}) {
            auto mesh=makeMesh(4);Case s;s.pureFirstJet=true;s.planar=true;
            s.planeOffset=gas?1-std::ldexp(1.,-40):std::ldexp(1.,-40);
            const auto g=HydroPhaseGeometry::build(*mesh,metricsFor(*mesh,false),publishAt(*mesh,s,{}));
            inspect(*mesh,s,g,true);
            int tiny=0;for(const auto& m:g.faceMeasures()) {
                const auto small=gas?m.gas:m.liquid;
                if(small>0 && small<1e-10){close(small,std::ldexp(1.L,-40),1e-24,"retained tiny interval");++tiny;}
            }
            check(tiny>0,"tiny face phase disappeared");std::cout<<"TINY,"<<gas<<','<<tiny<<'\n';
        }
        {
            auto mesh=makeMesh(4);Case s;s.pureFirstJet=true;s.cx=.42;s.cy=.39;s.radius=.105;
            const auto g=HydroPhaseGeometry::build(*mesh,metricsFor(*mesh,false),publishAt(*mesh,s,{}));
            int doubleRootFaces=0;
            for(int i=0;i<g.layout().velocityCount();++i) {
                const int native=g.layout().nativeFace(i);
                const auto& partitions=g.publication().facePartitions();
                if(partitions.halfFace(native,0).roots().size()+partitions.halfFace(native,1).roots().size()>=2)++doubleRootFaces;
            }
            check(doubleRootFaces>0,"two-root native face was not exercised");
            inspect(*mesh,s,g,true);std::cout<<"DOUBLE_ROOT,"<<doubleRootFaces<<'\n';
        }
        for(bool radial : {false,true})for(int amr : {0,1}) {
            std::vector<Error> errors;
            for(int n : {4,8,16}) {
                auto mesh=makeMesh(n,amr);Case s;s.radial=radial;s.amr=amr!=0;s.nonpolynomial=true;s.pureFirstJet=true;
                const auto g=HydroPhaseGeometry::build(*mesh,metricsFor(*mesh,radial),publishAt(*mesh,s,{}));
                const auto error=inspect(*mesh,s,g,false);errors.push_back(error);++convergenceChecks;
                std::cout<<"CONVERGENCE,"<<radial<<','<<amr<<','<<n<<','<<error.cell<<','<<error.face<<'\n';
            }
            for(int i=1;i<3;++i)check(errors[i].cell<errors[i-1].cell && errors[i].face<errors[i-1].face,"phase measurement error did not decrease");
            const auto cellOrder=std::log2(errors[1].cell/errors[2].cell),faceOrder=std::log2(errors[1].face/errors[2].face);
            check(cellOrder>=2 && faceOrder>=2,"phase measurement convergence order below gate");
            std::cout<<"ORDER,"<<radial<<','<<amr<<','<<cellOrder<<','<<faceOrder<<'\n';
        }
        rejections();
        std::cout<<"SUMMARY,"<<scenarios<<','<<cellChecks<<','<<faceChecks<<','<<coarseChecks<<','<<rejectionChecks<<','
            <<convergenceChecks<<','<<maximumCell<<','<<maximumFace<<','<<maximumRate<<','<<maximumClosure<<'\n';
        std::cout<<"H07a HYDRO same-geometry phase rates PASS\n";
        nativePhaseTraces();
        spatialPhaseContract();
        spatialPhaseShapes();
        spatialPhaseGuards();
        completeNormalRateContract();
        vofMaterialMassContract();
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
