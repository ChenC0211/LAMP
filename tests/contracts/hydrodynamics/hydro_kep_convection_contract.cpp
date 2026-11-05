#include "Modules/Hydrodynamics/HydroKepConvection.h" // Public header stands alone.
#include "Modules/Hydrodynamics/Detail/HydroFiniteContraction.h"
#include "interface_first_jet_fixture.h"
#include <boost/multiprecision/cpp_int.hpp>
#include <sstream>
#include <set>
#include <map>
#include <tuple>
#include <cfenv>
#include <cerrno>
#include <cstring>
#include <exception>
#if defined(_M_X64) || defined(__SSE2__)
#include <xmmintrin.h>
#endif

namespace
{
using namespace LAMP::Hydrodynamics;
int scenarios=0,traceChecks=0,coarseFineChecks=0,alignedChecks=0,loadChecks=0,actualCrossings=0,rejected=0,regrids=0;
long double peakTrace=0,peakLoad=0,peakMomentum=0,peakEnergy=0,peakConstant=0,peakLinearity=0;
long double amrMeanGap=0,amrRawPowerGap=0,actualDefect=0,actualTransfer=0;
int momentChecks=0,pairChecks=0;long double peakMoment=0;
int migrationFailures=0;
int controlledCases=0,controlledConnections=0,controlledActive=0,controlledQuiet=0,controlledCoarseFine=0;
long double controlledPeak=0,controlledPeakDefect=0;
int materialEnergyCases=0,materialEnergyDofs=0;
int impermeableCases=0,impermeableFreeRows=0,impermeableWallTangentRows=0,frozenOperatorCases=0;
int nativeMeanP1AnalyticCases=0,nativeMeanP1OldPointNegative=0,nativeMeanP1AffineChecks=0;
int physicalPairSpatialFamilies=0,physicalPairSpatialRows=0,physicalPairWallRows=0;
std::array<int,4> nativeMeanP2Negative{}; // Truth zero/sign/xy-factor and translated-bound controls; separate from old reject counts.
std::array<int,3> materialEnergyNegative{};
std::array<long double,5> materialEnergyError{};
std::array<long double,3> materialEnergyMutation{};
static_assert(!std::is_default_constructible_v<HydroKepConvection>);
static_assert(!std::is_default_constructible_v<HydroConvectionBoundary>);
static_assert(std::is_same_v<decltype(std::declval<HydroKepConvectionResult&>().momentumRates()),const std::vector<long double>&>);
static_assert(std::is_same_v<decltype(std::declval<HydroKepConvection&>().couplings()),const std::vector<HydroConvectionCoupling>&>);
long double error(long double a,long double b,long double scale,long double limit,const std::string& why)
{
    const auto e=std::abs(a-b)/(scale>0?scale:1);std::ostringstream message;
    message<<std::setprecision(17)<<why<<" actual="<<a<<" expected="<<b<<" scale="<<scale<<" error="<<e;
    check(std::isfinite(a) && std::isfinite(b) && std::isfinite(e) && e<=limit,message.str());return e;
}
int axis(const HydroDofLayout& l,int u){return l.normalAxis(u)==HydroDofLayout::NormalAxis::First?0:1;}
bool inside(const HydroMomentumSegment& g){return g.negativeVelocity>=0 && g.positiveVelocity>=0;}
int owner(const HydroMomentumSegment& g){return g.negativeVelocity>=0?g.negativeVelocity:g.positiveVelocity;}
bool coarseFine(const QuadtreeMesh& mesh,const HydroMomentumGeometry& geometry,const HydroMomentumSegment& g)
{
    int lo=100,hi=-1;for(int u:{g.negativeVelocity,g.positiveVelocity})for(const auto& overlap:geometry.volumes()[u].overlaps){
        const int level=mesh.cells(mesh.cellIndexToLabel[geometry.layout().nativeCell(overlap.pressure)]).level();lo=std::min(lo,level);hi=std::max(hi,level);}
    return lo<hi;
}
HydroPhysicalDivergence metrics(const QuadtreeMesh& mesh,bool rz,int selection=0)
{
    FaceSegmentGraph graph;graph.build(mesh,[selection](const QuadtreeMesh&,const Cell& c){const auto p=c.centerCoord();
        if(selection==1)return p.x<.5;if(selection==2)return !(p.x>.25 && p.x<.75 && p.y>.25 && p.y<.75);
        if(selection==3)return p.x<.25 || p.x>.75;return true;});
    return HydroPhysicalDivergence::build(mesh,HydroIntegralDivergence::build(mesh,HydroDofLayout::build(mesh,NativeLeafDomain::build(mesh,graph,rz))));
}
HydroNativeValues values(const QuadtreeMesh& mesh,const HydroDofLayout& l,const std::function<Vector2D(Vector2D)>& f)
{
    std::vector<double> u;for(int i=0;i<l.velocityCount();++i){const auto face=l.domain().face(mesh,l.nativeFace(i));const auto v=f(face.center);u.push_back(face.horizontal?v.y:v.x);}
    return HydroNativeValues(l,u,std::vector<double>(l.pressureCount()));
}
Vector2D centroid(const HydroMomentumSegment& g,bool rz)
{
    auto p=g.begin+(g.end-g.begin)*.5;
    if(rz && g.normalAxis==1){const long double a=g.begin.x,b=g.end.x;p.x=static_cast<double>((2.L/3)*(b*b*b-a*a*a)/(b*b-a*a));}
    return p;
}
std::vector<double> traces(const HydroDualMaterial& m,const std::function<Vector2D(Vector2D)>& field)
{
    std::vector<double> b(m.faceRates().size());for(std::size_t g=0;g<b.size();++g)if(!inside(m.geometry().segments()[g])){
        const auto& s=m.geometry().segments()[g];const auto p=centroid(s,m.layout().domain().axisymmetric());const auto v=field(p);b[g]=axis(m.layout(),owner(s))==0?v.x:v.y;
    }return b;
}
HydroConvectionBoundary bindBoundary(const QuadtreeMesh& mesh,const HydroDualMaterial& m,const std::vector<double>& t)
{return HydroConvectionBoundary::bind(mesh,m,t,m.sources(),m.timeLayer());}
HermiteGalsPublication publishAt(const QuadtreeMesh& mesh,const Case& s,HermiteGalsTimeLayer time)
{
    const auto f=s.profile(time.time());const auto seed=stateFor(mesh,f,s.radial,time,[](int,Vector2D,Vector2D){return 0.;},s.symmetry(mesh),HermiteGalsTransport::FirstJet);
    HermiteGalsCandidate candidate(mesh,time,seed.primary());const auto mixed=closeHermiteGalsMixedAmr(mesh,candidate,s.symmetry(mesh));
    HermiteGalsPublicationConfig config;config.volume.symmetryPlanes=s.symmetry(mesh);const auto id=candidate.identity();
    return publishHermiteGalsFirstJetGeometry(mesh,std::move(candidate),mixed,{id,centers(mesh,f)},config);
}
struct Frames {HydroPhaseGeometry a,m,b;HydroMomentumMass ma,mm,mb;};
Frames frames(const QuadtreeMesh& mesh,const HydroMomentumGeometry& g,HermiteGalsPublication a,HermiteGalsPublication m,HermiteGalsPublication b,double rho)
{
    const auto pa=HydroPhaseGeometry::build(mesh,g.metrics(),a),pm=HydroPhaseGeometry::build(mesh,g.metrics(),m),pb=HydroPhaseGeometry::build(mesh,g.metrics(),b);
    return {pa,pm,pb,HydroMomentumMass::build(mesh,g,pa,rho,1),HydroMomentumMass::build(mesh,g,pm,rho,1),HydroMomentumMass::build(mesh,g,pb,rho,1)};
}
HydroDualMaterial material(const QuadtreeMesh& mesh,const Frames& f,const HydroNativeValues& u,const HydroAdvectorBoundary* boundary=nullptr)
{
    const auto reconstruction=HydroAdvectorReconstruction::prepare(mesh,f.mm.geometry());
    const auto actual=boundary?reconstruction.evaluate(mesh,u,f.m.geometryIdentity(),f.m.timeLayer(),*boundary):
        reconstruction.evaluate(mesh,u,f.m.geometryIdentity(),f.m.timeLayer());
    const auto q=HydroPhaseRate::prepareSpatial(mesh,f.m,actual);
    return HydroDualMaterial::prepare(mesh,HydroMaterialTransport::prepare(mesh,f.a,q,f.b,f.ma.liquidDensity(),1),f.ma,f.mm,f.mb);
}
std::vector<long double> outflow(const HydroDualMaterial& m)
{
    std::vector<long double> q(m.layout().velocityCount());for(std::size_t g=0;g<m.faceRates().size();++g){const auto& s=m.geometry().segments()[g];
        if(s.negativeVelocity>=0)q[s.negativeVelocity]+=m.faceRates()[g];if(s.positiveVelocity>=0)q[s.positiveVelocity]-=m.faceRates()[g];}return q;
}
void nativeMoments(const QuadtreeMesh& mesh,const HydroKepConvection& op)
{
    const auto& m=op.material();const auto& l=op.layout();const int n=l.velocityCount();
    std::vector<Vector2D> points;std::vector<std::vector<int>> cellFaces(l.pressureCount()),neighbors(l.pressureCount());
    for(int u=0;u<n;++u){const auto face=l.domain().face(mesh,l.nativeFace(u));points.push_back(face.center);
        const int a=face.negativeCell<0?-1:l.pressureDof(face.negativeCell),b=face.positiveCell<0?-1:l.pressureDof(face.positiveCell);
        if(a>=0)cellFaces[a].push_back(u);if(b>=0)cellFaces[b].push_back(u);
        if(a>=0 && b>=0){neighbors[a].push_back(b);neighbors[b].push_back(a);}}
    std::set<std::pair<int,int>> allowed;
    for(int u=0;u<n;++u){std::set<int> roots,patch,support;
        for(const auto& p:m.geometry().volumes()[u].overlaps)roots.insert(p.pressure);patch=roots;
        for(int p:roots)patch.insert(neighbors[p].begin(),neighbors[p].end());
        for(int p:patch)for(int v:cellFaces[p])if(axis(l,v)==axis(l,u))support.insert(v);
        for(int i:support)for(int j:support)if(i<j)allowed.emplace(i,j);}
    std::vector<long double> kernel(static_cast<std::size_t>(n)*n),correction(kernel.size());std::vector<bool> exterior(n);
    for(std::size_t g=0;g<m.faceRates().size();++g){const auto& s=m.geometry().segments()[g];if(!inside(s)){exterior[owner(s)]=true;continue;}
        for(const auto& w:op.weights()[g]){kernel[static_cast<std::size_t>(s.negativeVelocity)*n+w.velocity]+=.5L*m.faceRates()[g]*w.weight;kernel[static_cast<std::size_t>(s.positiveVelocity)*n+w.velocity]-=.5L*m.faceRates()[g]*w.weight;}}
    for(int i=0;i<n;++i)for(int j=i;j<n;++j){const auto value=kernel[static_cast<std::size_t>(i)*n+j]-kernel[static_cast<std::size_t>(j)*n+i];kernel[static_cast<std::size_t>(i)*n+j]=value;kernel[static_cast<std::size_t>(j)*n+i]=-value;}
    for(const auto& e:op.couplings()){check(allowed.count({e.first,e.second})==1,"native correction left geometry one-ring support");++pairChecks;
        correction[static_cast<std::size_t>(e.first)*n+e.second]+=e.coefficient;correction[static_cast<std::size_t>(e.second)*n+e.first]-=e.coefficient;}
    for(int i=0;i<n;++i){const auto p=points[i];const auto& box=m.geometry().volumes()[i];const double h=std::max(box.upper.x-box.lower.x,box.upper.y-box.lower.y);
        std::array<long double,6> moments{},gross{},target{};
        for(int j=0;j<n;++j){const auto x=(points[j].x-p.x)/h,y=(points[j].y-p.y)/h;const std::array<double,6> q{1,x,y,x*x,x*y,y*y};
            for(int c=0;c<6;++c){const auto k=correction[static_cast<std::size_t>(i)*n+j]+(c?kernel[static_cast<std::size_t>(i)*n+j]:0);moments[c]+=k*q[c];gross[c]+=std::abs(k*q[c]);}}
        if(!exterior[i])for(std::size_t g=0;g<m.faceRates().size();++g){const auto& s=m.geometry().segments()[g];const int sign=s.negativeVelocity==i?1:(s.positiveVelocity==i?-1:0);if(!sign)continue;
            const auto at=centroid(s,l.domain().axisymmetric());target[1]+=sign*m.faceRates()[g]*(at.x-p.x)/h;target[2]+=sign*m.faceRates()[g]*(at.y-p.y)/h;}
        for(int c=0;c<(exterior[i]?1:6);++c){peakMoment=std::max(peakMoment,error(moments[c],target[c],gross[c]+std::abs(target[c]),1e-12,"independent final native spatial moment"));++momentChecks;}
    }
}
void affineTraces(const QuadtreeMesh& mesh,const HydroKepConvection& op,bool amr)
{
    const auto& l=op.layout();const bool rz=l.domain().axisymmetric();
    for(std::size_t g=0;g<op.weights().size();++g){const auto& segment=op.material().geometry().segments()[g];if(!inside(segment))continue;
        const auto p=centroid(segment,rz);long double one=0,x=0,y=0,gross=0;
        for(const auto& w:op.weights()[g]){const auto f=l.domain().face(mesh,l.nativeFace(w.velocity));check(axis(l,w.velocity)==axis(l,segment.negativeVelocity),"mixed transported components");
            one+=w.weight;x+=w.weight*f.center.x;y+=w.weight*f.center.y;gross+=std::abs(w.weight);}
        peakTrace=std::max({peakTrace,error(one,1,gross+1,1e-12,"constant trace"),error(x,p.x,gross+1,1e-12,"x affine trace"),error(y,p.y,gross+1,1e-12,"y affine trace")});
        error(op.points()[g].x,p.x,1,1e-14,"independent physical centroid");error(op.points()[g].y,p.y,1,1e-14,"independent physical centroid");
        const auto a=l.domain().face(mesh,l.nativeFace(segment.negativeVelocity)).center,b=l.domain().face(mesh,l.nativeFace(segment.positiveVelocity)).center;
        const auto mid=a+(b-a)*.5;const auto& volume=op.material().geometry().volumes()[segment.negativeVelocity];
        const double h=std::max(volume.upper.x-volume.lower.x,volume.upper.y-volume.lower.y);
        if(amr && coarseFine(mesh,op.material().geometry(),segment)){amrMeanGap=std::max(amrMeanGap,static_cast<long double>(std::hypot(p.x-mid.x,p.y-mid.y)/h));++coarseFineChecks;}
        if(!rz && p.x==mid.x && p.y==mid.y){for(const auto& w:op.weights()[g])check((w.velocity==segment.negativeVelocity || w.velocity==segment.positiveVelocity) && w.weight==.5,"aligned XY trace changed arithmetic center");check(op.weights()[g].size()==2,"aligned XY trace support changed");++alignedChecks;}++traceChecks;
    }
}
// Physical mass differences and dual-side incidence, independent of the
// production defect/report and of the dense convection reference below.
void materialEnergyReference(const HydroDualMaterial& m,const HydroNativeValues& u,
    const HydroKepConvectionResult& result,const std::vector<double>& trace)
{
    std::array<long double,2> convection{},outwardPower{},kinetic{},massPower{},materialPower{},boundaryPower{},scale{};
    for(int i=0;i<m.layout().velocityCount();++i){const int a=axis(m.layout(),i);const long double v=u.normalVelocity()[i];
        const auto begin=m.beginMass().phaseVolumes()[i],end=m.endMass().phaseVolumes()[i];
        const auto dm=m.material().liquidDensity()*(end.liquid-begin.liquid)+m.material().gasDensity()*(end.gas-begin.gas);
        long double q=0,gross=0;for(int side=0;side<4;++side)for(int id:m.geometry().volumes()[i].sides[side]){
            const auto f=m.faceRates()[id];q+=(side%2?1:-1)*f;gross+=std::abs(f);}
        const auto mdot=dm/m.deltaTime(),r=mdot+q,rateScale=std::max(1.L,std::abs(mdot)+gross);
        materialEnergyError[0]=std::max(materialEnergyError[0],error(m.cells()[i].massChange,dm,rateScale*m.deltaTime(),1e-12,"independent geometry mass change"));
        materialEnergyError[1]=std::max({materialEnergyError[1],error(m.cells()[i].outwardRate,q,rateScale,1e-12,"independent dual-side outflow"),
            error(m.cells()[i].defectRate,r,rateScale,1e-12,"independent actual material R")});
        const auto c=result.momentumRates()[i],load=result.momentumLoads()[i];
        materialEnergyError[2]=std::max(materialEnergyError[2],error(load,-c+v*r,std::max(1.L,std::abs(c)+std::abs(v)*(std::abs(mdot)+gross)),1e-12,"independent actual uR load"));
        convection[a]+=v*c;outwardPower[a]+=.5L*v*v*q;kinetic[a]+=v*load;
        massPower[a]+=.5L*v*v*mdot;materialPower[a]+=.5L*v*v*r;
        scale[a]+=std::abs(v*c)+std::abs(v*load)+v*v*(std::abs(mdot)+gross);++materialEnergyDofs;
    }
    for(std::size_t id=0;id<m.faceRates().size();++id){const auto& s=m.geometry().segments()[id];if(inside(s))continue;
        const int i=owner(s),a=axis(m.layout(),i);const long double v=u.normalVelocity()[i],f=(s.negativeVelocity>=0?1:-1)*m.faceRates()[id];
        boundaryPower[a]+=f*(v*trace[id]-.5L*v*v);scale[a]+=std::abs(f*v*trace[id])+std::abs(.5L*f*v*v);
    }
    for(int a=0;a<2;++a){const auto budget=std::max(1.L,scale[a]),expected=materialPower[a]-boundaryPower[a],actual=kinetic[a]-massPower[a];
        materialEnergyError[3]=std::max(materialEnergyError[3],error(convection[a],outwardPower[a]+boundaryPower[a],budget,1e-12,"independent actual q convection power"));
        materialEnergyError[4]=std::max(materialEnergyError[4],error(actual,expected,budget,1e-12,"independent actual R kinetic budget"));
        const std::array<long double,3> wrong{actual-2*materialPower[a],actual-4*materialPower[a],convection[a]};
        const char* names[3]={"missing uR energy mutation","reversed uR energy mutation","zero q energy mutation"};
        for(int k=0;k<3;++k){const auto target=k==2?boundaryPower[a]:expected;const auto gap=std::abs(wrong[k]-target)/budget;
            materialEnergyMutation[k]=std::max(materialEnergyMutation[k],gap);
            if(!materialEnergyNegative[k]&&gap>1e-5){rejects([&]{error(wrong[k],target,budget,1e-12,names[k]);},names[k]);++materialEnergyNegative[k];}}
    }
    ++materialEnergyCases;
}

// Independent dense raw matrix and final shared skew coefficients. This
// does not reuse the production per-segment/pair load scatter or its report.
HydroKepConvectionResult inspect(const QuadtreeMesh& mesh,const HydroKepConvection& op,
    const HydroNativeValues& u,const std::vector<double>& t,bool amr=false,bool ordinary=false)
{
    const auto& m=op.material();const auto bound=bindBoundary(mesh,m,t);const auto result=op.apply(mesh,u,bound,m.sources(),m.timeLayer());const int n=op.layout().velocityCount();
    std::vector<long double> raw(static_cast<std::size_t>(n)*n),rawAbs(raw.size()),rawCf(raw.size()),qInternal(n),qCf(n),boundary(n),expected(n);std::array<long double,2> bMomentum{},bPower{},power{},massPower{},materialPower{},momentum{},scale{},momentumScale{};
    for(std::size_t g=0;g<m.faceRates().size();++g){const auto& s=m.geometry().segments()[g];const auto f=m.faceRates()[g];long double trace=0;
        if(inside(s)){
            const int i=s.negativeVelocity,j=s.positiveVelocity;qInternal[i]+=f;qInternal[j]-=f;
            const bool cf=amr && coarseFine(mesh,m.geometry(),s);if(cf){qCf[i]+=f;qCf[j]-=f;}
            for(const auto& w:op.weights()[g]){raw[static_cast<std::size_t>(i)*n+w.velocity]+=f*w.weight;raw[static_cast<std::size_t>(j)*n+w.velocity]-=f*w.weight;rawAbs[static_cast<std::size_t>(i)*n+w.velocity]+=std::abs(f*w.weight);rawAbs[static_cast<std::size_t>(j)*n+w.velocity]+=std::abs(f*w.weight);trace+=w.weight*u.normalVelocity()[w.velocity];
                if(cf){rawCf[static_cast<std::size_t>(i)*n+w.velocity]+=f*w.weight;rawCf[static_cast<std::size_t>(j)*n+w.velocity]-=f*w.weight;}}
            if(ordinary){const auto p=centroid(s,false),a=op.layout().domain().face(mesh,op.layout().nativeFace(i)).center,b=op.layout().domain().face(mesh,op.layout().nativeFace(j)).center;const auto mid=a+(b-a)*.5;
                if(p.x==mid.x && p.y==mid.y)error(result.centerMomentumRates()[g],f*(static_cast<long double>(u.normalVelocity()[i])+u.normalVelocity()[j])/2,std::max(1.L,std::abs(f*trace)),1e-12,"aligned ordinary XY center flux");}
        }else{
            const int i=owner(s),a=axis(op.layout(),i);const long double sign=s.negativeVelocity>=0?1:-1;trace=t[g];boundary[i]+=sign*f*trace;
            bMomentum[a]+=sign*f*trace;bPower[a]+=sign*f*(u.normalVelocity()[i]*trace-.5L*u.normalVelocity()[i]*u.normalVelocity()[i]);
        }
        error(result.centerMomentumRates()[g],f*trace,std::max(1.L,std::abs(f*trace)),1e-12,"raw shared center momentum diagnostic");
    }
    // A typed physical-pair operator is independent of the historical
    // common-trace skew matrix. Recover its expected action from unique F.
    std::vector<long double> physicalK(raw.size());
    if(op.isImpermeable())for(std::size_t g=0;g<m.faceRates().size();++g){const auto& face=m.geometry().segments()[g];if(!inside(face))continue;
        physicalK[static_cast<std::size_t>(face.negativeVelocity)*n+face.positiveVelocity]+=m.faceRates()[g]/2;
        physicalK[static_cast<std::size_t>(face.positiveVelocity)*n+face.negativeVelocity]-=m.faceRates()[g]/2;}
    std::vector<long double> correction(raw.size());
    for(const auto& e:op.couplings()){
        check(e.first>=0 && e.first<e.second && e.second<n && axis(op.layout(),e.first)==axis(op.layout(),e.second),"invalid shared native pair");
        correction[static_cast<std::size_t>(e.first)*n+e.second]+=e.coefficient;
        correction[static_cast<std::size_t>(e.second)*n+e.first]-=e.coefficient;
    }
    const auto q=outflow(m);long double rawPower=0,desiredRawPower=0,grossRawPower=0;
    for(int i=0;i<n;++i){long double rawValue=0,rowScale=0;for(int j=0;j<n;++j){expected[i]+=.5L*(raw[static_cast<std::size_t>(i)*n+j]-raw[static_cast<std::size_t>(j)*n+i])*u.normalVelocity()[j];rawValue+=rawCf[static_cast<std::size_t>(i)*n+j]*u.normalVelocity()[j];rowScale+=.5L*(rawAbs[static_cast<std::size_t>(i)*n+j]+rawAbs[static_cast<std::size_t>(j)*n+i])*std::abs(u.normalVelocity()[j]);}
        for(int j=0;j<n;++j){const auto k=correction[static_cast<std::size_t>(i)*n+j];expected[i]+=k*u.normalVelocity()[j];rowScale+=std::abs(k*u.normalVelocity()[j]);}
        const int a=axis(op.layout(),i);const long double v=u.normalVelocity()[i];
        if(op.isImpermeable()){expected[i]=0;rowScale=0;for(int j=0;j<n;++j){const auto term=physicalK[static_cast<std::size_t>(i)*n+j]*u.normalVelocity()[j];expected[i]+=term;rowScale+=std::abs(term);}
            const auto diagonal=.5L*m.cells()[i].outwardIntegral/m.deltaTime()*v;expected[i]+=diagonal;rowScale+=std::abs(diagonal);}
        else {expected[i]+=.5L*qInternal[i]*v+boundary[i];rowScale+=std::abs(.5L*qInternal[i]*v)+std::abs(boundary[i]);}
        const auto load=-expected[i]+v*m.cells()[i].integratedDefect/m.deltaTime();
        peakLoad=std::max({peakLoad,error(result.momentumRates()[i],expected[i],std::max(1.L,rowScale),1e-12,"independent dense compatible load"),error(result.momentumLoads()[i],load,std::max(1.L,rowScale+std::abs(v*m.cells()[i].integratedDefect/m.deltaTime())),1e-12,"independent material momentum load")});
        momentum[a]+=result.momentumRates()[i];power[a]+=v*result.momentumLoads()[i];massPower[a]+=.5L*v*v*m.cells()[i].massChange/m.deltaTime();materialPower[a]+=.5L*v*v*m.cells()[i].integratedDefect/m.deltaTime();
        scale[a]+=std::abs(v)*rowScale+std::abs(v*load)+std::abs(.5L*v*v*m.cells()[i].massChange/m.deltaTime())+std::abs(v*v*m.cells()[i].integratedDefect/m.deltaTime());
        momentumScale[a]+=rowScale+std::abs(boundary[i]);rawPower+=v*rawValue;desiredRawPower+=.5L*v*v*qCf[i];grossRawPower+=std::abs(v*rawValue)+std::abs(.5L*v*v*qCf[i]);++loadChecks;
    }
    if(amr)amrRawPowerGap=std::max(amrRawPowerGap,std::abs(rawPower-desiredRawPower)/std::max(1.L,grossRawPower));
    for(int a=0;a<2;++a){peakMomentum=std::max(peakMomentum,error(momentum[a],bMomentum[a],std::max(1.L,momentumScale[a]),1e-12,"independent momentum boundary balance"));
        peakEnergy=std::max(peakEnergy,error(power[a]-massPower[a],materialPower[a]-bPower[a],std::max(1.L,scale[a]+std::abs(bPower[a])),1e-12,"independent material kinetic budget"));
        error(result.report().components[a].boundaryKineticPower,bPower[a],std::max(1.L,scale[a]),1e-12,"boundary KE report");error(result.report().components[a].materialKineticPower,materialPower[a],std::max(1.L,scale[a]),1e-12,"material KE report");}
    check(&result.material().faceRates()==&m.faceRates(),"result detached authoritative material");
    materialEnergyReference(m,u,result,t);++scenarios;return result;
}
void constant(const QuadtreeMesh& mesh,const HydroKepConvection& op,Vector2D constantVelocity,bool actual)
{
    const auto field=[constantVelocity](Vector2D){return constantVelocity;};const auto u=values(mesh,op.layout(),field);const auto t=traces(op.material(),field);const auto result=inspect(mesh,op,u,t);const auto q=outflow(op.material());
    std::vector<long double> gross(q.size());for(std::size_t g=0;g<op.material().faceRates().size();++g){const auto& segment=op.material().geometry().segments()[g];const auto f=std::abs(op.material().faceRates()[g]);
        if(segment.negativeVelocity>=0)gross[segment.negativeVelocity]+=f;if(segment.positiveVelocity>=0)gross[segment.positiveVelocity]+=f;for(const auto& w:op.weights()[g])gross[w.velocity]+=f*std::abs(w.weight);}
    for(const auto& e:op.couplings()){gross[e.first]+=std::abs(e.coefficient);gross[e.second]+=std::abs(e.coefficient);}
    for(int i=0;i<op.layout().velocityCount();++i){const auto v=static_cast<long double>(u.normalVelocity()[i]);error(result.momentumRates()[i],v*q[i],std::max(1.L,std::abs(v)*gross[i]),1e-12,"constant consistent convection");
        const auto recovered=(op.material().beginMass().masses()[i]*v+op.material().deltaTime()*result.momentumLoads()[i])/op.material().endMass().masses()[i];
        peakConstant=std::max(peakConstant,error(recovered,v,std::max(1.L,std::abs(v)),1e-10,"actual no-force constant material update"));
        if(actual){actualDefect=std::max(actualDefect,std::abs(op.material().cells()[i].integratedDefect));actualTransfer=std::max(actualTransfer,std::abs(op.material().cells()[i].geometryTransferChange));}}
}
void controlled(const QuadtreeMesh& mesh,const HydroKepConvection& center,bool quadratic=true,bool paired=false)
{
    const auto& m=center.material();const auto& l=m.layout();const int n=l.velocityCount();const bool rz=l.domain().axisymmetric();
    const auto op=paired?HydroKepConvection::buildControlled(mesh,m,center.periodicPairs(),{}):HydroKepConvection::buildControlled(mesh,m,{});
    const auto saved=m.faceIntegrals();
    const auto apply=[&](const HydroKepConvection& c,HydroNativeValues u,const std::function<Vector2D(Vector2D)>& field){
        auto v=u.normalVelocity();for(int i=0;i<n;++i)v[i]=v[c.representatives()[i]];u=HydroNativeValues(l,v,u.pressure());
        return paired?c.applyPaired(mesh,u,m.sources(),m.timeLayer()):c.apply(mesh,u,bindBoundary(mesh,m,traces(m,field)),m.sources(),m.timeLayer());};
    const auto constantField=[rz](Vector2D){return Vector2D{rz?0:.37,-.21};};
    const auto constantResult=apply(op,values(mesh,l,constantField),constantField);
    for(const auto& e:constantResult.dissipationSamples())check(e.sensor==0 && e.chi==0 && e.momentumRate==0,"constant field activated dissipation");
    const auto polynomialField=[rz,paired](Vector2D p){
        if(paired){const double t=rz?p.x:p.y;return Vector2D{rz?0:.2+.3*t+.4*t*t,.3-.2*t+.1*t*t};}
        return Vector2D{rz?p.x*(.2+.3*p.y):.2+.3*p.x+.4*p.y*p.y,.3-.2*p.x+.1*p.x*p.y+.2*p.y*p.y};};
    if(quadratic){const auto exact=apply(op,values(mesh,l,polynomialField),polynomialField);
        for(const auto& e:exact.dissipationSamples())if(e.chi!=0 || e.momentumRate!=0){
            const auto a=l.domain().face(mesh,l.nativeFace(e.first)).center,b=l.domain().face(mesh,l.nativeFace(e.second)).center;
            std::cout<<"QUADRATIC_SENSOR_WITNESS,"<<rz<<','<<paired<<','<<e.firstSegment<<','<<e.first<<','<<e.second<<','<<axis(l,e.first)<<','<<a.x<<','<<a.y<<','<<b.x<<','<<b.y<<','<<e.sensor<<','<<e.chi<<','<<e.halfJump<<std::endl;
            check(false,"resolved quadratic activated dissipation");}}
    const auto rough=[rz](Vector2D p){return Vector2D{(rz?p.x:1)*(p.y>.25 && p.y<.625?std::sin(71*p.x+53*p.y):0),.4+(p.x>.25 && p.x<.625?std::cos(61*p.x-47*p.y):0)};};
    auto u=values(mesh,l,rough);auto velocity=u.normalVelocity();for(int i=0;i<n;++i)velocity[i]=velocity[op.representatives()[i]];u=HydroNativeValues(l,velocity,u.pressure());
    const auto base=apply(center,u,rough),result=apply(op,u,rough);
    std::vector<long double> L(n),classMass(n),classLoad(n);std::array<long double,2> D{},gross{},kinetic{},massPower{},materialPower{},boundaryPower{},momentum{};
    std::set<int> seen;std::vector<int> partner(m.geometry().segments().size(),-1);
    for(const auto& p:op.periodicPairs()){partner[p.firstSegment]=p.secondSegment;partner[p.secondSegment]=p.firstSegment;}
    for(const auto& e:result.dissipationSamples()){
        check(seen.insert(e.firstSegment).second,"repeated stabilization connection");const auto& s=m.geometry().segments()[e.firstSegment];
        const int i=inside(s)?s.negativeVelocity:owner(s),j=inside(s)?s.positiveVelocity:owner(m.geometry().segments()[partner[e.firstSegment]]);
        check(e.first==i && e.second==j && (inside(s)?e.secondSegment==-1:e.secondSegment==partner[e.firstSegment]),"stabilization lost actual shared topology");
        const auto jump=static_cast<long double>(velocity[i])-velocity[j];
        const auto ramp=e.sensor<=.2L?0.L:e.sensor>=.5L?1.L:(e.sensor-.2L)/.3L;
        error(e.chi,ramp*ramp,1,1e-12,"sensor activation law");check(e.chi>=0 && e.chi<=1,"unbounded dissipation coefficient");
        const auto d=.5L*e.chi*std::abs(m.faceRates()[e.firstSegment]),f=d*jump,power=f*jump;L[i]+=f;L[j]-=f;const int a=axis(l,i);
        D[a]+=power;gross[a]+=std::abs(velocity[i]*f)+std::abs(velocity[j]*f)+power;
        controlledPeak=std::max(controlledPeak,error(e.momentumRate,f,std::max(1.L,std::abs(f)),1e-12,"independent stabilization segment"));
        if(e.chi>0 && m.faceRates()[e.firstSegment]!=0)++controlledActive;else ++controlledQuiet;
        if(inside(s)&&coarseFine(mesh,m.geometry(),s)&&power>0)++controlledCoarseFine;++controlledConnections;
    }
    for(std::size_t k=0;k<m.geometry().segments().size();++k){const auto& s=m.geometry().segments()[k];int i=-1,j=-1;
        if(inside(s)){i=s.negativeVelocity;j=s.positiveVelocity;}else if(partner[k]>static_cast<int>(k)){i=owner(s);j=owner(m.geometry().segments()[partner[k]]);}else continue;
        if(op.representatives()[i]!=op.representatives()[j])check(seen.count(static_cast<int>(k))==1,"missing actual stabilization connection");}
    for(int i=0;i<n;++i){const int r=op.representatives()[i];classMass[r]+=m.midpointMass().masses()[i];classLoad[r]+=L[i];}
    for(int i=0;i<n;++i){const int r=op.representatives()[i];L[i]=classLoad[r]*m.midpointMass().masses()[i]/classMass[r];const int a=axis(l,i);const long double v=velocity[i];
        const long double scale=std::max(1.L,std::abs(base.momentumLoads()[i])+std::abs(L[i]));
        controlledPeak=std::max({controlledPeak,error(result.dissipationRates()[i],L[i],scale,1e-12,"independent native stabilization"),
            error(result.momentumLoads()[i],base.momentumLoads()[i]-L[i],scale,1e-12,"actual controlled load"),
            error(result.momentumRates()[i],base.momentumRates()[i]+L[i],scale,1e-12,"actual controlled rate")});
        long double q=0;for(int side=0;side<4;++side)for(int k:m.geometry().volumes()[i].sides[side])q+=(side%2?1:-1)*m.faceRates()[k];
        const auto va=m.beginMass().phaseVolumes()[i],vb=m.endMass().phaseVolumes()[i];
        const auto mdot=(m.material().liquidDensity()*(vb.liquid-va.liquid)+m.material().gasDensity()*(vb.gas-va.gas))/m.deltaTime();
        kinetic[a]+=v*result.momentumLoads()[i];massPower[a]+=.5L*v*v*mdot;materialPower[a]+=.5L*v*v*(mdot+q);
        gross[a]+=std::abs(v*result.momentumLoads()[i])+v*v*(std::abs(mdot)+std::abs(q));momentum[a]+=L[i];controlledPeakDefect=std::max(controlledPeakDefect,std::abs(mdot+q));}
    for(std::size_t k=0;k<m.geometry().segments().size();++k){const auto& s=m.geometry().segments()[k];if(inside(s))continue;const int i=owner(s),a=axis(l,i);
        const long double v=velocity[i],f=(s.negativeVelocity>=0?1:-1)*m.faceRates()[k],b=result.boundary().traces()[k];const auto power=f*(v*b-.5L*v*v);boundaryPower[a]+=power;gross[a]+=std::abs(power);}
    for(int a=0;a<2;++a){check(D[a]>=0,"independent negative dissipation");const auto scale=std::max(1.L,gross[a]);
        controlledPeak=std::max({controlledPeak,error(result.report().components[a].dissipationPower,D[a],scale,1e-12,"independent dissipation diagnostic"),
            error(kinetic[a]-massPower[a],materialPower[a]-boundaryPower[a]-D[a],scale,1e-12,"actual material KE with dissipation"),
            error(momentum[a],0,std::max(1.L,D[a]),1e-12,"internal dissipative momentum conservation")});}
    // Scale/sign and admissible constant shifts do not change smoothness.
    const auto transformed=[&](Vector2D p){const auto v=rough(p);return Vector2D{-2*v.x,3-2*v.y};};
    const auto transformedResult=apply(op,values(mesh,l,transformed),transformed);
    check(result.dissipationSamples().size()==transformedResult.dissipationSamples().size(),"sensor topology depends on velocity");
    for(std::size_t k=0;k<result.dissipationSamples().size();++k)controlledPeak=std::max(controlledPeak,error(result.dissipationSamples()[k].chi,transformedResult.dissipationSamples()[k].chi,1,1e-10,"sensor shift/scale invariance"));
    check(saved==m.faceIntegrals() && &op.material().faceRates()==&m.faceRates(),"controlled operator modified material source");
    ++controlledCases;std::cout<<"CONTROLLED,"<<rz<<','<<paired<<','<<n<<','<<m.material().liquidDensity()<<','<<result.dissipationSamples().size()<<','<<D[0]<<','<<D[1]<<','<<controlledPeak<<std::endl;
}
void algebra(bool rz,int amr,int selection)
{
    // A one-column isolated RZ slab cannot identify an off-column physical
    // centroid. Its explicit rejection is tested below; two columns suffice.
    auto mesh=makeMesh(rz && selection==3?8:4,amr);Case s;s.pureFirstJet=true;s.radial=rz;s.planar=true;s.planeOffset=2;
    const auto geometry=HydroMomentumGeometry::build(*mesh,metrics(*mesh,rz,selection));const auto f=frames(*mesh,geometry,publishAt(*mesh,s,{}),publishAt(*mesh,s,{1,0,.04,Position::Midpoint}),publishAt(*mesh,s,{1,0,.04,Position::End}),1000);
    const auto advector=values(*mesh,geometry.layout(),[rz](Vector2D p){return Vector2D{(rz?p.x:1)*(.3+std::sin(3*p.y)),.1+std::cos(4*p.x)};});
    const auto op=HydroKepConvection::build(*mesh,material(*mesh,f,advector));nativeMoments(*mesh,op);affineTraces(*mesh,op,amr!=0);constant(*mesh,op,{rz?0:.37,-.21},false);
    const auto fa=[rz](Vector2D p){return Vector2D{(rz?p.x:1)*std::sin(5*p.x+3*p.y),std::cos(7*p.x-2*p.y)};};
    const auto fb=[rz](Vector2D p){return Vector2D{(rz?p.x:1)*(.2+.3*p.y),.4-.1*p.x};};const auto fs=[&](Vector2D p){return fa(p)+fb(p);};
    const auto a=inspect(*mesh,op,values(*mesh,op.layout(),fa),traces(op.material(),fa),amr!=0,!rz && !amr);
    const auto b=inspect(*mesh,op,values(*mesh,op.layout(),fb),traces(op.material(),fb));const auto sum=inspect(*mesh,op,values(*mesh,op.layout(),fs),traces(op.material(),fs));
    for(std::size_t i=0;i<a.momentumRates().size();++i)peakLinearity=std::max(peakLinearity,error(sum.momentumRates()[i],a.momentumRates()[i]+b.momentumRates()[i],std::max(1.L,std::abs(a.momentumRates()[i])+std::abs(b.momentumRates()[i])),1e-12,"compatible convection linearity"));
    const auto copied=op;check(&copied.weights()==&op.weights(),"prepared copy detached reconstruction");std::cout<<"ALGEBRA,"<<rz<<','<<amr<<','<<selection<<'\n';
    controlled(*mesh,op,selection==0);
}
void crossing(bool rz,int amr,bool planar)
{
    auto mesh=makeMesh(4,amr);Case s;s.pureFirstJet=true;s.radial=rz;s.planar=planar;s.radius=rz?.34:.19;s.motion=1;s.speed=.4;s.planeOffset=.49;
    const auto geometry=HydroMomentumGeometry::build(*mesh,metrics(*mesh,rz));const auto a=initial(*mesh,s);std::string stage;
    const auto mid=advanceFirstJet(*mesh,a,s,{1,0,.1,Position::Midpoint},stage),end=advanceFirstJet(*mesh,a,s,{1,0,.1,Position::End},stage);
    const Vector2D v=rz?Vector2D{0,.4}:Vector2D{.4,0};const auto advector=values(*mesh,geometry.layout(),[v](Vector2D){return v;});
    for(double rho:{1.,1000.,1e6}){const auto op=HydroKepConvection::build(*mesh,material(*mesh,frames(*mesh,geometry,a,mid,end,rho),advector));constant(*mesh,op,v,true);
        const auto probe=[rz](Vector2D p){return Vector2D{(rz?p.x:1)*(.2+std::sin(p.y)),.1+std::cos(2*p.x)};};
        inspect(*mesh,op,values(*mesh,op.layout(),probe),traces(op.material(),probe));++actualCrossings;
        controlled(*mesh,op);
        std::cout<<"CROSSING,"<<rz<<','<<amr<<','<<planar<<','<<rho<<'\n';}
}
// Independent line integration of the public actual source. GL8 differs
// from the production rule; every public piece knot is split before sampling.
// This helper reads neither production mean preparation nor moment reports.
struct NativeMeanReference
{
    std::array<long double,2> mean{},absoluteMean{};
    std::array<long double,3> quadratic{},absoluteQuadratic{}; // Physical xx/xy/yy advective response, before M/h^2.
    int pieces=0;
};
NativeMeanReference nativeMeanReference(const QuadtreeMesh& mesh,const HydroSpatialAdvector& actual,int u)
{
    const auto& l=actual.momentumGeometry().layout();check(!l.domain().axisymmetric(),"mean reference requires XY");
    const auto face=l.domain().face(mesh,l.nativeFace(u));const int normal=axis(l,u),tangent=1-normal;
    const auto& edge=face.horizontal?mesh.hEdges(face.edgeLabel):mesh.vEdges(face.edgeLabel);
    auto a=mesh.verts(edge.edgeToVert[0]).coord(),b=mesh.verts(edge.edgeToVert[1]).coord();
    const auto coordinate=[](Vector2D p,int c){return c==0?p.x:p.y;};
    if(coordinate(a,tangent)>coordinate(b,tangent))std::swap(a,b);
    const double lo=coordinate(a,tangent),hi=coordinate(b,tangent),plane=coordinate(a,normal);
    check(hi>lo && coordinate(b,normal)==plane,"mean reference native face is not canonical");
    const int native=face.negativeCell>=0?face.negativeCell:face.positiveCell,p=l.pressureDof(native);
    check(p>=0,"mean reference lacks actual pressure owner");std::vector<double> cuts{lo,hi};
    for(const auto& piece:actual.reconstruction().pieces())if(piece.pressure==p &&
        plane>=coordinate(piece.lower,normal) && plane<=coordinate(piece.upper,normal)){
        const double x=std::max(lo,coordinate(piece.lower,tangent)),y=std::min(hi,coordinate(piece.upper,tangent));
        if(x<y){cuts.push_back(x);cuts.push_back(y);}}
    std::sort(cuts.begin(),cuts.end());cuts.erase(std::unique(cuts.begin(),cuts.end()),cuts.end());
    constexpr std::array<long double,8> nodes{{-.96028985649753623168L,-.79666647741362673959L,-.52553240991632898582L,-.18343464249564980494L,
        .18343464249564980494L,.52553240991632898582L,.79666647741362673959L,.96028985649753623168L}};
    constexpr std::array<long double,8> weights{{.10122853629037625915L,.22238103445337447054L,.31370664587788728734L,.36268378337836198297L,
        .36268378337836198297L,.31370664587788728734L,.22238103445337447054L,.10122853629037625915L}};
    NativeMeanReference result;long double measured=0;
    for(std::size_t k=1;k<cuts.size();++k){const long double x=cuts[k-1],y=cuts[k],half=(y-x)/2,center=(x+y)/2;
        check(half>0,"empty independent source piece");++result.pieces;
        for(std::size_t q=0;q<nodes.size();++q){const double t=static_cast<double>(center+half*nodes[q]);
            const auto point=normal==0?Vector2D{plane,t}:Vector2D{t,plane};const auto v=actual.cellVelocityJet(p,point).jet.value;
            const long double w=half*weights[q];check(w>0,"independent mean has nonpositive physical weight");measured+=w;
            result.mean[0]+=w*v.x;result.mean[1]+=w*v.y;result.absoluteMean[0]+=w*std::abs(v.x);result.absoluteMean[1]+=w*std::abs(v.y);
            const long double dx=static_cast<long double>(point.x)-face.center.x,dy=static_cast<long double>(point.y)-face.center.y;
            const std::array<long double,3> q2{{2*v.x*dx,v.x*dy+v.y*dx,2*v.y*dy}};
            for(int c=0;c<3;++c){result.quadratic[c]+=w*q2[c];result.absoluteQuadratic[c]+=w*std::abs(q2[c]);}}}
    const long double length=static_cast<long double>(hi)-lo;
    error(measured,length,length,1e-12,"independent native line capacity");
    error(actual.momentumGeometry().metrics().faceMeasures()[u],length,length,1e-12,"independent true XY face measure");
    for(int c=0;c<2;++c){result.mean[c]/=length;result.absoluteMean[c]/=length;}
    for(int c=0;c<3;++c){result.quadratic[c]/=length;result.absoluteQuadratic[c]/=length;}return result;
}
// The translated slab keeps the original material-gross half-width and
// tolerance. Its centre is independently integrated; no target/reaction norm
// is added to the allowable residual. This also handles exact zero scale.
void nativeP2Slab(long double response,long double target,long double halfWidth,long double physicalGross,const char* why)
{
    check(std::isfinite(response)&&std::isfinite(target)&&std::isfinite(halfWidth)&&std::isfinite(physicalGross)
        &&halfWidth>=0&&physicalGross>=0,std::string(why)+": invalid physical slab data");
    check(std::abs(response-target)<=halfWidth+1e-12L*physicalGross,why);
}
void closedAnalyticMeans(const QuadtreeMesh& mesh,const HydroKepConvection& op,
    const std::vector<NativeMeanReference>& means,const std::vector<long double>& k,
    const std::vector<Vector2D>& points,const std::vector<bool>& fixed)
{
    // The source truth and its negative controls survive the method change.
    // They do not certify exact P1/P2 moments of the physical-pair operator.
    const auto& m=op.material();const auto& l=m.layout();const auto& actual=m.material().midpointRate().spatialAdvector();
    const int n=l.velocityCount();long double peak=0,pointMeanDifference=0,quadraticPeak=0;int jets=0,meanRows=0,quadraticReferences=0;
    const auto analytic=[](Vector2D p){const long double x=p.x,y=p.y;return std::array<long double,6>{{
        x*(1-x)*(1-2*y),-(1-2*x)*y*(1-y),(1-2*x)*(1-2*y),-2*x*(1-x),2*y*(1-y),-(1-2*x)*(1-2*y)}};};
    for(const auto& piece:actual.reconstruction().pieces())for(double tx:{1./6,.5,5./6})for(double ty:{1./6,.5,5./6}){
        const Vector2D p{piece.lower.x+(piece.upper.x-piece.lower.x)*tx,piece.lower.y+(piece.upper.y-piece.lower.y)*ty};
        const auto exact=analytic(p);const auto v=actual.cellVelocityJet(piece.pressure,p).jet;
        const std::array<long double,6> observed{{v.value.x,v.value.y,v.derivative.xx,v.derivative.xy,v.derivative.yx,v.derivative.yy}};
        for(int c=0;c<6;++c)peak=std::max(peak,error(observed[c],exact[c],c<2?.25L:1.L,1e-12,"independent closed P3 actual value/J"));++jets;}
    for(int i=0;i<n;++i){const auto face=l.domain().face(mesh,l.nativeFace(i));const int normal=axis(l,i);
        const long double x=face.center.x,y=face.center.y,h=m.geometry().metrics().faceMeasures()[i],variance=h*h/12;
        const std::array<long double,2> exact{{(x*(1-x)-(normal==1?variance:0))*(1-2*y),
            -(1-2*x)*(y*(1-y)-(normal==0?variance:0))}};
        const auto at=analytic(face.center);const auto point=actual.nativeVelocity(i);
        const long double factor=h*h/6;
        const std::array<long double,3> exactQuadratic=normal==0
            ?std::array<long double,3>{{0,-factor*x*(1-x),-factor*(1-2*x)*(1-2*y)}}
            :std::array<long double,3>{{factor*(1-2*x)*(1-2*y),factor*y*(1-y),0}};
        for(int c=0;c<3;++c){
            quadraticPeak=std::max(quadraticPeak,error(means[i].quadratic[c],exactQuadratic[c],factor,1e-12,"independent closed P3 analytic quadratic mean"));
            ++quadraticReferences;}
        for(int c=0;c<2;++c){peak=std::max(peak,error(means[i].mean[c],exact[c],.25L,1e-12,"independent closed P3 analytic face mean"));
            peak=std::max(peak,error(c==0?point.x:point.y,at[c],.25L,1e-12,"independent closed P3 analytic point"));}
        peak=std::max(peak,error(actual.values().normalVelocity()[i],exact[normal],.25L,1e-12,"analytic input true native face mean"));
        if(fixed[i]){check(actual.values().normalVelocity()[i]==0,"analytic closed normal wall changed");continue;}++meanRows;
        for(int c=0;c<3;++c)if(std::abs(exactQuadratic[c])>1e-6L*factor){
            const long double truthScale=means[i].absoluteQuadratic[c]+std::abs(exactQuadratic[c]);
            if(!nativeMeanP2Negative[0]){rejects([&]{error(means[i].quadratic[c],0,truthScale,1e-12,"P2 zero-point centre mutation");},"P2 zero-point centre mutation");++nativeMeanP2Negative[0];}
            if(!nativeMeanP2Negative[1]){rejects([&]{error(means[i].quadratic[c],-exactQuadratic[c],truthScale,1e-12,"P2 reversed response mutation");},"P2 reversed response mutation");++nativeMeanP2Negative[1];}
            if(c==1&&!nativeMeanP2Negative[2]){rejects([&]{error(means[i].quadratic[c],2*exactQuadratic[c],truthScale,1e-12,"P2 xy extra-two mutation");},"P2 xy extra-two mutation");++nativeMeanP2Negative[2];}}
        for(int c=0;c<2;++c){pointMeanDifference=std::max(pointMeanDifference,std::abs(means[i].mean[c]-at[c]));
            if(!nativeMeanP1OldPointNegative && std::abs(exact[c]-at[c])>1e-4L){
                const long double scale=means[i].absoluteMean[c]+std::abs(exact[c]);
                rejects([&]{error(means[i].mean[c],at[c],scale,1e-12,"old point-P1 target mutation");},"old point-P1 target mutation");
                ++nativeMeanP1OldPointNegative;}}}
    check(meanRows>0 && jets>0 && pointMeanDifference>1e-4L,"analytic closed group did not distinguish point from mean");
    const auto& report=op.momentReport();
    check(report.physicalPairing==1 && report.meanP1Rows==0 && report.meanP2Rows==0,
        "physical pair path claimed retired exact mean-P1/P2 certificates");
    // Recover the linear response of this K from the independently recovered
    // F-incidence matrix, rather than asserting the retired exact-P1 target.
    for(int i=0;i<n;++i)if(!fixed[i])for(int c=0;c<2;++c){long double observed=0,expected=0,gross=0;
        for(int j=0;j<n;++j){const long double delta=c==0?points[j].x-points[i].x:points[j].y-points[i].y;
            const auto term=k[static_cast<std::size_t>(i)*n+j]*delta;observed+=term;gross+=std::abs(term);}
        for(std::size_t g=0;g<m.faceRates().size();++g){const auto& segment=m.geometry().segments()[g];if(!inside(segment))continue;
            const int other=segment.negativeVelocity==i?segment.positiveVelocity:(segment.positiveVelocity==i?segment.negativeVelocity:-1);
            if(other<0)continue;const long double delta=c==0?points[other].x-points[i].x:points[other].y-points[i].y;
            const auto term=(segment.negativeVelocity==i?.5L:-.5L)*m.faceRates()[g]*delta;expected+=term;gross+=std::abs(term);}
        // The segment path and the aggregated-pair path have different
        // cancellation order. Include both paths' absolute work; a near-zero
        // final affine response is not its floating-point arithmetic scale.
        peak=std::max(peak,error(observed,expected,gross+std::abs(expected),1e-12,"analytic source physical-pair linear response"));}
    ++nativeMeanP1AnalyticCases;
    std::cout<<"SOURCE_MEAN_P1_ANALYTIC,"<<n<<','<<jets<<','<<meanRows<<','<<pointMeanDifference<<','<<peak<<",PASS\n";
    std::cout<<"SOURCE_MEAN_P2_ANALYTIC,"<<n<<','<<quadraticReferences<<','<<quadraticPeak<<",PASS\n";
    std::cout<<"P2_REFERENCE_NEGATIVE";for(auto count:nativeMeanP2Negative){check(count==1,"missing independent P2 reference negative control");std::cout<<','<<count;}std::cout<<",PASS\n";
}
void impermeableNative(const QuadtreeMesh&,const HydroKepConvection&,const HydroAdvectorBoundary&);
void closed(bool rz,int amr)
{
    auto mesh=makeMesh(4,amr);Case s;s.pureFirstJet=true;s.radial=rz;s.planar=true;s.planeOffset=2;
    const auto geometry=HydroMomentumGeometry::build(*mesh,metrics(*mesh,rz));const auto& l=geometry.layout();
    const auto stream=[rz](Vector2D p){const double radial=rz?p.x*p.x*(1-p.x)*(1-p.x):p.x*(1-p.x);return radial*p.y*(1-p.y);};
    std::vector<double> transported;
    for(int i=0;i<l.velocityCount();++i){const auto face=l.domain().face(*mesh,l.nativeFace(i));const auto& edge=face.horizontal?mesh->hEdges(face.edgeLabel):mesh->vEdges(face.edgeLabel);
        auto a=mesh->verts(edge.edgeToVert[0]).coord(),b=mesh->verts(edge.edgeToVert[1]).coord();if(face.horizontal?a.x>b.x:a.y>b.y)std::swap(a,b);
        const auto area=geometry.metrics().faceMeasures()[i];const double volumeFlux=(face.horizontal?-1:1)*(rz?static_cast<double>(2*pi):1)*(stream(b)-stream(a));
        transported.push_back(area>0?volumeFlux/area:0);}
    const auto f=frames(*mesh,geometry,publishAt(*mesh,s,{}),publishAt(*mesh,s,{1,0,.04,Position::Midpoint}),publishAt(*mesh,s,{1,0,.04,Position::End}),1);
    std::vector<int> walls;for(int u=0;u<l.velocityCount();++u)if(l.domain().face(*mesh,l.nativeFace(u)).outwardSign!=0)walls.push_back(u);
    const auto boundary=HydroAdvectorBoundary::impermeable(*mesh,geometry,walls);
    const auto m=material(*mesh,f,HydroNativeValues(l,transported,std::vector<double>(l.pressureCount())),&boundary);const auto op=HydroKepConvection::build(*mesh,m);
    const auto& actual=m.material().midpointRate().spatialAdvector();
    check(actual.boundary() && actual.boundary()->sameIdentity(boundary),"closed material lost explicit normal-law source");
    check(&m.faceRates()==&m.material().dualFaceRates(),"closed dual detached unique material F");
    for(int u:walls)for(int id:actual.reconstruction().primalAtoms()[u]){
        const auto& polynomial=actual.atomPolynomials()[id];for(auto c:polynomial.beginCoefficients)check(c==0,"bound normal polynomial nonzero");
        for(auto c:polynomial.endCoefficients)check(c==0,"bound end normal polynomial nonzero");
        const auto q=m.material().midpointRate().atomRates()[id];check(q.liquid==0 && q.gas==0 && m.material().atomRates()[id]==0,"bound phase/material atom nonzero");}
    for(std::size_t g=0;g<m.faceRates().size();++g)if(!inside(m.geometry().segments()[g]) && m.faceRates()[g]!=0){const auto& s=m.geometry().segments()[g];
        std::cout<<"CLOSED_WITNESS,"<<rz<<','<<amr<<','<<g<<','<<s.primalVelocity<<','<<s.normalAxis<<','<<s.begin.x<<','<<s.begin.y<<','<<s.end.x<<','<<s.end.y<<','<<s.area<<','<<m.faceRates()[g]<<','<<transported[s.primalVelocity]<<','<<m.material().faceRates()[s.primalVelocity]<<std::endl;
        check(false,"closed streamfunction produced boundary mass flux");}
    constant(*mesh,op,{rz?0:.37,-.21},false);
    const auto probe=[rz](Vector2D p){return Vector2D{(rz?p.x:1)*std::sin(5*p.x+3*p.y),std::cos(7*p.x-2*p.y)};};const auto u=values(*mesh,l,probe);const auto result=inspect(*mesh,op,u,std::vector<double>(m.faceRates().size()));
    long double power=0,gross=0;for(int i=0;i<l.velocityCount();++i){const auto value=u.normalVelocity()[i]*result.momentumRates()[i];power+=value;gross+=std::abs(value);}
    long double reference=0,referenceGross=0,qMax=0;
    for(int i=0;i<l.velocityCount();++i){
        check(m.beginMass().masses()[i]==m.midpointMass().masses()[i]&&m.midpointMass().masses()[i]==m.endMass().masses()[i],"static closed geometry mass changed");
        check(m.cells()[i].massChange==0,"static closed material mass changed");
        long double q=0;for(int side=0;side<4;++side)for(int id:m.geometry().volumes()[i].sides[side])q+=(side%2?1:-1)*m.faceRates()[id];
        const auto term=.5L*q*u.normalVelocity()[i]*u.normalVelocity()[i];reference+=term;referenceGross+=std::abs(term);qMax=std::max(qMax,std::abs(q));
    }
    const auto e=error(power,reference,std::max(1.L,gross+referenceGross),1e-12,"closed actual material convection power");
    if(!rz && amr==0)impermeableNative(*mesh,op,boundary);
    // Preserve the old zero-power number as a diagnostic, not an extra R=0 gate.
    std::cout<<"CLOSED,"<<rz<<','<<amr<<','<<std::abs(power)/std::max(1.L,gross)<<'\n';
    std::cout<<"CLOSED_MATERIAL_POWER,"<<rz<<','<<amr<<','<<power<<','<<reference<<','<<gross<<','<<referenceGross<<','<<qMax<<','<<e<<'\n';
}
template<class F>void reject(F f,const std::string& why){rejects(f,why);++rejected;}
void impermeableNative(const QuadtreeMesh& mesh,const HydroKepConvection& historical,const HydroAdvectorBoundary& law)
{
    const auto& m=historical.material();const auto& l=m.layout();const int n=l.velocityCount();
    const auto saved=m.faceRates();const auto source=m.sources();const auto clock=m.timeLayer();
    const auto op=HydroKepConvection::buildImpermeable(mesh,m,law);
    check(op.isImpermeable() && !historical.isImpermeable() && historical.fixedVelocityDofs().empty(),"typed closed path confused the historical operator");
    check(op.periodicPairs().empty() && op.fixedVelocityDofs()==law.velocityDofs(),"impermeable path manufactured periodic aliases or lost normal constraints");
    check(&op.material().faceRates()==&m.faceRates(),"impermeable preparation replaced unique actual F");
    std::vector<bool> fixed(n);std::vector<Vector2D> points;long double slip=0;
    const auto& actual=m.material().midpointRate().spatialAdvector();
    for(int i=0;i<n;++i){const auto face=l.domain().face(mesh,l.nativeFace(i));points.push_back(face.center);
        fixed[i]=face.outwardSign!=0;check(op.representatives()[i]==i,"impermeable native row was aliased");
        check(fixed[i]==std::binary_search(op.fixedVelocityDofs().begin(),op.fixedVelocityDofs().end(),i),"fixed set differs from actual native normal walls");
        if(fixed[i]){const auto b=actual.nativeVelocity(i);check((axis(l,i)==0?b.x:b.y)==0,"actual source violates normal impermeability");
            slip=std::max(slip,static_cast<long double>(std::abs(axis(l,i)==0?b.y:b.x)));}}
    check(slip>1e-3,"closed test does not distinguish slip from a no-slip wall");
    std::vector<NativeMeanReference> means;means.reserve(n);
    for(int i=0;i<n;++i){means.push_back(nativeMeanReference(mesh,actual,i));const int normal=axis(l,i);
        const long double input=actual.values().normalVelocity()[i],scale=means.back().absoluteMean[normal]+std::abs(input);
        check(scale>0 || means.back().mean[normal]==input,"normal mean discrepancy at exactly zero physical scale");
        error(means.back().mean[normal],input,scale,1e-12,"independent actual normal source mean/native U");}

    std::vector<long double> outward(n),gross(n);
    for(int i=0;i<n;++i)for(int side=0;side<4;++side)for(int g:m.geometry().volumes()[i].sides[side]){
        outward[i]+=(side%2?1:-1)*m.faceRates()[g];gross[i]+=std::abs(m.faceRates()[g]);}
    std::vector<long double> expectedK(static_cast<std::size_t>(n)*n);
    std::map<std::pair<int,int>,long double> physicalPairs;int physicalSegments=0;
    for(std::size_t g=0;g<saved.size();++g){const auto& segment=m.geometry().segments()[g];
        if(!inside(segment)){check(saved[g]==0,"typed source has actual external mass flux");continue;}++physicalSegments;
        const int i=segment.negativeVelocity,j=segment.positiveVelocity;const long double half=saved[g]/2;
        expectedK[static_cast<std::size_t>(i)*n+j]+=half;expectedK[static_cast<std::size_t>(j)*n+i]-=half;
        physicalPairs[{std::min(i,j),std::max(i,j)}]+=(i<j?half:-half);}
    int nonzeroPairs=0;for(const auto& pair:physicalPairs)if(pair.second!=0)++nonzeroPairs;
    check(op.momentReport().physicalPairing==1 && op.momentReport().physicalSegments==physicalSegments
        && op.momentReport().physicalPairs==nonzeroPairs,"physical segment/pair report does not match unique F incidence");
    const std::vector<double> boundaryTrace(saved.size()),pressure(l.pressureCount());const auto boundary=bindBoundary(mesh,m,boundaryTrace);
    // Recover the FINAL native operator by applying coordinate vectors. This
    // does not reuse the constraint assembler, couplings or reported moments.
    const auto frozen=op.frozenMomentumOperator(mesh,source,clock);
    check(frozen.rows()==n && frozen.cols()==n,"frozen momentum dimensions differ from native layout");
    long double frozenActionError=0,frozenCompatibilityError=0;
    std::vector<long double> k(static_cast<std::size_t>(n)*n);
    for(int j=0;j<n;++j){std::vector<double> unit(n);unit[j]=1;
        try {const auto r=op.apply(mesh,HydroNativeValues(l,unit,pressure),boundary,source,clock);
            for(int i=0;i<n;++i){k[static_cast<std::size_t>(i)*n+j]=r.momentumRates()[i]-(i==j?.5L*m.cells()[i].outwardIntegral/m.deltaTime():0);
                // Independent action oracle: the original apply, including its
                // actual u*r consistency impulse, not the matrix assembler.
                const long double h=frozen.coeff(i,j),expected=-r.momentumLoads()[i];
                const long double scale=std::abs(h)+std::abs(r.momentumRates()[i])
                    +(i==j?std::abs(m.cells()[i].integratedDefect/m.deltaTime()):0);
                check(scale>0 || h==expected,"frozen unit action differs at exact zero work");
                frozenActionError=std::max(frozenActionError,error(h,expected,scale,1e-12,"frozen unit action/original apply"));}
        } catch(const std::exception& failure) {
            const auto& cell=m.cells()[j];std::ostringstream context;context<<std::setprecision(17)
                <<"[impermeable unit basis column="<<j<<" massChange="<<cell.massChange
                <<" outwardIntegral="<<cell.outwardIntegral<<" outwardRate="<<cell.outwardRate
                <<" defectRate="<<cell.defectRate<<" independentFaceOutward="<<outward[j]
                <<" dt="<<m.deltaTime()<<"] "<<failure.what();throw std::runtime_error(context.str());
        }}
    long double peak=0,reactionNorm=0;int freeRows=0,tangentRows=0;
    const double length=std::max(mesh.maxCoord.x-mesh.minCoord.x,mesh.maxCoord.y-mesh.minCoord.y);
    const auto physicalError=[&](long double a,long double b,long double scale,const char* why){
        check(scale>0 || a==b,std::string(why)+": nonzero discrepancy at exact zero scale");return error(a,b,scale,1e-12,why);};
    for(int i=0;i<n;++i){const auto& box=m.geometry().volumes()[i];const double h=std::max(box.upper.x-box.lower.x,box.upper.y-box.lower.y);
        long double constantMoment=0,constantGross=0,norm=0;
        for(int j=0;j<n;++j){const auto coefficient=k[static_cast<std::size_t>(i)*n+j],reverse=k[static_cast<std::size_t>(j)*n+i];
            peak=std::max(peak,physicalError(coefficient,-reverse,std::abs(coefficient)+std::abs(reverse),"actual closed matrix skew"));
            if(axis(l,i)!=axis(l,j))peak=std::max(peak,physicalError(coefficient,0,gross[i]+gross[j],"closed matrix mixed transported components"));
            const auto expected=expectedK[static_cast<std::size_t>(i)*n+j];
            peak=std::max(peak,physicalError(coefficient,expected,gross[i]+gross[j]+std::abs(expected),"public unit-basis physical F pair"));
            constantMoment+=coefficient;constantGross+=std::abs(coefficient);norm+=std::abs(coefficient);}
        const auto canonical=m.cells()[i].outwardIntegral/m.deltaTime();
        peak=std::max(peak,physicalError(constantMoment,.5L*canonical,constantGross+std::abs(canonical),"actual all-row closed constant moment"));
        if(fixed[i]){reactionNorm=std::max(reactionNorm,norm);continue;}++freeRows;
        check(norm<=2*gross[i]+1e-12*gross[i],"actual free closed row violates original l1 bound");
        // Preserve the old translated-bound negative control as a reference
        // helper test only; the physical-pair K is not subject to that slab.
        if(!nativeMeanP2Negative[3]&&gross[i]>0){
            const long double target=m.midpointMass().masses()[i]*means[i].quadratic[0]/(h*h),halfWidth=gross[i]*h/length;
            nativeP2Slab(target+halfWidth,target,halfWidth,gross[i],"translated P2 upper boundary");
            nativeP2Slab(target-halfWidth,target,halfWidth,gross[i],"translated P2 lower boundary");
            rejects([&]{nativeP2Slab(target+halfWidth+gross[i],target,halfWidth,gross[i],"translated P2 overflow mutation");},"translated P2 overflow mutation");++nativeMeanP2Negative[3];}
        if(box.lower.x==mesh.minCoord.x || box.upper.x==mesh.maxCoord.x || box.lower.y==mesh.minCoord.y || box.upper.y==mesh.maxCoord.y)++tangentRows;
    }
    check(op.momentReport().boundedRows==freeRows && op.momentReport().coefficientRatio<=1+1e-12,"free stability certificate includes unbounded reaction rows");
    check(op.momentReport().meanP1Rows==0 && op.momentReport().meanP2Rows==0,"physical pairing claimed retired moment rows");
    peak=std::max(peak,physicalError(op.momentReport().maxReactionCoefficientL1,reactionNorm,reactionNorm,"independent finite reaction coefficient norm"));
    check(freeRows==n-static_cast<int>(law.velocityDofs().size()) && tangentRows>0,"physical-pair checks omitted wall-adjacent free tangential rows");
    std::cout<<"PHYSICAL_F_PAIR_MATRIX,"<<n<<','<<physicalSegments<<','<<nonzeroPairs<<','<<freeRows<<','<<tangentRows<<','<<peak<<",PASS\n";
    auto physical=values(mesh,l,[](Vector2D p){return Vector2D{std::sin(5*p.x+3*p.y),std::cos(7*p.x-2*p.y)};});
    auto velocity=physical.normalVelocity();for(int i:law.velocityDofs())velocity[i]=0;physical=HydroNativeValues(l,velocity,pressure);
    const auto r=inspect(mesh,op,physical,boundaryTrace);
    const Eigen::Map<const Eigen::VectorXd> vector(velocity.data(),n);
    const Eigen::VectorXd frozenAction=frozen*vector;
    long double maximumR=0,consistencyWork=0;
    for(int i=0;i<n;++i){long double actionGross=0;
        for(int j=0;j<n;++j)actionGross+=std::abs(static_cast<long double>(frozen.coeff(i,j))*velocity[j]);
        const auto consistency=velocity[i]*m.cells()[i].integratedDefect/m.deltaTime();
        frozenActionError=std::max(frozenActionError,physicalError(frozenAction[i],-r.momentumLoads()[i],
            actionGross+std::abs(r.momentumRates()[i])+std::abs(consistency),"frozen physical action/original apply"));
        maximumR=std::max(maximumR,std::abs(m.cells()[i].integratedDefect));consistencyWork+=std::abs(consistency);
        for(int j=0;j<n;++j){const long double forward=m.deltaTime()*static_cast<long double>(frozen.coeff(i,j));
            const long double reverse=m.deltaTime()*static_cast<long double>(frozen.coeff(j,i));
            const long double delta=i==j?static_cast<long double>(m.endMass().masses()[i])-m.beginMass().masses()[i]:0;
            const long double expected=i==j?-m.cells()[i].integratedDefect:0;
            const long double scale=std::abs(forward)+std::abs(reverse)+std::abs(expected)
                +(i==j?std::abs(static_cast<long double>(m.endMass().masses()[i]))+std::abs(static_cast<long double>(m.beginMass().masses()[i])):0);
            frozenCompatibilityError=std::max(frozenCompatibilityError,physicalError(forward+reverse+delta,expected,scale,
                "dt frozen symmetric part plus actual endpoint mass change"));}}
    check(maximumR>0 && consistencyWork>0,"frozen momentum witness erased actual nonzero material R");
    long double skewPower=0,powerGross=0;
    for(int i=0;i<n;++i){long double row=0,rowGross=0;for(int j=0;j<n;++j){const auto term=k[static_cast<std::size_t>(i)*n+j]*velocity[j];row+=term;rowGross+=std::abs(term);}
        skewPower+=velocity[i]*row;powerGross+=std::abs(velocity[i])*rowGross;}
    peak=std::max(peak,physicalError(skewPower,0,powerGross,"actual final closed skew kinetic power"));
    HydroConvectionDissipationConfig quiet;quiet.strength=0;
    const auto controlled=HydroKepConvection::buildImpermeableControlled(mesh,m,law,quiet);
    check(controlled.isImpermeable()&&controlled.fixedVelocityDofs()==law.velocityDofs(),"controlled branch lost typed impermeable source");
    const auto cr=controlled.apply(mesh,physical,boundary,source,clock);
    check(cr.momentumRates()==r.momentumRates()&&cr.momentumLoads()==r.momentumLoads(),"zero-strength controlled branch changed actual closed load");
    const auto quietFrozen=controlled.frozenMomentumOperator(mesh,source,clock);
    for(int i=0;i<n;++i)for(int j=0;j<n;++j)check(quietFrozen.coeff(i,j)==frozen.coeff(i,j),"zero-strength changed frozen operator");
    HydroConvectionDissipationConfig sensed;
    const auto nonlinear=HydroKepConvection::buildImpermeableControlled(mesh,m,law,sensed);
    reject([&]{nonlinear.frozenMomentumOperator(mesh,source,clock);},"frozen momentum rejects sensor dissipation");
    reject([&]{historical.frozenMomentumOperator(mesh,source,clock);},"frozen momentum requires typed impermeable source");
    auto partial=law.velocityDofs();partial.pop_back();const auto missing=HydroAdvectorBoundary::impermeable(mesh,m.geometry(),partial);
    reject([&]{HydroKepConvection::buildImpermeable(mesh,m,missing);},"impermeable");
    const auto foreign=HydroAdvectorBoundary::impermeable(mesh,HydroMomentumGeometry::build(mesh,m.geometry().metrics()),law.velocityDofs());
    reject([&]{HydroKepConvection::buildImpermeable(mesh,m,foreign);},"impermeable");
    reject([&]{op.applyPaired(mesh,physical,source,clock);},"paired apply");
    auto wrong=source;wrong.advector={};reject([&]{op.apply(mesh,physical,boundary,wrong,clock);},"different actual advector identity");
    auto wrongClock=clock;wrongClock.step+=1;reject([&]{op.apply(mesh,physical,boundary,source,wrongClock);},"different full time layer");
    reject([&]{op.frozenMomentumOperator(mesh,wrong,clock);},"different actual advector identity");
    reject([&]{op.frozenMomentumOperator(mesh,source,wrongClock);},"different full time layer");
    ++frozenOperatorCases;
    std::cout<<"FROZEN_MOMENTUM_OPERATOR,"<<n<<','<<frozen.nonZeros()<<','<<maximumR<<','<<consistencyWork<<','
        <<frozenActionError<<','<<frozenCompatibilityError<<",PASS\n";
    quiet.strength=-1;reject([&]{HydroKepConvection::buildImpermeableControlled(mesh,m,law,quiet);},"strength");
    closedAnalyticMeans(mesh,op,means,k,points,fixed);
    check(saved==m.faceRates() && &op.material().faceRates()==&m.faceRates(),"closed checks changed the authoritative F");
    ++impermeableCases;impermeableFreeRows+=freeRows;impermeableWallTangentRows+=tangentRows;
    std::cout<<"IMPERMEABLE_NATIVE,"<<n<<','<<law.velocityDofs().size()<<','<<freeRows<<','<<tangentRows<<','<<slip<<','<<peak<<",PASS\n";
}
// Independent continuum face functional for the smooth closed P3 source.
// It uses analytic value/J and physical line weights, never K or its report.
long double closedP3AdvectiveMean(Vector2D a,Vector2D b,int component)
{
    constexpr std::array<long double,3> nodes{{-.77459666924148337704L,0,.77459666924148337704L}};
    constexpr std::array<long double,3> weights{{5.L/9,8.L/9,5.L/9}};
    long double mean=0;
    for(int q=0;q<3;++q){const long double t=(1+nodes[q])/2,x=a.x+(static_cast<long double>(b.x)-a.x)*t,y=a.y+(static_cast<long double>(b.y)-a.y)*t;
        const long double ux=x*(1-x)*(1-2*y),uy=-(1-2*x)*y*(1-y);
        const long double ax=ux*(1-2*x)*(1-2*y)+uy*(-2*x*(1-x));
        const long double ay=ux*2*y*(1-y)+uy*(-(1-2*x)*(1-2*y));
        mean+=weights[q]*(component==0?ax:ay)/2;}
    return mean;
}
void physicalPairSpatial()
{
    std::array<long double,3> allError{},wallError{},allReference{},wallReference{};int family=0;
    for(int n:{4,8,16}){
        auto mesh=makeMesh(n);Case shape;shape.pureFirstJet=true;shape.planar=true;shape.planeOffset=2;
        const auto geometry=HydroMomentumGeometry::build(*mesh,metrics(*mesh,false));const auto& l=geometry.layout();
        const auto psi=[](Vector2D p)->long double{const long double x=p.x,y=p.y;return x*(1-x)*y*(1-y);};
        std::vector<double> normal(l.velocityCount());std::vector<int> walls;
        for(int u=0;u<l.velocityCount();++u){const auto face=l.domain().face(*mesh,l.nativeFace(u));const auto& edge=face.horizontal?mesh->hEdges(face.edgeLabel):mesh->vEdges(face.edgeLabel);
            auto a=mesh->verts(edge.edgeToVert[0]).coord(),b=mesh->verts(edge.edgeToVert[1]).coord();if(face.horizontal?a.x>b.x:a.y>b.y)std::swap(a,b);
            const long double measure=geometry.metrics().faceMeasures()[u];check(measure>0,"physical P3 native face has no measure");
            normal[u]=static_cast<double>((face.horizontal?-1.L:1.L)*(psi(b)-psi(a))/measure);
            if(face.outwardSign!=0){walls.push_back(u);check(normal[u]==0,"physical P3 input violates exact normal wall");}}
        const auto law=HydroAdvectorBoundary::impermeable(*mesh,geometry,walls);
        const auto f=frames(*mesh,geometry,publishAt(*mesh,shape,{}),publishAt(*mesh,shape,{1,0,.04,Position::Midpoint}),publishAt(*mesh,shape,{1,0,.04,Position::End}),1);
        const HydroNativeValues native(l,normal,std::vector<double>(l.pressureCount()));const auto m=material(*mesh,f,native,&law);const auto saved=m.faceRates();
        const auto op=HydroKepConvection::buildImpermeable(*mesh,m,law);const auto boundary=bindBoundary(*mesh,m,std::vector<double>(saved.size()));
        const auto result=op.apply(*mesh,native,boundary,m.sources(),m.timeLayer());const auto& actual=m.material().midpointRate().spatialAdvector();
        check(op.momentReport().physicalPairing==1 && op.momentReport().meanP1Rows==0 && op.momentReport().meanP2Rows==0,"spatial family did not use the physical-pair method");
        long double sourcePeak=0,allSquared=0,wallSquared=0,referenceSquared=0,wallReferenceSquared=0;int freeRows=0,wallRows=0;
        const auto& pieces=actual.reconstruction().pieces();
        for(std::size_t k=0;k<pieces.size();++k){const auto& piece=pieces[k];const auto p=piece.lower+(piece.upper-piece.lower)*.5;const long double x=p.x,y=p.y;
            const auto j=actual.pieceVelocityJet(static_cast<int>(k),p).jet;
            const std::array<long double,6> got{{j.value.x,j.value.y,j.derivative.xx,j.derivative.xy,j.derivative.yx,j.derivative.yy}};
            const std::array<long double,6> truth{{x*(1-x)*(1-2*y),-(1-2*x)*y*(1-y),(1-2*x)*(1-2*y),-2*x*(1-x),2*y*(1-y),-(1-2*x)*(1-2*y)}};
            for(int c=0;c<6;++c)sourcePeak=std::max(sourcePeak,error(got[c],truth[c],c<2?.25L:1.L,1e-12,"spatial family actual P3 source value/J"));}
        for(int u=0;u<l.velocityCount();++u){if(std::binary_search(walls.begin(),walls.end(),u))continue;++freeRows;
            const auto face=l.domain().face(*mesh,l.nativeFace(u));const auto& edge=face.horizontal?mesh->hEdges(face.edgeLabel):mesh->vEdges(face.edgeLabel);
            const auto a=mesh->verts(edge.edgeToVert[0]).coord(),b=mesh->verts(edge.edgeToVert[1]).coord();const auto& box=m.geometry().volumes()[u];
            const long double mass=m.midpointMass().masses()[u],outward=m.cells()[u].outwardIntegral/m.deltaTime();
            check(mass>0,"physical P3 dual mass is not positive");
            const long double response=(result.momentumRates()[u]-normal[u]*outward)/mass,truth=closedP3AdvectiveMean(a,b,axis(l,u)),difference=response-truth;
            allSquared+=mass*difference*difference;referenceSquared+=mass*truth*truth;
            if(box.lower.x==mesh->minCoord.x || box.upper.x==mesh->maxCoord.x || box.lower.y==mesh->minCoord.y || box.upper.y==mesh->maxCoord.y){
                ++wallRows;wallSquared+=mass*difference*difference;wallReferenceSquared+=mass*truth*truth;}}
        check(freeRows==2*n*(n-1) && wallRows==4*(n-1),"physical P3 family omitted free or wall-adjacent rows");
        check(saved==m.faceRates() && &op.material().faceRates()==&m.faceRates(),"physical spatial check replaced unique F");
        allError[family]=std::sqrt(allSquared);wallError[family]=std::sqrt(wallSquared);
        allReference[family]=std::sqrt(referenceSquared);wallReference[family]=std::sqrt(wallReferenceSquared);
        check(allReference[family]>0 && wallReference[family]>0,"physical P3 family has no independent reference scale");
        std::cout<<"PHYSICAL_PAIR_P3_SPACE,"<<n<<','<<freeRows<<','<<wallRows<<','<<allError[family]<<','<<wallError[family]<<','<<allReference[family]<<','<<wallReference[family]<<','<<sourcePeak<<'\n';
        physicalPairSpatialRows+=freeRows;physicalPairWallRows+=wallRows;++family;}
    check(family==3,"physical P3 family incomplete");
    for(int i=1;i<3;++i){check(allError[i-1]>1e-11L*allReference[i-1] && allError[i]>1e-11L*allReference[i]
            && wallError[i-1]>1e-11L*wallReference[i-1] && wallError[i]>1e-11L*wallReference[i],"physical P3 spatial error is unresolved");
        const long double allOrder=std::log2(allError[i-1]/allError[i]),wallOrder=std::log2(wallError[i-1]/wallError[i]);
        std::cout<<"PHYSICAL_PAIR_P3_ORDER,"<<(i==1?4:8)<<','<<(i==1?8:16)<<','<<allOrder<<','<<wallOrder<<'\n';
        check(allError[i]<allError[i-1] && wallError[i]<wallError[i-1],"physical-pair native face-mean advection is not monotone");
        // As in the original three-grid family, retain both raw orders and
        // apply the unchanged 1.8 gate to the finest resolved pair.
        if(i==2)check(allOrder>=1.8L && wallOrder>=1.8L,"physical-pair native face-mean advection fails original 1.8 spatial gate");}
    ++physicalPairSpatialFamilies;std::cout<<"PHYSICAL_PAIR_P3_SPACE_RESULT,PASS\n";
}
// Independent exact values are powers of two or integer multiples of the
// format quantum. This selector does not replace any legacy range negative.
void finiteRationalOracleOnly()
{
    using Exact=LAMP::Hydrodynamics::Detail::HydroFiniteContraction;
    using Integer=boost::multiprecision::cpp_int;
    using Rational=boost::multiprecision::cpp_rational;
    struct RestoreEnvironment {
        std::fenv_t state{};int error=errno;
#if defined(_M_X64) || defined(__SSE2__)
        unsigned csr=_mm_getcsr();
#endif
        RestoreEnvironment(){check(std::fegetenv(&state)==0,"rational oracle cannot capture native environment");}
        ~RestoreEnvironment(){if(std::fesetenv(&state)!=0)std::terminate();
#if defined(_M_X64) || defined(__SSE2__)
            _mm_setcsr(csr);
#endif
            errno=error;}
    } restoreEnvironment;
    unsigned roundedChecks=0,comparisonChecks=0,identityChecks=0,negativeChecks=0;
    // This conversion is exact in Boost 1.88 rational_adaptor::operator=(Float):
    // frexp/ldexp form an integral significand, then integer shifts and gcd.
    // It never rounds a rational reference back to a floating-point number.
    const auto native=[](long double value){check(std::isfinite(value),"oracle received nonfinite native value");return Rational(value);};
    const auto sign=[](const Rational& value){return value<0?-1:value>0?1:0;};
    const auto certificate=[&](long double actual,const Rational& truth,const std::string& label){
        check(std::isfinite(actual),label+": finite RN result required");
        if(truth==0)check(actual==0&&!std::signbit(actual),label+": exact zero must publish positive zero");
        else check(std::signbit(actual)==(truth<0),label+": rounded value/zero has wrong sign");
        const Rational value=native(actual);
        Rational lower,upper;
        const auto maximum=std::numeric_limits<long double>::max();
        if(actual==-maximum){upper=native(std::nextafter(actual,0.L));lower=value-(upper-value);}
        else if(actual==maximum){lower=native(std::nextafter(actual,0.L));upper=value+(value-lower);}
        else {lower=native(std::nextafter(actual,-std::numeric_limits<long double>::infinity()));
            upper=native(std::nextafter(actual,std::numeric_limits<long double>::infinity()));}
        // At +/-max the extrapolated neighbour is the exact unrepresentable
        // power of two. An overflow tie must not certify the odd max value.
        const Rational lowMid=(lower+value)/2,highMid=(value+upper)/2;
        check(truth>=lowMid&&truth<=highMid,label+": exact truth is outside the RN neighbour interval");
        const auto even=[&](const Rational& gap){const Rational units=value/gap;
            check(boost::multiprecision::denominator(units)==1,label+": neighbour quantum is not integral");
            Integer bits=boost::multiprecision::numerator(units);if(bits<0)bits=-bits;return (bits&1)==0;};
        if(truth==lowMid)check(even(value-lower),label+": lower halfway chose an odd significand");
        if(truth==highMid)check(even(upper-value),label+": upper halfway chose an odd significand");
        ++roundedChecks;
    };
    struct Witness {
        Exact value;Rational truth=0;
        void term(std::initializer_list<long double> factors,long double divisor=1.L){
            value.addQuotientProduct(factors,divisor);Rational term=1;
            for(const auto factor:factors)term*=Rational(factor);
            term/=Rational(divisor);truth+=term;
        }
    };
    const auto certify=[&](const Witness& witness,const std::string& label){const auto actual=witness.value.rounded();
        certificate(actual,witness.truth,label);Exact roundedNative;roundedNative.add(actual);
        check(witness.value.compare(roundedNative)==sign(witness.truth-native(actual)),
            label+": exact compare lost the rational remainder around RN");++comparisonChecks;};
    const auto equal=[&](const Witness& left,const Witness& right,const std::string& label){
        check(left.value.compare(right.value)==sign(left.truth-right.truth),label+": exact rational comparison differs");
        check(right.value.compare(left.value)==sign(right.truth-left.truth),label+": reverse exact comparison differs");
        comparisonChecks+=2;
    };
    const auto identical=[&](const Witness& left,const Witness& right,const std::string& label){
        check(left.truth==right.truth,label+": oracle identity fixture is incorrect");equal(left,right,label);++identityChecks;
    };
    const auto expectInvalid=[&](const char* message,const auto& operation){bool rejectedHere=false;
        try{operation();}catch(const std::invalid_argument& e){rejectedHere=true;
            check(std::string(e.what())==std::string("[HydroFiniteContraction] ")+message,
                std::string("rational oracle changed first error: ")+e.what());}
        check(rejectedHere,std::string("rational oracle missing rejection: ")+message);++negativeChecks;
    };
    const auto q=std::numeric_limits<long double>::denorm_min();
    const auto epsilon=std::numeric_limits<long double>::epsilon();
    const auto huge=std::scalbn(1.L,std::numeric_limits<long double>::max_exponent/2);
    check(q>0&&std::isfinite(huge),"rational oracle needs native binary range witnesses");
    Witness one;one.term({1.L});
    Witness zero;
    Witness dyadic;dyadic.term({3.L,std::scalbn(1.L,-17),std::scalbn(1.L,20)});
    dyadic.term({-5.L,4.L});dyadic.term({-4.L});
    identical(dyadic,zero,"dyadic exact cancellation");certify(dyadic,"dyadic exact cancellation");
    // 1/3+1/4+5/12-2/24=11/12; the large common-odd terms
    // differ by binary exponent only and cancel without primitive overflow.
    Witness sameOdd;sameOdd.term({1.L},3.L);sameOdd.term({.25L});sameOdd.term({5.L},12.L);
    sameOdd.term({-2.L},24.L);sameOdd.term({huge},3.L);sameOdd.term({-huge},3.L);
    Witness elevenTwelfths;elevenTwelfths.term({11.L},12.L);
    identical(sameOdd,elevenTwelfths,"same odd denominator mixed with D1");certify(sameOdd,"same odd denominator mixed with D1");
    Witness negativeDivisor;negativeDivisor.term({1.L},-6.L);negativeDivisor.term({-3.L},-12.L);
    negativeDivisor.term({.125L});negativeDivisor.term({.5L},-3.L);
    Witness oneTwentyFourth;oneTwentyFourth.term({1.L},24.L);
    identical(negativeDivisor,oneTwentyFourth,"negative native divisors");certify(negativeDivisor,"negative native divisors");
    Witness incompatible;incompatible.term({1.L},3.L);incompatible.term({1.L},5.L);incompatible.term({7.L},15.L);
    incompatible.term({.25L});incompatible.term({-.25L});
    identical(incompatible,one,"incompatible odd denominators 3 and 5");certify(incompatible,"incompatible odd denominators 3 and 5");
    Witness primeCancellation;primeCancellation.term({1.L});
    for(const long double divisor:{3.L,5.L,7.L,11.L,257.L,65537.L}){
        primeCancellation.term({1.L},divisor);certify(primeCancellation,"prime reciprocal before cancellation");
        primeCancellation.term({-1.L},divisor);identical(primeCancellation,one,"prime reciprocal cancellation");
    }
    certify(primeCancellation,"prime cancellation final");
    Witness wideSpan;wideSpan.term({huge},3.L);wideSpan.term({-huge},3.L);
    wideSpan.term({q},3.L);wideSpan.term({q},3.L);wideSpan.term({q},3.L);
    Witness quantum;quantum.term({q});identical(wideSpan,quantum,"large exponent cancellation with rational subquantum terms");
    certify(wideSpan,"large exponent cancellation with rational subquantum terms");
    // Compare values that intentionally have the same native RN result.
    Witness justAboveOne;justAboveOne.term({1.L});justAboveOne.term({q},3.L);
    Witness justBelowOne;justBelowOne.term({1.L});justBelowOne.term({-q},3.L);
    certify(justAboveOne,"exact compare above identical RN");certify(justBelowOne,"exact compare below identical RN");
    equal(justAboveOne,one,"positive sub-ULP exact compare");equal(justBelowOne,one,"negative sub-ULP exact compare");
    equal(justAboveOne,justBelowOne,"opposite exact remainders with identical RN");
    Witness original;original.term({-7.L},15.L);Witness copy(original);
    copy.value.add(copy.value);copy.truth*=2;certify(copy,"self add snapshot");
    copy.value.addAbsolute(copy.value);copy.truth=0;identical(copy,zero,"negative abs-self cancellation");certify(copy,"negative abs-self cancellation");
    copy=original;copy.value=copy.value;certify(copy,"self copy assignment");
    copy.value.add(copy.value,-1);copy.truth=0;identical(copy,zero,"self subtract snapshot");certify(copy,"self subtract snapshot");
    Witness positive;positive.term({7.L},15.L);positive.value.addAbsolute(positive.value);positive.truth*=2;
    Witness fourteenFifteenths;fourteenFifteenths.term({14.L},15.L);
    identical(positive,fourteenFifteenths,"positive abs-self doubles");certify(positive,"positive abs-self doubles");
    Witness independent(original);original.value.reset();original.truth=0;
    certify(original,"reset to exact positive zero");certify(independent,"copy remains independent after reset");
    Witness assigned;assigned=independent;assigned.value.reset();assigned.truth=0;
    certify(assigned,"copy-assigned reset is independent");certify(independent,"copy source survives assigned reset");
    Witness moved(std::move(independent));certify(moved,"move construction retains exact rational");
    Witness moveAssigned;moveAssigned=std::move(moved);certify(moveAssigned,"move assignment retains exact rational");
    // The query denominator is an accumulated rational, not a native dt.
    Witness denominator;denominator.term({1.L},3.L);denominator.term({1.L},5.L);
    Witness numerator;numerator.term({7.L},11.L);numerator.term({-.125L});
    certificate(numerator.value.scaledQuotient(.25L,denominator.value),
        numerator.truth*Rational(1,4)/denominator.truth,"scaled quotient of accumulated odd-denominator rationals");
    denominator.value.add(denominator.value,-1);denominator.truth=0;
    denominator.term({-1.L},3.L);denominator.term({-1.L},5.L);
    certificate(numerator.value.scaledQuotient(-.25L,denominator.value),
        numerator.truth*Rational(-1,4)/denominator.truth,"scaled quotient negative accumulated divisor");
    Witness largeOdd;largeOdd.term({std::scalbn(1.L,70)});largeOdd.term({1.L});
    check(boost::multiprecision::denominator(largeOdd.truth)==1
        &&boost::multiprecision::msb(boost::multiprecision::numerator(largeOdd.truth))>=64,
        "scaled quotient lacks a beyond-word accumulated odd divisor");
    certificate(numerator.value.scaledQuotient(1.5L,largeOdd.value),
        numerator.truth*Rational(3,2)/largeOdd.truth,"scaled quotient with beyond-word odd accumulated divisor");
    Witness arbitraryRational;
    arbitraryRational.term({std::scalbn(1.L,70)});arbitraryRational.term({1.L},3.L);arbitraryRational.term({q},5.L);
    check(boost::multiprecision::msb(boost::multiprecision::numerator(arbitraryRational.truth))>=64
        &&boost::multiprecision::msb(boost::multiprecision::denominator(arbitraryRational.truth))>=64,
        "scaled quotient lacks large accumulated rational numerator and denominator");
    certificate(numerator.value.scaledQuotient(.5L,arbitraryRational.value),
        numerator.truth*Rational(1,2)/arbitraryRational.truth,"scaled quotient with general large accumulated rational divisor");
    certificate(numerator.value.scaledQuotient(-0.L,largeOdd.value),Rational(0),"scaled quotient exact zero scale");
    certify(numerator,"scaled quotient leaves numerator unchanged");certify(largeOdd,"scaled quotient leaves denominator unchanged");
    const auto midpointCase=[&](long double base,std::initializer_list<long double> factors,const char* label){
        Witness value;value.term({base});value.term(factors);certify(value,label);};
    midpointCase(1.L,{.5L,epsilon},"normal even lower halfway");
    midpointCase(1.L,{1.5L,epsilon},"normal odd lower halfway rounds upward");
    midpointCase(1.L,{-.25L,epsilon},"normal binade-boundary halfway from below");
    midpointCase(-1.L,{-.5L,epsilon},"negative normal halfway");
    Witness halfwayBelow;halfwayBelow.term({1.L});halfwayBelow.term({.5L,epsilon});halfwayBelow.term({-q});
    Witness halfwayAbove;halfwayAbove.term({1.L});halfwayAbove.term({.5L,epsilon});halfwayAbove.term({q});
    certify(halfwayBelow,"normal halfway minus one exact minimum quantum");
    certify(halfwayAbove,"normal halfway plus one exact minimum quantum");
    const auto minimum=std::numeric_limits<long double>::min();
    midpointCase(minimum,{-.5L,q},"normal-subnormal boundary even halfway");
    midpointCase(minimum,{-1.5L,q},"normal-subnormal boundary odd halfway");
    for(const auto coefficient:{.25L,.5L,.75L,1.5L,2.5L,-.25L,-.5L,-.75L,-1.5L,-2.5L}){
        Witness value;value.term({coefficient,q});certify(value,"subnormal/zero signed neighbour certificate");}
    Witness subquantum;subquantum.term({q},3.L);certify(subquantum,"rational positive subquantum zero");
    Witness negativeSubquantum;negativeSubquantum.term({-q},3.L);certify(negativeSubquantum,"rational negative subquantum zero");
    const auto maximum=std::numeric_limits<long double>::max();
    const auto topSpacing=std::scalbn(1.L,std::numeric_limits<long double>::max_exponent-std::numeric_limits<long double>::digits);
    midpointCase(maximum,{.25L,topSpacing},"finite maximum below overflow halfway");
    midpointCase(-maximum,{-.25L,topSpacing},"negative finite maximum below overflow halfway");
    expectInvalid("final contraction overflow",[&]{Exact value;value.add(maximum);value.addProduct({.5L,topSpacing});(void)value.rounded();});
    expectInvalid("final contraction overflow",[&]{Exact value;value.add(-maximum);value.addProduct({-.5L,topSpacing});(void)value.rounded();});
    const auto nan=std::numeric_limits<long double>::quiet_NaN(),infinity=std::numeric_limits<long double>::infinity();
    Exact errors(sameOdd.value),before(errors),zeroDenominator;
    expectInvalid("invalid bounded contraction arity",[&]{errors.addQuotientProduct({},0.L);});
    expectInvalid("invalid bounded contraction arity",[&]{errors.addQuotientProduct({0.L,1.L,1.L,1.L,1.L,1.L,1.L,1.L,nan},0.L);});
    expectInvalid("zero contraction denominator",[&]{errors.addQuotientProduct({nan},0.L);});
    expectInvalid("zero contraction denominator",[&]{errors.addQuotientProduct({0.L},-0.L);});
    expectInvalid("nonfinite input factor",[&]{errors.addQuotientProduct({0.L,nan},3.L);});
    expectInvalid("nonfinite input factor",[&]{errors.addQuotientProduct({0.L,infinity},3.L);});
    expectInvalid("nonfinite input factor",[&]{errors.addQuotientProduct({0.L},nan);});
    expectInvalid("invalid contraction sum sign",[&]{errors.add(errors,0);});
    expectInvalid("zero exact contraction denominator",[&]{(void)errors.scaledQuotient(nan,zeroDenominator);});
    expectInvalid("zero exact contraction denominator",[&]{(void)errors.scaledQuotient(0.L,zeroDenominator);});
    expectInvalid("nonfinite input factor",[&]{(void)errors.scaledQuotient(nan,one.value);});
    check(errors.compare(before)==0,"rejected rational operation mutated the retained exact value");++comparisonChecks;
    certificate(errors.rounded(),sameOdd.truth,"rejected rational operations preserve original value");
    std::cout<<"FINITE_CONTRACTION_RATIONAL_ORACLE_PASS,rounded_certificates,"<<roundedChecks
        <<",exact_comparisons,"<<comparisonChecks<<",algebraic_identities,"<<identityChecks
        <<",exact_error_checks,"<<negativeChecks<<",boost_rational_to_float_conversion,0\n";
}

void finiteContractionOnly()
{
    using Exact=LAMP::Hydrodynamics::Detail::HydroFiniteContraction;
    const auto q=std::numeric_limits<long double>::denorm_min();
    check(q>0,"finite contraction contract requires gradual underflow");
    const auto product=[](std::initializer_list<long double> factors){Exact x;x.addProduct(factors);return x.rounded();};
    const auto pzero=product({.5L,q}),nzero=product({-.5L,q});
    check(pzero==0&&!std::signbit(pzero)&&nzero==0&&std::signbit(nzero),"exact half-quantum signed ties-to-even zero");
    Exact trace;trace.addProduct({.5L,q});trace.addProduct({.5L,q});
    check(trace.rounded()==q,"complete trace lost two half-quantum terms");
    check(product({1.5L,q})==2*q&&product({2.5L,q})==2*q,"subnormal ties-to-even integer parity");
    Exact rate;rate.addQuotientProduct({.25L,q},.125L);
    check(rate.rounded()==2*q,"complete u*R/dt lost representable consistency rate");
    check(product({4.L,.5L,q})==2*q,"complete F*w*U lost representable centre flux");
    Exact row;row.addProduct({q,.5L});row.addProduct({q,.5L});
    check(row.rounded()==q,"complete row rounded each product early");
    const int exponent=std::numeric_limits<long double>::max_exponent/2;
    const auto u=std::scalbn(1.L,exponent),small=std::scalbn(1.L,-exponent);
    check(std::isfinite(u)&&small>0&&product({u,u,small})==u,"large-u/small-F contraction formed u squared early");
    const auto maximum=std::numeric_limits<long double>::max();
    Exact cancellation;cancellation.addProduct({maximum,2.L});cancellation.add(-maximum);
    check(cancellation.rounded()==maximum,"complete signed contraction overflowed before cancellation");
    Exact thirds;thirds.addQuotientProduct({1.L},3.L);
    check(thirds.rounded()==1.L/3.L,"exact rational quotient rounded incorrectly");
    reject([&]{Exact overflow;overflow.addProduct({maximum,2.L});(void)overflow.rounded();},"final contraction overflow");
    // reset retains exact rational semantics; references below are original
    // thirds, powers of two, or integer multiples of the format quantum.
    static_assert(!noexcept(std::declval<Exact&>().reset()),"reset must retain its environment failure path");
    {
        struct NativeEnvironment {
            std::fenv_t state{};int error=errno;
#if defined(_M_X64) || defined(__SSE2__)
            unsigned csr=_mm_getcsr();
#endif
            NativeEnvironment(){check(std::fegetenv(&state)==0,"reset test cannot capture native environment");}
            bool restore() const noexcept {
                const bool ok=std::fesetenv(&state)==0;
#if defined(_M_X64) || defined(__SSE2__)
                _mm_setcsr(csr);
#endif
                errno=error;return ok;
            }
            ~NativeEnvironment(){if(!restore())std::terminate();}
            bool same(const NativeEnvironment& other) const {
                return error==other.error&&std::memcmp(&state,&other.state,sizeof state)==0
#if defined(_M_X64) || defined(__SSE2__)
                    &&csr==other.csr
#endif
                    ;
            }
        } entry;
        unsigned values=0,environmentPairs=0,ftzPairs=0,dazPairs=0;bool ftzRejected=false,dazRejected=false;
        const auto positiveZero=[&](const Exact& value){const auto z=value.rounded();
            check(z==0&&!std::signbit(z),"reset did not restore positive exact zero");++values;};
        Exact reused;reused.add(8.L);const Exact copied(reused);
        reused.reset();positiveZero(reused);reused.addProduct({.5L,q});reused.addProduct({.5L,q});
        check(reused.rounded()==q&&copied.rounded()==8.L,"reset/reuse changed an independent copy");
        Exact assigned;assigned=copied;assigned.reset();positiveZero(assigned);
        check(copied.rounded()==8.L&&reused.rounded()==q,"copy-assigned reset aliased another exact value");
        Exact odd(thirds);odd.reset();positiveZero(odd);odd.add(1.L);odd.addQuotientProduct({-2.L},3.L);
        check(odd.compare(thirds)==0&&odd.rounded()==1.L/3.L,"reset lost the original odd-denominator reference");
        Exact span;span.addProduct({maximum,2.L});span.addQuotientProduct({q},3.L);
        span.reset();positiveZero(span);span.addProduct({1.5L,q});span.addProduct({.5L,q});
        check(span.rounded()==2*q,"large-span rational reset retained stale terms");
        Exact signedZero;signedZero.addProduct({-.5L,q});const auto negative=signedZero.rounded();
        check(negative==0&&std::signbit(negative),"reset fixture lacks a negative half-quantum RN zero");
        signedZero.reset();positiveZero(signedZero);signedZero.add(-q);
        check(signedZero.rounded()==-q,"reset positive zero could not accept a new negative quantum");
        // Compare the real fresh-constructor check with reset under the SAME
        // primed flags/errno. Capture before fresh's eventual destruction.
        {
            NativeEnvironment original;
            check(std::feclearexcept(FE_ALL_EXCEPT)==0&&std::feraiseexcept(FE_INEXACT)==0,
                "reset test cannot prime native exception flags");
            NativeEnvironment primed;
            check(primed.restore(),"reset test cannot restore primed constructor state");errno=173;
            Exact fresh;(void)fresh;NativeEnvironment constructed;
            check(primed.restore(),"reset test cannot restore primed reset state");errno=173;
            reused.reset();NativeEnvironment reset;
            check(constructed.same(reset),"reset differs from fresh construction in primed FP/errno state");
        }
        // Restore controls AND errno before inspecting the retained exact
        // value. FTZ and DAZ follow actual original-constructor admission.
        // The original environment probe does not promise raw-MXCSR rejection.
        const auto rejectedEnvironment=[&](const char* mode,const auto& change,bool mustReject){
            NativeEnvironment original;Exact value(thirds);const Exact prior(value);
            change();NativeEnvironment unsupported;bool freshRejected=false,resetRejected=false;std::string freshReason,resetReason;
            check(unsupported.restore(),"reset test cannot restore unsupported constructor controls/errno");
            try{Exact fresh;(void)fresh;}catch(const std::invalid_argument& e){freshRejected=true;freshReason=e.what();}
            NativeEnvironment afterFresh;
            check(unsupported.restore(),"reset test cannot restore unsupported reset controls/errno");
            try{value.reset();}catch(const std::invalid_argument& e){resetRejected=true;resetReason=e.what();}
            NativeEnvironment afterReset;
            check(original.restore(),"reset test failed to restore native controls/errno");
            check(afterFresh.same(afterReset),"reset differs from constructor rejection FP/CSR/errno state");
            check(freshRejected==resetRejected&&freshReason==resetReason,
                "reset changed the fresh-constructor environment rejection");
            std::cout<<"FINITE_CONTRACTION_ENV_DIAGNOSTIC,mode,"<<mode
                <<",fresh_rejected,"<<freshRejected<<",reset_rejected,"<<resetRejected
                <<",FP_errno_parity,"<<afterFresh.same(afterReset)
#if defined(_M_X64) || defined(__SSE2__)
                <<",unsupported_CSR,"<<unsupported.csr<<",fresh_CSR,"<<afterFresh.csr<<",reset_CSR,"<<afterReset.csr
#endif
                <<",fresh_errno,"<<afterFresh.error<<",reset_errno,"<<afterReset.error<<std::endl;
            check(original.restore(),"reset diagnostic failed to restore native controls/errno");
            if(mustReject)check(freshRejected,"required real constructor environment negative was not exercised");
            if(resetRejected)check(value.compare(prior)==0&&value.rounded()==1.L/3.L,
                "rejected reset mutated the pre-existing exact rational");
            else{const auto z=value.rounded();check(z==0&&!std::signbit(z),"accepted constructor-equivalent reset is not positive zero");}
            ++environmentPairs;return freshRejected;
        };
        rejectedEnvironment("FE_UPWARD",[&]{check(std::fesetround(FE_UPWARD)==0&&std::fegetround()==FE_UPWARD,
            "reset test did not enter a real directed rounding mode");},true);
#if defined(_MSC_VER) && defined(_M_X64)
        ftzRejected=rejectedEnvironment("FTZ_ONLY",[&]{_mm_setcsr((_mm_getcsr()&~0x0040u)|0x8000u);
            check((_mm_getcsr()&0x8040u)==0x8000u,"reset FTZ-only mode witness is missing");},false);++ftzPairs;
        dazRejected=rejectedEnvironment("DAZ_ONLY",[&]{_mm_setcsr((_mm_getcsr()&~0x8000u)|0x0040u);
            check((_mm_getcsr()&0x8040u)==0x0040u,"reset DAZ-only mode witness is missing");},false);++dazPairs;
#endif
        check(values==5&&environmentPairs==1+ftzPairs+dazPairs,"reset value/environment coverage incomplete");
        std::cout<<"FINITE_CONTRACTION_RESET_PASS,value_cases,"<<values<<",copy_independent,1,primed_FP_errno_pair,1"
            <<",environment_pairs,"<<environmentPairs<<",FTZ_only_pairs,"<<ftzPairs<<",DAZ_only_pairs,"<<dazPairs
            <<",FTZ_original_constructor_rejects,"<<ftzRejected<<",DAZ_original_constructor_rejects,"<<dazRejected<<",no_allocation_or_performance_claim\n";
    }
    // Range detection leaves an ordinary multiplication's native result and
    // sign intact; ordinary consumer statements retain their original order.
    volatile long double left=.75L,right=.625L;
    const auto ordinary=left*right;
    {LAMP::Hydrodynamics::Detail::HydroFiniteRangeScope range(true);const auto checked=left*right;
        LAMP::Hydrodynamics::Detail::hydroFiniteProductRange(left,right,checked);
        check(checked==ordinary&&std::signbit(checked)==std::signbit(ordinary),"ordinary product changed native bits");}
    bool marked=false;try{LAMP::Hydrodynamics::Detail::HydroFiniteRangeScope range(true);LAMP::Hydrodynamics::Detail::hydroFiniteProductRange(.5L,q,.5L*q);}
    catch(const LAMP::Hydrodynamics::Detail::HydroFiniteRange&){marked=true;}
    check(marked&&!LAMP::Hydrodynamics::Detail::hydroFiniteRangeEnabled(),"typed range marker/scope restoration failed");
    finiteRationalOracleOnly();
    std::cout<<"FINITE_CONTRACTION_RESULT,PASS,"<<std::numeric_limits<long double>::digits<<'\n';
}
// D12 selector reuses the original guards() numerical-range fixture only.
// It adds no scenario to the default matrix and changes no production guard.
void finiteRangeOnly()
{
    auto mesh=makeMesh(4);Case s;s.pureFirstJet=true;s.planar=true;s.planeOffset=2;const auto geometry=HydroMomentumGeometry::build(*mesh,metrics(*mesh,false));
    const auto f=frames(*mesh,geometry,publishAt(*mesh,s,{}),publishAt(*mesh,s,{1,0,.04,Position::Midpoint}),publishAt(*mesh,s,{1,0,.04,Position::End}),1000);
    const auto field=[](Vector2D p){return Vector2D{.3+.1*p.x,.2+.4*p.y};};const auto u=values(*mesh,geometry.layout(),field);
    const auto m=material(*mesh,f,u);const auto op=HydroKepConvection::build(*mesh,m);const auto boundary=bindBoundary(*mesh,m,traces(m,field));
    // Original positive range path and independent physical-power oracle.
    const auto sf=frames(*mesh,geometry,publishAt(*mesh,s,{}),publishAt(*mesh,s,{1,0,.04,Position::Midpoint}),publishAt(*mesh,s,{1,0,.04,Position::End}),1e-290);
    const auto sm=material(*mesh,sf,u);const auto sop=HydroKepConvection::build(*mesh,sm);const auto hugeField=[](Vector2D){return Vector2D{1e160,1e160};};
    std::cout<<"KEP_RANGE_BEGIN,positive_large_velocity_small_material"<<std::endl;
    const auto range=sop.apply(*mesh,values(*mesh,sm.layout(),hugeField),bindBoundary(*mesh,sm,traces(sm,hugeField)),sm.sources(),sm.timeLayer());
    std::array<long double,2> rangePower{},rangeGross{};
    for(std::size_t g=0;g<sm.faceRates().size();++g){const auto& segment=sm.geometry().segments()[g];if(inside(segment))continue;
        const auto flux=sm.faceRates()[g]*(segment.negativeVelocity>=0?1:-1);const auto work=(flux*1e160L)*.5e160L;const int a=axis(sm.layout(),owner(segment));rangePower[a]+=work;rangeGross[a]+=std::abs(work);}
    for(int a=0;a<2;++a)error(range.report().components[a].boundaryKineticPower,rangePower[a],rangeGross[a],1e-12,"finite boundary power with large velocity/small flux");
    std::cout<<"KEP_RANGE_POSITIVE,velocity=1e160,material=1e-290,PASS"<<std::endl;
    // Original negative range declaration, with the actual unchanged failure
    // reason printed so product/impulse/caller can be identified directly.
    if(std::numeric_limits<long double>::digits==std::numeric_limits<double>::digits){
        bool rejectedRange=false;
        std::cout<<"KEP_RANGE_BEGIN,old_minimum_denorm_velocity_negative"<<std::endl;
        try{(void)op.apply(*mesh,HydroNativeValues(m.layout(),std::vector<double>(m.layout().velocityCount(),std::numeric_limits<double>::denorm_min()),u.pressure()),boundary,m.sources(),m.timeLayer());}
        catch(const std::exception& error){const std::string reason=error.what();
            check(reason.find("unrepresentable")!=std::string::npos,"old minimum-denorm range failed for a different reason");
            std::cout<<"KEP_RANGE_NEGATIVE_REASON,"<<reason<<std::endl;rejectedRange=true;}
        check(rejectedRange,"old minimum-denorm range declaration was not rejected");
    }else std::cout<<"KEP_RANGE_NEGATIVE,not_applicable_to_wider_long_double"<<std::endl;
    std::cout<<"KEP_RANGE_ONLY,original_positive_and_negative_declarations,PASS"<<std::endl;
}
void guards()
{
    auto mesh=makeMesh(4);Case s;s.pureFirstJet=true;s.planar=true;s.planeOffset=2;const auto geometry=HydroMomentumGeometry::build(*mesh,metrics(*mesh,false));
    const auto f=frames(*mesh,geometry,publishAt(*mesh,s,{}),publishAt(*mesh,s,{1,0,.04,Position::Midpoint}),publishAt(*mesh,s,{1,0,.04,Position::End}),1000);
    const auto field=[](Vector2D p){return Vector2D{.3+.1*p.x,.2+.4*p.y};};const auto u=values(*mesh,geometry.layout(),field);const auto m=material(*mesh,f,u);const auto op=HydroKepConvection::build(*mesh,m);const auto t=traces(m,field);const auto boundary=bindBoundary(*mesh,m,t);const auto saved=m.faceRates();
    reject([&]{bindBoundary(*mesh,m,std::vector<double>(1));},"boundary trace dimensions");auto bad=t;bad[0]=std::numeric_limits<double>::infinity();reject([&]{bindBoundary(*mesh,m,bad);},"nonfinite boundary trace");
    for(std::size_t g=0;g<t.size();++g)if(inside(m.geometry().segments()[g])){bad=t;bad[g]=1;reject([&]{bindBoundary(*mesh,m,bad);},"nonzero unused internal");break;}
    const auto& affineActual=m.material().midpointRate().spatialAdvector();
    for(int i=0;i<m.layout().velocityCount();++i){const auto mean=nativeMeanReference(*mesh,affineActual,i);const auto p=affineActual.nativeVelocity(i);
        for(int c=0;c<2;++c)error(mean.mean[c],c==0?p.x:p.y,mean.absoluteMean[c]+std::abs(c==0?p.x:p.y),1e-12,"independent affine mean equals point");
        ++nativeMeanP1AffineChecks;}
    const auto other=material(*mesh,f,u);reject([&]{op.apply(*mesh,u,bindBoundary(*mesh,other,t),m.sources(),m.timeLayer());},"different prepared material owner");
    for(double strength:{-1.,1.1,std::numeric_limits<double>::quiet_NaN()})reject([&]{HydroKepConvection::buildControlled(*mesh,m,{strength,.2,.5});},"dissipation strength");
    for(const auto cfg:{HydroConvectionDissipationConfig{1,-.1,.5},HydroConvectionDissipationConfig{1,.5,.5},HydroConvectionDissipationConfig{1,.2,std::numeric_limits<double>::infinity()}})
        reject([&]{HydroKepConvection::buildControlled(*mesh,m,cfg);},"dissipation thresholds");
    const auto disabled=HydroKepConvection::buildControlled(*mesh,m,{0,.2,.5}).apply(*mesh,u,boundary,m.sources(),m.timeLayer());
    const auto originalResult=op.apply(*mesh,u,boundary,m.sources(),m.timeLayer());
    check(disabled.momentumLoads()==originalResult.momentumLoads() && disabled.momentumRates()==originalResult.momentumRates() && disabled.dissipationSamples().empty(),"zero strength changed center path");
    const auto stabilized=HydroKepConvection::buildControlled(*mesh,m,{});
    auto foreignSource=m.sources();foreignSource.advector={};
    reject([&]{stabilized.apply(*mesh,u,boundary,foreignSource,m.timeLayer());},"different actual advector identity");
    auto foreignClock=m.timeLayer();++foreignClock.step;
    reject([&]{stabilized.apply(*mesh,u,boundary,m.sources(),foreignClock);},"different full time layer");
    reject([&]{stabilized.apply(*mesh,u,bindBoundary(*mesh,other,t),m.sources(),m.timeLayer());},"different prepared material owner");
    reject([&]{stabilized.apply(*mesh,values(*mesh,metrics(*mesh,false).layout(),field),boundary,m.sources(),m.timeLayer());},"different layout identity");
    const auto zeroMaterial=material(*mesh,f,values(*mesh,m.layout(),[](Vector2D){return Vector2D{};}));
    for(auto F:zeroMaterial.faceRates())check(F==0,"zero source has nonzero material F");
    controlled(*mesh,HydroKepConvection::build(*mesh,zeroMaterial));
    // Finite physical power despite u*u being outside the machine range.
    // This is a numerical-range contract, not a physical density benchmark.
    const auto sf=frames(*mesh,geometry,publishAt(*mesh,s,{}),publishAt(*mesh,s,{1,0,.04,Position::Midpoint}),publishAt(*mesh,s,{1,0,.04,Position::End}),1e-290);
    const auto sm=material(*mesh,sf,u);const auto sop=HydroKepConvection::build(*mesh,sm);const auto hugeField=[](Vector2D){return Vector2D{1e160,1e160};};
    const auto range=sop.apply(*mesh,values(*mesh,sm.layout(),hugeField),bindBoundary(*mesh,sm,traces(sm,hugeField)),sm.sources(),sm.timeLayer());
    std::array<long double,2> rangePower{},rangeGross{};
    for(std::size_t g=0;g<sm.faceRates().size();++g){const auto& segment=sm.geometry().segments()[g];if(inside(segment))continue;
        const auto flux=sm.faceRates()[g]*(segment.negativeVelocity>=0?1:-1);const auto work=(flux*1e160L)*.5e160L;const int a=axis(sm.layout(),owner(segment));rangePower[a]+=work;rangeGross[a]+=std::abs(work);}
    for(int a=0;a<2;++a)error(range.report().components[a].boundaryKineticPower,rangePower[a],rangeGross[a],1e-12,"finite boundary power with large velocity/small flux");
    std::cout<<"RANGE,1e160,1e-290\n";
    for(int k=0;k<3;++k){auto source=m.sources();if(k==0)source.begin={};if(k==1)source.midpoint={};if(k==2)source.end={};reject([&]{op.apply(*mesh,u,boundary,source,m.timeLayer());},"different geometry identity");}
    auto wrongAdvector=m.sources();wrongAdvector.advector={};reject([&]{op.apply(*mesh,u,boundary,wrongAdvector,m.timeLayer());},"different actual advector identity");
    auto clock=m.timeLayer();clock.begin=-.01;reject([&]{op.apply(*mesh,u,boundary,m.sources(),clock);},"different full time layer");clock=m.timeLayer();clock.step+=1;reject([&]{HydroConvectionBoundary::bind(*mesh,m,t,m.sources(),clock);},"different full time layer");
    reject([&]{op.apply(*mesh,values(*mesh,metrics(*mesh,false).layout(),field),boundary,m.sources(),m.timeLayer());},"different layout identity");
    reject([&]{op.apply(*mesh,HydroNativeValues(m.layout(),std::vector<double>(m.layout().velocityCount(),1e308),u.pressure()),boundary,m.sources(),m.timeLayer());},"nonfinite");
    if(std::numeric_limits<long double>::digits==std::numeric_limits<double>::digits){reject([&]{op.apply(*mesh,HydroNativeValues(m.layout(),std::vector<double>(m.layout().velocityCount(),std::numeric_limits<double>::denorm_min()),u.pressure()),boundary,m.sources(),m.timeLayer());},"unrepresentable");}
    auto radial=makeMesh(4);Case rs=s;rs.radial=true;const auto rg=HydroMomentumGeometry::build(*radial,metrics(*radial,true));
    const auto rf=frames(*radial,rg,publishAt(*radial,rs,{}),publishAt(*radial,rs,{1,0,.04,Position::Midpoint}),publishAt(*radial,rs,{1,0,.04,Position::End}),1);
    const auto rfield=[](Vector2D p){return Vector2D{.1*p.x,.2};};const auto ru=values(*radial,rg.layout(),rfield);const auto rm=material(*radial,rf,ru);const auto rop=HydroKepConvection::build(*radial,rm);const auto rt=traces(rm,rfield);const auto rb=bindBoundary(*radial,rm,rt);
    auto wrong=ru.normalVelocity();wrong[rg.metrics().axisVelocityDofs().front()]=1;reject([&]{rop.apply(*radial,HydroNativeValues(rg.layout(),wrong,ru.pressure()),rb,rm.sources(),rm.timeLayer());},"nonzero radial axis velocity");
    for(std::size_t g=0;g<rt.size();++g){const auto& segment=rm.geometry().segments()[g];if(!inside(segment) && segment.normalAxis==0 && segment.begin.x==0 && axis(rm.layout(),owner(segment))==0){bad=rt;bad[g]=1;reject([&]{bindBoundary(*radial,rm,bad);},"nonzero radial axis boundary trace");break;}}
    const auto narrow=HydroMomentumGeometry::build(*radial,metrics(*radial,true,3));
    const auto nf=frames(*radial,narrow,publishAt(*radial,rs,{}),publishAt(*radial,rs,{1,0,.04,Position::Midpoint}),publishAt(*radial,rs,{1,0,.04,Position::End}),1);
    const auto nm=material(*radial,nf,values(*radial,narrow.layout(),rfield));reject([&]{HydroKepConvection::build(*radial,nm);},"unidentifiable affine common trace");
    auto equivalent=makeMesh(4);mesh.reset();op.apply(*equivalent,u,boundary,m.sources(),m.timeLayer());auto& vertex=equivalent->verts(equivalent->vertIndexToLabel[0]);const auto original=vertex.coord();vertex.setCoord(original+Vector2D{.01,0});
    reject([&]{op.apply(*equivalent,u,boundary,m.sources(),m.timeLayer());},"NativeMeshIdentity");
    reject([&]{stabilized.apply(*equivalent,u,boundary,m.sources(),m.timeLayer());},"NativeMeshIdentity");vertex.setCoord(original);DictionaryMesh::RefineBox box;box.level=1;box.min={.25,.25};box.max={.75,.75};equivalent->applyRefineBoxes({box});reject([&]{op.apply(*equivalent,u,boundary,m.sources(),m.timeLayer());},"NativeMeshIdentity");
    reject([&]{stabilized.apply(*equivalent,u,boundary,m.sources(),m.timeLayer());},"NativeMeshIdentity");check(m.faceRates()==saved,"rejections changed authoritative flux");
}
void regrid(bool rz)
{
    auto coarse=makeMesh(4),fine=makeMesh(4,1);Case s;s.pureFirstJet=true;s.radial=rz;s.planar=true;s.planeOffset=.49;s.motion=1;s.speed=.4;
    const auto start=initial(*coarse,s);HermiteGalsPublicationConfig config;config.volume.symmetryPlanes=s.symmetry(*fine);
    const auto transferred=transferHermiteGalsFirstJetMesh(*coarse,start,*fine,HermiteGalsMeshTransferKind::ExactRefinement,config);const auto a=transferred.publication();std::string stage;
    check(transferred.sourceIdentity()==start.state().identity() && a.state().timeLayer().time()==0,"refinement lost source/time");
    const auto mid=advanceFirstJet(*fine,a,s,{1,0,.1,Position::Midpoint},stage),end=advanceFirstJet(*fine,a,s,{1,0,.1,Position::End},stage);const Vector2D v=rz?Vector2D{0,.4}:Vector2D{.4,0};const auto field=[v](Vector2D){return v;};
    const auto fg=HydroMomentumGeometry::build(*fine,metrics(*fine,rz));const auto fm=material(*fine,frames(*fine,fg,a,mid,end,1e6),values(*fine,fg.layout(),field));const auto op=HydroKepConvection::build(*fine,fm);constant(*fine,op,v,true);++regrids;
    controlled(*fine,op);
    config.volume.symmetryPlanes=s.symmetry(*coarse);const auto projected=transferHermiteGalsFirstJetMesh(*fine,end,*coarse,HermiteGalsMeshTransferKind::FirstJetProjection,config);const auto b=projected.publication();check(b.state().timeLayer().time()==.1,"coarsening changed physical time");
    const auto fineU=values(*fine,fg.layout(),field);const auto fineBoundary=bindBoundary(*fine,fm,traces(fm,field));
    reject([&]{op.apply(*coarse,fineU,fineBoundary,fm.sources(),fm.timeLayer());},"NativeMeshIdentity");const auto cm=advanceFirstJet(*coarse,b,s,{2,.1,.2,Position::Midpoint},stage),ce=advanceFirstJet(*coarse,b,s,{2,.1,.2,Position::End},stage);
    const auto cg=HydroMomentumGeometry::build(*coarse,metrics(*coarse,rz));const auto cmat=material(*coarse,frames(*coarse,cg,b,cm,ce,1e6),values(*coarse,cg.layout(),field));
    check(cmat.sources().begin==b.state().identity() && cmat.timeLayer().step==2 && !cmat.layout().sameIdentity(fm.layout()),"coarsening retagged old material interval/layout");
    constant(*coarse,HydroKepConvection::build(*coarse,cmat),v,true);++regrids;std::cout<<"REGRID,"<<rz<<','<<cmat.timeLayer().step<<','<<cmat.timeLayer().begin<<','<<cmat.timeLayer().end<<'\n';
    controlled(*coarse,HydroKepConvection::build(*coarse,cmat));
}
void boundedPorts(bool rz,bool amr,double rho)
{
    auto mesh=makeMesh(8,amr?1:0);Case s;s.pureFirstJet=true;s.radial=rz;s.planar=true;s.planeOffset=2;
    const auto g=HydroMomentumGeometry::build(*mesh,metrics(*mesh,rz));const auto& l=g.layout();
    const auto f=frames(*mesh,g,publishAt(*mesh,s,{}),publishAt(*mesh,s,{1,0,.04,Position::Midpoint}),publishAt(*mesh,s,{1,0,.04,Position::End}),rho);
    const auto m=material(*mesh,f,values(*mesh,l,[rz](Vector2D){return rz?Vector2D{0,1}:Vector2D{1,0};}));const auto saved=m.faceRates();
    std::map<std::tuple<int,int,double,double>,std::vector<int>> groups;
    for(std::size_t k=0;k<g.segments().size();++k){const auto& p=g.segments()[k];if(inside(p) || m.faceRates()[k]==0)continue;
        groups[{axis(l,owner(p)),p.normalAxis,p.normalAxis==0?p.begin.y:p.begin.x,p.normalAxis==0?p.end.y:p.end.x}].push_back(static_cast<int>(k));}
    std::vector<HydroConvectionPeriodicPair> ports;for(const auto& group:groups){check(group.second.size()==2,"test unmatched port");const int a=group.second[0],b=group.second[1];ports.push_back({a,b,g.segments()[a].begin-g.segments()[b].begin});}
    const auto op=HydroKepConvection::build(*mesh,m,ports);const auto& report=op.momentReport();
    check(report.boundedRows>0 && report.coefficientRatio<=1+1e-12 && report.quadraticBudgetRatio<=1+1e-10,"bounded native certificate");
    auto field=values(*mesh,l,[rz](Vector2D p){const double t=static_cast<double>(2*pi)*(rz?p.y:p.x);return rz?Vector2D{p.x*std::sin(t),.4+std::cos(t)}:Vector2D{std::cos(t),std::sin(t)};});
    auto u=field.normalVelocity();for(int i=0;i<l.velocityCount();++i)u[i]=u[op.representatives()[i]];field=HydroNativeValues(l,u,field.pressure());
    const auto result=op.applyPaired(*mesh,field,m.sources(),m.timeLayer());const int n=l.velocityCount();
    std::vector<std::map<int,long double>> matrix(n);std::vector<int> partners(g.segments().size(),-1);for(const auto& p:ports){partners[p.firstSegment]=p.secondSegment;partners[p.secondSegment]=p.firstSegment;}
    for(std::size_t k=0;k<g.segments().size();++k){const auto& p=g.segments()[k];const auto F=m.faceRates()[k];
        if(inside(p)){const int i=p.negativeVelocity,j=p.positiveVelocity;matrix[i][i]+=.5L*F;matrix[j][j]-=.5L*F;
            for(const auto& w:op.weights()[k]){matrix[i][w.velocity]+=.5L*F*w.weight;matrix[j][w.velocity]-=.5L*F*w.weight;matrix[w.velocity][i]-=.5L*F*w.weight;matrix[w.velocity][j]+=.5L*F*w.weight;}}
        else if(partners[k]>=0){const int i=owner(p),j=owner(g.segments()[partners[k]]);const auto v=.5L*(p.negativeVelocity>=0?1:-1)*F;matrix[i][i]+=v;matrix[i][j]+=v;}}
    for(const auto& p:op.couplings()){matrix[p.first][p.second]+=p.coefficient;matrix[p.second][p.first]-=p.coefficient;}
    std::vector<std::vector<int>> classes(n);for(int i=0;i<n;++i)classes[op.representatives()[i]].push_back(i);
    long double maxError=0;for(const auto& members:classes)if(!members.empty()){
        long double classMass=0,classRate=0;for(int i:members){classMass+=m.midpointMass().masses()[i];for(const auto& e:matrix[i])classRate+=e.second*u[e.first];classRate-=m.cells()[i].outwardRate*u[i];}
        const int root=members.front();std::map<int,long double> skew;long double gross=0,q=0;double h=0;
        for(int i:members){q+=m.cells()[i].outwardRate;const auto& box=g.volumes()[i];h=std::max({h,box.upper.x-box.lower.x,box.upper.y-box.lower.y});
            for(const auto& e:matrix[i])skew[op.representatives()[e.first]]+=e.second;}
        skew[root]-=.5L*q;for(std::size_t k=0;k<g.segments().size();++k)for(int i:{g.segments()[k].negativeVelocity,g.segments()[k].positiveVelocity})if(i>=0 && op.representatives()[i]==root)gross+=std::abs(m.faceRates()[k]);
        const auto p=l.domain().face(*mesh,l.nativeFace(root)).center;std::array<long double,6> moments{};long double norm=0;
        for(const auto& e:skew){auto d=l.domain().face(*mesh,l.nativeFace(e.first)).center-p;if(rz)d.y-=std::round(d.y);else d.x-=std::round(d.x);const double x=d.x/h,y=d.y/h;const std::array<double,6> polynomial{1,x,y,x*x,x*y,y*y};
            norm+=std::abs(e.second);for(int c=0;c<6;++c)moments[c]+=e.second*polynomial[c];}
        error(moments[0],.5L*q,gross,1e-12,"independent constant moment");
        error(moments[1],rz?0:classMass/h,gross,1e-12,"independent point PDE x moment");error(moments[2],rz?classMass/h:0,gross,1e-12,"independent point PDE y moment");
        check(norm<=2*gross*(1+1e-12),"independent final l1 bound");for(int c=3;c<6;++c)check(std::abs(moments[c])<=gross*h+1e-12*gross,"independent final quadratic bound");
        for(int i:members){const auto expected=-classRate*m.midpointMass().masses()[i]/classMass+u[i]*(m.cells()[i].defectRate-m.cells()[i].outwardRate);
            maxError=std::max(maxError,error(result.momentumLoads()[i],expected,rho,1e-12,"independent periodic native load"));}}
    for(const auto& r:result.report().components){check(r.momentumIdentity<=1e-12 && r.kineticIdentity<=1e-12,"periodic material KE");}
    auto constantField=values(*mesh,l,[rz](Vector2D){return Vector2D{rz?0:.7,-.3};});const auto constantResult=op.applyPaired(*mesh,constantField,m.sources(),m.timeLayer());
    for(int i=0;i<n;++i)error(constantResult.momentumLoads()[i],0,m.midpointMass().masses()[i],1e-10,"periodic constant acceleration");
    const auto transverse=values(*mesh,l,[rz](Vector2D p){const double z=rz?p.x:p.y;return Vector2D{z,z};});
    const auto affineResult=op.applyPaired(*mesh,transverse,m.sources(),m.timeLayer());
    for(int i=0;i<n;++i)error(affineResult.momentumLoads()[i],0,m.midpointMass().masses()[i],1e-10,"periodic transverse affine acceleration");
    auto missing=ports;missing.pop_back();reject([&]{HydroKepConvection::build(*mesh,m,missing);},"unpaired nonzero material boundary rate");
    auto wrong=ports;wrong.front().translation=wrong.front().translation*.5;reject([&]{HydroKepConvection::build(*mesh,m,wrong);},"periodic segment begin mismatch");
    check(&op.material().faceRates()==&m.faceRates() && saved==m.faceRates(),"bounded solve changed shared material flux");
    std::cout<<"BOUNDED_PORTS,"<<rz<<','<<amr<<','<<rho<<','<<report.boundedRows<<','<<report.iterations<<','<<report.coefficientRatio<<','<<maxError<<std::endl;
    controlled(*mesh,op,true,true);
}
template<class F>void migrationCheck(const std::string& name,F f)
{try{f();}catch(const std::exception& e){++migrationFailures;std::cout<<"MIGRATION_FAILURE,"<<name<<','<<e.what()<<std::endl;}}
}
int main(int argc,char** argv)
{
    try{
        std::cout<<std::setprecision(17);
        if(argc==2&&std::string(argv[1])=="--finite-contraction-only"){finiteContractionOnly();return 0;}
        if(argc==2&&std::string(argv[1])=="--range-only"){finiteRangeOnly();return 0;}
        if(argc==2&&std::string(argv[1])=="--impermeable-only"){closed(false,0);physicalPairSpatial();check(impermeableCases==1&&frozenOperatorCases==1,"missing actual impermeable/frozen operator contract");check(nativeMeanP1AnalyticCases==1&&nativeMeanP1OldPointNegative==1,"missing independent source mean group");return 0;}
        if(argc==2&&std::string(argv[1])=="--physical-pair-spatial"){physicalPairSpatial();return 0;}
        check(argc==1,"unexpected convection contract argument");
        for(bool rz:{false,true})for(int amr:{0,1}){const auto name=std::to_string(rz)+","+std::to_string(amr);for(int selection=0;selection<4;++selection) migrationCheck("algebra,"+name+","+std::to_string(selection),[&]{algebra(rz,amr,selection);});
            for(double rho:{1.,1000.,1e6})migrationCheck("bounded,"+name+","+std::to_string(rho),[&]{boundedPorts(rz,amr,rho);});
            for(bool planar:{false,true})migrationCheck("crossing,"+name+","+std::to_string(planar),[&]{crossing(rz,amr,planar);});migrationCheck("closed,"+name,[&]{closed(rz,amr);});}migrationCheck("physical_pair_spatial",physicalPairSpatial);migrationCheck("guards",[&]{guards();});migrationCheck("regrid,0",[&]{regrid(false);});migrationCheck("regrid,1",[&]{regrid(true);});
        check(amrMeanGap>1e-3,"AMR arithmetic mean location control absent");check(amrRawPowerGap>1e-5,"AMR raw reconstructed power control absent");check(actualDefect>1e-5 && actualTransfer>1e-5,"actual material defect/geometry transfer erased");
        std::cout<<"CONTROL,"<<amrMeanGap<<','<<amrRawPowerGap<<','<<actualDefect<<','<<actualTransfer<<'\n';
        std::cout<<"GEOMETRY,"<<coarseFineChecks<<','<<alignedChecks<<'\n';
        std::cout<<"NATIVE_MOMENTS,"<<momentChecks<<','<<pairChecks<<','<<peakMoment<<'\n';
        std::cout<<"SUMMARY,"<<scenarios<<','<<traceChecks<<','<<loadChecks<<','<<actualCrossings<<','<<rejected<<','<<regrids<<','<<peakTrace<<','<<peakLoad<<','<<peakMomentum<<','<<peakEnergy<<','<<peakConstant<<','<<peakLinearity<<'\n';
        check(materialEnergyCases==scenarios&&materialEnergyDofs==loadChecks,"independent material energy coverage");
        for(int k=0;k<3;++k)check(materialEnergyNegative[k]==1&&materialEnergyMutation[k]>1e-5,"ineffective material energy mutation");
        std::cout<<"MATERIAL_ENERGY_SUMMARY,"<<materialEnergyCases<<','<<materialEnergyDofs;for(auto e:materialEnergyError)std::cout<<','<<e;std::cout<<",PASS\n";
        std::cout<<"MATERIAL_ENERGY_NEGATIVE";for(auto e:materialEnergyMutation)std::cout<<','<<e;std::cout<<",PASS\n";
        std::cout<<"CONTROLLED_SUMMARY,"<<controlledCases<<','<<controlledConnections<<','<<controlledActive<<','<<controlledQuiet<<','<<controlledCoarseFine<<','<<controlledPeakDefect<<','<<controlledPeak<<'\n';
        check(impermeableCases==1&&frozenOperatorCases==1&&impermeableFreeRows>0&&impermeableWallTangentRows>0,"missing impermeable/frozen free-row contract coverage");
        check(nativeMeanP1AnalyticCases==1&&nativeMeanP1OldPointNegative==1&&nativeMeanP1AffineChecks>0,"missing independent source mean/affine/mutation coverage");
        check(physicalPairSpatialFamilies==1 && physicalPairSpatialRows==616 && physicalPairWallRows==100,"missing physical-pair P3 spatial family coverage");
        if(migrationFailures==0)check(controlledCases==57 && controlledActive>0 && controlledQuiet>0 && controlledCoarseFine>0 && controlledPeakDefect>1e-5,"incomplete controlled dissipation coverage");
        std::cout<<"MIGRATION_RESULT,"<<migrationFailures<<','<<(migrationFailures?"FAIL":"PASS")<<'\n';return migrationFailures?1:0;
    }catch(const std::exception& e){std::cerr<<"H10c FAIL: "<<e.what()<<'\n';return 1;}
}
