#include "Modules/Hydrodynamics/CchKepConvection.h" // Public header stands alone.
#include "interface_first_jet_fixture.h"
#include <sstream>
#include <set>
#include <map>
#include <tuple>

namespace
{
using namespace Quantum::Hydrodynamics;
int scenarios=0,traceChecks=0,coarseFineChecks=0,alignedChecks=0,loadChecks=0,actualCrossings=0,rejected=0,regrids=0;
long double peakTrace=0,peakLoad=0,peakMomentum=0,peakEnergy=0,peakConstant=0,peakLinearity=0;
long double amrMeanGap=0,amrRawPowerGap=0,actualDefect=0,actualTransfer=0;
int momentChecks=0,pairChecks=0;long double peakMoment=0;
int migrationFailures=0;
int controlledCases=0,controlledConnections=0,controlledActive=0,controlledQuiet=0,controlledCoarseFine=0;
long double controlledPeak=0,controlledPeakDefect=0;
int materialEnergyCases=0,materialEnergyDofs=0;
std::array<int,3> materialEnergyNegative{};
std::array<long double,5> materialEnergyError{};
std::array<long double,3> materialEnergyMutation{};
static_assert(!std::is_default_constructible_v<CchKepConvection>);
static_assert(!std::is_default_constructible_v<CchConvectionBoundary>);
static_assert(std::is_same_v<decltype(std::declval<CchKepConvectionResult&>().momentumRates()),const std::vector<long double>&>);
static_assert(std::is_same_v<decltype(std::declval<CchKepConvection&>().couplings()),const std::vector<CchConvectionCoupling>&>);
long double error(long double a,long double b,long double scale,long double limit,const std::string& why)
{
    const auto e=std::abs(a-b)/(scale>0?scale:1);std::ostringstream message;
    message<<std::setprecision(17)<<why<<" actual="<<a<<" expected="<<b<<" scale="<<scale<<" error="<<e;
    check(std::isfinite(a) && std::isfinite(b) && std::isfinite(e) && e<=limit,message.str());return e;
}
int axis(const CchDofLayout& l,int u){return l.normalAxis(u)==CchDofLayout::NormalAxis::First?0:1;}
bool inside(const CchMomentumSegment& g){return g.negativeVelocity>=0 && g.positiveVelocity>=0;}
int owner(const CchMomentumSegment& g){return g.negativeVelocity>=0?g.negativeVelocity:g.positiveVelocity;}
bool coarseFine(const QuadtreeMesh& mesh,const CchMomentumGeometry& geometry,const CchMomentumSegment& g)
{
    int lo=100,hi=-1;for(int u:{g.negativeVelocity,g.positiveVelocity})for(const auto& overlap:geometry.volumes()[u].overlaps){
        const int level=mesh.cells(mesh.cellIndexToLabel[geometry.layout().nativeCell(overlap.pressure)]).level();lo=std::min(lo,level);hi=std::max(hi,level);}
    return lo<hi;
}
CchPhysicalDivergence metrics(const QuadtreeMesh& mesh,bool rz,int selection=0)
{
    FaceSegmentGraph graph;graph.build(mesh,[selection](const QuadtreeMesh&,const Cell& c){const auto p=c.centerCoord();
        if(selection==1)return p.x<.5;if(selection==2)return !(p.x>.25 && p.x<.75 && p.y>.25 && p.y<.75);
        if(selection==3)return p.x<.25 || p.x>.75;return true;});
    return CchPhysicalDivergence::build(mesh,CchIntegralDivergence::build(mesh,CchDofLayout::build(mesh,NativeLeafDomain::build(mesh,graph,rz))));
}
CchNativeValues values(const QuadtreeMesh& mesh,const CchDofLayout& l,const std::function<Vector2D(Vector2D)>& f)
{
    std::vector<double> u;for(int i=0;i<l.velocityCount();++i){const auto face=l.domain().face(mesh,l.nativeFace(i));const auto v=f(face.center);u.push_back(face.horizontal?v.y:v.x);}
    return CchNativeValues(l,u,std::vector<double>(l.pressureCount()));
}
Vector2D centroid(const CchMomentumSegment& g,bool rz)
{
    auto p=g.begin+(g.end-g.begin)*.5;
    if(rz && g.normalAxis==1){const long double a=g.begin.x,b=g.end.x;p.x=static_cast<double>((2.L/3)*(b*b*b-a*a*a)/(b*b-a*a));}
    return p;
}
std::vector<double> traces(const CchDualMaterial& m,const std::function<Vector2D(Vector2D)>& field)
{
    std::vector<double> b(m.faceRates().size());for(std::size_t g=0;g<b.size();++g)if(!inside(m.geometry().segments()[g])){
        const auto& s=m.geometry().segments()[g];const auto p=centroid(s,m.layout().domain().axisymmetric());const auto v=field(p);b[g]=axis(m.layout(),owner(s))==0?v.x:v.y;
    }return b;
}
CchConvectionBoundary bindBoundary(const QuadtreeMesh& mesh,const CchDualMaterial& m,const std::vector<double>& t)
{return CchConvectionBoundary::bind(mesh,m,t,m.sources(),m.timeLayer());}
HermiteGalsPublication publishAt(const QuadtreeMesh& mesh,const Case& s,HermiteGalsTimeLayer time)
{
    const auto f=s.profile(time.time());const auto seed=stateFor(mesh,f,s.radial,time,[](int,Vector2D,Vector2D){return 0.;},s.symmetry(mesh),HermiteGalsTransport::FirstJet);
    HermiteGalsCandidate candidate(mesh,time,seed.primary());const auto mixed=closeHermiteGalsMixedAmr(mesh,candidate,s.symmetry(mesh));
    HermiteGalsPublicationConfig config;config.volume.symmetryPlanes=s.symmetry(mesh);const auto id=candidate.identity();
    return publishHermiteGalsFirstJetGeometry(mesh,std::move(candidate),mixed,{id,centers(mesh,f)},config);
}
struct Frames {CchPhaseGeometry a,m,b;CchMomentumMass ma,mm,mb;};
Frames frames(const QuadtreeMesh& mesh,const CchMomentumGeometry& g,HermiteGalsPublication a,HermiteGalsPublication m,HermiteGalsPublication b,double rho)
{
    const auto pa=CchPhaseGeometry::build(mesh,g.metrics(),a),pm=CchPhaseGeometry::build(mesh,g.metrics(),m),pb=CchPhaseGeometry::build(mesh,g.metrics(),b);
    return {pa,pm,pb,CchMomentumMass::build(mesh,g,pa,rho,1),CchMomentumMass::build(mesh,g,pm,rho,1),CchMomentumMass::build(mesh,g,pb,rho,1)};
}
CchDualMaterial material(const QuadtreeMesh& mesh,const Frames& f,const CchNativeValues& u,const CchAdvectorBoundary* boundary=nullptr)
{
    const auto reconstruction=CchAdvectorReconstruction::prepare(mesh,f.mm.geometry());
    const auto actual=boundary?reconstruction.evaluate(mesh,u,f.m.geometryIdentity(),f.m.timeLayer(),*boundary):
        reconstruction.evaluate(mesh,u,f.m.geometryIdentity(),f.m.timeLayer());
    const auto q=CchPhaseRate::prepareSpatial(mesh,f.m,actual);
    return CchDualMaterial::prepare(mesh,CchMaterialTransport::prepare(mesh,f.a,q,f.b,f.ma.liquidDensity(),1),f.ma,f.mm,f.mb);
}
std::vector<long double> outflow(const CchDualMaterial& m)
{
    std::vector<long double> q(m.layout().velocityCount());for(std::size_t g=0;g<m.faceRates().size();++g){const auto& s=m.geometry().segments()[g];
        if(s.negativeVelocity>=0)q[s.negativeVelocity]+=m.faceRates()[g];if(s.positiveVelocity>=0)q[s.positiveVelocity]-=m.faceRates()[g];}return q;
}
void nativeMoments(const QuadtreeMesh& mesh,const CchKepConvection& op)
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
void affineTraces(const QuadtreeMesh& mesh,const CchKepConvection& op,bool amr)
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
void materialEnergyReference(const CchDualMaterial& m,const CchNativeValues& u,
    const CchKepConvectionResult& result,const std::vector<double>& trace)
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
CchKepConvectionResult inspect(const QuadtreeMesh& mesh,const CchKepConvection& op,
    const CchNativeValues& u,const std::vector<double>& t,bool amr=false,bool ordinary=false)
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
    std::vector<long double> correction(raw.size());
    for(const auto& e:op.couplings()){
        check(e.first>=0 && e.first<e.second && e.second<n && axis(op.layout(),e.first)==axis(op.layout(),e.second),"invalid shared native pair");
        correction[static_cast<std::size_t>(e.first)*n+e.second]+=e.coefficient;
        correction[static_cast<std::size_t>(e.second)*n+e.first]-=e.coefficient;
    }
    const auto q=outflow(m);long double rawPower=0,desiredRawPower=0,grossRawPower=0;
    for(int i=0;i<n;++i){long double rawValue=0,rowScale=0;for(int j=0;j<n;++j){expected[i]+=.5L*(raw[static_cast<std::size_t>(i)*n+j]-raw[static_cast<std::size_t>(j)*n+i])*u.normalVelocity()[j];rawValue+=rawCf[static_cast<std::size_t>(i)*n+j]*u.normalVelocity()[j];rowScale+=.5L*(rawAbs[static_cast<std::size_t>(i)*n+j]+rawAbs[static_cast<std::size_t>(j)*n+i])*std::abs(u.normalVelocity()[j]);}
        for(int j=0;j<n;++j){const auto k=correction[static_cast<std::size_t>(i)*n+j];expected[i]+=k*u.normalVelocity()[j];rowScale+=std::abs(k*u.normalVelocity()[j]);}
        const int a=axis(op.layout(),i);const long double v=u.normalVelocity()[i];expected[i]+=.5L*qInternal[i]*v+boundary[i];
        const auto load=-expected[i]+v*m.cells()[i].integratedDefect/m.deltaTime();
        rowScale+=std::abs(.5L*qInternal[i]*v)+std::abs(boundary[i]);
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
void constant(const QuadtreeMesh& mesh,const CchKepConvection& op,Vector2D constantVelocity,bool actual)
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
void controlled(const QuadtreeMesh& mesh,const CchKepConvection& center,bool quadratic=true,bool paired=false)
{
    const auto& m=center.material();const auto& l=m.layout();const int n=l.velocityCount();const bool rz=l.domain().axisymmetric();
    const auto op=paired?CchKepConvection::buildControlled(mesh,m,center.periodicPairs(),{}):CchKepConvection::buildControlled(mesh,m,{});
    const auto saved=m.faceIntegrals();
    const auto apply=[&](const CchKepConvection& c,CchNativeValues u,const std::function<Vector2D(Vector2D)>& field){
        auto v=u.normalVelocity();for(int i=0;i<n;++i)v[i]=v[c.representatives()[i]];u=CchNativeValues(l,v,u.pressure());
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
    auto u=values(mesh,l,rough);auto velocity=u.normalVelocity();for(int i=0;i<n;++i)velocity[i]=velocity[op.representatives()[i]];u=CchNativeValues(l,velocity,u.pressure());
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
    const auto geometry=CchMomentumGeometry::build(*mesh,metrics(*mesh,rz,selection));const auto f=frames(*mesh,geometry,publishAt(*mesh,s,{}),publishAt(*mesh,s,{1,0,.04,Position::Midpoint}),publishAt(*mesh,s,{1,0,.04,Position::End}),1000);
    const auto advector=values(*mesh,geometry.layout(),[rz](Vector2D p){return Vector2D{(rz?p.x:1)*(.3+std::sin(3*p.y)),.1+std::cos(4*p.x)};});
    const auto op=CchKepConvection::build(*mesh,material(*mesh,f,advector));nativeMoments(*mesh,op);affineTraces(*mesh,op,amr!=0);constant(*mesh,op,{rz?0:.37,-.21},false);
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
    const auto geometry=CchMomentumGeometry::build(*mesh,metrics(*mesh,rz));const auto a=initial(*mesh,s);std::string stage;
    const auto mid=advanceFirstJet(*mesh,a,s,{1,0,.1,Position::Midpoint},stage),end=advanceFirstJet(*mesh,a,s,{1,0,.1,Position::End},stage);
    const Vector2D v=rz?Vector2D{0,.4}:Vector2D{.4,0};const auto advector=values(*mesh,geometry.layout(),[v](Vector2D){return v;});
    for(double rho:{1.,1000.,1e6}){const auto op=CchKepConvection::build(*mesh,material(*mesh,frames(*mesh,geometry,a,mid,end,rho),advector));constant(*mesh,op,v,true);
        const auto probe=[rz](Vector2D p){return Vector2D{(rz?p.x:1)*(.2+std::sin(p.y)),.1+std::cos(2*p.x)};};
        inspect(*mesh,op,values(*mesh,op.layout(),probe),traces(op.material(),probe));++actualCrossings;
        controlled(*mesh,op);
        std::cout<<"CROSSING,"<<rz<<','<<amr<<','<<planar<<','<<rho<<'\n';}
}
void closed(bool rz,int amr)
{
    auto mesh=makeMesh(4,amr);Case s;s.pureFirstJet=true;s.radial=rz;s.planar=true;s.planeOffset=2;
    const auto geometry=CchMomentumGeometry::build(*mesh,metrics(*mesh,rz));const auto& l=geometry.layout();
    const auto stream=[rz](Vector2D p){const double radial=rz?p.x*p.x*(1-p.x)*(1-p.x):p.x*(1-p.x);return radial*p.y*(1-p.y);};
    std::vector<double> transported;
    for(int i=0;i<l.velocityCount();++i){const auto face=l.domain().face(*mesh,l.nativeFace(i));const auto& edge=face.horizontal?mesh->hEdges(face.edgeLabel):mesh->vEdges(face.edgeLabel);
        auto a=mesh->verts(edge.edgeToVert[0]).coord(),b=mesh->verts(edge.edgeToVert[1]).coord();if(face.horizontal?a.x>b.x:a.y>b.y)std::swap(a,b);
        const auto area=geometry.metrics().faceMeasures()[i];const double volumeFlux=(face.horizontal?-1:1)*(rz?static_cast<double>(2*pi):1)*(stream(b)-stream(a));
        transported.push_back(area>0?volumeFlux/area:0);}
    const auto f=frames(*mesh,geometry,publishAt(*mesh,s,{}),publishAt(*mesh,s,{1,0,.04,Position::Midpoint}),publishAt(*mesh,s,{1,0,.04,Position::End}),1);
    std::vector<int> walls;for(int u=0;u<l.velocityCount();++u)if(l.domain().face(*mesh,l.nativeFace(u)).outwardSign!=0)walls.push_back(u);
    const auto boundary=CchAdvectorBoundary::impermeable(*mesh,geometry,walls);
    const auto m=material(*mesh,f,CchNativeValues(l,transported,std::vector<double>(l.pressureCount())),&boundary);const auto op=CchKepConvection::build(*mesh,m);
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
    // Preserve the old zero-power number as a diagnostic, not an extra R=0 gate.
    std::cout<<"CLOSED,"<<rz<<','<<amr<<','<<std::abs(power)/std::max(1.L,gross)<<'\n';
    std::cout<<"CLOSED_MATERIAL_POWER,"<<rz<<','<<amr<<','<<power<<','<<reference<<','<<gross<<','<<referenceGross<<','<<qMax<<','<<e<<'\n';
}
template<class F>void reject(F f,const std::string& why){rejects(f,why);++rejected;}
void guards()
{
    auto mesh=makeMesh(4);Case s;s.pureFirstJet=true;s.planar=true;s.planeOffset=2;const auto geometry=CchMomentumGeometry::build(*mesh,metrics(*mesh,false));
    const auto f=frames(*mesh,geometry,publishAt(*mesh,s,{}),publishAt(*mesh,s,{1,0,.04,Position::Midpoint}),publishAt(*mesh,s,{1,0,.04,Position::End}),1000);
    const auto field=[](Vector2D p){return Vector2D{.3+.1*p.x,.2+.4*p.y};};const auto u=values(*mesh,geometry.layout(),field);const auto m=material(*mesh,f,u);const auto op=CchKepConvection::build(*mesh,m);const auto t=traces(m,field);const auto boundary=bindBoundary(*mesh,m,t);const auto saved=m.faceRates();
    reject([&]{bindBoundary(*mesh,m,std::vector<double>(1));},"boundary trace dimensions");auto bad=t;bad[0]=std::numeric_limits<double>::infinity();reject([&]{bindBoundary(*mesh,m,bad);},"nonfinite boundary trace");
    for(std::size_t g=0;g<t.size();++g)if(inside(m.geometry().segments()[g])){bad=t;bad[g]=1;reject([&]{bindBoundary(*mesh,m,bad);},"nonzero unused internal");break;}
    const auto other=material(*mesh,f,u);reject([&]{op.apply(*mesh,u,bindBoundary(*mesh,other,t),m.sources(),m.timeLayer());},"different prepared material owner");
    for(double strength:{-1.,1.1,std::numeric_limits<double>::quiet_NaN()})reject([&]{CchKepConvection::buildControlled(*mesh,m,{strength,.2,.5});},"dissipation strength");
    for(const auto cfg:{CchConvectionDissipationConfig{1,-.1,.5},CchConvectionDissipationConfig{1,.5,.5},CchConvectionDissipationConfig{1,.2,std::numeric_limits<double>::infinity()}})
        reject([&]{CchKepConvection::buildControlled(*mesh,m,cfg);},"dissipation thresholds");
    const auto disabled=CchKepConvection::buildControlled(*mesh,m,{0,.2,.5}).apply(*mesh,u,boundary,m.sources(),m.timeLayer());
    const auto originalResult=op.apply(*mesh,u,boundary,m.sources(),m.timeLayer());
    check(disabled.momentumLoads()==originalResult.momentumLoads() && disabled.momentumRates()==originalResult.momentumRates() && disabled.dissipationSamples().empty(),"zero strength changed center path");
    const auto stabilized=CchKepConvection::buildControlled(*mesh,m,{});
    auto foreignSource=m.sources();foreignSource.advector={};
    reject([&]{stabilized.apply(*mesh,u,boundary,foreignSource,m.timeLayer());},"different actual advector identity");
    auto foreignClock=m.timeLayer();++foreignClock.step;
    reject([&]{stabilized.apply(*mesh,u,boundary,m.sources(),foreignClock);},"different full time layer");
    reject([&]{stabilized.apply(*mesh,u,bindBoundary(*mesh,other,t),m.sources(),m.timeLayer());},"different prepared material owner");
    reject([&]{stabilized.apply(*mesh,values(*mesh,metrics(*mesh,false).layout(),field),boundary,m.sources(),m.timeLayer());},"different layout identity");
    const auto zeroMaterial=material(*mesh,f,values(*mesh,m.layout(),[](Vector2D){return Vector2D{};}));
    for(auto F:zeroMaterial.faceRates())check(F==0,"zero source has nonzero material F");
    controlled(*mesh,CchKepConvection::build(*mesh,zeroMaterial));
    // Finite physical power despite u*u being outside the machine range.
    // This is a numerical-range contract, not a physical density benchmark.
    const auto sf=frames(*mesh,geometry,publishAt(*mesh,s,{}),publishAt(*mesh,s,{1,0,.04,Position::Midpoint}),publishAt(*mesh,s,{1,0,.04,Position::End}),1e-290);
    const auto sm=material(*mesh,sf,u);const auto sop=CchKepConvection::build(*mesh,sm);const auto hugeField=[](Vector2D){return Vector2D{1e160,1e160};};
    const auto range=sop.apply(*mesh,values(*mesh,sm.layout(),hugeField),bindBoundary(*mesh,sm,traces(sm,hugeField)),sm.sources(),sm.timeLayer());
    std::array<long double,2> rangePower{},rangeGross{};
    for(std::size_t g=0;g<sm.faceRates().size();++g){const auto& segment=sm.geometry().segments()[g];if(inside(segment))continue;
        const auto flux=sm.faceRates()[g]*(segment.negativeVelocity>=0?1:-1);const auto work=(flux*1e160L)*.5e160L;const int a=axis(sm.layout(),owner(segment));rangePower[a]+=work;rangeGross[a]+=std::abs(work);}
    for(int a=0;a<2;++a)error(range.report().components[a].boundaryKineticPower,rangePower[a],rangeGross[a],1e-12,"finite boundary power with large velocity/small flux");
    std::cout<<"RANGE,1e160,1e-290\n";
    for(int k=0;k<3;++k){auto source=m.sources();if(k==0)source.begin={};if(k==1)source.midpoint={};if(k==2)source.end={};reject([&]{op.apply(*mesh,u,boundary,source,m.timeLayer());},"different geometry identity");}
    auto wrongAdvector=m.sources();wrongAdvector.advector={};reject([&]{op.apply(*mesh,u,boundary,wrongAdvector,m.timeLayer());},"different actual advector identity");
    auto clock=m.timeLayer();clock.begin=-.01;reject([&]{op.apply(*mesh,u,boundary,m.sources(),clock);},"different full time layer");clock=m.timeLayer();clock.step+=1;reject([&]{CchConvectionBoundary::bind(*mesh,m,t,m.sources(),clock);},"different full time layer");
    reject([&]{op.apply(*mesh,values(*mesh,metrics(*mesh,false).layout(),field),boundary,m.sources(),m.timeLayer());},"different layout identity");
    reject([&]{op.apply(*mesh,CchNativeValues(m.layout(),std::vector<double>(m.layout().velocityCount(),1e308),u.pressure()),boundary,m.sources(),m.timeLayer());},"nonfinite");
    if(std::numeric_limits<long double>::digits==std::numeric_limits<double>::digits){reject([&]{op.apply(*mesh,CchNativeValues(m.layout(),std::vector<double>(m.layout().velocityCount(),std::numeric_limits<double>::denorm_min()),u.pressure()),boundary,m.sources(),m.timeLayer());},"unrepresentable");}
    auto radial=makeMesh(4);Case rs=s;rs.radial=true;const auto rg=CchMomentumGeometry::build(*radial,metrics(*radial,true));
    const auto rf=frames(*radial,rg,publishAt(*radial,rs,{}),publishAt(*radial,rs,{1,0,.04,Position::Midpoint}),publishAt(*radial,rs,{1,0,.04,Position::End}),1);
    const auto rfield=[](Vector2D p){return Vector2D{.1*p.x,.2};};const auto ru=values(*radial,rg.layout(),rfield);const auto rm=material(*radial,rf,ru);const auto rop=CchKepConvection::build(*radial,rm);const auto rt=traces(rm,rfield);const auto rb=bindBoundary(*radial,rm,rt);
    auto wrong=ru.normalVelocity();wrong[rg.metrics().axisVelocityDofs().front()]=1;reject([&]{rop.apply(*radial,CchNativeValues(rg.layout(),wrong,ru.pressure()),rb,rm.sources(),rm.timeLayer());},"nonzero radial axis velocity");
    for(std::size_t g=0;g<rt.size();++g){const auto& segment=rm.geometry().segments()[g];if(!inside(segment) && segment.normalAxis==0 && segment.begin.x==0 && axis(rm.layout(),owner(segment))==0){bad=rt;bad[g]=1;reject([&]{bindBoundary(*radial,rm,bad);},"nonzero radial axis boundary trace");break;}}
    const auto narrow=CchMomentumGeometry::build(*radial,metrics(*radial,true,3));
    const auto nf=frames(*radial,narrow,publishAt(*radial,rs,{}),publishAt(*radial,rs,{1,0,.04,Position::Midpoint}),publishAt(*radial,rs,{1,0,.04,Position::End}),1);
    const auto nm=material(*radial,nf,values(*radial,narrow.layout(),rfield));reject([&]{CchKepConvection::build(*radial,nm);},"unidentifiable affine common trace");
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
    const auto fg=CchMomentumGeometry::build(*fine,metrics(*fine,rz));const auto fm=material(*fine,frames(*fine,fg,a,mid,end,1e6),values(*fine,fg.layout(),field));const auto op=CchKepConvection::build(*fine,fm);constant(*fine,op,v,true);++regrids;
    controlled(*fine,op);
    config.volume.symmetryPlanes=s.symmetry(*coarse);const auto projected=transferHermiteGalsFirstJetMesh(*fine,end,*coarse,HermiteGalsMeshTransferKind::FirstJetProjection,config);const auto b=projected.publication();check(b.state().timeLayer().time()==.1,"coarsening changed physical time");
    const auto fineU=values(*fine,fg.layout(),field);const auto fineBoundary=bindBoundary(*fine,fm,traces(fm,field));
    reject([&]{op.apply(*coarse,fineU,fineBoundary,fm.sources(),fm.timeLayer());},"NativeMeshIdentity");const auto cm=advanceFirstJet(*coarse,b,s,{2,.1,.2,Position::Midpoint},stage),ce=advanceFirstJet(*coarse,b,s,{2,.1,.2,Position::End},stage);
    const auto cg=CchMomentumGeometry::build(*coarse,metrics(*coarse,rz));const auto cmat=material(*coarse,frames(*coarse,cg,b,cm,ce,1e6),values(*coarse,cg.layout(),field));
    check(cmat.sources().begin==b.state().identity() && cmat.timeLayer().step==2 && !cmat.layout().sameIdentity(fm.layout()),"coarsening retagged old material interval/layout");
    constant(*coarse,CchKepConvection::build(*coarse,cmat),v,true);++regrids;std::cout<<"REGRID,"<<rz<<','<<cmat.timeLayer().step<<','<<cmat.timeLayer().begin<<','<<cmat.timeLayer().end<<'\n';
    controlled(*coarse,CchKepConvection::build(*coarse,cmat));
}
void boundedPorts(bool rz,bool amr,double rho)
{
    auto mesh=makeMesh(8,amr?1:0);Case s;s.pureFirstJet=true;s.radial=rz;s.planar=true;s.planeOffset=2;
    const auto g=CchMomentumGeometry::build(*mesh,metrics(*mesh,rz));const auto& l=g.layout();
    const auto f=frames(*mesh,g,publishAt(*mesh,s,{}),publishAt(*mesh,s,{1,0,.04,Position::Midpoint}),publishAt(*mesh,s,{1,0,.04,Position::End}),rho);
    const auto m=material(*mesh,f,values(*mesh,l,[rz](Vector2D){return rz?Vector2D{0,1}:Vector2D{1,0};}));const auto saved=m.faceRates();
    std::map<std::tuple<int,int,double,double>,std::vector<int>> groups;
    for(std::size_t k=0;k<g.segments().size();++k){const auto& p=g.segments()[k];if(inside(p) || m.faceRates()[k]==0)continue;
        groups[{axis(l,owner(p)),p.normalAxis,p.normalAxis==0?p.begin.y:p.begin.x,p.normalAxis==0?p.end.y:p.end.x}].push_back(static_cast<int>(k));}
    std::vector<CchConvectionPeriodicPair> ports;for(const auto& group:groups){check(group.second.size()==2,"test unmatched port");const int a=group.second[0],b=group.second[1];ports.push_back({a,b,g.segments()[a].begin-g.segments()[b].begin});}
    const auto op=CchKepConvection::build(*mesh,m,ports);const auto& report=op.momentReport();
    check(report.boundedRows>0 && report.coefficientRatio<=1+1e-12 && report.quadraticBudgetRatio<=1+1e-10,"bounded native certificate");
    auto field=values(*mesh,l,[rz](Vector2D p){const double t=static_cast<double>(2*pi)*(rz?p.y:p.x);return rz?Vector2D{p.x*std::sin(t),.4+std::cos(t)}:Vector2D{std::cos(t),std::sin(t)};});
    auto u=field.normalVelocity();for(int i=0;i<l.velocityCount();++i)u[i]=u[op.representatives()[i]];field=CchNativeValues(l,u,field.pressure());
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
    auto missing=ports;missing.pop_back();reject([&]{CchKepConvection::build(*mesh,m,missing);},"unpaired nonzero material boundary rate");
    auto wrong=ports;wrong.front().translation=wrong.front().translation*.5;reject([&]{CchKepConvection::build(*mesh,m,wrong);},"periodic segment begin mismatch");
    check(&op.material().faceRates()==&m.faceRates() && saved==m.faceRates(),"bounded solve changed shared material flux");
    std::cout<<"BOUNDED_PORTS,"<<rz<<','<<amr<<','<<rho<<','<<report.boundedRows<<','<<report.iterations<<','<<report.coefficientRatio<<','<<maxError<<std::endl;
    controlled(*mesh,op,true,true);
}
template<class F>void migrationCheck(const std::string& name,F f)
{try{f();}catch(const std::exception& e){++migrationFailures;std::cout<<"MIGRATION_FAILURE,"<<name<<','<<e.what()<<std::endl;}}
}
int main()
{
    try{
        std::cout<<std::setprecision(17);for(bool rz:{false,true})for(int amr:{0,1}){const auto name=std::to_string(rz)+","+std::to_string(amr);for(int selection=0;selection<4;++selection) migrationCheck("algebra,"+name+","+std::to_string(selection),[&]{algebra(rz,amr,selection);});
            for(double rho:{1.,1000.,1e6})migrationCheck("bounded,"+name+","+std::to_string(rho),[&]{boundedPorts(rz,amr,rho);});
            for(bool planar:{false,true})migrationCheck("crossing,"+name+","+std::to_string(planar),[&]{crossing(rz,amr,planar);});migrationCheck("closed,"+name,[&]{closed(rz,amr);});}migrationCheck("guards",[&]{guards();});migrationCheck("regrid,0",[&]{regrid(false);});migrationCheck("regrid,1",[&]{regrid(true);});
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
        if(migrationFailures==0)check(controlledCases==57 && controlledActive>0 && controlledQuiet>0 && controlledCoarseFine>0 && controlledPeakDefect>1e-5,"incomplete controlled dissipation coverage");
        std::cout<<"MIGRATION_RESULT,"<<migrationFailures<<','<<(migrationFailures?"FAIL":"PASS")<<'\n';return migrationFailures?1:0;
    }catch(const std::exception& e){std::cerr<<"H10c FAIL: "<<e.what()<<'\n';return 1;}
}
