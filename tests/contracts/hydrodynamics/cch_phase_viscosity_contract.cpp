#include "Modules/Hydrodynamics/CchPhaseViscosity.h" // Public headers stand alone.
#include "Modules/Hydrodynamics/CchSinglePhaseSaddle.h"
#include "interface_first_jet_fixture.h"
#include "Core/Mesh/CartesianMeasure.h"
#include <Eigen/Eigenvalues>
#include <optional>

namespace
{
using namespace Quantum::Hydrodynamics;
int scenarios=0,pointChecks=0,spectra=0,scalings=0,saddles=0,rejected=0,mms=0;
long double peakVolume=0,peakMoment=0,peakEnergy=0,peakIdentity=0,peakSymmetry=0,peakScaling=0,peakSaddle=0;
static_assert(!std::is_default_constructible_v<HermiteGalsPhaseQuadrature>);
static_assert(!std::is_default_constructible_v<CchPhaseViscosity>);
struct Gauss
{
    std::vector<long double> x,w;
    explicit Gauss(int n):x(n),w(n) {
        for(int i=0;i<(n+1)/2;++i) {
            long double z=std::cos(pi*(i+.75L)/(n+.5L)),derivative=0;
            for(int k=0;k<30;++k) {
                long double p=1,previous=0;
                for(int j=1;j<=n;++j){const long double saved=p;p=((2*j-1)*z*p-(j-1)*previous)/j;previous=saved;}
                derivative=n*(z*p-previous)/(z*z-1);const long double next=z-p/derivative;
                if(next==z)break;z=next;
            }
            const long double weight=2/((1-z*z)*derivative*derivative);
            x[i]=-z;x[n-1-i]=z;w[i]=w[n-1-i]=weight;
        }
    }
};
const Gauss outerReference(24),innerReference(12);
// Independent continuum disk/half-plane slicing. Angular endpoint mapping
// removes the disk square-root endpoints; no production scalar/roots/weights.
template<class F>long double reference(Vector2D a,Vector2D b,const Case& s,F f)
{
    std::vector<long double> cuts{a.x,b.x};
    const long double cx=s.radial?0:s.cx,cy=s.cy,r=s.radius;
    const auto add=[&](long double x){if(x>a.x && x<b.x)cuts.push_back(x);};
    if(s.planar){if(!s.axial())add(s.planeOffset);}
    else {
        add(cx-r);add(cx+r);
        for(long double y:{static_cast<long double>(a.y),static_cast<long double>(b.y)}) {
            const long double d=y-cy;if(std::abs(d)<r){const auto width=std::sqrt(r*r-d*d);add(cx-width);add(cx+width);}
        }
    }
    std::sort(cuts.begin(),cuts.end());cuts.erase(std::unique(cuts.begin(),cuts.end()),cuts.end());
    long double result=0;
    for(std::size_t k=1;k<cuts.size();++k)for(std::size_t i=0;i<outerReference.x.size();++i) {
        const long double t=(1+outerReference.x[i])/2,angle=pi*t/2;
        const long double x=cuts[k-1]+(cuts[k]-cuts[k-1])*std::sin(angle)*std::sin(angle);
        const long double wx=(cuts[k]-cuts[k-1])*pi*std::sin(pi*t)*outerReference.w[i]/4;
        long double liquidLo=a.y,liquidHi=b.y;
        if(s.planar) {
            if(s.axial()) {if(s.planeSign>0)liquidHi=std::min(liquidHi,static_cast<long double>(s.planeOffset));else liquidLo=std::max(liquidLo,static_cast<long double>(s.planeOffset));}
            else if(s.planeSign*(x-s.planeOffset)>=0)liquidHi=liquidLo;
        }else {
            const long double q=r*r-(x-cx)*(x-cx);
            if(q<=0)liquidHi=liquidLo;
            else {const auto height=std::sqrt(q);liquidLo=std::max(liquidLo,cy-height);liquidHi=std::min(liquidHi,cy+height);if(liquidHi<liquidLo)liquidHi=liquidLo;}
        }
        const auto interval=[&](long double lo,long double hi,bool liquid) {
            if(hi<=lo)return;
            for(std::size_t j=0;j<innerReference.x.size();++j) {
                const long double y=(lo+hi)/2+(hi-lo)*innerReference.x[j]/2;
                const long double weight=wx*(hi-lo)*innerReference.w[j]/2*(s.radial?CartesianMeasure::twoPi*x:1.L);
                result+=weight*f(Vector2D{static_cast<double>(x),static_cast<double>(y)},liquid);
            }
        };
        if(liquidHi>liquidLo){interval(a.y,liquidLo,false);interval(liquidLo,liquidHi,true);interval(liquidHi,b.y,false);}
        else interval(a.y,b.y,false);
    }
    return result;
}
CchPhaseGeometry phaseFor(const QuadtreeMesh& mesh,const Case& s,bool subset=false)
{
    FaceSegmentGraph graph;graph.build(mesh,[subset](const QuadtreeMesh&,const Cell& c){return !subset || c.centerCoord().x<.5;});
    const auto layout=CchDofLayout::build(mesh,NativeLeafDomain::build(mesh,graph,s.radial));
    return CchPhaseGeometry::build(mesh,CchPhysicalDivergence::build(mesh,CchIntegralDivergence::build(mesh,layout)),initial(mesh,s));
}
struct FluidField{Vector2D u;std::array<double,4> strain;};
FluidField polynomial(Vector2D p,bool rz)
{
    const double x=p.x,y=p.y;
    if(rz)return {{x*(1+2*y),3+x*x+y+2*y*y},{1+2*y,1+4*y,1+2*y,4*x/std::sqrt(2.)}};
    return {{1+x+2*y+3*x*x+4*x*y+5*y*y,-1+2*x-y+4*x*x+3*x*y+2*y*y},
        {1+6*x+4*y,-1+3*x+4*y,(4+12*x+13*y)/std::sqrt(2.),0}};
}
FluidField smooth(Vector2D p,bool rz)
{
    const double x=p.x,y=p.y,k=static_cast<double>(pi),sx=std::sin(k*x),cx=std::cos(k*x),sy=std::sin(k*y),cy=std::cos(k*y);
    if(rz)return {{x*sy,cx*sy},{sy,k*cx*cy,sy,(k*x*cy-k*sx*sy)/std::sqrt(2.)}};
    return {{sx*cy,cx*sy},{k*cx*cy,k*cx*cy,-std::sqrt(2.)*k*sx*sy,0}};
}
template<class F>CchNativeValues values(const QuadtreeMesh& mesh,const CchDofLayout& layout,F field)
{
    std::vector<double> u;
    for(int i=0;i<layout.velocityCount();++i){const auto face=layout.domain().face(mesh,layout.nativeFace(i));const auto v=field(face.center).u;u.push_back(face.horizontal?v.y:v.x);}
    return CchNativeValues(layout,u,std::vector<double>(layout.pressureCount()));
}
template<class F>long double referenceEnergy(const QuadtreeMesh& mesh,const CchPhaseGeometry& geometry,const Case& s,double liquid,double gas,F field)
{
    long double energy=0;
    for(int p=0;p<geometry.layout().pressureCount();++p) {
        const auto bounds=geometry.publication().state().geometry().snapshot().nativeCellBounds(geometry.layout().nativeCell(p));
        energy+=reference(bounds[0],bounds[1],s,[&](Vector2D q,bool l){const auto f=field(q);long double sum=0;for(double strain:f.strain)sum+=strain*strain;return 2*(l?liquid:gas)*sum;});
    }
    return energy;
}
long double gate(long double actual,long double expected,long double scale,long double tolerance,const std::string& why)
{return close(actual,expected,tolerance*scale,why)/scale;}
void inspect(const QuadtreeMesh& mesh,const Case& s,const CchPhaseGeometry& geometry,const HermiteGalsPhaseQuadrature& quad)
{
    const auto& scalar=geometry.publication().state().geometry().snapshot();
    std::vector<std::array<long double,6>> moments(mesh.cellIndexNum());
    for(const auto& q:quad.points()) {
        check(q.physicalWeight>0 && std::isfinite(q.physicalWeight),"nonpositive phase quadrature weight");
        check(scalar.nativeQuarticBoxSign(q.cell,q.point,q.point)==(q.liquid?-1:1),"quadrature lies in wrong phase");++pointChecks;
        if(!q.liquid)continue;
        const auto b=scalar.nativeCellBounds(q.cell);
        const long double x=(q.point.x-b[0].x)/(b[1].x-b[0].x),y=(q.point.y-b[0].y)/(b[1].y-b[0].y),w=q.physicalWeight;
        const std::array<long double,6> m{{1,x,y,x*x,x*y,y*y}};
        for(int i=0;i<6;++i)moments[q.cell][i]+=w*m[i];
    }
    for(const auto& cell:quad.cells()) {
        const auto b=scalar.nativeCellBounds(cell.cell);const int p=geometry.layout().pressureDof(cell.cell);
        if(p<0)continue;const long double cap=geometry.metrics().cellVolumes()[p];
        peakVolume=std::max(peakVolume,gate(cell.liquid+cell.gas,cap,cap,2e-8,"phase quadrature capacity"));
        const long double analytic=s.planar?reference(b[0],b[1],s,[](Vector2D,bool l){return l?1.L:0.L;}):circleCell(b[0],b[1],{s.radial?0:s.cx,s.cy},s.radius,s.radial);
        peakVolume=std::max(peakVolume,gate(cell.liquid,analytic,cap,2e-8,"analytic phase volume"));
        for(int i=0;i<6;++i) {
            const long double exact=reference(b[0],b[1],s,[&](Vector2D q,bool l){if(!l)return 0.L;
                const long double x=(q.x-b[0].x)/(b[1].x-b[0].x),y=(q.y-b[0].y)/(b[1].y-b[0].y);
                return std::array<long double,6>{{1,x,y,x*x,x*y,y*y}}[i];});
            peakMoment=std::max(peakMoment,gate(moments[cell.cell][i],exact,cap,2e-8,"analytic phase moment"));
        }
    }
    check(geometry.publication().state().primary().phaseVolumes.empty(),"restored material inventory");++scenarios;
    std::cout<<"QUADRATURE,"<<s.radial<<','<<s.amr<<','<<s.planar<<','<<s.nonpolynomial<<','<<quad.points().size()<<','
        <<quad.report().sliceQueries<<','<<quad.report().maxRelativeVolumeDifference<<','<<quad.report().maxRelativeMomentEstimate<<'\n';
}
void saddle(const QuadtreeMesh& mesh,const CchPhaseViscosity& phase,const CchNativeValues& v)
{
    const auto& layout=phase.geometry().layout();CchEnergyViscosity::Matrix a=phase.viscosity().matrix();
    for(int i=0;i<layout.velocityCount();++i)a.coeffRef(i,i)+=1;
    const auto pressure=CchPressureLoad::build(phase.geometry().metrics().divergence());
    const auto solver=CchSinglePhaseSaddle::buildClosed(pressure,CchSaddleVelocityBlock::build(layout,a));
    const auto div=phase.geometry().metrics().divergence().apply(mesh,v);
    for(double x:div)check(std::abs(x)<1e-12,"saddle manufactured velocity not solenoidal");
    std::vector<double> p(layout.pressureCount());
    for(int c=0;c<layout.pressureCount();++c)
        p[c]=mesh.cells(mesh.cellIndexToLabel[layout.nativeCell(c)]).centerCoord().y-.5;
    const Eigen::Map<const Eigen::VectorXd> u(v.normalVelocity().data(),layout.velocityCount()),ep(p.data(),layout.pressureCount());
    const Eigen::VectorXd rhs=a*u+pressure.matrix()*ep;
    auto prescribed=v.normalVelocity();
    for(int free:solver.freeVelocityDofs())prescribed[free]=0;
    const auto result=solver.solve(mesh,CchSaddleRhs(layout,std::vector<double>(rhs.data(),rhs.data()+rhs.size()),prescribed));
    for(int i=0;i<layout.velocityCount();++i)peakSaddle=std::max(peakSaddle,gate(result.values().normalVelocity()[i],u[i],1,1e-10,"phase viscosity saddle"));++saddles;
    for(int i=0;i<layout.pressureCount();++i)peakSaddle=std::max(peakSaddle,gate(result.values().pressure()[i],p[i],1,1e-10,"phase viscosity saddle pressure"));
}
void algebra(const QuadtreeMesh& mesh,const Case& s,const CchPhaseGeometry& geometry,const HermiteGalsPhaseQuadrature& quad,bool spectrum)
{
    const auto field=[&](Vector2D p){return polynomial(p,s.radial);};const auto v=values(mesh,geometry.layout(),field);
    const auto phase=CchPhaseViscosity::build(mesh,geometry,quad,1,.01);
    const long double referenceValue=referenceEnergy(mesh,geometry,s,1,.01,field);
    const double energy=phase.dissipation(mesh,v,geometry.geometryIdentity(),geometry.timeLayer());
    peakEnergy=std::max(peakEnergy,gate(energy,referenceValue,std::max(1.L,referenceValue),5e-8,"independent phase energy"));
    const auto& k=phase.viscosity().matrix();const Eigen::Map<const Eigen::VectorXd> u(v.normalVelocity().data(),geometry.layout().velocityCount());
    peakIdentity=std::max(peakIdentity,gate(u.dot(k*u),energy,std::max(1.,energy),1e-10,"matrix and phase integral"));
    const CchEnergyViscosity::Matrix difference=k-CchEnergyViscosity::Matrix(k.transpose());
    peakSymmetry=std::max(peakSymmetry,gate(difference.norm(),0,std::max(1.,k.norm()),1e-12,"phase matrix symmetry"));
    if(spectrum) {
        // Remove exactly constrained radial axis entries before spectral check.
        std::vector<int> free;for(int i=0;i<geometry.layout().velocityCount();++i)if(std::find(geometry.metrics().axisVelocityDofs().begin(),geometry.metrics().axisVelocityDofs().end(),i)==geometry.metrics().axisVelocityDofs().end())free.push_back(i);
        Eigen::MatrixXd dense(free.size(),free.size());for(std::size_t i=0;i<free.size();++i)for(std::size_t j=0;j<free.size();++j)dense(i,j)=k.coeff(free[i],free[j]);
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigen(dense);check(eigen.info()==Eigen::Success,"phase spectrum solve");
        check(eigen.eigenvalues().minCoeff()>=-1e-10*std::max(1.,eigen.eigenvalues().maxCoeff()),"negative phase eigenvalue");
        std::cout<<"SPECTRUM,"<<s.radial<<','<<s.amr<<','<<eigen.eigenvalues().minCoeff()<<','<<eigen.eigenvalues().maxCoeff()<<'\n';++spectra;
        for(const auto& material:std::array<std::array<double,2>,3>{{{{1,0}},{{0,1}},{{1,1}}}}) {
            const auto one=CchPhaseViscosity::build(mesh,geometry,quad,material[0],material[1]);
            const auto expected=referenceEnergy(mesh,geometry,s,material[0],material[1],field);
            peakEnergy=std::max(peakEnergy,gate(one.dissipation(mesh,v,geometry.geometryIdentity(),geometry.timeLayer()),
                expected,std::max(1.L,expected),5e-8,"one-sided/equal phase viscosity"));++scalings;
        }
        const auto solenoidal=values(mesh,geometry.layout(),[&](Vector2D p){return s.radial?
            FluidField{{p.x,-2*p.y},{1,-2,1,0}}:FluidField{{p.x,-p.y},{1,-1,0,0}};});
        saddle(mesh,phase,solenoidal);
    }
    for(double factor:{0.,1e-6,1e-12}) {
        const auto scaled=CchPhaseViscosity::build(mesh,geometry,quad,factor,.01*factor);
        const double actual=scaled.dissipation(mesh,v,geometry.geometryIdentity(),geometry.timeLayer());
        peakScaling=std::max(peakScaling,gate(actual,energy*factor,std::max(1.L,static_cast<long double>(energy)*factor),1e-10,"phase viscosity scaling"));
        const CchEnergyViscosity::Matrix matrixError=scaled.viscosity().matrix()-factor*k;
        peakScaling=std::max(peakScaling,gate(matrixError.norm(),0,std::max(1.,factor*k.norm()),1e-10,"phase matrix scaling"));
        if(factor==0) {
            check(actual==0 && scaled.viscosity().matrix().nonZeros()==0,"zero phase viscosity retained dissipation");
            for(double value:scaled.apply(mesh,v,geometry.geometryIdentity(),geometry.timeLayer()))check(value==0,"zero phase viscosity retained load");
        }
        ++scalings;
    }
}
template<class F>void reject(F f,const std::string& message){rejects(f,message);++rejected;}
void guards()
{
    auto mesh=makeMesh(4);Case s;s.pureFirstJet=true;s.planar=true;s.planeOffset=.43;
    const auto geometry=phaseFor(*mesh,s);const auto quad=HermiteGalsPhaseQuadrature::build(*mesh,geometry.publication());
    const auto phase=CchPhaseViscosity::build(*mesh,geometry,quad,1,.01);
    const auto v=values(*mesh,geometry.layout(),[](Vector2D p){return polynomial(p,false);});
    const auto before=phase.apply(*mesh,v,geometry.geometryIdentity(),geometry.timeLayer());
    for(double bad:{-1.,std::numeric_limits<double>::quiet_NaN(),std::numeric_limits<double>::infinity()}) {
        reject([&]{CchPhaseViscosity::build(*mesh,geometry,quad,bad,1);},"negative/nonfinite phase");
        reject([&]{CchPhaseViscosity::build(*mesh,geometry,quad,1,bad);},"negative/nonfinite phase");
    }
    const auto other=phaseFor(*mesh,s);
    reject([&]{CchPhaseViscosity::build(*mesh,other,quad,1,1);},"different geometry identity");
    reject([&]{phase.validateFor(*mesh,{},geometry.timeLayer());},"different geometry identity");
    auto time=geometry.timeLayer();time.step=1;
    reject([&]{phase.validateFor(*mesh,geometry.geometryIdentity(),time);},"different full time layer");
    time=geometry.timeLayer();time.begin-=1;
    reject([&]{phase.validateFor(*mesh,geometry.geometryIdentity(),time);},"initial/restart layer must have zero interval");
    time=geometry.timeLayer();time.end+=1;
    reject([&]{phase.validateFor(*mesh,geometry.geometryIdentity(),time);},"initial/restart layer must have zero interval");
    time=geometry.timeLayer();time.begin=time.end=1;
    reject([&]{phase.validateFor(*mesh,geometry.geometryIdentity(),time);},"different full time layer");
    reject([&]{phase.apply(*mesh,values(*mesh,other.layout(),[](Vector2D p){return polynomial(p,false);}),geometry.geometryIdentity(),geometry.timeLayer());},"different layout identity");
    for(int kind=0;kind<9;++kind) {
        HermiteGalsPhaseQuadratureConfig c;
        if(kind==0)c.relativeMomentTolerance=0;if(kind==1)c.relativeVolumeTolerance=-1;
        if(kind==2)c.relativeVolumeTolerance=std::numeric_limits<double>::quiet_NaN();
        if(kind==3)c.relativeMomentTolerance=1e-3;if(kind==4)c.maximumDepth=-1;if(kind==5)c.maximumDepth=41;
        if(kind==6)c.maximumSliceQueriesPerCell=0;if(kind==7)c.maximumPoints=0;if(kind==8)c.relativeMomentTolerance=std::numeric_limits<double>::infinity();
        reject([&]{HermiteGalsPhaseQuadrature::build(*mesh,geometry.publication(),c);},"invalid quadrature");
    }
    HermiteGalsPhaseQuadratureConfig low;low.maximumSliceQueriesPerCell=1;
    reject([&]{HermiteGalsPhaseQuadrature::build(*mesh,geometry.publication(),low);},"slice work limit");
    low={};low.maximumPoints=1;
    reject([&]{HermiteGalsPhaseQuadrature::build(*mesh,geometry.publication(),low);},"point work limit");
    low={};low.maximumDepth=0;Case circle=s;circle.planar=false;const auto curved=phaseFor(*mesh,circle);
    reject([&]{HermiteGalsPhaseQuadrature::build(*mesh,curved.publication(),low);},"quadrature depth limit");
    Case historical=s;historical.pureFirstJet=false;
    reject([&]{HermiteGalsPhaseQuadrature::build(*mesh,initial(*mesh,historical));},"requires FirstJet");
    reject([&]{CchPhaseViscosity::build(*mesh,geometry,quad,std::numeric_limits<double>::max(),std::numeric_limits<double>::max());},"nonfinite");
    const auto after=phase.apply(*mesh,v,geometry.geometryIdentity(),geometry.timeLayer());check(before==after,"failure altered accepted phase operator");
    Case rz=s;rz.radial=true;const auto radial=phaseFor(*mesh,rz);
    const auto rq=HermiteGalsPhaseQuadrature::build(*mesh,radial.publication());const auto rk=CchPhaseViscosity::build(*mesh,radial,rq,0,0);
    std::vector<double> axisU(radial.layout().velocityCount());axisU[radial.metrics().axisVelocityDofs().front()]=1e-300;
    reject([&]{rk.apply(*mesh,CchNativeValues(radial.layout(),axisU,std::vector<double>(radial.layout().pressureCount())),radial.geometryIdentity(),radial.timeLayer());},"nonzero radial axis");
    auto& vertex=mesh->verts(mesh->vertIndexToLabel[0]);const auto old=vertex.coord();vertex.setCoord(old+Vector2D{.01,0});
    reject([&]{phase.validateFor(*mesh,geometry.geometryIdentity(),geometry.timeLayer());},"NativeMeshIdentity");vertex.setCoord(old);
    DictionaryMesh::RefineBox box;box.level=1;box.min={.25,.25};box.max={.75,.75};mesh->applyRefineBoxes({box});
    reject([&]{phase.validateFor(*mesh,geometry.geometryIdentity(),geometry.timeLayer());},"NativeMeshIdentity");
    const auto retained=[&]{const auto g=phaseFor(*mesh,s);const auto q=HermiteGalsPhaseQuadrature::build(*mesh,g.publication());return CchPhaseViscosity::build(*mesh,g,q,1,.01);}();
    const auto copy=retained;
    check(&copy.quadrature().publication().state().geometry().snapshot()==&retained.geometry().publication().state().geometry().snapshot(),"copied a second scalar");
    const auto rvalues=values(*mesh,copy.geometry().layout(),[](Vector2D p){return polynomial(p,false);});
    check(copy.dissipation(*mesh,rvalues,copy.geometry().geometryIdentity(),copy.geometry().timeLayer())>0,"retained phase operator lost source lifetime");
}
void convergence(bool rz,int amr)
{
    double previous=0,order=0;
    for(int n:{4,8,16}) {
        auto mesh=makeMesh(n,amr);Case s;s.pureFirstJet=true;s.radial=rz;s.amr=amr!=0;s.nonpolynomial=true;
        const auto geometry=phaseFor(*mesh,s);const auto quad=HermiteGalsPhaseQuadrature::build(*mesh,geometry.publication());
        const auto field=[&](Vector2D p){return smooth(p,rz);};const auto v=values(*mesh,geometry.layout(),field);
        const auto phase=CchPhaseViscosity::build(*mesh,geometry,quad,1,.01);
        const long double exact=referenceEnergy(*mesh,geometry,s,1,.01,field);
        const double actual=phase.dissipation(*mesh,v,geometry.geometryIdentity(),geometry.timeLayer());
        const double error=static_cast<double>(std::abs(actual-exact)/exact);
        if(n>4){check(error<previous,"phase energy MMS did not decrease");order=std::log2(previous/error);}
        std::cout<<"MMS,"<<rz<<','<<amr<<','<<n<<','<<actual<<','<<exact<<','<<error<<','<<quad.points().size()<<'\n';previous=error;++mms;
    }
    check(order>=1.8,"phase viscosity MMS order below gate");std::cout<<"ORDER,"<<rz<<','<<amr<<','<<order<<'\n';
}
}
int main()
{
    try {
        std::cout<<std::setprecision(17);
        for(bool rz:{false,true})for(int amr:{0,1}) {
            auto mesh=makeMesh(4,amr);
            for(int kind=0;kind<4;++kind) {
                Case s;s.radial=rz;s.amr=amr!=0;s.pureFirstJet=true;s.planar=kind!=0;s.planeOffset=kind==1?.43:kind==2?2:-1;
                const auto geometry=phaseFor(*mesh,s);const auto quad=HermiteGalsPhaseQuadrature::build(*mesh,geometry.publication());
                inspect(*mesh,s,geometry,quad);algebra(*mesh,s,geometry,quad,kind==0);
            }
        }
        for(bool radial:{false,true})for(bool gas:{false,true}) {
            auto mesh=makeMesh(4);Case s;s.pureFirstJet=true;s.radial=radial;s.planar=true;s.planeOffset=gas?1-std::ldexp(1.,-40):std::ldexp(1.,-40);
            const auto geometry=phaseFor(*mesh,s);const auto quad=HermiteGalsPhaseQuadrature::build(*mesh,geometry.publication());
            long double volume=0;for(const auto& q:quad.points())if(q.liquid!=gas)volume+=q.physicalWeight;
            const long double expected=std::ldexp(1.L,-40)*(radial?CartesianMeasure::twoPi/2:1.L);
            gate(volume,expected,expected,1e-5,"tiny phase volume survived");
            std::cout<<"TINY,"<<radial<<','<<gas<<','<<volume<<'\n';inspect(*mesh,s,geometry,quad);
        }
        for(bool rz:{false,true}) {
            auto mesh=makeMesh(8,1);Case s;s.pureFirstJet=true;s.radial=rz;s.amr=true;
            const auto geometry=phaseFor(*mesh,s,true);const auto quad=HermiteGalsPhaseQuadrature::build(*mesh,geometry.publication());
            check(geometry.layout().pressureCount()<static_cast<int>(quad.cells().size()),"missing explicit selected subset");
            inspect(*mesh,s,geometry,quad);algebra(*mesh,s,geometry,quad,false);
            std::cout<<"SUBSET,"<<rz<<','<<geometry.layout().pressureCount()<<','<<quad.cells().size()<<'\n';
        }
        guards();for(bool rz:{false,true})for(int amr:{0,1})convergence(rz,amr);
        std::cout<<"SUMMARY,"<<scenarios<<','<<pointChecks<<','<<spectra<<','<<scalings<<','<<saddles<<','<<rejected<<','<<mms<<','
            <<peakVolume<<','<<peakMoment<<','<<peakEnergy<<','<<peakIdentity<<','<<peakSymmetry<<','<<peakScaling<<','<<peakSaddle<<'\n';
        std::cout<<"H07b CCH phase quadrature viscosity PASS\n";
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
