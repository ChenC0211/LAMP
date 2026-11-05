#include "Modules/Hydrodynamics/HydroInterfaceNormalLoad.h"
#include "Modules/Hydrodynamics/HydroMomentumGeometry.h"
#include "Modules/Hydrodynamics/HydroPhaseViscosity.h"
#include "Modules/Hydrodynamics/HydroSinglePhaseSaddle.h"
#include "interface_first_jet_fixture.h"
#include <map>
#include <numeric>
namespace {
using namespace LAMP::Hydrodynamics;
using Matrix=HydroSaddleVelocityBlock::Matrix;
int scenes=0,solves=0,guards=0,negativeControls=0,coarseFaces=0,refinements=0;
double peakLoad=0,peakVelocity=0,peakPressure=0,peakBalance=0,peakContinuity=0;
double relative(double e,double s){check(std::isfinite(e)&&std::isfinite(s)&&s>=0,"invalid independent scale");if(s==0){check(e==0,"nonzero unscaled error");return 0;}return std::abs(e)/s;}
template<class F>void reject(F f,const char* message){rejects(f,message);++guards;}
HydroPhysicalDivergence metricsFor(const QuadtreeMesh& mesh,bool rz,int selected=0){
    FaceSegmentGraph graph;graph.build(mesh,[selected](const QuadtreeMesh&,const Cell& cell){const auto p=cell.centerCoord();
        return selected==0 || (selected==1?!(p.x>.375&&p.x<.625&&p.y>.375&&p.y<.625):selected==2?p.x<.375||p.x>.625:p.x>.25);});
    return HydroPhysicalDivergence::build(mesh,HydroIntegralDivergence::build(mesh,HydroDofLayout::build(mesh,NativeLeafDomain::build(mesh,graph,rz))));
}
HermiteGalsPublication pure(const QuadtreeMesh& mesh,bool rz,double sign){
    const Profile f=[sign](Vector2D){return Jet{sign,0,0,0};};
    const auto state=stateFor(mesh,f,rz,{},[](int,Vector2D,Vector2D){return 0.;},{},HermiteGalsTransport::FirstJet);
    HermiteGalsCandidate candidate(mesh,{},state.primary());const auto mixed=closeHermiteGalsMixedAmr(mesh,candidate);
    const auto id=candidate.identity();return publishHermiteGalsFirstJetGeometry(mesh,std::move(candidate),mixed,{id,centers(mesh,f)});
}
struct Reference {
    std::vector<int> minus,plus,component;std::vector<double> area,volume,alpha;int components=0;
    Reference(const QuadtreeMesh& mesh,const HydroPhaseGeometry& phase){
        const auto& l=phase.layout();const int np=l.pressureCount(),nu=l.velocityCount();
        component.resize(np);std::iota(component.begin(),component.end(),0);
        const auto root=[&](int p){while(component[p]!=p)p=component[p];return p;};
        std::map<int,long double> liquid;for(const auto& cell:phase.publication().volumes().cells())liquid[cell.cell]=cell.rawLiquid;
        for(int p=0;p<np;++p){const auto& c=mesh.cells(mesh.cellIndexToLabel[l.nativeCell(p)]);
            const auto a=mesh.verts(c.cellToVert[0]).coord(),b=mesh.verts(c.cellToVert[2]).coord();
            const double v=(b.x-a.x)*(b.y-a.y)*(l.domain().axisymmetric()?double(pi)*(b.x+a.x):1);volume.push_back(v);
            close(phase.metrics().cellVolumes()[p],v,1e-14*v,"independent physical volume");
            close(phase.cellVolumes()[p].liquid,liquid.at(l.nativeCell(p)),1e-12*v,"same publication volume");
            alpha.push_back(double(phase.cellVolumes()[p].liquid/v));close(phase.liquidFractions()[p],alpha.back(),1e-14,"source alpha");}
        for(int u=0;u<nu;++u){const auto f=l.domain().face(mesh,l.nativeFace(u));const auto& e=f.horizontal?mesh.hEdges(f.edgeLabel):mesh.vEdges(f.edgeLabel);
            const auto a=mesh.verts(e.edgeToVert[0]).coord(),b=mesh.verts(e.edgeToVert[1]).coord();const double length=f.horizontal?std::abs(b.x-a.x):std::abs(b.y-a.y);
            area.push_back(length*(l.domain().axisymmetric()?double(2*pi)*(f.horizontal?(a.x+b.x)/2:a.x):1));
            close(area.back(),phase.metrics().faceMeasures()[u],1e-14*area.back(),"independent face area");
            minus.push_back(f.negativeCell<0?-1:l.pressureDof(f.negativeCell));plus.push_back(f.positiveCell<0?-1:l.pressureDof(f.positiveCell));
            if(minus.back()>=0&&plus.back()>=0){component[root(plus.back())]=root(minus.back());
                if(mesh.cells(mesh.cellIndexToLabel[f.negativeCell]).level()!=mesh.cells(mesh.cellIndexToLabel[f.positiveCell]).level())++coarseFaces;}}
        std::map<int,int> numbered;for(int p=0;p<np;++p){const int r=root(p);if(!numbered.count(r))numbered[r]=components++;}
        std::vector<int> final(np);for(int p=0;p<np;++p)final[p]=numbered[root(p)];component=std::move(final);
    }
    std::vector<double> gradient(const std::vector<double>& p)const{std::vector<double> f(area.size());for(std::size_t u=0;u<f.size();++u)f[u]=area[u]*((plus[u]<0?0:p[plus[u]])-(minus[u]<0?0:p[minus[u]]));return f;}
};
void inspect(const QuadtreeMesh& mesh,const HydroPhaseGeometry& phase,double density,double mu){
    ++scenes;const auto& l=phase.layout();const Reference ref(mesh,phase);
    check(phase.publication().state().primary().phaseVolumes.empty(),"Interface inventory introduced");
    const auto pressure=HydroPressureLoad::build(phase.metrics().divergence());
    const auto mass=HydroMomentumMass::build(mesh,HydroMomentumGeometry::build(mesh,phase.metrics()),phase,density,1);
    const auto viscosity=HydroPhaseViscosity::build(mesh,phase,HermiteGalsPhaseQuadrature::build(mesh,phase.publication()),mu,.01*mu);
    Matrix a=viscosity.viscosity().matrix();if(mu==0)check(a.nonZeros()==0,"zero viscosity must disappear");
    for(int u=0;u<l.velocityCount();++u)a.coeffRef(u,u)+=mass.masses()[u]/.125;
    const auto saddle=HydroSinglePhaseSaddle::buildClosed(pressure,HydroSaddleVelocityBlock::build(l,a));
    check(saddle.pressureGaugeCount()==ref.components,"independent connected components");
    std::vector<bool> free(l.velocityCount(),false);for(int u:saddle.freeVelocityDofs())free[u]=true;
    for(int u=0;u<l.velocityCount();++u)check(free[u]==(ref.minus[u]>=0&&ref.plus[u]>=0&&ref.area[u]>0),"free face classification");
    for(double jump:{-1.75,0.,2.5}){
        const auto normal=HydroInterfaceNormalLoad::constant(mesh,pressure,phase,jump);normal.validateFor(mesh,pressure,phase.geometryIdentity(),phase.timeLayer());
        auto potential=ref.alpha;for(double& v:potential)v*=jump;const auto expected=ref.gradient(potential);
        double loadError=0,loadScale=0;for(int u=0;u<l.velocityCount();++u){loadError=std::hypot(loadError,normal.balanceLoads()[u]-expected[u]);loadScale=std::hypot(loadScale,expected[u]);}
        const double load=relative(loadError,loadScale);check(load<=1e-12,"independent normal load");peakLoad=std::max(peakLoad,load);
        auto p0=potential;for(double& v:p0)v=-v;const auto cp=pressure.apply(mesh,HydroNativeValues(l,std::vector<double>(l.velocityCount()),p0));
        double full=0;for(int u=0;u<l.velocityCount();++u)full=std::hypot(full,cp[u]+normal.balanceLoads()[u]);check(relative(full,loadScale)<=1e-12,"full row nongauge balance");
        // Virtual work from physical incidence, never C coefficients.
        std::vector<double> div(l.pressureCount());double power=0,workScale=0;
        for(int u=0;u<l.velocityCount();++u){const double v=std::sin(.31*(u+1));power+=v*normal.balanceLoads()[u];workScale+=std::abs(v*normal.balanceLoads()[u]);
            if(ref.minus[u]>=0)div[ref.minus[u]]+=ref.area[u]*v;if(ref.plus[u]>=0)div[ref.plus[u]]-=ref.area[u]*v;}
        for(int p=0;p<l.pressureCount();++p){power+=potential[p]*div[p];workScale+=std::abs(potential[p]*div[p]);}check(relative(power,workScale)<=1e-12,"pressure work identity");
        auto rhs=normal.balanceLoads();for(double& v:rhs)v=-v;std::cout<<"STATIC_BEGIN,"<<scenes<<','<<jump<<','<<density<<','<<mu<<std::endl;
        const auto result=saddle.solve(mesh,HydroSaddleRhs(l,rhs,std::vector<double>(l.velocityCount())));++solves;refinements+=result.report().linearRefinements;
        std::vector<double> sum(ref.components),volume(ref.components),pExpected(l.pressureCount());
        for(int p=0;p<l.pressureCount();++p){sum[ref.component[p]]+=ref.volume[p]*ref.alpha[p];volume[ref.component[p]]+=ref.volume[p];}
        double pressureError=0,velocityMax=0,velocityScale=0,balanceError=0,balanceScale=0;
        for(int p=0;p<l.pressureCount();++p){pExpected[p]=jump*(sum[ref.component[p]]/volume[ref.component[p]]-ref.alpha[p]);pressureError=std::max(pressureError,std::abs(result.values().pressure()[p]-pExpected[p]));}
        const auto solvedLoad=ref.gradient(result.values().pressure());
        for(int u=0;u<l.velocityCount();++u){velocityMax=std::max(velocityMax,std::abs(result.values().normalVelocity()[u]));velocityScale=std::max(velocityScale,std::abs(normal.balanceLoads()[u])/a.coeff(u,u));
            if(!free[u])check(result.values().normalVelocity()[u]==0,"boundary velocity");else{balanceError=std::hypot(balanceError,solvedLoad[u]+expected[u]);balanceScale=std::hypot(balanceScale,std::abs(solvedLoad[u])+std::abs(expected[u]));}
            if(ref.area[u]==0)check(normal.balanceLoads()[u]==0,"axis normal load");}
        const double v=relative(velocityMax,velocityScale),p=relative(pressureError,std::abs(jump)),balance=relative(balanceError,balanceScale);
        check(v<=1e-11,"static velocity");check(p<=1e-10,"physical pressure jump and volume gauge");check(balance<=1e-11,"free stress balance");
        peakVelocity=std::max(peakVelocity,v);peakPressure=std::max(peakPressure,p);peakBalance=std::max(peakBalance,balance);peakContinuity=std::max(peakContinuity,result.report().maxCellContinuityDefect);
        if(jump!=0&&phase.report().cutCells>0){double wrongPressure=0;for(int k=0;k<l.pressureCount();++k)wrongPressure=std::max(wrongPressure,std::abs(-result.values().pressure()[k]-pExpected[k]));
            check(wrongPressure/std::abs(jump)>1e-10,"sign negative control invisible");++negativeControls;
            auto permuted=potential;std::rotate(permuted.begin(),permuted.begin()+1,permuted.end());const auto wrong=ref.gradient(permuted);double error=0;
            for(int u=0;u<l.velocityCount();++u)error=std::hypot(error,wrong[u]-normal.balanceLoads()[u]);check(relative(error,loadScale)>1e-12,"alpha permutation negative control invisible");++negativeControls;}
        std::cout<<"STATIC_PASS,"<<scenes<<','<<jump<<','<<v<<','<<p<<','<<balance<<','<<result.report().linearRefinements<<'\n';
    }
}
void provenance(){
    auto mesh=makeMesh(4);Case s;s.pureFirstJet=true;s.planar=true;s.planeOffset=.43;
    const auto phase=HydroPhaseGeometry::build(*mesh,metricsFor(*mesh,false),initial(*mesh,s));const auto pressure=HydroPressureLoad::build(phase.metrics().divergence());
    const auto load=HydroInterfaceNormalLoad::constant(*mesh,pressure,phase,2.5),copy=load;const auto saved=load.balanceLoads();
    check(&load.balanceLoads()==&copy.balanceLoads(),"copy duplicated immutable load");
    check(&load.geometry().publication().state().geometry().snapshot()==&phase.publication().state().geometry().snapshot(),"copied scalar source");
    const auto validate=[&](auto id,auto time){load.validateFor(*mesh,pressure,id,time);};
    reject([&]{validate(HermiteGalsIdentity{},phase.timeLayer());},"identity");reject([&]{validate(initial(*mesh,s).state().identity(),phase.timeLayer());},"identity");
    for(auto time:{HermiteGalsTimeLayer{1,0,0,Position::Initial},HermiteGalsTimeLayer{0,0,0,Position::End},HermiteGalsTimeLayer{0,0,1,Position::Initial},HermiteGalsTimeLayer{1,-1,1,Position::Midpoint}})reject([&]{validate(phase.geometryIdentity(),time);},"time");
    auto bad=phase.timeLayer();bad.end=std::numeric_limits<double>::quiet_NaN();reject([&]{validate(phase.geometryIdentity(),bad);},"time");
    for(double j:{std::numeric_limits<double>::infinity(),-std::numeric_limits<double>::infinity(),std::numeric_limits<double>::quiet_NaN()})reject([&]{HydroInterfaceNormalLoad::constant(*mesh,pressure,phase,j);},"nonfinite");
    const auto other=HydroPressureLoad::build(HydroIntegralDivergence::build(*mesh,phase.layout()));reject([&]{HydroInterfaceNormalLoad::constant(*mesh,other,phase,1);},"divergence owner");
    reject([&]{load.validateFor(*mesh,other,phase.geometryIdentity(),phase.timeLayer());},"pressure owner");const auto foreign=HydroPressureLoad::build(metricsFor(*mesh,false).divergence());reject([&]{HydroInterfaceNormalLoad::constant(*mesh,foreign,phase,1);},"layout");
    const auto kept=[&]{auto m=makeMesh(4);const auto g=HydroPhaseGeometry::build(*m,metricsFor(*m,false),initial(*m,s));return HydroInterfaceNormalLoad::constant(*m,HydroPressureLoad::build(g.metrics().divergence()),g,2.5);}();
    kept.validateFor(*mesh,kept.pressureLoad(),kept.geometry().geometryIdentity(),kept.geometry().timeLayer());
    auto target=makeMesh(4,1);reject([&]{load.validateFor(*target,pressure,phase.geometryIdentity(),phase.timeLayer());},"mesh");
    for(bool rz:{false,true}){s.radial=rz;const auto original=initial(*mesh,s);HermiteGalsPublicationConfig config;config.volume.symmetryPlanes=s.symmetry(*mesh);
        const auto transfer=transferHermiteGalsFirstJetMesh(*mesh,original,*target,HermiteGalsMeshTransferKind::ExactRefinement,config);
        check(transfer.report().refinedCells>0&&transfer.publication().state().primary().phaseVolumes.empty(),"actual FirstJet transfer");
        const auto rebuilt=HydroPhaseGeometry::build(*target,metricsFor(*target,rz),transfer.publication());inspect(*target,rebuilt,1000,0);
        const auto nextPressure=HydroPressureLoad::build(rebuilt.metrics().divergence());const auto next=HydroInterfaceNormalLoad::constant(*target,nextPressure,rebuilt,2.5);
        reject([&]{next.validateFor(*target,nextPressure,original.state().identity(),rebuilt.timeLayer());},"identity");}
    DictionaryMesh::RefineBox box;box.min={.25,.25};box.max={.75,.75};box.level=1;mesh->applyRefineBoxes({box});reject([&]{validate(phase.geometryIdentity(),phase.timeLayer());},"mesh");check(saved==load.balanceLoads(),"rejection changed retained result");
}
}
int main(){try{
    std::cout<<std::setprecision(17);
    // Reproduced zero-viscosity witness: an ordinary residual loses continuity
    // in a corner with O(1e-32) flux scale despite an accurate force solve.
    {auto mesh=makeMesh(16);Case shape;shape.pureFirstJet=true;shape.cx=.4;shape.cy=.5;
        inspect(*mesh,HydroPhaseGeometry::build(*mesh,metricsFor(*mesh,false),initial(*mesh,shape)),1,0);}
    for(bool rz:{false,true})for(int amr:{0,1,2}){auto mesh=makeMesh(4,amr);Case s;s.pureFirstJet=true;s.radial=rz;const auto metrics=metricsFor(*mesh,rz);const auto phase=HydroPhaseGeometry::build(*mesh,metrics,initial(*mesh,s));
        inspect(*mesh,phase,1,0);inspect(*mesh,phase,1000,1e-6);inspect(*mesh,phase,1e6,1);s.planar=true;s.planeOffset=.43;inspect(*mesh,HydroPhaseGeometry::build(*mesh,metrics,initial(*mesh,s)),1000,0);}
    for(bool rz:{false,true})for(int selected:{1,2}){auto mesh=makeMesh(8);Case s;s.pureFirstJet=true;s.radial=rz;s.planar=true;s.planeOffset=.43;inspect(*mesh,HydroPhaseGeometry::build(*mesh,metricsFor(*mesh,rz,selected),initial(*mesh,s)),1000,1e-3);}
    for(bool rz:{false,true})for(double sign:{-1.,1.}){auto mesh=makeMesh(4,1);inspect(*mesh,HydroPhaseGeometry::build(*mesh,metricsFor(*mesh,rz),pure(*mesh,rz,sign)),1e6,0);}
    // The Interface scalar retains axis support; Hydro selects a genuine annular domain.
    auto annular=makeMesh(4);Case plane;plane.pureFirstJet=true;plane.radial=true;plane.planar=true;plane.planeOffset=.43;
    inspect(*annular,HydroPhaseGeometry::build(*annular,metricsFor(*annular,true,3),initial(*annular,plane)),1000,1e-6);
    provenance();check(coarseFaces>0&&refinements>0&&negativeControls>0,"coverage witness missing");
    std::cout<<"H12_SUMMARY,"<<scenes<<','<<solves<<','<<guards<<','<<negativeControls<<','<<coarseFaces<<','<<refinements<<','<<peakLoad<<','<<peakVelocity<<','<<peakPressure<<','<<peakBalance<<','<<peakContinuity<<",PASS\n";return 0;
}catch(const std::exception& e){std::cerr<<"H12 constant normal stress FAIL: "<<e.what()<<'\n';return 1;}}
