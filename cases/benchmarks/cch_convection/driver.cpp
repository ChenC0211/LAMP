// Research owner: cases/benchmarks. Calls the actual H07-H10 CCH operators.
// Autonomous vector-transport experiment, not a production NS time stepper.
#include "Modules/Hydrodynamics/CchKepConvection.h"
// H10j0 observer begin
#include "Modules/Hydrodynamics/CchSpatialAdvector.h"
// H10j0 observer end
#include "Modules/Interface/HermiteGalsPhaseRegion.h"
#include "interface_first_jet_fixture.h"
#include <filesystem>
#include <fstream>
#include <map>
#include <set>
#include <Eigen/SVD>
#include <complex>
#include <tuple>

namespace
{
using namespace Quantum::Hydrodynamics;
using Field=std::vector<double>;
// Read-only CSV schema follows the complete coefficient inventory of its
// source type; zero high-order entries are still named and written.
template<class Coefficients> void writeCoefficientHeader(std::ostream& out,const char* stem)
{for(std::size_t k=0;k<std::tuple_size<Coefficients>::value;++k)out<<','<<stem<<k;}
template<class Coefficients> void writeCoefficients(std::ostream& out,const Coefficients& coefficients)
{for(auto value:coefficients)out<<','<<value;}
constexpr double duration=1;
int failed=0,windows=0;
std::ofstream summaries,history,gates;
// H10n: explicit versioned diagnostics; no numerical state.
std::ofstream legacyGates,donorPower,donorSegments,donorNodes;
int legacyFailed=0;
void legacyGate(const std::string& name,bool pass,long double value,long double limit)
{
    legacyGates<<name<<','<<(pass?"PASS":"FAIL")<<','<<value<<','<<limit<<'\n';
    std::cout<<"LEGACY_GATE,"<<name<<','<<(pass?"PASS":"FAIL")<<','<<value<<','<<limit<<std::endl;
    if(!pass)++legacyFailed;
}
void gate(const std::string& name,bool pass,long double value,long double limit)
{
    gates<<name<<','<<(pass?"PASS":"FAIL")<<','<<value<<','<<limit<<'\n';
    std::cout<<"GATE,"<<name<<','<<(pass?"PASS":"FAIL")<<','<<value<<','<<limit<<std::endl;
    if(!pass)++failed;
}
int component(const CchDofLayout& l,int u)
{return l.normalAxis(u)==CchDofLayout::NormalAxis::First?0:1;}
bool internal(const CchMomentumSegment& s)
{return s.negativeVelocity>=0 && s.positiveVelocity>=0;}
int owner(const CchMomentumSegment& s)
{return s.negativeVelocity>=0?s.negativeVelocity:s.positiveVelocity;}
long double relative(long double value,long double scale)
{return std::abs(value)/std::max(1.L,scale);}
HermiteGalsPublication publication(const QuadtreeMesh& mesh,bool rz,HermiteGalsTimeLayer clock)
{
    // Constant negative FirstJet: actual pure liquid geometry, no inventory.
    const Profile f=[](Vector2D){return Jet{-1,0,0,0};};
    Case s;s.pureFirstJet=true;s.radial=rz;
    const auto planes=s.symmetry(mesh);
    const auto state=stateFor(mesh,f,rz,clock,[](int,Vector2D,Vector2D){return 0.;},planes,HermiteGalsTransport::FirstJet);
    HermiteGalsCandidate candidate(mesh,clock,state.primary());
    const auto mixed=closeHermiteGalsMixedAmr(mesh,candidate,planes);
    HermiteGalsPublicationConfig cfg;cfg.volume.symmetryPlanes=planes;
    const auto id=candidate.identity();
    return publishHermiteGalsFirstJetGeometry(mesh,std::move(candidate),mixed,{id,centers(mesh,f)},cfg);
}
CchDualMaterial prepare(const QuadtreeMesh& mesh,bool rz)
{
    FaceSegmentGraph graph;graph.build(mesh);
    const auto l=CchDofLayout::build(mesh,NativeLeafDomain::build(mesh,graph,rz));
    const auto metrics=CchPhysicalDivergence::build(mesh,CchIntegralDivergence::build(mesh,l));
    const auto geometry=CchMomentumGeometry::build(mesh,metrics);
    const auto a=CchPhaseGeometry::build(mesh,metrics,publication(mesh,rz,{}));
    const auto m=CchPhaseGeometry::build(mesh,metrics,publication(mesh,rz,{1,0,duration,Position::Midpoint}));
    const auto b=CchPhaseGeometry::build(mesh,metrics,publication(mesh,rz,{1,0,duration,Position::End}));
    const auto ma=CchMomentumMass::build(mesh,geometry,a,1,1);
    const auto mm=CchMomentumMass::build(mesh,geometry,m,1,1);
    const auto mb=CchMomentumMass::build(mesh,geometry,b,1,1);
    Field v(l.velocityCount());for(int i=0;i<l.velocityCount();++i)v[i]=component(l,i)==(rz?1:0)?1:0;
    const CchNativeValues advector(l,v,Field(l.pressureCount()));
    // H10k migration begin
    const auto actual=CchAdvectorReconstruction::prepare(mesh,geometry).evaluate(mesh,advector,m.geometryIdentity(),m.timeLayer());
    const auto rate=CchPhaseRate::prepareSpatial(mesh,m,actual);
    // H10k migration end
    return CchDualMaterial::prepare(mesh,CchMaterialTransport::prepare(mesh,a,rate,b,1,1),ma,mm,mb);
}
std::vector<Vector2D> materialAdvector(const QuadtreeMesh& mesh,const CchDualMaterial& material)
{
    // H10k migration begin
    material.material().validateFor(mesh,material.geometry(),material.sources(),material.timeLayer());
    const auto& actual=material.material().midpointRate().spatialAdvector();std::vector<Vector2D> result;
    for(int i=0;i<material.layout().velocityCount();++i)result.push_back(actual.nativeVelocity(i));
    // H10k migration end
    return result;
}
struct Study
{
    std::unique_ptr<QuadtreeMesh> mesh;
    CchKepConvection op;
    int n;bool rz,amr;
    std::vector<int> paired;
    std::vector<Vector2D> points;
    std::vector<bool> fixedAxis;
    long double maxPairError=0,maxBoundaryPower=0,maxMaterialPower=0,maxAxisLoad=0;
    long double maxDefect=0,maxMassChange=0,maxGeometryChange=0;
    int coarseFine=0;
    static CchKepConvection build(const QuadtreeMesh& mesh,bool radial,bool closed)
    {
        const auto m=prepare(mesh,radial);if(!closed)return CchKepConvection::build(mesh,m);
        std::map<std::tuple<int,int,double,double>,std::vector<int>> groups;const auto& segments=m.geometry().segments();
        for(std::size_t g=0;g<segments.size();++g){const auto& s=segments[g];if(internal(s) || m.faceRates()[g]==0)continue;
            groups[{component(m.layout(),owner(s)),s.normalAxis,s.normalAxis==0?s.begin.y:s.begin.x,s.normalAxis==0?s.end.y:s.end.x}].push_back(static_cast<int>(g));}
        std::vector<CchConvectionPeriodicPair> pairs;
        for(const auto& entry:groups){check(entry.second.size()==2,"declared periodic interval is not paired");const int g=entry.second[0],h=entry.second[1];pairs.push_back({g,h,segments[g].begin-segments[h].begin});}
        return CchKepConvection::build(mesh,m,pairs);
    }
    Study(int resolution,bool radial,bool refined,bool closed=true)
        :mesh(makeMesh(resolution,refined?1:0)),op(build(*mesh,radial,closed)),
        n(resolution),rz(radial),amr(refined)
    {
        const auto& l=op.layout();const auto& m=op.material();
        for(int i=0;i<l.velocityCount();++i){
            const auto p=l.domain().face(*mesh,l.nativeFace(i)).center;points.push_back(p);
            fixedAxis.push_back(rz && component(l,i)==0 && p.x==0);
            maxDefect=std::max(maxDefect,relative(m.cells()[i].defectRate,m.midpointMass().masses()[i]));
            maxMassChange=std::max(maxMassChange,relative(m.cells()[i].massChange,m.midpointMass().masses()[i]));
            maxGeometryChange=std::max(maxGeometryChange,relative(m.cells()[i].geometryTransferChange,m.midpointMass().masses()[i]));
        }
        const auto& segments=m.geometry().segments();paired.assign(segments.size(),-1);
        using Key=std::tuple<int,int,double,double>;
        std::map<Key,std::vector<int>> boundaries;
        for(std::size_t g=0;g<segments.size();++g){
            const auto& s=segments[g];
            if(internal(s)){
                int low=100,high=-1;
                for(int u:{s.negativeVelocity,s.positiveVelocity})for(const auto& p:m.geometry().volumes()[u].overlaps){
                    const int level=mesh->cells(mesh->cellIndexToLabel[l.nativeCell(p.pressure)]).level();low=std::min(low,level);high=std::max(high,level);}
                if(low<high)++coarseFine;continue;
            }
            if(m.faceRates()[g]==0)continue;
            const double lo=s.normalAxis==0?s.begin.y:s.begin.x,hi=s.normalAxis==0?s.end.y:s.end.x;
            boundaries[{component(l,owner(s)),s.normalAxis,lo,hi}].push_back(static_cast<int>(g));
        }
        for(const auto& entry:boundaries){
            check(entry.second.size()==2,"periodic boundary interval is not uniquely paired");
            const int g=entry.second[0],h=entry.second[1];const auto& a=segments[g];const auto& b=segments[h];
            check((a.negativeVelocity>=0)!=(b.negativeVelocity>=0),"periodic pair has same orientation");
            const double pa=a.normalAxis==0?a.begin.x:a.begin.y,pb=b.normalAxis==0?b.begin.x:b.begin.y;
            check(std::abs(pa-pb)==1,"periodic pair does not join opposite exterior planes");
            maxPairError=std::max({maxPairError,relative(a.area-b.area,std::abs(a.area)+std::abs(b.area)),
                relative(m.faceRates()[g]-m.faceRates()[h],std::abs(m.faceRates()[g])+std::abs(m.faceRates()[h]))});
            paired[g]=h;paired[h]=g;
        }
    }
    Field exact(double time,bool derivative=false)const
    {
        Field w(points.size());for(std::size_t i=0;i<w.size();++i){
            const int c=component(op.layout(),static_cast<int>(i));const auto p=points[i];
            const double theta=static_cast<double>(2*pi)*((rz?p.y:p.x)-time);
            if(rz)w[i]=c==0?(derivative?-static_cast<double>(2*pi)*p.x*std::cos(theta):p.x*std::sin(theta)):
                (derivative?static_cast<double>(2*pi)*std::sin(theta):.4+std::cos(theta));
            else w[i]=c==0?(derivative?static_cast<double>(2*pi)*std::sin(theta):std::cos(theta)):
                (derivative?-static_cast<double>(2*pi)*std::cos(theta):std::sin(theta));
        }return w;
    }
    std::pair<long double,long double> error(const Field& w,const Field& reference)const
    {
        long double numerator=0,denominator=0,linf=0;
        for(std::size_t i=0;i<w.size();++i){const long double d=static_cast<long double>(w[i])-reference[i],mass=op.material().midpointMass().masses()[i];
            numerator+=mass*d*d;denominator+=mass*reference[i]*reference[i];linf=std::max(linf,std::abs(d));}
        return {std::sqrt(numerator/std::max(denominator,1e-30L)),linf};
    }
    long double energy(const Field& w)const
    {long double result=0;for(std::size_t i=0;i<w.size();++i)result+=.5L*op.material().midpointMass().masses()[i]*w[i]*w[i];return result;}
    Field rhs(const Field& w,bool donor=false)
    {
        const auto& m=op.material();const auto& segments=m.geometry().segments();
        Field result(w.size()),trace(segments.size());
        if(donor){
            std::vector<long double> rates(w.size());
            for(std::size_t g=0;g<segments.size();++g){const auto& s=segments[g];const auto f=m.faceRates()[g];long double t=0;
                if(internal(s)){t=f>=0?w[s.negativeVelocity]:w[s.positiveVelocity];rates[s.negativeVelocity]+=f*t;rates[s.positiveVelocity]-=f*t;}
                else{const int i=owner(s);const long double outward=(s.negativeVelocity>=0?1:-1)*f;
                    t=outward>=0?w[i]:(paired[g]>=0?w[owner(segments[paired[g]])]:w[i]);rates[i]+=outward*t;}}
            for(std::size_t i=0;i<w.size();++i)result[i]=static_cast<double>((-rates[i]+w[i]*m.cells()[i].defectRate)/m.midpointMass().masses()[i]);
        }else{
            for(std::size_t g=0;g<trace.size();++g)if(!internal(segments[g])){const int i=owner(segments[g]);
                trace[g]=paired[g]>=0?.5*w[i]+.5*w[owner(segments[paired[g]])]:w[i];}
            const auto boundary=CchConvectionBoundary::bind(*mesh,m,trace,m.sources(),m.timeLayer());
            const auto applied=op.apply(*mesh,CchNativeValues(op.layout(),w,Field(op.layout().pressureCount())),boundary,m.sources(),m.timeLayer());
            const long double scale=energy(w)+1;
            for(const auto& c:applied.report().components){maxBoundaryPower=std::max(maxBoundaryPower,relative(c.boundaryKineticPower,scale));
                maxMaterialPower=std::max(maxMaterialPower,relative(c.materialKineticPower,scale));}
            for(std::size_t i=0;i<w.size();++i)result[i]=static_cast<double>(applied.momentumLoads()[i]/m.midpointMass().masses()[i]);
        }
        for(std::size_t i=0;i<w.size();++i)if(fixedAxis[i]){if(!donor)maxAxisLoad=std::max(maxAxisLoad,std::abs(static_cast<long double>(result[i])));result[i]=0;}
        return result;
    }
    Field integrate(int steps,bool donor,const std::string& label)
    {
        Field w=exact(0);const long double e0=energy(w);const double dt=duration/steps;
        for(int k=0;k<=steps;++k){
            const double time=static_cast<double>(k)/steps;
            if(k%std::max(1,steps/32)==0 || k==steps){const auto e=error(w,exact(time));
                history<<label<<','<<rz<<','<<amr<<','<<n<<','<<steps<<','<<k<<','<<time<<','<<e.first<<','<<e.second<<','<<(energy(w)-e0)/e0<<'\n';}
            if(k==steps)break;
            const auto a=rhs(w,donor);Field stage(w.size());
            for(std::size_t i=0;i<w.size();++i)stage[i]=w[i]+.5*dt*a[i];const auto b=rhs(stage,donor);
            for(std::size_t i=0;i<w.size();++i)stage[i]=w[i]+.5*dt*b[i];const auto c=rhs(stage,donor);
            for(std::size_t i=0;i<w.size();++i)stage[i]=w[i]+dt*c[i];const auto d=rhs(stage,donor);
            for(std::size_t i=0;i<w.size();++i){w[i]+=dt*(a[i]+2*b[i]+2*c[i]+d[i])/6;check(std::isfinite(w[i]),"nonfinite RK4 trajectory");}
        }
        ++windows;return w;
    }
};
// H11: optional production stabilization; original H10b modes stay unchanged.
// This remains prescribed-advector vector transport, not an Euler/NS solver.
void controlledStudy(bool probe)
{
    std::filesystem::create_directories("output/series");
    std::ofstream summary("output/series/controlled_summary.csv"),trajectory("output/series/controlled_trajectory.csv"),
        edges("output/series/controlled_edges.csv"),nodes("output/series/controlled_nodes.csv"),checks("output/series/controlled_gates.csv"),
        layouts("output/series/controlled_layout.csv"),overlaps("output/series/controlled_overlaps.csv"),materials("output/series/controlled_material.csv"),
        stepOutput("output/series/controlled_steps.csv");
    check(summary.good()&&trajectory.good()&&edges.good()&&nodes.good()&&checks.good()&&layouts.good()&&overlaps.good()&&materials.good()&&stepOutput.good(),"cannot open controlled convection study");
    summary<<std::setprecision(17)<<"rz,amr,n,profile,method,steps,energy_initial,energy_final,relative_l2,linf,initial_D_over_E,integrated_D,integrated_material_power,integrated_boundary_power,energy_balance_error,initial_active,initial_connections,max_active,max_sensor,minimum_stage_D,maximum_identity,maximum_axis_reaction,layout_elapsed_seconds\n";
    trajectory<<std::setprecision(17)<<"rz,amr,n,profile,method,step,time,energy,relative_l2,linf,D,material_power,boundary_power,active,connections,max_sensor\n";
    edges<<std::setprecision(17)<<"rz,amr,n,profile,method,sample,first,second,segment,partner,F,u_first,u_second,sensor,chi,L_first,D\n";
    nodes<<std::setprecision(17)<<"rz,amr,n,profile,method,sample,dof,component,representative,fixed_axis,mass,u,L,load,R,mass_change_rate\n";
    checks<<std::setprecision(17)<<"gate,status,value,limit\n";
    layouts<<std::setprecision(17)<<"rz,amr,n,dof,component,x,y,h,representative,mass,negative_pressure,positive_pressure\n";
    overlaps<<"rz,amr,n,dof,pressure\n";
    materials<<std::setprecision(17)<<"rz,amr,n,segment,negative,positive,partner,axis,F,I,translation_x,translation_y\n";
    stepOutput<<std::setprecision(17)<<"rz,amr,n,profile,method,step,dt,energy_begin,energy_end,D_a,D_b,D_c,D_d,M_a,M_b,M_c,M_d,B_a,B_b,B_c,B_d\n";
    const auto accept=[&](const std::string& name,bool ok,long double value,long double limit){checks<<name<<','<<(ok?"PASS":"FAIL")<<','<<value<<','<<limit<<'\n';checks.flush();
        std::cout<<"CONTROLLED_GATE,"<<name<<','<<(ok?"PASS":"FAIL")<<','<<value<<','<<limit<<std::endl;if(!ok)++failed;};
    struct Row {long double l2=0,linf=0,decay=0;};
    for(bool rz:{false,true})for(bool amr:{false,true}){
        std::vector<Row> rows;
        for(int n:probe?std::vector<int>{8}:std::vector<int>{8,16,32}){
            const auto started=std::chrono::steady_clock::now();Study s(n,rz,amr);
            const auto controlled=CchKepConvection::buildControlled(*s.mesh,s.op.material(),s.op.periodicPairs(),{});
            const auto& m=s.op.material();const auto& l=s.op.layout();
            std::vector<Vector2D> translation(m.faceRates().size());
            for(const auto& p:s.op.periodicPairs()){translation[p.firstSegment]=p.translation;translation[p.secondSegment]=p.translation*(-1);}
            for(int i=0;i<l.velocityCount();++i){const auto f=l.domain().face(*s.mesh,l.nativeFace(i));const auto& b=m.geometry().volumes()[i];
                layouts<<rz<<','<<amr<<','<<n<<','<<i<<','<<component(l,i)<<','<<f.center.x<<','<<f.center.y<<','<<std::max(b.upper.x-b.lower.x,b.upper.y-b.lower.y)<<','<<s.op.representatives()[i]<<','<<m.midpointMass().masses()[i]<<','<<(f.negativeCell<0?-1:l.pressureDof(f.negativeCell))<<','<<(f.positiveCell<0?-1:l.pressureDof(f.positiveCell))<<'\n';
                for(const auto& p:b.overlaps)overlaps<<rz<<','<<amr<<','<<n<<','<<i<<','<<p.pressure<<'\n';}
            for(std::size_t k=0;k<m.faceRates().size();++k){const auto& g=m.geometry().segments()[k];materials<<rz<<','<<amr<<','<<n<<','<<k<<','<<g.negativeVelocity<<','<<g.positiveVelocity<<','<<s.paired[k]<<','<<g.normalAxis<<','<<m.faceRates()[k]<<','<<m.faceIntegrals()[k]<<','<<translation[k].x<<','<<translation[k].y<<'\n';}
            layouts.flush();overlaps.flush();materials.flush();
            const auto packet=[&](double time){Field w(l.velocityCount());
                for(int i=0;i<l.velocityCount();++i){const auto p=s.points[i];double t=(rz?p.y:p.x)-time;t-=std::floor(t);const double z=(t-.5)/.2;
                    const double value=std::abs(z)<1?std::pow(1-z*z,3)*std::sin(static_cast<double>(2*pi)*.35*n*(t-.5)):0;
                    w[i]=component(l,i)==0?(rz?p.x:1)*value:.4+.7*value;}
                return w;};
            for(int profile=0;profile<2;++profile)for(int method=profile==0?1:0;method<2;++method){
                const auto& op=method?controlled:s.op;const std::string name=std::to_string(rz)+"_"+std::to_string(amr)+"_"+std::to_string(n)+"_"+std::to_string(profile)+"_"+std::to_string(method);
                const int steps=(profile==0?32:128)*n;const double dt=duration/steps;Field w=profile==0?s.exact(0):packet(0);
                for(int i=0;i<l.velocityCount();++i)w[i]=w[op.representatives()[i]];
                const auto initialEnergy=s.energy(w);long double sumD=0,sumMaterial=0,sumBoundary=0,minimumD=std::numeric_limits<long double>::infinity(),maxSensor=0,maxIdentity=0,maxAxisReaction=0;
                int initialActive=0,initialConnections=0,maxActive=0;long double initialDecay=0;
                struct Evaluation {Field rhs;long double D=0,material=0,boundary=0,sensor=0;int active=0,connections=0;};
                const auto evaluate=[&](const Field& input,const std::string& sample){
                    const auto result=op.applyPaired(*s.mesh,CchNativeValues(l,input,Field(l.pressureCount())),m.sources(),m.timeLayer());
                    Evaluation e;e.rhs.resize(input.size());
                    for(const auto& c:result.report().components){e.D+=c.dissipationPower;e.material+=c.materialKineticPower;e.boundary+=c.boundaryKineticPower;
                        e.sensor=std::max(e.sensor,c.maximumSensor);e.active+=c.activeDissipationConnections;e.connections+=c.dissipationConnections;maxIdentity=std::max(maxIdentity,c.kineticIdentity);}
                    check(e.D>=0,"negative controlled stage dissipation");minimumD=std::min(minimumD,e.D);maxActive=std::max(maxActive,e.active);maxSensor=std::max(maxSensor,e.sensor);
                    for(int i=0;i<l.velocityCount();++i){e.rhs[i]=static_cast<double>(result.momentumLoads()[i]/m.midpointMass().masses()[i]);
                        if(s.fixedAxis[i]){check(input[i]==0,"nonzero prescribed radial axis velocity");maxAxisReaction=std::max(maxAxisReaction,std::abs(result.momentumLoads()[i]));e.rhs[i]=0;}}
                    if(!sample.empty()){
                        for(const auto& edge:result.dissipationSamples())edges<<rz<<','<<amr<<','<<n<<','<<profile<<','<<method<<','<<sample<<','<<edge.first<<','<<edge.second<<','<<edge.firstSegment<<','<<edge.secondSegment<<','<<m.faceRates()[edge.firstSegment]<<','<<input[edge.first]<<','<<input[edge.second]<<','<<edge.sensor<<','<<edge.chi<<','<<edge.momentumRate<<','<<edge.power<<'\n';
                        for(int i=0;i<l.velocityCount();++i)nodes<<rz<<','<<amr<<','<<n<<','<<profile<<','<<method<<','<<sample<<','<<i<<','<<component(l,i)<<','<<op.representatives()[i]<<','<<s.fixedAxis[i]<<','<<m.midpointMass().masses()[i]<<','<<input[i]<<','<<result.dissipationRates()[i]<<','<<result.momentumLoads()[i]<<','<<m.cells()[i].defectRate<<','<<m.cells()[i].massChange/m.deltaTime()<<'\n';
                    }
                    return e;};
                for(int step=0;step<=steps;++step){const double time=static_cast<double>(step)/steps;
                    const auto a=evaluate(w,step==0?"initial":step==steps?"terminal":"");
                    if(step==0){initialActive=a.active;initialConnections=a.connections;initialDecay=a.D/initialEnergy;}
                    if(step%std::max(1,steps/32)==0 || step==steps){const auto error=s.error(w,profile==0?s.exact(time):packet(time));
                        trajectory<<rz<<','<<amr<<','<<n<<','<<profile<<','<<method<<','<<step<<','<<time<<','<<s.energy(w)<<','<<error.first<<','<<error.second<<','<<a.D<<','<<a.material<<','<<a.boundary<<','<<a.active<<','<<a.connections<<','<<a.sensor<<'\n';trajectory.flush();}
                    if(step==steps)break;
                    const auto energyBefore=s.energy(w);
                    Field stage(w.size());for(std::size_t i=0;i<w.size();++i)stage[i]=w[i]+.5*dt*a.rhs[i];const auto b=evaluate(stage,"");
                    for(std::size_t i=0;i<w.size();++i)stage[i]=w[i]+.5*dt*b.rhs[i];const auto c=evaluate(stage,"");
                    for(std::size_t i=0;i<w.size();++i)stage[i]=w[i]+dt*c.rhs[i];const auto d=evaluate(stage,"");
                    for(std::size_t i=0;i<w.size();++i){w[i]+=dt*(a.rhs[i]+2*b.rhs[i]+2*c.rhs[i]+d.rhs[i])/6;check(std::isfinite(w[i]),"nonfinite controlled RK4 trajectory");}
                    sumD+=dt*(a.D+2*b.D+2*c.D+d.D)/6;sumMaterial+=dt*(a.material+2*b.material+2*c.material+d.material)/6;sumBoundary+=dt*(a.boundary+2*b.boundary+2*c.boundary+d.boundary)/6;
                    stepOutput<<rz<<','<<amr<<','<<n<<','<<profile<<','<<method<<','<<step<<','<<dt<<','<<energyBefore<<','<<s.energy(w)<<','<<a.D<<','<<b.D<<','<<c.D<<','<<d.D<<','<<a.material<<','<<b.material<<','<<c.material<<','<<d.material<<','<<a.boundary<<','<<b.boundary<<','<<c.boundary<<','<<d.boundary<<'\n';
                }
                const auto energy=s.energy(w);const auto err=s.error(w,profile==0?s.exact(duration):packet(duration));
                const auto balance=std::abs(energy-initialEnergy+sumD-sumMaterial+sumBoundary)/initialEnergy;
                const auto seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
                summary<<rz<<','<<amr<<','<<n<<','<<profile<<','<<method<<','<<steps<<','<<initialEnergy<<','<<energy<<','<<err.first<<','<<err.second<<','<<initialDecay<<','<<sumD<<','<<sumMaterial<<','<<sumBoundary<<','<<balance<<','<<initialActive<<','<<initialConnections<<','<<maxActive<<','<<maxSensor<<','<<minimumD<<','<<maxIdentity<<','<<maxAxisReaction<<','<<seconds<<'\n';summary.flush();edges.flush();nodes.flush();stepOutput.flush();++windows;
                accept(name+"_identity",maxIdentity<=1e-12L,maxIdentity,1e-12L);
                if(profile==0)rows.push_back({err.first,err.second,initialDecay});
                else{
                    accept(name+"_energy_budget",balance<=1e-6L,balance,1e-6L);
                    accept(name+"_energy_growth",energy<=initialEnergy*(1+1e-6L),energy/initialEnergy,1+1e-6L);
                    accept(name+"_fixed_material",std::abs(sumMaterial)+std::abs(sumBoundary)<=1e-12L*initialEnergy,(std::abs(sumMaterial)+std::abs(sumBoundary))/initialEnergy,1e-12L);
                    if(method){accept(name+"_active",sumD>0&&initialActive>0,sumD,0);
                        if(n==32)accept(name+"_local",initialActive<initialConnections,initialActive,initialConnections);}}
                std::cout<<"CONTROLLED_WINDOW,"<<rz<<','<<amr<<','<<n<<','<<profile<<','<<method<<','<<err.first<<','<<balance<<','<<seconds<<std::endl;
            }
        }
        if(!probe){const auto& fine=rows[2];const auto& medium=rows[1];const auto prefix=std::to_string(rz)+"_"+std::to_string(amr);
            const auto p=std::log(medium.l2/fine.l2)/std::log(2.L);accept(prefix+"_smooth_order",p>=1.7L,p,1.7L);
            accept(prefix+"_smooth_l2",fine.l2<=.05L,fine.l2,.05L);accept(prefix+"_smooth_linf",fine.linf<medium.linf,fine.linf,medium.linf);
            const auto decayOrder=fine.decay<=1e-14L?0.L:std::log(medium.decay/fine.decay)/std::log(2.L);
            accept(prefix+"_smooth_D",fine.decay<=1e-14L || (medium.decay>fine.decay && decayOrder>=1.7L),decayOrder,1.7L);}
    }
    std::cout<<"CONTROLLED_COMPLETE,"<<windows<<','<<failed<<std::endl;
}
// H10n observer begin
long double donorDissipation(Study& s,const Field& w,const std::string& sample)
{
    const auto& m=s.op.material();const auto& segments=m.geometry().segments();
    const auto rhs=s.rhs(w,true);std::vector<long double> q(w.size());
    long double inside=0,periodic=0,openPower=0,correction=0,power=0,gross=0;
    const auto ratio=[](long double x,long double scale){return scale>0?std::abs(x)/scale:(x==0?0:std::numeric_limits<long double>::infinity());};
    check(std::max({s.maxDefect,s.maxMassChange,s.maxGeometryChange,s.maxPairError})<=1e-12L,"donor dissipation control is not fixed material with matching ports");
    for(std::size_t g=0;g<segments.size();++g){const auto& segment=segments[g];const auto f=m.faceRates()[g];
        if(segment.negativeVelocity>=0)q[segment.negativeVelocity]+=f;
        if(segment.positiveVelocity>=0)q[segment.positiveVelocity]-=f;
        int i=-1,j=-1,kind=0;long double rate=f,other=-f,term=0;
        if(internal(segment)){i=segment.negativeVelocity;j=segment.positiveVelocity;}
        else{
            i=owner(segment);rate=(segment.negativeVelocity>=0?1:-1)*f;
            if(s.paired[g]>=0){
                if(static_cast<int>(g)>s.paired[g])continue;
                const auto& p=segments[s.paired[g]];j=owner(p);kind=1;other=(p.negativeVelocity>=0?1:-1)*m.faceRates()[s.paired[g]];
                check(ratio(rate+other,std::abs(rate)+std::abs(other))<=1e-12L,"donor port antisymmetry failed");
            }else{kind=2;other=0;openPower+=.5L*rate*w[i]*w[i];gross+=std::abs(.5L*rate*w[i]*w[i]);}
        }
        if(j>=0){const long double jump=static_cast<long double>(w[i])-w[j];term=.5L*std::abs(rate)*jump*jump;if(kind==0)inside+=term;else periodic+=term;}
        donorSegments<<s.rz<<','<<s.amr<<','<<s.n<<','<<sample<<','<<g<<','<<s.paired[g]<<','<<kind<<','<<i<<','<<j<<','<<rate<<','<<other<<','<<term<<'\n';
    }
    for(std::size_t i=0;i<w.size();++i){const long double mass=m.midpointMass().masses()[i],R=m.cells()[i].defectRate;
        const long double u=w[i],work=mass*u*rhs[i];power+=work;correction+=(R-.5L*q[i])*u*u;
        gross+=std::abs(work)+(std::abs(R)+.5L*std::abs(q[i]))*u*u;
        if(s.fixedAxis[i])check(w[i]==0 && rhs[i]==0,"nonzero donor fixed-axis work");
        donorNodes<<s.rz<<','<<s.amr<<','<<s.n<<','<<sample<<','<<i<<','<<s.fixedAxis[i]<<','<<mass<<','<<u<<','<<rhs[i]<<','<<q[i]<<','<<R<<','<<m.cells()[i].massChange/m.deltaTime()<<'\n';
    }
    const auto D=inside+periodic,E=s.energy(w);gross+=D;
    const auto identity=ratio(power+D+openPower-correction,gross),fixed=ratio(power+D,gross);
    check(E>0 && std::isfinite(D) && D>=0 && identity<=1e-12L && fixed<=1e-12L && ratio(correction,gross)<=1e-12L && ratio(openPower,gross)<=1e-12L,"independent donor energy identity failed");
    donorPower<<s.rz<<','<<s.amr<<','<<s.n<<','<<sample<<','<<E<<','<<inside<<','<<periodic<<','<<openPower<<','<<correction<<','<<power<<','<<gross<<','<<identity<<','<<fixed<<','<<D/E<<'\n';
    donorPower.flush();donorSegments.flush();donorNodes.flush();return D/E;
}
// H10n observer end
struct Row{long double truncation,l2,linf,drift,donorLoss,donorDecay;};
Row spatial(Study& s)
{
    const auto initial=s.exact(0);const auto trunc=s.error(s.rhs(initial),s.exact(0,true));
    const auto decay=donorDissipation(s,initial,"initial");
    const int steps=32*s.n;const auto center=s.integrate(steps,false,"center"),donor=s.integrate(steps,true,"donor");
    const auto e=s.error(center,s.exact(duration));const long double e0=s.energy(initial);
    donorDissipation(s,donor,"terminal");
    const Row r{trunc.first,e.first,e.second,(s.energy(center)-e0)/e0,(e0-s.energy(donor))/e0,decay};
    summaries<<s.rz<<','<<s.amr<<','<<s.n<<','<<s.op.layout().velocityCount()<<','<<s.coarseFine<<','<<steps<<','<<r.truncation<<','<<r.l2<<','<<r.linf<<','<<r.drift<<','<<r.donorLoss<<','<<s.maxBoundaryPower<<','<<s.maxMaterialPower<<','<<s.maxDefect<<','<<s.maxMassChange<<','<<s.maxGeometryChange<<','<<s.maxAxisLoad<<'\n';summaries.flush();history.flush();
    std::cout<<"SPACE,"<<s.rz<<','<<s.amr<<','<<s.n<<','<<r.truncation<<','<<r.l2<<','<<r.linf<<','<<r.drift<<','<<r.donorLoss<<std::endl;
    const std::string name=std::string(s.rz?"rz":"xy")+(s.amr?"_amr":"_uniform")+"_n"+std::to_string(s.n);
    gate(name+"_paired_flux",s.maxPairError<=1e-12,s.maxPairError,1e-12);
    gate(name+"_boundary_power",s.maxBoundaryPower<=1e-12,s.maxBoundaryPower,1e-12);
    gate(name+"_material_power",s.maxMaterialPower<=1e-12,s.maxMaterialPower,1e-12);
    const auto defect=std::max({s.maxDefect,s.maxMassChange,s.maxGeometryChange});gate(name+"_fixed_material",defect<=1e-12,defect,1e-12);
    if(s.amr)gate(name+"_actual_coarse_fine",s.coarseFine>0,s.coarseFine,0);
    return r;
}
long double order(long double coarse,long double fine)
{return coarse>0 && fine>0?std::log(coarse/fine)/std::log(2.L):std::numeric_limits<long double>::quiet_NaN();}
void full()
{
    for(bool rz:{false,true})for(bool amr:{false,true}){
        std::vector<Row> rows;for(int n:{8,16,32}){Study s(n,rz,amr);rows.push_back(spatial(s));}
        const std::string name=std::string(rz?"rz":"xy")+(amr?"_amr":"_uniform");
        for(int j=1;j<3;++j)std::cout<<"SPACE_ORDER,"<<rz<<','<<amr<<','<<(j==1?16:32)<<','<<order(rows[j-1].truncation,rows[j].truncation)<<','<<order(rows[j-1].l2,rows[j].l2)<<','<<order(rows[j-1].donorLoss,rows[j].donorLoss)<<std::endl;
        const auto& r=rows.back();const auto truncOrder=order(rows[1].truncation,r.truncation),solutionOrder=order(rows[1].l2,r.l2),lossOrder=order(rows[1].donorLoss,r.donorLoss);
        gate(name+"_truncation_order",truncOrder>=1.7,truncOrder,1.7);
        gate(name+"_solution_order",solutionOrder>=1.7,solutionOrder,1.7);
        gate(name+"_fine_l2",r.l2<=.05,r.l2,.05);
        gate(name+"_linf_decreases",rows[0].linf>rows[1].linf && rows[1].linf>r.linf,r.linf,rows[1].linf);
        gate(name+"_energy_drift",std::abs(r.drift)<=1e-6,std::abs(r.drift),1e-6);
        gate(name+"_energy_vs_donor",r.donorLoss>0 && std::abs(r.drift)<=.01*r.donorLoss,std::abs(r.drift),.01*r.donorLoss);
        legacyGate(name+"_donor_loss_order",lossOrder>=.7 && lossOrder<=1.4,lossOrder,1.4);
        const auto decayOrder=order(rows[1].donorDecay,r.donorDecay);
        gate(name+"_donor_decay_order",decayOrder>=.7 && decayOrder<=1.4,decayOrder,1.4);
    }
    Study timeStudy(8,false,false);const auto reference=timeStudy.integrate(640,false,"time_reference");
    long double previous=0;for(int steps:{20,40,80}){
        const auto w=timeStudy.integrate(steps,false,"time");const auto e=timeStudy.error(w,reference).first;
        std::cout<<"TIME,"<<steps<<','<<e<<','<<timeStudy.error(w,timeStudy.exact(duration)).first<<std::endl;
        if(previous>0){const auto p=order(previous,e);gate("rk4_time_order_"+std::to_string(steps),p>=3.5,p,3.5);}previous=e;
    }
    std::cout<<"SUMMARY,"<<windows<<','<<failed<<','<<(failed?"FAIL":"PASS")<<std::endl;
    std::cout<<"LEGACY_SUMMARY,4,"<<legacyFailed<<','<<(legacyFailed?"FAIL":"PASS")<<std::endl;
}
// H10m observer begin: reuse the original donor, never the center operator.
void donorAudit()
{
    std::ofstream out("output/series/donor_audit.csv"),modes("output/series/donor_fourier.csv");
    check(out.good() && modes.good(),"cannot open donor gate audit");
    out<<std::setprecision(17)<<"rz,amr,n,steps,loss,loss_order,initial_decay_rate,decay_order,fourier_loss,fourier_actual_loss,fourier_max_error\n";
    modes<<std::setprecision(17)<<"rz,n,dof,actual,independent_fourier,error\n";
    for(bool rz:{false,true})for(bool amr:{false,true}){
        long double previousLoss=0,previousDecay=0,previousDissipation=0;
        for(int n:{8,16,32}){
            Study s(n,rz,amr);const auto w0=s.exact(0),rhs=s.rhs(w0,true);
            const auto independentDecay=donorDissipation(s,w0,"initial");
            const auto& mass=s.op.material().midpointMass().masses();
            const long double e0=s.energy(w0);long double power=0;
            for(std::size_t i=0;i<w0.size();++i)power+=mass[i]*w0[i]*rhs[i];
            const long double decay=-power/e0;const int steps=32*n;
            const auto w=s.integrate(steps,true,"donor_audit");
            const long double loss=(e0-s.energy(w))/e0;
            donorDissipation(s,w,"terminal");check(loss>0,"donor has no detectable energy loss");
            const auto p=previousLoss>0?order(previousLoss,loss):0;
            const auto pd=previousDecay>0?order(previousDecay,decay):0;
            long double fourierLoss=0,actualLoss=0,maxError=0;
            if(!amr){
                // Regular periodic transverse samples; true half-dual masses retained.
                using Complex=std::complex<long double>;
                const Complex imaginary(0,1),lambda=-static_cast<long double>(n)*(1.L-std::exp(-imaginary*(2*pi/n)));
                const auto z=lambda*static_cast<long double>(duration/steps),z2=z*z;
                const auto rk=1.L+z+z2/2.L+z2*z/6.L+z2*z2/24.L;
                const auto amplification=std::pow(rk,steps);
                fourierLoss=1-std::norm(amplification);
                long double start=0,end=0;
                for(std::size_t i=0;i<w.size();++i)if(component(s.op.layout(),static_cast<int>(i))==(rz?0:1)){
                    const auto point=s.points[i];const auto phase=std::exp(imaginary*(2*pi*(rz?point.y:point.x)));
                    const auto expected=static_cast<long double>(rz?point.x:1)*std::imag(amplification*phase);
                    const auto error=std::abs(w[i]-expected);maxError=std::max(maxError,error);
                    start+=mass[i]*w0[i]*w0[i];end+=mass[i]*w[i]*w[i];
                    modes<<rz<<','<<n<<','<<i<<','<<w[i]<<','<<expected<<','<<error<<'\n';
                }
                actualLoss=(start-end)/start;
                check(maxError<=1e-12L && std::abs(actualLoss-fourierLoss)<=1e-12L,"original donor differs from independent RK4 Fourier mode");
            }
            out<<rz<<','<<amr<<','<<n<<','<<steps<<','<<loss<<','<<p<<','<<decay<<','<<pd<<','<<fourierLoss<<','<<actualLoss<<','<<maxError<<'\n';out.flush();modes.flush();history.flush();
            std::cout<<"DONOR_AUDIT,"<<rz<<','<<amr<<','<<n<<','<<loss<<','<<p<<','<<decay<<','<<pd<<','<<maxError<<std::endl;
            if(n==32){const std::string name=std::string(rz?"rz":"xy")+(amr?"_amr":"_uniform");
                legacyGate(name+"_donor_loss_order",p>=.7 && p<=1.4,p,1.4);
                const auto pD=order(previousDissipation,independentDecay);
                gate(name+"_donor_decay_order",pD>=.7 && pD<=1.4,pD,1.4);}
            previousLoss=loss;previousDecay=decay;previousDissipation=independentDecay;
        }
    }
    std::cout<<"DONOR_SUMMARY,"<<windows<<','<<failed<<','<<(failed?"FAIL":"PASS")<<std::endl;
    std::cout<<"LEGACY_SUMMARY,4,"<<legacyFailed<<','<<(legacyFailed?"FAIL":"PASS")<<std::endl;
}
// H10m observer end
void rowAudit(const std::vector<int>& grids)
{
    // Independent sparse native matrix, including the actual research port.
    // Observer only: no coefficients or reference values enter the trajectory.
    std::ofstream out("output/series/rows.csv");check(out.good(),"cannot open native row audit");
    out<<std::setprecision(17)<<"rz,amr,n,dof,component,x,y,exterior,fixed_axis,mass,h,actual,error,squared_contribution,row_scale,moment1,moment2_over_h,moment3_over_h2,reference_error\n";
    for(bool rz:{false,true})for(bool amr:{false,true})for(int n:grids){
        Study s(n,rz,amr);const auto& m=s.op.material();const auto& segments=m.geometry().segments();const int count=s.op.layout().velocityCount();
        std::vector<std::map<int,long double>> matrix(count);std::vector<bool> exterior(count);
        for(std::size_t g=0;g<segments.size();++g){const auto& p=segments[g];const auto f=m.faceRates()[g];
            if(internal(p)){const int i=p.negativeVelocity,j=p.positiveVelocity;
                for(const auto& w:s.op.weights()[g]){matrix[i][w.velocity]+=.5L*f*w.weight;matrix[j][w.velocity]-=.5L*f*w.weight;matrix[w.velocity][i]-=.5L*f*w.weight;matrix[w.velocity][j]+=.5L*f*w.weight;}
                matrix[i][i]+=.5L*f;matrix[j][j]-=.5L*f;
            }else{const int i=owner(p);exterior[i]=true;const long double outward=(p.negativeVelocity>=0?1:-1)*f;
                if(s.paired[g]>=0){matrix[i][i]+=.5L*outward;matrix[i][owner(segments[s.paired[g]])]+=.5L*outward;}else matrix[i][i]+=outward;}
        }
        for(const auto& p:s.op.couplings()){matrix[p.first][p.second]+=p.coefficient;matrix[p.second][p.first]-=p.coefficient;}
        std::vector<std::vector<int>> classes(count);for(int i=0;i<count;++i)classes[s.op.representatives()[i]].push_back(i);
        for(const auto& members:classes)if(members.size()>1){std::map<int,long double> sum;long double mass=0;
            for(int i:members){mass+=m.midpointMass().masses()[i];for(const auto& e:matrix[i])sum[e.first]+=e.second;sum[i]-=m.cells()[i].outwardRate;}
            for(int i:members){matrix[i].clear();for(const auto& e:sum)matrix[i][e.first]=m.midpointMass().masses()[i]/mass*e.second;matrix[i][i]+=m.cells()[i].outwardRate;}}
        const auto& moment=s.op.momentReport();std::cout<<"BOUNDED_MOMENTS,"<<rz<<','<<amr<<','<<n<<','<<moment.iterations<<','<<moment.projectionCalls<<','<<moment.constantError<<','<<moment.affineError<<','<<moment.quadraticBoundError<<','<<moment.coefficientRatio<<','<<moment.maximumSupportRatio<<std::endl;
        const auto w=s.exact(0),expected=s.exact(0,true),actual=s.rhs(w);long double denominator=0;
        for(int i=0;i<count;++i)denominator+=m.midpointMass().masses()[i]*expected[i]*expected[i];
        std::array<long double,3> squared{},maximumError{},maximumScale{},maximumThird{};long double referenceError=0;
        for(int i=0;i<count;++i){const auto mass=m.midpointMass().masses()[i];const auto& box=m.geometry().volumes()[i];const double h=std::max(box.upper.x-box.lower.x,box.upper.y-box.lower.y);
            long double independent=w[i]*m.cells()[i].defectRate/mass,gross=std::abs(independent),scale=0;std::array<long double,3> moments{};
            for(const auto& entry:matrix[i]){const auto k=-entry.second/mass;independent+=k*w[entry.first];gross+=std::abs(k*w[entry.first]);scale+=std::abs(k)*h;
                double delta=rz?s.points[entry.first].y-s.points[i].y:s.points[entry.first].x-s.points[i].x;delta-=std::round(delta); // declared period is one
                moments[0]+=k*delta;moments[1]+=k*delta*delta/h;moments[2]+=k*delta*delta*delta/(h*h);}
            if(s.fixedAxis[i])independent=0;
            const auto discrepancy=relative(independent-actual[i],gross);referenceError=std::max(referenceError,discrepancy);check(discrepancy<=1e-12L,"independent native row matrix differs from actual apply");
            const long double e=static_cast<long double>(actual[i])-expected[i],contribution=mass*e*e/denominator;const int group=s.fixedAxis[i]?2:(exterior[i]?1:0);
            squared[group]+=contribution;maximumError[group]=std::max(maximumError[group],std::abs(e));maximumScale[group]=std::max(maximumScale[group],scale);maximumThird[group]=std::max(maximumThird[group],std::abs(moments[2]));
            out<<rz<<','<<amr<<','<<n<<','<<i<<','<<component(s.op.layout(),i)<<','<<s.points[i].x<<','<<s.points[i].y<<','<<exterior[i]<<','<<s.fixedAxis[i]<<','<<mass<<','<<h<<','<<actual[i]<<','<<e<<','<<contribution<<','<<scale<<','<<moments[0]<<','<<moments[1]<<','<<moments[2]<<','<<discrepancy<<'\n';
        }
        const auto full=s.error(actual,expected).first;check(relative(squared[0]+squared[1]+squared[2]-full*full,full*full)<=1e-12L,"native row contributions do not sum to full error");
        for(int group=0;group<3;++group)std::cout<<"ROW_AUDIT,"<<rz<<','<<amr<<','<<n<<','<<group<<','<<std::sqrt(squared[group])<<','<<maximumError[group]<<','<<maximumScale[group]<<','<<maximumThird[group]<<','<<referenceError<<std::endl;
    }
}
void affineAudit(const std::vector<int>& grids)
{
    std::ofstream out("output/series/affine.csv");check(out.good(),"cannot open affine audit");
    out<<std::setprecision(17)<<"rz,amr,n,dof,component,x,y,exterior,coarse_fine,actual_derivative,exact_derivative,absolute_error\n";
    for(bool rz:{false,true})for(bool amr:{false,true})for(int n:grids){
        Study s(n,rz,amr,false);const auto& m=s.op.material();const auto& geometry=m.geometry();
        const auto& moments=s.op.momentReport();
        std::cout<<"NATIVE_MOMENTS,"<<rz<<','<<amr<<','<<n<<','<<s.op.couplings().size()<<','<<moments.constraints<<','<<moments.iterations<<','<<moments.constantError<<','<<moments.affineError<<','<<moments.quadraticError<<std::endl;
        Field w(s.points.size()),t(geometry.segments().size());std::vector<bool> cf(w.size());
        // XY: both transported components x. RZ: radial r, axial z.
        for(std::size_t i=0;i<w.size();++i)w[i]=rz?(component(s.op.layout(),static_cast<int>(i))==0?s.points[i].x:s.points[i].y):s.points[i].x;
        for(std::size_t g=0;g<t.size();++g){const auto& segment=geometry.segments()[g];
            if(internal(segment)){
                int low=100,high=-1;for(int u:{segment.negativeVelocity,segment.positiveVelocity})for(const auto& p:geometry.volumes()[u].overlaps){
                    const int level=s.mesh->cells(s.mesh->cellIndexToLabel[s.op.layout().nativeCell(p.pressure)]).level();low=std::min(low,level);high=std::max(high,level);}
                if(low<high){cf[segment.negativeVelocity]=cf[segment.positiveVelocity]=true;for(const auto& weight:s.op.weights()[g])cf[weight.velocity]=true;}
            }else{
                auto p=segment.begin+(segment.end-segment.begin)*.5;
                if(rz && segment.normalAxis==1){const long double a=segment.begin.x,b=segment.end.x;p.x=static_cast<double>((2.L/3)*(b*b*b-a*a*a)/(b*b-a*a));}
                t[g]=rz?(component(s.op.layout(),owner(segment))==0?p.x:p.y):p.x;
            }
        }
        const auto boundary=CchConvectionBoundary::bind(*s.mesh,m,t,m.sources(),m.timeLayer());
        const auto value=s.op.apply(*s.mesh,CchNativeValues(s.op.layout(),w,Field(s.op.layout().pressureCount())),boundary,m.sources(),m.timeLayer());
        long double interior=0,coarseFine=0,exterior=0;
        for(std::size_t i=0;i<w.size();++i){const auto& box=geometry.volumes()[i];const auto p=s.points[i];
            const bool ext=box.lower.x==0 || box.lower.y==0 || box.upper.x==1 || box.upper.y==1;
            const int c=component(s.op.layout(),static_cast<int>(i));const long double expected=rz && c==0?0:-1;
            const long double actual=value.momentumLoads()[i]/m.midpointMass().masses()[i],e=std::abs(actual-expected);
            if(ext)exterior=std::max(exterior,e);else{interior=std::max(interior,e);if(cf[i])coarseFine=std::max(coarseFine,e);}
            out<<rz<<','<<amr<<','<<n<<','<<i<<','<<c<<','<<p.x<<','<<p.y<<','<<ext<<','<<cf[i]<<','<<actual<<','<<expected<<','<<e<<'\n';
        }
        std::cout<<"AFFINE_AUDIT,"<<rz<<','<<amr<<','<<n<<','<<interior<<','<<coarseFine<<','<<exterior<<std::endl;
        check(interior<=1e-10L,"final native interior affine counterexample failed");
    }
}
// Necessary-condition certificate only. No corrected target, convection pair,
// material rate or state is returned. The SVD is an independent norm check.
void materialMomentBounds(std::ofstream& out,std::ofstream& sources,const QuadtreeMesh& mesh,bool rz,bool amr,int n,double rho,
    const CchDualMaterial& material,const std::vector<Vector2D>& points,
    const std::vector<long double>& q,const std::vector<Vector2D>& advector)
{
    const auto& l=material.layout();const auto& masses=material.midpointMass().masses();
    long double c0=0,c1=0,gross0=0,gross1=0,totalMass=0,bSum=0,zbSum=0;
    std::vector<long double> xi(l.pressureCount()),eta(l.pressureCount());
    std::vector<long double> z,m,h;
    for(int i=0;i<l.velocityCount();++i){totalMass+=masses[i];if(component(l,i)!=0)continue;
        const long double coordinate=rz?points[i].x:points[i].y,b=masses[i]*(rz?advector[i].x:advector[i].y);
        c0-=b+coordinate*q[i];c1-=coordinate*b+.5L*q[i]*coordinate*coordinate;
        bSum+=b;zbSum+=coordinate*b;
        for(const auto& overlap:material.geometry().volumes()[i].overlaps){xi[overlap.pressure]+=overlap.cellFraction*coordinate;eta[overlap.pressure]+=overlap.cellFraction*coordinate*coordinate;}
        gross0+=std::abs(b)+std::abs(coordinate*q[i]);gross1+=std::abs(coordinate*b)+std::abs(.5L*q[i]*coordinate*coordinate);
        if(rz && coordinate==0)continue; // Fixed axis retained in the above sums.
        const auto& box=material.geometry().volumes()[i];z.push_back(coordinate);m.push_back(masses[i]);
        h.push_back(std::max(box.upper.x-box.lower.x,box.upper.y-box.lower.y));}
    check(totalMass>0 && !z.empty(),"empty material moment certificate");
    long double bound[2]{},variance[2]{},reference=0;
    for(int scaled=0;scaled<2;++scaled){
        long double sumW=0,sumWz=0;std::vector<long double> weights(z.size());
        for(std::size_t i=0;i<z.size();++i){weights[i]=m[i]*(scaled?h[i]*h[i]*h[i]*h[i]:1);sumW+=weights[i];sumWz+=weights[i]*z[i];}
        check(sumW>0 && std::isfinite(sumW),"invalid material certificate weight");const long double mean=sumWz/sumW;
        long double d=0;for(std::size_t i=0;i<z.size();++i)d+=weights[i]*(z[i]-mean)*(z[i]-mean);
        check(d>0 && std::isfinite(d),"unidentifiable material moment Gram");variance[scaled]=d/sumW;
        const long double centered=c1-mean*c0,minimum=c0*c0/sumW+centered*centered/d;
        bound[scaled]=std::sqrt(minimum/totalMass);
        Eigen::MatrixXd a(2,z.size());for(std::size_t i=0;i<z.size();++i){a(0,i)=std::sqrt(weights[i]/sumW);a(1,i)=std::sqrt(weights[i]/d)*(z[i]-mean);}
        const Eigen::Vector2d rhs(static_cast<double>(c0/std::sqrt(sumW)),static_cast<double>(centered/std::sqrt(d)));
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(a,Eigen::ComputeThinU|Eigen::ComputeThinV);svd.setThreshold(1e-12);
        check(svd.rank()==2,"material certificate SVD rank differs");const Eigen::VectorXd witness=svd.solve(rhs);
        reference=std::max(reference,relative(witness.squaredNorm()-minimum,minimum));
        reference=std::max(reference,static_cast<long double>((a*witness-rhs).norm()/std::max(1.,rhs.norm())));
        // Dual Cauchy--Schwarz check on a deterministic unrelated error field.
        Eigen::VectorXd trial(z.size());for(std::size_t i=0;i<z.size();++i)trial[i]=std::sin(3*z[i]+static_cast<double>(i%7));
        const Eigen::Vector2d sampled=a*trial;
        check(sampled.squaredNorm()<=trial.squaredNorm()*(1+1e-12),"material certificate exceeds actual error norm");}
    check(reference<=1e-12L,"independent material certificate SVD differs");
    long double mapped0=-bSum,mapped1=-zbSum,face0=-bSum,face1=-zbSum,firstTransfer=0,secondTransfer=0;
    const auto& primal=material.material();const auto& di=primal.midpointRate().geometry().metrics().divergence().matrix();
    // D_I entries already contain face area; primal faceRates are integrated
    // material rates, so their scatter uses orientation only, not area twice.
    for(int p=0;p<l.pressureCount();++p){long double primalQ=0;for(CchIntegralDivergence::Matrix::InnerIterator it(di,p);it;++it)primalQ+=(it.value()>0?1.L:-1.L)*primal.faceRates()[it.col()];
        const long double declared=primal.cells()[p].outwardIntegral/primal.deltaTime();reference=std::max(reference,relative(primalQ-declared,std::abs(primalQ)+std::abs(declared)));
        const auto point=mesh.cells(mesh.cellIndexToLabel[l.nativeCell(p)]).centerCoord();const long double coordinate=rz?point.x:point.y;
        firstTransfer=std::max(firstTransfer,std::abs(xi[p]-coordinate));secondTransfer=std::max(secondTransfer,std::abs(eta[p]-coordinate*coordinate));
        mapped0-=xi[p]*primalQ;mapped1-=.5L*eta[p]*primalQ;
        sources<<rz<<','<<amr<<','<<n<<','<<rho<<','<<p<<','<<point.x<<','<<point.y<<','<<xi[p]<<','<<eta[p]<<','<<primalQ<<'\n';}
    // Read actual primal owner faces directly, independently of D_I rows.
    for(int i=0;i<l.velocityCount();++i){const auto face=l.domain().face(mesh,l.nativeFace(i));const int a=face.negativeCell<0?-1:l.pressureDof(face.negativeCell),b=face.positiveCell<0?-1:l.pressureDof(face.positiveCell);
        const long double f=primal.faceRates()[i];face0-=f*((a<0?0:xi[a])-(b<0?0:xi[b]));face1-=.5L*f*((a<0?0:eta[a])-(b<0?0:eta[b]));}
    const long double sourceError=std::max({relative(mapped0-c0,gross0),relative(mapped1-c1,gross1),relative(face0-c0,gross0),relative(face1-c1,gross1)});
    check(sourceError<=1e-12L && reference<=1e-12L,"material source moment decomposition differs");
    out<<rz<<','<<amr<<','<<n<<','<<rho<<','<<c0<<','<<c1<<','<<gross0<<','<<gross1<<','<<totalMass<<','<<variance[0]<<','<<variance[1]<<','<<bound[0]<<','<<bound[1]<<','<<reference<<','<<mapped0<<','<<mapped1<<','<<face0<<','<<face1<<','<<firstTransfer<<','<<secondTransfer<<','<<sourceError<<'\n';
    std::cout<<"MATERIAL_MOMENT_BOUND,"<<rz<<','<<amr<<','<<n<<','<<rho<<','<<bound[0]<<','<<bound[1]<<','<<reference<<std::endl;
}
// Read-only H10h dual certificate for the fixed final-operator moment box.
// No pair coefficients, fitted target, second material flux or state is made.
void materialIntervalCertificate(std::ofstream& summary,std::ofstream& normals,std::ofstream& rows,
    bool rz,bool amr,int n,double rho,const CchDualMaterial& material,
    const std::vector<Vector2D>& points,const std::vector<long double>& q,
    const std::vector<Vector2D>& advector,const std::vector<int>& partners)
{
    const auto& l=material.layout();const auto& geometry=material.geometry();const auto& segments=geometry.segments();
    const long double chi=std::abs(material.material().liquidDensity()-material.material().gasDensity())/std::max(material.material().liquidDensity(),material.material().gasDensity());
    std::map<double,long double> budgetAt;long double c0=0,c1=0,gross=0,sumBudget=0,forced=0,gasForced=0,pairError=0;
    int worst=-1,gasWorst=-1;
    for(int i=0;i<l.velocityCount();++i){const long double mass=material.midpointMass().masses()[i];const auto& box=geometry.volumes()[i];const long double h=std::max(box.upper.x-box.lower.x,box.upper.y-box.lower.y);
        long double s=0;for(const auto& side:box.sides)for(int g:side)s+=std::abs(material.faceRates()[g]);const long double speed=std::abs(advector[i].x)+std::abs(advector[i].y),scale=s+mass*speed/h;
        const bool fixed=rz && component(l,i)==0 && points[i].x==0,exterior=box.lower.x==0 || box.lower.y==0 || box.upper.x==1 || box.upper.y==1;
        const auto phase=material.midpointMass().phaseVolumes()[i];const long double lower=std::abs(q[i])*h/mass;
        if(!fixed && !exterior && lower>forced){forced=lower;worst=i;}if(!fixed && !exterior && phase.liquid==0 && lower>gasForced){gasForced=lower;gasWorst=i;}
        long double gasReference=0,gasGross=0;for(int side=0;side<4;++side)for(int g:box.sides[side]){const auto& segment=segments[g];const long double value=material.material().gasDensity()*segment.area*(segment.normalAxis==0?advector[i].x:advector[i].y);gasReference+=(side%2==0?-1:1)*value;gasGross+=std::abs(value);}
        check(relative(gasReference,gasGross)<=1e-12L,"constant actual gas volume-rate closure differs");
        rows<<rz<<','<<amr<<','<<n<<','<<rho<<','<<i<<','<<component(l,i)<<','<<points[i].x<<','<<points[i].y<<','<<fixed<<','<<exterior<<','<<mass<<','<<box.volume<<','<<phase.liquid<<','<<phase.gas<<','<<h<<','<<q[i]<<','<<s<<','<<speed<<','<<scale<<','<<lower<<','<<material.cells()[i].defectRate<<','<<material.beginMass().phaseVolumes()[i].liquid<<','<<material.endMass().phaseVolumes()[i].liquid<<','<<material.cells()[i].massChange/material.deltaTime()<<','<<advector[i].x<<','<<advector[i].y<<','<<gasReference<<'\n';
        if(component(l,i)!=0)continue;const long double z=rz?points[i].x:points[i].y,b=mass*(rz?advector[i].x:advector[i].y),budget=fixed?0:chi*scale*h*h*h;
        c0-=b+z*q[i];c1-=z*b+.5L*z*z*q[i];gross+=std::abs(b)+std::abs(z*q[i])+s*h+mass*speed;budgetAt[static_cast<double>(z)]+=budget;sumBudget+=budget;}
    struct Normal {long double a,b;};std::vector<Normal> directions{{1,0}};for(const auto& entry:budgetAt)directions.push_back({-entry.first,1});
    long double maximumRatio=0,maximumViolation=0;Normal witness{};long double witnessDefect=0,witnessSupport=0;
    for(const auto& lambda:directions){const long double value=std::abs(lambda.a*c0+lambda.b*c1);long double support=0;
        for(const auto& entry:budgetAt)support+=entry.second*std::abs(lambda.a+lambda.b*entry.first);
        const long double roundoff=gross*(std::abs(lambda.a)+std::abs(lambda.b)),excess=std::max(0.L,value-support),violation=excess/std::max(1.L,support+roundoff);
        const long double ratio=support>0?value/support:(violation<=256*std::numeric_limits<double>::epsilon()?0:std::numeric_limits<long double>::infinity());
        if(violation>maximumViolation || (violation==maximumViolation && ratio>maximumRatio)){witness=lambda;witnessDefect=value;witnessSupport=support;}
        maximumRatio=std::max(maximumRatio,ratio);maximumViolation=std::max(maximumViolation,violation);
        normals<<rz<<','<<amr<<','<<n<<','<<rho<<','<<lambda.a<<','<<lambda.b<<','<<c0<<','<<c1<<','<<value<<','<<support<<','<<roundoff<<','<<ratio<<','<<violation<<'\n';}
    // The symbolic identity annihilates ANY shared pair. Check real internal
    // segments and actual declared ports as a separate finite witness.
    for(std::size_t g=0;g<segments.size();++g){const auto& segment=segments[g];if(component(l,owner(segment))!=0)continue;int i=-1,j=-1;
        if(internal(segment)){i=segment.negativeVelocity;j=segment.positiveVelocity;}else if(partners[g]>=0 && static_cast<int>(g)<partners[g]){i=owner(segment);j=owner(segments[partners[g]]);}else continue;
        const long double a=rz?points[i].x:points[i].y,b=rz?points[j].x:points[j].y,delta=b-a;
        for(const auto& lambda:directions){const long double term0=2*lambda.a*(a-b)+lambda.b*(a*a-b*b),term1=(2*lambda.a+lambda.b*(a+b))*delta;
            pairError=std::max(pairError,relative(term0+term1,std::abs(term0)+std::abs(term1)));}}
    check(pairError<=1e-12L,"final pair Farkas multiplier identity differs");
    const bool feasible=maximumViolation<=256*std::numeric_limits<double>::epsilon();
    summary<<rz<<','<<amr<<','<<n<<','<<rho<<','<<chi<<','<<c0<<','<<c1<<','<<maximumRatio<<','<<maximumViolation<<','<<witness.a<<','<<witness.b<<','<<witnessDefect<<','<<witnessSupport<<','<<feasible<<','<<pairError<<','<<forced<<','<<worst<<','<<gasForced<<','<<gasWorst<<'\n';
    std::cout<<"MATERIAL_INTERVAL_CERTIFICATE,"<<rz<<','<<amr<<','<<n<<','<<rho<<','<<feasible<<','<<maximumRatio<<','<<maximumViolation<<','<<pairError<<','<<forced<<','<<gasForced<<std::endl;
}
// Read-only H10i source/functional audit. No velocity coefficients, phase
// rates, material flux, final operator or accepted state is constructed.
struct SourceSegmentMeasure {std::array<long double,5> liquid{},gas{};long double closure=0,quadrature=0;std::size_t ranges=0,roots=0;};
SourceSegmentMeasure sourceSegmentMeasure(const std::shared_ptr<const LevelSetHermiteGeometry::Snapshot>& scalar,Vector2D a,Vector2D b,bool rz)
{
    LevelSetMappedSegment query(scalar,a,b);check(query.usesSource(*scalar),"spatial segment changed scalar source");
    const long double dx=static_cast<long double>(b.x)-a.x,dy=static_cast<long double>(b.y)-a.y,length=std::hypot(dx,dy);
    static constexpr int choose[5][5]={{1,0,0,0,0},{1,1,0,0,0},{1,2,1,0,0},{1,3,3,1,0},{1,4,6,4,1}};
    SourceSegmentMeasure result;result.ranges=query.ranges().size();result.roots=query.roots().size();
    const long double node=std::sqrt(.6L),nodes[3]={-node,0,node},weights[3]={5.L/9,8.L/9,5.L/9};
    for(const auto& range:query.ranges()){
        const long double width=range.length,sign=range.fromEnd?-1:1,t=range.fromEnd?1-range.offset:range.offset;
        const long double radialOrigin=(range.fromEnd?b.x:a.x)+(range.fromEnd?-dx:dx)*range.offset,radialSlope=sign*dx;
        for(int k=0;k<5;++k){long double moment=0,gauss=0;
            for(int j=0;j<=k;++j){const long double coefficient=choose[k][j]*std::pow(t,k-j)*std::pow(sign,j);
                moment+=coefficient*(rz?radialOrigin*std::pow(width,j+1)/(j+1)+radialSlope*std::pow(width,j+2)/(j+2):std::pow(width,j+1)/(j+1));}
            moment*=length*(rz?2*pi:1);for(int j=0;j<3;++j){const long double z=.5L*width*(1+nodes[j]);
                gauss+=weights[j]*std::pow(t+sign*z,k)*(rz?radialOrigin+radialSlope*z:1.L);}
            gauss*=.5L*width*length*(rz?2*pi:1);result.quadrature=std::max(result.quadrature,std::abs(moment-gauss));
            (range.negative?result.liquid:result.gas)[k]+=moment;}}
    const long double capacity=length*(rz?pi*(static_cast<long double>(a.x)+b.x):1.L);
    for(int k=0;k<5;++k){const long double expected=length*(rz?2*pi*(static_cast<long double>(a.x)/(k+1)+dx/(k+2)):1.L/(k+1));
        result.closure=std::max(result.closure,std::abs(result.liquid[k]+result.gas[k]-expected));}
    result.closure=capacity>0?result.closure/capacity:result.closure;result.quadrature=capacity>0?result.quadrature/capacity:result.quadrature;
    check(result.closure<=1e-12L && result.quadrature<=1e-12L,"source phase moment or GL3 identity failed");return result;
}
void spatialSourceAudit(std::ofstream& rows,std::ofstream& patches,std::ofstream& summary,std::ofstream& witness,
    const QuadtreeMesh& mesh,bool rz,bool amr,int n,const CchPhaseGeometry& begin,const CchPhaseGeometry& middle,const CchPhaseGeometry& end)
{
    const auto view=std::make_shared<const HermiteGalsGeometryView>(middle.publication().state().geometry());
    const std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> scalar(view,&view->snapshot());
    const auto include=[](std::vector<int>& list,int i){if(std::find(list.begin(),list.end(),i)==list.end())list.push_back(i);};
    for(int selection=0;selection<4;++selection){FaceSegmentGraph graph;graph.build(mesh,[selection](const QuadtreeMesh&,const Cell& c){const auto p=c.centerCoord();
            if(selection==1)return p.x<.5;if(selection==2)return !(p.x>.25 && p.x<.75 && p.y>.25 && p.y<.75);if(selection==3)return p.x<.25 || p.x>.75;return true;});
        const auto l=CchDofLayout::build(mesh,NativeLeafDomain::build(mesh,graph,rz));const auto metric=CchPhysicalDivergence::build(mesh,CchIntegralDivergence::build(mesh,l));
        const auto geometry=CchMomentumGeometry::build(mesh,metric);const auto phase=CchPhaseGeometry::build(mesh,metric,middle.publication());
        std::vector<Vector2D> lower(l.velocityCount()),upper(l.velocityCount()),means(l.velocityCount());
        std::vector<std::vector<int>> faces(l.pressureCount()),neighbors(l.pressureCount());
        for(int i=0;i<l.velocityCount();++i){const auto f=l.domain().face(mesh,l.nativeFace(i));const auto& edge=f.horizontal?mesh.hEdges(f.edgeLabel):mesh.vEdges(f.edgeLabel);
            auto a=mesh.verts(edge.edgeToVert[0]).coord(),b=mesh.verts(edge.edgeToVert[1]).coord();const int axis=component(l,i);if(axis==0?a.y>b.y:a.x>b.x)std::swap(a,b);lower[i]=a;upper[i]=b;
            means[i]={rz?(axis==0?.5*f.center.x*f.center.x:.25*(a.x*a.x+b.x*b.x)):f.center.x,f.center.y};
            const int p=f.negativeCell<0?-1:l.pressureDof(f.negativeCell),q=f.positiveCell<0?-1:l.pressureDof(f.positiveCell);
            if(p>=0)faces[p].push_back(i);if(q>=0)faces[q].push_back(i);if(p>=0 && q>=0){include(neighbors[p],q);include(neighbors[q],p);}}
        using Key=std::tuple<int,double,double,double>;std::map<Key,SourceSegmentMeasure> cache;long double peak=0,quadrature=0,missingPeak=0;int deficient=0,missing=0;
        const auto observe=[&](const char* kind,int id,int comp,int axis,Vector2D a,Vector2D b,long double area){
            const Key key{axis,axis==0?a.x:a.y,axis==0?a.y:a.x,axis==0?b.y:b.x};auto found=cache.find(key);
            if(found==cache.end())found=cache.emplace(key,sourceSegmentMeasure(scalar,a,b,rz)).first;const auto& m=found->second;
            const long double areaError=area>0?std::abs(m.liquid[0]+m.gas[0]-area)/area:std::abs(m.liquid[0]+m.gas[0]-area);
            check(areaError<=1e-12L,"retained phase source/physical segment measure differs");peak=std::max({peak,m.closure,areaError});quadrature=std::max(quadrature,m.quadrature);
            rows<<rz<<','<<amr<<','<<n<<','<<selection<<','<<kind<<','<<id<<','<<comp<<','<<axis<<','<<a.x<<','<<a.y<<','<<b.x<<','<<b.y<<','<<area<<','<<middle.geometryIdentity().value()<<','<<m.ranges<<','<<m.roots;
            for(int k=0;k<5;++k)rows<<','<<m.liquid[k]<<','<<m.gas[k];rows<<','<<m.closure<<','<<m.quadrature<<'\n';};
        for(int i=0;i<l.velocityCount();++i){observe("primal",i,component(l,i),component(l,i),lower[i],upper[i],metric.faceMeasures()[i]);const auto face=l.domain().face(mesh,l.nativeFace(i));
            std::vector<int> root,patch,support;for(int p:{face.negativeCell,face.positiveCell})if(p>=0)include(root,l.pressureDof(p));patch=root;for(int p:root)for(int q:neighbors[p])include(patch,q);
            for(int p:patch)for(int j:faces[p])if(component(l,j)==component(l,i) && !(rz && component(l,j)==0 && metric.faceMeasures()[j]==0))include(support,j);
            const double hx=rz?.5*(upper[i].x*upper[i].x-lower[i].x*lower[i].x):face.normalSpacing,hy=face.tangentSpacing;
            const double scaleX=rz?std::max(hx,.5*face.normalSpacing*face.normalSpacing+face.center.x*face.normalSpacing):std::max(face.normalSpacing,face.tangentSpacing),scaleY=std::max(hy,face.normalSpacing);
            if(rz && component(l,i)==0 && metric.faceMeasures()[i]==0){patches<<rz<<','<<amr<<','<<n<<','<<selection<<','<<i<<','<<component(l,i)<<','<<support.size()<<",0,0,0,1\n";continue;}
            check(!support.empty(),"source has no actual local face data");Eigen::MatrixXd a(3,support.size());for(std::size_t j=0;j<support.size();++j){const auto d=means[support[j]]-means[i];a(0,j)=1;a(1,j)=d.x/scaleX;a(2,j)=d.y/scaleY;}
            Eigen::JacobiSVD<Eigen::MatrixXd> svd(a,Eigen::ComputeThinU|Eigen::ComputeThinV);svd.setThreshold(1e-12);const Eigen::Vector3d target=component(l,i)==0?Eigen::Vector3d{0,0,1}:Eigen::Vector3d{0,1,0};
            const auto weights=svd.solve(target).eval();const double deficit=(a*weights-target).norm();if(svd.rank()<3)++deficient;if(deficit>1e-12)++missing;missingPeak=std::max(missingPeak,static_cast<long double>(deficit));
            patches<<rz<<','<<amr<<','<<n<<','<<selection<<','<<i<<','<<component(l,i)<<','<<support.size()<<','<<svd.rank()<<','<<deficit<<','<<weights.norm()<<",0\n";}
        for(std::size_t g=0;g<geometry.segments().size();++g){const auto& s=geometry.segments()[g];observe("dual",static_cast<int>(g),component(l,owner(s)),s.normalAxis,s.begin,s.end,s.area);}
        if(!rz && amr && n==16 && selection==0){const int u=103;for(int layer=0;layer<3;++layer){const auto& stage=layer==0?begin:layer==1?middle:end;const auto retained=std::make_shared<const HermiteGalsGeometryView>(stage.publication().state().geometry());
                const std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> actual(retained,&retained->snapshot());for(int side=0;side<4;++side)for(int g:geometry.volumes()[u].sides[side]){const auto& s=geometry.segments()[g];const auto m=sourceSegmentMeasure(actual,s.begin,s.end,rz);
                    check(m.liquid[0]==0,"three-layer pure gas witness has a liquid boundary interval");witness<<layer<<','<<stage.geometryIdentity().value()<<','<<u<<','<<side<<','<<g<<','<<m.liquid[0]<<','<<m.gas[0]<<','<<m.closure<<'\n';}}}
        summary<<rz<<','<<amr<<','<<n<<','<<selection<<','<<l.velocityCount()<<','<<geometry.segments().size()<<','<<cache.size()<<','<<peak<<','<<quadrature<<','<<deficient<<','<<missing<<','<<missingPeak<<'\n';
        std::cout<<"SPATIAL_SOURCE_AUDIT,"<<rz<<','<<amr<<','<<n<<','<<selection<<','<<l.velocityCount()<<','<<geometry.segments().size()<<','<<cache.size()<<','<<peak<<','<<quadrature<<','<<deficient<<','<<missing<<std::endl;}
}

// Passive H10i exact-plane source capability check, not a flux convention.
void coincidentSourceAudit()
{
    std::ofstream out("output/series/coincident_source.csv");check(out.good(),"cannot open coincident source audit");
    out<<std::setprecision(17)<<"rz,amr,geometry_identity,publication_ok,dual_mass_ok,queried_segments,unsupported_zero_traces,offset,mass_reason,segment_reason\n";
    for(bool rz:{false,true})for(bool amr:{false,true}){auto mesh=makeMesh(4,amr?1:0);Case shape;shape.pureFirstJet=true;shape.radial=rz;shape.planar=true;shape.planeOffset=amr?.125:.375;
        const auto pub=initial(*mesh,shape);FaceSegmentGraph graph;graph.build(*mesh);const auto l=CchDofLayout::build(*mesh,NativeLeafDomain::build(*mesh,graph,rz));
        const auto metric=CchPhysicalDivergence::build(*mesh,CchIntegralDivergence::build(*mesh,l));const auto geometry=CchMomentumGeometry::build(*mesh,metric);
        const auto phase=CchPhaseGeometry::build(*mesh,metric,pub);bool massOK=true;std::string massReason;
        try{const auto mass=CchMomentumMass::build(*mesh,geometry,phase,1000,1);(void)mass;}catch(const std::exception& e){massOK=false;massReason=e.what();}
        const auto view=std::make_shared<const HermiteGalsGeometryView>(pub.state().geometry());const std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> scalar(view,&view->snapshot());
        int zero=0;std::string reason;for(const auto& g:geometry.segments()){
            try{LevelSetMappedSegment query(scalar,g.begin,g.end);check(query.usesSource(*scalar),"coincident query changed scalar");}
            catch(const std::exception& e){const std::string detail=e.what();check(detail.find("positive-length coincident zero trace is unsupported")!=std::string::npos,"unexpected exact-plane segment rejection");++zero;reason=detail;}}
        check(zero>0,"exact-plane source did not exercise a coincident dual segment");
        std::replace(massReason.begin(),massReason.end(),',',';');std::replace(reason.begin(),reason.end(),',',';');
        out<<rz<<','<<amr<<','<<pub.state().identity().value()<<",1,"<<massOK<<','<<geometry.segments().size()<<','<<zero<<','<<shape.planeOffset<<','<<massReason<<','<<reason<<'\n';
        std::cout<<"COINCIDENT_SOURCE,"<<rz<<','<<amr<<','<<massOK<<','<<geometry.segments().size()<<','<<zero<<','<<reason<<std::endl;}
}

// Read-only H10i1a port observation on the actual published native scalar.
// Closed-form planar references are diagnostics, never phase inputs.
void tracePortObservation(std::ofstream& out,std::ofstream& ranges,
    const std::shared_ptr<const LevelSetHermiteGeometry::Snapshot>& scalar,
    const CchMomentumSegment& g,bool rz,bool amr,int sign,double plane,std::size_t segment,std::uint64_t identity)
{
    auto q=LevelSetMappedSegment::nativeQuarticTrace(scalar,g.begin,g.end,g.normalAxis);
    check(q.usesSource(*scalar),"new trace port changed scalar");
    const auto area=[&](const auto& r){
        if(!rz || g.normalAxis==0)return g.area*r.length;
        const long double delta=static_cast<long double>(g.end.x)-g.begin.x;
        const long double center=r.fromEnd?g.end.x-delta*(r.offset+r.length/2):g.begin.x+delta*(r.offset+r.length/2);
        return 2*pi*std::abs(delta)*r.length*center;
    };
    long double ordinary=0,liquid=0,zero=0,minus=0,plus=0;
    for(const auto& r:q.ordinaryRanges()){ordinary+=area(r);if(r.negative)liquid+=area(r);}
    for(const auto& z:q.coincidentRanges()){
        zero+=area(z.interval);
        for(int side:{-1,1}){const auto& s=side<0?z.minus:z.plus;
            check(s.available && scalar->nativeCellActive(s.nativeCellIndex) && s.normalOrder==1,"actual planar normal trace provenance");
            for(const auto& r:s.ranges){if(r.negative)(side<0?minus:plus)+=area(r);
                ranges<<rz<<','<<amr<<','<<sign<<','<<segment<<','<<identity<<','<<z.chart<<','<<side<<','<<s.nativeCellIndex<<','<<s.normalOrder<<','<<r.fromEnd<<','<<r.offset<<','<<r.length<<','<<r.negative<<'\n';}}
    }
    const double first=rz?g.begin.y:g.begin.x,last=rz?g.end.y:g.end.x;
    const bool coincident=first==plane && last==plane;
    std::size_t pointOnly=0;
    for(std::size_t k=0;k<q.events().size();++k){const auto& e=q.events()[k];if(e.root.pointOnly)++pointOnly;
        if(coincident)check(e.normalOrder==0 && e.root.pointOnly && (q.rootAtInputEndpoint(k,false) || q.rootAtInputEndpoint(k,true)),
            "zero trace invented a finite value root instead of a point-only native endpoint incidence");}
    long double reference=0;
    if(!coincident){
        if(first==last){if(sign*(first-plane)<0)reference=g.area;}
        else {const long double root=std::clamp((static_cast<long double>(plane)-first)/(last-first),0.L,1.L);
            const bool prefix=sign*(last-first)>0;const long double lo=prefix?0:root,hi=prefix?root:1;
            if(rz && g.normalAxis==1){const long double dr=static_cast<long double>(g.end.x)-g.begin.x;
                reference=2*pi*std::abs(dr)*(hi-lo)*(g.begin.x+dr*(lo+hi)/2);}
            else reference=g.area*(hi-lo);}
        LevelSetMappedSegment old(scalar,g.begin,g.end);
        check(q.ordinaryRanges().size()==old.ranges().size() && q.ordinary().roots().size()==old.roots().size(),"ordinary native port changed inventory");
        for(std::size_t k=0;k<old.ranges().size();++k){const auto& a=old.ranges()[k];const auto& b=q.ordinaryRanges()[k];
            check(a.fromEnd==b.fromEnd && a.offset==b.offset && a.length==b.length && a.negative==b.negative,"ordinary native port changed range");}
    }
    const long double expectedMinus=coincident && sign>0?g.area:0,expectedPlus=coincident && sign<0?g.area:0;
    const auto error=std::max({std::abs(ordinary+zero-g.area),std::abs(liquid-reference),std::abs(minus-expectedMinus),std::abs(plus-expectedPlus)});
    check(error<=1e-12L*g.area,"actual planar trace port measure mismatch");
    check(coincident==!q.coincidentRanges().empty(),"actual planar zero kind mismatch");
    out<<rz<<','<<amr<<','<<sign<<','<<segment<<','<<identity<<','<<g.normalAxis<<','<<g.begin.x<<','<<g.begin.y<<','<<g.end.x<<','<<g.end.y<<','<<g.area
        <<','<<plane<<','<<q.ordinaryRanges().size()<<','<<q.coincidentRanges().size()<<','<<q.events().size()<<','<<pointOnly<<','<<ordinary<<','<<liquid<<','<<zero<<','<<minus<<','<<plus<<','<<error<<'\n';
}

// Passive H10i1 capability inventory. References are closed-form rectangle
// volumes only; they never supply a phase trace, material rate or geometry.
void traceContractAudit()
{
    // H10j observer begin
    std::ofstream spatialZero("output/series/spatial_zero_rate.csv");
    spatialZero<<std::setprecision(17)<<"rz,amr,location,phase_sign,offset,geometry,mode,advector,atom,axis,ax,ay,bx,by,area,liquid,gas";
    writeCoefficientHeader<decltype(CchNormalPolynomial::beginCoefficients)>(spatialZero,"c");spatialZero<<",total,coincident\n";
    // H10j observer end
    // H10i2 observer begin
    std::ofstream sourceFaces("output/series/source_chain_faces.csv"),sourceDuals("output/series/source_chain_duals.csv"),sourceOverlaps("output/series/source_chain_overlaps.csv");
    check(sourceFaces.good() && sourceDuals.good() && sourceOverlaps.good(),"cannot open independent source chain coordinates");
    const auto provenanceHeader="rz,amr,location,phase_sign,offset,geometry_identity,step,clock_begin,clock_end,position,dof,";
    sourceFaces<<std::setprecision(17)<<provenanceHeader<<"native_face,normal_axis,begin_x,begin_y,end_x,end_y,area,liquid,gas,kind,negative_cell,positive_cell,minus_liquid,minus_gas,plus_liquid,plus_gas,zero_halves,same_snapshot,no_inventory\n";
    sourceDuals<<std::setprecision(17)<<provenanceHeader<<"lower_x,lower_y,upper_x,upper_y,capacity,liquid,gas,mass,same_snapshot,no_inventory\n";
    sourceOverlaps<<std::setprecision(17)<<provenanceHeader<<"pressure,native_cell,lower_x,lower_y,upper_x,upper_y,capacity,cell_fraction,liquid,gas\n";
    // H10i2 observer end
    std::ofstream spatialFaces("output/series/spatial_face_trace.csv"),regionSource("output/series/native_region_source.csv");
    check(spatialFaces.good() && regionSource.good(),"cannot open native phase consumers audit");
    spatialFaces<<std::setprecision(17)<<"rz,amr,location,phase_sign,geometry_identity,step,time,dof,native_face,kind,negative_cell,positive_cell,zero_halves,area,liquid,gas,minus_liquid,minus_gas,plus_liquid,plus_gas,expected_liquid,relative_error\n";
    regionSource<<std::setprecision(17)<<"rz,amr,location,phase_sign,geometry_identity,step,time,components,side_graphs,capacity,liquid,gas,axial_liquid,radial_liquid,expected_liquid,integrated_arcs,relative_error,same_source\n";
    std::ofstream publicationRows("output/series/publication_trace_cells.csv"),publicationSummary("output/series/publication_trace_summary.csv");
    check(publicationRows.good() && publicationSummary.good(),"cannot open actual trace publication audit");
    publicationRows<<std::setprecision(17)<<"rz,amr,location,phase_sign,geometry_identity,step,time,cell,capacity,liquid,expected_liquid,axial_liquid,radial_liquid,relative_error\n";
    publicationSummary<<std::setprecision(17)<<"rz,amr,location,phase_sign,geometry_identity,step,time,primary_cells,components,graphs,side_graphs,side_liquid_owners,zero_halves,total_liquid,expected_liquid,max_relative_error,max_green_defect,max_capacity_defect,same_source,no_inventory\n";
    std::ofstream out("output/series/trace_contract.csv"),boxes("output/series/trace_clip.csv");
    std::ofstream masses("output/series/trace_mass.csv");
    std::ofstream port("output/series/native_trace_port.csv"),ranges("output/series/native_trace_ranges.csv");
    check(out.good() && boxes.good() && masses.good(),"cannot open trace contract audit");
    check(port.good() && ranges.good(),"cannot open native trace port audit");
    out<<std::setprecision(17)<<"rz,amr,location,phase_sign,offset,publication_ok,phase_ok,dual_mass_ok,queried_segments,unsupported_zero_traces,failure_stage,reason\n";
    boxes<<std::setprecision(17)<<"rz,amr,phase_sign,side,offset,capacity,expected_liquid,expected_gas,measure_ok,liquid,gas,axial_liquid,radial_liquid,axial_gas,radial_gas,error_estimate,geometry_identity,reason\n";
    masses<<std::setprecision(17)<<"rz,amr,phase_sign,dof,geometry_identity,capacity,liquid,gas,expected_liquid,expected_gas,mass,expected_mass\n";
    port<<std::setprecision(17)<<"rz,amr,phase_sign,segment,geometry_identity,normal_axis,begin_x,begin_y,end_x,end_y,area,plane,ordinary_ranges,zero_charts,events,point_only_events,ordinary_area,ordinary_liquid,zero_area,minus_liquid,plus_liquid,error\n";
    ranges<<std::setprecision(17)<<"rz,amr,phase_sign,segment,geometry_identity,chart,side,native_cell,normal_order,from_end,offset,length,negative\n";
    const auto clean=[](std::string s){for(char& c:s)if(c==',' || c=='\n' || c=='\r')c=';';return s;};
    for(bool rz:{false,true})for(bool amr:{false,true})for(int location=0;location<2;++location)for(int sign:{-1,1}){
        auto mesh=makeMesh(4,amr?1:0);Case shape;shape.pureFirstJet=true;shape.radial=rz;shape.planar=true;shape.planeSign=sign;
        shape.planeOffset=location==0?(amr?.125:.375):(amr?.5:.25);
        bool publicationOK=false,phaseOK=false,massOK=false;std::size_t queried=0;int zero=0;
        std::string stage="publication",reason;
        try{
            const auto pub=initial(*mesh,shape);publicationOK=true;stage="phase";
            const auto& actualScalar=pub.state().geometry().snapshot();const auto& actualLedger=pub.volumes().nativeTraceLedger();
            const auto actualTime=pub.state().timeLayer();long double actualPeak=0;
            const bool sameSource=&pub.facePartitions().source()==&actualScalar && &pub.volumes().geometry().snapshot()==&actualScalar
                && pub.volumes().sourceIdentity()==pub.state().identity() && pub.volumes().timeLayer().time()==actualTime.time();
            check(sameSource && pub.state().primary().phaseVolumes.empty(),"actual publication scalar/revision/time/inventory mismatch");
            for(const auto& cell:pub.volumes().cells()) {
                const auto bounds=actualScalar.nativeCellBounds(cell.cell);const auto lo=bounds[0],hi=bounds[1];
                const double lower=rz?lo.y:lo.x,upper=rz?hi.y:hi.x;
                const long double length=std::max(0.,std::min(shape.planeOffset,upper)-lower);
                const auto cross=rz?pi*(static_cast<long double>(hi.x)-lo.x)*(static_cast<long double>(hi.x)+lo.x):static_cast<long double>(hi.y)-lo.y;
                const auto expected=sign>0?length*cross:cell.capacity-length*cross;
                const auto axial=actualLedger.axialVolumeByCell[cell.cell],radial=actualLedger.radialVolumeByCell[cell.cell];
                const auto trial=std::max({std::abs(cell.rawLiquid-expected),std::abs(axial-expected),std::abs(radial-expected)})/cell.capacity;
                check(trial<=1e-12L,"actual trace publication independent planar volume");actualPeak=std::max(actualPeak,trial);
                publicationRows<<rz<<','<<amr<<','<<location<<','<<sign<<','<<pub.state().identity().value()<<','<<actualTime.step<<','<<actualTime.time()<<','<<cell.cell<<','<<cell.capacity
                    <<','<<cell.rawLiquid<<','<<expected<<','<<axial<<','<<radial<<','<<trial<<'\n';
            }
            const auto expectedTotal=(sign>0?shape.planeOffset:1-shape.planeOffset)*(rz?pi:1.L);
            publicationSummary<<rz<<','<<amr<<','<<location<<','<<sign<<','<<pub.state().identity().value()<<','<<actualTime.step<<','<<actualTime.time()<<','<<pub.state().primary().geometryCells.size()
                <<','<<pub.report().components<<','<<pub.report().graphs<<','<<actualLedger.sideGraphs<<','<<actualLedger.sideLiquidOwners<<','<<actualLedger.zeroHalfFaces
                <<','<<pub.volumes().report().totalRawLiquid<<','<<expectedTotal<<','<<actualPeak<<','<<actualLedger.maxRelativeGreenDefect<<','<<actualLedger.maxRelativeCapacityDefect<<','<<sameSource<<",1\n";
            std::cout<<"TRACE_PUBLICATION,"<<rz<<','<<amr<<','<<location<<','<<sign<<','<<actualLedger.sideGraphs<<','<<actualLedger.zeroHalfFaces<<','<<actualPeak<<",PASS\n";
            FaceSegmentGraph graph;graph.build(*mesh);const auto l=CchDofLayout::build(*mesh,NativeLeafDomain::build(*mesh,graph,rz));
            const auto metric=CchPhysicalDivergence::build(*mesh,CchIntegralDivergence::build(*mesh,l));
            const auto phase=CchPhaseGeometry::build(*mesh,metric,pub);phaseOK=true;stage="mass";
            for(int u=0;u<phase.layout().velocityCount();++u) {
                const auto face=phase.layout().domain().face(*mesh,phase.layout().nativeFace(u));const auto& t=phase.faceTraces()[u];const auto m=phase.faceMeasures()[u];
                const auto& edge=face.horizontal?mesh->hEdges(face.edgeLabel):mesh->vEdges(face.edgeLabel);
                const auto a=mesh->verts(edge.edgeToVert[0]).coord(),b=mesh->verts(edge.edgeToVert[1]).coord();
                const double x=rz?a.y:a.x,dx=rz?b.y-a.y:b.x-a.x;const bool coincident=dx==0 && x==shape.planeOffset;
                long double expected=0;
                if(coincident){check(t.kind==CchPhaseFaceTraceKind::InteriorCoincident && t.minusAvailable() && t.plusAvailable(),"actual native face lost both spatial sides");
                    check(std::abs(t.minus.liquid-(sign>0?face.area:0))<=1e-12L*face.area && std::abs(t.plus.liquid-(sign<0?face.area:0))<=1e-12L*face.area,"actual native side physical reference");expected=face.area*.5L;}
                else{check(t.kind==CchPhaseFaceTraceKind::Ordinary,"ordinary native face mislabeled");
                    if(dx==0)expected=sign*(x-shape.planeOffset)<0?face.area:0;
                    else{const long double crossing=std::clamp((static_cast<long double>(shape.planeOffset)-x)/dx,0.L,1.L),lo=sign>0?0:crossing,hi=sign>0?crossing:1;
                        if(rz && face.horizontal){const long double dr=static_cast<long double>(b.x)-a.x;expected=2*pi*dr*(hi-lo)*(a.x+dr*(lo+hi)*.5L);}
                        else expected=face.area*(hi-lo);}}
                const auto error=face.area>0?std::abs(m.liquid-expected)/face.area:std::abs(m.liquid-expected);check(error<=1e-12L,"actual native spatial face independent reference");
                spatialFaces<<rz<<','<<amr<<','<<location<<','<<sign<<','<<phase.geometryIdentity().value()<<','<<phase.timeLayer().step<<','<<phase.timeLayer().time()<<','<<u<<','<<phase.layout().nativeFace(u)
                    <<','<<static_cast<int>(t.kind)<<','<<t.negativeCell<<','<<t.positiveCell<<','<<t.coincidentHalfFaces<<','<<face.area<<','<<m.liquid<<','<<m.gas
                    <<','<<t.minus.liquid<<','<<t.minus.gas<<','<<t.plus.liquid<<','<<t.plus.gas<<','<<expected<<','<<error<<'\n';
                // H10i2 observer begin
                sourceFaces<<rz<<','<<amr<<','<<location<<','<<sign<<','<<shape.planeOffset<<','<<phase.geometryIdentity().value()<<','<<actualTime.step<<','<<actualTime.begin<<','<<actualTime.end<<','<<static_cast<int>(actualTime.position)<<','<<u
                    <<','<<phase.layout().nativeFace(u)<<','<<(face.horizontal?1:0)<<','<<a.x<<','<<a.y<<','<<b.x<<','<<b.y<<','<<face.area<<','<<m.liquid<<','<<m.gas
                    <<','<<static_cast<int>(t.kind)<<','<<t.negativeCell<<','<<t.positiveCell<<','<<t.minus.liquid<<','<<t.minus.gas<<','<<t.plus.liquid<<','<<t.plus.gas<<','<<t.coincidentHalfFaces
                    <<','<<(&phase.publication().state().geometry().snapshot()==&pub.state().geometry().snapshot())<<','<<pub.state().primary().phaseVolumes.empty()<<'\n';
                // H10i2 observer end
            }
            const auto geometry=CchMomentumGeometry::build(*mesh,metric);
            // H10j observer begin
            const auto spatialReconstruction=CchAdvectorReconstruction::prepare(*mesh,geometry);
            for(int mode:{-1,0,1,2}){
                Field velocity(l.velocityCount());for(int u=0;u<l.velocityCount();++u){const auto f=l.domain().face(*mesh,l.nativeFace(u));
                    if(component(l,u)!=(rz?1:0))continue;
                    if(mode!=2)velocity[u]=.4*mode;
                    else if(!rz)velocity[u]=f.center.y-.375;
                    else {const auto& e=mesh->hEdges(f.edgeLabel);const auto a=mesh->verts(e.edgeToVert[0]).coord(),b=mesh->verts(e.edgeToVert[1]).coord();velocity[u]=.25*(a.x*a.x+b.x*b.x)-.015625;}}
                const auto advector=spatialReconstruction.evaluate(*mesh,CchNativeValues(l,velocity,Field(l.pressureCount())),phase.geometryIdentity(),phase.timeLayer());
                const auto rate=CchPhaseRate::prepareSpatial(*mesh,phase,advector);
                for(std::size_t i=0;i<spatialReconstruction.atoms().size();++i){const auto& a=spatialReconstruction.atoms()[i];const auto& p=advector.atomPolynomials()[i];const auto r=rate.atomRates()[i];
                    const bool zero=a.normalAxis==(rz?1:0) && (rz?a.begin.y:a.begin.x)==shape.planeOffset;
                    spatialZero<<rz<<','<<amr<<','<<location<<','<<sign<<','<<shape.planeOffset<<','<<phase.geometryIdentity().value()<<','<<mode<<','<<advector.identity().value()<<','<<i<<','<<a.normalAxis<<','<<a.begin.x<<','<<a.begin.y<<','<<a.end.x<<','<<a.end.y<<','<<a.area<<','<<r.liquid<<','<<r.gas;writeCoefficients(spatialZero,p.beginCoefficients);spatialZero<<','<<p.integral(false,0,p.length)<<','<<zero<<'\n';}
            }
            // H10j observer end
            try{const auto mass=CchMomentumMass::build(*mesh,geometry,phase,1000,1);
                for(int u=0;u<mass.layout().velocityCount();++u){const auto& volume=geometry.volumes()[u];
                    auto lo=volume.lower,hi=volume.upper;
                    const long double capacity=(static_cast<long double>(hi.x)-lo.x)*(static_cast<long double>(hi.y)-lo.y)*(rz?pi*(static_cast<long double>(hi.x)+lo.x):1.L);
                    if(rz){if(sign>0)hi.y=std::min(hi.y,shape.planeOffset);else lo.y=std::max(lo.y,shape.planeOffset);}
                    else{if(sign>0)hi.x=std::min(hi.x,shape.planeOffset);else lo.x=std::max(lo.x,shape.planeOffset);}
                    const long double vl=lo.x>=hi.x || lo.y>=hi.y?0:(static_cast<long double>(hi.x)-lo.x)*(static_cast<long double>(hi.y)-lo.y)*(rz?pi*(static_cast<long double>(hi.x)+lo.x):1.L);
                    const auto actual=mass.phaseVolumes()[u];const auto expected=1000*vl+capacity-vl;
                    check(actual.liquid>=0 && actual.gas>=0 && std::abs(actual.liquid-vl)<=1e-12L*capacity
                        && std::abs(actual.gas-(capacity-vl))<=1e-12L*capacity
                        && std::abs(mass.masses()[u]-expected)<=1e-12L*1001*capacity,"actual coincident dual mass closed reference");
                    // H10i2 observer begin
                    const auto massClock=mass.timeLayer();
                    sourceDuals<<rz<<','<<amr<<','<<location<<','<<sign<<','<<shape.planeOffset<<','<<mass.geometryIdentity().value()<<','<<massClock.step<<','<<massClock.begin<<','<<massClock.end<<','<<static_cast<int>(massClock.position)<<','<<u
                        <<','<<volume.lower.x<<','<<volume.lower.y<<','<<volume.upper.x<<','<<volume.upper.y<<','<<volume.volume<<','<<actual.liquid<<','<<actual.gas<<','<<mass.masses()[u]
                        <<','<<(&mass.phaseGeometry().publication().state().geometry().snapshot()==&pub.state().geometry().snapshot())<<','<<pub.state().primary().phaseVolumes.empty()<<'\n';
                    for(std::size_t k=0;k<volume.overlaps.size();++k){const auto& piece=volume.overlaps[k];const auto value=mass.overlapPhaseVolumes()[u][k];
                        sourceOverlaps<<rz<<','<<amr<<','<<location<<','<<sign<<','<<shape.planeOffset<<','<<mass.geometryIdentity().value()<<','<<massClock.step<<','<<massClock.begin<<','<<massClock.end<<','<<static_cast<int>(massClock.position)<<','<<u
                            <<','<<piece.pressure<<','<<mass.layout().nativeCell(piece.pressure)<<','<<piece.lower.x<<','<<piece.lower.y<<','<<piece.upper.x<<','<<piece.upper.y<<','<<piece.volume<<','<<piece.cellFraction<<','<<value.liquid<<','<<value.gas<<'\n';}
                    // H10i2 observer end
                    masses<<rz<<','<<amr<<','<<sign<<','<<u<<','<<phase.geometryIdentity().value()<<','<<capacity<<','<<actual.liquid<<','<<actual.gas<<','<<vl<<','<<(capacity-vl)<<','<<mass.masses()[u]<<','<<expected<<'\n';}
                massOK=true;stage="none";}
            catch(const std::exception& e){reason=clean(e.what());}
            const auto view=std::make_shared<const HermiteGalsGeometryView>(pub.state().geometry());
            const std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> scalar(view,&view->snapshot());
            for(const auto& g:geometry.segments()){
                tracePortObservation(port,ranges,scalar,g,rz,amr,sign,shape.planeOffset,queried,phase.geometryIdentity().value());
                ++queried;try{LevelSetMappedSegment query(scalar,g.begin,g.end);check(query.usesSource(*scalar),"trace inventory changed scalar");}
                catch(const std::exception& e){check(std::string(e.what()).find("positive-length coincident zero trace is unsupported")!=std::string::npos,"unexpected trace inventory rejection");++zero;}}
            HermiteGalsPhaseRegionConfig config;config.symmetryPlanes=shape.symmetry(*mesh);
            HermiteGalsPhaseRegionSource source(*mesh,pub.state(),config);
            const auto whole=source.measure(*mesh,{0,0},{1,1});const auto expectedWhole=(sign>0?shape.planeOffset:1-shape.planeOffset)*(rz?pi:1.L);
            const bool regionSame=&source.source().geometry().snapshot()==&pub.state().geometry().snapshot() && whole.sourceIdentity==phase.geometryIdentity();
            const auto regionError=std::max({std::abs(whole.liquid-expectedWhole),std::abs(whole.rawAxialLiquid-expectedWhole),std::abs(whole.rawRadialLiquid-expectedWhole)})/whole.capacity;
            check(regionSame && regionError<=1e-12L,"actual native region source independent reference");
            regionSource<<rz<<','<<amr<<','<<location<<','<<sign<<','<<whole.sourceIdentity.value()<<','<<source.source().timeLayer().step<<','<<source.source().timeLayer().time()
                <<','<<source.components()<<','<<source.nativeSideGraphs()<<','<<whole.capacity<<','<<whole.liquid<<','<<whole.gas<<','<<whole.rawAxialLiquid<<','<<whole.rawRadialLiquid
                <<','<<expectedWhole<<','<<whole.integratedArcs<<','<<regionError<<','<<regionSame<<'\n';
            for(int side:{-1,1}){
                Vector2D lo{0,0},hi{1,1};const double a=shape.planeOffset;
                if(rz){lo.y=side<0?a-.125:a;hi.y=side<0?a:a+.125;}
                else{lo.x=side<0?a-.125:a;hi.x=side<0?a:a+.125;}
                const long double capacity=(static_cast<long double>(hi.x)-lo.x)*(static_cast<long double>(hi.y)-lo.y)*(rz?pi*(static_cast<long double>(hi.x)+lo.x):1.L);
                const bool liquid=sign*side<0;bool measured=false;long double vl=0,vg=0,al=0,rl=0,ag=0,rg=0,error=0;std::string detail;
                try{const auto value=source.measure(*mesh,lo,hi);measured=true;vl=value.liquid;vg=value.gas;
                    al=value.rawAxialLiquid;rl=value.rawRadialLiquid;ag=value.rawAxialGas;rg=value.rawRadialGas;error=value.momentErrorEstimate;}
                catch(const std::exception& e){detail=clean(e.what());}
                boxes<<rz<<','<<amr<<','<<sign<<','<<side<<','<<a<<','<<capacity<<','<<(liquid?capacity:0)<<','<<(liquid?0:capacity)<<','<<measured<<','<<vl<<','<<vg<<','<<al<<','<<rl<<','<<ag<<','<<rg<<','<<error<<','<<phase.geometryIdentity().value()<<','<<detail<<'\n';
            }
        }catch(const std::exception& e){reason=clean(e.what());}
        out<<rz<<','<<amr<<','<<location<<','<<sign<<','<<shape.planeOffset<<','<<publicationOK<<','<<phaseOK<<','<<massOK<<','<<queried<<','<<zero<<','<<stage<<','<<reason<<'\n';
        std::cout<<"TRACE_CONTRACT,"<<rz<<','<<amr<<','<<location<<','<<sign<<','<<publicationOK<<','<<phaseOK<<','<<massOK<<','<<queried<<','<<zero<<','<<stage<<','<<reason<<std::endl;
    }
    // Independent, explicit Core capability after the original observations:
    // a seed source is not an accepted Interface/Hydro publication.
    std::ofstream topo("output/series/trace_topology.csv"),sides("output/series/native_side_support.csv");
    std::ofstream cacheRows("output/series/half_trace_cache.csv"),cacheSummary("output/series/half_trace_cache_summary.csv");
    check(topo.good() && sides.good(),"cannot open native topology observation");
    check(cacheRows.good() && cacheSummary.good(),"cannot open native half trace cache observation");
    topo<<std::setprecision(17)<<"rz,amr,location,phase_sign,offset,seed_identity,core_ok,curves,graphs,side_graphs,length,physical_measure,radial_moment,axial_moment,error,reason\n";
    sides<<std::setprecision(17)<<"rz,amr,phase_sign,seed_identity,graph,owner,normal_axis,edge_label,begin_vertex,end_vertex,minus_cell,minus_phase,plus_cell,plus_phase,begin_x,begin_y,end_x,end_y\n";
    cacheRows<<std::setprecision(17)<<"rz,amr,location,phase_sign,seed_identity,face,half,normal_axis,begin_x,begin_y,end_x,end_y,area,ordinary_ranges,zero_charts,value_events,normal_events,liquid,expected_liquid,minus_liquid,plus_liquid,error\n";
    cacheSummary<<std::setprecision(17)<<"rz,amr,location,phase_sign,seed_identity,halves,zero_halves,value_events,normal_events,maximum_area_error\n";
    for(bool rz:{false,true})for(bool amr:{false,true})for(int location=0;location<2;++location)for(int sign:{-1,1}){
        auto mesh=makeMesh(4,amr?1:0);Case shape;shape.pureFirstJet=true;shape.radial=rz;shape.planar=true;shape.planeSign=sign;
        shape.planeOffset=location==0?(amr?.125:.375):(amr?.5:.25);
        const auto seed=stateFor(*mesh,shape.profile(0),rz,{},[](int,Vector2D,Vector2D){return 0.;},shape.symmetry(*mesh),HermiteGalsTransport::FirstJet);
        check(seed.primary().phaseVolumes.empty(),"Core observer minted a material inventory");
        const auto view=std::make_shared<const HermiteGalsGeometryView>(seed.geometry());
        const std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> scalar(view,&view->snapshot());
        bool ok=false;std::size_t curves=0,graphs=0,sideCount=0;long double length=0,measure=0,radial=0,axial=0,err=0;std::string detail;
        try{const auto result=LevelSetInterfacePanels::discoverNativeTraceComponents(*mesh,scalar,rz,{},shape.symmetry(*mesh));curves=result.curves.size();
            if(location==0){const auto old=LevelSetInterfacePanels::discoverNativeComponents(*mesh,scalar,rz,{},shape.symmetry(*mesh));
                check(old.curves.size()==result.curves.size() && old.cellPhase==result.cellPhase,"ordinary Core topology inventory changed");
                for(std::size_t c=0;c<old.curves.size();++c){const auto& a=old.curves[c];const auto& b=result.curves[c];check(a.graphs().size()==b.graphs().size(),"ordinary Core graph count changed");
                    for(int k=0;k<static_cast<int>(a.graphs().size());++k){const auto& x=a.graphs()[k];const auto& y=b.graphs()[k];
                        check(!y.nativeSide && x.nativeCellIndex==y.nativeCellIndex && x.coordinateIsX==y.coordinateIsX && x.begin.side==y.begin.side && x.end.side==y.end.side
                            && x.begin.rootOrdinal==y.begin.rootOrdinal && x.end.rootOrdinal==y.end.rootOrdinal && x.localLower.x==y.localLower.x && x.localLower.y==y.localLower.y
                            && x.localUpper.x==y.localUpper.x && x.localUpper.y==y.localUpper.y,"ordinary retained graph identity changed");
                        for(double fraction:{0.,.375,1.}){const auto p=a.frameAt(k,fraction),q=b.frameAt(k,fraction);
                            check(p.point.x==q.point.x && p.point.y==q.point.y && p.normal.x==q.normal.x && p.normal.y==q.normal.y && p.levelSetResidual==q.levelSetResidual,"ordinary Core frame changed");}
                        const auto p=LevelSetMappedInterface::nativeMoments(*scalar,x,rz),q=LevelSetMappedInterface::nativeMoments(*scalar,y,rz);
                        check(p.radial==q.radial && p.axial==q.axial && p.errorEstimate==q.errorEstimate,"ordinary Core geometry moment changed");}}}
            for(const auto& chain:result.curves)for(int k=0;k<static_cast<int>(chain.graphs().size());++k){const auto& g=chain.graphs()[k];++graphs;
                const auto m=LevelSetMappedInterface::nativeMoments(*scalar,g,rz);radial+=m.radial;axial+=m.axial;
                for(const auto& q:chain.quadrature(k)){length+=q.weight;measure+=q.physicalWeight;}
                if(g.nativeSide){++sideCount;const auto& s=*g.nativeSide;check(s.usesSource(*scalar),"side observation wrong source");
                    sides<<rz<<','<<amr<<','<<sign<<','<<seed.identity().value()<<','<<k<<','<<s.ownerCell()<<','<<s.normalAxis()<<','<<s.edgeLabel()<<','<<s.beginVertex()<<','<<s.endVertex()<<','<<s.minus().nativeCellIndex<<','<<s.minus().phase<<','<<s.plus().nativeCellIndex<<','<<s.plus().phase<<','<<s.begin().x<<','<<s.begin().y<<','<<s.end().x<<','<<s.end().y<<'\n';}}
            const long double expected=rz?pi:1,green=sign*(rz?pi:1)*shape.planeOffset;
            err=std::max({std::abs(length-1.L),std::abs(measure-expected)/expected,std::abs((rz?axial:radial)-green)/expected,std::abs(rz?radial:axial)/expected});
            check(curves==1 && (location==1?sideCount==graphs:sideCount==0) && err<=1e-12L,"native Core topology closed reference failed");ok=true;
        }catch(const std::exception& e){detail=clean(e.what());}
        topo<<rz<<','<<amr<<','<<location<<','<<sign<<','<<shape.planeOffset<<','<<seed.identity().value()<<','<<ok<<','<<curves<<','<<graphs<<','<<sideCount<<','<<length<<','<<measure<<','<<radial<<','<<axial<<','<<err<<','<<detail<<'\n';
        check(ok,"actual native topology source failed: "+detail);
        std::cout<<"TRACE_TOPOLOGY,"<<rz<<','<<amr<<','<<location<<','<<sign<<','<<sideCount<<','<<err<<",PASS\n";
        const auto domain=domainFor(*mesh,rz);const auto cache=InterfacePhaseGeometry::MappedFacePartitions::prepareNativeTrace(*mesh,domain,scalar);cache.validate(*mesh,domain);
        check(&cache.source()==scalar.get(),"native cache observer changed seed source");
        std::size_t halfCount=0,zeroHalves=0,valueEvents=0,normalEvents=0;long double peak=0;
        for(const auto face:domain.faces(*mesh))for(int half=0;half<2;++half){const auto& q=cache.traceHalfFace(face.nativeRow,half);++halfCount;
            const auto a=q.begin(),b=q.end();const long double dx=static_cast<long double>(b.x)-a.x,dy=static_cast<long double>(b.y)-a.y,ds=std::hypot(dx,dy);
            const auto measure=[&](long double lo,long double hi){return ds*(hi-lo)*(rz?2*pi*(a.x+dx*(hi+lo)/2):1);};
            const auto negative=[&](const std::vector<LevelSetMappedSegment::Range>& ranges){long double v=0;
                for(const auto& r:ranges)if(r.negative){const long double t=r.fromEnd?1-r.offset-r.length/2:r.offset+r.length/2;
                    v+=ds*r.length*(rz?2*pi*(a.x+dx*t):1);}return v;};
            const long double area=measure(0,1);std::size_t ve=0,ne=0;for(const auto& e:q.events())if(e.normalOrder==0)++ve;else ++ne;
            valueEvents+=ve;normalEvents+=ne;long double liquid=0,expectedLiquid=0,minus=0,plus=0,trialError=0;
            check(q.usesSource(*scalar) && q.normalAxis()==(face.horizontal?1:0),"actual cached half provenance");
            if(q.coincidentRanges().empty()){
                check(&cache.halfFace(face.nativeRow,half)==&q.ordinary(),"cache observer ordinary source fork");liquid=negative(q.ordinaryRanges());
                const long double x=rz?a.y:a.x,h=rz?dy:dx;
                if(h==0)expectedLiquid=sign*(x-shape.planeOffset)<0?area:0;
                else{const auto crossing=std::clamp((shape.planeOffset-x)/h,0.L,1.L);expectedLiquid=sign>0?measure(0,crossing):measure(crossing,1);}
                trialError=std::abs(liquid-expectedLiquid);
            }else{++zeroHalves;
                for(const auto& z:q.coincidentRanges()){check(z.minus.available && z.plus.available && z.minus.normalOrder==1 && z.plus.normalOrder==1,"actual cached coincident normal source");
                    minus+=negative(z.minus.ranges);plus+=negative(z.plus.ranges);}
                bool rejected=false;try{cache.halfFace(face.nativeRow,half);}catch(const std::exception& e){rejected=std::string(e.what()).find("no legacy phase view")!=std::string::npos;}
                check(rejected,"cached zero support leaked a legacy material phase");
                trialError=std::max(std::abs(minus-(sign>0?area:0)),std::abs(plus-(sign<0?area:0)));
            }
            check(trialError<=1e-12L*area,"actual half trace cache closed physical reference failed");if(area>0)peak=std::max(peak,trialError/area);
            cacheRows<<rz<<','<<amr<<','<<location<<','<<sign<<','<<seed.identity().value()<<','<<face.nativeRow<<','<<half<<','<<q.normalAxis()<<','<<a.x<<','<<a.y<<','<<b.x<<','<<b.y<<','<<area<<','<<q.ordinaryRanges().size()<<','<<q.coincidentRanges().size()<<','<<ve<<','<<ne<<','<<liquid<<','<<expectedLiquid<<','<<minus<<','<<plus<<','<<trialError<<'\n';
        }
        check(cache.rootIncidenceCount()==valueEvents && cache.traceEventIncidenceCount()==valueEvents+normalEvents,"actual cache event counts confused");
        cacheSummary<<rz<<','<<amr<<','<<location<<','<<sign<<','<<seed.identity().value()<<','<<halfCount<<','<<zeroHalves<<','<<valueEvents<<','<<normalEvents<<','<<peak<<'\n';
        std::cout<<"HALF_TRACE_CACHE,"<<rz<<','<<amr<<','<<location<<','<<sign<<','<<halfCount<<','<<zeroHalves<<','<<peak<<",PASS\n";
    }
}

void materialCompatibility(const std::vector<int>& grids,bool bounds=false,bool intervals=false,bool spatialSources=false)
{
    // Passive necessary-condition observer. No convection coefficients, fitted
    // flux, accepted state or reference correction is constructed here.
    std::ofstream out("output/series/material_compatibility.csv");check(out.good(),"cannot open material compatibility audit");
    std::ofstream affine("output/series/material_affine.csv");check(affine.good(),"cannot open material affine audit");
    out<<std::setprecision(17)<<"rz,amr,n,rho,dt,energy_residual,energy_gross,affine_covariance,covariance_gross,acceleration_rms_lower_bound,R,G,mass_change,boundary_power,reference_error,raw_affine_rms,raw_affine_linf,point_target_rms,pde_constant_pair_defect,pde_quadratic_pair_defect,advector_error\n";
    affine<<std::setprecision(17)<<"rz,amr,n,rho,dof,component,x,y,fixed_axis,mass,h,first_moment,raw_derivative,q,gross,advector_x,advector_y,pde_first_moment\n";
    std::ofstream boundOut,sourceOut;if(bounds){boundOut.open("output/series/material_moment_bounds.csv");sourceOut.open("output/series/material_moment_sources.csv");check(boundOut.good() && sourceOut.good(),"cannot open material moment certificate");
        boundOut<<std::setprecision(17)<<"rz,amr,n,rho,c0,c1,gross0,gross1,total_mass,physical_variance,scaled_variance,rms_lower_bound,scaled_rms_lower_bound,reference_error,mapped_c0,mapped_c1,face_c0,face_c1,first_transfer_error,second_transfer_error,source_error\n";
        sourceOut<<std::setprecision(17)<<"rz,amr,n,rho,pressure,x,y,xi,eta,primal_q\n";}
    std::ofstream intervalOut,normalOut,scaleOut;if(intervals){intervalOut.open("output/series/material_interval_summary.csv");normalOut.open("output/series/material_interval_normals.csv");scaleOut.open("output/series/material_transport_scale.csv");check(intervalOut.good() && normalOut.good() && scaleOut.good(),"cannot open interval certificate");
        intervalOut<<std::setprecision(17)<<"rz,amr,n,rho,chi,c0,c1,maximum_ratio,maximum_violation,lambda0,lambda1,defect,support,necessary_feasible,pair_identity_error,forced_scale,worst_dof,pure_gas_forced_scale,pure_gas_worst_dof\n";
        normalOut<<std::setprecision(17)<<"rz,amr,n,rho,lambda0,lambda1,c0,c1,defect,support,gross,ratio,violation\n";
        scaleOut<<std::setprecision(17)<<"rz,amr,n,rho,dof,component,x,y,fixed_axis,exterior,mass,dual_volume,liquid_volume,gas_volume,h,q,gross,advector_scale,S,forced_scale,R,begin_liquid,end_liquid,mass_change_rate,advector_x,advector_y,constant_gas_rate\n";}
    std::ofstream segmentOut,patchOut,sourceSummary,gasWitness;if(spatialSources){segmentOut.open("output/series/spatial_source_segments.csv");patchOut.open("output/series/spatial_source_patches.csv");sourceSummary.open("output/series/spatial_source_summary.csv");gasWitness.open("output/series/spatial_source_witness.csv");check(segmentOut.good() && patchOut.good() && sourceSummary.good() && gasWitness.good(),"cannot open spatial source audit");
        segmentOut<<std::setprecision(17)<<"rz,amr,n,selection,kind,id,component,normal_axis,begin_x,begin_y,end_x,end_y,area,geometry_identity,ranges,roots,liquid0,gas0,liquid1,gas1,liquid2,gas2,liquid3,gas3,liquid4,gas4,closure,quadrature_error\n";
        patchOut<<std::setprecision(17)<<"rz,amr,n,selection,face,component,support,rank,tangent_deficit,weight_norm,fixed_axis\n";
        sourceSummary<<std::setprecision(17)<<"rz,amr,n,selection,native_faces,dual_segments,unique_queries,closure,quadrature_error,rank_deficient,tangent_unidentified,maximum_tangent_deficit\n";
        gasWitness<<std::setprecision(17)<<"layer,geometry_identity,dof,side,segment,liquid,gas,closure\n";}
    for(bool rz:{false,true})for(bool amr:{false,true})for(int n:grids){
        auto mesh=makeMesh(n,amr?1:0);FaceSegmentGraph graph;graph.build(*mesh);
        const auto l=CchDofLayout::build(*mesh,NativeLeafDomain::build(*mesh,graph,rz));
        const auto metric=CchPhysicalDivergence::build(*mesh,CchIntegralDivergence::build(*mesh,l));
        const auto geometry=CchMomentumGeometry::build(*mesh,metric);
        Case shape;shape.pureFirstJet=true;shape.radial=rz;shape.motion=1;shape.speed=.4;shape.radius=rz?.34:.19;
        const double dt=.1/n; // Fine-cell CFL=0.08, independent of material density.
        const auto begin=initial(*mesh,shape);std::string stage;
        const auto mid=advanceFirstJet(*mesh,begin,shape,{1,0,dt,Position::Midpoint},stage);
        const auto end=advanceFirstJet(*mesh,begin,shape,{1,0,dt,Position::End},stage);
        const auto a=CchPhaseGeometry::build(*mesh,metric,begin),m=CchPhaseGeometry::build(*mesh,metric,mid),b=CchPhaseGeometry::build(*mesh,metric,end);
        if(spatialSources)spatialSourceAudit(segmentOut,patchOut,sourceSummary,gasWitness,*mesh,rz,amr,n,a,m,b);
        Field v(l.velocityCount()),w(l.velocityCount());std::vector<Vector2D> points;
        for(int i=0;i<l.velocityCount();++i){const auto p=l.domain().face(*mesh,l.nativeFace(i)).center;points.push_back(p);
            v[i]=component(l,i)==(rz?1:0)?.4:0;w[i]=component(l,i)==0?(rz?p.x:1+p.y):0;}
        // H10k migration begin
        const auto actual=CchAdvectorReconstruction::prepare(*mesh,geometry).evaluate(*mesh,CchNativeValues(l,v,Field(l.pressureCount())),m.geometryIdentity(),m.timeLayer());
        const auto rate=CchPhaseRate::prepareSpatial(*mesh,m,actual);
        // H10k migration end
        for(double rho:{1.,1000.,1e6}){
            const auto ma=CchMomentumMass::build(*mesh,geometry,a,rho,1),mm=CchMomentumMass::build(*mesh,geometry,m,rho,1),mb=CchMomentumMass::build(*mesh,geometry,b,rho,1);
            const auto material=CchDualMaterial::prepare(*mesh,CchMaterialTransport::prepare(*mesh,a,rate,b,rho,1),ma,mm,mb);
            // Read-only explicit research ports; no production pairing API.
            const auto& boundarySegments=geometry.segments();using Key=std::tuple<int,int,double,double>;
            std::map<Key,std::vector<int>> groups;std::vector<int> partners(material.faceRates().size(),-1);
            for(std::size_t g=0;g<boundarySegments.size();++g){const auto& s=boundarySegments[g];if(internal(s) || material.faceRates()[g]==0)continue;
                groups[{component(l,owner(s)),s.normalAxis,s.normalAxis==0?s.begin.y:s.begin.x,s.normalAxis==0?s.end.y:s.end.x}].push_back(static_cast<int>(g));}
            for(const auto& group:groups){check(group.second.size()==2,"material observer has unmatched periodic interval");const int g=group.second[0],h=group.second[1];
                const auto& left=boundarySegments[g];const auto& right=boundarySegments[h];
                check((left.negativeVelocity>=0)!=(right.negativeVelocity>=0),"material observer port orientation differs");
                check(std::abs(left.normalAxis==0?left.begin.x-right.begin.x:left.begin.y-right.begin.y)==1,"material observer port translation differs");
                check(relative(left.area-right.area,std::abs(left.area)+std::abs(right.area))<=256*std::numeric_limits<double>::epsilon(),"material observer port measure differs");
                check(relative(material.faceRates()[g]-material.faceRates()[h],std::abs(material.faceRates()[g])+std::abs(material.faceRates()[h]))<=256*std::numeric_limits<double>::epsilon(),"material observer port flux differs");
                partners[g]=h;partners[h]=g;}
            std::vector<long double> q(w.size()),first(w.size()),gross(w.size());long double power=0,powerGross=0,facePhi=0,facePsi=0;
            const auto& segments=geometry.segments();
            for(std::size_t g=0;g<segments.size();++g){const auto& s=segments[g];const auto f=material.faceRates()[g];auto c=s.begin+(s.end-s.begin)*.5;
                if(rz && s.normalAxis==1){const long double x=s.begin.x,y=s.end.x;c.x=static_cast<double>((2.L/3)*(y*y*y-x*x*x)/(y*y-x*x));}
                for(int i:{s.negativeVelocity,s.positiveVelocity})if(i>=0){const long double flux=(i==s.negativeVelocity?1:-1)*f;q[i]+=flux;gross[i]+=std::abs(f);first[i]+=flux*(rz?c.x-points[i].x:c.y-points[i].y);}
                if(internal(s)){
                    const int i=s.negativeVelocity,j=s.positiveVelocity;facePhi+=.5L*f*(static_cast<long double>(w[i])*w[i]-static_cast<long double>(w[j])*w[j]);
                    if(component(l,i)==0){const long double zi=rz?points[i].x:points[i].y,zj=rz?points[j].x:points[j].y;
                        facePsi+=f*(zi-zj)*((rz?c.x:c.y)-.5L*(zi+zj));}
                }
                if(!internal(s)){const int i=owner(s);const long double flux=(s.negativeVelocity>=0?1:-1)*f;
                    const long double trace=partners[g]>=0?.5L*(w[i]+w[owner(segments[partners[g]])]):w[i];
                    const auto term=flux*(w[i]*trace-.5L*w[i]*w[i]);power+=term;powerGross+=std::abs(term);facePhi+=.5L*flux*w[i]*w[i]-term;
                    if(component(l,i)==0){const long double z=rz?points[i].x:points[i].y;facePsi+=flux*(z*(rz?c.x:c.y)-.5L*z*z);}}
            }
            long double phi=-power,phiGross=powerGross,psi=0,psiGross=0,mass=0,norm=0,maxR=0,maxG=0,maxChange=0;
            for(std::size_t i=0;i<w.size();++i){const auto term=.5L*q[i]*w[i]*w[i];phi+=term;phiGross+=std::abs(term);
                if(component(l,static_cast<int>(i))==0){const long double z=rz?points[i].x:points[i].y;
                    const long double covariance=z*first[i]+.5L*q[i]*z*z;psi+=covariance;psiGross+=std::abs(z*first[i])+std::abs(.5L*q[i]*z*z);}
                mass+=mm.masses()[i];norm+=mm.masses()[i]*w[i]*w[i];const auto& cell=material.cells()[i];
                maxR=std::max(maxR,std::abs(cell.defectRate));maxG=std::max(maxG,std::abs(cell.geometryTransferChange));maxChange=std::max(maxChange,std::abs(cell.massChange));}
            const auto lower=std::abs(phi)/std::sqrt(mass*norm);
            check(std::isfinite(lower) && norm>0,"invalid material compatibility norm");
            const auto reference=std::max(relative(phi-facePhi,phiGross),relative(psi-facePsi,psiGross));
            check(reference<=1e-12L,"independent material moment expression differs from native sums");
            if(rho==1)check(relative(phi,phiGross)<=1e-12L,"homogeneous affine material-KE control failed");
            const auto advector=materialAdvector(*mesh,material);long double squared=0,linf=0,pointSquared=0,advectorError=0,c0=0,c1=0;
            for(std::size_t i=0;i<w.size();++i){const bool used=component(l,static_cast<int>(i))==0,fixed=rz && used && points[i].x==0;const long double z=rz?points[i].x:points[i].y;
                const auto& box=geometry.volumes()[i];const double h=std::max(box.upper.x-box.lower.x,box.upper.y-box.lower.y);
                const long double raw=used?-first[i]/mm.masses()[i]:0,pdeTarget=used?mm.masses()[i]*(rz?advector[i].x:advector[i].y):0;
                if(!fixed){squared+=mm.masses()[i]*raw*raw;linf=std::max(linf,std::abs(raw));pointSquared+=pdeTarget*pdeTarget/mm.masses()[i];}
                if(used){c0-=pdeTarget+z*q[i];c1-=z*pdeTarget+.5L*q[i]*z*z;}
                advectorError=std::max(advectorError,static_cast<long double>(std::hypot(advector[i].x-(rz?0:.4),advector[i].y-(rz?.4:0))));
                affine<<rz<<','<<amr<<','<<n<<','<<rho<<','<<i<<','<<component(l,static_cast<int>(i))<<','<<points[i].x<<','<<points[i].y<<','<<fixed<<','<<mm.masses()[i]<<','<<h<<','<<first[i]<<','<<raw<<','<<q[i]<<','<<gross[i]<<','<<advector[i].x<<','<<advector[i].y<<','<<pdeTarget<<'\n';
            }
            check(advectorError<=1e-12L,"actual native advector reconstruction failed constant reference");
            if(bounds)materialMomentBounds(boundOut,sourceOut,*mesh,rz,amr,n,rho,material,points,q,advector);
            if(intervals)materialIntervalCertificate(intervalOut,normalOut,scaleOut,rz,amr,n,rho,material,points,q,advector,partners);
            out<<rz<<','<<amr<<','<<n<<','<<rho<<','<<dt<<','<<phi<<','<<phiGross<<','<<psi<<','<<psiGross<<','<<lower<<','<<maxR<<','<<maxG<<','<<maxChange<<','<<power<<','<<reference<<','<<std::sqrt(squared/mass)<<','<<linf<<','<<std::sqrt(pointSquared/mass)<<','<<c0<<','<<c1<<','<<advectorError<<'\n';
            std::cout<<"MATERIAL_COMPATIBILITY,"<<rz<<','<<amr<<','<<n<<','<<rho<<','<<phi<<','<<psi<<','<<lower<<','<<maxR<<','<<maxG<<std::endl;
            std::cout<<"MATERIAL_AFFINE,"<<rz<<','<<amr<<','<<n<<','<<rho<<','<<std::sqrt(squared/mass)<<','<<linf<<','<<std::sqrt(pointSquared/mass)<<','<<c0<<','<<c1<<','<<advectorError<<std::endl;
        }
    }
}
void centerCertificate(const std::vector<int>& grids)
{
    // A feasibility witness, not a minimum-change solver or a trajectory.
    std::ofstream rows("output/series/center_certificate_rows.csv"),coefficients("output/series/center_certificate_pairs.csv");
    check(rows.good() && coefficients.good(),"cannot open center coefficient certificate");
    rows<<std::setprecision(17)<<"rz,amr,n,dof,x,y,mass,h,gross,q_half,first_x,first_y,value_constant,value_x,value_y,value_xx,value_xy,value_yy,quadratic_bound,l1,constant_error,affine_error,quadratic_error,l1_ratio,pde_x,pde_y\n";
    coefficients<<std::setprecision(17)<<"rz,amr,n,first,second,coefficient\n";
    const auto include=[](std::vector<int>& list,int i){if(std::find(list.begin(),list.end(),i)==list.end())list.push_back(i);};
    for(bool rz:{false,true})for(bool amr:{false,true})for(int n:grids){auto mesh=makeMesh(n,amr?1:0);const auto material=prepare(*mesh,rz);const auto& l=material.layout();const auto& segments=material.geometry().segments();const int count=l.velocityCount();
        std::vector<Vector2D> points;std::vector<int> representative(count),partners(segments.size(),-1);std::vector<std::vector<int>> cellFaces(l.pressureCount()),neighbors(l.pressureCount()),roots(count);
        for(int i=0;i<count;++i){representative[i]=i;const auto face=l.domain().face(*mesh,l.nativeFace(i));points.push_back(face.center);
            const int a=face.negativeCell<0?-1:l.pressureDof(face.negativeCell),b=face.positiveCell<0?-1:l.pressureDof(face.positiveCell);
            if(a>=0)cellFaces[a].push_back(i);if(b>=0)cellFaces[b].push_back(i);if(a>=0 && b>=0){include(neighbors[a],b);include(neighbors[b],a);}}
        const auto root=[&](int i){while(representative[i]!=i)i=representative[i];return i;};
        Vector2D period{};using Key=std::tuple<int,int,double,double>;std::map<Key,std::vector<int>> ports;
        for(std::size_t g=0;g<segments.size();++g){const auto& s=segments[g];if(internal(s) || material.faceRates()[g]==0)continue;
            ports[{component(l,owner(s)),s.normalAxis,s.normalAxis==0?s.begin.y:s.begin.x,s.normalAxis==0?s.end.y:s.end.x}].push_back(static_cast<int>(g));}
        for(const auto& group:ports){check(group.second.size()==2,"certificate has unpaired interval");const int g=group.second[0],h=group.second[1];const auto& a=segments[g];const auto& b=segments[h];const int i=owner(a),j=owner(b);
            check((a.negativeVelocity>=0)!=(b.negativeVelocity>=0),"certificate paired orientations differ");
            check(relative(a.area-b.area,a.area+b.area)<=256*std::numeric_limits<double>::epsilon(),"certificate paired measures differ");
            check(relative(material.faceRates()[g]-material.faceRates()[h],std::abs(material.faceRates()[g])+std::abs(material.faceRates()[h]))<=256*std::numeric_limits<double>::epsilon(),"certificate paired rates differ");
            const double length=std::abs(a.normalAxis==0?a.begin.x-b.begin.x:a.begin.y-b.begin.y);check(length>0,"zero certificate period");if(a.normalAxis==0)period.x=length;else period.y=length;partners[g]=h;partners[h]=g;
            if(component(l,i)==a.normalAxis){const int ri=root(i),rj=root(j);representative[std::max(ri,rj)]=std::min(ri,rj);
                for(const auto& x:material.geometry().volumes()[i].overlaps)for(const auto& y:material.geometry().volumes()[j].overlaps){include(neighbors[x.pressure],y.pressure);include(neighbors[y.pressure],x.pressure);}}}
        const auto delta=[&](Vector2D a,Vector2D b){auto d=b-a;if(period.x>0)d.x-=period.x*std::round(d.x/period.x);if(period.y>0)d.y-=period.y*std::round(d.y/period.y);return d;};
        const auto actualVelocity=materialAdvector(*mesh,material);
        std::vector<double> widths(count);std::vector<long double> mass(count),gross(count),norm(count),pdeX(count),pdeY(count);std::vector<std::array<long double,6>> target(count),value(count),termGross(count);
        Vector2D lower=points.front(),upper=lower;
        for(int i=0;i<count;++i){const int r=root(i);representative[i]=r;for(const auto& p:material.geometry().volumes()[i].overlaps)include(roots[r],p.pressure);
            const auto& box=material.geometry().volumes()[i];widths[r]=std::max(widths[r],std::max(box.upper.x-box.lower.x,box.upper.y-box.lower.y));mass[r]+=material.midpointMass().masses()[i];
            pdeX[r]+=material.midpointMass().masses()[i]*actualVelocity[i].x;pdeY[r]+=material.midpointMass().masses()[i]*actualVelocity[i].y;
            lower.x=std::min(lower.x,points[i].x);lower.y=std::min(lower.y,points[i].y);upper.x=std::max(upper.x,points[i].x);upper.y=std::max(upper.y,points[i].y);}
        const double length=std::max({upper.x-lower.x,upper.y-lower.y,period.x,period.y});
        using Pair=std::pair<int,int>;std::map<Pair,long double> center;
        const auto add=[&](int i,int j,long double k){i=representative[i];j=representative[j];if(i==j)return;if(i>j){std::swap(i,j);k=-k;}center[{i,j}]+=k;};
        for(std::size_t g=0;g<segments.size();++g){const auto& s=segments[g];const auto f=material.faceRates()[g];auto p=s.begin+(s.end-s.begin)*.5;
            if(rz && s.normalAxis==1){const long double a=s.begin.x,b=s.end.x;p.x=static_cast<double>((2.L/3)*(b*b*b-a*a*a)/(b*b-a*a));}
            for(int i:{s.negativeVelocity,s.positiveVelocity})if(i>=0){const int r=representative[i];const long double signedRate=(i==s.negativeVelocity?1:-1)*f;gross[r]+=std::abs(f);target[r][0]+=.5L*signedRate;
                const auto d=delta(points[r],p);target[r][1]+=signedRate*d.x/widths[r];target[r][2]+=signedRate*d.y/widths[r];}
            if(internal(s))add(s.negativeVelocity,s.positiveVelocity,.5L*f);else if(static_cast<int>(g)<partners[g])add(owner(s),owner(segments[partners[g]]),.5L*(s.negativeVelocity>=0?1:-1)*f);
        }
        std::set<Pair> allowed;
        for(int i=0;i<count;++i)if(representative[i]==i){auto patch=roots[i];std::vector<int> support;for(int p:roots[i])for(int j:neighbors[p])include(patch,j);
            for(int p:patch)for(int j:cellFaces[p])if(component(l,i)==component(l,j))include(support,representative[j]);std::sort(support.begin(),support.end());
            for(std::size_t a=0;a<support.size();++a)for(std::size_t b=a+1;b<support.size();++b)allowed.emplace(support[a],support[b]);}
        for(const auto& entry:center)check(entry.second==0 || allowed.count(entry.first)==1,"center witness escapes true closed one-ring");
        for(const auto& pair:allowed){const auto [i,j]=pair;const auto found=center.find(pair);const long double k=found==center.end()?0:found->second;coefficients<<rz<<','<<amr<<','<<n<<','<<i<<','<<j<<','<<k<<'\n';
            for(int r:{i,j}){const auto d=delta(points[r],points[r==i?j:i]);const long double x=d.x/widths[r],y=d.y/widths[r],signedK=r==i?k:-k;const std::array<long double,6> p{1,x,y,x*x,x*y,y*y};
                norm[r]+=std::abs(k);for(int c=0;c<6;++c){value[r][c]+=signedK*p[c];termGross[r][c]+=std::abs(signedK*p[c]);}}}
        long double maximumConstant=0,maximumAffine=0,maximumQuadratic=0,maximumNorm=0;int activeRows=0;
        for(int i=0;i<count;++i)if(representative[i]==i){++activeRows;long double ec=0,ea=0,eq=0;const auto bound=gross[i]*widths[i]/length;
            for(int c=0;c<3;++c){const auto scale=gross[i]+std::abs(target[i][c])+termGross[i][c];const auto e=scale>0?std::abs(value[i][c]-target[i][c])/scale:std::abs(value[i][c]-target[i][c]);if(c==0)ec=std::max(ec,e);else ea=std::max(ea,e);}
            for(int c=3;c<6;++c){const auto scale=bound+gross[i]+termGross[i][c],excess=std::max(0.L,std::abs(value[i][c])-bound);eq=std::max(eq,scale>0?excess/scale:excess);}
            const auto ratio=gross[i]>0?norm[i]/(2*gross[i]):(norm[i]==0?0:std::numeric_limits<long double>::infinity());
            maximumConstant=std::max(maximumConstant,ec);maximumAffine=std::max(maximumAffine,ea);maximumQuadratic=std::max(maximumQuadratic,eq);maximumNorm=std::max(maximumNorm,ratio);
            rows<<rz<<','<<amr<<','<<n<<','<<i<<','<<points[i].x<<','<<points[i].y<<','<<mass[i]<<','<<widths[i]<<','<<gross[i];for(int c=0;c<3;++c)rows<<','<<target[i][c];for(int c=0;c<6;++c)rows<<','<<value[i][c];rows<<','<<bound<<','<<norm[i]<<','<<ec<<','<<ea<<','<<eq<<','<<ratio<<','<<pdeX[i]/widths[i]<<','<<pdeY[i]/widths[i]<<'\n';}
        const long double tolerance=256*std::numeric_limits<double>::epsilon();const bool feasible=std::max({maximumConstant,maximumAffine,maximumQuadratic})<=tolerance && maximumNorm<=1+tolerance;
        std::cout<<"CENTER_CERTIFICATE,"<<rz<<','<<amr<<','<<n<<','<<activeRows<<','<<center.size()<<','<<feasible<<','<<maximumConstant<<','<<maximumAffine<<','<<maximumQuadratic<<','<<maximumNorm<<std::endl;
        if(!rz && !amr && n==16)check(feasible,"uniform n=16 center witness did not satisfy fixed H10e bounds");
    }
}
}
// H10j0 observer begin
// Passive actual owner study; analytic data are inputs/reference only here.
// No phase/material rate, convection or accepted state is changed.
namespace {
void spatialAdvectorAudit()
{
    std::filesystem::create_directories("output/series");
    std::ofstream faces("output/series/advector_faces.csv"),atoms("output/series/advector_atoms.csv"),cells("output/series/advector_cells.csv"),summary("output/series/advector_summary.csv");
    check(faces.good() && atoms.good() && cells.good() && summary.good(),"cannot open actual advector observations");
    faces<<std::setprecision(17)<<"rz,amr,n,selection,revision,geometry,advector,step,clock_begin,clock_end,position,u,axis,x,y,input,point_x,point_y,integral,area\n";
    atoms<<std::setprecision(17)<<"rz,amr,n,selection,revision,atom,axis,primal,owner,ax,ay,bx,by,area,length";
    writeCoefficientHeader<decltype(CchNormalPolynomial::beginCoefficients)>(atoms,"c");
    writeCoefficientHeader<decltype(CchNormalPolynomial::endCoefficients)>(atoms,"e");atoms<<",integral,end_integral\n";
    cells<<std::setprecision(17)<<"rz,amr,n,selection,revision,p,lo_x,lo_y,hi_x,hi_y,x,y,velocity_x,velocity_y,divergence,original_divergence,volume\n";
    summary<<std::setprecision(17)<<"rz,amr,n,selection,faces,segments,atoms,weights,cell_pieces,view_entries,svd,rank_deficient,unidentified,geometry_shared,updates_shared,fresh_identity,prepare_seconds,evaluate_seconds\n";
    for(bool rz:{false,true})for(bool amr:{false,true})for(int n:{4,8,16,32}){
        auto mesh=makeMesh(n,amr?1:0);const auto pub=publication(*mesh,rz,{9,.2,.3,Position::Midpoint});
        for(int selection=0;selection<4;++selection){FaceSegmentGraph graph;graph.build(*mesh,[selection](const QuadtreeMesh&,const Cell& c){const auto p=c.centerCoord();
            if(selection==1)return p.x<.5;if(selection==2)return !(p.x>.25 && p.x<.75 && p.y>.25 && p.y<.75);if(selection==3)return p.x<.25 || p.x>.75;return true;});
            const auto l=CchDofLayout::build(*mesh,NativeLeafDomain::build(*mesh,graph,rz));const auto metric=CchPhysicalDivergence::build(*mesh,CchIntegralDivergence::build(*mesh,l));const auto geometry=CchMomentumGeometry::build(*mesh,metric);
            const auto start=std::chrono::steady_clock::now();const auto reconstruction=CchAdvectorReconstruction::prepare(*mesh,geometry);const auto prepared=std::chrono::steady_clock::now();double updateSeconds=0;bool same=true;std::uint64_t previous=0;
            for(int revision=0;revision<2;++revision){Field velocity(l.velocityCount());
                for(int u=0;u<l.velocityCount();++u){const auto f=l.domain().face(*mesh,l.nativeFace(u));const int axis=component(l,u);const auto& edge=f.horizontal?mesh->hEdges(f.edgeLabel):mesh->vEdges(f.edgeLabel);const auto a=mesh->verts(edge.edgeToVert[0]).coord(),b=mesh->verts(edge.edgeToVert[1]).coord();
                    if(rz)velocity[u]=axis==0?(.13+.05*revision)*f.center.x:.2+.11*revision+.31*.25*(a.x*a.x+b.x*b.x)+.19*f.center.y;
                    else velocity[u]=axis==0?.3+.11*revision+.17*f.center.x-.23*f.center.y:-.2-.07*revision+.31*f.center.x+.19*f.center.y;}
                const CchNativeValues values(l,velocity,Field(l.pressureCount()));const auto update=std::chrono::steady_clock::now();const auto advector=reconstruction.evaluate(*mesh,values,pub.state().identity(),pub.state().timeLayer());updateSeconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-update).count();
                advector.validateFor(*mesh,geometry,pub.state().identity(),pub.state().timeLayer());same=same && advector.reconstruction().sameIdentity(reconstruction);check(previous!=advector.identity().value(),"advector revision reused identity");previous=advector.identity().value();
                for(int u=0;u<l.velocityCount();++u){const auto f=l.domain().face(*mesh,l.nativeFace(u));const auto v=advector.nativeVelocity(u);long double total=0;for(int atom:reconstruction.primalAtoms()[u]){const auto& poly=advector.atomPolynomials()[atom];total+=poly.integral(false,0,poly.length);}const auto clock=advector.timeLayer();
                    faces<<rz<<','<<amr<<','<<n<<','<<selection<<','<<revision<<','<<advector.geometryIdentity().value()<<','<<advector.identity().value()<<','<<clock.step<<','<<clock.begin<<','<<clock.end<<','<<static_cast<int>(clock.position)<<','<<u<<','<<component(l,u)<<','<<f.center.x<<','<<f.center.y<<','<<velocity[u]<<','<<v.x<<','<<v.y<<','<<total<<','<<metric.faceMeasures()[u]<<'\n';}
                for(std::size_t i=0;i<reconstruction.atoms().size();++i){const auto& a=reconstruction.atoms()[i];const auto& poly=advector.atomPolynomials()[i];atoms<<rz<<','<<amr<<','<<n<<','<<selection<<','<<revision<<','<<i<<','<<a.normalAxis<<','<<a.primalVelocity<<','<<a.ownerPressure<<','<<a.begin.x<<','<<a.begin.y<<','<<a.end.x<<','<<a.end.y<<','<<a.area<<','<<poly.length;writeCoefficients(atoms,poly.beginCoefficients);writeCoefficients(atoms,poly.endCoefficients);atoms<<','<<poly.integral(false,0,poly.length)<<','<<poly.integral(true,0,poly.length)<<'\n';}
                const auto div=metric.apply(*mesh,values);for(int p=0;p<l.pressureCount();++p){const auto& cell=mesh->cells(mesh->cellIndexToLabel[l.nativeCell(p)]);const auto lo=mesh->verts(cell.cellToVert[0]).coord(),hi=mesh->verts(cell.cellToVert[2]).coord();const auto point=lo+Vector2D{.37*(hi.x-lo.x),.61*(hi.y-lo.y)};const auto v=advector.cellVelocity(p,point);
                    cells<<rz<<','<<amr<<','<<n<<','<<selection<<','<<revision<<','<<p<<','<<lo.x<<','<<lo.y<<','<<hi.x<<','<<hi.y<<','<<point.x<<','<<point.y<<','<<v.x<<','<<v.y<<','<<advector.cellDivergence()[p]<<','<<div[p]<<','<<metric.cellVolumes()[p]<<'\n';}
            }
            const auto& r=reconstruction.report();summary<<rz<<','<<amr<<','<<n<<','<<selection<<','<<l.velocityCount()<<','<<geometry.segments().size()<<','<<r.atoms<<','<<r.weights<<','<<r.cellPieces<<','<<r.viewEntries<<','<<r.svdCount<<','<<r.rankDeficient<<','<<r.tangentUnidentified<<','<<reconstruction.momentumGeometry().sameIdentity(geometry)<<','<<same<<",1,"<<std::chrono::duration<double>(prepared-start).count()<<','<<updateSeconds<<'\n';
        }
    }
    std::cout<<"H10j0 actual advector audit: 64 layouts / 128 revisions PASS\n";
}
}
// H10j0 observer end
// H10j observer begin
namespace {
void spatialPhaseAudit()
{
    std::filesystem::create_directories("output/series");
    std::ofstream rows("output/series/spatial_phase_atoms.csv"),views("output/series/spatial_phase_views.csv"),summary("output/series/spatial_phase_summary.csv"),gas("output/series/spatial_gas_witness.csv");
    check(rows.good() && views.good() && summary.good() && gas.good(),"cannot open spatial phase observation");
    rows<<std::setprecision(17)<<"rz,amr,n,selection,shape,geometry,advector,step,begin,end,position,atom,axis,ax,ay,bx,by,area,liquid,gas,total,cx,cy,radius,plane,speed,tangent_unidentified\n";
    views<<std::setprecision(17)<<"rz,amr,n,selection,shape,kind,id,axis,negative,positive,ax,ay,bx,by,area,liquid,gas,atoms,total\n";
    summary<<std::setprecision(17)<<"rz,amr,n,selection,shape,faces,segments,atoms,queries,ordinary,coincident,events,refinements,integration_error,refinement_change,closure,tangent_unidentified,same_advector,same_scalar,seconds\n";
    gas<<std::setprecision(17)<<"layer,geometry,advector,dof,side,segment,axis,ax,ay,bx,by,area,liquid,gas\n";
    for(bool rz:{false,true})for(bool amr:{false,true})for(int n:{4,8,16,32}){
        auto mesh=makeMesh(n,amr?1:0);Case circle;circle.pureFirstJet=true;circle.radial=rz;circle.motion=1;circle.speed=.4;circle.radius=rz?.34:.19;const double dt=.1/n;std::string stage;
        const auto initialCircle=initial(*mesh,circle);const auto midpoint=advanceFirstJet(*mesh,initialCircle,circle,{1,0,dt,Position::Midpoint},stage);const auto ending=advanceFirstJet(*mesh,initialCircle,circle,{1,0,dt,Position::End},stage);
        Case plane;plane.pureFirstJet=true;plane.radial=rz;plane.planar=true;plane.planeOffset=.43;const auto planar=initial(*mesh,plane);
        for(int selection=0;selection<4;++selection){FaceSegmentGraph graph;graph.build(*mesh,[selection](const QuadtreeMesh&,const Cell& c){const auto p=c.centerCoord();if(selection==1)return p.x<.5;if(selection==2)return !(p.x>.25 && p.x<.75 && p.y>.25 && p.y<.75);if(selection==3)return p.x<.25 || p.x>.75;return true;});
            const auto l=CchDofLayout::build(*mesh,NativeLeafDomain::build(*mesh,graph,rz));const auto metrics=CchPhysicalDivergence::build(*mesh,CchIntegralDivergence::build(*mesh,l));const auto geometry=CchMomentumGeometry::build(*mesh,metrics);const auto reconstruction=CchAdvectorReconstruction::prepare(*mesh,geometry);
            Field velocity(l.velocityCount());for(int u=0;u<l.velocityCount();++u){const auto f=l.domain().face(*mesh,l.nativeFace(u));if(!rz)velocity[u]=component(l,u)==0?.3+.17*f.center.x-.23*f.center.y:-.2+.31*f.center.x+.19*f.center.y;
                else if(component(l,u)==0)velocity[u]=.13*f.center.x;else{const auto& e=mesh->hEdges(f.edgeLabel);const auto a=mesh->verts(e.edgeToVert[0]).coord(),b=mesh->verts(e.edgeToVert[1]).coord();velocity[u]=.2+.155*.5*(a.x*a.x+b.x*b.x)+.19*f.center.y;}}
            for(int shape=0;shape<2;++shape){const auto& pub=shape==0?planar:midpoint;const auto phase=CchPhaseGeometry::build(*mesh,metrics,pub);const auto advector=reconstruction.evaluate(*mesh,CchNativeValues(l,velocity,Field(l.pressureCount())),phase.geometryIdentity(),phase.timeLayer());const auto started=std::chrono::steady_clock::now();const auto rate=CchPhaseRate::prepareSpatial(*mesh,phase,advector);const auto seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
                const auto clock=phase.timeLayer();for(std::size_t i=0;i<reconstruction.atoms().size();++i){const auto& a=reconstruction.atoms()[i];const auto value=rate.atomRates()[i];
                    rows<<rz<<','<<amr<<','<<n<<','<<selection<<','<<shape<<','<<phase.geometryIdentity().value()<<','<<rate.advectorIdentity().value()<<','<<clock.step<<','<<clock.begin<<','<<clock.end<<','<<static_cast<int>(clock.position)<<','<<i<<','<<a.normalAxis<<','<<a.begin.x<<','<<a.begin.y<<','<<a.end.x<<','<<a.end.y<<','<<a.area<<','<<value.liquid<<','<<value.gas<<','<<advector.atomPolynomials()[i].integral(false,0,advector.atomPolynomials()[i].length)<<','<<(rz?0:circle.cx+circle.speed*clock.time())<<','<<circle.cy+(rz?circle.speed*clock.time():0)<<','<<circle.radius<<','<<plane.planeOffset<<','<<circle.speed<<','<<reconstruction.report().tangentUnidentified<<'\n';}
                for(int u=0;u<l.velocityCount();++u){const auto f=l.domain().face(*mesh,l.nativeFace(u));const auto& e=f.horizontal?mesh->hEdges(f.edgeLabel):mesh->vEdges(f.edgeLabel);auto a=mesh->verts(e.edgeToVert[0]).coord(),b=mesh->verts(e.edgeToVert[1]).coord();if(component(l,u)==0?a.y>b.y:a.x>b.x)std::swap(a,b);const auto v=rate.faceRates()[u];
                    views<<rz<<','<<amr<<','<<n<<','<<selection<<','<<shape<<",primal,"<<u<<','<<component(l,u)<<','<<f.negativeCell<<','<<f.positiveCell<<','<<a.x<<','<<a.y<<','<<b.x<<','<<b.y<<','<<metrics.faceMeasures()[u]<<','<<v.liquid<<','<<v.gas<<','<<reconstruction.primalAtoms()[u].size()<<','<<rate.totalFaceRates()[u]<<'\n';}
                for(std::size_t i=0;i<geometry.segments().size();++i){const auto& a=geometry.segments()[i];const auto v=rate.dualFaceRates()[i];long double total=0;for(int atom:reconstruction.dualAtoms()[i])total+=advector.atomPolynomials()[atom].integral(false,0,advector.atomPolynomials()[atom].length);
                    views<<rz<<','<<amr<<','<<n<<','<<selection<<','<<shape<<",dual,"<<i<<','<<a.normalAxis<<','<<a.negativeVelocity<<','<<a.positiveVelocity<<','<<a.begin.x<<','<<a.begin.y<<','<<a.end.x<<','<<a.end.y<<','<<a.area<<','<<v.liquid<<','<<v.gas<<','<<reconstruction.dualAtoms()[i].size()<<','<<total<<'\n';}
                const auto& r=rate.spatialReport();summary<<rz<<','<<amr<<','<<n<<','<<selection<<','<<shape<<','<<l.velocityCount()<<','<<geometry.segments().size()<<','<<r.atoms<<','<<r.sourceQueries<<','<<r.ordinaryRanges<<','<<r.coincidentRanges<<','<<r.events<<','<<r.refinements<<','<<r.positiveIntegrationErrorEstimate<<','<<r.maxRelativeRefinementChange<<','<<r.maxRelativeAtomClosure<<','<<r.tangentUnidentified<<','<<rate.spatialAdvector().reconstruction().sameIdentity(reconstruction)<<','<<(&rate.geometry().publication().state().geometry().snapshot()==&pub.state().geometry().snapshot())<<','<<seconds<<'\n';
            }
            if(!rz && amr && n==16 && selection==0){Field uniform(l.velocityCount());for(int u=0;u<l.velocityCount();++u)uniform[u]=component(l,u)==0?.4:0;
                for(int layer=0;layer<3;++layer){const auto& pub=layer==0?initialCircle:layer==1?midpoint:ending;const auto phase=CchPhaseGeometry::build(*mesh,metrics,pub);const auto advector=reconstruction.evaluate(*mesh,CchNativeValues(l,uniform,Field(l.pressureCount())),phase.geometryIdentity(),phase.timeLayer());const auto rate=CchPhaseRate::prepareSpatial(*mesh,phase,advector);
                    for(int side=0;side<4;++side)for(int i:geometry.volumes()[103].sides[side]){const auto& a=geometry.segments()[i];const auto v=rate.dualFaceRates()[i];check(v.liquid==0,"new spatial pure-gas witness leaked liquid");gas<<layer<<','<<phase.geometryIdentity().value()<<','<<advector.identity().value()<<",103,"<<side<<','<<i<<','<<a.normalAxis<<','<<a.begin.x<<','<<a.begin.y<<','<<a.end.x<<','<<a.end.y<<','<<a.area<<','<<v.liquid<<','<<v.gas<<'\n';}}
            }
        }
        std::cout<<"SPATIAL_PHASE_LAYOUT,"<<rz<<','<<amr<<','<<n<<",PASS\n";
    }
}
}
// H10j observer end
// H10k observer begin
namespace {
void spatialMaterialAudit()
{
    std::filesystem::create_directories("output/series");
    std::ofstream atoms("output/series/material_spatial_atoms.csv"),views("output/series/material_spatial_views.csv"),cells("output/series/material_spatial_cells.csv"),summary("output/series/material_spatial_summary.csv"),gas("output/series/material_gas_witness.csv");
    check(atoms.good() && views.good() && cells.good() && summary.good() && gas.good(),"cannot open spatial material observation");
    atoms<<std::setprecision(17)<<"rz,amr,n,selection,shape,rho,atom,axis,ax,ay,bx,by,area,liquid,gas,F,total\n";
    views<<std::setprecision(17)<<"rz,amr,n,selection,shape,rho,kind,id,axis,negative,positive,ax,ay,bx,by,area,F,atoms\n";
    cells<<std::setprecision(17)<<"rz,amr,n,selection,shape,rho,p,native,lo_x,lo_y,hi_x,hi_y,volume,liquid_a,gas_a,liquid_m,gas_m,liquid_b,gas_b,mass_a,mass_m,mass_b,change,outward,defect,R,relative\n";
    summary<<std::setprecision(17)<<"rz,amr,n,selection,shape,rho,a,m,b,advector,step,begin,end,position,faces,segments,atoms,phase_queries,unidentified,mass_a,mass_m,mass_b,change,boundary,defect,l1defect,balance,same_F,same_phase,same_geometry,seconds,cxa,cya,cxm,cym,cxb,cyb,radius,plane\n";
    gas<<std::setprecision(17)<<"layer,rho,geometry,advector,step,begin,end,position,dof,side,segment,axis,ax,ay,bx,by,area,liquid,gas,F\n";
    for(bool rz:{false,true})for(bool amr:{false,true})for(int n:{4,8,16,32}){
        auto mesh=makeMesh(n,amr?1:0);Case circle;circle.pureFirstJet=true;circle.radial=rz;circle.motion=1;circle.speed=.4;circle.radius=rz?.34:.19;const double dt=.1/n;std::string stage;
        const auto ca=initial(*mesh,circle),cm=advanceFirstJet(*mesh,ca,circle,{1,0,dt,Position::Midpoint},stage),cb=advanceFirstJet(*mesh,ca,circle,{1,0,dt,Position::End},stage);
        Case plane;plane.pureFirstJet=true;plane.radial=rz;plane.planar=true;plane.planeOffset=.43;
        const auto pa=initial(*mesh,plane),pm=advanceFirstJet(*mesh,pa,plane,{1,0,dt,Position::Midpoint},stage),pb=advanceFirstJet(*mesh,pa,plane,{1,0,dt,Position::End},stage);
        for(int selection=0;selection<4;++selection){FaceSegmentGraph graph;graph.build(*mesh,[selection](const QuadtreeMesh&,const Cell& c){const auto p=c.centerCoord();if(selection==1)return p.x<.5;if(selection==2)return !(p.x>.25 && p.x<.75 && p.y>.25 && p.y<.75);if(selection==3)return p.x<.25 || p.x>.75;return true;});
            const auto l=CchDofLayout::build(*mesh,NativeLeafDomain::build(*mesh,graph,rz));const auto metrics=CchPhysicalDivergence::build(*mesh,CchIntegralDivergence::build(*mesh,l));const auto geometry=CchMomentumGeometry::build(*mesh,metrics);const auto reconstruction=CchAdvectorReconstruction::prepare(*mesh,geometry);
            Field velocity(l.velocityCount());for(int u=0;u<l.velocityCount();++u){const auto f=l.domain().face(*mesh,l.nativeFace(u));if(!rz)velocity[u]=component(l,u)==0?.3+.17*f.center.x-.23*f.center.y:-.2+.31*f.center.x+.19*f.center.y;
                else if(component(l,u)==0)velocity[u]=.13*f.center.x;else{const auto& e=mesh->hEdges(f.edgeLabel);const auto a=mesh->verts(e.edgeToVert[0]).coord(),b=mesh->verts(e.edgeToVert[1]).coord();velocity[u]=.2+.155*.5*(a.x*a.x+b.x*b.x)+.19*f.center.y;}}
            for(int shape=0;shape<2;++shape){const auto a=CchPhaseGeometry::build(*mesh,metrics,shape==0?pa:ca),m=CchPhaseGeometry::build(*mesh,metrics,shape==0?pm:cm),b=CchPhaseGeometry::build(*mesh,metrics,shape==0?pb:cb);
                const auto advector=reconstruction.evaluate(*mesh,CchNativeValues(l,velocity,Field(l.pressureCount())),m.geometryIdentity(),m.timeLayer());const auto phase=CchPhaseRate::prepareSpatial(*mesh,m,advector);
                for(double rho:{1.,1000.,1e6}){const auto started=std::chrono::steady_clock::now();const auto flux=CchMaterialFlux::prepare(*mesh,phase,rho,1);const auto material=CchMaterialTransport::prepare(*mesh,a,flux,b);const auto seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
                    material.validateFor(*mesh,geometry,material.sources(),material.timeLayer());const auto copied=material;
                    const bool sameF=&flux.atomRates()==&material.atomRates() && &flux.faceRates()==&material.faceRates() && &flux.dualFaceRates()==&material.dualFaceRates() && &copied.atomRates()==&material.atomRates();check(sameF,"material study interval recomputed unique F");
                    for(std::size_t i=0;i<reconstruction.atoms().size();++i){const auto& atom=reconstruction.atoms()[i];const auto q=phase.atomRates()[i];atoms<<rz<<','<<amr<<','<<n<<','<<selection<<','<<shape<<','<<rho<<','<<i<<','<<atom.normalAxis<<','<<atom.begin.x<<','<<atom.begin.y<<','<<atom.end.x<<','<<atom.end.y<<','<<atom.area<<','<<q.liquid<<','<<q.gas<<','<<flux.atomRates()[i]<<','<<advector.atomPolynomials()[i].integral(false,0,advector.atomPolynomials()[i].length)<<'\n';}
                    for(int u=0;u<l.velocityCount();++u){const auto f=l.domain().face(*mesh,l.nativeFace(u));const auto& e=f.horizontal?mesh->hEdges(f.edgeLabel):mesh->vEdges(f.edgeLabel);auto lo=mesh->verts(e.edgeToVert[0]).coord(),hi=mesh->verts(e.edgeToVert[1]).coord();if(component(l,u)==0?lo.y>hi.y:lo.x>hi.x)std::swap(lo,hi);
                        views<<rz<<','<<amr<<','<<n<<','<<selection<<','<<shape<<','<<rho<<",primal,"<<u<<','<<component(l,u)<<','<<f.negativeCell<<','<<f.positiveCell<<','<<lo.x<<','<<lo.y<<','<<hi.x<<','<<hi.y<<','<<metrics.faceMeasures()[u]<<','<<material.faceRates()[u]<<','<<reconstruction.primalAtoms()[u].size()<<'\n';}
                    for(std::size_t i=0;i<geometry.segments().size();++i){const auto& atom=geometry.segments()[i];views<<rz<<','<<amr<<','<<n<<','<<selection<<','<<shape<<','<<rho<<",dual,"<<i<<','<<atom.normalAxis<<','<<atom.negativeVelocity<<','<<atom.positiveVelocity<<','<<atom.begin.x<<','<<atom.begin.y<<','<<atom.end.x<<','<<atom.end.y<<','<<atom.area<<','<<material.dualFaceRates()[i]<<','<<reconstruction.dualAtoms()[i].size()<<'\n';}
                    for(int p=0;p<l.pressureCount();++p){const auto box=m.publication().state().geometry().snapshot().nativeCellBounds(l.nativeCell(p));const auto va=a.cellVolumes()[p],vm=m.cellVolumes()[p],vb=b.cellVolumes()[p];const auto c=material.cells()[p];
                        cells<<rz<<','<<amr<<','<<n<<','<<selection<<','<<shape<<','<<rho<<','<<p<<','<<l.nativeCell(p)<<','<<box[0].x<<','<<box[0].y<<','<<box[1].x<<','<<box[1].y<<','<<metrics.cellVolumes()[p]<<','<<va.liquid<<','<<va.gas<<','<<vm.liquid<<','<<vm.gas<<','<<vb.liquid<<','<<vb.gas<<','<<c.massBegin<<','<<c.massMidpoint<<','<<c.massEnd<<','<<c.massChange<<','<<c.outwardIntegral<<','<<c.integratedDefect<<','<<c.defectRate<<','<<c.relativeDefect<<'\n';}
                    const auto source=material.sources();const auto clock=material.timeLayer();const auto& r=material.report();
                    summary<<rz<<','<<amr<<','<<n<<','<<selection<<','<<shape<<','<<rho<<','<<source.begin.value()<<','<<source.midpoint.value()<<','<<source.end.value()<<','<<source.advector.value()<<','<<clock.step<<','<<clock.begin<<','<<clock.end<<','<<static_cast<int>(clock.position)<<','<<l.velocityCount()<<','<<geometry.segments().size()<<','<<flux.atomRates().size()<<','<<phase.spatialReport().sourceQueries<<','<<reconstruction.report().tangentUnidentified<<','<<r.geometryMass[0]<<','<<r.geometryMass[1]<<','<<r.geometryMass[2]<<','<<r.massChange<<','<<r.boundaryIntegral<<','<<r.integratedDefect<<','<<r.l1IntegratedDefect<<','<<r.balanceIdentity<<','<<sameF<<','<<(&flux.phaseRate().atomRates()==&phase.atomRates())<<','<<flux.momentumGeometry().sameIdentity(geometry)<<','<<seconds<<','<<(rz?0:circle.cx)<<','<<circle.cy<<','<<(rz?0:circle.cx+circle.speed*clock.time())<<','<<circle.cy+(rz?circle.speed*clock.time():0)<<','<<(rz?0:circle.cx+circle.speed*dt)<<','<<circle.cy+(rz?circle.speed*dt:0)<<','<<circle.radius<<','<<plane.planeOffset<<'\n';
                }
                if(!rz && amr && n==16 && selection==0 && shape==1){Field uniform(l.velocityCount());for(int u=0;u<l.velocityCount();++u)uniform[u]=component(l,u)==0?.4:0;
                    for(int layer=0;layer<3;++layer){const auto& p=layer==0?a:layer==1?m:b;const auto v=reconstruction.evaluate(*mesh,CchNativeValues(l,uniform,Field(l.pressureCount())),p.geometryIdentity(),p.timeLayer());const auto rates=CchPhaseRate::prepareSpatial(*mesh,p,v);
                        for(double rho:{1.,1000.,1e6}){const auto flux=CchMaterialFlux::prepare(*mesh,rates,rho,1);const auto clock=flux.timeLayer();
                            for(int side=0;side<4;++side)for(int i:geometry.volumes()[103].sides[side]){const auto& atom=geometry.segments()[i];const auto q=rates.dualFaceRates()[i];check(q.liquid==0,"spatial material gas witness leaked liquid");
                                gas<<layer<<','<<rho<<','<<flux.geometryIdentity().value()<<','<<flux.advectorIdentity().value()<<','<<clock.step<<','<<clock.begin<<','<<clock.end<<','<<static_cast<int>(clock.position)<<",103,"<<side<<','<<i<<','<<atom.normalAxis<<','<<atom.begin.x<<','<<atom.begin.y<<','<<atom.end.x<<','<<atom.end.y<<','<<atom.area<<','<<q.liquid<<','<<q.gas<<','<<flux.dualFaceRates()[i]<<'\n';}}
                    }
                }
            }
        }
        std::cout<<"SPATIAL_MATERIAL_LAYOUT,"<<rz<<','<<amr<<','<<n<<",PASS\n";
    }
}
}
// H10k observer end
// H10l1 observer begin
namespace {
int spatialClosureAudit()
{
    std::filesystem::create_directories("output/series");
    std::ofstream faces("output/series/closure_faces.csv"),cells("output/series/closure_cells.csv"),status("output/series/closure_status.csv");
    std::ofstream atoms("output/series/closure_atoms.csv"),views("output/series/closure_views.csv"),samples("output/series/closure_samples.csv"),sources("output/series/closure_sources.csv");
    check(faces.good() && cells.good() && status.good() && atoms.good() && views.good() && samples.good() && sources.good(),"cannot open actual closure observations");
    faces<<std::setprecision(17)<<"rz,amr,basis,u,axis,ax,ay,bx,by,cx,cy,input,negative,positive,area\n";
    cells<<std::setprecision(17)<<"rz,amr,basis,p,ax,ay,bx,by,volume,original_divergence\n";
    atoms<<std::setprecision(17)<<"rz,amr,basis,atom,axis,primal,owner,ax,ay,bx,by,area,length";
    writeCoefficientHeader<decltype(CchNormalPolynomial::beginCoefficients)>(atoms,"c");
    writeCoefficientHeader<decltype(CchNormalPolynomial::endCoefficients)>(atoms,"e");atoms<<",integral,end_integral\n";
    views<<std::setprecision(17)<<"rz,amr,basis,kind,id,atom,axis,negative,positive,ax,ay,bx,by,area\n";
    samples<<std::setprecision(17)<<"rz,amr,basis,atom,node,pressure,x,y,normal\n";
    sources<<std::setprecision(17)<<"rz,amr,basis,geometry,advector,step,begin,end,position,faces,cells,atoms,segments,view_entries\n";
    status<<"rz,amr,basis,accepted\n";int failures=0,scenes=0;
    for(bool rz:{false,true})for(int amr:{0,1}){auto mesh=makeMesh(4,amr);const auto pub=publication(*mesh,rz,{});FaceSegmentGraph graph;graph.build(*mesh);
        const auto l=CchDofLayout::build(*mesh,NativeLeafDomain::build(*mesh,graph,rz));const auto metric=CchPhysicalDivergence::build(*mesh,CchIntegralDivergence::build(*mesh,l));const auto geometry=CchMomentumGeometry::build(*mesh,metric);const auto reconstruction=CchAdvectorReconstruction::prepare(*mesh,geometry);
        for(int chosen=0;chosen<l.velocityCount();++chosen){if(metric.faceMeasures()[chosen]==0)continue;Field velocity(l.velocityCount());velocity[chosen]=1/metric.faceMeasures()[chosen];const CchNativeValues values(l,velocity,Field(l.pressureCount()));
            for(int u=0;u<l.velocityCount();++u){const auto f=l.domain().face(*mesh,l.nativeFace(u));const auto& edge=f.horizontal?mesh->hEdges(f.edgeLabel):mesh->vEdges(f.edgeLabel);auto a=mesh->verts(edge.edgeToVert[0]).coord(),b=mesh->verts(edge.edgeToVert[1]).coord();if(f.horizontal?a.x>b.x:a.y>b.y)std::swap(a,b);
                faces<<rz<<','<<amr<<','<<chosen<<','<<u<<','<<component(l,u)<<','<<a.x<<','<<a.y<<','<<b.x<<','<<b.y<<','<<f.center.x<<','<<f.center.y<<','<<velocity[u]<<','<<f.negativeCell<<','<<f.positiveCell<<','<<metric.faceMeasures()[u]<<'\n';}
            const auto div=metric.apply(*mesh,values);for(int p=0;p<l.pressureCount();++p){const auto& cell=mesh->cells(mesh->cellIndexToLabel[l.nativeCell(p)]);const auto a=mesh->verts(cell.cellToVert[0]).coord(),b=mesh->verts(cell.cellToVert[2]).coord();cells<<rz<<','<<amr<<','<<chosen<<','<<l.nativeCell(p)<<','<<a.x<<','<<a.y<<','<<b.x<<','<<b.y<<','<<metric.cellVolumes()[p]<<','<<div[p]<<'\n';}
            bool accepted=true;try{const auto actual=reconstruction.evaluate(*mesh,values,pub.state().identity(),pub.state().timeLayer());actual.validateFor(*mesh,geometry,pub.state().identity(),pub.state().timeLayer());
                const auto clock=actual.timeLayer();sources<<rz<<','<<amr<<','<<chosen<<','<<actual.geometryIdentity().value()<<','<<actual.identity().value()<<','<<clock.step<<','<<clock.begin<<','<<clock.end<<','<<static_cast<int>(clock.position)<<','<<l.velocityCount()<<','<<l.pressureCount()<<','<<reconstruction.atoms().size()<<','<<geometry.segments().size()<<','<<reconstruction.report().viewEntries<<'\n';
                // Observe complete public source data. No reconstruction/patch weights are consulted.
                for(std::size_t i=0;i<reconstruction.atoms().size();++i){const auto& a=reconstruction.atoms()[i];const auto& poly=actual.atomPolynomials()[i];
                    atoms<<rz<<','<<amr<<','<<chosen<<','<<i<<','<<a.normalAxis<<','<<a.primalVelocity<<','<<a.ownerPressure<<','<<a.begin.x<<','<<a.begin.y<<','<<a.end.x<<','<<a.end.y<<','<<a.area<<','<<poly.length;writeCoefficients(atoms,poly.beginCoefficients);writeCoefficients(atoms,poly.endCoefficients);atoms<<','<<poly.integral(false,0,poly.length)<<','<<poly.integral(true,0,poly.length)<<'\n';
                    int p=a.ownerPressure;if(p<0){const auto f=l.domain().face(*mesh,l.nativeFace(a.primalVelocity));p=l.pressureDof(f.negativeCell>=0?f.negativeCell:f.positiveCell);}
                    // GL3 integrates the complete metric trace; two extra
                    // points make the value comparison unisolvent for degree 4.
                    const double nodes[5]={-std::sqrt(.6),0,std::sqrt(.6),-.5,.5};
                    for(int j=0;j<5;++j){const auto point=a.begin+(a.end-a.begin)*((1+nodes[j])/2);const auto v=actual.cellVelocity(p,point);samples<<rz<<','<<amr<<','<<chosen<<','<<i<<','<<j<<','<<p<<','<<point.x<<','<<point.y<<','<<(a.normalAxis==0?v.x:v.y)<<'\n';}
                }
                const auto writeView=[&](const char* kind,int id,int axis,int negative,int positive,Vector2D a,Vector2D b,long double area,const std::vector<int>& ids){
                    for(int atom:ids)views<<rz<<','<<amr<<','<<chosen<<','<<kind<<','<<id<<','<<atom<<','<<axis<<','<<negative<<','<<positive<<','<<a.x<<','<<a.y<<','<<b.x<<','<<b.y<<','<<area<<'\n';};
                for(int u=0;u<l.velocityCount();++u){const auto f=l.domain().face(*mesh,l.nativeFace(u));const auto& e=f.horizontal?mesh->hEdges(f.edgeLabel):mesh->vEdges(f.edgeLabel);auto a=mesh->verts(e.edgeToVert[0]).coord(),b=mesh->verts(e.edgeToVert[1]).coord();if(f.horizontal?a.x>b.x:a.y>b.y)std::swap(a,b);
                    writeView("primal",u,component(l,u),f.negativeCell,f.positiveCell,a,b,metric.faceMeasures()[u],reconstruction.primalAtoms()[u]);}
                for(std::size_t i=0;i<geometry.segments().size();++i){const auto& s=geometry.segments()[i];writeView("dual",static_cast<int>(i),s.normalAxis,s.negativeVelocity,s.positiveVelocity,s.begin,s.end,s.area,reconstruction.dualAtoms()[i]);}
            }
            catch(const std::exception& e){accepted=false;++failures;std::cout<<"CLOSURE_REJECT,"<<rz<<','<<amr<<','<<chosen<<','<<e.what()<<std::endl;}
            status<<rz<<','<<amr<<','<<chosen<<','<<accepted<<'\n';++scenes;
        }
    }
    std::cout<<"CLOSURE_AUDIT,"<<scenes<<','<<failures<<'\n';return failures?1:0;
}
void spatialBoundaryAudit()
{
    std::filesystem::create_directories("output/series");std::ofstream atoms("output/series/boundary_atoms.csv"),faces("output/series/boundary_faces.csv"),segments("output/series/boundary_segments.csv"),cells("output/series/boundary_cells.csv"),summary("output/series/boundary_summary.csv");
    check(atoms.good() && faces.good() && segments.good() && cells.good() && summary.good(),"cannot open boundary observations");
    atoms<<std::setprecision(17)<<"rz,amr,n,selection,bound,id,axis,primal,owner,ax,ay,bx,by,area";
    writeCoefficientHeader<decltype(CchNormalPolynomial::beginCoefficients)>(atoms,"c");atoms<<",ql,qg,F\n";
    faces<<std::setprecision(17)<<"rz,amr,n,selection,bound,u,axis,negative,positive,outward,area,mean,F\n";
    segments<<std::setprecision(17)<<"rz,amr,n,selection,bound,id,component,axis,negative,positive,ax,ay,bx,by,area,F\n";
    cells<<std::setprecision(17)<<"rz,amr,n,selection,bound,p,ax,ay,bx,by,volume,DI_div,B_div,Bx_minus,Bx_plus,By_minus,By_plus,hx,hy\n";
    summary<<"rz,amr,n,selection,bound,geometry,advector,step,begin,end,position,faces,cells,atoms,segments,walls,retained_law\n";
    int count=0;for(bool rz:{false,true})for(int amr:{0,1})for(int n:{4,8,16,32}){auto mesh=makeMesh(n,amr);const auto pub=publication(*mesh,rz,{1,0,.1,Position::Midpoint});
        for(int selection=0;selection<2;++selection){FaceSegmentGraph graph;graph.build(*mesh,[selection](const QuadtreeMesh&,const Cell& c){return selection==0 || c.centerCoord().x<.5;});const auto l=CchDofLayout::build(*mesh,NativeLeafDomain::build(*mesh,graph,rz));const auto metric=CchPhysicalDivergence::build(*mesh,CchIntegralDivergence::build(*mesh,l));const auto g=CchMomentumGeometry::build(*mesh,metric);const auto r=CchAdvectorReconstruction::prepare(*mesh,g);const auto phase=CchPhaseGeometry::build(*mesh,metric,pub);std::vector<int> walls;Field values(l.velocityCount());const double right=selection?.5:1;
            const auto psi=[&](Vector2D p){return (rz?p.x*p.x*(right-p.x)*(right-p.x):p.x*(right-p.x))*p.y*(1-p.y);};
            for(int u=0;u<l.velocityCount();++u){const auto f=l.domain().face(*mesh,l.nativeFace(u));const auto& edge=f.horizontal?mesh->hEdges(f.edgeLabel):mesh->vEdges(f.edgeLabel);auto a=mesh->verts(edge.edgeToVert[0]).coord(),b=mesh->verts(edge.edgeToVert[1]).coord();if(f.horizontal?a.x>b.x:a.y>b.y)std::swap(a,b);const auto area=metric.faceMeasures()[u];values[u]=area>0?(f.horizontal?-1:1)*(rz?static_cast<double>(2*pi):1)*(psi(b)-psi(a))/area:0;if(f.outwardSign!=0)walls.push_back(u);}
            const auto law=CchAdvectorBoundary::impermeable(*mesh,g,walls);const CchNativeValues native(l,values,Field(l.pressureCount()));const auto original=metric.apply(*mesh,native);
            for(int bound=0;bound<2;++bound){const auto actual=bound?r.evaluate(*mesh,native,phase.geometryIdentity(),phase.timeLayer(),law):r.evaluate(*mesh,native,phase.geometryIdentity(),phase.timeLayer());const auto rates=CchPhaseRate::prepareSpatial(*mesh,phase,actual);const auto flux=CchMaterialFlux::prepare(*mesh,rates,1,1);
                for(std::size_t i=0;i<r.atoms().size();++i){const auto& a=r.atoms()[i];const auto& c=actual.atomPolynomials()[i].beginCoefficients;const auto q=rates.atomRates()[i];atoms<<rz<<','<<amr<<','<<n<<','<<selection<<','<<bound<<','<<i<<','<<a.normalAxis<<','<<a.primalVelocity<<','<<a.ownerPressure<<','<<a.begin.x<<','<<a.begin.y<<','<<a.end.x<<','<<a.end.y<<','<<a.area;writeCoefficients(atoms,c);atoms<<','<<q.liquid<<','<<q.gas<<','<<flux.atomRates()[i]<<'\n';}
                for(int u=0;u<l.velocityCount();++u){const auto f=l.domain().face(*mesh,l.nativeFace(u));faces<<rz<<','<<amr<<','<<n<<','<<selection<<','<<bound<<','<<u<<','<<component(l,u)<<','<<f.negativeCell<<','<<f.positiveCell<<','<<f.outwardSign<<','<<metric.faceMeasures()[u]<<','<<values[u]<<','<<flux.faceRates()[u]<<'\n';}
                for(std::size_t i=0;i<g.segments().size();++i){const auto& s=g.segments()[i];segments<<rz<<','<<amr<<','<<n<<','<<selection<<','<<bound<<','<<i<<','<<component(l,owner(s))<<','<<s.normalAxis<<','<<s.negativeVelocity<<','<<s.positiveVelocity<<','<<s.begin.x<<','<<s.begin.y<<','<<s.end.x<<','<<s.end.y<<','<<s.area<<','<<flux.dualFaceRates()[i]<<'\n';}
                for(int p=0;p<l.pressureCount();++p){const auto& cell=mesh->cells(mesh->cellIndexToLabel[l.nativeCell(p)]);const auto a=mesh->verts(cell.cellToVert[0]).coord(),b=mesh->verts(cell.cellToVert[2]).coord();const long double xl=rz?static_cast<long double>(a.x)*a.x/2:a.x,xh=rz?static_cast<long double>(b.x)*b.x/2:b.x,dx=xh-xl,dy=b.y-a.y;const auto x=xl+dx/4,y=a.y+dy/4,hx=dx/8,hy=dy/8;
                    const auto sample=[&](long double xx,long double yy,int c){const Vector2D q{static_cast<double>(rz?std::sqrt(2*xx):xx),static_cast<double>(yy)};const auto v=actual.cellVelocity(p,q);return c==0?(rz?q.x:1)*static_cast<long double>(v.x):static_cast<long double>(v.y);};
                    cells<<rz<<','<<amr<<','<<n<<','<<selection<<','<<bound<<','<<l.nativeCell(p)<<','<<a.x<<','<<a.y<<','<<b.x<<','<<b.y<<','<<metric.cellVolumes()[p]<<','<<original[p]<<','<<actual.cellDivergence()[p]<<','<<sample(x-hx,y,0)<<','<<sample(x+hx,y,0)<<','<<sample(x,y-hy,1)<<','<<sample(x,y+hy,1)<<','<<hx<<','<<hy<<'\n';}
                const auto clock=actual.timeLayer();summary<<rz<<','<<amr<<','<<n<<','<<selection<<','<<bound<<','<<actual.geometryIdentity().value()<<','<<actual.identity().value()<<','<<clock.step<<','<<clock.begin<<','<<clock.end<<','<<static_cast<int>(clock.position)<<','<<l.velocityCount()<<','<<l.pressureCount()<<','<<r.atoms().size()<<','<<g.segments().size()<<','<<walls.size()<<','<<(bound?(actual.boundary() && actual.boundary()->sameIdentity(law)):actual.boundary()==nullptr)<<'\n';++count;
            }
        }
    }std::cout<<"BOUNDARY_AUDIT,"<<count<<" actual revisions\n";
}
}
// H10l1 observer end
// H10l2 observer begin
namespace {
void spatialDualMaterialAudit()
{
    std::filesystem::create_directories("output/series");
    std::ofstream segments("output/series/dual_segments.csv"),duals("output/series/dual_cells.csv"),overlaps("output/series/dual_overlaps.csv"),faces("output/series/dual_primal_faces.csv"),cells("output/series/dual_primal_cells.csv"),summary("output/series/dual_summary.csv"),gas("output/series/dual_gas_witness.csv");
    check(segments.good()&&duals.good()&&overlaps.good()&&faces.good()&&cells.good()&&summary.good()&&gas.good(),"cannot open same-F dual observer");
    const auto head="rz,amr,n,selection,mode,rho,";
    segments<<std::setprecision(17)<<head<<"id,component,axis,negative,positive,ax,ay,bx,by,area,liquid,gas,F,I\n";
    duals<<std::setprecision(17)<<head<<"u,component,ax,ay,bx,by,volume,liquid_a,gas_a,liquid_m,gas_m,liquid_b,gas_b,mass_a,mass_m,mass_b,change,outward,defect,R,relative,mapped_change,G,mapped_defect,q,mapped_q,E_sp,impulse,kinetic\n";
    overlaps<<std::setprecision(17)<<head<<"u,p,native,cx0,cy0,cx1,cy1,x0,y0,x1,y1,volume,fraction\n";
    faces<<std::setprecision(17)<<head<<"u,axis,negative,positive,ax,ay,bx,by,area,F,I\n";
    cells<<std::setprecision(17)<<head<<"p,native,volume,change,outward,defect\n";
    summary<<std::setprecision(17)<<head<<"a,m,b,advector,step,begin,end,position,duals,segments,unidentified,same_F,same_geometry,max_E,max_relative_E,max_decomposition,mass_a_0,mass_m_0,mass_b_0,change_0,boundary_0,defect_0,mass_a_1,mass_m_1,mass_b_1,change_1,boundary_1,defect_1,kinetic_0,kinetic_1,cxa,cya,cxm,cym,cxb,cyb,radius,plane\n";
    gas<<std::setprecision(17)<<"layer,rho,geometry,advector,step,begin,end,position,dof,side,segment,axis,ax,ay,bx,by,area,liquid,gas,F\n";
    for(bool rz:{false,true})for(bool amr:{false,true})for(int n:{4,8,16,32}){auto mesh=makeMesh(n,amr?1:0);const double dt=.1/n;std::string stage;
        Case circle;circle.pureFirstJet=true;circle.radial=rz;circle.motion=1;circle.speed=.4;circle.radius=rz?.34:.19;
        const auto ca=initial(*mesh,circle),cm=advanceFirstJet(*mesh,ca,circle,{1,0,dt,Position::Midpoint},stage),cb=advanceFirstJet(*mesh,ca,circle,{1,0,dt,Position::End},stage);
        Case plane;plane.pureFirstJet=true;plane.radial=rz;plane.planar=true;plane.planeOffset=.43;const auto pa=initial(*mesh,plane),pm=advanceFirstJet(*mesh,pa,plane,{1,0,dt,Position::Midpoint},stage),pb=advanceFirstJet(*mesh,pa,plane,{1,0,dt,Position::End},stage);
        for(int selection=0;selection<4;++selection){FaceSegmentGraph graph;graph.build(*mesh,[selection](const QuadtreeMesh&,const Cell& c){const auto p=c.centerCoord();if(selection==1)return p.x<.5;if(selection==2)return !(p.x>.25&&p.x<.75&&p.y>.25&&p.y<.75);if(selection==3)return p.x<.25||p.x>.75;return true;});
            const auto l=CchDofLayout::build(*mesh,NativeLeafDomain::build(*mesh,graph,rz));const auto metrics=CchPhysicalDivergence::build(*mesh,CchIntegralDivergence::build(*mesh,l));const auto geometry=CchMomentumGeometry::build(*mesh,metrics);const auto reconstruction=CchAdvectorReconstruction::prepare(*mesh,geometry);
            for(int mode=0;mode<3;++mode){const auto a=CchPhaseGeometry::build(*mesh,metrics,mode==0?pa:ca),m=CchPhaseGeometry::build(*mesh,metrics,mode==0?pm:cm),b=CchPhaseGeometry::build(*mesh,metrics,mode==0?pb:cb);
                Field v(l.velocityCount());for(int u=0;u<l.velocityCount();++u){const auto f=l.domain().face(*mesh,l.nativeFace(u));if(mode==2)v[u]=component(l,u)==(rz?1:0)?.4:0;
                    else if(!rz)v[u]=component(l,u)==0?.3+.17*f.center.x-.23*f.center.y:-.2+.31*f.center.x+.19*f.center.y;
                    else if(component(l,u)==0)v[u]=.13*f.center.x;else{const auto& e=mesh->hEdges(f.edgeLabel);const auto lo=mesh->verts(e.edgeToVert[0]).coord(),hi=mesh->verts(e.edgeToVert[1]).coord();v[u]=.2+.155*.5*(lo.x*lo.x+hi.x*hi.x)+.19*f.center.y;}}
                const auto values=CchNativeValues(l,v,Field(l.pressureCount()));const auto actual=reconstruction.evaluate(*mesh,values,m.geometryIdentity(),m.timeLayer());const auto phase=CchPhaseRate::prepareSpatial(*mesh,m,actual);
                for(double rho:{1.,1000.,1e6}){const auto flux=CchMaterialFlux::prepare(*mesh,phase,rho,1);const auto material=CchMaterialTransport::prepare(*mesh,a,flux,b);
                    const auto d=CchDualMaterial::prepare(*mesh,material,CchMomentumMass::build(*mesh,geometry,a,rho,1),CchMomentumMass::build(*mesh,geometry,m,rho,1),CchMomentumMass::build(*mesh,geometry,b,rho,1));d.validateFor(*mesh,d.sources(),d.timeLayer());const auto impulse=d.consistencyImpulse(*mesh,values,d.sources(),d.timeLayer());const auto kinetic=d.materialKineticChange(*mesh,values,d.sources(),d.timeLayer());
                    check(&d.faceRates()==&flux.dualFaceRates()&&d.geometry().sameIdentity(flux.momentumGeometry()),"same-F dual observer detached owner/source");
                    const auto prefix=std::to_string(rz)+","+std::to_string(amr)+","+std::to_string(n)+","+std::to_string(selection)+","+std::to_string(mode)+","+std::to_string(static_cast<int>(rho))+",";
                    for(std::size_t id=0;id<geometry.segments().size();++id){const auto& s=geometry.segments()[id];const auto q=phase.dualFaceRates()[id];segments<<prefix<<id<<','<<component(l,owner(s))<<','<<s.normalAxis<<','<<s.negativeVelocity<<','<<s.positiveVelocity<<','<<s.begin.x<<','<<s.begin.y<<','<<s.end.x<<','<<s.end.y<<','<<s.area<<','<<q.liquid<<','<<q.gas<<','<<d.faceRates()[id]<<','<<d.faceIntegrals()[id]<<'\n';}
                    for(int u=0;u<l.velocityCount();++u){const auto f=l.domain().face(*mesh,l.nativeFace(u));const auto& e=f.horizontal?mesh->hEdges(f.edgeLabel):mesh->vEdges(f.edgeLabel);auto lo=mesh->verts(e.edgeToVert[0]).coord(),hi=mesh->verts(e.edgeToVert[1]).coord();if(f.horizontal?lo.x>hi.x:lo.y>hi.y)std::swap(lo,hi);
                        faces<<prefix<<u<<','<<component(l,u)<<','<<f.negativeCell<<','<<f.positiveCell<<','<<lo.x<<','<<lo.y<<','<<hi.x<<','<<hi.y<<','<<metrics.faceMeasures()[u]<<','<<material.faceRates()[u]<<','<<material.faceIntegrals()[u]<<'\n';
                        const auto va=d.beginMass().phaseVolumes()[u],vm=d.midpointMass().phaseVolumes()[u],vb=d.endMass().phaseVolumes()[u];const auto c=d.cells()[u];
                        duals<<prefix<<u<<','<<component(l,u)<<','<<lo.x<<','<<lo.y<<','<<hi.x<<','<<hi.y<<','<<geometry.volumes()[u].volume<<','<<va.liquid<<','<<va.gas<<','<<vm.liquid<<','<<vm.gas<<','<<vb.liquid<<','<<vb.gas<<','<<d.beginMass().masses()[u]<<','<<d.midpointMass().masses()[u]<<','<<d.endMass().masses()[u]<<','<<c.massChange<<','<<c.outwardIntegral<<','<<c.integratedDefect<<','<<c.defectRate<<','<<c.relativeDefect<<','<<c.mappedPrimalMassChange<<','<<c.geometryTransferChange<<','<<c.mappedPrimalIntegratedDefect<<','<<c.outwardRate<<','<<c.mappedPrimalOutwardRate<<','<<c.spatialExchangeRate<<','<<impulse[u]<<','<<(.5L*v[u])*impulse[u]<<'\n';
                        for(const auto& p:geometry.volumes()[u].overlaps){const auto box=m.publication().state().geometry().snapshot().nativeCellBounds(l.nativeCell(p.pressure));overlaps<<prefix<<u<<','<<p.pressure<<','<<l.nativeCell(p.pressure)<<','<<box[0].x<<','<<box[0].y<<','<<box[1].x<<','<<box[1].y<<','<<p.lower.x<<','<<p.lower.y<<','<<p.upper.x<<','<<p.upper.y<<','<<p.volume<<','<<p.cellFraction<<'\n';}}
                    for(int p=0;p<l.pressureCount();++p){const auto c=material.cells()[p];cells<<prefix<<p<<','<<l.nativeCell(p)<<','<<metrics.cellVolumes()[p]<<','<<c.massChange<<','<<c.outwardIntegral<<','<<c.integratedDefect<<'\n';}
                    const auto source=d.sources();const auto clock=d.timeLayer();const auto report=d.report();summary<<prefix<<source.begin.value()<<','<<source.midpoint.value()<<','<<source.end.value()<<','<<source.advector.value()<<','<<clock.step<<','<<clock.begin<<','<<clock.end<<','<<static_cast<int>(clock.position)<<','<<l.velocityCount()<<','<<geometry.segments().size()<<','<<reconstruction.report().tangentUnidentified<<",1,1,"<<report.maxAbsSpatialExchange<<','<<report.maxRelativeSpatialExchange<<','<<report.maxDefectDecomposition;
                    for(int c=0;c<2;++c){const auto r=report.components[c];summary<<','<<r.geometryMass[0]<<','<<r.geometryMass[1]<<','<<r.geometryMass[2]<<','<<r.massChange<<','<<r.boundaryIntegral<<','<<r.integratedDefect;}
                    summary<<','<<kinetic[0]<<','<<kinetic[1]<<','<<(rz?0:circle.cx)<<','<<circle.cy<<','<<(rz?0:circle.cx+circle.speed*clock.time())<<','<<circle.cy+(rz?circle.speed*clock.time():0)<<','<<(rz?0:circle.cx+circle.speed*dt)<<','<<circle.cy+(rz?circle.speed*dt:0)<<','<<circle.radius<<','<<plane.planeOffset<<'\n';
                }
                if(!rz&&amr&&n==16&&selection==0&&mode==2)for(int layer=0;layer<3;++layer){const auto p=layer==0?a:layer==1?m:b;const auto vLayer=reconstruction.evaluate(*mesh,values,p.geometryIdentity(),p.timeLayer());const auto q=CchPhaseRate::prepareSpatial(*mesh,p,vLayer);
                    for(double rho:{1.,1000.,1e6}){const auto f=CchMaterialFlux::prepare(*mesh,q,rho,1);const auto clock=f.timeLayer();for(int side=0;side<4;++side)for(int id:geometry.volumes()[103].sides[side]){const auto& s=geometry.segments()[id];gas<<layer<<','<<rho<<','<<f.geometryIdentity().value()<<','<<f.advectorIdentity().value()<<','<<clock.step<<','<<clock.begin<<','<<clock.end<<','<<static_cast<int>(clock.position)<<",103,"<<side<<','<<id<<','<<s.normalAxis<<','<<s.begin.x<<','<<s.begin.y<<','<<s.end.x<<','<<s.end.y<<','<<s.area<<','<<q.dualFaceRates()[id].liquid<<','<<q.dualFaceRates()[id].gas<<','<<f.dualFaceRates()[id]<<'\n';}}}
            }
        }
        std::cout<<"SAME_F_DUAL_LAYOUT,"<<rz<<','<<amr<<','<<n<<",PASS\n";
    }
}
}
// H10l2 observer end
int main(int argc,char** argv)
{
    try{
        if(argc==2 && (std::string(argv[1])=="controlled" || std::string(argv[1])=="controlled-probe")){
            std::cout<<std::setprecision(17);controlledStudy(std::string(argv[1])=="controlled-probe");return failed?1:0;}
        if(argc==2 && std::string(argv[1])=="--help")std::cout<<"controlled: H11 full 36-window local-dissipation study; controlled-probe: development n8 only\n";
        // H10l2 observer begin
        if(argc==2&&std::string(argv[1])=="dual"){spatialDualMaterialAudit();return 0;}
        // H10l2 observer end
        // H10k observer begin
        // H10l1 observer begin
        if(argc==2 && std::string(argv[1])=="closure"){return spatialClosureAudit();}
        if(argc==2 && std::string(argv[1])=="boundary"){spatialBoundaryAudit();return 0;}
        // H10l1 observer end
        if(argc==2 && std::string(argv[1])=="material"){spatialMaterialAudit();return 0;}
        if(argc==2 && std::string(argv[1])=="--help")std::cout<<"material: read-only H10k actual unique spatial F/source study\n";
        // H10k observer end
        // H10j observer begin
        if(argc==2 && std::string(argv[1])=="spatial"){spatialPhaseAudit();return 0;}
        if(argc==2 && std::string(argv[1])=="--help")std::cout<<"spatial: read-only H10j actual causal phase owner study\n";
        // H10j observer end
        // H10j0 observer begin
        if(argc==2 && std::string(argv[1])=="advector"){spatialAdvectorAudit();return 0;}
        if(argc==2 && std::string(argv[1])=="--help")std::cout<<"advector: read-only H10j0 compatible velocity owner study\n";
        // H10j0 observer end
        if(argc==2 && std::string(argv[1])=="--help"){std::cout<<"LAMP_cch_convection_benchmark donor-audit|full|probe|audit|audit32|rows|rows32|compatibility|compatibility32|certificate|certificate32|bounds|bounds32|intervals|intervals32|segments|segments32|coincident|trace-contract\nUse tools/run_case.cmd benchmarks/cch_convection for the fixed full study.\n";return 0;}
        check(argc==2 && (std::string(argv[1])=="donor-audit" || std::string(argv[1])=="full" || std::string(argv[1])=="probe" || std::string(argv[1])=="audit" || std::string(argv[1])=="audit32" || std::string(argv[1])=="rows" || std::string(argv[1])=="rows32" || std::string(argv[1])=="compatibility" || std::string(argv[1])=="compatibility32" || std::string(argv[1])=="certificate" || std::string(argv[1])=="certificate32" || std::string(argv[1])=="bounds" || std::string(argv[1])=="bounds32" || std::string(argv[1])=="intervals" || std::string(argv[1])=="intervals32" || std::string(argv[1])=="segments" || std::string(argv[1])=="segments32" || std::string(argv[1])=="coincident" || std::string(argv[1])=="trace-contract"),"usage: LAMP_cch_convection_benchmark donor-audit|full|probe|audit|audit32|rows|rows32|compatibility|compatibility32|certificate|certificate32|bounds|bounds32|intervals|intervals32|segments|segments32|coincident|trace-contract");
        std::filesystem::create_directories("output/series");
        summaries.open("output/series/space.csv");history.open("output/series/trajectory.csv");gates.open("output/series/gates.csv");
        check(summaries.good() && history.good() && gates.good(),"cannot open benchmark outputs");
        std::cout<<std::setprecision(17);summaries<<std::setprecision(17);history<<std::setprecision(17);gates<<std::setprecision(17);
        summaries<<"rz,amr,n,velocity_dofs,coarse_fine_segments,steps,truncation_l2,terminal_l2,terminal_linf,relative_energy_drift,donor_energy_loss,boundary_power,material_power,R,mass_change,G,raw_axis_acceleration\n";
        history<<"method,rz,amr,n,steps,step,tau,relative_l2,linf,relative_energy_drift\n";gates<<"gate,status,value,limit\n";
        // H10n: one explicit gate version, with legacy failures kept separately.
        legacyGates.open("output/series/legacy_gates.csv");donorPower.open("output/series/donor_power.csv");
        donorSegments.open("output/series/donor_energy_segments.csv");donorNodes.open("output/series/donor_energy_nodes.csv");
        std::ofstream contract("output/series/gate_contract.csv");
        check(legacyGates.good() && donorPower.good() && donorSegments.good() && donorNodes.good() && contract.good(),"cannot open versioned donor diagnostics");
        legacyGates<<std::setprecision(17)<<"gate,status,value,limit\n";
        donorPower<<std::setprecision(17)<<"rz,amr,n,sample,energy,internal_D,periodic_D,open_power,material_correction,power,gross,identity_error,fixed_error,decay_rate\n";
        donorSegments<<std::setprecision(17)<<"rz,amr,n,sample,segment,partner,kind,first,second,rate,other_rate,D\n";
        donorNodes<<std::setprecision(17)<<"rz,amr,n,sample,dof,fixed_axis,mass,u,rhs,q,R,mass_change_rate\n";
        contract<<"version,legacy_gate,current_gate,minimum,maximum,preserved_other_full_gates\n";
        for(bool rz:{false,true})for(bool amr:{false,true}){const auto name=std::string(rz?"rz":"xy")+(amr?"_amr":"_uniform");
            contract<<"H10b-dissipation-v2,"<<name<<"_donor_loss_order,"<<name<<"_donor_decay_order,0.7,1.4,80\n";}
        contract.close();std::cout<<"BENCHMARK_CONTRACT,H10b-dissipation-v2,84,80,4\n";
        if(std::string(argv[1])=="donor-audit")donorAudit();
        else if(std::string(argv[1])=="audit")affineAudit({8,16});
        else if(std::string(argv[1])=="audit32")affineAudit({32});
        else if(std::string(argv[1])=="rows")rowAudit({8,16});
        else if(std::string(argv[1])=="rows32")rowAudit({32});
        else if(std::string(argv[1])=="compatibility")materialCompatibility({4,8,16});
        else if(std::string(argv[1])=="compatibility32")materialCompatibility({32});
        else if(std::string(argv[1])=="bounds")materialCompatibility({4,8,16},true);
        else if(std::string(argv[1])=="bounds32")materialCompatibility({32},true);
        else if(std::string(argv[1])=="intervals")materialCompatibility({4,8,16},true,true);
        else if(std::string(argv[1])=="intervals32")materialCompatibility({32},true,true);
        else if(std::string(argv[1])=="segments")materialCompatibility({4,8,16},true,false,true);
        else if(std::string(argv[1])=="segments32")materialCompatibility({32},true,false,true);
        else if(std::string(argv[1])=="coincident")coincidentSourceAudit();
        else if(std::string(argv[1])=="trace-contract")traceContractAudit();
        else if(std::string(argv[1])=="certificate")centerCertificate({8,16});
        else if(std::string(argv[1])=="certificate32")centerCertificate({32});
        else if(std::string(argv[1])=="probe"){Study s(8,false,true);spatial(s);}
        else full();
        return failed?1:0;
    }catch(const std::exception& e){std::cerr<<"H10b research ERROR: "<<e.what()<<std::endl;return 2;}
}
