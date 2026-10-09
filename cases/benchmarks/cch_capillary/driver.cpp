// H13 research: same-source capillary force and actual native/FirstJet motion.
// The capillary-only Strang split is a research integrator, not H14 acceptance.
#include "Modules/Hydrodynamics/CchInterfaceNormalLoad.h"
#include "Modules/Hydrodynamics/CchDualMaterial.h"
#include "Modules/Hydrodynamics/CchPhaseViscosity.h"
#include "Modules/Hydrodynamics/CchSinglePhaseSaddle.h"
#include "interface_first_jet_fixture.h"
#include <filesystem>
#include <fstream>
#include <chrono>

namespace {
using namespace Quantum::Hydrodynamics;
using Matrix=CchSaddleVelocityBlock::Matrix;
using Clock=std::chrono::steady_clock;
CchPhysicalDivergence metric(const QuadtreeMesh& mesh,bool rz=false){FaceSegmentGraph graph;graph.build(mesh);return CchPhysicalDivergence::build(mesh,CchIntegralDivergence::build(mesh,CchDofLayout::build(mesh,NativeLeafDomain::build(mesh,graph,rz))));}
CchAdvectorBoundary closedBoundary(const QuadtreeMesh& mesh,const CchMomentumGeometry& geometry){const auto& l=geometry.layout();std::vector<int> walls;
    for(int u=0;u<l.velocityCount();++u){const auto f=l.domain().face(mesh,l.nativeFace(u));if((f.negativeCell<0||f.positiveCell<0)&&geometry.metrics().faceMeasures()[u]>0)walls.push_back(u);}return CchAdvectorBoundary::impermeable(mesh,geometry,walls);}
HermiteGalsPublication waveInitial(const QuadtreeMesh& mesh,double amplitude){
    const double k=double(2*pi);const Profile f=[=](Vector2D p){return Jet{p.y-.5-amplitude*std::cos(k*p.x),(p.x==0||p.x==1)?0:amplitude*k*std::sin(k*p.x),1,0};};
    const std::vector<LevelSetInterfacePanels::SymmetryPlane> planes{{0,0},{0,1}};
    const auto seed=stateFor(mesh,f,false,{},[](int,Vector2D,Vector2D){return 0.;},planes,HermiteGalsTransport::FirstJet);
    HermiteGalsCandidate candidate(mesh,{},seed.primary());const auto id=candidate.identity();const auto mixed=closeHermiteGalsMixedAmr(mesh,candidate,planes);
    HermiteGalsPublicationConfig cfg;cfg.volume.symmetryPlanes=planes;return publishHermiteGalsFirstJetGeometry(mesh,std::move(candidate),mixed,{id,centers(mesh,f)},cfg);
}
struct Surface {double length=0,amplitude=0,mean=0,modeNorm=0,curvatureL2=0;};
// Passive raw exports. These observations never enter a force or update.
struct WaveOutput {
    std::ofstream nodes{"output/series/wave_surface.csv"},faces{"output/series/wave_faces.csv"},
        cells{"output/series/wave_cells.csv"},exchange{"output/series/wave_exchange.csv"},
        jets{"output/series/wave_primary.csv"};
    WaveOutput(){
        for(auto* f:{&nodes,&faces,&cells,&exchange,&jets})*f<<std::setprecision(17);
        nodes<<"n,amr,steps,mu,step,time,component,graph,node,cell,x,y,nx,ny,curvature,weight,vx,vy,jxx,jxy,jyx,jyy\n";
        faces<<"n,amr,steps,mu,step,time,u,minus,plus,horizontal,x0,y0,x1,y1,area,velocity,mass,force,support,curvature\n";
        cells<<"n,amr,steps,mu,step,time,p,native,x0,y0,x1,y1,volume,liquid,pressure\n";
        exchange<<"n,amr,steps,mu,step,time,surface_power,mechanical_power,exchange_defect,wall_tangent_normal_derivative\n";
        jets<<"n,amr,steps,mu,step,time,kind,index,a,b,phi,dx,dy,mixed\n";
    }
    void write(int n,bool amr,int steps,double mu,int step,const QuadtreeMesh& mesh,
        const CchPhaseGeometry& phase,const CchPressureLoad& pressure,const CchMomentumMass& mass,
        const CchSpatialAdvector& field){
        const auto& l=phase.layout();const auto& pub=phase.publication();const double time=phase.timeLayer().time();
        const auto prefix=[&](std::ostream& out){out<<n<<','<<amr<<','<<steps<<','<<mu<<','<<step<<','<<time<<',';};
        const auto load=CchInterfaceNormalLoad::capillary(mesh,pressure,phase,1,field.reconstruction(),field.boundary());
        double surfacePower=0,mechanicalPower=0,wallDerivative=0;
        for(std::size_t c=0;c<pub.topology().curves.size();++c)for(std::size_t g=0;g<pub.topology().curves[c].graphs().size();++g){int node=0;for(const auto& q:pub.quadrature(c,g)){
            const auto& f=q.sample.frame;const auto v=field.cellVelocityJet(l.pressureDof(f.nativeCellIndex),f.point).jet;
            surfacePower+=q.physicalWeight*q.sample.curvature*(f.normal.x*v.value.x+f.normal.y*v.value.y);
            prefix(nodes);nodes<<c<<','<<g<<','<<node++<<','<<f.nativeCellIndex<<','<<f.point.x<<','<<f.point.y<<','<<f.normal.x<<','<<f.normal.y<<','<<q.sample.curvature<<','<<q.physicalWeight<<','<<v.value.x<<','<<v.value.y<<','<<v.derivative.xx<<','<<v.derivative.xy<<','<<v.derivative.yx<<','<<v.derivative.yy<<'\n';}}
        for(int u=0;u<l.velocityCount();++u){const auto f=l.domain().face(mesh,l.nativeFace(u));const auto& edge=f.horizontal?mesh.hEdges(f.edgeLabel):mesh.vEdges(f.edgeLabel);const auto a=mesh.verts(edge.edgeToVert[0]).coord(),b=mesh.verts(edge.edgeToVert[1]).coord();
            mechanicalPower-=load.balanceLoads()[u]*field.values().normalVelocity()[u];
            if(!f.horizontal&&(f.negativeCell<0||f.positiveCell<0))wallDerivative=std::max(wallDerivative,std::abs(field.cellVelocityJet(l.pressureDof(f.negativeCell<0?f.positiveCell:f.negativeCell),(a+b)*.5).jet.derivative.yx));
            prefix(faces);faces<<u<<','<<(f.negativeCell<0?-1:l.pressureDof(f.negativeCell))<<','<<(f.positiveCell<0?-1:l.pressureDof(f.positiveCell))<<','<<f.horizontal<<','<<a.x<<','<<a.y<<','<<b.x<<','<<b.y<<','<<phase.metrics().faceMeasures()[u]<<','<<field.values().normalVelocity()[u]<<','<<mass.masses()[u]<<','<<load.balanceLoads()[u]<<','<<load.supportMeasures()[u]<<','<<load.curvatureAverages()[u]<<'\n';}
        for(int p=0;p<l.pressureCount();++p){const auto& cell=mesh.cells(mesh.cellIndexToLabel[l.nativeCell(p)]);const auto a=mesh.verts(cell.cellToVert[0]).coord(),b=mesh.verts(cell.cellToVert[2]).coord();prefix(cells);cells<<p<<','<<l.nativeCell(p)<<','<<a.x<<','<<a.y<<','<<b.x<<','<<b.y<<','<<phase.metrics().cellVolumes()[p]<<','<<phase.cellVolumes()[p].liquid<<','<<field.values().pressure()[p]<<'\n';}
        prefix(exchange);exchange<<surfacePower<<','<<mechanicalPower<<','<<surfacePower+mechanicalPower<<','<<wallDerivative<<'\n';
        const auto& state=pub.state();const auto& primary=state.primary();const auto& closure=state.closure();
        for(std::size_t i=0;i<primary.vertices.size();++i){const auto& v=primary.vertices[i];prefix(jets);jets<<0<<','<<i<<','<<v.vertex<<",-1,"<<v.phi<<','<<v.dx<<','<<v.dy<<','<<closure.vertexMixedDerivatives[i]<<'\n';}
        for(std::size_t i=0;i<primary.edges.size();++i){const auto& e=primary.edges[i];prefix(jets);jets<<1<<','<<i<<','<<e.vertices[0]<<','<<e.vertices[1]<<','<<e.phi<<','<<e.dn<<",0,0\n";}
        for(std::size_t i=0;i<primary.geometryCells.size();++i){prefix(jets);jets<<2<<','<<i<<','<<primary.geometryCells[i]<<",-1,"<<closure.cellBubbleCoefficients[i]<<",0,0,0\n";}
        for(auto* f:{&nodes,&faces,&cells,&exchange,&jets}){f->flush();check(bool(*f),"wave raw export write failure");}
    }
};
Surface surface(const HermiteGalsPublication& pub){Surface s;const double k=double(2*pi);double curvature=0;
    for(std::size_t c=0;c<pub.topology().curves.size();++c)for(std::size_t g=0;g<pub.topology().curves[c].graphs().size();++g)for(const auto& q:pub.quadrature(c,g)){
        const auto p=q.sample.frame.point;const double dx=q.physicalWeight*q.sample.frame.normal.y,mode=std::cos(k*p.x);s.length+=q.physicalWeight;s.amplitude+=(p.y-.5)*mode*dx;s.mean+=(p.y-.5)*dx;s.modeNorm+=mode*mode*dx;curvature+=q.physicalWeight*q.sample.curvature*q.sample.curvature;}
    s.amplitude/=s.modeNorm;s.curvatureL2=std::sqrt(curvature/s.length);return s;
}
struct Kick {CchSaddleResult result;double capillaryWork=0,viscousWork=0,kineticChange=0;};
Kick kick(const QuadtreeMesh& mesh,const CchPhaseGeometry& phase,const CchPressureLoad& pressure,const CchMomentumMass& mass,
    const std::vector<double>& old,double halfStep,double mu,const CchAdvectorReconstruction& source,const CchAdvectorBoundary& law){
    const auto load=CchInterfaceNormalLoad::capillary(mesh,pressure,phase,1,source,&law);const int nu=phase.layout().velocityCount();
    Matrix viscosity(nu,nu);if(mu>0)viscosity=CchPhaseViscosity::build(mesh,phase,HermiteGalsPhaseQuadrature::build(mesh,phase.publication()),mu,mu).viscosity().matrix();
    Matrix block=.5*viscosity;for(int u=0;u<nu;++u)block.coeffRef(u,u)+=mass.masses()[u]/halfStep;
    const Eigen::Map<const Eigen::VectorXd> before(old.data(),nu);const Eigen::VectorXd previousViscosity=viscosity*before;
    std::vector<double> rhs(nu);for(int u=0;u<nu;++u)rhs[u]=mass.masses()[u]*old[u]/halfStep-.5*previousViscosity[u]-load.balanceLoads()[u];
    const auto saddle=CchSinglePhaseSaddle::buildClosed(pressure,CchSaddleVelocityBlock::build(phase.layout(),block));
    Kick out{saddle.solve(mesh,CchSaddleRhs(phase.layout(),rhs,std::vector<double>(nu)))};
    Eigen::VectorXd average(nu);for(int u=0;u<nu;++u){average[u]=.5*(old[u]+out.result.values().normalVelocity()[u]);out.capillaryWork-=halfStep*load.balanceLoads()[u]*average[u];out.kineticChange+=.5*mass.masses()[u]*(out.result.values().normalVelocity()[u]*out.result.values().normalVelocity()[u]-old[u]*old[u]);}
    out.viscousWork=-halfStep*average.dot(viscosity*average);return out;
}
void wave(int n,bool amr,int steps,double periods,double amplitude,double mu,std::ofstream& history,std::ofstream& summary,WaveOutput& raw){
    const auto started=Clock::now();auto mesh=makeMesh(n,amr?1:0);auto original=waveInitial(*mesh,amplitude);
    const auto metrics=metric(*mesh);const auto geometry=CchMomentumGeometry::build(*mesh,metrics);const auto pressure=CchPressureLoad::build(metrics.divergence());
    const auto reconstruction=CchAdvectorReconstruction::prepare(*mesh,geometry);const auto search=CchCharacteristicGeometry::prepare(*mesh,reconstruction);
    const auto& l=geometry.layout();std::vector<int> walls;for(int u=0;u<l.velocityCount();++u){const auto f=l.domain().face(*mesh,l.nativeFace(u));if(f.negativeCell<0||f.positiveCell<0)walls.push_back(u);}
    const auto law=CchAdvectorBoundary::impermeable(*mesh,geometry,walls);
    const double k=double(2*pi),omega=std::sqrt(k*k*k*std::tanh(k*.5)/2),duration=periods*double(2*pi)/omega,dt=duration/steps;
    std::vector<double> velocity(l.velocityCount());auto phase=CchPhaseGeometry::build(*mesh,metrics,original);auto mass=CchMomentumMass::build(*mesh,geometry,phase,1,1);
    const auto initialSurface=surface(original);const double initialEnergy=initialSurface.length-1;double work=0,dissipation=0,peakEnergy=0,peakContinuity=0,peakR=0,peakUR=0;
    double previousAmplitude=initialSurface.amplitude,crossing=-1;HermiteGalsPublicationConfig publication;publication.volume.symmetryPlanes={{0,0},{0,1}};
    CchSpaceTimeConfig config;config.initialPanels=1;config.characteristic.maxStep=dt;
    raw.write(n,amr,steps,mu,0,*mesh,phase,pressure,mass,reconstruction.evaluate(*mesh,CchNativeValues(l,velocity,std::vector<double>(l.pressureCount())),phase.geometryIdentity(),phase.timeLayer(),law));
    std::cout<<"WAVE_BEGIN,"<<n<<','<<amr<<','<<steps<<','<<periods<<','<<amplitude<<','<<mu<<','<<duration<<std::endl;
    for(int step=0;step<steps;++step){const double begin=duration*step/steps,end=duration*(step+1)/steps;
        const auto first=kick(*mesh,phase,pressure,mass,velocity,(end-begin)/2,mu,reconstruction,law);const auto& middleVelocity=first.result.values();
        const auto beginField=reconstruction.evaluate(*mesh,middleVelocity,original.state().identity(),original.state().timeLayer(),law);
        const auto endField=reconstruction.evaluate(*mesh,middleVelocity,original.state().identity(),{step+1,begin,end,Position::End},law);
        const auto motion=search.linear(*mesh,beginField,endField);
        // Boundary values are traced from the SAME old scalar and actual
        // velocity. The declared side-wall even symmetry supplies dphi/dx=0;
        // no analytic wave amplitude or future interface enters this data.
        const HermiteGalsBoundaryFirstJetQuery boundary=[&](Vector2D p,HermiteGalsTimeLayer t,Boundary){
            const auto trace=motion.traceBackward(p,t.time(),begin,config.characteristic);check(trace.jacobian().has_value(),"wall drift lacks actual derivative");
            const auto j=original.state().geometry().snapshot().scalarJetNative(trace.point());auto gradient=trace.jacobian()->transposeTimes({j.dx,j.dy});
            if(p.x==mesh->minCoord.x||p.x==mesh->maxCoord.x)gradient.x=0;return HermiteGalsBoundaryFirstJet{j.value,gradient};};
        const auto material=CchMaterialTransport::prepareSpaceTime(*mesh,original,motion,1,1,config,boundary,publication);
        const auto midMass=CchMomentumMass::build(*mesh,geometry,material.midpointRate().geometry(),1,1);const auto nextMass=CchMomentumMass::build(*mesh,geometry,material.endGeometry(),1,1);
        const auto dual=CchDualMaterial::prepare(*mesh,material,mass,midMass,nextMass);
        const auto second=kick(*mesh,material.endGeometry(),pressure,nextMass,middleVelocity.normalVelocity(),(end-begin)/2,mu,reconstruction,law);
        double kinetic=0,materialPower=0,l1R=0;for(int u=0;u<l.velocityCount();++u){const auto value=second.result.values().normalVelocity()[u];kinetic+=.5*nextMass.masses()[u]*value*value;
            const double r=double(dual.cells()[u].defectRate);l1R+=std::abs(r);materialPower+=.5*r*middleVelocity.normalVelocity()[u]*middleVelocity.normalVelocity()[u];}
        velocity=second.result.values().normalVelocity();original=material.endGeometry().publication();phase=material.endGeometry();mass=nextMass;
        raw.write(n,amr,steps,mu,step+1,*mesh,phase,pressure,mass,reconstruction.evaluate(*mesh,second.result.values(),phase.geometryIdentity(),phase.timeLayer(),law));
        const auto s=surface(original);work+=first.capillaryWork+second.capillaryWork;dissipation+=first.viscousWork+second.viscousWork;
        const double energy=kinetic+s.length-1,energyDefect=(energy-initialEnergy-dissipation)/initialEnergy;
        peakEnergy=std::max(peakEnergy,std::abs(energyDefect));peakContinuity=std::max({peakContinuity,first.result.report().maxCellContinuityDefect,second.result.report().maxCellContinuityDefect});
        peakR=std::max(peakR,l1R);peakUR=std::max(peakUR,std::abs(materialPower));
        if(crossing<0&&previousAmplitude>0&&s.amplitude<=0)crossing=begin+(end-begin)*previousAmplitude/(previousAmplitude-s.amplitude);previousAmplitude=s.amplitude;
        history<<n<<','<<amr<<','<<steps<<','<<periods<<','<<amplitude<<','<<mu<<','<<step+1<<','<<end<<','<<s.amplitude<<','<<s.mean<<','<<s.length<<','<<kinetic<<','<<work<<','<<dissipation<<','<<energyDefect<<','<<s.curvatureL2<<','<<material.report().l1IntegratedDefect/(end-begin)<<','<<l1R<<','<<materialPower<<','<<peakContinuity<<'\n';history.flush();
        if(step%8==0||step+1==steps)std::cout<<"WAVE_STEP,"<<n<<','<<amr<<','<<step+1<<','<<end<<','<<s.amplitude<<','<<energyDefect<<','<<l1R<<std::endl;
    }
    const double measured=crossing>0?double(pi)/(2*crossing):0;const auto finalSurface=surface(original);
    summary<<n<<','<<amr<<','<<steps<<','<<periods<<','<<amplitude<<','<<mu<<','<<omega<<','<<measured<<','<<finalSurface.amplitude<<','<<peakEnergy<<','<<peakContinuity<<','<<peakR<<','<<peakUR<<','<<std::chrono::duration<double>(Clock::now()-started).count()<<'\n';summary.flush();
    std::cout<<"WAVE_COMPLETE,"<<n<<','<<amr<<','<<steps<<','<<measured<<','<<omega<<','<<peakEnergy<<std::endl;
}
std::unique_ptr<QuadtreeMesh> dropMesh(int n,bool rz,bool amr){
    DictionaryMesh d;d.Nx=d.Ny=n;d.maxCellNum=32768;d.maxInitialIterations=0;auto mesh=std::make_unique<QuadtreeMesh>(d);
    if(amr){DictionaryMesh::RefineBox b;b.min={rz?0:.25,.25};b.max={rz?.5:.75,.75};b.level=1;mesh->applyRefineBoxes({b});}return mesh;
}
void statics(){
    std::ofstream geometries("output/series/drop_geometries.csv"),samples("output/series/drop_surface.csv"),faces("output/series/drop_faces.csv"),cells("output/series/drop_cells.csv"),summary("output/series/drop_solves.csv"),fields("output/series/drop_fields.csv"),velocities("output/series/drop_velocity.csv");
    for(auto* f:{&geometries,&samples,&faces,&cells,&summary,&fields,&velocities})*f<<std::setprecision(17);
    geometries<<"geometry,n,rz,amr,nonpolynomial,center,radius,curvature_l2,curvature_linf,area,liquid_volume,cells,faces\n";
    samples<<"geometry,component,graph,node,cell,x,y,nx,ny,curvature,weight\n";
    cells<<"geometry,p,native,x0,y0,x1,y1,volume,liquid,alpha\n";
    faces<<"geometry,u,minus,plus,horizontal,x0,y0,x1,y1,area,surface_measure,curvature,force\n";
    summary<<"solve,geometry,rho_l,mu_l,analytic,velocity_max,velocity_load_relative,velocity_physical_relative,pressure_error,algebraic,momentum,continuity,pressure_mean,linear_refinements,seconds\n";
    fields<<"solve,p,pressure\n";velocities<<"solve,u,velocity,mass,diagonal,force,pressure_load,inertia,reaction\n";
    int gid=0,sid=0;
    for(bool rz:{false,true})for(bool amr:{false,true})for(int n:{8,16,32})for(double center:{.4,.44,.48})for(bool nonpolynomial:{false,true}){
        ++gid;auto mesh=dropMesh(n,rz,amr);Case shape;shape.pureFirstJet=true;shape.radial=rz;shape.nonpolynomial=nonpolynomial;shape.cx=center;shape.cy=rz?center:.5;
        const auto metrics=metric(*mesh,rz);const auto phase=CchPhaseGeometry::build(*mesh,metrics,initial(*mesh,shape));const auto& l=phase.layout();
        const auto pressure=CchPressureLoad::build(metrics.divergence());const auto geometry=CchMomentumGeometry::build(*mesh,metrics);
        const auto reconstruction=CchAdvectorReconstruction::prepare(*mesh,geometry);const auto law=closedBoundary(*mesh,geometry);
        const auto actual=CchInterfaceNormalLoad::capillary(*mesh,pressure,phase,1,reconstruction,&law);const double exact=(rz?2:1)/shape.radius;
        const auto analytic=CchInterfaceNormalLoad::constant(*mesh,pressure,phase,-exact);const auto& pub=phase.publication();
        double area=0,error=0,linf=0,liquid=0,volume=0;
        for(std::size_t c=0;c<pub.topology().curves.size();++c)for(std::size_t g=0;g<pub.topology().curves[c].graphs().size();++g){int node=0;for(const auto& q:pub.quadrature(c,g)){
            const double d=q.sample.curvature-exact;area+=q.physicalWeight;error+=q.physicalWeight*d*d;linf=std::max(linf,std::abs(d));const auto& f=q.sample.frame;
            samples<<gid<<','<<c<<','<<g<<','<<node++<<','<<f.nativeCellIndex<<','<<f.point.x<<','<<f.point.y<<','<<f.normal.x<<','<<f.normal.y<<','<<q.sample.curvature<<','<<q.physicalWeight<<'\n';}}
        for(int p=0;p<l.pressureCount();++p){const auto& cell=mesh->cells(mesh->cellIndexToLabel[l.nativeCell(p)]);const auto a=mesh->verts(cell.cellToVert[0]).coord(),b=mesh->verts(cell.cellToVert[2]).coord();liquid+=phase.cellVolumes()[p].liquid;volume+=metrics.cellVolumes()[p];
            cells<<gid<<','<<p<<','<<l.nativeCell(p)<<','<<a.x<<','<<a.y<<','<<b.x<<','<<b.y<<','<<metrics.cellVolumes()[p]<<','<<phase.cellVolumes()[p].liquid<<','<<phase.liquidFractions()[p]<<'\n';}
        for(int u=0;u<l.velocityCount();++u){const auto f=l.domain().face(*mesh,l.nativeFace(u));const auto& edge=f.horizontal?mesh->hEdges(f.edgeLabel):mesh->vEdges(f.edgeLabel);const auto a=mesh->verts(edge.edgeToVert[0]).coord(),b=mesh->verts(edge.edgeToVert[1]).coord();
            faces<<gid<<','<<u<<','<<(f.negativeCell<0?-1:l.pressureDof(f.negativeCell))<<','<<(f.positiveCell<0?-1:l.pressureDof(f.positiveCell))<<','<<f.horizontal<<','<<a.x<<','<<a.y<<','<<b.x<<','<<b.y<<','<<metrics.faceMeasures()[u]<<','<<actual.supportMeasures()[u]<<','<<actual.curvatureAverages()[u]<<','<<actual.balanceLoads()[u]<<'\n';}
        geometries<<gid<<','<<n<<','<<rz<<','<<amr<<','<<nonpolynomial<<','<<center<<','<<shape.radius<<','<<std::sqrt(error/area)/exact<<','<<linf/exact<<','<<area<<','<<liquid<<','<<l.pressureCount()<<','<<l.velocityCount()<<'\n';geometries.flush();
        const auto quadrature=HermiteGalsPhaseQuadrature::build(*mesh,pub);std::vector<CchPhaseViscosity> viscosities;
        for(double mu:{1.,1e-3,1e-6,0.})viscosities.push_back(CchPhaseViscosity::build(*mesh,phase,quadrature,mu,.01*mu));
        for(double rho:{1.,1e3,1e6}){const auto mass=CchMomentumMass::build(*mesh,geometry,phase,rho,1);
            for(const auto& viscosity:viscosities){Matrix block=viscosity.viscosity().matrix();for(int u=0;u<l.velocityCount();++u)block.coeffRef(u,u)+=mass.masses()[u]/.125;
                const auto saddle=CchSinglePhaseSaddle::buildClosed(pressure,CchSaddleVelocityBlock::build(l,block));
                for(bool useAnalytic:{true,false}){++sid;const auto started=Clock::now();const auto& force=(useAnalytic?analytic:actual).balanceLoads();auto rhs=force;for(double& value:rhs)value=-value;
                    std::cout<<"DROP_BEGIN,"<<sid<<','<<gid<<','<<rho<<','<<viscosity.liquidViscosity()<<','<<useAnalytic<<std::endl;
                    const auto result=saddle.solve(*mesh,CchSaddleRhs(l,rhs,std::vector<double>(l.velocityCount())));const auto& u=result.values().normalVelocity();const auto& p=result.values().pressure();
                    const auto cp=pressure.apply(*mesh,result.values());const Eigen::Map<const Eigen::VectorXd> ev(u.data(),u.size());const Eigen::VectorXd inertia=block*ev;
                    double maxVelocity=0,loadScale=0,maxPressure=0;
                    for(int k=0;k<l.pressureCount();++k){maxPressure=std::max(maxPressure,std::abs(p[k]-exact*(phase.liquidFractions()[k]-liquid/volume)));fields<<sid<<','<<k<<','<<p[k]<<'\n';}
                    for(int k=0;k<l.velocityCount();++k){maxVelocity=std::max(maxVelocity,std::abs(u[k]));loadScale=std::max(loadScale,std::abs(force[k])/block.coeff(k,k));
                        velocities<<sid<<','<<k<<','<<u[k]<<','<<mass.masses()[k]<<','<<block.coeff(k,k)<<','<<force[k]<<','<<cp[k]<<','<<inertia[k]<<','<<inertia[k]+cp[k]+force[k]<<'\n';}
                    const double normalized=loadScale>0?maxVelocity/loadScale:0;const auto& report=result.report();
                    summary<<sid<<','<<gid<<','<<rho<<','<<viscosity.liquidViscosity()<<','<<useAnalytic<<','<<maxVelocity<<','<<normalized<<','<<maxVelocity/(.125/(shape.radius*shape.radius))<<','<<maxPressure/exact<<','<<report.algebraicRelativeResidual<<','<<report.freeMomentumDefect<<','<<report.maxCellContinuityDefect<<','<<report.maxPressureMeanDefect<<','<<report.linearRefinements<<','<<std::chrono::duration<double>(Clock::now()-started).count()<<'\n';summary.flush();
                    if(useAnalytic)check(normalized<=1e-11&&maxPressure/exact<=1e-10,"original H12 analytic static gate");
                    std::cout<<"DROP_RESULT,"<<sid<<','<<normalized<<','<<maxPressure/exact<<std::endl;
                }
            }
        }
        for(auto* f:{&samples,&faces,&cells,&fields,&velocities}){f->flush();check(bool(*f),"static export write failure");}
    }
    check(gid==72&&sid==1728,"incomplete static parameter product");std::cout<<"DROP_DATA_COMPLETE,"<<gid<<','<<sid<<std::endl;
}
// Read-only load/kinematic range audit. The dense basis enumeration is a
// research reference, never a replacement reconstruction or production force.
void workWitness(int n){
    auto mesh=makeMesh(n);const auto pub=waveInitial(*mesh,1e-4);const auto metrics=metric(*mesh);
    const auto phase=CchPhaseGeometry::build(*mesh,metrics,pub);const auto geometry=CchMomentumGeometry::build(*mesh,metrics);
    const auto pressure=CchPressureLoad::build(metrics.divergence());
    const auto mass=CchMomentumMass::build(*mesh,geometry,phase,1,1);const auto reconstruction=CchAdvectorReconstruction::prepare(*mesh,geometry);
    const auto& l=geometry.layout();const int nu=l.velocityCount(),np=l.pressureCount();std::vector<int> walls;
    for(int u=0;u<nu;++u){const auto f=l.domain().face(*mesh,l.nativeFace(u));if(f.negativeCell<0||f.positiveCell<0)walls.push_back(u);}
    const auto law=CchAdvectorBoundary::impermeable(*mesh,geometry,walls);
    const auto cap=CchInterfaceNormalLoad::capillary(*mesh,pressure,phase,1,reconstruction,&law);
    std::vector<HermiteGalsSurfaceQuadrature> quadrature;
    for(std::size_t c=0;c<pub.topology().curves.size();++c)for(std::size_t g=0;g<pub.topology().curves[c].graphs().size();++g)for(const auto& q:pub.quadrature(c,g))quadrature.push_back(q);
    const auto incidence=pressure.apply(*mesh,CchNativeValues(l,std::vector<double>(nu),phase.liquidFractions()));
    std::vector<double> surfaceGradient(nu),volumeGradient(nu),basis(nu);double curvatureMean=0,area=0;
    for(const auto& q:quadrature){curvatureMean+=q.physicalWeight*q.sample.curvature;area+=q.physicalWeight;}curvatureMean/=area;
    const std::string suffix=std::to_string(n)+".csv";
    std::ofstream nodes("output/series/work_nodes_"+suffix),response("output/series/work_response_"+suffix),loads("output/series/work_loads_"+suffix),cycles("output/series/work_cycles_"+suffix),cells("output/series/work_cells_"+suffix);
    for(auto* f:{&nodes,&response,&loads,&cycles,&cells})*f<<std::setprecision(17);
    nodes<<"q,cell,x,y,nx,ny,weight,curvature\n";response<<"q,u,normal_velocity\n";
    loads<<"u,minus,plus,horizontal,x0,y0,x1,y1,area,mass,C_alpha,balance_load,surface_gradient,volume_gradient\n";
    cycles<<"i,j,x,y,incidence_motion,divergence,surface_power,mechanical_power,volume_rate\n";
    cells<<"p,native,volume,alpha\n";
    for(int p=0;p<np;++p)cells<<p<<','<<l.nativeCell(p)<<','<<metrics.cellVolumes()[p]<<','<<phase.liquidFractions()[p]<<'\n';
    for(std::size_t q=0;q<quadrature.size();++q){const auto& s=quadrature[q];const auto& f=s.sample.frame;nodes<<q<<','<<f.nativeCellIndex<<','<<f.point.x<<','<<f.point.y<<','<<f.normal.x<<','<<f.normal.y<<','<<s.physicalWeight<<','<<s.sample.curvature<<'\n';}
    for(int u=0;u<nu;++u){const auto f=l.domain().face(*mesh,l.nativeFace(u));if(f.negativeCell<0||f.positiveCell<0)continue;
        basis[u]=1;const auto field=reconstruction.evaluate(*mesh,CchNativeValues(l,basis,std::vector<double>(np)),phase.geometryIdentity(),phase.timeLayer(),law);basis[u]=0;
        for(std::size_t q=0;q<quadrature.size();++q){const auto& s=quadrature[q];const auto& frame=s.sample.frame;const auto velocity=field.cellVelocity(l.pressureDof(frame.nativeCellIndex),frame.point);
            const double normal=velocity.x*frame.normal.x+velocity.y*frame.normal.y;surfaceGradient[u]+=s.physicalWeight*s.sample.curvature*normal;volumeGradient[u]+=s.physicalWeight*normal;
            response<<q<<','<<u<<','<<normal<<'\n';}}
    for(int u=0;u<nu;++u){const auto f=l.domain().face(*mesh,l.nativeFace(u));const auto& edge=f.horizontal?mesh->hEdges(f.edgeLabel):mesh->vEdges(f.edgeLabel);const auto a=mesh->verts(edge.edgeToVert[0]).coord(),b=mesh->verts(edge.edgeToVert[1]).coord();
        loads<<u<<','<<(f.negativeCell<0?-1:l.pressureDof(f.negativeCell))<<','<<(f.positiveCell<0?-1:l.pressureDof(f.positiveCell))<<','<<f.horizontal<<','<<a.x<<','<<a.y<<','<<b.x<<','<<b.y<<','<<metrics.faceMeasures()[u]<<','<<mass.masses()[u]<<','<<incidence[u]<<','<<cap.balanceLoads()[u]<<','<<surfaceGradient[u]<<','<<volumeGradient[u]<<'\n';}
    double invisiblePeak=0;int invisibleCount=0;
    for(int i=1;i<n;++i)for(int j=1;j<n;++j){const double x=double(i)/n,y=double(j)/n;double seen=0,power=0,mechanical=0,volume=0;
        const auto psi=[&](Vector2D p){return p.x==x&&p.y==y?1.:0.;};std::vector<double> velocity(nu);
        for(int u=0;u<nu;++u){const auto f=l.domain().face(*mesh,l.nativeFace(u));const auto& edge=f.horizontal?mesh->hEdges(f.edgeLabel):mesh->vEdges(f.edgeLabel);auto a=mesh->verts(edge.edgeToVert[0]).coord(),b=mesh->verts(edge.edgeToVert[1]).coord();if(f.horizontal?a.x>b.x:a.y>b.y)std::swap(a,b);
            velocity[u]=f.horizontal?-(psi(b)-psi(a))/(b.x-a.x):(psi(b)-psi(a))/(b.y-a.y);
            if(incidence[u]!=0)seen=std::max(seen,std::abs(velocity[u]));power+=surfaceGradient[u]*velocity[u];mechanical-=cap.balanceLoads()[u]*velocity[u];volume+=volumeGradient[u]*velocity[u];}
        const auto div=metrics.apply(*mesh,CchNativeValues(l,velocity,std::vector<double>(np)));double maxDiv=0;for(double v:div)maxDiv=std::max(maxDiv,std::abs(v));check(maxDiv==0,"native streamfunction cycle is not exactly divergence free");
        if(seen==0){++invisibleCount;invisiblePeak=std::max(invisiblePeak,std::abs(power));}
        cycles<<i<<','<<j<<','<<x<<','<<y<<','<<seen<<','<<maxDiv<<','<<power<<','<<mechanical<<','<<volume<<'\n';}
    Matrix block(nu,nu);for(int u=0;u<nu;++u)block.coeffRef(u,u)=mass.masses()[u];const auto saddle=CchSinglePhaseSaddle::buildClosed(pressure,CchSaddleVelocityBlock::build(l,block));
    for(bool workDual:{false,true}){std::vector<double> rhs(nu);for(int u=0;u<nu;++u)rhs[u]=workDual?curvatureMean*incidence[u]-(surfaceGradient[u]-curvatureMean*volumeGradient[u]):-cap.balanceLoads()[u];
        const auto result=saddle.solve(*mesh,CchSaddleRhs(l,rhs,std::vector<double>(nu)));const auto field=reconstruction.evaluate(*mesh,result.values(),phase.geometryIdentity(),phase.timeLayer(),law);
        double modal=0,normalizer=0,surfacePower=0,mechanicalPower=0;const double k=double(2*pi);
        for(const auto& q:quadrature){const auto& f=q.sample.frame;const auto v=field.cellVelocity(l.pressureDof(f.nativeCellIndex),f.point);const double normal=v.x*f.normal.x+v.y*f.normal.y,mode=std::cos(k*f.point.x);
            modal+=q.physicalWeight*mode*normal;normalizer+=q.physicalWeight*f.normal.y*mode*mode;surfacePower+=q.physicalWeight*q.sample.curvature*normal;}
        for(int u=0;u<nu;++u)mechanicalPower+=rhs[u]*result.values().normalVelocity()[u];
        std::cout<<"WORK_ACCELERATION,"<<n<<','<<workDual<<','<<-modal/(1e-4*normalizer)<<','<<surfacePower<<','<<mechanicalPower<<','<<surfacePower+mechanicalPower<<std::endl;}
    for(auto* f:{&nodes,&response,&loads,&cycles,&cells}){f->flush();check(bool(*f),"work witness export failure");}
    std::cout<<"WORK_RANGE,"<<n<<','<<invisibleCount<<','<<invisiblePeak<<','<<curvatureMean<<std::endl;
}
}
int main(int argc,char** argv){try{
    const std::string mode=argc>1?argv[1]:"probe";check(argc<=5&&(mode=="probe"||mode=="wave"||mode=="static"||mode=="work")&&(argc<3||mode=="probe"||(mode=="work"&&argc==3)),"usage: capillary benchmark [probe [n [steps [periods]]]|work [n]|wave|static]");
    const int probeResolution=argc>=3?std::stoi(argv[2]):8;check(probeResolution>=4&&probeResolution<=64,"probe resolution range");
    const int probeSteps=argc>=4?std::stoi(argv[3]):probeResolution;const double probePeriods=argc>=5?std::stod(argv[4]):1./16;
    check(probeSteps>0&&std::isfinite(probePeriods)&&probePeriods>0,"invalid probe duration/steps");
    std::filesystem::create_directories("output/series");std::cout<<std::setprecision(17);
    std::ofstream history("output/series/wave_history.csv"),summary("output/series/wave_summary.csv");history<<std::setprecision(17);summary<<std::setprecision(17);
    history<<"n,amr,steps,periods,amplitude,mu,step,time,mode,mean,length,kinetic,capillary_work,viscous_work,energy_defect,curvature_l2,R_primal_l1,R_dual_l1,R_kinetic_power,continuity\n";
    summary<<"n,amr,steps,periods,amplitude,mu,reference_omega,measured_omega,final_mode,peak_energy,peak_continuity,peak_R,peak_R_power,seconds\n";
    if(mode=="static")statics();
    else if(mode=="work")workWitness(probeResolution);
    else {WaveOutput raw;
        if(mode=="probe")wave(probeResolution,false,probeSteps,probePeriods,1e-4,0,history,summary,raw);
        else for(int n:{8,16,32})for(bool amr:{false,true})wave(n,amr,32*n,1,1e-4,0,history,summary,raw);
    }
    return 0;
}catch(const std::exception& e){std::cerr<<"H13 RESEARCH FAIL: "<<e.what()<<'\n';return 1;}}
