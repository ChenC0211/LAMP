#include "Modules/Electricity/DielectricEvolution.h"
#include "Modules/Electricity/IO/DictionaryDielectric.h"
#include "Core/Geometry/InterfaceGeometry.h"
#include "Core/Field/StaggeredVelocitySampler.h"
#include "Core/Field/VEdgeField.h"
#include "Core/Field/HEdgeField.h"
#include "Core/IO/DictionaryMesh.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <filesystem>
#include <fstream>
#include <chrono>

namespace {
using namespace Quantum::Electricity;
void checkContract(bool value, const char* message)
{ if (!value) throw std::runtime_error(message); }
template<class F> void rejects(F action)
{
    try { action(); } catch (const std::logic_error&) { return; }
    throw std::runtime_error("invalid electric lifecycle input accepted");
}
void equalField(const std::vector<double>& a, const std::vector<double>& b)
{
    checkContract(a.size() == b.size(), "field size differs");
    for (int i=0; i<a.size(); ++i)
        checkContract(std::isfinite(a[i]) && a[i] == b[i], "electric field changed from original solver");
}
const DielectricFields& fields(const DielectricEvolution& service)
{ return *service.evaluated()->electric.fields; }
static_assert(!std::is_assignable_v<decltype(std::declval<const DielectricFields&>().potential[0]),double>);
static_assert(!std::is_assignable_v<decltype(std::declval<const DielectricChargeSnapshot&>().integratedCharge[0]),double>);
static_assert(std::is_same_v<decltype(std::declval<DielectricSolveResult>().fields),std::shared_ptr<const DielectricFields>>);

void domainInputContract()
{
    namespace fs=std::filesystem;
    const auto cwd=fs::current_path();
    const auto directory=cwd/("electric_domain_input_"+std::to_string(std::chrono::high_resolution_clock::now().time_since_epoch().count()));
    fs::create_directories(directory);
    const auto electric=directory/"physical.dict", geometry=directory/"shape.dict", mesh=directory/"grid.dict";
    const auto write=[](const fs::path& path,const std::string& text) {
        std::ofstream stream(path);stream<<text;checkContract(bool(stream),"cannot write electric input fixture");
    };
    write(geometry,R"(Domain
{
 min (0.2, -0.5)
 max (1.2, 0.5)
 dimension 2
}
Point
{
 (0.4, -0.2)
 (0.8, -0.2)
 (0.8, 0.2)
 (0.4, 0.2)
}
Line
{
 (1, 2)
 {
  type straight
 }
 (2, 3)
 {
  type straight
 }
 (3, 4)
 {
  type straight
 }
 (4, 1)
 {
  type straight
 }
}
LevelSet
{
 "liquid"
 {
  lines (1, 2, 3, 4)
  innerPoint (0.6, 0)
  moving false
 }
}
)");
    write(mesh,"adaptiveMesh false\nBackgroundMesh\n{\n divide (8, 10)\n}\n");
    const auto input=[](const char* method,const std::string& binding) {
        return std::string("Electricity { model perfectDielectric; interfaceDiscretization ")+method+
            "; permittivityNegative 3; permittivityPositive 1; electrodeAxis y; lowerPotential 0; upperPotential 0.01; }\n"
            "Domain { coordinateSystem cartesianXY; interfaceLevelSet liquid; "+binding+" }\n";
    };
    const auto read=[&]{return readDielectricDomainInput(.001,electric.string(),geometry.string(),mesh.string());};
    write(electric,input("diffuseTanh","interfaceHalfWidth 0.07; motion fixed;"));
    const auto diffuse=read();
    checkContract(diffuse.electric.interfaceDiscretization=="diffuseTanh" && diffuse.electric.permittivityNegative==3.
        && diffuse.diffuseThickness==.07 && !diffuse.ignoredSharpHalfWidth && !diffuse.electric.axisymmetric,
        "module input lost actual material/domain values");
    checkContract(diffuse.geometry.hasDomain && diffuse.geometry.levelSetName.size()==1
        && diffuse.geometry.levelSetName[0]=="liquid" && !diffuse.geometry.moving[0]
        && diffuse.mesh.Nx==8 && diffuse.mesh.Ny==10 && diffuse.mesh.minCoord.x==.2
        && diffuse.mesh.minCoord.y==-.5 && diffuse.mesh.maxCoord.x==1.2 && diffuse.mesh.maxCoord.y==.5,
        "explicit mesh reader did not consume the same named geometry file");
    write(electric,input("sharpHermiteP3","interfaceHalfWidth 0.08; motion fixed;"));
    const auto compatibility=read();
    checkContract(compatibility.diffuseThickness==0. && compatibility.ignoredSharpHalfWidth==.08,
        "sharp compatibility width entered physical material smoothing");
    write(electric,input("sharpHermiteP3","motion fixed;"));
    const auto sharp=read();
    checkContract(sharp.diffuseThickness==0. && !sharp.ignoredSharpHalfWidth,"module input retained another read's half-width");
    const auto rejected=[&](const std::string& text,const char* expected) {
        write(electric,text);bool caught=false;
        try {read();}catch(const std::exception& error) {caught=std::string(error.what()).find(expected)!=std::string::npos;}
        checkContract(caught,"physical-domain rejection changed");
    };
    rejected(input("sharpHermiteP3","motion moving;"),"fixed");
    rejected(input("diffuseTanh","motion fixed;"),"interfaceHalfWidth");
    rejected(input("sharpHermiteP3","interfaceHalfWidth -1; motion fixed;"),"positive finite legacy half-width");
    write(electric,input("sharpHermiteP3","motion fixed;"));
    write(mesh,"adaptiveMesh false\nBackgroundMesh\n{\n divide (4, 10)\n}\n");
    bool small=false;try {read();}catch(const std::invalid_argument& error) {
        small=std::string(error.what()).find("Nx/Ny >= 8")!=std::string::npos;
    }
    checkContract(small && fs::current_path()==cwd && !fs::exists(directory/"output") && !fs::exists(directory/"cache"),
        "domain reader changed cwd, created output, or accepted an unsupported grid");

    const std::string common="permittivityNegative 2; permittivityPositive 1; electrodeAxis y; "
        "lowerPotential 0; upperPotential 0.01; linearResidualTolerance 2e-5; ";
    const auto physical=[&](bool charging,const std::string& extra) {
        return "Electricity { "+common+(charging
            ? "model leakyDielectric; interfaceDiscretization sharpHermiteP3; conductivityNegative .1; "
              "conductivityPositive .01; initialSurfaceChargeDensity .001; chargeConvection false; "
            : "model perfectDielectric; interfaceDiscretization diffuseTanh; ")+extra+" }";
    };
    write(electric,physical(false,""));
    const auto legacyDiffuse=readDielectricConfig(electric.string());
    const auto formalDiffuse=readModularDielectricConfig(electric.string());
    checkContract(legacyDiffuse.boundaryFluxTolerance==2e-5 && !formalDiffuse.boundaryFluxTolerance,
        "formal diffuse conservation inherited linear precision");
    write(electric,physical(true,""));
    const auto legacySharp=readDielectricConfig(electric.string());
    const auto formalSharp=readModularDielectricConfig(electric.string());
    checkContract(legacySharp.sharpFluxChecks->reconstructedPoisson==2e-5
        && legacySharp.sharpFluxChecks->gaussClosure==2e-5 && legacySharp.chargeClosureTolerance==16.*2e-5
        && !formalSharp.sharpFluxChecks && !formalSharp.chargeClosureTolerance,
        "old uniform expansion or independent formal defaults changed");
    write(electric,physical(true,"Conservation { reconstructedPoisson 3e-9; gaussClosure 4e-9; chargeClosure 5e-9; }"));
    const auto configured=readModularDielectricConfig(electric.string());
    checkContract(configured.linearResidualTolerance==2e-5 && configured.sharpFluxChecks->reconstructedPoisson==3e-9
        && configured.sharpFluxChecks->gaussClosure==4e-9 && configured.chargeClosureTolerance==5e-9,
        "independent sharp/charge settings lost their consumers");
    bool legacyRejected=false;
    try {readDielectricConfig(electric.string());}catch(const std::exception& error) {
        legacyRejected=std::string(error.what()).find("requires formal")!=std::string::npos;
    }
    checkContract(legacyRejected,"legacy input silently accepted formal Conservation");
    write(electric,physical(false,"Conservation { boundaryFlux 6e-9; }"));
    checkContract(readModularDielectricConfig(electric.string()).boundaryFluxTolerance==6e-9,
        "diffuse physical flux setting was ignored");
    const auto badSetting=[&](bool charging,const std::string& extra,const char* expected) {
        write(electric,physical(charging,extra));bool rejected=false;
        try {readModularDielectricConfig(electric.string());}catch(const std::exception& error) {
            rejected=std::string(error.what()).find(expected)!=std::string::npos;
        }
        checkContract(rejected,"invalid/unconsumed electric Conservation setting accepted");
    };
    badSetting(true,"Conservation { boundaryFlux 1e-9; }","only consumed by diffuseTanh");
    badSetting(false,"Conservation { gaussClosure 1e-9; }","not consumed by diffuse");
    badSetting(false,"Conservation { chargeClosure 1e-9; }","requires leakyDielectric");
    badSetting(true,"Conservation { chargeClosure 0; }","finite and in (0,1)");
    badSetting(true,"Conservation { reconstructedPoisson 1; }","finite and in (0,1)");
    badSetting(false,"Conservation { boundaryFlux nan; }","finite number");
    badSetting(true,"Conservation { chargeClosure 1e-9; chargeClosure 2e-9; }","duplicate");
    badSetting(true,"Conservation { unknown 1e-9; }","unknown");
}

void setCircle(const QuadtreeMesh& mesh, VertField<double>& scalar, double x, double y)
{
    for (int i=0; i<mesh.vertIndexNum(); ++i) {
        const auto p=mesh.verts(mesh.vertIndexToLabel[i]).coord();
        scalar.value[i]=std::hypot(p.x-x,p.y-y)-.23;
    }
}
DielectricConfig config()
{
    DielectricConfig c;
    c.model="leakyDielectric"; c.interfaceDiscretization="sharpHermiteP3";
    c.conductivityNegative=.1; c.conductivityPositive=.01;
    c.initialSurfaceChargeDensity=.001; c.upperPotential=.01; c.linearResidualTolerance=1e-8;
    // Keep the historical solve/transport reference's explicit physical gates.
    c.sharpFluxChecks=DielectricFluxChecks(c.linearResidualTolerance);
    c.chargeClosureTolerance=16.*c.linearResidualTolerance;
    return c;
}
void conservationContract(QuadtreeMesh& mesh)
{
    auto c=config();c.linearResidualTolerance=2e-5;
    c.sharpFluxChecks.reset();c.chargeClosureTolerance.reset();
    VertField<double> phi(mesh,0.);setCircle(mesh,phi,.5,.5);
    InterfaceGeometry geometry(mesh,phi);
    SurfaceChargeField charge(mesh,geometry.snapshot(),c.initialSurfaceChargeDensity);
    const auto source=charge.source();
    DielectricElectrostatics solver(mesh);
    const auto solve=[&](const DielectricConfig& settings) {
        return solver.solve(mesh,phi,0.,settings,&geometry,nullptr,&source);
    };
    const auto baseline=solve(c);
    checkContract(baseline.report.success && baseline.report.usedSuperLU && baseline.report.relativeResidual>0.,
        "conservation contract requires a real sharp solve with nonzero reconstructed residual");
    auto tight=c;tight.sharpFluxChecks.emplace();
    tight.sharpFluxChecks->reconstructedPoisson=.5*baseline.report.relativeResidual;
    const auto physicalFailure=solve(tight);
    checkContract(!physicalFailure.report.success && !physicalFailure.fields && physicalFailure.report.usedSuperLU
        && physicalFailure.report.reason.find("reconstructedPoisson tolerance")!=std::string::npos
        && physicalFailure.report.relativeResidual==baseline.report.relativeResidual,
        "physical reconstruction tolerance did not independently reject the same solved field");
    tight=c;tight.linearResidualTolerance=1e-30;
    const auto linearFailure=solve(tight);
    checkContract(!linearFailure.report.success && !linearFailure.fields
        && linearFailure.report.stage=="sharp-solve", "linear failure was replaced by physical acceptance");
    DielectricEvolution ordinary(mesh,c,0.);
    const auto evaluated=ordinary.evaluate(mesh,phi,{1,0.,.001},&geometry);
    const auto completed=ordinary.completeCharge({1,0.,.001});
    checkContract(completed.transport.transport.success && completed.charge->completedStep==1,
        "default independent charge closure failed");
    tight=c;tight.chargeClosureTolerance=1e-30;
    DielectricEvolution strict(mesh,tight,0.);
    const auto retained=strict.evaluate(mesh,phi,{1,0.,.001},&geometry);
    equalField(retained->electric.fields->potential,evaluated->electric.fields->potential);
    bool rejected=false;
    try {strict.completeCharge({1,0.,.001});}catch(const std::runtime_error& error) {
        rejected=std::string(error.what()).find("non-interface conduction current")!=std::string::npos
            || std::string(error.what()).find("close with electrode")!=std::string::npos;
    }
    checkContract(rejected && strict.failed() && retained->charge->completedStep==0,
        "chargeClosure did not independently reject physical current/charge imbalance");
    rejects([&]{strict.completeCharge({1,0.,.001});});
}
std::shared_ptr<const InterfaceTraceGeometry> traceQueries(const InterfaceGeometry& geometry)
{
    std::vector<InterfaceTraceGeometry::Point> points;
    const auto& snapshot=geometry.snapshot();
    for(int i=0;i<8;++i) {
        const double angle=i*(2.*std::acos(-1.))/8.;
        const auto frame=snapshot.projectedFrame({.5+.25*std::cos(angle),.5+.25*std::sin(angle)});
        points.push_back({frame,snapshot.reconstructedCurvatureTangentialDerivative(frame.point),0.});
    }
    return std::make_shared<const InterfaceTraceGeometry>(geometry.source(),std::move(points));
}
void runTransport(QuadtreeMesh& mesh, int mode)
{
    auto c=config(); c.chargeConvection=mode!=0;
    DielectricEvolution service(mesh,c,.08);
    DielectricElectrostatics original(mesh);
    DielectricElectrostatics withoutTraces(mesh);
    VertField<double> scalar(mesh,0.);
    setCircle(mesh,scalar,.5,.5);
    std::unique_ptr<SurfaceChargeField> reference;
    rejects([&] { fields(service); });
    VEdgeField<double> ux(mesh,.02);
    HEdgeField<double> uy(mesh,-.01);
    StaggeredVelocitySampler velocity(mesh,ux,uy,ux,uy);
    constexpr double dt=.001;
    for(int step=1; step<=3; ++step) {
        const DielectricStepTime time{step,(step-1)*dt,dt};
        InterfaceGeometry load(mesh,scalar);
        const auto queries=traceQueries(load);
        rejects([&] { service.evaluate(mesh,scalar,{step+1,time.forceTime,dt},&load); });
        checkContract(!service.failed(),"bad preflight poisoned electric state");
        auto evaluated=service.evaluate(mesh,scalar,time,&load,nullptr,queries);
        const auto& report=evaluated->electric.report;
        const auto oldCharge=evaluated->charge;
        const auto oldIntegrals=oldCharge->integratedCharge;
        const auto oldPotential=evaluated->electric.fields->potential;
        const auto oldScalar=evaluated->electric.fields->levelSet;
        checkContract(oldCharge==service.charge() && oldCharge->completedStep==step-1
            && oldCharge->time==time.forceTime,"field publication lost its actual charge time");
        checkContract(report.usedSuperLU && report.usedSharpInterface && report.usedSurfaceCharge,
            "production sharp/SuperLU path was not used");
        rejects([&] { service.evaluate(mesh,scalar,{step+1,step*dt,dt},&load); });
        rejects([&] { service.completeCharge({step,time.forceTime+dt,dt}); });
        if(!reference) reference=std::make_unique<SurfaceChargeField>(mesh,load.snapshot(),c.initialSurfaceChargeDensity);
        auto source=reference->source();
        DielectricConductivity conductivity{c.conductivityNegative,c.conductivityPositive};
        const auto oldSolve=original.solve(mesh,scalar,.08,c,&load,nullptr,&source,&conductivity,queries);
        const auto& oldReport=oldSolve.report;
        checkContract(oldReport.success && oldSolve.fields && oldReport.conductionCurrent,"reference electric solve failed");
        const auto& originalFields=*oldSolve.fields;
        const auto fieldOnly=withoutTraces.solve(mesh,scalar,0.,c,&load,nullptr,&source,&conductivity);
        checkContract(fieldOnly.report.success && fieldOnly.report.localTraceFits==0 && !fieldOnly.fields->traceGeometry
            && fieldOnly.fields->interfaceTraces.empty(),"field-only sharp solve created load sampling");
        equalField(fieldOnly.fields->potential,originalFields.potential);
        checkContract(fieldOnly.report.conductionCurrent->cellIntegratedCurrent==oldReport.conductionCurrent->cellIntegratedCurrent,
            "omitting traction queries changed actual conduction currents");
        equalField(fields(service).potential,originalFields.potential);
        equalField(fields(service).ex,originalFields.ex);
        equalField(fields(service).ey,originalFields.ey);
        checkContract(fields(service).traceGeometry==queries, "sharp trace source was lost");
        const auto actual=fields(service).interfaceTractions(), expected=originalFields.interfaceTractions();
        checkContract(actual.size()==expected.size(), "traction count changed");
        for(std::size_t i=0;i<actual.size();++i)
            checkContract(actual[i].x==expected[i].x && actual[i].y==expected[i].y,"physical traction changed");
        const auto& current=*oldReport.conductionCurrent;
        SurfaceChargeFieldStep oldUpdate;
        DielectricChargeUpdate update;
        if(mode==2) {
            setCircle(mesh,scalar,.5+step*dt*.02,.5-step*dt*.01);
            LevelSetHermiteGeometry::Snapshot next(mesh,scalar);
            const std::function<Vector2D(Vector2D)> forward=[dt](Vector2D p) { return p+Vector2D(dt*.02,-dt*.01); };
            const DielectricChargeMotion motion{mesh,next,forward};
            update=service.completeCharge(time,nullptr,&motion);
            reference=reference->advanceMoving(mesh,next,forward,current.cellIntegratedCurrent,
                current.absoluteFaceCurrent,16.0*c.linearResidualTolerance,dt,oldUpdate);
        } else if(mode==1) {
            update=service.completeCharge(time,&velocity);
            oldUpdate=reference->advance(current.cellIntegratedCurrent,reference->advectiveMeasureFlux(velocity),
                current.absoluteFaceCurrent,16.0*c.linearResidualTolerance,dt);
        } else {
            update=service.completeCharge(time);
            oldUpdate=reference->advanceNormalCurrent(current.cellIntegratedCurrent,current.absoluteFaceCurrent,
                16.0*c.linearResidualTolerance,dt);
        }
        checkContract(update.current==report.conductionCurrent,"charge did not consume the actual staged current");
        checkContract(service.evaluated()==evaluated && update.charge==service.charge()
            && update.charge!=oldCharge && update.charge->completedStep==step && update.charge->time==step*dt,
            "charge completion replaced the old field or published the wrong charge time");
        checkContract(oldCharge->integratedCharge==oldIntegrals && evaluated->electric.fields->potential==oldPotential,
            "advancing charge mutated a retained field/charge publication");
        checkContract(evaluated->electric.fields->levelSet==oldScalar,"material motion changed the retained source scalar");
        if(mode!=2)checkContract(update.charge->geometry==oldCharge->geometry,"fixed charge step rebuilt immutable geometry");
        checkContract(oldUpdate.transport.success && update.transport.transport.success,"charge transport failed");
        checkContract(service.charge()->integratedCharge==reference->transport().panelCharge(),
            "charge integrals changed from original transport");
        const auto& q=oldUpdate.transport;
        checkContract(update.electrodeResidual==q.totalChargeAfter-q.totalChargeBefore+dt*current.netBoundaryCurrent,
            "electrode accounting changed");
        equalField(fields(service).potential,originalFields.potential); // Still the t^n field, no hidden solve.
        rejects([&] { service.completeCharge(time); });
    }
}
void failureContract(QuadtreeMesh& mesh)
{
    auto c=config(); c.chargeConvection=true;
    DielectricEvolution service(mesh,c,.08);
    VertField<double> phi(mesh,0.); setCircle(mesh,phi,.5,.5);
    InterfaceGeometry load(mesh,phi);
    rejects([&] { DielectricEvolution::explicitRelaxationCfl(c,100.); });
    const auto retained=service.evaluate(mesh,phi,{1,0.,.001},&load);
    const auto oldPotential=retained->electric.fields->potential;
    const auto oldCharge=retained->charge->integratedCharge;
    rejects([&] { service.completeCharge({1,0.,.001}); });
    checkContract(!service.failed(),"missing material map consumed current");
    LevelSetHermiteGeometry::Snapshot next(mesh,phi);
    const std::function<Vector2D(Vector2D)> broken=[](Vector2D)->Vector2D {
        throw std::runtime_error("injected material trajectory failure");
    };
    DielectricChargeMotion motion{mesh,next,broken};
    bool threw=false;
    try { service.completeCharge({1,0.,.001},nullptr,&motion); }
    catch(const std::runtime_error&) { threw=true; }
    checkContract(threw && service.failed(),"numerical failure did not stop state");
    checkContract(retained->electric.fields->potential==oldPotential && retained->charge->integratedCharge==oldCharge,
        "failed update invalidated an earlier successful publication");
    rejects([&] { fields(service); }); rejects([&] { service.charge(); });
    rejects([&] { service.evaluate(mesh,phi,{2,.001,.001},&load); });
}
void dielectricReference(QuadtreeMesh& mesh)
{
    DielectricConfig c; c.permittivityNegative=c.permittivityPositive=1.; c.upperPotential=.01;
    DielectricEvolution service(mesh,c,.08);
    VertField<double> phi(mesh,0.); setCircle(mesh,phi,.5,.5);
    for(int step=1; step<=2; ++step) {
        const auto evaluated=service.evaluate(mesh,phi,{step,(step-1)*.001,.001});
        const auto& report=evaluated->electric.report;
        checkContract(report.success && report.usedSuperLU && !service.charge(),"perfect dielectric created charge");
        const auto& potential=fields(service).potential;
        for(int i=0; i<mesh.cellIndexNum(); ++i) {
            const auto y=mesh.cells(mesh.cellIndexToLabel[i]).centerCoord().y;
            checkContract(std::abs(potential[i]-.01*y)<1e-12,"uniform dielectric potential differs from analytic solution");
        }
    }
}
void polarizationReference(QuadtreeMesh& mesh)
{
    DielectricConfig c;c.upperPotential=.01;
    DielectricEvolution service(mesh,c,.08);
    VertField<double> phi(mesh,0.);
    for(int i=0;i<mesh.vertIndexNum();++i)
        phi.value[i]=mesh.verts(mesh.vertIndexToLabel[i]).coord().x-.5;
    const auto evaluated=service.evaluate(mesh,phi,{1,0.,.001});
    const auto& report=evaluated->electric.report;
    const auto& published=fields(service);
    checkContract(report.usedSuperLU && published.polarization && !published.traceGeometry,"diffuse polarization was not published");
    const auto& polarization=*published.polarization;
    polarization.meshIdentity.validate(mesh);
    FaceSegmentGraph graph;graph.build(mesh);
    for(int i=0;i<mesh.cellIndexNum();++i) {
        const auto p=mesh.cells(mesh.cellIndexToLabel[i]).centerCoord();
        checkContract(std::abs(published.potential[i]-.01*p.y)<1e-12,"layered diffuse potential differs from analytic reference");
        checkContract(published.permittivity[i]==polarization.permittivityByCell[i],"polarization lost its actual material");
    }
    for(const auto& f:graph.segments()) {
        const int row=f.isHorizontal?mesh.vEdgeIndexNum()+mesh.hEdges(f.edgeLabel).index():mesh.vEdges(f.edgeLabel).index();
        const double expected=f.boundary==Boundary::non?-.5*.01*.01:0.;
        checkContract(std::abs(polarization.coefficientByNative[row]-expected)<1e-14,"polarization coefficient differs from Maxwell reference");
    }
}
void geometryQueryContract(QuadtreeMesh& mesh)
{
    auto c=config();
    VertField<double> phi(mesh,0.);setCircle(mesh,phi,.5,.5);
    InterfaceGeometry geometry(mesh,phi), other(mesh,phi);
    const auto queries=traceQueries(geometry);
    DielectricEvolution service(mesh,c,0.);
    rejects([&]{service.evaluate(mesh,phi,{1,0.,.001},&geometry,nullptr,traceQueries(other));});
    checkContract(!service.failed(),"foreign query source poisoned preflight");
    const double saved=phi.value[0];phi.value[0]+=.01;
    rejects([&]{service.evaluate(mesh,phi,{1,0.,.001},&geometry);});
    phi.value[0]=saved;
    const auto evaluated=service.evaluate(mesh,phi,{1,0.,.001},&geometry,nullptr,queries);
    const auto& report=evaluated->electric.report;
    checkContract(report.localTraceFits==queries->points().size(),"physical query collection was not sampled completely");
    auto points=queries->points();points.front().frame.normal.x+=.1;
    const auto wrong=std::make_shared<const InterfaceTraceGeometry>(geometry.source(),std::move(points));
    DielectricElectrostatics direct(mesh);
    SurfaceChargeField reference(mesh,geometry.snapshot(),c.initialSurfaceChargeDensity);
    const auto source=reference.source();
    const auto before=direct.solve(mesh,phi,0.,c,&geometry,nullptr,&source,nullptr,queries);
    checkContract(before.report.success && before.fields,"valid trace query failed");
    const auto potential=before.fields->potential;
    const auto failed=direct.solve(mesh,phi,0.,c,&geometry,nullptr,&source,nullptr,wrong);
    checkContract(!failed.report.success && !failed.fields && before.fields->traceGeometry==queries,"invalid scalar jet was committed");
    equalField(before.fields->potential,potential);
}
void publicationLifetime(QuadtreeMesh& mesh)
{
    std::shared_ptr<const DielectricEvaluation> retained;
    std::shared_ptr<const DielectricChargeSnapshot> advanced;
    std::vector<double> oldPotential,oldCharge;
    {
        auto c=config();
        DielectricEvolution service(mesh,c,0.);
        VertField<double> phi(mesh,0.);setCircle(mesh,phi,.5,.5);
        InterfaceGeometry geometry(mesh,phi);
        retained=service.evaluate(mesh,phi,{1,0.,.001},&geometry);
        oldPotential=retained->electric.fields->potential;oldCharge=retained->charge->integratedCharge;
        auto completed=service.completeCharge({1,0.,.001});advanced=completed.charge;
        const auto next=service.evaluate(mesh,phi,{2,.001,.001},&geometry);
        checkContract(next!=retained && next->electric.fields!=retained->electric.fields
            && next->charge==advanced && advanced->integratedCharge!=oldCharge,
            "later field/charge was not published as a distinct actual state");
        checkContract(retained->electric.fields->potential==oldPotential && retained->charge->integratedCharge==oldCharge,
            "later solve overwrote an earlier publication");
    }
    retained->electric.fields->meshIdentity.validate(mesh);
    advanced->meshIdentity.validate(mesh);
    checkContract(retained->electric.fields->potential==oldPotential && retained->charge->integratedCharge==oldCharge
        && advanced->completedStep==1 && advanced->time==.001 && !advanced->geometry->panels().empty(),
        "publications borrowed destroyed service/charge storage");
    checkContract(advanced->density().size()==advanced->geometry->panels().size(),"retained charge density lost its owned geometry");
}
}
int main()
try {
    domainInputContract();
    DictionaryMesh input; input.Nx=input.Ny=16; input.maxCellNum=1024; input.maxInitialIterations=0;
    QuadtreeMesh mesh(input);
    conservationContract(mesh);
    for(int mode=0; mode<3; ++mode) runTransport(mesh,mode);
    failureContract(mesh); dielectricReference(mesh); polarizationReference(mesh); geometryQueryContract(mesh); publicationLifetime(mesh);
    std::cout<<"Electricity alone: three transport policies match original solve/charge; lifecycle and analytic potential passed\n";
    return 0;
} catch(const std::exception& error) { std::cerr<<error.what()<<'\n'; return 1; }
