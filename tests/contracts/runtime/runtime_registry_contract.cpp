#include "AppRegistry.h"
#include "Electrohydrodynamics/IO/VofElectrohydrodynamicsInput.h"
#include "Electrohydrodynamics/VofElectrohydrodynamicsRun.h"
#include "Electrohydrodynamics/IO/VofElectrohydrodynamicsOutput.h"
#include "Hydrodynamics/HydrodynamicsVofRegrid.h"
#include "Modules/Hydrodynamics/HydroVofStep.h"
#include "Modules/Hydrodynamics/HydroVofMeshTransfer.h"
#include "Core/Geometry/Detail/FiniteVolumeClosure.h" // F26 exact independent Initial oracle.
#include "Modules/Electric/ElectricVofSharpField.h"
#include "Modules/Interface/VofSmoothGeometry.h"
#include "Modules/Interface/GeometricVofInitialization.h"
#include "Modules/Interface/GeometricVofMeshTransfer.h"
#include "Core/IO/DictionaryMesh.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Geometry/FittedDomainGeometry.h"
#include "Core/Mesh/FaceSegmentGraph.h"
#include <cfenv>
#include <cerrno>
#include <cstring>
#include <set>
#if defined(_M_X64) || defined(__SSE2__)
#include <xmmintrin.h>
#endif
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iterator>
#include <memory>
#include <type_traits>
#include <array>
#include <chrono>
#include <filesystem>
#include <fstream>
#include "RunControls.h"
#include "Hydrodynamics/HydrodynamicsTimeStep.h"

#include <iomanip>
#include <limits>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <utility>

static void require(bool ok, const char* message) {
    if (!ok) throw std::runtime_error(message);
}
struct Capture {
    std::ostringstream out, err;
    std::streambuf* coutBuffer = std::cout.rdbuf(out.rdbuf());
    std::streambuf* cerrBuffer = std::cerr.rdbuf(err.rdbuf());
    ~Capture() { std::cout.rdbuf(coutBuffer); std::cerr.rdbuf(cerrBuffer); }
};
struct Probe : IApp {
    std::vector<std::string>& events;
    RunConfig config;
    int result;
    bool fail;
    Probe(std::vector<std::string>& e, const RunConfig& c, int r, bool f = false)
        : events(e), config(c), result(r), fail(f) { events.push_back("create"); }
    ~Probe() override { events.push_back("destroy"); }
    std::string name() const override { return "Probe"; }
    int run() override {
        events.push_back("run");
        require(config.checkCase,
            "dispatcher changed prepared configuration");
        if (fail) throw std::runtime_error("probe failure");
        return result;
    }
};
static void checkFixedSteps() {
    RunControls run;
    require(!run.fixedStepCount(),"empty time settings accepted");
    run.dt=.01;
    for(double end:{.03,std::nextafter(.03,0.),std::nextafter(.03,1.)}) {
        run.endTime=end;
        require(run.fixedStepCount()==3,"roundoff changed the final step");
    }
    for(double end:{0.,-.03,.031,std::numeric_limits<double>::infinity(),
            std::numeric_limits<double>::quiet_NaN()}) {
        run.endTime=end;require(!run.fixedStepCount(),"invalid end time accepted");
    }
    run.endTime=.03;
    for(double dt:{0.,-.01,1.,std::numeric_limits<double>::infinity(),
            std::numeric_limits<double>::quiet_NaN(),std::numeric_limits<double>::min()}) {
        run.dt=dt;require(!run.fixedStepCount(),"invalid step size or unrepresentable count accepted");
    }
    run.dt=1.;run.endTime=std::numeric_limits<int>::max();
    require(run.fixedStepCount()==std::numeric_limits<int>::max(),"maximum integer step lost");
    run.endTime+=1.;require(!run.fixedStepCount(),"step count overflow accepted");
    run.dt=2.;run.endTime=6.;require(run.fixedStepCount()==3,"step count depends on time units");
}
static void checkAccumulatedEndpoint() {
    RunControls run; run.endTime=.5;
    double begin=0., dt=.001, smallest=1.; int accepted=0;
    while(begin<run.endTime) {
        require(accepted<901,"rounded clock created an extra tail step");
        const auto end=run.stepEnd(begin,dt,accepted);
        const auto actual=end-begin;
        smallest=std::min(smallest,actual);
        begin=end; ++accepted; dt=actual;
        if(accepted==274 || accepted==552) dt*=.5;
    }
    // Independently specified physical schedule: .274+.139+.087 = .5.
    require(accepted==900 && begin==.5 && smallest>.000249,
        "retry-sized physical schedule lost exact endTime or produced an ulp step");
    require(run.stepEnd(.25,.01,500)==.25+.01,"ordinary step was lengthened");
    require(run.stepEnd(.49,.02,500)==.5,"real final interval was not truncated");
}
static void checkAdaptiveTimeStep() {
    RunControls run;run.dt=.125;run.endTime=1.;
    run.adaptive=RunControls::Adaptive{.01,.5,.4,.5,2.};
    const auto select=[&](double begin,double previous,int accepted,long double flow,
        bool capillary=false,long double capillaryRate=0,double h=0){
        return selectHydrodynamicsTimeStep(run,begin,previous,accepted,flow,capillary,capillaryRate,h,.5);
    };
    const auto near=[](double a,double b){return std::abs(a-b)<=8*std::numeric_limits<double>::epsilon()*std::max(std::abs(a),std::abs(b));};
    const auto seed=select(0.,0.,0,0);
    require(seed.proposed_dt==.125&&seed.actual_dt==.125&&std::string(seed.limiter)=="initialDt",
        "adaptive first step lost the Time.dt seed");
    require(seed.flow_rate_available&&seed.flow_rate==0&&!seed.flow_bound_available&&std::isinf(seed.flow_bound)
        &&!seed.capillary_available&&!seed.capillary_bound_available,
        "zero measured flow/sigma must be explicit unconstrained bounds");
    const auto growth=select(.125,.125,1,0);
    require(growth.actual_dt==.25&&growth.growth_bound==.25L&&std::string(growth.limiter)=="growth",
        "adaptive growth does not use the previous actual accepted interval");
    const auto afterRetry=select(.25,.03125,8,0);
    require(afterRetry.actual_dt==.0625,"retry-reduced accepted dt was not the next growth source");
    const auto flow=select(0.,0.,0,4);
    require(flow.flow_bound_available&&near(flow.actual_dt,.05)&&std::string(flow.limiter)=="flow"
        &&static_cast<long double>(flow.actual_dt)<=flow.flow_bound,"actual flow bound violated");
    const auto capillary=select(0.,0.,0,0,true,8,.0625);
    require(capillary.actual_dt==.03125&&capillary.capillary_bound_available&&capillary.interface_spacing_available
        &&capillary.interface_minimum_spacing==.0625&&std::string(capillary.limiter)=="capillary",
        "existing capillary Courant and accepted interface h were not retained");
    const auto maximum=select(0.,.5,5,0);
    require(maximum.actual_dt==.5&&std::string(maximum.limiter)=="maxDt","adaptive maximum bound ignored");
    const auto rejected=[&](auto call){bool failed=false;try{call();}catch(const std::invalid_argument&){failed=true;}return failed;};
    require(rejected([&]{select(0.,0.,0,100);}),"stability below minDt was clamped upward");
    const auto final=select(.99609375,.125,8,4);
    require(final.final&&final.below_min_final&&final.end==1.&&final.actual_dt==.00390625
        &&std::string(final.limiter)=="endTime","stable true final interval below minDt rejected");
    require(rejected([&]{select(.99609375,.125,8,100);}),"unstable true final interval below minDt accepted");
    require(rejected([&]{select(0.,0.,0,-1);})
        &&rejected([&]{select(0.,0.,0,std::numeric_limits<long double>::quiet_NaN());})
        &&rejected([&]{select(0.,0.,0,0,true,0,.125);}),"invalid/unknown rates became unconstrained");
    run.stepRetry=RunControls::StepRetry{.01,4,1};
    require(rejected([&]{select(0.,0.,0,0);}),"retry and adaptive both owned growth");
    run.stepRetry.reset();
    run.dt=.001;run.endTime=.5;run.adaptive=RunControls::Adaptive{.00001,.001,.4,.9,1.25};
    double begin=0.,previous=0.,smallest=1.;int accepted=0;bool aligned=false;
    while(begin<run.endTime){
        require(accepted<501,"adaptive endpoint generated an ulp tail step");
        const auto next=select(begin,previous,accepted,0);
        previous=next.actual_dt;smallest=std::min(smallest,previous);begin=next.end;
        aligned=aligned||next.endpoint_aligned;++accepted;
    }
    require(accepted==500&&begin==.5&&smallest>.00099&&aligned,
        "adaptive nominal maximum lost cumulative endpoint alignment");
    // min=max=seed is a legitimate nominal policy, despite the accepted
    // endpoint subtraction occasionally lying one or two clock ulps below it.
    run.adaptive=RunControls::Adaptive{.001,.001,.4,.9,1.25};
    begin=previous=0.;smallest=1.;accepted=0;bool minimumRounded=false;
    while(begin<run.endTime){
        require(accepted<501,"equal min/max policy generated an ulp tail step");
        const auto next=select(begin,previous,accepted,0);
        require(next.proposed_dt>=run.adaptive->minDt,
            "equal min/max policy selected a genuinely below-min proposal");
        if(next.actual_dt<run.adaptive->minDt)
            require(next.minimum_clock_rounded||next.below_min_final,
                "short actual interval has no explicit bounded clock provenance");
        minimumRounded=minimumRounded||next.minimum_clock_rounded;
        previous=next.actual_dt;smallest=std::min(smallest,previous);begin=next.end;++accepted;
    }
    require(accepted==500&&begin==.5&&smallest>.00099&&minimumRounded,
        "equal min/max policy mistook clock subtraction ulps for a below-min choice");
}
static void check() {
    Capture capture;
    AppRegistry registry;
    RunConfig cfg;
    std::vector<std::string> events;
    require(registry.run(cfg) == 1 && capture.out.str().empty() && capture.err.str().empty(), "empty selection");
    cfg.app = "missing";
    require(registry.run(cfg) == 2 && capture.err.str().find("Unknown app: 'missing'") != std::string::npos,
        "unknown selection");
    cfg.checkCase = true;
    const auto factory = [&](const RunConfig& input) -> std::unique_ptr<IApp> {
        require(&input == &cfg, "factory did not receive caller configuration");
        return std::make_unique<Probe>(events, input, 0);
    };
    registry.registerApp("Zulu", factory);
    registry.registerApp("Alpha", factory);
    require(registry.registeredNames() == std::vector<std::string>({"alpha", "zulu"}), "stable normalized names");
    cfg.app = "aLpHa";
    require(registry.run(cfg) == 0 && registry.run(cfg) == 0, "case-insensitive factory execution");
    require(events == std::vector<std::string>({"create", "run", "destroy", "create", "run", "destroy"}),
        "each execution must own a fresh App until completion");
    require(cfg.app == "aLpHa" && cfg.checkCase, "caller configuration changed");
    events.clear();
    registry.registerApp("ALPHA", [&](const RunConfig& input) { return std::make_unique<Probe>(events, input, 37); });
    require(registry.run(cfg) == 37 && registry.registeredNames().size() == 2, "replacement/nonzero return code");
    require(capture.err.str().find("App 'Probe' returned rc=37") != std::string::npos, "failure diagnostic");
    require(events == std::vector<std::string>({"create", "run", "destroy"}), "nonzero App destruction");
    registry.registerApp("Null", [](const RunConfig&) -> std::unique_ptr<IApp> { return {}; });
    cfg.app = "null";
    require(registry.run(cfg) == 2 && capture.err.str().find("Available: alpha null zulu") != std::string::npos,
        "null factory result and sorted available names");
    events.clear();
    registry.registerApp("Throws", [&](const RunConfig& input) { return std::make_unique<Probe>(events, input, 0, true); });
    cfg.app = "throws";
    bool rejected = false;
    try { registry.run(cfg); } catch (const std::runtime_error& e) { rejected = std::string(e.what()) == "probe failure"; }
    require(rejected && events == std::vector<std::string>({"create", "run", "destroy"}), "App exception propagation/destruction");
    registry.registerApp("FactoryThrows", [](const RunConfig&) -> std::unique_ptr<IApp> { throw std::runtime_error("factory failure"); });
    cfg.app = "factorythrows"; rejected = false;
    try { registry.run(cfg); } catch (const std::runtime_error& e) { rejected = std::string(e.what()) == "factory failure"; }
    require(rejected, "factory exception propagation");
    AppRegistry independent;
    cfg.app = "alpha";
    require(independent.run(cfg) == 2, "registry instances must not share factories");
}
static void checkFieldSchedule() {
    RunControls run;
    run.dt=.03;run.endTime=.27;run.writeInterval=.08;
    run.writeInitial=false;run.writeFirst=false;run.writeFinal=false;
    const std::vector<int> expected{3,6,8};
    std::vector<int> selected;
    for(int step=0;step<=9;++step) if(run.writesFieldsAtStep(step,9)) selected.push_back(step);
    require(selected==expected,"time intervals did not select the first completed step at each boundary");
    run.writeFirst=true;
    require(run.writesFieldsAtStep(1,9) && !run.writesFieldsAtStep(2,9)
        && run.writesFieldsAtStep(3,9),"forced first frame shifted regular cadence");
    require(run.writesFieldsAtStep(8,9) && run.writesFieldsAtStep(3,9),"field selection depends on query order");
    run.writeFirst=false;run.dt=.01;run.writeInterval=std::nextafter(.03,1.);
    require(!run.writesFieldsAtStep(2,3) && run.writesFieldsAtStep(3,3),"rounded interval boundary lost its frame");
    run.writeEverySteps=4;run.writeInitial=true;run.writeFirst=true;run.writeFinal=true;
    selected.clear();
    for(int step=0;step<=9;++step) if(run.writesFieldsAtStep(step,9)) selected.push_back(step);
    require(selected==std::vector<int>({0,1,4,8,9}),"step interval or initial/first/final flags changed");
    require(!run.writesFieldsAtStep(-1,9) && !run.writesFieldsAtStep(10,9)
        && !run.writesFieldsAtStep(0,0),"out-of-run step selected");
    run.writeFields=false;
    for(int step=0;step<=9;++step) require(!run.writesFieldsAtStep(step,9),"disabled fields selected a forced frame");
    run.writeFields=true;run.writeInitial=false;run.writeFirst=false;run.writeFinal=false;run.writeEverySteps=0;
    run.dt=1e100;run.endTime=3e100;run.writeInterval=1e-301;
    require(run.fixedStepCount()==3,"overflow probe is not a valid run");
    for(int step=1;step<=3;++step)
        require(run.writesFieldsAtStep(step,3),"substep interval overflow lost an accepted frame");
    run.dt=1e-100;run.writeInterval=2e-100;
    const int last=std::numeric_limits<int>::max();
    require(run.writesFieldsAtStep(last-1,last) && !run.writesFieldsAtStep(last,last),
        "large step index overflowed or changed the interval");
}
namespace VofEhdInputContract
{
namespace fs=std::filesystem;
constexpr std::array<const char*,7> names{{"run.dict","geometry.dict","mesh.dict","interface.dict","materials.dict","hydro.dict","electric.dict"}};
unsigned checks=0,guards=0;
void check(bool ok,const char* why){++checks;if(!ok)throw std::runtime_error(why);}
void write(const fs::path& path,const std::string& text)
{std::ofstream out(path,std::ios::binary);out<<text;if(!out)throw std::runtime_error("cannot save reader fixture");}
struct Fixture
{
    std::array<std::string,7> text;
    // Edit the real Geometry document explicitly. save() writes seven current
    // documents verbatim; there is no old Fields slot or schema conversion.
    void setConstants(const std::string& body)
    {
        auto& geometry=text[1];const auto old=geometry.find("Constants {");
        if(old!=std::string::npos){const auto end=geometry.find('}',old);
            check(end!=std::string::npos,"fixture Constants block is incomplete");geometry.erase(old,end-old+1);}
        const auto begin=geometry.find('{');check(begin!=std::string::npos,"fixture Geometry block is missing");
        geometry.insert(begin+1," Constants { "+body+" } ");
    }

    explicit Fixture(bool rz=false){
        text[0]="Run { app VofElectrohydrodynamics; } Parallel { threads 1; } Time { units code; dt .001; endTime .002; Adaptive { minDt .00025; maxDt .002; targetCourant .2; } } Output { Fields { enabled false; } }";
        text[1]=std::string("Geometry { coordinates ")+(rz?"rz":"xy")+"; Domain { min (0 0); max (1 1); } Phases { negative 7; positive 3; } Initial { phi \""+
            (rz?"r^2+(z-.5)^2-radius^2":"(x-.5)^2+(y-.5)^2-radius^2")+"\"; } }";
        text[2]="Mesh { resolution (8 8); initialMesh interfaceBand; Refinement { level 1; bandCells 2; } Regrid { minLevel 0; maxLevel 1; interfaceBandLayers 2; vorticityTolerance .5; coarsenRatio .25; } }";
        setConstants("radius .125; chargeScale 2;");
        text[3]="Interface { method vof; Boundary { type none; } Initialization { subcellsPerAxis 8; } }";
        // Source order is deliberately gas,liquid and IDs are deliberately not 1,2.
        text[4]="Materials { Phases { Gas { phaseId 3; density 1; dynamicViscosity 0; permittivity 2; conductivity 0; } Liquid { phaseId 7; density 2; dynamicViscosity 0; permittivity 4; conductivity .25; } } Interfaces { Shared { phases (7 3); surfaceTension 0; } } }";
        text[5]="Hydrodynamics { method staggered; densityRepresentation cellAverageProjection; Initial { velocity (0 0); } Boundary { type impermeable; } BodyForce { acceleration (0 0); } TimeIntegration { method adaptiveLinear; } }";
        text[6]=std::string("Electric { method vofFiniteVolume; Initial { chargeDensity \"")+(rz?"chargeScale*(z-.5)":"chargeScale*(x-.5)")+"\"; } Boundary { W { type "+(rz?"Neumann; value 0":"Dirichlet; value -.25")+"; } S { type Dirichlet; value -.5; } E { type Dirichlet; value 1; } N { type Neumann; value 0; } } }";
    }
    fs::path save(const fs::path& root,const std::string& name) const
    {const auto path=root/name;fs::create_directories(path);for(std::size_t i=0;i<names.size();++i)write(path/names[i],text[i]);return path;}
};
void replace(std::string& text,const std::string& from,const std::string& to)
{const auto p=text.find(from);check(p!=std::string::npos,"reader negative fixture anchor missing");check(text.find(from,p+from.size())==std::string::npos,"reader negative fixture anchor not unique");text.replace(p,from.size(),to);}
template<class Action> void reject(const fs::path& root,Fixture fixture,const std::string& name,const char* field,Action change)
{
    change(fixture);const auto path=fixture.save(root,name);bool caught=false;
    try{readVofElectrohydrodynamicsInput(path.string());}catch(const std::exception& e){caught=true;check(std::string(e.what()).find(field)!=std::string::npos,"reader guard lost actual field/reason");}
    check(caught,"invalid VOF-EHD input was accepted");++guards;
}
void positive(const fs::path& root,bool rz)
{
    Fixture fixture(rz);if(rz)replace(fixture.text[6],"Electric {","Electric { Solver { relativeLinearTolerance 1e-12; relativeConservationTolerance 1e-11; }");
    const auto path=fixture.save(root,rz?"rz":"xy");const auto input=readVofElectrohydrodynamicsInput(path.string(),3);const auto& fluid=input.fluid;
    check(fluid.interfaceMethod==HydrodynamicsInput::InterfaceMethod::Vof&&fluid.radial==rz,"reader selected wrong method/measure");
    check(fluid.parallelThreads==3&&fluid.run.adaptive&&fluid.vofRegrid&&fluid.initialInterfaceBand,"actual shared reader/adaptive/refinement were bypassed");
    check(!fluid.run.diagnosticSeries,"omitted series profile must default to physical histories");
    check(fluid.negativePhase==7&&fluid.positivePhase==3&&input.materials.epsilonLiquid==4&&input.materials.epsilonGas==2
        &&input.materials.conductivityLiquid==.25&&input.materials.conductivityGas==0,"material ID/source order or legal zero conductivity changed");
    check(input.initialChargeDensity({.25,.25})==-.5&&input.initialChargeDensity({.75,.75})==.5,"signed volume-charge expression was clipped, complemented or unsigned");
    check(input.boundaries[0].type==(rz?LAMP::Electric::ElectricVofBoundary::Type::HomogeneousNeumann:LAMP::Electric::ElectricVofBoundary::Type::Dirichlet)
        &&input.boundaries[3].type==LAMP::Electric::ElectricVofBoundary::Type::HomogeneousNeumann,"static native boundary/true axis mapping changed");
    check(input.controls.relativeLinearTolerance==(rz?1e-12:1e-11)&&input.controls.relativeConservationTolerance==(rz?1e-11:1e-10),"actual original solver controls/defaults changed");
    check(fluid.documents.size()==7,"VOF-EHD must retain exactly seven real documents");
    const auto effective=root/(rz?"rz_effective":"xy_effective");fs::create_directories(effective);
    for(std::size_t i=0;i<names.size();++i){check(fluid.documents[i].sourceText==fixture.text[i],"reader changed source dictionary text");
        check(fs::path(fluid.documents[i].file).filename()==names[i],"document provenance/order changed");write(effective/names[i],fluid.documents[i].effectiveText);}
    const auto roundtrip=readVofElectrohydrodynamicsInput(effective.string());
    check(roundtrip.fluid.documents.size()==7&&roundtrip.fluid.parallelThreads==3&&roundtrip.fluid.radial==rz
        &&roundtrip.materials.epsilonLiquid==4&&roundtrip.materials.conductivityGas==0&&roundtrip.controls.relativeLinearTolerance==input.controls.relativeLinearTolerance
        &&roundtrip.initialChargeDensity({.25,.25})==-.5,"seven-document effective text does not roundtrip the actual typed reader");
    bool chargeTrace=false,solverTrace=false;
    for(const auto& entry:fluid.trace.entries){if(entry.field=="Electric.Initial.chargeDensity")chargeTrace=entry.units=="code_charge/code_volume"&&entry.used&&!entry.defaulted;
        if(entry.field=="Electric.Solver.relativeLinearTolerance")solverTrace=entry.defaulted==!rz&&entry.used;}
    check(chargeTrace&&solverTrace,"charge units/default-vs-explicit solver provenance missing");
}
void seriesProfiles(const fs::path& root)
{
    for(const auto* profile:{"physical","diagnostic"}){
        Fixture fixture;replace(fixture.text[0],"Output {",std::string("Output { Series { profile ")+profile+"; }");
        const auto input=readVofElectrohydrodynamicsInput(fixture.save(root,std::string("series_")+profile).string());
        const bool diagnostic=std::string(profile)=="diagnostic";
        check(input.fluid.run.diagnosticSeries==diagnostic,"explicit series profile was not consumed");
        bool traced=false;for(const auto& entry:input.fluid.trace.entries)if(entry.field=="Output.Series.profile")
            traced=entry.used&&!entry.defaulted&&fs::path(entry.file).filename()=="run.dict";
        check(traced,"explicit series profile lost actual source/use provenance");
        const auto effective=root/(std::string("series_effective_")+profile);fs::create_directories(effective);
        for(std::size_t i=0;i<names.size();++i)write(effective/names[i],input.fluid.documents[i].effectiveText);
        const auto again=readVofElectrohydrodynamicsInput(effective.string());
        check(again.fluid.run.diagnosticSeries==diagnostic,"effective series profile did not roundtrip");
    }
    reject(root,Fixture{},"series_bad_profile","Output.Series.profile",[](Fixture& f){
        replace(f.text[0],"Output {","Output { Series { profile verbose; }");});
    reject(root,Fixture{},"series_unknown_control","Output.Series.futureFlag",[](Fixture& f){
        replace(f.text[0],"Output {","Output { Series { futureFlag true; }");});
}
void sharpAndCoupling(const fs::path& root,const Fixture& xy)
{
    Fixture sharp;sharp.text[2]="Mesh { resolution (8 8); initialMesh uniform; }";
    replace(sharp.text[6],"vofFiniteVolume","vofSharpSurface");replace(sharp.text[6],"chargeDensity","surfaceChargeDensity");
    const auto b1=readVofElectrohydrodynamicsInput(sharp.save(root,"b1_xy").string());
    check(b1.method==VofElectrohydrodynamicsInput::Method::VofSharpSurface&&b1.initialSurfaceChargeDensity({.25,.25})==-.5,
        "B1 did not consume a signed surface expression");
    check(b1.fluid.vof.requireMaxwellCorrection&&b1.fluid.vof.maximumMaxwellLoadChange==.25&&b1.fluid.vof.maximumMaxwellWorkChange==.25,
        "B1 default actual correction or original force/work budget changed");
    bool wrongAccessor=false;try{b1.initialChargeDensity({.25,.25});}catch(const std::logic_error&){wrongAccessor=true;}
    check(wrongAccessor,"B1 surface density was exposed as volume density");
    bool surfaceTrace=false;for(const auto& entry:b1.fluid.trace.entries)if(entry.field=="Electric.Initial.surfaceChargeDensity")
        surfaceTrace=entry.units=="code_charge/code_interface_measure"&&entry.used&&!entry.defaulted
            &&fs::path(entry.file).filename()=="electric.dict";
    check(surfaceTrace,"B1 surface units or actual consumption trace missing");
    auto explicitSharp=sharp;replace(explicitSharp.text[6],"Electric {","Electric { SmoothGeometry { relativeVolumeTolerance 5e-13; } SurfaceTransport { maximumDisplacementInCells .4; maximumProjectionDistanceInCells .15; relativeChargeTolerance 5e-11; } Coupling { maximumForceChange .2; maximumWorkChange .125; }");
    const auto declared=readVofElectrohydrodynamicsInput(explicitSharp.save(root,"b1_explicit_budgets").string());
    check(declared.smoothGeometry.relativeVolumeTolerance==static_cast<long double>(5e-13)
        &&declared.surfaceTransport.maximumDisplacementInCells==static_cast<long double>(.4)
        &&declared.surfaceTransport.maximumProjectionDistanceInCells==static_cast<long double>(.15)
        &&declared.surfaceTransport.relativeChargeTolerance==static_cast<long double>(5e-11)
        &&declared.fluid.vof.maximumMaxwellLoadChange==.2&&declared.fluid.vof.maximumMaxwellWorkChange==.125,
        "B1 declared geometry/transport/correction budgets were not consumed");
    const auto sharpEffective=root/"b1_effective";fs::create_directories(sharpEffective);
    for(std::size_t i=0;i<names.size();++i){check(declared.fluid.documents[i].sourceText==explicitSharp.text[i],"B1 reader changed actual source text");
        write(sharpEffective/names[i],declared.fluid.documents[i].effectiveText);}
    const auto b1Roundtrip=readVofElectrohydrodynamicsInput(sharpEffective.string());
    check(b1Roundtrip.method==declared.method&&b1Roundtrip.initialSurfaceChargeDensity({.25,.25})==-.5
        &&b1Roundtrip.smoothGeometry.relativeVolumeTolerance==declared.smoothGeometry.relativeVolumeTolerance
        &&b1Roundtrip.surfaceTransport.maximumDisplacementInCells==declared.surfaceTransport.maximumDisplacementInCells
        &&b1Roundtrip.surfaceTransport.maximumProjectionDistanceInCells==declared.surfaceTransport.maximumProjectionDistanceInCells
        &&b1Roundtrip.surfaceTransport.relativeChargeTolerance==declared.surfaceTransport.relativeChargeTolerance
        &&b1Roundtrip.fluid.vof.requireMaxwellCorrection&&b1Roundtrip.fluid.vof.maximumMaxwellLoadChange==.2
        &&b1Roundtrip.fluid.vof.maximumMaxwellWorkChange==.125,"B1 effective input did not roundtrip");
    const auto a=readVofElectrohydrodynamicsInput(xy.save(root,"a_original_coupling_policy").string());
    check(a.method==VofElectrohydrodynamicsInput::Method::VofFiniteVolume&&!a.fluid.vof.requireMaxwellCorrection
        &&a.fluid.documents[6].effectiveText.find("Coupling")==std::string::npos,"A without Coupling changed its default correction policy");
    wrongAccessor=false;try{a.initialSurfaceChargeDensity({.25,.25});}catch(const std::logic_error&){wrongAccessor=true;}
    check(wrongAccessor,"A volume density was exposed as surface density");
    auto coupled=xy;replace(coupled.text[6],"Electric {","Electric { Coupling { maximumForceChange .25; maximumWorkChange .125; }");
    const auto aCoupled=readVofElectrohydrodynamicsInput(coupled.save(root,"a_explicit_coupling").string());
    check(aCoupled.fluid.vof.requireMaxwellCorrection&&aCoupled.fluid.vof.maximumMaxwellLoadChange==.25
        &&aCoupled.fluid.vof.maximumMaxwellWorkChange==.125,"A explicit Coupling block was not consumed");
    const auto aEffective=root/"a_coupling_effective";fs::create_directories(aEffective);
    for(std::size_t i=0;i<names.size();++i)write(aEffective/names[i],aCoupled.fluid.documents[i].effectiveText);
    const auto aRoundtrip=readVofElectrohydrodynamicsInput(aEffective.string());
    check(aRoundtrip.fluid.vof.requireMaxwellCorrection&&aRoundtrip.fluid.vof.maximumMaxwellLoadChange==.25
        &&aRoundtrip.fluid.vof.maximumMaxwellWorkChange==.125,"A Coupling effective input did not roundtrip");
    reject(root,sharp,"b1_bulk_key","Electric.Initial.chargeDensity",[](Fixture& f){replace(f.text[6],"surfaceChargeDensity","chargeDensity");});
    reject(root,sharp,"b1_both_keys","Electric.Initial.chargeDensity",[](Fixture& f){replace(f.text[6],"Initial {","Initial { chargeDensity 0;");});
    reject(root,sharp,"b1_rz","Electric.method",[](Fixture& f){replace(f.text[1],"coordinates xy","coordinates rz");replace(f.text[1],"(x-.5)^2+(y-.5)^2","(r-.5)^2+(z-.5)^2");});
    reject(root,sharp,"b1_regrid","Electric.method",[&](Fixture& f){f.text[2]=xy.text[2];});
    reject(root,sharp,"b1_initial_refinement","Electric.method",[](Fixture& f){f.text[2]="Mesh { resolution (8 8); initialMesh interfaceBand; Refinement { level 1; bandCells 2; } }";});
    reject(root,sharp,"b1_prescribed_flow","Electric.method",[](Fixture& f){replace(f.text[5],"Boundary { type impermeable; }",
        "Boundary { type prescribedFlow; Patches { feed { side S; interval (0 .5); type prescribedVelocity; phase liquid; velocity (.03125 .125); } out { side N; interval (0 1); type pressureOutlet; phase gas; pressure 2; } } }");});
    reject(root,sharp,"b1_symmetry_contact","Electric.method",[](Fixture& f){replace(f.text[3],"Boundary { type none; }","Boundary { type symmetry; sides (bottom); }");});
    reject(root,sharp,"b1_south_even90_contact","Electric.method",[](Fixture& f){replace(f.text[3],"Boundary { type none; }","Boundary { type southEven90; interval (0 1); }");});
    reject(root,sharp,"b1_timed_surface","Electric.Initial.surfaceChargeDensity",[](Fixture& f){replace(f.text[6],"chargeScale*(x-.5)","chargeScale*(x-.5)+t");});
    reject(root,sharp,"b1_looser_volume","relativeVolumeTolerance",[](Fixture& f){replace(f.text[6],"Electric {","Electric { SmoothGeometry { relativeVolumeTolerance 1e-11; }");});
    reject(root,sharp,"b1_looser_displacement","maximumDisplacementInCells",[](Fixture& f){replace(f.text[6],"Electric {","Electric { SurfaceTransport { maximumDisplacementInCells .46; }");});
    reject(root,sharp,"b1_looser_projection","maximumProjectionDistanceInCells",[](Fixture& f){replace(f.text[6],"Electric {","Electric { SurfaceTransport { maximumProjectionDistanceInCells .21; }");});
    reject(root,sharp,"b1_looser_surface","relativeChargeTolerance",[](Fixture& f){replace(f.text[6],"Electric {","Electric { SurfaceTransport { relativeChargeTolerance 1e-9; }");});
    reject(root,sharp,"b1_looser_force","maximumForceChange",[](Fixture& f){replace(f.text[6],"Electric {","Electric { Coupling { maximumForceChange .3; }");});
    reject(root,sharp,"b1_zero_work","maximumWorkChange",[](Fixture& f){replace(f.text[6],"Electric {","Electric { Coupling { maximumWorkChange 0; }");});
    reject(root,sharp,"b1_unused_advection","Electric.Advection",[](Fixture& f){replace(f.text[6],"Electric {","Electric { Advection { Inlets {} }");});
    reject(root,xy,"a_both_keys","Electric.Initial.surfaceChargeDensity",[](Fixture& f){replace(f.text[6],"Initial {","Initial { surfaceChargeDensity 0;");});
    reject(root,xy,"a_unused_smooth","Electric.SmoothGeometry",[](Fixture& f){replace(f.text[6],"Electric {","Electric { SmoothGeometry { relativeVolumeTolerance 1e-12; }");});
    reject(root,xy,"a_unused_surface","Electric.SurfaceTransport",[](Fixture& f){replace(f.text[6],"Electric {","Electric { SurfaceTransport { relativeChargeTolerance 1e-10; }");});
    reject(root,xy,"a_looser_force","maximumForceChange",[](Fixture& f){replace(f.text[6],"Electric {","Electric { Coupling { maximumForceChange .3; }");});
    reject(root,xy,"a_looser_work","maximumWorkChange",[](Fixture& f){replace(f.text[6],"Electric {","Electric { Coupling { maximumWorkChange .3; }");});
}
void meshPhasePolicy(const fs::path& root)
{
    Fixture legacy;
    const auto original=readVofElectrohydrodynamicsInput(legacy.save(root,"mesh_legacy_defaults").string());
    check(original.fluid.initialBandScale==HydrodynamicsInput::InitialBandScale::BaseCell
        &&!original.fluid.vofRegrid->liquidMinLevel&&!original.fluid.vofRegrid->gasMinLevel,
        "omitted mesh policy changed legacy base-cell band or inherited common floors");
    bool defaultScale=false;
    for(const auto& entry:original.fluid.trace.entries)if(entry.field=="Mesh.Refinement.bandScale")
        defaultScale=entry.used&&entry.defaulted;
    check(defaultScale,"legacy band scale lost default provenance");
    Fixture supplied;
    replace(supplied.text[2],"level 1; bandCells 2;","level 2; bandCells 2; bandScale currentCell;");
    replace(supplied.text[2],"minLevel 0; maxLevel 1;","minLevel 1; maxLevel 3; liquidMinLevel 2; gasMinLevel 1;");
    const auto input=readVofElectrohydrodynamicsInput(supplied.save(root,"mesh_phase_policy").string());
    const auto& f=input.fluid;const auto& policy=*f.vofRegrid;
    check(f.initialBandScale==HydrodynamicsInput::InitialBandScale::CurrentCell
        &&policy.minLevel==1&&policy.maxLevel==3&&policy.liquidMinLevel==2&&policy.gasMinLevel==1,
        "explicit current-cell band or separate phase floors were not consumed");
    bool scale=false,liquid=false,gas=false;
    for(const auto& entry:f.trace.entries){
        if(entry.field=="Mesh.Refinement.bandScale")scale=entry.used&&!entry.defaulted;
        if(entry.field=="Mesh.Regrid.liquidMinLevel")liquid=entry.used&&!entry.defaulted;
        if(entry.field=="Mesh.Regrid.gasMinLevel")gas=entry.used&&!entry.defaulted;
    }
    check(scale&&liquid&&gas,"explicit mesh policy lost source/use provenance");
    const auto effective=root/"mesh_policy_effective";fs::create_directories(effective);
    for(std::size_t i=0;i<names.size();++i)write(effective/names[i],f.documents[i].effectiveText);
    const auto reread=readVofElectrohydrodynamicsInput(effective.string());const auto& again=reread.fluid;
    check(again.initialBandScale==f.initialBandScale&&again.vofRegrid->minLevel==policy.minLevel
        &&again.vofRegrid->maxLevel==policy.maxLevel&&again.vofRegrid->liquidMinLevel==policy.liquidMinLevel
        &&again.vofRegrid->gasMinLevel==policy.gasMinLevel,
        "effective mesh policy did not roundtrip actual typed scale and floors");
    Fixture oneFloor=legacy;replace(oneFloor.text[2],"minLevel 0; maxLevel 1;","minLevel 1; maxLevel 2; liquidMinLevel 2;");
    const auto inherited=readVofElectrohydrodynamicsInput(oneFloor.save(root,"mesh_one_phase_floor").string());
    check(inherited.fluid.vofRegrid->liquidMinLevel==2&&!inherited.fluid.vofRegrid->gasMinLevel
        &&inherited.fluid.vofRegrid->gasMinLevel.value_or(inherited.fluid.vofRegrid->minLevel)==1,
        "one explicit phase floor changed omitted phase inheritance");
    reject(root,supplied,"mesh_bad_scale","bandScale",[](Fixture& f){replace(f.text[2],"currentCell","allCells");});
    reject(root,supplied,"mesh_liquid_below_common","liquidMinLevel",[](Fixture& f){replace(f.text[2],"liquidMinLevel 2;","liquidMinLevel 0;");});
    reject(root,supplied,"mesh_liquid_above_max","liquidMinLevel",[](Fixture& f){replace(f.text[2],"liquidMinLevel 2;","liquidMinLevel 4;");});
    reject(root,supplied,"mesh_gas_below_common","gasMinLevel",[](Fixture& f){replace(f.text[2],"gasMinLevel 1;","gasMinLevel 0;");});
    reject(root,supplied,"mesh_gas_above_max","gasMinLevel",[](Fixture& f){replace(f.text[2],"gasMinLevel 1;","gasMinLevel 4;");});
    reject(root,supplied,"mesh_nonintegral_floor","liquidMinLevel",[](Fixture& f){replace(f.text[2],"liquidMinLevel 2;","liquidMinLevel 1.5;");});
}

void geometryConstantsSchema(const fs::path& root)
{
    Fixture shared;replace(shared.text[6],"chargeScale*(x-.5)","chargeScale*radius*(x-.5)");
    const auto input=readVofElectrohydrodynamicsInput(shared.save(root,"geometry_shared_constant").string());
    check(input.fluid.documents.size()==7&&input.fluid.constants.at("radius")==.125&&input.fluid.constants.at("chargeScale")==2,
        "Geometry.Constants did not produce one shared typed constant map");
    check(input.fluid.initialJet({.5,.5}).value==-.015625&&input.initialChargeDensity({.25,.25})==-.0625,
        "one Geometry.Constants radius was not shared by geometry and electric initialization");
    unsigned traces=0;for(const auto& e:input.fluid.trace.entries){
        check(e.field.find("Fields.Constants")==std::string::npos,"retired Fields constant trace owner");
        if(e.field=="Geometry.Constants.radius"||e.field=="Geometry.Constants.chargeScale"){
            check(e.used&&!e.defaulted&&e.units=="code"&&!e.consumer.empty()&&fs::path(e.file).filename()=="geometry.dict",
                "shared constant trace does not identify actual geometry.dict source");++traces;}}
    check(traces==2,"shared Geometry constants have duplicate/missing trace entries");
    const auto effective=root/"geometry_shared_constant_effective";fs::create_directories(effective);
    for(std::size_t i=0;i<names.size();++i){check(fs::path(input.fluid.documents[i].file).filename()==names[i]
        &&input.fluid.documents[i].sourceText==shared.text[i],"seven-document constant provenance/order changed");
        write(effective/names[i],input.fluid.documents[i].effectiveText);}
    const auto again=readVofElectrohydrodynamicsInput(effective.string());
    check(again.fluid.documents.size()==7&&again.fluid.constants==input.fluid.constants&&again.initialChargeDensity({.25,.25})==-.0625,
        "Geometry.Constants effective seven-document roundtrip failed");
    reject(root,shared,"unused_geometry_constant","Geometry.Constants.unused",[](Fixture& f){f.setConstants("radius .125; chargeScale 2; unused 1;");});
    for(const auto name:{"x","t","sin"})reject(root,shared,std::string("reserved_geometry_constant_")+name,
        (std::string("Geometry.Constants.")+name).c_str(),[&](Fixture& f){f.setConstants(std::string("radius .125; chargeScale 2; ")+name+" 1;");});
    {const auto path=shared.save(root,"obsolete_fields_extra");write(path/"fields.dict","Fields { Constants { stale 1; } }");bool caught=false;
        try{readVofElectrohydrodynamicsInput(path.string());}catch(const std::exception& e){caught=std::string(e.what()).find("fields.dict")!=std::string::npos;}
        check(caught,"obsolete fields.dict was accepted with current seven documents");++guards;}
    for(bool explicitEmpty:{false,true}){
        Fixture empty;empty.text[1]="Geometry { coordinates xy; Domain { min (0 0); max (1 1); } Phases { negative 7; positive 3; } Initial { phi \"(x-.5)^2+(y-.5)^2-.125^2\"; } }";
        replace(empty.text[6],"chargeScale*(x-.5)","0");if(explicitEmpty)empty.setConstants("");
        if(!explicitEmpty)for(const auto raw:{"123","(1 2)"})
            reject(root,empty,std::string("geometry_constants_wrong_type_")+(raw[0]=='1'?"scalar":"tuple"),"Geometry.Constants",
                [raw](Fixture& f){f.text[1].insert(f.text[1].find('{')+1,std::string(" Constants ")+raw+"; ");});
        const auto value=readVofElectrohydrodynamicsInput(empty.save(root,explicitEmpty?"geometry_empty_constants":"geometry_no_constants").string());
        check(value.fluid.constants.empty()&&value.initialChargeDensity({.25,.25})==0&&value.fluid.documents.size()==7,
            "empty/omitted Geometry.Constants rejected a parameter-free input");
        for(const auto& e:value.fluid.trace.entries)check(e.field.find("Geometry.Constants.")==std::string::npos,
            "empty/omitted Constants fabricated an expression parameter");}
}

void conservativeMomentumSchema(const fs::path& root)
{
    for(bool rz:{false,true}){
        Fixture f(rz);replace(f.text[5],"densityRepresentation cellAverageProjection; ","");
        const auto input=readVofElectrohydrodynamicsInput(f.save(root,rz?"rz_default_conservative":"xy_default_conservative").string());
        check(input.fluid.vof.momentumDensityRepresentation==LAMP::Hydrodynamics::HydroMomentumMass::DensityRepresentation::CellAverageProjection,
            "VOF default did not select conservative inventory P0 mass");
        bool traced=false;for(const auto& e:input.fluid.trace.entries)if(e.field=="Hydrodynamics.densityRepresentation")traced=e.defaulted&&e.used;
        check(traced&&input.fluid.documents[5].effectiveText.find("densityRepresentation cellAverageProjection")!=std::string::npos,
            "conservative default missing from actual trace/effective dictionary");
        const auto effective=root/(rz?"rz_default_conservative_effective":"xy_default_conservative_effective");fs::create_directories(effective);
        for(std::size_t i=0;i<names.size();++i)write(effective/names[i],input.fluid.documents[i].effectiveText);
        const auto roundtrip=readVofElectrohydrodynamicsInput(effective.string());
        check(roundtrip.fluid.vof.momentumDensityRepresentation==input.fluid.vof.momentumDensityRepresentation,
            "conservative default did not survive effective roundtrip");
        reject(root,Fixture(rz),rz?"rz_sharp_momentum_reference":"xy_sharp_momentum_reference","densityRepresentation",
            [](Fixture& bad){replace(bad.text[5],"cellAverageProjection","sharpGeometry");});
    }
}

void run(const fs::path& requested)
{
    const auto root=requested/("reader_"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));fs::create_directories(root);
    positive(root,false);positive(root,true);Fixture xy,rz(true);seriesProfiles(root);sharpAndCoupling(root,xy);meshPhasePolicy(root);geometryConstantsSchema(root);conservativeMomentumSchema(root);
    reject(root,xy,"wrong_app","Run.app",[](Fixture& f){replace(f.text[0],"VofElectrohydrodynamics","Hydrodynamics");});
    reject(root,xy,"hg","Interface.method",[](Fixture& f){replace(f.text[3],"method vof","method hermiteGals");});
    reject(root,xy,"sevof","Interface.method",[](Fixture& f){replace(f.text[3],"method vof","method sevof");});
    reject(root,xy,"old_electric","Electric.method",[](Fixture& f){replace(f.text[6],"vofFiniteVolume","HHO");});
    reject(root,xy,"old_surface_charge","Electric.Initial.surfaceChargeDensity",[](Fixture& f){replace(f.text[6],"chargeDensity","surfaceChargeDensity");});
    reject(root,xy,"timed_charge","Electric.Initial.chargeDensity",[](Fixture& f){replace(f.text[6],"chargeScale*(x-.5)","chargeScale*(x-.5)+t");});
    reject(root,xy,"negative_conductivity","conductivity",[](Fixture& f){replace(f.text[4],"conductivity .25","conductivity -.25");});
    reject(root,xy,"zero_epsilon","permittivity",[](Fixture& f){replace(f.text[4],"permittivity 4","permittivity 0");});
    reject(root,xy,"neumann_flux","value",[](Fixture& f){replace(f.text[6],"N { type Neumann; value 0","N { type Neumann; value 1");});
    reject(root,rz,"axis_dirichlet","Boundary.W.type",[](Fixture& f){replace(f.text[6],"W { type Neumann","W { type Dirichlet");});
    reject(root,xy,"no_anchor","Boundary",[](Fixture& f){f.text[6]="Electric { method vofFiniteVolume; Initial { chargeDensity \"chargeScale*(x-.5)\"; } Boundary { W { type Neumann; value 0; } S { type Neumann; value 0; } E { type Neumann; value 0; } N { type Neumann; value 0; } } }";});
    reject(root,xy,"linear_looser","relativeLinearTolerance",[](Fixture& f){replace(f.text[6],"Electric {","Electric { Solver { relativeLinearTolerance 1e-10; }");});
    reject(root,xy,"charge_looser","relativeConservationTolerance",[](Fixture& f){replace(f.text[6],"Electric {","Electric { Solver { relativeConservationTolerance 1e-9; }");});
    reject(root,xy,"unknown_electric","Electric.futureFlag",[](Fixture& f){replace(f.text[6],"Electric {","Electric { futureFlag true;");});
    {const auto path=xy.save(root,"missing_electric");fs::remove(path/"electric.dict");bool caught=false;
        try{readVofElectrohydrodynamicsInput(path.string());}catch(const std::exception& e){caught=std::string(e.what()).find("electric.dict")!=std::string::npos;}
        check(caught,"missing seventh input did not reject at its actual source");++guards;}
    {const auto path=xy.save(root,"extra_eighth");write(path/"extra.dict","Unused {}");bool caught=false;
        try{readVofElectrohydrodynamicsInput(path.string());}catch(const std::exception& e){caught=std::string(e.what()).find("extra.dict")!=std::string::npos;}
        check(caught,"eighth unused document was silently accepted");++guards;}
    std::cout<<"VOF-EHD actual reader: xy/rz/effective-roundtrip, guards="<<guards<<" checks="<<checks<<" fixtures="<<root.string()<<'\n';
}
}
namespace VofEhdOpenInputContract
{
namespace fs=std::filesystem;
using VofEhdInputContract::Fixture;
using VofEhdInputContract::replace;
using VofEhdInputContract::write;
using Patch=HydrodynamicsInput::VofBoundaryPatch;
using Mode=LAMP::Interface::GeometricVofAdvectorSlab::Boundary::Mode;
unsigned checks=0,guards=0;
void check(bool ok,const char* why){++checks;if(!ok)throw std::runtime_error(why);}
Fixture open(bool rz=false)
{
    Fixture fixture(rz);
    replace(fixture.text[5],"Boundary { type impermeable; }",
        "Boundary { type prescribedFlow; Patches { feed { side S; interval (0 .5); type prescribedVelocity; phase liquid; velocity (.03125 .125); } out { side N; interval (0 1); type pressureOutlet; phase gas; pressure 2; } } }");
    return fixture;
}
void density(Fixture& fixture,const std::string& patch,const std::string& value)
{replace(fixture.text[6],"Electric {","Electric { Advection { Inlets { "+patch+" { chargeDensity "+value+"; } } }");}
template<class Action> void rejectOpen(const fs::path& root,Fixture fixture,const std::string& name,const char* field,Action change)
{
    change(fixture);const auto path=fixture.save(root,name);bool caught=false;
    try{readVofElectrohydrodynamicsInput(path.string());}catch(const std::exception& e){caught=true;check(std::string(e.what()).find(field)!=std::string::npos,"open reader guard lost actual source field");}
    check(caught,"invalid open VOF-EHD input accepted");++guards;
}
void positive(const fs::path& root,bool rz,const std::string& name,const char* supplied,double expected)
{
    auto fixture=open(rz);if(supplied)density(fixture,"feed",supplied);
    const auto path=fixture.save(root,name);const auto input=readVofElectrohydrodynamicsInput(path.string());
    const auto& patches=input.fluid.vofBoundaryPatches;
    check(patches.size()==2&&input.inletChargeDensityByPatch.size()==2,"open config order/charge mapping missing");
    const auto& inlet=patches[0];const auto& outlet=patches[1];
    check(inlet.name=="feed"&&inlet.side==Patch::Side::S&&inlet.kind==Patch::Kind::PrescribedVelocity
        &&inlet.interval.x==0&&inlet.interval.y==.5&&inlet.velocity.x==.03125&&inlet.velocity.y==.125
        &&inlet.reservoir==Mode::PureLiquidReservoir,"physical complete inlet vector/interval/phase changed");
    check(outlet.name=="out"&&outlet.side==Patch::Side::N&&outlet.kind==Patch::Kind::PressureOutlet
        &&outlet.interval.x==0&&outlet.interval.y==1&&outlet.pressure==2&&outlet.reservoir==Mode::PureGasReservoir,
        "physical free outlet/pressure/phase changed");
    check(input.inletChargeDensityByPatch[0]==expected&&input.inletChargeDensityByPatch[1]==0,
        "signed/default charge density is not in actual config-patch order");
    check(input.fluid.initialVelocity.x==0&&input.fluid.initialVelocity.y==0&&input.fluid.vof.surfaceTension==0,
        "open reader changed physical initial seed or first-scope sigma");
    bool trace=false;
    for(const auto& entry:input.fluid.trace.entries)if(entry.field=="Electric.Advection.Inlets.feed.chargeDensity")
        trace=entry.used&&entry.defaulted==!supplied&&entry.units=="code_charge/code_volume"
            &&fs::path(entry.file).filename()=="electric.dict";
    check(trace,"incoming signed/default0 charge lacks actual electric.dict provenance");
    const auto effective=root/(name+"_effective");fs::create_directories(effective);
    for(std::size_t i=0;i<VofEhdInputContract::names.size();++i){
        check(input.fluid.documents[i].sourceText==fixture.text[i],"open reader rewrote source text");
        write(effective/VofEhdInputContract::names[i],input.fluid.documents[i].effectiveText);}
    const auto again=readVofElectrohydrodynamicsInput(effective.string());
    check(again.fluid.documents.size()==7&&again.fluid.radial==rz&&again.fluid.vofBoundaryPatches.size()==2
        &&again.fluid.vofBoundaryPatches[0].name==inlet.name&&again.fluid.vofBoundaryPatches[0].side==inlet.side
        &&again.fluid.vofBoundaryPatches[0].kind==inlet.kind&&again.fluid.vofBoundaryPatches[0].interval.x==inlet.interval.x
        &&again.fluid.vofBoundaryPatches[0].interval.y==inlet.interval.y&&again.fluid.vofBoundaryPatches[0].velocity.x==inlet.velocity.x
        &&again.fluid.vofBoundaryPatches[0].velocity.y==inlet.velocity.y&&again.fluid.vofBoundaryPatches[0].reservoir==inlet.reservoir
        &&again.fluid.vofBoundaryPatches[1].pressure==2&&again.inletChargeDensityByPatch==input.inletChargeDensityByPatch,
        "seven-document open effective input failed typed roundtrip");
}
void run(const fs::path& requested)
{
    const auto root=requested/("open_reader_"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));fs::create_directories(root);
    positive(root,false,"xy_negative","-.25",-.25);positive(root,true,"rz_default",nullptr,0);
    positive(root,false,"xy_positive",".5",.5);positive(root,false,"xy_explicit_zero","0",0);
    const auto xy=open(),rz=open(true);
    {Fixture fixture;const auto input=readVofElectrohydrodynamicsInput(fixture.save(root,"closed_unchanged").string());
        check(input.fluid.vofBoundaryPatches.empty()&&input.inletChargeDensityByPatch.empty(),"closed input acquired open state");}
    rejectOpen(root,xy,"sigma","surfaceTension",[](Fixture& f){replace(f.text[4],"surfaceTension 0","surfaceTension 1");});
    rejectOpen(root,xy,"unaligned","Boundary.Patches.feed.interval",[](Fixture& f){replace(f.text[5],"interval (0 .5)","interval (0 .3)");});
    rejectOpen(root,xy,"reversed","Boundary.Patches.feed.interval",[](Fixture& f){replace(f.text[5],"interval (0 .5)","interval (.5 0)");});
    rejectOpen(root,xy,"overlap","Boundary.Patches.out.interval",[](Fixture& f){replace(f.text[5],"out { side N","out { side S");});
    rejectOpen(root,rz,"true_axis","Boundary.Patches.feed.side",[](Fixture& f){replace(f.text[5],"feed { side S","feed { side W");});
    rejectOpen(root,xy,"zero_normal","Boundary.Patches.feed.velocity",[](Fixture& f){replace(f.text[5],"velocity (.03125 .125)","velocity (.03125 0)");});
    rejectOpen(root,xy,"no_inlet","Boundary.Patches",[](Fixture& f){replace(f.text[5],"velocity (.03125 .125)","velocity (.03125 -.125)");});
    rejectOpen(root,xy,"no_outlet","Boundary.Patches",[](Fixture& f){replace(f.text[5],"type pressureOutlet; phase gas; pressure 2","type prescribedVelocity; phase gas; velocity (0 -.125)");});
    rejectOpen(root,xy,"pressure_payload","Boundary.Patches.out.velocity",[](Fixture& f){replace(f.text[5],"pressure 2;","pressure 2; velocity (0 0);");});
    rejectOpen(root,xy,"unknown_injection","Advection.Inlets.wall",[](Fixture& f){density(f,"wall","1");});
    rejectOpen(root,xy,"free_outlet_injection","Advection.Inlets.out",[](Fixture& f){density(f,"out","1");});
    rejectOpen(root,xy,"outgoing_injection","Advection.Inlets.out",[](Fixture& f){replace(f.text[5],"type pressureOutlet; phase gas; pressure 2","type prescribedVelocity; phase gas; velocity (0 .125)");density(f,"out","1");});
    rejectOpen(root,xy,"nonfinite_injection","Advection.Inlets.feed.chargeDensity",[](Fixture& f){density(f,"feed","nan");});
    rejectOpen(root,Fixture{},"closed_advection","Electric.Advection",[](Fixture& f){density(f,"feed","1");});
    // A legal six-document Hydro fixture must reach the actual Boundary.type
    // choice, rather than failing on the seventh dictionary or electric material keys.
    {auto fixture=xy;replace(fixture.text[0],"VofElectrohydrodynamics","Hydrodynamics");
        replace(fixture.text[4]," permittivity 2; conductivity 0;","");replace(fixture.text[4]," permittivity 4; conductivity .25;","");
        replace(fixture.text[1]," chargeScale 2;","");const auto path=fixture.save(root,"hydro_optin_rejected");fs::remove(path/"electric.dict");bool caught=false;
        try{readHydrodynamicsInput(path.string());}catch(const std::exception& e){caught=std::string(e.what()).find("Hydrodynamics.Boundary.type")!=std::string::npos;}
        check(caught,"standalone Hydro silently opted into open VOF-EHD input");++guards;}
    std::cout<<"VOF-EHD open reader: xy/rz signed/default0/effective-roundtrip, guards="<<guards<<" checks="<<checks<<" fixtures="<<root.string()<<'\n';
}
}
namespace VofEhdJointContract
{
namespace fs=std::filesystem;
namespace H=LAMP::Hydrodynamics;
namespace I=LAMP::Interface;
using Time=LAMP::TimeLayer;
using Wall=HydrodynamicsTransaction::Clock;
unsigned checks=0,guards=0,acceptedSteps=0,changedSteps=0;
bool response=false,workResponse=false;
void check(bool ok,const char* why){++checks;if(!ok)throw std::runtime_error(why);}
bool sameClock(Time a,Time b)
{return a.step==b.step&&a.begin==b.begin&&a.end==b.end&&a.position==b.position;}
std::string read(const fs::path& p)
{
    std::ifstream in(p,std::ios::binary);check(bool(in),"joint contract cannot read actual output prefix");
    return {std::istreambuf_iterator<char>(in),std::istreambuf_iterator<char>()};
}
template<class T> void append(std::string& out,const T& value)
{
    static_assert(std::is_trivially_copyable<T>::value,"raw snapshot requires scalar value");
    out.append(reinterpret_cast<const char*>(&value),sizeof(value));
}
void append(std::string& out,const std::vector<double>& values)
{
    const auto size=values.size();append(out,size);
    if(size)out.append(reinterpret_cast<const char*>(values.data()),size*sizeof(double));
}
struct Snapshot
{
    const void *fluid=nullptr,*charge=nullptr,*field=nullptr;
    Time clock;
    std::uint64_t geometry=0,chargeId=0;
    std::string payload,manifest,frames;
    Snapshot(const QuadtreeMesh& mesh,const VofElectrohydrodynamicsTransaction& tx,const fs::path& run)
    {
        const auto* h=tx.accepted();const auto* q=tx.acceptedCharge();const auto* e=tx.acceptedField();
        check(h&&q&&e,"joint accepted getter is missing");fluid=h;charge=q;field=e;clock=h->timeLayer();
        const auto& pub=h->phaseGeometry().vofPublication();geometry=pub.identity().value();chargeId=q->identity();
        append(payload,h->values().normalVelocity());append(payload,h->values().pressure());
        append(payload,q->integratedCharge());append(payload,e->potential());
        for(const auto& c:pub.plicCells()){
            for(const auto& value:c.inventory)append(payload,value);
            for(const auto& value:c.measured)append(payload,value);
            append(payload,c.normal.x);append(payload,c.normal.y);append(payload,c.offset);append(payload,c.resolved);
        }
        append(payload,mesh.cellIndexNum());
        for(int ci=0;ci<mesh.cellIndexNum();++ci){const auto label=mesh.cellIndexToLabel[ci];
            const auto& cell=mesh.cells(label);append(payload,label);append(payload,cell.level());
            for(int vertex:cell.cellToVert){const auto point=mesh.verts(vertex).coord();append(payload,point.x);append(payload,point.y);}}
        manifest=read(run/"output/meta/manifest.json");frames=read(run/"output/meta/frames.csv");
    }
    void unchanged(const QuadtreeMesh& mesh,const VofElectrohydrodynamicsTransaction& tx,const fs::path& run) const
    {
        const Snapshot now(mesh,tx,run);
        check(fluid==now.fluid&&charge==now.charge&&field==now.field&&sameClock(clock,now.clock)
            &&geometry==now.geometry&&chargeId==now.chargeId&&payload==now.payload,
            "failed/discarded transaction changed accepted mesh/VOF/U/p/Q/phi bytes or owner/clock");
        check(manifest==now.manifest&&frames==now.frames,"uncommitted transaction advanced actual output prefix");
    }
};
template<class Action> void rejects(Action&& action,const char* why)
{bool caught=false;try{action();}catch(const std::exception&){caught=true;}check(caught,why);++guards;}
void validate(const QuadtreeMesh& mesh,const VofElectrohydrodynamicsTransaction& tx,Time expected)
{
    const auto* h=tx.accepted();const auto* q=tx.acceptedCharge();const auto* e=tx.acceptedField();
    check(h&&q&&e,"joint accepted fluid/charge/field incomplete");h->validateFor(mesh);
    const auto& pub=h->phaseGeometry().vofPublication();q->validateFor(mesh,expected);
    e->validateFor(mesh,pub,e->sourceCharge(),expected);
    check(sameClock(h->timeLayer(),expected)&&sameClock(pub.timeLayer(),expected)&&sameClock(q->timeLayer(),expected)
        &&sameClock(e->geometryTime(),expected)&&sameClock(e->candidateTime(),expected),"joint complete accepted clocks differ");
    check(e->geometry().identity()==pub.identity()&&e->chargeCandidate().sameInstance(*q),"joint geometry/Q candidate owner mismatch");
    if(e->isElectrostaticSnapshot())check(e->effectiveDt()==0&&e->sourceCharge().sameInstance(*q)
        &&sameClock(e->chargeSourceTime(),expected),"actual accepted snapshot replaced source/Q/clock");
    else check(expected.position==Time::Position::End&&e->effectiveDt()>0,"ordinary End diagnostic lacks true Ohmic interval");
    for(double u:h->values().normalVelocity()){check(std::isfinite(u),"nonfinite accepted U");response=response||u!=0;}
    bool positive=false,negative=false;
    for(double value:q->integratedCharge()){check(std::isfinite(value),"nonfinite accepted Q");positive=positive||value>0;negative=negative||value<0;}
    check(positive&&negative,"joint fixture lost actual positive/negative charge");
}
void checkReceipt(const QuadtreeMesh& mesh,const VofElectrohydrodynamicsTransaction& tx,Time previous,Time target)
{
    const auto* r=tx.remapReceipt();check(r&&r->decision&&r->sourceMesh&&r->physicalEnd&&r->physicalStep,"actual regrid/physical-step receipt absent");
    check(sameClock(r->decision->sourceClock,target)&&sameClock(r->physicalEnd->timeLayer(),target),"regrid changed actual physical End clock");
    r->physicalEnd->validateFor(*r->sourceMesh);
    const auto& step=*r->physicalStep;check(sameClock(step.beginState().timeLayer(),previous)
        &&step.candidateState().sameIdentity(*r->physicalEnd),"physical step source/candidate owner mismatch");
    check(step.frozenMidpoint().hasMaxwellLoad()&&!step.report().stages.empty(),"actual AC step has no coeval Maxwell load");
    const auto& stage=step.report().stages.back();check(std::isfinite(stage.maxwellWork),"actual final-stage Maxwell work is nonfinite");
    workResponse=workResponse||stage.maxwellWork!=0;
    const auto* h=tx.accepted();const auto* q=tx.acceptedCharge();const auto* e=tx.acceptedField();
    const auto& pub=h->phaseGeometry().vofPublication();
    if(r->decision->changed){
        ++changedSteps;check(r->interfaceTransfer&&r->fluidTransfer,"changed receipt lacks actual phase/fluid transfers");
        check(r->decision->sourceLeaves==static_cast<std::size_t>(r->sourceMesh->cellIndexNum())
            &&r->decision->targetLeaves==static_cast<std::size_t>(mesh.cellIndexNum())
            &&r->sourceMesh->cellIndexNum()!=mesh.cellIndexNum(),"changed fixture did not produce a real native leaf-count change");
        check(r->interfaceTransfer->sourceIdentity()==r->physicalEnd->phaseGeometry().vofPublication().identity()
            &&r->interfaceTransfer->candidate().identity()==pub.identity(),"same-End phase transfer identities are not actual source/target");
        check(e->isElectrostaticSnapshot()&&e->effectiveDt()==0&&e->chargeCandidate().sameInstance(*q)
            &&sameClock(e->chargeSourceTime(),target),"changed accepted field is not the actual same-End target snapshot");
        check(!h->velocityHistory(),"changed native mesh retained incompatible velocity history");
        unsigned hanging=0;
        for(const auto face:pub.domain().faces(mesh))if(face.negativeCell>=0&&face.positiveCell>=0){
            const int a=mesh.cells(mesh.cellIndexToLabel[face.negativeCell]).level();
            const int b=mesh.cells(mesh.cellIndexToLabel[face.positiveCell]).level();
            check(std::abs(a-b)<=1,"actual target shared native face violates 2:1");hanging+=a!=b;}
        check(hanging>0,"changed fixture has no true 2:1 native face");
    }else{
        check(!r->interfaceTransfer&&!r->fluidTransfer,"no-op fabricated a remap receipt");
        check(!e->isElectrostaticSnapshot()&&sameClock(e->chargeSourceTime(),previous),"unchanged End lost original accepted source clock");
    }
}
void run(const fs::path& system,const fs::path& requested)
{
    const auto input=readVofElectrohydrodynamicsInput(fs::absolute(system).string());const auto& f=input.fluid;
    check(f.vofRegrid&&f.vofRegrid->maxLevel==1&&f.nx==8&&f.ny==8&&f.run.fixedStepCount()==2
        &&!f.run.adaptive&&f.run.writeFields&&f.run.writeEverySteps==1,"joint selector requires the frozen small two-step actual seven-document fixture");
    const auto run=fs::absolute(requested)/(std::string("joint_")+(f.radial?"rz_":"xy_")
        +std::to_string(Wall::now().time_since_epoch().count()));fs::create_directories(run);
    DictionaryMesh d;d.Nx=f.nx;d.Ny=f.ny;d.minCoord=f.domainMin;d.maxCoord=f.domainMax;
    d.maxInitialIterations=0;d.maxCellNum=16384;d.maxLevel_local={f.vofRegrid->maxLevel};
    auto mesh=std::make_unique<QuadtreeMesh>(d);
    // A small public geometry fixture, not a copy of the App's band planner:
    // initialize once at allocated finest resolution, then let actual policy coarsen.
    DictionaryMesh::RefineBox fine;fine.level=f.vofRegrid->maxLevel;fine.min=f.domainMin;fine.max=f.domainMax;
    mesh->applyRefineBoxes({fine});
    VofElectrohydrodynamicsOutput writer(input,{run});VofElectrohydrodynamicsTransaction tx(input,writer);
    writer.begin();const char* stage="joint test initialization";
    {
        const auto pub=I::GeometricVofInitialization::prepare(*mesh,f.radial,[&](Vector2D point){return f.initialJet(point).value;},f.vofInitialization);
        const auto layout=H::HydroDofLayout::build(*mesh,pub.domain());
        const auto metrics=H::HydroPhysicalDivergence::build(*mesh,H::HydroIntegralDivergence::build(*mesh,layout));
        const auto geometry=H::HydroPhaseGeometry::build(*mesh,metrics,pub);
        tx.initializeVof(*mesh,H::HydroTwoPhaseState::createVof(*mesh,geometry,
            H::HydroNativeValues(layout,std::vector<double>(layout.velocityCount()),std::vector<double>(layout.pressureCount()))),Wall::now(),stage);
    }
    validate(*mesh,tx,{0,0,0,Time::Position::Initial});
    const Snapshot initial(*mesh,tx,run);
    rejects([&]{tx.prepare(*mesh,{2,0,f.run.dt,Time::Position::End},f.step,false,stage);},"noncontinuous target preparation was accepted");
    initial.unchanged(*mesh,tx,run);
    const Time first{1,0,f.run.stepEnd(0,f.run.dt,0),Time::Position::End};
    tx.prepare(*mesh,first,f.step,false,stage);initial.unchanged(*mesh,tx,run);
    check(!tx.remapReceipt(),"uncommitted preparation exposed a committed receipt");
    tx.releaseStep();initial.unchanged(*mesh,tx,run);
    rejects([&]{tx.commit(*mesh);},"discarded valid candidate could still commit");initial.unchanged(*mesh,tx,run);
    for(int index=0;index<2;++index){
        const auto previous=tx.accepted()->timeLayer();const double begin=previous.time();
        const Time target{previous.step+1,begin,f.run.stepEnd(begin,f.run.dt,static_cast<int>(previous.step)),Time::Position::End};
        const Snapshot source(*mesh,tx,run);const auto started=Wall::now();
        tx.prepare(*mesh,target,f.step,false,stage);source.unchanged(*mesh,tx,run);
        tx.commit(*mesh);++acceptedSteps;validate(*mesh,tx,target);checkReceipt(*mesh,tx,previous,target);
        if(index==0){const Snapshot committed(*mesh,tx,run);
            rejects([&]{tx.commit(*mesh);},"joint candidate committed twice");committed.unchanged(*mesh,tx,run);}
        tx.output(*mesh,started,std::chrono::duration<double>(Wall::now()-started).count(),index==1);tx.releaseStep();
    }
    check(acceptedSteps==2&&tx.accepted()->timeLayer().time()==f.run.endTime,"joint final accepted time/count differs");
    check(changedSteps>0,"joint fixture did not exercise actual changed AMR");
    check(response&&workResponse,"joint two-step fixture did not produce nonzero native response and actual Maxwell work");
    writer.complete();
    const auto manifest=read(run/"output/meta/manifest.json");check(manifest.find("\"status\":\"COMPLETED\"")!=std::string::npos,"joint completed writer lacks terminal manifest");
    const auto frames=read(run/"output/meta/frames.csv");check(std::count(frames.begin(),frames.end(),'\n')==4,"joint writer did not publish exactly initial plus two unique accepted frames");
    tx.clear();check(!tx.accepted()&&!tx.acceptedCharge()&&!tx.acceptedField(),"joint clear retained accepted subset");
    std::cout<<"VOF-EHD actual joint transaction: shape="<<(f.radial?"rz":"xy")<<" accepted="<<acceptedSteps
        <<" changed="<<changedSteps<<" guards="<<guards<<" checks="<<checks<<" output="<<run.string()<<'\n';
}
}
namespace VofB1ZeroWiringContract
{
namespace fs=std::filesystem;
namespace H=LAMP::Hydrodynamics;
namespace I=LAMP::Interface;
namespace E=LAMP::Electric;
using Time=LAMP::TimeLayer;
using Wall=HydrodynamicsTransaction::Clock;
using VofEhdJointContract::check;
using VofEhdJointContract::read;
using VofEhdJointContract::append;
using VofEhdJointContract::rejects;
using VofEhdJointContract::sameClock;
struct Snapshot
{
    const void *fluid=nullptr,*surface=nullptr,*field=nullptr;
    Time clock;std::uint64_t geometry=0;std::string payload,manifest,frames;
    Snapshot(const QuadtreeMesh& mesh,const VofElectrohydrodynamicsTransaction& tx,const fs::path& run)
    {
        const auto* h=tx.accepted();const auto* q=tx.acceptedSurfaceCharge();const auto* e=tx.acceptedSharpField();
        check(h&&q&&e&&!tx.acceptedCharge()&&!tx.acceptedField(),"B1 accepted getters mixed surface and volume charge");
        fluid=h;surface=q;field=e;clock=h->timeLayer();geometry=h->phaseGeometry().vofIdentity().value();
        append(payload,h->values().normalVelocity());append(payload,h->values().pressure());
        for(const auto& moment:q->moments())for(double value:moment)append(payload,value);
        for(const auto& cell:h->phaseGeometry().vofPublication().plicCells()){
            for(auto value:cell.inventory)append(payload,value);for(auto value:cell.measured)append(payload,value);
            append(payload,cell.normal.x);append(payload,cell.normal.y);append(payload,cell.offset);append(payload,cell.resolved);}
        append(payload,mesh.cellIndexNum());
        manifest=read(run/"output/meta/manifest.json");frames=read(run/"output/meta/frames.csv");
    }
    void unchanged(const QuadtreeMesh& mesh,const VofElectrohydrodynamicsTransaction& tx,const fs::path& run) const
    {
        const Snapshot now(mesh,tx,run);
        check(fluid==now.fluid&&surface==now.surface&&field==now.field&&sameClock(clock,now.clock)
            &&geometry==now.geometry&&payload==now.payload,"failed/discarded B1 candidate changed accepted owners/clock/VOF/U/p/Q moments");
        check(manifest==now.manifest&&frames==now.frames,"uncommitted B1 candidate changed actual accepted output prefix");
    }
};
void validate(const QuadtreeMesh& mesh,const VofElectrohydrodynamicsTransaction& tx,Time target)
{
    const auto* h=tx.accepted();const auto* q=tx.acceptedSurfaceCharge();const auto* e=tx.acceptedSharpField();
    check(h&&q&&e,"B1 accepted joint bundle is incomplete");h->validateFor(mesh);
    const auto& pub=h->phaseGeometry().vofPublication();const auto& smooth=q->geometry();smooth.validateFor(mesh,pub);q->validateFor(mesh,smooth,target);
    check(sameClock(h->timeLayer(),target)&&sameClock(pub.timeLayer(),target)&&sameClock(smooth.sourceTime(),target)
        &&sameClock(e->candidateTime(),target)&&e->geometry().sameInstance(smooth)&&e->chargeCandidate().sameInstance(*q),"B1 accepted full clocks or actual geometry/charge owners differ");
    if(target.position==Time::Position::Initial)check(e->effectiveDt()==0&&e->sourceCharge().sameInstance(*q)&&!e->advection(),"B1 initial field is not a genuine electrostatic snapshot");
    else check(e->advection()!=nullptr&&e->effectiveDt()==target.time()-target.begin,"B1 End field lost full-interval material/Ohmic receipt");
    for(const auto& moment:q->moments())for(double value:moment)check(value==0,"B1 zero wiring generated a surface moment");
}
void checkStep(const QuadtreeMesh& mesh,const VofElectrohydrodynamicsTransaction& tx,const E::ElectricVofSurfaceCharge& oldQ,Time target,const H::HydroVofAcConfig& config)
{
    const auto* step=tx.preparedStep();check(step&&!step->report().stages.empty(),"B1 did not expose a complete actual physical step");
    const auto& report=step->report();
    check(report.stages.size()==2&&report.maxwellCorrectionChecked&&std::isfinite(report.maxwellLoadChange)&&std::isfinite(report.maxwellWorkChange)
        &&report.maxwellLoadChange<=config.maximumMaxwellLoadChange&&report.maxwellWorkChange<=config.maximumMaxwellWorkChange,"B1 did not apply the actual two-stage electric load/work correction gate");
    check(step->frozenMidpoint().hasMaxwellLoad(),"B1 actual Hydro solve did not consume Maxwell load");
    const auto& load=*step->frozenMidpoint().forces().maxwellLoad;check(load.isSharp(),"B1 actual Hydro solve consumed A stress");
    const auto& field=load.sharpStress().field();const auto mid=step->transport().midpoint().timeLayer();
    load.validateSharpFor(mesh,load.pressureLoad(),step->frozenMidpoint().material().material().midpointGeometry(),load.sharpStress(),mid);
    check(field.sourceCharge().sameInstance(oldQ)&&field.advection()!=nullptr
        &&field.advection()->sourceCharge().sameInstance(oldQ)&&field.advection()->interfaceTransport()!=nullptr
        &&field.advection()->interfaceTransport()->sameInstance(step->transport()),"B1 final Midpoint restarted from a candidate Q or another transport");
    check(sameClock(field.candidateTime(),mid)&&field.geometry().sourceIdentity()==step->transport().midpoint().identity()
        &&field.effectiveDt()==mid.time()-target.begin,"B1 force clock was manufactured from an End diagnostic");
    check(step->transport().smoothSource()&&step->transport().smoothMidpoint()
        &&step->transport().smoothSource()->sameInstance(oldQ.geometry())
        &&step->transport().smoothMidpoint()->sameInstance(field.geometry()),"B1 actual phase flux and electric field lost their identical smooth geometry owners");
    check(std::all_of(load.momentumRhsLoad().begin(),load.momentumRhsLoad().end(),[](double x){return x==0;}),"B1 zero fixture produced a spurious native electric force");
    check(step->report().stages.back().maxwellWork==0,"B1 actual zero RHS produced nonzero physical electric work");
}
void run(const fs::path& system,const fs::path& requested)
{
    const auto input=readVofElectrohydrodynamicsInput(fs::absolute(system).string());const auto& f=input.fluid;
    check(input.method==VofElectrohydrodynamicsInput::Method::VofSharpSurface&&!f.radial&&!f.vofRegrid&&!f.initialInterfaceBand
        &&f.run.fixedStepCount()==2&&!f.run.adaptive&&f.run.writeFields&&f.run.writeEverySteps==1,
        "B1 joint mode requires its actual uniform closed XY two-step fixture");
    check(f.vof.surfaceTension==0&&f.vof.bodyAcceleration.x==0&&f.vof.bodyAcceleration.y==0
        &&f.initialVelocity.x==0&&f.initialVelocity.y==0,"B1 zero wiring fixture includes a mechanical forcing");
    check(input.initialSurfaceChargeDensity({.487,.513})==0&&std::all_of(input.boundaries.begin(),input.boundaries.end(),[](const auto& b){return b.value==0;}),"B1 zero wiring fixture has nonzero electrical input");
    check(requested.is_absolute(),"B1 joint output parent must be an explicit absolute path");
    const auto run=requested/("b1_zero_wiring_"+std::to_string(Wall::now().time_since_epoch().count()));fs::create_directories(run);
    DictionaryMesh d;d.Nx=f.nx;d.Ny=f.ny;d.minCoord=f.domainMin;d.maxCoord=f.domainMax;d.maxInitialIterations=0;d.maxCellNum=4096;
    auto mesh=std::make_unique<QuadtreeMesh>(d);VofElectrohydrodynamicsOutput writer(input,{run});VofElectrohydrodynamicsTransaction tx(input,writer);
    writer.begin();const char* stage="B1 joint actual initialization";
    {
        const auto publication=I::GeometricVofInitialization::prepare(*mesh,false,[&](Vector2D p){return f.initialJet(p).value;},f.vofInitialization);
        const auto layout=H::HydroDofLayout::build(*mesh,publication.domain());
        const auto metrics=H::HydroPhysicalDivergence::build(*mesh,H::HydroIntegralDivergence::build(*mesh,layout));
        tx.initializeVof(*mesh,H::HydroTwoPhaseState::createVof(*mesh,H::HydroPhaseGeometry::build(*mesh,metrics,publication),
            H::HydroNativeValues(layout,std::vector<double>(static_cast<std::size_t>(layout.velocityCount())),std::vector<double>(static_cast<std::size_t>(layout.pressureCount())))),Wall::now(),stage);
    }
    validate(*mesh,tx,{0,0,0,Time::Position::Initial});const Snapshot initial(*mesh,tx,run);
    rejects([&]{tx.prepare(*mesh,{2,0,f.run.dt,Time::Position::End},f.step,false,stage);},"B1 wrong next-step clock was accepted");initial.unchanged(*mesh,tx,run);
    const Time first{1,0,f.run.stepEnd(0,f.run.dt,0),Time::Position::End};
    tx.prepare(*mesh,first,f.step,false,stage);initial.unchanged(*mesh,tx,run);checkStep(*mesh,tx,*tx.acceptedSurfaceCharge(),first,f.vof);
    const auto* proposed=tx.preparedStep();
    rejects([&]{writer.accepted(*mesh,proposed->candidateState(),*tx.acceptedSurfaceCharge(),*tx.acceptedSharpField(),proposed,Wall::now(),0,false);},
        "B1 writer accepted a physical candidate with unrelated old Q/field");initial.unchanged(*mesh,tx,run);
    rejects([&]{tx.output(*mesh,Wall::now(),0,false);},"B1 uncommitted candidate reached accepted output");initial.unchanged(*mesh,tx,run);
    tx.releaseStep();check(!tx.preparedStep(),"B1 release retained a complete candidate");initial.unchanged(*mesh,tx,run);
    rejects([&]{tx.commit(*mesh);},"B1 discarded candidate could still commit");initial.unchanged(*mesh,tx,run);
    for(int index=0;index<2;++index){const auto previous=tx.accepted()->timeLayer();const double begin=previous.time();
        const Time target{previous.step+1,begin,f.run.stepEnd(begin,f.run.dt,static_cast<int>(previous.step)),Time::Position::End};
        const auto oldQ=*tx.acceptedSurfaceCharge();const Snapshot source(*mesh,tx,run);const auto started=Wall::now();
        tx.prepare(*mesh,target,f.step,false,stage);source.unchanged(*mesh,tx,run);checkStep(*mesh,tx,oldQ,target,f.vof);
        tx.commit(*mesh);validate(*mesh,tx,target);check(!tx.remapReceipt(),"B1 same-mesh step fabricated an AMR receipt");
        const auto* endpoint=tx.acceptedSharpField();const auto* step=tx.preparedStep();
        check(endpoint->sourceCharge().sameInstance(oldQ)&&endpoint->advection()->sourceCharge().sameInstance(oldQ)
            &&endpoint->advection()->interfaceTransport()->sameInstance(step->transport()),"B1 End accumulated Midpoint charge instead of restarting accepted Q");
        check(step->transport().smoothEndpoint()&&step->transport().smoothEndpoint()->sameInstance(endpoint->geometry()),
            "B1 End field differs from the actual final phase-flux smooth geometry");
        check(!endpoint->implicitInputCharge().sameInstance(step->frozenMidpoint().forces().maxwellLoad->sharpStress().field().implicitInputCharge()),"B1 reused Midpoint material input as End input");
        check(std::all_of(tx.accepted()->values().normalVelocity().begin(),tx.accepted()->values().normalVelocity().end(),[](double u){return u==0;}),"B1 zero wiring generated a spurious native velocity");
        const Snapshot committed(*mesh,tx,run);rejects([&]{tx.commit(*mesh);},"B1 joint candidate committed twice");committed.unchanged(*mesh,tx,run);
        tx.output(*mesh,started,std::chrono::duration<double>(Wall::now()-started).count(),index==1);tx.releaseStep();
    }
    writer.complete();const auto frames=read(run/"output/meta/frames.csv");
    check(std::count(frames.begin(),frames.end(),'\n')==4,"B1 did not publish exactly Initial plus two accepted frames");
    check(read(run/"output/meta/manifest.json").find("\"status\":\"COMPLETED\"")!=std::string::npos,"B1 actual output did not complete");
    tx.clear();check(!tx.accepted()&&!tx.acceptedSurfaceCharge()&&!tx.acceptedSharpField()&&!tx.preparedStep(),"B1 clear kept a subset alive");
    std::cout<<"B1 ZERO-FORCE actual two-step typed joint/output wiring PASS; charged motion unqualified; output="<<run.string()<<'\n';
}
}
namespace VofEhdOpenJointContract
{
namespace fs=std::filesystem;
namespace H=LAMP::Hydrodynamics;
namespace E=LAMP::Electric;
namespace I=LAMP::Interface;
namespace J=VofEhdJointContract;
using Clock=LAMP::TimeLayer;
using Wall=HydrodynamicsTransaction::Clock;
unsigned checks=0,accepted=0,changed=0,phaseChecks=0,chargeChecks=0,axisFaces=0;
long double peakDiv=0,peakPhase=0,peakCharge=0,peakGauss=0;
bool realInlet=false,realOutlet=false,realSignedQt=false,realOhmic=false;
void check(bool ok,const char* why){++checks;if(!ok)throw std::runtime_error(why);}
void bound(long double error,long double scale,long double tolerance,const char* why)
{
    check(std::isfinite(error)&&std::isfinite(scale)&&scale>=0
        &&(scale==0?error==0:std::abs(error)<=tolerance*scale),why);
}
std::array<long double,2> stocks(const I::GeometricVofPublication& p)
{
    std::array<long double,2> result{};
    for(const auto& cell:p.plicCells())for(int k=0;k<2;++k)result[k]+=cell.inventory[k];
    return result;
}
void nativeFlow(const QuadtreeMesh& mesh,const H::HydroTwoPhaseState& state)
{
    state.validateFor(mesh);const auto* flow=state.vofPrescribedFlow();
    check(flow&&!flow->pressureOutlets.empty(),"F6B open state lacks actual free pressure outlet");
    const auto& layout=state.values().layout();const auto& domain=layout.domain();
    const auto& u=state.values().normalVelocity();const auto& D=state.phaseGeometry().metrics().divergence().matrix();
    for(int p=0;p<D.rows();++p){
        long double residual=0,scale=0;
        for(H::HydroIntegralDivergence::Matrix::InnerIterator it(D,p);it;++it){
            const long double term=static_cast<long double>(it.value())*u[it.col()];residual+=term;scale+=std::abs(term);}
        bound(residual,scale,1e-11L,"F6B accepted Initial/End original native D_I U");
        peakDiv=std::max(peakDiv,scale?std::abs(residual)/scale:0.L);
    }
    long double total=0,activity=0;unsigned fixed=0,outlets=0,axis=0;
    for(int row=0;row<domain.nativeFaceCount();++row){
        const auto f=domain.face(mesh,row);if(f.negativeCell>=0&&f.positiveCell>=0)continue;
        const int dof=layout.velocityDof(row);const auto rule=flow->rulesByNativeRow.at(row);
        if(f.area==0){
            check(domain.axisymmetric()&&f.center.x==0&&rule==I::GeometricVofAdvectorSlab::Boundary::Mode::Closed,"F6B true axis cannot become open");
            if(dof>=0)check(u[dof]==0,"F6B true axis accepted normal nonzero");++axis;continue;
        }
        check(dof>=0,"F6B positive-area exterior native face omitted");
        const long double flux=static_cast<long double>(f.outwardSign)*f.area*u[dof];total+=flux;activity+=std::abs(flux);
        if(rule==I::GeometricVofAdvectorSlab::Boundary::Mode::Closed)check(u[dof]==0,"F6B wall changed accepted normal");
        else if(f.horizontal&&f.negativeCell<0){
            check(rule==I::GeometricVofAdvectorSlab::Boundary::Mode::PureLiquidReservoir
                &&u[dof]==flow->normalVelocity[dof]
                &&(state.timeLayer().position!=Clock::Position::Initial||flow->normalVelocity[dof]==.02)
                &&std::find(flow->prescribedVelocityDofs.begin(),flow->prescribedVelocityDofs.end(),dof)!=flow->prescribedVelocityDofs.end(),
                "F6B actual inlet fixed normal/phase/DOF mapping");
            const auto vector=flow->endpointVelocity(row,state.timeLayer());
            check(vector.x==0&&vector.y==flow->normalVelocity[dof],"F6B native-row full endpoint vector equals actual fixed normal");++fixed;
        }else{
            check(f.horizontal&&f.positiveCell<0&&flux>0
                &&rule==I::GeometricVofAdvectorSlab::Boundary::Mode::PureGasReservoir,"F6B free pressure outlet reversed or wrong phase");
            check(std::find(flow->prescribedVelocityDofs.begin(),flow->prescribedVelocityDofs.end(),dof)==flow->prescribedVelocityDofs.end(),"F6B pressure outlet was fixed");
            bool found=false;for(const auto& port:flow->pressureOutlets)if(port.velocityDof==dof){check(port.pressure==0,"F6B actual physical outlet pressure");found=true;}
            check(found,"F6B free outlet absent from actual Saddle specification");++outlets;
        }
    }
    bound(total,activity,1e-12L,"F6B independent physical side total flow balance");
    check(fixed>0&&outlets>0&&axis>0,"F6B true-axis/inlet/free-outlet roles missing");axisFaces+=axis;
}
void phaseLedger(const QuadtreeMesh& mesh,const I::GeometricVofTransport& transport,bool half)
{
    const auto source=stocks(transport.source());const auto target=stocks(half?transport.midpoint():transport.endpoint());
    std::array<long double,2> outward{},activity{};
    for(const auto& f:transport.faces())if(f.negativeCell<0||f.positiveCell<0){
        const auto& flux=half?f.half:f.full;const int sign=f.negativeCell>=0?1:-1;
        const std::array<long double,2> phase{{flux.liquid,flux.gas}};
        for(int k=0;k<2;++k){outward[k]+=sign*phase[k];activity[k]+=std::abs(phase[k]);}
        const auto actual=transport.source().domain().face(mesh,f.nativeRow);
        if(actual.area==0)check(flux.total==0&&flux.liquid==0&&flux.gas==0,"F6B axis cannot transfer physical phase");
    }
    for(int k=0;k<2;++k){
        const long double scale=std::abs(source[k])+std::abs(target[k])+activity[k];
        const long double error=target[k]-source[k]+outward[k];bound(error,scale,1e-12L,"F6B independent physical phase delta plus actual exterior phase flux");
        peakPhase=std::max(peakPhase,scale?std::abs(error)/scale:0.L);++phaseChecks;
    }
    const auto& ledger=half?transport.report().boundaryHalf:transport.report().boundaryFull;
    check(ledger.enabled&&ledger.physicalRelativeConservation<=1e-12L
        &&ledger.algebraicRelativeConservation<=1e-12L&&ledger.massRelativeConservation<=1e-12L,"F6B original F1 physical/WY/mass gates");
}
void physicalCharge(const QuadtreeMesh& mesh,const E::ElectricVofCharge& source,
    const I::GeometricVofTransport& transport,const E::ElectricVofChargeTransport& moved,const E::ElectricVofField& field)
{
    const auto& adv=moved.endpoint();adv.validateFor(mesh,source,transport,transport.endTime());
    check(moved.sourceCharge().sameInstance(source)&&moved.transport().sameInstance(transport)
        &&adv.boundary()&&adv.boundary()->sourceCharge().sameInstance(source)
        &&adv.boundary()->transport().sameInstance(transport),"F6B actual source-bound F4 final transport receipt");
    field.validateFor(mesh,transport.endpoint(),source,transport.endTime());
    check(!field.isElectrostaticSnapshot()&&field.advection()&&field.advection()->sameInstance(adv)
        &&field.sourceCharge().sameInstance(source)&&field.geometry().identity()==transport.endpoint().identity()
        &&J::sameClock(field.chargeSourceTime(),source.timeLayer())
        &&J::sameClock(field.candidateTime(),transport.endTime()),"F6B physical End field lost actual Qt/source/full clock");
    const double dt=field.effectiveDt();check(dt==transport.endTime().end-transport.endTime().begin&&dt>0,"F6B physical End actual dt");
    const auto& old=source.integratedCharge();const auto& next=field.chargeCandidate().integratedCharge();
    std::vector<long double> qt(old.size()),j(old.size()),d(old.size()),scaleQt(old.size()),scaleJ(old.size()),scaleD(old.size());
    long double exteriorQt=0,exteriorJ=0,exteriorScale=0;unsigned incoming=0,outgoing=0;
    for(const auto& f:adv.faces()){
        for(const auto incident:{std::pair<int,int>{f.negativeCell,1},{f.positiveCell,-1}})if(incident.first>=0){
            qt[incident.first]+=incident.second*static_cast<long double>(f.integratedCharge);scaleQt[incident.first]+=std::abs(f.integratedCharge);}
        const bool boundary=f.negativeCell<0||f.positiveCell<0;
        if(boundary){
            const int sign=f.negativeCell>=0?1:-1;const auto outward=sign*f.totalVolume;
            const double rho=adv.boundary()->inflowChargeDensities().at(f.nativeRow);
            const auto actual=transport.source().domain().face(mesh,f.nativeRow);
            if(outward<0){
                check(rho==-.01,"F6B signed inlet rho mapped by actual configured patch");
                const long double model=static_cast<long double>(rho)*f.totalVolume;
                bound(static_cast<long double>(f.integratedCharge)-model,std::abs(model)+std::abs(f.integratedCharge),1e-10L,"F6B actual inlet signed rho times real T");
                check(f.integratedCharge<0,"F6B negative reservoir charge flux sign");++incoming;realInlet=true;
            }else{check(rho==0,"F6B unused/outlet/axis density payload not zero");if(outward>0){++outgoing;realOutlet=true;}}
            if(actual.area==0)check(f.totalVolume==0&&f.integratedCharge==0,"F6B true axis Qt zero");
            exteriorQt+=sign*static_cast<long double>(f.integratedCharge);exteriorScale+=std::abs(f.integratedCharge);
        }
        realSignedQt=realSignedQt||f.integratedCharge<0;
    }
    for(const auto& f:field.faces()){
        for(const auto incident:{std::pair<int,int>{f.negativeCell,1},{f.positiveCell,-1}})if(incident.first>=0){
            j[incident.first]+=incident.second*static_cast<long double>(f.current);scaleJ[incident.first]+=std::abs(f.current);
            d[incident.first]+=incident.second*static_cast<long double>(f.displacement);scaleD[incident.first]+=std::abs(f.displacement);}
        if(f.negativeCell<0||f.positiveCell<0){const int sign=f.negativeCell>=0?1:-1;exteriorJ+=sign*static_cast<long double>(f.current);exteriorScale+=dt*std::abs(f.current);}
        realOhmic=realOhmic||f.current!=0;
    }
    long double delta=0,totalScale=exteriorScale;
    for(std::size_t ci=0;ci<old.size();++ci){
        const long double error=static_cast<long double>(next[ci])-old[ci]+qt[ci]+dt*j[ci];
        const long double scale=std::abs(next[ci])+std::abs(old[ci])+scaleQt[ci]+dt*scaleJ[ci];
        bound(error,scale,1e-10L,"F6B independent actual physical cell Qend-Qbegin+Qt+dtJ");
        peakCharge=std::max(peakCharge,scale?std::abs(error)/scale:0.L);
        const long double gauss=d[ci]-next[ci],gs=std::abs(next[ci])+scaleD[ci];
        bound(gauss,gs,1e-10L,"F6B independent physical End native displacement Gauss");
        peakGauss=std::max(peakGauss,gs?std::abs(gauss)/gs:0.L);
        delta+=static_cast<long double>(next[ci])-old[ci];totalScale+=std::abs(next[ci])+std::abs(old[ci]);++chargeChecks;
    }
    bound(delta+exteriorQt+dt*exteriorJ,totalScale,1e-10L,"F6B global physical boundary Qt and actual End J");
    check(incoming>0&&outgoing>0&&adv.report().boundaryEnabled&&adv.report().relativeChargeResidual<=1e-10
        &&adv.report().relativeBoundaryModelChargeResidual<=1e-10,"F6B actual signed reservoir/accepted native Qt original gates");
    check(field.report().relativeLinearResidual<=1e-11&&field.report().maximumRelativeGaussResidual<=1e-10,"F6B retained physical End linear/Gauss gates");
}
void run(const fs::path& requested)
{
    const auto run=fs::absolute(requested)/("open_joint_rz_"+std::to_string(Wall::now().time_since_epoch().count()));fs::create_directories(run);
    VofEhdInputContract::Fixture fixture(true);
    fixture.text[0]="Run { app VofElectrohydrodynamics; } Parallel { threads 1; } Time { units code; dt .00001; endTime .00002; } Output { Fields { enabled true; writeEverySteps 1; } }";
    fixture.text[1]="Geometry { coordinates rz; Domain { min (0 0); max (1 1); } Phases { negative 7; positive 3; } Initial { phi \"z-.5\"; } }";
    fixture.setConstants("chargeScale .02;");
    fixture.text[4]="Materials { Phases { Gas { phaseId 3; density 1; dynamicViscosity 0; permittivity 2; conductivity .1; } Liquid { phaseId 7; density 2; dynamicViscosity 0; permittivity 4; conductivity .25; } } Interfaces { Shared { phases (7 3); surfaceTension 0; } } }";
    fixture.text[5]="Hydrodynamics { method staggered; densityRepresentation cellAverageProjection; capillaryModel heightFunction; Initial { velocity (0 0); } Boundary { type prescribedFlow; Patches { feed { side S; interval (0 1); type prescribedVelocity; phase liquid; velocity (0 .02); } drain { side N; interval (0 1); type pressureOutlet; phase gas; pressure 0; } } } BodyForce { acceleration (0 0); } TimeIntegration { method adaptiveLinear; } }";
    fixture.text[6]="Electric { method vofFiniteVolume; Initial { chargeDensity \"chargeScale*(z-.5)\"; } Boundary { W { type Neumann; value 0; } S { type Dirichlet; value 1; } E { type Neumann; value 0; } N { type Dirichlet; value 0; } } Advection { Inlets { feed { chargeDensity -.01; } } } }";
    const auto system=fixture.save(run,"system");const auto input=readVofElectrohydrodynamicsInput(system.string());const auto& f=input.fluid;
    check(f.radial&&f.vofRegrid&&f.vofRegrid->maxLevel==1&&f.run.fixedStepCount()==2
        &&f.vofBoundaryPatches.size()==2&&input.inletChargeDensityByPatch.size()==2,"F6B actual reader-typed open fixture");
    for(std::size_t i=0;i<f.vofBoundaryPatches.size();++i)check(input.inletChargeDensityByPatch[i]==(f.vofBoundaryPatches[i].name=="feed"?-.01:0),"F6B charge config order follows named physical patches");
    DictionaryMesh dictionary;dictionary.Nx=f.nx;dictionary.Ny=f.ny;dictionary.minCoord=f.domainMin;dictionary.maxCoord=f.domainMax;
    dictionary.maxInitialIterations=0;dictionary.maxCellNum=16384;dictionary.maxLevel_local={1};
    auto mesh=std::make_unique<QuadtreeMesh>(dictionary);
    DictionaryMesh::RefineBox fine;fine.level=1;fine.min=f.domainMin;fine.max=f.domainMax;mesh->applyRefineBoxes({fine});
    VofElectrohydrodynamicsOutput writer(input,{run});VofElectrohydrodynamicsTransaction tx(input,writer);writer.begin();
    const char* stage="F6B actual native Initial projection";
    {
        const auto pub=I::GeometricVofInitialization::prepare(*mesh,true,[&](Vector2D p){return f.initialJet(p).value;},f.vofInitialization);
        const auto layout=H::HydroDofLayout::build(*mesh,pub.domain());
        const auto metrics=H::HydroPhysicalDivergence::build(*mesh,H::HydroIntegralDivergence::build(*mesh,layout));
        const auto geometry=H::HydroPhaseGeometry::build(*mesh,metrics,pub);
        const H::HydroNativeValues zeroSeed(layout,std::vector<double>(layout.velocityCount()),std::vector<double>(layout.pressureCount()));
        auto initial=tx.makeInitialVofState(*mesh,geometry,zeroSeed);
        check(initial.vofPrescribedFlow()&&initial.timeLayer().position==Clock::Position::Initial,"F6B Initial hook returned an unsupported seed");
        nativeFlow(*mesh,initial);
        // One necessary request-mode probe on this same real Initial owner:
        // positive mu without no-slip must retain Full/canonical raw for K.
        auto positiveInput=input;positiveInput.fluid.vof.liquidViscosity=.1;positiveInput.fluid.vof.gasViscosity=.1;
        VofElectrohydrodynamicsOutput positiveWriter(positiveInput,{run/"positive_mu_probe"});
        VofElectrohydrodynamicsTransaction positive(positiveInput,positiveWriter);positiveWriter.begin();
        auto positiveInitial=positive.makeInitialVofState(*mesh,geometry,zeroSeed);
        check(positiveInitial.vofPrescribedFlow()&&positiveInitial.vofPrescribedFlow()->stationaryNoSlipVelocityDofs.empty(),"positive-mu request probe invented no-slip");
        positive.initializeVof(*mesh,std::move(positiveInitial),Wall::now(),stage);
        positive.prepare(*mesh,{1,0,f.run.dt,Clock::Position::End},positiveInput.fluid.step,false,stage);
        const auto* positiveStep=positive.preparedStep();check(positiveStep&&positiveStep->frozenMidpoint().hasMaxwellLoad(),"positive-mu actual Main probe has no Frozen/Maxwell");
        const auto& positiveFrozen=positiveStep->frozenMidpoint();const auto& full=positiveFrozen.forces().maxwellLoad->geometry();
        check(full.isRt0GreenCompatible()&&full.isStressCompatible()&&full.hasCanonicalRawStrain()
            &&positiveFrozen.viscosity().liquidViscosity()==.1&&positiveFrozen.viscosity().gasViscosity()==.1
            &&positiveFrozen.viscosity().viscosity().matrix().nonZeros()>0,"no-wall positive-mu Main skipped Full/canonical strain or real viscosity");
        positive.releaseStep();check(!positive.preparedStep(),"positive-mu request release retained candidate");positive.clear();
        std::cout<<"RT0_GREEN_MAIN_POSITIVE_MU_FULL_PASS prepares=1 releases=1 same_Initial_phase_layout_seed=1 original_open_two_steps_unchanged=1\n";
        tx.initializeVof(*mesh,std::move(initial),Wall::now(),stage);
    }
    J::validate(*mesh,tx,{0,0,0,Clock::Position::Initial});
    check(!tx.physicalChargeTransport()&&!tx.physicalEndField(),"F6B initial state fabricated physical step receipts");
    const J::Snapshot initial(*mesh,tx,run);
    J::rejects([&]{tx.prepare(*mesh,{2,0,f.run.dt,Clock::Position::End},f.step,false,stage);},"F6B noncontinuous target accepted");
    initial.unchanged(*mesh,tx,run);
    const Clock first{1,0,f.run.dt,Clock::Position::End};
    tx.prepare(*mesh,first,f.step,false,stage);initial.unchanged(*mesh,tx,run);
    check(tx.physicalChargeTransport()&&tx.physicalEndField()&&!tx.remapReceipt(),"F6B uncommitted final physical receipts absent or exposed committed regrid");
    tx.releaseStep();check(!tx.physicalChargeTransport()&&!tx.physicalEndField(),"F6B release retained old physical receipts");initial.unchanged(*mesh,tx,run);
    J::rejects([&]{tx.commit(*mesh);},"F6B discarded candidate committed");initial.unchanged(*mesh,tx,run);
    for(int index=0;index<2;++index){
        const auto previous=tx.accepted()->timeLayer();const auto sourceQ=*tx.acceptedCharge();
        const J::Snapshot before(*mesh,tx,run);
        const Clock target{previous.step+1,previous.time(),f.run.stepEnd(previous.time(),f.run.dt,static_cast<int>(previous.step)),Clock::Position::End};
        const auto start=Wall::now();tx.prepare(*mesh,target,f.step,false,stage);before.unchanged(*mesh,tx,run);
        check(tx.physicalChargeTransport()&&tx.physicalEndField(),"F6B prepare did not retain actual final transport/End field");
        tx.commit(*mesh);++accepted;J::validate(*mesh,tx,target);nativeFlow(*mesh,*tx.accepted());
        const auto* receipt=tx.remapReceipt();check(receipt&&receipt->sourceMesh&&receipt->physicalEnd&&receipt->physicalStep&&receipt->decision,"F6B actual physical/regrid receipt missing");
        const auto& sourceMesh=*receipt->sourceMesh;const auto& physical=*receipt->physicalEnd;const auto& step=*receipt->physicalStep;
        const auto& actualLoad=step.frozenMidpoint().forces().maxwellLoad;check(actualLoad&&actualLoad->geometry().isRt0GreenCompatible()
            &&!actualLoad->geometry().isStressCompatible()&&!actualLoad->geometry().hasCanonicalRawStrain(),"no-wall zero-mu Main did not request Green-only");
        check(step.beginState().phaseGeometry().vofPublication().identity()==tx.physicalChargeTransport()->transport().sourceIdentity()
            &&J::sameClock(step.beginState().timeLayer(),previous)&&step.candidateState().sameIdentity(physical),"F6B actual physical source/candidate mismatch");
        phaseLedger(sourceMesh,step.transport(),true);phaseLedger(sourceMesh,step.transport(),false);
        physicalCharge(sourceMesh,sourceQ,step.transport(),*tx.physicalChargeTransport(),*tx.physicalEndField());
        check(tx.physicalEndField()->geometry().identity()==physical.phaseGeometry().vofPublication().identity()
            &&J::sameClock(tx.physicalEndField()->candidateTime(),physical.timeLayer()),"F6B physical End field relabeled with target geometry");
        J::checkReceipt(*mesh,tx,previous,target);
        if(receipt->decision->changed){
            ++changed;const auto& remap=*receipt->fluidTransfer;
            const auto& targetLayout=tx.accepted()->values().layout();
            const auto expectedMeans=H::HydroVofMeshTransfer::boundaryNormalMeans(sourceMesh,physical,*mesh,targetLayout.domain());
            const auto& actualFlow=*tx.accepted()->vofPrescribedFlow();
            for(const int dof:actualFlow.prescribedVelocityDofs)
                check(actualFlow.normalVelocity[dof]==expectedMeans.at(targetLayout.nativeFace(dof)),
                    "F6B changed inlet g must equal actual source-polynomial target-face mean exactly");
            check(remap.report().prescribedFlow&&remap.report().pressureOutletNormals>0&&!remap.report().pressureGaugeApplied
                &&remap.report().projectionEnergyIdentityRelativeResidual<=2e-11L,"F6B actual open same-End F5 projection");
            const auto oldStocks=stocks(physical.phaseGeometry().vofPublication());const auto newStocks=stocks(tx.accepted()->phaseGeometry().vofPublication());
            for(int k=0;k<2;++k)bound(newStocks[k]-oldStocks[k],std::abs(newStocks[k])+std::abs(oldStocks[k]),1e-12L,"F6B changed target independent same-End phase stocks");
            long double q0=0,q1=0,scale=0;for(double q:tx.physicalEndField()->chargeCandidate().integratedCharge()){q0+=q;scale+=std::abs(q);}
            for(double q:tx.acceptedCharge()->integratedCharge()){q1+=q;scale+=std::abs(q);}
            bound(q1-q0,scale,1e-10L,"F6B actual same-End signed Q remap");
            check(tx.acceptedField()->isElectrostaticSnapshot()&&tx.acceptedField()->effectiveDt()==0
                &&!tx.acceptedField()->advection(),"F6B target snapshot is not physical End current");
        }else check(!tx.acceptedField()->isElectrostaticSnapshot()
            &&tx.acceptedCharge()->sameInstance(tx.physicalEndField()->chargeCandidate()),"F6B ordinary End lost actual Ohmic candidate");
        tx.output(*mesh,start,std::chrono::duration<double>(Wall::now()-start).count(),index==1);tx.releaseStep();
        check(!tx.physicalChargeTransport()&&!tx.physicalEndField(),"F6B released physical owner still exposed");
    }
    check(accepted==2&&changed>0&&tx.accepted()->timeLayer().time()==f.run.endTime,"F6B two accepted steps with changed AMR missing");
    check(realInlet&&realOutlet&&realSignedQt&&realOhmic&&axisFaces>0,"F6B actual open T/signed Qt/Ohmic J/true axis not exercised");
    writer.complete();check(J::read(run/"output/meta/manifest.json").find("\"status\":\"COMPLETED\"")!=std::string::npos,"F6B native terminal not completed");
    const auto frames=J::read(run/"output/meta/frames.csv");check(std::count(frames.begin(),frames.end(),'\n')==4,"F6B complete index not Initial plus two accepted frames");
    tx.clear();check(!tx.accepted()&&!tx.acceptedCharge()&&!tx.acceptedField(),"F6B clear retained accepted subset");
    std::cout<<"VOF-EHD actual open joint: RZ true_axis accepted="<<accepted<<" changed="<<changed
        <<" phase_checks="<<phaseChecks<<" charge_cells="<<chargeChecks<<" checks="<<checks
        <<" div_rel="<<peakDiv<<" phase_rel="<<peakPhase<<" physical_Qt_J_rel="<<peakCharge
        <<" physical_Gauss_rel="<<peakGauss<<" actual_physical_End_separate_from_target_snapshot output="<<run.string()<<'\n';
}
}

namespace VofEhdOpenInputContract { namespace SouthEven90
{
unsigned readerChecks=0,readerGuards=0;
void checkSouth(bool ok,const char* why){++readerChecks;if(!ok)throw std::runtime_error(why);}
Fixture inputFixture(bool opened)
{
    auto fixture=opened?VofEhdOpenInputContract::open(true):Fixture(true);
    fixture.text[1]="Geometry { coordinates rz; Domain { min (0 0); max (1 1); } Phases { negative 7; positive 3; } Initial { phi \"r^2+z^2-.25^2\"; } }";
    fixture.setConstants("chargeScale .0002;");
    replace(fixture.text[4],"surfaceTension 0","surfaceTension .01");
    replace(fixture.text[4],"conductivity 0;","conductivity .1;");
    replace(fixture.text[3],"Boundary { type none; }","Boundary { type southEven90; interval (0 1); }");
    if(opened){replace(fixture.text[5],"interval (0 .5)","interval (0 .25)");
        replace(fixture.text[5],"velocity (.03125 .125)","velocity (0 .02)");replace(fixture.text[5],"pressure 2;","pressure 0;");}
    fixture.text[6]="Electric { method vofFiniteVolume; Initial { chargeDensity \"chargeScale*(z-.5)\"; } Boundary { W { type Neumann; value 0; } S { type Dirichlet; value .01; } E { type Neumann; value 0; } N { type Dirichlet; value 0; } } }";
    if(opened)VofEhdOpenInputContract::density(fixture,"feed","-.0001");
    return fixture;
}
template<class Action> void rejectSouth(const fs::path& root,Fixture fixture,const char* name,const char* field,Action action)
{
    action(fixture);const auto path=fixture.save(root,name);bool caught=false;
    try{readVofElectrohydrodynamicsInput(path.string());}catch(const std::exception& e){caught=true;checkSouth(std::string(e.what()).find(field)!=std::string::npos,"southEven90 reader guard lost actual field provenance");}
    checkSouth(caught,"invalid southEven90 reader input accepted");++readerGuards;
}
void positiveSouth(const fs::path& root,bool opened)
{
    auto fixture=inputFixture(opened);
    if(!opened)replace(fixture.text[3],"interval (0 1)","interval (.125 .875)");
    const auto path=fixture.save(root,opened?"open_south_even90":"closed_south_even90");
    const auto input=readVofElectrohydrodynamicsInput(path.string());const auto& f=input.fluid;
    const Vector2D expected=opened?Vector2D{0,1}:Vector2D{.125,.875};
    checkSouth(f.radial&&f.vof.surfaceTension==.01&&f.vofSouthEven90Interval
        &&f.vofSouthEven90Interval->x==expected.x&&f.vofSouthEven90Interval->y==expected.y,
        "explicit physical southEven90 interval or positive sigma changed");
    checkSouth(f.vofSymmetryPlanes.empty(),"southEven90 was silently converted to the old symmetry model");
    if(opened)checkSouth(f.vofBoundaryPatches.size()==2&&f.vofBoundaryPatches[0].interval.y==.25
        &&f.vofSouthEven90Interval->y==1&&input.inletChargeDensityByPatch[0]==-.0001,
        "geometry interval was inferred from hydraulic inlet patch or signed charge lost");
    else checkSouth(f.vofBoundaryPatches.empty()&&input.inletChargeDensityByPatch.empty(),"closed explicit geometry acquired an open flow");
    unsigned intervalTrace=0,typeTrace=0;
    for(const auto& entry:f.trace.entries){
        if(entry.field=="Interface.Boundary.interval"){
            ++intervalTrace;checkSouth(entry.used&&!entry.defaulted&&entry.units=="code_length"
                &&entry.consumer=="HydroVofStep.southEven90Boundary.actualStageGeometry"
                &&fs::path(entry.file).filename()=="interface.dict"&&entry.line>0,
                "southEven90 interval trace does not identify actual stage geometry/file/units");}
        if(entry.field=="Interface.Boundary.type"){
            ++typeTrace;checkSouth(entry.used&&!entry.defaulted&&entry.units=="dimensionless"
                &&entry.consumer=="HydroVofStep.southEven90Boundary.actualStageGeometry"
                &&fs::path(entry.file).filename()=="interface.dict","southEven90 boundary type trace missing");}}
    checkSouth(intervalTrace==1&&typeTrace==1,"southEven90 has absent or duplicated field ownership");
    const auto jet=f.initialJet({.25,0});
    checkSouth(jet.value==0&&jet.derivative==.5&&jet.partialSecond==0,
        "reader did not retain the generic hemisphere expression/contact first jet");
    const auto effective=root/(opened?"open_effective":"closed_effective");fs::create_directories(effective);
    for(std::size_t i=0;i<VofEhdInputContract::names.size();++i){checkSouth(f.documents[i].sourceText==fixture.text[i],"southEven90 reader changed original document");write(effective/VofEhdInputContract::names[i],f.documents[i].effectiveText);}
    const auto again=readVofElectrohydrodynamicsInput(effective.string());
    checkSouth(again.fluid.documents.size()==7&&again.fluid.vofSouthEven90Interval
        &&again.fluid.vofSouthEven90Interval->x==expected.x&&again.fluid.vofSouthEven90Interval->y==expected.y
        &&again.fluid.vof.surfaceTension==.01&&again.fluid.vofSymmetryPlanes.empty()
        &&again.fluid.vofBoundaryPatches.size()==f.vofBoundaryPatches.size()
        &&again.inletChargeDensityByPatch==input.inletChargeDensityByPatch,
        "southEven90 seven-dictionary effective input failed actual typed roundtrip");
}
void run(const fs::path& requested)
{
    const auto root=requested/("south_even90_reader_"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));fs::create_directories(root);
    positiveSouth(root,true);positiveSouth(root,false);
    // Original none/symmetry do not synthesize geometry or interval provenance.
    for(bool symmetry:{false,true}){Fixture fixture(true);if(symmetry)replace(fixture.text[3],"Boundary { type none; }","Boundary { type symmetry; sides (bottom); }");
        const auto input=readVofElectrohydrodynamicsInput(fixture.save(root,symmetry?"old_symmetry":"old_none").string());
        checkSouth(!input.fluid.vofSouthEven90Interval&&input.fluid.vofSymmetryPlanes.size()==(symmetry?1u:0u),"old none/symmetry geometry behavior changed");
        for(const auto& entry:input.fluid.trace.entries)checkSouth(entry.field!="Interface.Boundary.interval","absent declaration invented interval/default trace");
        for(const auto& entry:input.fluid.trace.entries)if(entry.field=="Interface.Boundary.type")
            checkSouth(entry.consumer=="Interface.VofInitialization.boundary","old none/symmetry trace consumer changed");}
    const auto fixture=inputFixture(true);
    rejectSouth(root,fixture,"unaligned","Interface.Boundary.interval",[](Fixture& f){replace(f.text[3],"interval (0 1)","interval (0 .3)");});
    rejectSouth(root,fixture,"reversed","Interface.Boundary.interval",[](Fixture& f){replace(f.text[3],"interval (0 1)","interval (1 0)");});
    rejectSouth(root,fixture,"empty","Interface.Boundary.interval",[](Fixture& f){replace(f.text[3],"interval (0 1)","interval (.5 .5)");});
    rejectSouth(root,fixture,"outside","Interface.Boundary.interval",[](Fixture& f){replace(f.text[3],"interval (0 1)","interval (-.125 1)");});
    rejectSouth(root,fixture,"nonfinite","Interface.Boundary.interval",[](Fixture& f){replace(f.text[3],"interval (0 1)","interval (nan 1)");});
    rejectSouth(root,fixture,"arity","Interface.Boundary.interval",[](Fixture& f){replace(f.text[3],"interval (0 1)","interval (0 .5 1)");});
    rejectSouth(root,fixture,"missing_interval","interval",[](Fixture& f){replace(f.text[3]," interval (0 1);","");});
    rejectSouth(root,fixture,"extra_side","Interface.Boundary.sides",[](Fixture& f){replace(f.text[3],"type southEven90;","type southEven90; sides (bottom);");});
    rejectSouth(root,fixture,"none_payload","Interface.Boundary.interval",[](Fixture& f){replace(f.text[3],"type southEven90;","type none;");});
    rejectSouth(root,fixture,"missing_geometry","surfaceTension",[](Fixture& f){replace(f.text[3],"Boundary { type southEven90; interval (0 1); }","Boundary { type none; }");});
    rejectSouth(root,fixture,"old_symmetry_positive_sigma","Interface.Boundary",[](Fixture& f){replace(f.text[3],"Boundary { type southEven90; interval (0 1); }","Boundary { type symmetry; sides (bottom); }");});
    {auto f=inputFixture(false);replace(f.text[0],"VofElectrohydrodynamics","Hydrodynamics");
        replace(f.text[4]," permittivity 2; conductivity .1;","");replace(f.text[4]," permittivity 4; conductivity .25;","");
        f.setConstants("");const auto path=f.save(root,"standalone_hydro");fs::remove(path/"electric.dict");bool caught=false;
        try{readHydrodynamicsInput(path.string());}catch(const std::exception& e){caught=std::string(e.what()).find("Interface.Boundary.type")!=std::string::npos;}
        checkSouth(caught,"standalone Hydro silently enabled EHD-only southEven90 stage geometry");++readerGuards;}
    std::cout<<"VOF_EHD_SOUTH_EVEN90_INPUT_PASS,reader_only,open_closed_positive_sigma,effective_trace,old_none_symmetry,guards,"<<readerGuards
        <<",checks,"<<readerChecks<<",no_Step_or_contact_accuracy_claim,fixtures,"<<root.string()<<'\n';
}
}}

// BEGIN F19 EXPLICIT PRESSURE-OUTLET BACKFLOW READER CONTRACT
namespace F19BackflowInputContract
{
namespace fs=std::filesystem;
using Fixture=VofEhdInputContract::Fixture;
using Model=LAMP::Hydrodynamics::HydroVofPressureOutletMomentumModel;
unsigned checks=0,guards=0,positiveCases=0;
void checkF19(bool ok,const char* why){++checks;if(!ok)throw std::runtime_error(why);}
void configureBackflow(Fixture& fixture,const std::string& body)
{
    VofEhdInputContract::replace(fixture.text[5],"Boundary { type prescribedFlow;",
        "Boundary { type prescribedFlow; Backflow { "+body+" }");
}
bool hasBackflowTrace(const HydrodynamicsInput& fluid)
{
    return std::any_of(fluid.trace.entries.begin(),fluid.trace.entries.end(),[](const auto& e){
        return e.field.find("Hydrodynamics.Boundary.Backflow")==0;});
}
void positiveF19(const fs::path& root,bool rz,bool ambient,const std::string& name,Vector2D velocity)
{
    auto fixture=VofEhdOpenInputContract::open(rz);
    const std::string body=ambient?(velocity.x==0&&velocity.y==0?"model ambient; velocity (0 0);":"model ambient; velocity (.0625 -.03125);"):
        "model interiorContinuationZeroShear;";
    configureBackflow(fixture,body);const auto path=fixture.save(root,name);
    const auto input=readVofElectrohydrodynamicsInput(path.string());const auto& option=input.fluid.vofPressureOutletBackflow;
    checkF19(option.has_value()&&option->model==(ambient?Model::Ambient:Model::InteriorContinuationZeroShear),"F19 explicit model did not reach typed input");
    checkF19(option->velocity.has_value()==ambient,"F19 Interior acquired an unused Ambient vector or Ambient lost its vector");
    if(ambient)checkF19(option->velocity->x==velocity.x&&option->velocity->y==velocity.y,"F19 complete signed Ambient vector changed");
    checkF19(input.fluid.vofBoundaryPatches.size()==2&&input.fluid.vofBoundaryPatches[1].pressure==2
        &&input.fluid.vofBoundaryPatches[0].velocity.x==.03125&&input.fluid.vofBoundaryPatches[0].velocity.y==.125,
        "F19 backflow changed prescribed g/physical pressure/full inlet vector");
    bool modelTrace=false,velocityTrace=false;
    for(const auto& e:input.fluid.trace.entries){
        if(e.field=="Hydrodynamics.Boundary.Backflow.model")modelTrace=e.used&&!e.defaulted&&e.units=="dimensionless"
            &&e.consumer=="HydroVofStep.prescribedOpenMomentum"&&fs::path(e.file).filename()=="hydro.dict";
        if(e.field=="Hydrodynamics.Boundary.Backflow.velocity")velocityTrace=e.used&&!e.defaulted&&e.units=="code_velocity"
            &&e.consumer=="HydroVofStep.pressureOutletBackflowMomentum"&&fs::path(e.file).filename()=="hydro.dict";
    }
    checkF19(modelTrace&&velocityTrace==ambient,"F19 model/vector provenance is missing or marks an unused consumer");
    const auto effective=root/(name+"_effective");fs::create_directories(effective);
    for(std::size_t i=0;i<VofEhdInputContract::names.size();++i){
        checkF19(input.fluid.documents[i].sourceText==fixture.text[i],"F19 reader rewrote actual input bytes");
        VofEhdInputContract::write(effective/VofEhdInputContract::names[i],input.fluid.documents[i].effectiveText);
    }
    const auto again=readVofElectrohydrodynamicsInput(effective.string());const auto& second=again.fluid.vofPressureOutletBackflow;
    checkF19(second.has_value()&&second->model==option->model&&second->velocity.has_value()==ambient,"F19 effective roundtrip lost explicit backflow");
    if(ambient)checkF19(second->velocity->x==velocity.x&&second->velocity->y==velocity.y,"F19 effective roundtrip changed Ambient vector");
    checkF19(again.fluid.vofBoundaryPatches[1].pressure==2&&again.fluid.vofBoundaryPatches[0].velocity.y==.125,
        "F19 effective roundtrip altered physical flow settings");++positiveCases;
}
template<class Change> void rejectF19(const fs::path& root,Fixture fixture,const std::string& name,const char* field,Change change)
{
    change(fixture);const auto path=fixture.save(root,name);bool caught=false;
    try{(void)readVofElectrohydrodynamicsInput(path.string());}
    catch(const std::exception& e){caught=true;checkF19(std::string(e.what()).find(field)!=std::string::npos,"F19 rejection lost actual input field");}
    checkF19(caught,"F19 invalid or unused backflow input accepted");++guards;
}
void run(const fs::path& requested)
{
    const auto root=requested/("f19_backflow_reader_"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));fs::create_directories(root);
    positiveF19(root,false,true,"ambient_xy_nonzero",{.0625,-.03125});positiveF19(root,true,true,"ambient_rz_zero",{0,0});
    positiveF19(root,true,false,"interior_rz",{0,0});
    for(bool opened:{false,true}){
        auto fixture=opened?VofEhdOpenInputContract::open():Fixture{};
        const auto input=readVofElectrohydrodynamicsInput(fixture.save(root,opened?"legacy_open":"legacy_closed").string());
        checkF19(!input.fluid.vofPressureOutletBackflow&&!hasBackflowTrace(input.fluid),"F19 omitted option created backflow state or fake active trace");
        checkF19(input.fluid.documents[5].effectiveText.find("Backflow")==std::string::npos,"F19 omitted option appeared active in effective hydro.dict");
    }
    const auto open=VofEhdOpenInputContract::open();
    rejectF19(root,open,"unknown_model","Backflow.model",[](Fixture& f){configureBackflow(f,"model arbitrary;");});
    rejectF19(root,open,"missing_vector","Backflow.velocity",[](Fixture& f){configureBackflow(f,"model ambient;");});
    rejectF19(root,open,"nonfinite_vector","Backflow.velocity",[](Fixture& f){configureBackflow(f,"model ambient; velocity (nan 0);");});
    rejectF19(root,open,"scalar_vector","Backflow.velocity",[](Fixture& f){configureBackflow(f,"model ambient; velocity 0;");});
    rejectF19(root,open,"interior_unused_vector","Backflow.velocity",[](Fixture& f){configureBackflow(f,"model interiorContinuationZeroShear; velocity (0 0);");});
    rejectF19(root,open,"unknown_payload","Backflow.extra",[](Fixture& f){configureBackflow(f,"model ambient; velocity (0 0); extra 1;");});
    rejectF19(root,Fixture{},"closed_unused","Backflow",[](Fixture& f){VofEhdInputContract::replace(f.text[5],"Boundary { type impermeable; }",
        "Boundary { type impermeable; Backflow { model ambient; velocity (0 0); } }");});
    rejectF19(root,open,"no_pressure_port","Backflow",[](Fixture& f){configureBackflow(f,"model ambient; velocity (0 0);");
        VofEhdInputContract::replace(f.text[5],"type pressureOutlet; phase gas; pressure 2","type prescribedVelocity; phase gas; velocity (0 .125)");});
    // A legal six-document standalone Hydro fixture reaches the Backflow
    // guard, rather than failing on electric-only fields or document count.
    {Fixture fixture;VofEhdInputContract::replace(fixture.text[5],"Boundary { type impermeable; }",
        "Boundary { type impermeable; Backflow { model ambient; velocity (0 0); } }");
        VofEhdInputContract::replace(fixture.text[0],"VofElectrohydrodynamics","Hydrodynamics");
        VofEhdInputContract::replace(fixture.text[4]," permittivity 2; conductivity 0;","");
        VofEhdInputContract::replace(fixture.text[4]," permittivity 4; conductivity .25;","");
        VofEhdInputContract::replace(fixture.text[1]," chargeScale 2;","");const auto path=fixture.save(root,"standalone_hydro_rejects_backflow");fs::remove(path/"electric.dict");
        bool caught=false;try{(void)readHydrodynamicsInput(path.string());}
        catch(const std::exception& e){caught=std::string(e.what()).find("Hydrodynamics.Boundary.Backflow")!=std::string::npos;}
        checkF19(caught,"F19 standalone Hydro accepted an unconsumed VofEHD backflow option");++guards;
    }
    std::cout<<"F19_BACKFLOW_INPUT_PASS positive="<<positiveCases<<" guards="<<guards<<" checks="<<checks
        <<" reader_effective_trace_only_no_physical_backflow_claim fixtures="<<root.string()<<'\n';
}
}
// END F19 EXPLICIT PRESSURE-OUTLET BACKFLOW READER CONTRACT

// BEGIN F22 PRESSURE-BACKFLOW INDEPENDENT CHARGE CONTRACT
namespace F22PressureBackflowChargeContract
{
namespace fs=std::filesystem;
namespace H=LAMP::Hydrodynamics;
namespace I=LAMP::Interface;
namespace E=LAMP::Electric;
namespace J=VofEhdJointContract;
namespace O=VofEhdOpenJointContract;
using Fixture=VofEhdInputContract::Fixture;
using Clock=LAMP::TimeLayer;
using Wall=HydrodynamicsTransaction::Clock;
unsigned checks=0,guards=0,readers=0,accepted=0,changed=0,pressureIncoming=0,negativeInitial=0;
void check(bool ok,const char* why){++checks;if(!ok)throw std::runtime_error(why);}
void bound(long double error,long double scale,long double tolerance,const char* why)
{check(std::isfinite(error)&&std::isfinite(scale)&&scale>=0&&(scale==0?error==0:std::abs(error)<=tolerance*scale),why);}
Fixture fixture(bool rz,const char* drainDensity)
{
    Fixture f(rz);
    f.text[0]="Run { app VofElectrohydrodynamics; } Parallel { threads 1; } Time { units code; dt .00001; endTime .00002; } Output { Fields { enabled true; writeEverySteps 1; } }";
    f.text[1]=std::string("Geometry { coordinates ")+(rz?"rz":"xy")+"; Domain { min (0 0); max (1 1); } Phases { negative 7; positive 3; } Initial { phi \""+(rz?"z-.5":"y-.5")+"\"; } }";
    f.setConstants("chargeScale .02;");
    f.text[4]="Materials { Phases { Gas { phaseId 3; density 1; dynamicViscosity 0; permittivity 2; conductivity .1; } Liquid { phaseId 7; density 2; dynamicViscosity 0; permittivity 4; conductivity .25; } } Interfaces { Shared { phases (7 3); surfaceTension 0; } } }";
    f.text[5]="Hydrodynamics { method staggered; densityRepresentation cellAverageProjection; capillaryModel heightFunction; Initial { velocity (0 0); } Boundary { type prescribedFlow; Backflow { model ambient; velocity (0 0); } Patches { feed { side S; interval (0 1); type prescribedVelocity; phase liquid; velocity (0 .02); } drain { side N; interval (0 1); type pressureOutlet; phase gas; pressure 0; } } } BodyForce { acceleration (0 0); } TimeIntegration { method adaptiveLinear; } }";
    f.text[6]=std::string("Electric { method vofFiniteVolume; Initial { chargeDensity \"chargeScale*(")+(rz?"z":"y")+"-.5)\"; } Boundary { W { type Neumann; value 0; } S { type Dirichlet; value 1; } E { type Neumann; value 0; } N { type Dirichlet; value 0; } } Advection { Inlets { feed { chargeDensity -.01; } "+(drainDensity?std::string("drain { chargeDensity ")+drainDensity+"; } ":"")+"} } }";
    return f;
}
void reader(const fs::path& root,bool rz,const char* supplied,double expected,const std::string& name)
{
    const auto source=fixture(rz,supplied);const auto input=readVofElectrohydrodynamicsInput(source.save(root,name).string());
    check(input.fluid.vofPressureOutletBackflow.has_value()&&input.fluid.vofPressureOutletBackflow->model==H::HydroVofPressureOutletMomentumModel::Ambient,
        "F22 explicit fluid backflow declaration missing");
    check(input.fluid.vofBoundaryPatches.size()==2&&input.fluid.vofBoundaryPatches[1].name=="drain"
        &&input.fluid.vofBoundaryPatches[1].kind==HydrodynamicsInput::VofBoundaryPatch::Kind::PressureOutlet
        &&input.fluid.vofBoundaryPatches[1].reservoir==I::GeometricVofAdvectorSlab::Boundary::Mode::PureGasReservoir
        &&input.inletChargeDensityByPatch.size()==2&&input.inletChargeDensityByPatch[0]==-.01&&input.inletChargeDensityByPatch[1]==expected,
        "F22 independent signed/default drain rho must follow physical patch order");
    bool trace=false;
    for(const auto& e:input.fluid.trace.entries)if(e.field=="Electric.Advection.Inlets.drain.chargeDensity")
        trace=e.used&&e.defaulted==!supplied&&e.units=="code_charge/code_volume"&&fs::path(e.file).filename()=="electric.dict"
            &&e.consumer=="ElectricVofChargeBoundary.signedInflowDensity";
    check(trace,"F22 pressure reservoir charge lacks signed/default electric provenance");
    const auto effective=root/(name+"_effective");fs::create_directories(effective);
    for(std::size_t i=0;i<VofEhdInputContract::names.size();++i){
        check(input.fluid.documents[i].sourceText==source.text[i],"F22 reader changed actual source document");
        VofEhdInputContract::write(effective/VofEhdInputContract::names[i],input.fluid.documents[i].effectiveText);}
    const auto again=readVofElectrohydrodynamicsInput(effective.string());
    check(again.inletChargeDensityByPatch==input.inletChargeDensityByPatch&&again.fluid.radial==rz
        &&again.fluid.vofPressureOutletBackflow.has_value()&&again.fluid.vofBoundaryPatches[1].pressure==0,
        "F22 effective roundtrip changed independent charge/pressure declaration");++readers;
}
template<class Change> void reject(const fs::path& root,const std::string& name,const char* field,Change change)
{
    auto f=fixture(false,".03125");change(f);bool caught=false;
    try{(void)readVofElectrohydrodynamicsInput(f.save(root,name).string());}
    catch(const std::exception& e){caught=true;check(std::string(e.what()).find(field)!=std::string::npos,"F22 reader rejection lost source field");}
    check(caught,"F22 invalid or undeclared reservoir charge accepted");++guards;
}
H::HydroNativeValues seed(const QuadtreeMesh& mesh,const H::HydroDofLayout& layout,bool rz)
{
    const auto psi=[rz](Vector2D p){const long double x=p.x,z=p.y;
        return rz?.01L*x*x+.05L*x*x*(1-x*x)*z*z:.02L*x+.05L*x*(1-x)*z*z;};
    std::vector<double> u(layout.velocityCount());
    for(int v=0;v<layout.velocityCount();++v){const auto face=layout.domain().face(mesh,layout.nativeFace(v));
        if(face.area==0){u[v]=0;continue;}
        const auto& edge=face.horizontal?mesh.hEdges(face.edgeLabel):mesh.vEdges(face.edgeLabel);
        auto a=mesh.verts(edge.edgeToVert[0]).coord(),b=mesh.verts(edge.edgeToVert[1]).coord();
        if((face.horizontal?a.x:a.y)>(face.horizontal?b.x:b.y))std::swap(a,b);
        const long double integrated=(face.horizontal?1.L:-1.L)*(psi(b)-psi(a))*(rz?2*std::acos(-1.L):1.L);
        u[v]=static_cast<double>(integrated/face.area);
    }
    return H::HydroNativeValues(layout,u,std::vector<double>(layout.pressureCount()));
}
void native(const QuadtreeMesh& mesh,const H::HydroTwoPhaseState& state,bool initial)
{
    state.validateFor(mesh);const auto* flow=state.vofPrescribedFlow();check(flow&&!flow->pressureOutlets.empty(),"F22 missing typed physical pressure flow");
    const auto& layout=state.values().layout();const auto& u=state.values().normalVelocity();
    const auto& D=state.phaseGeometry().metrics().divergence().matrix();
    for(int row=0;row<D.rows();++row){long double error=0,scale=0;
        for(H::HydroIntegralDivergence::Matrix::InnerIterator it(D,row);it;++it){const long double term=static_cast<long double>(it.value())*u[it.col()];error+=term;scale+=std::abs(term);}
        bound(error,scale,1e-11L,"F22 actual projected native divergence");}
    unsigned incoming=0;
    for(const auto& port:flow->pressureOutlets){check(port.pressure==0,"F22 physical pressure changed");
        const auto face=layout.domain().face(mesh,layout.nativeFace(port.velocityDof));
        check(face.horizontal&&face.positiveCell<0&&std::find(flow->prescribedVelocityDofs.begin(),flow->prescribedVelocityDofs.end(),port.velocityDof)==flow->prescribedVelocityDofs.end(),
            "F22 pressure unknown was changed to prescribed velocity");
        if(u[port.velocityDof]<0)++incoming;}
    if(initial){check(incoming>0,"F22 public Initial M projection lost actual negative-begin pressure normals");negativeInitial+=incoming;}
}
void charge(const QuadtreeMesh& mesh,const E::ElectricVofCharge& source,const E::ElectricVofChargeTransport& moved,const E::ElectricVofField& field,double drainRho)
{
    const auto& t=moved.transport();const auto& domain=t.source().domain();
    check(moved.sourceCharge().sameInstance(source),"F22 physical charge transport has foreign accepted source");
    unsigned actualPressureIn=0;
    for(bool half:{true,false}){const auto& adv=half?moved.midpoint():moved.endpoint();adv.validateFor(mesh,source,t,adv.timeLayer());
        check(adv.boundary()&&adv.boundary()->sourceCharge().sameInstance(source)&&adv.boundary()->transport().sameInstance(t),"F22 source-bound electric reservoir owner lost");
        for(const auto& f:adv.faces())if(f.negativeCell<0||f.positiveCell<0){
            const auto actual=domain.face(mesh,f.nativeRow);const int sign=f.negativeCell>=0?1:-1;
            const auto& tf=t.faces().at(f.nativeRow);const auto& volume=half?tf.half:tf.full;
            const long double in=sign>0?volume.negativeTotal:volume.positiveTotal,out=sign>0?volume.positiveTotal:volume.negativeTotal;
            const long double fullIn=sign>0?tf.full.negativeTotal:tf.full.positiveTotal;
            const bool pressure=actual.horizontal&&actual.positiveCell<0;
            const double rho=adv.boundary()->inflowChargeDensities().at(f.nativeRow),expected=fullIn>0?(pressure?drainRho:-.01):0;
            check(rho==expected,"F22 native-row reservoir rho must use full incoming gross and independent configured patch charge");
            if(actual.area==0)check(in==0&&out==0&&rho==0&&f.integratedCharge==0,"F22 true axis carried volume or charge");
            if(in>0){
                const int inside=f.negativeCell>=0?f.negativeCell:f.positiveCell;
                const double capacity=domain.cellVolumes().at(domain.cellRowByIndex().at(inside));
                const long double outgoing=static_cast<long double>(source.integratedCharge().at(inside))*out/capacity;
                const long double incoming=static_cast<long double>(rho)*in;
                const long double model=sign*(outgoing-incoming);
                bound(static_cast<long double>(f.integratedCharge)-model,std::abs(f.integratedCharge)+std::abs(outgoing)+std::abs(incoming),1e-10L,
                    "F22 actual signed Qt must contain configured rho times actual incoming gross");
                if(pressure&&!half){check(rho==drainRho&&rho!=0&&f.integratedCharge!=0,"F22 actual pressure-backflow signed Q absent");++actualPressureIn;}
            }
        }
        check(adv.report().boundaryEnabled&&adv.report().relativeChargeResidual<=1e-10&&adv.report().relativeBoundaryModelChargeResidual<=1e-10,
            "F22 original signed Q/shared Qt/boundary model gates");
    }
    check(actualPressureIn>0,"F22 real transaction did not retain pressure-port incoming gross");pressureIncoming+=actualPressureIn;
    const auto& adv=moved.endpoint();field.validateFor(mesh,t.endpoint(),source,t.endTime());
    check(!field.isElectrostaticSnapshot()&&field.advection()&&field.advection()->sameInstance(adv)
        &&field.sourceCharge().sameInstance(source)&&field.geometry().identity()==t.endpoint().identity(),"F22 physical End field owner changed");
    const double dt=field.effectiveDt();check(dt==t.endTime().end-t.endTime().begin&&dt>0,"F22 physical End dt changed");
    std::vector<long double> qt(source.integratedCharge().size()),current(qt.size()),displacement(qt.size()),activity(qt.size()),gaussActivity(qt.size());
    for(const auto& f:adv.faces())for(const auto side:{std::pair<int,int>{f.negativeCell,1},{f.positiveCell,-1}})if(side.first>=0){qt[side.first]+=side.second*static_cast<long double>(f.integratedCharge);activity[side.first]+=std::abs(f.integratedCharge);}
    for(const auto& f:field.faces())for(const auto side:{std::pair<int,int>{f.negativeCell,1},{f.positiveCell,-1}})if(side.first>=0){
        current[side.first]+=side.second*static_cast<long double>(f.current);activity[side.first]+=dt*std::abs(f.current);
        displacement[side.first]+=side.second*static_cast<long double>(f.displacement);gaussActivity[side.first]+=std::abs(f.displacement);}
    for(std::size_t c=0;c<qt.size();++c){const double old=source.integratedCharge()[c],next=field.chargeCandidate().integratedCharge()[c];
        bound(static_cast<long double>(next)-old+qt[c]+dt*current[c],std::abs(old)+std::abs(next)+activity[c],1e-10L,"F22 independent actual cell Qend-Qbegin+Qt+dtJ");
        bound(displacement[c]-next,std::abs(next)+gaussActivity[c],1e-10L,"F22 independent physical End native Gauss");}
    check(field.report().relativeLinearResidual<=1e-11&&field.report().maximumRelativeGaussResidual<=1e-10,"F22 retained physical electric linear/Gauss gates");
}
void transaction(const fs::path& root,bool rz,const char* density,double rho)
{
    const auto run=root/(rz?"transaction_rz":"transaction_xy");fs::create_directories(run);
    const auto input=readVofElectrohydrodynamicsInput(fixture(rz,density).save(run,"system").string());const auto& f=input.fluid;
    DictionaryMesh d;d.Nx=f.nx;d.Ny=f.ny;d.minCoord=f.domainMin;d.maxCoord=f.domainMax;d.maxInitialIterations=0;d.maxCellNum=16384;d.maxLevel_local={1};
    auto mesh=std::make_unique<QuadtreeMesh>(d);DictionaryMesh::RefineBox fine;fine.level=1;fine.min=f.domainMin;fine.max=f.domainMax;mesh->applyRefineBoxes({fine});
    VofElectrohydrodynamicsOutput writer(input,{run});VofElectrohydrodynamicsTransaction tx(input,writer);writer.begin();const char* stage="F22 actual Initial projection";
    {
        const auto pub=I::GeometricVofInitialization::prepare(*mesh,rz,[&](Vector2D p){return f.initialJet(p).value;},f.vofInitialization);
        const auto layout=H::HydroDofLayout::build(*mesh,pub.domain());
        const auto metrics=H::HydroPhysicalDivergence::build(*mesh,H::HydroIntegralDivergence::build(*mesh,layout));
        const auto geometry=H::HydroPhaseGeometry::build(*mesh,metrics,pub);
        auto initial=tx.makeInitialVofState(*mesh,geometry,seed(*mesh,layout,rz));native(*mesh,initial,true);
        tx.initializeVof(*mesh,std::move(initial),Wall::now(),stage);
    }
    J::validate(*mesh,tx,{0,0,0,Clock::Position::Initial});unsigned localChanged=0;
    for(int index=0;index<2;++index){const auto previous=tx.accepted()->timeLayer();const auto source=*tx.acceptedCharge();const J::Snapshot before(*mesh,tx,run);
        const Clock target{previous.step+1,previous.time(),f.run.stepEnd(previous.time(),f.run.dt,static_cast<int>(previous.step)),Clock::Position::End};
        const auto start=Wall::now();tx.prepare(*mesh,target,f.step,false,stage);before.unchanged(*mesh,tx,run);
        check(tx.physicalChargeTransport()&&tx.physicalEndField()&&!tx.remapReceipt(),"F22 prepared real physical receipts missing or exposed regrid before commit");
        charge(*mesh,source,*tx.physicalChargeTransport(),*tx.physicalEndField(),rho);
        O::phaseLedger(*mesh,tx.physicalChargeTransport()->transport(),true);O::phaseLedger(*mesh,tx.physicalChargeTransport()->transport(),false);
        tx.commit(*mesh);++accepted;J::validate(*mesh,tx,target);J::checkReceipt(*mesh,tx,previous,target);native(*mesh,*tx.accepted(),false);
        const auto* receipt=tx.remapReceipt();check(receipt&&receipt->decision,"F22 actual regrid decision missing");
        if(receipt->decision->changed){++localChanged;++changed;}
        tx.output(*mesh,start,std::chrono::duration<double>(Wall::now()-start).count(),index==1);tx.releaseStep();
        check(!tx.physicalChargeTransport()&&!tx.physicalEndField(),"F22 released physical receipt retained");
    }
    check(localChanged>0&&tx.accepted()->timeLayer().time()==f.run.endTime,"F22 two-step actual changed-AMR fixture missing");
    writer.complete();check(J::read(run/"output/meta/manifest.json").find("\"status\":\"COMPLETED\"")!=std::string::npos,"F22 terminal writer incomplete");
    tx.clear();check(!tx.accepted()&&!tx.acceptedCharge()&&!tx.acceptedField(),"F22 clear retained accepted subset");
}
void run(const fs::path& requested)
{
    const auto root=fs::absolute(requested)/("f22_pressure_backflow_charge_"+std::to_string(Wall::now().time_since_epoch().count()));fs::create_directories(root);
    reader(root,false,".03125",.03125,"reader_positive_xy");reader(root,true,"-.015625",-.015625,"reader_negative_rz");reader(root,false,nullptr,0,"reader_default_zero");
    reject(root,"no_backflow_optin","Advection.Inlets.drain",[](Fixture& f){VofEhdInputContract::replace(f.text[5]," Backflow { model ambient; velocity (0 0); }","");});
    reject(root,"unknown_charge_patch","Advection.Inlets.ghost",[](Fixture& f){VofEhdInputContract::replace(f.text[6],"drain { chargeDensity .03125; }","ghost { chargeDensity .03125; }");});
    reject(root,"nonfinite_reservoir","Advection.Inlets.drain.chargeDensity",[](Fixture& f){VofEhdInputContract::replace(f.text[6],"drain { chargeDensity .03125; }","drain { chargeDensity nan; }");});
    transaction(root,false,"-.015625",-.015625);transaction(root,true,".03125",.03125);
    check(readers==3&&guards==3&&accepted==4&&changed>=2&&negativeInitial>0&&pressureIncoming>0,"F22 required actual reader/transaction coverage missing");
    std::cout<<"F22_PRESSURE_BACKFLOW_CHARGE_PASS xy_rz accepted="<<accepted<<" changed="<<changed<<" reader_cases="<<readers<<" guards="<<guards
        <<" negative_begin_pressure_rows="<<negativeInitial<<" pressure_incoming_gross_rows="<<pressureIncoming<<" checks="<<checks
        <<" actual_transaction_no_forced_temporal_reversal_or_net0_claim output="<<root.string()<<'\n';
}
}
// END F22 PRESSURE-BACKFLOW INDEPENDENT CHARGE CONTRACT

namespace VofB1TopologyDiagnostic
{
namespace fs=std::filesystem;
namespace H=LAMP::Hydrodynamics;
namespace I=LAMP::Interface;
namespace E=LAMP::Electric;
using Time=LAMP::TimeLayer;
using Wall=HydrodynamicsTransaction::Clock;
void dump(const fs::path& path,const I::GeometricVofPublication& actual,const I::GeometricVofPublication& initial)
{
    require(!fs::exists(path),"topology diagnostic refuses existing inventory output");std::ofstream file(path);file.exceptions(std::ios::failbit|std::ios::badbit);
    file<<std::scientific<<std::setprecision(std::numeric_limits<long double>::max_digits10);
    file<<"publication_identity,step,begin,end,position,time,native_cell,lower_x,lower_y,upper_x,upper_y,liquid,gas,capacity,minor_relative_capacity,source_measured_liquid,source_measured_gas,source_resolved,source_normal_x,source_normal_y,source_offset,initial_liquid,initial_gas,initial_mixed,actual_mixed,new_mixed,delta_liquid,delta_gas\n";
    const auto clock=actual.timeLayer();
    for(int native:actual.domain().cellIndexByRow()){const auto cell=static_cast<std::size_t>(native);const auto& p=actual.plicCells().at(cell);const auto& old=initial.plicCells().at(cell);const auto capacity=actual.cells().at(cell).capacity;
        const bool oldMixed=old.inventory[0]>0&&old.inventory[1]>0,mixed=p.inventory[0]>0&&p.inventory[1]>0;
        file<<actual.identity().value()<<','<<clock.step<<','<<clock.begin<<','<<clock.end<<','<<static_cast<int>(clock.position)<<','<<clock.time()<<','<<native<<','
            <<p.rectangle.lower.x<<','<<p.rectangle.lower.y<<','<<p.rectangle.upper.x<<','<<p.rectangle.upper.y<<','<<p.inventory[0]<<','<<p.inventory[1]<<','<<capacity<<','<<std::min(p.inventory[0],p.inventory[1])/capacity<<','
            <<p.measured[0]<<','<<p.measured[1]<<','<<(p.resolved?1:0)<<','<<p.normal.x<<','<<p.normal.y<<','<<p.offset<<','<<old.inventory[0]<<','<<old.inventory[1]<<','<<(oldMixed?1:0)<<','<<(mixed?1:0)<<','<<(mixed&&!oldMixed?1:0)<<','<<p.inventory[0]-old.inventory[0]<<','<<p.inventory[1]-old.inventory[1]<<'\n';}
    file.close();
}
void run(const fs::path& system,const fs::path& requested)
{
    require(requested.is_absolute()&&(requested.root_name()=="D:"||requested.root_name()=="d:"),"topology diagnostic needs an absolute D: output parent");
    const auto input=readVofElectrohydrodynamicsInput(fs::absolute(system).string());const auto& f=input.fluid;
    require(input.method==VofElectrohydrodynamicsInput::Method::VofSharpSurface&&!f.radial&&!f.vofRegrid&&!f.initialInterfaceBand&&f.vofBoundaryPatches.empty()&&f.vofSymmetryPlanes.empty()&&f.initialVelocity.x==0&&f.initialVelocity.y==0&&!f.vof.southEven90Boundary&&!input.smoothGeometry.allowBoundaryEndpoints,
        "topology diagnostic requires the original same-mesh closed XY B1 physical input");
    const auto root=requested/("b1_topology_diagnostic_"+std::to_string(Wall::now().time_since_epoch().count()));fs::create_directories(root);
    DictionaryMesh specification;specification.Nx=f.nx;specification.Ny=f.ny;specification.minCoord=f.domainMin;specification.maxCoord=f.domainMax;specification.maxInitialIterations=0;specification.maxCellNum=16384;
    QuadtreeMesh mesh(specification);const auto publication=I::GeometricVofInitialization::prepare(mesh,false,[&](Vector2D p){return f.initialJet(p).value;},f.vofInitialization);
    const auto layout=H::HydroDofLayout::build(mesh,publication.domain());const auto metrics=H::HydroPhysicalDivergence::build(mesh,H::HydroIntegralDivergence::build(mesh,layout));
    const auto fluid=H::HydroTwoPhaseState::createVof(mesh,H::HydroPhaseGeometry::build(mesh,metrics,publication),H::HydroNativeValues(layout,std::vector<double>(static_cast<std::size_t>(layout.velocityCount())),std::vector<double>(static_cast<std::size_t>(layout.pressureCount()))));
    VofElectrohydrodynamicsOutput writer(input,{root});VofElectrohydrodynamicsTransaction tx(input,writer);writer.begin();const char* stage="topology diagnostic actual Initial";tx.initializeVof(mesh,fluid,Wall::now(),stage);
    const VofB1ZeroWiringContract::Snapshot preserved(mesh,tx,root);const auto acceptedQ=*tx.acceptedSurfaceCharge();dump(root/"initial.csv",publication,publication);
    {
        std::ofstream faces(root/"native_faces.csv");faces.exceptions(std::ios::failbit|std::ios::badbit);faces<<std::scientific<<std::setprecision(std::numeric_limits<long double>::max_digits10)
            <<"source_identity,native_row,negative_cell,positive_cell,horizontal,center_x,center_y,physical_measure,tangent_spacing\n";
        for(const auto face:publication.domain().faces(mesh))faces<<publication.identity().value()<<','<<face.nativeRow<<','<<face.negativeCell<<','<<face.positiveCell<<','<<(face.horizontal?1:0)<<','<<face.center.x<<','<<face.center.y<<','<<face.area<<','<<face.tangentSpacing<<'\n';faces.close();
    }
    E::ElectricVofSharpProblem problem;problem.epsilon={{input.materials.epsilonLiquid,input.materials.epsilonGas}};problem.conductivity={{input.materials.conductivityLiquid,input.materials.conductivityGas}};
    for(std::size_t side=0;side<4;++side){const auto boundary=input.boundaries[side];problem.boundaries[side].type=boundary.type==E::ElectricVofBoundary::Type::Dirichlet?E::ElectricVofSharpBoundary::Type::Dirichlet:E::ElectricVofSharpBoundary::Type::HomogeneousNeumann;problem.boundaries[side].value=[value=boundary.value](Vector2D){return value;};}
    unsigned callbacks=0,dumped=0;auto config=f.vof;config.requireMaxwellCorrection=true;
    config.maxwellLoad=[&](const QuadtreeMesh& grid,const I::GeometricVofTransport& transport,const H::HydroPhaseGeometry& mid,const H::HydroPressureLoad& pressure){
        ++callbacks;const auto prefix="stage_"+std::to_string(callbacks);require(transport.source().identity()==publication.identity()&&mid.vofPublication().identity()==transport.midpoint().identity(),"diagnostic callback lost actual source or Midpoint publication");
        dump(root/(prefix+"_source.csv"),transport.source(),publication);dump(root/(prefix+"_mid.csv"),transport.midpoint(),publication);dump(root/(prefix+"_end.csv"),transport.endpoint(),publication);
        {std::ofstream flux(root/(prefix+"_actual_flux.csv"));flux.exceptions(std::ios::failbit|std::ios::badbit);flux<<std::scientific<<std::setprecision(std::numeric_limits<long double>::max_digits10)
                <<"source_identity,midpoint_identity,endpoint_identity,step,begin,end,position,native_row,negative_cell,positive_cell,begin_normal,end_normal,half_total,half_raw_liquid,half_liquid,half_gas,half_limiter,half_mass,half_positive_total,half_negative_total,half_positive_liquid,half_positive_gas,half_negative_liquid,half_negative_gas,full_total,full_raw_liquid,full_liquid,full_gas,full_limiter,full_mass,full_positive_total,full_negative_total,full_positive_liquid,full_positive_gas,full_negative_liquid,full_negative_gas\n";
            const auto target=transport.endTime();const auto& slab=transport.advector();for(const auto& face:transport.faces())flux<<transport.source().identity().value()<<','<<transport.midpoint().identity().value()<<','<<transport.endpoint().identity().value()<<','<<target.step<<','<<target.begin<<','<<target.end<<','<<static_cast<int>(target.position)<<','<<face.nativeRow<<','<<face.negativeCell<<','<<face.positiveCell<<','<<slab.beginNormal().at(static_cast<std::size_t>(face.nativeRow))<<','<<slab.endNormal().at(static_cast<std::size_t>(face.nativeRow))<<','<<face.half.total<<','<<face.half.rawLiquid<<','<<face.half.liquid<<','<<face.half.gas<<','<<face.half.limiter<<','<<face.half.mass<<','<<face.half.positiveTotal<<','<<face.half.negativeTotal<<','<<face.half.positiveLiquid<<','<<face.half.positiveGas<<','<<face.half.negativeLiquid<<','<<face.half.negativeGas<<','<<face.full.total<<','<<face.full.rawLiquid<<','<<face.full.liquid<<','<<face.full.gas<<','<<face.full.limiter<<','<<face.full.mass<<','<<face.full.positiveTotal<<','<<face.full.negativeTotal<<','<<face.full.positiveLiquid<<','<<face.full.positiveGas<<','<<face.full.negativeLiquid<<','<<face.full.negativeGas<<'\n';flux.close();}
        ++dumped;
        // Diagnostic-only copy of the App's public module chain. The original
        // HydroVofStep owns all transport, native NS stages and physical gates.
        const auto smooth=I::VofSmoothGeometry::prepare(grid,mid.vofPublication(),input.smoothGeometry);
        const auto moved=[&]{try{return E::ElectricSurfaceChargeTransport::prepareStage(grid,acceptedQ,transport,smooth,input.surfaceTransport);}
            catch(const E::ElectricSurfaceChargeTransferError& error){if(error.reason()==E::ElectricSurfaceChargeTransferFailure::ExcessiveDisplacement||error.reason()==E::ElectricSurfaceChargeTransferFailure::ProjectionDistance)throw H::HydroTwoPhaseStepRetry(H::HydroTwoPhaseRetryReason::Transport,error.what());throw;}}();
        const auto field=E::ElectricVofSharpField::prepareOhmic(grid,moved,problem,{input.controls.relativeLinearTolerance,input.controls.relativeConservationTolerance});const auto stress=E::ElectricVofSharpMaxwellStress::prepare(grid,field);return H::HydroMaxwellLoad::prepareSharp(grid,pressure,mid,stress);
    };
    bool prepared=false;std::string reason;int geometryCell=-1,geometryReason=-1;
    const Time first{1,0,f.run.stepEnd(0,f.run.dt,0),Time::Position::End};
    try{const auto candidate=H::HydroVofStep::prepare(mesh,*tx.accepted(),first,config);prepared=true;require(candidate.beginState().sameIdentity(*tx.accepted()),"diagnostic prepare changed its source owner");}
    catch(const I::VofSmoothGeometryError& error){reason=error.what();geometryCell=error.nativeCell();geometryReason=static_cast<int>(error.reason());}
    catch(const std::exception& error){reason=error.what();}
    require(callbacks>0,"diagnostic never reached an actual transport callback");require(dumped==callbacks,"an actual callback inventory/flux export failed");preserved.unchanged(mesh,tx,root);
    std::ofstream result(root/"diagnostic_status.txt");result.exceptions(std::ios::failbit|std::ios::badbit);result<<std::setprecision(std::numeric_limits<long double>::max_digits10)
        <<"diagnostic_only=1\naccepted_commit=0\naccepted_snapshot_and_output_prefix_unchanged=1\nactual_callbacks="<<callbacks<<"\noriginal_prepare_completed="<<(prepared?1:0)<<"\ngeometry_failure_cell="<<geometryCell<<"\ngeometry_failure_reason="<<geometryReason<<"\naccepted_source_identity="<<publication.identity().value()<<"\naccepted_step="<<tx.accepted()->timeLayer().step<<"\naccepted_begin="<<tx.accepted()->timeLayer().begin<<"\naccepted_end="<<tx.accepted()->timeLayer().end<<"\naccepted_position="<<static_cast<int>(tx.accepted()->timeLayer().position)<<"\nobserved_exception="<<reason<<'\n';result.close();tx.clear();
    std::cout<<"B1 topology diagnostic CAPTURED callbacks="<<callbacks<<" prepare_completed="<<(prepared?1:0)<<" no_commit=1 output="<<root.string()<<'\n';
}
}


// BEGIN F26 ACTUAL A SURFACE-PROJECTION INITIAL CONTRACT
namespace F26SurfaceProjectionInitialContract {
namespace fs=std::filesystem;
namespace H=LAMP::Hydrodynamics;
namespace I=LAMP::Interface;
namespace E=LAMP::Electric;
namespace FV=LAMP::Geometry::Detail::FiniteVolume;
namespace Old=VofEhdInputContract;
using Fixture=Old::Fixture;
using Type=VofElectrohydrodynamicsInput::InitialChargeType;
using Time=LAMP::TimeLayer;
using Wall=HydrodynamicsTransaction::Clock;
unsigned readers=0,initials=0,projected=0,guards=0,cellChecks=0,arcsChecked=0,ordinarySteps=0;
void check(bool ok,const char* why){if(!ok)throw std::runtime_error(why);}
void exact(double a,double b,const char* why){check(std::isfinite(a)&&std::isfinite(b)&&a==b&&(a!=0||std::signbit(a)==std::signbit(b)),why);}
template<class F>void logicReject(F action,const char* why){bool caught=false;try{action();}catch(const std::logic_error&){caught=true;}check(caught,why);++guards;}
void validateAccepted(const QuadtreeMesh& mesh,const VofElectrohydrodynamicsTransaction& tx,Time expected,const E::ElectricVofControls& controls)
{
    const auto* fluid=tx.accepted();const auto* charge=tx.acceptedCharge();const auto* field=tx.acceptedField();
    check(fluid&&charge&&field&&!tx.acceptedSurfaceCharge()&&!tx.acceptedSharpField(),"F26 accepted ordinary A owner set is incomplete");
    fluid->validateFor(mesh);const auto& pub=fluid->phaseGeometry().vofPublication();
    pub.validateFor(mesh,pub.identity(),expected);charge->validateFor(mesh,expected);field->validateFor(mesh,pub,field->sourceCharge(),expected);
    check(VofEhdJointContract::sameClock(fluid->timeLayer(),expected)&&VofEhdJointContract::sameClock(pub.timeLayer(),expected)
        &&VofEhdJointContract::sameClock(charge->timeLayer(),expected)&&VofEhdJointContract::sameClock(field->geometryTime(),expected)
        &&VofEhdJointContract::sameClock(field->candidateTime(),expected)&&field->geometry().identity()==pub.identity()
        &&field->chargeCandidate().sameInstance(*charge),"F26 accepted ordinary A owner/fullclock mismatch");
    if(expected.position==Time::Position::Initial)
        check(field->isElectrostaticSnapshot()&&field->effectiveDt()==0&&field->sourceCharge().sameInstance(*charge)
            &&!field->advection()&&VofEhdJointContract::sameClock(field->chargeSourceTime(),expected),
            "F26 actual Initial snapshot changed its same-Q owner or fabricated an interval");
    else check(expected.position==Time::Position::End&&!field->isElectrostaticSnapshot()
        &&field->effectiveDt()==expected.end-expected.begin&&field->effectiveDt()>0&&field->advection(),
        "F26 ordinary A End lacks its true Ohmic/advection interval");
    const auto& report=field->report();
    check(report.relativeLinearResidual<=controls.relativeLinearTolerance&&report.maximumRelativeGaussResidual<=controls.relativeConservationTolerance,
        "F26 accepted ordinary A field changed retained linear/Gauss gates");
    for(double value:fluid->values().normalVelocity())check(std::isfinite(value),"F26 accepted ordinary A has nonfinite velocity");
    for(double value:fluid->values().pressure())check(std::isfinite(value),"F26 accepted ordinary A has nonfinite pressure");
    for(double value:charge->integratedCharge())check(std::isfinite(value),"F26 accepted ordinary A has nonfinite signed charge");
    for(double value:field->potential())check(std::isfinite(value),"F26 accepted ordinary A has nonfinite potential");
}
Fixture fixture(bool projection,bool explicitCell=false,bool spatial=false)
{
    Fixture f;
    f.text[0]="Run { app VofElectrohydrodynamics; } Parallel { threads 1; } Time { units code; dt .000001; endTime .000001; } Output { Fields { enabled false; } }";
    f.text[1]="Geometry { coordinates xy; Domain { min (0 0); max (1 1); } Phases { negative 7; positive 3; } Initial { phi \"(x-centerX)^2+(y-centerY)^2-radius^2\"; } }";
    f.text[2]="Mesh { resolution (32 32); initialMesh uniform; }";
    f.setConstants("centerX .487; centerY .513; radius .18; chargeScale .5;");
    const std::string initial=projection?std::string("type surfaceChargeProjection; surfaceChargeDensity \"")+(spatial?"chargeScale*(x-centerX)":"chargeScale")+"\";":
        std::string(explicitCell?"type cellDensity; ":"")+"chargeDensity \"chargeScale*(x-.5)\";";
    f.text[6]="Electric { method vofFiniteVolume; Initial { "+initial+" } Boundary { W { type Dirichlet; value -.25; } S { type Dirichlet; value -.5; } E { type Dirichlet; value 1; } N { type Neumann; value 0; } } }";
    return f;
}
VofElectrohydrodynamicsInput readChecked(const fs::path& root,const Fixture& f,const std::string& name,bool projection,bool explicitCell)
{
    const auto input=readVofElectrohydrodynamicsInput(f.save(root,name).string());++readers;
    check(input.method==VofElectrohydrodynamicsInput::Method::VofFiniteVolume&&input.initialChargeType==(projection?Type::SurfaceChargeProjection:Type::CellDensity),
        "F26 reader selected another electric route or Initial type");
    const auto& fluid=input.fluid;check(fluid.documents.size()==7,"F26 changed seven-document ownership");
    const char* field=projection?"Electric.Initial.surfaceChargeDensity":"Electric.Initial.chargeDensity";
    bool typeTrace=false,densityTrace=false;
    for(const auto& e:fluid.trace.entries){
        if(e.field=="Electric.Initial.type")typeTrace=e.units=="dimensionless"&&e.consumer=="ElectricVofCharge.initialization"
            &&e.used&&e.defaulted==(!projection&&!explicitCell)&&std::holds_alternative<std::string>(e.value)
            &&std::get<std::string>(e.value)==(projection?"surfaceChargeProjection":"cellDensity")&&fs::path(e.file).filename()=="electric.dict";
        if(e.field==field)densityTrace=e.units==(projection?"code_charge/code_interface_measure":"code_charge/code_volume")
            &&e.consumer==(projection?"ElectricVofCharge.initialSurfaceMomentProjection":"ElectricVofCharge.initialIntegratedCellCharge")
            &&e.used&&!e.defaulted&&fs::path(e.file).filename()=="electric.dict";
    }
    check(typeTrace&&densityTrace,"F26 Initial type/expression units, consumer, explicit/default provenance missing");
    if(projection)logicReject([&]{input.initialChargeDensity({.25,.25});},"F26 projection surface density leaked as volume density");
    else logicReject([&]{input.initialSurfaceChargeDensity({.25,.25});},"F26 old volume density leaked as surface density");
    const auto effective=root/(name+"_effective");fs::create_directories(effective);
    for(std::size_t k=0;k<Old::names.size();++k){check(fluid.documents[k].sourceText==f.text[k],"F26 reader modified source dictionary text");
        Old::write(effective/Old::names[k],fluid.documents[k].effectiveText);}
    const auto& electric=fluid.documents[6].effectiveText;
    if(projection)check(electric.find("Initial { type surfaceChargeProjection; surfaceChargeDensity ")!=std::string::npos,"F26 explicit surface projection did not survive effective text");
    else if(explicitCell)check(electric.find("Initial { type cellDensity; chargeDensity ")!=std::string::npos,"F26 explicit cellDensity type omitted from effective text");
    else check(electric.find("Initial { chargeDensity ")!=std::string::npos&&electric.find("Initial { type ")==std::string::npos,
        "F26 omitted type changed the old effective Initial block format");
    const auto again=readVofElectrohydrodynamicsInput(effective.string());
    check(again.method==input.method&&again.initialChargeType==input.initialChargeType&&again.fluid.documents.size()==7,
        "F26 effective roundtrip changed type/method/document ownership");
    for(Vector2D p:std::array<Vector2D,3>{{{.25,.25},{.487,.513},{.75,.75}}})
        exact(projection?input.initialSurfaceChargeDensity(p):input.initialChargeDensity(p),
            projection?again.initialSurfaceChargeDensity(p):again.initialChargeDensity(p),"F26 effective roundtrip changed actual signed spatial expression");
    return input;
}
std::string payload(const I::GeometricVofPublication& pub)
{
    std::string out;
    for(const auto& c:pub.plicCells()){for(auto x:c.inventory)VofEhdJointContract::append(out,x);for(auto x:c.measured)VofEhdJointContract::append(out,x);
        VofEhdJointContract::append(out,c.normal.x);VofEhdJointContract::append(out,c.normal.y);
        VofEhdJointContract::append(out,c.offset);VofEhdJointContract::append(out,c.resolved);}
    return out;
}
std::unique_ptr<QuadtreeMesh> meshFor(const VofElectrohydrodynamicsInput& input,bool refined=false)
{
    const auto& f=input.fluid;DictionaryMesh d;d.Nx=f.nx;d.Ny=f.ny;d.minCoord=f.domainMin;d.maxCoord=f.domainMax;d.maxInitialIterations=0;d.maxCellNum=8192;
    auto mesh=std::make_unique<QuadtreeMesh>(d);
    if(refined){DictionaryMesh::RefineBox b;b.level=1;b.min={.5,0};b.max={1,1};mesh->applyRefineBoxes({b});}
    return mesh;
}
std::vector<double> actualInitial(const fs::path& root,const VofElectrohydrodynamicsInput& input,const std::string& name,bool projection,bool spatial,bool advance)
{
    auto mesh=meshFor(input);const auto& f=input.fluid;const auto run=root/(name+"_transaction");fs::create_directories(run);
    VofElectrohydrodynamicsOutput writer(input,{run});VofElectrohydrodynamicsTransaction tx(input,writer);writer.begin();const char* stage="F26 actual Initial";
    const auto pub=I::GeometricVofInitialization::prepare(*mesh,false,[&](Vector2D p){return f.initialJet(p).value;},f.vofInitialization);
    const auto original=payload(pub);const auto layout=H::HydroDofLayout::build(*mesh,pub.domain());
    const auto metrics=H::HydroPhysicalDivergence::build(*mesh,H::HydroIntegralDivergence::build(*mesh,layout));
    const auto phase=H::HydroPhaseGeometry::build(*mesh,metrics,pub);
    const auto fluid=H::HydroTwoPhaseState::createVof(*mesh,phase,H::HydroNativeValues(layout,
        std::vector<double>(layout.velocityCount()),std::vector<double>(layout.pressureCount())));
    tx.initializeVof(*mesh,fluid,Wall::now(),stage);++initials;
    if(spatial)VofEhdJointContract::validate(*mesh,tx,{0,0,0,Time::Position::Initial});
    else validateAccepted(*mesh,tx,{0,0,0,Time::Position::Initial},input.controls);
    check(tx.acceptedCharge()&&tx.acceptedField()&&!tx.acceptedSurfaceCharge()&&!tx.acceptedSharpField()
        &&tx.accepted()->phaseGeometry().vofIdentity()==pub.identity(),"F26 Initial retained B1 state or another actual VOF owner");
    std::vector<double> expected(mesh->cellIndexNum());FV::Dyadic momentTotal,cellTotal;unsigned nonzero=0,positive=0,negative=0,empty=0;
    if(projection){
        const auto smooth=I::VofSmoothGeometry::prepare(*mesh,pub,input.smoothGeometry);
        smooth.validateFor(*mesh,pub);check(smooth.sourceIdentity()==pub.identity()&&!smooth.arcs().empty()&&smooth.report().boundaryEndpoints==0,
            "F26 oracle needs actual same-source closed smooth geometry");
        const auto surface=E::ElectricVofSurfaceCharge::fromDensity(*mesh,smooth,[&](std::size_t arc,double t){const auto p=smooth.arcs()[arc].sample(t);
            return input.initialSurfaceChargeDensity({static_cast<double>(p.x),static_cast<double>(p.y)});});
        surface.validateFor(*mesh,smooth,pub.timeLayer());std::vector<FV::Dyadic> byCell(mesh->cellIndexNum());std::vector<unsigned> arcs(mesh->cellIndexNum());
        for(std::size_t arc=0;arc<smooth.arcs().size();++arc){const int ci=smooth.arcs()[arc].nativeCell;
            check(ci>=0&&ci<mesh->cellIndexNum(),"F26 typed oracle arc has no actual native-cell mapping");
            exact(surface.moments()[arc][0],surface.integratedCharge()[arc],"F26 typed M0 is not its integrated surface charge");
            byCell[ci].add(surface.moments()[arc][0]);momentTotal.add(surface.moments()[arc][0]);++arcs[ci];++arcsChecked;}
        for(int ci=0;ci<mesh->cellIndexNum();++ci){expected[ci]=FV::roundDouble(byCell[ci]).value;if(!arcs[ci]){check(expected[ci]==0,"F26 empty arc-cell oracle is nonzero");++empty;}}
        check(empty>0,"F26 projection oracle must distinguish interface and non-interface native cells");++projected;
    }else for(int ci=0;ci<mesh->cellIndexNum();++ci){const auto center=mesh->cells(mesh->cellIndexToLabel[ci]).centerCoord();
        expected[ci]=FV::roundDouble(FV::Dyadic::product({input.initialChargeDensity(center),pub.domain().cellVolumes()[pub.domain().cellRowByIndex()[ci]]})).value;}
    const auto actual=tx.acceptedCharge()->integratedCharge();check(actual.size()==expected.size(),"F26 actual Initial Q vector has another native layout");
    for(std::size_t ci=0;ci<actual.size();++ci){exact(actual[ci],expected[ci],"F26 actual Initial per-cell Q differs from original typed-source exact oracle");
        nonzero+=actual[ci]!=0;positive+=actual[ci]>0;negative+=actual[ci]<0;cellTotal.add(actual[ci]);++cellChecks;}
    check(nonzero>0&&(spatial?positive>0&&negative>0:positive>0&&negative==0),"F26 actual signed Initial charge witness missing");
    check(payload(pub)==original,"F26 transient projection mutated original VOF inventory/plane bytes");
    if(advance){
        const auto oldFluid=*tx.accepted();const auto oldCharge=*tx.acceptedCharge();const auto savedQ=oldCharge.integratedCharge();
        const auto savedU=oldFluid.values().normalVelocity(),savedP=oldFluid.values().pressure();const auto saved=payload(pub);
        const Time end{1,0,f.run.dt,Time::Position::End};tx.prepare(*mesh,end,f.step,false,stage);
        check(tx.accepted()->sameIdentity(oldFluid)&&tx.acceptedCharge()->sameInstance(oldCharge)
            &&tx.acceptedCharge()->integratedCharge()==savedQ&&tx.accepted()->values().normalVelocity()==savedU
            &&tx.accepted()->values().pressure()==savedP&&payload(pub)==saved,"F26 ordinary A prepare changed accepted Initial source");
        check(tx.preparedStep()&&tx.physicalChargeTransport()&&tx.physicalEndField()
            &&!tx.preparedStep()->transport().smoothMidpoint()&&!tx.preparedStep()->transport().smoothEndpoint()
            &&!tx.acceptedSurfaceCharge()&&!tx.acceptedSharpField(),"F26 projection incorrectly retained sharp transport after Initial");
        tx.commit(*mesh);validateAccepted(*mesh,tx,end,input.controls);
        check(tx.acceptedCharge()&&tx.acceptedField()&&!tx.acceptedSurfaceCharge()&&!tx.acceptedSharpField(),
            "F26 ordinary End switched route after Initial projection");
        tx.output(*mesh,Wall::now(),0,true);tx.releaseStep();++ordinarySteps;
    }
    if(advance)writer.complete();
    else writer.failure(0,0,0,"F26 finite Initial-only qualification",
        "Expected bounded qualification prefix: no physical time window was requested or executed by this contract; retained Initial only.");
    tx.clear();
    std::cout<<std::setprecision(17)<<"F26_ACTUAL_INITIAL,name,"<<name<<",projection,"<<projection<<",cells,"<<actual.size()<<",nonzero,"<<nonzero
        <<",positive,"<<positive<<",negative,"<<negative<<",oracle_total_surface_M0,"<<momentTotal.rounded()<<",actual_total_cell_Q,"<<cellTotal.rounded()
        <<",normalization,none,subsequent_A_steps,"<<(advance?1:0)<<'\n';
    return actual;
}
template<class F>void readerReject(const fs::path& root,Fixture f,const std::string& name,const char* token,F change)
{
    change(f);bool caught=false;try{(void)readVofElectrohydrodynamicsInput(f.save(root,name).string());}
    catch(const std::exception& e){caught=std::string(e.what()).find(token)!=std::string::npos;}
    check(caught,"F26 selected-field/type/scope rejection lost its actual input provenance");++guards;
}
void negativeInputs(const fs::path& root)
{
    const auto p=fixture(true);
    readerReject(root,p,"unknown_type","Electric.Initial.type",[](Fixture& f){Old::replace(f.text[6],"surfaceChargeProjection","surfaceAlias");});
    readerReject(root,p,"projection_bulk_only","Electric.Initial.chargeDensity",[](Fixture& f){Old::replace(f.text[6],"surfaceChargeDensity","chargeDensity");});
    readerReject(root,p,"both_charge_keys","Electric.Initial.chargeDensity",[](Fixture& f){Old::replace(f.text[6],"Initial {","Initial { chargeDensity 0;");});
    readerReject(root,p,"cell_with_surface","Electric.Initial.surfaceChargeDensity",[](Fixture& f){Old::replace(f.text[6],"surfaceChargeProjection","cellDensity");});
    readerReject(root,p,"missing_surface","Electric.Initial.surfaceChargeDensity",[](Fixture& f){Old::replace(f.text[6],"surfaceChargeDensity \"chargeScale\";","");});
    readerReject(root,p,"timed_surface","Electric.Initial.surfaceChargeDensity",[](Fixture& f){Old::replace(f.text[6],"\"chargeScale\"","\"chargeScale+t\"");});
    readerReject(root,p,"nonfinite_surface","Electric.Initial.surfaceChargeDensity",[](Fixture& f){Old::replace(f.text[6],"\"chargeScale\"","\"chargeScale/(x-x)\"");});
    readerReject(root,p,"rz_projection","Electric.Initial.type",[](Fixture& f){Old::replace(f.text[1],"coordinates xy","coordinates rz");
        Old::replace(f.text[1],"(x-centerX)^2+(y-centerY)^2","(r-centerX)^2+(z-centerY)^2");});
    readerReject(root,p,"open_projection","Electric.Initial.type",[](Fixture& f){Old::replace(f.text[5],"Boundary { type impermeable; }",
        "Boundary { type prescribedFlow; Patches { feed { side S; interval (0 .5); type prescribedVelocity; phase liquid; velocity (0 .125); } out { side N; interval (0 1); type pressureOutlet; phase gas; pressure 0; } } }");});
    readerReject(root,p,"contact_projection","Electric.Initial.type",[](Fixture& f){Old::replace(f.text[3],"Boundary { type none; }","Boundary { type symmetry; sides (bottom); }");});
    auto b1=p;Old::replace(b1.text[6],"vofFiniteVolume","vofSharpSurface");Old::replace(b1.text[6],"type surfaceChargeProjection; ","");
    const auto oldB1=readVofElectrohydrodynamicsInput(b1.save(root,"unchanged_b1_reader").string());++readers;
    check(oldB1.method==VofElectrohydrodynamicsInput::Method::VofSharpSurface&&oldB1.initialSurfaceChargeDensity({.25,.25})==.5,
        "F26 changed original B1 initial surface accessor");
    readerReject(root,b1,"b1_type_is_not_a_projection_port","Electric.Initial.type",[](Fixture& f){Old::replace(f.text[6],"Initial {","Initial { type surfaceChargeProjection;");});
    auto future=p;future.text[2]="Mesh { resolution (32 32); initialMesh uniform; Regrid { minLevel 0; maxLevel 1; interfaceBandLayers 2; vorticityTolerance .5; coarsenRatio .25; } }";
    const auto amr=readVofElectrohydrodynamicsInput(future.save(root,"future_amr_allowed").string());++readers;
    check(amr.initialChargeType==Type::SurfaceChargeProjection&&amr.fluid.vofRegrid.has_value(),"F26 projection incorrectly forbids future A AMR policy");
}
void actualScopeReject(const fs::path& root,bool nonuniform)
{
    const auto input=readVofElectrohydrodynamicsInput(fixture(true).save(root,nonuniform?"nonuniform_input":"empty_input").string());
    auto mesh=meshFor(input,nonuniform);const auto run=root/(nonuniform?"nonuniform_actual_scope":"empty_actual_scope");fs::create_directories(run);
    VofElectrohydrodynamicsOutput writer(input,{run});VofElectrohydrodynamicsTransaction tx(input,writer);writer.begin();const char* stage="F26 expected actual Initial scope";
    const auto pub=I::GeometricVofInitialization::prepare(*mesh,false,[&](Vector2D p){return nonuniform?input.fluid.initialJet(p).value:1.;},input.fluid.vofInitialization);
    check(nonuniform?!pub.domain().uniformLeaves():pub.report().mixedCells==0,"F26 runtime negative did not produce its declared actual source geometry");
    const auto layout=H::HydroDofLayout::build(*mesh,pub.domain());const auto metrics=H::HydroPhysicalDivergence::build(*mesh,H::HydroIntegralDivergence::build(*mesh,layout));
    const auto phase=H::HydroPhaseGeometry::build(*mesh,metrics,pub);const auto before=payload(pub);bool caught=false;std::string reason;
    try{tx.initializeVof(*mesh,H::HydroTwoPhaseState::createVof(*mesh,phase,H::HydroNativeValues(layout,
        std::vector<double>(layout.velocityCount()),std::vector<double>(layout.pressureCount()))),Wall::now(),stage);}
    catch(const std::exception& e){reason=e.what();caught=nonuniform?reason.find("surfaceChargeProjection requires actual closed uniform XY Initial")!=std::string::npos:
        std::string(stage)=="A actual Initial surface-charge projection geometry";}
    check(caught&&!tx.accepted()&&!tx.acceptedCharge()&&!tx.acceptedField()&&!tx.acceptedSurfaceCharge()&&!tx.acceptedSharpField()
        &&payload(pub)==before,"F26 actual geometry scope failure published/mutated a partial Initial owner");
    writer.failure(0,0,0,stage,"expected F26 qualification rejection: "+reason);++guards;
}
void run(const fs::path& root)
{
    check(root.is_absolute()&&root.root_name().string()=="D:"&&!fs::exists(root),"F26 selector needs a unique explicit absolute D: work directory");
    fs::create_directories(root);
    const auto cell=readChecked(root,fixture(false),"old_default_density",false,false);
    const auto explicitCell=readChecked(root,fixture(false,true),"explicit_cell_density",false,true);
    const auto oldQ=actualInitial(root,cell,"old_default_density",false,true,false);
    const auto explicitQ=actualInitial(root,explicitCell,"explicit_cell_density",false,true,false);
    check(oldQ.size()==explicitQ.size(),"F26 explicit/default volume Q layouts differ");
    for(std::size_t k=0;k<oldQ.size();++k)exact(oldQ[k],explicitQ[k],"F26 explicit cellDensity changed legacy default initial Q");
    const auto constant=readChecked(root,fixture(true),"constant_surface_projection",true,false);
    actualInitial(root,constant,"constant_surface_projection",true,false,true);
    const auto signedInput=readChecked(root,fixture(true,false,true),"signed_surface_projection",true,false);
    actualInitial(root,signedInput,"signed_surface_projection",true,true,false);
    negativeInputs(root);actualScopeReject(root,true);actualScopeReject(root,false);
    check(readers==6&&initials==4&&projected==2&&ordinarySteps==1&&guards==17&&cellChecks==4096&&arcsChecked>0,
        "F26 required finite Initial/reader/signed/legacy/ordinary-A coverage incomplete");
    std::cout<<"F26_A_SURFACE_PROJECTION_INITIAL_PASS,readers,"<<readers<<",actual_initials,"<<initials<<",projections,"<<projected
        <<",per_cell_Q_checks,"<<cellChecks<<",typed_arcs,"<<arcsChecked<<",guards,"<<guards<<",ordinary_A_steps,"<<ordinarySteps
        <<",exact_value_and_signed_zero,no_total_Q_normalization,no_B1_transport_retained,NO_LONGRUN_OR_A_B1_ACCURACY_CLAIM\n";
}
} // namespace F26SurfaceProjectionInitialContract
// END F26 ACTUAL A SURFACE-PROJECTION INITIAL CONTRACT

namespace VofChargeAdvectionInputContract {
namespace fs=std::filesystem;
using Fixture=VofEhdInputContract::Fixture;
using Scheme=LAMP::Electric::ElectricVofChargeAdvectionScheme;
using VofEhdInputContract::replace;
using VofEhdInputContract::write;
unsigned readers=0,guards=0;
void scheme(Fixture& f,const std::string& value){replace(f.text[6],"Electric {","Electric { Advection { scheme "+value+"; }");}
VofElectrohydrodynamicsInput read(const fs::path& root,const Fixture& f,const std::string& name,Scheme expected,bool isDefault){
    const auto input=readVofElectrohydrodynamicsInput(f.save(root,name).string());++readers;
    require(input.controls.chargeAdvection==expected,"charge advection reader produced wrong typed scheme");
    bool trace=false;for(const auto& e:input.fluid.trace.entries)if(e.field=="Electric.Advection.scheme")
        trace=e.used&&e.defaulted==isDefault&&e.units=="dimensionless"&&fs::path(e.file).filename()=="electric.dict";
    require(trace,"charge advection scheme trace lacks default/explicit electric.dict provenance");
    const auto effective=root/(name+"_effective");fs::create_directories(effective);
    for(std::size_t k=0;k<VofEhdInputContract::names.size();++k)write(effective/VofEhdInputContract::names[k],input.fluid.documents[k].effectiveText);
    const auto again=readVofElectrohydrodynamicsInput(effective.string());++readers;
    require(again.controls.chargeAdvection==expected&&again.inletChargeDensityByPatch==input.inletChargeDensityByPatch,"charge advection effective text lost scheme or physical reservoir state");
    return input;
}
template<class Action>void rejectScheme(const fs::path& root,Fixture f,const std::string& name,const char* field,Action action){
    action(f);bool caught=false;try{readVofElectrohydrodynamicsInput(f.save(root,name).string());}
    catch(const std::exception& e){caught=std::string(e.what()).find(field)!=std::string::npos;}
    require(caught,"charge advection invalid combination accepted or real field diagnostic lost");++guards;
}
void run(const fs::path& requested){
    require(requested.is_absolute()&&requested.root_name().string()=="D:","charge advection selector requires an absolute D: output directory");
    const auto root=requested/("reader_"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    require(!fs::exists(root),"charge advection selector output must be fresh");fs::create_directories(root);
    require(LAMP::Electric::ElectricVofControls{}.chargeAdvection==Scheme::DonorCell,"charge advection changed direct module P0 default");
    const auto defaultXY=read(root,Fixture{},"closed_xy_default",Scheme::LimitedLinearFct,true);
    read(root,Fixture{true},"closed_rz_default",Scheme::LimitedLinearFct,true);
    require(defaultXY.fluid.vofBoundaryPatches.empty()&&defaultXY.inletChargeDensityByPatch.empty(),"closed FCT default manufactured a reservoir");
    auto explicitLow=Fixture{};scheme(explicitLow,"donorCell");read(root,explicitLow,"closed_scheme_only_p0",Scheme::DonorCell,false);
    auto explicitHigh=Fixture{true};scheme(explicitHigh,"limitedLinearFct");read(root,explicitHigh,"closed_scheme_only_fct",Scheme::LimitedLinearFct,false);
    auto open=VofEhdOpenInputContract::open();VofEhdOpenInputContract::density(open,"feed","-.25");
    replace(open.text[6],"Advection {","Advection { scheme donorCell;");
    const auto opened=read(root,open,"open_p0_signed_inlet",Scheme::DonorCell,false);
    require(opened.inletChargeDensityByPatch.size()==2&&opened.inletChargeDensityByPatch.front()==-.25,"scheme changed signed physical reservoir");
    rejectScheme(root,Fixture{},"invalid_scheme","Electric.Advection.scheme",[](Fixture& f){scheme(f,"unlimitedLinear");});
    rejectScheme(root,Fixture{},"closed_inlets_forbidden","Electric.Advection",[](Fixture& f){VofEhdOpenInputContract::density(f,"feed","1");});
    rejectScheme(root,VofEhdOpenInputContract::open(),"unknown_inlet","Advection.Inlets.ghost",[](Fixture& f){VofEhdOpenInputContract::density(f,"ghost","1");});
    auto sharp=Fixture{};sharp.text[2]="Mesh { resolution (8 8); initialMesh uniform; }";
    replace(sharp.text[6],"vofFiniteVolume","vofSharpSurface");replace(sharp.text[6],"chargeDensity","surfaceChargeDensity");
    const auto b1=readVofElectrohydrodynamicsInput(sharp.save(root,"b1_unchanged").string());++readers;
    require(b1.method==VofElectrohydrodynamicsInput::Method::VofSharpSurface&&b1.controls.chargeAdvection==Scheme::DonorCell,"A FCT default changed B1 controls");
    for(const auto& e:b1.fluid.trace.entries)require(e.field!="Electric.Advection.scheme","B1 acquired unused A advection default");
    rejectScheme(root,sharp,"b1_advection_forbidden","Electric.Advection",[](Fixture& f){scheme(f,"limitedLinearFct");});
    require(readers==11&&guards==4,"charge advection reader coverage incomplete");
    std::cout<<"VOF_CHARGE_ADVECTION_INPUT_PASS,readers,"<<readers<<",guards,"<<guards<<",A_default_FCT,module_default_P0,B1_unchanged,closed_scheme_only,explicit_roundtrip_trace,signed_inlets\n";
}
} // namespace VofChargeAdvectionInputContract
namespace VofChargeRegridContract {
namespace fs=std::filesystem;
namespace H=LAMP::Hydrodynamics;
namespace I=LAMP::Interface;
using Fixture=VofEhdInputContract::Fixture;
using VofEhdInputContract::replace;
using VofEhdInputContract::write;
unsigned readers=0,guards=0,steps=0;
void option(Fixture& f,const std::string& value){
    replace(f.text[6],"Electric {","Electric { Regrid { chargeDensityJumpTolerance "+value+"; }");
}
void schema(const fs::path& root){
    for(bool rz:{false,true}){
        Fixture f(rz);option(f,".0001");
        const auto input=readVofElectrohydrodynamicsInput(f.save(root,rz?"reader_rz":"reader_xy").string());++readers;
        require(input.chargeDensityJumpTolerance==.0001,"charge regrid typed threshold lost");
        bool traced=false;for(const auto& e:input.fluid.trace.entries)if(e.field=="Electric.Regrid.chargeDensityJumpTolerance")
            traced=e.used&&!e.defaulted&&e.units=="code_charge/code_volume"&&e.consumer=="ElectricVofRefinement.signedDensityJump";
        require(traced,"charge regrid threshold lost actual source/units/consumer trace");
        const auto effective=root/(rz?"effective_rz":"effective_xy");fs::create_directories(effective);
        for(std::size_t k=0;k<VofEhdInputContract::names.size();++k)write(effective/VofEhdInputContract::names[k],input.fluid.documents[k].effectiveText);
        const auto again=readVofElectrohydrodynamicsInput(effective.string());++readers;
        require(again.chargeDensityJumpTolerance==input.chargeDensityJumpTolerance,"charge regrid effective dictionary does not roundtrip");
    }
    const auto legacy=readVofElectrohydrodynamicsInput(Fixture{}.save(root,"no_electric_policy").string());++readers;
    require(!legacy.chargeDensityJumpTolerance&&legacy.fluid.documents[6].effectiveText.find("Regrid")==std::string::npos,
        "omitted charge policy fabricated an unused setting");
    auto rejects=[&](Fixture f,const std::string& name,const char* field){
        bool caught=false;try{readVofElectrohydrodynamicsInput(f.save(root,name).string());}
        catch(const std::exception& e){caught=std::string(e.what()).find(field)!=std::string::npos;}
        require(caught,"invalid charge regrid input accepted or field provenance lost");++guards;
    };
    for(const auto* value:{"0","-.1","nan","inf"}){Fixture f;option(f,value);rejects(f,std::string("invalid_")+value,"Electric.Regrid.chargeDensityJumpTolerance");}
    {Fixture f;option(f,".1");replace(f.text[6],"chargeDensityJumpTolerance","unknownTolerance");rejects(f,"unknown","Electric.Regrid.unknownTolerance");}
    {Fixture f;option(f,".1");f.text[2]="Mesh { resolution (8 8); initialMesh uniform; }";rejects(f,"no_amr","Electric.Regrid");}
    {Fixture f;option(f,".1");f.text[2]="Mesh { resolution (8 8); initialMesh uniform; }";
        replace(f.text[6],"vofFiniteVolume","vofSharpSurface");replace(f.text[6],"chargeDensity \"","surfaceChargeDensity \"");rejects(f,"b1_rejected","Electric.Regrid");}
}
enum class Edge { None,BelowMinimum,ZeroUnderflow };
std::size_t actual(const fs::path& root,bool rz,bool protectedCharge,Edge edge=Edge::None){
    const bool belowMinimum=edge==Edge::BelowMinimum,zeroUnderflow=edge==Edge::ZeroUnderflow;
    const int resolution=belowMinimum?4:8,maximumLevel=belowMinimum?3:1;
    Fixture f(rz);f.text[0]="Run { app VofElectrohydrodynamics; } Parallel { threads 1; } Time { units code; dt .00001; endTime .00001; } Output { Fields { enabled false; } }";
    f.text[1]=std::string("Geometry { coordinates ")+(rz?"rz":"xy")+"; Domain { min (0 0); max (1 1); } Phases { negative 7; positive 3; } Constants { chargeScale .03125; } Initial { phi \"1\"; } }";
    f.text[2]=belowMinimum?
        "Mesh { resolution (4 4); initialMesh uniform; Regrid { minLevel 2; maxLevel 3; interfaceBandLayers 0; vorticityTolerance .5; coarsenRatio .25; } }":
        "Mesh { resolution (8 8); initialMesh uniform; Regrid { minLevel 0; maxLevel 1; interfaceBandLayers 0; vorticityTolerance .5; coarsenRatio .25; } }";
    f.text[6]=std::string("Electric { method vofFiniteVolume; Initial { chargeDensity \"chargeScale*(")+(rz?"r":"x")+"-.75)*exp(-128*("+(rz?"r":"x")+"-.75)^2)\"; } Boundary { W { type Neumann; value 0; } S { type Dirichlet; value 0; } E { type Neumann; value 0; } N { type Dirichlet; value 0; } } }";
    if(zeroUnderflow){
        replace(f.text[6],std::string("chargeScale*(")+(rz?"r":"x")+
            "-.75)*exp(-128*("+(rz?"r":"x")+"-.75)^2)","0");
        replace(f.text[1],"Constants { chargeScale .03125; } ","");
        replace(f.text[2],"coarsenRatio .25;","coarsenRatio 1e-20;");
    }
    if(protectedCharge)option(f,zeroUnderflow?"2.2250738585072014e-308":".0001");
    const std::string name=std::string(rz?"rz_":"xy_")+(belowMinimum?"below_minimum":zeroUnderflow?"zero_underflow":protectedCharge?"protected":"baseline");
    const auto run=root/name;fs::create_directories(run);
    const auto input=readVofElectrohydrodynamicsInput(f.save(run,"system").string());
    if(zeroUnderflow){
        require(input.chargeDensityJumpTolerance==std::numeric_limits<double>::min()
            &&input.fluid.vofRegrid->coarsenRatio==1e-20,
            "zero-charge underflow fixture lost its explicit minimum-normal threshold and ratio");
        volatile double nativeLow=*input.chargeDensityJumpTolerance*input.fluid.vofRegrid->coarsenRatio;
        require(nativeLow==0,"zero-charge fixture did not exercise native hysteresis-product underflow");
    }
    DictionaryMesh d;d.Nx=resolution;d.Ny=resolution;d.minCoord={0,0};d.maxCoord={1,1};d.maxInitialIterations=0;d.maxCellNum=4096;d.maxLevel_local={maximumLevel};
    QuadtreeMesh mesh(d);
    if(!belowMinimum){DictionaryMesh::RefineBox box;box.level=1;box.min={0,0};box.max={1,1};mesh.applyRefineBoxes({box});}
    for(int ci=0;ci<mesh.cellIndexNum();++ci)require(mesh.cells(mesh.cellIndexToLabel[ci]).level()==(belowMinimum?0:1),
        "charge regrid boundary fixture has the wrong actual source level");
    VofElectrohydrodynamicsOutput writer(input,{run});VofElectrohydrodynamicsTransaction tx(input,writer);writer.begin();
    const auto pub=I::GeometricVofInitialization::prepare(mesh,rz,[](Vector2D){return 1.;},input.fluid.vofInitialization);
    const auto layout=H::HydroDofLayout::build(mesh,pub.domain());
    const auto metric=H::HydroPhysicalDivergence::build(mesh,H::HydroIntegralDivergence::build(mesh,layout));
    const auto geometry=H::HydroPhaseGeometry::build(mesh,metric,pub);
    const char* stage="charge-layer initial";
    tx.initializeVof(mesh,H::HydroTwoPhaseState::createVof(mesh,geometry,H::HydroNativeValues(layout,
        std::vector<double>(layout.velocityCount()),std::vector<double>(layout.pressureCount()))),HydrodynamicsTransaction::Clock::now(),stage);
    const auto oldQ=*tx.acceptedCharge();const auto oldValues=oldQ.integratedCharge();
    const LAMP::TimeLayer target{1,0,.00001,LAMP::TimeLayer::Position::End};
    tx.prepare(mesh,target,input.fluid.step,false,stage);
    require(tx.acceptedCharge()->sameInstance(oldQ)&&tx.acceptedCharge()->integratedCharge()==oldValues,
        "charge AMR preparation mutated accepted Q");
    if(belowMinimum){
        const auto* physical=tx.physicalEndField();require(physical,"below-minimum fixture lacks its actual End field");
        const auto& actualDomain=physical->geometry().domain();const auto& actualQ=physical->chargeCandidate().integratedCharge();
        bool strongJump=false;
        for(const auto face:actualDomain.faces(mesh))if(face.negativeCell>=0&&face.positiveCell>=0){
            const int a=face.negativeCell,b=face.positiveCell;
            const long double rhoA=static_cast<long double>(actualQ[a])/actualDomain.cellVolumes()[actualDomain.cellRowByIndex()[a]],
                rhoB=static_cast<long double>(actualQ[b])/actualDomain.cellVolumes()[actualDomain.cellRowByIndex()[b]];
            strongJump=strongJump||std::abs(rhoB-rhoA)>2*(*input.chargeDensityJumpTolerance);
        }
        require(strongJump,"below-minimum fixture did not exercise a strong actual same-End charge jump");
    }
    tx.commit(mesh);++steps;
    const auto* receipt=tx.remapReceipt();require(receipt&&receipt->decision&&receipt->decision->changed,"charge AMR fixture did not actually adapt");
    const auto& report=*receipt->decision;
    require(report.interfaceSeedCells==0&&report.flowRefineCells==0,"charge layer fixture accidentally refined for interface/vorticity");
    if(edge==Edge::None)
        require(protectedCharge?report.additionalHoldCells>0:report.additionalHoldCells==0,"actual End charge protection request was not consumed");
    bool fineCharge=false,coarseQuiet=false;
    for(int ci:tx.accepted()->phaseGeometry().vofPublication().domain().cellIndexByRow()){
        const auto& c=mesh.cells(mesh.cellIndexToLabel[ci]);const auto p=mesh.verts(c.cellToVert[0]).coord();
        if(p.x>=.625&&p.x<.875&&c.level()==1)fineCharge=true;
        if(p.x<.25&&c.level()==0)coarseQuiet=true;
    }
    if(edge==Edge::None)
        require(coarseQuiet&&(!protectedCharge||fineCharge),"charge protection lost the remote signed layer or refined the quiet domain globally");
    if(belowMinimum)for(int ci=0;ci<mesh.cellIndexNum();++ci)
        require(mesh.cells(mesh.cellIndexToLabel[ci]).level()>=2,
            "same-End charge request failed to clamp a below-minimum source into active policy levels");
    if(zeroUnderflow){
        require(report.additionalRefineCells==0&&report.additionalHoldCells==0&&mesh.cellIndexNum()==resolution*resolution,
            "an underflowed positive hysteresis threshold protected the exact-zero charge domain");
        const auto* accepted=tx.accepted();const auto* acceptedQ=tx.acceptedCharge();const auto* acceptedField=tx.acceptedField();
        require(accepted&&acceptedQ&&acceptedField,"zero-charge accepted transaction is incomplete");
        accepted->validateFor(mesh);const auto& acceptedVof=accepted->phaseGeometry().vofPublication();
        acceptedQ->validateFor(mesh,acceptedVof,target);
        acceptedField->validateFor(mesh,acceptedVof,acceptedField->sourceCharge(),target);
        require(VofEhdJointContract::sameClock(accepted->timeLayer(),target)
            &&VofEhdJointContract::sameClock(acceptedField->candidateTime(),target)
            &&acceptedField->isElectrostaticSnapshot()&&acceptedField->effectiveDt()==0
            &&acceptedField->chargeCandidate().sameInstance(*acceptedQ),
            "zero-charge coarse candidate lost actual same-End geometry/Q/snapshot owners");
        for(int ci=0;ci<mesh.cellIndexNum();++ci)
            require(mesh.cells(mesh.cellIndexToLabel[ci]).level()==0&&acceptedQ->integratedCharge()[ci]==0,
                "zero-charge underflow regression did not truly coarsen all native cells to L0");
    }else VofEhdJointContract::validate(mesh,tx,target);
    const auto count=static_cast<std::size_t>(mesh.cellIndexNum());
    tx.output(mesh,HydrodynamicsTransaction::Clock::now(),0,true);tx.releaseStep();writer.complete();tx.clear();return count;
}
void run(const fs::path& requested){
    require(requested.is_absolute()&&requested.root_name().string()=="D:","charge regrid selector requires an absolute D: output directory");
    const auto root=requested/("charge_regrid_"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    require(!fs::exists(root),"charge regrid output must be fresh");fs::create_directories(root);schema(root);
    for(bool rz:{false,true}){const auto baseline=actual(root,rz,false),protectedLeaves=actual(root,rz,true);
        require(protectedLeaves>baseline&&protectedLeaves<256,"charge criterion failed to preserve a local layer at bounded leaf count");
        actual(root,rz,true,Edge::BelowMinimum);actual(root,rz,true,Edge::ZeroUnderflow);}
    require(readers==5&&guards==7&&steps==8,"charge regrid coverage incomplete");
    std::cout<<"VOF_CHARGE_REGRID_PASS,readers,"<<readers<<",guards,"<<guards<<",actual_End_steps,"<<steps<<",XY_RZ,signed_layer_protected,quiet_domain_coarsened,below_minimum_clamped,zero_underflow_coarsened,no_additional_field_solve\n";
}
} // namespace VofChargeRegridContract
// BEGIN PRESSURE-DRIVEN RESERVOIR READER AND APP CONTRACT
namespace VofPressureDrivenReservoirContract
{
namespace fs=std::filesystem;
namespace H=LAMP::Hydrodynamics;
namespace I=LAMP::Interface;
namespace E=LAMP::Electric;
namespace J=VofEhdJointContract;
namespace O=VofEhdOpenJointContract;
using Fixture=VofEhdInputContract::Fixture;
using Patch=HydrodynamicsInput::VofBoundaryPatch;
using Clock=LAMP::TimeLayer;
using Wall=HydrodynamicsTransaction::Clock;
using VofEhdInputContract::replace;
unsigned checks=0,readers=0,guards=0,accepted=0,changed=0;
void check(bool ok,const char* why){++checks;if(!ok)throw std::runtime_error(why);}
fs::path fresh(const fs::path& requested,const char* name)
{
    check(requested.is_absolute()&&requested.root_name().string()=="D:","pressure-driven selectors require an absolute D: output directory");
    const auto root=requested/(std::string(name)+std::to_string(Wall::now().time_since_epoch().count()));
    check(!fs::exists(root),"pressure-driven output must be fresh");fs::create_directories(root);return root;
}
Fixture fixture(bool rz=false,bool continuation=false)
{
    Fixture f(rz);
    f.text[0]="Run { app VofElectrohydrodynamics; } Parallel { threads 1; } Time { units code; dt .00001; endTime .00004; } Output { Fields { enabled true; writeEverySteps 1; } }";
    f.text[1]=std::string("Geometry { coordinates ")+(rz?"rz":"xy")+"; Domain { min (0 0); max (1 1); } Phases { negative 7; positive 3; } Initial { phi \""+(rz?"z-.5":"y-.5")+"\"; } }";
    f.setConstants("chargeScale .0002;");
    f.text[4]="Materials { Phases { Gas { phaseId 3; density 1; dynamicViscosity 0; permittivity 2; conductivity .1; } Liquid { phaseId 7; density 2; dynamicViscosity 0; permittivity 4; conductivity .25; } } Interfaces { Shared { phases (7 3); surfaceTension 0; } } }";
    f.text[5]=std::string("Hydrodynamics { method staggered; densityRepresentation cellAverageProjection; capillaryModel heightFunction; Initial { velocity (0 0); } Boundary { type pressureDrivenReservoirs; Backflow { ")+
        (continuation?"model interiorContinuationZeroShear;":"model ambient; velocity (0 0);")+" } Patches { supply { side S; interval (0 1); type pressureOutlet; phase liquid; pressure 1; } drain { side N; interval (0 1); type pressureOutlet; phase gas; pressure 0; } } } BodyForce { acceleration (0 0); } TimeIntegration { method adaptiveLinear; } }";
    f.text[6]=std::string("Electric { method vofFiniteVolume; Initial { chargeDensity \"chargeScale*(")+(rz?"z":"y")+"-.5)\"; } Boundary { W { type Neumann; value 0; } S { type Dirichlet; value .01; } E { type Neumann; value 0; } N { type Dirichlet; value 0; } } Advection { Inlets { supply { chargeDensity -.01; } drain { chargeDensity .005; } } } }";
    return f;
}
void readerCase(const fs::path& root,bool rz,bool continuation)
{
    const auto source=fixture(rz,continuation);const std::string name=std::string(rz?"rz_":"xy_")+(continuation?"continuation":"ambient");
    const auto input=readVofElectrohydrodynamicsInput(source.save(root,name).string());const auto& f=input.fluid;
    check(f.vofPressureDrivenReservoirs&&f.radial==rz&&f.vofBoundaryPatches.size()==2&&f.vofPressureOutletBackflow.has_value(),"pressure-driven reader lost explicit mode, ports or Backflow");
    check(f.initialVelocity.x==0&&f.initialVelocity.y==0&&f.vof.surfaceTension==0&&f.run.fixedStepCount()==4,"pressure-driven reader changed resting sigma=0 four-step data");
    check(f.vofBoundaryPatches[0].name=="supply"&&f.vofBoundaryPatches[0].side==Patch::Side::S&&f.vofBoundaryPatches[0].pressure==1
        &&f.vofBoundaryPatches[0].reservoir==I::GeometricVofAdvectorSlab::Boundary::Mode::PureLiquidReservoir
        &&f.vofBoundaryPatches[1].name=="drain"&&f.vofBoundaryPatches[1].side==Patch::Side::N&&f.vofBoundaryPatches[1].pressure==0
        &&f.vofBoundaryPatches[1].reservoir==I::GeometricVofAdvectorSlab::Boundary::Mode::PureGasReservoir,"pressure-driven physical pressure, side, phase or patch order changed");
    for(const auto& p:f.vofBoundaryPatches)check(p.kind==Patch::Kind::PressureOutlet&&p.velocity.x==0&&p.velocity.y==0,"pressure-driven reader created a prescribed feed");
    check(input.inletChargeDensityByPatch==std::vector<double>({-.01,.005})&&input.materials.conductivityLiquid>0&&input.materials.conductivityGas>0,"pressure reservoirs lost independent signed charge or positive conductivity");
    const auto expected=continuation?H::HydroVofPressureOutletMomentumModel::InteriorContinuationZeroShear:H::HydroVofPressureOutletMomentumModel::Ambient;
    check(f.vofPressureOutletBackflow->model==expected&&f.vofPressureOutletBackflow->velocity.has_value()==!continuation,"pressure-driven explicit momentum model changed");
    bool trace=false;for(const auto& e:f.trace.entries)if(e.field=="Hydrodynamics.Boundary.type"){
        const auto* value=std::get_if<std::string>(&e.value);trace=e.used&&!e.defaulted&&e.consumer=="Hydrodynamics.nativeBoundary"
            &&e.units=="dimensionless"&&fs::path(e.file).filename()=="hydro.dict"&&value&&*value=="pressureDrivenReservoirs";}
    check(trace,"pressure-driven type lacks actual native-boundary consumption trace");
    const auto effective=root/(name+"_effective");fs::create_directories(effective);
    for(std::size_t n=0;n<VofEhdInputContract::names.size();++n){check(f.documents[n].sourceText==source.text[n],"pressure-driven reader changed source document bytes");
        VofEhdInputContract::write(effective/VofEhdInputContract::names[n],f.documents[n].effectiveText);}
    check(f.documents[5].effectiveText.find("pressureDrivenReservoirs")!=std::string::npos,"pressure-driven effective input lost selected type");
    const auto again=readVofElectrohydrodynamicsInput(effective.string());
    check(again.fluid.vofPressureDrivenReservoirs&&again.fluid.radial==rz&&again.fluid.vofBoundaryPatches.size()==2
        &&again.fluid.vofBoundaryPatches[0].pressure==1&&again.fluid.vofBoundaryPatches[1].pressure==0
        &&again.fluid.vofPressureOutletBackflow&&again.fluid.vofPressureOutletBackflow->model==expected
        &&again.inletChargeDensityByPatch==input.inletChargeDensityByPatch,"pressure-driven effective roundtrip changed physical mode or reservoirs");++readers;
}
template<class Change>void reject(const fs::path& root,const char* name,const char* field,Change change,bool rz=false)
{
    auto f=fixture(rz);change(f);bool caught=false;
    try{(void)readVofElectrohydrodynamicsInput(f.save(root,name).string());}
    catch(const std::exception& e){caught=true;check(std::string(e.what()).find(field)!=std::string::npos,"pressure-driven guard lost actual source field");}
    check(caught,"invalid pressure-driven reservoir input accepted");++guards;
}
void reader(const fs::path& requested)
{
    const auto root=fresh(requested,"pressure_driven_reader_");for(bool rz:{false,true})for(bool continuation:{false,true})readerCase(root,rz,continuation);
    check(!HydrodynamicsInput{}.vofPressureDrivenReservoirs,"pressure-driven mode does not default false");
    {Fixture f;const auto input=readVofElectrohydrodynamicsInput(f.save(root,"legacy_closed").string());check(!input.fluid.vofPressureDrivenReservoirs&&input.fluid.vofBoundaryPatches.empty(),"legacy closed input acquired pressure supply");}
    {auto f=VofEhdOpenInputContract::open();const auto input=readVofElectrohydrodynamicsInput(f.save(root,"legacy_prescribed").string());check(!input.fluid.vofPressureDrivenReservoirs&&input.fluid.vofBoundaryPatches[0].kind==Patch::Kind::PrescribedVelocity,"legacy prescribedFlow changed mode");}
    reject(root,"missing_backflow","Boundary",[](Fixture& f){replace(f.text[5]," Backflow { model ambient; velocity (0 0); }","");});
    reject(root,"equal_pressure","Boundary.Patches",[](Fixture& f){replace(f.text[5],"pressure 1;","pressure 0;");});
    reject(root,"one_port","Boundary.Patches",[](Fixture& f){replace(f.text[5]," drain { side N; interval (0 1); type pressureOutlet; phase gas; pressure 0; }","");});
    reject(root,"fixed_feed","Boundary.Patches.supply",[](Fixture& f){replace(f.text[5],"type pressureOutlet; phase liquid; pressure 1;","type prescribedVelocity; phase liquid; velocity (0 .02);");});
    reject(root,"velocity_payload","Boundary.Patches.supply.velocity",[](Fixture& f){replace(f.text[5],"pressure 1;","pressure 1; velocity (0 .02);");});
    reject(root,"missing_phase","Boundary.Patches.supply",[](Fixture& f){replace(f.text[5]," phase liquid;","");});
    reject(root,"missing_pressure","Boundary.Patches.supply",[](Fixture& f){replace(f.text[5]," pressure 1;","");});
    reject(root,"nonfinite_pressure","Boundary.Patches.supply.pressure",[](Fixture& f){replace(f.text[5],"pressure 1;","pressure nan;");});
    reject(root,"unaligned","Boundary.Patches.supply.interval",[](Fixture& f){replace(f.text[5],"supply { side S; interval (0 1)","supply { side S; interval (0 .3)");});
    reject(root,"reversed","Boundary.Patches.supply.interval",[](Fixture& f){replace(f.text[5],"supply { side S; interval (0 1)","supply { side S; interval (1 0)");});
    reject(root,"overlap","Boundary.Patches.drain.interval",[](Fixture& f){replace(f.text[5],"drain { side N;","drain { side S;");});
    reject(root,"axis","Boundary.Patches.supply.side",[](Fixture& f){replace(f.text[5],"supply { side S;","supply { side W;");},true);
    reject(root,"nonzero_initial","Initial.velocity",[](Fixture& f){replace(f.text[5],"Initial { velocity (0 0); }","Initial { velocity (.01 0); }");});
    reject(root,"ambient_missing_vector","Backflow",[](Fixture& f){replace(f.text[5],"model ambient; velocity (0 0);","model ambient;");});
    reject(root,"continuation_unused_vector","Backflow.velocity",[](Fixture& f){replace(f.text[5],"model ambient;","model interiorContinuationZeroShear;");});
    check(readers==4&&guards==15,"pressure-driven reader coverage incomplete");
    std::cout<<"VOF_PRESSURE_DRIVEN_READER_PASS readers="<<readers<<" guards="<<guards<<" checks="<<checks<<" reader_effective_trace_only output="<<root.string()<<'\n';
}
void native(const QuadtreeMesh& mesh,const H::HydroTwoPhaseState& state,bool initial)
{
    state.validateFor(mesh);const auto* flow=state.vofPrescribedFlow();
    check(flow&&flow->pressureDrivenReservoirs&&flow->prescribedVelocityDofs.empty()&&!flow->endpointVelocity
        &&flow->pressureOutletBackflowModel&&flow->stationaryNoSlipVelocityDofs.empty(),"pressure-driven App fabricated a fixed feed/vector or lost explicit momentum model");
    check(std::all_of(flow->normalVelocity.begin(),flow->normalVelocity.end(),[](double x){return x==0;}),"pressure-driven free-port normal payload is not exactly zero");
    const auto& layout=state.values().layout();const auto& domain=layout.domain();const auto& u=state.values().normalVelocity();
    const auto& D=state.phaseGeometry().metrics().divergence().matrix();
    for(int c=0;c<D.rows();++c){long double error=0,scale=0;for(H::HydroIntegralDivergence::Matrix::InnerIterator it(D,c);it;++it){const long double term=static_cast<long double>(it.value())*u[it.col()];error+=term;scale+=std::abs(term);}O::bound(error,scale,1e-11L,"pressure-driven actual native D_I U");}
    long double total=0,activity=0,supply=0,drain=0;unsigned ports=0;
    for(int row=0;row<domain.nativeFaceCount();++row){const auto face=domain.face(mesh,row);if(face.negativeCell>=0&&face.positiveCell>=0)continue;
        const int dof=layout.velocityDof(row);const auto mode=flow->rulesByNativeRow.at(row);
        if(face.area==0){check(domain.axisymmetric()&&mode==I::GeometricVofAdvectorSlab::Boundary::Mode::Closed,"pressure-driven flow opened true axis");if(dof>=0)check(u[dof]==0,"pressure-driven axis normal nonzero");continue;}
        check(dof>=0,"pressure-driven positive-area exterior face lacks native DOF");const long double flux=static_cast<long double>(face.outwardSign)*face.area*u[dof];total+=flux;activity+=std::abs(flux);
        if(mode==I::GeometricVofAdvectorSlab::Boundary::Mode::Closed){check(u[dof]==0,"pressure-driven wall normal changed");continue;}
        check(face.horizontal,"pressure-driven test port is not an actual S/N face");const bool south=face.negativeCell<0;
        check(mode==(south?I::GeometricVofAdvectorSlab::Boundary::Mode::PureLiquidReservoir:I::GeometricVofAdvectorSlab::Boundary::Mode::PureGasReservoir),"pressure-driven port phase changed");
        bool found=false;for(const auto& port:flow->pressureOutlets)if(port.velocityDof==dof){check(port.pressure==(south?1:0),"pressure-driven physical pressure changed or became remap multiplier");found=true;}
        check(found,"pressure-driven free port absent from actual Saddle");++ports;(south?supply:drain)+=flux;
    }
    check(ports==flow->pressureOutlets.size()&&ports>1,"pressure-driven actual free-port coverage incomplete");
    O::bound(total,activity,1e-12L,"pressure-driven independent exterior flux balance");
    if(initial){check(std::all_of(u.begin(),u.end(),[](double x){return x==0;}),"pressure-driven Initial hook fabricated nonzero velocity");
        check(std::all_of(state.values().pressure().begin(),state.values().pressure().end(),[](double p){return p==0;}),"pressure-driven Initial projection multiplier became physical pressure");}
    else check(supply<0&&drain>0,"pressure-driven actual End did not supply from high to low physical pressure");
}
void actual(const fs::path& root,bool rz)
{
    const auto run=root/(rz?"actual_rz":"actual_xy");fs::create_directories(run);const auto input=readVofElectrohydrodynamicsInput(fixture(rz).save(run,"system").string());const auto& f=input.fluid;
    check(f.vofPressureDrivenReservoirs&&f.vofRegrid&&f.vofRegrid->maxLevel==1&&f.run.fixedStepCount()==4
        &&f.vof.surfaceTension==0&&f.vof.bodyAcceleration.x==0&&f.vof.bodyAcceleration.y==0&&!f.vofSouthEven90Interval,"pressure-driven actual App fixture lost four-step sigma=0 scope");
    DictionaryMesh d;d.Nx=f.nx;d.Ny=f.ny;d.minCoord=f.domainMin;d.maxCoord=f.domainMax;d.maxInitialIterations=0;d.maxCellNum=16384;d.maxLevel_local={1};
    auto mesh=std::make_unique<QuadtreeMesh>(d);DictionaryMesh::RefineBox fine;fine.level=1;fine.min=f.domainMin;fine.max=f.domainMax;mesh->applyRefineBoxes({fine});
    VofElectrohydrodynamicsOutput writer(input,{run});VofElectrohydrodynamicsTransaction tx(input,writer);writer.begin();const char* stage="pressure-driven actual Initial";
    {const auto pub=I::GeometricVofInitialization::prepare(*mesh,rz,[&](Vector2D p){return f.initialJet(p).value;},f.vofInitialization);
        const auto layout=H::HydroDofLayout::build(*mesh,pub.domain());const auto metric=H::HydroPhysicalDivergence::build(*mesh,H::HydroIntegralDivergence::build(*mesh,layout));
        auto initial=tx.makeInitialVofState(*mesh,H::HydroPhaseGeometry::build(*mesh,metric,pub),H::HydroNativeValues(layout,std::vector<double>(layout.velocityCount()),std::vector<double>(layout.pressureCount())));
        native(*mesh,initial,true);tx.initializeVof(*mesh,std::move(initial),Wall::now(),stage);}
    J::validate(*mesh,tx,{0,0,0,Clock::Position::Initial});unsigned localChanged=0;bool realE=false,realJ=false;long double pressureWork=0;
    for(int index=0;index<4;++index){const auto previous=tx.accepted()->timeLayer();const auto sourceQ=*tx.acceptedCharge();const J::Snapshot before(*mesh,tx,run);
        const Clock target{previous.step+1,previous.time(),f.run.stepEnd(previous.time(),f.run.dt,static_cast<int>(previous.step)),Clock::Position::End};const auto start=Wall::now();
        tx.prepare(*mesh,target,f.step,false,stage);before.unchanged(*mesh,tx,run);
        check(tx.physicalChargeTransport()&&tx.physicalEndField()&&!tx.remapReceipt(),"pressure-driven prepare lost genuine physical receipts or exposed uncommitted remap");
        tx.commit(*mesh);++accepted;J::validate(*mesh,tx,target);native(*mesh,*tx.accepted(),false);J::checkReceipt(*mesh,tx,previous,target);
        const auto* receipt=tx.remapReceipt();check(receipt&&receipt->physicalStep&&receipt->physicalEnd&&receipt->sourceMesh&&receipt->decision,"pressure-driven physical/same-End receipt missing");
        const auto& oldMesh=*receipt->sourceMesh;const auto& step=*receipt->physicalStep;const auto& physical=*receipt->physicalEnd;
        O::phaseLedger(oldMesh,step.transport(),true);O::phaseLedger(oldMesh,step.transport(),false);
        O::physicalCharge(oldMesh,sourceQ,step.transport(),*tx.physicalChargeTransport(),*tx.physicalEndField());
        check(!step.report().stages.empty()&&step.frozenMidpoint().hasMaxwellLoad()&&!step.hasCapillaryLoad(),"pressure-driven actual A lost Maxwell or inserted capillary force");
        const auto& final=step.report().stages.back();const auto* flow=step.beginState().vofPrescribedFlow();const auto& layout=step.beginState().values().layout();const auto& mid=step.solvedMidpointValues().normalVelocity();long double work=0,workScale=0;
        for(const auto& port:flow->pressureOutlets){const auto face=layout.domain().face(oldMesh,layout.nativeFace(port.velocityDof));const long double term=-static_cast<long double>(target.end-target.begin)*port.pressure*face.outwardSign*face.area*mid[port.velocityDof];work+=term;workScale+=std::abs(term);}
        O::bound(final.pressureOutletWork-work,std::abs(final.pressureOutletWork)+workScale,1e-11L,"pressure-driven actual pressure traction work disagrees with physical free-port flux");
        check(work>0&&final.prescribedReactionWork==0,"pressure-driven actual flow lacks positive physical pressure work or includes prescribed-feed work");pressureWork+=work;
        for(const auto& face:tx.physicalEndField()->faces())if(face.physicalMeasure>0){check(face.conductivity>0,"pressure-driven actual Ohmic face lost positive K");realE=realE||face.integratedField!=0;realJ=realJ||face.current!=0;}
        if(receipt->decision->changed){++localChanged;++changed;check(receipt->fluidTransfer&&receipt->interfaceTransfer,"pressure-driven changed End missing actual transfers");const auto& remap=receipt->fluidTransfer->report();
            check(remap.prescribedFlow&&remap.fixedNormals==0&&remap.pressureOutletNormals>0&&!remap.pressureGaugeApplied&&remap.projectionEnergyIdentityRelativeResidual<=2e-11L,"pressure-driven same-End remap fixed a reservoir, gauged physical p or broke original projection gate");
            const auto a=O::stocks(physical.phaseGeometry().vofPublication()),b=O::stocks(tx.accepted()->phaseGeometry().vofPublication());for(int k=0;k<2;++k)O::bound(b[k]-a[k],std::abs(a[k])+std::abs(b[k]),1e-12L,"pressure-driven same-End independent phase stock changed");
            long double q0=0,q1=0,scale=0;for(double q:tx.physicalEndField()->chargeCandidate().integratedCharge()){q0+=q;scale+=std::abs(q);}for(double q:tx.acceptedCharge()->integratedCharge()){q1+=q;scale+=std::abs(q);}O::bound(q1-q0,scale,1e-10L,"pressure-driven same-End signed Q changed");
            check(tx.acceptedField()->isElectrostaticSnapshot()&&tx.acceptedField()->effectiveDt()==0&&!tx.acceptedField()->advection(),"pressure-driven target snapshot masquerades as physical End current");}
        tx.output(*mesh,start,std::chrono::duration<double>(Wall::now()-start).count(),index==3);tx.releaseStep();check(!tx.physicalChargeTransport()&&!tx.physicalEndField(),"pressure-driven release retained old physical receipts");
    }
    check(localChanged>0&&realE&&realJ&&pressureWork>0&&tx.accepted()->timeLayer().time()==f.run.endTime,"pressure-driven four-step actual E/J/pressure/AMR scope incomplete");
    writer.complete();check(J::read(run/"output/meta/manifest.json").find("\"status\":\"COMPLETED\"")!=std::string::npos,"pressure-driven completed manifest absent");const auto frames=J::read(run/"output/meta/frames.csv");check(std::count(frames.begin(),frames.end(),'\n')==6,"pressure-driven frame index is not Initial plus four unique accepted frames");
    tx.clear();check(!tx.accepted()&&!tx.acceptedCharge()&&!tx.acceptedField(),"pressure-driven clear retained accepted subset");
}
void app(const fs::path& requested)
{
    const auto root=fresh(requested,"pressure_driven_app_");actual(root,false);actual(root,true);
    check(accepted==8&&changed>=2,"pressure-driven actual XY/RZ four-step changed-AMR coverage incomplete");
    std::cout<<"VOF_PRESSURE_DRIVEN_APP_PASS XY_RZ accepted="<<accepted<<" changed="<<changed<<" checks="<<checks<<" actual_A_sigma0_positive_K_E_J_pressure_work_same_End_AMR_no_contact_or_noSlip_qualification output="<<root.string()<<'\n';
}
}
// END PRESSURE-DRIVEN RESERVOIR READER AND APP CONTRACT

// BEGIN CSG INITIAL PUBLIC RUN CONTRACT
#include "Core/IO/DictionaryGeometry.h"
namespace VofCsgInitialAppContract
{
namespace fs=std::filesystem;
namespace H=LAMP::Hydrodynamics;
namespace I=LAMP::Interface;
namespace J=VofEhdJointContract;
using Fixture=VofEhdInputContract::Fixture;
using Time=LAMP::TimeLayer;
using VofEhdInputContract::replace;
unsigned checks=0,guards=0,observedSteps=0;
void check(bool ok,const char* why){++checks;if(!ok)throw std::runtime_error(why);}
fs::path originalCsg()
{
    auto root=fs::path(__FILE__);for(int k=0;k<4;++k)root=root.parent_path();
    const auto source=root/"cases/electrospray/system/nozzle.csg";
    check(fs::is_regular_file(source),"CSG contract cannot locate original electrospray geometry asset");return source;
}
Fixture fixture()
{
    Fixture f(true);
    f.text[0]="Run { app VofElectrohydrodynamics; } Parallel { threads 1; } Time { units code; dt .00001; endTime .00002; } Output { Fields { enabled true; writeEverySteps 1; } }";
    f.text[1]="Geometry { coordinates rz; Domain { min (0 0); max (5 16); } Phases { negative 2; positive 3; } Initial { csgFile \"nozzle.csg\"; lengthScaleSI 1e-6; excludedPhases (1); Pin { pointId 6; wetBoundaryId 3; dryBoundaryId 2; } } }";
    f.text[2]="Mesh { resolution (40 80); initialMesh uniform; }";
    f.text[4]="Materials { Phases { Gas { phaseId 3; density 1; dynamicViscosity 0; permittivity 1; conductivity .1; } Liquid { phaseId 2; density 2; dynamicViscosity 0; permittivity 1; conductivity .1; } } Interfaces { Shared { phases (2 3); surfaceTension 0; } } }";
    f.text[5]="Hydrodynamics { method staggered; densityRepresentation cellAverageProjection; capillaryModel heightFunction; Initial { velocity (0 0); } Boundary { type impermeable; } BodyForce { acceleration (0 0); } TimeIntegration { method adaptiveLinear; } }";
    f.text[6]="Electric { method vofFiniteVolume; Initial { chargeDensity 0; } Boundary { Patches { electrodeDry { boundaryId 2; type Dirichlet; value 0; } electrodeWet { boundaryId 3; type Dirichlet; value 0; } vacuumTop { boundaryId 9; type Dirichlet; value 0; } liquidInlet { boundaryId 5; type Neumann; value 0; } vacuumOuter { boundaryId 8; type Neumann; value 0; } vacuumAxis { boundaryId 10; type Neumann; value 0; } liquidAxisUpper { boundaryId 11; type Neumann; value 0; } liquidAxisLower { boundaryId 12; type Neumann; value 0; } } } }";
    return f;
}
fs::path save(const fs::path& root,const std::string& name,const Fixture& f)
{
    const auto system=f.save(root,name);fs::copy_file(originalCsg(),system/"nozzle.csg");return system;
}
void typed(const VofElectrohydrodynamicsInput& input,const std::string& raw,bool expectedCheckpointDefault=true)
{
    const auto& f=input.fluid;check(f.csgGeometry&&f.csgFile=="nozzle.csg"&&f.csgLengthScaleSI==1e-6&&f.csgExcludedPhases==std::vector<int>({1}),"CSG reader lost typed raw asset, L0 or excluded physical phase IDs");
    check(f.radial&&f.nx==40&&f.ny==80&&f.domainMin.x==0&&f.domainMin.y==0&&f.domainMax.x==5&&f.domainMax.y==16
        &&f.negativePhase==2&&f.positivePhase==3&&!f.initialInterfaceBand&&!f.vofRegrid&&f.run.checkpointEverySteps==0&&f.run.restartDirectory.empty(),"CSG reader changed original units, physical phases or uniform/no-checkpoint scope");
    check(f.csgPin&&f.csgPin->pointId==6&&f.csgPin->wetBoundaryId==3&&f.csgPin->dryBoundaryId==2,"CSG reader lost explicit original point/wet/dry pin IDs");
    check(f.csgGeometry->points.size()==9&&f.csgGeometry->points[5].x==1&&f.csgGeometry->points[5].y==3
        &&f.csgGeometry->domainMax.x==5&&f.csgGeometry->domainMax.y==16&&f.csgGeometry->phaseName.size()==3,"CSG L0 conversion changed original nozzle lip, domain or phases");
    bool noPhi=false;try{(void)f.initialJet({.5,2});}catch(const std::exception&){noPhi=true;}check(noPhi,"CSG input exposed a fabricated phi expression/jet");
    check(f.documents.size()==8,"CSG case must retain seven main documents and one raw asset");unsigned assets=0,main=0;
    for(const auto& doc:f.documents){const auto name=fs::path(doc.file).filename().string();if(name=="nozzle.csg"){++assets;check(doc.sourceText==raw&&doc.effectiveText==raw,"CSG source/effective asset bytes changed");}
        else{check(std::find(VofEhdInputContract::names.begin(),VofEhdInputContract::names.end(),name)!=VofEhdInputContract::names.end(),"CSG input retained an undeclared extra document");++main;}}
    check(assets==1&&main==7,"CSG asset/main provenance is missing or duplicated");
    bool fileTrace=false,lengthTrace=false,excludedTrace=false,checkpointTrace=false;
    for(const auto& e:f.trace.entries){if(e.field=="Geometry.Initial.csgFile")fileTrace=e.used&&!e.defaulted&&fs::path(e.file).filename()=="geometry.dict";
        if(e.field=="Geometry.Initial.lengthScaleSI")lengthTrace=e.used&&!e.defaulted;
        if(e.field=="Geometry.Initial.excludedPhases")excludedTrace=e.used&&!e.defaulted;
        if(e.field=="Checkpoint.everySteps"){const auto* value=std::get_if<int>(&e.value);checkpointTrace=e.used&&(e.defaulted==expectedCheckpointDefault)&&value&&*value==0;}}
    check(fileTrace&&lengthTrace&&excludedTrace&&checkpointTrace,"CSG actual asset/L0/solid/default-zero-checkpoint trace missing");
}
template<class Change>void reject(const fs::path& root,const char* name,const char* field,Change change)
{
    auto f=fixture();change(f);const auto path=save(root,name,f);bool caught=false;
    try{(void)readVofElectrohydrodynamicsInput(path.string());}catch(const std::exception& e){caught=true;const std::string why=e.what();
        check(why.find(field)!=std::string::npos||why.find("Geometry.Initial")!=std::string::npos,"CSG reader guard lost actual source field");}
    check(caught,"unsupported or malformed CSG Initial input accepted");++guards;
}
void readerGuards(const fs::path& root)
{
    reject(root,"phi_and_csg","Initial",[](Fixture& f){replace(f.text[1],"Initial {","Initial { phi \"z-3\";");});
    reject(root,"missing_length","lengthScaleSI",[](Fixture& f){replace(f.text[1]," lengthScaleSI 1e-6;","");});
    reject(root,"zero_length","lengthScaleSI",[](Fixture& f){replace(f.text[1],"lengthScaleSI 1e-6","lengthScaleSI 0");});
    reject(root,"traversal","csgFile",[](Fixture& f){replace(f.text[1],"nozzle.csg","../nozzle.csg");});
    reject(root,"missing_excluded","excludedPhases",[](Fixture& f){replace(f.text[1]," excludedPhases (1);","");});
    reject(root,"duplicate_excluded","excludedPhases",[](Fixture& f){replace(f.text[1],"excludedPhases (1)","excludedPhases (1 1)");});
    reject(root,"invalid_excluded","excludedPhases",[](Fixture& f){replace(f.text[1],"excludedPhases (1)","excludedPhases (0)");});
    reject(root,"initial_band","Mesh",[](Fixture& f){f.text[2]="Mesh { resolution (40 80); initialMesh interfaceBand; Refinement { level 1; bandCells 2; } }";});
    reject(root,"unfitted_base","selected solid cuts a leaf",[](Fixture& f){f.text[2]="Mesh { resolution (20 40); initialMesh uniform; }";});
    reject(root,"checkpoint","Checkpoint",[](Fixture& f){f.text[0]+=" Checkpoint { everySteps 1; }";});
    reject(root,"restart","Restart",[](Fixture& f){f.text[0]+=" Restart { directory \"D:/not-used-csg-checkpoint\"; }";});
    reject(root,"missing_pin_leg","Pin",[](Fixture& f){replace(f.text[1]," dryBoundaryId 2;","");});
    {const auto path=save(root,"missing_asset",fixture());fs::remove(path/"nozzle.csg");bool caught=false;try{(void)readVofElectrohydrodynamicsInput(path.string());}
        catch(const std::exception& e){caught=true;check(std::string(e.what()).find("nozzle.csg")!=std::string::npos,"missing CSG asset lost actual filename");}check(caught,"CSG reader accepted missing explicit raw asset");++guards;}
    {Fixture old;const auto input=readVofElectrohydrodynamicsInput(old.save(root,"legacy_phi").string());check(!input.fluid.csgGeometry&&input.fluid.csgExcludedPhases.empty()&&!input.fluid.csgPin&&input.fluid.initialJet({.5,.5}).value<0,"legacy phi reader path changed under CSG opt-in");}
}
std::vector<std::array<long double,2>> stocks;
void accepted(const VofElectrohydrodynamicsTransaction& tx,const Time& clock)
{
    const auto* h=tx.accepted();const auto* q=tx.acceptedCharge();const auto* e=tx.acceptedField();check(h&&q&&e,"CSG actual public run accepted an incomplete joint bundle");
    const auto& p=h->phaseGeometry().vofPublication();const auto& domain=p.domain();const auto& map=domain.cellRowByIndex();const auto& cells=p.plicCells();
    check(domain.axisymmetric()&&domain.uniformLeaves()&&map.size()==3200&&domain.cellIndexByRow().size()==2720&&cells.size()==3200
        &&q->integratedCharge().size()==3200&&e->potential().size()==3200,"CSG actual public run did not build original 3200-native/2720-fluid RZ domain");
    check(J::sameClock(clock,h->timeLayer())&&J::sameClock(clock,p.timeLayer())&&J::sameClock(clock,q->timeLayer())&&J::sameClock(clock,e->candidateTime())
        &&e->geometry().identity()==p.identity()&&e->chargeCandidate().sameInstance(*q),"CSG accepted owners or complete End clock differ");
    const auto pin=p.pinnedCorner();check(pin&&pin->sourceIdentity()==p.identity()&&J::sameClock(pin->sourceTime(),clock)&&pin->physicalPoint().x==1&&pin->physicalPoint().y==3
        &&pin->nativeCell()>=0&&map.at(pin->nativeCell())>=0&&pin->wetNativeRow()>=0&&pin->dryNativeRow()>=0&&pin->wetNativeRow()!=pin->dryNativeRow(),"CSG actual accepted publication lacks original typed pin certificate");
    check(!h->vofPrescribedFlow()&&std::all_of(h->values().normalVelocity().begin(),h->values().normalVelocity().end(),[](double x){return x==0;})
        &&std::all_of(h->values().pressure().begin(),h->values().pressure().end(),[](double x){return x==0;}),"CSG sigma=0 closed zero-forcing run generated mechanical motion");
    unsigned solid=0,fluid=0;long double volume=0,liquid=0,gas=0;std::vector<std::array<long double,2>> current;current.reserve(cells.size());
    for(std::size_t ci=0;ci<cells.size();++ci){const auto& cell=cells[ci];current.push_back(cell.inventory);check(q->integratedCharge()[ci]==0&&e->potential()[ci]==0,"CSG zero source/voltage generated Q or phi");
        if(map[ci]<0){++solid;check(cell.inventory[0]==0&&cell.inventory[1]==0&&cell.measured[0]==0&&cell.measured[1]==0&&!cell.resolved
            &&p.cells()[ci].nativeCell==-1&&p.cells()[ci].capacity==0&&p.alpha()[ci]==0&&e->epsilonByCell()[ci]==0&&e->conductivityByCell()[ci]==0,"CSG solid slot acquired fake fluid/gas inventory, Q, field material or geometry");continue;}
        ++fluid;const auto row=static_cast<std::size_t>(map[ci]);const long double cap=domain.cellVolumes().at(row);const auto center=(cell.rectangle.lower+cell.rectangle.upper)*.5;
        check(!(center.x>1&&center.y<3),"CSG active fluid contains an original electrode-solid leaf");
        check(std::isfinite(cell.inventory[0])&&std::isfinite(cell.inventory[1])&&cell.inventory[0]>=0&&cell.inventory[1]>=0&&cap>0,"CSG actual independent phase stock is invalid");
        VofEhdOpenJointContract::bound(cell.inventory[0]+cell.inventory[1]-cap,cap,1e-12L,"CSG authoritative liquid/gas capacity closure");volume+=cap;liquid+=cell.inventory[0];gas+=cell.inventory[1];
        if(center.x<.8&&center.y<2.8)check(cell.inventory[0]>0&&cell.inventory[1]==0,"CSG original liquid column lost its phase");
        if(center.y>4.2)check(cell.inventory[0]==0&&cell.inventory[1]>0,"CSG original upper vacuum lost its phase");}
    check(solid==480&&fluid==2720&&liquid>0&&gas>0,"CSG original solid/fluid or two-stock count missing");
    VofEhdOpenJointContract::bound(volume-328*std::acos(-1.L),328*std::acos(-1.L),1e-12L,"CSG original true revolved fluid-domain volume");
    for(const auto& face:e->faces())check(face.integratedField==0&&face.displacement==0&&face.current==0,"CSG zero electrical data generated E/D/J");
    if(stocks.empty())stocks=std::move(current);else check(current==stocks,"CSG zero-motion actual public steps changed independent L/G inventories");++observedSteps;
}
void run(const fs::path& requested)
{
    check(requested.is_absolute()&&requested.root_name().string()=="D:","CSG selector requires absolute D: output root");
    const auto root=requested/("csg_initial_"+std::to_string(HydrodynamicsTransaction::Clock::now().time_since_epoch().count()));check(!fs::exists(root),"CSG output must be fresh");fs::create_directories(root);
    const auto raw=J::read(originalCsg());const auto source=fixture();const auto system=save(root,"system",source);const auto input=readVofElectrohydrodynamicsInput(system.string());typed(input,raw);
    const auto effective=root/"effective";fs::create_directories(effective);for(const auto& doc:input.fluid.documents){const auto name=fs::path(doc.file).filename();VofEhdInputContract::write(effective/name,doc.effectiveText);
        if(name!="nozzle.csg"){const auto found=std::find(VofEhdInputContract::names.begin(),VofEhdInputContract::names.end(),name.string());check(found!=VofEhdInputContract::names.end()&&doc.sourceText==source.text[std::distance(VofEhdInputContract::names.begin(),found)],"CSG reader rewrote source main document");}}
    const auto again=readVofElectrohydrodynamicsInput(effective.string());typed(again,raw,false);readerGuards(root);
    VofElectrohydrodynamicsOutput writer(input,{root});VofElectrohydrodynamicsTransaction tx(input,writer);unsigned failures=0;std::string failure;
    const HydrodynamicsRunIO io{[&]{writer.begin();},[&](auto attempted,auto step,double time,double oldDt,double newDt,const std::string& reason){writer.retry(attempted,step,time,oldDt,newDt,reason);},
        [&](auto attempted,auto step,double time,const char* stage,const std::string& reason){++failures;failure=std::string(stage)+": "+reason;writer.failure(attempted,step,time,stage,reason);},
        [&]{check(observedSteps==2,"CSG actual public run skipped accepted step observations");writer.complete();},
        [&](const auto& report,const auto& clock,int retries){writer.timeStep(report,clock,retries);},
        [&](const auto& clock,double seconds){accepted(tx,clock);writer.stepCompleted(clock,seconds);}};
    RunConfig config;config.app="VofElectrohydrodynamics";const auto result=runHydrodynamics(input.fluid,config,tx,io);
    if(result!=0)throw std::runtime_error("CSG actual public run failed: "+failure);
    check(failures==0&&observedSteps==2&&!tx.accepted()&&!tx.acceptedCharge()&&!tx.acceptedField(),"CSG public run did not complete/clear its two accepted steps");
    check(J::read(root/"output/meta/manifest.json").find("\"status\":\"COMPLETED\"")!=std::string::npos,"CSG real writer lacks completed manifest");
    check(J::read(root/"output/meta/effective/nozzle.csg")==raw,"CSG real writer did not preserve raw effective asset");
    (void)readVofElectrohydrodynamicsInput((root/"output/meta/effective").string());
    const auto index=J::read(root/"output/meta/frames.csv");check(std::count(index.begin(),index.end(),'\n')==4,"CSG frame index is not actual Initial plus two accepted End frames");
    std::istringstream lines(index);std::string line;std::getline(lines,line);unsigned frames=0;while(std::getline(lines,line)){if(line.empty())continue;std::istringstream values(line);std::vector<std::string> columns;std::string value;while(std::getline(values,value,','))columns.push_back(value);
        check(columns.size()==9&&std::stoul(columns[1])==frames,"CSG frame index accepted steps are not unique and ordered");const auto volume=root/"output"/columns[5];check(volume.extension()==".dat","CSG volume output is not DAT");const auto text=J::read(volume);
        check(text.find("VARIABLES")!=std::string::npos&&text.find("ELEMENTS=2720")!=std::string::npos,"CSG real DAT does not contain exactly selected fluid cells");++frames;}
    check(frames==3&&guards==13,"CSG direct reader/public-run coverage incomplete");
    std::cout<<"VOF_CSG_INITIAL_APP_PASS RZ native=3200 active=2720 solid=480 accepted="<<observedSteps<<" guards="<<guards<<" checks="<<checks<<" original_raw_asset_effective_roundtrip_typed_pin_authoritative_LG_no_solid_Q_actual_runHydrodynamics_sigma0_closed_DAT_only output="<<root.string()<<'\n';
}
}
// END CSG INITIAL PUBLIC RUN CONTRACT


// BEGIN ACTUAL CSG ELECTRIC BOUNDARY-ID CONTRACT
namespace VofCsgElectricBoundaryContract {
namespace fs=std::filesystem;namespace C=VofCsgInitialAppContract;namespace J=VofEhdJointContract;
using Fixture=VofEhdInputContract::Fixture;using VofEhdInputContract::replace;
unsigned checks=0,guards=0,steps=0;long double voltageError=0,gauss=0,maxField=0,maxCurrent=0;
void check(bool ok,const char* why){++checks;if(!ok)throw std::runtime_error(why);}
Fixture fixture(){auto f=C::fixture();
    replace(f.text[6],"boundaryId 2; type Dirichlet; value 0","boundaryId 2; type Dirichlet; value 6");
    replace(f.text[6],"boundaryId 3; type Dirichlet; value 0","boundaryId 3; type Dirichlet; value 6");return f;}
void field(const VofElectrohydrodynamicsTransaction& tx){
    const auto* e=tx.acceptedField();const auto* h=tx.accepted();check(e&&h,"CSG electric actual transaction has no coeval field/fluid");
    const auto& p=h->phaseGeometry().vofPublication();const auto& d=p.domain();const auto& cells=p.plicCells();
    check(e->geometry().identity()==p.identity()&&(e->geometryTime().step==h->timeLayer().step&&e->geometryTime().begin==h->timeLayer().begin&&e->geometryTime().end==h->timeLayer().end&&e->geometryTime().position==h->timeLayer().position),"CSG electric field lost actual publication/clock");
    double lo=6,hi=0;for(std::size_t ci=0;ci<cells.size();++ci){
        if(d.cellRowByIndex()[ci]<0){check(e->potential()[ci]==0&&e->epsilonByCell()[ci]==0&&e->chargeCandidate().integratedCharge()[ci]==0,"excluded electrode acquired phi/material/charge unknown");continue;}
        const double value=e->potential()[ci];check(std::isfinite(value)&&value>=-1e-10&&value<=6+1e-10,"CSG electric boundary potentials violate fixed 0..6 data");lo=std::min(lo,value);hi=std::max(hi,value);}
    check(hi-lo>1,"CSG 300V/code6 electrode and 0 collector collapsed to a constant field");
    unsigned wet=0,dry=0,top=0;
    for(int row=0;row<d.nativeFaceCount();++row){const auto& f=p.faces().at(row);const auto& ef=e->faces().at(row);
        const auto origin=p.frameOrigin();auto a=p.vertices().at(f.vertices[0]),b=p.vertices().at(f.vertices[1]);
        a.x+=origin.x;a.y+=origin.y;b.x+=origin.x;b.y+=origin.y;
        const bool horizontal=a.y==b.y;const long double centerX=(a.x+b.x)/2,centerY=(a.y+b.y)/2;
        if(f.negativeCell<0&&f.positiveCell<0){check(ef.physicalMeasure==0&&ef.integratedField==0&&ef.current==0,"excluded face acquired electric flux");continue;}
        maxField=std::max(maxField,std::abs(static_cast<long double>(ef.integratedField)));maxCurrent=std::max(maxCurrent,std::abs(static_cast<long double>(ef.current)));
        if(f.negativeCell>=0&&f.positiveCell>=0)continue;
        if(ef.physicalMeasure==0){check(ef.integratedField==0&&ef.displacement==0&&ef.current==0,"RZ zero-area axis acquired flux");continue;}
        const bool isWet=!horizontal&&a.x==1&&b.x==1&&centerY<3;
        const bool isDry=horizontal&&a.y==3&&b.y==3&&centerX>1;
        const bool isTop=horizontal&&a.y==16&&b.y==16;
        if(isWet||isDry||isTop){long double coefficient=0;for(const auto& c:ef.gradientCoefficients)coefficient+=c.value;
            check(coefficient!=0,"actual Dirichlet boundary-ID face has no affine boundary coefficient");
            const long double actual=-static_cast<long double>(ef.gradientConstant)/coefficient,expected=isTop?0:6;
            voltageError=std::max(voltageError,std::abs(actual-expected));check(std::abs(actual-expected)<=1e-12L,"actual face affine trace does not encode original voltage by physical ID");
            wet+=isWet;dry+=isDry;top+=isTop;
        }else check(ef.integratedField==0&&ef.current==0,"explicit homogeneous Neumann boundary-ID face acquired flux");
    }
    check(wet>0&&dry>0&&top>0&&maxField>0&&maxCurrent>0,"actual internal wet/dry electrodes, collector or nonzero E/J missing");
    gauss=std::max(gauss,static_cast<long double>(e->report().maximumRelativeGaussResidual));check(gauss<=1e-10L,"CSG actual electric field failed original Gauss gate");
}

void run(const fs::path& requested){
    check(requested.is_absolute()&&requested.root_name().string()=="D:","CSG electric test requires absolute D: output");
    const auto root=requested/("csg_electric_"+std::to_string(HydrodynamicsTransaction::Clock::now().time_since_epoch().count()));fs::create_directories(root);
    const auto source=fixture();const auto system=C::save(root,"system",source);const auto input=readVofElectrohydrodynamicsInput(system.string());
    check(input.csgElectricBoundaries.size()==8,"CSG reader lost exact named boundary-ID table");
    const auto effective=root/"effective";fs::create_directories(effective);for(const auto& doc:input.fluid.documents)VofEhdInputContract::write(effective/fs::path(doc.file).filename(),doc.effectiveText);
    const auto roundtrip=readVofElectrohydrodynamicsInput(effective.string());check(roundtrip.csgElectricBoundaries.size()==8,"CSG electric effective roundtrip lost table");
    for(std::size_t k=0;k<8;++k){const auto& a=input.csgElectricBoundaries[k];const auto& b=roundtrip.csgElectricBoundaries[k];check(a.name==b.name&&a.boundaryId==b.boundaryId&&a.condition.type==b.condition.type&&a.condition.value==b.condition.value,"CSG electric table changed on effective roundtrip");}
    const auto reject=[&](const char* name,auto edit){auto f=fixture();edit(f);const auto path=C::save(root,name,f);bool caught=false;
        try{(void)readVofElectrohydrodynamicsInput(path.string());}catch(const std::exception&){caught=true;}check(caught,"invalid CSG electric patch table accepted");++guards;};
    reject("unknown_id",[](Fixture& f){replace(f.text[6],"boundaryId 2;","boundaryId 999;");});
    reject("duplicate_id",[](Fixture& f){replace(f.text[6],"boundaryId 3;","boundaryId 2;");});
    reject("nonzero_neumann",[](Fixture& f){replace(f.text[6],"boundaryId 5; type Neumann; value 0","boundaryId 5; type Neumann; value 1");});
    VofElectrohydrodynamicsOutput writer(input,{root});VofElectrohydrodynamicsTransaction tx(input,writer);std::string failure;
    const HydrodynamicsRunIO io{[&]{writer.begin();},[&](auto attempted,auto step,double time,double oldDt,double newDt,const std::string& reason){writer.retry(attempted,step,time,oldDt,newDt,reason);},
        [&](auto attempted,auto step,double time,const char* stage,const std::string& reason){failure=std::string(stage)+": "+reason;writer.failure(attempted,step,time,stage,reason);},
        [&]{writer.complete();},[&](const auto& report,const auto& clock,int retries){writer.timeStep(report,clock,retries);},
        [&](const auto& clock,double seconds){check(tx.accepted()!=nullptr,"CSG electric real committed bundle missing");
            field(tx);writer.stepCompleted(clock,seconds);++steps;}};
    RunConfig config;config.app="VofElectrohydrodynamics";const int result=runHydrodynamics(input.fluid,config,tx,io);
    if(result!=0)throw std::runtime_error("CSG electric actual workflow failed: "+failure);
    check(steps==2&&tx.accepted()==nullptr,"CSG electric real workflow did not finish/clear two steps");
    check(J::read(root/"output/meta/manifest.json").find("\"status\":\"COMPLETED\"")!=std::string::npos,"CSG electric workflow did not publish completed manifest");
    std::cout<<std::setprecision(17)<<"VOF_CSG_ELECTRIC_BOUNDARY_PASS accepted="<<steps<<" guards="<<guards<<" checks="<<checks<<" voltage_error="<<voltageError<<" gauss="<<gauss<<" max_field="<<maxField<<" max_current="<<maxCurrent<<" actual_original_CSG_code6_electrode_collector0_equal_epsilon_fixture_no_original_material_or_pressure_claim output="<<root.string()<<'\n';
}
}
// END ACTUAL CSG ELECTRIC BOUNDARY-ID CONTRACT

// BEGIN ACTUAL CSG FLUID BOUNDARY-ID CONTRACT
namespace VofCsgFluidBoundaryContract {
namespace fs=std::filesystem;namespace C=VofCsgInitialAppContract;namespace E=VofCsgElectricBoundaryContract;
namespace J=VofEhdJointContract;namespace H=LAMP::Hydrodynamics;namespace I=LAMP::Interface;namespace O=VofEhdOpenJointContract;
using Fixture=VofEhdInputContract::Fixture;using VofEhdInputContract::replace;
unsigned checks=0,guards=0,steps=0;long double pressureWork=0,maximumPhaseResidual=0;
void report(const VofEhdCsgBoundaryReport*,std::size_t);
void check(bool ok,const char* why){++checks;if(!ok)throw std::runtime_error(why);}
void fieldReceipt(const LAMP::Electric::ElectricVofField& field){const auto& r=field.report();
    check(r.integratedExpandedResidualAvailable&&r.integratedExpandedResidualRows==field.geometry().domain().cellIndexByRow().size(),
        "V56 selected CSG field lacks complete-A actual rows");
    check(r.integratedExpandedResidualScratchBytes<=64u*1024u*1024u&&r.postprocessTimingAvailable&&
        r.integratedExpandedResidualTimingAvailable,"V56 selected CSG receipt/timing admission");
    for(double t:{r.postprocessFaceSeconds,r.postprocessCellSeconds,r.postprocessFinalizeSeconds,r.postprocessSeconds,
        r.integratedExpandedResidualRouteSeconds,r.integratedExpandedResidualSeconds})check(std::isfinite(t)&&t>=0,"V56 CSG actual stage time");
    if(field.isElectrostaticSnapshot())check(field.effectiveDt()==0&&field.chargeCandidate().sameInstance(field.sourceCharge())&&
        field.chargeCandidate().integratedCharge()==field.implicitInputCharge(),"V56 selected snapshot changed same accepted Q");
}
void fieldOutputTail(const fs::path& root){std::istringstream input(J::read(root/"output/series/electric_preparation.dat"));
    std::string line;std::vector<std::string> columns;unsigned rows=0;bool variableHeader=false;
    const std::vector<std::string> tail={"field_postprocess_timing_available","field_postprocess_face_seconds","field_postprocess_cell_seconds",
        "field_postprocess_finalize_seconds","field_postprocess_seconds","solver_expanded_residual_available","solver_expanded_residual_parallel_used",
        "solver_expanded_residual_parallel_fallback","solver_expanded_residual_workers","solver_expanded_residual_rows","solver_expanded_residual_scratch_bytes",
        "solver_expanded_residual_timing_available","solver_expanded_residual_route_seconds","solver_expanded_residual_seconds"};
    while(std::getline(input,line)){
        const auto firstNonSpace=line.find_first_not_of(" \t\r");
        const bool startsVariables=firstNonSpace!=std::string::npos&&line.compare(firstNonSpace,9,"VARIABLES")==0;
        if(startsVariables||(variableHeader&&firstNonSpace!=std::string::npos&&line[firstNonSpace]=='"')){
            variableHeader=true;for(std::size_t at=line.find('"');at!=std::string::npos;){const auto end=line.find('"',at+1);
                check(end!=std::string::npos,"V56 actual VARIABLES quote");columns.push_back(line.substr(at+1,end-at-1));at=line.find('"',end+1);}continue;}
        if(firstNonSpace!=std::string::npos)variableHeader=false;
        std::istringstream row(line);double first=0;if(!(row>>first))continue;std::vector<double> values{first};for(double x;row>>x;)values.push_back(x);
        check(!columns.empty()&&values.size()==columns.size(),"V56 actual numeric row width differs from VARIABLES (stage belongs to ZONE)");
        const auto get=[&](const std::string& key){const auto it=std::find(columns.begin(),columns.end(),key);check(it!=columns.end(),"V56 actual tail column missing");return values[static_cast<std::size_t>(it-columns.begin())];};
        for(const auto& key:tail){const double available=get(key+"Available");check(available==0||available==1,"V56 tail Available is not boolean");
            if(available==1)check(std::isfinite(get(key))&&get(key)>=0,"V56 actual tail numeric payload");}
        check(get("solver_expanded_residual_availableAvailable")==1&&get("solver_expanded_residual_available")==1&&
            get("solver_expanded_residual_rowsAvailable")==1&&get("solver_expanded_residual_rows")>0,
            "V56 actual selected Field residual receipt not forwarded");++rows;
    }
    check(rows>=3,"V56 actual retained Midpoint/physical End field rows absent");
}
Fixture fixture(){auto f=E::fixture();
    replace(f.text[0],"Output {","Output { Series { profile diagnostic; }");
    replace(f.text[0],"endTime .00002;","endTime .00002; Adaptive { minDt .00001; maxDt .00001; targetCourant .4; }");
    f.text[4]="Materials { Phases { Gas { phaseId 3; density 1; dynamicViscosity .1; permittivity 1; conductivity .1; } Liquid { phaseId 2; density 2; dynamicViscosity .1; permittivity 1; conductivity .1; } } Interfaces { Shared { phases (2 3); surfaceTension 0; } } }";
    replace(f.text[5],"Boundary { type impermeable; }",R"(Boundary { type pressureDrivenReservoirs;
        Patches { supply { boundaryId 5; type pressureOutlet; phase liquid; pressure .001; }
            environmentOuter { boundaryId 8; type pressureOutlet; phase gas; pressure 0; }
            collector { boundaryId 9; type pressureOutlet; phase gas; pressure 0; } }
        Backflow { model interiorContinuationZeroShear; } StationaryNoSlip { boundaryIds (2 3); } })");
    const auto closing=f.text[6].find_last_not_of(" \r\n\t");check(closing!=std::string::npos&&f.text[6][closing]=='}',"CSG fluid fixture Electric closing brace absent");f.text[6].erase(closing,1);
    f.text[6]+=" Advection { scheme limitedLinearFct; Inlets { supply { chargeDensity .02; } environmentOuter { chargeDensity -.01; } collector { chargeDensity -.03; } } } }";
    return f;
}
void typed(const VofElectrohydrodynamicsInput& input){const auto& f=input.fluid;
    check(f.csgGeometry&&f.vofPressureDrivenReservoirs&&f.vofBoundaryPatches.size()==3&&f.vofStationaryNoSlipBoundaryIds==std::vector<int>({2,3})
        &&f.vofPressureOutletBackflow&&f.vofPressureOutletBackflow->model==H::HydroVofPressureOutletMomentumModel::InteriorContinuationZeroShear
        &&!f.vofPressureOutletBackflow->velocity&&f.vof.surfaceTension==0&&f.vof.liquidViscosity==.1&&f.vof.gasViscosity==.1,
        "CSG fluid typed pressure/no-slip/positive physical viscosity or sigma0 scope changed");
    const std::array<int,3> ids{{5,8,9}};const std::array<double,3> pressure{{.001,0,0}},rho{{.02,-.01,-.03}};
    for(std::size_t k=0;k<ids.size();++k){const auto& p=f.vofBoundaryPatches[k];check(p.boundaryId&&*p.boundaryId==ids[k]&&p.kind==HydrodynamicsInput::VofBoundaryPatch::Kind::PressureOutlet
        &&p.pressure==pressure[k]&&p.velocity.x==0&&p.velocity.y==0&&input.inletChargeDensityByPatch[k]==rho[k],"CSG physical port pressure or independent signed-Q named reservoir changed");}
    unsigned physical=0,wall=0;for(const auto& e:f.trace.entries){if(e.field.find(".boundaryId")!=std::string::npos&&e.consumer=="Hydrodynamics.VOF.physicalPatchBoundaryId"){check(e.used&&!e.defaulted,"CSG pressure boundary ID trace lost explicit source");++physical;}
        if(e.field=="Hydrodynamics.Boundary.StationaryNoSlip.boundaryIds"){check(e.used&&!e.defaulted,"CSG no-slip ID trace lost explicit source");++wall;}}
    check(physical==3&&wall==1,"CSG physical boundary-ID trace missing or duplicated");
}
void observe(const VofElectrohydrodynamicsTransaction& tx,const LAMP::TimeLayer& clock){const auto* state=tx.accepted();const auto* step=tx.preparedStep();const auto* ef=tx.acceptedField();
    check(state&&step&&ef&&J::sameClock(clock,state->timeLayer()),"CSG actual accepted/prepared fluid/field or clock absent");
    fieldReceipt(*ef);
    const auto& p=state->phaseGeometry().vofPublication();const auto& d=p.domain();const auto& layout=state->values().layout();const auto* flow=state->vofPrescribedFlow();
    check(flow&&flow->pressureDrivenReservoirs&&flow->prescribedVelocityDofs.empty()&&!flow->endpointVelocity&&flow->pressureOutletBackflowModel
        &&*flow->pressureOutletBackflowModel==H::HydroVofPressureOutletMomentumModel::InteriorContinuationZeroShear,
        "CSG actual pressure-driven flow fabricated fixed feed/vector or lost momentum model");
    check(flow->stationaryNoSlipVelocityDofs.size()==47&&std::is_sorted(flow->stationaryNoSlipVelocityDofs.begin(),flow->stationaryNoSlipVelocityDofs.end())
        &&std::adjacent_find(flow->stationaryNoSlipVelocityDofs.begin(),flow->stationaryNoSlipVelocityDofs.end())==flow->stationaryNoSlipVelocityDofs.end(),"CSG actual wet/dry wall DOFs are not complete sorted unique");
    check(std::all_of(flow->normalVelocity.begin(),flow->normalVelocity.end(),[](double u){return u==0;}),"CSG pressure reservoirs carry prescribed normal payload");
    const auto& greenLoad=step->frozenMidpoint().forces().maxwellLoad;check(greenLoad&&greenLoad->geometry().isRt0GreenCompatible()
        &&!greenLoad->geometry().isStressCompatible()&&!greenLoad->geometry().hasCanonicalRawStrain(),"actual no-slip positive-mu Main did not request Green-only Maxwell geometry");
    const auto& viscosity=step->frozenMidpoint().viscosity().viscosity();check(viscosity.strain().stationaryNoSlipVelocityDofs()==flow->stationaryNoSlipVelocityDofs
        &&viscosity.matrix().nonZeros()>0&&!viscosity.rowWeights().empty(),"CSG actual Frozen viscosity did not consume original stationary full-trace no-slip rows");
    const auto& pin=p.pinnedCorner();check(pin&&pin->sourceIdentity()==p.identity()&&J::sameClock(pin->sourceTime(),clock)&&pin->physicalPoint().x==1&&pin->physicalPoint().y==3,"CSG actual open End lost original stock-derived pin certificate");
    const auto& u=state->values().normalVelocity();const auto& matrix=state->phaseGeometry().metrics().divergence().matrix();
    for(int c=0;c<matrix.rows();++c){long double value=0,scale=0;for(H::HydroIntegralDivergence::Matrix::InnerIterator it(matrix,c);it;++it){const long double term=static_cast<long double>(it.value())*u[it.col()];value+=term;scale+=std::abs(term);}O::bound(value,scale,1e-11L,"CSG actual D_I U failed original relative gate");}
    long double supply=0,drain=0,total=0,activity=0,work=0,workScale=0;unsigned wet=0,dry=0,axis=0,ports=0;
    const auto& mid=step->solvedMidpointValues().normalVelocity();const double dt=clock.end-clock.begin;
    for(int row=0;row<d.nativeFaceCount();++row){const auto& face=p.faces().at(row);if((face.negativeCell<0&&face.positiveCell<0)||(face.negativeCell>=0&&face.positiveCell>=0))continue;
        auto a=p.vertices().at(face.vertices[0]),b=p.vertices().at(face.vertices[1]);const auto origin=p.frameOrigin();a.x+=origin.x;a.y+=origin.y;b.x+=origin.x;b.y+=origin.y;
        const bool horizontal=a.y==b.y;const long double x=(a.x+b.x)/2,y=(a.y+b.y)/2;const auto mode=flow->rulesByNativeRow.at(row);const int dof=layout.velocityDof(row);
        const long double measure=ef->faces().at(row).physicalMeasure;const int sign=face.negativeCell<0?-1:1;
        if(measure==0){const auto& axes=state->phaseGeometry().metrics().axisVelocityDofs();
            check(!horizontal&&x==0&&mode==I::GeometricVofAdvectorSlab::Boundary::Mode::Closed&&dof>=0&&u[dof]==0
                &&std::find(axes.begin(),axes.end(),dof)!=axes.end()
                &&!std::binary_search(flow->stationaryNoSlipVelocityDofs.begin(),flow->stationaryNoSlipVelocityDofs.end(),dof),
                "CSG genuine axis lost its original constrained zero-velocity DOF or became a no-slip wall");++axis;continue;}
        check(dof>=0,"CSG positive-area physical face lacks native velocity DOF");const long double flux=sign*measure*u[dof];total+=flux;activity+=std::abs(flux);
        const bool wallWet=!horizontal&&x==1&&y<3,wallDry=horizontal&&y==3&&x>1;
        if(wallWet||wallDry){check(mode==I::GeometricVofAdvectorSlab::Boundary::Mode::Closed&&u[dof]==0
            &&std::binary_search(flow->stationaryNoSlipVelocityDofs.begin(),flow->stationaryNoSlipVelocityDofs.end(),dof),"CSG internal wet/dry wall was inferred as exterior port or lost no-slip");wet+=wallWet;dry+=wallDry;continue;}
        const bool inlet=horizontal&&y==0&&x<1,outer=!horizontal&&x==5&&y>3,top=horizontal&&y==16;
        check(inlet||outer||top,"CSG positive-area physical boundary missing independent original geometry role");
        check(mode==(inlet?I::GeometricVofAdvectorSlab::Boundary::Mode::PureLiquidReservoir:I::GeometricVofAdvectorSlab::Boundary::Mode::PureGasReservoir),"CSG actual native reservoir phase differs from physical ID");
        const auto found=std::find_if(flow->pressureOutlets.begin(),flow->pressureOutlets.end(),[&](const auto& q){return q.velocityDof==dof;});
        check(found!=flow->pressureOutlets.end()&&found->pressure==(inlet?.001:0),"CSG actual physical pressure was not bound to its Core port ID");++ports;(inlet?supply:drain)+=flux;
        const long double term=-dt*found->pressure*sign*measure*mid[dof];work+=term;workScale+=std::abs(term);
    }
    check(wet==15&&dry==32&&axis==80&&ports==113&&flow->pressureOutlets.size()==ports,"CSG actual wet/dry/axis/pressure-port coverage counts changed");
    O::bound(total,activity,1e-12L,"CSG independently summed physical total volume flux");check(activity>0&&supply!=0&&drain!=0&&supply*drain<0,"CSG solved pressure-reservoir flow lacks balanced nonzero physical exchange");
    check(!step->report().stages.empty()&&step->frozenMidpoint().hasMaxwellLoad()&&!step->hasCapillaryLoad(),"CSG sigma0 A transaction missing actual Maxwell or inserted capillary load");
    const auto& final=step->report().stages.back();O::bound(final.pressureOutletWork-work,std::abs(final.pressureOutletWork)+workScale,1e-11L,"CSG physical pressure traction work disagrees with actual port flux");
    check(work!=0&&final.prescribedReactionWork==0,"CSG actual signed pressure work missing or replaced by fixed feed work");pressureWork+=work;std::cout<<std::setprecision(17)<<"CSG_FLUID_SIGNED_WORK step="<<clock.step<<" supply_outward="<<supply<<" environment_outward="<<drain<<" pressure_work="<<work<<std::endl;
    const auto& transport=step->transport();const auto& ledger=transport.report().boundaryFull;check(ledger.enabled&&ledger.physicalRelativeConservation<=1e-12L&&ledger.algebraicRelativeConservation<=1e-12L,"CSG original physical/algebraic boundary phase gates failed");
    const auto source=O::stocks(transport.source()),target=O::stocks(transport.endpoint());std::array<long double,2> outward{},gross{};
    for(const auto& face:transport.faces()){if((face.negativeCell<0&&face.positiveCell<0)||(face.negativeCell>=0&&face.positiveCell>=0))continue;const int sign=face.negativeCell<0?-1:1;
        outward[0]+=sign*face.full.liquid;outward[1]+=sign*face.full.gas;gross[0]+=face.full.positiveLiquid+face.full.negativeLiquid;gross[1]+=face.full.positiveGas+face.full.negativeGas;}
    for(int k=0;k<2;++k){const long double residual=target[k]-source[k]+outward[k];O::bound(residual,std::abs(source[k])+std::abs(target[k])+gross[k],1e-12L,"CSG independent L/G stock plus actual native boundary flux");maximumPhaseResidual=std::max(maximumPhaseResidual,std::abs(residual));}
    check(tx.physicalChargeTransport()&&tx.physicalEndField()&&tx.physicalEndField()->advection(),"CSG actual physical End charge/Ohmic receipts absent");E::field(tx);++steps;
    report(tx.csgBoundaryPreparation(),1);
}
// Use the original complete Core classification as a test-only uncached
// reference. The real phase/seed owner is shared by cold and warm consumers.
struct BindingEnvironment {
    std::fenv_t native{};int error=errno;
#if defined(_M_X64) || defined(__SSE2__)
    unsigned csr=_mm_getcsr();
#endif
    BindingEnvironment(){check(std::fegetenv(&native)==0,"CSG binding FP capture failed");errno=error;}
    void restore() const {check(std::fesetenv(&native)==0,"CSG binding FP restore failed");
#if defined(_M_X64) || defined(__SSE2__)
        _mm_setcsr(csr);
#endif
        errno=error;}
    bool same(const BindingEnvironment& b)const{return std::memcmp(&native,&b.native,sizeof(native))==0&&error==b.error
#if defined(_M_X64) || defined(__SSE2__)
        &&csr==b.csr
#endif
        ;}
    ~BindingEnvironment()noexcept{std::fesetenv(&native);
#if defined(_M_X64) || defined(__SSE2__)
        _mm_setcsr(csr);
#endif
        errno=error;}
};
void report(const VofEhdCsgBoundaryReport* r,std::size_t scope){
    check(r&&r->byteLimit==32u*1024u*1024u,"actual CSG preparation report absent or wrong budget");
    const auto& s=r->stages.at(scope);check(s.requests>0&&s.requests==s.hits+s.rebuilds&&s.selectorCalls==s.rebuilds
        &&s.completeGraphBuilds==s.rebuilds&&s.retainedBuilds<=s.rebuilds&&s.peakRetainedBytes<=r->byteLimit,
        "actual CSG request/rebuild/graph counters or retained budget inconsistent");
    if(s.timingAvailable)check(std::isfinite(s.totalSeconds)&&std::isfinite(s.rebuildSeconds)
        &&s.totalSeconds>=0&&s.rebuildSeconds>=0&&s.rebuildSeconds<=s.totalSeconds,
        "actual CSG nested timing is not finite/bounded");
}
std::string flowPayload(const H::HydroTwoPhaseState& state){
    const auto* f=state.vofPrescribedFlow();check(f,"CSG binding reference lost actual open flow");std::string out;
    J::append(out,state.values().normalVelocity());J::append(out,state.values().pressure());
    J::append(out,f->rulesByNativeRow.size());for(auto x:f->rulesByNativeRow)J::append(out,x);
    J::append(out,f->prescribedVelocityDofs.size());for(int x:f->prescribedVelocityDofs)J::append(out,x);
    J::append(out,f->normalVelocity);J::append(out,f->stationaryNoSlipVelocityDofs.size());for(int x:f->stationaryNoSlipVelocityDofs)J::append(out,x);
    J::append(out,f->pressureOutlets.size());for(const auto& p:f->pressureOutlets){J::append(out,p.velocityDof);J::append(out,p.pressure);}
    J::append(out,f->pressureDrivenReservoirs);J::append(out,bool(f->endpointVelocity));J::append(out,bool(f->pressureOutletBackflowMomentum));
    J::append(out,bool(f->pressureOutletBackflowModel));if(f->pressureOutletBackflowModel)J::append(out,*f->pressureOutletBackflowModel);return out;
}
void bindingContracts(const fs::path& root,const VofElectrohydrodynamicsInput& input){
    const BindingEnvironment caller;const auto& f=input.fluid;
    DictionaryMesh spec;spec.Nx=f.nx;spec.Ny=f.ny;spec.minCoord=f.domainMin;spec.maxCoord=f.domainMax;
    spec.levelSetCount=0;spec.maxInitialIterations=0;spec.maxCellNum=16384;QuadtreeMesh mesh(spec);
    const auto solids=FittedDomainGeometry::selectFixedSolids(mesh,*f.csgGeometry,*f.csgSampler,f.csgExcludedPhases);
    FaceSegmentGraph complete;complete.build(mesh);FaceSegmentGraph selected;
    selected.build(mesh,[&](const QuadtreeMesh&,const Cell& c){return solids.excludedPhaseByCell.at(c.index())==0;});
    const auto domain=NativeLeafDomain::build(mesh,selected,f.radial);
    auto pub=I::GeometricVofInitialization::prepare(mesh,domain,[&](Vector2D p){return f.csgSampler->phaseValue(f.negativePhase,p);},f.vofInitialization);
    check(f.csgPin&&solids.boundaryIdsBySegment.size()==complete.segments().size(),"CSG binding independent Initial source absent");
    const auto& pin=*f.csgPin;const auto point=f.csgGeometry->points.at(pin.pointId-1);std::set<int> wet,dry;
    std::vector<int> ids(static_cast<std::size_t>(domain.nativeFaceCount()));
    for(std::size_t i=0;i<complete.segments().size();++i){const auto& segment=complete.segments()[i];
        const auto& edge=segment.isHorizontal?mesh.hEdges(segment.edgeLabel):mesh.vEdges(segment.edgeLabel);
        const int row=edge.index()+(segment.isHorizontal?mesh.vEdgeIndexNum():0),id=solids.boundaryIdsBySegment.at(i);
        if(id==0)continue;check(ids.at(row)==0||ids.at(row)==id,"uncached Core reference has conflicting nonzero native IDs");ids.at(row)=id;
        if(id!=pin.wetBoundaryId&&id!=pin.dryBoundaryId)continue;
        bool endpoint=false;for(int v:edge.edgeToVert){const auto p=mesh.verts(v).coord();endpoint=endpoint||(p.x==point.x&&p.y==point.y);}
        if(endpoint)(id==pin.wetBoundaryId?wet:dry).insert(row);
    }
    check(wet.size()==1&&dry.size()==1,"CSG binding original exact pin support changed");
    pub=I::GeometricVofPublication::withPinnedCorner(mesh,pub,point,*wet.begin(),*dry.begin(),pin.wallTrace);
    const auto layout=H::HydroDofLayout::build(mesh,domain);
    const auto metric=H::HydroPhysicalDivergence::build(mesh,H::HydroIntegralDivergence::build(mesh,layout));
    const auto phase=H::HydroPhaseGeometry::build(mesh,metric,pub);
    const H::HydroNativeValues seed(layout,std::vector<double>(layout.velocityCount()),std::vector<double>(layout.pressureCount()));
    const auto probe=root/"binding_probe";VofElectrohydrodynamicsOutput writer(input,{probe}),freshWriter(input,{root/"binding_fresh"});
    VofElectrohydrodynamicsTransaction tx(input,writer),fresh(input,freshWriter);
    check(std::fesetround(FE_TONEAREST)==0&&std::feraiseexcept(FE_ALL_EXCEPT)==0,"CSG binding cannot establish finite sticky-state tape");
#if defined(_M_X64) || defined(__SSE2__)
    _mm_setcsr((_mm_getcsr()&~(0x6000u|0x8040u))|0x1f80u|0x3fu);
#endif
    errno=173;const BindingEnvironment entry;
    auto cold=tx.makeInitialVofState(mesh,phase,seed);const BindingEnvironment coldFp;
    report(tx.csgBoundaryPreparation(),0);
#if defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
    check(tx.csgBoundaryPreparation()->stages[0].retainedBuilds==1,"CSG binding tape did not retain a quiet actual cold classification");
#endif
    entry.restore();auto warm=tx.makeInitialVofState(mesh,phase,seed);const BindingEnvironment warmFp;
    entry.restore();auto reference=fresh.makeInitialVofState(mesh,phase,seed);const BindingEnvironment referenceFp;
    check(flowPayload(cold)==flowPayload(warm)&&flowPayload(warm)==flowPayload(reference)
        &&coldFp.same(warmFp)&&warmFp.same(referenceFp),"CSG cached/fresh original flow/U/p or complete native FP/rawCSR/errno differ");
    check(warm.phaseGeometry().vofIdentity()==pub.identity()&&reference.phaseGeometry().vofIdentity()==pub.identity()
        &&J::sameClock(warm.timeLayer(),pub.timeLayer()),"CSG flow cache changed genuine phase owner or Initial clock");
#if defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
    check(tx.csgBoundaryPreparation()->stages[0].hits==1&&fresh.csgBoundaryPreparation()->stages[0].hits==0,
        "CSG warm/fresh comparison did not execute a real hit and uncached classification");
#endif
    const auto* flow=warm.vofPrescribedFlow();unsigned walls=0,ports=0;
    for(const auto face:domain.faces(mesh)){if((face.negativeCell<0)==(face.positiveCell<0)||face.area==0)continue;
        const int row=face.nativeRow,u=layout.velocityDof(row),id=ids.at(row);check(id>0&&u>=0,"uncached Core reference lacks physical boundary ID");
        if(id==2||id==3){check(flow->rulesByNativeRow.at(row)==I::GeometricVofAdvectorSlab::Boundary::Mode::Closed
            &&std::binary_search(flow->stationaryNoSlipVelocityDofs.begin(),flow->stationaryNoSlipVelocityDofs.end(),u),"cached wall disagrees with original uncached Core ID");++walls;}
        else{check(id==5||id==8||id==9,"uncached Core selected boundary has unexpected physical role");
            check(flow->rulesByNativeRow.at(row)==(id==5?I::GeometricVofAdvectorSlab::Boundary::Mode::PureLiquidReservoir:I::GeometricVofAdvectorSlab::Boundary::Mode::PureGasReservoir),"cached phase rule differs from uncached Core ID");
            const auto p=std::find_if(flow->pressureOutlets.begin(),flow->pressureOutlets.end(),[&](const auto& p){return p.velocityDof==u;});
            check(p!=flow->pressureOutlets.end()&&p->pressure==(id==5?.001:0),"cached pressure differs from uncached physical table");++ports;}
    }
    check(walls==47&&ports==113,"CSG uncached physical coverage changed");caller.restore();writer.begin();const char* stage="CSG cache probe";
    tx.initializeVof(mesh,std::move(warm),HydrodynamicsTransaction::Clock::now(),stage);const J::Snapshot accepted(mesh,tx,probe);
    const LAMP::TimeLayer target{1,0,f.run.dt,LAMP::TimeLayer::Position::End};
    for(int attempt=0;attempt<2;++attempt){tx.prepare(mesh,target,f.step,false,stage);accepted.unchanged(mesh,tx,probe);
        const auto* prepared=tx.preparedStep();check(prepared&&prepared->frozenMidpoint().hasMaxwellLoad(),"CSG actual release/retry omitted Frozen/Maxwell");
        const auto& frozen=prepared->frozenMidpoint();const auto& green=frozen.forces().maxwellLoad;const auto& actual=frozen.viscosity().viscosity();
        check(green&&green->geometry().isRt0GreenCompatible()&&!green->geometry().isStressCompatible()&&!green->geometry().hasCanonicalRawStrain()
            &&frozen.viscosity().liquidViscosity()==.1&&frozen.viscosity().gasViscosity()==.1,"CSG positive-mu retry lost actual Green/no-slip request");
        const auto raw=H::HydroStrainOperator::buildStationaryNoSlip(mesh,actual.strain().metrics(),actual.strain().points(),flow->stationaryNoSlipVelocityDofs);
        const auto independent=H::HydroEnergyViscosity::build(raw,actual.viscosity());
        const auto& ak=actual.matrix();const auto& bk=independent.matrix();check(ak.rows()==bk.rows()&&ak.cols()==bk.cols()&&ak.nonZeros()==bk.nonZeros(),"Green-only Maxwell changed independent no-slip K shape");
        for(Eigen::Index row=0;row<=ak.outerSize();++row)check(ak.outerIndexPtr()[row]==bk.outerIndexPtr()[row],"Green-only Maxwell changed independent no-slip K row offsets");
        for(Eigen::Index k=0;k<ak.nonZeros();++k)check(ak.innerIndexPtr()[k]==bk.innerIndexPtr()[k]&&std::memcmp(ak.valuePtr()+k,bk.valuePtr()+k,sizeof(double))==0,"Green-only Maxwell changed independent no-slip K coefficient/sign bits");
        std::cout<<"RT0_GREEN_MAIN_NOSLIP_K_PASS attempt="<<attempt<<" original_positive_mu=.1 independent_public_noSlip_K=1\n";
        report(tx.csgBoundaryPreparation(),1);check(tx.csgBoundaryPreparation()->stages[0].requests==0
            &&(attempt==0||tx.csgBoundaryPreparation()->stages[1].rebuilds>=1),"CSG retry reused a consumed/discarded attempt or accumulated Initial report");
        tx.releaseStep();check(!tx.csgBoundaryPreparation()&&!tx.preparedStep(),"CSG release retained attempt report/candidate");accepted.unchanged(mesh,tx,probe);
        bool discarded=false;try{tx.commit(mesh);}catch(const std::exception&){discarded=true;}
        check(discarded,"CSG discarded candidate remained committable");accepted.unchanged(mesh,tx,probe);}
    // ID 6 is declared and nonignored, but lies entirely on the excluded solid.
    // A warm ID map must still run the current electric table's coverage gate.
    auto invalidFixture=fixture();replace(invalidFixture.text[6],"Patches {","Patches { unusedSolid { boundaryId 6; type Neumann; value 0; }");
    const auto invalid=readVofElectrohydrodynamicsInput(C::save(root,"uncovered_actual_id",invalidFixture).string());
    VofElectrohydrodynamicsOutput badWriter(invalid,{root/"binding_uncovered"});VofElectrohydrodynamicsTransaction bad(invalid,badWriter);
    auto badState=bad.makeInitialVofState(mesh,phase,seed);badWriter.begin();const auto before=J::read(root/"binding_uncovered/output/meta/frames.csv");
    bool rejected=false;try{bad.initializeVof(mesh,std::move(badState),HydrodynamicsTransaction::Clock::now(),stage);}
    catch(const std::exception& e){rejected=std::string(e.what()).find("no actual selected-fluid coverage")!=std::string::npos;}
    check(rejected&&!bad.accepted()&&J::read(root/"binding_uncovered/output/meta/frames.csv")==before,
        "warm CSG map hid a declared-but-uncovered electric condition or published failed Initial state");
    std::cout<<"CSG_BINDING_CACHE_PASS actual_hit_vs_fresh_original_Core_map complete_flow_U_p_FP_errno release_retry_2 uncovered_declared_id_6\n";
}
void run(const fs::path& requested){check(requested.is_absolute()&&requested.root_name().string()=="D:","CSG fluid selector requires absolute D: root");
    const auto root=requested/("csg_fluid_"+std::to_string(HydrodynamicsTransaction::Clock::now().time_since_epoch().count()));fs::create_directories(root);
    const auto system=C::save(root,"system",fixture());const auto input=readVofElectrohydrodynamicsInput(system.string());typed(input);
    const auto effective=root/"effective";fs::create_directories(effective);for(const auto& doc:input.fluid.documents)VofEhdInputContract::write(effective/fs::path(doc.file).filename(),doc.effectiveText);
    const auto again=readVofElectrohydrodynamicsInput(effective.string());typed(again);
    const auto reject=[&](const char* name,auto edit){auto f=fixture();edit(f);const auto path=C::save(root,name,f);bool caught=false;try{(void)readVofElectrohydrodynamicsInput(path.string());}catch(const std::exception&){caught=true;}check(caught,"invalid CSG fluid physical-ID input accepted");++guards;};
    reject("id_and_side",[](Fixture& f){replace(f.text[5],"boundaryId 5;","boundaryId 5; side S; interval (0 1);");});
    reject("missing_id",[](Fixture& f){replace(f.text[5],"boundaryId 5;","side S; interval (0 1);");});
    reject("duplicate_port",[](Fixture& f){replace(f.text[5],"boundaryId 8;","boundaryId 5;");});
    reject("unknown_port",[](Fixture& f){replace(f.text[5],"boundaryId 5;","boundaryId 999;");});
    reject("wall_port_overlap",[](Fixture& f){replace(f.text[5],"boundaryIds (2 3)","boundaryIds (2 3 5)");});
    reject("duplicate_wall",[](Fixture& f){replace(f.text[5],"boundaryIds (2 3)","boundaryIds (2 3 3)");});
    reject("missing_walls",[](Fixture& f){replace(f.text[5],"StationaryNoSlip { boundaryIds (2 3); }","");});
    bindingContracts(root,input);
    VofElectrohydrodynamicsOutput writer(input,{root});VofElectrohydrodynamicsTransaction tx(input,writer);std::string failure;
    const HydrodynamicsRunIO io{[&]{writer.begin();},[&](auto attempted,auto step,double time,double oldDt,double newDt,const std::string& reason){writer.retry(attempted,step,time,oldDt,newDt,reason);},
        [&](auto attempted,auto step,double time,const char* stage,const std::string& reason){failure=std::string(stage)+": "+reason;writer.failure(attempted,step,time,stage,reason);},
        [&]{writer.complete();},[&](const auto& report,const auto& clock,int retries){observe(tx,clock);writer.timeStep(report,clock,retries);},
        [&](const auto& clock,double seconds){writer.stepCompleted(clock,seconds);}};
    RunConfig config;config.app="VofElectrohydrodynamics";const int result=runHydrodynamics(input.fluid,config,tx,io);
    if(result!=0)throw std::runtime_error("CSG fluid actual public workflow failed: "+failure);
    check(steps==2&&!tx.accepted()&&pressureWork!=0,"CSG fluid actual two-step pressure/no-slip workflow incomplete");
    check(J::read(root/"output/meta/manifest.json").find("\"status\":\"COMPLETED\"")!=std::string::npos,"CSG fluid completed manifest missing");
    const auto index=J::read(root/"output/meta/frames.csv");check(std::count(index.begin(),index.end(),'\n')==4,"CSG fluid real DAT index lacks Initial plus two End frames");
    const auto bindingTable=J::read(root/"output/series/csg_boundary_preparation.dat");
    check(bindingTable.find("physical_geometry_identity")!=std::string::npos&&bindingTable.find("accepted_geometry_identity")!=std::string::npos
        &&bindingTable.find("complete_graph_builds")!=std::string::npos,"CSG DAT lost declared identity/scope/count schema");
    std::istringstream lines(bindingTable);std::string line;unsigned bindingRows=0;
    while(std::getline(lines,line)){std::istringstream row(line);double first=0;if(!(row>>first))continue;
        std::vector<double> values{first};for(double x;row>>x;)values.push_back(x);
        check(values.size()==32,"CSG DAT does not have 16 logical numeric columns and 16 availability columns");
        check(values[2]==values[3]&&values[4]==(bindingRows?1:0)&&values[5]>0&&values[5]==values[6]+values[7]
            &&values[8]<=values[7]&&values[9]==values[7]&&values[10]==values[7]&&values[11]<=values[12]
            &&values[12]==32u*1024u*1024u,"CSG DAT provenance/scopes/counters disagree with actual fixed-domain workflow");++bindingRows;
    }
    check(bindingRows==3,"CSG DAT repeated an Initial/physical scope or omitted an actual Initial/End output");
    fieldOutputTail(root);
    std::cout<<std::setprecision(17)<<"VOF_CSG_FLUID_BOUNDARY_PASS accepted="<<steps<<" guards="<<guards<<" checks="<<checks<<" pressure_work="<<pressureWork<<" phase_residual="<<maximumPhaseResidual
        <<" actual_original_CSG_fixed_uniform_sigma0_wet_dry_fulltrace_noSlip_pressure_ports_signed_Q_E_J_no_original_material_or_positiveSigma_claim output="<<root.string()<<'\n';
}
}
// END ACTUAL CSG FLUID BOUNDARY-ID CONTRACT


// BEGIN CSG SELECTED TARGET AMR PUBLIC RUN CONTRACT
#include "Core/Mesh/NativeLeafOverlap.h"
namespace VofCsgTargetAmrContract {
namespace fs=std::filesystem;namespace H=LAMP::Hydrodynamics;namespace I=LAMP::Interface;
namespace C=VofCsgInitialAppContract;namespace J=VofEhdJointContract;using VofEhdInputContract::replace;
unsigned checks=0,steps=0,changed=0,released=0;long double coverage=0;
std::array<long double,2> original{};std::vector<std::size_t> outputCells{2720};std::uint64_t previousQ=0;
void check(bool ok,const char* why){++checks;if(!ok)throw std::runtime_error(why);}
std::array<long double,2> totals(const I::GeometricVofPublication& p){std::array<long double,2> result{};
    for(int ci:p.domain().cellIndexByRow())for(int k=0;k<2;++k)result[k]+=p.plicCells()[ci].inventory[k];return result;}
void observe(const VofElectrohydrodynamicsTransaction& tx,const LAMP::TimeLayer& clock,const VofElectrohydrodynamicsInput& input){
    const auto* h=tx.accepted();const auto* q=tx.acceptedCharge();const auto* e=tx.acceptedField();const auto* step=tx.preparedStep();const auto* receipt=tx.remapReceipt();
    check(h&&q&&e&&step&&receipt&&receipt->decision&&receipt->sourceMesh&&receipt->physicalEnd&&receipt->physicalStep==step,
        "CSG public AMR observation must occur after commit and before release of actual physical/tree receipts");
    VofCsgFluidBoundaryContract::fieldReceipt(*e);
    const auto& green=step->frozenMidpoint().forces().maxwellLoad;check(green&&green->geometry().isRt0GreenCompatible()
        &&!green->geometry().isStressCompatible()&&!green->geometry().hasCanonicalRawStrain(),"CSG changed/noop AMR actual A did not retain Green-only request");
    check(clock.step==static_cast<std::int64_t>(steps+1)&&clock.position==LAMP::TimeLayer::Position::End,"CSG public AMR observations skipped an actual End");
    receipt->physicalEnd->validateFor(*receipt->sourceMesh);
    const auto& begin=step->beginState().phaseGeometry().vofPublication();const auto& physical=receipt->physicalEnd->phaseGeometry().vofPublication();
    const auto* physicalQ=tx.physicalChargeTransport();const auto* physicalField=tx.physicalEndField();
    check(physicalQ&&physicalField&&J::sameClock(physicalField->candidateTime(),clock)&&!physicalField->isElectrostaticSnapshot()
        &&physicalField->effectiveDt()==clock.end-clock.begin,"CSG AMR missing one real physical End transport/Ohmic receipt");
    physicalQ->sourceCharge().validateFor(*receipt->sourceMesh,begin,begin.timeLayer());
    check(physicalQ->transport().sameInstance(step->transport())&&physicalQ->transport().sourceIdentity()==begin.identity()
        &&physicalQ->transport().endpoint().identity()==physical.identity(),"CSG physical Q did not consume the same actual Interface End transport");
    if(steps==0){check(begin.domain().cellRowByIndex().size()==3200&&begin.domain().cellIndexByRow().size()==2720&&begin.timeLayer().position==LAMP::TimeLayer::Position::Initial,
            "CSG public Run did not initialize original 3200-native/2720-fluid fitted nozzle");original=totals(begin);}
    else check(physicalQ->sourceCharge().identity()==previousQ,"CSG physical Q advance restarted from a different accepted Q owner");
    const auto& p=h->phaseGeometry().vofPublication();const auto& d=p.domain();const auto& layout=h->values().layout();const auto n=d.cellRowByIndex().size();
    check(J::sameClock(h->timeLayer(),clock)&&J::sameClock(p.timeLayer(),clock)&&J::sameClock(q->timeLayer(),clock)&&J::sameClock(e->candidateTime(),clock)
        &&e->geometry().identity()==p.identity()&&e->chargeCandidate().sameInstance(*q),"CSG target AMR lost same-End actual fluid/Q/field owners");
    check(d.axisymmetric()&&d.cellIndexByRow().size()<n&&layout.domain().cellRowByIndex()==d.cellRowByIndex()
        &&layout.domain().cellIndexByRow()==d.cellIndexByRow()&&p.plicCells().size()==n&&q->integratedCharge().size()==n&&e->potential().size()==n,
        "CSG AMR target became complete domain or lost exact selected native maps/slots");
    unsigned inactive=0,fine=0,coarse=0;long double capacity=0;const double baseWidth=(input.fluid.domainMax.x-input.fluid.domainMin.x)/input.fluid.nx;
    for(std::size_t ci=0;ci<n;++ci){const int row=d.cellRowByIndex()[ci];const auto& cell=p.plicCells()[ci];
        check(layout.pressureDof(static_cast<int>(ci))==row,"CSG target pressure map differs from actual selected geometry");
        if(row<0){++inactive;check(cell.inventory[0]==0&&cell.inventory[1]==0&&cell.measured[0]==0&&cell.measured[1]==0&&!cell.resolved
            &&p.cells()[ci].capacity==0&&q->integratedCharge()[ci]==0&&e->potential()[ci]==0&&e->epsilonByCell()[ci]==0,
            "CSG target excluded solid acquired fake phase/Q/field material");continue;}
        check(d.cellIndexByRow().at(row)==static_cast<int>(ci)&&layout.nativeCell(row)==static_cast<int>(ci),"CSG target selected native maps are not invertible");
        const auto center=(cell.rectangle.lower+cell.rectangle.upper)*.5;check(!(center.x>1&&center.y<3),"CSG target active fluid overlaps original excluded electrode interior");
        const double width=cell.rectangle.upper.x-cell.rectangle.lower.x;fine+=width==baseWidth/2;coarse+=width==baseWidth;capacity+=d.cellVolumes()[row];
        check(q->integratedCharge()[ci]==0&&e->potential()[ci]==0,"CSG zero charge/voltage target acquired Q or phi");
    }
    check(inactive>0&&fine>0&&coarse>0&&n>3200&&n<12800,"CSG AMR failed to refine locally while retaining coarse selected fluid and excluded solids");
    VofEhdOpenJointContract::bound(capacity-328*std::acos(-1.L),328*std::acos(-1.L),1e-12L,"CSG target selected physical-domain RZ coverage");
    const auto after=totals(p);for(int k=0;k<2;++k)VofEhdOpenJointContract::bound(after[k]-original[k],original[k],1e-12L,"CSG same-End target changed independent global phase stock");
    check(std::all_of(h->values().normalVelocity().begin(),h->values().normalVelocity().end(),[](double value){return value==0;}),"CSG quiet AMR created mechanical flow");
    const auto pin=p.pinnedCorner();check(pin&&pin->sourceIdentity()==p.identity()&&J::sameClock(pin->sourceTime(),clock)&&pin->physicalPoint().x==1&&pin->physicalPoint().y==3,
        "CSG target pin was not requalified at actual owner/full clock");
    for(const auto& face:p.faces())if(face.negativeCell>=0&&face.positiveCell>=0){
        const auto& a=p.plicCells()[face.negativeCell].rectangle;const auto& b=p.plicCells()[face.positiveCell].rectangle;
        const double ha=a.upper.x-a.lower.x,hb=b.upper.x-b.lower.x;check(ha<=2*hb&&hb<=2*ha,"CSG target actual native-face adjacency violates quadtree 2:1");}
    check(receipt->decision->additionalRefineCells==0&&receipt->decision->additionalHoldCells==0,"zero-charge CSG AMR invented charge refinement requirements");
    if(receipt->decision->changed){++changed;check(receipt->interfaceTransfer&&receipt->fluidTransfer,"changed CSG AMR omitted separate same-End Interface/Hydro transfers");
        const auto& transfer=*receipt->interfaceTransfer;const auto& overlap=transfer.nativeOverlaps();coverage=std::max(coverage,overlap.maximumRelativeCoverageResidual());
        check(transfer.sourceIdentity()==physical.identity()&&transfer.sourceIdentity()!=begin.identity()&&transfer.candidate().identity()==p.identity(),
            "CSG AMR transfer source/target is not actual physical End/accepted target publication");
        check(e->isElectrostaticSnapshot()&&e->effectiveDt()==0&&!physicalField->chargeCandidate().sameInstance(*q),"CSG AMR target repeated Ohmic time instead of same-End Q remap/snapshot");
        const auto& sourceMap=physical.domain().cellRowByIndex();for(std::size_t ci=0;ci<sourceMap.size();++ci)
            check((sourceMap[ci]>=0)==!overlap.bySource()[ci].empty(),"CSG inactive source acquired physical transfer overlap");
        for(std::size_t ci=0;ci<n;++ci)check((d.cellRowByIndex()[ci]>=0)==!overlap.byTarget()[ci].empty(),"CSG inactive target acquired physical transfer overlap");
        const auto& report=transfer.report();for(int k=0;k<2;++k)check(report.relativeGlobalPhaseResidual[k]<=1e-12L,"CSG independent same-End phase remap failed original relative gate");
    }else check(!receipt->interfaceTransfer&&!receipt->fluidTransfer&&!e->isElectrostaticSnapshot()&&e->chargeCandidate().sameInstance(physicalField->chargeCandidate()),
        "unchanged CSG AMR fabricated physical remap or repeated Q advancement");
    VofCsgFluidBoundaryContract::report(tx.csgBoundaryPreparation(),1);
    if(receipt->decision->changed){VofCsgFluidBoundaryContract::report(tx.csgBoundaryPreparation(),2);
        check(tx.csgBoundaryPreparation()->stages[2].rebuilds>0,"CSG target owner change did not rebuild actual selected boundary IDs");}
    previousQ=q->identity();outputCells.push_back(d.cellIndexByRow().size());++steps;
}
void run(const fs::path& requested){
    check(requested.is_absolute()&&requested.root_name().string()=="D:","CSG target-AMR selector requires absolute D: output root");
    const auto root=requested/("csg_target_amr_"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));check(!fs::exists(root),"CSG AMR output must be fresh");fs::create_directories(root);
    auto f=C::fixture();f.text[0]="Run { app VofElectrohydrodynamics; } Parallel { threads 1; } Time { units code; dt .00001; endTime .00002; Adaptive { minDt .00001; maxDt .00001; targetCourant .2; safety .9; maxGrowth 1; } } Output { Fields { enabled true; writeEverySteps 1; } Series { profile diagnostic; } }";
    f.text[2]="Mesh { resolution (40 80); initialMesh uniform; Regrid { minLevel 0; maxLevel 1; interfaceBandLayers 1; vorticityTolerance 1000; coarsenRatio .25; } }";
    replace(f.text[6],"Electric {","Electric { Regrid { chargeDensityJumpTolerance .01; }");
    const auto input=readVofElectrohydrodynamicsInput(C::save(root,"system",f).string());
    check(input.fluid.csgGeometry&&input.fluid.csgSampler&&input.fluid.vofRegrid&&input.fluid.run.adaptive&&input.chargeDensityJumpTolerance==.01,
        "CSG actual reader did not preserve explicit AMR/charge/adaptive policy");
    const auto effective=root/"effective";fs::create_directories(effective);for(const auto& doc:input.fluid.documents)VofEhdInputContract::write(effective/fs::path(doc.file).filename(),doc.effectiveText);
    const auto again=readVofElectrohydrodynamicsInput(effective.string());check(again.fluid.csgGeometry&&again.fluid.vofRegrid&&again.fluid.run.adaptive&&again.chargeDensityJumpTolerance==.01,
        "CSG AMR effective configuration lost immutable asset or actual consumed policies");
    VofElectrohydrodynamicsOutput writer(input,{root});VofElectrohydrodynamicsTransaction tx(input,writer);std::string failure;unsigned failures=0;
    const HydrodynamicsRunIO io{[&]{writer.begin();},[&](auto attempted,auto step,double time,double oldDt,double newDt,const std::string& reason){writer.retry(attempted,step,time,oldDt,newDt,reason);},
        [&](auto attempted,auto step,double time,const char* stage,const std::string& reason){++failures;failure=std::string(stage)+": "+reason;writer.failure(attempted,step,time,stage,reason);},
        [&]{check(steps==2&&released==2,"CSG public workflow skipped pre-release observations or post-release callbacks");writer.complete();},
        [&](const auto& report,const auto& clock,int retries){check(retries==0,"CSG quiet AMR unexpectedly retried");observe(tx,clock,input);writer.timeStep(report,clock,retries);},
        [&](const auto& clock,double seconds){check(tx.accepted()&&J::sameClock(tx.accepted()->timeLayer(),clock)&&!tx.preparedStep()&&!tx.remapReceipt()&&!tx.physicalChargeTransport()&&!tx.physicalEndField(),
            "CSG completed-step callback did not follow release of physical/tree receipts");++released;writer.stepCompleted(clock,seconds);}};
    RunConfig config;config.app="VofElectrohydrodynamics";const int result=runHydrodynamics(input.fluid,config,tx,io);
    if(result!=0)throw std::runtime_error("CSG target-AMR actual public workflow failed: "+failure);
    check(failures==0&&steps==2&&changed>=1&&!tx.accepted()&&!tx.acceptedCharge()&&!tx.acceptedField(),"CSG public two-step workflow did not complete/clear changed selected AMR");
    check(J::read(root/"output/meta/manifest.json").find("\"status\":\"COMPLETED\"")!=std::string::npos,"CSG AMR real writer lacks completed transaction manifest");
    std::istringstream lines(J::read(root/"output/meta/frames.csv"));std::string line;std::getline(lines,line);unsigned frames=0;
    while(std::getline(lines,line)){if(line.empty())continue;std::istringstream values(line);std::vector<std::string> columns;std::string value;while(std::getline(values,value,','))columns.push_back(value);
        check(columns.size()==9&&frames<outputCells.size()&&std::stoul(columns[1])==frames,"CSG AMR DAT index lost unique actual Initial/End ordering");
        const auto volume=root/"output"/columns[5];check(volume.extension()==".dat"&&J::read(volume).find("ELEMENTS="+std::to_string(outputCells[frames]))!=std::string::npos,"CSG actual DAT contains a different selected cell count");++frames;}
    check(frames==3,"CSG public AMR real DAT writer did not publish Initial plus both End frames");
    VofCsgFluidBoundaryContract::fieldOutputTail(root);
    std::cout<<"VOF_CSG_TARGET_AMR_PASS accepted="<<steps<<" changed="<<changed<<" checks="<<checks<<" coverage="<<static_cast<double>(coverage)
        <<" actual_runHydrodynamics_original_CSG_selected_domain_typed_pin_independent_LG_zero_Q_requirements_once_End_Ohmic_target_snapshot_DAT output="<<root.string()<<'\n';
}
}
// END CSG SELECTED TARGET AMR PUBLIC RUN CONTRACT


// BEGIN ELECTROSPRAY EXPLICIT READER PROFILE CONTRACT
namespace ElectrosprayProfileInputContract {
namespace fs=std::filesystem;namespace C=VofCsgInitialAppContract;namespace F=VofCsgFluidBoundaryContract;namespace J=VofEhdJointContract;
using Fixture=VofEhdInputContract::Fixture;using Application=VofElectrohydrodynamicsInput::Application;using VofEhdInputContract::replace;
unsigned checks=0,readers=0,guards=0;
void check(bool ok,const char* why){++checks;if(!ok)throw std::runtime_error(why);}
Fixture fixture(){auto f=F::fixture();replace(f.text[0],"app VofElectrohydrodynamics;","app Electrospray;");return f;}
void provenance(const VofElectrohydrodynamicsInput& input,const char* expected){
    check(std::string(vofElectrohydrodynamicsApplicationName(input.application))==expected,"typed application profile lost its canonical name");
    unsigned traces=0,runDocuments=0,assets=0;
    for(const auto& e:input.fluid.trace.entries)if(e.field=="Run.app"){
        const auto* name=std::get_if<std::string>(&e.value);check(e.used&&!e.defaulted&&e.consumer=="Runtime.App.dispatch"&&name&&*name==expected&&fs::path(e.file).filename()=="run.dict",
            "actual Run.app trace was replaced by VOF-EHD engine name or lost explicit source");++traces;}
    for(const auto& document:input.fluid.documents){const auto name=fs::path(document.file).filename();
        if(name=="run.dict"){++runDocuments;const auto marker=std::string("app ")+expected+";";
            check(document.sourceText.find(marker)!=std::string::npos&&document.effectiveText.find(marker)!=std::string::npos,"application source/effective Run.app no longer matches explicit profile");}
        if(name=="nozzle.csg"){++assets;check(document.sourceText==J::read(C::originalCsg())&&document.effectiveText==document.sourceText,"profile lost declared original Core CSG raw asset");}}
    check(traces==1&&runDocuments==1&&assets==1,"profile source/effective/trace has missing or duplicate App/raw-asset provenance");
}
void samePhysics(const VofElectrohydrodynamicsInput& a,const VofElectrohydrodynamicsInput& b){
    const auto& x=a.fluid;const auto& y=b.fluid;
    check(a.method==b.method&&a.initialChargeType==b.initialChargeType&&a.materials.epsilonLiquid==b.materials.epsilonLiquid&&a.materials.epsilonGas==b.materials.epsilonGas
        &&a.materials.conductivityLiquid==b.materials.conductivityLiquid&&a.materials.conductivityGas==b.materials.conductivityGas
        &&a.controls.chargeAdvection==b.controls.chargeAdvection&&a.inletChargeDensityByPatch==b.inletChargeDensityByPatch&&a.csgElectricBoundaries.size()==b.csgElectricBoundaries.size(),
        "Electrospray profile changed declared physical electric parameters or charge transport");
    check(x.csgGeometry&&y.csgGeometry&&x.csgSampler&&y.csgSampler&&x.csgPin&&y.csgPin&&x.csgPin->pointId==y.csgPin->pointId
        &&x.csgPin->wetBoundaryId==y.csgPin->wetBoundaryId&&x.csgPin->dryBoundaryId==y.csgPin->dryBoundaryId&&x.csgPin->wallTrace==y.csgPin->wallTrace&&x.csgExcludedPhases==y.csgExcludedPhases
        &&x.vofPressureDrivenReservoirs&&y.vofPressureDrivenReservoirs&&x.vofStationaryNoSlipBoundaryIds==y.vofStationaryNoSlipBoundaryIds
        &&x.vofBoundaryPatches.size()==y.vofBoundaryPatches.size()&&x.vof.surfaceTension==y.vof.surfaceTension&&x.vof.liquidViscosity==y.vof.liquidViscosity&&x.vof.gasViscosity==y.vof.gasViscosity
        &&x.vof.liquidDensity==y.vof.liquidDensity&&x.vof.gasDensity==y.vof.gasDensity&&x.nx==y.nx&&x.ny==y.ny&&x.run.dt==y.run.dt&&x.run.endTime==y.run.endTime,
        "Electrospray profile changed CSG/pin/pressure/no-slip/Hydro physical configuration");
    for(std::size_t k=0;k<x.vofBoundaryPatches.size();++k){const auto& p=x.vofBoundaryPatches[k];const auto& q=y.vofBoundaryPatches[k];
        check(p.name==q.name&&p.boundaryId==q.boundaryId&&p.kind==q.kind&&p.reservoir==q.reservoir&&p.pressure==q.pressure&&p.velocity.x==q.velocity.x&&p.velocity.y==q.velocity.y,
            "profile inserted different physical reservoir or prescribed feed values");}
    for(std::size_t k=0;k<a.csgElectricBoundaries.size();++k){const auto& p=a.csgElectricBoundaries[k];const auto& q=b.csgElectricBoundaries[k];
        check(p.name==q.name&&p.boundaryId==q.boundaryId&&p.condition.type==q.condition.type&&p.condition.value==q.condition.value,"profile inserted different electrode voltage/boundary parameters");}
}
using WallTraceLaw=LAMP::Interface::GeometricVofPublication::PinCertificate::WallTraceLaw;
unsigned wallTraceChecks=0;
void wallTraceCheck(bool ok,const char* why){++wallTraceChecks;if(!ok)throw std::runtime_error(why);}
const char* wallTraceName(WallTraceLaw law){
    switch(law){case WallTraceLaw::StrictWholeCell:return "strictWholeCell";
        case WallTraceLaw::PrescribedWetDry:return "prescribedWetDry";}
    throw std::invalid_argument("unknown pin wall trace law in reader contract");
}
void wallTraceEvidence(const VofElectrohydrodynamicsInput& input,WallTraceLaw expected,bool defaulted){
    const auto word=std::string(wallTraceName(expected));
    wallTraceCheck(input.fluid.csgPin&&input.fluid.csgPin->wallTrace==expected,"typed CSG pin wallTrace differs from its source declaration");
    unsigned traces=0,documents=0;
    for(const auto& entry:input.fluid.trace.entries)if(entry.field=="Geometry.Initial.Pin.wallTrace"){
        const auto* value=std::get_if<std::string>(&entry.value);
        wallTraceCheck(value&&*value==word&&entry.used&&entry.defaulted==defaulted&&entry.consumer=="Interface.PinnedCorner.wallTrace"&&
            fs::path(entry.file).filename()=="geometry.dict"&&entry.rawValue==(defaulted?"<omitted>":word),"pin wallTrace trace lost canonical value, source or strict default provenance");++traces;}
    for(const auto& document:input.fluid.documents)if(fs::path(document.file).filename()=="geometry.dict"){
        const auto marker="wallTrace "+word+";";
        wallTraceCheck(document.effectiveText.find(marker)!=std::string::npos&&
            (defaulted?document.sourceText.find("wallTrace")==std::string::npos:document.sourceText.find(marker)!=std::string::npos),
            "pin wallTrace source/effective dictionary lost the explicit law or changed omitted source bytes");++documents;}
    wallTraceCheck(traces==1&&documents==1,"pin wallTrace provenance is missing or duplicated");
}
VofElectrohydrodynamicsInput wallTraceRoundtrip(const fs::path& system,const VofElectrohydrodynamicsInput& input){
    fs::create_directories(system);
    for(const auto& document:input.fluid.documents)VofEhdInputContract::write(system/fs::path(document.file).filename(),document.effectiveText);
    return readVofElectrohydrodynamicsInput(system.string(),input.fluid.parallelThreads,Application::Electrospray);
}
void wallTraceReaders(const fs::path& root,const VofElectrohydrodynamicsInput& omitted,const VofElectrohydrodynamicsInput& omittedAgain){
    wallTraceCheck(HydrodynamicsInput::CsgPin{}.wallTrace==WallTraceLaw::StrictWholeCell,"CsgPin native typed default is not strictWholeCell");
    wallTraceEvidence(omitted,WallTraceLaw::StrictWholeCell,true);wallTraceEvidence(omittedAgain,WallTraceLaw::StrictWholeCell,false);
    unsigned sources=1,roundtrips=1,refused=0;
    for(const auto law:{WallTraceLaw::StrictWholeCell,WallTraceLaw::PrescribedWetDry}){
        auto f=fixture();const auto word=std::string(wallTraceName(law));
        replace(f.text[1],"dryBoundaryId 2; }","dryBoundaryId 2; wallTrace "+word+"; }");
        const auto system=C::save(root,"wall_trace_"+word,f);
        const auto input=readVofElectrohydrodynamicsInput(system.string(),3,Application::Electrospray);wallTraceEvidence(input,law,false);++sources;
        const auto again=wallTraceRoundtrip(root/("wall_trace_"+word+"_effective"),input);wallTraceEvidence(again,law,false);++roundtrips;
    }
    auto invalid=fixture();replace(invalid.text[1],"dryBoundaryId 2; }","dryBoundaryId 2; wallTrace prescribedwetdry; }");
    const auto bad=C::save(root,"wall_trace_unknown_exact",invalid);bool caught=false;
    try{(void)readVofElectrohydrodynamicsInput(bad.string(),3,Application::Electrospray);}
    catch(const std::exception& e){const auto text=std::string(e.what());caught=text.find("geometry.dict:")!=std::string::npos&&
        text.find("Geometry.Initial.Pin.wallTrace = 'prescribedwetdry'; expected strictWholeCell | prescribedWetDry")!=std::string::npos;}
    wallTraceCheck(caught,"unknown pin wallTrace was accepted or lost its exact word/source/allowed-value refusal");++refused;
    wallTraceCheck(sources==3&&roundtrips==3&&refused==1,"pin wallTrace reader coverage incomplete");
    std::cout<<"PIN_WALL_TRACE_INPUT_PASS source_cases="<<sources<<" effective_roundtrips="<<roundtrips<<" invalid_guards="<<refused<<" checks="<<wallTraceChecks
        <<" strict_omitted_explicit_prescribed_canonical_trace_source_effective_exact_unknown_refusal\n";
}
template<class Change>void reject(const fs::path& root,const char* name,const char* field,Change change){auto f=fixture();change(f);
    const auto system=C::save(root,name,f);bool caught=false;try{(void)readVofElectrohydrodynamicsInput(system.string(),{},Application::Electrospray);}
    catch(const std::exception& e){caught=true;check(std::string(e.what()).find(field)!=std::string::npos,"Electrospray profile refusal lost actual offending source field");}
    check(caught,"unsupported Electrospray profile or legacy App block was accepted");++guards;}
void run(const fs::path& requested){
    check(requested.is_absolute()&&requested.root_name().string()=="D:","Electrospray profile selector requires absolute D: output root");
    const auto root=requested/("electrospray_profile_"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));check(!fs::exists(root),"profile reader fixtures must be fresh");fs::create_directories(root);
    const auto generalSystem=C::save(root/"dispatch_general","system",F::fixture());const auto general=readVofElectrohydrodynamicsInput(generalSystem.string(),3);++readers;
    const auto electrospraySystem=C::save(root/"dispatch_electrospray","system",fixture());const auto input=readVofElectrohydrodynamicsInput(electrospraySystem.string(),3,Application::Electrospray);++readers;
    check(general.application==Application::General&&input.application==Application::Electrospray&&input.fluid.parallelThreads==3,"default General or explicit Electrospray/threads override changed");
    provenance(general,"VofElectrohydrodynamics");provenance(input,"Electrospray");samePhysics(general,input);
    const auto effective=root/"effective";fs::create_directories(effective);for(const auto& doc:input.fluid.documents)VofEhdInputContract::write(effective/fs::path(doc.file).filename(),doc.effectiveText);
    const auto again=readVofElectrohydrodynamicsInput(effective.string(),3,Application::Electrospray);++readers;provenance(again,"Electrospray");samePhysics(input,again);
    {bool caught=false;try{(void)readVofElectrohydrodynamicsInput(electrospraySystem.string());}catch(const std::exception& e){caught=std::string(e.what()).find("Run.app")!=std::string::npos;}
        check(caught,"default General reader silently aliased Run.app Electrospray");++guards;}
    {bool caught=false;try{(void)readVofElectrohydrodynamicsInput(generalSystem.string(),{},Application::Electrospray);}catch(const std::exception& e){caught=std::string(e.what()).find("Run.app")!=std::string::npos;}
        check(caught,"Electrospray reader silently aliased General Run.app");++guards;}
    reject(root,"missing_pin","Geometry.Initial",[](Fixture& f){replace(f.text[1]," Pin { pointId 6; wetBoundaryId 3; dryBoundaryId 2; }","");});
    reject(root,"missing_noSlip","Hydrodynamics.Boundary",[](Fixture& f){replace(f.text[5],"StationaryNoSlip { boundaryIds (2 3); }","");});
    reject(root,"legacy_app_block","Electrospray",[](Fixture& f){f.text[0]+=" Electrospray { enabled true; }";});
    reject(root,"hg_interface","Interface.method",[](Fixture& f){replace(f.text[3],"method vof;","method hermiteGals;");});
    reject(root,"unsupported_B1","Electric.method",[](Fixture& f){replace(f.text[6],"method vofFiniteVolume;","method vofSharpSurface;");});
    {auto f=C::fixture();replace(f.text[0],"app VofElectrohydrodynamics;","app Electrospray;");const auto system=C::save(root,"impermeable_not_spray",f);bool caught=false;
        try{(void)readVofElectrohydrodynamicsInput(system.string(),{},Application::Electrospray);}catch(const std::exception& e){caught=std::string(e.what()).find("Hydrodynamics.Boundary.type")!=std::string::npos;}
        check(caught,"Electrospray profile accepted closed flow without physical pressure reservoirs");++guards;}
    {Fixture f(true);replace(f.text[0],"app VofElectrohydrodynamics;","app Electrospray;");const auto system=f.save(root,"phi_not_CSG");bool caught=false;
        try{(void)readVofElectrohydrodynamicsInput(system.string(),{},Application::Electrospray);}catch(const std::exception& e){caught=std::string(e.what()).find("Geometry.Initial")!=std::string::npos;}
        check(caught,"Electrospray profile accepted a free phi droplet without fixed solid CSG");++guards;}
    {bool caught=false;try{(void)vofElectrohydrodynamicsApplicationName(static_cast<Application>(-1));}catch(const std::invalid_argument&){caught=true;}
        check(caught,"unknown typed production application profile had a silent fallback name");++guards;}
    check(readers==3&&guards==10,"Electrospray reader profile direct coverage incomplete");
    wallTraceReaders(root,input,again);
    std::ostringstream dispatch;dispatch<<"{\"Electrospray\":"<<std::quoted((root/"dispatch_electrospray").generic_string())
        <<",\"VofElectrohydrodynamics\":"<<std::quoted((root/"dispatch_general").generic_string())<<"}\n";
    VofEhdInputContract::write(root/"dispatch_cases.json",dispatch.str());
    std::cout<<"ELECTROSPRAY_PROFILE_INPUT_PASS readers="<<readers<<" guards="<<guards<<" checks="<<checks
        <<" explicit_General_default_Electrospray_A_CSG_pin_pressure_noSlip_source_effective_trace_no_HG_alias_no_solver_qualification dispatch_cases="<<(root/"dispatch_cases.json").string()<<'\n';
}
}
// END ELECTROSPRAY EXPLICIT READER PROFILE CONTRACT


// BEGIN ACTUAL ELECTROSPRAY HF PEAK OBSERVER CONTRACT
#include "Core/Geometry/FittedDomainGeometry.h"
#include "Core/Mesh/FaceSegmentGraph.h"
#include "Core/IO/StrictDictionary.h"
#include "IO/DictionaryParallel.h"
#include <cfenv>
#include <cerrno>
#include <exception>
#include <map>
#include <set>
#if defined(_M_X64) || defined(__SSE2__)
#include <xmmintrin.h>
#endif
namespace ElectrosprayHfPeakContract
{
namespace fs=std::filesystem;namespace H=LAMP::Hydrodynamics;namespace I=LAMP::Interface;
namespace J=VofEhdJointContract;using Time=LAMP::TimeLayer;
unsigned checks=0,released=0;
void check(bool value,const char* why){++checks;if(!value)throw std::runtime_error(why);}
struct Environment {
    std::fenv_t native{};int error=errno;
#if defined(_M_X64) || defined(__SSE2__)
    unsigned csr=_mm_getcsr();
#endif
    Environment(){check(std::fegetenv(&native)==0,"HF observer cannot capture caller floating environment");errno=error;}
    ~Environment() noexcept {
        if(std::fesetenv(&native)!=0)std::terminate();
#if defined(_M_X64) || defined(__SSE2__)
        _mm_setcsr(csr);
#endif
        errno=error;
    }
};
struct Observation {std::uint64_t identity=0;double time=0,minimum=0,maximum=0;
    Time acceptedClock;std::uint64_t acceptedIdentity=0;};
using WallTraceLaw=I::GeometricVofPublication::PinCertificate::WallTraceLaw;
unsigned pinLawChecks=0,pinInitial=0,pinMid=0,pinPhysicalEnd=0,pinAcceptedEnd=0;
void pinLawCheck(bool ok,const char* why){++pinLawChecks;if(!ok)throw std::runtime_error(why);}
void observePinLaw(const I::GeometricVofPublication& publication,const Time& expected,WallTraceLaw law,
    const char* label,const QuadtreeMesh* mesh,std::ostream& out)
{
    const auto pin=publication.pinnedCorner();
    pinLawCheck(pin&&pin->wallTraceLaw()==law&&pin->sourceIdentity()==publication.identity()&&
        J::sameClock(publication.timeLayer(),expected)&&J::sameClock(pin->sourceTime(),expected)&&
        pin->physicalPoint().x==1&&pin->physicalPoint().y==3,"actual App pin lost typed law, physical endpoint or same-owner full clock");
    if(mesh)pin->validateFor(*mesh,publication);
    const bool prescribed=law==WallTraceLaw::PrescribedWetDry&&pin->kind()==I::GeometricVofPublication::PinCertificate::Kind::ResolvedCorner;
    pinLawCheck(pin->prescribedWetWallTrace()==prescribed&&pin->prescribedDryWallTrace()==prescribed,
        "actual App pin wall receipts do not match declared law/current resolved contact (not bulk purity)");
    const auto clock=publication.timeLayer(),pinClock=pin->sourceTime();
    out<<"PIN_WALL_TRACE layer="<<label<<",law="<<ElectrosprayProfileInputContract::wallTraceName(pin->wallTraceLaw())
        <<",kind="<<(pin->kind()==I::GeometricVofPublication::PinCertificate::Kind::ResolvedCorner?"ResolvedCorner":"NativeCoincident")
        <<",wet_receipt="<<pin->prescribedWetWallTrace()<<",dry_receipt="<<pin->prescribedDryWallTrace()
        <<",gid="<<publication.identity().value()<<",pin_gid="<<pin->sourceIdentity().value()
        <<",clock=("<<clock.step<<','<<clock.begin<<','<<clock.end<<','<<static_cast<int>(clock.position)
        <<"),pin_clock=("<<pinClock.step<<','<<pinClock.begin<<','<<pinClock.end<<','<<static_cast<int>(pinClock.position)<<")\n";
}
std::string prescribedWallTrace(const VofElectrohydrodynamicsInput& formal,const std::string& source)
{
    pinLawCheck(formal.fluid.csgPin.has_value(),"target-band pin fixture requires the original declared CSG pin");
    if(formal.fluid.csgPin->wallTrace==WallTraceLaw::PrescribedWetDry)return source;
    pinLawCheck(formal.fluid.csgPin->wallTrace==WallTraceLaw::StrictWholeCell,"target-band fixture cannot translate an unknown wallTrace law");
    bool omitted=false;unsigned traces=0;
    for(const auto& entry:formal.fluid.trace.entries)if(entry.field=="Geometry.Initial.Pin.wallTrace"){omitted=entry.defaulted;++traces;}
    pinLawCheck(traces==1,"target-band fixture lacks one actual pin wallTrace declaration/default");
    auto result=source;
    if(omitted){const auto marker="dryBoundaryId "+std::to_string(formal.fluid.csgPin->dryBoundaryId)+";";
        VofEhdInputContract::replace(result,marker,marker+" wallTrace prescribedWetDry;");}
    else VofEhdInputContract::replace(result,"wallTrace strictWholeCell;","wallTrace prescribedWetDry;");
    return result;
}
std::optional<std::uint64_t> initialGeometryIdentity;
std::map<std::int64_t,Observation> observations;
bool targetBandPolicy(const VofElectrohydrodynamicsInput& input)
{
    return input.fluid.initialInterfaceBand&&input.fluid.initialBandScale==HydrodynamicsInput::InitialBandScale::TargetCell&&
        input.fluid.initialRefinementLevel==4&&input.fluid.initialInterfaceBandCells==3;
}
std::string regridSuffix(const std::string& mesh)
{
    const auto begin=mesh.find("Regrid");
    check(begin!=std::string::npos&&mesh.find("Regrid",begin+1)==std::string::npos,"HF input fixture needs one actual Regrid block");
    return mesh.substr(begin);
}
std::string prepareTargetMesh(const VofElectrohydrodynamicsInput& formal,const std::string& source)
{
    check(!formal.fluid.initialInterfaceBand||targetBandPolicy(formal),
        "HF target-band observer requires the formal L4/B3/targetCell policy; refusing a mesh override");
    if(formal.fluid.initialInterfaceBand)return source;
    auto result=source;VofEhdInputContract::replace(result,"initialMesh uniform;",
        "initialMesh interfaceBand; Refinement { level 4; bandCells 3; bandScale targetCell; }");
    check(regridSuffix(result)==regridSuffix(source),"HF target-band preparation changed dynamic Regrid bytes");
    return result;
}
std::string prepareUniformControlMesh(const VofElectrohydrodynamicsInput& input,const std::string& source)
{
    if(!input.fluid.initialInterfaceBand)return source;
    check(targetBandPolicy(input),"HF uniform control requires the qualified target-band source");
    auto result=source;VofEhdInputContract::replace(result,"initialMesh interfaceBand;","initialMesh uniform;");
    const auto begin=result.find("Refinement"),open=result.find('{',begin),end=result.find('}',open);
    check(begin!=std::string::npos&&open!=std::string::npos&&end!=std::string::npos&&
        result.find("Refinement",begin+1)==std::string::npos&&result.find('{',open+1)>end&&end<result.find("Regrid"),
        "HF uniform control cannot isolate the actual flat Refinement block");
    result.erase(begin,end-begin+1);
    check(regridSuffix(result)==regridSuffix(source),"HF uniform control changed dynamic Regrid bytes");
    return result;
}
void observe(const VofElectrohydrodynamicsTransaction& tx,const Time& accepted,const VofElectrohydrodynamicsInput& input,std::ostream& out)
{
    const Environment environment;
    const auto* step=tx.preparedStep();const auto* receipt=tx.remapReceipt();const auto* target=tx.accepted();
    check(step&&step->hasCapillaryLoad()&&receipt&&receipt->sourceMesh&&receipt->physicalEnd&&receipt->physicalStep==step&&target,
        "HF observer needs the actual pre-release physical source-mesh receipt");
    const auto& mesh=*receipt->sourceMesh;const auto& load=step->capillaryLoad();const auto& phase=load.phaseGeometry();
    const auto& publication=phase.vofPublication();const auto& layout=phase.layout();const auto clock=publication.timeLayer();
    receipt->physicalEnd->validateFor(mesh);
    load.validateVofFor(mesh,step->frozenMidpoint().pressureLoad(),publication.identity(),clock);
    check(clock.position==Time::Position::Midpoint&&clock.step==accepted.step&&clock.begin==accepted.begin&&clock.end==accepted.end&&
        J::sameClock(receipt->physicalEnd->timeLayer(),accepted)&&J::sameClock(target->timeLayer(),accepted),
        "HF observer conflated actual Mid force and accepted/remapped End clocks");
    pinLawCheck(input.fluid.csgPin.has_value(),"actual HF observer lacks the typed input pin wall law");
    const auto law=input.fluid.csgPin->wallTrace;
    observePinLaw(publication,clock,law,"Midpoint",&mesh,out);++pinMid;
    observePinLaw(receipt->physicalEnd->phaseGeometry().vofPublication(),accepted,law,"PhysicalEnd",&mesh,out);++pinPhysicalEnd;
    // The public transaction has no target-mesh getter. The actual accepted
    // publication still supplies its own identity/full clock and pin receipt.
    observePinLaw(target->phaseGeometry().vofPublication(),accepted,law,"AcceptedEnd",nullptr,out);++pinAcceptedEnd;
    const auto& faceK=load.heightFunctionFaceCurvatures();const auto& support=load.supportReports();
    const auto& physical=load.heightFunctionCellCurvatures();const auto& effective=load.heightFunctionCellEffectiveCurvatures();
    const auto& methods=load.heightFunctionCellCurvatureMethods();const auto& jets=load.heightFunctionCellJets();
    int faceMin=-1,faceMax=-1,rhsPeak=-1;
    for(int u=0;u<layout.velocityCount();++u)if(support[u].curvatureDefined&&std::isfinite(faceK[u])){
        if(faceMin<0||faceK[u]<faceK[faceMin])faceMin=u;
        if(faceMax<0||faceK[u]>faceK[faceMax])faceMax=u;}
    for(int u=0;u<layout.velocityCount();++u)if(!support[u].boundary&&std::isfinite(load.momentumRhsLoad()[u])&&
        (rhsPeak<0||std::abs(load.momentumRhsLoad()[u])>std::abs(load.momentumRhsLoad()[rhsPeak])))rhsPeak=u;
    check(faceMin>=0&&faceMax>=0&&rhsPeak>=0,"HF observer lacks actual supported curvature extrema/interior RHS peak");
    if(accepted.step==1){const auto& initial=step->beginState();const auto initialClock=initial.timeLayer();
        check(initialClock.step==0&&initialClock.begin==0&&initialClock.end==0&&initialClock.position==Time::Position::Initial,
            "HF observer first accepted step has no original Initial clock");
        initialGeometryIdentity=initial.phaseGeometry().vofIdentity().value();
        observePinLaw(initial.phaseGeometry().vofPublication(),initialClock,law,"Initial",&mesh,out);++pinInitial;}
    check(observations.emplace(accepted.step,Observation{publication.identity().value(),clock.time(),faceK[faceMin],faceK[faceMax],
        accepted,target->phaseGeometry().vofIdentity().value()}).second,
        "HF observer duplicated an accepted physical step");
    const auto& report=load.curvatureReport();const auto pin=publication.pinnedCorner();
    out<<"HF_FORCE_PEAK step="<<accepted.step<<",accepted_gid="<<target->phaseGeometry().vofIdentity().value()
        <<",physical_end_gid="<<receipt->physicalEnd->phaseGeometry().vofIdentity().value()<<",force_gid="<<publication.identity().value()
        <<",force_time="<<clock.time()<<",force_clock=("<<clock.step<<','<<clock.begin<<','<<clock.end<<','<<static_cast<int>(clock.position)
        <<"),face_range=("<<faceK[faceMin]<<','<<faceK[faceMax]<<"),height="<<report.heightFunctionCells
        <<",extended="<<report.heightFunctionExtendedCells<<",refined="<<report.heightFunctionRefinedCells<<",fit="<<report.heightFunctionFitCells
        <<",pinned="<<report.pinnedCornerHeightCells<<",coarse_aligned="<<report.heightFunctionCoarseAlignedCells;
    if(pin)out<<",pin=(native:"<<pin->nativeCell()<<",point:"<<pin->physicalPoint().x<<','<<pin->physicalPoint().y<<')';
    out<<'\n';std::set<int> cells,printed;
    for(int u:{faceMin,faceMax,rhsPeak}){if(!printed.insert(u).second)continue;const auto f=layout.domain().face(mesh,layout.nativeFace(u));const auto& r=support[u];
        out<<" FACE min_kappa="<<(u==faceMin)<<",max_kappa="<<(u==faceMax)<<",max_interior_abs_rhs="<<(u==rhsPeak)
            <<",velocity="<<u<<",native_row="<<f.nativeRow<<",axis="<<int(f.horizontal)
            <<",center=("<<f.center.x<<','<<f.center.y<<"),native_cells=("<<f.negativeCell<<','<<f.positiveCell<<"),kappa="<<faceK[u]
            <<",gradient="<<load.colorGradient()[u]<<",rhs="<<load.momentumRhsLoad()[u]<<",coincident="<<r.coincidentPhaseHeight
            <<",amr="<<r.crossesAmr<<",support_count="<<r.supportCount<<",boundary="<<r.boundary
            <<",subresolution_gradient="<<r.subResolutionGradient<<",gradient_roundoff_scale="<<r.gradientRoundoffScale
            <<",support_status="<<static_cast<int>(r.status)<<'\n';
        for(int native:{f.negativeCell,f.positiveCell})if(native>=0){const int p=layout.pressureDof(native);if(p>=0)cells.insert(p);}}
    for(int p:cells){const int native=layout.nativeCell(p);const auto& cell=publication.plicCells()[native];const auto& b=cell.rectangle;
        out<<" CELL pressure="<<p<<",native="<<native<<",level="<<mesh.cells(mesh.cellIndexToLabel[native]).level()
            <<",rectangle=[("<<b.lower.x<<','<<b.lower.y<<"),("<<b.upper.x<<','<<b.upper.y<<")],inventory=("<<cell.inventory[0]<<','<<cell.inventory[1]
            <<"),measured=("<<cell.measured[0]<<','<<cell.measured[1]<<"),normal=("<<cell.normal.x<<','<<cell.normal.y<<"),offset="<<cell.offset
            <<",resolved="<<cell.resolved<<",method="<<static_cast<int>(methods[p])<<",physical_available="<<std::isfinite(physical[p])
            <<",physical_kappa="<<physical[p]<<",effective_kappa="<<effective[p];
        if(pin){const auto point=pin->physicalPoint(),center=(b.lower+b.upper)*.5;const auto& contact=publication.plicCells()[pin->nativeCell()].rectangle;
            const auto size=contact.upper-contact.lower;out<<",center_minus_pin=("<<center.x-point.x<<','<<center.y-point.y
                <<"),distance_in_contact_cells=("<<(center.x-point.x)/size.x<<','<<(center.y-point.y)/size.y<<')';}
        if(cell.resolved&&cell.measured[0]>0&&cell.measured[1]>0){const auto segment=LAMP::Geometry::PlanarCellCut::interfaceSegment(b,cell.normal,cell.offset);
            const long double cx=(static_cast<long double>(b.lower.x)+b.upper.x)/2,cy=(static_cast<long double>(b.lower.y)+b.upper.y)/2;
            if(segment)out<<",actual_mixed_segment=[("<<cx+(*segment)[0].x<<','<<cy+(*segment)[0].y<<"),("
                <<cx+(*segment)[1].x<<','<<cy+(*segment)[1].y<<")]";}
        if(const auto& selected=jets[p]){const auto& jet=*selected;const auto parts=H::HydroBalancedCsfLoad::axisymmetricHeightCurvatures(jet);
            out<<",jet={axis:"<<jet.heightAxis<<",outward:"<<jet.outward<<",origin:"<<jet.origin<<",height:"<<jet.height<<",first:"<<jet.first
                <<",second:"<<jet.second<<",tangent:"<<jet.tangentCoordinate<<",normal_spacing:"<<jet.normalSpacing<<",tangent_spacing:"<<jet.tangentSpacing
                <<",coarse_aligned:"<<jet.coarseAligned<<",axis_regular:"<<jet.axisRegular<<",axis_reflected:"<<jet.axisReflected
                <<",south_reflected:"<<jet.southEven90Reflected<<",pinned:"<<jet.pinnedCornerConstrained<<",pin_bin:"<<jet.pinnedCornerBin
                <<",meridional:"<<parts[0]<<",hoop:"<<parts[1]<<",radius:"
                <<(jet.heightAxis==0?jet.origin+jet.outward*jet.height:jet.tangentCoordinate)<<'}';}
        out<<'\n';
    }
    out.flush();
    const auto peak=layout.domain().face(mesh,layout.nativeFace(faceMax));
    std::cout<<"ELECTROSPRAY_HF_PEAK step="<<accepted.step<<" force_gid="<<publication.identity().value()<<" force_time="<<clock.time()
        <<" face_range=("<<faceK[faceMin]<<','<<faceK[faceMax]<<") max_native_face="<<peak.nativeRow
        <<" max_center=("<<peak.center.x<<','<<peak.center.y<<") max_interior_abs_rhs="<<std::abs(load.momentumRhsLoad()[rhsPeak])
        <<" rhs_peak_native_face="<<layout.nativeFace(rhsPeak)<<'\n';
}

void initialTargetReaders(const fs::path& root,const VofElectrohydrodynamicsInput& input,
    const std::array<const char*,8>& names)
{
    using Scale=HydrodynamicsInput::InitialBandScale;
    const auto sourceMesh=J::read(root/"system/mesh.dict");
    const auto uniformMesh=prepareUniformControlMesh(input,sourceMesh);const auto uniformSystem=root/"uniform_control";
    fs::create_directories(uniformSystem);
    for(const char* name:names)VofEhdInputContract::write(uniformSystem/name,
        std::string(name)=="mesh.dict"?uniformMesh:J::read(root/"system"/name));
    const auto uniform=readVofElectrohydrodynamicsInput(uniformSystem.string(),2,VofElectrohydrodynamicsInput::Application::Electrospray);
    check(!uniform.fluid.initialInterfaceBand&&regridSuffix(uniformMesh)==regridSuffix(sourceMesh),
        "explicit uniform control retained Refinement or changed Regrid");
    const auto restoredMesh=prepareTargetMesh(uniform,uniformMesh);const auto restoredSystem=root/"uniform_target_control";
    fs::create_directories(restoredSystem);
    for(const char* name:names)VofEhdInputContract::write(restoredSystem/name,
        std::string(name)=="mesh.dict"?restoredMesh:J::read(uniformSystem/name));
    const auto restored=readVofElectrohydrodynamicsInput(restoredSystem.string(),2,VofElectrohydrodynamicsInput::Application::Electrospray);
    check(targetBandPolicy(restored)&&prepareTargetMesh(input,sourceMesh)==sourceMesh,
        "uniform opt-in or already-target mesh byte preservation failed");
    ElectrosprayProfileInputContract::samePhysics(input,restored);
    check(!uniform.fluid.initialInterfaceBand&&input.fluid.initialInterfaceBand&&input.fluid.initialBandScale==Scale::TargetCell&&
        input.fluid.initialRefinementLevel==4&&input.fluid.initialInterfaceBandCells==3,
        "targetCell opt-in changed the uniform default or lost its explicit level/width/scale");
    ElectrosprayProfileInputContract::samePhysics(uniform,input);
    unsigned traces=0;
    for(const auto& entry:input.fluid.trace.entries)if(entry.field=="Mesh.Refinement.bandScale"){
        const auto* value=std::get_if<std::string>(&entry.value);
        check(value&&*value=="targetCell"&&entry.used&&!entry.defaulted,"targetCell trace lost explicit consumed source");++traces;}
    check(traces==1,"targetCell trace is missing or duplicated");
    fs::create_directories(root/"roundtrip");
    for(const auto& doc:input.fluid.documents)VofEhdInputContract::write(root/"roundtrip"/fs::path(doc.file).filename(),doc.effectiveText);
    const auto again=readVofElectrohydrodynamicsInput((root/"roundtrip").string(),2,VofElectrohydrodynamicsInput::Application::Electrospray);
    check(again.fluid.initialInterfaceBand&&again.fluid.initialBandScale==Scale::TargetCell&&again.fluid.initialRefinementLevel==4&&
        again.fluid.initialInterfaceBandCells==3,"targetCell source/effective roundtrip lost its actual policy");
    ElectrosprayProfileInputContract::samePhysics(input,again);
    unsigned refused=0;
    for(const char* scale:{"baseCell","currentCell"}){const auto system=root/("guard_"+std::string(scale));fs::create_directories(system);
        for(const char* name:names){auto text=J::read(root/"system"/name);
            if(std::string(name)=="mesh.dict")VofEhdInputContract::replace(text,"bandScale targetCell;",std::string("bandScale ")+scale+";");
            VofEhdInputContract::write(system/name,text);}
        bool rejected=false;try{(void)readVofElectrohydrodynamicsInput(system.string(),2,VofElectrohydrodynamicsInput::Application::Electrospray);}
        catch(const std::exception&){rejected=true;}check(rejected,"CSG base/current band scale silently entered targetCell");++refused;}
    VofEhdInputContract::Fixture phi(true);const auto baseline=phi.save(root,"phi_default");
    (void)readVofElectrohydrodynamicsInput(baseline.string());
    phi.text[2]="Mesh { resolution (8 8); initialMesh interfaceBand; Refinement { level 4; bandCells 3; bandScale targetCell; } }";
    const auto foreign=phi.save(root,"guard_non_csg");bool rejected=false;
    try{(void)readVofElectrohydrodynamicsInput(foreign.string());}catch(const std::exception&){rejected=true;}
    check(rejected,"non-CSG targetCell silently selected another initialization algorithm");++refused;
    std::cout<<"CSG_INITIAL_TARGET_READERS roundtrip=1 guards="<<refused<<" default_uniform=1 target_level=4 target_band=3\n";
}
void inspectTargetInitial(const VofElectrohydrodynamicsTransaction& tx,const Time& accepted,
    const VofElectrohydrodynamicsInput& input,const fs::path& root,std::ostream& out)
{
    const Environment environment;const auto* step=tx.preparedStep();const auto* receipt=tx.remapReceipt();
    check(accepted.step==1&&step&&receipt&&receipt->sourceMesh,"targetCell Initial check lost its first actual source mesh");
    const auto& mesh=*receipt->sourceMesh;const auto& begin=step->beginState();const auto& pub=begin.phaseGeometry().vofPublication();
    const auto clock=pub.timeLayer();const auto& domain=pub.domain();const auto& layout=begin.values().layout();
    begin.validateFor(mesh);
    check(clock.step==0&&clock.begin==0&&clock.end==0&&clock.position==Time::Position::Initial&&J::sameClock(clock,begin.timeLayer())&&
        step->transport().sourceIdentity()==pub.identity(),"targetCell did not expose one true Initial source before physical transport");
    check(mesh.cellIndexNum()>input.fluid.nx*input.fluid.ny&&!domain.uniformLeaves()&&domain.axisymmetric()&&
        domain.cellRowByIndex()==layout.domain().cellRowByIndex()&&domain.cellIndexByRow()==layout.domain().cellIndexByRow(),
        "targetCell Initial mesh is still coarse/uniform or lost exact selected pressure maps");
    const auto solids=FittedDomainGeometry::selectFixedSolids(mesh,*input.fluid.csgGeometry,*input.fluid.csgSampler,input.fluid.csgExcludedPhases);
    FaceSegmentGraph all;all.build(mesh);const auto complete=NativeLeafDomain::build(mesh,all,true);
    unsigned fine=0,coarse=0,inactive=0,mixed=0,pureJumps=0;long double liquid=0,gas=0,capacity=0;
    for(int native=0;native<mesh.cellIndexNum();++native){const auto& cell=mesh.cells(mesh.cellIndexToLabel[native]);const auto& phase=pub.plicCells()[native];
        check(cell.level()<=input.fluid.initialRefinementLevel,"targetCell Initial exceeded its declared target level");
        fine+=cell.level()==input.fluid.initialRefinementLevel;coarse+=cell.level()<input.fluid.initialRefinementLevel;
        const bool active=domain.cellRowByIndex()[native]>=0;
        check(active==(solids.excludedPhaseByCell[native]==0),"targetCell final Initial selected mask differs from original fitted solid");
        if(!active){++inactive;check(phase.inventory[0]==0&&phase.inventory[1]==0&&phase.measured[0]==0&&phase.measured[1]==0&&
            !phase.resolved&&pub.cells()[native].capacity==0&&pub.alpha()[native]==0,"targetCell excluded solid acquired a phantom phase");continue;}
        check(layout.nativeCell(layout.pressureDof(native))==native&&phase.inventory[0]>=0&&phase.inventory[1]>=0,"targetCell native map or independent stock changed");
        liquid+=phase.inventory[0];gas+=phase.inventory[1];capacity+=pub.cells()[native].capacity;
        if(phase.measured[0]>0&&phase.measured[1]>0){++mixed;check(cell.level()==input.fluid.initialRefinementLevel,"targetCell left an actual mixed Initial PLIC at coarse resolution");}}
    for(const auto f:complete.faces(mesh))if(f.negativeCell>=0&&f.positiveCell>=0){
        const auto& a=mesh.cells(mesh.cellIndexToLabel[f.negativeCell]);const auto& b=mesh.cells(mesh.cellIndexToLabel[f.positiveCell]);
        check(std::abs(a.level()-b.level())<=1,"targetCell final complete native mesh is not locally 2:1 across solid/fluid mask");}
    const auto purePhase=[&](int native){const auto& c=pub.plicCells()[native];return c.inventory[0]>0&&c.inventory[1]==0?0:c.inventory[0]==0&&c.inventory[1]>0?1:-1;};
    for(const auto f:domain.faces(mesh))if(f.negativeCell>=0&&f.positiveCell>=0){const int a=purePhase(f.negativeCell),b=purePhase(f.positiveCell);
        if(a>=0&&b>=0&&a!=b){++pureJumps;check(mesh.cells(mesh.cellIndexToLabel[f.negativeCell]).level()==input.fluid.initialRefinementLevel&&
            mesh.cells(mesh.cellIndexToLabel[f.positiveCell]).level()==input.fluid.initialRefinementLevel,"targetCell left actual pure-jump Initial seeds below target level");}}
    check(fine>0&&coarse>0&&inactive>0&&mixed>0&&liquid>0&&gas>0,"targetCell Initial lost local coarse/fine/selected/two-phase witnesses");
    VofEhdOpenJointContract::bound(capacity-328*std::acos(-1.L),328*std::acos(-1.L),1e-12L,"targetCell changed the original physical RZ fluid-domain volume");
    const auto pin=pub.pinnedCorner();check(pin&&pin->sourceIdentity()==pub.identity()&&J::sameClock(pin->sourceTime(),clock)&&
        pin->physicalPoint().x==1&&pin->physicalPoint().y==3,"targetCell Initial lost the original same-source pin/full Initial clock");
    pin->validateFor(mesh,pub);check(mesh.cells(mesh.cellIndexToLabel[pin->nativeCell()]).level()==input.fluid.initialRefinementLevel,"targetCell contact source stayed coarse");
    const auto wet=domain.face(mesh,pin->wetNativeRow()),dry=domain.face(mesh,pin->dryNativeRow());
    check(!wet.horizontal&&wet.center.x==1&&wet.center.y<3&&dry.horizontal&&dry.center.y==3&&dry.center.x>1,
        "targetCell pin wall legs differ from original nozzle geometry");
    // Test-only direct public initialization on the final source mesh. Exact
    // stock equality distinguishes one fine CSG initialization from an
    // already-published coarse PLIC followed by conservative prolongation.
    const auto reference=I::GeometricVofInitialization::prepare(mesh,domain,[&](Vector2D p){
        return input.fluid.csgSampler->phaseValue(input.fluid.negativePhase,p);},input.fluid.vofInitialization);
    for(int native=0;native<mesh.cellIndexNum();++native){check(pub.plicCells()[native].inventory==reference.plicCells()[native].inventory,
        "targetCell Initial stocks are not the sole same-physical-source initializer on its final mesh");
        if(native!=pin->nativeCell())check(pub.plicCells()[native].normal.x==reference.plicCells()[native].normal.x&&
            pub.plicCells()[native].normal.y==reference.plicCells()[native].normal.y&&pub.plicCells()[native].offset==reference.plicCells()[native].offset,
            "targetCell Initial geometry inherited a different coarse plane before publication");}
    std::istringstream index(J::read(root/"output/meta/frames.csv"));std::string line;unsigned initialFrames=0,records=0;
    std::getline(index,line);
    while(std::getline(index,line)){if(line.empty())continue;++records;std::istringstream row(line);std::vector<std::string> values;std::string value;
        while(std::getline(row,value,','))values.push_back(value);check(values.size()==9,"targetCell actual writer index columns changed");
        if(std::stoll(values[1])==0){++initialFrames;check(std::stod(values[2])==0&&std::stoull(values[3])==pub.identity().value(),
            "targetCell Initial output identity differs from the original transported owner");}}
    check(initialFrames==1&&records==1,"targetCell pre-output observer lacks exactly one published actual Initial frame");
    out<<"INITIAL_TARGET_CELL native="<<mesh.cellIndexNum()<<",active="<<domain.cellIndexByRow().size()<<",fine="<<fine<<",coarse="<<coarse
        <<",inactive="<<inactive<<",mixed="<<mixed<<",pure_jump_seeds="<<pureJumps<<",initial_gid="<<pub.identity().value()
        <<",initial_clock=(0,0,0,0),liquid="<<liquid<<",gas="<<gas<<",same_final_mesh_initializer_stock_exact=1\n";out.flush();
}
void run(const fs::path& requested,bool targetBand=false,bool initialOnly=false,bool inputOnly=false)
{
    check(requested.is_absolute()&&requested.root_name().string()=="D:","HF peak observer requires absolute D: output root");
    const unsigned expectedSteps=initialOnly?1:4;
    auto repo=fs::path(__FILE__);for(int parent=0;parent<4;++parent)repo=repo.parent_path();const auto source=repo/"cases/electrospray/system";
    const auto root=requested/(std::string(initialOnly?"electrospray_initial_target_band_":targetBand?"electrospray_hf_peaks_target_band_":"electrospray_hf_peaks_")+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    check(!fs::exists(root),"HF peak observer output must be fresh");fs::create_directories(root/"source_system");fs::create_directories(root/"system");
    const std::array<const char*,8> names{{"run.dict","geometry.dict","mesh.dict","interface.dict","materials.dict","hydro.dict","electric.dict","nozzle.csg"}};
    std::array<std::string,8> originals;
    for(std::size_t k=0;k<names.size();++k){originals[k]=J::read(source/names[k]);VofEhdInputContract::write(root/"source_system"/names[k],originals[k]);
        VofEhdInputContract::write(root/"system"/names[k],originals[k]);}
    const auto formal=readVofElectrohydrodynamicsInput((root/"source_system").string(),2,VofElectrohydrodynamicsInput::Application::Electrospray);
    auto observedRun=originals[0];VofEhdInputContract::replace(observedRun,"endTime 20;",initialOnly?"endTime 5e-6;":"endTime 2e-5;");
    VofEhdInputContract::replace(observedRun,"writeEverySteps "+std::to_string(formal.fluid.run.writeEverySteps)+";","writeEverySteps 1;");
    VofEhdInputContract::replace(observedRun,"everySteps "+std::to_string(formal.fluid.run.logEverySteps)+";","everySteps 1;");
    VofEhdInputContract::write(root/"system/run.dict",observedRun);
    const auto observedMesh=targetBand?prepareTargetMesh(formal,originals[2]):originals[2];
    if(observedMesh!=originals[2])VofEhdInputContract::write(root/"system/mesh.dict",observedMesh);
    const auto observedGeometry=targetBand?prescribedWallTrace(formal,originals[1]):originals[1];
    if(observedGeometry!=originals[1])VofEhdInputContract::write(root/"system/geometry.dict",observedGeometry);
    std::ostringstream overrides;overrides<<"Eight formal inputs copied byte-for-byte into source_system. No formal input is edited.\n"
        <<"Time.endTime:20->"<<(initialOnly?"5e-6":"2e-5")<<"; Output.Fields.writeEverySteps:"<<formal.fluid.run.writeEverySteps
        <<"->1; Log.everySteps:"<<formal.fluid.run.logEverySteps<<"->1.\n"
        <<"Explicit process resource threadsOverride=2 matches qualified v7. No dt/physical coefficient/stock/HF gate change.\n";
    if(targetBand){if(observedMesh==originals[2])overrides<<"Copied mesh already L4/B3/targetCell; no mesh override; original mesh bytes preserved.\n";
        else overrides<<"Copied mesh initialMesh:uniform->interfaceBand; Refinement { level4; bandCells3; bandScale targetCell; }. Dynamic Regrid unchanged.\n";}
    if(targetBand)overrides<<"Copied Geometry.Initial.Pin.wallTrace explicitly prescribedWetDry (retained if already declared); formal/source_system bytes unchanged.\n";
    if(inputOnly)overrides<<"Input-only preparation/reader contract; no writer, transaction or physical step is run.\n";
    overrides<<"hf_peaks.txt reads the actual retained Mid/sourceMesh, never a recomputed End force.\n";
    VofEhdInputContract::write(root/"observation_overrides.txt",overrides.str());
    const auto input=readVofElectrohydrodynamicsInput((root/"system").string(),2,VofElectrohydrodynamicsInput::Application::Electrospray);
    if(initialOnly)initialTargetReaders(root,input,names);
    if(targetBand){
        ElectrosprayProfileInputContract::wallTraceEvidence(input,WallTraceLaw::PrescribedWetDry,false);
        const auto again=ElectrosprayProfileInputContract::wallTraceRoundtrip(root/"pin_wall_trace_effective",input);
        ElectrosprayProfileInputContract::wallTraceEvidence(again,WallTraceLaw::PrescribedWetDry,false);
    }
    check(input.application==VofElectrohydrodynamicsInput::Application::Electrospray&&input.fluid.parallelThreads==2&&input.fluid.run.endTime==(initialOnly?5e-6:2e-5)&&
        input.fluid.run.dt==5e-6&&input.fluid.run.adaptive&&input.fluid.vofRegrid,"HF observer bypassed the actual Electrospray adaptive/AMR input");
    const auto checkInputBytes=[&]{
        for(std::size_t k=0;k<names.size();++k){check(J::read(source/names[k])==originals[k]&&J::read(root/"source_system"/names[k])==originals[k],
            "HF observer changed the formal case or its frozen source copy");check(J::read(root/"system"/names[k])==(k==0?observedRun:k==1?observedGeometry:k==2?observedMesh:originals[k]),
            "HF observer changed an unapproved copied input");}
    };
    if(inputOnly){checkInputBytes();std::cout<<"ELECTROSPRAY_INITIAL_TARGET_BAND_INPUT_PASS checks="<<checks
        <<" formal_mesh_preserved="<<(observedMesh==originals[2])<<" independent_uniform_control=1 source_effective_roundtrip=1 reader_guards=3 no_physical_steps output="<<root.string()<<'\n';return;}
    AppRuntime::configureParallelThreads(input.fluid.parallelThreads,"command_line");
    VofElectrohydrodynamicsOutput writer(input,{root});VofElectrohydrodynamicsTransaction tx(input,writer);std::string failure;
    std::ofstream peaks(root/"hf_peaks.txt");peaks.exceptions(std::ios::failbit|std::ios::badbit);peaks<<std::setprecision(std::numeric_limits<long double>::max_digits10);
    const HydrodynamicsRunIO io{[&]{writer.begin();},
        [&](auto attempted,auto step,double time,double oldDt,double newDt,const std::string& why){writer.retry(attempted,step,time,oldDt,newDt,why);},
        [&](auto attempted,auto step,double time,const char* stage,const std::string& why){failure=std::string(stage)+": "+why;writer.failure(attempted,step,time,stage,why);},
        [&]{check(observations.size()==expectedSteps&&released==expectedSteps,"HF observer short window missed an actual accepted/released step");
            pinLawCheck(pinInitial==1&&pinMid==expectedSteps&&pinPhysicalEnd==expectedSteps&&pinAcceptedEnd==expectedSteps,
                "pin wall law observer missed an actual Initial/Midpoint/physical End/accepted End owner");writer.complete();},
        [&](const auto& report,const auto& clock,int retries){writer.timeStep(report,clock,retries);
            if(targetBand&&clock.step==1)inspectTargetInitial(tx,clock,input,root,peaks);observe(tx,clock,input,peaks);},
        [&](const auto& clock,double seconds){check(tx.accepted()&&!tx.preparedStep()&&!tx.remapReceipt(),"HF observer physical receipt survived release");
            ++released;writer.stepCompleted(clock,seconds);}};
    RunConfig config;config.app="Electrospray";config.threadsOverride=2;
    const int result=runHydrodynamics(input.fluid,config,tx,io);
    if(result!=0)throw std::runtime_error("HF actual Electrospray workflow failed: "+failure);
    check(!tx.accepted()&&!tx.acceptedCharge()&&!tx.acceptedField(),"HF observer retained an accepted owner after workflow clear");
    std::istringstream rows(J::read(root/"output/series/vof_csf.dat"));std::string line;unsigned matched=0;
    while(std::getline(rows,line)){std::istringstream values(line);std::int64_t step=0;double accepted=0,time=0,lo=0,hi=0;std::uint64_t id=0;
        if(!(values>>step>>accepted>>time>>id>>lo>>hi))continue;const auto found=observations.find(step);
        check(found!=observations.end()&&found->second.identity==id&&found->second.time==time&&found->second.minimum==lo&&found->second.maximum==hi,
            "HF observer source identity/time/face extrema differ from actual writer vof_csf");++matched;}
    check(matched==expectedSteps,"HF observer did not match every actual writer force row");
    // timeStep precedes tx.output. Final frame-count/End validation belongs
    // here, after all accepted output and release callbacks have completed.
    std::istringstream frameRows(J::read(root/"output/meta/frames.csv"));unsigned frames=0,initialFrames=0;double finalFrameTime=0;
    std::getline(frameRows,line);
    while(std::getline(frameRows,line)){if(line.empty())continue;std::istringstream row(line);std::vector<std::string> values;std::string value;
        while(std::getline(row,value,','))values.push_back(value);check(values.size()==9,"HF observer final frame-index columns changed");
        const auto frame=std::stoul(values[0]);const auto step=std::stoll(values[1]);const auto time=std::stod(values[2]);const auto id=std::stoull(values[3]);
        check(frame==frames&&step==static_cast<std::int64_t>(frames),"HF observer final accepted frames are missing, duplicated or unordered");
        if(step==0){++initialFrames;check(initialGeometryIdentity&&id==*initialGeometryIdentity&&time==0,
            "HF observer final Initial frame lost the original source identity/time");}
        else{const auto found=observations.find(step);check(found!=observations.end()&&found->second.acceptedClock.position==Time::Position::End&&
            time==found->second.acceptedClock.time()&&id==found->second.acceptedIdentity,
            "HF observer final End frame differs from its actual accepted full clock/geometry");}
        finalFrameTime=time;++frames;
    }
    check(initialFrames==1&&frames==expectedSteps+1&&finalFrameTime==input.fluid.run.endTime,
        "HF observer final publication lacks one Initial plus every actual End frame");
    checkInputBytes();
    peaks.close();std::cout<<(initialOnly?"ELECTROSPRAY_INITIAL_TARGET_BAND_PASS actual_steps=":"ELECTROSPRAY_HF_PEAKS_PASS actual_steps=")<<observations.size()<<" writer_force_rows="<<matched<<" checks="<<checks
        <<" pin_initial="<<pinInitial<<" pin_mid="<<pinMid<<" pin_physical_end="<<pinPhysicalEnd<<" pin_accepted_end="<<pinAcceptedEnd<<" pin_law_checks="<<pinLawChecks
        <<" original_Mid_source_mesh_public_getters_no_production_change_no_End_HF output="<<root.string()<<'\n';
}
}
// END ACTUAL ELECTROSPRAY HF PEAK OBSERVER CONTRACT


// BEGIN TYPED A ELECTRIC BACKEND INPUT CONTRACT
namespace VofEhdBackendInputContract {
namespace fs=std::filesystem;using Fixture=VofEhdInputContract::Fixture;
using Backend=LAMP::Electric::ElectricVofControls::Backend;
unsigned checks=0,guards=0,readers=0,roundtrips=0;
void check(bool ok,const char* why){++checks;if(!ok)throw std::runtime_error(why);}
void evidence(const VofElectrohydrodynamicsInput& input,Backend expected,bool defaulted){
    const std::string word=expected==Backend::SuperLU?"superlu":"gmres";
    check(input.method==VofElectrohydrodynamicsInput::Method::VofFiniteVolume&&input.controls.backend==expected,
        "A backend reader lost its typed explicit/default selection");
    check(input.controls.relativeLinearTolerance==1e-11&&input.controls.relativeConservationTolerance==1e-10,
        "A backend choice changed original strict tolerances");
    unsigned traces=0,documents=0;
    for(const auto& e:input.fluid.trace.entries)if(e.field=="Electric.Solver.backend"){
        const auto* value=std::get_if<std::string>(&e.value);
        check(value&&*value==word&&e.used&&e.defaulted==defaulted&&e.consumer=="ElectricVofField.strictIntegratedBackend"
            &&fs::path(e.file).filename()=="electric.dict"&&e.rawValue==(defaulted?"<omitted>":word),
            "A backend trace lost exact source/default/value/consumer provenance");++traces;}
    for(const auto& d:input.fluid.documents)if(fs::path(d.file).filename()=="electric.dict"){
        const auto marker="backend "+word+";";
        check(d.effectiveText.find(marker)!=std::string::npos&&
            (defaulted?d.sourceText.find("backend")==std::string::npos:d.sourceText.find(marker)!=std::string::npos),
            "A backend source/effective selection changed omitted source or lost declared word");++documents;}
    check(traces==1&&documents==1,"A backend provenance missing or duplicated");
}
VofElectrohydrodynamicsInput reread(const fs::path& system,const VofElectrohydrodynamicsInput& input){
    fs::create_directories(system);
    for(const auto& d:input.fluid.documents)VofEhdInputContract::write(system/fs::path(d.file).filename(),d.effectiveText);
    return readVofElectrohydrodynamicsInput(system.string());
}
void refused(const fs::path& root,const char* name,Fixture f){
    const auto path=f.save(root,name);bool caught=false;
    try{(void)readVofElectrohydrodynamicsInput(path.string());}
    catch(const std::exception& e){const auto message=std::string(e.what());
        caught=message.find("electric.dict:")!=std::string::npos&&message.find("Electric.Solver.backend")!=std::string::npos;}
    check(caught,"backend reader accepted unsupported word or widened B1 vocabulary");++guards;
}
void run(const fs::path& requested){
    check(requested.is_absolute()&&requested.root_name().string()=="D:","backend reader contract requires absolute D: output root");
    const auto root=requested/("field_backend_input_"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    check(!fs::exists(root),"backend reader fixture path must be fresh");fs::create_directories(root);
    check(LAMP::Electric::ElectricVofControls{}.backend==Backend::Gmres,"module backend default changed from GMRES");
    for(bool rz:{false,true})for(const std::string word:{"","gmres","superlu"}){
        Fixture f(rz);if(!word.empty())VofEhdInputContract::replace(f.text[6],"Electric {","Electric { Solver { backend "+word+"; }");
        const auto name=std::string(rz?"rz_":"xy_")+(word.empty()?"omitted":word);
        const auto input=readVofElectrohydrodynamicsInput(f.save(root,name).string(),3);
        const auto expected=word=="superlu"?Backend::SuperLU:Backend::Gmres;evidence(input,expected,word.empty());++readers;
        for(std::size_t k=0;k<VofEhdInputContract::names.size();++k)
            check(input.fluid.documents[k].sourceText==f.text[k],"backend reader changed one of seven actual source documents");
        const auto again=reread(root/(name+"_effective"),input);evidence(again,expected,false);++roundtrips;
        check(again.fluid.radial==rz&&again.controls.chargeAdvection==input.controls.chargeAdvection&&
            again.materials.epsilonLiquid==input.materials.epsilonLiquid&&again.initialChargeDensity({.25,.25})==-.5,
            "backend effective roundtrip changed material/charge/metric");
    }
    Fixture unknown;VofEhdInputContract::replace(unknown.text[6],"Electric {","Electric { Solver { backend auto; }");
    refused(root,"unknown_backend",unknown);
    Fixture b1;b1.text[2]="Mesh { resolution (8 8); initialMesh uniform; }";
    VofEhdInputContract::replace(b1.text[6],"vofFiniteVolume","vofSharpSurface");
    VofEhdInputContract::replace(b1.text[6],"chargeDensity","surfaceChargeDensity");
    const auto sharp=readVofElectrohydrodynamicsInput(b1.save(root,"b1_original").string());
    const auto sharpAgain=reread(root/"b1_original_effective",sharp);
    for(const auto* input:{&sharp,&sharpAgain}){
        check(input->method==VofElectrohydrodynamicsInput::Method::VofSharpSurface&&input->initialSurfaceChargeDensity({.25,.25})==-.5,
            "B1 original surface reader or effective roundtrip changed");
        for(const auto& d:input->fluid.documents)if(fs::path(d.file).filename()=="electric.dict")
            check(d.effectiveText.find("backend")==std::string::npos,"B1 effective vocabulary gained A backend");
        for(const auto& e:input->fluid.trace.entries)check(e.field!="Electric.Solver.backend","B1 gained an unused A backend trace");
    }
    for(const std::string word:{"gmres","superlu"}){auto wrong=b1;
        VofEhdInputContract::replace(wrong.text[6],"Electric {","Electric { Solver { backend "+word+"; }");
        refused(root,("b1_"+word).c_str(),wrong);}
    check(readers==6&&roundtrips==6&&guards==3,"backend reader coverage incomplete");
    std::cout<<"VOF_EHD_FIELD_BACKEND_INPUT_PASS source="<<readers<<" roundtrip="<<roundtrips<<" guards="<<guards
        <<" checks="<<checks<<" default_GMRES_explicit_strict_SuperLU_A_only_B1_unchanged_no_solve output="<<root.string()<<'\n';
}
}
// END TYPED A ELECTRIC BACKEND INPUT CONTRACT

namespace VofRegridIntervalInputContract {
namespace fs=std::filesystem;
using VofEhdInputContract::Fixture;
using VofEhdInputContract::replace;
using VofEhdInputContract::write;
void run(const fs::path& requested){
    require(requested.is_absolute()&&requested.root_name().string()=="D:","regrid interval reader requires absolute D: output root");
    const auto root=requested/("regrid_interval_"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    require(!fs::exists(root),"regrid interval reader output must be fresh");fs::create_directories(root);
    const auto traces=[](const HydrodynamicsInput& input){
        int count=0;for(const auto& e:input.trace.entries)if(e.field=="Mesh.Regrid.interval"){
            require(e.used&&!e.defaulted&&e.units=="accepted physical steps","explicit regrid interval lost use/source/units provenance");++count;}
        return count;
    };
    Fixture omitted;const auto legacy=readVofElectrohydrodynamicsInput(omitted.save(root,"omitted").string());
    require(legacy.fluid.vofRegrid&&legacy.fluid.vofRegrid->interval==1&&traces(legacy.fluid)==0,
        "omitted interval changed original every-End behavior or old trace");
    const std::string oldMesh="Mesh { resolution (8 8); initialMesh interfaceBand; Refinement { level 1; bandCells 2; } Regrid { minLevel 0; maxLevel 1; interfaceBandLayers 2; vorticityTolerance 0.5; coarsenRatio 0.25; } }\n";
    require(legacy.fluid.documents[2].effectiveText==oldMesh,"omitted interval changed legacy effective mesh bytes");
    Fixture one;replace(one.text[2],"Regrid {","Regrid { interval 1;");
    const auto explicitOne=readVofElectrohydrodynamicsInput(one.save(root,"one").string());
    require(explicitOne.fluid.vofRegrid->interval==1&&traces(explicitOne.fluid)==1,"explicit interval1 was not consumed");
    require(explicitOne.fluid.documents.size()==legacy.fluid.documents.size(),"interval1 changed actual document inventory");
    for(std::size_t i=0;i<legacy.fluid.documents.size();++i)
        require(explicitOne.fluid.documents[i].effectiveText==legacy.fluid.documents[i].effectiveText,
            "explicit interval1 changed old effective document bytes");
    for(bool rz:{false,true}){
        Fixture supplied(rz);replace(supplied.text[2],"Regrid {","Regrid { interval 50;");
        const auto input=readVofElectrohydrodynamicsInput(supplied.save(root,rz?"rz50":"xy50").string());
        require(input.fluid.vofRegrid->interval==50&&traces(input.fluid)==1
            &&input.fluid.documents[2].effectiveText.find("interval 50;")!=std::string::npos,
            "explicit positive interval missing from typed reader/effective text");
        const auto effective=root/(rz?"rz_effective50":"xy_effective50");fs::create_directories(effective);
        for(std::size_t i=0;i<VofEhdInputContract::names.size();++i)
            write(effective/VofEhdInputContract::names[i],input.fluid.documents[i].effectiveText);
        const auto again=readVofElectrohydrodynamicsInput(effective.string());
        require(again.fluid.vofRegrid->interval==50&&traces(again.fluid)==1&&again.fluid.radial==rz
            &&again.chargeDensityJumpTolerance==input.chargeDensityJumpTolerance
            &&again.materials.epsilonLiquid==input.materials.epsilonLiquid,
            "effective interval50 roundtrip changed actual policy/material/metric");
    }
    const auto overflow=std::to_string(static_cast<long long>(std::numeric_limits<int>::max())+1);
    for(const auto& value:std::vector<std::string>{"0","-1","1.5",overflow}){
        Fixture bad;replace(bad.text[2],"Regrid {","Regrid { interval "+value+";");
        const auto path=bad.save(root,"invalid_"+value);bool caught=false;
        try{(void)readVofElectrohydrodynamicsInput(path.string());}
        catch(const std::exception& error){const std::string reason=error.what();
            caught=reason.find("mesh.dict")!=std::string::npos&&reason.find("Mesh.Regrid.interval")!=std::string::npos;}
        require(caught,"invalid regrid interval accepted or actual dictionary field provenance lost");
    }
    std::cout<<"VOF_REGRID_INTERVAL_INPUT_PASS default1_explicit1_legacy_effective_bytes XY_RZ_interval50_roundtrip guards=4 output="<<root.string()<<'\n';
}
}

int main(int argc,char** argv) {
    if(argc==3&&std::string(argv[1])=="--regrid-interval-only"){try{VofRegridIntervalInputContract::run(argv[2]);return 0;}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
    if(argc==3&&std::string(argv[1])=="--vof-ehd-field-backend-input-only"){try{VofEhdBackendInputContract::run(argv[2]);return 0;}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
    if(argc==3&&std::string(argv[1])=="--electrospray-initial-target-band-input-only"){try{ElectrosprayHfPeakContract::run(argv[2],true,true,true);return 0;}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
    if(argc==3&&std::string(argv[1])=="--electrospray-initial-target-band-only"){try{ElectrosprayHfPeakContract::run(argv[2],true,true);return 0;}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
    if(argc==3&&std::string(argv[1])=="--electrospray-hf-peaks-target-band-only"){try{ElectrosprayHfPeakContract::run(argv[2],true,false);return 0;}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
    if(argc==3&&std::string(argv[1])=="--electrospray-hf-peaks-only"){try{ElectrosprayHfPeakContract::run(argv[2]);return 0;}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}

    if(argc==3&&std::string(argv[1])=="--electrospray-profile-input-only"){try{ElectrosprayProfileInputContract::run(argv[2]);return 0;}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}

    if(argc==3&&std::string(argv[1])=="--vof-csg-target-amr-only"){try{VofCsgTargetAmrContract::run(argv[2]);return 0;}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}

if(argc==3&&std::string(argv[1])=="--vof-csg-fluid-boundary-only"){try{VofCsgFluidBoundaryContract::run(argv[2]);return 0;}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
    if(argc==3&&std::string(argv[1])=="--vof-csg-electric-boundary-only"){try{VofCsgElectricBoundaryContract::run(argv[2]);return 0;}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
    if(argc==3&&std::string(argv[1])=="--vof-csg-initial-only"){try{VofCsgInitialAppContract::run(argv[2]);return 0;}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
    if(argc==3&&std::string(argv[1])=="--vof-pressure-driven-reader-only"){try{VofPressureDrivenReservoirContract::reader(argv[2]);return 0;}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
    if(argc==3&&std::string(argv[1])=="--vof-pressure-driven-app-only"){try{VofPressureDrivenReservoirContract::app(argv[2]);return 0;}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
    if(argc==3&&std::string(argv[1])=="--vof-charge-regrid-only"){try{VofChargeRegridContract::run(argv[2]);return 0;}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
    if(argc==3&&std::string(argv[1])=="--vof-charge-advection-input-only"){try{VofChargeAdvectionInputContract::run(argv[2]);return 0;}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
    if(argc==3&&std::string(argv[1])=="--vof-ehd-surface-projection-initial-only"){try{F26SurfaceProjectionInitialContract::run(argv[2]);return 0;}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
    if(argc==4&&std::string(argv[1])=="--vof-b1-topology-dump-only"){try{VofB1TopologyDiagnostic::run(argv[2],argv[3]);return 0;}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
    if(argc==4&&std::string(argv[1])=="--vof-b1-zero-wiring-only"){try{VofB1ZeroWiringContract::run(argv[2],argv[3]);return 0;}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
    if(argc==3&&std::string(argv[1])=="--vof-ehd-pressure-backflow-charge-only"){try{F22PressureBackflowChargeContract::run(argv[2]);return 0;}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
    if(argc==3&&std::string(argv[1])=="--vof-ehd-backflow-input-only"){try{F19BackflowInputContract::run(argv[2]);return 0;}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
    if(argc==3&&std::string(argv[1])=="--vof-ehd-south-even90-input-only"){try{VofEhdOpenInputContract::SouthEven90::run(argv[2]);return 0;}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
    if(argc==3&&std::string(argv[1])=="--vof-ehd-open-joint-only"){try{VofEhdOpenJointContract::run(argv[2]);return 0;}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
    if(argc==3&&std::string(argv[1])=="--vof-ehd-open-input-only"){try{VofEhdOpenInputContract::run(argv[2]);return 0;}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
    if(argc==4&&std::string(argv[1])=="--vof-ehd-joint-only"){try{VofEhdJointContract::run(argv[2],argv[3]);return 0;}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
    if(argc==3&&std::string(argv[1])=="--vof-ehd-input-only"){try{VofEhdInputContract::run(argv[2]);return 0;}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
    try { checkFixedSteps(); checkAccumulatedEndpoint(); checkAdaptiveTimeStep(); checkFieldSchedule(); check(); std::cout << "App entry facilities: fixed steps, field schedule and AppRegistry execution/lifetime passed\n"; return 0; }
    catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
