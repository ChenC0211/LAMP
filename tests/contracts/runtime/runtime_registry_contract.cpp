#include "Runtime/AppRegistry.h"
#include "Runtime/RunControls.h"

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
        require(config.checkCase && config.startMode == RunStartMode::Resume && config.schemaVersion == 2,
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
static void check() {
    Capture capture;
    AppRegistry registry;
    RunConfig cfg;
    std::vector<std::string> events;
    require(registry.run(cfg) == 1 && capture.out.str().empty() && capture.err.str().empty(), "empty selection");
    cfg.app = "missing";
    require(registry.run(cfg) == 2 && capture.err.str().find("Unknown app: 'missing'") != std::string::npos,
        "unknown selection");
    cfg.checkCase = true; cfg.startMode = RunStartMode::Resume; cfg.schemaVersion = 2;
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
    require(cfg.app == "aLpHa" && cfg.checkCase && cfg.schemaVersion == 2, "caller configuration changed");
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
int main() {
    try { checkFixedSteps(); checkFieldSchedule(); check(); std::cout << "Independent Runtime: fixed steps, field schedule and AppRegistry execution/lifetime passed\n"; return 0; }
    catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
