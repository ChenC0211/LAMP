#include "Modules/Plasma/IO/DictionaryMotionSchedule.h"
#include "Modules/Plasma/Boundary/LegacyGridMotionEvolution.h"
#include "Core/IO/DictionaryGeometry.h"
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>

namespace {
void check(bool value,const char* why) { if (!value) throw std::runtime_error(why); }
void near(double actual,double expected,const char* why)
{ check(std::isfinite(actual) && std::abs(actual-expected)<1e-14,why); }
struct DirectoryScope {
    std::filesystem::path previous=std::filesystem::current_path();
    DirectoryScope() {
        const auto path=previous/("grid_motion_"+std::to_string(std::chrono::high_resolution_clock::now().time_since_epoch().count()));
        std::filesystem::create_directories(path/"system"); std::filesystem::current_path(path);
    }
    ~DirectoryScope() { std::filesystem::current_path(previous); }
};
}
int main() try
{
    DirectoryScope directory;
    {
        std::ofstream input("system/geometry.dict");
        input<<R"(MotionSchedule
{
 updateEvery 3
 outputEvery 2
 metricsPath "grid.csv"
 startupHoldSteps 1
 stage1Steps 2
 stage2Steps 2
 gsaInitial 0.7 mm
 gsaHot 0.72 mm
 gadInitial 1 mm
 gadHot 0.97 mm
 daInitial 0.3 mm
 daFinal 0.32 mm
 accelShiftStage1 0.02 mm
 decelShiftStage1 -0.01 mm
 accelRecessionStage2 0.01 mm
 accelRegionXMin 1.5 mm
 accelRegionXMax 1.8 mm
 decelRegionXMin 2.8 mm
 decelRegionXMax 3.1 mm
 apertureHalpGapInitial 0.15 mm
 accelDownstreamFaceX 1.8 mm
 recessionBandXMin 1.5 mm
 recessionBandXMax 1.8 mm
 recessionBandYMin 0.15 mm
 recessionBandYMax 0.8 mm
}
Domain
{
 min (0, 0)
 max (1, 1)
}
)";
    }
    DictionaryGeometry geometry; geometry.readDict();
    auto config=readMotionSchedule(geometry);
    check(config.enabled && geometry.extensionEnabled("MotionSchedule") && config.updateEvery==3
        && config.outputEvery==2 && config.metricsPath=="grid.csv","legacy defaults/cadence/output");
    near(config.apertureHalfGapInitial_m,.00015,"legacy misspelled aperture alias");
    near(config.accelRegionXMin_m,.0015,"accel selector SI"); near(config.accelRegionXMax_m,.0018,"accel selector SI");
    near(config.decelRegionXMin_m,.0028,"decel selector SI"); near(config.decelRegionXMax_m,.0031,"decel selector SI");
    std::vector<Line> reference(6);
    for (int i=0;i<6;++i)
    {
        const double left=i<2?.0005:(i<4?.0015:.0028);
        const double y=(i%2==0)? .00015:-.00015;
        reference[i].points[0]={left,y}; reference[i].points[1]={left+.0003,y};
        reference[i].setType(LineType::straight);
    }
    LegacyGridMotion::Bindings bindings{{0},{1},{2},{3},{4},{5}};
    const LegacyGridMotion motion(config,bindings);
    bindings = {};
    config.accelShiftStage1_m=99.; // Construction owns the original config/binding.
    auto hold=motion.state(0,true); check(hold.stageName=="stage0_hold","hold stage"); near(hold.accelShift_m,0.,"hold position");
    const auto half=motion.state(1,true); check(half.stageName=="stage1_hot_gap","gap stage");
    near(half.stage1Alpha,.5,"first gap fraction"); near(half.accelShift_m,.00001,"half displacement");
    near(motion.state(2,true).stage2Alpha,0.,"last gap step has no recession");
    near(motion.state(3,true).stage2Alpha,.5,"first recession fraction");
    auto final=motion.state(4,true); check(final.stageName=="stage2_recession","recession stage");
    near(final.stage2Alpha,1.,"final recession fraction"); near(motion.state(99,true).stage2Alpha,1.,"saturated schedule");
    check(motion.state(3,false).stageName=="static","inactive workflow");
    auto lines=reference;
    motion.deform(lines,final);
    near(lines[0].points[1].x,.0008,"screen fixed");
    near(lines[2].points[0].x,.00152,"accel upstream translation");
    near(lines[2].points[1].x,.00183,"downstream face recession");
    near(lines[2].points[1].y,.00016,"upper recession direction");
    near(lines[3].points[1].y,-.00016,"lower recession direction");
    near(lines[4].points[0].x,.00279,"decel translation");
    motion.measure(lines,final);
    near(final.gsa_m,.00072,"measured screen-accel gap"); near(final.gad_m,.00096,"measured accel-decel gap");
    near(final.da_m,.0003,"measured aperture includes unchanged upstream endpoint");
    // Restoring the same reference must reproduce the same absolute geometry.
    auto repeated=reference; motion.deform(repeated,motion.state(4,true));
    for (int i=0;i<6;++i) for (int j=0;j<2;++j)
    { near(repeated[i].points[j].x,lines[i].points[j].x,"absolute x repeat"); near(repeated[i].points[j].y,lines[i].points[j].y,"absolute y repeat"); }
    auto zero=readMotionSchedule(geometry); zero.startupHoldSteps=-1; zero.stage1Steps=zero.stage2Steps=0;
    LegacyGridMotion instantaneous(zero,{});
    near(instantaneous.state(0,true).stage1Alpha,1.,"zero-duration legacy gap");
    near(instantaneous.state(0,true).stage2Alpha,1.,"zero-duration legacy recession");
    // Recover every accepted state, including skipped cadence hits and stage transitions.
    for (int every : {1,3,7}) for (int holdSteps : {0,1,3})
        for (int gapSteps : {0,2}) for (int recessionSteps : {0,2})
    {
        auto settings=readMotionSchedule(geometry);
        settings.updateEvery=every; settings.startupHoldSteps=holdSteps;
        settings.stage1Steps=gapSteps; settings.stage2Steps=recessionSteps;
        const LegacyGridMotion::Bindings lineIds{{0},{1},{2},{3},{4},{5}};
        LegacyGridMotionEvolution::Reference initial{reference,{{0,0},{1,0},{2,0},{3,0},{4,0},{5,0}}};
        auto evolution=LegacyGridMotionEvolution(settings,lineIds,initial,{{2,3},{4,5}});
        initial.lines[2].points[0].x=99.; initial.innerPoints[2].x=99.;
        auto live=reference;
        auto candidate=evolution.prepare(0,true);
        evolution.deformReference(live,candidate); evolution.accept(candidate,live);
        for (int step=0;step<12;++step)
        {
            candidate=evolution.prepare(step);
            if(candidate.changed) { evolution.deformReference(live,candidate); evolution.accept(candidate,live); }
            LegacyGridMotionEvolution restored(settings,lineIds,{reference,{{0,0},{1,0},{2,0},{3,0},{4,0},{5,0}}},{{2,3},{4,5}});
            restored.restore(step+1,live);
            const auto& actual=restored.current(); const auto& expected=evolution.current();
            check(actual.stageName==expected.stageName,"restored accepted stage");
            near(actual.stage1Alpha,expected.stage1Alpha,"restored gap progress");
            near(actual.stage2Alpha,expected.stage2Alpha,"restored recession progress");
            near(actual.gsa_m,expected.gsa_m,"restored measured gap");
            const auto next=restored.prepare(step+1),fullNext=evolution.prepare(step+1);
            check(next.changed==fullNext.changed,"resumed next update cadence");
            near(next.stage1Alpha,fullNext.stage1Alpha,"resumed next gap progress");
            near(next.stage2Alpha,fullNext.stage2Alpha,"resumed next recession progress");
            near(evolution.innerPoint(2,expected).x,2+expected.accelShift_m,"owned initial inner point");
            auto incompatible=live; incompatible[2].points[0].x+=.001;
            const double savedGap=actual.gsa_m;
            bool rejected=false;
            try { restored.restore(step+1,incompatible); } catch(const std::runtime_error&) { rejected=true; }
            check(rejected,"incompatible reference checkpoint rejected");
            near(restored.current().gsa_m,savedGap,"failed restore preserves accepted state");
        }
    }
    { std::ofstream input("system/geometry.dict"); input<<"MotionSchedule\n{\n enabled false\n}\n"; }
    geometry.readDict(); check(!readMotionSchedule(geometry).enabled && !geometry.extensionEnabled("MotionSchedule"),"disabled input");
    { std::ofstream input("system/geometry.dict"); input<<"Domain\n{\n min (0, 0)\n max (1, 1)\n}\n"; }
    geometry.readDict(); check(!readMotionSchedule(geometry).enabled,"missing input defaults");
    std::cout<<"Legacy grid schedule, SI input, binding ownership, deformation and measurement checks passed\n";
    return 0;
}
catch (const std::exception& error) { std::cerr<<error.what()<<'\n'; return 1; }
