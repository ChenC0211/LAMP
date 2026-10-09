#include "Modules/Plasma/IO/DictionaryMovingEB.h"
#include "Modules/Plasma/Boundary/MovingEBController.h"
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>

namespace {
void check(bool valid,const char* reason) { if (!valid) throw std::runtime_error(reason); }
void near(double actual,double expected,const char* reason)
{ check(std::isfinite(actual) && std::abs(actual-expected)<1e-12,reason); }
struct DirectoryScope {
    std::filesystem::path previous=std::filesystem::current_path();
    DirectoryScope() {
        const auto path=previous/("movingeb_input_"+std::to_string(std::chrono::high_resolution_clock::now().time_since_epoch().count()));
        std::filesystem::create_directories(path/"system"); std::filesystem::current_path(path);
    }
    ~DirectoryScope() { std::filesystem::current_path(previous); }
};
const char* domain="Domain\n{\n min (0, 0)\n max (1, 1)\n dimension 2\n}\n";
void input(const std::string& body) { std::ofstream out("system/geometry.dict"); out<<body<<domain; }
}
int main() try
{
    DirectoryScope directory;
    DictionaryGeometry geometry;
    input(""); geometry.readDict();
    check(!geometry.extensionEnabled("MovingEB") && !readMovingEBConfig(geometry).enabled,"absent input defaults");
    input(R"(
MovingEB
{
 updateEvery 2
 outputEvery 3
 restart off
 restartPath "wall_cache.bin"
 metricsPath "wall_metrics.csv"
 ImpactModel
 {
  ionsOnly no
  minSamples 3
  smoothing exponential
  smoothingLength 2 mm
 }
 Region
 {
  name "wall {lower"
  boundaryId 7
  levelSetIndex 2
  selector lineIds (1, 3)
  xmin 2 cm
  xmax 30 mm
  ymin -1 mm
  ymax 1 mm
  motionMode normalErosion
  model ionErosion
  startStep 2
  endStep 4
  velocityX 1.25
  velocityY -2.5
  normalSpeed 3.75
  c0 4.5
  alpha 0.5
  energyModel power
  energyRef 25
  angleModel cosine
  angleExponent 2
  smoothingLength 0.5 mm
  vmax 10
 }
}
)");
    geometry.readDict();
    auto config=readMovingEBConfig(geometry);
    check(config.enabled && geometry.extensionEnabled("MovingEB") && config.updateEvery==2
        && config.outputEvery==3 && !config.restart && config.restartPath=="wall_cache.bin"
        && config.metricsPath=="wall_metrics.csv","root settings and default enabled");
    check(!config.impactModel.ionsOnly && config.impactModel.minSamples==3
        && config.impactModel.smoothing=="exponential","impact model");
    near(config.impactModel.smoothingLength_m,.002,"SI smoothing length");
    check(config.regions.size()==1,"one complete region"); const auto& region=config.regions[0];
    check(region.name=="wall {lower" && region.boundaryId==7 && region.levelSetIndex==1
        && region.selector.lineIds==std::vector<int>({0,2}) && !region.selector.wholeBoundary,"names and legacy indices");
    check(region.selector.hasXRange && region.selector.hasYRange && region.startStep==2 && region.endStep==4,"selector and active window");
    near(region.selector.xMin_m,.02,"xmin units"); near(region.selector.xMax_m,.03,"xmax units");
    near(region.selector.yMin_m,-.001,"ymin units"); near(region.selector.yMax_m,.001,"ymax units");
    near(region.velocityX_mps,1.25,"velocity x"); near(region.velocityY_mps,-2.5,"velocity y");
    near(region.normalSpeed_mps,3.75,"normal speed"); near(region.C0,4.5,"erosion coefficient");
    near(region.alpha,.5,"flux exponent"); near(region.energyRef_eV,25.,"reference energy");
    near(region.angleExponent,2.,"angle exponent"); near(region.smoothingLength_m,.0005,"region smoothing");
    near(region.vmax_mps,10.,"speed limit");
    check(region.motionMode=="normalErosion" && region.model=="ionErosion"
        && region.energyModel=="power" && region.angleModel=="cosine","model selection");
    MovingEBController controller(geometry,config);
    config.enabled=false; config.regions.clear();
    check(controller.enabled() && controller.activeRegionDiagnostics(3).size()==1
        && controller.activeRegionDiagnostics(3)[0].name=="wall {lower"
        && controller.activeRegionDiagnostics(3)[0].normalErosion
        && controller.activeRegionDiagnostics(1).empty(),"controller owns immutable input and active window");

    input("MovingEB\n{\n enabled false\n Region\n{\n enabled true\n}\n}\n"); geometry.readDict();
    check(!readMovingEBConfig(geometry).enabled && !geometry.extensionEnabled("MovingEB"),"nested enable cannot override disabled section");
    input("MovingEB\n{\n enabled false\n}\nMovingEB\n{\n outputEvery 9\n}\n"); geometry.readDict();
    check(readMovingEBConfig(geometry).enabled && geometry.extensionEnabled("MovingEB")
        && readMovingEBConfig(geometry).outputEvery==9,"legacy repeated-section order");
    input(""); geometry.readDict();
    check(!geometry.extensionBlock("MovingEB") && !readMovingEBConfig(geometry).enabled,"re-reading clears extensions");
    std::cout<<"MovingEB schema, SI values, module ownership and capability checks passed\n";
    return 0;
}
catch (const std::exception& error) { std::cerr<<error.what()<<'\n'; return 1; }
