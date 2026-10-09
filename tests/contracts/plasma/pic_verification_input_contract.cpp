#include "Verification/Pic/DictionaryPicVerification.h"
#include "Verification/Pic/RegressionHarness.h"
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>

namespace fs = std::filesystem;
static void require(bool value, const char* message) {
    if (!value) throw std::runtime_error(message);
}
int main() {
    const auto scratch = fs::current_path() / ("verification_input_" +
        std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    try {
        fs::create_directories(scratch);
        const auto path = scratch / "pic.dict";
        DictionaryPicVerification input;
        require(!input.readDict(path.string()), "missing optional dictionary");
        require(!input.config().selfForceTest3.enabled && input.config().regressionHarness.mode == "off",
                "verification defaults");
        std::ofstream(path) << R"(PIC
{
 dt 2 ns
 steps 7
 regressionTest yes
 regressionUnitBoxPeriodicY on
 regressionElectricFieldPeriodicY 1
}
RegressionHarness
{
 mode CHECK
 baseline "case_baseline.dict"
 report "case_report.txt"
 absTol 3e-12
 relTol 4e-10
}
RegressionCutCircle
{
 enabled true
 levelSetIndex 2
 center (.4,.6)
 radius .3
 alphaSamplingN 5
}
RegressionCutCirclePoisson
{
 enabled true
 coeff 2.5
 computeEError yes
}
RegressionCutCirclePolarizationPoisson
{
 enabled true
 kappaMode sin2d
 kappa0 2
 kappaAmp .1
 kappaKx 3
 kappaKy 4
}
RegressionEmbeddedConductor
{
 enabled true
 particle (.2,.3)
 charge .25
 sigma .02
 sampleMode HOSTCELL
}
SelfForceTest3
{
 enabled true
 subtests (uniform traverse periodicCross)
 ablations (A0 A1)
 Particle
 {
  q 1e-12 C
  m 2 kg
 }
 Neutrality
 {
  enabled false
  mode custom
 }
 GaugeFix
 {
  enabled off
 }
 Uniform
 {
  Nlist (8 12 bad)
  samplePoints ((.37,.41) (.51,.63))
 }
 Traverse
 {
  path (.2,.3) (.7,.8) 9
 }
}
Species
{
 "electron"
 {
  charge -1.602176634e-19
  mass 9.1093837015e-31
  numParticles 4
 }
}
Output
{
 Fields
 {
  every 3
  list (phi rho)
 }
}
)";
        require(input.readDict(path.string()), "verification reader");
        const auto& c = input.config(); const auto& p = input.physics();
        require(c.regressionUnitBoxSinSin && c.regressionUnitBoxSinCosPeriodicY && c.regressionElectricFieldSinCosPeriodicY,
                "PIC-level verification aliases");
        require(c.regressionHarness.mode == "check" && c.regressionHarness.baselinePath == "case_baseline.dict"
            && c.regressionHarness.reportPath == "case_report.txt" && c.regressionHarness.absTol == 3e-12
            && c.regressionHarness.relTol == 4e-10, "baseline settings");
        require(c.regressionCutCircle.enabled && c.regressionCutCircle.center.x == .4
            && c.regressionCutCircle.radius == .3 && c.regressionCutCircle.alphaSamplingN == 5, "cut metrics settings");
        require(c.regressionCutCirclePoisson.enabled && c.regressionCutCirclePoisson.coeff == 2.5
            && c.regressionCutCirclePoisson.computeEError, "manufactured Poisson settings");
        require(c.regressionCutCirclePolarizationPoisson.enabled && c.regressionCutCirclePolarizationPoisson.kappaMode == "sin2d"
            && c.regressionCutCirclePolarizationPoisson.kappaKy == 4, "polarization settings");
        require(c.regressionEmbeddedConductor.enabled && c.regressionEmbeddedConductor.sampleMode == "hostcell"
            && c.regressionEmbeddedConductor.particle.y == .3 && c.regressionEmbeddedConductor.charge == .25, "conductor settings");
        const auto& s = c.selfForceTest3;
        require(s.enabled && s.subtests.size() == 3 && s.ablations.size() == 2 && s.particle.q == 1e-12 && s.particle.m == 2,
                "self force lists and SI units");
        require(!s.neutrality.enabled && s.neutrality.mode == "custom" && !s.gaugeFix.enabled
            && s.uniform.Nlist == std::vector<int>({8,12}) && s.uniform.samplePoints.size() == 2
            && s.traverse.start.x == .2 && s.traverse.end.y == .8 && s.traverse.numSamples == 9, "self force nested records");
        require(std::abs(p.dt-2e-9) < 1e-23 && p.steps == 7 && p.species.size() == 1
            && p.species.front().numParticles == 4 && p.output.every == 3, "physics reader continues after App records");
        DictionaryPIC ordinary;
        require(ordinary.readDict(path.string()) && ordinary.config().steps == p.steps
            && ordinary.config().species.front().numParticles == 4, "ordinary PIC requires no App reader");
        auto harness = c.regressionHarness;
        harness.baselinePath = (scratch / "baseline.dict").string();
        harness.reportPath = (scratch / "report.txt").string();
        harness.mode = "record";
        require(RegressionHarness(harness).process({{"metric",2.0}}) == 0, "baseline recording");
        harness.mode = "check";
        require(RegressionHarness(harness).process({{"metric",2.0}}) == 0, "matching baseline");
        require(RegressionHarness(harness).process({{"metric",3.0}}) != 0, "mismatch must fail");
        require(RegressionHarness(harness).process({{"different",2.0}}) != 0, "missing metric must fail");
        fs::remove(scratch / "baseline.dict");
        require(RegressionHarness(harness).process({{"metric",2.0}}) != 0, "missing baseline must fail");
        fs::remove(path); fs::remove(scratch / "report.txt"); fs::remove(scratch);
        std::cout << "PIC verification input ownership and baseline checks passed\n";
        return 0;
    } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
