#include "Pic/PicApp.h"

#include <iostream>

#include "Pic/PicRunAdapter.h"
#include "Pic/Setup/PicRunSetup.h"
#include "Pic/IO/PicInput.h"
#include "Pic/IO/PicOutput.h"
#include "Modules/Plasma/ExplicitPIC/ExplicitPICModule.h"
#include "Modules/Plasma/ImplicitPIC/ImplicitPICModule.h"

PicApp::PicApp(const RunConfig& config) : config_(config) {}

std::string PicApp::name() const
{
    return "Pic";
}

PicAppConfig PicApp::readConfig(const std::string& runDictPath)
{
    return readPicAppConfig(runDictPath);
}

int PicApp::run()
{
    if (config_.checkCase)
    {
        std::cerr << "[Pic] --check-case is not supported by the Pic App.\n";
        return 2;
    }

    const auto selection = readConfig();
    const auto input = readPicSharedInput();
    if (!input)
        return 2;
    const auto units = preparePicUnits(*input);
    if (!units.warning.empty())
        std::cerr << "[Pic] WARNING: " << units.warning << std::endl;
    const auto geometry = readPicGeometryClosure(*input,*units.units);
    auto setup = preparePicRun(*input,units,geometry);
    writePicRunReports(setup);
    setup.plasmaCollisions = readPicRunCollisions(*input);
    if (!setup.plasmaCollisions)
        return 2;

    if (selection.integrator == PicIntegrator::Implicit)
        return PicRunAdapter<ImplicitPICModule>(setup).run();
    return PicRunAdapter<ExplicitPICModule>(setup).run();
}