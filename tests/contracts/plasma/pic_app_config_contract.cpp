#include "Pic/PicApp.h"
#include "Pic/IO/PicInput.h"
#include "Pic/IO/PicOutput.h"
#include "Pic/Setup/PicRunSetup.h"
#include "Core/Units/UnitsContext.h"
#include "Modules/Plasma/PlasmaRunInputs.h"

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

namespace fs = std::filesystem;

static void require(bool condition, const char* message)
{
    if (!condition) throw std::runtime_error(message);
}

// Test the actual responsibility boundary in a private working directory.
struct ScopedCaseDirectory
{
    const fs::path previous = fs::current_path();
    explicit ScopedCaseDirectory(const fs::path& path) { fs::current_path(path); }
    ~ScopedCaseDirectory() { std::error_code ignored; fs::current_path(previous, ignored); }
};

int main()
{
    const auto scratch = fs::current_path() / ("pic_app_config_" +
        std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    try
    {
        fs::create_directories(scratch);
        const auto path = scratch / "run.dict";
        std::ofstream(path) << "Run { app Pic; }\n";
        require(PicApp::readConfig(path.string()).integrator == PicIntegrator::Explicit,
            "default PIC integrator");

        std::ofstream(path) << "Run { app Pic; }\nPic { integrator implicit; }\nParallel { threads 4; }\n";
        require(PicApp::readConfig(path.string()).integrator == PicIntegrator::Implicit,
            "implicit integrator selection");

        auto rejects = [&](const std::string& contents) {
            std::ofstream(path) << contents;
            try
            {
                (void)PicApp::readConfig(path.string());
                return false;
            }
            catch (const std::exception&)
            {
                return true;
            }
        };
        require(rejects("Run { app Pic; }\nPic { integrator unknown; }\n"), "unknown integrator rejected");
        require(rejects("Run { app Pic; }\nPic { integrator explicit; modules ( ExplicitPIC ); }\n"),
            "generic module list rejected");
        require(rejects("Run { app Pic; }\nPic { integrator explicit; integrator implicit; }\n"),
            "duplicate integrator rejected");
        require(rejects("Run { app Pic; }\nSingleModule { modules ( ExplicitPIC ); }\n"),
            "retired SingleModule block rejected");
        require(rejects("Run { app Electric; }\n"), "another App selection rejected");
        require(rejects("Pic { integrator explicit; }\n"), "missing Run block rejected");

        // The public App helper forwards to the real IO selection service.
        std::ofstream(path) << "Run { app Pic; }\nPic { integrator explicit; }\n";
        require(readPicAppConfig(path.string()).integrator == PicApp::readConfig(path.string()).integrator,
            "App selection delegates to Pic IO");
        {
            ScopedCaseDirectory directory(scratch);
            PicSharedInput parsed;
            parsed.units = std::make_shared<UnitsContext>(UnitsContext::MakeCodeUnits());
            parsed.normalization.config.strategy = NormalizationStrategy::kManual;
            parsed.normalizationSource = "default";
            parsed.writeSIOutputs = true;
            const auto units = preparePicUnits(parsed);
            PicGeometryClosure geometry;
            geometry.hasDimension = true; geometry.dimension = 3;
            geometry.hasSlabThickness = true; geometry.slabThickness_m = .125;
            const auto prepared = preparePicRun(parsed, units, geometry);
            require(prepared.units == units.units && prepared.normalization->units() == units.units,
                "typed Setup shares its actual Units owner");
            require(prepared.normalization->config().dimension == 3
                && prepared.normalization->config().slabThickness_m == .125,
                "typed Setup applies supplied geometry closure");
            require(prepared.normalizationExplicit && prepared.normalizationSource == "geometry.dict/Domain"
                && prepared.dimensionSource == "geometry.dict/Domain" && prepared.slabSource == "geometry.dict/Domain",
                "typed Setup retains geometry provenance");
            const auto bound = bindPicRunInputs(prepared);
            require(bound.units == prepared.units && bound.normalization == prepared.normalization && bound.writeSIOutputs,
                "runtime binding uses the prepared contexts");
            require(!fs::exists("output") && !fs::exists("cache"),
                "typed Setup performs no directory or report IO");
            writePicRunReports(prepared);
            require(fs::is_regular_file("output/meta/units_report.txt")
                && fs::is_regular_file("output/meta/normalization_report.txt") && !fs::exists("cache"),
                "Pic Output owns the actual reports, without cache input IO");
            fs::remove("output/meta/units_report.txt");
            fs::remove("output/meta/normalization_report.txt");
            fs::remove("output/meta"); fs::remove("output");
        }

        fs::remove(path);
        fs::remove(scratch);
        std::cout << "PIC App integrator input contract passed\n";
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        std::error_code ignored;
        fs::remove_all(scratch, ignored);
        return 1;
    }
}