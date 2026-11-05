#include "Pic/IO/PicOutput.h"
#include "Pic/Setup/PicRunSetup.h"
#include "Core/Units/UnitsContext.h"
#include "Core/Units/NormalizationContext.h"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

namespace
{
void writeUnitsReportOrWarn(const std::string& outDir,
                            const UnitsContext& units,
                            bool writeSIOutputs)
{
    namespace fs = std::filesystem;
    std::error_code ec;
    fs::create_directories(outDir, ec);

    const std::string filePath = (fs::path(outDir) / "units_report.txt").string();
    std::ofstream out(filePath, std::ios::out);
    if (!out)
    {
        std::cerr << "[Pic] WARNING: failed to write units report: " << filePath << std::endl;
        return;
    }

    out << units.TableString();
    out << "\n";
    out << "OutputOptions\n";
    out << "  writeSIOutputs : " << (writeSIOutputs ? "true" : "false") << "\n";
    out << "\n";
    out << "Conversions\n";
    out << "  1 code length   = " << std::scientific << std::setprecision(8) << units.scales().L0_m << " m\n";
    out << "  1 code time     = " << std::scientific << std::setprecision(8) << units.scales().T0_s << " s\n";
    out << "  1 code potential= " << std::scientific << std::setprecision(8) << units.scales().Phi0_V << " V\n";
    out << "  1 code density  = " << std::scientific << std::setprecision(8) << units.scales().n0_m3 << " 1/m^3\n";
    out << "  1 code velocity = " << std::scientific << std::setprecision(8) << units.scales().v0_mps << " m/s\n";
    out << "  1 code field    = " << std::scientific << std::setprecision(8) << units.scales().E0_Vpm << " V/m\n";
    out << "  1 code rho      = " << std::scientific << std::setprecision(8) << units.scales().rho0_Cpm3 << " C/m^3\n";
}

void writeNormalizationReportOrWarn(const std::string& outDir,
                                   const NormalizationContext& norm,
                                   bool explicitConfig,
                                   const std::string& source,
                                   const std::string& dimensionSource,
                                   const std::string& slabSource)
{
    namespace fs = std::filesystem;
    std::error_code ec;
    fs::create_directories(outDir, ec);

    const std::string filePath = (fs::path(outDir) / "normalization_report.txt").string();
    std::ofstream out(filePath, std::ios::out);
    if (!out)
    {
        std::cerr << "[Pic] WARNING: failed to write normalization report: " << filePath << std::endl;
        return;
    }

    out << norm.TableString();
    out << "\n";
    out << "NormalizationConfig\n";
    out << "  explicit : " << (explicitConfig ? "true" : "false") << "\n";
    out << "  source   : " << (source.empty() ? "(unknown)" : source) << "\n";
    out << "  dimension_source : " << (dimensionSource.empty() ? "(unknown)" : dimensionSource) << "\n";
    out << "  slabThickness_source : " << (slabSource.empty() ? "(unknown)" : slabSource) << "\n";
}

}

void writePicRunReports(const PicRunSetup& setup)
{
    const auto& normLoad = setup.normalizationInput;
    const bool normExplicit = setup.normalizationExplicit;
    const auto& normSource = setup.normalizationSource;
    const auto& normDimSource = setup.dimensionSource;
    const auto& normSlabSource = setup.slabSource;
    std::cout << setup.units->TableString();
    writeUnitsReportOrWarn("./output/meta", *setup.units, setup.writeSIOutputs);

    std::cout << setup.normalization->TableString();
    {
        auto fmtSci = [](double v, int p = 6) {
            std::ostringstream oss;
            oss << std::scientific << std::setprecision(p) << v;
            return oss.str();
        };
        std::cout << "[Pic] Normalization source: " << normSource << "\n";
        std::cout << "[Pic] Normalization closure: dimension=" << normLoad.config.dimension
                  << " (" << normDimSource << ")";
        if (normLoad.config.dimension == 2)
        {
            std::cout << " slabThickness=" << fmtSci(normLoad.config.slabThickness_m) << " m"
                      << " (" << normSlabSource << ")";
        }
        std::cout << "\n";
    }
    if (!normExplicit)
    {
        std::cerr << "[Pic] WARNING: Normalization not specified; using defaults (dimension=2, slabThickness=1 m).\n";
    }
    if (normLoad.config.dimension == 2 && normSlabSource == "default")
    {
        std::cerr << "[Pic] WARNING: slabThickness not specified for 2D; defaulting to 1 m. "
                     "Set Normalization.slabThickness or Geometry.Domain.slabThickness.\n";
    }
    writeNormalizationReportOrWarn("./output/meta", *setup.normalization, normExplicit,
                                   normSource, normDimSource, normSlabSource);

}
