#include "Modules/Plasma/IO/PlasmaCollisionLoader.h"
#include "Modules/Plasma/IO/DictionaryGas.h"
#include "Modules/Plasma/IO/DictionaryCrossSections.h"
#include "Modules/Plasma/IO/DictionaryReactions.h"
#include "Modules/Plasma/IO/DictionaryCollisions.h"
#include "Core/IO/DictionaryBlockFile.h"
#include "Core/IO/DictionaryBlockText.h"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

namespace
{
void writeGasReportOrWarn(const std::string& outDir, const GasState& gas)
{
    namespace fs = std::filesystem;
    std::error_code ec;
    fs::create_directories(outDir, ec);

    const std::string filePath = (fs::path(outDir) / "gas_report.txt").string();
    std::ofstream out(filePath, std::ios::out);
    if (!out)
    {
        std::cerr << "[RunController] WARNING: failed to write gas report: " << filePath << std::endl;
        return;
    }

    out << "BackgroundGas\n";
    out << "  name        : " << gas.name << "\n";
    out << "  pressure_Pa : " << std::scientific << std::setprecision(10) << gas.pressure_Pa << "\n";
    out << "  temperature_K : " << std::scientific << std::setprecision(10) << gas.temperature_K << "\n";
    out << "  n_gas_1pm3  : " << std::scientific << std::setprecision(10) << gas.number_density_m3 << "\n";
    out << "  n_gas_1pcm3 : " << std::scientific << std::setprecision(10) << (gas.number_density_m3 * 1e-6) << "\n";
    if (gas.particle_mass_kg > 0.0)
    {
        out << "  particle_mass_kg : " << std::scientific << std::setprecision(10) << gas.particle_mass_kg << "\n";
        out << "  rho_kgpm3        : " << std::scientific << std::setprecision(10) << gas.mass_density_kgpm3 << "\n";
        out << "  vth_mps          : " << std::scientific << std::setprecision(10) << gas.thermal_speed_mps << "\n";
    }
}

void writeCrossSectionsReportOrWarn(const std::string& outDir, const CrossSectionDB& db)
{
    namespace fs = std::filesystem;
    std::error_code ec;
    fs::create_directories(outDir, ec);

    const std::string filePath = (fs::path(outDir) / "crossSections_report.txt").string();
    std::ofstream out(filePath, std::ios::out);
    if (!out)
    {
        std::cerr << "[RunController] WARNING: failed to write cross-sections report: " << filePath << std::endl;
        return;
    }
    out << db.TableString();
    out << "\n";
}

void writeReactionsReportOrWarn(const std::string& outDir, const ReactionDatabase& db)
{
    namespace fs = std::filesystem;
    std::error_code ec;
    fs::create_directories(outDir, ec);

    const std::string filePath = (fs::path(outDir) / "reactions_report.txt").string();
    std::ofstream out(filePath, std::ios::out);
    if (!out)
    {
        std::cerr << "[RunController] WARNING: failed to write reactions report: " << filePath << std::endl;
        return;
    }
    out << db.TableString();
    out << "\n";
}

} // namespace

std::shared_ptr<const PlasmaCollisionInputs> loadPlasmaCollisionInputs(const std::string& picPath)
{
    namespace fs = std::filesystem;
    DictionaryBlockText parser;
    std::string picText, mccText;
    if (!picPath.empty() && fs::exists(picPath) && !parser.loadFile(picPath, picText))
    {
        std::cerr << "[Plasma] ERROR: cannot read " << picPath << '\n';
        return nullptr;
    }
    if (fs::exists("./system/mcc.dict"))
        parser.loadFile("./system/mcc.dict", mccText);

    // Each block has its own fallback, so mixed legacy sources remain valid.
    auto select = [&](const std::string& label, const std::string& name,
                      const std::vector<std::string>& mccNames,
                      const std::vector<std::string>& picNames, std::string& path) {
        std::string block;
        if (parser.extractFirstBlockText(mccText, mccNames, block, DictionaryBlockText::Scope::TopLevel) ||
            parser.extractFirstBlockText(picText, picNames, block, DictionaryBlockText::Scope::TopLevel))
        {
            path = (fs::path("./cache") / ("pic_" + label + ".dict")).string();
            if (!writeDictionaryBlockFile(path, name, block))
            {
                std::cerr << "[Plasma] ERROR: failed to write " << name << " block to " << path << '\n';
                return false;
            }
        }
        else
            path = "./system/" + label + ".dict";
        return true;
    };
    auto inputs = std::make_shared<PlasmaCollisionInputs>();
    std::string path;
    if (!select("gas", "Gas", {"Gas", "gas", "BackgroundGas", "backgroundGas"}, {"Gas", "gas"}, path))
        return nullptr;
    DictionaryGas gas;
    if (gas.readDict(path)) inputs->gas = std::make_shared<GasState>(gas.state());
    if (inputs->gas)
    {
        const auto& g = *inputs->gas;
        std::cout << "[RunController] BackgroundGas: name=" << g.name
                  << " p=" << std::scientific << std::setprecision(6) << g.pressure_Pa << " Pa"
                  << " T=" << std::scientific << std::setprecision(6) << g.temperature_K << " K"
                  << " n_gas=" << std::scientific << std::setprecision(6) << g.number_density_m3 << " 1/m^3"
                  << " (" << std::scientific << std::setprecision(6) << (g.number_density_m3 * 1e-6) << " 1/cm^3)"
                  << std::endl;
        writeGasReportOrWarn("./output/meta", g);
    }
    else
    {
        std::cout << "[RunController] BackgroundGas: (none)" << std::endl;
    }

    if (!select("crossSections", "CrossSections", {"CrossSections", "crossSections", "Crosssections", "crosssections"}, {"CrossSections", "crossSections", "Crosssections", "crosssections"}, path))
        return nullptr;
    if (fs::exists(path))
    {
        DictionaryCrossSections dictXS;
        if (!dictXS.readDict(path))
        {
            std::cerr << "[RunController] ERROR: failed to parse cross sections dictionary" << std::endl;
            return nullptr;
        }
        inputs->crossSections = std::make_shared<CrossSectionDB>(dictXS.db());
        std::cout << "[RunController] CrossSections\n";
        std::cout << inputs->crossSections->TableString() << std::endl;
        writeCrossSectionsReportOrWarn("./output/meta", *inputs->crossSections);
    }
    else
    {
        std::cout << "[RunController] CrossSections: (none)" << std::endl;
    }
    if (!select("reactions", "Reactions", {"Reactions", "reactions"}, {"Reactions", "reactions"}, path))
        return nullptr;
    if (fs::exists(path))
    {
        if (!inputs->crossSections)
        {
            std::cerr << "[RunController] ERROR: reactions.dict provided but crossSections.dict is missing" << std::endl;
            return nullptr;
        }
        DictionaryReactions dictR;
        if (!dictR.readDict(path, inputs->crossSections.get()))
        {
            std::cerr << "[RunController] ERROR: failed to parse reactions dictionary" << std::endl;
            return nullptr;
        }
        inputs->reactions = std::make_shared<ReactionDatabase>(dictR.db());
        std::cout << "[RunController] Reactions\n";
        std::cout << inputs->reactions->TableString() << std::endl;
        writeReactionsReportOrWarn("./output/meta", *inputs->reactions);
    }
    else
    {
        std::cout << "[RunController] Reactions: (none)" << std::endl;
    }
    if (!select("collisions", "Collisions", {"Collisions", "collisions"}, {"Collisions", "collisions"}, path))
        return nullptr;
    if (fs::exists(path))
    {
        if (!inputs->gas)
        {
            std::cerr << "[RunController] ERROR: collisions.dict provided but gas.dict is missing" << std::endl;
            return nullptr;
        }
        if (!inputs->crossSections)
        {
            std::cerr << "[RunController] ERROR: collisions.dict provided but crossSections.dict is missing" << std::endl;
            return nullptr;
        }
        if (!inputs->reactions)
        {
            std::cerr << "[RunController] ERROR: collisions.dict provided but reactions.dict is missing" << std::endl;
            return nullptr;
        }
        DictionaryCollisions dictC;
        if (!dictC.readDict(path))
        {
            std::cerr << "[RunController] ERROR: failed to parse collisions dictionary" << std::endl;
            return nullptr;
        }
        inputs->collisions = std::make_shared<CollisionsConfig>(dictC.config());
        std::cout << "[RunController] Collisions\n";
        std::cout << inputs->collisions->TableString() << std::endl;
    }
    else
    {
        std::cout << "[RunController] Collisions: (none)" << std::endl;
    }

    return inputs;
}
