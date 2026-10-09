#include "Modules/Plasma/IO/DictionaryNormalization.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include "Core/Units/UnitsContext.h"
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "Modules/Plasma/Physics/PhysicalConstants.h"
#include "Core/Units/UnitParser.h"

using namespace std;

namespace
{
string joinTokens(const vector<string>& tokens, size_t begin_index)
{
    string out;
    for (size_t i = begin_index; i < tokens.size(); ++i)
    {
        if (!out.empty())
            out += " ";
        out += tokens[i];
    }
    return out;
}

static std::string toLowerCopyLocal(std::string s)
{
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return s;
}

NormalizationStrategy parseStrategyToken(const std::string& token)
{
    const std::string t = toLowerCopyLocal(token);
    if (t == "manual")
        return NormalizationStrategy::kManual;
    if (t == "plasmafrequency" || t == "plasma" || t == "omegap" || t == "wp")
        return NormalizationStrategy::kPlasmaFrequency;
    if (t == "thermal" || t == "te" || t == "temperature")
        return NormalizationStrategy::kThermal;
    return NormalizationStrategy::kManual;
}

bool parseValuePlusUnitToSI(const std::vector<std::string>& tokens,
                           size_t value_index,
                           UnitParser::Result& out)
{
    if (tokens.size() <= value_index)
        return false;

    const string joined = joinTokens(tokens, value_index);
    if (!UnitParser::TryParseToSI(joined, out))
        return false;
    return true;
}

bool parseChargeToSI(const std::vector<std::string>& tokens,
                     size_t value_index,
                     double& out_C)
{
    UnitParser::Result r;
    if (!parseValuePlusUnitToSI(tokens, value_index, r))
        return false;

    if (r.has_unit && r.unit != "c" && r.unit != "e")
        return false;

    out_C = r.value_si;
    return true;
}

bool parseMassToSI(const std::vector<std::string>& tokens,
                   size_t value_index,
                   double& out_kg)
{
    UnitParser::Result r;
    if (!parseValuePlusUnitToSI(tokens, value_index, r))
        return false;

    if (r.has_unit && r.unit != "kg" && r.unit != "g" && r.unit != "amu" && r.unit != "u")
        return false;

    out_kg = r.value_si;
    return true;
}

struct RefSpeciesInfo
{
    double charge_C = 0.0;
    double mass_kg = 0.0;
};

RefSpeciesInfo resolveRefSpeciesInfo(const NormalizationConfig& cfg, std::string& warn)
{
    constexpr double kElementaryCharge_C = PhysicalConstants::qe_C;
    constexpr double kElectronMass_kg = PhysicalConstants::me_kg;
    constexpr double kProtonMass_kg = PhysicalConstants::mp_kg;

    auto lowerCopy = [](std::string s) -> std::string {
        std::transform(s.begin(), s.end(), s.begin(),
                       [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
        return s;
    };

    RefSpeciesInfo out;
    if (cfg.refCharge_C > 0.0 && cfg.refMass_kg > 0.0)
    {
        out.charge_C = cfg.refCharge_C;
        out.mass_kg = cfg.refMass_kg;
        return out;
    }

    std::string name = cfg.refSpecies;
    if (name.empty())
        name = "electron";
    std::string n = lowerCopy(name);

    if (n == "electron" || n == "e")
    {
        out.charge_C = -kElementaryCharge_C;
        out.mass_kg = kElectronMass_kg;
        return out;
    }
    if (n == "proton" || n == "p" || n == "h+" || n == "hplus")
    {
        out.charge_C = kElementaryCharge_C;
        out.mass_kg = kProtonMass_kg;
        return out;
    }

    warn = "unknown refSpecies '" + name + "', fallback to electron";
    out.charge_C = -kElementaryCharge_C;
    out.mass_kg = kElectronMass_kg;
    return out;
}

} // namespace

std::string DictionaryNormalization::toLowerCopy(std::string s)
{
    transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return static_cast<char>(tolower(c)); });
    return s;
}

bool DictionaryNormalization::readDict(const std::string& dictPath)
{
    ifstream file(dictPath);
    if (!file.is_open())
        return false;

    cfg_ = NormalizationConfig{};
    hasDimension_ = false;
    hasSlabThickness_ = false;
    hasStrategy_ = false;
    hasAny_ = false;

    bool inBlock = false;
    string line;

    while (getline(file, line))
    {
        line = removeInlineComments(line);
        line = trim(line);

        if (line.empty() || line == "{")
            continue;

        if (line == "Normalization" || line == "normalization")
        {
            inBlock = true;
            continue;
        }

        if (line == "}")
        {
            if (inBlock)
                inBlock = false;
            continue;
        }

        // Accept both inside Normalization{} and at top-level.
        vector<string> tokens = split(line);
        if (tokens.empty())
            continue;

        const string key = toLowerCopy(tokens[0]);

        if (key == "strategy" && tokens.size() >= 2)
        {
            cfg_.strategy = parseStrategyToken(tokens[1]);
            hasStrategy_ = true;
            hasAny_ = true;
            continue;
        }

        if (key == "dimension" && tokens.size() >= 2)
        {
            try
            {
                cfg_.dimension = stoi(tokens[1]);
                hasDimension_ = true;
                hasAny_ = true;
            }
            catch (...)
            {
                cerr << "[DictionaryNormalization] WARNING: failed to parse dimension from line: " << line << endl;
            }
            continue;
        }

        if (key == "slabthickness" && tokens.size() >= 2)
        {
            UnitParser::Result r;
            if (parseValuePlusUnitToSI(tokens, 1, r))
            {
                cfg_.slabThickness_m = r.value_si;
                hasSlabThickness_ = true;
                hasAny_ = true;
            }
            else
            {
                cerr << "[DictionaryNormalization] WARNING: failed to parse slabThickness from line: " << line << endl;
            }
            continue;
        }

        if ((key == "refspecies" || key == "ref" || key == "reference") && tokens.size() >= 2)
        {
            cfg_.refSpecies = tokens[1];
            hasAny_ = true;
            continue;
        }

        if ((key == "refcharge" || key == "refq" || key == "refcharge_c") && tokens.size() >= 2)
        {
            double q = 0.0;
            if (parseChargeToSI(tokens, 1, q))
            {
                cfg_.refCharge_C = q;
                hasAny_ = true;
            }
            else
            {
                cerr << "[DictionaryNormalization] WARNING: failed to parse refCharge from line: " << line << endl;
            }
            continue;
        }

        if ((key == "refmass" || key == "refm" || key == "refmass_kg") && tokens.size() >= 2)
        {
            double m = 0.0;
            if (parseMassToSI(tokens, 1, m))
            {
                cfg_.refMass_kg = m;
                hasAny_ = true;
            }
            else
            {
                cerr << "[DictionaryNormalization] WARNING: failed to parse refMass from line: " << line << endl;
            }
            continue;
        }

        if ((key == "te0" || key == "t0e" || key == "te") && tokens.size() >= 2)
        {
            UnitParser::Result r;
            if (parseValuePlusUnitToSI(tokens, 1, r))
            {
                // UnitParser returns SI base: eV->J, K->K.
                const string u = toLowerCopy(r.unit);
                if (u == "ev")
                {
                    cfg_.Te0_eV = r.value_si / PhysicalConstants::eV_J;
                    hasAny_ = true;
                }
                else if (u == "k")
                {
                    cfg_.Te0_eV = (PhysicalConstants::kB_JpK * r.value_si) / PhysicalConstants::eV_J;
                    hasAny_ = true;
                }
                else
                {
                    // Interpret as energy in Joules by default.
                    cfg_.Te0_eV = r.value_si / PhysicalConstants::eV_J;
                    hasAny_ = true;
                }
            }
            else
            {
                cerr << "[DictionaryNormalization] WARNING: failed to parse Te0 from line: " << line << endl;
            }
            continue;
        }

        // Silently ignore unknown keys.
        (void)inBlock;
    }

    file.close();
    return true;
}

UnitsScales resolvePlasmaReferenceScales(const UnitsContext& baseUnits,
                                       const NormalizationConfig& cfg,
                                       std::string& warn)
{
    UnitsScales out = baseUnits.scales();

    if (baseUnits.mode() != UnitsMode::kSI)
    {
        if (cfg.strategy != NormalizationStrategy::kManual)
            warn = "Units mode is CodeUnits; normalization strategy ignored.";
        return out;
    }

    if (cfg.strategy == NormalizationStrategy::kManual)
        return out;

    constexpr double kEpsilon0_Fpm = PhysicalConstants::epsilon0_Fpm;
    constexpr double kElementaryCharge_C = PhysicalConstants::qe_C;

    std::string refWarn;
    const RefSpeciesInfo ref = resolveRefSpeciesInfo(cfg, refWarn);
    if (!refWarn.empty())
        warn = refWarn;

    if (cfg.strategy == NormalizationStrategy::kPlasmaFrequency)
    {
        const double n0 = out.n0_m3;
        const double q = std::abs(ref.charge_C);
        const double m = ref.mass_kg;
        if (n0 <= 0.0 || q <= 0.0 || m <= 0.0)
        {
            warn = "plasmaFrequency strategy: invalid n0/q/m; fallback to manual.";
            return out;
        }
        const double omega_p = std::sqrt(n0 * q * q / (kEpsilon0_Fpm * m));
        if (!std::isfinite(omega_p) || omega_p <= 0.0)
        {
            warn = "plasmaFrequency strategy: omega_p invalid; fallback to manual.";
            return out;
        }
        out.T0_s = 1.0 / omega_p;
        out.v0_mps = -1.0;
        out.E0_Vpm = -1.0;
        out.rho0_Cpm3 = -1.0;
        return out;
    }

    if (cfg.strategy == NormalizationStrategy::kThermal)
    {
        if (cfg.Te0_eV <= 0.0)
        {
            warn = "thermal strategy: Te0 missing or non-positive; fallback to manual.";
            return out;
        }

        const double qabs = std::max(std::abs(ref.charge_C), kElementaryCharge_C);
        const double kTe_J = cfg.Te0_eV * kElementaryCharge_C;
        out.Phi0_V = kTe_J / qabs;

        if (ref.mass_kg <= 0.0 || out.L0_m <= 0.0)
        {
            warn = "thermal strategy: invalid refMass or L0; fallback to manual.";
            return out;
        }

        const double vth = std::sqrt(kTe_J / ref.mass_kg);
        if (!std::isfinite(vth) || vth <= 0.0)
        {
            warn = "thermal strategy: vth invalid; fallback to manual.";
            return out;
        }
        out.T0_s = out.L0_m / vth;
        out.v0_mps = -1.0;
        out.E0_Vpm = -1.0;
        out.rho0_Cpm3 = -1.0;
        return out;
    }

    return out;
}

