#include "Core/Units/UnitParser.h"

#include <cctype>
#include <cstdlib>
#include <cmath>
#include <stdexcept>
#include <unordered_map>

namespace
{
// 1 Torr = 133.32236842105263 Pa (exact definition uses 1 atm / 760).
constexpr double kTorrToPa = 101325.0 / 760.0;

// Fundamental constants (SI). We keep them local to avoid cross-module coupling.
constexpr double kElementaryCharge_C = 1.602176634e-19; // exact by SI definition
constexpr double kEVToJ = kElementaryCharge_C; // 1 eV = e joules
constexpr double kAMUToKg = 1.66053906660e-27; // 1 u (amu) in kg
constexpr double kElectronMassKg = 9.1093837015e-31; // electron mass in kg


bool IsNumberChar(char c)
{
    return (std::isdigit(static_cast<unsigned char>(c)) != 0) || c == '+' || c == '-' || c == '.' || c == 'e' || c == 'E';
}
} // namespace

std::string UnitParser::Trim(const std::string& s)
{
    size_t b = 0;
    while(b < s.size() && std::isspace(static_cast<unsigned char>(s[b])) != 0)
        ++b;
    size_t e = s.size();
    while(e > b && std::isspace(static_cast<unsigned char>(s[e - 1])) != 0)
        --e;
    return s.substr(b, e - b);
}

std::string UnitParser::ToLower(std::string s)
{
    for(char& c : s)
    {
        c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    }
    return s;
}

std::string UnitParser::NormalizeMicro(std::string s)
{
    // Replace UTF-8 micro signs with 'u'. Both "µ" (U+00B5) and "μ" (U+03BC) appear in the wild.
    const std::string micro1 = u8"µ";
    const std::string micro2 = u8"μ";

    for(const auto& micro : {micro1, micro2})
    {
        size_t pos = 0;
        while((pos = s.find(micro, pos)) != std::string::npos)
        {
            s.replace(pos, micro.size(), "u");
            pos += 1;
        }
    }
    return s;
}

bool UnitParser::SplitNumberAndUnit(const std::string& s, std::string& number, std::string& unit)
{
    // Strategy:
    // 1) Trim.
    // 2) Scan for first non-number char; split there.
    // Supports both "1e-3 s" and "1e-3s".
    const std::string t = Trim(s);
    if(t.empty())
        return false;

    // Find the start of the unit token.
    size_t i = 0;
    while(i < t.size() && IsNumberChar(t[i]))
        ++i;

    number = Trim(t.substr(0, i));
    unit = Trim(t.substr(i));

    // If the unit is separated by whitespace, Trim already handled it.
    if(number.empty())
        return false;

    return true;
}

bool UnitParser::UnitToSIMultiplier(const std::string& unit_norm, double& multiplier)
{
    // All multipliers map the unit to the SI base unit for the corresponding dimension:
    // length -> meter, time -> second, voltage -> volt, pressure -> pascal.
    static const std::unordered_map<std::string, double> kUnitToSI = {
        // Length
        {"m", 1.0},
        {"cm", 1e-2},
        {"mm", 1e-3},
        {"um", 1e-6},
        {"nm", 1e-9},
        {"km", 1e3},

        // Time
        {"s", 1.0},
        {"ms", 1e-3},
        {"us", 1e-6},
        {"ns", 1e-9},
        {"min", 60.0},
        {"h", 3600.0},
        {"hr", 3600.0},

        // Frequency
        {"hz", 1.0},
        {"1/s", 1.0},
        {"s^-1", 1.0},
        {"khz", 1e3},
        {"mhz", 1e6},
        {"ghz", 1e9},

        // Velocity
        {"m/s", 1.0},
        {"mps", 1.0},
        {"cm/s", 1e-2},
        {"cmps", 1e-2},
        {"mm/s", 1e-3},
        {"mmps", 1e-3},
        {"km/s", 1e3},
        {"kmps", 1e3},

        // Voltage / potential
        {"v", 1.0},
        {"mv", 1e-3},
        {"kv", 1e3},

        // Electric field
        {"v/m", 1.0},
        {"mv/m", 1e-3},
        {"kv/m", 1e3},

        // Temperature
        {"k", 1.0},

        // Mass
        {"kg", 1.0},
        {"g", 1e-3},
        {"amu", kAMUToKg},
        {"u", kAMUToKg},
        {"me", kElectronMassKg},
        {"m_e", kElectronMassKg},

        // Charge
        {"c", 1.0},
        {"e", kElementaryCharge_C},

        // Energy
        {"j", 1.0},
        {"ev", kEVToJ},

        // Area (cross section)
        {"m2", 1.0},
        {"m^2", 1.0},
        {"cm2", 1e-4},
        {"cm^2", 1e-4},

        // Pressure
        {"pa", 1.0},
        {"kpa", 1e3},
        {"mpa", 1e6},
        {"bar", 1e5},
        {"atm", 101325.0},
        {"torr", kTorrToPa},
        {"mtorr", 1e-3 * kTorrToPa},

        // Number density (convenience)
        {"1/m3", 1.0},
        {"1/m^3", 1.0},
        {"m-3", 1.0},
        {"m^-3", 1.0},
        {"1/cm3", 1e6},
        {"1/cm^3", 1e6},
        {"cm-3", 1e6},
        {"cm^-3", 1e6},

        // Charge density (convenience)
        {"c/m3", 1.0},
        {"c/m^3", 1.0},
        {"c/cm3", 1e6},
        {"c/cm^3", 1e6},
    };

    auto it = kUnitToSI.find(unit_norm);
    if(it == kUnitToSI.end())
        return false;
    multiplier = it->second;
    return true;
}

bool UnitParser::TryParseToSI(const std::string& text, Result& out)
{
    out = Result{};

    std::string number_token;
    std::string unit_token;
    if(!SplitNumberAndUnit(text, number_token, unit_token))
        return false;

    char* end_ptr = nullptr;
    const double value = std::strtod(number_token.c_str(), &end_ptr);
    if(end_ptr == number_token.c_str())
        return false;

    // If there is no unit suffix, treat it as SI already.
    if(unit_token.empty())
    {
        out.value_si = value;
        out.unit.clear();
        out.has_unit = false;
        return true;
    }

    std::string unit_norm = NormalizeMicro(ToLower(unit_token));
    double mult = 1.0;
    if(!UnitToSIMultiplier(unit_norm, mult))
        return false;

    out.value_si = value * mult;
    out.unit = unit_norm;
    out.has_unit = true;
    return true;
}

UnitParser::Result UnitParser::ParseToSIOrThrow(const std::string& text)
{
    Result r;
    if(!TryParseToSI(text, r))
        throw std::runtime_error("UnitParser: failed to parse '" + text + "'");
    return r;
}

bool UnitParser::TryParseQuantityToSI(const std::string& text, Dimension expected, Result& out)
{
    out = Result{};
    std::string numeric, unit;
    if (!SplitNumberAndUnit(text, numeric, unit) || unit.empty()) return false;
    std::size_t used = 0;
    double value = 0.0;
    try { value = std::stod(numeric, &used); }
    catch (...) { return false; }
    if (used != numeric.size() || !std::isfinite(value)) return false;
    const auto normalized = NormalizeMicro(ToLower(unit));
    static const std::unordered_map<std::string, Dimension> dimensions{
        {"m/s",Dimension::Velocity},{"cm/s",Dimension::Velocity},{"mm/s",Dimension::Velocity},
        {"km/s",Dimension::Velocity},{"1/s",Dimension::InverseTime},{"s^-1",Dimension::InverseTime},
        {"m",Dimension::Length},{"cm",Dimension::Length},{"mm",Dimension::Length},
        {"um",Dimension::Length},{"nm",Dimension::Length},{"km",Dimension::Length},
        {"s",Dimension::Time},{"ms",Dimension::Time},{"us",Dimension::Time},
        {"ns",Dimension::Time},{"min",Dimension::Time},{"h",Dimension::Time},{"hr",Dimension::Time},
        {"v",Dimension::Voltage},{"mv",Dimension::Voltage},{"kv",Dimension::Voltage},
        {"pa",Dimension::Pressure},{"kpa",Dimension::Pressure},{"mpa",Dimension::Pressure},
        {"bar",Dimension::Pressure},{"atm",Dimension::Pressure},{"torr",Dimension::Pressure},{"mtorr",Dimension::Pressure},
        {"kg/m^3",Dimension::MassDensity},{"kg/m3",Dimension::MassDensity},{"g/cm^3",Dimension::MassDensity},
        {"pa*s",Dimension::DynamicViscosity},{"n*s/m^2",Dimension::DynamicViscosity},
        {"s/m",Dimension::Conductivity},{"n/m",Dimension::SurfaceTension},
        {"c/m^2",Dimension::SurfaceChargeDensity},{"c/m2",Dimension::SurfaceChargeDensity},
        {"c/m^3",Dimension::VolumeChargeDensity},{"c/cm^3",Dimension::VolumeChargeDensity},
        {"1/m^3",Dimension::NumberDensity},{"1/cm^3",Dimension::NumberDensity}
    };
    const auto found = dimensions.find(normalized);
    double factor = 0.0;
    if (found == dimensions.end() || found->second != expected) return false;
    if (!UnitToSIMultiplier(normalized, factor))
    {
        // New material spellings belong to the dimension-checked API only;
        // do not broaden accepted units in legacy geometry/PIC readers.
        static const std::unordered_map<std::string, double> materialFactors{
            {"kg/m^3",1.0},{"kg/m3",1.0},{"g/cm^3",1e3},
            {"pa*s",1.0},{"n*s/m^2",1.0},{"s/m",1.0},{"n/m",1.0},
            {"c/m^2",1.0},{"c/m2",1.0}
        };
        const auto material = materialFactors.find(normalized);
        if (material == materialFactors.end()) return false;
        factor = material->second;
    }
    const double si = value * factor;
    if (!std::isfinite(si) || (value != 0.0 && si == 0.0)) return false;
    out.value_si = si; out.unit = normalized; out.has_unit = true;
    return true;
}
