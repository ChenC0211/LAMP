#include "Modules/Plasma/ExplicitPIC/IO/DictionaryPIC.h"

#include "Core/Units/UnitParser.h"
#include "Core/Units/UnitsContext.h"
#include "Modules/Plasma/Physics/PhysicalConstants.h"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stack>

using namespace std;

string DictionaryPIC::toLowerCopy(string s)
{
    transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return static_cast<char>(tolower(c)); });
    return s;
}

bool DictionaryPIC::parseBoolToken(const string& s, bool& out)
{
    const string v = toLowerCopy(s);
    if (v == "true" || v == "1" || v == "on" || v == "yes")
    {
        out = true;
        return true;
    }
    if (v == "false" || v == "0" || v == "off" || v == "no")
    {
        out = false;
        return true;
    }
    return false;
}

PicBCType DictionaryPIC::parseBCTypeToken(const string& s, bool& ok)
{
    const string v = toLowerCopy(s);
    ok = true;
    if (v == "dirichlet")
        return PicBCType::Dirichlet;
    if (v == "neumann")
        return PicBCType::Neumann;
    if (v == "robin")
        return PicBCType::Robin;
    if (v == "periodic")
        return PicBCType::Periodic;

    ok = false;
    return PicBCType::Dirichlet;
}

static bool isQuotedName(const string& s)
{
    return s.size() >= 2 && s.front() == '"' && s.back() == '"';
}

static string unquote(const string& s)
{
    if (isQuotedName(s))
        return s.substr(1, s.size() - 2);
    return s;
}

namespace
{
std::string trim_ws(const std::string& s)
{
    size_t b = 0;
    while (b < s.size() && std::isspace(static_cast<unsigned char>(s[b])) != 0)
        ++b;
    size_t e = s.size();
    while (e > b && std::isspace(static_cast<unsigned char>(s[e - 1])) != 0)
        --e;
    return s.substr(b, e - b);
}

std::string lowerCopy(const std::string& s)
{
    std::string out = s;
    std::transform(out.begin(), out.end(), out.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return out;
}

std::string stripListPunct(std::string s)
{
    auto isPunct = [](char c) {
        return c == '(' || c == ')' || c == '[' || c == ']' || c == '{' || c == '}' || c == ',' || c == ';';
    };
    while (!s.empty() && isPunct(s.front()))
        s.erase(s.begin());
    while (!s.empty() && isPunct(s.back()))
        s.pop_back();
    return s;
}

void parseIntListTokens(const std::vector<std::string>& tokens,
                         size_t start,
                         std::vector<int>& out,
                         bool oneBasedToZero)
{
    out.clear();
    if (start >= tokens.size())
        return;

    std::string text;
    for (size_t i = start; i < tokens.size(); ++i)
    {
        if (!text.empty())
            text.push_back(' ');
        text += unquote(tokens[i]);
    }

    text = stripListPunct(text);
    for (char& c : text)
    {
        if (c == ',' || c == ';' || c == '|' || c == '/')
            c = ' ';
    }

    std::istringstream iss(text);
    std::string part;
    while (iss >> part)
    {
        part = stripListPunct(part);
        if (part.empty())
            continue;
        try
        {
            int v = std::stoi(part);
            if (oneBasedToZero && v > 0)
                v -= 1;
            out.push_back(v);
        }
        catch (...)
        {
            // ignore invalid entries
        }
    }
}

void parseBoundaryGroupTokens(const std::vector<std::string>& tokens,
                              size_t start,
                              std::vector<int>& out)
{
    parseIntListTokens(tokens, start, out, false);
}

// Forward declaration (implementation below).
bool tryParseStrictDouble(const std::string& s, double& out);

// Parse a 2D vector written like:
//   (1e3, 0) m/s
//   (1e3 m/s, 0 m/s)
// If no unit is present, numeric values are interpreted as SI.
bool parseVector2ToSIWithUnits(const std::string& line, Vector2D& out)
{
    const size_t l = line.find('(');
    const size_t r = line.find(')', l == std::string::npos ? 0 : l + 1);
    if (l == std::string::npos || r == std::string::npos || r <= l)
        return false;

    const std::string inside = trim_ws(line.substr(l + 1, r - l - 1));
    const size_t comma = inside.find(',');
    if (comma == std::string::npos)
        return false;

    const std::string xs = trim_ws(inside.substr(0, comma));
    const std::string ys = trim_ws(inside.substr(comma + 1));

    // Optional vector-level unit after ')'. Take the first token, stripping punctuation.
    std::string unit;
    {
        std::string tail = trim_ws(line.substr(r + 1));
        tail = stripListPunct(tail);
        std::istringstream iss(tail);
        iss >> unit;
        unit = stripListPunct(unit);
    }

    auto parseComp = [&](const std::string& s, double& out_si) -> bool {
        // 1) Try as-is (component may already include a unit).
        UnitParser::Result rr;
        if (UnitParser::TryParseToSI(s, rr))
        {
            out_si = rr.value_si;
            return true;
        }

        // 2) If there is a trailing vector unit, try appending it.
        if (!unit.empty())
        {
            if (UnitParser::TryParseToSI(s + " " + unit, rr))
            {
                out_si = rr.value_si;
                return true;
            }
        }

        // 3) Bare number -> SI.
        double numeric = 0.0;
        if (tryParseStrictDouble(s, numeric))
        {
            out_si = numeric;
            return true;
        }

        return false;
    };

    double x = 0.0;
    double y = 0.0;
    if (!parseComp(xs, x) || !parseComp(ys, y))
        return false;

    out = Vector2D(x, y);
    return true;
}

void appendUnique(std::vector<std::string>& out, const std::string& v)
{
    if (v.empty())
        return;
    if (std::find(out.begin(), out.end(), v) == out.end())
        out.push_back(v);
}

void appendSeeSideToken(const std::string& token, std::vector<std::string>& out)
{
    std::string t = lowerCopy(unquote(token));
    for (char& c : t)
    {
        if (c == ',' || c == ';' || c == '|' || c == '/')
            c = ' ';
    }

    std::istringstream iss(t);
    std::string part;
    while (iss >> part)
    {
        part = stripListPunct(part);
        if (part.empty())
            continue;
        if (part == "both" || part == "lr")
        {
            appendUnique(out, "left");
            appendUnique(out, "right");
            continue;
        }
        if (part == "all")
        {
            appendUnique(out, "left");
            appendUnique(out, "right");
            appendUnique(out, "bottom");
            appendUnique(out, "top");
            continue;
        }
        if (part == "tb")
        {
            appendUnique(out, "top");
            appendUnique(out, "bottom");
            continue;
        }
        appendUnique(out, part);
    }
}

bool tryParseStrictDouble(const std::string& s, double& out)
{
    size_t idx = 0;
    try
    {
        out = std::stod(s, &idx);
    }
    catch (...)
    {
        return false;
    }
    return idx == s.size();
}

bool parseTimeToCode(const std::vector<std::string>& tokens, size_t start,
                     const UnitsContext* units, double& out_code)
{
    if (start >= tokens.size())
        return false;

    // Join the remainder tokens into a single "value unit" string.
    std::string text;
    for (size_t i = start; i < tokens.size(); ++i)
    {
        if (!text.empty())
            text += " ";
        text += tokens[i];
    }

    // Legacy pure-number format -> code value.
    if (tokens.size() == start + 1)
    {
        double numeric = 0.0;
        if (tryParseStrictDouble(tokens[start], numeric))
        {
            out_code = numeric;
            return true;
        }
    }

    UnitParser::Result r;
    if (!UnitParser::TryParseToSI(text, r))
        return false;

    // Require time units when a unit suffix is present.
    if (r.has_unit)
    {
        const std::string& u = r.unit;
        if (!(u == "s" || u == "ms" || u == "us" || u == "ns" || u == "min" || u == "h" || u == "hr"))
            return false;
    }

    // No unit suffix: treat as code value (legacy behavior).
    if (!r.has_unit)
    {
        out_code = r.value_si;
        return true;
    }

    out_code = units ? units->TimeSIToCode(r.value_si) : r.value_si;
    return true;
}

bool parsePotentialToCode(const std::vector<std::string>& tokens, size_t start,
                          const UnitsContext* units, double& out_code)
{
    if (start >= tokens.size())
        return false;

    std::string text;
    for (size_t i = start; i < tokens.size(); ++i)
    {
        if (!text.empty())
            text += " ";
        text += tokens[i];
    }

    // Legacy pure-number format -> code value.
    if (tokens.size() == start + 1)
    {
        double numeric = 0.0;
        if (tryParseStrictDouble(tokens[start], numeric))
        {
            out_code = numeric;
            return true;
        }
    }

    UnitParser::Result r;
    if (!UnitParser::TryParseToSI(text, r))
        return false;

    if (r.has_unit)
    {
        const std::string& u = r.unit;
        if (!(u == "v" || u == "mv" || u == "kv"))
            return false;
    }

    if (!r.has_unit)
    {
        out_code = r.value_si;
        return true;
    }

    out_code = units ? units->PotentialSIToCode(r.value_si) : r.value_si;
    return true;
}

bool parseScalarToSI(const std::vector<std::string>& tokens, size_t start, double& out_si)
{
    if (start >= tokens.size())
        return false;

    std::string text;
    for (size_t i = start; i < tokens.size(); ++i)
    {
        if (!text.empty())
            text += " ";
        text += tokens[i];
    }

    // Pure-number: interpret as SI.
    if (tokens.size() == start + 1)
    {
        double numeric = 0.0;
        if (tryParseStrictDouble(tokens[start], numeric))
        {
            out_si = numeric;
            return true;
        }
    }

    UnitParser::Result r;
    if (!UnitParser::TryParseToSI(text, r))
        return false;
    out_si = r.value_si;
    return true;
}

bool parseNumberDensityToSI(const std::vector<std::string>& tokens, size_t start, double& out_si)
{
    if (start >= tokens.size())
        return false;

    // Pure-number: interpret as SI (1/m^3).
    if (tokens.size() == start + 1)
    {
        double numeric = 0.0;
        if (tryParseStrictDouble(tokens[start], numeric))
        {
            out_si = numeric;
            return true;
        }
    }

    if (tokens.size() >= start + 2)
    {
        double numeric = 0.0;
        if (tryParseStrictDouble(tokens[start], numeric))
        {
            std::string unit = stripListPunct(lowerCopy(tokens[start + 1]));
            if (unit == "1/m3" || unit == "1/m^3" || unit == "m-3" || unit == "m^-3")
            {
                out_si = numeric;
                return true;
            }
            if (unit == "1/cm3" || unit == "1/cm^3" || unit == "cm-3" || unit == "cm^-3")
            {
                out_si = numeric * 1e6;
                return true;
            }
        }
    }

    // Fallback: treat as plain scalar in SI.
    return parseScalarToSI(tokens, start, out_si);
}

bool parseEnergyToEV(const std::vector<std::string>& tokens, size_t start, double& out_eV)
{
    if (start >= tokens.size())
        return false;

    // Pure-number: interpret as eV.
    if (tokens.size() == start + 1)
    {
        double numeric = 0.0;
        if (tryParseStrictDouble(tokens[start], numeric))
        {
            out_eV = numeric;
            return true;
        }
    }

    double val_si = 0.0;
    if (!parseScalarToSI(tokens, start, val_si))
        return false;

    // UnitParser returns SI Joules for eV inputs.
    out_eV = val_si / std::max(1e-300, PhysicalConstants::eV_J);
    return true;
}

// Parse Te where the preferred input is in eV.
// - If no unit is present, treat it as eV.
// - If unit is eV, UnitParser returns Joules; we convert back to eV.
bool parseTeToEV(const std::vector<std::string>& tokens, size_t start, double& out_eV)
{
    if (start >= tokens.size())
        return false;

    // Pure numeric -> assume eV.
    if (tokens.size() == start + 1)
    {
        double numeric = 0.0;
        if (tryParseStrictDouble(tokens[start], numeric))
        {
            out_eV = numeric;
            return true;
        }
    }

    std::string text;
    for (size_t i = start; i < tokens.size(); ++i)
    {
        if (!text.empty())
            text += " ";
        text += tokens[i];
    }

    UnitParser::Result r;
    if (!UnitParser::TryParseToSI(text, r))
        return false;

    if (!r.has_unit)
    {
        out_eV = r.value_si; // interpret as eV (legacy)
        return true;
    }

    // UnitParser returns SI for the specific unit's base dimension.
    // - For energy-like inputs (eV, J), value_si is in Joules.
    // - For temperature-like inputs (K), value_si is in Kelvin.
    // Convert to eV accordingly.
    if (r.unit == "k")
    {
        out_eV = (PhysicalConstants::kB_JpK * r.value_si) / PhysicalConstants::eV_J;
        return true;
    }

    // Default: interpret as energy in Joules.
    out_eV = r.value_si / PhysicalConstants::eV_J;
    return true;
}

bool parseEnergyToEVAt(const std::vector<std::string>& tokens, size_t start, double& out_eV, size_t& consumed)
{
    if (start >= tokens.size())
        return false;

    if (start + 1 < tokens.size())
    {
        std::vector<std::string> sub{tokens[start], tokens[start + 1]};
        if (parseEnergyToEV(sub, 0, out_eV))
        {
            consumed = 2;
            return true;
        }
    }

    std::vector<std::string> sub{tokens[start]};
    if (parseEnergyToEV(sub, 0, out_eV))
    {
        consumed = 1;
        return true;
    }
    return false;
}

bool parseScalarToSIAt(const std::vector<std::string>& tokens, size_t start, double& out_si, size_t& consumed)
{
    if (start >= tokens.size())
        return false;

    if (start + 1 < tokens.size())
    {
        std::vector<std::string> sub{tokens[start], tokens[start + 1]};
        if (parseScalarToSI(sub, 0, out_si))
        {
            consumed = 2;
            return true;
        }
    }

    std::vector<std::string> sub{tokens[start]};
    if (parseScalarToSI(sub, 0, out_si))
    {
        consumed = 1;
        return true;
    }
    return false;
}

bool parsePotentialValueAt(const std::vector<std::string>& tokens, size_t start,
                           const UnitsContext* units, double& out_code, size_t& consumed)
{
    if (start >= tokens.size())
        return false;

    if (start + 1 < tokens.size())
    {
        std::vector<std::string> sub{tokens[start], tokens[start + 1]};
        if (parsePotentialToCode(sub, 0, units, out_code))
        {
            consumed = 2;
            return true;
        }
    }

    std::vector<std::string> sub{tokens[start]};
    if (parsePotentialToCode(sub, 0, units, out_code))
    {
        consumed = 1;
        return true;
    }
    return false;
}

bool parseUnitValueAt(const std::vector<std::string>& tokens, size_t start,
                      const std::vector<std::string>& allowed_units,
                      double& out_si, size_t& consumed, bool allow_unitless = true)
{
    if (start >= tokens.size())
        return false;

    auto unitAllowed = [&](const std::string& u) -> bool {
        if (allowed_units.empty())
            return true;
        const std::string ul = lowerCopy(u);
        for (const auto& au : allowed_units)
        {
            if (ul == au)
                return true;
        }
        return false;
    };

    auto tryParse = [&](size_t count) -> bool {
        if (start + count > tokens.size())
            return false;
        std::string text = tokens[start];
        for (size_t i = 1; i < count; ++i)
        {
            text += " ";
            text += tokens[start + i];
        }
        UnitParser::Result r;
        if (!UnitParser::TryParseToSI(text, r))
            return false;
        if (r.has_unit)
        {
            if (!unitAllowed(r.unit))
                return false;
        }
        else if (!allow_unitless)
        {
            return false;
        }
        out_si = r.value_si;
        consumed = count;
        return true;
    };

    if (start + 1 < tokens.size() && tryParse(2))
        return true;
    if (tryParse(1))
        return true;
    return false;
}

bool parseFrequencyToHzAt(const std::vector<std::string>& tokens, size_t start,
                          double& out_hz, size_t& consumed)
{
    static const std::vector<std::string> kUnits = {"hz", "khz", "mhz", "ghz"};
    return parseUnitValueAt(tokens, start, kUnits, out_hz, consumed, /*allow_unitless=*/true);
}

bool parseTimeToSIAt(const std::vector<std::string>& tokens, size_t start,
                     double& out_s, size_t& consumed)
{
    static const std::vector<std::string> kUnits = {"s", "ms", "us", "ns", "min", "h", "hr"};
    return parseUnitValueAt(tokens, start, kUnits, out_s, consumed, /*allow_unitless=*/true);
}

bool parseAngleToRadAt(const std::vector<std::string>& tokens, size_t start,
                       double& out_rad, size_t& consumed)
{
    if (start >= tokens.size())
        return false;

    constexpr double kPi = 3.141592653589793238462643383279502884;

    if (start + 1 < tokens.size())
    {
        const std::string u = lowerCopy(tokens[start + 1]);
        if (u == "deg" || u == "rad")
        {
            double v = 0.0;
            if (!tryParseStrictDouble(tokens[start], v))
                return false;
            consumed = 2;
            out_rad = (u == "deg") ? (v * kPi / 180.0) : v;
            return true;
        }
    }

    double v = 0.0;
    if (!tryParseStrictDouble(tokens[start], v))
        return false;
    consumed = 1;
    out_rad = v;
    return true;
}

bool parseTimeDirichletTokens(const std::vector<std::string>& tokens, size_t start,
                              const UnitsContext* units,
                              PicBoundarySide::TimeFunction& out,
                              double fallback_dc_code)
{
    if (start >= tokens.size())
        return false;

    out = PicBoundarySide::TimeFunction{};
    out.dc = fallback_dc_code;
    out.enabled = true;
    out.type = "sine";

    std::string head = lowerCopy(tokens[start]);
    if (head == "off" || head == "none" || head == "disable" || head == "false")
    {
        out.enabled = false;
        out.type = "constant";
        return true;
    }
    if (head == "sine" || head == "sin")
    {
        out.type = "sine";
        ++start;
    }
    else if (head == "cosine" || head == "cos")
    {
        out.type = "cosine";
        ++start;
    }
    else if (head == "constant" || head == "const")
    {
        out.enabled = false;
        out.type = "constant";
        return true;
    }

    size_t i = start;
    while (i < tokens.size())
    {
        const std::string key = lowerCopy(tokens[i]);
        if (key == "vdc" || key == "dc" || key == "offset")
        {
            double v = 0.0;
            size_t n = 0;
            if (parsePotentialValueAt(tokens, i + 1, units, v, n))
            {
                out.dc = v;
                i += 1 + n;
                continue;
            }
        }
        else if (key == "vac" || key == "ac" || key == "amp" || key == "amplitude")
        {
            double v = 0.0;
            size_t n = 0;
            if (parsePotentialValueAt(tokens, i + 1, units, v, n))
            {
                out.ac = v;
                i += 1 + n;
                continue;
            }
        }
        else if (key == "vpp" || key == "pp")
        {
            double vpp = 0.0;
            size_t n = 0;
            if (parsePotentialValueAt(tokens, i + 1, units, vpp, n))
            {
                out.ac = 0.5 * vpp;
                i += 1 + n;
                continue;
            }
        }
        else if (key == "freq" || key == "f" || key == "frequency")
        {
            double f = 0.0;
            size_t n = 0;
            if (parseFrequencyToHzAt(tokens, i + 1, f, n))
            {
                out.freq_hz = f;
                i += 1 + n;
                continue;
            }
        }
        else if (key == "phase" || key == "phi")
        {
            double rad = 0.0;
            size_t n = 0;
            if (parseAngleToRadAt(tokens, i + 1, rad, n))
            {
                out.phase_rad = rad;
                i += 1 + n;
                continue;
            }
        }
        else if (key == "t0" || key == "delay" || key == "start")
        {
            double t = 0.0;
            size_t n = 0;
            if (parseTimeToSIAt(tokens, i + 1, t, n))
            {
                out.t0_s = t;
                i += 1 + n;
                continue;
            }
        }

        // Unknown or malformed token: skip.
        ++i;
    }

    return true;
}
} // namespace

bool DictionaryPIC::readDict(const string& picDictPath, const UnitsContext* units, const AdditionalRecord& additionalRecord)
{
    ifstream file(picDictPath);
    if (!file.is_open())
    {
        // pic.dict is optional at runtime; caller decides whether this is fatal.
        return false;
    }

    stack<string> blocks;
    string pendingBlock;
    int activeSpeciesIndex = -1;
    int pendingFieldBoundaryId = -1;
    int pendingParticleBoundaryId = -1;
    int activeFieldBoundaryId = -1;
    int activeParticleBoundaryId = -1;
    bool sawPicBlock = false;
    bool sawSpeciesBlock = false;
    bool sawSpeciesEntry = false;
    bool sawDt = false;
    bool sawSteps = false;
    bool sawOutputEvery = false;
    bool sawPeriodicDiagnosticsEnabled = false;
    bool sawCurrentWaveformEnabled = false;
    bool sawProfile1DEnabled = false;
    bool sawProfile1DAvgEnabled = false;
    bool sawCurrentBalanceEnabled = false;
    bool sawDiagnosticsEnabled = false;
    bool sawElectrodeCurrentEnabled = false;
    bool sawReactionDiagnosticsEnabled = false;
    bool sawEnergyHistogramEnabled = false;
    bool warnedDeprecatedEmbeddedGeometryBlock = false;

    auto getOrCreateFieldBC = [&](int id) -> PicConfig::BoundaryIdFieldBC& {
        for (auto& b : config_.fieldBCById)
        {
            if (b.id == id)
                return b;
        }
        PicConfig::BoundaryIdFieldBC b;
        b.id = id;
        config_.fieldBCById.push_back(b);
        return config_.fieldBCById.back();
    };

    auto getOrCreateParticleBC = [&](int id) -> PicConfig::BoundaryIdParticleBC& {
        for (auto& b : config_.particleBCById)
        {
            if (b.id == id)
                return b;
        }
        PicConfig::BoundaryIdParticleBC b;
        b.id = id;
        config_.particleBCById.push_back(b);
        return config_.particleBCById.back();
    };

    auto openBlock = [&](const string& name) {
        const string lower = toLowerCopy(name);
        if (lower == "pic")
            sawPicBlock = true;
        else if (lower == "species")
            sawSpeciesBlock = true;
        else if (lower == "fieldboundaryitem")
        {
            activeFieldBoundaryId = pendingFieldBoundaryId;
            if (activeFieldBoundaryId >= 0)
                (void)getOrCreateFieldBC(activeFieldBoundaryId);
        }
        else if (lower == "particleboundaryitem")
        {
            activeParticleBoundaryId = pendingParticleBoundaryId;
            if (activeParticleBoundaryId >= 0)
                (void)getOrCreateParticleBC(activeParticleBoundaryId);
        }
        blocks.push(name);
    };

    auto closeBlock = [&]() {
        if (blocks.empty())
            return;
        if (blocks.top() == "SpeciesItem")
            activeSpeciesIndex = -1;
        if (toLowerCopy(blocks.top()) == "fieldboundaryitem")
            activeFieldBoundaryId = -1;
        if (toLowerCopy(blocks.top()) == "particleboundaryitem")
            activeParticleBoundaryId = -1;
        blocks.pop();
    };

    auto currentBlock = [&]() -> string {
        return blocks.empty() ? string() : blocks.top();
    };
    auto parentBlock = [&]() -> string {
        if (blocks.size() < 2)
            return string();
        auto temp = blocks;
        temp.pop();
        return temp.top();
    };

    auto setDomainSide = [&](const string& side, PicBCType type, double value) {
        const string s = toLowerCopy(side);
        PicBoundarySide* bc = nullptr;
        if (s == "left")
            bc = &config_.domainBC.left;
        else if (s == "right")
            bc = &config_.domainBC.right;
        else if (s == "bottom")
            bc = &config_.domainBC.bottom;
        else if (s == "top")
            bc = &config_.domainBC.top;

        if (bc)
        {
            bc->type = type;
            bc->value = value;
        }
    };

    auto getDomainSide = [&](const string& side) -> PicBoundarySide* {
        const string s = toLowerCopy(side);
        if (s == "left")
            return &config_.domainBC.left;
        if (s == "right")
            return &config_.domainBC.right;
        if (s == "bottom")
            return &config_.domainBC.bottom;
        if (s == "top")
            return &config_.domainBC.top;
        return nullptr;
    };

    string line;
    while (getline(file, line))
    {
        line = removeInlineComments(line);
        line = trim(line);

        if (line.empty())
            continue;

        // If we previously saw a "block title" line (stored in pendingBlock), the next
        // non-empty line must be the opening brace "{". If it isn't, the title was not
        // actually a block header (e.g. OpenFOAM-style list keys like "fields" followed by
        // "("), so we must drop it to avoid corrupting the brace stack later.
        if (!pendingBlock.empty() && line != "{" && line != "}")
        {
            pendingBlock.clear();
        }

        // -----------------------------
        // Braces
        // -----------------------------
        if (line == "{")
        {
            if (!pendingBlock.empty())
            {
                openBlock(pendingBlock);
                pendingBlock.clear();
            }
            continue;
        }

        if (line == "}")
        {
            closeBlock();
            continue;
        }

        // Inline "NAME {" opening
        const size_t bracePos = line.find('{');
        if (bracePos != string::npos)
        {
            // Special handling for one-line inline blocks that both open and close on the
            // same line, e.g.:
            //   boundary 3 { type Periodic }
            //   boundary 3 { model absorb SEE { enabled false } }
            // The generic brace-stack logic below intentionally stops at '{' and would
            // discard the inline content and miss the closing '}', corrupting parsing.
            const size_t closePos = line.rfind('}');
            if (closePos != string::npos && closePos > bracePos)
            {
                const string before = trim(line.substr(0, bracePos));
                const string inside = trim(line.substr(bracePos + 1, closePos - (bracePos + 1)));

                const string parent = toLowerCopy(currentBlock());
                if (toLowerCopy(before) == "see")
                {
                    SEEConfig* seeCfg = &config_.see;
                    if (toLowerCopy(currentBlock()) == "particleboundaryitem" && activeParticleBoundaryId >= 0)
                    {
                        auto& entry = getOrCreateParticleBC(activeParticleBoundaryId);
                        entry.seeSpecified = true;
                        seeCfg = &entry.see;
                    }

                    const vector<string> innerToks = split(inside);
                    for (size_t iTok = 0; iTok < innerToks.size();)
                    {
                        const string k = toLowerCopy(innerToks[iTok]);
                        if (iTok + 1 >= innerToks.size())
                            break;
                        if (k == "enabled")
                        {
                            parseBoolToken(innerToks[iTok + 1], seeCfg->enabled);
                            iTok += 2;
                        }
                        else if (k == "cathodeside" || k == "cathodesides")
                        {
                            seeCfg->cathodeSides.clear();
                            for (size_t i = iTok + 1; i < innerToks.size(); ++i)
                                appendSeeSideToken(innerToks[i], seeCfg->cathodeSides);
                            if (!seeCfg->cathodeSides.empty())
                                seeCfg->cathodeSide = seeCfg->cathodeSides.front();
                            else if (iTok + 1 < innerToks.size())
                                seeCfg->cathodeSide = unquote(innerToks[iTok + 1]);
                            break;
                        }
                        else if (k == "emittedspecies")
                        {
                            seeCfg->emittedSpecies = unquote(innerToks[iTok + 1]);
                            iTok += 2;
                        }
                        else if (k == "gamma")
                        {
                            seeCfg->gamma = stod(innerToks[iTok + 1]);
                            iTok += 2;
                        }
                        else if (k == "gammamodel" || k == "gamma_model" || k == "model")
                        {
                            seeCfg->gammaModel = unquote(innerToks[iTok + 1]);
                            iTok += 2;
                        }
                        else if (k == "gamma0")
                        {
                            seeCfg->gamma0 = stod(innerToks[iTok + 1]);
                            iTok += 2;
                        }
                        else if (k == "eref" || k == "eref_ev" || k == "eref_eV" || k == "erefenergy" || k == "erefenergy_ev")
                        {
                            double e_eV = 0.0;
                            size_t consumed = 0;
                            if (parseEnergyToEVAt(innerToks, iTok + 1, e_eV, consumed))
                                seeCfg->Eref_eV = e_eV;
                            iTok += 1 + std::max<size_t>(1, consumed);
                        }
                        else if (k == "alpha")
                        {
                            seeCfg->alpha = stod(innerToks[iTok + 1]);
                            iTok += 2;
                        }
                        else if (k == "gammamax" || k == "gamma_max")
                        {
                            seeCfg->gammaMax = stod(innerToks[iTok + 1]);
                            iTok += 2;
                        }
                        else if (k == "emin" || k == "emin_ev" || k == "emin_eV")
                        {
                            double e_eV = 0.0;
                            size_t consumed = 0;
                            if (parseEnergyToEVAt(innerToks, iTok + 1, e_eV, consumed))
                                seeCfg->Emin_eV = e_eV;
                            iTok += 1 + std::max<size_t>(1, consumed);
                        }
                        else if (k == "piecewise" || k == "gammapiecewise" || k == "gammatable" || k == "table")
                        {
                            seeCfg->piecewiseE_eV.clear();
                            seeCfg->piecewiseGamma.clear();
                            const size_t n = innerToks.size() - iTok;
                            if (n >= 5 && ((n - 1) % 2 == 0))
                            {
                                for (size_t i = iTok + 1; i + 1 < innerToks.size(); i += 2)
                                {
                                    double E = 0.0, g = 0.0;
                                    if (!tryParseStrictDouble(innerToks[i], E) || !tryParseStrictDouble(innerToks[i + 1], g))
                                        continue;
                                    seeCfg->piecewiseE_eV.push_back(E);
                                    seeCfg->piecewiseGamma.push_back(g);
                                }
                            }
                            break;
                        }
                        else if (k == "emittedenergy" || k == "emittedenergy_ev" || k == "emittedenergy_eV")
                        {
                            double e_eV = 0.0;
                            size_t consumed = 0;
                            if (parseEnergyToEVAt(innerToks, iTok + 1, e_eV, consumed))
                                seeCfg->emittedEnergy_eV = e_eV;
                            iTok += 1 + std::max<size_t>(1, consumed);
                        }
                        else if (k == "emittedenergymodel" || k == "emitted_energy_model" || k == "energymodel" || k == "emittedmodel")
                        {
                            seeCfg->emittedEnergyModel = unquote(innerToks[iTok + 1]);
                            iTok += 2;
                        }
                        else if (k == "emittedte" || k == "emittedte_ev" || k == "emittedte_eV" || k == "emittedtemperature" || k == "te_emit" || k == "teemit")
                        {
                            double te_eV = 0.0;
                            size_t consumed = 0;
                            if (parseEnergyToEVAt(innerToks, iTok + 1, te_eV, consumed))
                                seeCfg->emittedTe_eV = te_eV;
                            iTok += 1 + std::max<size_t>(1, consumed);
                        }
                        else if (k == "emittedemean" || k == "emittedmean" || k == "emean" || k == "emean_ev" || k == "emean_eV" || k == "emittedemean_ev" || k == "emittedemean_eV")
                        {
                            double em_eV = 0.0;
                            size_t consumed = 0;
                            if (parseEnergyToEVAt(innerToks, iTok + 1, em_eV, consumed))
                                seeCfg->emittedEmean_eV = em_eV;
                            iTok += 1 + std::max<size_t>(1, consumed);
                        }
                        else if (k == "injecteps" || k == "injecteps_m")
                        {
                            double eps_si = 0.0;
                            size_t consumed = 0;
                            if (parseScalarToSIAt(innerToks, iTok + 1, eps_si, consumed))
                                seeCfg->injectEps_m = eps_si;
                            iTok += 1 + std::max<size_t>(1, consumed);
                        }
                        else if (k == "maxemitmacro")
                        {
                            seeCfg->maxEmitMacro = stoi(innerToks[iTok + 1]);
                            iTok += 2;
                        }
                        else if (k == "angular" || k == "angularmodel" || k == "distribution")
                        {
                            seeCfg->angular = unquote(innerToks[iTok + 1]);
                            iTok += 2;
                        }
                        else
                        {
                            ++iTok;
                        }
                    }
                    continue;
                }

                if ((parent == "field" || parent == "particle") && !before.empty())
                {
                    const vector<string> headerToks = split(before);
                    if (headerToks.size() >= 2 && toLowerCopy(headerToks[0]) == "boundary")
                    {
                        const int id = stoi(headerToks[1]);
                        const vector<string> innerToks = split(inside);

                        if (parent == "field")
                        {
                            auto& entry = getOrCreateFieldBC(id);
                            for (size_t iTok = 0; iTok + 1 < innerToks.size(); ++iTok)
                            {
                                const string k = toLowerCopy(innerToks[iTok]);
                                if (k == "type")
                                {
                                    bool ok = false;
                                    entry.bc.type = parseBCTypeToken(innerToks[iTok + 1], ok);
                                    ++iTok;
                                }
                                else if (k == "value")
                                {
                                    // Minimal inline support: expect value in the form "value <num> <unit>".
                                    if (iTok + 2 < innerToks.size())
                                    {
                                        vector<string> tmp;
                                        tmp.push_back("value");
                                        tmp.push_back(innerToks[iTok + 1]);
                                        tmp.push_back(innerToks[iTok + 2]);
                                        double value_code = 0.0;
                                        if (!parsePotentialToCode(tmp, 1, units, value_code))
                                        {
                                            std::cerr << "[DictionaryPIC] ERROR: failed to parse boundary BC value: '" << line << "'" << std::endl;
                                            exit(1);
                                        }
                                        entry.bc.value = value_code;
                                        iTok += 2;
                                    }
                                }
                            }
                        }
                        else
                        {
                            auto& entry = getOrCreateParticleBC(id);
                            for (size_t iTok = 0; iTok + 1 < innerToks.size(); ++iTok)
                            {
                                const string k = toLowerCopy(innerToks[iTok]);
                                if (k == "model")
                                {
                                    entry.model = unquote(innerToks[iTok + 1]);
                                    break;
                                }
                            }
                        }

                        continue;
                    }
                }
            }

            string before = trim(line.substr(0, bracePos));
            if (before.empty())
            {
                if (!pendingBlock.empty())
                {
                    openBlock(pendingBlock);
                    pendingBlock.clear();
                }
                continue;
            }

            // Special case: species item name may appear as: "electron" {
            if (toLowerCopy(currentBlock()) == "species" && isQuotedName(before))
            {
                SpeciesConfig sp;
                sp.name = unquote(before);
                config_.species.push_back(sp);
                activeSpeciesIndex = static_cast<int>(config_.species.size()) - 1;
                sawSpeciesEntry = true;
                openBlock("SpeciesItem");
                continue;
            }

            // Inline boundary item under BC/Field or BC/Particle
            const string parent = toLowerCopy(currentBlock());
            if (parent == "field" || parent == "particle")
            {
                vector<string> toks = split(before);
                if (toks.size() >= 2 && toLowerCopy(toks[0]) == "boundary")
                {
                    const int id = stoi(toks[1]);
                    if (parent == "field")
                    {
                        pendingFieldBoundaryId = id;
                        openBlock("FieldBoundaryItem");
                    }
                    else
                    {
                        pendingParticleBoundaryId = id;
                        openBlock("ParticleBoundaryItem");
                    }
                    continue;
                }
            }

            openBlock(before);
            continue;
        }

        // -----------------------------
        // Block titles (waiting for '{')
        // -----------------------------
        {
            const string key = toLowerCopy(line);
            if (key == "pic" || key == "species" || key == "bc" || key == "deposition" || key == "output" || key == "domain" || key == "embedded" ||
                key == "embeddedgeometry" ||
                key == "see" || key == "particlemanagement" || key == "amr" || key == "debye" || key == "indicator" || key == "particlerelabel" ||
                key == "refinementindicators" || key == "debyelength" || key == "electricfield" ||
                key == "report" || key == "console")
            {
                pendingBlock = line;
                continue;
            }
        }

        // Species name line (waiting for '{')
        if (toLowerCopy(currentBlock()) == "species" && isQuotedName(line))
        {
            SpeciesConfig sp;
            sp.name = unquote(line);
            config_.species.push_back(sp);
            activeSpeciesIndex = static_cast<int>(config_.species.size()) - 1;
            sawSpeciesEntry = true;
            pendingBlock = "SpeciesItem";
            continue;
        }

        // Generic block-title support: allow arbitrary nested blocks written as
        //
        //   someBlock
        //   {
        //     ...
        //   }
        //
        // This is needed for user-authored sub-blocks such as initBox in the
        // Species section. Even if we don't parse the inner keys, we must keep
        // brace nesting consistent so later blocks (e.g. other species) are
        // parsed correctly.
        if (pendingBlock.empty())
        {
            const vector<string> tok = split(line);
            if (tok.size() == 1 && !isQuotedName(line))
            {
                pendingBlock = line;
                continue;
            }
        }

        // Boundary item under BC/Field or BC/Particle (block-style)
        {
            const string parent = toLowerCopy(currentBlock());
            const vector<string> tok = split(line);
            if ((parent == "field" || parent == "particle") && tok.size() >= 2 && toLowerCopy(tok[0]) == "boundary")
            {
                const int id = stoi(tok[1]);
                if (parent == "field")
                {
                    pendingFieldBoundaryId = id;
                    pendingBlock = "FieldBoundaryItem";
                }
                else
                {
                    pendingParticleBoundaryId = id;
                    pendingBlock = "ParticleBoundaryItem";
                }
                continue;
            }
        }

        // -----------------------------
        // Key-value parsing within blocks
        // -----------------------------
        const string b = currentBlock();
        const vector<string> tokens = split(line);
        if (tokens.empty())
            continue;
        if (additionalRecord && additionalRecord(b, parentBlock(), line, tokens))
            continue;

        // PIC block
        if (toLowerCopy(b) == "pic")
        {
            if (tokens.size() >= 2)
            {
                const string k = toLowerCopy(tokens[0]);
                if (k == "dt")
                {
                    double dt_code = 0.0;
                    if (!parseTimeToCode(tokens, 1, units, dt_code))
                    {
                        std::cerr << "[DictionaryPIC] ERROR: failed to parse dt: '" << line << "'" << std::endl;
                        exit(1);
                    }
                    config_.dt = dt_code;
                    sawDt = true;
                }
                else if (k == "steps" && tokens.size() == 2)
                {
                    config_.steps = stoi(tokens[1]);
                    sawSteps = true;
                }
                else if ((k == "randomseed" || k == "seed") && tokens.size() == 2)
                    config_.randomSeed = static_cast<unsigned int>(stoul(tokens[1]));
                else if ((k == "omp" || k == "ompthreads" || k == "omp_threads") && tokens.size() == 2)
                    config_.ompThreads = stoi(tokens[1]);
                else if ((k == "segmentermaxlevel" || k == "predictorsegmentermaxlevel") && tokens.size() == 2)
                    config_.segmenterMaxLevel = stoi(tokens[1]);
                else if ((k == "orbitavgpatch.enabled" || k == "orbitavgpatch") && tokens.size() == 2)
                    parseBoolToken(tokens[1], config_.orbitAvgPatch.enabled);
                else if ((k == "orbitavgpatch.tol" || k == "orbitavgpatch.abstol") && tokens.size() == 2)
                    config_.orbitAvgPatch.tol = stod(tokens[1]);
                else if (k == "orbitavgpatch.maxrefinedepth" && tokens.size() == 2)
                    config_.orbitAvgPatch.maxRefineDepth = stoi(tokens[1]);
            }
            continue;
        }

        // Deposition block
        if (toLowerCopy(b) == "deposition")
        {
            if (tokens.size() == 2)
            {
                const string k = toLowerCopy(tokens[0]);
                if (k == "shape")
                    config_.deposition.shape = unquote(tokens[1]);
                else if (k == "ebregularization")
                    parseBoolToken(tokens[1], config_.deposition.ebRegularization);
                else if (k == "debugparticleid" || k == "deposit.debugparticleid")
                    config_.deposition.debugParticleId = stoi(tokens[1]);
                else if (k == "dumpweightsonce" || k == "deposit.dumpweightsonce")
                    parseBoolToken(tokens[1], config_.deposition.dumpWeightsOnce);
                else if (k == "dumpweightspath" || k == "deposit.dumpweightspath")
                    config_.deposition.dumpWeightsPath = unquote(tokens[1]);
            }
            continue;
        }

        // Output block
        if (toLowerCopy(b) == "output")
        {
            if (tokens.size() >= 2)
            {
                const string k = toLowerCopy(tokens[0]);
                if (k == "fields.every" || k == "output.fields.every")
                {
                    config_.output.every = stoi(tokens[1]);
                    sawOutputEvery = true;
                }
                // Deprecated: legacy toggles. They are accepted
                // for backward compatibility but ignored.
                else if (k == "writeparticles" || k == "writefields")
                {
                    (void)tokens;
                }
                else if (k == "fields" || k == "output.fields")
                {
                    // Parse a parenthesized list, allowing flexible spacing:
                    //   fields ( phi rho n Te E )
                    //   fields (phi rho)
                    //   fields phi rho
                    std::vector<std::string> out;
                    out.reserve(tokens.size());
                    for (size_t i = 1; i < tokens.size(); ++i)
                    {
                        std::string t = tokens[i];
                        // Strip leading '(' and trailing ')' and commas.
                        while (!t.empty() && (t.front() == '(' || t.front() == '['))
                            t.erase(t.begin());
                        while (!t.empty() && (t.back() == ')' || t.back() == ']' || t.back() == ','))
                            t.pop_back();
                        if (t.empty() || t == "(" || t == ")")
                            continue;
                        out.push_back(unquote(t));
                    }
                    config_.output.fields = std::move(out);
                }
                else if (k == "console.level" || k == "output.console.level")
                    config_.output.console.level = unquote(tokens[1]);
                else if (k == "console.every" || k == "output.console.every")
                    config_.output.console.every = stoi(tokens[1]);
                else if (k == "console.dqrelwarn" || k == "output.console.dqrelwarn")
                    config_.output.console.dQrelWarn = stod(tokens[1]);
                else if (k == "console.dqrelrecover" || k == "output.console.dqrelrecover")
                    config_.output.console.dQrelRecover = stod(tokens[1]);
                else if (k == "console.dqreldebouncesteps" || k == "output.console.dqreldebouncesteps")
                    config_.output.console.dQrelDebounceSteps = stoi(tokens[1]);
                else if (k == "console.dqreleventcooldownsteps" || k == "output.console.dqreleventcooldownsteps")
                    config_.output.console.dQrelEventCooldownSteps = stoi(tokens[1]);
                else if (k == "console.dqrelprintrecover" || k == "output.console.dqrelprintrecover")
                    parseBoolToken(tokens[1], config_.output.console.dQrelPrintRecover);
                else if (k == "console.physics.enabled" || k == "output.console.physics.enabled")
                    parseBoolToken(tokens[1], config_.output.console.physicsEnabled);
                else if (k == "console.physics.detail" || k == "output.console.physics.detail")
                    config_.output.console.physicsDetail = unquote(tokens[1]);
                else if (k == "console.physics.showcharge" || k == "output.console.physics.showcharge")
                    parseBoolToken(tokens[1], config_.output.console.physicsShowCharge);
                else if (k == "console.physics.showcurrent" || k == "output.console.physics.showcurrent")
                    parseBoolToken(tokens[1], config_.output.console.physicsShowCurrent);
                else if (k == "console.physics.showwall" || k == "output.console.physics.showwall")
                    parseBoolToken(tokens[1], config_.output.console.physicsShowWall);
                else if (k == "console.physics.showmcc" || k == "output.console.physics.showmcc")
                    parseBoolToken(tokens[1], config_.output.console.physicsShowMCC);
                else if (k == "console.physics.showperf" || k == "output.console.physics.showperf")
                    parseBoolToken(tokens[1], config_.output.console.physicsShowPerf);
                else if (k == "console.format" || k == "output.console.format")
                    config_.output.v2.console.format = unquote(tokens[1]);
                // Structured series outputs: preferred keys
                else if (k == "series.enabled" || k == "output.series.enabled")
                    parseBoolToken(tokens[1], config_.output.v2.enabled);
                else if (k == "series.every" || k == "output.series.every" ||
                         k == "series.seriesevery" || k == "output.series.seriesevery")
                {
                    config_.output.v2.seriesEvery = stoi(tokens[1]);
                }
                else if (k == "amount.enabled" || k == "output.amount.enabled")
                    parseBoolToken(tokens[1], config_.output.amount.enabled);
                else if (k == "periodicdiagnostics.enabled" || k == "output.periodicdiagnostics.enabled")
                {
                    sawPeriodicDiagnosticsEnabled = true;
                    parseBoolToken(tokens[1], config_.output.periodicDiagnostics.enabled);
                }
                else if (k == "periodicdiagnostics.console" || k == "output.periodicdiagnostics.console")
                {
                    if (!sawPeriodicDiagnosticsEnabled)
                        config_.output.periodicDiagnostics.enabled = true;
                    parseBoolToken(tokens[1], config_.output.periodicDiagnostics.printConsole);
                }

                // Output::CurrentWaveform flat keys (optional)
                else if (k == "currentwaveform.enabled" || k == "output.currentwaveform.enabled")
                {
                    sawCurrentWaveformEnabled = true;
                    parseBoolToken(tokens[1], config_.output.currentWaveform.enabled);
                }
                else if (k == "currentwaveform.side" || k == "currentwaveform.electrode" ||
                         k == "output.currentwaveform.side" || k == "output.currentwaveform.electrode")
                {
                    if (!sawCurrentWaveformEnabled)
                        config_.output.currentWaveform.enabled = true;
                    config_.output.currentWaveform.side = unquote(tokens[1]);
                }

                // Output::Profile1D flat keys (optional)
                else if (k == "profile1d.enabled" || k == "output.profile1d.enabled")
                {
                    sawProfile1DEnabled = true;
                    parseBoolToken(tokens[1], config_.output.profile1D.enabled);
                }
                else if (k == "profile1d.nbinsx" || k == "profile1d.nbins" || k == "output.profile1d.nbinsx")
                {
                    if (!sawProfile1DEnabled)
                        config_.output.profile1D.enabled = true;
                    config_.output.profile1D.nbinsX = stoi(tokens[1]);
                }

                // Output::Profile1DAvg flat keys (optional)
                else if (k == "profile1davg.enabled" || k == "output.profile1davg.enabled" ||
                         k == "profilesavg.enabled" || k == "output.profilesavg.enabled" ||
                         k == "profiles_avg.enabled" || k == "output.profiles_avg.enabled")
                {
                    sawProfile1DAvgEnabled = true;
                    parseBoolToken(tokens[1], config_.output.profile1DAvg.enabled);
                }
                else if (k == "profile1davg.avgwindowsteps" || k == "output.profile1davg.avgwindowsteps" ||
                         k == "profilesavg.avgwindowsteps" || k == "output.profilesavg.avgwindowsteps" ||
                         k == "profile1davg.window" || k == "profile1davg.windowsteps")
                {
                    if (!sawProfile1DAvgEnabled)
                        config_.output.profile1DAvg.enabled = true;
                    config_.output.profile1DAvg.avgWindowSteps = stoi(tokens[1]);
                }
                // Output::CurrentBalance flat keys (optional)
                else if (k == "currentbalance.enabled" || k == "output.currentbalance.enabled")
                {
                    sawCurrentBalanceEnabled = true;
                    parseBoolToken(tokens[1], config_.output.currentBalance.enabled);
                }

                // Output::EBBandRMS flat keys (optional)
                else if (k == "ebbandrms.enabled" || k == "output.ebbandrms.enabled" ||
                         k == "ebrms.enabled" || k == "output.ebrms.enabled")
                    parseBoolToken(tokens[1], config_.output.ebBandRMS.enabled);
                else if (k == "ebbandrms.bandlayers" || k == "output.ebbandrms.bandlayers" ||
                         k == "ebbandrms.layers" || k == "output.ebbandrms.layers")
                    config_.output.ebBandRMS.bandLayers = stoi(tokens[1]);
                else if (k == "diagnostics.enabled" || k == "output.diagnostics.enabled")
                {
                    sawDiagnosticsEnabled = true;
                    parseBoolToken(tokens[1], config_.output.diagnostics.enabled);
                }
                else if (k == "diagnostics.every" || k == "output.diagnostics.every")
                    config_.output.diagnostics.every = stoi(tokens[1]);
                else if (k == "diagnostics.basedir" || k == "output.diagnostics.basedir")
                    config_.output.diagnostics.baseDir = unquote(tokens[1]);
                else if (k == "diagnostics.electrodecurrent.enabled" || k == "output.diagnostics.electrodecurrent.enabled")
                {
                    sawElectrodeCurrentEnabled = true;
                    parseBoolToken(tokens[1], config_.output.diagnostics.electrodeCurrent.enabled);
                }
                else if (k == "diagnostics.reactions.enabled" || k == "output.diagnostics.reactions.enabled")
                {
                    sawReactionDiagnosticsEnabled = true;
                    parseBoolToken(tokens[1], config_.output.diagnostics.reactions.enabled);
                }
                else if (k == "diagnostics.energyhistogram.enabled" || k == "output.diagnostics.energyhistogram.enabled")
                {
                    sawEnergyHistogramEnabled = true;
                    parseBoolToken(tokens[1], config_.output.diagnostics.energyHistogram.enabled);
                }
            }
            continue;
        }

        // Output::Series nested block (optional)
        // Example:
        //   Output {
        //     Series {
        //       enabled true
        //       seriesEvery 50
        //     }
        //   }
        if (toLowerCopy(b) == "series" && toLowerCopy(parentBlock()) == "output")
        {
            if (tokens.size() >= 2)
            {
                const string k = toLowerCopy(tokens[0]);
                if (k == "enabled" || k == "enable")
                    parseBoolToken(tokens[1], config_.output.v2.enabled);
                else if (k == "every" || k == "seriesevery")
                {
                    config_.output.v2.seriesEvery = stoi(tokens[1]);
                }
            }
            continue;
        }

        if (toLowerCopy(b) == "fields" && toLowerCopy(parentBlock()) == "output")
        {
            if (tokens.size() >= 2)
            {
                const string k = toLowerCopy(tokens[0]);
                if (k == "every")
                    config_.output.every = stoi(tokens[1]);
                else if (k == "list" || k == "fields")
                {
                    std::vector<std::string> out;
                    out.reserve(tokens.size());
                    for (size_t i = 1; i < tokens.size(); ++i)
                    {
                        std::string t = unquote(tokens[i]);
                        while (!t.empty() && (t.front() == '(' || t.front() == ','))
                            t.erase(t.begin());
                        while (!t.empty() && (t.back() == ')' || t.back() == ','))
                            t.pop_back();
                        if (t.empty() || t == "(" || t == ")")
                            continue;
                        out.push_back(unquote(t));
                    }
                    config_.output.fields = std::move(out);
                }
            }
            continue;
        }

        if (toLowerCopy(b) == "amount" &&
            (toLowerCopy(parentBlock()) == "output" || toLowerCopy(parentBlock()) == "series"))
        {
            if (tokens.size() >= 2)
            {
                const string k = toLowerCopy(tokens[0]);
                if (k == "enabled" || k == "enable")
                    parseBoolToken(tokens[1], config_.output.amount.enabled);
            }
            continue;
        }

        // Output::Profile1D nested block (optional)
        // Example:
        //   Output {
        //     Profile1D {
        //       enabled true
        //       every 50
        //       nbinsX 64
        //     }
        //   }
        if (toLowerCopy(b) == "profile1d" &&
            (toLowerCopy(parentBlock()) == "output" || toLowerCopy(parentBlock()) == "fields"))
        {
            if (tokens.size() >= 2)
            {

                const string k = toLowerCopy(tokens[0]);
                if (k == "enabled" || k == "enable")
                {
                    sawProfile1DEnabled = true;
                    parseBoolToken(tokens[1], config_.output.profile1D.enabled);
                }
                else if (k == "nbinsx" || k == "nbins")
                {
                    if (!sawProfile1DEnabled)
                        config_.output.profile1D.enabled = true;
                    config_.output.profile1D.nbinsX = stoi(tokens[1]);
                }
            }
            continue;
        }

        // Output::Profile1DAvg nested block (optional)
        // Example:
        //   Output {
        //     Profile1DAvg {
        //       enabled true
        //       every 50
        //       avgWindowSteps 200
        //       path "profiles_avg.csv"
        //     }
        //   }
        if ((toLowerCopy(b) == "profile1davg" || toLowerCopy(b) == "profilesavg" || toLowerCopy(b) == "profiles_avg") &&
            (toLowerCopy(parentBlock()) == "output" || toLowerCopy(parentBlock()) == "fields"))
        {
            if (tokens.size() >= 2)
            {
                const string k = toLowerCopy(tokens[0]);
                if (k == "enabled" || k == "enable")
                {
                    sawProfile1DAvgEnabled = true;
                    parseBoolToken(tokens[1], config_.output.profile1DAvg.enabled);
                }
                else if (k == "avgwindowsteps" || k == "window" || k == "windowsteps")
                {
                    if (!sawProfile1DAvgEnabled)
                        config_.output.profile1DAvg.enabled = true;
                    config_.output.profile1DAvg.avgWindowSteps = stoi(tokens[1]);
                }
            }
            continue;
        }

        // Output::PeriodicDiagnostics nested block (optional)
        if (toLowerCopy(b) == "periodicdiagnostics" &&
            (toLowerCopy(parentBlock()) == "output" || toLowerCopy(parentBlock()) == "console"))
        {
            if (tokens.size() >= 2)
            {
                const string k = toLowerCopy(tokens[0]);
                if (k == "enabled" || k == "enable")
                {
                    sawPeriodicDiagnosticsEnabled = true;
                    parseBoolToken(tokens[1], config_.output.periodicDiagnostics.enabled);
                }
                else if (k == "console" || k == "printconsole")
                {
                    if (!sawPeriodicDiagnosticsEnabled)
                        config_.output.periodicDiagnostics.enabled = true;
                    parseBoolToken(tokens[1], config_.output.periodicDiagnostics.printConsole);
                }
            }
            continue;
        }

        // Output::CurrentBalance nested block (optional)
        if (toLowerCopy(b) == "currentbalance" &&
            (toLowerCopy(parentBlock()) == "output" || toLowerCopy(parentBlock()) == "series"))
        {
            if (tokens.size() >= 2)
            {
                const string k = toLowerCopy(tokens[0]);
                if (k == "enabled" || k == "enable")
                {
                    sawCurrentBalanceEnabled = true;
                    parseBoolToken(tokens[1], config_.output.currentBalance.enabled);
                }
            }
            continue;
        }

        // Output::CurrentWaveform nested block (optional)
        // Example:
        //   CurrentWaveform {
        //     enabled true
        //     every 20
        //     side auto
        //     path "current_waveform.csv"
        //   }
        if ((toLowerCopy(b) == "currentwaveform" || toLowerCopy(b) == "waveform") &&
            (toLowerCopy(parentBlock()) == "output" || toLowerCopy(parentBlock()) == "series"))
        {
            if (tokens.size() >= 2)
            {
                const string k = toLowerCopy(tokens[0]);
                if (k == "enabled" || k == "enable")
                {
                    sawCurrentWaveformEnabled = true;
                    parseBoolToken(tokens[1], config_.output.currentWaveform.enabled);
                }
                else if (k == "side" || k == "electrode")
                {
                    if (!sawCurrentWaveformEnabled)
                        config_.output.currentWaveform.enabled = true;
                    config_.output.currentWaveform.side = unquote(tokens[1]);
                }
            }
            continue;
        }

        // Output::EBBandRMS nested block (optional)
        if ((toLowerCopy(b) == "ebbandrms" || toLowerCopy(b) == "ebrms") &&
            (toLowerCopy(parentBlock()) == "output" || toLowerCopy(parentBlock()) == "series"))
        {
            if (tokens.size() >= 2)
            {
                const string k = toLowerCopy(tokens[0]);
                if (k == "enabled" || k == "enable")
                    parseBoolToken(tokens[1], config_.output.ebBandRMS.enabled);
                else if (k == "bandlayers" || k == "layers")
                    config_.output.ebBandRMS.bandLayers = stoi(tokens[1]);
            }
            continue;
        }

        if (toLowerCopy(b) == "diagnostics" &&
            (toLowerCopy(parentBlock()) == "output" || toLowerCopy(parentBlock()) == "series"))
        {
            if (tokens.size() >= 2)
            {
                const string k = toLowerCopy(tokens[0]);
                if (k == "enabled" || k == "enable")
                {
                    sawDiagnosticsEnabled = true;
                    parseBoolToken(tokens[1], config_.output.diagnostics.enabled);
                }
                else if (k == "every")
                    config_.output.diagnostics.every = stoi(tokens[1]);
                else if (k == "basedir" || k == "base" || k == "dir")
                    config_.output.diagnostics.baseDir = unquote(tokens[1]);
            }
            continue;
        }

        if (toLowerCopy(b) == "electrodecurrent" &&
            (toLowerCopy(parentBlock()) == "diagnostics" || toLowerCopy(parentBlock()) == "series"))
        {
            if (tokens.size() >= 2)
            {
                const string k = toLowerCopy(tokens[0]);
                if (k == "enabled" || k == "enable")
                {
                    sawElectrodeCurrentEnabled = true;
                    parseBoolToken(tokens[1], config_.output.diagnostics.electrodeCurrent.enabled);
                }
                else if (k == "every")
                    config_.output.diagnostics.electrodeCurrent.every = stoi(tokens[1]);
                else if (k == "perboundary")
                    parseBoolToken(tokens[1], config_.output.diagnostics.electrodeCurrent.perBoundary);
                else if (k == "includedomain")
                    parseBoolToken(tokens[1], config_.output.diagnostics.electrodeCurrent.includeDomain);
                else if (k == "group" && tokens.size() >= 3)
                {
                    OutputConfig::DiagnosticsConfig::ElectrodeGroup group;
                    group.name = unquote(tokens[1]);
                    parseBoundaryGroupTokens(tokens, 2, group.boundaryIds);
                    if (!group.name.empty() && !group.boundaryIds.empty())
                        config_.output.diagnostics.electrodeCurrent.groups.push_back(std::move(group));
                }
            }
            continue;
        }

        if ((toLowerCopy(b) == "reactions" || toLowerCopy(b) == "reactiondiagnostics") &&
            (toLowerCopy(parentBlock()) == "diagnostics" || toLowerCopy(parentBlock()) == "series"))
        {
            if (tokens.size() >= 2)
            {
                const string k = toLowerCopy(tokens[0]);
                if (k == "enabled" || k == "enable")
                {
                    sawReactionDiagnosticsEnabled = true;
                    parseBoolToken(tokens[1], config_.output.diagnostics.reactions.enabled);
                }
                else if (k == "every")
                    config_.output.diagnostics.reactions.every = stoi(tokens[1]);
                else if (k == "summary")
                    parseBoolToken(tokens[1], config_.output.diagnostics.reactions.summary);
                else if (k == "perchannel")
                    parseBoolToken(tokens[1], config_.output.diagnostics.reactions.perChannel);
                else if (k == "perprojectile")
                    parseBoolToken(tokens[1], config_.output.diagnostics.reactions.perProjectile);
                else if (k == "spawnedspecies")
                    parseBoolToken(tokens[1], config_.output.diagnostics.reactions.spawnedSpecies);
            }
            continue;
        }

        if ((toLowerCopy(b) == "energyhistogram" || toLowerCopy(b) == "energy") &&
            (toLowerCopy(parentBlock()) == "diagnostics" || toLowerCopy(parentBlock()) == "series"))
        {
            if (tokens.size() >= 2)
            {
                const string k = toLowerCopy(tokens[0]);
                if (k == "enabled" || k == "enable")
                {
                    sawEnergyHistogramEnabled = true;
                    parseBoolToken(tokens[1], config_.output.diagnostics.energyHistogram.enabled);
                }
                else if (k == "every")
                    config_.output.diagnostics.energyHistogram.every = stoi(tokens[1]);
                else if (k == "bins")
                    config_.output.diagnostics.energyHistogram.bins = stoi(tokens[1]);
                else if (k == "maxenergyev" || k == "maxe_ev")
                    config_.output.diagnostics.energyHistogram.maxEnergyEV = stod(tokens[1]);
                else if (k == "sample" && tokens.size() >= 4)
                {
                    OutputConfig::DiagnosticsConfig::EnergySample sample;
                    sample.name = unquote(tokens[1]);
                    sample.species = unquote(tokens[2]);
                    sample.location = unquote(tokens[3]);
                    if (!sample.name.empty())
                        config_.output.diagnostics.energyHistogram.samples.push_back(std::move(sample));
                }
            }
            continue;
        }

        if ((toLowerCopy(b) == "particletrajectory" || toLowerCopy(b) == "particletracks" || toLowerCopy(b) == "trajectories") &&
            (toLowerCopy(parentBlock()) == "output" || toLowerCopy(parentBlock()) == "series" || toLowerCopy(parentBlock()) == "diagnostics"))
        {
            if (tokens.size() >= 2)
            {
                const string k = toLowerCopy(tokens[0]);
                if (k == "enabled" || k == "enable")
                    parseBoolToken(tokens[1], config_.output.particleTrajectory.enabled);
                else if (k == "every")
                    config_.output.particleTrajectory.every = stoi(tokens[1]);
                else if (k == "electronevery" || k == "electroninterval" || k == "electronoutputevery")
                    config_.output.particleTrajectory.electronEvery = stoi(tokens[1]);
                else if (k == "ionevery" || k == "ioninterval" || k == "ionoutputevery")
                    config_.output.particleTrajectory.ionEvery = stoi(tokens[1]);
                else if (k == "electronsampleevery" || k == "electronsampleinterval" || k == "electrondumpevery")
                    config_.output.particleTrajectory.electronSampleEvery = stoi(tokens[1]);
                else if (k == "ionsampleevery" || k == "ionsampleinterval" || k == "iondumpevery")
                    config_.output.particleTrajectory.ionSampleEvery = stoi(tokens[1]);
                else if (k == "electroncount" || k == "electrons" || k == "electrontracks")
                    config_.output.particleTrajectory.electronCount = stoi(tokens[1]);
                else if (k == "ioncount" || k == "ions" || k == "iontracks")
                    config_.output.particleTrajectory.ionCount = stoi(tokens[1]);
                else if (k == "electronmaxsamples" || k == "electronsteps" || k == "electronrecords")
                    config_.output.particleTrajectory.electronMaxSamples = stoi(tokens[1]);
                else if (k == "ionmaxsamples" || k == "ionsteps" || k == "ionrecords")
                    config_.output.particleTrajectory.ionMaxSamples = stoi(tokens[1]);
                else if (k == "basedir" || k == "base" || k == "dir")
                    config_.output.particleTrajectory.baseDir = unquote(tokens[1]);
            }
            continue;
        }

        // Checkpoint / restart block (optional)
        if (toLowerCopy(b) == "checkpoint" || toLowerCopy(b) == "cache" || toLowerCopy(b) == "restart")
        {
            if (tokens.size() >= 2)
            {
                const string k = toLowerCopy(tokens[0]);
                if (k == "enabled" || k == "enable")
                    parseBoolToken(tokens[1], config_.checkpoint.enabled);
                else if (k == "every" || k == "everysteps")
                    config_.checkpoint.every = stoi(tokens[1]);
                else if (k == "keep" || k == "keeplast" || k == "keep_last")
                    config_.checkpoint.keep = stoi(tokens[1]);
                else if (k == "path" || k == "dir" || k == "cachedir" || k == "folder")
                    config_.checkpoint.path = unquote(tokens[1]);
                else if (k == "restart" || k == "resume")
                    parseBoolToken(tokens[1], config_.checkpoint.restart);
                else if (k == "restartstep" || k == "step" || k == "resumestep")
                    config_.checkpoint.restartStep = stoi(tokens[1]);
                else if (k == "saverng" || k == "rng" || k == "rngstate")
                    parseBoolToken(tokens[1], config_.checkpoint.saveRng);
            }
            continue;
        }

        // Output::DebugFields nested block removed.

        // Output::Console nested block (optional)
        if (toLowerCopy(b) == "console")
        {
            if (tokens.size() >= 2)
            {
                const string k = toLowerCopy(tokens[0]);
                if (k == "level" || k == "loglevel" || k == "verbosity")
                    config_.output.console.level = unquote(tokens[1]);
                else if (k == "every")
                    config_.output.console.every = stoi(tokens[1]);
                else if (k == "dqrelwarn")
                    config_.output.console.dQrelWarn = stod(tokens[1]);
                else if (k == "dqrelrecover")
                    config_.output.console.dQrelRecover = stod(tokens[1]);
                else if (k == "dqreldebouncesteps")
                    config_.output.console.dQrelDebounceSteps = stoi(tokens[1]);
                else if (k == "dqreleventcooldownsteps")
                    config_.output.console.dQrelEventCooldownSteps = stoi(tokens[1]);
                else if (k == "dqrelprintrecover")
                    parseBoolToken(tokens[1], config_.output.console.dQrelPrintRecover);
                else if (k == "physics.enabled")
                    parseBoolToken(tokens[1], config_.output.console.physicsEnabled);
                else if (k == "physics.detail")
                    config_.output.console.physicsDetail = unquote(tokens[1]);
                else if (k == "physics.showcharge")
                    parseBoolToken(tokens[1], config_.output.console.physicsShowCharge);
                else if (k == "physics.showcurrent")
                    parseBoolToken(tokens[1], config_.output.console.physicsShowCurrent);
                else if (k == "physics.showwall")
                    parseBoolToken(tokens[1], config_.output.console.physicsShowWall);
                else if (k == "physics.showmcc")
                    parseBoolToken(tokens[1], config_.output.console.physicsShowMCC);
                else if (k == "physics.showperf")
                    parseBoolToken(tokens[1], config_.output.console.physicsShowPerf);
                else if (k == "format")
                    config_.output.v2.console.format = unquote(tokens[1]);
            }
            continue;
        }

        // Output::Console::Physics nested block (optional)
        if (toLowerCopy(b) == "physics" && toLowerCopy(parentBlock()) == "console")
        {
            if (tokens.size() >= 2)
            {
                const string k = toLowerCopy(tokens[0]);
                if (k == "enabled" || k == "enable")
                    parseBoolToken(tokens[1], config_.output.console.physicsEnabled);
                else if (k == "detail")
                    config_.output.console.physicsDetail = unquote(tokens[1]);
                else if (k == "showcharge")
                    parseBoolToken(tokens[1], config_.output.console.physicsShowCharge);
                else if (k == "showcurrent")
                    parseBoolToken(tokens[1], config_.output.console.physicsShowCurrent);
                else if (k == "showwall")
                    parseBoolToken(tokens[1], config_.output.console.physicsShowWall);
                else if (k == "showmcc")
                    parseBoolToken(tokens[1], config_.output.console.physicsShowMCC);
                else if (k == "showperf")
                    parseBoolToken(tokens[1], config_.output.console.physicsShowPerf);
            }
            continue;
        }

        // Secondary electron emission (SEE) block (optional)
        // Example:
        //   SEE {
        //     enabled true
        //     cathodeSides left right
        //     emittedSpecies electron
        //     gamma 0.1
        //     emittedEnergy 2 eV
        //     injectEps 1e-6 m
        //     maxEmitMacro 100
        //   }
        if (toLowerCopy(b) == "see")
        {
            SEEConfig* seeCfg = &config_.see;
            if (toLowerCopy(parentBlock()) == "particleboundaryitem" && activeParticleBoundaryId >= 0)
            {
                auto& entry = getOrCreateParticleBC(activeParticleBoundaryId);
                entry.seeSpecified = true;
                seeCfg = &entry.see;
            }

            if (tokens.size() >= 2)
            {
                const string k = toLowerCopy(tokens[0]);
                if (k == "enabled")
                {
                    parseBoolToken(tokens[1], seeCfg->enabled);
                }
                else if (k == "cathodeside" || k == "cathodesides")
                {
                    seeCfg->cathodeSides.clear();
                    for (size_t i = 1; i < tokens.size(); ++i)
                        appendSeeSideToken(tokens[i], seeCfg->cathodeSides);
                    if (!seeCfg->cathodeSides.empty())
                        seeCfg->cathodeSide = seeCfg->cathodeSides.front();
                    else if (tokens.size() >= 2)
                        seeCfg->cathodeSide = unquote(tokens[1]);
                }
                else if (k == "emittedspecies")
                {
                    seeCfg->emittedSpecies = unquote(tokens[1]);
                }
                else if (k == "gamma")
                {
                    seeCfg->gamma = stod(tokens[1]);
                }
                else if (k == "gammamodel" || k == "gamma_model" || k == "model")
                {
                    seeCfg->gammaModel = unquote(tokens[1]);
                }
                else if (k == "gamma0")
                {
                    seeCfg->gamma0 = stod(tokens[1]);
                }
                else if (k == "eref" || k == "eref_ev" || k == "eref_eV" || k == "erefenergy" || k == "erefenergy_ev")
                {
                    double e_eV = 0.0;
                    if (parseEnergyToEV(tokens, 1, e_eV))
                        seeCfg->Eref_eV = e_eV;
                }
                else if (k == "alpha")
                {
                    seeCfg->alpha = stod(tokens[1]);
                }
                else if (k == "gammamax" || k == "gamma_max")
                {
                    seeCfg->gammaMax = stod(tokens[1]);
                }
                else if (k == "emin" || k == "emin_ev" || k == "emin_eV")
                {
                    double e_eV = 0.0;
                    if (parseEnergyToEV(tokens, 1, e_eV))
                        seeCfg->Emin_eV = e_eV;
                }
                else if (k == "piecewise" || k == "gammapiecewise" || k == "gammatable" || k == "table")
                {
                    // Format: piecewise E1 gamma1 E2 gamma2 ... (E in eV)
                    // Example: piecewise 0 0.0 50 0.05 100 0.10
                    seeCfg->piecewiseE_eV.clear();
                    seeCfg->piecewiseGamma.clear();
                    const size_t n = tokens.size();
                    if (n >= 5 && ((n - 1) % 2 == 0))
                    {
                        for (size_t i = 1; i + 1 < n; i += 2)
                        {
                            double E = 0.0, g = 0.0;
                            if (!tryParseStrictDouble(tokens[i], E) || !tryParseStrictDouble(tokens[i + 1], g))
                                continue;
                            seeCfg->piecewiseE_eV.push_back(E);
                            seeCfg->piecewiseGamma.push_back(g);
                        }
                    }
                }
                else if (k == "emittedenergy" || k == "emittedenergy_ev" || k == "emittedenergy_eV")
                {
                    double e_eV = 0.0;
                    if (parseEnergyToEV(tokens, 1, e_eV))
                        seeCfg->emittedEnergy_eV = e_eV;
                }
                else if (k == "emittedenergymodel" || k == "emitted_energy_model" || k == "energymodel" || k == "emittedmodel")
                {
                    seeCfg->emittedEnergyModel = unquote(tokens[1]);
                }
                else if (k == "emittedte" || k == "emittedte_ev" || k == "emittedte_eV" || k == "emittedtemperature" || k == "te_emit" || k == "teemit")
                {
                    double te_eV = 0.0;
                    if (parseEnergyToEV(tokens, 1, te_eV))
                        seeCfg->emittedTe_eV = te_eV;
                }
                else if (k == "emittedemean" || k == "emittedmean" || k == "emean" || k == "emean_ev" || k == "emean_eV" || k == "emittedemean_ev" || k == "emittedemean_eV")
                {
                    double em_eV = 0.0;
                    if (parseEnergyToEV(tokens, 1, em_eV))
                        seeCfg->emittedEmean_eV = em_eV;
                }
                else if (k == "injecteps" || k == "injecteps_m")
                {
                    double eps_si = 0.0;
                    if (parseScalarToSI(tokens, 1, eps_si))
                        seeCfg->injectEps_m = eps_si;
                }
                else if (k == "maxemitmacro")
                {
                    seeCfg->maxEmitMacro = stoi(tokens[1]);
                }
                else if (k == "angular" || k == "angularmodel" || k == "distribution")
                {
                    seeCfg->angular = unquote(tokens[1]);
                }
            }
            continue;
        }

        // Domain BC block
        if (toLowerCopy(b) == "domain")
        {
            if (tokens.size() >= 2)
            {
                const string key = toLowerCopy(tokens[0]);
                const size_t dot = key.find('.');
                if (dot != string::npos)
                {
                    const string side = key.substr(0, dot);
                    const string field = key.substr(dot + 1);
                    if (field == "waveform" || field == "time" || field == "timedirichlet" || field == "rf")
                    {
                        PicBoundarySide* bc = getDomainSide(side);
                        if (bc)
                        {
                            parseTimeDirichletTokens(tokens, 1, units, bc->timeDirichlet, bc->value);
                        }
                        continue;
                    }
                }
            }

            if (tokens.size() >= 3)
            {
                bool ok = false;
                PicBCType type = parseBCTypeToken(tokens[1], ok);
                if (ok)
                {
                    double value_code = 0.0;
                    if (!parsePotentialToCode(tokens, 2, units, value_code))
                    {
                        std::cerr << "[DictionaryPIC] ERROR: failed to parse domain BC value: '" << line << "'" << std::endl;
                        exit(1);
                    }
                    setDomainSide(tokens[0], type, value_code);
                }
            }
            continue;
        }

        // Field BC by boundaryId (BC/Field/boundary <id> { ... })
        if (toLowerCopy(b) == "fieldboundaryitem")
        {
            if (activeFieldBoundaryId >= 0)
            {
                auto& entry = getOrCreateFieldBC(activeFieldBoundaryId);
                if (tokens.size() >= 2)
                {
                    const string k = toLowerCopy(tokens[0]);
                    if (k == "type")
                    {
                        bool ok = false;
                        entry.bc.type = parseBCTypeToken(tokens[1], ok);
                    }
                    else if (k == "value")
                    {
                        double value_code = 0.0;
                        if (!parsePotentialToCode(tokens, 1, units, value_code))
                        {
                            std::cerr << "[DictionaryPIC] ERROR: failed to parse boundary BC value: '" << line << "'" << std::endl;
                            exit(1);
                        }
                        entry.bc.value = value_code;
                    }
                    else if (k == "waveform" || k == "time" || k == "timedirichlet" || k == "rf")
                    {
                        parseTimeDirichletTokens(tokens, 1, units, entry.bc.timeDirichlet, entry.bc.value);
                    }
                }
            }
            continue;
        }

        // Particle BC by boundaryId (BC/Particle/boundary <id> { ... })
        if (toLowerCopy(b) == "particleboundaryitem")
        {
            if (activeParticleBoundaryId >= 0)
            {
                auto& entry = getOrCreateParticleBC(activeParticleBoundaryId);
                if (tokens.size() >= 2)
                {
                    const string k = toLowerCopy(tokens[0]);
                    if (k == "model")
                        entry.model = unquote(tokens[1]);
                }
            }
            continue;
        }

        // AMR block (PIC-side policy)
        if (toLowerCopy(b) == "amr")
        {
            if (tokens.size() >= 2)
            {
                const string k = toLowerCopy(tokens[0]);
                if (k == "adapteverysteps" || k == "adaptevery" || k == "every")
                    config_.amr.adaptEverySteps = stoi(tokens[1]);
                else if (k == "indicatoreverysteps" || k == "indicatorevery" || k == "indicatoreverystep")
                    config_.amr.indicatorEverySteps = stoi(tokens[1]);
                else if (k == "persistencesteps" || k == "persistence")
                    config_.amr.persistenceSteps = stoi(tokens[1]);
                else if (k == "cellbudgetperadapt" || k == "cellbudget")
                    config_.amr.cellBudgetPerAdapt = stoi(tokens[1]);
                else if (k == "hysteresis")
                    config_.amr.hysteresis = stod(tokens[1]);
                else if (k == "geometryrefreshmode" || k == "geomrefreshmode")
                {
                    config_.amr.geometryRefresh.present = true;
                    config_.amr.geometryRefresh.mode = unquote(tokens[1]);
                }
                else if (k == "geometryrefreshchangerel" || k == "geomrefreshchangerel")
                {
                    config_.amr.geometryRefresh.present = true;
                    config_.amr.geometryRefresh.changeRel = stod(tokens[1]);
                }
                else if (k == "geometryrefreshverbose" || k == "geomrefreshverbose")
                {
                    config_.amr.geometryRefresh.present = true;
                    parseBoolToken(tokens[1], config_.amr.geometryRefresh.verbose);
                }
                else if (k == "geometryrefreshallowwithmovingeb" || k == "geomrefreshallowwithmovingeb")
                {
                    config_.amr.geometryRefresh.present = true;
                    parseBoolToken(tokens[1], config_.amr.geometryRefresh.allowWithMovingEb);
                }
                else if (k == "movingebrefreshmode")
                {
                    config_.amr.geometryRefresh.present = true;
                    config_.amr.geometryRefresh.movingEbRefreshMode = unquote(tokens[1]);
                }
            }
            continue;
        }

        if (toLowerCopy(b) == "geometryrefresh" && toLowerCopy(parentBlock()) == "amr")
        {
            if (tokens.size() >= 2)
            {
                const string k = toLowerCopy(tokens[0]);
                config_.amr.geometryRefresh.present = true;
                if (k == "mode")
                    config_.amr.geometryRefresh.mode = unquote(tokens[1]);
                else if (k == "changerel" || k == "relchange" || k == "relth")
                    config_.amr.geometryRefresh.changeRel = stod(tokens[1]);
                else if (k == "verbose")
                    parseBoolToken(tokens[1], config_.amr.geometryRefresh.verbose);
                else if (k == "allowwithmovingeb" || k == "allowwithmovingebrefresh")
                    parseBoolToken(tokens[1], config_.amr.geometryRefresh.allowWithMovingEb);
                else if (k == "movingebrefreshmode" || k == "movingebmode")
                    config_.amr.geometryRefresh.movingEbRefreshMode = unquote(tokens[1]);
            }
            continue;
        }

        if (toLowerCopy(b) == "debye")
        {
            if (tokens.size() >= 2)
            {
                const string k = toLowerCopy(tokens[0]);
                if (k == "enabled")
                    parseBoolToken(tokens[1], config_.amr.debye.enabled);
                if (k == "cd")
                    config_.amr.debye.CD = stod(tokens[1]);
                else if (k == "te" || k == "te_ev")
                {
                    double e_eV = 0.0;
                    if (parseEnergyToEV(tokens, 1, e_eV))
                        config_.amr.debye.Te_eV = e_eV;
                }
                else if (k == "nefloor" || k == "ne_floor")
                {
                    double v = 0.0;
                    if (parseNumberDensityToSI(tokens, 1, v))
                        config_.amr.debye.neFloor = v;
                }
                else if (k == "avgwindowsteps" || k == "avgwindow")
                    config_.amr.debye.avgWindowSteps = stoi(tokens[1]);
                else if (k == "smoothwindowcells" || k == "smoothwindowradius" ||
                         k == "spatialsmoothcells" || k == "spatialwindowcells")
                    config_.amr.debye.smoothWindowCells = stoi(tokens[1]);
            }
            continue;
        }

        if (toLowerCopy(b) == "indicator")
        {
            if (tokens.size() >= 2)
            {
                const string k = toLowerCopy(tokens[0]);
                if (k == "useebband" ||
                    k == "ebbandwidth" || k == "eb_bandwidth" || k == "ebbandwidthcells" ||
                    k == "ebbandcoeff" || k == "eb_bandcoeff" || k == "ebbandc" || k == "eb_bandc" ||
                    k == "ebbandtargetlevel" || k == "eb_bandtargetlevel" || k == "ebbandlevel" || k == "eb_bandlevel" ||
                    k == "ebbandboost" || k == "eb_bandboost" ||
                    k == "ebbandminlevel" || k == "eb_bandminlevel")
                {
                    std::cerr << "[DictionaryPIC] WARNING: EBBand refinement is no longer supported; "
                                 "remove '" << tokens[0] << "' from AMR/Indicator.\n";
                }
                else if (k == "useefield")
                    parseBoolToken(tokens[1], config_.amr.indicator.useEField);
                else if (k == "efieldpercentile" || k == "percentile")
                    config_.amr.indicator.eFieldPercentile = stod(tokens[1]);
                else if (k == "efieldavgwindow" || k == "efieldavgwindowsteps")
                    config_.amr.indicator.eFieldAvgWindow = stoi(tokens[1]);
                else if (k == "smoothwindowcells" || k == "smoothwindowradius" ||
                         k == "efieldsmoothwindowcells" || k == "spatialsmoothcells" ||
                         k == "spatialwindowcells")
                    config_.amr.indicator.eFieldSmoothWindowCells = stoi(tokens[1]);
                else if (k == "efieldabsfloor" || k == "absfloor" || k == "efieldabs")
                {
                    double v = 0.0;
                    if (parseScalarToSI(tokens, 1, v))
                        config_.amr.indicator.eFieldAbsFloor = v;
                }
                else if (k == "efieldreffrac" || k == "reffrac")
                    config_.amr.indicator.eFieldRefFrac = stod(tokens[1]);
                else if (k == "efieldrefavgwindow" || k == "efieldrefavgwindowsteps" ||
                         k == "refavgwindow" || k == "refavgwindowsteps")
                    config_.amr.indicator.eFieldRefAvgWindow = stoi(tokens[1]);
            }
            continue;
        }

        // V3 RefinementIndicators block
        if (toLowerCopy(b) == "debyelength" && toLowerCopy(parentBlock()) == "refinementindicators")
        {
            if (tokens.size() >= 2)
            {
                const string k = toLowerCopy(tokens[0]);
                if (k == "enabled")
                    parseBoolToken(tokens[1], config_.amr.debye.enabled);
                else if (k == "cd")
                    config_.amr.debye.CD = stod(tokens[1]);
                else if (k == "te" || k == "te_ev")
                {
                    double e_eV = 0.0;
                    if (parseEnergyToEV(tokens, 1, e_eV))
                        config_.amr.debye.Te_eV = e_eV;
                }
                else if (k == "nefloor" || k == "ne_floor")
                {
                    double v = 0.0;
                    if (parseNumberDensityToSI(tokens, 1, v))
                        config_.amr.debye.neFloor = v;
                }
                else if (k == "avgwindowsteps" || k == "avgwindow")
                    config_.amr.debye.avgWindowSteps = stoi(tokens[1]);
                else if (k == "smoothwindowcells" || k == "smoothwindowradius" ||
                         k == "spatialsmoothcells" || k == "spatialwindowcells")
                    config_.amr.debye.smoothWindowCells = stoi(tokens[1]);
            }
            continue;
        }

        if (toLowerCopy(b) == "electricfield" && toLowerCopy(parentBlock()) == "refinementindicators")
        {
            if (tokens.size() >= 2)
            {
                const string k = toLowerCopy(tokens[0]);
                if (k == "enabled")
                    parseBoolToken(tokens[1], config_.amr.indicator.useEField);
                else if (k == "percentile")
                    config_.amr.indicator.eFieldPercentile = stod(tokens[1]);
                else if (k == "smoothwindowsteps" || k == "avgwindowsteps" || k == "efieldavgwindow" ||
                         k == "efieldavgwindowsteps")
                    config_.amr.indicator.eFieldAvgWindow = stoi(tokens[1]);
                else if (k == "smoothwindowcells" || k == "smoothwindowradius" ||
                         k == "efieldsmoothwindowcells" || k == "spatialsmoothcells" ||
                         k == "spatialwindowcells")
                    config_.amr.indicator.eFieldSmoothWindowCells = stoi(tokens[1]);
                else if (k == "absfloor" || k == "efieldabsfloor" || k == "efieldabs")
                {
                    double v = 0.0;
                    if (parseScalarToSI(tokens, 1, v))
                        config_.amr.indicator.eFieldAbsFloor = v;
                }
                else if (k == "reffrac" || k == "efieldreffrac")
                    config_.amr.indicator.eFieldRefFrac = stod(tokens[1]);
                else if (k == "refavgwindow" || k == "refavgwindowsteps" ||
                         k == "efieldrefavgwindow" || k == "efieldrefavgwindowsteps")
                    config_.amr.indicator.eFieldRefAvgWindow = stoi(tokens[1]);
            }
            continue;
        }

        if (toLowerCopy(b) == "particlerelabel")
        {
            if (tokens.size() >= 2)
            {
                const string k = toLowerCopy(tokens[0]);
                if (k == "mode")
                    config_.amr.particleRelabel.mode = unquote(tokens[1]);
                else if (k == "fallback")
                    config_.amr.particleRelabel.fallback = unquote(tokens[1]);
            }
            continue;
        }

        // ParticleManagement block (optional, simplified)
        if (toLowerCopy(b) == "particlemanagement")
        {
            if (tokens.size() >= 2)
            {
                auto& pm = config_.particleManagement;
                const string k = toLowerCopy(tokens[0]);
                if (k == "enabled")
                {
                    parseBoolToken(tokens[1], pm.enabled);
                }
                else if (k == "every" || k == "everysteps")
                {
                    pm.every = stoi(tokens[1]);
                }
                else if (k == "softcaptotal" || k == "softcap")
                {
                    pm.softCapTotal = stoll(tokens[1]);
                }
                else if (k == "hardcaptotal" || k == "hardcap")
                {
                    pm.hardCapTotal = stoll(tokens[1]);
                }
                else if (k == "percellmax" || k == "cellmax")
                {
                    pm.perCellMax = stoi(tokens[1]);
                }
                else if (k == "percellmin" || k == "cellmin")
                {
                    pm.perCellMin = stoi(tokens[1]);
                }
                else if (k == "merge" || k == "mergeenabled")
                {
                    parseBoolToken(tokens[1], pm.mergeEnabled);
                }
                else if (k == "mergemaxdeltae" || k == "mergemaxdeltae_ev" || k == "merge.maxdeltae")
                {
                    double e_eV = 0.0;
                    if (parseEnergyToEV(tokens, 1, e_eV))
                        pm.mergeMaxDeltaE_eV = e_eV;
                }
                else if (k == "mergemaxdeltaebyspecies" || k == "mergemaxdeltae_list" ||
                         k == "merge.maxdeltae.byspecies" || k == "merge.maxdeltae.list")
                {
                    std::vector<std::string> listTokens;
                    listTokens.reserve(tokens.size());
                    for (size_t i = 1; i < tokens.size(); ++i)
                    {
                        std::string t = stripListPunct(unquote(tokens[i]));
                        if (t.empty())
                            continue;
                        listTokens.push_back(t);
                    }

                    pm.mergeMaxDeltaE_bySpecies_eV.clear();
                    size_t idx = 0;
                    while (idx < listTokens.size())
                    {
                        double e_eV = 0.0;
                        size_t consumed = 0;
                        if (parseEnergyToEVAt(listTokens, idx, e_eV, consumed))
                        {
                            pm.mergeMaxDeltaE_bySpecies_eV.push_back(e_eV);
                            idx += consumed;
                        }
                        else
                        {
                            // Skip unrecognized token.
                            idx += 1;
                        }
                    }
                }
                else if (k == "mergemaxpairspercell" || k == "mergemaxpairs" || k == "merge.maxpairs")
                {
                    pm.mergeMaxPairsPerCell = stoi(tokens[1]);
                }
                else if (k == "split" || k == "splitenabled")
                {
                    parseBoolToken(tokens[1], pm.splitEnabled);
                }
                else if (k == "splitbycount" || k == "splitcount" || k == "split.count")
                {
                    parseBoolToken(tokens[1], pm.splitByCountEnabled);
                }
                else if (k == "splitwmax" || k == "wmax" || k == "split.maxweight")
                {
                    pm.splitWMax = stod(tokens[1]);
                }
                else if (k == "splitmaxpercell" || k == "splitmax" || k == "split.maxpercell")
                {
                    pm.splitMaxPerCell = stoi(tokens[1]);
                }
                else if (k == "splitjitterfrac" || k == "splitjitter" || k == "split.jitter")
                {
                    pm.splitJitterFrac = stod(tokens[1]);
                }
                else if (k == "splitebsafe" || k == "split.ebsafe")
                {
                    parseBoolToken(tokens[1], pm.splitEbSafe);
                }
                else if (k == "splitmaxretries" || k == "split.retries")
                {
                    pm.splitMaxRetries = stoi(tokens[1]);
                }
                else if (k == "thinning" || k == "thinningenabled")
                {
                    parseBoolToken(tokens[1], pm.thinningEnabled);
                }
                else if (k == "thinningmccspawn" || k == "thinningmcc" || k == "thinning.mcc")
                {
                    parseBoolToken(tokens[1], pm.thinningMccSpawn);
                }
                else if (k == "thinningseespawn" || k == "thinningsee" || k == "thinning.see")
                {
                    parseBoolToken(tokens[1], pm.thinningSeeSpawn);
                }
                else if (k == "minkeepprob" || k == "minkeepprobability")
                {
                    pm.minKeepProb = stod(tokens[1]);
                }
                else if (k == "report" || k == "reportenabled")
                {
                    parseBoolToken(tokens[1], pm.reportEnabled);
                }
                else if (k == "reportevery" || k == "reporteverysteps")
                {
                    pm.reportEvery = stoi(tokens[1]);
                }
            }
            continue;
        }

        // Embedded BC block (optional)
        if (toLowerCopy(b) == "embedded")
        {
            if (tokens.size() >= 2)
            {
                const string k = toLowerCopy(tokens[0]);
                if (k == "enabled")
                {
                    parseBoolToken(tokens[1], config_.embeddedBC.enabled);
                }
                else if (k == "type")
                {
                    bool ok = false;
                    config_.embeddedBC.type = parseBCTypeToken(tokens[1], ok);
                    config_.embeddedBC.enabled = ok;
                }
                else if (k == "a")
                {
                    config_.embeddedBC.a = stod(tokens[1]);
                }
                else if (k == "b")
                {
                    config_.embeddedBC.b = stod(tokens[1]);
                }
                else if (k == "g")
                {
                    config_.embeddedBC.g = stod(tokens[1]);
                }
                else if (k == "waveform" || k == "time" || k == "timedirichlet" || k == "rf")
                {
                    parseTimeDirichletTokens(tokens, 1, units, config_.embeddedBC.timeDirichlet, config_.embeddedBC.g);
                }
            }
            continue;
        }

        // Embedded geometry block (optional)
        if (toLowerCopy(b) == "embeddedgeometry")
        {
            if (!warnedDeprecatedEmbeddedGeometryBlock)
            {
                std::cerr << "[DictionaryPIC] WARNING: pic.dict/EmbeddedGeometry is deprecated and ignored. "
                             "Embedded boundary activation now auto-detects immersed non-domain boundaries from geometry.dict.\n";
                warnedDeprecatedEmbeddedGeometryBlock = true;
            }
            continue;
        }

        // Species item block
        if (toLowerCopy(b) == "speciesitem" && activeSpeciesIndex >= 0)
        {
            SpeciesConfig& sp = config_.species[activeSpeciesIndex];

            // Vector2D values
            if (toLowerCopy(tokens[0]) == "drift")
            {
                Vector2D d;
                // Drift is stored in SI (m/s). Accept both legacy bare numbers and unit-bearing inputs.
                if (parseVector2ToSIWithUnits(line, d) || parseCoord(line, d))
                    sp.drift = d;
                continue;
            }
            if (toLowerCopy(tokens[0]) == "sourcedrift")
            {
                Vector2D d;
                if (parseVector2ToSIWithUnits(line, d) || parseCoord(line, d))
                {
                    sp.source.drift = d;
                    sp.source.driftSet = true;
                }
                continue;
            }

            if (tokens.size() >= 2)
            {
                const string k = toLowerCopy(tokens[0]);
                if (k == "q" || k == "charge")
                {
                    double q_si = 0.0;
                    if (parseScalarToSI(tokens, 1, q_si))
                        sp.q = q_si;
                }
                else if (k == "m" || k == "mass")
                {
                    // Mass is stored in SI kg.
                    // Support common plasma shorthands explicitly to avoid silent fallback.
                    // Prefer parsing from the raw RHS text to be resilient to odd whitespace.
                    {
                        const size_t pos = line.find(tokens[0]);
                        if (pos != std::string::npos)
                        {
                            const std::string rhs = trim_ws(line.substr(pos + tokens[0].size()));
                            UnitParser::Result rr;
                            if (UnitParser::TryParseToSI(rhs, rr))
                            {
                                sp.m = rr.value_si;
                                continue;
                            }
                        }
                    }

                    if (tokens.size() >= 3)
                    {
                        const std::string unit = stripListPunct(lowerCopy(tokens[2]));
                        double numeric = 0.0;
                        if (tryParseStrictDouble(tokens[1], numeric))
                        {
                            if (unit == "me")
                            {
                                sp.m = numeric * PhysicalConstants::me_kg;
                                continue;
                            }
                            if (unit == "mp")
                            {
                                sp.m = numeric * PhysicalConstants::mp_kg;
                                continue;
                            }
                        }
                    }

                    double m_si = 0.0;
                    if (parseScalarToSI(tokens, 1, m_si))
                    {
                        sp.m = m_si;
                    }
                    else
                    {
                        std::cerr << "[DictionaryPIC] WARNING: failed to parse species mass (using default): '" << line << "'" << std::endl;
                    }
                }
                else if (k == "weight" || k == "wp")
                    sp.weight = stod(tokens[1]);
                else if (k == "numparticles" || k == "np")
                    sp.numParticles = stoi(tokens[1]);
                else if (k == "energycutoff" || k == "energycutoff_ev" || k == "energycutoffev")
                {
                    double cutoff_eV = 0.0;
                    if (parseEnergyToEV(tokens, 1, cutoff_eV))
                        sp.energyCutoff_eV = cutoff_eV;
                }
                else if (k == "init")
                    sp.init = unquote(tokens[1]);
                else if (k == "temperature" || k == "temp" || k == "te")
                {
                    double te_eV = 0.0;
                    if (parseTeToEV(tokens, 1, te_eV))
                        sp.Te_eV = te_eV;
                }
                else if (k == "sourceenabled" || k == "source")
                    parseBoolToken(tokens[1], sp.source.enabled);
                else if (k == "sourceboundary" || k == "sourceboundaryid" || k == "sourceboundaryids")
                    sp.source.boundaryId = stoi(tokens[1]);
                else if (k == "sourcen" || k == "sourcenumberdensity" || k == "sourcedensity")
                {
                    double v = 0.0;
                    if (parseNumberDensityToSI(tokens, 1, v))
                        sp.source.numberDensity_m3 = v;
                }
                else if (k == "sourcetemperature" || k == "sourcete")
                {
                    double te_eV = 0.0;
                    if (parseTeToEV(tokens, 1, te_eV))
                        sp.source.temperature_eV = te_eV;
                }
                else if (k == "sourcecurrent" || k == "sourcei")
                {
                    double v = 0.0;
                    if (parseScalarToSI(tokens, 1, v))
                        sp.source.current_A = v;
                }
                else if (k == "sourceevery")
                    sp.source.every = stoi(tokens[1]);
                else if (k == "sourcemaxmacro" || k == "sourcemaxemit")
                    sp.source.maxEmitMacro = stoi(tokens[1]);
            }
            continue;
        }

        // Ignore unknown lines.
    }

    file.close();

    if (!sawPicBlock)
    {
        std::cerr << "[DictionaryPIC] ERROR: missing PIC { } block in " << picDictPath << std::endl;
        return false;
    }
    if (!sawSpeciesBlock || !sawSpeciesEntry || config_.species.empty())
    {
        std::cerr << "[DictionaryPIC] ERROR: missing Species { } block or species entries in " << picDictPath << std::endl;
        return false;
    }
    if (!sawDt)
    {
        std::cerr << "[DictionaryPIC] WARNING: missing dt in PIC block; using default dt=" << config_.dt << std::endl;
    }
    if (!sawSteps)
    {
        std::cerr << "[DictionaryPIC] WARNING: missing steps in PIC block; using default steps=" << config_.steps << std::endl;
    }

    // Defaults / sanitation
    if (config_.dt <= 0.0)
        config_.dt = 1e-8;
    if (config_.steps <= 0)
        config_.steps = 10;
    if (config_.output.every <= 0)
        config_.output.every = 10;
    if (config_.output.console.every <= 0)
        config_.output.console.every = 10;
    if (config_.output.ebBandRMS.bandLayers < 0)
        config_.output.ebBandRMS.bandLayers = 0;
    if (config_.output.diagnostics.every < 0)
        config_.output.diagnostics.every = 0;
    if (config_.output.diagnostics.baseDir.empty())
        config_.output.diagnostics.baseDir = "diagnostics";
    if (config_.output.diagnostics.energyHistogram.bins <= 0)
        config_.output.diagnostics.energyHistogram.bins = 80;
    if (!(config_.output.diagnostics.energyHistogram.maxEnergyEV > 0.0))
        config_.output.diagnostics.energyHistogram.maxEnergyEV = 200.0;
    if (config_.output.particleTrajectory.every < 0)
        config_.output.particleTrajectory.every = 0;
    if (config_.output.particleTrajectory.electronCount < 0)
        config_.output.particleTrajectory.electronCount = 0;
    if (config_.output.particleTrajectory.ionCount < 0)
        config_.output.particleTrajectory.ionCount = 0;
    if (config_.output.particleTrajectory.electronMaxSamples < 0)
        config_.output.particleTrajectory.electronMaxSamples = 0;
    if (config_.output.particleTrajectory.ionMaxSamples < 0)
        config_.output.particleTrajectory.ionMaxSamples = 0;
    if (config_.output.particleTrajectory.baseDir.empty())
        config_.output.particleTrajectory.baseDir = "diagnostics/particle_trajectories";
    config_.output.v2.fieldsEvery = std::max(1, config_.output.every);
    if (config_.output.v2.seriesEvery <= 0)
        config_.output.v2.seriesEvery = 100;
    config_.output.v2.perfEvery = config_.output.v2.seriesEvery;
    config_.output.v2.physicsEvery = config_.output.v2.seriesEvery;
    config_.output.v2.speciesEvery = config_.output.v2.seriesEvery;
    config_.output.v2.boundaryEvery = config_.output.v2.seriesEvery;
    config_.output.v2.flushEvery = config_.output.v2.seriesEvery;
    config_.output.v2.amountEvery = config_.output.amount.enabled ? config_.output.v2.seriesEvery : 0;
    {
        std::string fmt = toLowerCopy(config_.output.v2.console.format);
        if (fmt != "oneline" && fmt != "multiline")
            fmt = "oneline";
        config_.output.v2.console.format = fmt;
    }
    config_.output.diagnostics.enabled = config_.output.diagnostics.enabled ||
        config_.output.diagnostics.electrodeCurrent.enabled ||
        config_.output.diagnostics.reactions.enabled ||
        config_.output.diagnostics.energyHistogram.enabled;

    for (auto& sp : config_.species)
    {
        if (sp.m <= 0.0)
            sp.m = 1.0;
        if (sp.weight <= 0.0)
            sp.weight = 1.0;
        if (sp.numParticles < 0)
            sp.numParticles = 0;
        if (sp.energyCutoff_eV < 0.0)
            sp.energyCutoff_eV = 0.0;
    }

    // SEE sanitation
    if (config_.see.gamma < 0.0)
        config_.see.gamma = 0.0;
    if (config_.see.gamma0 < 0.0)
        config_.see.gamma0 = 0.0;
    if (config_.see.Eref_eV <= 0.0)
        config_.see.Eref_eV = 100.0;
    if (config_.see.gammaMax < 0.0)
        config_.see.gammaMax = 0.0;
    if (config_.see.Emin_eV < 0.0)
        config_.see.Emin_eV = 0.0;
    if (config_.see.emittedEnergy_eV < 0.0)
        config_.see.emittedEnergy_eV = 0.0;
    if (config_.see.emittedTe_eV < 0.0)
        config_.see.emittedTe_eV = 0.0;
    if (config_.see.emittedEmean_eV < 0.0)
        config_.see.emittedEmean_eV = 0.0;
    if (config_.see.injectEps_m < 0.0)
        config_.see.injectEps_m = 0.0;
    if (config_.see.maxEmitMacro < 0)
        config_.see.maxEmitMacro = 0;

    // Normalize SEE gamma model string + table.
    {
        std::string m = toLowerCopy(config_.see.gammaModel);
        if (m.empty())
            m = "constant";
        if (m != "constant" && m != "powerlaw" && m != "piecewise")
            m = "constant";
        config_.see.gammaModel = m;

        if (m == "piecewise")
        {
            const size_t n = std::min(config_.see.piecewiseE_eV.size(), config_.see.piecewiseGamma.size());
            config_.see.piecewiseE_eV.resize(n);
            config_.see.piecewiseGamma.resize(n);

            std::vector<std::pair<double, double>> pts;
            pts.reserve(n);
            for (size_t i = 0; i < n; ++i)
            {
                const double E = config_.see.piecewiseE_eV[i];
                const double g = config_.see.piecewiseGamma[i];
                if (E < 0.0)
                    continue;
                pts.emplace_back(E, std::max(0.0, g));
            }
            std::sort(pts.begin(), pts.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
            // De-duplicate energies by keeping the last value.
            std::vector<std::pair<double, double>> uniq;
            for (const auto& p : pts)
            {
                if (!uniq.empty() && std::abs(p.first - uniq.back().first) < 1e-12)
                    uniq.back().second = p.second;
                else
                    uniq.push_back(p);
            }

            config_.see.piecewiseE_eV.clear();
            config_.see.piecewiseGamma.clear();
            for (const auto& p : uniq)
            {
                config_.see.piecewiseE_eV.push_back(p.first);
                config_.see.piecewiseGamma.push_back(p.second);
            }

            // If the table is not usable, fall back to constant.
            if (config_.see.piecewiseE_eV.size() < 2)
                config_.see.gammaModel = "constant";
        }
    }

    // Normalize SEE emitted-energy model string, and tie defaults to legacy emittedEnergy_eV.
    {
        std::string em = toLowerCopy(config_.see.emittedEnergyModel);
        if (em.empty())
            em = "mono";
        if (em == "mono" || em == "monochromatic" || em == "constant")
            em = "mono";
        else if (em == "maxwell" || em == "maxwellian")
            em = "maxwell";
        else if (em == "exponential" || em == "exp")
            em = "exponential";
        else
            em = "mono";
        config_.see.emittedEnergyModel = em;

        // If the new parameters are not specified, fall back to legacy emittedEnergy_eV.
        if (config_.see.emittedTe_eV <= 0.0)
            config_.see.emittedTe_eV = std::max(0.0, config_.see.emittedEnergy_eV);
        if (config_.see.emittedEmean_eV <= 0.0)
            config_.see.emittedEmean_eV = std::max(0.0, config_.see.emittedEnergy_eV);
    }

    // Normalize SEE angular model string.
    {
        const std::string a = toLowerCopy(config_.see.angular);
        if (a == "cosine" || a == "isotropic" || a == "fixed")
            config_.see.angular = a;
        else
            config_.see.angular = "cosine";
    }

    return true;
}
