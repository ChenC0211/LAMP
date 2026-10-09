#include "Core/IO/DictionaryUnits.h"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>

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
} // namespace

string DictionaryUnits::toLowerCopy(string s)
{
    transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return static_cast<char>(tolower(c)); });
    return s;
}

bool DictionaryUnits::parseBoolToken(const std::string& s, bool& out)
{
    const std::string v = toLowerCopy(s);
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

static bool parseValuePlusUnitToSI(const std::vector<std::string>& tokens,
                                  size_t value_index,
                                  double& out_si)
{
    if (tokens.size() <= value_index)
        return false;

    const string joined = joinTokens(tokens, value_index);

    UnitParser::Result r;
    if (UnitParser::TryParseToSI(joined, r))
    {
        out_si = r.value_si;
        return true;
    }

    // Fallback: assume the value is already SI.
    try
    {
        out_si = stod(tokens[value_index]);
        return true;
    }
    catch (...)
    {
        return false;
    }
}

bool DictionaryUnits::readDict(const string& unitsDictPath)
{
    ifstream file(unitsDictPath);
    if (!file.is_open())
        return false;

    UnitsScales scales; // defaults to 1s
    UnitsMode mode = UnitsMode::kCodeUnits;
    writeSIOutputs_ = false;

    bool inUnitsBlock = false;
    string line;

    while (getline(file, line))
    {
        line = removeInlineComments(line);
        line = trim(line);

        if (line.empty() || line == "{")
            continue;

        if (line == "Units")
        {
            inUnitsBlock = true;
            continue;
        }

        if (line == "}")
        {
            if (inUnitsBlock)
                inUnitsBlock = false;
            continue;
        }

        // Accept both inside Units{} and at top-level.
        vector<string> tokens = split(line);
        if (tokens.empty())
            continue;

        const string key = toLowerCopy(tokens[0]);

        if (key == "mode" && tokens.size() >= 2)
        {
            const string m = toLowerCopy(tokens[1]);
            if (m == "si")
                mode = UnitsMode::kSI;
            else if (m == "codeunits" || m == "code" || m == "legacy")
                mode = UnitsMode::kCodeUnits;
            else
                cerr << "[DictionaryUnits] WARNING: unknown mode '" << tokens[1] << "' (use SI|CodeUnits)" << endl;
            continue;
        }

        if ((key == "writesi" || key == "writesios" || key == "writesioutputs") && tokens.size() >= 2)
        {
            bool b = false;
            if (parseBoolToken(tokens[1], b))
                writeSIOutputs_ = b;
            else
                cerr << "[DictionaryUnits] WARNING: failed to parse writeSI from line: " << line << endl;
            continue;
        }

        // If a Units{} block exists, ignore keys outside the block except 'mode'.
        if (!inUnitsBlock)
        {
            // Allow single-line syntax without blocks by treating the whole file as implicit Units block.
            // If the file uses block form, inUnitsBlock will be true and this path won't be used.
        }

        double val_si = 0.0;
        if (key == "l0")
        {
            if (parseValuePlusUnitToSI(tokens, 1, val_si))
                scales.L0_m = val_si;
            else
                cerr << "[DictionaryUnits] WARNING: failed to parse L0 from line: " << line << endl;
        }
        else if (key == "t0")
        {
            if (parseValuePlusUnitToSI(tokens, 1, val_si))
                scales.T0_s = val_si;
            else
                cerr << "[DictionaryUnits] WARNING: failed to parse T0 from line: " << line << endl;
        }
        else if (key == "phi0")
        {
            if (parseValuePlusUnitToSI(tokens, 1, val_si))
                scales.Phi0_V = val_si;
            else
                cerr << "[DictionaryUnits] WARNING: failed to parse Phi0 from line: " << line << endl;
        }
        else if (key == "n0")
        {
            if (parseValuePlusUnitToSI(tokens, 1, val_si))
                scales.n0_m3 = val_si;
            else
                cerr << "[DictionaryUnits] WARNING: failed to parse n0 from line: " << line << endl;
        }
        else if (key == "v0")
        {
            if (parseValuePlusUnitToSI(tokens, 1, val_si))
                scales.v0_mps = val_si;
            else
                cerr << "[DictionaryUnits] WARNING: failed to parse v0 from line: " << line << endl;
        }
        else if (key == "e0")
        {
            if (parseValuePlusUnitToSI(tokens, 1, val_si))
                scales.E0_Vpm = val_si;
            else
                cerr << "[DictionaryUnits] WARNING: failed to parse E0 from line: " << line << endl;
        }
        else if (key == "rho0")
        {
            if (parseValuePlusUnitToSI(tokens, 1, val_si))
                scales.rho0_Cpm3 = val_si;
            else
                cerr << "[DictionaryUnits] WARNING: failed to parse rho0 from line: " << line << endl;
        }
        else
        {
            // Silently ignore unknown keys for forward compatibility.
        }
    }

    file.close();

    try
    {
        if (mode == UnitsMode::kSI)
            context_ = UnitsContext::MakeSI(scales);
        else
            context_ = UnitsContext::MakeCodeUnits();
    }
    catch (const std::exception& e)
    {
        cerr << "[DictionaryUnits] ERROR: failed to build UnitsContext: " << e.what() << endl;
        context_ = UnitsContext::MakeCodeUnits();
    }

    return true;
}
