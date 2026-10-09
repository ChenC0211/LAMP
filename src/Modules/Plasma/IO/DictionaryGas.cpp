#include "Modules/Plasma/IO/DictionaryGas.h"

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

bool parseValuePlusUnitToSI(const std::vector<std::string>& tokens,
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

    // Fallback: assume already SI.
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

} // namespace

std::string DictionaryGas::toLowerCopy(std::string s)
{
    transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return static_cast<char>(tolower(c)); });
    return s;
}

bool DictionaryGas::readDict(const std::string& gasDictPath)
{
    ifstream file(gasDictPath);
    if (!file.is_open())
        return false;

    string name = "";
    double pressure_Pa = 0.0;
    double temperature_K = 0.0;

    bool inGasBlock = false;
    string line;

    while (getline(file, line))
    {
        line = removeInlineComments(line);
        line = trim(line);

        if (line.empty() || line == "{")
            continue;

        if (line == "Gas")
        {
            inGasBlock = true;
            continue;
        }

        if (line == "}")
        {
            if (inGasBlock)
                inGasBlock = false;
            continue;
        }

        // Accept both inside Gas{} and at top-level.
        vector<string> tokens = split(line);
        if (tokens.empty())
            continue;

        const string key = toLowerCopy(tokens[0]);

        if (key == "name" || key == "gas" || key == "species")
        {
            if (tokens.size() >= 2)
                name = tokens[1];
            continue;
        }

        if (key == "p" || key == "pressure")
        {
            double v = 0.0;
            if (parseValuePlusUnitToSI(tokens, 1, v))
                pressure_Pa = v;
            else
                cerr << "[DictionaryGas] WARNING: failed to parse pressure from line: " << line << endl;
            continue;
        }

        if (key == "t" || key == "temperature")
        {
            double v = 0.0;
            if (parseValuePlusUnitToSI(tokens, 1, v))
                temperature_K = v;
            else
                cerr << "[DictionaryGas] WARNING: failed to parse temperature from line: " << line << endl;
            continue;
        }

        // If a Gas{} block exists, ignore unknown keys outside the block except known ones.
        if (!inGasBlock)
        {
            // allow implicit block
        }
    }

    file.close();

    if (name.empty())
        name = "Ar";
    if (temperature_K <= 0.0)
        temperature_K = 300.0;
    if (pressure_Pa < 0.0)
        pressure_Pa = 0.0;

    BackgroundGas gas(name, pressure_Pa, temperature_K);
    double m = 0.0;
    if (BackgroundGas::LookupParticleMassKg(name, &m))
        gas.SetParticleMassKg(m);

    string err;
    if (!gas.Validate(&err))
    {
        cerr << "[DictionaryGas] WARNING: invalid gas specification: " << err << " (using vacuum)" << endl;
        BackgroundGas vacuum(name, 0.0, temperature_K);
        if (BackgroundGas::LookupParticleMassKg(name, &m))
            vacuum.SetParticleMassKg(m);
        state_ = vacuum.Evaluate();
        return true;
    }

    state_ = gas.Evaluate();
    return true;
}
