#include "Modules/Plasma/IO/DictionaryCollisions.h"

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

#include "Modules/Plasma/Physics/PhysicalConstants.h"
#include "Core/Units/UnitParser.h"

using namespace std;

namespace
{
string joinTokens(const vector<string>& tokens, size_t begin)
{
    string out;
    for (size_t i = begin; i < tokens.size(); ++i)
    {
        if (!out.empty())
            out += " ";
        out += tokens[i];
    }
    return out;
}

bool parseBool(const std::string& s, bool& out)
{
    std::string v = s;
    transform(v.begin(), v.end(), v.begin(), [](unsigned char c) { return static_cast<char>(tolower(c)); });
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

bool parseSizeT(const std::string& s, std::size_t& out)
{
    try
    {
        size_t idx = 0;
        const unsigned long long v = stoull(s, &idx);
        if (idx != s.size())
            return false;
        out = static_cast<size_t>(v);
        return true;
    }
    catch (...)
    {
        return false;
    }
}

bool parseUint64(const std::string& s, std::uint64_t& out)
{
    try
    {
        size_t idx = 0;
        const unsigned long long v = stoull(s, &idx);
        if (idx != s.size())
            return false;
        out = static_cast<std::uint64_t>(v);
        return true;
    }
    catch (...)
    {
        return false;
    }
}

bool parseInt(const std::string& s, int& out)
{
    try
    {
        size_t idx = 0;
        const int v = stoi(s, &idx);
        if (idx != s.size())
            return false;
        out = v;
        return true;
    }
    catch (...)
    {
        return false;
    }
}

}

string DictionaryCollisions::toLowerCopy(string s)
{
    transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return static_cast<char>(tolower(c)); });
    return s;
}

string DictionaryCollisions::stripQuotesAndSemicolon(const std::string& s)
{
    std::string out = s;
    if (!out.empty() && out.back() == ';')
        out.pop_back();
    if (out.size() >= 2 && ((out.front() == '"' && out.back() == '"') || (out.front() == '\'' && out.back() == '\'')))
        out = out.substr(1, out.size() - 2);
    return out;
}

bool DictionaryCollisions::readDict(const std::string& path)
{
    ifstream file(path);
    if (!file.is_open())
        return false;

    cfg_ = CollisionsConfig{};

    struct Entry
    {
        std::string label;
        std::string projectile;
        std::string reactionSet;
    };

    bool foundCollisions = false;
    vector<string> blocks;
    string pendingBlock;
    bool inEntry = false;
    Entry cur;
    vector<CollisionMapping> mappings;

    auto currentBlock = [&]() -> std::string { return blocks.empty() ? std::string() : blocks.back(); };
    auto parentBlock = [&]() -> std::string {
        if (blocks.size() < 2) return std::string();
        return blocks[blocks.size() - 2];
    };

    auto openBlock = [&](const std::string& name) {
        blocks.push_back(name);
        if (toLowerCopy(name) == "collisions")
            foundCollisions = true;
        if (toLowerCopy(parentBlock()) == "collisions" && toLowerCopy(name) != "collisions")
        {
            inEntry = true;
            cur = Entry{};
            cur.label = stripQuotesAndSemicolon(name);
        }
    };

    auto finalizeEntry = [&]() -> bool {
        if (!inEntry)
            return true;
        if (cur.projectile.empty())
        {
            cerr << "[DictionaryCollisions] ERROR: entry '" << cur.label << "' missing projectile" << endl;
            return false;
        }
        if (cur.reactionSet.empty())
        {
            cerr << "[DictionaryCollisions] ERROR: entry '" << cur.label << "' missing reactionSet" << endl;
            return false;
        }
        CollisionMapping m;
        m.projectile = cur.projectile;
        m.reaction_set = cur.reactionSet;
        mappings.push_back(std::move(m));
        inEntry = false;
        return true;
    };

    auto closeBlock = [&]() -> bool {
        if (inEntry && toLowerCopy(parentBlock()) == "collisions" && !blocks.empty())
        {
            if (!finalizeEntry())
                return false;
        }
        if (!blocks.empty())
            blocks.pop_back();
        return true;
    };

    string line;
    while (getline(file, line))
    {
        line = removeInlineComments(line);
        line = trim(line);
        if (line.empty())
            continue;

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
            if (!closeBlock())
                return false;
            continue;
        }

        const size_t bracePos = line.find('{');
        if (bracePos != string::npos)
        {
            string before = trim(line.substr(0, bracePos));
            if (!before.empty())
                openBlock(before);
            else if (!pendingBlock.empty())
            {
                openBlock(pendingBlock);
                pendingBlock.clear();
            }
            continue;
        }

        if (toLowerCopy(line) == "collisions")
        {
            pendingBlock = line;
            continue;
        }

        const vector<string> tokens = split(line);
        
        // Single token could be a block name (e.g., "electron") waiting for next "{"
        if (tokens.size() == 1)
        {
            string stripped = stripQuotesAndSemicolon(tokens[0]);
            if (!stripped.empty() && toLowerCopy(stripped) != "collisions")
            {
                pendingBlock = stripped;
            }
            continue;
        }
        
        if (tokens.size() < 2)
            continue;

        const string k = toLowerCopy(tokens[0]);
        const string v0 = stripQuotesAndSemicolon(tokens[1]);

        // Top-level keys under Collisions
        if (!inEntry && toLowerCopy(currentBlock()) == "collisions")
        {
            if (k == "enabled")
            {
                bool b = false;
                if (!parseBool(v0, b))
                {
                    cerr << "[DictionaryCollisions] ERROR: invalid enabled: " << v0 << endl;
                    return false;
                }
                cfg_.enabled = b;
            }
            else if (k == "parallel")
            {
                bool b = true;
                if (!parseBool(v0, b))
                {
                    cerr << "[DictionaryCollisions] ERROR: invalid parallel: " << v0 << endl;
                    return false;
                }
                cfg_.parallel = b;
            }
            else if (k == "rngmode" || k == "rng_mode")
            {
                CollisionsRngMode m;
                if (!ParseCollisionsRngMode(v0, &m))
                {
                    cerr << "[DictionaryCollisions] ERROR: invalid rngMode: " << v0 << endl;
                    return false;
                }
                cfg_.rng_mode = m;
            }
            else if (k == "rngseed" || k == "rng_seed")
            {
                std::uint64_t s = 0;
                if (!parseUint64(v0, s))
                {
                    cerr << "[DictionaryCollisions] ERROR: invalid rngSeed: " << v0 << endl;
                    return false;
                }
                cfg_.rng_seed = s;
            }

            else if (k == "stats.every" || k == "statsevery" || k == "stats_every")
            {
                int n = 0;
                if (!parseInt(v0, n))
                {
                    cerr << "[DictionaryCollisions] ERROR: invalid stats.every: " << v0 << endl;
                    return false;
                }
                cfg_.stats_every = n;
            }

            else if (k == "trackfastneutrals" || k == "track_fast_neutrals" || k == "trackfastneutral")
            {
                bool b = false;
                if (!parseBool(v0, b))
                {
                    cerr << "[DictionaryCollisions] ERROR: invalid trackFastNeutrals: " << v0 << endl;
                    return false;
                }
                cfg_.track_fast_neutrals = b;
            }
            else if (k == "usegasthermalmotion" || k == "use_gas_thermal_motion" || k == "gasthermalmotion")
            {
                bool b = false;
                if (!parseBool(v0, b))
                {
                    cerr << "[DictionaryCollisions] ERROR: invalid useGasThermalMotion: " << v0 << endl;
                    return false;
                }
                cfg_.use_gas_thermal_motion = b;
            }
            else if (k == "method")
            {
                string m = toLowerCopy(v0);
                if (m != "nullcollision" && m != "null" && m != "majorant")
                {
                    cerr << "[DictionaryCollisions] ERROR: unsupported method: " << v0 << endl;
                    return false;
                }
            }
            else if (k == "majorantemax" || k == "emax")
            {
                // value + unit; default unit is eV when unit is omitted.
                const string expr = stripQuotesAndSemicolon(joinTokens(tokens, 1));
                UnitParser::Result r;
                if (UnitParser::TryParseToSI(expr, r) && r.has_unit)
                {
                    // UnitParser converts eV -> Joules.
                    if (r.unit == "ev")
                    {
                        cfg_.null_collision.majorant_Emax_eV = r.value_si / PhysicalConstants::eV_J;
                    }
                    else
                    {
                        cerr << "[DictionaryCollisions] ERROR: majorantEmax expects eV (or unitless), got unit='" << r.unit
                             << "'" << endl;
                        return false;
                    }
                }
                else
                {
                    // Unitless: interpret as eV.
                    try
                    {
                        cfg_.null_collision.majorant_Emax_eV = stod(expr);
                    }
                    catch (...)
                    {
                        cerr << "[DictionaryCollisions] ERROR: majorantEmax parse failed: " << expr << endl;
                        return false;
                    }
                }
            }
            else if (k == "safetyfactor" || k == "safety")
            {
                try
                {
                    cfg_.null_collision.safety_factor = stod(v0);
                }
                catch (...)
                {
                    cerr << "[DictionaryCollisions] ERROR: invalid safetyFactor: " << v0 << endl;
                    return false;
                }
            }
            else if (k == "energysamples" || k == "samples")
            {
                size_t n = 0;
                if (!parseSizeT(v0, n))
                {
                    cerr << "[DictionaryCollisions] ERROR: invalid energySamples: " << v0 << endl;
                    return false;
                }
                cfg_.null_collision.energy_samples = n;
            }
            else if (k == "clampenergy" || k == "clamp")
            {
                bool b = true;
                if (!parseBool(v0, b))
                {
                    cerr << "[DictionaryCollisions] ERROR: invalid clampEnergy: " << v0 << endl;
                    return false;
                }
                cfg_.null_collision.clamp_energy = b;
            }
            else if (k == "usefastlookup" || k == "fastlookup" || k == "use_lookup" || k == "uselookup")
            {
                bool b = true;
                if (!parseBool(v0, b))
                {
                    cerr << "[DictionaryCollisions] ERROR: invalid useFastLookup: " << v0 << endl;
                    return false;
                }
                cfg_.null_collision.use_fast_lookup = b;
            }
            else if (k == "lookupbins" || k == "lookup_bins" || k == "bins")
            {
                size_t n = 0;
                if (!parseSizeT(v0, n))
                {
                    cerr << "[DictionaryCollisions] ERROR: invalid lookupBins: " << v0 << endl;
                    return false;
                }
                cfg_.null_collision.lookup_bins = n;
            }
            continue;
        }

        // Nested stats block: Collisions { stats { every N; } }
        if (!inEntry && toLowerCopy(currentBlock()) == "stats" && toLowerCopy(parentBlock()) == "collisions")
        {
            if (k == "every")
            {
                int n = 0;
                if (!parseInt(v0, n))
                {
                    cerr << "[DictionaryCollisions] ERROR: invalid stats.every: " << v0 << endl;
                    return false;
                }
                cfg_.stats_every = n;
            }
            continue;
        }

        // Entry keys
        if (inEntry && toLowerCopy(parentBlock()) == "collisions")
        {
            if (k == "projectile")
                cur.projectile = stripQuotesAndSemicolon(v0);
            else if (k == "reactionset" || k == "reaction" || k == "set")
                cur.reactionSet = stripQuotesAndSemicolon(joinTokens(tokens, 1));
        }
    }

    if (!foundCollisions)
    {
        cerr << "[DictionaryCollisions] ERROR: missing 'Collisions' block" << endl;
        return false;
    }

    cfg_.mappings = std::move(mappings);
    string err;
    if (!cfg_.Validate(&err))
    {
        cerr << "[DictionaryCollisions] ERROR: " << err << endl;
        return false;
    }

    return true;
}
