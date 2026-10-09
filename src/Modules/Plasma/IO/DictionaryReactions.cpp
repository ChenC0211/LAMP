#include "Modules/Plasma/IO/DictionaryReactions.h"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

#include "Modules/Plasma/Physics/CrossSectionDB.h"
#include "Modules/Plasma/Physics/PhysicalConstants.h"
#include "Core/Units/UnitParser.h"

using namespace std;

namespace
{
string normalizeBlockName(const string& s)
{
    string out = s;
    if (out.size() >= 2 && ((out.front() == '"' && out.back() == '"') || (out.front() == '\'' && out.back() == '\'')))
        out = out.substr(1, out.size() - 2);
    transform(out.begin(), out.end(), out.begin(), [](unsigned char c) { return static_cast<char>(tolower(c)); });
    return out;
}

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

bool normalizeEnergyFrameToken(const std::string& text, std::string& out)
{
    std::string v = text;
    if (v.size() >= 2 && ((v.front() == '"' && v.back() == '"') || (v.front() == '\'' && v.back() == '\'')))
        v = v.substr(1, v.size() - 2);
    transform(v.begin(), v.end(), v.begin(), [](unsigned char c) { return static_cast<char>(tolower(c)); });

    if (v == "lab")
    {
        out = "lab";
        return true;
    }
    if (v == "centerofmass")
    {
        out = "centerOfMass";
        return true;
    }
    return false;
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

// Parse an energy value to eV. Accepts:
//   - "15.76 eV"
//   - "2.5 J"
//   - "15.76" (interpreted as eV)
bool parseEnergyToEV(const std::string& text, double& out_eV)
{
    UnitParser::Result r;
    if (UnitParser::TryParseToSI(text, r))
    {
        if (!r.has_unit)
        {
            out_eV = r.value_si;
            return true;
        }
        // value_si is in Joules for eV/J units.
        if (r.unit == "ev" || r.unit == "j")
        {
            out_eV = r.value_si / PhysicalConstants::eV_J;
            return true;
        }
        return false;
    }

    try
    {
        out_eV = stod(text);
        return true;
    }
    catch (...)
    {
        return false;
    }
}

bool parseStringList(const std::string& line, vector<string>& out)
{
    // Expect something like:
    //   products ( e Ar+ e )
    // or: products ("e" "Ar+" "e")
    const size_t l = line.find('(');
    const size_t r = line.find(')');
    if (l == string::npos || r == string::npos || r <= l)
        return false;

    string inside = line.substr(l + 1, r - l - 1);
    // Normalize separators
    for (char& c : inside)
    {
        if (c == ',' || c == ';')
            c = ' ';
    }

    std::istringstream iss(inside);
    string tok;
    out.clear();
    while (iss >> tok)
    {
        if (tok.size() >= 2 && ((tok.front() == '"' && tok.back() == '"') || (tok.front() == '\'' && tok.back() == '\'')))
            tok = tok.substr(1, tok.size() - 2);
        out.push_back(tok);
    }
    return true;
}
}

string DictionaryReactions::toLowerCopy(string s)
{
    transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return static_cast<char>(tolower(c)); });
    return s;
}

string DictionaryReactions::stripQuotesAndSemicolon(const std::string& s)
{
    std::string out = s;
    if (!out.empty() && out.back() == ';')
        out.pop_back();
    if (out.size() >= 2 && ((out.front() == '"' && out.back() == '"') || (out.front() == '\'' && out.back() == '\'')))
        out = out.substr(1, out.size() - 2);
    return out;
}

bool DictionaryReactions::readDict(const std::string& path, const CrossSectionDB* xs_db)
{
    ifstream file(path);
    if (!file.is_open())
        return false;

    bool foundReactions = false;
    vector<string> blocks;
    string pendingBlock;

    bool inSet = false;
    bool inChannels = false;
    bool inChannel = false;

    ReactionSet curSet;
    ReactionChannel curCh;

    auto currentBlock = [&]() -> string { return blocks.empty() ? string() : blocks.back(); };
    auto parentBlock = [&]() -> string {
        if (blocks.size() < 2) return string();
        return blocks[blocks.size() - 2];
    };
    auto grandParentBlock = [&]() -> string {
        if (blocks.size() < 3) return string();
        return blocks[blocks.size() - 3];
    };

    auto openBlock = [&](const string& name) {
        blocks.push_back(name);
        const string low = normalizeBlockName(name);
        if (low == "reactions")
            foundReactions = true;

        // A set block is directly under Reactions.
        if (normalizeBlockName(parentBlock()) == "reactions" && low != "reactions")
        {
            inSet = true;
            curSet = ReactionSet{};
            curSet.name = stripQuotesAndSemicolon(name);
            return;
        }

        if (inSet && normalizeBlockName(parentBlock()) == normalizeBlockName(curSet.name) && (low == "channels" || low == "channel"))
        {
            inChannels = true;
            return;
        }

        const string parentLow = normalizeBlockName(parentBlock());
        if (inSet && inChannels && (parentLow == "channels" || parentLow == "channel"))
        {
            inChannel = true;
            curCh = ReactionChannel{};
            curCh.name = stripQuotesAndSemicolon(name);
            return;
        }
    };

    auto finalizeChannel = [&]() -> bool {
        if (!inChannel)
            return true;
        // Basic validation.
        if (curCh.name.empty())
        {
            cerr << "[DictionaryReactions] ERROR: channel name is empty" << endl;
            return false;
        }
        // Allow disabled channels to omit crossSection entirely (placeholder/feature flag).
        if (curCh.enabled && curCh.cross_section_key.empty())
        {
            cerr << "[DictionaryReactions] ERROR: channel '" << curCh.name << "' missing crossSection" << endl;
            return false;
        }
        curSet.channels.push_back(curCh);
        inChannel = false;
        return true;
    };

    auto finalizeSet = [&]() -> bool {
        if (!inSet)
            return true;
        string err;
        if (!curSet.Validate(xs_db, &err))
        {
            cerr << "[DictionaryReactions] ERROR: invalid ReactionSet '" << curSet.name << "': " << err << endl;
            return false;
        }
        if (!db_.Register(curSet, false, &err))
        {
            cerr << "[DictionaryReactions] ERROR: failed to register ReactionSet '" << curSet.name << "': " << err << endl;
            return false;
        }
        inSet = false;
        inChannels = false;
        inChannel = false;
        return true;
    };

    auto closeBlock = [&]() -> bool {
        // Closing order matters: channel -> channels -> set.
        const string parentLow = normalizeBlockName(parentBlock());
        if (inChannel && (parentLow == "channels" || parentLow == "channel"))
        {
            if (!finalizeChannel())
                return false;
        }
        else if (inChannels && (normalizeBlockName(currentBlock()) == "channels" || normalizeBlockName(currentBlock()) == "channel"))
        {
            inChannels = false;
        }
        else if (inSet && normalizeBlockName(parentBlock()) == "reactions")
        {
            // closing set block
            if (!finalizeSet())
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
            const string before = trim(line.substr(0, bracePos));
            if (!before.empty())
                openBlock(before);
            else if (!pendingBlock.empty())
            {
                openBlock(pendingBlock);
                pendingBlock.clear();
            }
            continue;
        }

        // Block titles waiting for '{'
        const string lowLine = toLowerCopy(line);
        if (lowLine == "reactions")
        {
            pendingBlock = line;
            continue;
        }
        if (inSet && !inChannels && (lowLine == "channels" || lowLine == "channel"))
        {
            pendingBlock = line;
            continue;
        }

        // Key-value lines
        const vector<string> tokens = split(line);
        
        // Single token could be a block name (e.g., "electron-Ar", "ELASTIC") waiting for next "{"
        if (tokens.size() == 1)
        {
            string stripped = stripQuotesAndSemicolon(tokens[0]);
            if (!stripped.empty())
            {
                pendingBlock = stripped;
            }
            continue;
        }
        
        if (tokens.size() < 2)
            continue;

        const string k = toLowerCopy(tokens[0]);
        if (inChannel)
        {
            if (k == "crosssection" || k == "xs")
                curCh.cross_section_key = stripQuotesAndSemicolon(joinTokens(tokens, 1));
            else if (k == "threshold" || k == "eth" || k == "ethreshold")
            {
                double ev = 0.0;
                if (!parseEnergyToEV(joinTokens(tokens, 1), ev))
                {
                    cerr << "[DictionaryReactions] ERROR: failed to parse threshold in line: " << line << endl;
                    return false;
                }
                curCh.threshold_eV = ev;
            }
            else if (k == "products" || k == "product")
            {
                vector<string> list;
                if (!parseStringList(line, list))
                {
                    cerr << "[DictionaryReactions] ERROR: failed to parse products list in line: " << line << endl;
                    return false;
                }
                curCh.products = list;
            }
            else if (k == "scatter" || k == "scattermodel")
                curCh.scatter_model = stripQuotesAndSemicolon(joinTokens(tokens, 1));
            else if (k == "energyloss" || k == "loss")
                curCh.energy_loss_model = stripQuotesAndSemicolon(joinTokens(tokens, 1));
            else if (k == "energyshare" || k == "share")
                curCh.energy_share_model = stripQuotesAndSemicolon(joinTokens(tokens, 1));
            else if (k == "enabled")
            {
                bool b = true;
                if (!parseBool(stripQuotesAndSemicolon(tokens[1]), b))
                {
                    cerr << "[DictionaryReactions] ERROR: invalid bool for enabled in line: " << line << endl;
                    return false;
                }
                curCh.enabled = b;
            }
            continue;
        }

        if (inSet && !inChannels)
        {
            if (k == "projectile" || k == "proj")
                curSet.projectile = stripQuotesAndSemicolon(tokens[1]);
            else if (k == "target" || k == "gas" || k == "neutral")
                curSet.target = stripQuotesAndSemicolon(tokens[1]);
            else if (k == "energyframe" || k == "energy_frame")
            {
                std::string normalized;
                if (!normalizeEnergyFrameToken(tokens[1], normalized))
                {
                    cerr << "[DictionaryReactions] ERROR: invalid energyFrame in line: " << line << endl;
                    return false;
                }
                curSet.energy_frame = normalized;
            }
            continue;
        }
    }

    while (!blocks.empty())
    {
        if (!closeBlock())
            return false;
    }

    if (!foundReactions)
    {
        cerr << "[DictionaryReactions] ERROR: missing Reactions{ } block" << endl;
        return false;
    }
    return true;
}
