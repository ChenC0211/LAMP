#include "Modules/Plasma/IO/DictionaryCrossSections.h"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

#include "Modules/Plasma/Physics/BolsigCrossSectionFile.h"
#include "Modules/Plasma/Physics/BuiltInCrossSections.h"
#include "Modules/Plasma/Physics/TabulatedCrossSection.h"

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

TabulatedCrossSection::Interpolation parseInterp(const std::string& s, bool* ok)
{
    std::string v = s;
    transform(v.begin(), v.end(), v.begin(), [](unsigned char c) { return static_cast<char>(tolower(c)); });
    if (v == "lin" || v == "linear")
    {
        if (ok) *ok = true;
        return TabulatedCrossSection::Interpolation::kLinear;
    }
    if (v == "loglog" || v == "log-log" || v == "log")
    {
        if (ok) *ok = true;
        return TabulatedCrossSection::Interpolation::kLogLog;
    }
    if (ok) *ok = false;
    return TabulatedCrossSection::Interpolation::kLinear;
}
}

string DictionaryCrossSections::toLowerCopy(string s)
{
    transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return static_cast<char>(tolower(c)); });
    return s;
}

string DictionaryCrossSections::stripQuotesAndSemicolon(const std::string& s)
{
    std::string out = s;
    if (!out.empty() && out.back() == ';')
        out.pop_back();
    if (out.size() >= 2 && ((out.front() == '"' && out.back() == '"') || (out.front() == '\'' && out.back() == '\'')))
        out = out.substr(1, out.size() - 2);
    return out;
}

bool DictionaryCrossSections::readDict(const std::string& path)
{
    ifstream file(path);
    if (!file.is_open())
        return false;

    struct Entry
    {
        std::string key;
        std::string type;
        std::string library;
        std::string file;
        std::string projectile = "e";
        std::string prefix;
        std::string energyUnit = "eV";
        std::string sigmaUnit = "m2";
        std::string interp = "linear";
        bool overwrite = false;
        bool strict = true;
    };

    bool foundCrossSections = false;
    std::vector<std::string> blocks;
    std::string pendingBlock;

    bool inEntry = false;
    Entry cur;

    auto currentBlock = [&]() -> std::string { return blocks.empty() ? std::string() : blocks.back(); };
    auto parentBlock = [&]() -> std::string {
        if (blocks.size() < 2) return std::string();
        return blocks[blocks.size() - 2];
    };

    auto openBlock = [&](const std::string& name) {
        blocks.push_back(name);
        if (toLowerCopy(name) == "crosssections")
            foundCrossSections = true;
        if (toLowerCopy(parentBlock()) == "crosssections" && toLowerCopy(name) != "crosssections")
        {
            // Start new entry
            inEntry = true;
            cur = Entry{};
            cur.key = stripQuotesAndSemicolon(name);
        }
    };

    auto finalizeEntry = [&]() -> bool {
        if (!inEntry)
            return true;

        const std::string t = toLowerCopy(cur.type);
        if (cur.key.empty())
        {
            std::cerr << "[DictionaryCrossSections] ERROR: empty cross-section key" << std::endl;
            return false;
        }

        if (t.empty())
        {
            std::cerr << "[DictionaryCrossSections] ERROR: CrossSection '" << cur.key << "' missing 'type'" << std::endl;
            return false;
        }

        if (t == "table")
        {
            if (cur.file.empty())
            {
                std::cerr << "[DictionaryCrossSections] ERROR: CrossSection '" << cur.key << "' missing 'file'" << std::endl;
                return false;
            }

            TabulatedCrossSection::Options opt;
            opt.default_energy_unit = cur.energyUnit;
            opt.default_sigma_unit = cur.sigmaUnit;
            bool ok = false;
            opt.interp = parseInterp(cur.interp, &ok);
            if (!ok)
            {
                std::cerr << "[DictionaryCrossSections] ERROR: CrossSection '" << cur.key << "' unknown interp: " << cur.interp << std::endl;
                return false;
            }

            auto xs = std::make_shared<TabulatedCrossSection>(opt);
            std::string err;
            if (!xs->LoadFromFile(cur.file, &err))
            {
                std::cerr << "[DictionaryCrossSections] ERROR: failed to load table '" << cur.file << "' for key '" << cur.key << "': " << err << std::endl;
                return false;
            }
            if (!xs->Validate(&err))
            {
                std::cerr << "[DictionaryCrossSections] ERROR: invalid table for key '" << cur.key << "': " << err << std::endl;
                return false;
            }

            if (!db_.Register(cur.key, xs, cur.overwrite, &err))
            {
                std::cerr << "[DictionaryCrossSections] ERROR: failed to register key '" << cur.key << "': " << err << std::endl;
                return false;
            }
        }
        else if (t == "lxcatfile" || t == "bolsigfile" || t == "lxcat" || t == "bolsig")
        {
            if (cur.file.empty())
            {
                std::cerr << "[DictionaryCrossSections] ERROR: CrossSection source '" << cur.key << "' missing 'file'" << std::endl;
                return false;
            }

            BolsigCrossSectionFile::Options opt;
            opt.projectile = cur.projectile.empty() ? "e" : cur.projectile;
            opt.overwrite = cur.overwrite;
            opt.strict = cur.strict;
            opt.xs_options.default_energy_unit = cur.energyUnit;
            opt.xs_options.default_sigma_unit = cur.sigmaUnit;
            bool ok = false;
            opt.xs_options.interp = parseInterp(cur.interp, &ok);
            if (!ok)
            {
                std::cerr << "[DictionaryCrossSections] ERROR: CrossSection source '" << cur.key << "' unknown interp: " << cur.interp << std::endl;
                return false;
            }
            opt.key_prefix = cur.prefix;

            std::string err;
            std::vector<BolsigCrossSectionFile::BlockInfo> blocks_info;
            if (!BolsigCrossSectionFile::LoadAndRegister(cur.file, &db_, opt, &blocks_info, &err))
            {
                std::cerr << "[DictionaryCrossSections] ERROR: failed to load LXCat/BOLSIG file '" << cur.file << "': " << err << std::endl;
                return false;
            }
        }
        else if (t == "builtin" || t == "built-in" || t == "built_in")
        {
            if (cur.library.empty())
            {
                std::cerr << "[DictionaryCrossSections] ERROR: CrossSection source '" << cur.key << "' missing 'library'" << std::endl;
                return false;
            }

            std::string lib = toLowerCopy(cur.library);
            std::string err;

            if (lib == "arv1" || lib == "ar_v1" || lib == "argonv1" || lib == "argon_v1")
            {
                if (!BuiltInCrossSections::RegisterArV1(&db_, cur.prefix.empty() ? "builtin/Ar/v1:" : cur.prefix,
                                                       cur.overwrite, &err))
                {
                    std::cerr << "[DictionaryCrossSections] ERROR: failed to register built-in library '" << cur.library
                              << "': " << err << std::endl;
                    return false;
                }
            }
            else
            {
                std::cerr << "[DictionaryCrossSections] ERROR: unknown built-in library: " << cur.library << std::endl;
                return false;
            }
        }
        else
        {
            std::cerr << "[DictionaryCrossSections] ERROR: CrossSection '" << cur.key << "' unknown type: " << cur.type << std::endl;
            return false;
        }

        inEntry = false;
        return true;
    };

    auto closeBlock = [&]() -> bool {
        // If we are about to pop an entry directly under CrossSections, finalize it.
        if (inEntry && toLowerCopy(parentBlock()) == "crosssections" && !blocks.empty())
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

        // Block titles waiting for '{'
        if (toLowerCopy(line) == "crosssections")
        {
            pendingBlock = line;
            continue;
        }

        // Key-value inside an entry.
        if (!inEntry)
        {
            // Single token could be a block name (entry name) waiting for next "{"
            const vector<string> singleTokenCheck = split(line);
            if (singleTokenCheck.size() == 1)
            {
                string stripped = stripQuotesAndSemicolon(singleTokenCheck[0]);
                if (!stripped.empty())
                {
                    pendingBlock = stripped;
                }
            }
            continue;
        }
        if (toLowerCopy(parentBlock()) != "crosssections")
            continue;

        const vector<string> tokens = split(line);
        if (tokens.size() < 2)
            continue;

        const string k = toLowerCopy(tokens[0]);
        const string v = stripQuotesAndSemicolon(tokens[1]);

        if (k == "type")
            cur.type = v;
        else if (k == "file" || k == "path")
            cur.file = stripQuotesAndSemicolon(joinTokens(tokens, 1));
        else if (k == "interp" || k == "interpolation")
            cur.interp = v;
        else if (k == "energyunit" || k == "eunit")
            cur.energyUnit = v;
        else if (k == "sigmaunit" || k == "sunit")
            cur.sigmaUnit = v;
        else if (k == "projectile")
            cur.projectile = v;
        else if (k == "library")
            cur.library = v;
        else if (k == "prefix")
            cur.prefix = stripQuotesAndSemicolon(tokens[1]);
        else if (k == "library")
            cur.library = stripQuotesAndSemicolon(joinTokens(tokens, 1));
        else if (k == "overwrite")
        {
            bool b = false;
            if (parseBool(v, b))
                cur.overwrite = b;
            else
                cerr << "[DictionaryCrossSections] WARNING: invalid bool for overwrite: " << line << endl;
        }
        else if (k == "strict")
        {
            bool b = true;
            if (parseBool(v, b))
                cur.strict = b;
            else
                cerr << "[DictionaryCrossSections] WARNING: invalid bool for strict: " << line << endl;
        }
    }

    // Close remaining blocks (in case file ends without braces).
    while (!blocks.empty())
    {
        if (!closeBlock())
            return false;
    }

    if (!foundCrossSections)
    {
        cerr << "[DictionaryCrossSections] ERROR: missing CrossSections{ } block" << endl;
        return false;
    }
    return true;
}
