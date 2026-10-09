#include "Modules/Plasma/Physics/BolsigCrossSectionFile.h"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <sstream>
#include <string>

namespace
{

std::string Trim(const std::string& s)
{
    std::size_t b = 0;
    while (b < s.size() && std::isspace(static_cast<unsigned char>(s[b])))
    {
        ++b;
    }
    std::size_t e = s.size();
    while (e > b && std::isspace(static_cast<unsigned char>(s[e - 1])))
    {
        --e;
    }
    return s.substr(b, e - b);
}

bool StartsWith(const std::string& s, const std::string& prefix)
{
    if (s.size() < prefix.size())
    {
        return false;
    }
    return s.compare(0, prefix.size(), prefix) == 0;
}

std::string ToLower(std::string s)
{
    for (char& c : s)
    {
        c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    }
    return s;
}

std::string StripCR(const std::string& s)
{
    if (!s.empty() && s.back() == '\r')
    {
        return s.substr(0, s.size() - 1);
    }
    return s;
}

std::string StripComment(const std::string& line)
{
    std::size_t cut = std::string::npos;
    auto consider = [&](std::size_t p)
    {
        if (p != std::string::npos)
        {
            cut = (cut == std::string::npos) ? p : std::min(cut, p);
        }
    };
    consider(line.find('#'));
    consider(line.find(';'));
    consider(line.find("//"));
    if (cut == std::string::npos)
    {
        return line;
    }
    return line.substr(0, cut);
}

bool IsDashLine(const std::string& line)
{
    const std::string t = Trim(line);
    if (t.size() < 5)
    {
        return false;
    }
    for (char c : t)
    {
        if (c != '-')
        {
            return false;
        }
    }
    return true;
}

bool IsKeyword(const std::string& s)
{
    static const std::vector<std::string> k = {
        "ELASTIC",
        "EFFECTIVE",
        "MOMENTUM",
        "IONIZATION",
        "ATTACHMENT",
        "EXCITATION",
        "ROTATION",
    };
    return std::find(k.begin(), k.end(), s) != k.end();
}

std::vector<std::string> Tokenize(const std::string& line)
{
    std::vector<std::string> toks;
    std::istringstream iss(line);
    std::string t;
    while (iss >> t)
    {
        toks.push_back(t);
    }
    return toks;
}

bool TryParseDouble(const std::string& token, double& v)
{
    try
    {
        size_t idx = 0;
        v = std::stod(token, &idx);
        return idx == token.size();
    }
    catch (...)
    {
        return false;
    }
}

bool TryParseTwoDoubles(const std::string& line, double& a, double& b)
{
    const auto toks = Tokenize(line);
    if (toks.size() < 2)
    {
        return false;
    }
    return TryParseDouble(toks[0], a) && TryParseDouble(toks[1], b);
}

void ParseTargetAndState(const std::string& line, std::string& target, std::string& state)
{
    // Handle optional "->" or "<->". BOLSIG+ ignores arrow and everything after,
    // but we keep the RHS as an optional state label for DB keys.
    const std::string s = Trim(line);
    const std::size_t p2 = s.find("<->");
    const std::size_t p1 = s.find("->");
    std::size_t p = std::string::npos;
    std::size_t arrow_len = 0;
    if (p2 != std::string::npos)
    {
        p = p2;
        arrow_len = 3;
    }
    else if (p1 != std::string::npos)
    {
        p = p1;
        arrow_len = 2;
    }

    if (p == std::string::npos)
    {
        target = Trim(s);
        state.clear();
        return;
    }

    target = Trim(s.substr(0, p));
    state = Trim(s.substr(p + arrow_len));
}

std::string SanitizeKeyToken(const std::string& s)
{
    std::string out;
    out.reserve(s.size());
    for (char c : s)
    {
        if (std::isalnum(static_cast<unsigned char>(c)) || c == '_' || c == '-' || c == '.')
        {
            out.push_back(c);
        }
        else if (std::isspace(static_cast<unsigned char>(c)))
        {
            out.push_back('_');
        }
        else
        {
            out.push_back('_');
        }
    }
    // Trim underscores.
    while (!out.empty() && out.front() == '_') out.erase(out.begin());
    while (!out.empty() && out.back() == '_') out.pop_back();
    return out;
}

std::string BuildKey(const std::string& projectile,
                     const std::string& target,
                     const std::string& type,
                     const std::string& state)
{
    std::string key = SanitizeKeyToken(projectile) + "-" + SanitizeKeyToken(target) + ":" + type;
    if (!state.empty())
    {
        key += ":" + SanitizeKeyToken(state);
    }
    return key;
}

bool ParsePhelpsSpeciesLine(const std::string& line, std::string& projectile, std::string& target)
{
    if (!StartsWith(line, "SPECIES:"))
    {
        return false;
    }

    const std::string rest = Trim(line.substr(8));
    projectile.clear();
    target.clear();

    const std::size_t slash = rest.find('/');
    if (slash != std::string::npos)
    {
        projectile = Trim(rest.substr(0, slash));
        target = Trim(rest.substr(slash + 1));
    }
    else
    {
        target = rest;
    }
    return !target.empty();
}

std::string ExtractProcessState(const std::string& processText)
{
    const std::size_t comma = processText.rfind(',');
    if (comma != std::string::npos)
    {
        const std::string tail = Trim(processText.substr(comma + 1));
        if (!tail.empty())
        {
            return tail;
        }
    }

    const std::string low = ToLower(processText);
    if (low.find("backscat") != std::string::npos)
        return "Backscat";
    if (low.find("isotropic") != std::string::npos)
        return "Isotropic";

    return std::string();
}

std::string GuessProcessType(const std::string& processText)
{
    const std::string low = ToLower(processText);
    if (low.find("ionization") != std::string::npos)
        return "IONIZATION";
    if (low.find("excitation") != std::string::npos)
        return "EXCITATION";
    if (low.find("attachment") != std::string::npos)
        return "ATTACHMENT";
    if (low.find("effective") != std::string::npos)
        return "EFFECTIVE";
    if (low.find("elastic") != std::string::npos)
        return "ELASTIC";
    return "ELASTIC";
}

} // namespace

bool BolsigCrossSectionFile::LoadAndRegister(const std::string& file_path,
                                             CrossSectionDB* db,
                                             const Options& options,
                                             std::vector<BlockInfo>* out_blocks,
                                             std::string* err)
{
    if (!db)
    {
        if (err) *err = "CrossSectionDB pointer is null";
        return false;
    }

    std::ifstream ifs(file_path);
    if (!ifs.is_open())
    {
        if (err) *err = "Failed to open file: " + file_path;
        return false;
    }

    std::vector<BlockInfo> blocks;
    std::string line;
    std::size_t lineNo = 0;
    int phelpsIndex = 0;

    auto fail_or_skip = [&](const std::string& msg) -> bool
    {
        if (options.strict)
        {
            if (err) *err = msg;
            return false;
        }
        // non-strict: skip and continue
        return true;
    };

    auto parsePhelpsBlock = [&](const std::string& speciesLine) -> bool
    {
        std::string projectile;
        std::string target;
        if (!ParsePhelpsSpeciesLine(speciesLine, projectile, target))
        {
            return fail_or_skip("Line " + std::to_string(lineNo) + ": failed to parse SPECIES line");
        }

        BlockInfo info;
        info.type = "ELASTIC";
        info.target = target;
        info.state.clear();
        info.comments.push_back(speciesLine);

        // Find PROCESS line.
        std::string processText;
        std::string line;
        while (std::getline(ifs, line))
        {
            ++lineNo;
            line = Trim(StripCR(line));
            if (line.empty())
                continue;
            if (StartsWith(line, "PROCESS:"))
            {
                processText = Trim(line.substr(8));
                info.comments.push_back(line);
                break;
            }
            info.comments.push_back(line);
        }
        if (processText.empty())
        {
            return fail_or_skip("Line " + std::to_string(lineNo) + ": missing PROCESS line for Phelps block");
        }

        info.type = GuessProcessType(processText);
        info.state = ExtractProcessState(processText);
        if (info.state.empty())
        {
            ++phelpsIndex;
            info.state = "PhelpsBlock" + std::to_string(phelpsIndex);
        }

        // Scan forward to the opening dash line.
        bool foundDashStart = false;
        while (std::getline(ifs, line))
        {
            ++lineNo;
            line = Trim(StripCR(line));
            if (line.empty())
                continue;
            if (IsDashLine(line))
            {
                foundDashStart = true;
                break;
            }
            info.comments.push_back(line);
        }
        if (!foundDashStart)
        {
            return fail_or_skip("EOF while searching for table start after Phelps block");
        }

        // Read table lines until closing dash line.
        std::vector<std::string> tableLines;
        bool foundDashEnd = false;
        while (std::getline(ifs, line))
        {
            ++lineNo;
            line = Trim(StripCR(line));
            if (line.empty())
                continue;
            if (IsDashLine(line))
            {
                foundDashEnd = true;
                break;
            }
            line = Trim(StripComment(line));
            if (!line.empty())
            {
                tableLines.push_back(line);
            }
        }
        if (!foundDashEnd)
        {
            return fail_or_skip("EOF while reading table for Phelps block");
        }
        if (tableLines.empty())
        {
            return fail_or_skip("Empty cross-section table for Phelps block");
        }

        auto xs = std::make_shared<TabulatedCrossSection>(options.xs_options);
        std::string localErr;

        std::vector<std::pair<double, double>> pairs;
        pairs.reserve(tableLines.size());
        for (const auto& l : tableLines)
        {
            double E = 0.0;
            double S = 0.0;
            if (!TryParseTwoDoubles(l, E, S))
            {
                if (!fail_or_skip("Failed to parse table line '" + l + "'"))
                {
                    return false;
                }
                pairs.clear();
                break;
            }
            pairs.emplace_back(E, S);
        }
        if (pairs.empty())
        {
            return true;
        }
        std::sort(pairs.begin(), pairs.end(), [](auto a, auto b) { return a.first < b.first; });

        std::vector<double> Es;
        std::vector<double> Ss;
        Es.reserve(pairs.size());
        Ss.reserve(pairs.size());
        for (const auto& pr : pairs)
        {
            Es.push_back(pr.first);
            Ss.push_back(pr.second);
        }

        if (!xs->LoadFromColumns(Es, Ss, &localErr))
        {
            return fail_or_skip("Invalid table in Phelps block: " + localErr);
        }

        const std::string key = options.key_prefix + BuildKey(options.projectile, info.target, info.type, info.state);
        info.key = key;

        if (!db->Register(key, xs, options.overwrite, &localErr))
        {
            return fail_or_skip("Failed to register key '" + key + "': " + localErr);
        }

        blocks.push_back(info);
        return true;
    };

    while (std::getline(ifs, line))
    {
        ++lineNo;
        line = Trim(StripCR(line));
        if (line.empty())
        {
            continue;
        }

        const std::string type = line;
        if (!IsKeyword(type))
        {
            // Support LXCat "Phelps" ion blocks (SPECIES / PROCESS / dashed table).
            if (StartsWith(type, "SPECIES:"))
            {
                if (!parsePhelpsBlock(type))
                {
                    return false;
                }
                continue;
            }
            // Skip non-keyword lines (headers/comments).
            continue;
        }

        // Line 2: target (optionally with arrow/state)
        std::string targetLine;
        while (std::getline(ifs, targetLine))
        {
            ++lineNo;
            targetLine = Trim(StripCR(targetLine));
            if (!targetLine.empty())
            {
                break;
            }
        }
        if (targetLine.empty())
        {
            if (err) *err = "Unexpected EOF after keyword '" + type + "'";
            return false;
        }

        std::string target, state;
        ParseTargetAndState(targetLine, target, state);
        if (target.empty())
        {
            if (err) *err = "Line " + std::to_string(lineNo) + ": empty target name";
            return false;
        }

        // Line 3: parameter line (type-dependent). May be omitted for ATTACHMENT.
        std::string paramLine;
        std::streampos posAfterTarget = ifs.tellg();
        std::size_t lineNoAfterTarget = lineNo;
        while (std::getline(ifs, paramLine))
        {
            ++lineNo;
            paramLine = Trim(StripCR(paramLine));
            if (!paramLine.empty())
            {
                break;
            }
        }
        if (paramLine.empty())
        {
            if (err) *err = "Unexpected EOF after target '" + target + "'";
            return false;
        }

        BlockInfo info;
        info.type = type;
        info.target = target;
        info.state = state;

        double paramVal = 0.0;
        // If the next non-empty line is a dash line, then param line is absent.
        if (IsDashLine(paramLine))
        {
            // param line absent; dash line is the table start
            // param line absent
            // Rewind to treat this as the dash line.
            ifs.seekg(posAfterTarget);
            lineNo = lineNoAfterTarget;
        }
        else
        {
            // Try parse first token as double.
            const auto toks = Tokenize(StripComment(paramLine));
            if (!toks.empty() && TryParseDouble(toks[0], paramVal))
            {
                info.parameter = paramVal;
                // Threshold for inelastic processes.
                if (type == "IONIZATION" || type == "EXCITATION" || type == "ROTATION")
                {
                    info.threshold_eV = paramVal;
                }
            }
            else
            {
                // Not a numeric param; treat as comment and rewind so we don't lose it.
                ifs.seekg(posAfterTarget);
                lineNo = lineNoAfterTarget;
                // treat as comment
            }
        }

        // Now scan forward to the opening dash line, collecting comments.
        bool foundDashStart = false;
        while (std::getline(ifs, line))
        {
            ++lineNo;
            line = Trim(StripCR(line));
            if (line.empty())
            {
                continue;
            }
            if (IsDashLine(line))
            {
                foundDashStart = true;
                break;
            }
            info.comments.push_back(line);
        }
        if (!foundDashStart)
        {
            if (err) *err = "EOF while searching for table start after '" + type + " " + target + "'";
            return false;
        }

        // Read table lines until closing dash line.
        std::vector<std::string> tableLines;
        bool foundDashEnd = false;
        while (std::getline(ifs, line))
        {
            ++lineNo;
            line = Trim(StripCR(line));
            if (line.empty())
            {
                continue;
            }
            if (IsDashLine(line))
            {
                foundDashEnd = true;
                break;
            }
            line = Trim(StripComment(line));
            if (!line.empty())
            {
                tableLines.push_back(line);
            }
        }
        if (!foundDashEnd)
        {
            if (err) *err = "EOF while reading table for '" + type + " " + target + "'";
            return false;
        }
        if (tableLines.empty())
        {
            if (err) *err = "Empty cross-section table for '" + type + " " + target + "'";
            return false;
        }

        auto xs = std::make_shared<TabulatedCrossSection>(options.xs_options);
        std::string localErr;

        // Power-law shorthand: single data line "p sigma1".
        if (tableLines.size() == 1)
        {
            double p = 0.0;
            double sigma1 = 0.0;
            if (!TryParseTwoDoubles(tableLines[0], p, sigma1))
            {
                if (!fail_or_skip("Line " + std::to_string(lineNo) + ": failed to parse power-law line"))
                {
                    return false;
                }
                continue;
            }
            if (!xs->SetPowerLaw(p, sigma1, &localErr))
            {
                if (!fail_or_skip("Invalid power-law in '" + type + " " + target + "': " + localErr))
                {
                    return false;
                }
                continue;
            }
        }
        else
        {
            std::vector<std::pair<double, double>> pairs;
            pairs.reserve(tableLines.size());
            for (const auto& l : tableLines)
            {
                double E = 0.0;
                double S = 0.0;
                if (!TryParseTwoDoubles(l, E, S))
                {
                    if (!fail_or_skip("Failed to parse table line '" + l + "'"))
                    {
                        return false;
                    }
                    pairs.clear();
                    break;
                }
                pairs.emplace_back(E, S);
            }
            if (pairs.empty())
            {
                continue;
            }
            std::sort(pairs.begin(), pairs.end(), [](auto a, auto b) { return a.first < b.first; });

            std::vector<double> Es;
            std::vector<double> Ss;
            Es.reserve(pairs.size());
            Ss.reserve(pairs.size());
            for (const auto& pr : pairs)
            {
                Es.push_back(pr.first);
                Ss.push_back(pr.second);
            }

            if (!xs->LoadFromColumns(Es, Ss, &localErr))
            {
                if (!fail_or_skip("Invalid table in '" + type + " " + target + "': " + localErr))
                {
                    return false;
                }
                continue;
            }
        }

        const std::string key = options.key_prefix + BuildKey(options.projectile, target, type, state);
        info.key = key;

        if (!db->Register(key, xs, options.overwrite, &localErr))
        {
            if (!fail_or_skip("Failed to register key '" + key + "': " + localErr))
            {
                return false;
            }
            continue;
        }

        blocks.push_back(info);
    }

    if (out_blocks)
    {
        *out_blocks = blocks;
    }
    return true;
}
