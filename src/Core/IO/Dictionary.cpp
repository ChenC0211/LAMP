#include "Core/IO/Dictionary.h"

#include "Core/Units/UnitParser.h"
#include "Core/Units/UnitsContext.h"

#include <algorithm>
#include <cctype>

using namespace std;

Dictionary::Dictionary() = default;
Dictionary::~Dictionary() = default;

string Dictionary::trim(const string& line) const
{
    const size_t start = line.find_first_not_of(" \t\r\n");
    const size_t end = line.find_last_not_of(" \t\r\n");
    return (start == string::npos) ? string() : line.substr(start, end - start + 1);
}

string Dictionary::removeInlineComments(const string& line) const
{
    // Support both '#' and '//' style inline comments.
    //
    // We only treat '//' as a comment starter if it is at the beginning of
    // the line or preceded by whitespace. This avoids cutting URLs like
    // "http://...".

    const size_t hashPos = line.find('#');

    size_t slashPos = std::string::npos;
    {
        const size_t p = line.find("//");
        if (p != std::string::npos)
        {
            if (p == 0 || std::isspace(static_cast<unsigned char>(line[p - 1])) != 0)
                slashPos = p;
        }
    }

    size_t commentPos = std::string::npos;
    if (hashPos != std::string::npos && slashPos != std::string::npos)
        commentPos = std::min(hashPos, slashPos);
    else if (hashPos != std::string::npos)
        commentPos = hashPos;
    else
        commentPos = slashPos;

    if (commentPos != std::string::npos)
        return line.substr(0, commentPos);
    return line;
}

vector<string> Dictionary::split(const string& str) const
{
    istringstream stream(str);
    string token;
    vector<string> tokens;
    while (stream >> token)
        tokens.push_back(token);
    return tokens;
}

vector<string> Dictionary::split(const string& str, const string& delimiter) const
{
    vector<string> tokens;
    size_t start = 0;
    size_t end = str.find(delimiter);

    while (end != string::npos)
    {
        tokens.push_back(str.substr(start, end - start));
        start = end + delimiter.length();
        end = str.find(delimiter, start);
    }

    tokens.push_back(str.substr(start));
    return tokens;
}

bool Dictionary::parseCoord(const string& line, Vector2D& coord) const
{
    const regex coordRegex(R"(\(\s*(-?[\d.]+(?:[eE][+-]?[\d]+)?)\s*,\s*(-?[\d.]+(?:[eE][+-]?[\d]+)?)\s*\))");
    smatch match;
    if (regex_search(line, match, coordRegex) && match.size() == 3)
    {
        coord = Vector2D(stod(match[1].str()), stod(match[2].str()));
        return true;
    }
    return false;
}

namespace
{
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

    // Ensure there are no trailing non-space characters.
    while (idx < s.size() && std::isspace(static_cast<unsigned char>(s[idx])) != 0)
        ++idx;
    return idx == s.size();
}

bool parseLengthToCode(const std::string& token, const std::string& trailing_unit, const UnitsContext* units, double& out)
{
    // 1) Strict numeric: use trailing unit if provided.
    double numeric = 0.0;
    if (tryParseStrictDouble(token, numeric))
    {
        if (!trailing_unit.empty())
        {
            UnitParser::Result r;
            if (!UnitParser::TryParseToSI(std::to_string(numeric) + " " + trailing_unit, r))
                return false;
            const std::string u = r.unit;
            if (!(u == "m" || u == "cm" || u == "mm" || u == "um" || u == "nm" || u == "km"))
                return false;
            if (units)
                out = units->LengthSIToCode(r.value_si);
            else
                out = r.value_si;
            return true;
        }

        out = numeric;
        return true;
    }

    // 2) token includes unit.
    UnitParser::Result r;
    if (!UnitParser::TryParseToSI(token, r))
        return false;

    if (!r.has_unit)
    {
        // Treat as raw numeric code value.
        out = r.value_si;
        return true;
    }

    // Length unit check (best-effort; parser already validated unit).
    const std::string u = r.unit;
    if (!(u == "m" || u == "cm" || u == "mm" || u == "um" || u == "nm" || u == "km"))
        return false;

    if (units)
        out = units->LengthSIToCode(r.value_si);
    else
        out = r.value_si;
    return true;
}
} // namespace

bool Dictionary::parseCoordWithUnits(const std::string& line, Vector2D& coord, const UnitsContext* units) const
{
    // Match: (x, y) [unit]
    const std::regex coordRegex(R"(\(\s*([^,]+?)\s*,\s*([^\)]+?)\s*\)\s*([A-Za-zµμu][A-Za-z0-9µμu]*)?)");
    std::smatch m;
    if (!std::regex_search(line, m, coordRegex) || m.size() < 3)
        return false;

    const std::string xs = trim(m[1].str());
    const std::string ys = trim(m[2].str());
    const std::string unit = (m.size() >= 4) ? trim(m[3].str()) : std::string();

    double x = 0.0;
    double y = 0.0;
    if (!parseLengthToCode(xs, unit, units, x))
        return false;
    if (!parseLengthToCode(ys, unit, units, y))
        return false;

    coord = Vector2D(x, y);
    return true;
}

bool Dictionary::parseIndex(const string& line, vector<SignInt>& index) const
{
    // Simple parser without regex to avoid memory issues with long lines (e.g., 96+ numbers)
    size_t start = line.find('(');
    size_t end = line.rfind(')');
    
    if (start == string::npos || end == string::npos || end <= start)
        return false;
    
    string content = line.substr(start + 1, end - start - 1);
    
    size_t pos = 0;
    while (pos < content.length())
    {
        // Skip whitespace
        while (pos < content.length() && (std::isspace(content[pos]) || content[pos] == ','))
            ++pos;
        
        if (pos >= content.length())
            break;
        
        // Parse number
        size_t numStart = pos;
        bool isNegative = false;
        
        if (content[pos] == '-')
        {
            isNegative = true;
            ++pos;
        }
        
        while (pos < content.length() && std::isdigit(content[pos]))
            ++pos;
        
        if (pos == numStart + (isNegative ? 1 : 0))
            return false; // Invalid token: do not loop forever without consuming input.
        if (pos > numStart + (isNegative ? 1 : 0))
        {
            string numStr = content.substr(numStart, pos - numStart);
            try {
                int value = stoi(numStr);
                index.push_back({abs(value), isNegative});
            } catch (...) {
                return false;
            }
        }
    }
    
    return !index.empty();
}

bool Dictionary::parseIndex_brace(const string& line, vector<SignInt>& index) const
{
    // Simple parser without regex to avoid memory issues with long lines
    size_t start = line.find('{');
    size_t end = line.rfind('}');
    
    if (start == string::npos || end == string::npos || end <= start)
        return false;
    
    string content = line.substr(start + 1, end - start - 1);
    
    size_t pos = 0;
    while (pos < content.length())
    {
        // Skip whitespace
        while (pos < content.length() && (std::isspace(content[pos]) || content[pos] == ','))
            ++pos;
        
        if (pos >= content.length())
            break;
        
        // Parse number
        size_t numStart = pos;
        bool isNegative = false;
        
        if (content[pos] == '-')
        {
            isNegative = true;
            ++pos;
        }
        
        while (pos < content.length() && std::isdigit(content[pos]))
            ++pos;
        
        if (pos == numStart + (isNegative ? 1 : 0))
            return false; // Same progress guarantee for brace-delimited indices.
        if (pos > numStart + (isNegative ? 1 : 0))
        {
            string numStr = content.substr(numStart, pos - numStart);
            try {
                int value = stoi(numStr);
                index.push_back({abs(value), isNegative});
            } catch (...) {
                return false;
            }
        }
    }
    
    return !index.empty();
}

bool Dictionary::parseCoords(const string& line, vector<Vector2D>& coords) const
{
    const regex coordRegex(R"(\(\s*(-?[\d.]+(?:[eE][+-]?[\d]+)?)\s*,\s*(-?[\d.]+(?:[eE][+-]?[\d]+)?)\s*\))");
    auto begin = sregex_iterator(line.begin(), line.end(), coordRegex);
    auto end = sregex_iterator();

    for (auto it = begin; it != end; ++it)
    {
        const smatch m = *it;
        if (m.size() == 3)
            coords.push_back(Vector2D(stod(m[1].str()), stod(m[2].str())));
    }

    return !coords.empty();
}

bool Dictionary::parseCoordsWithUnits(const std::string& line, std::vector<Vector2D>& coords, const UnitsContext* units) const
{
    coords.clear();

    const std::regex coordRegex(R"(\(\s*([^,]+?)\s*,\s*([^\)]+?)\s*\)\s*([A-Za-zµμu][A-Za-z0-9µμu]*)?)");
    auto begin = std::sregex_iterator(line.begin(), line.end(), coordRegex);
    auto end = std::sregex_iterator();

    for (auto it = begin; it != end; ++it)
    {
        const std::smatch m = *it;
        if (m.size() < 3)
            continue;

        const std::string xs = trim(m[1].str());
        const std::string ys = trim(m[2].str());
        const std::string unit = (m.size() >= 4) ? trim(m[3].str()) : std::string();

        double x = 0.0;
        double y = 0.0;
        if (!parseLengthToCode(xs, unit, units, x))
            return false;
        if (!parseLengthToCode(ys, unit, units, y))
            return false;

        coords.push_back(Vector2D(x, y));
    }

    return !coords.empty();
}
