#include "Core/IO/DictionaryBlockText.h"
#include <cctype>

namespace
{
bool space(char ch) { return std::isspace(static_cast<unsigned char>(ch)) != 0; }

// Advance over one quoted literal or a legacy line comment. Keep the original
// offsets so extracted bodies retain their whitespace and nested structure.
bool skipLiteral(const std::string& text, std::size_t& pos)
{
    const char ch = text[pos];
    if (ch == '"' || ch == '\'')
    {
        ++pos;
        while (pos < text.size())
        {
            if (text[pos] == '\\' && pos + 1 < text.size()) pos += 2;
            else if (text[pos++] == ch) break;
        }
        return true;
    }
    if (ch == '#' || (ch == '/' && pos + 1 < text.size() && text[pos + 1] == '/'
                     && (pos == 0 || space(text[pos - 1]))))
    {
        const auto end = text.find('\n', pos);
        pos = end == std::string::npos ? text.size() : end;
        return true;
    }
    return false;
}

bool bodyAt(const std::string& text, std::size_t open, std::string& body)
{
    int depth = 1;
    for (std::size_t pos = open + 1; pos < text.size();)
    {
        if (skipLiteral(text, pos)) continue;
        if (text[pos] == '{') ++depth;
        else if (text[pos] == '}' && --depth == 0)
        {
            body = text.substr(open + 1, pos - open - 1);
            return true;
        }
        ++pos;
    }
    return false;
}
}

bool DictionaryBlockText::loadFile(const std::string& path, std::string& outText)
{
    std::ifstream file(path);
    if (!file.is_open()) return false;
    std::ostringstream buffer;
    std::string line;
    while (std::getline(file, line)) buffer << removeInlineComments(line) << '\n';
    outText = buffer.str();
    return true;
}

bool DictionaryBlockText::extractBlockText(const std::string& text, const std::string& name,
    std::string& body, Scope scope) const
{
    body.clear();
    if (name.empty()) return false;
    int depth = 0, listDepth = 0;
    bool recordStart = true;
    // A UTF-8 BOM is a file marker, not part of the first block name.
    const std::size_t first = text.compare(0, 3, "\xEF\xBB\xBF") == 0 ? 3 : 0;
    for (std::size_t pos = first; pos < text.size();)
    {
        const char ch = text[pos];
        if (space(ch))
        {
            if (ch == '\n' || ch == '\r') recordStart = true;
            ++pos;
            continue;
        }
        if (skipLiteral(text, pos)) { recordStart = false; continue; }
        if (ch == '{') { ++depth; recordStart = true; ++pos; continue; }
        if (ch == '}')
        {
            if (depth == 0) return false;
            --depth; recordStart = true; ++pos; continue;
        }
        if (ch == '(') { ++listDepth; recordStart = false; ++pos; continue; }
        if (ch == ')') { if (listDepth > 0) --listDepth; recordStart = false; ++pos; continue; }
        if (ch == ';') { recordStart = true; ++pos; continue; }
        const auto start = pos;
        while (pos < text.size() && !space(text[pos]) &&
               std::string("{}();\"'").find(text[pos]) == std::string::npos) ++pos;
        if (recordStart && listDepth == 0 && (scope == Scope::AnyDepth || depth == 0) &&
            text.compare(start, pos - start, name) == 0)
        {
            auto open = pos;
            while (open < text.size())
            {
                if (space(text[open])) { ++open; continue; }
                // Only comments may intervene between a name and its brace.
                if (text[open] != '"' && text[open] != '\'' && skipLiteral(text, open)) continue;
                break;
            }
            if (open < text.size() && text[open] == '{') return bodyAt(text, open, body);
        }
        recordStart = false;
    }
    return false;
}

bool DictionaryBlockText::extractFirstBlockText(const std::string& text,
    const std::vector<std::string>& names, std::string& body, Scope scope) const
{
    body.clear();
    for (const auto& name : names)
        if (extractBlockText(text, name, body, scope)) return true;
    return false;
}
