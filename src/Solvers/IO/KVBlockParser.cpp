#include "Solvers/IO/KVBlockParser.h"

#include <algorithm>
#include <cctype>
#include <sstream>

using namespace std;

namespace Solvers
{
namespace IO
{

bool KVBlock::has(const string& keyLower) const
{
    return kv.find(keyLower) != kv.end();
}

const vector<string>* KVBlock::get(const string& keyLower) const
{
    auto it = kv.find(keyLower);
    if (it == kv.end())
        return nullptr;
    return &it->second;
}

string KVBlockParser::toLowerCopy(string s)
{
    transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return static_cast<char>(tolower(c)); });
    return s;
}

bool KVBlockParser::parseBoolToken(const string& s, bool& out)
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

string KVBlockParser::stripTrailingPunct(string s)
{
    while (!s.empty())
    {
        const char c = s.back();
        if (c == ';' || c == ',' )
            s.pop_back();
        else
            break;
    }
    return s;
}

vector<string> KVBlockParser::stripParenTokens(const vector<string>& toks)
{
    vector<string> out;
    out.reserve(toks.size());
    for (auto t : toks)
    {
        t = stripTrailingPunct(t);
        if (t == "(" || t == ")")
            continue;
        // Handle tokens like "(A" or "B)".
        if (!t.empty() && t.front() == '(')
            t.erase(t.begin());
        if (!t.empty() && t.back() == ')')
            t.pop_back();
        if (!t.empty())
            out.push_back(t);
    }
    return out;
}

KVBlock KVBlockParser::parseKeyValues(const string& blockText) const
{
    KVBlock out;
    istringstream iss(blockText);
    string line;

    while (getline(iss, line))
    {
        line = trim(removeInlineComments(line));
        if (line.empty())
            continue;

        vector<string> toks = split(line);
        if (toks.size() < 2)
            continue;

        const string key = toLowerCopy(stripTrailingPunct(toks[0]));
        vector<string> vals(toks.begin() + 1, toks.end());
        vals = stripParenTokens(vals);
        if (vals.empty())
            continue;

        out.kv[key] = vals;
    }

    return out;
}

bool KVBlockParser::getBool(const KVBlock& block, const string& key, bool defaultVal) const
{
    const string k = toLowerCopy(key);
    const vector<string>* vals = block.get(k);
    if (!vals || vals->empty())
        return defaultVal;

    bool v = defaultVal;
    if (parseBoolToken((*vals)[0], v))
        return v;
    return defaultVal;
}

int KVBlockParser::getInt(const KVBlock& block, const string& key, int defaultVal) const
{
    const string k = toLowerCopy(key);
    const vector<string>* vals = block.get(k);
    if (!vals || vals->empty())
        return defaultVal;
    try
    {
        return stoi((*vals)[0]);
    }
    catch (...) { return defaultVal; }
}

double KVBlockParser::getDouble(const KVBlock& block, const string& key, double defaultVal) const
{
    const string k = toLowerCopy(key);
    const vector<string>* vals = block.get(k);
    if (!vals || vals->empty())
        return defaultVal;
    try
    {
        return stod((*vals)[0]);
    }
    catch (...) { return defaultVal; }
}

string KVBlockParser::getString(const KVBlock& block, const string& key, const string& defaultVal) const
{
    const string k = toLowerCopy(key);
    const vector<string>* vals = block.get(k);
    if (!vals || vals->empty())
        return defaultVal;
    return (*vals)[0];
}

vector<string> KVBlockParser::getList(const KVBlock& block, const string& key) const
{
    const string k = toLowerCopy(key);
    const vector<string>* vals = block.get(k);
    if (!vals)
        return {};
    return *vals;
}

} // namespace IO
} // namespace Solvers
