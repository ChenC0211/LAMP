#include "Electricity/Compatibility/DictionaryElectricity.h"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <stack>

using namespace std;

string DictionaryElectricity::toLowerCopy(string s)
{
    transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return static_cast<char>(tolower(c)); });
    return s;
}

bool DictionaryElectricity::parseBoolToken(const string& s, bool& out)
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

static bool isQuoted(const string& s)
{
    return s.size() >= 2 && ((s.front() == '"' && s.back() == '"') || (s.front() == '\'' && s.back() == '\''));
}

static string unquote(const string& s)
{
    return isQuoted(s) ? s.substr(1, s.size() - 2) : s;
}

bool DictionaryElectricity::readDict(const string& path)
{
    ifstream file(path);
    if (!file.is_open())
    {
        return false;
    }

    stack<string> blocks;
    string pendingBlock;

    auto openBlock = [&](const string& name) { blocks.push(name); };
    auto closeBlock = [&]() {
        if (!blocks.empty())
            blocks.pop();
    };
    auto currentBlock = [&]() -> string { return blocks.empty() ? string() : blocks.top(); };

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
            closeBlock();
            continue;
        }

        // Inline "NAME {" opening
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

        // Block title waiting for '{'
        {
            const string key = toLowerCopy(line);
            if (key == "electricity" || key == "output")
            {
                pendingBlock = line;
                continue;
            }
        }

        const string b = toLowerCopy(currentBlock());
        const vector<string> tokens = split(line);
        if (tokens.size() < 2)
            continue;

        const string k = toLowerCopy(tokens[0]);
        const string v = tokens[1];

        // Accept keys both in Electricity{} and Output{}
        if (b == "electricity" || b == "output" || b.empty())
        {
            if (k == "outputname" || k == "output" || k == "name")
            {
                cfg_.outputName = unquote(v);
            }
            else if (k == "writefields" || k == "write" || k == "dump")
            {
                parseBoolToken(v, cfg_.writeFields);
            }
            else if (k == "computehydroproperties" || k == "computehydro" || k == "hydroproperties")
            {
                parseBoolToken(v, cfg_.computeHydroProperties);
            }
            else if (k == "propertylevelsetindex" || k == "levelsetindex" || k == "levelset")
            {
                cfg_.propertyLevelSetIndex = stoi(v);
            }
        }
    }

    return true;
}
