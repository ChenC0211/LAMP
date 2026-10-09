#include "Electrospray/Compatibility/DictionaryElectrospray.h"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <iostream>
#include <stack>
#include <stdexcept>

using namespace std;

string DictionaryElectrospray::toLowerCopy(string s)
{
    transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return static_cast<char>(tolower(c)); });
    return s;
}

bool DictionaryElectrospray::parseBoolToken(const string& s, bool& out)
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

bool DictionaryElectrospray::readDict(const string& path)
{
    ifstream file(path);
    if (!file.is_open())
    {
        // run.dict not found: caller decides if fatal.
        return false;
    }

    bool foundElectrosprayBlock = false;

    stack<string> blocks;
    string pendingBlock;

    auto openBlock = [&](const string& name) {
        blocks.push(name);
        if (toLowerCopy(name) == "electrospray")
            foundElectrosprayBlock = true;
    };
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
        {
            const string key = toLowerCopy(line);
            if (key == "electrospray")
            {
                pendingBlock = line;
                continue;
            }
        }

        const string b = toLowerCopy(currentBlock());
        const vector<string> tokens = split(line);
        if (tokens.size() < 2)
            continue;

        if (b == "electrospray")
        {
            const string k = toLowerCopy(tokens[0]);
            const string v = tokens[1];

            if (k == "outputname")
                cfg_.outputName = v;
            else if (k == "physicsmode")
                cfg_.physicsMode = v;
            else if (k == "capillaryhalfwidthcells")
                cfg_.capillaryHalfWidthCells = stod(v);
            else if (k == "writefields")
                parseBoolToken(v, cfg_.writeFields);
            else if (k == "logevery")
            {
                cfg_.logEvery = stoi(v);
                if (cfg_.logEvery < 1) throw std::invalid_argument("Electrospray logEvery must be positive");
            }
            else if (k == "verbosediagnostics")
            {
                if (!parseBoolToken(v, cfg_.verboseDiagnostics))
                    throw std::invalid_argument("Electrospray verboseDiagnostics must be boolean");
            }
            else if (k == "restartloop")
                cfg_.restartLoop = stoi(v);
            else if (k == "dt")
                cfg_.dt = stod(v);
            else if (k == "finaltime")
                cfg_.finalTime = stod(v);
            else if (k == "writedt")
                cfg_.writeDt = stod(v);
            else if (k == "initialadaptiterations" || k == "initialiterations")
                cfg_.initialAdaptIterations = stoi(v);
            else if (k == "initializelevelsets")
                parseBoolToken(v, cfg_.initializeLevelSets);
        }
    }

    return foundElectrosprayBlock;
}
