#include "Runtime/IO/DictionaryRun.h"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <regex>
#include "Core/IO/StrictDictionary.h"
#include "Solvers/IO/KVBlockParser.h"

using namespace std;

string DictionaryRun::toLowerCopy(string s)
{
    transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return static_cast<char>(tolower(c)); });
    return s;
}

bool DictionaryRun::readDict(const string& runDictPath)
{
    Solvers::IO::KVBlockParser parser;
    string text;
    if (!parser.loadFile(runDictPath, text))
        return false;

    config_.app.clear();
    config_.schemaVersion = 1;
    if (std::regex_search(text, std::regex(R"(\bschemaVersion\b)")))
    {
        DictionaryInput::Document versioned(runDictPath);
        const auto& header = versioned.root().block("Run");
        config_.schemaVersion = header.get("schemaVersion").integer();
        if (config_.schemaVersion != 1 && config_.schemaVersion != 2)
            header.get("schemaVersion").fail("1 (legacy) or 2 (versioned Electrospray)");
        if (config_.schemaVersion == 2)
        {
            header.allow({"schemaVersion", "app"});
            config_.app = header.get("app").choice({"Electrospray"});
            return true;
        }
    }
    string runBlock;
    size_t runPos = 0;
    while ((runPos = text.find("Run", runPos)) != string::npos)
    {
        const size_t afterName = runPos + 3;
        const bool leftBoundary = runPos == 0 || std::isspace(static_cast<unsigned char>(text[runPos - 1]))
            || text[runPos - 1] == '}';
        size_t brace = afterName;
        while (brace < text.size() && std::isspace(static_cast<unsigned char>(text[brace]))) ++brace;
        if (leftBoundary && brace < text.size() && text[brace] == '{'
            && parser.extractBlockText(text.substr(runPos), "Run", runBlock))
        {
            // Handles both Run { app SingleModule } and the usual multiline block.
            const auto values = parser.parseKeyValues(runBlock);
            config_.app = parser.getString(values, "app", "");
            return true;
        }
        runPos = afterName;
    }
    // Preserve the historical brace-less Run/app records and top-level app form.
    istringstream file(text);

    string line;
    bool inRunBlock = false;

    while (getline(file, line))
    {
        line = removeInlineComments(line);
        line = trim(line);

        if (line.empty() || line == "{")
            continue;

        if (line == "Run")
        {
            inRunBlock = true;
            continue;
        }

        if (line == "}")
        {
            if (inRunBlock)
            {
                inRunBlock = false;
            }
            continue;
        }

        // --------------------------------------------------------------------
        // Run block
        // --------------------------------------------------------------------
        if (inRunBlock)
        {
            vector<string> tokens = split(line);
            if (tokens.size() >= 2 && toLowerCopy(tokens[0]) == "app")
                config_.app = tokens[1];
            continue;
        }
        // Fallback: allow single-line syntax without blocks.
        vector<string> tokens = split(line);
        if (tokens.size() >= 2 && toLowerCopy(tokens[0]) == "app")
            config_.app = tokens[1];
    }

    return true;
}
