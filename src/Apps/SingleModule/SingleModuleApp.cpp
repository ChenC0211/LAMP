#include "SingleModule/SingleModuleApp.h"
#include "SingleModule/PicModuleAdapter.h"

#include <algorithm>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <vector>

#include "Runtime/ModuleRegistry.h"

#include "Solvers/IO/KVBlockParser.h"

#include "Modules/Plasma/ExplicitPIC/ExplicitPICModule.h"
#include "Modules/Plasma/ImplicitPIC/ImplicitPICModule.h"
#include "Electricity/Compatibility/LegacyElectricityModule.h"
#include "Hydrodynamics/HydrodynamicsApp.h"

using namespace std;

namespace
{
string toLowerCopy(string s)
{
    transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return static_cast<char>(tolower(c)); });
    return s;
}

ModuleRegistry makeDefaultModuleRegistry()
{
    ModuleRegistry registry;

    // Built-in modules
    registry.registerModule("ExplicitPIC", [](const RunConfig& c) { return make_unique<PicModuleAdapter<ExplicitPICModule>>(c); });
    registry.registerModule("ImplicitPIC", [](const RunConfig& c) { return make_unique<PicModuleAdapter<ImplicitPICModule>>(c); });
    registry.registerModule("Electricity", [](const RunConfig& c) { return make_unique<ElectricityModule>(c); });
    registry.registerModule("Hydrodynamics", [](const RunConfig& c) { return make_unique<HydrodynamicsApp>(c); });

    // Optional aliases
    registry.registerModule("PIC", [](const RunConfig& c) { return make_unique<PicModuleAdapter<ExplicitPICModule>>(c); });
    registry.registerModule("iPIC", [](const RunConfig& c) { return make_unique<PicModuleAdapter<ImplicitPICModule>>(c); });
    registry.registerModule("Hydro", [](const RunConfig& c) { return make_unique<HydrodynamicsApp>(c); });

    return registry;
}

} // namespace

std::vector<std::string> SingleModuleApp::readModulesFromRunDict(const std::string& path)
{
    Solvers::IO::KVBlockParser parser;
    string text;
    if (!parser.loadFile(path, text)) return {};

    // Locate the block at top level: the app value in "Run { app SingleModule }"
    // is not a block name. Reject duplicate blocks and unbalanced braces.
    string body;
    int depth = 0;
    int blocks = 0;
    for (size_t pos = 0; pos < text.size();)
    {
        const char c = text[pos];
        if (c == '{') { ++depth; ++pos; continue; }
        if (c == '}')
        {
            if (--depth < 0) return {};
            ++pos;
            continue;
        }
        if (depth == 0 && (std::isalpha(static_cast<unsigned char>(c)) || c == '_'))
        {
            const size_t start = pos++;
            while (pos < text.size() && (std::isalnum(static_cast<unsigned char>(text[pos]))
                || text[pos] == '_')) ++pos;
            if (text.substr(start, pos - start) != "SingleModule") continue;
            size_t brace = pos;
            while (brace < text.size() && std::isspace(static_cast<unsigned char>(text[brace]))) ++brace;
            // Also preserve the old top-level "app SingleModule" record.
            if (brace == text.size() || text[brace] != '{') continue;
            if (++blocks != 1
                || !parser.extractBlockText(text.substr(start), "SingleModule", body)) return {};
            continue;
        }
        ++pos;
    }
    if (depth != 0 || blocks != 1) return {};

    size_t pos = 0;
    auto skipSpace = [&]()
    {
        while (pos < body.size() && std::isspace(static_cast<unsigned char>(body[pos]))) ++pos;
    };
    skipSpace();
    const size_t keyStart = pos;
    while (pos < body.size() && std::isalpha(static_cast<unsigned char>(body[pos]))) ++pos;
    if (toLowerCopy(body.substr(keyStart, pos - keyStart)) != "modules") return {};
    if (pos < body.size() && !std::isspace(static_cast<unsigned char>(body[pos])) && body[pos] != '(')
        return {};
    skipSpace();
    if (pos == body.size()) return {};

    string list;
    if (body[pos] == '(')
    {
        const size_t start = ++pos;
        const size_t end = body.find(')', start);
        if (end == string::npos) return {};
        list = body.substr(start, end - start);
        pos = end + 1;
    }
    else
    {
        // The historical unparenthesized syntax is one record on a single line.
        const size_t start = pos;
        while (pos < body.size() && body[pos] != '\n' && body[pos] != '\r' && body[pos] != ';') ++pos;
        list = body.substr(start, pos - start);
    }
    skipSpace();
    if (pos < body.size() && body[pos] == ';') { ++pos; skipSpace(); }
    if (pos != body.size()) return {}; // no duplicate key, nested block, or trailing fragment

    vector<string> modules;
    size_t tokenStart = 0;
    while (tokenStart < list.size())
    {
        while (tokenStart < list.size() && (std::isspace(static_cast<unsigned char>(list[tokenStart]))
            || list[tokenStart] == ',')) ++tokenStart;
        if (tokenStart == list.size()) break;
        size_t tokenEnd = tokenStart;
        while (tokenEnd < list.size() && !std::isspace(static_cast<unsigned char>(list[tokenEnd]))
            && list[tokenEnd] != ',') ++tokenEnd;
        const string token = list.substr(tokenStart, tokenEnd - tokenStart);
        if (token.empty() || !(std::isalpha(static_cast<unsigned char>(token.front())) || token.front() == '_'))
            return {};
        for (const unsigned char letter : token)
            if (!std::isalnum(letter) && letter != '_' && letter != '-') return {};
        modules.push_back(token);
        tokenStart = tokenEnd;
    }
    return modules;
}

SingleModuleApp::SingleModuleApp(const RunConfig& cfg) : cfg_(cfg) {}

std::string SingleModuleApp::name() const { return "SingleModule"; }

int SingleModuleApp::run()
{
    // Determine execution list from ./system/run.dict (optional).
    namespace fs = std::filesystem;
    const bool hasRunDict = fs::exists("./system/run.dict");
    vector<string> selected;
    if (hasRunDict)
    {
        selected = readModulesFromRunDict("./system/run.dict");
    }
    else
    {
        selected = {"ExplicitPIC"};
        cout << "[SingleModuleApp] run.dict not found; defaulting to module: ExplicitPIC" << endl;
    }

    if (selected.empty())
    {
        cerr << "[SingleModuleApp] ERROR: cannot find module list in ./system/run.dict.\n"
             << "  Expected:\n"
             << "    Run { app SingleModule }\n"
             << "    SingleModule { modules ( Electricity ) }\n";
        return 2;
    }

    ModuleRegistry registry = makeDefaultModuleRegistry();

    cout << "[SingleModuleApp] Module selection:";
    for (const auto& m : selected)
        cout << " " << m;
    cout << endl;

    // Validate/construct the complete selection before running any module.
    vector<unique_ptr<IModule>> modules;
    for (const auto& modName : selected)
    {
        unique_ptr<IModule> module = registry.create(modName, cfg_);
        if (!module)
        {
            cerr << "[SingleModuleApp] ERROR: Unknown module: '" << modName << "'. Available:";
            for (const auto& k : registry.registeredNames())
                cerr << " " << k;
            cerr << endl;
            return 2;
        }
        modules.push_back(std::move(module));
    }
    for (const auto& module : modules)
    {
        cout << "[SingleModuleApp] Running module: " << module->name() << endl;
        const int rc = module->run();
        if (rc != 0)
        {
            cerr << "[SingleModuleApp] Module '" << module->name() << "' returned rc=" << rc << endl;
            return rc;
        }
    }

    return 0;
}
