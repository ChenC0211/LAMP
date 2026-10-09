#include "Core/IO/DictionaryMesh.h"

#include "Core/Base/Base.h"
#include "Core/Units/UnitParser.h"
#include "Core/Units/UnitsContext.h"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <sstream>
#include <set>
#include <stdexcept>

using namespace std;

namespace
{
string lowerCopy(string s)
{
    transform(s.begin(), s.end(), s.begin(),
              [](unsigned char c) { return static_cast<char>(tolower(c)); });
    return s;
}

string stripSymbols(string s)
{
    auto isStrip = [](char c) {
        return c == '{' || c == '}' || c == '(' || c == ')' || c == ',' || c == ';';
    };
    while (!s.empty() && isStrip(s.front()))
        s.erase(s.begin());
    while (!s.empty() && isStrip(s.back()))
        s.pop_back();
    return s;
}

bool parseBoolToken(const string& s, bool& out)
{
    const string v = lowerCopy(s);
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

bool isMeshBlockName(const string& name)
{
    const string n = lowerCopy(name);
    return (n == "runtimeadaptation" || n == "geometryrefinement" || n == "physicsrefinement" || n == "geometryrefresh" ||
            n == "refineboundary" || n == "refinephase" || n == "refineregion" || n == "caps");
}
} // namespace

DictionaryMesh::DictionaryMesh()
{
    // Provide safe defaults so that a minimal mesh.dict (domain+divide only)
    // can still instantiate a QuadtreeMesh without requiring adaptive blocks.
    maxCellNum = 200000;
    levelSetCount = 1;
    minLevel_neg.assign(levelSetCount, 0);
    minLevel_pos.assign(levelSetCount, 0);
    maxLevel_local.assign(levelSetCount, 0);
}

DictionaryMesh::~DictionaryMesh() {}

void DictionaryMesh::ensureLevelSetCount(int n)
{
    levelSetCount = std::max(1, n);
    if (minLevel_neg.size() < static_cast<size_t>(levelSetCount))
        minLevel_neg.resize(levelSetCount, 0);
    if (minLevel_pos.size() < static_cast<size_t>(levelSetCount))
        minLevel_pos.resize(levelSetCount, 0);
    if (maxLevel_local.size() < static_cast<size_t>(levelSetCount))
        maxLevel_local.resize(levelSetCount, 0);
}

void DictionaryMesh::readDict(const UnitsContext* units, bool acceptsRuntimeAdaptation)
{
    readDict(meshDictPath, geometryDictPath, units, acceptsRuntimeAdaptation);
}

void DictionaryMesh::readDict(const std::string& meshPath, const std::string& geometryPath,
    const UnitsContext* units, bool acceptsRuntimeAdaptation)
{
    ifstream file(meshPath);

    if (!file.is_open())
    {
        throw std::runtime_error(meshPath + ": cannot open mesh dictionary");
    }

    string line;

    bool isInDomain = false;
    bool isAdaptive = false;
    bool foundDomainMin = false;
    bool foundDomainMax = false;
    bool foundDivide = false;

    vector<string> blockStack;
    string pendingBlock;
    string pendingCapIndicator;
    string activeCapIndicator;
    bool pendingRegionEntry = false;
    bool inRegionEntry = false;
    DictionaryMesh::GeometryRefinementConfig::RegionRule regionEntry;
    bool regionHasName = false;
    bool regionHasShape = false;
    bool regionHasXMin = false;
    bool regionHasXMax = false;
    bool regionHasYMin = false;
    bool regionHasYMax = false;
    bool regionHasMaxLevel = false;
    bool sawLegacyMaxLevelLocal = false;
    bool sawLegacyMinLevel = false;
    std::set<std::string> runtimeKeys;
    std::set<std::string> flowPhysicsKeys;
    bool sawFlowPhysicsBlock = false;
    int braceDepth = 0;

    auto pushBlock = [&](const string& name) {
        if (lowerCopy(name) == "runtimeadaptation")
        {
            if (runtimeAdaptation.present || !blockStack.empty())
                throw std::invalid_argument("mesh.dict: RuntimeAdaptation must be one top-level block");
            runtimeAdaptation.present = true;
        }
        if (acceptsRuntimeAdaptation && lowerCopy(name) == "physicsrefinement")
        {
            if (sawFlowPhysicsBlock || !blockStack.empty())
                throw std::invalid_argument("mesh.dict: PhysicsRefinement must be one top-level block");
            sawFlowPhysicsBlock = true;
        }
        blockStack.push_back(lowerCopy(name));
    };
    auto popBlock = [&]() {
        if (!blockStack.empty())
            blockStack.pop_back();
    };
    auto inBlock = [&](const string& name) -> bool {
        const string key = lowerCopy(name);
        return std::find(blockStack.begin(), blockStack.end(), key) != blockStack.end();
    };

    auto resetRegionEntry = [&]() {
        regionEntry = DictionaryMesh::GeometryRefinementConfig::RegionRule{};
        regionHasName = false;
        regionHasShape = false;
        regionHasXMin = false;
        regionHasXMax = false;
        regionHasYMin = false;
        regionHasYMax = false;
        regionHasMaxLevel = false;
    };

    auto regionEntryLabel = [&]() -> std::string {
        if (regionHasName && !regionEntry.name.empty())
            return "'" + regionEntry.name + "'";
        return "<unnamed>";
    };

    auto parseLengthValue = [&](const std::string& text, double& out) -> bool {
        const std::string expr = trim(text);
        if (expr.empty())
            return false;

        UnitParser::Result r;
        if (UnitParser::TryParseToSI(expr, r))
        {
            out = units ? units->LengthSIToCode(r.value_si) : r.value_si;
            return true;
        }

        try
        {
            out = std::stod(expr);
            return true;
        }
        catch (...)
        {
            return false;
        }
    };

    auto failRegionEntry = [&](const std::string& msg) -> void {
        throw std::invalid_argument(meshPath + ": GeometryRefinement.refineRegion entry "
            + regionEntryLabel() + " " + msg);
    };

    auto finalizeRegionEntry = [&]() {
        if (!regionHasName || regionEntry.name.empty())
            failRegionEntry("is missing required field 'name'.");
        if (!regionHasShape)
            failRegionEntry("is missing required field 'shape'.");
        if (regionEntry.shape != DictionaryMesh::GeometryRefinementConfig::RegionShape::kBox)
            failRegionEntry("uses unsupported shape; only 'box' is implemented.");
        if (!regionHasXMin || !regionHasXMax || !regionHasYMin || !regionHasYMax)
            failRegionEntry("is missing one or more required bounds: xmin/xmax/ymin/ymax.");
        if (!regionHasMaxLevel)
            failRegionEntry("is missing required field 'maxLevel'.");
        if (regionEntry.maxLevel < 0)
            failRegionEntry("has invalid maxLevel (must be >= 0).");
        if (regionEntry.min.x > regionEntry.max.x)
            failRegionEntry("has xmin > xmax.");
        if (regionEntry.min.y > regionEntry.max.y)
            failRegionEntry("has ymin > ymax.");

        geometryRefinement.regionRules.push_back(regionEntry);
        resetRegionEntry();
        inRegionEntry = false;
        pendingRegionEntry = false;
    };

    auto parseIndicatorList = [&](const string& lineText, vector<string>& out) {
        const size_t l = lineText.find('(');
        const size_t r = lineText.find(')');
        if (l == string::npos || r == string::npos || r <= l + 1)
            return;
        string inside = lineText.substr(l + 1, r - l - 1);
        for (char& c : inside)
        {
            if (c == ',' || c == ';')
                c = ' ';
        }
        istringstream iss(inside);
        string token;
        while (iss >> token)
        {
            token = stripSymbols(token);
            if (token.empty())
                continue;
            if (find(out.begin(), out.end(), token) == out.end())
                out.push_back(token);
        }
    };

    while (getline(file, line)) 
    {
        line = removeInlineComments(line);
        line = trim(line);

        if (line.empty()) 
            continue;

        if (lowerCopy(pendingBlock) == "runtimeadaptation" && line != "{")
            throw std::invalid_argument("mesh.dict/RuntimeAdaptation: expected opening brace");
        if (acceptsRuntimeAdaptation && lowerCopy(pendingBlock) == "physicsrefinement" && line != "{")
            throw std::invalid_argument("mesh.dict/PhysicsRefinement: expected opening brace");
        const auto opening = line.find('{');
        if (lowerCopy(trim(line.substr(0, opening))) == "runtimeadaptation")
        {
            if (braceDepth != 0 || (opening != std::string::npos && !trim(line.substr(opening + 1)).empty()))
                throw std::invalid_argument("mesh.dict/RuntimeAdaptation: use a top-level multiline block");
        }
        if (acceptsRuntimeAdaptation && lowerCopy(trim(line.substr(0, opening))) == "physicsrefinement"
            && (braceDepth != 0 || (opening != std::string::npos && !trim(line.substr(opening + 1)).empty())))
            throw std::invalid_argument("mesh.dict/PhysicsRefinement: use a top-level multiline block");
        braceDepth += static_cast<int>(std::count(line.begin(), line.end(), '{'))
            - static_cast<int>(std::count(line.begin(), line.end(), '}'));

        if (inBlock("runtimeadaptation") && line != "}")
        {
            auto tokens = split(line);
            if (tokens.size() != 2 || line.find_first_of("{}") != std::string::npos)
                throw std::invalid_argument("mesh.dict/RuntimeAdaptation: use one key and value per line");
            const auto key = tokens[0];
            const auto value = stripSymbols(tokens[1]);
            if (!runtimeKeys.insert(key).second)
                throw std::invalid_argument("mesh.dict/RuntimeAdaptation: duplicate key " + key);
            if (key == "enabled")
            {
                if (value != "true" && value != "false")
                    throw std::invalid_argument("mesh.dict/RuntimeAdaptation/enabled: expected true or false");
                runtimeAdaptation.enabled = value == "true";
            }
            else if (key == "everySteps" || key == "maxSweeps")
            {
                std::size_t consumed = 0;
                const int count = std::stoi(value, &consumed);
                if (consumed != value.size() || count <= 0 || (key == "maxSweeps" && count > 32))
                    throw std::invalid_argument("mesh.dict/RuntimeAdaptation: invalid positive count " + key);
                (key == "everySteps" ? runtimeAdaptation.everySteps : runtimeAdaptation.maxSweeps) = count;
            }
            else throw std::invalid_argument("mesh.dict/RuntimeAdaptation: unsupported key " + key);
            continue;
        }

        if (inRegionEntry)
        {
            if (line == "{")
                continue;

            if (line == "}")
            {
                finalizeRegionEntry();
                continue;
            }

            vector<string> tokens = split(line);
            if (tokens.empty())
                continue;

            const string key = lowerCopy(stripSymbols(tokens[0]));
            if (key == "name")
            {
                const size_t q1 = line.find('"');
                const size_t q2 = (q1 == string::npos) ? string::npos : line.find('"', q1 + 1);
                if (q1 != string::npos && q2 != string::npos && q2 > q1 + 1)
                    regionEntry.name = line.substr(q1 + 1, q2 - q1 - 1);
                else if (tokens.size() >= 2)
                    regionEntry.name = stripSymbols(tokens[1]);
                regionHasName = !regionEntry.name.empty();
                continue;
            }

            if (key == "shape")
            {
                if (tokens.size() < 2)
                    failRegionEntry("is missing value for 'shape'.");
                const string shape = lowerCopy(stripSymbols(tokens[1]));
                if (shape == "box")
                {
                    regionEntry.shape = DictionaryMesh::GeometryRefinementConfig::RegionShape::kBox;
                }
                else if (shape == "circle")
                {
                    regionEntry.shape = DictionaryMesh::GeometryRefinementConfig::RegionShape::kCircle;
                    failRegionEntry("uses unsupported shape 'circle'; only 'box' is implemented.");
                }
                else if (shape == "polygon")
                {
                    regionEntry.shape = DictionaryMesh::GeometryRefinementConfig::RegionShape::kPolygon;
                    failRegionEntry("uses unsupported shape 'polygon'; only 'box' is implemented.");
                }
                else
                {
                    failRegionEntry("uses unsupported shape '" + shape + "'; only 'box' is implemented.");
                }
                regionHasShape = true;
                continue;
            }

            auto parseBoundField = [&](double& dst, bool& seen, const char* fieldName) {
                const size_t keyPos = line.find(tokens[0]);
                const std::string rhs = (keyPos == std::string::npos)
                    ? std::string{}
                    : trim(line.substr(keyPos + tokens[0].size()));
                if (!parseLengthValue(rhs, dst))
                    failRegionEntry(std::string("has invalid value for '") + fieldName + "'.");
                seen = true;
            };

            if (key == "xmin")
            {
                parseBoundField(regionEntry.min.x, regionHasXMin, "xmin");
                continue;
            }
            if (key == "xmax")
            {
                parseBoundField(regionEntry.max.x, regionHasXMax, "xmax");
                continue;
            }
            if (key == "ymin")
            {
                parseBoundField(regionEntry.min.y, regionHasYMin, "ymin");
                continue;
            }
            if (key == "ymax")
            {
                parseBoundField(regionEntry.max.y, regionHasYMax, "ymax");
                continue;
            }
            if (key == "maxlevel")
            {
                if (tokens.size() < 2)
                    failRegionEntry("is missing value for 'maxLevel'.");
                try
                {
                    regionEntry.maxLevel = std::stoi(tokens[1]);
                }
                catch (...)
                {
                    failRegionEntry("has invalid integer value for 'maxLevel'.");
                }
                regionHasMaxLevel = true;
                continue;
            }

            continue;
        }

        if (inBlock("refineregion"))
        {
            const string lower = lowerCopy(line);
            if (pendingRegionEntry && line == "{")
            {
                inRegionEntry = true;
                pendingRegionEntry = false;
                continue;
            }
            if (lower == "entry")
            {
                resetRegionEntry();
                pendingRegionEntry = true;
                continue;
            }
            if (lower.rfind("entry", 0) == 0 && line.find('{') != string::npos)
            {
                resetRegionEntry();
                inRegionEntry = true;
                pendingRegionEntry = false;
                continue;
            }
        }

        if (line == "{")
        {
            if (!pendingBlock.empty())
            {
                pushBlock(pendingBlock);
                pendingBlock.clear();
            }
            else if (!pendingCapIndicator.empty())
            {
                activeCapIndicator = pendingCapIndicator;
                pendingCapIndicator.clear();
            }
            continue;
        }
        if (line == "}")
        {
            if (!activeCapIndicator.empty())
            {
                activeCapIndicator.clear();
                continue;
            }
            popBlock();
            continue;
        }

        {
            const size_t bracePos = line.find('{');
            if (bracePos != string::npos)
            {
                const string before = trim(line.substr(0, bracePos));
                if (!before.empty() && isMeshBlockName(before))
                {
                    const string lower = lowerCopy(before);
                    if (lower == "geometryrefinement")
                        geometryRefinement.present = true;
                    if (lower == "physicsrefinement")
                        physicsRefinement.present = true;
                    if (lower == "geometryrefresh")
                        geometryRefresh.present = true;
                    pushBlock(before);
                    continue;
                }
            }
        }

        if (isMeshBlockName(line))
        {
            pendingBlock = line;
            if (lowerCopy(line) == "geometryrefinement")
                geometryRefinement.present = true;
            if (lowerCopy(line) == "physicsrefinement")
                physicsRefinement.present = true;
            if (lowerCopy(line) == "geometryrefresh")
                geometryRefresh.present = true;
            continue;
        }

        if (line == "domain") 
        {
            isInDomain = true;
            continue;
        } 

        if (isInDomain) 
        {
            if (line.find("min") != string::npos) 
            {
                // Support both legacy (code) and SI (value+unit) formats.
                if (!parseCoordWithUnits(line, minCoord, units))
                    parseCoord(line, minCoord);
                foundDomainMin = true;
            } 
            else if (line.find("max") != std::string::npos) 
            {
                if (!parseCoordWithUnits(line, maxCoord, units))
                    parseCoord(line, maxCoord);
                foundDomainMax = true;
                isInDomain = false;
            }
            continue;
        }

        // -----------------------------------------------------------------
        // Static refinement boxes (optional)
        // Syntax:
        //   refineBox (xmin, ymin) (xmax, ymax) <targetLevel>
        // -----------------------------------------------------------------
        {
            if (line.rfind("refineBox", 0) == 0 || line.rfind("refinebox", 0) == 0)
            {
                vector<Vector2D> coords;
                // Support both legacy and SI formats.
                if (!parseCoordsWithUnits(line, coords, units))
                    parseCoords(line, coords);
                vector<string> tokens = split(line);
                if (coords.size() >= 2 && !tokens.empty())
                {
                    int lvl = 0;
                    // Expect last token to be an integer level.
                    try
                    {
                        lvl = stoi(tokens.back());
                    }
                    catch (...)
                    {
                        lvl = 0;
                    }

                    RefineBox rb;
                    rb.min = coords[0];
                    rb.max = coords[1];
                    rb.level = lvl;
                    refineBoxes.push_back(rb);
                }
                continue;
            }
        }

        if (line.find("divide") != string::npos) 
        {
            Vector2D divide;
            parseCoord(line, divide);
            Nx = static_cast<int>(divide.x);
            Ny = static_cast<int>(divide.y);
            if (Nx <= 0 || Ny <= 0) 
            {
                throw std::invalid_argument(meshPath + ": invalid grid division (Nx/Ny must be > 0).");
            }
            foundDivide = true;
            continue;
        } 

        if (line.find("adaptiveMesh") != string::npos)
        {
            vector<string> tokens = split(line);
            if (tokens.size() == 2) 
            {
                isAdaptive = (tokens[1] == "true" || tokens[1] == "True");
                if (isAdaptive)
                {
                    // Switch to explicit adaptive settings from file.
                    minLevel_neg.clear();
                    minLevel_pos.clear();
                    maxLevel_local.clear();
                }
            }
            continue;
        }

        if (line.find("maxLevel_local") != string::npos)
        {
            sawLegacyMaxLevelLocal = true;
            continue;
        }

        // -----------------------------------------------------------------
        // V3 GeometryRefinement / PhysicsRefinement parsing
        // -----------------------------------------------------------------
        if (inBlock("refineboundary"))
        {
            string normLine = line;
            for (char& c : normLine)
            {
                if (c == '{' || c == '}' || c == '(' || c == ')' || c == ',' || c == ';')
                    c = ' ';
            }
            vector<string> tokens = split(normLine);
            if (!tokens.empty())
            {
                DictionaryMesh::GeometryRefinementConfig::BoundaryRule rule;
                bool hasBoundary = false;
                bool hasLevel = false;
                for (size_t i = 0; i < tokens.size(); ++i)
                {
                    const string t = lowerCopy(stripSymbols(tokens[i]));
                    if (t == "entry")
                        continue;
                    if ((t == "boundaryid" || t == "boundary" || t == "id") && i + 1 < tokens.size())
                    {
                        rule.boundaryId = stoi(tokens[i + 1]);
                        hasBoundary = true;
                    }
                    else if ((t == "c_gamma" || t == "cgamma" || t == "c-gamma" || t == "cgamma") && i + 1 < tokens.size())
                    {
                        rule.cGamma = stod(tokens[i + 1]);
                    }
                    else if ((t == "maxlevel" || t == "max") && i + 1 < tokens.size())
                    {
                        rule.maxLevel = stoi(tokens[i + 1]);
                        hasLevel = true;
                    }
                }
                if (hasBoundary)
                {
                    if (!hasLevel)
                        rule.maxLevel = -1;
                    if (rule.cGamma <= 0.0)
                        rule.cGamma = 1.0;
                    geometryRefinement.boundaryRules.push_back(rule);
                }
            }
            continue;
        }

        if (inBlock("refinephase"))
        {
            string normLine = line;
            for (char& c : normLine)
            {
                if (c == '{' || c == '}' || c == '(' || c == ')' || c == ',' || c == ';')
                    c = ' ';
            }
            vector<string> tokens = split(normLine);
            if (!tokens.empty())
            {
                DictionaryMesh::GeometryRefinementConfig::PhaseRule rule;
                bool hasPhase = false;
                for (size_t i = 0; i < tokens.size(); ++i)
                {
                    const string t = lowerCopy(stripSymbols(tokens[i]));
                    if (t == "phase" && i + 1 < tokens.size())
                    {
                        rule.phase = stoi(tokens[i + 1]);
                        hasPhase = true;
                    }
                    else if (t == "minlevel" && i + 1 < tokens.size())
                    {
                        rule.minLevel = stoi(tokens[i + 1]);
                    }
                    else if (t == "maxlevel" && i + 1 < tokens.size())
                    {
                        rule.maxLevel = stoi(tokens[i + 1]);
                    }
                }
                if (hasPhase)
                    geometryRefinement.phaseRules.push_back(rule);
            }
            continue;
        }

        if (inBlock("refineregion"))
        {
            continue;
        }

        if (inBlock("geometryrefinement"))
        {
            vector<string> tokens = split(line);
            if (tokens.size() >= 2)
            {
                const string k = lowerCopy(tokens[0]);
                if (k == "enabled")
                    parseBoolToken(tokens[1], geometryRefinement.enabled);
            }
            continue;
        }

        if (inBlock("geometryrefresh"))
        {
            vector<string> tokens = split(line);
            if (tokens.size() >= 2)
            {
                const string k = lowerCopy(tokens[0]);
                geometryRefresh.present = true;
                if (k == "mode")
                    geometryRefresh.mode = stripSymbols(tokens[1]);
                else if (k == "changerel" || k == "relchange" || k == "relth")
                    geometryRefresh.changeRel = stod(tokens[1]);
                else if (k == "verbose")
                    parseBoolToken(tokens[1], geometryRefresh.verbose);
                else if (k == "allowwithmovingeb" || k == "allowwithmovingebrefresh")
                    parseBoolToken(tokens[1], geometryRefresh.allowWithMovingEb);
                else if (k == "movingebrefreshmode" || k == "movingebmode")
                    geometryRefresh.movingEbRefreshMode = stripSymbols(tokens[1]);
            }
            continue;
        }

        if (inBlock("caps"))
        {
            if (acceptsRuntimeAdaptation && inBlock("physicsrefinement"))
            {
                const auto fields = split(line);
                if (fields.size() != 3 || fields[1] != "maxLevel" || !activeCapIndicator.empty())
                    throw std::invalid_argument("mesh.dict/PhysicsRefinement/caps: use Indicator maxLevel N on one line");
                std::size_t consumed = 0;
                const int level = std::stoi(fields[2], &consumed);
                if (level < 0 || consumed != fields[2].size()
                    || !physicsRefinement.maxLevelCaps.emplace(lowerCopy(fields[0]), level).second)
                    throw std::invalid_argument("mesh.dict/PhysicsRefinement/caps: invalid or duplicate cap");
                continue;
            }
            string normLine = line;
            for (char& c : normLine)
            {
                if (c == '{' || c == '}' || c == '(' || c == ')' || c == ',' || c == ';')
                    c = ' ';
            }
            vector<string> tokens = split(normLine);
            if (tokens.empty())
                continue;

            string indicator = activeCapIndicator;
            size_t startIdx = 0;
            if (indicator.empty())
            {
                const string t0 = stripSymbols(tokens[0]);
                if (!t0.empty() && lowerCopy(t0) != "maxlevel")
                {
                    indicator = t0;
                    startIdx = 1;
                }
                else
                {
                    startIdx = 0;
                }
            }

            bool sawMax = false;
            for (size_t i = startIdx; i < tokens.size(); ++i)
            {
                const string t = lowerCopy(stripSymbols(tokens[i]));
                if (t == "maxlevel" && i + 1 < tokens.size())
                {
                    int lvl = stoi(tokens[i + 1]);
                    // Support both per-indicator caps (e.g. `DebyeLength { maxLevel 2 }`)
                    // and a global cap (e.g. `caps { maxLevel 2 }`).
                    if (indicator.empty())
                        indicator = "global";
                    physicsRefinement.maxLevelCaps[lowerCopy(indicator)] = lvl;
                    sawMax = true;
                    break;
                }
            }

            if (!sawMax && indicator.empty() && tokens.size() == 1)
            {
                pendingCapIndicator = stripSymbols(tokens[0]);
            }
            continue;
        }

        if (inBlock("physicsrefinement"))
        {
            vector<string> tokens = split(line);
            if (!tokens.empty())
            {
                const string k = lowerCopy(tokens[0]);
                if (acceptsRuntimeAdaptation)
                {
                    if (!flowPhysicsKeys.insert(k).second)
                        throw std::invalid_argument("mesh.dict/PhysicsRefinement: duplicate key " + tokens[0]);
                    if (k != "enabled" && k != "indicators" && k != "hysteresis"
                        && k != "velocitythreshold" && k != "velocitygradientthreshold")
                        throw std::invalid_argument("mesh.dict/PhysicsRefinement: unsupported flow control " + tokens[0]);
                    if (k == "enabled" && (tokens.size() != 2 || (tokens[1] != "true" && tokens[1] != "false")))
                        throw std::invalid_argument("mesh.dict/PhysicsRefinement/enabled: expected true or false");
                    if (k == "hysteresis")
                    {
                        std::size_t consumed = 0;
                        if (tokens.size() != 2) throw std::invalid_argument("mesh.dict/PhysicsRefinement/hysteresis: expected one number");
                        const double value = std::stod(tokens[1], &consumed);
                        if (consumed != tokens[1].size() || !std::isfinite(value) || value <= 0 || value >= 1)
                            throw std::invalid_argument("mesh.dict/PhysicsRefinement/hysteresis: expected 0 < value < 1");
                    }
                    if (k == "indicators")
                    {
                        const auto begin = line.find('('), end = line.find(')');
                        if (begin == std::string::npos || end == std::string::npos || end <= begin + 1
                            || !trim(line.substr(end + 1)).empty() || trim(line.substr(0, begin)) != "indicators")
                            throw std::invalid_argument("mesh.dict/PhysicsRefinement/indicators: expected a nonempty parenthesized list");
                        auto names = split(line.substr(begin + 1, end - begin - 1));
                        std::set<std::string> unique(names.begin(), names.end());
                        if (names.empty() || names.size() != unique.size())
                            throw std::invalid_argument("mesh.dict/PhysicsRefinement/indicators: empty or duplicate indicator");
                    }
                }
                if (k == "velocitythreshold" || k == "velocitygradientthreshold")
                {
                    UnitParser::Result quantity;
                    const auto dimension = k == "velocitythreshold" ? UnitParser::Dimension::Velocity : UnitParser::Dimension::InverseTime;
                    if (!acceptsRuntimeAdaptation || !UnitParser::TryParseQuantityToSI(trim(line.substr(tokens[0].size())), dimension, quantity)
                        || !(quantity.value_si > 0))
                        throw std::invalid_argument("mesh.dict/PhysicsRefinement: positive, dimensioned Hydro threshold requires a runtime consumer: " + tokens[0]);
                    (k == "velocitythreshold" ? physicsRefinement.velocityThresholdSI : physicsRefinement.velocityGradientThresholdSI) = quantity.value_si;
                    continue;
                }
                if (k == "enabled" && tokens.size() >= 2)
                    parseBoolToken(tokens[1], physicsRefinement.enabled);
                else if (k == "indicators")
                    parseIndicatorList(line, physicsRefinement.indicators);
                else if (k == "adapteverysteps" && tokens.size() >= 2)
                    physicsRefinement.adaptEverySteps = stoi(tokens[1]);
                else if ((k == "indicatoreverysteps" || k == "indicatorevery" || k == "indicatoreverystep") && tokens.size() >= 2)
                    physicsRefinement.indicatorEverySteps = stoi(tokens[1]);
                else if ((k == "reportevery" || k == "reporteverysteps" || k == "amrreportevery") &&
                         tokens.size() >= 2)
                    physicsRefinement.reportEvery = stoi(tokens[1]);
                else if (k == "hysteresis" && tokens.size() >= 2)
                    physicsRefinement.hysteresis = stod(tokens[1]);
                else if (k == "persistencesteps" && tokens.size() >= 2)
                    physicsRefinement.persistenceSteps = stoi(tokens[1]);
            }
            continue;
        }

        if (isAdaptive) 
        {
            if (line.find("maxCellNum") != string::npos) 
            {
                vector<string> tokens = split(line);
                if (tokens.size() == 2) 
                {
                    maxCellNum = stoi(tokens[1]);
                }
                continue;
            }

            if (line.find("maxInitialIterations") != string::npos) 
            {
                vector<string> tokens = split(line);
                if (tokens.size() == 2) 
                {
                    maxInitialIterations = std::stoi(tokens[1]);
                }
                continue;
            }

            if (line.find("C_Lipschitz") != string::npos) 
            {
                vector<string> tokens = split(line);
                if (tokens.size() == 2) 
                {
                    C_Lipschitz = stod(tokens[1]);
                }
                continue;
            }

            if (line.front() == '"' && line.back() == '"') 
            {
                string name = line.substr(1, line.size() - 2); // 去除引号

                while (getline(file, line)) 
                {
                    line = removeInlineComments(line);
                    line = trim(line);

                    if (line.empty()) 
                        continue;

                    if (line.find("minLevel_neg") != string::npos) 
                    {
                        vector<string> tokens = split(line);
                        if (tokens.size() == 2) 
                        {
                            sawLegacyMinLevel = true;
                        }
                    }
                    if (line.find("minLevel_pos") != string::npos) 
                    {
                        vector<string> tokens = split(line);
                        if (tokens.size() == 2) 
                        {
                            sawLegacyMinLevel = true;
                        }
                    }
                    if (line.find("maxLevel_local") != string::npos) 
                    {
                        vector<string> tokens = split(line);
                        if (tokens.size() == 2) 
                        {
                            sawLegacyMaxLevelLocal = true;
                        }
                    }
                }
                continue;
            }
        }
    }

    file.close();

    if (lowerCopy(pendingBlock) == "runtimeadaptation")
        throw std::invalid_argument("mesh.dict/RuntimeAdaptation: missing block body");

    if (runtimeAdaptation.present)
    {
        if (inBlock("runtimeadaptation") || !runtimeKeys.count("enabled"))
            throw std::invalid_argument("mesh.dict/RuntimeAdaptation: incomplete block or missing enabled");
        if (!runtimeAdaptation.enabled && runtimeKeys.size() != 1)
            throw std::invalid_argument("mesh.dict/RuntimeAdaptation: disabled block must not carry inactive controls");
        if (runtimeAdaptation.enabled && (!acceptsRuntimeAdaptation || !isAdaptive
            || !geometryRefinement.enabled))
            throw std::invalid_argument("mesh.dict/RuntimeAdaptation requires an explicit workflow consumer, adaptiveMesh and GeometryRefinement");
    }

    if (acceptsRuntimeAdaptation && physicsRefinement.present)
    {
        if (!sawFlowPhysicsBlock || inBlock("physicsrefinement") || lowerCopy(pendingBlock) == "physicsrefinement")
            throw std::invalid_argument("mesh.dict/PhysicsRefinement: missing or incomplete block");
        if (!flowPhysicsKeys.count("enabled"))
            throw std::invalid_argument("mesh.dict/PhysicsRefinement: missing enabled");
        if (!physicsRefinement.enabled && (flowPhysicsKeys.size() != 1 || !physicsRefinement.maxLevelCaps.empty()))
            throw std::invalid_argument("mesh.dict/PhysicsRefinement: disabled block carries inactive controls");
        if (physicsRefinement.enabled && (!runtimeAdaptation.enabled || !flowPhysicsKeys.count("hysteresis")))
            throw std::invalid_argument("mesh.dict/PhysicsRefinement: flow indicators require RuntimeAdaptation and explicit hysteresis");
    }

    if (pendingRegionEntry || inRegionEntry)
    {
        throw std::invalid_argument(meshPath + ": unterminated GeometryRefinement.refineRegion entry.");
    }

    if (sawLegacyMinLevel)
    {
        std::cerr << "[DictionaryMesh] WARNING: legacy minLevel_* settings are deprecated and ignored under v3 mesh.dict.\n";
    }
    if (sawLegacyMaxLevelLocal)
    {
        std::cerr << "[DictionaryMesh] WARNING: maxLevel_local is deprecated in v3 mesh.dict; "
                     "using requiredMax from refinement rules instead.\n";
    }

    // For adaptive meshes, derive maxLevel_local from v3 refinement caps/rules.
    if (isAdaptive)
    {
        int requiredMax = 0;
        if (geometryRefinement.enabled)
        {
            for (const auto& r : geometryRefinement.boundaryRules)
            {
                if (r.maxLevel > 0)
                    requiredMax = std::max(requiredMax, r.maxLevel);
            }
            for (const auto& r : geometryRefinement.phaseRules)
            {
                if (r.maxLevel > 0)
                    requiredMax = std::max(requiredMax, r.maxLevel);
                if (r.minLevel > 0)
                    requiredMax = std::max(requiredMax, r.minLevel);
            }
            for (const auto& r : geometryRefinement.regionRules)
            {
                if (r.maxLevel > 0)
                    requiredMax = std::max(requiredMax, r.maxLevel);
            }
        }
        if (physicsRefinement.enabled)
        {
            for (const auto& kv : physicsRefinement.maxLevelCaps)
                requiredMax = std::max(requiredMax, kv.second);
        }

        if (requiredMax > 0)
        {
            maxLevel_local.assign(std::max(1, levelSetCount), requiredMax);
        }
        else if (geometryRefinement.enabled || physicsRefinement.enabled)
        {
            std::cerr << "[DictionaryMesh] WARNING: adaptiveMesh enabled but maxLevel_local is missing; "
                         "no refinement caps detected. AMR may be disabled." << std::endl;
        }
    }

    // -----------------------------------------------------------------
    // Domain ownership (v3): the computational domain (box) is defined in
    // geometry.dict. mesh.dict may still contain a legacy 'domain' block,
    // but it is deprecated and will be ignored when geometry.dict provides
    // Domain/min/max.
    // -----------------------------------------------------------------
    bool geomMin = false;
    bool geomMax = false;
    bool inDomain = false;

    std::ifstream geomFile(geometryPath);
    if (!geomFile.is_open())
    {
        throw std::runtime_error(geometryPath + ": cannot open geometry dictionary (Domain must be defined there).");
    }

    while (std::getline(geomFile, line))
    {
        line = removeInlineComments(line);
        line = trim(line);
        if (line.empty())
            continue;

        const std::string lower = lowerCopy(line);
        if (lower == "domain")
        {
            inDomain = true;
            continue;
        }

        if (inDomain)
        {
            if (line.find("min") != std::string::npos)
            {
                if (!parseCoordWithUnits(line, minCoord, units))
                    parseCoord(line, minCoord);
                geomMin = true;
            }
            else if (line.find("max") != std::string::npos)
            {
                if (!parseCoordWithUnits(line, maxCoord, units))
                    parseCoord(line, maxCoord);
                geomMax = true;
                inDomain = false;
            }
        }
    }
    geomFile.close();

    if (!geomMin || !geomMax)
    {
        throw std::invalid_argument(geometryPath + ": missing Domain/min/max.");
    }
    if (foundDomainMin || foundDomainMax)
    {
        std::cerr << "[DictionaryMesh] WARNING: 'domain' in mesh.dict is deprecated; using geometry.dict Domain instead.\n";
    }
    foundDomainMin = true;
    foundDomainMax = true;
    if (maxCoord.x <= minCoord.x || maxCoord.y <= minCoord.y)
    {
        throw std::invalid_argument(geometryPath + ": invalid domain bounds (max must be greater than min).");
    }
    if (!foundDivide)
    {
        throw std::invalid_argument(meshPath + ": missing divide (Nx, Ny).");
    }

    // Safety: ensure adaptive level arrays are at least levelSetCount long
    // to avoid out-of-bounds access in QuadtreeMesh construction.
    ensureLevelSetCount(levelSetCount);
}
