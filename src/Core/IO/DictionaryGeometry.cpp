#include "Core/IO/DictionaryGeometry.h"

#include <algorithm>
#include <cassert>
#include <limits>
#include <sstream>
#include <stdexcept>

#include "Core/Units/UnitParser.h"

using namespace std;

namespace
{

// Braces in quoted names/paths do not change section nesting.
int braceDelta(const std::string& line)
{
    int delta = 0;
    char quote = 0;
    bool escaped = false;
    for (char c : line)
    {
        if (escaped) { escaped = false; continue; }
        if (quote && c == '\\') { escaped = true; continue; }
        if (quote) { if (c == quote) quote = 0; continue; }
        if (c == '"' || c == '\'') { quote = c; continue; }
        if (c == '{') ++delta;
        if (c == '}') --delta;
    }
    return delta;
}
string toLowerCopy(string s)
{
    transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return static_cast<char>(tolower(c)); });
    return s;
}

bool parseBoolToken(const string& s, bool& out)
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

int toZeroBased(int v)
{
    if (v > 0)
        return v - 1;
    return v;
}

void adjustIndex(SignInt& s)
{
    s.value = toZeroBased(s.value);
}

bool extractQuotedName(const string& line, string& out)
{
    const size_t q1 = line.find('"');
    if (q1 == string::npos)
        return false;
    const size_t q2 = line.find('"', q1 + 1);
    if (q2 == string::npos || q2 <= q1 + 1)
        return false;
    out = line.substr(q1 + 1, q2 - q1 - 1);
    return true;
}

bool parseConstraint(const string& line, string& lsName, bool& greaterThanZero)
{
    string name;
    if (!extractQuotedName(line, name))
        return false;

    const size_t gt = line.find('>');
    const size_t lt = line.find('<');
    if (gt == string::npos && lt == string::npos)
        return false;

    greaterThanZero = (gt != string::npos && (lt == string::npos || gt < lt));
    lsName = name;
    return true;
}

bool inferPhaseLevelSetInsideIsNegative(const DictionaryGeometry& dictGeometry,
                                        const int phaseIndex,
                                        const int levelSetIndex,
                                        bool& insideIsNegative)
{
    if (phaseIndex < 0 || phaseIndex >= static_cast<int>(dictGeometry.levelSetIndexA.size()))
        return false;
    if (levelSetIndex < 0 || levelSetIndex >= static_cast<int>(dictGeometry.levelSetName.size()))
        return false;

    bool found = false;
    bool candidate = true;
    auto collect = [&](const vector<SignInt>& refs) -> bool {
        for (const auto& ref : refs)
        {
            if (ref.value != levelSetIndex)
                continue;
            const bool localInsideIsNegative = !ref.isNegative;
            if (!found)
            {
                candidate = localInsideIsNegative;
                found = true;
                continue;
            }
            if (candidate != localInsideIsNegative)
                return false;
        }
        return true;
    };

    if (!collect(dictGeometry.levelSetIndexA[static_cast<size_t>(phaseIndex)]))
        return false;
    if (phaseIndex < static_cast<int>(dictGeometry.levelSetIndexB.size()) &&
        !collect(dictGeometry.levelSetIndexB[static_cast<size_t>(phaseIndex)]))
        return false;

    if (!found)
        return false;

    insideIsNegative = candidate;
    return true;
}
}

std::vector<int> DictionaryGeometry::immersedBoundaryIds() const
{
    std::vector<int> out;
    out.reserve(boundaries.size());
    for (const auto& ent : boundaries)
    {
        if (ent.id <= 0 || ent.ignored || !ent.immersed || ent.kind == GeometricBoundary::domain)
            continue;
        out.push_back(ent.id);
    }
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
    return out;
}

bool DictionaryGeometry::inferLevelSetInsideIsNegative(const int levelSetIndex,
                                                       bool& insideIsNegative) const
{
    bool found = false;
    bool candidate = true;
    for (const auto& ent : boundaries)
    {
        if (ent.ignored || !ent.immersed || ent.kind == GeometricBoundary::domain)
            continue;
        if (ent.levelSetIndex != levelSetIndex)
            continue;

        bool localInsideIsNegative = true;
        bool localFound = false;
        if (ent.phase >= 0 &&
            inferPhaseLevelSetInsideIsNegative(*this, ent.phase, ent.levelSetIndex, localInsideIsNegative))
        {
            localFound = true;
        }
        else if (ent.opposedPhase >= 0)
        {
            bool opposedInsideIsNegative = true;
            if (inferPhaseLevelSetInsideIsNegative(*this, ent.opposedPhase, ent.levelSetIndex, opposedInsideIsNegative))
            {
                localInsideIsNegative = !opposedInsideIsNegative;
                localFound = true;
            }
        }
        if (!localFound)
            continue;

        if (!found)
        {
            candidate = localInsideIsNegative;
            found = true;
            continue;
        }
        if (candidate != localInsideIsNegative)
            return false;
    }

    if (!found)
    {
        for (int phaseIndex = 0; phaseIndex < static_cast<int>(phaseName.size()); ++phaseIndex)
        {
            bool localInsideIsNegative = true;
            if (!inferPhaseLevelSetInsideIsNegative(*this, phaseIndex, levelSetIndex, localInsideIsNegative))
                continue;

            if (!found)
            {
                candidate = localInsideIsNegative;
                found = true;
                continue;
            }
            if (candidate != localInsideIsNegative)
                return false;
        }
    }

    if (!found)
        return false;

    insideIsNegative = candidate;
    return true;
}

bool DictionaryGeometry::inferBoundaryInsideIsNegative(const int boundaryId,
                                                       bool& insideIsNegative) const
{
    const auto it = boundaryById.find(boundaryId);
    if (it == boundaryById.end())
        return false;

    const auto& ent = it->second;
    if (ent.levelSetIndex < 0 || ent.levelSetIndex >= static_cast<int>(levelSetName.size()))
        return false;

    if (ent.phase >= 0 &&
        inferPhaseLevelSetInsideIsNegative(*this, ent.phase, ent.levelSetIndex, insideIsNegative))
    {
        return true;
    }

    if (ent.opposedPhase >= 0)
    {
        bool opposedInsideIsNegative = true;
        if (inferPhaseLevelSetInsideIsNegative(*this, ent.opposedPhase, ent.levelSetIndex, opposedInsideIsNegative))
        {
            insideIsNegative = !opposedInsideIsNegative;
            return true;
        }
    }

    return inferLevelSetInsideIsNegative(ent.levelSetIndex, insideIsNegative);
}

bool DictionaryGeometry::inferEmbeddedInsideIsNegative(bool& insideIsNegative) const
{
    bool found = false;
    bool candidate = true;
    for (const int boundaryId : immersedBoundaryIds())
    {
        bool localInsideIsNegative = true;
        if (!inferBoundaryInsideIsNegative(boundaryId, localInsideIsNegative))
            continue;

        if (!found)
        {
            candidate = localInsideIsNegative;
            found = true;
            continue;
        }
        if (candidate != localInsideIsNegative)
            return false;
    }

    if (!found)
        return false;

    insideIsNegative = candidate;
    return true;
}


void DictionaryGeometry::captureExtension(std::istream& file, const std::string& name)
{
    const auto begin = file.tellg();
    std::string line;
    while (std::getline(file, line))
    {
        line = trim(removeInlineComments(line));
        if (!line.empty()) break;
    }
    if (line != "{")
    {
        file.clear();
        file.seekg(begin);
        return;
    }
    std::string block = name + "\n{\n";
    int depth = 1;
    while (depth > 0 && std::getline(file, line))
    {
        block += line + "\n";
        depth += braceDelta(removeInlineComments(line));
    }
    extensionBlocks_[name] += block;
}

bool DictionaryGeometry::extensionEnabled(const std::string& name, bool defaultWhenPresent) const
{
    const auto* text = extensionBlock(name);
    if (!text) return false;
    std::istringstream input(*text);
    bool enabled = defaultWhenPresent;
    int depth = 0;
    std::string line;
    while (std::getline(input, line))
    {
        line = trim(removeInlineComments(line));
        if (depth == 0 && line == name) enabled = defaultWhenPresent;
        if (depth == 1)
        {
            const auto tokens = split(line);
            if (tokens.size() >= 2 && toLowerCopy(tokens[0]) == "enabled")
                parseBoolToken(tokens[1], enabled);
        }
        depth += braceDelta(line);
    }
    return enabled;
}

DictionaryGeometry::DictionaryGeometry() {}

DictionaryGeometry::~DictionaryGeometry() {}

void DictionaryGeometry::readDict(const UnitsContext* units)
{
    readDict(geometryDictPath, units);
}

void DictionaryGeometry::readDict(const std::string& path, const UnitsContext* units)
{
    ifstream file(path);

    if (!file.is_open())
    {
        throw std::runtime_error(path + ": cannot open geometry dictionary");
    }

    points.clear();
    lines.clear();
    levelSetName.clear();
    innerPoints.clear();
    moving.clear();
    linesOnInterface.clear();
    signClosureLines.clear();
    phaseName.clear();
    levelSetOperation.clear();
    levelSetIndexA.clear();
    levelSetIndexB.clear();
    boundaryLines.clear();
    boundaries.clear();
    boundaryById.clear();
    periodicPairs.clear();
    periodicMateByBoundary.clear();
    periodicShiftToMate.clear();
    hasDomain = false;
    hasDimension = false;
    hasSlabThickness = false;
    dimension = 2;
    slabThickness_m = 1.0;
    extensionBlocks_.clear();

    string line;
    while (getline(file, line))
    {
        line = removeInlineComments(line);
        line = trim(line);

        if (line.empty() || line == "{")
            continue;

        if (line == "Domain")
        {
            while (getline(file, line))
            {
                line = removeInlineComments(line);
                line = trim(line);

                if (line.empty() || line == "{")
                    continue;
                if (line == "}")
                    break;

                if (line.find("min") != string::npos)
                {
                    if (!parseCoordWithUnits(line, domainMin, units))
                        parseCoord(line, domainMin);
                    hasDomain = true;
                }
                else if (line.find("max") != string::npos)
                {
                    if (!parseCoordWithUnits(line, domainMax, units))
                        parseCoord(line, domainMax);
                    hasDomain = true;
                }
                else
                {
                    // Optional domain metadata.
                    // Supported keys: dimension, slabThickness.
                    auto tokens = split(line);
                    if (tokens.size() >= 2)
                    {
                        const std::string key = toLowerCopy(tokens[0]);
                        if (key == "dimension" || key == "dim")
                        {
                            try
                            {
                                dimension = std::stoi(tokens[1]);
                                hasDimension = true;
                            }
                            catch (...) {}
                        }
                        else if (key == "slabthickness" || key == "thickness" || key == "h0")
                        {
                            // Parse as length to SI meters. Accept both "<val> <unit>" and "<val><unit>".
                            std::string joined;
                            for (size_t i = 1; i < tokens.size(); ++i)
                            {
                                if (!joined.empty())
                                    joined += " ";
                                joined += tokens[i];
                            }
                            UnitParser::Result r;
                            if (UnitParser::TryParseToSI(joined, r))
                            {
                                slabThickness_m = r.value_si;
                                hasSlabThickness = true;
                            }
                            else
                            {
                                // Fallback: treat as meters.
                                try
                                {
                                    slabThickness_m = std::stod(tokens[1]);
                                    hasSlabThickness = true;
                                }
                                catch (...) {}
                            }
                        }
                    }
                }
            }
            continue;
        }

        if (line == "PeriodicPairs")
        {
            // Format:
            // PeriodicPairs
            // {
            //     pair <id>
            //     {
            //         boundaryA <int>
            //         boundaryB <int>
            //         shift (dx, dy) [unit]
            //         mode strict
            //     }
            // }
            while (getline(file, line))
            {
                line = removeInlineComments(line);
                line = trim(line);

                if (line.empty() || line == "{")
                    continue;

                if (line == "}")
                    break;

                vector<string> tokens = split(line);
                if (tokens.empty())
                    continue;

                const string key0 = toLowerCopy(tokens[0]);
                if (key0 == "pair")
                {
                    PeriodicPair pp;

                    // Enter pair block
                    while (getline(file, line))
                    {
                        line = removeInlineComments(line);
                        line = trim(line);

                        if (line.empty() || line == "{")
                            continue;
                        if (line == "}")
                            break;

                        vector<string> t = split(line);
                        if (t.empty())
                            continue;

                        const string k = toLowerCopy(t[0]);
                        if (k == "boundarya" && t.size() >= 2)
                        {
                            pp.boundaryA = stoi(t[1]);
                        }
                        else if (k == "boundaryb" && t.size() >= 2)
                        {
                            pp.boundaryB = stoi(t[1]);
                        }
                        else if (k == "shift")
                        {
                            Vector2D sh;
                            if (parseCoordWithUnits(line, sh, units) || parseCoord(line, sh))
                                pp.shiftAtoB = sh;
                        }
                        else if (k == "mode" && t.size() >= 2)
                        {
                            pp.mode = toLowerCopy(t[1]);
                        }
                    }

                    if (pp.boundaryA <= 0 || pp.boundaryB <= 0)
                    {
                        cerr << "[DictionaryGeometry] WARNING: PeriodicPairs has invalid boundary ids; skipping pair." << endl;
                        continue;
                    }
                    if (pp.boundaryA == pp.boundaryB)
                    {
                        cerr << "[DictionaryGeometry] WARNING: PeriodicPairs uses same boundary for A and B (" << pp.boundaryA
                             << "); skipping pair." << endl;
                        continue;
                    }
                    periodicPairs.push_back(pp);
                }
            }
            continue;
        }

        if (line == "Point")
        {
            while (getline(file, line))
            {
                line = removeInlineComments(line);
                line = trim(line);

                if (line.empty() || line == "{")
                    continue;

                if (line == "}")
                    break;

                Vector2D coord;
                if (parseCoordWithUnits(line, coord, units) || parseCoord(line, coord))
                    points.push_back(coord);
            }
            continue;
        }

        if (line == "Line")
        {
            while (getline(file, line))
            {
                line = removeInlineComments(line);
                line = trim(line);

                if (line.empty() || line == "{")
                    continue;

                if (line == "}")
                    break;

                vector<SignInt> pointIndex;
                if (parseIndex(line, pointIndex))
                {
                    if (pointIndex.size() < 2)
                        continue;

                    adjustIndex(pointIndex[0]);
                    adjustIndex(pointIndex[1]);

                    const int idx0 = pointIndex[0].value;
                    const int idx1 = pointIndex[1].value;
                    if (idx0 < 0 || idx0 >= static_cast<int>(points.size()) ||
                        idx1 < 0 || idx1 >= static_cast<int>(points.size()))
                    {
                        cerr << "[DictionaryGeometry] WARNING: line uses invalid point index: "
                             << pointIndex[0].value + 1 << ", " << pointIndex[1].value + 1 << endl;
                        continue;
                    }

                    Line l;
                    l.points[0] = points[idx0];
                    l.points[1] = points[idx1];
                    lines.push_back(l);

                    while (getline(file, line))
                    {
                        line = removeInlineComments(line);
                        line = trim(line);

                        if (line.empty() || line == "{")
                            continue;

                        if (line == "}")
                            break;

                        if (line.find("type") != string::npos)
                        {
                            vector<string> tokens = split(line);
                            if (tokens.size() >= 2)
                            {
                                if (tokens[1] == "straight")
                                    lines.back().setType(LineType::straight);
                                else if (tokens[1] == "arc")
                                    lines.back().setType(LineType::arc);
                                else if (tokens[1] == "spline")
                                    lines.back().setType(LineType::spline);
                            }
                        }
                        else if (line.find("center") != string::npos)
                        {
                            Vector2D center;
                            if (parseCoordWithUnits(line, center, units) || parseCoord(line, center))
                                lines.back().setArcCenter(center);
                        }
                        else if (line.find("radius") != string::npos)
                        {
                            vector<string> tokens = split(line);
                            if (tokens.size() >= 2)
                                lines.back().setRadius(stod(tokens[1]));
                        }
                        else if (line.find("radian") != string::npos)
                        {
                            vector<string> tokens = split(line);
                            if (tokens.size() >= 2)
                                lines.back().setRadian(stod(tokens[1]));
                        }
                        else if (line.find("featurePoints") != string::npos)
                        {
                            vector<Vector2D> featurePoints;
                            if (parseCoordsWithUnits(line, featurePoints, units) || parseCoords(line, featurePoints))
                            {
                                for (Vector2D& point : featurePoints)
                                    lines.back().addFeaturePoint(point);
                            }
                        }
                    }
                }
            }
            continue;
        }

        if (line == "LevelSet")
        {
            while (getline(file, line))
            {
                line = removeInlineComments(line);
                line = trim(line);

                if (line.empty() || line == "{")
                    continue;

                if (line == "}")
                    break;

                string lsName;
                if (extractQuotedName(line, lsName))
                {
                    vector<SignInt> lsLines;
                    vector<SignInt> signClosure;
                    bool haveSignClosure = false;
                    Vector2D innerPoint(0.0, 0.0);
                    bool haveInner = false;
                    bool isMoving = false;

                    while (getline(file, line))
                    {
                        line = removeInlineComments(line);
                        line = trim(line);

                        if (line.empty() || line == "{")
                            continue;

                        if (line == "}")
                            break;

                        if (line.rfind("signClosure", 0) == 0)
                        {
                            const auto begin = line.find('('), end = line.rfind(')');
                            if (haveSignClosure || begin == string::npos || end == string::npos || end <= begin
                                || trim(line.substr(0, begin)) != "signClosure"
                                || (!trim(line.substr(end + 1)).empty() && trim(line.substr(end + 1)) != ";"))
                                throw invalid_argument("[DictionaryGeometry] invalid or repeated LevelSet signClosure");
                            haveSignClosure = true;
                            string contents = line.substr(begin + 1, end - begin - 1);
                            replace(contents.begin(), contents.end(), ',', ' ');
                            istringstream indices(contents);
                            string token;
                            while (indices >> token)
                            {
                                size_t consumed = 0;
                                const long long id = stoll(token, &consumed);
                                if (consumed != token.size() || id == 0
                                    || id < -static_cast<long long>(numeric_limits<int>::max())
                                    || id > numeric_limits<int>::max())
                                    throw invalid_argument("[DictionaryGeometry] signClosure requires nonzero signed 1-based Line IDs");
                                signClosure.push_back({static_cast<int>(id < 0 ? -id : id) - 1, id < 0});
                            }
                            if (signClosure.empty())
                                throw invalid_argument("[DictionaryGeometry] explicit signClosure must not be empty");
                        }
                        else if (line.rfind("lines", 0) == 0)
                        {
                            vector<SignInt> lineIndex;
                            if (parseIndex(line, lineIndex))
                            {
                                for (auto& s : lineIndex)
                                {
                                    adjustIndex(s);
                                    lsLines.push_back(s);
                                }
                            }
                        }
                        else if (!line.empty() && line.front() == '(')
                        {
                            vector<SignInt> lineIndex;
                            if (parseIndex(line, lineIndex))
                            {
                                for (auto& s : lineIndex)
                                {
                                    adjustIndex(s);
                                    lsLines.push_back(s);
                                }
                            }
                        }
                        else if (line.find("innerPoint") != string::npos)
                        {
                            Vector2D coord;
                            if (parseCoordWithUnits(line, coord, units) || parseCoord(line, coord))
                            {
                                innerPoint = coord;
                                haveInner = true;
                            }
                        }
                        else if (line.find("moving") != string::npos)
                        {
                            vector<string> tokens = split(line);
                            if (tokens.size() >= 2)
                                parseBoolToken(tokens[1], isMoving);
                        }
                    }

                    levelSetName.push_back(lsName);
                    linesOnInterface.push_back(lsLines);
                    signClosureLines.push_back(signClosure);
                    innerPoints.push_back(innerPoint);
                    moving.push_back(isMoving);

                    if (!haveInner)
                    {
                        cerr << "[DictionaryGeometry] WARNING: LevelSet '" << lsName << "' missing innerPoint; using (0,0)." << endl;
                    }
                }
            }
            continue;
        }

        if (line == "Phase")
        {
            while (getline(file, line))
            {
                line = removeInlineComments(line);
                line = trim(line);

                if (line.empty() || line == "{")
                    continue;

                if (line == "}")
                    break;

                string pName;
                if (!extractQuotedName(line, pName))
                    continue;

                vector<LevelSetOperation> ops;
                vector<SignInt> idxA;
                vector<SignInt> idxB;
                vector<SignInt> bLines;
                vector<pair<int, bool>> constraints; // (lsIndex, gt0)

                bool inConstraints = false;
                while (getline(file, line))
                {
                    line = removeInlineComments(line);
                    line = trim(line);

                    if (line.empty() || line == "{")
                        continue;

                    if (line == "}")
                    {
                        if (inConstraints)
                        {
                            inConstraints = false;
                            continue;
                        }
                        break;
                    }

                    const string lower = toLowerCopy(line);
                    if (lower.rfind("constraints", 0) == 0)
                    {
                        inConstraints = true;
                        continue;
                    }

                    if (inConstraints)
                    {
                        string lsName;
                        bool gt0 = false;
                        if (parseConstraint(line, lsName, gt0))
                        {
                            auto it = find(levelSetName.begin(), levelSetName.end(), lsName);
                            if (it == levelSetName.end())
                            {
                                cerr << "[DictionaryGeometry] WARNING: Phase '" << pName
                                     << "' references unknown level set '" << lsName << "'." << endl;
                                continue;
                            }
                            int idx = static_cast<int>(distance(levelSetName.begin(), it));
                            constraints.push_back({idx, gt0});
                        }
                        continue;
                    }

                    if (lower.rfind("boundaries", 0) == 0)
                    {
                        vector<SignInt> lineIndex;
                        if (parseIndex(line, lineIndex))
                        {
                            for (auto& s : lineIndex)
                            {
                                adjustIndex(s);
                                bLines.push_back(s);
                            }
                        }
                    }
                }

                if (constraints.empty())
                {
                    cerr << "[DictionaryGeometry] WARNING: Phase '" << pName << "' has no constraints; skipping." << endl;
                    continue;
                }
                const int baseLevelSetCount = static_cast<int>(levelSetName.size());
                for (std::size_t ci = 0; ci < constraints.size(); ++ci)
                {
                    const auto& c = constraints[ci];

                    SignInt rhs;
                    rhs.value = c.first;
                    rhs.isNegative = c.second; // >0 means flip sign so signed LS < 0 matches the constraint

                    if (ci == 0)
                    {
                        ops.push_back(LevelSetOperation::non);
                        idxA.push_back(rhs);
                        idxB.push_back(rhs);
                        continue;
                    }

                    SignInt lhs;
                    lhs.value = baseLevelSetCount + static_cast<int>(ci - 1); // previous appLevelSet index
                    lhs.isNegative = false;

                    ops.push_back(LevelSetOperation::max); // intersection: max(signed constraints) < 0
                    idxA.push_back(lhs);
                    idxB.push_back(rhs);
                }

                phaseName.push_back(pName);
                levelSetOperation.push_back(ops);
                levelSetIndexA.push_back(idxA);
                levelSetIndexB.push_back(idxB);
                boundaryLines.push_back(bLines);
            }
            continue;
        }

        if (line == "Boundary")
        {
            while (getline(file, line))
            {
                line = removeInlineComments(line);
                line = trim(line);

                if (line.empty() || line == "{")
                    continue;

                if (line == "}")
                    break;

                GeometricBoundary kind = GeometricBoundary::domain;
                if (line.find("domainBoundary") != string::npos)
                    kind = GeometricBoundary::domain;
                else if (line.find("fixedBoundary") != string::npos)
                    kind = GeometricBoundary::fixed;
                else if (line.find("movingBoundary") != string::npos)
                    kind = GeometricBoundary::moving;
                else
                    continue;

                while (getline(file, line))
                {
                    line = removeInlineComments(line);
                    line = trim(line);

                    if (line.empty() || line == "{")
                        continue;

                    if (line == "}")
                        break;

                    vector<string> tokens = split(line);
                    if (tokens.size() >= 2 && toLowerCopy(tokens[0]) == "boundary")
                    {
                        int id = stoi(tokens[1]);
                        BoundaryEntity ent;
                        ent.id = id;
                        ent.kind = kind;

                        while (getline(file, line))
                        {
                            line = removeInlineComments(line);
                            line = trim(line);

                            if (line.empty() || line == "{")
                                continue;

                            if (line == "}")
                                break;

                            vector<string> t = split(line);
                            if (t.empty())
                                continue;

                            const string key = toLowerCopy(t[0]);
                            if (key == "name" && t.size() >= 2)
                            {
                                ent.name = t[1];
                                if (ent.name.size() >= 2 && ent.name.front() == '"' && ent.name.back() == '"')
                                    ent.name = ent.name.substr(1, ent.name.size() - 2);
                            }
                            else if (key == "lines")
                            {
                                vector<SignInt> lineIndex;
                                if (parseIndex(line, lineIndex))
                                {
                                    ent.lineIds.clear();
                                    for (auto& s : lineIndex)
                                    {
                                        adjustIndex(s);
                                        ent.lineIds.push_back(s.value);
                                    }
                                }
                            }
                            else if (key == "levelsetindex" && t.size() >= 2)
                            {
                                ent.levelSetIndex = toZeroBased(stoi(t[1]));
                            }
                            else if (key == "phase" && t.size() >= 2)
                            {
                                ent.phase = toZeroBased(stoi(t[1]));
                            }
                            else if (key == "opposedphase" && t.size() >= 2)
                            {
                                ent.opposedPhase = toZeroBased(stoi(t[1]));
                            }
                            else if (key == "adjacentphase" && t.size() >= 2)
                            {
                                ent.adjacentPhase = toZeroBased(stoi(t[1]));
                            }
                            else if (key == "immersed" && t.size() >= 2)
                            {
                                bool v = true;
                                if (parseBoolToken(t[1], v))
                                    ent.immersed = v;
                            }
                            else if (key == "ignored" && t.size() >= 2)
                            {
                                bool v = false;
                                if (parseBoolToken(t[1], v))
                                    ent.ignored = v;
                            }
                        }

                        boundaries.push_back(ent);
                        boundaryById[ent.id] = ent;

                        for (int lineId : ent.lineIds)
                        {
                            if (lineId < 0 || lineId >= static_cast<int>(lines.size()))
                                continue;
                            lines[lineId].setGeoBoundary(kind);
                            if (kind != GeometricBoundary::domain)
                            {
                                lines[lineId].levelSetIndex = ent.levelSetIndex;
                                lines[lineId].phaseIndex = ent.phase;
                                lines[lineId].opposedPhaseIndex = ent.opposedPhase;
                                lines[lineId].setImmersed(ent.immersed);
                                lines[lineId].setIgnored(ent.ignored);
                            }
                            else
                            {
                                lines[lineId].phaseIndex = ent.adjacentPhase;
                            }
                        }
                    }
                }
            }
            continue;
        }

        // Retain unknown top-level blocks without interpreting their contents.
        // A nested geometry-looking name must not enter the geometry parser.
        if (line != "}") captureExtension(file, line);

    }

    // Build periodic lookup maps (boundaryId -> mate / shift).
    periodicMateByBoundary.clear();
    periodicShiftToMate.clear();

    for (const auto& pp : periodicPairs)
    {
        const int a = pp.boundaryA;
        const int b = pp.boundaryB;

        // Detect conflicting mappings.
        auto ita = periodicMateByBoundary.find(a);
        if (ita != periodicMateByBoundary.end() && ita->second != b)
        {
            cerr << "[DictionaryGeometry] WARNING: boundary " << a << " is already mapped to " << ita->second
                 << "; ignoring additional periodic mapping to " << b << "." << endl;
            continue;
        }
        auto itb = periodicMateByBoundary.find(b);
        if (itb != periodicMateByBoundary.end() && itb->second != a)
        {
            cerr << "[DictionaryGeometry] WARNING: boundary " << b << " is already mapped to " << itb->second
                 << "; ignoring additional periodic mapping to " << a << "." << endl;
            continue;
        }

        periodicMateByBoundary[a] = b;
        periodicMateByBoundary[b] = a;
        periodicShiftToMate[a] = pp.shiftAtoB;
        periodicShiftToMate[b] = pp.shiftAtoB * (-1.0);
    }

    file.close();
}
