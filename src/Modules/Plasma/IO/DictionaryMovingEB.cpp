#include "Modules/Plasma/IO/DictionaryMovingEB.h"
#include "Core/IO/DictionaryGeometry.h"
#include "Core/Units/UnitParser.h"
#include <algorithm>
#include <sstream>

using namespace std;
namespace {
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

bool parseLengthValueToSI(const std::string& text, double& out)
{
    UnitParser::Result r;
    if (!UnitParser::TryParseToSI(text, r))
        return false;
    out = r.value_si;
    return true;
}


class Reader : public Dictionary
{
public:
    MovingEBConfig read(const std::string& text)
    {
        std::istringstream file(text);
        MovingEBConfig movingEB;
        std::string line;
        while (getline(file, line))
        {
            line = trim(removeInlineComments(line));
            if (line == "MovingEB")
            {
                movingEB.enabled = true;
                while (getline(file, line))
                {
                    line = removeInlineComments(line);
                    line = trim(line);

                    if (line.empty() || line == "{")
                        continue;
                    if (line == "}")
                        break;

                    if (line == "ImpactModel")
                    {
                        while (getline(file, line))
                        {
                            line = removeInlineComments(line);
                            line = trim(line);

                            if (line.empty() || line == "{")
                                continue;
                            if (line == "}")
                                break;

                            auto tokens = split(line);
                            if (tokens.size() < 2)
                                continue;

                            const string key = toLowerCopy(tokens[0]);
                            if (key == "ionsonly")
                            {
                                parseBoolToken(tokens[1], movingEB.impactModel.ionsOnly);
                            }
                            else if (key == "minsamples")
                            {
                                try { movingEB.impactModel.minSamples = std::stoi(tokens[1]); } catch (...) {}
                            }
                            else if (key == "smoothing")
                            {
                                movingEB.impactModel.smoothing = tokens[1];
                            }
                            else if (key == "smoothinglength")
                            {
                                std::string expr;
                                for (size_t i = 1; i < tokens.size(); ++i)
                                {
                                    if (!expr.empty())
                                        expr += " ";
                                    expr += tokens[i];
                                }
                                parseLengthValueToSI(expr, movingEB.impactModel.smoothingLength_m);
                            }
                        }
                        continue;
                    }

                    if (line.rfind("Region", 0) == 0)
                    {
                        MovingEBRegion region;
                        while (getline(file, line))
                        {
                            line = removeInlineComments(line);
                            line = trim(line);

                            if (line.empty() || line == "{")
                                continue;
                            if (line == "}")
                                break;

                            auto tokens = split(line);
                            if (tokens.empty())
                                continue;

                            const string key = toLowerCopy(tokens[0]);
                            auto parseIntValue = [&](int& dst) {
                                if (tokens.size() >= 2)
                                {
                                    try { dst = std::stoi(tokens[1]); } catch (...) {}
                                }
                            };
                            auto parseDoubleValue = [&](double& dst) {
                                if (tokens.size() >= 2)
                                {
                                    try { dst = std::stod(tokens[1]); } catch (...) {}
                                }
                            };
                            auto parseStringValue = [&](std::string& dst) {
                                if (tokens.size() >= 2)
                                {
                                    dst = tokens[1];
                                    if (dst.size() >= 2 && dst.front() == '"' && dst.back() == '"')
                                        dst = dst.substr(1, dst.size() - 2);
                                }
                            };
                            auto parseLengthSi = [&](double& dst) {
                                if (tokens.size() < 2)
                                    return;
                                std::string expr;
                                for (size_t i = 1; i < tokens.size(); ++i)
                                {
                                    if (!expr.empty())
                                        expr += " ";
                                    expr += tokens[i];
                                }
                                parseLengthValueToSI(expr, dst);
                            };
                            auto parseLineIdList = [&](std::vector<int>& dst) {
                                dst.clear();
                                const auto l = line.find('(');
                                const auto r = line.find(')', l == string::npos ? 0 : l + 1);
                                if (l == string::npos || r == string::npos || r <= l + 1)
                                    return;
                                const string body = line.substr(l + 1, r - l - 1);
                                auto raw = split(body, ",");
                                for (auto& t : raw)
                                {
                                    t = trim(t);
                                    if (t.empty())
                                        continue;
                                    try
                                    {
                                        dst.push_back(toZeroBased(std::stoi(t)));
                                    }
                                    catch (...) {}
                                }
                            };

                            if (key == "name")
                            {
                                if (!extractQuotedName(line, region.name))
                                    parseStringValue(region.name);
                            }
                            else if (key == "boundaryid")
                            {
                                parseIntValue(region.boundaryId);
                            }
                            else if (key == "levelsetindex")
                            {
                                parseIntValue(region.levelSetIndex);
                                region.levelSetIndex = toZeroBased(region.levelSetIndex);
                            }
                            else if (key == "selector")
                            {
                                if (tokens.size() >= 2)
                                {
                                    const std::string selectorMode = toLowerCopy(tokens[1]);
                                    region.selector.wholeBoundary = (selectorMode == "wholeboundary");
                                    if (selectorMode == "lineids")
                                        parseLineIdList(region.selector.lineIds);
                                }
                            }
                            else if (key == "lineids")
                            {
                                parseLineIdList(region.selector.lineIds);
                            }
                            else if (key == "xmin")
                            {
                                parseLengthSi(region.selector.xMin_m);
                                region.selector.hasXRange = true;
                            }
                            else if (key == "xmax")
                            {
                                parseLengthSi(region.selector.xMax_m);
                                region.selector.hasXRange = true;
                            }
                            else if (key == "ymin")
                            {
                                parseLengthSi(region.selector.yMin_m);
                                region.selector.hasYRange = true;
                            }
                            else if (key == "ymax")
                            {
                                parseLengthSi(region.selector.yMax_m);
                                region.selector.hasYRange = true;
                            }
                            else if (key == "motionmode")
                            {
                                parseStringValue(region.motionMode);
                            }
                            else if (key == "model")
                            {
                                parseStringValue(region.model);
                            }
                            else if (key == "startstep")
                            {
                                parseIntValue(region.startStep);
                            }
                            else if (key == "endstep")
                            {
                                parseIntValue(region.endStep);
                            }
                            else if (key == "velocityx")
                            {
                                parseDoubleValue(region.velocityX_mps);
                            }
                            else if (key == "velocityy")
                            {
                                parseDoubleValue(region.velocityY_mps);
                            }
                            else if (key == "normalspeed")
                            {
                                parseDoubleValue(region.normalSpeed_mps);
                            }
                            else if (key == "c0")
                            {
                                parseDoubleValue(region.C0);
                            }
                            else if (key == "alpha")
                            {
                                parseDoubleValue(region.alpha);
                            }
                            else if (key == "energymodel")
                            {
                                parseStringValue(region.energyModel);
                            }
                            else if (key == "anglemodel")
                            {
                                parseStringValue(region.angleModel);
                            }
                            else if (key == "energyref")
                            {
                                parseDoubleValue(region.energyRef_eV);
                            }
                            else if (key == "angleexponent")
                            {
                                parseDoubleValue(region.angleExponent);
                            }
                            else if (key == "smoothinglength")
                            {
                                parseLengthSi(region.smoothingLength_m);
                            }
                            else if (key == "vmax")
                            {
                                parseDoubleValue(region.vmax_mps);
                            }
                        }
                        movingEB.regions.push_back(region);
                        continue;
                    }

                    auto tokens = split(line);
                    if (tokens.size() < 2)
                        continue;

                    const string key = toLowerCopy(tokens[0]);
                    auto parseIntValue = [&](int& dst) {
                        try { dst = std::stoi(tokens[1]); } catch (...) {}
                    };
                    auto parseStringValue = [&](std::string& dst) {
                        dst = tokens[1];
                        if (dst.size() >= 2 && dst.front() == '"' && dst.back() == '"')
                            dst = dst.substr(1, dst.size() - 2);
                    };
                    auto parseBoolValue = [&](bool& dst) {
                        parseBoolToken(tokens[1], dst);
                    };

                    if (key == "enabled")
                    {
                        parseBoolValue(movingEB.enabled);
                    }
                    else if (key == "updateevery")
                    {
                        parseIntValue(movingEB.updateEvery);
                    }
                    else if (key == "outputevery")
                    {
                        parseIntValue(movingEB.outputEvery);
                    }
                    else if (key == "restart")
                    {
                        parseBoolValue(movingEB.restart);
                    }
                    else if (key == "metricspath")
                    {
                        parseStringValue(movingEB.metricsPath);
                    }
                    else if (key == "consistencypath")
                    {
                        parseStringValue(movingEB.consistencyPath);
                    }
                    else if (key == "erosionmetricspath")
                    {
                        parseStringValue(movingEB.erosionMetricsPath);
                    }
                    else if (key == "impactstatspath")
                    {
                        parseStringValue(movingEB.impactStatsPath);
                    }
                    else if (key == "restartpath")
                    {
                        parseStringValue(movingEB.restartPath);
                    }
                }
                continue;
            }
        }
        return movingEB;
    }
};
}

MovingEBConfig readMovingEBConfig(const DictionaryGeometry& geometry)
{
    const auto* section = geometry.extensionBlock("MovingEB");
    return section ? Reader{}.read(*section) : MovingEBConfig{};
}
