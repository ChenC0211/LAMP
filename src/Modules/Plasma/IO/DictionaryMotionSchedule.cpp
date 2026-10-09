#include "Modules/Plasma/IO/DictionaryMotionSchedule.h"
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
    LegacyGridMotionConfig read(const std::string& text)
    {
        std::istringstream file(text);
        std::string line;
        LegacyGridMotionConfig motionSchedule;
        while (getline(file, line))
        {
            line = trim(removeInlineComments(line));
            if (line == "MotionSchedule")
            {
                motionSchedule.enabled = true;
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
                    auto parseIntValue = [&](int& dst) {
                        try { dst = std::stoi(tokens[1]); } catch (...) {}
                    };
                    auto parseStringValue = [&](std::string& dst) {
                        dst = tokens[1];
                        if (dst.size() >= 2 && dst.front() == '"' && dst.back() == '"')
                            dst = dst.substr(1, dst.size() - 2);
                    };
                    auto parseLengthSi = [&](double& dst) {
                        std::string expr;
                        for (size_t i = 1; i < tokens.size(); ++i)
                        {
                            if (!expr.empty())
                                expr += " ";
                            expr += tokens[i];
                        }
                        parseLengthValueToSI(expr, dst);
                    };

                    if (key == "enabled")
                    {
                        parseBoolToken(tokens[1], motionSchedule.enabled);
                    }
                    else if (key == "updateevery")
                    {
                        parseIntValue(motionSchedule.updateEvery);
                    }
                    else if (key == "outputevery")
                    {
                        parseIntValue(motionSchedule.outputEvery);
                    }
                    else if (key == "metricspath")
                    {
                        parseStringValue(motionSchedule.metricsPath);
                    }
                    else if (key == "startupholdsteps")
                    {
                        parseIntValue(motionSchedule.startupHoldSteps);
                    }
                    else if (key == "stage1steps")
                    {
                        parseIntValue(motionSchedule.stage1Steps);
                    }
                    else if (key == "stage2steps")
                    {
                        parseIntValue(motionSchedule.stage2Steps);
                    }
                    else if (key == "gsainitial")
                    {
                        parseLengthSi(motionSchedule.gsaInitial_m);
                    }
                    else if (key == "gsahot")
                    {
                        parseLengthSi(motionSchedule.gsaHot_m);
                    }
                    else if (key == "gadinitial")
                    {
                        parseLengthSi(motionSchedule.gadInitial_m);
                    }
                    else if (key == "gadhot")
                    {
                        parseLengthSi(motionSchedule.gadHot_m);
                    }
                    else if (key == "dainitial")
                    {
                        parseLengthSi(motionSchedule.daInitial_m);
                    }
                    else if (key == "dafinal")
                    {
                        parseLengthSi(motionSchedule.daFinal_m);
                    }
                    else if (key == "accelshiftstage1")
                    {
                        parseLengthSi(motionSchedule.accelShiftStage1_m);
                    }
                    else if (key == "decelshiftstage1")
                    {
                        parseLengthSi(motionSchedule.decelShiftStage1_m);
                    }
                    else if (key == "accelrecessionstage2")
                    {
                        parseLengthSi(motionSchedule.accelRecessionStage2_m);
                    }
                    else if (key == "accelregionxmin")
                    {
                        parseLengthSi(motionSchedule.accelRegionXMin_m);
                    }
                    else if (key == "accelregionxmax")
                    {
                        parseLengthSi(motionSchedule.accelRegionXMax_m);
                    }
                    else if (key == "decelregionxmin")
                    {
                        parseLengthSi(motionSchedule.decelRegionXMin_m);
                    }
                    else if (key == "decelregionxmax")
                    {
                        parseLengthSi(motionSchedule.decelRegionXMax_m);
                    }
                    else if (key == "aperturehalpgapinitial" || key == "aperturehalfgapinitial")
                    {
                        parseLengthSi(motionSchedule.apertureHalfGapInitial_m);
                    }
                    else if (key == "acceldownstreamfacex")
                    {
                        parseLengthSi(motionSchedule.accelDownstreamFaceX_m);
                    }
                    else if (key == "recessionbandxmin")
                    {
                        parseLengthSi(motionSchedule.recessionBandXMin_m);
                    }
                    else if (key == "recessionbandxmax")
                    {
                        parseLengthSi(motionSchedule.recessionBandXMax_m);
                    }
                    else if (key == "recessionbandymin")
                    {
                        parseLengthSi(motionSchedule.recessionBandYMin_m);
                    }
                    else if (key == "recessionbandymax")
                    {
                        parseLengthSi(motionSchedule.recessionBandYMax_m);
                    }
                }
                continue;
            }


        }
        return motionSchedule;
    }
};
}
LegacyGridMotionConfig readMotionSchedule(const DictionaryGeometry& geometry)
{
    const auto* block=geometry.extensionBlock("MotionSchedule");
    return block ? Reader{}.read(*block) : LegacyGridMotionConfig{};
}
