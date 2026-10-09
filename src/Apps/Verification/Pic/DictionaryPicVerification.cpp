#include "Verification/Pic/DictionaryPicVerification.h"
#include "Core/Units/UnitParser.h"
#include <algorithm>
#include <cctype>
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

static bool isQuotedName(const string& s)
{
    return s.size() >= 2 && s.front() == '"' && s.back() == '"';
}

static string unquote(const string& s)
{
    if (isQuotedName(s))
        return s.substr(1, s.size() - 2);
    return s;
}

std::string stripListPunct(std::string s)
{
    auto isPunct = [](char c) {
        return c == '(' || c == ')' || c == '[' || c == ']' || c == '{' || c == '}' || c == ',' || c == ';';
    };
    while (!s.empty() && isPunct(s.front()))
        s.erase(s.begin());
    while (!s.empty() && isPunct(s.back()))
        s.pop_back();
    return s;
}

bool tryParseStrictDouble(const std::string& s, double& out)
{
    size_t idx = 0;
    try
    {
        out = std::stod(s, &idx);
    }
    catch (...)
    {
        return false;
    }
    return idx == s.size();
}

bool parseScalarToSI(const std::vector<std::string>& tokens, size_t start, double& out_si)
{
    if (start >= tokens.size())
        return false;

    std::string text;
    for (size_t i = start; i < tokens.size(); ++i)
    {
        if (!text.empty())
            text += " ";
        text += tokens[i];
    }

    // Pure-number: interpret as SI.
    if (tokens.size() == start + 1)
    {
        double numeric = 0.0;
        if (tryParseStrictDouble(tokens[start], numeric))
        {
            out_si = numeric;
            return true;
        }
    }

    UnitParser::Result r;
    if (!UnitParser::TryParseToSI(text, r))
        return false;
    out_si = r.value_si;
    return true;
}
}

bool DictionaryPicVerification::readDict(const std::string& path)
{
    return pic_.readDict(path, nullptr, [this](const auto& block, const auto& parent,
        const auto& line, const auto& tokens) { return readRecord(block, parent, line, tokens); });
}

bool DictionaryPicVerification::readRecord(const std::string& b, const std::string& parent,
    const std::string& line, const std::vector<std::string>& tokens)
{
        // RegressionUnitBoxSinSin block (compat: allow either
        //   PIC{ regressionUnitBoxSinSin true }
        // or
        //   RegressionUnitBoxSinSin { enabled true }
        // ).
        if (toLowerCopy(b) == "regressionunitboxsinsin")
        {
            if (tokens.size() >= 2)
            {
                const string k = toLowerCopy(tokens[0]);
                if (k == "enabled")
                    parseBoolToken(tokens[1], config_.regressionUnitBoxSinSin);
            }
            return true;
        }

        // RegressionUnitBoxSinCosPeriodicY block (Periodic Y wrap).
        //   PIC{ regressionUnitBoxSinCosPeriodicY true }
        // or
        //   RegressionUnitBoxSinCosPeriodicY { enabled true }.
        if (toLowerCopy(b) == "regressionunitboxsincosperiodicy")
        {
            if (tokens.size() >= 2)
            {
                const string k = toLowerCopy(tokens[0]);
                if (k == "enabled")
                    parseBoolToken(tokens[1], config_.regressionUnitBoxSinCosPeriodicY);
            }
            return true;
        }

        // RegressionElectricFieldSinCosPeriodicY block.
        //   PIC{ regressionElectricFieldSinCosPeriodicY true }
        // or
        //   RegressionElectricFieldSinCosPeriodicY { enabled true }.
        if (toLowerCopy(b) == "regressionelectricfieldsincosperiodicy")
        {
            if (tokens.size() >= 2)
            {
                const string k = toLowerCopy(tokens[0]);
                if (k == "enabled")
                    parseBoolToken(tokens[1], config_.regressionElectricFieldSinCosPeriodicY);
            }
            return true;
        }

        // RegressionHarness block
        if (toLowerCopy(b) == "regressionharness")
        {
            if (tokens.size() >= 2)
            {
                const string k = toLowerCopy(tokens[0]);
                if (k == "mode")
                    config_.regressionHarness.mode = toLowerCopy(tokens[1]);
                else if (k == "baseline" || k == "baselinepath")
                    config_.regressionHarness.baselinePath = unquote(tokens[1]);
                else if (k == "report" || k == "reportpath")
                    config_.regressionHarness.reportPath = unquote(tokens[1]);
                else if (k == "abstol" || k == "abs")
                    config_.regressionHarness.absTol = stod(tokens[1]);
                else if (k == "reltol" || k == "rel")
                    config_.regressionHarness.relTol = stod(tokens[1]);
            }
            return true;
        }

        // RegressionCutCircle block
        if (toLowerCopy(b) == "regressioncutcircle")
        {
            if (tokens.size() >= 2)
            {
                const string k = toLowerCopy(tokens[0]);

                if (k == "enabled")
                    parseBoolToken(tokens[1], config_.regressionCutCircle.enabled);
                else if (k == "levelsetindex" || k == "levelset")
                    config_.regressionCutCircle.levelSetIndex = stoi(tokens[1]);
                else if (k == "radius")
                    config_.regressionCutCircle.radius = stod(tokens[1]);
                else if (k == "dumpcsv")
                    parseBoolToken(tokens[1], config_.regressionCutCircle.dumpCSV);
                else if (k == "cellscsv")
                    config_.regressionCutCircle.cellsCSV = unquote(tokens[1]);
                else if (k == "segmentscsv")
                    config_.regressionCutCircle.segmentsCSV = unquote(tokens[1]);
                else if (k == "insideisnegative" || k == "insideneg")
                    parseBoolToken(tokens[1], config_.regressionCutCircle.insideIsNegative);
                else if (k == "alphasamplingn" || k == "alphasamplen")
                    config_.regressionCutCircle.alphaSamplingN = stoi(tokens[1]);
                else if (k == "cutlengthclampeps" || k == "cutlengthclamp")
                    config_.regressionCutCircle.cutLengthClampEps = stod(tokens[1]);
                else if (k == "center")
                {
                    // Expect: center (x, y)
                    parseCoord(line, config_.regressionCutCircle.center);
                }
            }
            return true;
        }



        // RegressionCutCirclePoisson block
        if (toLowerCopy(b) == "regressioncutcirclepoisson")
        {
            if (tokens.size() >= 2)
            {
                const string k = toLowerCopy(tokens[0]);

                if (k == "enabled")
                    parseBoolToken(tokens[1], config_.regressionCutCirclePoisson.enabled);
                else if (k == "levelsetindex" || k == "levelset")
                    config_.regressionCutCirclePoisson.levelSetIndex = stoi(tokens[1]);
                else if (k == "radius")
                    config_.regressionCutCirclePoisson.radius = stod(tokens[1]);
                else if (k == "insideisnegative" || k == "insideneg")
                    parseBoolToken(tokens[1], config_.regressionCutCirclePoisson.insideIsNegative);
                else if (k == "alphasamplingn" || k == "alphasamplen")
                    config_.regressionCutCirclePoisson.alphaSamplingN = stoi(tokens[1]);
                else if (k == "cutlengthclampeps" || k == "cutlengthclamp")
                    config_.regressionCutCirclePoisson.cutLengthClampEps = stod(tokens[1]);
                else if (k == "coeff" || k == "coefficient")
                    config_.regressionCutCirclePoisson.coeff = stod(tokens[1]);
                else if (k == "useembeddedbc" || k == "useebbc")
                    parseBoolToken(tokens[1], config_.regressionCutCirclePoisson.useEmbeddedBC);
                else if (k == "computeeerror" || k == "computeefielderror" || k == "computee")
                    parseBoolToken(tokens[1], config_.regressionCutCirclePoisson.computeEError);
                else if (k == "dumpphicsv")
                    parseBoolToken(tokens[1], config_.regressionCutCirclePoisson.dumpPhiCSV);
                else if (k == "phicsv")
                    config_.regressionCutCirclePoisson.phiCSV = unquote(tokens[1]);
                else if (k == "dumpmetricscsv")
                    parseBoolToken(tokens[1], config_.regressionCutCirclePoisson.dumpMetricsCSV);
                else if (k == "cellscsv")
                    config_.regressionCutCirclePoisson.cellsCSV = unquote(tokens[1]);
                else if (k == "segmentscsv")
                    config_.regressionCutCirclePoisson.segmentsCSV = unquote(tokens[1]);
                else if (k == "center")
                {
                    // Expect: center (x, y)
                    parseCoord(line, config_.regressionCutCirclePoisson.center);
                }
            }
            return true;
        }

        // RegressionCutCirclePolarizationPoisson block
        if (toLowerCopy(b) == "regressioncutcirclepolarizationpoisson" ||
            toLowerCopy(b) == "regressioncutcirclepolpoisson" ||
            toLowerCopy(b) == "regressioncutcirclepolarization")
        {
            if (tokens.size() >= 2)
            {
                const string k = toLowerCopy(tokens[0]);

                if (k == "enabled")
                    parseBoolToken(tokens[1], config_.regressionCutCirclePolarizationPoisson.enabled);
                else if (k == "levelsetindex" || k == "levelset")
                    config_.regressionCutCirclePolarizationPoisson.levelSetIndex = stoi(tokens[1]);
                else if (k == "radius")
                    config_.regressionCutCirclePolarizationPoisson.radius = stod(tokens[1]);
                else if (k == "insideisnegative" || k == "insideneg")
                    parseBoolToken(tokens[1], config_.regressionCutCirclePolarizationPoisson.insideIsNegative);
                else if (k == "alphasamplingn" || k == "alphasamplen")
                    config_.regressionCutCirclePolarizationPoisson.alphaSamplingN = stoi(tokens[1]);
                else if (k == "cutlengthclampeps" || k == "cutlengthclamp")
                    config_.regressionCutCirclePolarizationPoisson.cutLengthClampEps = stod(tokens[1]);
                else if (k == "coeff" || k == "coefficient")
                    config_.regressionCutCirclePolarizationPoisson.coeff = stod(tokens[1]);
                else if (k == "kappa")
                    config_.regressionCutCirclePolarizationPoisson.kappa = stod(tokens[1]);
                else if (k == "kappamode")
                    config_.regressionCutCirclePolarizationPoisson.kappaMode = unquote(tokens[1]);
                else if (k == "kappa0")
                    config_.regressionCutCirclePolarizationPoisson.kappa0 = stod(tokens[1]);
                else if (k == "kappaamp")
                    config_.regressionCutCirclePolarizationPoisson.kappaAmp = stod(tokens[1]);
                else if (k == "kappakx")
                    config_.regressionCutCirclePolarizationPoisson.kappaKx = stod(tokens[1]);
                else if (k == "kappaky")
                    config_.regressionCutCirclePolarizationPoisson.kappaKy = stod(tokens[1]);
                else if (k == "useembeddedbc" || k == "useebbc")
                    parseBoolToken(tokens[1], config_.regressionCutCirclePolarizationPoisson.useEmbeddedBC);
                else if (k == "computeeerror" || k == "computeefielderror" || k == "computee")
                    parseBoolToken(tokens[1], config_.regressionCutCirclePolarizationPoisson.computeEError);
                else if (k == "dumpphicsv")
                    parseBoolToken(tokens[1], config_.regressionCutCirclePolarizationPoisson.dumpPhiCSV);
                else if (k == "phicsv")
                    config_.regressionCutCirclePolarizationPoisson.phiCSV = unquote(tokens[1]);
                else if (k == "dumpmetricscsv")
                    parseBoolToken(tokens[1], config_.regressionCutCirclePolarizationPoisson.dumpMetricsCSV);
                else if (k == "cellscsv")
                    config_.regressionCutCirclePolarizationPoisson.cellsCSV = unquote(tokens[1]);
                else if (k == "segmentscsv")
                    config_.regressionCutCirclePolarizationPoisson.segmentsCSV = unquote(tokens[1]);
                else if (k == "center")
                {
                    // Expect: center (x, y)
                    parseCoord(line, config_.regressionCutCirclePolarizationPoisson.center);
                }
            }
            return true;
        }

        // RegressionEmbeddedConductor block (Test 5)
        if (toLowerCopy(b) == "regressionembeddedconductor")
        {
            if (tokens.size() >= 2)
            {
                const string k = toLowerCopy(tokens[0]);

                if (k == "enabled")
                    parseBoolToken(tokens[1], config_.regressionEmbeddedConductor.enabled);
                else if (k == "levelsetindex" || k == "levelset")
                    config_.regressionEmbeddedConductor.levelSetIndex = stoi(tokens[1]);
                else if (k == "radius")
                    config_.regressionEmbeddedConductor.radius = stod(tokens[1]);
                else if (k == "insideisnegative" || k == "insideneg")
                    parseBoolToken(tokens[1], config_.regressionEmbeddedConductor.insideIsNegative);
                else if (k == "alphasamplingn" || k == "alphasamplen")
                    config_.regressionEmbeddedConductor.alphaSamplingN = stoi(tokens[1]);
                else if (k == "cutlengthclampeps" || k == "cutlengthclamp")
                    config_.regressionEmbeddedConductor.cutLengthClampEps = stod(tokens[1]);
                else if (k == "charge" || k == "q")
                    config_.regressionEmbeddedConductor.charge = stod(tokens[1]);
                else if (k == "chargeunit")
                    config_.regressionEmbeddedConductor.chargeUnit = toLowerCopy(unquote(tokens[1]));
                else if (k == "sigma")
                    config_.regressionEmbeddedConductor.sigma = stod(tokens[1]);
                else if (k == "sourceshape")
                    config_.regressionEmbeddedConductor.sourceShape = toLowerCopy(unquote(tokens[1]));
                else if (k == "samplemode")
                    config_.regressionEmbeddedConductor.sampleMode = toLowerCopy(unquote(tokens[1]));
                else if (k == "sigmaautofromfine")
                    parseBoolToken(tokens[1], config_.regressionEmbeddedConductor.sigmaAutoFromFine);
                else if (k == "sigmafinefactor")
                    config_.regressionEmbeddedConductor.sigmaFineFactor = stod(tokens[1]);
                else if (k == "dumpphicsv")
                    parseBoolToken(tokens[1], config_.regressionEmbeddedConductor.dumpPhiCSV);
                else if (k == "phicsv")
                    config_.regressionEmbeddedConductor.phiCSV = unquote(tokens[1]);
                else if (k == "dumpencsv")
                    parseBoolToken(tokens[1], config_.regressionEmbeddedConductor.dumpEnCSV);
                else if (k == "encsv")
                    config_.regressionEmbeddedConductor.enCSV = unquote(tokens[1]);
                else if (k == "dumpmetricscsv")
                    parseBoolToken(tokens[1], config_.regressionEmbeddedConductor.dumpMetricsCSV);
                else if (k == "cellscsv")
                    config_.regressionEmbeddedConductor.cellsCSV = unquote(tokens[1]);
                else if (k == "segmentscsv")
                    config_.regressionEmbeddedConductor.segmentsCSV = unquote(tokens[1]);
                else if (k == "center")
                {
                    parseCoord(line, config_.regressionEmbeddedConductor.center);
                }
                else if (k == "particle")
                {
                    parseCoord(line, config_.regressionEmbeddedConductor.particle);
                }
            }
            return true;
        }

        // -----------------------------------------------------------------
        // SelfForceTest3 block (paperA Test 3)
        // -----------------------------------------------------------------
        if (toLowerCopy(b) == "selfforcetest3")
        {
            if (tokens.size() >= 2)
            {
                const string k = toLowerCopy(tokens[0]);
                if (k == "enabled" || k == "enable")
                {
                    parseBoolToken(tokens[1], config_.selfForceTest3.enabled);
                }
                else if (k == "subtests" || k == "tests")
                {
                    config_.selfForceTest3.subtests.clear();
                    for (size_t i = 1; i < tokens.size(); ++i)
                    {
                        std::string t = stripListPunct(unquote(tokens[i]));
                        if (t.empty())
                            continue;
                        config_.selfForceTest3.subtests.push_back(toLowerCopy(t));
                    }
                }
                else if (k == "ablations" || k == "ablation")
                {
                    config_.selfForceTest3.ablations.clear();
                    for (size_t i = 1; i < tokens.size(); ++i)
                    {
                        std::string t = stripListPunct(unquote(tokens[i]));
                        if (t.empty())
                            continue;
                        config_.selfForceTest3.ablations.push_back(t);
                    }
                }
                else if (k == "outputuniform")
                    config_.selfForceTest3.outputUniform = unquote(tokens[1]);
                else if (k == "outputtraverse")
                    config_.selfForceTest3.outputTraverse = unquote(tokens[1]);
                else if (k == "outputperiodiccross")
                    config_.selfForceTest3.outputPeriodicCross = unquote(tokens[1]);
                else if (k == "outputsummary")
                    config_.selfForceTest3.outputSummary = unquote(tokens[1]);
                else if (k == "neutrality.enabled")
                    parseBoolToken(tokens[1], config_.selfForceTest3.neutrality.enabled);
                else if (k == "gaugefix.enabled")
                    parseBoolToken(tokens[1], config_.selfForceTest3.gaugeFix.enabled);
                else if (k == "neutrality.mode")
                    config_.selfForceTest3.neutrality.mode = unquote(tokens[1]);
                else if (k == "gaugefix.mode")
                    config_.selfForceTest3.gaugeFix.mode = unquote(tokens[1]);
            }
            return true;
        }

        // SelfForceTest3::Particle block
        if (toLowerCopy(b) == "particle" && toLowerCopy(parent) == "selfforcetest3")
        {
            if (tokens.size() >= 2)
            {
                const string k = toLowerCopy(tokens[0]);
                if (k == "q" || k == "charge")
                {
                    double q_si = 0.0;
                    if (parseScalarToSI(tokens, 1, q_si))
                        config_.selfForceTest3.particle.q = q_si;
                    else
                        config_.selfForceTest3.particle.q = stod(tokens[1]);
                }
                else if (k == "m" || k == "mass")
                {
                    double m_si = 0.0;
                    if (parseScalarToSI(tokens, 1, m_si))
                        config_.selfForceTest3.particle.m = m_si;
                    else
                        config_.selfForceTest3.particle.m = stod(tokens[1]);
                }
            }
            return true;
        }

        // SelfForceTest3::Neutrality block
        if (toLowerCopy(b) == "neutrality" && toLowerCopy(parent) == "selfforcetest3")
        {
            if (tokens.size() >= 2)
            {
                const string k = toLowerCopy(tokens[0]);
                if (k == "enabled" || k == "enable")
                    parseBoolToken(tokens[1], config_.selfForceTest3.neutrality.enabled);
                else if (k == "mode")
                    config_.selfForceTest3.neutrality.mode = unquote(tokens[1]);
            }
            return true;
        }

        // SelfForceTest3::GaugeFix block
        if (toLowerCopy(b) == "gaugefix" && toLowerCopy(parent) == "selfforcetest3")
        {
            if (tokens.size() >= 2)
            {
                const string k = toLowerCopy(tokens[0]);
                if (k == "enabled" || k == "enable")
                    parseBoolToken(tokens[1], config_.selfForceTest3.gaugeFix.enabled);
                else if (k == "mode")
                    config_.selfForceTest3.gaugeFix.mode = unquote(tokens[1]);
            }
            return true;
        }

        // SelfForceTest3::Uniform block
        if (toLowerCopy(b) == "uniform" && toLowerCopy(parent) == "selfforcetest3")
        {
            if (tokens.size() >= 2)
            {
                const string k = toLowerCopy(tokens[0]);
                if (k == "nlist" || k == "n")
                {
                    config_.selfForceTest3.uniform.Nlist.clear();
                    for (size_t i = 1; i < tokens.size(); ++i)
                    {
                        std::string t = stripListPunct(tokens[i]);
                        if (t.empty())
                            continue;
                        try
                        {
                            config_.selfForceTest3.uniform.Nlist.push_back(stoi(t));
                        }
                        catch (...)
                        {
                            // Ignore malformed entries.
                        }
                    }
                }
                else if (k == "samplepoints" || k == "points")
                {
                    std::vector<Vector2D> pts;
                    if (parseCoords(line, pts))
                        config_.selfForceTest3.uniform.samplePoints = std::move(pts);
                }
            }
            return true;
        }

        // SelfForceTest3::Traverse block
        if (toLowerCopy(b) == "traverse" && toLowerCopy(parent) == "selfforcetest3")
        {
            if (tokens.size() >= 2)
            {
                const string k = toLowerCopy(tokens[0]);
                if (k == "path")
                {
                    std::vector<Vector2D> pts;
                    if (parseCoords(line, pts) && pts.size() >= 2)
                    {
                        config_.selfForceTest3.traverse.start = pts[0];
                        config_.selfForceTest3.traverse.end = pts[1];
                    }
                    std::string last = stripListPunct(tokens.back());
                    try
                    {
                        config_.selfForceTest3.traverse.numSamples = stoi(last);
                    }
                    catch (...)
                    {
                        // Ignore if not provided.
                    }
                }
            }
            return true;
        }

        // SelfForceTest3::TraverseCut block
        if (toLowerCopy(b) == "traversecut" && toLowerCopy(parent) == "selfforcetest3")
        {
            if (tokens.size() >= 2)
            {
                const string k = toLowerCopy(tokens[0]);
                if (k == "path")
                {
                    std::vector<Vector2D> pts;
                    if (parseCoords(line, pts) && pts.size() >= 2)
                    {
                        config_.selfForceTest3.traverseCut.start = pts[0];
                        config_.selfForceTest3.traverseCut.end = pts[1];
                    }
                    std::string last = stripListPunct(tokens.back());
                    try
                    {
                        config_.selfForceTest3.traverseCut.numSamples = stoi(last);
                    }
                    catch (...)
                    {
                        // Ignore if not provided.
                    }
                }
            }
            return true;
        }

        // SelfForceTest3::CutCross block
        if (toLowerCopy(b) == "cutcross" && toLowerCopy(parent) == "selfforcetest3")
        {
            if (tokens.size() >= 2)
            {
                const string k = toLowerCopy(tokens[0]);
                if (k == "path")
                {
                    std::vector<Vector2D> pts;
                    if (parseCoords(line, pts) && pts.size() >= 2)
                    {
                        config_.selfForceTest3.cutCross.start = pts[0];
                        config_.selfForceTest3.cutCross.end = pts[1];
                    }
                    std::string last = stripListPunct(tokens.back());
                    try
                    {
                        config_.selfForceTest3.cutCross.numSamples = stoi(last);
                    }
                    catch (...)
                    {
                        // Ignore if not provided.
                    }
                }
            }
            return true;
        }

        // SelfForceTest3::PeriodicCross block
        if (toLowerCopy(b) == "periodiccross" && toLowerCopy(parent) == "selfforcetest3")
        {
            if (tokens.size() >= 2)
            {
                const string k = toLowerCopy(tokens[0]);
                if (k == "path")
                {
                    std::vector<Vector2D> pts;
                    if (parseCoords(line, pts) && pts.size() >= 2)
                    {
                        config_.selfForceTest3.periodicCross.start = pts[0];
                        config_.selfForceTest3.periodicCross.end = pts[1];
                    }
                    std::string last = stripListPunct(tokens.back());
                    try
                    {
                        config_.selfForceTest3.periodicCross.numSamples = stoi(last);
                    }
                    catch (...)
                    {
                        // Ignore if not provided.
                    }
                }
            }
            return true;
        }

        if (toLowerCopy(b) == "pic" && tokens.size() == 2)
        {
            const string k = toLowerCopy(tokens[0]);
            if ((k == "regressionunitboxsinsin" || k == "regressiontest") && tokens.size() == 2)
                    parseBoolToken(tokens[1], config_.regressionUnitBoxSinSin);
                else if ((k == "regressionunitboxsincosperiodicy" || k == "regressionunitboxperiodicy") && tokens.size() == 2)
                    parseBoolToken(tokens[1], config_.regressionUnitBoxSinCosPeriodicY);
                else if ((k == "regressionelectricfieldsincosperiodicy" || k == "regressionelectricfieldperiodicy") && tokens.size() == 2)
                    parseBoolToken(tokens[1], config_.regressionElectricFieldSinCosPeriodicY);
        }
        return false;
}
