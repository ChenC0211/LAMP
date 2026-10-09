#include "Verification/Pic/RegressionHarness.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <vector>

namespace
{
static std::string toLowerCopy(std::string s)
{
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return s;
}

static std::string trimCopy(const std::string& s)
{
    size_t b = 0;
    while (b < s.size() && std::isspace(static_cast<unsigned char>(s[b])))
        ++b;
    size_t e = s.size();
    while (e > b && std::isspace(static_cast<unsigned char>(s[e - 1])))
        --e;
    return s.substr(b, e - b);
}

static std::string removeInlineCommentsCopy(const std::string& line)
{
    const auto posHash = line.find('#');
    const auto posSlash = line.find("//");
    size_t pos = std::string::npos;
    if (posHash != std::string::npos)
        pos = posHash;
    if (posSlash != std::string::npos)
        pos = std::min(pos, posSlash);
    if (pos == std::string::npos)
        return line;
    return line.substr(0, pos);
}

static std::vector<std::string> splitWs(const std::string& line)
{
    std::istringstream iss(line);
    std::vector<std::string> t;
    std::string w;
    while (iss >> w)
        t.push_back(w);
    return t;
}

static std::string unquoteCopy(const std::string& s)
{
    if (s.size() >= 2 && s.front() == '"' && s.back() == '"')
        return s.substr(1, s.size() - 2);
    return s;
}

static std::string parentDir(const std::string& path)
{
    try
    {
        std::filesystem::path p(path);
        if (p.has_parent_path())
            return p.parent_path().string();
    }
    catch (...)
    {
    }
    return std::string();
}
} // namespace

RegressionHarness::RegressionHarness(const PicVerificationConfig::RegressionHarnessConfig& cfg)
    : mode_(parseMode(cfg.mode))
    , baselinePath_(cfg.baselinePath)
    , reportPath_(cfg.reportPath)
    , absTol_(cfg.absTol)
    , relTol_(cfg.relTol)
{
}

RegressionHarness::Mode RegressionHarness::parseMode(const std::string& s)
{
    const std::string v = toLowerCopy(s);
    if (v == "record" || v == "write")
        return Mode::Record;
    if (v == "check" || v == "verify")
        return Mode::Check;
    return Mode::Off;
}

bool RegressionHarness::nearlyEqual(double a, double b) const
{
    const double diff = std::abs(a - b);
    const double scale = std::abs(b);
    return diff <= absTol_ + relTol_ * scale;
}

bool RegressionHarness::writeBaseline(const std::string& path, const Metrics& metrics)
{
    const std::string dir = parentDir(path);
    if (!dir.empty())
        std::filesystem::create_directories(dir);

    std::ofstream out(path);
    if (!out)
        return false;

    out << "# regression baseline (auto-generated)\n";
    out << "RegressionBaseline\n";
    out << "{\n";
    out << "  version 1\n";
    out << "  metrics\n";
    out << "  {\n";
    out << std::setprecision(17);
    for (const auto& kv : metrics)
        out << "    " << kv.first << " " << kv.second << "\n";
    out << "  }\n";
    out << "}\n";
    return true;
}

bool RegressionHarness::readBaseline(const std::string& path, Metrics& out)
{
    std::ifstream in(path);
    if (!in)
        return false;

    bool inMetrics = false;
    int braceDepth = 0;
    std::string line;
    while (std::getline(in, line))
    {
        line = removeInlineCommentsCopy(line);
        line = trimCopy(line);
        if (line.empty())
            continue;

        if (line.find('{') != std::string::npos)
            ++braceDepth;
        if (line.find('}') != std::string::npos)
        {
            --braceDepth;
            if (inMetrics && braceDepth <= 2)
                inMetrics = false;
        }

        const std::vector<std::string> t = splitWs(line);
        if (t.empty())
            continue;

        if (!inMetrics)
        {
            if (toLowerCopy(t[0]) == "metrics")
            {
                inMetrics = true;
                continue;
            }
            continue;
        }

        if (t.size() >= 2)
        {
            const std::string key = unquoteCopy(t[0]);
            const double val = std::stod(t[1]);
            out[key] = val;
        }
    }
    return !out.empty();
}

bool RegressionHarness::writeReport(const std::string& path,
                                   const Metrics& baseline,
                                   const Metrics& current,
                                   const std::string& summary) const
{
    const std::string dir = parentDir(path);
    if (!dir.empty())
        std::filesystem::create_directories(dir);

    std::ofstream out(path);
    if (!out)
        return false;

    out << summary << "\n\n";
    out << std::setprecision(17);
    out << "metric,baseline,current,absDiff\n";
    for (const auto& kv : baseline)
    {
        const auto it = current.find(kv.first);
        const double cur = (it != current.end()) ? it->second : std::numeric_limits<double>::quiet_NaN();
        out << kv.first << "," << kv.second << "," << cur << "," << std::abs(cur - kv.second) << "\n";
    }
    return true;
}

int RegressionHarness::process(const Metrics& metrics) const
{
    if (mode_ == Mode::Off)
        return 0;

    if (metrics.empty())
    {
        std::cerr << "[RegressionHarness] ERROR: no metrics provided." << std::endl;
        return 2;
    }

    if (mode_ == Mode::Record)
    {
        if (!writeBaseline(baselinePath_, metrics))
        {
            std::cerr << "[RegressionHarness] ERROR: failed to write baseline: " << baselinePath_ << std::endl;
            return 2;
        }
        std::cout << "[RegressionHarness] Baseline recorded: " << baselinePath_ << std::endl;
        return 0;
    }

    Metrics baseline;
    if (!readBaseline(baselinePath_, baseline))
    {
        std::cerr << "[RegressionHarness] ERROR: missing or invalid baseline: " << baselinePath_ << std::endl;
        return 3;
    }

    int mismatches = 0;
    std::ostringstream summary;
    summary << "[RegressionHarness] mode=check  baseline=\"" << baselinePath_ << "\"";
    summary << "  absTol=" << absTol_ << "  relTol=" << relTol_ << "\n";

    for (const auto& kv : baseline)
    {
        const auto it = metrics.find(kv.first);
        if (it == metrics.end())
        {
            ++mismatches;
            summary << "  MISSING metric: " << kv.first << "\n";
            continue;
        }

        const double ref = kv.second;
        const double cur = it->second;
        if (!nearlyEqual(cur, ref))
        {
            ++mismatches;
            summary << "  MISMATCH " << kv.first << ": ref=" << ref << " cur=" << cur
                    << " absDiff=" << std::abs(cur - ref) << "\n";
        }
    }

    if (mismatches == 0)
    {
        std::cout << "[RegressionHarness] PASS (" << baseline.size() << " metrics)." << std::endl;
        return 0;
    }

    const std::string sum = summary.str();
    std::cerr << sum;
    (void)writeReport(reportPath_, baseline, metrics, sum);
    std::cerr << "[RegressionHarness] FAIL: " << mismatches << " mismatches (report: " << reportPath_ << ")" << std::endl;
    return 4;
}
