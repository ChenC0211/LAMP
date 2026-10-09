#include "Modules/Plasma/Physics/TabulatedCrossSection.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <sstream>

#include "Modules/Plasma/Physics/PhysicalConstants.h"
#include "Core/Units/UnitParser.h"

namespace
{
std::string Trim(const std::string& s)
{
    std::size_t b = 0;
    while (b < s.size() && std::isspace(static_cast<unsigned char>(s[b])))
    {
        ++b;
    }
    std::size_t e = s.size();
    while (e > b && std::isspace(static_cast<unsigned char>(s[e - 1])))
    {
        --e;
    }
    return s.substr(b, e - b);
}
std::string StripComment(const std::string& line)
{
    std::size_t cut = std::string::npos;

    auto consider = [&](std::size_t p)
    {
        if (p != std::string::npos)
        {
            cut = (cut == std::string::npos) ? p : std::min(cut, p);
        }
    };

    consider(line.find('#'));
    consider(line.find(';'));
    consider(line.find("//"));

    if (cut == std::string::npos)
    {
        return line;
    }
    return line.substr(0, cut);
}
std::vector<std::string> Tokenize(const std::string& line)
{
    std::vector<std::string> toks;
    std::istringstream iss(line);
    std::string t;
    while (iss >> t)
    {
        toks.push_back(t);
    }
    return toks;
}

bool LooksLikeUnitToken(const std::string& t)
{
    for (char c : t)
    {
        if (std::isalpha(static_cast<unsigned char>(c)) || c == '^')
        {
            return true;
        }
    }
    return false;
}

bool ParseScalarWithOptionalUnit(const std::vector<std::string>& tokens,
                                std::size_t& idx,
                                const std::string& defaultUnit,
                                UnitParser::Result& out,
                                std::string* err)
{
    if (idx >= tokens.size())
    {
        if (err) *err = "Expected scalar token";
        return false;
    }

    // Try inline token first (e.g. "1e-19m2").
    {
        UnitParser::Result r;
        if (UnitParser::TryParseToSI(tokens[idx], r) && r.has_unit)
        {
            out = r;
            ++idx;
            return true;
        }
    }

    // Try "value unit" pair.
    if (idx + 1 < tokens.size() && LooksLikeUnitToken(tokens[idx + 1]))
    {
        const std::string text = tokens[idx] + " " + tokens[idx + 1];
        UnitParser::Result r;
        if (!UnitParser::TryParseToSI(text, r))
        {
            if (err) *err = "Failed to parse '" + text + "'";
            return false;
        }
        out = r;
        idx += 2;
        return true;
    }

    // Fallback: attach default unit.
    {
        const std::string text = tokens[idx] + " " + defaultUnit;
        UnitParser::Result r;
        if (!UnitParser::TryParseToSI(text, r))
        {
            if (err) *err = "Failed to parse '" + text + "'";
            return false;
        }
        out = r;
        ++idx;
        return true;
    }
}

bool UnitIsEnergy(const std::string& u)
{
    return (u == "ev" || u == "j");
}

bool UnitIsArea(const std::string& u)
{
    return (u == "m2" || u == "m^2" || u == "cm2" || u == "cm^2");
}

} // namespace

TabulatedCrossSection::TabulatedCrossSection() : options_() {}

TabulatedCrossSection::TabulatedCrossSection(const Options& options) : options_(options) {}

bool TabulatedCrossSection::LoadFromFile(const std::string& path, std::string* err)
{
    energies_eV_.clear();
    sigmas_m2_.clear();
    has_power_law_ = false;
    power_p_ = 0.0;
    power_sigma1_m2_ = 0.0;

    std::ifstream ifs(path);
    if (!ifs.is_open())
    {
        if (err) *err = "Failed to open file: " + path;
        return false;
    }

    std::string line;
    std::size_t lineNo = 0;
    while (std::getline(ifs, line))
    {
        ++lineNo;
        line = Trim(StripComment(line));
        if (line.empty())
        {
            continue;
        }

        const auto tokens = Tokenize(line);
        if (tokens.empty())
        {
            continue;
        }

        std::size_t idx = 0;
        UnitParser::Result Eres;
        UnitParser::Result Sres;
        std::string localErr;
        if (!ParseScalarWithOptionalUnit(tokens, idx, options_.default_energy_unit, Eres, &localErr))
        {
            if (err) *err = "Line " + std::to_string(lineNo) + ": " + localErr;
            return false;
        }
        if (!ParseScalarWithOptionalUnit(tokens, idx, options_.default_sigma_unit, Sres, &localErr))
        {
            if (err) *err = "Line " + std::to_string(lineNo) + ": " + localErr;
            return false;
        }
        if (idx != tokens.size())
        {
            // Allow trailing tokens only if they are comments (already stripped).
            if (err)
            {
                *err = "Line " + std::to_string(lineNo) + ": Expected exactly two scalars";
            }
            return false;
        }

        if (!Eres.has_unit || !UnitIsEnergy(Eres.unit))
        {
            if (err) *err = "Line " + std::to_string(lineNo) + ": Energy must use unit eV or J";
            return false;
        }
        if (!Sres.has_unit || !UnitIsArea(Sres.unit))
        {
            if (err) *err = "Line " + std::to_string(lineNo) + ": Sigma must use area unit (m2/cm2)";
            return false;
        }

        const double E_eV = Eres.value_si / PhysicalConstants::eV_J;
        const double sigma_m2 = Sres.value_si;

        energies_eV_.push_back(E_eV);
        sigmas_m2_.push_back(sigma_m2);
    }

    return Validate(err);
}

bool TabulatedCrossSection::LoadFromColumns(const std::vector<double>& energies_eV,
                                           const std::vector<double>& sigmas_m2,
                                           std::string* err)
{
    energies_eV_ = energies_eV;
    sigmas_m2_ = sigmas_m2;
    has_power_law_ = false;
    power_p_ = 0.0;
    power_sigma1_m2_ = 0.0;
    return Validate(err);
}

bool TabulatedCrossSection::SetPowerLaw(double p, double sigma1_m2, std::string* err)
{
    energies_eV_.clear();
    sigmas_m2_.clear();
    has_power_law_ = true;
    power_p_ = p;
    power_sigma1_m2_ = sigma1_m2;
    return Validate(err);
}

bool TabulatedCrossSection::Validate(std::string* err) const
{
    if (has_power_law_)
    {
        if (!(power_sigma1_m2_ >= 0.0))
        {
            if (err) *err = "Power-law sigma1 must be non-negative";
            return false;
        }
        if (options_.interp == Interpolation::kLogLog)
        {
            if (!(power_sigma1_m2_ > 0.0))
            {
                if (err) *err = "LogLog interpolation requires sigma1>0 for power-law";
                return false;
            }
        }
        // BOLSIG+ power-law option requires p > -2, but we keep it permissive.
        return true;
    }

    if (energies_eV_.empty())
    {
        if (err) *err = "No data points";
        return false;
    }
    if (energies_eV_.size() != sigmas_m2_.size())
    {
        if (err) *err = "Internal error: size mismatch";
        return false;
    }
    if (energies_eV_.size() < 2)
    {
        if (err) *err = "Need at least 2 data points";
        return false;
    }
    for (std::size_t i = 0; i < energies_eV_.size(); ++i)
    {
        const double E = energies_eV_[i];
        const double S = sigmas_m2_[i];
        if (!(E >= 0.0))
        {
            if (err) *err = "Energy must be non-negative";
            return false;
        }
        if (!(S >= 0.0))
        {
            if (err) *err = "Sigma must be non-negative";
            return false;
        }
        if (options_.interp == Interpolation::kLogLog)
        {
            if (!(E > 0.0))
            {
                if (err) *err = "LogLog interpolation requires E>0";
                return false;
            }
            if (!(S > 0.0))
            {
                if (err) *err = "LogLog interpolation requires sigma>0";
                return false;
            }
        }
        if (i > 0)
        {
            if (!(energies_eV_[i] > energies_eV_[i - 1]))
            {
                if (err) *err = "Energies must be strictly increasing";
                return false;
            }
        }
    }
    return true;
}

double TabulatedCrossSection::SigmaM2(double energy_eV) const
{
    return SigmaM2WithHint(energy_eV, nullptr);
}

double TabulatedCrossSection::SigmaM2WithHint(double energy_eV, std::size_t* i1_hint) const
{
    if (has_power_law_)
    {
        if (energy_eV <= 0.0)
        {
            return 0.0;
        }
        return power_sigma1_m2_ * std::pow(energy_eV, power_p_);
    }

    if (energies_eV_.empty())
    {
        return 0.0;
    }

    const double Emin = energies_eV_.front();
    const double Emax = energies_eV_.back();

    if (energy_eV <= Emin)
    {
        if (options_.clamp_outside_range)
        {
            if (i1_hint)
            {
                *i1_hint = 1;
            }
            return sigmas_m2_.front();
        }
        // Extrapolate with first segment.
    }
    if (energy_eV >= Emax)
    {
        if (options_.clamp_outside_range)
        {
            if (i1_hint)
            {
                *i1_hint = energies_eV_.size() - 1;
            }
            return sigmas_m2_.back();
        }
        // Extrapolate with last segment.
    }

    const std::size_t n = energies_eV_.size();
    // Find i1 such that energies[i1-1] < E <= energies[i1]. (lower_bound index)
    std::size_t i1 = 0;
    if (i1_hint && (*i1_hint > 0) && (*i1_hint < n))
    {
        i1 = *i1_hint;
    }
    else
    {
        auto it = std::lower_bound(energies_eV_.begin(), energies_eV_.end(), energy_eV);
        i1 = static_cast<std::size_t>(it - energies_eV_.begin());
    }

    // Local walk to correct i1 if the hint was off.
    // Ensure i1 is in [1, n-1] for interpolation/extrapolation using neighboring points.
    if (i1 == 0) i1 = 1;
    if (i1 >= n) i1 = n - 1;

    while (i1 < n && energies_eV_[i1] < energy_eV)
    {
        ++i1;
        if (i1 >= n)
        {
            i1 = n;
            break;
        }
    }
    while (i1 > 0 && energies_eV_[i1 - 1] >= energy_eV)
    {
        --i1;
    }

    if (i1_hint)
    {
        // Clamp to a usable range.
        if (i1 == 0) i1 = 1;
        if (i1 >= n) i1 = n - 1;
        *i1_hint = i1;
    }

    if (i1 == 0)
    {
        // E <= Emin and no clamp: extrapolate with first interval.
        const double E0 = energies_eV_[0];
        const double E1 = energies_eV_[1];
        const double S0 = sigmas_m2_[0];
        const double S1 = sigmas_m2_[1];
        const double t = (energy_eV - E0) / (E1 - E0);
        if (options_.interp == Interpolation::kLogLog)
        {
            const double le = std::log(std::max(energy_eV, 1e-300));
            const double le0 = std::log(E0);
            const double le1 = std::log(E1);
            const double ls0 = std::log(S0);
            const double ls1 = std::log(S1);
            const double tt = (le - le0) / (le1 - le0);
            return std::exp(ls0 + (ls1 - ls0) * tt);
        }
        return S0 + (S1 - S0) * t;
    }
    if (i1 >= n)
    {
        // E >= Emax and no clamp: extrapolate with last interval.
        const double E0 = energies_eV_[n - 2];
        const double E1 = energies_eV_[n - 1];
        const double S0 = sigmas_m2_[n - 2];
        const double S1 = sigmas_m2_[n - 1];
        const double t = (energy_eV - E0) / (E1 - E0);
        if (options_.interp == Interpolation::kLogLog)
        {
            const double le = std::log(std::max(energy_eV, 1e-300));
            const double le0 = std::log(E0);
            const double le1 = std::log(E1);
            const double ls0 = std::log(S0);
            const double ls1 = std::log(S1);
            const double tt = (le - le0) / (le1 - le0);
            return std::exp(ls0 + (ls1 - ls0) * tt);
        }
        return S0 + (S1 - S0) * t;
    }

    const std::size_t i0 = i1 - 1;

    const double E0 = energies_eV_[i0];
    const double E1 = energies_eV_[i1];
    const double S0 = sigmas_m2_[i0];
    const double S1 = sigmas_m2_[i1];

    if (energy_eV == E1)
    {
        return S1;
    }

    if (options_.interp == Interpolation::kLogLog)
    {
        const double le = std::log(energy_eV);
        const double le0 = std::log(E0);
        const double le1 = std::log(E1);
        const double ls0 = std::log(S0);
        const double ls1 = std::log(S1);
        const double t = (le - le0) / (le1 - le0);
        return std::exp(ls0 + (ls1 - ls0) * t);
    }

    const double t = (energy_eV - E0) / (E1 - E0);
    return S0 + (S1 - S0) * t;
}

std::string TabulatedCrossSection::SummaryString() const
{
    std::ostringstream oss;
    if (has_power_law_)
    {
        oss << "TabulatedCrossSection: powerlaw sigma(E)=sigma1*E^p";
        oss << " sigma1=" << power_sigma1_m2_ << " m2";
        oss << " p=" << power_p_;
        oss << " interp=" << (options_.interp == Interpolation::kLogLog ? "loglog" : "lin");
        return oss.str();
    }

    oss << "TabulatedCrossSection: n=" << size();
    if (!energies_eV_.empty())
    {
        oss << " Emin=" << energies_eV_.front() << " eV";
        oss << " Emax=" << energies_eV_.back() << " eV";
    }
    oss << " interp=" << (options_.interp == Interpolation::kLogLog ? "loglog" : "lin");
    return oss.str();
}
