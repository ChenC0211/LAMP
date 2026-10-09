#include "Modules/Plasma/Physics/Species.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>

namespace
{
inline std::string ToLowerCopy(const std::string& s)
{
    std::string out = s;
    std::transform(out.begin(), out.end(), out.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return out;
}

inline bool IsNearlyZero(double x)
{
    return std::abs(x) <= 1e-30;
}
} // namespace

std::string SpeciesRoleToString(SpeciesRole role)
{
    switch (role)
    {
    case SpeciesRole::kElectron:
        return "electron";
    case SpeciesRole::kIon:
        return "ion";
    case SpeciesRole::kNeutral:
        return "neutral";
    case SpeciesRole::kUnknown:
    default:
        return "unknown";
    }
}

bool ParseSpeciesRole(const std::string& token, SpeciesRole* out_role)
{
    if (!out_role)
    {
        return false;
    }

    const std::string s = ToLowerCopy(token);
    if (s == "electron" || s == "e" || s == "elec")
    {
        *out_role = SpeciesRole::kElectron;
        return true;
    }
    if (s == "ion" || s == "i")
    {
        *out_role = SpeciesRole::kIon;
        return true;
    }
    if (s == "neutral" || s == "n" || s == "gas")
    {
        *out_role = SpeciesRole::kNeutral;
        return true;
    }
    if (s == "unknown")
    {
        *out_role = SpeciesRole::kUnknown;
        return true;
    }

    return false;
}

bool Species::Validate(std::string* err) const
{
    if (name.empty())
    {
        if (err)
        {
            *err = "Species name is empty";
        }
        return false;
    }
    if (!(mass_kg > 0.0))
    {
        if (err)
        {
            *err = "Species '" + name + "' has non-positive mass";
        }
        return false;
    }
    if (!(macro_weight > 0.0))
    {
        if (err)
        {
            *err = "Species '" + name + "' has non-positive macro_weight";
        }
        return false;
    }

    // Role-specific sanity checks.
    if (role == SpeciesRole::kNeutral)
    {
        if (!IsNearlyZero(charge_C))
        {
            if (err)
            {
                *err = "Neutral species '" + name + "' should have zero charge";
            }
            return false;
        }
    }
    if (role == SpeciesRole::kElectron)
    {
        if (!(charge_C < 0.0))
        {
            if (err)
            {
                *err = "Electron species '" + name + "' should have negative charge";
            }
            return false;
        }
    }
    if (role == SpeciesRole::kIon)
    {
        if (IsNearlyZero(charge_C))
        {
            if (err)
            {
                *err = "Ion species '" + name + "' should have non-zero charge";
            }
            return false;
        }
    }

    return true;
}

void Species::PrintRow(std::ostream& os) const
{
    // Fixed-width, readable for console logs.
    os << std::left << std::setw(12) << name << "  "
       << std::left << std::setw(8) << SpeciesRoleToString(role) << "  "
       << std::scientific << std::setprecision(6)
       << std::setw(14) << mass_kg << "  "
       << std::setw(14) << charge_C << "  "
       << std::setw(12) << macro_weight << "  "
       << std::setw(14) << init.kT_J << "  "
       << "(" << init.drift_mps.x << ", " << init.drift_mps.y << ")";
}

std::string Species::RowString() const
{
    std::ostringstream oss;
    PrintRow(oss);
    return oss.str();
}
