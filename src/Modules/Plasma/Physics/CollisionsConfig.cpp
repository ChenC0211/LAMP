#include "Modules/Plasma/Physics/CollisionsConfig.h"

#include <algorithm>
#include <cctype>
#include <iomanip>
#include <ostream>
#include <sstream>

using namespace std;

namespace
{
string toLowerCopy(string s)
{
    transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return static_cast<char>(tolower(c)); });
    return s;
}
}

const char* ToString(CollisionsRngMode m)
{
    switch (m)
    {
        case CollisionsRngMode::kGlobal:
            return "global";
        case CollisionsRngMode::kPerParticle:
            return "perParticle";
        case CollisionsRngMode::kCounterBased:
            return "counterBased";
        default:
            return "global";
    }
}

bool ParseCollisionsRngMode(const std::string& s, CollisionsRngMode* out)
{
    if (!out)
        return false;
    const string v = toLowerCopy(s);
    if (v == "global" || v == "legacy" || v == "shared")
    {
        *out = CollisionsRngMode::kGlobal;
        return true;
    }
    if (v == "perparticle" || v == "per_particle" || v == "particle")
    {
        *out = CollisionsRngMode::kPerParticle;
        return true;
    }
    if (v == "counterbased" || v == "counter_based" || v == "counter")
    {
        *out = CollisionsRngMode::kCounterBased;
        return true;
    }
    return false;
}

bool CollisionsConfig::Validate(std::string* err) const
{
    if (!enabled)
        return true;

    if (stats_every < 0)
    {
        if (err) *err = "stats_every must be >= 0";
        return false;
    }

    if (mappings.empty())
    {
        if (err) *err = "Collisions are enabled but no mappings are provided";
        return false;
    }

    if (null_collision.majorant_Emax_eV <= 0.0)
    {
        if (err) *err = "null_collision.majorant_Emax_eV must be > 0";
        return false;
    }

    for (const auto& m : mappings)
    {
        if (m.projectile.empty())
        {
            if (err) *err = "Collision mapping has empty projectile";
            return false;
        }
        if (m.reaction_set.empty())
        {
            if (err) *err = "Collision mapping for projectile '" + m.projectile + "' has empty reaction_set";
            return false;
        }
    }
    return true;
}

void CollisionsConfig::PrintTable(std::ostream& os) const
{
    os << "CollisionsConfig\n";
    os << "  enabled                : " << (enabled ? "true" : "false") << "\n";
    os << "  parallel               : " << (parallel ? "true" : "false") << "\n";
    os << "  rng_mode               : " << ToString(rng_mode) << "\n";
    os << "  rng_seed               : " << rng_seed << "\n";
    os << "  stats_every            : " << stats_every << "\n";
    os << "  track_fast_neutrals    : " << (track_fast_neutrals ? "true" : "false") << "\n";
    os << "  use_gas_thermal_motion : " << (use_gas_thermal_motion ? "true" : "false") << "\n";
    os << "  method                 : nullCollision\n";
    os << "  majorant_Emax_eV        : " << null_collision.majorant_Emax_eV << "\n";
    os << "  safety_factor           : " << null_collision.safety_factor << "\n";
    os << "  energy_samples          : " << null_collision.energy_samples << "\n";
    os << "  clamp_energy            : " << (null_collision.clamp_energy ? "true" : "false") << "\n";
    os << "  use_fast_lookup         : " << (null_collision.use_fast_lookup ? "true" : "false") << "\n";
    os << "  lookup_bins             : " << null_collision.lookup_bins << "\n";
    os << "  mappings                : " << mappings.size() << "\n";
    for (const auto& m : mappings)
    {
        os << "    - projectile=" << m.projectile << "  reactionSet=" << m.reaction_set << "\n";
    }
}

std::string CollisionsConfig::TableString() const
{
    ostringstream ss;
    PrintTable(ss);
    return ss.str();
}
