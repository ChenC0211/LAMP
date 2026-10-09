#include "Modules/Plasma/Physics/ReactionSet.h"

#include <algorithm>
#include <iomanip>
#include <ostream>
#include <sstream>
#include <set>

#include "Modules/Plasma/Physics/CrossSectionDB.h"

using namespace std;

namespace
{
string toLowerCopy(string s)
{
    transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return static_cast<char>(tolower(c)); });
    return s;
}
}

bool ReactionSet::Validate(const CrossSectionDB* xs_db, std::string* err) const
{
    if (projectile.empty())
    {
        if (err)
            *err = "ReactionSet projectile is empty";
        return false;
    }
    if (target.empty())
    {
        if (err)
            *err = "ReactionSet target is empty";
        return false;
    }
    if (channels.empty())
    {
        if (err)
            *err = "ReactionSet has no channels";
        return false;
    }
    const string energyFrame = toLowerCopy(energy_frame);
    if (energyFrame != "lab" && energyFrame != "centerofmass")
    {
        if (err)
            *err = "ReactionSet energyFrame must be 'lab' or 'centerOfMass'";
        return false;
    }

    set<string> seen;
    for (const auto& ch : channels)
    {
        if (ch.name.empty())
        {
            if (err)
                *err = "ReactionChannel name is empty";
            return false;
        }
        const string key = toLowerCopy(ch.name);
        if (seen.count(key) != 0)
        {
            if (err)
                *err = "Duplicate channel name: " + ch.name;
            return false;
        }
        seen.insert(key);

        // Allow disabled channels to omit crossSection entirely (placeholder/feature flag).
        // When disabled, do not validate the crossSection key even if present.
        if (ch.enabled)
        {
            if (ch.cross_section_key.empty())
            {
                if (err)
                    *err = "Channel '" + ch.name + "' missing crossSection key";
                return false;
            }
        }
        if (ch.threshold_eV < 0.0)
        {
            if (err)
                *err = "Channel '" + ch.name + "' has negative threshold";
            return false;
        }

        if (ch.enabled && xs_db && !xs_db->Has(ch.cross_section_key))
        {
            if (err)
                *err = "Channel '" + ch.name + "' references unknown crossSection key: " + ch.cross_section_key;
            return false;
        }
    }

    return true;
}

void ReactionSet::PrintSummary(std::ostream& os) const
{
    os << "ReactionSet: " << (name.empty() ? "(unnamed)" : name) << "\n";
    os << "  projectile : " << projectile << "\n";
    os << "  target     : " << target << "\n";
    os << "  energyFrame: " << energy_frame << "\n";
    os << "  channels   : " << channels.size() << "\n";
    for (const auto& ch : channels)
    {
        os << "    - " << ch.name;
        os << "  xs=" << ch.cross_section_key;
        if (ch.threshold_eV > 0.0)
            os << "  Eth=" << scientific << setprecision(6) << ch.threshold_eV << " eV";
        if (!ch.products.empty())
        {
            os << "  products=";
            for (size_t i = 0; i < ch.products.size(); ++i)
            {
                if (i) os << ",";
                os << ch.products[i];
            }
        }
        if (!ch.scatter_model.empty())
            os << "  scatter=" << ch.scatter_model;
        if (!ch.energy_loss_model.empty())
            os << "  loss=" << ch.energy_loss_model;
        if (!ch.energy_share_model.empty())
            os << "  share=" << ch.energy_share_model;
        if (!ch.enabled)
            os << "  (disabled)";
        os << "\n";
    }
}

std::string ReactionSet::SummaryString() const
{
    std::ostringstream ss;
    PrintSummary(ss);
    return ss.str();
}
