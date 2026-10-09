#include "Modules/Plasma/Physics/SpeciesDatabase.h"

#include <sstream>

bool SpeciesDatabase::Add(const Species& species, std::string* err)
{
    std::string localErr;
    if (!species.Validate(&localErr))
    {
        if (err)
        {
            *err = localErr;
        }
        return false;
    }

    if (index_.find(species.name) != index_.end())
    {
        if (err)
        {
            *err = "Duplicate species name: " + species.name;
        }
        return false;
    }

    const std::size_t idx = species_.size();
    species_.push_back(species);
    index_[species.name] = idx;
    order_.push_back(species.name);
    return true;
}

bool SpeciesDatabase::Has(const std::string& name) const
{
    return index_.find(name) != index_.end();
}

const Species* SpeciesDatabase::Find(const std::string& name) const
{
    auto it = index_.find(name);
    if (it == index_.end())
    {
        return nullptr;
    }
    return &species_[it->second];
}

Species* SpeciesDatabase::FindMutable(const std::string& name)
{
    auto it = index_.find(name);
    if (it == index_.end())
    {
        return nullptr;
    }
    return &species_[it->second];
}

void SpeciesDatabase::PrintTable(std::ostream& os) const
{
    os << "SpeciesDatabase: " << size() << " species\n";
    os << "Name          Role      m[kg]          q[C]          w_macro      kT[J]         drift[m/s]\n";
    os << "-----------------------------------------------------------------------------------------------\n";
    for (const auto& name : order_)
    {
        const Species* s = Find(name);
        if (!s)
        {
            continue;
        }
        s->PrintRow(os);
        os << "\n";
    }
}

std::string SpeciesDatabase::TableString() const
{
    std::ostringstream oss;
    PrintTable(oss);
    return oss.str();
}

bool SpeciesDatabase::ValidateAll(std::string* err) const
{
    for (const auto& name : order_)
    {
        const Species* s = Find(name);
        if (!s)
        {
            if (err)
            {
                *err = "Internal error: missing species '" + name + "'";
            }
            return false;
        }
        std::string localErr;
        if (!s->Validate(&localErr))
        {
            if (err)
            {
                *err = localErr;
            }
            return false;
        }
    }
    return true;
}
