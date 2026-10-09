#include "Modules/Plasma/Physics/ReactionDatabase.h"

#include <algorithm>
#include <ostream>
#include <sstream>

using namespace std;

bool ReactionDatabase::Register(const ReactionSet& set, bool overwrite, std::string* err)
{
    if (set.name.empty())
    {
        if (err)
            *err = "ReactionSet name is empty";
        return false;
    }

    const auto it = sets_.find(set.name);
    if (it != sets_.end() && !overwrite)
    {
        if (err)
            *err = "Duplicate ReactionSet name: " + set.name;
        return false;
    }
    sets_[set.name] = set;
    return true;
}

bool ReactionDatabase::Has(const std::string& name) const
{
    return sets_.find(name) != sets_.end();
}

const ReactionSet* ReactionDatabase::Get(const std::string& name) const
{
    auto it = sets_.find(name);
    if (it == sets_.end())
        return nullptr;
    return &it->second;
}

std::vector<std::string> ReactionDatabase::Keys() const
{
    std::vector<std::string> keys;
    keys.reserve(sets_.size());
    for (const auto& kv : sets_)
        keys.push_back(kv.first);
    sort(keys.begin(), keys.end());
    return keys;
}

void ReactionDatabase::PrintTable(std::ostream& os) const
{
    os << "Reactions: " << sets_.size() << " set(s)\n";
    for (const auto& k : Keys())
    {
        const auto* s = Get(k);
        if (!s)
            continue;
        os << "\n";
        s->PrintSummary(os);
    }
}

std::string ReactionDatabase::TableString() const
{
    std::ostringstream ss;
    PrintTable(ss);
    return ss.str();
}
