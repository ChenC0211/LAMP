#include "Modules/Plasma/Physics/CrossSectionDB.h"

#include <algorithm>
#include <sstream>

bool CrossSectionDB::Register(const std::string& key,
                             CrossSectionModelPtr model,
                             bool overwrite,
                             std::string* err)
{
    if (key.empty())
    {
        if (err) *err = "Empty key";
        return false;
    }
    if (!model)
    {
        if (err) *err = "Null model";
        return false;
    }

    auto it = models_.find(key);
    if (it != models_.end() && !overwrite)
    {
        if (err) *err = "Duplicate key: " + key;
        return false;
    }

    models_[key] = std::move(model);
    return true;
}

bool CrossSectionDB::Has(const std::string& key) const
{
    return models_.find(key) != models_.end();
}

CrossSectionModelPtr CrossSectionDB::Get(const std::string& key) const
{
    auto it = models_.find(key);
    if (it == models_.end())
    {
        return nullptr;
    }
    return it->second;
}

std::vector<std::string> CrossSectionDB::Keys() const
{
    std::vector<std::string> keys;
    keys.reserve(models_.size());
    for (const auto& kv : models_)
    {
        keys.push_back(kv.first);
    }
    std::sort(keys.begin(), keys.end());
    return keys;
}

void CrossSectionDB::PrintTable(std::ostream& os) const
{
    os << "CrossSectionDB: " << models_.size() << " models\n";
    const auto keys = Keys();
    for (const auto& k : keys)
    {
        os << "  - " << k << "\n";
    }
}

std::string CrossSectionDB::TableString() const
{
    std::ostringstream oss;
    PrintTable(oss);
    return oss.str();
}
