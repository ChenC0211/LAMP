#pragma once

#include <iosfwd>
#include <string>
#include <unordered_map>
#include <vector>

#include "Modules/Plasma/Physics/ReactionSet.h"

class ReactionDatabase
{
public:
    ReactionDatabase() = default;

    // Register a reaction set by name.
    // Returns false if name exists unless overwrite==true.
    bool Register(const ReactionSet& set, bool overwrite = false, std::string* err = nullptr);

    bool Has(const std::string& name) const;
    const ReactionSet* Get(const std::string& name) const;

    std::vector<std::string> Keys() const;

    void PrintTable(std::ostream& os) const;
    std::string TableString() const;

private:
    std::unordered_map<std::string, ReactionSet> sets_;
};
