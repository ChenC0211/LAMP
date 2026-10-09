#pragma once

#include <iosfwd>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include "Modules/Plasma/Physics/CrossSectionModel.h"

class CrossSectionDB
{
public:
    CrossSectionDB() = default;

    // Register a model with a unique key.
    // Returns false on duplicate key unless overwrite==true.
    bool Register(const std::string& key, CrossSectionModelPtr model,
                  bool overwrite = false, std::string* err = nullptr);

    bool Has(const std::string& key) const;

    // Returns nullptr if not found.
    CrossSectionModelPtr Get(const std::string& key) const;

    std::vector<std::string> Keys() const;

    void PrintTable(std::ostream& os) const;
    std::string TableString() const;

private:
    std::unordered_map<std::string, CrossSectionModelPtr> models_;
};
