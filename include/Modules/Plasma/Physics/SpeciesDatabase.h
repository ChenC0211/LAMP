#pragma once

#include <iosfwd>
#include <string>
#include <unordered_map>
#include <vector>

#include "Modules/Plasma/Physics/Species.h"

// A small, deterministic species registry.
//
// - Names are unique keys.
// - Insertion order is preserved for printing.
// - All quantities stored in SI.
class SpeciesDatabase
{
public:
    SpeciesDatabase() = default;

    // Add a species definition.
    // Returns false if validation fails or name duplicates.
    bool Add(const Species& species, std::string* err = nullptr);

    bool Has(const std::string& name) const;

    // Find a species by name. Returns nullptr if not found.
    const Species* Find(const std::string& name) const;
    Species* FindMutable(const std::string& name);

    // Number of registered species.
    std::size_t size() const { return order_.size(); }

    // Get names in insertion order.
    const std::vector<std::string>& Names() const { return order_; }

    // Pretty-print as a table.
    void PrintTable(std::ostream& os) const;
    std::string TableString() const;

    // Validate the full database.
    // Returns false if any species is invalid.
    bool ValidateAll(std::string* err = nullptr) const;

private:
    // Insertion order for stable printing.
    std::vector<std::string> order_;

    // Name -> index into species_ vector.
    std::unordered_map<std::string, std::size_t> index_;

    // Stored definitions.
    std::vector<Species> species_;
};
