#pragma once

#include <iosfwd>
#include <string>
#include <vector>

#include "Modules/Plasma/Physics/ReactionChannel.h"

class CrossSectionDB;

// A set of reaction channels for one projectile-target pair.
// Example: projectile="e", target="Ar".
struct ReactionSet
{
    std::string name;       // dictionary name (unique in ReactionDatabase)
    std::string projectile; // species name
    std::string target;     // species name
    std::string energy_frame = "centerOfMass"; // lookup energy frame: "lab" | "centerOfMass"

    std::vector<ReactionChannel> channels;

    bool Validate(const CrossSectionDB* xs_db, std::string* err = nullptr) const;
    void PrintSummary(std::ostream& os) const;
    std::string SummaryString() const;
};
