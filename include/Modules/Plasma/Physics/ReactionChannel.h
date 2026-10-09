#pragma once

#include <string>
#include <vector>

// A single collision/reaction channel for a given projectile-target pair.
//
// This is a lightweight data model used by MCC and by transport/chemistry tooling.
// The execution models (scattering, energy partition, secondary generation, etc.)
// are implemented elsewhere; here we only store the parsed configuration.

struct ReactionChannel
{
    // Human-readable / dictionary name of the channel (unique within a ReactionSet).
    std::string name;

    // Lookup key into CrossSectionDB.
    std::string cross_section_key;

    // Threshold energy for inelastic channels, in eV.
    // For elastic channels this should be 0.
    double threshold_eV = 0.0;

    // Product species names (may be empty for elastic).
    std::vector<std::string> products;

    // Scattering model identifier (e.g. "isotropic", "backscatter").
    std::string scatter_model;

    // Optional model identifiers for energy loss and energy share.
    std::string energy_loss_model;
    std::string energy_share_model;

    bool enabled = true;
};
