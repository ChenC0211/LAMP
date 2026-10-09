#pragma once

#include <iosfwd>
#include <string>

#include "Core/Base/Vector2D.h"

// Species metadata for plasma simulations.
//
// Design goals:
// - Lightweight, Core-level (reusable by PIC/MCC/fluids).
// - Store physical quantities in SI.
// - Provide basic validation and human-readable tables.

enum class SpeciesRole
{
    kUnknown = 0,
    kElectron,
    kIon,
    kNeutral
};

// String conversions for dictionaries / logs.
std::string SpeciesRoleToString(SpeciesRole role);
bool ParseSpeciesRole(const std::string& token, SpeciesRole* out_role);

// Default kinematics for initializing a species distribution.
//
// We store thermal energy as kT in Joules:
// - For neutral gases, kT = k_B T.
// - For plasmas, kT is often specified in eV; convert by kT[J] = Te[eV] * e.
struct SpeciesKinematics
{
    double kT_J = 0.0;          // thermal energy [J]
    Vector2D drift_mps{0.0, 0.0}; // drift velocity [m/s]
};

struct Species
{
    std::string name;     // unique identifier (e.g., "e", "Ar+", "Ar")
    SpeciesRole role = SpeciesRole::kUnknown;

    double mass_kg = 0.0;    // particle mass [kg]
    double charge_C = 0.0;  // particle charge [C]

    // Macro-particle weight: number of real particles represented by one macro-particle.
    // Dimensionless, must be > 0.
    double macro_weight = 1.0;

    SpeciesKinematics init;

    bool IsElectron() const { return role == SpeciesRole::kElectron; }
    bool IsIon() const { return role == SpeciesRole::kIon; }
    bool IsNeutral() const { return role == SpeciesRole::kNeutral; }

    // Basic sanity checks. Returns true if valid, false otherwise.
    // If err is not null, it will be filled with a human-readable message.
    bool Validate(std::string* err = nullptr) const;

    // Pretty-print a single-row representation.
    void PrintRow(std::ostream& os) const;
    std::string RowString() const;
};
