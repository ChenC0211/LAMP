#pragma once

#include <iosfwd>
#include <memory>
#include <string>

class UnitsContext;

// -----------------------------------------------------------------------------
// NormalizationContext
//
// UnitsContext answers: "how to interpret inputs" (SI->code mapping).
// NormalizationContext answers: "what physical reference scales are used" and
// exposes a coherent set of derived scales for non-dimensionalization.
//
// This initial version introduces the context + dictionary + report. Physics/solver
// chains will start to use this in later extensions.
// -----------------------------------------------------------------------------

enum class NormalizationStrategy
{
    kManual = 0,
    kPlasmaFrequency = 1,
    kThermal = 2,
};

struct NormalizationConfig
{
    NormalizationStrategy strategy = NormalizationStrategy::kManual;

    // Problem dimension used by the PIC/Poisson closure.
    // 2: 2D slab (requires slabThickness_m)
    // 3: fully 3D.
    int dimension = 2;

    // Out-of-plane thickness used to close 2D volume as V = A * h0.
    // Only used when dimension==2.
    double slabThickness_m = 1.0;

    // Optional knobs for future steps.
    std::string refSpecies; // e.g. "electron"
    double Te0_eV = -1.0;   // thermal strategy
    double refCharge_C = 0.0;
    double refMass_kg = 0.0;
};

struct NormalizationScales
{
    // Primary scales (SI).
    double L0_m = 1.0;
    double T0_s = 1.0;
    double Phi0_V = 1.0;
    double n0_m3 = 1.0;

    // Derived scales (SI).
    double v0_mps = 1.0;     // L0/T0
    double E0_Vpm = 1.0;     // Phi0/L0
    double rho0_Cpm3 = 1.0;  // eps0*Phi0/L0^2
    double J0_Apm2 = 1.0;    // rho0*v0

    int dimension = 2;
    double slabThickness_m = 1.0;
};

class NormalizationContext
{
public:
    NormalizationContext();

    // Default "legacy" normalization: all reference scales are 1, dimension=2,
    // slabThickness=1. This is intentionally compatible with the baseline
    // behavior where the context is not yet used by physics kernels.
    static NormalizationContext MakeLegacy();

    // Build from an existing UnitsContext plus an explicit normalization config.
    // The primary scales are taken from UnitsContext scales.
    static NormalizationContext MakeFromUnits(const std::shared_ptr<const UnitsContext>& units,
                                              const NormalizationConfig& cfg);

    const std::shared_ptr<const UnitsContext>& units() const { return units_; }
    const NormalizationConfig& config() const { return cfg_; }
    const NormalizationScales& scales() const { return scales_; }

    // --- Common derived conversions (SI <-> dimensionless) ---
    double RhoSIToTilde(double rho_Cpm3) const;
    double RhoTildeToSI(double rho_tilde) const;

    double NumberDensitySIToTilde(double n_per_m3) const;
    double NumberDensityTildeToSI(double n_tilde) const;

    // Close 2D volume as V = A*h0. For dimension==3, this is identity.
    // The input 'volume_code' is the mesh cell volume in code units.
    // For 2D meshes in this codebase, cell "volume" is typically an area.
    double VolumeCodeToSI(double volume_code) const;

    // Pretty-print a scale table.
    void PrintTable(std::ostream& os) const;
    std::string TableString() const;

private:
    std::shared_ptr<const UnitsContext> units_; // may be null for legacy
    NormalizationConfig cfg_;
    NormalizationScales scales_;

    void finalizeFromUnits();
    void computeDerived();
};
