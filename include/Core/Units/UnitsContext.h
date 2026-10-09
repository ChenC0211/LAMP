#pragma once

#include <iosfwd>
#include <string>

// Units/normalization context.
//
// The codebase historically takes "code" values as its native representation.
// This provides a way to interpret real SI inputs and consistently convert
// them into code values using reference scales.
//
// This header intentionally stays lightweight and framework-free.

enum class UnitsMode
{
    kSI = 0,       // Inputs are interpreted as SI and converted using reference scales.
    kCodeUnits = 1 // Inputs are already in code units (identity conversion).
};

struct UnitsScales
{
    // Primary reference scales.
    double L0_m = 1.0;   // length [m]
    double T0_s = 1.0;   // time [s]
    double Phi0_V = 1.0; // potential [V]
    double n0_m3 = 1.0;  // number density [1/m^3]

    // Derived or user-provided scales.
    // If set to <= 0, UnitsContext will derive the value from the primaries.
    double v0_mps = -1.0;     // velocity [m/s]
    double E0_Vpm = -1.0;     // electric field [V/m]
    double rho0_Cpm3 = -1.0;  // charge density [C/m^3]
};

class UnitsContext
{
public:
    UnitsContext();

    // Default "legacy" behavior: interpret all values as code units.
    static UnitsContext MakeCodeUnits();

    // Interpret values as SI and convert to code using the provided scales.
    static UnitsContext MakeSI(const UnitsScales& scales);

    UnitsMode mode() const { return mode_; }
    const UnitsScales& scales() const { return scales_; }

    // --- Conversions between SI and code ---
    // These are always defined, regardless of mode. When mode==kCodeUnits,
    // the reference scales are 1 by default so conversions are identity.
    double LengthSIToCode(double meters) const;
    double LengthCodeToSI(double code_value) const;

    double TimeSIToCode(double seconds) const;
    double TimeCodeToSI(double code_value) const;

    double PotentialSIToCode(double volts) const;
    double PotentialCodeToSI(double code_value) const;

    double ElectricFieldSIToCode(double volts_per_meter) const;
    double ElectricFieldCodeToSI(double code_value) const;

    double NumberDensitySIToCode(double per_m3) const;
    double NumberDensityCodeToSI(double code_value) const;

    double ChargeDensitySIToCode(double coulomb_per_m3) const;
    double ChargeDensityCodeToSI(double code_value) const;

    // Pretty-print a scale table.
    void PrintTable(std::ostream& os) const;
    std::string TableString() const;

private:
    UnitsMode mode_ = UnitsMode::kCodeUnits;
    UnitsScales scales_;

    void finalizeDerivedScales();
};
