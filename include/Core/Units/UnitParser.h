#pragma once

#include <string>

// Minimal "value + unit" parser.
//
// Examples:
//   "50 mTorr" -> 6.666... Pa
//   "3 cm"     -> 0.03 m
//   "500 V"    -> 500 V
//   "1e-2s"    -> 0.01 s
//   "3 eV"     -> 4.806e-19 J
//   "39.9 amu" -> 6.63e-26 kg
//   "-1 e"     -> -1.602e-19 C
//   "1e-19 m2" -> 1e-19 m^2
//
// The parser converts known units to their SI base representation.

class UnitParser
{
public:
    // Strict versioned physical fields require an explicit unit of the right
    // dimension. The existing untyped parser remains the legacy API.
    enum class Dimension { Length, Time, Voltage, Pressure, MassDensity,
        DynamicViscosity, Conductivity, SurfaceTension, SurfaceChargeDensity,
        VolumeChargeDensity, NumberDensity, Velocity, InverseTime };
    struct Result
    {
        double value_si = 0.0;
        std::string unit;   // normalized unit token (lowercase, micro normalized to 'u')
        bool has_unit = false;
    };

    // Parse a string into an SI value. Returns false if parsing fails or the
    // unit is unknown.
    static bool TryParseToSI(const std::string& text, Result& out);

    static bool TryParseQuantityToSI(const std::string& text, Dimension expected, Result& out);

    // Convenience: throws std::runtime_error on failure.
    static Result ParseToSIOrThrow(const std::string& text);

private:
    static std::string Trim(const std::string& s);
    static std::string ToLower(std::string s);
    static std::string NormalizeMicro(std::string s);
    static bool SplitNumberAndUnit(const std::string& s, std::string& number, std::string& unit);
    static bool UnitToSIMultiplier(const std::string& unit_norm, double& multiplier);
};
