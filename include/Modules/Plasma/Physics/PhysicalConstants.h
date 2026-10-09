#pragma once

#include <iosfwd>
#include <string>

// Centralized SI physical constants.
//
// Rationale:
// - MCC / reactions need consistent m, q, k_B, etc.
// - Normalization and derived scales should never hard-code constants in scattered
//   locations.
// - All values in this file are SI.
//
// Notes:
// - kB and qe are exact in the 2019 SI.
// - mu0 is no longer exact post-2019; we keep a conventional value widely used in plasma simulation.
// - Prefer the provided derived constants (eV in Joules, amu in kg) to avoid repeated conversions.

class PhysicalConstants
{
public:
    // Vacuum permittivity [F/m].
    static constexpr double epsilon0_Fpm = 8.8541878128e-12;

    // Vacuum permeability [H/m]. Conventional value used in many codes.
    static constexpr double mu0_Hpm = 1.25663706212e-6;

    // Speed of light in vacuum [m/s]. Exact.
    static constexpr double c0_mps = 299792458.0;

    // Boltzmann constant [J/K]. Exact.
    static constexpr double kB_JpK = 1.380649e-23;

    // Elementary charge [C]. Exact.
    static constexpr double qe_C = 1.602176634e-19;

    // Electron mass [kg].
    static constexpr double me_kg = 9.1093837015e-31;

    // Proton mass [kg].
    static constexpr double mp_kg = 1.67262192369e-27;

    // Atomic mass unit (Dalton) [kg].
    static constexpr double amu_kg = 1.66053906660e-27;

    // 1 electron-volt in Joules [J].
    static constexpr double eV_J = qe_C;

    // Avogadro constant [1/mol]. Exact.
    static constexpr double NA_1pmol = 6.02214076e23;

    // Universal gas constant [J/(mol*K)]. Exact via kB*NA.
    static constexpr double R_JpMolK = kB_JpK * NA_1pmol;

    // Pretty-print a constant table.
    static void PrintTable(std::ostream& os);
    static std::string TableString();
};
