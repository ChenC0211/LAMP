#pragma once

#include <iosfwd>
#include <string>

// Background gas description (single-species ideal gas).
//
// Purpose:
// - Provide a standardized "gas / background" model for MCC and reactions.
// - Given gas name, pressure p, temperature T, compute number density
//   n = p / (k_B T) (ideal gas) and other derived properties.
//
// All stored quantities are SI.

struct GasState
{
    std::string name;

    double pressure_Pa = 0.0;      // p [Pa]
    double temperature_K = 0.0;    // T [K]

    double number_density_m3 = 0.0;    // n [1/m^3]

    // Optional: if particle_mass_kg is known (>0), we also compute mass density
    // and a representative thermal speed.
    double particle_mass_kg = 0.0;     // [kg]
    double mass_density_kgpm3 = 0.0;   // rho = n*m [kg/m^3]
    double thermal_speed_mps = 0.0;    // mean thermal speed ~ sqrt(8 kB T / (pi m)) [m/s]
};

class BackgroundGas
{
public:
    BackgroundGas() = default;
    BackgroundGas(std::string name, double pressure_Pa, double temperature_K);

    // If known, set the particle mass [kg] (e.g., Ar: 39.948 amu).
    void SetParticleMassKg(double particle_mass_kg);

    const std::string& Name() const { return name_; }
    double PressurePa() const { return pressure_Pa_; }
    double TemperatureK() const { return temperature_K_; }
    double ParticleMassKg() const { return particle_mass_kg_; }

    // Ideal gas number density n = p / (k_B T). Returns 0 if p<=0.
    double NumberDensity_m3() const;

    // Evaluate all derived properties (n, rho, vth). Always fills n.
    GasState Evaluate() const;

    // Basic validation. Returns true if valid, else false and fills err.
    bool Validate(std::string* err = nullptr) const;

    // Pretty-print.
    void Print(std::ostream& os) const;
    std::string TableString() const;

    // Minimal built-in lookup for common gases. Returns true if found.
    // Currently covers at least Argon ("Ar" / "argon").
    static bool LookupParticleMassKg(const std::string& gas_name, double* out_particle_mass_kg);

private:
    std::string name_;
    double pressure_Pa_ = 0.0;
    double temperature_K_ = 0.0;
    double particle_mass_kg_ = 0.0;
};
