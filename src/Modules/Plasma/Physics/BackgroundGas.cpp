#include "Modules/Plasma/Physics/BackgroundGas.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>

#include "Modules/Plasma/Physics/PhysicalConstants.h"

namespace
{
std::string ToLower(std::string s)
{
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return s;
}
} // namespace

BackgroundGas::BackgroundGas(std::string name, double pressure_Pa, double temperature_K)
    : name_(std::move(name)), pressure_Pa_(pressure_Pa), temperature_K_(temperature_K)
{
}

void BackgroundGas::SetParticleMassKg(double particle_mass_kg)
{
    particle_mass_kg_ = particle_mass_kg;
}

double BackgroundGas::NumberDensity_m3() const
{
    if(pressure_Pa_ <= 0.0)
        return 0.0;
    // For robustness, avoid division by zero.
    if(temperature_K_ <= 0.0)
        return 0.0;
    return pressure_Pa_ / (PhysicalConstants::kB_JpK * temperature_K_);
}

GasState BackgroundGas::Evaluate() const
{
    GasState s;
    s.name = name_;
    s.pressure_Pa = pressure_Pa_;
    s.temperature_K = temperature_K_;
    s.particle_mass_kg = particle_mass_kg_;

    s.number_density_m3 = NumberDensity_m3();

    if(particle_mass_kg_ > 0.0 && s.number_density_m3 > 0.0 && temperature_K_ > 0.0)
    {
        s.mass_density_kgpm3 = s.number_density_m3 * particle_mass_kg_;
        // Mean speed for Maxwellian: <|v|> = sqrt(8 kB T / (pi m)).
        constexpr double kPi = 3.141592653589793238462643383279502884;
        s.thermal_speed_mps = std::sqrt(8.0 * PhysicalConstants::kB_JpK * temperature_K_ / (kPi * particle_mass_kg_));
    }
    else
    {
        s.mass_density_kgpm3 = 0.0;
        s.thermal_speed_mps = 0.0;
    }

    return s;
}

bool BackgroundGas::Validate(std::string* err) const
{
    if(name_.empty())
    {
        if(err)
            *err = "BackgroundGas: name is empty.";
        return false;
    }
    if(pressure_Pa_ < 0.0)
    {
        if(err)
            *err = "BackgroundGas: pressure must be >= 0.";
        return false;
    }
    if(pressure_Pa_ > 0.0 && temperature_K_ <= 0.0)
    {
        if(err)
            *err = "BackgroundGas: temperature must be > 0 when pressure > 0.";
        return false;
    }
    if(particle_mass_kg_ < 0.0)
    {
        if(err)
            *err = "BackgroundGas: particle mass must be >= 0.";
        return false;
    }
    return true;
}

void BackgroundGas::Print(std::ostream& os) const
{
    os << TableString();
}

std::string BackgroundGas::TableString() const
{
    const GasState s = Evaluate();
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(8);
    oss << "BackgroundGas\n";
    oss << "  name            : " << s.name << "\n";
    oss << "  pressure [Pa]    : " << s.pressure_Pa << "\n";
    oss << "  temperature [K]  : " << s.temperature_K << "\n";
    oss << "  n_gas [1/m^3]    : " << s.number_density_m3 << "\n";
    if(s.particle_mass_kg > 0.0)
    {
        oss << "  m_particle [kg]  : " << s.particle_mass_kg << "\n";
        oss << "  rho_gas [kg/m^3] : " << s.mass_density_kgpm3 << "\n";
        oss << "  v_th_mean [m/s]  : " << s.thermal_speed_mps << "\n";
    }
    return oss.str();
}

bool BackgroundGas::LookupParticleMassKg(const std::string& gas_name, double* out_particle_mass_kg)
{
    if(out_particle_mass_kg == nullptr)
        return false;
    const std::string key = ToLower(gas_name);
    // Minimal built-in neutral masses for common MCC gases.
    if(key == "he" || key == "helium")
    {
        // 4.002602 amu.
        *out_particle_mass_kg = 4.002602 * PhysicalConstants::amu_kg;
        return true;
    }
    if(key == "ar" || key == "argon")
    {
        // 39.948 amu.
        *out_particle_mass_kg = 39.948 * PhysicalConstants::amu_kg;
        return true;
    }
    return false;
}
