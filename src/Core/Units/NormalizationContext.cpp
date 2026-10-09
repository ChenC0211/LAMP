#include "Core/Units/NormalizationContext.h"

#include <cmath>
#include <iomanip>
#include <sstream>
#include <stdexcept>

#include "Core/Units/UnitsContext.h"

NormalizationContext::NormalizationContext()
{
    // Legacy default.
    units_.reset();
    cfg_ = NormalizationConfig{};
    scales_ = NormalizationScales{};
    computeDerived();
}

NormalizationContext NormalizationContext::MakeLegacy()
{
    NormalizationContext ctx;
    ctx.units_.reset();
    ctx.cfg_ = NormalizationConfig{};
    ctx.scales_ = NormalizationScales{};
    ctx.computeDerived();
    return ctx;
}

NormalizationContext NormalizationContext::MakeFromUnits(const std::shared_ptr<const UnitsContext>& units,
                                                         const NormalizationConfig& cfg)
{
    NormalizationContext ctx;
    ctx.units_ = units;
    ctx.cfg_ = cfg;
    ctx.finalizeFromUnits();
    ctx.computeDerived();
    return ctx;
}

void NormalizationContext::finalizeFromUnits()
{
    if (!units_)
    {
        scales_ = NormalizationScales{};
        scales_.dimension = cfg_.dimension;
        scales_.slabThickness_m = cfg_.slabThickness_m;
        return;
    }

    const auto& us = units_->scales();
    scales_.L0_m = us.L0_m;
    scales_.T0_s = us.T0_s;
    scales_.Phi0_V = us.Phi0_V;
    scales_.n0_m3 = us.n0_m3;
    scales_.dimension = cfg_.dimension;
    scales_.slabThickness_m = cfg_.slabThickness_m;
}

void NormalizationContext::computeDerived()
{
    // Sanity.
    if (scales_.L0_m <= 0.0 || scales_.T0_s <= 0.0 || scales_.Phi0_V <= 0.0)
        throw std::runtime_error("NormalizationContext: reference scales must be positive (L0,T0,Phi0)");

    if (cfg_.dimension != 2 && cfg_.dimension != 3)
        throw std::runtime_error("NormalizationContext: dimension must be 2 or 3");

    if (cfg_.dimension == 2 && cfg_.slabThickness_m <= 0.0)
        throw std::runtime_error("NormalizationContext: slabThickness must be positive for 2D");

    scales_.v0_mps = scales_.L0_m / scales_.T0_s;
    scales_.E0_Vpm = scales_.Phi0_V / scales_.L0_m;
    // Local epsilon0 to keep Core independent of Plasma physics.
    constexpr double kEpsilon0_Fpm = 8.8541878128e-12;
    scales_.rho0_Cpm3 = kEpsilon0_Fpm * scales_.Phi0_V / (scales_.L0_m * scales_.L0_m);
    scales_.J0_Apm2 = scales_.rho0_Cpm3 * scales_.v0_mps;

    scales_.dimension = cfg_.dimension;
    scales_.slabThickness_m = cfg_.slabThickness_m;
}

double NormalizationContext::RhoSIToTilde(double rho_Cpm3) const
{
    return rho_Cpm3 / scales_.rho0_Cpm3;
}

double NormalizationContext::RhoTildeToSI(double rho_tilde) const
{
    return rho_tilde * scales_.rho0_Cpm3;
}

double NormalizationContext::NumberDensitySIToTilde(double n_per_m3) const
{
    return n_per_m3 / scales_.n0_m3;
}

double NormalizationContext::NumberDensityTildeToSI(double n_tilde) const
{
    return n_tilde * scales_.n0_m3;
}

double NormalizationContext::VolumeCodeToSI(double volume_code) const
{
    // Mesh geometric quantities are stored in code length units.
    // For 2D, the "cell volume" is an area in code units.
    const double L0 = scales_.L0_m;
    if (scales_.dimension == 2)
    {
        return volume_code * (L0 * L0) * scales_.slabThickness_m;
    }
    return volume_code * (L0 * L0 * L0);
}

void NormalizationContext::PrintTable(std::ostream& os) const
{
    os << std::scientific << std::setprecision(8);
    os << "NormalizationContext\n";
    os << "  strategy        : ";
    switch (cfg_.strategy)
    {
    case NormalizationStrategy::kManual: os << "manual"; break;
    case NormalizationStrategy::kPlasmaFrequency: os << "plasmaFrequency"; break;
    case NormalizationStrategy::kThermal: os << "thermal"; break;
    default: os << "unknown"; break;
    }
    os << "\n";
    os << "  dimension       : " << scales_.dimension << "\n";
    if (scales_.dimension == 2)
        os << "  slabThickness[m]: " << scales_.slabThickness_m << "\n";

    os << "Scales (SI)\n";
    os << "  L0   [m]      : " << scales_.L0_m << "\n";
    os << "  T0   [s]      : " << scales_.T0_s << "\n";
    os << "  Phi0 [V]      : " << scales_.Phi0_V << "\n";
    os << "  n0   [1/m^3]  : " << scales_.n0_m3 << "\n";
    os << "  v0   [m/s]    : " << scales_.v0_mps << "\n";
    os << "  E0   [V/m]    : " << scales_.E0_Vpm << "\n";
    os << "  rho0 [C/m^3]  : " << scales_.rho0_Cpm3 << "\n";
    os << "  J0   [A/m^2]  : " << scales_.J0_Apm2 << "\n";

    if (!cfg_.refSpecies.empty())
        os << "  refSpecies     : " << cfg_.refSpecies << "\n";
    if (cfg_.refCharge_C > 0.0)
        os << "  refCharge [C]  : " << cfg_.refCharge_C << "\n";
    if (cfg_.refMass_kg > 0.0)
        os << "  refMass [kg]   : " << cfg_.refMass_kg << "\n";
    if (cfg_.Te0_eV > 0.0)
        os << "  Te0   [eV]     : " << cfg_.Te0_eV << "\n";
}

std::string NormalizationContext::TableString() const
{
    std::ostringstream oss;
    PrintTable(oss);
    return oss.str();
}
