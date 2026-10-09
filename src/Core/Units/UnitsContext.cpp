#include "Core/Units/UnitsContext.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <stdexcept>

UnitsContext::UnitsContext()
{
    // Legacy default: code units with identity scales.
    scales_ = UnitsScales{};
    mode_ = UnitsMode::kCodeUnits;
    finalizeDerivedScales();
}

UnitsContext UnitsContext::MakeCodeUnits()
{
    UnitsContext ctx;
    ctx.mode_ = UnitsMode::kCodeUnits;
    ctx.scales_ = UnitsScales{};
    ctx.finalizeDerivedScales();
    return ctx;
}

UnitsContext UnitsContext::MakeSI(const UnitsScales& scales)
{
    UnitsContext ctx;
    ctx.mode_ = UnitsMode::kSI;
    ctx.scales_ = scales;
    ctx.finalizeDerivedScales();
    return ctx;
}

void UnitsContext::finalizeDerivedScales()
{
    // Basic sanity.
    if(scales_.L0_m <= 0.0 || scales_.T0_s <= 0.0 || scales_.Phi0_V <= 0.0)
    {
        throw std::runtime_error("UnitsContext: reference scales must be positive (L0,T0,Phi0)");
    }

    if(scales_.v0_mps <= 0.0)
    {
        scales_.v0_mps = scales_.L0_m / scales_.T0_s;
    }
    if(scales_.E0_Vpm <= 0.0)
    {
        scales_.E0_Vpm = scales_.Phi0_V / scales_.L0_m;
    }
    // rho0 depends on a charge scale (not yet defined here).
    // Keep it explicit: if not provided, default to 1 to avoid accidental
    // physical assumptions.
    if(scales_.rho0_Cpm3 <= 0.0)
    {
        scales_.rho0_Cpm3 = 1.0;
    }
}

double UnitsContext::LengthSIToCode(double meters) const { return meters / scales_.L0_m; }
double UnitsContext::LengthCodeToSI(double code_value) const { return code_value * scales_.L0_m; }

double UnitsContext::TimeSIToCode(double seconds) const { return seconds / scales_.T0_s; }
double UnitsContext::TimeCodeToSI(double code_value) const { return code_value * scales_.T0_s; }

double UnitsContext::PotentialSIToCode(double volts) const { return volts / scales_.Phi0_V; }
double UnitsContext::PotentialCodeToSI(double code_value) const { return code_value * scales_.Phi0_V; }

double UnitsContext::ElectricFieldSIToCode(double volts_per_meter) const { return volts_per_meter / scales_.E0_Vpm; }
double UnitsContext::ElectricFieldCodeToSI(double code_value) const { return code_value * scales_.E0_Vpm; }

double UnitsContext::NumberDensitySIToCode(double per_m3) const { return per_m3 / scales_.n0_m3; }
double UnitsContext::NumberDensityCodeToSI(double code_value) const { return code_value * scales_.n0_m3; }

double UnitsContext::ChargeDensitySIToCode(double coulomb_per_m3) const { return coulomb_per_m3 / scales_.rho0_Cpm3; }
double UnitsContext::ChargeDensityCodeToSI(double code_value) const { return code_value * scales_.rho0_Cpm3; }

void UnitsContext::PrintTable(std::ostream& os) const
{
    os << std::scientific << std::setprecision(8);
    os << "UnitsContext (mode=" << (mode_ == UnitsMode::kSI ? "SI" : "CodeUnits") << ")\n";
    os << "  L0   [m]      : " << scales_.L0_m << "\n";
    os << "  T0   [s]      : " << scales_.T0_s << "\n";
    os << "  Phi0 [V]      : " << scales_.Phi0_V << "\n";
    os << "  n0   [1/m^3]  : " << scales_.n0_m3 << "\n";
    os << "  v0   [m/s]    : " << scales_.v0_mps << "\n";
    os << "  E0   [V/m]    : " << scales_.E0_Vpm << "\n";
    os << "  rho0 [C/m^3]  : " << scales_.rho0_Cpm3 << "\n";
}

std::string UnitsContext::TableString() const
{
    std::ostringstream oss;
    PrintTable(oss);
    return oss.str();
}
