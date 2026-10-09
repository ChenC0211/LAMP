#pragma once

#include <cstddef>
#include <string>
#include <vector>

#include "Modules/Plasma/Physics/CrossSectionModel.h"

class TabulatedCrossSection final : public ICrossSectionModel
{
public:
    enum class Interpolation
    {
        kLinear,
        kLogLog,
    };

    struct Options
    {
        Interpolation interp = Interpolation::kLinear;
        std::string default_energy_unit = "eV";
        std::string default_sigma_unit = "m2";
        bool clamp_outside_range = true;
    };
    TabulatedCrossSection();
    explicit TabulatedCrossSection(const Options& options);

    bool LoadFromFile(const std::string& path, std::string* err = nullptr);

    // Load from already-parsed (E[eV], sigma[m^2]) columns.
    bool LoadFromColumns(const std::vector<double>& energies_eV,
                         const std::vector<double>& sigmas_m2,
                         std::string* err = nullptr);

    // Special-case support for the BOLSIG+ power-law shorthand: sigma(E) = sigma1 * E^p
    // where sigma1 is in m^2 at E=1 eV.
    bool SetPowerLaw(double p, double sigma1_m2, std::string* err = nullptr);

    double SigmaM2(double energy_eV) const override;

    // Fast-path evaluation for hot loops.
    //
    // This function evaluates the same tabulated/interpolated cross section as SigmaM2(),
    // but avoids a full binary search by using a caller-provided hint for the lower_bound
    // index (i1): the smallest index such that energies_eV()[i1] >= energy_eV.
    //
    // The hint is treated as an initial guess and will be adjusted (by a small local walk)
    // to the correct lower_bound index; as long as the returned i1 is correct, this
    // method returns bitwise-identical results to SigmaM2() for the same energy.
    //
    // If i1_hint is nullptr, this method falls back to a binary search.
    double SigmaM2WithHint(double energy_eV, std::size_t* i1_hint) const;

    std::size_t size() const { return energies_eV_.size(); }
    const std::vector<double>& energies_eV() const { return energies_eV_; }
    const std::vector<double>& sigmas_m2() const { return sigmas_m2_; }

    bool Validate(std::string* err = nullptr) const;
    std::string SummaryString() const;

private:
    Options options_;
    std::vector<double> energies_eV_;
    std::vector<double> sigmas_m2_;

    bool has_power_law_ = false;
    double power_p_ = 0.0;
    double power_sigma1_m2_ = 0.0;
};
