#pragma once

#include <cmath>
#include <functional>
#include <string>

#include "Modules/Plasma/Physics/CrossSectionModel.h"

// A cross-section model backed by a user-provided callable.
//
// Notes:
//  - The callable must be well-defined for energy_eV >= 0.
//  - This wrapper applies basic sanitation (NaN/negative -> 0).
class FunctionalCrossSection final : public ICrossSectionModel
{
public:
    using Func = std::function<double(double)>;

    explicit FunctionalCrossSection(Func f, std::string name = "")
        : f_(std::move(f)), name_(std::move(name))
    {
    }

    double SigmaM2(double energy_eV) const override
    {
        if (!(energy_eV >= 0.0) || !std::isfinite(energy_eV))
            return 0.0;
        const double s = f_(energy_eV);
        if (!std::isfinite(s) || s < 0.0)
            return 0.0;
        return s;
    }

    const std::string& name() const { return name_; }

private:
    Func f_;
    std::string name_;
};
