#pragma once

#include <memory>

// Cross section interface.
//
// Convention:
//   - Input energy is electron-volts (eV).
//   - Output sigma is square meters (m^2).
class ICrossSectionModel
{
public:
    virtual ~ICrossSectionModel() = default;
    virtual double SigmaM2(double energy_eV) const = 0;
};

using CrossSectionModelPtr = std::shared_ptr<ICrossSectionModel>;
