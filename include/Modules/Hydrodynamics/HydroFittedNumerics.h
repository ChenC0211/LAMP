#pragma once

#include <cmath>
#include <initializer_list>
#include <stdexcept>

// Acceptance settings for the fitted-fluid path, not material properties.
// Keep the three actual linear systems distinct from physical balance checks.
struct HydroFittedNumerics
{
    double momentumLinear;
    double pressureLinear;
    double coupledLinear;
    double momentumBalance;
    double localContinuity;
    double outletAdmissibility;

    // Uniform values preserve the pre-migration HydroSolveOptions contract.
    explicit HydroFittedNumerics(double uniformTolerance = 1e-8)
        : momentumLinear(uniformTolerance), pressureLinear(uniformTolerance),
          coupledLinear(uniformTolerance), momentumBalance(uniformTolerance),
          localContinuity(uniformTolerance), outletAdmissibility(uniformTolerance) {}

    void validate() const
    {
        for (double value : {momentumLinear, pressureLinear, coupledLinear,
                 momentumBalance, localContinuity, outletAdmissibility})
            if (!std::isfinite(value) || value <= 0.0 || value >= 1.0)
                throw std::invalid_argument("fitted Hydro numerical tolerances must be finite and in (0,1)");
    }
};
