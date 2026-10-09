#pragma once

#include <cmath>
#include <initializer_list>
#include <stdexcept>

// Physical acceptance for complete-domain flow, independent of linear solves.
// Each consumer keeps its original norm, roundoff allowance and comparison.
// Initial continuity and boundary/outlet checks apply only with boundary patches.
struct HydroCompleteConservation
{
    double momentumBalance = 1e-8;
    double initialContinuity = 1e-8;
    double outletAdmissibility = 1e-8;
    double boundaryFluxBalance = 1e-8;

    // Only the legacy input adapter expands its former single tolerance here.
    // New inputs and direct module calls have independent fixed defaults.
    HydroCompleteConservation() = default;
    explicit HydroCompleteConservation(double value)
        : momentumBalance(value), initialContinuity(value),
          outletAdmissibility(value), boundaryFluxBalance(value) {}

    void validate() const
    {
        for (double value : {momentumBalance, initialContinuity, outletAdmissibility, boundaryFluxBalance})
            if (!std::isfinite(value) || value <= 0.0 || value >= 1.0)
                throw std::invalid_argument("complete Hydro conservation tolerances must be finite and in (0,1)");
    }
};
