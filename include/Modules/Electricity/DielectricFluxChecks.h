#pragma once

#include <cmath>
#include <initializer_list>
#include <stdexcept>

// Independent flux reconstruction/charge closure checks for the sharp path.
// These do not select or configure the electric linear solver.
struct DielectricFluxChecks
{
    double reconstructedPoisson;
    double gaussClosure;

    explicit DielectricFluxChecks(double uniformTolerance = 1e-9)
        : reconstructedPoisson(uniformTolerance), gaussClosure(uniformTolerance) {}

    void validate() const
    {
        for (double value : {reconstructedPoisson, gaussClosure})
            if (!std::isfinite(value) || value <= 0.0 || value >= 1.0)
                throw std::invalid_argument("dielectric flux checks must be finite and in (0,1)");
    }
};
