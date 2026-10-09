#pragma once

#include "Modules/Hydrodynamics/HydroPressureSystem.h"
#include <vector>

class HydroFittedDomain;

// Algebra assembled from one active-fluid domain. No solver, mesh copy, gauge
// or alternative time integrator is owned here. Hydrodynamics' checked SuperLU
// path consumes these systems after its CN/AB2 momentum predictor.
namespace HydroFittedProjection
{
class System
{
public:
    System(const System&) = delete;
    System(System&&) = default;
    const HydroPressureSystem& equations() const { return equations_; }
    // Change only the velocity-dependent RHS. A rejected input leaves the
    // previous equations intact; matrix storage and inverse mass stay fixed.
    void updateRightHandSide(const Eigen::VectorXd& predictor);
    Eigen::VectorXd correctVelocity(const Eigen::VectorXd& predictor,
        const Eigen::VectorXd& pressureIncrement) const;
private:
    explicit System(const HydroFittedDomain& domain) : domain_(&domain) {}
    friend System assemble(const HydroFittedDomain&, double,
        const std::vector<double>&, const Eigen::VectorXd&);
    // Synchronous borrow: domain must outlive this system and must not be
    // moved, assigned or changed. No mesh/operator copy or cross-step cache.
    const HydroFittedDomain* domain_;
    HydroPressureSystem equations_;
    Eigen::VectorXd inverseMass_;
};

// Face densities follow the domain's compact free velocity order. All cell
// equations are retained; physical open traction faces anchor the multiplier.
System assemble(const HydroFittedDomain& domain, double dt,
    const std::vector<double>& densityByVelocityRow,
    const Eigen::VectorXd& predictor);
System assemble(const HydroFittedDomain&&, double,
    const std::vector<double>&, const Eigen::VectorXd&) = delete;

// Unsplit momentum/continuity correction belongs to HydroCoupledBalance;
// this service only prepares and applies the fitted mass-weighted projection.
}
