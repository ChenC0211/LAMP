#pragma once
#include "Modules/Hydrodynamics/HydroCompatibleMomentSpace.h"
#include <memory>
#include <array>
#include <vector>

namespace Quantum::Hydrodynamics {
// Exact block factorization of the complete midpoint saddle system, ordered
// as (q, eta, pressure, boundary multiplier, optional pressure gauge).
// Full-rank G in D R=G implies eta=0. The q/boundary equations are solved
// first, then the remaining momentum equations recover physical pressure.
// This is a solver for the original system, not a velocity projection.
// Changing matrices may reuse the preceding velocity factor as a GMRES
// preconditioner. Acceptance uses the current, unpreconditioned equation.
class HydroCompatibleMomentBlockSolve {
public:
    using Matrix=HydroCompatibleMomentSpace::Matrix;
    struct Statistics {int velocityFactorizations=0,pressureFactorizations=0,krylovIterations=0,retainedVelocityUnknowns=0;};
    HydroCompatibleMomentBlockSolve();
    ~HydroCompatibleMomentBlockSolve();
    HydroCompatibleMomentBlockSolve(const HydroCompatibleMomentBlockSolve&)=delete;
    HydroCompatibleMomentBlockSolve& operator=(const HydroCompatibleMomentBlockSolve&)=delete;
    void compute(const Matrix&,int solenoidal,int complement,int pressure,int boundary,int gauge,
        const std::vector<std::array<int,9>>& cellInteriors);
    // Zero budget requests direct factorization of the current matrix.
    Eigen::VectorXd solve(const Eigen::VectorXd&,double absoluteResidualBudget=0);
    const Statistics& statistics() const noexcept;
private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};
}
