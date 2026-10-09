#pragma once

#include <memory>
#include <string>

#include "Solvers/Poisson/PoissonFVSolverCellCentered.h"

namespace Solvers
{

// ----------------------------------------------------------------------------
// SolversConfig
//
// Global solver-layer configuration parsed from ./system/solvers.dict.
// - This file contains *solver algorithm/linear solver knobs*.
// - Physics modules still own *equation coefficients / RHS / BC values* and
//   pass them into solver calls.
//
// Two intended usage modes:
//   (A) Primary: the case entry loads solvers.dict once (runCase) and passes the
//       parsed config to apps/modules via RunConfig.
//   (B) Optional: Advanced user code may construct a DictionarySolvers (or a
//       SolversConfig) directly and call solver APIs without any Module.
// ----------------------------------------------------------------------------
struct SolversConfig
{
    // Options for the cell-centered FV Poisson solver.
    PoissonFVSolverCellCentered::Options poissonCellCentered;
};

} // namespace Solvers
