#pragma once

#include <functional>
#include <string>
#include <vector>
#include "Modules/Electricity/SurfaceChargeTransport.h"

namespace Quantum::Electricity
{
struct SurfaceChargeRelaxationOptions
{
    double underRelaxation = 0.1; // passive_fed's charge fixed-point iteration.
    int maxIterations = 500;
    // In the caller's consistent surface-charge-density units, not coulombs.
    double densityChangeTolerance = 1e-6;
    double relativeEquationTolerance = 1e-10;
};

struct SurfaceChargeRelaxationReport
{
    bool success = false;
    std::string stage, reason;
    int iterations = 0, temporalOrder = 0;
    double alpha = 0.0, beta = 0.0;
    double maxDensityChange = 0.0;
    double maxEquationDensityResidual = 0.0;
    double equationDensityScale = 0.0, relativeEquationResidual = 0.0;
    // Sum of the actual discrete charge equation in integral units:
    // alpha*(Q_new-Q_old) + beta*dt/dt_old*(Q_old-Q_older) - dt*I_new.
    double integratedEquationResidual = 0.0;
};

struct SurfaceChargeRelaxationResult
{
    SurfaceChargeRelaxationReport report;
    // Published only on success. Current is evaluated at EXACTLY this charge,
    // not at the pre-relaxation iterate of a subsequently updated charge.
    std::vector<double> integratedCharge, integratedNormalCurrent;
};

// passive_fed/rhoe_compt time coefficients and under-relaxed fixed point,
// applied to conservative material-panel integrals. All histories must first
// be remapped onto history.panelMeasure()/panel order by the geometry owner.
// First step: backward Euler. Later: actual variable-step two-level formula.
// Surface stretching is already in the material remap: do NOT add its density
// expansion term a second time. No emission or tangential flux is inferred.
//
// evaluateCurrent(Q) solves/reconstructs the actual electric state for Q and
// returns panel-integrated normal current INTO the surface. It is synchronous;
// captures must remain valid for this call. The caller owns field and electrode
// closure checks, factor reuse and cleanup of callback state after failure.
// Both adjacent-iterate change and full time-equation residual must pass.
// Inputs are unchanged and result arrays empty on failure. This kernel neither
// advances a geometry nor commits a Field/App state or changes a time-step gate.
SurfaceChargeRelaxationResult relaxSurfaceCharge(
    const SurfaceChargeTransport& history, const std::vector<double>* previousCharge,
    double dt, double previousDt,
    const std::function<std::vector<double>(const std::vector<double>&)>& evaluateCurrent,
    const SurfaceChargeRelaxationOptions& options = {});
}
