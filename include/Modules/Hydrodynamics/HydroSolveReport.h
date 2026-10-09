#pragma once

#include <limits>
#include <optional>
#include <string>
#include "Modules/Hydrodynamics/HydroBoundary.h"
#include "Modules/Hydrodynamics/HydroFittedNumerics.h"
#include "Modules/Hydrodynamics/HydroCompleteConservation.h"

// Fixed-mesh Hydro options. Material fields may be constant or supplied from a
// level set. The three coupling flags select the legacy geometry-aware path;
// managed callers pass explicit force densities through advance/prepareStep.
struct HydroSolveOptions
{
    bool surfaceForces = false;
    bool evaporation = false;
    bool directForce = false;
    bool closedDomainProjection = true;
    int maxDirectForceIterations = 600;
    double directForceTolerance = 1e-4;
    double linearResidualTolerance = 1e-8;
    // Only fitted Hydro steps consume this explicit split. Absent: expand the
    // legacy linearResidualTolerance uniformly, preserving existing callers.
    std::optional<HydroFittedNumerics> fittedNumerics;
    // Complete-domain physical checks. Absent means independent 1e-8
    // defaults, never inheritance from linearResidualTolerance. Fitted
    // steps reject this option and keep their own explicit numerics.
    std::optional<HydroCompleteConservation> completeConservation;
    double accelerationX = 0.0;
    double accelerationY = 0.0;
    // Empty: unchanged closed domain. Explicit patches retain the checked
    // fixed-dt projection/coupled path, not the legacy spray projection.
    std::vector<HydroBoundaryPatch> boundaryPatches;
    // Explicit CCH pressure/force skeleton; same native conservative mass
    // and coupled solver. Uniform complete XY/RZ closed domain only, selected
    // before the first step. Native viscosity uses E^T H_mu E with the owned
    // cell measure; midpoint uses its actual phase geometry instead. This
    // alone does not select KEP or midpoint EHD.
    bool compatibleNativePressure = false;
    // One Picard iterate: centered momentum from the actual shared mass flux,
    // phase-integrated E^T H_mu E viscosity and arithmetic midpoint velocity.
    // Requires HydroMidpointData and compatible uniform XY/RZ pressure.
    // RZ presently accepts frozen single-phase material; parameter geometry
    // retains its own XY capability until RZ interface transport is connected.
    bool conservativeMidpoint = false;
};

enum class HydroStateChange { Unchanged, Advanced, MayHaveChanged };

struct HydroSolveReport
{
    // Actual complete mixed conservative step. These are discrete powers;
    // pressure power need not vanish for the historical non-adjoint C/D pair.
    // The explicit compatible variant reports its selection separately.
    bool physicalNativeMomentumUsed=false;
    bool compatibleNativePressureUsed=false;
    // Physical pressure is p^(n+1/2) for this scheme; velocity is u^(n+1).
    bool midpointNativeMomentumUsed=false;
    double physicalConvectionPower=std::numeric_limits<double>::quiet_NaN();
    double physicalExternalPower=std::numeric_limits<double>::quiet_NaN();
    double midpointMassEnergyDefect=std::numeric_limits<double>::quiet_NaN();
    double viscousMeasureFractionAdjustment=0;
    double physicalPressurePower=std::numeric_limits<double>::quiet_NaN();
    double physicalViscousPower=std::numeric_limits<double>::quiet_NaN();
    double physicalEnergyRate=std::numeric_limits<double>::quiet_NaN();
    bool success = false;
    // Physical/history fields only; scratch caches may change even when the
    // physical state is unchanged. A successful Hydro substep is not coupled
    // App acceptance. The managed advance port stops after any kernel failure.
    HydroStateChange stateChange = HydroStateChange::Unchanged;
    bool requiresReinitialization = false;
    std::string stage = "validation";
    std::string reason;
    bool momentumUsedSuperLU = false;
    bool pressureUsedSuperLU = false;
    bool coupledCorrectionUsed = false;
    bool conservativeMomentumUsed = false;
    double conservativeMassRelativeClosure = std::numeric_limits<double>::quiet_NaN();
    bool coupledUsedSuperLU = false;
    bool momentumFactorReused = false;
    bool pressureFactorReused = false;
    bool coupledFactorReused = false;
    // Preserve the original predictor and projection linear-residual meanings.
    double momentumRelativeResidual = std::numeric_limits<double>::quiet_NaN();
    double pressureRelativeResidual = std::numeric_limits<double>::quiet_NaN();
    // Separate coupled-system and full discrete momentum-balance diagnostics.
    // Prescribed normal-velocity rows are constraints, not physical momentum
    // rows. Their zero identity residual does NOT measure a wall reaction or
    // contact-line support force; these diagnostics cannot certify either.
    // Native nonuniform complete domains measure the splitting defect after
    // their existing projection, with coupledCorrectionUsed=false; that value
    // is diagnostic, not an added acceptance gate or a claim of CN closure.
    // NaN means not evaluated (including legacy paths and early failures).
    double coupledRelativeResidual = std::numeric_limits<double>::quiet_NaN();
    double momentumBalanceRelativeResidual = std::numeric_limits<double>::quiet_NaN();
    double momentumBalanceMaxForceDensity = std::numeric_limits<double>::quiet_NaN();
    double divergenceBeforeProjection = std::numeric_limits<double>::quiet_NaN();
    double divergenceAfterProjection = std::numeric_limits<double>::quiet_NaN();
    double maxSpeed = std::numeric_limits<double>::quiet_NaN();
    // Populated only by the explicit fitted-fluid port. These count assembly
    // rows, not validation cases; old complete-domain report meanings remain.
    int fittedRegularMomentumRows = 0;
    int fittedBoundaryMomentumRows = 0;
    int fittedAdaptiveMomentumRows = 0;
    // Residual solves with the SAME coupled operator after cancellation.
    // No dt/history changes; original coupled diagnostics stay unchanged.
    int fittedBalanceRefinements = 0;
    double fittedLocalContinuityBeforeRefinement = std::numeric_limits<double>::quiet_NaN();
    double fittedLocalContinuityAfterRefinement = std::numeric_limits<double>::quiet_NaN();
};

// Independent, zero-time initialization result; the projection multiplier is
// not physical pressure and is never installed in the time-step state.
struct HydroVelocityInitializationReport
{
    bool success = false;
    std::string stage = "validation";
    std::string reason;
    bool usedSuperLU = false;
    double linearRelativeResidual = std::numeric_limits<double>::quiet_NaN();
    // Before means after imposing the prescribed normal boundary velocities.
    double divergenceBeforeProjection = std::numeric_limits<double>::quiet_NaN();
    double divergenceAfterProjection = std::numeric_limits<double>::quiet_NaN();
    double maxSpeed = std::numeric_limits<double>::quiet_NaN();
};
