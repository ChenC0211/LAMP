#pragma once

#include <string>
#include <vector>

// SI, geometry-ID input for the restored fitted two-phase flow capabilities.
// No material constants, clock, electric boundary or solver backend lives here.
struct FittedPressureBoundarySI
{
    int boundaryId = -1;
    double pressureSI = 0.0;
};

struct FittedFlowCaseConfig
{
    int liquidPhase = -1, gasPhase = -1;
    FittedPressureBoundarySI supply, outlet;
    int outletIncomingPhaseId = -1;
    std::vector<FittedPressureBoundarySI> ambient;
    std::string incomingMomentumModel, solidWallModel;
    double capillaryHalfWidthCells = 0.0;
    double momentumBalance = 0.0, localContinuity = 0.0;
    double outletAdmissibility = 0.0, phaseFluxBalance = 0.0;
};
