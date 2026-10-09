#pragma once

#include <string>
#include <vector>
#include <optional>
#include "Modules/Electricity/DielectricFluxChecks.h"
#include "Modules/Electricity/SurfaceChargeRelaxation.h"

struct FittedElectricBoundarySI
{
    int boundaryId = -1;
    std::string type;
    double potentialSI = 0.0;
};

// SI electric case binding, separate from flow/material/time configuration.
struct FittedElectricCaseConfig
{
    std::string model, interfaceDiscretization, surfaceChargeModel;
    int conductorPhaseId = -1;
    double conductorPotentialSI = 0.0;
    double initialSurfaceChargeSI = 0.0;
    std::vector<FittedElectricBoundarySI> domainBoundaries;
    double potentialBounds = 0.0;
    DielectricFluxChecks fluxChecks;
    double polarizationConsistency = 0.0;
    double chargeClosure = 0.0;
    // Present only for implicitConductionMaterialTransport. Density tolerance
    // is SI C/m^2 here; the App explicitly converts it for the module kernel.
    std::optional<Quantum::Electricity::SurfaceChargeRelaxationOptions> relaxationSI;
};
