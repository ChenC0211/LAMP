#pragma once

#include <string>
#include <vector>
#include "Electrospray/Setup/ElectrosprayNumericalControls.h"
#include "Core/Materials/MaterialTable.h"
#include "Modules/Hydrodynamics/FittedFlowCaseConfig.h"
#include "Modules/Electricity/FittedElectricCaseConfig.h"

// App input only: physical numbers are SI; phase/boundary IDs are 1-based.
// No mesh, field, normalization or solver state is owned by this reader.
struct ElectrosprayPhysicalConfig
{
    // Shared/module-owned data definitions. This flat legacy binding remains
    // the migration adapter around the existing App setup/solver calls, not a
    // versioned all-physics dictionary reader.
    using PhaseProperties = PhaseMaterial;
    using DomainElectricBoundary = FittedElectricBoundarySI;
    using PressureBoundary = FittedPressureBoundarySI;

    std::string model;
    std::string interfaceDiscretization;
    std::string fluidInterfaceLevelSet;
    int movingBoundaryId = -1;
    PhaseProperties liquid;
    PhaseProperties gas;
    int conductorPhaseId = -1;
    double conductorPotentialSI = 0.0; // V
    double surfaceTensionSI = 0.0; // N/m
    double initialSurfaceChargeSI = 0.0; // C/m^2
    std::vector<DomainElectricBoundary> domainElectricBoundaries;
    PressureBoundary supply;
    PressureBoundary outlet;
    // Optional 1-based liquid/gas phase ID for outlet inflow. -1 retains the
    // original PressureOutlet with no incoming state (backflow is rejected).
    int outletIncomingPhaseId = -1;
    // Ambient gas reservoirs may exchange mass in either direction. Supply
    // incoming material is liquid, ambient incoming material is gas; neither
    // prescribes normal velocity or flow rate. Outlet inflow is separately opt-in.
    std::vector<PressureBoundary> ambient;
    std::string incomingMomentumModel;
    std::string solidWallModel;
    std::string contactLineModel; // pinned: positions fixed, no prescribed angle.
    std::vector<int> pinnedContactPointIds; // 1-based geometry Point IDs.
    ElectrosprayNumericalControls numerics;
    // Provenance only; numerical consumers use the independently named values.
    std::string legacyToleranceSource;

    // Requires one complete ElectrosprayPhysical block. Missing/duplicate or
    // unknown keys, malformed numbers and invalid values throw; no defaults
    // silently select physical conditions. Geometry ID existence is App-owned.
    static ElectrosprayPhysicalConfig read(const std::string& path = "./system/electricity.dict");
};
