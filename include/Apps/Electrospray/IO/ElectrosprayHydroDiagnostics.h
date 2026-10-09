#pragma once

struct FittedFlowCaseConfig;

#include "Core/IO/TecplotWriter.h"
#include <optional>
#include <string>
#include <vector>
#include "Modules/Interface/Compatibility/InterfacePhaseFlux.h"

class QuadtreeMesh;
class HydroFittedDomain;
class Hydrodynamics;
struct HydroSolveReport;
struct ElectrosprayHydroScales;
struct InterfaceMotionReport;

void reportElectrosprayHydroDomain(const HydroFittedDomain& domain);
void writeElectrosprayHydroFields(const QuadtreeMesh& mesh, const HydroFittedDomain& domain,
    const Hydrodynamics& hydro, const std::vector<int>& phaseByCell,
    const ElectrosprayHydroScales& scale, const std::string& prefix);
void writeElectrosprayPhaseFluxBudget(const HydroFittedDomain& domain,
    const InterfacePhaseFlux::FluxBudget& budget, const ElectrosprayHydroScales& scale,
    const std::string& prefix);

// SI m^3/s, positive outward at every pressure port. Console supply convention
// is converted explicitly by the driver; it does not change the native ledger.
struct ElectrosprayBoundaryFlowRatesSI
{
    double supplyOutward = 0.0, ambientOutward = 0.0, outletOutward = 0.0, netOutward = 0.0;
};
ElectrosprayBoundaryFlowRatesSI evaluateElectrosprayBoundaryFlow(const QuadtreeMesh& mesh,
    const HydroFittedDomain& domain, const Hydrodynamics& hydro,
    const FittedFlowCaseConfig& flowConfig, const ElectrosprayHydroScales& scale);

// Values already evaluated/checked by the coupling driver. This writer never
// accepts a candidate or recomputes volume, material identity or electric fields.
struct ElectrosprayMotionStepValues
{
    long double volumeSI = 0.0L, initialVolumeSI = 0.0L, expectedVolumeChangeSI = 0.0L;
    double liquidFluxSI = 0.0;
    int changedMaterials = 0;
    bool electricSuperLU = false;
    double electricResidual = 0.0;
    long double ledgerBoundarySI = 0.0L, fluxFactorSI = 0.0L, pressureDifferenceSI = 0.0L;
};

// Owns the three XY histories. Calls remain at their original
// initialization/accepted-step stages; all input reports are borrowed per call.
// Every row is flushed and closed before the accepted-step call returns.
class ElectrosprayHydroDiagnostics
{
public:
    ElectrosprayHydroDiagnostics(const std::string& prefix, bool charged, bool moving, bool ehd, bool capillary);
    void initializeMotion(const std::string& prefix);
    void writeStep(int step, double time, const HydroSolveReport& report,
        const ElectrosprayBoundaryFlowRatesSI& boundaryFlow, const ElectrosprayHydroScales& scale);
    void writeAcceptedMotion(int step, double time, const InterfaceMotionReport& report,
        const InterfacePhaseFlux::FluxBudget& phaseBudget,
        const ElectrosprayMotionStepValues& values, const ElectrosprayHydroScales& scale);
private:
    std::optional<TecplotWriter::History> diagnostics, motionDiagnostics, phaseFluxDiagnostics;
};
