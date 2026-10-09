#pragma once

#include <string>
#include <optional>
#include "Core/IO/TecplotWriter.h"
#include "Electrospray/Coupling/ElectrosprayInterfaceMotion.h"

class HydroFittedDomain;
namespace InterfacePhaseFlux { struct FluxBudget; }
namespace InterfaceVolumeCorrection { struct Report; struct Progress; }

// App-level output of an already computed local liquid-volume balance.
// This writer converts code volume/phi to SI; it does not correct geometry,
// recompute phase fluxes, or impose a conservation acceptance threshold.
// Each call checks both flush and close, including the appended step series.
// An undefined predictor ratio (zero transport, positive defect) has an
// explicit validity flag and a zero display placeholder, not a zero error.
class ElectrosprayVolumeDiagnostics
{
public:
    ElectrosprayVolumeDiagnostics(const std::string& seriesPrefix, double lengthSI);
    static void writeProgress(const InterfaceVolumeCorrection::Progress& progress);
    void writeStep(int step, double timeSI, const InterfaceVolumeCorrection::Report& report);
    void writeFields(const std::string& fieldPrefix, const HydroFittedDomain& domain,
        const InterfaceVolumeMeasure& oldMeasure,
        const InterfacePhaseFlux::FluxBudget& oldFlux,
        const InterfaceVolumeMeasure& newMeasure,
        const InterfacePhaseFlux::FluxBudget& newFlux, double dtCode);

private:
    std::string seriesPath_;
    std::optional<TecplotWriter::History> history_;
    double lengthSI_;
    long double volumeFactor_;
};
