#pragma once

#include "Electrospray/Setup/ElectrosprayMaterials.h"
#include <string>
#include <optional>
#include "Core/IO/TecplotWriter.h"
class UnitsContext;
struct ElectrosprayMeshReport;
enum class ElectrosprayMeshStage;

// One compact row per accepted topology event; no duplicate spatial output.
class ElectrosprayMeshDiagnostics
{
public:
    ElectrosprayMeshDiagnostics(const std::string& prefix,
        const ElectrosprayHydroScales& scale, const UnitsContext& units);
    void write(int acceptedStep, double timeCode, const ElectrosprayMeshReport& report);
    static void progress(ElectrosprayMeshStage stage, int cells);
    // A negative threshold denotes an unselected indicator (not printed).
    static void flowPolicy(double speedThresholdSI, double gradientThresholdSI);
private:
    std::string path_;
    std::optional<TecplotWriter::History> history_;
    ElectrosprayHydroScales scale_;
    double chargeScale_;
};
