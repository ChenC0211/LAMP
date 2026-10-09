#pragma once

#include <string>
#include <optional>
#include "Core/IO/TecplotWriter.h"

class UnitsContext;
struct DielectricConductionReport;
struct ElectrosprayConductingState;
struct ElectrosprayConductingStepReport;
namespace Quantum::Electricity { struct SurfaceChargeFieldStep; }

// SI reporting of an already accepted conducting state. The source current
// belongs to the previous time level in explicit mode, or the accepted time
// level in implicit mode. Original numerical columns are retained in XY data;
// the constant charge-model name is stored once in zone AUXDATA.
// No charge update, normalization correction or numerical gate is added here.
class ElectrosprayChargeDiagnostics
{
public:
    ElectrosprayChargeDiagnostics(const std::string& seriesPrefix, const UnitsContext& units, bool implicit = false);
    void writeStep(int step, double timeCode, const ElectrosprayConductingState& state,
        const DielectricConductionReport& previousCurrent,
        const Quantum::Electricity::SurfaceChargeFieldStep& chargeStep,
        const ElectrosprayConductingStepReport& report);
private:
    std::string path_;
    std::optional<TecplotWriter::History> history_;
    double length_, time_, charge_, density_, current_;
    bool implicit_;
};

// Sparse console-only progress; never writes an iteration's candidate field.
void reportElectrosprayChargeIteration(int iteration, int limit, bool geometryReused, bool factorReused);
