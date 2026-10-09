#pragma once

#include <fstream>

struct DielectricConductionReport;
namespace Quantum::Electricity
{
struct DielectricChargeSnapshot;
struct SurfaceChargeFieldStep;
}

// Electric-module output only: reads the committed charge state and its completed current
// update. force_time labels the electric load used during this step;
// charge_time=force_time+dt labels the NEW charge, not that electric load.
// This writer neither evolves state nor imposes a physical residual threshold.
class SurfaceChargeDiagnostics
{
public:
    void initialize();
    void writeStep(const Quantum::Electricity::DielectricChargeSnapshot& charge,
        const DielectricConductionReport& current,
        const Quantum::Electricity::SurfaceChargeFieldStep& update,
        int step, double forceTime, double dt, bool writeFields);

private:
    std::ofstream summary_;
};
