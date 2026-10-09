#pragma once

#include <memory>
#include "Runtime/RunControls.h"
class HydroCaseCoupling;

// Immutable controls shared by the Hydro App and its EHD adapter. Fixed
// callers retain the original exact step count. With conservative StepRetry,
// steps is only the nominal count; the App's accepted time ends the loop.
struct HydroCaseSchedule
{
    int steps = 0;
    RunControls controls;

    bool writeStepFields(int step, bool includeInitialSource = false) const
    {
        // EHD evaluates the initial electric source during its first step.
        return controls.writesFieldsAtStep(step, steps) ||
            (steps > 0 && step == 1 && includeInitialSource && controls.writeFields && controls.writeInitial);
    }
    bool writeVariableStepFields(int step,double begin,double end) const
    {
        if(!controls.writeFields || step<0) return false;
        if(step==0) return controls.writeInitial;
        if((step==1 && controls.writeFirst) || (end==controls.endTime && controls.writeFinal)) return true;
        if(controls.writeEverySteps>0) return step%controls.writeEverySteps==0;
        if(controls.writeInterval<=end-begin) return true;
        return std::floor(end/controls.writeInterval+1e-12)>std::floor(begin/controls.writeInterval+1e-12);
    }
};

struct HydroCaseProgress
{
    int acceptedStep = 0; // All physical substeps/checks completed.
    int outputStep = 0;   // Scheduled fields and diagnostics also completed.
    double acceptedTime = 0.0, outputTime = 0.0;
    int nativeLeafCells=0,meshChanges=0;
};

// App-owned Hydro case state, shared by the standalone and EHD
// drivers. This is not a solver or a main loop. The caller selects every step
// and its source time; the numerical state/NS advance belong to Hydrodynamics.
// Interface owns the nodal liquid state, histories, transport and derived
// boundary report. This adapter supplies velocity layers and reads that state.
class HydroCase
{
public:
    // Read-only input/capability checks; no mesh or output yet. The optional
    // App coupling provider must outlive the case and all synchronous calls.
    explicit HydroCase(HydroCaseCoupling* coupling = nullptr);
    ~HydroCase();
    HydroCase(const HydroCase&) = delete;
    HydroCase& operator=(const HydroCase&) = delete;

    int steps() const;
    double dt() const;
    // Standalone driver entry: accepted time is authoritative when StepRetry
    // is enabled. EHD and legacy callers retain advance(step, forceTime).
    bool hasNextStep() const;
    int advanceNext();
    // Same App clock, two completion stages. Output failure can leave an
    // accepted physical step whose files are incomplete; the case then stops.
    // A failed conservative physical preparation restores its previous Hydro
    // and Interface state. Optional StepRetry is limited to the conservative
    // one-step method. Closed conservative Hydro can stage AMR with the next
    // physical step; rejected candidates retain the old mesh and both histories.
    // Legacy/EHD coupling has no rollback.
    // Other numerical paths may retain changed internal work on failure.
    HydroCaseProgress progress() const noexcept;
    int initialize();
    int advance(int step, double forceTime);
    void finish();

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};
