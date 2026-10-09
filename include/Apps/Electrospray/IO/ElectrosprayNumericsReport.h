#pragma once

struct ElectrosprayPhysicalConfig;
struct ElectrosprayConfig;
struct ElectrosprayWorkflow;
namespace Solvers { struct SolversConfig; }

// Pre-solve numerical configuration, not evidence of an accepted time step.
// Write failures propagate to the App. No field or physics state is changed.
void writeElectrosprayNumericsReport(const ElectrosprayPhysicalConfig& physical,
    const ElectrosprayConfig& run, const ElectrosprayWorkflow& workflow,
    const Solvers::SolversConfig* legacySharedSolvers);
