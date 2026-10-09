#pragma once

#include "Electrospray/Compatibility/ElectrosprayConfig.h"
#include "Electrospray/Compatibility/ElectrosprayPhysicalConfig.h"
#include "Electrospray/Setup/ElectrosprayWorkflow.h"
struct ElectrosprayModuleInputs;

// One-way conversion at the old input boundary. No module settings are ever
// reconstructed as a flat legacy configuration by the formal workflow.
ElectrosprayModuleInputs adaptLegacyElectrosprayInputs(
    const ElectrosprayPhysicalConfig& physical, const ElectrosprayConfig& run,
    const ElectrosprayWorkflow& workflow);
ElectrosprayWorkflow adaptLegacyElectrosprayWorkflow(const ElectrosprayConfig& run);
