#pragma once

#include "Electrospray/Setup/ElectrosprayModuleInputs.h"
#include "Electrospray/Setup/ElectrosprayWorkflow.h"
#include "Core/IO/ConfigurationTrace.h"
#include <memory>
#include <optional>
#include <string>
#include "Modules/Electricity/SurfaceChargeRelaxation.h"
class UnitsContext;
class NormalizationContext;

// Read-only migration binding. Readers are owned by Core/materials, module IO,
// solver IO and Runtime respectively; only this App joins their identities.
struct ElectrosprayVersionedInput
{
    ElectrosprayModuleInputs modules;
    // Selected once by the App reader; execution and resolved IO share it.
    ElectrosprayProblem problem = ElectrosprayProblem::Conducting;
    std::shared_ptr<const UnitsContext> units;
    std::shared_ptr<const NormalizationContext> normalization;
    DictionaryInput::ConfigurationTrace trace;
    // Code-unit iteration controls, distinct from the old all-physics adapter.
    std::optional<Quantum::Electricity::SurfaceChargeRelaxationOptions> chargeRelaxation;
};
ElectrosprayVersionedInput readElectrosprayVersionedInput(const std::string& systemDirectory = "./system");
