#pragma once

#include <string>
#include "Modules/Interface/Compatibility/InterfaceCaseConfig.h"
#include "Modules/Interface/Compatibility/InterfaceNodalConfig.h"
#include "Modules/Interface/Compatibility/InterfaceResolutionPlan.h"

namespace DictionaryInput { class Reader; struct Node; struct ConfigurationTrace; }

// Shared conservative-interface refinement grammar. The containing App owns
// mesh placement, timing, allocation and coupled publication, not this reader.
struct InterfaceAdaptationControls
{
    struct CoarseningControls {
        double hysteresis=.5;
        double relativeVelocityError=1e-3;
    };
    InterfaceResolutionPlan::Options policy;
    int cellBudget=0;
    double remapScalarTolerance=1e-12,remapRelativeGeometryTolerance=1e-3;
    std::optional<CoarseningControls> coarsening;
};
InterfaceAdaptationControls readInterfaceAdaptation(const DictionaryInput::Node& block,
    double distanceBandCells,DictionaryInput::ConfigurationTrace& trace,bool allowCoarsening=false);

// Shared grammar for the ordinary nodal interface, including legacy inline
// entries. Returns false without consuming input for an unrelated key.
bool readNodalInterfaceEntry(const std::string& key, DictionaryInput::Reader& input,
    InterfaceNodalConfig& config);
InterfaceNodalConfig readNodalInterfaceCase(const std::string& path);

// The containing dictionary owns block placement; Interface owns its grammar.
FittedInterfaceCaseConfig readFittedInterfaceCase(const std::string& path,
    DictionaryInput::ConfigurationTrace& trace);
FittedInterfaceCaseConfig readFittedInterfaceCase(const DictionaryInput::Node& interface,
    DictionaryInput::ConfigurationTrace& trace);

InterfaceCaseConfig readInterfaceCase(const std::string& path,
    DictionaryInput::ConfigurationTrace& trace);
