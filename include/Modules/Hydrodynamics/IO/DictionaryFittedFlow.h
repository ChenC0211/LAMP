#pragma once

#include "Modules/Hydrodynamics/FittedFlowCaseConfig.h"

namespace DictionaryInput { struct Node; struct ConfigurationTrace; }

// Parse only Hydro's physical block. The App selects files and adapts the
// legacy inline Interface block before calling the two module readers.
FittedFlowCaseConfig readFittedFlowCase(const DictionaryInput::Node& hydrodynamics,
    DictionaryInput::ConfigurationTrace& trace);
