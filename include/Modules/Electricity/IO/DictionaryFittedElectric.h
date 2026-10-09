#pragma once

#include <string>
#include "Modules/Electricity/FittedElectricCaseConfig.h"

namespace DictionaryInput { struct ConfigurationTrace; }

FittedElectricCaseConfig readFittedElectricCase(const std::string& path, DictionaryInput::ConfigurationTrace& trace);
