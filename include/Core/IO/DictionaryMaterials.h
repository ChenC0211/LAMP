#pragma once

#include <string>
#include "Core/Materials/MaterialTable.h"

namespace DictionaryInput { struct ConfigurationTrace; }

MaterialTable readMaterialTable(const std::string& path, DictionaryInput::ConfigurationTrace& trace);
