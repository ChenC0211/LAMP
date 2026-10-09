#pragma once

#include "Core/IO/DictionaryTokenReader.h"

// Preserve the first-batch Hydro reader API and strict token semantics.
namespace HydroInput { using Reader = DictionaryInput::Reader; }
