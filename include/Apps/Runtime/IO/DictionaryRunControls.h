#pragma once

#include "Runtime/RunControls.h"
class UnitsContext;
namespace DictionaryInput { struct Node; struct ConfigurationTrace; }

// The selected App validates allowed top-level/workflow blocks. This reader
// handles only the common Time/Output/Log/Restart blocks of schema 2.
RunControls readRunControls(const DictionaryInput::Node& root,
    const UnitsContext& units, DictionaryInput::ConfigurationTrace& trace,bool allowHydroControls=false);
