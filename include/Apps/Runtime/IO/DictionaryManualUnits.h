#pragma once
#include <memory>
#include <string>
class UnitsContext;
namespace DictionaryInput { struct ConfigurationTrace; }

// Strict schema-2 SI references. Conversion and derived scales still use the
// common UnitsContext, not a second normalization implementation.
std::shared_ptr<const UnitsContext> readManualSIUnits(const std::string& path,
    DictionaryInput::ConfigurationTrace& trace);
