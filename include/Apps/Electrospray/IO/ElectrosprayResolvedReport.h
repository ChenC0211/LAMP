#pragma once
#include "Core/IO/ConfigurationTrace.h"
#include <map>
#include <string>
struct ElectrosprayConfig;
struct ElectrosprayPhysicalConfig;
struct ElectrosprayVersionedInput;
class UnitsContext;
class DictionaryMesh;

// Canonical actual numerical inputs, shared by reports and migration checks.
// Raw requested legacy model names and unused controls stay in provenance.
std::map<std::string, DictionaryInput::ConfigurationValue> electrosprayResolvedValues(
    const ElectrosprayVersionedInput& input);
void writeElectrosprayResolvedReport(const ElectrosprayVersionedInput& input,
    const DictionaryMesh* mesh = nullptr);

// Compatibility entry points, not used by the versioned production report.
std::map<std::string, DictionaryInput::ConfigurationValue> electrosprayResolvedValues(
    const ElectrosprayPhysicalConfig& physical, const ElectrosprayConfig& run, const UnitsContext& units);
void writeElectrosprayResolvedReport(const ElectrosprayPhysicalConfig& physical,
    const ElectrosprayConfig& run, const UnitsContext& units, int schemaVersion,
    const DictionaryInput::ConfigurationTrace* trace = nullptr,
    const DictionaryMesh* mesh = nullptr);
