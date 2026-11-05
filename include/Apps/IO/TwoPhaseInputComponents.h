#pragma once
#include "Hydrodynamics/IO/HydrodynamicsInput.h"
#include <array>
#include <string>

namespace AppIO
{
// Common fluid schema and its effective text. Material prefixes retain the
// original source order; the caller may serialize its own additional values.
struct TwoPhaseInputComponents
{
    std::array<std::string,6> effective;
    struct Material {int phaseId=0;std::string prefix;};
    std::array<Material,2> materials;
    std::string materialInterfaces;
};
TwoPhaseInputComponents readTwoPhaseInputComponents(HydrodynamicsInput&,const CaseDocuments&,
    ConstantsUsage&,std::optional<int>,int checkpointDefault=10);
void recordTwoPhaseInputDefaults(const CaseDocuments&,const HydrodynamicsInput&,
    DictionaryInput::ConfigurationTrace&);
double boundedCaseValue(const DictionaryInput::Node&,DictionaryInput::ConfigurationTrace&,
    double lower,bool strictLower,const char* consumer,const char* units="code");
}
