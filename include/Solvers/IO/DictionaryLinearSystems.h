#pragma once

#include <map>
#include <string>
#include <vector>
#include "Solvers/DirectLinearSystemInput.h"

namespace DictionaryInput { struct ConfigurationTrace; }

namespace Solvers
{
// App supplies actual system identities; this IO adapter owns no physics and
// the numerical kernels never reopen a dictionary. Direct backends have no
// unconsumed CG/preconditioner/iteration options.
std::map<std::string, DirectLinearSystemInput> readDirectLinearSystems(
    const std::string& path, const std::vector<std::string>& requiredSystems,
    DictionaryInput::ConfigurationTrace& trace);
}
