#pragma once

#include <string>

namespace Solvers
{
// Values owned by the linear-system configuration, independent of its reader.
struct DirectLinearSystemInput
{
    std::string backend;
    double relativeTolerance = 0.0;
};
}
