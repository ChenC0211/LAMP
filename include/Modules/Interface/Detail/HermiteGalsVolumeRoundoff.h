#pragma once

#include <limits>

namespace LAMP::Interface::Detail
{
// One relative volume budget shared by inventory assembly and the working
// accuracy of geometry closure. This is a floating roundoff allowance, not
// an estimate/bound on quadrature or discretization error.
inline constexpr long double hermiteInventoryRoundoffBudget=64*std::numeric_limits<double>::epsilon();
}
