#pragma once

#include "Core/Geometry/UnitPolynomialRoots.h"

namespace LAMP::Hydrodynamics
{
struct HydroNormalPolynomial;
namespace Detail
{
// Metric density on the canonical unit chart, from the original physical
// normal polynomial. Used by HydroPhaseRate; no velocity/phase/root owner.
std::vector<Geometry::UnitPolynomialTerm> normalRatePolynomial(const HydroNormalPolynomial&);
}
}
