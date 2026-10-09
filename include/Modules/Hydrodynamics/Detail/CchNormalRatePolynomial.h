#pragma once

#include "Core/Geometry/UnitPolynomialRoots.h"

namespace Quantum::Hydrodynamics
{
struct CchNormalPolynomial;
namespace Detail
{
// Metric density on the canonical unit chart, from the original physical
// normal polynomial. Used by CchPhaseRate; no velocity/phase/root owner.
std::vector<Geometry::UnitPolynomialTerm> normalRatePolynomial(const CchNormalPolynomial&);
}
}
