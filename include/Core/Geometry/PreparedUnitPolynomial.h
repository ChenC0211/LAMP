#pragma once
#include "Core/Geometry/UnitPolynomialRoots.h"
#include <memory>

namespace LAMP::Geometry
{
// Immutable exact expansion of one complete weighted polynomial. Reuse only
// for this same polynomial; no coefficient deletion or root approximation.
// The sole implementation uses the existing strict Core polynomial engine.
class PreparedUnitPolynomial
{
public:
    explicit PreparedUnitPolynomial(const std::vector<UnitPolynomialTerm>&);
    // Copy all retained exact components, then append these terms in their
    // original order. This is not a separately summed polynomial merge.
    PreparedUnitPolynomial(const PreparedUnitPolynomial& base,const std::vector<UnitPolynomialTerm>& additionalTerms);
    UnitPolynomialRoots isolate() const;
    // Same complete exact-root engine, with an explicit absolute resolution
    // policy. Every nonendpoint enclosure has outward width <=64*binary64
    // epsilon; a representative rounded to 0/1 does not make it an endpoint.
    // isolate() retains its original endpoint-relative policy.
    UnitPolynomialRoots isolateNativeAbsolute() const;
    std::array<double,2> enclose(double a=0,double b=1) const;
private:
    struct Data;
    std::shared_ptr<const Data> data_;
};
}
