#pragma once
#include <initializer_list>
#include <memory>
#include <stdexcept>

namespace LAMP::Hydrodynamics { class HydroKepConvection; }

namespace LAMP::Hydrodynamics::Detail
{
// Derived VOF closed-KEP scalars and explicitly declared finite P0 material
// moments only: no source-inventory/accepted-flux conversion.
// The exceptional path retains the complete signed rational contraction and
// rounds once to native long double, round-to-nearest/ties-to-even. RN zero
// is allowed; final overflow and unsupported floating environments reject.
class HydroFiniteContraction
{
public:
    HydroFiniteContraction();
    ~HydroFiniteContraction();
    // Restore exact 0/1 in the retained Data after the original environment check.
    void reset();
    HydroFiniteContraction(const HydroFiniteContraction&);
    HydroFiniteContraction& operator=(const HydroFiniteContraction&);
    HydroFiniteContraction(HydroFiniteContraction&&) noexcept;
    HydroFiniteContraction& operator=(HydroFiniteContraction&&) noexcept;
    void add(long double value);
    void addProduct(std::initializer_list<long double> factors);
    void addQuotientProduct(std::initializer_list<long double> factors,long double denominator);
    void add(const HydroFiniteContraction& other,int sign=1);
    void addAbsolute(const HydroFiniteContraction& other);
    // Compare exact retained rationals without quantizing either operand.
    int compare(const HydroFiniteContraction& other) const;
    // RN(scale * this / denominator), retaining the entire rational ratio.
    // No intermediate primitive/geometry ratio is rounded.
    long double scaledQuotient(long double scale,const HydroFiniteContraction& denominator) const;
    long double rounded() const;
private:
    friend class LAMP::Hydrodynamics::HydroKepConvection;
    // Only the sealed production request may reuse its stable native qualification.
    void resetQualifiedStorage();
    void addBinary64Product(long double,long double,int binaryScale);
    long double roundedQualified() const;
    long double roundedImpl(bool checkEnvironment) const;
    struct Data;
    std::unique_ptr<Data> data_;
};
// Scoped range detection makes the original ordinary arithmetic an unchanged
// fast path. A marker restarts only an explicitly enabled typed consumer.
struct HydroFiniteRange final:std::exception
{ const char* what() const noexcept override{return "derived finite contraction requires range fallback";} };
class HydroFiniteRangeScope
{
public:
    explicit HydroFiniteRangeScope(bool enabled) noexcept;
    ~HydroFiniteRangeScope();
    HydroFiniteRangeScope(const HydroFiniteRangeScope&)=delete;
    HydroFiniteRangeScope& operator=(const HydroFiniteRangeScope&)=delete;
private:bool previous_;
};
bool hydroFiniteRangeEnabled() noexcept;
void hydroFiniteProductRange(long double a,long double b,long double result);
void hydroFiniteQuotientRange(long double a,long double b,long double result);
}
