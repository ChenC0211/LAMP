#pragma once
#include "Modules/Electric/ElectricVofField.h"
#include "Core/Base/Vector2D.h"
#include <array>
#include <memory>
#include <vector>

namespace LAMP::Electric
{
struct ElectricVofMaxwellTensor { double xx=0,yy=0,xy=0,hoop=0; };
struct ElectricVofMaxwellSample
{
    int nativeCell=-1;
    Vector2D point{0,0},displacement{0,0},electricField{0,0};
    ElectricVofMaxwellTensor tensor;
    // Analytic cell RT0 tensor divergence at the represented electric field.
    // Native width/radius reciprocals define its finite coefficients; exact
    // sums/products end in binary64 RNE. RZ includes the hoop and axis limit.
    Vector2D tensorDivergence{0,0};
};
struct ElectricVofMaxwellCell
{
    int nativeCell=-1;
    Vector2D lower{0,0},upper{0,0};
    double epsilon=0;
    // Canonical +x/r or +y/z D integrals, order left,right,bottom,top.
    // These are whole cell-side sums, not individual fine-face point traces.
    std::array<long double,4> sideDisplacement{},sidePhysicalMeasure{};
};
struct ElectricVofMaxwellReport
{
    std::size_t cells=0,axisCells=0;
    double maximumRelativeSideMeasureDefect=0;
    bool wholeCellSidesOnly=true;
};
// Readonly coeval Maxwell stress from ONE actual electric field. No Hydro,
// force deposition, extra solve, charge update, clock relabel or commit.
// XY: Cartesian RT0 D. RZ: r*Dr is affine in r^2, Dz affine in z.
// Both body and exterior-traction sampling use this SAME reconstruction.
class ElectricVofMaxwellStress
{
public:
    static ElectricVofMaxwellStress prepare(const QuadtreeMesh&,const ElectricVofField&);
    const ElectricVofField& field() const;
    LAMP::TimeLayer timeLayer() const;
    bool axisymmetric() const;
    const std::vector<ElectricVofMaxwellCell>& cells() const;
    const ElectricVofMaxwellReport& report() const;
    ElectricVofMaxwellSample sample(int nativeCell,Vector2D point) const;
    void validateFor(const QuadtreeMesh&,const ElectricVofField&,LAMP::TimeLayer) const;
private:
    struct Data;
    std::shared_ptr<const Data> data_;
    explicit ElectricVofMaxwellStress(std::shared_ptr<const Data>);
};
}
