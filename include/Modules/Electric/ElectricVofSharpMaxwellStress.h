#pragma once
#include "Modules/Electric/ElectricVofSharpField.h"
#include "Core/Geometry/PlanarCellCut.h"
#include <memory>

namespace LAMP::Electric
{
struct ElectricVofSharpMaxwellTensor
{
    double xx=0,xy=0,yy=0,hoop=0;
};
struct ElectricVofSharpMaxwellSample
{
    Vector2D electricField;
    ElectricVofSharpMaxwellTensor tensor;
};
struct ElectricVofSharpMaxwellInterface
{
    ElectricVofSharpInterfaceField field;
    std::array<ElectricVofSharpMaxwellTensor,2> tensor;
    double normalJump=0,tangentialJump=0;
    double expectedNormalJump=0,expectedTangentialJump=0;
    double relativeNormalIdentityResidual=0,relativeTangentialJumpResidual=0;
};
// Read-only frozen sharp-field stress. Interface samples use the actual
// equilibrated normal flux and ONE trace tangential derivative from this field.
// It does not supply a Hydro load or a moving-geometry energy certificate.
class ElectricVofSharpMaxwellStress
{
public:
    static ElectricVofSharpMaxwellStress prepare(const QuadtreeMesh&,
        const ElectricVofSharpField&);
    const ElectricVofSharpField& field() const;
    bool axisymmetric() const;
    ElectricVofSharpMaxwellSample sample(int nativeCell,unsigned phase,Vector2D) const;
    ElectricVofSharpMaxwellSample sampleAt(int nativeCell,unsigned phase,
        Geometry::PlanarCellCut::AnchoredPoint) const;
    ElectricVofSharpMaxwellInterface interfaceSample(std::size_t arc,double parameter) const;
    bool sameInstance(const ElectricVofSharpMaxwellStress&) const noexcept;
    void validateFor(const QuadtreeMesh&,const ElectricVofSharpField&) const;
private:
    struct Data;
    std::shared_ptr<const Data> data_;
    explicit ElectricVofSharpMaxwellStress(std::shared_ptr<const Data>);
};
}
