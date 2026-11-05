#include "Modules/Hydrodynamics/HydroVofAdvector.h"

#include "Core/Mesh/QuadtreeMesh.h"
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

namespace LAMP::Hydrodynamics
{
Interface::GeometricVofAdvectorSlab prepareHydroVofAdvectorSlab(
    const QuadtreeMesh& mesh,const Interface::GeometricVofPublication& oldVof,
    Interface::GeometricVofTimeLayer nextEnd,
    const HydroNativeValues& beginValues,const HydroNativeValues& endGuessValues)
{
    oldVof.validateMesh(mesh);
    const auto source=oldVof.identity();
    const auto require=[source](bool ok,const char* reason){
        if(!ok)throw std::invalid_argument("[HydroVofAdvector/prepare source="
            +std::to_string(source.value())+"] "+reason);
    };
    const auto& layout=beginValues.layout();
    beginValues.validateFor(layout,mesh);
    endGuessValues.validateFor(layout,mesh);
    const auto& domain=layout.domain();
    const auto& vofDomain=oldVof.domain();
    require(domain.axisymmetric()==vofDomain.axisymmetric(),"Hydro and VOF measures differ");
    require(domain.meshIdentity()==vofDomain.meshIdentity(),"different native mesh binding; explicit transfer required");
    require(static_cast<std::size_t>(layout.pressureCount())==vofDomain.cellIndexByRow().size()
        &&domain.cellIndexByRow()==vofDomain.cellIndexByRow()
        &&domain.cellRowByIndex()==vofDomain.cellRowByIndex(),"Hydro layout does not cover the actual selected VOF leaf domain");
    require(domain.nativeFaceCount()==vofDomain.nativeFaceCount(),"Hydro native face layouts differ");

    std::vector<double> beginNormal(static_cast<std::size_t>(domain.nativeFaceCount()));
    std::vector<double> endNormal(beginNormal.size());
    for(int nativeRow=0;nativeRow<domain.nativeFaceCount();++nativeRow){
        const auto face=vofDomain.face(mesh,nativeRow);
        const int velocity=layout.velocityDof(nativeRow);
        require(face.nativeRow==nativeRow,"invalid VOF native subface row");
        if(face.negativeCell<0&&face.positiveCell<0){
            require(velocity<0,"excluded VOF face cannot hold a Hydro velocity DOF");continue;}
        require(velocity>=0&&layout.nativeFace(velocity)==nativeRow,"invalid native subface velocity mapping");
        const auto axis=face.horizontal?HydroDofLayout::NormalAxis::Second:HydroDofLayout::NormalAxis::First;
        require(layout.normalAxis(velocity)==axis,"Hydro normal axis differs from the canonical native subface direction");
        const double a=beginValues.normalVelocity()[static_cast<std::size_t>(velocity)];
        const double b=endGuessValues.normalVelocity()[static_cast<std::size_t>(velocity)];
        require(std::isfinite(a)&&std::isfinite(b),"nonfinite native endpoint normal mean");
        // Hydro stores +x/r or +y/z means, exactly the VOF slab convention.
        // No outward sign, area multiplication, or reconstructed point value.
        beginNormal[static_cast<std::size_t>(nativeRow)]=a;
        endNormal[static_cast<std::size_t>(nativeRow)]=b;
    }
    return Interface::GeometricVofAdvectorSlab::prepare(mesh,oldVof,nextEnd,beginNormal,endNormal);
}
}