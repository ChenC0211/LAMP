#include "Modules/Hydrodynamics/HydroVofRefinement.h"
#include "Modules/Hydrodynamics/HydroTwoPhaseStep.h"
#include "Modules/Hydrodynamics/Detail/HydroRegularAdvectorGeometry.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
namespace LAMP::Hydrodynamics
{
namespace
{
void require(bool ok,const char* why)
{if(!ok)throw std::invalid_argument(std::string("[HydroVofRefinement] ")+why);}
}
HydroVofRefinementIndicators HydroVofRefinement::sample(
    const QuadtreeMesh& mesh,const HydroTwoPhaseState& accepted)
{
    require(accepted.model()==HydroPhaseGeometry::Model::GeometricVof&&!accepted.prescribedFlow(),
        "indicator requires a genuine closed VOF accepted state");
    const auto time=accepted.timeLayer();
    require(time.position==LAMP::TimeLayer::Position::End,"indicator requires the accepted End clock");
    accepted.validateVofFor(mesh,accepted.phaseGeometry().vofIdentity(),time);
    const auto& publication=accepted.phaseGeometry().vofPublication();
    require(publication.geometryKind()==Interface::GeometricVofGeometryKind::LocalPlic,
        "indicator requires LocalPLIC inventory");
    const auto& layout=accepted.values().layout();
    const auto& domain=publication.domain();
    require(layout.pressureCount()==static_cast<int>(domain.cellIndexByRow().size())&&
        domain.cellRowByIndex().size()==static_cast<std::size_t>(mesh.cellIndexNum())&&
        layout.domain().cellRowByIndex()==domain.cellRowByIndex()&&
        layout.domain().cellIndexByRow()==domain.cellIndexByRow(),"indicator requires the exact selected native maps");
    for(int p=0;p<layout.pressureCount();++p){const int native=domain.cellIndexByRow()[p];
        require(native>=0&&native<mesh.cellIndexNum()&&domain.cellRowByIndex()[native]==p&&
            layout.nativeCell(p)==native&&layout.pressureDof(native)==p,"indicator selected pressure map is not invertible");}
    const auto& regular=accepted.reconstruction().regularGeometry();
    // Existing sparse coefficient application only: no fit/SVD/solve or
    // separate velocity identity is prepared for this passive query.
    const auto coefficients=regular.apply(mesh,accepted.values(),&accepted.boundary());
    HydroVofRefinementIndicators result;
    result.vorticity.resize(mesh.cellIndexNum());result.hAbsVorticity.resize(mesh.cellIndexNum());
    for(int p=0;p<layout.pressureCount();++p){
        const int native=layout.nativeCell(p);const auto& cell=mesh.cells(mesh.cellIndexToLabel[native]);
        const auto jet=regular.cellJet(coefficients,p,cell.centerCoord()).jet;
        const long double omega=static_cast<long double>(jet.derivative.yx)-jet.derivative.xy;
        const auto size=cell.size();const double h=std::min(size.x,size.y);
        const long double indicator=h*std::abs(omega);
        const double value=static_cast<double>(omega),scaled=static_cast<double>(indicator);
        require(std::isfinite(h)&&h>0&&std::isfinite(value)&&std::isfinite(scaled),
            "nonfinite accepted regular vorticity or native spacing");
        result.vorticity[native]=value;result.hAbsVorticity[native]=scaled;
        result.maximumHAbsVorticity=std::max(result.maximumHAbsVorticity,scaled);
    }
    return result;
}
}
