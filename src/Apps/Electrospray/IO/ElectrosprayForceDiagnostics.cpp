#include "Electrospray/IO/ElectrosprayForceDiagnostics.h"
#include "Runtime/RunControls.h"
#include "Electrospray/Coupling/ElectrosprayForces.h"
#include "Modules/Hydrodynamics/HydroFittedDomain.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Units/UnitsContext.h"
#include "Core/IO/OutputLayout.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>

// Legacy CSV names and SI operation order intentionally remain unchanged.
void writeElectrosprayCapillaryDiagnostics(const QuadtreeMesh& mesh,
    const HydroFittedCapillaryState& state, double referenceDensitySI, double capillaryHalfWidthCells,
    const RunControls& run, const UnitsContext& units)
{
    const double length = units.scales().L0_m;
    const double velocity = length / units.scales().T0_s;
    const double pressure = referenceDensitySI * velocity * velocity;
    const double sigma = state.sigma, halfWidth = state.halfWidth;
    if (run.detailedLog)
        std::cout << std::setprecision(17) << "[Electrospray] Fitted capillary load: sigma_code=" << sigma
        << ", halfWidth_code=" << halfWidth << ", halfWidthCells=" << capillaryHalfWidthCells
        << ", loadedFaces=" << state.faces.size() << ", endpointProjections=" << state.endpointProjections
        << ", C_H_identity_error=" << state.maxGradientIdentityError
        << ", maxComponent_Npm3=" << state.maxAbsComponent * pressure / length
        << ", activeRadialLoad_N=" << state.integralX * pressure * length * length
        << ", activeAxialLoad_N=" << state.integralY * pressure * length * length
        << "; actual fitted pressure gradient/duals, no electric force or contact-line reaction.\n";
    if (run.writeFields)
    {
        const std::string prefix = (OutputLayout::fields() / run.outputPrefix).string();
        std::ofstream faces(prefix + "_capillary_faces_initial.csv"), cells(prefix + "_capillary_indicator_initial.csv");
        if (!faces || !cells) throw std::runtime_error("cannot open fitted capillary initial output");
        faces << std::setprecision(17)
            << "horizontal,edge_index,low_cell,high_cell,r_m,z_m,projected_r_m,projected_z_m,n_r,n_z,kappa_meridional_per_m,kappa_RZ_per_m,G_H_per_m,dual_volume_m3,force_density_Npm3,integrated_component_N,endpoint_projection,kkt_residual_m,levelset_residual_m,stationarity_residual_m\n";
        for (const auto& face : state.faces)
        {
            const double force = face.horizontal ? state.forceY.value[face.edgeIndex] : state.forceX.value[face.edgeIndex];
            faces << face.horizontal << ',' << face.edgeIndex << ',' << face.lowCellIndex << ',' << face.highCellIndex
                << ',' << face.point.x * length << ',' << face.point.y * length
                << ',' << face.frame.point.x * length << ',' << face.frame.point.y * length
                << ',' << face.frame.normal.x << ',' << face.frame.normal.y
                << ',' << face.frame.curvature / length << ',' << face.surfaceCurvature / length
                << ',' << face.gradient / length << ',' << face.dualMeasure * length * length * length
                << ',' << force * pressure / length << ',' << force * face.dualMeasure * pressure * length * length
                << ',' << face.boundaryEndpointProjection << ',' << face.endpointKktResidual * length
                << ',' << face.frame.levelSetResidual * length << ',' << face.frame.stationarityResidual * length << '\n';
        }
        cells << std::setprecision(17) << "cell_index,r_m,z_m,active,H_positive\n";
        for (int i = 0; i < mesh.cellIndexNum(); ++i)
        {
            const auto point = mesh.cells(mesh.cellIndexToLabel[i]).centerCoord();
            const bool active = std::isfinite(state.indicator[i]);
            cells << i << ',' << point.x * length << ',' << point.y * length << ',' << active << ',';
            if (active) cells << state.indicator[i];
            cells << '\n';
        }
        faces.close(); cells.close();
        if (!faces || !cells) throw std::runtime_error("incomplete fitted capillary initial output");
    }
}

void writeElectrosprayDielectricForceDiagnostics(const QuadtreeMesh& mesh,
    const HydroFittedDomain& domain, const ElectrosprayDielectricForces& result,
    const HydroFittedCapillaryState& capillary, double referenceDensitySI,
    const RunControls& run, const UnitsContext& units)
{
    constexpr double eps0 = 8.8541878128e-12;
    const double length = units.scales().L0_m;
    const double velocity = length / units.scales().T0_s;
    const double pressure = referenceDensitySI * velocity * velocity;
    const double fieldScale = units.scales().Phi0_V / length;
    const double electricToHydro = eps0 * fieldScale * fieldScale / pressure;
    const auto& electric = result.electric;
    const auto& load = result.polarization;
    if (run.detailedLog)
        std::cout << std::setprecision(17) << "[Electrospray] Fitted cellHarmonic polarization: loadedFaces=" << load.faces.size()
        << ", electricToHydro=" << electricToHydro << ", actualGaussResidual=" << electric.chargeFreeResidual
        << ", boundaryDisplacementClosure=" << electric.boundaryFluxRelativeClosure
        << ", C_epsilon_identity_error=" << load.maxGradientIdentityError
        << ", maxComponent_Npm3=" << load.maxAbsComponent * pressure / length
        << ", activeRadialLoad_N=" << load.integralX * pressure * length * length
        << ", activeAxialLoad_N=" << load.integralY * pressure * length * length
        << "; q=0 FV polarization, not sharp Maxwell traces or an energy-variation certificate.\n";
    if (run.writeFields)
    {
        const std::string prefix = (OutputLayout::fields() / run.outputPrefix).string();
        std::ofstream faces(prefix + "_polarization_faces_initial.csv"), cells(prefix + "_reconstructed_electric_initial.csv"),
            total(prefix + "_total_face_forces_initial.csv");
        if (!faces || !cells || !total) throw std::runtime_error("cannot open spray electric-force output");
        faces << std::setprecision(17) << "native_row,horizontal,low_cell,high_cell,r_m,z_m,epsilon_low,epsilon_high,D_cartesian_Cpm2,tangent_E_Vpm,G_epsilon_relative_per_m,dual_volume_m3,polarization_Npm3,integrated_component_N\n";
        for (const auto& face : load.faces)
        {
            const auto& native = domain.faces()[face.nativeRow];
            faces << face.nativeRow << ',' << native.horizontal << ',' << face.lowCell << ',' << face.highCell
                << ',' << native.center.x * length << ',' << native.center.y * length
                << ',' << electric.permittivityByCell[face.lowCell] << ',' << electric.permittivityByCell[face.highCell]
                << ',' << electric.displacementByNative[face.nativeRow] * eps0 * fieldScale << ',' << electric.tangentByNative[face.nativeRow] * fieldScale
                << ',' << face.gradient / length << ',' << face.dualVolume * length * length * length
                << ',' << face.forceDensity * pressure / length
                << ',' << face.forceDensity * face.dualVolume * pressure * length * length << '\n';
        }
        cells << std::setprecision(17) << "cell_index,active,r_m,z_m,E_r_Vpm,E_z_Vpm\n";
        for (int k = 0; k < mesh.cellIndexNum(); ++k)
        {
            const auto point = mesh.cells(mesh.cellIndexToLabel[k]).centerCoord();
            const bool active = domain.cellRowByIndex()[k] >= 0;
            cells << k << ',' << active << ',' << point.x * length << ',' << point.y * length << ',';
            if (active) cells << electric.fieldX[k] * fieldScale << ',' << electric.fieldY[k] * fieldScale;
            else cells << ',';
            cells << '\n';
        }
        total << std::setprecision(17) << "native_row,horizontal,velocity_row,r_m,z_m,dual_volume_m3,capillary_Npm3,polarization_Npm3,total_Npm3\n";
        for (const auto& face : domain.faces())
        {
            const int i = face.horizontal ? face.nativeRow - mesh.vEdgeIndexNum() : face.nativeRow;
            const double cap = face.horizontal ? capillary.forceY.value[i] : capillary.forceX.value[i];
            const double elec = face.horizontal ? load.forceY.value[i] : load.forceX.value[i];
            if (cap == 0.0 && elec == 0.0) continue;
            const double sum = face.horizontal ? result.total.y()[i] : result.total.x()[i];
            total << face.nativeRow << ',' << face.horizontal << ',' << face.velocityRow
                << ',' << face.center.x * length << ',' << face.center.y * length
                << ',' << face.dualVolume * length * length * length
                << ',' << cap * pressure / length << ',' << elec * pressure / length << ',' << sum * pressure / length << '\n';
        }
        faces.close(); cells.close(); total.close();
        if (!faces || !cells || !total) throw std::runtime_error("incomplete spray electric-force output");
    }
}

void reportElectrosprayConductingForces(const ElectrosprayConductingForces& result,
    const UnitsContext& units)
{
    constexpr double eps0 = 8.8541878128e-12;
    const double fieldScale = units.scales().Phi0_V / units.scales().L0_m;
    std::ostringstream message;
    message << std::scientific << std::setprecision(6)
        << "[电喷雾] 完整界面受力完成：载荷面=" << result.geometry.faces().size()
        << "，仅横向梯度非零的面=" << result.geometry.transverseOnlyFaces
        << "，非零固壁标量梯度来源=" << result.geometry.nonzeroWallGradientSources
        << "，端点投影=" << result.geometry.endpointProjections << "。\n"
        << "[电喷雾] Maxwell 电应力峰值：法向=" << result.electric.diagnostics.maxNormalTraction
            * eps0 * fieldScale * fieldScale
        << " Pa，切向=" << result.electric.diagnostics.maxTangentialTraction * eps0 * fieldScale * fieldScale
        << " Pa；主梯度与实际压力矩阵的一致性误差=" << result.geometry.maxGradientIdentityError
        << "；新载荷供下一次 NS 使用。\n";
    std::cout << message.str() << std::flush;
}
