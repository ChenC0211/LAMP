#include "Electrospray/IO/ElectrosprayElectricDiagnostics.h"
#include "Electrospray/Setup/ElectrosprayElectrostatics.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Units/UnitsContext.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>

namespace
{
using BCType = PolarizationPoissonFVSolverCellCentered::FittedBoundaryType;
constexpr double vacuumPermittivitySI = 8.8541878128e-12;
constexpr double twoPi = 6.283185307179586476925286766559;
}

void reportElectrosprayElectricInitialization(const QuadtreeMesh& mesh,
    const ElectrosprayElectricState& state, double conductorPotentialSI,
    const UnitsContext& units)
{
    double minPhi = conductorPotentialSI, maxPhi = minPhi;
    std::map<int, int> phaseCounts;
    for (int i = 0; i < mesh.cellIndexNum(); ++i)
    {
        const double voltage = units.PotentialCodeToSI(state.dielectric.potential()[i]);
        minPhi = std::min(minPhi, voltage); maxPhi = std::max(maxPhi, voltage);
        ++phaseCounts[state.phaseByCell[i]];
    }
    std::ostringstream message;
    message << std::scientific << std::setprecision(4)
            << "[电喷雾] 初始电势求解完成：SuperLU，未知量=" << state.dielectric.report().N
            << "，相对残差=" << state.dielectric.report().finalResidual << "，电势范围=[" << minPhi << ',' << maxPhi << "] V。\n";
    std::cout << message.str();
    for (const auto& count : phaseCounts)
        std::cout << "[电喷雾] 初始材料相 " << count.first << "：单元数=" << count.second << '\n';
}

void writeElectrosprayElectrostatics(const QuadtreeMesh& mesh,
    const ElectrosprayElectricState& state, const UnitsContext& units,
    const std::string& outputPrefix, bool initialReport, bool verboseDiagnostics)
{
    std::ofstream cells(outputPrefix + "_electric_cells.csv"), faces(outputPrefix + "_electric_faces.csv");
    if (!cells || !faces) throw std::runtime_error("cannot open physical electric output");
    cells << std::setprecision(17) << "cell_index,cell_label,level,r_m,z_m,phase_id,active,epsilon_r,potential_V,net_outward_displacement_C\n";
    faces << std::setprecision(17) << "segment,owner_label,neighbor_label,boundary_id,r_m,z_m,normal_r,normal_z,area_m2,epsilon_r,distance_m,effective_normal_E_Vpm,outward_displacement_C\n";
    const double length = units.scales().L0_m;
    const double voltage = units.scales().Phi0_V;
    const double chargeScale = vacuumPermittivitySI * voltage * length;
    std::vector<double> balance(mesh.cellIndexNum(), 0.0);
    std::map<int, double> boundaryCharge;
    double energy = 0.0, maxField = 0.0;
    const auto faceValues=state.dielectric.faceValues(mesh);
    for (std::size_t j = 0; j < state.dielectric.graph().segments().size(); ++j)
    {
        const auto& segment = state.dielectric.graph().segments()[j];
        const auto value = faceValues[j];
        const double charge = value.chargeCode * chargeScale;
        const double field = value.voltageDropCode * voltage / (value.distanceCode * length);
        balance[mesh.cells(segment.K).index()] += charge;
        if (state.dielectric.boundaries()[j].type != BCType::Interior) boundaryCharge[state.boundaryIdBySegment[j]] += charge;
        if (state.dielectric.boundaries()[j].type != BCType::Interior || segment.K < segment.L)
            energy += 0.5 * charge * value.voltageDropCode * voltage;
        maxField = std::max(maxField, std::abs(field));
        faces << j << ',' << segment.K << ',' << segment.L << ',' << state.boundaryIdBySegment[j]
              << ',' << segment.faceCenter.x * length << ',' << segment.faceCenter.y * length
              << ',' << segment.n.x << ',' << segment.n.y << ',' << twoPi * segment.faceCenter.x * segment.Af * length * length
              << ',' << state.dielectric.epsilonBySegment()[j] << ',' << value.distanceCode * length << ',' << field << ',' << charge << '\n';
    }
    double maxCellBalance = 0.0;
    for (int i = 0; i < mesh.cellIndexNum(); ++i)
    {
        const auto& cell = mesh.cells(mesh.cellIndexToLabel[i]);
        const auto range = state.dielectric.graph().segmentRangeForCellIndex(i);
        maxCellBalance = std::max(maxCellBalance, std::abs(balance[i]));
        cells << i << ',' << cell.label() << ',' << cell.level() << ',' << cell.centerCoord().x * length
              << ',' << cell.centerCoord().y * length << ',' << state.phaseByCell[i] << ',' << (range.first != range.second)
              << ',' << state.dielectric.epsilonByCell()[i] << ',' << state.dielectric.potential()[i] * voltage << ',' << balance[i] << '\n';
    }
    cells.close(); faces.close();
    if (!cells || !faces) throw std::runtime_error("incomplete physical electric output");
    double netBoundaryCharge = 0.0, totalAbsoluteBoundaryCharge = 0.0;
    for (const auto& entry : boundaryCharge)
    {
        netBoundaryCharge += entry.second;
        totalAbsoluteBoundaryCharge += std::abs(entry.second);
        if (verboseDiagnostics)
            std::cout << "[Electrospray] Boundary " << entry.first << " outward displacement_C=" << entry.second << '\n';
    }
    if (verboseDiagnostics)
        std::cout << "[Electrospray] Electric diagnostics: discreteEnergy_J=" << energy
              << ", maxEffectiveFaceE_Vpm=" << maxField << ", maxCellFluxBalance_C=" << maxCellBalance
              << ", boundaryFluxRelativeClosure=" << (totalAbsoluteBoundaryCharge > 0.0 ? std::abs(netBoundaryCharge) / totalAbsoluteBoundaryCharge : 0.0)
              << (initialReport ? ". t=0 cellHarmonic only; conductivity/charge/NS have not advanced.\n"
                  : ". Updated moving-interface q=0 cellHarmonic state; conducting charge is not advanced.\n");
    else if (initialReport)
    {
        std::ostringstream message;
        message << std::scientific << std::setprecision(4)
            << "[电喷雾] 初始静电场完成：最大面电场=" << maxField << " V/m，电场能=" << energy
            << " J；当前为无自由电荷的介电模型。\n";
        std::cout << message.str();
    }
}
