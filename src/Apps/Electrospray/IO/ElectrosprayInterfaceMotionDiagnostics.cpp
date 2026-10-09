#include "Electrospray/IO/ElectrosprayInterfaceMotionDiagnostics.h"
#include "Electrospray/Coupling/ElectrosprayInterfaceMotion.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/UniformLeafGrid.h"
#include "Core/Mesh/CartesianLeafLocator.h"
#include "Core/Units/UnitsContext.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <locale>
#include <sstream>
#include <stdexcept>
#include <optional>

namespace
{
void require(bool valid, const char* reason)
{
    if (!valid) throw std::runtime_error(std::string("[ElectrosprayInterfaceMotion] ") + reason);
}
}

void writeElectrosprayInterfaceCurvature(int step, InterfaceMotionStage stage,
    const LevelSetInterfacePanels& panels, const UnitsContext& units)
{
    require(panels.axisymmetric() && !panels.panels().empty(), "curvature observation requires RZ panels");
    const double length = units.scales().L0_m;
    require(std::isfinite(length) && length > 0.0, "invalid curvature observation length scale");
    const char* stageName = nullptr;
    switch (stage)
    {
    case InterfaceMotionStage::Previous: stageName = "输运前"; break;
    case InterfaceMotionStage::Transported: stageName = "输运后、体积修正前"; break;
    case InterfaceMotionStage::VolumeCorrected: stageName = "体积修正后、耦合尚未接受"; break;
    }
    require(stageName != nullptr, "unknown curvature observation stage");
    double minimum = 0.0, maximum = 0.0, polePanel = 0.0;
    int minimumCell = -1, maximumCell = -1;
    for (const auto& panel : panels.panels())
    {
        const auto& frame = panel.frame;
        require(frame.point.x > 0.0, "curvature representative is not at positive radius");
        const double curvature = (frame.curvature + frame.normal.x / frame.point.x) / length;
        require(std::isfinite(curvature), "nonfinite observed RZ curvature");
        if (minimumCell < 0 || curvature < minimum) { minimum = curvature; minimumCell = panel.cellIndex; }
        if (maximumCell < 0 || curvature > maximum) { maximum = curvature; maximumCell = panel.cellIndex; }
        polePanel = curvature;
    }
    const auto& pole = panels.panels().back();
    std::ostringstream message;
    message.imbue(std::locale::classic());
    message << std::scientific << std::setprecision(9)
        << "[电喷雾] 步 " << step << " 液面曲率（" << stageName << "，面元代表点）：范围=["
        << minimum << ", " << maximum << "] 1/m，极值单元=" << minimumCell << '/' << maximumCell
        << "；近轴面元=" << pole.cellIndex << "，曲率=" << polePanel
        << " 1/m，尖端 z=" << pole.end.y * length << " m。\n";
    std::cout << message.str() << std::flush;
    require(static_cast<bool>(std::cout), "failed to print stage curvature observation");
}

void writeElectrosprayInterfaceMotion(const QuadtreeMesh& mesh,
    const ElectrosprayInterfaceState& state, const InterfaceMotionReport& report,
    const UnitsContext& units, const std::string& prefix)
{
    require(state.surface.geometry() && state.surface.panels() && state.surface.pinnedCoordinates().size() == 1
        && state.surface.pinnedVertices().size() == 1
        && state.surface.nodalData().size() == static_cast<std::size_t>(mesh.vertIndexNum())
        && state.surface.activeVertices().size() == static_cast<std::size_t>(mesh.vertIndexNum()), "invalid motion output state");
    InterfaceMotion::validatePinnedToAxisChain(*state.surface.panels(), state.surface.pinnedCoordinates().front());
    const bool native = state.surface.geometry()->usesNativeLeaves();
    std::optional<UniformLeafGrid> grid;
    std::optional<CartesianLeafLocator> locator;
    if (native) locator.emplace(mesh);
    else grid.emplace(mesh,CartesianGridScope::UniformLeaves);
    const double length = units.scales().L0_m;
    require(std::isfinite(length) && length > 0.0, "invalid output length scale");
    std::ofstream nodes(prefix + "_interface_nodes.csv"), panels(prefix + "_interface_panels.csv"),
        summary(prefix + "_interface_summary.csv"), sectors(prefix + "_interface_corner_sectors.csv"),
        quadrature(prefix + "_interface_quadrature.csv");
    require(nodes && panels && summary && sectors && quadrature, "cannot open interface motion output");
    nodes << std::setprecision(17)
        << "vertex_index,vertex_label,actual_r_m,actual_z_m,canonical_r_m,canonical_z_m,active,pinned,phi_m,phi_r,phi_z,phi_rz_per_m\n";
    for (int position = 0; position < mesh.vertIndexNum(); ++position)
    {
        const int label = native ? mesh.vertIndexToLabel[position] : grid->vertices()[position];
        const int k = mesh.verts(label).index();
        const auto point = mesh.verts(label).coord();
        const int i = native ? 0 : position % (grid->nx()+1), j = native ? 0 : position / (grid->nx()+1);
        const double x = native ? point.x : i == grid->nx() ? grid->maxCoord().x : grid->minCoord().x+i*grid->spacing().x;
        const double y = native ? point.y : j == grid->ny() ? grid->maxCoord().y : grid->minCoord().y+j*grid->spacing().y;
        nodes << k << ',' << label << ',' << point.x * length << ',' << point.y * length
            << ',' << x * length << ',' << y * length << ',' << static_cast<int>(state.surface.activeVertices()[k])
            << ',' << (k == state.surface.pinnedVertices().front());
        if (state.surface.activeVertices()[k])
        {
            auto node = state.surface.nodalData()[k];
            if (native && mesh.verts(label).supCell >= 0)
            {
                const auto incident = locator->incident(point);
                bool found = false;
                for (int n = 0; n < incident.count; ++n)
                    if (state.surface.geometry()->nativeCellActive(incident.indices[n]))
                    {
                        const auto jet = state.surface.geometry()->scalarJetInNativeCell(incident.indices[n],point);
                        node = {state.surface.displayValue(k),jet.dx,jet.dy,jet.dxy};
                        found = true; break;
                    }
                require(found,"hanging output node has no actual active scalar leaf");
            }
            nodes << ',' << node.value * length << ',' << node.dx << ',' << node.dy << ',' << node.dxy / length;
        }
        else nodes << ",,,,";
        nodes << '\n';
    }
    panels << std::setprecision(17)
        << "panel,cell_index,begin_r_m,begin_z_m,end_r_m,end_z_m,r_m,z_m,n_r,n_z,kappa_meridional_per_m,kappa_RZ_per_m,area_m2\n";
    // Base-node output alone does not define a sector-enabled scalar. Keep
    // the explicit cell-local extension and all sub-global-ULP curve positions.
    sectors << std::setprecision(17)
        << "cell_index,vertex_index,phase_sign,phi_m,phi_r,phi_z,phi_rz_per_m\n";
    for (const auto& sector : state.surface.cornerSectors())
        sectors << sector.cellIndex << ',' << sector.vertexIndex << ',' << sector.phaseSign << ','
            << sector.dof.value*length << ',' << sector.dof.dx << ',' << sector.dof.dy << ',' << sector.dof.dxy/length << '\n';
    quadrature << std::setprecision(17)
        << "panel,node,r_m,z_m,has_cell_local,cell_i,cell_j,local_r,local_z,n_r,n_z,kappa_meridional_per_m,weight_m,physical_weight_m2\n";
    int index = 0;
    for (const auto& panel : state.surface.panels()->panels())
    {
        const auto& frame = panel.frame;
        require(frame.point.x > 0.0, "free-panel interior lies on or behind the axis");
        panels << index++ << ',' << panel.cellIndex << ',' << panel.begin.x * length << ',' << panel.begin.y * length
            << ',' << panel.end.x * length << ',' << panel.end.y * length << ',' << frame.point.x * length
            << ',' << frame.point.y * length << ',' << frame.normal.x << ',' << frame.normal.y
            << ',' << frame.curvature / length << ',' << (frame.curvature + frame.normal.x / frame.point.x) / length
            << ',' << panel.physicalMeasure * length * length << '\n';
        int q = 0;
        for (const auto& node : panel.quadratureNodes)
        {
            const auto& f = node.frame;
            quadrature << index-1 << ',' << q++ << ',' << f.point.x*length << ',' << f.point.y*length
                << ',' << f.hasCellLocalPoint;
            if (f.hasCellLocalPoint)
                quadrature << ',' << f.cellI << ',' << f.cellJ << ',' << f.cellLocalPoint.x << ',' << f.cellLocalPoint.y;
            else quadrature << ",,,,";
            quadrature << ',' << f.normal.x << ',' << f.normal.y << ',' << f.curvature/length
                << ',' << node.weight*length << ',' << node.physicalWeight*length*length << '\n';
        }
    }
    summary << std::setprecision(17)
        << "dt_code,active_nodes,changed_values,stationary_nodes,tapered_nodes,support_cells,pressure_cells,full_phi_width_m,zero_phi_width_m,max_phi_increment_m,max_trace_ratio,panels,free_area_before_m2,free_area_after_m2,pin_r_m,pin_z_m,pin_normal_before_r,pin_normal_before_z,pin_normal_after_r,pin_normal_after_z,pole_before_r_m,pole_before_z_m,pole_after_r_m,pole_after_z_m\n"
        << report.dt << ',' << report.activeNodes << ',' << report.changedValues << ',' << report.stationaryNodes
        << ',' << report.taperedNodes << ',' << report.supportCells << ',' << report.pressureCells
        << ',' << report.fullWidth * length << ',' << report.zeroWidth * length
        << ',' << report.maxValueIncrement * length << ',' << report.maxTraceRatio << ',' << state.surface.panels()->panels().size()
        << ',' << report.areaBefore * length * length << ',' << report.areaAfter * length * length
        << ',' << state.surface.pinnedCoordinates().front().x * length << ',' << state.surface.pinnedCoordinates().front().y * length
        << ',' << report.pinNormalBefore.x << ',' << report.pinNormalBefore.y
        << ',' << report.pinNormalAfter.x << ',' << report.pinNormalAfter.y
        << ',' << report.poleBefore.x * length << ',' << report.poleBefore.y * length
        << ',' << report.poleAfter.x * length << ',' << report.poleAfter.y * length << '\n';
    nodes.close(); panels.close(); summary.close(); sectors.close(); quadrature.close();
    require(nodes && panels && summary && sectors && quadrature, "incomplete interface motion output");
}
