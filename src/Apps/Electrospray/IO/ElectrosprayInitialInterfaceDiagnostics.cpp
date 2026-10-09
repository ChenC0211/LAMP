#include "Electrospray/IO/ElectrosprayInitialInterfaceDiagnostics.h"
#include "Modules/Interface/Compatibility/InterfaceCaseConfig.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/UniformLeafGrid.h"
#include "Core/Units/UnitsContext.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <type_traits>

namespace
{
void require(bool valid, const char* message)
{
    if (!valid) throw std::invalid_argument(std::string("[ElectrosprayInterface] ") + message);
}
}

void reportElectrosprayInterfaceInitialization(const QuadtreeMesh& mesh,
    const ElectrosprayInterfaceState& state, const FittedInterfaceCaseConfig& interfaceConfig)
{
    require(state.surface.pinnedCoordinates().size() == interfaceConfig.pinnedContactPointIds.size()
        && state.surface.pinnedVertices().size() == state.surface.pinnedCoordinates().size()
        && state.surface.initialContactJets().size() == state.surface.pinnedCoordinates().size(), "initial contact report size mismatch");
    for (std::size_t i = 0; i < state.surface.pinnedCoordinates().size(); ++i)
    {
        const int pointId = interfaceConfig.pinnedContactPointIds[i];
        const auto point = state.surface.pinnedCoordinates()[i];
        const int label = mesh.vertIndexToLabel[state.surface.pinnedVertices()[i]];
        const auto actual = mesh.verts(label).coord();
        std::cout << std::setprecision(17) << "[Electrospray] Pinned Point " << pointId
            << " canonical=(" << point.x << ',' << point.y << "), actualVertex=(" << actual.x << ',' << actual.y
            << "), vertexLabel=" << label << "; mesh coordinates unchanged.\n";
    }
    for (std::size_t i = 0; i < state.surface.pinnedCoordinates().size(); ++i)
    {
        const auto point = state.surface.pinnedCoordinates()[i];
        const auto& jet = state.surface.initialContactJets()[i];
        std::cout << std::setprecision(17) << "[Electrospray] Initial one-sided contact jet at (" << point.x << ',' << point.y
            << "): n=(" << jet.normal.x << ',' << jet.normal.y << "), kappaRZ=" << jet.axisymmetricCurvature
            << ", admittedSamples=" << jet.admittedSamples << "; no analytic curvature or dynamic angle imposed.\n";
    }
    const auto activeCount = std::count(state.surface.activeVertices().begin(), state.surface.activeVertices().end(), 1);
    std::cout << "[电喷雾] 自由液面初始化：叶单元=" << mesh.cellIndexNum() << "，有效流体顶点=" << activeCount
        << ", freeLines=" << state.freeLineIds.size() << ", fixedLines=" << state.fixedLineIds.size()
        << ", pinnedContacts=" << state.surface.pinnedVertices().size()
        << ". Fixed solid/composite LS unchanged; active-fluid scalar only, not a whole-domain smooth extension.\n";
}

void writeElectrosprayInitialInterface(const QuadtreeMesh& mesh, const ElectrosprayInterfaceState& state,
    const UnitsContext& units, const std::string& prefix)
{
    require(state.surface.vertexCount() == mesh.vertIndexNum() && state.surface.activeVertices().size() == static_cast<std::size_t>(mesh.vertIndexNum()),
        "free-interface output size mismatch");
    std::ofstream file(prefix + "_free_interface_initial.csv");
    if (!file) throw std::runtime_error("cannot open free-interface initial output");
    file << std::setprecision(17) << "vertex_index,vertex_label,r_m,z_m,active_fluid_vertex,pinned_contact,free_phi_m,pinned_r_m,pinned_z_m\n";
    const double length = units.scales().L0_m;
    for (int i = 0; i < mesh.vertIndexNum(); ++i)
    {
        const int label = mesh.vertIndexToLabel[i];
        const auto point = mesh.verts(label).coord();
        const auto contact = std::find(state.surface.pinnedVertices().begin(), state.surface.pinnedVertices().end(), i);
        const bool pinned = contact != state.surface.pinnedVertices().end();
        file << i << ',' << label << ',' << point.x * length << ',' << point.y * length << ','
            << static_cast<int>(state.surface.activeVertices()[i]) << ',' << pinned << ',';
        if (state.surface.activeVertices()[i]) file << state.surface.displayValue(i) * length;
        if (pinned)
        {
            const auto pin = state.surface.pinnedCoordinates()[static_cast<std::size_t>(contact - state.surface.pinnedVertices().begin())];
            file << ',' << pin.x * length << ',' << pin.y * length;
        }
        else file << ",,";
        file << '\n'; // Inactive is missing, never a plotted zero interface.
    }
    file.close();
    if (!file) throw std::runtime_error("incomplete free-interface initial output");
    require(state.surface.initialContactJets().size() == state.surface.pinnedCoordinates().size()
        && state.surface.pinnedVertices().size() == state.surface.pinnedCoordinates().size(), "contact jet output size mismatch");
    std::ofstream jets(prefix + "_pinned_initial_jets.csv");
    if (!jets) throw std::runtime_error("cannot open pinned initial jets output");
    jets << std::setprecision(17) << "vertex_index,pinned_r_m,pinned_z_m,inward_t_r,inward_t_z,n_r,n_z,phi_r,phi_z,phi_rr_per_m,phi_rz_per_m,phi_zz_per_m,kappa_meridional_per_m,kappa_RZ_per_m,admitted_samples\n";
    for (std::size_t i = 0; i < state.surface.initialContactJets().size(); ++i)
    {
        const auto& jet = state.surface.initialContactJets()[i];
        jets << state.surface.pinnedVertices()[i] << ',' << state.surface.pinnedCoordinates()[i].x * length << ',' << state.surface.pinnedCoordinates()[i].y * length
            << ',' << jet.inwardTangent.x << ',' << jet.inwardTangent.y << ',' << jet.normal.x << ',' << jet.normal.y
            << ',' << jet.dx << ',' << jet.dy << ',' << jet.dxx / length << ',' << jet.dxy / length << ',' << jet.dyy / length
            << ',' << jet.meridionalCurvature / length << ',' << jet.axisymmetricCurvature / length << ',' << jet.admittedSamples << '\n';
    }
    jets.close();
    if (!jets) throw std::runtime_error("incomplete pinned initial jets output");
}


InterfaceInitialGeometryObserver makeElectrosprayInitialGeometryObserver(const QuadtreeMesh& mesh,
    const ElectrosprayInterfaceState& state, const UnitsContext& units, const std::string& outputPrefix)
{
    InterfaceInitialGeometryObserver observer;
    observer.initialTangencyAdjusted = [](int vertexIndex, double tangential, double normal)
    {
        std::cout << std::setprecision(17) << "[Electrospray] Initial curve-tangent compatibility: vertex="
            << vertexIndex << ", rawTangentialDerivative=" << tangential
            << ", retainedNormalDerivative=" << normal << "; t=0 only, raw Q4 report unchanged.\n";
    };
    observer.patchPrepared = [&mesh, &state, &units, outputPrefix](const InterfaceScalarPatch& stored)
    {
        std::visit([&](const auto& patch) {
            if (!outputPrefix.empty())
            {
                std::ofstream file(outputPrefix + "_active_nodal_dofs.csv");
                if (!file) throw std::runtime_error("cannot open active nodal DOF output");
                file << std::setprecision(17) << "vertex_index,active,phi_m,phi_r,phi_z,phi_rz_per_m\n";
                const double length = units.scales().L0_m;
                for (int i = 0; i < mesh.vertIndexNum(); ++i)
                {
                    file << i << ',' << static_cast<int>(state.surface.activeVertices()[i]);
                    if (state.surface.activeVertices()[i])
                    {
                        const auto& node = patch.nodes[i];
                        file << ',' << node.value * length << ',' << node.dx << ',' << node.dy << ',' << node.dxy / length;
                    }
                    else file << ",,,,";
                    file << '\n';
                }
                file.close();
                if (!file) throw std::runtime_error("incomplete active nodal DOF output");
                if constexpr (std::is_same_v<std::decay_t<decltype(patch)>, Quantum::Geometry::ActiveQuarticHermitePatch>)
                {
                    // Diagnostic profile only: these additional independent data
                    // are not contained in the nodal file. No duplicate phi grid.
                    std::ofstream data(outputPrefix + "_active_quartic_data.csv");
                    if (!data) throw std::runtime_error("cannot open quartic scalar data output");
                    data << std::setprecision(17) << "kind,index,vertex_begin,vertex_end,phi_m,normal_derivative\n";
                    for (std::size_t i = 0; i < patch.edgeMidpoints.size(); ++i)
                    {
                        const auto& edge = patch.edgeMidpoints[i];
                        if (!std::isfinite(edge.value)) continue;
                        data << "edge," << i << ',' << edge.vertices[0] << ',' << edge.vertices[1]
                            << ',' << edge.value*length << ',' << edge.normalDerivative << '\n';
                    }
                    for (std::size_t c = 0; c < patch.cellCenters.size(); ++c)
                        if (patch.activeCells[c]) data << "center," << c << ",,," << patch.cellCenters[c]*length << ",\n";
                    data.close();
                    if (!data) throw std::runtime_error("incomplete quartic scalar data output");
                }
            }
        }, stored);
    };
    observer.scalarConstructed = [](const LevelSetHermiteGeometry::Snapshot& geometry)
    {
        if (geometry.usesNativeLeaves())
        {
            std::size_t active = 0;
            for (std::size_t c = 0; c < geometry.nativeCellCount(); ++c)
                active += geometry.nativeCellActive(static_cast<int>(c));
            std::cout << "[电喷雾] 粗细网格双四次液面已构造（角点、边中点和中心共享数据）：活动流体叶单元="
                << active << "；悬挂点使用共享粗细边约束，固体不参与液面构造。\n" << std::flush;
            return;
        }
        std::cout << "[Electrospray] Active Hermite scalar constructed: candidateCells=" << geometry.candidateCellCount()
            << ", seeds=" << geometry.seedCount() << "; inactive solid is not a zero polynomial.\n";
    };
    observer.panelsConstructed = [&units, outputPrefix](const LevelSetInterfacePanels& panels)
    {
        if (!outputPrefix.empty())
        {
            std::ofstream file(outputPrefix + "_free_panels_initial.csv");
            if (!file) throw std::runtime_error("cannot open free-panel output");
            file << std::setprecision(17) << "panel,cell_index,r_m,z_m,n_r,n_z,kappa_meridional_per_m,kappa_RZ_per_m,area_m2\n";
            const double length = units.scales().L0_m;
            int id = 0;
            for (const auto& panel : panels.panels())
            {
                const auto& frame = panel.frame;
                require(frame.point.x > 0.0, "free-panel interior must stay off the axis");
                file << id++ << ',' << panel.cellIndex << ',' << frame.point.x * length << ',' << frame.point.y * length
                    << ',' << frame.normal.x << ',' << frame.normal.y << ',' << frame.curvature / length
                    << ',' << (frame.curvature + frame.normal.x / frame.point.x) / length
                    << ',' << panel.physicalMeasure * length * length << '\n';
            }
            file.close();
            if (!file) throw std::runtime_error("incomplete free-panel output");
        }
        std::cout << "[Electrospray] Active free cap: panels=" << panels.panels().size()
            << ", pinned/axis endpoints retained; no force or motion step is claimed.\n";
    };
    return observer;
}
