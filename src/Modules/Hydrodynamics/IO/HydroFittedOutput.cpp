#include "Modules/Hydrodynamics/IO/HydroFittedOutput.h"
#include "Modules/Hydrodynamics/HydroFittedDomain.h"
#include "Modules/Hydrodynamics/HydroFittedVelocity.h"
#include "Modules/Hydrodynamics/Hydrodynamics.h"
#include "Core/Mesh/QuadtreeMesh.h"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

namespace HydroFittedOutput
{
namespace
{
void require(bool condition, const char* message)
{
    if (!condition) throw std::invalid_argument(std::string("[HydroFittedOutput] ")+message);
}

double dimensional(double value, double scale, bool divide = false)
{
    require(std::isfinite(value), "nonfinite consumed physical field or coordinate");
    const double result = divide ? value/scale : value*scale;
    require(std::isfinite(result) && (value == 0.0 || result != 0.0),
        "SI conversion is nonfinite or underflowed a nonzero value");
    return result;
}
}

Frame cells(const QuadtreeMesh& mesh, const HydroFittedDomain& domain,
    const Hydrodynamics& hydro, const Scales& scales, int step, double timeSI, int strandId)
{
    const auto state = hydro.state();
    for (double scale : {scales.length, scales.time, scales.density,
        scales.viscosity, scales.pressure, scales.velocity})
        require(std::isfinite(scale) && scale > 0.0, "all SI scales must be positive and finite");
    require(step >= 0 && std::isfinite(timeSI) && timeSI >= 0.0 && strandId > 0,
        "invalid step, SI time or time-series strand identifier");
    domain.validateMeshIdentity(mesh);
    require(domain.axisymmetric(), "RZFluxCompatible output requires an axisymmetric fitted domain");
    require(!domain.cellIndexByRow().empty()
        && domain.cellRowByIndex().size() == static_cast<std::size_t>(mesh.cellIndexNum()),
        "empty or inconsistent active-cell layout");
    for (const auto* field : {&state.density, &state.viscosity, &state.physicalPressure, &state.divergence})
        require(field->size() == mesh.cellIndexNum() && field->hasStorage(),
            "cell field storage does not match the native mesh");
    const auto velocity = hydro.fittedVelocity(mesh, domain, HydroVelocityTimeLevel::Current,
        HydroFittedVelocitySampler::Reconstruction::RZFluxCompatible);

    Frame result;
    auto& zone = result.zone;
    zone.title = "active_fluid_step_"+std::to_string(step);
    zone.coordinateNames = {{"r_m", "z_m"}};
    zone.element = TecplotWriter::Element::Quadrilateral;
    zone.solutionTime = timeSI;
    zone.strandId = strandId;
    zone.auxiliary = {{"step", std::to_string(step)}, {"units", "SI"},
        {"time_unit", "s"}, {"geometry", "axisymmetric_RZ_meridional_actual_leaf_vertices"},
        {"domain", "active_fluid_only"},
        {"velocity_reconstruction", "RZFluxCompatible sampled at actual cell centers; not native face values"},
        {"divergence_source", "stored native FV divergence; not a derivative of plotted cell vectors"}};
    for (const char* name : {"cell_index", "mesh_level", "rho_kg_per_m3", "mu_Pa_s", "p_Pa",
        "div_u_per_s", "u_r_reconstructed_m_per_s", "u_z_reconstructed_m_per_s", "speed_reconstructed_m_per_s"})
        zone.fields.push_back({name, TecplotWriter::Location::CellCentered, {}});
    for (auto& field : zone.fields) field.values.reserve(domain.cellIndexByRow().size());
    result.nativeCellIndices.reserve(domain.cellIndexByRow().size());
    zone.connectivity.reserve(domain.cellIndexByRow().size());
    std::vector<int> compactByNativeVertex(mesh.vertIndexNum(), -1);
    for (int nativeCell : domain.cellIndexByRow())
    {
        require(nativeCell >= 0 && nativeCell < mesh.cellIndexNum()
            && domain.cellRowByIndex()[nativeCell] == static_cast<int>(result.nativeCellIndices.size()),
            "active compact/native cell maps disagree");
        const int label = mesh.cellIndexToLabel[nativeCell];
        const auto& cell = mesh.cells(label);
        require(cell.index() == nativeCell, "native cell label/index mismatch");
        for (const auto* field : {&state.density, &state.viscosity, &state.physicalPressure, &state.divergence})
            require(field->label(nativeCell) == label, "stale consumed cell-field labels");
        std::vector<int> corners;
        corners.reserve(4);
        for (int corner = 0; corner < 4; ++corner)
        {
            const auto& vertex = mesh.verts(cell.cellToVert[corner]);
            const int nativeVertex = vertex.index();
            require(nativeVertex >= 0 && nativeVertex < mesh.vertIndexNum()
                && mesh.vertIndexToLabel[nativeVertex] == cell.cellToVert[corner],
                "active cell corner has invalid native vertex identity");
            int& compact = compactByNativeVertex[nativeVertex];
            if (compact < 0)
            {
                require(result.nativeVertexIndices.size() < static_cast<std::size_t>(std::numeric_limits<int>::max()),
                    "compact vertex count exceeds index range");
                compact = static_cast<int>(result.nativeVertexIndices.size());
                result.nativeVertexIndices.push_back(nativeVertex);
                const auto point = vertex.coord();
                zone.points.push_back({dimensional(point.x, scales.length), dimensional(point.y, scales.length)});
            }
            corners.push_back(compact);
        }
        result.nativeCellIndices.push_back(nativeCell);
        zone.connectivity.push_back(std::move(corners)); // Original cellToVert order.
        const auto u = velocity.sample(cell.centerCoord());
        const double ur = dimensional(u.x, scales.velocity), uz = dimensional(u.y, scales.velocity);
        const double speed = std::hypot(ur, uz);
        require(std::isfinite(speed), "nonfinite reconstructed SI speed");
        const double values[] = {static_cast<double>(nativeCell), static_cast<double>(cell.level()),
            dimensional(state.density[nativeCell], scales.density),
            dimensional(state.viscosity[nativeCell], scales.viscosity),
            dimensional(state.physicalPressure[nativeCell], scales.pressure),
            dimensional(state.divergence[nativeCell], scales.time, true), ur, uz, speed};
        for (std::size_t field = 0; field < zone.fields.size(); ++field)
            zone.fields[field].values.push_back(values[field]);
    }
    return result;
}
}
