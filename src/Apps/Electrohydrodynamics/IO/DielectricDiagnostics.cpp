#include "Modules/Hydrodynamics/HydroUniformLoad.h"
#include "Electrohydrodynamics/IO/DielectricDiagnostics.h"
#include "Modules/Electricity/DielectricElectrostatics.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/IO/OutputLayout.h"
#include "Core/IO/TecplotMeshOutput.h"

#include <cmath>
#include <filesystem>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace {
VertField<double> nodalField(const QuadtreeMesh& mesh, const std::vector<double>& values)
{
    if (values.size() != static_cast<std::size_t>(mesh.cellIndexNum()))
        throw std::runtime_error("electric output storage does not match the mesh");
    // Explicit output copy preserves the existing interpolation; no mutable
    // alias to the electric publication escapes through the legacy Field API.
    CellField<double> cell(mesh, 0.0);
    for (int i = 0; i < mesh.cellIndexNum(); ++i) cell.value[i] = values[i];
    return *cell.toVertField(mesh);
}
}

void DielectricDiagnostics::initialize()
{
    const auto diagnosticsPath = OutputLayout::series() / "electric_diagnostics.csv";
    summary_.open(diagnosticsPath);
    if (!summary_) throw std::runtime_error("cannot open " + diagnosticsPath.string());
    summary_ << "step,force_time,dt,electric_relative_residual,boundary_flux_residual,"
                    "max_electric_field,max_electric_force_density,electric_force_x,electric_force_y,"
                    "electric_energy,permittivity_min,permittivity_max,potential_min,potential_max\n"
                 << std::setprecision(17);
}

void DielectricDiagnostics::writeStep(const DielectricReport& report, const Quantum::Hydrodynamics::HydroUniformLoad& load, int step, double forceTime, double dt)
{
    summary_ << step << ',' << forceTime << ',' << dt << ',' << report.relativeResidual
                 << ',' << report.boundaryFluxImbalance << ',' << report.maxField << ',' << load.maxForce
                 << ',' << load.totalForceX << ',' << load.totalForceY << ',' << report.fieldEnergy
                 << ',' << report.minPermittivity << ',' << report.maxPermittivity
                 << ',' << report.minPotential << ',' << report.maxPotential << '\n';
    summary_.flush();
    if (!summary_) throw std::runtime_error("failed to write electric diagnostics");
}

void writeDielectricLoadFields(const QuadtreeMesh& mesh,
    int step, double forceTime, const DielectricFields& electric, const Quantum::Hydrodynamics::HydroUniformLoad& load,
    bool axisymmetric, const std::string& outputName)
{
    electric.meshIdentity.validate(mesh);
    if (electric.levelSet.size() != static_cast<std::size_t>(mesh.vertIndexNum()))
        throw std::runtime_error("electric scalar publication does not match mesh storage");
    VertField<double> levelSet(mesh, 0.0);
    for (int i = 0; i < mesh.vertIndexNum(); ++i) levelSet.value[i] = electric.levelSet[i];
    std::ostringstream name;
    name << outputName << "_electric_load_step" << std::setw(4) << std::setfill('0') << step;
    const std::vector<std::string> names{
        "epsilon", "potential", "Ex", "Ey", "electricForceX", "electricForceY", "levelSet"};
    const std::vector<VertField<double>> fields{
        nodalField(mesh, electric.permittivity), nodalField(mesh, electric.potential),
        nodalField(mesh, electric.ex), nodalField(mesh, electric.ey),
        *load.forceX.toVertField(mesh), *load.forceY.toVertField(mesh), levelSet};
    for (std::size_t f = 0; f < fields.size(); ++f)
    {
        if (fields[f].size() != mesh.vertIndexNum() || !fields[f].value)
            throw std::runtime_error("electric field output size mismatch: " + names[f]);
        for (int i = 0; i < fields[f].size(); ++i)
            if (!std::isfinite(fields[f].value[i]))
                throw std::runtime_error("electric field output is non-finite: " + names[f]);
    }
    const auto path = OutputLayout::fields() / (name.str() + ".dat");
    auto zone = TecplotMeshOutput::nodalZone(mesh, names, fields);
    zone.solutionTime = forceTime;
    zone.strandId = 2;
    zone.auxiliary = {{"CoordinateUnits", "code_length"}, {"FieldUnits", "code_units"},
        {"TimeUnits", "code_time"}, {"TimeRole", "force_time"},
        {"FieldRole", "explicit_electric_load_nodal_diagnostics"},
        {"CoordinateSystem", axisymmetric ? "axisymmetricRZ" : "planar2D"}};
    TecplotWriter::write(path, zone);
    if (!std::filesystem::exists(path) || std::filesystem::file_size(path) == 0)
        throw std::runtime_error("failed to write electric load frame: " + path.string());
}
