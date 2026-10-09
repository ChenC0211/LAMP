#include "Modules/Hydrodynamics/HydroUniformLoad.h"
#include "Electrohydrodynamics/IO/SharpDielectricDiagnostics.h"

#include <filesystem>
#include <iomanip>
#include <locale>
#include <sstream>
#include <stdexcept>
#include "Modules/Hydrodynamics/HydroInterfaceLoadGeometry.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/CartesianMeasure.h"
#include "Core/IO/OutputLayout.h"
#include "Modules/Electricity/DielectricElectrostatics.h"
#include "Modules/Electricity/DielectricAssemblyDiagnostics.h"

namespace
{
std::ofstream newOutput(const std::filesystem::path& path)
{
    if (std::filesystem::exists(path))
        throw std::runtime_error("refusing to overwrite sharp electric audit: " + path.string());
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot create sharp electric audit: " + path.string());
    out.imbue(std::locale::classic());
    out << std::scientific << std::setprecision(17);
    return out;
}
void finish(std::ofstream& out)
{
    out.close();
    if (!out) throw std::runtime_error("failed to write sharp electric audit");
}

void writeAssembly(const std::filesystem::path& prefix, const DielectricSharpAssemblyDiagnostics& diagnostics,
    bool axisymmetric)
{
    auto faces = newOutput(prefix.string() + "_global_faces.tsv");
    auto coefficients = newOutput(prefix.string() + "_global_coefficients.tsv");
    auto pieces = newOutput(prefix.string() + "_global_pieces.tsv");
    faces << "graph_segment\tcomponent\tedge_index\towner_cell\tneighbor_cell\towner_negative\tneighbor_negative\t"
             "x\ty\toutward_nx\toutward_ny\tface_area\tcenter_distance\tepsilon_face\tlink_negative_fraction\t"
             "face_negative_fraction\tinterface_x\tinterface_y\tnx\tny\tkappa\tkappa_s\tcondition\treproduction\t"
             "correction_constant\tcoefficient_count\tpiece_count\tsample_cell_indices\tsample_negative";
    if (axisymmetric) faces << "\tphysical_face_measure";
    faces << '\n';
    coefficients << "graph_segment\tordinal\tcell_index\tvalue\n";
    pieces << "graph_segment\tordinal\tbegin_x\tbegin_y\tend_x\tend_y\tnegative\n";
    for (const auto& face : diagnostics.correctedFaces)
    {
        faces << face.graphSegmentIndex << '\t' << (face.horizontal ? 'y' : 'x') << '\t' << face.edgeIndex
            << '\t' << face.ownerCellIndex << '\t' << face.neighborCellIndex
            << '\t' << face.ownerNegative << '\t' << face.neighborNegative
            << '\t' << face.faceCenter.x << '\t' << face.faceCenter.y
            << '\t' << face.outwardNormal.x << '\t' << face.outwardNormal.y
            << '\t' << face.faceArea << '\t' << face.centerDistance << '\t' << face.epsilonFace
            << '\t' << face.linkNegativeFraction << '\t' << face.faceNegativeFraction
            << '\t' << face.interfacePoint.x << '\t' << face.interfacePoint.y
            << '\t' << face.interfaceNormal.x << '\t' << face.interfaceNormal.y
            << '\t' << face.curvature << '\t' << face.curvatureTangentialDerivative
            << '\t' << face.condition << '\t' << face.reproduction << '\t' << face.correctionConstant
            << '\t' << face.correction.size() << '\t' << face.facePieces.size() << '\t';
        for (std::size_t k = 0; k < face.sampleCellIndices.size(); ++k)
            faces << (k ? "," : "") << face.sampleCellIndices[k];
        faces << '\t';
        for (std::size_t k = 0; k < face.sampleNegative.size(); ++k)
            faces << (k ? "," : "") << static_cast<int>(face.sampleNegative[k]);
        if (axisymmetric) faces << '\t' << face.physicalFaceMeasure;
        faces << '\n';
        for (std::size_t k = 0; k < face.correction.size(); ++k)
            coefficients << face.graphSegmentIndex << '\t' << k << '\t' << face.correction[k].cellIndex
                << '\t' << face.correction[k].value << '\n';
        for (std::size_t k = 0; k < face.facePieces.size(); ++k)
        {
            const auto& piece = face.facePieces[k];
            pieces << face.graphSegmentIndex << '\t' << k << '\t' << piece.begin.x << '\t' << piece.begin.y
                << '\t' << piece.end.x << '\t' << piece.end.y << '\t' << piece.negative << '\n';
        }
    }
    finish(faces); finish(coefficients); finish(pieces);
}
}

void SharpDielectricDiagnostics::initialize()
{
    summary_ = newOutput(OutputLayout::series() / "electric_sharp_diagnostics.csv");
    summary_ << "step,force_time,dt,used_superlu,corrected_physical_faces,certified_boundary_faces,"
                "zero_neumann_faces,local_trace_fits,active_load_faces,cell_integrated_flux_residual,"
                "max_trace_condition,max_trace_reproduction_error,zero_area_axis_faces\n";
}

void SharpDielectricDiagnostics::writeStep(const QuadtreeMesh& mesh, const VertField<double>& phi,
    const HydroInterfaceLoadGeometry& geometry, const DielectricFields& electric, const Quantum::Hydrodynamics::HydroUniformLoad& load,
    const DielectricConfig& config, const DielectricReport& report,
    int step, double forceTime, double dt, double materialThickness,
    const std::string& outputName, bool writeFields,
    const DielectricSharpAssemblyDiagnostics* assemblyDiagnostics)
{
    geometry.validateSource(mesh, phi);
    electric.meshIdentity.validate(mesh);
    if (electric.traceGeometry != geometry.traceGeometry() ||
        !report.success || !report.usedSharpInterface || !report.usedSuperLU ||
        report.localTraceFits != geometry.faces().size() || electric.interfaceTraces.size() != geometry.faces().size())
        throw std::runtime_error("sharp audit requires the complete actual solve and shared trace set");
    summary_ << step << ',' << forceTime << ',' << dt << ',' << report.usedSuperLU << ','
        << report.correctedPhysicalFaces << ',' << report.certifiedBoundaryFaces << ','
        << report.zeroNeumannFaces << ',' << report.localTraceFits << ',' << geometry.faces().size() << ','
        << report.cellIntegratedFluxResidual << ',' << report.maxTraceCondition << ','
        << report.maxTraceReproductionError << ',' << report.zeroAreaAxisFaces << '\n';
    summary_.flush();
    if (!summary_) throw std::runtime_error("failed to write sharp electric step audit");
    if (!writeFields) return;
    if (!assemblyDiagnostics || assemblyDiagnostics->correctedFaces.size() != report.correctedPhysicalFaces)
        throw std::runtime_error("sharp audit requires every actual corrected physical face");
    std::ostringstream stem;
    stem << outputName << "_electric_sharp_step" << std::setw(4) << std::setfill('0') << step;
    const auto prefix = OutputLayout::fields() / stem.str();
    const auto metadataPrefix = OutputLayout::meta() / stem.str();

    // Exactly the old quantitative analyzer's CELL columns. Metadata below
    // identifies an actual App/shared level-set source, not analyticCircle.
    auto cells = newOutput(prefix.string() + "_cells.tsv");
    cells << "x\ty\tarea\tV\tepsilon\tEx\tEy\tfx\tfy";
    if (config.axisymmetric) cells << "\tphysical_cell_volume";
    cells << '\n';
    for (int i = 0; i < mesh.cellIndexNum(); ++i)
    {
        const auto& cell = mesh.cells(mesh.cellIndexToLabel[i]);
        const auto p = cell.centerCoord();
        cells << p.x << '\t' << p.y << '\t' << cell.area() << '\t' << electric.potential[i]
            << '\t' << electric.permittivity[i] << '\t' << electric.ex[i] << '\t' << electric.ey[i]
            << '\t' << load.forceX.value[i] << '\t' << load.forceY.value[i];
        if (config.axisymmetric) cells << '\t' << CartesianMeasure::cell(cell.area(), p.x, true);
        cells << '\n';
    }
    finish(cells);
    auto vertices = newOutput(prefix.string() + "_phi.tsv");
    vertices << "vertex_index\tx\ty\tphi\n";
    for (int i = 0; i < mesh.vertIndexNum(); ++i)
    {
        const auto p = mesh.verts(mesh.vertIndexToLabel[i]).coord();
        vertices << i << '\t' << p.x << '\t' << p.y << '\t' << electric.levelSet[i] << '\n';
    }
    finish(vertices);
    auto faces = newOutput(prefix.string() + "_faces.tsv");
    faces << "component\tedge_index\tx\ty\tdual_area\tlow_cell\thigh_cell\tH_low\tH_high\tGH\t"
             "interface_x\tinterface_y\tnx\tny\tkappa\tkappa_s\tV_Gamma\tEt\tDn\tdsEt\tdsDn\tdssEt\tdssDn\t"
             "tau\telectric_force\tcondition\treproduction\tsample_cell_indices";
    if (report.usedSurfaceCharge) faces << "\ttransverse_GH\tDn_positive\tsurface_charge_density\ttau_t";
    if (config.axisymmetric) faces << "\tphysical_dual_volume\tsurface_curvature";
    faces << '\n';
    const auto& fraction = geometry.loadIndicator();
    for (std::size_t i = 0; i < geometry.faces().size(); ++i)
    {
        const auto& face = geometry.faces()[i];
        const auto& trace = electric.interfaceTraces[i];
        const auto& frame = face.frame;
        const double force = face.horizontal ? load.faceForceY.value[face.edgeIndex] : load.faceForceX.value[face.edgeIndex];
        faces << (face.horizontal ? 'y' : 'x') << '\t' << face.edgeIndex << '\t' << face.point.x << '\t' << face.point.y
            << '\t' << face.dualArea << '\t' << face.lowCellIndex << '\t' << face.highCellIndex
            << '\t' << fraction[face.lowCellIndex] << '\t' << fraction[face.highCellIndex] << '\t' << face.gradient
            << '\t' << frame.point.x << '\t' << frame.point.y << '\t' << frame.normal.x << '\t' << frame.normal.y
            << '\t' << frame.curvature << '\t' << face.curvatureTangentialDerivative << '\t' << trace.interfacePotential
            << '\t' << trace.tangentialElectricField << '\t' << trace.normalDisplacement
            << '\t' << trace.tangentialFieldDerivative << '\t' << trace.normalDisplacementDerivative
            << '\t' << trace.tangentialFieldSecondDerivative << '\t' << trace.normalDisplacementSecondDerivative
            << '\t' << trace.normalTraction << '\t' << force << '\t' << trace.condition << '\t' << trace.reproduction << '\t';
        for (std::size_t k = 0; k < trace.sampleCellIndices.size(); ++k)
            faces << (k ? "," : "") << trace.sampleCellIndices[k];
        if (report.usedSurfaceCharge)
            faces << '\t' << face.transverseGradient << '\t' << trace.positiveNormalDisplacement
                  << '\t' << trace.surfaceChargeDensity << '\t' << trace.tangentialTraction;
        if (config.axisymmetric) faces << '\t' << face.dualMeasure << '\t' << face.surfaceCurvature;
        faces << '\n';
    }
    finish(faces);
    auto meta = newOutput(metadataPrefix.string() + "_metadata.json");
    meta << "{\"source\":\"ElectrohydrodynamicsApp\",\"geometry\":\"sharedLevelSetHermite\","
            "\"interfaceDiscretization\":\"sharpHermiteP3\",\"model\":\"" << config.model << "\",\"usedSuperLU\":true,"
            "\"step\":" << step << ",\"force_time\":" << forceTime << ",\"dt\":" << dt
        << ",\"nx\":" << mesh.Nx << ",\"ny\":" << mesh.Ny << ",\"min\":[" << mesh.minCoord.x << ',' << mesh.minCoord.y
        << "],\"max\":[" << mesh.maxCoord.x << ',' << mesh.maxCoord.y << "],\"eps_negative\":" << config.permittivityNegative
        << ",\"eps_positive\":" << config.permittivityPositive << ",\"axis\":\"" << config.electrodeAxis
        << "\",\"lowerPotential\":" << config.lowerPotential << ",\"upperPotential\":" << config.upperPotential
        << ",\"materialThickness\":" << materialThickness
        << ",\"load_layer\":\"" << (geometry.options().layer == HydroInterfaceLoadGeometry::Layer::GeometricPL ? "geometricPL" : "regularized")
        << "\",\"load_half_width\":" << geometry.options().halfWidth
        << ",\"cell_force_scope\":\"electric face averages only; excludes capillarity\"";
    if (report.usedSurfaceCharge)
        meta << ",\"surface_charge\":true,\"Dn_side\":\"negative\",\"conductivity_negative\":" << config.conductivityNegative
             << ",\"conductivity_positive\":" << config.conductivityPositive
             << ",\"charge_convection\":" << (config.chargeConvection ? "true" : "false");
    if (config.axisymmetric)
        meta << ",\"coordinateSystem\":\"axisymmetricRZ\",\"x_component\":\"radial\","
                "\"integral_measure\":\"revolved 2*pi*r\","
                "\"totalForceX_scope\":\"integral of radial component, not net 3D Cartesian force\","
                "\"cell_epsilon_scope\":\"meridional PL mixture for diagnostics only\""
             << ",\"radial_axis_symmetry\":" << (geometry.snapshot().radialAxisSymmetry() ? "true" : "false")
             << ",\"zero_area_axis_faces\":" << report.zeroAreaAxisFaces
             << ",\"certified_boundary_faces_scope\":\"strict no-contact, excludes zero-area symmetry-axis faces\"";
    meta << "}\n";
    finish(meta);
    writeAssembly(metadataPrefix, *assemblyDiagnostics, config.axisymmetric);
}
