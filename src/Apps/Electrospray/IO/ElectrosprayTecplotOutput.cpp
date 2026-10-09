#include "Electrospray/IO/ElectrosprayTecplotOutput.h"
#include "Electrospray/Setup/ElectrosprayWorkflow.h"
#include "Core/Materials/MaterialTable.h"
#include "Electrospray/Setup/ElectrosprayElectrostatics.h"
#include "Electrospray/Setup/ElectrosprayInterface.h"
#include "Modules/Hydrodynamics/IO/HydroFittedOutput.h"
#include "Modules/Hydrodynamics/HydroFittedDomain.h"
#include "Modules/Electricity/FittedDielectricPolarization.h"
#include "Modules/Electricity/FittedSharpDielectricElectrostatics.h"
#include "Modules/Electricity/SurfaceChargeSource.h"
#include "Modules/Electricity/SurfaceChargeField.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/CartesianMeasure.h"
#include "Core/Units/UnitsContext.h"
#include "Core/IO/OutputLayout.h"
#include "Core/IO/TecplotWriter.h"
#include "Electrospray/Coupling/ElectrosprayForces.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <locale>
#include <map>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace
{
using TecplotField = TecplotWriter::Field;
using Location = TecplotWriter::Location;
void require(bool ok,const char* message)
{ if(!ok) throw std::runtime_error(std::string("[ElectrosprayTecplotOutput] ")+message); }
TecplotField cellField(std::string name) { return {std::move(name),Location::CellCentered,{}}; }
void modelMetadata(TecplotWriter::Zone& zone,std::string_view mode)
{
    zone.auxiliary["PhysicsMode"]=mode;
    zone.auxiliary["Coordinates"]="axisymmetric r-z meridian; coordinates in m";
    zone.auxiliary["TimeUnit"]="s";
    zone.auxiliary["ChargeModel"]="q=0 dielectric approximation; conductive charge NOT advanced";
    zone.auxiliary["UnavailableFields"]="surface free charge, conduction current, emission current";
}
}

ElectrosprayTecplotOutput::ElectrosprayTecplotOutput(std::string name, bool diagnosticFields)
    :name_(std::move(name)),diagnosticFields_(diagnosticFields)
{
    require(!name_.empty() && name_.find_first_not_of("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-.")
        ==std::string::npos && name_!="." && name_!="..","output name must be a plain file stem");
    OutputLayout::fields();
    series_.open(OutputLayout::meta()/(name_+"_series.csv"));
    require(static_cast<bool>(series_),"cannot open accepted-frame index");
    series_<<"step,time_s,flow_file,interface_file,force_file,physics_mode,charge_advanced\n";
    series_.flush();
    require(static_cast<bool>(series_),"cannot write accepted-frame index header");
}

void ElectrosprayTecplotOutput::write(const QuadtreeMesh& mesh,const HydroFittedDomain& domain,
    const Hydrodynamics& hydro,const ElectrosprayElectricState& electric,
    const PhaseMaterial& liquid, double surfaceTensionSI,const UnitsContext& units,
    const ElectrosprayInterfaceState* interfaceState,const ElectrosprayDielectricForces* polarization,
    const VEdgeField<double>* totalForceX,const HEdgeField<double>* totalForceY,
    int step,double timeCode,ElectrosprayProblem problem)
{
    const auto physicsMode = electrosprayProblemLabel(problem);
    require(units.mode()==UnitsMode::kSI && step>=0 && std::isfinite(timeCode) && timeCode>=0.0,
        "SI state and finite nonnegative output time required");
    require((totalForceX==nullptr)==(totalForceY==nullptr),"both total-force components are required");
    electric.dielectric.validateMesh(mesh);
    require(electric.dielectric.report().success,"dielectric output requires a solved electric state");
    if(electric.dielectric.potential().source())
        require(electric.dielectric.time()==timeCode && interfaceState
            && electric.dielectric.potential().source()==interfaceState->surface.geometry(),
            "dielectric output requires the matching field/interface time layer");
    HydroFittedOutput::Scales scales;
    scales.length=units.scales().L0_m; scales.time=units.scales().T0_s;
    scales.velocity=scales.length/scales.time; scales.density=liquid.densitySI;
    scales.pressure=scales.density*scales.velocity*scales.velocity;
    scales.viscosity=scales.density*scales.length*scales.velocity;
    const double voltage=units.scales().Phi0_V,fieldScale=voltage/scales.length;
    const double timeSI=timeCode*scales.time,forceScale=scales.pressure/scales.length;
    std::ostringstream stepName;
    stepName<<name_<<"_step"<<std::setw(6)<<std::setfill('0')<<step;
    const std::string base=stepName.str(),flowName=base+"_flow.dat";
    std::string interfaceName,forceName;
    auto flow=HydroFittedOutput::cells(mesh,domain,hydro,scales,step,timeSI,1);
    modelMetadata(flow.zone,physicsMode);
    auto phase=cellField("phase_id"),permittivity=cellField("epsilon_relative"),potential=cellField("potential_V");
    auto ex=cellField("E_r_reconstructed_V_per_m"),ey=cellField("E_z_reconstructed_V_per_m");
    for(int cell:flow.nativeCellIndices)
    {
        phase.values.push_back(electric.phaseByCell.at(cell));
        permittivity.values.push_back(electric.dielectric.epsilonByCell().at(cell));
        potential.values.push_back(electric.dielectric.potential()[cell]*voltage);
        if(polarization)
        {
            ex.values.push_back(polarization->electric.fieldX[cell]*fieldScale);
            ey.values.push_back(polarization->electric.fieldY[cell]*fieldScale);
        }
    }
    for(auto* field:{&phase,&permittivity,&potential}) flow.zone.fields.push_back(std::move(*field));
    if(polarization)
    {
        flow.zone.fields.push_back(std::move(ex));flow.zone.fields.push_back(std::move(ey));
        flow.zone.auxiliary["ElectricFieldMeaning"]="FV displacement reconstruction divided by cell epsilon; not sharp two-sided traces";
    }
    if(interfaceState)
    {
        TecplotField phi{"phi_node_m",Location::Nodal,{}};
        for(int vertex:flow.nativeVertexIndices)
        {
            require(interfaceState->surface.activeVertices().at(vertex)!=0,"fluid output includes an inactive scalar vertex");
            phi.values.push_back((interfaceState->surface.geometry()->usesNativeLeaves()
                ? interfaceState->surface.displayValue(vertex)
                : interfaceState->surface.nodalData().at(vertex).value)*scales.length);
        }
        flow.zone.fields.push_back(std::move(phi));
        flow.zone.auxiliary["InterfaceDisplay"]="Use separate interface zone; a linear nodal phi contour does not retain Hermite derivatives/corner sectors";
    }
    TecplotWriter::write(OutputLayout::fields()/flowName,flow.zone);

    if(interfaceState)
    {
        require(static_cast<bool>(interfaceState->surface.panels()),"missing retained interface panels");
        TecplotWriter::Zone curve;
        curve.title="Free surface";curve.coordinateNames={{"r_m","z_m"}};
        curve.element=TecplotWriter::Element::Line;curve.solutionTime=timeSI;curve.strandId=2;
        modelMetadata(curve,physicsMode);
        curve.auxiliary["GeometryMeaning"]="Polyline through original Hermite panel endpoints and GL8 nodes; display only, not a replacement solver interface";
        curve.auxiliary["FieldLocation"]="Normal, curvature and sigma*kappa are original panel midpoint values, constant on its display segments";
        auto panelId=cellField("panel_id"),normalR=cellField("normal_r"),normalZ=cellField("normal_z");
        auto curvature=cellField("kappa_RZ_per_m"),capillaryPressure=cellField("sigma_kappa_Pa");
        int id=0;
        for(const auto& panel:interfaceState->surface.panels()->panels())
        {
            const int first=static_cast<int>(curve.points.size());
            curve.points.push_back(panel.begin*scales.length);
            for(const auto& q:panel.quadratureNodes) curve.points.push_back(q.frame.point*scales.length);
            curve.points.push_back(panel.end*scales.length);
            require(panel.frame.point.x>0.0,"RZ panel midpoint must be off-axis");
            const double k=(panel.frame.curvature+panel.frame.normal.x/panel.frame.point.x)/scales.length;
            for(int node=first+1;node<static_cast<int>(curve.points.size());++node)
            {
                curve.connectivity.push_back({node-1,node});
                panelId.values.push_back(id);normalR.values.push_back(panel.frame.normal.x);normalZ.values.push_back(panel.frame.normal.y);
                curvature.values.push_back(k);capillaryPressure.values.push_back(surfaceTensionSI*k);
            }
            ++id;
        }
        for(auto* field:{&panelId,&normalR,&normalZ,&curvature,&capillaryPressure}) curve.fields.push_back(std::move(*field));
        interfaceName=base+"_interface.dat";
        TecplotWriter::write(OutputLayout::fields()/interfaceName,curve);
    }
    if(totalForceX)
    {
        TecplotWriter::Zone force;
        force.title="Native face force components";force.coordinateNames={{"r_m","z_m"}};
        force.element=TecplotWriter::Element::Line;force.solutionTime=timeSI;force.strandId=3;
        modelMetadata(force,physicsMode);
        force.auxiliary["ForceLocation"]="Native staggered component per volume; horizontal face=+z, vertical face=+r";
        force.auxiliary["ForceTime"]="Rebuilt at this output time for the next NS step in movingDielectric mode";
        if(polarization)
            force.auxiliary["CapillaryComponentMeaning"]="Derived as dimensional total minus polarization; may include cancellation rounding";
        auto native=cellField("native_face"),nr=cellField("component_r"),nz=cellField("component_z");
        auto total=cellField("total_force_component_N_per_m3"),ep=cellField("polarization_force_component_N_per_m3");
        auto cap=cellField("capillary_force_component_N_per_m3");
        for(const auto& face:domain.faces())
        {
            if(face.negativeCell<0 && face.positiveCell<0) continue;
            const auto& edge=face.horizontal?mesh.hEdges(face.edgeLabel):mesh.vEdges(face.edgeLabel);
            const int row=face.horizontal?face.nativeRow-mesh.vEdgeIndexNum():face.nativeRow;
            const double f=(face.horizontal?totalForceY->value[row]:totalForceX->value[row])*forceScale;
            const double e=polarization?(face.horizontal?polarization->polarization.forceY.value[row]:polarization->polarization.forceX.value[row])*forceScale:0.0;
            const int first=static_cast<int>(force.points.size());
            force.points.push_back(mesh.verts(edge.edgeToVert[0]).coord()*scales.length);
            force.points.push_back(face.center*scales.length);
            force.points.push_back(mesh.verts(edge.edgeToVert[1]).coord()*scales.length);
            for(int half=0;half<2;++half)
            {
                force.connectivity.push_back({first+half,first+half+1});
                native.values.push_back(face.nativeRow);nr.values.push_back(face.horizontal?0.0:1.0);nz.values.push_back(face.horizontal?1.0:0.0);
                total.values.push_back(f);
                if(polarization){ep.values.push_back(e);cap.values.push_back(f-e);}
            }
        }
        for(auto* field:{&native,&nr,&nz,&total}) force.fields.push_back(std::move(*field));
        if(polarization){force.fields.push_back(std::move(ep));force.fields.push_back(std::move(cap));}
        forceName=base+"_forces.dat";
        TecplotWriter::write(OutputLayout::fields()/forceName,force);
    }
    // A failed field write leaves no completed-frame entry. Existing CSV
    // conservation diagnostics remain independent and are never rewritten.
    // File references are relative to output/, not the meta/ manifest folder.
    const auto reference=[](const std::string& file) { return file.empty()?std::string{}:"fields/"+file; };
    series_<<std::setprecision(17)<<step<<','<<timeSI<<','<<reference(flowName)<<','<<reference(interfaceName)
        <<','<<reference(forceName)<<','<<physicsMode<<",false\n";
    series_.flush();require(static_cast<bool>(series_),"failed to register a completed Tecplot frame");
}

namespace
{
// Match the SI normalization used by the App's conducting solve, not the
// unrelated legacy Base.h reference scales.
constexpr double conductingVacuumPermittivitySI = 8.8541878128e-12;

double conductingFinite(long double value)
{
    const double result = static_cast<double>(value);
    require(std::isfinite(value) && std::isfinite(result) && (value == 0.0L || result != 0.0),
        "conducting initial output has a nonfinite or underflowed SI value");
    return result;
}

double conductingSI(double value, double scale)
{
    require(std::isfinite(value), "conducting initial output consumed a nonfinite value");
    return conductingFinite(static_cast<long double>(value)*scale);
}

void conductingAgreement(long double a, long double b, long double scale, const char* message)
{
    // The same 4096eps roundoff-level accounting gate as the phase-current
    // consumer. This checks exported data identity, not a new physical test.
    require(std::isfinite(a) && std::isfinite(b) && std::isfinite(scale) && scale >= 0.0L
        && std::abs(a-b) <= 4096.0L*std::numeric_limits<double>::epsilon()*scale, message);
}

void conductingFrameMetadata(TecplotWriter::Zone& zone, int step, double timeCode, double timeSI)
{
    zone.auxiliary["PhysicsMode"] = "movingConducting_sharpHermiteP3";
    zone.auxiliary["Coordinates"] = "axisymmetric r-z meridian; coordinates in m";
    zone.auxiliary["TimeUnit"] = "s";
    zone.auxiliary["charge_advanced"] = step > 0 ? "true" : "false";
    zone.auxiliary["Step"] = std::to_string(step);
    zone.auxiliary["Emission"] = "not implemented; conduction current is not emission current";
    std::ostringstream time;
    time.imbue(std::locale::classic()); time << std::setprecision(17) << timeCode;
    zone.auxiliary["ForcesTimeCode"] = time.str();
    time.str(""); time.clear(); time << std::setprecision(17) << timeSI;
    zone.auxiliary["ForcesTime_s"] = time.str();
    zone.auxiliary["ForceTimeMeaning"] = "same-time forces rebuilt for the NEXT NS step";
}

struct ConductingBoundaryOutput
{
    std::size_t segments = 0;
    long double area = 0.0L, displacement = 0.0L, current = 0.0L; // Code-unit integrals.
};

const char* conductingBoundaryName(PolarizationPoissonFVSolverCellCentered::FittedBoundaryType type)
{
    using Type = PolarizationPoissonFVSolverCellCentered::FittedBoundaryType;
    if (type == Type::Interior) return "interior";
    if (type == Type::Dirichlet) return "dirichlet";
    require(type == Type::HomogeneousNeumann, "invalid conducting output boundary type");
    return "homogeneous_neumann";
}

template<class Write>
void writeConductingCSV(const std::filesystem::path& path, const Write& write)
{
    std::ofstream output;
    output.imbue(std::locale::classic());
    output.exceptions(std::ios::failbit|std::ios::badbit);
    try
    {
        output.open(path, std::ios::out|std::ios::trunc);
        output << std::setprecision(17) << std::defaultfloat;
        write(output);
        output.flush();
        output.close();
    }
    catch (const std::ios_base::failure& error)
    { throw std::runtime_error("[ElectrosprayTecplotOutput] incomplete conducting output "+path.string()+": "+error.what()); }
}
} // namespace

namespace
{
void writeConductingElectricOutput(const QuadtreeMesh& mesh,
    const ElectrosprayElectricState& initialElectricBinding,
    const FittedSharpDielectricResult& result, const SurfaceChargeSource& charge,
    const UnitsContext& units, const std::string& outputName, double timeSI, bool initial,
    int step = 0, double timeCode = 0.0, bool diagnosticFields = true)
{
    using Type = PolarizationPoissonFVSolverCellCentered::FittedBoundaryType;
    require(!outputName.empty() && outputName.find_first_not_of("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-.")
        == std::string::npos && outputName != "." && outputName != "..", "output name must be a plain file stem");
    require(units.mode() == UnitsMode::kSI && charge.axisymmetric,
        "conducting initial output requires explicit SI and axisymmetric surface-charge units");
    require(std::isfinite(timeSI) && timeSI >= 0.0, "invalid conducting output physical time");
    const auto& scales = units.scales();
    for (double scale : {scales.L0_m, scales.T0_s, scales.Phi0_V})
        require(std::isfinite(scale) && scale > 0.0, "conducting initial SI scales must be finite and positive");
    const double length = scales.L0_m, voltage = scales.Phi0_V;
    const double areaScale = conductingFinite(static_cast<long double>(length)*length);
    const double chargeScale = conductingFinite(static_cast<long double>(conductingVacuumPermittivitySI)*voltage*length);
    const double currentScale = conductingFinite(static_cast<long double>(chargeScale)/scales.T0_s);
    require(result.report.success && result.report.usedSharpInterface && result.report.usedSuperLU
        && std::isfinite(result.report.relativeResidual) && result.report.relativeResidual >= 0.0
        && result.report.conductionCurrent, "conducting output requires a successful P3/SuperLU solve and its phase-current report");
    const auto& currentReport = *result.report.conductionCurrent;
    const auto& graph = initialElectricBinding.dielectric.graph();
    const auto& segments = graph.segments();
    const auto& offsets = graph.cellOffsets();
    const int count = mesh.cellIndexNum();
    require(count > 0 && mesh.cellIndexToLabel && mesh.vertIndexToLabel
        && offsets.size() == static_cast<std::size_t>(count)+1 && offsets.front() == 0
        && offsets.back() == segments.size() && !segments.empty(), "invalid initial active graph/index layout");
    require(initialElectricBinding.dielectric.boundaries().size() == segments.size()
        && initialElectricBinding.boundaryIdBySegment.size() == segments.size()
        && result.displacementFlux.size() == segments.size()
        && result.conductionCurrentBySegment.size() == segments.size()
        && result.epsilonBySegment.size() == segments.size() && result.corrections.size() == segments.size()
        && charge.cellIntegratedCharge.size() == static_cast<std::size_t>(count)
        && currentReport.cellIntegratedCurrent.size() == static_cast<std::size_t>(count)
        && result.potential.size() == count, "conducting output arrays do not match the initial active binding");
    std::vector<unsigned char> active(static_cast<std::size_t>(count), 0);
    for (int index = 0; index < count; ++index)
    {
        require(offsets[index] <= offsets[index+1] && offsets[index+1] <= segments.size(), "invalid active graph cell offsets");
        const int label = mesh.cellIndexToLabel[index];
        require(label >= 0 && label < mesh.cellLabelCapacity() && mesh.cells(label).index() == index
            && mesh.cells(label).label() == label && result.potential.label(index) == label,
            "conducting output has stale native cell/field identities");
        active[index] = offsets[index] != offsets[index+1];
        if (!active[index])
        {
            // These are inactive-mask sentinels of integrated source/current,
            // not solid potential values. Never read potential.value[index].
            require(charge.cellIntegratedCharge[index] == 0.0 && currentReport.cellIntegratedCurrent[index] == 0.0,
                "surface charge or conduction report disagrees with the inactive mask");
        }
    }
    FaceSegmentGraph actual;
    actual.build(mesh, [&active](const QuadtreeMesh&, const Cell& cell) { return active[cell.index()] != 0; });
    require(actual.cellOffsets() == offsets && actual.segments().size() == segments.size(),
        "initial graph does not cover the actual active leaf faces");

    std::vector<double> areaSI(segments.size()), displacementSI(segments.size()), currentSI(segments.size());
    std::vector<int> owners(segments.size()), neighbors(segments.size(), -1);
    std::vector<long double> accumulation(count, 0.0L), absoluteCellFaces(count, 0.0L);
    std::map<int, ConductingBoundaryOutput> boundaries;
    long double netCurrent = 0.0L, absoluteBoundaryCurrent = 0.0L, netDisplacement = 0.0L, absoluteBoundaryD = 0.0L;
    for (std::size_t sid = 0; sid < segments.size(); ++sid)
    {
        const auto& face = segments[sid]; const auto& expected = actual.segments()[sid];
        require(face.K == expected.K && face.L == expected.L && face.side == expected.side
            && face.boundary == expected.boundary && face.isHorizontal == expected.isHorizontal
            && face.edgeLabel == expected.edgeLabel && face.Af == expected.Af && face.df == expected.df
            && face.n.x == expected.n.x && face.n.y == expected.n.y
            && face.faceCenter.x == expected.faceCenter.x && face.faceCenter.y == expected.faceCenter.y,
            "conducting output graph has stale or altered actual face identities/metrics");
        owners[sid] = mesh.cells(face.K).index();
        if (face.L >= 0) neighbors[sid] = mesh.cells(face.L).index();
        const bool interior = neighbors[sid] >= 0 && active[neighbors[sid]];
        const auto& boundary = initialElectricBinding.dielectric.boundaries()[sid];
        const int id = initialElectricBinding.boundaryIdBySegment[sid];
        require(active[owners[sid]] && std::isfinite(boundary.value)
            && (interior ? boundary.type == Type::Interior && boundary.value == 0.0 && id == -1
                : id > 0 && (boundary.type == Type::Dirichlet
                    || (boundary.type == Type::HomogeneousNeumann && boundary.value == 0.0))),
            "conducting output boundary ID/type disagrees with actual activity");
        const double area = CartesianMeasure::face(face.Af, face.faceCenter.x, true);
        const double displacement = result.displacementFlux[sid], current = result.conductionCurrentBySegment[sid];
        areaSI[sid] = conductingSI(area, areaScale);
        displacementSI[sid] = conductingSI(displacement, chargeScale);
        currentSI[sid] = conductingSI(current, currentScale);
        // Validate coordinate conversions before opening any output file.
        conductingSI(face.faceCenter.x, length); conductingSI(face.faceCenter.y, length);
        require(std::isfinite(face.n.x) && std::isfinite(face.n.y), "nonfinite actual face normal");
        accumulation[owners[sid]] -= current;
        absoluteCellFaces[owners[sid]] += std::abs(static_cast<long double>(current));
        if (!interior)
        {
            auto& total = boundaries[id]; ++total.segments;
            total.area += area; total.displacement += displacement; total.current += current;
            netCurrent += current; absoluteBoundaryCurrent += std::abs(static_cast<long double>(current));
            netDisplacement += displacement; absoluteBoundaryD += std::abs(static_cast<long double>(displacement));
        }
    }
    require(!boundaries.empty(), "conducting initial output has no actual boundaries");
    conductingAgreement(netCurrent, currentReport.netBoundaryCurrent, absoluteBoundaryCurrent,
        "oriented boundary current sum disagrees with the same solve's net current");
    conductingAgreement(netDisplacement, result.report.netBoundaryFlux, absoluteBoundaryD,
        "oriented boundary displacement sum disagrees with the same solve's net displacement");

    TecplotWriter::Zone zone;
    zone.title = "Conducting P3 initial electric state"; zone.coordinateNames = {{"r_m", "z_m"}};
    zone.element = TecplotWriter::Element::Quadrilateral; zone.solutionTime = timeSI; zone.strandId = 4;
    zone.auxiliary = {{"PhysicsMode", "conducting_initial_sharpHermiteP3"}, {"TimeUnit", "s"},
        {"Coordinates", "actual active native leaf vertices; axisymmetric r-z meridian in m"},
        {"FieldLocation", "all non-coordinate fields are cell centered; no interpolated nodal potential"},
        {"ChargeMeaning", "supplied integral(q dA) per cell in C; NOT pointwise surface-charge density"},
        {"CurrentMeaning", "same-solve conduction accumulation=-sum(outward oriented J) per cell in A; NOT current density"},
        {"CurrentAvailable", "conductivity-weighted phase current from this P3 solve is available"},
        {"NotAdvanced", "t=0 only; no charge or interface time advancement; no Maxwell load or emission computed"},
        {"CurrentSign", "outward current positive; current into active domain negative; cell accumulation positive"}};
    for (const char* name : {"cell_index", "potential_V", "cell_surface_charge_C", "conduction_accumulation_A"})
        zone.fields.push_back(cellField(name));
    std::vector<int> compactVertex(mesh.vertIndexNum(), -1);
    bool zeroIntegratedCharge = true;
    long double totalCharge = 0.0L, totalAccumulation = 0.0L, absoluteAccumulation = 0.0L;
    for (int index = 0; index < count; ++index)
    {
        if (!active[index]) continue;
        const auto& cell = mesh.cells(mesh.cellIndexToLabel[index]);
        conductingAgreement(accumulation[index], currentReport.cellIntegratedCurrent[index], absoluteCellFaces[index],
            "cell conduction accumulation disagrees with minus the same oriented current sum");
        std::vector<int> corners;
        for (int vertexLabel : cell.cellToVert)
        {
            require(vertexLabel >= 0 && vertexLabel < mesh.vertLabelCapacity(), "invalid active-cell vertex label");
            const auto& vertex = mesh.verts(vertexLabel); const int native = vertex.index();
            require(native >= 0 && native < mesh.vertIndexNum() && mesh.vertIndexToLabel[native] == vertexLabel,
                "conducting active-cell vertex has stale native identity");
            int& compact = compactVertex[native];
            if (compact < 0)
            {
                require(zone.points.size() < static_cast<std::size_t>(std::numeric_limits<int>::max()), "too many output vertices");
                compact = static_cast<int>(zone.points.size());
                const auto point = vertex.coord();
                zone.points.push_back({conductingSI(point.x,length), conductingSI(point.y,length)});
            }
            corners.push_back(compact);
        }
        zone.connectivity.push_back(std::move(corners)); // Preserve actual cellToVert orientation.
        const double q = charge.cellIntegratedCharge[index], current = currentReport.cellIntegratedCurrent[index];
        const double values[] = {static_cast<double>(index), conductingSI(result.potential[index],voltage),
            conductingSI(q,chargeScale), conductingSI(current,currentScale)};
        for (std::size_t k = 0; k < zone.fields.size(); ++k) zone.fields[k].values.push_back(values[k]);
        zeroIntegratedCharge = zeroIntegratedCharge && q == 0.0;
        totalCharge += q; totalAccumulation += current; absoluteAccumulation += std::abs(static_cast<long double>(current));
    }
    conductingAgreement(totalAccumulation, currentReport.totalCellCurrent, absoluteAccumulation,
        "exported cell accumulation sum disagrees with the same solve's total cell current");
    zone.auxiliary["InitialChargeState"] = zeroIntegratedCharge
        ? "all supplied initial cell-integrated charges are zero; no additional pointwise q inference"
        : "nonzero supplied initial cell-integrated charge; the writer does not impose q=0";
    if (!initial)
    {
        zone.title = "Accepted conducting P3 electric state";
        zone.auxiliary["PhysicsMode"] = "movingConducting_sharpHermiteP3";
        zone.auxiliary.erase("NotAdvanced");
        zone.auxiliary.erase("InitialChargeState");
        zone.auxiliary["State"] = "accepted same-time electric and surface-charge state; no emission model";
        conductingFrameMetadata(zone, step, timeCode, timeSI);
    }
    long double groupedCurrent = 0.0L, groupedDisplacement = 0.0L;
    for (const auto& entry : boundaries)
    {
        conductingSI(conductingFinite(entry.second.area),areaScale);
        conductingSI(conductingFinite(entry.second.displacement),chargeScale);
        conductingSI(conductingFinite(entry.second.current),currentScale);
        groupedCurrent += entry.second.current; groupedDisplacement += entry.second.displacement;
    }
    conductingAgreement(groupedCurrent, currentReport.netBoundaryCurrent, absoluteBoundaryCurrent,
        "boundary-ID current totals disagree with the same solve's net current");
    conductingAgreement(groupedDisplacement, result.report.netBoundaryFlux, absoluteBoundaryD,
        "boundary-ID displacement totals disagree with the same solve's net displacement");
    const double totalChargeSI = conductingSI(conductingFinite(totalCharge),chargeScale);
    const double totalAccumulationSI = conductingSI(conductingFinite(totalAccumulation),currentScale);
    const double netCurrentSI = conductingSI(conductingFinite(groupedCurrent),currentScale);
    const double currentClosureSI = conductingSI(conductingFinite(totalAccumulation+groupedCurrent),currentScale);

    const std::string& stem = outputName;
    const auto electricPath = OutputLayout::fields()/(stem+"_electric.dat");
    const auto facePath = OutputLayout::fields()/(stem+"_faces.csv");
    const auto boundaryPath = OutputLayout::series()/(stem+"_boundaries.dat");
    require(!std::filesystem::exists(boundaryPath), "refusing to overwrite conducting boundary output");
    // All same-solve identity and D/J accounting checks above remain active.
    // The physical flow zone already owns these cell electric quantities;
    // diagnostic profile alone writes their duplicate mesh and native faces.
    if (diagnosticFields)
    {
        for (const auto& path : {electricPath, facePath})
            require(!std::filesystem::exists(path), "refusing to overwrite conducting diagnostic output");
        TecplotWriter::write(electricPath,zone);
        writeConductingCSV(facePath,[&](std::ostream& output)
        {
            output << "time_s,segment,owner_label,owner_cell_index,neighbor_label,neighbor_cell_index,boundary_id,boundary_type,"
                "r_m,z_m,normal_r,normal_z,area_m2,outward_D_C,outward_J_A\n";
            for (std::size_t sid = 0; sid < segments.size(); ++sid)
            {
                const auto& face = segments[sid];
                output << timeSI << ',' << sid << ',' << face.K << ',' << owners[sid] << ',' << face.L << ',' << neighbors[sid]
                    << ',' << initialElectricBinding.boundaryIdBySegment[sid] << ','
                    << conductingBoundaryName(initialElectricBinding.dielectric.boundaries()[sid].type)
                    << ',' << conductingSI(face.faceCenter.x,length) << ',' << conductingSI(face.faceCenter.y,length)
                    << ',' << face.n.x << ',' << face.n.y << ',' << areaSI[sid] << ',' << displacementSI[sid] << ',' << currentSI[sid] << '\n';
            }
        });
    }
    TecplotWriter::History boundaryTable(boundaryPath, "Boundary groups at one saved physical time",
        {"time_s","boundary_id","segment_count","area_m2","outward_D_C","outward_J_A"},
        {{"CurrentMeaning", "positive outward conduction current, not emission current"},
         {"Abscissa", "boundary_id; each file is one saved-time boundary distribution, not a time history"}});
    for (const auto& entry : boundaries)
        boundaryTable.append(TecplotWriter::row(timeSI, entry.first, entry.second.segments,
            conductingSI(conductingFinite(entry.second.area),areaScale),
            conductingSI(conductingFinite(entry.second.displacement),chargeScale),
            conductingSI(conductingFinite(entry.second.current),currentScale)));
    if (!initial) return; // The accepted-frame caller supplies its own concise status line.
    std::ostringstream message;
    message.imbue(std::locale::classic());
    message << std::scientific << std::setprecision(6)
        << "[电喷雾导电初态] 已输出 t=0 的 P3 电势、单元积分电荷和导电电流；活动单元=" << zone.connectivity.size()
        << "，去重节点=" << zone.points.size() << "。\n"
        << "[电喷雾导电初态] 电流换算 I0=" << currentScale << " A；外向为正，流入计算域为负，单元电荷累积为正。\n";
    for (const auto& entry : boundaries)
        message << "[电喷雾导电初态] 边界 " << entry.first << "：面积="
            << conductingSI(conductingFinite(entry.second.area),areaScale) << " m^2，外向电位移积分="
            << conductingSI(conductingFinite(entry.second.displacement),chargeScale) << " C，外向导电电流="
            << conductingSI(conductingFinite(entry.second.current),currentScale) << " A。\n";
    message << "[电喷雾导电初态] 初始表面电荷总积分=" << totalChargeSI << " C"
        << (zeroIntegratedCharge ? "（各活动单元积分均为零）" : "（保留输入的非零初始电荷）")
        << "；总边界外向电流=" << netCurrentSI << " A，单元总累积率=" << totalAccumulationSI
        << " A，电流收支误差=" << currentClosureSI << " A。\n"
        << "[电喷雾导电初态] 以上为初始化结果，尚未推进电荷、界面或 NS；未计算 Maxwell 力与发射。\n";
    std::cout << message.str() << std::flush;
    require(static_cast<bool>(std::cout), "failed to print conducting initial output summary");
}
} // namespace

void writeElectrosprayConductingInitialOutput(const QuadtreeMesh& mesh,
    const ElectrosprayElectricState& initialElectricBinding,
    const FittedSharpDielectricResult& result, const SurfaceChargeSource& charge,
    const UnitsContext& units, const std::string& outputName)
{
    require(!outputName.empty() && outputName.find_first_not_of("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-.")
        == std::string::npos && outputName != "." && outputName != "..", "output name must be a plain file stem");
    writeConductingElectricOutput(mesh, initialElectricBinding, result, charge, units,
        outputName+"_conducting_initial", 0.0, true);
}

namespace
{
void writeConductingLoads(const QuadtreeMesh& mesh,
    const HydroFittedDomain& domain, const ElectrosprayConductingForces& forces,
    const PhaseMaterial& liquid, const PhaseMaterial& gas, const UnitsContext& units,
    const std::string& outputName, double timeSI, double timeCode, bool initial, int step = 0)
{
    require(!outputName.empty() && outputName.find_first_not_of("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-.")
        == std::string::npos && outputName != "." && outputName != "..", "output name must be a plain file stem");
    require(units.mode() == UnitsMode::kSI && domain.axisymmetric(), "charged load output requires SI/RZ");
    domain.validateMeshIdentity(mesh);
    const double length = units.scales().L0_m, voltage = units.scales().Phi0_V;
    const double velocity = length / units.scales().T0_s;
    const double pressure = liquid.densitySI * velocity * velocity;
    const double fieldScale = voltage / length;
    const double tractionScale = conductingVacuumPermittivitySI * fieldScale * fieldScale;
    const double qScale = conductingVacuumPermittivitySI * fieldScale;
    for (double scale : {length, voltage, velocity, pressure, fieldScale, tractionScale, qScale})
        require(std::isfinite(scale) && scale > 0.0, "invalid charged load output SI scale");
    const auto& faces = forces.geometry.faces();
    const auto& traces = forces.electric.interfaceTraces;
    require(!faces.empty() && traces.size() == faces.size()
        && forces.electric.diagnostics.samplePoints == faces.size(), "incomplete charged load-face trace output");
    const auto nativeSizes = [&](const auto& x, const auto& y) {
        return x.value && y.value && x.size() == mesh.vEdgeIndexNum() && y.size() == mesh.hEdgeIndexNum();
    };
    const auto totalX = forces.total.x();
    const auto totalY = forces.total.y();
    require(totalX.hasStorage() && totalY.hasStorage()
        && totalX.size() == mesh.vEdgeIndexNum() && totalY.size() == mesh.hEdgeIndexNum()
        && nativeSizes(forces.capillary.forceX,forces.capillary.forceY)
        && nativeSizes(forces.maxwell.forceX,forces.maxwell.forceY), "invalid charged load native storage");
    TecplotWriter::Zone zone;
    zone.title = "Conducting initial full Maxwell and capillary native-face loads";
    zone.coordinateNames = {"r_m", "z_m"};
    zone.element = TecplotWriter::Element::Line;
    zone.solutionTime = timeSI;
    zone.strandId = 5;
    zone.auxiliary["State"] = "t=0; no accepted conducting time step; velocity/interface/charge not advanced";
    zone.auxiliary["FieldLocation"] = "native staggered face components; projected interface traces are reconstructed from actual P3 potential";
    zone.auxiliary["WallClosure"] = "explicit ActiveHermiteWallGradient; no wall reaction/contact-line force in sum";
    zone.auxiliary["Curvature"] = "same cell-local Hermite frame for Maxwell and capillarity; kappa_s is same-cell D3";
    zone.auxiliary["Current"] = "no emission or jet current inferred from these loads";
    for (const char* name : {"native_row","horizontal","velocity_row","projected_r_m","projected_z_m",
        "n_r","n_z","kappa_meridional_per_m","kappa_RZ_per_m","kappa_s_per_m2",
        "G_H_per_m","transverse_G_H_per_m","dual_volume_m3","interface_potential_V",
        "tangent_E_Vpm","negative_normal_D_Cpm2","positive_normal_D_Cpm2","surface_charge_Cpm2",
        "normal_Maxwell_Pa","tangential_Maxwell_Pa","capillary_component_Npm3","Maxwell_component_Npm3",
        "total_component_Npm3","integrated_total_component_N","trace_condition","trace_reproduction",
        "endpoint_projection","endpoint_KKT_residual_m"})
        zone.fields.push_back({name, Location::CellCentered, {}});
    if (!initial)
    {
        require(std::isfinite(timeSI) && timeSI >= 0.0 && std::isfinite(timeCode) && timeCode >= 0.0,
            "invalid accepted conducting load time");
        require(std::isfinite(liquid.relativePermittivity) && liquid.relativePermittivity > 0.0
            && std::isfinite(gas.relativePermittivity) && gas.relativePermittivity > 0.0,
            "two-sided interface E requires actual positive phase permittivities");
        zone.title = "Accepted full Maxwell and capillary native-face loads";
        zone.auxiliary["PhysicsMode"] = "movingConducting_sharpHermiteP3";
        zone.auxiliary["State"] = "accepted same-time forces; rebuilt for next NS step, not the old step's applied force";
        zone.auxiliary["TwoSidedElectricField"] = "normal E=D/(epsilon0*epsilon_relative), shared tangent E; negative=liquid, positive=gas; reconstructed at projected interface, not at the native face or cell center";
        zone.auxiliary["ProjectionSourceCoordinate"] = "cell_chart=1: actual coordinate is offset or 1-offset according to from_upper; cell_chart=0 retains legacy local offsets. Do not materialize 1-offset for geometry; rounded projected physical display points may coincide";
        conductingFrameMetadata(zone, step, timeCode, timeSI);
        for (const char* name : {"negative_normal_E_Vpm","positive_normal_E_Vpm",
            "source_cell_i","source_cell_j","cell_chart","local_x_offset","local_x_from_upper",
            "local_y_offset","local_y_from_upper"})
            zone.fields.push_back(cellField(name));
    }
    std::vector<bool> seen(domain.faces().size(), false);
    for (std::size_t i = 0; i < faces.size(); ++i)
    {
        const auto& face = faces[i]; const auto& trace = traces[i];
        const int nativeRow = face.edgeIndex + (face.horizontal ? mesh.vEdgeIndexNum() : 0);
        require(nativeRow >= 0 && nativeRow < static_cast<int>(domain.faces().size()) && !seen[nativeRow],
            "invalid or repeated charged output native face");
        seen[nativeRow] = true;
        const auto& native = domain.faces()[nativeRow];
        require(native.velocityRow >= 0 && native.horizontal == face.horizontal
            && native.center.x == face.point.x && native.center.y == face.point.y
            && native.dualVolume == face.dualMeasure, "charged output face disagrees with the fitted NS dual");
        const auto& edge = face.horizontal ? mesh.hEdges(native.edgeLabel) : mesh.vEdges(native.edgeLabel);
        const int begin = static_cast<int>(zone.points.size());
        for (int label : edge.edgeToVert)
        {
            const auto point = mesh.verts(label).coord();
            zone.points.push_back({conductingSI(point.x,length),conductingSI(point.y,length)});
        }
        zone.connectivity.push_back({begin,begin+1});
        const double capillary = face.horizontal ? forces.capillary.forceY.value[face.edgeIndex]
            : forces.capillary.forceX.value[face.edgeIndex];
        const double electric = face.horizontal ? forces.maxwell.forceY.value[face.edgeIndex]
            : forces.maxwell.forceX.value[face.edgeIndex];
        const double total = face.horizontal ? totalY[face.edgeIndex] : totalX[face.edgeIndex];
        conductingAgreement(total,static_cast<long double>(capillary)+electric*forces.electricToHydro,
            std::abs(capillary)+std::abs(electric*forces.electricToHydro),"output total load differs from capillary plus Maxwell");
        const double values[] = {static_cast<double>(nativeRow), face.horizontal ? 1.0 : 0.0,
            static_cast<double>(native.velocityRow), conductingSI(face.frame.point.x,length),conductingSI(face.frame.point.y,length),
            face.frame.normal.x,face.frame.normal.y,conductingSI(face.frame.curvature,1.0/length),
            conductingSI(face.surfaceCurvature,1.0/length),conductingSI(face.curvatureTangentialDerivative,1.0/(length*length)),
            conductingSI(face.gradient,1.0/length),conductingSI(face.transverseGradient,1.0/length),
            conductingSI(face.dualMeasure,length*length*length),conductingSI(trace.interfacePotential,voltage),
            conductingSI(trace.tangentialElectricField,fieldScale),conductingSI(trace.normalDisplacement,qScale),
            conductingSI(trace.positiveNormalDisplacement,qScale),conductingSI(trace.surfaceChargeDensity,qScale),
            conductingSI(trace.normalTraction,tractionScale),conductingSI(trace.tangentialTraction,tractionScale),
            conductingSI(capillary,pressure/length),conductingSI(electric,tractionScale/length),
            conductingSI(total,pressure/length),conductingFinite(static_cast<long double>(total)*face.dualMeasure*pressure*length*length),
            trace.condition,trace.reproduction,face.boundaryEndpointProjection ? 1.0 : 0.0,
            conductingSI(face.endpointKktResidual,length)};
        constexpr std::size_t valueCount = sizeof(values)/sizeof(values[0]);
        require(valueCount+(initial ? 0 : 9) == zone.fields.size(), "charged load output column mismatch");
        for (std::size_t f = 0; f < valueCount; ++f) zone.fields[f].values.push_back(values[f]);
        if (!initial)
        {
            zone.fields[valueCount].values.push_back(conductingSI(trace.normalDisplacement,
                fieldScale/liquid.relativePermittivity));
            zone.fields[valueCount+1].values.push_back(conductingSI(trace.positiveNormalDisplacement,
                fieldScale/gas.relativePermittivity));
            const auto point = LevelSetHermiteGeometry::cellChartPoint(face.frame);
            const double sourceValues[] = {static_cast<double>(face.frame.cellI),static_cast<double>(face.frame.cellJ),
                face.frame.hasCellChartPoint ? 1.0 : 0.0,point.x.offset,point.x.fromUpper ? 1.0 : 0.0,
                point.y.offset,point.y.fromUpper ? 1.0 : 0.0};
            for (std::size_t f = 0; f < sizeof(sourceValues)/sizeof(sourceValues[0]); ++f)
                zone.fields[valueCount+2+f].values.push_back(sourceValues[f]);
        }
    }
    const auto path = OutputLayout::fields()/(outputName+".dat");
    require(!std::filesystem::exists(path), "refusing to overwrite conducting initial loads");
    TecplotWriter::write(path,zone);
    if (initial) std::cout << "[电喷雾] 初始完整受力及液面两侧电场重构量已输出：" << path.string() << '\n' << std::flush;
}
} // namespace

void writeElectrosprayConductingInitialLoads(const QuadtreeMesh& mesh,
    const HydroFittedDomain& domain, const ElectrosprayConductingForces& forces,
    const PhaseMaterial& liquid, const PhaseMaterial& gas, const UnitsContext& units,
    const std::string& outputName)
{
    require(!outputName.empty() && outputName.find_first_not_of("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-.")
        == std::string::npos && outputName != "." && outputName != "..", "output name must be a plain file stem");
    writeConductingLoads(mesh, domain, forces, liquid, gas, units,
        outputName+"_conducting_initial_loads", 0.0, 0.0, true);
}

void ElectrosprayTecplotOutput::writeConducting(const QuadtreeMesh& mesh,
    const HydroFittedDomain& domain, const Hydrodynamics& hydro,
    const std::vector<int>& phaseByCell, const ElectrosprayElectricState& binding,
    const Quantum::Electricity::FittedChargeState& electrical,
    const PhaseMaterial& liquid, const PhaseMaterial& gas, double surfaceTensionSI, const UnitsContext& units,
    const ElectrosprayInterfaceState& interfaceState, const ElectrosprayConductingForces& forces,
    int step, double timeCode)
{
    electrical.validateMesh(mesh);
    require(electrical.completedStep()==step && electrical.time()==timeCode,"conducting frame time does not match its electrical state");
    const auto& electric=electrical.electric();
    const auto& charge=*electrical.charge();
    require(units.mode() == UnitsMode::kSI && domain.axisymmetric() && step >= 0
        && std::isfinite(timeCode) && timeCode >= 0.0, "conducting frame requires accepted SI/RZ state and finite time");
    domain.validateMeshIdentity(mesh);
    require(interfaceState.surface.geometry() && interfaceState.surface.panels()
        && interfaceState.surface.panels()->usesSource(*interfaceState.surface.geometry())
        && charge.geometry().usesSource(*interfaceState.surface.geometry()) && charge.geometry().axisymmetric(),
        "conducting charge/interface output does not retain the identical scalar");
    const int count = mesh.cellIndexNum();
    const auto& offsets = binding.dielectric.graph().cellOffsets();
    require(phaseByCell.size() == static_cast<std::size_t>(count)
        && domain.cellRowByIndex().size() == static_cast<std::size_t>(count)
        && offsets.size() == static_cast<std::size_t>(count)+1
        && electric.potential.size() == count
        && interfaceState.surface.activeVertices().size() == static_cast<std::size_t>(mesh.vertIndexNum())
        && interfaceState.surface.nodalData().size() == static_cast<std::size_t>(mesh.vertIndexNum()),
        "incompatible conducting frame mesh/state arrays");
    for (int cell = 0; cell < count; ++cell)
        require((domain.cellRowByIndex()[cell] >= 0) == (offsets[cell] != offsets[cell+1]),
            "conducting output electric and NS active masks differ");
    const auto& panels = charge.geometry().panels();
    const auto& interfacePanels = interfaceState.surface.panels()->panels();
    const auto& transport = charge.transport();
    const auto& panelCharge = transport.panelCharge();
    const auto& panelArea = transport.panelMeasure();
    require(!panels.empty() && panels.size() == interfacePanels.size()
        && panels.size() == panelCharge.size() && panels.size() == panelArea.size()
        && charge.geometry().cellToPanel() == interfaceState.surface.panels()->cellToPanel()
        && transport.closed() == charge.geometry().closed(),
        "conducting charge panels and interface output topology disagree");
    for (std::size_t p = 0; p < panels.size(); ++p)
    {
        if (interfaceState.surface.geometry()->usesNativeLeaves())
        {
            const auto& a = charge.geometry().nativeChain().graphs().at(p);
            const auto& b = interfaceState.surface.panels()->nativeChain().graphs().at(p);
            require(a.nativeCellIndex == b.nativeCellIndex && a.nativeCellIndex == panels[p].cellIndex
                && a.nativeCellIndex == interfacePanels[p].cellIndex && a.coordinateIsX == b.coordinateIsX
                && a.begin.side == b.begin.side && a.begin.rootOrdinal == b.begin.rootOrdinal
                && a.end.side == b.end.side && a.end.rootOrdinal == b.end.rootOrdinal
                && panelArea[p] == panels[p].physicalMeasure && panelArea[p] == interfacePanels[p].physicalMeasure
                && std::isfinite(panelArea[p]) && panelArea[p] > 0.0 && std::isfinite(panelCharge[p]),
                "native conducting panel identity, area or charge is stale");
            continue;
        }
        const auto& a = charge.geometry().panelGraph(static_cast<int>(p));
        const auto& b = interfaceState.surface.panels()->panelGraph(static_cast<int>(p));
        require(panels[p].cellIndex == interfacePanels[p].cellIndex
            && a.cellI == b.cellI && a.cellJ == b.cellJ && a.coordinateIsX == b.coordinateIsX
            && a.localBegin.x == b.localBegin.x && a.localBegin.y == b.localBegin.y
            && a.localEnd.x == b.localEnd.x && a.localEnd.y == b.localEnd.y
            && a.beginRoot.side == b.beginRoot.side && a.beginRoot.rootOrdinal == b.beginRoot.rootOrdinal
            && a.endRoot.side == b.endRoot.side && a.endRoot.rootOrdinal == b.endRoot.rootOrdinal
            && panelArea[p] == panels[p].physicalMeasure && panelArea[p] == interfacePanels[p].physicalMeasure
            && std::isfinite(panelArea[p]) && panelArea[p] > 0.0 && std::isfinite(panelCharge[p]),
            "conducting panel identity, physical area or integrated charge is stale");
    }
    // Verify that the handed-off load traces refer to this time's scalar.
    // The writer does not project points, reconstruct E, or recompute forces.
    for (const auto& face : forces.geometry.faces())
    {
        const auto& frame = face.frame;
        require(frame.hasCellLocalPoint != frame.hasCellChartPoint,
            "conducting force output requires exactly one source-local or boundary-chart frame");
        const auto actual = interfaceState.surface.geometry()->usesNativeLeaves()
            ? interfaceState.surface.geometry()->jetAtNativeCellChart(frame.nativeCellIndex,frame.cellChartPoint)
            : frame.hasCellChartPoint
            ? interfaceState.surface.geometry()->jetAtCellChart(frame.cellI, frame.cellJ, frame.cellChartPoint)
            : interfaceState.surface.geometry()->jetAtCellLocal(frame.cellI, frame.cellJ, frame.cellLocalPoint);
        require(actual.point.x == frame.point.x && actual.point.y == frame.point.y
            && actual.normal.x == frame.normal.x && actual.normal.y == frame.normal.y
            && actual.curvature == frame.curvature && actual.levelSetResidual == frame.levelSetResidual,
            "conducting force output uses a stale local scalar frame");
    }

    HydroFittedOutput::Scales scales;
    scales.length = units.scales().L0_m; scales.time = units.scales().T0_s;
    scales.velocity = conductingFinite(static_cast<long double>(scales.length)/scales.time);
    scales.density = liquid.densitySI;
    scales.pressure = conductingFinite(static_cast<long double>(scales.density)*scales.velocity*scales.velocity);
    scales.viscosity = conductingFinite(static_cast<long double>(scales.density)*scales.length*scales.velocity);
    const double voltage = units.scales().Phi0_V;
    for (double scale : {scales.length,scales.time,scales.velocity,scales.density,scales.pressure,scales.viscosity,voltage})
        require(std::isfinite(scale) && scale > 0.0, "invalid conducting frame SI scale");
    const double timeSI = conductingSI(timeCode,scales.time);
    const double areaScale = conductingFinite(static_cast<long double>(scales.length)*scales.length);
    const double chargeScale = conductingFinite(static_cast<long double>(conductingVacuumPermittivitySI)*voltage*scales.length);
    const double currentScale = conductingFinite(static_cast<long double>(chargeScale)/scales.time);
    const auto source = charge.source(); // Integrated cell ledger; no extra pointwise-q evaluation.
    require(source.axisymmetric && source.cellIntegratedCharge.size() == static_cast<std::size_t>(count)
        && electric.report.conductionCurrent
        && electric.report.conductionCurrent->cellIntegratedCurrent.size() == static_cast<std::size_t>(count),
        "conducting frame lacks actual cell charge or same-solve current");
    const auto& current = *electric.report.conductionCurrent;

    std::ostringstream stepName;
    stepName.imbue(std::locale::classic());
    stepName << name_ << "_step" << std::setw(6) << std::setfill('0') << step;
    const std::string base = stepName.str();
    const std::string flowName = base+"_flow.dat", interfaceName = base+"_interface.dat";
    const std::string forceName = base+"_forces.dat";
    const std::string electricName = diagnosticFields_ ? base+"_electric.dat" : flowName;
    const std::string faceName = base+"_faces.csv", boundaryName = base+"_boundaries.dat";
    for (const auto& name : {flowName,interfaceName,forceName})
        require(!std::filesystem::exists(OutputLayout::fields()/name), "refusing to overwrite an existing conducting frame");
    if (diagnosticFields_) for (const auto& name : {electricName,faceName})
        require(!std::filesystem::exists(OutputLayout::fields()/name), "refusing to overwrite an existing conducting diagnostic frame");
    require(!std::filesystem::exists(OutputLayout::series()/boundaryName), "refusing to overwrite conducting boundary output");

    auto flow = HydroFittedOutput::cells(mesh,domain,hydro,scales,step,timeSI,1);
    conductingFrameMetadata(flow.zone,step,timeCode,timeSI);
    flow.zone.auxiliary["ElectricFieldMeaning"] = diagnosticFields_
        ? "no cell E supplied; actual integrated FV D/J are in face CSV; two-sided P3 interface E/D are in forces DAT"
        : "potential and integrated cell charge/current included; two-sided P3 interface E/D in forces DAT; native face CSV omitted in physical profile";
    flow.zone.auxiliary["ChargeMeaning"] = "integral(q dA) per native cell; not a point charge density";
    flow.zone.auxiliary["CurrentMeaning"] = "same-solve conduction accumulation=-sum(outward J) in A; not a current density";
    flow.zone.auxiliary["InterfaceDisplay"] = "separate original Hermite panel/GL8 polyline; nodal phi contours do not retain Hermite derivatives or sectors";
    auto phase = cellField("phase_id"), epsilon = cellField("epsilon_relative");
    auto potential = cellField("potential_V"), cellCharge = cellField("cell_surface_charge_C");
    auto accumulation = cellField("conduction_accumulation_A");
    for (int cell : flow.nativeCellIndices)
    {
        const int id = phaseByCell[cell];
        require(id == liquid.phaseId || id == gas.phaseId, "conducting active output cell has no liquid/gas phase");
        const bool negative = interfaceState.surface.geometry()->negative(mesh.cells(mesh.cellIndexToLabel[cell]).centerCoord());
        require(id == (negative ? liquid.phaseId : gas.phaseId),
            "conducting output phase IDs do not match the current scalar");
        phase.values.push_back(id);
        epsilon.values.push_back(negative ? liquid.relativePermittivity : gas.relativePermittivity);
        potential.values.push_back(conductingSI(electric.potential[cell],voltage));
        cellCharge.values.push_back(conductingSI(source.cellIntegratedCharge[cell],chargeScale));
        accumulation.values.push_back(conductingSI(current.cellIntegratedCurrent[cell],currentScale));
    }
    for (auto* field : {&phase,&epsilon,&potential,&cellCharge,&accumulation}) flow.zone.fields.push_back(std::move(*field));
    TecplotField phi{"phi_node_m",Location::Nodal,{}};
    for (int vertex : flow.nativeVertexIndices)
    {
        require(interfaceState.surface.activeVertices()[vertex] != 0, "conducting flow zone contains an inactive scalar vertex");
        phi.values.push_back(conductingSI(interfaceState.surface.geometry()->usesNativeLeaves()
            ? interfaceState.surface.displayValue(vertex) : interfaceState.surface.nodalData()[vertex].value,scales.length));
    }
    flow.zone.fields.push_back(std::move(phi));

    TecplotWriter::Zone curve;
    curve.title = "Accepted conducting free surface and conserved panel charge";
    curve.coordinateNames = {{"r_m","z_m"}}; curve.element = TecplotWriter::Element::Line;
    curve.solutionTime = timeSI; curve.strandId = 2;
    conductingFrameMetadata(curve,step,timeCode,timeSI);
    curve.auxiliary["GeometryMeaning"] = "polyline through original retained panel endpoints and GL8 nodes; display only, not a solver replacement";
    curve.auxiliary["ChargeMeaning"] = "panel_Q_C is the stored integral(q dA); panel_average_q_Cpm2=Q/panel_area, not a pointwise reconstructed density";
    curve.auxiliary["DisplayDuplication"] = "each panel's Q and area repeat on its display segments; count each panel_id ONCE for sums";
    curve.auxiliary["FieldLocation"] = "normal/kappa/sigma*kappa use the original panel midpoint; Q and area belong to the full original panel";
    for (const char* name : {"panel_id","cell_index","normal_r","normal_z","kappa_RZ_per_m",
        "sigma_kappa_Pa","panel_area_m2","panel_Q_C","panel_average_q_Cpm2"})
        curve.fields.push_back(cellField(name));
    long double totalCharge = 0.0L;
    double maximumDensitySI = 0.0;
    for (std::size_t p = 0; p < panels.size(); ++p)
    {
        const auto& panel = panels[p];
        const int first = static_cast<int>(curve.points.size());
        const auto appendPoint = [&](Vector2D point) {
            curve.points.push_back({conductingSI(point.x,scales.length),conductingSI(point.y,scales.length)});
        };
        appendPoint(panel.begin);
        for (const auto& node : panel.quadratureNodes) appendPoint(node.frame.point);
        appendPoint(panel.end);
        require(panel.frame.point.x > 0.0, "conducting panel midpoint must lie off the axis");
        const double curvature = conductingSI(panel.frame.curvature+panel.frame.normal.x/panel.frame.point.x,1.0/scales.length);
        const double qSI = conductingSI(panelCharge[p],chargeScale);
        const double areaSI = conductingSI(panelArea[p],areaScale);
        const double densitySI = conductingFinite(static_cast<long double>(qSI)/areaSI);
        totalCharge += panelCharge[p]; maximumDensitySI = std::max(maximumDensitySI,std::abs(densitySI));
        const double values[] = {static_cast<double>(p),static_cast<double>(panel.cellIndex),panel.frame.normal.x,panel.frame.normal.y,
            curvature,conductingSI(curvature,surfaceTensionSI),areaSI,qSI,densitySI};
        require(sizeof(values)/sizeof(values[0]) == curve.fields.size(), "conducting panel output column mismatch");
        for (int node = first+1; node < static_cast<int>(curve.points.size()); ++node)
        {
            curve.connectivity.push_back({node-1,node});
            for (std::size_t f = 0; f < curve.fields.size(); ++f) curve.fields[f].values.push_back(values[f]);
        }
    }

    // The reused electric formatter validates actual native identities and
    // the complete same-solve cell/boundary D/J accounting at the old 4096eps
    // gate. A later failure can leave partial files but NEVER a completed-frame
    // entry; the initial wrappers keep their old names, data and t=0 metadata.
    writeConductingElectricOutput(mesh,binding,electric,source,units,base,timeSI,false,step,timeCode,diagnosticFields_);
    writeConductingLoads(mesh,domain,forces,liquid, gas,units,base+"_forces",timeSI,timeCode,false,step);
    TecplotWriter::write(OutputLayout::fields()/flowName,flow.zone);
    TecplotWriter::write(OutputLayout::fields()/interfaceName,curve);

    if (!conductingSeries_.is_open())
    {
        const auto path = OutputLayout::meta()/(name_+"_coupling_series.csv");
        require(!std::filesystem::exists(path), "refusing to overwrite conducting accepted-frame index");
        conductingSeries_.imbue(std::locale::classic());
        conductingSeries_.open(path);
        require(static_cast<bool>(conductingSeries_), "cannot open conducting accepted-frame index");
        conductingSeries_ << "step,time_s,forces_time_code,forces_time_s,flow_file,interface_file,force_file,electric_file,face_flux_file,boundary_file,charge_advanced\n";
        conductingSeries_.flush();
        require(static_cast<bool>(conductingSeries_), "cannot write conducting accepted-frame index header");
    }
    conductingSeries_ << std::setprecision(17) << step << ',' << timeSI << ',' << timeCode << ',' << timeSI
        << ",fields/" << flowName << ",fields/" << interfaceName << ",fields/" << forceName
        << ",fields/" << electricName << ',' << (diagnosticFields_ ? "fields/"+faceName : std::string{})
        << ",series/" << boundaryName
        << ',' << (step > 0 ? "true" : "false") << '\n';
    conductingSeries_.flush();
    require(static_cast<bool>(conductingSeries_), "failed to register the complete conducting frame");
    std::ostringstream frameRecord;
    frameRecord.imbue(std::locale::classic());
    frameRecord << std::setprecision(17) << step << ',' << timeSI << ",fields/" << flowName
        << ",fields/" << interfaceName << ",fields/" << forceName << ",movingConducting_sharpHermiteP3,"
        << (step > 0 ? "true" : "false") << '\n';
    series_ << frameRecord.str();
    series_.flush(); require(static_cast<bool>(series_), "failed to register a completed conducting Tecplot frame");
    std::ostringstream message;
    message.imbue(std::locale::classic());
    message << std::scientific << std::setprecision(6) << "[电喷雾导电输出] 已接受步=" << step
        << "，t=" << timeSI << " s；活动单元=" << flow.nativeCellIndices.size() << "，液面面元=" << panels.size()
        << "，总面电荷=" << conductingSI(conductingFinite(totalCharge),chargeScale)
        << " C，最大|面元平均q|=" << maximumDensitySI << " C/m^2；边界总外向电流="
        << conductingSI(current.netBoundaryCurrent,currentScale) << " A（流入为负）；流场、电势、液面及载荷已输出"
        << (diagnosticFields_ ? "，含原生面诊断。\n" : "；physical 模式不写重复电场和原生面明细。\n");
    std::cout << message.str() << std::flush;
    require(static_cast<bool>(std::cout), "failed to print accepted conducting frame summary");
}
