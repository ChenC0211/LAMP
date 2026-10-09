#include "Modules/Electricity/DielectricElectrostatics.h"
#include "Modules/Electricity/DielectricPolarizationForce.h"
#include "Core/Geometry/InterfaceGeometry.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <vector>

#include "Core/Mesh/QuadtreeMesh.h"

namespace
{

void require(bool condition, const std::string& message)
{
    if (!condition) throw std::runtime_error(message);
}


std::uint64_t meshSignature(const QuadtreeMesh& mesh)
{
    require(mesh.Nx >= 2 && mesh.Ny >= 2 &&
            static_cast<long long>(mesh.Nx) * mesh.Ny == mesh.cellIndexNum() &&
            (static_cast<long long>(mesh.Nx) + 1) * (static_cast<long long>(mesh.Ny) + 1) == mesh.vertIndexNum() &&
            static_cast<long long>(mesh.Nx) * (static_cast<long long>(mesh.Ny) + 1) == mesh.hEdgeIndexNum() &&
            (static_cast<long long>(mesh.Nx) + 1) * mesh.Ny == mesh.vEdgeIndexNum(),
            "only a fixed uniform Cartesian mesh with at least 2 by 2 cells is supported");
    const double hx = mesh.cellSize.x, hy = mesh.cellSize.y;
    require(std::isfinite(hx) && hx > 0.0 && std::isfinite(hy) && hy > 0.0 &&
            std::isfinite(hx * hy) && hx * hy > 0.0, "invalid mesh spacing/area");
    std::uint64_t signature = 1469598103934665603ull;
    auto mix = [&](std::uint64_t value) { signature = (signature ^ value) * 1099511628211ull; };
    auto mixDouble = [&](double value)
    {
        require(std::isfinite(value), "mesh contains a non-finite coordinate");
        std::uint64_t bits = 0;
        std::memcpy(&bits, &value, sizeof(value));
        mix(bits);
    };
    mix(mesh.Nx); mix(mesh.Ny); mixDouble(hx); mixDouble(hy);
    mixDouble(mesh.minCoord.x); mixDouble(mesh.minCoord.y);
    mixDouble(mesh.maxCoord.x); mixDouble(mesh.maxCoord.y);
    const double tolerance = 128.0 * std::numeric_limits<double>::epsilon();
    for (int i = 0; i < mesh.vertIndexNum(); ++i)
    {
        const int label = mesh.vertIndexToLabel[i];
        require(label >= 0 && label < mesh.vertLabelCapacity() && mesh.verts(label).index() == i,
                "invalid vertex mapping");
        mix(label); mixDouble(mesh.verts(label).coord().x); mixDouble(mesh.verts(label).coord().y);
    }
    for (int i = 0; i < mesh.cellIndexNum(); ++i)
    {
        const int label = mesh.cellIndexToLabel[i];
        require(label >= 0 && label < mesh.cellLabelCapacity() && mesh.cells(label).index() == i,
                "invalid cell mapping");
        const auto& cell = mesh.cells(label);
        require(cell.parent == -1 && cell.level() == 0 &&
                std::abs(cell.size().x - hx) <= tolerance * hx &&
                std::abs(cell.size().y - hy) <= tolerance * hy,
                "refined/coarsened or nonuniform cells are unsupported");
        mix(label); mixDouble(cell.centerCoord().x); mixDouble(cell.centerCoord().y);
        mixDouble(cell.size().x); mixDouble(cell.size().y); mix(static_cast<std::uint64_t>(cell.boundary()));
        for (int corner = 0; corner < 4; ++corner)
        {
            require(cell.son[corner] == -1, "AMR parent cells are unsupported");
            const int vertex = cell.cellToVert[corner];
            require(vertex >= 0 && vertex < mesh.vertLabelCapacity(), "invalid cell vertex label");
            const int index = mesh.verts(vertex).index();
            require(index >= 0 && index < mesh.vertIndexNum() && mesh.vertIndexToLabel[index] == vertex,
                    "invalid cell vertex mapping");
            mix(vertex);
            for (int part = 0; part < 2; ++part) mix(static_cast<std::uint64_t>(cell.neiCell[corner][part]));
        }
        for (int side = 0; side < 2; ++side)
        {
            require(cell.cellToHEdge[side][1] == -1 && cell.cellToVEdge[side][1] == -1,
                    "split faces are unsupported");
            require(cell.cellToHEdge[side][0] >= 0 && cell.cellToHEdge[side][0] < mesh.hEdgeLabelCapacity() &&
                    cell.cellToVEdge[side][0] >= 0 && cell.cellToVEdge[side][0] < mesh.vEdgeLabelCapacity(),
                    "invalid cell-to-edge label");
            const int hIndex = mesh.hEdges(cell.cellToHEdge[side][0]).index();
            const int vIndex = mesh.vEdges(cell.cellToVEdge[side][0]).index();
            require(hIndex >= 0 && hIndex < mesh.hEdgeIndexNum() &&
                    vIndex >= 0 && vIndex < mesh.vEdgeIndexNum() &&
                    mesh.hEdgeIndexToLabel[hIndex] == cell.cellToHEdge[side][0] &&
                    mesh.vEdgeIndexToLabel[vIndex] == cell.cellToVEdge[side][0], "inactive cell-to-edge mapping");
            mix(cell.cellToHEdge[side][0]); mix(cell.cellToVEdge[side][0]);
        }
    }
    // FaceSegmentGraph caches these edge metrics and owner/neighbor labels.
    // Deliberately exclude cut(), cutLength(), isCut and all phase/LS data:
    // this operator solves both fluids and is unaffected by interface motion.
    for (bool horizontal : {true, false})
    {
        const int edgeCount = horizontal ? mesh.hEdgeIndexNum() : mesh.vEdgeIndexNum();
        const int capacity = horizontal ? mesh.hEdgeLabelCapacity() : mesh.vEdgeLabelCapacity();
        const auto* labels = horizontal ? mesh.hEdgeIndexToLabel : mesh.vEdgeIndexToLabel;
        mix(horizontal); mix(edgeCount);
        for (int i = 0; i < edgeCount; ++i)
        {
            const int label = labels[i];
            require(label >= 0 && label < capacity, "invalid edge index-to-label mapping");
            const auto& edge = horizontal ? mesh.hEdges(label) : mesh.vEdges(label);
            require(edge.index() == i && edge.isHorizontal() == horizontal &&
                    edge.parent == -1 && edge.son[0] == -1 && edge.son[1] == -1 &&
                    std::isfinite(edge.length()) && edge.length() > 0.0,
                    "invalid/refined edge topology or metric");
            mix(label); mixDouble(edge.length()); mixDouble(edge.midpointCoord().x); mixDouble(edge.midpointCoord().y);
            mix(static_cast<std::uint64_t>(edge.boundary()));
            for (int vertex : edge.edgeToVert)
            {
                require(vertex >= 0 && vertex < mesh.vertLabelCapacity(), "invalid edge vertex label");
                const int index = mesh.verts(vertex).index();
                require(index >= 0 && index < mesh.vertIndexNum() && mesh.vertIndexToLabel[index] == vertex,
                        "invalid edge vertex mapping");
                mix(vertex);
            }
            for (int neighbor : {edge.masterCell, edge.slaveCell})
            {
                require(neighbor == -1 || (neighbor >= 0 && neighbor < mesh.cellLabelCapacity()), "invalid edge cell label");
                if (neighbor != -1)
                {
                    const int index = mesh.cells(neighbor).index();
                    require(index >= 0 && index < mesh.cellIndexNum() && mesh.cellIndexToLabel[index] == neighbor,
                            "invalid edge cell mapping");
                }
                mix(static_cast<std::uint64_t>(neighbor));
            }
        }
    }
    return signature;
}

double harmonicMean(double a, double b)
{
    const double lower = std::min(a, b), upper = std::max(a, b);
    return lower / (0.5 + 0.5 * (lower / upper));
}

bool isLowerElectrode(Boundary side, bool xAxis)
{
    return side == (xAxis ? Boundary::west : Boundary::south);
}
bool isUpperElectrode(Boundary side, bool xAxis)
{
    return side == (xAxis ? Boundary::east : Boundary::north);
}
} // namespace


DielectricElectrostatics::DielectricElectrostatics(const QuadtreeMesh& mesh) : mesh_(&mesh)
{
    meshSignature_ = meshSignature(mesh);
    graph_.build(mesh); // full domain: neither phase is treated as an excluded solid
    for (int i = 0; i < mesh.cellIndexNum(); ++i)
    {
        const auto range = graph_.segmentRangeForCellIndex(i);
        require(range.second - range.first == 4, "[DielectricElectrostatics] uniform cell must have four face segments");
    }
    for (const auto& face : graph_.segments())
    {
        require(face.K >= 0 && face.K < mesh.cellLabelCapacity() &&
                std::isfinite(face.Af) && face.Af > 0.0 && std::isfinite(face.df) && face.df > 0.0,
                "[DielectricElectrostatics] invalid face metric");
        require((std::abs(face.n.x) == 1.0 && face.n.y == 0.0) ||
                (std::abs(face.n.y) == 1.0 && face.n.x == 0.0),
                "[DielectricElectrostatics] only axis-aligned faces are supported");
        if (face.boundary == Boundary::non)
            require(face.L >= 0 && face.L < mesh.cellLabelCapacity(), "[DielectricElectrostatics] invalid neighbor face");
    }
}

DielectricSolveResult DielectricElectrostatics::solve(const QuadtreeMesh& mesh,
    const VertField<double>& levelSet, double thickness, const DielectricConfig& config,
    const InterfaceGeometry* sharedGeometry,
    DielectricSharpAssemblyDiagnostics* assemblyDiagnostics, const SurfaceChargeSource* surfaceCharge,
    const DielectricConductivity* conductivity, std::shared_ptr<const InterfaceTraceGeometry> traceQueries)
{
    DielectricSolveResult result;
    auto& report = result.report;
    report.stage = "input";
    try
    {
        require(config.model != "leakyDielectric" || surfaceCharge,
            "leakyDielectric requires an explicit surface-charge state from the App");
        validateDielectricConfig(config);
        require(!config.nativeAxisNeckHarmonics,"axis-neck harmonics are available only through the native fitted solver");
        require(!sharedGeometry || sharedGeometry->axisymmetric() == config.axisymmetric,
            "electric coordinate system must match the retained interface geometry");
        if (config.axisymmetric)
        {
            require(mesh.minCoord.x >= 0.0, "RZ electric domain cannot contain negative radius");
            require(mesh.minCoord.x != 0.0 || config.electrodeAxis == "y",
                "the symmetry axis requires homogeneous radial Neumann data, not an electrode");
        }
        require(&mesh == mesh_ && meshSignature(mesh) == meshSignature_,
                "mesh identity/topology/coordinates changed; reconstruct the fixed-grid electrostatics model");
        require(levelSet.value && levelSet.size() == mesh.vertIndexNum(), "level-set size must match mesh vertices");
        require(config.interfaceDiscretization == "sharpHermiteP3" ||
                (std::isfinite(thickness) && thickness > 0.0 && std::isfinite(2.0 * thickness)),
                "interface thickness must be finite, positive and representable when doubled");
        const int count = mesh.cellIndexNum();
        if (config.interfaceDiscretization == "sharpHermiteP3")
        {
            require(sharedGeometry != nullptr, "sharpHermiteP3 requires the caller's retained interface geometry");
            solveSharp(mesh, levelSet, config, *sharedGeometry, report, assemblyDiagnostics, surfaceCharge, conductivity, traceQueries, result.fields);
            return result;
        }
        require(!traceQueries, "physical interface traces require sharpHermiteP3");
        require(assemblyDiagnostics == nullptr, "assembly diagnostics require the sharp interface discretization");
        require(surfaceCharge == nullptr, "prescribed surface charge requires the sharp interface discretization");
        require(conductivity == nullptr, "conduction current requires the sharp interface discretization");
        const bool xAxis = config.electrodeAxis == "x";
        const double voltage = config.upperPotential - config.lowerPotential;

        report.stage = "material";
        std::vector<double> epsilonVertex(static_cast<std::size_t>(mesh.vertIndexNum()));
        for (int i = 0; i < mesh.vertIndexNum(); ++i)
        {
            const double phi = levelSet.value[i];
            require(std::isfinite(phi), "level-set contains a non-finite value");
            const double H = 0.5 * (1.0 + std::tanh(phi / (2.0 * thickness)));
            // Algebraically epsilon_- + H*(epsilon_+ - epsilon_-), evaluated
            // from the closer endpoint to avoid cancellation at large ratios.
            // Equal permittivities are therefore exactly constant.
            epsilonVertex[i] = H <= 0.5
                ? config.permittivityNegative + H * (config.permittivityPositive - config.permittivityNegative)
                : config.permittivityPositive + (1.0 - H) * (config.permittivityNegative - config.permittivityPositive);
            require(std::isfinite(epsilonVertex[i]) && epsilonVertex[i] > 0.0, "non-finite/nonpositive vertex permittivity");
        }

        std::vector<double> epsilonCell(count);
        report.minPermittivity = std::numeric_limits<double>::infinity();
        report.maxPermittivity = 0.0;
        for (int i = 0; i < count; ++i)
        {
            const auto& cell = mesh.cells(mesh.cellIndexToLabel[i]);
            double e[4];
            for (int corner = 0; corner < 4; ++corner)
                e[corner] = epsilonVertex[mesh.verts(cell.cellToVert[corner]).index()];
            epsilonCell[i] = 0.25 * e[0] + 0.25 * e[1] + 0.25 * e[2] + 0.25 * e[3];
            require(std::isfinite(epsilonCell[i]) && epsilonCell[i] > 0.0, "invalid cell permittivity");
            report.minPermittivity = std::min(report.minPermittivity, epsilonCell[i]);
            report.maxPermittivity = std::max(report.maxPermittivity, epsilonCell[i]);
        }

        report.stage = "assembly";
        const auto& faces = graph_.segments();
        std::vector<double> epsilonFace(faces.size()), diagonal(count, 0.0), boundaryRhs(count, 0.0);
        for (std::size_t i = 0; i < faces.size(); ++i)
        {
            const auto& face = faces[i];
            const int owner = mesh.cells(face.K).index();
            require(owner >= 0 && owner < count, "invalid owner index in fixed face graph");
            double coefficient = epsilonCell[owner];
            if (face.boundary == Boundary::non)
            {
                const int neighbor = mesh.cells(face.L).index();
                require(neighbor >= 0 && neighbor < count, "invalid neighbor index in fixed face graph");
                coefficient = harmonicMean(coefficient, epsilonCell[neighbor]);
            }
            epsilonFace[i] = coefficient;
            const double conductance = coefficient * face.Af / face.df;
            require(std::isfinite(coefficient) && coefficient > 0.0 &&
                    std::isfinite(conductance) && conductance > 0.0, "invalid/overflowed face conductance");
            if (face.boundary == Boundary::non || isLowerElectrode(face.boundary, xAxis) || isUpperElectrode(face.boundary, xAxis))
                diagonal[owner] += conductance;
            if (isUpperElectrode(face.boundary, xAxis))
                boundaryRhs[owner] += conductance * voltage;
        }
        double rhsNorm = 0.0;
        for (int i = 0; i < count; ++i)
        {
            require(std::isfinite(diagonal[i]) && diagonal[i] > 0.0 && std::isfinite(boundaryRhs[i]),
                    "non-finite assembled diagonal or electrode RHS");
            rhsNorm = std::hypot(rhsNorm, boundaryRhs[i]);
        }
        require(std::isfinite(rhsNorm), "electrode RHS norm overflow");

        PicDomainBC bc;
        bc.left.type = bc.right.type = bc.bottom.type = bc.top.type = PicBCType::Neumann;
        if (xAxis)
        {
            bc.left.type = bc.right.type = PicBCType::Dirichlet;
            bc.left.value = 0.0; bc.right.value = voltage;
        }
        else
        {
            bc.bottom.type = bc.top.type = PicBCType::Dirichlet;
            bc.bottom.value = 0.0; bc.top.value = voltage;
        }
        // Solve the voltage relative to the low electrode. E and residuals
        // never subtract a potentially large common potential offset.
        CellField<double> rhs, shiftedPotential;
        rhs.initialize(mesh, 0.0);
        shiftedPotential.initialize(mesh, 0.0);
        PolarizationPoissonFVSolverCellCentered::Options options;
        options.useSuperLU = true;
        options.requireSuperLU = true;
        options.useCG = false;
        options.useIC = false;
        options.relTol = config.linearResidualTolerance;
        options.freezeDirectMatrixUntilPatternChange = false;
        options.useTemporalExtrapolatedGuess = false;
        PolarizationPoissonFVSolverCellCentered::SolveReport linearReport;
        report.stage = "solve";
        const bool solved = solver_.solve(mesh, graph_, bc, rhs, shiftedPotential, &epsilonFace, &epsilonCell,
                                         &linearReport, options, nullptr, nullptr, nullptr);
        report.usedSuperLU = linearReport.usedSuperLU;
        require(solved && linearReport.success,
                "variable-permittivity Poisson solver failed: " + linearReport.message);
        require(report.usedSuperLU, "perfect-dielectric electrostatics requires actual SuperLU execution; alternative backend rejected");
        for (int i = 0; i < count; ++i)
            require(std::isfinite(shiftedPotential.value[i]), "Poisson solution contains a non-finite value");

        // Independently rebuild the actual FV equation residual. The public
        // direct solver's finalResidual is not a verified algebraic residual.
        std::vector<double> fieldX(count, 0.0), fieldY(count, 0.0), weightX(count, 0.0), weightY(count, 0.0);
        std::vector<double> equationResidual(count, 0.0);
        std::vector<double> displacementNormal(faces.size(), 0.0);
        double absoluteBoundaryFlux = 0.0;
        for (std::size_t i = 0; i < faces.size(); ++i)
        {
            const auto& face = faces[i];
            const int owner = mesh.cells(face.K).index();
            double normalField = 0.0;
            if (face.boundary == Boundary::non)
                normalField = -(shiftedPotential.value[mesh.cells(face.L).index()] - shiftedPotential.value[owner]) / face.df;
            else if (isLowerElectrode(face.boundary, xAxis))
                normalField = shiftedPotential.value[owner] / face.df;
            else if (isUpperElectrode(face.boundary, xAxis))
                normalField = -(voltage - shiftedPotential.value[owner]) / face.df;
            const double flux = epsilonFace[i] * face.Af * normalField;
            require(std::isfinite(normalField) && std::isfinite(flux), "electric field or displacement flux overflow");
            equationResidual[owner] += flux;
            displacementNormal[i] = epsilonFace[i] * normalField;
            // Reconstruct D first, then E=D/epsilon_cell. This exactly
            // preserves the normal field of a charge-free planar layer.
            fieldX[owner] += face.n.x * flux;
            fieldY[owner] += face.n.y * flux;
            weightX[owner] += face.Af * face.n.x * face.n.x;
            weightY[owner] += face.Af * face.n.y * face.n.y;
            if (face.boundary != Boundary::non)
            {
                report.netBoundaryFlux += flux;
                absoluteBoundaryFlux += std::abs(flux);
                if (isLowerElectrode(face.boundary, xAxis)) report.lowerBoundaryFlux += flux;
                if (isUpperElectrode(face.boundary, xAxis)) report.upperBoundaryFlux += flux;
            }
        }
        double residualNorm = 0.0;
        for (double residual : equationResidual)
        {
            require(std::isfinite(residual), "non-finite discrete Gauss-law residual");
            residualNorm = std::hypot(residualNorm, residual);
        }
        report.relativeResidual = rhsNorm > 0.0 ? residualNorm / rhsNorm : residualNorm;
        report.boundaryFluxImbalance = absoluteBoundaryFlux > 0.0
            ? std::abs(report.netBoundaryFlux) / absoluteBoundaryFlux : std::abs(report.netBoundaryFlux);
        require(std::isfinite(report.relativeResidual) && report.relativeResidual <= config.linearResidualTolerance,
                "independently evaluated Poisson relative residual exceeds linearResidualTolerance");
        // The global flux balance sums all cell residuals. Its roundoff can
        // grow with cell count; retain an explicit size-aware conservation
        // check instead of treating a direct-solver success flag as proof.
        const double fluxTolerance = std::max(config.boundaryFluxTolerance.value_or(1e-9),
            128.0 * std::numeric_limits<double>::epsilon() * count);
        require(std::isfinite(report.netBoundaryFlux) && std::isfinite(report.lowerBoundaryFlux) &&
                std::isfinite(report.upperBoundaryFlux) && std::isfinite(absoluteBoundaryFlux) &&
                std::isfinite(report.boundaryFluxImbalance) && report.boundaryFluxImbalance <= fluxTolerance,
                "integrated boundary displacement flux fails charge-free closure");

        report.stage = "polarization";
        std::vector<double> nextPotential(count);
        auto nextPolarization = std::make_shared<DielectricPolarization>(mesh);
        nextPolarization->permittivityByCell = epsilonCell;
        nextPolarization->coefficientByNative.assign(mesh.vEdgeIndexNum() + mesh.hEdgeIndexNum(), 0.0);
        for (int i = 0; i < count; ++i)
        {
            require(std::isfinite(weightX[i]) && weightX[i] > 0.0 && std::isfinite(weightY[i]) && weightY[i] > 0.0,
                    "invalid cell electric-field reconstruction weights");
            fieldX[i] = (fieldX[i] / weightX[i]) / epsilonCell[i];
            fieldY[i] = (fieldY[i] / weightY[i]) / epsilonCell[i];
            require(std::isfinite(fieldX[i]) && std::isfinite(fieldY[i]), "non-finite reconstructed electric field");
        }
        for (std::size_t i = 0; i < faces.size(); ++i)
        {
            const auto& face = faces[i];
            // Each interior edge occurs twice in the oriented graph. Write
            // once, using the positive Cartesian normal; never double count.
            if (face.boundary != Boundary::non || face.n.x + face.n.y < 0.0) continue;
            const int k = mesh.cells(face.K).index(), l = mesh.cells(face.L).index();
            const double tangent = face.isHorizontal
                ? 0.5 * fieldX[k] + 0.5 * fieldX[l]
                : 0.5 * fieldY[k] + 0.5 * fieldY[l];
            const double dn = displacementNormal[i];
            // For constant Dn and Et this is exactly Delta(Tnn)/d, with
            // Tnn=.5*(Dn^2/epsilon-epsilon*Et^2). For general 2D fields it
            // approximates -.5*E^2 G(epsilon), NOT G of the whole stress:
            // genuine non-gradient deformation loads are retained.
            const double coefficient = DielectricPolarizationForce::gradientCoefficient(dn, tangent,
                epsilonCell[k], epsilonCell[l]);
            require(std::isfinite(coefficient), "non-finite dielectric polarization coefficient");
            const int row = face.isHorizontal ? mesh.vEdgeIndexNum() + mesh.hEdges(face.edgeLabel).index()
                : mesh.vEdges(face.edgeLabel).index();
            nextPolarization->coefficientByNative[row] = coefficient;
        }
        report.minPotential = std::numeric_limits<double>::infinity();
        report.maxPotential = -std::numeric_limits<double>::infinity();
        for (int i = 0; i < count; ++i)
        {
            const double field2 = fieldX[i] * fieldX[i] + fieldY[i] * fieldY[i];
            nextPotential[i] = shiftedPotential.value[i] + config.lowerPotential;
            require(std::isfinite(fieldX[i]) && std::isfinite(fieldY[i]) && std::isfinite(field2) &&
                    std::isfinite(nextPotential[i]),
                    "electric field or restored potential is non-finite");
            const double area = mesh.cells(mesh.cellIndexToLabel[i]).area();
            report.maxField = std::max(report.maxField, std::sqrt(field2));
            report.fieldEnergy += 0.5 * epsilonCell[i] * field2 * area;
            report.minPotential = std::min(report.minPotential, nextPotential[i]);
            report.maxPotential = std::max(report.maxPotential, nextPotential[i]);
        }
        require(std::isfinite(report.maxField) && std::isfinite(report.fieldEnergy),
                "integrated electric diagnostics overflow");

        auto next = std::make_shared<DielectricFields>(mesh);
        next->levelSet.assign(levelSet.value.get(), levelSet.value.get() + levelSet.size());
        next->axisymmetric = config.axisymmetric;
        next->potential = std::move(nextPotential);
        next->permittivity = std::move(epsilonCell);
        next->ex = std::move(fieldX); next->ey = std::move(fieldY);
        next->polarization = std::move(nextPolarization);
        report.stage = "complete";
        report.reason = "charge-free perfect-dielectric electrostatics and polarization coefficient evaluated";
        result.fields = std::move(next);
        report.success = true;
    }
    catch (const std::exception& error)
    {
        result.fields.reset();
        report.success = false;
        report.reason = error.what();
    }
    return result;
}
