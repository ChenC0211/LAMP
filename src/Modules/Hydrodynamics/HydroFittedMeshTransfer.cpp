#include "Modules/Hydrodynamics/HydroFittedMeshTransfer.h"

#include "Modules/Hydrodynamics/Hydrodynamics.h"
#include "Modules/Hydrodynamics/HydroFittedDomain.h"
#include "Modules/Hydrodynamics/HydroFittedMomentum.h"
#include "Modules/Hydrodynamics/HydroFittedProjection.h"
#include "Modules/Hydrodynamics/HydroContinuityBalance.h"
#include "Core/Mesh/QuadtreeMesh.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <optional>
#include <stdexcept>
#include <tuple>

namespace
{
using Domain = HydroFittedDomain;
using Type = Domain::BoundaryType;
using VelocityReport = HydroFittedMeshTransfer::VelocityReport;

void require(bool condition, const char* message)
{
    if (!condition) throw std::invalid_argument(std::string("[HydroFittedMeshTransfer] ") + message);
}

template<class Field>
void checkField(const Field& field, int count, const int* labels)
{
    require(field.size() == count && field.value, "field size does not match mesh");
    for (int i = 0; i < count; ++i)
        require(field.label(i) == labels[i] && std::isfinite(field.value[i]),
            "stale field label or non-finite transfer value");
}

void checkMaterials(const QuadtreeMesh& mesh, const Domain& domain,
    const CellField<double>& density, const CellField<double>& viscosity,
    const VertField<double>& vertexDensity, const VertField<double>& vertexViscosity)
{
    for (int i : domain.cellIndexByRow())
    {
        require(density.value[i] > 0 && viscosity.value[i] >= 0, "invalid active cell material");
        for (int label : mesh.cells(mesh.cellIndexToLabel[i]).cellToVert)
        {
            const int j = mesh.verts(label).index();
            require(vertexDensity.value[j] > 0 && vertexViscosity.value[j] >= 0,
                "invalid active vertex material; rebuild closure on transferred geometry");
        }
    }
}

Eigen::VectorXd compact(const QuadtreeMesh& mesh, const Domain& domain,
    const VEdgeField<double>& radial, const HEdgeField<double>& axial)
{
    checkField(radial, mesh.vEdgeIndexNum(), mesh.vEdgeIndexToLabel);
    checkField(axial, mesh.hEdgeIndexNum(), mesh.hEdgeIndexToLabel);
    Eigen::VectorXd result(domain.nativeByVelocityRow().size());
    for (int i = 0; i < result.size(); ++i)
    {
        const int row = domain.nativeByVelocityRow()[i];
        result[i] = row < radial.size() ? radial.value[row] : axial.value[row - radial.size()];
    }
    return result;
}

// Index-independent boundary identity. The first runtime policy keeps these
// faces AND their normal half-duals unchanged; boundary coarsening is rejected.
using BoundaryKey = std::tuple<bool, double, double, double, double, double,
    double, double, int, int, double, double, double>;
std::vector<BoundaryKey> boundaryKeys(const Domain& domain)
{
    std::vector<BoundaryKey> keys;
    for (const auto& f : domain.faces())
    {
        if ((f.negativeCell < 0 && f.positiveCell < 0) || f.boundary.type == Type::Interior) continue;
        require(f.boundary.type != Type::PressureOutlet,
            "mesh transfer requires explicit incoming material at open pressure boundaries");
        keys.emplace_back(f.horizontal, f.center.x, f.center.y, f.area, f.dualVolume,
            f.normalSpacing, f.tangentSpacing, f.outwardSign, static_cast<int>(f.boundary.type),
            f.boundary.id, f.boundary.pressure, f.boundary.incomingDensity, f.boundary.incomingViscosity);
    }
    std::sort(keys.begin(), keys.end());
    return keys;
}

void checkCoverage(const QuadtreeMesh& a, const Domain& ad, const QuadtreeMesh& b, const Domain& bd)
{
    // Dyadic cells on one background lattice are disjoint or nested. Checking
    // the centres in BOTH directions detects every changed active/solid leaf.
    for (int i = 0; i < a.cellIndexNum(); ++i)
    {
        const int label = b.findCellLabel(a.cells(a.cellIndexToLabel[i]).centerCoord());
        require(label >= 0 && label < b.cellLabelCapacity(), "uncovered source leaf");
        const int j = b.cells(label).index();
        require(j >= 0 && j < b.cellIndexNum(), "coverage lookup is not a target leaf");
        require((ad.cellRowByIndex()[i] >= 0) == (bd.cellRowByIndex()[j] >= 0),
            "mesh transfer changed the physical fluid/solid domain");
    }
}

// Diagnostics only: native half-dual rho*u sums and signed net boundary flux.
// Use the same original face density as the momentum/projection coefficient.
std::array<double, 3> integrals(const Domain& domain, const std::vector<double>& faceDensity,
    const Eigen::VectorXd& velocity)
{
    std::array<long double, 3> sums{};
    for (int native : domain.nativeByVelocityRow())
    {
        const auto& f = domain.faces()[native];
        const double u = velocity[f.velocityRow];
        const double r = faceDensity[f.velocityRow];
        sums[f.horizontal ? 1 : 0] += static_cast<long double>(r) * f.dualVolume * u;
        if (f.boundary.type != Type::Interior)
            sums[2] += static_cast<long double>(f.outwardSign) * f.area * u;
    }
    const std::array<double, 3> result{static_cast<double>(sums[0]), static_cast<double>(sums[1]),
        static_cast<double>(sums[2])};
    for (double value : result) require(std::isfinite(value), "non-finite transferred momentum/flux diagnostic");
    return result;
}

double pressureIntegral(const Domain& domain, const CellField<double>& pressure)
{
    long double sum = 0;
    for (int i = 0; i < static_cast<int>(domain.cellIndexByRow().size()); ++i)
        sum += static_cast<long double>(domain.cellVolumes()[i]) * pressure.value[domain.cellIndexByRow()[i]];
    const double result = static_cast<double>(sum);
    require(std::isfinite(result), "non-finite pressure integral");
    return result;
}
}

HydroFittedMeshTransfer::Result::Result() = default;
HydroFittedMeshTransfer::Result::~Result() = default;
HydroFittedMeshTransfer::Result::Result(Result&&) noexcept = default;
HydroFittedMeshTransfer::Result& HydroFittedMeshTransfer::Result::operator=(Result&&) noexcept = default;
HydroFittedMeshTransfer::~HydroFittedMeshTransfer() = default;

std::vector<double> HydroFittedMeshTransfer::nativeFaceDensities(const QuadtreeMesh& mesh,
    const HydroFittedDomain& domain, const Hydrodynamics& solver)
{
    std::vector<double> result(domain.nativeByVelocityRow().size());
    for (int row : domain.nativeByVelocityRow())
    {
        const auto& f = domain.faces()[row];
        const auto& edge = f.horizontal ? mesh.hEdges(f.edgeLabel) : mesh.vEdges(f.edgeLabel);
        result[f.velocityRow] = f.negativeCell >= 0 && f.positiveCell >= 0
            ? solver.pressureFaceData(mesh, edge, f.horizontal, solver.p).density
            : solver.rho.value[f.negativeCell >= 0 ? f.negativeCell : f.positiveCell];
        require(std::isfinite(result[f.velocityRow]) && result[f.velocityRow] > 0,
            "invalid native face density");
    }
    return result;
}

const Hydrodynamics& HydroFittedMeshTransfer::checkedSource(const Hydrodynamics& source)
{
    require(!source.stepTransaction_ && !source.stepProgress_.stopped && !source.stepProgress_.pendingAcceptance,
        "cannot capture a pending or stopped Hydro step");
    require(source.stepEntry_ == Hydrodynamics::StepEntry::Managed && source.initialStateSet_
        && source.materialState_ && source.stepProgress_.completedSteps > 0
        && source.stepProgress_.acceptedSteps == source.stepProgress_.completedSteps
        && source.stepProgress_.acceptedEndTime == source.stepProgress_.endTime,
        "fitted mesh transfer requires accepted Hydro time and material history");
    return source;
}

HydroFittedMeshTransfer::HydroFittedMeshTransfer(const QuadtreeMesh& mesh, const Hydrodynamics& source)
    : pressure(checkedSource(source).p), density(source.rho), viscosity(source.mu),
      vertexDensity(source.vertRho), vertexViscosity(source.vertMu),
      velocityR(source.ux), previousVelocityR(source.uxPre), predictorR(source.uxMid),
      velocityZ(source.uy), previousVelocityZ(source.uyPre), predictorZ(source.uyMid),
      sourceMesh_(&mesh), sourceDomain_(source.fittedDomain_), fixedDt_(source.fixedDt_),
      stepProgress_(source.stepProgress_), materialState_(*source.materialState_)
{
    require(sourceDomain_ && fixedDt_ > 0 && source.activeReport_.success
        && source.pressureMode_ == Hydrodynamics::PressureMode::Incremental,
        "capture requires a successfully accepted fitted CN/AB2 step");
    const double allowance = 64.0 * std::numeric_limits<double>::epsilon()
        * std::max({std::abs(stepProgress_.endTime), fixedDt_, std::abs(materialState_.sourceTime)});
    require(materialState_.evolution == HydroMaterialEvolution::Frozen
        || std::abs(materialState_.sourceTime - stepProgress_.endTime) <= allowance,
        "AMR requires materials updated to the completed Hydro state time");
    sourceDomain_->validateMeshIdentity(mesh);
    boundaryKeys(*sourceDomain_);
    for (const auto* f : {&pressure, &density, &viscosity})
        checkField(*f, mesh.cellIndexNum(), mesh.cellIndexToLabel);
    for (const auto* f : {&vertexDensity, &vertexViscosity})
        checkField(*f, mesh.vertIndexNum(), mesh.vertIndexToLabel);
    compact(mesh, *sourceDomain_, predictorR, predictorZ);
    checkMaterials(mesh, *sourceDomain_, density, viscosity, vertexDensity, vertexViscosity);
    const auto sourceFaceDensity = nativeFaceDensities(mesh, *sourceDomain_, source);
    const auto capture = [&](const auto& r, const auto& z, VelocityReport& report) {
        const auto u = compact(mesh, *sourceDomain_, r, z);
        for (int i = 0; i < static_cast<int>(sourceDomain_->velocityRowByNative().size()); ++i)
            if (sourceDomain_->velocityRowByNative()[i] < 0)
                require((i < r.size() ? r.value[i] : z.value[i - r.size()]) == 0,
                    "source constrained velocity/history is nonzero");
        const auto sums = integrals(*sourceDomain_, sourceFaceDensity, u);
        report.sourceMomentumR = sums[0]; report.sourceMomentumZ = sums[1];
        report.sourceBoundaryFlux = sums[2];
    };
    capture(velocityR, velocityZ, sourceReport_.current);
    capture(previousVelocityR, previousVelocityZ, sourceReport_.previous);
    sourceReport_.sourcePressureIntegral = pressureIntegral(*sourceDomain_, pressure);
}

void HydroFittedMeshTransfer::appendTransferFields(std::vector<VertField<double>*>& vertices,
    std::vector<HEdgeField<double>*>& horizontal, std::vector<VEdgeField<double>*>& vertical,
    std::vector<CellField<double>*>& cells)
{
    const auto absent = [](const auto& list, auto* field) {
        require(std::find(list.begin(), list.end(), field) == list.end(), "duplicate transfer field registration");
    };
    for (auto* f : {&vertexDensity, &vertexViscosity}) absent(vertices, f);
    for (auto* f : {&velocityZ, &previousVelocityZ, &predictorZ}) absent(horizontal, f);
    for (auto* f : {&velocityR, &previousVelocityR, &predictorR}) absent(vertical, f);
    for (auto* f : {&pressure, &density, &viscosity}) absent(cells, f);
    vertices.insert(vertices.end(), {&vertexDensity, &vertexViscosity});
    horizontal.insert(horizontal.end(), {&velocityZ, &previousVelocityZ, &predictorZ});
    vertical.insert(vertical.end(), {&velocityR, &previousVelocityR, &predictorR});
    cells.insert(cells.end(), {&pressure, &density, &viscosity});
}

HydroFittedMeshTransfer::Result HydroFittedMeshTransfer::prepare(QuadtreeMesh& mesh,
    const HydroFittedDomain& domain, const HydroFittedNumerics& numerics) const
{
    numerics.validate();
    sourceDomain_->validateMeshIdentity(*sourceMesh_);
    domain.validateMeshIdentity(mesh);
    require(domain.axisymmetric() && sourceDomain_->axisymmetric(), "only fitted RZ transfer is supported");
    require(mesh.Nx == sourceMesh_->Nx && mesh.Ny == sourceMesh_->Ny
        && mesh.minCoord.x == sourceMesh_->minCoord.x && mesh.minCoord.y == sourceMesh_->minCoord.y
        && mesh.maxCoord.x == sourceMesh_->maxCoord.x && mesh.maxCoord.y == sourceMesh_->maxCoord.y,
        "different background domain/lattice");
    require(boundaryKeys(domain) == boundaryKeys(*sourceDomain_), "physical boundary faces or data changed");
    checkCoverage(*sourceMesh_, *sourceDomain_, mesh, domain);
    checkCoverage(mesh, domain, *sourceMesh_, *sourceDomain_);
    for (const auto* f : {&pressure, &density, &viscosity}) checkField(*f, mesh.cellIndexNum(), mesh.cellIndexToLabel);
    for (const auto* f : {&vertexDensity, &vertexViscosity}) checkField(*f, mesh.vertIndexNum(), mesh.vertIndexToLabel);
    checkMaterials(mesh, domain, density, viscosity, vertexDensity, vertexViscosity);
    const auto current = compact(mesh, domain, velocityR, velocityZ);
    const auto previous = compact(mesh, domain, previousVelocityR, previousVelocityZ);
    compact(mesh, domain, predictorR, predictorZ);

    Result result;
    result.report = sourceReport_;
    auto candidate = std::make_unique<Hydrodynamics>(mesh);
    candidate->rho = density; candidate->mu = viscosity;
    candidate->vertRho = vertexDensity; candidate->vertMu = vertexViscosity;
    candidate->p = pressure;
    candidate->uxMid = predictorR; candidate->uyMid = predictorZ;
    candidate->activeOptions_.fittedNumerics = numerics;
    // Original native halo gate and virtual-child density, exactly as in the
    // fitted momentum/projection path. No new face-material interpolation.
    if (!domain.uniformLeaves())
    {
        HydroFittedMomentumLayout::build(mesh, domain);
        candidate->searchSupCoeffs(mesh);
    }
    const auto faceDensity = nativeFaceDensities(mesh, domain, *candidate);
    // Both history levels use the same domain and masses. Prepare lazily so
    // an already compatible transfer still performs no matrix assembly/solve.
    std::optional<HydroFittedProjection::System> projection;
    const auto project = [&](const Eigen::VectorXd& input, VelocityReport& report) {
        Eigen::VectorXd u = input;
        auto continuity = HydroContinuityBalance::LocalFlux(domain.divergence(), u);
        report.maxDivergenceBefore = continuity.maxDivergence(domain.cellVolumes());
        report.localContinuityResidual = continuity.relativeResidual();
        // Unit algebraic scale only; this is not an NS time step. The old
        // local continuity gate has no absolute floor or predictor scale.
        for (int pass = 0; report.localContinuityResidual > numerics.localContinuity && pass < 4; ++pass)
        {
            if (!projection) projection.emplace(HydroFittedProjection::assemble(domain, 1.0, faceDensity, u));
            else projection->updateRightHandSide(u);
            const auto& system = projection->equations();
            double residual = 0;
            const auto multiplier = candidate->checkedLinearSolve(system.matrix, system.rhs,
                "pressure", residual, numerics.pressureLinear);
            u = projection->correctVelocity(u, multiplier);
            ++report.projectionSolves;
            continuity = HydroContinuityBalance::LocalFlux(domain.divergence(), u);
            report.localContinuityResidual = continuity.relativeResidual();
        }
        require(report.localContinuityResidual <= numerics.localContinuity,
            "transferred velocity/history fails original local continuity tolerance");
        report.maxDivergenceAfter = continuity.maxDivergence(domain.cellVolumes());
        report.maxCorrection = (u - input).cwiseAbs().maxCoeff();
        const auto sums = integrals(domain, faceDensity, u);
        report.targetMomentumR = sums[0]; report.targetMomentumZ = sums[1]; report.targetBoundaryFlux = sums[2];
        return u;
    };
    const Eigen::VectorXd now = project(current, result.report.current);
    const Eigen::VectorXd pre = project(previous, result.report.previous);
    result.report.usedSuperLU = candidate->activeReport_.pressureUsedSuperLU;
    result.report.targetPressureIntegral = pressureIntegral(domain, candidate->p);
    const Eigen::VectorXd flux = domain.divergence() * now;
    for (int i = 0; i < static_cast<int>(domain.faces().size()); ++i)
    {
        const int j = domain.velocityRowByNative()[i];
        const int k = i < mesh.vEdgeIndexNum() ? i : i - mesh.vEdgeIndexNum();
        auto& u = i < mesh.vEdgeIndexNum() ? candidate->ux.value[k] : candidate->uy.value[k];
        auto& old = i < mesh.vEdgeIndexNum() ? candidate->uxPre.value[k] : candidate->uyPre.value[k];
        auto& mid = i < mesh.vEdgeIndexNum() ? candidate->uxMid.value[k] : candidate->uyMid.value[k];
        u = j >= 0 ? now[j] : 0; old = j >= 0 ? pre[j] : 0;
        if (j < 0) mid = 0;
    }
    for (int i = 0; i < static_cast<int>(domain.cellIndexByRow().size()); ++i)
    {
        const int index = domain.cellIndexByRow()[i];
        candidate->divU.value[index] = flux[i] / domain.cellVolumes()[i];
        const auto& cell = mesh.cells(mesh.cellIndexToLabel[index]);
        const double ur = .5 * (candidate->ux.value[mesh.vEdges(cell.cellToVEdge[0][0]).index()]
            + candidate->ux.value[mesh.vEdges(cell.cellToVEdge[1][0]).index()]);
        const double uz = .5 * (candidate->uy.value[mesh.hEdges(cell.cellToHEdge[0][0]).index()]
            + candidate->uy.value[mesh.hEdges(cell.cellToHEdge[1][0]).index()]);
        candidate->uMag.value[index] = std::hypot(ur, uz);
        candidate->uMagMax_ = std::max(candidate->uMagMax_, candidate->uMag.value[index]);
    }
    candidate->uMagMax_ = std::max(candidate->uMagMax_, now.cwiseAbs().maxCoeff());
    candidate->maxDivURes_ = result.report.current.maxDivergenceAfter;
    candidate->fittedDomain_ = std::make_shared<const Domain>(domain);
    candidate->fixedDt_ = fixedDt_;
    candidate->stepProgress_ = stepProgress_; // AMR is not a physical time step.
    candidate->initialStateSet_ = true;
    candidate->materialState_ = materialState_;
    candidate->stepEntry_ = Hydrodynamics::StepEntry::Managed;
    candidate->pressureMode_ = Hydrodynamics::PressureMode::Incremental;
    // No source matrices, transient pressure potential, or forces survive.
    candidate->pressureLinearCache_.reset();
    candidate->activeReport_ = HydroSolveReport{};
    candidate->activeReport_.success = true;
    candidate->activeReport_.stage = "fitted-mesh-transfer";
    result.solver = std::move(candidate);
    return result;
}
