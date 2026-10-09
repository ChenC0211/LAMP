#pragma once

#include "Core/Base/Base.h"
#include "Core/Base/Enumeration.h"
#include "Modules/Hydrodynamics/Hydrodynamics.h"
#include "Modules/Hydrodynamics/LegacyNodalForcing.h"
#include "Modules/Hydrodynamics/HydroContinuityBalance.h"

using namespace std;
using namespace Eigen;
typedef Eigen::Triplet<double> Trip;

template<typename Mesh>
Hydrodynamics::Hydrodynamics(Mesh& mesh)
{
    // vertRho.initialize(mesh, 0.0);
    // vertMu.initialize(mesh, 0.0);
    p.initialize(mesh, 0.0);
    pressurePotential_.initialize(mesh, 0.0);
    rho.initialize(mesh, 0.0);
    mu.initialize(mesh, 0.0);
    externalFx.initialize(mesh, 0.0);
    externalFy.initialize(mesh, 0.0);
    directFx.initialize(mesh, 0.0);
    directFy.initialize(mesh, 0.0);
    evapUx.initialize(mesh, 0.0);
    evapUy.initialize(mesh, 0.0);
    ux.initialize(mesh, 0.0);
    uy.initialize(mesh, 0.0);
    uxPre.initialize(mesh, 0.0);
    uyPre.initialize(mesh, 0.0);
    uxMid.initialize(mesh, 0.0);
    uyMid.initialize(mesh, 0.0);

    uMag.initialize(mesh, 0.0);
    divU.initialize(mesh, 0.0);
    magGradU.initialize(mesh, 0.0);
}

template<typename Mesh>
void Hydrodynamics::solve(Mesh& mesh, const double& dt, LegacyNodalState& sharpInterface, ImmersedGeometry<Mesh>& iGeometry)
{
    HydroSolveOptions options;
    options.surfaceForces = options.evaporation = options.directForce = true;
    options.closedDomainProjection = false;
    const auto report = solveStep(mesh, dt, false, options, &sharpInterface, &iGeometry, false);
    if (!report.success)
        throw std::runtime_error("[Hydrodynamics] " + report.stage + ": " + report.reason);
}

template<typename Mesh>
void Hydrodynamics::solveRZ(Mesh& mesh, const double& dt, LegacyNodalState& sharpInterface, ImmersedGeometry<Mesh>& iGeometry)
{
    HydroSolveOptions options;
    options.surfaceForces = options.evaporation = options.directForce = true;
    options.closedDomainProjection = false;
    const auto report = solveStep(mesh, dt, true, options, &sharpInterface, &iGeometry, false);
    if (!report.success)
        throw std::runtime_error("[Hydrodynamics] " + report.stage + ": " + report.reason);
}

template<typename Mesh>
void Hydrodynamics::assignConstantMaterialFields(Mesh& mesh, double density, double viscosity)
{
    if (!std::isfinite(density) || density <= 0.0 || !std::isfinite(viscosity) || viscosity < 0.0)
        throw std::invalid_argument("Hydro requires finite density > 0 and viscosity >= 0");
    rho.initialize(mesh, density);
    mu.initialize(mesh, viscosity);
    vertRho.initialize(mesh, density);
    vertMu.initialize(mesh, viscosity);
}

template<typename Mesh>
HydroSolveReport Hydrodynamics::solveStep(Mesh& mesh, double dt, bool axisymmetric,
    const HydroSolveOptions& options, LegacyNodalState* sharpInterface,
    ImmersedGeometry<Mesh>* iGeometry, bool fixedTimeStep,
    const CellField<double>* forceX, const CellField<double>* forceY,
    const VEdgeField<double>* faceForceX, const HEdgeField<double>* faceForceY, StepEntry entry,
    const Quantum::Hydrodynamics::HydroConservativeTransport* conservativeTransport,bool variableConservativeStep,const HydroMidpointData* midpoint)
{
    if (const auto rejected = claimStepEntry(entry)) return *rejected;
    activeOptions_ = options;
    activeReport_ = HydroSolveReport{};
    try
    {
        if (fittedDomain_)
            throw std::invalid_argument("cannot switch a fitted-fluid Hydro history to a complete-domain or legacy step");
        if (options.fittedNumerics)
            throw std::invalid_argument("explicit fittedNumerics require a fitted Hydro step");
        activeOptions_.completeConservation = options.completeConservation.value_or(HydroCompleteConservation{});
        activeOptions_.completeConservation->validate();
        const bool incrementalPressure = fixedTimeStep && options.closedDomainProjection;
        const bool mixedLeaves=static_cast<long long>(mesh.Nx)*mesh.Ny!=mesh.cellIndexNum();
        if(options.compatibleNativePressure && (entry!=StepEntry::Managed || !conservativeTransport || mixedLeaves || !options.boundaryPatches.empty()))
            throw std::invalid_argument("compatible native pressure currently requires managed uniform closed conservative transport");
        const bool physicalNative=conservativeTransport && (mixedLeaves || options.compatibleNativePressure);
        if(options.conservativeMidpoint && (!physicalNative || !midpoint || !options.compatibleNativePressure))
            throw std::invalid_argument("midpoint momentum requires its compatible native data");
        if(physicalNative && (entry!=StepEntry::Managed || !incrementalPressure || !options.boundaryPatches.empty()))
            throw std::invalid_argument("mixed conservative physical NS currently requires a managed closed domain");
        const PressureMode requestedPressureMode = incrementalPressure
            ? PressureMode::Incremental : PressureMode::Nonincremental;
        if (pressureMode_ != PressureMode::Unset && pressureMode_ != requestedPressureMode)
            throw std::invalid_argument("cannot switch Hydro pressure projection mode or use legacy cache with incremental pressure");
        if (!std::isfinite(dt) || dt <= 0.0)
            throw std::invalid_argument("dt must be finite and positive");
        if(variableConservativeStep && (!fixedTimeStep || !conservativeTransport || entry!=StepEntry::Managed))
            throw std::invalid_argument("variable Hydro steps require managed conservative mass and implicit advection");
        if (fixedTimeStep && !variableConservativeStep && fixedDt_ > 0.0 && std::abs(dt - fixedDt_) > 1e-12 * fixedDt_)
            throw std::invalid_argument("the standalone Hydro path requires a fixed dt");
        if (!std::isfinite(options.linearResidualTolerance) || options.linearResidualTolerance <= 0.0)
            throw std::invalid_argument("linearResidualTolerance must be finite and positive");
        if (!std::isfinite(options.accelerationX) || !std::isfinite(options.accelerationY))
            throw std::invalid_argument("body acceleration components must be finite");
        if (options.directForce && (options.maxDirectForceIterations <= 0
            || !std::isfinite(options.directForceTolerance) || options.directForceTolerance <= 0.0))
            throw std::invalid_argument("direct-force iteration limit and tolerance must be positive");
        if ((options.surfaceForces || options.evaporation || options.directForce) && (!sharpInterface || !iGeometry))
            throw std::invalid_argument("interface coupling requires LegacyNodalState and ImmersedGeometry");
        if (options.evaporation && iGeometry->levelSets.size() < 2)
            throw std::invalid_argument("legacy evaporation coupling requires level set 1");
        if (mesh.cellIndexNum() == 0)
            throw std::invalid_argument("mesh is empty");
        if (!options.boundaryPatches.empty() && (!incrementalPressure || !fixedTimeStep
            || options.surfaceForces || options.evaporation || options.directForce))
            throw std::invalid_argument("boundary patches require the checked fixed-dt incremental Hydro path");
        if (boundaryConfigurationFixed_)
        {
            const auto same = [](const HydroBoundaryPatch& a, const HydroBoundaryPatch& b) {
                return a.name == b.name && a.side == b.side && a.type == b.type
                    && a.start == b.start && a.end == b.end && a.value == b.value;
            };
            if (frozenBoundaryPatches_.size() != options.boundaryPatches.size()
                || !std::equal(frozenBoundaryPatches_.begin(), frozenBoundaryPatches_.end(),
                    options.boundaryPatches.begin(), same))
                throw std::invalid_argument("static boundary patches cannot change after a successful Hydro step");
        }
        auto validField = [](const auto& field, int size, bool positive, bool nonnegative)
        {
            if (field.size() != size || !field.value) return false;
            for (int i = 0; i < size; ++i)
                if (!std::isfinite(field.value[i]) || (positive && field.value[i] <= 0.0)
                    || (nonnegative && field.value[i] < 0.0)) return false;
            return true;
        };
        if ((forceX == nullptr) != (forceY == nullptr))
            throw std::invalid_argument("cell forcing requires both components");
        if (forceX)
        {
            if (!validField(*forceX, mesh.cellIndexNum(), false, false)
                || !validField(*forceY, mesh.cellIndexNum(), false, false))
                throw std::invalid_argument("cell force fields must match the mesh and be finite");
            if (forceX->value == externalFx.value || forceX->value == externalFy.value
                || forceY->value == externalFx.value || forceY->value == externalFy.value)
                throw std::invalid_argument("cell force inputs must not alias Hydro externalFx/Fy outputs");
        }
        if ((faceForceX == nullptr) != (faceForceY == nullptr))
            throw std::invalid_argument("face forcing requires both components");
        if (faceForceX)
        {
            if (!incrementalPressure)
                throw std::invalid_argument("face forcing requires fixed-dt closed-domain incremental pressure");
            if (!validField(*faceForceX, mesh.vEdgeIndexNum(), false, false)
                || !validField(*faceForceY, mesh.hEdgeIndexNum(), false, false))
                throw std::invalid_argument("face force fields must match the mesh and be finite");
            auto aliasesMutableState = [&](const double* values)
            {
                return values == ux.value.get() || values == uy.value.get() || values == uxPre.value.get()
                    || values == uyPre.value.get() || values == uxMid.value.get() || values == uyMid.value.get()
                    || values == p.value.get() || values == pressurePotential_.value.get()
                    || values == externalFx.value.get() || values == externalFy.value.get();
            };
            if (faceForceX->value == faceForceY->value
                || aliasesMutableState(faceForceX->value.get()) || aliasesMutableState(faceForceY->value.get()))
                throw std::invalid_argument("face force inputs must not alias each other or mutable Hydro state");
        }
        if (!validField(rho, mesh.cellIndexNum(), true, false)
            || !validField(mu, mesh.cellIndexNum(), false, true)
            || !validField(vertRho, mesh.vertIndexNum(), true, false)
            || !validField(vertMu, mesh.vertIndexNum(), false, true))
            throw std::invalid_argument("material fields must match the mesh, with finite density > 0 and viscosity >= 0");
        if (!validField(ux, mesh.vEdgeIndexNum(), false, false)
            || !validField(uxPre, mesh.vEdgeIndexNum(), false, false)
            || !validField(uy, mesh.hEdgeIndexNum(), false, false)
            || !validField(uyPre, mesh.hEdgeIndexNum(), false, false))
            throw std::invalid_argument("velocity and history fields must match the mesh and be finite");
        if (incrementalPressure && !validField(p, mesh.cellIndexNum(), false, false))
            throw std::invalid_argument("pressure history must match the mesh and be finite");
        for (int i = 0; i < mesh.cellIndexNum(); ++i)
        {
            const auto& cell = mesh.cells(mesh.cellIndexToLabel[i]);
            if (!(cell.size().x > 0.0) || !(cell.size().y > 0.0)
                || !std::isfinite(cell.size().x) || !std::isfinite(cell.size().y)
                || !std::isfinite(cell.centerCoord().x) || !std::isfinite(cell.centerCoord().y)
                || !std::isfinite(cell.area()) || (axisymmetric && !(cell.centerCoord().x > 0.0)))
                throw std::invalid_argument("invalid cell dimensions or nonpositive RZ cell-center radius");
            if (faceForceX && !physicalNative)
            {
                const auto& first = mesh.cells(mesh.cellIndexToLabel[0]);
                if (cell.level() != first.level() || !(cell.area() > 0.0)
                    || std::abs(cell.size().x - first.size().x) > 1e-12 * first.size().x
                    || std::abs(cell.size().y - first.size().y) > 1e-12 * first.size().y
                    || cell.cellToVEdge[0][1] != -1 || cell.cellToVEdge[1][1] != -1
                    || cell.cellToHEdge[0][1] != -1 || cell.cellToHEdge[1][1] != -1)
                    throw std::invalid_argument("face forcing requires a uniform fixed mesh without coarse/fine faces");
            }
        }
        if (axisymmetric)
            for (int i = 0; i < mesh.vertIndexNum(); ++i)
                if (!(mesh.verts(mesh.vertIndexToLabel[i]).coord().x >= 0.0)
                    || !std::isfinite(mesh.verts(mesh.vertIndexToLabel[i]).coord().x))
                    throw std::invalid_argument("RZ mesh must have finite nonnegative radius");
        // Complete-domain boundary preparation may write current/history
        // velocities. No rollback is promised after this point.
        activeReport_.stateChange = HydroStateChange::MayHaveChanged;
        prepareBoundaryPatches(mesh, axisymmetric);
        checkOutletVelocity(mesh, dt);
        if (options.closedDomainProjection && options.boundaryPatches.empty())
        {
            // Impose the same zero normal flux on predictor/history and corrected velocity.
            for (int i = 0; i < mesh.vEdgeIndexNum(); ++i)
            {
                const auto& edge = mesh.vEdges(mesh.vEdgeIndexToLabel[i]);
                if (edge.masterCell == -1 || edge.slaveCell == -1)
                    ux.value[i] = uxPre.value[i] = 0.0;
            }
            for (int i = 0; i < mesh.hEdgeIndexNum(); ++i)
            {
                const auto& edge = mesh.hEdges(mesh.hEdgeIndexToLabel[i]);
                if (edge.masterCell == -1 || edge.slaveCell == -1)
                    uy.value[i] = uyPre.value[i] = 0.0;
            }
        }
        activeReport_.stage = "surface-forces";
        externalFx.initialize(mesh, 0.0);
        externalFy.initialize(mesh, 0.0);
#pragma omp parallel for num_threads(thread_num)
        for (int i = 0; i < mesh.cellIndexNum(); ++i)
        {
            externalFx.value[i] = rho.value[i] * options.accelerationX + (forceX ? forceX->value[i] : 0.0);
            externalFy.value[i] = rho.value[i] * options.accelerationY + (forceY ? forceY->value[i] : 0.0);
        }
        if (options.surfaceForces)
            for (int levelSetI = 0; levelSetI < iGeometry->levelSets.size(); ++levelSetI)
                if (iGeometry->levelSets[levelSetI].moving())
                    calculateExternalForce(mesh, levelSetI, *sharpInterface);

        if (!validField(externalFx, mesh.cellIndexNum(), false, false)
            || !validField(externalFy, mesh.cellIndexNum(), false, false))
            throw std::runtime_error("assembled external force density is non-finite");

        if (incrementalPressure && options.boundaryPatches.empty())
        {
            // Consume explicit force inputs before modifying p: a caller may
            // supply p itself as a force input, just as with other cell fields.
            const double gauge = p.value[0];
            for (int i = 0; i < mesh.cellIndexNum(); ++i)
            {
                p.value[i] -= gauge;
                if (!std::isfinite(p.value[i]))
                    throw std::invalid_argument("pressure history gauge subtraction overflowed");
            }
        }
        activeReport_.stage = "momentum-assembly";
        if(physicalNative)solveConservativeNativeStep(mesh,dt,axisymmetric,*conservativeTransport,faceForceX,faceForceY,midpoint);
        else {
        searchSupCoeffs(mesh);
        HydroMomentumAssembly assembly;
        if (axisymmetric)
        {
            buildNativeMomentumMatrix<true>(mesh, dt, assembly);
            buildNativeMomentumRHS<true>(mesh, dt, assembly,nullptr,nullptr,nullptr,!conservativeTransport);
        }
        else
        {
            buildNativeMomentumMatrix<false>(mesh, dt, assembly);
            buildNativeMomentumRHS<false>(mesh, dt, assembly,nullptr,nullptr,nullptr,!conservativeTransport);
        }
        addBoundaryMomentum(mesh, dt, axisymmetric, assembly,!conservativeTransport);
        if(conservativeTransport) addConservativeMomentum(mesh,dt,axisymmetric,*conservativeTransport,assembly);
        if (fixedTimeStep && axisymmetric) addNativeRZHoop(mesh, dt, assembly);
        if (faceForceX)
        {
            activeReport_.stage = "face-forces";
            addFaceForce(mesh, axisymmetric, *faceForceX, *faceForceY, assembly.rhs);
        }
        if (incrementalPressure) addHistoricalPressure(mesh, axisymmetric, assembly.rhs);
        SparseMatrix<double> matrix;
        const VectorXd momentumRhs = std::move(assembly).finish(matrix);
        const VectorXd solution = checkedLinearSolve(matrix, momentumRhs, "momentum",
            activeReport_.momentumRelativeResidual);
        uxMid.initialize(mesh, 0.0);
        uyMid.initialize(mesh, 0.0);
#pragma omp parallel for num_threads(thread_num)
        for (int i = 0; i < mesh.vEdgeIndexNum(); ++i)
        {
            uxPre.value[i] = ux.value[i];
            uxMid.value[i] = solution[i];
        }
#pragma omp parallel for num_threads(thread_num)
        for (int i = 0; i < mesh.hEdgeIndexNum(); ++i)
        {
            uyPre.value[i] = uy.value[i];
            uyMid.value[i] = solution[i + mesh.vEdgeIndexNum()];
        }
        activeReport_.stage = "predictor-diagnostics";
        activeReport_.divergenceBeforeProjection = velocityDiagnostics(mesh, uxMid, uyMid, axisymmetric, false).first;
        activeReport_.stage = "pressure-assembly";
        bool nativeAdaptiveProjection = false;
        if (options.closedDomainProjection && options.boundaryPatches.empty())
            for (int i = 0; i < mesh.cellIndexNum(); ++i)
                nativeAdaptiveProjection |= mesh.cells(mesh.cellIndexToLabel[i]).supCellIndex >= 0;
        if (!options.boundaryPatches.empty())
        {
            projectBoundaryVelocity(mesh, dt, axisymmetric);
        }
        else
        {
            // Row strategies prepare equations only. The common pressure solve
            // and potential publication finish before the original correction.
            {
                const auto pressureSystem = nativeAdaptiveProjection
                    ? assembleNativeAdaptivePressure(mesh, dt, axisymmetric)
                    : (axisymmetric ? assembleCellPressure<true>(mesh, dt) : assembleCellPressure<false>(mesh, dt));
                const auto pressure = checkedLinearSolve(pressureSystem.matrix, pressureSystem.rhs,
                    "pressure", activeReport_.pressureRelativeResidual);
                pressurePotential_.initialize(mesh, 0.0);
                for (int cellI = 0; cellI < mesh.cellIndexNum(); cellI++)
                    pressurePotential_.value[cellI] = pressure[cellI];
            }
            if (axisymmetric) correctNativeVelocity<true>(mesh, dt);
            else correctNativeVelocity<false>(mesh, dt);
        }
        activeReport_.stage = "pressure-update";
        // Legacy callers need no historical p mapping after mesh adaptation.
        if (!incrementalPressure) p.initialize(mesh, 0.0);
        // Projection corrects velocity with phi, not with p_old + phi.
        // No constant-viscosity rotational formula is applied to variable mu/rho.
        for (int i = 0; i < mesh.cellIndexNum(); ++i)
        {
            p.value[i] = (incrementalPressure ? p.value[i] : 0.0) + pressurePotential_.value[i];
            if (!std::isfinite(p.value[i]))
                throw std::runtime_error("physical pressure update is non-finite");
        }
        if (incrementalPressure && options.boundaryPatches.empty())
        {
            const double gauge = p.value[0];
            for (int i = 0; i < mesh.cellIndexNum(); ++i)
            {
                p.value[i] -= gauge;
                if (!std::isfinite(p.value[i]))
                    throw std::runtime_error("physical pressure gauge subtraction overflowed");
            }
        }
        if (incrementalPressure) correctCoupledBalance(mesh, axisymmetric, matrix, momentumRhs);
        }
        checkOutletVelocity(mesh, dt, true);
        activeReport_.stage = "interface-correction";
        if (options.evaporation) applyEvapVel(mesh, *sharpInterface);
        if (options.directForce) applyDirectForce(mesh, dt, *sharpInterface, *iGeometry);
        activeReport_.stage = "final-diagnostics";
        const auto diagnostic = velocityDiagnostics(mesh, ux, uy, axisymmetric, true);
        activeReport_.divergenceAfterProjection = diagnostic.first;
        activeReport_.maxSpeed = diagnostic.second;
        maxDivURes_ = diagnostic.first;
        uMagMax_ = diagnostic.second;
        if (fixedTimeStep)
        {
            fixedDt_ = dt;
            frozenBoundaryPatches_ = options.boundaryPatches;
            boundaryConfigurationFixed_ = true;
        }
        pressureMode_ = requestedPressureMode;
        activeReport_.success = true;
        activeReport_.stateChange = HydroStateChange::Advanced;
        activeReport_.stage = "complete";
    }
    catch (const std::exception& error)
    {
        activeReport_.reason = error.what();
    }
    return activeReport_;
}

template<typename Mesh>
std::pair<double, double> Hydrodynamics::velocityDiagnostics(Mesh& mesh,
    const VEdgeField<double>& vx, const HEdgeField<double>& vy,
    bool axisymmetric, bool storeFields, CellField<double>* divergenceTarget,
    CellField<double>* magnitudeTarget)
{
    auto& storedDivergence = divergenceTarget ? *divergenceTarget : divU;
    auto& storedMagnitude = magnitudeTarget ? *magnitudeTarget : uMag;
    return HydroContinuityBalance::nativeVelocity(mesh, vx, vy, axisymmetric, thread_num,
        storeFields ? &storedDivergence : nullptr, storeFields ? &storedMagnitude : nullptr);
}

template<typename Mesh>
void Hydrodynamics::assignLevelSetMaterialFields(Mesh& mesh, const LevelSet& levelSet,
    double densityNegative, double densityPositive, double viscosityNegative,
    double viscosityPositive, double thickness)
{
    if (!std::isfinite(densityNegative) || densityNegative <= 0.0
        || !std::isfinite(densityPositive) || densityPositive <= 0.0
        || !std::isfinite(viscosityNegative) || viscosityNegative < 0.0
        || !std::isfinite(viscosityPositive) || viscosityPositive < 0.0
        || !std::isfinite(thickness) || thickness <= 0.0)
        throw std::invalid_argument("two-material Hydro requires finite positive densities/thickness and nonnegative viscosities");
    if (levelSet.size() != mesh.vertIndexNum() || !levelSet.value)
        throw std::invalid_argument("material level-set field must match the mesh vertices");
    for (int i = 0; i < levelSet.size(); ++i)
        if (!std::isfinite(levelSet.value[i]))
            throw std::invalid_argument("material level-set field contains non-finite values");

    hyperTanLS.setThickness(thickness);
    hyperTanLS.initialize(mesh, 0.0);
    hyperTanLS.value = levelSet.convertToHyperTan(thickness);
    vertRho.initialize(mesh, 0.0);
    vertMu.initialize(mesh, 0.0);
#pragma omp parallel for num_threads(thread_num)
    for (int i = 0; i < mesh.vertIndexNum(); ++i)
    {
        const double positiveFraction = hyperTanLS.value[i];
        vertRho.value[i] = densityNegative * (1.0 - positiveFraction) + densityPositive * positiveFraction;
        vertMu.value[i] = viscosityNegative * (1.0 - positiveFraction) + viscosityPositive * positiveFraction;
    }
    rho = std::move(*vertRho.toCellField(mesh));
    mu = std::move(*vertMu.toCellField(mesh));
}

template<typename Mesh>
void Hydrodynamics::updateProperties(Mesh& mesh, LevelSet& levelSet)
{
    requireLegacyMaterialEntry();
    hyperTanLS.initialize(mesh, 0.0);
    hyperTanLS.value = levelSet.convertToHyperTan(hyperTanLS.thickness());

    vertRho.initialize(mesh, 0.0);
    vertMu.initialize(mesh, 0.0);
#pragma omp parallel for num_threads(thread_num)
    for (int i = 0; i < mesh.vertIndexNum(); i++)
    {
        vertRho.value[i] = nonDim_rho_l * (1.0 - hyperTanLS.value[i]) + nonDim_rho_g * hyperTanLS.value[i];
        // 此处乘以常数无量纲Re数的倒数，减小后面的计算量
        vertMu.value[i] = number_Re_inv * (nonDim_mu_l * (1.0 - hyperTanLS.value[i]) + nonDim_mu_g * hyperTanLS.value[i]);
    }
    rho = std::move(*vertRho.toCellField(mesh));
    mu = std::move(*vertMu.toCellField(mesh));
}

template<typename Mesh>
void Hydrodynamics::calculateExternalForce(Mesh& mesh, const int levelSetI, LegacyNodalState& sharpInterface)
{
    externalFx.initialize(mesh, 0.0);
    externalFy.initialize(mesh, 0.0);

#pragma omp parallel for num_threads(thread_num)
    for (int i = 0; i < mesh.cellIndexNum(); i++)
    {
        externalFy.value[i] = rho.value[i] * nonDim_g;
    }

    LegacyNodalForcing::spreadSurface(mesh, sharpInterface.geometry().nodeIndexToLabel[levelSetI], sharpInterface.geometry().allNodes(levelSetI),
        sharpInterface.Fxs[levelSetI], sharpInterface.Fys[levelSetI], externalFx, externalFy);
}

template<typename Mesh>
void Hydrodynamics::interpolationVelToNode(Mesh& mesh, const int levelSetI, LegacyNodalState& sharpInterface)
{
    CellField<double> cellUx(mesh, 0.0);
    CellField<double> cellUy(mesh, 0.0);
#pragma omp parallel for num_threads(thread_num)
    for (int i = 0; i < mesh.cellIndexNum(); i++)
    {
        Cell& cell = mesh.cells(mesh.cellIndexToLabel[i]);
        // if (cell.supCellIndex == -1)
        {
            cellUx.value[i] = 0.5 * (ux.value[mesh.vEdges(cell.cellToVEdge[0][0]).index()] + ux.value[mesh.vEdges(cell.cellToVEdge[1][0]).index()]);
            cellUy.value[i] = 0.5 * (uy.value[mesh.hEdges(cell.cellToHEdge[0][0]).index()] + uy.value[mesh.hEdges(cell.cellToHEdge[1][0]).index()]);
        }
        // else
        // {
        //     int supCellIndex = cell.supCellIndex;
        //     double uWest = 0.0;
        //     for (int j = 0; j < supCells_[supCellIndex].vEdges[0].ind.size(); j++)
        //     {
        //         uWest += supCells_[supCellIndex].vEdges[0].coe[j] * ux.value[supCells_[supCellIndex].vEdges[0].ind[j]];
        //     }
        //     double uEast = 0.0;
        //     for (int j = 0; j < supCells_[supCellIndex].vEdges[1].ind.size(); j++)
        //     {
        //         uEast += supCells_[supCellIndex].vEdges[1].coe[j] * ux.value[supCells_[supCellIndex].vEdges[1].ind[j]];
        //     }
        //     cellUx.value[i] = 0.5 * (uWest + uEast);

        //     double vSouth = 0.0;
        //     for (int j = 0; j < supCells_[supCellIndex].hEdges[0].ind.size(); j++)
        //     {
        //         vSouth += supCells_[supCellIndex].hEdges[0].coe[j] * uy.value[supCells_[supCellIndex].hEdges[0].ind[j]];
        //     }
        //     double vNorth = 0.0;
        //     for (int j = 0; j < supCells_[supCellIndex].hEdges[1].ind.size(); j++)
        //     {
        //         vNorth += supCells_[supCellIndex].hEdges[1].coe[j] * uy.value[supCells_[supCellIndex].hEdges[1].ind[j]];
        //     }
        //     cellUy.value[i] = 0.5 * (vSouth + vNorth);
        // }
    }
    LegacyNodalForcing::interpolateNoSlip(sharpInterface.geometry().noSlipNodes(), sharpInterface.geometry().noSlipNodeNum,
        cellUx, cellUy, sharpInterface.noSlipUx, sharpInterface.noSlipUy);
}

template<typename Mesh>
void Hydrodynamics::searchSupCoeffs(Mesh& mesh)
{
#pragma omp parallel for num_threads(thread_num)
    for (int cellI = 0; cellI < mesh.cellIndexNum(); cellI++)
    {
        Cell& cell = mesh.cells(mesh.cellIndexToLabel[cellI]);
        cell.supCellIndex = -1; // 初始化辅助单元索引
    }

    int supCellNum = 0;
    for (int i = 0; i < mesh.vertIndexNum(); i++)
    {
        Vert& vert = mesh.verts(mesh.vertIndexToLabel[i]);
        if (vert.supCell != -1)
        {
            if (vert.supCell < 0 || vert.supCell >= mesh.cellLabelCapacity())
                throw std::invalid_argument("native Hydro hanging vertex has an invalid coarse cell");
            Cell& supCell = mesh.cells(vert.supCell);
            if (supCell.index() < 0 || supCell.index() >= mesh.cellIndexNum())
                throw std::invalid_argument("native Hydro hanging vertex does not reference a leaf");
            if (supCell.supCellIndex != -1)
            {
                continue;
            }
            
            supCell.supCellIndex = supCellNum;
            supCellNum++;
        }
    }

    // Catch missing derived mesh identities before any parallel stencil use.
    // A coarse/fine face must be 2:1 and its coarse leaf must own a stencil.
    for (int i = 0; i < mesh.cellIndexNum(); ++i)
    {
        const auto& cell = mesh.cells(mesh.cellIndexToLabel[i]);
        for (int side = 0; side < 4; ++side)
            for (int part = 0; part < 2; ++part)
            {
                const int label = cell.neiCell[side][part];
                if (label == -1) continue;
                if (label < 0 || label >= mesh.cellLabelCapacity())
                    throw std::invalid_argument("native Hydro has an invalid neighbor label");
                const auto& neighbor = mesh.cells(label);
                if (std::abs(cell.level() - neighbor.level()) > 1)
                    throw std::invalid_argument("native Hydro requires 2:1 coarse/fine faces");
                if (neighbor.level() > cell.level() && cell.supCellIndex < 0)
                    throw std::invalid_argument("native Hydro coarse/fine face lacks hanging-vertex metadata");
            }
    }

    supCells_ = make_unique<virtualSupCell[]>(supCellNum);
    supSonCells_ = make_unique<virtualSupCell[]>(supCellNum * 4);

    // Checked cell stencils are prepared serially: no exception may cross
    // the parallel velocity-stencil construction below.
    for (int i = 0; i < mesh.cellIndexNum(); ++i)
    {
        const auto& cell = mesh.cells(mesh.cellIndexToLabel[i]);
        if (cell.supCellIndex < 0) continue;
        for (int child = 0; child < 4; ++child)
            supSonCells_[cell.supCellIndex * 4 + child].center =
                HydroNativePressureStencil::childCenter(mesh, cell, child);
    }

#pragma omp parallel for num_threads(thread_num)
    for (int cellI = 0; cellI < mesh.cellIndexNum(); cellI++)
    {
        Cell& cell = mesh.cells(mesh.cellIndexToLabel[cellI]);
        if (cell.supCellIndex != -1)
        {
            // 垂直边-西
            if (cell.cellToVEdge[0][1] == -1)
            {
                Edge& vEdge00 = mesh.vEdges(cell.cellToVEdge[0][0]);
                supCells_[cell.supCellIndex].vEdges[0].ind.push_back(vEdge00.index());
                supCells_[cell.supCellIndex].vEdges[0].coe.push_back(1.0);
            }
            else
            {
                Edge& vEdge00 = mesh.vEdges(cell.cellToVEdge[0][0]);
                Edge& vEdge01 = mesh.vEdges(cell.cellToVEdge[0][1]);
                supCells_[cell.supCellIndex].vEdges[0].ind.push_back(vEdge00.index());
                supCells_[cell.supCellIndex].vEdges[0].ind.push_back(vEdge01.index());
                Cell& neiCell0 = mesh.cells(cell.neiCell[0][0]);
                Edge& supVEdge3 = mesh.vEdges(neiCell0.cellToVEdge[0][0]);
                supCells_[cell.supCellIndex].vEdges[0].ind.push_back(supVEdge3.index());
                vector<Vector2D> _coords(3);
                for (int i = 0; i < 3; i++)
                {
                    _coords[i] = mesh.vEdges(mesh.vEdgeIndexToLabel[supCells_[cell.supCellIndex].vEdges[0].ind[i]]).midpointCoord();
                }
                Vector2D _aimCoord(vEdge00.midpointCoord().x, cell.centerCoord().y);
                vector<double> _coeff = Interpolation::leastSquareInterpolationCoeff(_coords, _aimCoord);
                for (int i = 0; i < 3; i++)
                {
                    supCells_[cell.supCellIndex].vEdges[0].coe.push_back(_coeff[i]);
                }
            }
            // 垂直边-东
            if (cell.cellToVEdge[1][1] == -1)
            {
                Edge& vEdge10 = mesh.vEdges(cell.cellToVEdge[1][0]);
                supCells_[cell.supCellIndex].vEdges[1].ind.push_back(vEdge10.index());
                supCells_[cell.supCellIndex].vEdges[1].coe.push_back(1.0);
            }
            else
            {
                Edge& vEdge10 = mesh.vEdges(cell.cellToVEdge[1][0]);
                Edge& vEdge11 = mesh.vEdges(cell.cellToVEdge[1][1]);
                supCells_[cell.supCellIndex].vEdges[1].ind.push_back(vEdge10.index());
                supCells_[cell.supCellIndex].vEdges[1].ind.push_back(vEdge11.index());
                Cell& neiCell2 = mesh.cells(cell.neiCell[2][0]);
                Edge& supVEdge3 = mesh.vEdges(neiCell2.cellToVEdge[1][0]);
                supCells_[cell.supCellIndex].vEdges[1].ind.push_back(supVEdge3.index());
                vector<Vector2D> _coords(3);
                for (int i = 0; i < 3; i++)
                {
                    _coords[i] = mesh.vEdges(mesh.vEdgeIndexToLabel[supCells_[cell.supCellIndex].vEdges[1].ind[i]]).midpointCoord();
                }
                Vector2D _aimCoord(vEdge10.midpointCoord().x, cell.centerCoord().y);
                vector<double> _coeff = Interpolation::leastSquareInterpolationCoeff(_coords, _aimCoord);
                for (int i = 0; i < 3; i++)
                {
                    supCells_[cell.supCellIndex].vEdges[1].coe.push_back(_coeff[i]);
                }
            }
            // 水平边-南
            if (cell.cellToHEdge[0][1] == -1)
            {
                Edge& hEdge00 = mesh.hEdges(cell.cellToHEdge[0][0]);
                supCells_[cell.supCellIndex].hEdges[0].ind.push_back(hEdge00.index());
                supCells_[cell.supCellIndex].hEdges[0].coe.push_back(1.0);
            }
            else
            {
                Edge& hEdge00 = mesh.hEdges(cell.cellToHEdge[0][0]);
                Edge& hEdge01 = mesh.hEdges(cell.cellToHEdge[0][1]);
                supCells_[cell.supCellIndex].hEdges[0].ind.push_back(hEdge00.index());
                supCells_[cell.supCellIndex].hEdges[0].ind.push_back(hEdge01.index());
                Cell& neiCell1 = mesh.cells(cell.neiCell[1][0]);
                Edge& supHEdge3 = mesh.hEdges(neiCell1.cellToHEdge[0][0]);
                supCells_[cell.supCellIndex].hEdges[0].ind.push_back(supHEdge3.index());
                vector<Vector2D> _coords(3);
                for (int i = 0; i < 3; i++)
                {
                    _coords[i] = mesh.hEdges(mesh.hEdgeIndexToLabel[supCells_[cell.supCellIndex].hEdges[0].ind[i]]).midpointCoord();
                }
                Vector2D _aimCoord(cell.centerCoord().x, hEdge00.midpointCoord().y);
                vector<double> _coeff = Interpolation::leastSquareInterpolationCoeff(_coords, _aimCoord);
                for (int i = 0; i < 3; i++)
                {
                    supCells_[cell.supCellIndex].hEdges[0].coe.push_back(_coeff[i]);
                }
            }
            // 水平边-北
            if (cell.cellToHEdge[1][1] == -1)
            {
                Edge& hEdge10 = mesh.hEdges(cell.cellToHEdge[1][0]);
                supCells_[cell.supCellIndex].hEdges[1].ind.push_back(hEdge10.index());
                supCells_[cell.supCellIndex].hEdges[1].coe.push_back(1.0);
            }
            else
            {
                Edge& hEdge10 = mesh.hEdges(cell.cellToHEdge[1][0]);
                Edge& hEdge11 = mesh.hEdges(cell.cellToHEdge[1][1]);
                supCells_[cell.supCellIndex].hEdges[1].ind.push_back(hEdge10.index());
                supCells_[cell.supCellIndex].hEdges[1].ind.push_back(hEdge11.index());
                Cell& neiCell3 = mesh.cells(cell.neiCell[1][0]);
                Edge& supHEdge3 = mesh.hEdges(neiCell3.cellToHEdge[1][0]);
                supCells_[cell.supCellIndex].hEdges[1].ind.push_back(supHEdge3.index());
                vector<Vector2D> _coords(3);
                for (int i = 0; i < 3; i++)
                {
                    _coords[i] = mesh.hEdges(mesh.hEdgeIndexToLabel[supCells_[cell.supCellIndex].hEdges[1].ind[i]]).midpointCoord();
                }
                Vector2D _aimCoord(cell.centerCoord().x, hEdge10.midpointCoord().y);
                vector<double> _coeff = Interpolation::leastSquareInterpolationCoeff(_coords, _aimCoord);
                for (int i = 0; i < 3; i++)
                {
                    supCells_[cell.supCellIndex].hEdges[1].coe.push_back(_coeff[i]);
                }
            }

            // 子网格-0
            int supSonCellIndex0 = cell.supCellIndex * 4;
            // 子网格-0 // 垂直边-西
            if (cell.cellToVEdge[0][1] == -1)
            {
                Edge& vEdge00 = mesh.vEdges(cell.cellToVEdge[0][0]);
                supSonCells_[supSonCellIndex0].vEdges[0].ind.push_back(vEdge00.index());
                if (cell.cellToHEdge[0][1] != -1)
                {
                    Cell& neiCell = mesh.cells(cell.neiCell[1][1]);
                    Edge& supVEdge2 = mesh.vEdges(neiCell.cellToVEdge[0][0]);
                    supSonCells_[supSonCellIndex0].vEdges[0].ind.push_back(supVEdge2.index());
                }
                else
                {
                    Edge& vEdge10 = mesh.vEdges(cell.cellToVEdge[1][0]);
                    supSonCells_[supSonCellIndex0].vEdges[0].ind.push_back(vEdge10.index());
                }
                if (mesh.verts(vEdge00.edgeToVert[0]).vertToVEdge[0] != -1)
                {
                    Edge& supVEdge3 = mesh.vEdges(mesh.verts(vEdge00.edgeToVert[0]).vertToVEdge[0]);
                    supSonCells_[supSonCellIndex0].vEdges[0].ind.push_back(supVEdge3.index());
                }
                else
                {
                    Edge& supVEdge3 = mesh.vEdges(mesh.verts(vEdge00.edgeToVert[1]).vertToVEdge[1]);
                    supSonCells_[supSonCellIndex0].vEdges[0].ind.push_back(supVEdge3.index());
                }
                vector<Vector2D> _coords(3);
                for (int i = 0; i < 3; i++)
                {
                    _coords[i] = mesh.vEdges(mesh.vEdgeIndexToLabel[supSonCells_[supSonCellIndex0].vEdges[0].ind[i]]).midpointCoord();
                }
                Vector2D _aimCoord(vEdge00.midpointCoord().x, cell.centerCoord().y - cell.size().y * 0.25);
                vector<double> _coeff = Interpolation::leastSquareInterpolationCoeff(_coords, _aimCoord);
                for (int i = 0; i < 3; i++)
                {
                    supSonCells_[supSonCellIndex0].vEdges[0].coe.push_back(_coeff[i]);
                }
            }
            else
            {
                Edge& vEdge00 = mesh.vEdges(cell.cellToVEdge[0][0]);
                supSonCells_[supSonCellIndex0].vEdges[0].ind.push_back(vEdge00.index());
                supSonCells_[supSonCellIndex0].vEdges[0].coe.push_back(1.0);
            }
            // 子网格-0 // 垂直边-东
            {
                Edge& vEdge00 = mesh.vEdges(cell.cellToVEdge[0][0]);
                Edge& vEdge10 = mesh.vEdges(cell.cellToVEdge[1][0]);
                supSonCells_[supSonCellIndex0].vEdges[1].ind.push_back(vEdge00.index());
                supSonCells_[supSonCellIndex0].vEdges[1].ind.push_back(vEdge10.index());
                if (cell.cellToVEdge[0][1] != -1)
                {
                    Edge& vEdge01 = mesh.vEdges(cell.cellToVEdge[0][1]);
                    supSonCells_[supSonCellIndex0].vEdges[1].ind.push_back(vEdge01.index());
                }
                else if (cell.cellToVEdge[1][1] != -1)
                {
                    Edge& vEdge11 = mesh.vEdges(cell.cellToVEdge[1][1]);
                    supSonCells_[supSonCellIndex0].vEdges[1].ind.push_back(vEdge11.index());
                }
                else if (cell.neiCell[1][1] != -1)
                {
                    Cell& neiCell = mesh.cells(cell.neiCell[1][1]);
                    Edge& supVEdge3 = mesh.vEdges(neiCell.cellToVEdge[0][0]);
                    supSonCells_[supSonCellIndex0].vEdges[1].ind.push_back(supVEdge3.index());
                }
                else if (cell.neiCell[1][0] != -1)
                {
                    Cell& neiCell = mesh.cells(cell.neiCell[1][0]);
                    Edge& supVEdge3 = mesh.vEdges(neiCell.cellToVEdge[0][0]);
                    supSonCells_[supSonCellIndex0].vEdges[1].ind.push_back(supVEdge3.index());
                }
                else
                {
                    Edge& supVEdge3 = mesh.vEdges(mesh.verts(vEdge00.edgeToVert[1]).vertToVEdge[1]);
                    supSonCells_[supSonCellIndex0].vEdges[1].ind.push_back(supVEdge3.index());
                }
                vector<Vector2D> _coords(3);
                for (int i = 0; i < 3; i++)
                {
                    _coords[i] = mesh.vEdges(mesh.vEdgeIndexToLabel[supSonCells_[supSonCellIndex0].vEdges[1].ind[i]]).midpointCoord();
                }
                Vector2D _aimCoord(cell.centerCoord().x, cell.centerCoord().y - cell.size().y * 0.25);
                vector<double> _coeff = Interpolation::leastSquareInterpolationCoeff(_coords, _aimCoord);
                for (int i = 0; i < 3; i++)
                {
                    supSonCells_[supSonCellIndex0].vEdges[1].coe.push_back(_coeff[i]);
                }
            }
            // 子网格-0 // 水平边-南
            if (cell.cellToHEdge[0][1] == -1)
            {
                Edge& hEdge00 = mesh.hEdges(cell.cellToHEdge[0][0]);
                supSonCells_[supSonCellIndex0].hEdges[0].ind.push_back(hEdge00.index());
                if (cell.cellToVEdge[0][1] != -1)
                {
                    Cell& neiCell = mesh.cells(cell.neiCell[0][1]);
                    Edge& supHEdge2 = mesh.hEdges(neiCell.cellToHEdge[0][0]);
                    supSonCells_[supSonCellIndex0].hEdges[0].ind.push_back(supHEdge2.index());
                }
                else
                {
                    Edge& hEdge10 = mesh.hEdges(cell.cellToHEdge[1][0]);
                    supSonCells_[supSonCellIndex0].hEdges[0].ind.push_back(hEdge10.index());
                }
                if (mesh.verts(hEdge00.edgeToVert[0]).vertToHEdge[0] != -1)
                {
                    Edge& supHEdge3 = mesh.hEdges(mesh.verts(hEdge00.edgeToVert[0]).vertToHEdge[0]);
                    supSonCells_[supSonCellIndex0].hEdges[0].ind.push_back(supHEdge3.index());
                }
                else
                {
                    Edge& supHEdge3 = mesh.hEdges(mesh.verts(hEdge00.edgeToVert[1]).vertToHEdge[1]);
                    supSonCells_[supSonCellIndex0].hEdges[0].ind.push_back(supHEdge3.index());
                }
                vector<Vector2D> _coords(3);
                for (int i = 0; i < 3; i++)
                {
                    _coords[i] = mesh.hEdges(mesh.hEdgeIndexToLabel[supSonCells_[supSonCellIndex0].hEdges[0].ind[i]]).midpointCoord();
                }
                Vector2D _aimCoord(cell.centerCoord().x - cell.size().x * 0.25, hEdge00.midpointCoord().y);
                vector<double> _coeff = Interpolation::leastSquareInterpolationCoeff(_coords, _aimCoord);
                for (int i = 0; i < 3; i++)
                {
                    supSonCells_[supSonCellIndex0].hEdges[0].coe.push_back(_coeff[i]);
                }
            }
            else
            {
                Edge& hEdge00 = mesh.hEdges(cell.cellToHEdge[0][0]);
                supSonCells_[supSonCellIndex0].hEdges[0].ind.push_back(hEdge00.index());
                supSonCells_[supSonCellIndex0].hEdges[0].coe.push_back(1.0);
            }
            // 子网格-0 // 水平边-北
            {
                Edge& hEdge00 = mesh.hEdges(cell.cellToHEdge[0][0]);
                Edge& hEdge10 = mesh.hEdges(cell.cellToHEdge[1][0]);
                supSonCells_[supSonCellIndex0].hEdges[1].ind.push_back(hEdge00.index());
                supSonCells_[supSonCellIndex0].hEdges[1].ind.push_back(hEdge10.index());
                if (cell.cellToHEdge[0][1] != -1)
                {
                    Edge& hEdge01 = mesh.hEdges(cell.cellToHEdge[0][1]);
                    supSonCells_[supSonCellIndex0].hEdges[1].ind.push_back(hEdge01.index());
                }
                else if (cell.cellToHEdge[1][1] != -1)
                {
                    Edge& hEdge11 = mesh.hEdges(cell.cellToHEdge[1][1]);
                    supSonCells_[supSonCellIndex0].hEdges[1].ind.push_back(hEdge11.index());
                }
                else if (cell.neiCell[0][1] != -1)
                {
                    Cell& neiCell = mesh.cells(cell.neiCell[0][1]);
                    Edge& supHEdge3 = mesh.hEdges(neiCell.cellToHEdge[0][0]);
                    supSonCells_[supSonCellIndex0].hEdges[1].ind.push_back(supHEdge3.index());
                }
                else if (cell.neiCell[0][0] != -1)
                {
                    Cell& neiCell = mesh.cells(cell.neiCell[0][0]);
                    Edge& supHEdge3 = mesh.hEdges(neiCell.cellToHEdge[0][0]);
                    supSonCells_[supSonCellIndex0].hEdges[1].ind.push_back(supHEdge3.index());
                }
                else
                {
                    Edge& supHEdge3 = mesh.hEdges(mesh.verts(hEdge00.edgeToVert[1]).vertToHEdge[1]);
                    supSonCells_[supSonCellIndex0].hEdges[1].ind.push_back(supHEdge3.index());
                }
                vector<Vector2D> _coords(3);
                for (int i = 0; i < 3; i++)
                {
                    _coords[i] = mesh.hEdges(mesh.hEdgeIndexToLabel[supSonCells_[supSonCellIndex0].hEdges[1].ind[i]]).midpointCoord();
                }
                Vector2D _aimCoord(cell.centerCoord().x - cell.size().x * 0.25, cell.centerCoord().y);
                vector<double> _coeff = Interpolation::leastSquareInterpolationCoeff(_coords, _aimCoord);
                for (int i = 0; i < 3; i++)
                {
                    supSonCells_[supSonCellIndex0].hEdges[1].coe.push_back(_coeff[i]);
                }
            }

            // 子网格-1
            int supSonCellIndex1 = cell.supCellIndex * 4 + 1;
            // 子网格-1 // 垂直边-西
            supSonCells_[supSonCellIndex1].vEdges[0] = supSonCells_[supSonCellIndex0].vEdges[1];
            // 子网格-1 // 垂直边-东
            if (cell.cellToVEdge[1][1] == -1)
            {
                Edge& vEdge10 = mesh.vEdges(cell.cellToVEdge[1][0]);
                supSonCells_[supSonCellIndex1].vEdges[1].ind.push_back(vEdge10.index());
                if (cell.cellToHEdge[0][1] != -1)
                {
                    Cell& neiCell = mesh.cells(cell.neiCell[1][1]);
                    Edge& supVEdge2 = mesh.vEdges(neiCell.cellToVEdge[0][0]);
                    supSonCells_[supSonCellIndex1].vEdges[1].ind.push_back(supVEdge2.index());
                }
                else
                {
                    Edge& vEdge00 = mesh.vEdges(cell.cellToVEdge[0][0]);
                    supSonCells_[supSonCellIndex1].vEdges[1].ind.push_back(vEdge00.index());
                }
                if (mesh.verts(vEdge10.edgeToVert[0]).vertToVEdge[0] != -1)
                {
                    Edge& supVEdge3 = mesh.vEdges(mesh.verts(vEdge10.edgeToVert[0]).vertToVEdge[0]);
                    supSonCells_[supSonCellIndex1].vEdges[1].ind.push_back(supVEdge3.index());
                }
                else
                {
                    Edge& supVEdge3 = mesh.vEdges(mesh.verts(vEdge10.edgeToVert[1]).vertToVEdge[1]);
                    supSonCells_[supSonCellIndex1].vEdges[1].ind.push_back(supVEdge3.index());
                }
                vector<Vector2D> _coords(3);
                for (int i = 0; i < 3; i++)
                {
                    _coords[i] = mesh.vEdges(mesh.vEdgeIndexToLabel[supSonCells_[supSonCellIndex1].vEdges[1].ind[i]]).midpointCoord();
                }
                Vector2D _aimCoord(vEdge10.midpointCoord().x, cell.centerCoord().y - cell.size().y * 0.25);
                vector<double> _coeff = Interpolation::leastSquareInterpolationCoeff(_coords, _aimCoord);
                for (int i = 0; i < 3; i++)
                {
                    supSonCells_[supSonCellIndex1].vEdges[1].coe.push_back(_coeff[i]);
                }
            }
            else
            {
                Edge& vEdge10 = mesh.vEdges(cell.cellToVEdge[1][0]);
                supSonCells_[supSonCellIndex1].vEdges[1].ind.push_back(vEdge10.index());
                supSonCells_[supSonCellIndex1].vEdges[1].coe.push_back(1.0);
            }
            // 子网格-1 // 水平边-南
            if (cell.cellToHEdge[0][1] == -1)
            {
                Edge& hEdge00 = mesh.hEdges(cell.cellToHEdge[0][0]);
                supSonCells_[supSonCellIndex1].hEdges[0].ind.push_back(hEdge00.index());
                if (cell.cellToVEdge[1][1] != -1)
                {
                    Cell& neiCell = mesh.cells(cell.neiCell[2][1]);
                    Edge& supHEdge2 = mesh.hEdges(neiCell.cellToHEdge[0][0]);
                    supSonCells_[supSonCellIndex1].hEdges[0].ind.push_back(supHEdge2.index());
                }
                else
                {
                    if (cell.cellToHEdge[1][1] != -1)
                    {
                        Edge& hEdge10 = mesh.hEdges(cell.cellToHEdge[1][1]);
                        supSonCells_[supSonCellIndex1].hEdges[0].ind.push_back(hEdge10.index());
                    }
                    else
                    {
                        Edge& hEdge10 = mesh.hEdges(cell.cellToHEdge[1][0]);
                        supSonCells_[supSonCellIndex1].hEdges[0].ind.push_back(hEdge10.index());
                    }
                }
                if (mesh.verts(hEdge00.edgeToVert[1]).vertToHEdge[1] != -1)
                {
                    Edge& supHEdge3 = mesh.hEdges(mesh.verts(hEdge00.edgeToVert[1]).vertToHEdge[1]);
                    supSonCells_[supSonCellIndex1].hEdges[0].ind.push_back(supHEdge3.index());
                }
                else
                {
                    Edge& supHEdge3 = mesh.hEdges(mesh.verts(hEdge00.edgeToVert[0]).vertToHEdge[0]);
                    supSonCells_[supSonCellIndex1].hEdges[0].ind.push_back(supHEdge3.index());
                }
                vector<Vector2D> _coords(3);
                for (int i = 0; i < 3; i++)
                {
                    _coords[i] = mesh.hEdges(mesh.hEdgeIndexToLabel[supSonCells_[supSonCellIndex1].hEdges[0].ind[i]]).midpointCoord();
                }
                Vector2D _aimCoord(cell.centerCoord().x + cell.size().x * 0.25, hEdge00.midpointCoord().y);
                vector<double> _coeff = Interpolation::leastSquareInterpolationCoeff(_coords, _aimCoord);
                for (int i = 0; i < 3; i++)
                {
                    supSonCells_[supSonCellIndex1].hEdges[0].coe.push_back(_coeff[i]);
                }
            }
            else
            {
                Edge& hEdge01 = mesh.hEdges(cell.cellToHEdge[0][1]);
                supSonCells_[supSonCellIndex1].hEdges[0].ind.push_back(hEdge01.index());
                supSonCells_[supSonCellIndex1].hEdges[0].coe.push_back(1.0);
            }
            // 子网格-1 // 水平边-北
            {
                Edge& hEdge00 = mesh.hEdges(cell.cellToHEdge[0][0]);
                Edge& hEdge10 = mesh.hEdges(cell.cellToHEdge[1][0]);
                supSonCells_[supSonCellIndex1].hEdges[1].ind.push_back(hEdge00.index());
                supSonCells_[supSonCellIndex1].hEdges[1].ind.push_back(hEdge10.index());
                if (cell.cellToHEdge[0][1] != -1)
                {
                    Edge& hEdge01 = mesh.hEdges(cell.cellToHEdge[0][1]);
                    supSonCells_[supSonCellIndex1].hEdges[1].ind.push_back(hEdge01.index());
                }
                else if (cell.cellToHEdge[1][1] != -1)
                {
                    Edge& hEdge11 = mesh.hEdges(cell.cellToHEdge[1][1]);
                    supSonCells_[supSonCellIndex1].hEdges[1].ind.push_back(hEdge11.index());
                }
                else if (cell.neiCell[2][1] != -1)
                {
                    Cell& neiCell = mesh.cells(cell.neiCell[2][1]);
                    Edge& supHEdge3 = mesh.hEdges(neiCell.cellToHEdge[0][0]);
                    supSonCells_[supSonCellIndex1].hEdges[1].ind.push_back(supHEdge3.index());
                }
                else if (cell.neiCell[2][0] != -1)
                {
                    Cell& neiCell = mesh.cells(cell.neiCell[2][0]);
                    Edge& supHEdge3 = mesh.hEdges(neiCell.cellToHEdge[0][0]);
                    supSonCells_[supSonCellIndex1].hEdges[1].ind.push_back(supHEdge3.index());
                }
                else
                {
                    Edge& supHEdge3 = mesh.hEdges(mesh.verts(hEdge00.edgeToVert[0]).vertToHEdge[0]);
                    supSonCells_[supSonCellIndex1].hEdges[1].ind.push_back(supHEdge3.index());
                }
                vector<Vector2D> _coords(3);
                for (int i = 0; i < 3; i++)
                {
                    _coords[i] = mesh.hEdges(mesh.hEdgeIndexToLabel[supSonCells_[supSonCellIndex1].hEdges[1].ind[i]]).midpointCoord();
                }
                Vector2D _aimCoord(cell.centerCoord().x + cell.size().x * 0.25, cell.centerCoord().y);
                vector<double> _coeff = Interpolation::leastSquareInterpolationCoeff(_coords, _aimCoord);
                for (int i = 0; i < 3; i++)
                {
                    supSonCells_[supSonCellIndex1].hEdges[1].coe.push_back(_coeff[i]);
                }
            }
            // 子网格-2
            int supSonCellIndex2 = cell.supCellIndex * 4 + 2;
            // 子网格-2 // 垂直边-西
            {
                Edge& vEdge00 = mesh.vEdges(cell.cellToVEdge[0][0]);
                Edge& vEdge10 = mesh.vEdges(cell.cellToVEdge[1][0]);
                supSonCells_[supSonCellIndex2].vEdges[0].ind.push_back(vEdge00.index());
                supSonCells_[supSonCellIndex2].vEdges[0].ind.push_back(vEdge10.index());
                if (cell.cellToVEdge[0][1] != -1)
                {
                    Edge& vEdge01 = mesh.vEdges(cell.cellToVEdge[0][1]);
                    supSonCells_[supSonCellIndex2].vEdges[0].ind.push_back(vEdge01.index());
                }
                else if (cell.cellToVEdge[1][1] != -1)
                {
                    Edge& vEdge11 = mesh.vEdges(cell.cellToVEdge[1][1]);
                    supSonCells_[supSonCellIndex2].vEdges[0].ind.push_back(vEdge11.index());
                }
                else if (cell.neiCell[3][1] != -1)
                {
                    Cell& neiCell = mesh.cells(cell.neiCell[3][1]);
                    Edge& supVEdge3 = mesh.vEdges(neiCell.cellToVEdge[0][0]);
                    supSonCells_[supSonCellIndex2].vEdges[0].ind.push_back(supVEdge3.index());
                }
                else if (cell.neiCell[3][0] != -1)
                {
                    Cell& neiCell = mesh.cells(cell.neiCell[1][0]);
                    Edge& supVEdge3 = mesh.vEdges(neiCell.cellToVEdge[0][0]);
                    supSonCells_[supSonCellIndex2].vEdges[0].ind.push_back(supVEdge3.index());
                }
                else
                {
                    Edge& supVEdge3 = mesh.vEdges(mesh.verts(vEdge00.edgeToVert[0]).vertToVEdge[0]);
                    supSonCells_[supSonCellIndex2].vEdges[0].ind.push_back(supVEdge3.index());
                }
                vector<Vector2D> _coords(3);
                for (int i = 0; i < 3; i++)
                {
                    _coords[i] = mesh.vEdges(mesh.vEdgeIndexToLabel[supSonCells_[supSonCellIndex2].vEdges[0].ind[i]]).midpointCoord();
                }
                Vector2D _aimCoord(cell.centerCoord().x, cell.centerCoord().y + cell.size().y * 0.25);
                vector<double> _coeff = Interpolation::leastSquareInterpolationCoeff(_coords, _aimCoord);
                for (int i = 0; i < 3; i++)
                {
                    supSonCells_[supSonCellIndex2].vEdges[0].coe.push_back(_coeff[i]);
                }
            }
            // 子网格-2 // 垂直边-东
            if (cell.cellToVEdge[1][1] == -1)
            {
                Edge& vEdge10 = mesh.vEdges(cell.cellToVEdge[1][0]);
                supSonCells_[supSonCellIndex2].vEdges[1].ind.push_back(vEdge10.index());
                if (cell.cellToHEdge[1][1] != -1)
                {
                    Cell& neiCell = mesh.cells(cell.neiCell[3][1]);
                    Edge& supVEdge2 = mesh.vEdges(neiCell.cellToVEdge[0][0]);
                    supSonCells_[supSonCellIndex2].vEdges[1].ind.push_back(supVEdge2.index());
                }
                else
                {
                    if (cell.cellToVEdge[0][1] != -1)
                    {
                        Edge& vEdge01 = mesh.vEdges(cell.cellToVEdge[0][1]);
                        supSonCells_[supSonCellIndex2].vEdges[1].ind.push_back(vEdge01.index());
                    }
                    else
                    {
                        Edge& vEdge00 = mesh.vEdges(cell.cellToVEdge[0][0]);
                        supSonCells_[supSonCellIndex2].vEdges[1].ind.push_back(vEdge00.index());
                    }
                }
                if (mesh.verts(vEdge10.edgeToVert[1]).vertToVEdge[1] != -1)
                {
                    Edge& supVEdge3 = mesh.vEdges(mesh.verts(vEdge10.edgeToVert[1]).vertToVEdge[1]);
                    supSonCells_[supSonCellIndex2].vEdges[1].ind.push_back(supVEdge3.index());
                }
                else
                {
                    Edge& supVEdge3 = mesh.vEdges(mesh.verts(vEdge10.edgeToVert[0]).vertToVEdge[0]);
                    supSonCells_[supSonCellIndex2].vEdges[1].ind.push_back(supVEdge3.index());
                }
                vector<Vector2D> _coords(3);
                for (int i = 0; i < 3; i++)
                {
                    _coords[i] = mesh.vEdges(mesh.vEdgeIndexToLabel[supSonCells_[supSonCellIndex2].vEdges[1].ind[i]]).midpointCoord();
                }
                Vector2D _aimCoord(vEdge10.midpointCoord().x, cell.centerCoord().y + cell.size().y * 0.25);
                vector<double> _coeff = Interpolation::leastSquareInterpolationCoeff(_coords, _aimCoord);
                for (int i = 0; i < 3; i++)
                {
                    supSonCells_[supSonCellIndex2].vEdges[1].coe.push_back(_coeff[i]);
                }
            }
            else
            {
                Edge& vEdge11 = mesh.vEdges(cell.cellToVEdge[1][1]);
                supSonCells_[supSonCellIndex2].vEdges[1].ind.push_back(vEdge11.index());
                supSonCells_[supSonCellIndex2].vEdges[1].coe.push_back(1.0);
            }
            // 子网格-2 // 水平边-南
            supSonCells_[supSonCellIndex2].hEdges[0] = supSonCells_[supSonCellIndex1].hEdges[1];
            // 子网格-2 // 水平边-北
            if (cell.cellToHEdge[1][1] == -1)
            {
                Edge& hEdge10 = mesh.hEdges(cell.cellToHEdge[1][0]);
                supSonCells_[supSonCellIndex2].hEdges[1].ind.push_back(hEdge10.index());
                if (cell.cellToVEdge[1][1] != -1)
                {
                    Cell& neiCell = mesh.cells(cell.neiCell[2][1]);
                    Edge& supHEdge2 = mesh.hEdges(neiCell.cellToHEdge[0][0]);
                    supSonCells_[supSonCellIndex2].hEdges[1].ind.push_back(supHEdge2.index());
                }
                else
                {
                    if (cell.cellToHEdge[0][1] != -1)
                    {
                        Edge& hEdge01 = mesh.hEdges(cell.cellToHEdge[0][1]);
                        supSonCells_[supSonCellIndex2].hEdges[1].ind.push_back(hEdge01.index());
                    }
                    else
                    {
                        Edge& hEdge00 = mesh.hEdges(cell.cellToHEdge[0][0]);
                        supSonCells_[supSonCellIndex2].hEdges[1].ind.push_back(hEdge00.index());
                    }
                }
                if (mesh.verts(hEdge10.edgeToVert[1]).vertToHEdge[1] != -1)
                {
                    Edge& supHEdge3 = mesh.hEdges(mesh.verts(hEdge10.edgeToVert[1]).vertToHEdge[1]);
                    supSonCells_[supSonCellIndex2].hEdges[1].ind.push_back(supHEdge3.index());
                }
                else
                {
                    Edge& supHEdge3 = mesh.hEdges(mesh.verts(hEdge10.edgeToVert[0]).vertToHEdge[0]);
                    supSonCells_[supSonCellIndex2].hEdges[1].ind.push_back(supHEdge3.index());
                }
                vector<Vector2D> _coords(3);
                for (int i = 0; i < 3; i++)
                {
                    _coords[i] = mesh.hEdges(mesh.hEdgeIndexToLabel[supSonCells_[supSonCellIndex2].hEdges[1].ind[i]]).midpointCoord();
                }
                Vector2D _aimCoord(cell.centerCoord().x + cell.size().x * 0.25, hEdge10.midpointCoord().y);
                vector<double> _coeff = Interpolation::leastSquareInterpolationCoeff(_coords, _aimCoord);
                for (int i = 0; i < 3; i++)
                {
                    supSonCells_[supSonCellIndex2].hEdges[1].coe.push_back(_coeff[i]);
                }
            }
            else
            {
                Edge& hEdge11 = mesh.hEdges(cell.cellToHEdge[1][1]);
                supSonCells_[supSonCellIndex2].hEdges[1].ind.push_back(hEdge11.index());
                supSonCells_[supSonCellIndex2].hEdges[1].coe.push_back(1.0);
            }
            // 子网格-3
            int supSonCellIndex3 = cell.supCellIndex * 4 + 3;
            // 子网格-3 // 垂直边-西
            if (cell.cellToVEdge[0][1] == -1)
            {
                Edge& vEdge00 = mesh.vEdges(cell.cellToVEdge[0][0]);
                supSonCells_[supSonCellIndex3].vEdges[0].ind.push_back(vEdge00.index());
                if (cell.cellToHEdge[1][1] != -1)
                {
                    Cell& neiCell = mesh.cells(cell.neiCell[3][1]);
                    Edge& supVEdge2 = mesh.vEdges(neiCell.cellToVEdge[0][0]);
                    supSonCells_[supSonCellIndex3].vEdges[0].ind.push_back(supVEdge2.index());
                }
                else
                {
                    if (cell.cellToVEdge[1][1] != -1)
                    {
                        Edge& vEdge10 = mesh.vEdges(cell.cellToVEdge[1][1]);
                        supSonCells_[supSonCellIndex3].vEdges[0].ind.push_back(vEdge10.index());
                    }
                    else
                    {
                        Edge& vEdge10 = mesh.vEdges(cell.cellToVEdge[1][0]);
                        supSonCells_[supSonCellIndex3].vEdges[0].ind.push_back(vEdge10.index());
                    }
                }
                if (mesh.verts(vEdge00.edgeToVert[1]).vertToVEdge[1] != -1)
                {
                    Edge& supVEdge3 = mesh.vEdges(mesh.verts(vEdge00.edgeToVert[1]).vertToVEdge[1]);
                    supSonCells_[supSonCellIndex3].vEdges[0].ind.push_back(supVEdge3.index());
                }
                else
                {
                    Edge& supVEdge3 = mesh.vEdges(mesh.verts(vEdge00.edgeToVert[0]).vertToVEdge[0]);
                    supSonCells_[supSonCellIndex3].vEdges[0].ind.push_back(supVEdge3.index());
                }
                vector<Vector2D> _coords(3);
                for (int i = 0; i < 3; i++)
                {
                    _coords[i] = mesh.vEdges(mesh.vEdgeIndexToLabel[supSonCells_[supSonCellIndex3].vEdges[0].ind[i]]).midpointCoord();
                }
                Vector2D _aimCoord(vEdge00.midpointCoord().x, cell.centerCoord().y + cell.size().y * 0.25);
                vector<double> _coeff = Interpolation::leastSquareInterpolationCoeff(_coords, _aimCoord);
                for (int i = 0; i < 3; i++)
                {
                    supSonCells_[supSonCellIndex3].vEdges[0].coe.push_back(_coeff[i]);
                }
            }
            else
            {
                Edge& vEdge01 = mesh.vEdges(cell.cellToVEdge[0][1]);
                supSonCells_[supSonCellIndex3].vEdges[0].ind.push_back(vEdge01.index());
                supSonCells_[supSonCellIndex3].vEdges[0].coe.push_back(1.0);
            }
            // 子网格-3 // 垂直边-东
            supSonCells_[supSonCellIndex3].vEdges[1] = supSonCells_[supSonCellIndex2].vEdges[0];
            // 子网格-3 // 水平边-南
            supSonCells_[supSonCellIndex3].hEdges[0] = supSonCells_[supSonCellIndex0].hEdges[1];
            // 子网格-3 // 水平边-北
            if (cell.cellToHEdge[1][1] == -1)
            {
                Edge& hEdge10 = mesh.hEdges(cell.cellToHEdge[1][0]);
                supSonCells_[supSonCellIndex3].hEdges[1].ind.push_back(hEdge10.index());
                if (cell.cellToVEdge[0][1] != -1)
                {
                    Cell& neiCell = mesh.cells(cell.neiCell[0][1]);
                    Edge& supHEdge2 = mesh.hEdges(neiCell.cellToHEdge[0][0]);
                    supSonCells_[supSonCellIndex3].hEdges[1].ind.push_back(supHEdge2.index());
                }
                else
                {
                    Edge& hEdge00 = mesh.hEdges(cell.cellToHEdge[0][0]);
                    supSonCells_[supSonCellIndex3].hEdges[1].ind.push_back(hEdge00.index());
                }
                if (mesh.verts(hEdge10.edgeToVert[0]).vertToHEdge[0] != -1)
                {
                    Edge& supHEdge3 = mesh.hEdges(mesh.verts(hEdge10.edgeToVert[0]).vertToHEdge[0]);
                    supSonCells_[supSonCellIndex3].hEdges[1].ind.push_back(supHEdge3.index());
                }
                else
                {
                    Edge& supHEdge3 = mesh.hEdges(mesh.verts(hEdge10.edgeToVert[1]).vertToHEdge[1]);
                    supSonCells_[supSonCellIndex3].hEdges[1].ind.push_back(supHEdge3.index());
                }
                vector<Vector2D> _coords(3);
                for (int i = 0; i < 3; i++)
                {
                    _coords[i] = mesh.hEdges(mesh.hEdgeIndexToLabel[supSonCells_[supSonCellIndex3].hEdges[1].ind[i]]).midpointCoord();
                }
                Vector2D _aimCoord(cell.centerCoord().x - cell.size().x * 0.25, hEdge10.midpointCoord().y);
                vector<double> _coeff = Interpolation::leastSquareInterpolationCoeff(_coords, _aimCoord);
                for (int i = 0; i < 3; i++)
                {
                    supSonCells_[supSonCellIndex3].hEdges[1].coe.push_back(_coeff[i]);
                }
            }
            else
            {
                Edge& hEdge10 = mesh.hEdges(cell.cellToHEdge[1][0]);
                supSonCells_[supSonCellIndex3].hEdges[1].ind.push_back(hEdge10.index());
                supSonCells_[supSonCellIndex3].hEdges[1].coe.push_back(1.0);
            }
        }
    }
}









template<typename Mesh>
void Hydrodynamics::applyEvapVel(Mesh& mesh, LegacyNodalState& sharpInterface)
{
    evapUx.initialize(mesh, 0.0);
    evapUy.initialize(mesh, 0.0);

    LegacyNodalForcing::spreadEvaporation(sharpInterface.geometry().nodeIndexToLabel[1], sharpInterface.geometry().allNodes(1), sharpInterface.evaporationVelocity(1), evapUx, evapUy);

    VEdgeField<double> edgeEvapUx = std::move(*std::move(*evapUx.toVertField(mesh)).toVEdgeField(mesh));
    HEdgeField<double> edgeEvapUy = std::move(*std::move(*evapUy.toVertField(mesh)).toHEdgeField(mesh));

#pragma omp parallel for num_threads(thread_num)
    for (int i = 0; i < mesh.vEdgeIndexNum(); i++)
    {
        ux.value[i] += edgeEvapUx.value[i];
    }
#pragma omp parallel for num_threads(thread_num)
    for (int i = 0; i < mesh.hEdgeIndexNum(); i++)
    {
        uy.value[i] += edgeEvapUy.value[i];
    }
}

template<typename Mesh>
void Hydrodynamics::applyDirectForce(Mesh& mesh, const double& dt, LegacyNodalState& sharpInterface, ImmersedGeometry<Mesh>& iGeometry)
{
    directFx.initialize(mesh, 0.0);
    directFy.initialize(mesh, 0.0);

    double _alpha = nonDim_rho_l / dt;
    double _omega = iterCoeff_ / _alpha;
    const double aimRes = activeOptions_.directForceTolerance;

    int iter = 0;
    double _residual = 0.0;
    do
    {
        directFx.initialize(mesh, 0.0);
        directFy.initialize(mesh, 0.0);

        // for (int levelSetI = 0; levelSetI < iGeometry.levelSets.size(); levelSetI++)
        // {
        //     if (!iGeometry.levelSets[levelSetI].moving())
        //     {
        //         interpolationVelToNode(mesh, levelSetI, sharpInterface);
        //     }
        // }
        interpolationVelToNode(mesh, 1, sharpInterface);
        _residual = LegacyNodalForcing::fixedVelocityForce(sharpInterface.geometry().noSlipNodeNum, sharpInterface.noSlipUx, sharpInterface.noSlipUy,
            Vector2D(0.0, 0.0), _alpha, sharpInterface.noSlipFx, sharpInterface.noSlipFy);
        LegacyNodalForcing::spreadVolume(sharpInterface.geometry().noSlipNodes(), sharpInterface.geometry().noSlipNodeNum,
            sharpInterface.noSlipFx, sharpInterface.noSlipFy, directFx, directFy);
        
#pragma omp parallel for num_threads(thread_num)
        for (int i = 0; i < mesh.vEdgeIndexNum(); i++)
        {
            Edge& vEdge = mesh.vEdges(mesh.vEdgeIndexToLabel[i]);
            if (vEdge.masterCell != -1 && vEdge.slaveCell != -1)
            {
                Cell& masterCell = mesh.cells(vEdge.masterCell);
                Cell& slaveCell = mesh.cells(vEdge.slaveCell);
                double _forceX = 0.0;
                if (masterCell.level() < slaveCell.level())
                {
                    _forceX = (directFx.value[masterCell.index()] + 2.0 * directFx.value[slaveCell.index()]) / 3.0;
                }
                else if (masterCell.level() > slaveCell.level())
                {
                    _forceX = (2.0 * directFx.value[masterCell.index()] + directFx.value[slaveCell.index()]) / 3.0;
                }
                else
                {
                    _forceX = 0.5 * (directFx.value[masterCell.index()] + directFx.value[slaveCell.index()]);
                }
                ux.value[i] += _forceX * _omega;
            }
            else if (vEdge.masterCell != -1)
            {
                Cell& masterCell = mesh.cells(vEdge.masterCell);
                ux.value[i] += directFx.value[masterCell.index()] * _omega;
            }
            else if (vEdge.slaveCell != -1)
            {
                Cell& slaveCell = mesh.cells(vEdge.slaveCell);
                ux.value[i] += directFx.value[slaveCell.index()] * _omega;
            }
        }
#pragma omp parallel for num_threads(thread_num)
        for (int i = 0; i < mesh.hEdgeIndexNum(); i++)
        {
            Edge& hEdge = mesh.hEdges(mesh.hEdgeIndexToLabel[i]);
            if (hEdge.masterCell != -1 && hEdge.slaveCell != -1)
            {
                Cell& masterCell = mesh.cells(hEdge.masterCell);
                Cell& slaveCell = mesh.cells(hEdge.slaveCell);
                double _forceY = 0.0;
                if (masterCell.level() < slaveCell.level())
                {
                    _forceY = (directFy.value[masterCell.index()] + 2.0 * directFy.value[slaveCell.index()]) / 3.0;
                }
                else if (masterCell.level() > slaveCell.level())
                {
                    _forceY = (2.0 * directFy.value[masterCell.index()] + directFy.value[slaveCell.index()]) / 3.0;
                }
                else
                {
                    _forceY = 0.5 * (directFy.value[masterCell.index()] + directFy.value[slaveCell.index()]);
                }
                uy.value[i] += _forceY * _omega;
            }
            else if (hEdge.masterCell != -1)
            {
                Cell& masterCell = mesh.cells(hEdge.masterCell);
                uy.value[i] += directFy.value[masterCell.index()] * _omega;
            }
            else if (hEdge.slaveCell != -1)
            {
                Cell& slaveCell = mesh.cells(hEdge.slaveCell);
                uy.value[i] += directFy.value[slaveCell.index()] * _omega;
            }
        }

        if (iter % 200 == 0)
        {
            cout << "[Hydrodynamics] Iter=" << iter << ", target velocity residual=" << _residual << endl;
        }

        iter++;

        if (!std::isfinite(_residual))
            throw std::runtime_error("direct-force iteration returned a non-finite residual");
    } while (_residual > aimRes && iter < activeOptions_.maxDirectForceIterations);
    if (_residual > aimRes)
        throw std::runtime_error("direct-force iteration limit reached, residual=" + std::to_string(_residual));
    cout << "[Hydrodynamics] Iter=" << iter << ", target velocity residual=" << _residual << endl;
}

template<typename Mesh>
void Hydrodynamics::writeCache(const Mesh& mesh, const std::string& filePath) const
{
    if (stepEntry_ == StepEntry::Managed)
        throw std::logic_error("managed Hydro state cannot be written as an unversioned legacy cache");
    if (pressureMode_ == PressureMode::Incremental)
        throw std::logic_error("legacy Hydro cache cannot store incremental-pressure state; versioned restart is not implemented");
    cout << "[Hydrodynamics] Writing flow-field cache..." << endl;

    ostringstream fileName;
	fileName << filePath << "hydrodynamics_cache.bin";

    ofstream cache(fileName.str(), ios::binary | ios::out);
	if (cache.fail())
	{
		cout << "Writing cache \"" << fileName.str() << "\" failed!" << endl;
		exit(-1);
	}

    for (int i = 0; i < mesh.vertIndexNum(); i++)
    {
        cache.write(reinterpret_cast<const char*>(&vertRho.value[i]), sizeof(double));
        cache.write(reinterpret_cast<const char*>(&vertMu.value[i]), sizeof(double));
    }
    for (int i = 0; i < mesh.vEdgeIndexNum(); i++)
    {
        cache.write(reinterpret_cast<const char*>(&ux.value[i]), sizeof(double));
        cache.write(reinterpret_cast<const char*>(&uxPre.value[i]), sizeof(double));
        cache.write(reinterpret_cast<const char*>(&uxMid.value[i]), sizeof(double));
    }
    for (int i = 0; i < mesh.hEdgeIndexNum(); i++)
    {
        cache.write(reinterpret_cast<const char*>(&uy.value[i]), sizeof(double));
        cache.write(reinterpret_cast<const char*>(&uyPre.value[i]), sizeof(double));
        cache.write(reinterpret_cast<const char*>(&uyMid.value[i]), sizeof(double));
    }
    for (int i = 0; i < mesh.cellIndexNum(); i++)
    {
        cache.write(reinterpret_cast<const char*>(&p.value[i]), sizeof(double));
        cache.write(reinterpret_cast<const char*>(&rho.value[i]), sizeof(double));
        cache.write(reinterpret_cast<const char*>(&mu.value[i]), sizeof(double));
        cache.write(reinterpret_cast<const char*>(&externalFx.value[i]), sizeof(double));
        cache.write(reinterpret_cast<const char*>(&externalFy.value[i]), sizeof(double));
        cache.write(reinterpret_cast<const char*>(&directFx.value[i]), sizeof(double));
        cache.write(reinterpret_cast<const char*>(&directFy.value[i]), sizeof(double));
        cache.write(reinterpret_cast<const char*>(&evapUx.value[i]), sizeof(double));
        cache.write(reinterpret_cast<const char*>(&evapUy.value[i]), sizeof(double));
        cache.write(reinterpret_cast<const char*>(&uMag.value[i]), sizeof(double));
        cache.write(reinterpret_cast<const char*>(&divU.value[i]), sizeof(double));
        cache.write(reinterpret_cast<const char*>(&magGradU.value[i]), sizeof(double));
    }

    cache.write(reinterpret_cast<const char*>(&maxDivURes_), sizeof(double));
    cache.write(reinterpret_cast<const char*>(&uMagMax_), sizeof(double));

    cache.close();
    cout << "[Hydrodynamics] Flow-field cache written." << endl;
}

template<typename Mesh>
void Hydrodynamics::readCache(Mesh& mesh, const std::string& filePath)
{
    if (const auto rejected = claimStepEntry(StepEntry::Legacy))
        throw std::logic_error(rejected->reason);
    if (pressureMode_ == PressureMode::Incremental)
        throw std::logic_error("cannot replace incremental-pressure state with a legacy Hydro cache");
    // The unversioned format encodes legacy pressure semantics only.
    pressureMode_ = PressureMode::Nonincremental;
    cout << "[Hydrodynamics] Reading flow-field cache..." << endl;

    ostringstream fileName;
    fileName << filePath << "hydrodynamics_cache.bin";

    ifstream cache(fileName.str(), ios::binary | ios::in);
    if (cache.fail())
    {
        cout << "Reading cache \"" << fileName.str() << "\" failed!" << endl;
        exit(-1);
    }

    vertRho.initialize(mesh, 0.0);
    vertMu.initialize(mesh, 0.0);
    for (int i = 0; i < mesh.vertIndexNum(); i++)
    {
        cache.read(reinterpret_cast<char*>(&vertRho.value[i]), sizeof(double));
        cache.read(reinterpret_cast<char*>(&vertMu.value[i]), sizeof(double));
    }
    ux.initialize(mesh, 0.0);
    uxPre.initialize(mesh, 0.0);
    uxMid.initialize(mesh, 0.0);
    for (int i = 0; i < mesh.vEdgeIndexNum(); i++)
    {
        cache.read(reinterpret_cast<char*>(&ux.value[i]), sizeof(double));
        cache.read(reinterpret_cast<char*>(&uxPre.value[i]), sizeof(double));
        cache.read(reinterpret_cast<char*>(&uxMid.value[i]), sizeof(double));
    }
    uy.initialize(mesh, 0.0);
    uyPre.initialize(mesh, 0.0);
    uyMid.initialize(mesh, 0.0);
    for (int i = 0; i < mesh.hEdgeIndexNum(); i++)
    {
        cache.read(reinterpret_cast<char*>(&uy.value[i]), sizeof(double));
        cache.read(reinterpret_cast<char*>(&uyPre.value[i]), sizeof(double));
        cache.read(reinterpret_cast<char*>(&uyMid.value[i]), sizeof(double));
    }
    p.initialize(mesh, 0.0);
    rho.initialize(mesh, 0.0);
    mu.initialize(mesh, 0.0);
    externalFx.initialize(mesh, 0.0);
    externalFy.initialize(mesh, 0.0);
    directFx.initialize(mesh, 0.0);
    directFy.initialize(mesh, 0.0);
    evapUx.initialize(mesh, 0.0);
    evapUy.initialize(mesh, 0.0);
    uMag.initialize(mesh, 0.0);
    divU.initialize(mesh, 0.0);
    magGradU.initialize(mesh, 0.0);
    for (int i = 0; i < mesh.cellIndexNum(); i++)
    {
        cache.read(reinterpret_cast<char*>(&p.value[i]), sizeof(double));
        cache.read(reinterpret_cast<char*>(&rho.value[i]), sizeof(double));
        cache.read(reinterpret_cast<char*>(&mu.value[i]), sizeof(double));
        cache.read(reinterpret_cast<char*>(&externalFx.value[i]), sizeof(double));
        cache.read(reinterpret_cast<char*>(&externalFy.value[i]), sizeof(double));
        cache.read(reinterpret_cast<char*>(&directFx.value[i]), sizeof(double));
        cache.read(reinterpret_cast<char*>(&directFy.value[i]), sizeof(double));
        cache.read(reinterpret_cast<char*>(&evapUx.value[i]), sizeof(double));
        cache.read(reinterpret_cast<char*>(&evapUy.value[i]), sizeof(double));
        cache.read(reinterpret_cast<char*>(&uMag.value[i]), sizeof(double));
        cache.read(reinterpret_cast<char*>(&divU.value[i]), sizeof(double));
        cache.read(reinterpret_cast<char*>(&magGradU.value[i]), sizeof(double));
    }

    cache.read(reinterpret_cast<char*>(&maxDivURes_), sizeof(double));
    cache.read(reinterpret_cast<char*>(&uMagMax_), sizeof(double));

    cache.close();
    cout << "[Hydrodynamics] Flow-field cache read." << endl;
}
