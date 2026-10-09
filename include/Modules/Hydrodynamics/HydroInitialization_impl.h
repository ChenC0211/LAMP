#pragma once

#include <type_traits>
#include "Core/Mesh/NativeMeshIdentity.h"
#include "Modules/Hydrodynamics/HydroContinuityBalance.h"

// Zero-time velocity compatibility, not a Navier--Stokes predictor. The same
// boundary metrics and Schur operator as the ordinary projection are used with
// inverseMassScale=1: min 1/2 sum_f rho_f V_f (u_f-u_hat_f)^2, B u=0.
template<typename Mesh>
HydroVelocityInitializationReport Hydrodynamics::initializeCompatibleVelocity(Mesh& mesh,
    const HydroSolveOptions& options, bool axisymmetric)
{
    HydroVelocityInitializationReport report;
    try
    {
        if (options.fittedNumerics)
            throw std::invalid_argument("compatible velocity initialization does not consume fittedNumerics");
        const auto conservation = options.completeConservation.value_or(HydroCompleteConservation{});
        conservation.validate();
        if (stepProgress_.stopped || stepEntry_ == StepEntry::Legacy
            || stepProgress_.completedSteps != 0 || fixedDt_ != 0.0
            || pressureMode_ != PressureMode::Unset || boundaryConfigurationFixed_)
            throw std::invalid_argument("compatible velocity initialization must precede the first Hydro time step");
        if (initialStateSet_ || materialState_)
        {
            if (stepProgress_.geometry != HydroStepGeometry::complete(axisymmetric).kind)
                throw std::invalid_argument("initial projection must use the initialized complete-domain capability");
            if (!materialState_)
                throw std::invalid_argument("initial projection requires initialized materials");
            if constexpr (std::is_same_v<Mesh, QuadtreeMesh>)
                if (stepMeshIdentity_) stepMeshIdentity_->validate(mesh);
        }
        if (!options.closedDomainProjection
            || options.surfaceForces || options.evaporation || options.directForce)
            throw std::invalid_argument("compatible velocity initialization requires checked complete-domain projection");
        if (!std::isfinite(options.linearResidualTolerance) || options.linearResidualTolerance <= 0.0)
            throw std::invalid_argument("linearResidualTolerance must be finite and positive");
        HydroBoundaryDetail::uniform(mesh);
        const int nx = mesh.vEdgeIndexNum(), ny = mesh.hEdgeIndexNum(), cells = mesh.cellIndexNum();
        const auto validField = [](const auto& field, int size, bool positive) {
            if (field.size() != size || !field.value) return false;
            for (int i = 0; i < size; ++i)
                if (!std::isfinite(field.value[i]) || (positive && !(field.value[i] > 0.0))) return false;
            return true;
        };
        if (!validField(rho, cells, true) || !validField(ux, nx, false) || !validField(uy, ny, false)
            || !validField(uxPre, nx, false) || !validField(uyPre, ny, false)
            || !validField(uxMid, nx, false) || !validField(uyMid, ny, false)
            || !validField(divU, cells, false) || !validField(uMag, cells, false))
            throw std::invalid_argument("initial velocity, histories, diagnostics and positive density must match the mesh and be finite");
        for (int i = 0; i < cells; ++i)
        {
            const auto& cell = mesh.cells(mesh.cellIndexToLabel[i]);
            if (cell.index() != i || !(cell.size().x > 0.0) || !(cell.size().y > 0.0)
                || !std::isfinite(cell.size().x) || !std::isfinite(cell.size().y)
                || !(cell.area() > 0.0) || !std::isfinite(cell.area())
                || !std::isfinite(cell.centerCoord().x) || !std::isfinite(cell.centerCoord().y)
                || (axisymmetric && !(cell.centerCoord().x > 0.0)))
                throw std::invalid_argument("invalid initial projection cell metric or indexing");
        }
        if (axisymmetric && !(mesh.minCoord.x >= 0.0))
            throw std::invalid_argument("RZ initial projection requires a nonnegative radial domain");
        const auto boundaries = makeBoundaryFaces(mesh, options.boundaryPatches, axisymmetric, true);
        VEdgeField<double> nextX(mesh, 0.0);
        HEdgeField<double> nextY(mesh, 0.0);
        for (int i = 0; i < nx; ++i) nextX.value[i] = ux.value[i];
        for (int i = 0; i < ny; ++i) nextY.value[i] = uy.value[i];
        for (int row = 0; row < nx + ny; ++row)
        {
            const auto& face = boundaries[row];
            if (face.owner < 0 || face.outlet) continue;
            const double component = face.sign * face.value;
            if (row < nx) nextX.value[row] = component;
            else nextY.value[row - nx] = component;
        }
        const auto before = velocityDiagnostics(mesh, nextX, nextY, axisymmetric, false);
        report.divergenceBeforeProjection = before.first;
        report.stage = "projection-assembly";
        const auto projection = assembleBoundaryProjection(mesh, 1.0, axisymmetric, boundaries, nextX, nextY);
        const auto& system = projection.equations();
        // Closed-box multiplier uses the same cell-zero gauge as complete NS.
        // Keep the original ungauged RHS for the EVERY-cell continuity check,
        // including the replaced row; the multiplier is never physical pressure.
        HydroPressureSystem gauged;
        if(options.boundaryPatches.empty()) {
            gauged=system;
            gauged.matrix.prune([](int row,int,double){return row!=0;});
            gauged.matrix.coeffRef(0,0)=1.;gauged.matrix.makeCompressed();gauged.rhs[0]=0;
        }
        const auto& solvedSystem=options.boundaryPatches.empty()?gauged:system;
        for (int row = 0; row < nx + ny; ++row)
        {
            const auto& face = boundaries[row];
            if (!std::isfinite(projection.face(row).area) || projection.face(row).area < 0.0
                || ((face.owner < 0 || face.outlet) && (!(projection.face(row).area > 0.0)
                    || !(projection.face(row).inverseMass > 0.0) || !std::isfinite(projection.face(row).inverseMass))))
                throw std::runtime_error("unrepresentable initial projection face metric or inverse mass");
        }
        Eigen::VectorXd multiplier;
        {
            // checkedLinearSolve uses the existing SuperLU backend. An unknown
            // stage gets a one-shot factor, not a physical pressure cache.
            // Restore even on a failed solve; no initialization state leaks
            // into the NS report, options, pressure, dt or history locks.
            struct RestoreSolveContext
            {
                HydroSolveReport& target;
                HydroSolveReport saved;
                double& tolerance;
                double previousTolerance;
                RestoreSolveContext(HydroSolveReport& r, double& t, double next)
                    : target(r), tolerance(t), previousTolerance(t)
                { std::swap(target, saved); tolerance = next; }
                ~RestoreSolveContext()
                { std::swap(target, saved); tolerance = previousTolerance; }
            } restore(activeReport_, activeOptions_.linearResidualTolerance, options.linearResidualTolerance);
            try
            {
                multiplier = checkedLinearSolve(solvedSystem.matrix, solvedSystem.rhs, "initialization",
                    report.linearRelativeResidual);
                report.usedSuperLU = true;
            }
            catch (...)
            {
                report.stage = activeReport_.stage;
                throw;
            }
        }
        report.stage = "projection-validation";
        projection.apply(multiplier,
            [&](std::size_t row) { return row < nx ? nextX.value[row] : nextY.value[row - nx]; },
            [&](std::size_t row, double value) { if (row < nx) nextX.value[row] = value; else nextY.value[row - nx] = value; },
            true, "nonfinite compatible initial velocity");
        CellField<double> nextDivergence, nextMagnitude;
        const auto after = velocityDiagnostics(mesh, nextX, nextY, axisymmetric, true,
            &nextDivergence, &nextMagnitude);
        report.divergenceAfterProjection = after.first;
        report.maxSpeed = after.second;
        // Check actual FV fluxes after forming velocities, not merely the
        // multiplier residual. Use dimensional integrated-flux roundoff scales.
        const auto flux = projection.fluxes(
            [&](std::size_t row) { return row < nx ? nextX.value[row] : nextY.value[row - nx]; });
        if (!HydroContinuityBalance::compatibleProjection(flux.integrated, system.rhs, flux.absolute,
            conservation.initialContinuity))
            throw std::runtime_error("compatible initial velocity does not satisfy every-cell FV continuity");
        const double outletTolerance = (conservation.outletAdmissibility
            + 128.0 * std::numeric_limits<double>::epsilon()) * std::max(before.second, after.second);
        for (int row = 0; row < nx + ny; ++row)
            if (boundaries[row].outlet)
            {
                const double value = row < nx ? nextX.value[row] : nextY.value[row - nx];
                if (boundaries[row].sign * value < -outletTolerance)
                    throw std::runtime_error("compatible initial projection requires unsupported pressure-outlet backflow");
            }
        // All allocations/checks finish before publication. Keep existing field
        // allocations and labels, and never install the multiplier as pressure.
        report.stage = "complete";
        report.success = true;
        for (int i = 0; i < nx; ++i)
            ux.value[i] = uxPre.value[i] = uxMid.value[i] = nextX.value[i];
        for (int i = 0; i < ny; ++i)
            uy.value[i] = uyPre.value[i] = uyMid.value[i] = nextY.value[i];
        for (int i = 0; i < cells; ++i)
        {
            divU.value[i] = nextDivergence.value[i];
            uMag.value[i] = nextMagnitude.value[i];
        }
        maxDivURes_ = after.first;
        uMagMax_ = after.second;
    }
    catch (const std::exception& error)
    {
        report.reason = error.what();
    }
    return report;
}
