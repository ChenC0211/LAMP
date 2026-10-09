#include "Modules/Hydrodynamics/Hydrodynamics.h"
#include "Modules/Hydrodynamics/HydroFittedDomain.h"
#include "Modules/Hydrodynamics/HydroFittedMomentum.h"
#include "Modules/Hydrodynamics/HydroFittedProjection.h"
#include "Modules/Hydrodynamics/HydroCoupledBalance.h"
#include "Modules/Hydrodynamics/HydroContinuityBalance.h"
#include "Core/Mesh/QuadtreeMesh.h"

#include <limits>
#include <iomanip>
#include <sstream>

namespace
{
using Fitted = HydroFittedDomain;

template<class Field, class LabelAt>
void checkField(const Field& field, int count, LabelAt labelAt)
{
    if (field.size() != count || !field.value) throw std::invalid_argument("fitted Hydro field does not match mesh");
    for (int i = 0; i < count; ++i)
        if (field.label(i) != labelAt(i) || !std::isfinite(field.value[i]))
            throw std::invalid_argument("fitted Hydro field has stale labels or non-finite values");
}

Eigen::VectorXd compact(const Fitted& domain, const Eigen::VectorXd& native)
{
    Eigen::VectorXd result(domain.nativeByVelocityRow().size());
    for (int i = 0; i < result.size(); ++i) result[i] = native[domain.nativeByVelocityRow()[i]];
    return result;
}

double outletTolerance(const Fitted::NativeFace& face, double velocity, double old,
    double previous, double predictor, double pressure, double density, double dt, double relativeTolerance)
{
    const double pressureScale = std::abs(pressure) + std::abs(face.boundary.pressure);
    const double impulse = dt * pressureScale * face.area / (density * face.dualVolume);
    const double advectiveScale = std::max({std::abs(velocity), std::abs(old), std::abs(previous)});
    // Uncancelled predictor/pressure impulses affect roundoff, not the allowed
    // physical backflow. Use the actual RZ half-dual area/volume ratio.
    return relativeTolerance * advectiveScale
        + 64.0 * std::numeric_limits<double>::epsilon() * (impulse + std::abs(predictor));
}

void checkOutlet(const Fitted& domain, const Eigen::VectorXd& velocity,
    const Eigen::VectorXd& old, const Eigen::VectorXd& previous,
    const Eigen::VectorXd& predictor, const Eigen::VectorXd& pressure,
    const std::vector<double>& density, double dt, double relativeTolerance)
{
    for (const auto& face : domain.faces())
    {
        if (face.boundary.type != Fitted::BoundaryType::PressureOutlet) continue;
        const int row = face.velocityRow;
        const int owner = face.negativeCell >= 0 ? face.negativeCell : face.positiveCell;
        const double tolerance = outletTolerance(face, velocity[row], old[row], previous[row], predictor[row],
            pressure[domain.cellRowByIndex()[owner]], density[row], dt, relativeTolerance);
        if (!std::isfinite(tolerance) || face.outwardSign * velocity[row] < -tolerance)
        {
            std::ostringstream message;
            message << std::setprecision(17)
                << "fitted pressure outlet has backflow without incoming material/momentum data: boundary "
                << face.boundary.id << ", nativeRow=" << face.nativeRow
                << ", outwardVelocity_code=" << face.outwardSign * velocity[row]
                << ", roundoffTolerance_code=" << tolerance;
            throw std::runtime_error(message.str());
        }
    }
}
}

HydroSolveReport Hydrodynamics::solveFittedStep(QuadtreeMesh& mesh, double dt,
    const HydroSolveOptions& options, const HydroFittedDomain& domain,
    const CellField<double>* forceX, const CellField<double>* forceY,
    const VEdgeField<double>* faceForceX, const HEdgeField<double>* faceForceY)
{
    if (const auto rejected = claimStepEntry(StepEntry::Managed)) return *rejected;
    activeOptions_ = options;
    activeReport_ = HydroSolveReport{};
    try
    {
        if (options.completeConservation)
            throw std::invalid_argument("fitted Hydro does not consume completeConservation");
        if (options.compatibleNativePressure)
            throw std::invalid_argument("compatible native pressure is not implemented on a fitted domain");
        const auto numerics = options.fittedNumerics.value_or(HydroFittedNumerics(options.linearResidualTolerance));
        numerics.validate();
        if (!domain.axisymmetric() || !options.closedDomainProjection || !options.boundaryPatches.empty()
            || options.surfaceForces || options.directForce || options.evaporation)
            throw std::invalid_argument("fitted RZ requires explicit domain, incremental pressure and supplied forces, not legacy coupling flags");
        if (!std::isfinite(dt) || dt <= 0.0
            || !std::isfinite(options.accelerationX) || !std::isfinite(options.accelerationY))
            throw std::invalid_argument("invalid fitted dt, residual tolerance or body acceleration");
        if (pressureMode_ == PressureMode::Nonincremental || (fixedDt_ > 0.0 && !fittedDomain_))
            throw std::invalid_argument("cannot reuse complete-domain/legacy Hydro histories in a fitted-fluid step");
        if (fixedDt_ > 0.0 && std::abs(dt - fixedDt_) > 1e-12 * fixedDt_)
            throw std::invalid_argument("fitted CN/AB2 Hydro requires a fixed dt");
        domain.validateMeshIdentity(mesh);
        if (fittedDomain_) fittedDomain_->validateMeshIdentity(mesh);
        if (fittedDomain_ && !fittedDomain_->sameDefinition(domain))
            throw std::invalid_argument("fitted fluid topology/boundaries cannot change after an accepted step");
        if ((forceX == nullptr) != (forceY == nullptr) || (faceForceX == nullptr) != (faceForceY == nullptr))
            throw std::invalid_argument("fitted forcing requires both components");
        const int cellCount = mesh.cellIndexNum(), vxCount = mesh.vEdgeIndexNum(), vyCount = mesh.hEdgeIndexNum();
        const auto cellLabel = [&](int i) { return mesh.cellIndexToLabel[i]; };
        const auto vertexLabel = [&](int i) { return mesh.vertIndexToLabel[i]; };
        const auto vxLabel = [&](int i) { return mesh.vEdgeIndexToLabel[i]; };
        const auto vyLabel = [&](int i) { return mesh.hEdgeIndexToLabel[i]; };
        for (const auto* field : {&rho, &mu, &p, &pressurePotential_, &externalFx, &externalFy, &divU, &uMag})
            checkField(*field, cellCount, cellLabel);
        checkField(vertRho, mesh.vertIndexNum(), vertexLabel);
        checkField(vertMu, mesh.vertIndexNum(), vertexLabel);
        for (const auto* field : {&ux, &uxPre, &uxMid}) checkField(*field, vxCount, vxLabel);
        for (const auto* field : {&uy, &uyPre, &uyMid}) checkField(*field, vyCount, vyLabel);
        if (forceX) { checkField(*forceX, cellCount, cellLabel); checkField(*forceY, cellCount, cellLabel); }
        if (faceForceX) { checkField(*faceForceX, vxCount, vxLabel); checkField(*faceForceY, vyCount, vyLabel); }

        Eigen::VectorXd native(vxCount + vyCount), previousNative(vxCount + vyCount);
        for (int i = 0; i < vxCount + vyCount; ++i)
        {
            native[i] = i < vxCount ? ux.value[i] : uy.value[i - vxCount];
            previousNative[i] = i < vxCount ? uxPre.value[i] : uyPre.value[i - vxCount];
            if (domain.velocityRowByNative()[i] < 0 && (native[i] != 0.0 || previousNative[i] != 0.0))
                throw std::invalid_argument("solid or zero-normal fitted face has nonzero velocity/history");
        }
        CellField<double> bodyX(mesh, 0.0), bodyY(mesh, 0.0);
        Eigen::VectorXd oldPressure(domain.cellIndexByRow().size());
        for (int i = 0; i < oldPressure.size(); ++i)
        {
            const int cell = domain.cellIndexByRow()[i];
            if (!(rho.value[cell] > 0.0) || mu.value[cell] < 0.0)
                throw std::invalid_argument("active fluid requires positive density and nonnegative viscosity");
            oldPressure[i] = p.value[cell];
            bodyX.value[cell] = rho.value[cell] * options.accelerationX + (forceX ? forceX->value[cell] : 0.0);
            bodyY.value[cell] = rho.value[cell] * options.accelerationY + (forceY ? forceY->value[cell] : 0.0);
        }
        activeReport_.stage = "fitted-momentum-assembly";
        const Eigen::VectorXd oldVelocity = compact(domain, native);
        const Eigen::VectorXd previousVelocity = compact(domain, previousNative);
        std::vector<double> outletTolerances(domain.nativeByVelocityRow().size(), 0.0);
        for (const auto& face : domain.faces())
        {
            if (face.boundary.type != Fitted::BoundaryType::PressureOutlet) continue;
            const int owner = face.negativeCell >= 0 ? face.negativeCell : face.positiveCell;
            const int row = face.velocityRow;
            outletTolerances[row] = outletTolerance(face, oldVelocity[row], oldVelocity[row], previousVelocity[row],
                oldVelocity[row], oldPressure[domain.cellRowByIndex()[owner]], rho.value[owner], dt,
                numerics.outletAdmissibility);
        }
        HydroFittedNativeInterior nativeInterior;
        auto momentumLayout = fittedMomentumLayout_;
        if (!domain.uniformLeaves())
        {
            if (!momentumLayout) momentumLayout = HydroFittedMomentumLayout::build(mesh, domain);
            momentumLayout->validate(mesh, domain);
            nativeInterior.layout = momentumLayout;
            const auto& enabled = momentumLayout->nativeInteriorRows();
            searchSupCoeffs(mesh);
            HydroMomentumAssembly assembly;
            buildNativeMomentumMatrix<true>(mesh,dt,assembly,&enabled);
            buildNativeMomentumRHS<true>(mesh,dt,assembly,&enabled,&bodyX,&bodyY);
            nativeInterior.rhs=std::move(assembly).finish(nativeInterior.matrix);
            nativeInterior.density.resize(native.size());
            nativeInterior.viscosity.resize(native.size());
            const auto centers=[&](const Cell& cell,int child) -> const supCoeff& {
                return supSonCells_[cell.supCellIndex*4+child].center;
            };
            for (int row:domain.nativeByVelocityRow()) if (enabled[row])
            {
                const auto& face=domain.faces()[row];
                const auto& edge=face.horizontal?mesh.hEdges(face.edgeLabel):mesh.vEdges(face.edgeLabel);
                nativeInterior.density[row]=pressureFaceData(mesh,edge,face.horizontal,p).density;
                nativeInterior.viscosity[row]=HydroNativePressureStencil::faceData(mesh,edge,face.horizontal,mu,mu,centers).density;
            }
        }
        auto momentum = assembleHydroFittedMomentum(mesh, domain, dt, rho, mu, vertRho, vertMu,
            native, previousNative, bodyX, bodyY, outletTolerances,domain.uniformLeaves()?nullptr:&nativeInterior);
        activeReport_.fittedRegularMomentumRows = static_cast<int>(momentum.regularNativeRows);
        activeReport_.fittedBoundaryMomentumRows = static_cast<int>(momentum.boundaryPieceRows);
        activeReport_.fittedAdaptiveMomentumRows = static_cast<int>(momentum.adaptiveNativeRows);
        if (faceForceX)
            for (const auto& face : domain.faces())
                if (face.velocityRow >= 0)
                {
                    const double force = face.horizontal ? faceForceY->value[face.nativeRow - vxCount]
                        : faceForceX->value[face.nativeRow];
                    momentum.rhs[face.velocityRow] += face.dualVolume * force;
                }
        checkOutlet(domain, oldVelocity, oldVelocity, previousVelocity, oldVelocity, oldPressure,
            momentum.densityByVelocityRow, dt, numerics.outletAdmissibility);
        // The predictor consumes historical physical p through the SAME C as
        // projection and the coupled residual. No cell-to-face pressure lift.
        const Eigen::VectorXd rhs = momentum.rhs - domain.pressureGradientLoad() * oldPressure;
        const Eigen::VectorXd predictor = checkedLinearSolve(momentum.matrix, rhs, "momentum",
            activeReport_.momentumRelativeResidual, numerics.momentumLinear);
        activeReport_.divergenceBeforeProjection = HydroContinuityBalance::LocalFlux::maximumDivergence(
            domain.divergence() * predictor, domain.cellVolumes());
        const auto projection = HydroFittedProjection::assemble(domain, dt, momentum.densityByVelocityRow, predictor);
        const auto& pressureSystem = projection.equations();
        Eigen::VectorXd increment = checkedLinearSolve(pressureSystem.matrix, pressureSystem.rhs, "pressure",
            activeReport_.pressureRelativeResidual, numerics.pressureLinear);
        const Eigen::VectorXd projected = projection.correctVelocity(predictor, increment);
        const auto coupled = HydroCoupledBalance::assemble(momentum.matrix,
            domain.pressureGradientLoad(), domain.divergence());
        const Eigen::VectorXd coupledRhs = coupled.rightHandSide(rhs, projected, increment);
        activeReport_.coupledCorrectionUsed = true;
        const Eigen::VectorXd correction = checkedLinearSolve(coupled.matrix(), coupledRhs, "coupled",
            activeReport_.coupledRelativeResidual, numerics.coupledLinear);
        const int velocities = static_cast<int>(domain.nativeByVelocityRow().size());
        Eigen::VectorXd finalVelocity = projected + correction.head(velocities);
        increment += correction.tail(increment.size());
        auto continuity = HydroContinuityBalance::LocalFlux(domain.divergence(), finalVelocity);
        bool continuous = continuity.satisfiesScaledFlux(numerics.localContinuity);
        activeReport_.fittedLocalContinuityBeforeRefinement = continuity.relativeResidual();
        activeReport_.fittedLocalContinuityAfterRefinement = activeReport_.fittedLocalContinuityBeforeRefinement;
        // Near static balance, projected + correction may lose digits compared
        // with the tiny final velocity. Refine the SAME CN/continuity equations;
        // do not zero small velocities or relax the physical gate. Already
        // accepted steps retain the original single-solve arithmetic.
        for (int refinement = 0; !continuous && refinement < 3; ++refinement)
        {
            // A, C, D and the pressure-reference policy are fixed for this
            // step. Re-evaluate the residual against the same frozen operator.
            const Eigen::VectorXd residual = coupled.rightHandSide(rhs, finalVelocity, increment);
            double linearResidual = 0.0;
            const bool originalReuse = activeReport_.coupledFactorReused;
            const Eigen::VectorXd delta = checkedLinearSolve(coupled.matrix(), residual, "coupled",
                linearResidual, numerics.coupledLinear);
            activeReport_.coupledFactorReused = originalReuse;
            finalVelocity += delta.head(velocities);
            increment += delta.tail(increment.size());
            ++activeReport_.fittedBalanceRefinements;
            continuity = HydroContinuityBalance::LocalFlux(domain.divergence(), finalVelocity);
            continuous = continuity.satisfiesScaledFlux(numerics.localContinuity);
            activeReport_.fittedLocalContinuityAfterRefinement = continuity.relativeResidual();
        }
        // Existing full momentum, outlet and local-continuity checks below
        // still decide acceptance. Failed refinements publish no state.
        const Eigen::VectorXd finalPressure = oldPressure + increment;
        const Eigen::VectorXd pressureLoad = domain.pressureGradientLoad() * increment;
        const auto balance = HydroCoupledBalance::momentumResidual(momentum.matrix, finalVelocity,
            pressureLoad, rhs, [&](int row) { return domain.faces()[domain.nativeByVelocityRow()[row]].dualVolume; });
        activeReport_.momentumBalanceRelativeResidual = balance.relative;
        activeReport_.momentumBalanceMaxForceDensity = balance.maxForceDensity;
        if (!increment.allFinite() || !finalPressure.allFinite() || balance.relative > numerics.momentumBalance)
            throw std::runtime_error("fitted full CN momentum balance failed");
        activeReport_.divergenceAfterProjection = continuity.maxDivergence(domain.cellVolumes());
        activeReport_.maxSpeed = continuity.maxFaceSpeed();
        activeReport_.stage = "fitted-outlet-admissibility";
        checkOutlet(domain, finalVelocity, oldVelocity, previousVelocity, predictor, finalPressure,
            momentum.densityByVelocityRow, dt, numerics.outletAdmissibility);

        // Precompute cell diagnostics and the frozen domain before touching
        // physical/history fields. Inactive values are not fluid unknowns.
        // The final velocity has not changed since its last residual evaluation.
        // Preserve the original scaled-flux gate, then publish those same fluxes.
        const auto& flux = continuity.netFlux();
        if (!continuous)
            throw std::runtime_error("fitted local fluid continuity exceeds localContinuity tolerance");
        std::vector<double> divergence(cellCount, 0.0), magnitude(cellCount, 0.0);
        Eigen::VectorXd finalNative = Eigen::VectorXd::Zero(vxCount + vyCount);
        Eigen::VectorXd predictorNative = finalNative;
        for (int i = 0; i < velocities; ++i)
        {
            finalNative[domain.nativeByVelocityRow()[i]] = finalVelocity[i];
            predictorNative[domain.nativeByVelocityRow()[i]] = predictor[i];
        }
        for (int i = 0; i < flux.size(); ++i)
        {
            const int index = domain.cellIndexByRow()[i];
            const auto& cell = mesh.cells(mesh.cellIndexToLabel[index]);
            divergence[index] = flux[i] / domain.cellVolumes()[i];
            const double ur = 0.5 * (finalNative[mesh.vEdges(cell.cellToVEdge[0][0]).index()]
                + finalNative[mesh.vEdges(cell.cellToVEdge[1][0]).index()]);
            const double uz = 0.5 * (finalNative[vxCount + mesh.hEdges(cell.cellToHEdge[0][0]).index()]
                + finalNative[vxCount + mesh.hEdges(cell.cellToHEdge[1][0]).index()]);
            magnitude[index] = std::hypot(ur, uz);
            activeReport_.maxSpeed = std::max(activeReport_.maxSpeed, magnitude[index]);
        }
        const auto frozen = fittedDomain_ ? fittedDomain_ : std::make_shared<const HydroFittedDomain>(domain);
        activeReport_.stage = "complete";
        activeReport_.success = true;
        activeReport_.stateChange = HydroStateChange::Advanced;
        for (int i = 0; i < vxCount; ++i)
        {
            uxPre.value[i] = native[i]; uxMid.value[i] = predictorNative[i]; ux.value[i] = finalNative[i];
        }
        for (int i = 0; i < vyCount; ++i)
        {
            uyPre.value[i] = native[vxCount + i]; uyMid.value[i] = predictorNative[vxCount + i]; uy.value[i] = finalNative[vxCount + i];
        }
        for (int i = 0; i < cellCount; ++i)
        {
            const int row = domain.cellRowByIndex()[i];
            p.value[i] = row >= 0 ? finalPressure[row] : 0.0;
            pressurePotential_.value[i] = row >= 0 ? increment[row] : 0.0;
            externalFx.value[i] = bodyX.value[i]; externalFy.value[i] = bodyY.value[i];
            divU.value[i] = divergence[i]; uMag.value[i] = magnitude[i];
        }
        fittedDomain_ = frozen;
        fittedMomentumLayout_ = std::move(momentumLayout);
        fixedDt_ = dt;
        pressureMode_ = PressureMode::Incremental;
        maxDivURes_ = activeReport_.divergenceAfterProjection;
        uMagMax_ = activeReport_.maxSpeed;
    }
    catch (const std::exception& error)
    {
        activeReport_.success = false;
        activeReport_.reason = error.what();
    }
    return activeReport_;
}
