#pragma once

#include "Modules/Hydrodynamics/HydroCoupledBalance.h"

// Close the *assembled* CN momentum equation and FV continuity constraint.
// The existing predictor and Poisson projection remain separately observable;
// this correction removes their viscous splitting defect using the same A,
// pressure face gradient and geometric fluxes. No constant-mu identity,
// analytic initial pressure, dense Schur complement or physical time substep.
template<typename Mesh>
void Hydrodynamics::correctCoupledBalance(Mesh& mesh, bool axisymmetric,
    const Eigen::SparseMatrix<double>& momentum, const Eigen::VectorXd& momentumRhs)
{
    const int cells = mesh.cellIndexNum();
    const int nxFaces = mesh.vEdgeIndexNum(), nyFaces = mesh.hEdgeIndexNum();
    const int velocities = nxFaces + nyFaces;
    // The native nonuniform strategy retains its existing split projection.
    // Still measure A*U + C*pi - rhs against the actual assembled equation:
    // a small predictor/Poisson residual does not bound this splitting defect.
    const auto reportProjectedBalance = [&] {
        activeReport_.stage = "projected-momentum-balance";
        Eigen::VectorXd velocity(velocities);
        Eigen::VectorXd pressureLoad = Eigen::VectorXd::Zero(velocities);
        std::vector<double> dualVolumes(velocities, 0.0);
        const auto sample = [&](const Edge& edge, bool horizontal, int row, double value) {
            velocity[row] = value;
            // This strategy has closed prescribed-zero normal velocities.
            // Their rows are constraints, not wall-reaction measurements.
            if (edge.masterCell == -1 || edge.slaveCell == -1) return;
            const auto face = pressureFaceData(mesh, edge, horizontal, pressurePotential_);
            const double measure = face.dualArea
                * (axisymmetric ? 2.0 * Pi * edge.midpointCoord().x : 1.0);
            dualVolumes[row] = measure;
            pressureLoad[row] = measure * face.gradient;
        };
        for (int i = 0; i < nxFaces; ++i)
            sample(mesh.vEdges(mesh.vEdgeIndexToLabel[i]), false, i, ux.value[i]);
        for (int i = 0; i < nyFaces; ++i)
            sample(mesh.hEdges(mesh.hEdgeIndexToLabel[i]), true, i + nxFaces, uy.value[i]);
        const auto balance = HydroCoupledBalance::momentumResidual(momentum, velocity,
            pressureLoad, momentumRhs, [&](int row) { return dualVolumes[row]; });
        activeReport_.momentumBalanceRelativeResidual = balance.relative;
        activeReport_.momentumBalanceMaxForceDensity = balance.maxForceDensity;
        // Diagnostic only: keep the original projection, gates and histories.
        // No coupled matrix, factorization, velocity or pressure update here.
    };
    // Coupled correction remains restricted to its existing uniform topology.
    if (mesh.Nx <= 0 || mesh.Ny <= 0
        || static_cast<long long>(mesh.Nx) * mesh.Ny != cells) { reportProjectedBalance(); return; }
    const auto spacing = mesh.cells(mesh.cellIndexToLabel[0]).size();
    for (int i = 0; i < cells; ++i)
    {
        const auto& cell = mesh.cells(mesh.cellIndexToLabel[i]);
        if (cell.level() != mesh.cells(mesh.cellIndexToLabel[0]).level()
            || std::abs(cell.size().x - spacing.x) > 1e-12 * spacing.x
            || std::abs(cell.size().y - spacing.y) > 1e-12 * spacing.y
            || cell.cellToVEdge[0][1] != -1 || cell.cellToVEdge[1][1] != -1
            || cell.cellToHEdge[0][1] != -1 || cell.cellToHEdge[1][1] != -1) { reportProjectedBalance(); return; }
    }
    activeReport_.stage = "coupled-assembly";
    activeReport_.coupledCorrectionUsed = true;
    if (momentum.rows() != velocities || momentum.cols() != velocities
        || momentumRhs.size() != velocities
        || cells > std::numeric_limits<int>::max() - velocities)
        throw std::runtime_error("invalid coupled momentum dimensions");

    Eigen::VectorXd projected(velocities);
    for (int i = 0; i < nxFaces; ++i) projected[i] = ux.value[i];
    for (int i = 0; i < nyFaces; ++i) projected[i + nxFaces] = uy.value[i];
    if (!projected.allFinite()) throw std::runtime_error("non-finite projected velocity before coupled correction");
    Eigen::VectorXd increment(cells);
    for (int i = 0; i < cells; ++i) increment[i] = pressurePotential_.value[i];
    Eigen::VectorXd pressureLoad = Eigen::VectorXd::Zero(velocities);
    Eigen::VectorXd continuityDefect = Eigen::VectorXd::Zero(cells);
    HydroCoupledBalance::Assembly assembly(momentum, cells, 4ULL * velocities,
        boundaryFaces_.empty() ? std::optional<int>(0) : std::nullopt);

    auto assemble = [&](const Edge& edge, bool horizontal, int row)
    {
        if (edge.masterCell == -1 || edge.slaveCell == -1)
        {
            if (!boundaryFaces_.empty())
            {
                const auto& boundary = boundaryFaces_[row];
                const double outward = boundary.sign * boundary.area;
                // Keep every continuity row: external traction anchors p.
                assembly.addDivergence(boundary.owner, row, outward);
                continuityDefect[boundary.owner] += outward * projected[row];
                if (boundary.outlet)
                {
                    assembly.addPressureGradient(row, boundary.owner, -outward);
                    pressureLoad[row] = -outward * increment[boundary.owner];
                }
            }
            return;
        }
        const auto& master = mesh.cells(edge.masterCell);
        const auto& slave = mesh.cells(edge.slaveCell);
        const int m = master.index(), s = slave.index();
        const auto face = pressureFaceData(mesh, edge, horizontal, pressurePotential_);
        const double radialMeasure = axisymmetric ? 2.0 * Pi * edge.midpointCoord().x : 1.0;
        const double dualVolume = face.dualArea * radialMeasure;
        const double gradientLoad = dualVolume / (horizontal ? spacing.y : spacing.x);
        const double fluxArea = edge.length() * radialMeasure;
        if (m < 0 || m >= cells || s < 0 || s >= cells || !std::isfinite(gradientLoad)
            || !(gradientLoad > 0.0) || !std::isfinite(fluxArea) || !(fluxArea > 0.0))
            throw std::runtime_error("invalid coupled pressure/continuity face metric");
        // C = Vdual G, exactly the historical pressure load. Master is the
        // positive-coordinate cell. Boundary normal velocity rows get no C.
        assembly.addPressureGradient(row, m, gradientLoad);
        assembly.addPressureGradient(row, s, -gradientLoad);
        pressureLoad[row] = gradientLoad * (increment[m] - increment[s]);
        // Lower block is -D (integrated flux, no density): a harmless sign
        // change makes the interior off-diagonal blocks adjoints on this grid.
        // Replace only cell 0's redundant continuity equation by delta_p0=0.
        if (m != 0 || !boundaryFaces_.empty())
        {
            assembly.addDivergence(m, row, -fluxArea);
            continuityDefect[m] -= fluxArea * projected[row];
        }
        if (s != 0 || !boundaryFaces_.empty())
        {
            assembly.addDivergence(s, row, fluxArea);
            continuityDefect[s] += fluxArea * projected[row];
        }
    };
    for (int i = 0; i < nxFaces; ++i)
        assemble(mesh.vEdges(mesh.vEdgeIndexToLabel[i]), false, i);
    for (int i = 0; i < nyFaces; ++i)
        assemble(mesh.hEdges(mesh.hEdgeIndexToLabel[i]), true, i + nxFaces);
    const auto coupled = assembly.finish();
    const Eigen::VectorXd coupledRhs = coupled.rightHandSide(momentumRhs, projected, pressureLoad, continuityDefect);
    const Eigen::VectorXd correction = checkedLinearSolve(coupled.matrix(), coupledRhs, "coupled",
        activeReport_.coupledRelativeResidual);

    activeReport_.stage = "coupled-balance-check";
    Eigen::VectorXd finalVelocity = projected + correction.head(velocities);
    // Retain the existing exact homogeneous normal-velocity boundary contract,
    // including r=0. Validate the resulting velocity, not a pre-clamp residual.
    for (int i = 0; i < nxFaces; ++i)
    {
        const auto& edge = mesh.vEdges(mesh.vEdgeIndexToLabel[i]);
        if (edge.masterCell == -1 || edge.slaveCell == -1)
        {
            if (boundaryFaces_.empty()) finalVelocity[i] = 0.0;
            else if (!boundaryFaces_[i].outlet)
                finalVelocity[i] = boundaryFaces_[i].sign * boundaryFaces_[i].value;
        }
    }
    for (int i = 0; i < nyFaces; ++i)
    {
        const auto& edge = mesh.hEdges(mesh.hEdgeIndexToLabel[i]);
        if (edge.masterCell == -1 || edge.slaveCell == -1)
        {
            if (boundaryFaces_.empty()) finalVelocity[i + nxFaces] = 0.0;
            else if (!boundaryFaces_[i + nxFaces].outlet)
                finalVelocity[i + nxFaces] = boundaryFaces_[i + nxFaces].sign * boundaryFaces_[i + nxFaces].value;
        }
    }
    increment += correction.tail(cells);
    pressureLoad.setZero();
    std::vector<double> dualVolumes(velocities, 0.0);
    auto finalLoad = [&](const Edge& edge, bool horizontal, int row)
    {
        if (edge.masterCell == -1 || edge.slaveCell == -1)
        {
            if (!boundaryFaces_.empty() && boundaryFaces_[row].outlet)
            {
                const auto& boundary = boundaryFaces_[row];
                dualVolumes[row] = boundary.volume;
                pressureLoad[row] = -boundary.sign * boundary.area * increment[boundary.owner];
            }
            return;
        }
        const auto face = pressureFaceData(mesh, edge, horizontal, pressurePotential_);
        const double volume = face.dualArea * (axisymmetric ? 2.0 * Pi * edge.midpointCoord().x : 1.0);
        const double gradient = (increment[mesh.cells(edge.masterCell).index()]
            - increment[mesh.cells(edge.slaveCell).index()]) / (horizontal ? spacing.y : spacing.x);
        dualVolumes[row] = volume;
        pressureLoad[row] = volume * gradient;
    };
    for (int i = 0; i < nxFaces; ++i)
        finalLoad(mesh.vEdges(mesh.vEdgeIndexToLabel[i]), false, i);
    for (int i = 0; i < nyFaces; ++i)
        finalLoad(mesh.hEdges(mesh.hEdgeIndexToLabel[i]), true, i + nxFaces);
    const auto balance = HydroCoupledBalance::momentumResidual(momentum, finalVelocity,
        pressureLoad, momentumRhs, [&](int row) { return dualVolumes[row]; });
    activeReport_.momentumBalanceRelativeResidual = balance.relative;
    activeReport_.momentumBalanceMaxForceDensity = balance.maxForceDensity;
    if (!increment.allFinite() || balance.relative > activeOptions_.completeConservation->momentumBalance)
        throw std::runtime_error("corrected unsplit momentum balance exceeds momentumBalance or is non-finite");
    // No velocity history, AB2 forcing, material or geometry changes occur
    // inside this correction. Retain the original phi/predictor diagnostics.
    for (int i = 0; i < nxFaces; ++i) ux.value[i] = finalVelocity[i];
    for (int i = 0; i < nyFaces; ++i) uy.value[i] = finalVelocity[i + nxFaces];
    for (int i = 0; i < cells; ++i)
    {
        p.value[i] += correction[velocities + i];
        if (!std::isfinite(p.value[i])) throw std::runtime_error("non-finite coupled physical pressure update");
    }
    if (boundaryFaces_.empty())
    {
        const double gauge = p.value[0];
        for (int i = 0; i < cells; ++i)
        {
            p.value[i] -= gauge;
            if (!std::isfinite(p.value[i])) throw std::runtime_error("non-finite coupled pressure gauge subtraction");
        }
    }
}
