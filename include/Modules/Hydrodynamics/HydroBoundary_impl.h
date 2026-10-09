#pragma once

#include "Modules/Hydrodynamics/HydroBoundary.h"
#include <iomanip>
#include <sstream>

namespace HydroBoundaryDetail
{
template<typename Mesh> void uniform(const Mesh& mesh)
{
    if (mesh.Nx < 2 || mesh.Ny < 2 || mesh.cellIndexNum() != mesh.Nx * mesh.Ny
        || mesh.vEdgeIndexNum() != (mesh.Nx + 1) * mesh.Ny
        || mesh.hEdgeIndexNum() != mesh.Nx * (mesh.Ny + 1))
        throw std::invalid_argument("Hydro boundary/RZ stress requires a complete fixed uniform grid");
    const auto h = mesh.cells(mesh.cellIndexToLabel[0]).size();
    for (int k = 0; k < mesh.cellIndexNum(); ++k)
    {
        const auto& c = mesh.cells(mesh.cellIndexToLabel[k]);
        if (c.level() != 0 || std::abs(c.size().x - h.x) > 1e-12 * h.x
            || std::abs(c.size().y - h.y) > 1e-12 * h.y
            || c.cellToVEdge[0][1] != -1 || c.cellToVEdge[1][1] != -1
            || c.cellToHEdge[0][1] != -1 || c.cellToHEdge[1][1] != -1)
            throw std::invalid_argument("Hydro boundary/RZ stress does not support AMR or nonuniform cells");
    }
}
}

template<typename Mesh>
std::vector<Hydrodynamics::BoundaryFaceData> Hydrodynamics::makeBoundaryFaces(Mesh& mesh,
    const std::vector<HydroBoundaryPatch>& patches, bool axisymmetric, bool includeClosed) const
{
    std::vector<BoundaryFaceData> faces;
    if (patches.empty() && !includeClosed) return faces;
    HydroBoundaryDetail::uniform(mesh);
    validateHydroBoundaryPatches(patches, mesh.minCoord, mesh.maxCoord,
        mesh.cells(mesh.cellIndexToLabel[0]).size(), axisymmetric);
    const int nx = mesh.vEdgeIndexNum(), count = nx + mesh.hEdgeIndexNum();
    faces.resize(count);
    int outletCount = 0;
    for (int row = 0; row < count; ++row)
    {
        const bool horizontal = row >= nx;
        const int k = horizontal ? row - nx : row;
        const auto& edge = horizontal ? mesh.hEdges(mesh.hEdgeIndexToLabel[k]) : mesh.vEdges(mesh.vEdgeIndexToLabel[k]);
        if (edge.masterCell != -1 && edge.slaveCell != -1) continue;
        auto& face = faces[row];
        face.horizontal = horizontal;
        face.sign = edge.masterCell == -1 ? 1.0 : -1.0;
        const int ownerLabel = edge.masterCell == -1 ? edge.slaveCell : edge.masterCell;
        if (ownerLabel < 0) throw std::invalid_argument("boundary face has no owner");
        const auto& owner = mesh.cells(ownerLabel);
        face.owner = owner.index();
        const auto point = edge.midpointCoord();
        const HydroBoundarySide side = horizontal
            ? (face.sign < 0.0 ? HydroBoundarySide::South : HydroBoundarySide::North)
            : (face.sign < 0.0 ? HydroBoundarySide::West : HydroBoundarySide::East);
        const double coordinate = horizontal ? point.x : point.y;
        for (const auto& patch : patches)
            if (patch.side == side && coordinate > patch.start && coordinate < patch.end)
            {
                face.outlet = patch.type == HydroBoundaryType::PressureOutlet;
                face.value = patch.value;
            }
        face.normalSpacing = horizontal ? owner.size().y : owner.size().x;
        face.tangentSpacing = horizontal ? owner.size().x : owner.size().y;
        const int innerSide = face.sign < 0.0 ? 1 : 0;
        face.inner = horizontal ? nx + mesh.hEdges(owner.cellToHEdge[innerSide][0]).index()
            : mesh.vEdges(owner.cellToVEdge[innerSide][0]).index();
        face.area = edge.length() * (axisymmetric ? 2.0 * Pi * point.x : 1.0);
        face.innerArea = edge.length() * (axisymmetric
            ? 2.0 * Pi * (horizontal ? point.x : owner.centerCoord().x) : 1.0);
        face.volume = axisymmetric && !horizontal
            ? Pi * (point.x + owner.centerCoord().x) * std::abs(point.x - owner.centerCoord().x) * edge.length()
            : 0.5 * owner.area() * (axisymmetric ? 2.0 * Pi * point.x : 1.0);
        if (!(face.volume > 0.0) || !std::isfinite(face.volume) || !std::isfinite(face.area)
            || !(face.innerArea > 0.0) || (face.outlet && !(face.area > 0.0)))
            throw std::invalid_argument("invalid boundary half-control-volume metric");
        if (face.outlet) ++outletCount;
    }
    if (outletCount == 0 && !patches.empty()) throw std::invalid_argument("boundary patches select no physical pressure outlet face");
    return faces;
}

template<typename Mesh>
void Hydrodynamics::prepareBoundaryPatches(Mesh& mesh, bool axisymmetric)
{
    boundaryFaces_ = makeBoundaryFaces(mesh, activeOptions_.boundaryPatches, axisymmetric);
    const int nx = mesh.vEdgeIndexNum();
    for (int row = 0; row < static_cast<int>(boundaryFaces_.size()); ++row)
    {
        const auto& face = boundaryFaces_[row];
        if (face.owner < 0 || face.outlet) continue;
        const double component = face.sign * face.value;
        if (face.horizontal) uy.value[row - nx] = uyPre.value[row - nx] = component;
        else ux.value[row] = uxPre.value[row] = component;
    }
}

template<typename Mesh>
void Hydrodynamics::addBoundaryMomentum(Mesh& mesh, double dt, bool axisymmetric, HydroMomentumAssembly& assembly,bool includeAdvection)
{
    if (boundaryFaces_.empty()) return;
    const int nx = mesh.vEdgeIndexNum(), count = nx + mesh.hEdgeIndexNum();
    const auto velocity = [&](int row, bool previous) {
        return row < nx ? (previous ? uxPre.value[row] : ux.value[row])
            : (previous ? uyPre.value[row - nx] : uy.value[row - nx]);
    };
    // Existing boundary rows are identities. Do not replace any interior row.
    assembly.entries.erase(std::remove_if(assembly.entries.begin(), assembly.entries.end(), [&](const auto& entry) {
        return boundaryFaces_[entry.row()].owner >= 0;
    }), assembly.entries.end());
    for (int row = 0; row < count; ++row)
    {
        const auto& face = boundaryFaces_[row];
        if (face.owner < 0) continue;
        if (!face.outlet)
        {
            assembly.entries.emplace_back(row, row, 1.0);
            assembly.rhs[row] = face.sign * face.value;
            continue;
        }
        const int k = face.horizontal ? row - nx : row;
        const auto& edge = face.horizontal ? mesh.hEdges(mesh.hEdgeIndexToLabel[k]) : mesh.vEdges(mesh.vEdgeIndexToLabel[k]);
        const auto& owner = mesh.cells(mesh.cellIndexToLabel[face.owner]);
        const double mass = rho.value[face.owner] * face.volume / dt;
        // Interior normal stress is 2 mu du_n/dn. Its half-CN coefficient is
        // mu*A_inner/h; the outer full traction is prescribed, not u_b=u_i.
        const double viscous = mu.value[face.owner] * face.innerArea / face.normalSpacing;
        assembly.entries.emplace_back(row, row, mass + viscous);
        assembly.entries.emplace_back(row, face.inner, -viscous);
        const auto convection = [&](bool previous) {
            if(!includeAdvection) return 0.;
            const double ub = velocity(row, previous), ui = velocity(face.inner, previous);
            const double qInner = rho.value[face.owner] * 0.5 * (ub + ui) * face.innerArea;
            const double donor = face.sign * qInner > 0.0 ? ui : ub;
            double result = face.sign * (rho.value[face.owner] * ub * ub * face.area - qInner * donor);
            for (int end = 0; end < 2; ++end)
            {
                const auto& vertex = mesh.verts(edge.edgeToVert[end]);
                const int adjacentLabel = face.horizontal ? vertex.vertToHEdge[end] : vertex.vertToVEdge[end];
                const int tangentRow = face.horizontal ? mesh.vEdges(owner.cellToVEdge[end][0]).index()
                    : nx + mesh.hEdges(owner.cellToHEdge[end][0]).index();
                const double area = !axisymmetric ? 0.5 * face.normalSpacing
                    : face.horizontal ? Pi * vertex.coord().x * face.normalSpacing
                    : face.volume / edge.length();
                if (adjacentLabel < 0)
                {
                    // This half-dual side is a real patch on the perpendicular
                    // physical boundary, not an absent/zero-flux side. Use its
                    // native normal velocity, including an inlet/outlet value.
                    const auto& corner = boundaryFaces_[tangentRow];
                    const double sideSign = end ? 1.0 : -1.0;
                    if (corner.owner < 0 || corner.sign != sideSign || corner.horizontal == face.horizontal)
                        throw std::runtime_error("corner half-volume lacks its perpendicular boundary face");
                    const double qOut = sideSign * vertRho.value[vertex.index()]
                        * velocity(tangentRow, previous) * area;
                    double cornerDonor = ub;
                    if (qOut < 0.0)
                    {
                        // The incoming tangential trace on that perpendicular
                        // side obeys tau_nt=0. Estimate d_t u_n from the actual
                        // two inward native boundary faces, not an exterior zero.
                        const auto& cornerEdge = face.horizontal
                            ? mesh.vEdges(mesh.vEdgeIndexToLabel[tangentRow])
                            : mesh.hEdges(mesh.hEdgeIndexToLabel[tangentRow - nx]);
                        const int inwardEnd = face.sign > 0.0 ? 0 : 1;
                        const auto& inwardVertex = mesh.verts(cornerEdge.edgeToVert[inwardEnd]);
                        const int inwardLabel = face.horizontal ? inwardVertex.vertToVEdge[inwardEnd]
                            : inwardVertex.vertToHEdge[inwardEnd];
                        if (inwardLabel < 0)
                            throw std::runtime_error("corner inflow lacks an inward boundary trace stencil");
                        const auto& inwardEdge = face.horizontal ? mesh.vEdges(inwardLabel) : mesh.hEdges(inwardLabel);
                        const int inwardRow = (face.horizontal ? 0 : nx) + inwardEdge.index();
                        if (boundaryFaces_[inwardRow].owner < 0 || boundaryFaces_[inwardRow].sign != sideSign)
                            throw std::runtime_error("corner inflow trace leaves its physical boundary side");
                        const double derivative = face.sign * (velocity(tangentRow, previous)
                            - velocity(inwardRow, previous)) / face.normalSpacing;
                        cornerDonor = HydroBoundaryDetail::zeroShearTrace(cornerDonor,
                            face.tangentSpacing, sideSign, derivative);
                    }
                    result += qOut * cornerDonor;
                    continue;
                }
                const auto& adjacent = face.horizontal ? mesh.hEdges(adjacentLabel) : mesh.vEdges(adjacentLabel);
                const int adjacentRow = (face.horizontal ? nx : 0) + adjacent.index();
                if (boundaryFaces_[adjacentRow].owner < 0)
                    throw std::runtime_error("boundary half-volume has a nonboundary tangential neighbor");
                const double other = velocity(adjacentRow, previous);
                const double derivative = (end ? other - ub : ub - other) / face.tangentSpacing;
                // tau_nt=0 => d_n(out) u_t = -d_t u_n(out).
                const double tangent = HydroBoundaryDetail::zeroShearTrace(velocity(tangentRow, previous),
                    face.normalSpacing, face.sign, derivative);
                const double q = vertRho.value[vertex.index()] * tangent * area;
                const double transported = q > 0.0 ? (end ? ub : other) : (end ? other : ub);
                result += (end ? 1.0 : -1.0) * q * transported;
            }
            return result;
        };
        const double body = face.horizontal ? externalFy.value[face.owner] : externalFx.value[face.owner];
        assembly.rhs[row] = mass * velocity(row, false)
            - viscous * (velocity(row, false) - velocity(face.inner, false))
            - 1.5 * convection(false) + 0.5 * convection(true) + face.volume * body
            - face.sign * face.area * face.value;
    }
    // A wall-adjacent tangential row previously omitted d_n u_t but retained
    // mu*d_t u_n. Remove that remaining shear on BOTH CN sides. Complete its
    // explicit inflow donor with the same zero-shear one-sided boundary trace.
    for (int row = 0; row < count; ++row)
    {
        if (boundaryFaces_[row].owner >= 0) continue;
        const bool horizontal = row >= nx;
        const int k = horizontal ? row - nx : row;
        const auto& edge = horizontal ? mesh.hEdges(mesh.hEdgeIndexToLabel[k]) : mesh.vEdges(mesh.vEdgeIndexToLabel[k]);
        const auto& high = mesh.cells(edge.masterCell);
        const auto& low = mesh.cells(edge.slaveCell);
        for (int end = 0; end < 2; ++end)
        {
            const auto& vertex = mesh.verts(edge.edgeToVert[end]);
            const double coordinate = horizontal ? vertex.coord().x : vertex.coord().y;
            const double wall = horizontal ? (end ? mesh.maxCoord.x : mesh.minCoord.x)
                : (end ? mesh.maxCoord.y : mesh.minCoord.y);
            if (coordinate != wall) continue;
            const double sign = end ? 1.0 : -1.0;
            const double hn = horizontal ? high.size().x : high.size().y;
            const double ht = horizontal ? high.size().y : high.size().x;
            const int highRow = horizontal ? mesh.vEdges(high.cellToVEdge[end][0]).index()
                : nx + mesh.hEdges(high.cellToHEdge[end][0]).index();
            const int lowRow = horizontal ? mesh.vEdges(low.cellToVEdge[end][0]).index()
                : nx + mesh.hEdges(low.cellToHEdge[end][0]).index();
            const double radius = horizontal ? vertex.coord().x : edge.midpointCoord().x;
            const double area = ht * (axisymmetric ? 2.0 * Pi * radius : 1.0);
            const double coefficient = 0.5 * sign * vertMu.value[vertex.index()] * area / ht;
            assembly.entries.emplace_back(row, highRow, coefficient);
            assembly.entries.emplace_back(row, lowRow, -coefficient);
            assembly.rhs[row] -= coefficient * (velocity(highRow, false) - velocity(lowRow, false));
            const auto addedConvection = [&](bool previous) {
                const double a = velocity(lowRow, previous), b = velocity(highRow, previous);
                const double q = sign * vertRho.value[vertex.index()] * 0.5 * (a + b) * area;
                const double trace = HydroBoundaryDetail::zeroShearTrace(velocity(row, previous),
                    hn, sign, b - a, ht);
                return q < 0.0 ? q * trace : 0.0; // Existing outflow donor is already the interior value.
            };
            if(includeAdvection) assembly.rhs[row] += -1.5 * addedConvection(false) + 0.5 * addedConvection(true);
        }
    }
    if (!assembly.rhs.allFinite()) throw std::runtime_error("nonfinite boundary momentum assembly");
}

template<typename Mesh>
void Hydrodynamics::addNativeRZHoop(Mesh& mesh, double dt, HydroMomentumAssembly& assembly)
{
    (void)dt;
    // div(2mu D)_r has -2mu*u_r/r^2 in addition to the existing r-weighted
    // radial/shear face divergence. This was absent from the old checked RZ
    // operator. Add it to all fixed-grid RZ runs, not just open-boundary runs.
    for (int i = 0; i < mesh.vEdgeIndexNum(); ++i)
    {
        const auto& edge = mesh.vEdges(mesh.vEdgeIndexToLabel[i]);
        double volume = 0.0, viscosity = 0.0;
        if (edge.masterCell != -1 && edge.slaveCell != -1)
        {
            const auto& high = mesh.cells(edge.masterCell);
            const auto& low = mesh.cells(edge.slaveCell);
            volume = 2.0 * Pi * edge.midpointCoord().x * high.area();
            viscosity = 0.5 * (mu.value[high.index()] + mu.value[low.index()]);
            if (high.level() != low.level())
            {
                // Use the same virtual-child material and fine dual area as
                // the existing native RZ momentum row, not the coarse area.
                const auto material = HydroNativePressureStencil::faceData(mesh, edge, false, mu, mu,
                    [&](const Cell& cell, int child) -> const supCoeff& {
                        return supSonCells_[cell.supCellIndex * 4 + child].center;
                    });
                volume = 2.0 * Pi * edge.midpointCoord().x * material.dualArea;
                viscosity = material.density;
            }
        }
        else if (!boundaryFaces_.empty() && boundaryFaces_[i].outlet)
        {
            const auto& boundary = boundaryFaces_[i];
            const long double r = edge.midpointCoord().x;
            const long double delta = -boundary.sign * 0.5L * boundary.normalSpacing;
            const long double dr = -boundary.sign * boundary.normalSpacing;
            const long double logarithm = std::log1p(delta / r);
            // Exact integral of the native linear u(r) over the radial half
            // dual: (1/2)*integral 2mu*u/r^2 dV = 2pi*mu*dz*integral u/r dr.
            // In particular u=a*r cancels its normal-stress face divergence.
            const long double innerWeight = -boundary.sign * (delta - r * logarithm) / dr;
            const long double boundaryWeight = -boundary.sign * logarithm - innerWeight;
            const long double multiplier = 2.0L * Pi * edge.length() * mu.value[boundary.owner];
            const double diagonal = static_cast<double>(multiplier * boundaryWeight);
            const double neighbor = static_cast<double>(multiplier * innerWeight);
            if (!(r > 0.0L) || !(r + delta > 0.0L) || !std::isfinite(diagonal)
                || !std::isfinite(neighbor) || boundaryWeight < 0.0L || innerWeight < 0.0L)
                throw std::runtime_error("invalid radial outlet integrated hoop stress");
            assembly.entries.emplace_back(i, i, diagonal);
            assembly.entries.emplace_back(i, boundary.inner, neighbor);
            assembly.rhs[i] -= diagonal * ux.value[i] + neighbor * ux.value[boundary.inner];
            continue;
        }
        else continue; // Includes exactly prescribed u_r=0 at r=0.
        const double r = edge.midpointCoord().x;
        const double coefficient = viscosity * volume / (r * r);
        if (!(r > 0.0) || !std::isfinite(coefficient) || coefficient < 0.0)
            throw std::runtime_error("invalid axisymmetric hoop-stress coefficient");
        assembly.entries.emplace_back(i, i, coefficient);
        assembly.rhs[i] -= coefficient * ux.value[i];
    }
}

template<typename Mesh>
HydroBoundaryProjection Hydrodynamics::assembleBoundaryProjection(Mesh& mesh,
    double inverseMassScale, bool axisymmetric, const std::vector<BoundaryFaceData>& boundaries,
    const VEdgeField<double>& vx, const HEdgeField<double>& vy) const
{
    const int nx = mesh.vEdgeIndexNum(), count = nx + mesh.hEdgeIndexNum();
    std::vector<HydroBoundaryProjection::Face> faces(count);
    for (int row = 0; row < count; ++row)
    {
        const bool horizontal = row >= nx;
        const int k = horizontal ? row - nx : row;
        const auto& edge = horizontal ? mesh.hEdges(mesh.hEdgeIndexToLabel[k]) : mesh.vEdges(mesh.vEdgeIndexToLabel[k]);
        const auto& boundary = boundaries[row];
        auto& face = faces[row];
        face.area = edge.length() * (axisymmetric ? 2.0 * Pi * edge.midpointCoord().x : 1.0);
        if (boundary.owner >= 0)
        {
            face.high = boundary.owner;
            face.sign = boundary.sign;
            face.fixed = !boundary.outlet;
            if (face.fixed) face.fixedVelocity = boundary.sign * boundary.value;
            else face.inverseMass = inverseMassScale / (rho.value[boundary.owner] * boundary.volume);
        }
        else
        {
            const auto& high = mesh.cells(edge.masterCell);
            const auto& low = mesh.cells(edge.slaveCell);
            const double volume = high.area() * (axisymmetric ? 2.0 * Pi * edge.midpointCoord().x : 1.0);
            face.inverseMass = inverseMassScale / (0.5 * (rho.value[high.index()] + rho.value[low.index()]) * volume);
            face.high = high.index(); face.low = low.index();
        }
    }
    return HydroBoundaryProjection(mesh.cellIndexNum(), std::move(faces),
        [&](std::size_t row) { return row < nx ? vx.value[row] : vy.value[row - nx]; });
}

template<typename Mesh>
void Hydrodynamics::projectBoundaryVelocity(Mesh& mesh, double dt, bool axisymmetric)
{
    const int nx = mesh.vEdgeIndexNum(), cells = mesh.cellIndexNum();
    const auto projection = assembleBoundaryProjection(mesh, dt, axisymmetric, boundaryFaces_, uxMid, uyMid);
    const auto& system = projection.equations();
    const auto correction = checkedLinearSolve(system.matrix, system.rhs, "pressure", activeReport_.pressureRelativeResidual);
    pressurePotential_.initialize(mesh, 0.0);
    for (int i = 0; i < cells; ++i) pressurePotential_.value[i] = correction[i];
    projection.apply(correction,
        [&](std::size_t row) { return row < nx ? uxMid.value[row] : uyMid.value[row - nx]; },
        [&](std::size_t row, double value) { if (row < nx) ux.value[row] = value; else uy.value[row - nx] = value; },
        false, "nonfinite boundary-compatible projected velocity");
}

template<typename Mesh>
void Hydrodynamics::checkOutletVelocity(Mesh& mesh, double dt, bool withPredictor) const
{
    if (boundaryFaces_.empty()) return;
    double scale = 0.0;
    for (int i = 0; i < ux.size(); ++i)
        scale = std::max(scale, std::max(std::abs(ux.value[i]), std::abs(uxPre.value[i])));
    for (int i = 0; i < uy.size(); ++i)
        scale = std::max(scale, std::max(std::abs(uy.value[i]), std::abs(uyPre.value[i])));
    const double relativeTolerance = (activeOptions_.completeConservation->outletAdmissibility
        + 128.0 * std::numeric_limits<double>::epsilon()) * scale;
    // A static pressure/force balance leaves velocities at roundoff. Scaling
    // roundoff by that already-cancelled result would reject every tiny signed
    // outlet fluctuation. Use uncancelled pressure impulses (velocity units),
    // and only at the post-solve check the actual current predictor. Pressure
    // anywhere in the domain matters because projection is globally coupled.
    // This is a machine-epsilon term, NOT linearTolerance times a new scale;
    // it introduces no unit-dependent floor and never clips a velocity.
    long double cancellationScale = 0.0L;
    for (int i = 0; i < mesh.cellIndexNum(); ++i)
    {
        const auto& cell = mesh.cells(mesh.cellIndexToLabel[i]);
        const long double impulse = static_cast<long double>(dt) * 2.0L * std::abs(p.value[i])
            / rho.value[i] / std::min(cell.size().x, cell.size().y);
        cancellationScale = std::max(cancellationScale, impulse);
    }
    if (withPredictor)
    {
        for (int i = 0; i < uxMid.size(); ++i)
            cancellationScale = std::max(cancellationScale, static_cast<long double>(std::abs(uxMid.value[i])));
        for (int i = 0; i < uyMid.size(); ++i)
            cancellationScale = std::max(cancellationScale, static_cast<long double>(std::abs(uyMid.value[i])));
    }
    for (std::size_t row = 0; row < boundaryFaces_.size(); ++row)
        if (boundaryFaces_[row].outlet)
        {
            const auto& face = boundaryFaces_[row];
            const long double boundaryImpulse = static_cast<long double>(dt)*face.area
                * (std::abs(static_cast<long double>(p.value[face.owner]))+std::abs(static_cast<long double>(face.value)))
                / rho.value[face.owner] / face.volume;
            const long double tolerance = relativeTolerance + 128.0L*std::numeric_limits<double>::epsilon()
                * std::max(cancellationScale, boundaryImpulse);
            if (!std::isfinite(tolerance) || !std::isfinite(cancellationScale) || !std::isfinite(boundaryImpulse))
                throw std::runtime_error("pressure outlet velocity resolution scale is nonfinite");
            const double component = row < static_cast<std::size_t>(mesh.vEdgeIndexNum())
                ? ux.value[row] : uy.value[row - mesh.vEdgeIndexNum()];
            const double previous = row < static_cast<std::size_t>(mesh.vEdgeIndexNum())
                ? uxPre.value[row] : uyPre.value[row - mesh.vEdgeIndexNum()];
            if (boundaryFaces_[row].sign * component < -tolerance
                || boundaryFaces_[row].sign * previous < -tolerance)
            {
                std::ostringstream reason;
                reason << std::setprecision(17)
                    << "pressure outlet backflow requires explicit incoming velocity/material data; outward="
                    << face.sign*component << " previous=" << face.sign*previous
                    << " tolerance=" << tolerance;
                throw std::runtime_error(reason.str());
            }
        }
}
