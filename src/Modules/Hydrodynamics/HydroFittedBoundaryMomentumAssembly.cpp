#include "Modules/Hydrodynamics/Detail/HydroFittedBoundaryMomentum.h"
#include "Modules/Hydrodynamics/HydroFittedMomentum.h"
#include "Modules/Hydrodynamics/Detail/HydroMomentumAssembly.h"
#include "Modules/Hydrodynamics/HydroFittedDomain.h"
#include "Core/Field/CellField.h"
#include "Core/Field/VertField.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/CartesianMeasure.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace HydroFittedDetail
{
void BoundaryMomentum::append(const QuadtreeMesh& mesh, const HydroFittedDomain& domain, double dt,
    const CellField<double>& rho, const CellField<double>& mu,
    const VertField<double>& vertexRho, const VertField<double>& vertexMu,
    const Eigen::VectorXd& current, const Eigen::VectorXd& previous,
    const CellField<double>& bodyX, const CellField<double>& bodyY,
    HydroFittedMomentumSystem& system, HydroMomentumAssembly& assembly) const
{
    using Type = HydroFittedDomain::BoundaryType;
    const auto require = [](bool good, const char* reason) {
        if (!good) throw std::invalid_argument(reason);
    };
    const auto value = [&](int native, const Eigen::VectorXd& history) {
        require(native >= 0 && native < history.size(), "missing boundary momentum native sample");
        const auto& face = domain.faces()[native];
        require((face.negativeCell >= 0 || face.positiveCell >= 0) && std::isfinite(history[native]),
            "boundary momentum attempted a solid-only/nonfinite velocity read");
        return history[native];
    };
    const auto evaluate = [&](const MomentumTrace& trace, const Eigen::VectorXd& history) {
        double result = 0.0;
        for (const auto& entry : trace) result += entry.second * value(entry.first, history);
        require(std::isfinite(result), "nonfinite boundary momentum trace");
        return result;
    };
    const auto coefficient = [&](int row, int column, double number) {
        const auto& target = domain.faces()[row]; const auto& source = domain.faces()[column];
        require(target.velocityRow >= 0 && (source.negativeCell >= 0 || source.positiveCell >= 0)
            && std::isfinite(number), "invalid shared momentum coefficient");
        if (source.velocityRow >= 0) assembly.entries.emplace_back(target.velocityRow, source.velocityRow, number);
        else require((source.boundary.type == Type::NoSlipWall || source.boundary.type == Type::Symmetry)
            && current[column] == 0.0 && previous[column] == 0.0, "invalid eliminated shared momentum constraint");
    };
    std::vector<double> scale(selected_.size(), 0.0);
    for (int native : domain.nativeByVelocityRow()) if (selected_[native])
    {
        const auto& face = domain.faces()[native];
        const auto& volume = duals_.volumes()[native];
        const auto& body = face.horizontal ? bodyY : bodyX;
        double mass = 0.0, force = 0.0;
        for (int side = 0; side < 2; ++side)
        {
            const int owner = side == 0 ? face.negativeCell : face.positiveCell;
            if (owner < 0) continue;
            const double part = ownerVolumes_[native][side];
            mass += part * rho.value[owner]; force += part * body.value[owner];
        }
        // Whole-equation normalization: C and its projection inverse mass
        // retain the existing pressure metric. Physical force/flux balance is
        // obtained by dividing this row by scale, not by mixing dual volumes.
        scale[native] = face.dualVolume / volume.measure;
        const double density = mass / volume.measure;
        const double diagonal = face.dualVolume * density / dt;
        require(std::isfinite(scale[native]) && scale[native] > 0.0 && std::isfinite(density) && density > 0.0,
            "invalid normalized physical momentum mass");
        system.densityByVelocityRow[face.velocityRow] = density;
        coefficient(native, native, diagonal);
        assembly.rhs[face.velocityRow] = diagonal * current[native] + scale[native] * force;
        ++system.boundaryPieceRows;
        if (face.horizontal) continue;

        // Original RZ -2*mu*u_r/r^2 term, integrated over each real owner half
        // using a linear native trace. This exactly cancels constant radial
        // normal stress for u_r=a*r, including an asymmetric coarse/fine dual.
        long double constant = 0.0L, slope = 0.0L;
        for (int side = 0; side < 2; ++side)
        {
            const int owner = side == 0 ? face.negativeCell : face.positiveCell;
            if (owner < 0) continue;
            const long double r = face.center.x;
            const long double c = mesh.cells(mesh.cellIndexToLabel[owner]).centerCoord().x;
            const long double low = std::min(c, r), high = std::max(c, r);
            require(low > 0.0L && high > low, "free radial momentum volume intersects a singular axis");
            const long double logarithm = std::log1p((high - low) / low);
            const long double multiplier = CartesianMeasure::twoPi * mu.value[owner]
                * (static_cast<long double>(volume.upper.y) - volume.lower.y);
            constant += multiplier * logarithm;
            slope += multiplier * ((high - low) - r * logarithm);
        }
        const double diagonalHoop = static_cast<double>(scale[native] * constant);
        const double slopeHoop = static_cast<double>(scale[native] * slope);
        coefficient(native, native, diagonalHoop);
        for (const auto& term : radialSlope_[native]) coefficient(native, term.first, slopeHoop * term.second);
        assembly.rhs[face.velocityRow] -= diagonalHoop * current[native]
            + slopeHoop * evaluate(radialSlope_[native], current);
    }
    const auto material = [&](const Connection& connection, bool viscosity) {
        const auto& cells = viscosity ? mu : rho;
        const auto& vertices = viscosity ? vertexMu : vertexRho;
        double result = 0.0;
        if (connection.materialCell >= 0) result = cells.value[connection.materialCell];
        else if (connection.materialVertex >= 0) result = vertices.value[connection.materialVertex];
        else
        {
            for (int id : connection.parts)
            {
                const auto& part = duals_.segments()[id];
                double sample = 0.0;
                if (part.ownerCell >= 0) sample = cells.value[part.ownerCell];
                else
                {
                    const auto& face = domain.faces()[part.primalNative];
                    if (face.negativeCell < 0 || face.positiveCell < 0)
                        sample = cells.value[face.negativeCell >= 0 ? face.negativeCell : face.positiveCell];
                    else
                    {
                        const auto& edge = face.horizontal ? mesh.hEdges(face.edgeLabel) : mesh.vEdges(face.edgeLabel);
                        const auto& start = mesh.verts(edge.edgeToVert[0]); const auto& end = mesh.verts(edge.edgeToVert[1]);
                        double center = 0.5 * (part.begin.y + part.end.y);
                        if (face.horizontal)
                            center = (2.0 / 3.0) * (part.begin.x * part.begin.x + part.begin.x * part.end.x + part.end.x * part.end.x)
                                / (part.begin.x + part.end.x);
                        const double lo = face.horizontal ? start.coord().x : start.coord().y;
                        const double hi = face.horizontal ? end.coord().x : end.coord().y;
                        const double fraction = (center - lo) / (hi - lo);
                        require(std::isfinite(fraction) && fraction >= 0.0 && fraction <= 1.0,
                            "shared material trace leaves its actual native face");
                        sample = (1.0 - fraction) * vertices.value[start.index()] + fraction * vertices.value[end.index()];
                    }
                }
                result += part.measure * sample;
            }
            result /= connection.area;
        }
        require(std::isfinite(result) && (viscosity ? result >= 0.0 : result > 0.0),
            "invalid active shared momentum material");
        return result;
    };
    for (const auto& connection : connections_)
    {
        if (connection.area == 0.0) continue; // Only the true RZ symmetry axis.
        double viscosity = material(connection, true);
        double flux[2]{};
        double traction = 0.0;
        const bool exterior = connection.kind == Kind::Exterior;
        if (exterior)
        {
            const auto& boundary = domain.faces()[duals_.segments()[connection.parts[0]].primalNative];
            const int native = connection.negative >= 0 ? connection.negative : connection.positive;
            const auto& face = domain.faces()[native];
            const bool open = boundary.boundary.type == Type::PressureReservoir || boundary.boundary.type == Type::PressureOutlet;
            if (open)
            {
                viscosity = 0.0; // Complete prescribed traction, not another viscous boundary force.
                if (face.horizontal == boundary.horizontal) traction = -boundary.boundary.pressure * connection.area;
                for (int history = 0; history < 2; ++history)
                {
                    const auto& state = history == 0 ? current : previous;
                    const double normal = value(boundary.nativeRow, state);
                    const double outward = boundary.outwardSign * normal;
                    const int owner = boundary.negativeCell >= 0 ? boundary.negativeCell : boundary.positiveCell;
                    const double density = boundary.boundary.type == Type::PressureReservoir && outward < 0.0
                        ? boundary.boundary.incomingDensity : rho.value[owner];
                    // Preserve the original outlet rule even for tolerated
                    // roundoff-sized inward normal velocity. Only reservoirs
                    // continue the normal component from the interior.
                    const bool continueIncoming = outward < 0.0
                        && (face.horizontal != boundary.horizontal || boundary.boundary.type == Type::PressureReservoir);
                    const double donor = continueIncoming ? evaluate(connection.incoming, state) : value(native, state);
                    flux[history] = density * normal * connection.area * donor;
                }
            }
            // Wall/symmetry mass flux stays exactly zero, including backflow histories.
        }
        else
        {
            const double density = material(connection, false);
            for (int history = 0; history < 2; ++history)
            {
                const auto& state = history == 0 ? current : previous;
                const double q = density * evaluate(connection.advector, state) * connection.area;
                flux[history] = q * value(q > 0.0 ? connection.negative : connection.positive, state);
            }
        }
        const auto& explicitStrain = connection.currentStrain.empty() ? connection.strain : connection.currentStrain;
        const double currentStress = viscosity * connection.area * evaluate(explicitStrain, current);
        for (const auto side : {std::make_pair(connection.negative, 1.0), std::make_pair(connection.positive, -1.0)})
        {
            const int native = side.first;
            if (native < 0 || !selected_[native]) continue;
            const double orientation = side.second * scale[native];
            const int compact = domain.faces()[native].velocityRow;
            for (const auto& term : connection.strain)
                coefficient(native, term.first, -0.5 * orientation * viscosity * connection.area * term.second);
            assembly.rhs[compact] += orientation * (0.5 * currentStress + traction - 1.5 * flux[0] + 0.5 * flux[1]);
        }
    }
}
}
