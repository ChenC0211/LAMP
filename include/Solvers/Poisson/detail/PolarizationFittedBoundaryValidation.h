#pragma once

// Private, topology-preserving validation for the face-fitted active-leaf port.
// It deliberately knows neither material phase nor solid/interface geometry.
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/CartesianMeasure.h"
#include "Solvers/Poisson/PolarizationPoissonFVSolverCellCentered.h"

namespace PolarizationFittedBoundaryValidation
{
inline void require(bool condition, const char* message)
{
    if (!condition) throw std::runtime_error(std::string("fitted boundary: ") + message);
}

inline bool near(double a, double b, double scale)
{
    return std::isfinite(a) && std::isfinite(b) && std::isfinite(scale) &&
        std::abs(a - b) <= 256.0 * std::numeric_limits<double>::epsilon() * scale;
}

inline void validate(const QuadtreeMesh& mesh, const FaceSegmentGraph& graph,
    const CellField<double>& rhs, const CellField<double>& phi,
    const std::vector<double>& kappa,
    const std::vector<PolarizationPoissonFVSolverCellCentered::FittedBoundary>& boundaries,
    const PolarizationPoissonFVSolverCellCentered::Options& opt)
{
    using Type = PolarizationPoissonFVSolverCellCentered::FittedBoundaryType;
    require(opt.requireSuperLU && opt.useSuperLU && !opt.useCG && !opt.useIC,
            "requires requireSuperLU=true, useSuperLU=true, useCG=false, useIC=false");
#ifndef EIGEN_SUPERLU_SUPPORT
    require(false, "required SuperLU support is unavailable; no fallback is allowed");
#endif
    require(std::isfinite(opt.relTol) && opt.relTol > 0.0 && opt.relTol < 1.0,
            "relTol must be finite and in (0,1)");
    const int count = mesh.cellIndexNum();
    require(count > 0 && mesh.cellIndexToLabel, "missing mesh cell index mapping");
    require(std::isfinite(mesh.minCoord.x) && std::isfinite(mesh.minCoord.y) &&
            std::isfinite(mesh.maxCoord.x) && std::isfinite(mesh.maxCoord.y) &&
            mesh.maxCoord.x > mesh.minCoord.x && mesh.maxCoord.y > mesh.minCoord.y,
            "invalid domain bounds");
    if (opt.axisymmetric)
        require(mesh.minCoord.x >= 0.0, "axisymmetric domain requires r>=0");
    require(rhs.size() == count && rhs.value && phi.size() == count && phi.value,
            "rhs and phi must already be mesh-initialized");
    require(&rhs != &phi && rhs.value.get() != phi.value.get(), "rhs/phi aliasing is unsupported");
    const auto& segments = graph.segments();
    const auto& offsets = graph.cellOffsets();
    require(offsets.size() == static_cast<std::size_t>(count) + 1 &&
            offsets.front() == 0 && offsets.back() == segments.size() && !segments.empty(),
            "invalid or empty active graph offsets");
    require(kappa.size() == segments.size() && boundaries.size() == segments.size(),
            "coefficient/boundary arrays must cover every oriented segment");

    std::vector<unsigned char> active(static_cast<std::size_t>(count), 0);
    for (int ci = 0; ci < count; ++ci)
    {
        require(offsets[ci] <= offsets[ci + 1] && offsets[ci + 1] <= segments.size(),
                "invalid graph cell range");
        const int label = mesh.cellIndexToLabel[ci];
        require(label >= 0 && label < mesh.cellLabelCapacity(), "cell label out of range");
        const auto& cell = mesh.cells(label);
        require(cell.index() == ci && cell.label() == label &&
                rhs.label(ci) == label && phi.label(ci) == label, "stale mesh/field cell identity");
        require(std::isfinite(rhs.value[ci]) && std::isfinite(phi.value[ci]),
                "non-finite rhs or initial phi");
        if (offsets[ci] == offsets[ci + 1]) continue;
        active[ci] = 1;
        for (int son : cell.son) require(son == -1, "active unknown must be an actual leaf cell");
        const auto size = cell.size(), center = cell.centerCoord();
        require(std::isfinite(size.x) && size.x > 0.0 && std::isfinite(size.y) && size.y > 0.0 &&
                std::isfinite(cell.area()) && cell.area() > 0.0 &&
                std::isfinite(center.x) && std::isfinite(center.y), "invalid active cell metric");
        require(center.x > mesh.minCoord.x && center.x < mesh.maxCoord.x &&
                center.y > mesh.minCoord.y && center.y < mesh.maxCoord.y, "cell center outside domain");
        CartesianMeasure::cell(cell.area(), center.x, opt.axisymmetric);
    }

    // Derive only the active predicate from offsets. Rebuilding retains all
    // actual labels, levels, parent/child identities and split-edge incidence.
    // Exact comparison prevents a caller from deleting an individual face.
    FaceSegmentGraph actual;
    actual.build(mesh, [&active](const QuadtreeMesh&, const Cell& cell) {
        return active[static_cast<std::size_t>(cell.index())] != 0;
    });
    require(actual.cellOffsets() == offsets && actual.segments().size() == segments.size(),
            "graph does not cover all faces of the selected actual leaf cells");
    // Native labels already identify the same (orientation, edge) key used by
    // the former ordered map. Dense lookup avoids one tree allocation per
    // physical face on EVERY electric iteration; all pairing checks remain.
    const auto unseen = std::numeric_limits<std::size_t>::max();
    std::vector<std::size_t> firstVertical(mesh.vEdgeLabelCapacity(), unseen);
    std::vector<std::size_t> firstHorizontal(mesh.hEdgeLabelCapacity(), unseen);
    for (int ci = 0; ci < count; ++ci)
    {
        if (!active[ci]) continue;
        const auto& cell = mesh.cells(mesh.cellIndexToLabel[ci]);
        const auto center = cell.centerCoord(), size = cell.size();
        std::array<double, 4> covered{};
        for (std::size_t sid = offsets[ci]; sid < offsets[ci + 1]; ++sid)
        {
            const auto& face = segments[sid];
            const auto& expected = actual.segments()[sid];
            require(face.K == expected.K && face.L == expected.L && face.side == expected.side &&
                    face.boundary == expected.boundary && face.isHorizontal == expected.isHorizontal &&
                    face.edgeLabel == expected.edgeLabel && face.Af == expected.Af && face.df == expected.df &&
                    face.n.x == expected.n.x && face.n.y == expected.n.y &&
                    face.faceCenter.x == expected.faceCenter.x && face.faceCenter.y == expected.faceCenter.y,
                    "stale or altered graph geometry/identity");
            const int side = static_cast<int>(face.side);
            require(side >= 0 && side < 4 && face.K == mesh.cellIndexToLabel[ci], "invalid owner side");
            require(std::isfinite(face.Af) && face.Af > 0.0 && std::isfinite(face.df) && face.df > 0.0 &&
                    std::isfinite(face.faceCenter.x) && std::isfinite(face.faceCenter.y), "invalid segment metric");
            covered[side] += face.Af;
            const double halfDistance = (face.faceCenter - center).dotProduct(face.n);
            const double normalSize = face.isHorizontal ? size.y : size.x;
            const double normalCenter = face.isHorizontal ? center.y : center.x;
            require(halfDistance > 0.0 && near(halfDistance, 0.5 * normalSize,
                    std::max(std::abs(normalCenter), normalSize)), "face is not on owner cell boundary");
            int neighbor = -1;
            if (face.L != -1)
            {
                require(face.L >= 0 && face.L < mesh.cellLabelCapacity(), "neighbor label out of range");
                neighbor = mesh.cells(face.L).index();
                require(neighbor >= 0 && neighbor < count && mesh.cellIndexToLabel[neighbor] == face.L,
                        "neighbor is not a current indexed cell");
            }
            else
            {
                const double wall = side == 0 ? mesh.minCoord.x : side == 1 ? mesh.minCoord.y :
                                    side == 2 ? mesh.maxCoord.x : mesh.maxCoord.y;
                const double coordinate = face.isHorizontal ? face.faceCenter.y : face.faceCenter.x;
                require(near(coordinate, wall, std::max(std::abs(wall), normalSize)),
                        "missing neighbor away from a physical domain wall");
            }
            const auto& boundary = boundaries[sid];
            require(std::isfinite(boundary.value), "non-finite boundary value");
            const bool interior = neighbor >= 0 && active[neighbor];
            if (interior)
                require(boundary.type == Type::Interior && boundary.value == 0.0,
                        "active-active faces must be Interior with zero unused value");
            else
                require(boundary.type == Type::Dirichlet ||
                        (boundary.type == Type::HomogeneousNeumann && boundary.value == 0.0),
                        "domain/inactive faces require static Dirichlet or exactly zero Neumann");
            require(std::isfinite(kappa[sid]) && kappa[sid] > 0.0, "kappa must be finite and positive");
            const double area = CartesianMeasure::face(face.Af, face.faceCenter.x, opt.axisymmetric);
            if (area == 0.0)
                require(opt.axisymmetric && face.L == -1 && side == 0 && mesh.minCoord.x == 0.0 &&
                        face.faceCenter.x == 0.0 && boundary.type == Type::HomogeneousNeumann,
                        "zero-area axis faces require exactly zero Neumann");
            else if (boundary.type != Type::HomogeneousNeumann)
            {
                const double conductance = kappa[sid] * area / (interior ? face.df : halfDistance);
                require(std::isfinite(conductance) && conductance > 0.0, "unrepresentable face conductance");
            }
            auto& firstByEdge = face.isHorizontal ? firstHorizontal : firstVertical;
            require(face.edgeLabel >= 0 && static_cast<std::size_t>(face.edgeLabel) < firstByEdge.size(),
                    "physical face label is outside native edge storage");
            auto& first = firstByEdge[static_cast<std::size_t>(face.edgeLabel)];
            if (first == unseen) first = sid;
            else
            {
                const auto otherSid = first;
                require(otherSid < segments.size(), "more than two owners of a physical face");
                const auto& other = segments[otherSid];
                require(interior && other.K == face.L && other.L == face.K &&
                        other.n.x == -face.n.x && other.n.y == -face.n.y &&
                        other.Af == face.Af && other.df == face.df && kappa[otherSid] == kappa[sid] &&
                        boundaries[otherSid].type == Type::Interior, "inconsistent reverse active face");
                // A third incidence must fail rather than reuse a completed pair.
                first = segments.size();
            }
        }
        for (int side = 0; side < 4; ++side)
        {
            const double length = side == 0 || side == 2 ? size.y : size.x;
            require(near(covered[side], length, length), "incomplete/duplicate active cell face coverage");
        }
    }
    for (const auto* firstByEdge : {&firstVertical, &firstHorizontal})
        for (std::size_t first : *firstByEdge)
        {
            if (first == unseen || first == segments.size()) continue;
            require(boundaries[first].type != Type::Interior, "missing reverse active face");
        }
    // Component anchoring is checked by the shared assembly before factorization.
}
} // namespace PolarizationFittedBoundaryValidation
