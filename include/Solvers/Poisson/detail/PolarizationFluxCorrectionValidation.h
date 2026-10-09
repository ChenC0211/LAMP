#pragma once

// Private validation for the generic integrated-flux port. No material/interface
// model is assumed here; phase and cut-edge diagnostics are deliberately ignored.
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/CartesianMeasure.h"
#include "Solvers/Poisson/PolarizationPoissonFVSolverCellCentered.h"

namespace PolarizationFluxCorrectionValidation
{
inline void require(bool condition, const std::string& message)
{
    if (!condition) throw std::runtime_error("fluxCorrections: " + message);
}

inline bool near(double a, double b, double scale)
{
    return std::isfinite(a) && std::isfinite(b) && std::isfinite(scale) &&
        std::abs(a - b) <= 256.0 * std::numeric_limits<double>::epsilon() * scale;
}

inline void validate(
    const QuadtreeMesh& mesh, const FaceSegmentGraph& graph, const PicDomainBC& bc,
    const CellField<double>& rhs, const CellField<double>& phi,
    const std::vector<double>* kappaBySegment, const std::vector<double>* kappaByCell,
    const PolarizationPoissonFVSolverCellCentered::Options& opt,
    const Solvers::Poisson::CutCellMetrics* cutMetrics,
    const PicEmbeddedBC* embeddedBC, const PicEmbeddedBCMap* embeddedBCById,
    const std::vector<PolarizationPoissonFVSolverCellCentered::SegmentFluxCorrection>& corrections)
{
    require(opt.requireSuperLU && opt.useSuperLU && !opt.useCG && !opt.useIC,
            "requires requireSuperLU=true, useSuperLU=true, useCG=false, useIC=false");
#ifndef EIGEN_SUPERLU_SUPPORT
    require(false, "required SuperLU support is unavailable; no fallback is allowed");
#endif
    require(std::isfinite(opt.relTol) && opt.relTol > 0.0 && opt.relTol < 1.0,
            "relTol must be finite and in (0,1)");
    require(!cutMetrics && !embeddedBC && !embeddedBCById,
            "cutMetrics and embedded BC pointers are unsupported");
    bool hasDirichlet = false;
    for (const auto* side : {&bc.left, &bc.right, &bc.bottom, &bc.top})
    {
        if (side->type == PicBCType::Neumann)
            require(std::isfinite(side->value) && side->value == 0.0 &&
                    !side->timeDirichlet.enabled,
                    "Neumann domain sides must be static and exactly homogeneous (value=0)");
        else
        {
            // Keep the existing Dirichlet predicate and diagnostic order.
            require(side->type == PicBCType::Dirichlet && std::isfinite(side->value) &&
                    !side->timeDirichlet.enabled, "all domain sides must be static finite Dirichlet BCs");
            hasDirichlet = true;
        }
    }
    require(hasDirichlet, "at least one static Dirichlet domain side is required; all-Neumann gauge pinning is forbidden");

    const int count = mesh.cellIndexNum();
    require(mesh.Nx >= 1 && mesh.Ny >= 1 &&
            static_cast<long long>(mesh.Nx) * mesh.Ny == count &&
            (static_cast<long long>(mesh.Nx) + 1) * (static_cast<long long>(mesh.Ny) + 1) == mesh.vertIndexNum() &&
            static_cast<long long>(mesh.Nx) * (static_cast<long long>(mesh.Ny) + 1) == mesh.hEdgeIndexNum() &&
            (static_cast<long long>(mesh.Nx) + 1) * mesh.Ny == mesh.vEdgeIndexNum(),
            "requires a complete uniform Cartesian mesh (no AMR or omitted cells)");
    require(mesh.cellIndexToLabel && mesh.vertIndexToLabel && mesh.hEdgeIndexToLabel && mesh.vEdgeIndexToLabel,
            "missing mesh index mappings");
    const double hx = mesh.cellSize.x, hy = mesh.cellSize.y;
    require(std::isfinite(hx) && hx > 0.0 && std::isfinite(hy) && hy > 0.0 &&
            std::isfinite(hx * hy) && hx * hy > 0.0 &&
            std::isfinite(mesh.minCoord.x) && std::isfinite(mesh.minCoord.y) &&
            std::isfinite(mesh.maxCoord.x) && std::isfinite(mesh.maxCoord.y), "invalid mesh geometry");
    if (opt.axisymmetric)
        require(mesh.minCoord.x >= 0.0, "axisymmetric radial domain must start at r>=0");
    const double sx = std::max({std::abs(mesh.minCoord.x), std::abs(mesh.maxCoord.x), hx});
    const double sy = std::max({std::abs(mesh.minCoord.y), std::abs(mesh.maxCoord.y), hy});
    require(hx > 1024.0 * std::numeric_limits<double>::epsilon() * sx &&
            hy > 1024.0 * std::numeric_limits<double>::epsilon() * sy &&
            near(mesh.maxCoord.x - mesh.minCoord.x, mesh.Nx * hx, sx) &&
            near(mesh.maxCoord.y - mesh.minCoord.y, mesh.Ny * hy, sy),
            "mesh extents/spacing are inconsistent or numerically unresolved");
    require(rhs.size() == count && rhs.value, "rhs must be initialized with one finite value per cell index");
    require(&rhs != &phi && (!phi.value || phi.value.get() != rhs.value.get()), "rhs/phi aliasing is unsupported");
    for (int ci = 0; ci < count; ++ci)
        require(std::isfinite(rhs.value[ci]), "rhs contains a non-finite value");

    const auto& segments = graph.segments();
    const auto& offsets = graph.cellOffsets();
    require(segments.size() == 4u * static_cast<std::size_t>(count) &&
            offsets.size() == static_cast<std::size_t>(count) + 1u &&
            offsets.front() == 0u && offsets.back() == segments.size(),
            "graph must cover exactly four unsplit faces of every cell");
    require(corrections.size() == segments.size(), "correction count must equal oriented segment count");
    require(kappaBySegment && kappaBySegment->size() == segments.size(),
            "kappaBySegment must explicitly contain one value per oriented segment");
    for (double value : *kappaBySegment)
        require(std::isfinite(value) && value > 0.0, "segment kappa must be finite and positive");
    if (kappaByCell)
    {
        require(kappaByCell->size() == static_cast<std::size_t>(count), "cell kappa size mismatch");
        for (double value : *kappaByCell)
            require(std::isfinite(value) && value > 0.0, "cell kappa must be finite and positive");
    }

    std::vector<int> logicalCell(static_cast<std::size_t>(count), -1);
    std::vector<std::array<int, 2>> logicalPosition(static_cast<std::size_t>(count));
    for (int ci = 0; ci < count; ++ci)
    {
        const int label = mesh.cellIndexToLabel[ci];
        require(label >= 0 && label < mesh.cellLabelCapacity(), "cell label out of range");
        const auto& cell = mesh.cells(label);
        require(cell.index() == ci && cell.label() == label && cell.parent == -1 && cell.level() == 0,
                "invalid cell mapping or nonuniform/refined cell");
        for (int child : cell.son) require(child == -1, "AMR parent cells are unsupported");
        require(near(cell.size().x, hx, hx) && near(cell.size().y, hy, hy), "nonuniform cell sizes");
        const auto center = cell.centerCoord();
        if (opt.axisymmetric)
            CartesianMeasure::cell(cell.area(), center.x, true);
        const double ix = (center.x - mesh.minCoord.x) / hx - 0.5;
        const double iy = (center.y - mesh.minCoord.y) / hy - 0.5;
        require(std::isfinite(ix) && std::isfinite(iy) && ix > -0.5 && ix < mesh.Nx - 0.5 &&
                iy > -0.5 && iy < mesh.Ny - 0.5, "cell center is outside the Cartesian lattice");
        const int x = static_cast<int>(std::round(ix)), y = static_cast<int>(std::round(iy));
        require(near(center.x, mesh.minCoord.x + (x + 0.5) * hx, sx) &&
                near(center.y, mesh.minCoord.y + (y + 0.5) * hy, sy), "cell center is off the Cartesian lattice");
        const auto slot = static_cast<std::size_t>(y) * mesh.Nx + x;
        require(logicalCell[slot] == -1, "duplicate Cartesian cell location");
        logicalCell[slot] = ci;
        logicalPosition[static_cast<std::size_t>(ci)] = {x, y};
        require(offsets[ci] <= segments.size() && offsets[ci + 1] <= segments.size() &&
                offsets[ci + 1] == offsets[ci] + 4u, "invalid or incomplete graph cell offsets");
    }

    // FaceSide order is west, south, east, north, not DomainBC member order.
    const std::array<const PicBoundarySide*, 4> boundaryBySide{{&bc.left, &bc.bottom, &bc.right, &bc.top}};
    bool hasPositiveMeasureDirichlet = false;
    std::vector<std::array<std::size_t, 4>> segmentBySide(static_cast<std::size_t>(count));
    for (int ci = 0; ci < count; ++ci)
    {
        std::array<bool, 4> seen{};
        const int label = mesh.cellIndexToLabel[ci];
        const auto& cell = mesh.cells(label);
        const auto center = cell.centerCoord();
        const auto position = logicalPosition[static_cast<std::size_t>(ci)];
        for (std::size_t sid = offsets[ci]; sid < offsets[ci + 1]; ++sid)
        {
            const auto& face = segments[sid];
            const int side = static_cast<int>(face.side);
            require(side >= 0 && side < 4 && !seen[side] && face.K == label, "invalid/duplicate oriented owner face");
            seen[side] = true;
            segmentBySide[static_cast<std::size_t>(ci)][side] = sid;
            const auto normal = faceSideToNormal(face.side);
            const bool horizontal = side == 1 || side == 3;
            const int x = position[0] + (side == 0 ? -1 : (side == 2 ? 1 : 0));
            const int y = position[1] + (side == 1 ? -1 : (side == 3 ? 1 : 0));
            const bool boundary = x < 0 || x >= mesh.Nx || y < 0 || y >= mesh.Ny;
            const int neighborIndex = boundary ? -1 : logicalCell[static_cast<std::size_t>(y) * mesh.Nx + x];
            const int expectedNeighbor = boundary ? -1 : mesh.cellIndexToLabel[neighborIndex];
            require(face.L == expectedNeighbor && face.boundary == (boundary ? faceSideToBoundary(face.side) : Boundary::non),
                    "graph neighbor/boundary does not match the complete Cartesian mesh");
            const double distance = (horizontal ? hy : hx) * (boundary ? 0.5 : 1.0);
            require(face.isHorizontal == horizontal && face.n.x == normal.x && face.n.y == normal.y &&
                    near(face.Af, horizontal ? hx : hy, horizontal ? hx : hy) &&
                    near(face.df, distance, horizontal ? sy : sx) && face.df > 0.0 && face.Af > 0.0 &&
                    near(face.faceCenter.x, center.x + 0.5 * hx * normal.x, sx) &&
                    near(face.faceCenter.y, center.y + 0.5 * hy * normal.y, sy), "invalid Cartesian segment metric/normal");
            const int edgeSide = side >= 2 ? 1 : 0;
            const int* cellEdges = horizontal ? cell.cellToHEdge[edgeSide] : cell.cellToVEdge[edgeSide];
            const int capacity = horizontal ? mesh.hEdgeLabelCapacity() : mesh.vEdgeLabelCapacity();
            require(cellEdges[1] == -1 && face.edgeLabel == cellEdges[0] &&
                    face.edgeLabel >= 0 && face.edgeLabel < capacity, "invalid or split face edge mapping");
            const auto& edge = horizontal ? mesh.hEdges(face.edgeLabel) : mesh.vEdges(face.edgeLabel);
            const int edgeCount = horizontal ? mesh.hEdgeIndexNum() : mesh.vEdgeIndexNum();
            const auto* edgeLabels = horizontal ? mesh.hEdgeIndexToLabel : mesh.vEdgeIndexToLabel;
            require(edge.index() >= 0 && edge.index() < edgeCount && edgeLabels[edge.index()] == face.edgeLabel &&
                    edge.parent == -1 && edge.son[0] == -1 && edge.son[1] == -1 && edge.isHorizontal() == horizontal &&
                    edge.length() == face.Af && edge.midpointCoord().x == face.faceCenter.x &&
                    edge.midpointCoord().y == face.faceCenter.y, "stale graph or invalid edge geometry");
            if (opt.axisymmetric)
            {
                if (boundary && side == 0 && mesh.minCoord.x == 0.0)
                    require(face.faceCenter.x == 0.0, "axis face must retain its exact zero radius");
                const double measure = CartesianMeasure::face(face.Af, face.faceCenter.x, true);
                if (measure == 0.0)
                    require(boundary && side == 0 &&
                            boundaryBySide[static_cast<std::size_t>(side)]->type == PicBCType::Neumann,
                            "zero-area axis faces require homogeneous Neumann BCs");
                if (boundary && measure > 0.0 &&
                    boundaryBySide[static_cast<std::size_t>(side)]->type == PicBCType::Dirichlet)
                    hasPositiveMeasureDirichlet = true;
            }
            const auto& correction = corrections[sid];
            require(std::isfinite(correction.constant), "non-finite correction constant");
            int previous = -1;
            for (const auto& coefficient : correction.coefficients)
            {
                require(coefficient.cellIndex > previous && coefficient.cellIndex < count &&
                        std::isfinite(coefficient.value), "coefficients must have sorted unique active cell indices and finite values");
                previous = coefficient.cellIndex;
            }
            if (boundary && boundaryBySide[static_cast<std::size_t>(side)]->type == PicBCType::Neumann)
            {
                // A prescribed zero flux must not receive an added operator.
                // Empty or explicit signed-zero entries are legal, but their
                // index/finite checks above still apply. Do not erase a tiny
                // nonzero value or accept cancellation between coefficients.
                require(correction.constant == 0.0,
                        "homogeneous Neumann boundary correction constant must be exactly zero");
                for (const auto& coefficient : correction.coefficients)
                    require(coefficient.value == 0.0,
                            "homogeneous Neumann boundary correction coefficients must each be exactly zero");
            }
        }
    }
    if (opt.axisymmetric)
        require(hasPositiveMeasureDirichlet,
                "at least one positive-area Dirichlet face is required; the axis cannot anchor the solve");
    for (std::size_t sid = 0; sid < segments.size(); ++sid)
    {
        const auto& face = segments[sid];
        if (face.L == -1) continue;
        const int neighborIndex = mesh.cells(face.L).index();
        const auto reverseId = segmentBySide[static_cast<std::size_t>(neighborIndex)][(static_cast<int>(face.side) + 2) % 4];
        const auto& reverse = segments[reverseId];
        require(reverse.L == face.K && reverse.edgeLabel == face.edgeLabel && reverse.Af == face.Af &&
                reverse.df == face.df && reverse.faceCenter.x == face.faceCenter.x &&
                reverse.faceCenter.y == face.faceCenter.y && reverse.n.x == -face.n.x && reverse.n.y == -face.n.y &&
                (*kappaBySegment)[reverseId] == (*kappaBySegment)[sid],
                "interior oriented pair must have identical kappa/geometry and opposite normals");
        const auto& forwardCorrection = corrections[sid];
        const auto& reverseCorrection = corrections[reverseId];
        require(reverseCorrection.constant == -forwardCorrection.constant &&
                reverseCorrection.coefficients.size() == forwardCorrection.coefficients.size(),
                "interior correction pair must be exactly antisymmetric");
        for (std::size_t j = 0; j < forwardCorrection.coefficients.size(); ++j)
            require(reverseCorrection.coefficients[j].cellIndex == forwardCorrection.coefficients[j].cellIndex &&
                    reverseCorrection.coefficients[j].value == -forwardCorrection.coefficients[j].value,
                    "interior correction coefficients must be exactly antisymmetric");
    }
}
} // namespace PolarizationFluxCorrectionValidation
