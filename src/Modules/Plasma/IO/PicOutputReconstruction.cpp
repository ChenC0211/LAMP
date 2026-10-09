#include "Modules/Plasma/IO/PicOutputReconstruction.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/FaceSegmentGraph.h"
#include "Solvers/Poisson/CutCellMetrics.h"
#include "Core/Field/CellField.h"
#include "Core/Base/Interpolation.h"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <unordered_map>

namespace
{
std::int64_t periodicCoordKey(double coord, double eps)
{
    const double safeEps = (eps > 0.0 && std::isfinite(eps)) ? eps : 1e-12;
    return static_cast<std::int64_t>(std::llround(coord / safeEps));
}
}

namespace PicOutputReconstruction
{
std::vector<unsigned char> buildActiveCellMask(const QuadtreeMesh& mesh,
    const FaceSegmentGraph& faceGraph, const Solvers::Poisson::CutCellMetrics* cutMetricsPtr,
    const std::vector<double>* alphaByCellIndex)
{
    std::vector<unsigned char> outputCellValid(static_cast<std::size_t>(mesh.cellIndexNum()), 0u);
    const bool outputHasAlphaMask = alphaByCellIndex &&
        alphaByCellIndex->size() == static_cast<std::size_t>(mesh.cellIndexNum());
    for (int ci = 0; ci < mesh.cellIndexNum(); ++ci)
    {
        const int cellLabel = mesh.cellIndexToLabel[ci];
        const Cell& cell = mesh.cells(cellLabel);
        const bool isLeaf =
            (cell.son[0] == -1) && (cell.son[1] == -1) && (cell.son[2] == -1) && (cell.son[3] == -1);
        if (!isLeaf)
            continue;

        // For moving EB + AMR, visualization should follow the active fluid mask
        // from cut metrics. faceGraph can transiently miss segments on otherwise
        // valid cells, which creates artificial holes in Field plots.
        if (cutMetricsPtr != nullptr)
        {
            if (!cutMetricsPtr->isActiveCellIndex(ci))
                continue;
        }
        else
        {
            const auto range = faceGraph.segmentRangeForCellIndex(ci);
            if (!(range.first < range.second))
                continue;
        }

        if (outputHasAlphaMask)
        {
            double alpha = (*alphaByCellIndex)[static_cast<std::size_t>(ci)];
            if (!std::isfinite(alpha))
                alpha = 1.0;
            alpha = std::clamp(alpha, 0.0, 1.0);
            if (!(alpha > 0.0))
                continue;
        }

        outputCellValid[static_cast<std::size_t>(ci)] = 1u;
    }
    return outputCellValid;
}

std::vector<bool> markZeroPotentialCells(const QuadtreeMesh& mesh,
    const CellField<double>& phi, const std::vector<unsigned char>& outputCellValid)
{
    std::vector<bool> phiZeroArtifact(static_cast<std::size_t>(mesh.cellIndexNum()), false);
    {
        // Seed: valid cells with phi==0 that have at least one valid
        // neighbour with |phi|>0  OR  border an invalid cell.
        std::vector<int> seedQueue;
        for (int ci = 0; ci < mesh.cellIndexNum(); ++ci)
        {
            if (outputCellValid[static_cast<std::size_t>(ci)] == 0u)
                continue;
            if (phi.value[ci] != 0.0)
                continue;
            const int cLabel = mesh.cellIndexToLabel[ci];
            const Cell& cell = mesh.cells(cLabel);
            bool bordersNonZeroOrInvalid = false;
            for (int j = 0; j < 4 && !bordersNonZeroOrInvalid; ++j)
            {
                const Vert& v = mesh.verts(cell.cellToVert[j]);
                for (int k = 0; k < 4; ++k)
                {
                    const int nLabel = v.vertToCell[k];
                    if (nLabel < 0) { bordersNonZeroOrInvalid = true; break; }
                    const int nci = mesh.cells(nLabel).index();
                    if (nci < 0 || nci >= mesh.cellIndexNum() || nci == ci) continue;
                    if (outputCellValid[static_cast<std::size_t>(nci)] == 0u)
                    { bordersNonZeroOrInvalid = true; break; }
                    if (phi.value[nci] != 0.0)
                    { bordersNonZeroOrInvalid = true; break; }
                }
            }
            if (bordersNonZeroOrInvalid)
            {
                phiZeroArtifact[static_cast<std::size_t>(ci)] = true;
                seedQueue.push_back(ci);
            }
        }
        // Flood-fill: mark remaining phi=0 valid cells adjacent to
        // already-marked cells (handles interior of phi=0 clusters).
        for (std::size_t qi = 0; qi < seedQueue.size(); ++qi)
        {
            const int sci = seedQueue[qi];
            const int scLabel = mesh.cellIndexToLabel[sci];
            const Cell& scell = mesh.cells(scLabel);
            for (int j = 0; j < 4; ++j)
            {
                const Vert& v = mesh.verts(scell.cellToVert[j]);
                for (int k = 0; k < 4; ++k)
                {
                    const int nLabel = v.vertToCell[k];
                    if (nLabel < 0) continue;
                    const int nci = mesh.cells(nLabel).index();
                    if (nci < 0 || nci >= mesh.cellIndexNum()) continue;
                    if (phiZeroArtifact[static_cast<std::size_t>(nci)]) continue;
                    if (outputCellValid[static_cast<std::size_t>(nci)] == 0u) continue;
                    if (phi.value[nci] != 0.0) continue;
                    phiZeroArtifact[static_cast<std::size_t>(nci)] = true;
                    seedQueue.push_back(nci);
                }
            }
        }
    }
    return phiZeroArtifact;
}

void fillMarkedCells(const QuadtreeMesh& mesh, const std::vector<unsigned char>& outputCellValid,
    const std::vector<bool>& phiZeroArtifact, CellField<double>& cOut)
{
    const int N = mesh.cellIndexNum();
    std::vector<bool> filled(static_cast<std::size_t>(N), false);
    for (int pass = 0; pass < 20; ++pass)
    {
        int filledThisPass = 0;
        for (int ci = 0; ci < N; ++ci)
        {
            if (!phiZeroArtifact[static_cast<std::size_t>(ci)])
                continue;
            if (filled[static_cast<std::size_t>(ci)])
                continue;
            const int cLabel = mesh.cellIndexToLabel[ci];
            const Cell& cell = mesh.cells(cLabel);
            double sum = 0.0;
            int cnt = 0;
            for (int j = 0; j < 4; ++j)
            {
                const Vert& v = mesh.verts(cell.cellToVert[j]);
                for (int k = 0; k < 4; ++k)
                {
                    const int nLabel = v.vertToCell[k];
                    if (nLabel < 0) continue;
                    const int nci = mesh.cells(nLabel).index();
                    if (nci < 0 || nci >= N || nci == ci) continue;
                    if (outputCellValid[static_cast<std::size_t>(nci)] == 0u) continue;
                    // Accept non-artifact cells, or artifact cells already filled
                    if (phiZeroArtifact[static_cast<std::size_t>(nci)] &&
                        !filled[static_cast<std::size_t>(nci)])
                        continue;
                    sum += cOut.value[nci];
                    ++cnt;
                }
            }
            if (cnt > 0)
            {
                cOut.value[ci] = sum / static_cast<double>(cnt);
                filled[static_cast<std::size_t>(ci)] = true;
                ++filledThisPass;
            }
        }
        if (filledThisPass == 0)
            break;
    }
}

std::unique_ptr<VertField<double>> buildMaskedVertField(const QuadtreeMesh& mesh,
                                                               const CellField<double>& cellField,
                                                               const std::vector<unsigned char>* validCellMask)
{
    auto vertField = std::make_unique<VertField<double>>();
    constexpr double kSentinel = std::numeric_limits<double>::quiet_NaN();
    vertField->initialize(mesh, kSentinel);

    auto cellIndexIsValid = [&](int ci) -> bool {
        if (ci < 0 || ci >= mesh.cellIndexNum())
            return false;
        if (validCellMask == nullptr)
            return true;
        if (validCellMask->size() != static_cast<std::size_t>(mesh.cellIndexNum()))
            return true;
        return ((*validCellMask)[static_cast<std::size_t>(ci)] != 0u);
    };

    auto appendCellLabel = [&](const int cellLabel,
                               std::vector<int>& usedCellLabels,
                               std::vector<double>& values,
                               std::vector<Vector2D>& coords) {
        if (cellLabel < 0)
            return;
        const int ci = mesh.cells(cellLabel).index();
        if (!cellIndexIsValid(ci))
            return;
        for (const int usedCellLabel : usedCellLabels)
        {
            if (usedCellLabel == cellLabel)
                return;
        }
        usedCellLabels.push_back(cellLabel);
        values.push_back(cellField.value[ci]);
        coords.push_back(mesh.cells(cellLabel).centerCoord());
    };

    auto appendStencilVert = [&](const int vertLabel,
                                 std::vector<int>& usedCellLabels,
                                 std::vector<double>& values,
                                 std::vector<Vector2D>& coords) {
        if (vertLabel < 0)
            return;
        const Vert& stencilVert = mesh.verts(vertLabel);
        for (int j = 0; j < 4; ++j)
            appendCellLabel(stencilVert.vertToCell[j], usedCellLabels, values, coords);
    };

    for (int i = 0; i < mesh.vertIndexNum(); ++i)
    {
        const Vert& vert = mesh.verts(mesh.vertIndexToLabel[i]);
        std::vector<int> usedCellLabels;
        std::vector<double> values;
        std::vector<Vector2D> coords;

        auto appendOwnCells = [&]() {
            for (int j = 0; j < 4; ++j)
                appendCellLabel(vert.vertToCell[j], usedCellLabels, values, coords);
        };

        appendOwnCells();

        if (values.size() < 3 && vert.boundary() == Boundary::west)
        {
            appendStencilVert(vert.neiVertices[2], usedCellLabels, values, coords);
        }
        else if (values.size() < 3 && vert.boundary() == Boundary::south)
        {
            appendStencilVert(vert.neiVertices[3], usedCellLabels, values, coords);
        }
        else if (values.size() < 3 && vert.boundary() == Boundary::east)
        {
            appendStencilVert(vert.neiVertices[0], usedCellLabels, values, coords);
        }
        else if (values.size() < 3 && vert.boundary() == Boundary::north)
        {
            appendStencilVert(vert.neiVertices[1], usedCellLabels, values, coords);
        }
        else if (values.size() < 3 && vert.boundary() == Boundary::southWest)
        {
            if (vert.vertToCell[2] != -1)
                appendStencilVert(mesh.cells(vert.vertToCell[2]).cellToVert[2], usedCellLabels, values, coords);
        }
        else if (values.size() < 3 && vert.boundary() == Boundary::southEast)
        {
            if (vert.vertToCell[3] != -1)
                appendStencilVert(mesh.cells(vert.vertToCell[3]).cellToVert[3], usedCellLabels, values, coords);
        }
        else if (values.size() < 3 && vert.boundary() == Boundary::northEast)
        {
            if (vert.vertToCell[0] != -1)
                appendStencilVert(mesh.cells(vert.vertToCell[0]).cellToVert[0], usedCellLabels, values, coords);
        }
        else if (values.size() < 3 && vert.boundary() == Boundary::northWest)
        {
            if (vert.vertToCell[1] != -1)
                appendStencilVert(mesh.cells(vert.vertToCell[1]).cellToVert[1], usedCellLabels, values, coords);
        }

        if (values.empty() && vert.supCell != -1)
        {
            appendCellLabel(vert.supCell, usedCellLabels, values, coords);
        }

        if (coords.size() < 3 || coords.size() != values.size())
        {
            if (values.empty())
            {
                // No valid active-fluid stencil at this vertex. Search for
                // the nearest valid (masked) cells to avoid reading solid-
                // interior values (phi=0) that create visible patches near
                // EB aperture corners.
                double sum = 0.0;
                int cnt = 0;

                // 1-ring: direct neighbors, masked
                for (int j = 0; j < 4; ++j)
                {
                    const int cLabel = vert.vertToCell[j];
                    if (cLabel < 0)
                        continue;
                    const int ci = mesh.cells(cLabel).index();
                    if (!cellIndexIsValid(ci))
                        continue;
                    sum += cellField.value[ci];
                    ++cnt;
                }

                // 2-ring: neighbor-vertex cells, masked
                if (cnt == 0)
                {
                    for (int k = 0; k < 4; ++k)
                    {
                        const int nv = vert.neiVertices[k];
                        if (nv < 0)
                            continue;
                        const Vert& nVert = mesh.verts(nv);
                        for (int j = 0; j < 4; ++j)
                        {
                            const int cLabel = nVert.vertToCell[j];
                            if (cLabel < 0)
                                continue;
                            const int ci = mesh.cells(cLabel).index();
                            if (!cellIndexIsValid(ci))
                                continue;
                            sum += cellField.value[ci];
                            ++cnt;
                        }
                    }
                }

                // supCell masked fallback
                if (cnt == 0 && vert.supCell >= 0)
                {
                    const int ci = mesh.cells(vert.supCell).index();
                    if (cellIndexIsValid(ci))
                    {
                        sum += cellField.value[ci];
                        ++cnt;
                    }
                }

                // If still nothing, leave as NaN sentinel. The post-pass
                // below will fill it from valid cells that reference this
                // vertex via cellToVert (handles quadtree T-junctions where
                // vertToCell does not reach any active cell).
                vertField->value[i] = (cnt > 0) ? (sum / static_cast<double>(cnt)) : kSentinel;
            }
            else
            {
                double sum = 0.0;
                for (const double v : values)
                    sum += v;
                vertField->value[i] = sum / static_cast<double>(values.size());
            }
        }
        else
        {
            double interp = Interpolation::leastSquareInterpolation(coords, values, vert.coord(), false);
            double vMin = values.front();
            double vMax = values.front();
            double sum = 0.0;
            for (const double v : values)
            {
                vMin = std::min(vMin, v);
                vMax = std::max(vMax, v);
                sum += v;
            }
            if (!std::isfinite(interp))
                interp = sum / static_cast<double>(values.size());
            // For visualization output, keep nodal reconstruction bounded by the
            // contributing cell values. Moving-EB + AMR stencils can make the
            // unconstrained least-squares fit overshoot badly, which shows up as
            // spurious phi/Te spikes in packed Field plots.
            vertField->value[i] = std::clamp(interp, vMin, vMax);
        }
    }

    // Post-pass: fill remaining NaN sentinel vertices from valid cells
    // that reference them via cellToVert. This handles the case where a
    // vertex's vertToCell neighbourhood contains only inactive cells (common
    // at quadtree T-junctions near EB aperture corners) while an active cut
    // cell still uses it as a corner vertex for FEPOINT connectivity.
    for (int ci = 0; ci < mesh.cellIndexNum(); ++ci)
    {
        if (!cellIndexIsValid(ci))
            continue;
        const int cLabel = mesh.cellIndexToLabel[ci];
        const Cell& cell = mesh.cells(cLabel);
        for (int j = 0; j < 4; ++j)
        {
            const Vert& v = mesh.verts(cell.cellToVert[j]);
            const int vi = v.index();
            if (vi < 0 || vi >= mesh.vertIndexNum())
                continue;
            if (std::isnan(vertField->value[vi]))
                vertField->value[vi] = cellField.value[ci];
        }
    }

    // Final cleanup: any vertex still NaN (not referenced by any valid cell)
    // is deep inside solid — set to 0.
    for (int i = 0; i < mesh.vertIndexNum(); ++i)
    {
        if (std::isnan(vertField->value[i]))
            vertField->value[i] = 0.0;
    }

    return vertField;
}

void canonicalizePeriodicBoundaryVertices(const QuadtreeMesh& mesh,
                                                 const bool periodicX,
                                                 const bool periodicY,
                                                 std::vector<VertField<double>>& vars)
{
    if (vars.empty() || !(periodicX || periodicY))
        return;

    const double baseTol = std::max(1e-12, 1e-9 * std::min(mesh.cellSize.x, mesh.cellSize.y));

    if (periodicY)
    {
        std::unordered_map<std::int64_t, int> southByX;
        std::unordered_map<std::int64_t, int> northByX;
        southByX.reserve(static_cast<std::size_t>(mesh.vertIndexNum()));
        northByX.reserve(static_cast<std::size_t>(mesh.vertIndexNum()));

        for (int vi = 0; vi < mesh.vertIndexNum(); ++vi)
        {
            const Vert& vert = mesh.verts(mesh.vertIndexToLabel[vi]);
            const Vector2D p = vert.coord();
            const std::int64_t key = periodicCoordKey(p.x, baseTol);
            if (std::abs(p.y - mesh.minCoord.y) <= baseTol)
                southByX.emplace(key, vi);
            if (std::abs(p.y - mesh.maxCoord.y) <= baseTol)
                northByX.emplace(key, vi);
        }

        for (const auto& it : southByX)
        {
            const auto mate = northByX.find(it.first);
            if (mate == northByX.end())
                continue;
            const int south = it.second;
            const int north = mate->second;
            for (auto& vf : vars)
            {
                const double avg = 0.5 * (vf.value[south] + vf.value[north]);
                vf.value[south] = avg;
                vf.value[north] = avg;
            }
        }
    }

    if (periodicX)
    {
        std::unordered_map<std::int64_t, int> westByY;
        std::unordered_map<std::int64_t, int> eastByY;
        westByY.reserve(static_cast<std::size_t>(mesh.vertIndexNum()));
        eastByY.reserve(static_cast<std::size_t>(mesh.vertIndexNum()));

        for (int vi = 0; vi < mesh.vertIndexNum(); ++vi)
        {
            const Vert& vert = mesh.verts(mesh.vertIndexToLabel[vi]);
            const Vector2D p = vert.coord();
            const std::int64_t key = periodicCoordKey(p.y, baseTol);
            if (std::abs(p.x - mesh.minCoord.x) <= baseTol)
                westByY.emplace(key, vi);
            if (std::abs(p.x - mesh.maxCoord.x) <= baseTol)
                eastByY.emplace(key, vi);
        }

        for (const auto& it : westByY)
        {
            const auto mate = eastByY.find(it.first);
            if (mate == eastByY.end())
                continue;
            const int west = it.second;
            const int east = mate->second;
            for (auto& vf : vars)
            {
                const double avg = 0.5 * (vf.value[west] + vf.value[east]);
                vf.value[west] = avg;
                vf.value[east] = avg;
            }
        }
    }
}
}
