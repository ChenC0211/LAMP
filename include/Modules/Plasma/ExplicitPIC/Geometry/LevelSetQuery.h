#pragma once

#include <cmath>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include "Core/Mesh/Edge.h"

#include "Core/Geometry/ImmersedGeometry.h"
#include "Core/Base/Interpolation.h"

// ============================================================================
// LevelSetQuery
// ----------------------------------------------------------------------------
// Minimal helper utilities for querying the level-set value and its gradient
// at an arbitrary coordinate using *existing* data structures:
//   - ImmersedGeometry<Mesh>::levelSets (psi on vertices)
//   - ImmersedGeometry<Mesh>::LS_dx/LS_dy (dpsi/dx, dpsi/dy on vertices)
//   - Mesh::findCellLabel + Cell::cellToVert for bilinear interpolation
//
// This file intentionally does NOT modify existing geometry/mesh core classes.
// ============================================================================

namespace ExplicitPIC
{
namespace Geometry
{

inline const char* boundaryName(const Boundary b)
{
    switch (b)
    {
    case Boundary::non:       return "non";
    case Boundary::west:      return "west";
    case Boundary::south:     return "south";
    case Boundary::east:      return "east";
    case Boundary::north:     return "north";
    case Boundary::southWest: return "southWest";
    case Boundary::southEast: return "southEast";
    case Boundary::northEast: return "northEast";
    case Boundary::northWest: return "northWest";
    default:                  return "unknown";
    }
}

struct LevelSetQuery
{
    template<typename Mesh>
    static double psiAt(const Mesh& mesh, ImmersedGeometry<Mesh>& iGeometry, const int levelSetI, const Vector2D& x)
    {
        const int cLabel = mesh.findCellLabel(x);
        const Cell& c = mesh.cells(cLabel);
        auto& psi = iGeometry.levelSets.at(levelSetI);
        return Interpolation::biLinearInterpolation(mesh, psi, c, x);
    }

    template<typename Mesh>
    static Vector2D gradAt(const Mesh& mesh, ImmersedGeometry<Mesh>& iGeometry, const int levelSetI, const Vector2D& x)
    {
        const int cLabel = mesh.findCellLabel(x);
        const Cell& c = mesh.cells(cLabel);
        auto& dpsidx = iGeometry.LS_dx(levelSetI);
        auto& dpsidy = iGeometry.LS_dy(levelSetI);
        const double gx = Interpolation::biLinearInterpolation(mesh, dpsidx, c, x);
        const double gy = Interpolation::biLinearInterpolation(mesh, dpsidy, c, x);
        return Vector2D(gx, gy);
    }

    template<typename Mesh>
    static double gradMagAt(const Mesh& mesh, ImmersedGeometry<Mesh>& iGeometry, const int levelSetI, const Vector2D& x)
    {
        const Vector2D g = gradAt(mesh, iGeometry, levelSetI, x);
        return std::sqrt(g.x * g.x + g.y * g.y);
    }

    template<typename Mesh>
    static void debugSampleLeafCellsAndEdges(
        const Mesh& mesh,
        ImmersedGeometry<Mesh>& iGeometry,
        const int levelSetI,
        const int nCellSamples = 5,
        const int nEdgeSamples = 5,
        const unsigned seed = 12345u)
    {
        std::vector<int> leafCells;
        leafCells.reserve(mesh.cellIndexNum());
        for (int i = 0; i < mesh.cellIndexNum(); ++i)
        {
            const int cLabel = mesh.cellIndexToLabel[i];
            const Cell& c = mesh.cells(cLabel);
            if (c.son[0] == -1) leafCells.push_back(cLabel);
        }

        std::vector<int> leafHEdges;
        leafHEdges.reserve(mesh.hEdgeIndexNum());
        for (int i = 0; i < mesh.hEdgeIndexNum(); ++i)
        {
            const int eLabel = mesh.hEdgeIndexToLabel[i];
            const Edge& e = mesh.hEdges(eLabel);
            if (e.son[0] == -1) leafHEdges.push_back(eLabel);
        }

        std::vector<int> leafVEdges;
        leafVEdges.reserve(mesh.vEdgeIndexNum());
        for (int i = 0; i < mesh.vEdgeIndexNum(); ++i)
        {
            const int eLabel = mesh.vEdgeIndexToLabel[i];
            const Edge& e = mesh.vEdges(eLabel);
            if (e.son[0] == -1) leafVEdges.push_back(eLabel);
        }

        std::mt19937 rng(seed);

        std::cout << "\n[LevelSetQuery::debug] levelSetI=" << levelSetI
                  << ", leafCells=" << leafCells.size()
                  << ", leafHEdges=" << leafHEdges.size()
                  << ", leafVEdges=" << leafVEdges.size() << "\n";

        auto printHEdge = [&](const int eLabel, const char* tag)
        {
            if (eLabel < 0) return;
            const Edge& e = mesh.hEdges(eLabel);
            std::cout << "  " << tag
                      << " HEdge(label=" << e.label()
                      << ", idx=" << e.index()
                      << ", bnd=" << boundaryName(e.boundary())
                      << ", len=" << e.length()
                      << ", cutLen=" << e.cutLength()
                      << ", cut=" << (e.cut() ? 1 : 0)
                      << ", isCut[ls]=" << (e.isCut[levelSetI] ? 1 : 0)
                      << ", master=" << e.masterCell
                      << ", slave=" << e.slaveCell
                      << ")\n";
        };

        auto printVEdge = [&](const int eLabel, const char* tag)
        {
            if (eLabel < 0) return;
            const Edge& e = mesh.vEdges(eLabel);
            std::cout << "  " << tag
                      << " VEdge(label=" << e.label()
                      << ", idx=" << e.index()
                      << ", bnd=" << boundaryName(e.boundary())
                      << ", len=" << e.length()
                      << ", cutLen=" << e.cutLength()
                      << ", cut=" << (e.cut() ? 1 : 0)
                      << ", isCut[ls]=" << (e.isCut[levelSetI] ? 1 : 0)
                      << ", master=" << e.masterCell
                      << ", slave=" << e.slaveCell
                      << ")\n";
        };

        // ---- Sample leaf cells (print psi(center), |grad|, and its 4 sides' edges) ----
        if (!leafCells.empty())
        {
            std::uniform_int_distribution<int> distCell(0, static_cast<int>(leafCells.size()) - 1);
            const int n = std::min(nCellSamples, static_cast<int>(leafCells.size()));
            std::cout << "\n[LevelSetQuery::debug] Sample leaf cells (n=" << n << "):\n";
            for (int s = 0; s < n; ++s)
            {
                const int cLabel = leafCells[distCell(rng)];
                const Cell& c = mesh.cells(cLabel);
                const Vector2D xc = c.centerCoord();
                const double psiC = psiAt(mesh, iGeometry, levelSetI, xc);
                const double gmag = gradMagAt(mesh, iGeometry, levelSetI, xc);
                std::cout << "\nCell(label=" << c.label() << ", idx=" << c.index() << ", level=" << c.level()
                          << ", bnd=" << boundaryName(c.boundary())
                          << ", center=(" << xc.x << "," << xc.y << ")"
                          << ", psi=" << psiC
                          << ", |grad|=" << gmag
                          << ", isCut[ls]=" << (c.isCut[levelSetI] ? 1 : 0)
                          << ")\n";

                // bottom/top H edges (row 0: south, row 1: north)
                printHEdge(c.cellToHEdge[0][0], "south[0]");
                printHEdge(c.cellToHEdge[0][1], "south[1]");
                printHEdge(c.cellToHEdge[1][0], "north[0]");
                printHEdge(c.cellToHEdge[1][1], "north[1]");

                // left/right V edges (row 0: west, row 1: east)
                printVEdge(c.cellToVEdge[0][0], "west[0] ");
                printVEdge(c.cellToVEdge[0][1], "west[1] ");
                printVEdge(c.cellToVEdge[1][0], "east[0] ");
                printVEdge(c.cellToVEdge[1][1], "east[1] ");
            }
        }

        // ---- Sample leaf edges (print length/cutLength) ----
        const int nH = std::min(nEdgeSamples, static_cast<int>(leafHEdges.size()));
        if (nH > 0)
        {
            std::uniform_int_distribution<int> distE(0, static_cast<int>(leafHEdges.size()) - 1);
            std::cout << "\n[LevelSetQuery::debug] Sample leaf HEdges (n=" << nH << "):\n";
            for (int s = 0; s < nH; ++s)
            {
                const int eLabel = leafHEdges[distE(rng)];
                printHEdge(eLabel, "sample");
            }
        }

        const int nV = std::min(nEdgeSamples, static_cast<int>(leafVEdges.size()));
        if (nV > 0)
        {
            std::uniform_int_distribution<int> distE(0, static_cast<int>(leafVEdges.size()) - 1);
            std::cout << "\n[LevelSetQuery::debug] Sample leaf VEdges (n=" << nV << "):\n";
            for (int s = 0; s < nV; ++s)
            {
                const int eLabel = leafVEdges[distE(rng)];
                printVEdge(eLabel, "sample");
            }
        }

        std::cout << "\n[LevelSetQuery::debug] Done.\n";
    }
};

} // namespace Geometry
} // namespace ExplicitPIC
