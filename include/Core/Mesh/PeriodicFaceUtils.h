#pragma once

#include <algorithm>
#include <cmath>

#include "Core/Mesh/FaceSegmentGraph.h"
#include "Core/Mesh/QuadtreeMesh.h"

namespace PeriodicFaceUtils
{

inline Boundary mateBoundary(const Boundary b)
{
    switch (b)
    {
    case Boundary::west:  return Boundary::east;
    case Boundary::east:  return Boundary::west;
    case Boundary::south: return Boundary::north;
    case Boundary::north: return Boundary::south;
    default:              return Boundary::non;
    }
}

inline Vector2D shiftToMate(const QuadtreeMesh& mesh, const Boundary b)
{
    const double Lx = mesh.maxCoord.x - mesh.minCoord.x;
    const double Ly = mesh.maxCoord.y - mesh.minCoord.y;

    switch (b)
    {
    case Boundary::west:  return Vector2D(+Lx, 0.0);
    case Boundary::east:  return Vector2D(-Lx, 0.0);
    case Boundary::south: return Vector2D(0.0, +Ly);
    case Boundary::north: return Vector2D(0.0, -Ly);
    default:              return Vector2D(0.0, 0.0);
    }
}

inline Vector2D shiftToMate(const QuadtreeMesh& mesh, const FaceSide side)
{
    return shiftToMate(mesh, faceSideToBoundary(side));
}

inline double probeEpsilon(const QuadtreeMesh& mesh)
{
    return 1e-9 * std::max(1e-12, std::min(mesh.cellSize.x, mesh.cellSize.y));
}

inline Vector2D clampProbeToDomain(const QuadtreeMesh& mesh, Vector2D p)
{
    const double eps = probeEpsilon(mesh);
    p.x = std::clamp(p.x, mesh.minCoord.x + eps, mesh.maxCoord.x - eps);
    p.y = std::clamp(p.y, mesh.minCoord.y + eps, mesh.maxCoord.y - eps);
    return p;
}

inline Vector2D mateProbePoint(const QuadtreeMesh& mesh,
                               const FaceSegment& seg,
                               const Vector2D& shiftToMateVec)
{
    const double eps = probeEpsilon(mesh);
    return clampProbeToDomain(mesh, seg.faceCenter + shiftToMateVec + seg.n * eps);
}

inline bool findMateCellLabel(const QuadtreeMesh& mesh,
                              const FaceSegment& seg,
                              const Vector2D& shiftToMateVec,
                              int& mateLabel)
{
    mateLabel = mesh.findCellLabel(mateProbePoint(mesh, seg, shiftToMateVec));
    if (mateLabel >= 0)
        return true;

    const double eps = probeEpsilon(mesh);
    if (seg.K >= 0 && std::isfinite(seg.df) && (seg.df > eps))
    {
        const Vector2D cK = mesh.cells(seg.K).centerCoord();
        mateLabel = mesh.findCellLabel(clampProbeToDomain(mesh, cK + shiftToMateVec + seg.n * (2.0 * seg.df)));
    }
    return (mateLabel >= 0);
}

inline int findMateCellLabel(const QuadtreeMesh& mesh,
                             const FaceSegment& seg,
                             const Vector2D& shiftToMateVec)
{
    int mateLabel = -1;
    findMateCellLabel(mesh, seg, shiftToMateVec, mateLabel);
    return mateLabel;
}

inline double centerDistance(const QuadtreeMesh& mesh,
                             const FaceSegment& seg,
                             const int mateLabel,
                             const Vector2D& shiftToMateVec)
{
    const Vector2D cK = mesh.cells(seg.K).centerCoord();
    const Vector2D cL = mesh.cells(mateLabel).centerCoord();
    const Vector2D cL_inK = cL - shiftToMateVec;
    return std::abs((cL_inK - cK).dotProduct(seg.n));
}

inline double tangentialCenterOffset(const QuadtreeMesh& mesh,
                                     const FaceSegment& seg,
                                     const int mateLabel,
                                     const Vector2D& shiftToMateVec)
{
    const Vector2D cK = mesh.cells(seg.K).centerCoord();
    const Vector2D cL = mesh.cells(mateLabel).centerCoord();
    const Vector2D cL_inK = cL - shiftToMateVec;
    const Vector2D t(-seg.n.y, seg.n.x);
    return std::abs((cL_inK - cK).dotProduct(t));
}

inline bool isCleanCellMate(const QuadtreeMesh& mesh,
                            const FaceSegment& seg,
                            const int mateLabel,
                            const Vector2D& shiftToMateVec)
{
    if (seg.K < 0 || mateLabel < 0)
        return false;

    const Cell& cellK = mesh.cells(seg.K);
    const Cell& cellL = mesh.cells(mateLabel);
    const Vector2D sizeK = cellK.size();
    const Vector2D sizeL = cellL.size();
    const bool periodicY = (seg.boundary == Boundary::south || seg.boundary == Boundary::north);
    const double tangentialSizeK = periodicY ? sizeK.x : sizeK.y;
    const double tangentialSizeL = periodicY ? sizeL.x : sizeL.y;
    const double tangentialTol = 1e-9 * std::max(1e-12, std::min(tangentialSizeK, tangentialSizeL));
    return tangentialCenterOffset(mesh, seg, mateLabel, shiftToMateVec) <= tangentialTol;
}

} // namespace PeriodicFaceUtils
