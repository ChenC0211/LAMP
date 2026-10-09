#pragma once

#include <vector>

#include "Core/Base/Vector2D.h"
#include "Core/Field/VertField.h"

class BackgroundMesh;
class Cell;

// ----------------------------------------------------------------------------
// CICMomentExactCut
//
// Moment-exact cut weighting for CIC deposition/gather on cut cells.
//
// Goal:
//   Replace the crude truncation W <- alpha_K * W on cut cells by a geometric
//   integration of the CIC kernel over the *fluid* portion of the cell.
//
// Model used here (minimal, robust):
//   - For a cut cell, we build a planar approximation of the level set
//     phi(x,y) inside the cell using vertex samples.
//   - The fluid region is approximated by a half-plane (phi >= 0 or phi <= 0,
//     controlled by insideIsNegative).
//   - We clip the (cell ∩ CIC support) rectangle by that half-plane.
//   - We then integrate the CIC kernel exactly over the clipped polygon by
//     splitting at x=xp and y=yp (to handle |.| in the hat function) and using
//     polygon moments (A, \int x, \int y, \int xy).
//
// This is "moment-exact" w.r.t. the planar cut and the bilinear CIC kernel.
// ----------------------------------------------------------------------------

namespace CICMomentExactCut {

struct Plane
{
    // plane: a*x + b*y + c = 0
    double a = 0.0;
    double b = 0.0;
    double c = 0.0;
};

// Build a planar approximation of the level set in a cell using vertex samples.
// Returns false on failure (invalid indices, degenerate geometry, etc.).
bool buildCellPlaneFromVertexLevelSet(const BackgroundMesh& mesh,
                                      const Cell& cell,
                                      const VertField<double>& levelSetVert,
                                      Plane* outPlane);

// Compute moment-exact CIC weight over an arbitrary axis-aligned rectangle
// clipped by a half-plane defined by a planar level-set approximation.
//
//   W = \int_{rect \cap support \cap fluid(pl)} S(x,y)\,dA,
//
// where S is the CIC kernel (this module). The rectangle is typically a sub-rectangle
// of a (possibly coarse) leaf cell when using virtual partitioning.
//
// Returns true on success (including the case of empty overlap, where outW=0).
// Returns false on failure; callers should fall back to alpha-truncation.
bool cutRectWeight2DFromPlane(double rect_xL,
                              double rect_xR,
                              double rect_yL,
                              double rect_yR,
                              const Plane& pl,
                              double xp,
                              double yp,
                              double hp,
                              bool insideIsNegative,
                              double* outW);

// Convenience wrapper: builds the cell plane from a vertex level-set and then
// evaluates cutRectWeight2DFromPlane on the provided rectangle.
bool cutRectWeight2D(const BackgroundMesh& mesh,
                     int cellLabel,
                     const VertField<double>& levelSetVert,
                     double rect_xL,
                     double rect_xR,
                     double rect_yL,
                     double rect_yR,
                     double xp,
                     double yp,
                     double hp,
                     bool insideIsNegative,
                     double* outW);

// Compute moment-exact CIC weight on a cut cell:
//   W = \int_{K \cap support \cap fluid} S(x,y) dA
// where S is the this module CIC kernel.
//
// Returns true and writes outW on success. Returns false if the computation is
// not applicable (degenerate hp, empty overlap) or fails; callers should fall
// back to alpha-truncation.
bool cutCellWeight2D(const BackgroundMesh& mesh,
                     int cellLabel,
                     const VertField<double>& levelSetVert,
                     double xp,
                     double yp,
                     double hp,
                     bool insideIsNegative,
                     double* outW);

} // namespace CICMomentExactCut
