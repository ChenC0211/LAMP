#include "Modules/Plasma/ExplicitPIC/Deposit/CICMomentExactCut.h"

#include "Core/Mesh/BackgroundMesh.h"
#include "Core/Mesh/Cell.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace CICMomentExactCut {

namespace {

inline double cross2(const Vector2D& a, const Vector2D& b)
{
    return a.x * b.y - a.y * b.x;
}

struct HalfPlane
{
    // Keep points where (nx * x + ny * y + d) >= -eps.
    double nx = 0.0;
    double ny = 0.0;
    double d = 0.0;
    double eps = 1e-14;

    inline double eval(const Vector2D& p) const { return nx * p.x + ny * p.y + d; }
    inline bool inside(const Vector2D& p) const { return eval(p) >= -eps; }
};

std::vector<Vector2D> clipPolygonHalfPlane(const std::vector<Vector2D>& poly, const HalfPlane& hp)
{
    std::vector<Vector2D> out;
    if (poly.size() < 3)
        return out;

    out.reserve(poly.size() + 2);
    Vector2D S = poly.back();
    bool inS = hp.inside(S);
    double fS = hp.eval(S);

    for (const Vector2D& E : poly)
    {
        const bool inE = hp.inside(E);
        const double fE = hp.eval(E);

        if (inE)
        {
            if (!inS)
            {
                const double denom = (fS - fE);
                if (std::abs(denom) > 1e-30)
                {
                    const double t = fS / denom; // in [0,1] for proper crossings
                    out.emplace_back(Vector2D{S.x + t * (E.x - S.x), S.y + t * (E.y - S.y)});
                }
                else
                {
                    // Nearly parallel; best-effort: push midpoint.
                    out.emplace_back(Vector2D{0.5 * (S.x + E.x), 0.5 * (S.y + E.y)});
                }
            }
            out.emplace_back(E);
        }
        else if (inS)
        {
            const double denom = (fS - fE);
            if (std::abs(denom) > 1e-30)
            {
                const double t = fS / denom;
                out.emplace_back(Vector2D{S.x + t * (E.x - S.x), S.y + t * (E.y - S.y)});
            }
            else
            {
                out.emplace_back(Vector2D{0.5 * (S.x + E.x), 0.5 * (S.y + E.y)});
            }
        }

        S = E;
        inS = inE;
        fS = fE;
    }

    // Remove tiny duplicate endpoints (can occur when clipping hits a vertex).
    if (out.size() >= 2)
    {
        const Vector2D& a = out.front();
        const Vector2D& b = out.back();
        if (std::abs(a.x - b.x) <= 1e-14 && std::abs(a.y - b.y) <= 1e-14)
            out.pop_back();
    }

    return out;
}

inline double signedArea(const std::vector<Vector2D>& poly)
{
    if (poly.size() < 3)
        return 0.0;
    double s = 0.0;
    for (std::size_t i = 0; i < poly.size(); ++i)
    {
        const Vector2D& p = poly[i];
        const Vector2D& q = poly[(i + 1) % poly.size()];
        s += (p.x * q.y - q.x * p.y);
    }
    return 0.5 * s;
}

struct Moments
{
    double A = 0.0;   // area
    double Mx = 0.0;  // \int x dA
    double My = 0.0;  // \int y dA
    double Mxy = 0.0; // \int x y dA
};

Moments polygonMoments(std::vector<Vector2D> poly)
{
    Moments m;
    if (poly.size() < 3)
        return m;

    // Ensure CCW for positive signed triangle areas.
    const double A0 = signedArea(poly);
    if (A0 < 0.0)
        std::reverse(poly.begin(), poly.end());

    const Vector2D p0 = poly.front();
    for (std::size_t i = 1; i + 1 < poly.size(); ++i)
    {
        const Vector2D p1 = poly[i];
        const Vector2D p2 = poly[i + 1];
        const double tri2 = cross2(Vector2D{p1.x - p0.x, p1.y - p0.y},
                                   Vector2D{p2.x - p0.x, p2.y - p0.y});
        const double A = 0.5 * tri2;
        if (!(A > 0.0))
            continue;

        m.A += A;
        m.Mx += A * (p0.x + p1.x + p2.x) / 3.0;
        m.My += A * (p0.y + p1.y + p2.y) / 3.0;

        const double diag = p0.x * p0.y + p1.x * p1.y + p2.x * p2.y;
        const double off = p0.x * p1.y + p0.x * p2.y + p1.x * p0.y + p1.x * p2.y + p2.x * p0.y + p2.x * p1.y;
        m.Mxy += A * (diag / 6.0 + off / 12.0);
    }

    return m;
}

inline double integrateCICOverPolygonQuadrant(const std::vector<Vector2D>& poly,
                                             double xp,
                                             double yp,
                                             double hp,
                                             bool isRight,
                                             bool isTop)
{
    if (poly.size() < 3 || !(hp > 0.0))
        return 0.0;

    const Moments m = polygonMoments(poly);
    if (!(m.A > 0.0))
        return 0.0;

    // s_x and s_y are linear on each quadrant.
    const double invh = 1.0 / hp;

    double ax = 0.0, bx = 0.0;
    if (isRight)
    {
        // x >= xp: s = (xp + hp - x) / hp
        ax = -invh;
        bx = (xp + hp) * invh;
    }
    else
    {
        // x <= xp: s = (x - (xp - hp)) / hp
        ax = invh;
        bx = -(xp - hp) * invh;
    }

    double ay = 0.0, by = 0.0;
    if (isTop)
    {
        // y >= yp: s = (yp + hp - y) / hp
        ay = -invh;
        by = (yp + hp) * invh;
    }
    else
    {
        // y <= yp: s = (y - (yp - hp)) / hp
        ay = invh;
        by = -(yp - hp) * invh;
    }

    const double I = ax * ay * m.Mxy + ax * by * m.Mx + ay * bx * m.My + bx * by * m.A;
    const double W = I * invh * invh; // divide by hp^2

    // Numerical noise can create tiny negative values.
    return (W > 0.0) ? W : 0.0;
}

} // namespace

bool cutRectWeight2DFromPlane(double rect_xL,
                              double rect_xR,
                              double rect_yL,
                              double rect_yR,
                              const Plane& pl,
                              double xp,
                              double yp,
                              double hp,
                              bool insideIsNegative,
                              double* outW)
{
    if (!outW)
        return false;
    *outW = 0.0;

    if (!(hp > 0.0))
        return false;
    if (!(rect_xR > rect_xL) || !(rect_yR > rect_yL))
        return true; // empty rect

    // Intersect rect with CIC support.
    const double sL = xp - hp;
    const double sR = xp + hp;
    const double tL = yp - hp;
    const double tR = yp + hp;

    const double ixL = std::max(rect_xL, sL);
    const double ixR = std::min(rect_xR, sR);
    const double iyL = std::max(rect_yL, tL);
    const double iyR = std::min(rect_yR, tR);

    if (!(ixR > ixL) || !(iyR > iyL))
        return true; // empty overlap

    // Start polygon = (rect ∩support) rectangle, CCW.
    std::vector<Vector2D> poly;
    poly.reserve(8);
    poly.emplace_back(Vector2D{ixL, iyL});
    poly.emplace_back(Vector2D{ixR, iyL});
    poly.emplace_back(Vector2D{ixR, iyR});
    poly.emplace_back(Vector2D{ixL, iyR});

    // Clip by fluid half-plane.
    // Convention: insideIsNegative==true means "solid" is phi<0, so fluid is phi>=0.
    const double sgn = insideIsNegative ? 1.0 : -1.0;
    HalfPlane hpFluid;
    hpFluid.nx = sgn * pl.a;
    hpFluid.ny = sgn * pl.b;
    hpFluid.d = sgn * pl.c;
    hpFluid.eps = 1e-14;

    poly = clipPolygonHalfPlane(poly, hpFluid);
    if (poly.size() < 3)
        return true; // fully excluded

    auto clipXLE = [&](const std::vector<Vector2D>& in, double xMax) {
        HalfPlane hp2;
        hp2.nx = -1.0;
        hp2.ny = 0.0;
        hp2.d = xMax;
        hp2.eps = 1e-14;
        return clipPolygonHalfPlane(in, hp2);
    };
    auto clipXGE = [&](const std::vector<Vector2D>& in, double xMin) {
        HalfPlane hp2;
        hp2.nx = 1.0;
        hp2.ny = 0.0;
        hp2.d = -xMin;
        hp2.eps = 1e-14;
        return clipPolygonHalfPlane(in, hp2);
    };
    auto clipYLE = [&](const std::vector<Vector2D>& in, double yMax) {
        HalfPlane hp2;
        hp2.nx = 0.0;
        hp2.ny = -1.0;
        hp2.d = yMax;
        hp2.eps = 1e-14;
        return clipPolygonHalfPlane(in, hp2);
    };
    auto clipYGE = [&](const std::vector<Vector2D>& in, double yMin) {
        HalfPlane hp2;
        hp2.nx = 0.0;
        hp2.ny = 1.0;
        hp2.d = -yMin;
        hp2.eps = 1e-14;
        return clipPolygonHalfPlane(in, hp2);
    };

    // LL: x<=xp, y<=yp
    {
        std::vector<Vector2D> q = clipXLE(poly, xp);
        q = clipYLE(q, yp);
        *outW += integrateCICOverPolygonQuadrant(q, xp, yp, hp, false, false);
    }
    // LR: x>=xp, y<=yp
    {
        std::vector<Vector2D> q = clipXGE(poly, xp);
        q = clipYLE(q, yp);
        *outW += integrateCICOverPolygonQuadrant(q, xp, yp, hp, true, false);
    }
    // UL: x<=xp, y>=yp
    {
        std::vector<Vector2D> q = clipXLE(poly, xp);
        q = clipYGE(q, yp);
        *outW += integrateCICOverPolygonQuadrant(q, xp, yp, hp, false, true);
    }
    // UR: x>=xp, y>=yp
    {
        std::vector<Vector2D> q = clipXGE(poly, xp);
        q = clipYGE(q, yp);
        *outW += integrateCICOverPolygonQuadrant(q, xp, yp, hp, true, true);
    }

    if (!std::isfinite(*outW))
        return false;

    // Cap to an upper bound based on uncut overlap of this rectangle.
    const double wMax = std::max(0.0, std::min(1.0, (ixR - ixL) * (iyR - iyL) / (hp * hp)));
    if (*outW < 0.0)
        *outW = 0.0;
    if (*outW > wMax * (1.0 + 1e-10))
        *outW = std::min(*outW, wMax);

    return true;
}

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
                     double* outW)
{
    if (!outW)
        return false;
    *outW = 0.0;
    if (cellLabel < 0)
        return false;

    const Cell& cell = mesh.cells(cellLabel);
    Plane pl;
    if (!buildCellPlaneFromVertexLevelSet(mesh, cell, levelSetVert, &pl))
        return false;

    return cutRectWeight2DFromPlane(rect_xL, rect_xR, rect_yL, rect_yR, pl, xp, yp, hp, insideIsNegative, outW);
}

bool buildCellPlaneFromVertexLevelSet(const BackgroundMesh& mesh,
                                      const Cell& cell,
                                      const VertField<double>& levelSetVert,
                                      Plane* outPlane)
{
    if (!outPlane)
        return false;

    // Vertex ordering: 0=BL,1=BR,2=TR,3=TL
    const int v0 = cell.cellToVert[0];
    const int v1 = cell.cellToVert[1];
    const int v2 = cell.cellToVert[2];
    const int v3 = cell.cellToVert[3];

    if (v0 < 0 || v1 < 0 || v2 < 0 || v3 < 0)
        return false;

    const int i0 = mesh.verts(v0).index();
    const int i1 = mesh.verts(v1).index();
    const int i2 = mesh.verts(v2).index();
    const int i3 = mesh.verts(v3).index();

    if (i0 < 0 || i1 < 0 || i2 < 0 || i3 < 0)
        return false;
    if (i0 >= levelSetVert.size() || i1 >= levelSetVert.size() || i2 >= levelSetVert.size() || i3 >= levelSetVert.size())
        return false;

    const double phi0 = levelSetVert.value[i0];
    const double phi1 = levelSetVert.value[i1];
    const double phi2 = levelSetVert.value[i2];
    const double phi3 = levelSetVert.value[i3];

    if (!std::isfinite(phi0) || !std::isfinite(phi1) || !std::isfinite(phi2) || !std::isfinite(phi3))
        return false;

    const Vector2D cc = cell.centerCoord();
    const Vector2D sz = cell.size();
    if (!(sz.x > 0.0) || !(sz.y > 0.0))
        return false;

    // Cell-center value from bilinear midpoint (same as vertex average).
    const double phiC = 0.25 * (phi0 + phi1 + phi2 + phi3);

    // Central-difference gradient estimate.
    const double gx = ((phi1 + phi2) - (phi0 + phi3)) / (2.0 * sz.x);
    const double gy = ((phi3 + phi2) - (phi0 + phi1)) / (2.0 * sz.y);

    if (!std::isfinite(gx) || !std::isfinite(gy))
        return false;

    const double gnorm = std::sqrt(gx * gx + gy * gy);
    if (!(gnorm > 0.0) || gnorm < 1e-14)
        return false;

    // Plane: gx*(x-cc.x) + gy*(y-cc.y) + phiC = 0
    outPlane->a = gx;
    outPlane->b = gy;
    outPlane->c = phiC - gx * cc.x - gy * cc.y;

    return true;
}

bool cutCellWeight2D(const BackgroundMesh& mesh,
                     int cellLabel,
                     const VertField<double>& levelSetVert,
                     double xp,
                     double yp,
                     double hp,
                     bool insideIsNegative,
                     double* outW)
{
    if (!outW)
        return false;
    *outW = 0.0;

    if (!(hp > 0.0) || cellLabel < 0)
        return false;

    const Cell& cell = mesh.cells(cellLabel);

    // Build plane approximation.
    Plane pl;
    if (!buildCellPlaneFromVertexLevelSet(mesh, cell, levelSetVert, &pl))
        return false;

    // Cell bounds.
    const Vector2D cc = cell.centerCoord();
    const Vector2D sz = cell.size();
    const double xL = cc.x - 0.5 * sz.x;
    const double xR = cc.x + 0.5 * sz.x;
    const double yL = cc.y - 0.5 * sz.y;
    const double yR = cc.y + 0.5 * sz.y;

    // Intersect with CIC support.
    const double sL = xp - hp;
    const double sR = xp + hp;
    const double tL = yp - hp;
    const double tR = yp + hp;

    const double ixL = std::max(xL, sL);
    const double ixR = std::min(xR, sR);
    const double iyL = std::max(yL, tL);
    const double iyR = std::min(yR, tR);

    if (!(ixR > ixL) || !(iyR > iyL))
        return true; // empty overlap -> weight 0

    // Start polygon = (cell ∩support) rectangle, CCW.
    std::vector<Vector2D> poly;
    poly.reserve(8);
    poly.emplace_back(Vector2D{ixL, iyL});
    poly.emplace_back(Vector2D{ixR, iyL});
    poly.emplace_back(Vector2D{ixR, iyR});
    poly.emplace_back(Vector2D{ixL, iyR});

    // Clip by fluid half-plane.
    // Convention: insideIsNegative==true means "solid" is phi<0, so fluid is phi>=0.
    const double sign = insideIsNegative ? 1.0 : -1.0;
    HalfPlane hpFluid;
    hpFluid.nx = sign * pl.a;
    hpFluid.ny = sign * pl.b;
    hpFluid.d = sign * pl.c;
    hpFluid.eps = 1e-14;

    poly = clipPolygonHalfPlane(poly, hpFluid);
    if (poly.size() < 3)
        return true; // fully excluded

    // Split by x=xp and y=yp (quadrants) and integrate exactly.
    auto clipXLE = [&](const std::vector<Vector2D>& in, double xMax) {
        HalfPlane hp;
        hp.nx = -1.0;
        hp.ny = 0.0;
        hp.d = xMax;
        hp.eps = 1e-14;
        return clipPolygonHalfPlane(in, hp);
    };
    auto clipXGE = [&](const std::vector<Vector2D>& in, double xMin) {
        HalfPlane hp;
        hp.nx = 1.0;
        hp.ny = 0.0;
        hp.d = -xMin;
        hp.eps = 1e-14;
        return clipPolygonHalfPlane(in, hp);
    };
    auto clipYLE = [&](const std::vector<Vector2D>& in, double yMax) {
        HalfPlane hp;
        hp.nx = 0.0;
        hp.ny = -1.0;
        hp.d = yMax;
        hp.eps = 1e-14;
        return clipPolygonHalfPlane(in, hp);
    };
    auto clipYGE = [&](const std::vector<Vector2D>& in, double yMin) {
        HalfPlane hp;
        hp.nx = 0.0;
        hp.ny = 1.0;
        hp.d = -yMin;
        hp.eps = 1e-14;
        return clipPolygonHalfPlane(in, hp);
    };

    // LL: x<=xp, y<=yp
    {
        std::vector<Vector2D> q = clipXLE(poly, xp);
        q = clipYLE(q, yp);
        *outW += integrateCICOverPolygonQuadrant(q, xp, yp, hp, false, false);
    }
    // LR: x>=xp, y<=yp
    {
        std::vector<Vector2D> q = clipXGE(poly, xp);
        q = clipYLE(q, yp);
        *outW += integrateCICOverPolygonQuadrant(q, xp, yp, hp, true, false);
    }
    // UL: x<=xp, y>=yp
    {
        std::vector<Vector2D> q = clipXLE(poly, xp);
        q = clipYGE(q, yp);
        *outW += integrateCICOverPolygonQuadrant(q, xp, yp, hp, false, true);
    }
    // UR: x>=xp, y>=yp
    {
        std::vector<Vector2D> q = clipXGE(poly, xp);
        q = clipYGE(q, yp);
        *outW += integrateCICOverPolygonQuadrant(q, xp, yp, hp, true, true);
    }

    if (!std::isfinite(*outW))
        return false;

    // Guard: tiny numerical overshoot above the uncut overlap is possible; cap it.
    const double wMax = std::max(0.0, std::min(1.0, (ixR - ixL) * (iyR - iyL) / (hp * hp)));
    if (*outW < 0.0)
        *outW = 0.0;
    if (*outW > wMax * (1.0 + 1e-10))
        *outW = std::min(*outW, wMax);

    return true;
}

} // namespace CICMomentExactCut
