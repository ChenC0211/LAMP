#pragma once

// Shared arithmetic for one complete, same-level native RZ momentum stencil.
// No mesh numbering, solid mask, boundary closure, pressure, or hoop correction
// is inferred here. Callers supply actual indices, physical RZ measures, and
// both velocity histories. Missing/coarse-fine stencils keep their old path.
namespace HydroNativeRZRow
{
struct Indices
{
    int center = -1, west = -1, east = -1, south = -1, north = -1;
    // The OTHER velocity component at the four geometric stencil corners.
    int crossSW = -1, crossSE = -1, crossNW = -1, crossNE = -1;
};

struct Metrics
{
    double volume = 0.0, areaWest = 0.0, areaEast = 0.0, areaNS = 0.0;
    double dx = 0.0, dy = 0.0;
    // Actual native midpoint differences used by the old explicit shear.
    // Radial row: south/north distances; axial row: west/east distances.
    // Matrix arithmetic intentionally retains its original dx/dy quotients.
    double shearLowDistance = 0.0, shearHighDistance = 0.0;
};

struct Materials
{
    double rho = 0.0;
    double rhoWest = 0.0, rhoEast = 0.0, rhoSouth = 0.0, rhoNorth = 0.0;
    double muWest = 0.0, muEast = 0.0, muSouth = 0.0, muNorth = 0.0;
};

struct Row
{
    Indices indices;
    Metrics metrics;
    Materials materials;
};

struct Values
{
    double center = 0.0, west = 0.0, east = 0.0, south = 0.0, north = 0.0;
    double crossSW = 0.0, crossSE = 0.0, crossNW = 0.0, crossNE = 0.0;
};

template<typename Read>
Values values(const Indices& i, Read read)
{
    return {read(i.center), read(i.west), read(i.east), read(i.south), read(i.north),
        read(i.crossSW), read(i.crossSE), read(i.crossNW), read(i.crossNE)};
}

// emit(row, column, coefficient) is called in the historical triplet order.
// Cross shear is part of this CN row; it must not be omitted by a fitted caller.
template<typename Emit>
void appendRadialMatrix(const Row& row, double dt, Emit emit)
{
    const auto& i = row.indices;
    const auto& g = row.metrics;
    const auto& m = row.materials;
    const double AreaNS_divDy = g.areaNS / g.dy;
    const double AreaE_divDx = g.areaEast / g.dx;
    const double AreaW_divDx = g.areaWest / g.dx;
    emit(i.center, i.south, -0.5 * m.muSouth * AreaNS_divDy);
    emit(i.center, i.north, -0.5 * m.muNorth * AreaNS_divDy);
    emit(i.center, i.east, -m.muEast * AreaE_divDx);
    emit(i.center, i.crossSE, 0.5 * m.muSouth * g.areaNS / g.dx);
    emit(i.center, i.crossNE, -0.5 * m.muNorth * g.areaNS / g.dx);
    emit(i.center, i.west, -m.muWest * AreaW_divDx);
    emit(i.center, i.crossSW, -0.5 * m.muSouth * g.areaNS / g.dx);
    emit(i.center, i.crossNW, 0.5 * m.muNorth * g.areaNS / g.dx);
    emit(i.center, i.center, g.volume * m.rho / dt
        + (m.muWest * AreaW_divDx + m.muEast * AreaE_divDx)
        + 0.5 * (m.muSouth * AreaNS_divDy + m.muNorth * AreaNS_divDy));
}

template<typename Emit>
void appendAxialMatrix(const Row& row, double dt, Emit emit)
{
    const auto& i = row.indices;
    const auto& g = row.metrics;
    const auto& m = row.materials;
    const double AreaNS_divDy = g.areaNS / g.dy;
    const double AreaE_divDx = g.areaEast / g.dx;
    const double AreaW_divDx = g.areaWest / g.dx;
    emit(i.center, i.west, -0.5 * m.muWest * AreaW_divDx);
    emit(i.center, i.east, -0.5 * m.muEast * AreaE_divDx);
    emit(i.center, i.north, -m.muNorth * AreaNS_divDy);
    emit(i.center, i.crossNW, 0.5 * m.muWest * g.areaWest / g.dy);
    emit(i.center, i.crossNE, -0.5 * m.muEast * g.areaEast / g.dy);
    emit(i.center, i.south, -m.muSouth * AreaNS_divDy);
    emit(i.center, i.crossSW, -0.5 * m.muWest * g.areaWest / g.dy);
    emit(i.center, i.crossSE, 0.5 * m.muEast * g.areaEast / g.dy);
    emit(i.center, i.center, g.volume * m.rho / dt
        + (m.muSouth + m.muNorth) * AreaNS_divDy
        + 0.5 * (m.muWest * AreaW_divDx + m.muEast * AreaE_divDx));
}

inline double radialRhs(const Row& row, const Values& u, const Values& pre,
    double forceMaster, double forceSlave, double dt,bool includeAdvection=true)
{
    const auto& g = row.metrics;
    const auto& m = row.materials;
    const double qWest = m.rhoWest * 0.5 * (u.center + u.west) * g.areaWest;
    const double qWestPre = m.rhoWest * 0.5 * (pre.center + pre.west) * g.areaWest;
    const double westUpwind = qWest > 0.0 ? u.west : u.center;
    const double westUpwindPre = qWestPre > 0.0 ? pre.west : pre.center;
    const double diffWest = 2.0 * m.muWest * (u.center - u.west) / g.dx;

    const double qEast = m.rhoEast * 0.5 * (u.center + u.east) * g.areaEast;
    const double qEastPre = m.rhoEast * 0.5 * (pre.center + pre.east) * g.areaEast;
    const double eastUpwind = qEast > 0.0 ? u.center : u.east;
    const double eastUpwindPre = qEastPre > 0.0 ? pre.center : pre.east;
    const double diffEast = 2.0 * m.muEast * (u.east - u.center) / g.dx;

    const double vSouth = 0.5 * (u.crossSE + u.crossSW);
    const double vSouthPre = 0.5 * (pre.crossSE + pre.crossSW);
    const double qSouth = m.rhoSouth * vSouth * g.areaNS;
    const double qSouthPre = m.rhoSouth * vSouthPre * g.areaNS;
    const double diffSouth = m.muSouth * ((u.center - u.south) / g.shearLowDistance
        + (u.crossSE - u.crossSW) / g.dx);
    const double southUpwind = qSouth > 0.0 ? u.south : u.center;
    const double southUpwindPre = qSouthPre > 0.0 ? pre.south : pre.center;

    const double vNorth = 0.5 * (u.crossNE + u.crossNW);
    const double vNorthPre = 0.5 * (pre.crossNE + pre.crossNW);
    const double qNorth = m.rhoNorth * vNorth * g.areaNS;
    const double qNorthPre = m.rhoNorth * vNorthPre * g.areaNS;
    const double diffNorth = m.muNorth * ((u.north - u.center) / g.shearHighDistance
        + (u.crossNE - u.crossNW) / g.dx);
    const double northUpwind = qNorth > 0.0 ? u.center : u.north;
    const double northUpwindPre = qNorthPre > 0.0 ? pre.center : pre.north;

    // Preserve E-W+N-S and the original CN/AB2 grouping verbatim.
    const double n = qEast * eastUpwind - qWest * westUpwind
        + qNorth * northUpwind - qSouth * southUpwind;
    const double nPre = qEastPre * eastUpwindPre - qWestPre * westUpwindPre
        + qNorthPre * northUpwindPre - qSouthPre * southUpwindPre;
    const double nDiff = (diffEast * g.areaEast - diffWest * g.areaWest)
        + (diffNorth - diffSouth) * g.areaNS;
    if(!includeAdvection) return g.volume*m.rho/dt*u.center+.5*nDiff+g.volume*.5*(forceMaster+forceSlave);
    return g.volume * m.rho / dt * u.center - 1.5 * n + 0.5 * (nPre + nDiff)
        + g.volume * 0.5 * (forceMaster + forceSlave);
}

inline double axialRhs(const Row& row, const Values& u, const Values& pre,
    double forceMaster, double forceSlave, double dt,bool includeAdvection=true)
{
    const auto& g = row.metrics;
    const auto& m = row.materials;
    const double qSouth = m.rhoSouth * 0.5 * (u.center + u.south) * g.areaNS;
    const double qSouthPre = m.rhoSouth * 0.5 * (pre.center + pre.south) * g.areaNS;
    const double southUpwind = qSouth > 0.0 ? u.south : u.center;
    const double southUpwindPre = qSouthPre > 0.0 ? pre.south : pre.center;
    const double diffSouth = 2.0 * m.muSouth * (u.center - u.south) / g.dy;

    const double qNorth = m.rhoNorth * 0.5 * (u.center + u.north) * g.areaNS;
    const double qNorthPre = m.rhoNorth * 0.5 * (pre.center + pre.north) * g.areaNS;
    const double northUpwind = qNorth > 0.0 ? u.center : u.north;
    const double northUpwindPre = qNorthPre > 0.0 ? pre.center : pre.north;
    const double diffNorth = 2.0 * m.muNorth * (u.north - u.center) / g.dy;

    const double uWest = 0.5 * (u.crossNW + u.crossSW);
    const double uWestPre = 0.5 * (pre.crossNW + pre.crossSW);
    const double qWest = m.rhoWest * uWest * g.areaWest;
    const double qWestPre = m.rhoWest * uWestPre * g.areaWest;
    const double diffWest = m.muWest * ((u.center - u.west) / g.shearLowDistance
        + (u.crossNW - u.crossSW) / g.dy);
    const double westUpwind = qWest > 0.0 ? u.west : u.center;
    const double westUpwindPre = qWestPre > 0.0 ? pre.west : pre.center;

    const double uEast = 0.5 * (u.crossNE + u.crossSE);
    const double uEastPre = 0.5 * (pre.crossNE + pre.crossSE);
    const double qEast = m.rhoEast * uEast * g.areaEast;
    const double qEastPre = m.rhoEast * uEastPre * g.areaEast;
    const double diffEast = m.muEast * ((u.east - u.center) / g.shearHighDistance
        + (u.crossNE - u.crossSE) / g.dy);
    const double eastUpwind = qEast > 0.0 ? u.center : u.east;
    const double eastUpwindPre = qEastPre > 0.0 ? pre.center : pre.east;

    const double n = qEast * eastUpwind - qWest * westUpwind
        + qNorth * northUpwind - qSouth * southUpwind;
    const double nPre = qEastPre * eastUpwindPre - qWestPre * westUpwindPre
        + qNorthPre * northUpwindPre - qSouthPre * southUpwindPre;
    const double nDiff = (diffEast * g.areaEast - diffWest * g.areaWest)
        + (diffNorth - diffSouth) * g.areaNS;
    if(!includeAdvection) return g.volume*m.rho/dt*u.center+.5*nDiff+g.volume*.5*(forceMaster+forceSlave);
    return g.volume * m.rho / dt * u.center - 1.5 * n + 0.5 * (nPre + nDiff)
        + g.volume * 0.5 * (forceMaster + forceSlave);
}
}
