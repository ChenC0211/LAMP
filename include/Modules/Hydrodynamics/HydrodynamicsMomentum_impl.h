#pragma once

#include "Modules/Hydrodynamics/Detail/HydroNativeRZRowPreparation.h"

namespace HydroNativeMomentumDetail
{
enum class FaceRefinement { SameLevel, MasterCoarse, SlaveCoarse };

// These small values are local to one row. Preserve reciprocal XY ratios and
// RZ multiply-then-divide cross stresses; do not precompute cross scale factors.
struct CartesianRadial
{
    double ratio;
    double shear() const { return 1.0 / ratio; }
    double east() const { return ratio; }
    double west() const { return ratio; }
    double cross(double value) const { return value; }
    double normalViscosity(double westMu, double eastMu) const { return (westMu + eastMu) * ratio; }
};
struct CartesianAxial
{
    double ratio;
    double shearWest() const { return 1.0 / ratio; }
    double shearEast() const { return 1.0 / ratio; }
    double normal() const { return ratio; }
    double crossWest(double value) const { return value; }
    double crossEast(double value) const { return value; }
};
struct AxisymmetricRadial
{
    double AreaNS, AreaNS_divDy, AreaE_divDx, AreaW_divDx, dx;
    double shear() const { return AreaNS_divDy; }
    double east() const { return AreaE_divDx; }
    double west() const { return AreaW_divDx; }
    double cross(double value) const { return value * AreaNS / dx; }
    double normalViscosity(double westMu, double eastMu) const
    { return westMu * AreaW_divDx + eastMu * AreaE_divDx; }
};
struct AxisymmetricAxial
{
    double AreaE, AreaW, AreaNS_divDy, AreaE_divDx, AreaW_divDx, dy;
    double shearWest() const { return AreaW_divDx; }
    double shearEast() const { return AreaE_divDx; }
    double normal() const { return AreaNS_divDy; }
    double crossWest(double value) const { return value * AreaW / dy; }
    double crossEast(double value) const { return value * AreaE / dy; }
};

template<bool Axisymmetric>
double volume(const Edge& edge, const Cell& cell)
{
    if constexpr (Axisymmetric) return 2.0 * Pi * edge.midpointCoord().x * cell.area();
    else return cell.area();
}

struct RZFaceAreas { double AreaNS, AreaE, AreaW; };

// Matrix coefficients and RHS fluxes use the same native face areas. The
// coarse-side radii are the original virtual-child radii, not cell averages.
template<FaceRefinement Refinement>
RZFaceAreas radialAreas(const Edge& edge, const Cell& master, const Cell& slave)
{
    const Cell& cell = Refinement == FaceRefinement::MasterCoarse ? slave : master;
    double AreaNS = 2.0 * Pi * edge.midpointCoord().x * cell.size().x;
    double AreaE, AreaW;
    if constexpr (Refinement == FaceRefinement::MasterCoarse)
        AreaE = 2.0 * Pi * (edge.midpointCoord().x + 0.5 * cell.size().x) * cell.size().y;
    else AreaE = 2.0 * Pi * master.centerCoord().x * cell.size().y;
    if constexpr (Refinement == FaceRefinement::SlaveCoarse)
        AreaW = 2.0 * Pi * (edge.midpointCoord().x - 0.5 * cell.size().x) * cell.size().y;
    else AreaW = 2.0 * Pi * slave.centerCoord().x * cell.size().y;
    return {AreaNS, AreaE, AreaW};
}

inline RZFaceAreas axialAreas(const Edge& edge, const Cell& cell, const Vert& start, const Vert& end)
{
    double AreaNS = 2.0 * Pi * edge.midpointCoord().x * cell.size().x;
    double AreaE = 2.0 * Pi * end.coord().x * cell.size().y;
    double AreaW = 2.0 * Pi * start.coord().x * cell.size().y;
    return {AreaNS, AreaE, AreaW};
}


template<bool Axisymmetric, FaceRefinement Refinement>
auto radialMetrics(const Edge& edge, const Cell& master, const Cell& slave)
{
    const Cell& cell = Refinement == FaceRefinement::MasterCoarse ? slave : master;
    if constexpr (Axisymmetric)
    {
        const auto areas = radialAreas<Refinement>(edge, master, slave);
        double AreaNS = areas.AreaNS;
        double AreaE = areas.AreaE;
        double AreaW = areas.AreaW;
        double AreaNS_divDy = AreaNS / cell.size().y;
        double AreaE_divDx = AreaE / cell.size().x;
        double AreaW_divDx = AreaW / cell.size().x;
        return AxisymmetricRadial{AreaNS, AreaNS_divDy, AreaE_divDx, AreaW_divDx, cell.size().x};
    }
    else return CartesianRadial{cell.size().y / cell.size().x};
}

template<bool Axisymmetric>
auto axialMetrics(const Edge& edge, const Cell& cell, const Vert& start, const Vert& end)
{
    if constexpr (Axisymmetric)
    {
        const auto areas = axialAreas(edge, cell, start, end);
        double AreaNS = areas.AreaNS;
        double AreaE = areas.AreaE;
        double AreaW = areas.AreaW;
        double AreaNS_divDy = AreaNS / cell.size().y;
        double AreaE_divDx = AreaE / cell.size().x;
        double AreaW_divDx = AreaW / cell.size().x;
        return AxisymmetricAxial{AreaE, AreaW, AreaNS_divDy, AreaE_divDx, AreaW_divDx, cell.size().y};
    }
    else return CartesianAxial{cell.size().x / cell.size().y};
}

// Keep the XY difference-before-area and RZ separate radial products. These
// values live only while assembling a row; no mesh field/cache is retained.
struct CartesianFlux
{
    double dx, dy;
    double westArea() const { return dy; }
    double eastArea() const { return dy; }
    double axialArea() const { return dx; }
    double diffusionNet(double east, double west, double north, double south) const
    { return (east - west) * dy + (north - south) * dx; }
};
struct AxisymmetricFlux
{
    RZFaceAreas areas;
    double dx, dy;
    double westArea() const { return areas.AreaW; }
    double eastArea() const { return areas.AreaE; }
    double axialArea() const { return areas.AreaNS; }
    double diffusionNet(double east, double west, double north, double south) const
    { return (east * areas.AreaE - west * areas.AreaW) + (north - south) * areas.AreaNS; }
};

template<bool Axisymmetric, FaceRefinement Refinement>
auto radialFluxMetrics(const Edge& edge, const Cell& master, const Cell& slave)
{
    const Cell& cell = Refinement == FaceRefinement::MasterCoarse ? slave : master;
    if constexpr (Axisymmetric)
        return AxisymmetricFlux{radialAreas<Refinement>(edge, master, slave), cell.size().x, cell.size().y};
    else return CartesianFlux{cell.size().x, cell.size().y};
}

template<bool Axisymmetric>
auto axialFluxMetrics(const Edge& edge, const Cell& cell, const Vert& start, const Vert& end)
{
    if constexpr (Axisymmetric)
        return AxisymmetricFlux{axialAreas(edge, cell, start, end), cell.size().x, cell.size().y};
    else return CartesianFlux{cell.size().x, cell.size().y};
}

}


// Shared native matrix topology and interpolation. Coordinate-specific row
// arithmetic, RZ regular dispatch and leaf-side policies retain their order.
template<bool Axisymmetric, typename Mesh>
void Hydrodynamics::buildNativeMomentumMatrix(Mesh& mesh, const double& dt, HydroMomentumAssembly& assembly,
    const std::vector<unsigned char>* selectedRows)
{
    if constexpr (Axisymmetric)
    {
        if (selectedRows && selectedRows->size()!=static_cast<std::size_t>(mesh.vEdgeIndexNum()+mesh.hEdgeIndexNum()))
            throw std::invalid_argument("native RZ matrix row selection has wrong size");
    }
    assembly.entries.clear();
    vector<vector<Eigen::Triplet<double>>> momCoeff(mesh.vEdgeIndexNum() + mesh.hEdgeIndexNum());
#pragma omp parallel for num_threads(thread_num)
    for (int i = 0; i < mesh.vEdgeIndexNum(); i++)
    {
        Edge& vEdge = mesh.vEdges(mesh.vEdgeIndexToLabel[i]);
        int row = i;
        if constexpr (Axisymmetric)
        {
            if (selectedRows && !(*selectedRows)[row]) continue;
        }
        if (vEdge.masterCell != -1 && vEdge.slaveCell != -1)
        {
            Cell& masterCell = mesh.cells(vEdge.masterCell);
            Cell& slaveCell = mesh.cells(vEdge.slaveCell);
            Vert& startVert = mesh.verts(vEdge.edgeToVert[0]);
            Vert& endVert = mesh.verts(vEdge.edgeToVert[1]);
            if (masterCell.level() == slaveCell.level())
            {
                if constexpr (Axisymmetric)
                {
                    HydroNativeRZRow::Row nativeRow;
                    if (HydroNativeRZRowDetail::regularRow(mesh, vEdge, true, rho, mu, vertMu, nativeRow))
                    {
                        HydroNativeRZRow::appendRadialMatrix(nativeRow, dt,
                            [&](int nativeIndex, int column, double coefficient) {
                                momCoeff[row].push_back(Trip(nativeIndex, column, coefficient));
                            });
                        continue;
                    }
                }
                double _dA = HydroNativeMomentumDetail::volume<Axisymmetric>(vEdge, masterCell);
                double _rho = 0.5 * (rho.value[masterCell.index()] + rho.value[slaveCell.index()]);
                double _muWest = mu.value[slaveCell.index()];
                double _muEast = mu.value[masterCell.index()];
                double _muSouth = vertMu.value[startVert.index()];
                double _muNorth = vertMu.value[endVert.index()];
                const auto metrics = HydroNativeMomentumDetail::radialMetrics<Axisymmetric, HydroNativeMomentumDetail::FaceRefinement::SameLevel>(vEdge, masterCell, slaveCell);
                double _dxDivdySouth = 0.0;
                double _dxDivdyNorth = 0.0;
                if (startVert.vertToVEdge[0] != -1)
                {
                    Edge& southVEdge = mesh.vEdges(startVert.vertToVEdge[0]);
                    _dxDivdySouth = metrics.shear();
                    if (mesh.cells(southVEdge.masterCell).level() == masterCell.level() && mesh.cells(southVEdge.slaveCell).level() == slaveCell.level())
                    {
                        momCoeff[row].push_back(Trip(row, southVEdge.index(), -0.5 * _muSouth * _dxDivdySouth));
                    }
                    else if (mesh.cells(southVEdge.masterCell).level() < masterCell.level())
                    {
                        int _supSonCellIndex = mesh.cells(southVEdge.masterCell).supCellIndex * 4 + 3;
                        for (int j = 0; j < supSonCells_[_supSonCellIndex].vEdges[0].ind.size(); j++)
                        {
                            momCoeff[row].push_back(Trip(row, supSonCells_[_supSonCellIndex].vEdges[0].ind[j], -0.5 * _muSouth * _dxDivdySouth * supSonCells_[_supSonCellIndex].vEdges[0].coe[j]));
                        }
                    }
                    else if (mesh.cells(southVEdge.slaveCell).level() < slaveCell.level())
                    {
                        int _supSonCellIndex = mesh.cells(southVEdge.slaveCell).supCellIndex * 4 + 2;
                        for (int j = 0; j < supSonCells_[_supSonCellIndex].vEdges[1].ind.size(); j++)
                        {
                            momCoeff[row].push_back(Trip(row, supSonCells_[_supSonCellIndex].vEdges[1].ind[j], -0.5 * _muSouth * _dxDivdySouth * supSonCells_[_supSonCellIndex].vEdges[1].coe[j]));
                        }
                    }
                    else if (mesh.cells(southVEdge.masterCell).level() > masterCell.level() && mesh.cells(southVEdge.slaveCell).level() > slaveCell.level())
                    {
                        Cell& neiParentCell = mesh.cells(mesh.cells(southVEdge.masterCell).parent);
                        Cell& sonCell0 = mesh.cells(neiParentCell.son[0]);
                        Cell& sonCell3 = mesh.cells(neiParentCell.son[3]);
                        if constexpr (Axisymmetric)
                        {
                            const auto emit = [&](int column, double coefficient) {
                                momCoeff[row].emplace_back(row, column, coefficient);
                            };
                            HydroNativeRZRowDetail::visitLeafLowSide(mesh, sonCell0, true, -0.5 * 0.5 * _muSouth * _dxDivdySouth, emit);
                            HydroNativeRZRowDetail::visitLeafLowSide(mesh, sonCell3, true, -0.5 * 0.5 * _muSouth * _dxDivdySouth, emit);
                        }
                        else
                        {
                            momCoeff[row].push_back(Trip(row, mesh.vEdges(sonCell0.cellToVEdge[0][0]).index(), -0.5 * 0.5 * _muSouth * _dxDivdySouth));
                            momCoeff[row].push_back(Trip(row, mesh.vEdges(sonCell3.cellToVEdge[0][0]).index(), -0.5 * 0.5 * _muSouth * _dxDivdySouth));
                        }
                    }
                    else if (mesh.cells(southVEdge.masterCell).level() > masterCell.level())
                    {
                        int _supCellIndex = mesh.cells(southVEdge.slaveCell).supCellIndex;
                        for (int j = 0; j < supCells_[_supCellIndex].vEdges[1].ind.size(); j++)
                        {
                            momCoeff[row].push_back(Trip(row, supCells_[_supCellIndex].vEdges[1].ind[j], -0.5 * _muSouth * _dxDivdySouth * supCells_[_supCellIndex].vEdges[1].coe[j]));
                        }
                    }
                    else if (mesh.cells(southVEdge.slaveCell).level() > slaveCell.level())
                    {
                        int _supCellIndex = mesh.cells(southVEdge.masterCell).supCellIndex;
                        for (int j = 0; j < supCells_[_supCellIndex].vEdges[0].ind.size(); j++)
                        {
                            momCoeff[row].push_back(Trip(row, supCells_[_supCellIndex].vEdges[0].ind[j], -0.5 * _muSouth * _dxDivdySouth * supCells_[_supCellIndex].vEdges[0].coe[j]));
                        }
                    }
                }
                else if (startVert.supCell != -1)
                {
                    _dxDivdySouth = metrics.shear();
                    Cell& _supCell = mesh.cells(startVert.supCell);
                    int _supSonCellIndex = _supCell.supCellIndex * 4 + 3;
                    for (int j = 0; j < supSonCells_[_supSonCellIndex].vEdges[1].ind.size(); j++)
                    {
                        momCoeff[row].push_back(Trip(row, supSonCells_[_supSonCellIndex].vEdges[1].ind[j], -0.5 * _muSouth * _dxDivdySouth * supSonCells_[_supSonCellIndex].vEdges[1].coe[j]));
                    }
                }
                if (endVert.vertToVEdge[1] != -1)
                {
                    Edge& northVEdge = mesh.vEdges(endVert.vertToVEdge[1]);
                    _dxDivdyNorth = metrics.shear();
                    if (mesh.cells(northVEdge.masterCell).level() == masterCell.level() && mesh.cells(northVEdge.slaveCell).level() == slaveCell.level())
                    {
                        momCoeff[row].push_back(Trip(row, northVEdge.index(), -0.5 * _muNorth * _dxDivdyNorth));
                    }
                    else if (mesh.cells(northVEdge.masterCell).level() < masterCell.level())
                    {
                        int _supSonCellIndex = mesh.cells(northVEdge.masterCell).supCellIndex * 4;
                        for (int j = 0; j < supSonCells_[_supSonCellIndex].vEdges[0].ind.size(); j++)
                        {
                            momCoeff[row].push_back(Trip(row, supSonCells_[_supSonCellIndex].vEdges[0].ind[j], -0.5 * _muNorth * _dxDivdyNorth * supSonCells_[_supSonCellIndex].vEdges[0].coe[j]));
                        }
                    }
                    else if (mesh.cells(northVEdge.slaveCell).level() < slaveCell.level())
                    {
                        int _supSonCellIndex = mesh.cells(northVEdge.slaveCell).supCellIndex * 4 + 1;
                        for (int j = 0; j < supSonCells_[_supSonCellIndex].vEdges[1].ind.size(); j++)
                        {
                            momCoeff[row].push_back(Trip(row, supSonCells_[_supSonCellIndex].vEdges[1].ind[j], -0.5 * _muNorth * _dxDivdyNorth * supSonCells_[_supSonCellIndex].vEdges[1].coe[j]));
                        }
                    }
                    else if (mesh.cells(northVEdge.masterCell).level() > masterCell.level() && mesh.cells(northVEdge.slaveCell).level() > slaveCell.level())
                    {
                        Cell& neiParentCell = mesh.cells(mesh.cells(northVEdge.masterCell).parent);
                        Cell& sonCell0 = mesh.cells(neiParentCell.son[0]);
                        Cell& sonCell3 = mesh.cells(neiParentCell.son[3]);
                        if constexpr (Axisymmetric)
                        {
                            const auto emit = [&](int column, double coefficient) {
                                momCoeff[row].emplace_back(row, column, coefficient);
                            };
                            HydroNativeRZRowDetail::visitLeafLowSide(mesh, sonCell0, true, -0.5 * 0.5 * _muNorth * _dxDivdyNorth, emit);
                            HydroNativeRZRowDetail::visitLeafLowSide(mesh, sonCell3, true, -0.5 * 0.5 * _muNorth * _dxDivdyNorth, emit);
                        }
                        else
                        {
                            momCoeff[row].push_back(Trip(row, mesh.vEdges(sonCell0.cellToVEdge[0][0]).index(), -0.5 * 0.5 * _muNorth * _dxDivdyNorth));
                            momCoeff[row].push_back(Trip(row, mesh.vEdges(sonCell3.cellToVEdge[0][0]).index(), -0.5 * 0.5 * _muNorth * _dxDivdyNorth));
                        }
                    }
                    else if (mesh.cells(northVEdge.masterCell).level() > masterCell.level())
                    {
                        int _supCellIndex = mesh.cells(northVEdge.slaveCell).supCellIndex;
                        for (int j = 0; j < supCells_[_supCellIndex].vEdges[1].ind.size(); j++)
                        {
                            momCoeff[row].push_back(Trip(row, supCells_[_supCellIndex].vEdges[1].ind[j], -0.5 * _muNorth * _dxDivdyNorth * supCells_[_supCellIndex].vEdges[1].coe[j]));
                        }
                    }
                    else if (mesh.cells(northVEdge.slaveCell).level() > slaveCell.level())
                    {
                        int _supCellIndex = mesh.cells(northVEdge.masterCell).supCellIndex;
                        for (int j = 0; j < supCells_[_supCellIndex].vEdges[0].ind.size(); j++)
                        {
                            momCoeff[row].push_back(Trip(row, supCells_[_supCellIndex].vEdges[0].ind[j], -0.5 * _muNorth * _dxDivdyNorth * supCells_[_supCellIndex].vEdges[0].coe[j]));
                        }
                    }
                }
                else if (endVert.supCell != -1)
                {
                    _dxDivdyNorth = metrics.shear();
                    Cell& _supCell = mesh.cells(endVert.supCell);
                    int _supSonCellIndex = _supCell.supCellIndex * 4;
                    for (int j = 0; j < supSonCells_[_supSonCellIndex].vEdges[1].ind.size(); j++)
                    {
                        momCoeff[row].push_back(Trip(row, supSonCells_[_supSonCellIndex].vEdges[1].ind[j], -0.5 * _muNorth * _dxDivdyNorth * supSonCells_[_supSonCellIndex].vEdges[1].coe[j]));
                    }
                }
                if (masterCell.supCellIndex == -1)
                {
                    momCoeff[row].push_back(Trip(row, mesh.vEdges(masterCell.cellToVEdge[1][0]).index(), -_muEast * metrics.east()));
                    momCoeff[row].push_back(Trip(row, mesh.vEdgeIndexNum() + mesh.hEdges(masterCell.cellToHEdge[0][0]).index(), metrics.cross(0.5 * _muSouth)));
                    momCoeff[row].push_back(Trip(row, mesh.vEdgeIndexNum() + mesh.hEdges(masterCell.cellToHEdge[1][0]).index(), metrics.cross(-0.5 * _muNorth)));
                }
                else
                {
                    int _supCellIndex = masterCell.supCellIndex;
                    for (int j = 0; j < supCells_[_supCellIndex].vEdges[1].ind.size(); j++)
                    {
                        momCoeff[row].push_back(Trip(row, supCells_[_supCellIndex].vEdges[1].ind[j], -_muEast * metrics.east() * supCells_[_supCellIndex].vEdges[1].coe[j]));
                    }
                    for (int j = 0; j < supCells_[_supCellIndex].hEdges[0].ind.size(); j++)
                    {
                        momCoeff[row].push_back(Trip(row, mesh.vEdgeIndexNum() + supCells_[_supCellIndex].hEdges[0].ind[j], metrics.cross(0.5 * _muSouth) * supCells_[_supCellIndex].hEdges[0].coe[j]));
                    }
                    for (int j = 0; j < supCells_[_supCellIndex].hEdges[1].ind.size(); j++)
                    {
                        momCoeff[row].push_back(Trip(row, mesh.vEdgeIndexNum() + supCells_[_supCellIndex].hEdges[1].ind[j], metrics.cross(-0.5 * _muNorth) * supCells_[_supCellIndex].hEdges[1].coe[j]));
                    }
                }
                if (slaveCell.supCellIndex == -1)
                {
                    momCoeff[row].push_back(Trip(row, mesh.vEdges(slaveCell.cellToVEdge[0][0]).index(), -_muWest * metrics.west()));
                    momCoeff[row].push_back(Trip(row, mesh.vEdgeIndexNum() + mesh.hEdges(slaveCell.cellToHEdge[0][0]).index(), metrics.cross(-0.5 * _muSouth)));
                    momCoeff[row].push_back(Trip(row, mesh.vEdgeIndexNum() + mesh.hEdges(slaveCell.cellToHEdge[1][0]).index(), metrics.cross(0.5 * _muNorth)));
                }
                else
                {
                    int _supCellIndex = slaveCell.supCellIndex;
                    for (int j = 0; j < supCells_[_supCellIndex].vEdges[0].ind.size(); j++)
                    {
                        momCoeff[row].push_back(Trip(row, supCells_[_supCellIndex].vEdges[0].ind[j], -_muWest * metrics.west() * supCells_[_supCellIndex].vEdges[0].coe[j]));
                    }
                    for (int j = 0; j < supCells_[_supCellIndex].hEdges[0].ind.size(); j++)
                    {
                        momCoeff[row].push_back(Trip(row, mesh.vEdgeIndexNum() + supCells_[_supCellIndex].hEdges[0].ind[j], metrics.cross(-0.5 * _muSouth) * supCells_[_supCellIndex].hEdges[0].coe[j]));
                    }
                    for (int j = 0; j < supCells_[_supCellIndex].hEdges[1].ind.size(); j++)
                    {
                        momCoeff[row].push_back(Trip(row, mesh.vEdgeIndexNum() + supCells_[_supCellIndex].hEdges[1].ind[j], metrics.cross(0.5 * _muNorth) * supCells_[_supCellIndex].hEdges[1].coe[j]));
                    }
                }
                momCoeff[row].push_back(Trip(row, vEdge.index(), _dA * _rho / dt + metrics.normalViscosity(_muWest, _muEast) + 0.5 * (_muSouth * _dxDivdySouth + _muNorth * _dxDivdyNorth)));
            }
            else if (masterCell.level() < slaveCell.level())
            {
                double _dA = HydroNativeMomentumDetail::volume<Axisymmetric>(vEdge, slaveCell);
                double _rhoMaster = 0.0;
                int supSonCellIndex = masterCell.centerCoord().y < slaveCell.centerCoord().y? masterCell.supCellIndex * 4 + 3 : masterCell.supCellIndex * 4;
                for (int j = 0; j < supSonCells_[supSonCellIndex].center.ind.size(); j++)
                {
                    _rhoMaster += rho.value[supSonCells_[supSonCellIndex].center.ind[j]] * supSonCells_[supSonCellIndex].center.coe[j];
                }
                double _rho = 0.5 * (rho.value[slaveCell.index()] + _rhoMaster);
                double _muWest = mu.value[slaveCell.index()];
                double _muEast = 0.0;
                for (int j = 0; j < supSonCells_[supSonCellIndex].center.ind.size(); j++)
                {
                    _muEast += mu.value[supSonCells_[supSonCellIndex].center.ind[j]] * supSonCells_[supSonCellIndex].center.coe[j];
                }
                double _muSouth = vertMu.value[startVert.index()];
                double _muNorth = vertMu.value[endVert.index()];
                const auto metrics = HydroNativeMomentumDetail::radialMetrics<Axisymmetric, HydroNativeMomentumDetail::FaceRefinement::MasterCoarse>(vEdge, masterCell, slaveCell);
                double _dxDivdySouth = 0.0;
                double _dxDivdyNorth = 0.0;
                if (startVert.vertToVEdge[0] != -1)
                {
                    Edge& southVEdge = mesh.vEdges(startVert.vertToVEdge[0]);
                    _dxDivdySouth = metrics.shear();
                    if (mesh.cells(southVEdge.slaveCell).level() == slaveCell.level())
                    {
                        momCoeff[row].push_back(Trip(row, southVEdge.index(), -0.5 * _muSouth * _dxDivdySouth));
                    }
                    else if (mesh.cells(southVEdge.slaveCell).level() < slaveCell.level())
                    {
                        int _supSonCellIndex = mesh.cells(southVEdge.slaveCell).supCellIndex * 4 + 2;
                        for (int j = 0; j < supSonCells_[_supSonCellIndex].vEdges[1].ind.size(); j++)
                        {
                            momCoeff[row].push_back(Trip(row, supSonCells_[_supSonCellIndex].vEdges[1].ind[j], -0.5 * _muSouth * _dxDivdySouth * supSonCells_[_supSonCellIndex].vEdges[1].coe[j]));
                        }
                    }
                }
                if (endVert.vertToVEdge[1] != -1)
                {
                    Edge& northVEdge = mesh.vEdges(endVert.vertToVEdge[1]);
                    _dxDivdyNorth = metrics.shear();
                    if (mesh.cells(northVEdge.slaveCell).level() == slaveCell.level())
                    {
                        momCoeff[row].push_back(Trip(row, northVEdge.index(), -0.5 * _muNorth * _dxDivdyNorth));
                    }
                    else if (mesh.cells(northVEdge.slaveCell).level() < slaveCell.level())
                    {
                        int _supSonCellIndex = mesh.cells(northVEdge.slaveCell).supCellIndex * 4 + 1;
                        for (int j = 0; j < supSonCells_[_supSonCellIndex].vEdges[1].ind.size(); j++)
                        {
                            momCoeff[row].push_back(Trip(row, supSonCells_[_supSonCellIndex].vEdges[1].ind[j], -0.5 * _muNorth * _dxDivdyNorth * supSonCells_[_supSonCellIndex].vEdges[1].coe[j]));
                        }
                    }
                }
                {
                    for (int j = 0; j < supSonCells_[supSonCellIndex].vEdges[1].ind.size(); j++)
                    {
                        momCoeff[row].push_back(Trip(row, supSonCells_[supSonCellIndex].vEdges[1].ind[j], -_muEast * metrics.east() * supSonCells_[supSonCellIndex].vEdges[1].coe[j]));
                    }
                    for (int j = 0; j < supSonCells_[supSonCellIndex].hEdges[0].ind.size(); j++)
                    {
                        momCoeff[row].push_back(Trip(row, mesh.vEdgeIndexNum() + supSonCells_[supSonCellIndex].hEdges[0].ind[j], metrics.cross(0.5 * _muSouth) * supSonCells_[supSonCellIndex].hEdges[0].coe[j]));
                    }
                    for (int j = 0; j < supSonCells_[supSonCellIndex].hEdges[1].ind.size(); j++)
                    {
                        momCoeff[row].push_back(Trip(row, mesh.vEdgeIndexNum() + supSonCells_[supSonCellIndex].hEdges[1].ind[j], metrics.cross(-0.5 * _muNorth) * supSonCells_[supSonCellIndex].hEdges[1].coe[j]));
                    }
                }
                {
                    momCoeff[row].push_back(Trip(row, mesh.vEdges(slaveCell.cellToVEdge[0][0]).index(), -_muWest * metrics.west()));
                    momCoeff[row].push_back(Trip(row, mesh.vEdgeIndexNum() + mesh.hEdges(slaveCell.cellToHEdge[0][0]).index(), metrics.cross(-0.5 * _muSouth)));
                    momCoeff[row].push_back(Trip(row, mesh.vEdgeIndexNum() + mesh.hEdges(slaveCell.cellToHEdge[1][0]).index(), metrics.cross(0.5 * _muNorth)));
                }
                momCoeff[row].push_back(Trip(row, vEdge.index(), _dA * _rho / dt + metrics.normalViscosity(_muWest, _muEast) + 0.5 * (_muSouth * _dxDivdySouth + _muNorth * _dxDivdyNorth)));
            }
            else
            {
                double _dA = HydroNativeMomentumDetail::volume<Axisymmetric>(vEdge, masterCell);
                double _rhoSlave = 0.0;
                int supSonCellIndex = slaveCell.centerCoord().y < masterCell.centerCoord().y? slaveCell.supCellIndex * 4 + 2 : slaveCell.supCellIndex * 4 + 1;
                for (int j = 0; j < supSonCells_[supSonCellIndex].center.ind.size(); j++)
                {
                    _rhoSlave += rho.value[supSonCells_[supSonCellIndex].center.ind[j]] * supSonCells_[supSonCellIndex].center.coe[j];
                }
                double _rho = 0.5 * (rho.value[masterCell.index()] + _rhoSlave);
                double _muWest = 0.0;
                for (int j = 0; j < supSonCells_[supSonCellIndex].center.ind.size(); j++)
                {
                    _muWest += mu.value[supSonCells_[supSonCellIndex].center.ind[j]] * supSonCells_[supSonCellIndex].center.coe[j];
                }
                double _muEast = mu.value[masterCell.index()];
                double _muSouth = vertMu.value[startVert.index()];
                double _muNorth = vertMu.value[endVert.index()];
                const auto metrics = HydroNativeMomentumDetail::radialMetrics<Axisymmetric, HydroNativeMomentumDetail::FaceRefinement::SlaveCoarse>(vEdge, masterCell, slaveCell);
                double _dxDivdySouth = 0.0;
                double _dxDivdyNorth = 0.0;
                if (startVert.vertToVEdge[0] != -1)
                {
                    Edge& southVEdge = mesh.vEdges(startVert.vertToVEdge[0]);
                    _dxDivdySouth = metrics.shear();
                    if (mesh.cells(southVEdge.masterCell).level() == masterCell.level())
                    {
                        momCoeff[row].push_back(Trip(row, southVEdge.index(), -0.5 * _muSouth * _dxDivdySouth));
                    }
                    else if (mesh.cells(southVEdge.masterCell).level() < masterCell.level())
                    {
                        int _supSonCellIndex = mesh.cells(southVEdge.masterCell).supCellIndex * 4 + 3;
                        for (int j = 0; j < supSonCells_[_supSonCellIndex].vEdges[0].ind.size(); j++)
                        {
                            momCoeff[row].push_back(Trip(row, supSonCells_[_supSonCellIndex].vEdges[0].ind[j], -0.5 * _muSouth * _dxDivdySouth * supSonCells_[_supSonCellIndex].vEdges[0].coe[j]));
                        }
                    }
                }
                if (endVert.vertToVEdge[1] != -1)
                {
                    Edge& northVEdge = mesh.vEdges(endVert.vertToVEdge[1]);
                    _dxDivdyNorth = metrics.shear();
                    if (mesh.cells(northVEdge.masterCell).level() == masterCell.level())
                    {
                        momCoeff[row].push_back(Trip(row, northVEdge.index(), -0.5 * _muNorth * _dxDivdyNorth));
                    }
                    else if (mesh.cells(northVEdge.masterCell).level() < masterCell.level())
                    {
                        int _supSonCellIndex = mesh.cells(northVEdge.masterCell).supCellIndex * 4;
                        for (int j = 0; j < supSonCells_[_supSonCellIndex].vEdges[0].ind.size(); j++)
                        {
                            momCoeff[row].push_back(Trip(row, supSonCells_[_supSonCellIndex].vEdges[0].ind[j], -0.5 * _muNorth * _dxDivdyNorth * supSonCells_[_supSonCellIndex].vEdges[0].coe[j]));
                        }
                    }
                }
                {
                    for (int j = 0; j < supSonCells_[supSonCellIndex].vEdges[0].ind.size(); j++)
                    {
                        momCoeff[row].push_back(Trip(row, supSonCells_[supSonCellIndex].vEdges[0].ind[j], -_muWest * metrics.west() * supSonCells_[supSonCellIndex].vEdges[0].coe[j]));
                    }
                    for (int j = 0; j < supSonCells_[supSonCellIndex].hEdges[0].ind.size(); j++)
                    {
                        momCoeff[row].push_back(Trip(row, mesh.vEdgeIndexNum() + supSonCells_[supSonCellIndex].hEdges[0].ind[j], metrics.cross(-0.5 * _muSouth) * supSonCells_[supSonCellIndex].hEdges[0].coe[j]));
                    }
                    for (int j = 0; j < supSonCells_[supSonCellIndex].hEdges[1].ind.size(); j++)
                    {
                        momCoeff[row].push_back(Trip(row, mesh.vEdgeIndexNum() + supSonCells_[supSonCellIndex].hEdges[1].ind[j], metrics.cross(0.5 * _muNorth) * supSonCells_[supSonCellIndex].hEdges[1].coe[j]));
                    }
                }
                {
                    momCoeff[row].push_back(Trip(row, mesh.vEdges(masterCell.cellToVEdge[1][0]).index(), -_muEast * metrics.east()));
                    momCoeff[row].push_back(Trip(row, mesh.vEdgeIndexNum() + mesh.hEdges(masterCell.cellToHEdge[0][0]).index(), metrics.cross(0.5 * _muSouth)));
                    momCoeff[row].push_back(Trip(row, mesh.vEdgeIndexNum() + mesh.hEdges(masterCell.cellToHEdge[1][0]).index(), metrics.cross(-0.5 * _muNorth)));
                }
                momCoeff[row].push_back(Trip(row, vEdge.index(), _dA * _rho / dt + metrics.normalViscosity(_muWest, _muEast) + 0.5 * (_muSouth * _dxDivdySouth + _muNorth * _dxDivdyNorth)));
            }
        }
        else
        {
            momCoeff[row].push_back(Trip(row, vEdge.index(), 1.0));
        }
    }
#pragma omp parallel for num_threads(thread_num)
for (int i = 0; i < mesh.hEdgeIndexNum(); i++)
    {
        Edge& hEdge = mesh.hEdges(mesh.hEdgeIndexToLabel[i]);
        int row = mesh.vEdgeIndexNum() + i;
        if constexpr (Axisymmetric)
        {
            if (selectedRows && !(*selectedRows)[row]) continue;
        }
        if (!activeOptions_.closedDomainProjection && fabs(hEdge.midpointCoord().y) < 1.0e-4)
        {
            momCoeff[row].push_back(Trip(row, hEdge.index() + mesh.vEdgeIndexNum(), 1.0));
            continue;
        }
        if (hEdge.masterCell != -1 && hEdge.slaveCell != -1)
        {
            Cell& masterCell = mesh.cells(hEdge.masterCell);
            Cell& slaveCell = mesh.cells(hEdge.slaveCell);
            Vert& startVert = mesh.verts(hEdge.edgeToVert[0]);
            Vert& endVert = mesh.verts(hEdge.edgeToVert[1]);
            if (masterCell.level() == slaveCell.level())
            {
                if constexpr (Axisymmetric)
                {
                    HydroNativeRZRow::Row nativeRow;
                    if (HydroNativeRZRowDetail::regularRow(mesh, hEdge, false, rho, mu, vertMu, nativeRow))
                    {
                        HydroNativeRZRow::appendAxialMatrix(nativeRow, dt,
                            [&](int nativeIndex, int column, double coefficient) {
                                momCoeff[row].push_back(Trip(nativeIndex, column, coefficient));
                            });
                        continue;
                    }
                }
                double _dA = HydroNativeMomentumDetail::volume<Axisymmetric>(hEdge, masterCell);
                double _rho = 0.5 * (rho.value[masterCell.index()] + rho.value[slaveCell.index()]);
                double _muSouth = mu.value[slaveCell.index()];
                double _muNorth = mu.value[masterCell.index()];
                double _muWest = vertMu.value[startVert.index()];
                double _muEast = vertMu.value[endVert.index()];
                const auto metrics = HydroNativeMomentumDetail::axialMetrics<Axisymmetric>(hEdge, masterCell, startVert, endVert);
                double _dyDivdxWest = 0.0;
                double _dyDivdxEast = 0.0;
                if (startVert.vertToHEdge[0] != -1)
                {
                    Edge& westHEdge = mesh.hEdges(startVert.vertToHEdge[0]);
                    _dyDivdxWest = metrics.shearWest();
                    if (mesh.cells(westHEdge.masterCell).level() == masterCell.level() && mesh.cells(westHEdge.slaveCell).level() == slaveCell.level())
                    {
                        momCoeff[row].push_back(Trip(row, mesh.vEdgeIndexNum() + westHEdge.index(), -0.5 * _muWest * _dyDivdxWest));
                    }
                    else if (mesh.cells(westHEdge.masterCell).level() < masterCell.level())
                    {
                        int _supSonCellIndex = mesh.cells(westHEdge.masterCell).supCellIndex * 4 + 1;
                        for (int j = 0; j < supSonCells_[_supSonCellIndex].hEdges[0].ind.size(); j++)
                        {
                            momCoeff[row].push_back(Trip(row, mesh.vEdgeIndexNum() + supSonCells_[_supSonCellIndex].hEdges[0].ind[j], -0.5 * _muWest * _dyDivdxWest * supSonCells_[_supSonCellIndex].hEdges[0].coe[j]));
                        }
                    }
                    else if (mesh.cells(westHEdge.slaveCell).level() < slaveCell.level())
                    {
                        int _supSonCellIndex = mesh.cells(westHEdge.slaveCell).supCellIndex * 4 + 2;
                        for (int j = 0; j < supSonCells_[_supSonCellIndex].hEdges[1].ind.size(); j++)
                        {
                            momCoeff[row].push_back(Trip(row, mesh.vEdgeIndexNum() + supSonCells_[_supSonCellIndex].hEdges[1].ind[j], -0.5 * _muWest * _dyDivdxWest * supSonCells_[_supSonCellIndex].hEdges[1].coe[j]));
                        }
                    }
                    else if (mesh.cells(westHEdge.masterCell).level() > masterCell.level() && mesh.cells(westHEdge.slaveCell).level() > slaveCell.level())
                    {
                        Cell& neiParentCell = mesh.cells(mesh.cells(westHEdge.masterCell).parent);
                        Cell& sonCell0 = mesh.cells(neiParentCell.son[0]);
                        Cell& sonCell1 = mesh.cells(neiParentCell.son[1]);
                        if constexpr (Axisymmetric)
                        {
                            const auto emit = [&](int column, double coefficient) {
                                momCoeff[row].emplace_back(row, column, coefficient);
                            };
                            HydroNativeRZRowDetail::visitLeafLowSide(mesh, sonCell0, false, -0.5 * 0.5 * _muWest * _dyDivdxWest, emit);
                            HydroNativeRZRowDetail::visitLeafLowSide(mesh, sonCell1, false, -0.5 * 0.5 * _muWest * _dyDivdxWest, emit);
                        }
                        else
                        {
                            momCoeff[row].push_back(Trip(row, mesh.vEdgeIndexNum() + mesh.hEdges(sonCell0.cellToHEdge[0][0]).index(), -0.5 * 0.5 * _muWest * _dyDivdxWest));
                            momCoeff[row].push_back(Trip(row, mesh.vEdgeIndexNum() + mesh.hEdges(sonCell1.cellToHEdge[0][0]).index(), -0.5 * 0.5 * _muWest * _dyDivdxWest));
                        }
                    }
                    else if (mesh.cells(westHEdge.masterCell).level() > masterCell.level())
                    {
                        int _supCellIndex = mesh.cells(westHEdge.slaveCell).supCellIndex;
                        for (int j = 0; j < supCells_[_supCellIndex].hEdges[1].ind.size(); j++)
                        {
                            momCoeff[row].push_back(Trip(row, mesh.vEdgeIndexNum() + supCells_[_supCellIndex].hEdges[1].ind[j], -0.5 * _muWest * _dyDivdxWest * supCells_[_supCellIndex].hEdges[1].coe[j]));
                        }
                    }
                    else if (mesh.cells(westHEdge.slaveCell).level() > slaveCell.level())
                    {
                        int _supCellIndex = mesh.cells(westHEdge.masterCell).supCellIndex;
                        for (int j = 0; j < supCells_[_supCellIndex].hEdges[0].ind.size(); j++)
                        {
                            momCoeff[row].push_back(Trip(row, mesh.vEdgeIndexNum() + supCells_[_supCellIndex].hEdges[0].ind[j], -0.5 * _muWest * _dyDivdxWest * supCells_[_supCellIndex].hEdges[0].coe[j]));
                        }
                    }
                }
                else if (startVert.supCell != -1)
                {
                    _dyDivdxWest = metrics.shearWest();
                    Cell& _supCell = mesh.cells(startVert.supCell);
                    int _supSonCellIndex = _supCell.supCellIndex * 4 + 2;
                    for (int j = 0; j < supSonCells_[_supSonCellIndex].hEdges[0].ind.size(); j++)
                    {
                        momCoeff[row].push_back(Trip(row, mesh.vEdgeIndexNum() + supSonCells_[_supSonCellIndex].hEdges[0].ind[j], -0.5 * _muWest * _dyDivdxWest * supSonCells_[_supSonCellIndex].hEdges[0].coe[j]));
                    }
                }
                if (endVert.vertToHEdge[1] != -1)
                {
                    Edge& eastHEdge = mesh.hEdges(endVert.vertToHEdge[1]);
                    _dyDivdxEast = metrics.shearEast();
                    if (mesh.cells(eastHEdge.masterCell).level() == masterCell.level() && mesh.cells(eastHEdge.slaveCell).level() == slaveCell.level())
                    {
                        momCoeff[row].push_back(Trip(row, mesh.vEdgeIndexNum() + eastHEdge.index(), -0.5 * _muEast * _dyDivdxEast));
                    }
                    else if (mesh.cells(eastHEdge.masterCell).level() < masterCell.level())
                    {
                        int _supSonCellIndex = mesh.cells(eastHEdge.masterCell).supCellIndex * 4;
                        for (int j = 0; j < supSonCells_[_supSonCellIndex].hEdges[0].ind.size(); j++)
                        {
                            momCoeff[row].push_back(Trip(row, mesh.vEdgeIndexNum() + supSonCells_[_supSonCellIndex].hEdges[0].ind[j], -0.5 * _muEast * _dyDivdxEast * supSonCells_[_supSonCellIndex].hEdges[0].coe[j]));
                        }
                    }
                    else if (mesh.cells(eastHEdge.slaveCell).level() < slaveCell.level())
                    {
                        int _supSonCellIndex = mesh.cells(eastHEdge.slaveCell).supCellIndex * 4 + 3;
                        for (int j = 0; j < supSonCells_[_supSonCellIndex].hEdges[1].ind.size(); j++)
                        {
                            momCoeff[row].push_back(Trip(row, mesh.vEdgeIndexNum() + supSonCells_[_supSonCellIndex].hEdges[1].ind[j], -0.5 * _muEast * _dyDivdxEast * supSonCells_[_supSonCellIndex].hEdges[1].coe[j]));
                        }
                    }
                    else if (mesh.cells(eastHEdge.masterCell).level() > masterCell.level() && mesh.cells(eastHEdge.slaveCell).level() > slaveCell.level())
                    {
                        Cell& neiParentCell = mesh.cells(mesh.cells(eastHEdge.masterCell).parent);
                        Cell& sonCell0 = mesh.cells(neiParentCell.son[0]);
                        Cell& sonCell1 = mesh.cells(neiParentCell.son[1]);
                        if constexpr (Axisymmetric)
                        {
                            const auto emit = [&](int column, double coefficient) {
                                momCoeff[row].emplace_back(row, column, coefficient);
                            };
                            HydroNativeRZRowDetail::visitLeafLowSide(mesh, sonCell0, false, -0.5 * 0.5 * _muEast * _dyDivdxEast, emit);
                            HydroNativeRZRowDetail::visitLeafLowSide(mesh, sonCell1, false, -0.5 * 0.5 * _muEast * _dyDivdxEast, emit);
                        }
                        else
                        {
                            momCoeff[row].push_back(Trip(row, mesh.vEdgeIndexNum() + mesh.hEdges(sonCell0.cellToHEdge[0][0]).index(), -0.5 * 0.5 * _muEast * _dyDivdxEast));
                            momCoeff[row].push_back(Trip(row, mesh.vEdgeIndexNum() + mesh.hEdges(sonCell1.cellToHEdge[0][0]).index(), -0.5 * 0.5 * _muEast * _dyDivdxEast));
                        }
                    }
                    else if (mesh.cells(eastHEdge.masterCell).level() > masterCell.level())
                    {
                        int _supCellIndex = mesh.cells(eastHEdge.slaveCell).supCellIndex;
                        for (int j = 0; j < supCells_[_supCellIndex].hEdges[1].ind.size(); j++)
                        {
                            momCoeff[row].push_back(Trip(row, mesh.vEdgeIndexNum() + supCells_[_supCellIndex].hEdges[1].ind[j], -0.5 * _muEast * _dyDivdxEast * supCells_[_supCellIndex].hEdges[1].coe[j]));
                        }
                    }
                    else if (mesh.cells(eastHEdge.slaveCell).level() > slaveCell.level())
                    {
                        int _supCellIndex = mesh.cells(eastHEdge.masterCell).supCellIndex;
                        for (int j = 0; j < supCells_[_supCellIndex].hEdges[0].ind.size(); j++)
                        {
                            momCoeff[row].push_back(Trip(row, mesh.vEdgeIndexNum() + supCells_[_supCellIndex].hEdges[0].ind[j], -0.5 * _muEast * _dyDivdxEast * supCells_[_supCellIndex].hEdges[0].coe[j]));
                        }
                    }
                }
                else if (endVert.supCell != -1)
                {
                    _dyDivdxEast = metrics.shearEast();
                    Cell& _supCell = mesh.cells(endVert.supCell);
                    int _supSonCellIndex = _supCell.supCellIndex * 4;
                    for (int j = 0; j < supSonCells_[_supSonCellIndex].hEdges[1].ind.size(); j++)
                    {
                        momCoeff[row].push_back(Trip(row, mesh.vEdgeIndexNum() + supSonCells_[_supSonCellIndex].hEdges[1].ind[j], -0.5 * _muEast * _dyDivdxEast * supSonCells_[_supSonCellIndex].hEdges[1].coe[j]));
                    }
                }
                if (masterCell.supCellIndex == -1)
                {
                    momCoeff[row].push_back(Trip(row, mesh.vEdgeIndexNum() + mesh.hEdges(masterCell.cellToHEdge[1][0]).index(), -_muNorth * metrics.normal()));
                    momCoeff[row].push_back(Trip(row, mesh.vEdges(masterCell.cellToVEdge[0][0]).index(), metrics.crossWest(0.5 * _muWest)));
                    momCoeff[row].push_back(Trip(row, mesh.vEdges(masterCell.cellToVEdge[1][0]).index(), metrics.crossEast(-0.5 * _muEast)));
                }
                else
                {
                    int _supCellIndex = masterCell.supCellIndex;
                    for (int j = 0; j < supCells_[_supCellIndex].hEdges[1].ind.size(); j++)
                    {
                        momCoeff[row].push_back(Trip(row, mesh.vEdgeIndexNum() + supCells_[_supCellIndex].hEdges[1].ind[j], -_muNorth * metrics.normal() * supCells_[_supCellIndex].hEdges[1].coe[j]));
                    }
                    for (int j = 0; j < supCells_[_supCellIndex].vEdges[0].ind.size(); j++)
                    {
                        momCoeff[row].push_back(Trip(row, supCells_[_supCellIndex].vEdges[0].ind[j], metrics.crossWest(0.5 * _muWest) * supCells_[_supCellIndex].vEdges[0].coe[j]));
                    }
                    for (int j = 0; j < supCells_[_supCellIndex].vEdges[1].ind.size(); j++)
                    {
                        momCoeff[row].push_back(Trip(row, supCells_[_supCellIndex].vEdges[1].ind[j], metrics.crossEast(-0.5 * _muEast) * supCells_[_supCellIndex].vEdges[1].coe[j]));
                    }
                }
                if (slaveCell.supCellIndex == -1)
                {
                    momCoeff[row].push_back(Trip(row, mesh.vEdgeIndexNum() + mesh.hEdges(slaveCell.cellToHEdge[0][0]).index(), -_muSouth * metrics.normal()));
                    momCoeff[row].push_back(Trip(row, mesh.vEdges(slaveCell.cellToVEdge[0][0]).index(), metrics.crossWest(-0.5 * _muWest)));
                    momCoeff[row].push_back(Trip(row, mesh.vEdges(slaveCell.cellToVEdge[1][0]).index(), metrics.crossEast(0.5 * _muEast)));
                }
                else
                {
                    int _supCellIndex = slaveCell.supCellIndex;
                    for (int j = 0; j < supCells_[_supCellIndex].hEdges[0].ind.size(); j++)
                    {
                        momCoeff[row].push_back(Trip(row, mesh.vEdgeIndexNum() + supCells_[_supCellIndex].hEdges[0].ind[j], -_muSouth * metrics.normal() * supCells_[_supCellIndex].hEdges[0].coe[j]));
                    }
                    for (int j = 0; j < supCells_[_supCellIndex].vEdges[0].ind.size(); j++)
                    {
                        momCoeff[row].push_back(Trip(row, supCells_[_supCellIndex].vEdges[0].ind[j], metrics.crossWest(-0.5 * _muWest) * supCells_[_supCellIndex].vEdges[0].coe[j]));
                    }
                    for (int j = 0; j < supCells_[_supCellIndex].vEdges[1].ind.size(); j++)
                    {
                        momCoeff[row].push_back(Trip(row, supCells_[_supCellIndex].vEdges[1].ind[j], metrics.crossEast(0.5 * _muEast) * supCells_[_supCellIndex].vEdges[1].coe[j]));
                    }
                }
                momCoeff[row].push_back(Trip(row, mesh.vEdgeIndexNum() + hEdge.index(), _dA * _rho / dt + (_muSouth + _muNorth) * metrics.normal() + 0.5 * (_muWest * _dyDivdxWest + _muEast * _dyDivdxEast)));
            }
            else if (masterCell.level() < slaveCell.level())
            {
                double _dA = HydroNativeMomentumDetail::volume<Axisymmetric>(hEdge, slaveCell);
                double _rhoMaster = 0.0;
                int supSonCellIndex = masterCell.centerCoord().x < slaveCell.centerCoord().x? masterCell.supCellIndex * 4 + 1 : masterCell.supCellIndex * 4;
                for (int j = 0; j < supSonCells_[supSonCellIndex].center.ind.size(); j++)
                {
                    _rhoMaster += rho.value[supSonCells_[supSonCellIndex].center.ind[j]] * supSonCells_[supSonCellIndex].center.coe[j];
                }
                double _rho = 0.5 * (rho.value[slaveCell.index()] + _rhoMaster);
                double _muSouth = mu.value[slaveCell.index()];
                double _muNorth = 0.0;
                for (int j = 0; j < supSonCells_[supSonCellIndex].center.ind.size(); j++)
                {
                    _muNorth += mu.value[supSonCells_[supSonCellIndex].center.ind[j]] * supSonCells_[supSonCellIndex].center.coe[j];
                }
                double _muWest = vertMu.value[startVert.index()];
                double _muEast = vertMu.value[endVert.index()];
                const auto metrics = HydroNativeMomentumDetail::axialMetrics<Axisymmetric>(hEdge, slaveCell, startVert, endVert);
                double _dyDivdxWest = 0.0;
                double _dyDivdxEast = 0.0;
                if (startVert.vertToHEdge[0] != -1)
                {
                    Edge& westHEdge = mesh.hEdges(startVert.vertToHEdge[0]);
                    _dyDivdxWest = metrics.shearWest();
                    if (mesh.cells(westHEdge.slaveCell).level() == slaveCell.level())
                    {
                        momCoeff[row].push_back(Trip(row, mesh.vEdgeIndexNum() + westHEdge.index(), -0.5 * _muWest * _dyDivdxWest));
                    }
                    else if (mesh.cells(westHEdge.slaveCell).level() < slaveCell.level())
                    {
                        int _supSonCellIndex = mesh.cells(westHEdge.slaveCell).supCellIndex * 4 + 2;
                        for (int j = 0; j < supSonCells_[_supSonCellIndex].hEdges[1].ind.size(); j++)
                        {
                            momCoeff[row].push_back(Trip(row, mesh.vEdgeIndexNum() + supSonCells_[_supSonCellIndex].hEdges[1].ind[j], -0.5 * _muWest * _dyDivdxWest * supSonCells_[_supSonCellIndex].hEdges[1].coe[j]));
                        }
                    }
                }
                if (endVert.vertToHEdge[1] != -1)
                {
                    Edge& eastHEdge = mesh.hEdges(endVert.vertToHEdge[1]);
                    _dyDivdxEast = metrics.shearEast();
                    if (mesh.cells(eastHEdge.slaveCell).level() == slaveCell.level())
                    {
                        momCoeff[row].push_back(Trip(row, mesh.vEdgeIndexNum() + eastHEdge.index(), -0.5 * _muEast * _dyDivdxEast));
                    }
                    else if (mesh.cells(eastHEdge.slaveCell).level() < slaveCell.level())
                    {
                        int _supSonCellIndex = mesh.cells(eastHEdge.slaveCell).supCellIndex * 4 + 3;
                        for (int j = 0; j < supSonCells_[_supSonCellIndex].hEdges[1].ind.size(); j++)
                        {
                            momCoeff[row].push_back(Trip(row, mesh.vEdgeIndexNum() + supSonCells_[_supSonCellIndex].hEdges[1].ind[j], -0.5 * _muEast * _dyDivdxEast * supSonCells_[_supSonCellIndex].hEdges[1].coe[j]));
                        }
                    }
                }
                {
                    for (int j = 0; j < supSonCells_[supSonCellIndex].hEdges[1].ind.size(); j++)
                    {
                        momCoeff[row].push_back(Trip(row, mesh.vEdgeIndexNum() + supSonCells_[supSonCellIndex].hEdges[1].ind[j], -_muNorth * metrics.normal() * supSonCells_[supSonCellIndex].hEdges[1].coe[j]));
                    }
                    for (int j = 0; j < supSonCells_[supSonCellIndex].vEdges[0].ind.size(); j++)
                    {
                        momCoeff[row].push_back(Trip(row, supSonCells_[supSonCellIndex].vEdges[0].ind[j], metrics.crossWest(0.5 * _muWest) * supSonCells_[supSonCellIndex].vEdges[0].coe[j]));
                    }
                    for (int j = 0; j < supSonCells_[supSonCellIndex].vEdges[1].ind.size(); j++)
                    {
                        momCoeff[row].push_back(Trip(row, supSonCells_[supSonCellIndex].vEdges[1].ind[j], metrics.crossEast(-0.5 * _muEast) * supSonCells_[supSonCellIndex].vEdges[1].coe[j]));
                    }
                }
                {
                    momCoeff[row].push_back(Trip(row, mesh.vEdgeIndexNum() + mesh.hEdges(slaveCell.cellToHEdge[0][0]).index(), -_muSouth * metrics.normal()));
                    momCoeff[row].push_back(Trip(row, mesh.vEdges(slaveCell.cellToVEdge[0][0]).index(), metrics.crossWest(-0.5 * _muWest)));
                    momCoeff[row].push_back(Trip(row, mesh.vEdges(slaveCell.cellToVEdge[1][0]).index(), metrics.crossEast(0.5 * _muEast)));
                }
                momCoeff[row].push_back(Trip(row, mesh.vEdgeIndexNum() + hEdge.index(), _dA * _rho / dt + (_muSouth + _muNorth) * metrics.normal() + 0.5 * (_muWest * _dyDivdxWest + _muEast * _dyDivdxEast)));
            }
            else
            {
                double _dA = HydroNativeMomentumDetail::volume<Axisymmetric>(hEdge, masterCell);
                double _rhoSlave = 0.0;
                int supSonCellIndex = slaveCell.centerCoord().x < masterCell.centerCoord().x? slaveCell.supCellIndex * 4 + 2 : slaveCell.supCellIndex * 4 + 3;
                for (int j = 0; j < supSonCells_[supSonCellIndex].center.ind.size(); j++)
                {
                    _rhoSlave += rho.value[supSonCells_[supSonCellIndex].center.ind[j]] * supSonCells_[supSonCellIndex].center.coe[j];
                }
                double _rho = 0.5 * (rho.value[masterCell.index()] + _rhoSlave);
                double _muSouth = 0.0;
                for (int j = 0; j < supSonCells_[supSonCellIndex].center.ind.size(); j++)
                {
                    _muSouth += mu.value[supSonCells_[supSonCellIndex].center.ind[j]] * supSonCells_[supSonCellIndex].center.coe[j];
                }
                double _muNorth = mu.value[masterCell.index()];
                double _muWest = vertMu.value[startVert.index()];
                double _muEast = vertMu.value[endVert.index()];
                const auto metrics = HydroNativeMomentumDetail::axialMetrics<Axisymmetric>(hEdge, masterCell, startVert, endVert);
                double _dyDivdxWest = 0.0;
                double _dyDivdxEast = 0.0;
                if (startVert.vertToHEdge[0] != -1)
                {
                    Edge& westHEdge = mesh.hEdges(startVert.vertToHEdge[0]);
                    _dyDivdxWest = metrics.shearWest();
                    if (mesh.cells(westHEdge.masterCell).level() == masterCell.level())
                    {
                        momCoeff[row].push_back(Trip(row, mesh.vEdgeIndexNum() + westHEdge.index(), -0.5 * _muWest * _dyDivdxWest));
                    }
                    else if (mesh.cells(westHEdge.masterCell).level() < masterCell.level())
                    {
                        int _supSonCellIndex = mesh.cells(westHEdge.masterCell).supCellIndex * 4 + 1;
                        for (int j = 0; j < supSonCells_[_supSonCellIndex].hEdges[0].ind.size(); j++)
                        {
                            momCoeff[row].push_back(Trip(row, mesh.vEdgeIndexNum() + supSonCells_[_supSonCellIndex].hEdges[0].ind[j], -0.5 * _muWest * _dyDivdxWest * supSonCells_[_supSonCellIndex].hEdges[0].coe[j]));
                        }
                    }
                }
                if (endVert.vertToHEdge[1] != -1)
                {
                    Edge& eastHEdge = mesh.hEdges(endVert.vertToHEdge[1]);
                    _dyDivdxEast = metrics.shearEast();
                    if (mesh.cells(eastHEdge.masterCell).level() == masterCell.level())
                    {
                        momCoeff[row].push_back(Trip(row, mesh.vEdgeIndexNum() + eastHEdge.index(), -0.5 * _muEast * _dyDivdxEast));
                    }
                    else if (mesh.cells(eastHEdge.masterCell).level() < masterCell.level())
                    {
                        int _supSonCellIndex = mesh.cells(eastHEdge.masterCell).supCellIndex * 4;
                        for (int j = 0; j < supSonCells_[_supSonCellIndex].hEdges[0].ind.size(); j++)
                        {
                            momCoeff[row].push_back(Trip(row, mesh.vEdgeIndexNum() + supSonCells_[_supSonCellIndex].hEdges[0].ind[j], -0.5 * _muEast * _dyDivdxEast * supSonCells_[_supSonCellIndex].hEdges[0].coe[j]));
                        }
                    }
                }
                {
                    for (int j = 0; j < supSonCells_[supSonCellIndex].hEdges[0].ind.size(); j++)
                    {
                        momCoeff[row].push_back(Trip(row, mesh.vEdgeIndexNum() + supSonCells_[supSonCellIndex].hEdges[0].ind[j], -_muSouth * metrics.normal() * supSonCells_[supSonCellIndex].hEdges[0].coe[j]));
                    }
                    for (int j = 0; j < supSonCells_[supSonCellIndex].vEdges[0].ind.size(); j++)
                    {
                        momCoeff[row].push_back(Trip(row, supSonCells_[supSonCellIndex].vEdges[0].ind[j], metrics.crossWest(-0.5 * _muWest) * supSonCells_[supSonCellIndex].vEdges[0].coe[j]));
                    }
                    for (int j = 0; j < supSonCells_[supSonCellIndex].vEdges[1].ind.size(); j++)
                    {
                        momCoeff[row].push_back(Trip(row, supSonCells_[supSonCellIndex].vEdges[1].ind[j], metrics.crossEast(0.5 * _muEast) * supSonCells_[supSonCellIndex].vEdges[1].coe[j]));
                    }
                }
                {
                    momCoeff[row].push_back(Trip(row, mesh.vEdgeIndexNum() + mesh.hEdges(masterCell.cellToHEdge[1][0]).index(), -_muNorth * metrics.normal()));
                    momCoeff[row].push_back(Trip(row, mesh.vEdges(masterCell.cellToVEdge[0][0]).index(), metrics.crossWest(0.5 * _muWest)));
                    momCoeff[row].push_back(Trip(row, mesh.vEdges(masterCell.cellToVEdge[1][0]).index(), metrics.crossEast(-0.5 * _muEast)));
                }
                momCoeff[row].push_back(Trip(row, mesh.vEdgeIndexNum() + hEdge.index(), _dA * _rho / dt + (_muSouth + _muNorth) * metrics.normal() + 0.5 * (_muWest * _dyDivdxWest + _muEast * _dyDivdxEast)));
            }
        }
        else
        {
            momCoeff[row].push_back(Trip(row, mesh.vEdgeIndexNum() + hEdge.index(), 1.0));
        }
    }
    for (int i = 0; i < mesh.vEdgeIndexNum() + mesh.hEdgeIndexNum(); i++)
    {
        assembly.entries.insert(assembly.entries.end(), momCoeff[i].begin(), momCoeff[i].end());
    }
}

// Shared native CN/AB2 RHS. Upwind decisions, historical interpolation and
// force accumulation keep their original sequence in each supported policy.
template<bool Axisymmetric, typename Mesh>
void Hydrodynamics::buildNativeMomentumRHS(Mesh& mesh, const double& dt, HydroMomentumAssembly& assembly,
    const std::vector<unsigned char>* selectedRows,
    const CellField<double>* bodyX, const CellField<double>* bodyY,bool includeAdvection)
{
    if constexpr (Axisymmetric)
    {
        if (selectedRows && selectedRows->size()!=static_cast<std::size_t>(mesh.vEdgeIndexNum()+mesh.hEdgeIndexNum()))
            throw std::invalid_argument("native RZ RHS row selection has wrong size");
    }
    const auto& nativeBodyX=bodyX?*bodyX:externalFx;
    const auto& nativeBodyY=bodyY?*bodyY:externalFy;
    assembly.rhs.resize(mesh.vEdgeIndexNum() + mesh.hEdgeIndexNum());
    if constexpr (Axisymmetric)
    {
        if (selectedRows) assembly.rhs.setZero();
    }
#pragma omp parallel for num_threads(thread_num)
    for (int i = 0; i < mesh.vEdgeIndexNum(); i++)
    {
        Edge& vEdge = mesh.vEdges(mesh.vEdgeIndexToLabel[i]);
        int row = i;
        if constexpr (Axisymmetric)
        {
            if (selectedRows && !(*selectedRows)[row]) continue;
        }
        if (vEdge.masterCell != -1 && vEdge.slaveCell != -1)
        {
            Cell& masterCell = mesh.cells(vEdge.masterCell);
            Cell& slaveCell = mesh.cells(vEdge.slaveCell);
            Vert& startVert = mesh.verts(vEdge.edgeToVert[0]);
            Vert& endVert = mesh.verts(vEdge.edgeToVert[1]);
            if (masterCell.level() == slaveCell.level())
            {
                if constexpr (Axisymmetric)
                {
                    HydroNativeRZRow::Row nativeRow;
                    if (HydroNativeRZRowDetail::regularRow(mesh, vEdge, true, rho, mu, vertMu, nativeRow))
                    {
                        HydroNativeRZRowDetail::prepareTransverseDensity(nativeRow.materials, true,
                            startVert.index(), endVert.index(), [&](int index) { return vertRho.value[index]; });
                        const int radialCount = mesh.vEdgeIndexNum();
                        const auto current = HydroNativeRZRow::values(nativeRow.indices, [&](int k) {
                            return k < radialCount ? ux.value[k] : uy.value[k - radialCount];
                        });
                        const auto previous = HydroNativeRZRow::values(nativeRow.indices, [&](int k) {
                            return k < radialCount ? uxPre.value[k] : uyPre.value[k - radialCount];
                        });
                        assembly.rhs[row] = HydroNativeRZRow::radialRhs(nativeRow, current, previous,
                            nativeBodyX.value[masterCell.index()], nativeBodyX.value[slaveCell.index()], dt,includeAdvection);
                        continue;
                    }
                }
                double _dA = HydroNativeMomentumDetail::volume<Axisymmetric>(vEdge, masterCell);
                double _rho = 0.5 * (rho.value[masterCell.index()] + rho.value[slaveCell.index()]);
                double _muWest = mu.value[slaveCell.index()];
                double _muEast = mu.value[masterCell.index()];
                double _muSouth = vertMu.value[startVert.index()];
                double _muNorth = vertMu.value[endVert.index()];
                double _mu = 2.0 * _muEast * _muWest / (_muEast + _muWest);
                double _rhoWest = rho.value[slaveCell.index()];
                double _rhoEast = rho.value[masterCell.index()];
                double _rhoSouth = vertRho.value[startVert.index()];
                double _rhoNorth = vertRho.value[endVert.index()];
                const auto metrics = HydroNativeMomentumDetail::radialFluxMetrics<Axisymmetric, HydroNativeMomentumDetail::FaceRefinement::SameLevel>(vEdge, masterCell, slaveCell);
                double _vSouthMaster = 0.0;
                double _vSouthSlave = 0.0;
                double _vNorthMaster = 0.0;
                double _vNorthSlave = 0.0;
                double _vSouthMasterPre = 0.0;
                double _vSouthSlavePre = 0.0;
                double _vNorthMasterPre = 0.0;
                double _vNorthSlavePre = 0.0;
                double _uWestUpwind = 0.0;
                double _uWestUpwindPre = 0.0;
                double _diffWest = 0.0;
                double _uWest = 0.0;
                double _uWestPre = 0.0;
                if (slaveCell.supCellIndex == -1)
                {
                    Edge& westVEdge00 = mesh.vEdges(slaveCell.cellToVEdge[0][0]);
                    _uWest = ux.value[westVEdge00.index()];
                    _uWestPre = uxPre.value[westVEdge00.index()];
                    Edge& southHEdge00 = mesh.hEdges(slaveCell.cellToHEdge[0][0]);
                    Edge& northHEdge10 = mesh.hEdges(slaveCell.cellToHEdge[1][0]);
                    _vSouthSlave = uy.value[southHEdge00.index()];
                    _vSouthSlavePre = uyPre.value[southHEdge00.index()];
                    _vNorthSlave = uy.value[northHEdge10.index()];
                    _vNorthSlavePre = uyPre.value[northHEdge10.index()];
                }
                else
                {
                    int _supCellIndex = slaveCell.supCellIndex;
                    for (int j = 0; j < supCells_[_supCellIndex].vEdges[0].ind.size(); j++)
                    {
                        _uWest += ux.value[supCells_[_supCellIndex].vEdges[0].ind[j]] * supCells_[_supCellIndex].vEdges[0].coe[j];
                        _uWestPre += uxPre.value[supCells_[_supCellIndex].vEdges[0].ind[j]] * supCells_[_supCellIndex].vEdges[0].coe[j];
                    }
                    for (int j = 0; j < supCells_[_supCellIndex].hEdges[0].ind.size(); j++)
                    {
                        _vSouthSlave += uy.value[supCells_[_supCellIndex].hEdges[0].ind[j]] * supCells_[_supCellIndex].hEdges[0].coe[j];
                        _vSouthSlavePre += uyPre.value[supCells_[_supCellIndex].hEdges[0].ind[j]] * supCells_[_supCellIndex].hEdges[0].coe[j];
                    }
                    for (int j = 0; j < supCells_[_supCellIndex].hEdges[1].ind.size(); j++)
                    {
                        _vNorthSlave += uy.value[supCells_[_supCellIndex].hEdges[1].ind[j]] * supCells_[_supCellIndex].hEdges[1].coe[j];
                        _vNorthSlavePre += uyPre.value[supCells_[_supCellIndex].hEdges[1].ind[j]] * supCells_[_supCellIndex].hEdges[1].coe[j];
                    }
                }
                double _qWest = _rhoWest * 0.5 * (ux.value[vEdge.index()] + _uWest) * metrics.westArea();
                double _qWestPre = _rhoWest * 0.5 * (uxPre.value[vEdge.index()] + _uWestPre) * metrics.westArea();
                {
                    _uWestUpwind = _qWest > 0.0 ? _uWest : ux.value[vEdge.index()];
                }
                {
                    _uWestUpwindPre = _qWestPre > 0.0 ? _uWestPre : uxPre.value[vEdge.index()];
                }
                _diffWest = 2.0 * _muWest * (ux.value[vEdge.index()] - _uWest) / metrics.dx;
                double _uEastUpwind = 0.0;
                double _uEastUpwindPre = 0.0;
                double _diffEast = 0.0;
                double _uEast = 0.0;
                double _uEastPre = 0.0;
                if (masterCell.supCellIndex == -1)
                {
                    Edge& eastEdge10 = mesh.vEdges(masterCell.cellToVEdge[1][0]);
                    _uEast = ux.value[eastEdge10.index()];
                    _uEastPre = uxPre.value[eastEdge10.index()];
                    Edge& southHEdge00 = mesh.hEdges(masterCell.cellToHEdge[0][0]);
                    Edge& northHEdge10 = mesh.hEdges(masterCell.cellToHEdge[1][0]);
                    _vSouthMaster = uy.value[southHEdge00.index()];
                    _vSouthMasterPre = uyPre.value[southHEdge00.index()];
                    _vNorthMaster = uy.value[northHEdge10.index()];
                    _vNorthMasterPre = uyPre.value[northHEdge10.index()];
                }
                else
                {
                    int _supCellIndex = masterCell.supCellIndex;
                    for (int j = 0; j < supCells_[_supCellIndex].vEdges[1].ind.size(); j++)
                    {
                        _uEast += ux.value[supCells_[_supCellIndex].vEdges[1].ind[j]] * supCells_[_supCellIndex].vEdges[1].coe[j];
                        _uEastPre += uxPre.value[supCells_[_supCellIndex].vEdges[1].ind[j]] * supCells_[_supCellIndex].vEdges[1].coe[j];
                    }
                    for (int j = 0; j < supCells_[_supCellIndex].hEdges[0].ind.size(); j++)
                    {
                        _vSouthMaster += uy.value[supCells_[_supCellIndex].hEdges[0].ind[j]] * supCells_[_supCellIndex].hEdges[0].coe[j];
                        _vSouthMasterPre += uyPre.value[supCells_[_supCellIndex].hEdges[0].ind[j]] * supCells_[_supCellIndex].hEdges[0].coe[j];
                    }
                    for (int j = 0; j < supCells_[_supCellIndex].hEdges[1].ind.size(); j++)
                    {
                        _vNorthMaster += uy.value[supCells_[_supCellIndex].hEdges[1].ind[j]] * supCells_[_supCellIndex].hEdges[1].coe[j];
                        _vNorthMasterPre += uyPre.value[supCells_[_supCellIndex].hEdges[1].ind[j]] * supCells_[_supCellIndex].hEdges[1].coe[j];
                    }
                }
                double _qEast = _rhoEast * 0.5 * (ux.value[vEdge.index()] + _uEast) * metrics.eastArea();
                double _qEastPre = _rhoEast * 0.5 * (uxPre.value[vEdge.index()] + _uEastPre) * metrics.eastArea();
                {
                    _uEastUpwind = _qEast > 0.0 ? ux.value[vEdge.index()] : _uEast;
                }
                {
                    _uEastUpwindPre = _qEastPre > 0.0 ? uxPre.value[vEdge.index()] : _uEastPre;
                }
                _diffEast = 2.0 * _muEast * (_uEast - ux.value[vEdge.index()]) / metrics.dx;
                double _vSouth = 0.5 * (_vSouthMaster + _vSouthSlave);
                double _vSouthPre = 0.5 * (_vSouthMasterPre + _vSouthSlavePre);
                double _qSouth = _rhoSouth * _vSouth * metrics.axialArea();
                double _qSouthPre = _rhoSouth * _vSouthPre * metrics.axialArea();
                double _uSouthUpwind = 0.0;
                double _uSouthUpwindPre = 0.0;
                double _diffSouth = 0.0;
                double _uSouth = 0.0;
                double _uSouthPre = 0.0;
                if (startVert.vertToVEdge[0] != -1)
                {
                    Edge& southVEdge = mesh.vEdges(startVert.vertToVEdge[0]);
                    _uSouth = ux.value[southVEdge.index()];
                    _uSouthPre = uxPre.value[southVEdge.index()];
                    _diffSouth = _muSouth * ((ux.value[vEdge.index()] - _uSouth) / (vEdge.midpointCoord().y - southVEdge.midpointCoord().y) + (_vSouthMaster - _vSouthSlave) / metrics.dx);
                }
                else if (startVert.supCell != -1)
                {
                    Cell _supCell = mesh.cells(startVert.supCell);
                    int _supSonCellIndex = _supCell.supCellIndex * 4 + 2;
                    for (int j = 0; j < supSonCells_[_supSonCellIndex].vEdges[0].ind.size(); j++)
                    {
                        _uSouth += ux.value[supSonCells_[_supSonCellIndex].vEdges[0].ind[j]] * supSonCells_[_supSonCellIndex].vEdges[0].coe[j];
                        _uSouthPre += uxPre.value[supSonCells_[_supSonCellIndex].vEdges[0].ind[j]] * supSonCells_[_supSonCellIndex].vEdges[0].coe[j];
                    }
                    _diffSouth = _muSouth * ((ux.value[vEdge.index()] - _uSouth) / metrics.dy + (_vSouthMaster - _vSouthSlave) / metrics.dx);
                }
                else
                {
                    _diffSouth = _muSouth * (_vSouthMaster - _vSouthSlave) / metrics.dx;
                }
                {
                    _uSouthUpwind = _qSouth > 0.0 ? _uSouth : ux.value[vEdge.index()];
                }
                {
                    _uSouthUpwindPre = _qSouthPre > 0.0 ? _uSouthPre : uxPre.value[vEdge.index()];
                }
                double _vNorth = 0.5 * (_vNorthMaster + _vNorthSlave);
                double _vNorthPre = 0.5 * (_vNorthMasterPre + _vNorthSlavePre);
                double _qNorth = _rhoNorth * _vNorth * metrics.axialArea();
                double _qNorthPre = _rhoNorth * _vNorthPre * metrics.axialArea();
                double _uNorthUpwind = 0.0;
                double _uNorthUpwindPre = 0.0;
                double _diffNorth = 0.0;
                double _uNorth = 0.0;
                double _uNorthPre = 0.0;
                if (endVert.vertToVEdge[1] != -1)
                {
                    Edge& northVEdge = mesh.vEdges(endVert.vertToVEdge[1]);
                    _uNorth = ux.value[northVEdge.index()];
                    _uNorthPre = uxPre.value[northVEdge.index()];
                    _diffNorth = _muNorth * ((_uNorth - ux.value[vEdge.index()]) / (northVEdge.midpointCoord().y - vEdge.midpointCoord().y) + (_vNorthMaster - _vNorthSlave) / metrics.dx);
                }
                else if (endVert.supCell != -1)
                {
                    Cell _supCell = mesh.cells(endVert.supCell);
                    int _supSonCellIndex = _supCell.supCellIndex * 4;
                    for (int j = 0; j < supSonCells_[_supSonCellIndex].vEdges[1].ind.size(); j++)
                    {
                        _uNorth += ux.value[supSonCells_[_supSonCellIndex].vEdges[1].ind[j]] * supSonCells_[_supSonCellIndex].vEdges[1].coe[j];
                        _uNorthPre += uxPre.value[supSonCells_[_supSonCellIndex].vEdges[1].ind[j]] * supSonCells_[_supSonCellIndex].vEdges[1].coe[j];
                    }
                    _diffNorth = _muNorth * ((_uNorth - ux.value[vEdge.index()]) / metrics.dy + (_vNorthMaster - _vNorthSlave) / metrics.dx);
                }
                else
                {
                    _diffNorth = _muNorth * (_vNorthMaster - _vNorthSlave) / metrics.dx;
                }
                {
                    _uNorthUpwind = _qNorth > 0.0 ? ux.value[vEdge.index()] : _uNorth;
                }
                {
                    _uNorthUpwindPre = _qNorthPre > 0.0 ? uxPre.value[vEdge.index()] : _uNorthPre;
                }
                double _n = _qEast * _uEastUpwind - _qWest * _uWestUpwind + _qNorth * _uNorthUpwind - _qSouth * _uSouthUpwind;
                double _nPre = _qEastPre * _uEastUpwindPre - _qWestPre * _uWestUpwindPre + _qNorthPre * _uNorthUpwindPre - _qSouthPre * _uSouthUpwindPre;
                if(!includeAdvection) _n=_nPre=0;
                double _nDiff = metrics.diffusionNet(_diffEast, _diffWest, _diffNorth, _diffSouth);
                assembly.rhs[row] = _dA * _rho / dt * ux.value[vEdge.index()] - 1.5 * _n + 0.5 * (_nPre + _nDiff) + _dA * 0.5 * (nativeBodyX.value[masterCell.index()] + nativeBodyX.value[slaveCell.index()]);
            }
            else if (masterCell.level() < slaveCell.level())
            {
                double _dA = HydroNativeMomentumDetail::volume<Axisymmetric>(vEdge, slaveCell);
                double _rhoMaster = 0.0;
                int supSonCellIndex = masterCell.centerCoord().y < slaveCell.centerCoord().y? masterCell.supCellIndex * 4 + 3 : masterCell.supCellIndex * 4;
                for (int j = 0; j < supSonCells_[supSonCellIndex].center.ind.size(); j++)
                {
                    _rhoMaster += rho.value[supSonCells_[supSonCellIndex].center.ind[j]] * supSonCells_[supSonCellIndex].center.coe[j];
                }
                double _rho = 0.5 * (rho.value[slaveCell.index()] + _rhoMaster);
                double _muWest = mu.value[slaveCell.index()];
                double _muEast = 0.0;
                for (int j = 0; j < supSonCells_[supSonCellIndex].center.ind.size(); j++)
                {
                    _muEast += mu.value[supSonCells_[supSonCellIndex].center.ind[j]] * supSonCells_[supSonCellIndex].center.coe[j];
                }
                double _muSouth = vertMu.value[startVert.index()];
                double _muNorth = vertMu.value[endVert.index()];
                double _mu = 2.0 * _muEast * _muWest / (_muEast + _muWest);
                double _rhoWest = rho.value[slaveCell.index()];
                double _rhoEast = 0.0;
                for (int j = 0; j < supSonCells_[supSonCellIndex].center.ind.size(); j++)
                {
                    _rhoEast += rho.value[supSonCells_[supSonCellIndex].center.ind[j]] * supSonCells_[supSonCellIndex].center.coe[j];
                }
                double _rhoSouth = vertRho.value[startVert.index()];
                double _rhoNorth = vertRho.value[endVert.index()];
                const auto metrics = HydroNativeMomentumDetail::radialFluxMetrics<Axisymmetric, HydroNativeMomentumDetail::FaceRefinement::MasterCoarse>(vEdge, masterCell, slaveCell);
                double _vSouthMaster = 0.0;
                double _vSouthSlave = 0.0;
                double _vNorthMaster = 0.0;
                double _vNorthSlave = 0.0;
                double _vSouthMasterPre = 0.0;
                double _vSouthSlavePre = 0.0;
                double _vNorthMasterPre = 0.0;
                double _vNorthSlavePre = 0.0;
                double _uWestUpwind = 0.0;
                double _uWestUpwindPre = 0.0;
                double _diffWest = 0.0;
                double _uWest = 0.0;
                double _uWestPre = 0.0;
                {
                    Edge& westVEdge00 = mesh.vEdges(slaveCell.cellToVEdge[0][0]);
                    _uWest = ux.value[westVEdge00.index()];
                    _uWestPre = uxPre.value[westVEdge00.index()];
                    Edge& southHEdge00 = mesh.hEdges(slaveCell.cellToHEdge[0][0]);
                    Edge& northHEdge10 = mesh.hEdges(slaveCell.cellToHEdge[1][0]);
                    _vSouthSlave = uy.value[southHEdge00.index()];
                    _vSouthSlavePre = uyPre.value[southHEdge00.index()];
                    _vNorthSlave = uy.value[northHEdge10.index()];
                    _vNorthSlavePre = uyPre.value[northHEdge10.index()];
                }
                double _qWest = _rhoWest * 0.5 * (ux.value[vEdge.index()] + _uWest) * metrics.westArea();
                double _qWestPre = _rhoWest * 0.5 * (uxPre.value[vEdge.index()] + _uWestPre) * metrics.westArea();
                {
                    _uWestUpwind = _qWest > 0.0 ? _uWest : ux.value[vEdge.index()];
                }
                {
                    _uWestUpwindPre = _qWestPre > 0.0 ? _uWestPre : uxPre.value[vEdge.index()];
                }
                _diffWest = 2.0 * _muWest * (ux.value[vEdge.index()] - _uWest) / metrics.dx;
                double _uEastUpwind = 0.0;
                double _uEastUpwindPre = 0.0;
                double _diffEast = 0.0;
                double _uEast = 0.0;
                double _uEastPre = 0.0;
                {
                    for (int j = 0; j < supSonCells_[supSonCellIndex].vEdges[1].ind.size(); j++)
                    {
                        _uEast += ux.value[supSonCells_[supSonCellIndex].vEdges[1].ind[j]] * supSonCells_[supSonCellIndex].vEdges[1].coe[j];
                        _uEastPre += uxPre.value[supSonCells_[supSonCellIndex].vEdges[1].ind[j]] * supSonCells_[supSonCellIndex].vEdges[1].coe[j];
                    }
                    for (int j = 0; j < supSonCells_[supSonCellIndex].hEdges[0].ind.size(); j++)
                    {
                        _vSouthMaster += uy.value[supSonCells_[supSonCellIndex].hEdges[0].ind[j]] * supSonCells_[supSonCellIndex].hEdges[0].coe[j];
                        _vSouthMasterPre += uyPre.value[supSonCells_[supSonCellIndex].hEdges[0].ind[j]] * supSonCells_[supSonCellIndex].hEdges[0].coe[j];
                    }
                    for (int j = 0; j < supSonCells_[supSonCellIndex].hEdges[1].ind.size(); j++)
                    {
                        _vNorthMaster += uy.value[supSonCells_[supSonCellIndex].hEdges[1].ind[j]] * supSonCells_[supSonCellIndex].hEdges[1].coe[j];
                        _vNorthMasterPre += uyPre.value[supSonCells_[supSonCellIndex].hEdges[1].ind[j]] * supSonCells_[supSonCellIndex].hEdges[1].coe[j];
                    }
                }
                double _qEast = _rhoEast * 0.5 * (ux.value[vEdge.index()] + _uEast) * metrics.eastArea();
                double _qEastPre = _rhoEast * 0.5 * (uxPre.value[vEdge.index()] + _uEastPre) * metrics.eastArea();
                {
                    _uEastUpwind = _qEast > 0.0 ? ux.value[vEdge.index()] : _uEast;
                }
                {
                    _uEastUpwindPre = _qEastPre > 0.0 ? uxPre.value[vEdge.index()] : _uEastPre;
                }
                _diffEast = 2.0 * _muEast * (_uEast - ux.value[vEdge.index()]) / metrics.dx;
                double _vSouth = 0.5 * (_vSouthMaster + _vSouthSlave);
                double _vSouthPre = 0.5 * (_vSouthMasterPre + _vSouthSlavePre);
                double _qSouth = _rhoSouth * _vSouth * metrics.axialArea();
                double _qSouthPre = _rhoSouth * _vSouthPre * metrics.axialArea();
                double _uSouthUpwind = 0.0;
                double _uSouthUpwindPre = 0.0;
                double _diffSouth = 0.0;
                double _uSouth = 0.0;
                double _uSouthPre = 0.0;
                if (startVert.vertToVEdge[0] != -1)
                {
                    Edge& southVEdge = mesh.vEdges(startVert.vertToVEdge[0]);
                    _uSouth = ux.value[southVEdge.index()];
                    _uSouthPre = uxPre.value[southVEdge.index()];
                    _diffSouth = _muSouth * ((ux.value[vEdge.index()] - _uSouth) / (vEdge.midpointCoord().y - southVEdge.midpointCoord().y) + (_vSouthMaster - _vSouthSlave) / metrics.dx);
                }
                else
                {
                    _diffSouth = _muSouth * (_vSouthMaster - _vSouthSlave) / metrics.dx;
                }
                {
                    _uSouthUpwind = _qSouth > 0.0 ? _uSouth : ux.value[vEdge.index()];
                }
                {
                    _uSouthUpwindPre = _qSouthPre > 0.0 ? _uSouthPre : uxPre.value[vEdge.index()];
                }
                double _vNorth = 0.5 * (_vNorthMaster + _vNorthSlave);
                double _vNorthPre = 0.5 * (_vNorthMasterPre + _vNorthSlavePre);
                double _qNorth = _rhoNorth * _vNorth * metrics.axialArea();
                double _qNorthPre = _rhoNorth * _vNorthPre * metrics.axialArea();
                double _uNorthUpwind = 0.0;
                double _uNorthUpwindPre = 0.0;
                double _diffNorth = 0.0;
                double _uNorth = 0.0;
                double _uNorthPre = 0.0;
                if (endVert.vertToVEdge[1] != -1)
                {
                    Edge& northVEdge = mesh.vEdges(endVert.vertToVEdge[1]);
                    _uNorth = ux.value[northVEdge.index()];
                    _uNorthPre = uxPre.value[northVEdge.index()];
                    _diffNorth = _muNorth * ((_uNorth - ux.value[vEdge.index()]) / (northVEdge.midpointCoord().y - vEdge.midpointCoord().y) + (_vNorthMaster - _vNorthSlave) / metrics.dx);
                }
                else
                {
                    _diffNorth = _muNorth * (_vNorthMaster - _vNorthSlave) / metrics.dx;
                }
                {
                    _uNorthUpwind = _qNorth > 0.0 ? ux.value[vEdge.index()] : _uNorth;
                }
                {
                    _uNorthUpwindPre = _qNorthPre > 0.0 ? uxPre.value[vEdge.index()] : _uNorthPre;
                }
                double _n = _qEast * _uEastUpwind - _qWest * _uWestUpwind + _qNorth * _uNorthUpwind - _qSouth * _uSouthUpwind;
                double _nPre = _qEastPre * _uEastUpwindPre - _qWestPre * _uWestUpwindPre + _qNorthPre * _uNorthUpwindPre - _qSouthPre * _uSouthUpwindPre;
                if(!includeAdvection) _n=_nPre=0;
                double _nDiff = metrics.diffusionNet(_diffEast, _diffWest, _diffNorth, _diffSouth);
                double _FxMaster = 0.0;
                for (int j = 0; j < supSonCells_[supSonCellIndex].center.ind.size(); j++)
                {
                    _FxMaster += nativeBodyX.value[supSonCells_[supSonCellIndex].center.ind[j]] * supSonCells_[supSonCellIndex].center.coe[j];
                }
                assembly.rhs[row] = _dA * _rho / dt * ux.value[vEdge.index()] - 1.5 * _n + 0.5 * (_nPre + _nDiff) + _dA * 0.5 * (nativeBodyX.value[slaveCell.index()] + _FxMaster);
            }
            else
            {
                double _dA = HydroNativeMomentumDetail::volume<Axisymmetric>(vEdge, masterCell);
                double _rhoSlave = 0.0;
                int supSonCellIndex = slaveCell.centerCoord().y < masterCell.centerCoord().y? slaveCell.supCellIndex * 4 + 2 : slaveCell.supCellIndex * 4 + 1;
                for (int j = 0; j < supSonCells_[supSonCellIndex].center.ind.size(); j++)
                {
                    _rhoSlave += rho.value[supSonCells_[supSonCellIndex].center.ind[j]] * supSonCells_[supSonCellIndex].center.coe[j];
                }
                double _rho = 0.5 * (rho.value[masterCell.index()] + _rhoSlave);
                double _muWest = 0.0;
                for (int j = 0; j < supSonCells_[supSonCellIndex].center.ind.size(); j++)
                {
                    _muWest += mu.value[supSonCells_[supSonCellIndex].center.ind[j]] * supSonCells_[supSonCellIndex].center.coe[j];
                }
                double _muEast = mu.value[masterCell.index()];
                double _muSouth = vertMu.value[startVert.index()];
                double _muNorth = vertMu.value[endVert.index()];
                double _mu = 2.0 * _muEast * _muWest / (_muEast + _muWest);
                double _rhoWest = 0.0;
                for (int j = 0; j < supSonCells_[supSonCellIndex].center.ind.size(); j++)
                {
                    _rhoWest += rho.value[supSonCells_[supSonCellIndex].center.ind[j]] * supSonCells_[supSonCellIndex].center.coe[j];
                }
                double _rhoEast = rho.value[masterCell.index()];
                double _rhoSouth = vertRho.value[startVert.index()];
                double _rhoNorth = vertRho.value[endVert.index()];
                const auto metrics = HydroNativeMomentumDetail::radialFluxMetrics<Axisymmetric, HydroNativeMomentumDetail::FaceRefinement::SlaveCoarse>(vEdge, masterCell, slaveCell);
                double _vSouthMaster = 0.0;
                double _vSouthSlave = 0.0;
                double _vNorthMaster = 0.0;
                double _vNorthSlave = 0.0;
                double _vSouthMasterPre = 0.0;
                double _vSouthSlavePre = 0.0;
                double _vNorthMasterPre = 0.0;
                double _vNorthSlavePre = 0.0;
                double _uWestUpwind = 0.0;
                double _uWestUpwindPre = 0.0;
                double _diffWest = 0.0;
                double _uWest = 0.0;
                double _uWestPre = 0.0;
                {
                    for (int j = 0; j < supSonCells_[supSonCellIndex].vEdges[0].ind.size(); j++)
                    {
                        _uWest += ux.value[supSonCells_[supSonCellIndex].vEdges[0].ind[j]] * supSonCells_[supSonCellIndex].vEdges[0].coe[j];
                        _uWestPre += uxPre.value[supSonCells_[supSonCellIndex].vEdges[0].ind[j]] * supSonCells_[supSonCellIndex].vEdges[0].coe[j];
                    }
                    for (int j = 0; j < supSonCells_[supSonCellIndex].hEdges[0].ind.size(); j++)
                    {
                        _vSouthSlave += uy.value[supSonCells_[supSonCellIndex].hEdges[0].ind[j]] * supSonCells_[supSonCellIndex].hEdges[0].coe[j];
                        _vSouthSlavePre += uyPre.value[supSonCells_[supSonCellIndex].hEdges[0].ind[j]] * supSonCells_[supSonCellIndex].hEdges[0].coe[j];
                    }
                    for (int j = 0; j < supSonCells_[supSonCellIndex].hEdges[1].ind.size(); j++)
                    {
                        _vNorthSlave += uy.value[supSonCells_[supSonCellIndex].hEdges[1].ind[j]] * supSonCells_[supSonCellIndex].hEdges[1].coe[j];
                        _vNorthSlavePre += uyPre.value[supSonCells_[supSonCellIndex].hEdges[1].ind[j]] * supSonCells_[supSonCellIndex].hEdges[1].coe[j];
                    }
                }
                double _qWest = _rhoWest * 0.5 * (ux.value[vEdge.index()] + _uWest) * metrics.westArea();
                double _qWestPre = _rhoWest * 0.5 * (uxPre.value[vEdge.index()] + _uWestPre) * metrics.westArea();
                {
                    _uWestUpwind = _qWest > 0.0 ? _uWest : ux.value[vEdge.index()];
                }
                {
                    _uWestUpwindPre = _qWestPre > 0.0 ? _uWestPre : uxPre.value[vEdge.index()];
                }
                _diffWest = 2.0 * _muWest * (ux.value[vEdge.index()] - _uWest) / metrics.dx;
                double _uEastUpwind = 0.0;
                double _uEastUpwindPre = 0.0;
                double _diffEast = 0.0;
                double _uEast = 0.0;
                double _uEastPre = 0.0;
                {
                    Edge& eastEdge10 = mesh.vEdges(masterCell.cellToVEdge[1][0]);
                    _uEast = ux.value[eastEdge10.index()];
                    _uEastPre = uxPre.value[eastEdge10.index()];
                    Edge& southHEdge00 = mesh.hEdges(masterCell.cellToHEdge[0][0]);
                    Edge& northHEdge10 = mesh.hEdges(masterCell.cellToHEdge[1][0]);
                    _vSouthMaster = uy.value[southHEdge00.index()];
                    _vSouthMasterPre = uyPre.value[southHEdge00.index()];
                    _vNorthMaster = uy.value[northHEdge10.index()];
                    _vNorthMasterPre = uyPre.value[northHEdge10.index()];
                }
                double _qEast = _rhoEast * 0.5 * (ux.value[vEdge.index()] + _uEast) * metrics.eastArea();
                double _qEastPre = _rhoEast * 0.5 * (uxPre.value[vEdge.index()] + _uEastPre) * metrics.eastArea();
                {
                    _uEastUpwind = _qEast > 0.0 ? ux.value[vEdge.index()] : _uEast;
                }
                {
                    _uEastUpwindPre = _qEastPre > 0.0 ? uxPre.value[vEdge.index()] : _uEastPre;
                }
                _diffEast = 2.0 * _muEast * (_uEast - ux.value[vEdge.index()]) / metrics.dx;
                double _vSouth = 0.5 * (_vSouthMaster + _vSouthSlave);
                double _vSouthPre = 0.5 * (_vSouthMasterPre + _vSouthSlavePre);
                double _qSouth = _rhoSouth * _vSouth * metrics.axialArea();
                double _qSouthPre = _rhoSouth * _vSouthPre * metrics.axialArea();
                double _uSouthUpwind = 0.0;
                double _uSouthUpwindPre = 0.0;
                double _diffSouth = 0.0;
                double _uSouth = 0.0;
                double _uSouthPre = 0.0;
                if (startVert.vertToVEdge[0] != -1)
                {
                    Edge& southVEdge = mesh.vEdges(startVert.vertToVEdge[0]);
                    _uSouth = ux.value[southVEdge.index()];
                    _uSouthPre = uxPre.value[southVEdge.index()];
                    _diffSouth = _muSouth * ((ux.value[vEdge.index()] - _uSouth) / (vEdge.midpointCoord().y - southVEdge.midpointCoord().y) + (_vSouthMaster - _vSouthSlave) / metrics.dx);
                }
                else
                {
                    _diffSouth = _muSouth * (_vSouthMaster - _vSouthSlave) / metrics.dx;
                }
                {
                    _uSouthUpwind = _qSouth > 0.0 ? _uSouth : ux.value[vEdge.index()];
                }
                {
                    _uSouthUpwindPre = _qSouthPre > 0.0 ? _uSouthPre : uxPre.value[vEdge.index()];
                }
                double _vNorth = 0.5 * (_vNorthMaster + _vNorthSlave);
                double _vNorthPre = 0.5 * (_vNorthMasterPre + _vNorthSlavePre);
                double _qNorth = _rhoNorth * _vNorth * metrics.axialArea();
                double _qNorthPre = _rhoNorth * _vNorthPre * metrics.axialArea();
                double _uNorthUpwind = 0.0;
                double _uNorthUpwindPre = 0.0;
                double _diffNorth = 0.0;
                double _uNorth = 0.0;
                double _uNorthPre = 0.0;
                if (endVert.vertToVEdge[1] != -1)
                {
                    Edge& northVEdge = mesh.vEdges(endVert.vertToVEdge[1]);
                    _uNorth = ux.value[northVEdge.index()];
                    _uNorthPre = uxPre.value[northVEdge.index()];
                    _diffNorth = _muNorth * ((_uNorth - ux.value[vEdge.index()]) / (northVEdge.midpointCoord().y - vEdge.midpointCoord().y) + (_vNorthMaster - _vNorthSlave) / metrics.dx);
                }
                else
                {
                    _diffNorth = _muNorth * (_vNorthMaster - _vNorthSlave) / metrics.dx;
                }
                {
                    _uNorthUpwind = _qNorth > 0.0 ? ux.value[vEdge.index()] : _uNorth;
                }
                {
                    _uNorthUpwindPre = _qNorthPre > 0.0 ? uxPre.value[vEdge.index()] : _uNorthPre;
                }
                double _n = _qEast * _uEastUpwind - _qWest * _uWestUpwind + _qNorth * _uNorthUpwind - _qSouth * _uSouthUpwind;
                double _nPre = _qEastPre * _uEastUpwindPre - _qWestPre * _uWestUpwindPre + _qNorthPre * _uNorthUpwindPre - _qSouthPre * _uSouthUpwindPre;
                if(!includeAdvection) _n=_nPre=0;
                double _nDiff = metrics.diffusionNet(_diffEast, _diffWest, _diffNorth, _diffSouth);
                double _FxSlave = 0.0;
                for (int j = 0; j < supSonCells_[supSonCellIndex].center.ind.size(); j++)
                {
                    _FxSlave += nativeBodyX.value[supSonCells_[supSonCellIndex].center.ind[j]] * supSonCells_[supSonCellIndex].center.coe[j];
                }
                assembly.rhs[row] = _dA * _rho / dt * ux.value[vEdge.index()] - 1.5 * _n + 0.5 * (_nPre + _nDiff) + _dA * 0.5 * (nativeBodyX.value[masterCell.index()] + _FxSlave);
            }
        }
        else if (vEdge.masterCell == -1)
        {
            if constexpr (Axisymmetric)
            {
                assembly.rhs[row] = ux.value[i];
            }
            else
            {
                assembly.rhs[row] = 0.0;
            }
        }
        else
        {
            if constexpr (Axisymmetric)
            {
                assembly.rhs[row] = ux.value[i];
            }
            else
            {
                assembly.rhs[row] = 0.0;
            }
        }
    }
#pragma omp parallel for num_threads(thread_num)
    for (int i = 0; i < mesh.hEdgeIndexNum(); i++)
    {
        Edge& hEdge = mesh.hEdges(mesh.hEdgeIndexToLabel[i]);
        int row = mesh.vEdgeIndexNum() + i;
        if constexpr (Axisymmetric)
        {
            if (selectedRows && !(*selectedRows)[row]) continue;
        }
        if (!activeOptions_.closedDomainProjection && fabs(hEdge.midpointCoord().y) < 1.0e-4)
        {
            assembly.rhs[row] = 0.0;
            continue;
        }
        if (hEdge.masterCell != -1 && hEdge.slaveCell != -1)
        {
            Cell& masterCell = mesh.cells(hEdge.masterCell);
            Cell& slaveCell = mesh.cells(hEdge.slaveCell);
            Vert& startVert = mesh.verts(hEdge.edgeToVert[0]);
            Vert& endVert = mesh.verts(hEdge.edgeToVert[1]);
            if (masterCell.level() == slaveCell.level())
            {
                if constexpr (Axisymmetric)
                {
                    HydroNativeRZRow::Row nativeRow;
                    if (HydroNativeRZRowDetail::regularRow(mesh, hEdge, false, rho, mu, vertMu, nativeRow))
                    {
                        HydroNativeRZRowDetail::prepareTransverseDensity(nativeRow.materials, false,
                            startVert.index(), endVert.index(), [&](int index) { return vertRho.value[index]; });
                        const int radialCount = mesh.vEdgeIndexNum();
                        const auto current = HydroNativeRZRow::values(nativeRow.indices, [&](int k) {
                            return k < radialCount ? ux.value[k] : uy.value[k - radialCount];
                        });
                        const auto previous = HydroNativeRZRow::values(nativeRow.indices, [&](int k) {
                            return k < radialCount ? uxPre.value[k] : uyPre.value[k - radialCount];
                        });
                        assembly.rhs[row] = HydroNativeRZRow::axialRhs(nativeRow, current, previous,
                            nativeBodyY.value[masterCell.index()], nativeBodyY.value[slaveCell.index()], dt,includeAdvection);
                        continue;
                    }
                }
                double _dA = HydroNativeMomentumDetail::volume<Axisymmetric>(hEdge, masterCell);
                double _rho = 0.5 * (rho.value[masterCell.index()] + rho.value[slaveCell.index()]);
                double _dxDivdy = masterCell.size().x / masterCell.size().y;
                double _muSouth = mu.value[slaveCell.index()];
                double _muNorth = mu.value[masterCell.index()];
                double _muWest = vertMu.value[startVert.index()];
                double _muEast = vertMu.value[endVert.index()];
                double _mu = 2.0 * _muNorth * _muSouth / (_muNorth + _muSouth);
                double _rhoSouth = rho.value[slaveCell.index()];
                double _rhoNorth = rho.value[masterCell.index()];
                double _rhoWest = vertRho.value[startVert.index()];
                double _rhoEast = vertRho.value[endVert.index()];
                const auto metrics = HydroNativeMomentumDetail::axialFluxMetrics<Axisymmetric>(hEdge, masterCell, startVert, endVert);
                double _uWestMaster = 0.0;
                double _uWestSlave = 0.0;
                double _uEastMaster = 0.0;
                double _uEastSlave = 0.0;
                double _uWestMasterPre = 0.0;
                double _uWestSlavePre = 0.0;
                double _uEastMasterPre = 0.0;
                double _uEastSlavePre = 0.0;
                double _vSouthUpwind = 0.0;
                double _vSouthUpwindPre = 0.0;
                double _diffSouth = 0.0;
                double _vSouth = 0.0;
                double _vSouthPre = 0.0;
                if (slaveCell.supCellIndex == -1)
                {
                    Edge& southHEdge00 = mesh.hEdges(slaveCell.cellToHEdge[0][0]);
                    _vSouth = uy.value[southHEdge00.index()];
                    _vSouthPre = uyPre.value[southHEdge00.index()];
                    Edge& westVEdge00 = mesh.vEdges(slaveCell.cellToVEdge[0][0]);
                    Edge& eastVEdge10 = mesh.vEdges(slaveCell.cellToVEdge[1][0]);
                    _uWestSlave = ux.value[westVEdge00.index()];
                    _uWestSlavePre = uxPre.value[westVEdge00.index()];
                    _uEastSlave = ux.value[eastVEdge10.index()];
                    _uEastSlavePre = uxPre.value[eastVEdge10.index()];
                }
                else
                {
                    int _supCellIndex = slaveCell.supCellIndex;
                    for (int j = 0; j < supCells_[_supCellIndex].hEdges[0].ind.size(); j++)
                    {
                        _vSouth += uy.value[supCells_[_supCellIndex].hEdges[0].ind[j]] * supCells_[_supCellIndex].hEdges[0].coe[j];
                        _vSouthPre += uyPre.value[supCells_[_supCellIndex].hEdges[0].ind[j]] * supCells_[_supCellIndex].hEdges[0].coe[j];
                    }
                    for (int j = 0; j < supCells_[_supCellIndex].vEdges[0].ind.size(); j++)
                    {
                        _uWestSlave += ux.value[supCells_[_supCellIndex].vEdges[0].ind[j]] * supCells_[_supCellIndex].vEdges[0].coe[j];
                        _uWestSlavePre += uxPre.value[supCells_[_supCellIndex].vEdges[0].ind[j]] * supCells_[_supCellIndex].vEdges[0].coe[j];
                    }
                    for (int j = 0; j < supCells_[_supCellIndex].vEdges[1].ind.size(); j++)
                    {
                        _uEastSlave += ux.value[supCells_[_supCellIndex].vEdges[1].ind[j]] * supCells_[_supCellIndex].vEdges[1].coe[j];
                        _uEastSlavePre += uxPre.value[supCells_[_supCellIndex].vEdges[1].ind[j]] * supCells_[_supCellIndex].vEdges[1].coe[j];
                    }
                }
                double _qSouth = _rhoSouth * 0.5 * (uy.value[hEdge.index()] + _vSouth) * metrics.axialArea();
                double _qSouthPre = _rhoSouth * 0.5 * (uyPre.value[hEdge.index()] + _vSouthPre) * metrics.axialArea();
                {
                    _vSouthUpwind = _qSouth > 0.0 ? _vSouth : uy.value[hEdge.index()];
                }
                {
                    _vSouthUpwindPre = _qSouthPre > 0.0 ? _vSouthPre : uyPre.value[hEdge.index()];
                }
                _diffSouth = 2.0 * _muSouth * (uy.value[hEdge.index()] - _vSouth) / metrics.dy;
                double _vNorthUpwind = 0.0;
                double _vNorthUpwindPre = 0.0;
                double _diffNorth = 0.0;
                double _vNorth = 0.0;
                double _vNorthPre = 0.0;
                if (masterCell.supCellIndex == -1)
                {
                    Edge& northHEdge10 = mesh.hEdges(masterCell.cellToHEdge[1][0]);
                    _vNorth = uy.value[northHEdge10.index()];
                    _vNorthPre = uyPre.value[northHEdge10.index()];
                    Edge& westVEdge00 = mesh.vEdges(masterCell.cellToVEdge[0][0]);
                    Edge& eastVEdge10 = mesh.vEdges(masterCell.cellToVEdge[1][0]);
                    _uWestMaster = ux.value[westVEdge00.index()];
                    _uWestMasterPre = uxPre.value[westVEdge00.index()];
                    _uEastMaster = ux.value[eastVEdge10.index()];
                    _uEastMasterPre = uxPre.value[eastVEdge10.index()];
                }
                else
                {
                    int _supCellIndex = masterCell.supCellIndex;
                    for (int j = 0; j < supCells_[_supCellIndex].hEdges[1].ind.size(); j++)
                    {
                        _vNorth += uy.value[supCells_[_supCellIndex].hEdges[1].ind[j]] * supCells_[_supCellIndex].hEdges[1].coe[j];
                        _vNorthPre += uyPre.value[supCells_[_supCellIndex].hEdges[1].ind[j]] * supCells_[_supCellIndex].hEdges[1].coe[j];
                    }
                    for (int j = 0; j < supCells_[_supCellIndex].vEdges[0].ind.size(); j++)
                    {
                        _uWestMaster += ux.value[supCells_[_supCellIndex].vEdges[0].ind[j]] * supCells_[_supCellIndex].vEdges[0].coe[j];
                        _uWestMasterPre += uxPre.value[supCells_[_supCellIndex].vEdges[0].ind[j]] * supCells_[_supCellIndex].vEdges[0].coe[j];
                    }
                    for (int j = 0; j < supCells_[_supCellIndex].vEdges[1].ind.size(); j++)
                    {
                        _uEastMaster += ux.value[supCells_[_supCellIndex].vEdges[1].ind[j]] * supCells_[_supCellIndex].vEdges[1].coe[j];
                        _uEastMasterPre += uxPre.value[supCells_[_supCellIndex].vEdges[1].ind[j]] * supCells_[_supCellIndex].vEdges[1].coe[j];
                    }
                }
                double _qNorth = _rhoNorth * 0.5 * (uy.value[hEdge.index()] + _vNorth) * metrics.axialArea();
                double _qNorthPre = _rhoNorth * 0.5 * (uyPre.value[hEdge.index()] + _vNorthPre) * metrics.axialArea();
                {
                    _vNorthUpwind = _qNorth > 0.0 ? uy.value[hEdge.index()] : _vNorth;
                }
                {
                    _vNorthUpwindPre = _qNorthPre > 0.0 ? uyPre.value[hEdge.index()] : _vNorthPre;
                }
                _diffNorth = 2.0 * _muNorth * (_vNorth - uy.value[hEdge.index()]) / metrics.dy;
                double _uWest = 0.5 * (_uWestMaster + _uWestSlave);
                double _uWestPre = 0.5 * (_uWestMasterPre + _uWestSlavePre);
                double _qWest = _rhoWest * _uWest * metrics.westArea();
                double _qWestPre = _rhoWest * _uWestPre * metrics.westArea();
                double _vWestUpwind = 0.0;
                double _vWestUpwindPre = 0.0;
                double _diffWest = 0.0;
                double _vWest = 0.0;
                double _vWestPre = 0.0;
                if (startVert.vertToHEdge[0] != -1)
                {
                    Edge& westHEdge = mesh.hEdges(startVert.vertToHEdge[0]);
                    _vWest = uy.value[westHEdge.index()];
                    _vWestPre = uyPre.value[westHEdge.index()];
                    _diffWest = _muWest * ((uy.value[hEdge.index()] - _vWest) / (hEdge.midpointCoord().x - westHEdge.midpointCoord().x) + (_uWestMaster - _uWestSlave) / metrics.dy);
                }
                else if (startVert.supCell != -1)
                {
                    Cell& supCell = mesh.cells(startVert.supCell);
                    int supSonCellIndex = supCell.supCellIndex * 4 + 1;
                    for (int j = 0; j < supSonCells_[supSonCellIndex].hEdges[1].ind.size(); j++)
                    {
                        _vWest += uy.value[supSonCells_[supSonCellIndex].hEdges[1].ind[j]] * supSonCells_[supSonCellIndex].hEdges[1].coe[j];
                        _vWestPre += uyPre.value[supSonCells_[supSonCellIndex].hEdges[1].ind[j]] * supSonCells_[supSonCellIndex].hEdges[1].coe[j];
                    }
                    _diffWest = _muWest * ((uy.value[hEdge.index()] - _vWest) / metrics.dx + (_uWestMaster - _uWestSlave) / metrics.dy);
                }
                else
                {
                    _diffWest = _muWest * (_uWestMaster - _uWestSlave) / metrics.dy;
                }
                {
                    _vWestUpwind = _qWest > 0.0 ? _vWest : uy.value[hEdge.index()];
                }
                {
                    _vWestUpwindPre = _qWestPre > 0.0 ? _vWestPre : uyPre.value[hEdge.index()];
                }
                double _uEast = 0.5 * (_uEastMaster + _uEastSlave);
                double _uEastPre = 0.5 * (_uEastMasterPre + _uEastSlavePre);
                double _qEast = _rhoEast * _uEast * metrics.eastArea();
                double _qEastPre = _rhoEast * _uEastPre * metrics.eastArea();
                double _vEastUpwind = 0.0;
                double _vEastUpwindPre = 0.0;
                double _diffEast = 0.0;
                double _vEast = 0.0;
                double _vEastPre = 0.0;
                if (endVert.vertToHEdge[1] != -1)
                {
                    Edge& eastHEdge = mesh.hEdges(endVert.vertToHEdge[1]);
                    _vEast = uy.value[eastHEdge.index()];
                    _vEastPre = uyPre.value[eastHEdge.index()];
                    _diffEast = _muEast * ((_vEast - uy.value[hEdge.index()]) / (eastHEdge.midpointCoord().x - hEdge.midpointCoord().x) + (_uEastMaster - _uEastSlave) / metrics.dy);
                }
                else if (endVert.supCell != -1)
                {
                    Cell& supCell = mesh.cells(endVert.supCell);
                    int supSonCellIndex = supCell.supCellIndex * 4;
                    for (int j = 0; j < supSonCells_[supSonCellIndex].hEdges[1].ind.size(); j++)
                    {
                        _vEast += uy.value[supSonCells_[supSonCellIndex].hEdges[1].ind[j]] * supSonCells_[supSonCellIndex].hEdges[1].coe[j];
                        _vEastPre += uyPre.value[supSonCells_[supSonCellIndex].hEdges[1].ind[j]] * supSonCells_[supSonCellIndex].hEdges[1].coe[j];
                    }
                    _diffEast = _muEast * ((_vEast - uy.value[hEdge.index()]) / metrics.dx + (_uEastMaster - _uEastSlave) / metrics.dy);
                }
                else
                {
                    _diffEast = _muEast * (_uEastMaster - _uEastSlave) / metrics.dy;
                }
                {
                    _vEastUpwind = _qEast > 0.0 ? uy.value[hEdge.index()] : _vEast;
                }
                {
                    _vEastUpwindPre = _qEastPre > 0.0 ? uyPre.value[hEdge.index()] : _vEastPre;
                }
                double _n = _qEast * _vEastUpwind - _qWest * _vWestUpwind + _qNorth * _vNorthUpwind - _qSouth * _vSouthUpwind;
                double _nPre = _qEastPre * _vEastUpwindPre - _qWestPre * _vWestUpwindPre + _qNorthPre * _vNorthUpwindPre - _qSouthPre * _vSouthUpwindPre;
                if(!includeAdvection) _n=_nPre=0;
                double _nDiff = metrics.diffusionNet(_diffEast, _diffWest, _diffNorth, _diffSouth);
                assembly.rhs[row] = _dA * _rho / dt * uy.value[hEdge.index()] - 1.5 * _n + 0.5 * (_nPre + _nDiff) + _dA * 0.5 * (nativeBodyY.value[masterCell.index()] + nativeBodyY.value[slaveCell.index()]);
            }
            else if (masterCell.level() < slaveCell.level())
            {
                double _dA = HydroNativeMomentumDetail::volume<Axisymmetric>(hEdge, slaveCell);
                double _rhoMaster = 0.0;
                int supSonCellIndex = masterCell.centerCoord().x < slaveCell.centerCoord().x ? masterCell.supCellIndex * 4 + 1 : masterCell.supCellIndex * 4;
                for (int j = 0; j < supSonCells_[supSonCellIndex].center.ind.size(); j++)
                {
                    _rhoMaster += rho.value[supSonCells_[supSonCellIndex].center.ind[j]] * supSonCells_[supSonCellIndex].center.coe[j];
                }
                double _rhoSlave = rho.value[slaveCell.index()];
                double _rho = 0.5 * (_rhoMaster + _rhoSlave);
                double _muSouth = mu.value[slaveCell.index()];
                double _muNorth = 0.0;
                for (int j = 0; j < supSonCells_[supSonCellIndex].center.ind.size(); j++)
                {
                    _muNorth += mu.value[supSonCells_[supSonCellIndex].center.ind[j]] * supSonCells_[supSonCellIndex].center.coe[j];
                }
                double _muWest = vertMu.value[startVert.index()];
                double _muEast = vertMu.value[endVert.index()];
                double _rhoSouth = rho.value[slaveCell.index()];
                double _rhoNorth = 0.0;
                for (int j = 0; j < supSonCells_[supSonCellIndex].center.ind.size(); j++)
                {
                    _rhoNorth += rho.value[supSonCells_[supSonCellIndex].center.ind[j]] * supSonCells_[supSonCellIndex].center.coe[j];
                }
                double _rhoWest = vertRho.value[startVert.index()];
                double _rhoEast = vertRho.value[endVert.index()];
                const auto metrics = HydroNativeMomentumDetail::axialFluxMetrics<Axisymmetric>(hEdge, slaveCell, startVert, endVert);
                double _uWestMaster = 0.0;
                double _uWestSlave = 0.0;
                double _uEastMaster = 0.0;
                double _uEastSlave = 0.0;
                double _uWestMasterPre = 0.0;
                double _uWestSlavePre = 0.0;
                double _uEastMasterPre = 0.0;
                double _uEastSlavePre = 0.0;
                double _vSouthUpwind = 0.0;
                double _vSouthUpwindPre = 0.0;
                double _diffSouth = 0.0;
                double _vSouth = 0.0;
                double _vSouthPre = 0.0;
                {
                    Edge& southHEdge00 = mesh.hEdges(slaveCell.cellToHEdge[0][0]);
                    _vSouth = uy.value[southHEdge00.index()];
                    _vSouthPre = uyPre.value[southHEdge00.index()];
                    Edge& westVEdge00 = mesh.vEdges(slaveCell.cellToVEdge[0][0]);
                    Edge& eastVEdge10 = mesh.vEdges(slaveCell.cellToVEdge[1][0]);
                    _uWestSlave = ux.value[westVEdge00.index()];
                    _uWestSlavePre = uxPre.value[westVEdge00.index()];
                    _uEastSlave = ux.value[eastVEdge10.index()];
                    _uEastSlavePre = uxPre.value[eastVEdge10.index()];
                }
                double _qSouth = _rhoSouth * 0.5 * (uy.value[hEdge.index()] + _vSouth) * metrics.axialArea();
                double _qSouthPre = _rhoSouth * 0.5 * (uyPre.value[hEdge.index()] + _vSouthPre) * metrics.axialArea();
                {
                    _vSouthUpwind = _qSouth > 0.0 ? _vSouth : uy.value[hEdge.index()];
                }
                {
                    _vSouthUpwindPre = _qSouthPre > 0.0 ? _vSouthPre : uyPre.value[hEdge.index()];
                }
                _diffSouth = 2.0 * _muSouth * (uy.value[hEdge.index()] - _vSouth) / metrics.dy;
                double _vNorthUpwind = 0.0;
                double _vNorthUpwindPre = 0.0;
                double _diffNorth = 0.0;
                double _vNorth = 0.0;
                double _vNorthPre = 0.0;
                {
                    for (int j = 0; j < supSonCells_[supSonCellIndex].hEdges[1].ind.size(); j++)
                    {
                        _vNorth += uy.value[supSonCells_[supSonCellIndex].hEdges[1].ind[j]] * supSonCells_[supSonCellIndex].hEdges[1].coe[j];
                        _vNorthPre += uyPre.value[supSonCells_[supSonCellIndex].hEdges[1].ind[j]] * supSonCells_[supSonCellIndex].hEdges[1].coe[j];
                    }
                    for (int j = 0; j < supSonCells_[supSonCellIndex].vEdges[0].ind.size(); j++)
                    {
                        _uWestMaster += ux.value[supSonCells_[supSonCellIndex].vEdges[0].ind[j]] * supSonCells_[supSonCellIndex].vEdges[0].coe[j];
                        _uWestMasterPre += uxPre.value[supSonCells_[supSonCellIndex].vEdges[0].ind[j]] * supSonCells_[supSonCellIndex].vEdges[0].coe[j];
                    }
                    for (int j = 0; j < supSonCells_[supSonCellIndex].vEdges[1].ind.size(); j++)
                    {
                        _uEastMaster += ux.value[supSonCells_[supSonCellIndex].vEdges[1].ind[j]] * supSonCells_[supSonCellIndex].vEdges[1].coe[j];
                        _uEastMasterPre += uxPre.value[supSonCells_[supSonCellIndex].vEdges[1].ind[j]] * supSonCells_[supSonCellIndex].vEdges[1].coe[j];
                    }
                }
                double _qNorth = _rhoNorth * 0.5 * (uy.value[hEdge.index()] + _vNorth) * metrics.axialArea();
                double _qNorthPre = _rhoNorth * 0.5 * (uyPre.value[hEdge.index()] + _vNorthPre) * metrics.axialArea();
                {
                    _vNorthUpwind = _qNorth > 0.0 ? uy.value[hEdge.index()] : _vNorth;
                }
                {
                    _vNorthUpwindPre = _qNorthPre > 0.0 ? uyPre.value[hEdge.index()] : _vNorthPre;
                }
                _diffNorth = 2.0 * _muNorth * (_vNorth - uy.value[hEdge.index()]) / metrics.dy;
                double _uWest = 0.5 * (_uWestMaster + _uWestSlave);
                double _uWestPre = 0.5 * (_uWestMasterPre + _uWestSlavePre);
                double _qWest = _rhoWest * _uWest * metrics.westArea();
                double _qWestPre = _rhoWest * _uWestPre * metrics.westArea();
                double _vWestUpwind = 0.0;
                double _vWestUpwindPre = 0.0;
                double _diffWest = 0.0;
                double _vWest = 0.0;
                double _vWestPre = 0.0;
                if (startVert.vertToHEdge[0] != -1)
                {
                    Edge& westHEdge = mesh.hEdges(startVert.vertToHEdge[0]);
                    _vWest = uy.value[westHEdge.index()];
                    _vWestPre = uyPre.value[westHEdge.index()];
                    _diffWest = _muWest * ((uy.value[hEdge.index()] - _vWest) / (hEdge.midpointCoord().x - westHEdge.midpointCoord().x) + (_uWestMaster - _uWestSlave) / metrics.dy);
                }
                else
                {
                    _diffWest = _muWest * (_uWestMaster - _uWestSlave) / metrics.dy;
                }
                {
                    _vWestUpwind = _qWest > 0.0 ? _vWest : uy.value[hEdge.index()];
                }
                {
                    _vWestUpwindPre = _qWestPre > 0.0 ? _vWestPre : uyPre.value[hEdge.index()];
                }
                double _uEast = 0.5 * (_uEastMaster + _uEastSlave);
                double _uEastPre = 0.5 * (_uEastMasterPre + _uEastSlavePre);
                double _qEast = _rhoEast * _uEast * metrics.eastArea();
                double _qEastPre = _rhoEast * _uEastPre * metrics.eastArea();
                double _vEastUpwind = 0.0;
                double _vEastUpwindPre = 0.0;
                double _diffEast = 0.0;
                double _vEast = 0.0;
                double _vEastPre = 0.0;
                if (endVert.vertToHEdge[1] != -1)
                {
                    Edge& eastHEdge = mesh.hEdges(endVert.vertToHEdge[1]);
                    _vEast = uy.value[eastHEdge.index()];
                    _vEastPre = uyPre.value[eastHEdge.index()];
                    _diffEast = _muEast * ((_vEast - uy.value[hEdge.index()]) / (eastHEdge.midpointCoord().x - hEdge.midpointCoord().x) + (_uEastMaster - _uEastSlave) / metrics.dy);
                }
                else
                {
                    _diffEast = _muEast * (_uEastMaster - _uEastSlave) / metrics.dy;
                }
                {
                    _vEastUpwind = _qEast > 0.0 ? uy.value[hEdge.index()] : _vEast;
                }
                {
                    _vEastUpwindPre = _qEastPre > 0.0 ? uyPre.value[hEdge.index()] : _vEastPre;
                }
                double _n = _qEast * _vEastUpwind - _qWest * _vWestUpwind + _qNorth * _vNorthUpwind - _qSouth * _vSouthUpwind;
                double _nPre = _qEastPre * _vEastUpwindPre - _qWestPre * _vWestUpwindPre + _qNorthPre * _vNorthUpwindPre - _qSouthPre * _vSouthUpwindPre;
                if(!includeAdvection) _n=_nPre=0;
                double _nDiff = metrics.diffusionNet(_diffEast, _diffWest, _diffNorth, _diffSouth);
                double _FyMaster = 0.0;
                for (int j = 0; j < supSonCells_[supSonCellIndex].center.ind.size(); j++)
                {
                    _FyMaster += nativeBodyY.value[supSonCells_[supSonCellIndex].center.ind[j]] * supSonCells_[supSonCellIndex].center.coe[j];
                }
                assembly.rhs[row] = _dA * _rho / dt * uy.value[hEdge.index()] - 1.5 * _n + 0.5 * (_nPre + _nDiff) + _dA * 0.5 * (_FyMaster + nativeBodyY.value[slaveCell.index()]);
            }
            else
            {
                double _dA = HydroNativeMomentumDetail::volume<Axisymmetric>(hEdge, masterCell);
                double _rhoMaster = rho.value[masterCell.index()];
                double _rhoSlave = 0.0;
                int supSonCellIndex = slaveCell.centerCoord().x < masterCell.centerCoord().x ? slaveCell.supCellIndex * 4 + 2 : slaveCell.supCellIndex * 4 + 3;
                for (int j = 0; j < supSonCells_[supSonCellIndex].center.ind.size(); j++)
                {
                    _rhoSlave += rho.value[supSonCells_[supSonCellIndex].center.ind[j]] * supSonCells_[supSonCellIndex].center.coe[j];
                }
                double _rho = 0.5 * (_rhoMaster + _rhoSlave);
                double _muNorth = mu.value[masterCell.index()];
                double _muSouth = 0.0;
                for (int j = 0; j < supSonCells_[supSonCellIndex].center.ind.size(); j++)
                {
                    _muSouth += mu.value[supSonCells_[supSonCellIndex].center.ind[j]] * supSonCells_[supSonCellIndex].center.coe[j];
                }
                double _muWest = vertMu.value[startVert.index()];
                double _muEast = vertMu.value[endVert.index()];
                double _rhoNorth = rho.value[masterCell.index()];
                double _rhoSouth = 0.0;
                for (int j = 0; j < supSonCells_[supSonCellIndex].center.ind.size(); j++)
                {
                    _rhoSouth += rho.value[supSonCells_[supSonCellIndex].center.ind[j]] * supSonCells_[supSonCellIndex].center.coe[j];
                }
                double _rhoWest = vertRho.value[startVert.index()];
                double _rhoEast = vertRho.value[endVert.index()];
                const auto metrics = HydroNativeMomentumDetail::axialFluxMetrics<Axisymmetric>(hEdge, masterCell, startVert, endVert);
                double _uWestMaster = 0.0;
                double _uWestSlave = 0.0;
                double _uEastMaster = 0.0;
                double _uEastSlave = 0.0;
                double _uWestMasterPre = 0.0;
                double _uWestSlavePre = 0.0;
                double _uEastMasterPre = 0.0;
                double _uEastSlavePre = 0.0;
                double _vSouthUpwind = 0.0;
                double _vSouthUpwindPre = 0.0;
                double _diffSouth = 0.0;
                double _vSouth = 0.0;
                double _vSouthPre = 0.0;
                {
                    for (int j = 0; j < supSonCells_[supSonCellIndex].hEdges[0].ind.size(); j++)
                    {
                        _vSouth += uy.value[supSonCells_[supSonCellIndex].hEdges[0].ind[j]] * supSonCells_[supSonCellIndex].hEdges[0].coe[j];
                        _vSouthPre += uyPre.value[supSonCells_[supSonCellIndex].hEdges[0].ind[j]] * supSonCells_[supSonCellIndex].hEdges[0].coe[j];
                    }
                    for (int j = 0; j < supSonCells_[supSonCellIndex].vEdges[0].ind.size(); j++)
                    {
                        _uWestSlave += ux.value[supSonCells_[supSonCellIndex].vEdges[0].ind[j]] * supSonCells_[supSonCellIndex].vEdges[0].coe[j];
                        _uWestSlavePre += uxPre.value[supSonCells_[supSonCellIndex].vEdges[0].ind[j]] * supSonCells_[supSonCellIndex].vEdges[0].coe[j];
                    }
                    for (int j = 0; j < supSonCells_[supSonCellIndex].vEdges[1].ind.size(); j++)
                    {
                        _uEastSlave += ux.value[supSonCells_[supSonCellIndex].vEdges[1].ind[j]] * supSonCells_[supSonCellIndex].vEdges[1].coe[j];
                        _uEastSlavePre += uxPre.value[supSonCells_[supSonCellIndex].vEdges[1].ind[j]] * supSonCells_[supSonCellIndex].vEdges[1].coe[j];
                    }
                }
                double _qSouth = _rhoSouth * 0.5 * (uy.value[hEdge.index()] + _vSouth) * metrics.axialArea();
                double _qSouthPre = _rhoSouth * 0.5 * (uyPre.value[hEdge.index()] + _vSouthPre) * metrics.axialArea();
                {
                    _vSouthUpwind = _qSouth > 0.0 ? _vSouth : uy.value[hEdge.index()];
                }
                {
                    _vSouthUpwindPre = _qSouthPre > 0.0 ? _vSouthPre : uyPre.value[hEdge.index()];
                }
                _diffSouth = 2.0 * _muSouth * (uy.value[hEdge.index()] - _vSouth) / metrics.dy;
                double _vNorthUpwind = 0.0;
                double _vNorthUpwindPre = 0.0;
                double _diffNorth = 0.0;
                double _vNorth = 0.0;
                double _vNorthPre = 0.0;
                {
                    Edge& northHEdge10 = mesh.hEdges(masterCell.cellToHEdge[1][0]);
                    _vNorth = uy.value[northHEdge10.index()];
                    _vNorthPre = uyPre.value[northHEdge10.index()];
                    Edge& westVEdge00 = mesh.vEdges(masterCell.cellToVEdge[0][0]);
                    Edge& eastVEdge10 = mesh.vEdges(masterCell.cellToVEdge[1][0]);
                    _uWestMaster = ux.value[westVEdge00.index()];
                    _uWestMasterPre = uxPre.value[westVEdge00.index()];
                    _uEastMaster = ux.value[eastVEdge10.index()];
                    _uEastMasterPre = uxPre.value[eastVEdge10.index()];
                }
                double _qNorth = _rhoNorth * 0.5 * (uy.value[hEdge.index()] + _vNorth) * metrics.axialArea();
                double _qNorthPre = _rhoNorth * 0.5 * (uyPre.value[hEdge.index()] + _vNorthPre) * metrics.axialArea();
                {
                    _vNorthUpwind = _qNorth > 0.0 ? uy.value[hEdge.index()] : _vNorth;
                }
                {
                    _vNorthUpwindPre = _qNorthPre > 0.0 ? uyPre.value[hEdge.index()] : _vNorthPre;
                }
                _diffNorth = 2.0 * _muNorth * (_vNorth - uy.value[hEdge.index()]) / metrics.dy;
                double _uWest = 0.5 * (_uWestMaster + _uWestSlave);
                double _uWestPre = 0.5 * (_uWestMasterPre + _uWestSlavePre);
                double _qWest = _rhoWest * _uWest * metrics.westArea();
                double _qWestPre = _rhoWest * _uWestPre * metrics.westArea();
                double _vWestUpwind = 0.0;
                double _vWestUpwindPre = 0.0;
                double _diffWest = 0.0;
                double _vWest = 0.0;
                double _vWestPre = 0.0;
                if (startVert.vertToHEdge[0] != -1)
                {
                    Edge& westHEdge = mesh.hEdges(startVert.vertToHEdge[0]);
                    _vWest = uy.value[westHEdge.index()];
                    _vWestPre = uyPre.value[westHEdge.index()];
                    _diffWest = _muWest * ((uy.value[hEdge.index()] - _vWest) / (hEdge.midpointCoord().x - westHEdge.midpointCoord().x) + (_uWestMaster - _uWestSlave) / metrics.dy);
                }
                else
                {
                    _diffWest = _muWest * (_uWestMaster - _uWestSlave) / metrics.dy;
                }
                {
                    _vWestUpwind = _qWest > 0.0 ? _vWest : uy.value[hEdge.index()];
                }
                {
                    _vWestUpwindPre = _qWestPre > 0.0 ? _vWestPre : uyPre.value[hEdge.index()];
                }
                double _uEast = 0.5 * (_uEastMaster + _uEastSlave);
                double _uEastPre = 0.5 * (_uEastMasterPre + _uEastSlavePre);
                double _qEast = _rhoEast * _uEast * metrics.eastArea();
                double _qEastPre = _rhoEast * _uEastPre * metrics.eastArea();
                double _vEastUpwind = 0.0;
                double _vEastUpwindPre = 0.0;
                double _diffEast = 0.0;
                double _vEast = 0.0;
                double _vEastPre = 0.0;
                if (endVert.vertToHEdge[1] != -1)
                {
                    Edge& eastHEdge = mesh.hEdges(endVert.vertToHEdge[1]);
                    _vEast = uy.value[eastHEdge.index()];
                    _vEastPre = uyPre.value[eastHEdge.index()];
                    _diffEast = _muEast * ((_vEast - uy.value[hEdge.index()]) / (eastHEdge.midpointCoord().x - hEdge.midpointCoord().x) + (_uEastMaster - _uEastSlave) / metrics.dy);
                }
                else
                {
                    _diffEast = _muEast * (_uEastMaster - _uEastSlave) / metrics.dy;
                }
                {
                    _vEastUpwind = _qEast > 0.0 ? uy.value[hEdge.index()] : _vEast;
                }
                {
                    _vEastUpwindPre = _qEastPre > 0.0 ? uyPre.value[hEdge.index()] : _vEastPre;
                }
                double _n = _qEast * _vEastUpwind - _qWest * _vWestUpwind + _qNorth * _vNorthUpwind - _qSouth * _vSouthUpwind;
                double _nPre = _qEastPre * _vEastUpwindPre - _qWestPre * _vWestUpwindPre + _qNorthPre * _vNorthUpwindPre - _qSouthPre * _vSouthUpwindPre;
                if(!includeAdvection) _n=_nPre=0;
                double _nDiff = metrics.diffusionNet(_diffEast, _diffWest, _diffNorth, _diffSouth);
                double _FySlave = 0.0;
                for (int j = 0; j < supSonCells_[supSonCellIndex].center.ind.size(); j++)
                {
                    _FySlave += nativeBodyY.value[supSonCells_[supSonCellIndex].center.ind[j]] * supSonCells_[supSonCellIndex].center.coe[j];
                }
                assembly.rhs[row] = _dA * _rho / dt * uy.value[hEdge.index()] - 1.5 * _n + 0.5 * (_nPre + _nDiff) + _dA * 0.5 * (nativeBodyY.value[masterCell.index()] + _FySlave);
            }
        }
        else if (hEdge.masterCell == -1)
        {
            if constexpr (Axisymmetric)
            {
                assembly.rhs[row] = uy.value[i];
            }
            else
            {
                assembly.rhs[row] = 0.0;
            }
        }
        else
        {
            if constexpr (Axisymmetric)
            {
                assembly.rhs[row] = uy.value[i];
            }
            else
            {
                assembly.rhs[row] = 0.0;
            }
        }
    }
}
