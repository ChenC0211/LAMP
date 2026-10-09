#pragma once

// The pressure-face correction is identical in XY and RZ. Preserve the
// original diagnostic boundary policy and coordinate-specific divergence.
template<bool Axisymmetric, typename Mesh>
void Hydrodynamics::correctNativeVelocity(Mesh& mesh, const double& dt)
{
#pragma omp parallel for num_threads(thread_num)
    for (int i = 0; i < mesh.vEdgeIndexNum(); i++)
    {
        Edge& vEdge = mesh.vEdges(mesh.vEdgeIndexToLabel[i]);
        if (vEdge.masterCell != -1 && vEdge.slaveCell != -1)
        {
            const auto face = pressureFaceData(mesh, vEdge, false, pressurePotential_);
            ux.value[i] = uxMid.value[i] - dt / face.density * face.gradient;
        }
        else if (vEdge.masterCell != -1)
        {
            ux.value[i] = 0.0;
        }
        else
        {
            ux.value[i] = 0.0;
        }
    }

#pragma omp parallel for num_threads(thread_num)
    for (int i = 0; i < mesh.hEdgeIndexNum(); i++)
    {
        Edge& hEdge = mesh.hEdges(mesh.hEdgeIndexToLabel[i]);
        if (hEdge.masterCell != -1 && hEdge.slaveCell != -1)
        {
            const auto face = pressureFaceData(mesh, hEdge, true, pressurePotential_);
            uy.value[i] = uyMid.value[i] - dt / face.density * face.gradient;
        }
        else if (hEdge.masterCell != -1)
        {
            uy.value[i] = 0.0;
        }
        else
        {
            uy.value[i] = 0.0;
        }
    }

    uMag.initialize(mesh, 0.0);
    divU.initialize(mesh, 0.0);
    dxUx.initialize(mesh, 0.0);
    dyUx.initialize(mesh, 0.0);
    dxUy.initialize(mesh, 0.0);
    dyUy.initialize(mesh, 0.0);

    VertField<double> vertUx = std::move(*ux.toVertField(mesh));
    VertField<double> vertUy = std::move(*uy.toVertField(mesh));
#pragma omp parallel for num_threads(thread_num)
    for (int i = 0; i < mesh.cellIndexNum(); i++)
    {
        Cell& cell = mesh.cells(mesh.cellIndexToLabel[i]);
        if constexpr (!Axisymmetric)
        {
            if (cell.boundary() != Boundary::non)
            {
                continue;
            }
        }
        double _dx = cell.size().x;
        double _dy = cell.size().y;
        double _uWest = 0.0;
        double _vWest = 0.0;
        if (cell.cellToVEdge[0][1] != -1)
        {
            _uWest = vertUx.value[mesh.verts(mesh.vEdges(cell.cellToVEdge[0][0]).edgeToVert[1]).index()];
            _vWest = vertUy.value[mesh.verts(mesh.vEdges(cell.cellToVEdge[0][0]).edgeToVert[1]).index()];
        }
        else
        {
            _uWest = ux.value[mesh.vEdges(cell.cellToVEdge[0][0]).index()];
            _vWest = 0.5 * (vertUy.value[mesh.verts(mesh.vEdges(cell.cellToVEdge[0][0]).edgeToVert[0]).index()] + vertUy.value[mesh.verts(mesh.vEdges(cell.cellToVEdge[0][0]).edgeToVert[1]).index()]);
        }
        double _uEast = 0.0;
        double _vEast = 0.0;
        if (cell.cellToVEdge[1][1] != -1)
        {
            _uEast = vertUx.value[mesh.verts(mesh.vEdges(cell.cellToVEdge[1][0]).edgeToVert[1]).index()];
            _vEast = vertUy.value[mesh.verts(mesh.vEdges(cell.cellToVEdge[1][0]).edgeToVert[1]).index()];
        }
        else
        {
            _uEast = ux.value[mesh.vEdges(cell.cellToVEdge[1][0]).index()];
            _vEast = 0.5 * (vertUy.value[mesh.verts(mesh.vEdges(cell.cellToVEdge[1][0]).edgeToVert[0]).index()] + vertUy.value[mesh.verts(mesh.vEdges(cell.cellToVEdge[1][0]).edgeToVert[1]).index()]);
        }
        double _vSouth = 0.0;
        double _uSouth = 0.0;
        if (cell.cellToHEdge[0][1] != -1)
        {
            _vSouth = vertUy.value[mesh.verts(mesh.hEdges(cell.cellToHEdge[0][0]).edgeToVert[1]).index()];
            _uSouth = vertUx.value[mesh.verts(mesh.hEdges(cell.cellToHEdge[0][0]).edgeToVert[1]).index()];
        }
        else
        {
            _vSouth = uy.value[mesh.hEdges(cell.cellToHEdge[0][0]).index()];
            _uSouth = 0.5 * (vertUx.value[mesh.verts(mesh.hEdges(cell.cellToHEdge[0][0]).edgeToVert[0]).index()] + vertUx.value[mesh.verts(mesh.hEdges(cell.cellToHEdge[0][0]).edgeToVert[1]).index()]);
        }
        double _vNorth = 0.0;
        double _uNorth = 0.0;
        if (cell.cellToHEdge[1][1] != -1)
        {
            _vNorth = vertUy.value[mesh.verts(mesh.hEdges(cell.cellToHEdge[1][0]).edgeToVert[1]).index()];
            _uNorth = vertUx.value[mesh.verts(mesh.hEdges(cell.cellToHEdge[1][0]).edgeToVert[1]).index()];
        }
        else
        {
            _vNorth = uy.value[mesh.hEdges(cell.cellToHEdge[1][0]).index()];
            _uNorth = 0.5 * (vertUx.value[mesh.verts(mesh.hEdges(cell.cellToHEdge[1][0]).edgeToVert[0]).index()] + vertUx.value[mesh.verts(mesh.hEdges(cell.cellToHEdge[1][0]).edgeToVert[1]).index()]);
        }
        double _ux = 0.25 * (_uWest + _uEast + _uSouth + _uNorth);
        double _uy = 0.25 * (_vSouth + _vNorth + _vWest + _vEast);
        uMag.value[i] = sqrt(_ux * _ux + _uy * _uy);
        dxUx.value[i] = (_uEast - _uWest) / _dx;
        dyUx.value[i] = (_uNorth - _uSouth) / _dy;
        dxUy.value[i] = (_vEast - _vWest) / _dx;
        dyUy.value[i] = (_vNorth - _vSouth) / _dy;
        if constexpr (Axisymmetric)
        {
            double _dA = 2.0 * Pi * cell.centerCoord().x * cell.area();
            double AreaNS = 2.0 * Pi * cell.centerCoord().x * cell.size().x;
            double AreaE = 2.0 * Pi * mesh.verts(cell.cellToVert[1]).coord().x * cell.size().y;
            double AreaW = 2.0 * Pi * mesh.verts(cell.cellToVert[0]).coord().x * cell.size().y;
            divU.value[i] = (_uEast * AreaE - _uWest * AreaW + (_vNorth - _vSouth) * AreaNS) / _dA;
        }
        else
        {
            divU.value[i] = dxUx.value[i] + dyUy.value[i];
        }
        magGradU.value[i] = sqrt(dxUx.value[i] * dxUx.value[i] + dyUx.value[i] * dyUx.value[i] + dxUy.value[i] * dxUy.value[i] + dyUy.value[i] * dyUy.value[i]);
    }
}


namespace HydroNativePressureDetail
{
struct CartesianMetrics
{
    double _dx, _dy, _dxDivdy, _dyDivdx;
    double west() const { return _dyDivdx; }
    double east() const { return _dyDivdx; }
    double axial() const { return _dxDivdy; }
    double fluxRate(double dt, double _uEast, double _uWest, double _vNorth, double _vSouth) const
    {
        return 1.0 / dt * ((_uEast - _uWest) * _dy + (_vNorth - _vSouth) * _dx);
    }
};
struct AxisymmetricMetrics
{
    double AreaNS, AreaE, AreaW, AreaNS_divDy, AreaE_divDx, AreaW_divDx;
    double west() const { return AreaW_divDx; }
    double east() const { return AreaE_divDx; }
    double axial() const { return AreaNS_divDy; }
    double fluxRate(double dt, double _uEast, double _uWest, double _vNorth, double _vSouth) const
    {
        return 1.0 / dt * ((_uEast * AreaE - _uWest * AreaW) + (_vNorth - _vSouth) * AreaNS);
    }
};

// Keep the original operation order for both geometry strategies. In particular,
// XY uses the reciprocal aspect ratio; RZ keeps separate east/west flux products.
template<bool Axisymmetric, typename Mesh>
auto metrics(Mesh& mesh, Cell& cell)
{
    if constexpr (Axisymmetric)
    {
        double AreaNS = 2.0 * Pi * cell.centerCoord().x * cell.size().x;
        double AreaE = 2.0 * Pi * mesh.verts(cell.cellToVert[1]).coord().x * cell.size().y;
        double AreaW = 2.0 * Pi * mesh.verts(cell.cellToVert[0]).coord().x * cell.size().y;
        double AreaNS_divDy = AreaNS / cell.size().y;
        double AreaE_divDx = AreaE / cell.size().x;
        double AreaW_divDx = AreaW / cell.size().x;
        return AxisymmetricMetrics{AreaNS, AreaE, AreaW, AreaNS_divDy, AreaE_divDx, AreaW_divDx};
    }
    else
    {
        double _dx = cell.size().x;
        double _dy = cell.size().y;
        double _dxDivdy = _dx / _dy;
        double _dyDivdx = 1.0 / _dxDivdy;
        return CartesianMetrics{_dx, _dy, _dxDivdy, _dyDivdx};
    }
}
}

// Shared cell-row topology, interpolation, pressure gauge and triplet ordering.
// The conservative native-adaptive face assembly below remains a distinct policy.
template<bool Axisymmetric, typename Mesh>
HydroPressureSystem Hydrodynamics::assembleCellPressure(Mesh& mesh, const double& dt)
{
    vector<vector<Eigen::Triplet<double>>> _pressureCoeff(mesh.cellIndexNum());
    Eigen::VectorXd _pressureRHS(mesh.cellIndexNum());
#pragma omp parallel for num_threads(thread_num)
    for (int cellI = 0; cellI < mesh.cellIndexNum(); cellI++)
    {
        Cell& cell = mesh.cells(mesh.cellIndexToLabel[cellI]);

        if (activeOptions_.closedDomainProjection ? cellI == 0 : (
            cell.boundary() == Boundary::east || cell.boundary() == Boundary::northEast || cell.boundary() == Boundary::southEast
            // || 
            // cell.boundary() == Boundary::west || cell.boundary() == Boundary::northWest || cell.boundary() == Boundary::southWest
            // || 
            // cell.boundary() == Boundary::south || cell.boundary() == Boundary::southEast || cell.boundary() == Boundary::southWest
            // || 
            // cell.boundary() == Boundary::north || cell.boundary() == Boundary::northEast || cell.boundary() == Boundary::northWest
        ))
        {
            _pressureCoeff[cellI].push_back(Trip(cellI, cellI, 1.0));
            _pressureRHS[cellI] = 0.0;
            continue;
        }

        double _rho = rho.value[cell.index()];
        const auto metrics = HydroNativePressureDetail::metrics<Axisymmetric>(mesh, cell);

        double _pCoeff[4] = {0.0, 0.0, 0.0, 0.0};
        // west
        if (cell.neiCell[0][1] != -1)
        {
            double _rhoFace = vertRho.value[mesh.verts(mesh.vEdges(cell.cellToVEdge[0][0]).edgeToVert[1]).index()];
            _pCoeff[0] = metrics.west() / _rhoFace;

            Cell& neiCell00 = mesh.cells(cell.neiCell[0][0]);
            Cell& neiCell01 = mesh.cells(cell.neiCell[0][1]);
            if (neiCell00.neiCell[0][1] != -1 && neiCell01.neiCell[0][1] != -1)
            {
                vector<Vector2D> _coords(3);
                _coords[0] = neiCell00.centerCoord();
                Cell& supCell2 = mesh.cells(neiCell00.neiCell[0][1]);
                Cell& supCell3 = mesh.cells(neiCell01.neiCell[0][0]);
                _coords[1] = supCell2.centerCoord();
                _coords[2] = supCell3.centerCoord();
                Vector2D _aimCoord = mesh.verts(neiCell00.cellToVert[3]).coord();
                vector<double> _coeff = Interpolation::leastSquareInterpolationCoeff(_coords, _aimCoord);

                _pressureCoeff[cellI].push_back(Trip(cellI, neiCell00.index(), _coeff[0] * _pCoeff[0]));
                _pressureCoeff[cellI].push_back(Trip(cellI, supCell2.index(), _coeff[1] * _pCoeff[0]));
                _pressureCoeff[cellI].push_back(Trip(cellI, supCell3.index(), _coeff[2] * _pCoeff[0]));
            }
            else if (neiCell01.neiCell[0][1] != -1)
            {
                vector<Vector2D> _coords(3);
                _coords[0] = neiCell00.centerCoord();
                _coords[1] = neiCell01.centerCoord();
                Cell& supCell3 = mesh.cells(neiCell01.neiCell[0][0]);
                _coords[2] = supCell3.centerCoord();
                Vector2D _aimCoord = mesh.verts(neiCell00.cellToVert[3]).coord();
                vector<double> _coeff = Interpolation::leastSquareInterpolationCoeff(_coords, _aimCoord);

                _pressureCoeff[cellI].push_back(Trip(cellI, neiCell00.index(), _coeff[0] * _pCoeff[0]));
                _pressureCoeff[cellI].push_back(Trip(cellI, neiCell01.index(), _coeff[1] * _pCoeff[0]));
                _pressureCoeff[cellI].push_back(Trip(cellI, supCell3.index(), _coeff[2] * _pCoeff[0]));
            }
            else if (neiCell00.neiCell[0][1] != -1)
            {
                vector<Vector2D> _coords(3);
                _coords[0] = neiCell00.centerCoord();
                _coords[1] = neiCell01.centerCoord();
                Cell& supCell3 = mesh.cells(neiCell00.neiCell[0][1]);
                _coords[2] = supCell3.centerCoord();
                Vector2D _aimCoord = mesh.verts(neiCell00.cellToVert[3]).coord();
                vector<double> _coeff = Interpolation::leastSquareInterpolationCoeff(_coords, _aimCoord);

                _pressureCoeff[cellI].push_back(Trip(cellI, neiCell00.index(), _coeff[0] * _pCoeff[0]));
                _pressureCoeff[cellI].push_back(Trip(cellI, neiCell01.index(), _coeff[1] * _pCoeff[0]));
                _pressureCoeff[cellI].push_back(Trip(cellI, supCell3.index(), _coeff[2] * _pCoeff[0]));
            }
            else
            {
                Cell& neiParentCell = mesh.cells(neiCell00.parent);

                for (int i = 0; i < 4; i++)
                {
                    Cell& sonCell = mesh.cells(neiParentCell.son[i]);
                    _pressureCoeff[cellI].push_back(Trip(cellI, sonCell.index(), 0.25 * _pCoeff[0]));
                }
            }
        }
        else if (cell.neiCell[0][0] != -1)
        {
            Cell& neiCell = mesh.cells(cell.neiCell[0][0]);
            
            if (cell.level() > neiCell.level())
            {
                double _rhoNei = 0.0;
                int supSonCellIndex = neiCell.centerCoord().y < cell.centerCoord().y ? neiCell.supCellIndex * 4 + 2 : neiCell.supCellIndex * 4 + 1;
                for (int i = 0; i < supSonCells_[supSonCellIndex].center.ind.size(); i++)
                {
                    _rhoNei += rho.value[supSonCells_[supSonCellIndex].center.ind[i]] * supSonCells_[supSonCellIndex].center.coe[i];
                }
                double _rhoFace = 0.5 * (_rho + _rhoNei);
                _pCoeff[0] = metrics.west() / _rhoFace;
                for (int i = 0; i < supSonCells_[supSonCellIndex].center.ind.size(); i++)
                {
                    _pressureCoeff[cellI].push_back(Trip(cellI, supSonCells_[supSonCellIndex].center.ind[i], supSonCells_[supSonCellIndex].center.coe[i] * _pCoeff[0]));
                }
            }
            else
            {
                double _rhoFace = 0.5 * (_rho + rho.value[neiCell.index()]);
                _pCoeff[0] = metrics.west() / _rhoFace;

                _pressureCoeff[cellI].push_back(Trip(cellI, neiCell.index(), _pCoeff[0]));
            }
        }
        else
        {
            _pCoeff[0] = 0.0;
        }
        // east
        if (cell.neiCell[2][1] != -1)
        {
            double _rhoFace = vertRho.value[mesh.verts(mesh.vEdges(cell.cellToVEdge[1][0]).edgeToVert[1]).index()];
            _pCoeff[2] = metrics.east() / _rhoFace;

            Cell& neiCell20 = mesh.cells(cell.neiCell[2][0]);
            Cell& neiCell21 = mesh.cells(cell.neiCell[2][1]);
            if (neiCell20.neiCell[2][1] != -1 && neiCell21.neiCell[2][1] != -1)
            {
                vector<Vector2D> _coords(3);
                _coords[0] = neiCell20.centerCoord();
                Cell& supCell2 = mesh.cells(neiCell20.neiCell[2][1]);
                Cell& supCell3 = mesh.cells(neiCell21.neiCell[2][0]);
                _coords[1] = supCell2.centerCoord();
                _coords[2] = supCell3.centerCoord();
                Vector2D _aimCoord = mesh.verts(neiCell20.cellToVert[2]).coord();
                vector<double> _coeff = Interpolation::leastSquareInterpolationCoeff(_coords, _aimCoord);

                _pressureCoeff[cellI].push_back(Trip(cellI, neiCell20.index(), _coeff[0] * _pCoeff[2]));
                _pressureCoeff[cellI].push_back(Trip(cellI, supCell2.index(), _coeff[1] * _pCoeff[2]));
                _pressureCoeff[cellI].push_back(Trip(cellI, supCell3.index(), _coeff[2] * _pCoeff[2]));
            }
            else if (neiCell21.neiCell[2][1] != -1)
            {
                vector<Vector2D> _coords(3);
                _coords[0] = neiCell20.centerCoord();
                _coords[1] = neiCell21.centerCoord();
                Cell& supCell3 = mesh.cells(neiCell21.neiCell[2][0]);
                _coords[2] = supCell3.centerCoord();
                Vector2D _aimCoord = mesh.verts(neiCell20.cellToVert[2]).coord();
                vector<double> _coeff = Interpolation::leastSquareInterpolationCoeff(_coords, _aimCoord);

                _pressureCoeff[cellI].push_back(Trip(cellI, neiCell20.index(), _coeff[0] * _pCoeff[2]));
                _pressureCoeff[cellI].push_back(Trip(cellI, neiCell21.index(), _coeff[1] * _pCoeff[2]));
                _pressureCoeff[cellI].push_back(Trip(cellI, supCell3.index(), _coeff[2] * _pCoeff[2]));
            }
            else if (neiCell20.neiCell[2][1] != -1)
            {
                vector<Vector2D> _coords(3);
                _coords[0] = neiCell20.centerCoord();
                _coords[1] = neiCell21.centerCoord();
                Cell& supCell3 = mesh.cells(neiCell20.neiCell[2][1]);
                _coords[2] = supCell3.centerCoord();
                Vector2D _aimCoord = mesh.verts(neiCell20.cellToVert[2]).coord();
                vector<double> _coeff = Interpolation::leastSquareInterpolationCoeff(_coords, _aimCoord);

                _pressureCoeff[cellI].push_back(Trip(cellI, neiCell20.index(), _coeff[0] * _pCoeff[2]));
                _pressureCoeff[cellI].push_back(Trip(cellI, neiCell21.index(), _coeff[1] * _pCoeff[2]));
                _pressureCoeff[cellI].push_back(Trip(cellI, supCell3.index(), _coeff[2] * _pCoeff[2]));
            }
            else
            {
                Cell& neiParentCell = mesh.cells(neiCell20.parent);

                for (int i = 0; i < 4; i++)
                {
                    Cell& sonCell = mesh.cells(neiParentCell.son[i]);
                    _pressureCoeff[cellI].push_back(Trip(cellI, sonCell.index(), 0.25 * _pCoeff[2]));
                }
            }
        }
        else if (cell.neiCell[2][0] != -1)
        {
            Cell& neiCell = mesh.cells(cell.neiCell[2][0]);

            if (cell.level() > neiCell.level())
            {
                double _rhoNei = 0.0;
                int supSonCellIndex = neiCell.centerCoord().y < cell.centerCoord().y ? neiCell.supCellIndex * 4 + 3 : neiCell.supCellIndex * 4;
                for (int i = 0; i < supSonCells_[supSonCellIndex].center.ind.size(); i++)
                {
                    _rhoNei += rho.value[supSonCells_[supSonCellIndex].center.ind[i]] * supSonCells_[supSonCellIndex].center.coe[i];
                }
                double _rhoFace = 0.5 * (_rho + _rhoNei);
                _pCoeff[2] = metrics.east() / _rhoFace;
                for (int i = 0; i < supSonCells_[supSonCellIndex].center.ind.size(); i++)
                {
                    _pressureCoeff[cellI].push_back(Trip(cellI, supSonCells_[supSonCellIndex].center.ind[i], supSonCells_[supSonCellIndex].center.coe[i] * _pCoeff[2]));
                }
            }
            else
            {
                double _rhoFace = 0.5 * (_rho + rho.value[neiCell.index()]);
                _pCoeff[2] = metrics.east() / _rhoFace;

                _pressureCoeff[cellI].push_back(Trip(cellI, neiCell.index(), _pCoeff[2]));
            }
        }
        else
        {
            _pCoeff[2] = 0.0;
        }
        // south
        if (cell.neiCell[1][1] != -1)
        {
            double _rhoFace = vertRho.value[mesh.verts(mesh.hEdges(cell.cellToHEdge[0][0]).edgeToVert[1]).index()];
            _pCoeff[1] = metrics.axial() / _rhoFace;

            Cell& neiCell10 = mesh.cells(cell.neiCell[1][0]);
            Cell& neiCell11 = mesh.cells(cell.neiCell[1][1]);
            if (neiCell10.neiCell[1][1] != -1 && neiCell11.neiCell[1][1] != -1)
            {
                vector<Vector2D> _coords(3);
                _coords[0] = neiCell10.centerCoord();
                Cell& supCell2 = mesh.cells(neiCell10.neiCell[1][1]);
                Cell& supCell3 = mesh.cells(neiCell11.neiCell[1][0]);
                _coords[1] = supCell2.centerCoord();
                _coords[2] = supCell3.centerCoord();
                Vector2D _aimCoord = mesh.verts(neiCell10.cellToVert[1]).coord();
                vector<double> _coeff = Interpolation::leastSquareInterpolationCoeff(_coords, _aimCoord);

                _pressureCoeff[cellI].push_back(Trip(cellI, neiCell10.index(), _coeff[0] * _pCoeff[1]));
                _pressureCoeff[cellI].push_back(Trip(cellI, supCell2.index(), _coeff[1] * _pCoeff[1]));
                _pressureCoeff[cellI].push_back(Trip(cellI, supCell3.index(), _coeff[2] * _pCoeff[1]));
            }
            else if (neiCell11.neiCell[1][1] != -1)
            {
                vector<Vector2D> _coords(3);
                _coords[0] = neiCell10.centerCoord();
                _coords[1] = neiCell11.centerCoord();
                Cell& supCell3 = mesh.cells(neiCell11.neiCell[1][0]);
                _coords[2] = supCell3.centerCoord();
                Vector2D _aimCoord = mesh.verts(neiCell10.cellToVert[1]).coord();
                vector<double> _coeff = Interpolation::leastSquareInterpolationCoeff(_coords, _aimCoord);

                _pressureCoeff[cellI].push_back(Trip(cellI, neiCell10.index(), _coeff[0] * _pCoeff[1]));
                _pressureCoeff[cellI].push_back(Trip(cellI, neiCell11.index(), _coeff[1] * _pCoeff[1]));
                _pressureCoeff[cellI].push_back(Trip(cellI, supCell3.index(), _coeff[2] * _pCoeff[1]));
            }
            else if (neiCell10.neiCell[1][1] != -1)
            {
                vector<Vector2D> _coords(3);
                _coords[0] = neiCell10.centerCoord();
                _coords[1] = neiCell11.centerCoord();
                Cell& supCell3 = mesh.cells(neiCell10.neiCell[1][1]);
                _coords[2] = supCell3.centerCoord();
                Vector2D _aimCoord = mesh.verts(neiCell10.cellToVert[1]).coord();
                vector<double> _coeff = Interpolation::leastSquareInterpolationCoeff(_coords, _aimCoord);

                _pressureCoeff[cellI].push_back(Trip(cellI, neiCell10.index(), _coeff[0] * _pCoeff[1]));
                _pressureCoeff[cellI].push_back(Trip(cellI, neiCell11.index(), _coeff[1] * _pCoeff[1]));
                _pressureCoeff[cellI].push_back(Trip(cellI, supCell3.index(), _coeff[2] * _pCoeff[1]));
            }
            else
            {
                Cell& neiParentCell = mesh.cells(neiCell10.parent);

                for (int i = 0; i < 4; i++)
                {
                    Cell& sonCell = mesh.cells(neiParentCell.son[i]);
                    _pressureCoeff[cellI].push_back(Trip(cellI, sonCell.index(), 0.25 * _pCoeff[1]));
                }
            }
        }
        else if (cell.neiCell[1][0] != -1)
        {
            Cell& neiCell = mesh.cells(cell.neiCell[1][0]);

            if (cell.level() > neiCell.level())
            {
                double _rhoNei = 0.0;
                int supSonCellIndex = neiCell.centerCoord().x < cell.centerCoord().x ? neiCell.supCellIndex * 4 + 2 : neiCell.supCellIndex * 4 + 3;
                for (int i = 0; i < supSonCells_[supSonCellIndex].center.ind.size(); i++)
                {
                    _rhoNei += rho.value[supSonCells_[supSonCellIndex].center.ind[i]] * supSonCells_[supSonCellIndex].center.coe[i];
                }
                double _rhoFace = 0.5 * (_rho + _rhoNei);
                _pCoeff[1] = metrics.axial() / _rhoFace;
                for (int i = 0; i < supSonCells_[supSonCellIndex].center.ind.size(); i++)
                {
                    _pressureCoeff[cellI].push_back(Trip(cellI, supSonCells_[supSonCellIndex].center.ind[i], supSonCells_[supSonCellIndex].center.coe[i] * _pCoeff[1]));
                }
            }
            else
            {
                double _rhoFace = 0.5 * (_rho + rho.value[neiCell.index()]);
                _pCoeff[1] = metrics.axial() / _rhoFace;

                _pressureCoeff[cellI].push_back(Trip(cellI, neiCell.index(), _pCoeff[1]));
            }
        }
        else
        {
            _pCoeff[1] = 0.0;
        }
        // north
        if (cell.neiCell[3][1] != -1)
        {
            double _rhoFace = vertRho.value[mesh.verts(mesh.hEdges(cell.cellToHEdge[1][0]).edgeToVert[1]).index()];
            _pCoeff[3] = metrics.axial() / _rhoFace;
            Cell& neiCell30 = mesh.cells(cell.neiCell[3][0]);
            Cell& neiCell31 = mesh.cells(cell.neiCell[3][1]);
            if (neiCell30.neiCell[3][1] != -1 && neiCell31.neiCell[3][1] != -1)
            {
                vector<Vector2D> _coords(3);
                _coords[0] = neiCell30.centerCoord();
                Cell& supCell2 = mesh.cells(neiCell30.neiCell[3][1]);
                Cell& supCell3 = mesh.cells(neiCell31.neiCell[3][0]);
                _coords[1] = supCell2.centerCoord();
                _coords[2] = supCell3.centerCoord();
                Vector2D _aimCoord = mesh.verts(neiCell30.cellToVert[2]).coord();
                vector<double> _coeff = Interpolation::leastSquareInterpolationCoeff(_coords, _aimCoord);

                _pressureCoeff[cellI].push_back(Trip(cellI, neiCell30.index(), _coeff[0] * _pCoeff[3]));
                _pressureCoeff[cellI].push_back(Trip(cellI, supCell2.index(), _coeff[1] * _pCoeff[3]));
                _pressureCoeff[cellI].push_back(Trip(cellI, supCell3.index(), _coeff[2] * _pCoeff[3]));
            }
            else if (neiCell31.neiCell[3][1] != -1)
            {
                vector<Vector2D> _coords(3);
                _coords[0] = neiCell30.centerCoord();
                _coords[1] = neiCell31.centerCoord();
                Cell& supCell3 = mesh.cells(neiCell31.neiCell[3][0]);
                _coords[2] = supCell3.centerCoord();
                Vector2D _aimCoord = mesh.verts(neiCell30.cellToVert[2]).coord();
                vector<double> _coeff = Interpolation::leastSquareInterpolationCoeff(_coords, _aimCoord);

                _pressureCoeff[cellI].push_back(Trip(cellI, neiCell30.index(), _coeff[0] * _pCoeff[3]));
                _pressureCoeff[cellI].push_back(Trip(cellI, neiCell31.index(), _coeff[1] * _pCoeff[3]));
                _pressureCoeff[cellI].push_back(Trip(cellI, supCell3.index(), _coeff[2] * _pCoeff[3]));
            }
            else if (neiCell30.neiCell[3][1] != -1)
            {
                vector<Vector2D> _coords(3);
                _coords[0] = neiCell30.centerCoord();
                _coords[1] = neiCell31.centerCoord();
                Cell& supCell3 = mesh.cells(neiCell30.neiCell[3][1]);
                _coords[2] = supCell3.centerCoord();
                Vector2D _aimCoord = mesh.verts(neiCell30.cellToVert[2]).coord();
                vector<double> _coeff = Interpolation::leastSquareInterpolationCoeff(_coords, _aimCoord);

                _pressureCoeff[cellI].push_back(Trip(cellI, neiCell30.index(), _coeff[0] * _pCoeff[3]));
                _pressureCoeff[cellI].push_back(Trip(cellI, neiCell31.index(), _coeff[1] * _pCoeff[3]));
                _pressureCoeff[cellI].push_back(Trip(cellI, supCell3.index(), _coeff[2] * _pCoeff[3]));
            }
            else
            {
                Cell& neiParentCell = mesh.cells(neiCell30.parent);

                for (int i = 0; i < 4; i++)
                {
                    Cell& sonCell = mesh.cells(neiParentCell.son[i]);
                    _pressureCoeff[cellI].push_back(Trip(cellI, sonCell.index(), 0.25 * _pCoeff[3]));
                }
            }
        }
        else if (cell.neiCell[3][0] != -1)
        {
            Cell& neiCell = mesh.cells(cell.neiCell[3][0]);

            if (cell.level() > neiCell.level())
            {
                double _rhoNei = 0.0;
                int supSonCellIndex = neiCell.centerCoord().x < cell.centerCoord().x ? neiCell.supCellIndex * 4 + 1 : neiCell.supCellIndex * 4;
                for (int i = 0; i < supSonCells_[supSonCellIndex].center.ind.size(); i++)
                {
                    _rhoNei += rho.value[supSonCells_[supSonCellIndex].center.ind[i]] * supSonCells_[supSonCellIndex].center.coe[i];
                }
                double _rhoFace = 0.5 * (_rho + _rhoNei);
                _pCoeff[3] = metrics.axial() / _rhoFace;
                for (int i = 0; i < supSonCells_[supSonCellIndex].center.ind.size(); i++)
                {
                    _pressureCoeff[cellI].push_back(Trip(cellI, supSonCells_[supSonCellIndex].center.ind[i], supSonCells_[supSonCellIndex].center.coe[i] * _pCoeff[3]));
                }
            }
            else
            {
                double _rhoFace = 0.5 * (_rho + rho.value[neiCell.index()]);
                _pCoeff[3] = metrics.axial() / _rhoFace;

                _pressureCoeff[cellI].push_back(Trip(cellI, neiCell.index(), _pCoeff[3]));
            }
        }
        else
        {
            _pCoeff[3] = 0.0;
        }

        _pressureCoeff[cellI].push_back(Trip(cellI, cellI, -(_pCoeff[0] + _pCoeff[1] + _pCoeff[2] + _pCoeff[3])));

        // RHS
        double _uWest = 0.0;
        if (cell.cellToVEdge[0][1] != -1)
        {
            for (int i = 0; i < supCells_[cell.supCellIndex].vEdges[0].ind.size(); i++)
            {
                _uWest += uxMid.value[supCells_[cell.supCellIndex].vEdges[0].ind[i]] * supCells_[cell.supCellIndex].vEdges[0].coe[i];
            }
        }
        else
        {
            _uWest = uxMid.value[mesh.vEdges(cell.cellToVEdge[0][0]).index()];
        }
        double _uEast = 0.0;
        if (cell.cellToVEdge[1][1] != -1)
        {
            for (int i = 0; i < supCells_[cell.supCellIndex].vEdges[1].ind.size(); i++)
            {
                _uEast += uxMid.value[supCells_[cell.supCellIndex].vEdges[1].ind[i]] * supCells_[cell.supCellIndex].vEdges[1].coe[i];
            }
        }
        else
        {
            _uEast = uxMid.value[mesh.vEdges(cell.cellToVEdge[1][0]).index()];
        }
        double _vSouth = 0.0;
        if (cell.cellToHEdge[0][1] != -1)
        {
            for (int i = 0; i < supCells_[cell.supCellIndex].hEdges[0].ind.size(); i++)
            {
                _vSouth += uyMid.value[supCells_[cell.supCellIndex].hEdges[0].ind[i]] * supCells_[cell.supCellIndex].hEdges[0].coe[i];
            }
        }
        else
        {
            _vSouth = uyMid.value[mesh.hEdges(cell.cellToHEdge[0][0]).index()];
        }
        double _vNorth = 0.0;
        if (cell.cellToHEdge[1][1] != -1)
        {
            for (int i = 0; i < supCells_[cell.supCellIndex].hEdges[1].ind.size(); i++)
            {
                _vNorth += uyMid.value[supCells_[cell.supCellIndex].hEdges[1].ind[i]] * supCells_[cell.supCellIndex].hEdges[1].coe[i];
            }
        }
        else
        {
            _vNorth = uyMid.value[mesh.hEdges(cell.cellToHEdge[1][0]).index()];
        }

        _pressureRHS[cellI] = metrics.fluxRate(dt, _uEast, _uWest, _vNorth, _vSouth);
    }

    vector<Eigen::Triplet<double>> _pressureTrip;
    for (int i = 0; i < mesh.cellIndexNum(); i++)
    {
        _pressureTrip.insert(_pressureTrip.end(), _pressureCoeff[i].begin(), _pressureCoeff[i].end());
    }

    SparseMatrix<double> _pressureMatrix(mesh.cellIndexNum(), mesh.cellIndexNum());
    _pressureMatrix.setFromTriplets(_pressureTrip.begin(), _pressureTrip.end());

    return {std::move(_pressureMatrix), std::move(_pressureRHS)};
}


// Closed native adaptive projection: assemble EACH real face once and insert
// its identical corrected flux, with opposite signs, in both neighboring FV
// rows. Reuse the old virtual-child gradient and the existing final-velocity
// update. Same-level and open/fitted paths retain their original assembly.
template<typename Mesh>
HydroPressureSystem Hydrodynamics::assembleNativeAdaptivePressure(Mesh& mesh, double dt, bool axisymmetric)
{
    const int count = mesh.cellIndexNum();
    Eigen::VectorXd rhs = Eigen::VectorXd::Zero(count);
    std::vector<Eigen::Triplet<double>> entries;
    const auto centers = [&](const Cell& cell, int child) -> const supCoeff& {
        return supSonCells_[cell.supCellIndex * 4 + child].center;
    };
    const auto face = [&](const Edge& edge, bool horizontal, double predictor) {
        if (edge.masterCell < 0 || edge.slaveCell < 0) return; // prescribed zero normal flux
        const auto data = pressureFaceData(mesh, edge, horizontal, p);
        const double area = edge.length() * (axisymmetric ? 2.0 * Pi * edge.midpointCoord().x : 1.0);
        if (!(area > 0.0) || !(data.density > 0.0) || !std::isfinite(data.density))
            throw std::runtime_error("native adaptive pressure has invalid face area/density");
        const int high = mesh.cells(edge.masterCell).index(), low = mesh.cells(edge.slaveCell).index();
        const double fluxRate = area * predictor / dt;
        if (low != 0) rhs[low] += fluxRate;
        if (high != 0) rhs[high] -= fluxRate;
        HydroNativePressureStencil::visitGradient(mesh, edge, horizontal, centers,
            [&](int column, double gradient) {
                const double coefficient = area / data.density * gradient;
                if (low != 0) entries.emplace_back(low, column, coefficient);
                if (high != 0) entries.emplace_back(high, column, -coefficient);
            });
    };
    for (int i = 0; i < mesh.vEdgeIndexNum(); ++i)
        face(mesh.vEdges(mesh.vEdgeIndexToLabel[i]), false, uxMid.value[i]);
    for (int i = 0; i < mesh.hEdgeIndexNum(); ++i)
        face(mesh.hEdges(mesh.hEdgeIndexToLabel[i]), true, uyMid.value[i]);
    entries.emplace_back(0, 0, 1.0); // same closed-domain pressure gauge
    Eigen::SparseMatrix<double> matrix(count, count);
    matrix.setFromTriplets(entries.begin(), entries.end());
    matrix.makeCompressed();
    return {std::move(matrix), std::move(rhs)};
}

// One face stencil for the historical pressure in the momentum predictor and
// the pressure potential in the velocity projection. In particular, never
// compute a cell gradient and interpolate it back to the staggered face.
template<typename Mesh>
Hydrodynamics::PressureFaceData Hydrodynamics::pressureFaceData(Mesh& mesh,
    const Edge& edge, bool horizontal, const CellField<double>& pressure) const
{
    return HydroNativePressureStencil::faceData(mesh, edge, horizontal, pressure, rho,
        [&](const Cell& cell, int child) -> const supCoeff& {
            return supSonCells_[cell.supCellIndex * 4 + child].center;
        });
}

template<typename Mesh>
void Hydrodynamics::addFaceForce(Mesh& mesh, bool axisymmetric,
    const VEdgeField<double>& forceX, const HEdgeField<double>& forceY, Eigen::VectorXd& rhs)
{
    auto add = [&](const Edge& edge, bool horizontal, int row, double force)
    {
        if (edge.masterCell == -1 || edge.slaveCell == -1)
        {
            if (!boundaryFaces_.empty() && boundaryFaces_[row].outlet)
            {
                const double load = boundaryFaces_[row].volume * force;
                if (!std::isfinite(load) || !std::isfinite(rhs[row] + load))
                    throw std::runtime_error("invalid outlet half-control-volume face load");
                rhs[row] += load;
            }
            return;
        }
        const auto face = pressureFaceData(mesh, edge, horizontal, p);
        const double measure = face.dualArea
            * (axisymmetric ? 2.0 * Pi * edge.midpointCoord().x : 1.0);
        const double load = measure * force;
        if (!std::isfinite(measure) || measure <= 0.0 || !std::isfinite(load)
            || !std::isfinite(face.density) || face.density <= 0.0
            || !std::isfinite(rhs[row] + load))
            throw std::runtime_error("explicit face force load, dual volume or density is invalid");
        rhs[row] += load;
    };
    // Match the pressure-load quadrature exactly; density is consumed once by
    // the momentum mass matrix, not divided out again in this force density.
    for (int i = 0; i < mesh.vEdgeIndexNum(); ++i)
        add(mesh.vEdges(mesh.vEdgeIndexToLabel[i]), false, i, forceX.value[i]);
    for (int i = 0; i < mesh.hEdgeIndexNum(); ++i)
        add(mesh.hEdges(mesh.hEdgeIndexToLabel[i]), true, i + mesh.vEdgeIndexNum(), forceY.value[i]);
}

template<typename Mesh>
void Hydrodynamics::addHistoricalPressure(Mesh& mesh, bool axisymmetric, Eigen::VectorXd& rhs)
{
    auto subtract = [&](const Edge& edge, bool horizontal, int row)
    {
        if (edge.masterCell == -1 || edge.slaveCell == -1)
        {
            if (!boundaryFaces_.empty() && boundaryFaces_[row].outlet)
            {
                const auto& boundary = boundaryFaces_[row];
                rhs[row] += boundary.sign * boundary.area * p.value[boundary.owner];
                if (!std::isfinite(rhs[row]))
                    throw std::runtime_error("invalid outlet historical pressure load");
            }
            return;
        }
        const auto face = pressureFaceData(mesh, edge, horizontal, p);
        const double measure = face.dualArea
            * (axisymmetric ? 2.0 * Pi * edge.midpointCoord().x : 1.0);
        const double load = measure * face.gradient;
        if (!std::isfinite(load) || !std::isfinite(face.density) || face.density <= 0.0)
            throw std::runtime_error("historical pressure face load or density is invalid");
        rhs[row] -= load;
    };
    // This checked pass stays serial: exceptions must not leave an OpenMP loop.
    for (int i = 0; i < mesh.vEdgeIndexNum(); ++i)
        subtract(mesh.vEdges(mesh.vEdgeIndexToLabel[i]), false, i);
    for (int i = 0; i < mesh.hEdgeIndexNum(); ++i)
        subtract(mesh.hEdges(mesh.hEdgeIndexToLabel[i]), true, i + mesh.vEdgeIndexNum());
}
