#pragma once

using namespace std;

template<typename T>
VertField<T>::VertField() : Field<T>() {}

template<typename T>
VertField<T>::VertField(const VertField<T>& field) : Field<T>(field) {}

template<typename T>
template<typename Mesh>
VertField<T>::VertField(Mesh& mesh, const T& value) : Field<T>()
{
    this->size_ = mesh.vertIndexNum();
    this->value = std::make_unique<T[]>(this->size_);
    this->label_ = std::make_unique<int[]>(this->size_);
#pragma omp parallel for num_threads(thread_num)
    for (int i = 0; i < this->size_; i++)
    {
        Vert& vert = mesh.verts(mesh.vertIndexToLabel[i]);

        this->value[i] = value;
        this->label_[i] = vert.label();
    }
}

template<typename T>
VertField<T>::~VertField() {}

template<typename T>
template<typename Mesh>
void VertField<T>::initialize(const Mesh& mesh, const T& value)
{
    this->size_ = mesh.vertIndexNum();

    this->value = std::make_unique<T[]>(this->size_);
    this->label_ = std::make_unique<int[]>(this->size_);
#pragma omp parallel for num_threads(thread_num)
    for (int i = 0; i < this->size_; i++)
    {
        Vert& vert = mesh.verts(mesh.vertIndexToLabel[i]);

        this->value[i] = value;
        this->label_[i] = vert.label();
    }
}

template<typename T>
template<typename Mesh>
void VertField<T>::updateWithMesh(const Mesh& mesh)
{
    auto tempValue = std::make_unique<T[]>(mesh.vertIndexNum());
    auto tempLabel = std::make_unique<int[]>(mesh.vertIndexNum());

    int _ptr1 = 0;
    int _ptr2 = 0;
    while (_ptr1 < this->size_ && _ptr2 < mesh.vertIndexNum())
    {
        if (this->label_[_ptr1] == mesh.verts(mesh.vertIndexToLabel[_ptr2]).label())
        {
            mesh.verts(mesh.vertIndexToLabel[_ptr2]).setNeedUpdate(false);
            tempValue[_ptr2] = this->value[_ptr1];
            tempLabel[_ptr2] = this->label_[_ptr1];
            _ptr1++;
            _ptr2++;
        }
        else if (this->label_[_ptr1] < mesh.verts(mesh.vertIndexToLabel[_ptr2]).label())
        {
            _ptr1++;
        }
        else
        {
            mesh.verts(mesh.vertIndexToLabel[_ptr2]).setNeedUpdate(true);
            tempLabel[_ptr2] = mesh.verts(mesh.vertIndexToLabel[_ptr2]).label();
            _ptr2++;
        }
    }

    while (_ptr2 < mesh.vertIndexNum())
    {
        mesh.verts(mesh.vertIndexToLabel[_ptr2]).setNeedUpdate(true);
        tempLabel[_ptr2] = mesh.verts(mesh.vertIndexToLabel[_ptr2]).label();
        _ptr2++;
    }

    this->value = std::move(tempValue);
    this->label_ = std::move(tempLabel);
    this->size_ = mesh.vertIndexNum();

#pragma omp parallel for num_threads(thread_num)
    for (int cellI = 0; cellI < mesh.Nx * mesh.Ny; cellI++)
    {
        Cell& cell = mesh.cells(cellI);

        if (cell.son[0] != -1)
        {
            for (int i = 0; i < 4; i++)
            {
                subInterpolation(mesh, mesh.cells(cell.son[i]));
            }
        }
    }
}

template<typename T>
template<typename Mesh>
void VertField<T>::subInterpolation(const Mesh& mesh, const Cell& cell)
{
    Cell& parentCell = mesh.cells(cell.parent);

    for (int i = 0; i < 4; i++)
    {
        Vert& vert = mesh.verts(cell.cellToVert[i]);
        if (vert.needUpdate())
        {
            vert.setNeedUpdate(false);

            this->value[vert.index()] = Interpolation::biLinearInterpolation
                (
                    mesh, 
                    *this, 
                    parentCell, 
                    mesh.verts(cell.cellToVert[i]).coord()
                );
        }
    }

    if (cell.son[0] != -1)
    {
        for (int i = 0; i < 4; i++)
        {
            subInterpolation(mesh, mesh.cells(cell.son[i]));
        }
    }
}

template<typename T>
template<typename Mesh>
unique_ptr<CellField<T>> VertField<T>::toCellField(const Mesh& mesh) const
{
    auto cellField = make_unique<CellField<T>>();
    cellField->initialize(mesh, 0.0);

#pragma omp parallel for num_threads(thread_num)
    for (int i = 0; i < mesh.cellIndexNum(); i++)
    {
        Cell& cell = mesh.cells(mesh.cellIndexToLabel[i]);

        // cellField->value[i] = Interpolation::biLinearInterpolation
        //     (
        //         mesh, 
        //         *this, 
        //         cell, 
        //         cell.centerCoord()
        //     );

        vector<double> values(4);
        vector<Vector2D> coords(4);
        for (int i = 0; i < 4; i++)
        {
            Vert& vert = mesh.verts(cell.cellToVert[i]);
            values[i] = this->value[vert.index()];
            coords[i] = vert.coord();
        }

        cellField->value[i] = Interpolation::leastSquareInterpolation(coords, values, cell.centerCoord(), false);
    }

    return cellField;
}

template<typename T>
template<typename Mesh>
std::unique_ptr<HEdgeField<T>> VertField<T>::toHEdgeField(const Mesh& mesh) const
{
    auto hEdgeField = make_unique<HEdgeField<T>>();
    hEdgeField->initialize(mesh, 0.0);

#pragma omp parallel for num_threads(thread_num)
    for (int i = 0; i < mesh.hEdgeIndexNum(); i++)
    {
        Edge& hEdge = mesh.hEdges(mesh.hEdgeIndexToLabel[i]);
        Cell& cell = mesh.cells(hEdge.masterCell == -1 ? hEdge.slaveCell : hEdge.masterCell);

        // hEdgeField->value[i] = Interpolation::biLinearInterpolation
        //     (
        //         mesh, 
        //         *this, 
        //         cell,
        //         hEdge.midpointCoord()
        //     );

        vector<double> values(4);
        vector<Vector2D> coords(4);
        for (int i = 0; i < 4; i++)
        {
            Vert& vert = mesh.verts(cell.cellToVert[i]);
            values[i] = this->value[vert.index()];
            coords[i] = vert.coord();
        }

        hEdgeField->value[i] = Interpolation::leastSquareInterpolation(coords, values, hEdge.midpointCoord(), false);
    }

    return hEdgeField;
}

template<typename T>
template<typename Mesh>
std::unique_ptr<VEdgeField<T>> VertField<T>::toVEdgeField(const Mesh& mesh) const
{
    auto vEdgeField = make_unique<VEdgeField<T>>();
    vEdgeField->initialize(mesh, 0.0);

#pragma omp parallel for num_threads(thread_num)
    for (int i = 0; i < mesh.vEdgeIndexNum(); i++)
    {
        Edge& vEdge = mesh.vEdges(mesh.vEdgeIndexToLabel[i]);
        Cell& cell = mesh.cells(vEdge.masterCell == -1 ? vEdge.slaveCell : vEdge.masterCell);

        // vEdgeField->value[i] = Interpolation::biLinearInterpolation
        //     (
        //         mesh, 
        //         *this, 
        //         cell,
        //         vEdge.midpointCoord()
        //     );

        vector<double> values(4);
        vector<Vector2D> coords(4);
        for (int i = 0; i < 4; i++)
        {
            Vert& vert = mesh.verts(cell.cellToVert[i]);
            values[i] = this->value[vert.index()];
            coords[i] = vert.coord();
        }

        vEdgeField->value[i] = Interpolation::leastSquareInterpolation(coords, values, vEdge.midpointCoord(), false);
    }

    return vEdgeField;
}
