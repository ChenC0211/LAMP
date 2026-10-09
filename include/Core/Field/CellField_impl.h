#pragma once

#include "Core/Base/Interpolation.h"
#include "Core/Field/FieldTransferValidation.h"
using namespace std;

template<typename T>
CellField<T>::CellField() : Field<T>() {}

template<typename T>
CellField<T>::CellField(const CellField<T>& field) : Field<T>(field) {}

template<typename T>
template<typename Mesh>
CellField<T>::CellField(Mesh& mesh, const T& value) : Field<T>()
{
    this->size_ = mesh.cellIndexNum();
    this->value = std::make_unique<T[]>(this->size_);
    this->label_ = std::make_unique<int[]>(this->size_);
#pragma omp parallel for num_threads(thread_num)
    for (int i = 0; i < this->size_; i++)
    {
        Cell& cell = mesh.cells(mesh.cellIndexToLabel[i]);

        this->value[i] = value;
        this->label_[i] = cell.label();
    }
}

template<typename T>
CellField<T>::~CellField() {}

template<typename T>
template<typename Mesh>
void CellField<T>::initialize(const Mesh& mesh, const T& value)
{
    this->size_ = mesh.cellIndexNum();
    this->value = std::make_unique<T[]>(this->size_);
    this->label_ = std::make_unique<int[]>(this->size_);
#pragma omp parallel for num_threads(thread_num)
    for (int i = 0; i < this->size_; i++)
    {
        Cell& cell = mesh.cells(mesh.cellIndexToLabel[i]);

        this->value[i] = value;
        this->label_[i] = cell.label();
    }
}

template<typename T>
template<typename Mesh>
void CellField<T>::updateWithMesh(const Mesh& mesh)
{
    auto tempValue = std::make_unique<T[]>(mesh.cellIndexNum());
    auto tempLabel = std::make_unique<int[]>(mesh.cellIndexNum());

    int _ptr1 = 0;
    int _ptr2 = 0;
    while (_ptr1 < this->size_ && _ptr2 < mesh.cellIndexNum())
    {
        if (this->label_[_ptr1] == mesh.cells(mesh.cellIndexToLabel[_ptr2]).label())
        {
            // This field already owns the surviving value. A cached/shared
            // work flag is not evidence that THIS field needs interpolation.
            mesh.cells(mesh.cellIndexToLabel[_ptr2]).setNeedUpdate(false);
            tempValue[_ptr2] = this->value[_ptr1];
            tempLabel[_ptr2] = this->label_[_ptr1];
            _ptr1++;
            _ptr2++;
        }
        else if (this->label_[_ptr1] < mesh.cells(mesh.cellIndexToLabel[_ptr2]).label())
        {
            _ptr1++;
        }
        else
        {
            mesh.cells(mesh.cellIndexToLabel[_ptr2]).setNeedUpdate(true);
            tempLabel[_ptr2] = mesh.cells(mesh.cellIndexToLabel[_ptr2]).label();
            _ptr2++;
        }
    }

    while (_ptr2 < mesh.cellIndexNum())
    {
        mesh.cells(mesh.cellIndexToLabel[_ptr2]).setNeedUpdate(true);
        tempLabel[_ptr2] = mesh.cells(mesh.cellIndexToLabel[_ptr2]).label();
        _ptr2++;
    }

    for (int i = 0; i < mesh.cellIndexNum(); ++i)
    {
        const Cell& cell = mesh.cells(mesh.cellIndexToLabel[i]);
        if (cell.needUpdate()) Quantum::FieldTransfer::validateSources(cell, this->size_, "cell");
    }
#pragma omp parallel for num_threads(thread_num)
    for (int i = 0; i < mesh.cellIndexNum(); i++)
    {
        Cell& cell = mesh.cells(mesh.cellIndexToLabel[i]);
        if (cell.needUpdate())
        {
            cell.setNeedUpdate(false);

            vector<double> values(cell.relatedPreIndex.size());
            for (int j = 0; j < cell.relatedPreIndex.size(); j++)
            {
                values[j] = this->value[cell.relatedPreIndex[j]];
            }

            // Robustness: leastSquareInterpolation (first-order) requires >=3 points.
            // Coarse-fine updates can yield small stencils.
            if (cell.relatedPreCoord.size() < 3 || cell.relatedPreCoord.size() != values.size())
            {
                if (values.empty())
                {
                    tempValue[i] = T{};
                }
                else
                {
                    double sum = 0.0;
                    for (const double v : values)
                        sum += v;
                    tempValue[i] = static_cast<T>(sum / static_cast<double>(values.size()));
                }
            }
            else
            {
                tempValue[i] = Interpolation::leastSquareInterpolation
                (
                    cell.relatedPreCoord,
                    values,
                    cell.centerCoord(),
                    false
                );
            }
        }
    }

    this->value = std::move(tempValue);
    this->label_ = std::move(tempLabel);
    this->size_ = mesh.cellIndexNum();
}

// template<typename T>
// template<typename Mesh>
// void CellField<T>::subInterpolation(const Mesh& mesh, std::unique_ptr<T[]>& values, unique_ptr<int[]>& labels, int preI, Cell& cell)
// {
//     if (cell.son[0] != -1)
//     {
//         for (int i = 0; i < 4; i++)
//         {
//             subInterpolation(mesh, values, labels, preI, mesh.cells(cell.son[i]));
//         }
//     }
//     else if (cell.index() != -1)
//     {
//         values[cell.index()] = Interpolation::averageInterpolation
//             (
//                 mesh, 
//                 *this, 
//                 preI
//             );
//         labels[cell.index()] = cell.label();
//     }
// }

template<typename T>
template<typename Mesh>
unique_ptr<VertField<T>> CellField<T>::toVertField(const Mesh& mesh) const
{
    auto vertField = make_unique<VertField<T>>();
    vertField->initialize(mesh, 0.0);

#pragma omp parallel for num_threads(thread_num)
    for (int i = 0; i < mesh.vertIndexNum(); i++)
    {
        Vert& vert = mesh.verts(mesh.vertIndexToLabel[i]);

        vector<double> values;
        vector<Vector2D> coords;
        if (vert.boundary() == Boundary::non)
        {
            for (int j = 0; j < 4; j++)
            {
                if (vert.vertToCell[j] != -1)
                {
                    values.push_back(this->value[mesh.cells(vert.vertToCell[j]).index()]);
                    coords.push_back(mesh.cells(vert.vertToCell[j]).centerCoord());
                }
            }
        }
        else
        {
            if (vert.boundary() == Boundary::west)
            {
                Vert& stencilVert = mesh.verts(vert.neiVertices[2]);
                for (int j = 0; j < 4; j++)
                {
                    if (stencilVert.vertToCell[j] != -1)
                    {
                        values.push_back(this->value[mesh.cells(stencilVert.vertToCell[j]).index()]);
                        coords.push_back(mesh.cells(stencilVert.vertToCell[j]).centerCoord());
                    }
                }
            }
            else if (vert.boundary() == Boundary::south)
            {
                Vert& stencilVert = mesh.verts(vert.neiVertices[3]);
                for (int j = 0; j < 4; j++)
                {
                    if (stencilVert.vertToCell[j] != -1)
                    {
                        values.push_back(this->value[mesh.cells(stencilVert.vertToCell[j]).index()]);
                        coords.push_back(mesh.cells(stencilVert.vertToCell[j]).centerCoord());
                    }
                }
            }
            else if (vert.boundary() == Boundary::east)
            {
                Vert& stencilVert = mesh.verts(vert.neiVertices[0]);
                for (int j = 0; j < 4; j++)
                {
                    if (stencilVert.vertToCell[j] != -1)
                    {
                        values.push_back(this->value[mesh.cells(stencilVert.vertToCell[j]).index()]);
                        coords.push_back(mesh.cells(stencilVert.vertToCell[j]).centerCoord());
                    }
                }
            }
            else if (vert.boundary() == Boundary::north)
            {
                Vert& stencilVert = mesh.verts(vert.neiVertices[1]);
                for (int j = 0; j < 4; j++)
                {
                    if (stencilVert.vertToCell[j] != -1)
                    {
                        values.push_back(this->value[mesh.cells(stencilVert.vertToCell[j]).index()]);
                        coords.push_back(mesh.cells(stencilVert.vertToCell[j]).centerCoord());
                    }
                }
            }
            else if (vert.boundary() == Boundary::southWest)
            {
                Vert& stencilVert = mesh.verts(mesh.cells(vert.vertToCell[2]).cellToVert[2]);
                for (int j = 0; j < 4; j++)
                {
                    if (stencilVert.vertToCell[j] != -1)
                    {
                        values.push_back(this->value[mesh.cells(stencilVert.vertToCell[j]).index()]);
                        coords.push_back(mesh.cells(stencilVert.vertToCell[j]).centerCoord());
                    }
                }
            }
            else if (vert.boundary() == Boundary::southEast)
            {
                Vert& stencilVert = mesh.verts(mesh.cells(vert.vertToCell[3]).cellToVert[3]);
                for (int j = 0; j < 4; j++)
                {
                    if (stencilVert.vertToCell[j] != -1)
                    {
                        values.push_back(this->value[mesh.cells(stencilVert.vertToCell[j]).index()]);
                        coords.push_back(mesh.cells(stencilVert.vertToCell[j]).centerCoord());
                    }
                }
            }
            else if (vert.boundary() == Boundary::northEast)
            {
                Vert& stencilVert = mesh.verts(mesh.cells(vert.vertToCell[0]).cellToVert[0]);
                for (int j = 0; j < 4; j++)
                {
                    if (stencilVert.vertToCell[j] != -1)
                    {
                        values.push_back(this->value[mesh.cells(stencilVert.vertToCell[j]).index()]);
                        coords.push_back(mesh.cells(stencilVert.vertToCell[j]).centerCoord());
                    }
                }
            }
            else if (vert.boundary() == Boundary::northWest)
            {
                Vert& stencilVert = mesh.verts(mesh.cells(vert.vertToCell[1]).cellToVert[1]);
                for (int j = 0; j < 4; j++)
                {
                    if (stencilVert.vertToCell[j] != -1)
                    {
                        values.push_back(this->value[mesh.cells(stencilVert.vertToCell[j]).index()]);
                        coords.push_back(mesh.cells(stencilVert.vertToCell[j]).centerCoord());
                    }
                }
            }
        }
        
        // Robustness: leastSquareInterpolation (first-order) requires >=3 points.
        // On adaptive/coarse-fine meshes some boundary/corner vertices may see fewer
        // adjacent leaf cells. Fall back to a simple average to avoid aborting.
        if (coords.size() < 3 || coords.size() != values.size())
        {
            if (values.empty())
            {
                vertField->value[i] = T{};
            }
            else
            {
                double sum = 0.0;
                for (const double v : values)
                    sum += v;
                vertField->value[i] = static_cast<T>(sum / static_cast<double>(values.size()));
            }
        }
        else
        {
            vertField->value[i] = Interpolation::leastSquareInterpolation
                (
                    coords,
                    values,
                    vert.coord(),
                    false
                );
        }
    }

    return vertField;
}

template<typename T>
template<typename Mesh>
void CellField<T>::output(const Mesh& mesh, const string& _name) const
{
    // CellField插值到VertField
    auto vertField = toVertField(mesh);

    // 输出到文件
    vertField->output(mesh, _name);
}
