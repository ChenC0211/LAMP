#pragma once

#include "Core/Base/Interpolation.h"
#include "Core/Field/FieldTransferValidation.h"
using namespace std;

template<typename T>
VEdgeField<T>::VEdgeField() : Field<T>() {}

template<typename T>
VEdgeField<T>::VEdgeField(const VEdgeField<T>& field) : Field<T>(field) {}

template<typename T>
template<typename Mesh>
VEdgeField<T>::VEdgeField(Mesh& mesh, const T& value) : Field<T>()
{
    this->size_ = mesh.vEdgeIndexNum();
    this->value = std::make_unique<T[]>(this->size_);
    this->label_ = std::make_unique<int[]>(this->size_);
#pragma omp parallel for
    for (int i = 0; i < this->size_; i++)
    {
        Edge& vEdge = mesh.vEdges(mesh.vEdgeIndexToLabel[i]);

        this->value[i] = value;
        this->label_[i] = vEdge.label();
    }
}

template<typename T>
VEdgeField<T>::~VEdgeField() {}

template<typename T>
template<typename Mesh>
void VEdgeField<T>::initialize(const Mesh& mesh, const T& value)
{
    this->size_ = mesh.vEdgeIndexNum();
    this->value = std::make_unique<T[]>(this->size_);
    this->label_ = std::make_unique<int[]>(this->size_);
#pragma omp parallel for
    for (int i = 0; i < this->size_; i++)
    {
        Edge& vEdge = mesh.vEdges(mesh.vEdgeIndexToLabel[i]);

        this->value[i] = value;
        this->label_[i] = vEdge.label();
    }
}

template<typename T>
template<typename Mesh>
void VEdgeField<T>::updateWithMesh(const Mesh& mesh)
{
    auto tempValue = std::make_unique<T[]>(mesh.vEdgeIndexNum());
    auto tempLabel = std::make_unique<int[]>(mesh.vEdgeIndexNum());

    int _ptr1 = 0;
    int _ptr2 = 0;
    while (_ptr1 < this->size_ && _ptr2 < mesh.vEdgeIndexNum())
    {
        if (this->label_[_ptr1] == mesh.vEdges(mesh.vEdgeIndexToLabel[_ptr2]).label())
        {
            mesh.vEdges(mesh.vEdgeIndexToLabel[_ptr2]).setNeedUpdate(false);
            tempValue[_ptr2] = this->value[_ptr1];
            tempLabel[_ptr2] = this->label_[_ptr1];
            _ptr1++;
            _ptr2++;
        }
        else if (this->label_[_ptr1] < mesh.vEdges(mesh.vEdgeIndexToLabel[_ptr2]).label())
        {
            _ptr1++;
        }
        else
        {
            mesh.vEdges(mesh.vEdgeIndexToLabel[_ptr2]).setNeedUpdate(true);
            tempLabel[_ptr2] = mesh.vEdges(mesh.vEdgeIndexToLabel[_ptr2]).label();
            _ptr2++;
        }
    }

    while (_ptr2 < mesh.vEdgeIndexNum())
    {
        mesh.vEdges(mesh.vEdgeIndexToLabel[_ptr2]).setNeedUpdate(true);
        tempLabel[_ptr2] = mesh.vEdges(mesh.vEdgeIndexToLabel[_ptr2]).label();
        _ptr2++;
    }

    for (int i = 0; i < mesh.vEdgeIndexNum(); ++i)
    {
        const Edge& edge = mesh.vEdges(mesh.vEdgeIndexToLabel[i]);
        if (edge.needUpdate()) LAMP::FieldTransfer::validateSources(edge, this->size_, "vertical edge");
    }
#pragma omp parallel for
    for (int i = 0; i < mesh.vEdgeIndexNum(); i++)
    {
        Edge& vEdge = mesh.vEdges(mesh.vEdgeIndexToLabel[i]);
        if (vEdge.needUpdate())
        {
            vEdge.setNeedUpdate(false);

            vector<double> values(vEdge.relatedPreIndex.size());
            for (int j = 0; j < vEdge.relatedPreIndex.size(); j++)
            {
                values[j] = this->value[vEdge.relatedPreIndex[j]];
            }

            // Robustness: leastSquareInterpolation (first-order) requires >=3 points.
            if (vEdge.relatedPreCoord.size() < 3 || vEdge.relatedPreCoord.size() != values.size())
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
                    vEdge.relatedPreCoord,
                    values,
                    vEdge.midpointCoord(),
                    false
                );
            }
        }
    }

    this->value = std::move(tempValue);
    this->label_ = std::move(tempLabel);
    this->size_ = mesh.vEdgeIndexNum();
}

// template<typename T>
// template<typename Mesh>
// void VEdgeField<T>::subInterpolation(const Mesh& mesh, std::unique_ptr<T[]>& values, unique_ptr<int[]>& labels, int preI, Edge& vEdge)
// {
//     if (vEdge.son[0] != -1)
//     {
//         for (int i = 0; i < 2; i++)
//         {
//             subInterpolation(mesh, values, labels, preI, mesh.vEdges(vEdge.son[i]));
//         }
//     }
//     else if (vEdge.index() != -1)
//     {
//         values[vEdge.index()] = Interpolation::averageInterpolation
//             (
//                 mesh, 
//                 *this, 
//                 preI
//             );
//         labels[vEdge.index()] = vEdge.label();
//     }
// }

template<typename T>
template<typename Mesh>
unique_ptr<VertField<T>> VEdgeField<T>::toVertField(const Mesh& mesh) const
{
    auto vertField = make_unique<VertField<T>>();
    vertField->initialize(mesh, 0.0);

#pragma omp parallel for
    for (int i = 0; i < mesh.vertIndexNum(); i++)
    {
        Vert& vert = mesh.verts(mesh.vertIndexToLabel[i]);
        
        vector<double> values;
        vector<Vector2D> coords;

        Vert* stencilVert;
        if (vert.boundary() == Boundary::non)
        {
            stencilVert = &vert;
        }
        else if (vert.boundary() == Boundary::west)
        {
            stencilVert = &mesh.verts(vert.neiVertices[2]);
        }
        else if (vert.boundary() == Boundary::south)
        {
            stencilVert = &mesh.verts(vert.neiVertices[3]);
        }
        else if (vert.boundary() == Boundary::east)
        {
            stencilVert = &mesh.verts(vert.neiVertices[0]);
        }
        else if (vert.boundary() == Boundary::north)
        {
            stencilVert = &mesh.verts(vert.neiVertices[1]);
        }
        else if (vert.boundary() == Boundary::southWest)
        {
            stencilVert = &mesh.verts(mesh.verts(vert.neiVertices[3]).neiVertices[2]);
        }
        else if (vert.boundary() == Boundary::southEast)
        {
            stencilVert = &mesh.verts(mesh.verts(vert.neiVertices[3]).neiVertices[0]);
        }
        else if (vert.boundary() == Boundary::northEast)
        {
            stencilVert = &mesh.verts(mesh.verts(vert.neiVertices[1]).neiVertices[0]);
        }
        else if (vert.boundary() == Boundary::northWest)
        {
            stencilVert = &mesh.verts(mesh.verts(vert.neiVertices[1]).neiVertices[2]);
        }
        Cell* cells[4];
        for (int j = 0; j < 4; j++)
        {
            if (stencilVert->vertToCell[j] != -1)
            {
                cells[j] = &mesh.cells(stencilVert->vertToCell[j]);
            }
        }

        findRelated(mesh, cells, *stencilVert, values, coords);

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
void VEdgeField<T>::findRelated(const Mesh& mesh, Cell** cells, const Vert& vert, vector<double>& values, vector<Vector2D>& coords) const
{
    if (vert.type() == TVertType::non)
    {
        values.push_back(this->value[mesh.vEdges(cells[0]->cellToVEdge[0][0]).index()]);
        coords.push_back(mesh.vEdges(cells[0]->cellToVEdge[0][0]).midpointCoord());
        values.push_back(this->value[mesh.vEdges(cells[0]->cellToVEdge[1][0]).index()]);
        coords.push_back(mesh.vEdges(cells[0]->cellToVEdge[1][0]).midpointCoord());

        values.push_back(this->value[mesh.vEdges(cells[2]->cellToVEdge[0][0]).index()]);
        coords.push_back(mesh.vEdges(cells[2]->cellToVEdge[0][0]).midpointCoord());
        values.push_back(this->value[mesh.vEdges(cells[2]->cellToVEdge[1][0]).index()]);
        coords.push_back(mesh.vEdges(cells[2]->cellToVEdge[1][0]).midpointCoord());
    }
    else if (vert.type() == TVertType::west)
    {
        values.push_back(this->value[mesh.vEdges(cells[0]->cellToVEdge[0][0]).index()]);
        coords.push_back(mesh.vEdges(cells[0]->cellToVEdge[0][0]).midpointCoord());

        values.push_back(this->value[mesh.vEdges(cells[1]->cellToVEdge[0][0]).index()]);
        coords.push_back(mesh.vEdges(cells[1]->cellToVEdge[0][0]).midpointCoord());
        values.push_back(this->value[mesh.vEdges(cells[1]->cellToVEdge[1][0]).index()]);
        coords.push_back(mesh.vEdges(cells[1]->cellToVEdge[1][0]).midpointCoord());

        values.push_back(this->value[mesh.vEdges(cells[2]->cellToVEdge[0][0]).index()]);
        coords.push_back(mesh.vEdges(cells[2]->cellToVEdge[0][0]).midpointCoord());
        values.push_back(this->value[mesh.vEdges(cells[2]->cellToVEdge[1][0]).index()]);
        coords.push_back(mesh.vEdges(cells[2]->cellToVEdge[1][0]).midpointCoord());
    }
    else if (vert.type() == TVertType::south)
    {
        values.push_back(this->value[mesh.vEdges(cells[1]->cellToVEdge[0][0]).index()]);
        coords.push_back(mesh.vEdges(cells[1]->cellToVEdge[0][0]).midpointCoord());
        values.push_back(this->value[mesh.vEdges(cells[1]->cellToVEdge[1][0]).index()]);
        coords.push_back(mesh.vEdges(cells[1]->cellToVEdge[1][0]).midpointCoord());

        values.push_back(this->value[mesh.vEdges(cells[2]->cellToVEdge[0][0]).index()]);
        coords.push_back(mesh.vEdges(cells[2]->cellToVEdge[0][0]).midpointCoord());
        values.push_back(this->value[mesh.vEdges(cells[2]->cellToVEdge[1][0]).index()]);
        coords.push_back(mesh.vEdges(cells[2]->cellToVEdge[1][0]).midpointCoord());

        values.push_back(this->value[mesh.vEdges(cells[3]->cellToVEdge[0][0]).index()]);
        coords.push_back(mesh.vEdges(cells[3]->cellToVEdge[0][0]).midpointCoord());
    }
    else if (vert.type() == TVertType::east)
    {
        values.push_back(this->value[mesh.vEdges(cells[0]->cellToVEdge[1][0]).index()]);
        coords.push_back(mesh.vEdges(cells[0]->cellToVEdge[1][0]).midpointCoord());
        values.push_back(this->value[mesh.vEdges(cells[0]->cellToVEdge[0][0]).index()]);
        coords.push_back(mesh.vEdges(cells[0]->cellToVEdge[0][0]).midpointCoord());

        values.push_back(this->value[mesh.vEdges(cells[3]->cellToVEdge[0][0]).index()]);
        coords.push_back(mesh.vEdges(cells[3]->cellToVEdge[0][0]).midpointCoord());
        values.push_back(this->value[mesh.vEdges(cells[3]->cellToVEdge[1][0]).index()]);
        coords.push_back(mesh.vEdges(cells[3]->cellToVEdge[1][0]).midpointCoord());

        values.push_back(this->value[mesh.vEdges(cells[2]->cellToVEdge[1][0]).index()]);
        coords.push_back(mesh.vEdges(cells[2]->cellToVEdge[1][0]).midpointCoord());
    }
    else if (vert.type() == TVertType::north)
    {
        values.push_back(this->value[mesh.vEdges(cells[0]->cellToVEdge[0][0]).index()]);
        coords.push_back(mesh.vEdges(cells[0]->cellToVEdge[0][0]).midpointCoord());
        values.push_back(this->value[mesh.vEdges(cells[0]->cellToVEdge[1][0]).index()]);
        coords.push_back(mesh.vEdges(cells[0]->cellToVEdge[1][0]).midpointCoord());

        values.push_back(this->value[mesh.vEdges(cells[1]->cellToVEdge[1][0]).index()]);
        coords.push_back(mesh.vEdges(cells[1]->cellToVEdge[1][0]).midpointCoord());

        values.push_back(this->value[mesh.vEdges(cells[2]->cellToVEdge[0][0]).index()]);
        coords.push_back(mesh.vEdges(cells[2]->cellToVEdge[0][0]).midpointCoord());
        values.push_back(this->value[mesh.vEdges(cells[3]->cellToVEdge[1][0]).index()]);
        coords.push_back(mesh.vEdges(cells[3]->cellToVEdge[1][0]).midpointCoord());
    }
}

template<typename T>
template<typename Mesh>
void VEdgeField<T>::output(const Mesh& mesh, const string& _name) const
{
    // VEdgeField插值到VertField
    auto vertField = toVertField(mesh);

    // 输出到文件
    vertField->output(mesh, _name);
}
