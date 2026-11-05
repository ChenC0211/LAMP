#pragma once

#include "Core/Base/Enumeration.h"
#include "Core/Base/Interpolation.h"
#include "Core/Field/FieldTransferValidation.h"
using namespace std;

template<typename T>
HEdgeField<T>::HEdgeField() : Field<T>() {}

template<typename T>
HEdgeField<T>::HEdgeField(const HEdgeField<T>& field) : Field<T>(field) {}

template<typename T>
template<typename Mesh>
HEdgeField<T>::HEdgeField(Mesh& mesh, const T& value) : Field<T>()
{
    this->size_ = mesh.hEdgeIndexNum();
    this->value = std::make_unique<T[]>(this->size_);
    this->label_ = std::make_unique<int[]>(this->size_);
#pragma omp parallel for
    for (int i = 0; i < this->size_; i++)
    {
        Edge& hEdge = mesh.hEdges(mesh.hEdgeIndexToLabel[i]);

        this->value[i] = value;
        this->label_[i] = hEdge.label();
    }
}

template<typename T>
HEdgeField<T>::~HEdgeField() {}

template<typename T>
template<typename Mesh>
void HEdgeField<T>::initialize(const Mesh& mesh, const T& value)
{
    this->size_ = mesh.hEdgeIndexNum();
    this->value = std::make_unique<T[]>(this->size_);
    this->label_ = std::make_unique<int[]>(this->size_);
#pragma omp parallel for
    for (int i = 0; i < this->size_; i++)
    {
        Edge& hEdge = mesh.hEdges(mesh.hEdgeIndexToLabel[i]);

        this->value[i] = value;
        this->label_[i] = hEdge.label();
    }
}

template<typename T>
template<typename Mesh>
void HEdgeField<T>::updateWithMesh(const Mesh& mesh)
{
    auto tempValue = std::make_unique<T[]>(mesh.hEdgeIndexNum());
    auto tempLabel = std::make_unique<int[]>(mesh.hEdgeIndexNum());

    int _ptr1 = 0;
    int _ptr2 = 0;
    while (_ptr1 < this->size_ && _ptr2 < mesh.hEdgeIndexNum())
    {
        if (this->label_[_ptr1] == mesh.hEdges(mesh.hEdgeIndexToLabel[_ptr2]).label())
        {
            mesh.hEdges(mesh.hEdgeIndexToLabel[_ptr2]).setNeedUpdate(false);
            tempValue[_ptr2] = this->value[_ptr1];
            tempLabel[_ptr2] = this->label_[_ptr1];
            _ptr1++;
            _ptr2++;
        }
        else if (this->label_[_ptr1] < mesh.hEdges(mesh.hEdgeIndexToLabel[_ptr2]).label())
        {
            _ptr1++;
        }
        else
        {
            mesh.hEdges(mesh.hEdgeIndexToLabel[_ptr2]).setNeedUpdate(true);
            tempLabel[_ptr2] = mesh.hEdges(mesh.hEdgeIndexToLabel[_ptr2]).label();
            _ptr2++;
        }
    }

    while (_ptr2 < mesh.hEdgeIndexNum())
    {
        mesh.hEdges(mesh.hEdgeIndexToLabel[_ptr2]).setNeedUpdate(true);
        tempLabel[_ptr2] = mesh.hEdges(mesh.hEdgeIndexToLabel[_ptr2]).label();
        _ptr2++;
    }

    for (int i = 0; i < mesh.hEdgeIndexNum(); ++i)
    {
        const Edge& edge = mesh.hEdges(mesh.hEdgeIndexToLabel[i]);
        if (edge.needUpdate()) LAMP::FieldTransfer::validateSources(edge, this->size_, "horizontal edge");
    }
#pragma omp parallel for
    for (int i = 0; i < mesh.hEdgeIndexNum(); i++)
    {
        Edge& hEdge = mesh.hEdges(mesh.hEdgeIndexToLabel[i]);
        if (hEdge.needUpdate())
        {
            hEdge.setNeedUpdate(false);

            vector<double> values(hEdge.relatedPreIndex.size());
            for (int j = 0; j < hEdge.relatedPreIndex.size(); j++)
            {
                values[j] = this->value[hEdge.relatedPreIndex[j]];
            }

            // Robustness: leastSquareInterpolation (first-order) requires >=3 points.
            // Near boundaries / coarse-fine transitions the stencil can be smaller.
            if (hEdge.relatedPreCoord.size() < 3 || hEdge.relatedPreCoord.size() != values.size())
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
                    hEdge.relatedPreCoord,
                    values,
                    hEdge.midpointCoord(),
                    false
                );
            }
        }
    }

    this->value = std::move(tempValue);
    this->label_ = std::move(tempLabel);
    this->size_ = mesh.hEdgeIndexNum();
}

// template<typename T>
// template<typename Mesh>
// void HEdgeField<T>::subInterpolation(const Mesh& mesh, std::unique_ptr<T[]>& values, unique_ptr<int[]>& labels, int preI, Edge& hEdge)
// {
//     if (hEdge.son[0] != -1)
//     {
//         for (int i = 0; i < 2; i++)
//         {
//             subInterpolation(mesh, values, labels, preI, mesh.hEdges(hEdge.son[i]));
//         }
//     }
//     else if (hEdge.index() != -1)
//     {
//         values[hEdge.index()] = Interpolation::averageInterpolation
//             (
//                 mesh, 
//                 *this, 
//                 preI
//             );
//         labels[hEdge.index()] = hEdge.label();
//     }
// }

template<typename T>
template<typename Mesh>
unique_ptr<VertField<T>> HEdgeField<T>::toVertField(const Mesh& mesh) const
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

        // Robustness: handle degenerate stencils (e.g., boundary/corner) without aborting.
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
void HEdgeField<T>::findRelated(const Mesh& mesh, Cell** cells, const Vert& vert, vector<double>& values, vector<Vector2D>& coords) const
{
    if (vert.type() == TVertType::non)
    {
        values.push_back(this->value[mesh.hEdges(cells[0]->cellToHEdge[0][0]).index()]);
        coords.push_back(mesh.hEdges(cells[0]->cellToHEdge[0][0]).midpointCoord());
        values.push_back(this->value[mesh.hEdges(cells[0]->cellToHEdge[1][0]).index()]);
        coords.push_back(mesh.hEdges(cells[0]->cellToHEdge[1][0]).midpointCoord());

        values.push_back(this->value[mesh.hEdges(cells[2]->cellToHEdge[0][0]).index()]);
        coords.push_back(mesh.hEdges(cells[2]->cellToHEdge[0][0]).midpointCoord());
        values.push_back(this->value[mesh.hEdges(cells[2]->cellToHEdge[1][0]).index()]);
        coords.push_back(mesh.hEdges(cells[2]->cellToHEdge[1][0]).midpointCoord());
    }
    else if (vert.type() == TVertType::west)
    {
        values.push_back(this->value[mesh.hEdges(cells[0]->cellToHEdge[0][0]).index()]);
        coords.push_back(mesh.hEdges(cells[0]->cellToHEdge[0][0]).midpointCoord());
        values.push_back(this->value[mesh.hEdges(cells[0]->cellToHEdge[1][0]).index()]);
        coords.push_back(mesh.hEdges(cells[0]->cellToHEdge[1][0]).midpointCoord());

        values.push_back(this->value[mesh.hEdges(cells[1]->cellToHEdge[0][0]).index()]);
        coords.push_back(mesh.hEdges(cells[1]->cellToHEdge[0][0]).midpointCoord());
        values.push_back(this->value[mesh.hEdges(cells[1]->cellToHEdge[1][0]).index()]);
        coords.push_back(mesh.hEdges(cells[1]->cellToHEdge[1][0]).midpointCoord());

        values.push_back(this->value[mesh.hEdges(cells[2]->cellToHEdge[1][0]).index()]);
        coords.push_back(mesh.hEdges(cells[2]->cellToHEdge[1][0]).midpointCoord());
    }
    else if (vert.type() == TVertType::south)
    {
        values.push_back(this->value[mesh.hEdges(cells[1]->cellToHEdge[0][0]).index()]);
        coords.push_back(mesh.hEdges(cells[1]->cellToHEdge[0][0]).midpointCoord());

        values.push_back(this->value[mesh.hEdges(cells[2]->cellToHEdge[0][0]).index()]);
        coords.push_back(mesh.hEdges(cells[2]->cellToHEdge[0][0]).midpointCoord());
        values.push_back(this->value[mesh.hEdges(cells[2]->cellToHEdge[1][0]).index()]);
        coords.push_back(mesh.hEdges(cells[2]->cellToHEdge[1][0]).midpointCoord());

        values.push_back(this->value[mesh.hEdges(cells[3]->cellToHEdge[0][0]).index()]);
        coords.push_back(mesh.hEdges(cells[3]->cellToHEdge[0][0]).midpointCoord());
        values.push_back(this->value[mesh.hEdges(cells[3]->cellToHEdge[1][0]).index()]);
        coords.push_back(mesh.hEdges(cells[3]->cellToHEdge[1][0]).midpointCoord());
    }
    else if (vert.type() == TVertType::east)
    {
        values.push_back(this->value[mesh.hEdges(cells[0]->cellToHEdge[1][0]).index()]);
        coords.push_back(mesh.hEdges(cells[0]->cellToHEdge[1][0]).midpointCoord());
        values.push_back(this->value[mesh.hEdges(cells[0]->cellToHEdge[0][0]).index()]);
        coords.push_back(mesh.hEdges(cells[0]->cellToHEdge[0][0]).midpointCoord());

        values.push_back(this->value[mesh.hEdges(cells[3]->cellToHEdge[1][0]).index()]);
        coords.push_back(mesh.hEdges(cells[3]->cellToHEdge[1][0]).midpointCoord());

        values.push_back(this->value[mesh.hEdges(cells[2]->cellToHEdge[0][0]).index()]);
        coords.push_back(mesh.hEdges(cells[2]->cellToHEdge[0][0]).midpointCoord());
        values.push_back(this->value[mesh.hEdges(cells[2]->cellToHEdge[1][0]).index()]);
        coords.push_back(mesh.hEdges(cells[2]->cellToHEdge[1][0]).midpointCoord());
    }
    else if (vert.type() == TVertType::north)
    {
        values.push_back(this->value[mesh.hEdges(cells[0]->cellToHEdge[0][0]).index()]);
        coords.push_back(mesh.hEdges(cells[0]->cellToHEdge[0][0]).midpointCoord());
        values.push_back(this->value[mesh.hEdges(cells[0]->cellToHEdge[1][0]).index()]);
        coords.push_back(mesh.hEdges(cells[0]->cellToHEdge[1][0]).midpointCoord());

        values.push_back(this->value[mesh.hEdges(cells[1]->cellToHEdge[1][0]).index()]);
        coords.push_back(mesh.hEdges(cells[1]->cellToHEdge[1][0]).midpointCoord());
        values.push_back(this->value[mesh.hEdges(cells[1]->cellToHEdge[0][0]).index()]);
        coords.push_back(mesh.hEdges(cells[1]->cellToHEdge[0][0]).midpointCoord());

        values.push_back(this->value[mesh.hEdges(cells[3]->cellToHEdge[1][0]).index()]);
        coords.push_back(mesh.hEdges(cells[3]->cellToHEdge[1][0]).midpointCoord());
    }
}

template<typename T>
template<typename Mesh>
void HEdgeField<T>::output(const Mesh& mesh, const string& _name) const
{
    // HEdgeField插值到VertField
    auto vertField = toVertField(mesh);

    // 输出到文件
    vertField->output(mesh, _name);
}
