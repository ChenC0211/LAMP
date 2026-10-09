#include "Core/IO/TecplotMeshOutput.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Field/VertField.h"

#include <stdexcept>
#include <utility>

namespace TecplotMeshOutput
{
TecplotWriter::Zone nodalZone(const QuadtreeMesh& mesh,
    const std::vector<std::string>& names, const std::vector<VertField<double>>& fields)
{
    if (names.size() != fields.size())
        throw std::invalid_argument("TecplotMeshOutput: field/name count mismatch");
    const int vertices = mesh.vertIndexNum(), cells = mesh.cellIndexNum();
    if (vertices <= 0 || cells <= 0 || !mesh.vertIndexToLabel || !mesh.cellIndexToLabel)
        throw std::invalid_argument("TecplotMeshOutput: missing indexed mesh");
    TecplotWriter::Zone zone;
    zone.title = "Quadtree Mesh";
    zone.coordinateNames = {"X", "Y"};
    zone.element = TecplotWriter::Element::Quadrilateral;
    zone.points.reserve(vertices);
    for (int i = 0; i < vertices; ++i)
    {
        const int label = mesh.vertIndexToLabel[i];
        if (label < 0 || label >= mesh.vertLabelCapacity() || mesh.verts(label).index() != i)
            throw std::invalid_argument("TecplotMeshOutput: invalid vertex index mapping");
        zone.points.push_back(mesh.verts(label).coord());
    }
    zone.connectivity.reserve(cells);
    for (int i = 0; i < cells; ++i)
    {
        const int label = mesh.cellIndexToLabel[i];
        if (label < 0 || label >= mesh.cellLabelCapacity() || mesh.cells(label).index() != i)
            throw std::invalid_argument("TecplotMeshOutput: invalid cell index mapping");
        const auto& cell = mesh.cells(label);
        std::vector<int> connection;
        connection.reserve(4);
        for (int corner = 0; corner < 4; ++corner)
        {
            const int vertexLabel = cell.cellToVert[corner];
            if (vertexLabel < 0 || vertexLabel >= mesh.vertLabelCapacity())
                throw std::invalid_argument("TecplotMeshOutput: invalid cell vertex label");
            const int index = mesh.verts(vertexLabel).index();
            if (index < 0 || index >= vertices || mesh.vertIndexToLabel[index] != vertexLabel)
                throw std::invalid_argument("TecplotMeshOutput: unindexed cell vertex");
            connection.push_back(index);
        }
        zone.connectivity.push_back(std::move(connection));
    }
    zone.fields.reserve(fields.size());
    for (std::size_t f = 0; f < fields.size(); ++f)
    {
        if (fields[f].size() != vertices || !fields[f].value)
            throw std::invalid_argument("TecplotMeshOutput: nodal field size mismatch");
        TecplotWriter::Field field;
        field.name = names[f];
        field.location = TecplotWriter::Location::Nodal;
        field.values.assign(fields[f].value.get(), fields[f].value.get() + vertices);
        zone.fields.push_back(std::move(field));
    }
    return zone;
}
}
