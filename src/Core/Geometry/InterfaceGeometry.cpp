#include "Core/Geometry/InterfaceGeometry.h"
#include "Core/Geometry/LevelSetPhaseFraction.h"
#include "Core/Field/VertField.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <array>
#include <cmath>
#include <cstring>
#include <stdexcept>

namespace {
void require(bool ok, const char* message)
{ if (!ok) throw std::invalid_argument(message); }
std::uint64_t bits(double value)
{
    static_assert(sizeof(double)==sizeof(std::uint64_t), "binary64 source words required");
    std::uint64_t result;std::memcpy(&result,&value,sizeof(result));return result;
}
std::vector<std::uint64_t> sourceWords(const QuadtreeMesh& mesh, const VertField<double>& phi,
    CartesianGridScope scope)
{
    const UniformLeafGrid grid(mesh,scope);
    require(phi.value && phi.size()==mesh.vertIndexNum(),"interface scalar storage must match mesh vertices");
    for(int i=0;i<mesh.vertIndexNum();++i)
        require(phi.label(i)==mesh.vertIndexToLabel[i] && std::isfinite(phi.value[i]),
            "interface scalar mapping or finite values do not match the mesh");
    std::vector<std::uint64_t> words;
    const auto integer = [&](int value) { words.push_back(static_cast<std::uint64_t>(static_cast<std::int64_t>(value))); };
    const auto real = [&](double value) { words.push_back(bits(value)); };
    const auto point = [&](Vector2D value) { real(value.x); real(value.y); };
    integer(mesh.Nx); integer(mesh.Ny); point(mesh.minCoord); point(mesh.maxCoord); point(mesh.cellSize);
    integer(mesh.cellLabelCapacity()); integer(mesh.vertLabelCapacity());
    integer(mesh.vEdgeLabelCapacity()); integer(mesh.hEdgeLabelCapacity());
    for (int k = 0; k < mesh.cellIndexNum(); ++k)
    {
        const int label = mesh.cellIndexToLabel[k]; const auto& cell = mesh.cells(label);
        integer(label); integer(cell.index()); integer(cell.label()); integer(cell.level());
        integer(cell.parent); for (int child : cell.son) integer(child);
        point(cell.centerCoord()); point(cell.size()); real(cell.area()); integer(static_cast<int>(cell.boundary()));
        for (int vertex : cell.cellToVert) integer(vertex);
        for (int side = 0; side < 2; ++side) for (int part = 0; part < 2; ++part)
        { integer(cell.cellToVEdge[side][part]); integer(cell.cellToHEdge[side][part]); }
        for (int side = 0; side < 4; ++side) for (int part = 0; part < 2; ++part) integer(cell.neiCell[side][part]);
    }
    for (int k = 0; k < mesh.vertIndexNum(); ++k)
    {
        const int label = mesh.vertIndexToLabel[k]; const auto& vertex = mesh.verts(label);
        integer(label); integer(vertex.index()); integer(vertex.label()); point(vertex.coord());
        integer(phi.label(k)); real(phi.value[k]);
    }
    const auto edgeWords = [&](const Edge& edge, int label)
    {
        integer(label); integer(edge.label()); integer(edge.index()); point(edge.midpointCoord()); real(edge.length());
        integer(static_cast<int>(edge.isHorizontal())); integer(static_cast<int>(edge.boundary()));
        integer(edge.masterCell); integer(edge.slaveCell); integer(edge.parent);
        for (int child : edge.son) integer(child);
        for (int vertex : edge.edgeToVert) integer(vertex);
    };
    for (int k = 0; k < mesh.vEdgeIndexNum(); ++k)
    { const int label = mesh.vEdgeIndexToLabel[k]; edgeWords(mesh.vEdges(label), label); }
    for (int k = 0; k < mesh.hEdgeIndexNum(); ++k)
    { const int label = mesh.hEdgeIndexToLabel[k]; edgeWords(mesh.hEdges(label), label); }
    return words;
}
}

InterfaceGeometry::InterfaceGeometry(const QuadtreeMesh& mesh, const VertField<double>& phi,
    bool axisymmetric, CartesianGridScope scope, LevelSetHermiteGeometry::Reconstruction reconstruction)
    : mesh_(&mesh), axisymmetric_(axisymmetric), scope_(scope),
      source_(std::make_shared<const LevelSetHermiteGeometry::Snapshot>(mesh,phi,
          axisymmetric && mesh.minCoord.x==0.0,scope,reconstruction)),
      words_(sourceWords(mesh,phi,scope))
{
    require(!axisymmetric || mesh.minCoord.x>=0.0,"axisymmetric interface requires nonnegative radius");
    const UniformLeafGrid grid(mesh,scope);
    fraction_.resize(mesh.cellIndexNum());
    for(int j=0;j<grid.ny();++j)for(int i=0;i<grid.nx();++i)
    {
        const auto& cell=mesh.cells(grid.cells()[i+grid.nx()*j]);
        const auto vertex=[&](int x,int y){return grid.vertices()[x+(grid.nx()+1)*y];};
        const int corners[4]={vertex(i,j),vertex(i+1,j),vertex(i+1,j+1),vertex(i,j+1)};
        std::array<double,4> values;
        for(int corner=0;corner<4;++corner) {
            require(cell.cellToVert[corner]==corners[corner],"interface cell corners disagree with uniform lattice");
            values[corner]=phi.value[mesh.verts(corners[corner]).index()];
        }
        fraction_[cell.index()]=LevelSetPhaseFraction::positiveFraction(values);
    }
}
void InterfaceGeometry::validateSource(const QuadtreeMesh& mesh, const VertField<double>& phi) const
{
    require(source_ && &mesh==mesh_,"interface geometry belongs to a different mesh or was moved from");
    source_->value(source_->minCoord()); // Preserve the guarded FP environment check.
    require(sourceWords(mesh,phi,scope_)==words_,"stale interface geometry: mesh or scalar changed");
}
