#include "Core/Mesh/CartesianLeafLevelIndex.h"
#include "Core/Mesh/CartesianLeafMesh.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <limits>
#include <stdexcept>

CartesianLeafLevelIndex::CartesianLeafLevelIndex(const QuadtreeMesh& mesh)
    :coordinates_(mesh.minCoord,mesh.maxCoord,mesh.Nx,mesh.Ny),nx_(mesh.Nx),ny_(mesh.Ny),addresses_(mesh.cellIndexNum())
{
    validateCartesianLeafMesh(mesh);
    for(int c=0;c<mesh.cellIndexNum();++c) {
        const auto a=LAMP::Mesh::cellAddress(mesh,mesh.cells(mesh.cellIndexToLabel[c]));
        const auto size=extent(a.level);
        if(a.i>=static_cast<unsigned>(size[0]) || a.j>=static_cast<unsigned>(size[1]))
            throw std::invalid_argument("leaf address lies outside background lattice");
        addresses_[c]={static_cast<int>(a.i),static_cast<int>(a.j),a.level};
        if(!cells_.emplace(std::array<int,3>{{static_cast<int>(a.level),static_cast<int>(a.i),static_cast<int>(a.j)}},c).second)
            throw std::invalid_argument("duplicate same-level leaf address");
    }
}
std::array<int,2> CartesianLeafLevelIndex::extent(unsigned level) const
{
    if(level>=31 || nx_>(std::numeric_limits<int>::max()>>level) || ny_>(std::numeric_limits<int>::max()>>level))
        throw std::invalid_argument("same-level lattice exceeds supported logical index range");
    return {{nx_<<level,ny_<<level}};
}
int CartesianLeafLevelIndex::cell(unsigned level,int i,int j) const
{
    const auto size=extent(level);
    if(i<0 || j<0 || i>=size[0] || j>=size[1]) return -1;
    const auto found=cells_.find({{static_cast<int>(level),i,j}});
    return found==cells_.end()?-1:found->second;
}
