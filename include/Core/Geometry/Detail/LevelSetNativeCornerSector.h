#pragma once

#include "Core/Geometry/Detail/LevelSetHermiteCell.h"
#include "Core/Mesh/CartesianLeafLocator.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <algorithm>
#include <array>
#include <iterator>

namespace Quantum::Geometry::Detail
{
namespace
{
namespace LevelSetHermiteScalar
{
// Original native reentrant-corner topology/trace rules, independent of
// polynomial degree. The scalar owner checks uniqueness, exact zero/value,
// marked tangential components and its OWN full-degree phase certificate.
// No contact-angle model, local remeshing or relaxed activity tolerance here.
template<class Active>
std::array<bool,4> nativeCornerSharedComponents(const QuadtreeMesh& mesh,
    const CartesianLeafLocator& locator,int cellIndex,int vertexIndex,const Active& active)
{
    Detail::require(cellIndex>=0 && cellIndex<mesh.cellIndexNum() && active(cellIndex)
        && vertexIndex>=0 && vertexIndex<mesh.vertIndexNum(),
        "native corner sector requires active cell and valid native vertex");
    const auto& vertex=mesh.verts(mesh.vertIndexToLabel[vertexIndex]);
    const auto point=vertex.coord();
    Detail::require(point.x>mesh.minCoord.x && point.x<mesh.maxCoord.x
        && point.y>mesh.minCoord.y && point.y<mesh.maxCoord.y && vertex.supCell<0,
        "native corner sector requires an internal regular vertex");
    const auto incident=locator.incident(point);
    Detail::require(incident.count==4,"native corner sector requires four actual corner cells");
    const auto& native=mesh.cells(mesh.cellIndexToLabel[cellIndex]);
    int inactive=-1,activeCount=0;
    for(int k=0;k<incident.count;++k)
    {
        const int other=incident.indices[k];
        const auto& neighbor=mesh.cells(mesh.cellIndexToLabel[other]);
        const bool corner=std::any_of(std::begin(neighbor.cellToVert),std::end(neighbor.cellToVert),
            [&](int label){return mesh.verts(label).index()==vertexIndex;});
        Detail::require(neighbor.level()==native.level() && corner,
            "native corner sector requires a locally same-level boundary band");
        if(active(other))++activeCount;else inactive=other;
    }
    Detail::require(activeCount==3,"native corner sector requires three active cells and one solid cell");
    bool touchesSolid=false;
    for(int side=0;side<4;++side)for(int part=0;part<2;++part)
        touchesSolid=touchesSolid || native.neiCell[side][part]==mesh.cellIndexToLabel[inactive];
    Detail::require(touchesSolid,"native corner sector must share an edge with the inactive cell");
    std::array<bool,4> shared{{true,false,false,false}};
    for(int side=0;side<4;++side)for(int part=0;part<2;++part)
    {
        const int label=native.neiCell[side][part];
        if(label<0)continue;
        const auto& neighbor=mesh.cells(label);
        if(!active(neighbor.index()))continue;
        const bool corner=std::any_of(std::begin(neighbor.cellToVert),std::end(neighbor.cellToVert),
            [&](int v){return mesh.verts(v).index()==vertexIndex;});
        if(corner)shared[side==0 || side==2?2:1]=true;
    }
    return shared;
}
}
}
}
