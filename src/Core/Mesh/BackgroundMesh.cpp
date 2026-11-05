#include "Core/Mesh/BackgroundMesh.h"
#include "Core/Base/Base.h"
#include "Core/Mesh/MeshCacheIO.h"
#include "Core/Mesh/DyadicGridCoordinates.h"
#include <algorithm>
#include <cstdint>

using namespace std;

BackgroundMesh::BackgroundMesh
(
    Vector2D _minCoord, 
    Vector2D _maxCoord, 
    int _Nx, 
    int _Ny
) :
    minCoord(_minCoord), 
    maxCoord(_maxCoord), 
    Nx(_Nx), 
    Ny(_Ny),
    cellSize
    (
        (_maxCoord.x - _minCoord.x) / static_cast<double>(_Nx), 
        (_maxCoord.y - _minCoord.y) / static_cast<double>(_Ny)
    )
{
    levelSetCount_ = 1;
    vertIndexNum_ = (Nx + 1) * (Ny + 1);
    hEdgeIndexNum_ = Nx * (Ny + 1);
    vEdgeIndexNum_ = Ny * (Nx + 1);
    cellIndexNum_ = Nx * Ny;
    vertLabelCapacity_ = vertIndexNum_;
    hEdgeLabelCapacity_ = hEdgeIndexNum_;
    vEdgeLabelCapacity_ = vEdgeIndexNum_;
    cellLabelCapacity_ = cellIndexNum_;

    verts_.resetDense(vertIndexNum_);
    hEdges_.resetDense(hEdgeIndexNum_);
    vEdges_.resetDense(vEdgeIndexNum_);
    cells_.resetDense(cellIndexNum_);

    vertIndexToLabel = new int[vertIndexNum_];
    hEdgeIndexToLabel = new int[hEdgeIndexNum_];
    vEdgeIndexToLabel = new int[vEdgeIndexNum_];
    cellIndexToLabel = new int[cellIndexNum_];

    buildMesh();
    ensureLevelSetStorage();
}

BackgroundMesh::BackgroundMesh(const DictionaryMesh& dictMesh) :
    minCoord(dictMesh.minCoord), 
    maxCoord(dictMesh.maxCoord), 
    Nx(dictMesh.Nx), 
    Ny(dictMesh.Ny),
    cellSize
    (
        (dictMesh.maxCoord.x - dictMesh.minCoord.x) / static_cast<double>(dictMesh.Nx), 
        (dictMesh.maxCoord.y - dictMesh.minCoord.y) / static_cast<double>(dictMesh.Ny)
    )
{
    levelSetCount_ = dictMesh.levelSetCount == 0 ? 0 : std::max(1, dictMesh.levelSetCount);
    vertIndexNum_ = (Nx + 1) * (Ny + 1);
    hEdgeIndexNum_ = Nx * (Ny + 1);
    vEdgeIndexNum_ = Ny * (Nx + 1);
    cellIndexNum_ = Nx * Ny;
    vertLabelCapacity_ = vertIndexNum_;
    hEdgeLabelCapacity_ = hEdgeIndexNum_;
    vEdgeLabelCapacity_ = vEdgeIndexNum_;
    cellLabelCapacity_ = cellIndexNum_;

    verts_.resetDense(vertIndexNum_);
    hEdges_.resetDense(hEdgeIndexNum_);
    vEdges_.resetDense(vEdgeIndexNum_);
    cells_.resetDense(cellIndexNum_);

    vertIndexToLabel = new int[vertIndexNum_];
    hEdgeIndexToLabel = new int[hEdgeIndexNum_];
    vEdgeIndexToLabel = new int[vEdgeIndexNum_];
    cellIndexToLabel = new int[cellIndexNum_];
    
    buildMesh();
    ensureLevelSetStorage();
}

BackgroundMesh::BackgroundMesh(const BackgroundMesh& source, CopyMetadataTag)
    : minCoord(source.minCoord), maxCoord(source.maxCoord), Nx(source.Nx), Ny(source.Ny),
      cellSize(source.cellSize), levelSetCount_(source.levelSetCount_)
{
}

BackgroundMesh::~BackgroundMesh() 
{

    delete[] vertIndexToLabel;
    delete[] hEdgeIndexToLabel;
    delete[] vEdgeIndexToLabel;
    delete[] cellIndexToLabel;
}

void BackgroundMesh::buildMesh()
{
    // 生成顶点
    createVertices();

    // 生成边
    createEdges();

    // 生成单元
    createCells();

    // 查找相邻单元
    findNeighbourCells();

    // 定义类型
    defineTypes();

    statistics(verts_, vertIndexNum_);
    statistics(hEdges_, hEdgeIndexNum_);
    statistics(vEdges_, vEdgeIndexNum_);
    statistics(cells_, cellIndexNum_);
}

void BackgroundMesh::ensureLevelSetStorage()
{
    const std::size_t n = static_cast<std::size_t>(levelSetCount_ == 0 ? 0 : std::max(1, levelSetCount_));
    for (int i = 0; i < vertIndexNum_; ++i)
        verts_[i].ensureLevelSetSize(n);
    for (int i = 0; i < hEdgeIndexNum_; ++i)
        hEdges_[i].ensureLevelSetSize(n);
    for (int i = 0; i < vEdgeIndexNum_; ++i)
        vEdges_[i].ensureLevelSetSize(n);
    for (int i = 0; i < cellIndexNum_; ++i)
        cells_[i].ensureLevelSetSize(n);
}

void BackgroundMesh::createVertices()
{
    const LAMP::Mesh::DyadicGridCoordinates coordinates(minCoord,maxCoord,Nx,Ny);
    int count = 0;
    for (int j = 0; j < Ny + 1; j++)
    {
        for (int i = 0; i < Nx + 1; i++)
        {
            vertIndexToLabel[count] = count;// 顶点编号，后续动网格顶点编号变化，仅存储存在的顶点编号
            const int vertI = vertIndexToLabel[count];

            const Vector2D coord=coordinates.point(i,j,0);

            Vert vert(count, coord);
            verts_[vertI] = vert;

            count++;
        }
    }

    if (count != vertIndexNum_)
    {
        cerr << "Vertex count mismatch!" << endl;
        exit(-1);
    }
}

void BackgroundMesh::createEdges()
{
    const LAMP::Mesh::DyadicGridCoordinates coordinates(minCoord,maxCoord,Nx,Ny);
    int hCount = 0;
    // 水平边
    for (int j = 0; j < Ny + 1; j++)
    {
        for (int i = 0; i < Nx; i++)
        {
            hEdgeIndexToLabel[hCount] = hCount;// 边编号，后续动网格边编号变化，仅存储存在的边编号
            const int hEdgeI = hEdgeIndexToLabel[hCount];
            
            int vertI0 = j * (Nx + 1) + i;
            int vertI1 = j * (Nx + 1) + i + 1;
            Vector2D midpointCoord;
            midpointCoord = coordinates.point(2ULL*i+1,2ULL*j,1);

            Edge hEdge(hCount, midpointCoord, cellSize.x, true);
            hEdges_[hEdgeI] = hEdge;

            // 边到点的关系
            hEdges_[hEdgeI].edgeToVert[0] = vertI0;
            hEdges_[hEdgeI].edgeToVert[1] = vertI1;
            // 点到边的关系
            verts_[vertI0].vertToHEdge[1] = hEdgeI;
            verts_[vertI1].vertToHEdge[0] = hEdgeI;
            // 点到点的关系
            verts_[vertI0].neiVertices[2] = vertI1;
            verts_[vertI1].neiVertices[0] = vertI0;

            hCount++;
        }
    }

    int vCount = 0;
    // 垂直边
    for (int j = 0; j < Ny; j++)
    {
        for (int i = 0; i < Nx + 1; i++)
        {
            vEdgeIndexToLabel[vCount] = vCount;// 边编号，后续动网格边编号变化，仅存储存在的边编号
            const int vEdgeI = vEdgeIndexToLabel[vCount];

            int vertI0 = j * (Nx + 1) + i;
            int vertI1 = j * (Nx + 1) + i + Nx + 1;
            Vector2D midpointCoord;
            midpointCoord = coordinates.point(2ULL*i,2ULL*j+1,1);

            Edge vEdge(vCount, midpointCoord, cellSize.y, false);
            vEdges_[vEdgeI] = vEdge;

            // 边到点的关系
            vEdges_[vEdgeI].edgeToVert[0] = vertI0;
            vEdges_[vEdgeI].edgeToVert[1] = vertI1;
            // 点到边的关系
            verts_[vertI0].vertToVEdge[1] = vEdgeI;
            verts_[vertI1].vertToVEdge[0] = vEdgeI;
            // 点到点的关系
            verts_[vertI0].neiVertices[3] = vertI1;
            verts_[vertI1].neiVertices[1] = vertI0;

            vCount++;
        }
    }

    if (hCount != hEdgeIndexNum_ || vCount != vEdgeIndexNum_)
    {
        cerr << "Edge count mismatch!" << endl;
        exit(-1);
    }
}

void BackgroundMesh::createCells()
{
    const LAMP::Mesh::DyadicGridCoordinates coordinates(minCoord,maxCoord,Nx,Ny);
    int count = 0;
    for (int j = 0; j < Ny; j++)
    {
        for (int i = 0; i < Nx; i++)
        {
            cellIndexToLabel[count] = count;// 单元编号，后续动网格单元编号变化，仅存储存在的单元编号
            const int cellI = cellIndexToLabel[count];
            
            int vertI0 = j * (Nx + 1) + i;
            int vertI1 = j * (Nx + 1) + i + 1;
            int vertI2 = (j + 1) * (Nx + 1) + i + 1;
            int vertI3 = (j + 1) * (Nx + 1) + i;

            Vector2D centerCoord;
            centerCoord = coordinates.point(2ULL*i+1,2ULL*j+1,1);

            Cell cell(count, centerCoord, 0, cellSize);
            cells_[cellI] = cell;

            // 单元到点的关系
            cells_[cellI].cellToVert[0] = vertI0;
            cells_[cellI].cellToVert[1] = vertI1;
            cells_[cellI].cellToVert[2] = vertI2;
            cells_[cellI].cellToVert[3] = vertI3;

            // 单元到边的关系
            int hEdgeIndex0 = j * Nx + i;
            int hEdgeIndex1 = (j + 1) * Nx + i;
            int vEdgeIndex0 = j * (Nx + 1) + i;
            int vEdgeIndex1 = j * (Nx + 1) + i + 1;

            cells_[cellI].cellToHEdge[0][0] = hEdgeIndex0;
            cells_[cellI].cellToHEdge[1][0] = hEdgeIndex1; 
            cells_[cellI].cellToVEdge[0][0] = vEdgeIndex0;
            cells_[cellI].cellToVEdge[1][0] = vEdgeIndex1;

            // 点到单元的关系
            verts_[vertI0].vertToCell[2] = cellI;
            verts_[vertI1].vertToCell[3] = cellI;
            verts_[vertI2].vertToCell[0] = cellI;
            verts_[vertI3].vertToCell[1] = cellI;

            // 边到单元的关系
            hEdges_[hEdgeIndex0].masterCell = cellI;
            hEdges_[hEdgeIndex1].slaveCell = cellI;
            vEdges_[vEdgeIndex0].masterCell = cellI;
            vEdges_[vEdgeIndex1].slaveCell = cellI;

            count++;
        }
    }
}

void BackgroundMesh::findNeighbourCells()
{
    for (int i = 0; i < hEdgeIndexNum_; i++)
    {
        const int edgeI = hEdgeIndexToLabel[i];
        Edge& edge = hEdges_[edgeI];
        if (edge.masterCell != -1 && edge.slaveCell != -1)
        {
            cells_[edge.masterCell].neiCell[1][0] = edge.slaveCell;
            cells_[edge.slaveCell].neiCell[3][0] = edge.masterCell;
        }
    }

    for (int i = 0; i < vEdgeIndexNum_; i++)
    {
        const int edgeI = vEdgeIndexToLabel[i];
        Edge& edge = vEdges_[edgeI];
        if (edge.masterCell != -1 && edge.slaveCell != -1)
        {
            cells_[edge.masterCell].neiCell[0][0] = edge.slaveCell;
            cells_[edge.slaveCell].neiCell[2][0] = edge.masterCell;
        }
    }
}

void BackgroundMesh::defineTypes()
{
    // cell
    for (int i = 0; i < cellIndexNum_; i++)
    {
        Cell& cell = cells_[cellIndexToLabel[i]];
        cell.setBoundary(Boundary::non);
        if (cell.neiCell[0][0] == -1)
        {
            cell.setBoundary(Boundary::west);
        }

        if (cell.neiCell[2][0] == -1)
        {
            cell.setBoundary(Boundary::east);
        }

        if (cell.neiCell[1][0] == -1)
        {
            if (cell.boundary() == Boundary::west)
            {
                cell.setBoundary(Boundary::southWest);
            }
            else if (cell.boundary() == Boundary::east)
            {
                cell.setBoundary(Boundary::southEast);
            }
            else
            {
                cell.setBoundary(Boundary::south);
            }
        }

        if (cell.neiCell[3][0] == -1)
        {
            if (cell.boundary() == Boundary::west)
            {
                cell.setBoundary(Boundary::northWest);
            }
            else if (cell.boundary() == Boundary::east)
            {
                cell.setBoundary(Boundary::northEast);
            }
            else
            {
                cell.setBoundary(Boundary::north);
            }
        }
    }

    // hEdge
    for (int i = 0; i < hEdgeIndexNum_; i++)
    {
        Edge& edge = hEdges_[hEdgeIndexToLabel[i]];
        edge.setBoundary(Boundary::non);
        if (edge.masterCell == -1)
        {
            edge.setBoundary(Boundary::north);
        }

        if (edge.slaveCell == -1)
        {
            edge.setBoundary(Boundary::south);
        }
    }

    // vEdge
    for (int i = 0; i < vEdgeIndexNum_; i++)
    {
        Edge& edge = vEdges_[vEdgeIndexToLabel[i]];
        edge.setBoundary(Boundary::non);
        if (edge.masterCell == -1)
        {
            edge.setBoundary(Boundary::east);
        }

        if (edge.slaveCell == -1)
        {
            edge.setBoundary(Boundary::west);
        }
    }

    // vert
    for (int i = 0; i < vertIndexNum_; i++)
    {
        Vert& vert = verts_[vertIndexToLabel[i]];
        vert.setBoundary(Boundary::non);
        vert.setType(TVertType::non);
        vert.supCell = -1;
        vert.supVert[0] = -1;
        vert.supVert[1] = -1;
        if (vert.vertToHEdge[0] == -1)
        {
            if (vert.vertToCell[0] == -1)
            {
                if (vert.vertToVEdge[0] == -1)
                {
                    vert.setBoundary(Boundary::southWest);
                }
                else if (vert.vertToVEdge[1] == -1)
                {
                    vert.setBoundary(Boundary::northWest);
                }
                else
                {
                    vert.setBoundary(Boundary::west);
                }
            }
            else
            {
                vert.setType(TVertType::west);
                vert.supCell = vert.vertToCell[0];
                vert.supVert[0] = cells_[vert.supCell].cellToVert[0];
                vert.supVert[1] = cells_[vert.supCell].cellToVert[3];
            }
        }

        if (vert.vertToHEdge[1] == -1)
        {
            if (vert.vertToCell[2] == -1)
            {
                if (vert.vertToVEdge[0] == -1)
                {
                    vert.setBoundary(Boundary::southEast);
                }
                else if (vert.vertToVEdge[1] == -1)
                {
                    vert.setBoundary(Boundary::northEast);
                }
                else
                {
                    vert.setBoundary(Boundary::east);
                }
            }
            else
            {
                vert.setType(TVertType::east);
                vert.supCell = vert.vertToCell[2];
                vert.supVert[0] = cells_[vert.supCell].cellToVert[1];
                vert.supVert[1] = cells_[vert.supCell].cellToVert[2];
            }
        }

        if (vert.vertToVEdge[0] == -1)
        {
            if (vert.vertToCell[1] == -1)
            {
                if (vert.boundary() == Boundary::non)
                {
                    vert.setBoundary(Boundary::south);
                }
            }
            else
            {
                vert.setType(TVertType::south);
                vert.supCell = vert.vertToCell[1];
                vert.supVert[0] = cells_[vert.supCell].cellToVert[0];
                vert.supVert[1] = cells_[vert.supCell].cellToVert[1];
            }
        }

        if (vert.vertToVEdge[1] == -1)
        {
            if (vert.vertToCell[3] == -1)
            {
                if (vert.boundary() == Boundary::non)
                {
                    vert.setBoundary(Boundary::north);
                }
            }
            else
            {
                vert.setType(TVertType::north);
                vert.supCell = vert.vertToCell[3];
                vert.supVert[0] = cells_[vert.supCell].cellToVert[3];
                vert.supVert[1] = cells_[vert.supCell].cellToVert[2];
            }
        }
    }
}

template<typename T>
void BackgroundMesh::statistics(LAMP::Mesh::Detail::EntityPool<T>& t, int num)
{
    int count = 0;
    for (int i = 0; i < num; i++)
    {
        if (t.peek(i).label() != -1)
        {
            t[i].setIndex(count);
            count++;
        }
    }
}

int BackgroundMesh::findCellLabel(const Vector2D& coord) const
{
    return findCellLabel(coord, MINVAL);
}

int BackgroundMesh::findCellLabel(const Vector2D& coord, const double boundaryTolerance) const
{
    // Safety: reject non-finite coordinates early to avoid undefined behavior
    // in comparisons/casts and to prevent out-of-range indexing.
    if (!std::isfinite(coord.x) || !std::isfinite(coord.y) ||
        !std::isfinite(boundaryTolerance) || boundaryTolerance < 0)
    {
        return -1;
    }

    if (coord.x < minCoord.x - boundaryTolerance || coord.x > maxCoord.x + boundaryTolerance ||
        coord.y < minCoord.y - boundaryTolerance || coord.y > maxCoord.y + boundaryTolerance)
    {
        return -1;
    }

    int i = coord.x < maxCoord.x - boundaryTolerance? static_cast<int>((coord.x - minCoord.x) / cellSize.x) : Nx - 1;
    int j = coord.y < maxCoord.y - boundaryTolerance? static_cast<int>((coord.y - minCoord.y) / cellSize.y) : Ny - 1;

    // Clamp to valid background range (defensive against roundoff).
    if (i < 0) i = 0;
    if (i >= Nx) i = Nx - 1;
    if (j < 0) j = 0;
    if (j >= Ny) j = Ny - 1;

    int cellLabel = j * Nx + i;

    if (cellLabel < 0)
        return -1;
    if (cellLabelCapacity_ > 0 && cellLabel >= cellLabelCapacity_)
        return -1;
    
    while (cells_[cellLabel].son[0] != -1)
    {
        const Vector2D cc = cells_[cellLabel].centerCoord();
        int next = -1;
        if (coord.x < cc.x)
        {
            next = (coord.y < cc.y) ? cells_[cellLabel].son[0] : cells_[cellLabel].son[3];
        }
        else
        {
            next = (coord.y < cc.y) ? cells_[cellLabel].son[1] : cells_[cellLabel].son[2];
        }

        if (next < 0)
            return -1;
        if (cellLabelCapacity_ > 0 && next >= cellLabelCapacity_)
            return -1;
        cellLabel = next;
    }

    return cellLabel;
}

int BackgroundMesh::findBackgroundCellLabel(const Vector2D& coord) const
{
    if (coord.x < minCoord.x - MINVAL || coord.x > maxCoord.x + MINVAL ||
        coord.y < minCoord.y - MINVAL || coord.y > maxCoord.y + MINVAL)
    {
        return -1;
    }

    int i = coord.x < maxCoord.x - MINVAL? static_cast<int>((coord.x - minCoord.x) / cellSize.x) : Nx - 1;
    int j = coord.y < maxCoord.y - MINVAL? static_cast<int>((coord.y - minCoord.y) / cellSize.y) : Ny - 1;

    int bgCellLabel = j * Nx + i;

    return bgCellLabel;
}


bool BackgroundMesh::isInsideCell(const int cellLabel, const Vector2D& coord, const double eps) const
{
    if (cellLabel < 0)
        return false;
    if (cellLabelCapacity_ > 0 && cellLabel >= cellLabelCapacity_)
        return false;

    const Cell& cell = cells_[cellLabel];
    if (cell.label() == -1)
        return false;

    const Vector2D c = cell.centerCoord();
    const Vector2D h(0.5 * cell.size().x, 0.5 * cell.size().y);

    return (coord.x >= c.x - h.x - eps) && (coord.x <= c.x + h.x + eps) &&
           (coord.y >= c.y - h.y - eps) && (coord.y <= c.y + h.y + eps);
}

int BackgroundMesh::findCellLabelFromHint(const int hintCellLabel, const Vector2D& coord, const double eps) const
{
    // Safety: avoid out-of-range indexing in findCellLabel for out-of-domain coords.
    if (coord.x < minCoord.x - eps || coord.x > maxCoord.x + eps ||
        coord.y < minCoord.y - eps || coord.y > maxCoord.y + eps)
    {
        return -1;
    }

    auto descendToLeaf = [&](int cellLabel) -> int
    {
        int c = cellLabel;
        if (c < 0)
            return -1;
        if (cellLabelCapacity_ > 0 && c >= cellLabelCapacity_)
            return -1;
        if (cells_[c].label() == -1)
            return -1;

        while (cells_[c].son[0] != -1)
        {
            const Vector2D cc = cells_[c].centerCoord();
            if (coord.x < cc.x)
            {
                c = (coord.y < cc.y) ? cells_[c].son[0] : cells_[c].son[3];
            }
            else
            {
                c = (coord.y < cc.y) ? cells_[c].son[1] : cells_[c].son[2];
            }

            if (c < 0)
                return -1;
            if (cellLabelCapacity_ > 0 && c >= cellLabelCapacity_)
                return -1;
            if (cells_[c].label() == -1)
                return -1;
        }
        return c;
    };

    // Hint must be within allocated cell label space.
    if (hintCellLabel >= 0 && cellLabelCapacity_ > 0 && hintCellLabel >= cellLabelCapacity_)
        return findCellLabel(coord);

    if (hintCellLabel >= 0 && cells_[hintCellLabel].label() != -1)
    {
        int hintLeaf = hintCellLabel;
        if (cells_[hintLeaf].son[0] != -1)
            hintLeaf = descendToLeaf(hintLeaf);

        if (isInsideCell(hintLeaf, coord, eps))
            return hintLeaf;

        // Enumerate a small candidate set (hint + 1-ring neighbors + diagonal corners),
        // then descend (if needed) and accept the first that contains coord.
        std::vector<int> cand;
        cand.reserve(16);
        auto pushUnique = [&](int v)
        {
            if (v < 0) return;
            if (cellLabelCapacity_ > 0 && v >= cellLabelCapacity_) return;
            if (cells_[v].label() == -1) return;
            for (int x : cand)
                if (x == v) return;
            cand.push_back(v);
        };

        pushUnique(hintLeaf);

        const Cell& cell = cells_[hintLeaf];
        for (int dir = 0; dir < 4; ++dir)
        {
            pushUnique(cell.neiCell[dir][0]);
            pushUnique(cell.neiCell[dir][1]);
        }

        // Diagonal corners via neighbors-of-neighbors.
        const int west[2]  = { cell.neiCell[0][0], cell.neiCell[0][1] };
        const int east[2]  = { cell.neiCell[2][0], cell.neiCell[2][1] };
        const int south[2] = { cell.neiCell[1][0], cell.neiCell[1][1] };
        const int north[2] = { cell.neiCell[3][0], cell.neiCell[3][1] };

        for (int k = 0; k < 2; ++k)
        {
            if (west[k] >= 0 && (cellLabelCapacity_ <= 0 || west[k] < cellLabelCapacity_) && cells_[west[k]].label() != -1)
            {
                pushUnique(cells_[west[k]].neiCell[3][0]);
                pushUnique(cells_[west[k]].neiCell[3][1]);
                pushUnique(cells_[west[k]].neiCell[1][0]);
                pushUnique(cells_[west[k]].neiCell[1][1]);
            }
            if (east[k] >= 0 && (cellLabelCapacity_ <= 0 || east[k] < cellLabelCapacity_) && cells_[east[k]].label() != -1)
            {
                pushUnique(cells_[east[k]].neiCell[3][0]);
                pushUnique(cells_[east[k]].neiCell[3][1]);
                pushUnique(cells_[east[k]].neiCell[1][0]);
                pushUnique(cells_[east[k]].neiCell[1][1]);
            }
            if (south[k] >= 0 && (cellLabelCapacity_ <= 0 || south[k] < cellLabelCapacity_) && cells_[south[k]].label() != -1)
            {
                pushUnique(cells_[south[k]].neiCell[0][0]);
                pushUnique(cells_[south[k]].neiCell[0][1]);
                pushUnique(cells_[south[k]].neiCell[2][0]);
                pushUnique(cells_[south[k]].neiCell[2][1]);
            }
            if (north[k] >= 0 && (cellLabelCapacity_ <= 0 || north[k] < cellLabelCapacity_) && cells_[north[k]].label() != -1)
            {
                pushUnique(cells_[north[k]].neiCell[0][0]);
                pushUnique(cells_[north[k]].neiCell[0][1]);
                pushUnique(cells_[north[k]].neiCell[2][0]);
                pushUnique(cells_[north[k]].neiCell[2][1]);
            }
        }

        for (int c : cand)
        {
            if (!isInsideCell(c, coord, eps))
                continue;
            return descendToLeaf(c);
        }
    }

    // Fallback to the generic descent starting from the background cell.
    return findCellLabel(coord);
}

void BackgroundMesh::output(const string& name) const
{
    cout << "Writing mesh..." << endl;

    ostringstream fileName;
	fileName << "mesh_" << name;
	string filePath = "./output/" + fileName.str() + ".dat";

	ofstream fp(filePath, ofstream::out);
	if (!fp.is_open())
	{
        cerr << filePath << " write mesh file failed!" << endl;
		exit(-1);
	}

	fp << "Title = \"Background Mesh" << "\"" << endl;
	fp << "Variables = \"X\",\"Y\"" << endl;
	fp << "ZONE N = " << vertIndexNum_ << ", E = " << cellIndexNum_
		<< ", F = FEPOINT, ET = QUADRILATERAL" << endl;

	for (int i = 0; i < vertIndexNum_; i++)
    {
        const int vertI = vertIndexToLabel[i];
        const Vert& vert = verts_[vertI];
        fp << vert.coord().x << " " << vert.coord().y << endl;
    }

	for (int i = 0; i < cellIndexNum_; i++)
    {
        const int cellI = cellIndexToLabel[i];
        const Cell& cell = cells_[cellI];
        for (int v = 0; v < 4; v++)
        {
            const Vert& vert = verts_[cell.cellToVert[v]];
            fp << vert.index() + 1 << " ";
        }

        fp << endl;
    }

	fp.close();
    cout << "Mesh output done." << endl << endl;
}

void BackgroundMesh::writeCache(const string& filePath) const
{
    cout << "Writing mesh cache..." << endl;

    ostringstream fileName;
	fileName << filePath << "mesh_cache.bin";

    ofstream cache(fileName.str(), ios::binary | ios::out);
	if (cache.fail())
	{
		cout << "Writing cache \"" << fileName.str() << "\" failed!" << endl;
		exit(-1);
	}

    const char magic[4] = {'B', 'M', 'C', '2'};
    const std::uint32_t version = 2;
    MeshCacheIO::writeRawArray(cache, magic, 4);
    MeshCacheIO::writeRaw(cache, version);

    const int levelSetCount = levelSetCount_;
    MeshCacheIO::writeRaw(cache, levelSetCount);
    MeshCacheIO::writeRaw(cache, minCoord);
    MeshCacheIO::writeRaw(cache, maxCoord);
    MeshCacheIO::writeRaw(cache, Nx);
    MeshCacheIO::writeRaw(cache, Ny);

    MeshCacheIO::writeRaw(cache, vertIndexNum_);
    MeshCacheIO::writeRaw(cache, hEdgeIndexNum_);
    MeshCacheIO::writeRaw(cache, vEdgeIndexNum_);
    MeshCacheIO::writeRaw(cache, cellIndexNum_);

    MeshCacheIO::writeRawArray(cache, vertIndexToLabel, static_cast<std::size_t>(vertIndexNum_));
    MeshCacheIO::writeRawArray(cache, hEdgeIndexToLabel, static_cast<std::size_t>(hEdgeIndexNum_));
    MeshCacheIO::writeRawArray(cache, vEdgeIndexToLabel, static_cast<std::size_t>(vEdgeIndexNum_));
    MeshCacheIO::writeRawArray(cache, cellIndexToLabel, static_cast<std::size_t>(cellIndexNum_));

    for (int i = 0; i < vertIndexNum_; ++i)
        MeshCacheIO::writeVert(cache, verts_[i], levelSetCount_);
    for (int i = 0; i < hEdgeIndexNum_; ++i)
        MeshCacheIO::writeEdge(cache, hEdges_[i], levelSetCount_);
    for (int i = 0; i < vEdgeIndexNum_; ++i)
        MeshCacheIO::writeEdge(cache, vEdges_[i], levelSetCount_);
    for (int i = 0; i < cellIndexNum_; ++i)
        MeshCacheIO::writeCell(cache, cells_[i], levelSetCount_);

    cache.close();
    cout << "Mesh cache written." << endl << endl;
}

void BackgroundMesh::readCache()
{
    cout << "Reading mesh cache..." << endl;

    ostringstream fileName;
    fileName << "./cache/mesh_cache.bin";

    ifstream cache(fileName.str(), ios::binary | ios::in);
    if (cache.fail())
    {
        cout << "Reading cache \"" << fileName.str() << "\" failed!" << endl;
        exit(-1);
    }

    char magic[4] = {0, 0, 0, 0};
    MeshCacheIO::readRawArray(cache, magic, 4);
    const char expected[4] = {'B', 'M', 'C', '2'};
    if (!MeshCacheIO::magicEquals(magic, expected))
    {
        cout << "Reading cache \"" << fileName.str() << "\" failed: incompatible mesh cache (delete cache)!" << endl;
        exit(-1);
    }
    std::uint32_t version = 0;
    MeshCacheIO::readRaw(cache, version);
    if (version != 2u)
    {
        cout << "Reading cache \"" << fileName.str() << "\" failed: unsupported cache version " << version << "!" << endl;
        exit(-1);
    }

    int levelSetCount = 1;
    MeshCacheIO::readRaw(cache, levelSetCount);
    if (levelSetCount == 0) levelSetCount_ = 0;
    else setLevelSetCount(levelSetCount);

    Vector2D minCoordCache;
    Vector2D maxCoordCache;
    int nxCache = 0;
    int nyCache = 0;
    MeshCacheIO::readRaw(cache, minCoordCache);
    MeshCacheIO::readRaw(cache, maxCoordCache);
    MeshCacheIO::readRaw(cache, nxCache);
    MeshCacheIO::readRaw(cache, nyCache);

    if (nxCache != Nx || nyCache != Ny ||
        minCoordCache.x != minCoord.x || minCoordCache.y != minCoord.y ||
        maxCoordCache.x != maxCoord.x || maxCoordCache.y != maxCoord.y)
    {
        cout << "Reading cache \"" << fileName.str() << "\" failed: mesh dimensions mismatch (delete cache)." << endl;
        exit(-1);
    }

    MeshCacheIO::readRaw(cache, vertIndexNum_);
    MeshCacheIO::readRaw(cache, hEdgeIndexNum_);
    MeshCacheIO::readRaw(cache, vEdgeIndexNum_);
    MeshCacheIO::readRaw(cache, cellIndexNum_);
    vertLabelCapacity_ = vertIndexNum_;
    hEdgeLabelCapacity_ = hEdgeIndexNum_;
    vEdgeLabelCapacity_ = vEdgeIndexNum_;
    cellLabelCapacity_ = cellIndexNum_;

    MeshCacheIO::readRawArray(cache, vertIndexToLabel, static_cast<std::size_t>(vertIndexNum_));
    MeshCacheIO::readRawArray(cache, hEdgeIndexToLabel, static_cast<std::size_t>(hEdgeIndexNum_));
    MeshCacheIO::readRawArray(cache, vEdgeIndexToLabel, static_cast<std::size_t>(vEdgeIndexNum_));
    MeshCacheIO::readRawArray(cache, cellIndexToLabel, static_cast<std::size_t>(cellIndexNum_));

    for (int i = 0; i < vertIndexNum_; ++i)
        MeshCacheIO::readVert(cache, verts_[i], levelSetCount_);
    for (int i = 0; i < hEdgeIndexNum_; ++i)
        MeshCacheIO::readEdge(cache, hEdges_[i], levelSetCount_);
    for (int i = 0; i < vEdgeIndexNum_; ++i)
        MeshCacheIO::readEdge(cache, vEdges_[i], levelSetCount_);
    for (int i = 0; i < cellIndexNum_; ++i)
        MeshCacheIO::readCell(cache, cells_[i], levelSetCount_);

    cache.close();
    cout << "Mesh cache read." << endl << endl;
}
