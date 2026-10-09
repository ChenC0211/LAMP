//                       网格编号信息
// ---------------------------------------------------------------------------
//      *-----*-----*             3
//      |  3  |  2  |           *---*
//      *-----*-----*          0|   |2
//      |  0  |  1  |           *---*
//      *-----*-----*             1
//    父网格与子网格关系    单元与相邻单元的关系
// ---------------------------------------------------------------------------
//      3-----------2        *-----1-----*
//      |           |        |           |
//      |           |        0           1
//      |           |        |           |
//      0-----------1        *-----0-----*
//      单元到点的关系         单元到边的关系
// ---------------------------------------------------------------------------
//            3                    *
//            |                    1
//      0-----*-----2        *--0--*--1--*
//            |                    0
//            1                    *
//     点与相邻点的关系         点到边的关系
// ---------------------------------------------------------------------------

#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <fstream>
#include "Core/Base/Vector2D.h"
#include "Core/Base/Base.h"
#include "Core/Mesh/Vert.h"
#include "Core/Mesh/Edge.h"
#include "Core/Mesh/Cell.h"
#include "Core/IO/DictionaryMesh.h"

class BackgroundMesh
{
public:
    const Vector2D minCoord;
    const Vector2D maxCoord;
    const int Nx;
    const int Ny;
    const Vector2D cellSize;

    int* vertIndexToLabel = nullptr;
    int* hEdgeIndexToLabel = nullptr;
    int* vEdgeIndexToLabel = nullptr;
    int* cellIndexToLabel = nullptr;

    std::vector<std::vector<int>> negFaceVerts;
    std::vector<std::vector<int>> posFaceVerts;

    int vertIndexNum() const { return vertIndexNum_; }
    int hEdgeIndexNum() const { return hEdgeIndexNum_; }
    int vEdgeIndexNum() const { return vEdgeIndexNum_; }
    int cellIndexNum() const { return cellIndexNum_; }

    int vertLabelCapacity() const { return vertLabelCapacity_; }
    int hEdgeLabelCapacity() const { return hEdgeLabelCapacity_; }
    int vEdgeLabelCapacity() const { return vEdgeLabelCapacity_; }

    // Allocated capacity of label-addressable arrays (cells_/verts_/edges).
    // For uniform meshes this equals cellIndexNum(); for quadtree meshes it can
    // be larger after internal buffer growth.
    int cellLabelCapacity() const { return cellLabelCapacity_; }

    Vert& verts(int label) const { return verts_[label]; }
    Edge& hEdges(int label) const { return hEdges_[label]; }
    Edge& vEdges(int label) const { return vEdges_[label]; }
    Cell& cells(int label) const { return cells_[label]; }

    int levelSetCount() const { return levelSetCount_; }

    BackgroundMesh(Vector2D _minCoord, Vector2D _maxCoord, int _Nx, int _Ny);
    BackgroundMesh(const DictionaryMesh& dictMesh);
    ~BackgroundMesh();

    int findCellLabel(const Vector2D& coord) const;
    int findBackgroundCellLabel(const Vector2D& coord) const;

    // Fast in-cell check and hint-based cell lookup (stay-in-cell + neighbor ring + fallback).
    // eps is a small tolerance for boundary inclusion.
    bool isInsideCell(int cellLabel, const Vector2D& coord, double eps = MINVAL) const;
    int findCellLabelFromHint(int hintCellLabel, const Vector2D& coord, double eps = MINVAL) const;

    void output(const std::string& name) const;
    void writeCache(const std::string& filePath) const;
    void readCache();

protected:
    // For an owning derived clone only: copy lattice metadata, leave storage
    // empty. The clone fills native records directly, without buildMesh().
    struct CopyMetadataTag {};
    BackgroundMesh(const BackgroundMesh& source, CopyMetadataTag);

    void setLevelSetCount(int n) { levelSetCount_ = (n > 0) ? n : 1; }
    void ensureLevelSetStorage();

protected:
    Vert* verts_ = nullptr;
    Edge* hEdges_ = nullptr;
    Edge* vEdges_ = nullptr;
    Cell* cells_ = nullptr;

    int vertIndexNum_ = 0;
    int hEdgeIndexNum_ = 0;
    int vEdgeIndexNum_ = 0;
    int cellIndexNum_ = 0;

    int vertLabelCapacity_ = 0;
    int hEdgeLabelCapacity_ = 0;
    int vEdgeLabelCapacity_ = 0;

    // Capacity (allocated length) of cells_ in label space.
    int cellLabelCapacity_ = 0;

    int levelSetCount_ = 1;
    
    void buildMesh();

    void createVertices();
    void createEdges();
    void createCells();

    void findNeighbourCells();

    void defineTypes();

    template<typename T>
    void statistics(T*& t, int num);

};
