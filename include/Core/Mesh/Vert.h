#pragma once

#include <vector>
#include "Core/Base/Base.h"
#include "Core/Base/Vector2D.h"
#include "Core/Base/Enumeration.h"
#include "Core/Base/Structures.h"

class Vert
{
public:
    std::vector<unsigned char> onCutEdge;
    std::vector<unsigned char> onCutCell;
    std::vector<Vector2D> nearestNodeCoord;
    std::vector<int> nearestNodeLabel;

    int vertToHEdge[2] = {-1, -1};
    int vertToVEdge[2] = {-1, -1};
    int vertToCell[4] = {-1, -1, -1, -1};

    int neiVertices[4] = {-1, -1, -1, -1};

    int supVert[2] = {-1, -1};
    int supCell = -1;

    int son[1] = {-1};
    int parent = -1;

    int inPhaseNum = 0;
    int relatedVertsNum = 0;
    int relatedVerts[6] = {-1, -1, -1, -1, -1, -1};

    DiscCoeff lapCoeff;
    DiscCoeff selfCoeff;
    DiscCoeff nCoeff;
    DiscCoeff tCoeff;
    DiscCoeff xCoeff;
    DiscCoeff yCoeff;
    DiscCoeff nnCoeff;
    DiscCoeff ttCoeff;

    Vert();
    Vert(int label, Vector2D coord);
    ~Vert();

    void ensureLevelSetSize(std::size_t n)
    {
        if (onCutEdge.size() != n)
        {
            onCutEdge.assign(n, 0u);
            onCutCell.assign(n, 0u);
            nearestNodeCoord.assign(n, Vector2D(0.0, 0.0));
            nearestNodeLabel.assign(n, -1);
        }
    }

    int label() const { return label_; }
    int index() const { return index_; }
    Vector2D coord() const { return coord_; }
    bool needUpdate() const { return needUpdate_; }
    TVertType type() const { return type_; }
    Boundary boundary() const { return boundary_; }
    int phase() const { return phase_; }

    void setLabel(int label) { label_ = label; }
    void setIndex(int index) { index_ = index; }
    void setCoord(Vector2D coord) { coord_ = coord; }
    void setNeedUpdate(bool needUpdate) { needUpdate_ = needUpdate; }
    void setType(TVertType type) { type_ = type; }
    void setBoundary(Boundary boundary) { boundary_ = boundary; }
    void setPhase(int phase) { phase_ = phase; }

private:
    int label_ = -1;
    int index_ = -1;
    Vector2D coord_;

    bool needUpdate_ = false;

    Boundary boundary_ = Boundary::non;
    TVertType type_ = TVertType::non;

    int phase_ = -1;
};
