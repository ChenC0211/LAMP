#pragma once

#include <iostream>
#include <vector>
#include "Core/Base/Base.h"
#include "Core/Base/Vector2D.h"
#include "Core/Base/Structures.h"

class InterfaceNode
{
public:
    int cellStencil[25];
    double cellStencilCoeff[25];
    int relatedNodes[3] = {-1, -1, -1};

    IIMCoeff iimCoeff;
    
    InterfaceNode();
    ~InterfaceNode();

    inline int index() const { return index_; }
    inline int label() const { return label_; }
    inline int lineIndex() const { return lineIndex_; }
    inline double ds() const { return ds_; }
    inline Vector2D coord() const { return coord_; }
    inline Vector2D preCoord() const { return preCoord_; }
    inline Vector2D prePreCoord() const { return prePreCoord_; }
    inline Vector2D normal() const { return normal_; }
    inline double curvature() const { return curvature_; }
    inline double curvatureInPlane() const { return curvatureInPlane_; }
    inline int relatedNegVertLabel() const { return relatedNegVertLabel_; }
    inline int relatedPosVertLabel() const { return relatedPosVertLabel_; }
    inline int upNodeLabel() const { return upNodeLabel_; }
    inline int downNodeLabel() const { return downNodeLabel_; }


    void setIndex(const int index) { index_ = index; }
    void setLabel(const int label) { label_ = label; }
    void setLineIndex(const int lineIndex) { lineIndex_ = lineIndex; }
    void setDs(const double ds) { ds_ = ds; }
    void setCoord(const Vector2D& coord) { coord_ = coord; }
    void setPreCoord(const Vector2D& preCoord) { preCoord_ = preCoord; }
    void setPrePreCoord(const Vector2D& prePreCoord) { prePreCoord_ = prePreCoord; }
    void setNormal(const Vector2D& normal) { normal_ = normal; }
    void setCurvature(const double curvature) { curvature_ = curvature; }
    void setCurvatureInPlane(const double curvatureInPlane) { curvatureInPlane_ = curvatureInPlane; }
    void setRelatedNegVertLabel(const int relatedNegVertLabel) { relatedNegVertLabel_ = relatedNegVertLabel; }
    void setRelatedPosVertLabel(const int relatedPosVertLabel) { relatedPosVertLabel_ = relatedPosVertLabel; }
    void setUpNodeLabel(const int upNodeLabel) { upNodeLabel_ = upNodeLabel; }
    void setDownNodeLabel(const int downNodeLabel) { downNodeLabel_ = downNodeLabel; }


    InterfaceNode& operator=(const InterfaceNode& node);

private:
    int index_ = -1;
    int label_ = -1;
    int lineIndex_ = -1;
    double ds_;
    Vector2D coord_;
    Vector2D preCoord_;
    Vector2D prePreCoord_;
    Vector2D normal_;
    double curvature_;
    double curvatureInPlane_;

    int relatedNegVertLabel_ = -1;
    int relatedPosVertLabel_ = -1;
    int upNodeLabel_ = -1;
    int downNodeLabel_ = -1;


};