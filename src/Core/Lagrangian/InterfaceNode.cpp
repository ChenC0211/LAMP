#include "Core/Lagrangian/InterfaceNode.h"

using namespace std;

InterfaceNode::InterfaceNode()
{
    for (int i = 0; i < 25; i++)
    {
        cellStencil[i] = -1;
        cellStencilCoeff[i] = 0.0;
    }
}

InterfaceNode::~InterfaceNode() {}

InterfaceNode& InterfaceNode::operator=(const InterfaceNode& node)
{
    for (int i = 0; i < 25; i++)
    {
        cellStencil[i] = node.cellStencil[i];
        cellStencilCoeff[i] = node.cellStencilCoeff[i];
    }

    index_ = node.index();
    label_ = node.label();
    lineIndex_ = node.lineIndex();
    ds_ = node.ds();
    coord_ = node.coord();
    preCoord_ = node.preCoord();
    prePreCoord_ = node.prePreCoord();
    normal_ = node.normal();
    curvature_ = node.curvature();
    curvatureInPlane_ = node.curvatureInPlane();
    relatedNegVertLabel_ = node.relatedNegVertLabel();
    relatedPosVertLabel_ = node.relatedPosVertLabel();
    upNodeLabel_ = node.upNodeLabel();
    downNodeLabel_ = node.downNodeLabel();


    for (int i = 0; i < 3; i++)
    {
        relatedNodes[i] = node.relatedNodes[i];
    }

    return *this;
}