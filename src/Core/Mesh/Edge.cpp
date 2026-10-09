#include "Core/Mesh/Edge.h"

using namespace std;

Edge::Edge() {}

Edge::Edge
(
    int label,
    Vector2D midpointCoord, 
    double length,
    bool isHorizontal
) :
    label_(label),
    midpointCoord_(midpointCoord), 
    length_(length),
    isHorizontal_(isHorizontal)
{}

Edge::~Edge() {}