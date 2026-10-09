#include "Core/Geometry/Line.h"

using namespace std;

Line::Line() {}

void Line::reverse()
{
    Vector2D temp = points[0];
    points[0] = points[1];
    points[1] = temp;
}

void Line::setFeaturePoint(int i, Vector2D _featurePoint)
{
    if (i > parameters_.featurePoints.size())
    {
        cerr << "[Line] ERROR: feature point index out of range." << endl;
        exit(1);
    }
    else
    {
        parameters_.featurePoints[i] = _featurePoint;
    }
}