#pragma once

#include "Core/Base/Vector2D.h"
#include "Core/Geometry/BSpline.h"

#include <vector>

// Helper function to check if a point is on a line segment
bool isPointOnSegment(const Vector2D &p1, const Vector2D &p2, const Vector2D &point);
// Helper function to check if two line segments intersect
bool doLineSegmentsIntersect(const Vector2D &p1, const Vector2D &p2, const Vector2D &q1, const Vector2D &q2, double epsilon);
// Function to check if a line segment intersects the B-spline curve
bool isIntersecting(const Vector2D &p1, const Vector2D &p2, BSplineCurve &bspline);
// Function to check if a line segment intersects a pre-sampled polyline
bool isIntersecting(const Vector2D &p1, const Vector2D &p2, const std::vector<Vector2D> &polyline);

double minmod(const double& a, const double& b);

double operatorHG(const double& a, const double& b, const double& c, const double& d, const double& s);

double sgn(const double& a);

double mod(const double& a, const double& b);

double mod2(const double& a, const double& b);