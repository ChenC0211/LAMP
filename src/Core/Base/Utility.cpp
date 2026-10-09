#include "Core/Base/Utility.h"
#include "Core/Base/Base.h"
#include <cmath> // For std::fabs


// Function to check if a line segment intersects the B-spline curve
bool isIntersecting(const Vector2D &p1, const Vector2D &p2, BSplineCurve &bspline) {
    // Define a small epsilon value for floating-point comparisons
    const double epsilon = 1e-6;

    // Iterate through the segments of the B-spline curve
    for (size_t i = 0; i < bspline.getNumSegments(); ++i) {
        Vector2D splineStart = bspline.evaluate(i / static_cast<double>(bspline.getNumSegments()));
        Vector2D splineEnd = bspline.evaluate((i + 1) / static_cast<double>(bspline.getNumSegments()));

        // Check if the line segment intersects with the current B-spline segment
        if (doLineSegmentsIntersect(p1, p2, splineStart, splineEnd, epsilon)) {
            return true;
        }
    }
    return false;
}

// Function to check if a line segment intersects a pre-sampled polyline
bool isIntersecting(const Vector2D &p1, const Vector2D &p2, const std::vector<Vector2D> &polyline) {
    const double epsilon = 1e-6;
    if (polyline.size() < 2)
        return false;

    for (size_t i = 1; i < polyline.size(); ++i) {
        const Vector2D& splineStart = polyline[i - 1];
        const Vector2D& splineEnd = polyline[i];
        if (doLineSegmentsIntersect(p1, p2, splineStart, splineEnd, epsilon)) {
            return true;
        }
    }
    return false;
}

// Helper function to check if two line segments intersect
bool doLineSegmentsIntersect(const Vector2D &p1, const Vector2D &p2, const Vector2D &q1, const Vector2D &q2, double epsilon) {
    auto orientation = [epsilon](const Vector2D &a, const Vector2D &b, const Vector2D &c) {
        double val = (b.y - a.y) * (c.x - b.x) - (b.x - a.x) * (c.y - b.y);
        if (fabs(val) < epsilon) return 0; // collinear
        return (val > 0) ? 1 : 2; // clock or counterclock wise
    };

    int o1 = orientation(p1, p2, q1);
    int o2 = orientation(p1, p2, q2);
    int o3 = orientation(q1, q2, p1);
    int o4 = orientation(q1, q2, p2);

    // General case
    if (o1 != o2 && o3 != o4) return true;

    // Special cases
    if (o1 == 0 && isPointOnSegment(p1, p2, q1)) return true;
    if (o2 == 0 && isPointOnSegment(p1, p2, q2)) return true;
    if (o3 == 0 && isPointOnSegment(q1, q2, p1)) return true;
    if (o4 == 0 && isPointOnSegment(q1, q2, p2)) return true;

    return false;
}

// Helper function to check if a point is on a line segment
bool isPointOnSegment(const Vector2D &p1, const Vector2D &p2, const Vector2D &point) {
    double crossProduct = (point.y - p1.y) * (p2.x - p1.x) - (point.x - p1.x) * (p2.y - p1.y);
    if (fabs(crossProduct) > 1e-16) {
        return false;
    }

    double dotProduct = (point.x - p1.x) * (p2.x - p1.x) + (point.y - p1.y) * (p2.y - p1.y);
    if (dotProduct < 0) {
        return false;
    }

    double squaredLength = (p2.x - p1.x) * (p2.x - p1.x) + (p2.y - p1.y) * (p2.y - p1.y);
    if (dotProduct > squaredLength) {
        return false;
    }

    return true;
}

double minmod(const double& a, const double& b)
{
	if (a * b <= MINVAL)
		return 0.0;

	double sgn = (a > 0.0) ? 1.0 : -1.0;
	double min_abs = (fabs(a) < fabs(b)) ? fabs(a) : fabs(b);

	return sgn * min_abs;
}

double operatorHG(const double& a, const double& b, const double& c, const double& d, const double& s)
{
    if (s <= 0.0)
	{
        return sqrt(std::max(pow(std::max(a, 0.0), 2.0), pow(std::min(b, 0.0), 2.0))
			+ std::max(pow(std::max(c, 0.0), 2.0), pow(std::min(d, 0.0), 2.0)));
    }
    else
	{
        return sqrt(std::max(pow(std::min(a, 0.0), 2.0), pow(std::max(b, 0.0), 2.0))
			+ std::max(pow(std::min(c, 0.0), 2.0), pow(std::max(d, 0.0), 2.0)));
    }
}

double sgn(const double& a)
{
    if (a > 0.0)
		return 1.0;
    else
        return -1.0;
}

double mod(const double& a, const double& b)
{
    return sqrt(a * a + b * b);
}

double mod2(const double& a, const double& b)
{
    return a * a + b * b;
}