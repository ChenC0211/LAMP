#ifndef BSPLINECURVE_H
#define BSPLINECURVE_H

#include <vector>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include "Core/Base/Vector2D.h"

class BSplineCurve {
    
public:
    BSplineCurve(const std::vector<Vector2D>& controlPoints, int degree)
        : controlPoints(controlPoints), degree(degree) {
        generateKnotVector();
    }

Vector2D evaluate(double t) const {
    const int n = static_cast<int>(controlPoints.size()) - 1;
    const int k = degree;

    // Find the knot span index
    int span = -1;
    for (int idx = 0; idx < static_cast<int>(knotVector.size()) - 1; ++idx) {
        if (t >= knotVector[idx] && t < knotVector[idx + 1]) {
            span = idx;
            break;
        }
    }
    if (span == -1) {
        if (t == 1.0) {
            span = n;
        } else {
            throw std::out_of_range("Parameter t is out of range.");
        }
    }

    // Compute basis functions
    std::vector<double> N(k + 1, 0.0);
    N[0] = 1.0;
    std::vector<double> left(k + 1), right(k + 1);

    for (int j = 1; j <= k; ++j) {
        left[j] = t - knotVector[span + 1 - j];
        right[j] = knotVector[span + j] - t;
        double saved = 0.0;

        for (int r = 0; r < j; ++r) {
            double temp = N[r] / (right[r + 1] + left[j - r]);
            N[r] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        N[j] = saved;
    }

    // Compute the point on the curve
    Vector2D C(0.0, 0.0);
    for (int j = 0; j <= k; ++j) {
        C = C + controlPoints[span - k + j] * N[j];
    }

    return C;
}

int getNumSegments() const {
    return static_cast<int>(controlPoints.size()) - degree;
}

private:
    std::vector<Vector2D> controlPoints;
    std::vector<double> knotVector;
    int degree;

    void generateKnotVector()
    {
        int n = static_cast<int>(controlPoints.size()) - 1;
        int k = degree;
        int m = n + k + 1; // 节点向量的总长度

        knotVector.resize(m + 1);

        // 前 k+1 个节点设为 0（左端点）
        for (int idx = 0; idx <= k; ++idx)
        {
            knotVector[idx] = 0.0;
        }

        // 中间的节点均匀分布
        for (int idx = k + 1; idx <= n; ++idx)
        {
            knotVector[idx] = (idx - k) / static_cast<double>(n - k + 1);
        }

        // 后 k+1 个节点设为 1（右端点）
        for (int idx = n + 1; idx <= m; ++idx)
        {
            knotVector[idx] = 1.0;
        }
    }
};

#endif // BSPLINECURVE_H