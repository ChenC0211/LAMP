#include "Core/Base/Interpolation.h"
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

namespace Interpolation
{

static double linearCoreFunction(const double& r)
{
    double fabs_r = fabs(r);

    if (fabs_r < 1.0)
    {
        return 1.0 - fabs_r;
    }
    else
    {
        return 0.0;
    }
}
// static double PesKinCoreFunction(const double& r)
// {
//     double fabs_r = fabs(r);

//     if (fabs_r < 2.0)
//     {
//         return 0.25 * (1.0 + cos(0.5 * Pi * fabs_r));
//     }
//     else
//     {
//         return 0.0;
//     }
// }
static double PesKinCoreFunction(const double& r)
{
    double fabs_r = std::fabs(r);

    if (fabs_r < 1.0)
    {
        return 0.125 * (3 - 2 * fabs_r + std::sqrt(1 + 4 * fabs_r - 4 * fabs_r * fabs_r));
    }
    else if (fabs_r < 2.0)
    {
        return 0.125 * (5 - 2 * fabs_r - std::sqrt(-7 + 12 * fabs_r - 4 * fabs_r * fabs_r));
    }
    else
    {
        return 0.0;
    }
}

}; // namespace Interpolation


namespace Interpolation
{

double leastSquareInterpolation
(
    const std::vector<Vector2D>& points,
    const std::vector<double>& values,
    const Vector2D& aimCoord,
    bool secondOrder = false
)
{
    int n = static_cast<int>(points.size());
    if (n != values.size())
    {
        cout << "[Interpolation] ERROR: points.size() != values.size()." << endl;
        exit(0);
    }

    if (secondOrder && n < 6)
    {
        cout << "[Interpolation] ERROR: second-order interpolation requires at least 6 points." << endl;
        exit(0);
    }
    if (!secondOrder && n < 3)
    {
        cout << "[Interpolation] ERROR: first-order interpolation requires at least 3 points." << endl;
        exit(0);
    }

    MatrixXd A(n, secondOrder ? 6 : 3);
    VectorXd b(n);

    for (int i = 0; i < n; i++)
    {
        A(i, 0) = 1.0;
        A(i, 1) = points[i].x;
        A(i, 2) = points[i].y;
        if (secondOrder)
        {
            A(i, 3) = points[i].x * points[i].x;
            A(i, 4) = points[i].y * points[i].y;
            A(i, 5) = points[i].x * points[i].y;
        }
        b(i) = values[i];
    }

    // 求解最小二乘问题 (A^T * A) * p = A^T * b
    VectorXd p = (A.transpose() * A).ldlt().solve(A.transpose() * b);

    // 插值：计算目标点的值
    double result = 0.0;
    if (secondOrder)
    {
        result = p(0) + p(1) * aimCoord.x + p(2) * aimCoord.y + p(3) * aimCoord.x * aimCoord.x + p(4) * aimCoord.y * aimCoord.y + p(5) * aimCoord.x * aimCoord.y;
    }
    else
    {
        result = p(0) + p(1) * aimCoord.x + p(2) * aimCoord.y;
    }

    return result;
}

vector<double> leastSquareInterpolationCoeff
(
    const std::vector<Vector2D>& points,
    const Vector2D& aimCoord,
    bool secondOrder
)
{
    int n = static_cast<int>(points.size());

    if (secondOrder && n < 6)
    {
        cout << "[Interpolation] ERROR: second-order interpolation requires at least 6 points." << endl;
        exit(0);
    }
    if (!secondOrder && n < 3)
    {
        cout << "[Interpolation] ERROR: first-order interpolation requires at least 3 points." << endl;
        exit(0);
    }

    MatrixXd A(n, secondOrder ? 6 : 3);

    for (int i = 0; i < n; i++)
    {
        A(i, 0) = 1.0;
        A(i, 1) = points[i].x;
        A(i, 2) = points[i].y;
        if (secondOrder)
        {
            A(i, 3) = points[i].x * points[i].x;
            A(i, 4) = points[i].y * points[i].y;
            A(i, 5) = points[i].x * points[i].y;
        }
    }

    // 构造目标点的基函数向量 B
    VectorXd B(secondOrder ? 6 : 3);
    B(0) = 1.0;
    B(1) = aimCoord.x;
    B(2) = aimCoord.y;
    if (secondOrder)
    {
        B(3) = aimCoord.x * aimCoord.x;
        B(4) = aimCoord.y * aimCoord.y;
        B(5) = aimCoord.x * aimCoord.y;
    }

    // 计算权重 w = B * (A^T * A)^{-1} * A^T
    MatrixXd AtA = A.transpose() * A;
    MatrixXd inv_AtA = AtA.inverse();
    MatrixXd w = B.transpose() * inv_AtA * A.transpose();

    // 返回权重列表
    std::vector<double> weights(w.data(), w.data() + w.size());

    return weights;
}

double diracInterpolation
(
    const Vector2D& LagrangianCoord,
    const Vector2D& EulerianCoord,
    const Vector2D& h
)
{
    const Vector2D distance = EulerianCoord - LagrangianCoord;
    const double rx = distance.x / h.x;
    const double ry = distance.y / h.y;
    return PesKinCoreFunction(rx) * PesKinCoreFunction(ry);
}

}; // namespace Interpolation