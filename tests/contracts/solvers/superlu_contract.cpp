#include <Eigen/Sparse>
#include <Eigen/SuperLUSupport>
#include <cmath>
#include <iostream>

int main()
{
    Eigen::SparseMatrix<double> a(3, 3);
    const std::vector<Eigen::Triplet<double>> entries{
        {0,0,4}, {0,1,1}, {1,0,1}, {1,1,3}, {1,2,1}, {2,1,1}, {2,2,2}};
    a.setFromTriplets(entries.begin(), entries.end());
    Eigen::Vector3d expected(1, -2, 3);
    Eigen::VectorXd rhs = a * expected;
    Eigen::SuperLU<Eigen::SparseMatrix<double>> solver;
    solver.compute(a);
    if (solver.info() != Eigen::Success) return 1;
    Eigen::VectorXd x = solver.solve(rhs);
    if (solver.info() != Eigen::Success || !x.allFinite() ||
        (x - expected).norm() > 1e-12 || (a*x - rhs).norm() > 1e-12) return 2;
    std::cout << "SuperLU factorization and solve passed\n";
}
