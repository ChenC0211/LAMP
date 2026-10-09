#pragma once

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <utility>
#include <vector>

// One assembly lifetime, in either native or compact velocity numbering.
// Row kernels and boundary/load contributions write here; the completed RHS
// then travels with the matrix through prediction and coupled correction.
struct HydroMomentumAssembly
{
    std::vector<Eigen::Triplet<double>> entries;
    Eigen::VectorXd rhs;

    // Keep the caller's storage order and original triplet reduction order.
    // Fitted native-interior rows need RowMajor for their subsequent gathering.
    template<int StorageOrder>
    Eigen::VectorXd finish(Eigen::SparseMatrix<double, StorageOrder>& matrix) &&
    {
        matrix.resize(rhs.size(), rhs.size());
        matrix.setFromTriplets(entries.begin(), entries.end());
        matrix.makeCompressed();
        return std::move(rhs);
    }
};
