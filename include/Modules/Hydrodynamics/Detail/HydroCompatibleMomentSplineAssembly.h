#pragma once
#include "Modules/Hydrodynamics/detail/HydroCompatibleMomentSplineNative.h"
#include "Modules/Hydrodynamics/HydroCompatibleMomentComplement.h"
#include <functional>
#include <optional>

namespace Quantum::Hydrodynamics {
// Layout-owned, phase-independent assembly of T^T B T. Only the pattern and
// coordinate maps are retained; every material/transport B is supplied anew.
class MomentSplineCoordinateAssembly {
public:
    using Matrix=HydroCompatibleMomentSpace::Matrix;
    struct Cell {
        Matrix map;
        std::vector<int> positions;
    };
    std::vector<int> offsets;
    Matrix transform;
    Eigen::SparseMatrix<double,Eigen::RowMajor> rowTransform;
    // A byte-valued structural mask avoids keeping a second full array of
    // double coefficients merely to initialize each fresh global matrix.
    Eigen::SparseMatrix<char> pattern;
    std::vector<Cell> cells;
    MomentSplineCoordinateAssembly(Matrix,std::vector<int>);
    Matrix assemble(const std::function<const Eigen::MatrixXd&(int)>&,
        const std::vector<double>* scales=nullptr) const;
protected:
    MomentSplineCoordinateAssembly()=default;
    void initialize(Matrix,std::vector<int>);
};
struct HydroCompatibleMomentSpline::NativeLayout::CoordinateAssembly:MomentSplineCoordinateAssembly {
    std::optional<HydroCompatibleMomentComplement> complement;
    // Lazily prepared under the layout mutex, one independently owned cell
    // per worker. This is the constant-mode observation metric, not inertia.
    std::vector<Eigen::MatrixXd> geometricGram;
    CoordinateAssembly(const NativeLayout&,const HydroCompatibleMomentComplement*);
};
}
