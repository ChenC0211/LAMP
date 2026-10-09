#include "Modules/Hydrodynamics/HydroFittedProjection.h"
#include "Modules/Hydrodynamics/HydroFittedDomain.h"

#include <cmath>
#include <stdexcept>

namespace HydroFittedProjection
{
namespace
{
void require(bool condition, const char* message)
{
    if (!condition) throw std::invalid_argument(message);
}

void finiteMatrix(const Eigen::SparseMatrix<double>& matrix)
{
    for (int column = 0; column < matrix.outerSize(); ++column)
        for (Eigen::SparseMatrix<double>::InnerIterator entry(matrix, column); entry; ++entry)
            require(std::isfinite(entry.value()), "non-finite fitted operator coefficient");
}
}

System assemble(const HydroFittedDomain& domain, double dt,
    const std::vector<double>& densityByVelocityRow, const Eigen::VectorXd& predictor)
{
    const auto& divergence = domain.divergence();
    const int velocityCount = static_cast<int>(domain.nativeByVelocityRow().size());
    require(velocityCount > 0 && divergence.cols() == velocityCount && divergence.rows() > 0,
        "invalid fitted projection domain");
    require(std::isfinite(dt) && dt > 0.0, "fitted projection requires finite positive dt");
    require(densityByVelocityRow.size() == static_cast<std::size_t>(velocityCount)
        && predictor.size() == velocityCount && predictor.allFinite(), "invalid fitted predictor/density size or values");
    System system(domain);
    auto& equations = system.equations_;
    system.inverseMass_.resize(velocityCount);
    Eigen::SparseMatrix<double> weightedTranspose = divergence.transpose();
    if (!domain.uniformLeaves()) weightedTranspose = -domain.pressureGradientLoad();
    for (int i = 0; i < velocityCount; ++i)
    {
        const double density = densityByVelocityRow[i];
        const double volume = domain.faces()[domain.nativeByVelocityRow()[i]].dualVolume;
        require(std::isfinite(density) && density > 0.0 && std::isfinite(volume) && volume > 0.0,
            "fitted inverse mass requires physical positive fluid density and dual volume");
        const double mass = density * volume;
        require(std::isfinite(mass) && mass > 0.0, "unrepresentable fitted face mass");
        system.inverseMass_[i] = dt / mass;
        require(std::isfinite(system.inverseMass_[i]) && system.inverseMass_[i] > 0.0,
            "unrepresentable fitted inverse mass");
    }
    for (int column = 0; column < weightedTranspose.outerSize(); ++column)
        for (Eigen::SparseMatrix<double>::InnerIterator entry(weightedTranspose, column); entry; ++entry)
            entry.valueRef() *= system.inverseMass_[entry.row()];
    equations.matrix = divergence * weightedTranspose;
    equations.matrix.makeCompressed();
    equations.rhs = -(divergence * predictor);
    finiteMatrix(equations.matrix);
    require(equations.rhs.allFinite(), "non-finite fitted projection rhs");
    for (int i = 0; i < equations.matrix.rows(); ++i)
        require(equations.matrix.coeff(i, i) > 0.0, "fitted projection has unanchored or degenerate continuity row");
    return system;
}

void System::updateRightHandSide(const Eigen::VectorXd& predictor)
{
    const auto& divergence = domain_->divergence();
    require(predictor.size() == divergence.cols() && predictor.allFinite(),
        "invalid fitted projection predictor size or values");
    Eigen::VectorXd rhs = -(divergence * predictor);
    require(rhs.allFinite(), "non-finite fitted projection rhs");
    equations_.rhs = std::move(rhs);
}

Eigen::VectorXd System::correctVelocity(const Eigen::VectorXd& predictor,
    const Eigen::VectorXd& pressureIncrement) const
{
    const auto& gradient = domain_->pressureGradientLoad();
    require(predictor.size() == gradient.rows() && inverseMass_.size() == gradient.rows()
        && pressureIncrement.size() == gradient.cols() && predictor.allFinite() && pressureIncrement.allFinite()
        && inverseMass_.allFinite(), "invalid fitted velocity correction operands");
    for (int i = 0; i < inverseMass_.size(); ++i)
        require(inverseMass_[i] > 0.0, "fitted velocity correction requires positive inverse mass");
    Eigen::VectorXd corrected = predictor
        - (inverseMass_.array() * (gradient * pressureIncrement).array()).matrix();
    require(corrected.allFinite(), "non-finite fitted corrected velocity");
    return corrected;
}
}
