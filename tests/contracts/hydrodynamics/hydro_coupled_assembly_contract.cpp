#include "Modules/Hydrodynamics/HydroCoupledBalance.h"
#include <Eigen/SuperLUSupport>
#include <cstring>
#include <iostream>
#include <type_traits>

namespace {
using Sparse = Eigen::SparseMatrix<double>;
using Vector = Eigen::VectorXd;
using HydroCoupledBalance::Assembly;

void check(bool condition, const char* message)
{
    if (!condition) throw std::runtime_error(message);
}
template<class F> void rejects(F&& action)
{
    try { action(); }
    catch (const std::invalid_argument&) { return; }
    throw std::runtime_error("invalid coupled assembly/load was accepted");
}
Sparse sparse(const Eigen::MatrixXd& values) { return values.sparseView(); }

void checkMatrix(const Sparse& actual, const Eigen::MatrixXd& expected)
{
    check(actual.isCompressed() && actual.rows() == expected.rows()
        && actual.cols() == expected.cols(), "coupled matrix shape/storage");
    check((Eigen::MatrixXd(actual).array() == expected.array()).all(), "coupled block coefficients/signs");
}

void closedReference()
{
    Eigen::Matrix2d a; a << 4, 1, 1, 3;
    const Sparse momentum = sparse(a);
    Assembly rows(momentum, 2, 8, 0);
    rows.addPressureGradient(0, 0, 1); rows.addPressureGradient(0, 1, -1);
    rows.addPressureGradient(1, 0, -1); rows.addPressureGradient(1, 1, 1);
    rows.addDivergence(0, 0, 1); rows.addDivergence(0, 1, -1);
    rows.addDivergence(1, 0, -1); rows.addDivergence(1, 1, 1);
    const auto equations = rows.finish();
    Eigen::Matrix4d expected;
    expected << 4,1,1,-1, 1,3,-1,1, 0,0,1,0, 1,-1,0,0;
    checkMatrix(equations.matrix(), expected);
    rejects([&] { rows.finish(); });
    rejects([&] { rows.addPressureGradient(0, 0, 1); });
    rejects([&] { rows.addDivergence(1, 0, 1); });
    Vector b(2), velocity(2), pressure(2), load(2), flux(2);
    b << 7,11; velocity << 1,-2; pressure << 0,0.5;
    load << -0.5,0.5; flux << 123,-3; // Reference row is replaced, not added.
    const auto rhs = equations.rightHandSide(b, velocity, load, flux);
    Vector expectedRhs(4); expectedRhs << 5.5,15.5,0,-3;
    check((rhs.array() == expectedRhs.array()).all(), "native load/gauge RHS");
    Eigen::SuperLU<Sparse> solver;
    solver.compute(equations.matrix());
    check(solver.info() == Eigen::Success, "closed coupled factorization");
    const Vector correction = solver.solve(rhs);
    check(solver.info() == Eigen::Success, "closed coupled solve");
    Vector expectedCorrection(4); expectedCorrection << 1,4,0,2.5;
    check((correction - expectedCorrection).norm() < 1e-12, "closed pressure/velocity correction");
    rejects([&] { equations.rightHandSide(b, velocity, pressure); }); // No sparse C/D retained.
}

void anchoredNonAdjoint()
{
    Eigen::Matrix3d a; a << 4,1,0, 1,3,0, 0,0,2;
    Eigen::Matrix<double,3,2> c; c << 2,0, -1,3, 1,-2;
    Eigen::Matrix<double,2,3> d; d << 1,-1,0, 0,2,-1;
    const Sparse momentum = sparse(a), gradient = sparse(c), divergence = sparse(d);
    const auto equations = HydroCoupledBalance::assemble(momentum, gradient, divergence);
    Eigen::Matrix<double,5,5> expected;
    expected << 4,1,0,2,0, 1,3,0,-1,3, 0,0,2,1,-2, -1,1,0,0,0, 0,-2,1,0,0;
    checkMatrix(equations.matrix(), expected); // Keep both continuity rows; C != D^T.
    const auto* matrixValues = equations.matrix().valuePtr();
    const std::vector<double> coefficients(matrixValues, matrixValues + equations.matrix().nonZeros());
    Eigen::SuperLU<Sparse> solver;
    solver.compute(equations.matrix());
    check(solver.info() == Eigen::Success, "non-adjoint coupled factorization");
    Vector b(3), velocity(3), pressure(2), expectedVelocity(3), expectedPressure(2);
    b << 16,2,13; velocity << -2,1,0; pressure << 0.5,-2;
    expectedVelocity << 2,2,4; expectedPressure << 3,-1;
    // Repeated RHS/corrections use one actual SuperLU factorization. These
    // solve the same momentum + continuity equations, not a time substep.
    for (int correctionStep = 0; correctionStep < 3; ++correctionStep)
    {
        const Vector rhs = equations.rightHandSide(b, velocity, pressure);
        const Vector correction = solver.solve(rhs);
        check(solver.info() == Eigen::Success && correction.allFinite(), "repeated coupled correction");
        velocity += correction.head(3); pressure += correction.tail(2);
        check((velocity - expectedVelocity).norm() < 1e-12
            && (pressure - expectedPressure).norm() < 1e-12, "non-adjoint physical solution");
        check((a * velocity + c * pressure - b).norm() < 1e-12 && (d * velocity).norm() < 1e-12,
            "full momentum/continuity equations");
        check(equations.matrix().valuePtr() == matrixValues && std::memcmp(matrixValues,
            coefficients.data(), coefficients.size() * sizeof(double)) == 0, "refinement rebuilt/changed operator");
    }
    const Vector goodRhs = equations.rightHandSide(b, velocity, pressure);
    auto invalidPressure = pressure;
    invalidPressure[0] = std::numeric_limits<double>::quiet_NaN();
    rejects([&] { equations.rightHandSide(b, velocity, invalidPressure); });
    rejects([&] { equations.rightHandSide(b, Vector::Zero(2), pressure); });
    rejects([&] { equations.rightHandSide(Vector::Zero(2), velocity, pressure); });
    const Vector repeatedRhs = equations.rightHandSide(b, velocity, pressure);
    check(std::memcmp(goodRhs.data(), repeatedRhs.data(), goodRhs.size() * sizeof(double)) == 0,
        "invalid RHS changed frozen operator");
    const Vector finalLoad = gradient * pressure;
    const auto balance = HydroCoupledBalance::momentumResidual(momentum, velocity, finalLoad, b,
        [](int) { return 1.0; });
    check(balance.relative < 1e-12 && balance.maxForceDensity < 1e-12, "shared momentum residual");

    // Native row construction and sparse strategies produce the same block;
    // the operator survives destruction of its transient triplet builder.
    const auto native = [&] {
        Assembly rows(momentum, 2, 9);
        rows.addPressureGradient(0,0,2); rows.addPressureGradient(1,0,-1);
        rows.addPressureGradient(1,1,3); rows.addPressureGradient(2,0,1); rows.addPressureGradient(2,1,-2);
        rows.addDivergence(0,0,1); rows.addDivergence(0,1,-1);
        rows.addDivergence(1,1,2); rows.addDivergence(1,2,-1);
        return rows.finish();
    }();
    checkMatrix(native.matrix(), expected);
    const Vector flux = divergence * velocity;
    check((native.rightHandSide(b, velocity, finalLoad, flux) - goodRhs).norm() < 1e-12,
        "native/sparse residual strategy disagreement");
}
}

int main()
{
    static_assert(!std::is_copy_constructible_v<HydroCoupledBalance::Operator>);
    static_assert(!std::is_constructible_v<Assembly, Sparse&&, int, std::size_t>);
    static_assert(std::is_same_v<decltype(std::declval<const HydroCoupledBalance::Operator&>().matrix()), const Sparse&>);
    try
    {
        // A nonzero physical defect must retain both rows in its relative
        // norm, even when one prescribed-velocity row has no force measure.
        Sparse mass(2,2); mass.insert(0,0)=2.; mass.insert(1,1)=4.;
        Vector velocity(2), pressureLoad(2), rhs(2);
        velocity<<3.,-1.; pressureLoad<<1.,5.; rhs<<5.,0.;
        const auto unclosed=HydroCoupledBalance::momentumResidual(mass,velocity,pressureLoad,rhs,
            [](int row){return row==0?.5:0.;});
        const double expected=std::sqrt(5.)/(std::sqrt(52.)+std::sqrt(26.));
        check(std::abs(unclosed.relative-expected)<1e-15 && unclosed.maxForceDensity==4.,
            "nonzero momentum defect norm or force-density measure changed");
        closedReference();
        anchoredNonAdjoint();
        std::cout << "Coupled operator: native/sparse blocks, pressure reference, retained operators and repeated SuperLU corrections passed\n";
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
