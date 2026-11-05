#pragma once

#include "Core/Base/Vector2D.h"

#include <Eigen/Core>
#include <array>
#include <cstddef>
#include <functional>
#include <memory>
#include <vector>

namespace LAMP::ElectricDiscretization {

using CellVector = Eigen::Matrix<double, 6, 1>;
using PotentialVector = Eigen::Matrix<double, 10, 1>;
using FacePolynomial = Eigen::Matrix<double, 3, 1>;
using PotentialGradient = Eigen::Matrix<double, 2, 10>;
using PotentialHessian = std::array<PotentialVector, 3>; // xx, xy, yy
using ScalarFunction = std::function<double(Vector2D)>;

struct CellQuadraturePoint {
    Vector2D point;
    double weight = 0.0; // Ordinary dA, before the optional 2*pi*r weight.
};
struct FaceQuadraturePoint {
    Vector2D point;
    Vector2D outwardNormal;
    double weight = 0.0; // Ordinary ds, before the optional 2*pi*r weight.
    double coordinate = 0.0; // Oriented normalized arclength zeta in [-1,1].
};
struct FaceQuadrature {
    double length = 0.0; // Actual h_F; the axis face must be omitted.
    std::vector<FaceQuadraturePoint> nodes;
};
struct LocalGeometry {
    // Center and every volume/face point use this SAME coordinate frame.
    Vector2D center;
    double scale = 0.0; // Positive h_T for the scaled monomial basis.
    bool axisymmetric = false;
    std::vector<CellQuadraturePoint> volume;
    std::vector<FaceQuadrature> faces;
    // World anchor; polynomial arithmetic never adds it back. Only external
    // ScalarFunction callbacks receive world points. RZ requires zero here.
    Vector2D coordinateOrigin{0.0, 0.0};
};

// xi=(x-center.x)/scale, eta=(y-center.y)/scale.
// P2: [1,xi,eta,xi^2,xi*eta,eta^2].
// P3 adds [xi^3,xi^2*eta,xi*eta^2,eta^3] in that order.
// Gradient and Hessian entries are derivatives in physical coordinates.
CellVector p2Basis(Vector2D point, Vector2D center, double scale);
PotentialVector p3Basis(Vector2D point, Vector2D center, double scale);
PotentialGradient p3Gradient(Vector2D point, Vector2D center, double scale);
PotentialHessian p3Hessian(Vector2D point, Vector2D center, double scale);
FacePolynomial faceBasis(double coordinate);
double physicalWeight(Vector2D point, bool axisymmetric);

// Immutable prepared k=2 operator. Local ordering is six cell coefficients,
// followed by three coefficients per face, in the supplied geometry order.
// Owns quadrature and factors, without mesh, physical state, clock or IO.
class LocalOperator {
public:
    LocalOperator(LocalGeometry geometry, double permittivity,
                  double stabilization = 1.0);

    const LocalGeometry& geometry() const;
    std::size_t faceCount() const;
    double permittivity() const;
    double stabilization() const;
    const Eigen::MatrixXd& reconstructionMatrix() const; // 10 x (6+3*N_F)
    const Eigen::MatrixXd& stiffnessMatrix() const; // K
    const Eigen::MatrixXd& condensedMatrix() const; // Cell-only Schur block S
    const Eigen::Matrix<double, 10, 10>& gradientGramMatrix() const;
    const Eigen::MatrixXd& stabilizationMatrix() const;

    CellVector recoverCell(const Eigen::VectorXd& faces,
                           const CellVector& cellLoad = CellVector::Zero()) const;
    PotentialVector reconstruct(const CellVector& cell,
                                const Eigen::VectorXd& faces) const;
    Eigen::VectorXd condensedRhs(const CellVector& cellLoad,
                                const Eigen::VectorXd& faceLoad) const;
    Eigen::VectorXd faceResidual(const CellVector& cell,
                                const Eigen::VectorXd& faces,
                                const Eigen::VectorXd& faceLoad) const;

    // Cell reduction uses ordinary dA; face reduction uses m*ds.
    CellVector reductionCell(const ScalarFunction& value) const;
    Eigen::VectorXd reductionFaces(const ScalarFunction& value) const;
    Eigen::VectorXd reduction(const ScalarFunction& value) const;

    // Physical RHS moments; callers choose the source/flux and its sign.
    CellVector loadCell(const ScalarFunction& source) const;
    FacePolynomial loadFace(std::size_t face, const ScalarFunction& flux) const;

private:
    struct Data;
    const Data& data() const;
    std::shared_ptr<const Data> data_;
};

} // namespace LAMP::ElectricDiscretization
