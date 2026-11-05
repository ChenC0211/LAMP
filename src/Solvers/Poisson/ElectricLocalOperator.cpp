#include "Solvers/Poisson/ElectricLocalOperator.h"

#include <Eigen/Cholesky>
#include <Eigen/QR>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

namespace LAMP::ElectricDiscretization {
namespace {

constexpr double twoPi = 6.283185307179586476925286766559;
using CellMatrix = Eigen::Matrix<double, 6, 6>;
using PotentialMatrix = Eigen::Matrix<double, 10, 10>;
using AugmentedMatrix = Eigen::Matrix<double, 11, 11>;
using FaceMatrix = Eigen::Matrix<double, 3, 3>;
using CellPotentialMatrix = Eigen::Matrix<double, 6, 10>;

void require(bool condition, const char* message)
{
    if (!condition)
        throw std::runtime_error(std::string("[ElectricDiscretization/LocalOperator] ") + message);
}

bool finite(Vector2D point)
{
    return std::isfinite(point.x) && std::isfinite(point.y);
}

std::array<double, 2> normalizedPoint(Vector2D point, Vector2D center, double scale)
{
    require(finite(point) && finite(center), "nonfinite polynomial coordinate");
    require(std::isfinite(scale) && scale > 0.0, "polynomial scale must be positive");
    const double x = (point.x - center.x) / scale;
    const double y = (point.y - center.y) / scale;
    require(std::isfinite(x) && std::isfinite(y), "unrepresentable normalized coordinate");
    return {{x, y}};
}

template<class Derived>
void requireFinite(const Eigen::MatrixBase<Derived>& matrix, const char* message)
{
    require(matrix.allFinite(), message);
}

template<class Matrix>
void validateCholesky(const Eigen::LLT<Matrix>& factor, const char* message)
{
    require(factor.info() == Eigen::Success, message);
    const auto diagonal = factor.matrixL().toDenseMatrix().diagonal().eval();
    require(diagonal.allFinite() && (diagonal.array() > 0.0).all(), message);
}

void checkSymmetric(const Eigen::MatrixXd& matrix, const char* message)
{
    require(matrix.rows() == matrix.cols() && matrix.allFinite(), message);
    const double magnitude = matrix.cwiseAbs().maxCoeff();
    const double defect = (matrix - matrix.transpose()).cwiseAbs().maxCoeff();
    const double budget = 256.0 * std::numeric_limits<double>::epsilon()
        * static_cast<double>(matrix.rows()) * magnitude;
    require(defect <= budget, message);
}

void validateFunction(const ScalarFunction& function)
{
    require(static_cast<bool>(function), "empty quadrature callback");
}

double sample(const ScalarFunction& function, Vector2D point)
{
    const double value = function(point);
    require(std::isfinite(value), "nonfinite quadrature callback value");
    return value;
}

Vector2D callbackPoint(const LocalGeometry& geometry, Vector2D point)
{
    if (geometry.coordinateOrigin.x == 0.0 && geometry.coordinateOrigin.y == 0.0)
        return point;
    const auto world = geometry.coordinateOrigin + point;
    require(finite(world), "unrepresentable world callback coordinate");
    return world;
}

} // namespace

CellVector p2Basis(Vector2D point, Vector2D center, double scale)
{
    const auto xy = normalizedPoint(point, center, scale);
    const double x = xy[0], y = xy[1];
    CellVector result;
    result << 1.0, x, y, x*x, x*y, y*y;
    requireFinite(result, "nonfinite P2 basis");
    return result;
}

PotentialVector p3Basis(Vector2D point, Vector2D center, double scale)
{
    const auto xy = normalizedPoint(point, center, scale);
    const double x = xy[0], y = xy[1];
    PotentialVector result;
    result << 1.0, x, y, x*x, x*y, y*y, x*x*x, x*x*y, x*y*y, y*y*y;
    requireFinite(result, "nonfinite P3 basis");
    return result;
}

PotentialGradient p3Gradient(Vector2D point, Vector2D center, double scale)
{
    const auto xy = normalizedPoint(point, center, scale);
    const double x = xy[0], y = xy[1];
    PotentialGradient result;
    result << 0.0, 1.0, 0.0, 2.0*x, y, 0.0, 3.0*x*x, 2.0*x*y, y*y, 0.0,
              0.0, 0.0, 1.0, 0.0, x, 2.0*y, 0.0, x*x, 2.0*x*y, 3.0*y*y;
    result /= scale;
    requireFinite(result, "nonfinite physical P3 gradient");
    return result;
}

PotentialHessian p3Hessian(Vector2D point, Vector2D center, double scale)
{
    const auto xy = normalizedPoint(point, center, scale);
    const double x = xy[0], y = xy[1];
    PotentialHessian result;
    result[0] << 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 6.0*x, 2.0*y, 0.0, 0.0;
    result[1] << 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 2.0*x, 2.0*y, 0.0;
    result[2] << 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 2.0*x, 6.0*y;
    for (auto& derivative : result) {
        derivative /= scale;
        derivative /= scale;
        requireFinite(derivative, "nonfinite physical P3 Hessian");
    }
    return result;
}

FacePolynomial faceBasis(double coordinate)
{
    require(std::isfinite(coordinate) && coordinate >= -1.0 && coordinate <= 1.0,
            "face coordinate must be finite and in [-1,1]");
    FacePolynomial result;
    result << 1.0, coordinate, 0.5 * (3.0 * coordinate * coordinate - 1.0);
    return result;
}

double physicalWeight(Vector2D point, bool axisymmetric)
{
    require(finite(point), "nonfinite physical quadrature coordinate");
    if (!axisymmetric)
        return 1.0;
    require(point.x >= 0.0, "negative radius in RZ quadrature");
    const double result = twoPi * point.x;
    require(std::isfinite(result), "unrepresentable RZ physical weight");
    return result;
}

struct LocalOperator::Data {
    LocalGeometry geometry;
    double epsilon = 0.0, eta = 0.0;
    Eigen::MatrixXd reconstruction, stiffness, condensed, stabilization;
    PotentialMatrix gradient = PotentialMatrix::Zero();
    Eigen::LLT<CellMatrix> cellMassFactor, cellStiffnessFactor;
    std::vector<Eigen::LLT<FaceMatrix>> faceMassFactors;
    Eigen::MatrixXd cellFromFaces;

    Data(LocalGeometry input, double permittivity, double penalty)
        : geometry(std::move(input)), epsilon(permittivity), eta(penalty)
    {
        require(std::isfinite(epsilon) && epsilon > 0.0, "permittivity must be positive");
        require(std::isfinite(eta) && eta > 0.0, "stabilization must be positive");
        require(finite(geometry.center) && std::isfinite(geometry.scale)
                    && geometry.scale > 0.0, "invalid local center or scale");
        require(finite(geometry.coordinateOrigin), "nonfinite coordinate origin");
        require(!geometry.axisymmetric || (geometry.coordinateOrigin.x == 0.0
                    && geometry.coordinateOrigin.y == 0.0),
                "RZ local operator requires a zero coordinate origin");
        require(!geometry.volume.empty(), "volume quadrature is empty");
        require(!geometry.faces.empty(), "skeletal face list is empty");
        require(geometry.faces.size() <= static_cast<std::size_t>(
                    (std::numeric_limits<int>::max() - 6) / 3), "too many skeletal faces");
        const Eigen::Index faceDofs = static_cast<Eigen::Index>(3 * geometry.faces.size());
        const Eigen::Index localDofs = 6 + faceDofs;

        CellMatrix cellMass = CellMatrix::Zero();
        CellPotentialMatrix cellCross = CellPotentialMatrix::Zero();
        PotentialMatrix geometricGradient = PotentialMatrix::Zero();
        PotentialVector mean = PotentialVector::Zero();
        Eigen::MatrixXd coupling = Eigen::MatrixXd::Zero(10, localDofs);
        double weightedVolume = 0.0;

        for (const auto& node : geometry.volume) {
            require(finite(node.point) && std::isfinite(node.weight) && node.weight > 0.0,
                    "volume nodes need finite coordinates and positive ordinary weights");
            const double m = physicalWeight(node.point, geometry.axisymmetric);
            const double weight = node.weight * m;
            require(std::isfinite(weight) && (m == 0.0 || weight > 0.0),
                    "invalid weighted volume node");
            const auto cell = p2Basis(node.point, geometry.center, geometry.scale);
            const auto potential = p3Basis(node.point, geometry.center, geometry.scale);
            const auto gradientBasis = p3Gradient(node.point, geometry.center, geometry.scale);
            const auto hessianBasis = p3Hessian(node.point, geometry.center, geometry.scale);
            cellMass.noalias() += node.weight * (cell * cell.transpose());
            cellCross.noalias() += node.weight * (cell * potential.transpose());
            geometricGradient.noalias() += weight
                * (gradientBasis.transpose() * gradientBasis);
            mean.noalias() += weight * potential;
            weightedVolume += weight;
            // div(m grad chi)=m*(chi_xx+chi_yy)+(dm/dx)*chi_x.
            // The last term is essential in RZ, including at r=0; no 1/r is used.
            PotentialVector divergence = m * (hessianBasis[0] + hessianBasis[2]);
            if (geometry.axisymmetric)
                divergence.noalias() += twoPi * gradientBasis.row(0).transpose();
            coupling.leftCols(6).noalias() -= node.weight
                * (divergence * cell.transpose());
        }
        require(std::isfinite(weightedVolume) && weightedVolume > 0.0,
                "physical volume must be positive");
        requireFinite(cellMass, "nonfinite ordinary cell mass");
        requireFinite(cellCross, "nonfinite cell potential cross mass");
        requireFinite(geometricGradient, "nonfinite reconstruction gradient Gram matrix");
        requireFinite(mean, "nonfinite reconstruction mean");
        cellMassFactor.compute(cellMass);
        validateCholesky(cellMassFactor, "ordinary cell P2 mass is rank deficient");
        const CellPotentialMatrix cellProjector = cellMassFactor.solve(cellCross);
        requireFinite(cellProjector, "failed ordinary cell projection");

        // Cache each face's weighted projector, preserving the supplied face
        // orientation and arclength coordinate. Axis faces have no unknowns.
        std::vector<Eigen::Matrix<double, 3, 10>> faceProjectors;
        faceProjectors.reserve(geometry.faces.size());
        faceMassFactors.reserve(geometry.faces.size());
        for (std::size_t face = 0; face < geometry.faces.size(); ++face) {
            const auto& quadrature = geometry.faces[face];
            require(std::isfinite(quadrature.length) && quadrature.length > 0.0
                        && !quadrature.nodes.empty(), "invalid face length or empty face quadrature");
            FaceMatrix mass = FaceMatrix::Zero();
            Eigen::Matrix<double, 3, 10> cross = Eigen::Matrix<double, 3, 10>::Zero();
            double ordinaryLength = 0.0, physicalLength = 0.0;
            bool axisFace = geometry.axisymmetric;
            for (const auto& node : quadrature.nodes) {
                require(finite(node.point) && finite(node.outwardNormal)
                            && std::isfinite(node.weight) && node.weight > 0.0,
                        "face nodes need finite coordinates, normal and positive ordinary weights");
                const double normalSquared = node.outwardNormal.x * node.outwardNormal.x
                    + node.outwardNormal.y * node.outwardNormal.y;
                require(std::isfinite(normalSquared) && std::abs(normalSquared - 1.0)
                            <= 256.0 * std::numeric_limits<double>::epsilon(),
                        "face outward normal must be unit length");
                const auto basis = faceBasis(node.coordinate);
                const auto potential = p3Basis(node.point, geometry.center, geometry.scale);
                const auto gradientBasis = p3Gradient(node.point, geometry.center, geometry.scale);
                const double m = physicalWeight(node.point, geometry.axisymmetric);
                const double weight = node.weight * m;
                require(std::isfinite(weight) && (m == 0.0 || weight > 0.0),
                        "invalid weighted face node");
                axisFace = axisFace && node.point.x == 0.0;
                ordinaryLength += node.weight;
                physicalLength += weight;
                mass.noalias() += weight * (basis * basis.transpose());
                cross.noalias() += weight * (basis * potential.transpose());
                const PotentialVector normalDerivative =
                    node.outwardNormal.x * gradientBasis.row(0).transpose()
                    + node.outwardNormal.y * gradientBasis.row(1).transpose();
                coupling.middleCols(static_cast<Eigen::Index>(6 + 3*face), 3).noalias()
                    += weight * (normalDerivative * basis.transpose());
            }
            require(!axisFace, "RZ axis face must be omitted from the skeleton");
            require(std::isfinite(physicalLength) && physicalLength > 0.0,
                    "face physical measure must be positive");
            const double lengthBudget = 4096.0 * std::numeric_limits<double>::epsilon()
                * static_cast<double>(quadrature.nodes.size())
                * std::max(quadrature.length, ordinaryLength);
            require(std::isfinite(ordinaryLength)
                        && std::abs(ordinaryLength - quadrature.length) <= lengthBudget,
                    "ordinary face quadrature does not integrate the supplied face length");
            requireFinite(mass, "nonfinite weighted face P2 mass");
            requireFinite(cross, "nonfinite weighted face potential cross mass");
            faceMassFactors.emplace_back(mass);
            validateCholesky(faceMassFactors.back(), "weighted face P2 mass is rank deficient");
            faceProjectors.push_back(faceMassFactors.back().solve(cross));
            requireFinite(faceProjectors.back(), "failed weighted face projection");
        }
        requireFinite(coupling, "nonfinite weak Green coupling");

        // Constant epsilon cancels from the reconstruction equations. Scale
        // the two homogeneous blocks for QR, without material-dependent basis
        // normalization or altering their weighted Green/mean definitions.
        const double gradientScale = geometricGradient.cwiseAbs().maxCoeff();
        require(std::isfinite(gradientScale) && gradientScale > 0.0,
                "zero reconstruction gradient Gram matrix");
        const Eigen::Matrix<double, 9, 9> nonconstant =
            geometricGradient.bottomRightCorner<9, 9>() / gradientScale;
        Eigen::LLT<Eigen::Matrix<double, 9, 9>> nonconstantFactor(nonconstant);
        validateCholesky(nonconstantFactor, "reconstruction gradient Gram rank is not nine");

        AugmentedMatrix augmented = AugmentedMatrix::Zero();
        augmented.topLeftCorner<10, 10>() = geometricGradient / gradientScale;
        augmented.topRightCorner<10, 1>() = mean / weightedVolume;
        augmented.bottomLeftCorner<1, 10>() = mean.transpose() / weightedVolume;
        Eigen::MatrixXd rhs = Eigen::MatrixXd::Zero(11, localDofs);
        rhs.topRows(10) = coupling / gradientScale;
        rhs.bottomLeftCorner(1, 6) = mean.head<6>().transpose() / weightedVolume;
        Eigen::ColPivHouseholderQR<AugmentedMatrix> reconstructionFactor(augmented);
        require(reconstructionFactor.rank() == 11, "augmented reconstruction is rank deficient");
        const Eigen::MatrixXd solved = reconstructionFactor.solve(rhs);
        requireFinite(solved, "failed augmented reconstruction solve");
        reconstruction = solved.topRows(10);
        gradient = epsilon * geometricGradient;
        requireFinite(gradient, "nonfinite permittivity weighted gradient Gram matrix");

        // delta_T = pi_T(p-V_T), with the ordinary cell projector.
        // delta_TF = pi_F,m(p-V_F). Its second term is the actual
        // pointwise trace of delta_T, not a second face projection.
        Eigen::MatrixXd cellCorrection = cellProjector * reconstruction;
        cellCorrection.leftCols(6) -= CellMatrix::Identity();
        stabilization = Eigen::MatrixXd::Zero(localDofs, localDofs);
        for (std::size_t face = 0; face < geometry.faces.size(); ++face) {
            const auto& quadrature = geometry.faces[face];
            Eigen::MatrixXd faceCorrection = faceProjectors[face] * reconstruction;
            faceCorrection.middleCols(static_cast<Eigen::Index>(6 + 3*face), 3)
                -= FaceMatrix::Identity();
            const double penaltyScale = eta * epsilon / quadrature.length;
            require(std::isfinite(penaltyScale) && penaltyScale > 0.0,
                    "unrepresentable face stabilization scale");
            for (const auto& node : quadrature.nodes) {
                const double weight = node.weight * physicalWeight(node.point, geometry.axisymmetric);
                const auto faceValues = faceBasis(node.coordinate);
                const auto cellValues = p2Basis(node.point, geometry.center, geometry.scale);
                const Eigen::RowVectorXd delta = faceValues.transpose() * faceCorrection
                    - cellValues.transpose() * cellCorrection;
                stabilization.noalias() += (penaltyScale * weight) * (delta.transpose() * delta);
            }
        }
        requireFinite(stabilization, "nonfinite stabilization matrix");
        stiffness = reconstruction.transpose() * gradient * reconstruction + stabilization;
        checkSymmetric(stiffness, "local energy matrix is not symmetric within arithmetic precision");
        const CellMatrix cellBlock = stiffness.topLeftCorner<6, 6>();
        cellStiffnessFactor.compute(cellBlock);
        validateCholesky(cellStiffnessFactor, "cell energy block is not positive definite");
        cellFromFaces = cellStiffnessFactor.solve(stiffness.topRightCorner(6, faceDofs));
        requireFinite(cellFromFaces, "failed cell static condensation");
        condensed = stiffness.bottomRightCorner(faceDofs, faceDofs)
            - stiffness.bottomLeftCorner(faceDofs, 6) * cellFromFaces;
        checkSymmetric(condensed, "cell condensed matrix is not symmetric within arithmetic precision");
    }

    void validateFaces(const Eigen::VectorXd& faces) const
    {
        require(faces.size() == static_cast<Eigen::Index>(3 * geometry.faces.size())
                    && faces.allFinite(), "invalid local face coefficient vector");
    }

    CellVector solveCell(const CellVector& rhs) const
    {
        requireFinite(rhs, "nonfinite cell load");
        const CellVector result = cellStiffnessFactor.solve(rhs);
        requireFinite(result, "nonfinite cell recovery");
        return result;
    }
};

LocalOperator::LocalOperator(LocalGeometry geometry, double permittivity, double stabilization)
    : data_(std::make_shared<const Data>(std::move(geometry), permittivity, stabilization))
{
}

const LocalOperator::Data& LocalOperator::data() const
{
    require(static_cast<bool>(data_), "moved-from local operator");
    return *data_;
}

const LocalGeometry& LocalOperator::geometry() const { return data().geometry; }
std::size_t LocalOperator::faceCount() const { return data().geometry.faces.size(); }
double LocalOperator::permittivity() const { return data().epsilon; }
double LocalOperator::stabilization() const { return data().eta; }
const Eigen::MatrixXd& LocalOperator::reconstructionMatrix() const { return data().reconstruction; }
const Eigen::MatrixXd& LocalOperator::stiffnessMatrix() const { return data().stiffness; }
const Eigen::MatrixXd& LocalOperator::condensedMatrix() const { return data().condensed; }
const Eigen::Matrix<double, 10, 10>& LocalOperator::gradientGramMatrix() const { return data().gradient; }
const Eigen::MatrixXd& LocalOperator::stabilizationMatrix() const { return data().stabilization; }

CellVector LocalOperator::recoverCell(const Eigen::VectorXd& faces, const CellVector& cellLoad) const
{
    const auto& prepared = data();
    prepared.validateFaces(faces);
    const CellVector result = prepared.solveCell(cellLoad) - prepared.cellFromFaces * faces;
    requireFinite(result, "nonfinite recovered cell polynomial");
    return result;
}

PotentialVector LocalOperator::reconstruct(const CellVector& cell, const Eigen::VectorXd& faces) const
{
    const auto& prepared = data();
    requireFinite(cell, "nonfinite cell coefficients");
    prepared.validateFaces(faces);
    const PotentialVector result = prepared.reconstruction.leftCols(6) * cell
        + prepared.reconstruction.rightCols(faces.size()) * faces;
    requireFinite(result, "nonfinite reconstructed potential");
    return result;
}

Eigen::VectorXd LocalOperator::condensedRhs(const CellVector& cellLoad,
                                           const Eigen::VectorXd& faceLoad) const
{
    const auto& prepared = data();
    prepared.validateFaces(faceLoad);
    const Eigen::VectorXd result = faceLoad
        - prepared.stiffness.bottomLeftCorner(faceLoad.size(), 6) * prepared.solveCell(cellLoad);
    requireFinite(result, "nonfinite condensed right hand side");
    return result;
}

Eigen::VectorXd LocalOperator::faceResidual(const CellVector& cell, const Eigen::VectorXd& faces,
                                           const Eigen::VectorXd& faceLoad) const
{
    const auto& prepared = data();
    requireFinite(cell, "nonfinite cell coefficients");
    prepared.validateFaces(faces);
    prepared.validateFaces(faceLoad);
    const Eigen::VectorXd result = prepared.stiffness.bottomLeftCorner(faces.size(), 6) * cell
        + prepared.stiffness.bottomRightCorner(faces.size(), faces.size()) * faces - faceLoad;
    requireFinite(result, "nonfinite face residual");
    return result;
}

CellVector LocalOperator::reductionCell(const ScalarFunction& value) const
{
    const auto& prepared = data();
    validateFunction(value);
    CellVector rhs = CellVector::Zero();
    for (const auto& node : prepared.geometry.volume)
        rhs.noalias() += (node.weight * sample(value, callbackPoint(prepared.geometry, node.point)))
            * p2Basis(node.point, prepared.geometry.center, prepared.geometry.scale);
    const CellVector result = prepared.cellMassFactor.solve(rhs);
    requireFinite(result, "nonfinite ordinary cell reduction");
    return result;
}

Eigen::VectorXd LocalOperator::reductionFaces(const ScalarFunction& value) const
{
    const auto& prepared = data();
    validateFunction(value);
    Eigen::VectorXd result(static_cast<Eigen::Index>(3 * prepared.geometry.faces.size()));
    for (std::size_t face = 0; face < prepared.geometry.faces.size(); ++face) {
        FacePolynomial rhs = FacePolynomial::Zero();
        for (const auto& node : prepared.geometry.faces[face].nodes)
            rhs.noalias() += (node.weight * physicalWeight(node.point, prepared.geometry.axisymmetric)
                * sample(value, callbackPoint(prepared.geometry, node.point))) * faceBasis(node.coordinate);
        result.segment<3>(static_cast<Eigen::Index>(3*face)) = prepared.faceMassFactors[face].solve(rhs);
    }
    requireFinite(result, "nonfinite weighted face reduction");
    return result;
}

Eigen::VectorXd LocalOperator::reduction(const ScalarFunction& value) const
{
    Eigen::VectorXd result(static_cast<Eigen::Index>(6 + 3 * faceCount()));
    result.head<6>() = reductionCell(value);
    result.tail(result.size() - 6) = reductionFaces(value);
    return result;
}

CellVector LocalOperator::loadCell(const ScalarFunction& source) const
{
    const auto& prepared = data();
    validateFunction(source);
    CellVector result = CellVector::Zero();
    for (const auto& node : prepared.geometry.volume)
        result.noalias() += (node.weight * physicalWeight(node.point, prepared.geometry.axisymmetric)
            * sample(source, callbackPoint(prepared.geometry, node.point)))
            * p2Basis(node.point, prepared.geometry.center, prepared.geometry.scale);
    requireFinite(result, "nonfinite physical cell load");
    return result;
}

FacePolynomial LocalOperator::loadFace(std::size_t face, const ScalarFunction& flux) const
{
    const auto& prepared = data();
    require(face < prepared.geometry.faces.size(), "face load index out of range");
    validateFunction(flux);
    FacePolynomial result = FacePolynomial::Zero();
    for (const auto& node : prepared.geometry.faces[face].nodes)
        result.noalias() += (node.weight * physicalWeight(node.point, prepared.geometry.axisymmetric)
            * sample(flux, callbackPoint(prepared.geometry, node.point))) * faceBasis(node.coordinate);
    requireFinite(result, "nonfinite physical face load");
    return result;
}

} // namespace LAMP::ElectricDiscretization
