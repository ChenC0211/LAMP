#pragma once

#include <Eigen/Core>
#include <Eigen/QR>
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

// Shared fixed 25-point polynomial kernel. The caller owns stencil selection,
// sample meaning and activity/geometry checks; this class owns no mesh or field.
// Arithmetic, QR rank/reproduction gates and coefficient ordering are the
// original LevelSetCurvature reconstruction, used by both grid layouts.
namespace LevelSetPolynomial
{
enum class Basis
{
    TotalDegree,
    TensorProduct
};

struct Derivatives
{
    double dx = 0.0, dy = 0.0;
    double dxx = 0.0, dxy = 0.0, dyy = 0.0;
};

// Value and derivatives of ONE selected polynomial. This is not a globally
// continuous scalar across different stencils and does not define phase signs.
struct ScalarJet
{
    double value = 0.0;
    Derivatives derivatives;
};

namespace Detail
{
inline void require(bool condition, const char *reason)
{
    if (!condition)
        throw std::runtime_error(std::string("[LevelSetPolynomial] ") + reason);
}

inline double checked(long double value, const char *reason)
{
    const double result = static_cast<double>(value);
    require(std::isfinite(value) && std::isfinite(result), reason);
    return result;
}
} // namespace Detail

class Reconstruction
{
  public:
    Reconstruction() = default;
    Reconstruction(int degree, Basis basis) : degree_(degree), basis_(basis)
    {
        Detail::require(basis == Basis::TotalDegree || basis == Basis::TensorProduct,
                        "unsupported polynomial basis");
        Detail::require(degree == 3 || degree == 4, "total polynomial degree must be 3 or 4");
        Detail::require(basis != Basis::TensorProduct || degree == 4,
                        "the fixed 5x5 tensor-product basis requires degree 4");
        prepareProjection();
    }
    double conditionEstimate() const
    {
        return condition_;
    }

    // Samples are y-major over [-2,2]^2, normalized query relative to the
    // chosen anchor. Physical derivative scaling is independent in x and y.
    Derivatives evaluateDerivatives(const std::array<double, 25> &samples, long double x, long double y,
                                    double hx, double hy, bool oddRadialZero, long double &variation) const
    {
        return compute(samples, x, y, hx, hy, oddRadialZero, variation, nullptr);
    }

    ScalarJet evaluateJet(const std::array<double, 25> &samples, long double x, long double y,
                          double hx, double hy, bool oddRadialZero, long double &variation) const
    {
        ScalarJet result;
        result.derivatives = compute(samples, x, y, hx, hy, oddRadialZero, variation, &result.value);
        return result;
    }

    // Interpolate a complete 7/8 by 7/8 tensor stencil by Neville extension of
    // THIS existing Q4 interpolant. Storage is y-major with fixed row stride 8;
    // entries outside [0,sizeX)x[0,sizeY) are not read. No new QR solve/basis or
    // fitted weights: adjoining degree-four interpolants are recursively
    // combined into degree sizeX-1/sizeY-1. Differentiate the coordinate weights
    // as well as the Q4 polynomials; averaging their derivatives is incorrect.
    ScalarJet evaluateExtendedJet(const std::array<double,64>& samples, int sizeX, int sizeY,
        long double x, long double y, double hx, double hy, bool oddRadialZero) const
    {
        using Detail::require;
        require(terms_ == 25 && degree_ == 4 && basis_ == Basis::TensorProduct,
            "extended interpolation requires the existing Q4 tensor basis");
        require((sizeX == 7 || sizeX == 8) && (sizeY == 7 || sizeY == 8)
            && std::isfinite(x) && std::isfinite(y) && x >= 0 && x <= sizeX-1 && y >= 0 && y <= sizeY-1
            && std::isfinite(hx) && hx > 0 && std::isfinite(hy) && hy > 0,
            "invalid complete extended tensor stencil/query");
        const long double reference = samples[0];
        std::array<double,64> centered{};
        for (int j=0;j<sizeY;++j)
            for (int i=0;i<sizeX;++i)
                centered[8*j+i] = Detail::checked(static_cast<long double>(samples[8*j+i])-reference,
                    "nonfinite extended scalar difference");
        // Remove one COMMON value offset before the Q4 evaluations. The
        // differentiated Neville weights must not act on independently rounded
        // large-offset values; restore that offset once in the terminal value.
        std::array<std::array<ExtendedJet,4>,4> rows{};
        for (int j=0;j<sizeY-4;++j)
            for (int i=0;i<sizeX-4;++i)
            {
                std::array<double,25> local{};
                for (int v=0;v<5;++v)
                    for (int u=0;u<5;++u)
                        local[5*v+u]=centered[8*(j+v)+i+u];
                long double variation=0.0L;
                const auto q4=evaluateJet(local,x-i-2,y-j-2,hx,hy,false,variation);
                const auto& d=q4.derivatives;
                rows[j][i]={{q4.value,d.dx,d.dy,d.dxx,d.dxy,d.dyy}};
            }
        for (int j=0;j<sizeY-4;++j)
            for (int count=6;count<=sizeX;++count)
                for (int i=0;i<=sizeX-count;++i)
                    rows[j][i]=extendJet(rows[j][i],rows[j][i+1],(x-i)/(count-1),
                        1.0L/(count-1)/hx,true);
        for (int count=6;count<=sizeY;++count)
            for (int j=0;j<=sizeY-count;++j)
                rows[j][0]=extendJet(rows[j][0],rows[j+1][0],(y-j)/(count-1),
                    1.0L/(count-1)/hy,false);
        auto result=rows[0][0];
        result[0]+=reference;
        ScalarJet jet;
        jet.value=Detail::checked(result[0],"nonfinite extended value");
        jet.derivatives={Detail::checked(result[1],"nonfinite extended dx"),
            Detail::checked(result[2],"nonfinite extended dy"),Detail::checked(result[3],"nonfinite extended dxx"),
            Detail::checked(result[4],"nonfinite extended dxy"),Detail::checked(result[5],"nonfinite extended dyy")};
        if (oddRadialZero) {jet.derivatives.dx=0.0;jet.derivatives.dxy=0.0;}
        return jet;
    }

  private:
    using ExtendedJet=std::array<long double,6>;
    static ExtendedJet extendJet(const ExtendedJet& left,const ExtendedJet& right,
        long double weight,long double weightDerivative,bool xDirection)
    {
        ExtendedJet result{};
        for (std::size_t k=0;k<result.size();++k)
            result[k]=left[k]+weight*(right[k]-left[k]);
        if (xDirection)
        {
            result[1]+=weightDerivative*(right[0]-left[0]);
            result[3]+=2.0L*weightDerivative*(right[1]-left[1]);
            result[4]+=weightDerivative*(right[2]-left[2]);
        }
        else
        {
            result[2]+=weightDerivative*(right[0]-left[0]);
            result[5]+=2.0L*weightDerivative*(right[2]-left[2]);
            result[4]+=weightDerivative*(right[1]-left[1]);
        }
        return result;
    }
    // Same projection, normalization and derivative accumulation as the old
    // derivative-only port. The optional value shares these coefficients; it
    // neither repeats the fit nor changes the old derivative operation order.
    Derivatives compute(const std::array<double, 25> &samples, long double x, long double y,
                        double hx, double hy, bool oddRadialZero, long double &variation, double *value) const
    {
        using Detail::require;
        require(terms_ > 0 && std::isfinite(x) && std::isfinite(y) && std::isfinite(hx) && hx > 0 &&
                    std::isfinite(hy) && hy > 0,
                "invalid polynomial query");
        const long double reference = samples[12];
        std::array<long double, 25> differences{};
        variation = 0.0L;
        for (int row = 0; row < 25; ++row)
        {
            const long double value = static_cast<long double>(samples[row]) - reference;
            require(std::isfinite(value), "unrepresentable local phi difference");
            differences[row] = value;
            variation = std::max(variation, std::abs(value));
        }
        if (variation == 0.0L)
        {
            if (value)
                *value = samples[12];
            return Derivatives{};
        }
        // Remove the offset and scale the local right-hand side before the
        // projection; restore physical units only in the differentiated result.
        // All 25 projection rows consume the SAME normalized samples. Divide
        // each once; retain the original multiplication/addition order and
        // precision within every row. This caches arithmetic inputs, not a
        // refitted operator or a different polynomial/discretization.
        std::array<long double, 25> normalized{};
        for (int sample = 0; sample < 25; ++sample)
            normalized[sample] = differences[sample] / variation;
        std::array<long double, 25> coefficients{};
        for (int term = 0; term < terms_; ++term)
            for (int sample = 0; sample < 25; ++sample)
                coefficients[term] += projection_[term * 25 + sample] * normalized[sample];
        std::array<long double, 5> xp{{1.0L}}, yp{{1.0L}};
        for (int power = 1; power <= degree_; ++power)
        {
            xp[power] = xp[power - 1] * x;
            yp[power] = yp[power - 1] * y;
        }
        long double dx = 0.0L, dy = 0.0L, dxx = 0.0L, dxy = 0.0L, dyy = 0.0L, scalar = 0.0L;
        for (int term = 0; term < terms_; ++term)
        {
            const int u = powers_[term][0], v = powers_[term][1];
            const long double coefficient = coefficients[term];
            if (value)
                scalar += coefficient * xp[u] * yp[v];
            if (u)
                dx += coefficient * u * xp[u - 1] * yp[v];
            if (v)
                dy += coefficient * v * xp[u] * yp[v - 1];
            if (u >= 2)
                dxx += coefficient * u * (u - 1) * xp[u - 2] * yp[v];
            if (u && v)
                dxy += coefficient * u * v * xp[u - 1] * yp[v - 1];
            if (v >= 2)
                dyy += coefficient * v * (v - 1) * xp[u] * yp[v - 2];
        }
        if (value)
            *value = Detail::checked(reference + variation * scalar, "non-finite reconstructed value");
        Derivatives result;
        result.dx = Detail::checked(dx * variation / hx, "non-finite reconstructed dx");
        result.dy = Detail::checked(dy * variation / hy, "non-finite reconstructed dy");
        result.dxx = Detail::checked(dxx * variation / hx / hx, "non-finite reconstructed dxx");
        result.dxy = Detail::checked(dxy * variation / hx / hy, "non-finite reconstructed dxy");
        result.dyy = Detail::checked(dyy * variation / hy / hy, "non-finite reconstructed dyy");
        if (oddRadialZero)
        {
            // Prescribed odd-derivative identities: do not retain QR roundoff
            // as a nonzero axis normal or divide it by a vanishing radius.
            result.dx = 0.0;
            result.dxy = 0.0;
        }
        return result;
    }

  private:
    int degree_ = 0, terms_ = 0;
    Basis basis_ = Basis::TotalDegree;
    double condition_ = 0.0;
    std::array<std::array<int, 2>, 25> powers_{};
    std::array<double, 25 * 25> projection_{};
    void prepareProjection()
    {
        using Detail::require;
        if (basis_ == Basis::TotalDegree)
        {
            // Preserve the original total-degree column ordering and solve.
            for (int total = 0; total <= degree_; ++total)
                for (int u = total; u >= 0; --u)
                    powers_[terms_++] = {{u, total - u}};
        }
        else
        {
            for (int u = 0; u <= degree_; ++u)
                for (int v = 0; v <= degree_; ++v)
                    powers_[terms_++] = {{u, v}};
        }
        Eigen::MatrixXd design(25, terms_);
        int row = 0;
        for (int y = -2; y <= 2; ++y)
            for (int x = -2; x <= 2; ++x)
            {
                std::array<double, 5> xp{{1.0}}, yp{{1.0}};
                for (int p = 1; p <= degree_; ++p)
                {
                    xp[p] = xp[p - 1] * x;
                    yp[p] = yp[p - 1] * y;
                }
                for (int term = 0; term < terms_; ++term)
                    design(row, term) = xp[powers_[term][0]] * yp[powers_[term][1]];
                ++row;
            }
        Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(design);
        qr.setThreshold(64.0 * std::numeric_limits<double>::epsilon());
        require(qr.rank() == terms_, "rank-deficient polynomial design matrix");
        const Eigen::MatrixXd projection = qr.solve(Eigen::MatrixXd::Identity(25, 25));
        require(projection.allFinite(), "non-finite QR projection matrix");
        condition_ =
            design.cwiseAbs().rowwise().sum().maxCoeff() * projection.cwiseAbs().rowwise().sum().maxCoeff();
        require(std::isfinite(condition_) && condition_ >= 1.0 && condition_ <= 1e8,
                "ill-conditioned normalized polynomial design matrix");
        const double defect =
            (projection * design - Eigen::MatrixXd::Identity(terms_, terms_)).cwiseAbs().maxCoeff();
        require(std::isfinite(defect) &&
                    defect <= 256.0 * std::numeric_limits<double>::epsilon() * condition_,
                "QR projection fails polynomial reproduction check");
        for (int term = 0; term < terms_; ++term)
            for (int sample = 0; sample < 25; ++sample)
                projection_[term * 25 + sample] = projection(term, sample);
    }
};
} // namespace LevelSetPolynomial
