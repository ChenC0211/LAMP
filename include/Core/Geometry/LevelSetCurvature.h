#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "Core/Geometry/LevelSetPolynomialReconstruction.h"

#include "Core/Base/Vector2D.h"
#include "Core/Field/VertField.h"
#include "Core/Mesh/UniformLeafGrid.h"

// Fixed-uniform local scalar reconstruction. Method selection belongs to the
// caller; evaluation never changes the level set or interface node positions.
namespace LevelSetCurvature
{
using Basis = LevelSetPolynomial::Basis;
using Derivatives = LevelSetPolynomial::Derivatives;

struct Sample
{
    double dx = 0.0, dy = 0.0;
    double dxx = 0.0, dxy = 0.0, dyy = 0.0;
    double nx = 0.0, ny = 0.0;
    // In-plane -div(grad(phi)/|grad(phi)|). Any RZ -nx/r term is the
    // caller's responsibility; no coordinate-system convention is inferred.
    double curvature = 0.0;
};

namespace Detail
{
inline void require(bool condition, const char* reason)
{
    if (!condition) throw std::runtime_error(std::string("[LevelSetCurvature] ") + reason);
}

inline double checked(long double value, const char* reason)
{
    const double result = static_cast<double>(value);
    require(std::isfinite(value) && std::isfinite(result), reason);
    return result;
}
}

// Immutable snapshot: construct again after phi or mesh changes. The complete
// mesh is validated once, and nodal values are copied into Cartesian order;
// evaluate() retains no references and does no mesh search or QR factorization.
//
// A fixed 5x5 stencil fits all monomials of total degree 3 (10 coefficients)
// or 4 (15 coefficients), with equal sample weights. The separate TensorProduct
// degree-4 basis uses all x^u*y^v, 0<=u,v<=4 (25 coefficients): the square
// system is biquartic interpolation of the same 25 nodal values, not the
// total-degree least-squares fit. All modes share the same QR and derivative
// path. Its anchor is the nearest grid vertex, clamped to keep all 25 vertices
// inside the domain by default. Explicit radialAxisSymmetry instead uses an
// even extension of actual nodal values at r=0 and sets odd radial derivatives
// there to zero; outer-radial and axial boundaries retain their old stencils.
// Physical coordinates are normalized independently by hx
// and hy. No mode, degree or stencil is adapted to curvature, a circle radius,
// or measured error. Tensor-product interpolation does not promise smoothing
// of noisy phi or continuity between different stencil anchors.
// RootOnly preserves the original complete root-grid contract. UniformLeaves
// explicitly permits a complete equal-level leaf grid, not hanging-node AMR
// or an active-fluid subset. Neither scope changes the mesh's identities.
// An explicit native-index sample mask restricts ONLY the derivative stencil:
// excluded phi payload is never read. The unchanged 5x5 polynomial is shifted
// to the nearest wholly admitted stencil containing the query (no extrapolation).
// This does not create an active-domain Hermite scalar or a boundary condition.
template<class Mesh>
class Reconstruction
{
public:
    Reconstruction(const Mesh& mesh, const VertField<double>& phi, int degree = 4,
        Basis basis = Basis::TotalDegree, bool radialAxisSymmetry = false,
        CartesianGridScope scope = CartesianGridScope::RootOnly,
        const std::vector<std::uint8_t>& sampleMask = {})
        : grid_(validatedGrid(mesh, scope)), degree_(degree), nx_(grid_.nx()), ny_(grid_.ny()), basis_(basis), lower_(grid_.minCoord()),
          upper_(grid_.maxCoord()), hx_(grid_.spacing().x), hy_(grid_.spacing().y),
          radialAxisSymmetry_(radialAxisSymmetry)
    {
        using Detail::require;
        require(!radialAxisSymmetry_ || lower_.x == 0.0,
            "radial even extension requires the actual r=0 domain boundary");
        require(basis_ == Basis::TotalDegree || basis_ == Basis::TensorProduct, "unsupported polynomial basis");
        require(degree_ == 3 || degree_ == 4, "total polynomial degree must be 3 or 4");
        require(basis_ != Basis::TensorProduct || degree_ == 4,
            "the fixed 5x5 tensor-product basis requires degree 4");
        require(nx_ >= 4 && ny_ >= 4, "a 5x5 stencil requires Nx and Ny >= 4");
        const std::int64_t vertexCount = (static_cast<std::int64_t>(nx_) + 1) * (static_cast<std::int64_t>(ny_) + 1);
        require(phi.value && phi.size() == vertexCount, "nodal phi size must match the mesh");
        require(sampleMask.empty() || sampleMask.size() == static_cast<std::size_t>(vertexCount),
            "derivative sample mask must use native vertex indices");
        values_.resize(static_cast<std::size_t>(vertexCount));
        if (!sampleMask.empty()) sampleAllowed_.resize(static_cast<std::size_t>(vertexCount));
        for (std::size_t position = 0; position < grid_.vertices().size(); ++position)
        {
            // The layout owns actual labels, whereas phi uses native indices.
            const int k = mesh.verts(grid_.vertices()[position]).index();
            if (!sampleMask.empty())
            {
                require(sampleMask[k] <= 1, "derivative sample mask must be binary");
                sampleAllowed_[position] = sampleMask[k];
                if (!sampleMask[k]) continue;
            }
            require(std::isfinite(phi.value[k]), "nodal phi contains a non-finite value");
            values_[position] = phi.value[k];
        }
        polynomial_ = LevelSetPolynomial::Reconstruction(degree_, basis_);
        condition_ = polynomial_.conditionEstimate();
    }

    // Retained metadata lets a scalar built from these SAME Q4 derivatives
    // reuse the validated grid rather than invent a second root-grid map.
    const UniformLeafGrid& grid() const { return grid_; }
    int degree() const { return degree_; }
    Basis basis() const { return basis_; }
    bool radialAxisSymmetry() const { return radialAxisSymmetry_; }
    // Induced infinity norms ||A||_inf ||A^+||_inf for the normalized,
    // rectangular monomial design matrix and its QR-derived left inverse.
    double conditionEstimate() const { return condition_; }

    // Raw local derivatives are meaningful at constants and stationary
    // points too. This does not expose a globally continuous scalar value;
    // the existing nearest-anchor polynomial/stencil selection is unchanged.
    Derivatives evaluateDerivatives(const Vector2D& point) const
    {
        long double variation = 0.0L;
        return computeDerivatives(point, variation);
    }

    // Value and derivatives of the SAME selected 5x5 polynomial, for shared
    // Hermite midpoint/center data generation. Not a continuous global phi.
    LevelSetPolynomial::ScalarJet evaluateJet(const Vector2D& point) const
    {
        long double variation = 0.0L;
        LevelSetPolynomial::ScalarJet result;
        result.derivatives = computeDerivatives(point, variation, &result.value);
        return result;
    }

    // Half-grid point with at least one odd doubled index (edge midpoint or
    // cell center). Average ALL equally nearest admitted Q4 anchors, including
    // derivatives from those same fits. Integer distances preserve exact ties
    // without a physical-coordinate rounding decision. The original nearest-
    // anchor ports above are unchanged. This is DOF generation, not phase or
    // a globally continuous value/derivative query. Boundary/mask admission
    // uses the existing support rules; there is no solid extrapolation.
    LevelSetPolynomial::ScalarJet evaluateMidpointJet(int twiceX, int twiceY) const
    {
        using Detail::require;
        require(twiceX >= 0 && twiceY >= 0 && static_cast<std::int64_t>(twiceX) <= 2LL*nx_
            && static_cast<std::int64_t>(twiceY) <= 2LL*ny_ && ((twiceX&1) || (twiceY&1)),
            "midpoint query requires an in-domain half-grid index");
        const double x = twiceX/2.0, y = twiceY/2.0;
        int nearest = 33;
        std::vector<std::array<int,2>> anchors;
        for (int j = std::max(2, static_cast<int>(std::ceil(y-2.0)));
            j <= std::min(ny_-2, static_cast<int>(std::floor(y+2.0))); ++j)
            for (int i = std::max(radialAxisSymmetry_ ? 0 : 2, static_cast<int>(std::ceil(x-2.0)));
                i <= std::min(nx_-2, static_cast<int>(std::floor(x+2.0))); ++i)
            {
                if (!sampleAllowed_.empty() && !stencilAllowed(i,j)) continue;
                const auto dx = 2LL*i-twiceX, dy = 2LL*j-twiceY;
                const int distance = static_cast<int>(dx*dx+dy*dy); // offsets bounded by the 5x5 support
                if (distance > nearest) continue;
                if (distance < nearest) { nearest = distance; anchors.clear(); }
                anchors.push_back({{i,j}});
            }
        require(!anchors.empty(), "no admitted 5x5 stencil contains the midpoint; extrapolation is unsupported");
        std::array<long double,6> sums{};
        std::array<double,6> reference{};
        for (std::size_t k = 0; k < anchors.size(); ++k)
        {
            const auto& a = anchors[k];
            long double variation = 0.0L;
            const auto jet = polynomial_.evaluateJet(samplesAt(a[0],a[1]), x-a[0], y-a[1], hx_, hy_,
                radialAxisSymmetry_ && twiceX == 0, variation);
            const auto& d = jet.derivatives;
            const std::array<double,6> values{{jet.value,d.dx,d.dy,d.dxx,d.dxy,d.dyy}};
            if (k == 0) reference = values;
            for (std::size_t c = 0; c < values.size(); ++c)
                sums[c] += static_cast<long double>(values[c])-reference[c];
        }
        std::array<double,6> averaged{};
        for (std::size_t c = 0; c < averaged.size(); ++c)
            averaged[c] = Detail::checked(reference[c]+sums[c]/anchors.size(), "nonfinite shared midpoint jet");
        return {averaged[0],{averaged[1],averaged[2],averaged[3],averaged[4],averaged[5]}};
    }

    // Higher-accuracy DATA for a retained Hermite scalar, not a new curve or
    // curvature query. Integer grid points use seven nodes, half points eight,
    // per direction. Complete support, including its radial even extension,
    // is checked against the SAME sample mask before any payload is read.
    // All equally closest admitted supports are combined; there is no reduced-
    // order fallback or error-dependent support/weight selection.
    LevelSetPolynomial::ScalarJet evaluateHermiteDataJet(int twiceX,int twiceY) const
    {
        using Detail::require;
        require(degree_ == 4 && basis_ == Basis::TensorProduct,
            "enriched Hermite data require Q4 tensor reconstruction");
        require(twiceX >= 0 && twiceY >= 0 && static_cast<std::int64_t>(twiceX) <= 2LL*nx_
            && static_cast<std::int64_t>(twiceY) <= 2LL*ny_,"invalid enriched Hermite half-grid index");
        const int sx=7+(twiceX&1),sy=7+(twiceY&1);
        require(nx_ >= sx-1 && ny_ >= sy-1,"grid too small for complete enriched Hermite stencil");
        const double x=twiceX/2.0,y=twiceY/2.0;
        int nearest=std::numeric_limits<int>::max();
        std::vector<std::array<int,2>> starts;
        for (int j=std::max(0,static_cast<int>(std::ceil(y-(sy-1))));
            j<=std::min(ny_-(sy-1),static_cast<int>(std::floor(y)));++j)
            for (int i=std::max(radialAxisSymmetry_?-nx_:0,static_cast<int>(std::ceil(x-(sx-1))));
                i<=std::min(nx_-(sx-1),static_cast<int>(std::floor(x)));++i)
            {
                bool allowed=true;
                if (!sampleAllowed_.empty())
                    for (int v=0;v<sy && allowed;++v)
                        for (int u=0;u<sx;++u)
                            if (!sampleAllowed_[index(radialAxisSymmetry_?std::abs(i+u):i+u,j+v)])
                            {allowed=false;break;}
                if (!allowed) continue;
                const auto dx=2LL*i+(sx-1)-twiceX,dy=2LL*j+(sy-1)-twiceY;
                const int distance=static_cast<int>(dx*dx+dy*dy); // offsets bounded by the support
                if (distance>nearest) continue;
                if (distance<nearest) {nearest=distance;starts.clear();}
                starts.push_back({{i,j}});
            }
        require(!starts.empty(),"no admitted complete enriched Hermite stencil; extrapolation is unsupported");
        std::array<long double,6> sum{};
        std::array<double,6> reference{};
        for (std::size_t k=0;k<starts.size();++k)
        {
            const auto& start=starts[k];
            std::array<double,64> samples{};
            for (int v=0;v<sy;++v)
                for (int u=0;u<sx;++u)
                    samples[8*v+u]=values_[index(radialAxisSymmetry_?std::abs(start[0]+u):start[0]+u,start[1]+v)];
            const auto jet=polynomial_.evaluateExtendedJet(samples,sx,sy,x-start[0],y-start[1],hx_,hy_,
                radialAxisSymmetry_ && twiceX==0);
            const auto& d=jet.derivatives;
            const std::array<double,6> values{{jet.value,d.dx,d.dy,d.dxx,d.dxy,d.dyy}};
            if (k==0) reference=values;
            for (std::size_t c=0;c<values.size();++c) sum[c]+=static_cast<long double>(values[c])-reference[c];
        }
        std::array<double,6> values{};
        for (std::size_t c=0;c<values.size();++c)
            values[c]=Detail::checked(reference[c]+sum[c]/starts.size(),"nonfinite enriched Hermite data");
        return {values[0],{values[1],values[2],values[3],values[4],values[5]}};
    }

    Sample evaluate(const Vector2D& point) const
    {
        using Detail::require;
        long double variation = 0.0L;
        const Derivatives derivatives = computeDerivatives(point, variation);
        require(variation > 0.0L, "degenerate level-set gradient: constant stencil");
        Sample result;
        result.dx = derivatives.dx; result.dy = derivatives.dy;
        result.dxx = derivatives.dxx; result.dxy = derivatives.dxy; result.dyy = derivatives.dyy;
        const double magnitude = std::hypot(result.dx, result.dy);
        const long double resolutionFloor = 64.0L * std::numeric_limits<double>::epsilon() * condition_
            * variation / std::min(hx_, hy_);
        require(std::isfinite(magnitude) && std::isfinite(resolutionFloor) && magnitude > resolutionFloor,
            "degenerate or roundoff-unresolved level-set gradient");
        result.nx = result.dx / magnitude; result.ny = result.dy / magnitude;
        result.curvature = Detail::checked(-(static_cast<long double>(result.ny) * result.ny * result.dxx
            - 2.0L * result.nx * result.ny * result.dxy
            + static_cast<long double>(result.nx) * result.nx * result.dyy) / magnitude,
            "non-finite reconstructed curvature");
        return result;
    }

private:
    static UniformLeafGrid validatedGrid(const Mesh& mesh, CartesianGridScope scope)
    {
        try { return UniformLeafGrid(mesh, scope); }
        catch (const std::invalid_argument& error)
        {
            // Preserve this evaluator's existing explicit runtime-error port
            // while retaining the shared layout validator's concrete reason.
            throw std::runtime_error(std::string("[LevelSetCurvature] ") + error.what());
        }
    }

    Derivatives computeDerivatives(const Vector2D& point, long double& variation, double* value = nullptr) const
    {
        using Detail::require;
        require(std::isfinite(point.x) && std::isfinite(point.y)
            && point.x >= lower_.x && point.x <= upper_.x
            && point.y >= lower_.y && point.y <= upper_.y, "evaluation point must lie in the finite mesh domain");
        const double xGrid = (point.x - lower_.x) / hx_;
        const double yGrid = (point.y - lower_.y) / hy_;
        require(std::isfinite(xGrid) && std::isfinite(yGrid), "unrepresentable normalized evaluation point");
        int i = std::clamp(static_cast<int>(std::llround(xGrid)), radialAxisSymmetry_ ? 0 : 2, nx_ - 2);
        int j = std::clamp(static_cast<int>(std::llround(yGrid)), 2, ny_ - 2);
        if (!sampleAllowed_.empty() && !stencilAllowed(i, j))
        {
            // Geometry-independent choice. Never replace missing samples by
            // zero/solid values or choose a stencil by its resulting curvature.
            double nearest = std::numeric_limits<double>::infinity();
            int selectedX = -1, selectedY = -1;
            for (int b = std::max(2, static_cast<int>(std::ceil(yGrid - 2.0)));
                b <= std::min(ny_ - 2, static_cast<int>(std::floor(yGrid + 2.0))); ++b)
                for (int a = std::max(radialAxisSymmetry_ ? 0 : 2, static_cast<int>(std::ceil(xGrid - 2.0)));
                    a <= std::min(nx_ - 2, static_cast<int>(std::floor(xGrid + 2.0))); ++a)
                {
                    if (!stencilAllowed(a, b)) continue;
                    const double distance = (a - xGrid) * (a - xGrid) + (b - yGrid) * (b - yGrid);
                    if (distance < nearest) { nearest = distance; selectedX = a; selectedY = b; }
                }
            require(selectedX >= 0, "no admitted 5x5 stencil contains the query; extrapolation is unsupported");
            i = selectedX; j = selectedY;
        }
        const long double x = xGrid - i, y = yGrid - j;
        const auto samples = samplesAt(i,j);
        if (value)
        {
            const auto result = polynomial_.evaluateJet(samples, x, y, hx_, hy_,
                radialAxisSymmetry_ && point.x == 0.0, variation);
            *value = result.value;
            return result.derivatives;
        }
        return polynomial_.evaluateDerivatives(samples, x, y, hx_, hy_,
            radialAxisSymmetry_ && point.x == 0.0, variation);
    }

    std::array<double,25> samplesAt(int i, int j) const
    {
        std::array<double,25> samples{};
        int row = 0;
        for (int v = -2; v <= 2; ++v)
            for (int u = -2; u <= 2; ++u)
                samples[row++] = values_[index(radialAxisSymmetry_ ? std::abs(i+u) : i+u, j+v)];
        return samples;
    }

    UniformLeafGrid grid_;
    int degree_, nx_, ny_;
    Basis basis_;
    Vector2D lower_, upper_;
    double hx_, hy_, condition_ = 0.0;
    bool radialAxisSymmetry_ = false;
    std::vector<double> values_;
    std::vector<std::uint8_t> sampleAllowed_;
    LevelSetPolynomial::Reconstruction polynomial_;

    std::size_t index(int i, int j) const
    { return static_cast<std::size_t>(i) + (static_cast<std::size_t>(nx_) + 1) * j; }

    bool stencilAllowed(int i, int j) const
    {
        for (int v = -2; v <= 2; ++v)
            for (int u = -2; u <= 2; ++u)
                if (!sampleAllowed_[index(radialAxisSymmetry_ ? std::abs(i + u) : i + u, j + v)]) return false;
        return true;
    }


};
}
