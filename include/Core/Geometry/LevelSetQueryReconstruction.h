#pragma once

#include "Core/Base/Vector2D.h"
#include "Core/Geometry/LevelSetPolynomialReconstruction.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace LevelSetPolynomial
{
// The existing native transport's fixed 5x5 physical-query stencil. This is
// neither a replacement mesh nor interpolation of hanging-node field payload.
// The caller supplies an actual scalar and certifies the ENTIRE support box
// against its active domain/chart. Selection depends only on that support,
// never on sampled values, curvature or a target answer. No extrapolation.
class QueryReconstruction
{
    struct Support { int x, y, distance; };

  public:
    // Geometry-only selection. No sampled scalar, callback or mesh is retained.
    // Only prepareHermiteDataStencil can construct one; the caller owns the
    // validity of the admitted domain/chart when reusing it.
    class HermiteDataStencil
    {
        friend class QueryReconstruction;
        Vector2D anchor_, spacing_;
        std::array<int, 2> halfIndex_;
        bool radial_;
        std::vector<Support> selected_;
        HermiteDataStencil(Vector2D anchor, Vector2D spacing, std::array<int, 2> halfIndex, bool radial)
            : anchor_(anchor), spacing_(spacing), halfIndex_(halfIndex), radial_(radial) {}
    };

    explicit QueryReconstruction(bool radial) : radial_(radial), polynomial_(4, Basis::TensorProduct)
    {
    }
    double conditionEstimate() const
    {
        return polynomial_.conditionEstimate();
    }

    template <class Admitted, class Scalar>
    Derivatives evaluateDerivatives(Vector2D point, Vector2D spacing, const Admitted &admitted,
                                    const Scalar &scalar, long double &variation) const
    {
        return compute<false>(point, spacing, admitted, scalar, variation);
    }

    // Native physical-query counterpart of the same polynomial jet port.
    // Samples and support checks are performed once, with the old ordering.
    template <class Admitted, class Scalar>
    ScalarJet evaluateJet(Vector2D point, Vector2D spacing, const Admitted &admitted,
                          const Scalar &scalar, long double &variation) const
    {
        return compute<true>(point, spacing, admitted, scalar, variation);
    }

    // Native DATA-generation counterpart of the enriched grid port. The
    // target is anchor + halfIndex*spacing/2, with each halfIndex 0 or 1.
    // Keep that half-grid identity: do not round a midpoint first and then
    // infer its normalized coordinates. Samples are physical queries of the
    // supplied scalar; no native field or hanging payload is interpolated here.
    // Uses the SAME Q4 Neville extension, complete 7/8-point supports and
    // equally-nearest support average as Reconstruction::evaluateHermiteDataJet.
    template <class Admitted, class Scalar>
    ScalarJet evaluateHermiteDataJet(Vector2D anchor, Vector2D spacing,
        std::array<int,2> halfIndex,const Admitted& admitted,const Scalar& scalar) const
    {
        return evaluateHermiteDataJet(prepareHermiteDataStencil(anchor, spacing, halfIndex, admitted), scalar);
    }

    template <class Admitted>
    HermiteDataStencil prepareHermiteDataStencil(Vector2D anchor, Vector2D spacing,
        std::array<int, 2> halfIndex, const Admitted& admitted) const
    {
        using Detail::require;
        require(std::isfinite(anchor.x) && std::isfinite(anchor.y)
            && std::isfinite(spacing.x) && std::isfinite(spacing.y) && spacing.x>0 && spacing.y>0
            && (!radial_ || anchor.x>=0) && halfIndex[0]>=0 && halfIndex[0]<=1
            && halfIndex[1]>=0 && halfIndex[1]<=1,"invalid native enriched Hermite query");
        const int sizeX=7+halfIndex[0],sizeY=7+halfIndex[1];
        // Test nearest boxes first to avoid visiting distant supports on
        // every native datum. Stable ordering preserves the original y/x
        // summation order for ties. Selection never inspects scalar values.
        HermiteDataStencil result(anchor, spacing, halfIndex, radial_);
        for(const auto& support:hermiteSupportOrder())
        {
            if(!result.selected_.empty() && support.distance>result.selected_.front().distance)break;
            Vector2D lo{anchor.x+support.x*spacing.x,anchor.y+support.y*spacing.y};
            Vector2D hi{anchor.x+(support.x+sizeX-1)*spacing.x,anchor.y+(support.y+sizeY-1)*spacing.y};
            if(radial_ && lo.x<0) {hi.x=std::max(hi.x,-lo.x);lo.x=0;}
            if(!std::isfinite(lo.x) || !std::isfinite(lo.y) || !std::isfinite(hi.x)
                || !std::isfinite(hi.y) || !(lo.x<hi.x && lo.y<hi.y))continue;
            if(admitted(lo,hi))result.selected_.push_back(support);
        }
        require(!result.selected_.empty(),"no admitted complete native enriched Hermite stencil; extrapolation is unsupported");
        return result;
    }

    template <class Scalar>
    ScalarJet evaluateHermiteDataJet(const HermiteDataStencil& stencil, const Scalar& scalar) const
    {
        return computeHermiteDataJet(stencil, [&](std::size_t k, int u, int v) {
            return scalar(hermiteSamplePoint(stencil, k, u, v));
        });
    }

    // Same coordinate expression and y/x/support order as the scalar port.
    // Enables an owning batch to prepare sample indices without fitting any
    // dummy scalar or duplicating the recursive interpolation kernel.
    template <class Visitor>
    void visitHermiteDataSamples(const HermiteDataStencil& stencil, const Visitor& visitor) const
    {
        Detail::require(radial_ == stencil.radial_, "prepared Hermite stencil radial policy mismatch");
        for (std::size_t k=0; k<stencil.selected_.size(); ++k)
            for (int v=0; v<7+stencil.halfIndex_[1]; ++v)
                for (int u=0; u<7+stencil.halfIndex_[0]; ++u)
                    visitor(hermiteSamplePoint(stencil, k, u, v));
    }

    // Values in visitHermiteDataSamples order; the caller owns/count-checks
    // its sample array. Arithmetic below is shared with the physical port.
    template <class IndexedScalar>
    ScalarJet evaluateHermiteDataSamples(const HermiteDataStencil& stencil, const IndexedScalar& scalar) const
    {
        std::size_t index=0;
        return computeHermiteDataJet(stencil, [&](std::size_t, int, int) { return scalar(index++); });
    }

  private:
    static Vector2D hermiteSamplePoint(const HermiteDataStencil& stencil, std::size_t k, int u, int v)
    {
        const auto& start=stencil.selected_[k];
        Vector2D point{stencil.anchor_.x+(start.x+u)*stencil.spacing_.x,
            stencil.anchor_.y+(start.y+v)*stencil.spacing_.y};
        if (stencil.radial_) point.x=std::abs(point.x);
        return point;
    }

    template <class Sample>
    ScalarJet computeHermiteDataJet(const HermiteDataStencil& stencil, const Sample& sample) const
    {
        Detail::require(radial_ == stencil.radial_, "prepared Hermite stencil radial policy mismatch");
        const auto anchor = stencil.anchor_, spacing = stencil.spacing_;
        const auto halfIndex = stencil.halfIndex_;
        const int sizeX=7+halfIndex[0],sizeY=7+halfIndex[1];
        const double x=halfIndex[0]*.5,y=halfIndex[1]*.5;
        const auto& selected = stencil.selected_;
        std::array<long double,6> sum{};
        std::array<double,6> reference{};
        for(std::size_t k=0;k<selected.size();++k)
        {
            const auto& start=selected[k];
            std::array<double,64> samples{};
            for(int v=0;v<sizeY;++v)for(int u=0;u<sizeX;++u)
            {
                samples[8*v+u]=sample(k,u,v);
            }
            const auto jet=polynomial_.evaluateExtendedJet(samples,sizeX,sizeY,x-start.x,y-start.y,
                spacing.x,spacing.y,radial_ && anchor.x==0 && halfIndex[0]==0);
            const auto& d=jet.derivatives;
            const std::array<double,6> values{{jet.value,d.dx,d.dy,d.dxx,d.dxy,d.dyy}};
            if(k==0)reference=values;
            for(std::size_t c=0;c<values.size();++c)sum[c]+=static_cast<long double>(values[c])-reference[c];
        }
        std::array<double,6> values{};
        for(std::size_t c=0;c<values.size();++c)
            values[c]=Detail::checked(reference[c]+sum[c]/selected.size(),"nonfinite native enriched Hermite data");
        return {values[0],{values[1],values[2],values[3],values[4],values[5]}};
    }

  private:
    static const std::vector<Support>& hermiteSupportOrder()
    {
        // For both halfIndex=0,size=7 and halfIndex=1,size=8, the old
        // candidate starts are [-6,0] and doubled center distance is 2*i+6.
        // This immutable ordering is independent of coordinates or fields.
        static const std::vector<Support> order = [] {
            std::vector<Support> result;
            result.reserve(49);
            for (int j=-6; j<=0; ++j)
                for (int i=-6; i<=0; ++i)
                    result.push_back({i, j, (2*i+6)*(2*i+6)+(2*j+6)*(2*j+6)});
            std::stable_sort(result.begin(), result.end(),
                [](const auto& a, const auto& b) { return a.distance < b.distance; });
            return result;
        }();
        return order;
    }

    template <bool WithValue, class Admitted, class Scalar>
    auto compute(Vector2D point, Vector2D spacing, const Admitted &admitted,
                 const Scalar &scalar, long double &variation) const
    {
        if (!std::isfinite(point.x) || !std::isfinite(point.y) || !std::isfinite(spacing.x) ||
            !std::isfinite(spacing.y) || spacing.x <= 0 || spacing.y <= 0 || (radial_ && point.x < 0))
            throw std::invalid_argument("[LevelSetQueryReconstruction] invalid point or physical spacing");
        int selectedX = 0, selectedY = 0, nearest = 9;
        bool found = false;
        for (int sy = -2; sy <= 2; ++sy)
            for (int sx = -2; sx <= 2; ++sx)
            {
                const int distance = sx * sx + sy * sy;
                if (distance >= nearest)
                    continue;
                const Vector2D center{point.x + sx * spacing.x, point.y + sy * spacing.y};
                if (radial_ && center.x < 0)
                    continue;
                Vector2D lo{point.x + (sx - 2) * spacing.x, point.y + (sy - 2) * spacing.y};
                Vector2D hi{point.x + (sx + 2) * spacing.x, point.y + (sy + 2) * spacing.y};
                if (radial_ && lo.x < 0)
                {
                    hi.x = std::max(hi.x, -lo.x);
                    lo.x = 0;
                }
                if (!admitted(lo, hi))
                    continue;
                nearest = distance;
                selectedX = sx;
                selectedY = sy;
                found = true;
            }
        if (!found)
            throw std::runtime_error("[LevelSetQueryReconstruction] no admitted native 5x5 query stencil "
                                     "contains the vertex; extrapolation is unsupported");
        std::array<double, 25> samples{};
        int row = 0;
        for (int v = -2; v <= 2; ++v)
            for (int u = -2; u <= 2; ++u)
            {
                Vector2D sample{point.x + (selectedX + u) * spacing.x, point.y + (selectedY + v) * spacing.y};
                if (radial_)
                    sample.x = std::abs(sample.x);
                samples[row++] = scalar(sample);
            }
        if constexpr (WithValue)
            return polynomial_.evaluateJet(samples, -selectedX, -selectedY, spacing.x, spacing.y,
                                           radial_ && point.x == 0.0, variation);
        else
            return polynomial_.evaluateDerivatives(samples, -selectedX, -selectedY, spacing.x, spacing.y,
                                                   radial_ && point.x == 0.0, variation);
    }

  private:
    bool radial_;
    Reconstruction polynomial_;
};
} // namespace LevelSetPolynomial
