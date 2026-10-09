#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <initializer_list>
#include <utility>
#include <vector>

#include "Core/Geometry/Detail/LevelSetHermiteCell.h"

// Single-cell extension of the common Hermite scalar, not a mesh/Snapshot or
// another reconstruction. The caller must supply canonical shared-edge data.
// No cubic Bernstein certificate may be used for this quartic representation.
// Keep the same TU-local, strict floating-point implementation boundary.
namespace Quantum::Geometry::Detail
{
namespace
{
namespace LevelSetHermiteScalar
{
class QuarticCell
{
  public:
    // Tensor index 5*y+x. Per axis: value(0), value(1), derivative(0),
    // derivative(1), value(1/2). Derivatives are in LOCAL cell coordinates:
    // hx*phi_x, hy*phi_y, hx*hy*phi_xy. Thus the 16 corner data are supplemented
    // by four edge values, four edge normal derivatives and one center value.
    // Own the supplied doubles; exact queries never use the normalized cache.
    using Dofs = std::array<double, 25>;
    using ExpandedPolynomial = std::array<ExactExpansion::Expansion, 9>;
    using Interval = ExactExpansion::Interval;
    using ExpandedControls = std::array<ExactExpansion::Expansion, 25>;

    const Dofs &localDofs() const noexcept { return raw_; }

    struct RationalBasisJet
    {
        std::array<ExactExpansion::Expansion,5> value, derivative, second;
        ExactExpansion::Expansion denominator;
    };
    // Same basis as scalar/root restriction, evaluated BEFORE a retained
    // material coordinate is rounded. Both orders share positive D^4.
    // This is a terminal-integration adapter, not another interpolation.
    static RationalBasisJet rationalBasisJet(const ExactExpansion::Expansion& numerator,
        const ExactExpansion::Expansion& denominator,bool secondDerivative=false)
    {
        ExactExpansion::validateEnvironment();
        RationalSegmentCoordinate point{numerator,numerator,denominator};
        validateCoordinate(point);
        // Clearing denominators are arbitrary. Normalize them by an exact
        // common power of two before D^4, as in the mapped restriction port;
        // an extreme but harmless chart scale must not overflow/underflow.
        int exponent=0;
        std::frexp(ExactExpansion::estimate(denominator),&exponent);
        point.numeratorBegin=point.numeratorEnd=ExactExpansion::scalePow2(numerator,-exponent);
        point.denominator=ExactExpansion::scalePow2(denominator,-exponent);
        const auto values=restrictedBasisCleared(point,0), derivatives=restrictedBasisCleared(point,1);
        const auto squared=ExactExpansion::multiply(point.denominator,point.denominator);
        RationalBasisJet result;
        result.denominator=ExactExpansion::multiply(squared,squared);
        Detail::require(ExactExpansion::sign(result.denominator)>0,"nonpositive quartic basis denominator");
        for (std::size_t k=0;k<5;++k)
        { result.value[k]=values[k][0]; result.derivative[k]=derivatives[k][0]; }
        if (secondDerivative) {
            const auto second=restrictedBasisCleared(point,2);
            for (std::size_t k=0;k<5;++k) result.second[k]=second[k][0];
        }
        return result;
    }

    QuarticCell(const Dofs &data, double hx, double hy) : raw_(data), hx_(hx), hy_(hy)
    {
        ExactExpansion::validateEnvironment();
        Detail::require(std::isfinite(hx_) && hx_ > 0.0 && std::isfinite(hy_) && hy_ > 0.0,
                        "quartic cell requires positive finite dimensions");
        for (double value : raw_)
            Detail::require(std::isfinite(value), "non-finite quartic Hermite datum");
        offset_ = raw_[0];
        for (std::size_t y = 0; y < 5; ++y)
            for (std::size_t x = 0; x < 5; ++x)
            {
                const auto k = 5 * y + x;
                const bool valueDof = (x < 2 || x == 4) && (y < 2 || y == 4);
                normalized_[k] = valueDof
                                     ? Detail::checked(static_cast<long double>(raw_[k]) - offset_,
                                                       "quartic cell offset subtraction overflow")
                                     : raw_[k];
                scale_ = std::max(scale_, std::abs(normalized_[k]));
            }
        if (scale_ != 0.0)
            for (double &value : normalized_)
            {
                const double normalized = value / scale_;
                Detail::require(std::isfinite(normalized) && (value == 0.0 || normalized != 0.0),
                                "quartic cell normalization underflow");
                value = normalized;
            }
    }

    // Floating approximations for metric evaluation, NOT phase-sign authority.
    // Upper-end chart offsets are evaluated by basis reflection, without 1-a.
    Sample approximateJetAtChart(CellChartPoint point) const
    {
        validateQuery(point);
        const auto bx = approximateBasis(point.x), by = approximateBasis(point.y);
        return {Detail::checked(static_cast<long double>(offset_) + localDerivative(bx, by, 0, 0),
                                "non-finite quartic value"),
                physicalDerivative(bx, by, 1, 0), physicalDerivative(bx, by, 0, 1),
                physicalDerivative(bx, by, 2, 0), physicalDerivative(bx, by, 1, 1),
                physicalDerivative(bx, by, 0, 2)};
    }

    ThirdDerivatives approximateThirdDerivativesAtChart(CellChartPoint point) const
    {
        validateQuery(point);
        const auto bx = approximateBasis(point.x), by = approximateBasis(point.y);
        return {physicalDerivative(bx, by, 3, 0), physicalDerivative(bx, by, 2, 1),
                physicalDerivative(bx, by, 1, 2), physicalDerivative(bx, by, 0, 3)};
    }

    // Exact polynomial of Dx^4*Dy^4 * partial_u^dx partial_v^dy phi(u(t),v(t)),
    // u(t)=Nx(t)/Dx, v(t)=Ny(t)/Dy. The positive denominator is the SAME for
    // every derivative order. These are LOCAL derivatives, not physical ones.
    // All nine coefficients go directly to the existing segment root isolator;
    // no rounded monomial conversion, fitting or small-coefficient deletion.
    ExpandedPolynomial segmentPolynomialCleared(const RationalSegmentCoordinate &x,
                                                 const RationalSegmentCoordinate &y,
                                                 unsigned dx = 0, unsigned dy = 0) const
    {
        ExpandedControls data;
        for (std::size_t k=0;k<data.size();++k)
            data[k]=ExactExpansion::fromDouble(raw_[k]);
        return segmentPolynomialCleared(data,x,y,dx,dy);
    }

    // Native shared constraints retain homogeneous DOF numerators. Apply the
    // SAME tensor restriction without rounding them back to cell doubles;
    // the owning data adapter retains the additional positive denominator.
    static ExpandedPolynomial segmentPolynomialCleared(const ExpandedControls& data,
        const RationalSegmentCoordinate& x,const RationalSegmentCoordinate& y,
        unsigned dx=0,unsigned dy=0)
    {
        namespace Exact = ExactExpansion;
        Exact::validateEnvironment();
        Detail::require(dx <= 4 && dy <= 4, "quartic derivative order exceeds degree four");
        validateCoordinate(x);
        validateCoordinate(y);
        const auto bx = restrictedBasisCleared(x, dx), by = restrictedBasisCleared(y, dy);
        ExpandedPolynomial result{};
        // Contract the y tensor index BEFORE multiplying by the x basis.
        // Distributivity is exact in expansion arithmetic: no rounded power
        // coefficients, fitting, tolerance or dropped term is introduced.
        // A dense restriction uses 125+125 products instead of 125+625;
        // scalar cancellations also occur before the second contraction.
        for(std::size_t i=0;i<5;++i) {
            std::array<Exact::Expansion,5> alongY{};
            for(std::size_t j=0;j<5;++j) {
                const auto& datum=data[5*j+i];if(Exact::sign(datum)==0)continue;
                for(std::size_t py=0;py<5;++py)if(Exact::sign(by[j][py])!=0)
                    alongY[py]=Exact::add(alongY[py],Exact::multiply(datum,by[j][py]));
            }
            for(std::size_t py=0;py<5;++py)if(Exact::sign(alongY[py])!=0)
                for(std::size_t px=0;px<5;++px)if(Exact::sign(bx[i][px])!=0)
                    result[px+py]=Exact::add(result[px+py],Exact::multiply(alongY[py],bx[i][px]));
        }
        return result;
    }

    struct HomogeneousRestriction
    {
        ExpandedPolynomial numerator;
        ExactExpansion::Expansion denominator;
    };
    // Restrict the same scalar using a balanced homogeneous representation.
    // Q4 composition raises chart factors to the eighth total power: setting
    // just their leading denominator to one can underflow valid low tails.
    // Balance complete input objects before products and the complete result
    // afterwards. The raw clearing kernel above keeps its D_x^4 D_y^4 contract.
    static HomogeneousRestriction homogeneousRestriction(ExpandedControls data,
        ExactExpansion::Expansion denominator,RationalSegmentCoordinate x,
        RationalSegmentCoordinate y,unsigned dx=0,unsigned dy=0)
    {
        namespace E=ExactExpansion;
        E::centerHomogeneous(data,denominator);
        for(auto* coordinate:{&x,&y}) {
            E::ExponentRange chartRange;
            chartRange.include(coordinate->numeratorBegin);
            chartRange.include(coordinate->numeratorEnd);
            chartRange.include(coordinate->denominator);
            const int chartShift=chartRange.centeredShift();
            coordinate->numeratorBegin=E::scalePow2(coordinate->numeratorBegin,chartShift);
            coordinate->numeratorEnd=E::scalePow2(coordinate->numeratorEnd,chartShift);
            coordinate->denominator=E::scalePow2(coordinate->denominator,chartShift);
        }
        HomogeneousRestriction result{segmentPolynomialCleared(data,x,y,dx,dy),denominator};
        const auto factor=E::multiply(x.denominator,y.denominator);
        for(int k=0;k<4;++k)result.denominator=E::multiply(result.denominator,factor);
        E::centerHomogeneous(result.numerator,result.denominator);
        return result;
    }
    HomogeneousRestriction homogeneousRestriction(const RationalSegmentCoordinate& x,
        const RationalSegmentCoordinate& y,unsigned dx=0,unsigned dy=0) const
    {
        ExpandedControls data;
        for(std::size_t k=0;k<data.size();++k)data[k]=ExactExpansion::fromDouble(raw_[k]);
        return homogeneousRestriction(std::move(data),ExactExpansion::fromDouble(1),x,y,dx,dy);
    }

    // Exact local scalar/derivative at the represented chart point, including
    // offsets smaller than one ULP of the upper endpoint. Sign tests use this
    // expansion, not approximateJetAtChart(). No division: denominators are 1.
    ExactExpansion::Expansion localDerivativeExactAtChart(CellChartPoint point,
                                                          unsigned dx = 0, unsigned dy = 0) const
    {
        validateQuery(point);
        const auto coordinate = [](CellChartCoordinate chart)
        {
            namespace Exact = ExactExpansion;
            const auto one = Exact::fromDouble(1.0), offset = Exact::fromDouble(chart.offset);
            const auto value = chart.fromUpper ? Exact::subtract(one, offset) : offset;
            return RationalSegmentCoordinate{value, value, one};
        };
        return segmentPolynomialCleared(coordinate(point.x), coordinate(point.y), dx, dy)[0];
    }

    // Same authoritative 25 DOFs in the tensor Bernstein basis. One-dimensional
    // controls are f0, f0+d0/4, (16*fm-5*f0-5*f1-d0+d1)/6, f1-d1/4, f1.
    // The exact port clears 12 in EACH direction, so every control has common
    // positive denominator 144. Do not round these controls before sign tests.
    ExpandedControls bernsteinControlsCleared() const
    {
        ExactExpansion::validateEnvironment();
        namespace Exact = ExactExpansion;
        ExpandedControls values;
        for (std::size_t k = 0; k < values.size(); ++k)
            values[k] = Exact::fromDouble(raw_[k]);
        return bernsteinControlsCleared(values);
    }

    static ExpandedControls bernsteinControlsCleared(const ExpandedControls& values)
    {
        ExactExpansion::validateEnvironment();
        namespace Exact=ExactExpansion;
        return transformAxes(values, [](const std::array<Exact::Expansion, 5>& line)
        {
            constexpr int transform[5][5] = {
                {12,0,0,0,0}, {12,0,3,0,0}, {-10,-10,-2,2,32}, {0,12,0,-3,0}, {0,12,0,0,0}};
            std::array<Exact::Expansion, 5> result{};
            for (int row = 0; row < 5; ++row)
                for (int k = 0; k < 5; ++k)
                    if (transform[row][k] != 0)
                        result[row] = Exact::add(result[row],
                            Exact::multiply(line[k], Exact::fromDouble(transform[row][k])));
            return result;
        });
    }

    struct RegularChartBox
    {
        Vector2D lower,upper;
        int phase=0;
        std::array<int,2> partial{{0,0}};
    };
    // Bounded graded tensor partition. Bisect only coordinate intervals of
    // uncertified boxes, keeping a conforming Cartesian product (no hanging
    // chart edges). Every box must have either
    // a strict phase hull or a strict Cartesian partial of the original Q4.
    // No polynomial is fitted on the subboxes. In particular a critical zero
    // cannot acquire a regular chart by merely increasing this work budget.
    static std::vector<RegularChartBox> regularChartPartition(const ExpandedControls& exact)
    {
        std::array<std::array<Interval,25>,3> nets;
        double scale=0;
        for (int k=0;k<25;++k) {
            nets[0][k]=ExactExpansion::enclosure(exact[k]);
            scale=std::max({scale,std::abs(nets[0][k].lo),std::abs(nets[0][k].hi)});
        }
        constexpr int weights[5][4]={{4,0,0,0},{1,3,0,0},{0,2,2,0},{0,0,3,1},{0,0,0,4}};
        for (int axis=0;axis<2;++axis) for (int line=0;line<5;++line) {
            const auto index=[&](int k){return axis?5*k+line:5*line+k;};
            for (int k=0;k<5;++k) {
                ExactExpansion::Expansion value;
                for (int j=0;j<4;++j) if (weights[k][j]) value=ExactExpansion::add(value,
                    ExactExpansion::multiply(ExactExpansion::subtract(exact[index(j+1)],exact[index(j)]),
                        ExactExpansion::fromDouble(weights[k][j])));
                nets[axis+1][index(k)]=ExactExpansion::enclosure(value);
            }
        }
        const double guard=4096*std::numeric_limits<double>::epsilon()*scale;
        std::vector<double> xCuts{0,.5,1},yCuts{0,.5,1};
        // Preserve the former 32^2 maximum retained boxes and the total
        // 2^2+4^2+...+32^2 classification budget. Small regular corner arcs
        // receive local resolution without increasing uniform work limits.
        constexpr std::size_t maximumBoxes=1024,maximumVisits=1364;
        std::size_t visits=0;
        for (;;) {
            const auto nx=xCuts.size()-1,ny=yCuts.size()-1;
            if(nx*ny>maximumBoxes || visits+nx*ny>maximumVisits) return {};
            std::vector<RegularChartBox> boxes;
            std::vector<bool> splitX(nx,false),splitY(ny,false);
            bool resolved=true;
            for (std::size_t y=0;y<ny;++y) for (std::size_t x=0;x<nx;++x) {
                ++visits;
                RegularChartBox box{{xCuts[x],yCuts[y]},{xCuts[x+1],yCuts[y+1]}};
                const auto f=Detail::tensorBernsteinBounds(nets[0],box.lower,box.upper);
                box.phase=f.lo>0?1:f.hi<0?-1:0;
                if (!box.phase) {
                    const auto du=Detail::tensorBernsteinBounds(nets[1],box.lower,box.upper);
                    const auto dv=Detail::tensorBernsteinBounds(nets[2],box.lower,box.upper);
                    box.partial[0]=du.lo>guard?1:du.hi<-guard?-1:0;
                    box.partial[1]=dv.lo>guard?1:dv.hi<-guard?-1:0;
                    if (!box.partial[0] && !box.partial[1]) {
                        resolved=false;splitX[x]=true;splitY[y]=true;
                    }
                    // Prefer the certified partial with the larger conservative
                    // transverse margin. Integration controls its own numerical
                    // error; do not refine topology solely to satisfy a chosen
                    // fixed-order quadrature aspect ratio.
                    const auto minAbs=[](const Interval& d)
                        { return d.lo>0.0?d.lo:(d.hi<0.0?-d.hi:0.0); };
                    const auto maxAbs=[&guard](const Interval& d)
                        { return std::max({std::abs(d.lo),std::abs(d.hi),guard}); };
                    const double rhoU=box.partial[0]?minAbs(du)/maxAbs(dv):-1.0;
                    const double rhoV=box.partial[1]?minAbs(dv)/maxAbs(du):-1.0;
                    if (rhoV>rhoU) box.partial[0]=0;
                    else box.partial[1]=0;
                }
                boxes.push_back(box);
            }
            if (resolved) return boxes;
            const auto refine=[](std::vector<double>& cuts,const std::vector<bool>& split) {
                std::vector<double> next{cuts.front()};
                for(std::size_t i=0;i<split.size();++i) {
                    if(split[i]) {
                        const auto middle=cuts[i]+(cuts[i+1]-cuts[i])/2;
                        if(!(middle>cuts[i] && middle<cuts[i+1])) return false;
                        next.push_back(middle);
                    }
                    next.push_back(cuts[i+1]);
                }
                cuts=std::move(next);return true;
            };
            if(!refine(xCuts,splitX) || !refine(yCuts,splitY)) return {};
        }
    }

    // Certify one signed Cartesian partial only on boxes that can contain a
    // zero. Off-interface turning points must not invalidate a regular arc.
    // Controls represent the exact homogeneous numerator, with a positive
    // common scale. Every accepted box uses outward interval hulls; bounded
    // work exhaustion is inconclusive. No scalar data or roots are changed.
    static std::array<bool,2> zeroSetCartesianPartials(const ExpandedControls& controls,
        Vector2D lower={0,0},Vector2D upper={1,1})
    {
        Detail::require(std::isfinite(lower.x) && std::isfinite(lower.y)
            && std::isfinite(upper.x) && std::isfinite(upper.y)
            && lower.x>=0 && lower.y>=0 && upper.x<=1 && upper.y<=1
            && lower.x<upper.x && lower.y<upper.y,
            "zero-set partial requires a nonempty closed local subbox");
        using Net=std::array<Interval,25>;
        using A=Detail::Outward;
        Net values;
        double scale=0;
        for (std::size_t k=0;k<25;++k) {
            values[k]=ExactExpansion::enclosure(controls[k]);
            scale=std::max({scale,std::abs(values[k].lo),std::abs(values[k].hi)});
        }
        const double guard=4096*std::numeric_limits<double>::epsilon()*scale;
        // Restrict the SAME scalar and its original-cell derivatives. Keep
        // the full-cell resolution guard; a small chart cannot lower it.
        values=Detail::tensorBernsteinRestriction(values,lower,upper);
        const auto hull=[](const Net& net) {
            Interval range{std::numeric_limits<double>::infinity(),-std::numeric_limits<double>::infinity()};
            for (auto x:net) {
                if (!A::finite(x)) return A::whole();
                range.lo=std::min(range.lo,x.lo);range.hi=std::max(range.hi,x.hi);
            }
            return range;
        };
        const auto split=[](const Net& net,int axis) {
            std::array<Net,2> halves;
            for (int line=0;line<5;++line) {
                std::array<Interval,5> work,left,right;
                const auto index=[&](int k){return axis?5*k+line:5*line+k;};
                for (int k=0;k<5;++k) work[k]=net[index(k)];
                left[0]=work[0];right[4]=work[4];
                for (int stage=1;stage<5;++stage) {
                    for (int k=0;k<5-stage;++k)
                        work[k]=A::multiply(A::add(work[k],work[k+1]),A::point(.5));
                    left[stage]=work[0];right[4-stage]=work[4-stage];
                }
                for (int k=0;k<5;++k) {halves[0][index(k)]=left[k];halves[1][index(k)]=right[k];}
            }
            return halves;
        };
        std::array<bool,2> result{{false,false}};
        for (int axis=0;axis<2;++axis) {
            Net derivative;
            // Differentiate degree four, then elevate degree three back to
            // four. Integer weights cancel the factor four exactly.
            constexpr int weights[5][4]={{4,0,0,0},{1,3,0,0},{0,2,2,0},{0,0,3,1},{0,0,0,4}};
            for (int line=0;line<5;++line) {
                const auto index=[&](int k){return axis?5*k+line:5*line+k;};
                for (int k=0;k<5;++k) {
                    ExactExpansion::Expansion value;
                    for (int j=0;j<4;++j) if (weights[k][j]) value=ExactExpansion::add(value,
                        ExactExpansion::multiply(ExactExpansion::subtract(controls[index(j+1)],controls[index(j)]),
                            ExactExpansion::fromDouble(weights[k][j])));
                    derivative[index(k)]=ExactExpansion::enclosure(value);
                }
            }
            derivative=Detail::tensorBernsteinRestriction(derivative,lower,upper);
            for (int sign:{-1,1}) {
                struct Box { Net value,partial;int depth; };
                std::vector<Box> pending{{values,derivative,0}};
                int visited=0;bool certified=true;
                while (!pending.empty()) {
                    if (++visited>16384) {certified=false;break;}
                    const auto box=std::move(pending.back());pending.pop_back();
                    const auto f=hull(box.value);
                    if (f.lo>0 || f.hi<0) continue;
                    const auto d=hull(box.partial);
                    if (sign>0?d.lo>guard:d.hi<-guard) continue;
                    if (box.depth==24) {certified=false;break;}
                    const auto fv=split(box.value,box.depth&1),dv=split(box.partial,box.depth&1);
                    pending.push_back({fv[1],dv[1],box.depth+1});
                    pending.push_back({fv[0],dv[0],box.depth+1});
                }
                if (certified) {result[axis]=true;break;}
            }
        }
        return result;
    }

    // Resolution-scale/Cartesian-partial test only; these approximate values
    // are not the authority for phase, zero corners or closed-cell exclusion.
    std::array<long double, 25> bernsteinControlsApproximate() const
    {
        ExactExpansion::validateEnvironment();
        std::array<long double, 25> values;
        std::copy(raw_.begin(), raw_.end(), values.begin());
        auto controls = transformAxes(values, [](const std::array<long double, 5>& line)
        {
            return std::array<long double, 5>{{line[0], line[0]+line[2]/4.0L,
                (16.0L*line[4]-5.0L*line[0]-5.0L*line[1]-line[2]+line[3])/6.0L,
                line[1]-line[3]/4.0L, line[1]}};
        });
        for (long double value : controls)
            Detail::require(std::isfinite(value), "nonfinite quartic Bernstein approximation");
        return controls;
    }

    // Outward hull of the same 25 controls, optionally restricted to a closed
    // subbox by the EXISTING tensor blossom. Infinite bounds are inconclusive,
    // never a pure-phase certificate. No sampling or cubic surrogate is used.
    Interval cellBounds(Vector2D lower = {0.0,0.0}, Vector2D upper = {1.0,1.0}) const
    {
        ExactExpansion::validateEnvironment();
        Detail::require(std::isfinite(lower.x) && std::isfinite(lower.y)
            && std::isfinite(upper.x) && std::isfinite(upper.y)
            && lower.x >= 0.0 && lower.y >= 0.0 && upper.x <= 1.0 && upper.y <= 1.0
            && lower.x <= upper.x && lower.y <= upper.y,
            "quartic bounds require an ordered closed local subbox");
        return Detail::tensorBernsteinBounds(bernsteinControlsOutward(), lower, upper);
    }

    int cellSignExceptVertices(unsigned allowedCornerMask) const
    {
        ExactExpansion::validateEnvironment();
        Detail::require((allowedCornerMask & ~15u) == 0u, "invalid quartic allowed-corner mask");
        const auto controls = bernsteinControlsCleared();
        std::array<int, 25> signs;
        for (std::size_t k = 0; k < signs.size(); ++k)
            signs[k] = ExactExpansion::sign(controls[k]);
        return CellView::certifyBernsteinSigns(signs, allowedCornerMask);
    }

    // Outward de Casteljau evaluation at the represented chart point. Reflect
    // control order for an upper-end chart; never first form rounded 1-offset.
    // If the enclosure straddles zero, use localDerivativeExactAtChart().
    Interval evaluateEnclosureAtChart(CellChartPoint point) const
    {
        validateQuery(point);
        using Arithmetic = Detail::Outward;
        const auto evaluateLine = [](std::array<Interval, 5> line, CellChartCoordinate coordinate)
        {
            if (coordinate.fromUpper) std::reverse(line.begin(), line.end());
            if (coordinate.offset == 0.0) return line[0];
            if (coordinate.offset == 1.0) return line[4];
            const auto argument = Arithmetic::point(coordinate.offset);
            const auto complement = Arithmetic::add(Arithmetic::point(1.0), Arithmetic::negate(argument));
            for (int stage = 4; stage > 0; --stage)
                for (int k = 0; k < stage; ++k)
                    line[k] = Arithmetic::add(Arithmetic::multiply(complement, line[k]),
                        Arithmetic::multiply(argument, line[k+1]));
            return line[0];
        };
        const auto controls = bernsteinControlsOutward();
        std::array<Interval, 5> yValues;
        for (int j = 0; j < 5; ++j)
        {
            std::array<Interval, 5> line;
            for (int i = 0; i < 5; ++i) line[i] = controls[5*j+i];
            yValues[j] = evaluateLine(line, point.x);
        }
        return evaluateLine(yValues, point.y);
    }

  private:
    Dofs raw_{}, normalized_{};
    double hx_, hy_, offset_ = 0.0, scale_ = 0.0;
    using ApproximateBasis = std::array<std::array<long double, 5>, 4>;
    using ExpandedBasis = std::array<std::array<ExactExpansion::Expansion, 5>, 5>;

    template<class Value, class Transform>
    static std::array<Value, 25> transformAxes(std::array<Value, 25> data, Transform transform)
    {
        for (int axis = 0; axis < 2; ++axis)
            for (int fixed = 0; fixed < 5; ++fixed)
            {
                std::array<Value, 5> line;
                for (int k = 0; k < 5; ++k) line[k] = data[axis ? 5*k+fixed : 5*fixed+k];
                const auto controls = transform(line);
                for (int k = 0; k < 5; ++k) data[axis ? 5*k+fixed : 5*fixed+k] = controls[k];
            }
        return data;
    }

    std::array<Interval, 25> bernsteinControlsOutward() const
    {
        using Arithmetic = Detail::Outward;
        std::array<Interval, 25> values;
        for (std::size_t k = 0; k < values.size(); ++k) values[k] = Arithmetic::point(raw_[k]);
        return transformAxes(values, [](const std::array<Interval, 5>& line)
        {
            const auto scaled = [](Interval value, double factor)
                { return Arithmetic::multiply(value, Arithmetic::point(factor)); };
            auto middle = Arithmetic::add(scaled(line[4],16.0), scaled(line[0],-5.0));
            middle = Arithmetic::add(middle, scaled(line[1],-5.0));
            middle = Arithmetic::add(middle, Arithmetic::negate(line[2]));
            middle = Arithmetic::add(middle, line[3]);
            return std::array<Interval, 5>{{line[0],
                Arithmetic::add(line[0], Arithmetic::dividePositive(line[2],4.0)),
                Arithmetic::dividePositive(middle,6.0),
                Arithmetic::add(line[1], Arithmetic::negate(Arithmetic::dividePositive(line[3],4.0))),line[1]}};
        });
    }

    static void validateQuery(CellChartPoint point)
    {
        ExactExpansion::validateEnvironment();
        validateCellChart(point);
    }

    static void validateCoordinate(const RationalSegmentCoordinate &coordinate)
    {
        namespace Exact = ExactExpansion;
        Detail::require(Exact::sign(coordinate.denominator) > 0, "quartic restriction denominator must be positive");
        for (const auto *n : {&coordinate.numeratorBegin, &coordinate.numeratorEnd})
            Detail::require(Exact::sign(*n) >= 0 && Exact::sign(Exact::subtract(coordinate.denominator, *n)) >= 0,
                            "quartic restriction endpoint outside cell");
    }

    static ApproximateBasis approximateBasis(CellChartCoordinate point)
    {
        // Qi=Hi-Hi(1/2)*B for the four old cubic Hermite functions; Q4=B,
        // B=16*u^2*(1-u)^2. Integer coefficients avoid a second fitting step.
        constexpr std::array<std::array<int, 5>, 5> coefficients{{
            {{1, 0, -11, 18, -8}}, {{0, 0, -5, 14, -8}}, {{0, 1, -4, 5, -2}},
            {{0, 0, 1, -3, 2}}, {{0, 0, 16, -32, 16}}}};
        constexpr std::array<std::size_t, 5> reflected{{1, 0, 3, 2, 4}};
        constexpr std::array<int, 5> parity{{1, 1, -1, -1, 1}};
        ApproximateBasis result{};
        for (unsigned derivative = 0; derivative <= 3; ++derivative)
            for (std::size_t row = 0; row < 5; ++row)
            {
                const auto source = point.fromUpper ? reflected[row] : row;
                long double value = 0.0L;
                for (int power = 4; power >= static_cast<int>(derivative); --power)
                {
                    int coefficient = coefficients[source][power];
                    for (unsigned order = 0; order < derivative; ++order)
                        coefficient *= power - static_cast<int>(order);
                    value = value * point.offset + coefficient;
                }
                result[derivative][row] = point.fromUpper
                                              ? parity[row] * (derivative % 2 ? -value : value)
                                              : value;
            }
        return result;
    }

    long double localDerivative(const ApproximateBasis &bx, const ApproximateBasis &by,
                                unsigned dx, unsigned dy) const
    {
        long double sum = 0.0L;
        for (std::size_t y = 0; y < 5; ++y)
            for (std::size_t x = 0; x < 5; ++x)
                sum += normalized_[5 * y + x] * bx[dx][x] * by[dy][y];
        const long double value = scale_ * sum;
        Detail::require(std::isfinite(value) && (sum == 0.0L || value != 0.0L),
                        "quartic metric scaling overflow or underflow");
        return value;
    }

    double physicalDerivative(const ApproximateBasis &bx, const ApproximateBasis &by,
                              unsigned dx, unsigned dy) const
    {
        long double value = localDerivative(bx, by, dx, dy);
        for (unsigned order = 0; order < dx + dy; ++order)
        {
            const long double next = value / (order < dx ? hx_ : hy_);
            Detail::require(std::isfinite(next) && (value == 0.0L || next != 0.0L),
                            "quartic physical derivative overflow or underflow");
            value = next;
        }
        const double result = Detail::checked(value, "non-finite quartic physical derivative");
        Detail::require(value == 0.0L || result != 0.0, "quartic derivative conversion underflow");
        return result;
    }

    static ExpandedBasis restrictedBasisCleared(const RationalSegmentCoordinate &coordinate,
                                                 unsigned derivativeOrder)
    {
        namespace Exact = ExactExpansion;
        Exact::validateEnvironment();
        // Pure basis data, independent of any scalar, leaf, chart owner or
        // time level. Identical expansion COMPONENTS (not merely equal values)
        // preserve the original exact-operation ordering on a cache hit.
        // Thread-local bounded storage needs no synchronization or mesh-life
        // coupling; failed computations never enter the table.
        struct Entry
        {
            RationalSegmentCoordinate coordinate;
            unsigned derivativeOrder = 0;
            ExpandedBasis basis;
            bool ready = false;
        };
        struct Cache
        {
            std::array<Entry, 8> entries;
            std::size_t next = 0;
        };
        static thread_local Cache cache;
        const auto sameCoordinate = [&](const RationalSegmentCoordinate &other) {
            return coordinate.denominator.components() == other.denominator.components() &&
                coordinate.numeratorBegin.components() == other.numeratorBegin.components() &&
                coordinate.numeratorEnd.components() == other.numeratorEnd.components();
        };
        for (const auto &entry : cache.entries)
            if (entry.ready && derivativeOrder == entry.derivativeOrder && sameCoordinate(entry.coordinate))
                return entry.basis;
        const auto result = computeRestrictedBasisCleared(coordinate, derivativeOrder);
        // Assemble before assignment: allocation failure cannot leave an old
        // ready entry with new coordinates and incomplete basis coefficients.
        Entry prepared{coordinate, derivativeOrder, result, true};
        cache.entries[cache.next] = std::move(prepared);
        cache.next = (cache.next + 1) % cache.entries.size();
        return result;
    }

    static ExpandedBasis computeRestrictedBasisCleared(const RationalSegmentCoordinate &coordinate,
                                                        unsigned derivativeOrder)
    {
        namespace Exact = ExactExpansion;
        // Reuse the old exact cubic restriction, lifting its D^3 denominator
        // to D^4, then add the differentiated midpoint bubble. No sampling.
        ExpandedBasis result{}, powers{};
        if (derivativeOrder <= 3)
        {
            const auto cubic = CellView::restrictedBasisCleared(coordinate, derivativeOrder);
            for (std::size_t row = 0; row < 4; ++row)
                for (std::size_t k = 0; k < 4; ++k)
                    result[row][k] = Exact::multiply(coordinate.denominator, cubic[row][k]);
        }
        std::array<int, 5> bubble{{0, 0, 16, -32, 16}};
        for (unsigned order = 0; order < derivativeOrder; ++order)
        {
            for (std::size_t k = 0; k < 4; ++k)
                bubble[k] = static_cast<int>(k + 1) * bubble[k + 1];
            bubble[4] = 0;
        }
        std::array<Exact::Expansion, 5> denominators{};
        denominators[0] = powers[0][0] = Exact::fromDouble(1.0);
        const auto delta = Exact::subtract(coordinate.numeratorEnd, coordinate.numeratorBegin);
        const std::size_t parameterDegree = Exact::sign(delta) == 0 ? 0 : 4;
        for (std::size_t power = 1; power <= 4; ++power)
        {
            denominators[power] = Exact::multiply(denominators[power - 1], coordinate.denominator);
            for (std::size_t k = 0; k <= std::min(power, parameterDegree); ++k)
            {
                if (k < power)
                    powers[power][k] = Exact::multiply(coordinate.numeratorBegin, powers[power - 1][k]);
                if (k > 0)
                    powers[power][k] = Exact::add(powers[power][k], Exact::multiply(delta, powers[power - 1][k - 1]));
            }
        }
        for (std::size_t power = 0; power <= 4; ++power)
        {
            if (bubble[power] == 0)
                continue;
            const auto factor = Exact::multiply(Exact::fromDouble(bubble[power]), denominators[4 - power]);
            for (std::size_t k = 0; k <= std::min(power, parameterDegree); ++k)
                result[4][k] = Exact::add(result[4][k], Exact::multiply(factor, powers[power][k]));
        }
        constexpr std::array<double, 4> midpoint{{0.5, 0.5, 0.125, -0.125}};
        for (std::size_t row = 0; row < 4; ++row)
            for (std::size_t k = 0; k < 5; ++k)
                result[row][k] = Exact::subtract(result[row][k], Exact::multiply(Exact::fromDouble(midpoint[row]), result[4][k]));
        return result;
    }
};
} // namespace LevelSetHermiteScalar
} // namespace
} // namespace Quantum::Geometry::Detail
