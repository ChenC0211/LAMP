#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <initializer_list>
#include <mutex>
#include <utility>
#include <vector>

#include "Core/Geometry/Detail/LevelSetHermiteCell.h"

// Single-cell extension of the common Hermite scalar, not a mesh/Snapshot or
// another reconstruction. The caller must supply canonical shared-edge data.
// No cubic Bernstein certificate may be used for this quartic representation.
// Keep the same TU-local, strict floating-point implementation boundary.
namespace LAMP::Geometry::Detail
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
    double approximateValueAtChart(CellChartPoint point) const
    {
        validateQuery(point);
        const auto bx=approximateBasis(point.x,0),by=approximateBasis(point.y,0);
        return approximateValue(bx,by);
    }

    Sample approximateJetAtChart(CellChartPoint point) const
    {
        validateQuery(point);
        const auto bx = approximateBasis(point.x), by = approximateBasis(point.y);
        return {approximateValue(bx,by),
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
        // A transverse line with constant x has only bx[i][0]. Contract that
        // axis first: sum_j by[j](t) * sum_i data[j,i]*bx[i][0]. This is the
        // SAME exact tensor polynomial and denominator, with no fitted powers
        // or skipped root. Dense work is 25+25 rather than 125+25 products.
        // Component equality is a sufficient exact-constant test, never a
        // rounded coordinate key. General and constant-y paths retain their
        // original operation order below.
        if(x.numeratorBegin.components()==x.numeratorEnd.components()
            && y.numeratorBegin.components()!=y.numeratorEnd.components()) {
            for(std::size_t j=0;j<5;++j) {
                Exact::Expansion alongX;
                for(std::size_t i=0;i<5;++i) {
                    const auto& datum=data[5*j+i];
                    if(Exact::sign(datum)!=0 && Exact::sign(bx[i][0])!=0)
                        Exact::accumulate(alongX,Exact::multiply(datum,bx[i][0]));
                }
                if(Exact::sign(alongX)!=0)
                    for(std::size_t py=0;py<5;++py)if(Exact::sign(by[j][py])!=0)
                        Exact::accumulate(result[py],Exact::multiply(alongX,by[j][py]));
            }
            return result;
        }
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
                    Exact::accumulate(alongY[py],Exact::multiply(datum,by[j][py]));
            }
            for(std::size_t py=0;py<5;++py)if(Exact::sign(alongY[py])!=0)
                for(std::size_t px=0;px<5;++px)if(Exact::sign(bx[i][px])!=0)
                    Exact::accumulate(result[px+py],Exact::multiply(alongY[py],bx[i][px]));
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
        return homogeneousRestrictionCentered(data,denominator,std::move(x),std::move(y),dx,dy);
    }
    // The supplied immutable data are the EXACT result of the original first
    // centerHomogeneous call. A leaf may retain that result once; coordinates,
    // output scaling and all restriction/FP checks still run for every query.
    static HomogeneousRestriction homogeneousRestrictionCentered(const ExpandedControls& data,
        const ExactExpansion::Expansion& denominator,RationalSegmentCoordinate x,
        RationalSegmentCoordinate y,unsigned dx=0,unsigned dy=0)
    {
        namespace E=ExactExpansion;
        E::validateEnvironment();
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

    static ApproximateBasis approximateBasis(CellChartCoordinate point,unsigned maximumDerivative=3)
    {
        // Qi=Hi-Hi(1/2)*B for the four old cubic Hermite functions; Q4=B,
        // B=16*u^2*(1-u)^2. Integer coefficients avoid a second fitting step.
        constexpr std::array<std::array<int, 5>, 5> coefficients{{
            {{1, 0, -11, 18, -8}}, {{0, 0, -5, 14, -8}}, {{0, 1, -4, 5, -2}},
            {{0, 0, 1, -3, 2}}, {{0, 0, 16, -32, 16}}}};
        constexpr std::array<std::size_t, 5> reflected{{1, 0, 3, 2, 4}};
        constexpr std::array<int, 5> parity{{1, 1, -1, -1, 1}};
        ApproximateBasis result{};
        for (unsigned derivative = 0; derivative <= maximumDerivative; ++derivative)
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

    double approximateValue(const ApproximateBasis& bx,const ApproximateBasis& by) const
    {
        return Detail::checked(static_cast<long double>(offset_) + localDerivative(bx,by,0,0),
            "non-finite quartic value");
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

// Four quartic pieces of ONE scalar on a native cell. The 49 physical datum
// are supplied by the immutable shared completion, not new transported DOFs.
// Original Q4 inputs/queries above remain a separate representation.
class C2MacroQuarticCell
{
    using Exact = ExactExpansion::Expansion;
  public:
    using Dofs = std::array<double,49>;
    using Interval = ExactExpansion::Interval;

    C2MacroQuarticCell(const Dofs& physical,Vector2D lower,Vector2D upper)
        : physical_(physical),lower_(lower),upper_(upper)
    {
        ExactExpansion::validateEnvironment();
        for(double coordinate:{lower.x,lower.y,upper.x,upper.y})
            Detail::require(std::isfinite(coordinate),"C2 macro requires finite native bounds");
        Detail::require(lower.x<upper.x && lower.y<upper.y,
            "C2 macro requires positive actual native dimensions");
        for(double value:physical_)
            Detail::require(std::isfinite(value),"nonfinite C2 macro physical datum");
        width_[0]=difference(upper.x,lower.x);
        width_[1]=difference(upper.y,lower.y);
        constexpr std::array<int,7> orders{{0,0,0,1,1,2,2}};
        std::array<std::array<Exact,3>,2> scale;
        for(int axis=0;axis<2;++axis) {
            scale[axis][0]=number(1);
            scale[axis][1]=width_[axis];
            scale[axis][2]=ExactExpansion::multiply(width_[axis],width_[axis]);
        }
        std::array<Exact,49> native;
        for(int y=0;y<7;++y)for(int x=0;x<7;++x)
            native[7*y+x]=ExactExpansion::multiply(
                ExactExpansion::multiply(number(physical_[7*y+x]),scale[0][orders[x]]),
                scale[1][orders[y]]);
        for(int py=0;py<2;++py)for(int px=0;px<2;++px) {
            auto& piece=pieces_[2*py+px];
            std::array<Exact,35> alongX;
            for(int y=0;y<7;++y)for(int x=0;x<5;++x)
                for(int k=0;k<7;++k)if(pieceDatumNumerator[px][x][k]!=0)
                    ExactExpansion::accumulate(alongX[5*y+x],
                        ExactExpansion::multiply(native[7*y+k],
                            number(pieceDatumNumerator[px][x][k])));
            for(int y=0;y<5;++y)for(int x=0;x<5;++x)
                for(int k=0;k<7;++k)if(pieceDatumNumerator[py][y][k]!=0)
                    ExactExpansion::accumulate(piece.local[5*y+x],
                        ExactExpansion::multiply(alongX[5*k+x],
                            number(pieceDatumNumerator[py][y][k])));
            // Fixed rational seven-basis transform, denominator 768 per axis.
            // Exact native widths normalize all derivative orders. No rounded
            // piece25 cache is subsequently used as scalar authority.
            piece.denominator=number(768.0*768.0);
            ExactExpansion::centerHomogeneous(piece.local,piece.denominator);
            auto controls=QuarticCell::bernsteinControlsCleared(piece.local);
            auto denominator=ExactExpansion::multiply(piece.denominator,number(144));
            ExactExpansion::centerHomogeneous(controls,denominator);
            Detail::require(denominator.components().size()==1,
                "C2 macro Bernstein denominator is not exactly represented");
            const double divisor=ExactExpansion::estimate(denominator);
            for(std::size_t k=0;k<controls.size();++k)
                piece.bounds[k]=Detail::Outward::dividePositive(
                    ExactExpansion::enclosure(controls[k]),divisor);
            // Remove the true native corner value before estimating the
            // variation scale, so a large scalar offset cannot cancel it.
            for(const auto& control:controls) {
                const auto relative=ExactExpansion::subtract(control,
                    ExactExpansion::multiply(number(physical_[0]),denominator));
                const long double value=ratio(relative,denominator);
                controlLower_=std::min(controlLower_,value);
                controlUpper_=std::max(controlUpper_,value);
            }
        }
        const long double minimumWidth=std::min(
            static_cast<long double>(ExactExpansion::estimate(width_[0])),
            static_cast<long double>(ExactExpansion::estimate(width_[1])));
        resolutionScale_=2*(controlUpper_-controlLower_)/minimumWidth;
        Detail::require(std::isfinite(resolutionScale_),
            "unrepresentable C2 macro derivative-resolution scale");
    }

    const Dofs& physicalDofs() const noexcept { return physical_; }
    long double derivativeResolutionScale() const noexcept { return resolutionScale_; }

    Sample approximateJetAtChart(CellChartPoint point) const
    { return jet(selectNative(point)); }
    double approximateValueAtChart(CellChartPoint point) const
    {
        const auto selected=selectNative(point);
        return derivative(selected,0,0);
    }
    Sample approximateJetAtPieceChart(int piece,CellChartPoint point) const
    { return jet(explicitPiece(piece,point)); }
    Sample approximateJetAtPhysicalPoint(Vector2D point) const
    { return jet(selectPhysical(point)); }
    double approximateValueAtPhysicalPoint(Vector2D point) const
    { return derivative(selectPhysical(point),0,0); }
    ThirdDerivatives approximateThirdDerivativesAtChart(CellChartPoint point) const
    {
        const auto selected=selectNative(point);
        return {derivative(selected,3,0),derivative(selected,2,1),
            derivative(selected,1,2),derivative(selected,0,3)};
    }

    // Full linear derivative of this piece scalar with respect to the 49
    // PHYSICAL completed data. The same rational seven-basis transform and
    // actual native widths used to form the net are transposed here.
    std::array<long double,49> differentialWeights(CellChartPoint point,
        const std::array<long double,6>& weights) const
    {
        for(long double weight:weights)
            Detail::require(std::isfinite(weight),"nonfinite C2 macro differential weight");
        const auto selected=selectNative(point);
        const int half[2]={selected.piece%2,selected.piece/2};
        const auto localX=QuarticCell::rationalBasisJet(
            selected.x.numeratorBegin,selected.x.denominator,true);
        const auto localY=QuarticCell::rationalBasisJet(
            selected.y.numeratorBegin,selected.y.denominator,true);
        const QuarticCell::RationalBasisJet* local[2]={&localX,&localY};
        std::array<std::array<std::array<Exact,7>,3>,2> basis;
        const auto channel=[](const QuarticCell::RationalBasisJet& jet,int order)
            ->const std::array<Exact,5>& {
            return order==0?jet.value:order==1?jet.derivative:jet.second;
        };
        for(int axis=0;axis<2;++axis)for(int derivative=0;derivative<3;++derivative)
            for(int datum=0;datum<7;++datum)for(int k=0;k<5;++k) {
                const int factor=pieceDatumNumerator[half[axis]][k][datum];
                if(factor!=0)ExactExpansion::accumulate(basis[axis][derivative][datum],
                    ExactExpansion::multiply(channel(*local[axis],derivative)[k],number(factor)));
            }
        const auto basisDenominator=ExactExpansion::multiply(number(768.0*768.0),
            ExactExpansion::multiply(localX.denominator,localY.denominator));
        constexpr int inputOrder[7]={0,0,0,1,1,2,2};
        constexpr int xOrder[6]={0,1,0,2,1,0},yOrder[6]={0,0,1,0,1,2};
        std::array<long double,49> result{};
        for(int y=0;y<7;++y)for(int x=0;x<7;++x)for(int c=0;c<6;++c) {
            if(weights[c]==0)continue;
            auto numerator=ExactExpansion::multiply(basis[0][xOrder[c]][x],
                basis[1][yOrder[c]][y]);
            if(ExactExpansion::sign(numerator)==0)continue;
            auto denominator=basisDenominator;
            numerator=ExactExpansion::scalePow2(numerator,xOrder[c]+yOrder[c]);
            for(int k=0;k<inputOrder[x];++k)numerator=ExactExpansion::multiply(numerator,width_[0]);
            for(int k=0;k<inputOrder[y];++k)numerator=ExactExpansion::multiply(numerator,width_[1]);
            for(int k=0;k<xOrder[c];++k)denominator=ExactExpansion::multiply(denominator,width_[0]);
            for(int k=0;k<yOrder[c];++k)denominator=ExactExpansion::multiply(denominator,width_[1]);
            std::array<Exact,1> normalized{{std::move(numerator)}};
            ExactExpansion::centerHomogeneous(normalized,denominator);
            const long double coefficient=static_cast<long double>(ExactExpansion::estimate(normalized[0]))
                /ExactExpansion::estimate(denominator);
            const long double term=weights[c]*coefficient,sum=result[7*y+x]+term;
            Detail::require(std::isfinite(coefficient)&&coefficient!=0
                &&std::isfinite(term)&&term!=0&&std::isfinite(sum),
                "unrepresentable C2 macro differential coefficient");
            result[7*y+x]=sum;
        }
        return result;
    }

    // Outward bounds use the exact four-piece net and the existing tensor
    // blossom. This introduces no root/region capability or surrogate field.
    Interval cellBounds(Vector2D lower={0,0},Vector2D upper={1,1}) const
    {
        ExactExpansion::validateEnvironment();
        Detail::require(std::isfinite(lower.x)&&std::isfinite(lower.y)
            &&std::isfinite(upper.x)&&std::isfinite(upper.y)
            &&lower.x>=0&&lower.y>=0&&upper.x<=1&&upper.y<=1
            &&lower.x<=upper.x&&lower.y<=upper.y,
            "C2 macro bounds require an ordered closed native subbox");
        Interval result{std::numeric_limits<double>::infinity(),
            -std::numeric_limits<double>::infinity()};
        for(int py=0;py<2;++py)for(int px=0;px<2;++px) {
            const Vector2D lo{std::max(lower.x,px*.5),std::max(lower.y,py*.5)};
            const Vector2D hi{std::min(upper.x,(px+1)*.5),std::min(upper.y,(py+1)*.5)};
            if(lo.x>hi.x||lo.y>hi.y)continue;
            const auto range=Detail::tensorBernsteinBounds(pieces_[2*py+px].bounds,
                {2*lo.x-px,2*lo.y-py},{2*hi.x-px,2*hi.y-py});
            result.lo=std::min(result.lo,range.lo);
            result.hi=std::max(result.hi,range.hi);
        }
        return result;
    }

    // Outward enclosure of the SAME scalar or its physical derivative on
    // an actual native box. Each intersecting piece is differentiated in its
    // own Bernstein net; reduced-degree nets are elevated algebraically back
    // to Q4 before the existing tensor blossom. No exterior extrapolation.
    Interval derivativeBounds(Vector2D lower,Vector2D upper,unsigned dx,unsigned dy) const
    {
        ExactExpansion::validateEnvironment();
        Detail::require(dx<=4&&dy<=4,"C2 macro bound derivative order exceeds Q4");
        Detail::require(std::isfinite(lower.x)&&std::isfinite(lower.y)&&std::isfinite(upper.x)&&std::isfinite(upper.y)
            &&lower.x>=0&&lower.y>=0&&upper.x<=1&&upper.y<=1&&lower.x<=upper.x&&lower.y<=upper.y,
            "C2 macro derivative bounds leave the actual native cell");
        Interval result{std::numeric_limits<double>::infinity(),-std::numeric_limits<double>::infinity()};
        for(int py=0;py<2;++py)for(int px=0;px<2;++px){
            const Vector2D lo{std::max(lower.x,px*.5),std::max(lower.y,py*.5)},
                hi{std::min(upper.x,(px+1)*.5),std::min(upper.y,(py+1)*.5)};
            if(lo.x>hi.x||lo.y>hi.y)continue;
            const auto controls=physicalDerivativeControls(pieces_[2*py+px],dx,dy);
            const auto range=Detail::tensorBernsteinBounds(controls,
                {2*lo.x-px,2*lo.y-py},{2*hi.x-px,2*hi.y-py});
            result.lo=std::min(result.lo,range.lo);result.hi=std::max(result.hi,range.hi);
        }
        return result;
    }

    // Equivalent native-chart Q4 columns of ONE named piece extension.
    // These rounded terminal coefficients never become the scalar/root owner.
    // The mapped graph must prove membership in that same actual piece.
    std::array<double,25> pieceNativeLocalDofs(int piece) const
    {
        checkedPiece(piece);const int px=piece%2,py=piece/2;
        constexpr int orders[7]={0,0,0,1,1,2,2};
        std::array<std::array<Exact,3>,2> scales;
        for(int axis=0;axis<2;++axis){scales[axis][0]=number(1);scales[axis][1]=width_[axis];
            scales[axis][2]=ExactExpansion::multiply(width_[axis],width_[axis]);}
        std::array<Exact,49> native;
        for(int y=0;y<7;++y)for(int x=0;x<7;++x)
            native[7*y+x]=ExactExpansion::multiply(number(physical_[7*y+x]),
                ExactExpansion::multiply(scales[0][orders[x]],scales[1][orders[y]]));
        std::array<double,25> result;
        for(int y=0;y<5;++y)for(int x=0;x<5;++x){Exact sum;
            for(int j=0;j<7;++j)for(int i=0;i<7;++i){
                const double factor=static_cast<double>(nativePieceDatumNumerator[px][x][i])*nativePieceDatumNumerator[py][y][j];
                if(factor)ExactExpansion::accumulate(sum,ExactExpansion::multiply(native[7*j+i],number(factor)));}
            result[5*y+x]=ratio(std::move(sum),number(768.0*768.0));}
        return result;
    }
    std::array<long double,49> pieceNativeDataWeights(int piece,const std::array<long double,25>& weights) const
    {
        checkedPiece(piece);const int px=piece%2,py=piece/2;
        constexpr int inputOrder[7]={0,0,0,1,1,2,2},outputOrder[5]={0,0,1,1,0};
        const long double hx=ExactExpansion::estimate(width_[0]),hy=ExactExpansion::estimate(width_[1]);
        std::array<long double,49> result{};
        for(long double weight:weights)Detail::require(std::isfinite(weight),"nonfinite C2 native-piece sensitivity");
        for(int y=0;y<7;++y)for(int x=0;x<7;++x)for(int j=0;j<5;++j)for(int i=0;i<5;++i){
            if(weights[5*j+i]==0)continue;
            const long double factor=static_cast<long double>(nativePieceDatumNumerator[px][i][x])
                *nativePieceDatumNumerator[py][j][y]/(768.0L*768.0L);
            const long double term=weights[5*j+i]*factor*std::pow(hx,inputOrder[x]-outputOrder[i])
                *std::pow(hy,inputOrder[y]-outputOrder[j]);
            const auto sum=result[7*y+x]+term;
            Detail::require(std::isfinite(term)&&std::isfinite(sum),"unrepresentable C2 native-piece transpose");
            result[7*y+x]=sum;
        }
        return result;
    }
    int pieceIndexForBox(Vector2D lower,Vector2D upper) const
    {
        Detail::require(lower.x>=0&&lower.y>=0&&upper.x<=1&&upper.y<=1
            &&lower.x<=upper.x&&lower.y<=upper.y,"invalid C2 piece native box");
        const auto half=[](double lo,double hi){
            Detail::require(hi<=.5||lo>=.5,"C2 graph box crosses a scalar half-line");
            return lo>=.5?1:0;
        };
        return half(lower.x,upper.x)+2*half(lower.y,upper.y);
    }

    // Internal exact piece data are homogeneous LOCAL QuarticCell DOFs,
    // not the old native25 physical columns. Later restriction adapters must
    // retain piece provenance and the positive denominator together.
    const QuarticCell::ExpandedControls& pieceLocalNumerators(int piece) const
    { return checkedPiece(piece).local; }
    const Exact& pieceLocalDenominator(int piece) const
    { return checkedPiece(piece).denominator; }

  private:
    struct Piece {
        QuarticCell::ExpandedControls local;
        Exact denominator;
        std::array<Interval,25> bounds;
        mutable std::array<std::once_flag,2> gradientBoundsReady;
        mutable std::array<std::array<Interval,25>,2> gradientBounds;
    };
    struct PiecePoint {
        int piece=0;
        RationalSegmentCoordinate x,y;
    };
    Dofs physical_;
    Vector2D lower_,upper_;
    std::array<Exact,2> width_;
    std::array<Piece,4> pieces_;
    long double controlLower_=0,controlUpper_=0,resolutionScale_=0;

    // Fixed seven-basis evaluated as piece [V0,V1,D10,D11,Vmid].
    // These integers are the prepared power basis evaluated at 0,1,1/2
    // and differentiated at 0,1; denominator 768 in EACH native axis.
    inline static constexpr int pieceDatumNumerator[2][5][7]={
        {{768,0,0,0,0,0,0},{0,768,0,0,0,0,0},{0,0,0,384,0,0,0},
         {-768,0,768,-192,-192,-16,16},{576,240,-48,108,12,7,-1}},
        {{0,768,0,0,0,0,0},{0,0,768,0,0,0,0},
         {-768,0,768,-192,-192,-16,16},{0,0,0,0,384,0,0},
         {-48,240,576,-12,-108,-1,7}}
    };
    std::array<Interval,25> physicalDerivativeControls(const Piece& piece,unsigned dx,unsigned dy) const
    {
        if(dx==0&&dy==0)return piece.bounds;
        if(dx+dy==1){const int axis=dy?1:0;
            std::call_once(piece.gradientBoundsReady[axis],[&]{piece.gradientBounds[axis]=buildPhysicalDerivativeControls(piece,dx,dy);});
            return piece.gradientBounds[axis];}
        return buildPhysicalDerivativeControls(piece,dx,dy);
    }
    std::array<Interval,25> buildPhysicalDerivativeControls(const Piece& piece,unsigned dx,unsigned dy) const
    {
        auto controls=QuarticCell::bernsteinControlsCleared(piece.local);
        auto denominator=ExactExpansion::multiply(piece.denominator,number(144));
        int nx=4,ny=4;
        for(unsigned order=0;order<dx;++order){
            std::array<Exact,25> next;
            for(int y=0;y<=ny;++y)for(int x=0;x<nx;++x)
                next[5*y+x]=ExactExpansion::multiply(number(nx),ExactExpansion::subtract(controls[5*y+x+1],controls[5*y+x]));
            controls=std::move(next);--nx;
        }
        for(unsigned order=0;order<dy;++order){
            std::array<Exact,25> next;
            for(int y=0;y<ny;++y)for(int x=0;x<=nx;++x)
                next[5*y+x]=ExactExpansion::multiply(number(ny),ExactExpansion::subtract(controls[5*(y+1)+x],controls[5*y+x]));
            controls=std::move(next);--ny;
        }
        while(nx<4){std::array<Exact,25> next;
            for(int y=0;y<=ny;++y)for(int x=0;x<=nx+1;++x){
                if(x<=nx)ExactExpansion::accumulate(next[5*y+x],ExactExpansion::multiply(number(nx+1-x),controls[5*y+x]));
                if(x>0)ExactExpansion::accumulate(next[5*y+x],ExactExpansion::multiply(number(x),controls[5*y+x-1]));}
            denominator=ExactExpansion::multiply(denominator,number(nx+1));controls=std::move(next);++nx;}
        while(ny<4){std::array<Exact,25> next;
            for(int y=0;y<=ny+1;++y)for(int x=0;x<=nx;++x){
                if(y<=ny)ExactExpansion::accumulate(next[5*y+x],ExactExpansion::multiply(number(ny+1-y),controls[5*y+x]));
                if(y>0)ExactExpansion::accumulate(next[5*y+x],ExactExpansion::multiply(number(y),controls[5*(y-1)+x]));}
            denominator=ExactExpansion::multiply(denominator,number(ny+1));controls=std::move(next);++ny;}
        for(auto& control:controls)control=ExactExpansion::scalePow2(control,dx+dy);
        for(unsigned k=0;k<dx;++k)denominator=ExactExpansion::multiply(denominator,width_[0]);
        for(unsigned k=0;k<dy;++k)denominator=ExactExpansion::multiply(denominator,width_[1]);
        ExactExpansion::centerHomogeneous(controls,denominator);
        const auto divisor=ExactExpansion::enclosure(denominator);
        Detail::require(divisor.lo>0&&std::isfinite(divisor.hi),"unresolved C2 macro derivative-bound scale");
        std::array<Interval,25> result;
        for(int k=0;k<25;++k){const auto value=ExactExpansion::enclosure(controls[k]);
            result[k]={Detail::Outward::down(value.lo/(value.lo<0?divisor.lo:divisor.hi)),
                Detail::Outward::up(value.hi/(value.hi<0?divisor.hi:divisor.lo))};}
        return result;
    }
    inline static constexpr int nativePieceDatumNumerator[2][5][7]={
        {{768,0,0,0,0,0,0},{6912,-12288,6144,2304,-1536,256,128},{0,0,0,768,0,0,0},{43008,-73728,30720,14592,-7680,1664,640},{0,768,0,0,0,0,0}},
        {{6144,-12288,6912,1536,-2304,128,256},{0,0,768,0,0,0,0},{-30720,73728,-43008,-7680,14592,-640,-1664},{0,0,0,0,768,0,0},{0,768,0,0,0,0,0}}
    };
    static Exact number(double value) { return ExactExpansion::fromDouble(value); }
    static Exact difference(double first,double second)
    { return ExactExpansion::subtract(number(first),number(second)); }
    const Piece& checkedPiece(int piece) const
    {
        Detail::require(piece>=0&&piece<4,"C2 macro piece index must be 0..3");
        return pieces_[piece];
    }
    static RationalSegmentCoordinate chart(CellChartCoordinate point)
    {
        const auto one=number(1),tail=number(point.offset);
        const auto value=point.fromUpper?ExactExpansion::subtract(one,tail):tail;
        return {value,value,one};
    }
    static std::pair<int,RationalSegmentCoordinate> selectAxis(
        const RationalSegmentCoordinate& native)
    {
        Detail::require(ExactExpansion::sign(native.denominator)>0
            &&ExactExpansion::sign(native.numeratorBegin)>=0
            &&ExactExpansion::sign(ExactExpansion::subtract(native.denominator,
                native.numeratorBegin))>=0,
            "C2 macro coordinate leaves its actual native chart");
        const auto twice=ExactExpansion::scalePow2(native.numeratorBegin,1);
        const int half=ExactExpansion::sign(
            ExactExpansion::subtract(twice,native.denominator))>0?1:0;
        const auto local=half?ExactExpansion::subtract(twice,native.denominator):twice;
        return {half,{local,local,native.denominator}};
    }
    static PiecePoint selectRational(const RationalSegmentCoordinate& x,
        const RationalSegmentCoordinate& y)
    {
        const auto px=selectAxis(x),py=selectAxis(y);
        return {2*py.first+px.first,px.second,py.second};
    }
    static PiecePoint selectNative(CellChartPoint point)
    {
        ExactExpansion::validateEnvironment();
        validateCellChart(point);
        return selectRational(chart(point.x),chart(point.y));
    }
    static PiecePoint explicitPiece(int piece,CellChartPoint point)
    {
        ExactExpansion::validateEnvironment();
        Detail::require(piece>=0&&piece<4,"C2 macro piece index must be 0..3");
        validateCellChart(point);
        return {piece,chart(point.x),chart(point.y)};
    }
    PiecePoint selectPhysical(Vector2D point) const
    {
        ExactExpansion::validateEnvironment();
        Detail::require(std::isfinite(point.x)&&std::isfinite(point.y)
            &&point.x>=lower_.x&&point.x<=upper_.x
            &&point.y>=lower_.y&&point.y<=upper_.y,
            "C2 macro physical point leaves its actual native cell");
        const auto nx=difference(point.x,lower_.x),ny=difference(point.y,lower_.y);
        return selectRational({nx,nx,width_[0]},{ny,ny,width_[1]});
    }
    static double ratio(Exact numerator,Exact denominator)
    {
        Detail::require(ExactExpansion::sign(denominator)>0,
            "nonpositive C2 macro homogeneous denominator");
        std::array<Exact,1> values{{std::move(numerator)}};
        ExactExpansion::centerHomogeneous(values,denominator);
        const long double value=static_cast<long double>(ExactExpansion::estimate(values[0]))
            /ExactExpansion::estimate(denominator);
        const double result=Detail::checked(value,"unrepresentable C2 macro data ratio");
        Detail::require(ExactExpansion::sign(values[0])==0||result!=0,
            "C2 macro data ratio underflow");
        return result;
    }
    double physicalRatio(Exact numerator,Exact denominator,unsigned dx,unsigned dy) const
    {
        numerator=ExactExpansion::scalePow2(numerator,static_cast<int>(dx+dy));
        for(unsigned k=0;k<dx;++k)denominator=ExactExpansion::multiply(denominator,width_[0]);
        for(unsigned k=0;k<dy;++k)denominator=ExactExpansion::multiply(denominator,width_[1]);
        return ratio(std::move(numerator),std::move(denominator));
    }
    Sample jet(const PiecePoint& point) const
    {
        const auto& piece=checkedPiece(point.piece);
        const auto bx=QuarticCell::rationalBasisJet(
            point.x.numeratorBegin,point.x.denominator,true);
        const auto by=QuarticCell::rationalBasisJet(
            point.y.numeratorBegin,point.y.denominator,true);
        const auto basis=[](const QuarticCell::RationalBasisJet& value,int order)
            ->const std::array<Exact,5>& {
            return order==0?value.value:order==1?value.derivative:value.second;
        };
        std::array<std::array<Exact,5>,3> alongY;
        for(int order=0;order<3;++order)for(int x=0;x<5;++x)
            for(int y=0;y<5;++y)if(ExactExpansion::sign(piece.local[5*y+x])!=0
                &&ExactExpansion::sign(basis(by,order)[y])!=0)
                ExactExpansion::accumulate(alongY[order][x],
                    ExactExpansion::multiply(piece.local[5*y+x],basis(by,order)[y]));
        const auto denominator=ExactExpansion::multiply(piece.denominator,
            ExactExpansion::multiply(bx.denominator,by.denominator));
        constexpr int xOrder[6]={0,1,0,2,1,0},yOrder[6]={0,0,1,0,1,2};
        std::array<double,6> values;
        for(int channel=0;channel<6;++channel) {
            Exact numerator;
            for(int x=0;x<5;++x)if(ExactExpansion::sign(alongY[yOrder[channel]][x])!=0
                &&ExactExpansion::sign(basis(bx,xOrder[channel])[x])!=0)
                ExactExpansion::accumulate(numerator,ExactExpansion::multiply(
                    alongY[yOrder[channel]][x],basis(bx,xOrder[channel])[x]));
            values[channel]=physicalRatio(std::move(numerator),denominator,
                xOrder[channel],yOrder[channel]);
        }
        return {values[0],values[1],values[2],values[3],values[4],values[5]};
    }
    double derivative(const PiecePoint& point,unsigned dx,unsigned dy) const
    {
        const auto& piece=checkedPiece(point.piece);
        auto result=QuarticCell::homogeneousRestrictionCentered(
            piece.local,piece.denominator,point.x,point.y,dx,dy);
        return physicalRatio(std::move(result.numerator[0]),
            std::move(result.denominator),dx,dy);
    }
};
} // namespace LevelSetHermiteScalar
} // namespace
} // namespace LAMP::Geometry::Detail
