#include "Modules/Interface/Compatibility/InterfaceTransportInput.h"
#include "Core/Mesh/NativeLeafDomain.h"

#include <stdexcept>
#include <utility>
#include <cmath>
#include <sstream>
#include <iomanip>
#include <limits>
namespace {
Quantum::Geometry::UnitPolynomialRoots checkedRoots(Quantum::Geometry::UnitPolynomialRoots roots)
{
    if (roots.identicallyZero && !roots.roots.empty()) throw std::invalid_argument("zero flux trace has isolated roots");
    double previous=-1;
    for (const auto& root:roots.roots) {
        if (!std::isfinite(root.point) || !std::isfinite(root.lower) || !std::isfinite(root.upper)
            || root.lower<0 || root.upper>1 || root.lower>root.point || root.point>root.upper || root.lower<=previous)
            throw std::invalid_argument("invalid, overlapping or unordered retained flux roots");
        previous=root.upper;
    }
    return roots;
}
}

InterfaceTransportDomain::InterfaceTransportDomain(const NativeLeafDomain& geometry, BoundaryQuery boundary)
    : geometry_(&geometry), boundary_(std::move(boundary))
{
    if (!boundary_) throw std::invalid_argument("Interface transport requires an explicit boundary query");
}

InterfaceTransportDomain::Boundary InterfaceTransportDomain::boundary(int nativeFace) const
{
    if (nativeFace < 0 || nativeFace >= geometry_->nativeFaceCount())
        throw std::out_of_range("Interface transport boundary native face is out of range");
    const auto result = boundary_(nativeFace);
    switch (result)
    {
    case Boundary::Interior: case Boundary::StationaryWall:
    case Boundary::Symmetry: case Boundary::SinglePhasePort: case Boundary::PrescribedScalar: return result;
    }
    throw std::invalid_argument("Unsupported interface transport boundary constraint");
}

InterfaceNormalFluxView::InterfaceNormalFluxView(Queries queries) : queries_(std::move(queries))
{
    if (!queries_.chart || !queries_.integral || !queries_.phaseIntegral || !queries_.density)
        throw std::invalid_argument("Interface normal flux requires one complete frozen primitive");
    if((queries_.retainedMomentDegree!=0 && queries_.retainedMomentDegree!=2)
        || ((queries_.retainedMomentDegree || queries_.retainedVelocityProducts) && !queries_.polynomial))
        throw std::invalid_argument("retained face moments require degree two and the actual polynomial trace");
}
InterfaceNormalFluxView InterfaceNormalFluxView::scaled(std::shared_ptr<const InterfaceNormalFluxView> spatial,double amplitude)
{
    if (!spatial || spatial->spatial_ || !spatial->hasNormalFluxRoots() || !std::isfinite(amplitude))
        throw std::invalid_argument("separable flux requires one unscaled retained spatial primitive with roots and finite amplitude");
    Queries q;q.retainedMomentDegree=spatial->retainedMomentDegree();q.retainedVelocityProducts=spatial->retainedVelocityProducts();
    q.chart=[spatial](int r){return spatial->faceChart(r);};
    q.integral=[spatial,amplitude](int r,int half,double a,double b,double s,double t) {
        return amplitude*spatial->integratedNormalFlux(r,half,a,b,s,t);
    };
    q.phaseIntegral=[spatial,amplitude](int r,int half,const LevelSetMappedSegment& segment) {
        return amplitude*spatial->integratedNormalFlux(r,half,segment);
    };
    q.density=[spatial,amplitude](int r,int half,bool end,long double s) {
        return amplitude*spatial->normalFluxDensity(r,half,end,s);
    };
    q.roots=[spatial,amplitude](int r,int half) {
        if (amplitude==0) return Quantum::Geometry::UnitPolynomialRoots{{},true};
        return spatial->normalFluxRoots(r,half);
    };
    if(spatial->hasNormalFluxPolynomial())q.polynomial=[spatial,amplitude](int r,int half) {
        auto terms=spatial->normalFluxPolynomial(r,half);for(auto& term:terms)term.weight*=amplitude;return terms;
    };
    if(spatial->hasRangeIntegral())q.rangeIntegral=[spatial,amplitude](int r,int half,double a,double b,bool end,long double offset,long double length) {
        return amplitude*spatial->integratedRange(r,half,a,b,end,offset,length);
    };
    InterfaceNormalFluxView result(std::move(q));result.spatial_=std::move(spatial);result.amplitude_=amplitude;return result;
}
long double InterfaceNormalFluxView::densityMagnitudeBound(int row,int half) const
{
    long double bound=0;
    if(spatial_)bound=std::abs(amplitude_)*spatial_->densityMagnitudeBound(row,half);
    else if(queries_.densityMagnitudeBound)bound=queries_.densityMagnitudeBound(row,half);
    else if(queries_.polynomial) {
        // |sum a_k s^k| <= sum |a_k| on [0,1], with all retained factors.
        // Deliberately an upper bound, not a sampled maximum.
        for(const auto& term:queries_.polynomial(row,half))for(double a:term.coefficients)
            bound+=std::abs(static_cast<long double>(a)*term.factor*term.weight);
        bound*=1+64*std::numeric_limits<long double>::epsilon();
    } else throw std::invalid_argument("normal-flux trace has no analytic magnitude bound");
    if(!(bound>=0 && std::isfinite(bound)))throw std::invalid_argument("invalid normal-flux magnitude bound");
    return bound;
}
bool InterfaceNormalFluxView::canInterpolateRoots(const InterfaceNormalFluxView& next) const noexcept
{
    return (hasNormalFluxPolynomial() && next.hasNormalFluxPolynomial()) || (spatial_ && spatial_==next.spatial_);
}
Quantum::Geometry::UnitPolynomialRoots InterfaceNormalFluxView::interpolatedRoots(
    const InterfaceNormalFluxView& next,int row,int half,double fraction) const
{
    if (!canInterpolateRoots(next) || !std::isfinite(fraction) || fraction<0 || fraction>1)
        throw std::invalid_argument("flux layers have no common interpolation/root authority");
    if (spatial_ && spatial_==next.spatial_) {
        const auto amplitude=(1-static_cast<long double>(fraction))*amplitude_+static_cast<long double>(fraction)*next.amplitude_;
        if (amplitude==0) return {{},true};
        return checkedRoots(spatial_->normalFluxRoots(row,half));
    }
    auto terms=normalFluxPolynomial(row,half);for (auto& term:terms) term.weight*=1-fraction;
    for (auto term:next.normalFluxPolynomial(row,half)) {term.weight*=fraction;terms.push_back(term);}
    try {return checkedRoots(Quantum::Geometry::isolateUnitPolynomial(terms));}
    catch(const std::exception& error) {
        std::ostringstream message;message<<std::setprecision(21)<<error.what()
            <<"; interpolated flux row="<<row<<" half="<<half<<" fraction="<<fraction;
        for(const auto& term:terms) {
            message<<"; factor="<<term.factor<<" weight="<<term.weight<<" coefficients=";
            for(double coefficient:term.coefficients)message<<coefficient<<',';
        }
        throw std::runtime_error(message.str());
    }
}
long double InterfaceNormalFluxView::interpolatedIntegral(const InterfaceNormalFluxView& next,int row,int half,
    double fraction,double lo,double hi,double partLo,double partHi) const
{
    if (!canInterpolateRoots(next) || !std::isfinite(fraction) || fraction<0 || fraction>1)
        throw std::invalid_argument("flux layers have no common interpolation authority");
    if (spatial_ && spatial_==next.spatial_) return ((1-static_cast<long double>(fraction))*amplitude_
        +static_cast<long double>(fraction)*next.amplitude_)*spatial_->integratedNormalFlux(row,half,lo,hi,partLo,partHi);
    return static_cast<long double>(1-fraction)*integratedNormalFlux(row,half,lo,hi,partLo,partHi)
        +static_cast<long double>(fraction)*next.integratedNormalFlux(row,half,lo,hi,partLo,partHi);
}
long double InterfaceNormalFluxView::integratedRange(int row,int half,double a0,double a1,bool fromEnd,long double offset,long double length) const
{
    if(queries_.rangeIntegral)return queries_.rangeIntegral(row,half,a0,a1,fromEnd,offset,length);
    if(queries_.polynomial)return Quantum::Geometry::integrateUnitPolynomialRange(queries_.polynomial(row,half),a0,a1,fromEnd,offset,length);
    throw std::logic_error("normal flux has no retained range primitive");
}
long double InterfaceNormalFluxView::interpolatedRange(const InterfaceNormalFluxView& next,int row,int half,double fraction,
    double a0,double a1,bool fromEnd,long double offset,long double length) const
{
    if(!canInterpolateRoots(next) || !std::isfinite(fraction) || fraction<0 || fraction>1)
        throw std::invalid_argument("flux layers have no common range interpolation authority");
    if(spatial_ && spatial_==next.spatial_)return ((1-static_cast<long double>(fraction))*amplitude_+
        static_cast<long double>(fraction)*next.amplitude_)*spatial_->integratedRange(row,half,a0,a1,fromEnd,offset,length);
    return static_cast<long double>(1-fraction)*integratedRange(row,half,a0,a1,fromEnd,offset,length)+
        static_cast<long double>(fraction)*next.integratedRange(row,half,a0,a1,fromEnd,offset,length);
}
Quantum::Geometry::UnitPolynomialRoots InterfaceNormalFluxView::interpolationTimeRoots(const InterfaceNormalFluxView& next) const
{
    if (!canInterpolateRoots(next)) throw std::invalid_argument("flux layers have no common time-root authority");
    if (!spatial_ || spatial_!=next.spatial_) return {};
    Quantum::Geometry::UnitPolynomialTerm old,newer;
    old.coefficients[0]=1;old.coefficients[1]=-1;old.factor=amplitude_;
    newer.coefficients[1]=1;newer.factor=next.amplitude_;
    return Quantum::Geometry::isolateUnitPolynomial({old,newer});
}
