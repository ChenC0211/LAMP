#pragma once

#include "Core/Geometry/Detail/LevelSetMappedSegment.h"
#include "Core/Geometry/Detail/LevelSetMappedCurveCell.h"
#include "Core/Geometry/Detail/LevelSetMappedRootRefinement.h"
#include "Core/Geometry/LevelSetMappedInterface.h"
#include <initializer_list>
#include <sstream>
#include <iomanip>
#include <map>
#include <functional>
#include <numeric>
#include <type_traits>
#include <optional>
#include <memory>
#include <chrono>

namespace LAMP::Geometry::Detail
{
namespace
{
namespace MappedHermiteInterface
{
namespace M = MappedHermiteSegment;
namespace Exact = ExactExpansion;
using Ratio = M::Ratio;
using Expansion = Exact::Expansion;
using Point = std::array<Ratio,2>;
using Public = ::LevelSetMappedInterface;
// Query-local passive clock. Null profiles never read the clock. Nested
// restrictions are subtracted from root engine/refinement spans; the total
// scope assigns the unclassified remainder after all child scopes commit.
class NativeWorkTimer
{
    using Profile=Public::NativeMomentWorkProfile;
    using Field=double Profile::*;
    Profile* profile_;Field field_,excluded_;bool total_;
    double before_=0;
    std::chrono::steady_clock::time_point started_{};
    static double classified(const Profile& p) {
        return p.restrictionConstruction+p.initialRootCacheLookup+p.initialRootCacheRemember
            +p.initialRootIsolation+p.rootMappingRefinement+p.fixedAndSetup;
    }
public:
    NativeWorkTimer(Profile* profile,Field field,Field excluded=nullptr,bool total=false)
        :profile_(profile),field_(field),excluded_(excluded),total_(total)
    {
        if(profile_) {
            before_=total_?classified(*profile_):excluded_?profile_->*excluded_:0.;
            started_=std::chrono::steady_clock::now();
        }
    }
    void stop() {
        if(!profile_)return;
        const double elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-started_).count();
        if(total_) {
            profile_->total+=elapsed;
            profile_->remainingIntegration+=elapsed-(classified(*profile_)-before_);
        } else profile_->*field_+=elapsed-(excluded_?profile_->*excluded_-before_:0.);
        profile_=nullptr;
    }
    ~NativeWorkTimer(){stop();}
    NativeWorkTimer(const NativeWorkTimer&)=delete;
    NativeWorkTimer& operator=(const NativeWorkTimer&)=delete;
};

using EdgeIdentity=std::pair<int,int>; // Canonical side and complete-root ordinal.
struct LocalRootResolutionError : std::runtime_error
{
    std::vector<EdgeIdentity> roots;
    LocalRootResolutionError(std::string message,std::vector<EdgeIdentity> requested)
        : std::runtime_error(std::move(message)),roots(std::move(requested)) {}
};

inline void require(bool condition, const char* message)
{
    if (!condition) throw std::runtime_error(std::string("[LevelSetMappedInterface] ")+message);
}
// Failure-only formatting preserves arithmetic and operation order. Retained
// expansion components, rather than rounded points, identify a failing scale.
inline std::string describe(const Ratio& value)
{
    std::ostringstream out;
    out << std::setprecision(17) << "N={";
    for (const auto component:value.numerator.components()) out << component << ',';
    out << "},D={";
    for (const auto component:value.denominator.components()) out << component << ',';
    return out.str()+"}";
}
inline std::string describe(const Point& value)
{ return "("+describe(value[0])+";"+describe(value[1])+")"; }
template<class Action,class Description>
inline auto contextual(Action&& action,Description&& description) -> decltype(action())
{
    try { return action(); }
    catch (const Public::RootResolutionError&) { throw; }
    catch (const LocalRootResolutionError&) { throw; }
    catch (const std::exception& error)
    { throw std::runtime_error(description()+": "+error.what()); }
}
template<class Action>
inline auto stage(const char* name,Action&& action) -> decltype(action())
{ return contextual(std::forward<Action>(action),[&] { return std::string("[MappedInterface ")+name+"]"; }); }
inline Ratio normalized(Ratio value)
{
    if (Exact::sign(value.numerator) == 0) return M::constant(0.0);
    int exponent = 0;
    std::frexp(Exact::estimate(value.denominator), &exponent);
    value.numerator = Exact::scalePow2(value.numerator,-exponent);
    value.denominator = Exact::scalePow2(value.denominator,-exponent);
    return value;
}
inline Ratio sum(const Ratio& a, const Ratio& b)
{
    if (Exact::sign(a.numerator) == 0) return b;
    if (Exact::sign(b.numerator) == 0) return a;
    if (Exact::sign(Exact::subtract(a.denominator,b.denominator)) == 0)
        return normalized(M::ratio(Exact::add(a.numerator,b.numerator),a.denominator));
    return normalized(M::add(a,b));
}
inline Ratio difference(const Ratio& a, const Ratio& b)
{
    if (Exact::sign(b.numerator) == 0) return a;
    if (Exact::sign(Exact::subtract(a.denominator,b.denominator)) == 0)
        return normalized(M::ratio(Exact::subtract(a.numerator,b.numerator),a.denominator));
    return normalized(M::subtract(a,b));
}
inline Ratio product(const Ratio& a, const Ratio& b)
{
    if (Exact::sign(a.numerator) == 0 || Exact::sign(b.numerator) == 0) return M::constant(0.0);
    return normalized(M::ratio(Exact::multiply(a.numerator,b.numerator),
        Exact::multiply(a.denominator,b.denominator)));
}
inline Ratio quotient(const Ratio& a, const Ratio& b)
{
    require(Exact::sign(b.numerator) != 0, "zero mapped coordinate divisor");
    return normalized(M::ratio(Exact::multiply(a.numerator,b.denominator),
        Exact::multiply(a.denominator,b.numerator)));
}
inline Ratio between(const Ratio& a, const Ratio& b, double t)
{
    if (t == 0.0) return a;
    if (t == 1.0) return b;
    return sum(a,product(difference(b,a),M::constant(t)));
}
inline bool same(const Point& a, const Point& b)
{ return M::compare(a[0],b[0]) == 0 && M::compare(a[1],b[1]) == 0; }
inline Point point(Vector2D p) { return {{M::constant(p.x),M::constant(p.y)}}; }
inline Expansion times(const Expansion& value, double factor)
{ return Exact::multiply(value,M::number(factor)); }
inline int orientation(Vector2D a, Vector2D b, const Point& p)
{
    const auto x = difference(p[0],M::constant(a.x)), y = difference(p[1],M::constant(a.y));
    return Exact::sign(difference(product(M::ratio(M::difference(b.x,a.x),M::number(1.0)),y),
        product(M::ratio(M::difference(b.y,a.y),M::number(1.0)),x)).numerator);
}
// Query-private coefficients for the SAME transverse Q4. The stamp is in
// the original cell chart; any mapping or exact root refinement invalidates
// the packet rather than reusing values from the previous representative.
struct FilteredGeometryPacket
{
    std::array<Exact::Interval,5> controls{};
    double parameter=0.;
};
struct Root
{
    Ratio parameter, lower, upper;
    // Sign of dP/dt in the original increasing local coordinate. Nonzero
    // only for certified simple INTERIOR roots; endpoint multiplicity is
    // intentionally not classified by this field.
    int derivativeSign=0;
    // Original scalar signs just before/after this root on its edge. Zero
    // means outside the requested segment, not a zero-length phase interval.
    int phaseBefore=0,phaseAfter=0;
    bool deflated=false;
    std::optional<FilteredGeometryPacket> geometryPacket;
};
// Geometry-only numerical certificate for one ordinary transverse root.
// Its interval encloses the root of the SAME Q4, not a fitted or rounded
// surrogate. Lack of a certificate explicitly selects the exact reference.
struct FilteredTransverseRoot
{
    double parameter=0.,lower=0.,upper=0.;
    int derivativeSign=0;
};
inline Exact::Interval filteredIntervalDivide(Exact::Interval n,Exact::Interval d)
{
    using A=LevelSetHermiteScalar::Detail::Outward;
    if(!A::finite(n)||!A::finite(d)||(d.lo<=0.&&d.hi>=0.))return A::whole();
    if(A::zero(n))return {};
    const std::array<double,4> values{{n.lo/d.lo,n.lo/d.hi,n.hi/d.lo,n.hi/d.hi}};
    for(double value:values)if(!std::isfinite(value))return A::whole();
    return {A::down(*std::min_element(values.begin(),values.end())),
        A::up(*std::max_element(values.begin(),values.end()))};
}
inline bool filteredTransverseRoot(const std::array<Exact::Interval,5>& controls,
    double lower,double upper,FilteredTransverseRoot& output)
{
    using A=LevelSetHermiteScalar::Detail::Outward;
    using Interval=Exact::Interval;
    Exact::validateEnvironment();
    require(std::isfinite(lower)&&std::isfinite(upper)&&lower>=0.&&upper<=1.&&lower<upper,
        "invalid filtered transverse chart interval");
    for(const auto control:controls)if(!A::finite(control))return false;
    const auto evaluate=[&](double t) {
        auto work=controls;
        if(t==0.)return work.front();
        if(t==1.)return work.back();
        const auto argument=A::point(t),complement=A::add(A::point(1.),A::negate(argument));
        for(int degree=4;degree>0;--degree)for(int k=0;k<degree;++k)
            work[k]=A::add(A::multiply(complement,work[k]),A::multiply(argument,work[k+1]));
        return work[0];
    };
    const auto sign=[](Interval value) {return value.lo>0.?1:value.hi<0.?-1:0;};
    Interval derivative{std::numeric_limits<double>::infinity(),-std::numeric_limits<double>::infinity()};
    for(int k=0;k<4;++k) {
        const auto slope=A::multiply(A::point(4.),A::add(controls[k+1],A::negate(controls[k])));
        if(!A::finite(slope))return false;
        derivative.lo=std::min(derivative.lo,slope.lo);
        derivative.hi=std::max(derivative.hi,slope.hi);
    }
    // Bernstein convex hull certifies a nonzero derivative throughout the
    // entire transverse interval. Sampling equal signs is never a no-root
    // proof, and the graph's orientation flag is not a numerical slope bound.
    const int direction=sign(derivative);
    if(direction==0)return false;
    const auto first=evaluate(lower),last=evaluate(upper);
    if(!A::finite(first)||!A::finite(last)||sign(first)!=-direction||sign(last)!=direction)return false;
    // Strict endpoint signs plus strict monotonicity prove one interior root.
    // The mean-value theorem contracts its retained bracket using only
    // outward arithmetic; a residual estimate never authorizes acceptance.
    for(int iteration=0;iteration<24;++iteration) {
        const double middle=.5*lower+.5*upper;
        if(!(middle>lower&&middle<upper))return false;
        const double width=upper-lower;
        const double distance=std::min({middle,1.-middle,std::abs(middle-.5)});
        if(lower>0.&&upper<1.&&!(lower<=.5&&upper>=.5)
            &&distance>0.&&width<=1e-13*distance) {
            output={middle,lower,upper,direction};
            Exact::validateEnvironment();
            return true;
        }
        const auto value=evaluate(middle);
        const auto contracted=A::add(A::point(middle),A::negate(filteredIntervalDivide(value,derivative)));
        if(!A::finite(contracted))return false;
        const double nextLower=std::max(lower,contracted.lo),nextUpper=std::min(upper,contracted.hi);
        // A contradictory or stagnating filter is unqualified. The original
        // exact path diagnoses difficult input; no approximate root is used.
        if(!(nextLower<nextUpper)||nextUpper-nextLower>=width)return false;
        lower=nextLower;upper=nextUpper;
    }
    return false;
}
// Point estimates and their outward error bounds, in scalar units and
// dimensionless native-cell derivative units. Root-enclosure width and GL8
// quadrature error are separate and are never added to these point bounds.
struct FilteredTransverseJet
{
    long double value=0.,valueError=0.,derivative=0.,derivativeError=0.;
};
inline bool filteredTransverseGeometryJet(const std::array<Exact::Interval,5>& controls,
    double parameter,int derivativeSign,FilteredTransverseJet& output)
{
    using A=LevelSetHermiteScalar::Detail::Outward;
    Exact::validateEnvironment();
    if(!std::isfinite(parameter)||parameter<0.||parameter>1.
        ||(derivativeSign!=1&&derivativeSign!=-1))return false;
    auto work=controls;
    for(const auto control:work)if(!A::finite(control))return false;
    const auto argument=A::point(parameter),complement=A::add(A::point(1.),A::negate(argument));
    // After three de Casteljau levels, both the quartic value and its
    // derivative use the same last two values: P'=4*(work[1]-work[0]).
    // The outward operations enclose the actual source coefficients and q;
    // no fitted polynomial or rounded scalar DOF enters this evaluation.
    for(int degree=4;degree>1;--degree)for(int k=0;k<degree;++k)
        work[k]=A::add(A::multiply(complement,work[k]),A::multiply(argument,work[k+1]));
    const auto derivative=A::multiply(A::point(4.),A::add(work[1],A::negate(work[0])));
    const auto value=A::add(A::multiply(complement,work[0]),A::multiply(argument,work[1]));
    if(!A::finite(value)||!A::finite(derivative)
        ||(derivativeSign>0?derivative.lo<=0.:derivative.hi>=0.))return false;
    const auto estimate=[](Exact::Interval bound,long double& center,long double& error) {
        center=.5L*bound.lo+.5L*bound.hi;
        const long double radius=std::max(center-static_cast<long double>(bound.lo),
            static_cast<long double>(bound.hi)-center);
        error=radius==0.?0.:A::up(static_cast<double>(radius));
        return std::isfinite(center)&&std::isfinite(error)&&error>=0.;
    };
    FilteredTransverseJet sample;
    if(!estimate(value,sample.value,sample.valueError)
        ||!estimate(derivative,sample.derivative,sample.derivativeError)
        ||!(std::abs(sample.derivative)>sample.derivativeError))return false;
    Exact::validateEnvironment();
    output=sample;
    return true;
}
struct Cut
{
    Point local;
    Ratio lower, upper; // Chosen coordinate enclosure.
    std::vector<Public::RootResolutionRequest> sources;
    int edgeSide=-1,edgeOrdinal=-1;
};
struct EdgeBinding
{
    Root root; // Tangential canonical-edge coordinate, independent of chosen_.
    std::vector<Public::RootResolutionRequest> sources;
};
inline long double finiteProduct(std::initializer_list<long double> factors,
    long double divisor=1.0L);

// Terminal floating quadrature basis. Both endpoint distances are converted
// independently from the retained Ratio; never compute 1-rounded(q). This
// preserves the opposite-end Hermite columns even when rounded(q)==1.
struct Basis
{
    std::array<long double,4> value, derivative;
    explicit Basis(const Ratio& coordinate)
    {
        const auto& n=coordinate.numerator; const auto& d=coordinate.denominator;
        const auto complement=Exact::subtract(d,n);
        const long double q=M::approximate(coordinate);
        const long double m=M::approximate(M::ratio(complement,d));
        require(q>=0.0L && q<=1.0L && m>=0.0L && m<=1.0L,
            "terminal Hermite basis coordinate left the unit interval");
        // The potentially cancelling linear derivative factors are also
        // formed before terminal conversion, without raising N or D powers.
        const long double mMinus2q=M::approximate(M::ratio(Exact::subtract(complement,times(n,2.0)),d));
        const long double qMinus2m=M::approximate(M::ratio(Exact::subtract(n,times(complement,2.0)),d));
        value={{finiteProduct({m,m,1.0L+2.0L*q}),finiteProduct({q,q,1.0L+2.0L*m}),
            finiteProduct({q,m,m}),finiteProduct({-q,q,m})}};
        derivative={{finiteProduct({-6.0L,q,m}),finiteProduct({6.0L,q,m}),
            finiteProduct({m,mMinus2q}),finiteProduct({q,qMinus2m})}};
    }
};
struct QuarticBasis
{
    std::array<long double,5> value, derivative;
    explicit QuarticBasis(const Ratio& coordinate)
    {
        const auto basis=LevelSetHermiteScalar::QuarticCell::rationalBasisJet(
            coordinate.numerator,coordinate.denominator);
        for (std::size_t k=0;k<5;++k)
        {
            value[k]=M::approximate(M::ratio(basis.value[k],basis.denominator));
            derivative[k]=M::approximate(M::ratio(basis.derivative[k],basis.denominator));
        }
    }
};
template<std::size_t Count> struct CellJet
{
    // dx,dy below are P_u,P_v in dimensionless CELL coordinates, not
    // physical derivatives. The integral uses the shared-bound hx,hy in
    // ds/|grad(phi)| = hx*hy/|P_transverse| dq. DOF columns instead retain
    // the source derivative scales: legacy nominal hx,hy on uniform grids,
    // actual leaf widths on native grids. Those scale its input derivatives.
    long double value=0.0L, dx=0.0L, dy=0.0L;
    // Numerical diagnostic estimates. The original jet has a floating
    // cancellation screen; a filtered geometry packet supplies a same-Q4
    // point error bound. Topology is supplied by the certified root, never
    // inferred from these diagnostic estimates.
    long double dxRoundoff=0.0L,dyRoundoff=0.0L;
    std::array<long double,Count> columns{};
};
using Jet = CellJet<16>;

// Integrands may have individually very small/large factors. Reconstruct
// once from separated exponents; never silently erase a nonzero product or
// overflow an intermediate whose final product would be representable.
inline long double finiteProduct(std::initializer_list<long double> factors,
    long double divisor)
{
    require(std::isfinite(divisor) && divisor!=0.0L, "invalid integral divisor");
    bool zero=false;
    for (const auto factor:factors)
    {
        require(std::isfinite(factor), "nonfinite integral factor");
        zero=zero || factor==0.0L;
    }
    if (zero) return 0.0L;
    int divisorExponent=0;
    long double mantissa=1.0L/std::frexp(divisor,&divisorExponent);
    int exponent=-divisorExponent;
    for (const auto factor:factors)
    {
        int factorExponent=0;
        mantissa*=std::frexp(factor,&factorExponent);
        exponent+=factorExponent;
    }
    const long double result=std::scalbn(mantissa,exponent);
    require(std::isfinite(result) && result!=0.0L, "unrepresentable nonzero integral product");
    return result;
}

// Layout-independent geometry result: no scalar DOF columns or backend types.
struct GraphSample
{
    CellChartPoint chart;
    long double physicalSpan=0;
    int coordinateDirection=0;
};

template<class CellLayout> class CellIntegrator
{
    static constexpr std::size_t dofCount=CellLayout::dofCount;
    // Exact certificates apply to each retained native Q4 polynomial,
    // including a named piece of the completed C2 source. Cache/filtered
    // acceleration and complete native-side supports have separate scopes.
    static constexpr bool nativeQuarticCertificates=
        std::is_same_v<CellLayout,NativeQuarticCurveCell>
        || std::is_same_v<CellLayout,NativeC2MacroCurveCell>;
    static_assert(dofCount==16 || dofCount==25,"unsupported mapped cell data layout");
    using Jet = CellJet<dofCount>;
    using Basis = std::conditional_t<dofCount==16,MappedHermiteInterface::Basis,QuarticBasis>;
    using Arc = Public::CellArc<dofCount>;
    using Result = Public::IntegrationResult<Arc>;
    bool axisymmetric_ = true;
    bool geometryMeasure_ = false;
    bool guardedTransverseNewton_ = false;
    Public::NativeMomentWorkProfile* work_ = nullptr;
    struct CoincidentSupport { Point begin,end; int normalAxis; };
    const std::vector<CoincidentSupport>* coincidentSupports_=nullptr;
    std::size_t coincidentArcs_=0,localRootRefinements_=0;
    mutable std::size_t transverseRootRefinements_=0;
    Vector2D greenOrigin_{0,0};
    struct BatchOwnership {
        struct Rectangle {Point lower,upper;std::vector<int> owners;};
        std::vector<Rectangle> rectangles;
    };
    // Optional, query-local certificate: every channel has this exact same
    // ordered axis-aligned partition. No rounded coordinates or root state.
    const BatchOwnership* batchOwnership_=nullptr;
    struct GeometryOnly {};
    CellIntegrator(const typename CellLayout::Source& source, const Public::Graph& graph,
        const std::vector<const typename CellLayout::Query*>& cuts, const std::vector<Public::Region>& regions,
        GeometryOnly)
        : layout_(source,graph), graph_(graph), cuts_(cuts), regions_(regions), chosen_(graph.coordinateIsX?0:1)
    {
        Exact::validateEnvironment();
        nominal_=layout_.derivativeScale();
        const auto bounds=layout_.bounds();
        for (int axis=0;axis<2;++axis)
        {
            const double low=axis==0?bounds[0].x:bounds[0].y, high=axis==0?bounds[1].x:bounds[1].y;
            origin_[axis]=M::constant(low);
            width_[axis]=M::ratio(M::difference(high,low),M::number(1.0));
            require(Exact::sign(width_[axis].numerator)>0, "nonpositive shared cell width");
        }
        dofs_=layout_.coefficients();
        require(std::isfinite(graph.localLower.x) && std::isfinite(graph.localLower.y)
            && std::isfinite(graph.localUpper.x) && std::isfinite(graph.localUpper.y)
            && graph.localLower.x>=0 && graph.localLower.y>=0
            && graph.localUpper.x<=1 && graph.localUpper.y<=1
            && graph.localLower.x<graph.localUpper.x && graph.localLower.y<graph.localUpper.y,
            "invalid original-cell graph subbox");
        if(graph.nativeSide)validateNativeSide(source);
    }
public:
    CellIntegrator(const typename CellLayout::Source& source, const Public::Graph& graph,
        const std::vector<const typename CellLayout::Query*>& cuts, const std::vector<Public::Region>& regions,
        bool axisymmetric = true, Vector2D greenOrigin = {0,0})
        : CellIntegrator(source,graph,cuts,regions,GeometryOnly{})
    {
        axisymmetric_=axisymmetric;
        require(std::isfinite(greenOrigin.x) && std::isfinite(greenOrigin.y), "nonfinite Green-potential origin");
        greenOrigin_=greenOrigin;
        require(!regions.empty(), "free arc requires actual ownership regions");
        stage("validate ownership-region boundaries",[&] { validateRegions(); });
    }

    // Geometry-only use cannot obtain an unvalidated integral object. In
    // particular planar graph orientation imposes no RZ radius/volume law.
    static int graphOrientation(const typename CellLayout::Source& source,const Public::Graph& graph)
    {
        const std::vector<const typename CellLayout::Query*> cuts;
        const std::vector<Public::Region> regions;
        return CellIntegrator(source,graph,cuts,regions,GeometryOnly{}).graphDirection();
    }
    static Public::GraphMoments graphMoments(const typename CellLayout::Source& source,
        const Public::Graph& graph,bool axisymmetric)
    {
        require(graph.controlledIntegration,"native moments require a controlled complete graph");
        const std::vector<const typename CellLayout::Query*> cuts;
        const std::vector<Public::Region> regions;
        CellIntegrator integrator(source,graph,cuts,regions,GeometryOnly{});
        integrator.axisymmetric_=axisymmetric;integrator.geometryMeasure_=true;
        if(graph.nativeSide)return integrator.sideMoments(source,axisymmetric);
        Public::GraphMoments moments;
        for (const auto& arc:integrator.run().arcs) {
            moments.radial+=arc.radialMoment;moments.axial+=arc.axialMoment;
            moments.errorEstimate+=arc.momentErrorEstimate;
        }
        // Integrate in the certified graph coordinate without a transverse
        // derivative ratio. Recover the other Green moment by integration
        // by parts on the SAME curve and retained endpoint roots:
        // XY: A=R-sign*[xy] for a y-graph, R=A-sign*[xy] for an x-graph.
        // RZ: A=R-sign*pi*[r^2*z], or R=A-sign*pi*[r^2*z].
        // This is not a planar replacement or a relaxed quadrature gate.
        auto first=integrator.endpoint(graph.begin),last=integrator.endpoint(graph.end);
        if(M::compare(first.local[integrator.chosen_],last.local[integrator.chosen_])>0)
            std::swap(first,last);
        const auto midpoint=MappedHermiteInterface::between(first.local[integrator.chosen_],last.local[integrator.chosen_],.5);
        const int sign=integrator.integrationRoot(midpoint).derivativeSign;
        require(sign!=0,"geometry moment primitive lost its graph orientation");
        const auto primitive=[&](const Cut& cut) {
            const auto point=integrator.physical(cut.local);
            return product(product(point[0],point[1]),axisymmetric?point[0]:M::constant(1.0));
        };
        const long double endpoint=sign*(axisymmetric?3.1415926535897932384626433832795L:1.L)
            *M::approximate(difference(primitive(last),primitive(first)));
        if(graph.coordinateIsX) moments.radial=moments.axial-endpoint;
        else moments.axial=moments.radial-endpoint;
        // Apply the original moment accuracy scale also to the recovered
        // component: cancellation in the endpoint primitive must not allow
        // an inaccurate small moment to pass on the larger integral's scale.
        const long double fraction=M::approximate(difference(last.local[integrator.chosen_],first.local[integrator.chosen_]));
        const long double capacity=M::approximate(integrator.width_[0])*M::approximate(integrator.width_[1])
            *(axisymmetric?6.283185307179586476925286766559L
                *M::approximate(sum(integrator.origin_[0],integrator.width_[0])):1.L);
        const long double limit=2e-12L*capacity*fraction
            +2e-10L*std::min(std::abs(moments.axial),std::abs(moments.radial));
        require(std::isfinite(moments.axial) && std::isfinite(moments.radial)
            && moments.errorEstimate<=limit,"recovered Green moment exceeds controlled quadrature accuracy");
        return moments;
    }
    struct RegionGeometryResult {Result integral;std::size_t coincidentArcs,localRootRefinements;};
    // Retained supports are internal strict-TU data, not a public flag that
    // can make an unproved hand-written graph coincide with a clip edge.
    static RegionGeometryResult regionGeometry(const typename CellLayout::Source& source,
        const Public::Graph& graph,const std::vector<const typename CellLayout::Query*>& cuts,
        const std::vector<Public::Region>& regions,bool axisymmetric,Vector2D greenOrigin)
    {
        require(!graph.nativeSide,"side-attached region consumption requires the explicit trace-volume consumer");
        require(graph.controlledIntegration,"region moments require a controlled complete graph");
        CellIntegrator integrator(source,graph,cuts,regions,axisymmetric,greenOrigin);
        integrator.geometryMeasure_=true;
        std::vector<CoincidentSupport> supports;
        if constexpr(nativeQuarticCertificates) {
            const int normal=1-integrator.chosen_;
            for(const auto* query:cuts) for(const auto& zero:query->zeroTraces) {
                const bool incident=std::find(zero.cells.begin(),zero.cells.end(),integrator.layout_.address())!=zero.cells.end();
                if(!incident)continue;
                Vector2D a=query->inputBegin,b=query->inputEnd;
                if(std::make_pair(b.x,b.y)<std::make_pair(a.x,a.y))std::swap(a,b);
                const int axis=a.x==b.x?0:1;
                if(axis!=normal)continue; // Not a graph-coordinate support.
                auto first=integrator.segmentLocal(a,b,zero.chart.begin),last=integrator.segmentLocal(a,b,zero.chart.end);
                const auto lower=point(graph.localLower),upper=point(graph.localUpper);
                if(M::compare(first[normal],lower[normal])<=0 || M::compare(first[normal],upper[normal])>=0)continue;
                if(M::compare(last[integrator.chosen_],lower[integrator.chosen_])<=0
                    || M::compare(first[integrator.chosen_],upper[integrator.chosen_])>=0)continue;
                // The ORIGINAL complete zero-set certificate excludes another
                // transverse root and a critical point anywhere on this chart.
                require(source.resolvedZeroSetPartial(graph.nativeCellIndex,normal,graph.localLower,graph.localUpper),
                    "coincident clip lacks a complete regular graph certificate");
                supports.push_back({std::move(first),std::move(last),normal});
            }
        }
        integrator.coincidentSupports_=&supports;
        auto result=integrator.run();
        return {std::move(result),integrator.coincidentArcs_,integrator.localRootRefinements_+integrator.transverseRootRefinements_};
    }
    // A batch shares graph topology and every sampled root/physical jet.
    // Target integrators below supply only original ownership and direct
    // Green gauges; their graph integrations are never run separately.
    static std::vector<RegionGeometryResult> regionGeometryBatch(
        const typename CellLayout::Source& source,const Public::Graph& graph,
        const std::vector<const typename CellLayout::Query*>& cuts,
        const std::vector<Public::NativeMomentTarget>& targets,bool axisymmetric,
        Public::GraphMoments* nativeGraphMoments=nullptr,
        const std::vector<std::vector<int>>* selectedRegionIds=nullptr,
        bool guardedTransverseNewton=false,Public::NativeMomentWorkProfile* work=nullptr)
    {
        require(!targets.empty(),"empty native region moment batch");
        require(!selectedRegionIds || (!nativeGraphMoments && selectedRegionIds->size()==targets.size()),
            "selected native batch has a different channel inventory");
        require(!graph.nativeSide,"side-attached batch requires the explicit trace-volume consumer");
        require(graph.controlledIntegration,"batch region moments require a controlled complete graph");
        CellIntegrator shared(source,graph,cuts,targets.front().regions,axisymmetric,targets.front().greenOrigin);
        shared.geometryMeasure_=true;
        shared.work_=work;
        if constexpr(std::is_same_v<CellLayout,NativeQuarticCurveCell>)
            shared.guardedTransverseNewton_=guardedTransverseNewton && (selectedRegionIds || nativeGraphMoments);
        std::vector<std::unique_ptr<CellIntegrator>> owned;
        std::vector<CellIntegrator*> channels;
        owned.reserve(targets.size()-1+(nativeGraphMoments?1:0));
        channels.reserve(targets.size()+(nativeGraphMoments?1:0));
        // shared already owns target 0's actual constructor validation,
        // geometry measure and Green origin; retain the original channel order.
        channels.push_back(&shared);
        for(std::size_t t=1;t<targets.size();++t) {
            const auto& target=targets[t];
            owned.push_back(std::unique_ptr<CellIntegrator>(new CellIntegrator(
                source,graph,cuts,target.regions,axisymmetric,target.greenOrigin)));
            owned.back()->geometryMeasure_=true;channels.push_back(owned.back().get());
        }
        BatchOwnership ownership;
        const auto& partition=targets.front().regions;
        bool common=true;
        for(const auto& target:targets) {
            if(target.regions.size()!=partition.size()){common=false;break;}
            for(std::size_t r=0;r<partition.size() && common;++r)for(int v=0;v<4;++v)
                if(target.regions[r].vertices[v].x!=partition[r].vertices[v].x
                    || target.regions[r].vertices[v].y!=partition[r].vertices[v].y){common=false;break;}
        }
        if(common)for(std::size_t r=0;r<partition.size();++r) {
            Vector2D lo=partition[r].vertices[0],hi=lo;
            for(const auto& v:partition[r].vertices) {
                lo.x=std::min(lo.x,v.x);lo.y=std::min(lo.y,v.y);
                hi.x=std::max(hi.x,v.x);hi.y=std::max(hi.y,v.y);
            }
            bool rectangular=lo.x<hi.x && lo.y<hi.y;
            for(int v=0;v<4;++v) {
                const auto a=partition[r].vertices[v],b=partition[r].vertices[(v+1)%4];
                rectangular=rectangular && (a.x==lo.x || a.x==hi.x) && (a.y==lo.y || a.y==hi.y)
                    && ((a.x==b.x)!=(a.y==b.y));
            }
            if(!rectangular){common=false;break;}
            typename BatchOwnership::Rectangle box{point(lo),point(hi),{}};
            box.owners.reserve(targets.size());
            for(const auto& target:targets)box.owners.push_back(target.regions[r].regionId);
            ownership.rectangles.push_back(std::move(box));
        }
        if(common)shared.batchOwnership_=&ownership;
        std::vector<CoincidentSupport> supports;
        if constexpr(nativeQuarticCertificates) {
            const int normal=1-shared.chosen_;
            for(const auto* query:cuts)for(const auto& zero:query->zeroTraces) {
                const bool incident=std::find(zero.cells.begin(),zero.cells.end(),shared.layout_.address())!=zero.cells.end();
                if(!incident)continue;
                Vector2D a=query->inputBegin,b=query->inputEnd;
                if(std::make_pair(b.x,b.y)<std::make_pair(a.x,a.y))std::swap(a,b);
                const int axis=a.x==b.x?0:1;
                if(axis!=normal)continue;
                auto first=shared.segmentLocal(a,b,zero.chart.begin),last=shared.segmentLocal(a,b,zero.chart.end);
                const auto lower=point(graph.localLower),upper=point(graph.localUpper);
                if(M::compare(first[normal],lower[normal])<=0 || M::compare(first[normal],upper[normal])>=0)continue;
                if(M::compare(last[shared.chosen_],lower[shared.chosen_])<=0
                    || M::compare(first[shared.chosen_],upper[shared.chosen_])>=0)continue;
                require(source.resolvedZeroSetPartial(graph.nativeCellIndex,normal,graph.localLower,graph.localUpper),
                    "coincident batch clip lacks a complete regular graph certificate");
                supports.push_back({std::move(first),std::move(last),normal});
            }
        }
        shared.coincidentSupports_=&supports;
        for(auto* channel:channels)channel->coincidentSupports_=&supports;
        if(!nativeGraphMoments)
            return shared.withLocalRootResolution([&] { return shared.runBatchOnce(channels,selectedRegionIds); });
        // This extra channel integrates the COMPLETE original graph. Actual
        // target cuts only supply shared split/root information; they never
        // remove a global arc or give it a local clip owner.
        const std::vector<Public::Region> fullGraphRegions;
        owned.push_back(std::unique_ptr<CellIntegrator>(new CellIntegrator(
            source,graph,cuts,fullGraphRegions,GeometryOnly{})));
        auto* global=owned.back().get();
        global->axisymmetric_=axisymmetric;global->geometryMeasure_=true;
        channels.push_back(global);
        auto output=shared.withLocalRootResolution([&] { return shared.runBatchOnce(channels); });
        Public::GraphMoments moments;
        for(const auto& arc:output.back().integral.arcs) {
            moments.radial+=arc.radialMoment;moments.axial+=arc.axialMoment;
            moments.errorEstimate+=arc.momentErrorEstimate;
        }
        output.pop_back();
        // Recover the global other component ONCE from the original complete
        // graph's retained endpoints, exactly as the independent nativeMoments
        // path. No local accumulated moment or rounded gauge shift is used.
        auto first=shared.endpoint(graph.begin),last=shared.endpoint(graph.end);
        if(M::compare(first.local[shared.chosen_],last.local[shared.chosen_])>0)
            std::swap(first,last);
        const auto midpoint=MappedHermiteInterface::between(first.local[shared.chosen_],last.local[shared.chosen_],.5);
        const int sign=shared.integrationRoot(midpoint).derivativeSign;
        require(sign!=0,"geometry moment primitive lost its graph orientation");
        const auto primitive=[&](const Cut& cut) {
            const auto point=shared.physical(cut.local);
            return product(product(point[0],point[1]),axisymmetric?point[0]:M::constant(1.0));
        };
        const long double endpoint=sign*(axisymmetric?3.1415926535897932384626433832795L:1.L)
            *M::approximate(difference(primitive(last),primitive(first)));
        if(graph.coordinateIsX)moments.radial=moments.axial-endpoint;
        else moments.axial=moments.radial-endpoint;
        const long double fraction=M::approximate(difference(last.local[shared.chosen_],first.local[shared.chosen_]));
        const long double capacity=M::approximate(shared.width_[0])*M::approximate(shared.width_[1])
            *(axisymmetric?6.283185307179586476925286766559L
                *M::approximate(sum(shared.origin_[0],shared.width_[0])):1.L);
        const long double limit=2e-12L*capacity*fraction
            +2e-10L*std::min(std::abs(moments.axial),std::abs(moments.radial));
        require(std::isfinite(moments.axial) && std::isfinite(moments.radial)
            && moments.errorEstimate<=limit,"recovered Green moment exceeds controlled quadrature accuracy");
        *nativeGraphMoments=moments;
        return output;
    }
    // Strict-TU discovery uses these SAME edge inventories in original-cell
    // coordinates. They are not public floating points or another scalar.
    static std::vector<Root> chartEdgeRoots(const typename CellLayout::Source& source,
        const Public::Graph& graph,int side)
    {
        const std::vector<const typename CellLayout::Query*> cuts;
        const std::vector<Public::Region> regions;
        return CellIntegrator(source,graph,cuts,regions,GeometryOnly{}).edgeRoots(side);
    }
    static std::vector<int> chartEdgeIncidences(const typename CellLayout::Source& source,
        const Public::Graph& graph,int side)
    {
        const std::vector<const typename CellLayout::Query*> cuts;
        const std::vector<Public::Region> regions;
        CellIntegrator sampler(source,graph,cuts,regions,GeometryOnly{});
        std::vector<int> counts;
        for(const auto& root:sampler.edgeRoots(side)) {
            if(root.phaseBefore!=0 && root.phaseAfter!=0 && root.phaseBefore!=root.phaseAfter) {
                counts.push_back(1);continue;
            }
            int inward=0;
            if constexpr(nativeQuarticCertificates) {
                const auto lo=sampler.edgePoint(side,root.lower),hi=sampler.edgePoint(side,root.upper);
                const int normal=(side==0 || side==2)?0:1;
                const auto p=sampler.layout_.polynomial(
                    sampler.coordinate(lo[0],hi[0]),sampler.coordinate(lo[1],hi[1]),normal==0,normal==1);
                // A nonzero derivative over the COMPLETE root bracket proves
                // the local implicit graph. No rounded root/gradient test.
                try {
                    const auto zeros=PolynomialSegmentRoots::isolateEndpointRelative(p);
                    if(!zeros.identicallyZero && zeros.roots.empty())
                        inward=(side<2?1:-1)*Exact::sign(p[0]);
                } catch(const std::runtime_error&) {} // no certificate; another corner side may certify
            }
            if(!inward){counts.push_back(-1);continue;}
            const int phase=root.phaseBefore?root.phaseBefore:root.phaseAfter;
            require(phase!=0,"edge root has no adjacent phase");
            counts.push_back(phase==inward?0:(root.phaseBefore && root.phaseAfter?2:1));
        }
        return counts;
    }
    using GraphSample=MappedHermiteInterface::GraphSample;
    // Complete endpoint-root inventories belong to the immutable source and
    // ORIGINAL graph, not to a query fraction. This object owns their selected
    // Cuts and graph identity; the source owner must outlive it. It is not an
    // integral/region object and cannot import privately refined face roots.
    class GraphSampling
    {
        friend class CellIntegrator;
        const typename CellLayout::Source* source_;
        Public::Graph graph_;
        Cut first_,last_;
        GraphSampling(const typename CellLayout::Source& source,const Public::Graph& graph,Cut first,Cut last)
            : source_(&source),graph_(graph),first_(std::move(first)),last_(std::move(last)) {}
    };
    static GraphSampling prepareGraphSampling(const typename CellLayout::Source& source,const Public::Graph& graph)
    {
        const std::vector<const typename CellLayout::Query*> cuts;
        const std::vector<Public::Region> regions;
        CellIntegrator sampler(source,graph,cuts,regions,GeometryOnly{});
        if(graph.nativeSide){auto ends=sampler.sideEndpoints();
            return GraphSampling(source,graph,std::move(ends[0]),std::move(ends[1]));}
        auto ends=sampler.resolvedSamplingEndpoints(0,1,.5);
        return GraphSampling(source,graph,std::move(ends[0]),std::move(ends[1]));
    }
    // Terminal local-chart sample of the ORIGINAL graph. Topology remains
    // its side/root identity; converting a Ratio to its nearest-end tail is
    // a floating evaluation point, not an exact root or a global inverse.
    static GraphSample graphSample(const typename CellLayout::Source& source,
        const Public::Graph& graph,double fraction)
    { return graphSample(source,graph,0.0,1.0,fraction); }

    static GraphSample graphSample(const typename CellLayout::Source& source,
        const Public::Graph& graph,double begin,double end,double fraction)
    {
        validateSamplingInterval(begin,end,fraction);
        const auto prepared=prepareGraphSampling(source,graph);
        return graphSample(source,prepared,begin,end,fraction);
    }

    static GraphSample graphSample(const typename CellLayout::Source& source,
        const GraphSampling& prepared,double begin,double end,double fraction)
    {
        Exact::validateEnvironment();
        validateSamplingInterval(begin,end,fraction);
        require(prepared.source_==&source,"prepared graph sampling belongs to another scalar source");
        const auto& graph=prepared.graph_;
        const std::vector<const typename CellLayout::Query*> cuts;
        const std::vector<Public::Region> regions;
        CellIntegrator sampler(source,graph,cuts,regions,GeometryOnly{});
        return contextual([&]() -> GraphSample {
        const Cut* first=&prepared.first_;
        const Cut* last=&prepared.last_;
        const int chosen=sampler.chosen_;
        auto firstChosen=MappedHermiteInterface::between(first->local[chosen],last->local[chosen],begin);
        auto lastChosen=MappedHermiteInterface::between(first->local[chosen],last->local[chosen],end);
        const bool originalBegin=begin==0 && fraction==0;
        const bool originalEnd=end==1 && fraction==1;
        // Prepared endpoints normally suffice. Only an unresolved interior
        // query refines its own canonical roots; the shared cache stays immutable.
        std::optional<std::array<Cut,2>> refined;
        if (!originalBegin && !originalEnd && !sampler.separatedSample(*first,*last,
                MappedHermiteInterface::between(firstChosen,lastChosen,fraction))) {
            refined=sampler.resolvedSamplingEndpoints(begin,end,fraction);
            first=&(*refined)[0];last=&(*refined)[1];
            firstChosen=MappedHermiteInterface::between(first->local[chosen],last->local[chosen],begin);
            lastChosen=MappedHermiteInterface::between(first->local[chosen],last->local[chosen],end);
        }
        const auto extent=difference(lastChosen,firstChosen);
        const int direction=Exact::sign(extent.numerator);
        require(direction!=0,"graph sample requires a nonzero chosen-coordinate span");
        Point p=originalEnd?last->local:first->local;
        if(!originalBegin && !originalEnd)
        {
            const auto q=MappedHermiteInterface::between(firstChosen,lastChosen,fraction);
            require(sampler.separatedSample(*first,*last,q),
                "graph sampling point remains inside an endpoint root enclosure");
            p[chosen]=q;
            if(!graph.nativeSide){const auto root=sampler.sectionRoot(q);p[1-chosen]=root.parameter;}
        }
        const auto coordinate=[](const Ratio& value) {
            require(M::compare(value,M::constant(0))>=0 && M::compare(value,M::constant(1))<=0,
                "sampled graph point left its closed cell");
            const bool upper=M::compare(value,M::constant(.5))>0;
            const auto tail=upper?difference(M::constant(1),value):value;
            const double offset=static_cast<double>(M::approximate(tail));
            require(std::isfinite(offset) && offset>=0 && offset<=.5 &&
                (Exact::sign(tail.numerator)==0 || offset>0),"native graph tail collapsed during conversion");
            return CellChartCoordinate{offset,upper};
        };
        GraphSample result;
        result.chart={coordinate(p[0]),coordinate(p[1])};
        const auto firstCoordinate=coordinate(firstChosen),lastCoordinate=coordinate(lastChosen);
        require(compareCellChartCoordinate(lastCoordinate,firstCoordinate)==direction,
            "native graph endpoints collapse in the terminal local chart");
        if(fraction>0 && fraction<1)
        {
            const auto current=chosen==0?result.chart.x:result.chart.y;
            require(compareCellChartCoordinate(current,firstCoordinate)*direction>0 &&
                compareCellChartCoordinate(lastCoordinate,current)*direction>0,
                "graph quadrature coordinate is not representably interior in its local chart");
        }
        result.physicalSpan=std::abs(M::approximate(product(extent,sampler.width_[chosen])));
        result.coordinateDirection=direction;
        require(std::isfinite(result.physicalSpan) && result.physicalSpan>0,"invalid native graph physical span");
        return result;
        },[&] { return sampler.samplingFailureContext(begin,end,fraction); });
    }

private:
    static void validateSamplingInterval(double begin,double end,double fraction)
    {
        require(std::isfinite(fraction) && fraction>=0 && fraction<=1,"invalid graph sampling fraction");
        require(std::isfinite(begin) && std::isfinite(end) && begin>=0 && end<=1 && begin<end,
            "invalid native subarc interval");
    }
public:

    // Invert only the retained chosen-coordinate chart, never a physical
    // display point. Exact endpoint charts have their original identities;
    // an interior quotient rounding to 0/1 must reject, not snap to an end.
    static double graphFraction(const typename CellLayout::Source& source,const Public::Graph& graph,
        const CellChartPoint& point)
    {
        validateCellChart(point);
        const auto beginChart=graphSample(source,graph,0).chart;
        const auto endChart=graphSample(source,graph,1).chart;
        if(cellChartEqual(point,beginChart)) return 0;
        if(cellChartEqual(point,endChart)) return 1;
        if(!graph.nativeSide)require(cellChartInterior(point),"native interior material point is outside its cell");
        else {
            const auto& side=*graph.nativeSide;const auto normal=side.normalAxis()==0?point.x:point.y;
            require(compareCellChartCoordinate(normal,{side.side()<2?0.:1.,false})==0,
                "native side graph point left its certified normal support");
        }
        const std::vector<const typename CellLayout::Query*> cuts;
        const std::vector<Public::Region> regions;
        CellIntegrator sampler(source,graph,cuts,regions,GeometryOnly{});
        const auto prepared=prepareGraphSampling(source,graph);
        const auto& first=prepared.first_;
        const auto& last=prepared.last_;
        const int chosen=sampler.chosen_;
        const auto local=chosen==0?point.x:point.y;
        const auto q=local.fromUpper?difference(M::constant(1),M::constant(local.offset)):M::constant(local.offset);
        const auto fraction=quotient(difference(q,first.local[chosen]),
            difference(last.local[chosen],first.local[chosen]));
        require(M::compare(fraction,M::constant(0))>0 && M::compare(fraction,M::constant(1))<0,
            "native material point is outside its retained graph");
        const double result=static_cast<double>(M::approximate(fraction));
        require(std::isfinite(result) && result>0 && result<1,
            "native interior material point has no represented graph parameter");
        return result;
    }

    // Orient the original graph using the SAME transverse-root and jet
    // evaluation as GL8. This samples orientation; run() still validates all
    // quadrature sections and ownership. It is not a projection surrogate.
private:
    void validateNativeSide(const typename CellLayout::Source& source) const
    {
        if constexpr(std::is_same_v<CellLayout,NativeQuarticCurveCell>) {
            const auto& side=*graph_.nativeSide;const int normal=side.normalAxis();
            require(side.usesNativeScalar(&source) && side.ownerCell()==graph_.nativeCellIndex,
                "native side graph belongs to another scalar/cell");
            require(graph_.cellI==-1 && graph_.cellJ==-1 && graph_.controlledIntegration
                && same(point(graph_.localLower),point({0,0})) && same(point(graph_.localUpper),point({1,1}))
                && graph_.coordinateIsX==(normal==1) && graph_.begin.side==-1 && graph_.begin.rootOrdinal==-1
                && graph_.end.side==-1 && graph_.end.rootOrdinal==-1,"native side graph has inconsistent public addressing");
            const auto ends=sideEndpoints();
            require(M::compare(ends[0].local[normal],M::constant(side.side()<2?0:1))==0
                && M::compare(ends[1].local[normal],ends[0].local[normal])==0,"native side endpoints left their actual leaf side");
            for(const auto& end:ends)for(int a=0;a<2;++a)require(M::compare(end.local[a],M::constant(0))>=0
                && M::compare(end.local[a],M::constant(1))<=0,"native side endpoint left its original cell");
            require(M::compare(ends[0].local[chosen_],ends[1].local[chosen_])!=0,"zero native side extent");
            const auto value=source.restrictPhysical(graph_.nativeCellIndex,side.begin(),side.end()).numerator;
            require(std::all_of(value.begin(),value.end(),[](const auto& c){return Exact::sign(c)==0;}),"native side value support is not exact zero");
            require(source.resolvedZeroSetPartial(graph_.nativeCellIndex,normal),"native side lost its complete signed normal zero-set certificate");
            const auto derivative=source.restrictPhysical(graph_.nativeCellIndex,side.begin(),side.end(),normal==0,normal==1).numerator;
            const auto roots=PolynomialSegmentRoots::isolateEndpointRelative(derivative);
            require(!roots.identicallyZero && roots.roots.empty(),"native side lost strict first normal derivative");
            const int sign=Exact::sign(derivative[0]);
            require((chosen_==0?-sign:sign)*M::compare(ends[1].local[chosen_],ends[0].local[chosen_])>0,
                "native side support reverses the original interface orientation");
        } else require(false,"native side sampling requires native Q4");
    }
    std::array<Cut,2> sideEndpoints() const
    {
        const auto& side=*graph_.nativeSide;
        const auto a=segmentLocal(side.begin(),side.end(),M::constant(0)),b=segmentLocal(side.begin(),side.end(),M::constant(1));
        return {{Cut{a,a[chosen_],a[chosen_]},Cut{b,b[chosen_],b[chosen_]}}};
    }
    Public::GraphMoments sideMoments(const typename CellLayout::Source& source,bool axisymmetric) const
    {
        Public::GraphMoments result;
        if constexpr(std::is_same_v<CellLayout,NativeQuarticCurveCell>) {
            const auto& side=*graph_.nativeSide;const int normal=side.normalAxis();
            const auto derivative=source.restrictPhysical(graph_.nativeCellIndex,side.begin(),side.end(),normal==0,normal==1).numerator;
            const int sign=Exact::sign(derivative[0]);
            auto a=point(side.begin()),b=point(side.end());if(M::compare(a[chosen_],b[chosen_])>0)std::swap(a,b);
            require(!axisymmetric || (Exact::sign(a[0].numerator)>=0 && Exact::sign(b[0].numerator)>=0),"negative native side RZ radius");
            constexpr long double pi=3.1415926535897932384626433832795L;
            if(normal==0)result.radial=sign*(axisymmetric?pi:1.L)*M::approximate(product(
                axisymmetric?product(a[0],a[0]):a[0],difference(b[1],a[1])));
            else result.axial=sign*(axisymmetric?pi:1.L)*M::approximate(product(a[1],axisymmetric?
                product(difference(b[0],a[0]),sum(b[0],a[0])):difference(b[0],a[0])));
            require(std::isfinite(result.radial) && std::isfinite(result.axial),"nonfinite native side Green moment");
        }
        return result;
    }
    int graphDirection()
    {
        return contextual([&] {
        const auto ends=resolvedSamplingEndpoints(0,1,.5);
        const auto& begin=ends[0];const auto& end=ends[1];
        const int order=M::compare(end.local[chosen_],begin.local[chosen_]);
        require(order!=0,"graph endpoints have no chosen-coordinate extent");
        const auto q=MappedHermiteInterface::between(begin.local[chosen_],end.local[chosen_],0.5);
        const auto root=sectionRoot(q);
        Point p; p[chosen_]=q; p[1-chosen_]=root.parameter;
        const auto sample=jet(p);
        const long double derivative=chosen_==0?sample.dy:sample.dx;
        const long double roundoff=chosen_==0?sample.dyRoundoff:sample.dxRoundoff;
        require(root.derivativeSign!=0 && std::abs(derivative)>roundoff &&
            (derivative>0?1:-1)==root.derivativeSign,"unresolved free-graph orientation");
        return order*(chosen_==0?-root.derivativeSign:root.derivativeSign);
        },[&] { return samplingFailureContext(0,1,.5); });
    }

public:
    Result run()
    { return withLocalRootResolution([&] { return runOnce(); }); }

private:
    template<class Action> auto withLocalRootResolution(Action action) -> decltype(action())
    {
        // A bounded numerical-resolution policy, not a physical tolerance.
        // Only requested, face-independent canonical endpoints are refined
        // here. Public face roots are always handled by the caller's cache.
        constexpr std::size_t maximumLocalRefinements=16;
        std::size_t refinements=0;
        for (;;)
        {
            try { return action(); }
            catch (const LocalRootResolutionError& request)
            {
                require(!request.roots.empty(),"empty local endpoint-resolution request");
                for (const auto& identity:request.roots)
                {
                    require(refinements<maximumLocalRefinements,
                        "canonical endpoint refinement limit exhausted; no arc was accepted");
                    refineEdgeRoot(identity);
                    ++refinements;
                    ++localRootRefinements_;
                }
            }
        }
    }

    bool separatedSample(const Cut& first,const Cut& last,const Ratio& q) const
    {
        const int order=M::compare(last.local[chosen_],first.local[chosen_]);
        if (order==0) return false;
        const auto& low=order>0?first:last;
        const auto& high=order>0?last:first;
        return M::compare(low.upper,q)<0 && M::compare(q,high.lower)<0;
    }
    std::array<Cut,2> resolvedSamplingEndpoints(double begin,double end,double fraction)
    {
        return withLocalRootResolution([&] {
            std::array<Cut,2> ends{{endpoint(graph_.begin),endpoint(graph_.end)}};
            const auto a=MappedHermiteInterface::between(ends[0].local[chosen_],ends[1].local[chosen_],begin);
            const auto b=MappedHermiteInterface::between(ends[0].local[chosen_],ends[1].local[chosen_],end);
            if (!separatedSample(ends[0],ends[1],MappedHermiteInterface::between(a,b,fraction)))
                requestResolution({&ends[0],&ends[1]},"graph interior sample overlaps an endpoint root enclosure");
            return ends;
        });
    }
    std::string scalarFailureContext() const
    {
        // Failure-only data allow a local scalar/chart reproducer; no diagnostic
        // file IO, mutable geometry or extra sampling on accepted queries.
        std::ostringstream out;
        out<<std::setprecision(17)<<"[MappedInterface cell="<<graph_.nativeCellIndex
            <<", chosen="<<chosen_<<", lower=("<<graph_.localLower.x<<','<<graph_.localLower.y
            <<"), upper=("<<graph_.localUpper.x<<','<<graph_.localUpper.y<<"), ends=("
            <<graph_.begin.side<<','<<graph_.begin.rootOrdinal<<';'<<graph_.end.side<<','<<graph_.end.rootOrdinal
            <<"), localDofs={";
        for (double value:dofs_) out<<value<<',';
        return out.str()+"}]";
    }
    std::string samplingFailureContext(double begin,double end,double fraction) const
    {
        std::ostringstream out;
        out<<std::setprecision(17)<<scalarFailureContext()<<" [sampling interval=("
            <<begin<<','<<end<<"), fraction="<<fraction<<']';
        return out.str();
    }

private:
    Result runOnce()
    {
        Result result;
        coincidentArcs_=0;
        edgeBindings_.clear();
        std::vector<Cut> rawCuts;
        for (std::size_t segmentIndex=0;segmentIndex<cuts_.size();++segmentIndex)
        {
            const auto* query=cuts_[segmentIndex];
            require(query!=nullptr, "null retained face cut");
            Vector2D a=query->inputBegin,b=query->inputEnd;
            const bool inputReversed=std::make_pair(b.x,b.y)<std::make_pair(a.x,a.y);
            if (inputReversed) std::swap(a,b);
            for (std::size_t rootIndex=0;rootIndex<query->roots.size();++rootIndex)
            {
                const auto& root=query->roots[rootIndex];
                if (!(layout_.address()==M::CellAddress{root.info.cellI,root.info.cellJ,root.info.nativeCellIndex})) continue;
                ++result.sourceRootIncidences;
                contextual([&] {
                const auto convert=[&](const Ratio& t) {
                    Point local;
                    for (int axis=0;axis<2;++axis)
                    {
                        const double av=axis==0?a.x:a.y,bv=axis==0?b.x:b.y;
                        local[axis]=quotient(difference(sum(M::constant(av),
                            product(M::ratio(M::difference(bv,av),M::number(1.0)),t)),origin_[axis]),width_[axis]);
                    }
                    return local;
                };
                Cut cut=stage("face root physical-to-local conversion",[&] {
                    return Cut{convert(root.parameter),convert(root.lower)[chosen_],convert(root.upper)[chosen_]};
                });
                // A physical FV face can cross other charts of this SAME
                // cell. Exclude its roots only when the retained enclosure
                // lies outside this closed subbox. Ambiguous ordering asks
                // the caller to refine its shared root, never snaps it.
                const auto boxLo=point(graph_.localLower),boxHi=point(graph_.localUpper);
                const auto enclosureA=convert(root.lower),enclosureB=convert(root.upper);
                bool outside=false,uncertain=false;
                for (int axis=0;axis<2;++axis) {
                    auto lo=enclosureA[axis],hi=enclosureB[axis];
                    if (M::compare(lo,hi)>0) std::swap(lo,hi);
                    outside=outside || M::compare(hi,boxLo[axis])<0 || M::compare(lo,boxHi[axis])>0;
                    uncertain=uncertain || (M::compare(lo,boxLo[axis])<0 && M::compare(hi,boxLo[axis])>=0)
                        || (M::compare(lo,boxHi[axis])<=0 && M::compare(hi,boxHi[axis])>0);
                }
                if (outside) return;
                cut.sources.push_back({segmentIndex,inputReversed?query->roots.size()-1-rootIndex:rootIndex});
                if (M::compare(cut.lower,cut.upper)>0) std::swap(cut.lower,cut.upper);
                for (int axis=0;axis<2;++axis)
                    require(M::compare(cut.local[axis],M::constant(0.0))>=0
                        && M::compare(cut.local[axis],M::constant(1.0))<=0,
                        "retained root incidence left its source cell");

                // A cut supported by a canonical cell edge is the SAME edge
                // polynomial. Its enclosure must intersect exactly one root
                // of the complete edge inventory; this establishes ordinal
                // identity without comparing rounded representatives.
                for (int side=0;side<4;++side)
                {
                    const int axis=(side==0 || side==2)?0:1;
                    const auto edgePosition=sum(origin_[axis],product(width_[axis],sidePosition(side)));
                    if (M::compare(M::constant(axis==0?a.x:a.y),edgePosition)!=0
                        || M::compare(M::constant(axis==0?b.x:b.y),edgePosition)!=0) continue;
                    const int tangent=1-axis;
                    auto lo=convert(root.lower)[tangent],hi=convert(root.upper)[tangent];
                    if (M::compare(lo,hi)>0) std::swap(lo,hi);
                    const auto& edge=stage("face/supporting-edge ordinal match",[&]() -> const std::vector<Root>& {
                        return edgeRoots(side);
                    });
                    int match=-1;
                    for (int k=0;k<static_cast<int>(edge.size());++k)
                        if (M::compare(lo,edge[k].upper)<=0 && M::compare(edge[k].lower,hi)<=0)
                        {
                            require(match<0, "face/edge root ordinal remains ambiguous");
                            match=k;
                        }
                    if (match<0 && uncertain) throw Public::RootResolutionError(
                        "face root enclosure straddles an internal graph boundary",cut.sources);
                    require(match>=0, "face root is missing from the complete supporting edge inventory");
                    // Intersect two enclosures of this SAME complete-edge
                    // root. This also recognizes exact chart-corner roots:
                    // endlessly refining a face bracket across that exact
                    // corner would never establish strict one-sided order.
                    if (M::compare(edge[match].lower,lo)>0) lo=edge[match].lower;
                    if (M::compare(edge[match].upper,hi)<0) hi=edge[match].upper;
                    if (M::compare(cut.local[tangent],lo)<0 || M::compare(cut.local[tangent],hi)>0)
                        cut.local[tangent]=M::compare(edge[match].parameter,lo)>=0 &&
                            M::compare(edge[match].parameter,hi)<=0?edge[match].parameter:
                            MappedHermiteInterface::between(lo,hi,.5);
                    cut.lower=chosen_==tangent?lo:cut.local[chosen_];
                    cut.upper=chosen_==tangent?hi:cut.local[chosen_];
                    cut.edgeSide=side; cut.edgeOrdinal=match;
                    bindEdgeRoot(cut,std::move(lo),std::move(hi));
                    break;
                }
                if (uncertain && cut.edgeSide<0) throw Public::RootResolutionError(
                    "face root enclosure straddles an internal graph boundary",cut.sources);
                rawCuts.push_back(std::move(cut));
                },[&] {
                    std::ostringstream out;
                    out << std::setprecision(17) << "[MappedInterface face root " << rootIndex
                        << ", cell=(" << layout_.address() << "), segment=("
                        << a.x << ',' << a.y << ")->(" << b.x << ',' << b.y << "), chart="
                        << root.info.chart << ", t=" << describe(root.parameter) << ']';
                    return out.str();
                });
            }
        }
        // Resolve graph endpoints only AFTER binding actual face roots to
        // canonical ordinals. A coarse edge inventory never overwrites an
        // already refined, same-identity face root representation.
        Cut first=stage("begin graph endpoint",[&] { return endpoint(graph_.begin); });
        Cut last=stage("end graph endpoint",[&] { return endpoint(graph_.end); });
        require(M::compare(first.local[chosen_],last.local[chosen_])!=0,
            "graph has no positive chosen-coordinate extent");
        if (M::compare(first.local[chosen_],last.local[chosen_])>0) std::swap(first,last);
        std::vector<Cut> split{first,last};
        if(coincidentSupports_)for(const auto& support:*coincidentSupports_)for(const auto* p:{&support.begin,&support.end}) {
            if(M::compare((*p)[chosen_],first.local[chosen_])<=0 || M::compare((*p)[chosen_],last.local[chosen_])>=0)continue;
            const auto lower=point(graph_.localLower),upper=point(graph_.localUpper);
            if(M::compare((*p)[chosen_],lower[chosen_])<0 || M::compare((*p)[chosen_],upper[chosen_])>0)continue;
            split.push_back(Cut{*p,(*p)[chosen_],(*p)[chosen_]});
        }
        for (auto cut:rawCuts)
        {
            if (cut.edgeSide>=0) cut=endpoint({cut.edgeSide,cut.edgeOrdinal});
            if (M::compare(cut.local[chosen_],first.local[chosen_])<0)
            {
                if (M::compare(cut.upper,first.lower)>=0)
                    requestResolution({&cut,&first},"face cut and graph begin have unresolved ordering");
                continue;
            }
            if (M::compare(cut.local[chosen_],last.local[chosen_])>0)
            {
                if (M::compare(cut.lower,last.upper)<=0)
                    requestResolution({&cut,&last},"face cut and graph end have unresolved ordering");
                continue;
            }
            split.push_back(std::move(cut));
        }
        stage("sort retained chosen-coordinate cuts",[&] {
            std::sort(split.begin(),split.end(),[&](const Cut& a,const Cut& b) {
                return M::compare(a.local[chosen_],b.local[chosen_])<0;
            });
        });
        std::vector<Cut> unique;
        for (const auto& cut:split)
        {
            if (!unique.empty() && M::compare(unique.back().local[chosen_],cut.local[chosen_])==0)
            {
                require(same(unique.back().local,cut.local), "different curve points share one graph coordinate");
                mergeSources(unique.back().sources,cut.sources);
                if (M::compare(cut.lower,unique.back().lower)>0) unique.back().lower=cut.lower;
                if (M::compare(cut.upper,unique.back().upper)<0) unique.back().upper=cut.upper;
                continue;
            }
            if (!unique.empty() && M::compare(unique.back().upper,cut.lower)>=0)
                requestResolution({&unique.back(),&cut},"independent cut-root enclosures cannot be strictly ordered");
            unique.push_back(cut);
        }
        result.retainedCuts=unique.size();
        for (std::size_t k=1;k<unique.size();++k)
        {
            if(coincidentBoundary(unique[k-1],unique[k])) {++coincidentArcs_;continue;}
            result.arcs.push_back(contextual([&] {
                return integrateArc(unique[k-1],unique[k]);
            },[&] {
                return "[MappedInterface subarc "+std::to_string(k-1)+", chosen="+std::to_string(chosen_)
                    +", begin="+describe(unique[k-1].local[chosen_])+", end="+describe(unique[k].local[chosen_])+"]";
            }));
        }
        return result;
    }

    std::vector<RegionGeometryResult> runBatchOnce(const std::vector<CellIntegrator*>& channels,
        const std::vector<std::vector<int>>* selectedRegionIds=nullptr)
    {
        Result result;
        coincidentArcs_=0;
        edgeBindings_.clear();
        std::vector<Cut> rawCuts;
        for (std::size_t segmentIndex=0;segmentIndex<cuts_.size();++segmentIndex)
        {
            const auto* query=cuts_[segmentIndex];
            require(query!=nullptr, "null retained face cut");
            Vector2D a=query->inputBegin,b=query->inputEnd;
            const bool inputReversed=std::make_pair(b.x,b.y)<std::make_pair(a.x,a.y);
            if (inputReversed) std::swap(a,b);
            for (std::size_t rootIndex=0;rootIndex<query->roots.size();++rootIndex)
            {
                const auto& root=query->roots[rootIndex];
                if (!(layout_.address()==M::CellAddress{root.info.cellI,root.info.cellJ,root.info.nativeCellIndex})) continue;
                ++result.sourceRootIncidences;
                contextual([&] {
                const auto convert=[&](const Ratio& t) {
                    Point local;
                    for (int axis=0;axis<2;++axis)
                    {
                        const double av=axis==0?a.x:a.y,bv=axis==0?b.x:b.y;
                        local[axis]=quotient(difference(sum(M::constant(av),
                            product(M::ratio(M::difference(bv,av),M::number(1.0)),t)),origin_[axis]),width_[axis]);
                    }
                    return local;
                };
                Cut cut=stage("face root physical-to-local conversion",[&] {
                    return Cut{convert(root.parameter),convert(root.lower)[chosen_],convert(root.upper)[chosen_]};
                });
                // A physical FV face can cross other charts of this SAME
                // cell. Exclude its roots only when the retained enclosure
                // lies outside this closed subbox. Ambiguous ordering asks
                // the caller to refine its shared root, never snaps it.
                const auto boxLo=point(graph_.localLower),boxHi=point(graph_.localUpper);
                const auto enclosureA=convert(root.lower),enclosureB=convert(root.upper);
                bool outside=false,uncertain=false;
                for (int axis=0;axis<2;++axis) {
                    auto lo=enclosureA[axis],hi=enclosureB[axis];
                    if (M::compare(lo,hi)>0) std::swap(lo,hi);
                    outside=outside || M::compare(hi,boxLo[axis])<0 || M::compare(lo,boxHi[axis])>0;
                    uncertain=uncertain || (M::compare(lo,boxLo[axis])<0 && M::compare(hi,boxLo[axis])>=0)
                        || (M::compare(lo,boxHi[axis])<=0 && M::compare(hi,boxHi[axis])>0);
                }
                if (outside) return;
                cut.sources.push_back({segmentIndex,inputReversed?query->roots.size()-1-rootIndex:rootIndex});
                if (M::compare(cut.lower,cut.upper)>0) std::swap(cut.lower,cut.upper);
                for (int axis=0;axis<2;++axis)
                    require(M::compare(cut.local[axis],M::constant(0.0))>=0
                        && M::compare(cut.local[axis],M::constant(1.0))<=0,
                        "retained root incidence left its source cell");

                // A cut supported by a canonical cell edge is the SAME edge
                // polynomial. Its enclosure must intersect exactly one root
                // of the complete edge inventory; this establishes ordinal
                // identity without comparing rounded representatives.
                for (int side=0;side<4;++side)
                {
                    const int axis=(side==0 || side==2)?0:1;
                    const auto edgePosition=sum(origin_[axis],product(width_[axis],sidePosition(side)));
                    if (M::compare(M::constant(axis==0?a.x:a.y),edgePosition)!=0
                        || M::compare(M::constant(axis==0?b.x:b.y),edgePosition)!=0) continue;
                    const int tangent=1-axis;
                    auto lo=convert(root.lower)[tangent],hi=convert(root.upper)[tangent];
                    if (M::compare(lo,hi)>0) std::swap(lo,hi);
                    const auto& edge=stage("face/supporting-edge ordinal match",[&]() -> const std::vector<Root>& {
                        return edgeRoots(side);
                    });
                    int match=-1;
                    for (int k=0;k<static_cast<int>(edge.size());++k)
                        if (M::compare(lo,edge[k].upper)<=0 && M::compare(edge[k].lower,hi)<=0)
                        {
                            require(match<0, "face/edge root ordinal remains ambiguous");
                            match=k;
                        }
                    if (match<0 && uncertain) throw Public::RootResolutionError(
                        "face root enclosure straddles an internal graph boundary",cut.sources);
                    require(match>=0, "face root is missing from the complete supporting edge inventory");
                    // Intersect two enclosures of this SAME complete-edge
                    // root. This also recognizes exact chart-corner roots:
                    // endlessly refining a face bracket across that exact
                    // corner would never establish strict one-sided order.
                    if (M::compare(edge[match].lower,lo)>0) lo=edge[match].lower;
                    if (M::compare(edge[match].upper,hi)<0) hi=edge[match].upper;
                    if (M::compare(cut.local[tangent],lo)<0 || M::compare(cut.local[tangent],hi)>0)
                        cut.local[tangent]=M::compare(edge[match].parameter,lo)>=0 &&
                            M::compare(edge[match].parameter,hi)<=0?edge[match].parameter:
                            MappedHermiteInterface::between(lo,hi,.5);
                    cut.lower=chosen_==tangent?lo:cut.local[chosen_];
                    cut.upper=chosen_==tangent?hi:cut.local[chosen_];
                    cut.edgeSide=side; cut.edgeOrdinal=match;
                    bindEdgeRoot(cut,std::move(lo),std::move(hi));
                    break;
                }
                if (uncertain && cut.edgeSide<0) throw Public::RootResolutionError(
                    "face root enclosure straddles an internal graph boundary",cut.sources);
                rawCuts.push_back(std::move(cut));
                },[&] {
                    std::ostringstream out;
                    out << std::setprecision(17) << "[MappedInterface face root " << rootIndex
                        << ", cell=(" << layout_.address() << "), segment=("
                        << a.x << ',' << a.y << ")->(" << b.x << ',' << b.y << "), chart="
                        << root.info.chart << ", t=" << describe(root.parameter) << ']';
                    return out.str();
                });
            }
        }
        // Resolve graph endpoints only AFTER binding actual face roots to
        // canonical ordinals. A coarse edge inventory never overwrites an
        // already refined, same-identity face root representation.
        Cut first=stage("begin graph endpoint",[&] { return endpoint(graph_.begin); });
        Cut last=stage("end graph endpoint",[&] { return endpoint(graph_.end); });
        require(M::compare(first.local[chosen_],last.local[chosen_])!=0,
            "graph has no positive chosen-coordinate extent");
        if (M::compare(first.local[chosen_],last.local[chosen_])>0) std::swap(first,last);
        std::vector<Cut> split{first,last};
        if(coincidentSupports_)for(const auto& support:*coincidentSupports_)for(const auto* p:{&support.begin,&support.end}) {
            if(M::compare((*p)[chosen_],first.local[chosen_])<=0 || M::compare((*p)[chosen_],last.local[chosen_])>=0)continue;
            const auto lower=point(graph_.localLower),upper=point(graph_.localUpper);
            if(M::compare((*p)[chosen_],lower[chosen_])<0 || M::compare((*p)[chosen_],upper[chosen_])>0)continue;
            split.push_back(Cut{*p,(*p)[chosen_],(*p)[chosen_]});
        }
        for (auto cut:rawCuts)
        {
            if (cut.edgeSide>=0) cut=endpoint({cut.edgeSide,cut.edgeOrdinal});
            if (M::compare(cut.local[chosen_],first.local[chosen_])<0)
            {
                if (M::compare(cut.upper,first.lower)>=0)
                    requestResolution({&cut,&first},"face cut and graph begin have unresolved ordering");
                continue;
            }
            if (M::compare(cut.local[chosen_],last.local[chosen_])>0)
            {
                if (M::compare(cut.lower,last.upper)<=0)
                    requestResolution({&cut,&last},"face cut and graph end have unresolved ordering");
                continue;
            }
            split.push_back(std::move(cut));
        }
        stage("sort retained chosen-coordinate cuts",[&] {
            std::sort(split.begin(),split.end(),[&](const Cut& a,const Cut& b) {
                return M::compare(a.local[chosen_],b.local[chosen_])<0;
            });
        });
        std::vector<Cut> unique;
        for (const auto& cut:split)
        {
            if (!unique.empty() && M::compare(unique.back().local[chosen_],cut.local[chosen_])==0)
            {
                require(same(unique.back().local,cut.local), "different curve points share one graph coordinate");
                mergeSources(unique.back().sources,cut.sources);
                if (M::compare(cut.lower,unique.back().lower)>0) unique.back().lower=cut.lower;
                if (M::compare(cut.upper,unique.back().upper)<0) unique.back().upper=cut.upper;
                continue;
            }
            if (!unique.empty() && M::compare(unique.back().upper,cut.lower)>=0)
                requestResolution({&unique.back(),&cut},"independent cut-root enclosures cannot be strictly ordered");
            unique.push_back(cut);
        }
        result.retainedCuts=unique.size();
        std::vector<RegionGeometryResult> output;
        output.reserve(channels.size());
        for(std::size_t t=0;t<channels.size();++t)output.push_back({Result{},0,0});
        for(std::size_t k=1;k<unique.size();++k) {
            std::vector<bool> active(channels.size());bool any=false;
            for(std::size_t t=0;t<channels.size();++t) {
                active[t]=!channels[t]->coincidentBoundary(unique[k-1],unique[k]);
                any=any || active[t];
                if(!active[t])++output[t].coincidentArcs;
            }
            if(selectedRegionIds && any) {
                // Complete retained face cuts prove one owner on this open
                // subarc. Locate it with an actual enclosed midpoint root;
                // selected GL8 nodes retain all original ownership checks.
                const auto q=MappedHermiteInterface::between(unique[k-1].local[chosen_],unique[k].local[chosen_],.5);
                const bool leftBlocked=M::compare(unique[k-1].upper,q)>=0;
                const bool rightBlocked=M::compare(q,unique[k].lower)>=0;
                if(leftBlocked || rightBlocked)
                    requestResolution({leftBlocked?&unique[k-1]:nullptr,rightBlocked?&unique[k]:nullptr},
                        "selected batch midpoint lies inside a cutting-root enclosure");
                const auto sample=batchRegionIntegrationRoot(q,channels,active);
                any=false;
                for(std::size_t t=0;t<channels.size();++t)if(active[t]) {
                    const int owner=sample.owners[t];
                    require(owner!=-2,"selected batch midpoint has overlapping target owners");
                    require(owner>=0,"selected batch midpoint has no actual target owner");
                    require(sample.lowerOwners[t]==owner && sample.upperOwners[t]==owner,
                        "selected batch midpoint root enclosure crosses an actual target boundary");
                    const auto& ids=(*selectedRegionIds)[t];
                    active[t]=std::find(ids.begin(),ids.end(),owner)!=ids.end();
                    any=any || active[t];
                }
            }
            if(!any)continue;
            auto arcs=contextual([&] {
                return integrateArcBatch(unique[k-1],unique[k],channels,active);
            },[&] {
                return "[MappedInterface batch subarc "+std::to_string(k-1)+", chosen="+std::to_string(chosen_)
                    +", begin="+describe(unique[k-1].local[chosen_])+", end="+describe(unique[k].local[chosen_])+"]";
            });
            for(std::size_t t=0;t<channels.size();++t)if(active[t])
                output[t].integral.arcs.push_back(std::move(arcs[t]));
        }
        for(auto& target:output) {
            target.integral.retainedCuts=result.retainedCuts;
            target.integral.sourceRootIncidences=result.sourceRootIncidences;
            target.localRootRefinements=localRootRefinements_+transverseRootRefinements_;
        }
        return output;
    }

private:
    Point segmentLocal(Vector2D a,Vector2D b,const Ratio& t) const
    {
        Point local;
        for(int axis=0;axis<2;++axis) {
            const double av=axis?a.y:a.x,bv=axis?b.y:b.x;
            local[axis]=quotient(difference(sum(M::constant(av),product(
                M::ratio(M::difference(bv,av),M::number(1)),t)),origin_[axis]),width_[axis]);
        }
        return local;
    }
    bool coincidentBoundary(const Cut& first,const Cut& last) const
    {
        if(!coincidentSupports_)return false;
        for(const auto& s:*coincidentSupports_) {
            if(M::compare(first.local[chosen_],s.begin[chosen_])<0 || M::compare(last.local[chosen_],s.end[chosen_])>0)continue;
            // Exact support + complete signed normal partial proved the
            // unique zero graph throughout this positive chosen interval.
            Point a=first.local,b=last.local;a[s.normalAxis]=b[s.normalAxis]=s.begin[s.normalAxis];
            const auto pa=physical(a),pb=physical(b);
            int matches=0,region=-1;bool different=false;
            for(const auto& r:regions_)for(int e=0;e<4;++e) {
                const auto va=r.vertices[e],vb=r.vertices[(e+1)%4];
                if(orientation(va,vb,pa)!=0 || orientation(va,vb,pb)!=0)continue;
                const auto lo=M::constant(std::min(chosen_?va.y:va.x,chosen_?vb.y:vb.x));
                const auto hi=M::constant(std::max(chosen_?va.y:va.x,chosen_?vb.y:vb.x));
                if(M::compare(pa[chosen_],lo)<0 || M::compare(pb[chosen_],hi)>0)continue;
                ++matches;different=different || (region>=0 && region!=r.regionId);region=r.regionId;
            }
            require(matches>0 && matches<=2,"coincident clip has missing/overlapping region boundary support");
            if(matches==1 || different)return true;
        }
        return false;
    }
    CellLayout layout_;
    const Public::Graph& graph_;
    const std::vector<const typename CellLayout::Query*>& cuts_;
    const std::vector<Public::Region>& regions_;
    int chosen_;
    Vector2D nominal_;
    Point origin_,width_;
    std::array<double,dofCount> dofs_;
    std::array<std::vector<Root>,4> edgeRoots_;
    std::array<bool,4> edgeReady_{{false,false,false,false}};
    std::map<EdgeIdentity,EdgeBinding> edgeBindings_;
    std::vector<int> regionOrientation_;

    static void mergeSources(std::vector<Public::RootResolutionRequest>& target,
        const std::vector<Public::RootResolutionRequest>& incoming)
    {
        target.insert(target.end(),incoming.begin(),incoming.end());
        const auto key=[](const Public::RootResolutionRequest& r) { return std::make_pair(r.segmentIndex,r.rootIndex); };
        std::sort(target.begin(),target.end(),[&](const auto& a,const auto& b) { return key(a)<key(b); });
        target.erase(std::unique(target.begin(),target.end(),[&](const auto& a,const auto& b) { return key(a)==key(b); }),target.end());
    }
    void bindEdgeRoot(const Cut& cut,Ratio lower,Ratio upper)
    {
        const int tangent=(cut.edgeSide==0 || cut.edgeSide==2)?1:0;
        const EdgeIdentity identity{cut.edgeSide,cut.edgeOrdinal};
        const Root candidate{cut.local[tangent],std::move(lower),std::move(upper)};
        const auto inserted=edgeBindings_.emplace(identity,EdgeBinding{candidate,cut.sources});
        if (inserted.second) return;
        auto& binding=inserted.first->second;
        const Ratio lo=M::compare(binding.root.lower,candidate.lower)>0?binding.root.lower:candidate.lower;
        const Ratio hi=M::compare(binding.root.upper,candidate.upper)<0?binding.root.upper:candidate.upper;
        require(M::compare(lo,hi)<=0,"same canonical-edge root has disjoint face-source enclosures");
        const auto inside=[&](const Ratio& value) { return M::compare(lo,value)<=0 && M::compare(value,hi)<=0; };
        if (!inside(binding.root.parameter))
            binding.root.parameter=inside(candidate.parameter)?candidate.parameter:MappedHermiteInterface::between(lo,hi,0.5);
        binding.root.lower=lo; binding.root.upper=hi;
        mergeSources(binding.sources,cut.sources);
        // Identity was proved by the same source's unique complete-edge
        // ordinal, not closeness of coordinates. The intersection chooses a
        // common integral representative of that algebraic root; individual
        // face measures remain floating approximations within their own
        // retained enclosures. This is not an exact-watertightness claim.
    }
    [[noreturn]] void requestResolution(std::initializer_list<const Cut*> blocked,
        const std::string& reason) const
    {
        std::vector<Public::RootResolutionRequest> requests;
        std::vector<EdgeIdentity> local;
        for (const auto* cut:blocked)
        {
            if (!cut || M::compare(cut->lower,cut->upper)==0) continue;
            if (!cut->sources.empty()) mergeSources(requests,cut->sources);
            else if (cut->edgeSide>=0) local.emplace_back(cut->edgeSide,cut->edgeOrdinal);
        }
        const std::string message="[MappedInterface boundary-root resolution] "+reason;
        // External roots go first so every dependent face measure sees the
        // caller's new shared cache. Any remaining local-only endpoint is
        // revisited when that complete integration is restarted.
        if (!requests.empty()) throw Public::RootResolutionError(message,std::move(requests));
        std::sort(local.begin(),local.end());
        local.erase(std::unique(local.begin(),local.end()),local.end());
        if (!local.empty()) throw LocalRootResolutionError(message,std::move(local));
        throw std::runtime_error(message+"; no positive-width refinable source enclosure");
    }
    Point edgePoint(int side,const Ratio& tangent) const
    {
        const bool vertical=side==0 || side==2;
        Point p;
        p[vertical?0:1]=sidePosition(side);
        p[vertical?1:0]=tangent;
        return p;
    }
    Ratio sidePosition(int side) const
    {
        return M::constant(side==0?graph_.localLower.x:side==1?graph_.localLower.y:
            side==2?graph_.localUpper.x:graph_.localUpper.y);
    }
    Root mapRoot(Root root,const Ratio& lower,const Ratio& upper) const
    {
        NativeWorkTimer timer(work_,&Public::NativeMomentWorkProfile::rootMappingRefinement);
        if(work_)++work_->rootMappings;
        const auto extent=difference(upper,lower);
        const auto map=[&](const Ratio& t){return sum(lower,product(extent,t));};
        root.parameter=map(root.parameter);root.lower=map(root.lower);root.upper=map(root.upper);
        root.geometryPacket.reset();
        return root;
    }
    void refineEdgeRoot(const EdgeIdentity& identity)
    {
        require(edgeBindings_.find(identity)==edgeBindings_.end(),
            "a face-backed canonical root cannot be privately refined");
        const auto& inventory=edgeRoots(identity.first);
        require(identity.second>=0 && identity.second<static_cast<int>(inventory.size()),
            "local endpoint refinement ordinal is outside its complete inventory");
        auto& root=edgeRoots_[identity.first][identity.second];
        NativeWorkTimer timer(work_,&Public::NativeMomentWorkProfile::rootMappingRefinement,
            &Public::NativeMomentWorkProfile::restrictionConstruction);
        if(work_)++work_->rootRefinements;
        const auto refined=MappedHermiteRootRefinement::refine(root.lower,root.upper,
            [&](const Ratio& lo,const Ratio& hi) {
                auto p=polynomial(edgePoint(identity.first,lo),edgePoint(identity.first,hi));
                if constexpr(nativeQuarticCertificates)
                    if(root.deflated)p=PolynomialSquareFree::squareFree(p).polynomial;
                return p;
            });
        require(root.deflated || root.derivativeSign==0 || refined.derivativeSign==root.derivativeSign,
            "canonical endpoint refinement changed its simple-root crossing direction");
        root.parameter=refined.originalParameter;root.lower=refined.originalLower;root.upper=refined.originalUpper;
        if(!root.deflated)root.derivativeSign=refined.derivativeSign;
    }

    M::Coordinate coordinate(const Ratio& a, const Ratio& b) const
    {
        // Identical denominators already define one common chart. Avoid
        // inserting an unnecessary D/D factor, especially on a transverse
        // line whose fixed coordinate has a==b. This is exact algebra, not
        // approximate division or cancellation of expansion tails.
        M::Coordinate c;
        if (Exact::sign(Exact::subtract(a.denominator,b.denominator))==0)
            c={a.numerator,b.numerator,a.denominator};
        else
            c={Exact::multiply(a.numerator,b.denominator),
                Exact::multiply(b.numerator,a.denominator),Exact::multiply(a.denominator,b.denominator)};
        int exponent=0;
        std::frexp(Exact::estimate(c.denominator),&exponent);
        c.numeratorBegin=Exact::scalePow2(c.numeratorBegin,-exponent);
        c.numeratorEnd=Exact::scalePow2(c.numeratorEnd,-exponent);
        c.denominator=Exact::scalePow2(c.denominator,-exponent);
        return c;
    }
    typename CellLayout::Polynomial sourcePolynomial(const M::Coordinate& x,const M::Coordinate& y) const
    {
        NativeWorkTimer timer(work_,&Public::NativeMomentWorkProfile::restrictionConstruction);
        if(work_)++work_->restrictionQueries;
        return layout_.polynomial(x,y);
    }
    typename CellLayout::Polynomial polynomial(const Point& a, const Point& b) const
    { return stage("cleared restriction construction",[&] {
        return sourcePolynomial(coordinate(a[0],b[0]),coordinate(a[1],b[1]));
    }); }
    std::vector<Root> isolate(Point a, Point b,bool allowMultiple=false,bool guardedNewton=false) const
    {
        return contextual([&] {
            if constexpr(std::is_same_v<CellLayout,NativeQuarticCurveCell>) {
                using Certificate=typename CellLayout::Source::InitialRootCertificate;
                // Retain ALL eight raw expansion vectors before endpoint
                // reversal/normalization. The source leaf owns this cache;
                // graph/region/face state never enters isolateImpl.
                const std::array<Expansion,8> input{{a[0].numerator,a[0].denominator,
                    a[1].numerator,a[1].denominator,b[0].numerator,b[0].denominator,
                    b[1].numerator,b[1].denominator}};
                std::vector<Certificate> certificates;
                if(work_)++work_->initialRootCacheQueries;
                const bool retained=[&] {
                    NativeWorkTimer timer(work_,&Public::NativeMomentWorkProfile::initialRootCacheLookup);
                    return layout_.lookupInitialRoots(input,allowMultiple,certificates);
                }();
                if(work_) {
                    if(retained)++work_->initialRootCacheHits;else ++work_->initialRootCacheMisses;
                }
                if(retained) {
                    std::vector<Root> result;result.reserve(certificates.size());
                    for(const auto& c:certificates)result.push_back({
                        {c.parameterNumerator,c.parameterDenominator},
                        {c.lowerNumerator,c.lowerDenominator},{c.upperNumerator,c.upperDenominator},
                        c.derivativeSign,c.phaseBefore,c.phaseAfter,c.deflated});
                    return result;
                }
                auto result=isolateImpl(a,b,allowMultiple,guardedNewton);
                certificates.reserve(result.size());
                for(const auto& r:result)certificates.push_back({
                    r.parameter.numerator,r.parameter.denominator,r.lower.numerator,r.lower.denominator,
                    r.upper.numerator,r.upper.denominator,r.derivativeSign,r.phaseBefore,r.phaseAfter,r.deflated});
                // Only complete successful initial inventories are retained.
                // mapRoot, bindings and all subsequent refinements stay private.
                {
                    NativeWorkTimer timer(work_,&Public::NativeMomentWorkProfile::initialRootCacheRemember);
                    if(work_)++work_->initialRootCacheRememberQueries;
                    layout_.rememberInitialRoots(input,allowMultiple,certificates);
                }
                return result;
            } else return isolateImpl(a,b,allowMultiple,guardedNewton);
        },[&] {
            return "[MappedInterface isolate, a="+describe(a)+", b="+describe(b)+"]";
        });
    }
    // P0.1 bounded dyadic fallback for a one-shot near-critical failure on ONE
    // native Q4 chart edge (a dyadic subbox side or a transverse GL8 section).
    // This is the chart-edge counterpart to the physical-face recovery in
    // LevelSetMappedSegment: when the interface runs almost collinear with the
    // sampled edge the degree-eight restriction can carry a tiny extremum that
    // does not cross zero yet defeats the single-shot interval isolator
    // (Zalesak, rotated slot wall grazing a subbox edge). Re-restricting the
    // SAME retained scalar along dyadic subedges is a Bernstein subdivision: it
    // only tightens hulls, rounds no coefficient and adds/removes no root. A
    // multiple/tangential root on a native Q4 boundary edge can instead use
    // its exact square-free certificate. Transverse integration sections
    // remain strict: their derivative must belong to the original scalar.
    // Other unresolved cases retain the bounded subdivision rejection.
    PolynomialSegmentRoots::Result isolateEdgePieceWithDyadicFallback(
        const Point& a,const Point& b,const typename CellLayout::Polynomial& p,bool allowMultiple,bool& deflated,
        bool guardedNewton=false) const
    {
        using PolynomialSegmentRoots::isolateEndpointRelative;
        try {
            if constexpr(std::is_same_v<CellLayout,NativeQuarticCurveCell>)
                if(guardedNewton && !allowMultiple)
                    return PolynomialSegmentRoots::isolateEndpointRelativeGuardedNewton(p);
            return isolateEndpointRelative(p);
        }
        catch (const std::exception& outer)
        {
            if constexpr(nativeQuarticCertificates) {
                if(allowMultiple) {
                    std::optional<typename CellLayout::Polynomial> factor;
                    try {const auto result=PolynomialSquareFree::squareFree(p);if(result.deflated)factor=result.polynomial;}
                    catch(const std::runtime_error&) {} // preserve bounded same-scalar subdivision below
                    if(factor){auto roots=isolateEndpointRelative(*factor);deflated=true;return roots;}
                }
            }
            if constexpr (std::tuple_size<typename CellLayout::Polynomial>::value != 9) throw;
            const int fixedAxis=(M::compare(a[0],b[0])==0)?0:1;
            const int varyAxis=1-fixedAxis;
            constexpr int maxDepth=6;
            PolynomialSegmentRoots::Result merged;
            std::vector<bool> exactSplitRoots;
            bool failed=false;
            std::function<void(const Ratio&,const Ratio&,double,double,int)> recurse=
                [&](const Ratio& ra,const Ratio& rb,double lo,double hi,int depth)
            {
                if(failed) return;
                Point pa=a,pb=b; pa[varyAxis]=ra; pb[varyAxis]=rb;
                try
                {
                    const auto q=polynomial(pa,pb);
                    auto sub=isolateEndpointRelative(q);
                    if(sub.identicallyZero){failed=true;return;}
                    for(std::size_t k=0;k<sub.roots.size();++k)
                    {
                        const bool endpoint=sub.roots[k]==0 || sub.roots[k]==1;
                        if (endpoint && ((sub.roots[k]==0 && lo>0) || (sub.roots[k]==1 && hi<1)))
                            require(MappedHermiteRootRefinement::endpointDerivativeSign(q,sub.roots[k]==1)!=0,
                                "dyadic subdivision exposed a multiple root; physical refinement is required");
                        const auto mapped=[&](double t) {
                            return Exact::add(M::number(lo),Exact::multiply(M::number(hi-lo),M::number(t)));
                        };
                        merged.roots.push_back(lo+(hi-lo)*sub.roots[k]);
                        merged.enclosures.push_back({std::max(lo,Exact::enclosure(mapped(sub.enclosures[k].lo)).lo),
                            std::min(hi,Exact::enclosure(mapped(sub.enclosures[k].hi)).hi)});
                        exactSplitRoots.push_back(endpoint);
                    }
                }
                catch(const std::exception&)
                {
                    if(depth>=maxDepth){failed=true;return;}
                    const Ratio mid=MappedHermiteInterface::between(ra,rb,0.5);
                    recurse(ra,mid,lo,0.5*(lo+hi),depth+1);
                    recurse(mid,rb,0.5*(lo+hi),hi,depth+1);
                }
            };
            recurse(a[varyAxis],b[varyAxis],0.0,1.0,0);
            if(failed)
                throw std::runtime_error(std::string("dyadic edge root-isolation budget exhausted: ")+outer.what());
            std::vector<std::size_t> order(merged.roots.size());
            std::iota(order.begin(),order.end(),std::size_t{0});
            std::sort(order.begin(),order.end(),
                [&](std::size_t x,std::size_t y){return merged.roots[x]<merged.roots[y];});
            PolynomialSegmentRoots::Result sorted;
            bool previousExact=false;
            for(std::size_t idx:order)
            {
                const double r=merged.roots[idx];
                auto e=merged.enclosures[idx];
                if(!sorted.enclosures.empty() && e.lo<=sorted.enclosures.back().hi)
                {
                    require(previousExact && exactSplitRoots[idx] && r==sorted.roots.back(),
                        "distinct dyadic root enclosures overlap; no root was merged");
                    continue;
                }
                sorted.roots.push_back(r);
                sorted.enclosures.push_back(e);
                previousExact=exactSplitRoots[idx];
            }
            sorted.identicallyZero=false;
            return sorted;
        }
    }

    std::vector<Root> isolateImpl(Point a, Point b,bool allowMultiple,bool guardedNewton) const
    {
        // Select the endpoint-relative direction before constructing all
        // segment coefficients. Endpoint input coordinates retain the forward
        // segment's cleared coordinate denominators. The native backend
        // independently centers each returned homogeneous numerator, so these
        // endpoint restrictions and the full restriction need not share a
        // numerator scale; full-polynomial endpoint values cannot replace them.
        // A whole transverse line has an endpoint basis in its varying axis,
        // so these constant restrictions are much cheaper than two full ones.
        const auto endpoints=stage("cleared restriction construction",[&] {
            const auto x=coordinate(a[0],b[0]),y=coordinate(a[1],b[1]);
            const auto atEndpoint=[&](bool last) {
                auto cx=x,cy=y;
                cx.numeratorBegin=cx.numeratorEnd=last?x.numeratorEnd:x.numeratorBegin;
                cy.numeratorBegin=cy.numeratorEnd=last?y.numeratorEnd:y.numeratorBegin;
                return sourcePolynomial(cx,cy)[0];
            };
            return std::array<Expansion,2>{{atEndpoint(false),atEndpoint(true)}};
        });
        const auto absolute=[](Expansion e) { return Exact::sign(e)<0?times(e,-1.0):e; };
        const bool reverse=Exact::sign(Exact::subtract(absolute(endpoints[1]),absolute(endpoints[0])))<0;
        if (reverse) std::swap(a,b);
        const auto p=polynomial(a,b);
        bool deflated=false;
        const auto isolated=stage("endpoint-relative polynomial root isolation",[&] {
            NativeWorkTimer timer(work_,&Public::NativeMomentWorkProfile::initialRootIsolation,
                &Public::NativeMomentWorkProfile::restrictionConstruction);
            if(work_)++work_->initialRootIsolations;
            return isolateEdgePieceWithDyadicFallback(a,b,p,allowMultiple,deflated,guardedNewton);
        });
        require(!isolated.identicallyZero, "coincident zero restriction is unsupported");
        require(isolated.roots.size()==isolated.enclosures.size(), "root/enclosure inventory mismatch");
        int phase=0;
        for (const auto& coefficient:p) if ((phase=Exact::sign(coefficient))!=0) break;
        require(phase!=0, "nonzero restriction has no one-sided endpoint phase");
        std::vector<Root> result;
        for (std::size_t k=0;k<isolated.roots.size();++k)
        {
            const double t=isolated.roots[k];
            const auto enclosure=isolated.enclosures[k];
            require(t>=0.0 && t<=1.0 && enclosure.lo>=0.0 && enclosure.hi<=1.0
                && enclosure.lo<=t && t<=enclosure.hi, "invalid transverse root enclosure");
            require((t!=0.0 || Exact::sign(p[0])==0)
                && (t!=1.0 || Exact::sign(Exact::evaluate(p,1.0))==0),
                "transverse root collapsed to a nonzero endpoint");
            Root root{M::constant(t),M::constant(enclosure.lo),M::constant(enclosure.hi)};
            root.deflated=deflated;
            root.phaseBefore=t==0?0:phase;
            if (t>0.0 && t<1.0)
            {
                // Only the original strict inventory certifies simplicity.
                // A square-free edge inventory supplies original phases but
                // never a surrogate derivative for the original scalar.
                if(!deflated){root.derivativeSign=reverse?phase:-phase;phase=-phase;}
                else {
                    const double next=k+1<isolated.roots.size()?isolated.enclosures[k+1].lo:1.;
                    const double witness=.5*enclosure.hi+.5*next;
                    require(witness>enclosure.hi && witness<next,"no original phase witness after multiple edge root");
                    phase=Exact::sign(Exact::evaluate(p,witness));
                    require(phase!=0,"edge phase witness is an original-polynomial zero");
                }
            }
            root.phaseAfter=t==1?0:phase;
            if (reverse)
            {
                root.parameter=difference(M::constant(1.0),root.parameter);
                const auto lower=difference(M::constant(1.0),root.upper);
                root.upper=difference(M::constant(1.0),root.lower); root.lower=lower;
                std::swap(root.phaseBefore,root.phaseAfter);
            }
            result.push_back(std::move(root));
        }
        const int lastPhase=Exact::sign(Exact::evaluate(p,1.0));
        require(lastPhase==0 || lastPhase==phase, "root phase inventory disagrees with the final endpoint");
        if (reverse) std::reverse(result.begin(),result.end());
        return result;
    }
    const std::vector<Root>& edgeRoots(int side)
    {
        require(side>=0 && side<4, "endpoint side must be W/S/E/N");
        if (!edgeReady_[side])
        {
            const bool vertical=side==0 || side==2;
            Point a=point(graph_.localLower),b=point(graph_.localUpper);
            a[vertical?0:1]=b[vertical?0:1]=sidePosition(side);
            edgeRoots_[side]=contextual([&] { return isolate(a,b,true); },[&] {
                return "[MappedInterface complete edge inventory, side="+std::to_string(side)+"]";
            });
            for (auto& root:edgeRoots_[side]) root=mapRoot(std::move(root),a[vertical?1:0],b[vertical?1:0]);
            edgeReady_[side]=true;
        }
        return edgeRoots_[side];
    }
    Cut endpoint(const Public::Endpoint& endpoint)
    {
        const auto& roots=edgeRoots(endpoint.side);
        require(endpoint.rootOrdinal>=0 && endpoint.rootOrdinal<static_cast<int>(roots.size()),
            "endpoint root ordinal is outside the complete edge inventory");
        const auto binding=edgeBindings_.find({endpoint.side,endpoint.rootOrdinal});
        const auto& root=binding==edgeBindings_.end()?roots[endpoint.rootOrdinal]:binding->second.root;
        const bool vertical=endpoint.side==0 || endpoint.side==2;
        Point p;
        p[vertical?0:1]=sidePosition(endpoint.side);
        p[vertical?1:0]=root.parameter;
        Cut result{p,chosen_==(vertical?1:0)?root.lower:p[chosen_],
            chosen_==(vertical?1:0)?root.upper:p[chosen_]};
        result.edgeSide=endpoint.side; result.edgeOrdinal=endpoint.rootOrdinal;
        if (binding!=edgeBindings_.end()) result.sources=binding->second.sources;
        return result;
    }
    bool nativeSideCovered(Vector2D a,Vector2D b) const
    {
        // Native coarse sides are supplied as their real vertex-midpoint
        // half faces, not a new whole-side query. Prove exact collinearity and
        // gap-free coverage using represented physical endpoints. No root,
        // scalar, subface or small interval is replaced or merged here.
        const bool x=a.x!=b.x;
        const auto coordinate=[&](Vector2D p) { return x?p.x:p.y; };
        const double low=std::min(coordinate(a),coordinate(b));
        const double high=std::max(coordinate(a),coordinate(b));
        std::vector<std::pair<double,double>> intervals;
        for(const auto* query:cuts_)
        {
            require(query!=nullptr,"null retained face cut");
            const auto c=query->inputBegin,d=query->inputEnd;
            if(orientation(a,b,point(c))!=0 || orientation(a,b,point(d))!=0) continue;
            const double first=std::min(coordinate(c),coordinate(d));
            const double last=std::max(coordinate(c),coordinate(d));
            if(first>=low && last<=high && first<last) intervals.emplace_back(first,last);
        }
        std::sort(intervals.begin(),intervals.end());
        double reach=low;
        for(const auto& interval:intervals)
        {
            if(interval.first>reach) return false;
            reach=std::max(reach,interval.second);
        }
        return reach==high;
    }
    void validateRegions()
    {
        for (const auto& region:regions_)
        {
            require(region.regionId>=0, "negative ownership region id");
            for (const auto& p:region.vertices)
                require(std::isfinite(p.x) && std::isfinite(p.y) && (!axisymmetric_ || p.x>=0.0),
                    "region has a nonfinite coordinate or negative radius");
            int sign=0;
            for (int k=0;k<4;++k)
            {
                const int turn=orientation(region.vertices[k],region.vertices[(k+1)%4],point(region.vertices[(k+2)%4]));
                require(turn!=0 && (sign==0 || sign==turn), "region is not a strict convex cyclic quadrilateral");
                sign=turn;
            }
            regionOrientation_.push_back(sign);
        }
        const auto equal=[](Vector2D a,Vector2D b) { return a.x==b.x && a.y==b.y; };
        // Same-owner subquad seams need no curve cuts. Every other side must
        // be an actual supplied segment. Native layouts can instead supply
        // an exactly covering sequence of their real collinear half faces.
        for (std::size_t r=0;r<regions_.size();++r) for (int e=0;e<4;++e)
        {
            const auto a=regions_[r].vertices[e],b=regions_[r].vertices[(e+1)%4];
            bool covered=false;
            for (std::size_t s=0;s<regions_.size() && !covered;++s)
                if (s!=r && regions_[s].regionId==regions_[r].regionId)
                    for (int f=0;f<4;++f)
                    {
                        const auto c=regions_[s].vertices[f],d=regions_[s].vertices[(f+1)%4];
                        if ((equal(a,c)&&equal(b,d)) || (equal(a,d)&&equal(b,c))) covered=true;
                    }
            for (const auto* query:cuts_) if (!covered)
            {
                require(query!=nullptr, "null retained face cut");
                if ((equal(a,query->inputBegin)&&equal(b,query->inputEnd))
                    || (equal(a,query->inputEnd)&&equal(b,query->inputBegin))) covered=true;
            }
            if(!covered && graph_.nativeCellIndex>=0) covered=nativeSideCovered(a,b);
            require(covered, "an external region side is missing its actual retained cut segment");
        }
    }
    int uniqueOwner(const Point& p) const
    {
        if (geometryMeasure_ && regions_.empty()) return 0;
        int owner=-1;
        for (std::size_t k=0;k<regions_.size();++k)
        {
            bool inside=true;
            for (int e=0;e<4 && inside;++e)
                inside=orientation(regions_[k].vertices[e],regions_[k].vertices[(e+1)%4],p)*regionOrientation_[k]>=0;
            if (!inside) continue;
            if(owner>=0 && owner!=regions_[k].regionId) return -2;
            owner=regions_[k].regionId;
        }
        return owner;
    }
    int owner(const Point& p) const
    {
        const int result=uniqueOwner(p);
        require(result!=-2,"free-curve point has different overlapping region owners");
        require(result>=0,"free-curve point has no actual region owner");
        return result;
    }
    Point physical(const Point& local) const
    { return {{sum(origin_[0],product(width_[0],local[0])),sum(origin_[1],product(width_[1],local[1]))}}; }
    Jet jet(const Point& p) const
    {
        return contextual([&] { return jetImpl(p); },[&] {
            return "[MappedInterface retained-coordinate jet, local="+describe(p)+"]";
        });
    }
    // Batch geometry consumes no correction columns and only the transverse
    // first derivative. All other consumers retain the original full jet.
    Jet geometryJet(const Point& p) const
    {
        return contextual([&] { return jetImpl<true>(p); },[&] {
            return "[MappedInterface retained-coordinate jet, local="+describe(p)+"]";
        });
    }
    Jet geometryJet(const Point& p,const Root& root) const
    {
        if constexpr(std::is_same_v<CellLayout,NativeQuarticCurveCell>) {
            if(geometryMeasure_&&guardedTransverseNewton_&&root.geometryPacket) {
                const auto& packet=*root.geometryPacket;
                // Only the final actual representative may consume its
                // query-private coefficients. All refined/mapped roots are
                // constructed without a packet and use the original jet.
                if(M::compare(p[1-chosen_],M::constant(packet.parameter))==0) {
                    FilteredTransverseJet bounded;
                    if(filteredTransverseGeometryJet(packet.controls,packet.parameter,root.derivativeSign,bounded)) {
                        Jet result;result.value=bounded.value;
                        if(chosen_==0) {result.dy=bounded.derivative;result.dyRoundoff=bounded.derivativeError;}
                        else {result.dx=bounded.derivative;result.dxRoundoff=bounded.derivativeError;}
                        return result;
                    }
                }
            }
        }
        return geometryJet(p);
    }
    template<bool GeometryOnly=false> Jet jetImpl(const Point& p) const
    {
        const Basis bx=stage("r Hermite basis",[&] { return Basis(p[0]); });
        const Basis by=stage("z Hermite basis",[&] { return Basis(p[1]); });
        Jet result;
        std::array<long double,3> sums{},compensations{},absoluteSums{};
        const auto accumulate=[&](int component,long double term) {
            const long double adjusted=term-compensations[component],next=sums[component]+adjusted;
            const long double absoluteNext=absoluteSums[component]+std::abs(term);
            require(std::isfinite(adjusted) && std::isfinite(next) && std::isfinite(absoluteNext),
                "nonfinite terminal Hermite jet accumulation");
            compensations[component]=(next-sums[component])-adjusted;
            sums[component]=next;
            absoluteSums[component]=absoluteNext;
        };
        // Keep the original SW/SE/NE/NW, value/dx/dy/dxy order for Cubic.
        // Quartic supplies all 25 tensor data in its explicit basis order.
        for (std::size_t k=0;k<dofCount;++k)
        {
            const auto index=CellLayout::tensorIndex(k);
            const int x=index[0],y=index[1];
            const double coefficient=dofs_[k];
            if (coefficient!=0.0)
            {
                accumulate(0,finiteProduct({bx.value[x],by.value[y],coefficient}));
                if constexpr(!GeometryOnly) {
                    accumulate(1,finiteProduct({bx.derivative[x],by.value[y],coefficient}));
                    accumulate(2,finiteProduct({bx.value[x],by.derivative[y],coefficient}));
                } else if(chosen_==0)
                    accumulate(2,finiteProduct({bx.value[x],by.derivative[y],coefficient}));
                else
                    accumulate(1,finiteProduct({bx.derivative[x],by.value[y],coefficient}));
            }
            if constexpr(!GeometryOnly)
                result.columns[k]=finiteProduct({bx.value[x],by.value[y],
                    (x==2 || x==3)?nominal_.x:1.0L,(y==2 || y==3)?nominal_.y:1.0L});
        }
        result.value=sums[0]; result.dx=sums[1]; result.dy=sums[2];
        // Conservative floating cancellation screen for coordinate terminal
        // conversion, basis products and the 16/25-term compensated
        // sum. Use binary64 epsilon even where long double is wider because
        // retained-expansion estimates are binary64. This is NOT a physical
        // tolerance or a certified jet enclosure; topology is supplied by
        // the independent exact-coefficient simple-root inventory.
        constexpr long double accumulatedRoundoff=8.0L*dofCount*std::numeric_limits<double>::epsilon();
        constexpr long double gamma=accumulatedRoundoff/(1.0L-accumulatedRoundoff);
        if constexpr(!GeometryOnly) {
            result.dxRoundoff=finiteProduct({gamma,absoluteSums[1]});
            result.dyRoundoff=finiteProduct({gamma,absoluteSums[2]});
        } else if(chosen_==0)
            result.dyRoundoff=finiteProduct({gamma,absoluteSums[2]});
        else
            result.dxRoundoff=finiteProduct({gamma,absoluteSums[1]});
        return result;
    }
    Root sectionRoot(const Ratio& q) const
    {
        Point a=point(graph_.localLower),b=point(graph_.localUpper);
        a[chosen_]=b[chosen_]=q;
        if constexpr(std::is_same_v<CellLayout,NativeQuarticCurveCell>) {
            if(geometryMeasure_&&guardedTransverseNewton_) {
                FilteredTransverseRoot filtered;
                std::array<Exact::Interval,5> controls;
                const bool qualified=stage("filtered native transverse root",[&] {
                    NativeWorkTimer timer(work_,&Public::NativeMomentWorkProfile::initialRootIsolation);
                    const auto coordinate=filteredIntervalDivide(Exact::enclosure(q.numerator),
                        Exact::enclosure(q.denominator));
                    using A=LevelSetHermiteScalar::Detail::Outward;
                    if(!A::finite(coordinate)||coordinate.lo<0.||coordinate.hi>1.)return false;
                    return layout_.transverseBernsteinControls(chosen_,coordinate,controls)
                        &&filteredTransverseRoot(controls,
                            chosen_==0?graph_.localLower.y:graph_.localLower.x,
                            chosen_==0?graph_.localUpper.y:graph_.localUpper.x,filtered);
                });
                if(qualified) {
                    if(work_)++work_->initialRootIsolations;
                    NativeWorkTimer timer(work_,&Public::NativeMomentWorkProfile::rootMappingRefinement);
                    if(work_)++work_->rootMappings;
                    // Already in the original native-cell transverse chart;
                    // retain each certified endpoint independently as a Ratio.
                    Root root{M::constant(filtered.parameter),M::constant(filtered.lower),
                        M::constant(filtered.upper),filtered.derivativeSign,
                        -filtered.derivativeSign,filtered.derivativeSign};
                    root.geometryPacket=FilteredGeometryPacket{controls,filtered.parameter};
                    return root;
                }
            }
        }
        const auto roots=stage("GL8 transverse roots",[&] { return isolate(a,b,false,guardedTransverseNewton_); });
        int selected=-1;
        for(int n=0;n<static_cast<int>(roots.size());++n)
            if(M::compare(roots[n].parameter,M::constant(0.0))>0 &&
                M::compare(roots[n].parameter,M::constant(1.0))<0)
            {
                require(selected<0,"transverse line has more than one interior free-graph root");
                selected=n;
            }
        require(selected>=0,"transverse line has no interior free-graph root");
        return mapRoot(roots[selected],a[1-chosen_],b[1-chosen_]);
    }
    Root refineIntegrationRoot(const Ratio& q,const Root& root) const
    {
        NativeWorkTimer timer(work_,&Public::NativeMomentWorkProfile::rootMappingRefinement,
            &Public::NativeMomentWorkProfile::restrictionConstruction);
        if(work_)++work_->rootRefinements;
        ++transverseRootRefinements_;
        const auto refined=MappedHermiteRootRefinement::refine(root.lower,root.upper,
            [&](const Ratio& lo,const Ratio& hi) {
                Point a,b;a[chosen_]=b[chosen_]=q;
                a[1-chosen_]=lo;b[1-chosen_]=hi;
                return polynomial(a,b);
            });
        require(refined.derivativeSign==root.derivativeSign,
            "transverse integration refinement changed its simple-root direction");
        return {refined.originalParameter,refined.originalLower,refined.originalUpper,refined.derivativeSign};
    }
    Root integrationRoot(const Ratio& q) const
    {
        auto root=sectionRoot(q);
        if (!graph_.controlledIntegration) return root;
        // A root EXACTLY at the interior basis knot cannot satisfy a
        // width/distance test: distance is zero even with an arbitrarily
        // narrow enclosure. Certify it by exact restriction of the same
        // scalar at that rational point, inside this isolated simple-root
        // enclosure. This changes neither an external root identity nor a
        // near-knot root; no floating snap or absolute error floor is used.
        const auto half=M::constant(.5);
        if (root.derivativeSign!=0 && M::compare(root.lower,half)<=0 && M::compare(half,root.upper)<=0) {
            Point at;at[chosen_]=q;at[1-chosen_]=half;
            if (Exact::sign(polynomial(at,at)[0])==0)
                return {half,half,half,root.derivativeSign};
        }
        // The Q4 basis has zeros at 0, 1/2 and 1. A binary64 root near
        // one of them can corrupt a small sensitivity column relatively,
        // even though the physical curve is perfectly regular. Refining
        // the quadrature interval then amplifies, rather than cures, that
        // sampling error. Resolve this internal transverse root in its
        // retained rational chart before converting the basis values.
        // This is below the 2e-10 quadrature relative tolerance, not a new
        // acceptance allowance. External face/endpoint roots are untouched.
        for (int refinement=0;;++refinement) {
            if (M::compare(root.lower,root.upper)==0) return root;
            const long double width=M::approximate(difference(root.upper,root.lower));
            require(width>0,"unrepresentable transverse root enclosure width");
            long double distance=1;
            for (double knot:{0.,.5,1.})
                distance=std::min(distance,std::abs(M::approximate(difference(root.parameter,M::constant(knot)))));
            if (width<=1e-13L*distance) return root;
            require(refinement<4,"transverse root could not resolve controlled Hermite basis sampling");
            root=refineIntegrationRoot(q,root);
        }
    }
    Root regionIntegrationRoot(const Ratio& q) const
    {
        auto root=integrationRoot(q);
        if(geometryMeasure_ && regions_.empty()) return root;
        // A binary64 transverse enclosure may straddle a clipping side even
        // when the SAME simple root lies strictly in one retained subarc.
        // Resolve that internal root; never assign an ambiguous point to an
        // owner, alter external face roots, or discard the thin subarc.
        for(int refinement=0;refinement<4;++refinement) {
            Point p,lo,hi;p[chosen_]=lo[chosen_]=hi[chosen_]=q;
            p[1-chosen_]=root.parameter;lo[1-chosen_]=root.lower;hi[1-chosen_]=root.upper;
            const int region=uniqueOwner(physical(p));
            if(region>=0 && uniqueOwner(physical(lo))==region && uniqueOwner(physical(hi))==region) return root;
            if(M::compare(root.lower,root.upper)==0) return root;
            root=refineIntegrationRoot(q,root);
        }
        // Original contextual owner/enclosure checks below remain mandatory
        // after the bounded work, and diagnose a genuinely unresolved case.
        return root;
    }
    std::vector<int> batchOwners(const Point& p,const std::vector<CellIntegrator*>& channels) const
    {
        std::vector<int> owners(channels.size(),-1);
        for(std::size_t t=0;t<channels.size();++t)
            if(channels[t]->geometryMeasure_ && channels[t]->regions_.empty())owners[t]=0;
        if(!batchOwnership_) {
            for(std::size_t t=0;t<channels.size();++t)owners[t]=channels[t]->uniqueOwner(p);
            return owners;
        }
        // Exact closed-rectangle membership is equivalent to the validated
        // convex rectangle's four signed orientation predicates. Evaluate it
        // once, then preserve each channel's SAME-owner boundary semantics.
        for(const auto& box:batchOwnership_->rectangles) {
            if(M::compare(p[0],box.lower[0])<0 || M::compare(p[0],box.upper[0])>0
                || M::compare(p[1],box.lower[1])<0 || M::compare(p[1],box.upper[1])>0)continue;
            for(std::size_t t=0;t<channels.size();++t) {
                if(channels[t]->geometryMeasure_ && channels[t]->regions_.empty())continue;
                if(owners[t]==-2)continue;
                if(owners[t]>=0 && owners[t]!=box.owners[t])owners[t]=-2;
                else owners[t]=box.owners[t];
            }
        }
        return owners;
    }
    // One query-private packet for the final actual root and the original
    // coordinate/ownership work. Never retained across nodes or queries.
    struct BatchRegionSample {
        Root root;
        Point local,location;
        std::vector<int> owners,lowerOwners,upperOwners;
    };
    BatchRegionSample batchRegionIntegrationRoot(const Ratio& q,const std::vector<CellIntegrator*>& channels,
        const std::vector<bool>& active) const
    {
        auto root=integrationRoot(q);
        const auto sample=[&](const Root& actual,bool terminal) {
            BatchRegionSample result;result.root=actual;
            Point lo,hi;result.local[chosen_]=lo[chosen_]=hi[chosen_]=q;
            result.local[1-chosen_]=actual.parameter;
            lo[1-chosen_]=actual.lower;hi[1-chosen_]=actual.upper;
            // Only after the fourth refinement did the old caller first
            // perform this physical conversion with its explicit stage.
            if(terminal)result.location=stage("GL8 batch shared physical point",[&] {return physical(result.local);});
            else result.location=physical(result.local);
            const auto lower=physical(lo),upper=physical(hi);
            result.owners=batchOwners(result.location,channels);
            result.lowerOwners=batchOwners(lower,channels);result.upperOwners=batchOwners(upper,channels);
            return result;
        };
        for(int refinement=0;refinement<4;++refinement) {
            auto measured=sample(root,false);
            bool resolved=true;
            for(std::size_t t=0;t<channels.size();++t)if(active[t])
                resolved=resolved && measured.owners[t]>=0 && measured.lowerOwners[t]==measured.owners[t]
                    && measured.upperOwners[t]==measured.owners[t];
            if(resolved || M::compare(root.lower,root.upper)==0)return measured;
            root=refineIntegrationRoot(q,root);
        }
        // The fourth refinement changed the root after the last ownership
        // attempt. Recompute from that final root before the mandatory checks.
        return sample(root,true);
    }
    using BatchArc=std::vector<Arc>;
    BatchArc integrateArcRuleBatch(const Cut& firstCut,const Cut& lastCut,
        const std::vector<CellIntegrator*>& channels,const std::vector<bool>& active) const
    {
        const auto& first=firstCut.local[chosen_];
        const auto& last=lastCut.local[chosen_];
        static constexpr std::array<double,8> x{{-0.96028985649753623168,-0.79666647741362673959,
            -0.52553240991632898582,-0.18343464249564980494,0.18343464249564980494,
            0.52553240991632898582,0.79666647741362673959,0.96028985649753623168}};
        static constexpr std::array<double,8> w{{0.10122853629037625915,0.22238103445337447054,
            0.31370664587788728734,0.36268378337836198297,0.36268378337836198297,
            0.31370664587788728734,0.22238103445337447054,0.10122853629037625915}};
        const auto extent=difference(last,first);
        require(Exact::sign(extent.numerator)>0,"nonpositive retained batch free-arc width");
        Arc prototype;
        prototype.nativeCellIndex=graph_.nativeCellIndex;
        prototype.requiresNativeDofConstraints=graph_.nativeCellIndex>=0;
        prototype.fromEnd=M::compare(sum(first,last),M::constant(1.0))>0;
        prototype.offset=M::approximate(prototype.fromEnd?difference(M::constant(1.0),last):first);
        prototype.width=M::approximate(extent);
        require(prototype.offset>=0.0L && prototype.offset<=1.0L && prototype.width>0.0L && prototype.width<=1.0L,
            "unrepresentable terminal batch free-arc offset or width");
        BatchArc output(channels.size(),prototype);
        std::vector<std::array<long double,2>> sums(channels.size()),compensations(channels.size());
        const long double hx=M::approximate(width_[0]),hy=M::approximate(width_[1]);
        constexpr long double twoPi=6.283185307179586476925286766559005768L;
        int derivativeSign=0;
        for(std::size_t k=0;k<x.size();++k)contextual([&] {
            const auto q=stage("GL8 batch chosen-coordinate interpolation",[&] {
                return MappedHermiteInterface::between(first,last,0.5*(1.0+x[k]));
            });
            require(M::compare(first,q)<0 && M::compare(q,last)<0,"GL8 batch coordinate is not strictly interior");
            const bool leftBlocked=M::compare(firstCut.upper,q)>=0;
            const bool rightBlocked=M::compare(q,lastCut.lower)>=0;
            if(leftBlocked || rightBlocked)
                requestResolution({leftBlocked?&firstCut:nullptr,rightBlocked?&lastCut:nullptr},
                    "GL8 batch coordinate="+describe(q)+", left root upper="+describe(firstCut.upper)
                    +", right root lower="+describe(lastCut.lower)+"; quadrature node is inside a cutting-root enclosure");
            const auto pointSample=batchRegionIntegrationRoot(q,channels,active);
            const auto& root=pointSample.root;const auto& p=pointSample.local;
            const auto& location=stage("GL8 batch shared physical point",[&]() -> const Point& {return pointSample.location;});
            const auto& regions=pointSample.owners;const auto& lowerOwners=pointSample.lowerOwners;
            const auto& upperOwners=pointSample.upperOwners;
            for(std::size_t t=0;t<channels.size();++t)if(active[t]) {
                const int region=regions[t];
                require(region!=-2,"batch free-curve point has different overlapping target owners");
                require(region>=0,"batch free-curve point has no actual target owner");
                require(lowerOwners[t]==region && upperOwners[t]==region,
                    "batch transverse root enclosure crosses an actual target boundary");
                require(output[t].regionId<0 || output[t].regionId==region,
                    "batch free subarc changes target owner without a retained face cut");
                output[t].regionId=region;
            }
            // One actual root and the consumed geometry jet per node,
            // never one per gauge or an unused correction-sensitivity query.
            const auto sample=geometryJet(p,root);
            const int sign=root.derivativeSign;
            const long double transverse=chosen_==0?sample.dy:sample.dx;
            const long double transverseRoundoff=chosen_==0?sample.dyRoundoff:sample.dxRoundoff;
            require(sign!=0 && std::isfinite(transverse) && std::abs(transverse)>transverseRoundoff
                && (transverse>0.0L?1:-1)==sign && (derivativeSign==0 || derivativeSign==sign),
                "batch free arc lost its regular transverse derivative");
            derivativeSign=sign;
            const long double radius=M::approximate(location[0]),z=M::approximate(location[1]);
            require(!axisymmetric_ || radius>0.0L,"RZ batch free-arc interior must have positive radius");
            const long double weight=finiteProduct({0.5L,w[k],prototype.width});
            for(std::size_t t=0;t<channels.size();++t)if(active[t]) {
                const auto gauge=channels[t]->greenOrigin_;
                const long double shiftedZ=gauge.y==0?z:M::approximate(difference(location[1],M::constant(gauge.y)));
                const long double shiftedRadius=gauge.x==0?radius:M::approximate(difference(location[0],M::constant(gauge.x)));
                const long double radialSum=gauge.x==0?radius:M::approximate(sum(location[0],M::constant(gauge.x)));
                const long double moment=chosen_!=0?0:axisymmetric_?finiteProduct({twoPi,radius,shiftedZ,hx,weight,
                    static_cast<long double>(sign)},1.0L):finiteProduct({shiftedZ,hx,weight,static_cast<long double>(sign)},1.0L);
                const long double radialMoment=chosen_!=1?0:finiteProduct({axisymmetric_?twoPi*.5L:1.L,shiftedRadius,
                    axisymmetric_?radialSum:1.L,hy,weight,static_cast<long double>(sign)},1.L);
                for(std::size_t c=0;c<2;++c) {
                    const long double term=c==0?moment:radialMoment;
                    require(std::isfinite(term),"nonfinite batch mapped free-arc integral");
                    const long double adjusted=term-compensations[t][c],next=sums[t][c]+adjusted;
                    require(std::isfinite(adjusted) && std::isfinite(next),"nonfinite batch mapped integral accumulation");
                    compensations[t][c]=(next-sums[t][c])-adjusted;sums[t][c]=next;
                }
                auto& result=output[t];
                result.minTransverseDerivative=k==0?std::abs(transverse):std::min(result.minTransverseDerivative,std::abs(transverse));
                result.maxScalarResidual=std::max(result.maxScalarResidual,std::abs(sample.value));
                ++result.quadraturePoints;
            }
        },[&] {return "[MappedInterface batch GL8 point "+std::to_string(k)+"]";});
        for(std::size_t t=0;t<channels.size();++t)if(active[t]) {
            output[t].axialMoment=sums[t][0];output[t].radialMoment=sums[t][1];
        }
        return output;
    }
    BatchArc integrateGeometryControlledBatch(const Cut& first,const Cut& last,BatchArc whole,
        std::size_t& evaluations,const std::vector<CellIntegrator*>& channels,const std::vector<bool>& active) const
    {
        struct Piece {Cut first,middle,last;BatchArc left,right;int depth;std::vector<long double> error;};
        const auto split=[&](const Cut& a,const Cut& b,const BatchArc& coarse,int depth) {
            require(evaluations+16<=16384,"batch geometry moment quadrature work budget exhausted");
            Cut middle;middle.local=a.local;
            middle.local[chosen_]=MappedHermiteInterface::between(a.local[chosen_],b.local[chosen_],.5);
            middle.lower=middle.upper=middle.local[chosen_];
            auto left=integrateArcRuleBatch(a,middle,channels,active);
            auto right=integrateArcRuleBatch(middle,b,channels,active);evaluations+=16;
            std::vector<long double> error(channels.size());
            for(std::size_t t=0;t<channels.size();++t)if(active[t])
                error[t]=std::max(std::abs(left[t].axialMoment+right[t].axialMoment-coarse[t].axialMoment),
                    std::abs(left[t].radialMoment+right[t].radialMoment-coarse[t].radialMoment));
            return Piece{a,middle,b,std::move(left),std::move(right),depth,std::move(error)};
        };
        const long double fraction=M::approximate(difference(last.local[chosen_],first.local[chosen_]));
        const long double capacity=M::approximate(width_[0])*M::approximate(width_[1])
            *(axisymmetric_?6.283185307179586476925286766559L*M::approximate(sum(origin_[0],width_[0])):1.L);
        const long double budget=2e-12L*capacity*fraction;
        std::vector<Piece> pieces;pieces.push_back(split(first,last,whole,0));
        for(;;) {
            std::vector<long double> error(channels.size());std::size_t worst=0;long double worstError=-1;
            for(std::size_t k=0;k<pieces.size();++k) {
                long double largest=0;
                for(std::size_t t=0;t<channels.size();++t)if(active[t]) {
                    error[t]+=pieces[k].error[t];largest=std::max(largest,pieces[k].error[t]);
                }
                if(largest>worstError){worst=k;worstError=largest;}
            }
            bool resolved=true;
            for(std::size_t t=0;t<channels.size();++t)if(active[t])resolved=resolved && error[t]<=budget;
            if(resolved) {
                for(std::size_t t=0;t<channels.size();++t)if(active[t]) {
                    auto& out=whole[t];out.axialMoment=out.radialMoment=0;out.quadraturePoints=0;
                    long double axialCompensation=0,radialCompensation=0;
                    for(const auto& piece:pieces)for(const auto* arc:{&piece.left[t],&piece.right[t]}) {
                        require(arc->regionId==out.regionId,"batch controlled graph quadrature changed target ownership");
                        const auto add=[](long double value,long double& sum,long double& compensation) {
                            const auto increment=value-compensation,next=sum+increment;
                            compensation=(next-sum)-increment;sum=next;
                        };
                        add(arc->axialMoment,out.axialMoment,axialCompensation);
                        add(arc->radialMoment,out.radialMoment,radialCompensation);
                        out.quadraturePoints+=arc->quadraturePoints;
                        out.minTransverseDerivative=std::min(out.minTransverseDerivative,arc->minTransverseDerivative);
                        out.maxScalarResidual=std::max(out.maxScalarResidual,arc->maxScalarResidual);
                    }
                    out.momentErrorEstimate=error[t];
                }
                return whole;
            }
            const auto piece=pieces[worst];
            if(piece.depth>=16)throw std::runtime_error(scalarFailureContext()+" batch geometry moment quadrature depth exhausted");
            auto left=split(piece.first,piece.middle,piece.left,piece.depth+1);
            auto right=split(piece.middle,piece.last,piece.right,piece.depth+1);
            pieces[worst]=std::move(left);pieces.push_back(std::move(right));
        }
    }
    BatchArc integrateArcBatch(const Cut& first,const Cut& last,const std::vector<CellIntegrator*>& channels,
        const std::vector<bool>& active) const
    {
        auto whole=integrateArcRuleBatch(first,last,channels,active);
        std::size_t evaluations=8;
        whole=integrateGeometryControlledBatch(first,last,std::move(whole),evaluations,channels,active);
        const auto midpoint=MappedHermiteInterface::between(first.local[chosen_],last.local[chosen_],.5);
        const int sign=integrationRoot(midpoint).derivativeSign;
        require(sign!=0,"batch region moment primitive lost its graph orientation");
        // Recover the second Green component with each ORIGINAL gauge in
        // exact retained endpoint coordinates, never a rounded gauge shift
        // of a previously accumulated universal moment.
        const auto a=physical(first.local),b=physical(last.local);
        const long double fraction=M::approximate(difference(last.local[chosen_],first.local[chosen_]));
        const long double capacity=M::approximate(width_[0])*M::approximate(width_[1])
            *(axisymmetric_?6.283185307179586476925286766559L*M::approximate(sum(origin_[0],width_[0])):1.L);
        for(std::size_t t=0;t<channels.size();++t)if(active[t]) {
            // Complete global graph recovery has the original whole-graph
            // endpoint primitive and error gate, after all shared subarcs.
            if(channels[t]->regions_.empty())continue;
            const auto gauge=channels[t]->greenOrigin_;
            const auto primitive=[&](const Point& p) {
                const auto shiftedY=difference(p[1],M::constant(gauge.y));
                const auto shiftedX=axisymmetric_?product(difference(p[0],M::constant(gauge.x)),
                    sum(p[0],M::constant(gauge.x))):difference(p[0],M::constant(gauge.x));
                return product(shiftedX,shiftedY);
            };
            const long double endpoint=sign*(axisymmetric_?3.1415926535897932384626433832795L:1.L)
                *M::approximate(difference(primitive(b),primitive(a)));
            auto& out=whole[t];
            if(chosen_==0)out.radialMoment=out.axialMoment-endpoint;
            else out.axialMoment=out.radialMoment-endpoint;
            require(std::isfinite(out.axialMoment) && std::isfinite(out.radialMoment)
                && out.momentErrorEstimate<=2e-12L*capacity*fraction
                    +2e-10L*std::min(std::abs(out.axialMoment),std::abs(out.radialMoment)),
                "recovered batch region Green moment exceeds controlled quadrature accuracy");
        }
        return whole;
    }
    Arc integrateArc(const Cut& firstCut,const Cut& lastCut) const
    {
        auto whole=integrateArcRule(firstCut,lastCut);
        if (graph_.controlledIntegration) {
            std::size_t evaluations=8;
            whole=geometryMeasure_?integrateGeometryControlled(firstCut,lastCut,std::move(whole),evaluations)
                :integrateControlled(firstCut,lastCut,std::move(whole),0,evaluations);
        }
        if (!geometryMeasure_) layout_.finalizeSensitivity(whole);
        else if(!regions_.empty()) {
            auto a=firstCut,b=lastCut;
            const auto midpoint=MappedHermiteInterface::between(a.local[chosen_],b.local[chosen_],.5);
            const int sign=integrationRoot(midpoint).derivativeSign;
            require(sign!=0,"region moment primitive lost its graph orientation");
            // Internal subdivision cuts have only a chosen coordinate. Obtain
            // the other endpoint from the same transverse root, never from a
            // terminal rounded phase range. Original external cuts retain both.
            const auto primitive=[&](const Cut& c) {
                const auto p=physical(c.local);
                const auto shiftedY=difference(p[1],M::constant(greenOrigin_.y));
                const auto shiftedX=axisymmetric_?product(difference(p[0],M::constant(greenOrigin_.x)),
                    sum(p[0],M::constant(greenOrigin_.x))):difference(p[0],M::constant(greenOrigin_.x));
                return product(shiftedX,shiftedY);
            };
            const long double endpoint=sign*(axisymmetric_?3.1415926535897932384626433832795L:1.L)
                *M::approximate(difference(primitive(b),primitive(a)));
            if(chosen_==0)whole.radialMoment=whole.axialMoment-endpoint;
            else whole.axialMoment=whole.radialMoment-endpoint;
            const long double fraction=M::approximate(difference(b.local[chosen_],a.local[chosen_]));
            const long double capacity=M::approximate(width_[0])*M::approximate(width_[1])
                *(axisymmetric_?6.283185307179586476925286766559L*M::approximate(sum(origin_[0],width_[0])):1.L);
            require(std::isfinite(whole.axialMoment) && std::isfinite(whole.radialMoment)
                && whole.momentErrorEstimate<=2e-12L*capacity*fraction
                    +2e-10L*std::min(std::abs(whole.axialMoment),std::abs(whole.radialMoment)),
                "recovered region Green moment exceeds controlled quadrature accuracy");
        }
        return whole;
    }
    Arc integrateGeometryControlled(const Cut& first,const Cut& last,Arc whole,std::size_t& evaluations) const
    {
        // Sum the error estimates over a complete graph interval. Assigning
        // a vanishing budget proportional to every tiny terminal interval
        // over-refines near-tangent endpoints even after their contribution
        // is far below the physical CV accuracy budget.
        struct Piece { Cut first,middle,last; Arc left,right; int depth; long double error; };
        const auto split=[&](const Cut& a,const Cut& b,const Arc& coarse,int depth) {
            require(evaluations+16<=16384,"geometry moment quadrature work budget exhausted");
            Cut middle;middle.local=a.local;
            middle.local[chosen_]=MappedHermiteInterface::between(a.local[chosen_],b.local[chosen_],.5);
            middle.lower=middle.upper=middle.local[chosen_];
            auto left=integrateArcRule(a,middle),right=integrateArcRule(middle,b);evaluations+=16;
            const long double error=std::max(std::abs(left.axialMoment+right.axialMoment-coarse.axialMoment),
                std::abs(left.radialMoment+right.radialMoment-coarse.radialMoment));
            return Piece{a,middle,b,std::move(left),std::move(right),depth,error};
        };
        const long double fraction=M::approximate(difference(last.local[chosen_],first.local[chosen_]));
        const long double capacity=M::approximate(width_[0])*M::approximate(width_[1])
            *(axisymmetric_?6.283185307179586476925286766559L*M::approximate(sum(origin_[0],width_[0])):1.L);
        const long double budget=2e-12L*capacity*fraction;
        std::vector<Piece> pieces;pieces.push_back(split(first,last,whole,0));
        for(;;) {
            long double error=0;std::size_t worst=0;
            for(std::size_t k=0;k<pieces.size();++k) {
                error+=pieces[k].error;
                if(pieces[k].error>pieces[worst].error) worst=k;
            }
            if(error<=budget) {
                whole.axialMoment=whole.radialMoment=0;whole.quadraturePoints=0;
                long double axialCompensation=0,radialCompensation=0;
                for(const auto& piece:pieces) for(const auto* arc:{&piece.left,&piece.right}) {
                    require(arc->regionId==whole.regionId,"geometry moment quadrature changed ownership");
                    const auto add=[](long double value,long double& sum,long double& compensation) {
                        const auto increment=value-compensation,next=sum+increment;
                        compensation=(next-sum)-increment;sum=next;
                    };
                    add(arc->axialMoment,whole.axialMoment,axialCompensation);
                    add(arc->radialMoment,whole.radialMoment,radialCompensation);
                    whole.quadraturePoints+=arc->quadraturePoints;
                    whole.minTransverseDerivative=std::min(whole.minTransverseDerivative,arc->minTransverseDerivative);
                    whole.maxScalarResidual=std::max(whole.maxScalarResidual,arc->maxScalarResidual);
                }
                whole.momentErrorEstimate=error;return whole;
            }
            const auto piece=pieces[worst];
            if(piece.depth>=16) throw std::runtime_error(scalarFailureContext()+" geometry moment quadrature depth exhausted");
            auto left=split(piece.first,piece.middle,piece.left,piece.depth+1);
            auto right=split(piece.middle,piece.last,piece.right,piece.depth+1);
            pieces[worst]=std::move(left);pieces.push_back(std::move(right));
        }
    }
    Arc integrateControlled(const Cut& first,const Cut& last,Arc whole,
        int depth,std::size_t& evaluations) const
    {
        require(evaluations+16<=16384,"controlled graph quadrature work budget exhausted");
        const auto q=MappedHermiteInterface::between(first.local[chosen_],last.local[chosen_],.5);
        Cut middle;middle.local=first.local;middle.local[chosen_]=q;middle.lower=middle.upper=q;
        auto left=integrateArcRule(first,middle),right=integrateArcRule(middle,last);
        evaluations+=16;
        const long double fraction=M::approximate(difference(last.local[chosen_],first.local[chosen_]));
        const long double hx=M::approximate(width_[0]),hy=M::approximate(width_[1]);
        const long double capacity=hx*hy*(axisymmetric_?
            6.283185307179586476925286766559L*M::approximate(sum(origin_[0],width_[0])):1.L);
        const auto tolerance=[&](long double a,long double b,long double scale) {
            return 2e-12L*scale*fraction+2e-10L*std::max(std::abs(a),std::abs(b));
        };
        const auto agrees=[&](long double a,long double b,long double scale) {
            return std::abs(a-b)<=tolerance(a,b,scale);
        };
        const long double momentError=std::max(std::abs(left.axialMoment+right.axialMoment-whole.axialMoment),
            std::abs(left.radialMoment+right.radialMoment-whole.radialMoment));
        bool resolved=agrees(whole.axialMoment,left.axialMoment+right.axialMoment,capacity)
            && agrees(whole.radialMoment,left.radialMoment+right.radialMoment,capacity);
        for (std::size_t k=0;!geometryMeasure_ && k<dofCount && resolved;++k)
            resolved=agrees(whole.sensitivity[k],left.sensitivity[k]+right.sensitivity[k],
                std::max(std::abs(whole.sensitivity[k]),std::abs(left.sensitivity[k])+std::abs(right.sensitivity[k])));
        if (!resolved) {
            if (depth>=16) {
                // Preserve the rejected local polynomial and the limiting
                // integral column so investigation needs no full-flow replay.
                long double worst=-1,coarse=0,fine=0,worstTolerance=0;
                int column=-3;
                const auto inspect=[&](int id,long double a,long double b,long double scale) {
                    const long double limit=tolerance(a,b,scale);
                    const long double error=std::abs(a-b);
                    const long double ratio=limit>0?error/limit:error==0?0:std::numeric_limits<long double>::infinity();
                    if (ratio>worst) {worst=ratio;column=id;coarse=a;fine=b;worstTolerance=limit;}
                };
                inspect(-2,whole.axialMoment,left.axialMoment+right.axialMoment,capacity);
                inspect(-1,whole.radialMoment,left.radialMoment+right.radialMoment,capacity);
                for (std::size_t k=0;!geometryMeasure_ && k<dofCount;++k)
                    inspect(static_cast<int>(k),whole.sensitivity[k],left.sensitivity[k]+right.sensitivity[k],
                        std::max(std::abs(whole.sensitivity[k]),std::abs(left.sensitivity[k])+std::abs(right.sensitivity[k])));
                std::ostringstream message;
                message<<std::setprecision(17)<<scalarFailureContext()
                    <<" controlled graph quadrature depth exhausted; depth="<<depth<<", evaluations="<<evaluations
                    <<", chosenInterval=("<<describe(first.local[chosen_])<<';'<<describe(last.local[chosen_])
                    <<"), worstColumn="<<column<<" (-2=axial,-1=radial,otherwise=localSensitivity)"
                    <<", coarse="<<coarse<<", fine="<<fine<<", tolerance="<<worstTolerance<<", errorRatio="<<worst;
                throw std::runtime_error(message.str());
            }
            left=integrateControlled(first,middle,std::move(left),depth+1,evaluations);
            right=integrateControlled(middle,last,std::move(right),depth+1,evaluations);
        }
        require(left.regionId==right.regionId,"controlled graph quadrature changed ownership");
        whole.axialMoment=left.axialMoment+right.axialMoment;
        whole.radialMoment=left.radialMoment+right.radialMoment;
        whole.momentErrorEstimate=resolved?momentError:left.momentErrorEstimate+right.momentErrorEstimate;
        for (std::size_t k=0;k<dofCount;++k) whole.sensitivity[k]=left.sensitivity[k]+right.sensitivity[k];
        whole.minTransverseDerivative=std::min(left.minTransverseDerivative,right.minTransverseDerivative);
        whole.maxScalarResidual=std::max(left.maxScalarResidual,right.maxScalarResidual);
        whole.quadraturePoints=left.quadraturePoints+right.quadraturePoints;
        return whole;
    }
    Arc integrateArcRule(const Cut& firstCut, const Cut& lastCut) const
    {
        const auto& first=firstCut.local[chosen_];
        const auto& last=lastCut.local[chosen_];
        static constexpr std::array<double,8> x{{-0.96028985649753623168,-0.79666647741362673959,
            -0.52553240991632898582,-0.18343464249564980494,0.18343464249564980494,
            0.52553240991632898582,0.79666647741362673959,0.96028985649753623168}};
        static constexpr std::array<double,8> w{{0.10122853629037625915,0.22238103445337447054,
            0.31370664587788728734,0.36268378337836198297,0.36268378337836198297,
            0.31370664587788728734,0.22238103445337447054,0.10122853629037625915}};
        const auto extent=difference(last,first);
        require(Exact::sign(extent.numerator)>0, "nonpositive retained free-arc width");
        Arc result;
        result.nativeCellIndex=graph_.nativeCellIndex;
        result.requiresNativeDofConstraints=graph_.nativeCellIndex>=0;
        result.fromEnd=M::compare(sum(first,last),M::constant(1.0))>0;
        result.offset=M::approximate(result.fromEnd?difference(M::constant(1.0),last):first);
        result.width=M::approximate(extent);
        require(result.offset>=0.0L && result.offset<=1.0L && result.width>0.0L && result.width<=1.0L,
            "unrepresentable terminal free-arc offset or width");
        const long double hx=M::approximate(width_[0]),hy=M::approximate(width_[1]);
        constexpr long double twoPi=6.283185307179586476925286766559005768L;
        std::array<long double,2+dofCount> sums{},compensations{};
        int derivativeSign=0;
        for (std::size_t k=0;k<x.size();++k)
        {
            contextual([&] {
            const auto q=stage("GL8 chosen-coordinate interpolation",[&] {
                return MappedHermiteInterface::between(first,last,0.5*(1.0+x[k]));
            });
            require(M::compare(first,q)<0 && M::compare(q,last)<0, "GL8 chosen coordinate is not strictly interior");
            // A representative is not the exact cutting root. A quadrature
            // node lying inside either cut's enclosure has no certified side
            // of that boundary, even when the representatives are distinct.
            // Refine the COMMON face roots before integrating such an arc;
            // never absorb it into a neighbor or silently erase its width.
            const bool leftBlocked=M::compare(firstCut.upper,q)>=0;
            const bool rightBlocked=M::compare(q,lastCut.lower)>=0;
            if(leftBlocked || rightBlocked)
                requestResolution({leftBlocked?&firstCut:nullptr,rightBlocked?&lastCut:nullptr},
                    "GL8 coordinate="+describe(q)+", left root upper="+describe(firstCut.upper)
                    +", right root lower="+describe(lastCut.lower)+"; quadrature node is inside a cutting-root enclosure");
            const auto root=regionIntegrationRoot(q);
            Point p; p[chosen_]=q; p[1-chosen_]=root.parameter;
            const auto location=stage("GL8 shared-bound physical point",[&] { return physical(p); });
            const int region=contextual([&] { return owner(location); },[&] {
                std::ostringstream message;
                message << std::setprecision(17) << "[MappedInterface GL8 point region ownership, cell=("
                    << layout_.address() << "), local=" << describe(p)
                    << ", physical=" << describe(location) << ", transverse lower=" << describe(root.lower)
                    << ", transverse upper=" << describe(root.upper) << ", regions=";
                for(const auto& candidate:regions_)
                {
                    message << '{' << candidate.regionId << ':';
                    for(const auto& vertex:candidate.vertices) message << '(' << vertex.x << ',' << vertex.y << ')';
                    message << '}';
                }
                return message.str()+']';
            });
            Point lo=p,hi=p;
            lo[1-chosen_]=root.lower; hi[1-chosen_]=root.upper;
            stage("GL8 root-enclosure ownership",[&] {
                require(owner(physical(lo))==region && owner(physical(hi))==region,
                    "transverse root enclosure crosses an actual region boundary");
            });
            require(result.regionId<0 || result.regionId==region, "free subarc changes owner without a retained face cut");
            result.regionId=region;
            const auto sample=jet(p);
            const int sign=root.derivativeSign;
            const long double transverse=chosen_==0?sample.dy:sample.dx;
            const long double transverseRoundoff=chosen_==0?sample.dyRoundoff:sample.dxRoundoff;
            require(sign!=0 && std::isfinite(transverse) && std::abs(transverse)>transverseRoundoff
                && (transverse>0.0L?1:-1)==sign
                && (derivativeSign==0 || derivativeSign==sign), "free arc lost its regular transverse derivative");
            derivativeSign=sign;
            const long double radius=M::approximate(location[0]),z=M::approximate(location[1]);
            // Translate the retained location, not its rounded global display.
            // Keep the zero-origin operations unchanged for all old consumers.
            const long double shiftedZ=greenOrigin_.y==0?z:
                M::approximate(difference(location[1],M::constant(greenOrigin_.y)));
            const long double shiftedRadius=greenOrigin_.x==0?radius:
                M::approximate(difference(location[0],M::constant(greenOrigin_.x)));
            const long double radialSum=greenOrigin_.x==0?radius:
                M::approximate(sum(location[0],M::constant(greenOrigin_.x)));
            require(!axisymmetric_ || radius>0.0L, "RZ free-arc interior must have positive radius");
            const long double weight=finiteProduct({0.5L,w[k],result.width});
            const long double moment=geometryMeasure_ && chosen_!=0?0:axisymmetric_ ? finiteProduct({twoPi,radius,shiftedZ,hx,weight,
                chosen_==0?static_cast<long double>(sign):sample.dy},
                chosen_==0?1.0L:std::abs(sample.dx)) : finiteProduct({shiftedZ,hx,weight,
                chosen_==0?static_cast<long double>(sign):sample.dy},
                chosen_==0?1.0L:std::abs(sample.dx));
            const long double radialMoment=geometryMeasure_ && chosen_!=1?0:finiteProduct({axisymmetric_?twoPi*.5L:1.L,shiftedRadius,
                axisymmetric_?radialSum:1.L,hy,weight,
                chosen_==1?static_cast<long double>(sign):sample.dx},
                chosen_==1?1.L:std::abs(sample.dy));
            for (std::size_t c=0;c<=dofCount+1;++c)
            {
                // nativeMoments exposes only Green moments. Do not compute
                // or impose accuracy on unrequested correction Jacobians.
                // The full region/sensitivity query keeps every column.
                if (geometryMeasure_ && c>0 && c<=dofCount) continue;
                const long double term=c==dofCount+1?radialMoment:c==0?moment:axisymmetric_ ? finiteProduct({-twoPi,radius,hx,hy,weight,
                    sample.columns[c-1]},std::abs(transverse)) : finiteProduct({-1.0L,hx,hy,weight,
                    sample.columns[c-1]},std::abs(transverse));
                require(std::isfinite(term), "nonfinite mapped free-arc integral");
                const long double adjusted=term-compensations[c],next=sums[c]+adjusted;
                require(std::isfinite(adjusted) && std::isfinite(next), "nonfinite mapped integral accumulation");
                compensations[c]=(next-sums[c])-adjusted; sums[c]=next;
            }
            result.minTransverseDerivative=k==0?std::abs(transverse):std::min(result.minTransverseDerivative,std::abs(transverse));
            result.maxScalarResidual=std::max(result.maxScalarResidual,std::abs(sample.value));
            ++result.quadraturePoints;
            },[&] { return "[MappedInterface GL8 point "+std::to_string(k)+"]"; });
        }
        result.axialMoment=sums[0];
        result.radialMoment=sums[dofCount+1];
        for (std::size_t c=0;c<dofCount;++c) result.sensitivity[c]=sums[c+1];
        return result;
    }
};
using Integrator = CellIntegrator<UniformCurveCell>;
using QuarticIntegrator = CellIntegrator<QuarticCurveCell>;
using NativeIntegrator = CellIntegrator<NativeCurveCell>;
using NativeQuarticIntegrator = CellIntegrator<NativeQuarticCurveCell>;
// C2 sources select their actual graph piece and complete source pullback.
// The integration embedding remains the original native-cell chart for all
// adapters; ordinary native cubic/Q4 instantiations keep their old layout.
template<class Scalar> using NativeIntegratorFor=CellIntegrator<
    typename NativeCurveCellSelection<Scalar>::Type>;
}
}
}
