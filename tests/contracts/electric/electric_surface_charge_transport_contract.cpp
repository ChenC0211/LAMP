#include "Modules/Electric/ElectricSurfaceChargeTransport.h"
#include "Modules/Electric/Detail/ElectricSurfaceChargeProjection.h"
#include "Core/Geometry/PreparedUnitPolynomial.h"
#include "Core/Geometry/NativeOutwardRounding.h"
#include "Modules/Interface/GeometricVofInitialization.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <algorithm>
#include <cmath>
#include <chrono>
#include <cerrno>
#include <cfenv>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>
#if defined(_M_X64) || (defined(__x86_64__) && defined(__SSE2__))
#include <xmmintrin.h>
#endif

namespace
{
using namespace LAMP::Electric;
using namespace LAMP::Interface;
using Clock=LAMP::TimeLayer;
unsigned assertions=0;
void check(bool yes,const char* message){++assertions;if(!yes)throw std::runtime_error(message);}
template<class F> void rejects(F&& f,const char* message)
{bool rejected=false;try{f();}catch(const std::exception&){rejected=true;}check(rejected,message);}
void projectionProfile(const ElectricSurfaceChargeTransferReport& r)
{
    check(r.coreHessianFallbacks+r.coreEndpointSignFallbacks+r.coreNonfiniteFallbacks+r.coreRootWidthFallbacks
        +r.coreForcedDistanceFallbacks==r.rootIsolationCalls,"Core fallback classifications do not cover actual completed calls");
    check(r.coreEndpointProjections+r.coreInteriorProjections==r.rootIsolationCalls,"Core result classifications lost endpoint/interior identity");
    check(r.geometryBasePreparations<=r.targetArcs&&r.geometryBasePreparations<=r.rootIsolationCalls,
        "geometry prefix was prepared beyond the actual target arc/call count");
    for(double seconds:{r.prepareSeconds,r.materialIntegrationSeconds,r.pointProjectionSeconds,
        r.nativeArcProjectionSeconds,r.coreArcProjectionSeconds,r.identityArcProjectionSeconds})
        check(std::isfinite(seconds)&&seconds>=0,"passive projection timer is invalid");
    std::cout<<"q1_profile prepare_seconds="<<r.prepareSeconds<<" integration_seconds="<<r.materialIntegrationSeconds
        <<" point_projection_seconds="<<r.pointProjectionSeconds<<" native_arc_seconds="<<r.nativeArcProjectionSeconds
        <<" core_arc_seconds="<<r.coreArcProjectionSeconds<<" identity_arc_seconds="<<r.identityArcProjectionSeconds
        <<" endpoint_requests="<<r.endpointProjectionRequests<<" endpoint_memo_hits="<<r.endpointProjectionMemoHits
        <<" core_hessian="<<r.coreHessianFallbacks<<" core_endpoint_sign="<<r.coreEndpointSignFallbacks
        <<" core_nonfinite="<<r.coreNonfiniteFallbacks<<" core_rootwidth="<<r.coreRootWidthFallbacks
        <<" core_forced_distance="<<r.coreForcedDistanceFallbacks<<" core_endpoints="<<r.coreEndpointProjections
        <<" core_interiors="<<r.coreInteriorProjections<<" geometry_bases="<<r.geometryBasePreparations<<'\n';
}
void preparedPolynomial()
{
    using namespace LAMP::Geometry;
    static_assert(!std::is_default_constructible_v<PreparedUnitPolynomial>);
    const auto sameDouble=[](double a,double b){return std::memcmp(&a,&b,sizeof(double))==0;};
    const auto sameRoots=[&](const UnitPolynomialRoots& a,const UnitPolynomialRoots& b){
        if(a.identicallyZero!=b.identicallyZero||a.roots.size()!=b.roots.size())return false;
        for(std::size_t i=0;i<a.roots.size();++i)if(!sameDouble(a.roots[i].point,b.roots[i].point)
            ||!sameDouble(a.roots[i].lower,b.roots[i].lower)||!sameDouble(a.roots[i].upper,b.roots[i].upper))return false;
        return true;};
    UnitPolynomialTerm endpoints;endpoints.factor=1;endpoints.coefficients[1]=-1;endpoints.coefficients[2]=1;
    UnitPolynomialTerm repeated;repeated.factor=1;repeated.coefficients[0]=.25;repeated.coefficients[1]=-1;repeated.coefficients[2]=1;
    UnitPolynomialTerm linear;linear.factor=1;linear.coefficients[0]=-.25;linear.coefficients[1]=1;
    UnitPolynomialTerm huge;huge.factor=std::ldexp(1.L,80);huge.coefficients[0]=1;
    auto negativeHuge=huge;negativeHuge.weight=-1;
    auto microscopic=linear;microscopic.coefficients[0]=-std::ldexp(1.,-50);
    const std::vector<std::vector<UnitPolynomialTerm>> cases{{endpoints},{repeated},{},{huge,linear,negativeHuge},{microscopic}};
    const std::vector<std::vector<double>> expected{{0,1},{},{},{.25},{std::ldexp(1.,-50)}};
    for(std::size_t c=0;c<cases.size();++c){const PreparedUnitPolynomial prepared(cases[c]);
        if(c==1){std::string originalFailure,preparedFailure;
            try{isolateUnitPolynomial(cases[c]);}catch(const std::exception& e){originalFailure=e.what();}
            try{prepared.isolate();}catch(const std::exception& e){preparedFailure=e.what();}
            check(!originalFailure.empty()&&originalFailure==preparedFailure
                &&originalFailure.find("[PolynomialSegmentRootsCandidate]")!=std::string::npos,
                "prepared polynomial changed the original unresolved tangential-root rejection");
        }else{const auto original=isolateUnitPolynomial(cases[c]),retained=prepared.isolate();
            check(sameRoots(original,retained),"prepared polynomial changed exact root points, brackets or endpoint inventory");
            check(retained.identicallyZero==(c==2)&&retained.roots.size()==expected[c].size(),"prepared polynomial lost a known zero/endpoint/cancellation root");
            for(std::size_t i=0;i<expected[c].size();++i)check(retained.roots[i].point==expected[c][i],"prepared polynomial changed a known exact root");
            const auto absolute=prepared.isolateNativeAbsolute();
            check(absolute.identicallyZero==retained.identicallyZero&&absolute.roots.size()==expected[c].size(),
                "absolute policy changed the complete ordinary/zero/endpoint root inventory");
            for(std::size_t i=0;i<expected[c].size();++i){const auto& root=absolute.roots[i];
                check(root.lower<=expected[c][i]&&root.upper>=expected[c][i],"absolute policy lost a known exact root");
                check(std::nextafter(root.upper-root.lower,std::numeric_limits<double>::infinity())<=64*std::numeric_limits<double>::epsilon(),
                    "absolute policy returned an enclosure wider than its declared native resolution");}}
        for(const auto range:std::array<std::array<double,2>,5>{{{{0,0}},{{1,1}},{{0,1}},{{.125,.875}},{{.25,.5}}}}){
            const auto a=encloseUnitPolynomial(cases[c],range[0],range[1]),b=prepared.enclose(range[0],range[1]);
            check(sameDouble(a[0],b[0])&&sameDouble(a[1],b[1]),"prepared polynomial changed an outward endpoint/range enclosure");}
        for(std::size_t split=0;split<=cases[c].size();++split){
            const std::vector<UnitPolynomialTerm> prefix(cases[c].begin(),cases[c].begin()+split),suffix(cases[c].begin()+split,cases[c].end());
            const PreparedUnitPolynomial base(prefix),appended(base,suffix);
            if(c==1){std::string originalFailure,appendedFailure;
                try{prepared.isolate();}catch(const std::exception& e){originalFailure=e.what();}
                try{appended.isolate();}catch(const std::exception& e){appendedFailure=e.what();}
                check(!originalFailure.empty()&&originalFailure==appendedFailure,"prefix append changed unresolved-root rejection");
            }else{check(sameRoots(prepared.isolate(),appended.isolate()),"prefix append changed exact root points or brackets");
                check(sameRoots(prepared.isolateNativeAbsolute(),appended.isolateNativeAbsolute()),
                    "prefix append changed absolute-policy root points or brackets");}
            for(const auto range:std::array<std::array<double,2>,5>{{{{0,0}},{{1,1}},{{0,1}},{{.125,.875}},{{.25,.5}}}}){
                const auto a=prepared.enclose(range[0],range[1]),b=appended.enclose(range[0],range[1]);
                check(sameDouble(a[0],b[0])&&sameDouble(a[1],b[1]),"prefix append changed outward bounds or cancellation sequence");}
            const auto before=encloseUnitPolynomial(prefix),after=base.enclose();
            check(sameDouble(before[0],after[0])&&sameDouble(before[1],after[1]),"appending changed immutable base coefficients");}}
    std::vector<UnitPolynomialTerm> mutableTerms{linear};const PreparedUnitPolynomial immutable(mutableTerms);const auto copy=immutable;
    mutableTerms[0].coefficients[0]=-.75;
    check(immutable.isolate().roots.front().point==.25&&copy.isolate().roots.front().point==.25
        &&isolateUnitPolynomial(mutableTerms).roots.front().point==.75,"prepared polynomial retained mutable caller coefficient storage");
    std::vector<UnitPolynomialTerm> mutablePrefix{huge},mutableSuffix{linear,negativeHuge};
    const PreparedUnitPolynomial frozenBase(mutablePrefix),frozenAppend(frozenBase,mutableSuffix);
    mutablePrefix[0].factor=1;mutableSuffix[0].coefficients[0]=-.75;
    check(frozenAppend.isolate().roots.front().point==.25,"prefix constructor retained mutable base/additional term storage");
    for(const auto range:std::array<std::array<double,2>,5>{{{{-.01,1}},{{0,1.01}},{{.75,.25}},{{std::numeric_limits<double>::quiet_NaN(),1}},{{0,std::numeric_limits<double>::infinity()}}}})
        rejects([&]{immutable.enclose(range[0],range[1]);},"prepared polynomial accepted an invalid unit range");
    auto invalid=linear;invalid.factor=std::numeric_limits<long double>::infinity();
    rejects([&]{PreparedUnitPolynomial bad({invalid});},"prepared polynomial accepted a nonfinite exact factor");
    PreparedUnitPolynomial moved({linear});const auto retainedMove=std::move(moved);
    check(retainedMove.isolate().roots.front().point==.25,"prepared polynomial move changed retained exact data");
    rejects([&]{moved.isolate();},"moved-from prepared polynomial was dereferenced");
    rejects([&]{moved.isolateNativeAbsolute();},"moved-from absolute-policy polynomial was dereferenced");
    rejects([&]{moved.enclose();},"moved-from prepared enclosure was dereferenced");
    rejects([&]{PreparedUnitPolynomial bad(moved,{});},"prefix constructor dereferenced a moved-from base");
    std::cout<<"prepared_unit_polynomial_contract PASS cases=5 roots_and_bounds_bitwise=1 original_multiple_rejection=1 immutable_input=1 invalid_range=1\n";
}
void preparedEnvironment()
{
    using namespace LAMP::Geometry;
    struct RestoreEnvironment
    {
        std::fenv_t saved;int savedErrno=errno;
#if defined(_M_X64) || (defined(__x86_64__) && defined(__SSE2__))
        unsigned savedMxcsr=_mm_getcsr();
#endif
        RestoreEnvironment(){std::fegetenv(&saved);}
        void restore() const
        {
            std::fesetenv(&saved);
#if defined(_M_X64) || (defined(__x86_64__) && defined(__SSE2__))
            _mm_setcsr(savedMxcsr);
#endif
            errno=savedErrno;
        }
        ~RestoreEnvironment(){restore();}
    } environment;
    UnitPolynomialTerm linear;linear.factor=1;linear.coefficients[0]=-.25;linear.coefficients[1]=1;
    const std::vector<UnitPolynomialTerm> terms{linear};const PreparedUnitPolynomial valid(terms);
    const auto rejectEveryEntry=[&]{
        rejects([&]{PreparedUnitPolynomial rejected(terms);},"prepared constructor accepted an invalid actual FP environment");
        rejects([&]{PreparedUnitPolynomial rejected(valid,terms);},"prepared prefix constructor accepted an invalid actual FP environment");
        rejects([&]{valid.isolate();},"prepared root query accepted an invalid actual FP environment");
        rejects([&]{valid.isolateNativeAbsolute();},"absolute root query accepted an invalid actual FP environment");
        rejects([&]{valid.enclose();},"prepared range query accepted an invalid actual FP environment");};
    for(int rounding:{FE_DOWNWARD,FE_UPWARD,FE_TOWARDZERO}){
        environment.restore();check(std::fesetround(rounding)==0,"cannot set requested negative-test rounding mode");
        rejectEveryEntry();check(std::fegetround()==rounding,"prepared rejection changed the caller's rounding mode");}
    unsigned underflowModes=0;
#if defined(_M_X64) || (defined(__x86_64__) && defined(__SSE2__))
    for(unsigned mode:{1u<<15,1u<<6}){
        environment.restore();_mm_setcsr(environment.savedMxcsr|mode);rejectEveryEntry();
        check((_mm_getcsr()&mode)!=0,"prepared rejection repaired the caller's FTZ/DAZ mode");++underflowModes;}
#endif
    environment.restore();
    const auto checkOutsidePolicy=[&]{
        check(std::fesetround(FE_DOWNWARD)==0,"cannot set rounding for prepared-scope leak check");
        rejects([&]{evaluateUnitPolynomial(terms,.5);},"prepared scope leaked permission to a later old free query");
        check(std::fegetround()==FE_DOWNWARD,"old free rejection repaired the caller's rounding mode");environment.restore();};
    {const PreparedUnitPolynomial success(terms);}checkOutsidePolicy();
    {const PreparedUnitPolynomial success(valid,{});}checkOutsidePolicy();
    valid.isolate();checkOutsidePolicy();valid.isolateNativeAbsolute();checkOutsidePolicy();valid.enclose();checkOutsidePolicy();
    auto invalid=linear;invalid.factor=std::numeric_limits<long double>::infinity();
    rejects([&]{PreparedUnitPolynomial bad({invalid});},"prepared constructor lost the nonfinite factor guard");checkOutsidePolicy();
    rejects([&]{PreparedUnitPolynomial bad(valid,{invalid});},"prepared prefix constructor lost the nonfinite factor guard");checkOutsidePolicy();
    UnitPolynomialTerm repeated;repeated.factor=1;repeated.coefficients[0]=.25;repeated.coefficients[1]=-1;repeated.coefficients[2]=1;
    const PreparedUnitPolynomial critical({repeated});
    rejects([&]{critical.isolate();},"prepared query lost the unresolved interior tangential-root guard");checkOutsidePolicy();
    rejects([&]{critical.isolateNativeAbsolute();},"absolute query lost the unresolved interior tangential-root guard");checkOutsidePolicy();
    UnitPolynomialTerm large;large.factor=std::numeric_limits<double>::max()/8;large.coefficients[11]=1;
    const PreparedUnitPolynomial overflowingRange({large});
    rejects([&]{overflowingRange.enclose();},"prepared query lost its original finite intermediate guard");checkOutsidePolicy();
    std::cout<<"prepared_environment_contract PASS rounding_modes=3 underflow_modes="<<underflowModes
        <<" every_entry_checked=1 success_and_exception_scope_restored=1 old_free_policy=1\n";
}
void nativeAbsolutePolynomial()
{
    using namespace LAMP::Geometry;
    // Exact root 1-2^-54 lies strictly between prev(1) and 1. Its rounded
    // representative need not itself be strictly inside the unit interval.
    UnitPolynomialTerm linear;linear.factor=1;linear.coefficients[0]=-1;linear.coefficients[1]=1;
    UnitPolynomialTerm offset;offset.factor=std::ldexp(1.L,-54);offset.coefficients[0]=1;
    const std::vector<UnitPolynomialTerm> terms{linear,offset};const PreparedUnitPolynomial prepared(terms);
    const auto previous=std::nextafter(1.,0.);
    check(prepared.enclose(previous,previous)[1]<0&&prepared.enclose(1,1)[0]>0,
        "near-one root does not have independently strict exact endpoint signs");
    std::string oldFree,oldPrepared;
    try{isolateUnitPolynomial(terms);}catch(const std::exception& e){oldFree=e.what();}
    try{prepared.isolate();}catch(const std::exception& e){oldPrepared=e.what();}
    check(oldFree==oldPrepared&&oldFree.find("endpoint-relative root resolution is unrepresentable")!=std::string::npos,
        "old free/prepared endpoint-relative resolution policy changed");
    const auto roots=prepared.isolateNativeAbsolute();
    check(!roots.identicallyZero&&roots.roots.size()==1,"absolute policy lost the unique near-one root");
    const auto& root=roots.roots.front();
    check(root.lower<root.upper&&root.lower<1&&root.upper>=previous
        &&root.point>=root.lower&&root.point<=root.upper,
        "near-one interior root was lost or collapsed to a rounded endpoint");
    check(std::nextafter(root.upper-root.lower,std::numeric_limits<double>::infinity())<=64*std::numeric_limits<double>::epsilon(),
        "near-one absolute root exceeded the original native-width gate");
    UnitPolynomialTerm twoRoots;twoRoots.factor=1;twoRoots.coefficients[0]=.1875;
    twoRoots.coefficients[1]=-1;twoRoots.coefficients[2]=1;
    const auto separated=PreparedUnitPolynomial({twoRoots}).isolateNativeAbsolute();
    check(separated.roots.size()==2&&!separated.identicallyZero,"absolute policy merged distinct ordinary roots");
    for(std::size_t i=0;i<2;++i){const auto& r=separated.roots[i];const double known=i?.75:.25;
        check(r.lower<=known&&known<=r.upper&&std::nextafter(r.upper-r.lower,std::numeric_limits<double>::infinity())
            <=64*std::numeric_limits<double>::epsilon(),"absolute multiple-root inventory lost containment or native resolution");}
    std::cout<<"native_absolute_polynomial_contract PASS old_relative_rejection_preserved=1 near_one_interior=1 complete_inventory=1\n";
}
template<class Real> constexpr bool nativeOutwardBinary64=sizeof(Real)==sizeof(std::uint64_t)
    &&std::numeric_limits<Real>::is_iec559&&std::numeric_limits<Real>::radix==2&&std::numeric_limits<Real>::digits==53;
template<class Real> Real nativeOutwardFromBits(std::uint64_t bits)
{Real value=0;std::memcpy(&value,&bits,std::min(sizeof bits,sizeof value));return value;}
template<class Real> bool nativeOutwardSameValue(Real a,Real b)
{
    if constexpr(nativeOutwardBinary64<Real>){std::uint64_t aa=0,bb=0;
        std::memcpy(&aa,&a,sizeof aa);std::memcpy(&bb,&b,sizeof bb);return aa==bb;}
    return (a==b&&std::signbit(a)==std::signbit(b))
        ||(std::isnan(a)&&std::isnan(b)&&std::signbit(a)==std::signbit(b));
}
template<class Real> struct NativeOutwardObservation {Real value=0;int error=0,flags=0;};
template<class Real,bool Shared> NativeOutwardObservation<Real> observeNativeOutward(Real input,bool upward,bool preset)
{
    std::feclearexcept(FE_ALL_EXCEPT);if(preset)std::feraiseexcept(FE_DIVBYZERO);errno=preset?EDOM:0;
    volatile Real runtimeInput=input;
    const auto direction=upward?std::numeric_limits<Real>::infinity():-std::numeric_limits<Real>::infinity();
    const auto value=Shared?LAMP::Geometry::nativeOutwardNeighbor(runtimeInput,upward):std::nextafter(runtimeInput,direction);
    const int error=errno,flags=std::fetestexcept(FE_ALL_EXCEPT);return {value,error,flags};
}
template<class Real> void nativeOutwardType(const char* name)
{
    const auto tiny=std::numeric_limits<Real>::denorm_min(),minimum=std::numeric_limits<Real>::min();
    const auto maximum=std::numeric_limits<Real>::max(),infinity=std::numeric_limits<Real>::infinity();
    std::vector<Real> fixed{Real(0),-Real(0),tiny,-tiny,minimum,-minimum,maximum,-maximum,
        std::nextafter(minimum,Real(0)),-std::nextafter(minimum,Real(0)),
        std::nextafter(maximum,Real(0)),-std::nextafter(maximum,Real(0)),
        Real(.5),Real(-.5),Real(1),Real(-1),std::ldexp(Real(1),-50),-std::ldexp(Real(1),-50),
        infinity,-infinity,std::numeric_limits<Real>::quiet_NaN(),-std::numeric_limits<Real>::quiet_NaN()};
    if constexpr(nativeOutwardBinary64<Real>)for(const std::uint64_t bits:{
        0x0010000000000001ull,0x001fffffffffffffull,0x0020000000000000ull,0x0020000000000001ull,
        0x7fdfffffffffffffull,0x7fe0000000000000ull,0x7fefffffffffffffull,
        0x8010000000000001ull,0x801fffffffffffffull,0x8020000000000000ull,
        0xffdfffffffffffffull,0xffe0000000000000ull,0xffefffffffffffffull,
        0x7ff8000000001234ull,0xfff8000000004321ull,0x7ff0000000000001ull,0xfff0000000000001ull})
        fixed.push_back(nativeOutwardFromBits<Real>(bits));
    std::vector<Real> normal;normal.reserve(4096);
    volatile std::uint64_t runtimeSeed=0x20261018a17e5eedull;std::uint64_t state=runtimeSeed;
    for(unsigned i=0;i<4096;++i){auto random=(state+=0x9e3779b97f4a7c15ull);
        random=(random^(random>>30))*0xbf58476d1ce4e5b9ull;
        random=(random^(random>>27))*0x94d049bb133111ebull;random^=random>>31;
        if constexpr(nativeOutwardBinary64<Real>){const auto exponent=2+(random%2044);
            normal.push_back(nativeOutwardFromBits<Real>((random&0x800fffffffffffffull)|(exponent<<52)));}
        else{const auto fraction=Real(.5)+static_cast<Real>(random&0xffffu)/Real(131072);
            normal.push_back(std::ldexp((random>>63)?-fraction:fraction,static_cast<int>(i%121)-60));}}
    std::uint64_t comparisons=0;
    for(const auto* inputs:{&fixed,&normal})for(std::size_t i=0;i<inputs->size();++i)
        for(bool upward:{false,true})for(bool preset:{false,true}){
            const auto expected=observeNativeOutward<Real,false>((*inputs)[i],upward,preset);
            const auto actual=observeNativeOutward<Real,true>((*inputs)[i],upward,preset);
            const bool same=nativeOutwardSameValue(actual.value,expected.value)
                &&actual.error==expected.error&&actual.flags==expected.flags;
            if(!same)std::cerr<<"native_outward_mismatch type="<<name<<" fixed="<<(inputs==&fixed)<<" index="<<i
                <<" upward="<<upward<<" preset="<<preset<<" std_errno="<<expected.error<<" shared_errno="<<actual.error
                <<" std_FE="<<expected.flags<<" shared_FE="<<actual.flags<<'\n';
            check(same,"shared native outward neighbor differs from same-type std::nextafter result/errno/FE");++comparisons;}
    std::cout<<"native_outward_contract PASS type="<<name<<" fixed_inputs="<<fixed.size()<<" normal_inputs="<<normal.size()
        <<" comparisons="<<comparisons<<" binary64_bits_checked="<<nativeOutwardBinary64<Real><<" errno_FE_equal=1\n";
}
void nativeOutward()
{
    struct RestoreEnvironment
    {
        std::fenv_t saved{};int savedErrno=errno;
#if defined(_M_X64) || (defined(__x86_64__) && defined(__SSE2__))
        unsigned savedMxcsr=_mm_getcsr();
#endif
        RestoreEnvironment(){check(std::fegetenv(&saved)==0,"cannot save outward-neighbor contract environment");}
        ~RestoreEnvironment(){std::fesetenv(&saved);
#if defined(_M_X64) || (defined(__x86_64__) && defined(__SSE2__))
            _mm_setcsr(savedMxcsr);
#endif
            errno=savedErrno;}
    } environment;
    nativeOutwardType<double>("double");nativeOutwardType<long double>("long_double");
}
void projectionCertificate()
{
    nativeOutward();
    preparedEnvironment();
    preparedPolynomial();
    nativeAbsolutePolynomial();
    namespace Projection=LAMP::Electric::Detail::SurfaceChargeProjection;
    using Point=Projection::Point;
    if constexpr(nativeOutwardBinary64<long double>){
        // Actual full07 arc42/cell689. Exact arithmetic independently places
        // its unique D' root in (prev(1),1), with strictly positive D''.
        const std::array<Point,5> nearOne{{
            {.5625,.67635605183793868},{.55497136507613476,.67982953728296858},
            {.54721403633392118,.68298272932516413},{.53928718942474696,.68542609742690752},
            {.53125,.68746301448570335}}};
        const auto result=Projection::select(nearOne,{.53125001032123809,.68746305521084983});
        check(result.selected&&result.interiorRoots==1&&result.uniqueInteriorCertificate
            &&result.minimum==Projection::MinimumCertificate::StrictInterior,
            "actual near-one minimum lost its complete-root StrictInterior certificate");
        const auto& selected=result.candidates[*result.selected];const auto& r=selected.parameter;
        check(selected.endpoint==-1&&r.lower<r.upper&&r.lower<1&&r.upper>=std::nextafter(1.,0.)
            &&r.point>=r.lower&&r.point<=r.upper,
            "actual nondegenerate interior root was mislabeled from a rounded representative");
        check(std::nextafter(static_cast<long double>(r.upper)-r.lower,std::numeric_limits<long double>::infinity())
            <=64*std::numeric_limits<long double>::epsilon(),"actual near-one root exceeded the original Q1 native-width gate");
    }
    // Exact binary64 receipt from prefix_probe_02: a normal-distance floor
    // hides the descent to a near-endpoint root in absolute D enclosures.
    const std::array<Point,5> actual{{
        {.34375,.40384792714177337},{.34982184840634817,.39585789939825888},
        {.35658876964184055,.38836377306059811},{.36393289058017481,.38139974197808157},
        {.37166277990834662,.375}}};
    const Point query{.34375002114156522,.40384794320784739};
    const auto fixed=Projection::select(actual,query);
    check(fixed.interiorRoots==1&&fixed.uniqueInteriorCertificate&&fixed.selected.has_value(),
        "actual near-endpoint minimum lacks the complete-root sign certificate");
    const auto& root=fixed.candidates[*fixed.selected];const auto& endpoint=fixed.candidates.front();
    check(fixed.roundedWinner==0&&*fixed.selected==2&&root.endpoint==-1&&root.parameter.point>0&&root.parameter.point<1,
        "actual projection kept the rounded endpoint winner instead of the certified interior root");
    check(root.parameter.lower<=4.6924705376556855e-15&&root.parameter.upper>=4.6924705376556855e-15,
        "actual projection changed the retained Core root enclosure");
    check(!(root.squaredDistance[1]<endpoint.squaredDistance[0])&&!(endpoint.squaredDistance[1]<root.squaredDistance[0]),
        "actual counterexample no longer exercises overlapping absolute distance bounds");
    check(fixed.minimum==Projection::MinimumCertificate::StrictInterior,"exact interior proof lost its typed certificate");
    // Actual full_02 neighboring arc11 shares owned node8 with arc12.
    // Its endpoint wins the rounded upper-distance ordering but is strictly
    // worse than arc12's certified interior minimum at the SAME query.
    const std::array<Point,5> previous{{
        {.34195866850117956,.40625},{.34240106100192286,.40564546275138164},
        {.34346446933818786,.40504356493652383},{.34329674888352396,.4044443664388474},
        {.34375,.40384792714177337}}};
    const auto neighbor=Projection::select(previous,query);
    check(neighbor.selected&&neighbor.minimum==Projection::MinimumCertificate::StrictEndpoint,
        "actual neighboring arc lacks its strict endpoint certificate");
    const auto& other=neighbor.candidates[*neighbor.selected];
    check(other.endpoint==1&&!(root.squaredDistance[1]<other.squaredDistance[0])
        &&!(other.squaredDistance[1]<root.squaredDistance[0])&&other.squaredDistance[1]<=root.squaredDistance[1],
        "actual cross-arc counterexample lost its rounded endpoint winner or overlapping bounds");
    const Projection::LocalMinimum interior{12,{8,26},{actual.front(),actual.back()},root.squaredDistance,root.endpoint,fixed.minimum};
    const Projection::LocalMinimum ending{11,{28,8},{previous.front(),previous.back()},other.squaredDistance,other.endpoint,neighbor.minimum};
    check(Projection::selectLocalMinimum({ending,interior})==std::optional<std::size_t>{1}
        &&Projection::selectLocalMinimum({interior,ending})==std::optional<std::size_t>{0},
        "shared actual endpoint dominance did not reselect the certified interior winner independently of candidate order");
    auto falseNode=ending;falseNode.nodes[1]=80;
    check(!Projection::selectLocalMinimum({falseNode,interior}),"coincident coordinates with a foreign node ID gained a dominance proof");
    auto falsePoint=ending;falsePoint.endpoints[1].x+=.01L;
    check(!Projection::selectLocalMinimum({falsePoint,interior}),"matching node ID with a different actual endpoint gained a dominance proof");
    auto noProof=interior;noProof.minimum=Projection::MinimumCertificate::None;
    check(!Projection::selectLocalMinimum({ending,noProof}),"interior label alone forged a strict global-minimum proof");
    auto secondInterior=ending;secondInterior.endpoint=-1;secondInterior.minimum=Projection::MinimumCertificate::StrictInterior;
    check(!Projection::selectLocalMinimum({secondInterior,interior}),"two overlapping interior minima were exempted by a shared node");
    check(!Projection::selectLocalMinimum({ending,interior,secondInterior}),"winner selection ignored another unresolved interior minimum");
    const std::array<Point,5> line{{{0,0},{.25,0},{.5,0},{.75,0},{1,0}}};
    for(unsigned side=0;side<2;++side){const auto result=Projection::select(line,{side?1.125:-.125,.25});
        check(result.interiorRoots==0&&!result.uniqueInteriorCertificate&&result.selected.has_value()
            &&result.candidates[*result.selected].endpoint==static_cast<int>(side)
            &&result.minimum==Projection::MinimumCertificate::StrictEndpoint,"true endpoint minimum changed");}
    // This symmetric graph has two equal interior minima separated by a
    // maximum. The complete inventory forbids the single-root certificate.
    const std::array<Point,5> symmetric{{{0,.25},{.25,0},{.5,-.125},{.75,0},{1,.25}}};
    const auto ambiguous=Projection::select(symmetric,{.5,.625});
    check(ambiguous.interiorRoots>=3&&!ambiguous.uniqueInteriorCertificate&&!ambiguous.selected&&ambiguous.overlap.has_value(),
        "multiple stationary points or equal minima were given a false uniqueness certificate");
    if constexpr(nativeOutwardBinary64<long double>)for(const auto& candidate:ambiguous.candidates)if(candidate.endpoint<0){
        const auto& r=candidate.parameter;
        check(std::nextafter(static_cast<long double>(r.upper)-r.lower,std::numeric_limits<long double>::infinity())
            <=64*std::numeric_limits<long double>::epsilon(),
            "MVT early return escaped the final native-absolute width gate");}
    const auto maximum=Projection::select(symmetric,{.5,1});
    check(maximum.interiorRoots==1&&!maximum.uniqueInteriorCertificate&&!maximum.selected,
        "a unique stationary maximum was mistaken for a unique minimum");
    std::cout<<"projection_certificate_contract PASS actual_counterexample=1 shared_endpoint_reselection=1 false_topology_and_unproved_rejected=1 endpoint_minima=2 multiple_stationary_ambiguity=1 unique_maximum_rejected=1 assertions="<<assertions<<'\n';
}
std::unique_ptr<QuadtreeMesh> meshFor(unsigned n)
{DictionaryMesh d;d.Nx=d.Ny=n;d.maxCellNum=4*n*n;d.maxInitialIterations=0;return std::make_unique<QuadtreeMesh>(d);}
VofSmoothGeometry smooth(const QuadtreeMesh& mesh,const GeometricVofPublication& source)
{VofSmoothGeometryConfig c;c.allowBoundaryEndpoints=true;return VofSmoothGeometry::prepare(mesh,source,c);}
std::vector<double> translationVelocity(const QuadtreeMesh&,const GeometricVofPublication&);
void identityReport(const ElectricSurfaceChargeTransferReport& r)
{
    check(r.mode==ElectricSurfaceChargeTransferMode::ExactPartitionIdentity&&!r.integralDiagnosticsAvailable,
        "direct identity did not declare its mode and unavailable integral diagnostics");
    check(r.identityArcTransfers==r.sourceArcs&&r.sourceArcs==r.targetArcs&&r.sourceDensityReconstructions==0,
        "direct identity reconstructed a density or failed to copy a complete partition");
    check(r.materialPanels==0&&r.materialSamples==0&&r.parameterSplits==0&&r.nativeArcProjections==0
        &&r.exactArcProjections==0&&r.localChainCertificates==0&&r.rootIsolationCalls==0
        &&r.endpointProjectionRequests==0&&r.endpointProjectionMemoHits==0&&r.geometryBasePreparations==0,
        "direct identity performed material sampling, projection or root work");
    check(r.materialIntegrationSeconds==0&&r.pointProjectionSeconds==0&&r.nativeArcProjectionSeconds==0
        &&r.coreArcProjectionSeconds==0&&r.identityArcProjectionSeconds==0,
        "direct identity entered the original integration/projection path");
    check(std::isnan(r.sourceAbsoluteCharge)&&std::isnan(r.sourceMomentQuadratureResidual),
        "direct identity forged an unmeasured absolute-charge or quadrature diagnostic");
    for(unsigned k=0;k<2;++k)check(std::isnan(r.sourcePositionMoment[k])&&std::isnan(r.advectedPositionMoment[k])
        &&std::isnan(r.projectedPositionMoment[k]),"direct identity forged an unmeasured spatial moment");
    check(r.chargeResidual==0&&r.relativeChargeResidual==0&&r.maximumDisplacement==0&&r.maximumProjectionDistance==0,
        "direct identity introduced motion, projection or a charge residual");
}
void sameMomentBits(const std::vector<std::array<double,3>>& a,const std::vector<std::array<double,3>>& b)
{
    check(a.size()==b.size(),"identity changed the moment partition");
    for(std::size_t i=0;i<a.size();++i)for(unsigned k=0;k<3;++k)
        check(std::memcmp(&a[i][k],&b[i][k],sizeof(double))==0,"identity did not preserve a signed moment bitwise");
}
void preserve(const ElectricVofSurfaceCharge& source,const ElectricSurfaceChargeTransport& result)
{
    check(result.sourceCharge().sameInstance(source),"receipt lost accepted charge owner");
    check(result.report().relativeChargeResidual<=1e-10L,"signed material charge did not conserve");
    if(result.report().integralDiagnosticsAvailable){
        check(result.report().sourceAbsoluteCharge>0,"signed cancellation erased material scale");
        check(result.report().sourceMomentQuadratureResidual<=1e-10L*result.report().sourceAbsoluteCharge,"source material moment quadrature failed");
    }else identityReport(result.report());
}
void stationary()
{
    auto mesh=meshFor(8);const auto vof=GeometricVofInitialization::prepare(*mesh,false,[](Vector2D p){return p.y-.423;});
    const auto geometry=smooth(*mesh,vof);
    const auto source=ElectricVofSurfaceCharge::fromDensity(*mesh,geometry,[&](std::size_t a,double t){
        const auto p=geometry.arcs()[a].sample(t);return static_cast<double>(2*p.x-1);});
    const auto snapshot=source.moments();
    const auto rebuilt=smooth(*mesh,vof);
    const auto remap=ElectricSurfaceChargeTransport::prepareSameEnd(*mesh,source,rebuilt);
    preserve(source,remap);check(!remap.interfaceTransport(),"same-End receipt forged VOF transport");
    check(remap.chargeCandidate().geometry().sameInstance(rebuilt),"identity did not bind the actual rebuilt target owner");
    sameMomentBits(snapshot,remap.chargeCandidate().moments());
    std::vector<std::array<double,3>> zeroNet(snapshot.size());
    for(std::size_t a=0;a<zeroNet.size();++a)zeroNet[a]={{a%2?-0.:0.,
        a%2?-std::numeric_limits<double>::denorm_min():std::numeric_limits<double>::denorm_min(),
        a%2?-.5:.5}};
    const auto higherModes=ElectricVofSurfaceCharge::fromMoments(*mesh,geometry,vof.timeLayer(),zeroNet);
    const auto higherIdentity=ElectricSurfaceChargeTransport::prepareSameEnd(*mesh,higherModes,rebuilt);
    preserve(higherModes,higherIdentity);sameMomentBits(zeroNet,higherIdentity.chargeCandidate().moments());
    check(higherIdentity.report().sourceCharge==0&&higherIdentity.report().targetCharge==0,
        "zero net/high-mode identity invented net charge");
    Clock end{1,0,.01,Clock::Position::End};
    std::vector<double> zero(static_cast<std::size_t>(vof.domain().nativeFaceCount()));
    const auto slab=GeometricVofAdvectorSlab::prepare(*mesh,vof,end,zero,zero);
    check(vof.geometryKind()==GeometricVofGeometryKind::LocalPlic
        &&slab.beginNormal()==zero&&slab.endNormal()==zero,
        "identity fixture lost the actual LocalPlic endpoint normal unknowns");
    rejects([&]{slab.beginVelocity();},"LocalPlic identity fixture unexpectedly exposes reconstructed begin vectors");
    rejects([&]{slab.endVelocity();},"LocalPlic identity fixture unexpectedly exposes reconstructed end vectors");
    VofSmoothGeometryConfig smoothConfig;smoothConfig.allowBoundaryEndpoints=true;
    const auto transport=GeometricVofTransport::prepareLocalSmoothSweep(*mesh,vof,geometry,slab,1,1,{},smoothConfig);
    for(const auto* stage:{&transport.midpoint(),&transport.endpoint()}){
        const auto& target=*(stage->timeLayer().position==Clock::Position::Midpoint?transport.smoothMidpoint():transport.smoothEndpoint());
        const auto attempt=ElectricSurfaceChargeTransport::prepareStage(*mesh,source,transport,target);
        const auto repeat=ElectricSurfaceChargeTransport::prepareStage(*mesh,source,transport,target);
        preserve(source,attempt);check(attempt.interfaceTransport()->sameInstance(transport),"receipt lost real VOF stage");
        sameMomentBits(snapshot,attempt.chargeCandidate().moments());
        check(attempt.chargeCandidate().geometry().sameInstance(target),"identity relabeled an old geometry as its actual target");
        const auto highStage=ElectricSurfaceChargeTransport::prepareStage(*mesh,higherModes,transport,target);
        preserve(higherModes,highStage);sameMomentBits(zeroNet,highStage.chargeCandidate().moments());
        check(attempt.chargeCandidate().timeLayer().position==stage->timeLayer().position,"candidate clock position was relabeled");
        check(attempt.chargeCandidate().timeLayer().begin==stage->timeLayer().begin&&attempt.chargeCandidate().timeLayer().end==stage->timeLayer().end,"candidate lost full interval");
        check(attempt.chargeCandidate().moments()==repeat.chargeCandidate().moments(),"repeat advanced an earlier candidate");
        check(attempt.report().materialSamples==repeat.report().materialSamples
            &&attempt.report().endpointProjectionRequests==repeat.report().endpointProjectionRequests
            &&attempt.report().endpointProjectionMemoHits==repeat.report().endpointProjectionMemoHits,
            "a new material receipt inherited another attempt's endpoint memo");
        rejects([&]{ElectricSurfaceChargeTransport::prepareStage(*mesh,attempt.chargeCandidate(),transport,target);},"candidate was accepted as the old material owner");
        rejects([&]{ElectricSurfaceChargeTransport::prepareSameEnd(*mesh,source,target);},"same-End accepted another clock");
        rejects([&]{ElectricSurfaceChargeTransport::prepareStage(*mesh,source,transport,smooth(*mesh,*stage));},
            "identity accepted a reconstructed target in place of the actual smooth stage owner");
        const auto foreignSmoothCharge=ElectricVofSurfaceCharge::fromMoments(*mesh,rebuilt,vof.timeLayer(),snapshot);
        rejects([&]{ElectricSurfaceChargeTransport::prepareStage(*mesh,foreignSmoothCharge,transport,target);},
            "identity accepted a rebuilt source instead of the actual smooth phase-flux owner");
    }
    const auto foreign=GeometricVofPublication::prepareLocalPlicInventory(*mesh,vof.domain(),{},[&]{std::vector<std::array<long double,2>> v;for(const auto& c:vof.plicCells())v.push_back(c.inventory);return v;}(),{});
    rejects([&]{ElectricSurfaceChargeTransport::prepareStage(*mesh,source,transport,smooth(*mesh,foreign));},"foreign target stage accepted");
    rejects([&]{ElectricSurfaceChargeTransport::prepareSameEnd(*mesh,source,smooth(*mesh,foreign));},"same-End accepted unrelated inventories with the same clock");
    const auto rz=GeometricVofInitialization::prepare(*mesh,true,[](Vector2D p){return p.y-.423;});
    const auto rzGeometry=smooth(*mesh,rz);const auto rzCharge=ElectricVofSurfaceCharge::fromDensity(*mesh,rzGeometry,[](std::size_t,double){return 1.;});
    rejects([&]{ElectricSurfaceChargeTransport::prepareSameEnd(*mesh,rzCharge,rzGeometry);},"RZ material map was silently claimed");
    const auto band=GeometricVofInitialization::prepare(*mesh,false,[](Vector2D p){return (p.y-.323)*(p.y-.681);});
    const auto bandGeometry=smooth(*mesh,band);const auto bandCharge=ElectricVofSurfaceCharge::fromDensity(*mesh,bandGeometry,[](std::size_t,double){return 1.;});
    rejects([&]{ElectricSurfaceChargeTransport::prepareSameEnd(*mesh,bandCharge,bandGeometry);},"separate components mixed material charge");
    // Same full-stage geometry and zero MEAN do not prove zero material
    // motion: both actual endpoint velocity snapshots are nonzero.
    const auto velocity=translationVelocity(*mesh,vof);auto reversed=velocity;
    for(auto& value:reversed)value=-value;
    const auto cancelledSlab=GeometricVofAdvectorSlab::prepare(*mesh,vof,end,velocity,reversed);
    const auto cancelled=GeometricVofTransport::prepareLocalPlicSweep(*mesh,vof,cancelledSlab,1,1);
    const auto cancelledTarget=smooth(*mesh,cancelled.endpoint());
    check(cancelledTarget.arcs().size()==geometry.arcs().size(),"cancelling-velocity negative changed the geometric partition");
    for(std::size_t a=0;a<geometry.arcs().size();++a)for(unsigned k=0;k<5;++k){
        const auto x=geometry.arcs()[a].controlPoints[k],y=cancelledTarget.arcs()[a].controlPoints[k];
        check(x.x==y.x&&x.y==y.y,"cancelling-velocity negative requires the same parameter geometry");}
    const auto cancellation=ElectricSurfaceChargeTransport::prepareStage(*mesh,source,cancelled,cancelledTarget);
    preserve(source,cancellation);
    check(cancellation.report().mode==ElectricSurfaceChargeTransferMode::MaterialProjection
        &&cancellation.report().integralDiagnosticsAvailable&&cancellation.report().materialSamples>0
        &&cancellation.report().sourceDensityReconstructions==source.moments().size(),
        "nonzero actual velocity was mistaken for identity because its mean/geometry was unchanged");
    check(source.moments()==snapshot,"rejected candidate modified accepted source");
    std::cout<<"stationary/same-End material samples="<<remap.report().materialSamples<<" signed_scale="<<static_cast<double>(remap.report().sourceAbsoluteCharge)
        <<" relative_Q="<<static_cast<double>(remap.report().relativeChargeResidual)<<'\n';
    projectionProfile(remap.report());
}
long double taper(long double x)
{
    if(x<=0||x>=1)return 0;if(x>=.2L&&x<=.8L)return 1;
    const auto z=x<.2L?x/.2L:(1-x)/.2L;return z*z*(3-2*z);
}
void localChainAndSeparateBranches(bool includeRemote=true)
{
    auto mesh=meshFor(16);const auto vof=GeometricVofInitialization::prepare(*mesh,false,[](Vector2D p){return std::hypot(p.x-.487,p.y-.513)-.18;});
    const auto geometry=smooth(*mesh,vof);const auto& arcs=geometry.arcs();
    std::vector<std::vector<std::size_t>> incident(geometry.nodes().size());
    for(std::size_t a=0;a<arcs.size();++a)for(const auto node:arcs[a].nodes)incident[node].push_back(a);
    std::size_t shortLinks=0;
    for(std::size_t a=0;a<arcs.size();++a){long double length=0;
        for(unsigned i=1;i<5;++i)length+=std::hypot(arcs[a].controlPoints[i].x-arcs[a].controlPoints[i-1].x,arcs[a].controlPoints[i].y-arcs[a].controlPoints[i-1].y);
        const auto& left=incident[arcs[a].nodes[0]];const auto& right=incident[arcs[a].nodes[1]];
        if(length>=.2L/16||left.size()!=2||right.size()!=2)continue;
        const auto before=left[0]==a?left[1]:left[0],after=right[0]==a?right[1]:right[0];bool adjacent=false;
        for(const auto x:arcs[before].nodes)for(const auto y:arcs[after].nodes)adjacent=adjacent||x==y;
        if(before!=after&&!adjacent)++shortLinks;}
    check(shortLinks>0,"closed regression contains no short arc between nonadjacent continuous neighbors");
    const auto source=ElectricVofSurfaceCharge::fromDensity(*mesh,geometry,[&](std::size_t a,double t){return static_cast<double>(.03L*(arcs[a].sample(t).x-.487L));});
    const auto snapshot=source.moments();const auto remap=ElectricSurfaceChargeTransport::prepareSameEnd(*mesh,source,smooth(*mesh,vof));
    preserve(source,remap);
    std::uint64_t digest=14695981039346656037ull;
    for(const auto& moment:remap.chargeCandidate().moments())for(double value:moment){std::uint64_t bits=0;std::memcpy(&bits,&value,sizeof(value));
        for(unsigned byte=0;byte<8;++byte){digest^=(bits>>(8*byte))&255u;digest*=1099511628211ull;}}
    long double maximumMomentDifference=0;
    for(std::size_t a=0;a<snapshot.size();++a)for(unsigned k=0;k<3;++k){
        const auto difference=std::abs(remap.chargeCandidate().moments()[a][k]-snapshot[a][k]);
        maximumMomentDifference=std::max(maximumMomentDifference,static_cast<long double>(difference));
        check(difference<2e-12,"local-chain same-End changed material moments");}
    std::cout<<"local continuous chain short_links="<<shortLinks<<" projection="<<static_cast<double>(remap.report().maximumProjectionDistance)
        <<" samples="<<remap.report().materialSamples<<" native_arcs="<<remap.report().nativeArcProjections
        <<" exact_arcs="<<remap.report().exactArcProjections<<" local_chains="<<remap.report().localChainCertificates
        <<" core_roots="<<remap.report().rootIsolationCalls<<" moment_digest="<<digest
        <<" maximum_moment_difference="<<maximumMomentDifference
        <<" position_moment_difference="<<std::hypot(remap.report().projectedPositionMoment[0]-remap.report().sourcePositionMoment[0],remap.report().projectedPositionMoment[1]-remap.report().sourcePositionMoment[1])<<'\n';
    projectionProfile(remap.report());
    if(!includeRemote)return;

    // A resolved broad lower cap connects two thin upper branches. Near the
    // upper neck the opposite branch is inside R, but its connecting curve is
    // far outside the actual local patch and the 2R arclength qualification.
    auto hairpinMesh=meshFor(32);const auto hairpin=GeometricVofInitialization::prepare(*hairpinMesh,false,[](Vector2D p){
        const auto z=std::max(0.,std::min(1.,(p.y-.4)/.3));const auto width=.002+.148*(1-z*z*(3-2*z));
        return std::hypot(p.x-.5,std::min(0.,p.y-.4))-width;});
    const auto hairpinGeometry=smooth(*hairpinMesh,hairpin);
    bool nearOppositeBranches=false;
    for(std::size_t a=0;a<hairpinGeometry.arcs().size();++a)for(std::size_t b=a+1;b<hairpinGeometry.arcs().size();++b){
        const auto x=hairpinGeometry.arcs()[a].sample(.5L),y=hairpinGeometry.arcs()[b].sample(.5L);
        const auto nx=hairpinGeometry.arcs()[a].normal(.5L),ny=hairpinGeometry.arcs()[b].normal(.5L);
        if(x.y>.8L&&y.y>.8L&&nx.x*ny.x+nx.y*ny.y<-.9L&&std::hypot(x.x-y.x,x.y-y.y)<.2L/32)nearOppositeBranches=true;}
    check(nearOppositeBranches,"remote-branch fixture has no actual opposite curves inside the declared radius");
    const auto hairpinCharge=ElectricVofSurfaceCharge::fromDensity(*hairpinMesh,hairpinGeometry,[](std::size_t,double){return 1.;});
    const auto hairpinSnapshot=hairpinCharge.moments();const auto acceptedOwner=hairpinCharge;
    // Exact same-End identity has an explicit parameter correspondence and
    // therefore needs no closest-point choice. Exercise the ORIGINAL remote
    // branch gate with nonzero actual velocity snapshots and zero full mean.
    const auto velocity=translationVelocity(*hairpinMesh,hairpin);auto opposite=velocity;
    for(auto& value:opposite)value=-value;
    const Clock end{1,0,.001,Clock::Position::End};
    const auto slab=GeometricVofAdvectorSlab::prepare(*hairpinMesh,hairpin,end,velocity,opposite);
    const auto transport=GeometricVofTransport::prepareLocalPlicSweep(*hairpinMesh,hairpin,slab,1,1);
    const auto target=smooth(*hairpinMesh,transport.endpoint());
    bool rejected=false;try{ElectricSurfaceChargeTransport::prepareStage(*hairpinMesh,hairpinCharge,transport,target);}
    catch(const ElectricSurfaceChargeTransferError& e){rejected=e.reason()==ElectricSurfaceChargeTransferFailure::AmbiguousProjection&&
        std::string(e.what()).find("material lies near a distinct nonadjacent surface branch")!=std::string::npos;}
    check(rejected,"nearby opposite branches were exempted by remote component connectivity");
    check(hairpinCharge.sameInstance(acceptedOwner)&&hairpinCharge.moments()==hairpinSnapshot,"rejected remote branch changed accepted material owner");
    std::cout<<"remote branch negative PASS arcs="<<hairpinGeometry.arcs().size()<<'\n';
}
std::vector<double> translationVelocity(const QuadtreeMesh& mesh,const GeometricVofPublication& vof)
{
    const auto psi=[](long double x,long double y){return .04L*y*taper(x)*taper(y);};
    std::vector<double> velocity(static_cast<std::size_t>(vof.domain().nativeFaceCount()));
    for(int row=0;row<vof.domain().nativeFaceCount();++row){const auto f=vof.domain().face(mesh,row);const auto length=f.area;
        const long double x=f.center.x,y=f.center.y,h=length/2;
        velocity[row]=static_cast<double>(f.horizontal?-(psi(x+h,y)-psi(x-h,y))/length:(psi(x,y+h)-psi(x,y-h))/length);}
    return velocity;
}
void coarseSupportRejects()
{
    auto mesh=meshFor(16);const auto vof=GeometricVofInitialization::prepare(*mesh,false,[](Vector2D p){return std::hypot(p.x-.487,p.y-.513)-.18;});
    const auto geometry=smooth(*mesh,vof);
    const auto source=ElectricVofSurfaceCharge::fromDensity(*mesh,geometry,[&](std::size_t a,double t){
        return static_cast<double>(1+2*(geometry.arcs()[a].sample(t).x-.487));});
    const auto accepted=source;const auto moments=source.moments();const auto identity=vof.identity();const auto clock=vof.timeLayer();
    std::vector<std::array<long double,2>> stocks;for(const auto& cell:vof.plicCells())stocks.push_back(cell.inventory);
    const auto velocity=translationVelocity(*mesh,vof);const Clock end{1,0,.025,Clock::Position::End};
    const auto slab=GeometricVofAdvectorSlab::prepare(*mesh,vof,end,velocity,velocity);
    const auto transport=GeometricVofTransport::prepareLocalPlicSweep(*mesh,vof,slab,1,1);
    bool rejected=false;try{smooth(*mesh,transport.endpoint());}
    catch(const VofSmoothGeometryError& e){rejected=e.reason()==VofSmoothGeometryFailure::InfeasibleVolume;}
    check(rejected,"known n16 curve-family support failure did not reject with its typed geometry reason");
    check(source.sameInstance(accepted)&&source.geometry().sameInstance(geometry)&&source.moments()==moments,
        "unsupported target changed accepted surface owner/moments");
    const auto now=vof.timeLayer();
    check(vof.identity()==identity&&now.step==clock.step&&now.begin==clock.begin&&now.end==clock.end&&now.position==clock.position,
        "unsupported target changed accepted VOF identity/full clock");
    for(std::size_t i=0;i<stocks.size();++i)check(vof.plicCells()[i].inventory==stocks[i],"unsupported target changed accepted independent stocks");
    std::cout<<"n16 coarse curve-family support typed rejection PASS; original positive motion remains unsupported\n";
}
void translation(unsigned n)
{
    auto mesh=meshFor(n);const auto vof=GeometricVofInitialization::prepare(*mesh,false,[](Vector2D p){return std::hypot(p.x-.487,p.y-.513)-.18;});
    const auto geometry=smooth(*mesh,vof);
    const auto source=ElectricVofSurfaceCharge::fromDensity(*mesh,geometry,[&](std::size_t a,double t){
        const auto p=geometry.arcs()[a].sample(t);return static_cast<double>(1+2*(p.x-.487));});
    const auto acceptedMoments=source.moments();
    long double previousError=1;
    for(double dt:{.025,.0125}){
        const auto velocity=translationVelocity(*mesh,vof);
        Clock end{1,0,dt,Clock::Position::End};const auto slab=GeometricVofAdvectorSlab::prepare(*mesh,vof,end,velocity,velocity);
        const auto transport=GeometricVofTransport::prepareLocalPlicSweep(*mesh,vof,slab,1,1);
        const auto target=smooth(*mesh,transport.endpoint());
        const auto result=ElectricSurfaceChargeTransport::prepareStage(*mesh,source,transport,target);preserve(source,result);
        check(result.report().maximumDisplacement>.039L*dt&&result.report().maximumDisplacement<.041L*dt,"real slab translation was not consumed");
        const auto& r=result.report();
        check(std::abs((r.advectedPositionMoment[0]-r.sourcePositionMoment[0])-.04L*dt*r.sourceCharge)<1e-11L*r.sourceAbsoluteCharge,"material charge did not follow the actual tangential/normal trajectory");
        check(std::abs(r.advectedPositionMoment[1]-r.sourcePositionMoment[1])<1e-11L*r.sourceAbsoluteCharge,"horizontal translation created vertical material drift");
        check(std::hypot(r.projectedPositionMoment[0]-r.advectedPositionMoment[0],r.projectedPositionMoment[1]-r.advectedPositionMoment[1])<=1.00001L*r.maximumProjectionDistance*r.sourceAbsoluteCharge,"projected charge location exceeds reported distance budget");
        check(result.report().maximumProjectionDistance<previousError,"projection error did not decrease on halving dt");previousError=result.report().maximumProjectionDistance;
        check(result.report().maximumProjectionDistance<.2L/n,"projection left local branch");
        check(std::isfinite(result.report().maximumProjectionDistance),"projection certificate returned a nonfinite distance");
        if(dt==.025){ElectricSurfaceChargeTransferConfig tight;tight.maximumProjectionDistanceInCells=(n/2.L)*result.report().maximumProjectionDistance;
            bool rejected=false;try{ElectricSurfaceChargeTransport::prepareStage(*mesh,source,transport,target,tight);}
            catch(const ElectricSurfaceChargeTransferError& e){rejected=e.reason()==ElectricSurfaceChargeTransferFailure::ProjectionDistance;}
            check(rejected,"a distance gate below the actual certified map was bypassed");
            check(source.moments()==acceptedMoments,"rejected projection changed accepted signed material charge");}
        std::cout<<"translation n="<<n<<" dt="<<dt<<" projection="<<static_cast<double>(result.report().maximumProjectionDistance)
            <<" relative_Q="<<static_cast<double>(result.report().relativeChargeResidual)<<" panels="<<result.report().materialPanels<<'\n';
    }
}
}
int main(int argc,char** argv)
{
    static_assert(!std::is_default_constructible_v<ElectricSurfaceChargeTransport>);
    try{
        if(argc==2&&std::string(argv[1])=="--native-outward-only"){nativeOutward();return 0;}
        if(argc==2&&std::string(argv[1])=="--projection-certificate-only"){projectionCertificate();return 0;}
        if(argc==2&&std::string(argv[1])=="--stationary-only"){stationary();return 0;}
        if(argc==2&&std::string(argv[1])=="--closed-projection-cost"){
            const auto start=std::chrono::steady_clock::now();localChainAndSeparateBranches(false);
            std::cout<<"closed_projection_contract_seconds="<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<'\n';return 0;}
        if(argc==2&&std::string(argv[1])=="--original-n16-motion"){translation(16);return 0;}
        if(argc==2&&std::string(argv[1])=="--resolved-n32-motion"){translation(32);return 0;}
        if(argc==2&&std::string(argv[1])=="--coarse-support-rejection"){coarseSupportRejects();return 0;}
        projectionCertificate();stationary();localChainAndSeparateBranches();coarseSupportRejects();translation(32);
        std::cout<<"electric_surface_charge_transport_contract PASS assertions="<<assertions<<'\n';return 0;}
    catch(const std::exception& e){std::cerr<<"electric_surface_charge_transport_contract FAIL "<<e.what()<<'\n';return 1;}
}


