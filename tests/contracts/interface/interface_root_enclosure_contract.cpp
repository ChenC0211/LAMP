#include "Core/Geometry/LevelSetMappedInterface.h"
#include "Core/Geometry/LevelSetInterfacePanels.h"
#include "Core/Geometry/LevelSetMappedSegment.h"
#include "Core/Geometry/LevelSetHermiteGeometry.h"
#include "Core/Geometry/Detail/ExactExpansion.h"
#include "Core/Geometry/Detail/LevelSetHermiteQuarticCell.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <algorithm>
#include <cmath>
#include <cfenv>
#include <cstdint>
#include <cstring>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace {
using I=LevelSetMappedInterface;
using Snapshot=LevelSetHermiteGeometry::Snapshot;
using Patch=LAMP::Geometry::ActiveQuarticHermitePatch;
constexpr long double pi=3.1415926535897932384626433832795L;
int cases=0,rejections=0;
void check(bool ok,const std::string& why){if(!ok)throw std::runtime_error(why);}
void near(long double a,long double b,long double scale,const std::string& why)
{if(!(std::isfinite(a) && std::abs(a-b)<=2e-10L*scale))
    std::cout<<"NUMERIC_FAILURE,"<<why<<','<<a<<','<<b<<','<<a-b<<','<<scale<<'\n';
    check(std::isfinite(a) && std::abs(a-b)<=2e-10L*scale,why);}
template<class F> void rejects(F f,const std::string& why){try{f();}catch(const std::exception& e){
    check(std::string(e.what()).find(why)!=std::string::npos,std::string("unexpected rejection: ")+e.what());++rejections;return;}
    throw std::runtime_error("missing rejection: "+why);}
// Binary64 exponent scaling must preserve EVERY expansion component bit,
// including the normal/subnormal boundaries and low representable tails.
void exactPowerOfTwoScalingContract()
{
    namespace E=LAMP::Geometry::Detail::ExactExpansion;
    const auto bits=[](double value){std::uint64_t result=0;std::memcpy(&result,&value,sizeof result);return result;};
    int scales=0;
    const auto compare=[&](const E::Expansion& value,int shift){
        const auto source=value.components();const auto scaled=E::scalePow2(value,shift);
        const auto actual=scaled.components();
        check(actual.size()==source.size(),"power-of-two scale changed exact component inventory");
        auto out=actual.begin();
        for(double component:source){
            const double reference=std::scalbn(component,shift);
            check(std::isfinite(reference) && reference!=0,"invalid bitwise scale reference");
            check(bits(*out++)==bits(reference),"power-of-two scale differs bitwise from original scalbn");
            check(bits(std::scalbn(reference,-shift))==bits(component),"power-of-two scale lost round-trip component bits");
        }
        ++scales;
    };
    const double minimum=std::numeric_limits<double>::min(),tiny=std::numeric_limits<double>::denorm_min();
    const std::vector<std::pair<double,int>> probes{
        {1.25,11},{0x1.23456789abcdep+17,-21},{std::numeric_limits<double>::max(),0},
        {minimum,-1},{1.5*minimum,-1},{std::nextafter(minimum,0.),0},
        {std::nextafter(minimum,0.),1},{tiny,52},{tiny,2097},{0x1p+1023,-2097}
    };
    for(const auto& probe:probes)for(double sign:{1.,-1.})compare(E::fromDouble(sign*probe.first),probe.second);
    const auto tail=E::add(E::fromDouble(1.),E::fromDouble(0x1p-60));
    check(tail.components().size()==2,"scaling fixture lost its independently representable tail");
    for(int shift:{-1014,0,1000})compare(tail,shift);
    // Expansion zero deliberately has no signed-zero distinction. Neither
    // extreme shift may attempt to negate INT_MIN or manufacture a component.
    for(double zero:{0.,-0.})for(int shift:{std::numeric_limits<int>::min(),std::numeric_limits<int>::max()})
        check(E::scalePow2(E::fromDouble(zero),shift).components().empty(),"power-of-two scale changed exact zero");
    const auto refused=[&](const E::Expansion& value,int shift){
        rejects([&]{E::scalePow2(value,shift);},"unsupported arithmetic scale");
    };
    refused(E::fromDouble(std::numeric_limits<double>::max()),1); // overflow
    refused(E::fromDouble(minimum),-53); // nonzero bit below denorm_min
    refused(E::fromDouble(std::nextafter(minimum,1.)),-1); // nonzero rounded result still loses its low bit
    refused(tail,-1015); // head remains normal, but the independent tail is unrepresentable
    for(int shift:{std::numeric_limits<int>::min(),std::numeric_limits<int>::max()})
        refused(E::fromDouble(1.),shift);
    std::cout<<"EXACT_POW2_SCALE,SCALBN_BITWISE,"<<scales<<",NORMAL_SUBNORMAL,TAIL,EXTREME_GATES,PASS\n";
}
struct Fixture {
    std::unique_ptr<QuadtreeMesh> mesh;
    std::shared_ptr<const Snapshot> source;
    std::vector<LevelSetMappedSegment> cuts;
    std::vector<I::Region> regions;
    I::Graph graph;
    const Vector2D lo{1,2},hi{2,3};
    Vector2D origin;
    bool xGraph,radial;
    int sign;
    double depth;
    Fixture(bool x,bool rz,int s,int exponent):xGraph(x),radial(rz),sign(s),depth(std::ldexp(1.,-exponent)) {
        DictionaryMesh d;d.Nx=d.Ny=1;d.minCoord=lo;d.maxCoord=hi;d.maxCellNum=64;d.maxInitialIterations=0;
        mesh=std::make_unique<QuadtreeMesh>(d);
        const auto jet=[&](Vector2D p){const double q=(x?p.x-lo.x:p.y-lo.y)-.5,v=x?p.y-lo.y:p.x-lo.x;
            const double g=.625+.125*q*q-depth,dg=.25*q;
            return Patch::NodalDof{sign*(v-g),sign*(x?-dg:1),sign*(x?1:-dg),0};};
        Patch p;p.activeCells={1};p.nodes.resize(mesh->vertIndexNum());p.cellCenters={jet(lo+Vector2D{.5,.5}).value};
        for(int v=0;v<mesh->vertIndexNum();++v)p.nodes[v]=jet(mesh->verts(mesh->vertIndexToLabel[v]).coord());
        const auto& cell=mesh->cells(mesh->cellIndexToLabel[0]);
        for(int k=0;k<4;++k){int a=mesh->verts(cell.cellToVert[k]).index(),b=mesh->verts(cell.cellToVert[(k+1)%4]).index();
            auto pa=mesh->verts(mesh->vertIndexToLabel[a]).coord(),pb=mesh->verts(mesh->vertIndexToLabel[b]).coord();
            const bool h=pa.y==pb.y;if(h?pa.x>pb.x:pa.y>pb.y){std::swap(a,b);std::swap(pa,pb);}
            const auto j=jet(pa+(pb-pa)*.5);p.edgeMidpoints.push_back({{{a,b}},j.value,h?j.dy:j.dx});}
        source=std::make_shared<const Snapshot>(Snapshot::fromNativeLeaves(*mesh,p,false));
        graph.nativeCellIndex=0;graph.coordinateIsX=x;graph.begin={x?0:1,0};graph.end={x?2:3,0};graph.controlledIntegration=true;
        if(x){const double y=lo.y+.625;origin={lo.x+.5,y};
            regions={{0,{{lo,{hi.x,lo.y},{hi.x,y},{lo.x,y}}}},
                {1,{{{lo.x,y},{hi.x,y},hi,{lo.x,hi.y}}}}};}
        else{const double r=lo.x+.625;origin={r,lo.y+.5};
            regions={{0,{{lo,{r,lo.y},{r,hi.y},{lo.x,hi.y}}}},
                {1,{{{r,lo.y},{hi.x,lo.y},hi,{r,hi.y}}}}};}
        for(const auto& region:regions)for(int k=0;k<4;++k){const auto a=region.vertices[k],b=region.vertices[(k+1)%4];
            const auto equal=[](Vector2D u,Vector2D v){return u.x==v.x && u.y==v.y;};bool found=false;
            for(const auto& cut:cuts)found=found || (equal(a,cut.begin()) && equal(b,cut.end())) || (equal(a,cut.end()) && equal(b,cut.begin()));
            if(!found)cuts.emplace_back(source,a,b);}
        // Resolve the EXTERNAL clipping roots through their public immutable
        // API before comparing a tiny arc width with an analytic value. Their
        // default representatives alone do not promise relative tiny widths.
        // The independent TRANSVERSE GL8 roots still exercise the production
        // ownership defect and its bounded same-root refinement.
        for(auto& cut:cuts)for(std::size_t k=0;k<cut.roots().size();++k)cut=cut.refined(k);
    }
    I::QuarticResult integrate(const std::vector<I::Region>& owner){
        for(int retry=0;;++retry){std::vector<const LevelSetMappedSegment*> refs;for(const auto& c:cuts)refs.push_back(&c);
            try{return I::integrateQuartic(*source,graph,refs,owner,radial,origin);}
            catch(const I::RootResolutionError& error){check(retry<256 && !error.requests.empty(),"retained root work limit");
                std::set<std::size_t> replaced;for(const auto& r:error.requests)if(replaced.insert(r.segmentIndex).second)
                    cuts.at(r.segmentIndex)=cuts.at(r.segmentIndex).refined(r.rootIndex);}}
    }
    long double capMoment()const{
        const long double e=depth,a=std::sqrt(e/.125L),area=4.L/3*e*a,square=16.L/15*e*e*a;
        return sign*(radial?(xGraph?-2*pi*origin.x*area:pi*(-2*origin.x*area+square)):-area);
    }
};
void nativeQuarticQueryGuard(const Fixture& f)
{
    const auto& scalar=*f.source;const auto p=f.origin;
    const auto chart=LAMP::Geometry::cellChartFromLocal(p-f.lo);
    const auto pack=[](const Snapshot::ScalarJet2& j){return std::array<double,6>{j.value,j.dx,j.dy,j.dxx,j.dxy,j.dyy};};
    const auto bits=[](double x){std::uint64_t value=0;std::memcpy(&value,&x,sizeof value);return value;};
    const auto original=pack(scalar.scalarJetNative(p));
    const auto same=[&](const Snapshot::ScalarJet2& j){const auto values=pack(j);
        for(std::size_t k=0;k<values.size();++k)check(bits(values[k])==bits(original[k]),"native Q4 cold/warm component changed");};
    same(scalar.scalarJetNative(p));same(scalar.scalarJetInNativeCell(0,p));
    check(bits(scalar.value(p))==bits(original[0]),"native exact value differs from actual scalar jet");
    // This unchanged public chart query supplies the original approximate
    // derivative channels; its value is not the exact physical-value port.
    const auto reference=pack(scalar.scalarJetAtNativeCellChart(0,chart));
    for(std::size_t k=1;k<reference.size();++k)
        check(bits(reference[k])==bits(original[k]),"native Q4 derivative differs from original chart query");
    const auto query=[&](int port){
        if(port==0)(void)scalar.scalarJetNative(p);
        else if(port==1)(void)scalar.scalarJetInNativeCell(0,p);
        else (void)scalar.value(p);
    };
    const auto refused=[](auto call,const char* reason){bool failed=false;
        try{call();}catch(const std::exception& error){failed=std::string(error.what()).find(reason)!=std::string::npos;}
        check(failed,std::string("actual native Q4 guard missing: ")+reason);};
    struct RestoreRound {int saved=std::fegetround();~RestoreRound(){std::fesetround(saved);}} restore;
    check(restore.saved==FE_TONEAREST,"native query guard fixture requires original nearest mode");
    check(std::fesetround(FE_DOWNWARD)==0,"cannot exercise native Q4 warm FP guard");
    // The unchanged chart entry uses the SAME geometry.cpp exact state,
    // so this rejects a permit leaked by a successful physical query.
    refused([&]{scalar.scalarJetAtNativeCellChart(0,chart);},"round-to-nearest");
    for(int port=0;port<3;++port)refused([&]{query(port);},"round-to-nearest");
    check(std::fesetround(FE_TONEAREST)==0,"cannot restore native Q4 nearest mode");
    bool invalid=false;try{scalar.scalarJetInNativeCell(-1,p);}catch(const std::exception&){invalid=true;}
    check(invalid,"native Q4 invalid-cell query did not fail inside prepared query");
    check(std::fesetround(FE_DOWNWARD)==0,"cannot exercise native Q4 exceptional FP guard");
    refused([&]{scalar.scalarJetAtNativeCellChart(0,chart);},"round-to-nearest");
    check(std::fesetround(FE_TONEAREST)==0,"cannot restore native Q4 exceptional mode");
#if defined(_M_X64) || (defined(__x86_64__) && defined(__SSE2__))
    {
        struct RestoreCsr {unsigned saved=_mm_getcsr();~RestoreCsr(){_mm_setcsr(saved);}} restoreCsr;
        _mm_setcsr(restoreCsr.saved|(1u<<15)|(1u<<6));
        for(int port=0;port<3;++port)refused([&]{query(port);},"gradual underflow");
    }
#endif
    same(scalar.scalarJetNative(p));same(scalar.scalarJetInNativeCell(0,p));
    check(bits(scalar.value(p))==bits(original[0]),"failed native Q4 queries changed exact value");
    std::cout<<"NATIVE_Q4_QUERY_GUARD,3,COMPONENTS,CORE_WARM,CORE_EXCEPTION,FP,FTZ,PASS\n";
}
void numericalQuarticTransportPoint(const Fixture& f)
{
    namespace E=LAMP::Geometry::Detail::ExactExpansion;
    using Q=LAMP::Geometry::Detail::LevelSetHermiteScalar::QuarticCell;
    // This existing unit-width fixture has exactly representable dyadic Q2
    // DOFs, including the tiny52 tail. The unchanged FULL restriction is an
    // independent exact reference; it never uses the new enclosing nets.
    const auto raw=f.source->nativeQuarticLocalDofs(0);Q::ExpandedControls data;
    for(int k=0;k<25;++k)data[k]=E::fromDouble(raw[k]);
    const std::array<Vector2D,3> points{{f.origin,{1.234567890123,2.789123456789},
        {std::nextafter(f.hi.x,f.lo.x),2.9375}}};
    int fallbacks=0;
    for(const auto p:points) {
        // Start with scalar-only cold preparation, then require the original
        // full jet to prepare its derivative nets independently on this leaf.
        const auto valueOnly=f.source->numericalScalarValueNative(p,1e-13);
        const auto valueWarm=f.source->numericalScalarValueNative(p,1e-13);
        check(valueOnly.boundsResolved&&valueOnly.absoluteErrorBound<=1e-13
            &&valueOnly.absoluteErrorBound==valueWarm.absoluteErrorBound
            &&valueOnly.exactValueFallback==valueWarm.exactValueFallback
            &&std::memcmp(&valueOnly.value,&valueWarm.value,sizeof(double))==0,
            "scalar-only Q4 cold/warm value certificate changed");
        const auto numeric=f.source->numericalScalarJetNative(p,1e-13,1e-11);
        check(std::memcmp(&valueOnly.value,&numeric.jet.value,sizeof(double))==0
            &&valueOnly.absoluteErrorBound==numeric.absoluteErrorBounds[0]
            &&valueOnly.exactValueFallback==numeric.exactValueFallback,
            "scalar-only Q4 changed the original value expression/certificate");
        const auto warm=f.source->numericalScalarJetNative(p,1e-13,1e-11);
        check(numeric.boundsResolved && numeric.absoluteErrorBounds==warm.absoluteErrorBounds,
            "numerical Q4 original-source bounds unresolved or changed warm");
        const std::array<double,3> values{{numeric.jet.value,numeric.jet.dx,numeric.jet.dy}};
        const auto chart=[](double p,double lower){const auto n=E::subtract(E::fromDouble(p),E::fromDouble(lower));
            return LAMP::Geometry::Detail::LevelSetHermiteScalar::RationalSegmentCoordinate{n,n,E::fromDouble(1)};};
        for(int k=0;k<3;++k) {
            const auto exact=Q::homogeneousRestriction(data,E::fromDouble(1),chart(p.x,f.lo.x),chart(p.y,f.lo.y),k==1?1:0,k==2?1:0);
            const auto residual=E::enclosure(E::subtract(exact.numerator[0],
                E::multiply(E::fromDouble(values[k]),exact.denominator)));
            const auto denominator=E::enclosure(exact.denominator);
            check(denominator.lo>0,"invalid independent point denominator");
            const long double error=std::max(std::abs(static_cast<long double>(residual.lo)),
                std::abs(static_cast<long double>(residual.hi)))/denominator.lo;
            check(error<=numeric.absoluteErrorBounds[k],"numerical Q4 bound missed original exact channel");
            if(k==0)check(error<=valueOnly.absoluteErrorBound,
                "scalar-only Q4 bound missed independent full exact restriction");
            check(numeric.absoluteErrorBounds[k]<=(k?1e-11:1e-13),"numerical Q4 channel exceeded independent budget");
        }
        check(numeric.jet.value==warm.jet.value && numeric.jet.dx==warm.jet.dx && numeric.jet.dy==warm.jet.dy,
            "numerical Q4 warm query changed first jet");
        if(numeric.exactValueFallback)++fallbacks;
    }
    check(fallbacks>0,"tiny52 near-zero point did not use exact value fallback");
    const auto tightValue=f.source->numericalScalarValueNative(points[1],1e-30);
    check(tightValue.boundsResolved&&tightValue.exactValueFallback,
        "over-budget scalar-only Q4 did not use original exact value restriction");
    const auto tight=f.source->numericalScalarJetNative(points[1],1e-30,1e-30);
    check(std::memcmp(&tightValue.value,&tight.jet.value,sizeof(double))==0
        &&tightValue.absoluteErrorBound==tight.absoluteErrorBounds[0],
        "scalar-only tight exact fallback changed original value");
    check(tight.boundsResolved&&tight.exactValueFallback&&tight.exactGradientFallback,
        "over-budget numerical Q4 channels did not use exact restrictions");
    const auto refused=[](auto query,const char* text){bool failed=false;try{query();}
        catch(const std::exception& e){failed=std::string(e.what()).find(text)!=std::string::npos;}
        check(failed,std::string("numerical Q4 missing rejection: ")+text);};
    refused([&]{f.source->numericalScalarJetNative(points[1],0,1e-11);},"point budget");
    for(double budget:{0.,-1.,-std::numeric_limits<double>::infinity(),std::numeric_limits<double>::quiet_NaN()})
        refused([&]{f.source->numericalScalarValueNative(points[1],budget);},"point budget");
    const auto invalidPoint=[&](Vector2D point){bool failed=false;
        try{f.source->numericalScalarValueNative(point,1e-13);}catch(const std::exception&){failed=true;}
        check(failed,"scalar-only Q4 accepted an invalid physical point");};
    invalidPoint({f.lo.x-1,f.lo.y});
    invalidPoint({std::numeric_limits<double>::quiet_NaN(),points[1].y});
    invalidPoint({points[1].x,std::numeric_limits<double>::infinity()});
    struct RestoreRound {int saved=std::fegetround();~RestoreRound(){std::fesetround(saved);}} restore;
    check(std::fesetround(FE_DOWNWARD)==0,"cannot set numerical Q4 FP guard fixture");
    refused([&]{f.source->numericalScalarJetNative(points[1]);},"round-to-nearest");
    refused([&]{f.source->numericalScalarValueNative(points[1]);},"round-to-nearest");
    check(std::fesetround(FE_TONEAREST)==0,"cannot restore numerical Q4 nearest mode");
#if defined(_M_X64) || (defined(__x86_64__) && defined(__SSE2__))
    {struct RestoreCsr {unsigned saved=_mm_getcsr();~RestoreCsr(){_mm_setcsr(saved);}} restoreCsr;
        _mm_setcsr(restoreCsr.saved|(1u<<15)|(1u<<6));
        refused([&]{f.source->numericalScalarJetNative(points[1]);},"gradual underflow");
        refused([&]{f.source->numericalScalarValueNative(points[1]);},"gradual underflow");}
#endif
    check(f.source->numericalScalarJetNative(points[1],1e-13,1e-11).boundsResolved,
        "failed numerical Q4 query changed the source or FP scope");
    const auto valueAfter=f.source->numericalScalarValueNative(points[1],1e-13);
    const auto fullAfter=f.source->numericalScalarJetNative(points[1],1e-13,1e-11);
    check(valueAfter.boundsResolved&&std::memcmp(&valueAfter.value,&fullAfter.jet.value,sizeof(double))==0
        &&valueAfter.absoluteErrorBound==fullAfter.absoluteErrorBounds[0],
        "failed scalar-only query changed source, original value or FP scope");
    std::cout<<"NUMERICAL_Q4_TRANSPORT_POINT,3,FULL_EXACT_VALUE_GRADIENT_BOUNDS,TINY52,REFLECTED,FP,FTZ,PASS\n";
    std::cout<<"NUMERICAL_Q4_TRANSPORT_VALUE,3,FULL_EXACT_VALUE_BOUND,COLD_VALUE_THEN_JET,TINY52,REFLECTED,TIGHT,INVALID,FP,FTZ,PASS\n";
}
void exactQuarticRepeat(const I::QuarticResult& a,const I::QuarticResult& b){
    check(a.sourceRootIncidences==b.sourceRootIncidences && a.retainedCuts==b.retainedCuts
        && a.arcs.size()==b.arcs.size(),"warmed integration changed root/cut/arc inventory");
    for(std::size_t k=0;k<a.arcs.size();++k){const auto& x=a.arcs[k];const auto& y=b.arcs[k];
        check(x.regionId==y.regionId && x.fromEnd==y.fromEnd && x.offset==y.offset && x.width==y.width
            && x.axialMoment==y.axialMoment && x.radialMoment==y.radialMoment
            && x.momentErrorEstimate==y.momentErrorEstimate && x.sensitivity==y.sensitivity
            && x.minTransverseDerivative==y.minTransverseDerivative && x.maxScalarResidual==y.maxScalarResidual
            && x.quadraturePoints==y.quadraturePoints && x.nativeCellIndex==y.nativeCellIndex
            && x.requiresNativeDofConstraints==y.requiresNativeDofConstraints
            && x.independentSensitivity.size()==y.independentSensitivity.size(),
            "warmed integration changed an actual arc result");
        for(std::size_t j=0;j<x.independentSensitivity.size();++j){const auto& p=x.independentSensitivity[j];const auto& q=y.independentSensitivity[j];
            check(p.kind==q.kind && p.index==q.index && p.component==q.component && p.coefficient==q.coefficient,
                "warmed integration changed independent source sensitivity");}
    }
}
void verify(Fixture& f){
    const auto split=f.integrate(f.regions);const auto expected=f.capMoment();
    const bool reuseWitness=f.xGraph && !f.radial && f.sign==1 && f.depth==std::ldexp(1.,-52);
    std::vector<std::vector<LevelSetMappedSegment::Root>> retainedRoots;
    if(reuseWitness){
        for(const auto& cut:f.cuts)retainedRoots.push_back(cut.roots());
        exactQuarticRepeat(split,f.integrate(f.regions));
    }
    check(split.arcs.size()==3,"two crossings must retain three nonempty arcs");
    long double axial=0,radial=0,capWidth=0,capAxial=0,capRadial=0,capAllowance=0;int capArcs=0;
    for(const auto& arc:split.arcs){axial+=arc.axialMoment;radial+=arc.radialMoment;
        if(arc.regionId==0){++capArcs;capWidth+=arc.width;
            // The unchanged transverse root solver initially has binary64
            // resolution. A tiny moment need not have fixed relative accuracy
            // before ownership demands more root resolution. Compare against
            // the independent polynomial integral at that floating scale;
            // this is a test allowance, not a certified production bound.
            const long double rootScale=16*std::numeric_limits<double>::epsilon()*arc.width*(f.radial?2*pi*3:1);
            const long double limit=2e-10L*std::abs(expected)+rootScale;
            check(std::abs(arc.axialMoment-expected)<=limit && std::abs(arc.radialMoment-expected)<=limit,"tiny cap analytic integral");
            capAxial=arc.axialMoment;capRadial=arc.radialMoment;capAllowance=limit;
            check(arc.axialMoment!=0 && arc.radialMoment!=0,"tiny cap erased");}}
    check(capArcs==1,"cap has no unique interval");
    const auto width=2*std::sqrt(static_cast<long double>(f.depth)/.125L);
    near(capWidth,width,width,"retained tiny arc width");
    auto joined=f.regions;for(auto& r:joined)r.regionId=0;
    const auto unsplit=f.integrate(joined);long double wholeA=0,wholeR=0;
    for(const auto& a:unsplit.arcs){wholeA+=a.axialMoment;wholeR+=a.radialMoment;}
    near(axial,wholeA,std::max(std::abs(wholeA),1.L),"axial partition additivity");
    near(radial,wholeR,std::max(std::abs(wholeR),1.L),"radial partition additivity");
    if(reuseWitness){
        // The same initial root certificates must not retain this changed owner context.
        for(const auto& arc:unsplit.arcs)check(arc.regionId==0,"warmed integration retained the previous region owner");
        exactQuarticRepeat(split,f.integrate(f.regions));
        bool fpRejected=false;
        {
            struct RestoreRound {int saved=std::fegetround();~RestoreRound(){std::fesetround(saved);}} restore;
            check(std::fesetround(FE_DOWNWARD)==0,"cannot exercise warmed integration FP guard");
            try{f.integrate(f.regions);}catch(const std::exception& error){
                fpRejected=std::string(error.what()).find("round-to-nearest")!=std::string::npos;}
        }
        check(fpRejected,"warmed integration ignored the actual FP environment");
        exactQuarticRepeat(split,f.integrate(f.regions));
        for(std::size_t k=0;k<f.cuts.size();++k){const auto& now=f.cuts[k].roots();const auto& before=retainedRoots[k];
            check(f.cuts[k].usesSource(*f.source) && now.size()==before.size(),"warmed context changed retained input roots/source");
            for(std::size_t j=0;j<now.size();++j){const auto& a=now[j];const auto& b=before[j];
                check(a.cellI==b.cellI && a.cellJ==b.cellJ && a.chart==b.chart && a.parameter==b.parameter
                    && a.lower==b.lower && a.upper==b.upper && a.pointOnly==b.pointOnly
                    && a.nativeCellIndex==b.nativeCellIndex,"warmed context shared a refined input root");}
        }
    }
    ++cases;
    std::cout<<"CAP,"<<f.xGraph<<','<<f.radial<<','<<f.sign<<','<<f.depth<<','<<expected<<','<<capWidth
        <<','<<capAxial<<','<<capRadial<<','<<std::max(std::abs(capAxial-expected),std::abs(capRadial-expected))<<','<<capAllowance<<'\n';
}
using Trace=LevelSetMappedSegment::NativeQuarticTrace;
std::shared_ptr<const Snapshot> traceSource(QuadtreeMesh& mesh,
    const std::function<Patch::NodalDof(Vector2D)>& f,bool rz=false,int selectedAxis=-1,double selectedLine=0)
{
    Patch p;p.activeCells.assign(mesh.cellIndexNum(),1);p.nodes.resize(mesh.vertIndexNum());p.cellCenters.resize(mesh.cellIndexNum());
    for(int v=0;v<mesh.vertIndexNum();++v)p.nodes[v]=f(mesh.verts(mesh.vertIndexToLabel[v]).coord());
    std::set<std::array<int,2>> edges;
    for(int c=0;c<mesh.cellIndexNum();++c){const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);
        const auto lo=mesh.verts(cell.cellToVert[0]).coord(),hi=mesh.verts(cell.cellToVert[2]).coord();
        if(selectedAxis>=0 && (selectedAxis==0?lo.x:lo.y)>=selectedLine)p.activeCells[c]=0;
        p.cellCenters[c]=f(lo+(hi-lo)*.5).value;
        if(!p.activeCells[c])continue;
        for(int k=0;k<4;++k){std::array<int,2> e{{mesh.verts(cell.cellToVert[k]).index(),mesh.verts(cell.cellToVert[(k+1)%4]).index()}};
            const auto a=mesh.verts(mesh.vertIndexToLabel[e[0]]).coord(),b=mesh.verts(mesh.vertIndexToLabel[e[1]]).coord();
            if(a.y==b.y?a.x>b.x:a.y>b.y)std::swap(e[0],e[1]);edges.insert(e);}}
    for(const auto& e:edges){const auto a=mesh.verts(mesh.vertIndexToLabel[e[0]]).coord(),b=mesh.verts(mesh.vertIndexToLabel[e[1]]).coord();
        const auto j=f(a+(b-a)*.5);p.edgeMidpoints.push_back({e,j.value,a.y==b.y?j.dy:j.dx});}
    return std::make_shared<const Snapshot>(Snapshot::fromNativeLeaves(mesh,p,rz));
}
auto traceProfile(int axis,double line,int order,int factor,int sign=1,double event=.375,double perturb=0,bool radial=false)
{
    return [=](Vector2D p){const double q=radial && axis==0?p.x*p.x-line*line:(axis==0?p.x:p.y)-line,t=axis==0?p.y:p.x;
        double g=1+t+t*t,dg=1+2*t;
        if(radial && axis==1){g=1+t*t;dg=2*t;}
        if(factor==1){g=t-event;dg=1;}
        if(factor==2){g=(t-event)*(t-event);dg=2*(t-event);}
        if(factor==3){g=t*t-event;dg=2*t;}
        const double power=std::pow(q,order),derivative=order*std::pow(q,order-1)*(radial && axis==0?2*p.x:1);
        const double normal=sign*derivative*g,tangent=sign*power*dg;
        return Patch::NodalDof{sign*power*g+perturb,axis==0?normal:tangent,axis==0?tangent:normal,sign*derivative*dg};};
}
void traceNear(long double actual,long double expected,long double scale,const std::string& why)
{check(std::isfinite(actual) && std::abs(actual-expected)<=1e-12L*scale,why);}
long double phaseLength(const std::vector<LevelSetMappedSegment::Range>& ranges)
{long double l=0;for(const auto& r:ranges)if(r.negative)l+=r.length;return l;}
void traceContract()
{
    int layouts=0,charts=0,events=0,refined=0,seams=0;
    for(bool rz:{false,true})for(bool amr:{false,true})for(int axis:{0,1}){
        DictionaryMesh d;d.Nx=d.Ny=4;d.maxCellNum=512;d.maxInitialIterations=0;
        QuadtreeMesh mesh(d);const double rlo=0;
        if(amr){DictionaryMesh::RefineBox b;b.level=1;b.min={rlo+.25,.25};b.max={rlo+.75,.75};mesh.applyRefineBoxes({b});}
        for(double offset: {amr?.125:.375,amr?.5:.25})for(int order:{1,2,3,4})for(int sign:{-1,1}){
            // (r^2-a^2)^k is a legal even-r Q4 for k<=2; axial
            // powers and Cartesian normal powers also cover k=3,4.
            if(rz && axis==0 && order>2)continue;
            const double line=offset+(axis==0?rlo:0);
            const auto source=traceSource(mesh,traceProfile(axis,line,order,0,sign,.375,0,rz),rz);
            const Vector2D a=axis==0?Vector2D{line,0}:Vector2D{rlo,line},b=axis==0?Vector2D{line,1}:Vector2D{rlo+1,line};
            auto q=LevelSetMappedSegment::nativeQuarticTrace(source,a,b,axis);
            auto reverse=LevelSetMappedSegment::nativeQuarticTrace(source,b,a,axis);
            check(q.usesSource(*source) && q.normalAxis()==axis && q.ordinaryRanges().empty(),"typed native zero source or range kind");
            check(q.coincidentRanges().size()==q.gridPieceCount() && q.coincidentRanges().size()==reverse.coincidentRanges().size(),"native zero chart coverage");
            long double total=0,ml=0,pl=0,rml=0,rpl=0,physical=0;
            for(const auto& z:q.coincidentRanges()){
                total+=z.interval.length;++charts;
                check(z.minus.available && z.plus.available,"internal zero side missing");
                for(const auto* s:{&z.minus,&z.plus}){
                    check(s->nativeCellIndex>=0 && source->nativeCellActive(s->nativeCellIndex) && s->normalOrder==order,"normal Taylor native provenance/order");
                    long double covered=0;for(const auto& r:s->ranges)covered+=r.length;
                    traceNear(covered,z.interval.length,z.interval.length,"side partition coverage");}
                ml+=phaseLength(z.minus.ranges);pl+=phaseLength(z.plus.ranges);
                const auto& r=z.interval;const long double mid=r.fromEnd?1-r.offset-r.length/2:r.offset+r.length/2;
                physical+=r.length*(rz?(axis==0?2*pi*line:2*pi*(rlo+mid)):1);
            }
            for(const auto& z:reverse.coincidentRanges()){rml+=phaseLength(z.minus.ranges);rpl+=phaseLength(z.plus.ranges);}
            traceNear(total,1,1,"zero trace measure closure");traceNear(physical,rz?(axis==0?2*pi*line:pi):1,physical,"XY/RZ physical measure");
            traceNear(ml,(sign<0)!=(order%2==1)?1:0,1,"minus normal Taylor phase");traceNear(pl,sign<0?1:0,1,"plus normal Taylor phase");
            traceNear(rml,ml,1,"reversed minus trace");traceNear(rpl,pl,1,"reversed plus trace");
            const auto middle=a+(b-a)*.375;
            auto left=LevelSetMappedSegment::nativeQuarticTrace(source,a,middle,axis),right=LevelSetMappedSegment::nativeQuarticTrace(source,middle,b,axis);
            long double subM=0,subP=0;
            for(const auto& z:left.coincidentRanges()){subM+=.375L*phaseLength(z.minus.ranges);subP+=.375L*phaseLength(z.plus.ranges);}
            for(const auto& z:right.coincidentRanges()){subM+=.625L*phaseLength(z.minus.ranges);subP+=.625L*phaseLength(z.plus.ranges);}
            traceNear(subM,ml,1,"actual child minus additivity");traceNear(subP,pl,1,"actual child plus additivity");
            rejects([&]{q.ordinary();},"no legacy phase view");
            rejects([&]{LevelSetMappedSegment old(source,a,b);},"positive-length coincident zero trace is unsupported");
            check(q.events().empty(),"zero chart invented finite value roots");++layouts;
        }
        if(amr){
            auto source=traceSource(mesh,traceProfile(axis,.25,1,0,1,.375,0,rz),rz);
            const Vector2D a=axis==0?Vector2D{.25,0}:Vector2D{0,.25},b=axis==0?Vector2D{.25,1}:Vector2D{1,.25};
            auto q=LevelSetMappedSegment::nativeQuarticTrace(source,a,b,axis);
            int localSeams=0;long double minus=0,plus=0;
            for(const auto& z:q.coincidentRanges()){
                check(z.minus.available && z.plus.available,"hanging side has no actual owner");
                const auto m=source->nativeCellBounds(z.minus.nativeCellIndex),p=source->nativeCellBounds(z.plus.nativeCellIndex);
                const double hm=m[1].x-m[0].x,hp=p[1].x-p[0].x;
                if(hm!=hp){check(hm==2*hp || hp==2*hm,"native hanging ratio is not 2:1");++localSeams;}
                minus+=phaseLength(z.minus.ranges);plus+=phaseLength(z.plus.ranges);
            }
            check(localSeams>0,"actual hanging subface was not exercised");seams+=localSeams;
            traceNear(minus,1,1,"hanging minus-side phase");traceNear(plus,0,1,"hanging plus-side phase");
        }
    }
    DictionaryMesh d;d.Nx=d.Ny=1;d.maxCellNum=64;d.maxInitialIterations=0;QuadtreeMesh mesh(d);
    for(int axis:{0,1})for(int factor:{1,2,3}){
        const double event=factor==3?.125:.375;
        auto source=traceSource(mesh,traceProfile(axis,.5,1,factor,1,event));
        Vector2D a=axis==0?Vector2D{.5,0}:Vector2D{0,.5},b=axis==0?Vector2D{.5,1}:Vector2D{1,.5};
        auto q=LevelSetMappedSegment::nativeQuarticTrace(source,a,b,axis);
        check(q.coincidentRanges().size()==1 && q.events().size()==2,"two-side normal root events");
        const double root=factor==3?std::sqrt(event):event;
        const auto expected=factor==2?0:root;
        traceNear(phaseLength(q.coincidentRanges()[0].plus.ranges),expected,1,"normal coefficient phase event");
        for(std::size_t k=0;k<q.events().size();++k){const auto& e=q.events()[k];++events;
            check(e.normalOrder==1 && std::abs(e.normalSide)==1 && e.root.nativeCellIndex==0 && !q.rootAtInputEndpoint(k,false) && !q.rootAtInputEndpoint(k,true),"typed normal event provenance");
            if(e.root.lower<e.root.upper){auto candidate=q.refined(k);++refined;
                check(candidate.events().size()==q.events().size(),"derivative refinement lost event");
                traceNear(phaseLength(candidate.coincidentRanges()[0].plus.ranges),expected,1,"same-order derivative refinement");}
        }
        rejects([&]{q.refined(q.events().size());},"index is out of range");
        traceNear(phaseLength(q.coincidentRanges()[0].plus.ranges),expected,1,"failed refinement changed original");
    }
    check(refined>0,"normal derivative refinement was not exercised");
    const double tail=std::ldexp(1.,-48);
    for(int axis:{0,1})for(bool reverse:{false,true}){
        auto source=traceSource(mesh,traceProfile(axis,.5,1,1,-1,1-tail));
        Vector2D a=axis==0?Vector2D{.5,0}:Vector2D{0,.5},b=axis==0?Vector2D{.5,1}:Vector2D{1,.5};if(reverse)std::swap(a,b);
        auto q=LevelSetMappedSegment::nativeQuarticTrace(source,a,b,axis);
        traceNear(phaseLength(q.coincidentRanges()[0].plus.ranges),tail,tail,"tiny normal phase erased or rounded");
        auto nonzero=traceSource(mesh,traceProfile(axis,.5,1,0,1,.375,tail));
        auto ordinary=LevelSetMappedSegment::nativeQuarticTrace(nonzero,a,b,axis);
        check(ordinary.coincidentRanges().empty() && ordinary.ordinaryRanges().size()==1 && !ordinary.ordinaryRanges()[0].negative,"near-zero nonzero classified zero");
        check(ordinary.ordinary().usesSource(*nonzero),"ordinary same-source view");
        Vector2D c=axis==0?Vector2D{0,.25}:Vector2D{.25,0},e=axis==0?Vector2D{1,.25}:Vector2D{.25,1};
        auto cross=LevelSetMappedSegment::nativeQuarticTrace(source,c,e,1-axis);LevelSetMappedSegment old(source,c,e);
        check(cross.ordinaryRanges().size()==old.ranges().size() && cross.events().size()==old.roots().size(),"ordinary legacy inventory changed");
        for(std::size_t k=0;k<old.ranges().size();++k){const auto& x=old.ranges()[k];const auto& y=cross.ordinaryRanges()[k];
            check(x.fromEnd==y.fromEnd && x.offset==y.offset && x.length==y.length && x.negative==y.negative,"ordinary ranges changed");}
    }
    for(double event:{0.,1.}){
        auto source=traceSource(mesh,traceProfile(0,.5,1,1,1,event));auto q=LevelSetMappedSegment::nativeQuarticTrace(source,{.5,0},{.5,1},0);
        for(std::size_t k=0;k<q.events().size();++k){check(q.rootAtInputEndpoint(k,event==1),"normal endpoint identity");
            rejects([&]{q.refined(k);},"exact point root needs no refinement");}
    }
    for(int order:{1,2,3,4}){
        const double line=order==1?1-tail:.5;
        auto scalar=traceSource(mesh,traceProfile(0,line,order,0));
        auto q=LevelSetMappedSegment::nativeQuarticTrace(scalar,{0,.25},{1,.25},1);
        LevelSetMappedSegment old(scalar,{0,.25},{1,.25});
        check(q.coincidentRanges().empty() && q.ordinaryRanges().size()==old.ranges().size(),"ordinary multiple-root view");
        long double gas=0;for(const auto& r:q.ordinaryRanges())if(!r.negative)gas+=r.length;
        traceNear(gas,order%2?1-line:1,order%2?1-line:1,"ordinary tiny or even phase");
        for(std::size_t k=0;k<old.ranges().size();++k){const auto& a=old.ranges()[k];const auto& b=q.ordinaryRanges()[k];
            check(a.fromEnd==b.fromEnd && a.offset==b.offset && a.length==b.length && a.negative==b.negative,"ordinary multiple ranges changed");}
    }
    // Independently rounded nodal samples of a near-endpoint squared profile
    // are not an exact squared Q4. The unchanged retained source can be
    // near-critical; both APIs must keep its bounded rejection.
    auto unresolved=traceSource(mesh,traceProfile(0,1-tail,2,0));
    rejects([&]{LevelSetMappedSegment old(unresolved,{0,.25},{1,.25});},"root-isolation budget exhausted");
    rejects([&]{LevelSetMappedSegment::nativeQuarticTrace(unresolved,{0,.25},{1,.25},1);},"root-isolation budget exhausted");
    DictionaryMesh md;md.Nx=md.Ny=2;md.maxCellNum=64;md.maxInitialIterations=0;QuadtreeMesh mixedMesh(md);
    auto mixedSource=traceSource(mixedMesh,[](Vector2D p){const double t=std::max(0.,p.y-.5);
        return Patch::NodalDof{p.x-.5+t*t,1,2*t,0};});
    auto mixed=LevelSetMappedSegment::nativeQuarticTrace(mixedSource,{.5,0},{.5,1},0);
    long double ordinaryWidth=0,zeroWidth=0;for(const auto& r:mixed.ordinaryRanges())ordinaryWidth+=r.length;
    for(const auto& z:mixed.coincidentRanges())zeroWidth+=z.interval.length;
    traceNear(ordinaryWidth,.5,1,"mixed ordinary piece");traceNear(zeroWidth,.5,1,"mixed zero piece");
    check(!mixed.events().empty(),"ordinary incidence at zero-piece endpoint lost");
    for(std::size_t k=0;k<mixed.events().size();++k)check(mixed.events()[k].normalOrder==0
        && !mixed.rootAtInputEndpoint(k,false) && !mixed.rootAtInputEndpoint(k,true),"mixed value root identity");
    auto source=traceSource(mesh,traceProfile(0,0,1,0));auto boundary=LevelSetMappedSegment::nativeQuarticTrace(source,{0,0},{0,1},0);
    check(!boundary.coincidentRanges()[0].minus.available && boundary.coincidentRanges()[0].plus.available,"missing physical outside side extrapolated");
    auto radial=traceSource(mesh,traceProfile(0,0,1,0,1,.375,0,true),true);
    auto axis=LevelSetMappedSegment::nativeQuarticTrace(radial,{0,0},{0,1},0);
    check(!axis.coincidentRanges()[0].minus.available && axis.coincidentRanges()[0].plus.normalOrder==2
        && phaseLength(axis.coincidentRanges()[0].plus.ranges)==0,"actual radial axis trace/normal order");
    auto zero=traceSource(mesh,[](Vector2D){return Patch::NodalDof{};});
    rejects([&]{LevelSetMappedSegment::nativeQuarticTrace(zero,{.5,0},{.5,1},0);},"all normal Taylor coefficients vanish");
    rejects([&]{LevelSetMappedSegment::nativeQuarticTrace(source,{0,0},{0,1},1);},"normal axis requires");
    rejects([&]{LevelSetMappedSegment::nativeQuarticTrace(source,{0,0},{1,1},0);},"normal axis requires");
    rejects([&]{LevelSetMappedSegment::nativeQuarticTrace(nullptr,{0,0},{0,1},0);},"null retained source");
    auto foreign=traceSource(mesh,traceProfile(0,0,1,0));check(!boundary.usesSource(*foreign),"wrong source accepted");
    std::weak_ptr<const Snapshot> lifetime=source;source.reset();check(!lifetime.expired(),"trace source lifetime lost");
    auto moved=std::move(boundary);rejects([&]{boundary.events();},"moved-from trace");check(moved.normalAxis()==0,"moved trace invalid");
    const int saved=std::fegetround();std::fesetround(FE_DOWNWARD);bool fpRejected=false;
    try{moved.normalAxis();}catch(const std::exception&){fpRejected=true;}std::fesetround(saved);check(fpRejected,"trace ignored FP environment");
    DictionaryMesh sd;sd.Nx=2;sd.Ny=1;sd.maxCellNum=64;sd.maxInitialIterations=0;QuadtreeMesh selected(sd);
    auto partial=traceSource(selected,traceProfile(0,.5,1,0),false,0,.5);
    auto cut=LevelSetMappedSegment::nativeQuarticTrace(partial,{.5,0},{.5,1},0);
    check(cut.coincidentRanges()[0].minus.available && !cut.coincidentRanges()[0].plus.available,"inactive source side extrapolated");
    rejects([&]{LevelSetMappedSegment::nativeQuarticTrace(partial,{.75,0},{.75,1},0);},"no active closed-cell coverage");
    std::cout<<"NATIVE_TRACE,"<<layouts<<','<<charts<<','<<events<<','<<refined<<','<<seams<<'\n';
}
using Panels=LevelSetInterfacePanels;
static_assert(!std::is_default_constructible_v<I::NativeSideSupport>);
void sideTopologyContract()
{
    int layouts=0,supports=0,hanging=0;
    for(bool rz:{false,true})for(bool amr:{false,true})for(int axis:{0,1})for(int sign:{-1,1}){
        DictionaryMesh d;d.Nx=d.Ny=4;d.maxCellNum=512;d.maxInitialIterations=0;QuadtreeMesh mesh(d);
        if(amr){DictionaryMesh::RefineBox b;b.level=1;b.min={.25,.25};b.max={.75,.75};mesh.applyRefineBoxes({b});}
        const double line=.25;auto source=traceSource(mesh,traceProfile(axis,line,1,0,sign,.375,0,rz),rz);
        const Vector2D a=axis==0?Vector2D{line,0}:Vector2D{0,line},b=axis==0?Vector2D{line,1}:Vector2D{1,line};
        const std::vector<Vector2D> contacts=rz && axis==1?std::vector<Vector2D>{b}:std::vector<Vector2D>{a,b};
        rejects([&]{Panels::discoverNativeComponents(mesh,source,rz,contacts);},"coincident zero trace");
        const auto topology=Panels::discoverNativeTraceComponents(mesh,source,rz,contacts);
        check(topology.curves.size()==1,"side topology component count");const auto& chain=topology.curves[0];
        check(!chain.closed() && chain.usesSource(*source),"side topology chain source");
        std::set<std::pair<bool,int>> ids;long double length=0,physical=0,radial=0,axial=0;
        for(int k=0;k<static_cast<int>(chain.graphs().size());++k){const auto& g=chain.graphs()[k];
            check(bool(g.nativeSide),"plane must use side support");const auto& s=*g.nativeSide;++supports;
            check(s.usesSource(*source) && s.ownerCell()==g.nativeCellIndex && s.normalAxis()==axis
                && g.begin.side==-1 && g.end.side==-1 && g.begin.rootOrdinal==-1 && g.end.rootOrdinal==-1,"side provenance/finite root invented");
            check(ids.insert({s.horizontal(),s.edgeLabel()}).second,"duplicate zero edge geometry");
            const auto& edge=s.horizontal()?mesh.hEdges(s.edgeLabel()):mesh.vEdges(s.edgeLabel());
            check((edge.edgeToVert[0]==s.beginVertex() && edge.edgeToVert[1]==s.endVertex())
                || (edge.edgeToVert[1]==s.beginVertex() && edge.edgeToVert[0]==s.endVertex()),"actual zero subedge vertices lost");
            check(s.minus().available() && s.plus().available() && s.minus().phase==-sign && s.plus().phase==sign,"two sided whole-leaf phase certificate");
            check(chain.zeroSetPartialResolved(k,axis,0,1),"side lost original complete normal certificate");
            const auto m=source->nativeCellBounds(s.minus().nativeCellIndex),p=source->nativeCellBounds(s.plus().nativeCellIndex);
            const double hm=m[1].x-m[0].x,hp=p[1].x-p[0].x;if(hm!=hp){check(hm==2*hp || hp==2*hm,"nonactual hanging ratio");++hanging;}
            for(double f:{0.,.375,1.}){const auto frame=chain.frameAt(k,f);
                traceNear(axis==0?frame.normal.x:frame.normal.y,sign,1,"side raw normal direction");
                traceNear(frame.levelSetResidual,0,1,"side raw scalar residual");
                check(frame.nativeCellIndex==g.nativeCellIndex && frame.hasCellChartPoint,"side frame chart identity");
                traceNear(chain.graphFraction(k,frame.cellChartPoint),f,1,"side chart round trip");
                check(chain.graphAt(g.nativeCellIndex,frame.cellChartPoint)>=0,"side frame graph lookup");}
            long double local=0,split=0;
            for(const auto& q:chain.quadrature(k)){local+=q.weight;physical+=q.physicalWeight;check(q.gradientMagnitude>0,"side original derivative lost");}
            for(const auto& q:chain.quadrature(k,0,.375))split+=q.weight;
            for(const auto& q:chain.quadrature(k,.375,1))split+=q.weight;
            traceNear(split,local,local,"side actual subarc additivity");length+=local;
            const auto green=I::nativeMoments(*source,g,rz);radial+=green.radial;axial+=green.axial;
            rejects([&]{I::integrateQuartic(*source,g,{}, {},rz);},"side-attached support");
            auto wrong=g;wrong.begin.rootOrdinal=0;rejects([&]{I::nativeMoments(*source,wrong,rz);},"public addressing");
            wrong=g;wrong.nativeCellIndex=s.ownerCell()==s.minus().nativeCellIndex?s.plus().nativeCellIndex:s.minus().nativeCellIndex;
            rejects([&]{I::nativeMoments(*source,wrong,rz);},"another scalar/cell");
        }
        traceNear(length,1,1,"side geometry full length");
        const long double measure=rz?(axis==0?2*pi*line:pi):1;
        traceNear(physical,measure,measure,"side XY/RZ physical measure");
        traceNear(radial,axis==0?sign*(rz?pi*line*line:line):0,1,"side closed radial Green reference");
        traceNear(axial,axis==1?sign*(rz?pi*line:line):0,1,"side closed axial Green reference");
        for(int c=0;c<mesh.cellIndexNum();++c){const auto bounds=source->nativeCellBounds(c);
            const double midpoint=axis==0?(bounds[0].x+bounds[1].x)/2:(bounds[0].y+bounds[1].y)/2;
            check(topology.cellPhase[c]==(midpoint<line?-sign:sign),"side full native leaf phase lost");}
        auto foreign=traceSource(mesh,traceProfile(axis,line,1,0,sign,.375,0,rz),rz);
        rejects([&]{I::nativeMoments(*foreign,chain.graphs()[0],rz);},"different snapshot");
        std::weak_ptr<const Snapshot> life=source;source.reset();check(!life.expired(),"side graph source lifetime lost");++layouts;
    }
    check(hanging>0,"native topology did not exercise actual hanging sides");
    DictionaryMesh d;d.Nx=d.Ny=4;d.maxCellNum=512;d.maxInitialIterations=0;QuadtreeMesh mesh(d);
    const auto mixed=traceSource(mesh,[](Vector2D p){const double t=std::max(0.,p.y-.5);return Patch::NodalDof{p.x-.5+t*t,1,2*t,0};});
    const auto join=Panels::discoverNativeTraceComponents(mesh,mixed,false,{{.5,0},{.25,1}});
    check(join.curves.size()==1,"ordinary/side join changed component count");int side=0,ordinary=0;
    for(int k=0;k<static_cast<int>(join.curves[0].graphs().size());++k){const auto& g=join.curves[0].graphs()[k];
        if(g.nativeSide)++side;else ++ordinary;join.curves[0].frameAt(k,.5);I::nativeMoments(*mixed,g,false);}
    check(side>0 && ordinary>0,"ordinary/side join not exercised");
    const auto twice=traceSource(mesh,[](Vector2D p){return Patch::NodalDof{(p.x-.25)*(p.x-.75),2*p.x-1,0,0};});
    const auto components=Panels::discoverNativeTraceComponents(mesh,twice,false,{{.25,0},{.25,1},{.75,0},{.75,1}});
    check(components.curves.size()==2 && components.negativeRegions.size()==1,"multiple side components/negative connectivity");
    for(double line:{0.,1.})for(int sign:{-1,1}){
        const auto source=traceSource(mesh,traceProfile(0,line,1,0,sign));
        const auto topo=Panels::discoverNativeTraceComponents(mesh,source,false,{{line,0},{line,1}});
        check(topo.curves.size()==1,"physical boundary zero chain lost");
        for(const auto& g:topo.curves[0].graphs())check(g.nativeSide->minus().available()==(line==1)
            && g.nativeSide->plus().available()==(line==0),"physical exterior phase extrapolated");}
    const auto partial=traceSource(mesh,traceProfile(0,.5,1,0),false,0,.5);
    const auto selected=Panels::discoverNativeTraceComponents(mesh,partial,false,{{.5,0},{.5,1}});
    for(const auto& g:selected.curves[0].graphs())check(g.nativeSide->minus().available() && !g.nativeSide->plus().available(),"inactive side phase extrapolated");
    for(int order:{2,3,4}){const auto source=traceSource(mesh,traceProfile(0,.5,order,0));
        rejects([&]{Panels::discoverNativeTraceComponents(mesh,source,false,{{.5,0},{.5,1}});},"refinement required");}
    const auto crossing=traceSource(mesh,traceProfile(0,.5,1,1));
    rejects([&]{Panels::discoverNativeTraceComponents(mesh,crossing,false,{{.5,0},{.5,1}});},"refinement required");
    const auto zero=traceSource(mesh,[](Vector2D){return Patch::NodalDof{};});
    rejects([&]{Panels::discoverNativeTraceComponents(mesh,zero,false);},"refinement required");
    std::cout<<"NATIVE_SIDE_TOPOLOGY,"<<layouts<<','<<supports<<','<<hanging<<','<<side<<','<<ordinary<<",PASS\n";
}
}
int main(int argc,char** argv){try{
    std::cout<<std::scientific<<std::setprecision(16);
    if(argc==2 && std::string(argv[1])=="--numerical-q4"){
        Fixture f(true,false,1,52);nativeQuarticQueryGuard(f);numericalQuarticTransportPoint(f);return 0;}
    if(argc==2 && std::string(argv[1])=="--side-topology"){sideTopologyContract();return 0;}
    if(argc==2 && std::string(argv[1])=="--reproducer"){
        Fixture f(true,false,1,52);const auto result=f.integrate(f.regions);
        for(const auto& arc:result.arcs)std::cout<<"REPRODUCER_ARC,"<<arc.regionId<<','<<arc.width<<','<<arc.axialMoment<<','<<arc.radialMoment<<'\n';
        return 0;}
    exactPowerOfTwoScalingContract();
    for(int exponent:{24,42,48,52})for(bool radial:{false,true})for(bool x:{false,true})for(int sign:{-1,1}){
        Fixture f(x,radial,sign,exponent);verify(f);}
    Fixture f(true,false,1,52);nativeQuarticQueryGuard(f);numericalQuarticTransportPoint(f);auto overlap=f.regions;auto duplicate=overlap[0];duplicate.regionId=2;overlap.push_back(duplicate);
    rejects([&]{f.integrate(overlap);},"different overlapping region owners");
    auto missing=f.regions;missing.pop_back();rejects([&]{f.integrate(missing);},"no actual region owner");
    Fixture absent(true,false,1,52);
    absent.cuts.erase(std::remove_if(absent.cuts.begin(),absent.cuts.end(),[&](const auto& c){return c.begin().y==absent.origin.y && c.end().y==absent.origin.y;}),absent.cuts.end());
    rejects([&]{absent.integrate(absent.regions);},"external region side is missing");
    Fixture other(true,false,1,52);f.cuts[0]=std::move(other.cuts[0]);
    rejects([&]{f.integrate(f.regions);},"different snapshot");
    traceContract();sideTopologyContract();
    std::cout<<"I12a9r1 root enclosure PASS cases="<<cases<<" rejections="<<rejections<<'\n';return 0;
}catch(const std::exception& e){std::cerr<<"I12a9r1 root enclosure FAIL: "<<e.what()<<'\n';return 1;}}
