// Core V2a contract: independent analytic moments, physical XY/RZ measure,
// affine cuts, exact zero events and a real 2:1 parent/subcell geometry patch.
#include "Core/Geometry/PlanarCellCut.h"
#include "Core/Mesh/CartesianMeasure.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cfenv>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
// BEGIN V4a1a ADDED CORE INCLUDES
#include <cstdint>
#include <cstring>
// END V4a1a ADDED CORE INCLUDES

namespace
{
namespace Cut = LAMP::Geometry::PlanarCellCut;
using Cut::Point;
using Polygon = std::vector<Point>;
constexpr long double tolerance = 256*std::numeric_limits<long double>::epsilon();
// Physical 2*pi is the existing shared measure coefficient (stored as double).
// Independent analytic integrals use that declared coefficient; no tolerance
// is relaxed when long double has more precision on another compiler.
const long double pi = static_cast<long double>(CartesianMeasure::twoPi)/2;
long double peakScaledError = 0;
unsigned assertions = 0;
void check(bool condition,const char* name)
{
    ++assertions;
    if (!condition) throw std::runtime_error(name);
}
void near(long double actual,long double expected,long double scale,const char* name,
    long double relative=tolerance)
{
    ++assertions;
    if (!std::isfinite(actual) || !std::isfinite(expected) || !(scale>0))
        throw std::runtime_error(std::string(name)+": invalid comparison");
    const auto error=std::abs(actual-expected)/scale;
    peakScaledError=std::max(peakScaledError,error);
    if (error>relative) throw std::runtime_error(std::string(name)+": analytic discrepancy");
}
template<class F> void rejected(F&& function,const char* name)
{
    bool caught=false;
    try { function(); } catch (const std::exception&) { caught=true; }
    check(caught,name);
}
Polygon reversed(Polygon polygon)
{ std::reverse(polygon.begin(),polygon.end());return polygon; }
Polygon rectangle(long double a,long double b,long double c,long double d)
{ return {{a,c},{b,c},{b,d},{a,d}}; }
long double measured(const Polygon& polygon,bool rz,Point origin={})
{ return Cut::polygonMeasureInFrame(polygon,origin,rz); }
void moments()
{
    const Polygon box=rectangle(-1,3,-2,1);
    for (const auto& polygon:{box,reversed(box)}) {
        const auto m=Cut::polygonMoments(polygon);
        near(m.area,12,12,"rectangle area");
        near(m.firstX,12,12,"rectangle x moment");
        near(m.firstY,-6,12,"rectangle y moment");
        near(measured(polygon,false,{7,-19}),12,12,"XY frame invariance");
        near(measured(polygon,true,{7,-19}),192*pi,192*pi,"RZ actual radius moment");
    }
    const Polygon triangle{{0,0},{2,0},{0,3}};
    const auto m=Cut::polygonMoments(triangle);
    near(m.area,3,3,"triangle area");
    near(m.firstX,2,2,"triangle x moment");
    near(m.firstY,3,3,"triangle y moment");
    near(measured(triangle,true),4*pi,4*pi,"axis triangle radial integral");
    near(measured(triangle,true,{4,8}),28*pi,28*pi,"triangle global radius");
    // A local frame may contain negative x while every physical radius is positive.
    const Polygon shifted{{-4,0},{-2,0},{-4,3}};
    near(measured(shifted,true,{4,-6}),4*pi,4*pi,"explicit translated frame");
    const Polygon collinear{{0,0},{.5L,0},{1,0},{0,0}};
    check(Cut::polygonMoments(collinear).area==0,"collinear polygon not zero");
    check(measured({},false)==0&&measured({{0,0}},true)==0,"lower-dimensional measure");
    const auto cut=Cut::clipConvexPolygon(rectangle(-1,1,-1,1),{1,1},0);
    const auto cm=Cut::polygonMoments(cut);
    near(cm.area,2,2,"diagonal cut area");
    near(cm.firstX,-2.L/3,2,"diagonal cut x moment");
    near(cm.firstY,-2.L/3,2,"diagonal cut y moment");
    const auto scaled=Cut::clipConvexPolygon(rectangle(-1,1,-1,1),{7,7},0);
    near(measured(scaled,false),2,2,"raw normal plane scaling");
    const auto offset=Cut::clipConvexPolygon(rectangle(0,2,0,1),{2,0},1);
    near(measured(offset,false),.5L,.5L,"plane offset uses unnormalized normal");
    near(measured(offset,true),pi/4,pi/4,"axis strip exact radial weight");
    const auto positive=Cut::clipConvexPolygon(rectangle(-1,1,-1,1),{-1,-1},0);
    near(measured(cut,true,{2,0})+measured(positive,true,{2,0}),16*pi,16*pi,
        "polygon complement physical closure");
}
void interface(const Cut::LinearTriangleCut& cut,Point gradient,const char* name)
{
    check(cut.interface.has_value(),name);
    const auto norm=std::hypot(gradient.x,gradient.y);
    near(cut.unitNormal.x,gradient.x/norm,1,"liquid-to-gas normal x");
    near(cut.unitNormal.y,gradient.y/norm,1,"liquid-to-gas normal y");
    const auto a=(*cut.interface)[0],b=(*cut.interface)[1];
    check((-cut.unitNormal.y)*(b.x-a.x)+cut.unitNormal.x*(b.y-a.y)>0,
        "interface orientation does not leave negative phase on left");
}
void triangles()
{
    const std::array<Point,3> triangle{{{0,0},{2,0},{0,3}}};
    for (bool clockwise:{false,true}) {
        auto points=triangle;
        std::array<long double,3> scalar{{-1,1,1}};
        if (clockwise) { std::swap(points[1],points[2]);std::swap(scalar[1],scalar[2]); }
        const auto cut=Cut::cutTriangle(points,scalar);
        check(cut.kind==Cut::TriangleRegion::Cut&&cut.negativeCount==3,"triangle cut classification");
        check(cut.zeroVertexMask==0&&cut.zeroEdgeMask==0,"ordinary cut falsely tagged zero");
        near(Cut::negativeMeasure(cut,{},false),.75L,.75L,"affine triangle negative area");
        near(Cut::negativeMeasure(cut,{},true),pi/2,pi/2,"affine triangle exact radial integral");
        interface(cut,{1,2.L/3},"triangle missing zero segment");
        near(Cut::negativeMeasure(points,scalar,{},false),.75L,.75L,"fast triangle analytic area");
        near(Cut::negativeMeasure(points,scalar,{},true),pi/2,pi/2,"fast triangle analytic RZ integral");
        for (auto& value:scalar) value=-value;
        const auto complement=Cut::cutTriangle(points,scalar);
        check(complement.negativeCount==4,"triangle complement needs four fixed vertices");
        near(Cut::negativeMeasure(cut,{4,0},true)+Cut::negativeMeasure(complement,{4,0},true),
            28*pi,28*pi,"triangle complement RZ closure");
        interface(complement,{-1,-2.L/3},"complement missing zero segment");
    }
    const auto maximum=std::numeric_limits<long double>::max();
    const auto large=Cut::cutTriangle(triangle,{{-maximum,maximum,maximum}});
    near(Cut::negativeMeasure(large,{},false),.75L,.75L,"opposite large scalars overflowed crossing");
    const auto vertex=Cut::cutTriangle(triangle,{{0,-1,1}});
    check(vertex.kind==Cut::TriangleRegion::Cut&&vertex.zeroVertexMask==1&&vertex.zeroEdgeMask==0,
        "exact vertex event not retained");
    near(Cut::negativeMeasure(vertex,{},false),1.5L,1.5L,"zero vertex analytic area");
    interface(vertex,{-.5L,1.L/3},"zero vertex lost full segment");
    for (long double sign:{-1.L,1.L}) {
        const auto edge=Cut::cutTriangle(triangle,{{0,0,sign}});
        check(edge.zeroVertexMask==3&&edge.zeroEdgeMask==1,"coincident edge event masks");
        check(edge.kind==(sign<0?Cut::TriangleRegion::Negative:Cut::TriangleRegion::Positive),
            "coincident edge pure classification");
        near(Cut::negativeMeasure(edge,{},false),sign<0?3.L:0.L,3,"coincident edge phase area");
        interface(edge,{0,sign},"coincident edge lost segment");
        const auto isolated=Cut::cutTriangle(triangle,{{0,sign,sign}});
        check(isolated.zeroVertexMask==1&&!isolated.interface,"isolated zero point became segment");
        near(Cut::negativeMeasure(isolated,{},false),sign<0?3.L:0.L,3,"isolated vertex phase area");
    }
    const auto allZero=Cut::cutTriangle(triangle,{{0,0,0}});
    check(allZero.kind==Cut::TriangleRegion::AllZeroScalar&&allZero.zeroVertexMask==7
        &&allZero.zeroEdgeMask==7&&!allZero.interface,"positive-area all-zero classification");
    rejected([&]{Cut::negativeMeasure(allZero,{},false);},"all-zero area arbitrarily assigned phase");
    rejected([&]{Cut::negativeMeasure(triangle,{{0,0,0}},{},true);},"fast all-zero area arbitrarily assigned phase");
    // Resolve a tiny physical liquid phase. Compare against its own size,
    // never against unit capacity; no fraction clipping/floor may hide it.
    const std::array<Point,3> unit{{{0,0},{1,0},{0,1}}};
    constexpr long double delta=1e-9L;
    const auto small=Cut::cutTriangle(unit,{{-delta,1,1}});
    const auto t=delta/(1+delta);
    const auto area=t*t/2, radial=pi*t*t*t/3;
    check(Cut::negativeMeasure(small,{},false)>0,"small liquid area erased");
    near(Cut::negativeMeasure(small,{},false),area,area,"small phase relative area");
    near(Cut::negativeMeasure(small,{},true),radial,radial,"small phase relative RZ volume");
    near(Cut::negativeMeasure(unit,{{-delta,1,1}},{},false),area,area,"fast small phase relative area");
    near(Cut::negativeMeasure(unit,{{-delta,1,1}},{},true),radial,radial,"fast small phase relative RZ volume");
    const auto scalarSmall=Cut::cutTriangle(triangle,{{-1e-200L,1e-200L,1e-200L}});
    near(Cut::negativeMeasure(scalarSmall,{},false),.75L,.75L,"scalar scale changed phase geometry");
}
// Independent analytic integration of vertical strips in a rectangle under
// x+y <= intercept. No polygon clipping, fan rule or production quadrature
// participates in this reference. Full and affine partial strips integrate
// exactly with quadratic/cubic antiderivatives.
Cut::PolygonMoments lineRectangle(long double a,long double b,long double c,long double d,long double intercept)
{
    Cut::PolygonMoments result;
    const auto fullEnd=std::clamp(intercept-d,a,b);
    const auto cutBegin=fullEnd,cutEnd=std::clamp(intercept-c,a,b);
    if (fullEnd>a) {
        result.area=(fullEnd-a)*(d-c);
        result.firstX=(fullEnd*fullEnd-a*a)*(d-c)/2;
        result.firstY=(fullEnd-a)*(d*d-c*c)/2;
    }
    if (cutEnd>cutBegin) {
        const auto h0=intercept-c;
        const auto width=cutEnd-cutBegin;
        const auto second=(cutEnd*cutEnd-cutBegin*cutBegin)/2;
        const auto third=(cutEnd*cutEnd*cutEnd-cutBegin*cutBegin*cutBegin)/3;
        result.area+=h0*width-second;
        result.firstX+=h0*second-third;
        result.firstY+=((intercept*intercept-c*c)*width-2*intercept*second+third)/2;
    }
    return result;
}
struct FanResult { long double volume=0;std::vector<Point> sharedEnds; };
FanResult fan(const Polygon& corners,Point center,bool rz,Point origin={})
{
    FanResult result;
    for (std::size_t k=0;k<corners.size();++k) {
        const std::array<Point,3> triangle{{center,corners[k],corners[(k+1)%corners.size()]}};
        std::array<long double,3> scalar{};
        for (int i=0;i<3;++i) scalar[i]=triangle[i].x+triangle[i].y-1.375L;
        const auto negative=Cut::cutTriangle(triangle,scalar);
        const auto fastLiquid=Cut::negativeMeasure(triangle,scalar,origin,rz);
        for (auto& value:scalar) value=-value;
        const auto positive=Cut::cutTriangle(triangle,scalar);
        const auto fullArea=std::abs((triangle[1].x-triangle[0].x)*(triangle[2].y-triangle[0].y)
            -(triangle[1].y-triangle[0].y)*(triangle[2].x-triangle[0].x))/2;
        const auto fullVolume=fullArea*(rz?2*pi*(origin.x+(triangle[0].x+triangle[1].x+triangle[2].x)/3):1);
        const auto liquid=Cut::negativeMeasure(negative,origin,rz);
        near(fastLiquid,liquid,fullVolume,"fan fast/full measure agreement");
        near(liquid+Cut::negativeMeasure(positive,origin,rz),fullVolume,fullVolume,"fan local complement");
        result.volume+=liquid;
        if (negative.interface) for (auto p:*negative.interface) if (p.x==1) result.sharedEnds.push_back(p);
    }
    return result;
}
void amrPatch()
{
    // One true coarse leaf has the hanging node (1,.5) on its shared edge.
    // The adjacent parent [1,2]x[0,1] is replaced by all four .5x.5 children.
    // Every scalar is the same globally affine x+y-1.375 on actual vertices.
    const Polygon coarse{{0,0},{1,0},{1,.5L},{1,1},{0,1}};
    for (bool rz:{false,true}) {
        const auto coarseResult=fan(coarse,{.5L,.5L},rz);
        auto exact=lineRectangle(0,1,0,1,1.375L);
        near(coarseResult.volume,rz?2*pi*exact.firstX:exact.area,rz?pi:1,"coarse hanging-node fan");
        const auto parent=fan(rectangle(1,2,0,1),{1.5L,.5L},rz);
        long double childSum=0;unsigned sharedCrossings=0;
        for (int ix=0;ix<2;++ix) for (int iy=0;iy<2;++iy) {
            const auto a=1+ix*.5L,b=a+.5L,c=iy*.5L,d=c+.5L;
            const auto child=fan(rectangle(a,b,c,d),{(a+b)/2,(c+d)/2},rz);
            const auto childExact=lineRectangle(a,b,c,d,1.375L);
            const auto childCapacity=(b-a)*(d-c)*(rz?pi*(a+b):1);
            near(child.volume,rz?2*pi*childExact.firstX:childExact.area,childCapacity,"actual child analytic cut");
            childSum+=child.volume;
            for (const auto p:child.sharedEnds) {
                check(p.y==.375L,"fine shared trace crossing inconsistent");++sharedCrossings;
            }
        }
        near(childSum,parent.volume,rz?3*pi:1,"parent/four-child volume sum");
        const auto total=lineRectangle(0,2,0,1,1.375L);
        near(coarseResult.volume+childSum,rz?2*pi*total.firstX:total.area,rz?4*pi:2,"2:1 patch global volume");
        check(!coarseResult.sharedEnds.empty()&&sharedCrossings>0,"2:1 patch did not cross shared subface");
        for (const auto p:coarseResult.sharedEnds) check(p.y==.375L,"coarse shared trace crossing inconsistent");
    }
}
void thinBinaryPolygons()
{
    // These represented coordinates enclose q^2 despite cross products of
    // order one. Elementary area/centroid formulas are independent of Core;
    // every check is relative to the actual thin phase, never unit capacity.
    const long double q=std::ldexp(1.L,-28),thinTolerance=1e-12L;
    const auto initialAssertions=assertions;
    const auto triangleArea=q*q/2,cutArea=q*q/8,quadArea=q*q;
    const Point triangleCentroid{(2+3*q)/3,(2+q)/3};
    const Point cutCentroid{(2+3*q)/6,(2+q)/6};
    const Point quadCentroid{.5L+q,.5L+.5L*q};
    const auto verifyPolygon=[&](const Polygon& polygon,long double area,Point centroid) {
        const auto actual=Cut::polygonMoments(polygon);
        near(actual.area,area,area,"binary thin polygon area",thinTolerance);
        near(actual.firstX,area*centroid.x,area*centroid.x,"binary thin polygon x moment",thinTolerance);
        near(actual.firstY,area*centroid.y,area*centroid.y,"binary thin polygon y moment",thinTolerance);
        for(long double radiusOrigin:{.5L,1.L})for(bool rz:{false,true}) {
            const auto expected=area*(rz?2*pi*(radiusOrigin+centroid.x):1);
            near(Cut::polygonMeasureInFrame(polygon,{radiusOrigin,0},rz),expected,expected,
                "binary thin polygon physical measure",thinTolerance);
        }
    };
    for(long double translation:{0.L,16.L})for(bool clockwise:{false,true}) {
        std::array<Point,3> triangle{{{translation,translation},
            {translation+1+q,translation+1},{translation+1+2*q,translation+1+q}}};
        Polygon quad{{translation,translation},{translation+1+q,translation+1},
            {translation+1+2*q,translation+1+q},{translation+q,translation+q}};
        if(clockwise) {std::swap(triangle[1],triangle[2]);quad=reversed(quad);}
        const Point translatedTriangle{translation+triangleCentroid.x,translation+triangleCentroid.y};
        const Point translatedCut{translation+cutCentroid.x,translation+cutCentroid.y};
        const Point translatedQuad{translation+quadCentroid.x,translation+quadCentroid.y};
        verifyPolygon(Polygon(triangle.begin(),triangle.end()),triangleArea,translatedTriangle);
        verifyPolygon(quad,quadArea,translatedQuad);
        const std::array<long double,3> pureScalar{{-1,-1,-1}},cutScalar{{-1,1,1}},positiveScalar{{1,1,1}};
        const auto pure=Cut::cutTriangle(triangle,pureScalar),cut=Cut::cutTriangle(triangle,cutScalar);
        const auto positive=Cut::cutTriangle(triangle,positiveScalar);
        check(pure.kind==Cut::TriangleRegion::Negative&&pure.negativeCount==3,
            "binary thin triangle pure classification");
        check(cut.kind==Cut::TriangleRegion::Cut&&cut.negativeCount==3&&cut.interface.has_value(),
            "binary thin triangle affine classification");
        check(positive.kind==Cut::TriangleRegion::Positive&&positive.negativeCount==0,
            "binary thin triangle positive classification");
        verifyPolygon(Polygon(cut.negativeVertices.begin(),cut.negativeVertices.begin()+cut.negativeCount),
            cutArea,translatedCut);
        for(long double radiusOrigin:{.5L,1.L})for(bool rz:{false,true}) {
            const Point origin{radiusOrigin,0};
            const auto pureExpected=triangleArea*(rz?2*pi*(radiusOrigin+translatedTriangle.x):1);
            const auto cutExpected=cutArea*(rz?2*pi*(radiusOrigin+translatedCut.x):1);
            near(Cut::negativeMeasure(pure,origin,rz),pureExpected,pureExpected,
                "binary thin pure triangle physical measure",thinTolerance);
            near(Cut::negativeMeasure(triangle,pureScalar,origin,rz),pureExpected,pureExpected,
                "binary thin fast pure triangle physical measure",thinTolerance);
            near(Cut::negativeMeasure(cut,origin,rz),cutExpected,cutExpected,
                "binary thin cut triangle physical measure",thinTolerance);
            near(Cut::negativeMeasure(triangle,cutScalar,origin,rz),cutExpected,cutExpected,
                "binary thin fast cut triangle physical measure",thinTolerance);
            check(Cut::negativeMeasure(positive,origin,rz)==0,
                "binary thin theoretical zero phase is nonzero");
            check(Cut::negativeMeasure(triangle,positiveScalar,origin,rz)==0,
                "binary thin fast theoretical zero phase is nonzero");
        }
    }
    std::cout<<"group=binary-thin-polygons PASS assertions="<<(assertions-initialAssertions)
        <<" q_exponent=-28 thin_relative_tolerance="<<thinTolerance<<'\n';
}
void inverseAndLegacy()
{
    const Cut::Rectangle cell{{0,-.25},{2,.75}};
    for (bool rz:{false,true}) for (Vector2D normal:{Vector2D{1,0},Vector2D{0,1},Vector2D{3,4},Vector2D{-3,4}})
        for (long double fraction:{0.L,1e-9L,.125L,.5L,.875L,1.L}) {
            const auto total=Cut::capacity(cell,rz);
            const auto reconstruction=Cut::reconstruct(cell,normal,fraction,rz);
            near(reconstruction.volume,fraction*total,total,"legacy PLIC inverse volume",Cut::relativeVolumeResolution);
            near(Cut::volume(cell,reconstruction.normal,reconstruction.offset,rz),reconstruction.volume,total,
                "legacy PLIC forward consistency",Cut::relativeVolumeResolution);
            const auto complement=Cut::volume(cell,{-reconstruction.normal.x,-reconstruction.normal.y},
                -reconstruction.offset,rz);
            near(reconstruction.volume+complement,total,total,"legacy PLIC complement");
        }
    // Existing rectangle kernel compared with the independent strip integral.
    const auto exact=lineRectangle(0,2,0,1,1.375L);
    const Cut::Rectangle patch{{0,0},{2,1}};
    const auto offset=(1.375L-1.5L)/std::sqrt(2.L);
    for (bool rz:{false,true}) {
        const auto expected=rz?2*pi*exact.firstX:exact.area;
        near(Cut::volume(patch,{1,1},offset,rz),expected,rz?4*pi:2,"legacy oblique physical integral");
        const auto overlap=Cut::overlapInFrame(patch,rectangle(0,2,0,1),{0,0},{1,1},offset,rz);
        near(overlap.liquid,expected,rz?4*pi:2,"legacy shared-frame overlap");
    }
}
void invalidInput()
{
    const auto nan=std::numeric_limits<long double>::quiet_NaN();
    const auto inf=std::numeric_limits<long double>::infinity();
    const Polygon concave{{0,0},{2,0},{.5L,.5L},{2,2},{0,2}};
    const Polygon crossed{{0,0},{1,1},{0,1},{1,0}};
    rejected([&]{Cut::polygonMoments(concave);},"concave polygon moments accepted");
    rejected([&]{Cut::clipConvexPolygon(concave,{1,0},.5L);},"concave polygon clipping accepted");
    rejected([&]{Cut::polygonMoments(crossed);},"crossed polygon accepted");
    rejected([&]{Cut::polygonMoments({{0,0},{1,0},{nan,1}});},"nonfinite polygon accepted");
    rejected([&]{measured(rectangle(-.1L,1,0,1),true);},"negative physical radius accepted");
    rejected([&]{measured(rectangle(0,1,0,1),false,{nan,0});},"nonfinite frame accepted");
    rejected([&]{Cut::clipConvexPolygon(rectangle(0,1,0,1),{0,0},0);},"zero plane normal accepted");
    rejected([&]{Cut::clipConvexPolygon(rectangle(0,1,0,1),{1,0},inf);},"nonfinite plane accepted");
    const std::array<Point,3> triangle{{{0,0},{1,0},{0,1}}};
    rejected([&]{Cut::cutTriangle(triangle,{{0,inf,1}});},"nonfinite triangle scalar accepted");
    rejected([&]{Cut::cutTriangle({{{0,0},{1,0},{2,0}}},{{-1,0,1}});},"degenerate triangle accepted");
    rejected([&]{Cut::cutTriangle({{{0,0},{0,0},{1,1}}},{{-1,0,1}});},"repeated triangle vertex accepted");
    const auto tiny=std::numeric_limits<long double>::min();
    rejected([&]{Cut::negativeMeasure(Cut::cutTriangle(triangle,{{-tiny,1,1}}),{},false);},
        "unrepresentable tiny area silently erased");
    rejected([&]{Cut::negativeMeasure(triangle,{{-tiny,1,1}},{},false);},
        "fast unrepresentable tiny area silently erased");
    rejected([&]{Cut::negativeMeasure(triangle,{{1,1,1}},{-.125L,0},true);},
        "fast positive triangle accepted negative physical radius");
    const Cut::Rectangle cell{{0,0},{1,1}};
    rejected([&]{Cut::reconstruct(cell,{1,0},-.01L,false);},"negative alpha clipped");
    rejected([&]{Cut::reconstruct(cell,{1,0},1.01L,false);},"alpha above one clipped");
    rejected([&]{Cut::reconstruct(cell,{1,0},nan,true);},"nonfinite alpha accepted");
}
}
// BEGIN V4a1a ADDED ANCHORED TEST BLOCK
namespace
{
// Independent twofold reference: supplied components remain separate until final scalar output.
struct AnchoredRefNumber {long double hi=0,lo=0;};
AnchoredRefNumber arSum(long double a,long double b)
{const auto s=a+b,v=s-a;return {s,(a-(s-v))+(b-v)};}
AnchoredRefNumber arAdd(AnchoredRefNumber a,AnchoredRefNumber b)
{const auto s=arSum(a.hi,b.hi),t=arSum(a.lo,b.lo),u=arSum(s.lo,t.hi),v=arSum(s.hi,u.hi);return arSum(v.hi,v.lo+u.lo+t.lo);}
AnchoredRefNumber arNeg(AnchoredRefNumber a){return {-a.hi,-a.lo};}
AnchoredRefNumber arSub(AnchoredRefNumber a,AnchoredRefNumber b){return arAdd(a,arNeg(b));}
AnchoredRefNumber arMul(AnchoredRefNumber a,AnchoredRefNumber b)
{const auto p=a.hi*b.hi;auto e=arAdd({std::fma(a.hi,b.hi,-p),0},{a.hi*b.lo,0});e=arAdd(e,{a.lo*b.hi,0});e=arAdd(e,{a.lo*b.lo,0});return arAdd({p,0},e);}
AnchoredRefNumber arDiv(AnchoredRefNumber a,AnchoredRefNumber b)
{const auto q=a.hi/b.hi;const auto r=arSub(a,arMul(b,{q,0}));return arAdd({q,0},{(r.hi+r.lo)/b.hi,0});}
long double arValue(AnchoredRefNumber a){return a.hi+a.lo;}
struct AnchoredRefPoint {AnchoredRefNumber x,y;};
using AnchoredSite=LAMP::Geometry::PlanarCellCut::AnchoredPoint;
using AnchoredSites=std::vector<AnchoredSite>;
using AnchoredRefPolygon=std::vector<AnchoredRefPoint>;
AnchoredRefPoint arPoint(AnchoredSite p){return {arSum(p.anchor.x,p.offset.x),arSum(p.anchor.y,p.offset.y)};}
AnchoredRefPoint arSubPoint(AnchoredRefPoint a,AnchoredRefPoint b){return {arSub(a.x,b.x),arSub(a.y,b.y)};}
AnchoredRefNumber arCross(AnchoredRefPoint a,AnchoredRefPoint b){return arSub(arMul(a.x,b.y),arMul(a.y,b.x));}
AnchoredRefPolygon arPoints(const AnchoredSites& p){AnchoredRefPolygon out;for(auto q:p)out.push_back(arPoint(q));return out;}
struct AnchoredRefMoments {AnchoredRefNumber area,x,y;};
AnchoredRefMoments arMoments(const AnchoredRefPolygon& p)
{
    if(p.size()<3)return {};const auto base=p.front();AnchoredRefNumber area,x,y;
    for(std::size_t i=0;i<p.size();++i){const auto a=arSubPoint(p[i],base),b=arSubPoint(p[(i+1)%p.size()],base);const auto c=arCross(a,b);
        area=arAdd(area,c);x=arAdd(x,arMul(arAdd(a.x,b.x),c));y=arAdd(y,arMul(arAdd(a.y,b.y),c));}
    if(arValue(area)<0){area=arNeg(area);x=arNeg(x);y=arNeg(y);}
    area=arDiv(area,{2,0});x=arAdd(arMul(base.x,area),arDiv(x,{6,0}));y=arAdd(arMul(base.y,area),arDiv(y,{6,0}));return {area,x,y};
}
AnchoredRefNumber arMeasure(const AnchoredRefPolygon& p,LAMP::Geometry::PlanarCellCut::Point origin,bool rz)
{const auto m=arMoments(p);return rz?arMul({static_cast<long double>(CartesianMeasure::twoPi),0},arAdd(arMul({origin.x,0},m.area),m.x)):m.area;}
AnchoredRefPoint arCrossing(AnchoredRefPoint a,AnchoredRefPoint b,long double va,long double vb)
{
    const auto wa=std::abs(va),wb=std::abs(vb);const auto t=arDiv({wa,0},arAdd({wa,0},{wb,0}));
    const auto d=arSubPoint(b,a);return {arAdd(a.x,arMul(t,d.x)),arAdd(a.y,arMul(t,d.y))};
}
AnchoredRefPolygon arTriangleNegative(const std::array<AnchoredSite,3>& triangle,const std::array<long double,3>& scalar)
{
    AnchoredRefPolygon p;for(unsigned i=0;i<3;++i){const auto j=(i+1)%3;const auto a=arPoint(triangle[i]),b=arPoint(triangle[j]);
        if(scalar[i]<=0)p.push_back(a);if((scalar[i]<0&&scalar[j]>0)||(scalar[i]>0&&scalar[j]<0))p.push_back(arCrossing(a,b,scalar[i],scalar[j]));}return p;
}
AnchoredRefPolygon arHalfPlane(const AnchoredRefPolygon& p,bool top)
{
    AnchoredRefPolygon out;if(p.empty())return out;
    for(std::size_t i=0;i<p.size();++i){const auto a=p[i],b=p[(i+1)%p.size()];auto va=arSub(a.y,{.5L,0}),vb=arSub(b.y,{.5L,0});if(top){va=arNeg(va);vb=arNeg(vb);}
        const bool ina=arValue(va)<=0,inb=arValue(vb)<=0;if(ina)out.push_back(a);
        if(ina!=inb){const auto t=arDiv(va,arSub(va,vb));const auto d=arSubPoint(b,a);out.push_back({arAdd(a.x,arMul(t,d.x)),arAdd(a.y,arMul(t,d.y))});}}return out;
}
bool arSame(AnchoredSite a,AnchoredSite b)
{const auto d=arSubPoint(arPoint(a),arPoint(b));return d.x.hi==0&&d.x.lo==0&&d.y.hi==0&&d.y.lo==0;}
unsigned anchoredCoreAssertions=0,anchoredCoreCases=0,anchoredCoreGuards=0;
long double anchoredCorePeak=0,anchoredFixedInventory=0,anchoredInverseInventory=0,anchoredPointCollapse=0;
unsigned anchoredInverseIterations=0,anchoredInverseExpansions=0;
constexpr long double anchoredOwnGate=1e-12L;
void anchoredCheck(bool ok,const char* reason){++anchoredCoreAssertions;if(!ok)throw std::runtime_error(std::string("anchored Core: ")+reason);}
void anchoredNear(long double a,long double e,const char* reason)
{
    anchoredCheck(std::isfinite(a)&&std::isfinite(e)&&e>0,reason);const auto error=std::abs(a-e)/e;anchoredCorePeak=std::max(anchoredCorePeak,error);
    if(error>anchoredOwnGate)throw std::runtime_error(std::string("anchored Core own-phase gate: ")+reason);
}
template<class F>void anchoredReject(F&& f,const char* reason)
{bool threw=false;try{f();}catch(const std::exception&){threw=true;}anchoredCheck(threw,reason);++anchoredCoreGuards;}
std::uint64_t anchoredDoubleBits(double d){std::uint64_t b=0;std::memcpy(&b,&d,sizeof(d));return b;}
const std::array<Point,4> anchoredFanBoundary{{{.5L,-.375L},{.625L,-.375L},{.625L,-.25L},{.5L,-.25L}}};
const std::array<double,4> anchoredFanScalars{{0.00032032603052658398,0.094017880689772709,0.17653099104249845,0.088857572330733547}};
constexpr double anchoredFanCentre=-0x1.7ba31be9009c0p-17;
constexpr double anchoredFanTarget=0x1.f5409cbd24a41p-24;
constexpr double anchoredFanCapacity=0x1.c463abeccb2bbp-5;
constexpr double anchoredFanAlpha=0x1.1ba69cc747717p-19;
const Point anchoredFanOrigin{0,.5L},anchoredFanHub{.5625L,-.3125L};
std::array<AnchoredSite,3> anchoredFanTriangle(unsigned i)
{return {{{anchoredFanBoundary[i],{}},{anchoredFanBoundary[(i+1)%4],{}},{anchoredFanHub,{}}}};}
std::array<long double,3> anchoredFanValues(unsigned i,long double centre)
{return {{static_cast<long double>(anchoredFanScalars[i]),static_cast<long double>(anchoredFanScalars[(i+1)%4]),centre}};}
struct AnchoredPositiveSum{long double sum=0,error=0;void add(long double v){const auto y=v-error,t=sum+y;error=(t-sum)-y;sum=t;}};
long double anchoredFanFast(long double centre)
{AnchoredPositiveSum sum;for(unsigned i=0;i<4;++i)sum.add(Cut::negativeAnchoredMeasure(anchoredFanTriangle(i),anchoredFanValues(i,centre),anchoredFanOrigin,true));return sum.sum;}
long double anchoredInspectFan(long double centre)
{
    AnchoredRefNumber independentL,independentG;AnchoredPositiveSum publishedL,publishedG;std::vector<AnchoredSite> gamma;
    for(unsigned i=0;i<4;++i){const auto tri=anchoredFanTriangle(i);const auto s=anchoredFanValues(i,centre);const auto cut=Cut::cutAnchoredTriangle(tri,s);auto positive=s;for(auto& v:positive)v=-v;
        const auto gas=Cut::cutAnchoredTriangle(tri,positive);anchoredCheck(cut.kind==Cut::TriangleRegion::Cut&&gas.kind==Cut::TriangleRegion::Cut,"fan cut classification");
        const AnchoredSites liquid(cut.negativeVertices.begin(),cut.negativeVertices.begin()+cut.negativeCount),air(gas.negativeVertices.begin(),gas.negativeVertices.begin()+gas.negativeCount);
        anchoredCheck(cut.interface.has_value()&&gas.interface.has_value(),"fan interface absent");
        for(bool rz:{false,true}){const auto il=arMeasure(arPoints(liquid),anchoredFanOrigin,rz),ig=arMeasure(arPoints(air),anchoredFanOrigin,rz);
            anchoredNear(Cut::negativeAnchoredMeasure(cut,anchoredFanOrigin,rz),arValue(il),"published negative actual sites");
            anchoredNear(Cut::negativeAnchoredMeasure(gas,anchoredFanOrigin,rz),arValue(ig),"published positive actual sites");
            anchoredNear(Cut::negativeAnchoredMeasure(tri,s,anchoredFanOrigin,rz),arValue(il),"fast and published same negative sites");
            anchoredNear(Cut::negativeAnchoredMeasure(tri,positive,anchoredFanOrigin,rz),arValue(ig),"fast and published same positive sites");
            anchoredNear(arValue(il),arValue(arMeasure(arTriangleNegative(tri,s),anchoredFanOrigin,rz)),"independent negative affine clip");
            anchoredNear(arValue(ig),arValue(arMeasure(arTriangleNegative(tri,positive),anchoredFanOrigin,rz)),"independent positive affine clip");
            const auto total=arMeasure(arPoints(AnchoredSites(tri.begin(),tri.end())),anchoredFanOrigin,rz);
            anchoredNear(arValue(arAdd(il,ig)),arValue(total),"actual independent mixed-phase closure");
            if(rz){independentL=arAdd(independentL,il);independentG=arAdd(independentG,ig);publishedL.add(Cut::negativeAnchoredMeasure(cut,anchoredFanOrigin,true));publishedG.add(Cut::negativeAnchoredMeasure(gas,anchoredFanOrigin,true));}}
        const auto a=tri[0].anchor,b=tri[1].anchor,c=tri[2].anchor;
        const auto determinant=(b.x-a.x)*(c.y-a.y)-(b.y-a.y)*(c.x-a.x);
        const Point grad{((s[1]-s[0])*(c.y-a.y)-(s[2]-s[0])*(b.y-a.y))/determinant,((b.x-a.x)*(s[2]-s[0])-(c.x-a.x)*(s[1]-s[0]))/determinant};
        const auto norm=std::hypot(grad.x,grad.y);anchoredCheck(norm>0&&cut.unitNormal.x*grad.x+cut.unitNormal.y*grad.y>0,"fan normal does not face positive scalar");
        const auto d=arSubPoint(arPoint((*cut.interface)[1]),arPoint((*cut.interface)[0]));
        const auto chordLength=std::hypot(arValue(d.x),arValue(d.y)),normalGate=1024*std::numeric_limits<long double>::epsilon();
        anchoredCheck(chordLength>0&&std::abs(cut.unitNormal.x*arValue(d.x)+cut.unitNormal.y*arValue(d.y))<=normalGate*chordLength,"actual fan chord normal orthogonality");
        anchoredCheck(std::abs(std::hypot(cut.unitNormal.x,cut.unitNormal.y)-1)<=normalGate,"actual fan unit normal length");
        anchoredCheck(-cut.unitNormal.y*arValue(d.x)+cut.unitNormal.x*arValue(d.y)>0,"fan interface orientation");
        for(auto site:*cut.interface){const auto p=arSubPoint(arPoint(site),arPoint(tri[0]));const auto value=arAdd({s[0],0},arAdd(arMul({grad.x,0},p.x),arMul({grad.y,0},p.y)));
            anchoredCheck(std::abs(arValue(value))<=anchoredOwnGate*std::max({std::abs(s[0]),std::abs(s[1]),std::abs(s[2])}),"actual interface not on original affine scalar");gamma.push_back(site);}
        ++anchoredCoreCases;
    }
    anchoredCheck(gamma.size()==8,"fan Gamma endpoint count");unsigned unique=0;
    for(std::size_t i=0;i<gamma.size();++i){bool earlier=false;unsigned degree=0;for(std::size_t j=0;j<gamma.size();++j)if(arSame(gamma[i],gamma[j])){anchoredCheck(Cut::sameAnchoredPoint(gamma[i],gamma[j]),"public exact site equality disagrees with oracle");++degree;earlier=earlier||j<i;}
        anchoredCheck(degree==2,"canonical shared fan root does not have degree two");if(!earlier)++unique;}
    anchoredCheck(unique==4,"fan Gamma closed cycle root count");
    anchoredNear(publishedL.sum,arValue(independentL),"whole published liquid actual-site integral");anchoredNear(publishedG.sum,arValue(independentG),"whole published gas actual-site integral");
    anchoredNear(arValue(arAdd(independentL,independentG)),anchoredFanCapacity,"whole independently integrated L/G capacity");
    const auto minor=std::min(static_cast<long double>(anchoredFanTarget),static_cast<long double>(anchoredFanCapacity)-anchoredFanTarget);
    const auto residual=std::abs(arValue(independentL)-anchoredFanTarget)/minor;anchoredCheck(residual<=anchoredOwnGate,"independently measured public fan inventory minor gate");return residual;
}
std::array<Point,4> anchoredFrozenOffsets()
{
    if constexpr(std::numeric_limits<long double>::digits==53)return {{{-0x45de593670f87p-59L,-0x45de593670f87p-59L},{0xfc57066e1f853p-69L,-0xfc57066e1f853p-69L},{0x433339d50dcbfp-68L,0x433339d50dcbfp-68L},{-0x10afe11140803p-65L,0x10afe11140803p-65L}}};
    else if constexpr(std::numeric_limits<long double>::digits==64)return {{{-0x8bbcb26ce1f0deedp-72L,-0x8bbcb26ce1f0deedp-72L},{0xfc57066e1f853349p-81L,-0xfc57066e1f853349p-81L},{0x10ccce754372fc25p-78L,0x10ccce754372fc25p-78L},{-0x857f088a040183bdp-80L,0x857f088a040183bdp-80L}}};
    else if constexpr(std::numeric_limits<long double>::digits==113)return {{{-0x1177964d9c3e1bdd97b28d121356dp-121L,-0x1177964d9c3e1bdd97b28d121356dp-121L},{0x1f8ae0cdc3f0a669109e046a3167bp-130L,-0x1f8ae0cdc3f0a669109e046a3167bp-130L},{0x10ccce754372fc25023131b0a53bp-126L,0x10ccce754372fc25023131b0a53bp-126L},{-0x10afe1114080307793ae86c058947p-129L,0x10afe1114080307793ae86c058947p-129L}}};
    else throw std::runtime_error("anchored Core: missing frozen native-lattice precision certificate");
}
void anchoredOriginalFan()
{
    anchoredCheck(anchoredDoubleBits(anchoredFanAlpha)==4522100464421205783ULL&&anchoredDoubleBits(anchoredFanTarget)==4503410553444780609ULL,"captured source double bits changed");
    anchoredCheck(anchoredDoubleBits(anchoredFanAlpha*anchoredFanCapacity)==anchoredDoubleBits(anchoredFanTarget),"captured original capacity times alpha target changed");
    const auto expected=anchoredFrozenOffsets();
    for(unsigned i=0;i<4;++i)for(int gauge:{-512,0,512})for(bool reverse:{false,true}){
        AnchoredSite a{anchoredFanHub,{}},b{anchoredFanBoundary[i],{}};long double va=std::ldexp(static_cast<long double>(anchoredFanCentre),gauge),vb=std::ldexp(static_cast<long double>(anchoredFanScalars[i]),gauge);if(reverse){std::swap(a,b);std::swap(va,vb);}
        const auto actual=Cut::anchoredZeroCrossing(a,b,va,vb);anchoredCheck(actual.anchor.x==anchoredFanHub.x&&actual.anchor.y==anchoredFanHub.y,"canonical smaller-scalar anchor changed");
        anchoredCheck(actual.offset.x==expected[i].x&&actual.offset.y==expected[i].y,"offset differs from frozen independent native lattice rounding");anchoredCheck(arSame(actual,{anchoredFanHub,expected[i]}),"actual public root sum differs from frozen lattice");++anchoredCoreCases;}
    anchoredFixedInventory=anchoredInspectFan(anchoredFanCentre);
    AnchoredPositiveSum ordinary;for(unsigned i=0;i<4;++i){std::array<Point,3> tri{{anchoredFanBoundary[i],anchoredFanBoundary[(i+1)%4],anchoredFanHub}};ordinary.add(Cut::negativeMeasure(tri,anchoredFanValues(i,anchoredFanCentre),anchoredFanOrigin,true));}
    anchoredPointCollapse=std::abs(ordinary.sum-anchoredFanTarget)/anchoredFanTarget;
    if constexpr(std::numeric_limits<long double>::digits==53)anchoredCheck(anchoredPointCollapse>anchoredOwnGate,"original Point fan negative control ceased to expose quantization");
    long double span=.125L;for(auto s:anchoredFanScalars)span=std::max(span,std::abs(static_cast<long double>(s)));
    long double lo=-span,hi=span,vl=anchoredFanFast(lo),vh=anchoredFanFast(hi);
    while((vl<anchoredFanTarget||vh>anchoredFanTarget)&&anchoredInverseExpansions<32){if(vl<anchoredFanTarget){lo*=2;vl=anchoredFanFast(lo);}if(vh>anchoredFanTarget){hi*=2;vh=anchoredFanFast(hi);}++anchoredInverseExpansions;}
    anchoredCheck(vl>=anchoredFanTarget&&vh<=anchoredFanTarget,"same32 inverse failed finite volume bracket");
    long double best=std::abs(vl-anchoredFanTarget)<std::abs(vh-anchoredFanTarget)?lo:hi,actual=anchoredFanFast(best);
    const auto minor=std::min(static_cast<long double>(anchoredFanTarget),static_cast<long double>(anchoredFanCapacity)-anchoredFanTarget),limit=anchoredOwnGate*minor;
    const auto resolved=[&]{return actual>0&&actual<anchoredFanCapacity&&std::abs(actual-anchoredFanTarget)<=limit;};
    for(unsigned trial=0;!resolved()&&trial<80;++trial){const auto mid=lo+(hi-lo)/2;if(mid==lo||mid==hi)break;++anchoredInverseIterations;const auto vm=anchoredFanFast(mid);
        if(std::abs(vm-anchoredFanTarget)<std::abs(actual-anchoredFanTarget)){best=mid;actual=vm;}if(vm==anchoredFanTarget)break;if(vm>anchoredFanTarget)lo=mid;else hi=mid;}
    anchoredCheck(resolved(),"same80 inverse did not satisfy actual own-minor gate");anchoredInverseInventory=anchoredInspectFan(best);
}
AnchoredSites anchoredBand(long double delta,bool axis=false)
{const Point base{axis?0.L:1.L,0};return {{base,{0,0}},{base,{delta,0}},{base,{delta,1}},{base,{0,1}}};}
void anchoredGeometryAndGuards()
{
    const auto delta=std::ldexp(1.L,-(std::numeric_limits<long double>::digits+16));
    auto band=anchoredBand(delta);
    const auto m=Cut::anchoredPolygonMoments(band);const auto ref=arMoments(arPoints(band));
    anchoredNear(m.area,arValue(ref.area),"sub-ULP actual band area");anchoredNear(m.firstX,arValue(ref.x),"sub-ULP actual band radial moment");anchoredNear(m.firstY,arValue(ref.y),"sub-ULP actual band axial moment");
    for(bool rz:{false,true})for(bool reverse:{false,true}){auto p=band;if(reverse)std::reverse(p.begin(),p.end());anchoredNear(Cut::anchoredPolygonMeasureInFrame(p,{},rz),arValue(arMeasure(arPoints(p),{},rz)),"supplied sub-ULP polygon measure");
        for(bool top:{false,true}){const auto clipped=Cut::clipAnchoredConvexPolygon(p,{0,top?-1.L:1.L},top?-.5L:.5L);anchoredNear(Cut::anchoredPolygonMeasureInFrame(clipped,{},rz),arValue(arMeasure(arHalfPlane(arPoints(p),top),{},rz)),"same-site thin plane clip");}++anchoredCoreCases;}
    const Cut::Rectangle full{{0,0},{2,1}},half{{0,0},{1,1}};
    anchoredCheck(Cut::anchoredPolygonInsideBox(band,{},full)&&!Cut::anchoredPolygonInsideBox(band,{},half),"exact containment lost coordinate tail");
    anchoredNear(Cut::anchoredPolygonMeasureInFrame(Cut::clipAnchoredPolygonInFrame(full,band,{}),{},false),delta,"full ordinary box actual clip");
    anchoredCheck(Cut::anchoredPolygonMeasureInFrame(Cut::clipAnchoredPolygonInFrame(half,band,{}),{},false)==0,"ordinary-box boundary strict zero");
    Polygon collapsed;for(auto p:band){const auto approximate=Cut::approximateAnchoredPoint(p);anchoredCheck(approximate.x==1,"explicit approximate display point did not collapse native sub-ULP tail");collapsed.push_back(approximate);}
    anchoredCheck(Cut::polygonMeasureInFrame(collapsed,{},false)==0,"old Point collapse control has invented positive area");
    const AnchoredSite q{{1,0},{delta,0}},equivalent{{1,delta},{delta,-delta}};anchoredCheck(Cut::sameAnchoredPoint(q,equivalent)&&arSame(q,equivalent),"equivalent supplied anchor gauges");
    const auto difference=Cut::anchoredDifference(q,{{1,0},{}});anchoredCheck(difference.x==delta&&difference.y==0,"exact site difference lost tail");
    const auto axis=anchoredBand(delta,true);anchoredNear(Cut::anchoredPolygonMeasureInFrame(axis,{},true),arValue(arMeasure(arPoints(axis),{},true)),"axis small phase actual radial integral");
    const std::array<AnchoredSite,3> tri{{{{0,0},{}},{{1,0},{}},{{0,1},{}}}};
    for(long double sign:{-1.L,1.L}){const auto pure=Cut::cutAnchoredTriangle(tri,{{sign,sign,sign}});anchoredCheck(pure.kind==(sign<0?Cut::TriangleRegion::Negative:Cut::TriangleRegion::Positive),"anchored pure kind");
        if(sign<0)anchoredNear(Cut::negativeAnchoredMeasure(pure,{},false),.5L,"anchored pure area");else anchoredCheck(Cut::negativeAnchoredMeasure(pure,{},true)==0,"anchored pure positive strict zero");
        const auto edge=Cut::cutAnchoredTriangle(tri,{{0,0,sign}});anchoredCheck(edge.zeroVertexMask==3&&edge.zeroEdgeMask==1&&edge.interface.has_value(),"anchored exact coincident edge tags");
        const auto vertex=Cut::cutAnchoredTriangle(tri,{{0,sign,sign}});anchoredCheck(vertex.zeroVertexMask==1&&!vertex.interface,"anchored exact isolated vertex tags");++anchoredCoreCases;}
    const auto cut=Cut::cutAnchoredTriangle(tri,{{0,-1,1}});anchoredCheck(cut.kind==Cut::TriangleRegion::Cut&&cut.zeroVertexMask==1&&cut.interface.has_value(),"anchored exact vertex interface");anchoredNear(Cut::negativeAnchoredMeasure(cut,{},false),.25L,"anchored zero vertex area");
    const auto zero=Cut::cutAnchoredTriangle(tri,{{0,0,0}});anchoredCheck(zero.kind==Cut::TriangleRegion::AllZeroScalar&&zero.zeroVertexMask==7&&zero.zeroEdgeMask==7&&!zero.interface,"anchored all-zero tags");
    anchoredReject([&]{Cut::negativeAnchoredMeasure(zero,{},false);},"anchored all-zero measure accepted");
    anchoredReject([&]{Cut::negativeAnchoredMeasure(tri,{{0,0,0}},{},false);},"anchored fast all-zero measure accepted");
    const auto nan=std::numeric_limits<long double>::quiet_NaN(),inf=std::numeric_limits<long double>::infinity();
    auto invalid=band;invalid[0].offset.x=nan;anchoredReject([&]{Cut::anchoredPolygonMoments(invalid);},"nonfinite actual offset accepted");invalid=band;invalid[0].anchor.y=inf;anchoredReject([&]{Cut::anchoredPolygonMoments(invalid);},"nonfinite actual anchor accepted");
    anchoredReject([&]{Cut::anchoredPolygonMeasureInFrame(band,{nan,0},false);},"nonfinite frame accepted");
    anchoredReject([&]{Cut::clipAnchoredConvexPolygon(band,{0,0},0);},"anchored zero plane normal accepted");
    anchoredReject([&]{Cut::cutAnchoredTriangle(tri,{{-1,inf,1}});},"anchored nonfinite scalar accepted");
    anchoredReject([&]{Cut::cutAnchoredTriangle({{tri[0],tri[0],tri[1]}},{{-1,0,1}});},"anchored repeated site accepted");
    anchoredReject([&]{Cut::cutAnchoredTriangle({{{{0,0},{}},{{1,0},{}},{{2,0},{}}}},{{-1,0,1}});},"anchored degenerate triangle accepted");
    const AnchoredSites crossed{{{0,0},{}},{{1,1},{}},{{0,1},{}},{{1,0},{}}},concave{{{0,0},{}},{{2,0},{}},{{.5L,.5L},{}},{{2,2},{}},{{0,2},{}}};
    anchoredReject([&]{Cut::anchoredPolygonMoments(crossed);},"anchored crossed polygon accepted");anchoredReject([&]{Cut::clipAnchoredConvexPolygon(concave,{1,0},.5L);},"anchored concave polygon accepted");
    const AnchoredSites negativeRadius{{{0,0},{-delta,0}},{{0,0},{delta,0}},{{0,1},{}}};anchoredReject([&]{Cut::anchoredPolygonMeasureInFrame(negativeRadius,{},true);},"negative physical radius hidden in tail accepted");
    anchoredReject([&]{Cut::negativeAnchoredMeasure(tri,{{1,1,1}},{-.125L,0},true);},"positive phase fast query skipped invalid physical radius");
    const std::array<AnchoredSite,3> mixedBadRadius{{{{0,0},{-delta,0}},{{0,0},{delta,0}},{{0,1},{}}}};
    anchoredReject([&]{Cut::negativeAnchoredMeasure(mixedBadRadius,{{1,-1,1}},{},true);},"positive source vertex invalid radius hidden by negative clip");
    anchoredReject([&]{Cut::anchoredZeroCrossing(tri[0],tri[1],1,1);},"nonopposite anchored crossing accepted");
    const auto tiny=std::numeric_limits<long double>::min();anchoredReject([&]{Cut::negativeAnchoredMeasure(tri,{{-tiny,1,1}},{},false);},"unsupported positive tiny measure silently erased");
    anchoredReject([&]{Cut::anchoredPolygonMeasureInFrame(AnchoredSites(tri.begin(),tri.end()),{std::numeric_limits<long double>::max()/2,0},true);},"unrepresentable positive RZ integral silently accepted");
}
void anchoredDoubleFrameContract()
{
    const auto before=assertions;
    const auto bits=[](double a,double b){return anchoredDoubleBits(a)==anchoredDoubleBits(b);};
    const auto result=[&](AnchoredSite query,AnchoredSite frame,Point center,Vector2D expected,const char* name){
        const auto actual=Cut::anchoredDifferenceDouble(query,frame,center);
        check(bits(actual.x,expected.x)&&bits(actual.y,expected.y),name);
    };
    const long double delta=std::ldexp(1.L,-160),height=std::ldexp(1.L,-80);
    const AnchoredSite query{{1,-1},{delta,-delta}},frame{{1,-1},{height,-height}};
    result(query,frame,{-height,height},{static_cast<double>(delta),-static_cast<double>(delta)},
        "complete frame contraction lost the low retained offset");
    const AnchoredSite querySwapped{query.offset,query.anchor},frameSwapped{frame.offset,frame.anchor};
    check(Cut::sameAnchoredPoint(query,querySwapped)&&Cut::sameAnchoredPoint(frame,frameSwapped),
        "frame equivalent-gauge fixture changed its actual sites");
    result(querySwapped,frameSwapped,{-height,height},{static_cast<double>(delta),-static_cast<double>(delta)},
        "different affine anchors changed complete binary64 contraction");
    const auto huge=std::numeric_limits<long double>::max();
    result({{huge,huge},{huge,-huge}},{{huge,0},{huge,0}},{-delta,-delta},
        {static_cast<double>(delta),static_cast<double>(delta)},
        "cancelled overflowing partial sums lost their finite complete result");
    const long double half=std::ldexp(1.L,-53),ulp=std::ldexp(1.L,-52);
    result({{1,-1},{half,-half}},{},{},{1,-1},"binary64 even-lower tie changed");
    result({{1+ulp,-(1+ulp)},{half,-half}},{},{},
        {1+2*static_cast<double>(ulp),-(1+2*static_cast<double>(ulp))},"binary64 even-upper tie changed");
    // The low third term decides the binary64 tie even when native LD first
    // rounding followed by a double cast would lose that decision.
    result({{1,-1},{half,-half}},{{0,0},{-delta,delta}},{},
        {1+static_cast<double>(ulp),-(1+static_cast<double>(ulp))},
        "complete binary64 contraction double-rounded its tie witness");
    result({{1,-1},{delta,-delta}},{{1,-1},{0,0}},{},
        {static_cast<double>(delta),-static_cast<double>(delta)},"single-subtraction frame path lost a finite offset");
    result({{1,-1},{-0.L,0.L}},{{1,-1},{0,0}},{},{0,0},"exact same-anchor cancellation is not canonical positive zero");
    result({{huge,huge},{-huge,-huge}},{},{},{0,0},"complete exact cancellation is not canonical positive zero");
    const double minimum=std::numeric_limits<double>::denorm_min(),maximum=std::numeric_limits<double>::max();
    result({{0,0},{minimum,-minimum}},{},{},{minimum,-minimum},"binary64 minimum subnormal was discarded");
    result({{1,1},{0,0}},{{1,1},{-minimum,minimum}},{},{minimum,-minimum},
        "complete multi-term frame discarded its binary64 subnormal");
    result({{maximum,-maximum},{delta,-delta}},{},{},{maximum,-maximum},"finite final extreme binary64 result rejected");
    if constexpr(std::numeric_limits<long double>::min_exponent-std::numeric_limits<long double>::digits<-1074){
        const auto halfMinimum=std::ldexp(1.L,-1075);
        check(halfMinimum>0,"wider-native half-subnormal fixture vanished");
        result({{0,0},{halfMinimum,-halfMinimum}},{},{},{0.,-0.},"binary64 half-subnormal tie did not round to signed zero");
        result({{0,0},{minimum,-minimum}},{{0,0},{-halfMinimum,halfMinimum}},{},
            {2*minimum,-2*minimum},"binary64 odd subnormal tie did not round to its even neighbor");
    }
    rejected([&]{Cut::anchoredDifferenceDouble({{maximum,-maximum},{maximum,-maximum}},{});},
        "final binary64 affine overflow was accepted");
    rejected([&]{Cut::anchoredDifferenceDouble({{1,1},{maximum,-maximum}},{{1,1},{}},{-maximum,maximum});},
        "single-subtraction binary64 overflow was accepted");
    for(const auto invalid:{std::numeric_limits<long double>::quiet_NaN(),
        std::numeric_limits<long double>::infinity(),-std::numeric_limits<long double>::infinity()}){
        for(unsigned component=0;component<10;++component){
            AnchoredSite q{},f{};Point center{};
            const std::array<long double*,10> fields{{&q.anchor.x,&q.anchor.y,&q.offset.x,&q.offset.y,
                &f.anchor.x,&f.anchor.y,&f.offset.x,&f.offset.y,&center.x,&center.y}};
            *fields[component]=invalid;
            rejected([&]{Cut::anchoredDifferenceDouble(q,f,center);},"nonfinite affine input was cancelled or ignored");
        }
    }
    {
        struct RestoreRounding {int mode=std::fegetround();~RestoreRounding(){std::fesetround(mode);}} restore;
        check(std::fesetround(FE_UPWARD)==0,"cannot set negative affine rounding fixture");
        rejected([&]{Cut::anchoredDifferenceDouble(query,frame,{-height,height});},
            "binary64 affine port accepted an invalid actual rounding mode");
    }
    std::cout<<"anchored_binary64_frame assertions="<<assertions-before<<" complete_RNE=PASS\n";
}
void anchoredCoreContract()
{
    anchoredOriginalFan();anchoredGeometryAndGuards();
    std::cout<<std::setprecision(17)<<"anchored_core cases="<<anchoredCoreCases<<" assertions="<<anchoredCoreAssertions<<" guards="<<anchoredCoreGuards
        <<" inverse_iterations="<<anchoredInverseIterations<<" inverse_expansions="<<anchoredInverseExpansions<<" peak_own_relative="<<anchoredCorePeak
        <<" fixed_inventory_relative="<<anchoredFixedInventory<<" inverse_inventory_relative="<<anchoredInverseInventory<<" point_collapse_relative="<<anchoredPointCollapse
        <<" delta="<<std::ldexp(1.L,-(std::numeric_limits<long double>::digits+16))<<" own_gate="<<anchoredOwnGate<<" PASS\n";
}
}
// END V4a1a ADDED ANCHORED TEST BLOCK
namespace
{
struct CornerAnalyticMoments { long double area=0,firstX=0,firstY=0; };
CornerAnalyticMoments cornerAnalytic(const Cut::Rectangle& cell,Cut::RectangleCorner corner,
    long double slope,bool horizontalNegative)
{
    const int ci=static_cast<int>(corner),sx=ci==0||ci==3?-1:1,sy=ci<2?-1:1;
    const long double w=static_cast<long double>(cell.upper.x)-cell.lower.x;
    const long double h=static_cast<long double>(cell.upper.y)-cell.lower.y;
    long double a=0,u=0,v=0;
    // Independent integrals of V <= slope*U in inward corner coordinates.
    // The two branches are a triangle and a rectangle minus a triangle.
    if (slope*w<=h) {
        a=slope*w*w/2;u=slope*w*w*w/3;v=slope*slope*w*w*w/6;
    } else {
        a=w*h-h*h/(2*slope);u=w*w*h/2-h*h*h/(6*slope*slope);
        v=w*h*h/2-h*h*h/(3*slope);
    }
    if (!horizontalNegative) {a=w*h-a;u=w*w*h/2-u;v=w*h*h/2-v;}
    return {a,sx*w*a/2-sx*u,sy*h*a/2-sy*v};
}
void cornerInverseContract()
{
    const auto before=assertions;unsigned cases=0;
    const std::array<Cut::Rectangle,5> boxes{{{{1,2},{3,6}},{{0,0},{2,4}},
        {{2,-1},{2.25,7}},{{8,-1},{40,0}},{{1024,-4096},{1026,-4092}}}};
    const std::array<Cut::RectangleCorner,4> corners{{Cut::RectangleCorner::LowerLeft,
        Cut::RectangleCorner::LowerRight,Cut::RectangleCorner::UpperRight,Cut::RectangleCorner::UpperLeft}};
    const auto verifyPlane=[&](const Cut::Rectangle& box,Cut::RectangleCorner corner,bool rz,
        bool branch,long double fraction) {
        ++cases;
        const auto result=Cut::reconstructThroughCorner(box,corner,fraction,rz,branch);
        const auto& plane=result.plane;const int ci=static_cast<int>(corner);
        const int sx=ci==0||ci==3?-1:1,sy=ci<2?-1:1;
        const long double hx=(static_cast<long double>(box.upper.x)-box.lower.x)/2;
        const long double hy=(static_cast<long double>(box.upper.y)-box.lower.y)/2;
        const auto c=Cut::capacity(box,rz),actual=Cut::volume(box,plane.normal,plane.offset,rz);
        check(result.corner==corner&&result.localCorner.x==sx*hx&&result.localCorner.y==sy*hy,
            "corner inverse changed the actual local corner identity");
        check(actual==plane.volume&&plane.residual==actual-fraction*c,
            "corner inverse returned a different plane volume/residual");
        near(actual,fraction*c,c,"corner inverse original volume gate",Cut::relativeVolumeResolution);
        const auto length=std::hypot(static_cast<long double>(plane.normal.x),static_cast<long double>(plane.normal.y));
        const Point n{plane.normal.x/length,plane.normal.y/length};
        check(n.x*result.localCorner.x+n.y*result.localCorner.y-plane.offset==0,
            "corner inverse lost its exact normalized corner zero");
        const auto segment=Cut::interfaceSegment(box,plane.normal,plane.offset);
        check(segment.has_value(),"corner inverse has no actual segment/edge");
        check(((*segment)[0].x==result.localCorner.x&&(*segment)[0].y==result.localCorner.y)||
            ((*segment)[1].x==result.localCorner.x&&(*segment)[1].y==result.localCorner.y),
            "original interfaceSegment did not retain the exact local corner");
        if (fraction>0&&fraction<1) {
            check(actual>0&&Cut::volume(box,{-plane.normal.x,-plane.normal.y},-plane.offset,rz)>0,
                "corner inverse erased a positive geometric phase");
            const auto horizontal=n.x*(-sx*hx)+n.y*(sy*hy)-plane.offset;
            const auto vertical=n.x*(sx*hx)+n.y*(-sy*hy)-plane.offset;
            check(branch?(horizontal<0&&vertical>0):(horizontal>0&&vertical<0),
                "corner inverse changed the declared adjacent phase sides");
        }
        return result;
    };
    for (const auto& box:boxes) for (bool rz:{false,true}) for (auto corner:corners)
        for (bool branch:{false,true}) for (int shapeKind:{0,1}) {
            const long double w=static_cast<long double>(box.upper.x)-box.lower.x;
            const long double h=static_cast<long double>(box.upper.y)-box.lower.y;
            const long double slope=shapeKind?2*h/w:h/(2*w);
            const auto expected=cornerAnalytic(box,corner,slope,branch);
            const long double radius=static_cast<long double>(box.lower.x)+w/2;
            const auto expectedVolume=rz?2*pi*(radius*expected.area+expected.firstX):expected.area;
            const auto fraction=expectedVolume/Cut::capacity(box,rz);
            const auto result=verifyPlane(box,corner,rz,branch,fraction);
            const auto& p=result.plane;const auto length=std::hypot(static_cast<long double>(p.normal.x),static_cast<long double>(p.normal.y));
            const auto polygon=Cut::clipConvexPolygon(rectangle(-w/2,w/2,-h/2,h/2),
                {p.normal.x/length,p.normal.y/length},p.offset);
            const auto actual=Cut::polygonMoments(polygon);
            near(actual.area,expected.area,w*h,"corner analytic triangle/trapezoid area");
            near(actual.firstX,expected.firstX,w*w*h,"corner analytic first x moment");
            near(actual.firstY,expected.firstY,w*h*h,"corner analytic first y moment");
            near(p.volume,expectedVolume,Cut::capacity(box,rz),"corner analytic physical XY/RZ volume");
        }
    const Cut::Rectangle nearPureBox{{1,2},{3,6}};
    const long double small=std::ldexp(1.L,-40);
    for (bool rz:{false,true}) for (auto corner:corners) for (bool branch:{false,true}) {
        for (long double fraction:{0.L,1.L,small,1-small,std::nextafter(1.L,0.L)})
            verifyPlane(nearPureBox,corner,rz,branch,fraction);
    }
    const Cut::Rectangle valid{{1,2},{3,6}};
    const auto lowerLeft=Cut::RectangleCorner::LowerLeft;
    for (long double invalid:{-1.L,2.L,std::numeric_limits<long double>::quiet_NaN(),
        std::numeric_limits<long double>::infinity(),-std::numeric_limits<long double>::infinity()})
        rejected([&]{Cut::reconstructThroughCorner(valid,lowerLeft,invalid,false);},"corner inverse accepted an invalid fraction");
    rejected([&]{Cut::reconstructThroughCorner(valid,static_cast<Cut::RectangleCorner>(4),.5L,false);},
        "corner inverse accepted an invalid corner");
    rejected([&]{Cut::reconstructThroughCorner({{1,2},{1,6}},lowerLeft,.5L,false);},
        "corner inverse accepted an empty rectangle");
    rejected([&]{Cut::reconstructThroughCorner({{-1,2},{3,6}},lowerLeft,.5L,true);},
        "corner inverse accepted a negative physical radius");
    rejected([&]{Cut::reconstructThroughCorner({{1,2},{std::numeric_limits<double>::infinity(),6}},lowerLeft,.5L,false);},
        "corner inverse accepted a nonfinite rectangle");
    const auto tiny=std::numeric_limits<long double>::denorm_min();
    if (tiny>0) {
        check(tiny*.25L==0,"corner inverse underflow fixture is not unrepresentable");
        rejected([&]{Cut::reconstructThroughCorner({{0,0},{.5,.5}},lowerLeft,tiny,false);},
            "corner inverse silently erased an unrepresentable positive mixed target");
    }
    std::cout<<"corner_inverse cases="<<cases<<" assertions="<<assertions-before
        <<" local_corner_exact=PASS analytic_XY_RZ=PASS\n";
}
}
int main(int argc,char** argv)
{
    try {
        if (argc==2&&std::string(argv[1])=="--corner-inverse-only") {
            cornerInverseContract();
            std::cout<<"planar_cell_cut_contract corner_inverse PASS\n";return 0;
        }
        cornerInverseContract();
        moments();triangles();amrPatch();inverseAndLegacy();invalidInput();
        thinBinaryPolygons();
        // BEGIN V4a1a ADDED CORE MAIN CALL
        anchoredCoreContract();
        anchoredDoubleFrameContract();
        // END V4a1a ADDED CORE MAIN CALL
        std::cout<<std::setprecision(17)<<"planar_cell_cut_contract PASS assertions="<<assertions
            <<" peak_scaled_error="<<peakScaledError<<" tolerance="<<tolerance<<'\n';
        return 0;
    } catch (const std::exception& error) { std::cerr<<error.what()<<'\n';return 1; }
}
