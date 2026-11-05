#include "Modules/Hydrodynamics/HydroBalancedCsfLoad.h"
#include "Modules/Hydrodynamics/HydroSpatialAdvector.h"
#include "Modules/Interface/HermiteGalsAmrClosure.h"
#include "Modules/Hydrodynamics/Detail/HydroRegularAdvectorGeometry.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/UniformLeafGrid.h"
#include "Core/Mesh/CartesianLeafLocator.h"
#include "Core/Mesh/DyadicGridCoordinates.h"
#include "Core/Geometry/Detail/ExactFaceTimeRational.h"
#include "Core/Geometry/Detail/FiniteVolumeClosure.h"
#include "Core/Geometry/PolynomialVelocitySlab.h"
#include "Core/Geometry/CompactQuarticDerivative.h"
#include "Core/Geometry/PlanarCellCut.h"
#include <Eigen/QR>
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <mutex>
#include <set>
#include <numeric>
#include <optional>
#include <tuple>
#include <sstream>
#include <iomanip>
#include <stdexcept>
#include <string>
#include <utility>

namespace LAMP::Hydrodynamics
{
namespace
{
void require(bool value,const char* reason)
{if(!value)throw std::invalid_argument(std::string("[HydroBalancedCsfLoad] ")+reason);}
void validateCurvatureConfig(const HydroCsfCurvatureConfig& cfg)
{
    require(cfg.forceMethod==HydroCsfForceMethod::MotionAdjoint||cfg.forceMethod==HydroCsfForceMethod::HeightFunction
        ||cfg.forceMethod==HydroCsfForceMethod::HeightFunctionIncrement,"invalid CSF force method");
    require(cfg.maxCandidateQueriesPerFace>0&&cfg.maxVisibilityCellQueriesPerFace>0&&cfg.maxSurfaceQueriesPerFace>0&&cfg.maxGeometryQueriesPerFace>0&&cfg.maxShadowWallQueriesPerFace>0&&
        cfg.maxQuadratureDepth>0&&cfg.maxQuadratureDepth<=40&&std::isfinite(cfg.relativeQuadratureTolerance)&&cfg.relativeQuadratureTolerance>0&&cfg.relativeQuadratureTolerance<1&&
        std::isfinite(cfg.relativeGradientRoundoff)&&cfg.relativeGradientRoundoff>0&&cfg.relativeGradientRoundoff<1,
        "invalid CSF curvature config");
}
struct Box {Vector2D lo,hi;};
namespace Exact=LAMP::Geometry::Detail::ExactExpansion;
using RayParameter=LAMP::Geometry::Detail::FaceTime::FaceTimeRational;
// Only the exceptional finite-inventory height certificate uses rational
// arithmetic. It bounds an unknown minor geometry; it never creates a plane,
// changes either stock, or changes the ordinary Kahan height path.
namespace HeightFinite=LAMP::Geometry::Detail::FiniteVolume;
using HeightInteger=boost::multiprecision::cpp_int;
using HeightRational=boost::multiprecision::cpp_rational;
HeightRational exactHeightNative(long double value)
{
    require(std::isfinite(value),"nonfinite finite-height factor");
    if(value==0)return 0;
    int exponent=0;auto fraction=std::frexp(std::abs(value),&exponent);
    HeightInteger integer=0;int remaining=std::numeric_limits<long double>::digits;
    while(remaining){const int take=std::min(remaining,32);fraction=std::scalbn(fraction,take);
        const auto chunk=static_cast<std::uint32_t>(fraction);integer<<=take;integer+=chunk;
        fraction-=chunk;remaining-=take;}
    exponent-=std::numeric_limits<long double>::digits;
    HeightRational result=integer;
    if(exponent>=0)result*=HeightInteger(1)<<exponent;else result/=HeightInteger(1)<<(-exponent);
    if(std::signbit(value))result=-result;return result;
}
// One terminal nearest-even rounding of an exact rational. In particular,
// sub-eta contributions are retained until the integer remainder comparison.
bool roundedHeightNative(const HeightRational& value,long double& result)
{
    if(value==0){result=0;return true;}
    const bool negative=value<0;
    HeightInteger numerator=boost::multiprecision::numerator(value);
    if(negative)numerator=-numerator;
    const HeightInteger denominator=boost::multiprecision::denominator(value);
    int top=static_cast<int>(boost::multiprecision::msb(numerator))-static_cast<int>(boost::multiprecision::msb(denominator));
    if(top>=0?numerator<(denominator<<top):(numerator<<(-top))<denominator)--top;
    constexpr int digits=std::numeric_limits<long double>::digits;
    if(top>=std::numeric_limits<long double>::max_exponent)return false;
    const int quantum=std::max(top-(digits-1),std::numeric_limits<long double>::min_exponent-digits);
    HeightInteger scaledNumerator=numerator,scaledDenominator=denominator;
    if(quantum>=0)scaledDenominator<<=quantum;else scaledNumerator<<=-quantum;
    HeightInteger integer=scaledNumerator/scaledDenominator;
    const HeightInteger remainder=scaledNumerator%scaledDenominator,twice=remainder<<1;
    if(twice>scaledDenominator||(twice==scaledDenominator&&(integer&1)!=0))++integer;
    if(integer==0){result=std::copysign(0.L,negative?-1.L:1.L);return true;}
    if(static_cast<int>(boost::multiprecision::msb(integer))+quantum>=std::numeric_limits<long double>::max_exponent)return false;
    const auto magnitude=std::scalbn(integer.convert_to<long double>(),quantum);
    if(!std::isfinite(magnitude)||magnitude<=0)return false;
    result=negative?-magnitude:magnitude;return true;
}
bool sameHeightNative(long double a,long double b)
{return a==b&&(a!=0||std::signbit(a)==std::signbit(b));}
struct FiniteHeightMinor
{
    bool liquid=false;
    HeightRational lower=0,upper=0,area=0;
};
template<class Cell>
std::optional<FiniteHeightMinor> finiteHeightMinor(const Cell& cell,
    const Geometry::PlanarCellCut::Rectangle& clip,bool rz,std::string& reason)
{
    const auto& rectangle=cell.rectangle;
    const bool liquid=cell.inventory[0]<cell.inventory[1];const int minor=liquid?0:1,major=1-minor;
    const long double hx=(static_cast<long double>(rectangle.upper.x)-rectangle.lower.x)/2;
    const auto capacity=Geometry::PlanarCellCut::capacity(rectangle,rz);
    if(cell.resolved||!(cell.inventory[minor]>0&&cell.inventory[major]>cell.inventory[minor])||
        cell.normal.x!=1||cell.normal.y!=0||cell.offset!=(major==0?hx:-hx)||
        cell.measured[minor]!=0||cell.measured[major]!=capacity){
        reason="finite-height support is not a published unresolved pure boundary-placeholder plane";return std::nullopt;}
    if(rz&&clip.lower.x<=0){
        reason="finite-height support has no positive actual radial lower bound";return std::nullopt;}
    FiniteHeightMinor result;result.liquid=liquid;
    result.area=(exactHeightNative(clip.upper.x)-exactHeightNative(clip.lower.x))*
        (exactHeightNative(clip.upper.y)-exactHeightNative(clip.lower.y));
    const bool whole=clip.lower.x==rectangle.lower.x&&clip.lower.y==rectangle.lower.y&&
        clip.upper.x==rectangle.upper.x&&clip.upper.y==rectangle.upper.y;
    const auto stock=exactHeightNative(cell.inventory[minor]);
    if(rz){const auto factor=exactHeightNative(CartesianMeasure::twoPi);
        result.upper=stock/(factor*exactHeightNative(clip.lower.x));
        if(whole)result.lower=stock/(factor*exactHeightNative(clip.upper.x));}
    else{result.upper=stock;if(whole)result.lower=stock;}
    // Intersection with the actual query area is a geometric upper bound,
    // not a changed source stock or a rounded minor-area estimate.
    if(result.upper>result.area)result.upper=result.area;
    if(result.area<=0||result.lower<0||result.upper<result.lower){
        reason="finite-height stock bounds do not intersect the actual query area";return std::nullopt;}
    return result;
}
int compare(const RayParameter& a,const RayParameter& b)
{return Exact::sign(Exact::subtract(Exact::multiply(a.n,b.d),Exact::multiply(b.n,a.d)));}
struct RayLine
{
    Vector2D begin;
    std::array<Exact::Expansion,2> delta;
    RayLine(Vector2D a,Vector2D b):begin(a),delta{{Exact::subtract(Exact::fromDouble(b.x),Exact::fromDouble(a.x)),
        Exact::subtract(Exact::fromDouble(b.y),Exact::fromDouble(a.y))}}{}
    bool clip(Box box,RayParameter& low,RayParameter& high) const
    {
        for(int axis=0;axis<2;++axis){const double start=axis?begin.y:begin.x,lo=axis?box.lo.y:box.lo.x,hi=axis?box.hi.y:box.hi.x;
            if(Exact::sign(delta[axis])==0){if(start<lo||start>hi)return false;}
            else{RayParameter a(Exact::subtract(Exact::fromDouble(lo),Exact::fromDouble(start)),delta[axis]),
                    b(Exact::subtract(Exact::fromDouble(hi),Exact::fromDouble(start)),delta[axis]);
                if(compare(a,b)>0)std::swap(a,b);if(compare(a,low)>0)low=std::move(a);if(compare(b,high)<0)high=std::move(b);}
            if(compare(low,high)>0)return false;}
        return true;
    }
};
using RayIntervals=std::vector<std::pair<RayParameter,RayParameter>>;
bool covered(RayIntervals& intervals)
{
    std::sort(intervals.begin(),intervals.end(),[](const auto& a,const auto& b){const int order=compare(a.first,b.first);
        return order<0||(order==0&&compare(a.second,b.second)<0);});
    RayParameter end(0),one(1);
    for(const auto& span:intervals){if(compare(span.first,end)>0)return false;if(compare(span.second,end)>0)end=span.second;
        if(compare(end,one)>=0)return true;}
    return false;
}
Box surrounding(Vector2D p,double radius)
{
    const double inf=std::numeric_limits<double>::infinity();
    return {{std::nextafter(p.x-radius,-inf),std::nextafter(p.y-radius,-inf)},
            {std::nextafter(p.x+radius,inf),std::nextafter(p.y+radius,inf)}};
}
// Bins only accelerate the fixed physical support. They do not truncate it
// to a ring count or alter either the native mesh or the source quadrature.
struct Buckets
{
    using Key=std::pair<long long,long long>;
    Vector2D origin{0,0};double width=1;
    std::map<Key,std::vector<int>> data;
    Key key(Vector2D p) const
    {
        const double x=std::floor((p.x-origin.x)/width),y=std::floor((p.y-origin.y)/width);
        require(std::isfinite(x)&&std::isfinite(y)&&std::abs(x)<1e15&&std::abs(y)<1e15,"unrepresentable CSF spatial index");
        return {static_cast<long long>(x),static_cast<long long>(y)};
    }
    template<class F>void visit(Box b,F&& action) const
    {
        const auto lo=key(b.lo),hi=key(b.hi);
        require(hi.first>=lo.first&&hi.second>=lo.second&&hi.first-lo.first<=8&&hi.second-lo.second<=8,"unresolved CSF local index extent");
        for(auto x=lo.first;x<=hi.first;++x)for(auto y=lo.second;y<=hi.second;++y){
            const auto it=data.find({x,y});if(it!=data.end())for(int id:it->second)action(id);}
    }
    void insert(Box b,int id)
    {
        const auto lo=key(b.lo),hi=key(b.hi);
        require(hi.first-lo.first<=2&&hi.second-lo.second<=2,"unresolved CSF index insertion");
        for(auto x=lo.first;x<=hi.first;++x)for(auto y=lo.second;y<=hi.second;++y)data[{x,y}].push_back(id);
    }
};
struct SurfacePoint {Vector2D point;double weight,curvature;int pressure,level;};
struct Range {double lo,hi;};
Range embedded(double lo,double hi,double fraction)
{
    const auto range=Exact::enclosure(Exact::add(Exact::fromDouble(lo),Exact::multiply(
        Exact::subtract(Exact::fromDouble(hi),Exact::fromDouble(lo)),Exact::fromDouble(fraction))));
    return {range.lo,range.hi};
}
struct SurfaceGraph {std::size_t component,graph;int pressure,level,cell;bool coordinateIsX;Box bounds;};
struct ProjectionIndex
{
    std::vector<Box> cells,walls;
    std::vector<double> diagonal;
    std::vector<int> component,level;
    std::vector<NativeLeafFace> faces;
    std::vector<SurfaceGraph> graphs;
    std::vector<int> wallAxis,wallDirection,wallPressure;
    Buckets cellBins,wallBins,graphBins;
    const HydroRegularizedColor& color;
    std::map<std::tuple<int,double,double>,std::vector<SurfacePoint>> sampleCache;
    std::map<std::pair<int,double>,Vector2D> endpointSamples;
    std::map<int,std::array<Range,2>> endpointBounds;
    std::vector<int> scratch;
    std::vector<std::size_t> cellVisitStamp;
    std::vector<std::size_t> graphVisitStamp,graphPositiveStamp,wallVisitStamp;
    std::size_t visitEpoch=0;
    std::size_t graphVisitEpoch=0,wallVisitEpoch=0;
    RayIntervals intervals;
    Box selectedBounds;bool convexRectangle=false;
    ProjectionIndex(const QuadtreeMesh& mesh,const HydroRegularizedColor& inputColor,HydroCsfCurvatureReport& report):color(inputColor)
    {
        const auto& phase=color.phaseGeometry();const auto& l=phase.layout();const int np=l.pressureCount();
        require(np>0,"empty CSF pressure domain");cellVisitStamp.resize(np,0);
        const auto& scalar=phase.publication().state().geometry().snapshot();
        Vector2D origin{std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity()};double maxDiagonal=0;
        for(int p=0;p<np;++p){const int cell=l.nativeCell(p);const auto b=scalar.nativeCellBounds(cell);
            cells.push_back({b[0],b[1]});const auto h=b[1]-b[0];const double d=std::hypot(h.x,h.y);
            require(std::isfinite(d)&&d>0,"unrepresentable CSF leaf diameter");diagonal.push_back(d);maxDiagonal=std::max(maxDiagonal,d);
            origin.x=std::min(origin.x,b[0].x);origin.y=std::min(origin.y,b[0].y);level.push_back(mesh.cells(mesh.cellIndexToLabel[cell]).level());}
        const double width=color.config().epsilon+2*maxDiagonal;
        require(std::isfinite(width)&&width>0,"unrepresentable CSF support radius");
        // Smaller bins reduce exact visibility leaf visits while retaining
        // every source candidate. The largest support spans at most six bins
        // per axis; native leaf boxes still fit inside one diagonal width.
        const double binWidth=std::max(maxDiagonal,width/3);
        for(auto* bins:{&cellBins,&wallBins,&graphBins}){bins->origin=origin;bins->width=binWidth;}
        for(int p=0;p<np;++p)cellBins.insert(cells[p],p);
        selectedBounds={origin,origin};auto area=Exact::fromDouble(0);
        for(const auto& box:cells){selectedBounds.hi.x=std::max(selectedBounds.hi.x,box.hi.x);selectedBounds.hi.y=std::max(selectedBounds.hi.y,box.hi.y);
            area=Exact::add(area,Exact::multiply(Exact::subtract(Exact::fromDouble(box.hi.x),Exact::fromDouble(box.lo.x)),
                Exact::subtract(Exact::fromDouble(box.hi.y),Exact::fromDouble(box.lo.y))));}
        // Native selected leaves are disjoint. Exact area equality proves
        // this particular union is a complete closed rectangle, allowing
        // a convex-domain shortcut. Arbitrary public boxes use exact slabs.
        convexRectangle=Exact::sign(Exact::subtract(area,Exact::multiply(
            Exact::subtract(Exact::fromDouble(selectedBounds.hi.x),Exact::fromDouble(origin.x)),
            Exact::subtract(Exact::fromDouble(selectedBounds.hi.y),Exact::fromDouble(origin.y)))))==0;
        std::vector<int> parent(np);std::iota(parent.begin(),parent.end(),0);
        const auto root=[&](int p){int r=p;while(parent[r]!=r)r=parent[r];while(parent[p]!=p){const int next=parent[p];parent[p]=r;p=next;}return r;};
        for(int u=0;u<l.velocityCount();++u){const auto f=l.domain().face(mesh,l.nativeFace(u));faces.push_back(f);
            const int a=f.negativeCell<0?-1:l.pressureDof(f.negativeCell),b=f.positiveCell<0?-1:l.pressureDof(f.positiveCell);
            if(a>=0&&b>=0&&f.area>0)parent[root(b)]=root(a);
            if(a<0||b<0){const auto& e=f.horizontal?mesh.hEdges(f.edgeLabel):mesh.vEdges(f.edgeLabel);
                const auto x=mesh.verts(e.edgeToVert[0]).coord(),y=mesh.verts(e.edgeToVert[1]).coord();
                walls.push_back({{std::min(x.x,y.x),std::min(x.y,y.y)},{std::max(x.x,y.x),std::max(x.y,y.y)}});
                wallAxis.push_back(f.horizontal?1:0);wallDirection.push_back(a>=0?1:-1);wallPressure.push_back(a>=0?a:b);
                wallBins.insert(walls.back(),static_cast<int>(walls.size()-1));}}
        std::map<int,int> numbered;
        for(int p=0;p<np;++p){const int r=root(p);if(!numbered.count(r))numbered[r]=static_cast<int>(numbered.size());component.push_back(numbered[r]);}
        const auto& publication=phase.publication();
        const auto clock=color.timeLayer();
        for(std::size_t c=0;c<publication.topology().curves.size();++c)for(std::size_t g=0;g<publication.topology().curves[c].graphs().size();++g){
            const auto& original=publication.topology().curves[c].graphs()[g];
            const int gp=l.pressureDof(original.nativeCellIndex);
            if(gp>=0){const auto& leaf=cells[gp];Box graphBox;
                const auto x0=embedded(leaf.lo.x,leaf.hi.x,original.localLower.x),x1=embedded(leaf.lo.x,leaf.hi.x,original.localUpper.x),
                    y0=embedded(leaf.lo.y,leaf.hi.y,original.localLower.y),y1=embedded(leaf.lo.y,leaf.hi.y,original.localUpper.y);
                graphBox={{x0.lo,y0.lo},{x1.hi,y1.hi}};
                if(original.nativeSide){const auto a=original.nativeSide->begin(),b=original.nativeSide->end();
                    require(original.nativeSide->usesSource(scalar)&&original.nativeSide->ownerCell()==original.nativeCellIndex,
                        "different native side source/cell");
                    graphBox={{std::min(a.x,b.x),std::min(a.y,b.y)},{std::max(a.x,b.x),std::max(a.y,b.y)}};}
                graphs.push_back({c,g,gp,level[gp],original.nativeCellIndex,original.coordinateIsX,graphBox});
                graphBins.insert(graphBox,static_cast<int>(graphs.size()-1));}
            for(const auto& q:publication.quadrature(c,g)){
                ++report.sourceNodes;const auto& t=q.sample.timeLayer;
                require(q.sample.geometryIdentity==color.geometryIdentity()&&t.step==clock.step&&t.begin==clock.begin&&t.end==clock.end&&t.position==clock.position,
                    "different surface geometry/full clock");
                require(std::isfinite(q.physicalWeight)&&q.physicalWeight>0&&std::isfinite(q.sample.curvature)&&
                    std::isfinite(q.sample.frame.point.x)&&std::isfinite(q.sample.frame.point.y),"invalid original CSF surface sample");
                const int p=l.pressureDof(q.sample.frame.nativeCellIndex);if(p<0)continue;
                ++report.selectedNodes;}}
        graphVisitStamp.resize(graphs.size(),0);graphPositiveStamp.resize(graphs.size(),0);wallVisitStamp.resize(walls.size(),0);
    }
    bool nearWall(Vector2D p,double radius) const
    {
        bool found=false;wallBins.visit(surrounding(p,radius),[&](int id){
            const auto& b=walls[id];const double dx=std::max({0.,b.lo.x-p.x,p.x-b.hi.x}),dy=std::max({0.,b.lo.y-p.y,p.y-b.hi.y});
            found=found||std::hypot(dx,dy)<radius;});return found;
    }
    bool visible(Vector2D a,Vector2D b,int pressureComponent,HydroCsfFaceSupportReport& report,std::size_t queryLimit)
    {
        const auto inside=[&](Vector2D p){return p.x>=selectedBounds.lo.x&&p.x<=selectedBounds.hi.x&&p.y>=selectedBounds.lo.y&&p.y<=selectedBounds.hi.y;};
        if(convexRectangle&&inside(a)&&inside(b))return true;
        scratch.clear();intervals.clear();
        if(++visitEpoch==0){std::fill(cellVisitStamp.begin(),cellVisitStamp.end(),0);visitEpoch=1;}
        cellBins.visit({{std::min(a.x,b.x),std::min(a.y,b.y)},{std::max(a.x,b.x),std::max(a.y,b.y)}},[&](int id){
            if(cellVisitStamp[id]==visitEpoch)return;
            require(report.visibilityCellQueries<queryLimit,"CSF visibility cell query limit");++report.visibilityCellQueries;
            cellVisitStamp[id]=visitEpoch;scratch.push_back(id);});
        const RayLine ray(a,b);
        for(int id:scratch){
            if(component[id]!=pressureComponent)continue;RayParameter low(0),high(1);
            if(ray.clip(cells[id],low,high))intervals.emplace_back(std::move(low),std::move(high));}
        return covered(intervals);
    }
    std::vector<SurfacePoint> samples(int graph,double begin,double end,HydroCsfFaceSupportReport& report,const HydroCsfCurvatureConfig& cfg)
    {
        require(report.surfaceQueries<=cfg.maxSurfaceQueriesPerFace&&cfg.maxSurfaceQueriesPerFace-report.surfaceQueries>=8,"CSF surface query limit");
        report.surfaceQueries+=8;
        const auto key=std::make_tuple(graph,begin,end);const auto found=sampleCache.find(key);
        if(found!=sampleCache.end())return found->second;
        const auto& g=graphs[graph];const auto& pub=color.phaseGeometry().publication();std::vector<SurfacePoint> result;
        for(const auto& q:pub.quadrature(g.component,g.graph,begin,end)){
            const auto& t=q.sample.timeLayer;const auto clock=color.timeLayer();
            require(q.sample.geometryIdentity==color.geometryIdentity()&&t.step==clock.step&&t.begin==clock.begin&&t.end==clock.end&&t.position==clock.position,
                "different subarc geometry/full clock");
            require(q.sample.frame.nativeCellIndex==g.cell&&q.physicalWeight>0&&std::isfinite(q.physicalWeight)&&std::isfinite(q.sample.curvature),"invalid original subarc sample");
            result.push_back({q.sample.frame.point,q.physicalWeight,q.sample.curvature,g.pressure,g.level});}
        require(result.size()==8,"unexpected original subarc rule");
        if(sampleCache.size()>=16384)sampleCache.clear(); // bounded, disposable source-query acceleration
        sampleCache.emplace(key,result);return result;
    }
    Vector2D endpoint(int id,double fraction,HydroCsfFaceSupportReport& report,const HydroCsfCurvatureConfig& cfg)
    {
        require(report.surfaceQueries<cfg.maxSurfaceQueriesPerFace,"CSF surface query limit");++report.surfaceQueries;
        const auto key=std::make_pair(id,fraction);const auto found=endpointSamples.find(key);if(found!=endpointSamples.end())return found->second;
        const auto& graph=graphs[id];const auto q=color.phaseGeometry().publication().sample(graph.component,graph.graph,fraction);
        const auto t=q.timeLayer,clock=color.timeLayer();
        require(q.geometryIdentity==color.geometryIdentity()&&t.step==clock.step&&t.begin==clock.begin&&t.end==clock.end&&t.position==clock.position&&
            q.frame.nativeCellIndex==graph.cell,"different endpoint geometry/full clock/cell");
        require(std::isfinite(q.frame.point.x)&&std::isfinite(q.frame.point.y),"nonfinite original endpoint");
        if(endpointSamples.size()>=16384)endpointSamples.clear();endpointSamples.emplace(key,q.frame.point);return q.frame.point;
    }
    static bool outsideDisk(Box box,Vector2D center,double radius)
    {
        const auto difference=[](double lo,double hi,double c){return c<lo?Exact::subtract(Exact::fromDouble(lo),Exact::fromDouble(c)):
            c>hi?Exact::subtract(Exact::fromDouble(c),Exact::fromDouble(hi)):Exact::fromDouble(0);};
        const auto x=difference(box.lo.x,box.hi.x,center.x),y=difference(box.lo.y,box.hi.y,center.y),r=Exact::fromDouble(radius);
        return Exact::sign(Exact::subtract(Exact::add(Exact::multiply(x,x),Exact::multiply(y,y)),Exact::multiply(r,r)))>=0;
    }
    bool shadowed(Box box,Vector2D from,int pressureComponent,HydroCsfFaceSupportReport& report,const HydroCsfCurvatureConfig& cfg)
    {
        if(convexRectangle)return false;
        if(++wallVisitEpoch==0){std::fill(wallVisitStamp.begin(),wallVisitStamp.end(),0);wallVisitEpoch=1;}
        bool found=false;
        wallBins.visit({{std::min(box.lo.x,from.x),std::min(box.lo.y,from.y)},
            {std::max(box.hi.x,from.x),std::max(box.hi.y,from.y)}},[&](int i){
            if(found||wallVisitStamp[i]==wallVisitEpoch)return;wallVisitStamp[i]=wallVisitEpoch;
            require(report.shadowWallQueries<cfg.maxShadowWallQueriesPerFace,"CSF shadow wall query limit");++report.shadowWallQueries;
            if(component[wallPressure[i]]!=pressureComponent)return;
            const int axis=wallAxis[i],direction=wallDirection[i];const auto& wall=walls[i];
            const double k=axis?wall.lo.y:wall.lo.x,fk=axis?from.y:from.x,ft=axis?from.x:from.y,
                lo=axis?box.lo.y:box.lo.x,hi=axis?box.hi.y:box.hi.x,
                wlo=axis?wall.lo.x:wall.lo.y,whi=axis?wall.hi.x:wall.hi.y;
            if(direction>0?(fk>k||lo<=k):(fk<k||hi>=k))return;
            if(fk==k){found=ft>wlo&&ft<whi;return;}
            bool all=true;
            for(double x:{box.lo.x,box.hi.x})for(double y:{box.lo.y,box.hi.y}){
                const double qk=axis?y:x,qt=axis?x:y;
                const auto den=Exact::subtract(Exact::fromDouble(qk),Exact::fromDouble(fk));
                const RayParameter hit(Exact::add(Exact::multiply(den,Exact::fromDouble(ft)),
                    Exact::multiply(Exact::subtract(Exact::fromDouble(k),Exact::fromDouble(fk)),Exact::subtract(Exact::fromDouble(qt),Exact::fromDouble(ft)))),den);
                all=all&&compare(hit,RayParameter(wlo))>0&&compare(hit,RayParameter(whi))<0;}
            found=all;
        });
        return found;
    }
    bool noZeroSet(const SurfaceGraph& graph,Box box,HydroCsfFaceSupportReport& report,const HydroCsfCurvatureConfig& cfg) const
    {
        require(report.geometryQueries<cfg.maxGeometryQueriesPerFace,"CSF geometry query limit");++report.geometryQueries;
        const auto& scalar=color.phaseGeometry().publication().state().geometry().snapshot();
        if(scalar.usesC2Macro()) {
            // This graph is inside its original native leaf. Intersect the
            // search enclosure with that exact support before taking bounds.
            const auto b=scalar.nativeCellBounds(graph.cell);
            Vector2D lower,upper;
            for(int axis=0;axis<2;++axis) {
                const double lo=axis?b[0].y:b[0].x,hi=axis?b[1].y:b[1].x;
                const double first=std::max(lo,axis?box.lo.y:box.lo.x),last=std::min(hi,axis?box.hi.y:box.hi.x);
                if(first>last)return true;
                const auto width=Exact::enclosure(Exact::subtract(Exact::fromDouble(hi),Exact::fromDouble(lo)));
                const auto a=Exact::enclosure(Exact::subtract(Exact::fromDouble(first),Exact::fromDouble(lo)));
                const auto z=Exact::enclosure(Exact::subtract(Exact::fromDouble(last),Exact::fromDouble(lo)));
                require(width.lo>0,"unresolved macro graph box width");
                (axis?lower.y:lower.x)=std::max(0.,std::nextafter(a.lo/width.hi,-std::numeric_limits<double>::infinity()));
                (axis?upper.y:upper.x)=std::min(1.,std::nextafter(z.hi/width.lo,std::numeric_limits<double>::infinity()));
            }
            const auto range=scalar.nativeC2MacroScalarBounds(graph.cell,lower,upper);
            return range[0]>0||range[1]<0;
        }
        std::array<LAMP::Geometry::PolynomialAffineCoordinate,2> coordinates;
        for(int axis=0;axis<2;++axis){const double lo=axis?box.lo.y:box.lo.x,hi=axis?box.hi.y:box.hi.x,center=lo+(hi-lo)/2;
            const double radius=std::nextafter(std::max(center-lo,hi-center),std::numeric_limits<double>::infinity());
            coordinates[axis]={center,axis?0:radius,axis?radius:0,{0,0}};}
        const auto range=color.phaseGeometry().publication().state().geometry().snapshot().nativeQuarticAffineEnclosure(graph.cell,coordinates).enclose({-1,1},{-1,1});
        return range.lower>0||range.upper<0;
    }
    bool wholeGraphPacket(int id) const
    {
        const auto& graph=graphs[id];const auto& original=color.phaseGeometry().publication().topology().curves[graph.component].graphs()[graph.graph];
        if(original.nativeSide)return false;const int chosen=graph.coordinateIsX?0:1;
        for(const auto endpoint:{original.begin,original.end})if((endpoint.side==0||endpoint.side==2?0:1)!=chosen)return true;
        return false;
    }
    Range endpointBound(const SurfaceGraph& graph,bool last) const
    {
        const auto& chain=color.phaseGeometry().publication().topology().curves[graph.component];
        const auto& original=chain.graphs()[graph.graph];const int chosen=graph.coordinateIsX?0:1;
        const auto coordinate=[&](Vector2D p){return chosen?p.y:p.x;};
        if(original.nativeSide){const double point=coordinate(last?original.nativeSide->end():original.nativeSide->begin());return {point,point};}
        const auto endpoint=last?original.end:original.begin;const auto bounds=chain.source()->nativeCellBounds(graph.cell);
        require((endpoint.side==0||endpoint.side==2?0:1)==chosen,"uncertified varying graph endpoint bound");
        return embedded(coordinate(bounds[0]),coordinate(bounds[1]),coordinate(endpoint.side<2?original.localLower:original.localUpper));
    }
    bool provedZero(int id,double begin,double end,HydroCsfFaceSupportReport& report,const HydroCsfCurvatureConfig& cfg)
    {
        const auto& graph=graphs[id];Box box=graph.bounds;
        if(!wholeGraphPacket(id)){auto found=endpointBounds.find(id);
        if(found==endpointBounds.end()){const auto first=endpointBound(graph,false),last=endpointBound(graph,true);
            found=endpointBounds.emplace(id,std::array<Range,2>{{first,last}}).first;}
        const auto interpolate=[&](double t,bool upper){const auto& e=found->second;
            const auto v=Exact::enclosure(Exact::add(Exact::multiply(Exact::subtract(Exact::fromDouble(1),Exact::fromDouble(t)),
                Exact::fromDouble(upper?e[0].hi:e[0].lo)),Exact::multiply(Exact::fromDouble(t),Exact::fromDouble(upper?e[1].hi:e[1].lo))));
            return upper?v.hi:v.lo;};
        const double lo=std::min(interpolate(begin,false),interpolate(end,false)),hi=std::max(interpolate(begin,true),interpolate(end,true));
        if(graph.coordinateIsX){box.lo.x=std::max(box.lo.x,lo);box.hi.x=std::min(box.hi.x,hi);}
        else{box.lo.y=std::max(box.lo.y,lo);box.hi.y=std::min(box.hi.y,hi);}
        require(box.lo.x<=box.hi.x&&box.lo.y<=box.hi.y,"invalid certified graph subarc box");}
        struct Pending{Box box;int depth;};std::vector<Pending> pending{{box,0}};
        while(!pending.empty()){const auto item=pending.back();pending.pop_back();const auto& current=item.box;
            if(outsideDisk(current,report.center,report.radius)||shadowed(current,report.center,report.pressureComponent,report,cfg)||noZeroSet(graph,current,report,cfg))continue;
            if(item.depth>=52)return false;
            const bool splitX=current.hi.x-current.lo.x>current.hi.y-current.lo.y;
            const double mid=splitX?current.lo.x+(current.hi.x-current.lo.x)/2:current.lo.y+(current.hi.y-current.lo.y)/2;
            if(splitX?!(mid>current.lo.x&&mid<current.hi.x):!(mid>current.lo.y&&mid<current.hi.y))return false;
            Box left=current,right=current;if(splitX)left.hi.x=right.lo.x=mid;else left.hi.y=right.lo.y=mid;
            pending.push_back({left,item.depth+1});pending.push_back({right,item.depth+1});}
        return true;
    }
};

struct QuadratureRule
{
    std::vector<double> weights,curvatures;
    long double weight=0,moment=0;
    double maxDistance=0;
    bool crossesAmr=false;
    std::size_t radialInside=0,radialOutside=0;
};
QuadratureRule integrateRule(ProjectionIndex& index,int graph,double begin,double end,HydroCsfFaceSupportReport& report,
    const HydroCsfCurvatureConfig& cfg,int faceLevel)
{
    QuadratureRule rule;
    for(const auto& q:index.samples(graph,begin,end,report,cfg)){
        const double distance=std::hypot(q.point.x-report.center.x,q.point.y-report.center.y);
        if(distance>=report.radius){if(distance>report.radius)++rule.radialOutside;continue;}++rule.radialInside;
        if(!index.visible(report.center,q.point,report.pressureComponent,report,cfg.maxVisibilityCellQueriesPerFace)){++report.visibilityRejected;continue;}
        const double ratio=distance/report.radius,base=(1-ratio)*(1+ratio),kernel=base*base*base,weight=q.weight*kernel;
        require(std::isfinite(weight)&&weight>0,"unrepresentable positive CSF weight");
        if(!report.quadratureAnchorDefined){report.quadratureAnchor=q.curvature;report.quadratureAnchorDefined=true;}
        index.graphPositiveStamp[graph]=index.graphVisitEpoch;
        rule.weights.push_back(weight);rule.curvatures.push_back(q.curvature);rule.weight+=weight;
        const long double delta=static_cast<long double>(q.curvature)-report.quadratureAnchor,term=static_cast<long double>(weight)*delta;
        require(std::isfinite(delta)&&std::isfinite(term)&&(term!=0||delta==0),"unrepresentable CSF quadrature anchored term");
        rule.moment+=term;
        require(std::isfinite(rule.weight)&&std::isfinite(rule.moment),"unrepresentable CSF quadrature moment");
        report.curvatureScale=std::max(report.curvatureScale,std::abs(q.curvature));
        rule.maxDistance=std::max(rule.maxDistance,distance);rule.crossesAmr=rule.crossesAmr||q.level!=faceLevel;}
    return rule;
}
QuadratureRule combine(QuadratureRule a,QuadratureRule b)
{
    a.weight+=b.weight;a.moment+=b.moment;a.maxDistance=std::max(a.maxDistance,b.maxDistance);a.crossesAmr=a.crossesAmr||b.crossesAmr;
    a.radialInside+=b.radialInside;a.radialOutside+=b.radialOutside;
    require(std::isfinite(a.weight)&&std::isfinite(a.moment),"unrepresentable combined CSF quadrature moment");
    a.weights.insert(a.weights.end(),b.weights.begin(),b.weights.end());a.curvatures.insert(a.curvatures.end(),b.curvatures.begin(),b.curvatures.end());return a;
}
struct AdaptivePanel
{
    int graph;double begin,end;int depth;std::size_t subdivisions;
    QuadratureRule fine;
    long double weightError,momentError;
    bool unresolved,radialUnresolved;
};
// A detected radial cutoff uses a composite rule: a single coarse subarc can
// underresolve its tail even when one fine node is positive.
AdaptivePanel evaluatePanel(ProjectionIndex& index,int graph,double begin,double end,int depth,
    HydroCsfFaceSupportReport& report,const HydroCsfCurvatureConfig& cfg,int faceLevel,std::size_t subdivisions=1)
{
    if(!(subdivisions>0&&subdivisions<=cfg.maxSurfaceQueriesPerFace/24&&report.surfaceQueries<=cfg.maxSurfaceQueriesPerFace-24*subdivisions)){
        std::ostringstream message;message<<"[HydroBalancedCsfLoad] CSF surface query limit, graph="<<graph<<", subdivisions="<<subdivisions;
        throw std::invalid_argument(message.str());}
    QuadratureRule coarse,fine;
    // A whole-graph packet is refined uniformly. Its zero tails have no
    // independent absence claim; the error is the complete N/2N rule change.
    for(std::size_t i=0;i<subdivisions;++i){const double a=begin+(end-begin)*(double(i)/subdivisions),b=begin+(end-begin)*(double(i+1)/subdivisions);
        coarse=combine(std::move(coarse),integrateRule(index,graph,a,b,report,cfg,faceLevel));}
    for(std::size_t i=0;i<2*subdivisions;++i){const double a=begin+(end-begin)*(double(i)/(2*subdivisions)),b=begin+(end-begin)*(double(i+1)/(2*subdivisions));
        require(b>a,"unrepresentable CSF quadrature subdivision");
        fine=combine(std::move(fine),integrateRule(index,graph,a,b,report,cfg,faceLevel));}
    bool endInside=false,endOutside=false;
    if(index.convexRectangle){
        for(double fraction:{begin,end}){const auto point=index.endpoint(graph,fraction,report,cfg);
            const double distance=std::hypot(point.x-report.center.x,point.y-report.center.y);
            endInside=endInside||distance<report.radius;endOutside=endOutside||distance>report.radius;}}
    if(subdivisions==1&&((coarse.radialInside+fine.radialInside>0&&coarse.radialOutside+fine.radialOutside>0)||(endInside&&endOutside)))
        return evaluatePanel(index,graph,begin,end,depth,report,cfg,faceLevel,2);
    const bool positiveEvidence=report.quadratureAnchorDefined;
    const bool unresolved=coarse.weight==0&&fine.weight==0&&!positiveEvidence&&!index.provedZero(graph,begin,end,report,cfg);
    bool radialUnresolved=false;
    if(index.convexRectangle){
        radialUnresolved=(endInside&&coarse.radialInside+fine.radialInside==0)||(endOutside&&coarse.radialOutside+fine.radialOutside==0);
        if(begin==0&&end==1&&coarse.weight==0&&fine.weight==0&&!index.provedZero(graph,begin,end,report,cfg))radialUnresolved=true;
    }
    const long double weightError=std::abs(fine.weight-coarse.weight),momentError=std::abs(fine.moment-coarse.moment);
    report.deepestPanel=std::max(report.deepestPanel,depth);
    return {graph,begin,end,depth,subdivisions,std::move(fine),weightError,momentError,unresolved,radialUnresolved};
}
void projectAdaptive(ProjectionIndex& index,HydroCsfFaceSupportReport& report,const HydroCsfCurvatureConfig& cfg,int faceLevel,
    std::vector<double>& weights,std::vector<double>& curvatures,std::vector<HydroCsfQuadraturePanel>& ledger)
{
    std::vector<AdaptivePanel> panels;std::vector<int> candidates;
    if(++index.graphVisitEpoch==0){std::fill(index.graphVisitStamp.begin(),index.graphVisitStamp.end(),0);
        std::fill(index.graphPositiveStamp.begin(),index.graphPositiveStamp.end(),0);index.graphVisitEpoch=1;}
    index.graphBins.visit(surrounding(report.center,report.radius),[&](int id){if(index.graphVisitStamp[id]==index.graphVisitEpoch)return;index.graphVisitStamp[id]=index.graphVisitEpoch;
        require(report.candidateQueries<cfg.maxCandidateQueriesPerFace,"CSF support query limit");++report.candidateQueries;
        const auto& graph=index.graphs[id];if(index.component[graph.pressure]!=report.pressureComponent||ProjectionIndex::outsideDisk(graph.bounds,report.center,report.radius))return;
        candidates.push_back(id);});
    for(int id:candidates){panels.push_back(evaluatePanel(index,id,0,1,0,report,cfg,faceLevel));}
    for(;;){long double total=0,moment=0,weightError=0,momentError=0;int worst=-1;long double worstScore=-1;bool unresolved=false;
        const double scale=std::max(report.curvatureScale,1/report.radius);
        require(std::isfinite(scale)&&scale>0,"unrepresentable CSF quadrature scale");
        for(const auto& panel:panels){total+=panel.fine.weight;moment+=panel.fine.moment;weightError+=panel.weightError;momentError+=panel.momentError;
        }
        const bool conditioned=total>0&&weightError<=.01L*total;
        const long double offset=total>0?std::abs(moment/total):0;
        for(std::size_t i=0;i<panels.size();++i){auto& panel=panels[i];
            const bool lost=total==0&&report.quadratureAnchorDefined;
            const bool unresolvedNow=panel.radialUnresolved||lost||(panel.unresolved&&!report.quadratureAnchorDefined);
            const long double score=unresolvedNow?std::numeric_limits<long double>::infinity():
                conditioned?(panel.momentError+offset*panel.weightError)/scale:panel.weightError;
            unresolved=unresolved||unresolvedNow;
            if(score>worstScore||(score==worstScore&&unresolvedNow&&worst>=0&&panel.depth<panels[worst].depth)){
                worstScore=score;worst=static_cast<int>(i);}}
        const long double curvatureError=total>0?(total>weightError?
            (momentError+offset*weightError)/(total-weightError):std::numeric_limits<double>::max()):0;
        require(std::isfinite(total)&&std::isfinite(moment)&&std::isfinite(weightError)&&std::isfinite(momentError)&&std::isfinite(curvatureError),
            "unrepresentable CSF quadrature error estimate");
        report.totalWeight=total;report.weightErrorEstimate=weightError;report.anchoredMomentErrorEstimate=momentError;
        report.curvatureErrorEstimate=static_cast<double>(curvatureError);
        if(!unresolved&&(total==0?(weightError==0&&momentError==0&&!report.quadratureAnchorDefined):
                (conditioned&&curvatureError<=cfg.relativeQuadratureTolerance*scale))){
            report.weightErrorEstimate=weightError;report.anchoredMomentErrorEstimate=momentError;report.curvatureErrorEstimate=static_cast<double>(curvatureError);break;}
        require(worst>=0,"CSF quadrature has no refinable panel");auto panel=std::move(panels[worst]);
        if(panel.depth>=cfg.maxQuadratureDepth){std::ostringstream message;message<<std::setprecision(17)<<
            (panel.unresolved?"[HydroBalancedCsfLoad] unsupported active CSF face without local selected/visible curvature support: unresolved quadrature depth":
                "[HydroBalancedCsfLoad] CSF quadrature depth limit")<<", center=("<<report.center.x<<','<<report.center.y<<"), graph="<<panel.graph<<
                ", fraction=("<<panel.begin<<','<<panel.end<<"), queries="<<report.surfaceQueries;
            throw std::invalid_argument(message.str());}
        if(!report.quadratureAnchorDefined||total==0||(panel.radialUnresolved&&panel.begin==0&&panel.end==1&&panel.fine.weight==0)){
            require(panel.subdivisions<=cfg.maxSurfaceQueriesPerFace/48,"CSF surface query limit");
            auto searched=evaluatePanel(index,panel.graph,0,1,panel.depth+1,report,cfg,faceLevel,2*panel.subdivisions);
            panels.erase(std::remove_if(panels.begin(),panels.end(),[&](const auto& p){return p.graph==panel.graph;}),panels.end());
            if(searched.fine.weight>0){const auto count=searched.subdivisions;
                // Minimum batch cost; any cutoff promotion is additionally
                // charged and checked before each actual surface query.
                const std::size_t minimumPanelQueries=index.convexRectangle?26:24;
                require(count<=cfg.maxSurfaceQueriesPerFace/minimumPanelQueries&&
                    report.surfaceQueries<=cfg.maxSurfaceQueriesPerFace-minimumPanelQueries*count,"CSF surface query limit");
                panels.push_back(evaluatePanel(index,panel.graph,0,1./count,searched.depth,report,cfg,faceLevel));
                for(std::size_t i=1;i<count;++i)panels.push_back(evaluatePanel(index,panel.graph,double(i)/count,double(i+1)/count,searched.depth,report,cfg,faceLevel));
            }else panels.push_back(std::move(searched));
        }else{const double mid=panel.begin+(panel.end-panel.begin)/2;
            panels[worst]=evaluatePanel(index,panel.graph,panel.begin,mid,panel.depth+1,report,cfg,faceLevel);
            panels.push_back(evaluatePanel(index,panel.graph,mid,panel.end,panel.depth+1,report,cfg,faceLevel));}}
    std::map<int,long double> finalGraphWeights;
    for(auto& panel:panels){const auto& graph=index.graphs[panel.graph];finalGraphWeights[panel.graph]+=panel.fine.weight;
        ledger.push_back({graph.component,graph.graph,panel.begin,panel.end,panel.depth,panel.subdivisions,
            report.quadratureAnchorDefined,panel.weightError,panel.momentError});
        weights.insert(weights.end(),panel.fine.weights.begin(),panel.fine.weights.end());curvatures.insert(curvatures.end(),panel.fine.curvatures.begin(),panel.fine.curvatures.end());
        report.maxSupportDistance=std::max(report.maxSupportDistance,panel.fine.maxDistance);report.crossesAmr=report.crossesAmr||panel.fine.crossesAmr;}
    for(const auto& entry:finalGraphWeights)if(entry.second==0&&index.graphPositiveStamp[entry.first]==index.graphVisitEpoch)
        ++report.trialPositiveGraphsWithoutFineSupport;
    report.acceptedPanels=ledger.size();report.supportCount=weights.size();report.totalWeight=0;
    for(std::size_t i=0;i<weights.size();++i){report.totalWeight+=weights[i];if(i==0)report.curvatureMin=report.curvatureMax=curvatures[i];
        else{report.curvatureMin=std::min(report.curvatureMin,curvatures[i]);report.curvatureMax=std::max(report.curvatureMax,curvatures[i]);}}
}

// Geometric alpha has no diffuse force support. Each original cell graph is
// integrated once for its own cell; faces only consume these local means.
struct GeometricRule
{
    std::vector<double> weights,curvatures;
    std::vector<Vector2D> points;
    long double weight=0,moment=0;
};
struct GeometricPanel
{
    std::size_t component,graph;
    double begin,end;int depth;
    GeometricRule fine;
    long double weightError,momentError;
};
struct GeometricCellCurvature
{
    bool defined=false,anchored=false;
    double anchor=0,scale=0,average=0,minimum=0,maximum=0;
    long double area=0,weightError=0,momentError=0,curvatureError=0;
    std::size_t queries=0,supportCount=0;
    int deepestPanel=0;
    std::vector<GeometricPanel> panels;
};
GeometricRule geometricRule(const HydroPhaseGeometry& phase,int pressure,std::size_t c,std::size_t g,
    double begin,double end,GeometricCellCurvature& cell,const HydroCsfCurvatureConfig& cfg)
{
    const auto samples=phase.publication().quadrature(c,g,begin,end);
    require(cell.queries<=cfg.maxSurfaceQueriesPerFace&&samples.size()<=cfg.maxSurfaceQueriesPerFace-cell.queries,
        "geometric cell curvature surface query limit");
    cell.queries+=samples.size();GeometricRule rule;
    const auto clock=phase.timeLayer();
    for(const auto& q:samples){
        const auto time=q.sample.timeLayer;
        require(q.sample.geometryIdentity==phase.geometryIdentity()&&time.step==clock.step&&time.begin==clock.begin
            &&time.end==clock.end&&time.position==clock.position
            &&q.sample.frame.nativeCellIndex==phase.layout().nativeCell(pressure),
            "geometric cell curvature lost local publication/full clock");
        require(std::isfinite(q.physicalWeight)&&q.physicalWeight>0&&std::isfinite(q.sample.curvature),
            "invalid geometric cell curvature positive physical sample");
        if(!cell.anchored){cell.anchor=q.sample.curvature;cell.anchored=true;}
        const long double delta=static_cast<long double>(q.sample.curvature)-cell.anchor;
        const long double term=static_cast<long double>(q.physicalWeight)*delta;
        require(std::isfinite(delta)&&std::isfinite(term)&&(term!=0||delta==0),
            "unrepresentable geometric cell curvature anchored moment");
        rule.weight+=q.physicalWeight;rule.moment+=term;
        rule.weights.push_back(q.physicalWeight);rule.curvatures.push_back(q.sample.curvature);rule.points.push_back(q.sample.frame.point);
        cell.scale=std::max(cell.scale,std::abs(q.sample.curvature));
    }
    require(std::isfinite(rule.weight)&&std::isfinite(rule.moment),"nonfinite geometric cell curvature integral");return rule;
}
GeometricPanel geometricPanel(const HydroPhaseGeometry& phase,int p,std::size_t c,std::size_t g,
    double begin,double end,int depth,GeometricCellCurvature& cell,const HydroCsfCurvatureConfig& cfg)
{
    const double mid=begin+(end-begin)/2;
    require(mid>begin&&mid<end,"unrepresentable geometric cell curvature subdivision");
    const auto coarse=geometricRule(phase,p,c,g,begin,end,cell,cfg);
    auto fine=geometricRule(phase,p,c,g,begin,mid,cell,cfg);
    auto right=geometricRule(phase,p,c,g,mid,end,cell,cfg);
    fine.weight+=right.weight;fine.moment+=right.moment;
    fine.weights.insert(fine.weights.end(),right.weights.begin(),right.weights.end());
    fine.curvatures.insert(fine.curvatures.end(),right.curvatures.begin(),right.curvatures.end());
    fine.points.insert(fine.points.end(),right.points.begin(),right.points.end());
    cell.deepestPanel=std::max(cell.deepestPanel,depth);
    const auto weightError=std::abs(fine.weight-coarse.weight),momentError=std::abs(fine.moment-coarse.moment);
    return {c,g,begin,end,depth,std::move(fine),weightError,momentError};
}
void geometricCell(const HydroPhaseGeometry& phase,int p,const std::vector<std::pair<std::size_t,std::size_t>>& graphs,
    GeometricCellCurvature& cell,const HydroCsfCurvatureConfig& cfg)
{
    for(const auto& graph:graphs)cell.panels.push_back(geometricPanel(phase,p,graph.first,graph.second,0,1,0,cell,cfg));
    for(;;){
        long double area=0,moment=0,weightError=0,momentError=0;
        for(const auto& panel:cell.panels){area+=panel.fine.weight;moment+=panel.fine.moment;
            weightError+=panel.weightError;momentError+=panel.momentError;}
        require(area>0&&std::isfinite(area)&&std::isfinite(moment),"geometric cell curvature has no positive interface area");
        const long double offset=std::abs(moment/area);
        const long double curvatureError=area>weightError?(momentError+offset*weightError)/(area-weightError):
            std::numeric_limits<double>::infinity();
        const auto bounds=phase.publication().state().geometry().snapshot().nativeCellBounds(phase.layout().nativeCell(p));
        const double h=std::min(bounds[1].x-bounds[0].x,bounds[1].y-bounds[0].y);
        const double scale=std::max(cell.scale,1/h);cell.scale=scale;
        require(std::isfinite(scale)&&scale>0,"unrepresentable geometric cell curvature error scale");
        if(weightError<=cfg.relativeQuadratureTolerance*area&&curvatureError<=cfg.relativeQuadratureTolerance*scale){
            cell.area=area;cell.weightError=weightError;cell.momentError=momentError;cell.curvatureError=curvatureError;break;
        }
        auto worst=std::max_element(cell.panels.begin(),cell.panels.end(),[&](const auto& a,const auto& b){
            return a.momentError+offset*a.weightError+scale*a.weightError
                <b.momentError+offset*b.weightError+scale*b.weightError;});
        require(worst!=cell.panels.end(),"geometric cell curvature has no refinable graph");
        const auto old=*worst;
        if(old.depth>=cfg.maxQuadratureDepth){std::ostringstream message;message<<std::setprecision(17)
            <<"[HydroBalancedCsfLoad] geometric cell curvature quadrature depth limit, pressure="<<p
            <<", native_cell="<<phase.layout().nativeCell(p)<<", component="<<old.component<<", graph="<<old.graph
            <<", fraction=("<<old.begin<<','<<old.end<<"), area_error="<<weightError<<", curvature_error="<<curvatureError;
            throw std::invalid_argument(message.str());}
        const double mid=old.begin+(old.end-old.begin)/2;
        *worst=geometricPanel(phase,p,old.component,old.graph,old.begin,mid,old.depth+1,cell,cfg);
        cell.panels.push_back(geometricPanel(phase,p,old.component,old.graph,mid,old.end,old.depth+1,cell,cfg));
    }
    std::vector<double> weights,curvatures;
    for(const auto& panel:cell.panels){weights.insert(weights.end(),panel.fine.weights.begin(),panel.fine.weights.end());
        curvatures.insert(curvatures.end(),panel.fine.curvatures.begin(),panel.fine.curvatures.end());}
    cell.average=HydroBalancedCsfLoad::positiveCurvatureAverage(weights,curvatures);
    const auto range=std::minmax_element(curvatures.begin(),curvatures.end());
    cell.minimum=*range.first;cell.maximum=*range.second;cell.supportCount=weights.size();cell.defined=true;
}


// Local tangent of the actual FirstJet extraction/compact closure/scalar map.
// The cached primary rates use the direction field solely to freeze each
// one-sided scalar Hessian. No native unit-field probes or carried Hessians.
struct ProjectedRateBasis
{
    using Sparse=std::map<int,long double>;
    using Key=std::tuple<int,int,int>; // Q4 sensitivity kind, index, component
    const QuadtreeMesh& mesh;
    const HydroPhaseGeometry& phase;
    const HydroSpatialAdvector& motion;
    const Interface::HermiteGalsAmrClosureSensitivity& mixed;
    const Detail::HydroRegularAdvectorGeometry& regular;
    std::vector<std::optional<Detail::HydroRegularCellStencil>> stencils;
    std::map<int,std::size_t> vertexOrdinal;
    std::map<Key,Sparse> rates;
    std::map<std::pair<double,double>,std::array<Sparse,3>> pointRates;
    bool compatibleIncrement=false;
    std::optional<UniformLeafGrid> incrementGrid;
    std::vector<Vector2D> incrementPoints;
    std::map<std::pair<double,double>,std::array<int,2>> incrementPositions;
    std::map<std::pair<double,double>,Sparse> phiValueRates;
    ProjectedRateBasis(const QuadtreeMesh& m,const HydroPhaseGeometry& p,const HydroSpatialAdvector& u,
        const Interface::HermiteGalsAmrClosureSensitivity& c,bool increment=false)
        :mesh(m),phase(p),motion(u),mixed(c),regular(u.reconstruction().regularGeometry()),
        stencils(p.layout().pressureCount()),compatibleIncrement(increment)
    {
        const auto& primary=p.publication().state().primary();
        for(std::size_t k=0;k<primary.vertices.size();++k)
            require(vertexOrdinal.emplace(primary.vertices[k].vertex,k).second,"duplicate projected primary vertex");
        if(compatibleIncrement){
            incrementGrid.emplace(mesh,CartesianGridScope::UniformLeaves);
            const auto& grid=*incrementGrid;const int nx=grid.nx(),ny=grid.ny();
            require(nx>=2&&ny>=2,"increment rate requires five physical half-grid nodes per axis");
            incrementPoints.resize((2*nx+1)*(2*ny+1));
            const auto vertex=[&](int i,int j){return mesh.verts(grid.vertices()[j*(nx+1)+i]).coord();};
            for(int j=0;j<=2*ny;++j)for(int i=0;i<=2*nx;++i){Vector2D q;
                if(i%2&&j%2)q=mesh.cells(grid.cells()[(j/2)*nx+i/2]).centerCoord();
                else if(i%2){const auto a=vertex(i/2,j/2),b=vertex(i/2+1,j/2);q=a+(b-a)*.5;}
                else if(j%2){const auto a=vertex(i/2,j/2),b=vertex(i/2,j/2+1);q=a+(b-a)*.5;}
                else q=vertex(i/2,j/2);
                incrementPoints[j*(2*nx+1)+i]=q;
                require(incrementPositions.emplace(std::make_pair(q.x,q.y),std::array<int,2>{i,j}).second,
                    "duplicate physical increment-rate half-grid point");
            }
            // Extract precisely the same primary/mixed/D4 dependencies of
            // interface leaves, then authenticate their source stencils in
            // ONE batch. No repeated whole-mesh identity walk per leaf.
            using PointKey=std::pair<double,double>;
            std::set<PointKey> jetPoints,valuePoints;std::set<int> nativeCells;
            const auto addJet=[&](Vector2D q){jetPoints.emplace(q.x,q.y);valuePoints.emplace(q.x,q.y);};
            const auto edgePoint=[&](const auto& e){const auto a=vertexPoint(e.vertices[0]),b=vertexPoint(e.vertices[1]);return a+(b-a)*.5;};
            for(const auto& curve:phase.publication().topology().curves)
                for(const auto& graph:curve.graphs())nativeCells.insert(graph.nativeCellIndex);
            for(int native:nativeCells){const auto& cell=mesh.cells(mesh.cellIndexToLabel[native]);
                std::array<Vector2D,4> corners{};int corner=0;
                for(int label:cell.cellToVert){const int v=mesh.verts(label).index();const auto q=vertexPoint(v);
                    corners[corner++]=q;addJet(q);const auto ordinal=vertexOrdinal.at(v);
                    for(const auto& term:mixed.rows().at(ordinal)){
                        if(term.component==Interface::HermiteGalsPrimaryDerivative::EdgeNormal)addJet(edgePoint(primary.edges.at(term.primaryOrdinal)));
                        else addJet(vertexPoint(primary.vertices.at(term.primaryOrdinal).vertex));
                    }
                }
                for(int a=0;a<4;++a)for(int b=a+1;b<4;++b)
                    if(corners[a].x==corners[b].x||corners[a].y==corners[b].y)addJet(corners[a]+(corners[b]-corners[a])*.5);
                const auto center=cell.centerCoord();valuePoints.emplace(center.x,center.y);
            }
            for(const auto& key:jetPoints){const auto ij=incrementPositions.at(key);const Vector2D q{key.first,key.second};
                for(int axis=0;axis<2;++axis){
                    // Edge primary carries ONLY its normal derivative;
                    // tangential and center gradient rates are not data.
                    if((ij[0]%2&&axis==0)||(ij[1]%2&&axis==1))continue;
                    bool constrained=false;
                    for(const auto& plane:mixed.symmetryPlanes())
                        if(plane.normalAxis==axis&&(axis==0?q.x:q.y)==plane.coordinate)constrained=true;
                    if(constrained)continue;
                    const int count=2*(axis==0?nx:ny)+1,begin=std::clamp(ij[axis]-2,0,count-5);
                    for(int k=begin;k<begin+5;++k){const auto p=axis==0?incrementPoints[ij[1]*(2*nx+1)+k]
                        :incrementPoints[k*(2*nx+1)+ij[0]];valuePoints.emplace(p.x,p.y);}
                }
            }
            std::set<int> required;
            const auto& scalar=phase.publication().state().geometry().snapshot();
            for(const auto& key:valuePoints){int native=-1;const Vector2D q{key.first,key.second};
                for(int c:scalar.nativeCellCandidates(q,q))if(scalar.nativeCellActive(c)&&phase.layout().pressureDof(c)>=0
                    &&(native<0||c<native))native=c;
                require(native>=0,"increment-rate halo has no active scalar/source cell");required.insert(phase.layout().pressureDof(native));}
            const std::vector<int> pressures(required.begin(),required.end());auto prepared=regular.cellStencils(mesh,pressures,motion.boundary());
            for(std::size_t k=0;k<pressures.size();++k)stencils[pressures[k]]=std::move(prepared[k]);
        }
    }
    Vector2D vertexPoint(int v) const {return mesh.verts(mesh.vertIndexToLabel[v]).coord();}
    // Only numerical rates on the requested interface/closure halo are
    // sampled. The retained full half-grid is immutable physical topology.
    const Sparse& phiRate(Vector2D q)
    {
        const auto key=std::make_pair(q.x,q.y);const auto found=phiValueRates.find(key);
        if(found!=phiValueRates.end())return found->second;
        const auto& scalar=phase.publication().state().geometry().snapshot();
        int native=-1;
        for(int c:scalar.nativeCellCandidates(q,q))if(scalar.nativeCellActive(c)&&phase.layout().pressureDof(c)>=0
            &&(native<0||c<native))native=c;
        require(native>=0,"increment rate point has no active scalar/source support");
        // phi/g are the shared C1 source. In this method no donor Hessian,
        // transported-Jacobian gradient rate or independently fitted field enters.
        const auto bounds=scalar.nativeCellBounds(native);
        // Reproduce the source's exact same endpoint-relative physical chart.
        // scalarJetNative would additionally compute an exact homogeneous
        // value that this gradient-only rate never consumes.
        const auto local=[](double point,double low,double high){
            require(std::isfinite(point)&&point>=low&&point<=high,"increment phi rate leaves actual closed leaf");
            const bool upper=point-low>high-point;
            const double offset=static_cast<double>(upper?(static_cast<long double>(high)-point)/(static_cast<long double>(high)-low)
                :(static_cast<long double>(point)-low)/(static_cast<long double>(high)-low));
            return LAMP::Geometry::CellChartCoordinate{offset,upper};
        };
        const LAMP::Geometry::CellChartPoint chart{local(q.x,bounds[0].x,bounds[1].x),local(q.y,bounds[0].y,bounds[1].y)};
        const auto j=scalar.scalarJetAtNativeCellChart(native,chart);const int p=phase.layout().pressureDof(native);
        require(stencils[p].has_value(),"increment phi rate left authenticated stencil halo");
        Sparse value;
        for(const auto& row:regular.cellJetWeights(*stencils[p],q)){
            const long double rate=-static_cast<long double>(row.jet[0])*j.dx-static_cast<long double>(row.jet[1])*j.dy;
            require(std::isfinite(rate),"nonfinite increment phi rate");
            if(rate!=0)value[row.velocityDof]+=rate;
        }
        return phiValueRates.emplace(key,std::move(value)).first->second;
    }
    std::array<Sparse,3> incrementPoint(Vector2D q)
    {
        const auto ij=incrementPositions.at({q.x,q.y});const auto& grid=*incrementGrid;
        std::array<Sparse,3> value;value[0]=phiRate(q);
        for(int axis=0;axis<2;++axis){
            if((ij[0]%2&&axis==0)||(ij[1]%2&&axis==1))continue;
            bool constrained=false;
            for(const auto& plane:mixed.symmetryPlanes())
                if(plane.normalAxis==axis&&(axis==0?q.x:q.y)==plane.coordinate)constrained=true;
            if(constrained)continue;
            const int count=2*(axis==0?grid.nx():grid.ny())+1,target=ij[axis],begin=std::clamp(target-2,0,count-5);
            const auto pointAt=[&](int k){return axis==0?incrementPoints[ij[1]*(2*grid.nx()+1)+k]
                :incrementPoints[k*(2*grid.nx()+1)+ij[0]];};
            std::array<const Sparse*,5> samples{};std::set<int> dofs;
            for(int k=0;k<5;++k){samples[k]=&phiRate(pointAt(begin+k));
                for(const auto& entry:*samples[k])dofs.insert(entry.first);}
            for(int u:dofs){
                const double rate=LAMP::Geometry::compactQuarticDerivativeAtNode(target,count,
                    [&](int k){const auto p=pointAt(k);return axis==0?p.x:p.y;},
                    [&](int k){const auto& sample=*samples[k-begin];const auto found=sample.find(u);
                        return found==sample.end()?0.L:found->second;},require);
                if(rate!=0)value[axis+1].emplace(u,rate);
            }
        }
        return value;
    }
    const std::array<Sparse,3>& point(Vector2D q)
    {
        const auto key=std::make_pair(q.x,q.y);const auto found=pointRates.find(key);
        if(found!=pointRates.end())return found->second;
        if(compatibleIncrement)return pointRates.emplace(key,incrementPoint(q)).first->second;
        const auto& scalar=phase.publication().state().geometry().snapshot();
        const auto candidates=scalar.nativeCellCandidates(q,q);
        int sample=-1;
        for(int c:candidates)if(scalar.nativeCellActive(c)&&phase.layout().pressureDof(c)>=0
            &&(sample<0||c<sample))sample=c;
        require(sample>=0,"projected primary point has no active scalar/source support");
        const auto actual=motion.cellVelocity(phase.layout().pressureDof(sample),q);
        const Vector2D backward{-actual.x,-actual.y};int donor=-1;
        for(int c:candidates){
            if(!scalar.nativeCellActive(c)||phase.layout().pressureDof(c)<0)continue;
            const auto b=scalar.nativeCellBounds(c);
            if(q.x<b[0].x||q.x>b[1].x||q.y<b[0].y||q.y>b[1].y)continue;
            // Exactly the characteristic inward rule: no offset, tolerance,
            // clamping or Hessian averaging. Zero-normal ties retain the
            // same minimum native-index convention as scalarJetNative.
            if((q.x==b[0].x&&backward.x<0)||(q.x==b[1].x&&backward.x>0)
                ||(q.y==b[0].y&&backward.y<0)||(q.y==b[1].y&&backward.y>0))continue;
            if(donor<0||c<donor)donor=c;
        }
        require(donor>=0,"projected backward donor exits active closed domain");
        const auto j=scalar.scalarJetInNativeCell(donor,q);
        const int p=phase.layout().pressureDof(donor);
        if(!stencils[p])stencils[p]=regular.cellStencil(mesh,p,motion.boundary());
        std::array<Sparse,3> value;
        for(const auto& row:regular.cellJetWeights(*stencils[p],q)){
            const auto& b=row.jet;
            value[0][row.velocityDof]=-static_cast<long double>(b[0])*j.dx-static_cast<long double>(b[1])*j.dy;
            value[1][row.velocityDof]=-static_cast<long double>(b[2])*j.dx-static_cast<long double>(b[4])*j.dy
                -static_cast<long double>(j.dxx)*b[0]-static_cast<long double>(j.dxy)*b[1];
            value[2][row.velocityDof]=-static_cast<long double>(b[3])*j.dx-static_cast<long double>(b[5])*j.dy
                -static_cast<long double>(j.dxy)*b[0]-static_cast<long double>(j.dyy)*b[1];
        }
        // This is the actual scalar Neumann FirstJet constraint, not a
        // velocity parity law: retain value and tangential gradient rates.
        for(const auto& plane:mixed.symmetryPlanes())
            if((plane.normalAxis==0?q.x:q.y)==plane.coordinate)value[plane.normalAxis+1].clear();
        for(const auto& component:value)for(const auto& entry:component)
            require(std::isfinite(entry.second),"nonfinite projected primary rate");
        return pointRates.emplace(key,std::move(value)).first->second;
    }
    const Sparse& primary(int kind,int index,int component)
    {
        const Key key{kind,index,component};const auto found=rates.find(key);
        if(found!=rates.end())return found->second;
        using S=LAMP::Geometry::ActiveQuarticHermiteSensitivity;
        const auto& state=phase.publication().state().primary();Sparse result;
        if(kind==static_cast<int>(S::Kind::Vertex)){
            const auto ordinal=vertexOrdinal.find(index);require(ordinal!=vertexOrdinal.end(),"missing projected master vertex");
            require(component>=0&&component<=3,"invalid projected vertex component");
            if(component<3)result=point(vertexPoint(index))[component];
            else for(const auto& term:mixed.rows().at(ordinal->second)){
                const Sparse* source=nullptr;
                if(term.component==Interface::HermiteGalsPrimaryDerivative::EdgeNormal){
                    require(term.primaryOrdinal<state.edges.size(),"projected edge closure ordinal out of range");
                    source=&primary(static_cast<int>(S::Kind::EdgeMidpoint),static_cast<int>(term.primaryOrdinal),1);
                }else{
                    require(term.primaryOrdinal<state.vertices.size(),"projected vertex closure ordinal out of range");
                    const int c=term.component==Interface::HermiteGalsPrimaryDerivative::VertexDx?1:2;
                    source=&primary(kind,state.vertices[term.primaryOrdinal].vertex,c);
                }
                for(const auto& entry:*source)result[entry.first]+=static_cast<long double>(term.weight)*entry.second;
            }
        }else if(kind==static_cast<int>(S::Kind::EdgeMidpoint)){
            require(index>=0&&static_cast<std::size_t>(index)<state.edges.size()&&(component==0||component==1),
                "invalid projected edge primary");
            const auto& edge=state.edges[index];const auto a=vertexPoint(edge.vertices[0]),b=vertexPoint(edge.vertices[1]);
            result=point(a+(b-a)*.5)[component==0?0:a.y==b.y?2:1];
        }else{
            require(kind==static_cast<int>(S::Kind::CellCenter)&&component==0,"invalid projected center primary");
            const auto b=phase.publication().state().geometry().snapshot().nativeCellBounds(index);
            if(compatibleIncrement)result=phiRate(mesh.cells(mesh.cellIndexToLabel[index]).centerCoord());
            else result=point(b[0]+(b[1]-b[0])*.5)[0];
        }
        return rates.emplace(key,std::move(result)).first->second;
    }
    Sparse normal(const LevelSetHermiteGeometry::Frame& frame)
    {
        const auto& scalar=phase.publication().state().geometry().snapshot();
        const auto chart=LevelSetHermiteGeometry::cellChartPoint(frame);
        const auto j=scalar.scalarJetAtNativeCellChart(frame.nativeCellIndex,chart);
        const long double norm=std::hypot(static_cast<long double>(j.dx),static_cast<long double>(j.dy));
        require(std::isfinite(norm)&&norm>0,"unresolved projected interface gradient");
        Sparse result;
        // Macro sensitivity already pulls all 49 completed data back to the
        // original shared formal sources, including neighboring cell centers.
        // primary() then applies the SAME actual mixed closure/rate generator.
        const auto sensitivity=scalar.usesC2Macro()
            ?scalar.nativeC2MacroDifferentialSensitivity(frame.nativeCellIndex,chart,{1,0,0,0,0,0})
            :scalar.nativeQuarticDifferentialSensitivity(frame.nativeCellIndex,chart,{1,0,0,0,0,0});
        for(const auto& term:sensitivity)
            for(const auto& entry:primary(static_cast<int>(term.kind),term.index,term.component))
                result[entry.first]-=term.coefficient*entry.second/norm;
        for(const auto& entry:result)require(std::isfinite(entry.second),"nonfinite projected normal basis");
        return result;
    }
};

// Same-source weak traction moments. Unlike the cell-average factory, the
// accepted quadrature certifies each actual native load and normal moment,
// including quarter-patch source knots and signed cancellation.
struct AdjointMoment
{
    long double normal=0,traction=0,normalScale=0,tractionScale=0;
    long double tangent=0,tangentScale=0; // positive physical tangent contraction
    long double sourceScale=0; // external only: integral |R_x|+|R_y|, plus |P| when projected
    double minimum=0,maximum=0,maxDistance=0;
    std::size_t nodes=0;
};
struct AdjointRule
{
    long double area=0,moment=0,stressMoment=0;
    std::map<int,AdjointMoment> native;
    std::size_t nodes=0;
};
struct AdjointPanel
{
    std::size_t component=0,graph=0;int pressure=0,depth=0;
    double begin=0,end=1;
    AdjointRule fine;
    // These exact two fine rules become the child coarse rules on refinement.
    std::array<AdjointRule,2> halves;
    long double areaError=0,momentError=0,stressMomentError=0;
    std::map<int,std::array<long double,3>> errors; // normal, anchored full traction, tangent
};
struct AdjointAnchor {bool defined=false;double curvature=0,normalStress=0;};
struct AdjointExternal
{
    double sigma=0;
    const Interface::HermiteGalsSurfaceVectorView* source=nullptr;
    double* stressScale=nullptr;
};
void addAdjointRule(AdjointRule& target,const AdjointRule& source)
{
    target.area+=source.area;target.moment+=source.moment;target.stressMoment+=source.stressMoment;target.nodes+=source.nodes;
    for(const auto& entry:source.native){auto& a=target.native[entry.first];const auto& b=entry.second;
        a.normal+=b.normal;a.traction+=b.traction;a.normalScale+=b.normalScale;a.tractionScale+=b.tractionScale;
        a.tangent+=b.tangent;a.tangentScale+=b.tangentScale;a.sourceScale+=b.sourceScale;
        if(a.nodes==0){a.minimum=b.minimum;a.maximum=b.maximum;}
        else{a.minimum=std::min(a.minimum,b.minimum);a.maximum=std::max(a.maximum,b.maximum);}
        a.nodes+=b.nodes;a.maxDistance=std::max(a.maxDistance,b.maxDistance);
    }
}
void chargeAdjointBudget(std::size_t nodes,std::size_t& requested,const HydroCsfCurvatureConfig& cfg)
{
    require(requested<=cfg.maxSurfaceQueriesPerFace&&nodes<=cfg.maxSurfaceQueriesPerFace-requested,
        "adjoint surface query limit");requested+=nodes;
}
AdjointRule adjointRule(const HydroPhaseGeometry& phase,const Detail::HydroRegularAdvectorGeometry& regular,
    const Detail::HydroRegularCellStencil& stencil,const std::vector<Vector2D>& centers,
    std::size_t c,std::size_t g,double begin,double end,AdjointAnchor& anchor,
    std::size_t& queries,std::size_t& requested,double& curvatureScale,const HydroCsfCurvatureConfig& cfg,
    ProjectedRateBasis* projected,const AdjointExternal* external=nullptr)
{
    const auto samples=phase.publication().quadrature(c,g,begin,end);
    chargeAdjointBudget(samples.size(),requested,cfg);queries+=samples.size();
    AdjointRule result;const auto clock=phase.timeLayer();
    for(const auto& q:samples){
        const auto& frame=q.sample.frame;const auto time=q.sample.timeLayer;
        require(q.sample.geometryIdentity==phase.geometryIdentity()&&time.step==clock.step&&time.begin==clock.begin
            &&time.end==clock.end&&time.position==clock.position&&frame.nativeCellIndex==phase.layout().nativeCell(stencil.pressure),
            "adjoint quadrature lost source publication/full clock");
        require(q.physicalWeight>0&&std::isfinite(q.physicalWeight)&&std::isfinite(q.sample.curvature)
            &&std::isfinite(frame.normal.x)&&std::isfinite(frame.normal.y),"invalid adjoint physical sample");
        const bool first=!anchor.defined;
        if(first){anchor.curvature=q.sample.curvature;anchor.defined=true;}
        const long double w=q.physicalWeight,delta=static_cast<long double>(q.sample.curvature)-anchor.curvature;
        result.area+=w;result.moment+=w*delta;++result.nodes;
        curvatureScale=std::max(curvatureScale,std::abs(q.sample.curvature));
        std::map<int,long double> normalWeights;
        std::map<int,std::array<long double,2>> tangentWeights; // contraction, arithmetic family
        std::map<int,long double> sourceWeights; // actual full vector reconstruction family
        if(projected)normalWeights=projected->normal(frame);
        if(!projected||external){
            for(const auto& row:regular.cellJetWeights(stencil,frame.point)){
                if(!projected)normalWeights[row.velocityDof]=static_cast<long double>(frame.normal.x)*row.jet[0]
                    +static_cast<long double>(frame.normal.y)*row.jet[1];
                if(external){
                    const long double x=-static_cast<long double>(frame.normal.y)*row.jet[0],
                        y=static_cast<long double>(frame.normal.x)*row.jet[1];
                    tangentWeights[row.velocityDof]={{x+y,std::abs(x)+std::abs(y)}};
                    sourceWeights[row.velocityDof]=std::abs(static_cast<long double>(row.jet[0]))
                        +std::abs(static_cast<long double>(row.jet[1]));
                }
            }
        }
        long double stressDelta=0,tangentialTraction=0;
        if(external){
            const auto vector=external->source->sample(c,g,q.sample);
            const double normalTraction=vector.x*frame.normal.x+vector.y*frame.normal.y;
            const double tangentTraction=-vector.x*frame.normal.y+vector.y*frame.normal.x;
            const double stress=external->sigma*q.sample.curvature-normalTraction;
            require(std::isfinite(normalTraction)&&std::isfinite(tangentTraction)&&std::isfinite(stress),
                "nonfinite coeval interface stress");
            if(first)anchor.normalStress=stress;
            stressDelta=static_cast<long double>(stress)-anchor.normalStress;tangentialTraction=tangentTraction;
            result.stressMoment+=w*stressDelta;
            *external->stressScale=std::max({*external->stressScale,std::abs(stress),std::abs(anchor.normalStress),
                std::abs(external->sigma*q.sample.curvature),std::abs(normalTraction),std::abs(tangentTraction)});
            // A tangent may act on a native column with zero normal motion.
            for(const auto& row:tangentWeights)normalWeights.emplace(row.first,0);
        }
        for(const auto& row:normalWeights){
            require(row.first>=0&&static_cast<std::size_t>(row.first)<centers.size(),"adjoint source DOF out of range");
            long double tangent=0,tangentScale=0;
            if(external){
                const auto it=tangentWeights.find(row.first);
                if(it!=tangentWeights.end()){
                    tangent=w*tangentialTraction*it->second[0];
                    tangentScale=w*std::abs(tangentialTraction)*it->second[1];
                }
            }
            long double sourceScale=0;
            if(external){
                const auto it=sourceWeights.find(row.first);
                sourceScale=w*((it==sourceWeights.end()?0:it->second)+(projected?std::abs(row.second):0));
                require(std::isfinite(sourceScale)&&sourceScale>=0,"nonfinite interface vector source family");
            }
            if(row.second==0&&tangent==0&&tangentScale==0&&sourceScale==0)continue;
            const long double moment=w*row.second;
            // With no external source this is the original capillary arithmetic.
            const long double traction=external?stressDelta*moment-tangent:delta*moment;
            require(std::isfinite(moment)&&std::isfinite(traction)&&std::isfinite(tangent)&&std::isfinite(tangentScale),
                "nonfinite adjoint native traction");
            auto& value=result.native[row.first];
            value.normal+=moment;value.traction+=traction;value.normalScale+=std::abs(moment);
            value.tractionScale+=external?std::abs(stressDelta*moment)+tangentScale:std::abs(traction);
            value.tangent+=tangent;value.tangentScale+=tangentScale;value.sourceScale+=sourceScale;
            if(value.nodes==0)value.minimum=value.maximum=q.sample.curvature;
            else{value.minimum=std::min(value.minimum,q.sample.curvature);value.maximum=std::max(value.maximum,q.sample.curvature);}
            ++value.nodes;const auto center=centers[row.first];
            value.maxDistance=std::max(value.maxDistance,std::hypot(frame.point.x-center.x,frame.point.y-center.y));
        }
    }
    require(result.area>0&&std::isfinite(result.area)&&std::isfinite(result.moment)&&std::isfinite(result.stressMoment),
        "invalid adjoint positive surface measure");
    return result;
}

AdjointPanel adjointPanel(const HydroPhaseGeometry& phase,const Detail::HydroRegularAdvectorGeometry& regular,
    const Detail::HydroRegularCellStencil& stencil,const std::vector<Vector2D>& centers,
    std::size_t c,std::size_t g,double begin,double end,int depth,AdjointAnchor& anchor,
    std::size_t& queries,std::size_t& requested,double& curvatureScale,const HydroCsfCurvatureConfig& cfg,
    const AdjointRule* retainedCoarse=nullptr,ProjectedRateBasis* projected=nullptr,const AdjointExternal* external=nullptr)
{
    const double middle=begin+(end-begin)/2;
    require(middle>begin&&middle<end,"unrepresentable adjoint subdivision");
    std::optional<AdjointRule> sampledCoarse;
    // Reuse is confined to the exact parent half in this same owner/anchor.
    // Charge the original logical budget, but report only actual new queries.
    if(retainedCoarse)chargeAdjointBudget(retainedCoarse->nodes,requested,cfg);
    else sampledCoarse=adjointRule(phase,regular,stencil,centers,c,g,begin,end,anchor,queries,requested,curvatureScale,cfg,projected,external);
    const auto& coarse=retainedCoarse?*retainedCoarse:*sampledCoarse;
    auto left=adjointRule(phase,regular,stencil,centers,c,g,begin,middle,anchor,queries,requested,curvatureScale,cfg,projected,external);
    auto right=adjointRule(phase,regular,stencil,centers,c,g,middle,end,anchor,queries,requested,curvatureScale,cfg,projected,external);
    auto fine=left;addAdjointRule(fine,right);
    AdjointPanel result;result.component=c;result.graph=g;result.pressure=stencil.pressure;
    result.begin=begin;result.end=end;result.depth=depth;
    result.areaError=std::abs(fine.area-coarse.area);result.momentError=std::abs(fine.moment-coarse.moment);
    result.stressMomentError=std::abs(fine.stressMoment-coarse.stressMoment);
    for(const auto& entry:coarse.native)result.errors[entry.first]={std::abs(entry.second.normal),std::abs(entry.second.traction),std::abs(entry.second.tangent)};
    for(const auto& entry:fine.native){const auto found=coarse.native.find(entry.first);
        const long double n=found==coarse.native.end()?0:found->second.normal,t=found==coarse.native.end()?0:found->second.traction,
            tangent=found==coarse.native.end()?0:found->second.tangent;
        result.errors[entry.first]={std::abs(entry.second.normal-n),std::abs(entry.second.traction-t),std::abs(entry.second.tangent-tangent)};
    }
    result.fine=std::move(fine);result.halves={std::move(left),std::move(right)};return result;
}
struct AdjointGraphBudget
{
    long double area=0,areaError=0,momentError=0,stressMomentError=0;
    std::map<int,std::array<long double,7>> native; // E_N,E_full,abs N,abs full,E_tangent,abs tangent,full source family
};
long double adjointRelativeError(long double error,long double scale)
{return scale>0?error/scale:error==0?0:std::numeric_limits<long double>::infinity();}
AdjointGraphBudget adjointGraphBudget(const std::vector<AdjointPanel>& panels)
{
    AdjointGraphBudget result;
    for(const auto& panel:panels){
        result.area+=panel.fine.area;result.areaError+=panel.areaError;result.momentError+=panel.momentError;
        result.stressMomentError+=panel.stressMomentError;
        for(const auto& entry:panel.fine.native){auto& sum=result.native[entry.first];
            sum[2]+=entry.second.normalScale;sum[3]+=entry.second.tractionScale;sum[5]+=entry.second.tangentScale;
            sum[6]+=entry.second.sourceScale;}
        for(const auto& entry:panel.errors){auto& sum=result.native[entry.first];sum[0]+=entry.second[0];sum[1]+=entry.second[1];sum[4]+=entry.second[2];}
    }return result;
}
long double adjointGraphScore(const AdjointGraphBudget& budget,double scale,double stressScale=-1)
{
    long double result=std::max(adjointRelativeError(budget.areaError,budget.area),
        adjointRelativeError(budget.momentError,scale*budget.area));
    const double nativeScale=stressScale>=0?stressScale:scale;
    if(stressScale>=0)result=std::max(result,adjointRelativeError(budget.stressMomentError,stressScale*budget.area));
    for(const auto& entry:budget.native){const auto& value=entry.second;
        // External traction is one vector family. A mathematically zero
        // normal column must not normalize its error by its own roundoff.
        // The legacy capillary path keeps its original normal-moment scale.
        const long double sourceScale=stressScale>=0?value[6]:value[2];
        result=std::max({result,adjointRelativeError(value[0],sourceScale),
            adjointRelativeError(value[1],value[3]+nativeScale*sourceScale)});
        if(stressScale>=0)result=std::max(result,adjointRelativeError(value[4],value[5]+stressScale*sourceScale));
    }
    return result;
}
long double adjointPanelPriority(const AdjointPanel& panel,const AdjointGraphBudget& budget,double scale,double stressScale=-1)
{
    // Absolute panel contributions normalized by the COMPLETE original graph
    // budget. A tiny near-zero subarc never acquires its own relative target.
    long double result=adjointRelativeError(panel.areaError,budget.area)
        +adjointRelativeError(panel.momentError,scale*budget.area);
    const double nativeScale=stressScale>=0?stressScale:scale;
    if(stressScale>=0)result+=adjointRelativeError(panel.stressMomentError,stressScale*budget.area);
    for(const auto& entry:panel.errors){const auto& value=budget.native.at(entry.first);
        const long double sourceScale=stressScale>=0?value[6]:value[2];
        result+=adjointRelativeError(entry.second[0],sourceScale)
            +adjointRelativeError(entry.second[1],value[3]+nativeScale*sourceScale);
        if(stressScale>=0)result+=adjointRelativeError(entry.second[2],value[5]+stressScale*sourceScale);
    }
    return result;
}

std::vector<double> adjointChosenCoordinateCuts(const HydroPhaseGeometry& phase,
    const Detail::HydroRegularCell& cell,std::size_t c,std::size_t g)
{
    const auto& curve=phase.publication().topology().curves[c];const auto& graph=curve.graphs()[g];
    const int axis=graph.coordinateIsX?0:1;
    const auto begin=LevelSetHermiteGeometry::cellChartPoint(curve.frameAt(static_cast<int>(g),0));
    const auto end=LevelSetHermiteGeometry::cellChartPoint(curve.frameAt(static_cast<int>(g),1));
    const auto first=axis==0?begin.x:begin.y,last=axis==0?end.x:end.y;
    const bool ascending=LAMP::Geometry::compareCellChartCoordinate(first,last)<0;
    std::vector<double> cuts{0,1};
    for(double knot:cell.supportKnots[axis]){
        const double lower=axis==0?cell.lower.x:cell.lower.y,upper=axis==0?cell.upper.x:cell.upper.y;
        if(!(knot>lower&&knot<upper))continue;
        const LAMP::Geometry::CellChartCoordinate coordinate{(knot-lower)/(upper-lower),false};
        const auto low=ascending?first:last,high=ascending?last:first;
        if(LAMP::Geometry::compareCellChartCoordinate(coordinate,low)<=0
            ||LAMP::Geometry::compareCellChartCoordinate(coordinate,high)>=0)continue;
        // This only inverts the retained chosen coordinate. The transverse
        // dummy chart is not a curve sample or a reconstructed interface.
        LAMP::Geometry::CellChartPoint chart{{.5,false},{.5,false}};
        if(graph.nativeSide)chart=begin;
        (axis==0?chart.x:chart.y)=coordinate;
        const double fraction=curve.graphFraction(static_cast<int>(g),chart);
        require(fraction>0&&fraction<1,"unrepresentable chosen-coordinate source knot");
        cuts.push_back(fraction);
    }
    std::sort(cuts.begin(),cuts.end());
    for(std::size_t i=1;i<cuts.size();++i)require(cuts[i]>cuts[i-1],"distinct chosen-coordinate source knots collapsed");
    return cuts;
}
std::vector<int> adjointPressureComponents(const QuadtreeMesh& mesh,const HydroPhaseGeometry& phase,int& count)
{
    const auto& l=phase.layout();std::vector<int> parent(l.pressureCount());std::iota(parent.begin(),parent.end(),0);
    const auto root=[&](int p){int r=p;while(parent[r]!=r)r=parent[r];
        while(parent[p]!=p){const int next=parent[p];parent[p]=r;p=next;}return r;};
    for(int u=0;u<l.velocityCount();++u){const auto f=l.domain().face(mesh,l.nativeFace(u));
        if(f.negativeCell>=0&&f.positiveCell>=0&&phase.metrics().faceMeasures()[u]>0)
            parent[root(l.pressureDof(f.positiveCell))]=root(l.pressureDof(f.negativeCell));
    }
    std::map<int,int> numbered;std::vector<int> result(parent.size());
    for(std::size_t p=0;p<parent.size();++p){const int r=root(static_cast<int>(p));
        if(!numbered.count(r))numbered[r]=static_cast<int>(numbered.size());result[p]=numbered[r];}
    count=static_cast<int>(numbered.size());return result;
}


// A height is a read-only column measure of this one retained scalar/phase
// owner. Cached publication normals choose directions; no raw curvature is
// consumed and no new curve sampling, root discovery or surface rule is made.
std::string localPlicSourceDescription(int native,const Interface::GeometricVofPlicCell& cell)
{
    std::ostringstream out;out<<std::setprecision(std::numeric_limits<long double>::max_digits10)
        <<"native="<<native<<", inventory=("<<cell.inventory[0]<<','<<cell.inventory[1]
        <<"), measured=("<<cell.measured[0]<<','<<cell.measured[1]
        <<"), normal=("<<cell.normal.x<<','<<cell.normal.y<<"), offset="<<cell.offset<<", resolved="<<cell.resolved;
    return out.str();
}
// Topology ordering only: raw column volumes/heights are never clipped.
constexpr long double localPlicHeightTopologyTolerance=1e-6L;
struct LocalPlicHeightGrid
{
    struct Slice {std::size_t first=0,last=0;long double capacity=0,liquid=0,gas=0,normalLower=0,normalUpper=0;bool reflected=false,southReflected=false,plannedComplete=false;};
    struct Direction {
        std::array<Slice,27> slices;
        std::array<std::array<Slice,12>,2> extension;
        std::array<int,2> rings{{0,0}};
        std::string reason,extendedReason;
        bool extendedReady=false,wallBlocked=false;
    };
    struct Frame {Vector2D center{0,0};double hx=0,hy=0;bool coarseAligned=false;
        std::optional<LAMP::Mesh::DyadicCellAddress> dyadic;};
    struct RefinedSupport {Frame frame;std::array<Direction,2> directions;};
    const HydroPhaseGeometry& phase;HydroCsfCurvatureReport& report;
    std::map<int,RefinedSupport> refined;
    struct FitResult {double curvature=0;std::size_t segments=0;std::string reason;bool ready=false,mixedHeightPoints=false,twoSegmentCircle=false;};
    std::map<int,FitResult> fits;
    std::vector<std::optional<HydroCsfHeightJet>> jets;
    std::vector<HydroCsfCellCurvatureMethod> methods;
    std::vector<bool> subgrid,geometryUnsupported;
    bool rz=false,axisSymmetry=false,nativeSupportChart=false;
    const HydroVofSouthEven90Boundary* south=nullptr;
    std::optional<Interface::GeometricVofPublication::PinCertificate> pin;
    std::set<int> pinCoincidentRows;
    std::unique_ptr<CartesianLeafLocator> declaredLocator;
    std::optional<LAMP::Mesh::DyadicGridCoordinates> dyadicCoordinates;
    std::optional<Geometry::PlanarCellCut::Rectangle> chartQueryBox(
        const QuadtreeMesh& mesh,const Frame& frame,int axis,int line,int column) const
    {
        if(!frame.dyadic||!dyadicCoordinates)return std::nullopt;const auto address=*frame.dyadic;
        constexpr auto limit=std::numeric_limits<std::uint64_t>::max();
        if(address.level>=64||static_cast<std::uint64_t>(mesh.Nx)>(limit>>address.level)||
            static_cast<std::uint64_t>(mesh.Ny)>(limit>>address.level))return std::nullopt;
        const auto nx=static_cast<std::uint64_t>(mesh.Nx)<<address.level,ny=static_cast<std::uint64_t>(mesh.Ny)<<address.level;
        const int tangent=1-axis;const double dh=axis?frame.hy:frame.hx,dt=tangent?frame.hy:frame.hx;
        Vector2D center=frame.center;(axis?center.y:center.x)+=line*dh;(tangent?center.y:center.x)+=column*dt;
        const Vector2D half=axis?Vector2D{dt/2,dh/2}:Vector2D{dh/2,dt/2};
        Geometry::PlanarCellCut::Rectangle result{center-half,center+half};bool any=false;
        const auto endpoint=[&](int direction,std::uint64_t base,int delta,std::uint64_t extent,double& value){
            std::uint64_t index=0;bool reflected=false;
            if(delta<0){const auto magnitude=static_cast<unsigned>(-delta);
                if(base<magnitude){
                    // Only the already certified real r=0 axis has an image.
                    if(direction!=0||!axisSymmetry||physicalLower.x!=0)return false;
                    index=magnitude-base;reflected=true;
                }else index=base-magnitude;
            }else{if(static_cast<unsigned>(delta)>extent||base>extent-static_cast<unsigned>(delta))return false;
                index=base+static_cast<unsigned>(delta);}
            if(index>extent)return false;
            const auto point=direction?dyadicCoordinates->point(0,index,address.level):dyadicCoordinates->point(index,0,address.level);
            value=direction?point.y:point.x;if(reflected)value=-value;return true;
        };
        // One invalid coordinate cannot discard a valid coordinate's exact
        // physical endpoints. Other exterior coordinates keep the old path.
        for(int direction=0;direction<2;++direction){const int delta=direction?(axis?line:column):(axis?column:line);
            const auto base=direction?address.j:address.i,extent=direction?ny:nx;double lo=0,hi=0;
            if(endpoint(direction,base,delta,extent,lo)&&endpoint(direction,base,delta+1,extent,hi)){
                (direction?result.lower.y:result.lower.x)=lo;(direction?result.upper.y:result.upper.x)=hi;any=true;}
        }
        return any?std::optional<Geometry::PlanarCellCut::Rectangle>{result}:std::nullopt;
    }
    Geometry::PlanarCellCut::Rectangle fixedQueryBox(
        const QuadtreeMesh& mesh,const Frame& frame,int axis,int line,int column) const
    {
        if(const auto actual=chartQueryBox(mesh,frame,axis,line,column))return *actual;
        // Unauthenticated frames retain their original arithmetic.
        const int tangent=1-axis;const double dh=axis?frame.hy:frame.hx,dt=tangent?frame.hy:frame.hx;
        Vector2D center=frame.center;(axis?center.y:center.x)+=line*dh;(tangent?center.y:center.x)+=column*dt;
        const Vector2D half=axis?Vector2D{dt/2,dh/2}:Vector2D{dh/2,dt/2};
        return {center-half,center+half};
    }
    Vector2D physicalLower{0,0},physicalUpper{0,0};
    static bool coversSouth(const HydroVofSouthEven90Boundary& b,double lo,double hi){
        for(const auto& span:b.intervals())if(lo>=span[0]&&hi<=span[1])return true;return false;
    }
    // Read-only images of actual leaf boxes: no exterior material state.
    void appendDeclaredSlice(const QuadtreeMesh&,int axis,int outward,
        Geometry::PlanarCellCut::Rectangle box,Slice& slice,std::string& reason,
        std::vector<HydroVofBox>& batch,bool coarseAligned)
    {
        const bool interiorNative=nativeSupportChart&&!south&&!pin;
        require((south||interiorNative)&&declaredLocator,"missing declared local height source");
        slice.first=slice.last=batch.size();slice.normalLower=axis?box.lower.y:box.lower.x;slice.normalUpper=axis?box.upper.y:box.upper.x;
        if(interiorNative&&(box.lower.x<physicalLower.x||box.lower.y<physicalLower.y||box.upper.x>physicalUpper.x||box.upper.y>physicalUpper.y)){
            reason="native coincident height support reaches physical boundary";return;}
        // Radial-height columns truncate at the true axis. A face-centred
        // chart can have a half-slice straddling r=0; its physical origin is 0.
        if(rz&&axisSymmetry&&axis==0)slice.normalLower=std::max(slice.normalLower,static_cast<long double>(physicalLower.x));
        struct Image {Geometry::PlanarCellCut::Rectangle box;bool south=false,axis=false;};std::vector<Image> images,positive;
        if(box.lower.y<physicalLower.y){
            require(coversSouth(*south,std::max(box.lower.x,physicalLower.x),box.upper.x),"S Even90 height query has partial/undeclared physical segment coverage");
            if(box.upper.y>physicalLower.y)images.push_back({{{box.lower.x,physicalLower.y},box.upper},false,false});
            const auto end=std::min(box.upper.y,physicalLower.y);
            images.push_back({{{box.lower.x,2*physicalLower.y-end},{box.upper.x,2*physicalLower.y-box.lower.y}},true,false});
            slice.southReflected=true;++report.southEven90ReflectedSlices;
        }else images.push_back({box,false,false});
        for(const auto& image:images){auto piece=image.box;
            if(piece.lower.x<physicalLower.x){
                if(!(rz&&axisSymmetry&&physicalLower.x==0)){reason="declared local height support reaches undeclared physical wall";return;}
                if(axis==0){if(piece.upper.x<=0)continue;piece.lower.x=0;positive.push_back({piece,image.south,false});}
                else if(piece.upper.x<=0){positive.push_back({{{-piece.upper.x,piece.lower.y},{-piece.lower.x,piece.upper.y}},image.south,true});slice.reflected=true;}
                else{positive.push_back({{{0,piece.lower.y},piece.upper},image.south,false});positive.push_back({{{0,piece.lower.y},{-piece.lower.x,piece.upper.y}},image.south,true});slice.reflected=true;}
            }else positive.push_back(image);
        }
        const auto& local=phase.vofPublication().plicCells();long double covered=0;const double roundoff=64*std::numeric_limits<double>::epsilon();
        for(const auto& image:positive){const auto& piece=image.box;
            if(piece.lower.x<physicalLower.x||piece.lower.y<physicalLower.y||piece.upper.x>physicalUpper.x||piece.upper.y>physicalUpper.y){reason="declared/reflected local height support reaches undeclared physical wall";return;}
            slice.capacity+=static_cast<long double>(piece.upper.x-piece.lower.x)*(piece.upper.y-piece.lower.y);
            auto ids=declaredLocator->intersectingBox(piece.lower,piece.upper);std::sort(ids.begin(),ids.end());ids.erase(std::unique(ids.begin(),ids.end()),ids.end());
            for(int ci:ids){if(phase.layout().pressureDof(ci)<0)continue;const auto& c=local.at(ci);const auto& cb=c.rectangle;
                Geometry::PlanarCellCut::Rectangle clip{{std::max(piece.lower.x,cb.lower.x),std::max(piece.lower.y,cb.lower.y)},{std::min(piece.upper.x,cb.upper.x),std::min(piece.upper.y,cb.upper.y)}};
                if(clip.lower.x>=clip.upper.x||clip.lower.y>=clip.upper.y)continue;
                if(c.measured[0]>0&&c.measured[1]>0&&!(interiorNative&&subgrid[ci])){double cn=axis?c.normal.y:c.normal.x;if((axis==1&&image.south)||(axis==0&&image.axis))cn=-cn;
                    if(rz){const double spacing=axis?box.upper.x-box.lower.x:box.upper.y-box.lower.y,sourceSpacing=axis?cb.upper.x-cb.lower.x:cb.upper.y-cb.lower.y;
                        const bool partial=axis?(clip.lower.x!=cb.lower.x||clip.upper.x!=cb.upper.x):(clip.lower.y!=cb.lower.y||clip.upper.y!=cb.upper.y);
                        if((sourceSpacing>spacing*(1+roundoff)||(coarseAligned&&sourceSpacing>=spacing*(1-roundoff)))&&partial){reason="RZ tangent height strip subdivides a coarser/unaligned mixed PLIC leaf: "+localPlicSourceDescription(ci,c);return;}}
                    if(!c.resolved||outward*cn<=0){reason="declared local column has unresolved/oppositely oriented PLIC branch: "+localPlicSourceDescription(ci,c);return;}}
                covered+=static_cast<long double>(clip.upper.x-clip.lower.x)*(clip.upper.y-clip.lower.y);batch.push_back({ci,clip});
            }
        }
        slice.last=batch.size();if(slice.capacity==0&&rz&&axisSymmetry&&axis==0&&box.upper.x<=0)return;
        if(!(slice.capacity>0&&std::abs(covered-slice.capacity)/slice.capacity<=1e-12L))reason="incomplete declared actual-leaf height-box coverage";
    }

    // Only measured mixed cells consume height directions. Pure-phase rows
    // retain their public NaN/unavailable outputs without a 102-Slice workspace.
    std::vector<int> directionRow;
    std::vector<std::array<Direction,2>> directions;
    std::array<Direction,2>& cellDirections(int p)
    {return directions[static_cast<std::size_t>(directionRow[static_cast<std::size_t>(p)])];}
    const std::array<Direction,2>& cellDirections(int p) const
    {return directions[static_cast<std::size_t>(directionRow[static_cast<std::size_t>(p)])];}
    std::vector<std::array<Frame,2>> primaryFrames;
    std::vector<Vector2D> normal;std::vector<bool> mixed;
    LocalPlicHeightGrid(const QuadtreeMesh& mesh,const HydroPhaseGeometry& source,HydroCsfCurvatureReport& r,
        const HydroVofSouthEven90Boundary* boundary=nullptr)
        :phase(source),report(r),south(boundary)
    {
        const auto& l=phase.layout();const auto& local=phase.vofPublication().plicCells();const int np=l.pressureCount();
        pin=phase.vofPublication().pinnedCorner();if(pin)pin->validateFor(mesh,phase.vofPublication());
        auto resolution=phase.vofPublication().localPlicResolutionMasks(mesh);
        subgrid=std::move(resolution.capillarySubgrid);
        geometryUnsupported=std::move(resolution.geometryUnsupported);
        directionRow.assign(np,-1);primaryFrames.resize(np);normal.resize(np);mixed.resize(np);jets.resize(np);methods.assign(np,HydroCsfCellCurvatureMethod::Unavailable);rz=l.domain().axisymmetric();
        nativeSupportChart=!rz&&phase.vofPublication().localPlicResolution()==Interface::GeometricVofLocalResolution::NativeInterfaceSupport;
        Vector2D lower{std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity()},upper{-lower.x,-lower.y};
        double hx=lower.x,hy=lower.y;int directionCount=0;
        for(int p=0;p<np;++p){const auto& c=local.at(l.nativeCell(p));const auto& b=c.rectangle;
            const double dx=b.upper.x-b.lower.x,dy=b.upper.y-b.lower.y;
            require(std::isfinite(dx)&&std::isfinite(dy)&&dx>0&&dy>0,"local PLIC height-function invalid native rectangle");
            hx=std::min(hx,dx);hy=std::min(hy,dy);
            lower.x=std::min(lower.x,b.lower.x);lower.y=std::min(lower.y,b.lower.y);upper.x=std::max(upper.x,b.upper.x);upper.y=std::max(upper.y,b.upper.y);
            mixed[p]=c.measured[0]>0&&c.measured[1]>0;normal[p]=c.normal;
            if(mixed[p])directionRow[p]=directionCount++;
            for(auto& frame:primaryFrames[p])frame={(b.lower+b.upper)*.5,dx,dy,false};
        }
        directions.resize(static_cast<std::size_t>(directionCount));
        physicalLower=lower;physicalUpper=upper;
        if(pin||np!=mesh.cellIndexNum())declaredLocator=std::make_unique<CartesianLeafLocator>(mesh);
        if(south){south->validateFor(mesh,phase.vofPublication());require(south->coordinate()==lower.y,"different actual south plane");
            declaredLocator=std::make_unique<CartesianLeafLocator>(mesh);}
        // Qualify the actual symmetry line, never an annular physical wall.
        if(rz&&lower.x==0){
            for(const auto face:l.domain().faces(mesh))if(!face.horizontal&&face.center.x==0&&
                (face.negativeCell<0||face.positiveCell<0)){
                require(face.area==0,"local PLIC axis must have exact zero physical face measure");axisSymmetry=true;}
            require(axisSymmetry,"local PLIC r=0 support lacks actual symmetry-axis native faces");
        }
        const double roundoff=64*std::numeric_limits<double>::epsilon();
        const auto canonicalSize=[&](const Cell& cell,const Geometry::PlanarCellCut::Rectangle& rectangle){
            if(!dyadicCoordinates)dyadicCoordinates.emplace(mesh.minCoord,mesh.maxCoord,mesh.Nx,mesh.Ny);
            const auto address=LAMP::Mesh::cellAddress(mesh,cell);
            require(address.level==static_cast<unsigned>(cell.level()),
                "local PLIC height-function dyadic address/actual level mismatch");
            constexpr unsigned dx[4]={0,1,1,0},dy[4]={0,0,1,1};
            for(unsigned corner=0;corner<4;++corner){
                const auto expected=dyadicCoordinates->point(address.i+dx[corner],address.j+dy[corner],address.level);
                const auto actual=mesh.verts(cell.cellToVert[corner]).coord();
                require(actual.x==expected.x&&actual.y==expected.y,
                    "local PLIC height-function requires canonical dyadic corners for rounded 2:1 widths");
                if(corner==0)require(rectangle.lower.x==actual.x&&rectangle.lower.y==actual.y,
                    "local PLIC height-function published lower corner differs from actual dyadic vertex");
                if(corner==2)require(rectangle.upper.x==actual.x&&rectangle.upper.y==actual.y,
                    "local PLIC height-function published upper corner differs from actual dyadic vertex");
            }
            return dyadicCoordinates->cellSize(address.level);
        };
        // Actual shared native subfaces certify LOCAL 2:1. Distant coarse
        // leaves need not have the globally finest size used by virtual frames.
        for(int u=0;u<l.velocityCount();++u){const auto face=l.domain().face(mesh,l.nativeFace(u));
            if(face.negativeCell<0||face.positiveCell<0)continue;
            const auto& a=mesh.cells(mesh.cellIndexToLabel[face.negativeCell]);
            const auto& b=mesh.cells(mesh.cellIndexToLabel[face.positiveCell]);
            require(std::abs(a.level()-b.level())<=1,"local PLIC height-function requires actual native-face 2:1 levels");
            const auto& ar=local[face.negativeCell].rectangle;const auto& br=local[face.positiveCell].rectangle;
            const double ax=ar.upper.x-ar.lower.x,ay=ar.upper.y-ar.lower.y;
            const double bx=br.upper.x-br.lower.x,by=br.upper.y-br.lower.y;
            if(!(std::max(ax,bx)<=2*std::min(ax,bx)*(1+roundoff)&&std::max(ay,by)<=2*std::min(ay,by)*(1+roundoff))){
                // Decimal base spacing rounds each canonical endpoint before
                // subtraction. Authenticate that exact geometry and its tree
                // chart, rather than enlarging a relative width tolerance.
                const auto as=canonicalSize(a,ar),bs=canonicalSize(b,br);
                require(std::max(as.x,bs.x)<=2*std::min(as.x,bs.x)&&std::max(as.y,bs.y)<=2*std::min(as.y,bs.y),
                    "local PLIC height-function requires actual native-face 2:1 widths");
            }
        }
        // A rounded centre/width is not an inverse of canonical endpoint
        // evaluation. Authenticate the actual native rectangle once, then name
        // fixed-query endpoints by the same integer chart as its real vertices.
        if(rz)for(int p=0;p<np;++p)if(mixed[p]){
            const auto native=l.nativeCell(p);const auto& cell=mesh.cells(mesh.cellIndexToLabel[native]);
            try{(void)canonicalSize(cell,local[native].rectangle);const auto address=LAMP::Mesh::cellAddress(mesh,cell);
                for(auto& frame:primaryFrames[p])frame.dyadic=address;
            }catch(const std::invalid_argument&){/* Noncanonical frames retain the original query path. */}
        }
        const CartesianLeafLocator locator(mesh);
        const auto candidatesFor=[&](Vector2D lo,Vector2D hi){
            auto candidates=locator.intersectingBox(lo,hi);
            // Same actual-clip accumulation order as the former std::set.
            std::sort(candidates.begin(),candidates.end());
            candidates.erase(std::unique(candidates.begin(),candidates.end()),candidates.end());
            candidates.erase(std::remove_if(candidates.begin(),candidates.end(),[&](int ci){return l.pressureDof(ci)<0;}),candidates.end());return candidates;
        };
        for(int p=0;p<np;++p){const auto& b=local[l.nativeCell(p)].rectangle;
            const bool declaredSouth=south&&b.lower.y==lower.y&&coversSouth(*south,b.lower.x,b.upper.x);
            if(mixed[p]&&declaredSouth)++report.southEven90BoundaryCells;
            if(mixed[p])require((b.lower.x>lower.x||(axisSymmetry&&b.lower.x==0))&&(b.lower.y>lower.y||declaredSouth)&&b.upper.x<upper.x&&b.upper.y<upper.y,
                "local PLIC height-function physical-contact/boundary-cell support unsupported");
        }
        // RZ-only virtual boxes. Axial-height tangent strips reflect at the
        // genuine axis; radial-height normal strips truncate there. All query
        // pieces still name real nonnegative native leaves and the same plane.
        const auto appendRzSlice=[&](int axis,int outward,Geometry::PlanarCellCut::Rectangle box,
            Slice& slice,std::string& reason,std::vector<HydroVofBox>& batch,bool coarseAligned){
            if(south){appendDeclaredSlice(mesh,axis,outward,box,slice,reason,batch,coarseAligned);return;}
            slice.first=slice.last=batch.size();
            if(box.lower.y<lower.y||box.upper.y>upper.y||box.upper.x>upper.x||
                (box.lower.x<lower.x&&!axisSymmetry)){
                reason="bounded local height support reaches physical boundary";return;}
            std::vector<Geometry::PlanarCellCut::Rectangle> pieces;
            if(axis==0){
                if(axisSymmetry&&box.upper.x<=0)return; // Empty part of a common axis-truncated column.
                box.lower.x=std::max(box.lower.x,lower.x);pieces.push_back(box);
            }else if(axisSymmetry&&box.lower.x<0){
                slice.reflected=true;
                if(box.upper.x<=0)pieces.push_back({{-box.upper.x,box.lower.y},{-box.lower.x,box.upper.y}});
                else{pieces.push_back({{0,box.lower.y},{box.upper.x,box.upper.y}});
                    pieces.push_back({{0,box.lower.y},{-box.lower.x,box.upper.y}});}
            }else pieces.push_back(box);
            slice.normalLower=axis?box.lower.y:box.lower.x;slice.normalUpper=axis?box.upper.y:box.upper.x;
            long double covered=0;
            for(const auto& piece:pieces){
                if(piece.upper.x>upper.x){reason="reflected local height support reaches physical boundary";return;}
                slice.capacity+=static_cast<long double>(piece.upper.x-piece.lower.x)*(piece.upper.y-piece.lower.y);
                for(int ci:candidatesFor(piece.lower,piece.upper)){const auto& c=local[ci];const auto& cb=c.rectangle;
                    Geometry::PlanarCellCut::Rectangle clipped{{std::max(piece.lower.x,cb.lower.x),std::max(piece.lower.y,cb.lower.y)},
                        {std::min(piece.upper.x,cb.upper.x),std::min(piece.upper.y,cb.upper.y)}};
                    if(clipped.lower.x>=clipped.upper.x||clipped.lower.y>=clipped.upper.y)continue;
                    if(c.measured[0]>0&&c.measured[1]>0){const double cn=axis?c.normal.y:c.normal.x;
                        // A fine tangent strip cannot acquire higher-order
                        // information by subdividing a coarser linear plane.
                        // This applies to primary AND supplementary frames.
                        const double spacing=axis?box.upper.x-box.lower.x:box.upper.y-box.lower.y;
                        const double sourceSpacing=axis?cb.upper.x-cb.lower.x:cb.upper.y-cb.lower.y;
                        const bool partial=axis?(clipped.lower.x!=cb.lower.x||clipped.upper.x!=cb.upper.x):
                            (clipped.lower.y!=cb.lower.y||clipped.upper.y!=cb.upper.y);
                        const bool requiresWhole=sourceSpacing>spacing*(1+roundoff)||
                            (coarseAligned&&sourceSpacing>=spacing*(1-roundoff));
                        if(requiresWhole&&partial){reason="RZ tangent height strip subdivides a coarser/unaligned mixed PLIC leaf: "+localPlicSourceDescription(ci,c);return;}
                        if(!c.resolved||outward*cn<=0){reason="bounded local column has unresolved/oppositely oriented PLIC branch: "+localPlicSourceDescription(ci,c);return;}}
                    covered+=static_cast<long double>(clipped.upper.x-clipped.lower.x)*(clipped.upper.y-clipped.lower.y);
                    batch.push_back({ci,clipped});
                }
            }
            slice.last=batch.size();
            if(!(slice.capacity>0&&std::abs(covered-slice.capacity)/slice.capacity<=1e-12L))
                reason="incomplete actual-leaf coverage of RZ meridional height box";
        };
        // One RZ geometric plan per chart, made before measuring heights or
        // computing curvature. The original native9x3 bounding box names only
        // an observed coarse level; there is no recursive stencil expansion.
        if(rz)for(int p=0;p<np;++p)if(mixed[p])for(int axis=0;axis<2;++axis){
            const auto native=l.nativeCell(p);const auto& sourceCell=mesh.cells(mesh.cellIndexToLabel[native]);
            const auto& frame=primaryFrames[p][axis];const double spacing=axis?frame.hx:frame.hy;
            const Vector2D half=axis?Vector2D{1.5*frame.hx,4.5*frame.hy}:Vector2D{4.5*frame.hx,1.5*frame.hy};
            auto lo=frame.center-half,hi=frame.center+half;
            const auto firstBox=chartQueryBox(mesh,frame,axis,-4,-1),lastBox=chartQueryBox(mesh,frame,axis,4,1);
            if(firstBox&&lastBox){lo=firstBox->lower;hi=lastBox->upper;}
            if(lo.x<0&&axisSymmetry){if(axis==1){hi.x=std::max(hi.x,-lo.x);lo.x=0;}else lo.x=0;}
            int targetLevel=sourceCell.level();
            for(int ci:candidatesFor(lo,hi)){const auto& candidate=local[ci];if(!(candidate.measured[0]>0&&candidate.measured[1]>0))continue;
                const auto& cb=candidate.rectangle;
                if(!(std::min(hi.x,cb.upper.x)>std::max(lo.x,cb.lower.x)&&std::min(hi.y,cb.upper.y)>std::max(lo.y,cb.lower.y)))continue;
                const auto& leaf=mesh.cells(mesh.cellIndexToLabel[ci]);const double width=axis?cb.upper.x-cb.lower.x:cb.upper.y-cb.lower.y;
                if(leaf.level()<sourceCell.level()&&width>spacing*(1+roundoff))targetLevel=std::min(targetLevel,leaf.level());
            }
            if(targetLevel==sourceCell.level())continue;
            const Cell* ancestor=&sourceCell;bool available=true;
            for(int step=sourceCell.level();step>targetLevel;--step){
                if(ancestor->parent<0){available=false;break;}ancestor=&mesh.cells(ancestor->parent);
                if(ancestor->level()!=step-1){available=false;break;}
            }
            if(!available||ancestor->level()!=targetLevel){cellDirections(p)[axis].reason="RZ coarse-aligned height chart lacks its actual ancestor";continue;}
            const auto a=mesh.verts(ancestor->cellToVert[0]).coord(),b=mesh.verts(ancestor->cellToVert[2]).coord();
            if(!(a.x<=frame.center.x&&frame.center.x<=b.x&&a.y<=frame.center.y&&frame.center.y<=b.y)){
                cellDirections(p)[axis].reason="RZ coarse-aligned ancestor does not contain source chart";continue;}
            primaryFrames[p][axis]={(a+b)*.5,b.x-a.x,b.y-a.y,true};
            try{(void)canonicalSize(*ancestor,{a,b});primaryFrames[p][axis].dyadic=LAMP::Mesh::cellAddress(mesh,*ancestor);}
            catch(const std::invalid_argument&){/* Noncanonical ancestor frames retain the original query path. */}
            ++report.heightFunctionCoarseAlignedCharts;
        }
        std::vector<HydroVofBox> queries;
        for(int p=0;p<np;++p)if(mixed[p]&&!subgrid[l.nativeCell(p)]){const auto& native=local[l.nativeCell(p)];const auto& b=native.rectangle;
            if(!native.resolved)throw std::invalid_argument("[HydroBalancedCsfLoad] local PLIC height-function requires a resolved mixed PLIC source: "+
                localPlicSourceDescription(l.nativeCell(p),native));
            for(int axis=0;axis<2;++axis){auto& direction=cellDirections(p)[axis];
                const double component=axis?normal[p].y:normal[p].x;
                if(component==0){direction.reason="zero PLIC normal in height direction";continue;}
                const int outward=component>0?1:-1;const auto& frame=primaryFrames[p][axis];
                for(int column=-1;column<=1&&direction.reason.empty();++column)for(int line=-4;line<=4&&direction.reason.empty();++line){
                    const auto box=fixedQueryBox(mesh,frame,axis,line,column);
                    if(declaredLocator&&boxHitsSolid(mesh,box))direction.wallBlocked=true;
                    if(rz){auto& slice=direction.slices[9*(column+1)+line+4];
                        appendRzSlice(axis,outward,box,slice,direction.reason,queries,frame.coarseAligned);continue;}
                    if(south){auto& slice=direction.slices[9*(column+1)+line+4];
                        appendDeclaredSlice(mesh,axis,outward,box,slice,direction.reason,queries,false);continue;}
                    if(box.lower.x<lower.x||box.lower.y<lower.y||box.upper.x>upper.x||box.upper.y>upper.y){direction.reason="local 9x3 height support reaches physical boundary";break;}
                    auto& slice=direction.slices[9*(column+1)+line+4];slice.first=queries.size();slice.capacity=static_cast<long double>(box.upper.x-box.lower.x)*(box.upper.y-box.lower.y);
                    const auto candidates=candidatesFor(box.lower,box.upper);
                    long double covered=0;
                    for(int ci:candidates){const auto& c=local[ci];const auto& cb=c.rectangle;
                        Geometry::PlanarCellCut::Rectangle clipped{{std::max(box.lower.x,cb.lower.x),std::max(box.lower.y,cb.lower.y)},
                            {std::min(box.upper.x,cb.upper.x),std::min(box.upper.y,cb.upper.y)}};
                        if(clipped.lower.x>=clipped.upper.x||clipped.lower.y>=clipped.upper.y)continue;
                        if(c.measured[0]>0&&c.measured[1]>0&&!(nativeSupportChart&&subgrid[ci])){const double cn=axis?c.normal.y:c.normal.x;
                            if(!c.resolved||outward*cn<=0){direction.reason="local column has unresolved/oppositely oriented PLIC branch: "+localPlicSourceDescription(ci,c);break;}}
                        covered+=static_cast<long double>(clipped.upper.x-clipped.lower.x)*(clipped.upper.y-clipped.lower.y);
                        queries.push_back({ci,clipped});
                    }
                    slice.last=queries.size();
                    if(direction.reason.empty()&&!(slice.capacity>0&&std::abs(covered-slice.capacity)/slice.capacity<=1e-12L))direction.reason="incomplete actual-leaf coverage of local height box";
                    if(nativeSupportChart&&direction.reason.empty())slice.plannedComplete=true;
                }
            }
        }
        // One public source/mesh validation for all actual clipped-leaf queries.
        const auto measured=rz?phase.measureVofMeridionalBoxes(mesh,queries):phase.measureVofBoxes(mesh,queries);
        for(int p=0;p<np;++p)if(mixed[p]&&!subgrid[l.nativeCell(p)])for(auto& direction:cellDirections(p))if(direction.reason.empty())for(auto& slice:direction.slices){
            if(rz&&slice.capacity==0)continue;
            long double lc=0,gc=0;
            for(std::size_t q=slice.first;q<slice.last;++q){const auto l=measured[q].liquid-lc,ln=slice.liquid+l;lc=(ln-slice.liquid)-l;slice.liquid=ln;
                const auto g=measured[q].gas-gc,gn=slice.gas+g;gc=(gn-slice.gas)-g;slice.gas=gn;++report.heightFunctionCellQueries;}
            if(!(std::isfinite(slice.liquid)&&std::isfinite(slice.gas)&&slice.liquid>=0&&slice.gas>=0&&
                std::abs(slice.liquid+slice.gas-slice.capacity)/slice.capacity<=1e-12L))direction.reason="local PLIC height-box phase capacity closure failed";
        }
        struct Pending {Slice* slice;std::string* reason;};
        struct Work {int p,axis;Direction* direction;Frame frame;};
        const auto nativeFrame=[&](int p,int axis){const auto& b=local[l.nativeCell(p)].rectangle;
            if(rz)return primaryFrames[p][axis];
            return Frame{(b.lower+b.upper)*.5,b.upper.x-b.lower.x,b.upper.y-b.lower.y};};
        // The original 9x3 construction above is unchanged. Supplementary
        // frames use the same actual clipped-leaf boxes and retained owner.
        const auto appendSlice=[&](const Work& work,int line,int column,Slice& slice,std::string& reason,
            std::vector<HydroVofBox>& batch,std::vector<Pending>& pending){
            const int axis=work.axis;const auto& frame=work.frame;
            const int outward=(axis?normal[work.p].y:normal[work.p].x)>0?1:-1;
            const auto box=fixedQueryBox(mesh,frame,axis,line,column);
            if(rz){appendRzSlice(axis,outward,box,slice,reason,batch,frame.coarseAligned);
                if(reason.empty()&&slice.capacity>0)pending.push_back({&slice,&reason});return;}
            if(south){appendDeclaredSlice(mesh,axis,outward,box,slice,reason,batch,frame.coarseAligned);
                if(reason.empty()&&slice.capacity>0)pending.push_back({&slice,&reason});return;}
            if(box.lower.x<lower.x||box.lower.y<lower.y||box.upper.x>upper.x||box.upper.y>upper.y){
                reason="bounded local height support reaches physical boundary";return;}
            slice.capacity=static_cast<long double>(box.upper.x-box.lower.x)*(box.upper.y-box.lower.y);
            const auto candidates=candidatesFor(box.lower,box.upper);
            long double covered=0;std::vector<HydroVofBox> clippedQueries;
            for(int ci:candidates){const auto& c=local[ci];const auto& cb=c.rectangle;
                Geometry::PlanarCellCut::Rectangle clipped{{std::max(box.lower.x,cb.lower.x),std::max(box.lower.y,cb.lower.y)},
                    {std::min(box.upper.x,cb.upper.x),std::min(box.upper.y,cb.upper.y)}};
                if(clipped.lower.x>=clipped.upper.x||clipped.lower.y>=clipped.upper.y)continue;
                // An owner-declared subgrid normal is not a resolved branch
                // witness. Its full actual phase measure still enters below.
                if(c.measured[0]>0&&c.measured[1]>0&&!(nativeSupportChart&&subgrid[ci])){const double cn=axis?c.normal.y:c.normal.x;
                    if(!c.resolved||outward*cn<=0){reason="bounded local column has unresolved/oppositely oriented PLIC branch: "+localPlicSourceDescription(ci,c);return;}}
                covered+=static_cast<long double>(clipped.upper.x-clipped.lower.x)*(clipped.upper.y-clipped.lower.y);
                clippedQueries.push_back({ci,clipped});
            }
            if(!(slice.capacity>0&&std::abs(covered-slice.capacity)/slice.capacity<=1e-12L)){
                reason="incomplete actual-leaf coverage of bounded local height box";return;}
            slice.first=batch.size();batch.insert(batch.end(),clippedQueries.begin(),clippedQueries.end());slice.last=batch.size();
            pending.push_back({&slice,&reason});
        };
        const auto measureBatch=[&](const std::vector<HydroVofBox>& batch,const std::vector<Pending>& pending,bool isRefined){
            if(batch.empty())return;
            // One live source/mesh validation for this entire supplementary round.
            const auto measured=rz?phase.measureVofMeridionalBoxes(mesh,batch):phase.measureVofBoxes(mesh,batch);
            for(const auto& item:pending){auto& slice=*item.slice;long double lc=0,gc=0;
                for(std::size_t q=slice.first;q<slice.last;++q){const auto l=measured[q].liquid-lc,ln=slice.liquid+l;lc=(ln-slice.liquid)-l;slice.liquid=ln;
                    const auto g=measured[q].gas-gc,gn=slice.gas+g;gc=(gn-slice.gas)-g;slice.gas=gn;++report.heightFunctionCellQueries;}
                if(isRefined)++report.heightFunctionRefinedSlices;else ++report.heightFunctionExtendedSlices;
                if(!(std::isfinite(slice.liquid)&&std::isfinite(slice.gas)&&slice.liquid>=0&&slice.gas>=0&&
                    std::abs(slice.liquid+slice.gas-slice.capacity)/slice.capacity<=1e-12L))*item.reason="bounded local PLIC height-box phase capacity closure failed";
            }
        };
        const auto extend=[&](std::vector<Work> active,bool isRefined){
            for(int round=0;round<4&&!active.empty();++round){
                std::vector<HydroVofBox> batch;std::vector<Pending> pending;std::vector<Work> measuredWork;
                for(const auto& work:active){auto& direction=*work.direction;bool needLower=false,needUpper=false;
                    const int outward=(work.axis?normal[work.p].y:normal[work.p].x)>0?1:-1;
                    if(!orderedSupport(direction,outward,needLower,needUpper,direction.extendedReason))continue;
                    for(int side=0;side<2&&direction.extendedReason.empty();++side){
                        if(!(side?needUpper:needLower))continue;
                        const int ring=direction.rings[side],line=(side?1:-1)*(5+ring);
                        for(int column=-1;column<=1&&direction.extendedReason.empty();++column)
                            appendSlice(work,line,column,direction.extension[side][3*ring+column+1],direction.extendedReason,batch,pending);
                        if(direction.extendedReason.empty())++direction.rings[side];
                    }
                    if(direction.extendedReason.empty())measuredWork.push_back(work);
                }
                measureBatch(batch,pending,isRefined);active.clear();
                for(const auto& work:measuredWork){auto& direction=*work.direction;bool needLower=false,needUpper=false;
                    const int outward=(work.axis?normal[work.p].y:normal[work.p].x)>0?1:-1;
                    if(!direction.extendedReason.empty()||!orderedSupport(direction,outward,needLower,needUpper,direction.extendedReason))continue;
                    if(!needLower&&!needUpper)direction.extendedReady=true;
                    else if(round==3)direction.extendedReason="bounded local height column lacks opposite phase ends after four rings per missing end";
                    else active.push_back(work);
                }
            }
        };
        // Cached qualification is uncounted; old preferred9 -> alternate9
        // always wins over any supplementary frame or normal-end extension.
        std::vector<bool> nativeNineUsable(np,false);std::vector<Work> nativeWork;
        for(int p=0;p<np;++p)if(mixed[p]&&!subgrid[l.nativeCell(p)]){double value=0;std::string reason;
            const int first=std::abs(normal[p].x)>=std::abs(normal[p].y)?0:1;
            if(estimate(p,first,value,reason,false)||estimate(p,1-first,value,reason,false)){nativeNineUsable[p]=true;continue;}
            for(int axis=0;axis<2;++axis){auto& direction=cellDirections(p)[axis];bool needLower=false,needUpper=false;
                const int outward=(axis?normal[p].y:normal[p].x)>0?1:-1;
                if(orderedSupport(direction,outward,needLower,needUpper,direction.extendedReason)&&(needLower||needUpper))
                    nativeWork.push_back({p,axis,&direction,nativeFrame(p,axis)});
            }
        }
        extend(std::move(nativeWork),false);
        // Only a coarse mixed source that also failed both native CT results
        // receives the original GLOBAL-fine virtual frame. Its hx/hy above
        // are measured from actual leaves across all levels; this does not
        // introduce a new locally adaptive frame or change HF arithmetic.
        std::vector<HydroVofBox> refinedQueries;std::vector<Pending> refinedPending;
        for(int p=0;p<np;++p)if(mixed[p]&&!subgrid[l.nativeCell(p)]&&!nativeNineUsable[p]){const auto& source=local[l.nativeCell(p)];const auto& b=source.rectangle;
            if(!(b.upper.x-b.lower.x>hx*(1+roundoff)||b.upper.y-b.lower.y>hy*(1+roundoff)))continue;
            const int first=std::abs(normal[p].x)>=std::abs(normal[p].y)?0:1;double value=0;std::string reason;
            if(estimateExtended(p,first,value,reason,false)||estimateExtended(p,1-first,value,reason,false))continue;
            auto& support=refined[p];const auto segment=Geometry::PlanarCellCut::interfaceSegment(b,source.normal,source.offset);
            if(!segment){for(auto& direction:support.directions)direction.reason="coarse mixed PLIC source has no representable actual interface segment";continue;}
            const long double cx=static_cast<long double>(b.lower.x)+(static_cast<long double>(b.upper.x)-b.lower.x)/2;
            const long double cy=static_cast<long double>(b.lower.y)+(static_cast<long double>(b.upper.y)-b.lower.y)/2;
            support.frame={{static_cast<double>(cx+((*segment)[0].x+(*segment)[1].x)/2),
                static_cast<double>(cy+((*segment)[0].y+(*segment)[1].y)/2)},hx,hy};
            if(!(std::isfinite(support.frame.center.x)&&std::isfinite(support.frame.center.y)&&
                support.frame.center.x>=b.lower.x&&support.frame.center.x<=b.upper.x&&support.frame.center.y>=b.lower.y&&support.frame.center.y<=b.upper.y)){
                for(auto& direction:support.directions)direction.reason="unrepresentable coarse PLIC segment-midpoint height anchor";continue;}
            for(int axis=0;axis<2;++axis){auto& direction=support.directions[axis];const double component=axis?normal[p].y:normal[p].x;
                if(component==0){direction.reason="zero PLIC normal in refined height direction";continue;}
                const Work work{p,axis,&direction,support.frame};
                for(int column=-1;column<=1&&direction.reason.empty();++column)for(int line=-4;line<=4&&direction.reason.empty();++line)
                    appendSlice(work,line,column,direction.slices[9*(column+1)+line+4],direction.reason,refinedQueries,refinedPending);
            }
        }
        measureBatch(refinedQueries,refinedPending,true);
        std::vector<Work> refinedWork;
        for(auto& entry:refined){const int p=entry.first;auto& support=entry.second;
            for(int axis=0;axis<2;++axis){auto& direction=support.directions[axis];bool needLower=false,needUpper=false;
                const int outward=(axis?normal[p].y:normal[p].x)>0?1:-1;
                if(orderedSupport(direction,outward,needLower,needUpper,direction.extendedReason)&&!needLower&&!needUpper)direction.extendedReady=true;
            }
            const int first=std::abs(normal[p].x)>=std::abs(normal[p].y)?0:1;double value=0;std::string reason;
            if(estimateRefined(p,first,value,reason,false)||estimateRefined(p,1-first,value,reason,false))continue;
            for(int axis=0;axis<2;++axis){auto& direction=support.directions[axis];bool needLower=false,needUpper=false;
                const int outward=(axis?normal[p].y:normal[p].x)>0?1:-1;
                if(direction.extendedReason.empty()&&orderedSupport(direction,outward,needLower,needUpper,direction.extendedReason)&&(needLower||needUpper))
                    refinedWork.push_back({p,axis,&direction,support.frame});
            }
        }
        extend(std::move(refinedWork),true);
        // Snapshot only independently qualified height curvatures. A failed
        // source may borrow a local constant fit, but that fit never becomes
        // a donor for another source in this stage.
        std::vector<double> qualifiedHeight(np,std::numeric_limits<double>::quiet_NaN());
        if(nativeSupportChart)for(int p=0;p<np;++p)if(mixed[p]&&!subgrid[l.nativeCell(p)]){
            const int first=std::abs(normal[p].x)>=std::abs(normal[p].y)?0:1;double value=0;std::string reason;
            bool ready=estimate(p,first,value,reason,false)||estimate(p,1-first,value,reason,false);
            if(!ready&&!cellDirections(p)[0].wallBlocked&&!cellDirections(p)[1].wallBlocked)
                ready=estimateExtended(p,first,value,reason,false)||estimateExtended(p,1-first,value,reason,false)||
                    estimateRefined(p,first,value,reason,false)||estimateRefined(p,1-first,value,reason,false);
            if(ready&&std::isfinite(value))qualifiedHeight[p]=value;
        }
        // A declared fixed-neighbourhood fit is the last estimator only.
        // All earlier successful native/fine HF arithmetic retains priority.
        std::vector<int> fitTargets;
        for(int p=0;p<np;++p)if(mixed[p]&&!subgrid[l.nativeCell(p)]&&!nativeNineUsable[p]){double value=0;std::string reason;
            const int first=std::abs(normal[p].x)>=std::abs(normal[p].y)?0:1;
            if(estimateExtended(p,first,value,reason,false)||estimateExtended(p,1-first,value,reason,false)||
                estimateRefined(p,first,value,reason,false)||estimateRefined(p,1-first,value,reason,false))continue;
            fitTargets.push_back(p);
        }
        if(!fitTargets.empty()){
            // Actual native/subface adjacency is built once for this stage.
            // PLIC segments need not have watertight endpoints across a face.
            std::vector<std::vector<int>> adjacent(np),qualifiedPairs(nativeSupportChart?np:0);
            for(int u=0;u<l.velocityCount();++u){const auto face=l.domain().face(mesh,l.nativeFace(u));
                if(face.negativeCell<0||face.positiveCell<0)continue;
                const int a=l.pressureDof(face.negativeCell),b=l.pressureDof(face.positiveCell);
                if(a>=0&&b>=0&&mixed[a]&&mixed[b]&&!geometryUnsupported[l.nativeCell(a)]&&!geometryUnsupported[l.nativeCell(b)]){
                    adjacent[a].push_back(b);adjacent[b].push_back(a);
                    if(nativeSupportChart){const auto& trace=phase.faceTraces()[u];
                        const long double ta=face.horizontal?normal[a].x:normal[a].y,tb=face.horizontal?normal[b].x:normal[b].y;
                        if(trace.minus.liquid>0&&trace.minus.gas>0&&trace.plus.liquid>0&&trace.plus.gas>0&&ta*tb>0){
                            qualifiedPairs[a].push_back(b);qualifiedPairs[b].push_back(a);}}
                }
            }
            using Point=Geometry::PlanarCellCut::Point;using Segment=std::array<Point,2>;
            std::map<int,std::optional<Segment>> segmentCache;
            const auto actualSegment=[&](int native)->const std::optional<Segment>&{
                const auto found=segmentCache.find(native);if(found!=segmentCache.end())return found->second;
                const auto& cell=local[native];auto segment=Geometry::PlanarCellCut::interfaceSegment(cell.rectangle,cell.normal,cell.offset);
                if(segment){const long double cx=static_cast<long double>(cell.rectangle.lower.x)+(static_cast<long double>(cell.rectangle.upper.x)-cell.rectangle.lower.x)/2;
                    const long double cy=static_cast<long double>(cell.rectangle.lower.y)+(static_cast<long double>(cell.rectangle.upper.y)-cell.rectangle.lower.y)/2;
                    for(auto& point:*segment){point.x+=cx;point.y+=cy;}}
                return segmentCache.emplace(native,std::move(segment)).first->second;
            };
            // Independent native columns are a bounded fallback only. Keep
            // the original whole-Direction qualification and complete-HF
            // values intact, including its supplementary/refined frames.
            struct NativeHeightSample {Slice slice;bool valid=false;};
            using NativeHeightSamples=std::array<std::array<NativeHeightSample,27>,2>;
            std::map<int,NativeHeightSamples> nativeHeightSamples;
            if(nativeSupportChart&&!south&&!pin){
                std::vector<HydroVofBox> supplementaryQueries;
                std::vector<NativeHeightSample*> supplementarySamples;
                const auto finishSample=[&](NativeHeightSample& sample,const auto& values,bool countQueries){
                    auto& slice=sample.slice;slice.liquid=0;slice.gas=0;long double lc=0,gc=0;
                    for(std::size_t q=slice.first;q<slice.last;++q){
                        const auto liquid=values[q].liquid-lc,ln=slice.liquid+liquid;lc=(ln-slice.liquid)-liquid;slice.liquid=ln;
                        const auto gas=values[q].gas-gc,gn=slice.gas+gas;gc=(gn-slice.gas)-gas;slice.gas=gn;
                        if(countQueries)++report.heightFunctionCellQueries;
                    }
                    sample.valid=std::isfinite(slice.capacity)&&slice.capacity>0&&std::isfinite(slice.liquid)&&std::isfinite(slice.gas)&&
                        slice.liquid>=0&&slice.gas>=0&&std::abs(slice.liquid+slice.gas-slice.capacity)/slice.capacity<=1e-12L;
                };
                for(int p:fitTargets){auto& samples=nativeHeightSamples[p];
                    for(int axis=0;axis<2;++axis){const auto component=axis?normal[p].y:normal[p].x;
                        if(component==0)continue;const int outward=component>0?1:-1;const auto& frame=primaryFrames[p][axis];
                        for(int column=-1;column<=1;++column)for(int line=-4;line<=4;++line){
                            const auto index=9*(column+1)+line+4;auto& sample=samples[axis][index];
                            const auto& original=cellDirections(p)[axis].slices[index];
                            const auto box=fixedQueryBox(mesh,frame,axis,line,column);
                            if(original.plannedComplete){
                                // The original batch measured these queries
                                // even when a different slice rejected the
                                // Direction. Reuse those actual results.
                                sample.slice=original;sample.slice.normalLower=axis?box.lower.y:box.lower.x;
                                sample.slice.normalUpper=axis?box.upper.y:box.upper.x;finishSample(sample,measured,false);continue;
                            }
                            if(box.lower.x<lower.x||box.lower.y<lower.y||box.upper.x>upper.x||box.upper.y>upper.y)continue;
                            auto& slice=sample.slice;slice.capacity=static_cast<long double>(box.upper.x-box.lower.x)*(box.upper.y-box.lower.y);
                            slice.normalLower=axis?box.lower.y:box.lower.x;slice.normalUpper=axis?box.upper.y:box.upper.x;
                            long double covered=0;bool qualified=true;std::vector<HydroVofBox> pieces;
                            for(int ci:candidatesFor(box.lower,box.upper)){if(l.pressureDof(ci)<0)continue;
                                const auto& cell=local[ci];const auto& cb=cell.rectangle;
                                Geometry::PlanarCellCut::Rectangle piece{{std::max(box.lower.x,cb.lower.x),std::max(box.lower.y,cb.lower.y)},
                                    {std::min(box.upper.x,cb.upper.x),std::min(box.upper.y,cb.upper.y)}};
                                if(piece.lower.x>=piece.upper.x||piece.lower.y>=piece.upper.y)continue;
                                if(cell.measured[0]>0&&cell.measured[1]>0&&!subgrid[ci]){
                                    const auto cn=axis?cell.normal.y:cell.normal.x;
                                    if(!cell.resolved||outward*cn<=0){qualified=false;break;}
                                }
                                covered+=static_cast<long double>(piece.upper.x-piece.lower.x)*(piece.upper.y-piece.lower.y);
                                pieces.push_back({ci,piece});
                            }
                            if(!qualified||!(slice.capacity>0&&std::abs(covered-slice.capacity)/slice.capacity<=1e-12L))continue;
                            slice.first=supplementaryQueries.size();
                            supplementaryQueries.insert(supplementaryQueries.end(),pieces.begin(),pieces.end());slice.last=supplementaryQueries.size();
                            supplementarySamples.push_back(&sample);
                        }
                    }
                }
                if(!supplementaryQueries.empty()){
                    const auto supplementaryMeasured=phase.measureVofBoxes(mesh,supplementaryQueries);
                    for(auto* sample:supplementarySamples)finishSample(*sample,supplementaryMeasured,true);
                }
            }
            const auto mixedHeightFit=[&](int p,Point origin,long double nx,long double ny,long double h,FitResult& fit){
                const auto found=nativeHeightSamples.find(p);if(found==nativeHeightSamples.end())return false;
                std::vector<Point> points;const auto& sourceBox=local[l.nativeCell(p)].rectangle;
                for(int axis=0;axis<2;++axis){const auto component=axis?normal[p].y:normal[p].x;if(component==0)continue;
                    const int outward=component>0?1:-1;const auto& frame=primaryFrames[p][axis];
                    const long double dt=axis?frame.hx:frame.hy;
                    std::array<std::optional<Point>,3> axisPoints;
                    std::array<long double,3> axisGraphs{};
                    for(int column=-1;column<=1;++column){
                        const auto at=[&](int line)->const NativeHeightSample&{return found->second[axis][9*(column+1)+line+4];};
                        std::optional<int> liquidEnd,gasEnd;
                        for(int step=0;step<=4;++step){const int line=-outward*step;const auto& sample=at(line);if(!sample.valid)break;
                            if(sample.slice.gas<=localPlicHeightTopologyTolerance*sample.slice.capacity){liquidEnd=line;break;}}
                        for(int step=0;step<=4;++step){const int line=outward*step;const auto& sample=at(line);if(!sample.valid)break;
                            if(sample.slice.liquid<=localPlicHeightTopologyTolerance*sample.slice.capacity){gasEnd=line;break;}}
                        if(!liquidEnd||!gasEnd)continue;
                        long double height=0,compensation=0,previous=1;bool qualified=true;
                        for(int line=*liquidEnd;;line+=outward){const auto& sample=at(line);if(!sample.valid){qualified=false;break;}
                            const auto fraction=sample.slice.liquid/(sample.slice.liquid+sample.slice.gas);
                            if(fraction>previous+localPlicHeightTopologyTolerance){qualified=false;break;}previous=fraction;
                            const auto term=sample.slice.liquid/dt-compensation,next=height+term;compensation=(next-height)-term;height=next;
                            if(line==*gasEnd)break;
                        }
                        if(!qualified)continue;const auto& first=at(*liquidEnd).slice;
                        const long double graph=(outward>0?first.normalLower:first.normalUpper)+outward*height;
                        const long double tangent=(axis?frame.center.x:frame.center.y)+column*dt;
                        if(!(std::isfinite(graph)&&std::isfinite(tangent)))continue;
                        if(column==0){const long double lo=axis?sourceBox.lower.y:sourceBox.lower.x,hi=axis?sourceBox.upper.y:sourceBox.upper.x;
                            if(graph<lo-dt/2||graph>hi+dt/2)continue;}
                        axisGraphs[column+1]=graph;
                        axisPoints[column+1]=axis?Point{tangent,graph}:Point{graph,tangent};
                    }
                    // Nearest phase ends may give each column a different
                    // origin. Qualify their world graphs before mixing axes,
                    // retaining the complete HF chart's opposed-secant rule.
                    if(axisPoints[0]&&axisPoints[1]&&axisPoints[2]){std::string reason;
                        if(!qualifiedNativeChart(-1,frame,axis,outward,0,0,axisGraphs,reason))continue;}
                    for(const auto& candidate:axisPoints){if(!candidate)continue;const auto point=*candidate;bool independent=true;
                        for(const auto other:points)if(std::hypot((point.x-other.x)/h,(point.y-other.y)/h)<.5L){independent=false;break;}
                        if(independent)points.push_back(point);
                    }
                }
                if(points.size()<3)return false;
                using Matrix=Eigen::Matrix<long double,Eigen::Dynamic,3>;using Vector=Eigen::Matrix<long double,Eigen::Dynamic,1>;
                Matrix matrix(static_cast<Eigen::Index>(points.size()),3);Vector values(matrix.rows());bool negative=false,positive=false;
                for(std::size_t i=0;i<points.size();++i){const auto x=points[i].x-origin.x,y=points[i].y-origin.y;
                    const long double s=(-ny*x+nx*y)/h,v=(nx*x+ny*y)/h;
                    matrix(static_cast<Eigen::Index>(i),0)=1;matrix(static_cast<Eigen::Index>(i),1)=s;
                    matrix(static_cast<Eigen::Index>(i),2)=s*s;values[static_cast<Eigen::Index>(i)]=v;
                    negative=negative||s<0;positive=positive||s>0;
                }
                if(!negative||!positive||!matrix.allFinite()||!values.allFinite())return false;
                // Mixed axes must also describe a noncancelling graph in the
                // source frame. Sort only this qualification copy; retain the
                // original point and QR row order for every accepted graph.
                std::array<Point,6> graph{};
                for(std::size_t i=0;i<points.size();++i)
                    graph[i]={matrix(static_cast<Eigen::Index>(i),1),values[static_cast<Eigen::Index>(i)]};
                std::sort(graph.begin(),graph.begin()+points.size(),[](Point a,Point b){return a.x<b.x;});
                for(std::size_t i=1;i<points.size();++i)if(!(graph[i-1].x<graph[i].x))return false;
                for(std::size_t i=2;i<points.size();++i){
                    const auto left=(graph[i-1].y-graph[i-2].y)/(graph[i-1].x-graph[i-2].x);
                    const auto right=(graph[i].y-graph[i-1].y)/(graph[i].x-graph[i-1].x);
                    if(!std::isfinite(left)||!std::isfinite(right))return false;
                    if(((left<0&&right>0)||(left>0&&right<0))&&std::max(std::abs(left),std::abs(right))>1)return false;
                }
                Eigen::ColPivHouseholderQR<Matrix> qr(matrix);if(qr.rank()!=3)return false;
                const Eigen::Matrix<long double,3,1> coefficients=qr.solve(values);if(!coefficients.allFinite())return false;
                const auto metric=std::hypot(1.L,coefficients[1]),curvature=-2*coefficients[2]/(h*metric*metric*metric);
                const auto value=static_cast<double>(curvature);if(!std::isfinite(value)||(value==0&&curvature!=0))return false;
                fit.curvature=value;fit.segments=0;fit.mixedHeightPoints=true;fit.ready=true;return true;
            };
            for(int p:fitTargets){auto& fit=fits[p];const int sourceNative=l.nativeCell(p);const auto& source=local[sourceNative];
                const auto& sourceSegment=actualSegment(sourceNative);long double spanLower=0,spanUpper=0;int rank=0;
                const auto fail=[&](const std::string& why){std::ostringstream detail;detail<<std::setprecision(17)
                    <<(nativeSupportChart?"local PLIC circle fit ":"local PLIC quadratic fit ")<<why<<", segments="<<fit.segments<<", span=["<<spanLower<<','<<spanUpper<<"], rank="<<rank;fit.reason=detail.str();};
                if(!sourceSegment){fail("has no representable source segment");continue;}
                const Point origin{((*sourceSegment)[0].x+(*sourceSegment)[1].x)/2,((*sourceSegment)[0].y+(*sourceSegment)[1].y)/2};
                const long double normalLength=std::hypot(static_cast<long double>(source.normal.x),static_cast<long double>(source.normal.y));
                if(!(std::isfinite(normalLength)&&normalLength>0)){fail("has invalid source normal");continue;}
                const long double nx=source.normal.x/normalLength,ny=source.normal.y/normalLength;
                const long double h=std::max(source.rectangle.upper.x-source.rectangle.lower.x,source.rectangle.upper.y-source.rectangle.lower.y),half=2.5L*h;
                // Closed-box broad phase conservatively encloses rounding of
                // these long-double endpoints. The actual fixed5h segment clip
                // below is unchanged; no support expansion enters the fit.
                const double negativeInfinity=-std::numeric_limits<double>::infinity(),positiveInfinity=-negativeInfinity;
                const Vector2D queryLower{std::nextafter(static_cast<double>(origin.x-half),negativeInfinity),
                    std::nextafter(static_cast<double>(origin.y-half),negativeInfinity)};
                const Vector2D queryUpper{std::nextafter(static_cast<double>(origin.x+half),positiveInfinity),
                    std::nextafter(static_cast<double>(origin.y+half),positiveInfinity)};
                const auto candidates=candidatesFor(queryLower,queryUpper);
                std::map<int,Segment> clipped;
                for(int native:candidates){const int candidate=l.pressureDof(native);if(candidate<0||!mixed[candidate])continue;
                    // Missing own curvature does not invalidate an actual
                    // resolved, orientation-qualified geometric position.
                    const auto& cell=local[native];if(geometryUnsupported[native]||!cell.resolved||static_cast<long double>(cell.normal.x)*nx+static_cast<long double>(cell.normal.y)*ny<=0)continue;
                    const auto& segment=actualSegment(native);if(!segment)continue;
                    const Point a{(*segment)[0].x-origin.x,(*segment)[0].y-origin.y},b{(*segment)[1].x-origin.x,(*segment)[1].y-origin.y};
                    const Point delta{b.x-a.x,b.y-a.y};long double begin=0,end=1;bool inside=true;
                    for(int axis=0;axis<2;++axis){const long double coordinate=axis?a.y:a.x,direction=axis?delta.y:delta.x;
                        if(direction==0){if(coordinate< -half||coordinate>half){inside=false;break;}}
                        else{long double lo=(-half-coordinate)/direction,hi=(half-coordinate)/direction;if(lo>hi)std::swap(lo,hi);
                            begin=std::max(begin,lo);end=std::min(end,hi);if(!(end>begin)){inside=false;break;}}
                    }
                    if(!inside||!(end>begin))continue;
                    const Segment cut{{{a.x+begin*delta.x,a.y+begin*delta.y},{a.x+end*delta.x,a.y+end*delta.y}}};
                    if(std::hypot(cut[1].x-cut[0].x,cut[1].y-cut[0].y)>0)clipped.emplace(candidate,cut);
                }
                if(!clipped.count(p)){fail("source segment is unavailable inside fixed box");continue;}
                if(nativeSupportChart){
                    // A positive segment-length measure gives the rank-one
                    // constant-curvature least-squares fit. Its immutable HF
                    // donors use this same fixed box and normal hemisphere;
                    // it cannot introduce an extremum outside their values.
                    std::vector<double> weights,curvatures;
                    for(const auto& entry:clipped)if(std::isfinite(qualifiedHeight[entry.first])){
                        const auto& segment=entry.second;
                        const long double length=std::hypot(segment[1].x-segment[0].x,segment[1].y-segment[0].y);
                        const double weight=static_cast<double>(length/h);
                        if(!(std::isfinite(weight)&&weight>0))continue;
                        weights.push_back(weight);curvatures.push_back(qualifiedHeight[entry.first]);
                    }
                    if(!weights.empty()){
                        fit.curvature=HydroBalancedCsfLoad::positiveCurvatureAverage(weights,curvatures);
                        fit.segments=weights.size();fit.ready=true;continue;
                    }
                }
                if(mixedHeightFit(p,origin,nx,ny,h,fit))continue;
                std::vector<int> component{p};std::set<int> visited{p};
                for(std::size_t i=0;i<component.size();++i)for(int neighbour:adjacent[component[i]])
                    if(clipped.count(neighbour)&&visited.insert(neighbour).second)component.push_back(neighbour);
                fit.segments=component.size();std::set<long double> tangentPositions;bool negative=false,positive=false;
                for(int cell:component){const auto& segment=clipped.at(cell);const long double sx=(segment[0].x+segment[1].x)/2,sy=(segment[0].y+segment[1].y)/2;
                    const long double tangent=(-ny*sx+nx*sy)/h;tangentPositions.insert(tangent);spanLower=std::min(spanLower,tangent);spanUpper=std::max(spanUpper,tangent);
                    negative=negative||tangent<0;positive=positive||tangent>0;
                }
                // The sparse estimator requires the entire eligible box to
                // contain this one qualified shared-face pair, with genuine
                // source-tangent separation. It does not alter larger fits.
                bool twoSegmentCircle=false;
                if(nativeSupportChart&&component.size()==2&&clipped.size()==2&&tangentPositions.size()==2&&
                    *tangentPositions.rbegin()-*tangentPositions.begin()>=.5L){
                    const auto& pairs=qualifiedPairs[p];
                    twoSegmentCircle=std::find(pairs.begin(),pairs.end(),component[1])!=pairs.end();
                }
                // Retain every qualified two-sided face-cluster fit and its
                // original sample order. Only insufficient topology support
                // uses the already clipped, same-oriented fixed-box segments.
                if(!twoSegmentCircle&&(fit.segments<3||tangentPositions.size()<3||!negative||!positive)){
                    component.assign(1,p);std::vector<int> others;
                    for(const auto& entry:clipped)if(entry.first!=p)others.push_back(entry.first);
                    std::sort(others.begin(),others.end(),[&](int a,int b){return l.nativeCell(a)<l.nativeCell(b);});
                    component.insert(component.end(),others.begin(),others.end());
                    fit.segments=component.size();tangentPositions.clear();spanLower=0;spanUpper=0;
                    for(int cell:component){const auto& segment=clipped.at(cell);const long double sx=(segment[0].x+segment[1].x)/2,sy=(segment[0].y+segment[1].y)/2;
                        const long double tangent=(-ny*sx+nx*sy)/h;tangentPositions.insert(tangent);spanLower=std::min(spanLower,tangent);spanUpper=std::max(spanUpper,tangent);}
                }
                if(fit.segments<3&&!twoSegmentCircle){fail("requires at least three native segments");continue;}
                if(tangentPositions.size()<3&&!twoSegmentCircle){fail("lacks independent tangent midpoints");continue;}
                using Matrix=Eigen::Matrix<long double,Eigen::Dynamic,3>;using Vector=Eigen::Matrix<long double,Eigen::Dynamic,1>;
                Matrix matrix(static_cast<Eigen::Index>(3*fit.segments),3);Vector values(matrix.rows());
                const long double gauss=std::sqrt(3.L/5);const std::array<long double,3> nodes{{-gauss,0,gauss}},weights{{5.L/9,8.L/9,5.L/9}};
                Eigen::Index row=0;bool finite=true;
                for(int cell:component){const auto& segment=clipped.at(cell);const long double dx=segment[1].x-segment[0].x,dy=segment[1].y-segment[0].y;
                    const long double length=std::hypot(dx,dy);
                    for(int q=0;q<3;++q){const long double parameter=(nodes[q]+1)/2,x=segment[0].x+parameter*dx,y=segment[0].y+parameter*dy;
                        const long double s=(-ny*x+nx*y)/h,v=(nx*x+ny*y)/h,weight=std::sqrt(length*weights[q]/(2*h));
                        matrix(row,0)=weight;matrix(row,1)=weight*s;matrix(row,2)=weight*(nativeSupportChart?s*s+v*v:s*s);values[row]=weight*v;
                        finite=finite&&std::isfinite(weight)&&weight>0&&std::isfinite(s)&&std::isfinite(v);++row;
                    }
                }
                if(!finite||!matrix.allFinite()||!values.allFinite()){fail("has nonfinite weighted samples");continue;}
                Eigen::ColPivHouseholderQR<Matrix> qr(matrix);rank=static_cast<int>(qr.rank());
                if(rank!=3){fail("has deficient numerical rank");continue;}
                const Eigen::Matrix<long double,3,1> coefficients=qr.solve(values);
                if(!coefficients.allFinite()){fail("has nonfinite coefficients");continue;}
                if(nativeSupportChart){
                    // v=a+b*s+c*(s*s+v*v) fits a circle, including its planar
                    // limit c=0. Normalize its constant signed curvature by
                    // the measured RMS gradient of that implicit function.
                    // The intercept normal equation makes this equivalent to
                    // sqrt(1+b*b-4*a*c), without subtracting nearly equal terms.
                    long double gradientSquared=0,totalWeight=0;bool gradientAligned=true;
                    for(Eigen::Index q=0;q<matrix.rows();++q){
                        const long double tangent=coefficients[1]*matrix(q,0)+2*coefficients[2]*matrix(q,1);
                        const long double normal=matrix(q,0)-2*coefficients[2]*values[q];
                        gradientSquared+=tangent*tangent+normal*normal;totalWeight+=matrix(q,0)*matrix(q,0);
                        if(twoSegmentCircle){const auto& actualNormal=this->normal[component[static_cast<std::size_t>(q/3)]];
                            const long double actualT=-ny*actualNormal.x+nx*actualNormal.y,actualN=nx*actualNormal.x+ny*actualNormal.y;
                            const long double agreement=-tangent*actualT+normal*actualN;
                            gradientAligned=gradientAligned&&std::isfinite(agreement)&&agreement>0;}
                    }
                    if(!gradientAligned){fail("two-segment circle gradient opposes an actual PLIC normal");continue;}
                    if(!(std::isfinite(gradientSquared)&&gradientSquared>0&&std::isfinite(totalWeight)&&totalWeight>0)){
                        fail("has invalid implicit-circle gradient");continue;}
                    const long double curvature=-2*coefficients[2]/(h*std::sqrt(gradientSquared/totalWeight));
                    fit.curvature=static_cast<double>(curvature);
                    if(!std::isfinite(fit.curvature)||(fit.curvature==0&&curvature!=0)){fail("has unrepresentable curvature");continue;}
                    fit.twoSegmentCircle=twoSegmentCircle;fit.ready=true;continue;
                }
                const long double metric=std::hypot(1.L,coefficients[1]),curvature=-2*coefficients[2]/(h*metric*metric*metric);
                long double total=curvature;
                if(rz){
                    // The same rotated quadratic supplies point, normal and kappa.
                    // A fit alone does not certify an even jet at the true axis.
                    const long double radius=origin.x+h*coefficients[0]*nx;
                    const long double radialNormal=(nx+coefficients[1]*ny)/metric;
                    if(!(std::isfinite(radius)&&radius>0)){fail("lacks a positive-radius qualified RZ fit point");continue;}
                    total+=radialNormal/radius;
                }
                fit.curvature=static_cast<double>(total);
                if(!std::isfinite(fit.curvature)||(fit.curvature==0&&total!=0)){fail("has unrepresentable curvature");continue;}
                fit.ready=true;
            }
        }
    }
    const Slice& extendedSlice(const Direction& direction,int column,int line) const
    {
        if(line>=-4&&line<=4)return direction.slices[9*column+line+4];
        const int side=line>4?1:0,ring=std::abs(line)-5;
        return direction.extension[side][3*ring+column];
    }
    bool orderedSupport(const Direction& direction,int outward,bool& needLower,bool& needUpper,std::string& reason) const
    {
        if(!direction.reason.empty()){reason=direction.reason;return false;}
        const int lower=-4-direction.rings[0],upper=4+direction.rings[1];
        for(int column=0;column<3;++column){
            int firstLine=outward>0?lower:upper,lastLine=outward>0?upper:lower;
            if(rz){while(firstLine>=lower&&firstLine<=upper&&extendedSlice(direction,column,firstLine).capacity==0)firstLine+=outward;
                while(lastLine>=lower&&lastLine<=upper&&extendedSlice(direction,column,lastLine).capacity==0)lastLine-=outward;
                if(firstLine<lower||firstLine>upper||lastLine<lower||lastLine>upper){reason="empty axis-truncated height column";return false;}}
            const auto& first=extendedSlice(direction,column,firstLine);
            const auto& last=extendedSlice(direction,column,lastLine);
            if(first.gas>localPlicHeightTopologyTolerance*first.capacity)(outward>0?needLower:needUpper)=true;
            if(last.liquid>localPlicHeightTopologyTolerance*last.capacity)(outward>0?needUpper:needLower)=true;
            long double previous=1;
            for(int step=0;step<=upper-lower;++step){const int line=outward>0?lower+step:upper-step;
                const auto& slice=extendedSlice(direction,column,line);if(rz&&slice.capacity==0)continue;
                const auto fraction=slice.liquid/(slice.liquid+slice.gas);
                if(fraction>previous+localPlicHeightTopologyTolerance){
                    std::ostringstream detail;detail<<std::setprecision(17)<<"bounded local height column has multiple liquid-to-gas transitions beyond topology tolerance, column="<<column<<", line="<<line<<", previous="<<previous<<", fraction="<<fraction<<", reversal="<<fraction-previous;reason=detail.str();return false;}
                previous=fraction;
            }
        }
        return true;
    }
    bool finishJet(int p,const Direction& direction,const Frame& frame,int axis,int outward,int lower,int upper,
        const std::array<long double,3>& height,long double first,long double second,double& curvature,std::string& reason,bool count)
    {
        HydroCsfHeightJet jet;jet.height=height[1];jet.first=first;jet.second=second;jet.heightAxis=axis;jet.outward=outward;
        for(int column=0;column<3;++column){
            int begin=outward>0?lower:upper;
            while(begin>=lower&&begin<=upper&&extendedSlice(direction,column,begin).capacity==0)begin+=outward;
            if(begin<lower||begin>upper){reason="empty qualified RZ height column";return false;}
            const auto& endpoint=extendedSlice(direction,column,begin);
            const auto origin=outward>0?endpoint.normalLower:endpoint.normalUpper;
            if(column==0)jet.origin=origin;
            else if(origin!=jet.origin){reason="RZ height columns have different actual physical origins";return false;}
        }
        jet.tangentCoordinate=axis?frame.center.x:frame.center.y;
        jet.normalSpacing=axis?frame.hy:frame.hx;jet.tangentSpacing=axis?frame.hx:frame.hy;jet.coarseAligned=frame.coarseAligned;
        jet.axisRegular=axisSymmetry&&axis==1&&jet.tangentCoordinate==0&&first==0;
        for(int column=0;column<3;++column)for(int line=lower;line<=upper;++line){
            jet.axisReflected=jet.axisReflected||extendedSlice(direction,column,line).reflected;
            jet.southEven90Reflected=jet.southEven90Reflected||extendedSlice(direction,column,line).southReflected;}
        try{const auto parts=HydroBalancedCsfLoad::axisymmetricHeightCurvatures(jet);
            const auto total=parts[0]+parts[1];curvature=static_cast<double>(total);
            if(!std::isfinite(curvature)||(curvature==0&&total!=0)){reason="unrepresentable local PLIC RZ total height curvature";return false;}}
        catch(const std::invalid_argument& e){reason=e.what();return false;}
        if(count)jets[p]=jet;return true;
    }
    bool qualifiedNativeChart(int p,const Frame& frame,int axis,int outward,int lower,int upper,
        const std::array<long double,3>& height,std::string& reason) const
    {
        if(!nativeSupportChart)return true;
        const long double spacing=axis?frame.hx:frame.hy,dh=axis?frame.hy:frame.hx;
        const long double left=height[1]-height[0],right=height[2]-height[1];
        const long double transverseRange=std::max(spacing,std::max(std::abs(left),std::abs(right)));
        // Opposite steep secants can cancel in the centred first derivative
        // across a cell-scale cusp. Smooth oblique graphs with same-sign
        // secants retain their HF chart, including curves crossing 45 degrees.
        if(((left<0&&right>0)||(left>0&&right<0))&&transverseRange>spacing){
            reason="native-support height chart has opposed steep transverse secants";return false;}
        if(p>=0){
            const auto& box=phase.vofPublication().plicCells()[phase.layout().nativeCell(p)].rectangle;
            const long double origin=(axis?frame.center.y:frame.center.x)+(outward>0?lower-.5L:upper+.5L)*dh;
            const long double point=origin+outward*height[1],lo=axis?box.lower.y:box.lower.x,hi=axis?box.upper.y:box.upper.x;
            // The quadratic height graph has these two secants as its
            // central-cell endpoint slopes. Their maximum bounds the mean
            // height's distance to any point in that transverse column.
            if(point<lo-transverseRange/2||point>hi+transverseRange/2){
                reason="native-support height chart does not contain its source interface";return false;}
        }
        return true;
    }
    bool nativeHeightCircle(const std::array<long double,3>& height,long double spacing,double& curvature,std::string& reason) const
    {
        // Signed circumcircle curvature of the same three qualified integral
        // heights. The outer chord is at least 2h: the geometric operator is
        // intrinsically bounded by 1/h, without clipping its result.
        const auto left=(height[1]-height[0])/spacing,right=(height[2]-height[1])/spacing;
        const auto turn=(right-left)/std::hypot(1.L,left)/std::hypot(1.L,right);
        const auto k=(-2*turn/std::hypot(2.L,left+right))/spacing;
        curvature=static_cast<double>(k);
        if(!std::isfinite(curvature)||(turn==0&&right!=left)||(k==0&&turn!=0)||(curvature==0&&k!=0)){
            reason="unrepresentable native finite-height circle curvature";return false;}
        return true;
    }
    bool estimateWindow(int p,const Direction& direction,const Frame& frame,int axis,int outward,double& curvature,std::string& reason,bool count)
    {
        if(!direction.extendedReady){reason=direction.extendedReason.empty()?direction.reason:direction.extendedReason;return false;}
        const int lower=-4-direction.rings[0],upper=4+direction.rings[1];
        const long double dh=axis?frame.hy:frame.hx,spacing=axis?frame.hx:frame.hy;
        std::array<long double,3> height{};
        for(int column=0;column<3;++column){if(count)++report.heightFunctionColumns;long double compensation=0;
            for(int step=0;step<=upper-lower;++step){const auto& slice=extendedSlice(direction,column,outward>0?lower+step:upper-step);
                if(rz&&slice.capacity==0)continue;
                const auto fraction=slice.liquid/(slice.liquid+slice.gas),term=(rz?slice.liquid/spacing:fraction*dh)-compensation,next=height[column]+term;
                compensation=(next-height[column])-term;height[column]=next;}
        }
        if(!qualifiedNativeChart(p,frame,axis,outward,lower,upper,height,reason))return false;
        if(nativeSupportChart)return nativeHeightCircle(height,spacing,curvature,reason);
        // All three columns include the same physical endpoints and origin.
        const auto left=height[1]-height[0],right=height[2]-height[1];
        const auto first=(left+right)/(2*spacing),second=(right-left)/(spacing*spacing),metric=std::hypot(1.L,first);
        if(rz)return finishJet(p,direction,frame,axis,outward,lower,upper,height,first,second,curvature,reason,count);
        const auto k=-second/(metric*metric*metric);curvature=static_cast<double>(k);
        if(!std::isfinite(curvature)||(curvature==0&&k!=0)){reason="unrepresentable bounded local PLIC height curvature";return false;}
        return true;
    }
    bool estimateExtended(int p,int axis,double& curvature,std::string& reason,bool count=true)
    {
        const auto& b=phase.vofPublication().plicCells()[phase.layout().nativeCell(p)].rectangle;
        const Frame frame=rz?primaryFrames[p][axis]:Frame{(b.lower+b.upper)*.5,b.upper.x-b.lower.x,b.upper.y-b.lower.y};
        return estimateWindow(p,cellDirections(p)[axis],frame,axis,(axis?normal[p].y:normal[p].x)>0?1:-1,curvature,reason,count);
    }
    bool estimateRefined(int p,int axis,double& curvature,std::string& reason,bool count=true)
    {
        const auto found=refined.find(p);if(found==refined.end())return false;
        return estimateWindow(p,found->second.directions[axis],found->second.frame,axis,(axis?normal[p].y:normal[p].x)>0?1:-1,curvature,reason,count);
    }
    bool estimate(int p,int axis,double& curvature,std::string& reason,bool count=true)
    {
        const auto& direction=cellDirections(p)[axis];if(!direction.reason.empty()){reason=direction.reason;return false;}
        const auto& box=phase.vofPublication().plicCells()[phase.layout().nativeCell(p)].rectangle;
        const int outward=(axis?normal[p].y:normal[p].x)>0?1:-1;
        const auto& selectedFrame=primaryFrames[p][axis];
        const long double dh=rz?(axis?selectedFrame.hy:selectedFrame.hx):(axis?box.upper.y-box.lower.y:box.upper.x-box.lower.x);
        const long double spacing=rz?(axis?selectedFrame.hx:selectedFrame.hy):(axis?box.upper.x-box.lower.x:box.upper.y-box.lower.y);
        std::array<long double,3> height{};
        for(int column=0;column<3;++column){if(count)++report.heightFunctionColumns;
            int firstIndex=outward>0?0:8,lastIndex=outward>0?8:0;
            if(rz){while(firstIndex>=0&&firstIndex<9&&direction.slices[9*column+firstIndex].capacity==0)firstIndex+=outward;
                while(lastIndex>=0&&lastIndex<9&&direction.slices[9*column+lastIndex].capacity==0)lastIndex-=outward;
                if(firstIndex<0||firstIndex>=9||lastIndex<0||lastIndex>=9){reason="empty axis-truncated height column";return false;}}
            const auto& first=direction.slices[9*column+firstIndex];
            const auto& last=direction.slices[9*column+lastIndex];
            if(first.gas>localPlicHeightTopologyTolerance*first.capacity||last.liquid>localPlicHeightTopologyTolerance*last.capacity){
                std::ostringstream detail;detail<<std::setprecision(17)<<"local 9-cell column lacks opposite phase ends within topology tolerance, column="<<column<<", firstGasFraction="<<first.gas/first.capacity<<", lastLiquidFraction="<<last.liquid/last.capacity;reason=detail.str();return false;}
            long double previous=1,compensation=0;
            for(int step=0;step<9;++step){const auto& slice=direction.slices[9*column+(outward>0?step:8-step)];
                if(rz&&slice.capacity==0)continue;const auto fraction=slice.liquid/(slice.liquid+slice.gas);
                if(fraction>previous+localPlicHeightTopologyTolerance){
                    std::ostringstream detail;detail<<std::setprecision(17)<<"local height column has multiple liquid-to-gas transitions beyond topology tolerance, column="<<column<<", step="<<step<<", previous="<<previous<<", fraction="<<fraction<<", reversal="<<fraction-previous;reason=detail.str();return false;}previous=fraction;
                const auto term=(rz?slice.liquid/spacing:fraction*dh)-compensation,next=height[column]+term;compensation=(next-height[column])-term;height[column]=next;}
        }
        if(!qualifiedNativeChart(p,selectedFrame,axis,outward,-4,4,height,reason))return false;
        if(nativeSupportChart)return nativeHeightCircle(height,spacing,curvature,reason);
        // The same physical height origin cancels between all three columns.
        const auto left=height[1]-height[0],right=height[2]-height[1];
        const auto first=(left+right)/(2*spacing),second=(right-left)/(spacing*spacing),metric=std::hypot(1.L,first);
        if(rz)return finishJet(p,direction,selectedFrame,axis,outward,-4,4,height,first,second,curvature,reason,count);
        const auto k=-second/(metric*metric*metric);curvature=static_cast<double>(k);
        if(!std::isfinite(curvature)||(curvature==0&&k!=0)){reason="unrepresentable local PLIC height curvature";return false;}
        return true;
    }
    // Failure-only evidence for a direct local fixture. This method is never
    // called by a successful estimator and never contributes numerical input.
    // At most 96 actual native cells + 32 actual pure-jump faces are printed.
    std::string failureSupportDiagnostic(const QuadtreeMesh& mesh,int p) const
    {
        using Point=Geometry::PlanarCellCut::Point;using Segment=std::array<Point,2>;
        std::ostringstream out;out<<std::setprecision(std::numeric_limits<long double>::max_digits10);
        try{
            const auto& l=phase.layout();const auto& local=phase.vofPublication().plicCells();
            const int sourceNative=l.nativeCell(p);const auto& source=local.at(sourceNative);
            const auto segmentAt=[&](int native)->std::optional<Segment>{
                const auto& cell=local.at(native);auto segment=Geometry::PlanarCellCut::interfaceSegment(cell.rectangle,cell.normal,cell.offset);
                if(segment){const long double cx=static_cast<long double>(cell.rectangle.lower.x)+(static_cast<long double>(cell.rectangle.upper.x)-cell.rectangle.lower.x)/2;
                    const long double cy=static_cast<long double>(cell.rectangle.lower.y)+(static_cast<long double>(cell.rectangle.upper.y)-cell.rectangle.lower.y)/2;
                    for(auto& point:*segment){point.x+=cx;point.y+=cy;}}
                return segment;
            };
            const auto printSegment=[&](const Segment& segment){out<<"[("<<segment[0].x<<','<<segment[0].y<<"),("<<segment[1].x<<','<<segment[1].y<<")]";};
            const auto sourceSegment=segmentAt(sourceNative);
            const auto clock=phase.timeLayer();out<<"fixed5h_support{record_limit=128, cell_limit=96, pure_face_limit=32, source="<<sourceNative
                <<", vof_id="<<phase.vofPublication().identity().value()<<", clock=("<<clock.step<<','<<clock.begin<<','<<clock.end<<','<<static_cast<int>(clock.position)<<')';
            if(!sourceSegment){out<<", source_segment=none, box_unavailable=1}";return out.str();}
            const Point origin{((*sourceSegment)[0].x+(*sourceSegment)[1].x)/2,((*sourceSegment)[0].y+(*sourceSegment)[1].y)/2};
            const long double normalLength=std::hypot(static_cast<long double>(source.normal.x),static_cast<long double>(source.normal.y));
            if(!(std::isfinite(normalLength)&&normalLength>0)){out<<", source_normal_invalid=1}";return out.str();}
            const long double nx=source.normal.x/normalLength,ny=source.normal.y/normalLength;
            const long double h=std::max(source.rectangle.upper.x-source.rectangle.lower.x,source.rectangle.upper.y-source.rectangle.lower.y),half=2.5L*h;
            out<<", origin=("<<origin.x<<','<<origin.y<<"), h="<<h<<", box=[("<<origin.x-half<<','<<origin.y-half<<"),("<<origin.x+half<<','<<origin.y+half<<")]";
            // Identical fixed-box segment clipping to the estimator; returned
            // evidence is not reused by that estimator.
            const auto clip=[&](const Segment& global)->std::optional<Segment>{
                const Point a{global[0].x-origin.x,global[0].y-origin.y},b{global[1].x-origin.x,global[1].y-origin.y};
                const Point delta{b.x-a.x,b.y-a.y};long double begin=0,end=1;
                for(int axis=0;axis<2;++axis){const long double coordinate=axis?a.y:a.x,direction=axis?delta.y:delta.x;
                    if(direction==0){if(coordinate< -half||coordinate>half)return std::nullopt;}
                    else{long double lo=(-half-coordinate)/direction,hi=(half-coordinate)/direction;if(lo>hi)std::swap(lo,hi);
                        begin=std::max(begin,lo);end=std::min(end,hi);if(!(end>begin))return std::nullopt;}}
                const Segment result{{{a.x+begin*delta.x,a.y+begin*delta.y},{a.x+end*delta.x,a.y+end*delta.y}}};
                if(!(std::hypot(result[1].x-result[0].x,result[1].y-result[0].y)>0))return std::nullopt;return result;
            };
            const double negativeInfinity=-std::numeric_limits<double>::infinity(),positiveInfinity=-negativeInfinity;
            const Vector2D queryLower{std::nextafter(static_cast<double>(origin.x-half),negativeInfinity),std::nextafter(static_cast<double>(origin.y-half),negativeInfinity)};
            const Vector2D queryUpper{std::nextafter(static_cast<double>(origin.x+half),positiveInfinity),std::nextafter(static_cast<double>(origin.y+half),positiveInfinity)};
            const CartesianLeafLocator locator(mesh);auto candidates=locator.intersectingBox(queryLower,queryUpper);
            std::sort(candidates.begin(),candidates.end());candidates.erase(std::unique(candidates.begin(),candidates.end()),candidates.end());
            candidates.erase(std::remove_if(candidates.begin(),candidates.end(),[&](int ci){return l.pressureDof(ci)<0;}),candidates.end());
            // Always retain the source even when the sorted box has >96 leaves.
            const auto sourcePosition=std::find(candidates.begin(),candidates.end(),sourceNative);
            if(sourcePosition!=candidates.end())std::rotate(candidates.begin(),sourcePosition,sourcePosition+1);
            std::size_t emitted=0,eligible=0;out<<", candidate_count="<<candidates.size()<<", cells=[";
            for(int native:candidates){const auto& cell=local.at(native);const int pressure=l.pressureDof(native);
                const auto segment=segmentAt(native);const auto cut=segment?clip(*segment):std::optional<Segment>{};
                const long double dot=static_cast<long double>(cell.normal.x)*nx+static_cast<long double>(cell.normal.y)*ny;
                const char* reason=pressure<0?"missing_pressure":!mixed[pressure]?"not_measured_mixed":!cell.resolved?"unresolved":
                    !(dot>0)?"opposite_or_tangent_normal":!segment?"no_actual_segment":!cut?"no_positive_length_fixed_box_clip":"eligible_same_oriented_segment";
                if(pressure>=0&&mixed[pressure]&&cell.resolved&&dot>0&&segment&&cut)++eligible;
                if(emitted>=96)continue;if(emitted++)out<<';';
                out<<"{native="<<native<<", rectangle=[("<<cell.rectangle.lower.x<<','<<cell.rectangle.lower.y<<"),("<<cell.rectangle.upper.x<<','<<cell.rectangle.upper.y
                    <<")], Q=("<<cell.inventory[0]<<','<<cell.inventory[1]<<"), measured=("<<cell.measured[0]<<','<<cell.measured[1]<<"), normal=("<<cell.normal.x<<','<<cell.normal.y
                    <<"), offset="<<cell.offset<<", resolved="<<cell.resolved<<", dot="<<dot<<", filter="<<reason<<", actual_plane_rectangle_segment=";
                if(segment)printSegment(*segment);else out<<"none";out<<", clipped_relative_segment=";if(cut)printSegment(*cut);else out<<"none";out<<'}';
            }
            out<<"], eligible_segment_count="<<eligible<<", cell_records="<<emitted<<", cells_truncated="<<(emitted<candidates.size());
            std::size_t faceCount=0,faceEmitted=0;out<<", pure_jump_faces=[";
            const auto purePhase=[&](int native){const auto& c=local.at(native);if(!c.resolved)return -1;
                return c.inventory[0]>0&&c.inventory[1]==0&&c.measured[0]>0&&c.measured[1]==0?0:
                    c.inventory[1]>0&&c.inventory[0]==0&&c.measured[1]>0&&c.measured[0]==0?1:-1;};
            for(int row=0;row<l.domain().nativeFaceCount();++row){const auto f=l.domain().face(mesh,row);
                if(f.negativeCell<0||f.positiveCell<0)continue;const int a=purePhase(f.negativeCell),b=purePhase(f.positiveCell);if(a<0||b<0||a==b)continue;
                const auto& edge=f.horizontal?mesh.hEdges(f.edgeLabel):mesh.vEdges(f.edgeLabel);
                const auto x=mesh.verts(edge.edgeToVert[0]).coord(),y=mesh.verts(edge.edgeToVert[1]).coord();const Segment endpoints{{{x.x,x.y},{y.x,y.y}}};
                const auto cut=clip(endpoints);if(!cut)continue;++faceCount;if(faceEmitted>=32)continue;if(faceEmitted++)out<<';';
                out<<"{native_row="<<row<<", negative="<<f.negativeCell<<", positive="<<f.positiveCell<<", liquid="<<(a==0?f.negativeCell:f.positiveCell)
                    <<", gas="<<(a==1?f.negativeCell:f.positiveCell)<<", endpoints=";printSegment(endpoints);out<<", clipped_relative_segment=";printSegment(*cut);out<<'}';
            }
            out<<"], pure_jump_count="<<faceCount<<", pure_face_records="<<faceEmitted<<", pure_faces_truncated="<<(faceEmitted<faceCount)
                <<", records="<<emitted+faceEmitted<<", any_truncated="<<(emitted<candidates.size()||faceEmitted<faceCount)<<'}';
        }catch(const std::exception& error){out<<", diagnostic_error="<<error.what()<<'}';}
        return out.str();
    }
    static double coordinate(Vector2D p,int axis){return axis?p.y:p.x;}
    int pinSign(int axis) const
    {
        const int corner=static_cast<int>(pin->corner());
        return axis?(corner<2?-1:1):(corner==0||corner==3?-1:1);
    }
    int pinBin(int axis,long double tangent) const
    {
        if(!pin)return -1;const auto& c=phase.vofPublication().plicCells().at(pin->nativeCell());
        const long double h=coordinate(c.rectangle.upper,1-axis)-coordinate(c.rectangle.lower,1-axis);
        const long double distance=-pinSign(1-axis)*(tangent-coordinate(pin->physicalPoint(),1-axis));
        if(!(distance>0&&distance<2*h))return -1;return distance<h?0:1;
    }
    bool pinNormalWindowContains(int axis,const Geometry::PlanarCellCut::Rectangle& box) const
    {
        if(!pin)return false;const auto& c=phase.vofPublication().plicCells().at(pin->nativeCell());
        const double spacing=coordinate(c.rectangle.upper,axis)-coordinate(c.rectangle.lower,axis);
        const double center=coordinate((c.rectangle.lower+c.rectangle.upper)*.5,axis);
        double lower=center-4.5*spacing;const double upper=center+4.5*spacing;
        if(rz&&axisSymmetry&&axis==0)lower=std::max(0.,lower);
        return coordinate(box.lower,axis)>=lower&&coordinate(box.upper,axis)<=upper;
    }
    bool boxHitsSolid(const QuadtreeMesh& mesh,const Geometry::PlanarCellCut::Rectangle& box) const
    {
        if(!declaredLocator)return false;
        for(int ci:declaredLocator->intersectingBox(box.lower,box.upper))if(phase.layout().pressureDof(ci)<0){
            const auto& cell=mesh.cells(mesh.cellIndexToLabel[ci]);
            const auto lo=mesh.verts(cell.cellToVert[0]).coord(),hi=mesh.verts(cell.cellToVert[2]).coord();
            if(std::min(hi.x,box.upper.x)>std::max(lo.x,box.lower.x)&&
                std::min(hi.y,box.upper.y)>std::max(lo.y,box.lower.y))return true;
        }
        return false;
    }
    // A fixed 9x3 real-box plan. Only the true RZ axis has its existing even
    // reflection/truncation; missing selected fluid coverage is never gas.
    bool realColumnHeights(const QuadtreeMesh& mesh,int axis,int outward,Vector2D center,
        double dh,double spacing,const std::array<double,3>& tangentCenters,bool reflectAxis,
        std::array<long double,3>& heights,long double& origin,bool& axisImage,bool& wall,
        bool& finiteInventorySupport,std::string& reason)
    {
        const auto& local=phase.vofPublication().plicCells();std::array<Slice,27> slices;
        std::vector<HydroVofBox> queries;const double roundoff=64*std::numeric_limits<double>::epsilon();
        std::map<std::size_t,FiniteHeightMinor> finiteQueries;
        finiteInventorySupport=false;
        const auto finiteFailure=[&](const std::string& message,std::size_t q,int column,int line,
            bool rounded=false,long double lo=0,long double hi=0){
            const auto& query=queries.at(q);const auto& cell=local.at(query.nativeCell);const auto clock=phase.timeLayer();
            std::ostringstream detail;detail<<std::setprecision(std::numeric_limits<long double>::max_digits10)
                <<message<<", {"<<localPlicSourceDescription(query.nativeCell,cell)
                <<", query=[("<<query.rectangle.lower.x<<','<<query.rectangle.lower.y<<"),("<<query.rectangle.upper.x<<','<<query.rectangle.upper.y
                <<")], axis="<<axis<<", outward="<<outward<<", column="<<column<<", line="<<line
                <<", geometry_owner="<<phase.vofPublication().identity().value()<<", fullclock=(step="<<clock.step
                <<", begin="<<clock.begin<<", end="<<clock.end<<", position="<<static_cast<int>(clock.position)<<')';
            if(rounded)detail<<", endpoint_RN=("<<lo<<','<<hi<<"), endpoint_signbit=("<<std::signbit(lo)<<','<<std::signbit(hi)<<')';
            detail<<'}';reason=detail.str();return false;
        };
        for(int column=0;column<3;++column)for(int line=-4;line<=4;++line){
            Vector2D q=center;(axis?q.x:q.y)=tangentCenters[column];(axis?q.y:q.x)+=line*dh;
            const Vector2D half=axis?Vector2D{spacing/2,dh/2}:Vector2D{dh/2,spacing/2};
            Geometry::PlanarCellCut::Rectangle box{q-half,q+half};auto& slice=slices[9*column+line+4];
            slice.first=queries.size();slice.normalLower=coordinate(box.lower,axis);slice.normalUpper=coordinate(box.upper,axis);
            if(rz&&axisSymmetry&&axis==0){
                if(box.upper.x<=0){slice.last=queries.size();continue;}
                box.lower.x=std::max(0.,box.lower.x);slice.normalLower=box.lower.x;
            }
            std::vector<Geometry::PlanarCellCut::Rectangle> images;
            if(rz&&axisSymmetry&&axis==1&&box.lower.x<0&&reflectAxis){
                axisImage=true;
                if(box.upper.x<=0)images.push_back({{-box.upper.x,box.lower.y},{-box.lower.x,box.upper.y}});
                else{images.push_back({{0,box.lower.y},box.upper});images.push_back({{0,box.lower.y},{-box.lower.x,box.upper.y}});}
            }else images.push_back(box);
            long double covered=0;
            for(const auto& image:images){
                if(image.lower.x<physicalLower.x||image.lower.y<physicalLower.y||image.upper.x>physicalUpper.x||image.upper.y>physicalUpper.y){
                    wall=true;reason="fixed real column reaches an undeclared physical boundary";return false;}
                slice.capacity+=static_cast<long double>(image.upper.x-image.lower.x)*(image.upper.y-image.lower.y);
                auto candidates=declaredLocator->intersectingBox(image.lower,image.upper);
                std::sort(candidates.begin(),candidates.end());candidates.erase(std::unique(candidates.begin(),candidates.end()),candidates.end());
                for(int ci:candidates){
                    if(phase.layout().pressureDof(ci)<0){if(boxHitsSolid(mesh,image))wall=true;continue;}
                    const auto& c=local.at(ci);const auto& cb=c.rectangle;
                    Geometry::PlanarCellCut::Rectangle clip{{std::max(image.lower.x,cb.lower.x),std::max(image.lower.y,cb.lower.y)},
                        {std::min(image.upper.x,cb.upper.x),std::min(image.upper.y,cb.upper.y)}};
                    if(clip.lower.x>=clip.upper.x||clip.lower.y>=clip.upper.y)continue;
                    if(c.inventory[0]>0&&c.inventory[1]>0){
                        if(!c.resolved||c.measured[0]<=0||c.measured[1]<=0||outward*coordinate(c.normal,axis)<=0){
                            std::string finiteReason;
                            const auto finite=finiteHeightMinor(c,clip,rz,finiteReason);
                            if(finite)finiteQueries.emplace(queries.size(),*finite);
                            else{
                            // Failure-only readonly context; retain the exact
                            // original mixed-support predicate and reason prefix.
                            const auto& publication=phase.vofPublication();const auto clock=publication.timeLayer();
                            std::ostringstream detail;detail<<std::setprecision(std::numeric_limits<long double>::max_digits10)
                                <<"fixed real column has unresolved/oppositely oriented mixed PLIC support"
                                <<", {"<<localPlicSourceDescription(ci,c)<<", level="<<mesh.cells(mesh.cellIndexToLabel[ci]).level()
                                <<", rectangle=[("<<cb.lower.x<<','<<cb.lower.y<<"),("<<cb.upper.x<<','<<cb.upper.y
                                <<")], clip=[("<<clip.lower.x<<','<<clip.lower.y<<"),("<<clip.upper.x<<','<<clip.upper.y
                                <<")], query_box=[("<<box.lower.x<<','<<box.lower.y<<"),("<<box.upper.x<<','<<box.upper.y
                                <<")], actual_image=[("<<image.lower.x<<','<<image.lower.y<<"),("<<image.upper.x<<','<<image.upper.y
                                <<")], axis="<<axis<<", outward="<<outward<<", normal_component="<<coordinate(c.normal,axis)
                                <<", oriented_normal="<<outward*coordinate(c.normal,axis)<<", column="<<column<<", line="<<line
                                <<", center=("<<center.x<<','<<center.y<<"), query_center=("<<q.x<<','<<q.y
                                <<"), dh="<<dh<<", spacing="<<spacing<<", tangent_centers=("<<tangentCenters[0]<<','<<tangentCenters[1]<<','<<tangentCenters[2]
                                <<"), reflect_axis="<<reflectAxis<<", axis_image="<<axisImage<<", wall="<<wall
                                <<", geometry_owner="<<publication.identity().value()<<", fullclock=(step="<<clock.step
                                <<", begin="<<clock.begin<<", end="<<clock.end<<", position="<<static_cast<int>(clock.position)
                                <<"), failed_predicates=(unresolved="<<!c.resolved<<", measured_liquid_nonpositive="<<(c.measured[0]<=0)
                                <<", measured_gas_nonpositive="<<(c.measured[1]<=0)<<", outward_normal_nonpositive="<<(outward*coordinate(c.normal,axis)<=0)<<')';
                            if(pin){const auto point=pin->physicalPoint();const auto pinClock=pin->sourceTime();
                                detail<<", pin=(native="<<pin->nativeCell()<<", vertex="<<pin->nativeVertex()<<", point=("<<point.x<<','<<point.y
                                    <<"), corner="<<static_cast<int>(pin->corner())<<", kind="<<static_cast<int>(pin->kind())
                                    <<", source_owner="<<pin->sourceIdentity().value()<<", fullclock=(step="<<pinClock.step
                                    <<", begin="<<pinClock.begin<<", end="<<pinClock.end<<", position="<<static_cast<int>(pinClock.position)<<"))";}
                            else detail<<", pin=none";
                            detail<<", finite_height_reason={"<<finiteReason<<"}}";reason=detail.str();return false;}}
                        const double sourceSpacing=coordinate(cb.upper,1-axis)-coordinate(cb.lower,1-axis);
                        const bool partial=coordinate(clip.lower,1-axis)!=coordinate(cb.lower,1-axis)||coordinate(clip.upper,1-axis)!=coordinate(cb.upper,1-axis);
                        if(rz&&sourceSpacing>spacing*(1+roundoff)&&partial){
                            reason="RZ pin strip subdivides a coarser mixed PLIC leaf";return false;}
                    }
                    covered+=static_cast<long double>(clip.upper.x-clip.lower.x)*(clip.upper.y-clip.lower.y);queries.push_back({ci,clip});
                }
            }
            slice.last=queries.size();
            if(!(slice.capacity>0&&std::abs(covered-slice.capacity)/slice.capacity<=1e-12L)){
                reason="incomplete actual selected-fluid coverage of fixed real column";return false;}
        }
        const auto measured=rz?phase.measureVofMeridionalBoxes(mesh,queries):phase.measureVofBoxes(mesh,queries);
        for(const auto& item:finiteQueries){const auto& value=measured[item.first];const auto& minor=item.second;
            const auto major=minor.liquid?value.gas:value.liquid,empty=minor.liquid?value.liquid:value.gas;
            if(empty!=0||exactHeightNative(major)!=minor.area){
                int column=0,line=0;for(int k=0;k<27;++k)if(item.first>=slices[k].first&&item.first<slices[k].last){column=k/9;line=k%9-4;break;}
                return finiteFailure("finite-height pure-major query measure is not the exact actual rectangle area",item.first,column,line);}}
        for(auto& slice:slices){long double lc=0,gc=0;
            for(std::size_t q=slice.first;q<slice.last;++q){const auto l=measured[q].liquid-lc,ln=slice.liquid+l;lc=(ln-slice.liquid)-l;slice.liquid=ln;
                const auto g=measured[q].gas-gc,gn=slice.gas+g;gc=(gn-slice.gas)-g;slice.gas=gn;++report.heightFunctionCellQueries;}
            if(slice.capacity>0&&!(std::isfinite(slice.liquid)&&std::isfinite(slice.gas)&&slice.liquid>=0&&slice.gas>=0&&
                std::abs(slice.liquid+slice.gas-slice.capacity)/slice.capacity<=1e-12L)){
                reason="fixed real column phase capacity closure failed";return false;}
        }
        for(int column=0;column<3;++column){int first=outward>0?0:8,last=outward>0?8:0;
            while(first>=0&&first<9&&slices[9*column+first].capacity==0)first+=outward;
            while(last>=0&&last<9&&slices[9*column+last].capacity==0)last-=outward;
            if(first<0||first>=9||last<0||last>=9){reason="empty real height column";return false;}
            const auto& begin=slices[9*column+first];const auto& end=slices[9*column+last];
            if(begin.gas>localPlicHeightTopologyTolerance*begin.capacity||end.liquid>localPlicHeightTopologyTolerance*end.capacity){
                reason="fixed real column lacks opposite pure ends";return false;}
            const auto physicalOrigin=outward>0?begin.normalLower:begin.normalUpper;
            if(column==0)origin=physicalOrigin;else if(origin!=physicalOrigin){reason="real columns have different physical origins";return false;}
            long double previous=1,sum=0,compensation=0;
            for(int step=0;step<9;++step){const auto& slice=slices[9*column+(outward>0?step:8-step)];if(slice.capacity==0)continue;
                const auto fraction=slice.liquid/(slice.liquid+slice.gas);
                if(fraction>previous+localPlicHeightTopologyTolerance){reason="fixed real column has multiple phase transitions";return false;}
                previous=fraction;const auto term=slice.liquid/spacing-compensation,next=sum+term;compensation=(next-sum)-term;sum=next;
            }
            heights[column]=origin+outward*sum;
        }
        if(!finiteQueries.empty()){
            const long double inverseSpacing=1.L/static_cast<long double>(spacing);
            if(!(std::isfinite(inverseSpacing)&&inverseSpacing>0)){
                reason="unrepresentable fixed native finite-height inverse spacing";return false;}
            const auto inverse=exactHeightNative(inverseSpacing),tolerance=exactHeightNative(localPlicHeightTopologyTolerance);
            std::size_t certifiedColumns=0;std::array<bool,3> certified{};
            for(int column=0;column<3;++column){
                bool uncertain=false;std::size_t witnessQuery=0;int witnessLine=0;
                for(int line=0;line<9;++line){const auto& slice=slices[9*column+line];
                    const auto first=finiteQueries.lower_bound(slice.first);
                    if(first!=finiteQueries.end()&&first->first<slice.last){uncertain=true;witnessQuery=first->first;witnessLine=line-4;}}
                if(!uncertain)continue;
                std::array<HeightRational,9> liquid{},gas{},deltaLower{},deltaUpper{};
                HeightFinite::Dyadic exactBase(origin);HeightRational base=exactHeightNative(origin);
                HeightRational lower=0,upper=0;
                for(int line=0;line<9;++line){const auto& slice=slices[9*column+line];
                    for(std::size_t q=slice.first;q<slice.last;++q){
                        liquid[line]+=exactHeightNative(measured[q].liquid);gas[line]+=exactHeightNative(measured[q].gas);
                        exactBase.add(HeightFinite::Dyadic::product({measured[q].liquid,inverseSpacing,static_cast<long double>(outward)}));
                        base+=exactHeightNative(measured[q].liquid)*inverse*outward;
                        const auto found=finiteQueries.find(q);if(found==finiteQueries.end())continue;
                        const auto& minor=found->second;
                        if(minor.liquid){deltaLower[line]+=minor.lower;deltaUpper[line]+=minor.upper;}
                        else{deltaLower[line]-=minor.upper;deltaUpper[line]-=minor.lower;}
                    }
                    lower+=deltaLower[line];upper+=deltaUpper[line];
                }
                int first=outward>0?0:8,last=outward>0?8:0;
                while(first>=0&&first<9&&slices[9*column+first].capacity==0)first+=outward;
                while(last>=0&&last<9&&slices[9*column+last].capacity==0)last-=outward;
                if(gas[first]-deltaLower[first]>tolerance*exactHeightNative(slices[9*column+first].capacity)||
                    liquid[last]+deltaUpper[last]>tolerance*exactHeightNative(slices[9*column+last].capacity)){
                    return finiteFailure("finite-height interval lacks the original opposite pure ends",witnessQuery,column,witnessLine);}
                HeightRational previousLower=1;
                for(int step=0;step<9;++step){const int line=outward>0?step:8-step;
                    if(slices[9*column+line].capacity==0)continue;
                    const HeightRational total=liquid[line]+gas[line],lo=liquid[line]+deltaLower[line],hi=liquid[line]+deltaUpper[line];
                    if(total<=0||lo<0||hi>total||hi<lo){
                        return finiteFailure("finite-height interval is outside the actual phase-area capacity",witnessQuery,column,witnessLine);}
                    if(hi/total>previousLower+tolerance){
                        return finiteFailure("finite-height interval has multiple phase transitions under the original topology tolerance",witnessQuery,column,witnessLine);}
                    previousLower=lo/total;
                }
                long double canonicalBase=0;
                if(!roundedHeightNative(base,canonicalBase)||!sameHeightNative(canonicalBase,exactBase.rounded())){
                    return finiteFailure("finite-height exact rational/dyadic canonical sum disagree",witnessQuery,column,witnessLine);}
                HeightRational correctionLower=lower,correctionUpper=upper;
                if(outward<0){correctionLower=-upper;correctionUpper=-lower;}
                const HeightRational lo=base+inverse*correctionLower,hi=base+inverse*correctionUpper;
                long double roundedLower=0,roundedUpper=0;
                const bool lowerFinite=roundedHeightNative(lo,roundedLower),upperFinite=roundedHeightNative(hi,roundedUpper);
                if(!lowerFinite||!upperFinite)
                    return finiteFailure("finite-height interval endpoint is outside the native finite range",witnessQuery,column,witnessLine);
                if(!sameHeightNative(roundedLower,roundedUpper))
                    return finiteFailure("finite-height interval does not have a unique native nearest-even height",
                        witnessQuery,column,witnessLine,true,roundedLower,roundedUpper);
                heights[column]=roundedLower;certified[column]=true;++certifiedColumns;
            }
            finiteInventorySupport=true;
            report.finiteInventoryHeightColumns+=certifiedColumns;
            report.finiteInventoryHeightQueries+=finiteQueries.size();
            HydroCsfFiniteHeightCertificate certificate;certificate.heightAxis=axis;certificate.outward=outward;
            certificate.center=center;certificate.normalSpacing=dh;certificate.tangentSpacing=spacing;
            certificate.tangentCenters=tangentCenters;certificate.origin=origin;certificate.nativeInverseSpacing=inverseSpacing;
            certificate.physicalHeights=heights;certificate.certified=certified;certificate.queryCount=finiteQueries.size();
            report.finiteInventoryHeightCertificates.push_back(certificate);
        }
        report.heightFunctionColumns+=3;return true;
    }
    double jetCurvature(const HydroCsfHeightJet& jet) const
    {
        long double result=0;
        if(rz){const auto parts=HydroBalancedCsfLoad::axisymmetricHeightCurvatures(jet);result=parts[0]+parts[1];}
        else{const auto metric=std::hypot(1.L,jet.first);result=-jet.second/(metric*metric*metric);}
        const double value=static_cast<double>(result);
        require(std::isfinite(value)&&(value!=0||result==0),"unrepresentable qualified pin height curvature");return value;
    }
    std::optional<HydroCsfHeightJet> pinnedJet(const QuadtreeMesh& mesh,int axis,int outward,long double tangent)
    {
        const int bin=pinBin(axis,tangent);if(bin<0)return std::nullopt;
        pin->validateFor(mesh,phase.vofPublication());const auto& c=phase.vofPublication().plicCells().at(pin->nativeCell());
        const double dh=coordinate(c.rectangle.upper,axis)-coordinate(c.rectangle.lower,axis);
        const double h=coordinate(c.rectangle.upper,1-axis)-coordinate(c.rectangle.lower,1-axis),direction=-pinSign(1-axis);
        Vector2D center=(c.rectangle.lower+c.rectangle.upper)*.5;const auto point=pin->physicalPoint();
        const std::array<double,3> centers{{coordinate(point,1-axis)+direction*.5*h,
            coordinate(point,1-axis)+direction*1.5*h,coordinate(point,1-axis)+direction*2.5*h}};
        std::array<long double,3> height{};long double origin=0;bool reflected=false,wall=false,finiteInventory=false;std::string reason;
        const bool supported=realColumnHeights(mesh,axis,outward,center,dh,h,centers,false,height,origin,reflected,wall,finiteInventory,reason);
        require(supported,("pinned P3 actual support unavailable: "+reason).c_str());
        const long double hp=coordinate(point,axis),d0=height[0]-hp,d1=height[1]-hp,d2=height[2]-hp;
        const long double a=(85*d0-23*d1+4*d2)/18,b=(5*d1-10*d0-d2)/2,cc=(2*d2-7*d1+11*d0)/9;
        const long double t=direction*(tangent-coordinate(point,1-axis))/h;
        const long double value=hp+t*(a+t*(b+t*cc)),first=direction*(a+t*(2*b+3*t*cc))/h,second=(2*b+6*t*cc)/(h*h);
        HydroCsfHeightJet jet;jet.height=(value-origin)/outward;jet.first=first/outward;jet.second=second/outward;
        jet.origin=origin;jet.tangentCoordinate=tangent;jet.heightAxis=axis;jet.outward=outward;
        jet.normalSpacing=dh;jet.tangentSpacing=h;jet.pinnedCornerConstrained=true;jet.pinnedCornerBin=bin;
        jet.finiteInventoryHeightSupport=finiteInventory;
        require(std::isfinite(jet.height)&&jet.height>=0&&std::isfinite(jet.first)&&std::isfinite(jet.second),"invalid pin P3 physical jet");
        report.pinnedCornerHeightColumns+=3;return jet;
    }
    void qualifyPinCoincidentRows(const QuadtreeMesh& mesh)
    {
        if(!pin||pin->kind()!=Interface::GeometricVofPublication::PinCertificate::Kind::NativeCoincident||!pinCoincidentRows.empty())return;
        const auto& l=phase.layout();const auto& local=phase.vofPublication().plicCells();const auto& contact=local.at(pin->nativeCell());
        const int axis=contact.normal.x!=0?0:1,outward=coordinate(contact.normal,axis)>0?1:-1,tangent=1-axis;
        const double direction=-pinSign(tangent);const auto point=pin->physicalPoint();
        struct Span{double lower=0,upper=0;int row=-1;};std::vector<Span> spans;
        const auto pure=[&](int ci){if(ci<0)return -1;const auto& c=local.at(ci);
            return c.resolved&&c.inventory[0]>0&&c.inventory[1]==0?0:c.resolved&&c.inventory[1]>0&&c.inventory[0]==0?1:-1;};
        for(const auto f:l.domain().faces(mesh)){
            if((f.horizontal?1:0)!=axis||coordinate(f.center,axis)!=coordinate(point,axis)||f.area<=0)continue;
            const int minus=pure(f.negativeCell),plus=pure(f.positiveCell);
            if(minus<0||plus<0||minus==plus||(minus==0?1:-1)!=outward)continue;
            const auto& edge=f.horizontal?mesh.hEdges(f.edgeLabel):mesh.vEdges(f.edgeLabel);
            const auto x=mesh.verts(edge.edgeToVert[0]).coord(),y=mesh.verts(edge.edgeToVert[1]).coord();
            double a=direction*(coordinate(x,tangent)-coordinate(point,tangent)),b=direction*(coordinate(y,tangent)-coordinate(point,tangent));
            if(a>b)std::swap(a,b);if(a>=0&&b>a)spans.push_back({a,b,f.nativeRow});
        }
        std::sort(spans.begin(),spans.end(),[](const Span& a,const Span& b){return a.lower<b.lower;});
        double covered=0;for(const auto& span:spans){if(span.lower>covered)break;
            require(span.lower==covered,"ambiguous overlapping certified pure pin branch");pinCoincidentRows.insert(span.row);covered=span.upper;}
    }
    double pinnedCoincidentCurvature(const QuadtreeMesh& mesh,int velocity)
    {
        require(bool(pin),"coincident phase height requires explicit geometric boundary capability");
        pin->validateFor(mesh,phase.vofPublication());qualifyPinCoincidentRows(mesh);
        const auto& l=phase.layout();const auto f=l.domain().face(mesh,l.nativeFace(velocity));
        require(pinCoincidentRows.count(f.nativeRow)!=0,"pure face is not connected to the certified native coincident pin branch");
        const auto& local=phase.vofPublication().plicCells();const auto& a=local.at(f.negativeCell).rectangle;const auto& b=local.at(f.positiveCell).rectangle;
        const int axis=f.horizontal?1:0,outward=local[f.negativeCell].inventory[1]==0?1:-1,tangent=1-axis;
        const double dh=axis?std::min(a.upper.y-a.lower.y,b.upper.y-b.lower.y):std::min(a.upper.x-a.lower.x,b.upper.x-b.lower.x);
        const double dt=axis?std::min(a.upper.x-a.lower.x,b.upper.x-b.lower.x):std::min(a.upper.y-a.lower.y,b.upper.y-b.lower.y);
        const std::array<double,3> centers{{coordinate(f.center,tangent)-dt,coordinate(f.center,tangent),coordinate(f.center,tangent)+dt}};
        std::array<long double,3> heights{};long double origin=0;bool reflected=false,wall=false,finiteInventory=false;std::string reason;
        HydroCsfHeightJet jet;
        if(realColumnHeights(mesh,axis,outward,f.center,dh,dt,centers,true,heights,origin,reflected,wall,finiteInventory,reason)){
            jet.height=(heights[1]-origin)/outward;jet.first=(heights[2]-heights[0])/(2*dt*outward);
            jet.second=(heights[2]-2*heights[1]+heights[0])/(dt*dt*outward);jet.origin=origin;
            jet.heightAxis=axis;jet.outward=outward;jet.tangentCoordinate=coordinate(f.center,tangent);jet.axisReflected=reflected;
            jet.finiteInventoryHeightSupport=finiteInventory;
        }else{
            require(wall,("certified pure pin face lacks ordinary single-branch/pure-end support: "+reason).c_str());
            const auto pinned=pinnedJet(mesh,axis,outward,coordinate(f.center,tangent));
            require(bool(pinned),"certified pure pin face outside the first/second blocked tangent bin");jet=*pinned;++report.pinnedCornerHeightFaces;
        }
        ++report.pureJumpHeightFaces;return jetCurvature(jet);
    }
    double coincidentFaceCurvature(const QuadtreeMesh& mesh,int velocity)
    {
        if(pin){qualifyPinCoincidentRows(mesh);
            if(pinCoincidentRows.count(phase.layout().nativeFace(velocity)))return pinnedCoincidentCurvature(mesh,velocity);}
        const bool interiorNative=nativeSupportChart&&!south&&!pin;
        if(!south&&!interiorNative)return pinnedCoincidentCurvature(mesh,velocity);
        require(south!=nullptr||interiorNative,"coincident phase height requires explicit geometric boundary capability");
        const auto& layout=phase.layout();const auto f=layout.domain().face(mesh,layout.nativeFace(velocity));const auto& local=phase.vofPublication().plicCells();
        require(f.negativeCell>=0&&f.positiveCell>=0&&f.area>0,"coincident phase height requires a real positive-area shared face");
        const auto pure=[&](int ci){const auto& c=local.at(ci);
            // Force color and this interface position use actual measured
            // geometry. An independent unrepresentable tail is still kept in
            // inventory; it does not remove a full geometric phase face.
            if(interiorNative)return c.measured[0]>0&&c.measured[1]==0?0:c.measured[1]>0&&c.measured[0]==0?1:-1;
            if(!c.resolved)return -1;
            return c.inventory[0]>0&&c.inventory[1]==0&&c.measured[0]>0&&c.measured[1]==0?0:c.inventory[1]>0&&c.inventory[0]==0&&c.measured[1]>0&&c.measured[0]==0?1:-1;};
        const int minus=pure(f.negativeCell),plus=pure(f.positiveCell);require(minus>=0&&plus>=0&&minus!=plus,"active CSF face is not an exact pure-L/G shared interface");
        const auto& a=local.at(f.negativeCell).rectangle;const auto& b=local.at(f.positiveCell).rectangle;
        Frame frame{f.center,std::min(a.upper.x-a.lower.x,b.upper.x-b.lower.x),std::min(a.upper.y-a.lower.y,b.upper.y-b.lower.y),false};
        const int axis=f.horizontal?1:0,outward=minus==0?1:-1,tangent=1-axis;
        if(interiorNative){
            if(!declaredLocator)declaredLocator=std::make_unique<CartesianLeafLocator>(mesh);
            (tangent?frame.hy:frame.hx)=f.area;
            // Align normal slices with the adjacent liquid leaf, rather than
            // averaging different pure branches into half-cell .5 plateaus.
            (axis?frame.center.y:frame.center.x)-=outward*(axis?frame.hy:frame.hx)/2;
        }
        const double dh=axis?frame.hy:frame.hx,spacing=tangent?frame.hy:frame.hx;
        Direction direction;std::vector<HydroVofBox> queries;
        for(int column=-1;column<=1&&direction.reason.empty();++column)for(int line=-4;line<=4&&direction.reason.empty();++line){
            Vector2D center=frame.center;(axis?center.y:center.x)+=line*dh;(tangent?center.y:center.x)+=column*spacing;
            const Vector2D half=axis?Vector2D{spacing/2,dh/2}:Vector2D{dh/2,spacing/2};
            appendDeclaredSlice(mesh,axis,outward,{center-half,center+half},direction.slices[9*(column+1)+line+4],direction.reason,queries,false);}
        require(direction.reason.empty(),("coincident phase height support unavailable: "+direction.reason).c_str());
        const auto measured=rz?phase.measureVofMeridionalBoxes(mesh,queries):phase.measureVofBoxes(mesh,queries);
        for(auto& slice:direction.slices){if(rz&&slice.capacity==0)continue;long double lc=0,gc=0;
            for(std::size_t q=slice.first;q<slice.last;++q){const auto l=measured[q].liquid-lc,ln=slice.liquid+l;lc=(ln-slice.liquid)-l;slice.liquid=ln;
                const auto g=measured[q].gas-gc,gn=slice.gas+g;gc=(gn-slice.gas)-g;slice.gas=gn;++report.heightFunctionCellQueries;}
            require(std::isfinite(slice.liquid)&&std::isfinite(slice.gas)&&slice.liquid>=0&&slice.gas>=0&&std::abs(slice.liquid+slice.gas-slice.capacity)/slice.capacity<=1e-12L,"coincident phase height-box closure failed");}
        bool needLower=false,needUpper=false;std::string reason;
        require(orderedSupport(direction,outward,needLower,needUpper,reason)&&!needLower&&!needUpper,("coincident phase height lacks original single-branch/pure-end support: "+reason).c_str());
        if(interiorNative){long double height=0,correction=0;
            for(int step=0;step<9;++step){const auto& slice=direction.slices[9+(outward>0?step:8-step)];
                const auto term=slice.liquid/(slice.liquid+slice.gas)-correction,next=height+term;correction=(next-height)-term;height=next;}
            require(std::abs(height-5)<=9*localPlicHeightTopologyTolerance+64*std::numeric_limits<double>::epsilon(),
                "native coincident height graph does not locate its actual shared phase face");}
        direction.extendedReady=true;double result=0;
        require(estimateWindow(-1,direction,frame,axis,outward,result,reason,false),("coincident phase height curvature unavailable: "+reason).c_str());
        report.heightFunctionColumns+=3;++report.pureJumpHeightFaces;return result;
    }
    double curvature(int p,const QuadtreeMesh& mesh)
    {
        // Explicit native-support resolution model, selected before fitting.
        // Effective capillary closure is zero for an unresolved stock feature;
        // its physical curvature remains unavailable and no isolation is claimed.
        if(subgrid[phase.layout().nativeCell(p)]){
            methods[p]=HydroCsfCellCurvatureMethod::Unavailable;return 0;
        }
        const int first=std::abs(normal[p].x)>=std::abs(normal[p].y)?0:1;int selected=first;double result=0;std::array<std::string,2> reasons;
        bool extended=false,usedRefined=false,usedFit=false;
        if(!estimate(p,first,result,reasons[first])){++report.heightFunctionAlternateDirections;selected=1-first;
            if(!estimate(p,selected,result,reasons[selected])){
                const auto& source=phase.vofPublication().plicCells().at(phase.layout().nativeCell(p));
                const auto center=(source.rectangle.lower+source.rectangle.upper)*.5;
                if(pin&&cellDirections(p)[first].wallBlocked&&pinBin(first,coordinate(center,1-first))>=0&&
                    pinNormalWindowContains(first,source.rectangle)){
                    const int outward=coordinate(source.normal,first)>0?1:-1;
                    const auto jet=pinnedJet(mesh,first,outward,coordinate(center,1-first));
                    require(bool(jet),"missing qualified first/second-bin pin jet");jets[p]=*jet;result=jetCurvature(*jet);
                    methods[p]=HydroCsfCellCurvatureMethod::HeightFunction;++report.pinnedCornerHeightCells;
                    ++report.heightFunctionCells;(first==0?report.heightFunctionXCells:report.heightFunctionYCells)++;return result;
                }
                if(cellDirections(p)[0].wallBlocked||cellDirections(p)[1].wallBlocked){
                    const int native=phase.layout().nativeCell(p);const auto clock=phase.vofPublication().timeLayer();
                    std::ostringstream detail;detail<<std::setprecision(std::numeric_limits<long double>::max_digits10)
                        <<"wall-blocked LocalPlic height support requires the actual nearby typed pin chart; "
                        <<localPlicSourceDescription(native,source)<<", pressure="<<p
                        <<", level="<<mesh.cells(mesh.cellIndexToLabel[native]).level()
                        <<", rectangle=[("<<source.rectangle.lower.x<<','<<source.rectangle.lower.y<<"),("<<source.rectangle.upper.x<<','<<source.rectangle.upper.y
                        <<")], first_axis="<<first<<", clock=("<<clock.step<<','<<clock.begin<<','<<clock.end<<','<<static_cast<int>(clock.position)
                        <<"), publication="<<phase.vofPublication().identity().value()<<", pin="<<bool(pin);
                    if(pin){const auto point=pin->physicalPoint();const auto& contact=phase.vofPublication().plicCells().at(pin->nativeCell());
                        const auto size=contact.rectangle.upper-contact.rectangle.lower;
                        detail<<", pin_native="<<pin->nativeCell()<<", pin_level="<<mesh.cells(mesh.cellIndexToLabel[pin->nativeCell()]).level()
                            <<", pin_kind="<<static_cast<int>(pin->kind())<<", pin_corner="<<static_cast<int>(pin->corner())
                            <<", pin_point=("<<point.x<<','<<point.y<<"), pin_rectangle=[("<<contact.rectangle.lower.x<<','<<contact.rectangle.lower.y
                            <<"),("<<contact.rectangle.upper.x<<','<<contact.rectangle.upper.y<<")], center_minus_pin=("<<center.x-point.x<<','<<center.y-point.y
                            <<"), center_distance_in_contact_cells=("<<(center.x-point.x)/size.x<<','<<(center.y-point.y)/size.y<<')';}
                    for(int axis=0;axis<2;++axis){const auto& frame=primaryFrames[p][axis];const int tangent=1-axis;
                        detail<<"; axis="<<axis<<", wall_blocked="<<cellDirections(p)[axis].wallBlocked
                            <<", pin_bin="<<pinBin(axis,coordinate(center,tangent))<<", pin_normal_window="<<pinNormalWindowContains(axis,source.rectangle)
                            <<", estimate_reason={"<<reasons[axis]<<"}, initial_reason={"<<cellDirections(p)[axis].reason
                            <<"}, frame=("<<frame.center.x<<','<<frame.center.y<<','<<frame.hx<<','<<frame.hy<<','<<frame.coarseAligned<<')';
                        if(!cellDirections(p)[axis].wallBlocked)continue;
                        bool found=false;
                        // Replay only fixed query geometry after failure; no phase/source mutation.
                        for(int column=-1;column<=1&&!found;++column)for(int line=-4;line<=4&&!found;++line){
                            const auto query=fixedQueryBox(mesh,frame,axis,line,column);if(!boxHitsSolid(mesh,query))continue;
                            found=true;detail<<", first_wall_query={column="<<column<<", line="<<line<<", box=[("<<query.lower.x<<','<<query.lower.y
                                <<"),("<<query.upper.x<<','<<query.upper.y<<")], solid_leaves=[";bool emitted=false;
                            for(int ci:declaredLocator->intersectingBox(query.lower,query.upper))if(phase.layout().pressureDof(ci)<0){
                                const auto& cell=mesh.cells(mesh.cellIndexToLabel[ci]);const auto lo=mesh.verts(cell.cellToVert[0]).coord(),hi=mesh.verts(cell.cellToVert[2]).coord();
                                if(!(std::min(hi.x,query.upper.x)>std::max(lo.x,query.lower.x)&&std::min(hi.y,query.upper.y)>std::max(lo.y,query.lower.y)))continue;
                                if(emitted)detail<<';';emitted=true;detail<<"native="<<ci<<",level="<<cell.level()<<",rectangle=[("<<lo.x<<','<<lo.y<<"),("<<hi.x<<','<<hi.y<<")]";
                            }
                            detail<<"]}";
                        }
                    }
                    require(false,detail.str().c_str());
                }
                std::array<std::string,2> boundedReasons,refinedReasons;selected=first;
                if(estimateExtended(p,selected,result,boundedReasons[selected]))extended=true;
                else{selected=1-first;if(estimateExtended(p,selected,result,boundedReasons[selected]))extended=true;}
                if(!extended){selected=first;
                    if(estimateRefined(p,selected,result,refinedReasons[selected]))usedRefined=true;
                    else{selected=1-first;if(estimateRefined(p,selected,result,refinedReasons[selected]))usedRefined=true;}}
                if(!extended&&!usedRefined){const auto fit=fits.find(p);
                    if(fit!=fits.end()&&fit->second.ready){result=fit->second.curvature;usedFit=true;}
                    else{const int native=phase.layout().nativeCell(p);
                        // An insufficient fit alone does not authorize zero. The
                        // Interface owner must prove a closed isolated component
                        // whose whole geometry cannot supply Delta-separated
                        // positions. This never overwrites the source geometry.
                        const auto fragment=phase.vofPublication().isolatedLocalPlicFragment(mesh,native);
                        if(fragment){fragment->validateFor(mesh,phase.vofPublication());
                            methods[p]=HydroCsfCellCurvatureMethod::UnderresolvedIsolatedFragmentZero;
                            ++report.unresolvedIsolatedFragmentCells;return 0;}
                        // A separate typed proof can cover a closed two-cell
                        // component whose two-bin packing excludes three
                        // Delta-separated positions. No singleton diameter
                        // assertion or physical curvature is inferred here.
                        const auto deficit=phase.vofPublication().twoCellLocalPlicPositionDeficit(mesh,native);
                        if(deficit){deficit->validateFor(mesh,phase.vofPublication());
                            methods[p]=HydroCsfCellCurvatureMethod::UnderresolvedIsolatedFragmentZero;
                            ++report.unresolvedIsolatedFragmentCells;return 0;}
                        const auto star=phase.vofPublication().vertexStarLocalPlicFragment(mesh,native);
                        if(star){star->validateFor(mesh,phase.vofPublication());
                            methods[p]=HydroCsfCellCurvatureMethod::UnderresolvedIsolatedFragmentZero;
                            ++report.unresolvedIsolatedFragmentCells;++report.unresolvedVertexStarFragmentCells;return 0;}
                        for(int axis=0;axis<2;++axis){if(!boundedReasons[axis].empty())reasons[axis]+=", bounded="+boundedReasons[axis];
                            if(!refinedReasons[axis].empty())reasons[axis]+=", refined="+refinedReasons[axis];}
                        const std::string fitReason=fit==fits.end()?"local PLIC quadratic fit was not prepared":fit->second.reason;
                        throw std::invalid_argument("[HydroBalancedCsfLoad] local PLIC height support unavailable, "+
                            localPlicSourceDescription(native,phase.vofPublication().plicCells()[native])+", x="+reasons[0]+", y="+reasons[1]+", fit="+fitReason+
                            ", diagnostic="+failureSupportDiagnostic(mesh,p));}}}}
        if(usedFit){const auto& fit=fits.at(p);methods[p]=fit.mixedHeightPoints?HydroCsfCellCurvatureMethod::MixedHeightFit:
                fit.twoSegmentCircle?HydroCsfCellCurvatureMethod::TwoSegmentCircleFit:HydroCsfCellCurvatureMethod::SegmentFit;
            ++report.heightFunctionFitCells;report.heightFunctionFitSegments+=fit.segments;return result;}
        if(extended)++report.heightFunctionExtendedCells;
        if(usedRefined)++report.heightFunctionRefinedCells;
        if(rz&&jets[p]&&jets[p]->coarseAligned)++report.heightFunctionCoarseAlignedCells;
        methods[p]=HydroCsfCellCurvatureMethod::HeightFunction;++report.heightFunctionCells;(selected==0?report.heightFunctionXCells:report.heightFunctionYCells)++;return result;
    }
};

struct HeightFunctionGrid
{
    const HydroPhaseGeometry& phase;
    const LevelSetHermiteGeometry::Snapshot& scalar;
    UniformLeafGrid grid;
    HydroCsfCurvatureReport& report;
    std::vector<int> pressureAt,graphs;
    std::vector<std::array<int,2>> position,partial;
    std::vector<std::array<Vector2D,2>> bounds;
    std::vector<Vector2D> normal;
    std::vector<long double> fractions;
    HeightFunctionGrid(const QuadtreeMesh& mesh,const HydroPhaseGeometry& source,HydroCsfCurvatureReport& r)
        :phase(source),scalar(source.publication().state().geometry().snapshot()),
         grid(mesh,CartesianGridScope::UniformLeaves),report(r)
    {
        const auto& l=phase.layout();const int np=l.pressureCount();
        require(grid.nx()>=3&&grid.ny()>=3,"height-function requires at least three uniform XY cells in each direction");
        require(static_cast<std::size_t>(np)==grid.cells().size(),"height-function requires complete pressure-cell support");
        pressureAt.resize(grid.cells().size());position.resize(np);partial.assign(np,{{0,0}});
        bounds.resize(np);fractions.resize(np);normal.resize(np);graphs.assign(np,0);
        std::vector<unsigned char> seen(np,0);
        for(int j=0;j<grid.ny();++j)for(int i=0;i<grid.nx();++i){
            const int native=mesh.cells(grid.cells()[j*grid.nx()+i]).index();const int p=l.pressureDof(native);
            require(p>=0&&!seen[p],"height-function requires complete unique native pressure-cell support");seen[p]=1;
            require(scalar.nativeCellActive(native),"height-function requires complete active SingleQuartic support");
            pressureAt[j*grid.nx()+i]=p;position[p]={{i,j}};bounds[p]=scalar.nativeCellBounds(native);
            const long double capacity=phase.metrics().cellVolumes()[p],liquid=phase.cellVolumes()[p].liquid;
            require(std::isfinite(capacity)&&capacity>0&&std::isfinite(liquid)&&liquid>=0&&liquid<=capacity,
                "height-function invalid actual phase volume/capacity");
            fractions[p]=liquid/capacity;
        }
        const auto& publication=phase.publication();
        for(std::size_t c=0;c<publication.topology().curves.size();++c)
            for(std::size_t g=0;g<publication.topology().curves[c].graphs().size();++g){
                const int p=l.pressureDof(publication.topology().curves[c].graphs()[g].nativeCellIndex);
                require(p>=0,"height-function interface outside complete pressure domain");
                if(graphs[p]++==0){const auto& cached=publication.quadrature(c,g);
                    require(!cached.empty(),"height-function lacks certified publication normal");
                    const auto& sample=cached.front().sample;
                    const auto sampleClock=sample.timeLayer,clock=phase.timeLayer();
                    require(sample.geometryIdentity==phase.geometryIdentity()&&sample.frame.nativeCellIndex==l.nativeCell(p)
                        &&sampleClock.step==clock.step&&sampleClock.begin==clock.begin&&sampleClock.end==clock.end&&sampleClock.position==clock.position,
                        "height-function normal has different retained geometry/full clock owner");
                    normal[p]=sample.frame.normal;
                    require(std::isfinite(normal[p].x)&&std::isfinite(normal[p].y)&&std::hypot(normal[p].x,normal[p].y)>0,
                        "height-function invalid retained interface normal");}
            }
    }
    int at(int i,int j) const
    {return i<0||i>=grid.nx()||j<0||j>=grid.ny()?-1:pressureAt[j*grid.nx()+i];}
    long double value(int p)
    {++report.heightFunctionCellQueries;return fractions[p];}
    bool containsInterface(int p,long double fraction) const
    {return graphs[p]>0||(fraction>0&&fraction<1);}
    bool regular(int p,int axis,int outward,std::string& reason)
    {
        if(graphs[p]==0){reason="partial-volume cell lacks retained interface branch";return false;}
        if(partial[p][axis]==0){++report.heightFunctionPartialCertificates;
            partial[p][axis]=scalar.nativeZeroSetPartialResolved(phase.layout().nativeCell(p),axis)?1:-1;}
        if(partial[p][axis]<0){reason="complete native zero set is not single-valued in height direction";return false;}
        const double component=axis==0?normal[p].x:normal[p].y;
        if(component==0||(component>0?1:-1)!=outward){reason="column contains an oppositely oriented interface branch";return false;}
        return true;
    }
    bool column(int i,int j,int axis,int outward,long double& height,std::string& reason)
    {
        ++report.heightFunctionColumns;
        const int initial=at(i,j);if(initial<0){reason="height transverse column leaves physical domain";return false;}
        const auto endpoint=[&](int direction,long double expected,int& endpointCell){
            int x=i,y=j,p=initial;long double fraction=value(p);
            bool found=(expected==0?fraction<1:fraction>0)||graphs[p]>0;
            while(!found||containsInterface(p,fraction)){
                (axis==0?x:y)+=direction;p=at(x,y);
                if(p<0){reason="height column reaches physical boundary before a pure-phase end";return false;}
                fraction=value(p);if(containsInterface(p,fraction))found=true;
            }
            if(fraction!=expected){reason="height column closes in the wrong pure phase";return false;}
            endpointCell=p;return true;
        };
        int top=-1,bottom=-1;
        if(!endpoint(outward,0,top)||!endpoint(-outward,1,bottom))return false;
        const int start=position[bottom][axis],finish=position[top][axis];
        if(outward*(finish-start)<=0){reason="height column has no complete liquid-to-gas transition";return false;}
        long double depth=0,compensation=0,previous=1;bool branch=false;
        for(int k=start;;k+=outward){const int p=axis==0?at(k,j):at(i,k);const long double fraction=value(p);
            if(fraction>previous){reason="height column has multiple liquid-to-gas transitions";return false;}
            previous=fraction;
            if(containsInterface(p,fraction)){branch=true;if(!regular(p,axis,outward,reason))return false;}
            const long double width=axis==0?static_cast<long double>(bounds[p][1].x)-bounds[p][0].x:
                static_cast<long double>(bounds[p][1].y)-bounds[p][0].y;
            const long double term=fraction*width-compensation,next=depth+term;
            compensation=(next-depth)-term;depth=next;
            if(k==finish)break;
        }
        if(!branch){reason="height column has no identified single interface transition";return false;}
        const double coordinate=axis==0?(outward>0?bounds[bottom][0].x:bounds[bottom][1].x):
            (outward>0?bounds[bottom][0].y:bounds[bottom][1].y);
        height=static_cast<long double>(outward)*coordinate+depth;
        if(!std::isfinite(height)){reason="height column is not representable";return false;}
        return true;
    }
    bool estimate(int p,int axis,double& curvature,std::string& reason)
    {
        const double component=axis==0?normal[p].x:normal[p].y;
        if(component==0){reason="height direction has zero retained normal component";return false;}
        const int outward=component>0?1:-1;
        if(!regular(p,axis,outward,reason))return false;
        std::array<long double,3> heights{};const int i=position[p][0],j=position[p][1];
        for(int transverse=-1;transverse<=1;++transverse)
            if(!column(i+(axis==1?transverse:0),j+(axis==0?transverse:0),axis,outward,heights[transverse+1],reason))return false;
        const long double spacing=axis==0?grid.spacing().y:grid.spacing().x;
        const long double left=heights[1]-heights[0],right=heights[2]-heights[1];
        const long double first=(left+right)/(2*spacing),second=(right-left)/(spacing*spacing);
        const long double metric=std::hypot(1.L,first),k=-second/(metric*metric*metric);
        curvature=static_cast<double>(k);
        if(!std::isfinite(curvature)||(curvature==0&&k!=0)){reason="height curvature is not representable";return false;}
        return true;
    }
    double curvature(int p)
    {
        const int first=std::abs(normal[p].x)>=std::abs(normal[p].y)?0:1;
        std::array<std::string,2> reasons;double result=0;int selected=first;
        if(!estimate(p,first,result,reasons[first])){
            ++report.heightFunctionAlternateDirections;selected=1-first;
            if(!estimate(p,selected,result,reasons[selected])){
                throw std::invalid_argument("[HydroBalancedCsfLoad] height-function support unavailable, native="+
                    std::to_string(phase.layout().nativeCell(p))+", x="+reasons[0]+", y="+reasons[1]);}
        }
        ++report.heightFunctionCells;
        (selected==0?report.heightFunctionXCells:report.heightFunctionYCells)++;
        return result;
    }
};

}
struct HydroVofSouthEven90Boundary::Data {
    Data(Interface::GeometricVofPublication value):publication(std::move(value)){}
    Interface::GeometricVofPublication publication;std::vector<int> rows;std::vector<std::array<double,2>> intervals;double south=0;
};
HydroVofSouthEven90Boundary::HydroVofSouthEven90Boundary(std::shared_ptr<const Data> data):data_(std::move(data)){}
HydroVofSouthEven90Boundary HydroVofSouthEven90Boundary::prepare(const QuadtreeMesh& mesh,const Interface::GeometricVofPublication& pub,const std::vector<int>& rows) {
    pub.validateFor(mesh,pub.identity(),pub.timeLayer());require(pub.geometryKind()==Interface::GeometricVofGeometryKind::LocalPlic,"S Even90 requires an actual LocalPlic source");
    require(pub.domain().cellRowByIndex().size()==static_cast<std::size_t>(mesh.cellIndexNum()),"S Even90 requires the complete native-indexed selected map");require(!rows.empty(),"S Even90 requires declared actual south native faces");
    auto data=std::make_shared<Data>(pub);data->south=mesh.minCoord.y;data->rows=rows;std::sort(data->rows.begin(),data->rows.end());
    require(std::adjacent_find(data->rows.begin(),data->rows.end())==data->rows.end(),"S Even90 has duplicate actual face rows");std::vector<std::array<double,2>> pieces;
    for(int row:data->rows){const auto f=pub.domain().face(mesh,row);
        require(f.horizontal&&f.negativeCell<0&&f.positiveCell>=0&&f.outwardSign==-1&&f.center.y==data->south&&f.area>0,"S Even90 row is not a whole positive-area actual south face");
        const auto& edge=mesh.hEdges(f.edgeLabel);const auto a=mesh.verts(edge.edgeToVert[0]).coord(),b=mesh.verts(edge.edgeToVert[1]).coord();
        require(a.y==data->south&&b.y==data->south&&a.x!=b.x,"S Even90 actual face endpoints disagree with south plane");pieces.push_back({std::min(a.x,b.x),std::max(a.x,b.x)});}
    std::sort(pieces.begin(),pieces.end());for(const auto& p:pieces){if(!data->intervals.empty()&&p[0]==data->intervals.back()[1])data->intervals.back()[1]=p[1];
        else{require(data->intervals.empty()||p[0]>data->intervals.back()[1],"S Even90 actual intervals overlap");data->intervals.push_back(p);}}
    return HydroVofSouthEven90Boundary(std::move(data));
}
const Interface::GeometricVofPublication& HydroVofSouthEven90Boundary::source() const {require(bool(data_),"consumed S Even90 view");return data_->publication;}
const std::vector<int>& HydroVofSouthEven90Boundary::nativeRows() const {require(bool(data_),"consumed S Even90 view");return data_->rows;}
const std::vector<std::array<double,2>>& HydroVofSouthEven90Boundary::intervals() const {require(bool(data_),"consumed S Even90 view");return data_->intervals;}
double HydroVofSouthEven90Boundary::coordinate() const {require(bool(data_),"consumed S Even90 view");return data_->south;}
bool HydroVofSouthEven90Boundary::sameIdentity(const HydroVofSouthEven90Boundary& other) const noexcept {return data_&&data_==other.data_;}
void HydroVofSouthEven90Boundary::validateFor(const QuadtreeMesh& mesh,const Interface::GeometricVofPublication& pub) const {
    require(bool(data_),"consumed S Even90 view");pub.validateFor(mesh,pub.identity(),pub.timeLayer());const auto a=pub.timeLayer(),b=data_->publication.timeLayer();
    require(pub.identity()==data_->publication.identity()&&pub.domain().meshIdentity()==data_->publication.domain().meshIdentity()&&a.step==b.step&&a.begin==b.begin&&a.end==b.end&&a.position==b.position,"S Even90 has a foreign publication/mesh/full clock");
}
struct HydroBalancedCsfLoad::Data
{
    Data(HydroPressureLoad p,HydroRegularizedColor c,double s,double k)
        :pressure(std::move(p)),phase(c.phaseGeometry()),color(std::move(c)),sigma(s),curvature(k){}
    Data(HydroPressureLoad p,HydroPhaseGeometry g,double s)
        :pressure(std::move(p)),phase(std::move(g)),sigma(s),curvature(0),kind(HydroCsfDiscretization::GeometricPhaseFraction){}
    HydroPressureLoad pressure;
    HydroPhaseGeometry phase;
    std::optional<HydroRegularizedColor> color;
    std::optional<HydroVofSouthEven90Boundary> southEven90;
    double sigma,curvature;
    HydroCsfDiscretization kind=HydroCsfDiscretization::RegularizedColor;
    std::vector<double> gradient,rhs;
    bool capillary=false;
    HydroCsfCurvatureConfig config;
    std::vector<double> averages,heightCurvature,heightEffectiveCurvature;
    std::vector<HydroCsfCellCurvatureMethod> heightCurvatureMethods;
    std::vector<std::optional<HydroCsfHeightJet>> heightJets;
    std::vector<double> incrementReferenceLoad,incrementResidualLoad,incrementDeltaMean;
    std::shared_ptr<ProjectedRateBasis> incrementRates;
    std::shared_ptr<std::mutex> incrementRateMutex;
    std::vector<HydroCsfFaceSupportReport> support;
    std::vector<std::vector<HydroCsfQuadraturePanel>> panels;
    HydroCsfCurvatureReport report;
    std::optional<HydroAdvectorReconstruction> source;
    std::optional<HydroAdvectorBoundary> boundary;
    std::optional<HydroSpatialAdvector> projectedMotion;
    std::optional<Interface::HermiteGalsAmrClosureSensitivity> mixed;
    std::vector<double> referenceCurvature,referenceStress,referenceLoad,residualLoad,tangentialLoad;
    std::optional<Interface::HermiteGalsSurfaceVectorView> externalSurface;
};
HydroBalancedCsfLoad::HydroBalancedCsfLoad(std::shared_ptr<const Data> d):data_(std::move(d)){}
HydroBalancedCsfLoad HydroBalancedCsfLoad::constantCurvature(const QuadtreeMesh& mesh,
    HydroPressureLoad pressure,HydroRegularizedColor color,double sigma,double curvature)
{
    require(std::isfinite(sigma)&&sigma>=0,"invalid surface tension");
    require(std::isfinite(curvature),"nonfinite constant curvature");
    require(pressure.layout().sameIdentity(color.phaseGeometry().layout()),"different pressure/color layout");
    require(&pressure.divergence().matrix()==&color.phaseGeometry().metrics().divergence().matrix(),
        "different authoritative divergence owner");
    color.validateFor(mesh,color.geometryIdentity(),color.timeLayer());
    const double coefficient=sigma*curvature;
    require(std::isfinite(coefficient)&&(coefficient!=0||sigma==0||curvature==0),
        "unrepresentable constant CSF coefficient");
    auto d=std::make_shared<Data>(std::move(pressure),std::move(color),sigma,curvature);
    const int count=d->pressure.layout().velocityCount();
    // Reuse the only authoritative gradient. No new stencil, area factor,
    // density factor, sharp alpha substitution, or virtual coarse-face split.
    d->gradient=d->pressure.apply(mesh,HydroNativeValues(d->pressure.layout(),std::vector<double>(count,0),d->color->colors()));
    d->rhs.assign(count,0);
    if(sigma!=0&&curvature!=0)for(int u=0;u<count;++u){
        d->rhs[u]=coefficient*d->gradient[u];
        require(std::isfinite(d->rhs[u]),"nonfinite constant CSF momentum RHS");}
    return HydroBalancedCsfLoad(std::move(d));
}
double HydroBalancedCsfLoad::positiveCurvatureAverage(const std::vector<double>& weights,const std::vector<double>& values)
{
    require(!weights.empty()&&weights.size()==values.size(),"invalid curvature average dimensions");
    double ws=0,ks=0;std::size_t anchor=weights.size();bool positive=false,negative=false,constant=true;
    for(std::size_t i=0;i<weights.size();++i){require(std::isfinite(weights[i])&&weights[i]>=0&&std::isfinite(values[i]),"invalid positive curvature average input");
        if(weights[i]>0){ws=std::max(ws,weights[i]);ks=std::max(ks,std::abs(values[i]));if(anchor==weights.size())anchor=i;
            constant=constant&&values[i]==values[anchor];
            positive=positive||values[i]>0;negative=negative||values[i]<0;}}
    require(ws>0,"no positive curvature support");if(constant)return values[anchor];
    const long double base=static_cast<long double>(values[anchor])/ks;long double total=0,offset=0;
    for(std::size_t i=0;i<weights.size();++i)if(weights[i]>0){const long double w=static_cast<long double>(weights[i])/ws,
            k=static_cast<long double>(values[i])/ks;
        require(w>0&&(k!=0||values[i]==0),"unrepresentable normalized curvature average input");
        const long double delta=k-base,term=w*delta;
        require(term!=0||delta==0,"unrepresentable normalized curvature average moment");total+=w;offset+=term;}
    const double mean=static_cast<double>((base+offset/total)*ks);
    require(std::isfinite(mean)&&(mean!=0||positive==negative),"unrepresentable positive curvature average");return mean;
}
bool HydroBalancedCsfLoad::selectedSegmentVisible(Vector2D a,Vector2D b,const std::vector<std::array<Vector2D,2>>& boxes)
{
    Exact::validateEnvironment();const RayLine ray(a,b);RayIntervals spans;
    for(const auto& box:boxes){require(std::isfinite(box[0].x)&&std::isfinite(box[0].y)&&std::isfinite(box[1].x)&&std::isfinite(box[1].y)&&
            box[0].x<=box[1].x&&box[0].y<=box[1].y,"invalid selected visibility box");
        RayParameter low(0),high(1);if(ray.clip({box[0],box[1]},low,high))spans.emplace_back(std::move(low),std::move(high));}
    return covered(spans);
}
HydroBalancedCsfLoad HydroBalancedCsfLoad::capillary(const QuadtreeMesh& mesh,HydroPressureLoad pressure,
    HydroRegularizedColor color,double sigma,HydroCsfCurvatureConfig cfg)
{
    require(std::isfinite(sigma)&&sigma>=0,"invalid surface tension");
    validateCurvatureConfig(cfg);
    require(pressure.layout().sameIdentity(color.phaseGeometry().layout()),"different pressure/color layout");
    require(&pressure.divergence().matrix()==&color.phaseGeometry().metrics().divergence().matrix(),"different authoritative divergence owner");
    color.validateFor(mesh,color.geometryIdentity(),color.timeLayer());
    auto d=std::make_shared<Data>(std::move(pressure),std::move(color),sigma,0);d->capillary=true;d->config=cfg;
    const auto& l=d->pressure.layout();const int nu=l.velocityCount();
    d->gradient=d->pressure.apply(mesh,HydroNativeValues(l,std::vector<double>(nu,0),d->color->colors()));
    d->rhs.assign(nu,0);d->averages.assign(nu,0);d->support.resize(nu);d->panels.resize(nu);
    if(sigma==0){d->report.projectionSkippedSigmaZero=true;for(auto& r:d->support)r.status=HydroCsfSupportStatus::SigmaZero;
        return HydroBalancedCsfLoad(std::move(d));}
    ProjectionIndex index(mesh,*d->color,d->report);std::vector<double> weights,curvatures;
    const auto& cells=d->color->cellReports();
    for(int u=0;u<nu;++u){auto& r=d->support[u];const auto& face=index.faces[u];
        const int a=face.negativeCell<0?-1:l.pressureDof(face.negativeCell),b=face.positiveCell<0?-1:l.pressureDof(face.positiveCell);
        const int p=a>=0?a:b;require(p>=0,"CSF face without native pressure cell");r.pressureComponent=index.component[p];
        r.center=face.center;r.boundary=a<0||b<0;
        const double diameter=std::max(a<0?0:index.diagonal[a],b<0?0:index.diagonal[b]);
        r.radius=d->color->config().epsilon+2*diameter;r.nearAxis=l.domain().axisymmetric()&&r.center.x<r.radius;
        r.nearBoundary=index.nearWall(r.center,r.radius);
        r.crossesAmr=a>=0&&b>=0&&index.level[a]!=index.level[b];
        long double scale=0,estimatedError=0;
        for(HydroPressureLoad::Matrix::InnerIterator it(d->pressure.matrix(),u);it;++it){
            scale+=std::abs(it.value()*d->color->colors()[it.col()]);
            estimatedError+=std::abs(it.value())*cells[it.col()].absoluteErrorEstimate/d->color->phaseGeometry().metrics().cellVolumes()[it.col()];}
        r.gradientRoundoffScale=static_cast<double>(cfg.relativeGradientRoundoff*scale);
        r.colorGradientErrorEstimate=static_cast<double>(estimatedError);
        r.subResolutionGradient=d->gradient[u]!=0&&std::abs(d->gradient[u])<=r.gradientRoundoffScale;
        weights.clear();curvatures.clear();
        const bool faceCrossesAmr=r.crossesAmr;
        d->report.nearBoundaryFaces+=r.nearBoundary;d->report.nearAxisFaces+=r.nearAxis;d->report.amrFaces+=r.crossesAmr;
        if(d->gradient[u]==0){r.status=HydroCsfSupportStatus::ZeroGradient;continue;}
        try{projectAdaptive(index,r,cfg,index.level[p],weights,curvatures,d->panels[u]);}
        catch(const std::exception& e){std::ostringstream message;message<<std::setprecision(17)<<e.what()<<", velocity="<<u<<
            ", center=("<<r.center.x<<','<<r.center.y<<"), surface="<<r.surfaceQueries<<", geometry="<<r.geometryQueries<<
            ", depth="<<r.deepestPanel<<", weight="<<r.totalWeight<<", weight_error="<<r.weightErrorEstimate<<
            ", curvature_error="<<r.curvatureErrorEstimate<<", scale="<<r.curvatureScale;
            throw std::invalid_argument(message.str());}
        if(!faceCrossesAmr&&r.crossesAmr)++d->report.amrFaces;
        d->report.candidateQueries+=r.candidateQueries;d->report.visibilityRejected+=r.visibilityRejected;
        if(r.supportCount){require(std::isfinite(r.totalWeight)&&r.totalWeight>0,"unrepresentable CSF support measure");
            d->averages[u]=positiveCurvatureAverage(weights,curvatures);r.curvatureDefined=true;r.status=HydroCsfSupportStatus::Projected;++d->report.projectedFaces;
            const double coefficient=sigma*d->averages[u];require(std::isfinite(coefficient)&&(coefficient!=0||d->averages[u]==0),"unrepresentable projected CSF coefficient");
            d->rhs[u]=coefficient*d->gradient[u];
            require(std::isfinite(d->rhs[u])&&(d->rhs[u]!=0||coefficient==0||d->gradient[u]==0),"unrepresentable projected CSF momentum RHS");}
        else if(d->gradient[u]==0){r.status=HydroCsfSupportStatus::ZeroGradient;}
        else if(r.boundary&&cells[p].transitionVolumeUpper==0){
            // C1 on a certified band-free prescribed wall is a pressure
            // reaction term, not evidence of a physical interface here.
            r.status=HydroCsfSupportStatus::PureBoundaryReaction;++d->report.pureBoundaryRows;}
        else throw std::invalid_argument("[HydroBalancedCsfLoad] unsupported active CSF face without local selected/visible curvature support, velocity="+std::to_string(u));
    }
    return HydroBalancedCsfLoad(std::move(d));
}

HydroBalancedCsfLoad HydroBalancedCsfLoad::geometricCapillary(const QuadtreeMesh& mesh,HydroPressureLoad pressure,
    HydroPhaseGeometry phase,double sigma,HydroCsfCurvatureConfig cfg)
{
    require(std::isfinite(sigma)&&sigma>=0,"invalid surface tension");validateCurvatureConfig(cfg);
    require(pressure.layout().sameIdentity(phase.layout()),"different pressure/geometric phase layout");
    require(&pressure.divergence().matrix()==&phase.metrics().divergence().matrix(),"different authoritative divergence owner");
    phase.validateFor(mesh,phase.geometryIdentity(),phase.timeLayer());
    require(phase.publication().state().primary().transport==Interface::HermiteGalsTransport::FirstJet,
        "geometric CSF requires Pure FirstJet publication");
    auto d=std::make_shared<Data>(std::move(pressure),std::move(phase),sigma);d->capillary=true;d->config=cfg;
    const auto& l=d->pressure.layout();const int nu=l.velocityCount(),np=l.pressureCount();
    const auto& alpha=d->phase.liquidFractions();
    for(double value:alpha)require(std::isfinite(value)&&value>=0&&value<=1,"unbounded geometric phase fraction");
    d->gradient=d->pressure.apply(mesh,HydroNativeValues(l,std::vector<double>(nu,0),alpha));
    d->rhs.assign(nu,0);d->averages.assign(nu,0);d->support.resize(nu);d->panels.resize(nu);
    if(sigma==0){d->report.projectionSkippedSigmaZero=true;for(auto& r:d->support)r.status=HydroCsfSupportStatus::SigmaZero;
        return HydroBalancedCsfLoad(std::move(d));}
    std::vector<std::vector<std::pair<std::size_t,std::size_t>>> graphs(np);
    const auto& publication=d->phase.publication();
    for(std::size_t c=0;c<publication.topology().curves.size();++c)
        for(std::size_t g=0;g<publication.topology().curves[c].graphs().size();++g){
            const int p=l.pressureDof(publication.topology().curves[c].graphs()[g].nativeCellIndex);
            if(p>=0)graphs[p].push_back({c,g});
        }
    std::vector<GeometricCellCurvature> cells(np);
    std::vector<int> components(np);std::iota(components.begin(),components.end(),0);
    const auto root=[&](int p){int r=p;while(components[r]!=r)r=components[r];
        while(components[p]!=p){const int next=components[p];components[p]=r;p=next;}return r;};
    for(int u=0;u<nu;++u){const auto f=l.domain().face(mesh,l.nativeFace(u));
        if(f.negativeCell>=0&&f.positiveCell>=0&&d->phase.metrics().faceMeasures()[u]>0)
            components[root(l.pressureDof(f.positiveCell))]=root(l.pressureDof(f.negativeCell));}
    std::map<int,int> numbered;std::vector<int> component(np);
    for(int p=0;p<np;++p){const int r=root(p);if(!numbered.count(r))numbered[r]=static_cast<int>(numbered.size());component[p]=numbered[r];}
    for(int u=0;u<nu;++u){auto& r=d->support[u];const auto face=l.domain().face(mesh,l.nativeFace(u));
        const int a=face.negativeCell<0?-1:l.pressureDof(face.negativeCell),b=face.positiveCell<0?-1:l.pressureDof(face.positiveCell);
        const int p=a>=0?a:b;require(p>=0,"geometric CSF face has no native pressure cell");
        r.center=face.center;r.boundary=a<0||b<0;r.nearBoundary=r.boundary;r.pressureComponent=component[p];
        r.nearAxis=l.domain().axisymmetric()&&r.center.x==0;
        r.crossesAmr=a>=0&&b>=0&&mesh.cells(mesh.cellIndexToLabel[face.negativeCell]).level()!=mesh.cells(mesh.cellIndexToLabel[face.positiveCell]).level();
        d->report.nearBoundaryFaces+=r.nearBoundary;d->report.nearAxisFaces+=r.nearAxis;d->report.amrFaces+=r.crossesAmr;
        long double gradientScale=0;
        for(HydroPressureLoad::Matrix::InnerIterator it(d->pressure.matrix(),u);it;++it)gradientScale+=std::abs(it.value()*alpha[it.col()]);
        r.gradientRoundoffScale=static_cast<double>(cfg.relativeGradientRoundoff*gradientScale);
        r.subResolutionGradient=d->gradient[u]!=0&&std::abs(d->gradient[u])<=r.gradientRoundoffScale;
        if(d->gradient[u]==0){r.status=HydroCsfSupportStatus::ZeroGradient;continue;}
        std::vector<int> adjacent;for(int cell:{a,b})if(cell>=0&&!graphs[cell].empty())adjacent.push_back(cell);
        if(adjacent.empty()){
            if(r.boundary){r.status=HydroCsfSupportStatus::PureBoundaryReaction;++d->report.pureBoundaryRows;continue;}
            throw std::invalid_argument("[HydroBalancedCsfLoad] active geometric CSF face without adjacent cell interface, velocity="+std::to_string(u));
        }
        std::vector<double> equalWeights,averages;
        for(int cell:adjacent){auto& local=cells[cell];
            if(!local.defined){try{geometricCell(d->phase,cell,graphs[cell],local,cfg);}
                catch(const std::exception& e){throw std::invalid_argument(std::string(e.what())+", velocity="+std::to_string(u)+
                    ", geometry="+std::to_string(d->phase.geometryIdentity().value()));}
                d->report.sourceNodes+=local.queries;d->report.selectedNodes+=local.supportCount;}
            r.candidateQueries+=graphs[cell].size();r.surfaceQueries+=local.queries;r.supportCount+=local.supportCount;
            r.deepestPanel=std::max(r.deepestPanel,local.deepestPanel);r.totalWeight+=1;
            r.weightErrorEstimate+=local.weightError/local.area;r.anchoredMomentErrorEstimate+=local.momentError/local.area;
            r.curvatureErrorEstimate+=static_cast<double>(local.curvatureError/adjacent.size());
            r.curvatureScale=std::max(r.curvatureScale,local.scale);
            if(equalWeights.empty()){r.curvatureMin=local.minimum;r.curvatureMax=local.maximum;}
            else{r.curvatureMin=std::min(r.curvatureMin,local.minimum);r.curvatureMax=std::max(r.curvatureMax,local.maximum);}
            equalWeights.push_back(1);averages.push_back(local.average);
            for(const auto& panel:local.panels){d->panels[u].push_back({panel.component,panel.graph,panel.begin,panel.end,
                    panel.depth,1,false,panel.weightError/local.area,panel.momentError/local.area});
                for(const auto point:panel.fine.points)r.maxSupportDistance=std::max(r.maxSupportDistance,std::hypot(point.x-r.center.x,point.y-r.center.y));}
        }
        require(r.surfaceQueries<=cfg.maxSurfaceQueriesPerFace,"geometric face curvature surface query limit");
        require(r.candidateQueries<=cfg.maxCandidateQueriesPerFace,"geometric face curvature graph query limit");
        // Equal cell weights implement the local two-cell arithmetic mean,
        // independent of the different surface areas inside the two cells.
        d->averages[u]=positiveCurvatureAverage(equalWeights,averages);
        r.quadratureAnchor=averages.front();r.quadratureAnchorDefined=true;
        r.curvatureDefined=true;r.status=HydroCsfSupportStatus::Projected;r.acceptedPanels=d->panels[u].size();
        ++d->report.projectedFaces;d->report.candidateQueries+=r.candidateQueries;
        const double coefficient=sigma*d->averages[u];
        require(std::isfinite(coefficient)&&(coefficient!=0||d->averages[u]==0),"unrepresentable geometric CSF coefficient");
        d->rhs[u]=coefficient*d->gradient[u];
        require(std::isfinite(d->rhs[u])&&(d->rhs[u]!=0||coefficient==0||d->gradient[u]==0),"unrepresentable geometric CSF momentum RHS");
    }
    return HydroBalancedCsfLoad(std::move(d));
}


HydroBalancedCsfLoad HydroBalancedCsfLoad::heightFunctionCapillary(const QuadtreeMesh& mesh,
    HydroPressureLoad pressure,HydroPhaseGeometry phase,double sigma,HydroCsfCurvatureConfig cfg)
{
    require(std::isfinite(sigma)&&sigma>=0,"invalid surface tension");
    require(cfg.forceMethod==HydroCsfForceMethod::MotionAdjoint||cfg.forceMethod==HydroCsfForceMethod::HeightFunction,
        "invalid CSF force method");
    require(std::isfinite(cfg.relativeGradientRoundoff)&&cfg.relativeGradientRoundoff>0&&cfg.relativeGradientRoundoff<1,
        "invalid CSF gradient roundoff diagnostic");
    require(pressure.layout().sameIdentity(phase.layout()),"different pressure/geometric phase layout");
    require(&pressure.divergence().matrix()==&phase.metrics().divergence().matrix(),"different authoritative divergence owner");
    phase.validateFor(mesh,phase.geometryIdentity(),phase.timeLayer());
    const auto& primary=phase.publication().state().primary();
    require(primary.transport==Interface::HermiteGalsTransport::FirstJet,"height-function requires Pure FirstJet publication");
    require(!phase.layout().domain().axisymmetric(),"height-function requires XY geometry");
    require(phase.layout().pressureCount()==mesh.cellIndexNum()&&primary.geometryCells.size()==static_cast<std::size_t>(mesh.cellIndexNum()),
        "height-function requires complete native geometry and pressure domain");
    require(phase.layout().domain().uniformLeaves(),"height-function requires uniform XY mesh");
    require(phase.publication().state().geometry().snapshot().reconstruction()==LevelSetHermiteGeometry::Reconstruction::EnrichedQuartic,
        "height-function requires SingleQuartic scalar representation");
    auto d=std::make_shared<Data>(std::move(pressure),std::move(phase),sigma);
    d->capillary=true;d->kind=HydroCsfDiscretization::GeometricHeightFunction;
    cfg.forceMethod=HydroCsfForceMethod::HeightFunction;d->config=cfg;
    const auto& l=d->pressure.layout();const int nu=l.velocityCount(),np=l.pressureCount();
    d->gradient=d->pressure.apply(mesh,HydroNativeValues(l,std::vector<double>(nu,0),d->phase.liquidFractions()));
    d->rhs.assign(nu,0);d->averages.assign(nu,0);d->support.resize(nu);d->panels.resize(nu);
    d->heightCurvature.assign(np,std::numeric_limits<double>::quiet_NaN());
    for(auto& support:d->support)support.curvatureErrorEstimateAvailable=false;
    if(sigma==0){d->report.projectionSkippedSigmaZero=true;for(auto& r:d->support)r.status=HydroCsfSupportStatus::SigmaZero;
        return HydroBalancedCsfLoad(std::move(d));}
    HeightFunctionGrid heights(mesh,d->phase,d->report);
    for(int p=0;p<np;++p)if(heights.graphs[p]>0)d->heightCurvature[p]=heights.curvature(p);
    const auto& alpha=d->phase.liquidFractions();
    for(int u=0;u<nu;++u){auto& r=d->support[u];const auto face=l.domain().face(mesh,l.nativeFace(u));
        const int a=face.negativeCell<0?-1:l.pressureDof(face.negativeCell),b=face.positiveCell<0?-1:l.pressureDof(face.positiveCell);
        r.center=face.center;r.boundary=a<0||b<0;r.nearBoundary=r.boundary;r.pressureComponent=0;
        d->report.nearBoundaryFaces+=r.nearBoundary;
        long double gradientScale=0;
        for(HydroPressureLoad::Matrix::InnerIterator it(d->pressure.matrix(),u);it;++it)gradientScale+=std::abs(it.value()*alpha[it.col()]);
        r.gradientRoundoffScale=static_cast<double>(cfg.relativeGradientRoundoff*gradientScale);
        r.subResolutionGradient=d->gradient[u]!=0&&std::abs(d->gradient[u])<=r.gradientRoundoffScale;
        if(d->gradient[u]==0){r.status=HydroCsfSupportStatus::ZeroGradient;continue;}
        std::vector<double> weights,curvatures;
        for(int p:{a,b})if(p>=0&&std::isfinite(d->heightCurvature[p])){weights.push_back(1);curvatures.push_back(d->heightCurvature[p]);}
        if(curvatures.empty()){
            if(r.boundary){r.status=HydroCsfSupportStatus::PureBoundaryReaction;++d->report.pureBoundaryRows;continue;}
            throw std::invalid_argument("[HydroBalancedCsfLoad] active height-function CSF face without adjacent interface curvature, velocity="+std::to_string(u));}
        d->averages[u]=positiveCurvatureAverage(weights,curvatures);
        const auto range=std::minmax_element(curvatures.begin(),curvatures.end());r.curvatureMin=*range.first;r.curvatureMax=*range.second;
        r.totalWeight=static_cast<long double>(weights.size());r.supportCount=curvatures.size();r.curvatureDefined=true;r.status=HydroCsfSupportStatus::Projected;
        r.curvatureScale=std::max(std::abs(r.curvatureMin),std::abs(r.curvatureMax));++d->report.projectedFaces;
        const double coefficient=sigma*d->averages[u];
        require(std::isfinite(coefficient)&&(coefficient!=0||d->averages[u]==0),"unrepresentable height-function CSF coefficient");
        d->rhs[u]=coefficient*d->gradient[u];
        require(std::isfinite(d->rhs[u])&&(d->rhs[u]!=0||coefficient==0||d->gradient[u]==0),"unrepresentable height-function CSF momentum RHS");
    }
    return HydroBalancedCsfLoad(std::move(d));
}


HydroBalancedCsfLoad HydroBalancedCsfLoad::localPlicHeightFunctionCapillary(const QuadtreeMesh& mesh,
    HydroPressureLoad pressure,HydroPhaseGeometry phase,double sigma,HydroCsfCurvatureConfig cfg)
{return localPlicHeightFunctionCapillaryImpl(mesh,std::move(pressure),std::move(phase),sigma,cfg,nullptr);}
HydroBalancedCsfLoad HydroBalancedCsfLoad::localPlicHeightFunctionCapillary(const QuadtreeMesh& mesh,
    HydroPressureLoad pressure,HydroPhaseGeometry phase,double sigma,const HydroVofSouthEven90Boundary& boundary,HydroCsfCurvatureConfig cfg)
{return localPlicHeightFunctionCapillaryImpl(mesh,std::move(pressure),std::move(phase),sigma,cfg,&boundary);}
const HydroVofSouthEven90Boundary* HydroBalancedCsfLoad::southEven90Boundary() const noexcept
{return data_->southEven90?&*data_->southEven90:nullptr;}
HydroBalancedCsfLoad HydroBalancedCsfLoad::localPlicHeightFunctionCapillaryImpl(const QuadtreeMesh& mesh,
    HydroPressureLoad pressure,HydroPhaseGeometry phase,double sigma,HydroCsfCurvatureConfig cfg,const HydroVofSouthEven90Boundary* boundary)
{
    require(std::isfinite(sigma)&&sigma>=0,"invalid surface tension");
    require(cfg.forceMethod==HydroCsfForceMethod::MotionAdjoint||cfg.forceMethod==HydroCsfForceMethod::HeightFunction,
        "invalid CSF force method");
    require(std::isfinite(cfg.relativeGradientRoundoff)&&cfg.relativeGradientRoundoff>0&&cfg.relativeGradientRoundoff<1,
        "invalid CSF gradient roundoff diagnostic");
    require(pressure.layout().sameIdentity(phase.layout()),"different pressure/geometric phase layout");
    require(&pressure.divergence().matrix()==&phase.metrics().divergence().matrix(),"different authoritative divergence owner");
    require(phase.model()==HydroPhaseGeometry::Model::GeometricVof,"local PLIC height-function requires a genuine VOF phase owner");
    phase.validateVofFor(mesh,phase.vofIdentity(),phase.timeLayer());
    require(phase.vofPublication().geometryKind()==Interface::GeometricVofGeometryKind::LocalPlic,
        "local PLIC height-function requires LocalPlic geometry");
    const auto& sourceDomain=phase.vofPublication().domain();
    require(phase.vofPublication().plicCells().size()==static_cast<std::size_t>(mesh.cellIndexNum())&&
        sourceDomain.cellRowByIndex().size()==static_cast<std::size_t>(mesh.cellIndexNum())&&
        phase.layout().domain().cellRowByIndex()==sourceDomain.cellRowByIndex()&&
        phase.layout().domain().cellIndexByRow()==sourceDomain.cellIndexByRow()&&
        phase.layout().domain().nativeFaceCount()==sourceDomain.nativeFaceCount(),
        "local PLIC height-function requires the exact selected native geometry/pressure maps");
    auto d=std::make_shared<Data>(std::move(pressure),std::move(phase),sigma);
    d->capillary=true;d->kind=HydroCsfDiscretization::GeometricHeightFunction;
    cfg.forceMethod=HydroCsfForceMethod::HeightFunction;d->config=cfg;
    if(boundary){boundary->validateFor(mesh,d->phase.vofPublication());d->southEven90=*boundary;}
    const auto& l=d->pressure.layout();const int nu=l.velocityCount(),np=l.pressureCount();
    // This force color is the same actual PLIC geometry used by the height
    // columns. Independent inventory alpha remains the transport/mass view.
    const auto& local=d->phase.vofPublication().plicCells();std::vector<double> color(np);
    for(int p=0;p<np;++p){const int native=l.nativeCell(p);const auto& cell=local[native];
        const auto capacity=d->phase.metrics().cellVolumes()[p];
        color[p]=cell.measured[0]==0?0:cell.measured[1]==0?1:static_cast<double>(cell.measured[0]/capacity);
        if(!(std::isfinite(color[p])&&color[p]>=0&&color[p]<=1))throw std::invalid_argument(
            "[HydroBalancedCsfLoad] invalid actual PLIC force color: "+localPlicSourceDescription(native,cell));}
    d->gradient=d->pressure.apply(mesh,HydroNativeValues(l,std::vector<double>(nu,0),color));
    d->rhs.assign(nu,0);d->averages.assign(nu,0);d->support.resize(nu);d->panels.resize(nu);
    d->heightCurvature.assign(np,std::numeric_limits<double>::quiet_NaN());d->heightJets.resize(np);
    d->heightEffectiveCurvature=d->heightCurvature;d->heightCurvatureMethods.assign(np,HydroCsfCellCurvatureMethod::Unavailable);
    for(auto& support:d->support)support.curvatureErrorEstimateAvailable=false;
    if(sigma==0){d->report.projectionSkippedSigmaZero=true;for(auto& r:d->support)r.status=HydroCsfSupportStatus::SigmaZero;
        return HydroBalancedCsfLoad(std::move(d));}
    LocalPlicHeightGrid heights(mesh,d->phase,d->report,boundary);
    for(int p=0;p<np;++p)if(heights.mixed[p]){
        d->heightEffectiveCurvature[p]=heights.curvature(p,mesh);
        if(heights.methods[p]!=HydroCsfCellCurvatureMethod::UnderresolvedIsolatedFragmentZero&&
            heights.methods[p]!=HydroCsfCellCurvatureMethod::Unavailable)
            d->heightCurvature[p]=d->heightEffectiveCurvature[p];}
    d->heightCurvatureMethods=std::move(heights.methods);d->heightJets=std::move(heights.jets);
    const auto& alpha=color; // The gradient roundoff scale uses this SAME force color.
    for(int u=0;u<nu;++u){auto& r=d->support[u];const auto face=l.domain().face(mesh,l.nativeFace(u));
        const int a=face.negativeCell<0?-1:l.pressureDof(face.negativeCell),b=face.positiveCell<0?-1:l.pressureDof(face.positiveCell);
        r.center=face.center;r.boundary=a<0||b<0;r.nearBoundary=r.boundary;r.pressureComponent=0;
        d->report.nearBoundaryFaces+=r.nearBoundary;
        r.crossesAmr=face.negativeCell>=0&&face.positiveCell>=0&&
            mesh.cells(mesh.cellIndexToLabel[face.negativeCell]).level()!=mesh.cells(mesh.cellIndexToLabel[face.positiveCell]).level();
        d->report.amrFaces+=r.crossesAmr;
        long double gradientScale=0;
        for(HydroPressureLoad::Matrix::InnerIterator it(d->pressure.matrix(),u);it;++it)gradientScale+=std::abs(it.value()*alpha[it.col()]);
        r.gradientRoundoffScale=static_cast<double>(cfg.relativeGradientRoundoff*gradientScale);
        r.subResolutionGradient=d->gradient[u]!=0&&std::abs(d->gradient[u])<=r.gradientRoundoffScale;
        if(d->gradient[u]==0){r.status=HydroCsfSupportStatus::ZeroGradient;continue;}
        std::vector<double> weights,curvatures;bool nativeMaskedNeighbor=false;
        for(int p:{a,b})if(p>=0&&std::isfinite(d->heightEffectiveCurvature[p])){
            r.underresolvedFragmentClosure=r.underresolvedFragmentClosure||
                d->heightCurvatureMethods[p]==HydroCsfCellCurvatureMethod::UnderresolvedIsolatedFragmentZero||
                d->heightCurvatureMethods[p]==HydroCsfCellCurvatureMethod::Unavailable;
            // The owner-declared Native2 P0 closure has no physical kappa.
            // Preserve its provenance without diluting a valid neighbour's
            // constant face potential by averaging with an unavailable zero.
            if(heights.nativeSupportChart&&heights.subgrid[l.nativeCell(p)]&&
                d->heightCurvatureMethods[p]==HydroCsfCellCurvatureMethod::Unavailable){
                nativeMaskedNeighbor=true;continue;}
            weights.push_back(1);curvatures.push_back(d->heightEffectiveCurvature[p]);}
        if(curvatures.empty()&&nativeMaskedNeighbor){
            // No physical donor: retain only the explicit owner-masked face
            // closure that the old finite-effective gather already supplied.
            curvatures.push_back(0);weights.push_back(1);
        }
        if(curvatures.empty()&&(boundary||d->phase.vofPublication().pinnedCorner()||heights.nativeSupportChart)&&!r.boundary){
            // An actual shared pure-phase face has its own chart position;
            // never label either pure pressure cell as a mixed HF source.
            curvatures.push_back(heights.coincidentFaceCurvature(mesh,u));weights.push_back(1);r.coincidentPhaseHeight=true;
        }
        if(curvatures.empty()){
            if(r.boundary){r.status=HydroCsfSupportStatus::PureBoundaryReaction;++d->report.pureBoundaryRows;continue;}
            std::string details;for(int p:{a,b})if(p>=0){const int native=l.nativeCell(p);details+="; "+localPlicSourceDescription(native,local[native]);}
            throw std::invalid_argument("[HydroBalancedCsfLoad] active height-function CSF face without adjacent interface curvature, velocity="+std::to_string(u)+details);}
        d->averages[u]=positiveCurvatureAverage(weights,curvatures);
        const auto range=std::minmax_element(curvatures.begin(),curvatures.end());r.curvatureMin=*range.first;r.curvatureMax=*range.second;
        r.totalWeight=static_cast<long double>(weights.size());r.supportCount=curvatures.size();r.curvatureDefined=true;r.status=HydroCsfSupportStatus::Projected;
        r.curvatureScale=std::max(std::abs(r.curvatureMin),std::abs(r.curvatureMax));++d->report.projectedFaces;
        const double coefficient=sigma*d->averages[u];
        require(std::isfinite(coefficient)&&(coefficient!=0||d->averages[u]==0),"unrepresentable height-function CSF coefficient");
        d->rhs[u]=coefficient*d->gradient[u];
        require(std::isfinite(d->rhs[u])&&(d->rhs[u]!=0||coefficient==0||d->gradient[u]==0),"unrepresentable height-function CSF momentum RHS");
    }
    return HydroBalancedCsfLoad(std::move(d));
}


HydroBalancedCsfLoad HydroBalancedCsfLoad::heightFunctionIncrementCapillary(const QuadtreeMesh& mesh,
    HydroPressureLoad pressure,HydroPhaseGeometry phase,double sigma,const HydroSpatialAdvector& motion,
    const Interface::HermiteGalsAmrClosureSensitivity* mixed,HydroCsfCurvatureConfig cfg)
{
    require(cfg.forceMethod==HydroCsfForceMethod::MotionAdjoint||cfg.forceMethod==HydroCsfForceMethod::HeightFunction
        ||cfg.forceMethod==HydroCsfForceMethod::HeightFunctionIncrement,"invalid CSF force method");
    require(phase.publication().state().primary().gradientUpdate==Interface::HermiteGalsGradientUpdate::CompatibleIncrement,
        "height-function increment requires actual compatible-increment publication");
    const auto& geometry=motion.momentumGeometry();geometry.validateMesh(mesh);
    require(geometry.layout().sameIdentity(pressure.layout()),"different increment motion/native pressure layout");
    require(&geometry.metrics().divergence().matrix()==&pressure.divergence().matrix(),"different increment divergence owner");
    motion.validateFor(mesh,geometry,phase.geometryIdentity(),phase.timeLayer());
    require(motion.reconstruction().tangentialModes().empty(),"height-function increment requires native-only motion");
    const auto& regular=motion.reconstruction().regularGeometry();regular.validateFor(mesh,geometry);
    require(regular.usesCompatibleTensor()&&!phase.layout().domain().axisymmetric(),"height-function increment requires complete uniform XY");
    const auto* boundary=motion.boundary();
    require(boundary&&boundary->completeBoundaryLines(),"height-function increment requires complete homogeneous closed boundary");
    boundary->validateFor(mesh,geometry);
    const auto& mask=boundary->velocityDofs();
    for(int u=0;u<phase.layout().velocityCount();++u){const auto face=phase.layout().domain().face(mesh,phase.layout().nativeFace(u));
        if(face.negativeCell<0||face.positiveCell<0)
            require(std::find(mask.begin(),mask.end(),u)!=mask.end(),"height-function increment requires all closed wall faces");}
    if(const auto* bound=motion.reconstruction().phaseGeometry())bound->validateFor(mesh,phase.geometryIdentity(),phase.timeLayer());
    auto baseConfig=cfg;baseConfig.forceMethod=HydroCsfForceMethod::HeightFunction;
    const auto base=heightFunctionCapillary(mesh,std::move(pressure),std::move(phase),sigma,baseConfig);
    auto d=std::make_shared<Data>(*base.data_);d->kind=HydroCsfDiscretization::GeometricHeightFunctionIncrement;
    cfg.forceMethod=HydroCsfForceMethod::HeightFunctionIncrement;d->config=cfg;
    d->source=motion.reconstruction();d->boundary=*boundary;d->projectedMotion=motion;
    const auto& layout=d->phase.layout();const int nu=layout.velocityCount(),np=layout.pressureCount();
    d->incrementReferenceLoad.assign(nu,0);d->incrementResidualLoad.assign(nu,0);d->incrementDeltaMean.assign(np,0);
    // This is a distinct fixed accepted-GL8 discretization. No old adaptive
    // native-traction certificate, HF truncation estimate or face multiplier
    // is claimed for the corrected force.
    for(auto& r:d->support)r.curvatureErrorEstimateAvailable=false;
    if(sigma==0)return HydroBalancedCsfLoad(std::move(d));
    require(mixed,"height-function increment lacks actual mixed-closure owner");
    mixed->validateFor(mesh,d->phase.publication().state());d->mixed=*mixed;
    d->incrementRates=std::make_shared<ProjectedRateBasis>(mesh,d->phase,*d->projectedMotion,*d->mixed,true);
    d->incrementRateMutex=std::make_shared<std::mutex>();
    int count=0;const auto component=adjointPressureComponents(mesh,d->phase,count);
    struct Node {const Interface::HermiteGalsSurfaceQuadrature* q;int pressure,component;long double delta;};
    std::vector<Node> nodes;std::vector<long double> areas(count),moments(count);
    const auto& publication=d->phase.publication();const auto clock=d->phase.timeLayer();
    for(std::size_t c=0;c<publication.topology().curves.size();++c)
        for(std::size_t g=0;g<publication.topology().curves[c].graphs().size();++g){
            const int p=layout.pressureDof(publication.topology().curves[c].graphs()[g].nativeCellIndex);
            require(p>=0&&std::isfinite(d->heightCurvature[p]),"increment interface lacks same-owner height curvature");
            for(const auto& q:publication.quadrature(c,g)){
                const auto time=q.sample.timeLayer;
                require(q.sample.geometryIdentity==d->phase.geometryIdentity()&&time.step==clock.step
                    &&time.begin==clock.begin&&time.end==clock.end&&time.position==clock.position
                    &&q.sample.frame.nativeCellIndex==layout.nativeCell(p),"increment cached quadrature has different geometry/full clock");
                require(q.physicalWeight>0&&std::isfinite(q.physicalWeight)&&std::isfinite(q.sample.curvature),"invalid increment cached physical sample");
                const long double delta=static_cast<long double>(q.sample.curvature)-d->heightCurvature[p];
                const int id=component[p];areas[id]+=q.physicalWeight;moments[id]+=q.physicalWeight*delta;
                nodes.push_back({&q,p,id,delta});
            }
        }
    std::vector<long double> means(count);std::vector<double> potential(np);
    for(int id=0;id<count;++id)if(areas[id]>0){means[id]=moments[id]/areas[id];
        require(std::isfinite(means[id]),"nonfinite increment mean curvature difference");}
    const auto& alpha=d->phase.liquidFractions();
    for(int p=0;p<np;++p){d->incrementDeltaMean[p]=static_cast<double>(means[component[p]]);
        potential[p]=sigma*d->incrementDeltaMean[p]*alpha[p];
        require(std::isfinite(potential[p]),"nonfinite increment reference pressure");}
    d->incrementReferenceLoad=d->pressure.apply(mesh,HydroNativeValues(layout,std::vector<double>(nu,0),potential));
    std::vector<long double> residual(nu);
    // One scatter per accepted node. No new curve, root discovery, facewise
    // quadrature or recursively resampled traction integral is introduced.
    for(const auto& node:nodes){const auto weights=d->incrementRates->normal(node.q->sample.frame);
        const long double traction=-static_cast<long double>(sigma)*node.q->physicalWeight*(node.delta-means[node.component]);
        require(std::isfinite(traction),"nonfinite increment residual traction");
        for(const auto& entry:weights){require(entry.first>=0&&entry.first<nu,"increment source native DOF out of range");
            residual[entry.first]+=traction*entry.second;}
    }
    d->report.incrementNodes=nodes.size();d->report.incrementRatePoints=d->incrementRates->phiValueRates.size();
    for(int u=0;u<nu;++u){d->incrementResidualLoad[u]=static_cast<double>(residual[u]);
        const double increment=d->incrementReferenceLoad[u]+d->incrementResidualLoad[u];d->rhs[u]+=increment;
        require(std::isfinite(d->incrementResidualLoad[u])&&std::isfinite(increment)&&std::isfinite(d->rhs[u]),"nonfinite increment native momentum RHS");
        if(increment!=0)++d->report.incrementLoadDofs;
    }
    return HydroBalancedCsfLoad(std::move(d));
}

HydroBalancedCsfLoad HydroBalancedCsfLoad::adjointCapillary(const QuadtreeMesh& mesh,HydroPressureLoad pressure,
    HydroPhaseGeometry phase,double sigma,const HydroAdvectorReconstruction& reconstruction,
    const HydroAdvectorBoundary* boundary,HydroCsfCurvatureConfig cfg)
{return weakCapillary(mesh,std::move(pressure),std::move(phase),sigma,reconstruction,boundary,cfg,nullptr,nullptr);}

HydroBalancedCsfLoad HydroBalancedCsfLoad::projectedMotionCapillary(const QuadtreeMesh& mesh,HydroPressureLoad pressure,
    HydroPhaseGeometry phase,double sigma,const HydroSpatialAdvector& motion,
    const Interface::HermiteGalsAmrClosureSensitivity* mixed,HydroCsfCurvatureConfig cfg)
{return weakCapillary(mesh,std::move(pressure),std::move(phase),sigma,motion.reconstruction(),motion.boundary(),cfg,&motion,mixed);}

HydroBalancedCsfLoad HydroBalancedCsfLoad::adjointInterfaceLoad(const QuadtreeMesh& mesh,HydroPressureLoad pressure,
    HydroPhaseGeometry phase,double sigma,const Interface::HermiteGalsSurfaceVectorView& external,
    const HydroAdvectorReconstruction& reconstruction,const HydroAdvectorBoundary* boundary,HydroCsfCurvatureConfig cfg)
{return weakCapillary(mesh,std::move(pressure),std::move(phase),sigma,reconstruction,boundary,cfg,nullptr,nullptr,&external);}
HydroBalancedCsfLoad HydroBalancedCsfLoad::projectedMotionInterfaceLoad(const QuadtreeMesh& mesh,HydroPressureLoad pressure,
    HydroPhaseGeometry phase,double sigma,const Interface::HermiteGalsSurfaceVectorView& external,
    const HydroSpatialAdvector& motion,const Interface::HermiteGalsAmrClosureSensitivity* mixed,HydroCsfCurvatureConfig cfg)
{return weakCapillary(mesh,std::move(pressure),std::move(phase),sigma,motion.reconstruction(),motion.boundary(),cfg,&motion,mixed,&external);}

HydroBalancedCsfLoad HydroBalancedCsfLoad::weakCapillary(const QuadtreeMesh& mesh,HydroPressureLoad pressure,
    HydroPhaseGeometry phase,double sigma,const HydroAdvectorReconstruction& reconstruction,
    const HydroAdvectorBoundary* boundary,HydroCsfCurvatureConfig cfg,const HydroSpatialAdvector* motion,
    const Interface::HermiteGalsAmrClosureSensitivity* mixed,const Interface::HermiteGalsSurfaceVectorView* external)
{
    require(std::isfinite(sigma)&&sigma>=0,"invalid surface tension");validateCurvatureConfig(cfg);
    if(external)external->validateFor(mesh,phase.publication(),phase.timeLayer());
    require(pressure.layout().sameIdentity(phase.layout()),"different pressure/adjoint phase layout");
    require(&pressure.divergence().matrix()==&phase.metrics().divergence().matrix(),"different authoritative divergence owner");
    phase.validateFor(mesh,phase.geometryIdentity(),phase.timeLayer());
    require(phase.publication().state().primary().transport==Interface::HermiteGalsTransport::FirstJet,
        "adjoint CSF requires Pure FirstJet publication");
    const auto& geometry=reconstruction.momentumGeometry();geometry.validateMesh(mesh);
    require(geometry.layout().sameIdentity(pressure.layout()),"different adjoint source layout");
    require(&geometry.metrics().divergence().matrix()==&pressure.divergence().matrix(),"different adjoint source divergence owner");
    require(reconstruction.tangentialModes().empty(),"adjoint native load does not support independent tangential modes");
    if(boundary)boundary->validateFor(mesh,geometry);
    if(const auto* bound=reconstruction.phaseGeometry())bound->validateFor(mesh,phase.geometryIdentity(),phase.timeLayer());
    const auto& regular=reconstruction.regularGeometry();regular.validateFor(mesh,geometry);
    auto d=std::make_shared<Data>(std::move(pressure),std::move(phase),sigma);d->capillary=true;d->config=cfg;
    d->kind=motion?HydroCsfDiscretization::GeometricProjectedAdjoint:HydroCsfDiscretization::GeometricAdjoint;
    d->source=reconstruction;if(boundary)d->boundary=*boundary;
    if(external)d->externalSurface=*external;
    if(motion){
        motion->validateFor(mesh,geometry,d->phase.geometryIdentity(),d->phase.timeLayer());
        d->projectedMotion=*motion;
        if(sigma>0||external){
            require(mixed!=nullptr,"projected load lacks mixed closure owner");
            mixed->validateFor(mesh,d->phase.publication().state());
            require(regular.usesCompatibleTensor()&&!d->phase.layout().domain().axisymmetric(),
                "projected load requires complete uniform XY source");
            require(boundary&&boundary->completeBoundaryLines(),"projected load requires complete homogeneous closed boundary");
            const auto& mask=boundary->velocityDofs();
            for(int u=0;u<d->phase.layout().velocityCount();++u){
                const auto f=d->phase.layout().domain().face(mesh,d->phase.layout().nativeFace(u));
                if(f.negativeCell<0||f.positiveCell<0)
                    require(std::find(mask.begin(),mask.end(),u)!=mask.end(),"projected load requires all closed wall faces");
            }
            d->mixed=*mixed;
        }
    }
    const auto& l=d->pressure.layout();const int nu=l.velocityCount(),np=l.pressureCount();
    const auto& alpha=d->phase.liquidFractions();for(double a:alpha)require(std::isfinite(a)&&a>=0&&a<=1,"unbounded geometric phase fraction");
    d->gradient=d->pressure.apply(mesh,HydroNativeValues(l,std::vector<double>(nu,0),alpha));
    d->rhs.assign(nu,0);d->referenceLoad.assign(nu,0);d->residualLoad.assign(nu,0);d->tangentialLoad.assign(nu,0);
    d->referenceCurvature.assign(np,std::numeric_limits<double>::quiet_NaN());
    d->referenceStress.assign(np,std::numeric_limits<double>::quiet_NaN());
    d->support.resize(nu);d->panels.resize(nu);
    if(sigma==0&&!external){d->report.projectionSkippedSigmaZero=true;for(auto& r:d->support)r.status=HydroCsfSupportStatus::SigmaZero;
        return HydroBalancedCsfLoad(std::move(d));}
    std::optional<ProjectedRateBasis> projected;
    if(motion)projected.emplace(mesh,d->phase,*d->projectedMotion,*d->mixed);
    int componentCount=0;const auto component=adjointPressureComponents(mesh,d->phase,componentCount);
    std::vector<AdjointAnchor> anchors(componentCount);
    std::vector<std::optional<Detail::HydroRegularCellStencil>> stencils(np);
    std::vector<std::size_t> queries(np),requested(np),graphs(np);std::vector<double> scales(np),stressScales(np);
    std::vector<Vector2D> centers(nu);
    for(int u=0;u<nu;++u){const auto f=l.domain().face(mesh,l.nativeFace(u));centers[u]=f.center;
        auto& r=d->support[u];r.center=f.center;r.boundary=f.negativeCell<0||f.positiveCell<0;
        const int p=l.pressureDof(f.negativeCell>=0?f.negativeCell:f.positiveCell);r.pressureComponent=component[p];
        r.nearBoundary=r.boundary;r.nearAxis=l.domain().axisymmetric()&&f.center.x==0;
        r.crossesAmr=f.negativeCell>=0&&f.positiveCell>=0&&
            mesh.cells(mesh.cellIndexToLabel[f.negativeCell]).level()!=mesh.cells(mesh.cellIndexToLabel[f.positiveCell]).level();
        d->report.nearBoundaryFaces+=r.nearBoundary;d->report.nearAxisFaces+=r.nearAxis;d->report.amrFaces+=r.crossesAmr;
    }
    std::vector<AdjointPanel> accepted;const auto& publication=d->phase.publication();
    for(std::size_t c=0;c<publication.topology().curves.size();++c)
        for(std::size_t g=0;g<publication.topology().curves[c].graphs().size();++g){
            const int p=l.pressureDof(publication.topology().curves[c].graphs()[g].nativeCellIndex);
            require(p>=0,"adjoint interface graph outside actual pressure domain");
            require(++graphs[p]<=cfg.maxCandidateQueriesPerFace,"adjoint candidate query limit");
            if(!stencils[p]){stencils[p]=regular.cellStencil(mesh,p,boundary);
                const auto& cell=regular.cells()[p];scales[p]=1/std::min(cell.upper.x-cell.lower.x,cell.upper.y-cell.lower.y);}
            std::vector<AdjointPanel> panels;
            // Split every actual source knot in the graph's monotone chosen
            // coordinate using its original cell chart. Transverse source
            // knots remain controlled by the actual native-moment estimate.
            // Projected rates use the scalar's actual piece sensitivity.
            // Macro graphs have already been split at their half-cell seams;
            // source quarter knots were resolved at the primary queries and
            // are not additional knots of Q_Gamma along a scalar graph.
            const auto cuts=projected?std::vector<double>{0,1}:adjointChosenCoordinateCuts(d->phase,regular.cells()[p],c,g);
            const AdjointExternal externalRule{sigma,external,&stressScales[p]};
            for(std::size_t part=1;part<cuts.size();++part)
                panels.push_back(adjointPanel(d->phase,regular,*stencils[p],centers,c,g,cuts[part-1],cuts[part],0,
                    anchors[component[p]],queries[p],requested[p],scales[p],cfg,nullptr,projected?&*projected:nullptr,
                    external?&externalRule:nullptr));
            for(;;){
                const auto budget=adjointGraphBudget(panels);const auto score=adjointGraphScore(budget,scales[p],external?stressScales[p]:-1);
                if(score<=cfg.relativeQuadratureTolerance){
                    for(auto& panel:panels){panel.halves={};accepted.push_back(std::move(panel));}break;}
                const auto worst=std::max_element(panels.begin(),panels.end(),[&](const auto& a,const auto& b){
                    return adjointPanelPriority(a,budget,scales[p],external?stressScales[p]:-1)
                        <adjointPanelPriority(b,budget,scales[p],external?stressScales[p]:-1);});
                require(worst!=panels.end(),"adjoint graph has no refinable panel");auto old=std::move(*worst);
                if(old.depth>=cfg.maxQuadratureDepth){std::ostringstream message;message<<std::setprecision(17)
                    <<"[HydroBalancedCsfLoad] adjoint traction quadrature depth limit, pressure="<<p<<", graph="<<g
                    <<", fraction=("<<old.begin<<','<<old.end<<"), complete_graph_score="<<score;
                    throw std::invalid_argument(message.str());}
                const double middle=old.begin+(old.end-old.begin)/2;
                *worst=adjointPanel(d->phase,regular,*stencils[p],centers,c,g,old.begin,middle,old.depth+1,
                    anchors[component[p]],queries[p],requested[p],scales[p],cfg,&old.halves[0],projected?&*projected:nullptr,
                    external?&externalRule:nullptr);
                panels.push_back(adjointPanel(d->phase,regular,*stencils[p],centers,c,g,middle,old.end,old.depth+1,
                    anchors[component[p]],queries[p],requested[p],scales[p],cfg,&old.halves[1],projected?&*projected:nullptr,
                    external?&externalRule:nullptr));
            }
        }
    std::vector<long double> areas(componentCount),moments(componentCount),stressMoments(componentCount);
    for(const auto& panel:accepted){const int id=component[panel.pressure];areas[id]+=panel.fine.area;moments[id]+=panel.fine.moment;
        stressMoments[id]+=panel.fine.stressMoment;}
    std::vector<double> potential(np);
    for(int p=0;p<np;++p){const int id=component[p];
        if(areas[id]>0){require(anchors[id].defined,"adjoint component reference lacks actual surface sample");
            const double k=static_cast<double>(static_cast<long double>(anchors[id].curvature)+moments[id]/areas[id]);
            require(std::isfinite(k),"nonfinite actual component reference curvature");d->referenceCurvature[p]=k;
            d->referenceStress[p]=external?static_cast<double>(static_cast<long double>(anchors[id].normalStress)
                +stressMoments[id]/areas[id]):sigma*k;
            require(std::isfinite(d->referenceStress[p]),"nonfinite actual component reference stress");
            potential[p]=external?d->referenceStress[p]*alpha[p]:sigma*k*alpha[p];
            require(std::isfinite(potential[p]),"nonfinite adjoint reference pressure potential");}
    }
    d->referenceLoad=d->pressure.apply(mesh,HydroNativeValues(l,std::vector<double>(nu,0),potential));
    std::vector<long double> residual(nu),tangential(nu);
    for(const auto& panel:accepted){const int p=panel.pressure,id=component[p];
        const long double offset=external?static_cast<long double>(d->referenceStress[p])-anchors[id].normalStress
            :static_cast<long double>(d->referenceCurvature[p])-anchors[id].curvature;
        for(const auto& entry:panel.fine.native){const int u=entry.first;const auto& value=entry.second;
            residual[u]-=external?value.traction-offset*value.normal
                :static_cast<long double>(sigma)*(value.traction-offset*value.normal);
            tangential[u]+=value.tangent;
            auto& r=d->support[u];const bool first=!r.curvatureDefined;
            r.curvatureDefined=true;r.status=HydroCsfSupportStatus::Projected;
            r.curvatureMin=first?value.minimum:std::min(r.curvatureMin,value.minimum);
            r.curvatureMax=first?value.maximum:std::max(r.curvatureMax,value.maximum);
            r.supportCount+=value.nodes;r.surfaceQueries+=value.nodes;r.totalWeight+=value.normalScale;
            r.maxSupportDistance=std::max(r.maxSupportDistance,value.maxDistance);
            r.deepestPanel=std::max(r.deepestPanel,panel.depth);++r.candidateQueries;
            r.quadratureAnchor=d->referenceCurvature[p];r.quadratureAnchorDefined=true;
            r.weightErrorEstimate+=panel.areaError;r.anchoredMomentErrorEstimate+=panel.momentError;r.curvatureScale=std::max(r.curvatureScale,scales[p]);
            const auto error=panel.errors.find(u);if(error!=panel.errors.end()){
                r.normalMomentErrorEstimate+=error->second[0];
                r.tractionMomentErrorEstimate+=error->second[1]+std::abs(offset)*error->second[0];}
            d->panels[u].push_back({panel.component,panel.graph,panel.begin,panel.end,panel.depth,1,false,panel.areaError,panel.momentError});
        }
    }
    for(auto count:queries)d->report.sourceNodes+=count;
    for(int u=0;u<nu;++u){d->residualLoad[u]=static_cast<double>(residual[u]);d->tangentialLoad[u]=static_cast<double>(tangential[u]);
        d->rhs[u]=d->referenceLoad[u]+d->residualLoad[u];
        require(std::isfinite(d->rhs[u])&&std::isfinite(d->residualLoad[u])&&std::isfinite(d->tangentialLoad[u]),"nonfinite adjoint momentum RHS");
        auto& r=d->support[u];r.acceptedPanels=d->panels[u].size();
        require(r.surfaceQueries<=cfg.maxSurfaceQueriesPerFace&&r.candidateQueries<=cfg.maxCandidateQueriesPerFace,
            "adjoint native support query limit");
        d->report.selectedNodes+=r.supportCount;d->report.candidateQueries+=r.candidateQueries;
        if(r.curvatureDefined)++d->report.projectedFaces;
        else if(r.boundary&&d->gradient[u]!=0){r.status=HydroCsfSupportStatus::PureBoundaryReaction;++d->report.pureBoundaryRows;}
    }
    return HydroBalancedCsfLoad(std::move(d));
}

const HydroPressureLoad& HydroBalancedCsfLoad::pressureLoad() const noexcept{return data_->pressure;}
const HydroRegularizedColor& HydroBalancedCsfLoad::color() const
{require(data_->color.has_value(),"geometric phase-fraction CSF load has no regularized color");return *data_->color;}
const HydroPhaseGeometry& HydroBalancedCsfLoad::phaseGeometry() const noexcept{return data_->phase;}
HydroCsfDiscretization HydroBalancedCsfLoad::discretization() const noexcept{return data_->kind;}
const HydroDofLayout& HydroBalancedCsfLoad::layout() const noexcept{return data_->pressure.layout();}
double HydroBalancedCsfLoad::surfaceTension() const noexcept{return data_->sigma;}
double HydroBalancedCsfLoad::constantCurvature() const{require(!data_->capillary,"projected load has no prescribed constant curvature");return data_->curvature;}
bool HydroBalancedCsfLoad::isCapillary() const noexcept{return data_->capillary;}
const HydroCsfCurvatureConfig& HydroBalancedCsfLoad::curvatureConfig() const{require(data_->capillary,"constant load has no curvature projection config");return data_->config;}
const std::vector<double>& HydroBalancedCsfLoad::curvatureAverages() const
{require(data_->capillary,"constant load has no projected curvature");require(data_->kind!=HydroCsfDiscretization::GeometricAdjoint
    &&data_->kind!=HydroCsfDiscretization::GeometricProjectedAdjoint
    &&data_->kind!=HydroCsfDiscretization::GeometricHeightFunctionIncrement,
    "load has no single face curvature multiplier");return data_->averages;}
const std::vector<double>& HydroBalancedCsfLoad::heightFunctionFaceCurvatures() const
{require(data_->kind==HydroCsfDiscretization::GeometricHeightFunction
    ||data_->kind==HydroCsfDiscretization::GeometricHeightFunctionIncrement,"load has no height-function base face curvature");return data_->averages;}
std::array<long double,2> HydroBalancedCsfLoad::axisymmetricHeightCurvatures(const HydroCsfHeightJet& jet)
{
    require((jet.heightAxis==0||jet.heightAxis==1)&&(jet.outward==1||jet.outward==-1)&&
        std::isfinite(jet.height)&&jet.height>=0&&std::isfinite(jet.first)&&std::isfinite(jet.second)&&
        std::isfinite(jet.origin)&&std::isfinite(jet.tangentCoordinate),"invalid RZ height jet");
    const auto metric=std::hypot(1.L,jet.first),meridional=-jet.second/(metric*metric*metric);
    const long double radius=jet.heightAxis==0?jet.origin+jet.outward*jet.height:jet.tangentCoordinate;
    require(std::isfinite(radius)&&radius>=0,"RZ height jet has invalid interface radius");
    long double hoop=0;
    if(radius==0){
        require(jet.heightAxis==1&&jet.axisRegular&&jet.first==0,"RZ axis curvature requires an even regular axial height jet");
        hoop=-jet.second;
    }else hoop=jet.heightAxis==0?jet.outward/(metric*radius):-jet.first/(metric*radius);
    require(std::isfinite(meridional)&&std::isfinite(hoop),"unrepresentable RZ height jet curvature");
    return {meridional,hoop};
}
const std::vector<std::optional<HydroCsfHeightJet>>& HydroBalancedCsfLoad::heightFunctionCellJets() const
{
    require(data_->phase.model()==HydroPhaseGeometry::Model::GeometricVof&&
        data_->phase.vofPublication().geometryKind()==Interface::GeometricVofGeometryKind::LocalPlic&&
        data_->kind==HydroCsfDiscretization::GeometricHeightFunction,"load has no LocalPlic column jets");
    return data_->heightJets;
}
const std::vector<double>& HydroBalancedCsfLoad::heightFunctionCellCurvatures() const
{require(data_->kind==HydroCsfDiscretization::GeometricHeightFunction
    ||data_->kind==HydroCsfDiscretization::GeometricHeightFunctionIncrement,"load has no height-function cell curvature");return data_->heightCurvature;}
const std::vector<double>& HydroBalancedCsfLoad::heightFunctionCellEffectiveCurvatures() const
{
    require(data_->phase.model()==HydroPhaseGeometry::Model::GeometricVof&&
        data_->phase.vofPublication().geometryKind()==Interface::GeometricVofGeometryKind::LocalPlic&&
        data_->kind==HydroCsfDiscretization::GeometricHeightFunction,"load has no LocalPlic effective cell curvature");
    return data_->heightEffectiveCurvature;
}
const std::vector<HydroCsfCellCurvatureMethod>& HydroBalancedCsfLoad::heightFunctionCellCurvatureMethods() const
{
    require(data_->phase.model()==HydroPhaseGeometry::Model::GeometricVof&&
        data_->phase.vofPublication().geometryKind()==Interface::GeometricVofGeometryKind::LocalPlic&&
        data_->kind==HydroCsfDiscretization::GeometricHeightFunction,"load has no LocalPlic curvature methods");
    return data_->heightCurvatureMethods;
}
const std::vector<double>& HydroBalancedCsfLoad::heightFunctionIncrementReferenceLoad() const
{require(data_->kind==HydroCsfDiscretization::GeometricHeightFunctionIncrement,"load has no height-function increment reference");return data_->incrementReferenceLoad;}
const std::vector<double>& HydroBalancedCsfLoad::heightFunctionIncrementResidualLoad() const
{require(data_->kind==HydroCsfDiscretization::GeometricHeightFunctionIncrement,"load has no height-function increment residual");return data_->incrementResidualLoad;}
const std::vector<double>& HydroBalancedCsfLoad::heightFunctionIncrementMeanCurvatureDifference() const
{require(data_->kind==HydroCsfDiscretization::GeometricHeightFunctionIncrement,"load has no height-function increment mean");return data_->incrementDeltaMean;}
std::vector<std::pair<int,long double>> HydroBalancedCsfLoad::incrementNormalRateAtSample(const Interface::HermiteGalsSurfaceSample& sample) const
{
    require(data_->kind==HydroCsfDiscretization::GeometricHeightFunctionIncrement&&data_->incrementRates&&data_->incrementRateMutex,
        "load has no sampled increment generator (sigma-zero skips it)");
    const auto clock=data_->phase.timeLayer(),time=sample.timeLayer;
    require(sample.geometryIdentity==data_->phase.geometryIdentity()&&time.step==clock.step&&time.begin==clock.begin
        &&time.end==clock.end&&time.position==clock.position,"increment normal query has different geometry/full clock");
    require(data_->phase.layout().pressureDof(sample.frame.nativeCellIndex)>=0,"increment normal query has no retained native support");
    const std::lock_guard<std::mutex> lock(*data_->incrementRateMutex);
    const auto rows=data_->incrementRates->normal(sample.frame);
    return {rows.begin(),rows.end()};
}
const std::vector<HydroCsfFaceSupportReport>& HydroBalancedCsfLoad::supportReports() const{require(data_->capillary,"constant load has no curvature support reports");return data_->support;}
const HydroCsfCurvatureReport& HydroBalancedCsfLoad::curvatureReport() const{require(data_->capillary,"constant load has no curvature projection report");return data_->report;}
const std::vector<HydroCsfQuadraturePanel>& HydroBalancedCsfLoad::quadraturePanels(std::size_t velocity) const
{require(data_->capillary,"constant load has no curvature quadrature panels");return data_->panels.at(velocity);}
const std::vector<double>& HydroBalancedCsfLoad::colorGradient() const noexcept{return data_->gradient;}
const std::vector<double>& HydroBalancedCsfLoad::momentumRhsLoad() const noexcept{return data_->rhs;}
Interface::HermiteGalsIdentity HydroBalancedCsfLoad::geometryIdentity() const{return data_->phase.geometryIdentity();}
Interface::GeometricVofIdentity HydroBalancedCsfLoad::vofIdentity() const{return data_->phase.vofIdentity();}
void HydroBalancedCsfLoad::validateVofFor(const QuadtreeMesh& mesh,const HydroPressureLoad& pressure,
    Interface::GeometricVofIdentity geometry,Interface::GeometricVofTimeLayer clock) const
{
    require(&pressure.matrix()==&data_->pressure.matrix(),"different authoritative pressure owner");
    require(data_->phase.model()==HydroPhaseGeometry::Model::GeometricVof&&
        data_->kind==HydroCsfDiscretization::GeometricHeightFunction&&
        data_->phase.vofPublication().geometryKind()==Interface::GeometricVofGeometryKind::LocalPlic,
        "load is not a genuine LocalPlic height-function owner");
    data_->phase.validateVofFor(mesh,geometry,clock);
}
LAMP::TimeLayer HydroBalancedCsfLoad::timeLayer() const{return data_->phase.timeLayer();}
void HydroBalancedCsfLoad::validateFor(const QuadtreeMesh& mesh,const HydroPressureLoad& pressure,
    Interface::HermiteGalsIdentity geometry,LAMP::TimeLayer clock) const
{
    require(&pressure.matrix()==&data_->pressure.matrix(),"different authoritative pressure owner");
    data_->phase.validateFor(mesh,geometry,clock);
}
const HydroAdvectorReconstruction& HydroBalancedCsfLoad::reconstruction() const
{require(data_->source.has_value(),"load has no adjoint reconstruction");return *data_->source;}
const HydroAdvectorBoundary* HydroBalancedCsfLoad::boundary() const
{require(data_->source.has_value(),"load has no adjoint reconstruction");return data_->boundary?&*data_->boundary:nullptr;}
const std::vector<double>& HydroBalancedCsfLoad::referenceCurvatures() const
{require((data_->kind==HydroCsfDiscretization::GeometricAdjoint||data_->kind==HydroCsfDiscretization::GeometricProjectedAdjoint),"load has no adjoint reference curvature");return data_->referenceCurvature;}
bool HydroBalancedCsfLoad::hasExternalSurfaceVector() const noexcept{return data_->externalSurface.has_value();}
const Interface::HermiteGalsSurfaceVectorView& HydroBalancedCsfLoad::externalSurfaceVector() const
{require(data_->externalSurface.has_value(),"load has no authenticated external surface vector");return *data_->externalSurface;}
const std::vector<double>& HydroBalancedCsfLoad::referenceNormalStresses() const
{require((data_->kind==HydroCsfDiscretization::GeometricAdjoint||data_->kind==HydroCsfDiscretization::GeometricProjectedAdjoint),"load has no adjoint reference normal stress");return data_->referenceStress;}
const std::vector<double>& HydroBalancedCsfLoad::tangentialRhsLoad() const
{require((data_->kind==HydroCsfDiscretization::GeometricAdjoint||data_->kind==HydroCsfDiscretization::GeometricProjectedAdjoint),"load has no adjoint tangential RHS");return data_->tangentialLoad;}
const std::vector<double>& HydroBalancedCsfLoad::referencePressureLoad() const
{require((data_->kind==HydroCsfDiscretization::GeometricAdjoint||data_->kind==HydroCsfDiscretization::GeometricProjectedAdjoint),"load has no adjoint reference pressure");return data_->referenceLoad;}
const std::vector<double>& HydroBalancedCsfLoad::adjointResidualLoad() const
{require((data_->kind==HydroCsfDiscretization::GeometricAdjoint||data_->kind==HydroCsfDiscretization::GeometricProjectedAdjoint),"load has no adjoint residual traction");return data_->residualLoad;}
void HydroBalancedCsfLoad::validateFor(const QuadtreeMesh& mesh,const HydroPressureLoad& pressure,
    Interface::HermiteGalsIdentity geometry,LAMP::TimeLayer clock,
    const HydroAdvectorReconstruction& reconstruction,const HydroAdvectorBoundary* boundary) const
{
    validateFor(mesh,pressure,geometry,clock);
    require(data_->kind!=HydroCsfDiscretization::GeometricProjectedAdjoint
        &&data_->kind!=HydroCsfDiscretization::GeometricHeightFunctionIncrement,
        "projected load requires actual motion direction owner validation");
    require(data_->source.has_value(),"load has no adjoint reconstruction");
    require(data_->source->sameIdentity(reconstruction),"different adjoint reconstruction owner");
    require(data_->boundary.has_value()==(boundary!=nullptr),"different adjoint boundary presence");
    if(boundary)require(data_->boundary->sameIdentity(*boundary),"different adjoint boundary owner");
}
const HydroSpatialAdvector& HydroBalancedCsfLoad::projectedMotionAdvector() const
{require(data_->projectedMotion.has_value(),"load has no projected motion owner");return *data_->projectedMotion;}
const Interface::HermiteGalsAmrClosureSensitivity& HydroBalancedCsfLoad::mixedSensitivity() const
{require(data_->mixed.has_value(),"load has no projected mixed owner");return *data_->mixed;}
void HydroBalancedCsfLoad::validateFor(const QuadtreeMesh& mesh,const HydroPressureLoad& pressure,
    Interface::HermiteGalsIdentity geometry,LAMP::TimeLayer clock,
    const HydroSpatialAdvector& motion,const Interface::HermiteGalsAmrClosureSensitivity* mixed) const
{
    validateFor(mesh,pressure,geometry,clock);
    require(data_->projectedMotion.has_value(),"load has no projected motion owner");
    motion.validateFor(mesh,motion.momentumGeometry(),geometry,clock);
    require(data_->source->sameIdentity(motion.reconstruction()),"different adjoint reconstruction owner");
    require(data_->boundary.has_value()==(motion.boundary()!=nullptr),"different adjoint boundary presence");
    if(motion.boundary())require(data_->boundary->sameIdentity(*motion.boundary()),"different adjoint boundary owner");
    require(data_->projectedMotion->identity()==motion.identity(),"different projected motion direction owner");
    if(data_->sigma>0||data_->externalSurface){
        require(mixed&&data_->mixed.has_value(),"projected load lacks mixed closure owner");
        mixed->validateFor(mesh,geometry,clock);
        require(data_->mixed->symmetryPlanes()==mixed->symmetryPlanes(),"different projected mixed symmetry constraints");
    }
}
void HydroBalancedCsfLoad::validateFor(const QuadtreeMesh& mesh,const HydroPressureLoad& pressure,
    Interface::HermiteGalsIdentity geometry,LAMP::TimeLayer clock,
    const Interface::HermiteGalsSurfaceVectorView& external) const
{
    validateFor(mesh,pressure,geometry,clock);
    require(data_->externalSurface.has_value(),"load has no authenticated external surface vector");
    external.validateFor(mesh,data_->phase.publication(),clock);
    require(data_->externalSurface->sameIdentity(external),"different authoritative external surface vector owner");
}
void HydroBalancedCsfLoad::validateFor(const QuadtreeMesh& mesh,const HydroPressureLoad& pressure,const HydroRegularizedColor& color) const
{
    validateFor(mesh,pressure,color.geometryIdentity(),color.timeLayer());
    require(data_->color.has_value(),"geometric phase-fraction CSF load has no regularized color");
    require(&color.colors()==&data_->color->colors(),"different authoritative color owner/epsilon");
    data_->color->validateFor(mesh,color.phaseGeometry());
}
}
