#include "Core/Geometry/QuadraticCellCut.h"
#include "Core/Geometry/QuadraticCornerBoundary.h"
#include "Core/Geometry/Detail/QuadraticGraphExact.h"
#include <boost/multiprecision/cpp_int.hpp>
#include "Core/Geometry/Detail/ExactExpansion.h"
#include "Core/Geometry/Detail/ExactRationalRound.h"
#include "Core/Mesh/CartesianMeasure.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>
#include <stdexcept>
#include <sstream>

namespace LAMP::Geometry::Detail {
// This proof belongs to ONE immutable payload, not just equal native views.
// Copies/recentering retain no usable proof because owner identity differs.
// All quantities are the original actual graph's independent measurements.
struct QuadraticFixedPortCertificate {
    std::weak_ptr<const QuadraticGraphExact> owner;
    PlanarCellCut::Rectangle box;
    std::array<PlanarCellCut::Point,2> ports{};
    int phaseSign=0;bool axisymmetric=false;long double targetPhaseVolume=0;
    std::vector<GraphArc> arcs;
    QuadraticExact::Scalar total=0,below=0,above=0;
};
}

namespace LAMP::Geometry::QuadraticCellCut
{
namespace
{
using Poly=std::array<long double,5>;
void require(bool ok,const char* message) {if (!ok) throw std::invalid_argument(message);}
long double coordinate(Point p,int axis) {return axis?p.y:p.x;}
Point center(const Rectangle& b)
{return {static_cast<long double>(b.lower.x)+(static_cast<long double>(b.upper.x)-b.lower.x)/2,
    static_cast<long double>(b.lower.y)+(static_cast<long double>(b.upper.y)-b.lower.y)/2};}
void validate(Shape s)
{require((s.axis==0 || s.axis==1) && (s.sign==-1 || s.sign==1)
    && std::isfinite(s.slope) && std::isfinite(s.second),"[QuadraticCellCut] invalid graph shape");
    if(s.frameNormal) {const auto n=*s.frameNormal;const auto norm=std::hypot(n.x,n.y);
        require(std::isfinite(norm) && norm>0,"[QuadraticCellCut] invalid normal frame");}}
long double value(const Poly& p,long double x)
{long double result=p[4];for (int k=3;k>=0;--k) result=result*x+p[k];return result;}
void roots(long double a,long double b,long double c,long double lo,long double hi,std::vector<long double>& cuts)
{
    require(std::isfinite(a) && std::isfinite(b) && std::isfinite(c),"[QuadraticCellCut] nonfinite root polynomial");
    const auto add=[&](long double x){if (x>lo && x<hi) cuts.push_back(x);};
    if (a==0) {if (b!=0) add(-c/b);return;}
    const auto scale=std::max({std::abs(a),std::abs(b),std::abs(c)});a/=scale;b/=scale;c/=scale;
    const auto product=4*a*c;
    const auto discriminant=std::fma(b,b,-product)-std::fma(4*a,c,-product);
    if (discriminant<0) return;
    const auto q=-.5L*(b+std::copysign(std::sqrt(discriminant),b));
    if (q==0) add(-b/(2*a));else {add(q/a);add(c/q);}
}
// Integrate in a local interval coordinate, avoiding subtraction of large
// primitive values. The physical XY/RZ integrands have degree at most four.
long double integrate(const Poly& p,long double lo,long double hi)
{
    const auto middle=lo+(hi-lo)/2,half=(hi-lo)/2;
    const auto second=p[2]+3*p[3]*middle+6*p[4]*middle*middle;
    return 2*half*(value(p,middle)+second*half*half/3+p[4]*half*half*half*half/5);
}
Poly measure(const Poly& upper,const Poly& lower,bool rz,int axis,long double radius)
{
    Poly difference{};for (int k=0;k<5;++k) difference[k]=upper[k]-lower[k];
    if (!rz) return difference;
    Poly result{};
    if (axis==1) {
        for (int k=0;k<5;++k) {result[k]+=CartesianMeasure::twoPi*radius*difference[k];
            if (k+1<5) result[k+1]+=CartesianMeasure::twoPi*difference[k];}
    } else {
        Poly sum{};for (int k=0;k<5;++k) sum[k]=upper[k]+lower[k];sum[0]+=2*radius;
        for (int i=0;i<5;++i) for (int j=0;i+j<5;++j)
            result[i+j]+=(CartesianMeasure::twoPi/2)*difference[i]*sum[j];
    }
    return result;
}
GraphMeasure measureGraphPolygon(const Rectangle&,const std::vector<Point>&,Point polygonOrigin,
    Point origin,Vector2D normal,long double scale,const std::array<long double,3>&,bool rz,int sign,bool* full);
long double polygonVolume(const Rectangle& box,const std::vector<Point>& polygon,Point origin,
    const Reconstruction& reconstruction,bool rz,bool* entirelyLiquid=nullptr)
{
    if (entirelyLiquid) *entirelyLiquid=true;
    const auto s=reconstruction.shape;const auto c=center(box);
    if(s.frameNormal) {
        const long double scale=std::max(static_cast<long double>(box.upper.x)-box.lower.x,
            static_cast<long double>(box.upper.y)-box.lower.y);
        return measureGraphPolygon(box,polygon,origin,c,*s.frameNormal,scale,
            {{reconstruction.offset/scale,s.slope,s.second*scale/2}},rz,s.sign,entirelyLiquid).volume;
    }
    std::vector<Point> local;std::vector<long double> slabs;
    for (auto p:polygon) {p.x-=c.x-origin.x;p.y-=c.y-origin.y;
        const Point q{coordinate(p,1-s.axis),coordinate(p,s.axis)};local.push_back(q);slabs.push_back(q.x);}
    std::sort(slabs.begin(),slabs.end());slabs.erase(std::unique(slabs.begin(),slabs.end()),slabs.end());
    const Poly graph{{reconstruction.offset,s.slope,s.second/2,0,0}};long double volume=0;
    for (std::size_t k=1;k<slabs.size();++k) {
        const auto lo=slabs[k-1],hi=slabs[k],middle=lo+(hi-lo)/2;
        std::vector<Poly> edges;
        for (std::size_t e=0;e<local.size();++e) {
            const auto a=local[e],b=local[(e+1)%local.size()];
            if (a.x==b.x || middle<=std::min(a.x,b.x) || middle>=std::max(a.x,b.x)) continue;
            const auto slope=(b.y-a.y)/(b.x-a.x);edges.push_back({{a.y-slope*a.x,slope,0,0,0}});
        }
        if (edges.empty()) continue;
        require(edges.size()==2,"[QuadraticCellCut] invalid convex slice");
        if (value(edges[0],middle)>value(edges[1],middle)) std::swap(edges[0],edges[1]);
        std::vector<long double> cuts{lo,hi};
        for (const auto& edge:edges) roots(graph[2],graph[1]-edge[1],graph[0]-edge[0],lo,hi,cuts);
        std::sort(cuts.begin(),cuts.end());cuts.erase(std::unique(cuts.begin(),cuts.end()),cuts.end());
        for (std::size_t i=1;i<cuts.size();++i) {
            const auto a=cuts[i-1],b=cuts[i],mid=a+(b-a)/2;
            auto lower=edges[0],upper=edges[1];
            if (s.sign>0) {
                if (value(graph,mid)<=value(lower,mid)) {if (entirelyLiquid) *entirelyLiquid=false;continue;}
                if (value(graph,mid)<value(upper,mid)) {upper=graph;if (entirelyLiquid) *entirelyLiquid=false;}
            } else {
                if (value(graph,mid)>=value(upper,mid)) {if (entirelyLiquid) *entirelyLiquid=false;continue;}
                if (value(graph,mid)>value(lower,mid)) {lower=graph;if (entirelyLiquid) *entirelyLiquid=false;}
            }
            volume+=integrate(measure(upper,lower,rz,s.axis,c.x),a,b);
        }
    }
    require(std::isfinite(volume) && volume>=0,"[QuadraticCellCut] invalid curved volume");return volume;
}
}
Basis basis(const Shape& s)
{
    validate(s);
    if(!s.frameNormal)return s.axis?Basis{{1,0},{0,1}}:Basis{{0,1},{1,0}};
    const auto normal=*s.frameNormal;const auto length=std::hypot(static_cast<long double>(normal.x),static_cast<long double>(normal.y));
    const Point n{normal.x/length,normal.y/length};return {{-n.y,n.x},n};
}
namespace {
namespace QE=Detail::QuadraticExact;
using QR=QE::Scalar;using QV=QE::Vector;using QP=QE::Polynomial;
QV exactCenter(const Rectangle&b){return {(QE::scalar(b.lower.x)+QE::scalar(b.upper.x))/2,(QE::scalar(b.lower.y)+QE::scalar(b.upper.y))/2};}
std::shared_ptr<Detail::QuadraticGraphExact> lifted(const Rectangle&box,const Reconstruction&g) {
    Detail::validateGraph(g);if(g.exact)return std::make_shared<Detail::QuadraticGraphExact>(*g.exact);
    auto e=std::make_shared<Detail::QuadraticGraphExact>();const auto b=basis(g.shape);e->origin=exactCenter(box);
    e->tangent=QE::vector(b.tangent);e->normal=QE::vector(b.normal);e->determinant=QE::dot(e->tangent,e->tangent);
    // Legacy implicit uses raw dot coordinates. This conversion retains that
    // phase curve, including the native frame's non-unit determinant.
    e->q={QE::scalar(g.offset)/e->determinant,QE::scalar(g.shape.slope),QE::scalar(g.shape.second)*e->determinant/2};
    e->nativeShape=g.shape;e->nativeOffset=g.offset;return e;
}
void setViews(Reconstruction&g,const std::shared_ptr<Detail::QuadraticGraphExact>&e) {
    g.offset=QE::rounded(e->q[0]);g.shape.slope=QE::rounded(e->q[1]);g.shape.second=QE::rounded(2*e->q[2]);
    e->nativeOffset=g.offset;e->nativeShape=g.shape;g.exact=e;
}
struct ExactMeasure {QR total=0,below=0,above=0;std::array<QR,3> derivative{};};
ExactMeasure exactPolygon(const Detail::QuadraticGraphExact&g,const std::vector<QV>&physical,bool rz) {
    ExactMeasure out;if(physical.size()<3)return out;struct V{QR t,n;};std::vector<V> vertices;std::vector<QR> slabs;
    QR crossSum=0,rMoment=0;
    for(std::size_t i=0;i<physical.size();++i){const auto&a=physical[i];const auto&b=physical[(i+1)%physical.size()];const QR cross=a.x*b.y-b.x*a.y;crossSum+=cross;rMoment+=(a.x+b.x)*cross;
        const auto d=QE::subtract(a,g.origin);const QR t=QE::dot(d,g.tangent)/g.determinant,n=QE::dot(d,g.normal)/g.determinant;vertices.push_back({t,n});slabs.push_back(t);if(rz)require(a.x>=0,"[QuadraticCellCut] encoded polygon has negative radius");}
    if(crossSum<0){crossSum=-crossSum;rMoment=-rMoment;}out.total=rz?QE::scalar(CartesianMeasure::twoPi)*rMoment/6:crossSum/2;
    std::sort(slabs.begin(),slabs.end());slabs.erase(std::unique(slabs.begin(),slabs.end()),slabs.end());
    const QP graph{{g.q[0],g.q[1],g.q[2],0,0}};
    const auto density=[&](const QP&lo,const QP&hi){QP height{},r{},result{};for(int i=0;i<5;++i){height[i]=hi[i]-lo[i];r[i]=g.normal.x*(hi[i]+lo[i])/2;}r[0]+=g.origin.x;r[1]+=g.tangent.x;
        result=rz?QE::product(height,r):height;const QR metric=g.determinant*(rz?QE::scalar(CartesianMeasure::twoPi):QR(1));for(auto&x:result)x*=metric;return result;};
    for(std::size_t k=1;k<slabs.size();++k){const QR a=slabs[k-1],b=slabs[k],mid=(a+b)/2;std::vector<QP> edges;
        for(std::size_t i=0;i<vertices.size();++i){const auto&u=vertices[i];const auto&v=vertices[(i+1)%vertices.size()];if(u.t==v.t||mid<=std::min(u.t,v.t)||mid>=std::max(u.t,v.t))continue;
            const QR m=(v.n-u.n)/(v.t-u.t);edges.push_back({{u.n-m*u.t,m,0,0,0}});}
        if(edges.empty())continue;require(edges.size()==2,"[QuadraticCellCut] encoded polygon is not convex");if(QE::value(edges[0],mid)>QE::value(edges[1],mid))std::swap(edges[0],edges[1]);
        std::vector<QR> cuts{a,b};for(const auto&edge:edges){QP p{};for(int i=0;i<5;++i)p[i]=graph[i]-edge[i];const auto roots=QE::roots(p,a,b,g.portParameters);cuts.insert(cuts.end(),roots.begin(),roots.end());}
        std::sort(cuts.begin(),cuts.end());cuts.erase(std::unique(cuts.begin(),cuts.end()),cuts.end());
        for(std::size_t j=1;j<cuts.size();++j){const QR lo=cuts[j-1],hi=cuts[j],middle=(lo+hi)/2;const auto h=QE::value(graph,middle),lower=QE::value(edges[0],middle),upper=QE::value(edges[1],middle);
            if(h<=lower)out.above+=QE::integral(density(edges[0],edges[1]),lo,hi);
            else if(h>=upper)out.below+=QE::integral(density(edges[0],edges[1]),lo,hi);
            else {out.below+=QE::integral(density(edges[0],graph),lo,hi);out.above+=QE::integral(density(graph,edges[1]),lo,hi);
                QP weight{{1,0,0,0,0}};if(rz){weight={};for(int i=0;i<3;++i)weight[i]=g.normal.x*graph[i];weight[0]+=g.origin.x;weight[1]+=g.tangent.x;for(auto&x:weight)x*=QE::scalar(CartesianMeasure::twoPi);}
                for(int d=0;d<3;++d){QP monomial{};monomial[d]=g.determinant;out.derivative[d]+=QE::integral(QE::product(weight,monomial),lo,hi);}}
        }
    }
    require(out.below>=0&&out.above>=0&&out.below+out.above==out.total,"[QuadraticCellCut] encoded independent phase partition failed");return out;
}
std::vector<QV> exactCorners(const Rectangle&b){return {QE::vector({b.lower.x,b.lower.y}),QE::vector({b.upper.x,b.lower.y}),QE::vector({b.upper.x,b.upper.y}),QE::vector({b.lower.x,b.upper.y})};}
// A fixed-port chord is a LINE, so its two physical phases require only
// rational half-plane clipping and the polygon's exact area/radial moment.
// No quadratic edge roots or general graph slab decomposition are needed.
QR exactLinearPolygonVolume(const std::vector<QV>&polygon,bool rz){
    if(polygon.size()<3)return 0;QR crossSum=0,rMoment=0;
    for(std::size_t i=0;i<polygon.size();++i){const auto&a=polygon[i];const auto&b=polygon[(i+1)%polygon.size()];
        const QR cross=a.x*b.y-b.x*a.y;crossSum+=cross;rMoment+=(a.x+b.x)*cross;}
    if(crossSum<0){crossSum=-crossSum;rMoment=-rMoment;}
    return rz?QE::scalar(CartesianMeasure::twoPi)*rMoment/6:crossSum/2;
}
ExactMeasure exactFixedPortChord(const Detail::QuadraticGraphExact&chord,const Rectangle&box,bool rz){
    const auto polygon=exactCorners(box);if(rz)for(const auto&p:polygon)
        require(p.x>=0,"[QuadraticCellCut] encoded polygon has negative radius");
    const auto implicit=[&](const QV&p){const auto d=QE::subtract(p,chord.origin);
        return QR(QE::dot(d,chord.normal)-chord.q[1]*QE::dot(d,chord.tangent)-chord.determinant*chord.q[0]);};
    const auto phase=[&](bool below){std::vector<QV> clipped;
        for(std::size_t i=0;i<polygon.size();++i){const auto&a=polygon[i];const auto&b=polygon[(i+1)%polygon.size()];
            const QR fa=implicit(a),fb=implicit(b);const bool ia=below?fa<=0:fa>=0,ib=below?fb<=0:fb>=0;
            if(ia!=ib)clipped.push_back(QE::add(a,QE::multiply(QE::subtract(b,a),fa/(fa-fb))));
            if(ib)clipped.push_back(b);}
        return exactLinearPolygonVolume(clipped,rz);};
    ExactMeasure result;result.total=exactLinearPolygonVolume(polygon,rz);
    result.below=phase(true);result.above=phase(false);
    require(result.below>=0&&result.above>=0&&result.below+result.above==result.total,
        "[QuadraticCellCut] encoded independent phase partition failed");return result;
}
ExactMeasure exactFixedPortArc(const Detail::QuadraticGraphExact&graph,const Detail::QuadraticGraphExact&chord,
    const ExactMeasure&chordMeasure,const QR&lo,const QR&hi,bool rz){
    QP height{},radialAverage{},radialGraph{};
    for(int i=0;i<3;++i){height[i]=graph.q[i]-chord.q[i];
        radialAverage[i]=graph.normal.x*(graph.q[i]+chord.q[i])/2;radialGraph[i]=graph.normal.x*graph.q[i];}
    radialAverage[0]+=graph.origin.x;radialAverage[1]+=graph.tangent.x;
    radialGraph[0]+=graph.origin.x;radialGraph[1]+=graph.tangent.x;
    const QR metric=graph.determinant*(rz?QE::scalar(CartesianMeasure::twoPi):QR(1));
    const QP strip=rz?QE::product(height,radialAverage):height;
    const QR change=metric*QE::integral(strip,lo,hi);
    ExactMeasure result;result.total=chordMeasure.total;
    // Each phase starts from its independently clipped CHORD phase and adds
    // its own signed physical strip. Neither target nor a rounded dominant
    // complement supplies a phase measurement.
    result.below=chordMeasure.below+change;result.above=chordMeasure.above-change;
    for(int k=0;k<3;++k){QP monomial{};monomial[k]=metric;
        result.derivative[k]=QE::integral(rz?QE::product(radialGraph,monomial):monomial,lo,hi);}
    require(result.below>=0&&result.above>=0&&result.below+result.above==result.total,
        "[QuadraticCellCut] encoded independent phase partition failed");return result;
}
StripPhaseMeasure exactStrip(const Rectangle&box,const Reconstruction&g,int normal,Point a,Point b,long double depth,bool positive,bool rz){
    Detail::validateGraph(g);require((normal==0||normal==1)&&std::isfinite(depth)&&depth>0&&coordinate(a,normal)==coordinate(b,normal),"[QuadraticCellCut] invalid encoded strip");
    // positive is the outward transport direction: its donor strip is
    // [near-depth,near], exactly as the original axis/framed kernels.
    const auto d=QE::scalar(depth)*(positive?-1:1);const QV va=QE::vector(a),vb=QE::vector(b),shift=normal?QV{0,d}:QV{d,0};
    const std::vector<QV> polygon{va,vb,QE::add(vb,shift),QE::add(va,shift)};
    for(const auto&p:polygon)require(p.x>=QE::scalar(box.lower.x)&&p.x<=QE::scalar(box.upper.x)&&p.y>=QE::scalar(box.lower.y)&&p.y<=QE::scalar(box.upper.y),"[QuadraticCellCut] encoded strip exceeds donor");
    const auto m=exactPolygon(*g.exact,polygon,rz);return {QE::rounded(m.total),QE::rounded(g.shape.sign>0?m.below:m.above),QE::rounded(g.shape.sign>0?m.above:m.below)};
}
// Retain face, displacement and query coordinates independently until this
// same exact-polygon measure. In particular, a thin selected triangle may have
// both its world zero and its back-edge below native coordinate resolution.
StripPhaseMeasure exactLinearSweep(const Rectangle&box,const Reconstruction&g,int normal,
    Point a,Point b,const QR& firstDepth,const QR& lastDepth,bool positive,bool rz,
    const std::optional<std::array<Point,2>>& query,bool identifyCoveredPhases)
{
    Detail::validateGraph(g);require(g.exact!=nullptr,"[QuadraticCellCut/LinearSweep] exact graph required");
    validate(g.shape);
    require((normal==0||normal==1)&&std::isfinite(a.x)&&std::isfinite(a.y)&&std::isfinite(b.x)&&std::isfinite(b.y)
        &&coordinate(a,normal)==coordinate(b,normal)&&coordinate(a,1-normal)!=coordinate(b,1-normal),
        "[QuadraticCellCut/LinearSweep] invalid original near-face segment");
    require(std::isfinite(box.lower.x)&&std::isfinite(box.lower.y)&&std::isfinite(box.upper.x)&&std::isfinite(box.upper.y)
        &&box.lower.x<box.upper.x&&box.lower.y<box.upper.y&&(!rz||box.lower.x>=0),
        "[QuadraticCellCut/LinearSweep] invalid donor rectangle");
    const QV va=QE::vector(a),vb=QE::vector(b),edge=QE::subtract(vb,va);
    const auto inside=[&](const QV&p){return p.x>=QE::scalar(box.lower.x)&&p.x<=QE::scalar(box.upper.x)
        &&p.y>=QE::scalar(box.lower.y)&&p.y<=QE::scalar(box.upper.y);};
    require(inside(va)&&inside(vb),"[QuadraticCellCut/LinearSweep] original near segment exceeds donor");
    QR lo=0,hi=1;
    if(query){std::array<QR,2> parameters{};
        for(int k=0;k<2;++k){const auto p=(*query)[k];require(std::isfinite(p.x)&&std::isfinite(p.y)
                &&coordinate(p,normal)==coordinate(a,normal),"[QuadraticCellCut/LinearSweep] query leaves original face plane");
            parameters[k]=(QE::scalar(coordinate(p,1-normal))-QE::scalar(coordinate(a,1-normal)))
                /(QE::scalar(coordinate(b,1-normal))-QE::scalar(coordinate(a,1-normal)));
            require(parameters[k]>=0&&parameters[k]<=1,"[QuadraticCellCut/LinearSweep] query exceeds original segment");}
        lo=std::min(parameters[0],parameters[1]);hi=std::max(parameters[0],parameters[1]);
        require(lo<hi,"[QuadraticCellCut/LinearSweep] zero query interval");}
    const QR rate=lastDepth-firstDepth,direction=positive?QR(1):QR(-1);
    const QR dl=direction*(firstDepth+rate*lo),dh=direction*(firstDepth+rate*hi);
    if(dl<=0&&dh<=0)return {};
    // The unique affine zero is retained as a ratio of ORIGINAL binary input
    // values; never construct a rounded world splitting point in the caller.
    if(dl<0)lo=-firstDepth/rate;if(dh<0)hi=-firstDepth/rate;
    require(lo<hi,"[QuadraticCellCut/LinearSweep] invalid selected direction interval");
    const QV nearLo=QE::add(va,QE::multiply(edge,lo)),nearHi=QE::add(va,QE::multiply(edge,hi));
    const QR depthLo=firstDepth+rate*lo,depthHi=firstDepth+rate*hi;
    const QV shiftLo=normal?QV{0,-depthLo}:QV{-depthLo,0};
    const QV shiftHi=normal?QV{0,-depthHi}:QV{-depthHi,0};
    std::vector<QV> polygon;
    for(const auto&p:{nearLo,nearHi,QE::add(nearHi,shiftHi),QE::add(nearLo,shiftLo)}){
        require(inside(p),"[QuadraticCellCut/LinearSweep] selected wedge exceeds donor");
        if(polygon.empty()||p.x!=polygon.back().x||p.y!=polygon.back().y)polygon.push_back(p);}
    if(polygon.size()>1&&polygon.front().x==polygon.back().x&&polygon.front().y==polygon.back().y)polygon.pop_back();
    require(polygon.size()>=3,"[QuadraticCellCut/LinearSweep] selected wedge collapsed");
    const auto m=exactPolygon(*g.exact,polygon,rz);
    require(m.total>0,"[QuadraticCellCut/LinearSweep] nonpositive selected wedge");
    const QR&liquid=g.shape.sign>0?m.below:m.above;const QR&gas=g.shape.sign>0?m.above:m.below;
    StripPhaseMeasure result{QE::rounded(m.total),QE::rounded(liquid),QE::rounded(gas)};
    if(identifyCoveredPhases){const auto donor=exactPolygon(*g.exact,exactCorners(box),rz);
        const QR&donorLiquid=g.shape.sign>0?donor.below:donor.above;const QR&donorGas=g.shape.sign>0?donor.above:donor.below;
        if(donorLiquid>0&&liquid==donorLiquid)result.coveredPhases|=1;
        if(donorGas>0&&gas==donorGas)result.coveredPhases|=2;}
    return result;
}
// Exact binary-input second-order algebra for one fixed-face path. Values
// are reduced rational numbers; native rounding occurs only at the API edge.
struct BoundaryRJet {
    QR value=0,first=0,second=0;
    BoundaryRJet()=default;
    BoundaryRJet(long long v):value(v){}
    BoundaryRJet(const QR&v,const QR&d=QR(0),const QR&dd=QR(0)):value(v),first(d),second(dd){}
};
BoundaryRJet operator+(const BoundaryRJet&a,const BoundaryRJet&b){return {a.value+b.value,a.first+b.first,a.second+b.second};}
BoundaryRJet operator-(const BoundaryRJet&a,const BoundaryRJet&b){return {a.value-b.value,a.first-b.first,a.second-b.second};}
BoundaryRJet operator-(const BoundaryRJet&a){return {-a.value,-a.first,-a.second};}
BoundaryRJet operator*(const BoundaryRJet&a,const BoundaryRJet&b){return {a.value*b.value,a.first*b.value+a.value*b.first,a.second*b.value+2*a.first*b.first+a.value*b.second};}
BoundaryRJet operator/(const BoundaryRJet&a,const BoundaryRJet&b){
    require(b.value!=0,"[QuadraticCellCut/BoundaryJet] singular exact quotient");
    const BoundaryRJet inverse{QR(1)/b.value,-b.first/(b.value*b.value),2*b.first*b.first/(b.value*b.value*b.value)-b.second/(b.value*b.value)};
    return a*inverse;
}
using BoundaryJPoly=std::array<BoundaryRJet,5>;
BoundaryJPoly boundaryProduct(const BoundaryJPoly&a,const BoundaryJPoly&b){BoundaryJPoly out{};for(int i=0;i<5;++i)for(int j=0;i+j<5;++j)out[i+j]=out[i+j]+a[i]*b[j];return out;}
BoundaryRJet boundaryIntegral(const BoundaryJPoly&p,const BoundaryRJet&lo,const BoundaryRJet&hi){BoundaryRJet a=1,b=1,out=0;for(int i=0;i<5;++i){a=a*lo;b=b*hi;out=out+p[i]*(b-a)/BoundaryRJet(i+1);}return out;}
struct BoundaryJPoint {BoundaryRJet x,y;};
BoundaryRJet boundaryPolygonVolume(const std::vector<BoundaryJPoint>&polygon,bool rz){
    BoundaryRJet crossSum=0,radial=0;if(polygon.size()<3)return {};
    for(std::size_t i=0;i<polygon.size();++i){const auto&a=polygon[i];const auto&b=polygon[(i+1)%polygon.size()];const auto cross=a.x*b.y-b.x*a.y;crossSum=crossSum+cross;radial=radial+(a.x+b.x)*cross;}
    require(crossSum.value>=0,"[QuadraticCellCut/BoundaryJet] chord polygon orientation changed");
    return rz?BoundaryRJet(QE::scalar(CartesianMeasure::twoPi))*radial/6:crossSum/2;
}
Jet2 nativeBoundaryJet(const BoundaryRJet&j){return {QE::rounded(j.value),QE::rounded(j.first),QE::rounded(j.second)};}
struct ExactStationaryFaceBoundary {
    QR parameter=0;QV position{};bool inOriginalArc=false;
    std::array<int,2> branchSign{};std::array<QR,2> slack{};
    std::array<std::array<QR,3>,2> gradient{};
};
std::array<std::optional<ExactStationaryFaceBoundary>,2> exactStationaryFaceBoundaries(const Rectangle&box,
    const Detail::QuadraticGraphExact&g,const std::vector<Detail::GraphArc>&arcs){
    require(!arcs.empty(),"[QuadraticCellCut/StationaryFace] current exact graph has no actual arc");
    const QR h=std::max(QE::scalar(box.upper.x)-QE::scalar(box.lower.x),QE::scalar(box.upper.y)-QE::scalar(box.lower.y));
    require(h>0&&g.determinant>0,"[QuadraticCellCut/StationaryFace] invalid current cell/frame");
    std::array<std::optional<ExactStationaryFaceBoundary>,2> result;
    for(int axis=0;axis<2;++axis){const QR normal=axis?g.normal.y:g.normal.x,tangent=axis?g.tangent.y:g.tangent.x,a=normal*g.q[2];if(a==0)continue;
        ExactStationaryFaceBoundary out;out.parameter=-(tangent+normal*g.q[1])/(2*a);out.position=Detail::graphPosition(g,out.parameter);
        const QR other=axis?out.position.x:out.position.y,otherLo=QE::scalar(axis?box.lower.x:box.lower.y),otherHi=QE::scalar(axis?box.upper.x:box.upper.y);
        if(other<otherLo||other>otherHi)continue; // Actual physical face tangent interval, exact comparison.
        out.inOriginalArc=std::any_of(arcs.begin(),arcs.end(),[&](const Detail::GraphArc&arc){
            return out.parameter>=arc.first.parameter&&out.parameter<=arc.last.parameter;});
        const QR value=axis?out.position.y:out.position.x,lower=QE::scalar(axis?box.lower.y:box.lower.x),upper=QE::scalar(axis?box.upper.y:box.upper.x),scaledT=out.parameter/h;
        const std::array<QR,3> derivative{{normal,normal*scaledT,normal*scaledT*scaledT}};
        for(int side=0;side<2;++side){const QR inward=side?upper-value:value-lower;
            out.branchSign[side]=inward>0?1:inward<0?-1:out.inOriginalArc?1:-1;
            out.slack[side]=out.branchSign[side]*inward/h;
            for(int k=0;k<3;++k)out.gradient[side][k]=out.branchSign[side]*(side?-1:1)*derivative[k];}
        result[axis]=std::move(out);
    }return result;
}
StationaryFaceBoundary nativeStationaryFaceBoundary(const ExactStationaryFaceBoundary&value,int axis){
    StationaryFaceBoundary result;result.parameter=QE::rounded(value.parameter);result.position=QE::native(value.position);result.inOriginalArc=value.inOriginalArc;
    for(int side=0;side<2;++side){result.side[side]=2*axis+side;result.branchSign[side]=value.branchSign[side];result.slack[side]=QE::rounded(value.slack[side]);
        for(int k=0;k<3;++k)result.normalizedCoeffGradient[side][k]=QE::rounded(value.gradient[side][k]);}return result;
}
}
Reconstruction fromFixedPorts(const Rectangle&box,Shape shape,Point first,Point last,long double eta){
    validate(shape);require(std::isfinite(eta),"[QuadraticCellCut] invalid fixed-port eta");const auto b=basis(shape);auto e=std::make_shared<Detail::QuadraticGraphExact>();
    e->origin=exactCenter(box);e->tangent=QE::vector(b.tangent);e->normal=QE::vector(b.normal);e->determinant=QE::dot(e->tangent,e->tangent);
    const auto u=QE::subtract(QE::vector(first),e->origin),v=QE::subtract(QE::vector(last),e->origin);const QR t0=QE::dot(u,e->tangent)/e->determinant,t1=QE::dot(v,e->tangent)/e->determinant,n0=QE::dot(u,e->normal)/e->determinant,n1=QE::dot(v,e->normal)/e->determinant;
    require(t0!=t1,"[QuadraticCellCut] fixed ports have equal tangent parameter");const QR h=std::max(QE::scalar(box.upper.x)-QE::scalar(box.lower.x),QE::scalar(box.upper.y)-QE::scalar(box.lower.y));require(h>0,"[QuadraticCellCut] invalid fixed-port box");
    const QR c2=QE::scalar(eta)/h,m=(n1-n0)/(t1-t0);e->q={n0-m*t0+c2*t0*t1,m-c2*(t0+t1),c2};e->portParameters={t0,t1};e->ports={first,last};
    Reconstruction result;result.shape=shape;setViews(result,e);return result;
}
ThreePortGraphResult fromThreePorts(const Rectangle&box,Shape shape,const std::array<Point,3>&ports){
    validate(shape);(void)PlanarCellCut::capacity(box,false);const auto b=basis(shape);
    auto e=std::make_shared<Detail::QuadraticGraphExact>();e->origin=exactCenter(box);
    e->tangent=QE::vector(b.tangent);e->normal=QE::vector(b.normal);e->determinant=QE::dot(e->tangent,e->tangent);
    const QR x0=QE::scalar(box.lower.x),x1=QE::scalar(box.upper.x),y0=QE::scalar(box.lower.y),y1=QE::scalar(box.upper.y);
    const QR h=std::max(x1-x0,y1-y0);require(h>0&&e->determinant>0,"[QuadraticCellCut/ThreePorts] invalid cell/frame");
    std::array<QR,3> t{},n{};std::array<std::array<QR,3>,3> lagrange{};
    for(int i=0;i<3;++i){
        require(std::isfinite(ports[i].x)&&std::isfinite(ports[i].y),"[QuadraticCellCut/ThreePorts] nonfinite physical port");
        const auto p=QE::vector(ports[i]);require(p.x>=x0&&p.x<=x1&&p.y>=y0&&p.y<=y1,"[QuadraticCellCut/ThreePorts] port outside actual cell");
        unsigned face=0;if(p.x==x0)face|=1;if(p.x==x1)face|=2;if(p.y==y0)face|=4;if(p.y==y1)face|=8;
        require(face!=0&&(face&(face-1))==0,"[QuadraticCellCut/ThreePorts] actual open face required; corner/interior unsupported");
        const auto d=QE::subtract(p,e->origin);t[i]=QE::dot(d,e->tangent)/e->determinant;n[i]=QE::dot(d,e->normal)/e->determinant;
    }
    for(int i=0;i<3;++i){const int j=(i+1)%3,k=(i+2)%3;const QR denominator=(t[i]-t[j])*(t[i]-t[k]);
        require(denominator!=0,"[QuadraticCellCut/ThreePorts] equal tangent parameters; quadratic chart is not unisolvent");
        lagrange[i]={{t[j]*t[k]/denominator,-(t[j]+t[k])/denominator,QR(1)/denominator}};
        for(int coefficient=0;coefficient<3;++coefficient)e->q[coefficient]+=n[i]*lagrange[i][coefficient];
    }
    e->portParameters.assign(t.begin(),t.end());e->ports.assign(ports.begin(),ports.end());
    ThreePortGraphResult result;result.graph.shape=shape;setViews(result.graph,e);
    const std::array<QR,3> normalize{{QR(1)/h,1,h}};
    for(int i=0;i<3;++i){const QR slope=e->q[1]+2*e->q[2]*t[i];
        for(int coordinate=0;coordinate<2;++coordinate){
            const QR tangent=coordinate?e->tangent.y:e->tangent.x,normal=coordinate?e->normal.y:e->normal.x;
            const QR heightDerivative=(normal-slope*tangent)/e->determinant;
            for(int coefficient=0;coefficient<3;++coefficient)
                result.normalizedCoefficientJacobian[i][coordinate][coefficient]=QE::rounded(lagrange[i][coefficient]*heightDerivative*normalize[coefficient]);
        }
    }
    return result;
}
SameFaceDividedIncidence sameFaceDividedIncidence(const Rectangle&box,const Reconstruction&current,
    Point knownRoot,Point fourthNode){
    require(current.exact!=nullptr,"[QuadraticCellCut/DividedIncidence] current exact graph required; legacy unsupported");
    Detail::validateGraph(current);(void)basis(current.shape);(void)PlanarCellCut::capacity(box,false);
    const auto&g=*current.exact;const QR x0=QE::scalar(box.lower.x),x1=QE::scalar(box.upper.x),y0=QE::scalar(box.lower.y),y1=QE::scalar(box.upper.y);
    const QR h=std::max(x1-x0,y1-y0);require(h>0&&g.determinant>0,"[QuadraticCellCut/DividedIncidence] invalid cell/frame");
    const std::array<Point,2> points{{knownRoot,fourthNode}};std::array<unsigned,2> faces{};std::array<QR,2> t{};
    QR knownHeight;
    for(int i=0;i<2;++i){
        require(std::isfinite(points[i].x)&&std::isfinite(points[i].y),"[QuadraticCellCut/DividedIncidence] nonfinite physical point");
        const auto p=QE::vector(points[i]);require(p.x>=x0&&p.x<=x1&&p.y>=y0&&p.y<=y1,"[QuadraticCellCut/DividedIncidence] point outside actual cell");
        if(p.x==x0)faces[i]|=1;if(p.x==x1)faces[i]|=2;if(p.y==y0)faces[i]|=4;if(p.y==y1)faces[i]|=8;
        require(faces[i]!=0&&(faces[i]&(faces[i]-1))==0,"[QuadraticCellCut/DividedIncidence] actual open face required; corner/interior unsupported");
        const auto delta=QE::subtract(p,g.origin);t[i]=QE::dot(delta,g.tangent)/g.determinant;
        if(i==0)knownHeight=QE::dot(delta,g.normal)/g.determinant;
    }
    require(faces[0]==faces[1],"[QuadraticCellCut/DividedIncidence] points must share the same actual open face");
    require(knownHeight==(g.q[2]*t[0]+g.q[1])*t[0]+g.q[0],"[QuadraticCellCut/DividedIncidence] known point is not an exact root of the current curve");
    const int axis=(faces[0]&3)?0:1;const QR normal=axis?g.normal.y:g.normal.x,tangent=axis?g.tangent.y:g.tangent.x;
    const QR otherTangent=axis?g.tangent.x:g.tangent.y,a=normal*g.q[2],b=tangent+normal*g.q[1];
    require(a!=0,"[QuadraticCellCut/DividedIncidence] side polynomial is not quadratic");
    require(otherTangent!=0,"[QuadraticCellCut/DividedIncidence] zero projection onto physical face tangent");
    const QR companionParameter=-b/a-t[0];const auto companion=Detail::graphPosition(g,companionParameter);
    const QR factor=g.determinant/(otherTangent*h);
    SameFaceDividedIncidence result;result.residual=QE::rounded(factor*(t[0]+t[1]+b/a));result.companion=QE::native(companion);
    result.face=faces[0]==1?0:faces[0]==2?1:faces[0]==4?2:3;
    result.normalizedCoefficientJacobian={{0,QE::rounded(g.determinant/(otherTangent*h*g.q[2])),
        QE::rounded(-g.determinant*b/(otherTangent*h*h*normal*g.q[2]*g.q[2]))}};
    for(int point=0;point<2;++point)for(int coordinate=0;coordinate<2;++coordinate)
        result.pointJacobian[point][coordinate]=QE::rounded((coordinate?g.tangent.y:g.tangent.x)/(otherTangent*h));
    return result;
}
FixedPortVolumeResult fixedPortVolumeGraph(const Rectangle&box,Shape shape,Point first,Point last,
    long double targetPhaseVolume,bool rz,long double branchEta){
    validate(shape);require(std::isfinite(targetPhaseVolume)&&std::isfinite(branchEta),"[QuadraticCellCut] nonfinite fixed-port volume request");
    if(rz)require(first.x>0&&last.x>0,"[QuadraticCellCut] regular physical-axis ports require the separate axis volume chart");
    const auto onBoundary=[&](Point p){const auto x=QE::scalar(p.x),y=QE::scalar(p.y);return x>=QE::scalar(box.lower.x)&&x<=QE::scalar(box.upper.x)&&y>=QE::scalar(box.lower.y)&&y<=QE::scalar(box.upper.y)&&
        (x==QE::scalar(box.lower.x)||x==QE::scalar(box.upper.x)||y==QE::scalar(box.lower.y)||y==QE::scalar(box.upper.y));};
    require(onBoundary(first)&&onBoundary(last),"[QuadraticCellCut] fixed-volume ports are not actual cell-boundary events");
    const auto openPort=[&](Point p){unsigned mask=0;if(p.x==box.lower.x)mask|=1;if(p.x==box.upper.x)mask|=2;
        if(p.y==box.lower.y)mask|=4;if(p.y==box.upper.y)mask|=8;return mask!=0&&(mask&(mask-1))==0;};
    const bool ordinaryPorts=openPort(first)&&openPort(last);
    auto chord=fromFixedPorts(box,shape,first,last,0);const auto&e=*chord.exact;
    const auto plain=ordinaryPorts?exactFixedPortChord(e,box,rz):exactPolygon(e,exactCorners(box),rz);
    const QR target=QE::scalar(targetPhaseVolume),capacity=plain.total;
    require(target>0&&target<capacity,"[QuadraticCellCut] fixed-port volume requires actual mixed target");
    const QR h=std::max(QE::scalar(box.upper.x)-QE::scalar(box.lower.x),QE::scalar(box.upper.y)-QE::scalar(box.lower.y));
    const QR t0=e.portParameters[0],t1=e.portParameters[1],lo=std::min(t0,t1),hi=std::max(t0,t1);
    const QP bubble{{t0*t1/h,-(t0+t1)/h,QR(1)/h,0,0}};
    const QR v0=shape.sign>0?plain.below:plain.above;QR a=0,b=0;const QR metric=shape.sign*e.determinant;
    if(!rz)b=metric*QE::integral(bubble,lo,hi);
    else {
        QP radius{{e.origin.x+e.normal.x*e.q[0],e.tangent.x+e.normal.x*e.q[1],0,0,0}};
        const QR revolved=metric*QE::scalar(CartesianMeasure::twoPi);
        b=revolved*QE::integral(QE::product(bubble,radius),lo,hi);
        a=revolved*e.normal.x*QE::integral(QE::product(bubble,bubble),lo,hi)/2;
    }
    const QR c=v0-target;std::vector<QR> roots;FixedPortVolumeResult result;result.quadratic=a!=0;
    if(a==0){require(b!=0,"[QuadraticCellCut] singular fixed-port linear volume equation");roots.push_back(-c/b);}
    else {
        const QR discriminant=b*b-4*a*c;require(discriminant>0,"[QuadraticCellCut] fixed-port volume has no regular real branch");
        const QR radical=QE::radical(discriminant),q=-(b+(b>=0?radical:-radical))/2;
        require(q!=0,"[QuadraticCellCut] fixed-port quadratic has no regular q root");roots={q/a,c/q};
    }
    result.realRoots=static_cast<unsigned>(roots.size());const QR branchDerivative=2*a*QE::scalar(branchEta)+b;
    require(branchDerivative!=0,"[QuadraticCellCut] declared fixed-port volume branch is its vertex");QR eta=0;
    for(const auto&r:roots){const QR derivative=2*a*r+b;if(derivative!=0&&((derivative>0)==(branchDerivative>0))){eta=r;++result.branchCandidates;}}
    require(result.branchCandidates==1,"[QuadraticCellCut] fixed-port volume branch has no unique root");
    auto solved=std::make_shared<Detail::QuadraticGraphExact>(e);
    solved->q={e.q[0]+eta*bubble[0],e.q[1]+eta*bubble[1],eta*bubble[2]};result.graph=chord;setViews(result.graph,solved);
    const auto arcs=Detail::graphArcs(*solved,box);
    if(!(arcs.size()==1&&arcs[0].first.parameter==lo&&arcs[0].last.parameter==hi)){
        // Error-path observability only. All quantities below belong to this
        // SAME rejected exact polynomial, never an alternate candidate/root.
        std::ostringstream why;why.precision(std::numeric_limits<long double>::max_digits10);
        why<<"[QuadraticCellCut] analytic fixed-port volume leaves its original one-arc branch"
            <<" reason="<<(arcs.empty()?"no-arcs":arcs.size()>1?"multiple-arcs":"parameter-mismatch")
            <<" arcCount="<<arcs.size()<<" expectedLo="<<QE::rounded(lo)<<" expectedHi="<<QE::rounded(hi)
            <<" p0=("<<first.x<<','<<first.y<<") p1=("<<last.x<<','<<last.y<<")"
            <<" box=("<<box.lower.x<<','<<box.lower.y<<";"<<box.upper.x<<','<<box.upper.y<<")"
            <<" rz="<<rz<<" phaseSign="<<shape.sign<<" target="<<targetPhaseVolume
            <<" eta="<<QE::rounded(eta)<<" q=("<<QE::rounded(solved->q[0])<<','<<QE::rounded(solved->q[1])<<','<<QE::rounded(solved->q[2])<<")"
            <<" O=("<<QE::rounded(solved->origin.x)<<','<<QE::rounded(solved->origin.y)<<")"
            <<" T=("<<QE::rounded(solved->tangent.x)<<','<<QE::rounded(solved->tangent.y)<<")"
            <<" N=("<<QE::rounded(solved->normal.x)<<','<<QE::rounded(solved->normal.y)<<") D="<<QE::rounded(solved->determinant);
        for(std::size_t i=0;i<arcs.size();++i){const auto&event=arcs[i];
            why<<" arc["<<i<<"]=("<<QE::rounded(event.first.parameter)<<','<<QE::rounded(event.last.parameter)
                <<";begin="<<event.first.position.x<<','<<event.first.position.y<<";end="<<event.last.position.x<<','<<event.last.position.y
                <<";masks="<<event.first.boundary<<','<<event.last.boundary
                <<";loGap="<<QE::rounded(event.first.parameter-lo)<<";hiGap="<<QE::rounded(event.last.parameter-hi)<<')';}
        const auto corners=exactCorners(box);QR boxLo=0,boxHi=0;bool initialized=false;
        for(std::size_t i=0;i<corners.size();++i){const auto d=QE::subtract(corners[i],solved->origin);const QR t=QE::dot(d,solved->tangent)/solved->determinant,n=QE::dot(d,solved->normal)/solved->determinant;
            const QR implicit=QR(shape.sign)*(n-(solved->q[2]*t+solved->q[1])*t-solved->q[0])/h;
            why<<" corner["<<i<<"]="<<QE::rounded(implicit);
            if(!initialized){boxLo=boxHi=t;initialized=true;}else{boxLo=std::min(boxLo,t);boxHi=std::max(boxHi,t);}}
        const std::array<const char*,4> names{{"Left","Right","Bottom","Top"}};
        for(int side=0;side<4;++side){const int axis=side/2;const long double border=axis?(side%2?box.upper.y:box.lower.y):(side%2?box.upper.x:box.lower.x);
            const auto equation=Detail::coordinatePolynomial(*solved,axis,border);const QR disc=equation[1]*equation[1]-4*equation[2]*equation[0];
            why<<" side["<<names[side]<<"]=(c="<<QE::rounded(equation[0])<<",b="<<QE::rounded(equation[1])<<",a="<<QE::rounded(equation[2])
                <<",disc="<<QE::rounded(disc)<<",coincident="<<(equation[0]==0&&equation[1]==0&&equation[2]==0);
            const auto recordPoint=[&](const char*label,const QR&t){const auto p=Detail::graphPosition(*solved,t);const QR gap=QE::value(equation,t),insideGap=(side%2)?-gap:gap;
                const QR other=axis?p.x:p.y,otherLo=QE::scalar(axis?box.lower.x:box.lower.y),otherHi=QE::scalar(axis?box.upper.x:box.upper.y);
                why<<';'<<label<<"=(t="<<QE::rounded(t)<<",x="<<QE::rounded(p.x)<<",y="<<QE::rounded(p.y)
                    <<",sideGap="<<QE::rounded(gap)<<",insideGap="<<QE::rounded(insideGap)
                    <<",inExpected="<<(t>=lo&&t<=hi)<<",inFaceTangent="<<(other>=otherLo&&other<=otherHi)<<')';};
            if(equation[2]!=0)recordPoint("stationary",-equation[1]/(2*equation[2]));else why<<";linear=true";
            const auto roots=QE::roots(equation,boxLo,boxHi,solved->portParameters);why<<";roots="<<roots.size();for(const auto&t:roots)recordPoint("root",t);why<<')';
        }
        throw std::invalid_argument(why.str());
    }
    // The already computed exact arcs are the SAME side events consumed by
    // S::evaluate. Reuse their endpoint certificate without integrating area.
    const auto&arc=arcs.front();const Point expectedFirst=t0<t1?first:last,expectedLast=t0<t1?last:first;
    require(arc.first.position.x==expectedFirst.x&&arc.first.position.y==expectedFirst.y&&
        arc.last.position.x==expectedLast.x&&arc.last.position.y==expectedLast.y&&
        arc.first.boundary!=0&&arc.last.boundary!=0,
        "[QuadraticCellCut] analytic volume changed original surface endpoint events");
    // Positive exact width and D certify positive ds; retain the original
    // native event-interval representability gate without a quadrature call.
    require(solved->determinant>0&&hi>lo&&QE::rounded(hi)>QE::rounded(lo),
        "[QuadraticCellCut] analytic volume has no representable regular arc");
    if(rz){
        const QP radius{{solved->origin.x+solved->normal.x*solved->q[0],
            solved->tangent.x+solved->normal.x*solved->q[1],solved->normal.x*solved->q[2],0,0}};
        QR least=std::min(QE::value(radius,lo),QE::value(radius,hi));
        if(radius[2]>0){const QR vertex=-radius[1]/(2*radius[2]);if(vertex>lo&&vertex<hi)least=std::min(least,QE::value(radius,vertex));}
        require(least>=0,"[QuadraticCellCut] analytic volume arc crosses negative physical radius");
        // The ordinary chart already requires both original ports r>0.
        // Continuity and positive ds therefore certify positive revolved area.
        require(expectedFirst.x>0&&expectedLast.x>0,
            "[QuadraticCellCut] analytic volume has no positive ordinary-chart RZ area");
    }
    const auto simplePort=[&](const Detail::GraphEvent&event){const unsigned mask=event.boundary;
        if(mask==0||(mask&(mask-1))!=0)return false;const int axis=(mask&3)?0:1;
        const QR slope=solved->q[1]+2*solved->q[2]*event.parameter;
        return (axis?solved->tangent.y:solved->tangent.x)+(axis?solved->normal.y:solved->normal.x)*slope!=0;};
    const bool regularArc=ordinaryPorts&&simplePort(arc.first)&&simplePort(arc.last);
    const auto measured=regularArc?exactFixedPortArc(*solved,e,plain,lo,hi,rz):exactPolygon(*solved,exactCorners(box),rz);
    const QR actual=shape.sign>0?measured.below:measured.above;
    const QR residual=actual-target,own=std::min(target,capacity-target),error=residual<0?-residual:residual;
    require(error<=QE::scalar(1e-12L)*own&&error<=QE::scalar(PlanarCellCut::relativeVolumeResolution)*capacity,
        "[QuadraticCellCut] analytic fixed-port volume failed original volume/own-phase certificate");
    result.eta=QE::rounded(eta);result.chordPhaseVolume=QE::rounded(v0);result.volumeResidual=QE::rounded(residual);
    result.ownPhaseResidual=QE::rounded(error/own);result.volumeDerivativeEta=QE::rounded(2*a*eta+b);
    const std::array<QR,3> scale{{h,1,QR(1)/h}};
    for(int i=0;i<3;++i)result.selectedPhaseDerivative[i]=QE::rounded(shape.sign*measured.derivative[i]*scale[i]);
    result.graph.volume=QE::rounded(actual);result.graph.residual=QE::rounded(error/capacity);
    auto certificate=std::make_shared<Detail::QuadraticFixedPortCertificate>();
    certificate->owner=solved;certificate->box=box;certificate->ports={{first,last}};
    certificate->phaseSign=shape.sign;certificate->axisymmetric=rz;certificate->targetPhaseVolume=targetPhaseVolume;
    certificate->arcs=arcs;certificate->total=measured.total;certificate->below=measured.below;certificate->above=measured.above;
    solved->fixedPortCertificate=std::move(certificate);return result;
}
std::array<std::optional<StationaryFaceBoundary>,2> stationaryFaceBoundaries(const Rectangle&box,const Reconstruction&current){
    require(current.exact!=nullptr,"[QuadraticCellCut/StationaryFace] current exact graph required; legacy unsupported");
    Detail::validateGraph(current);(void)basis(current.shape);(void)PlanarCellCut::capacity(box,false);
    const auto arcs=Detail::graphArcs(*current.exact,box);const auto exact=exactStationaryFaceBoundaries(box,*current.exact,arcs);
    std::array<std::optional<StationaryFaceBoundary>,2> result;for(int axis=0;axis<2;++axis)if(exact[axis])result[axis]=nativeStationaryFaceBoundary(*exact[axis],axis);return result;
}
FixedPortVolumeBoundaryJet fixedPortVolumeBoundaryJet(const Rectangle&box,const Reconstruction&current,
    Point p0,Point p1,Point v0,Point v1,long double targetPhaseVolume,bool rz){
    require(current.exact!=nullptr,"[QuadraticCellCut/BoundaryJet] current exact graph required; legacy unsupported");
    Detail::validateGraph(current);(void)basis(current.shape);(void)PlanarCellCut::capacity(box,rz);
    require(std::isfinite(targetPhaseVolume),"[QuadraticCellCut/BoundaryJet] nonfinite original target phase");
    const auto* proof=current.exact->fixedPortCertificate.get();
    const auto samePoint=[](const auto& a,const auto& b){return a.x==b.x&&a.y==b.y;};
    const bool retained=proof&&proof->owner.lock().get()==current.exact.get()&&proof->phaseSign==current.shape.sign&&
        proof->axisymmetric==rz&&proof->targetPhaseVolume==targetPhaseVolume&&
        samePoint(proof->box.lower,box.lower)&&samePoint(proof->box.upper,box.upper)&&
        samePoint(proof->ports[0],p0)&&samePoint(proof->ports[1],p1);
    const auto graph=lifted(box,current);const auto&g=*graph;const QR h=std::max(QE::scalar(box.upper.x)-QE::scalar(box.lower.x),QE::scalar(box.upper.y)-QE::scalar(box.lower.y));
    require(h>0&&g.determinant>0,"[QuadraticCellCut/BoundaryJet] invalid current frame/cell");
    const std::array<Point,2> points{{p0,p1}},velocities{{v0,v1}};std::array<unsigned,2> face{};
    std::array<BoundaryJPoint,2> ports;std::array<BoundaryRJet,2> t{},n{};
    for(int i=0;i<2;++i){const auto p=QE::vector(points[i]),v=QE::vector(velocities[i]);
        require(p.x>=QE::scalar(box.lower.x)&&p.x<=QE::scalar(box.upper.x)&&p.y>=QE::scalar(box.lower.y)&&p.y<=QE::scalar(box.upper.y),"[QuadraticCellCut/BoundaryJet] port outside current cell");
        if(p.x==QE::scalar(box.lower.x))face[i]|=1;if(p.x==QE::scalar(box.upper.x))face[i]|=2;
        if(p.y==QE::scalar(box.lower.y))face[i]|=4;if(p.y==QE::scalar(box.upper.y))face[i]|=8;
        require(face[i]!=0,"[QuadraticCellCut/BoundaryJet] actual boundary port required");
        if((face[i]&(face[i]-1))!=0){
            require((face[i]&3)!=0&&(face[i]&12)!=0&&v.x==0&&v.y==0,
                "[QuadraticCellCut/BoundaryJet] actual corner velocity must be exactly zero");
        }else require((face[i]&3)?v.x==0:v.y==0,"[QuadraticCellCut/BoundaryJet] direction leaves original physical face");
        if(rz)require(p.x>0,"[QuadraticCellCut/BoundaryJet] physical-axis port belongs to independent fixed axis chart");
        ports[i]={{p.x,v.x,0},{p.y,v.y,0}};
        const BoundaryRJet dx=ports[i].x-BoundaryRJet(g.origin.x),dy=ports[i].y-BoundaryRJet(g.origin.y);
        t[i]=(dx*BoundaryRJet(g.tangent.x)+dy*BoundaryRJet(g.tangent.y))/BoundaryRJet(g.determinant);
        n[i]=(dx*BoundaryRJet(g.normal.x)+dy*BoundaryRJet(g.normal.y))/BoundaryRJet(g.determinant);
        require(n[i].value==(g.q[2]*t[i].value+g.q[1])*t[i].value+g.q[0],"[QuadraticCellCut/BoundaryJet] input port is not on the current exact curve");
    }
    require(t[0].value!=t[1].value,"[QuadraticCellCut/BoundaryJet] zero current tangent span");
    const bool forward=t[0].value<t[1].value;const auto lo=forward?t[0]:t[1],hi=forward?t[1]:t[0];
    const auto arcs=retained?proof->arcs:Detail::graphArcs(g,box);require(arcs.size()==1&&arcs[0].first.parameter==lo.value&&arcs[0].last.parameter==hi.value,"[QuadraticCellCut/BoundaryJet] current curve is not its original one-arc branch");
    const auto expectedFirst=forward?p0:p1,expectedLast=forward?p1:p0;const auto&arc=arcs[0];
    require(arc.first.position.x==expectedFirst.x&&arc.first.position.y==expectedFirst.y&&arc.last.position.x==expectedLast.x&&arc.last.position.y==expectedLast.y,"[QuadraticCellCut/BoundaryJet] current endpoints lost exact event identity");
    // A cache hit still runs the unchanged mixed/own-phase certificate below;
    // only their repeated arc construction and polygon integral are skipped.
    const auto measured=retained?ExactMeasure{proof->total,proof->below,proof->above,{}}:exactPolygon(g,exactCorners(box),rz);
    const QR target=QE::scalar(targetPhaseVolume),actual=current.shape.sign>0?measured.below:measured.above;
    require(target>0&&target<measured.total,"[QuadraticCellCut/BoundaryJet] actual mixed target required");
    const QR own=std::min(target,measured.total-target),error=actual>=target?actual-target:target-actual;
    require(error<=QE::scalar(1e-12L)*own&&error<=QE::scalar(PlanarCellCut::relativeVolumeResolution)*measured.total,"[QuadraticCellCut/BoundaryJet] current graph fails original target-phase certificate");
    const auto chordSlope=(n[1]-n[0])/(t[1]-t[0]),chordOffset=n[0]-chordSlope*t[0];
    const BoundaryJPoly bubble{{t[0]*t[1]/BoundaryRJet(h),-(t[0]+t[1])/BoundaryRJet(h),BoundaryRJet(QR(1)/h),0,0}};
    BoundaryRJet a=0,b=0;const BoundaryRJet metric(QR(current.shape.sign)*g.determinant);
    if(!rz)b=metric*boundaryIntegral(bubble,lo,hi);
    else{const BoundaryJPoly radius{{BoundaryRJet(g.origin.x)+BoundaryRJet(g.normal.x)*chordOffset,
        BoundaryRJet(g.tangent.x)+BoundaryRJet(g.normal.x)*chordSlope,0,0,0}};
        const auto revolved=metric*BoundaryRJet(QE::scalar(CartesianMeasure::twoPi));
        b=revolved*boundaryIntegral(boundaryProduct(bubble,radius),lo,hi);
        a=revolved*BoundaryRJet(g.normal.x)*boundaryIntegral(boundaryProduct(bubble,bubble),lo,hi)/2;}
    const std::array<QV,4> corners{{QE::vector({box.lower.x,box.lower.y}),QE::vector({box.upper.x,box.lower.y}),QE::vector({box.upper.x,box.upper.y}),QE::vector({box.lower.x,box.upper.y})}};
    BoundaryRJet chordVolume;
    if((face[0]&face[1])!=0){
        // Same actual face: every moving chord is that EXACT fixed boundary,
        // so its selected volume and both derivatives are constant.
        const auto middle=QE::subtract(exactCenter(box),g.origin);const QR ct=QE::dot(middle,g.tangent)/g.determinant,cn=QE::dot(middle,g.normal)/g.determinant;
        const QR sign=QR(current.shape.sign)*(cn-chordOffset.value-chordSlope.value*ct);
        require(sign!=0,"[QuadraticCellCut/BoundaryJet] degenerate same-face chord interior");chordVolume=BoundaryRJet(sign<0?measured.total:QR(0));
    }else{
        std::array<bool,4> inside{};for(int i=0;i<4;++i){const auto d=QE::subtract(corners[i],g.origin);const QR ct=QE::dot(d,g.tangent)/g.determinant,cn=QE::dot(d,g.normal)/g.determinant;
            const QR sign=QR(current.shape.sign)*(cn-chordOffset.value-chordSlope.value*ct);
            bool fixedPort=false;for(int port=0;port<2;++port)fixedPort=fixedPort||
                ((face[port]&(face[port]-1))!=0&&ports[port].x.value==corners[i].x&&ports[port].y.value==corners[i].y&&
                 ports[port].x.first==0&&ports[port].y.first==0);
            require(sign!=0||fixedPort,"[QuadraticCellCut/BoundaryJet] chord crosses an undeclared actual corner");inside[i]=sign<0;}
        const std::array<unsigned,4> edge{{4,2,8,1}};std::vector<BoundaryJPoint> polygon;
        for(int i=0;i<4;++i){if(inside[i])polygon.push_back({BoundaryRJet(corners[i].x),BoundaryRJet(corners[i].y)});
            if(inside[i]!=inside[(i+1)%4]){require((face[0]&edge[i])!=0||(face[1]&edge[i])!=0,"[QuadraticCellCut/BoundaryJet] chord event disagrees with original port face");polygon.push_back(ports[(face[0]&edge[i])!=0?0:1]);}}
        chordVolume=boundaryPolygonVolume(polygon,rz);
    }
    const auto c=chordVolume-BoundaryRJet(target);const QR eta0=g.q[2]*h,fEta=2*a.value*eta0+b.value;
    require(fEta!=0,"[QuadraticCellCut/BoundaryJet] singular current implicit volume branch");
    const QR fLambda=a.first*eta0*eta0+b.first*eta0+c.first,etaFirst=-fLambda/fEta;
    const QR fLambdaLambda=a.second*eta0*eta0+b.second*eta0+c.second,fLambdaEta=2*a.first*eta0+b.first;
    const QR etaSecond=-(fLambdaLambda+2*fLambdaEta*etaFirst+2*a.value*etaFirst*etaFirst)/fEta;
    const BoundaryRJet eta(eta0,etaFirst,etaSecond),q1=chordSlope-eta*(t[0]+t[1])/BoundaryRJet(h),q0=chordOffset+eta*t[0]*t[1]/BoundaryRJet(h),q2=eta/BoundaryRJet(h);
    FixedPortVolumeBoundaryJet result;result.eta=nativeBoundaryJet(eta);result.slope=nativeBoundaryJet(q1);result.span=nativeBoundaryJet((hi-lo)/BoundaryRJet(h));
    for(int i=0;i<4;++i){const auto d=QE::subtract(corners[i],g.origin);const BoundaryRJet ct(QE::dot(d,g.tangent)/g.determinant),cn(QE::dot(d,g.normal)/g.determinant);
        result.cornerImplicit[i]=nativeBoundaryJet(BoundaryRJet(current.shape.sign)*(cn-q0-q1*ct-q2*ct*ct)/BoundaryRJet(h));}
    for(int i=0;i<2;++i){const auto slope=q1+2*q2*t[i];result.portTangent[i][0]=nativeBoundaryJet(BoundaryRJet(g.tangent.x)+BoundaryRJet(g.normal.x)*slope);
        result.portTangent[i][1]=nativeBoundaryJet(BoundaryRJet(g.tangent.y)+BoundaryRJet(g.normal.y)*slope);}
    const auto stationary=exactStationaryFaceBoundaries(box,g,arcs);
    for(int axis=0;axis<2;++axis)if(stationary[axis]){
        const auto&actual=*stationary[axis];StationaryFaceBoundaryJet out;out.current=nativeStationaryFaceBoundary(actual,axis);
        const BoundaryRJet normal(axis?g.normal.y:g.normal.x),tangent(axis?g.tangent.y:g.tangent.x),origin(axis?g.origin.y:g.origin.x);
        const auto a=normal*q2,b=tangent+normal*q1,parameter=-b/(2*a);
        const auto position=origin+tangent*parameter+normal*(q0+q1*parameter+q2*parameter*parameter);
        const BoundaryRJet lower(QE::scalar(axis?box.lower.y:box.lower.x)),upper(QE::scalar(axis?box.upper.y:box.upper.x));
        for(int side=0;side<2;++side){const auto inward=side?upper-position:position-lower;
            out.slack[side]=nativeBoundaryJet(BoundaryRJet(actual.branchSign[side])*inward/BoundaryRJet(h));}
        result.stationaryFace[axis]=std::move(out);
    }
    return result;
}
Reconstruction fromAxisPorts(const Rectangle&box,int sign,Point axis,Point other){
    require(axis.x==0&&other.x!=0&&(sign==-1||sign==1),"[QuadraticCellCut] invalid axis ports");auto e=std::make_shared<Detail::QuadraticGraphExact>();e->origin=exactCenter(box);e->tangent={1,0};e->normal={0,1};
    const QR k=(QE::scalar(other.y)-QE::scalar(axis.y))/(QE::scalar(other.x)*QE::scalar(other.x));e->q={QE::scalar(axis.y)-e->origin.y+k*e->origin.x*e->origin.x,2*k*e->origin.x,k};
    e->portParameters={-e->origin.x,QE::scalar(other.x)-e->origin.x};e->ports={axis,other};Reconstruction result;result.shape.axis=1;result.shape.sign=sign;setViews(result,e);return result;
}
AxisPortVolumeResult axisPortVolumeGraph(const Rectangle&box,int sign,Point axis,
    Point begin,Point end,long double targetPhaseVolume){
    const auto finite=[](Point p){return std::isfinite(p.x)&&std::isfinite(p.y);};
    require(std::isfinite(box.lower.x)&&std::isfinite(box.lower.y)&&std::isfinite(box.upper.x)&&std::isfinite(box.upper.y)&&
        box.lower.x==0&&box.upper.x>0&&box.lower.y<box.upper.y,
        "[QuadraticCellCut] axis volume requires a finite physical rlo=0 box");
    require((sign==-1||sign==1)&&finite(axis)&&finite(begin)&&finite(end)&&std::isfinite(targetPhaseVolume),
        "[QuadraticCellCut] nonfinite or invalid axis volume request");
    const QR r=QE::scalar(box.upper.x),zlo=QE::scalar(box.lower.y),zhi=QE::scalar(box.upper.y),za=QE::scalar(axis.y);
    require(axis.x==0&&za>zlo&&za<zhi,
        "[QuadraticCellCut] axis volume point must be on the actual open physical axis face");
    const bool horizontal=begin.y==end.y&&begin.x!=end.x;
    const bool vertical=begin.x==end.x&&begin.y!=end.y;
    require(horizontal||vertical,"[QuadraticCellCut] axis companion needs one nondegenerate axis-aligned actual face segment");
    if(horizontal){
        require((begin.y==box.lower.y||begin.y==box.upper.y)&&begin.x>=0&&begin.x<=box.upper.x&&end.x>=0&&end.x<=box.upper.x,
            "[QuadraticCellCut] horizontal axis companion segment must lie on the actual top or bottom cell face");
    }else{
        require(begin.x==box.upper.x&&begin.y>=box.lower.y&&begin.y<=box.upper.y&&end.y>=box.lower.y&&end.y<=box.upper.y,
            "[QuadraticCellCut] vertical axis companion segment must lie on the actual right cell face");
    }
    const QR pi=QE::scalar(CartesianMeasure::twoPi)/2,r2=r*r,capacity=pi*r2*(zhi-zlo),target=QE::scalar(targetPhaseVolume);
    require(target>0&&target<capacity,"[QuadraticCellCut] axis volume requires an actual mixed selected-phase target");
    QR k;unsigned boundary=2;
    if(!horizontal){
        // Full-r below volume = pi R^2 (za-zlo) + pi k R^4/2.
        const QR below=sign>0?target:capacity-target;
        k=2*(below-pi*r2*(za-zlo))/(pi*r2*r2);
    }else{
        const bool top=begin.y==box.upper.y;boundary=top?8u:4u;
        const QR delta=(top?zhi:zlo)-za;
        // The actual outside wedge is pi |delta| rb^2/2. Select it
        // directly from the requested sign; all complements remain rational.
        const QR wedge=top?(sign>0?capacity-target:target):(sign>0?target:capacity-target);
        const QR absoluteDelta=delta>0?delta:-delta,radialSquared=2*wedge/(pi*absoluteDelta);
        require(radialSquared>0,"[QuadraticCellCut] axis companion has no positive radial-square solution");
        k=delta/radialSquared;
    }
    auto exact=std::make_shared<Detail::QuadraticGraphExact>();exact->origin=exactCenter(box);exact->tangent={1,0};exact->normal={0,1};
    exact->q={za-exact->origin.y+k*exact->origin.x*exact->origin.x,2*k*exact->origin.x,k};
    const QR axisParameter=-exact->origin.x;
    exact->portParameters={axisParameter};exact->ports={axis};
    AxisPortVolumeResult result;result.horizontalCompanion=horizontal;result.graph.shape.axis=1;result.graph.shape.sign=sign;setViews(result.graph,exact);
    const auto arcs=Detail::graphArcs(*exact,box);
    require(arcs.size()==1&&arcs.front().first.parameter==axisParameter&&
        arcs.front().first.boundary==1&&arcs.front().last.boundary==boundary,
        "[QuadraticCellCut] analytic axis volume leaves its original one-arc/face branch");
    const auto& arc=arcs.front();const auto companion=arc.last.position;
    require(arc.first.position.x==axis.x&&arc.first.position.y==axis.y&&companion.x>0&&
        arc.last.parameter>arc.first.parameter&&QE::rounded(arc.last.parameter)>QE::rounded(arc.first.parameter),
        "[QuadraticCellCut] analytic axis volume has no representable positive-radius regular arc");
    require(Detail::graphSlope(*exact,axisParameter)==0,
        "[QuadraticCellCut] analytic axis volume violates exact dz/dr=0");
    const auto point=Detail::graphPosition(*exact,arc.last.parameter);
    const QR coordinate=horizontal?point.x:point.y,first=QE::scalar(horizontal?begin.x:begin.y),last=QE::scalar(horizontal?end.x:end.y);
    const QR fraction=(coordinate-first)/(last-first);
    result.companionFraction=QE::rounded(fraction);
    const auto nativeCoordinate=horizontal?companion.x:companion.y;
    require(fraction>0&&fraction<1&&result.companionFraction>0&&result.companionFraction<1&&
        nativeCoordinate>std::min(horizontal?begin.x:begin.y,horizontal?end.x:end.y)&&
        nativeCoordinate<std::max(horizontal?begin.x:begin.y,horizontal?end.x:end.y),
        "[QuadraticCellCut] analytic axis companion leaves its actual open subface");
    if(horizontal)require(k!=0&&Detail::graphSlope(*exact,arc.last.parameter)!=0,
        "[QuadraticCellCut] analytic axis companion is tangent to its horizontal face");
    // Reuse the same side event for all later surface/material exports. The
    // physical polynomial above stays unchanged; no rounded-root refit occurs.
    exact->portParameters.push_back(arc.last.parameter);exact->ports.push_back(companion);
    const auto measured=exactPolygon(*exact,exactCorners(box),true);
    const QR actual=sign>0?measured.below:measured.above,residual=actual-target;
    const QR own=std::min(target,capacity-target),error=residual<0?-residual:residual;
    require(error<=QE::scalar(1e-12L)*own&&error<=QE::scalar(PlanarCellCut::relativeVolumeResolution)*capacity,
        "[QuadraticCellCut] analytic axis volume failed original volume/own-phase certificate");
    result.companion=companion;result.phaseVolume=QE::rounded(actual);result.volumeResidual=QE::rounded(residual);result.ownPhaseResidual=QE::rounded(error/own);
    result.graph.volume=result.phaseVolume;result.graph.residual=QE::rounded(error/capacity);return result;
}
Reconstruction recenter(const Rectangle&source,const Rectangle&target,const Reconstruction&g){
    auto e=lifted(source,g);const auto origin=exactCenter(target),delta=QE::subtract(origin,e->origin);const QR t=QE::dot(delta,e->tangent)/e->determinant,n=QE::dot(delta,e->normal)/e->determinant;
    e->q={e->q[0]+e->q[1]*t+e->q[2]*t*t-n,e->q[1]+2*e->q[2]*t,e->q[2]};e->origin=origin;for(auto&port:e->portParameters)port-=t;
    auto result=g;setViews(result,e);return result;
}
GraphMeasure graphMeasure(const Rectangle&box,const Reconstruction&g,bool rz){
    const auto e=lifted(box,g);const auto m=exactPolygon(*e,exactCorners(box),rz);GraphMeasure result;result.volume=QE::rounded(g.shape.sign>0?m.below:m.above);
    const QR h=std::max(QE::scalar(box.upper.x)-QE::scalar(box.lower.x),QE::scalar(box.upper.y)-QE::scalar(box.lower.y));const std::array<QR,3> scale{{h,1,QR(1)/h}};
    for(int i=0;i<3;++i)result.derivative[i]=QE::rounded(g.shape.sign*m.derivative[i]*scale[i]);return result;
}
std::shared_ptr<const Detail::QuadraticParametricData> materialCurveData(const Rectangle&box,const Reconstruction&g,long double lower,long double upper){
    const auto e=lifted(box,g);const auto events=Detail::graphBoundaryEvents(*e,box);const Detail::GraphEvent *first=nullptr,*last=nullptr;
    for(const auto&event:events){if(!event.boundary)continue;
        if(QE::rounded(event.parameter)==lower){require(!first,"[QuadraticCellCut] ambiguous native material lower event");first=&event;}
        if(QE::rounded(event.parameter)==upper){require(!last,"[QuadraticCellCut] ambiguous native material upper event");last=&event;}}
    require(first&&last&&first->parameter<last->parameter,"[QuadraticCellCut] material interval has no original side-root identities");return Detail::graphArcMaterialData(*e,{*first,*last});
}
CornerBoundaryKind cornerBoundaryKind(const Rectangle& box,const Reconstruction& graph,unsigned corner)
{
    require(graph.exact!=nullptr,"[QuadraticCellCut/CornerBoundary] current exact graph required; legacy unsupported");
    Detail::validateGraph(graph);(void)basis(graph.shape);(void)PlanarCellCut::capacity(box,false);
    require(corner<4,"[QuadraticCellCut/CornerBoundary] invalid corner index");
    const auto& e=*graph.exact;
    require(e.determinant>0,"[QuadraticCellCut/CornerBoundary] invalid exact frame determinant");
    const std::array<Point,4> corners{{{box.lower.x,box.lower.y},{box.upper.x,box.lower.y},
        {box.upper.x,box.upper.y},{box.lower.x,box.upper.y}}};
    const auto delta=QE::subtract(QE::vector(corners[corner]),e.origin);
    const QR t=QE::dot(delta,e.tangent)/e.determinant,n=QE::dot(delta,e.normal)/e.determinant;
    const QR f=graph.shape.sign*(n-(e.q[2]*t+e.q[1])*t-e.q[0]);
    const QR slope=e.q[1]+2*e.q[2]*t;
    const int inwardX=(corner==0||corner==3)?1:-1,inwardY=corner<2?1:-1;
    const QR gx=graph.shape.sign*inwardX*(e.normal.x-slope*e.tangent.x)/e.determinant;
    const QR gy=graph.shape.sign*inwardY*(e.normal.y-slope*e.tangent.y)/e.determinant;
    const bool positive=gx>0&&gy>0,negative=gx<0&&gy<0;
    const bool closed=(positive&&(f>=0))||(negative&&(f<=0));
    const QR h=std::max(QE::scalar(box.upper.x)-QE::scalar(box.lower.x),QE::scalar(box.upper.y)-QE::scalar(box.lower.y));
    return {QE::rounded(f/h),{{QE::rounded(gx),QE::rounded(gy)}},closed};
}
long double implicitValue(const Rectangle& box,const Reconstruction& graph,Point point)
{
    require(std::isfinite(box.lower.x) && std::isfinite(box.lower.y) && std::isfinite(box.upper.x)
        && std::isfinite(box.upper.y) && box.lower.x<box.upper.x && box.lower.y<box.upper.y
        && std::isfinite(graph.offset) && std::isfinite(point.x) && std::isfinite(point.y),
        "[QuadraticCellCut] invalid phase query");
    if(graph.exact){Detail::validateGraph(graph);const auto&e=*graph.exact;const auto delta=QE::subtract(QE::vector(point),e.origin);const QR t=QE::dot(delta,e.tangent)/e.determinant,n=QE::dot(delta,e.normal)/e.determinant;return QE::rounded(graph.shape.sign*(n-(e.q[2]*t+e.q[1])*t-e.q[0]));}
    const auto b=basis(graph.shape);const auto c=center(box);point.x-=c.x;point.y-=c.y;
    const auto t=point.x*b.tangent.x+point.y*b.tangent.y,n=point.x*b.normal.x+point.y*b.normal.y;
    return graph.shape.sign*(n-graph.offset-graph.shape.slope*t-graph.shape.second*t*t/2);
}
namespace {
GraphMeasure measureGraphPolygon(const Rectangle& box,const std::vector<Point>& polygon,Point polygonOrigin,
    Point origin,Vector2D normal,long double scale,const std::array<long double,3>& coefficients,bool rz,int sign,bool* fullOutput)
{
    const auto capacity=PlanarCellCut::capacity(box,rz);
    const auto length=std::hypot(static_cast<long double>(normal.x),static_cast<long double>(normal.y));
    require(std::isfinite(scale) && scale>0 && std::isfinite(length) && length>0
        && std::isfinite(origin.x) && std::isfinite(origin.y),"[QuadraticCellCut] invalid rotated graph frame");
    for (auto c:coefficients) require(std::isfinite(c),"[QuadraticCellCut] nonfinite rotated graph coefficient");
    const Point n{normal.x/length,normal.y/length},t{-n.y,n.x};
    std::vector<Point> local;local.reserve(polygon.size());
    std::vector<long double> slabs;
    for (auto point:polygon) {
        const auto x=(point.x-(origin.x-polygonOrigin.x))/scale,y=(point.y-(origin.y-polygonOrigin.y))/scale;
        local.push_back({x*t.x+y*t.y,x*n.x+y*n.y});slabs.push_back(local.back().x);
    }
    std::sort(slabs.begin(),slabs.end());slabs.erase(std::unique(slabs.begin(),slabs.end()),slabs.end());
    GraphMeasure result;bool full=true,empty=true;
    for (std::size_t k=1;k<slabs.size();++k) {
        const auto lo=slabs[k-1],hi=slabs[k],width=hi-lo;
        // Integrate in the slab's own [0,1] parameter. Nearly aligned frames
        // can create sub-ULP-width slivers: no rounded interior abscissa or
        // enormous global edge intercept is needed to represent their area.
        const Poly graph{{(coefficients[2]*lo+coefficients[1])*lo+coefficients[0],
            (coefficients[1]+2*coefficients[2]*lo)*width,coefficients[2]*width*width,0,0}};
        std::vector<Poly> edges;
        for (std::size_t i=0;i<local.size();++i) {
            const auto a=local[i],b=local[(i+1)%local.size()];
            if (a.x==b.x || lo<std::min(a.x,b.x) || hi>std::max(a.x,b.x)) continue;
            edges.push_back({{a.y+(lo-a.x)/(b.x-a.x)*(b.y-a.y),width/(b.x-a.x)*(b.y-a.y),0,0,0}});
        }
        if (edges.empty()) continue;
        require(edges.size()==2,"[QuadraticCellCut] invalid rotated cell slice");
        if (value(edges[0],.5L)>value(edges[1],.5L)) std::swap(edges[0],edges[1]);
        std::vector<long double> cuts{0,1};
        for (const auto& edge:edges) roots(graph[2],graph[1]-edge[1],graph[0]-edge[0],0,1,cuts);
        std::sort(cuts.begin(),cuts.end());cuts.erase(std::unique(cuts.begin(),cuts.end()),cuts.end());
        for (std::size_t j=1;j<cuts.size();++j) {
            const auto a=cuts[j-1],b=cuts[j],mid=a+(b-a)/2;
            auto lower=edges[0],upper=edges[1];bool cut=false;
            if(sign>0) {
                if(value(graph,mid)<=value(lower,mid)){full=false;continue;}
                cut=value(graph,mid)<value(upper,mid);if(cut)upper=graph;
            } else {
                if(value(graph,mid)>=value(upper,mid)){full=false;continue;}
                cut=value(graph,mid)>value(lower,mid);if(cut)lower=graph;
            }
            empty=false;if(cut)full=false;
            Poly difference{},density{};
            for(int i=0;i<5;++i) difference[i]=upper[i]-lower[i];
            if (rz) {
                Poly meanRadius{};
                for(int i=0;i<5;++i) meanRadius[i]=scale*n.x*(upper[i]+lower[i])/2;
                meanRadius[0]+=origin.x+scale*t.x*lo;meanRadius[1]+=scale*t.x*width;
                for(int i=0;i<5;++i) for(int m=0;i+m<5;++m)
                    density[i+m]+=CartesianMeasure::twoPi*difference[i]*meanRadius[m];
            } else density=difference;
            result.volume+=scale*scale*width*integrate(density,a,b);
            if (cut) {
                // Endpoint variations cancel because the two slice formulas
                // agree where the graph meets the rectangle boundary.
                Poly radius{};radius[0]=1;
                if (rz) {
                    radius={};for(int i=0;i<3;++i) radius[i]=scale*n.x*graph[i];
                    radius[0]+=origin.x+scale*t.x*lo;radius[1]+=scale*t.x*width;
                    for(auto& q:radius) q*=CartesianMeasure::twoPi;
                }
                const std::array<Poly,3> basis{{Poly{{1,0,0,0,0}},Poly{{lo,width,0,0,0}},
                    Poly{{lo*lo,2*lo*width,width*width,0,0}}}};
                for(int m=0;m<3;++m) {
                    Poly sensitivity{};
                    for(int i=0;i<5;++i) for(int p=0;i+p<5;++p) sensitivity[i+p]+=radius[i]*basis[m][p];
                    result.derivative[m]+=sign*scale*scale*width*integrate(sensitivity,a,b);
                }
            }
        }
    }
    if(fullOutput)*fullOutput=full;if(empty)result.volume=0;
    const auto allowance=PlanarCellCut::relativeVolumeResolution*capacity;
    require(std::isfinite(result.volume) && result.volume>=-allowance && result.volume<=capacity+allowance,
        "[QuadraticCellCut] invalid rotated graph volume");
    for(auto d:result.derivative) require(std::isfinite(d),"[QuadraticCellCut] nonfinite graph sensitivity");
    return result;
}
}
GraphMeasure graphMeasure(const Rectangle& box,Point origin,Vector2D normal,long double scale,
    const std::array<long double,3>& coefficients,bool rz)
{
    const std::vector<Point> corners{{box.lower.x,box.lower.y},{box.upper.x,box.lower.y},
        {box.upper.x,box.upper.y},{box.lower.x,box.upper.y}};
    bool full=false;auto result=measureGraphPolygon(box,corners,{0,0},origin,normal,scale,coefficients,rz,1,&full);
    if(full)result.volume=PlanarCellCut::capacity(box,rz);return result;
}
PlanarCellCut::Overlap overlapInFrame(const Rectangle& box,const std::vector<Point>& polygon,Point origin,
    const Reconstruction& reconstruction,bool rz)
{
    if(reconstruction.exact){Detail::validateGraph(reconstruction);const auto clipped=PlanarCellCut::clipPolygonInFrame(box,polygon,origin);std::vector<QV> exact;for(const auto&p:clipped)exact.push_back(QE::add(QE::vector(p),QE::vector(origin)));const auto m=exactPolygon(*reconstruction.exact,exact,rz);return {QE::rounded(m.total),QE::rounded(reconstruction.shape.sign>0?m.below:m.above)};}
    validate(reconstruction.shape);require(std::isfinite(reconstruction.offset),"[QuadraticCellCut] invalid offset");
    const auto total=PlanarCellCut::overlapInFrame(box,polygon,origin,{1,0},0,rz).total;
    const auto clipped=PlanarCellCut::clipPolygonInFrame(box,polygon,origin);
    bool entirelyLiquid=false;
    const auto integrated=polygonVolume(box,clipped,origin,reconstruction,rz,&entirelyLiquid);
    // A region classified wholly liquid shares the exact total metric. This
    // avoids comparing different integration roundoff for the same domain;
    // genuinely cut regions retain their independently integrated volume.
    const auto liquid=entirelyLiquid?total:integrated;
    require(liquid<=total+PlanarCellCut::relativeVolumeResolution*PlanarCellCut::capacity(box,rz),
        "[QuadraticCellCut] curved subvolume exceeds donor overlap");
    // No phase clipping: the caller retains the same strict shared-face bounds.
    return {total,liquid};
}
namespace
{
namespace StripExact=LAMP::Geometry::Detail::ExactExpansion;
namespace StripRound=LAMP::Geometry::Detail::ExactRationalRound;
using SE=StripExact::Expansion;
using SP=std::array<SE,5>;
SE stripScalar(long double value){return StripRound::scalar(value);}
SE stripAdd(const SE& a,const SE& b){return StripExact::add(a,b);}
SE stripSub(const SE& a,const SE& b){return StripExact::subtract(a,b);}
SE stripMul(const SE& a,const SE& b){return StripExact::multiply(a,b);}
SE stripScale(const SE& a,double b){return stripMul(a,StripExact::fromDouble(b));}
SP stripProduct(const SP& a,const SP& b)
{
    SP result{};for(int i=0;i<5;++i)for(int j=0;i+j<5;++j)
        result[i+j]=stripAdd(result[i+j],stripMul(a[i],b[j]));return result;
}
SE stripValue(const SP& p,const SE& at)
{auto result=p[4];for(int i=3;i>=0;--i)result=stripAdd(p[i],stripMul(at,result));return result;}
// Common integer denominator preserves exact cancellation of close interval
// primitives and thin positive phases. No large-coordinate subtraction rounds.
SE stripIntegral60(const SP& p,const SE& x,const SE& y)
{
    auto a=StripExact::fromDouble(1),b=a;SE numerator;
    for(int i=0;i<5;++i){a=stripMul(a,x);b=stripMul(b,y);
        numerator=stripAdd(numerator,stripScale(stripMul(p[i],stripSub(b,a)),60./(i+1)));}
    return numerator;
}
SE stripIntegral60(const SP& p,long double lo,long double hi)
{return stripIntegral60(p,stripScalar(lo),stripScalar(hi));}
// Each event is retained as its distance from the nearer endpoint. In
// particular, 1-distance is an exact expansion, not a rounded unit coordinate.
// The root model remains the original floating quadratic root approximation.
void stripRoots(const SP& p,std::vector<SE>& cuts)
{
    const auto unit=StripExact::fromDouble(1);
    const auto rounded=[&](const SE& value){return StripRound::quotient(value,unit);};
    std::vector<long double> left,right;
    roots(rounded(p[2]),rounded(p[1]),rounded(p[0]),0,1,left);
    const auto atOne=stripAdd(stripAdd(p[0],p[1]),p[2]);
    const auto reverseLinear=stripSub(SE{},stripAdd(p[1],stripScale(p[2],2)));
    roots(rounded(p[2]),rounded(reverseLinear),rounded(atOne),0,1,right);
    for(auto root:left)if(root<=.5L)cuts.push_back(stripScalar(root));
    for(auto distance:right)if(distance<.5L)cuts.push_back(stripSub(unit,stripScalar(distance)));
}
void stripSortCuts(std::vector<SE>& cuts)
{
    std::sort(cuts.begin(),cuts.end(),[](const SE& a,const SE& b){return StripExact::sign(stripSub(a,b))<0;});
    cuts.erase(std::unique(cuts.begin(),cuts.end(),[](const SE& a,const SE& b){return StripExact::sign(stripSub(a,b))==0;}),cuts.end());
}
// Preselected exact-integer rational backend for a genuine frame graph.
// Reduced arbitrary-size integer ratios do not inherit the expansion engine's
// binary64 exponent window; no floating homogeneous product is accumulated.
using FrameRatio=boost::multiprecision::number<boost::multiprecision::cpp_rational_backend,boost::multiprecision::et_off>;
using FrameInteger=boost::multiprecision::cpp_int;
using FramePoly=std::array<FrameRatio,5>;
FrameRatio frameScalar(long double value)
{
    require(std::isfinite(value),"[QuadraticCellCut] nonfinite exact framed scalar");if(value==0)return 0;
    const bool negative=value<0;int exponent=0;auto fraction=std::frexp(std::abs(value),&exponent);FrameInteger bits=0;
    for(int k=0;k<std::numeric_limits<long double>::digits;++k){fraction*=2;bits<<=1;if(fraction>=1){++bits;fraction-=1;}}
    require(fraction==0,"[QuadraticCellCut] native framed scalar bits were not retained");
    FrameRatio result(bits);const int shift=exponent-std::numeric_limits<long double>::digits;
    if(shift>=0)result*=FrameRatio(FrameInteger(1)<<shift);else result/=FrameRatio(FrameInteger(1)<<-shift);
    return negative?-result:result;
}
// Boost conversion is only a seed. Exact rational midpoint comparisons select
// the actual native nearest/ties-even representative, including subnormals.
long double frameRounded(const FrameRatio& input)
{
    if(input==0)return 0;const bool negative=input<0;const FrameRatio value=negative?-input:input;
    auto candidate=value.convert_to<long double>();if(candidate==0)candidate=std::numeric_limits<long double>::denorm_min();
    if(!std::isfinite(candidate))candidate=std::numeric_limits<long double>::max();
    for(int step=0;step<16;++step){const auto lower=std::nextafter(candidate,0.L),upper=std::nextafter(candidate,std::numeric_limits<long double>::infinity());
        const bool even=StripRound::evenSignificand(candidate);const FrameRatio lowerMid=(frameScalar(lower)+frameScalar(candidate))/2;
        if(value<lowerMid||(value==lowerMid&&!even)){require(lower>0,"[QuadraticCellCut] nonzero framed measure rounds to zero");candidate=lower;continue;}
        const FrameRatio upperMid=std::isfinite(upper)?(frameScalar(candidate)+frameScalar(upper))/2:
            frameScalar(candidate)+frameScalar(std::scalbn(1.L,std::numeric_limits<long double>::max_exponent-std::numeric_limits<long double>::digits-1));
        if(value>upperMid||(value==upperMid&&!even)){require(std::isfinite(upper),"[QuadraticCellCut] framed native measure overflows");candidate=upper;continue;}
        return negative?-candidate:candidate;
    }
    throw std::invalid_argument("[QuadraticCellCut] framed native rounding seed is outside its certified neighborhood");
}
FramePoly frameProduct(const FramePoly& a,const FramePoly& b)
{FramePoly result{};for(int i=0;i<5;++i)for(int j=0;i+j<5;++j)result[i+j]+=a[i]*b[j];return result;}
FrameRatio frameValue(const FramePoly& p,const FrameRatio& at)
{auto result=p[4];for(int i=3;i>=0;--i)result=p[i]+at*result;return result;}
FrameRatio frameIntegral(const FramePoly& p,const FrameRatio& lo,const FrameRatio& hi)
{FrameRatio a=1,b=1,result=0;for(int i=0;i<5;++i){a*=lo;b*=hi;result+=p[i]*(b-a)/(i+1);}return result;}
void frameRoots(const FramePoly& p,std::vector<FrameRatio>& cuts)
{
    std::vector<long double> left,right;roots(frameRounded(p[2]),frameRounded(p[1]),frameRounded(p[0]),0,1,left);
    roots(frameRounded(p[2]),frameRounded(-p[1]-2*p[2]),frameRounded(p[0]+p[1]+p[2]),0,1,right);
    for(auto root:left)if(root<=.5L)cuts.push_back(frameScalar(root));
    for(auto distance:right)if(distance<.5L)cuts.push_back(FrameRatio(1)-frameScalar(distance));
}
SP stripDensity(const SP& height,const SP& startRadius,const SE& twoPi,bool rz)
{
    if(!rz)return height;
    SP radius=startRadius;for(int k=0;k<5;++k)radius[k]=stripAdd(radius[k],stripScale(height[k],.5));
    auto density=stripProduct(height,radius);for(auto& value:density)value=stripMul(twoPi,value);return density;
}
}
StripPhaseMeasure linearSignedNormalSweepMeasure(const Rectangle&box,const Reconstruction&g,int normal,
    Point a,Point b,long double meanDepth,long double halfVariationDepth,bool positive,bool rz,
    std::optional<std::array<Point,2>> query,bool identifyCoveredPhases)
{
    require(std::isfinite(meanDepth)&&std::isfinite(halfVariationDepth),"[QuadraticCellCut/LinearSweep] nonfinite factored depth");
    const QR mean=QE::scalar(meanDepth),variation=QE::scalar(halfVariationDepth);
    return exactLinearSweep(box,g,normal,a,b,mean-variation,mean+variation,positive,rz,query,identifyCoveredPhases);
}
StripPhaseMeasure linearNormalSweepMeasure(const Rectangle&box,const Reconstruction&g,int normal,
    Point a,Point b,long double depthBegin,long double depthEnd,bool positive,bool rz,
    std::optional<std::array<Point,2>> query,bool identifyCoveredPhases)
{
    require(std::isfinite(depthBegin)&&std::isfinite(depthEnd)&&depthBegin>=0&&depthEnd>=0,
        "[QuadraticCellCut/LinearSweep] endpoint depths must be nonnegative");
    const QR direction=positive?QR(1):QR(-1);
    return exactLinearSweep(box,g,normal,a,b,direction*QE::scalar(depthBegin),direction*QE::scalar(depthEnd),
        positive,rz,query,identifyCoveredPhases);
}
StripPhaseMeasure axisStripMeasure(const Rectangle& box,const Reconstruction& graph,int normal,
    Point a,Point b,long double depth,bool positive,bool rz)
{
    if(graph.exact)return exactStrip(box,graph,normal,a,b,depth,positive,rz);
    validate(graph.shape);require(!graph.shape.frameNormal,"[QuadraticCellCut] axis strip requires the genuine axis graph branch");
    require((normal==0||normal==1)&&std::isfinite(depth)&&depth>0&&std::isfinite(graph.offset),"[QuadraticCellCut] invalid separate strip depth/normal");
    require(std::isfinite(a.x)&&std::isfinite(a.y)&&std::isfinite(b.x)&&std::isfinite(b.y)
        &&coordinate(a,normal)==coordinate(b,normal)&&coordinate(a,1-normal)!=coordinate(b,1-normal),"[QuadraticCellCut] invalid axis strip near interval");
    StripExact::validateActualEnvironment();const auto unit=StripExact::fromDouble(1),sixty=StripExact::fromDouble(60);
    const auto near=stripScalar(coordinate(a,normal)),d=stripScalar(depth);
    std::array<SE,2> lower,width;lower[normal]=positive?stripSub(near,d):near;width[normal]=d;
    lower[1-normal]=stripScalar(std::min(coordinate(a,1-normal),coordinate(b,1-normal)));
    width[1-normal]=stripSub(stripScalar(std::max(coordinate(a,1-normal),coordinate(b,1-normal))),lower[1-normal]);
    for(int axis=0;axis<2;++axis){const auto boundLo=stripScalar(axis==0?box.lower.x:box.lower.y),boundHi=stripScalar(axis==0?box.upper.x:box.upper.y);
        require(StripExact::sign(stripSub(lower[axis],boundLo))>=0&&StripExact::sign(stripSub(boundHi,stripAdd(lower[axis],width[axis])))>=0,
            "[QuadraticCellCut] exact axis strip leaves its original donor");}
    if(rz)require(StripExact::sign(lower[0])>=0,"[QuadraticCellCut] exact strip crosses negative radius");
    const auto c=center(box);const int tangent=1-graph.shape.axis,n=graph.shape.axis;
    const auto t0=stripSub(lower[tangent],stripScalar(coordinate(c,tangent))),n0=stripSub(lower[n],stripScalar(coordinate(c,n)));
    const auto slope=stripScalar(graph.shape.slope),second=stripScalar(graph.shape.second),tw=width[tangent],nw=width[n];
    SP height{};height[0]=stripSub(stripAdd(stripScalar(graph.offset),stripAdd(stripMul(slope,t0),stripScale(stripMul(second,stripMul(t0,t0)),.5))),n0);
    height[1]=stripMul(stripAdd(slope,stripMul(second,t0)),tw);height[2]=stripScale(stripMul(second,stripMul(tw,tw)),.5);
    const auto rounded=[&](const SE& value){return StripRound::quotient(value,unit);};
    std::vector<SE> cuts{SE{},unit};stripRoots(height,cuts);auto upperEvent=height;upperEvent[0]=stripSub(upperEvent[0],nw);
    stripRoots(upperEvent,cuts);stripSortCuts(cuts);
    const auto twoPi=stripScalar(CartesianMeasure::twoPi);SE belowIntegral,aboveIntegral;
    for(std::size_t span=1;span<cuts.size();++span){const auto lo=cuts[span-1],hi=cuts[span];
        const auto mid=stripScale(stripAdd(lo,hi),.5);const auto h=stripValue(height,mid);
        SP below{},above{};
        if(StripExact::sign(h)<=0)above[0]=nw;
        else if(StripExact::sign(stripSub(h,nw))>=0)below[0]=nw;
        else{below=height;above[0]=stripSub(nw,height[0]);for(int k=1;k<5;++k)above[k]=stripSub(SE{},height[k]);}
        SP densityBelow{},densityAbove{};
        if(rz&&n==1){SP radius{};radius[0]=lower[0];radius[1]=width[0];densityBelow=stripProduct(below,radius);densityAbove=stripProduct(above,radius);
            for(auto* p:{&densityBelow,&densityAbove})for(auto& value:*p)value=stripMul(twoPi,value);}
        else{SP baseRadius{};baseRadius[0]=lower[0];densityBelow=stripDensity(below,baseRadius,twoPi,rz);
            for(int k=0;k<5;++k)baseRadius[k]=stripAdd(baseRadius[k],below[k]);densityAbove=stripDensity(above,baseRadius,twoPi,rz);}
        const auto lowerPart=stripIntegral60(densityBelow,lo,hi),upperPart=stripIntegral60(densityAbove,lo,hi);
        require(StripExact::sign(lowerPart)>=0&&StripExact::sign(upperPart)>=0,"[QuadraticCellCut] unresolved positive strip phase integral");
        belowIntegral=stripAdd(belowIntegral,lowerPart);aboveIntegral=stripAdd(aboveIntegral,upperPart);
    }
    belowIntegral=stripMul(tw,belowIntegral);aboveIntegral=stripMul(tw,aboveIntegral);
    auto total=stripMul(width[0],width[1]);if(rz)total=stripMul(total,stripMul(twoPi,stripAdd(lower[0],stripScale(width[0],.5))));
    require(StripExact::sign(stripSub(stripAdd(belowIntegral,aboveIntegral),stripScale(total,60)))==0,
        "[QuadraticCellCut] exact strip phase partition identity failed");
    StripPhaseMeasure result;result.total=rounded(total);
    const auto below=StripRound::quotient(belowIntegral,sixty),above=StripRound::quotient(aboveIntegral,sixty);
    result.liquid=graph.shape.sign>0?below:above;result.gas=graph.shape.sign>0?above:below;return result;
}
StripPhaseMeasure framedStripMeasure(const Rectangle& box,const Reconstruction& graph,int normal,
    Point a,Point b,long double depth,bool positive,bool rz)
{
    if(graph.exact)return exactStrip(box,graph,normal,a,b,depth,positive,rz);
    validate(graph.shape);require(graph.shape.frameNormal.has_value(),"[QuadraticCellCut] framed strip requires the genuine original frame graph");
    require((normal==0||normal==1)&&std::isfinite(depth)&&depth>0&&std::isfinite(graph.offset),"[QuadraticCellCut] invalid separate framed strip depth/normal");
    require(std::isfinite(a.x)&&std::isfinite(a.y)&&std::isfinite(b.x)&&std::isfinite(b.y)
        &&coordinate(a,normal)==coordinate(b,normal)&&coordinate(a,1-normal)!=coordinate(b,1-normal),"[QuadraticCellCut] invalid framed strip near interval");
    StripExact::validateActualEnvironment(); // Native source bits and root rounding only; no expansion framed arithmetic.
    const auto near=frameScalar(coordinate(a,normal)),d=frameScalar(depth);std::array<FrameRatio,2> lower,width;
    lower[normal]=positive?near-d:near;width[normal]=d;lower[1-normal]=frameScalar(std::min(coordinate(a,1-normal),coordinate(b,1-normal)));
    width[1-normal]=frameScalar(std::max(coordinate(a,1-normal),coordinate(b,1-normal)))-lower[1-normal];
    for(int axis=0;axis<2;++axis)require(lower[axis]>=frameScalar(axis==0?box.lower.x:box.lower.y)
        &&lower[axis]+width[axis]<=frameScalar(axis==0?box.upper.x:box.upper.y),"[QuadraticCellCut] exact framed strip leaves its original donor");
    if(rz)require(lower[0]>=0,"[QuadraticCellCut] exact framed strip crosses negative radius");
    const auto frame=basis(graph.shape);const auto c=center(box);const auto tx=frameScalar(frame.tangent.x),ty=frameScalar(frame.tangent.y),
        nx=frameScalar(frame.normal.x),ny=frameScalar(frame.normal.y),determinant=nx*nx+ny*ny;
    require(determinant>0,"[QuadraticCellCut] invalid original framed strip determinant");
    struct Vertex{FrameRatio t,n;};std::array<Vertex,4> vertices;std::vector<FrameRatio> slabs;
    for(int i=0;i<4;++i){const auto x=lower[0]+((i==1||i==2)?width[0]:FrameRatio(0))-frameScalar(c.x);
        const auto y=lower[1]+(i>=2?width[1]:FrameRatio(0))-frameScalar(c.y);vertices[i]={tx*x+ty*y,nx*x+ny*y};slabs.push_back(vertices[i].t);}
    std::sort(slabs.begin(),slabs.end());slabs.erase(std::unique(slabs.begin(),slabs.end()),slabs.end());
    const auto offset=frameScalar(graph.offset),slope=frameScalar(graph.shape.slope),halfSecond=frameScalar(graph.shape.second)/2,
        twoPi=frameScalar(CartesianMeasure::twoPi),radiusOrigin=frameScalar(c.x);FrameRatio belowTotal=0,aboveTotal=0;
    for(std::size_t slab=1;slab<slabs.size();++slab){const auto lo=slabs[slab-1],hi=slabs[slab],tw=hi-lo;std::vector<FramePoly> boundary;
        for(int i=0;i<4;++i){auto first=vertices[i],last=vertices[(i+1)%4];if(last.t<first.t)std::swap(first,last);
            const auto dt=last.t-first.t;if(dt==0||lo<first.t||last.t<hi)continue;const auto dn=last.n-first.n;
            FramePoly edge{};edge[0]=first.n+(lo-first.t)*dn/dt;edge[1]=tw*dn/dt;boundary.push_back(std::move(edge));}
        require(boundary.size()==2,"[QuadraticCellCut] invalid exact framed strip slice");
        if(frameValue(boundary[0],FrameRatio(1)/2)>frameValue(boundary[1],FrameRatio(1)/2))std::swap(boundary[0],boundary[1]);
        FramePoly polynomial{};polynomial[0]=offset+slope*lo+halfSecond*lo*lo;
        polynomial[1]=(slope+2*halfSecond*lo)*tw;polynomial[2]=halfSecond*tw*tw;
        std::vector<FrameRatio> cuts{0,1};for(const auto& edge:boundary){FramePoly event{};for(int k=0;k<5;++k)event[k]=polynomial[k]-edge[k];frameRoots(event,cuts);}
        std::sort(cuts.begin(),cuts.end());cuts.erase(std::unique(cuts.begin(),cuts.end()),cuts.end());
        const auto density=[&](const FramePoly& upper,const FramePoly& lowerBound){FramePoly difference{};
            for(int k=0;k<5;++k)difference[k]=(upper[k]-lowerBound[k])/determinant;if(!rz)return difference;
            FramePoly radius{};radius[0]=radiusOrigin+tx*lo/determinant;radius[1]=tx*tw/determinant;
            for(int k=0;k<5;++k)radius[k]+=nx*(upper[k]+lowerBound[k])/(2*determinant);
            auto result=frameProduct(difference,radius);for(auto& q:result)q*=twoPi;return result;};
        FrameRatio belowIntegral=0,aboveIntegral=0;
        for(std::size_t span=1;span<cuts.size();++span){const auto left=cuts[span-1],right=cuts[span],at=(left+right)/2,g=frameValue(polynomial,at);
            auto split=polynomial;if(g<=frameValue(boundary[0],at))split=boundary[0];else if(g>=frameValue(boundary[1],at))split=boundary[1];
            const auto below=frameIntegral(density(split,boundary[0]),left,right),above=frameIntegral(density(boundary[1],split),left,right);
            require(below>=0&&above>=0,"[QuadraticCellCut] unresolved positive framed strip phase integral");belowIntegral+=below;aboveIntegral+=above;}
        require(belowIntegral+aboveIntegral==frameIntegral(density(boundary[1],boundary[0]),0,1),"[QuadraticCellCut] exact framed phase slice identity failed");
        belowTotal+=tw*belowIntegral;aboveTotal+=tw*aboveIntegral;
    }
    auto total=width[0]*width[1];if(rz)total*=twoPi*(lower[0]+width[0]/2);
    require(belowTotal+aboveTotal==total,"[QuadraticCellCut] exact framed strip partition identity failed");
    StripPhaseMeasure result;result.total=frameRounded(total);const auto below=frameRounded(belowTotal),above=frameRounded(aboveTotal);
    result.liquid=graph.shape.sign>0?below:above;result.gas=graph.shape.sign>0?above:below;return result;
}
Reconstruction reconstruct(const Rectangle& box,Shape shape,long double fraction,bool rz)
{
    validate(shape);require(std::isfinite(fraction) && fraction>=0 && fraction<=1,"[QuadraticCellCut] invalid fraction");
    const auto capacity=PlanarCellCut::capacity(box,rz);const auto c=center(box);
    const long double hx=(static_cast<long double>(box.upper.x)-box.lower.x)/2,hy=(static_cast<long double>(box.upper.y)-box.lower.y)/2;
    const auto frame=basis(shape);
    const auto ht=std::abs(frame.tangent.x)*hx+std::abs(frame.tangent.y)*hy,
        hn=std::abs(frame.normal.x)*hx+std::abs(frame.normal.y)*hy;
    const auto extent=hn+std::abs(shape.slope)*ht+std::abs(shape.second)*ht*ht/2;
    const std::vector<Point> polygon{{-hx,-hy},{hx,-hy},{hx,hy},{-hx,hy}};
    Reconstruction result{shape};long double lo=-extent,hi=extent;const auto target=fraction*capacity;
    for (int i=0;i<96;++i) {
        result.offset=lo+(hi-lo)/2;bool full=false;result.volume=polygonVolume(box,polygon,c,result,rz,&full);
        if(full)result.volume=capacity;
        if (result.volume==target || result.offset==lo || result.offset==hi) break;
        if ((result.volume<target)==(shape.sign>0)) lo=result.offset;else hi=result.offset;
    }
    result.residual=result.volume-target;
    require(std::abs(result.residual)<=PlanarCellCut::relativeVolumeResolution*capacity,
        "[QuadraticCellCut] unresolved curved volume constraint");return result;
}
std::optional<Reconstruction> projectEndpointVolume(const Rectangle& box,int axis,int sign,
    long double lowerHeight,long double upperHeight,long double liquidVolume,bool rz)
{
    validate({axis,sign,0,0});const auto capacity=PlanarCellCut::capacity(box,rz);const auto c=center(box);
    require(std::isfinite(lowerHeight) && std::isfinite(upperHeight) && std::isfinite(liquidVolume)
        && liquidVolume>=0 && liquidVolume<=capacity,"[QuadraticCellCut] invalid projection data");
    const long double width=axis?static_cast<long double>(box.upper.x)-box.lower.x:
        static_cast<long double>(box.upper.y)-box.lower.y;
    const long double nlo=axis?box.lower.y:box.lower.x,nhi=axis?box.upper.y:box.upper.x;
    if(!(liquidVolume>0 && liquidVolume<capacity) || lowerHeight<nlo || lowerHeight>nhi
        || upperHeight<nlo || upperHeight>nhi) return std::nullopt;
    const auto average=(lowerHeight+upperHeight)/2,slope=(upperHeight-lowerHeight)/width;
    // H(t)=average+slope*t+k*(t*t/width^2-1/4). Endpoints stay fixed.
    long double k=0;
    if(rz && axis==0) {
        const auto meanSquare=sign>0?nlo*nlo+liquidVolume/((CartesianMeasure::twoPi/2)*width):
            nhi*nhi-liquidVolume/((CartesianMeasure::twoPi/2)*width);
        const auto difference=(meanSquare-average*average)-slope*slope*width*width/12;
        const auto discriminant=100*average*average+120*difference;
        if(!(discriminant>0) || !(average>0)) return std::nullopt;
        // Physical root continuous with k=0, without subtracting close roots.
        k=-60*difference/(10*average+std::sqrt(discriminant));
    } else {
        const auto measure=rz?CartesianMeasure::twoPi*c.x*width:width;
        const auto mean=sign>0?nlo+liquidVolume/measure:nhi-liquidVolume/measure;
        const auto linearMean=average+(rz?slope*width*width/(12*c.x):0);
        k=6*(linearMean-mean);
    }
    const auto middle=average-k/4,second=2*k/(width*width);
    if(!std::isfinite(middle) || !std::isfinite(second)) return std::nullopt;
    if(second!=0) {
        const auto stationary=-slope/second;
        if(std::abs(stationary)<width/2) {
            const auto value=middle+slope*stationary+second*stationary*stationary/2;
            if(value<nlo || value>nhi) return std::nullopt;
        }
    }
    Reconstruction result{{axis,sign,slope,second},middle-coordinate(c,axis)};
    const long double hx=(static_cast<long double>(box.upper.x)-box.lower.x)/2,
        hy=(static_cast<long double>(box.upper.y)-box.lower.y)/2;
    result.volume=polygonVolume(box,{{-hx,-hy},{hx,-hy},{hx,hy},{-hx,hy}},c,result,rz);
    result.residual=result.volume-liquidVolume;
    if(std::abs(result.residual)>PlanarCellCut::relativeVolumeResolution*capacity) return std::nullopt;
    return result;
}
std::vector<std::array<long double,2>> intervals(const Rectangle& box,const Reconstruction& graph,Point begin,Point end,
    long double lower,long double upper)
{
    if(graph.exact){Detail::validateGraph(graph);const auto&e=*graph.exact;const auto a=QE::subtract(QE::vector(begin),e.origin),v=QE::subtract(QE::vector(end),QE::vector(begin));
        const QR t0=QE::dot(a,e.tangent)/e.determinant,t1=QE::dot(v,e.tangent)/e.determinant,n0=QE::dot(a,e.normal)/e.determinant,n1=QE::dot(v,e.normal)/e.determinant;
        const QP p{{n0-e.q[0]-e.q[1]*t0-e.q[2]*t0*t0,n1-e.q[1]*t1-2*e.q[2]*t0*t1,-e.q[2]*t1*t1,0,0}};
        std::vector<QR> cuts{QE::scalar(lower),QE::scalar(upper)};const auto rr=QE::roots(p,cuts[0],cuts[1]);cuts.insert(cuts.end(),rr.begin(),rr.end());std::sort(cuts.begin(),cuts.end());cuts.erase(std::unique(cuts.begin(),cuts.end()),cuts.end());
        std::vector<std::array<long double,2>> result;for(std::size_t i=1;i<cuts.size();++i)if(graph.shape.sign*QE::value(p,(cuts[i-1]+cuts[i])/2)<=0)result.push_back({QE::rounded(cuts[i-1]),QE::rounded(cuts[i])});return result;}

    validate(graph.shape);require(lower>=0 && upper<=1 && lower<=upper,"[QuadraticCellCut] invalid chord interval");
    const auto c=center(box);begin.x-=c.x;begin.y-=c.y;end.x-=c.x;end.y-=c.y;
    const auto s=graph.shape;const auto frame=basis(s);
    const auto t=begin.x*frame.tangent.x+begin.y*frame.tangent.y,
        dt=(end.x-begin.x)*frame.tangent.x+(end.y-begin.y)*frame.tangent.y;
    const auto n=begin.x*frame.normal.x+begin.y*frame.normal.y,
        dn=(end.x-begin.x)*frame.normal.x+(end.y-begin.y)*frame.normal.y;
    const auto a=-s.second*dt*dt/2,b=dn-s.slope*dt-s.second*t*dt,d=n-graph.offset-s.slope*t-s.second*t*t/2;
    std::vector<long double> cuts{lower,upper};roots(a,b,d,lower,upper,cuts);
    std::sort(cuts.begin(),cuts.end());cuts.erase(std::unique(cuts.begin(),cuts.end()),cuts.end());
    std::vector<std::array<long double,2>> result;
    for (std::size_t i=1;i<cuts.size();++i) {
        const auto mid=cuts[i-1]+(cuts[i]-cuts[i-1])/2;
        if (s.sign*((a*mid+b)*mid+d)<=0) result.push_back({{cuts[i-1],cuts[i]}});
    }
    return result;
}
}



