#include "Modules/Electric/ElectricSurfaceChargeTransport.h"
#include "Modules/Electric/Detail/ElectricSurfaceChargeProjection.h"
#include "Core/Geometry/PreparedUnitPolynomial.h"
#include "Core/Geometry/NativeOutwardRounding.h"
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <limits>
#include <map>
#include <numeric>
#include <optional>
#include <set>
#include <sstream>
#include <string>
#include <utility>

namespace LAMP::Electric
{
namespace Detail::SurfaceChargeProjection
{
namespace
{
using Polynomial=std::array<double,12>;
using Term=Geometry::UnitPolynomialTerm;
void require(bool condition,const char* message)
{if(!condition)throw ElectricSurfaceChargeTransferError(ElectricSurfaceChargeTransferFailure::NonconvexProjection,message);}
long double choose(unsigned n,unsigned k)
{long double r=1;for(unsigned j=1;j<=k;++j)r=r*(n+1-j)/j;return r;}
Polynomial bernsteinPower(unsigned i)
{
    Polynomial result{};for(unsigned k=i;k<=4;++k)result[k]=static_cast<double>(choose(4,i)*choose(4-i,k-i)*((k-i)%2?-1:1));
    return result;
}
Polynomial product(const Polynomial& a,const Polynomial& b)
{
    Polynomial result{};for(unsigned i=0;i<=4;++i)for(unsigned j=0;j<=4;++j)result[i+j]+=a[i]*b[j];return result;
}
void appendProduct(std::vector<Term>& terms,const Polynomial& polynomial,long double factor,long double weight)
{
    // Keep actual coordinate factors separate until Core exact expansion.
    for(unsigned part=0;weight!=0&&part<8;++part){const auto head=static_cast<double>(weight);
        require(std::isfinite(head)&&head!=0,"unrepresentable exact projection factor");
        terms.push_back({polynomial,factor,head});weight-=static_cast<long double>(head);}
    require(weight==0,"unresolved exact projection factor expansion");
}
using DistanceBasis=std::array<Polynomial,5>;
DistanceBasis distanceBasis()
{
    DistanceBasis b{};for(unsigned i=0;i<5;++i)b[i]=bernsteinPower(i);return b;
}
void appendGeometryTerms(std::vector<Term>& result,const std::array<Point,5>& controls,const DistanceBasis& b)
{
    for(unsigned i=0;i<5;++i)for(unsigned j=0;j<5;++j){const auto polynomial=product(b[i],b[j]);
        appendProduct(result,polynomial,controls[i].x,controls[j].x);
        appendProduct(result,polynomial,controls[i].y,controls[j].y);}
}
void appendQueryTerms(std::vector<Term>& result,const std::array<Point,5>& controls,Point p,const DistanceBasis& b)
{
    for(unsigned i=0;i<5;++i){auto polynomial=b[i];for(auto& coefficient:polynomial)coefficient*=-2;
        appendProduct(result,polynomial,p.x,controls[i].x);appendProduct(result,polynomial,p.y,controls[i].y);}
    Polynomial constant{};constant[0]=1;appendProduct(result,constant,p.x,p.x);appendProduct(result,constant,p.y,p.y);
}
std::vector<Term> geometryDistanceTerms(const std::array<Point,5>& controls)
{
    const auto b=distanceBasis();std::vector<Term> result;result.reserve(64);appendGeometryTerms(result,controls,b);return result;
}
std::vector<Term> queryDistanceTerms(const std::array<Point,5>& controls,Point p)
{
    const auto b=distanceBasis();std::vector<Term> result;result.reserve(16);appendQueryTerms(result,controls,p,b);return result;
}
std::vector<Term> distanceTerms(const std::array<Point,5>& controls,Point p)
{
    const auto b=distanceBasis();std::vector<Term> result;result.reserve(64);
    appendGeometryTerms(result,controls,b);appendQueryTerms(result,controls,p,b);return result;
}
std::vector<Term> derivativeTerms(std::vector<Term> terms)
{
    for(auto& term:terms){for(unsigned i=1;i<term.coefficients.size();++i)term.coefficients[i-1]=term.coefficients[i]*i;term.coefficients.back()=0;}return terms;
}
struct ExactDistanceBases
{
    Geometry::PreparedUnitPolynomial distance,derivative;
    explicit ExactDistanceBases(const std::array<Point,5>& controls):ExactDistanceBases(geometryDistanceTerms(controls)){}
private:
    explicit ExactDistanceBases(const std::vector<Term>& terms):distance(terms),derivative(derivativeTerms(terms)){}
};
Selection selectPrepared(const Geometry::PreparedUnitPolynomial& exactDistance,const Geometry::PreparedUnitPolynomial& exactDerivative)
{
    constexpr bool nativeBinary64=sizeof(long double)==sizeof(double)&&sizeof(double)==8
        &&std::numeric_limits<long double>::is_iec559&&std::numeric_limits<long double>::radix==2
        &&std::numeric_limits<long double>::digits==53&&std::numeric_limits<long double>::max_exponent==1024;
    const auto roots=nativeBinary64?exactDerivative.isolateNativeAbsolute():exactDerivative.isolate();
    require(!roots.identicallyZero,"constant-distance arc has no unique projection");
    std::vector<Geometry::UnitPolynomialRoots::Root> points{{0,0,0},{1,1,1}};
    for(const auto& root:roots.roots)if(!(root.lower==root.upper&&(root.lower==0||root.lower==1))){
        if constexpr(nativeBinary64)require(Geometry::nativeOutwardNeighbor(static_cast<long double>(root.upper)-root.lower,true)
            <=64*std::numeric_limits<long double>::epsilon(),"exact projection root exceeds native absolute resolution");
        points.push_back(root);}
    Selection result;result.interiorRoots=points.size()-2;result.candidates.reserve(points.size());
    for(const auto& root:points){const auto bounds=exactDistance.enclose(root.lower,root.upper);
        require(std::isfinite(bounds[0])&&std::isfinite(bounds[1])&&bounds[1]>=0,"invalid exact projection distance enclosure");
        const auto endpoint=root.lower==root.upper&&(root.lower==0||root.lower==1)?(root.lower==0?0:1):-1;
        result.candidates.push_back({root,{std::max(0.L,static_cast<long double>(bounds[0])),static_cast<long double>(bounds[1])},endpoint});}
    const auto best=std::min_element(result.candidates.begin(),result.candidates.end(),
        [](const auto& a,const auto& b){return a.squaredDistance[1]<b.squaredDistance[1];});
    result.roundedWinner=static_cast<std::size_t>(best-result.candidates.begin());
    // This is the COMPLETE exact D' root inventory. With one nonendpoint
    // candidate and strict endpoint signs, D' cannot change sign
    // anywhere else: D decreases to this root and increases after it.
    // The common normal-distance floor need not discriminate D(root),D(0).
    if(result.interiorRoots==1){const auto& root=result.candidates[2].parameter;
        if(root.lower>=0&&root.upper<=1&&root.point>=root.lower&&root.point<=root.upper){
            const auto left=exactDerivative.enclose(0,0),right=exactDerivative.enclose(1,1);
            if(left[1]<0&&right[0]>0){result.selected=2;result.uniqueInteriorCertificate=true;
                result.minimum=MinimumCertificate::StrictInterior;return result;}}}
    // Preserve the original certificate for every remaining stationary
    // candidate and both endpoints. Equal/overlapping minima still reject.
    for(std::size_t k=0;k<result.candidates.size();++k)if(k!=result.roundedWinner
        &&!(best->squaredDistance[1]<result.candidates[k].squaredDistance[0])){
        result.overlap=k;return result;}
    result.selected=result.roundedWinner;
    result.minimum=best->endpoint<0?MinimumCertificate::StrictInterior:MinimumCertificate::StrictEndpoint;return result;
}
}
Selection select(const std::array<Point,5>& controls,Point p)
{
    const auto distance=distanceTerms(controls,p),derivative=derivativeTerms(distance);
    const Geometry::PreparedUnitPolynomial exactDistance(distance),exactDerivative(derivative);
    return selectPrepared(exactDistance,exactDerivative);
}
static Selection selectCached(const std::array<Point,5>& controls,Point p,const ExactDistanceBases& base)
{
    const auto suffix=queryDistanceTerms(controls,p),derivative=derivativeTerms(suffix);
    const Geometry::PreparedUnitPolynomial exactDistance(base.distance,suffix),exactDerivative(base.derivative,derivative);
    return selectPrepared(exactDistance,exactDerivative);
}
std::optional<std::size_t> selectLocalMinimum(const std::vector<LocalMinimum>& candidates)
{
    const auto samePoint=[](Point a,Point b){return a.x==b.x&&a.y==b.y;};
    const auto sameEndpoint=[&](const LocalMinimum& a,const LocalMinimum& b){
        return a.minimum==MinimumCertificate::StrictEndpoint&&b.minimum==MinimumCertificate::StrictEndpoint
            &&a.endpoint>=0&&a.endpoint<2&&b.endpoint>=0&&b.endpoint<2
            &&a.nodes[a.endpoint]==b.nodes[b.endpoint]&&samePoint(a.endpoints[a.endpoint],b.endpoints[b.endpoint]);};
    const auto strictSharedEndpoint=[&](const LocalMinimum& a,const LocalMinimum& b){
        if(a.minimum!=MinimumCertificate::StrictInterior||a.endpoint!=-1
            ||b.minimum!=MinimumCertificate::StrictEndpoint||b.endpoint<0||b.endpoint>1)return false;
        // A strictly certified interior global minimum beats BOTH of its
        // actual endpoints. At a shared owned node, B's endpoint is exactly
        // that same physical point: D_A(min)<D_A(node)=D_B(endpoint).
        for(unsigned side=0;side<2;++side)if(a.nodes[side]==b.nodes[b.endpoint]
            &&samePoint(a.endpoints[side],b.endpoints[b.endpoint]))return true;
        return false;};
    for(std::size_t i=0;i<candidates.size();++i){bool certified=true;
        for(std::size_t j=0;j<candidates.size();++j)if(i!=j&&!sameEndpoint(candidates[i],candidates[j])
            &&!(candidates[i].squaredDistance[1]<candidates[j].squaredDistance[0])
            &&!strictSharedEndpoint(candidates[i],candidates[j])){certified=false;break;}
        if(certified){auto selected=i;
            for(std::size_t j=0;j<candidates.size();++j)if(candidates[j].arc<candidates[selected].arc
                &&sameEndpoint(candidates[i],candidates[j]))selected=j;
            return selected;}}
    return {};
}
}
namespace
{
using Point=Geometry::PlanarCellCut::Point;
using Smooth=Interface::VofSmoothGeometry;
using Arc=Interface::VofSmoothArc;
using Failure=ElectricSurfaceChargeTransferFailure;
using Clock=LAMP::TimeLayer;
using ExactDistanceBases=Detail::SurfaceChargeProjection::ExactDistanceBases;
void require(bool yes,Failure reason,const char* message)
{if(!yes)throw ElectricSurfaceChargeTransferError(reason,message);}
bool sameClock(Clock a,Clock b)
{return a.step==b.step&&a.begin==b.begin&&a.end==b.end&&a.position==b.position;}
Point sub(Point a,Point b){return {a.x-b.x,a.y-b.y};}
long double dot(Point a,Point b){return a.x*b.x+a.y*b.y;}
long double norm(Point a){return std::hypot(a.x,a.y);}
std::array<long double,3> basis(long double t)
{const auto z=2*t-1;return {1,z,(3*z*z-1)/2};}
constexpr std::array<long double,8> gaussX{{
    -.960289856497536231683560868569L,-.796666477413626739591553936476L,
    -.525532409916328985817739049189L,-.183434642495649804939476142360L,
     .183434642495649804939476142360L, .525532409916328985817739049189L,
     .796666477413626739591553936476L, .960289856497536231683560868569L}};
constexpr std::array<long double,8> gaussW{{
    .101228536290376259152531354310L,.222381034453374470544355994426L,
    .313706645877887287337962201987L,.362683783378361982965150449277L,
    .362683783378361982965150449277L,.313706645877887287337962201987L,
    .222381034453374470544355994426L,.101228536290376259152531354310L}};
struct Sum
{
    long double value=0,correction=0;
    void add(long double x){const auto y=x-correction,t=value+y;correction=(t-value)-y;value=t;}
};
using WallClock=std::chrono::steady_clock;
struct PassiveTimer
{
    double& seconds;WallClock::time_point start=WallClock::now();
    ~PassiveTimer(){seconds+=std::chrono::duration<double>(WallClock::now()-start).count();}
};
struct Topology
{
    std::vector<std::set<std::size_t>> adjacent;
    std::vector<Point> endpoints;
};
Topology topology(const Smooth& geometry)
{
    const auto& arcs=geometry.arcs();const auto& nodes=geometry.nodes();
    require(!arcs.empty(),Failure::UnsupportedTopology,"empty surface component");
    std::vector<std::vector<std::size_t>> incident(nodes.size());
    for(std::size_t i=0;i<arcs.size();++i)for(auto n:arcs[i].nodes){
        require(n<nodes.size(),Failure::UnsupportedTopology,"foreign arc node");incident[n].push_back(i);}
    Topology result;result.adjacent.resize(arcs.size());
    for(std::size_t n=0;n<nodes.size();++n){
        require(incident[n].size()==1||incident[n].size()==2,Failure::UnsupportedTopology,"branched surface component");
        if(incident[n].size()==1)result.endpoints.push_back(nodes[n].point);
        else {result.adjacent[incident[n][0]].insert(incident[n][1]);result.adjacent[incident[n][1]].insert(incident[n][0]);}}
    std::set<std::size_t> seen{0};std::vector<std::size_t> pending{0};
    while(!pending.empty()){const auto i=pending.back();pending.pop_back();
        for(auto j:result.adjacent[i])if(seen.insert(j).second)pending.push_back(j);}
    require(seen.size()==arcs.size()&&(result.endpoints.empty()||result.endpoints.size()==2),
        Failure::UnsupportedTopology,"only one nonbranching surface component is supported");
    return result;
}
struct Grid
{
    long double x0=0,y0=0,hx=0,hy=0;
    std::map<std::pair<int,int>,int> cells;
    std::map<std::pair<int,int>,std::size_t> arcs;
    std::vector<std::array<long double,4>> velocity; // W,S,E,N.
    std::pair<int,int> index(Point p) const
    {return {static_cast<int>(std::floor((p.x-x0)/hx)),static_cast<int>(std::floor((p.y-y0)/hy))};}
};
Grid gridFor(const QuadtreeMesh& mesh,const Smooth& source,const Smooth& target,
    const Interface::GeometricVofTransport* transport,bool half)
{
    Grid result;const auto& domain=source.source().domain();
    require(domain.uniformLeaves()&&!domain.axisymmetric(),Failure::UnsupportedMesh,"surface transfer requires uniform XY leaves");
    require(domain.meshIdentity()==target.source().domain().meshIdentity(),Failure::UnsupportedMesh,"AMR or a rebuilt mesh owner is unsupported");
    const auto& indices=domain.cellIndexByRow();
    require(!indices.empty(),Failure::UnsupportedMesh,"empty native domain");
    const auto& first=source.cells().at(indices.front()).rectangle;
    result.hx=first.upper.x-first.lower.x;result.hy=first.upper.y-first.lower.y;
    result.x0=first.lower.x;result.y0=first.lower.y;
    for(int ci:indices){const auto& r=source.cells().at(ci).rectangle;
        result.x0=std::min(result.x0,static_cast<long double>(r.lower.x));result.y0=std::min(result.y0,static_cast<long double>(r.lower.y));}
    require(result.hx>0&&result.hy>0,Failure::UnsupportedMesh,"invalid native spacing");
    result.velocity.resize(source.cells().size());std::vector<std::array<unsigned,4>> counts(source.cells().size());
    for(int ci:indices){const auto& r=source.cells().at(ci).rectangle;
        require(r.upper.x-r.lower.x==result.hx&&r.upper.y-r.lower.y==result.hy,Failure::UnsupportedMesh,"mixed native spacing");
        const Point center{(r.lower.x+r.upper.x)/2,(r.lower.y+r.upper.y)/2};
        require(result.cells.emplace(result.index(center),ci).second,Failure::UnsupportedMesh,"overlapping native cells");}
    for(std::size_t a=0;a<target.arcs().size();++a){const auto ci=target.arcs()[a].nativeCell;const auto& r=target.cells().at(ci).rectangle;
        const Point center{(r.lower.x+r.upper.x)/2,(r.lower.y+r.upper.y)/2};
        require(result.arcs.emplace(result.index(center),a).second,Failure::UnsupportedTopology,"multiple arcs per native cell");}
    if(!transport)return result;
    const auto& slab=transport->advector();
    require(!slab.boundary()||slab.boundary()->allClosed(),Failure::UnsupportedBoundary,"open surface transport is unsupported");
    for(int row=0;row<domain.nativeFaceCount();++row){const auto f=domain.face(mesh,row);
        const auto a=static_cast<long double>(slab.beginNormal()[row]),b=static_cast<long double>(slab.endNormal()[row]);
        const auto mean=half?(3*a+b)/4:(a+b)/2;
        if(f.negativeCell>=0){const auto side=f.horizontal?3:2;result.velocity[f.negativeCell][side]=mean;++counts[f.negativeCell][side];}
        if(f.positiveCell>=0){const auto side=f.horizontal?1:0;result.velocity[f.positiveCell][side]=mean;++counts[f.positiveCell][side];}}
    for(int ci:indices)for(auto n:counts[ci])require(n==1,Failure::UnsupportedMesh,"split native faces are unsupported");
    return result;
}
bool samePoint(Point a,Point b){return a.x==b.x&&a.y==b.y;}
// An exact oriented PARAMETER identity, not proximity or a geometric-set
// match. Each target arc uses the same t and therefore the same P2 basis.
std::optional<std::vector<std::size_t>> identityArcs(const Smooth& source,const Smooth& target)
{
    if(source.arcs().size()!=target.arcs().size()||source.nodes().size()!=target.nodes().size()
        ||source.source().domain().cellIndexByRow()!=target.source().domain().cellIndexByRow())return std::nullopt;
    std::map<int,std::size_t> targetByCell;
    for(std::size_t a=0;a<target.arcs().size();++a)
        if(!targetByCell.emplace(target.arcs()[a].nativeCell,a).second)return std::nullopt;
    const auto missing=std::numeric_limits<std::size_t>::max();
    std::vector<std::size_t> result(source.arcs().size()),nodes(source.nodes().size(),missing),inverse(target.nodes().size(),missing);
    for(std::size_t a=0;a<source.arcs().size();++a){const auto& old=source.arcs()[a];
        const auto found=targetByCell.find(old.nativeCell);if(found==targetByCell.end())return std::nullopt;
        const auto& next=target.arcs()[found->second];
        if(old.graphAxis!=next.graphAxis)return std::nullopt;
        for(unsigned k=0;k<5;++k)if(!samePoint(old.controlPoints[k],next.controlPoints[k])
            ||!samePoint(old.power[k],next.power[k]))return std::nullopt;
        for(unsigned side=0;side<2;++side){const auto p=old.nodes[side],q=next.nodes[side];
            const auto& x=source.nodes()[p];const auto& y=target.nodes()[q];
            if(x.nativeRow!=y.nativeRow||!samePoint(x.point,y.point)||!samePoint(x.tangent,y.tangent)
                ||(nodes[p]!=missing&&nodes[p]!=q)||(inverse[q]!=missing&&inverse[q]!=p))return std::nullopt;
            nodes[p]=q;inverse[q]=p;}
        result[a]=found->second;}
    if(std::find(nodes.begin(),nodes.end(),missing)!=nodes.end()
        ||std::find(inverse.begin(),inverse.end(),missing)!=inverse.end())return std::nullopt;
    return result;
}
bool zeroNativeStage(const Interface::GeometricVofTransport& transport)
{
    const auto& slab=transport.advector();
    for(const auto* samples:{&slab.beginNormal(),&slab.endNormal()})
        for(const auto value:*samples)if(value!=0)return false;
    // LocalPlic's actual velocity unknowns are the native normal means.
    // It deliberately retains no reconstructed vector samples. Zero at BOTH
    // endpoint normal snapshots proves zero native material advector; do not
    // query or invent an unavailable tangential reconstruction.
    if(transport.source().geometryKind()!=Interface::GeometricVofGeometryKind::LocalPlic){
        for(const auto* samples:{&slab.beginVelocity(),&slab.endVelocity()})
            for(const auto value:*samples)if(value.x!=0||value.y!=0)return false;}
    // Never infer zero physical motion from a cancelling full-interval mean,
    // a net phase flux, or unchanged target geometry alone.
    const auto zeroFlux=[](const Interface::GeometricVofIntegratedFlux& f){
        return f.total==0&&f.positiveTotal==0&&f.negativeTotal==0&&f.rawLiquid==0
            &&f.liquid==0&&f.gas==0&&f.mass==0&&f.positiveLiquid==0&&f.negativeLiquid==0
            &&f.positiveGas==0&&f.negativeGas==0;};
    for(const auto& face:transport.faces())require(zeroFlux(face.half)&&zeroFlux(face.full),
        Failure::ForeignSource,"zero native velocity has nonzero actual phase flux");
    require(transport.report().maximumCompressionResidual==0&&transport.report().maximumTotalDivergence==0,
        Failure::ForeignSource,"zero native velocity has a nonzero actual WY/divergence ledger");
    const auto& old=transport.source().plicCells();
    for(const auto* stage:{&transport.midpoint(),&transport.endpoint()}){
        require(stage->plicCells().size()==old.size(),Failure::ForeignSource,"zero native stage changed stock layout");
        for(std::size_t i=0;i<old.size();++i)require(stage->plicCells()[i].inventory==old[i].inventory,
            Failure::ForeignSource,"zero native stage changed independent phase stocks");}
    return true;
}
// Exact characteristic of each frozen cellwise linear normal velocity.
// This shares the real split order/temporal means, but is a surface material
// approximation, not an exact inverse of the VOF strip/WY stock operator.
Point directional(Point p,unsigned axis,long double dt,const Grid& grid)
{
    auto cell=grid.index(p);long double remaining=dt;
    if(!grid.cells.count(cell)&&p.x==grid.x0+cell.first*grid.hx&&grid.cells.count({cell.first-1,cell.second}))--cell.first;
    if(!grid.cells.count(cell)&&p.y==grid.y0+cell.second*grid.hy&&grid.cells.count({cell.first,cell.second-1}))--cell.second;
    const auto h=axis?grid.hy:grid.hx,origin=axis?grid.y0:grid.x0;
    for(unsigned crossing=0;crossing<8&&remaining>0;++crossing){
        auto found=grid.cells.find(cell);require(found!=grid.cells.end(),Failure::UnsupportedBoundary,"material left the closed native domain");
        const auto& v=grid.velocity[found->second];const auto lower=origin+(axis?cell.second:cell.first)*h;
        const auto lo=v[axis?1:0],hi=v[axis?3:2],slope=(hi-lo)/h;
        auto& x=axis?p.y:p.x;const auto speed=lo+slope*(x-lower);
        if(speed==0)return p;
        const auto end=x+(slope==0?speed*remaining:speed*std::expm1(slope*remaining)/slope);
        require(std::isfinite(end),Failure::ExcessiveDisplacement,"nonfinite directional characteristic");
        if(end>=lower&&end<=lower+h){x=end;return p;}
        const auto boundary=speed>0?lower+h:lower;
        long double hit;
        if(slope==0)hit=(boundary-x)/speed;
        else {const auto ratio=slope*(boundary-x)/speed;
            require(ratio>-1,Failure::ExcessiveDisplacement,"invalid characteristic face crossing");hit=std::log1p(ratio)/slope;}
        require(std::isfinite(hit)&&hit>=0&&hit<=remaining,Failure::ExcessiveDisplacement,"invalid characteristic crossing time");
        remaining-=hit;x=boundary;
        if(axis)cell.second+=speed>0?1:-1;else cell.first+=speed>0?1:-1;
    }
    require(remaining==0,Failure::ExcessiveDisplacement,"directional characteristic crossed too many leaves");return p;
}
long double choose(unsigned n,unsigned k)
{long double r=1;for(unsigned j=1;j<=k;++j)r=r*(n+1-j)/j;return r;}
// Outward intervals keep the inexpensive Bernstein predicate a certificate.
// The exact Core filter below is for the same Euclidean map, not another map.
struct Bounds{long double lower=0,upper=0;};
long double down(long double x){return Geometry::nativeOutwardNeighbor(x,false);}
long double up(long double x){return Geometry::nativeOutwardNeighbor(x,true);}
Bounds exact(long double x){return {x,x};}
Bounds plus(Bounds a,Bounds b){return {down(a.lower+b.lower),up(a.upper+b.upper)};}
Bounds minus(Bounds a,Bounds b){return {down(a.lower-b.upper),up(a.upper-b.lower)};}
Bounds times(Bounds a,Bounds b)
{
    const std::array<long double,4> v{{a.lower*b.lower,a.lower*b.upper,a.upper*b.lower,a.upper*b.upper}};
    return {down(*std::min_element(v.begin(),v.end())),up(*std::max_element(v.begin(),v.end()))};
}
Bounds dividePositive(Bounds a,Bounds b)
{return times(a,{down(1/b.upper),up(1/b.lower)});}
Bounds square(Bounds a)
{
    const auto hi=up(std::max(a.lower*a.lower,a.upper*a.upper));
    return {a.lower<=0&&a.upper>=0?0:std::max(0.L,down(std::min(a.lower*a.lower,a.upper*a.upper))),hi};
}
struct IntervalPoint{Bounds x,y;};
IntervalPoint exact(Point p){return {exact(p.x),exact(p.y)};}
IntervalPoint plus(IntervalPoint a,IntervalPoint b){return {plus(a.x,b.x),plus(a.y,b.y)};}
IntervalPoint minus(IntervalPoint a,IntervalPoint b){return {minus(a.x,b.x),minus(a.y,b.y)};}
IntervalPoint times(IntervalPoint a,Bounds b){return {times(a.x,b),times(a.y,b)};}
Bounds inner(IntervalPoint a,IntervalPoint b){return plus(times(a.x,b.x),times(a.y,b.y));}
template<std::size_t N> IntervalPoint bezier(std::array<IntervalPoint,N> c,Bounds t)
{
    if(t.lower==0&&t.upper==0)return c.front();
    if(t.lower==1&&t.upper==1)return c.back();
    const auto oneMinus=minus(exact(1),t);
    for(std::size_t n=N-1;n>0;--n)for(std::size_t i=0;i<n;++i)c[i]=plus(times(c[i],oneMinus),times(c[i+1],t));
    return c[0];
}
struct NativeProductWeights
{
    std::array<Bounds,16> derivativeProduct;
    std::array<Bounds,15> curvatureProduct;
};
NativeProductWeights prepareNativeProductWeights()
{
    NativeProductWeights result;
    for(unsigned i=0;i<4;++i)for(unsigned j=0;j<4;++j)
        result.derivativeProduct[4*i+j]=dividePositive(exact(choose(3,i)*choose(3,j)),exact(choose(6,i+j)));
    for(unsigned i=0;i<5;++i)for(unsigned j=0;j<3;++j)
        result.curvatureProduct[3*i+j]=dividePositive(exact(choose(4,i)*choose(2,j)),exact(choose(6,i+j)));
    return result;
}
struct NativeArcConstants
{
    std::array<IntervalPoint,5> controls;
    std::array<IntervalPoint,4> derivative;
    std::array<IntervalPoint,3> secondDerivative;
    std::array<Bounds,7> hessianPrefix{};
};
NativeArcConstants prepareNativeArcConstants(const Arc& arc,const NativeProductWeights& weights)
{
    NativeArcConstants result;
    for(unsigned i=0;i<5;++i)result.controls[i]=exact(arc.controlPoints[i]);
    for(unsigned i=0;i<4;++i)result.derivative[i]=times(minus(result.controls[i+1],result.controls[i]),exact(4));
    for(unsigned i=0;i<3;++i)result.secondDerivative[i]=times(minus(result.derivative[i+1],result.derivative[i]),exact(3));
    // This is the completed ORIGINAL first-loop prefix, accumulated from
    // zero in the same i/j order. Query-dependent terms append to its copy.
    for(unsigned i=0;i<4;++i)for(unsigned j=0;j<4;++j){const auto weight=weights.derivativeProduct[4*i+j];
        result.hessianPrefix[i+j]=plus(result.hessianPrefix[i+j],times(weight,inner(result.derivative[i],result.derivative[j])));}
    return result;
}
enum class CoreProjectionReason { None,Hessian,EndpointSign,Nonfinite,RootWidth,ForcedDistance };
struct Projection
{
    std::size_t arc=0;long double t=0,distance=0;
    Bounds squaredDistance{};bool exactCertificate=false;int endpoint=-1;Bounds parameter{0,1};
    bool usedCoreRoots=false;
    Detail::SurfaceChargeProjection::MinimumCertificate minimum=Detail::SurfaceChargeProjection::MinimumCertificate::None;
    CoreProjectionReason coreReason=CoreProjectionReason::None;
    bool preparedGeometryBase=false;
};
Projection enclosedProjection(const Arc& arc,std::size_t id,Point p,Bounds parameter,bool certifiedExact=false)
{
    std::array<IntervalPoint,5> relativeControls{};
    for(unsigned i=0;i<5;++i)relativeControls[i]=minus(exact(arc.controlPoints[i]),exact(p));
    const auto delta=bezier(relativeControls,parameter);
    auto distance=plus(square(delta.x),square(delta.y));distance.lower=std::max(0.L,distance.lower);
    const auto endpoint=parameter.lower==parameter.upper&&(parameter.lower==0||parameter.lower==1)?(parameter.lower==0?0:1):-1;
    const auto t=(parameter.lower+parameter.upper)/2;
    const auto nominalDistance=norm(sub(arc.sample(t),p));
    // Only reached after the native global strict-Hessian certificate and
    // its endpoint signs/root enclosure have established a strict minimum.
    using Certificate=Detail::SurfaceChargeProjection::MinimumCertificate;
    return {id,t,up(std::max(std::sqrt(std::max(0.L,distance.upper)),nominalDistance)),distance,certifiedExact,endpoint,parameter,
        false,endpoint<0?Certificate::StrictInterior:Certificate::StrictEndpoint};
}
Projection projectArcExact(const Arc& arc,std::size_t id,Point p,
    std::optional<ExactDistanceBases>& base,CoreProjectionReason reason)
{
    try{
        const bool preparedNow=!base;if(preparedNow)base.emplace(arc.controlPoints);
        const auto selection=Detail::SurfaceChargeProjection::selectCached(arc.controlPoints,p,*base);
        if(!selection.selected){
            std::ostringstream context;context<<std::setprecision(std::numeric_limits<long double>::max_digits10)
                <<"exact arc minimum distance is not separated arc="<<id<<" native_cell="<<arc.nativeCell
                <<" query=("<<p.x<<','<<p.y<<") graph_axis="<<arc.graphAxis<<" nodes=("<<arc.nodes[0]<<','<<arc.nodes[1]<<')'
                <<" winner="<<selection.roundedWinner<<" overlap="<<*selection.overlap<<" controls=";
            for(const auto control:arc.controlPoints)context<<'('<<control.x<<','<<control.y<<')';
            for(std::size_t k=0;k<selection.candidates.size();++k){const auto& candidate=selection.candidates[k];
                const auto point=arc.sample(candidate.parameter.point);const auto delta=sub(point,p);
                context<<" candidate=(index="<<k<<",endpoint="<<candidate.endpoint<<",t="<<candidate.parameter.point
                    <<",parameter=["<<candidate.parameter.lower<<','<<candidate.parameter.upper<<']'
                    <<",squared_distance=["<<candidate.squaredDistance[0]<<','<<candidate.squaredDistance[1]<<']'
                    <<",nominal_point=("<<point.x<<','<<point.y<<"),nominal_squared_distance="<<dot(delta,delta)<<')';}
            throw ElectricSurfaceChargeTransferError(Failure::AmbiguousProjection,context.str().c_str());}
        const auto& selected=selection.candidates[*selection.selected];const auto& parameter=selected.parameter;
        return {id,parameter.point,up(std::sqrt(selected.squaredDistance[1])),
            {selected.squaredDistance[0],selected.squaredDistance[1]},true,selected.endpoint,{parameter.lower,parameter.upper},true,selection.minimum,reason,preparedNow};
    }catch(const ElectricSurfaceChargeTransferError&){throw;}
    catch(const std::exception& error){
        std::ostringstream context;context<<std::setprecision(std::numeric_limits<long double>::max_digits10)
            <<"bounded exact closest-point root certificate failed original_reason="<<error.what()
            <<" arc="<<id<<" native_cell="<<arc.nativeCell<<" query=("<<p.x<<','<<p.y<<") controls=";
        for(const auto control:arc.controlPoints)context<<'('<<control.x<<','<<control.y<<')';
        throw ElectricSurfaceChargeTransferError(Failure::NonconvexProjection,context.str().c_str());}
}
Projection projectArc(const Arc& arc,std::size_t id,Point p,
    std::vector<std::optional<NativeArcConstants>>& arcCache,std::optional<NativeProductWeights>& weightCache,
    std::optional<ExactDistanceBases>& base)
{
    for(unsigned endpoint:{0U,4U})if(p.x==arc.controlPoints[endpoint].x&&p.y==arc.controlPoints[endpoint].y)
        return {id,endpoint==0?0.L:1.L,0,{0,0},true,endpoint==0?0:1,exact(endpoint==0?0.L:1.L),false,
            Detail::SurfaceChargeProjection::MinimumCertificate::StrictEndpoint}; // Strict graph injectivity is a G1 prerequisite.
    const auto exactProjection=[&](CoreProjectionReason reason){return projectArcExact(arc,id,p,base,reason);};
    if(!weightCache)weightCache=prepareNativeProductWeights();
    if(!arcCache[id])arcCache[id]=prepareNativeArcConstants(arc,*weightCache);
    const auto& prepared=*arcCache[id];const auto& controls=prepared.controls;
    const auto& d=prepared.derivative;const auto& dd=prepared.secondDerivative;
    std::array<IntervalPoint,5> relativeControls{};
    for(unsigned i=0;i<5;++i)relativeControls[i]=minus(controls[i],exact(p));
    auto coefficients=prepared.hessianPrefix;
    for(unsigned i=0;i<5;++i)for(unsigned j=0;j<3;++j){const auto weight=weightCache->curvatureProduct[3*i+j];
        coefficients[i+j]=plus(coefficients[i+j],times(weight,inner(minus(controls[i],exact(p)),dd[j])));}
    Bounds hessian{std::numeric_limits<long double>::infinity(),0};
    for(const auto coefficient:coefficients){if(!std::isfinite(coefficient.lower)||!std::isfinite(coefficient.upper))return exactProjection(CoreProjectionReason::Nonfinite);
        if(coefficient.lower<=0)return exactProjection(CoreProjectionReason::Hessian);
        hessian.lower=std::min(hessian.lower,coefficient.lower);hessian.upper=std::max(hessian.upper,coefficient.upper);}
    const auto derivative=[&](long double t){return inner(bezier(relativeControls,exact(t)),bezier(d,exact(t)));};
    const auto left=derivative(0),right=derivative(1);
    if(left.lower>=0)return enclosedProjection(arc,id,p,exact(0));
    if(right.upper<=0)return enclosedProjection(arc,id,p,exact(1));
    if(!(left.upper<0&&right.lower>0))return exactProjection(
        !std::isfinite(left.lower)||!std::isfinite(left.upper)||!std::isfinite(right.lower)||!std::isfinite(right.upper)
            ?CoreProjectionReason::Nonfinite:CoreProjectionReason::EndpointSign);
    Bounds parameter{0,1};
    for(unsigned iteration=0;iteration<80;++iteration){const auto mid=(parameter.lower+parameter.upper)/2;
        if(mid==parameter.lower||mid==parameter.upper)break;
        const auto value=derivative(mid);const auto old=parameter;
        if(!std::isfinite(value.lower)||!std::isfinite(value.upper))return exactProjection(CoreProjectionReason::Nonfinite);
        if(value.upper<0)parameter.lower=mid;
        else if(value.lower>0)parameter.upper=mid;
        const auto corrected=minus(exact(mid),dividePositive(value,hessian));
        if(!std::isfinite(corrected.lower)||!std::isfinite(corrected.upper))return exactProjection(CoreProjectionReason::Nonfinite);
        parameter.lower=std::max(parameter.lower,corrected.lower);parameter.upper=std::min(parameter.upper,corrected.upper);
        require(parameter.lower<=parameter.upper,Failure::NonconvexProjection,"inconsistent certified projection root bracket");
        if(parameter.lower==old.lower&&parameter.upper==old.upper)break;}
    // A root-containing interval alone cannot justify using its midpoint as
    // the material coordinate. Require a machine-scale parameter enclosure;
    // a wider/stalled native interval goes to the same exact Core filter.
    if(up(parameter.upper-parameter.lower)>64*std::numeric_limits<long double>::epsilon())return exactProjection(CoreProjectionReason::RootWidth);
    return enclosedProjection(arc,id,p,parameter);
}
long double aabbDistanceLower(const Arc& arc,Point p)
{
    Point low=arc.controlPoints[0],high=low;
    for(const auto c:arc.controlPoints){low.x=std::min(low.x,c.x);low.y=std::min(low.y,c.y);high.x=std::max(high.x,c.x);high.y=std::max(high.y,c.y);}
    const auto dx=p.x<low.x?std::max(0.L,down(low.x-p.x)):(p.x>high.x?std::max(0.L,down(p.x-high.x)):0.L);
    const auto dy=p.y<low.y?std::max(0.L,down(low.y-p.y)):(p.y>high.y?std::max(0.L,down(p.y-high.y)):0.L);
    return std::max(0.L,down(down(dx*dx)+down(dy*dy)));
}
bool sameNode(const Projection& a,const Projection& b,const Smooth& geometry)
{
    if(a.endpoint<0||b.endpoint<0)return false;
    return geometry.arcs()[a.arc].nodes[a.endpoint]==geometry.arcs()[b.arc].nodes[b.endpoint];
}
using Controls=std::array<IntervalPoint,5>;
Controls intervalControls(const Arc& arc)
{Controls result{};for(unsigned i=0;i<5;++i)result[i]=exact(arc.controlPoints[i]);return result;}
long double controlPolygonLengthUpper(const Controls& controls)
{
    long double result=0;
    for(unsigned i=1;i<5;++i){const auto edge=minus(controls[i],controls[i-1]);const auto squared=plus(square(edge.x),square(edge.y));
        result=up(result+up(std::sqrt(std::max(0.L,squared.upper))));}
    return result;
}
long double partialLengthUpper(const Arc& arc,const Projection& projection,unsigned endpoint)
{
    if(projection.endpoint==static_cast<int>(endpoint))return 0;
    // The outer root bound includes the whole possible endpoint-to-projection
    // subarc. Its actual Bezier control polygon bounds physical curve length.
    const auto t=endpoint==0?projection.parameter.upper:projection.parameter.lower;
    require(t>=0&&t<=1,Failure::NonconvexProjection,"invalid local-chain projection parameter enclosure");
    auto row=intervalControls(arc);Controls left{},right{};left[0]=row[0];right[4]=row[4];
    const auto parameter=exact(t),complement=minus(exact(1),parameter);
    for(unsigned level=1;level<=4;++level){for(unsigned i=0;i<=4-level;++i)row[i]=plus(times(row[i],complement),times(row[i+1],parameter));
        left[level]=row[0];right[4-level]=row[4-level];}
    return controlPolygonLengthUpper(endpoint==0?left:right);
}
bool sameLocalContinuousChain(const Projection& best,const Projection& other,const Smooth& geometry,
    const Topology& graph,const std::set<std::size_t>& local,long double radius)
{
    // This qualifies a bounded continuous local patch only. It grants no
    // component-wide exemption and does not enlarge the material radius R.
    const auto budget=2*radius;
    for(unsigned direction=0;direction<2;++direction){auto length=partialLengthUpper(geometry.arcs()[best.arc],best,direction);
        if(length>budget)continue;
        auto current=best.arc;auto node=geometry.arcs()[current].nodes[direction];std::set<std::size_t> visited{current};
        for(std::size_t walked=0;walked<local.size();++walked){std::optional<std::size_t> next;
            for(const auto adjacent:graph.adjacent[current]){const auto& arc=geometry.arcs()[adjacent];
                if(arc.nodes[0]==node||arc.nodes[1]==node){require(!next,Failure::UnsupportedTopology,"branched local-chain node");next=adjacent;}}
            if(!next||!local.count(*next)||!visited.insert(*next).second)break;
            const auto& arc=geometry.arcs()[*next];const unsigned enter=arc.nodes[0]==node?0:1;
            if(*next==other.arc){const auto bound=up(length+partialLengthUpper(arc,other,enter));if(bound<=budget)return true;break;}
            length=up(length+controlPolygonLengthUpper(intervalControls(arc)));if(length>budget)break;
            current=*next;node=arc.nodes[1-enter];}
    }
    return false;
}
std::array<long double,3> densityCoefficients(const Arc& arc,const std::array<double,3>& moments)
{
    std::array<std::array<long double,3>,3> gram{};
    for(unsigned panel=0;panel<16;++panel)for(unsigned g=0;g<8;++g){
        const auto t=(panel+(gaussX[g]+1)/2)/16,weight=gaussW[g]/32*norm(arc.derivative(t));const auto b=basis(t);
        for(unsigned i=0;i<3;++i)for(unsigned j=0;j<3;++j)gram[i][j]+=weight*b[i]*b[j];}
    std::array<std::array<long double,3>,3> l{};
    for(unsigned i=0;i<3;++i)for(unsigned j=0;j<=i;++j){auto v=gram[i][j];for(unsigned k=0;k<j;++k)v-=l[i][k]*l[j][k];
        require(std::isfinite(v)&&(i!=j||v>0),Failure::InvalidMoments,"nonpositive surface moment Gram matrix");l[i][j]=i==j?std::sqrt(v):v/l[j][j];}
    std::array<long double,3> y{},result{};
    for(unsigned i=0;i<3;++i){long double v=moments[i];for(unsigned j=0;j<i;++j)v-=l[i][j]*y[j];y[i]=v/l[i][i];}
    for(int i=2;i>=0;--i){auto v=y[i];for(unsigned j=i+1;j<3;++j)v-=l[j][i]*result[j];result[i]=v/l[i][i];}
    for(auto x:result)require(std::isfinite(x),Failure::InvalidMoments,"nonfinite reconstructed material charge density");return result;
}
struct MaterialTransfer
{
    const Smooth& source;const Smooth& target;const Grid& grid;const Topology& targetTopology;
    const ElectricSurfaceChargeTransferConfig& config;ElectricSurfaceChargeTransferReport& report;
    long double dt=0;bool firstY=false;
    std::vector<std::array<Sum,3>> targetMoments;
    std::array<Sum,3> measuredSource;
    Sum absoluteCharge;
    std::array<Sum,2> sourcePosition,advectedPosition,projectedPosition;
    std::map<std::pair<std::size_t,long double>,Projection> endpointProjectionMemo;
    std::vector<std::optional<NativeArcConstants>> nativeArcCache;
    std::optional<NativeProductWeights> nativeWeightCache;
    std::vector<std::optional<ExactDistanceBases>> exactDistanceBases;
    Point map(Point p) const
    {if(dt==0)return p;p=directional(p,firstY?1:0,dt,grid);return directional(p,firstY?0:1,dt,grid);}
    Projection projectOnArc(std::size_t arc,Point point,bool exactOnly=false)
    {
        const auto start=WallClock::now();
        auto result=exactOnly?projectArcExact(target.arcs()[arc],arc,point,exactDistanceBases[arc],CoreProjectionReason::ForcedDistance)
            :projectArc(target.arcs()[arc],arc,point,nativeArcCache,nativeWeightCache,exactDistanceBases[arc]);
        const auto seconds=std::chrono::duration<double>(WallClock::now()-start).count();
        if(result.usedCoreRoots)report.coreArcProjectionSeconds+=seconds;
        else if(result.exactCertificate)report.identityArcProjectionSeconds+=seconds;
        else report.nativeArcProjectionSeconds+=seconds;
        if(result.exactCertificate)++report.exactArcProjections;else ++report.nativeArcProjections;
        if(result.usedCoreRoots){++report.rootIsolationCalls;
            switch(result.coreReason){
                case CoreProjectionReason::Hessian:++report.coreHessianFallbacks;break;
                case CoreProjectionReason::EndpointSign:++report.coreEndpointSignFallbacks;break;
                case CoreProjectionReason::Nonfinite:++report.coreNonfiniteFallbacks;break;
                case CoreProjectionReason::RootWidth:++report.coreRootWidthFallbacks;break;
                case CoreProjectionReason::ForcedDistance:++report.coreForcedDistanceFallbacks;break;
                case CoreProjectionReason::None:break;}
            if(result.endpoint<0)++report.coreInteriorProjections;else ++report.coreEndpointProjections;
            if(result.preparedGeometryBase)++report.geometryBasePreparations;}
        return result;
    }
    Projection project(Point original,Point normal)
    {
        const PassiveTimer timer{report.pointProjectionSeconds};
        const auto p=map(original);const auto h=std::min(grid.hx,grid.hy),displacement=norm(sub(p,original));
        report.maximumDisplacement=std::max(report.maximumDisplacement,displacement);
        require(displacement<=config.maximumDisplacementInCells*h,Failure::ExcessiveDisplacement,"surface material displacement exceeds short-step qualification");
        struct LocalArc{std::size_t arc;long double lower;std::optional<Projection> projection;};
        const auto cell=grid.index(p);std::vector<LocalArc> local;
        for(int y=-1;y<=1;++y)for(int x=-1;x<=1;++x){auto found=grid.arcs.find({cell.first+x,cell.second+y});
            if(found!=grid.arcs.end())local.push_back({found->second,aabbDistanceLower(target.arcs()[found->second],p),{}});}
        require(!local.empty(),Failure::ProjectionDistance,"material has no local target arc");
        std::sort(local.begin(),local.end(),[](const auto& a,const auto& b){return a.lower<b.lower||(a.lower==b.lower&&a.arc<b.arc);});
        auto bestUpper=std::numeric_limits<long double>::infinity();
        for(auto& entry:local){if(entry.lower>bestUpper)continue;
            entry.projection=projectOnArc(entry.arc,p);bestUpper=std::min(bestUpper,entry.projection->squaredDistance.upper);}
        // An overlapping native interval is refined against the SAME exact
        // distance polynomial before comparing it with any other candidate.
        std::size_t bestIndex=0;bool haveBest=false;
        for(std::size_t i=0;i<local.size();++i)if(local[i].projection&&(!haveBest||local[i].projection->squaredDistance.upper<local[bestIndex].projection->squaredDistance.upper))
            {bestIndex=i;haveBest=true;}
        bool unresolved=false;
        for(std::size_t i=0;i<local.size();++i)if(i!=bestIndex&&local[i].projection&&
            local[i].projection->squaredDistance.lower<=local[bestIndex].projection->squaredDistance.upper&&
            !sameNode(*local[i].projection,*local[bestIndex].projection,target))unresolved=true;
        if(unresolved){
            for(auto& entry:local)if(entry.projection&&entry.projection->squaredDistance.lower<=bestUpper&&!entry.projection->exactCertificate)
                entry.projection=projectOnArc(entry.arc,p,true);
            for(std::size_t i=0;i<local.size();++i)if(local[i].projection&&local[i].projection->squaredDistance.upper<local[bestIndex].projection->squaredDistance.upper)bestIndex=i;}
        // Choose by strict certificates against ALL computed candidates,
        // including a certified interior minimum that beats the rounded
        // endpoint winner at their shared actual target node.
        std::vector<Detail::SurfaceChargeProjection::LocalMinimum> minima;minima.reserve(local.size());
        for(std::size_t i=0;i<local.size();++i)if(local[i].projection){const auto& hit=*local[i].projection;
            const auto& arc=target.arcs()[hit.arc];std::array<Point,2> endpoints{};
            for(unsigned side=0;side<2;++side){endpoints[side]=target.nodes()[arc.nodes[side]].point;
                const auto control=arc.controlPoints[side==0?0:4];
                require(control.x==endpoints[side].x&&control.y==endpoints[side].y,Failure::ForeignSource,
                    "target arc endpoint differs from its actual owned node");}
            minima.push_back({hit.arc,arc.nodes,endpoints,{hit.squaredDistance.lower,hit.squaredDistance.upper},hit.endpoint,hit.minimum});}
        const auto selected=Detail::SurfaceChargeProjection::selectLocalMinimum(minima);
        if(selected)for(std::size_t i=0;i<local.size();++i)if(local[i].arc==minima[*selected].arc){bestIndex=i;break;}
        auto best=*local[bestIndex].projection;
        if(!selected)for(const auto& entry:local)if(entry.projection&&entry.arc!=best.arc){const auto& other=*entry.projection;
            if(!sameNode(best,other,target)&&!(best.squaredDistance.upper<other.squaredDistance.lower)){
                std::ostringstream context;context<<std::setprecision(std::numeric_limits<long double>::max_digits10)
                    <<"local minimum projection distances are not separated source_point=("<<original.x<<','<<original.y
                    <<") mapped_point=("<<p.x<<','<<p.y<<") source_normal=("<<normal.x<<','<<normal.y<<')'
                    <<" dt="<<dt<<" first_y="<<firstY<<" spacing=("<<grid.hx<<','<<grid.hy<<')'
                    <<" source_identity="<<source.sourceIdentity().value()<<" target_identity="<<target.sourceIdentity().value()
                    <<" source_clock=("<<source.sourceTime().step<<','<<source.sourceTime().begin<<','<<source.sourceTime().end<<','
                    <<static_cast<int>(source.sourceTime().position)<<") target_clock=("<<target.sourceTime().step<<','
                    <<target.sourceTime().begin<<','<<target.sourceTime().end<<','<<static_cast<int>(target.sourceTime().position)<<')'
                    <<" best_arc="<<best.arc<<" conflicting_arc="<<other.arc;
                for(const auto& candidate:local){const auto& arc=target.arcs()[candidate.arc];
                    context<<" local=(arc="<<candidate.arc<<",native_cell="<<arc.nativeCell<<",nodes=("<<arc.nodes[0]<<','<<arc.nodes[1]
                        <<"),aabb_squared_lower="<<candidate.lower<<",evaluated="<<candidate.projection.has_value()<<",controls=";
                    for(const auto control:arc.controlPoints)context<<'('<<control.x<<','<<control.y<<')';
                    if(candidate.projection){const auto& projection=*candidate.projection;
                        context<<",t="<<projection.t<<",parameter=["<<projection.parameter.lower<<','<<projection.parameter.upper
                            <<"],endpoint="<<projection.endpoint<<",squared_distance=["<<projection.squaredDistance.lower<<','
                            <<projection.squaredDistance.upper<<"],exact_certificate="<<projection.exactCertificate
                            <<",used_core_roots="<<projection.usedCoreRoots<<",minimum_certificate="<<static_cast<int>(projection.minimum);}
                    context<<')';}
                throw ElectricSurfaceChargeTransferError(Failure::AmbiguousProjection,context.str().c_str());}}
        require(selected.has_value(),Failure::AmbiguousProjection,"local minima lack a complete strict comparison certificate");
        const auto limit=config.maximumProjectionDistanceInCells*h;const auto limitSquared=square(exact(limit));
        if(best.squaredDistance.upper>limitSquared.lower&&!best.exactCertificate)best=projectOnArc(best.arc,p,true);
        require(best.squaredDistance.upper<=limitSquared.lower,Failure::ProjectionDistance,"material projection exceeds declared distance");
        // Closest-point pruning never suppresses the separate branch-radius
        // check. Nearby nonadjacent arcs need an actual, length-bounded local
        // curve chain; graph connectivity by itself is insufficient.
        std::set<std::size_t> localArcs;for(const auto& entry:local)localArcs.insert(entry.arc);
        for(auto& entry:local)if(entry.arc!=best.arc&&!targetTopology.adjacent[best.arc].count(entry.arc)&&entry.lower<=limitSquared.upper){
            if(!entry.projection)entry.projection=projectOnArc(entry.arc,p);
            bool chain=entry.projection->squaredDistance.lower<=limitSquared.upper&&
                sameLocalContinuousChain(best,*entry.projection,target,targetTopology,localArcs,limit);
            // The native parameter interval already certifies this closest
            // point. A bounded chain certifies the branch condition directly;
            // exact root refinement cannot strengthen a condition already met.
            if(!chain&&entry.projection->squaredDistance.lower<=limitSquared.upper&&!entry.projection->exactCertificate){
                entry.projection=projectOnArc(entry.arc,p,true);
                chain=entry.projection->squaredDistance.lower<=limitSquared.upper&&
                    sameLocalContinuousChain(best,*entry.projection,target,targetTopology,localArcs,limit);}
            if(chain)++report.localChainCertificates;
            require(entry.projection->squaredDistance.lower>limitSquared.upper||chain,
                Failure::AmbiguousProjection,"material lies near a distinct nonadjacent surface branch");}
        require(dot(normal,target.arcs()[best.arc].normal(best.t))>.25L,Failure::AmbiguousProjection,"surface material orientation changed excessively");
        report.maximumProjectionDistance=std::max(report.maximumProjectionDistance,best.distance);++report.materialSamples;return best;
    }
    Projection projectEndpoint(std::size_t sourceArc,long double t)
    {
        ++report.endpointProjectionRequests;
        const auto key=std::make_pair(sourceArc,t);const auto found=endpointProjectionMemo.find(key);
        if(found!=endpointProjectionMemo.end()){
            // Same immutable target/slab/config and identical source point
            // AND normal. Its complete qualification is already retained;
            // every logical request still counts and all maxima are idempotent.
            ++report.endpointProjectionMemoHits;++report.materialSamples;return found->second;}
        const auto& arc=source.arcs()[sourceArc];const auto result=project(arc.sample(t),arc.normal(t));
        endpointProjectionMemo.emplace(key,result);return result;
    }
    void integrate(std::size_t sourceArc,const std::array<long double,3>& density,long double low,long double high,unsigned depth)
    {
        const auto& arc=source.arcs()[sourceArc];
        std::array<Projection,8> targets;std::array<long double,8> charges{},parameters{};
        const auto center=(low+high)/2,radius=(high-low)/2;bool single=true;
        for(unsigned g=0;g<8;++g){const auto t=center+radius*gaussX[g];const auto b=basis(t);
            parameters[g]=t;charges[g]=(density[0]*b[0]+density[1]*b[1]+density[2]*b[2])*norm(arc.derivative(t))*radius*gaussW[g];
            targets[g]=project(arc.sample(t),arc.normal(t));single=single&&(g==0||targets[g].arc==targets[0].arc);}
        // Interior Gauss sites can miss a short terminal crossing. Both actual
        // parameter endpoints must belong to the same target arc as well. A
        // shared target node is equivalent for this measure-zero endpoint.
        if(single){const auto owner=targets[0].arc;
            const auto endpointMatches=[&](long double t){const auto hit=projectEndpoint(sourceArc,t);
                if(hit.arc==owner)return true;
                if(hit.endpoint<0)return false;
                const auto node=target.arcs()[hit.arc].nodes[hit.endpoint];
                return target.arcs()[owner].nodes[0]==node||target.arcs()[owner].nodes[1]==node;};
            single=endpointMatches(low)&&endpointMatches(high);}
        if(!single&&depth<32){++report.parameterSplits;integrate(sourceArc,density,low,center,depth+1);integrate(sourceArc,density,center,high,depth+1);return;}
        if(!single)report.maximumTerminalParameterWidth=std::max(report.maximumTerminalParameterWidth,high-low);
        ++report.materialPanels;
        for(unsigned g=0;g<8;++g){const auto targetBasis=basis(targets[g].t),sourceBasis=basis(parameters[g]);
            for(unsigned k=0;k<3;++k){targetMoments[targets[g].arc][k].add(charges[g]*targetBasis[k]);measuredSource[k].add(charges[g]*sourceBasis[k]);}
            const auto sourcePoint=arc.sample(parameters[g]),advectedPoint=map(sourcePoint),projectedPoint=target.arcs()[targets[g].arc].sample(targets[g].t);
            sourcePosition[0].add(charges[g]*sourcePoint.x);sourcePosition[1].add(charges[g]*sourcePoint.y);
            advectedPosition[0].add(charges[g]*advectedPoint.x);advectedPosition[1].add(charges[g]*advectedPoint.y);
            projectedPosition[0].add(charges[g]*projectedPoint.x);projectedPosition[1].add(charges[g]*projectedPoint.y);
            absoluteCharge.add(std::abs(charges[g]));}
    }
};
}
ElectricSurfaceChargeTransferError::ElectricSurfaceChargeTransferError(Failure reason,const char* message)
    :std::runtime_error(std::string("[Electric/SurfaceChargeTransfer] ")+message),reason_(reason){}
struct ElectricSurfaceChargeTransport::Data
{
    ElectricVofSurfaceCharge source,candidate;
    std::optional<Interface::GeometricVofTransport> transport;
    ElectricSurfaceChargeTransferReport report;
    Data(ElectricVofSurfaceCharge a,ElectricVofSurfaceCharge b):source(std::move(a)),candidate(std::move(b)){}
};
ElectricSurfaceChargeTransport::ElectricSurfaceChargeTransport(std::shared_ptr<const Data> data):data_(std::move(data)){}
ElectricSurfaceChargeTransport ElectricSurfaceChargeTransport::prepareStage(const QuadtreeMesh& mesh,
    const ElectricVofSurfaceCharge& source,const Interface::GeometricVofTransport& transport,
    const Smooth& target,const ElectricSurfaceChargeTransferConfig& config)
{return prepare(mesh,source,&transport,target,config);}
ElectricSurfaceChargeTransport ElectricSurfaceChargeTransport::prepareSameEnd(const QuadtreeMesh& mesh,
    const ElectricVofSurfaceCharge& source,const Smooth& target,const ElectricSurfaceChargeTransferConfig& config)
{return prepare(mesh,source,nullptr,target,config);}
ElectricSurfaceChargeTransport ElectricSurfaceChargeTransport::prepare(const QuadtreeMesh& mesh,
    const ElectricVofSurfaceCharge& charge,const Interface::GeometricVofTransport* transport,
    const Smooth& target,const ElectricSurfaceChargeTransferConfig& config)
{
    const auto prepareStart=WallClock::now();
    const auto& source=charge.geometry();source.validateFor(mesh,source.source());target.validateFor(mesh,target.source());
    charge.validateFor(mesh,source,source.sourceTime());
    require(std::isfinite(config.maximumDisplacementInCells)&&config.maximumDisplacementInCells>0&&config.maximumDisplacementInCells<=.5L&&
        std::isfinite(config.maximumProjectionDistanceInCells)&&config.maximumProjectionDistanceInCells>0&&config.maximumProjectionDistanceInCells<=.25L&&
        std::isfinite(config.relativeChargeTolerance)&&config.relativeChargeTolerance>0&&config.relativeChargeTolerance<=1e-10L,
        Failure::InvalidMoments,"invalid transfer qualification or loosened conservation gate");
    const auto old=charge.timeLayer(),next=target.sourceTime();bool half=false;long double dt=0;bool firstY=false;
    require(old.position==Clock::Position::Initial||old.position==Clock::Position::End,Failure::InvalidClock,"transfer requires an accepted source clock");
    if(transport){
        require(transport->source().identity()==source.sourceIdentity(),Failure::ForeignSource,"transfer source is not the accepted charge geometry owner");
        transport->validateFor(mesh,source.source(),transport->advector(),transport->liquidDensity(),transport->gasDensity());
        half=next.position==Clock::Position::Midpoint;
        const auto& stage=half?transport->midpoint():transport->endpoint();
        require(stage.identity()==target.sourceIdentity()&&sameClock(next,stage.timeLayer()),Failure::ForeignSource,"target is not the actual half/full VOF stage");
        if(const auto* geometry=transport->smoothSource()){
            const auto* retained=half?transport->smoothMidpoint():transport->smoothEndpoint();
            require(geometry->sameInstance(source)&&retained&&retained->sameInstance(target),
                Failure::ForeignSource,"surface charge and phase flux must consume the same actual smooth source/target owners");
        }
        dt=static_cast<long double>(next.time())-old.time();firstY=(next.step%2)==0;
    }else {require(sameClock(old,next),Failure::InvalidClock,"same-End remap requires identical complete clocks");
        require(source.sourceIdentity()==target.sourceIdentity(),Failure::ForeignSource,"same-End reconstruction requires the same actual VOF publication");}
    const auto sourceTopology=topology(source),targetTopology=topology(target);
    require(sourceTopology.endpoints.size()==targetTopology.endpoints.size(),Failure::UnsupportedTopology,"surface topology changed");
    const auto grid=gridFor(mesh,source,target,transport,half);const auto h=std::min(grid.hx,grid.hy);
    ElectricSurfaceChargeTransferReport report;report.sourceArcs=source.arcs().size();report.targetArcs=target.arcs().size();
    if(!transport||zeroNativeStage(*transport))if(const auto mapping=identityArcs(source,target)){
        std::vector<std::array<double,3>> moments(target.arcs().size());Sum total;
        for(std::size_t a=0;a<mapping->size();++a){moments[(*mapping)[a]]=charge.moments()[a];total.add(charge.moments()[a][0]);}
        require(std::isfinite(total.value),Failure::InvalidMoments,"unrepresentable identity total charge");
        report.mode=ElectricSurfaceChargeTransferMode::ExactPartitionIdentity;
        report.integralDiagnosticsAvailable=false;report.identityArcTransfers=mapping->size();
        report.sourceCharge=total.value;report.targetCharge=total.value;
        const auto unavailable=std::numeric_limits<long double>::quiet_NaN();
        report.sourceAbsoluteCharge=unavailable;report.sourceMomentQuadratureResidual=unavailable;
        report.sourcePositionMoment.fill(unavailable);report.advectedPositionMoment.fill(unavailable);
        report.projectedPositionMoment.fill(unavailable);
        auto candidate=ElectricVofSurfaceCharge::fromPrepared(mesh,target,next,std::move(moments));
        auto data=std::make_shared<Data>(charge,std::move(candidate));data->report=report;if(transport)data->transport=*transport;
        data->report.prepareSeconds=std::chrono::duration<double>(WallClock::now()-prepareStart).count();
        return ElectricSurfaceChargeTransport(std::move(data));}
    MaterialTransfer transfer{source,target,grid,targetTopology,config,report,dt,firstY,std::vector<std::array<Sum,3>>(target.arcs().size()),{},{},{},{},{}};
    transfer.nativeArcCache.resize(target.arcs().size());
    transfer.exactDistanceBases.resize(target.arcs().size());
    for(const auto& p:sourceTopology.endpoints){const auto mapped=transfer.map(p);
        require(norm(sub(mapped,p))<=64*std::numeric_limits<double>::epsilon()*h,Failure::UnsupportedBoundary,"moving boundary surface endpoints are unsupported");
        auto distance=std::numeric_limits<long double>::infinity();for(const auto& q:targetTopology.endpoints)distance=std::min(distance,norm(sub(p,q)));
        require(distance<=64*std::numeric_limits<double>::epsilon()*h,Failure::UnsupportedBoundary,"boundary surface endpoint changed");}
    Sum sourceTotal,targetTotal;long double absoluteMomentScale=0;
    const auto integrationStart=WallClock::now();
    for(std::size_t a=0;a<source.arcs().size();++a){const auto& moment=charge.moments()[a];const auto density=densityCoefficients(source.arcs()[a],moment);
        ++report.sourceDensityReconstructions;
        transfer.endpointProjectionMemo.clear();
        transfer.measuredSource={};
        for(unsigned panel=0;panel<16;++panel)transfer.integrate(a,density,panel/16.L,(panel+1)/16.L,0);
        for(unsigned k=0;k<3;++k){absoluteMomentScale+=std::abs(moment[k]);report.sourceMomentQuadratureResidual=std::max(report.sourceMomentQuadratureResidual,std::abs(transfer.measuredSource[k].value-moment[k]));}
        sourceTotal.add(moment[0]);}
    report.materialIntegrationSeconds=std::chrono::duration<double>(WallClock::now()-integrationStart).count();
    std::vector<std::array<double,3>> moments(target.arcs().size());
    for(std::size_t a=0;a<moments.size();++a)for(unsigned k=0;k<3;++k){moments[a][k]=static_cast<double>(transfer.targetMoments[a][k].value);
        require(std::isfinite(moments[a][k]),Failure::InvalidMoments,"unrepresentable transferred surface moment");if(k==0)targetTotal.add(moments[a][0]);}
    report.sourceCharge=sourceTotal.value;report.targetCharge=targetTotal.value;report.sourceAbsoluteCharge=transfer.absoluteCharge.value;
    for(unsigned k=0;k<2;++k){report.sourcePositionMoment[k]=transfer.sourcePosition[k].value;report.advectedPositionMoment[k]=transfer.advectedPosition[k].value;report.projectedPositionMoment[k]=transfer.projectedPosition[k].value;}
    report.chargeResidual=report.targetCharge-report.sourceCharge;
    const auto scale=std::max(report.sourceAbsoluteCharge,absoluteMomentScale);
    report.relativeChargeResidual=scale>0?std::abs(report.chargeResidual)/scale:std::abs(report.chargeResidual);
    require(report.relativeChargeResidual<=config.relativeChargeTolerance&&report.sourceMomentQuadratureResidual<=config.relativeChargeTolerance*scale,
        Failure::Conservation,"material quadrature/charge conservation failed; no charge correction is permitted");
    auto candidate=ElectricVofSurfaceCharge::fromPrepared(mesh,target,next,std::move(moments));
    auto data=std::make_shared<Data>(charge,std::move(candidate));data->report=report;if(transport)data->transport=*transport;
    data->report.prepareSeconds=std::chrono::duration<double>(WallClock::now()-prepareStart).count();
    return ElectricSurfaceChargeTransport(std::move(data));
}
const ElectricVofSurfaceCharge& ElectricSurfaceChargeTransport::sourceCharge() const{return data_->source;}
const ElectricVofSurfaceCharge& ElectricSurfaceChargeTransport::chargeCandidate() const{return data_->candidate;}
const Interface::GeometricVofTransport* ElectricSurfaceChargeTransport::interfaceTransport() const{return data_->transport?&*data_->transport:nullptr;}
const ElectricSurfaceChargeTransferReport& ElectricSurfaceChargeTransport::report() const{return data_->report;}
bool ElectricSurfaceChargeTransport::sameInstance(const ElectricSurfaceChargeTransport& other) const noexcept{return data_&&data_==other.data_;}
}







