#include "Modules/Interface/VofSmoothGeometry.h"
#include "Modules/Interface/Detail/VofSmoothPhaseMeasure.h"
#include "Core/Geometry/UnitPolynomialRoots.h"
#include "Core/Geometry/NativeOutwardRounding.h"
#include "Core/Geometry/Detail/FiniteVolumeClosure.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <iomanip>
#include <sstream>
#include <optional>
#include <string>
#include <utility>

namespace LAMP::Interface
{
namespace
{
namespace Cut=Geometry::PlanarCellCut;
using Point=Cut::Point;
using Failure=VofSmoothGeometryFailure;
Point add(Point a,Point b) { return {a.x+b.x,a.y+b.y}; }
Point sub(Point a,Point b) { return {a.x-b.x,a.y-b.y}; }
Point mul(Point a,long double b) { return {a.x*b,a.y*b}; }
long double dot(Point a,Point b) { return a.x*b.x+a.y*b.y; }
long double norm(Point a) { return std::hypot(a.x,a.y); }
long double coordinate(Point a,unsigned axis) { return axis?a.y:a.x; }
Point unit(Point a,int cell)
{
    const auto length=norm(a);
    if(!(length>0)||!std::isfinite(length))
        throw VofSmoothGeometryError(Failure::InvalidTangent,cell,"zero or nonfinite shared tangent");
    return mul(a,1/length);
}
void require(bool condition,Failure failure,int cell,const char* message)
{ if(!condition)throw VofSmoothGeometryError(failure,cell,message); }
template<std::size_t N> Point bezier(std::array<Point,N> value,long double t)
{
    for(std::size_t n=N-1;n>0;--n)for(std::size_t i=0;i<n;++i)
        value[i]=add(mul(value[i],1-t),mul(value[i+1],t));
    return value[0];
}
void powerProjection(VofSmoothArc& arc)
{
    constexpr long double choose[5][5]={{1,0,0,0,0},{1,1,0,0,0},{1,2,1,0,0},
        {1,3,3,1,0},{1,4,6,4,1}};
    for(unsigned k=0;k<5;++k){
        Point value{};
        for(unsigned j=0;j<=k;++j)
            value=add(value,mul(arc.controlPoints[j],choose[k][j]*((k-j)%2?-1:1)));
        arc.power[k]=mul(value,choose[4][k]);
    }
}
struct Rule { std::array<long double,16> x{},w{}; };
const Rule& gauss()
{
    // Fixed positive GL16 on [0,1]. No process-initialization root solve.
    static constexpr Rule rule{{{
        .005299532504175033701922913274834L,.027712488463383711961005792232696L,
        .067184398806084128059766051143803L,.122297795822498483052449402576279L,
        .191061877798678125776664117975604L,.270991611171386306828790278508211L,
        .359198224610370543384769749269251L,.452493745081181279907340332287521L,
        .547506254918818720092659667712479L,.640801775389629456615230250730749L,
        .729008388828613693171209721491789L,.808938122201321874223335882024396L,
        .877702204177501516947550597423721L,.932815601193915871940233948856197L,
        .972287511536616288038994207767304L,.994700467495824966298077086725166L}},{{
        .013576229705877047425890286228009L,.031126761969323946431421918497189L,
        .047579255841246392404962553801123L,.062314485627766936026238141096008L,
        .074797994408288366040750865273739L,.084578259697501269594656039515180L,
        .091301707522461794433381833984610L,.094725305227534248142698361604142L,
        .094725305227534248142698361604142L,.091301707522461794433381833984610L,
        .084578259697501269594656039515180L,.074797994408288366040750865273739L,
        .062314485627766936026238141096008L,.047579255841246392404962553801123L,
        .031126761969323946431421918497189L,.013576229705877047425890286228009L}}};
    return rule;
}
constexpr std::array<long double,4> eta{{.069431844202973712388L,.33000947820757186760L,
    .66999052179242813240L,.93056815579702628761L}};
constexpr std::array<long double,4> etaWeight{{.17392742256872692869L,.32607257743127307131L,
    .32607257743127307131L,.17392742256872692869L}};
Point midpoint(const Cut::Rectangle& box)
{ return {(static_cast<long double>(box.lower.x)+box.upper.x)/2,(static_cast<long double>(box.lower.y)+box.upper.y)/2}; }
struct Endpoint {int cell=-1,face=-1;Point point,normal;};
unsigned graphOuterPhase(const VofSmoothArc& arc,const Cut::Rectangle& box,long double u,bool liquidLower)
{
    const auto first=arc.controlPoints.front(),last=arc.controlPoints.back();
    const auto a=coordinate(first,arc.graphAxis),b=coordinate(last,arc.graphAxis);
    const auto endpoint=u<std::min(a,b)?(a<b?first:last):(a>b?first:last);
    const auto lower=arc.graphAxis?box.lower.x:box.lower.y;
    const auto upper=arc.graphAxis?box.upper.x:box.upper.y;
    const auto v=coordinate(endpoint,1-arc.graphAxis);
    require(v==lower||v==upper,Failure::InvalidCurve,arc.nativeCell,"outer graph strip has no transverse boundary endpoint");
    const bool below=v==upper;
    return below==liquidLower?0:1;
}
// A graph-strip parametrization of the actual curved phase. Its Jacobian is
// positive after curve certification. During volume-root search its signed
// polynomial extension is used, never an indicator or a clipped phase.
std::array<long double,2> phaseMeasure(const VofSmoothArc& arc,
    const GeometricVofPlicCell& plic,bool rz,
    std::array<std::vector<VofSmoothWeightedPoint>,2>* output=nullptr)
{
    const auto& box=plic.rectangle;
    const unsigned axis=arc.graphAxis;
    const long double lo=axis?box.lower.y:box.lower.x,hi=axis?box.upper.y:box.upper.x;
    const long double lower=axis?box.lower.x:box.lower.y,upper=axis?box.upper.x:box.upper.y;
    const auto a=coordinate(arc.controlPoints.front(),axis),b=coordinate(arc.controlPoints.back(),axis);
    const long double direction=b>a?1:-1;
    const bool liquidLower=axis?direction>0:direction<0;
    // Measure the SAME actual world-control polynomial relative to each
    // transverse face before evaluating it. World sample minus face can lose
    // a small phase width although these control-to-face differences are
    // representable. Keep both widths independent, including the signed
    // polynomial extension used while fitting; never clamp a phase.
    std::array<Point,5> phaseWidths;
    for(unsigned i=0;i<5;++i){const auto v=coordinate(arc.controlPoints[i],1-axis);
        phaseWidths[i]={v-lower,upper-v};}
    std::array<long double,2> measured{};
    const auto accumulate=[&](unsigned phase,long double u,long double v0,long double v1,long double width,long double jacobian){
        if(width==0||jacobian==0)return;
        for(unsigned j=0;j<4;++j){
            const auto v=(1-eta[j])*v0+eta[j]*v1;
            const Point point=axis?Point{v,u}:Point{u,v};
            const auto weight=jacobian*width*etaWeight[j]*(rz?2*std::acos(-1.L)*point.x:1);
            measured[phase]+=weight;
            if(output) {
                require(weight>0&&std::isfinite(weight),Failure::InvalidCurve,arc.nativeCell,"nonpositive phase quadrature weight");
                (*output)[phase].push_back({point,weight});
            }
        }
    };
    const auto& rule=gauss();
    for(unsigned q=0;q<16;++q){
        const auto p=arc.sample(rule.x[q]),d=arc.derivative(rule.x[q]);
        const auto widths=bezier(phaseWidths,rule.x[q]);
        const auto u=coordinate(p,axis),v=coordinate(p,1-axis);
        const auto jacobian=rule.w[q]*direction*coordinate(d,axis);
        accumulate(liquidLower?0:1,u,lower,v,widths.x,jacobian);
        accumulate(liquidLower?1:0,u,v,upper,widths.y,jacobian);
    }
    // Outside the graph's endpoint interval the whole transverse strip has
    // one phase. Its side is fixed by the actual endpoint/native-face topology.
    const auto outer=[&](long double begin,long double end){
        if(!(end>begin))return;
        const unsigned phase=graphOuterPhase(arc,box,(begin+end)/2,liquidLower);
        for(unsigned q=0;q<16;++q)
            accumulate(phase,(1-rule.x[q])*begin+rule.x[q]*end,lower,upper,upper-lower,(end-begin)*rule.w[q]);
    };
    outer(lo,std::min(a,b));outer(std::max(a,b),hi);
    return measured;
}
using StripPolynomial=std::vector<Geometry::UnitPolynomialTerm>;
StripPolynomial stripCoordinate(const VofSmoothArc& arc,unsigned axis,bool reverse,
    long double anchor=0,long double sign=1)
{
    // Retain the Bezier basis and anchor as separate exact Core terms. In
    // particular, face +/- depth is never rounded into a world coordinate.
    constexpr int basis[5][5]={{1,-4,6,-4,1},{0,4,-12,12,-4},{0,0,6,-12,6},
        {0,0,0,4,-4},{0,0,0,0,1}};
    StripPolynomial result;
    for(unsigned j=0;j<5;++j){Geometry::UnitPolynomialTerm term;
        term.factor=sign*coordinate(arc.controlPoints[reverse?4-j:j],axis);
        for(unsigned k=0;k<5;++k)term.coefficients[k]=basis[j][k];result.push_back(term);}
    if(anchor!=0){Geometry::UnitPolynomialTerm term;term.factor=-sign*anchor;term.coefficients[0]=1;result.push_back(term);}
    return result;
}
StripPolynomial stripConstant(long double value)
{Geometry::UnitPolynomialTerm term;term.factor=value;term.coefficients[0]=1;return {term};}
StripPolynomial stripDifference(StripPolynomial first,const StripPolynomial& second)
{
    for(auto term:second){for(auto& coefficient:term.coefficients)coefficient=-coefficient;first.push_back(term);}return first;
}
StripPolynomial stripDerivative(StripPolynomial polynomial)
{
    for(auto& term:polynomial){for(unsigned k=0;k+1<term.coefficients.size();++k)term.coefficients[k]=(k+1)*term.coefficients[k+1];
        term.coefficients.back()=0;}
    return polynomial;
}
std::vector<double> stripBinaryParts(long double value)
{
    std::vector<double> result;
    for(unsigned k=0;value!=0&&k<4;++k){const auto part=static_cast<double>(value);
        require(std::isfinite(part)&&part!=0,Failure::VolumeResidual,-1,"unrepresentable smooth strip coefficient");
        result.push_back(part);value-=part;}
    require(value==0,Failure::VolumeResidual,-1,"unresolved smooth strip coefficient expansion");return result;
}
StripPolynomial stripScale(const StripPolynomial& polynomial,long double scale)
{
    using Exact=Geometry::Detail::FiniteVolume::Dyadic;StripPolynomial result;
    for(const auto& term:polynomial){auto value=Exact::product({term.factor,term.weight,scale});
        for(unsigned part=0;value.compare(Exact{})!=0&&part<4;++part){auto copy=term;copy.factor=value.rounded();copy.weight=1;
            require(std::isfinite(copy.factor)&&copy.factor!=0,Failure::VolumeResidual,-1,"unrepresentable scaled strip coefficient");
            result.push_back(copy);value.add(copy.factor,-1);}
        require(value.compare(Exact{})==0,Failure::VolumeResidual,-1,"unresolved exact scaled strip coefficient");}
    return result;
}
StripPolynomial stripEndpoint(const StripPolynomial& polynomial,bool last)
{
    StripPolynomial result;
    for(auto term:polynomial){double value=term.coefficients[0];
        if(last)for(unsigned k=1;k<term.coefficients.size();++k)value+=term.coefficients[k];
        term.coefficients={};term.coefficients[0]=value;result.push_back(term);}
    return result;
}
StripPolynomial stripProduct(const StripPolynomial& first,const StripPolynomial& second)
{
    StripPolynomial result;
    for(const auto& a:first)for(const auto& b:second){
        // Inputs have integer Bezier/derivative coefficients. Their degree-7
        // convolution is exact; the two coordinate factors remain separate.
        std::array<double,12> coefficients{};
        for(unsigned i=0;i<12;++i)for(unsigned j=0;i+j<12;++j)coefficients[i+j]+=a.coefficients[i]*b.coefficients[j];
        if(std::all_of(coefficients.begin(),coefficients.end(),[](double value){return value==0;}))continue;
        for(const auto x:stripBinaryParts(a.factor))for(const auto y:stripBinaryParts(b.factor)){
            Geometry::UnitPolynomialTerm term;term.coefficients=coefficients;term.factor=x;term.weight=y;
            require(a.weight==1&&b.weight==1,Failure::InvalidCurve,-1,"nested strip polynomial product");result.push_back(term);}}
    return result;
}
struct StripAccumulator
{
    using Exact=Geometry::Detail::FiniteVolume::Dyadic;
    std::array<Exact,2> phase,error;
    void positive(unsigned p,const StripPolynomial& density,double lower,double upper,
        long double scale,long double tolerance,long double divisor=1)
    {
        if(lower==upper)return;
        require(std::isfinite(divisor)&&divisor>0,Failure::VolumeResidual,-1,"invalid strip integral divisor");
        auto budget=tolerance*scale*divisor/64;
        auto integral=Geometry::integratePositiveUnitPolynomial(density,lower,upper,budget);
        for(unsigned iteration=0;iteration<3&&integral.value>0&&integral.absoluteError>tolerance*integral.value/64;++iteration){
            budget=tolerance*integral.value/64;
            integral=Geometry::integratePositiveUnitPolynomial(density,lower,upper,budget);}
        require(integral.value>=0&&integral.absoluteError>=0&&
            (integral.value==0?integral.absoluteError==0:integral.absoluteError<=tolerance*integral.value/32),
            Failure::VolumeResidual,-1,"positive smooth strip integral lacks its own-phase error budget");
        const auto value=integral.value/divisor;
        require(std::isfinite(value)&&(integral.value==0||value>0),Failure::VolumeResidual,-1,"unrepresentable positive strip quotient");
        auto uncertainty=integral.absoluteError/divisor;
        if(divisor!=1)uncertainty+=std::numeric_limits<long double>::epsilon()*value;
        if(uncertainty>0)uncertainty=std::nextafter(uncertainty,std::numeric_limits<long double>::infinity());
        phase[p].add(value);error[p].add(uncertainty);
    }
    void rootError(unsigned p,const StripPolynomial& density,const Geometry::UnitPolynomialRoots::Root& root,long double divisor=1)
    {
        if(root.lower==root.upper)return;
        const auto range=Geometry::encloseUnitPolynomial(density,root.lower,root.upper);
        const auto bound=std::max(std::abs(range[0]),std::abs(range[1]));
        const auto value=(static_cast<long double>(root.upper-root.lower)*bound)/divisor;
        require(std::isfinite(value)&&(value>0||bound==0),Failure::VolumeResidual,-1,"unrepresentable smooth strip root uncertainty");
        if(value>0)error[p].add(std::nextafter(value,std::numeric_limits<long double>::infinity()));
    }
};
unsigned stripOuterPhase(const VofSmoothArc& arc,const Cut::Rectangle& box,bool reverse,bool far,bool liquidLower)
{
    const auto& endpoint=arc.controlPoints[(reverse!=far)?4:0];
    const auto lower=arc.graphAxis?box.lower.x:box.lower.y,upper=arc.graphAxis?box.upper.x:box.upper.y;
    const auto v=coordinate(endpoint,1-arc.graphAxis);
    require(v==lower||v==upper,Failure::InvalidCurve,arc.nativeCell,"outer smooth strip has no transverse native endpoint");
    return (v==upper)==liquidLower?0:1;
}
StripAccumulator smoothStrip(const VofSmoothArc& arc,const Cut::Rectangle& box,unsigned normalAxis,
    bool upperFace,long double tangentWidth,long double volume,long double tolerance,std::size_t& rootsCount)
{
    const unsigned u=arc.graphAxis,v=1-u;
    const auto ulo=u?box.lower.y:box.lower.x,uhi=u?box.upper.y:box.upper.x;
    const auto vlo=v?box.lower.y:box.lower.x,vhi=v?box.upper.y:box.upper.x;
    const bool increasing=coordinate(arc.controlPoints.back(),u)>coordinate(arc.controlPoints.front(),u);
    const bool liquidLower=u?increasing:!increasing;
    // For a graph-normal strip start at the endpoint nearest its face. A
    // sub-ulp inward width then produces a small parameter near zero, never
    // the rounded subtraction 1-t near the opposite parameter endpoint.
    const auto faceV=upperFace?vhi:vlo;
    const auto firstDistance=(coordinate(arc.controlPoints.front(),v)-faceV)*(upperFace?-1.L:1.L);
    const auto lastDistance=(coordinate(arc.controlPoints.back(),v)-faceV)*(upperFace?-1.L:1.L);
    const bool reverse=normalAxis==u?increasing==upperFace:
        firstDistance==lastDistance?!increasing:firstDistance>lastDistance;
    const long double direction=reverse?(increasing?-1.L:1.L):(increasing?1.L:-1.L);
    const auto derivative=stripDerivative(stripCoordinate(arc,u,reverse,0,direction));
    StripAccumulator result;
    if(normalAxis==u){
        const auto distance=stripCoordinate(arc,u,reverse,upperFace?uhi:ulo,upperFace?-1.L:1.L);
        const auto d0=Geometry::evaluateUnitPolynomial(distance,0),d1=Geometry::evaluateUnitPolynomial(distance,1);
        require(d0>=0&&d1>d0&&d1<=uhi-ulo,Failure::InvalidCurve,arc.nativeCell,"invalid face-anchored graph strip range");
        const auto threshold=stripDifference(stripConstant(volume),stripScale(distance,tangentWidth));
        const auto firstSign=Geometry::evaluateUnitPolynomial(threshold,0),lastSign=Geometry::evaluateUnitPolynomial(threshold,1);
        if(firstSign<=0){result.phase[stripOuterPhase(arc,box,reverse,false,liquidLower)].add(volume);return result;}
        if(d0>0)result.positive(stripOuterPhase(arc,box,reverse,false,liquidLower),
            stripScale(stripEndpoint(distance,false),tangentWidth),0,1,volume,tolerance);
        if(lastSign>0)result.positive(stripOuterPhase(arc,box,reverse,true,liquidLower),
            stripEndpoint(threshold,true),0,1,volume,tolerance);
        double end=1;std::optional<Geometry::UnitPolynomialRoots::Root> boundary;
        if(lastSign<0){const auto isolated=Geometry::isolateUnitPolynomial(threshold);
            require(!isolated.identicallyZero&&isolated.roots.size()==1,Failure::InvalidCurve,arc.nativeCell,"graph strip has no unique inward boundary");
            boundary=isolated.roots.front();end=boundary->point;++rootsCount;
            require(end>0&&end<1,Failure::VolumeResidual,arc.nativeCell,"unrepresentable positive graph strip parameter");}
        const auto below=stripProduct(stripCoordinate(arc,v,reverse,vlo),derivative);
        const auto above=stripProduct(stripCoordinate(arc,v,reverse,vhi,-1),derivative);
        const unsigned lowerPhase=liquidLower?0:1;
        result.positive(lowerPhase,below,0,end,volume,tolerance);result.positive(1-lowerPhase,above,0,end,volume,tolerance);
        if(boundary){result.rootError(lowerPhase,below,*boundary);result.rootError(1-lowerPhase,above,*boundary);}
    }else{
        const bool graphReverse=!increasing;
        const auto first=coordinate(arc.controlPoints[graphReverse?4:0],u),last=coordinate(arc.controlPoints[graphReverse?0:4],u);
        if(first>ulo)result.positive(stripOuterPhase(arc,box,graphReverse,false,liquidLower),
            stripScale(stripConstant(volume),first-ulo),0,1,volume,tolerance,tangentWidth);
        if(last<uhi)result.positive(stripOuterPhase(arc,box,graphReverse,true,liquidLower),
            stripScale(stripConstant(volume),uhi-last),0,1,volume,tolerance,tangentWidth);
        const auto inward=stripCoordinate(arc,v,reverse,upperFace?vhi:vlo,upperFace?-1.L:1.L);
        const auto inside=stripDifference(stripConstant(volume),stripScale(inward,tangentWidth));
        const auto nearDensity=stripProduct(inward,derivative),farDensity=stripProduct(inside,derivative);
        const auto fullDensity=stripProduct(stripConstant(volume),derivative);
        const unsigned nearPhase=(upperFace?!liquidLower:liquidLower)?0:1;
        const auto isolated=Geometry::isolateUnitPolynomial(inside);
        std::vector<Geometry::UnitPolynomialRoots::Root> cuts{{0,0,0},{1,1,1}};
        for(const auto& root:isolated.roots){
            require((root.point>0&&root.point<1)||root.lower==root.upper,Failure::VolumeResidual,arc.nativeCell,
                "smooth strip boundary lacks a representable endpoint-relative parameter");
            if(root.point<=0||root.point>=1)continue;cuts.push_back(root);++rootsCount;
            // The two neighboring formulas agree at the exact root. Only
            // their difference can contribute enclosure error, not a full
            // constant strip density on that tiny terminal range.
            result.rootError(0,farDensity,root,tangentWidth);result.rootError(1,farDensity,root,tangentWidth);}
        std::sort(cuts.begin(),cuts.end(),[](const auto& a,const auto& b){return a.point<b.point;});
        for(std::size_t k=1;k<cuts.size();++k){const auto& a=cuts[k-1];const auto& b=cuts[k];
            require(a.upper<b.lower,Failure::InvalidCurve,arc.nativeCell,"unresolved smooth strip boundary ordering");
            const auto middle=a.upper+(b.lower-a.upper)/2;
            require(middle>a.upper&&middle<b.lower,Failure::InvalidCurve,arc.nativeCell,"no certified smooth strip phase interval");
            if(Geometry::evaluateUnitPolynomial(inside,middle)>0){
                result.positive(nearPhase,nearDensity,a.point,b.point,volume,tolerance);
                result.positive(1-nearPhase,farDensity,a.point,b.point,volume,tolerance,tangentWidth);
            }else result.positive(nearPhase,fullDensity,a.point,b.point,volume,tolerance,tangentWidth);
        }
    }
    return result;
}
bool certified(const VofSmoothArc& arc,const Cut::Rectangle& box,unsigned& maximumDepth)
{
    // A whole-arc hull may exceed the box while its actual curve stays inside.
    // Fixed-depth de Casteljau subdivision tightens the SAME represented
    // curve, not its geometry policy. Every retained subarc must have a box
    // hull and a strictly signed derivative hull. At most 2^8 pieces are
    // examined; unresolved bounds reject. No sample test or box expansion.
    const auto sign=coordinate(arc.controlPoints.back(),arc.graphAxis)>
        coordinate(arc.controlPoints.front(),arc.graphAxis)?1.L:-1.L;
    struct Piece { std::array<Point,5> controls;unsigned depth=0; };
    std::array<Piece,9> stack{};unsigned count=1;stack[0]={arc.controlPoints,0};
    while(count){
        const auto piece=stack[--count];maximumDepth=std::max(maximumDepth,piece.depth);
        bool bounded=true,monotone=true;
        for(const auto p:piece.controls)
            bounded=bounded&&std::isfinite(p.x)&&std::isfinite(p.y)&&p.x>=box.lower.x&&p.x<=box.upper.x
                &&p.y>=box.lower.y&&p.y<=box.upper.y;
        for(unsigned i=1;i<5;++i)
            monotone=monotone&&sign*(coordinate(piece.controls[i],arc.graphAxis)-
                coordinate(piece.controls[i-1],arc.graphAxis))>0;
        if(bounded&&monotone)continue;
        if(piece.depth==8)return false;
        Piece left,right;left.depth=right.depth=piece.depth+1;
        auto q=piece.controls;left.controls[0]=q[0];right.controls[4]=q[4];
        for(unsigned level=1;level<=4;++level){
            for(unsigned i=0;i<5-level;++i)q[i]=mul(add(q[i],q[i+1]),.5L);
            left.controls[level]=q[0];right.controls[4-level]=q[4-level];
        }
        // This is an actual curve point, so an exterior midpoint cannot be
        // rescued by tighter enclosures. It is a rejection, never clipping.
        const auto middle=left.controls[4];
        if(!std::isfinite(middle.x)||!std::isfinite(middle.y)||middle.x<=box.lower.x||middle.x>=box.upper.x
            ||middle.y<=box.lower.y||middle.y>=box.upper.y)return false;
        stack[count++]=right;stack[count++]=left;
    }
    return true;
}
using DirectExact=Geometry::Detail::FiniteVolume::Dyadic;
long double directDown(long double x){return Geometry::nativeOutwardNeighbor(x,false);}
long double directUp(long double x){return Geometry::nativeOutwardNeighbor(x,true);}
struct DirectRange {long double lower=0,upper=0;};
DirectRange exactRange(const DirectExact& exact)
{
    const auto value=exact.rounded();const int order=DirectExact(value).compare(exact);
    return {order>0?directDown(value):value,order<0?directUp(value):value};
}
DirectRange addRange(DirectRange a,DirectRange b)
{
    if(a.lower==0&&a.upper==0)return b;if(b.lower==0&&b.upper==0)return a;
    return {directDown(a.lower+b.lower),directUp(a.upper+b.upper)};
}
DirectRange scaleRange(DirectRange a,long double scale)
{
    if(scale==0)return {};if(scale==1)return a;
    const auto x=a.lower*scale,y=a.upper*scale;
    return {directDown(std::min(x,y)),directUp(std::max(x,y))};
}
DirectRange differenceRange(DirectRange a,DirectRange b)
{return addRange(a,{-b.upper,-b.lower});}
DirectRange ratioRange(const DirectExact& numerator,const DirectExact& denominator,int native)
{
    const auto n=exactRange(numerator),d=exactRange(denominator);
    require(d.lower>0||d.upper<0,Failure::InfeasibleVolume,native,"direct XY quotient has no finite denominator enclosure");
    const std::array<long double,4> values{{n.lower/d.lower,n.lower/d.upper,n.upper/d.lower,n.upper/d.upper}};
    const auto bounds=std::minmax_element(values.begin(),values.end());
    require(std::isfinite(*bounds.first)&&std::isfinite(*bounds.second),Failure::InfeasibleVolume,native,
        "direct XY quotient exceeds its native finite domain");
    return {directDown(*bounds.first),directUp(*bounds.second)};
}
std::array<DirectExact,2> analyticArea70XY(const VofSmoothArc& arc,const GeometricVofPlicCell& cell)
{
    const auto& p=arc.controlPoints;const auto axis=arc.graphAxis,v=1-axis;
    require(axis<2,Failure::InvalidCurve,arc.nativeCell,"invalid direct XY graph axis");
    const auto low=axis?cell.rectangle.lower.x:cell.rectangle.lower.y;
    const auto high=axis?cell.rectangle.upper.x:cell.rectangle.upper.y;
    const auto graphLow=axis?cell.rectangle.lower.y:cell.rectangle.lower.x;
    const auto graphHigh=axis?cell.rectangle.upper.y:cell.rectangle.upper.x;
    const auto first=coordinate(p.front(),axis),last=coordinate(p.back(),axis);
    require(first!=last,Failure::InvalidCurve,arc.nativeCell,"direct XY curve has zero graph span");
    const auto direction=last>first?1.L:-1.L;const bool liquidLower=axis?direction>0:direction<0;
    const unsigned lowerPhase=liquidLower?0:1;
    constexpr unsigned c4[5]{1,4,6,4,1},c3[4]{1,3,3,1},c7[8]{1,7,21,35,35,21,7,1};
    std::array<DirectExact,2> result;
    for(unsigned i=0;i<5;++i)for(unsigned j=0;j<4;++j){
        // 70*integral(B_i^4 * 4 B_j^3) is an integer. Original
        // coordinate factors remain separate until the exact contraction.
        const auto k=direction*(35*c4[i]*c3[j]/c7[i+j]);
        const auto value=coordinate(p[i],v),a=coordinate(p[j],axis),b=coordinate(p[j+1],axis);
        auto& below=result[lowerPhase];auto& above=result[1-lowerPhase];
        below.add(DirectExact::product({k,value,b}));below.add(DirectExact::product({k,value,a}),-1);
        below.add(DirectExact::product({k,low,b}),-1);below.add(DirectExact::product({k,low,a}));
        above.add(DirectExact::product({k,high,b}));above.add(DirectExact::product({k,high,a}),-1);
        above.add(DirectExact::product({k,value,b}),-1);above.add(DirectExact::product({k,value,a}));
    }
    const auto outer=[&](long double begin,long double end){if(!(end>begin))return;
        auto& area=result[graphOuterPhase(arc,cell.rectangle,begin+(end-begin)/2,liquidLower)];
        area.add(DirectExact::product({70,end,high}));area.add(DirectExact::product({70,end,low}),-1);
        area.add(DirectExact::product({70,begin,high}),-1);area.add(DirectExact::product({70,begin,low}));};
    outer(graphLow,std::min(first,last));outer(std::max(first,last),graphHigh);
    return result;
}
std::array<DirectExact,2> liquidAreaGradient70(const VofSmoothArc& arc)
{
    std::array<DirectExact,2> k;
    for(unsigned a=0;a<2;++a){k[a].add(DirectExact::product({5,coordinate(arc.controlPoints[4],a)}));
        k[a].add(DirectExact::product({5,coordinate(arc.controlPoints[0],a)}),-1);
        k[a].add(DirectExact::product({4,coordinate(arc.controlPoints[3],a)}));
        k[a].add(DirectExact::product({4,coordinate(arc.controlPoints[1],a)}),-1);}
    return {{k[1].multiplied(2),k[0].multiplied(-2)}};
}
std::array<DirectRange,2> directCenterDomain(const VofSmoothArc& arc,const Cut::Rectangle& box,
    Point anchor,Detail::SmoothDirectXYWork& work)
{
    using Row=std::array<long double,5>;
    struct Piece {std::array<Row,5> weights{};unsigned depth=0;};
    std::array<std::array<DirectRange,2>,5> relative{};std::array<DirectRange,2> low{},high{};
    for(unsigned a=0;a<2;++a){DirectExact lower(a?box.lower.y:box.lower.x),upper(a?box.upper.y:box.upper.x);
        lower.add(coordinate(anchor,a),-1);upper.add(coordinate(anchor,a),-1);low[a]=exactRange(lower);high[a]=exactRange(upper);
        for(unsigned i:{0U,1U,3U,4U}){DirectExact value(coordinate(arc.controlPoints[i],a));
            value.add(coordinate(anchor,a),-1);relative[i][a]=exactRange(value);}}
    std::array<DirectRange,2> domain{{{-std::numeric_limits<long double>::infinity(),std::numeric_limits<long double>::infinity()},
        {-std::numeric_limits<long double>::infinity(),std::numeric_limits<long double>::infinity()}}};
    const auto constrain=[&](unsigned axis,DirectRange constant,long double slope,bool strict){
        require(std::isfinite(constant.lower)&&std::isfinite(constant.upper)&&std::isfinite(slope),
            Failure::InfeasibleVolume,arc.nativeCell,"nonfinite direct XY geometry constraint");
        if(slope==0){require(strict?constant.lower>0:constant.lower>=0,Failure::InfeasibleVolume,arc.nativeCell,
            "fixed endpoint jet has no direct XY hull certificate");return;}
        const auto boundary=-constant.lower/slope;require(std::isfinite(boundary),Failure::InfeasibleVolume,arc.nativeCell,
            "direct XY hull boundary exceeds the finite native domain");
        if(slope>0)domain[axis].lower=std::max(domain[axis].lower,directUp(boundary));
        else domain[axis].upper=std::min(domain[axis].upper,directDown(boundary));
        require(domain[axis].lower<domain[axis].upper,Failure::InfeasibleVolume,arc.nativeCell,
            "empty direct XY sufficient center rectangle");
    };
    const auto constant=[&](const Row& row,unsigned axis){DirectRange value;
        for(unsigned i:{0U,1U,3U,4U})value=addRange(value,scaleRange(relative[i][axis],row[i]));return value;};
    Piece initial;for(unsigned i=0;i<5;++i)initial.weights[i][i]=1;
    std::array<Piece,9> stack{};unsigned count=1;stack[0]=initial;
    const auto direction=coordinate(arc.controlPoints[4],arc.graphAxis)>coordinate(arc.controlPoints[0],arc.graphAxis)?1.L:-1.L;
    while(count){const auto piece=stack[--count];
        if(piece.depth==8){++work.domainLeaves;
            for(unsigned i=0;i<5;++i){const auto& row=piece.weights[i];const bool endpoint=row[0]==1||row[4]==1;
                if(endpoint){const auto& p=arc.controlPoints[row[0]==1?0:4];
                    require(p.x>=box.lower.x&&p.x<=box.upper.x&&p.y>=box.lower.y&&p.y<=box.upper.y,
                        Failure::InfeasibleVolume,arc.nativeCell,"direct XY endpoint is outside its original cell");continue;}
                const bool strict=i==0||i==4;
                for(unsigned a=0;a<2;++a){const auto c=constant(row,a);
                    constrain(a,differenceRange(c,low[a]),row[2],strict);
                    constrain(a,differenceRange(high[a],c),-row[2],strict);}}
            for(unsigned i=1;i<5;++i){Row difference;for(unsigned j=0;j<5;++j)
                    difference[j]=direction*(piece.weights[i][j]-piece.weights[i-1][j]);
                constrain(arc.graphAxis,constant(difference,arc.graphAxis),difference[2],true);}
            continue;}
        Piece left,right;left.depth=right.depth=piece.depth+1;auto rows=piece.weights;
        left.weights[0]=rows[0];right.weights[4]=rows[4];
        for(unsigned level=1;level<=4;++level){for(unsigned i=0;i<5-level;++i)for(unsigned j=0;j<5;++j)
                rows[i][j]=(rows[i][j]+rows[i+1][j])/2;
            left.weights[level]=rows[0];right.weights[4-level]=rows[4-level];}
        stack[count++]=right;stack[count++]=left;
    }
    for(const auto interval:domain)require(std::isfinite(interval.lower)&&std::isfinite(interval.upper)&&interval.lower<interval.upper,
        Failure::InfeasibleVolume,arc.nativeCell,"direct XY rectangle has no finite interior");
    return domain;
}
Detail::SmoothDirectXYResult directSmoothXY(const VofSmoothArc& base,const GeometricVofPlicCell& cell,
    long double tolerance,Detail::SmoothDirectXYWork& work)
{
    const auto native=base.nativeCell;
    require(base.graphAxis<2&&cell.inventory[0]>0&&cell.inventory[1]>0&&std::isfinite(tolerance)&&tolerance>0&&tolerance<=1e-12L,
        Failure::InfeasibleVolume,native,"invalid direct XY source or tolerance");
    const auto areas=analyticArea70XY(base,cell);const auto liquidGradient=liquidAreaGradient70(base);
    const unsigned minor=cell.inventory[0]<cell.inventory[1]?0:1;
    std::array<DirectExact,2> gradient=liquidGradient;if(minor==1)for(auto& g:gradient)g=g.multiplied(-1);
    Point anchor;for(unsigned a=0;a<2;++a){const auto lo=a?cell.rectangle.lower.y:cell.rectangle.lower.x;
        const auto hi=a?cell.rectangle.upper.y:cell.rectangle.upper.x;const auto value=coordinate(base.controlPoints[2],a);
        if(a==0)anchor.x=std::abs(value-lo)<=std::abs(hi-value)?lo:hi;
        else anchor.y=std::abs(value-lo)<=std::abs(hi-value)?lo:hi;}
    const auto domain=directCenterDomain(base,cell.rectangle,anchor,work);
    auto constant=areas[minor];for(unsigned a=0;a<2;++a){constant.add(gradient[a].multiplied(coordinate(anchor,a)));
        constant.add(gradient[a].multiplied(coordinate(base.controlPoints[2],a)),-1);}
    const auto areaAt=[&](long double x,long double y){auto value=constant;
        value.add(gradient[0].multiplied(x));value.add(gradient[1].multiplied(y));return value;};
    auto minimum=areaAt(domain[0].lower,domain[1].lower),maximum=minimum;
    for(const auto x:{domain[0].lower,domain[0].upper})for(const auto y:{domain[1].lower,domain[1].upper}){
        auto value=areaAt(x,y);if(value.compare(minimum)<0)minimum=value;if(value.compare(maximum)>0)maximum=value;}
    const auto minorStock=DirectExact::product({70,cell.inventory[minor]});
    const auto majorStock=DirectExact::product({70,cell.inventory[1-minor]});
    const auto minorSlack=DirectExact::product({70,tolerance,cell.inventory[minor]});
    const auto majorSlack=DirectExact::product({70,tolerance,cell.inventory[1-minor]});
    auto minorLower=minorStock,minorUpper=minorStock;minorLower.add(minorSlack,-1);minorUpper.add(minorSlack);
    auto majorLower=majorStock,majorUpper=majorStock;majorLower.add(majorSlack,-1);majorUpper.add(majorSlack);
    // The actual independent integrals define the exact total. Do not create
    // another total by rounding capacity or summing the independent stocks.
    auto capacity=areas[0];capacity.add(areas[1]);
    auto fromMajorLower=capacity,fromMajorUpper=capacity;
    fromMajorLower.add(majorUpper,-1);fromMajorUpper.add(majorLower,-1);
    auto lower=minimum,upper=maximum;
    for(const auto& bound:{minorLower,fromMajorLower})if(bound.compare(lower)>0)lower=bound;
    for(const auto& bound:{minorUpper,fromMajorUpper})if(bound.compare(upper)<0)upper=bound;
    require(lower.compare(upper)<=0,Failure::InfeasibleVolume,native,
        "direct XY sufficient geometry and independent stock bands do not intersect");
    unsigned pivot=2;DirectExact bestQuantum;
    for(unsigned a=0;a<2;++a){if(gradient[a].compare(DirectExact{})==0)continue;
        const auto magnitude=std::max(std::abs(coordinate(anchor,a)+domain[a].lower),std::abs(coordinate(anchor,a)+domain[a].upper));
        const auto quantum=directUp(magnitude)-magnitude;
        require(std::isfinite(quantum)&&quantum>0,Failure::InfeasibleVolume,native,"unrepresentable direct XY world quantum");
        const auto error=gradient[a].absolute().multiplied(quantum);
        if(pivot==2||error.compare(bestQuantum)<0){pivot=a;bestQuantum=error;}}
    require(pivot<2,Failure::InfeasibleVolume,native,"direct XY area has zero center gradient");
    auto target=minorStock;
    if(target.compare(lower)<0||target.compare(upper)>0){
        // Only an exact-minor target outside the JOINT band needs moving.
        // Reserve eight a-priori world-area quanta plus 64 native eps of the
        // smaller stock; this is an inward selection budget, not a looser
        // acceptance gate or a universal factor/quotient error theorem.
        // Insufficient reserve rejects before constructing any candidate.
        auto reserve=bestQuantum.multiplied(8);
        reserve.add(DirectExact::product({70,64,std::numeric_limits<long double>::epsilon(),
            std::min(cell.inventory[0],cell.inventory[1])}));
        auto safeLower=lower,safeUpper=upper;safeLower.add(reserve);safeUpper.add(reserve,-1);
        require(safeLower.compare(safeUpper)<=0,Failure::InfeasibleVolume,native,
            "direct XY joint stock band has insufficient fixed native selection reserve");
        if(target.compare(safeLower)<0)target=safeLower;else if(target.compare(safeUpper)>0)target=safeUpper;
    }
    // Verify both original bands even when their independent stocks do not
    // sum to the exact geometric total. Neither stock is normalized/repaired.
    auto selected=areas;auto change=target;change.add(areas[minor],-1);selected[minor]=target;selected[1-minor].add(change,-1);
    for(unsigned p=0;p<2;++p){auto error=selected[p];error.add(DirectExact::product({70,cell.inventory[p]}),-1);
        require(error.absolute().compare(DirectExact::product({70,tolerance,cell.inventory[p]}))<=0,
            Failure::InfeasibleVolume,native,"direct XY selected area is outside an original independent stock band");}
    const unsigned free=1-pivot;auto freeInterval=domain[free];
    auto intercept=target;intercept.add(constant,-1);
    if(gradient[free].compare(DirectExact{})!=0){auto n0=intercept,n1=intercept;
        n0.add(gradient[pivot].multiplied(domain[pivot].lower),-1);
        n1.add(gradient[pivot].multiplied(domain[pivot].upper),-1);
        auto r0=ratioRange(n0,gradient[free],native),r1=ratioRange(n1,gradient[free],native);
        if((n0.compare(n1)>0)==(gradient[free].compare(DirectExact{})>0))std::swap(r0,r1);
        freeInterval.lower=std::max(freeInterval.lower,r0.upper);freeInterval.upper=std::min(freeInterval.upper,r1.lower);}
    require(freeInterval.lower<=freeInterval.upper,Failure::InfeasibleVolume,native,"direct XY stock equality has no native-safe segment");
    const std::array<long double,2> g{{gradient[0].rounded(),gradient[1].rounded()}};
    const auto magnitude=std::max(std::abs(g[0]),std::abs(g[1]));
    require(magnitude>0&&std::isfinite(magnitude)&&g[pivot]!=0,Failure::InfeasibleVolume,native,"direct XY native area gradient is unavailable");
    const std::array<long double,2> h{{static_cast<long double>(cell.rectangle.upper.x)-cell.rectangle.lower.x,
        static_cast<long double>(cell.rectangle.upper.y)-cell.rectangle.lower.y}};
    const auto hmax=std::max(h[0],h[1]);const std::array<long double,2> w{{h[0]/hmax,h[1]/hmax}};
    const std::array<long double,2> scaled{{g[0]/magnitude,g[1]/magnitude}};
    const auto denominator=w[0]*w[0]*scaled[0]*scaled[0]+w[1]*w[1]*scaled[1]*scaled[1];
    auto residual=target;residual.add(areas[minor],-1);
    const auto nominal=coordinate(base.controlPoints[2],free)-coordinate(anchor,free);
    const auto projected=nominal+w[free]*w[free]*scaled[free]/denominator*(residual.rounded()/magnitude);
    require(std::isfinite(projected),Failure::InfeasibleVolume,native,"nonfinite direct XY metric projection");
    const auto localFree=std::clamp(projected,freeInterval.lower,freeInterval.upper);
    const auto worldFree=coordinate(anchor,free)+localFree;
    // The single pivot quotient targets the ACTUAL already-rounded free world
    // coordinate. There is no sample/measure correction or adjacent-value try.
    auto numerator=target;numerator.add(areas[minor],-1);
    for(unsigned a=0;a<2;++a)numerator.add(gradient[a].multiplied(coordinate(base.controlPoints[2],a)));
    numerator.add(gradient[free].multiplied(worldFree),-1);
    const auto worldPivot=numerator.rounded()/g[pivot];
    require(std::isfinite(worldFree)&&std::isfinite(worldPivot),Failure::InfeasibleVolume,native,"direct XY candidate is not world-native representable");
    Detail::SmoothDirectXYResult result;result.arc=base;
    if(pivot==0)result.arc.controlPoints[2]={worldPivot,worldFree};else result.arc.controlPoints[2]={worldFree,worldPivot};
    ++work.candidateConstructions;++work.geometryCertifications;
    const auto rejectCandidate=[&](const char* message){std::ostringstream detail;detail<<std::setprecision(17)<<message<<" controls=";
        for(const auto point:result.arc.controlPoints)detail<<'('<<point.x<<','<<point.y<<')';
        detail<<" measured="<<result.measured[0]<<','<<result.measured[1]<<" stock="<<cell.inventory[0]<<','<<cell.inventory[1];
        throw VofSmoothGeometryError(Failure::InfeasibleVolume,native,detail.str().c_str());};
    if(!certified(result.arc,cell.rectangle,result.maximumCertificateSubdivisionDepth))
        rejectCandidate("single direct XY world candidate fails the original geometry certificate");
    auto actualArea=areas;
    for(unsigned a=0;a<2;++a){auto change=liquidGradient[a].multiplied(coordinate(result.arc.controlPoints[2],a));
        change.add(liquidGradient[a].multiplied(coordinate(base.controlPoints[2],a)),-1);
        actualArea[0].add(change);actualArea[1].add(change,-1);}
    for(unsigned p=0;p<2;++p){auto error=actualArea[p];error.add(DirectExact::product({70,cell.inventory[p]}),-1);
        if(error.absolute().compare(DirectExact::product({70,tolerance,cell.inventory[p]}))>0){
            result.measured={{actualArea[0].rounded()/70,actualArea[1].rounded()/70}};
            rejectCandidate("single direct XY world candidate lacks an exact own-phase area certificate");}}
    ++work.phaseQuadratureBuilds;result.measured=phaseMeasure(result.arc,cell,false,&result.phaseQuadrature);
    for(unsigned p=0;p<2;++p)if(!(std::abs(result.measured[p]-cell.inventory[p])<=tolerance*cell.inventory[p]))
        rejectCandidate("single direct XY world candidate fails an original own-phase gate");
    const auto direction=coordinate(base.controlPoints[4],base.graphAxis)>coordinate(base.controlPoints[0],base.graphAxis)?1.L:-1.L;
    const auto projectionSign=base.graphAxis==0?-direction:direction;
    result.arc.volumeCorrection=base.volumeCorrection+6*projectionSign*(coordinate(result.arc.controlPoints[2],1-base.graphAxis)
        -coordinate(base.controlPoints[2],1-base.graphAxis));
    powerProjection(result.arc);return result;
}

using CellSides=std::vector<std::array<int,4>>;
void appendSourceReceipt(std::ostream& detail,const GeometricVofPublication& source)
{
    const auto clock=source.timeLayer();
    detail<<" source_identity="<<source.identity().value()<<" source_clock=("<<clock.step<<','
        <<clock.begin<<','<<clock.end<<','<<static_cast<int>(clock.position)<<')';
}
void appendCellReceipt(std::ostream& detail,const std::vector<GeometricVofPlicCell>& cells,int native,bool rz)
{
    const auto& cell=cells.at(static_cast<std::size_t>(native));
    const auto capacity=Cut::capacity(cell.rectangle,rz);
    detail<<" cell_receipt=(native="<<native<<",box="<<cell.rectangle.lower.x<<','<<cell.rectangle.lower.y
        <<':'<<cell.rectangle.upper.x<<','<<cell.rectangle.upper.y<<",stock="<<cell.inventory[0]<<','<<cell.inventory[1]
        <<",capacity="<<capacity<<",minor_relative_capacity="<<std::min(cell.inventory[0],cell.inventory[1])/capacity
        <<",source_measured="<<cell.measured[0]<<','<<cell.measured[1]<<",source_resolved="<<cell.resolved
        <<",source_normal="<<cell.normal.x<<','<<cell.normal.y<<",source_offset="<<cell.offset<<')';
}
int shiftCell(int cell,unsigned axis,int amount,const CellSides& neighbors)
{
    const int side=static_cast<int>(2*axis)+(amount>0?1:0);
    for(int count=std::abs(amount);count>0&&cell>=0;--count)cell=neighbors[static_cast<std::size_t>(cell)][static_cast<std::size_t>(side)];
    return cell;
}
int stencilCell(int cell,unsigned normalAxis,int normal,int tangent,const CellSides& neighbors)
{ return shiftCell(shiftCell(cell,normalAxis,normal,neighbors),1-normalAxis,tangent,neighbors); }
int seedFacePhase(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,int row,
    const std::vector<GeometricVofPlicCell>& seed)
{
    const auto f=domain.face(mesh,row);Point a{f.center.x,f.center.y},b=a;
    if(f.horizontal){a.x-=f.tangentSpacing/2;b.x+=f.tangentSpacing/2;}
    else {a.y-=f.tangentSpacing/2;b.y+=f.tangentSpacing/2;}
    int phase=-1;
    for(const auto ci:{f.negativeCell,f.positiveCell})if(ci>=0){
        const auto& cell=seed[static_cast<std::size_t>(ci)];int current=-1;
        if(cell.inventory[0]==0)current=1;
        else if(cell.inventory[1]==0)current=0;
        else {const auto c=midpoint(cell.rectangle);
            const auto value=[&](Point p){return cell.normal.x*(p.x-c.x)+cell.normal.y*(p.y-c.y)-cell.offset;};
            const auto va=value(a),vb=value(b);
            if(va<0&&vb<0)current=0;else if(va>0&&vb>0)current=1;
        }
        if(current<0||(phase>=0&&phase!=current))return -1;phase=current;
    }
    return phase;
}
struct HeightSeed
{
    unsigned normalAxis=0;
    bool liquidLower=false,radialSquared=false;
    long double normalOrigin=0,tangentOrigin=0,tangentScale=1;
    std::array<long double,5> polynomial{};
    unsigned strips=3;
    std::vector<std::pair<int,unsigned>> endFaces;
    int tangentOffset=0,normalLower=0,normalUpper=0;
    long double value(long double tangent) const
    {
        const auto u=(tangent-tangentOrigin)/tangentScale;
        auto f=polynomial[4];for(int k=3;k>=0;--k)f=f*u+polynomial[static_cast<std::size_t>(k)];
        if(!radialSquared)return normalOrigin+f;
        const auto radicand=normalOrigin*normalOrigin+2*normalOrigin*f;
        if(!(radicand>0))return std::numeric_limits<long double>::quiet_NaN();
        const auto r=std::sqrt(radicand);
        return normalOrigin+(2*normalOrigin*f)/(r+normalOrigin);
    }
    Point normal(long double tangent,int cell) const
    {
        const auto u=(tangent-tangentOrigin)/tangentScale;
        auto slope=4*polynomial[4];for(int k=3;k>=1;--k)slope=slope*u+k*polynomial[static_cast<std::size_t>(k)];
        slope/=tangentScale;
        if(radialSquared)slope*=normalOrigin/value(tangent);
        const auto sign=liquidLower?1.L:-1.L;
        return unit(normalAxis?Point{-sign*slope,sign}:Point{sign,-sign*slope},cell);
    }
};
struct HeightEndpoint {Point point,normal;int face=-1;};
std::array<HeightEndpoint,2> heightEndpoints(const HeightSeed&,const Cut::Rectangle&,
    int,const std::array<int,4>&);
std::array<long double,5> solveHeight(std::array<std::array<long double,6>,5> rows,unsigned count)
{
    for(unsigned k=0;k<count;++k){
        unsigned pivot=k;for(unsigned i=k+1;i<count;++i)if(std::abs(rows[i][k])>std::abs(rows[pivot][k]))pivot=i;
        require(rows[pivot][k]!=0,Failure::InvalidCurve,-1,"dependent cell-average height equations");
        if(pivot!=k)std::swap(rows[pivot],rows[k]);
        for(unsigned i=k+1;i<count;++i){const auto factor=rows[i][k]/rows[k][k];
            for(unsigned j=k;j<count;++j)rows[i][j]-=factor*rows[k][j];
            rows[i][5]-=factor*rows[k][5];}
    }
    std::array<long double,5> result{};
    for(int i=static_cast<int>(count)-1;i>=0;--i){auto value=rows[static_cast<std::size_t>(i)][5];
        for(unsigned j=static_cast<unsigned>(i)+1;j<count;++j)value-=rows[static_cast<std::size_t>(i)][j]*result[j];
        result[static_cast<std::size_t>(i)]=value/rows[static_cast<std::size_t>(i)][static_cast<std::size_t>(i)];
        require(std::isfinite(result[static_cast<std::size_t>(i)]),Failure::InvalidCurve,-1,"nonfinite cell-average height solve");}
    return result;
}
long double binomial(unsigned degree,unsigned index)
{
    long double value=1;for(unsigned k=1;k<=index;++k)value=value*(degree+1-k)/k;return value;
}
std::optional<HeightSeed> heightSeed(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,
    int native,unsigned axis,const std::vector<GeometricVofPlicCell>& seed,
    const CellSides& neighbors,const CellSides& faceRows,const std::vector<int>& faceStatus)
{
    const auto& box=seed[static_cast<std::size_t>(native)].rectangle;const auto center=midpoint(box);
    const bool rz=domain.axisymmetric();
    unsigned available=1;
    for(const int direction:{-1,1}){int ci=native;for(unsigned k=0;k<4;++k){ci=shiftCell(ci,1-axis,direction,neighbors);if(ci<0)break;++available;}}
    if(available<3)return std::nullopt;
    const unsigned strips=available>=5?5:3;
    const std::vector<int> starts=strips==5?std::vector<int>{-2,-3,-1,-4,0}:std::vector<int>{-1,-2,0};
    // All choices are finite, ordered seed qualifications, before any actual
    // Hermite/root attempt. There is no retry with another chart after failure.
    for(const int tangentStart:starts)for(int width=1;width<=7;++width)
        for(int low=-3;low<=0;++low){const int high=low+width-1;if(high<0||high>3)continue;
            HeightSeed result;result.normalAxis=axis;result.radialSquared=rz&&axis==0;
            result.strips=strips;result.endFaces.resize(2*strips);
            result.normalOrigin=coordinate(center,axis);result.tangentOrigin=coordinate(center,1-axis);
            result.tangentScale=axis?static_cast<long double>(box.upper.x)-box.lower.x:
                static_cast<long double>(box.upper.y)-box.lower.y;
            result.tangentOffset=tangentStart;result.normalLower=low;result.normalUpper=high;
            std::array<std::array<long double,6>,5> equations{};int lowerPhase=-1;bool qualified=true;
            for(unsigned column=0;column<strips&&qualified;++column){
                const int tangent=tangentStart+static_cast<int>(column);
                const int lower=stencilCell(native,axis,low,tangent,neighbors),upper=stencilCell(native,axis,high,tangent,neighbors);
                if(lower<0||upper<0){qualified=false;break;}
                const auto lowerRow=faceRows[static_cast<std::size_t>(lower)][2*axis];
                const auto upperRow=faceRows[static_cast<std::size_t>(upper)][2*axis+1];
                if(faceStatus[static_cast<std::size_t>(lowerRow)]==1||faceStatus[static_cast<std::size_t>(upperRow)]==1){qualified=false;break;}
                const auto lp=seedFacePhase(mesh,domain,lowerRow,seed),up=seedFacePhase(mesh,domain,upperRow,seed);
                if(lp<0||up<0||lp==up||(lowerPhase>=0&&lp!=lowerPhase)){qualified=false;break;}
                lowerPhase=lp;result.liquidLower=lp==0;
                result.endFaces[static_cast<std::size_t>(2*column)]={lowerRow,static_cast<unsigned>(lp)};
                result.endFaces[static_cast<std::size_t>(2*column+1)]={upperRow,static_cast<unsigned>(up)};
                std::array<long double,2> stocks{};long double previous=result.liquidLower?1.L:0.L;
                for(int normal=low;normal<=high;++normal){const int ci=stencilCell(native,axis,normal,tangent,neighbors);
                    if(ci<0){qualified=false;break;}const auto& cell=seed[static_cast<std::size_t>(ci)];
                    const auto fraction=cell.inventory[0]==0?0.L:cell.inventory[1]==0?1.L:
                        cell.inventory[0]/Cut::capacity(cell.rectangle,rz);
                    if(result.liquidLower?fraction>previous:fraction<previous){qualified=false;break;}
                    previous=fraction;for(unsigned p=0;p<2;++p)stocks[p]+=cell.inventory[p];
                }
                if(!qualified)break;
                const auto& first=seed[static_cast<std::size_t>(lower)].rectangle;
                const auto& last=seed[static_cast<std::size_t>(upper)].rectangle;
                const auto normalLow=axis?static_cast<long double>(first.lower.y):first.lower.x;
                const auto normalHigh=axis?static_cast<long double>(last.upper.y):last.upper.x;
                const auto tangentLow=axis?static_cast<long double>(first.lower.x):first.lower.y;
                const auto tangentHigh=axis?static_cast<long double>(first.upper.x):first.upper.y;
                const auto tangentWidth=tangentHigh-tangentLow,tangentCenter=(tangentHigh+tangentLow)/2;
                auto liquidHeight=0.L,gasHeight=0.L;
                if(result.radialSquared){
                    const auto denominator=std::acos(-1.L)*tangentWidth;
                    const auto lowMoment=(normalLow-result.normalOrigin)*(normalLow+result.normalOrigin);
                    const auto highMoment=(normalHigh-result.normalOrigin)*(normalHigh+result.normalOrigin);
                    liquidHeight=(result.liquidLower?lowMoment+stocks[0]/denominator:highMoment-stocks[0]/denominator)/(2*result.normalOrigin);
                    gasHeight=(result.liquidLower?highMoment-stocks[1]/denominator:lowMoment+stocks[1]/denominator)/(2*result.normalOrigin);
                }else {
                    const auto measure=rz?std::acos(-1.L)*tangentWidth*(tangentHigh+tangentLow):tangentWidth;
                    liquidHeight=(result.liquidLower?normalLow+stocks[0]/measure:normalHigh-stocks[0]/measure)-result.normalOrigin;
                    gasHeight=(result.liquidLower?normalHigh-stocks[1]/measure:normalLow+stocks[1]/measure)-result.normalOrigin;
                }
                const auto d=(tangentCenter-result.tangentOrigin)/result.tangentScale;
                const auto halfWidth=tangentWidth/(2*result.tangentScale);
                for(unsigned degree=0;degree<strips;++degree){long double moment=0;
                    for(unsigned j=0;j<=degree;++j){const auto factor=binomial(degree,j)*std::pow(d,static_cast<int>(degree-j));
                        if(j%2==0)moment+=factor*std::pow(halfWidth,static_cast<int>(j))/(j+1);
                        if(rz&&axis==1&&(j+1)%2==0)
                            moment+=factor*result.tangentScale/tangentCenter*std::pow(halfWidth,static_cast<int>(j+1))/(j+2);
                    }
                    equations[column][degree]=moment;
                }
                equations[column][5]=(liquidHeight+gasHeight)/2;
            }
            if(qualified){
                result.polynomial=solveHeight(equations,strips);
                try {
                    const auto ends=heightEndpoints(result,box,native,faceRows[static_cast<std::size_t>(native)]);
                    for(const auto& end:ends)qualified=qualified&&faceStatus[static_cast<std::size_t>(end.face)]>=0;
                    for(const auto row:faceRows[static_cast<std::size_t>(native)])
                        if(faceStatus[static_cast<std::size_t>(row)]==1)
                            qualified=qualified&&(ends[0].face==row||ends[1].face==row);
                    if(qualified)return result;
                }catch(const VofSmoothGeometryError&){
                    // A seed chart is qualified before selection, never used
                    // as a fallback for a failed final curve or volume root.
                }
            }
        }
    return std::nullopt;
}
std::array<HeightEndpoint,2> heightEndpoints(const HeightSeed& height,const Cut::Rectangle& box,
    int native,const std::array<int,4>& rows)
{
    const unsigned axis=height.normalAxis,tangentAxis=1-axis;
    const long double tangentLow=axis?box.lower.x:box.lower.y,tangentHigh=axis?box.upper.x:box.upper.y;
    const long double normalLow=axis?box.lower.y:box.lower.x,normalHigh=axis?box.upper.y:box.upper.x;
    std::vector<HeightEndpoint> points;
    const auto append=[&](long double tangent,long double normal,int face){
        require((tangent>tangentLow&&tangent<tangentHigh)||(normal>normalLow&&normal<normalHigh),
            Failure::CornerEndpoint,native,"height seed meets a native corner");
        const Point p=axis?Point{tangent,normal}:Point{normal,tangent};
        for(const auto& existing:points)if(existing.point.x==p.x&&existing.point.y==p.y)return;
        points.push_back({p,height.normal(tangent,native),face});
    };
    for(unsigned side=0;side<2;++side){const auto t=side?tangentHigh:tangentLow,z=height.value(t);
        require(std::isfinite(z),Failure::InvalidCurve,native,"height seed has nonpositive radial square");
        if(z>=normalLow&&z<=normalHigh)append(t,z,rows[2*tangentAxis+side]);
    }
    const auto uLow=(tangentLow-height.tangentOrigin)/height.tangentScale;
    const auto uHigh=(tangentHigh-height.tangentOrigin)/height.tangentScale;
    for(unsigned side=0;side<2;++side){const auto boundary=side?normalHigh:normalLow;
        const auto target=height.radialSquared?(boundary-height.normalOrigin)*(boundary+height.normalOrigin)/(2*height.normalOrigin):
            boundary-height.normalOrigin;
        std::vector<long double> roots;
        if(height.strips==3){
            const auto a=height.polynomial[2],b=height.polynomial[1],c=height.polynomial[0]-target;
            if(a==0){if(b!=0)roots.push_back(-c/b);else require(c!=0,Failure::CornerEndpoint,native,"height seed coincides with a native face");}
            else {const auto discriminant=b*b-4*a*c;if(discriminant>=0){
                const auto q=-.5L*(b+std::copysign(std::sqrt(discriminant),b));
                if(q!=0){roots.push_back(q/a);roots.push_back(c/q);}else roots.push_back(-b/(2*a));}}
        }else {
            // Root the same retained quartic through the public bounded Core
            // engine. Terms remain separate, with no small-coefficient deletion.
            std::vector<Geometry::UnitPolynomialTerm> terms;
            for(unsigned k=0;k<5;++k){Geometry::UnitPolynomialTerm term;term.factor=height.polynomial[k];
                for(unsigned j=0;j<=k;++j)term.coefficients[j]=static_cast<double>(binomial(k,j)*
                    std::pow(uLow,static_cast<int>(k-j))*std::pow(uHigh-uLow,static_cast<int>(j)));
                terms.push_back(term);
            }
            Geometry::UnitPolynomialTerm constant;constant.factor=-target;constant.coefficients[0]=1;terms.push_back(constant);
            const auto isolated=Geometry::isolateUnitPolynomial(terms);
            require(!isolated.identicallyZero,Failure::CornerEndpoint,native,"quartic height seed coincides with a native face");
            for(const auto& root:isolated.roots)roots.push_back(uLow+(uHigh-uLow)*root.point);
        }
        for(const auto u:roots)if(u>=uLow&&u<=uHigh)
            append(height.tangentOrigin+height.tangentScale*u,boundary,rows[2*axis+side]);
    }
    require(points.size()==2,Failure::UnmatchedEndpoint,native,"height seed does not define exactly one regular cell crossing");
    return {points[0],points[1]};
}
void pureQuadrature(VofSmoothCell& cell,bool rz,unsigned phase)
{
    const auto& box=cell.rectangle;
    const auto dx=static_cast<long double>(box.upper.x)-box.lower.x;
    const auto dy=static_cast<long double>(box.upper.y)-box.lower.y;
    for(unsigned i=0;i<4;++i)for(unsigned j=0;j<4;++j){
        const Point p{(1-eta[i])*box.lower.x+eta[i]*box.upper.x,
            (1-eta[j])*box.lower.y+eta[j]*box.upper.y};
        const auto weight=dx*dy*etaWeight[i]*etaWeight[j]*(rz?2*std::acos(-1.L)*p.x:1);
        require(weight>0&&std::isfinite(weight),Failure::InvalidCurve,cell.nativeCell,"nonpositive pure-cell quadrature");
        cell.phaseQuadrature[phase].push_back({p,weight});cell.measured[phase]+=weight;
    }
}
}

std::array<long double,2> Detail::measureSmoothArcPhases(const VofSmoothArc& arc,
    const GeometricVofPlicCell& cell,bool rz,
    std::array<std::vector<VofSmoothWeightedPoint>,2>* output)
{ return phaseMeasure(arc,cell,rz,output); }
bool Detail::certifySmoothArc(const VofSmoothArc& arc,const Cut::Rectangle& box,unsigned& maximumDepth)
{ return certified(arc,box,maximumDepth); }
std::array<long double,2> Detail::analyticSmoothArcPhasesXY(const VofSmoothArc& arc,const GeometricVofPlicCell& cell)
{ const auto exact=analyticArea70XY(arc,cell);return {{exact[0].rounded()/70,exact[1].rounded()/70}}; }
Detail::SmoothDirectXYResult Detail::constructDirectSmoothArcXY(const VofSmoothArc& arc,const GeometricVofPlicCell& cell,
    long double tolerance,SmoothDirectXYWork* output)
{
    SmoothDirectXYWork local;auto& work=output?*output:local;work={};
    try{return directSmoothXY(arc,cell,tolerance,work);}
    catch(const VofSmoothGeometryError&){throw;}
    catch(const std::exception& error){throw VofSmoothGeometryError(Failure::InfeasibleVolume,arc.nativeCell,error.what());}
}

VofSmoothGeometryError::VofSmoothGeometryError(VofSmoothGeometryFailure reason,int nativeCell,const char* message)
    :std::runtime_error("[Interface/VofSmoothGeometry cell="+std::to_string(nativeCell)+"] "+message),
    reason_(reason),nativeCell_(nativeCell) {}
VofSmoothArc::Point VofSmoothArc::sample(long double t) const
{
    require(std::isfinite(t)&&t>=0&&t<=1,Failure::InvalidCurve,nativeCell,"curve parameter outside [0,1]");
    return bezier(controlPoints,t);
}
VofSmoothArc::Point VofSmoothArc::derivative(long double t) const
{
    require(std::isfinite(t)&&t>=0&&t<=1,Failure::InvalidCurve,nativeCell,"curve parameter outside [0,1]");
    std::array<Point,4> d{};for(unsigned i=0;i<4;++i)d[i]=mul(sub(controlPoints[i+1],controlPoints[i]),4);
    return bezier(d,t);
}
VofSmoothArc::Point VofSmoothArc::normal(long double t) const
{ const auto d=derivative(t);return unit({d.y,-d.x},nativeCell); }
struct VofSmoothGeometry::Data
{
    explicit Data(const GeometricVofPublication& p):source(p){}
    GeometricVofPublication source;
    std::vector<VofSmoothCell> cells;
    std::vector<VofSmoothNode> nodes;
    std::vector<VofSmoothArc> arcs;
    VofSmoothGeometryReport report;
    long double stripTolerance=1e-12L;
};
VofSmoothGeometry::VofSmoothGeometry(std::shared_ptr<const Data> data):data_(std::move(data)){}
VofSmoothGeometry VofSmoothGeometry::prepare(const QuadtreeMesh& mesh,
    const GeometricVofPublication& source,const VofSmoothGeometryConfig& config)
{
    source.validateMesh(mesh);
    require(source.geometryKind()==GeometricVofGeometryKind::LocalPlic,Failure::SourceKind,-1,"requires actual LocalPlic source");
    require(source.domain().uniformLeaves(),Failure::UnsupportedMesh,-1,"first B1 geometry requires same-level native leaves");
    require(std::isfinite(config.relativeVolumeTolerance)&&config.relativeVolumeTolerance>0
        &&config.relativeVolumeTolerance<=1e-12L,Failure::VolumeResidual,-1,"volume tolerance must be in (0,1e-12]");
    auto data=std::make_shared<Data>(source);
    data->stripTolerance=config.relativeVolumeTolerance;
    const auto& plic=source.plicCells();const auto n=plic.size();
    data->cells.resize(n);
    std::vector<std::vector<int>> cellFaces(n),faceEnds(static_cast<std::size_t>(source.domain().nativeFaceCount()));
    CellSides neighbors(n),faceRows(n);for(auto& sides:neighbors)sides.fill(-1);for(auto& sides:faceRows)sides.fill(-1);
    for(const auto f:source.domain().faces(mesh)){
        const unsigned axis=f.horizontal?1:0;
        for(unsigned side=0;side<2;++side){const int ci=side?f.positiveCell:f.negativeCell;if(ci<0)continue;
            cellFaces[static_cast<std::size_t>(ci)].push_back(f.nativeRow);
            const auto slot=2*axis+(side?0:1);
            require(faceRows[static_cast<std::size_t>(ci)][slot]<0,Failure::UnsupportedMesh,ci,"same-level seed has split native sides");
            faceRows[static_cast<std::size_t>(ci)][slot]=f.nativeRow;
            neighbors[static_cast<std::size_t>(ci)][slot]=side?f.negativeCell:f.positiveCell;
        }
    }
    // B1 owns its reconstruction seed. The source's alpha node-ring LS can
    // acquire a transverse boundary bias even when every tangential row has
    // identical stocks. A fixed direct-axis inventory difference preserves
    // that invariance, without changing the source VOF plane or any stock.
    auto seed=plic;
    const bool rz=source.domain().axisymmetric();
    const auto fraction=[&](std::size_t ci){return plic[ci].inventory[0]/Cut::capacity(plic[ci].rectangle,rz);};
    for(std::size_t ci=0;ci<n;++ci){
        const auto& in=plic[ci];if(in.inventory[0]==0||in.inventory[1]==0)continue;
        if(!in.resolved){
            std::ostringstream detail;detail<<std::setprecision(std::numeric_limits<long double>::max_digits10)
                <<"unresolved PLIC cannot supply B1 smooth geometry";
            appendSourceReceipt(detail,source);appendCellReceipt(detail,plic,static_cast<int>(ci),rz);
            for(const auto row:faceRows[ci]){const auto f=source.domain().face(mesh,row);
                detail<<" adjacent_face="<<row<<" cells="<<f.negativeCell<<','<<f.positiveCell;
                const auto other=f.negativeCell==static_cast<int>(ci)?f.positiveCell:f.negativeCell;
                if(other>=0)appendCellReceipt(detail,plic,other,rz);}
            throw VofSmoothGeometryError(Failure::UnresolvedCell,static_cast<int>(ci),detail.str().c_str());
        }
        const auto gradient=[&](unsigned axis){
            int lower=-1,upper=-1;
            for(const auto row:cellFaces[ci]){const auto f=source.domain().face(mesh,row);
                if(static_cast<unsigned>(f.horizontal)!=axis)continue;
                if(f.positiveCell==static_cast<int>(ci))lower=f.negativeCell;
                if(f.negativeCell==static_cast<int>(ci))upper=f.positiveCell;}
            if(lower<0&&upper<0)return 0.L;
            const auto a=lower>=0?static_cast<std::size_t>(lower):ci;
            const auto b=upper>=0?static_cast<std::size_t>(upper):ci;
            const auto distance=coordinate(midpoint(plic[b].rectangle),axis)-coordinate(midpoint(plic[a].rectangle),axis);
            require(distance>0&&std::isfinite(distance),Failure::UnsupportedMesh,static_cast<int>(ci),"invalid direct-axis native seed stencil");
            return (fraction(b)-fraction(a))/distance;
        };
        const auto gx=gradient(0),gy=gradient(1);
        const auto normal=unit({-gx,-gy},static_cast<int>(ci));
        const bool liquid=in.inventory[0]<=in.inventory[1];
        const Vector2D selected{static_cast<double>(liquid?normal.x:-normal.x),static_cast<double>(liquid?normal.y:-normal.y)};
        const auto reconstruction=Cut::reconstruct(in.rectangle,selected,
            in.inventory[liquid?0:1]/Cut::capacity(in.rectangle,rz),rz);
        seed[ci].normal=liquid?reconstruction.normal:Vector2D{-reconstruction.normal.x,-reconstruction.normal.y};
        seed[ci].offset=liquid?reconstruction.offset:-reconstruction.offset;
    }
    // Single-arc inventory topology qualification, before chart selection.
    // A positive phase cannot enter a neighbor whose authoritative stock is
    // exactly zero. Every qualified mixed cell requires two actual exits.
    // Only forced implications are propagated; unresolved choices are not
    // optimized or assigned an invented connectivity.
    std::vector<bool> mixed(n),queued(n);std::vector<int> queue;
    for(std::size_t ci=0;ci<n;++ci){mixed[ci]=plic[ci].inventory[0]>0&&plic[ci].inventory[1]>0;
        if(mixed[ci]){queue.push_back(static_cast<int>(ci));queued[ci]=true;}}
    std::vector<int> faceStatus(static_cast<std::size_t>(source.domain().nativeFaceCount()),-1);
    // Passive failure provenance only; this array never selects a geometry.
    std::vector<int> faceDecisionCell(faceStatus.size(),-1);
    for(const auto f:source.domain().faces(mesh)){
        const bool negative=f.negativeCell<0||mixed[static_cast<std::size_t>(f.negativeCell)];
        const bool positive=f.positiveCell<0||mixed[static_cast<std::size_t>(f.positiveCell)];
        if(negative&&positive)faceStatus[static_cast<std::size_t>(f.nativeRow)]=0;
        else faceDecisionCell[static_cast<std::size_t>(f.nativeRow)]=negative?f.positiveCell:f.negativeCell;
    }
    for(std::size_t cursor=0;cursor<queue.size();++cursor){
        const int ci=queue[cursor];queued[static_cast<std::size_t>(ci)]=false;
        unsigned forced=0,possible=0;
        for(const auto row:faceRows[static_cast<std::size_t>(ci)]){
            const auto status=faceStatus[static_cast<std::size_t>(row)];if(status>=0)++possible;if(status==1)++forced;}
        if(forced>2||possible<2){
            const auto& cell=seed[static_cast<std::size_t>(ci)];
            const auto segment=Cut::interfaceSegment(cell.rectangle,cell.normal,cell.offset);
            if(segment){const auto hx=(static_cast<long double>(cell.rectangle.upper.x)-cell.rectangle.lower.x)/2;
                const auto hy=(static_cast<long double>(cell.rectangle.upper.y)-cell.rectangle.lower.y)/2;
                for(const auto p:*segment)require(!(std::abs(p.x)==hx&&std::abs(p.y)==hy),
                    Failure::CornerEndpoint,ci,"single-arc inventory topology meets a native corner");}
            std::ostringstream detail;detail<<std::setprecision(std::numeric_limits<long double>::max_digits10)
                <<"inventory-qualified native exits cannot form one cell arc; forced="<<forced<<" possible="<<possible;
            appendSourceReceipt(detail,source);appendCellReceipt(detail,plic,ci,rz);
            for(const auto row:faceRows[static_cast<std::size_t>(ci)]){
                const auto index=static_cast<std::size_t>(row);const auto f=source.domain().face(mesh,row);
                const auto decision=faceDecisionCell[index];
                detail<<" face_receipt=(row="<<row<<",status="<<faceStatus[index]<<",cells="
                    <<f.negativeCell<<','<<f.positiveCell<<",decision="<<decision<<",rule=";
                if(decision<0)detail<<"undecided";
                else if(!mixed[static_cast<std::size_t>(decision)])detail<<"initial_pure_neighbor";
                else detail<<(faceStatus[index]==1?"possible_two":"forced_two");
                detail<<')';
                const auto other=f.negativeCell==ci?f.positiveCell:f.negativeCell;
                if(other>=0)appendCellReceipt(detail,plic,other,rz);
                if(decision>=0&&decision!=ci&&decision!=other)appendCellReceipt(detail,plic,decision,rz);
            }
            throw VofSmoothGeometryError(Failure::UnmatchedEndpoint,ci,detail.str().c_str());
        }
        if(possible!=2&&forced!=2)continue;
        const int assigned=possible==2?1:-1;
        for(const auto row:faceRows[static_cast<std::size_t>(ci)])if(faceStatus[static_cast<std::size_t>(row)]==0){
            faceStatus[static_cast<std::size_t>(row)]=assigned;
            faceDecisionCell[static_cast<std::size_t>(row)]=ci;
            const auto f=source.domain().face(mesh,row);
            for(const auto other:{f.negativeCell,f.positiveCell})if(other>=0&&other!=ci&&mixed[static_cast<std::size_t>(other)]&&!queued[static_cast<std::size_t>(other)]){
                queue.push_back(other);queued[static_cast<std::size_t>(other)]=true;}
        }
    }
    std::vector<Endpoint> endpoints;std::vector<std::array<std::size_t,2>> cellEnds(n);
    std::vector<std::optional<HeightSeed>> heights(n);
    struct EndFace {int row=-1,owner=-1;unsigned phase=0;};std::vector<EndFace> qualifiedEndFaces;
    for(std::size_t ci=0;ci<n;++ci){
        auto& out=data->cells[ci];const auto& in=seed[ci];out.nativeCell=static_cast<int>(ci);
        out.rectangle=in.rectangle;out.inventory=in.inventory;
        require(in.resolved,Failure::UnresolvedCell,out.nativeCell,"unresolved PLIC cannot supply B1 smooth geometry");
        if(in.inventory[0]==0||in.inventory[1]==0){pureQuadrature(out,rz,in.inventory[0]==0?1:0);continue;}
        const std::array<std::optional<HeightSeed>,2> qualified{{
            heightSeed(mesh,source.domain(),out.nativeCell,0,seed,neighbors,faceRows,faceStatus),
            heightSeed(mesh,source.domain(),out.nativeCell,1,seed,neighbors,faceRows,faceStatus)}};
        const unsigned preferred=std::abs(in.normal.x)>=std::abs(in.normal.y)?0:1;
        if(!qualified[0]&&!qualified[1]){
            const auto segment=Cut::interfaceSegment(in.rectangle,in.normal,in.offset);
            if(segment){const auto hx=(static_cast<long double>(in.rectangle.upper.x)-in.rectangle.lower.x)/2;
                const auto hy=(static_cast<long double>(in.rectangle.upper.y)-in.rectangle.lower.y)/2;
                for(const auto p:*segment)require(!(std::abs(p.x)==hx&&std::abs(p.y)==hy),
                    Failure::CornerEndpoint,out.nativeCell,"unqualified height seed has a native-corner branch ambiguity");}
            std::ostringstream detail;detail<<std::setprecision(17)<<"no qualified finite fixed-order height seed; stock="
                <<in.inventory[0]<<','<<in.inventory[1]<<" direct_n="<<in.normal.x<<','<<in.normal.y;
            throw VofSmoothGeometryError(Failure::UnresolvedCell,out.nativeCell,detail.str().c_str());
        }
        const auto selected=qualified[preferred]?preferred:1-preferred;
        heights[ci]=qualified[selected];const auto& height=*heights[ci];
        const auto segment=heightEndpoints(height,in.rectangle,out.nativeCell,faceRows[ci]);
        for(const auto& item:height.endFaces)qualifiedEndFaces.push_back({item.first,out.nativeCell,item.second});
        for(unsigned end=0;end<2;++end){
            const auto p=segment[end].point;const auto face=segment[end].face;
            require(face>=0,Failure::UnmatchedEndpoint,out.nativeCell,"endpoint has no actual native face");
            cellEnds[ci][end]=endpoints.size();faceEnds[static_cast<std::size_t>(face)].push_back(static_cast<int>(endpoints.size()));
            endpoints.push_back({out.nativeCell,face,p,segment[end].normal});
        }
    }
    std::vector<std::size_t> endpointNode(endpoints.size());
    for(std::size_t row=0;row<faceEnds.size();++row){
        const auto& list=faceEnds[row];if(list.empty())continue;
        const auto f=source.domain().face(mesh,static_cast<int>(row));
        const bool exterior=f.negativeCell<0||f.positiveCell<0;
        if(list.size()!=(exterior?1U:2U)){
            std::ostringstream detail;detail<<std::setprecision(17)
                <<"native face does not have the qualified endpoint incidence; face="<<row<<" count="<<list.size()
                <<" horizontal="<<f.horizontal<<" center="<<f.center.x<<','<<f.center.y;
            for(const auto native:{f.negativeCell,f.positiveCell})if(native>=0){
                const auto& s=seed[static_cast<std::size_t>(native)];const auto& old=plic[static_cast<std::size_t>(native)];
                detail<<" neighbor="<<native<<" stock="<<s.inventory[0]<<','<<s.inventory[1]
                    <<" seed_n="<<s.normal.x<<','<<s.normal.y<<" offset="<<s.offset
                    <<" source_n="<<old.normal.x<<','<<old.normal.y;
                if(heights[static_cast<std::size_t>(native)]){const auto& h=*heights[static_cast<std::size_t>(native)];
                    detail<<" height_axis="<<h.normalAxis<<" stencil="<<h.tangentOffset<<','<<h.normalLower<<':'<<h.normalUpper
                        <<" strips="<<h.strips<<" polynomial=";
                    for(const auto coefficient:h.polynomial)detail<<coefficient<<',';}
                const auto segment=Cut::interfaceSegment(s.rectangle,s.normal,s.offset);
                if(segment)for(const auto p:*segment)detail<<" local_end="<<p.x<<','<<p.y;
            }
            throw VofSmoothGeometryError(Failure::UnmatchedEndpoint,endpoints[static_cast<std::size_t>(list[0])].cell,detail.str().c_str());
        }
        require(!exterior||config.allowBoundaryEndpoints,Failure::BoundaryEndpoint,endpoints[static_cast<std::size_t>(list[0])].cell,
            "boundary endpoint requires explicit frozen-open-curve scope");
        Point p{},normal{};
        for(const auto i:list){p=add(p,mul(endpoints[static_cast<std::size_t>(i)].point,1.L/list.size()));
            normal=add(normal,endpoints[static_cast<std::size_t>(i)].normal);}
        // The native face owns its normal coordinate. Only the two actual
        // Height-seed tangential intersections are averaged for the new node.
        if(f.horizontal)p.y=f.center.y;else p.x=f.center.x;
        if(list.size()==2)
            require(dot(endpoints[static_cast<std::size_t>(list[0])].normal,
                endpoints[static_cast<std::size_t>(list[1])].normal)>0,
                Failure::InvalidTangent,endpoints[static_cast<std::size_t>(list[0])].cell,
                "opposed height-seed phase normals do not define a qualified shared tangent");
        normal=unit(normal,-1);const auto tangent=Point{-normal.y,normal.x};
        const auto id=data->nodes.size();data->nodes.push_back({static_cast<int>(row),p,tangent});
        for(const auto i:list)endpointNode[static_cast<std::size_t>(i)]=id;
        if(exterior)++data->report.boundaryEndpoints;
    }
    for(std::size_t ci=0;ci<n;++ci){
        auto& cell=data->cells[ci];const auto& in=seed[ci];
        if(in.inventory[0]==0||in.inventory[1]==0)continue;
        VofSmoothArc base;base.nativeCell=static_cast<int>(ci);
        base.nodes={endpointNode[cellEnds[ci][0]],endpointNode[cellEnds[ci][1]]};
        auto p0=data->nodes[base.nodes[0]].point,p1=data->nodes[base.nodes[1]].point;
        if(dot(sub(p1,p0),data->nodes[base.nodes[0]].tangent)<0){
            std::swap(base.nodes[0],base.nodes[1]);std::swap(p0,p1);}
        const auto chord=sub(p1,p0);
        const auto length=norm(chord);
        require(length>0&&std::isfinite(length),Failure::InvalidTangent,base.nativeCell,"zero shared-endpoint chord");
        const auto d0=mul(data->nodes[base.nodes[0]].tangent,length),d1=mul(data->nodes[base.nodes[1]].tangent,length);
        require(dot(chord,d0)>0&&dot(chord,d1)>0,Failure::InvalidTangent,base.nativeCell,"shared tangents oppose cell traversal");
        base.graphAxis=std::abs(chord.x)>=std::abs(chord.y)?0:1;
        base.controlPoints={p0,add(p0,mul(d0,.25L)),
            add(mul(add(p0,p1),.5L),mul(sub(d0,d1),1.L/6)),sub(p1,mul(d1,.25L)),p1};
        if(!rz){
            Detail::SmoothDirectXYWork work;
            try {
                auto direct=Detail::constructDirectSmoothArcXY(base,in,config.relativeVolumeTolerance,&work);
                data->report.xyDirectCandidates+=work.candidateConstructions;
                data->report.xyGeometryCertifications+=work.geometryCertifications;
                data->report.xyPhaseQuadratureBuilds+=work.phaseQuadratureBuilds;
                data->report.xyFitPhaseMeasurements+=work.fitPhaseMeasurements;
                data->report.xyVolumeNewtonSteps+=work.volumeNewtonSteps;
                data->report.maximumCertificateSubdivisionDepth=std::max(data->report.maximumCertificateSubdivisionDepth,
                    direct.maximumCertificateSubdivisionDepth);
                cell.arc=static_cast<int>(data->arcs.size());data->arcs.push_back(std::move(direct.arc));
                cell.measured=direct.measured;cell.phaseQuadrature=std::move(direct.phaseQuadrature);
                ++data->report.mixedCells;
            }catch(const VofSmoothGeometryError& error){
                std::ostringstream detail;detail<<std::setprecision(17)<<error.what()
                    <<" direct_xy=(candidates="<<work.candidateConstructions<<",fit_measures="<<work.fitPhaseMeasurements
                    <<",newton="<<work.volumeNewtonSteps<<",certificates="<<work.geometryCertifications
                    <<",quadratures="<<work.phaseQuadratureBuilds<<",domain_leaves="<<work.domainLeaves<<") base_controls=";
                for(const auto point:base.controlPoints)detail<<'('<<point.x<<','<<point.y<<')';
                appendSourceReceipt(detail,source);appendCellReceipt(detail,plic,cell.nativeCell,false);
                throw VofSmoothGeometryError(error.reason(),cell.nativeCell,detail.str().c_str());
            }
            continue;
        }
        // RZ retains its original transverse polynomial/three-update policy.
        // Fixed graph-transverse height direction, NOT a surface normal.
        // The graph coordinate is independent of delta; actual surface normals
        // always come from the retained curve derivative.
        const auto correctionDirection=base.graphAxis==0?Point{0,chord.x>0?-1.L:1.L}:
            Point{chord.y>0?1.L:-1.L,0};
        const auto scale=std::min(static_cast<long double>(cell.rectangle.upper.x)-cell.rectangle.lower.x,
            static_cast<long double>(cell.rectangle.upper.y)-cell.rectangle.lower.y);
        const auto curve=[&](long double lambda){auto a=base;a.volumeCorrection=lambda*scale;
            a.controlPoints[2]=add(a.controlPoints[2],mul(correctionDirection,a.volumeCorrection/6));return a;};
        // Match the minority phase without subtracting it from capacity.
        const unsigned phase=in.inventory[0]<in.inventory[1]?0:1;
        const auto target=in.inventory[phase];
        const auto v0=Detail::measureSmoothArcPhases(base,in,rz)[phase];
        const auto vp=Detail::measureSmoothArcPhases(curve(1),in,rz)[phase],vm=Detail::measureSmoothArcPhases(curve(-1),in,rz)[phase];
        const auto c=v0-target,b=(vp-vm)/2;
        // Only a radial transverse displacement changes the RZ quadratic
        // coefficient. This is the analytic degree, never a small-value cut.
        const auto a=rz&&base.graphAxis==1?(vp+vm)/2-v0:0.L;
        std::vector<long double> roots;
        if(a==0){if(b!=0)roots.push_back(-c/b);else if(c==0)roots.push_back(0);}
        else {
            const auto discriminant=b*b-4*a*c;
            if(discriminant>=0&&std::isfinite(discriminant)){
                const auto q=-.5L*(b+std::copysign(std::sqrt(discriminant),b));
                if(q!=0){roots.push_back(q/a);roots.push_back(c/q);}else roots.push_back(-b/(2*a));
            }
        }
        std::sort(roots.begin(),roots.end(),[](long double x,long double y){return std::abs(x)<std::abs(y);});
        bool found=false,lastCertified=false;VofSmoothArc accepted,lastCandidate=base;
        auto lastRoot=std::numeric_limits<long double>::quiet_NaN();
        std::array<long double,2> lastMeasure{};
        unsigned certificateDepth=0;
        for(auto root:roots){
            if(!std::isfinite(root))continue;
            std::array<VofSmoothArc,4> actualNative;unsigned actualCount=0;
            // Refine the one physical polynomial volume root, preserving its
            // mathematical law and all geometric admissibility conditions.
            for(unsigned iteration=0;iteration<3;++iteration){
                actualNative[actualCount++]=curve(root);
                const auto residual=Detail::measureSmoothArcPhases(actualNative[actualCount-1],in,rz)[phase]-target;
                const auto slope=b+2*a*root;if(slope==0)break;
                const auto next=root-residual/slope;if(!std::isfinite(next))break;root=next;
            }
            auto candidate=curve(root);lastCandidate=candidate;lastRoot=root;
            actualNative[actualCount++]=candidate;
            auto measure=Detail::measureSmoothArcPhases(candidate,in,rz);lastMeasure=measure;
            unsigned depth=0;lastCertified=certified(candidate,cell.rectangle,depth);
            certificateDepth=std::max(certificateDepth,depth);
            if(!lastCertified)continue;
            bool accurate=true;
            for(unsigned p=0;p<2;++p)
                accurate=accurate&&std::abs(measure[p]-in.inventory[p])<=config.relativeVolumeTolerance*in.inventory[p];
            if(!accurate)continue;
            accepted=std::move(candidate);found=true;break;
        }
        if(!found){
            std::ostringstream detail;detail<<std::setprecision(17)
                <<"no certified Bernstein-bounded monotone curve matches both phase stocks; roots="<<roots.size()
                <<" last_lambda="<<lastRoot<<" delta="<<lastCandidate.volumeCorrection
                <<" certified="<<lastCertified<<" subdivision="<<certificateDepth
                <<" stock="<<in.inventory[0]<<','<<in.inventory[1]
                <<" signed_measure="<<lastMeasure[0]<<','<<lastMeasure[1]
                <<" box="<<cell.rectangle.lower.x<<','<<cell.rectangle.lower.y<<':'
                <<cell.rectangle.upper.x<<','<<cell.rectangle.upper.y<<" controls=";
            for(const auto p:lastCandidate.controlPoints)detail<<'('<<p.x<<','<<p.y<<')';
            for(const auto node:base.nodes){const auto row=data->nodes[node].nativeRow;
                detail<<" raw_node_face="<<row;
                for(const auto index:faceEnds[static_cast<std::size_t>(row)]){
                    const auto& raw=endpoints[static_cast<std::size_t>(index)];
                    detail<<" raw_cell="<<raw.cell<<" point="<<raw.point.x<<','<<raw.point.y
                        <<" normal="<<raw.normal.x<<','<<raw.normal.y;
                    if(heights[static_cast<std::size_t>(raw.cell)]){const auto& h=*heights[static_cast<std::size_t>(raw.cell)];
                        detail<<" hf_axis="<<h.normalAxis<<" liquidLower="<<h.liquidLower
                            <<" stencil="<<h.tangentOffset<<','<<h.normalLower<<':'<<h.normalUpper
                            <<" origin="<<h.normalOrigin<<','<<h.tangentOrigin<<" scale="<<h.tangentScale
                            <<" strips="<<h.strips<<" poly=";
                        for(const auto coefficient:h.polynomial)detail<<coefficient<<',';}
                }
            }
            Point lower=lastCandidate.sample(0),upper=lower;
            long double minimumGraphSpeed=std::numeric_limits<long double>::infinity();
            const auto direction=coordinate(lastCandidate.controlPoints.back(),base.graphAxis)>
                coordinate(lastCandidate.controlPoints.front(),base.graphAxis)?1.L:-1.L;
            for(unsigned i=0;i<=64;++i){const auto t=i/64.L;const auto p=lastCandidate.sample(t);
                lower.x=std::min(lower.x,p.x);lower.y=std::min(lower.y,p.y);
                upper.x=std::max(upper.x,p.x);upper.y=std::max(upper.y,p.y);
                minimumGraphSpeed=std::min(minimumGraphSpeed,direction*coordinate(lastCandidate.derivative(t),base.graphAxis));}
            detail<<" diagnostic_samples="<<lower.x<<','<<lower.y<<':'<<upper.x<<','<<upper.y
                <<" diagnostic_min_graph_speed="<<minimumGraphSpeed;
            appendSourceReceipt(detail,source);appendCellReceipt(detail,plic,cell.nativeCell,rz);
            throw VofSmoothGeometryError(Failure::InfeasibleVolume,cell.nativeCell,detail.str().c_str());
        }
        data->report.maximumCertificateSubdivisionDepth=std::max(data->report.maximumCertificateSubdivisionDepth,certificateDepth);
        powerProjection(accepted);cell.arc=static_cast<int>(data->arcs.size());data->arcs.push_back(accepted);
        cell.measured=Detail::measureSmoothArcPhases(data->arcs.back(),in,rz,&cell.phaseQuadrature);
        ++data->report.mixedCells;
    }
    for(const auto& cell:data->cells)for(unsigned phase=0;phase<2;++phase){
        const auto stock=cell.inventory[phase],measured=cell.measured[phase];
        if(stock==0)require(measured==0,Failure::VolumeResidual,cell.nativeCell,"nonzero geometry for exactly empty phase");
        else {
            const auto relative=std::abs(measured-stock)/stock;
            require(relative<=config.relativeVolumeTolerance,Failure::VolumeResidual,cell.nativeCell,"own-phase quadrature volume residual");
            data->report.maximumOwnPhaseRelativeResidual=std::max(data->report.maximumOwnPhaseRelativeResidual,relative);
        }
    }
    data->report.nodes=data->nodes.size();VofSmoothGeometry result(data);
    for(const auto& qualified:qualifiedEndFaces){
        require(faceEnds[static_cast<std::size_t>(qualified.row)].empty(),Failure::UnmatchedEndpoint,qualified.owner,
            "qualified height window end face is crossed by the actual final geometry");
        const auto f=source.domain().face(mesh,qualified.row);
        for(const auto native:{f.negativeCell,f.positiveCell})if(native>=0)
            require(result.phaseAt(native,{f.center.x,f.center.y})==qualified.phase,
                Failure::UnmatchedEndpoint,qualified.owner,"actual final phase disagrees with its height-window end-face seed qualification");
    }
    return result;
}
const GeometricVofPublication& VofSmoothGeometry::source() const { return data_->source; }
GeometricVofIdentity VofSmoothGeometry::sourceIdentity() const { return source().identity(); }
GeometricVofTimeLayer VofSmoothGeometry::sourceTime() const { return source().timeLayer(); }
const std::vector<VofSmoothCell>& VofSmoothGeometry::cells() const { return data_->cells; }
const std::vector<VofSmoothNode>& VofSmoothGeometry::nodes() const { return data_->nodes; }
const std::vector<VofSmoothArc>& VofSmoothGeometry::arcs() const { return data_->arcs; }
const VofSmoothGeometryReport& VofSmoothGeometry::report() const { return data_->report; }
unsigned VofSmoothGeometry::phaseAt(int nativeCell,Point point) const
{
    require(nativeCell>=0&&static_cast<std::size_t>(nativeCell)<data_->cells.size(),Failure::ForeignSource,nativeCell,
        "phase query has no actual native cell owner");
    const auto& cell=data_->cells[static_cast<std::size_t>(nativeCell)];const auto& box=cell.rectangle;
    require(std::isfinite(point.x)&&std::isfinite(point.y)&&point.x>=box.lower.x&&point.x<=box.upper.x
        &&point.y>=box.lower.y&&point.y<=box.upper.y,Failure::InvalidCurve,nativeCell,"phase query outside its actual native rectangle");
    if(cell.inventory[0]==0)return 1;
    if(cell.inventory[1]==0)return 0;
    require(cell.arc>=0&&static_cast<std::size_t>(cell.arc)<data_->arcs.size(),Failure::InvalidCurve,nativeCell,
        "mixed phase query has no actual smooth arc");
    const auto& arc=data_->arcs[static_cast<std::size_t>(cell.arc)];const auto axis=arc.graphAxis;
    const auto a=coordinate(arc.controlPoints.front(),axis),b=coordinate(arc.controlPoints.back(),axis);
    const auto u=coordinate(point,axis),v=coordinate(point,1-axis);
    const bool increasing=b>a,liquidLower=axis?increasing:!increasing;
    if(u<std::min(a,b)||u>std::max(a,b))return graphOuterPhase(arc,box,u,liquidLower);
    long double lo=0,hi=1;
    for(unsigned iteration=0;iteration<80;++iteration){
        const auto mid=(lo+hi)/2;if(mid==lo||mid==hi)break;
        const auto value=coordinate(arc.sample(mid),axis);
        if(value==u){lo=hi=mid;break;}
        if((value<u)==increasing)lo=mid;else hi=mid;
    }
    const auto t=u==a?0.L:u==b?1.L:(lo+hi)/2;
    const auto curveV=coordinate(arc.sample(t),1-axis);
    // This is a query qualification band for roundoff at the retained trace,
    // not a displaced interface, inventory threshold, or phase repair.
    const auto width=axis?static_cast<long double>(box.upper.x)-box.lower.x:
        static_cast<long double>(box.upper.y)-box.lower.y;
    const auto roundoff=64*std::numeric_limits<long double>::epsilon()*
        std::max({std::abs(v),std::abs(curveV),width});
    require(std::abs(v-curveV)>roundoff,Failure::InvalidCurve,nativeCell,"phase query is indistinguishable from the actual smooth trace");
    return (v<curveV)==liquidLower?0:1;
}
VofSmoothStripMeasure VofSmoothGeometry::measureNativeStrip(const QuadtreeMesh& mesh,
    const GeometricVofPublication& actualSource,int nativeCell,unsigned normalAxis,
    bool upperFace,long double actualSweptVolume) const
{
    validateFor(mesh,actualSource);
    require(source().domain().uniformLeaves()&&!source().domain().axisymmetric(),Failure::UnsupportedMesh,nativeCell,
        "smooth native strips require a same-level XY source");
    require(nativeCell>=0&&static_cast<std::size_t>(nativeCell)<data_->cells.size(),Failure::ForeignSource,nativeCell,
        "smooth native strip has no actual cell owner");
    const auto& cell=data_->cells[static_cast<std::size_t>(nativeCell)];const auto& box=cell.rectangle;
    const auto capacity=Cut::capacity(box,false);
    require(normalAxis<2&&std::isfinite(actualSweptVolume)&&actualSweptVolume>=0&&actualSweptVolume<=capacity,
        Failure::InvalidCurve,nativeCell,"invalid smooth native strip axis or requested volume");
    VofSmoothStripMeasure result;result.requestedVolume=actualSweptVolume;
    if(actualSweptVolume==0)return result;
    if(cell.inventory[0]==0||cell.inventory[1]==0){
        result.phase[cell.inventory[0]==0?1:0]=actualSweptVolume;result.rawPhase=result.phase;return result;}
    require(cell.arc>=0&&static_cast<std::size_t>(cell.arc)<data_->arcs.size(),Failure::InvalidCurve,nativeCell,
        "mixed smooth native strip has no retained arc");
    const auto tangentWidth=normalAxis?static_cast<long double>(box.upper.x)-box.lower.x:
        static_cast<long double>(box.upper.y)-box.lower.y;
    const auto normalWidth=normalAxis?static_cast<long double>(box.upper.y)-box.lower.y:
        static_cast<long double>(box.upper.x)-box.lower.x;
    require(std::isfinite(tangentWidth)&&tangentWidth>0&&std::isfinite(normalWidth)&&normalWidth>0,
        Failure::VolumeResidual,nativeCell,"invalid face-anchored smooth strip widths");
    try{
        const auto sums=smoothStrip(data_->arcs[static_cast<std::size_t>(cell.arc)],box,normalAxis,upperFace,
            tangentWidth,actualSweptVolume,data_->stripTolerance,result.polynomialRoots);
        using Exact=Geometry::Detail::FiniteVolume::Dyadic;
        Exact combined;
        for(unsigned p=0;p<2;++p){result.rawPhase[p]=sums.phase[p].rounded();combined.add(sums.phase[p]);
            require(std::isfinite(result.rawPhase[p])&&result.rawPhase[p]>=0,Failure::VolumeResidual,nativeCell,
                "invalid independent positive smooth strip phase");
            auto error=sums.error[p];error.add(Exact::product({std::numeric_limits<long double>::epsilon(),result.rawPhase[p]}));
            result.absoluteError[p]=error.rounded();
            require(result.absoluteError[p]>0||error.compare(Exact{})==0,Failure::VolumeResidual,nativeCell,
                "unrepresentable smooth strip own-phase error bound");
            require(result.rawPhase[p]==0?result.absoluteError[p]==0:
                result.absoluteError[p]<=data_->stripTolerance*result.rawPhase[p],Failure::VolumeResidual,nativeCell,
                "smooth strip root/integral uncertainty exceeds its own-phase budget");}
        auto defect=combined;defect.add(actualSweptVolume,-1);
        result.relativeCoverageResidual=defect.absolute().rounded()/capacity;
        require(result.relativeCoverageResidual<=data_->stripTolerance,Failure::VolumeResidual,nativeCell,
            "independent smooth strip phases do not cover the requested native volume");
        const auto sum=combined.rounded();require(std::isfinite(sum)&&sum>0,Failure::VolumeResidual,nativeCell,
            "positive smooth strip has no phase measure");
        const unsigned minor=result.rawPhase[0]<=result.rawPhase[1]?0:1;
        result.phase[minor]=actualSweptVolume*(result.rawPhase[minor]/sum);
        result.phase[1-minor]=actualSweptVolume-result.phase[minor];
        for(unsigned p=0;p<2;++p)require(std::isfinite(result.phase[p])&&result.phase[p]>=0&&
            (result.rawPhase[p]==0||result.phase[p]>0),Failure::VolumeResidual,nativeCell,
            "shared smooth strip normalization lost a positive phase");
        return result;
    }catch(const VofSmoothGeometryError&){throw;}
    catch(const std::exception& error){throw VofSmoothGeometryError(Failure::VolumeResidual,nativeCell,error.what());}
}
bool VofSmoothGeometry::sameInstance(const VofSmoothGeometry& other) const noexcept { return data_==other.data_; }
void VofSmoothGeometry::validateFor(const QuadtreeMesh& mesh,const GeometricVofPublication& source) const
{
    require(source.identity()==sourceIdentity(),Failure::ForeignSource,-1,"foreign VOF publication identity");
    source.validateFor(mesh,sourceIdentity(),sourceTime());
}
}
