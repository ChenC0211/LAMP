#include "Core/Geometry/QuadraticMaterialVolume.h"
#include "Core/Geometry/Detail/QuadraticGraphExact.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>

namespace LAMP::Geometry::QuadraticMaterialVolume
{
namespace
{
namespace E=LAMP::Geometry::Detail::QuadraticExact;
using V=std::array<long double,2>;
using Polynomial=std::array<V,3>;
using Gradient=std::array<V,3>;
void require(bool yes,const char* why)
{if(!yes)throw std::invalid_argument(std::string("[QuadraticMaterialVolume] ")+why);}
E::Vector exact(Point p){return E::vector({p.x,p.y});}
struct Sum {long double value=0,error=0;void add(long double x){const auto y=x-error,t=value+y;error=(t-value)-y;value=t;}};
std::array<E::Vector,3> authoritative(const C::Curve& curve)
{
    require(!curve.usesExactBasis()||curve.exactViewUnchanged(),"stale curve coefficient view");
    if(curve.usesExactBasis())return curve.exactData()->coefficients;
    if(const auto* e=curve.endpointBasis()){
        const auto a=exact(e->first),b=exact(e->last),q=exact(e->bubble);
        return {{a,E::add(E::subtract(b,a),q),E::multiply(q,-1)}};
    }
    return {{exact(curve.coefficients[0]),exact(curve.coefficients[1]),exact(curve.coefficients[2])}};
}
struct NativeArc
{
    V anchor{};
    Polynomial source{},displacement{};
    std::array<Point,3> samples{};
    V lower{},upper{};
    int orientation=1;
};
V value(const Polynomial& p,long double s)
{return {{std::fma(s,std::fma(s,p[2][0],p[1][0]),p[0][0]),std::fma(s,std::fma(s,p[2][1],p[1][1]),p[0][1])}};}
V derivative(const Polynomial& p,long double s)
{return {{std::fma(2*s,p[2][0],p[1][0]),std::fma(2*s,p[2][1],p[1][1])}};}
Polynomial at(const NativeArc& a,const Rectangle& box,long double theta)
{
    Polynomial p{};
    for(unsigned j=0;j<3;++j)for(unsigned k=0;k<2;++k)p[j][k]=std::fma(theta,a.displacement[j][k],a.source[j][k]);
    p[0][0]+=a.anchor[0]-box.lower.x;p[0][1]+=a.anchor[1]-box.lower.y;
    return p;
}
void roots(std::vector<long double>& cuts,long double a,long double b,long double c)
{
    const auto scale=std::max({std::abs(a),std::abs(b),std::abs(c)});if(scale==0)return;
    a/=scale;b/=scale;c/=scale;
    const auto put=[&](long double x){if(std::isfinite(x)&&x>0&&x<1)cuts.push_back(x);};
    if(a==0){if(b!=0)put(-c/b);return;}
    long double discriminant=std::fma(b,b,-4*a*c);
    const auto error=16*std::numeric_limits<long double>::epsilon()*(b*b+std::abs(4*a*c));
    // Near a discriminant event, preselect the original rational sign/root
    // arithmetic on these SAME native coefficients. No failed-query fallback.
    if(std::abs(discriminant)<=error){
        E::Polynomial p{};p[0]=E::scalar(c);p[1]=E::scalar(b);p[2]=E::scalar(a);
        for(const auto& root:E::roots(p,E::Scalar(0),E::Scalar(1)))put(E::rounded(root));return;
    }
    if(discriminant<0)return;
    const auto q=-(b+std::copysign(std::sqrt(discriminant),b))/2;
    if(q==0)put(-b/(2*a));else{put(q/a);put(c/q);}
}
std::vector<long double> cuts(const Polynomial& p,long double width,long double height)
{
    std::vector<long double> result{0,1};
    for(unsigned axis=0;axis<2;++axis)for(long double border:{0.L,axis==0?width:height})
        roots(result,p[2][axis],p[1][axis],p[0][axis]-border);
    std::sort(result.begin(),result.end());result.erase(std::unique(result.begin(),result.end()),result.end());return result;
}
long double clippedGreen(const Polynomial& p,int orientation,long double width,long double height)
{
    const auto partition=cuts(p,width,height);Sum sum;
    constexpr long double node=.2886751345948128822545743902509787278L;
    for(std::size_t i=1;i<partition.size();++i){const auto lo=partition[i-1],hi=partition[i],middle=(lo+hi)/2,span=hi-lo;
        const auto center=value(p,middle);if(!(center[1]>0&&center[1]<height)||center[0]<=0)continue;
        for(long double sign:{-1.L,1.L}){const auto s=middle+sign*node*span;
            const auto q=value(p,s),d=derivative(p,s);
            // This saturation is the defining Green test function, not stock
            // clipping. Its x derivative is exactly the rectangle indicator.
            const auto g=center[0]>=width?width:q[0];sum.add(orientation*span*g*d[1]/2);}
    }
    return sum.value;
}
Gradient instantGradient(const Polynomial& p,int orientation,long double width,long double height)
{
    for(unsigned axis=0;axis<2;++axis)if(p[1][axis]==0&&p[2][axis]==0)
        require(p[0][axis]!=0&&p[0][axis]!=(axis==0?width:height),"finite material arc coincides with a native box side");
    Gradient result{};const auto partition=cuts(p,width,height);
    constexpr long double node=.2886751345948128822545743902509787278L;
    for(std::size_t i=1;i<partition.size();++i){const auto lo=partition[i-1],hi=partition[i],middle=(lo+hi)/2,span=hi-lo;
        const auto center=value(p,middle);
        if(!(center[0]>0&&center[0]<width&&center[1]>0&&center[1]<height))continue;
        for(long double sign:{-1.L,1.L}){const auto s=middle+sign*node*span;const auto d=derivative(p,s);
            const std::array<long double,3> n{{std::fma(s,2*s-3,1),4*s*(1-s),s*(2*s-1)}};
            for(unsigned j=0;j<3;++j){result[j][0]+=orientation*span*n[j]*d[1]/2;result[j][1]-=orientation*span*n[j]*d[0]/2;}}
    }
    return result;
}
constexpr std::array<long double,4> nodes4{{-.8611363115940525752L,-.3399810435848562648L,.3399810435848562648L,.8611363115940525752L}};
constexpr std::array<long double,4> weights4{{.3478548451374538574L,.6521451548625461426L,.6521451548625461426L,.3478548451374538574L}};
constexpr std::array<long double,8> nodes8{{-.9602898564975362317L,-.7966664774136267396L,-.5255324099163289858L,-.1834346424956498049L,.1834346424956498049L,.5255324099163289858L,.7966664774136267396L,.9602898564975362317L}};
constexpr std::array<long double,8> weights8{{.1012285362903762592L,.2223810344533744705L,.3137066458778872873L,.3626837833783619830L,.3626837833783619830L,.3137066458778872873L,.2223810344533744705L,.1012285362903762592L}};
struct Bubble {long double area=0,volume=0,error=0;std::size_t samples=0;};
long double bubbleArea(const Polynomial& p,V normal,long double lo,long double hi,bool high,long double integrationLo,long double integrationHi)
{
    const auto width=hi-lo,middle=lo+width/2;
    const auto middleDerivative=derivative(p,middle);
    // cross(d,d') is constant for this SAME quadratic. A compensated 2x2
    // determinant avoids subtracting the tangent projection from d'.
    const auto cross=[](const V& a,const V& b){
        const auto left=a[0]*b[1],right=a[1]*b[0];
        return (left-right)+(std::fma(a[0],b[1],-left)-std::fma(a[1],b[0],-right));};
    const V second{{2*p[2][0],2*p[2][1]}};
    const auto turn=cross(p[1],second);
    const auto eval=[&](long double z){
        const V d{{std::fma((2*p[2][0])*width,z-.5L,middleDerivative[0]),
                   std::fma((2*p[2][1])*width,z-.5L,middleDerivative[1])}};
        const auto length=std::hypot(d[0],d[1]);require(length>0,"nonregular material tangent");
        const auto phi=4*z*(1-z);
        // -t'(s).normal = -cross(c1,2c2)*cross(d,normal)/|d|^3.
        // Integrate in z: no subtraction of nearly equal world s endpoints.
        return -turn*cross(d,normal)/(length*length*length)*phi*width;};
    Sum sum;const auto mid=(integrationLo+integrationHi)/2,half=(integrationHi-integrationLo)/2;
    if(high){for(std::size_t i=0;i<nodes8.size();++i)sum.add(half*weights8[i]*eval(mid+half*nodes8[i]));}
    else for(std::size_t i=0;i<nodes4.size();++i)sum.add(half*weights4[i]*eval(mid+half*nodes4[i]));
    return sum.value;
}
std::array<long double,2> adaptiveBubbleArea(const Polynomial& p,V normal,long double bubbleLo,long double bubbleHi,
    long double lo,long double hi,const Config& cfg,unsigned depth,std::size_t& samples)
{
    samples+=nodes4.size()+nodes8.size();
    const auto a=bubbleArea(p,normal,bubbleLo,bubbleHi,false,lo,hi),b=bubbleArea(p,normal,bubbleLo,bubbleHi,true,lo,hi);
    const auto error=std::abs(a-b),natural=std::max(std::abs(a),std::abs(b));
    if(error<=cfg.areaVolumeRelativeTolerance*natural+128*std::numeric_limits<long double>::epsilon()*natural)return {{b,error}};
    const auto mid=(lo+hi)/2;require(depth<cfg.maximumDepth&&mid>lo&&mid<hi,"physical area-volume quadrature budget exhausted");
    const auto left=adaptiveBubbleArea(p,normal,bubbleLo,bubbleHi,lo,mid,cfg,depth+1,samples),right=adaptiveBubbleArea(p,normal,bubbleLo,bubbleHi,mid,hi,cfg,depth+1,samples);
    return {{left[0]+right[0],left[1]+right[1]}};
}
Bubble bubble(const Polynomial& p,int orientation,long double width,long double height,const Config& cfg)
{
    Bubble out;const auto partition=cuts(p,width,height);
    for(std::size_t i=1;i<partition.size();++i){const auto lo=partition[i-1],hi=partition[i],middle=(lo+hi)/2;
        const auto q=value(p,middle);if(!(q[0]>0&&q[0]<width&&q[1]>0&&q[1]<height))continue;
        const auto d=derivative(p,middle);const auto length=std::hypot(d[0],d[1]);require(length>0,"nonregular bubble tangent");
        const V normal{{orientation*d[1]/length,-orientation*d[0]/length}};
        // Integral phi*d is exact since d is affine and phi symmetric.
        out.volume+=orientation*(2*(hi-lo)/3)*(normal[0]*d[1]-normal[1]*d[0]);
        const auto area=adaptiveBubbleArea(p,normal,lo,hi,0,1,cfg,0,out.samples);out.area+=area[0];out.error+=area[1];
    }
    return out;
}
}
struct Prepared::Data {std::vector<NativeArc> arcs;};
Prepared::Prepared(std::shared_ptr<const Data> data):data_(std::move(data)){}
Prepared Prepared::prepare(const std::vector<ArcMotion>& input)
{
    require(!input.empty(),"requires a nonempty closed material curve");auto data=std::make_shared<Data>();data->arcs.reserve(input.size());
    using Key=std::pair<double,double>;std::map<Key,std::array<unsigned,2>> sourceDegree,imageDegree;
    std::map<Key,E::Vector> authoritativeEndpoint;
    std::map<Key,Point> retainedEndpointDisplacement;
    for(const auto& item:input){require(item.liquidLeftOrientation==1||item.liquidLeftOrientation==-1,"invalid liquid-left orientation");
        const auto c=authoritative(item.source);const auto first=C::position(item.source,0),last=C::position(item.source,1);
        NativeArc a;a.anchor={{first.x,first.y}};a.samples=item.displacement;a.orientation=item.liquidLeftOrientation;
        for(const auto& d:item.displacement)require(std::isfinite(d.x)&&std::isfinite(d.y),"nonfinite retained displacement");
        // Rounded endpoint identity alone cannot certify a closed material
        // path. Both the authoritative source polynomial and its retained
        // displacement must agree at every shared native endpoint identity.
        const auto retainEndpoint=[&](Point p,const E::Vector& physicalPoint,Point displacement){
            const auto [source,sourceInserted]=authoritativeEndpoint.emplace(Key(p.x,p.y),physicalPoint);
            require(sourceInserted||(source->second.x==physicalPoint.x&&source->second.y==physicalPoint.y),
                "authoritative source polynomial endpoints break the closed C0 material path");
            const auto [found,inserted]=retainedEndpointDisplacement.emplace(Key(p.x,p.y),displacement);
            require(inserted||(found->second.x==displacement.x&&found->second.y==displacement.y),
                "retained endpoint displacement breaks the closed C0 material path");
        };
        retainEndpoint(first,c[0],item.displacement[0]);
        retainEndpoint(last,E::add(c[0],E::add(c[1],c[2])),item.displacement[2]);
        auto local=c;local[0]=E::subtract(local[0],exact(first));
        const auto d0=exact(item.displacement[0]),dm=exact(item.displacement[1]),d1=exact(item.displacement[2]);
        const auto db=E::subtract(E::multiply(dm,4),E::multiply(E::add(d0,d1),2));
        const std::array<E::Vector,3> dc{{d0,E::add(E::subtract(d1,d0),db),E::multiply(db,-1)}};
        for(unsigned j=0;j<3;++j){a.source[j]={{E::rounded(local[j].x),E::rounded(local[j].y)}};a.displacement[j]={{E::rounded(dc[j].x),E::rounded(dc[j].y)}};}
        const auto image=C::displacedSamples(item.source,item.displacement[0],item.displacement[1],item.displacement[2]);
        const auto sourceBounds=C::bounds(item.source),imageBounds=C::bounds(image);
        a.lower={{std::min(sourceBounds.lower.x,imageBounds.lower.x),std::min(sourceBounds.lower.y,imageBounds.lower.y)}};
        a.upper={{std::max(sourceBounds.upper.x,imageBounds.upper.x),std::max(sourceBounds.upper.y,imageBounds.upper.y)}};
        const auto add=[&](auto& degree,Point begin,Point end){if(a.orientation<0)std::swap(begin,end);++degree[Key(begin.x,begin.y)][1];++degree[Key(end.x,end.y)][0];};
        add(sourceDegree,first,last);add(imageDegree,C::position(image,0),C::position(image,1));data->arcs.push_back(std::move(a));
    }
    for(const auto* degree:{&sourceDegree,&imageDegree})for(const auto& [key,count]:*degree){(void)key;require(count[0]==1&&count[1]==1,"material arcs are not a closed directed C0 chain");}
    return Prepared(std::move(data));
}
std::size_t Prepared::arcCount()const noexcept{return data_->arcs.size();}
std::vector<CellMotion> Prepared::evaluate(const std::vector<Rectangle>& boxes,const Config& cfg)const
{
    require(std::isfinite(cfg.relativeTolerance)&&cfg.relativeTolerance>0&&cfg.maximumDepth>0&&cfg.maximumDepth<=32&&
        std::isfinite(cfg.areaVolumeRelativeTolerance)&&cfg.areaVolumeRelativeTolerance>0&&cfg.areaVolumeRelativeTolerance<1,
        "invalid quadrature configuration");
    std::vector<CellMotion> result;result.reserve(boxes.size());
    for(const auto& box:boxes){const long double width=box.upper.x-box.lower.x,height=box.upper.y-box.lower.y;
        require(width>0&&height>0&&std::isfinite(width*height),"invalid native rectangle");CellMotion out;Sum contraction,area,volume,independent;
        for(std::size_t index=0;index<data_->arcs.size();++index){const auto& arc=data_->arcs[index];
            if(arc.upper[0]<box.lower.x||arc.lower[0]>box.upper.x||arc.upper[1]<box.lower.y||arc.lower[1]>box.upper.y)continue;
            ArcCovector row;row.arc=index;row.sample=instantGradient(at(arc,box,.5L),arc.orientation,width,height);out.gradient.push_back(row);
            for(unsigned j=0;j<3;++j)for(unsigned axis=0;axis<2;++axis){const auto d=axis==0?arc.samples[j].x:arc.samples[j].y;contraction.add(row.sample[j][axis]*d);}
            const auto av=bubble(at(arc,box,.5L),arc.orientation,width,height,cfg);area.add(av.area);volume.add(av.volume);out.areaQuadratureErrorEstimate+=av.error;
            out.quadratureSamples+=1+av.samples;
        }
        if(!out.gradient.empty())for(const auto& arc:data_->arcs){
            independent.add(clippedGreen(at(arc,box,1),arc.orientation,width,height)-clippedGreen(at(arc,box,0),arc.orientation,width,height));}
        out.nativeRoundoffEstimate=1024*std::numeric_limits<long double>::epsilon()*width*height*data_->arcs.size();
        // Gonzalez discrete gradient in the original quadratic sample L2
        // metric. A single midpoint shape derivative replaces theta quadrature;
        // the independent closed/clipped volume difference supplies its finite
        // secant. This preserves the zero-displacement physical derivative.
        long double motionScale=0;
        for(const auto& row:out.gradient)for(const auto& d:data_->arcs[row.arc].samples)
            motionScale=std::max({motionScale,std::abs(static_cast<long double>(d.x)),std::abs(static_cast<long double>(d.y))});
        const auto remainder=independent.value-contraction.value;
        if(motionScale>0&&std::abs(remainder)>out.nativeRoundoffEstimate){
            constexpr long double gram[3][3]={{4.L/30,2.L/30,-1.L/30},{2.L/30,16.L/30,2.L/30},{-1.L/30,2.L/30,4.L/30}};
            std::vector<Gradient> metricDisplacement(out.gradient.size());Sum denominator;
            for(std::size_t ri=0;ri<out.gradient.size();++ri){const auto& samples=data_->arcs[out.gradient[ri].arc].samples;
                for(unsigned axis=0;axis<2;++axis){std::array<long double,3> d{};
                    for(unsigned k=0;k<3;++k)d[k]=(axis==0?samples[k].x:samples[k].y)/motionScale;
                    for(unsigned j=0;j<3;++j){Sum md;for(unsigned k=0;k<3;++k)md.add(gram[j][k]*d[k]);metricDisplacement[ri][j][axis]=md.value;denominator.add(d[j]*md.value);}}}
            require(std::isfinite(denominator.value)&&denominator.value>0,"degenerate physical quadratic displacement metric");
            const auto factor=(remainder/motionScale)/denominator.value;require(std::isfinite(factor),"unrepresentable finite volume discrete gradient");
            for(std::size_t ri=0;ri<out.gradient.size();++ri)for(unsigned j=0;j<3;++j)for(unsigned axis=0;axis<2;++axis)
                out.gradient[ri].sample[j][axis]+=factor*metricDisplacement[ri][j][axis];
            contraction=Sum{};
            for(const auto& row:out.gradient)for(unsigned j=0;j<3;++j)for(unsigned axis=0;axis<2;++axis)
                contraction.add(row.sample[j][axis]*(axis==0?data_->arcs[row.arc].samples[j].x:data_->arcs[row.arc].samples[j].y));
        }
        // An unresolved secant is retained as numerical error, not divided by a
        // tiny velocity to create an unphysical load or declared exactly zero.
        out.displacementVolumeChange=contraction.value;out.independentVolumeDifference=independent.value;
        out.volumeIdentityResidual=out.displacementVolumeChange-out.independentVolumeDifference;
        out.areaDerivative=area.value;out.volumeDerivative=volume.value;
        if(volume.value>0){out.areaVolumeAvailable=true;out.areaVolumeCovector=area.value/volume.value;require(std::isfinite(out.areaVolumeCovector),"nonfinite physical area-volume covector");}
        require(std::abs(out.volumeIdentityResidual)<=out.quadratureErrorEstimate+cfg.relativeTolerance*width*height+out.nativeRoundoffEstimate,"independent clipped-volume work identity exceeds numerical budget");
        result.push_back(std::move(out));
    }
    return result;
}
}
