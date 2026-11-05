#include "Core/Geometry/ParametricHermiteCurve.h"
#include <sstream>
#include <iomanip>

namespace LAMP::Geometry::Parametric {
void validateContainedSurface(const Curve& curve,const DomainRectangle& box,bool rz){
    box.validate();
    for(const auto& arc:curve.arcs)for(int axis=0;axis<2;axis++){
        const auto extent=range(arc.coordinate(axis));const double lo=axis?box.bottom:box.left,hi=axis?box.top:box.right;
        const double tolerance=2e-11*(hi-lo)+64*std::numeric_limits<double>::epsilon()*std::max(std::abs(lo),std::abs(hi));
        bool lower=rz&&axis==0&&lo==0&&((curve.boundary.axis(0)&&arc.a==0)||(curve.boundary.axis(1)&&arc.b==1)),upper=false;
        for(int end=0;end<2;end++)if(curve.boundary.fixed(end)&&(end?arc.b==1:arc.a==0)){
            const auto p=curve.boundary.position[end];const double coordinate=axis?p.y:p.x;lower=lower||coordinate==lo;upper=upper||coordinate==hi;}
        check((lower?extent[0]>=lo-tolerance:extent[0]>lo)&&(upper?extent[1]<=hi+tolerance:extent[1]<hi),"free surface has undeclared domain contact or exit");
    }
}
SurfaceRestriction restrictSurface(const Curve& curve,const DomainRectangle& box){
    box.validate();check(!curve.arcs.empty(),"empty free surface restriction");SurfaceRestriction result;
    for(const auto& arc:curve.arcs){std::vector<std::array<double,2>> pieces;
        for(const auto& ab:rectangleIntervals(arc,box)){
            if(!pieces.empty()&&pieces.back()[1]==ab[0])pieces.back()[1]=ab[1];else pieces.push_back(ab);}
        result.clipped=result.clipped||pieces.size()!=1||pieces.front()[0]!=0||pieces.front()[1]!=1;
        result.intervals.push_back(std::move(pieces));}
    const auto begins=[&](std::size_t i){return !result.intervals[i].empty()&&result.intervals[i].front()[0]==0;};
    const auto finishes=[&](std::size_t i){return !result.intervals[i].empty()&&result.intervals[i].back()[1]==1;};
    for(std::size_t i=0;i<curve.arcs.size();i++)for(const auto& ab:result.intervals[i]){
        const bool originalFirst=curve.boundary.open()&&i==0&&ab[0]==0;
        const bool originalLast=curve.boundary.open()&&i+1==curve.arcs.size()&&ab[1]==1;
        if(!originalFirst&&(ab[0]>0||!finishes(i?i-1:curve.arcs.size()-1)))result.ends.push_back({i,ab[0],-1});
        if(!originalLast&&(ab[1]<1||!begins((i+1)%curve.arcs.size())))result.ends.push_back({i,ab[1],1});
    }
    return result;
}
namespace {
using Scalar=std::array<double,13>;
Scalar derivative(const Scalar& p){Scalar out{};for(int k=1;k<13;k++)out[k-1]=k*p[k];return out;}
Scalar product(const Scalar& a,const Scalar& b){
    Scalar out{};for(int i=0;i<13;i++)if(a[i]!=0)for(int j=0;j<13;j++)if(b[j]!=0){
        check(i+j<13,"surface stationary polynomial degree overflow");out[i+j]+=a[i]*b[j];}return out;
}
Scalar withoutEndpointFactor(const Scalar& p,bool last){
    Scalar out{};
    if(!last){for(int i=0;i<12;i++)out[i]=p[i+1];}
    else {out[11]=-p[12];for(int i=10;i>=0;i--)out[i]=out[i+1]-p[i+1];}
    return out;
}
void appendRoots(std::vector<double>& points,Scalar polynomial){
    double scale=0;for(double v:polynomial){check(std::isfinite(v),"nonfinite surface stationary polynomial");scale=std::max(scale,std::abs(v));}
    if(scale==0)return;
    for(double& v:polynomial)v/=scale;
    const auto found=roots(polynomial,0,1);points.insert(points.end(),found.begin(),found.end());
}
}

std::vector<double> surfacePartition(const Arc& arc,const CurveBoundary& boundary,bool rz){
    check(arc.b>arc.a,"surface partition requires a positive material interval");
    double scale=0;for(int k=1;k<5;k++)scale=std::max({scale,std::abs(arc.p[k].x),std::abs(arc.p[k].y)});
    check(scale>0&&std::isfinite(scale),"surface partition requires a regular finite arc");
    std::array<Point,5> p{};for(int k=1;k<5;k++)p[k]=arc.p[k]/scale;
    Scalar dx{},dy{},numerator{};
    for(int k=1;k<5;k++){dx[k-1]=k*p[k].x;dy[k-1]=k*p[k].y;}
    // cross(R',R'') has degree four, not five. Pair antisymmetric terms
    // analytically so the leading cancellation does not create false roots.
    for(int i=1;i<5;i++)for(int j=i+1;j<5;j++)numerator[i+j-3]+=i*j*(j-i)*cross(p[i],p[j]);
    auto speed2=product(dx,dx);const auto yy=product(dy,dy);
    for(int k=0;k<13;k++)speed2[k]+=yy[k];
    const auto speedDerivative=derivative(speed2);
    auto meridional=product(derivative(numerator),speed2);const auto denominatorTerm=product(numerator,speedDerivative);
    for(int k=0;k<13;k++)meridional[k]=2*meridional[k]-3*denominatorTerm[k];
    std::vector<double> result{0,1};appendRoots(result,speedDerivative);appendRoots(result,numerator);appendRoots(result,meridional);
    if(rz){
        Scalar radius{};double radialScale=0;for(const auto& v:arc.p)radialScale=std::max(radialScale,std::abs(v.x));
        check(radialScale>0&&std::isfinite(radialScale),"free surface cannot lie wholly on the RZ axis");
        for(int k=0;k<5;k++)radius[k]=arc.p[k].x/radialScale;
        // k_theta = z'/(r sqrt(Q)). At a declared regular pole both r and
        // z' vanish. Use their analytic common factor, never an epsilon r.
        auto axial=dy;
        for(int end=0;end<2;end++)if(boundary.axis(end)&&(end?arc.b==1:arc.a==0)){
            sampleSurfaceRZ(arc,double(end),true);
            radius=withoutEndpointFactor(radius,end!=0);axial=withoutEndpointFactor(axial,end!=0);
        }
        const auto dr=derivative(radius);
        auto hoop=product(product(derivative(axial),radius),speed2);
        const auto radialTerm=product(product(axial,dr),speed2),speedTerm=product(product(axial,radius),speedDerivative);
        for(int k=0;k<13;k++)hoop[k]=2*hoop[k]-2*radialTerm[k]-speedTerm[k];
        appendRoots(result,hoop);
    }
    std::sort(result.begin(),result.end());result.erase(std::unique(result.begin(),result.end()),result.end());
    // Stationary speed points also expose a nonregular chart before it can
    // be hidden between a quadrature rule's fixed abscissae.
    for(double t:result)sampleSurface(arc,t);
    return result;
}

SurfaceSample sampleSurface(const Arc& arc,double fraction) {
    check(std::isfinite(fraction) && fraction>=0 && fraction<=1,"parameter surface sample outside its arc");
    const auto d=localValue(arc.p,fraction,1),dd=localValue(arc.p,fraction,2);const double speed=norm(d);
    check(speed>0 && std::isfinite(speed),"nonregular parameter surface sample");
    const double curvature=cross(d,dd)/(speed*speed*speed);
    check(std::isfinite(curvature),"nonfinite parameter surface curvature");
    return {arc.sampleFraction(fraction),d,{d.y/speed,-d.x/speed},speed,curvature};
}
SurfaceSample sampleSurfaceRZ(const Arc& arc,double fraction,bool pole) {
    auto q=sampleSurface(arc,fraction);
    if(pole){check(fraction==0||fraction==1,"regular pole must be a material endpoint");
        check(std::abs(q.point.x)<2e-11&&std::abs(q.derivative.y)<2e-8*q.speed,"nonregular axis surface endpoint");
        q.point.x=0;q.derivative.y=0;q.normal.x=0;q.curvature*=2;
    }else{
        if(!(q.point.x>0)){
            std::ostringstream message;message<<std::setprecision(17)<<"RZ surface has nonpositive radius away from a declared pole: radius="
                <<q.point.x<<" fraction="<<fraction<<" arc="<<arc.a<<":"<<arc.b<<" coefficients=";
            for(const auto& p:arc.p)message<<p.x<<","<<p.y<<";";
            throw std::runtime_error(message.str());
        }
        q.curvature+=q.normal.x/q.point.x;
    }
    check(std::isfinite(q.curvature),"nonfinite RZ surface curvature");return q;
}
SurfaceSample sampleSurfaceRZ(const Arc& arc,double fraction,const CurveBoundary& boundary) {
    const bool first=boundary.axis(0)&&arc.a==0,last=boundary.axis(1)&&arc.b==1;
    if(!first&&!last)return sampleSurfaceRZ(arc,fraction,false);
    // The endpoint conditions are owned by the boundary, not inferred from
    // a small radius. Validate them before canceling their exact factors.
    if(first)sampleSurfaceRZ(arc,0,true);
    if(last)sampleSurfaceRZ(arc,1,true);
    auto q=sampleSurface(arc,fraction);Scalar radius{},axial{};
    for(int k=0;k<5;k++)radius[k]=arc.p[k].x;
    for(int k=1;k<5;k++)axial[k-1]=k*arc.p[k].y;
    for(int end=0;end<2;end++)if(end?last:first){radius=withoutEndpointFactor(radius,end!=0);axial=withoutEndpointFactor(axial,end!=0);}
    const auto value=[&](const Scalar& p){double v=p.back();for(int k=11;k>=0;k--)v=v*fraction+p[k];return v;};
    const double r=value(radius),dr=value(derivative(radius)),ddr=value(derivative(derivative(radius)));
    const double z=value(axial),dz=value(derivative(axial));
    const double factor=first?(last?fraction*(1-fraction):fraction):(1-fraction);
    const double df=first?(last?1-2*fraction:1):-1,ddf=first&&last?-2:0;
    check(r>0&&std::isfinite(r),"RZ pole quotient has nonpositive radius in its declared arc");
    q.point.x=factor*r;q.derivative={df*r+factor*dr,factor*z};
    const Point second{ddf*r+2*df*dr+factor*ddr,df*z+factor*dz};
    q.speed=norm(q.derivative);check(q.speed>0&&std::isfinite(q.speed),"nonregular factored RZ pole jet");
    q.normal={q.derivative.y/q.speed,-q.derivative.x/q.speed};
    // Both curvature terms have the same finite pole limit. Away from the
    // endpoint, the quotient z/r is evaluated after analytic cancellation;
    // there is no radial cutoff or transition between geometric models.
    q.curvature=cross(q.derivative,second)/(q.speed*q.speed*q.speed)+z/(q.speed*r);
    check(std::isfinite(q.curvature),"nonfinite factored RZ surface curvature");return q;
}
DistanceSample signedDistance(const Curve& curve,Point query) {
    check(std::isfinite(query.x) && std::isfinite(query.y),"nonfinite parametric distance query");
    struct Foot {double square=0,t=0;Point point;std::size_t arc=0;};
    std::vector<Foot> candidates;double best=std::numeric_limits<double>::infinity();
    for(std::size_t k=0;k<curve.arcs.size();++k) {
        const auto& arc=curve.arcs[k];std::array<double,8> stationary{};
        for(int i=0;i<5;++i)for(int j=1;j<5;++j)stationary[i+j-1]+=j*dot(i==0?arc.p[i]-query:arc.p[i],arc.p[j]);
        auto points=roots(stationary,0,1);points.push_back(0);points.push_back(1);
        for(double t:points){const Point foot=arc.sampleFraction(t),delta=foot-query;const double square=dot(delta,delta);
            best=std::min(best,square);candidates.push_back({square,t,foot,k});}
    }
    check(std::isfinite(best),"parametric distance query has no boundary");
    const double distance=std::sqrt(best),positionBudget=128*std::numeric_limits<double>::epsilon()*std::max(1.,norm(query));
    const double squaredBudget=2*distance*positionBudget+positionBudget*positionBudget;
    const Foot* selected=nullptr;bool ambiguous=false,node=false;
    for(const auto& foot:candidates)if(foot.square<=best+squaredBudget) {
        if(!selected)selected=&foot;else if(norm(foot.point-selected->point)>positionBudget)ambiguous=true;
        node=node || foot.t==0 || foot.t==1;
    }
    PhaseMeasure phase(curve);bool inside=false;
    for(const auto& interval:phase.intervals(query.y))if(query.x>=interval[0] && query.x<interval[1])inside=true;
    DistanceSample result;result.value=(inside?-1.:1.)*distance;result.foot=selected->point;
    if(ambiguous || node)return result;
    const auto& arc=curve.arcs[selected->arc];const Point d=arc.sampleFraction(selected->t,1),dd=arc.sampleFraction(selected->t,2);
    const double speed=norm(d);if(!(speed>0 && std::isfinite(speed)))return result;
    const double curvature=cross(d,dd)/(speed*speed*speed),offsetFactor=1+result.value*curvature;
    if(!(offsetFactor>0 && std::isfinite(curvature)))return result;
    result.curvature=curvature;result.distanceCurvature=curvature/offsetFactor;
    result.normal={d.y/speed,-d.x/speed};result.tangent=d/speed;result.differentialAvailable=true;return result;
}
} // namespace LAMP::Geometry::Parametric
