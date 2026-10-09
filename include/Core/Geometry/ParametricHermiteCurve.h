#pragma once
#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
#include <vector>

// Polynomial builders for a single oriented, regular parametric boundary.
// A persistent owner publishes these only through const access. Phase integrals
// require the effective liquid boundary, not overlapping raw component curves.
// Polynomial roots are numerical subdivision roots, not topology certificates.
namespace Quantum::Geometry::Parametric {
inline void check(bool ok,const char* message){if(!ok)throw std::runtime_error(message);}
constexpr std::array<double,5> gx{{-.9061798459386639928,-.5384693101056830910,0,.5384693101056830910,.9061798459386639928}};
constexpr std::array<double,5> gw{{.2369268850561890875,.4786286704993664680,.5688888888888888889,.4786286704993664680,.2369268850561890875}};
struct GaussEight {
    std::array<double,8> x{},w{};
    GaussEight() {
        const double pi=std::acos(-1.);
        for(int i=0;i<8;++i) {
            double z=std::cos(pi*(i+.75)/8.5),derivative=0;
            for(int iteration=0;iteration<20;++iteration) {
                double previous=1,current=z;
                for(int n=2;n<=8;++n){const double next=((2*n-1)*z*current-(n-1)*previous)/n;previous=current;current=next;}
                derivative=8*(z*current-previous)/(z*z-1);const double next=z-current/derivative;
                if(next==z)break;z=next;
            }
            x[i]=z;w[i]=2/((1-z*z)*derivative*derivative);
        }
    }
};
inline const GaussEight& gaussEight(){static const GaussEight rule;return rule;}
template<std::size_t N>
double polynomial(const std::array<double,N>& p,double t,int d=0)
{
    double out=0;
    for(int n=static_cast<int>(N)-1;n>=d;--n){double c=p[n];for(int j=0;j<d;++j)c*=n-j;out=out*t+c;}
    return out;
}
// Derivative-root subdivision finds all simple roots on this graph control.
// It does not certify a degenerate contact/topology event in production.
template<std::size_t N>
std::vector<double> roots(const std::array<double,N>& p,double lo,double hi,int degree=static_cast<int>(N)-1)
{
    while(degree>0 && p[degree]==0)--degree;
    if(degree==0)return {};
    std::array<double,N> derivative{};for(int i=1;i<=degree;++i)derivative[i-1]=i*p[i];
    auto points=roots(derivative,lo,hi,degree-1);points.insert(points.begin(),lo);points.push_back(hi);
    std::vector<double> out;
    for(std::size_t i=1;i<points.size();++i) {
        const auto value=[&](double x){double out=0;for(int j=degree;j>=0;--j)out=out*x+p[j];return out;};
        double a=points[i-1],b=points[i],fa=value(a),fb=value(b);
        if(fa==0 && a>lo && a<hi)out.push_back(a);
        if((fa<0 && fb>0)||(fa>0 && fb<0)) {
            for(int n=0;n<64;++n){const double m=.5*(a+b);if(m==a || m==b)break;const double fm=value(m);
                if(fm==0){a=b=m;break;}if((fa<0)==(fm<0)){a=m;fa=fm;}else b=m;}
            out.push_back(.5*(a+b));
        }
    }
    std::sort(out.begin(),out.end());out.erase(std::unique(out.begin(),out.end()),out.end());return out;
}
struct Point {
    double x=0,y=0;
    Point operator+(Point b)const{return {x+b.x,y+b.y};}
    Point operator-(Point b)const{return {x-b.x,y-b.y};}
    Point operator*(double s)const{return {s*x,s*y};}
    Point operator/(double s)const{return {x/s,y/s};}
    Point& operator+=(Point b){x+=b.x;y+=b.y;return *this;}
};
inline double dot(Point a,Point b){return a.x*b.x+a.y*b.y;}
inline double cross(Point a,Point b){return a.x*b.y-a.y*b.x;}
inline double norm(Point a){return std::sqrt(dot(a,a));}
template<std::size_t N>
struct PolynomialArc {
    double a=0,b=0;std::array<Point,N> p{};
    Point sample(double s,int derivative=0)const {
        return sampleFraction((s-a)/(b-a),derivative);
    }
    Point sampleFraction(double t,int derivative=0)const {
        const double length=b-a;Point out;
        for(int i=static_cast<int>(N)-1;i>=derivative;--i){double factor=1;for(int j=0;j<derivative;++j)factor*=(i-j)/length;out=out*t+p[i]*factor;}return out;
    }
    Point increment(double lo,double hi)const {
        const double a0=(lo-a)/(b-a),b0=(hi-a)/(b-a),dt=(hi-lo)/(b-a);Point out;
        for(int n=1;n<static_cast<int>(N);++n){double sum=0;for(int j=0;j<n;++j)sum+=std::pow(a0,j)*std::pow(b0,n-1-j);out+=p[n]*(dt*sum);}return out;
    }
    std::array<double,N> coordinate(int axis)const {std::array<double,N> out{};for(std::size_t i=0;i<N;++i)out[i]=axis==0?p[i].x:p[i].y;return out;}
};
using Arc=PolynomialArc<5>;
// Physical free-curve endpoints, separate from the domain-boundary path
// closing its liquid region. Closure segments carry phase measure only;
// they are never free-surface arcs, curvature or capillary energy.
struct CurveBoundary {
    enum class Kind {Periodic,Axis,Fixed};
    std::array<Kind,2> kind{{Kind::Periodic,Kind::Periodic}};
    std::array<Point,2> position{};
    // Ordered from the final endpoint back to the initial endpoint.
    std::vector<Point> corners;
    bool open() const {return kind[0]!=Kind::Periodic;}
    bool axis(int end) const {return kind.at(end)==Kind::Axis;}
    bool fixed(int end) const {return kind.at(end)==Kind::Fixed;}
    bool axisPair() const {return axis(0)&&axis(1);}
    bool constantMode(int component) const {return !open()||(component==1&&axisPair());}
    static CurveBoundary axisMeridian(){CurveBoundary out;out.kind={Kind::Axis,Kind::Axis};return out;}
    static CurveBoundary suppliedMeridian(Point rim,Point corner){
        CurveBoundary out;out.kind={Kind::Fixed,Kind::Axis};out.position[0]=rim;out.corners={corner};return out;
    }
    CurveBoundary homogeneous() const {auto out=*this;out.position={};out.corners.clear();return out;}
};
template<std::size_t N>
struct PolynomialCurve {
    std::vector<PolynomialArc<N>> arcs;
    CurveBoundary boundary;
    std::size_t index(double s)const {
        const auto it=std::upper_bound(arcs.begin(),arcs.end(),s,[](double x,const PolynomialArc<N>& p){return x<p.b;});
        return it==arcs.end()?arcs.size()-1:static_cast<std::size_t>(it-arcs.begin());
    }
    Point sample(double s,int d=0)const{
        Point p=arcs[index(s)].sample(s,d);
        if(s==0||s==1){const int end=s==1;
            if(boundary.axis(end)){if(d==0)p.x=0;if(d==1)p.y=0;}
            if(boundary.fixed(end)&&d==0)p=boundary.position[end];}return p;
    }
    std::vector<PolynomialArc<N>> closure() const {
        std::vector<PolynomialArc<N>> out;if(!boundary.open())return out;check(!arcs.empty(),"empty open free curve");
        Point previous=sample(1);auto corners=boundary.corners;corners.push_back(sample(0));
        for(const auto p:corners){if(p.x!=previous.x||p.y!=previous.y){PolynomialArc<N> line;line.a=0;line.b=1;line.p[0]=previous;line.p[1]=p-previous;out.push_back(line);}previous=p;}
        return out;
    }
    void validateBoundary(double left,double right,double bottom,double top,bool rz) const {
        check(!arcs.empty(),"empty parametric boundary");
        check((boundary.kind[0]==CurveBoundary::Kind::Periodic)==(boundary.kind[1]==CurveBoundary::Kind::Periodic),"periodic free curve has one open endpoint");
        if(!boundary.open()){check(boundary.corners.empty(),"closed free curve has domain closure corners");return;}
        for(int end=0;end<2;end++){
            check(boundary.axis(end)||boundary.fixed(end),"unsupported free-curve endpoint");
            const auto p=sample(end),raw=end?arcs.back().sampleFraction(1):arcs.front().p[0];
            check(std::isfinite(p.x)&&std::isfinite(p.y)&&norm(raw-p)<2e-11,"free endpoint does not match its physical owner");
            if(boundary.axis(end))check(rz&&left==0&&p.x==0,"axis endpoint requires a domain containing the RZ axis");
        }
        const auto edges=closure();
        // A moving pole may change only the zero-area axial closure. Every
        // positive-area boundary trace is bounded by fixed points/corners.
        if(boundary.axis(1)&&!edges.empty())check(edges.front().p[0].x==0&&edges.front().p[1].x==0,"moving pole closure must follow the axis");
        if(boundary.axis(0)&&!edges.empty())check(edges.back().p[0].x==0&&edges.back().p[1].x==0,"moving pole closure must follow the axis");
        for(const auto& edge:edges){
            const Point a=edge.p[0],b=a+edge.p[1];
            // A virtual axis closure carries exactly zero RZ measure. Its
            // material poles can pass an axial outlet; the physical region
            // is still obtained by rectangle restriction of this one curve.
            if(rz&&left==0&&a.x==0&&b.x==0){check(std::isfinite(a.y)&&std::isfinite(b.y),"nonfinite axial phase closure");continue;}
            for(auto p:{a,b})check(std::isfinite(p.x)&&std::isfinite(p.y)&&p.x>=left&&p.x<=right&&p.y>=bottom&&p.y<=top,"phase closure leaves the physical domain");
            check((a.x==b.x&&(a.x==left||a.x==right))||(a.y==b.y&&(a.y==bottom||a.y==top)),"phase closure is not a physical domain edge");
        }
    }
    Point sampleInterval(double a,double b,double fraction,int derivative=0)const {
        std::size_t i=index(a);
        while(i+1<arcs.size() && fraction>(arcs[i].b-a)/(b-a))++i;
        const auto& arc=arcs[i];const double t=(a-arc.a)/(arc.b-arc.a)+(b-a)/(arc.b-arc.a)*fraction;
        return arc.sampleFraction(t,derivative);
    }
    Point increment(double a,double b)const {
        Point out;for(std::size_t i=index(a);i<arcs.size() && arcs[i].a<b;++i){const double lo=std::max(a,arcs[i].a),hi=std::min(b,arcs[i].b);if(lo<hi)out+=arcs[i].increment(lo,hi);}return out;
    }
    std::vector<double> partition(double a,double b)const {
        std::vector<double> out{a};for(std::size_t i=index(a);i<arcs.size() && arcs[i].b<b;++i)if(arcs[i].b>a)out.push_back(arcs[i].b);out.push_back(b);return out;
    }
    double area()const {
        static_assert(N==5,"This five-point whole-curve rule is only for quartic geometry; transported higher-degree areas use their own exact rule.");
        double out=0;const auto add=[&](const PolynomialArc<N>& a){for(int q=0;q<5;++q){const double t=.5*(1+gx[q]);out+=.25*(a.b-a.a)*gw[q]*cross(a.sampleFraction(t)-arcs.front().p[0],a.sampleFraction(t,1));}};
        for(const auto& a:arcs)add(a);for(const auto& a:closure())add(a);return out;
    }
    double energy()const {
        static_assert(N==5,"This five-point energy rule is only for quartic geometry.");
        double out=0;for(const auto& a:arcs)for(int q=0;q<5;++q){const Point dd=a.sampleFraction(.5*(1+gx[q]),2);out+=.5*(a.b-a.a)*gw[q]*dot(dd,dd);}return out;
    }
};
using Curve=PolynomialCurve<5>;
template<std::size_t N>
std::array<double,2> range(const std::array<double,N>& p) {
    std::array<double,N> derivative{};for(std::size_t i=1;i<N;++i)derivative[i-1]=i*p[i];
    auto candidates=roots(derivative,0,1);candidates.push_back(0);candidates.push_back(1);
    std::array<double,2> out{{polynomial(p,0),polynomial(p,0)}};
    for(double s:candidates){const double value=polynomial(p,s);out[0]=std::min(out[0],value);out[1]=std::max(out[1],value);}return out;
}
template<std::size_t N>
PolynomialArc<N> restrictedArc(const PolynomialArc<N>& source,double a,double b) {
    PolynomialArc<N> out;out.a=a;out.b=b;
    const double origin=(a-source.a)/(source.b-source.a),width=(b-a)/(source.b-source.a);
    // Compose in local coordinates, retaining coefficients even for a
    // subinterval narrower than the global parameter's internal point ULP.
    for(int k=0;k<static_cast<int>(N);++k) {
        double choose=1,power=1;
        for(int j=k;j>=0;--j) {
            if(j<k){choose*=double(j+1)/(k-j);power*=origin;}
            out.p[j]+=source.p[k]*(choose*power*std::pow(width,j));
        }
    }
    return out;
}
struct DomainRectangle {
    double left=0,right=1,bottom=0,top=1;
    bool contains(Point p) const {return p.x>=left&&p.x<=right&&p.y>=bottom&&p.y<=top;}
    void validate() const {check(std::isfinite(left)&&std::isfinite(right)&&std::isfinite(bottom)&&std::isfinite(top)
        &&right>left&&top>bottom,"invalid physical surface rectangle");}
};
// Parameter intervals of the ORIGINAL polynomial inside a physical rectangle.
// No new curve, snapped control points or sampled phase indicator is formed.
// Numerical polynomial roots retain their existing non-certified semantics.
template<std::size_t N>
std::vector<std::array<double,2>> rectangleIntervals(const PolynomialArc<N>& arc,const DomainRectangle& box) {
    box.validate();const auto x=range(arc.coordinate(0)),y=range(arc.coordinate(1));
    if(x[1]<box.left||x[0]>box.right||y[1]<box.bottom||y[0]>box.top)return {};
    std::vector<double> cuts{0,1};
    for(int axis=0;axis<2;axis++)for(double line:axis?std::array<double,2>{box.bottom,box.top}:std::array<double,2>{box.left,box.right}){
        auto equation=arc.coordinate(axis);equation[0]-=line;const auto found=roots(equation,0,1);cuts.insert(cuts.end(),found.begin(),found.end());}
    std::sort(cuts.begin(),cuts.end());cuts.erase(std::unique(cuts.begin(),cuts.end()),cuts.end());
    std::vector<std::array<double,2>> result;
    for(std::size_t i=1;i<cuts.size();i++)if(box.contains(arc.sampleFraction(.5*(cuts[i-1]+cuts[i]))))result.push_back({cuts[i-1],cuts[i]});
    return result;
}
struct SurfaceRestriction {
    struct End {std::size_t arc=0;double fraction=0;int sign=0;};
    std::vector<std::vector<std::array<double,2>>> intervals;
    // Only newly exposed ends: + at an outgoing parameter end, - at an
    // incoming one. Original fixed endpoints and RZ poles are not new ends.
    std::vector<End> ends;
    bool clipped=false;
};
SurfaceRestriction restrictSurface(const Curve&,const DomainRectangle&);
// Closed-domain admission with the declared endpoint incidence tolerance.
// Shared by force construction and both physical candidate time layers.
void validateContainedSurface(const Curve&,const DomainRectangle&,bool axisymmetric);
template<std::size_t N>
Point localValue(const std::array<Point,N>& coefficients,double t,int derivative=0) {
    Point out;for(int k=static_cast<int>(N)-1;k>=derivative;--k){double factor=1;for(int j=0;j<derivative;++j)factor*=k-j;out=out*t+coefficients[k]*factor;}return out;
}
struct CrossingEvent {double coordinate=0;int jump=0;};
enum class RectanglePhase {Empty,Full,Cut};
struct RectanglePhaseMoments {
    std::array<double,6> values{};
    RectanglePhase phase=RectanglePhase::Empty;
};
template<std::size_t N>
struct PolynomialPhaseMeasure {
    const PolynomialCurve<N>& curve;
    struct Bounds {std::array<double,2> x,y;};
    std::vector<Bounds> bounds;
    std::vector<Point> sharedEnds;
    std::vector<PolynomialArc<N>> closureEdges;
    std::size_t endpointIncidences=0;
    double maximumEndpointIncidenceWidth=0;
    std::array<std::map<double,std::vector<std::array<double,2>>>,4> lineIntervals;
    explicit PolynomialPhaseMeasure(const PolynomialCurve<N>& value):curve(value),closureEdges(value.closure()) {
        std::size_t component=0;
        for(std::size_t i=0;i<curve.arcs.size();++i) {
            const auto& arc=curve.arcs[i];if(arc.a==0)component=i;
            const std::size_t next=arc.b==1?component:i+1;check(next<curve.arcs.size(),"phase arc has no shared endpoint owner");
            const Point end=curve.boundary.open()&&arc.b==1?curve.sample(1):curve.arcs[next].p[0];sharedEnds.push_back(end);
            check(norm(arc.sampleFraction(1)-end)<2e-11,"phase boundary shared endpoint exceeds its declared position tolerance");
            auto x=range(arc.coordinate(0)),y=range(arc.coordinate(1));x[0]=std::min(x[0],end.x);x[1]=std::max(x[1],end.x);y[0]=std::min(y[0],end.y);y[1]=std::max(y[1],end.y);bounds.push_back({x,y});
        }
    }
    const std::vector<std::array<double,2>>& intervals(double y) {return intervals(1,y,false);}
    // Exact one-sided trace; coincident closures need a branch convention,
    // not an epsilon displacement of the physical face or phase geometry.
    const std::vector<std::array<double,2>>& intervals(int axis,double y,bool fromLower=false) {
        check((axis==0||axis==1)&&std::isfinite(y),"invalid phase trace plane");
        auto& cache=lineIntervals[2*axis+int(fromLower)];
        const auto cached=cache.find(y);if(cached!=cache.end())return cached->second;
        const auto normal=[&](Point p){return axis?p.y:p.x;};
        const auto tangent=[&](Point p){return axis?p.x:p.y;};
        const auto crossesZero=[&](double a,double b){return fromLower?(0>std::min(a,b)&&0<=std::max(a,b)):(0>=std::min(a,b)&&0<std::max(a,b));};
        const auto crosses=[&](double a,double b){return crossesZero(a-y,b-y);};
        std::vector<CrossingEvent> events;
        for(std::size_t i=0;i<curve.arcs.size();++i) {
            const auto& arc=curve.arcs[i];const auto extent=axis?bounds[i].y:bounds[i].x;if(y<extent[0] || y>extent[1])continue;
            const auto coordinate=arc.coordinate(axis);auto equation=coordinate;equation[0]-=y;
            const auto intersections=roots(equation,0,1);
            std::array<double,N> derivative{};for(int j=1;j<static_cast<int>(N);++j)derivative[j-1]=j*coordinate[j];
            auto partition=roots(derivative,0,1);partition.insert(partition.begin(),0);partition.push_back(1);
            const auto point=[&](double t){return t==0?arc.p[0]:(t==1?sharedEnds[i]:arc.sampleFraction(t));};
            // Shared nodes own incidence. Opposite one-sided traces use
            // [min,max) or (min,max], so a regular join is counted once.
            for(std::size_t p=1;p<partition.size();++p) {
                double a=partition[p-1],b=partition[p];const Point pa=point(a),pb=point(b);
                const double na=normal(pa),nb=normal(pb);
                // Evaluate distance to the plane directly, as in rectangle
                // clipping. Evaluating an absolute coordinate first loses
                // sub-ULP normal distances near a tangent and shifts the
                // apparent crossing by their square root along the curve.
                const double va=a==0?na-y:polynomial(equation,a),vb=b==1?nb-y:polynomial(equation,b);
                if(!crossesZero(va,vb))continue;
                double x=0;
                if(va==0)x=tangent(pa);else if(vb==0)x=tangent(pb);
                else if(b==1 && ((polynomial(equation,b)<=0&&vb>=0)||(polynomial(equation,b)>=0&&vb<=0))) {
                    // Shared-node incidence spans the stored node and its
                    // independently rounded local-chart endpoint. This owns
                    // that endpoint event; it is NOT an interior root with a
                    // relaxed residual threshold. Neither geometry is moved.
                    x=tangent(pb);++endpointIncidences;
                    maximumEndpointIncidenceWidth=std::max(maximumEndpointIncidenceWidth,std::abs(polynomial(equation,b)-vb));
                }
                else {
                    // The phase trace and rectangle clipping must use the
                    // SAME polynomial root. Independent bisections can pick
                    // different points in a rounded zero plateau near a
                    // tangent, leaving an unmatched Green boundary term.
                    const auto root=std::lower_bound(intersections.begin(),intersections.end(),a);
                    check(root!=intersections.end()&&*root<=b,"phase trace has no owned clipping root");
                    const double t=*root;x=tangent(arc.sampleFraction(t));
                    if(std::abs(polynomial(equation,t))>=128*std::numeric_limits<double>::epsilon()*std::max(1.,std::abs(y)))
                        std::cerr<<std::setprecision(17)<<"PHASE_ROOT_BUDGET axis="<<axis<<" plane="<<y<<" value="<<polynomial(coordinate,t)<<" t="<<t<<" arc="<<i<<" parameter="<<arc.a<<':'<<arc.b<<" ends="<<na<<':'<<nb<<'\n';
                    check(std::abs(polynomial(equation,t))<128*std::numeric_limits<double>::epsilon()*std::max(1.,std::abs(y)),"phase line root exceeds floating coordinate budget");
                }
                events.push_back({x,(axis?1:-1)*(vb>va?-1:1)});
            }
        }
        for(std::size_t i=0;i<closureEdges.size();++i){
            const auto a=closureEdges[i].p[0];
            // Closure incidence belongs to the SAME stored nodes as the
            // free curve and adjacent closure edges. Reconstructing the
            // endpoint as a+(b-a) can shift it by an ULP and leave a lone
            // crossing on a plane between those values, even at an RZ pole.
            const auto b=i+1<closureEdges.size()?closureEdges[i+1].p[0]:curve.sample(0);
            const double na=normal(a),nb=normal(b);
            if(crosses(na,nb)){
                const double x=tangent(a)+(tangent(b)-tangent(a))*(y-na)/(nb-na);events.push_back({x,(axis?1:-1)*(nb>na?-1:1)});}
        }
        std::sort(events.begin(),events.end(),[](const CrossingEvent& a,const CrossingEvent& b){return a.coordinate<b.coordinate;});
        std::vector<std::array<double,2>> pieces;int winding=0;double begin=0;
        for(const auto& event:events){const int before=winding;winding+=event.jump;if(before==0 && winding!=0)begin=event.coordinate;
            if(before!=0 && winding==0)pieces.push_back({begin,event.coordinate});}
        if(winding!=0){std::cerr<<std::setprecision(17)<<"PHASE_LINE_OPEN y="<<y<<" winding="<<winding<<" events=";
            for(const auto& event:events)std::cerr<<' '<<event.coordinate<<':'<<event.jump;std::cerr<<'\n';}
        check(winding==0,"phase line does not close: unresolved endpoint/branch incidence");
        return cache.emplace(y,std::move(pieces)).first->second;
    }
    double length(double y,double a,double b) {
        double out=0;for(const auto& p:intervals(y))out+=std::max(0.,std::min(b,p[1])-std::max(a,p[0]));return out;
    }
    // Green integral of the ACTUAL clipped simple region, including the
    // liquid portions of the rectangle edges. Not a sum of arc defects.
    // Requires an actual oriented phase boundary (outer and hole loops),
    // with liquid on the left; raw overlapping/internal arcs are invalid.
    template<class Visitor>
    bool visitRectangle(double left,double right,double bottom,double top,Visitor&& visit) const {
        bool intersects=false;
        const auto add=[&](const PolynomialArc<N>& arc,const Bounds& box){
            if(box.x[1]<left || box.x[0]>right || box.y[1]<bottom || box.y[0]>top)return;
            intersects=true;
            // intervals(y) takes its trace from ABOVE a horizontal edge.
            // Therefore that trace owns a bottom-coincident arc; retaining
            // it here too counts the same boundary twice. Top-coincident
            // arcs remain: together with the above trace they give the
            // interior (below) trace, for either phase orientation.
            bool onBottom=arc.p[0].y==bottom;
            for(std::size_t k=1;k<N;++k)onBottom=onBottom&&arc.p[k].y==0;
            if(onBottom)return;
            for(const auto& interval:rectangleIntervals(arc,{left,right,bottom,top}))visit(arc,interval[0],interval[1]);
        };
        for(std::size_t i=0;i<curve.arcs.size();++i)add(curve.arcs[i],bounds[i]);
        for(const auto& edge:closureEdges)add(edge,{range(edge.coordinate(0)),range(edge.coordinate(1))});
        return intersects;
    }
    double simpleRectangle(double left,double right,double bottom,double top) {
        const double middle=.5*(bottom+top);double out=0;bool interiorBoundary=false;
        visitRectangle(left,right,bottom,top,[&](const PolynomialArc<N>& arc,double a,double b){
                const auto point=arc.sampleFraction(.5*(a+b));
                interiorBoundary=interiorBoundary||(point.x>left&&point.x<right&&point.y>bottom&&point.y<top);
                if constexpr(N==5) {
                    for(int q=0;q<5;++q){const double t=.5*(a+b+(b-a)*gx[q]);out-=.5*(b-a)*(arc.b-arc.a)*gw[q]*(arc.sampleFraction(t).y-middle)*arc.sampleFraction(t,1).x;}
                } else {
                    static_assert(N<=9,"phase area quadrature degree exceeded");const auto& rule=gaussEight();
                    for(int q=0;q<8;++q){const double t=.5*(a+b+(b-a)*rule.x[q]);out-=.5*(b-a)*(arc.b-arc.a)*rule.w[q]*(arc.sampleFraction(t).y-middle)*arc.sampleFraction(t,1).x;}
                }
        });
        // Use the same pure-cell classification as rectanglePhaseMoments.
        // A bounding-box contact alone is not an interior phase boundary:
        // rounded tangencies must not leave an unmatched Green edge term.
        if(!interiorBoundary){const double inside=length(middle,left,right);return inside>.5*(right-left)?(right-left)*(top-bottom):0.;}
        out+=(middle-bottom)*length(bottom,left,right)+(top-middle)*length(top,left,right);return out;
    }
    // Exact polynomial Green integrals on the SAME clipped phase boundary.
    // Order: (1,X,Y,X^2,XY,Y^2), X=(x-cx)/hx, Y=(y-cy)/hy.
    // RZ uses dV=2*pi*x*dA, including axis cells. No sampled indicator,
    // independent interface or correction of an FV inventory is introduced.
    RectanglePhaseMoments rectanglePhaseMoments(double left,double right,double bottom,double top,bool rz=false) {
        static_assert(N==5,"phase material moments currently own quartic geometry");
        check(std::isfinite(left)&&std::isfinite(right)&&std::isfinite(bottom)&&std::isfinite(top)
            &&right>left&&top>bottom&&(!rz||left>=0),"invalid phase moment rectangle");
        const double hx=right-left,hy=top-bottom,cx=left+.5*hx,cy=bottom+.5*hy;
        constexpr int powers[6][2]={{0,0},{1,0},{0,1},{2,0},{1,1},{0,2}};
        // Local polynomial products, degree <=19 for quartic RZ moments.
        // Compose each clipped arc first, avoiding powers of global x,y.
        using Poly=std::array<double,21>;
        const auto product=[](const Poly& a,const Poly& b){Poly c{};for(int i=0;i<21;i++)for(int j=0;j<21-i;j++)c[i+j]+=a[i]*b[j];return c;};
        double planar[4][3]{};bool interiorBoundary=false;
        visitRectangle(left,right,bottom,top,[&](const PolynomialArc<N>& arc,double a,double b){
            const auto middle=arc.sampleFraction(.5*(a+b));
            interiorBoundary=interiorBoundary||(middle.x>left&&middle.x<right&&middle.y>bottom&&middle.y<top);
            const auto piece=restrictedArc(arc,arc.a+a*(arc.b-arc.a),arc.a+b*(arc.b-arc.a));
            Poly x{},y{},dx{};for(int k=0;k<5;k++){x[k]=piece.p[k].x/hx;y[k]=piece.p[k].y/hy;}
            x[0]=(piece.p[0].x-cx)/hx;y[0]=(piece.p[0].y-cy)/hy;
            for(int k=1;k<5;k++)dx[k-1]=k*x[k];
            std::array<Poly,4> xp{},yp{};xp[0][0]=yp[0][0]=1;
            for(int k=1;k<4;k++){xp[k]=product(xp[k-1],x);yp[k]=product(yp[k-1],y);}
            for(int i=0;i<4;i++)for(int j=0;j<3&&i+j<=3;j++){
                const auto integrand=product(product(xp[i],yp[j+1]),dx);double integral=0;
                for(int k=20;k>=0;k--)integral+=integrand[k]/(k+1);
                planar[i][j]-=hx*hy*integral/(j+1);
            }
        });
        for(int edge=0;edge<2;edge++){
            const double y=edge?top:bottom,Y=edge?.5:-.5,sign=edge?1.:-1.;
            for(const auto& interval:intervals(y)){
                const double a=std::max(left,interval[0]),b=std::min(right,interval[1]);if(!(b>a))continue;
                const double A=(a-cx)/hx,B=(b-cx)/hx;
                for(int i=0;i<4;i++)for(int j=0;j<3&&i+j<=3;j++)
                    planar[i][j]+=sign*hx*hy*std::pow(Y,j+1)/(j+1)*(std::pow(B,i+1)-std::pow(A,i+1))/(i+1);
            }
        }
        // Without an interior boundary the connected open rectangle has a
        // single phase. This is a geometric classification, not a small-volume
        // threshold. Publish exact empty/full integrals; roundoff must not
        // create a material sliver in a geometrically pure cell.
        if(!interiorBoundary){
            if(length(cy,left,right)<.5*hx)return {};
            const double volume=hx*hy*(rz?2*std::acos(-1.)*cx:1.);
            return {{{volume,rz?volume*hx/(12*cx):0,0,volume/12,0,volume/12}},RectanglePhase::Full};
        }
        RectanglePhaseMoments result;result.phase=RectanglePhase::Cut;
        for(int k=0;k<6;k++){const int i=powers[k][0],j=powers[k][1];
            result.values[k]=rz?2*std::acos(-1.)*(cx*planar[i][j]+hx*planar[i+1][j]):planar[i][j];}
        return result;
    }
    std::array<double,6> rectangleMoments(double left,double right,double bottom,double top,bool rz=false) {
        return rectanglePhaseMoments(left,right,bottom,top,rz).values;
    }
    double indicatorArea(double tolerance,std::vector<double> cuts={},double left=-2,double right=3,double bottom=-2,double top=3) {
        // Integrate nonzero-winding phase length. Overlap contributes once.
        // Adaptive quadrature is a numerical control, not a general proof
        // that arbitrarily small unresolved arrangement features are absent.
        for(std::size_t i=0;i<curve.arcs.size();++i){cuts.push_back(bounds[i].y[0]);cuts.push_back(bounds[i].y[1]);cuts.push_back(curve.arcs[i].p[0].y);cuts.push_back(sharedEnds[i].y);
            for(double edge:{left,right}){auto p=curve.arcs[i].coordinate(0);p[0]-=edge;for(double t:roots(p,0,1))cuts.push_back(curve.arcs[i].sampleFraction(t).y);}}
        cuts.push_back(bottom);cuts.push_back(top);cuts.erase(std::remove_if(cuts.begin(),cuts.end(),[&](double y){return y<bottom || y>top;}),cuts.end());
        std::sort(cuts.begin(),cuts.end());cuts.erase(std::unique(cuts.begin(),cuts.end()),cuts.end());
        double out=0;const double pi=std::acos(-1.);
        for(std::size_t i=1;i<cuts.size();++i) {
            const double lo=cuts[i-1],hi=cuts[i];if(!(hi>lo))continue;
            const auto value=[&](double t){const double y=lo+(hi-lo)*std::pow(std::sin(.5*pi*t),2);return length(y,left,right)*.5*pi*(hi-lo)*std::sin(pi*t);};
            const auto rule=[&](double a,double b){double sum=0;for(int q=0;q<5;++q)sum+=.5*(b-a)*gw[q]*value(.5*(a+b+(b-a)*gx[q]));return sum;};
            std::function<double(double,double,double,int)> integrate;
            integrate=[&](double a,double b,double budget,int depth){const double m=.5*(a+b),coarse=rule(a,b),fine=rule(a,m)+rule(m,b);
                if(std::abs(coarse-fine)<=budget)return fine;
                if(depth>=18)std::cerr<<std::setprecision(17)<<"PHASE_QUADRATURE_UNRESOLVED lo="<<lo<<" hi="<<hi<<" a="<<a<<" b="<<b<<" difference="<<std::abs(coarse-fine)<<" budget="<<budget<<'\n';
                check(depth<18,"phase-indicator quadrature unresolved");return integrate(a,m,.5*budget,depth+1)+integrate(m,b,.5*budget,depth+1);};
            out+=integrate(0,1,tolerance*(hi-lo),0);
        }
        return out;
    }
};
using PhaseMeasure=PolynomialPhaseMeasure<5>;

struct DistanceSample {
    double value=0,curvature=std::numeric_limits<double>::quiet_NaN();
    double distanceCurvature=std::numeric_limits<double>::quiet_NaN();
    Point foot,normal,tangent;
    // C1 joins need not have a unique second derivative. A node, nonunique
    // footpoint or focal point never publishes an ordinary curvature frame.
    bool differentialAvailable=false;
};
// Differential geometry on ONE retained polynomial arc, at a unit fraction.
// The derivative is with respect to that local fraction, preserving short
// arc scales. At a C1 join the supplied arc specifies the one-sided curvature;
// this does not assert a unique Hessian for a closest-point distance query.
struct SurfaceSample {Point point,derivative,normal;double speed=0,curvature=0;};
SurfaceSample sampleSurface(const Arc&,double fraction);
// Axisymmetric total curvature, including the meridional and hoop terms.
// A known regular pole uses its analytic limit, never an epsilon radius.
SurfaceSample sampleSurfaceRZ(const Arc&,double fraction,bool regularPole=false);
// Boundary-aware evaluation retains the declared common pole factors in r
// and z'. This resolves the whole neighboring arc, including samples much
// closer to a pole than direct power-basis subtraction can distinguish.
SurfaceSample sampleSurfaceRZ(const Arc&,double fraction,const CurveBoundary&);
// Endpoints, stationary principal curvatures, inflections and stationary
// parameter speed on the SAME quartic. Used for geometry observations and
// quadrature/source-chart partitions; it never fits or filters the curve.
// Roots use the existing numerical polynomial solver, not certified interval
// arithmetic. Declared RZ poles remove the common radial/axial-jet factor.
std::vector<double> surfacePartition(const Arc&,const CurveBoundary&,bool axisymmetric);
DistanceSample signedDistance(const Curve& effectiveBoundary,Point query);

template<std::size_t N>
std::vector<double> unitGridKnots(const PolynomialCurve<N>& curve,int cells) {
    std::vector<double> out;for(int i=0;i<=cells;++i)out.push_back(double(i)/cells);
    for(const auto& arc:curve.arcs)for(int axis=0;axis<2;++axis)for(int j=1;j<cells;++j) {
        auto p=arc.coordinate(axis);p[0]-=double(j)/cells;
        for(double t:roots(p,0,1))out.push_back(arc.a+(arc.b-arc.a)*t);
        // A shared join may straddle the face only through accumulated
        // coefficient roundoff, with no interior root in either piece.
        // Retain its existing node when the floating coordinate enclosure
        // meets the face. Do not move the node, merge roots or drop a piece.
        double scale=1;for(double coefficient:arc.coordinate(axis))scale+=std::abs(coefficient);
        const double coordinateBudget=64*std::numeric_limits<double>::epsilon()*scale;
        if(std::abs(polynomial(p,0))<=coordinateBudget)out.push_back(arc.a);
        if(std::abs(polynomial(p,1))<=coordinateBudget)out.push_back(arc.b);
    }
    std::sort(out.begin(),out.end());out.erase(std::unique(out.begin(),out.end()),out.end());return out;
}
} // namespace Quantum::Geometry::Parametric
