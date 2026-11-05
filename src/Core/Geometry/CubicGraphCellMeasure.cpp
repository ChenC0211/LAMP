#include "Core/Geometry/CubicGraphCellMeasure.h"
#include "Core/Geometry/UnitPolynomialRoots.h"
#include "Core/Mesh/CartesianMeasure.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace LAMP::Geometry::CubicGraphCellMeasure
{
namespace
{
using Poly=std::array<long double,7>;
using Point=PlanarCellCut::Point;
void require(bool ok,const char* message) {if(!ok) throw std::invalid_argument(message);}
long double value(const Poly& p,long double x)
{long double v=p[6];for(int k=5;k>=0;--k) v=v*x+p[k];return v;}
long double integral(const Poly& p,long double lo,long double hi)
{
    const auto m=lo+(hi-lo)/2,h=(hi-lo)/2,m2=m*m,h2=h*h;
    const auto q2=p[2]+3*p[3]*m+6*p[4]*m2+10*p[5]*m2*m+15*p[6]*m2*m2;
    const auto q4=p[4]+5*p[5]*m+15*p[6]*m2;
    return 2*h*(value(p,m)+q2*h2/3+q4*h2*h2/5+p[6]*h2*h2*h2/7);
}
}
Result evaluate(const PlanarCellCut::Rectangle& box,Point origin,Vector2D normal,long double scale,
    const std::array<long double,4>& c,bool rz)
{
    const auto capacity=PlanarCellCut::capacity(box,rz);
    const auto length=std::hypot(static_cast<long double>(normal.x),static_cast<long double>(normal.y));
    require(std::isfinite(scale) && scale>0 && std::isfinite(length) && length>0
        && std::isfinite(origin.x) && std::isfinite(origin.y),"[CubicGraphCellMeasure] invalid frame");
    for(auto v:c) require(std::isfinite(v),"[CubicGraphCellMeasure] nonfinite coefficient");
    const Point n{normal.x/length,normal.y/length},t{-n.y,n.x};
    const std::array<Point,4> corners{{{box.lower.x,box.lower.y},{box.upper.x,box.lower.y},
        {box.upper.x,box.upper.y},{box.lower.x,box.upper.y}}};
    std::array<Point,4> local;std::vector<long double> slabs;
    for(int k=0;k<4;++k) {
        const auto x=(corners[k].x-origin.x)/scale,y=(corners[k].y-origin.y)/scale;
        local[k]={x*t.x+y*t.y,x*n.x+y*n.y};slabs.push_back(local[k].x);
        require(std::isfinite(local[k].x) && std::isfinite(local[k].y),"[CubicGraphCellMeasure] unrepresentable local frame");
    }
    std::sort(slabs.begin(),slabs.end());slabs.erase(std::unique(slabs.begin(),slabs.end()),slabs.end());
    require(slabs.size()>=2,"[CubicGraphCellMeasure] collapsed local frame");
    Result result;bool full=true,empty=true,visited=false;
    for(std::size_t k=1;k<slabs.size();++k) {
        const auto lo=slabs[k-1],hi=slabs[k],w=hi-lo;
        const Poly graph{{((c[3]*lo+c[2])*lo+c[1])*lo+c[0],
            (c[1]+2*c[2]*lo+3*c[3]*lo*lo)*w,(c[2]+3*c[3]*lo)*w*w,c[3]*w*w*w,0,0,0}};
        std::vector<Poly> edges;
        for(int i=0;i<4;++i) {
            const auto a=local[i],b=local[(i+1)%4];
            if(a.x==b.x || lo<std::min(a.x,b.x) || hi>std::max(a.x,b.x)) continue;
            edges.push_back({{a.y+(lo-a.x)/(b.x-a.x)*(b.y-a.y),w/(b.x-a.x)*(b.y-a.y),0,0,0,0,0}});
        }
        if(edges.empty()) continue;
        require(edges.size()==2,"[CubicGraphCellMeasure] invalid cell slice");
        visited=true;
        if(value(edges[0],.5L)>value(edges[1],.5L)) std::swap(edges[0],edges[1]);
        std::vector<UnitPolynomialRoots::Root> cuts{{0,0,0},{1,1,1}};
        for(const auto& edge:edges) {
            std::vector<UnitPolynomialTerm> terms;
            for(int j=0;j<4;++j) if(graph[j]!=0) {
                UnitPolynomialTerm term;term.factor=graph[j];term.coefficients[j]=1;terms.push_back(term);
            }
            for(int j=0;j<2;++j) if(edge[j]!=0) {
                UnitPolynomialTerm term;term.factor=-edge[j];term.coefficients[j]=1;terms.push_back(term);
            }
            const auto roots=isolateUnitPolynomial(terms);
            for(const auto& root:roots.roots) if(root.point>0 && root.point<1) cuts.push_back(root);
        }
        std::sort(cuts.begin(),cuts.end(),[](const auto& a,const auto& b){return a.point<b.point;});
        for(std::size_t j=1;j<cuts.size();++j) {
            const auto& left=cuts[j-1];const auto& right=cuts[j];
            require(left.upper<right.lower,"[CubicGraphCellMeasure] unresolved slice root ordering");
            const long double mid=left.upper+(right.lower-left.upper)/2,a=left.point,b=right.point;
            require(mid>left.upper && mid<right.lower,"[CubicGraphCellMeasure] no representable phase sample between root enclosures");
            const auto lower=edges[0];auto upper=edges[1];
            if(value(graph,mid)<=value(lower,mid)) {full=false;continue;}
            empty=false;const bool cut=value(graph,mid)<value(upper,mid);
            if(cut) {upper=graph;full=false;}
            Poly difference{},density{};for(int i=0;i<7;++i) difference[i]=upper[i]-lower[i];
            if(rz) {
                Poly meanRadius{};for(int i=0;i<7;++i) meanRadius[i]=scale*n.x*(upper[i]+lower[i])/2;
                meanRadius[0]+=origin.x+scale*t.x*lo;meanRadius[1]+=scale*t.x*w;
                for(int i=0;i<7;++i) for(int m=0;i+m<7;++m)
                    density[i+m]+=CartesianMeasure::twoPi*difference[i]*meanRadius[m];
            } else density=difference;
            result.volume+=scale*scale*w*integral(density,a,b);
            if(cut) {
                // Boundary terms cancel where the two slice formulas meet.
                Poly radius{};radius[0]=1;
                if(rz) {
                    radius={};for(int i=0;i<4;++i) radius[i]=scale*n.x*graph[i];
                    radius[0]+=origin.x+scale*t.x*lo;radius[1]+=scale*t.x*w;
                    for(auto& v:radius) v*=CartesianMeasure::twoPi;
                }
                Poly basis{};basis[0]=1;
                for(int m=0;m<4;++m) {
                    Poly sensitivity{};
                    for(int i=0;i<7;++i) for(int p=0;i+p<7;++p) sensitivity[i+p]+=radius[i]*basis[p];
                    result.derivative[m]+=scale*scale*w*integral(sensitivity,a,b);
                    Poly next{};for(int p=0;p<6;++p) {next[p]+=lo*basis[p];next[p+1]+=w*basis[p];}basis=next;
                }
            }
        }
    }
    require(visited,"[CubicGraphCellMeasure] no represented cell slice");
    if(full) result.volume=capacity;else if(empty) result.volume=0;
    const auto allowance=PlanarCellCut::relativeVolumeResolution*capacity;
    require(std::isfinite(result.volume) && result.volume>=-allowance && result.volume<=capacity+allowance,
        "[CubicGraphCellMeasure] invalid volume");
    for(auto d:result.derivative) require(std::isfinite(d),"[CubicGraphCellMeasure] nonfinite sensitivity");
    return result;
}
}
