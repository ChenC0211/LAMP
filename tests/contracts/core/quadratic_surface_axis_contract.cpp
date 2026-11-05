#include "Core/Geometry/QuadraticSurfaceMeasure.h"
#include "Core/Mesh/CartesianMeasure.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

namespace
{
namespace S=LAMP::Geometry::QuadraticSurfaceMeasure;
namespace Q=LAMP::Geometry::QuadraticCellCut;
using Rectangle=LAMP::Geometry::PlanarCellCut::Rectangle;
unsigned assertions=0,guards=0;long double peak=0;
void check(bool ok,const char* why){++assertions;if(!ok)throw std::runtime_error(why);}
void near(long double a,long double b,long double scale,long double tolerance,const char* why)
{check(std::isfinite(a)&&std::isfinite(b)&&scale>0,"invalid comparison");const auto error=std::abs(a-b)/scale;peak=std::max(peak,error);check(error<=tolerance,why);}
template<class F> void rejects(F&& f,const char* why)
{bool caught=false;try{f();}catch(const std::exception&){caught=true;}check(caught,why);++guards;}
long double reference(const Rectangle& box,const Q::Reconstruction& graph,long double lo,long double hi,bool rz)
{
    // Independent positive composite GL5 on the graph formula. No clipping,
    // endpoint correction or call to SurfaceMeasure's integration is used.
    constexpr long double x[5]={-.90617984593866399280L,-.53846931010568309104L,0,.53846931010568309104L,.90617984593866399280L};
    constexpr long double w[5]={.23692688505618908751L,.47862867049936646804L,.56888888888888888889L,.47862867049936646804L,.23692688505618908751L};
    const auto basis=Q::basis(graph.shape);const long double cx=(static_cast<long double>(box.lower.x)+box.upper.x)/2;long double sum=0;
    for(int panel=0;panel<16;++panel){const auto a=lo+(hi-lo)*panel/16,b=lo+(hi-lo)*(panel+1)/16,mid=(a+b)/2,half=(b-a)/2;
        for(unsigned k=0;k<5;++k){const auto t=mid+half*x[k],height=graph.offset+graph.shape.slope*t+graph.shape.second*t*t/2;
            const auto r=cx+basis.tangent.x*t+basis.normal.x*height;check(!rz||r>=0,"all independent positive interior samples have physical radius");
            const auto ds=std::hypot(1.L,graph.shape.slope+graph.shape.second*t);sum+=half*w[k]*ds*(rz?CartesianMeasure::twoPi*r:1);}}
    return sum;
}
void rotatedAxis()
{
    const Rectangle box{{0,0},{1,1}};Q::Reconstruction graph;graph.shape={1,1,.123L,.37L,Vector2D{.6,.8}};graph.offset=-.1L;
    const auto xy=S::evaluate(box,graph,box,false),rz=S::evaluate(box,graph,box,true);check(xy.patches.size()==1&&rz.patches.size()==1,"one genuine rotated curved clipped arc");
    const auto& patch=rz.patches.front();check((patch.beginBoundary&S::Left)||(patch.endBoundary&S::Left),"clipping proves the physical r=0 endpoint");
    check(graph.shape.second!=0&&graph.shape.frameNormal.has_value(),"curved rotated frame regression, not a radial-coordinate line");
    near(rz.integral.length,reference(box,graph,patch.lower,patch.upper,false),1,1e-12L,"original physical length with tagged axis endpoint");
    near(rz.integral.physicalMeasure,reference(box,graph,patch.lower,patch.upper,true),1,1e-12L,"true revolved area from independent positive rule");
    near(xy.integral.length,rz.integral.length,1,1e-12L,"XY length unchanged by radial constraint metadata");
    rejects([&]{S::integrate(box,graph,(patch.beginBoundary&S::Left)?patch.lower-.1L:patch.lower,(patch.endBoundary&S::Left)?patch.upper+.1L:patch.upper,true);},"untagged public integration still rejects an actually negative-radius endpoint");
    const auto derivative=S::differentiate(box,graph,box,true);const long double h=1e-5L;
    for(unsigned j=0;j<3;++j){auto plus=graph,minus=graph;
        if(j==0){plus.offset+=h;minus.offset-=h;}else if(j==1){plus.shape.slope+=h;minus.shape.slope-=h;}else{plus.shape.second+=h;minus.shape.second-=h;}
        const auto numerical=(S::evaluate(box,plus,box,true).integral.physicalMeasure-S::evaluate(box,minus,box,true).integral.physicalMeasure)/(2*h);
        near(derivative.derivative[j],numerical,std::max(1.L,std::abs(numerical)),2e-8L,"axis-clipped physical endpoint derivative, finite-difference reference");}
    const unsigned rawOutside=static_cast<unsigned>(patch.begin.x<0)+static_cast<unsigned>(patch.end.x<0);
    if(rawOutside)rejects([&]{S::integrate(box,graph,patch.lower,patch.upper,true);},"raw no-tag API does not silently acquire clipping provenance");
    std::cout<<"ROTATED_AXIS,raw_negative_radius_endpoints="<<rawOutside<<",area="<<static_cast<double>(rz.integral.physicalMeasure)<<'\n';
}
void strictInteriorAndLine()
{
    const Rectangle box{{0,0},{1,1}};Q::Reconstruction bad;bad.shape={0,1,0,2,{}};bad.offset=-.51L;
    rejects([&]{S::integrate(box,bad,-.2L,.2L,true);},"negative interior stationary radius remains rejected despite positive endpoint radii");
    Q::Reconstruction line;line.shape={1,1,.2L,0,{}};line.offset=0;
    const auto measured=S::evaluate(box,line,box,true);near(measured.integral.physicalMeasure,CartesianMeasure::twoPi*.5L*std::hypot(1.L,.2L),1,1e-12L,"exact physical radial line integral");
    const auto d=S::differentiate(box,line,box,true);near(d.derivative[0],0,1,1e-12L,"radial clipping endpoints give zero offset derivative for this line");
    near(d.derivative[1],CartesianMeasure::twoPi*.5L*.2L/std::hypot(1.L,.2L),1,1e-12L,"line slope derivative keeps true radial measure");
}
}
int main()
{
    try{rotatedAxis();strictInteriorAndLine();std::cout<<"quadratic surface axis contract PASS assertions="<<assertions<<" guards="<<guards<<" peak="<<static_cast<double>(peak)<<'\n';return 0;}
    catch(const std::exception& error){std::cerr<<"quadratic surface axis contract FAIL "<<error.what()<<'\n';return 1;}
}