#include "Core/Geometry/QuadraticCellCut.h"
#include "Core/Geometry/QuadraticSurfaceMeasure.h"
#include "Core/Geometry/QuadraticCurveMeasure.h"
#include "Core/Mesh/CartesianMeasure.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
namespace Q=LAMP::Geometry::QuadraticCellCut;
namespace S=LAMP::Geometry::QuadraticSurfaceMeasure;
namespace C=LAMP::Geometry::QuadraticCurveMeasure;
void check(bool ok,const char*message){if(!ok)throw std::runtime_error(message);}
void near(long double a,long double b,long double scale,const char*message){check(std::isfinite(a)&&std::isfinite(b)&&scale>0&&std::abs(a-b)<=1e-12L*scale,message);}
template<class F>void rejects(F&&f,const char*message){bool bad=false;try{f();}catch(const std::exception&){bad=true;}check(bad,message);}
constexpr long double pi=static_cast<long double>(CartesianMeasure::twoPi)/2;
const Q::Rectangle box{{0,.25},{.25,.5}};const Q::Point axis{0,.34375};
void inspect(const Q::AxisPortVolumeResult&solved,int sign,Q::Point first,Q::Point last,long double target,bool horizontal){
    check(solved.graph.exact!=nullptr&&solved.graph.shape.sign==sign,"same exact axis polynomial and phase retained");
    check(solved.horizontalCompanion==horizontal&&solved.companionFraction>0&&solved.companionFraction<1,"actual open companion segment report");
    const auto represented=horizontal?first.x+(last.x-first.x)*solved.companionFraction:first.y+(last.y-first.y)*solved.companionFraction;
    near(horizontal?solved.companion.x:solved.companion.y,represented,.25L,"reported fraction follows original partial/reversed segment");
    check((horizontal?solved.companion.y:solved.companion.x)==(horizontal?first.y:first.x),"companion actual physical face identity");
    const auto volume=Q::graphMeasure(box,solved.graph,true).volume;const long double capacity=pi*.25L*.25L*.25L,own=std::min(target,capacity-target);
    near(volume,target,own,"actual own-phase Core volume");
    auto opposite=solved.graph;opposite.shape.sign=-sign;const auto other=Q::graphMeasure(box,opposite,true).volume;
    near(other,capacity-target,capacity-target,"independent opposite phase own-scale certificate");
    near(solved.phaseVolume,target,own,"reported actual phase volume");check(solved.ownPhaseResidual<=1e-12L,"axis volume original own-phase certificate");
    near(solved.volumeResidual,volume-target,target,"reported true signed volume residual");
    near(Q::implicitValue(box,solved.graph,axis),0,.25L,"same polynomial passes fixed axis point");
    const auto surface=S::evaluate(box,solved.graph,box,true);check(surface.patches.size()==1,"one actual regular axis arc");const auto&arc=surface.patches.front();
    check(arc.begin.x==axis.x&&arc.begin.y==axis.y&&arc.end.x==solved.companion.x&&arc.end.y==solved.companion.y,"same source root events preserve axis/companion bit identity");
    const auto sample=S::sample(box,solved.graph,arc.lower,true);check(sample.position.x==0&&sample.unitNormal.x==0&&std::isfinite(sample.curvature),"true regular axis derivative/curvature limit");
    check(surface.integral.physicalMeasure>0&&std::isfinite(surface.integral.physicalMeasure),"true revolved surface measure");
    const auto material=C::Curve::fromGraph(Q::materialCurveData(box,solved.graph,arc.lower,arc.upper));const auto begin=C::position(material,0),end=C::position(material,1);
    check(begin.x==static_cast<double>(axis.x)&&begin.y==static_cast<double>(axis.y)&&end.x==static_cast<double>(solved.companion.x)&&end.y==static_cast<double>(solved.companion.y),"one exact source material restriction; no native sqrt refit");
}
void right(){
    const Q::Point lower{.25,.375},upper{.25,.4375};const long double z=.3984375L;
    const long double liquid=pi*.25L*.25L*((axis.y+z)/2-.25L),gas=pi*.25L*.25L*(.5L-(axis.y+z)/2);
    for(int sign:{1,-1})for(bool reverse:{false,true}){const auto first=reverse?upper:lower,last=reverse?lower:upper;const auto target=sign>0?liquid:gas;
        const auto solved=Q::axisPortVolumeGraph(box,sign,axis,first,last,target);inspect(solved,sign,first,last,target,false);near(solved.companion.y,z,.25L,"right face direct analytic companion height");}
}
void horizontal(bool top,bool tiny){
    const auto z=top?.5L:.25L,dz=std::abs(z-axis.y);const long double radiusSquared=.25L*.25L*(tiny?std::ldexp(1.L,-30):1.L/3);
    // Exact RZ wedge: small stock = pi*|z_face-z_axis|*r_B^2/2.
    // The ordinary case requires a non-rational native positive sqrt view.
    const long double small=pi*dz*radiusSquared/2,capacity=pi*.25L*.25L*.25L;
    const Q::Point lower{tiny?0.L:.0625L,z},upper{tiny?.0625L:.1875L,z};
    for(int sign:{1,-1}){if(tiny&&sign!=(top?-1:1))continue;const long double liquid=top?capacity-small:small,gas=top?small:capacity-small,target=sign>0?liquid:gas;
        for(bool reverse:{false,true}){const auto first=reverse?upper:lower,last=reverse?lower:upper;const auto solved=Q::axisPortVolumeGraph(box,sign,axis,first,last,target);inspect(solved,sign,first,last,target,true);
            near(solved.companion.x*solved.companion.x,radiusSquared,radiusSquared,"analytic r-squared retained through root view");}}
}
void invalid(){
    const auto target=pi*.25L*.25L*.125L;const Q::Point lower{.25,.375},upper{.25,.4375};
    rejects([&]{(void)Q::axisPortVolumeGraph(box,0,axis,lower,upper,target);},"invalid phase sign accepted");
    rejects([&]{(void)Q::axisPortVolumeGraph(box,1,{.01,.34375},lower,upper,target);},"nonaxis fixed point accepted");
    rejects([&]{(void)Q::axisPortVolumeGraph(box,1,{0,.25},lower,upper,target);},"axis corner accepted as regular open endpoint");
    rejects([&]{(void)Q::axisPortVolumeGraph(box,1,axis,lower,lower,target);},"zero companion segment accepted");
    rejects([&]{(void)Q::axisPortVolumeGraph(box,1,axis,{.1,.375},{.1,.4375},target);},"interior segment masquerades as physical face");
    rejects([&]{(void)Q::axisPortVolumeGraph(box,1,axis,{.25,.375},{.1875,.4375},target);},"nonaxis-aligned companion segment accepted");
    rejects([&]{(void)Q::axisPortVolumeGraph(box,1,axis,{0,.375},{0,.4375},target);},"second radial-axis port accepted as companion");
    rejects([&]{(void)Q::axisPortVolumeGraph(box,1,axis,lower,upper,0);},"zero phase target accepted");
    rejects([&]{(void)Q::axisPortVolumeGraph(box,1,axis,lower,upper,std::numeric_limits<long double>::quiet_NaN());},"nonfinite target accepted");
    const auto outside=pi*.25L*.25L*((axis.y+.48L)/2-.25L);
    rejects([&]{(void)Q::axisPortVolumeGraph(box,1,axis,lower,upper,outside);},"solved companion outside actual partial face accepted");
    rejects([&]{(void)Q::axisPortVolumeGraph({{.01,.25},{.25,.5}},1,axis,lower,upper,target);},"nonaxis cell accepted by regular-axis chart");
}
int main(){try{right();horizontal(true,false);horizontal(false,false);horizontal(true,true);horizontal(false,true);invalid();
    std::cout<<"quadratic axis volume graph contract PASS; right/top/bottom, two phases, partial/reversed ports and own-scale small stocks\n";return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}

