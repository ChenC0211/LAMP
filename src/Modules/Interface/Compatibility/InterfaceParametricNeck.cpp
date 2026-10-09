#include "Modules/Interface/Compatibility/InterfaceParametricNeck.h"
#include <numeric>

namespace InterfaceParametricProjection {
namespace {
using namespace Quantum::Geometry::Parametric;
void require(bool ok,const char* why){if(!ok)throw std::invalid_argument(why);}
double crossing(const Curve& c,double z){
    std::vector<double> found;
    for(const auto& a:c.arcs){
        auto p=a.coordinate(1);p[0]-=z;
        if(p[0]==0)found.push_back(a.a);
        for(double t:roots(p,0,1))found.push_back(a.a+(a.b-a.a)*t);
    }
    if(c.sample(1).y==z)found.push_back(1);
    std::sort(found.begin(),found.end());found.erase(std::unique(found.begin(),found.end()),found.end());
    require(found.size()==1,"neck section must have one transverse free-curve crossing");return found.front();
}
void append(Curve& result,const Curve& source,double begin,double end){
    for(const auto& a:source.arcs){const double lo=std::max(begin,a.a),hi=std::min(end,a.b);
        if(lo<hi)result.arcs.push_back(lo==a.a&&hi==a.b?a:restrictedArc(a,lo,hi));}
}
void normalize(Curve& c){
    const double start=c.arcs.front().a,length=c.arcs.back().b-start;
    require(length>0&&std::isfinite(length),"invalid child material interval");
    for(auto& a:c.arcs){a.a=(a.a-start)/length;a.b=(a.b-start)/length;}
    c.arcs.front().a=0;c.arcs.back().b=1;
}
double curvature(Point d,Point dd){return cross(d,dd)/std::pow(norm(d),3);}
// A single parametric quartic cap. Its endpoint arc-length gauges are
// |x_t(0)|=|x_t(1)|=L, x_tt(0)=L^2*kappa*N and r_tt(1)=0.
// The axis owns r(1)=0 and z_t(1)=0. These give
// (kappa*N_r)L^2/6 - (1-T_r)L + 2R = 0; take the branch continuous
// at kappa=0. The remaining axial coefficient is fixed by pi*int r^2 dz.
// No inverse r chart, tunable gap, optimization restart or alternate cap.
std::optional<Arc> cap(const Curve& source,double seam,double volume,bool lower,double section,double radialLimit,NeckSplitReport& report){
    // H2 curves are C1, so second derivatives at material knots are one-sided.
    // Match the retained parent arc, never the side that this event removes.
    const auto retained=lower?std::lower_bound(source.arcs.begin(),source.arcs.end(),seam,[](const Arc& a,double s){return a.b<s;}):source.arcs.begin()+source.index(seam);
    require(retained!=source.arcs.end(),"neck seam has no retained parent arc");
    const auto p=retained->sample(seam),d=retained->sample(seam,1),dd=retained->sample(seam,2);const double speed=norm(d);
    require(p.x>0&&d.y>0&&speed>0&&std::isfinite(speed),"neck seam requires a regular transverse axial graph");
    const int side=lower?0:1;const double sign=lower?1:-1;
    const Point tangent{sign*d.x/speed,d.y/speed};const double kappa=curvature(d,dd),radialSecond=-kappa*tangent.y;
    const double b=1-tangent.x,discriminant=b*b-4*radialSecond*p.x/3;
    report.seamRadialTangent[side]=tangent.x;report.seamSignedCurvature[side]=kappa;report.capDiscriminant[side]=discriminant;
    require(std::isfinite(kappa)&&std::isfinite(discriminant),"nonfinite parametric neck data");
    if(discriminant<0)return {};
    const double length=4*p.x/(b+std::sqrt(discriminant)),rr=radialSecond*length*length;
    require(length>0&&std::isfinite(length),"invalid parametric neck length");report.capLength[side]=length;
    Arc local;local.a=0;local.b=1;local.p[0]={p.x,0};local.p[1]=tangent*length;
    local.p[2]={.5*rr,.5*length*length*kappa*tangent.x};
    local.p[3].x=-2*p.x-2*local.p[1].x-5*rr/6;local.p[4].x=p.x+local.p[1].x+rr/3;
    local.p[4].y=-(local.p[1].y+2*local.p[2].y)/4;
    // The free axial basis is t^3-3t^4/4; its derivative is nonnegative.
    // Eight-point Gauss integrates the degree-11 volume exactly, up to roundoff.
    double base=0,coefficient=0;const auto& rule=gaussEight();
    for(int i=0;i<8;i++){const double t=.5*(1+rule.x[i]),r=local.sampleFraction(t).x,weight=.5*rule.w[i]*r*r;
        base+=weight*local.sampleFraction(t,1).y;coefficient+=weight*3*t*t*(1-t);}
    require(coefficient>0&&std::isfinite(base)&&std::isfinite(coefficient),"invalid parametric cap volume functional");
    const double free=(volume/std::acos(-1.)-base)/coefficient;local.p[3].y=free;local.p[4].y-=.75*free;
    const double height=.75*local.p[1].y+.5*local.p[2].y+.25*free,pole=p.y+sign*height;
    if(lower)report.lowerPole=pole;else report.upperPole=pole;
    // Factor out the exact axis endpoint before range evaluation: r=(1-t)q,
    // z_t=(1-t)w. This tests regular interior geometry, not a rounded endpoint.
    const std::array<double,4> radial{{p.x,p.x+local.p[1].x,p.x+local.p[1].x+local.p[2].x,-local.p[4].x}};
    const std::array<double,3> axial{{local.p[1].y,local.p[1].y+2*local.p[2].y,local.p[1].y+2*local.p[2].y+3*free}};
    require(std::isfinite(height)&&std::isfinite(free),"nonfinite constrained parametric neck cap");
    report.capAdmissible[side]=height>0&&(lower?pole<section:pole>section)&&range(radial)[0]>0&&range(axial)[0]>0&&range(local.coordinate(0))[1]<radialLimit;
    if(!report.capAdmissible[side])return {};
    for(auto& q:local.p)q.y*=sign;local.p[0].y+=p.y;
    Arc result=lower?local:restrictedArc(local,1,0);const double materialLength=length/speed;
    result.a=lower?seam:seam-materialLength;result.b=lower?seam+materialLength:seam;
    const auto cp=result.sampleFraction(lower?0:1),cd=result.sampleFraction(lower?0:1,1),cdd=result.sampleFraction(lower?0:1,2);
    report.seamPosition=std::max(report.seamPosition,norm(cp-p));report.seamTangent=std::max(report.seamTangent,norm(cd-d));
    report.seamCurvature=std::max(report.seamCurvature,std::abs(curvature(cd,cdd)-kappa));
    return result;
}
}
std::vector<NeckSite> inspectNecks(const Candidate& source,double limit){
    const auto& grid=source.grid();require(grid.axisymmetric&&grid.left==0&&std::isfinite(limit)&&limit>0&&limit<=.5,"invalid RZ neck detection scope");
    struct Section {double z,h;};std::vector<Section> sections;
    std::vector<Grid::Rectangle> axisCells;
    for(std::size_t k=0;k<grid.size();k++){const auto a=grid.rectangle(k);if(a.left==0)axisCells.push_back(a);}
    std::sort(axisCells.begin(),axisCells.end(),[](const auto& a,const auto& b){return a.bottom<b.bottom;});
    for(std::size_t k=1;k<axisCells.size();k++){const auto a=axisCells[k-1],b=axisCells[k];if(a.top==b.bottom)
        sections.push_back({a.top,std::min({a.right,b.right,a.top-a.bottom,b.top-b.bottom})});}
    std::vector<NeckSite> result;
    for(const auto& component:source.boundaryHandle()->components()){
        const auto& curve=*component.curve;if(!curve.boundary.axis(1)||(!curve.boundary.axis(0)&&!curve.boundary.fixed(0)))continue;
        for(const auto& arc:curve.arcs){auto coordinate=arc.coordinate(0);std::array<double,4> derivative{};for(int k=1;k<5;k++)derivative[k-1]=k*coordinate[k];
            auto minima=roots(derivative,0,1);if(polynomial(derivative,0)==0)minima.push_back(0);
            for(double t:minima){const double s=arc.a+(arc.b-arc.a)*t;const auto p=arc.sampleFraction(t),d=arc.sampleFraction(t,1),dd=arc.sampleFraction(t,2);
                if(!(s>0&&s<1&&p.x>0&&d.y>0&&dd.x>0)||sections.empty())continue;
                const auto nearest=std::min_element(sections.begin(),sections.end(),[&](const Section& a,const Section& b){return std::abs(a.z-p.y)<std::abs(b.z-p.y);});
                bool repeated=false;for(auto& old:result)if(old.options.parent==component.id&&old.options.section==nearest->z){
                    if(p.x<old.radius){old.radius=p.x;old.position=p.y;}repeated=true;}
                if(!repeated)result.push_back({{component.id,nearest->z,limit,1e-10},p.x,p.y,nearest->h});
            }
        }
    }
    return result;
}
std::vector<NeckSplitOptions> detectNecks(const Candidate& source,double limit){
    std::vector<NeckSplitOptions> result;for(const auto& site:inspectNecks(source,limit))if(site.radius<=limit*site.localScale)result.push_back(site.options);return result;
}
NeckAssessment assessNeckSplit(const Candidate& source,const NeckSplitOptions& options){
    const auto& grid=source.grid();
    require(source.componentCount()>0&&grid.axisymmetric&&grid.left==0,"neck split requires an owned RZ boundary on the axis");
    require(std::isfinite(options.section)&&std::isfinite(options.maximumRadiusCells)&&options.maximumRadiusCells>0&&options.maximumRadiusCells<=.5
        &&std::isfinite(options.inventoryFractionTolerance)&&options.inventoryFractionTolerance>0,"invalid explicit neck model scale or budget");
    std::size_t parent=source.componentCount();std::uint64_t largest=0;
    for(const auto& c:source.retiredComponents())largest=std::max(largest,c.id);
    for(std::size_t i=0;i<source.componentCount();i++){const auto id=source.boundaryHandle()->components()[i].id;largest=std::max(largest,id);if(id==options.parent)parent=i;}
    require(parent<source.componentCount()&&largest<std::numeric_limits<std::uint64_t>::max()-1,"invalid neck parent or exhausted material identity");
    const auto old=source.component(parent);const auto& curve=old.geometry();
    require(curve.boundary.axis(1)&&(curve.boundary.axis(0)||curve.boundary.fixed(0)),"neck split requires an axis-ending meridian");
    NeckAssessment assessment;auto& report=assessment.report;report.parent=options.parent;report.children={largest+1,largest+2};report.section=options.section;
    for(std::size_t k=0;k<grid.size();k++){const auto cell=grid.rectangle(k);if(cell.left!=0)continue;
        if(cell.top==options.section){require(report.cells[0]<0,"ambiguous lower neck cell");report.cells[0]=int(k);}
        if(cell.bottom==options.section){require(report.cells[1]<0,"ambiguous upper neck cell");report.cells[1]=int(k);}}
    require(report.cells[0]>=0&&report.cells[1]>=0,"neck section must be a shared axis-cell plane");
    const auto lower=grid.rectangle(report.cells[0]),upper=grid.rectangle(report.cells[1]);
    report.lowerSection=lower.bottom;report.upperSection=upper.top;
    report.localScale=std::min({lower.right,upper.right,lower.top-lower.bottom,upper.top-upper.bottom});
    const double s0=crossing(curve,lower.bottom),sn=crossing(curve,options.section),s1=crossing(curve,upper.top);
    require(0<s0&&s0<sn&&sn<s1&&s1<1,"neck patch must lie in the interior of one meridian");
    report.neckRadius=std::numeric_limits<double>::infinity();report.sectionRadius=curve.sample(sn).x;
    Curve patch;append(patch,curve,s0,s1);
    for(const auto& a:patch.arcs){
        const auto radial=range(a.coordinate(0));std::array<double,4> dz{};for(int k=1;k<5;k++)dz[k-1]=k*a.p[k].y;
        require(range(dz)[0]>0&&radial[0]>0&&radial[1]<std::min(lower.right,upper.right),"source neck patch must be an interior axial graph in the two axis cells");
        auto x=a.coordinate(0);std::array<double,4> dx{};for(int k=1;k<5;k++)dx[k-1]=k*x[k];auto extrema=roots(dx,0,1);extrema.push_back(0);extrema.push_back(1);
        for(double t:extrema){const auto p=a.sampleFraction(t);if(p.x<report.neckRadius){report.neckRadius=p.x;report.neckPosition=p.y;}}
    }
    require(report.neckPosition>lower.bottom&&report.neckPosition<upper.top,"neck minimum is outside the selected axial patch");
    if(report.sectionRadius>options.maximumRadiusCells*report.localScale)return assessment;
    std::array<std::vector<double>,2> inventories;for(auto& v:inventories)v.assign(grid.size(),0);
    for(std::size_t k=0;k<grid.size();k++){
        const auto c=grid.rectangle(k);const double v=old.inventory()[k];
        require(std::isfinite(v),"invalid parent neck inventory");
        if(c.top<=options.section)inventories[0][k]=v;else if(c.bottom>=options.section)inventories[1][k]=v;
        else require(v==0,"occupied coarse cell straddles the event section; refine the event support first");
        report.parentVolume+=v;
    }
    for(int k=0;k<2;k++){const auto slot=report.cells[k];require(old.inventory()[slot]>0&&old.inventory()[slot]<grid.capacity(slot),"neck cell requires strictly mixed inventory");report.modelRegionVolume+=old.inventory()[slot];}
    Curve first,second;first.boundary=curve.boundary;first.boundary.kind[1]=CurveBoundary::Kind::Axis;
    second.boundary=CurveBoundary::axisMeridian();
    const auto lowerCap=cap(curve,s0,old.inventory()[report.cells[0]],true,options.section,lower.right,report);
    const auto upperCap=cap(curve,s1,old.inventory()[report.cells[1]],false,options.section,upper.right,report);
    if(!lowerCap||!upperCap){assessment.readiness=NeckReadiness::CapPending;return assessment;}
    append(first,curve,0,s0);first.arcs.push_back(*lowerCap);
    second.arcs.push_back(*upperCap);append(second,curve,s1,1);
    normalize(first);normalize(second);
    report.retainedParentIntervals={{{0,s0},{s1,1}}};
    report.capIntervals={{{first.arcs.back().a,1},{0,second.arcs.front().b}}};
    std::vector<Candidate> result;std::vector<std::uint64_t> ids;
    for(std::size_t i=0;i<source.componentCount();i++){
        if(i!=parent){result.push_back(source.component(i));ids.push_back(source.boundaryHandle()->components()[i].id);continue;}
        std::array<Curve,2> children{{std::move(first),std::move(second)}};
        for(int child=0;child<2;child++){
            // initialize validates the effective boundary and evaluates its
            // actual geometry; it is not the new inventory's source.
            auto candidate=initialize(std::move(children[child]),grid);double totalDefect=0;
            for(std::size_t k=0;k<grid.size();k++){
                const double difference=candidate.inventory()[k]-inventories[child][k],fraction=std::abs(difference)/grid.capacity(k);
                require(fraction<=options.inventoryFractionTolerance,"neck caps do not satisfy the independent local inventory");
                candidate.report_.maximumInventoryFractionError=std::max(candidate.report_.maximumInventoryFractionError,fraction);totalDefect+=difference;
            }
            candidate.inventory_=std::move(inventories[child]);candidate.report_.totalGeometryInventoryDifference=totalDefect;
            report.childVolumes[child]=std::accumulate(candidate.inventory_.begin(),candidate.inventory_.end(),0.);
            report.maximumInventoryFractionError=std::max(report.maximumInventoryFractionError,candidate.report_.maximumInventoryFractionError);
            result.push_back(std::move(candidate));ids.push_back(report.children[child]);
        }
    }
    auto candidate=combine(grid,std::move(result),ids,source.retiredComponents());
    for(std::size_t k=0;k<grid.size();k++){report.maximumCellBalance=std::max(report.maximumCellBalance,std::abs(candidate.inventory()[k]-source.inventory()[k]));
        require(std::abs(candidate.inventory()[k]-source.inventory()[k])<=8*std::numeric_limits<double>::epsilon()*grid.capacity(k),"neck event changed a native cell inventory");}
    assessment.readiness=NeckReadiness::Ready;assessment.split=NeckSplit{std::move(candidate),report};return assessment;
}
NeckSplit prepareNeckSplit(const Candidate& source,const NeckSplitOptions& options){
    auto assessment=assessNeckSplit(source,options);require(bool(assessment.split),"neck event has not entered the joint radius/cap feasible set");return std::move(*assessment.split);
}
}
