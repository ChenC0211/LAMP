#include "Modules/Hydrodynamics/HydroBalancedCsfLoad.h"
#include "interface_first_jet_fixture.h"
#include <cstdint>
#include <cstring>
#include <map>
#include <numeric>

namespace
{
using namespace LAMP::Hydrodynamics;
using Box=std::array<Vector2D,2>;
int scenes=0,guards=0,amrWitness=0,axisWitness=0,visibilityWitness=0;
double peakMean=0,peakLoad=0;
template<class F>void rejection(F f,const char* expected){rejects(f,expected);++guards;}
double rel(double error,double scale){check(std::isfinite(error)&&std::isfinite(scale)&&scale>=0,"nonfinite reference scale");
    if(scale==0){check(error==0,"nonzero error with zero scale");return 0;}return std::abs(error)/scale;}
struct Shape
{
    std::string kind="circle";bool rz=false;double cx=.43,cy=.51,a=.19,b=.19,amplitude=.015;
    Profile profile()const
    {
        const auto s=*this;return [s](Vector2D p){if(s.kind=="sinus"){
            const double k=double(2*pi),angle=k*p.x;
            return Jet{p.y-.5-s.amplitude*std::cos(angle),(p.x==0||p.x==1)?0:s.amplitude*k*std::sin(angle),1,0};}
            const double x=p.x-(s.rz?0:s.cx),y=p.y-s.cy;
            return Jet{x*x/(s.a*s.a)+y*y/(s.b*s.b)-1,2*x/(s.a*s.a),2*y/(s.b*s.b),0};};
    }
    std::vector<LevelSetInterfacePanels::SymmetryPlane> planes()const
    {return kind=="sinus"?(rz?std::vector<LevelSetInterfacePanels::SymmetryPlane>{{0,1}}:
        std::vector<LevelSetInterfacePanels::SymmetryPlane>{{0,0},{0,1}}):std::vector<LevelSetInterfacePanels::SymmetryPlane>{};}
    struct Arc{Vector2D point;long double curvature=0,weight=0;};
    Arc arc(long double t)const
    {
        if(kind=="sinus"){const auto k=2*pi,angle=k*t,slope=-amplitude*k*std::sin(angle),speed=std::sqrt(1+slope*slope);
            const auto meridional=amplitude*k*k*std::cos(angle)/(speed*speed*speed);
            const auto azimuthal=t==0?amplitude*k*k:amplitude*k*std::sin(angle)/(t*speed);
            return {{static_cast<double>(t),static_cast<double>(.5L+amplitude*std::cos(angle))},
                meridional+(rz?azimuthal:0),speed*(rz?2*pi*t:1)};}
        if(rz){const auto r=a*std::sin(t),z=cy+b*std::cos(t),speed=std::hypot(a*std::cos(t),b*std::sin(t));
            return {{static_cast<double>(r),static_cast<double>(z)},a*b/(speed*speed*speed)+b/(a*speed),2*pi*r*speed};}
        const auto speed=std::hypot(a*std::sin(t),b*std::cos(t));
        return {{static_cast<double>(cx+a*std::cos(t)),static_cast<double>(cy+b*std::sin(t))},a*b/(speed*speed*speed),speed};
    }
    long double end()const{return kind=="sinus"?1:rz?pi:2*pi;}
    long double sharpCurvature(Vector2D p)const
    {
        // Test-only closest-point search on the analytic continuous curve.
        // This value is never passed to the capillary factory.
        constexpr int count=256;long double best=0,value=std::numeric_limits<long double>::infinity();
        const auto distance=[&](long double t){const auto q=arc(t).point;const long double x=q.x-p.x,y=q.y-p.y;return x*x+y*y;};
        for(int i=0;i<=count;++i){const auto t=end()*i/count,v=distance(t);if(v<value){best=t;value=v;}}
        const bool periodic=kind!="sinus"&&!rz;
        long double lo=periodic?best-end()/count:std::max(0.L,best-end()/count),
            hi=periodic?best+end()/count:std::min(end(),best+end()/count);
        for(int i=0;i<64;++i){const auto left=(2*lo+hi)/3,right=(lo+2*hi)/3;if(distance(left)<distance(right))hi=right;else lo=left;}
        return arc((lo+hi)/2).curvature;
    }
};
Shape shapeFor(std::string kind,bool rz)
{
    check(kind=="circle"||kind=="ellipse"||kind=="sinus","invalid BF03 shape");Shape s;s.kind=std::move(kind);s.rz=rz;
    if(s.kind=="ellipse"){s.a=rz?.31:.27;s.b=rz?.23:.18;s.cx=.46;s.cy=.52;}
    else if(s.kind=="circle"){s.a=s.b=rz?.34:.19;}return s;
}
HermiteGalsPublication publish(const QuadtreeMesh& mesh,const Profile& f,bool rz,
    const std::vector<LevelSetInterfacePanels::SymmetryPlane>& planes={},HermiteGalsPublicationConfig cfg={})
{
    const auto seed=stateFor(mesh,f,rz,{},[](int,Vector2D,Vector2D){return 0.;},planes,HermiteGalsTransport::FirstJet);
    HermiteGalsCandidate candidate(mesh,{},seed.primary());const auto mixed=closeHermiteGalsMixedAmr(mesh,candidate,planes);
    cfg.volume.symmetryPlanes=planes;const auto id=candidate.identity();
    return publishHermiteGalsFirstJetGeometry(mesh,std::move(candidate),mixed,{id,centers(mesh,f)},cfg);
}
HydroPhysicalDivergence metrics(const QuadtreeMesh& mesh,bool rz,int selection=0)
{
    FaceSegmentGraph graph;graph.build(mesh,[selection](const QuadtreeMesh&,const Cell& c){const auto p=c.centerCoord();
        return selection==0||(selection==1?!(p.x>.375&&p.x<.625&&p.y>.375&&p.y<.625):selection==2?p.x<.5:p.x<.375||p.x>.625);});
    return HydroPhysicalDivergence::build(mesh,HydroIntegralDivergence::build(mesh,
        HydroDofLayout::build(mesh,NativeLeafDomain::build(mesh,graph,rz))));
}
HydroRegularizedColor colorFor(const QuadtreeMesh& mesh,const HydroPhaseGeometry& phase,double epsilon)
{HydroRegularizedColorConfig cfg;cfg.epsilon=epsilon;return HydroRegularizedColor::build(mesh,phase,cfg);}
struct Reference
{
    std::vector<Box> boxes;std::vector<int> component,minus,plus;std::vector<double> area,diagonal;
    std::vector<Vector2D> centers;int components=0;bool rectangle=false;
    Reference(const QuadtreeMesh& mesh,const HydroPhaseGeometry& phase)
    {
        const auto& l=phase.layout();component.resize(l.pressureCount());std::iota(component.begin(),component.end(),0);
        const auto root=[&](int p){while(component[p]!=p)p=component[p];return p;};
        long double areaSum=0;Vector2D lower{1,1},upper{0,0};
        for(int p=0;p<l.pressureCount();++p){const auto& c=mesh.cells(mesh.cellIndexToLabel[l.nativeCell(p)]);
            const auto a=mesh.verts(c.cellToVert[0]).coord(),b=mesh.verts(c.cellToVert[2]).coord();boxes.push_back({a,b});
            diagonal.push_back(std::hypot(b.x-a.x,b.y-a.y));areaSum+=(b.x-a.x)*(b.y-a.y);
            lower.x=std::min(lower.x,a.x);lower.y=std::min(lower.y,a.y);upper.x=std::max(upper.x,b.x);upper.y=std::max(upper.y,b.y);}
        rectangle=areaSum==static_cast<long double>(upper.x-lower.x)*(upper.y-lower.y);
        for(int u=0;u<l.velocityCount();++u){const auto f=l.domain().face(mesh,l.nativeFace(u));
            const auto& e=f.horizontal?mesh.hEdges(f.edgeLabel):mesh.vEdges(f.edgeLabel);
            const auto a=mesh.verts(e.edgeToVert[0]).coord(),b=mesh.verts(e.edgeToVert[1]).coord();centers.push_back((a+b)*.5);
            area.push_back((f.horizontal?std::abs(b.x-a.x):std::abs(b.y-a.y))*(l.domain().axisymmetric()?double(2*pi)*(f.horizontal?(a.x+b.x)/2:a.x):1));
            minus.push_back(f.negativeCell<0?-1:l.pressureDof(f.negativeCell));plus.push_back(f.positiveCell<0?-1:l.pressureDof(f.positiveCell));
            if(minus.back()>=0&&plus.back()>=0)component[root(plus.back())]=root(minus.back());}
        std::map<int,int> ids;for(int p=0;p<l.pressureCount();++p){const auto r=root(p);if(!ids.count(r))ids[r]=components++;}
        std::vector<int> final(component.size());for(std::size_t p=0;p<component.size();++p)final[p]=ids[root(static_cast<int>(p))];component=std::move(final);
    }
    bool visible(Vector2D a,Vector2D b,int selectedComponent)const
    {
        if(rectangle&&components==1)return true;
        // Independent ordinary-scale ray/slab interval union. Exact sub-ULP
        // synthetic holes are tested directly against the production kernel.
        std::vector<std::array<long double,2>> intervals;
        for(std::size_t p=0;p<boxes.size();++p)if(component[p]==selectedComponent){long double lo=0,hi=1;bool hit=true;
            for(int axis=0;axis<2;++axis){const long double start=axis?a.y:a.x,end=axis?b.y:b.x,
                lower=axis?boxes[p][0].y:boxes[p][0].x,upper=axis?boxes[p][1].y:boxes[p][1].x,d=end-start;
                if(d==0){if(start<lower||start>upper){hit=false;break;}}else{auto t0=(lower-start)/d,t1=(upper-start)/d;if(t0>t1)std::swap(t0,t1);
                    lo=std::max(lo,t0);hi=std::min(hi,t1);if(lo>hi){hit=false;break;}}}
            if(hit)intervals.push_back({lo,hi});}
        std::sort(intervals.begin(),intervals.end());long double covered=0;
        for(const auto& i:intervals){if(i[0]>covered)return false;covered=std::max(covered,i[1]);if(covered>=1)return true;}return false;
    }
    double gradient(int u,const std::vector<double>& c)const
    {return area[u]*((plus[u]<0?0:c[plus[u]])-(minus[u]<0?0:c[minus[u]]));}
};
struct Node{Vector2D point;double weight=0,curvature=0;int pressure=-1;};
std::vector<Node> sourceNodes(const HydroPhaseGeometry& phase,const Reference& ref,const std::string& name,bool emit)
{
    std::vector<Node> selected;const auto& pub=phase.publication();std::size_t number=0;
    for(std::size_t c=0;c<pub.topology().curves.size();++c)for(std::size_t g=0;g<pub.topology().curves[c].graphs().size();++g)
        for(const auto& q:pub.quadrature(c,g)){const auto& f=q.sample.frame;const int p=phase.layout().pressureDof(f.nativeCellIndex);
            const double weight=q.weight*(phase.layout().domain().axisymmetric()?double(2*pi)*f.point.x:1);
            close(q.physicalWeight,weight,1e-14*weight,"independent RZ surface weight counted once");
            check(q.sample.geometryIdentity==phase.geometryIdentity(),"surface q identity mismatch");
            if(emit)std::cout<<"BF03_Q,"<<name<<','<<number<<','<<c<<','<<g<<','<<f.nativeCellIndex<<','<<p<<','
                <<(p<0?-1:ref.component[p])<<','<<f.point.x<<','<<f.point.y<<','<<f.normal.x<<','<<f.normal.y<<','<<q.weight<<','
                <<q.physicalWeight<<','<<q.sample.curvature<<','<<q.sample.azimuthalCurvature<<','<<f.levelSetResidual<<'\n';
            if(p>=0)selected.push_back({f.point,weight,q.sample.curvature,p});++number;}
    return selected;
}
struct RebuiltRule
{
    long double weight=0,moment=0,numerator=0;
    std::size_t count=0,rejectedVisibility=0;
    double minimum=0,maximum=0,maxDistance=0;
};
RebuiltRule rebuildSubarc(const HermiteGalsPublication& pub,const Reference& ref,const HydroCsfFaceSupportReport& report,
    const HydroCsfQuadraturePanel& panel,double begin,double end,const std::string& name,int velocity,std::size_t panelIndex,int half,bool emit)
{
    RebuiltRule result;std::size_t ordinal=0;
    for(const auto& q:pub.quadrature(panel.component,panel.graph,begin,end)){
        const auto& f=q.sample.frame;const auto physical=q.weight*(pub.domain().axisymmetric()?double(2*pi)*f.point.x:1);
        close(q.physicalWeight,physical,1e-14*physical,"subarc physical weight counted once");
        check(q.sample.geometryIdentity==pub.state().identity(),"adaptive source identity changed");
        const auto clock=pub.state().timeLayer(),actual=q.sample.timeLayer;
        check(actual.step==clock.step&&actual.begin==clock.begin&&actual.end==clock.end&&actual.position==clock.position,"adaptive full clock changed");
        const auto distance=std::hypot(f.point.x-report.center.x,f.point.y-report.center.y);
        double kernel=0,weight=0;bool visible=false;
        if(distance<report.radius){visible=ref.visible(report.center,f.point,report.pressureComponent);
            if(visible){const auto ratio=distance/report.radius,base=(1-ratio)*(1+ratio);kernel=base*base*base;weight=q.physicalWeight*kernel;
                check(std::isfinite(weight)&&weight>0,"adaptive positive kernel/weight became unrepresentable");
                result.weight+=weight;result.moment+=static_cast<long double>(weight)*(static_cast<long double>(q.sample.curvature)-report.quadratureAnchor);
                result.numerator+=static_cast<long double>(weight)*q.sample.curvature;
                if(result.count==0)result.minimum=result.maximum=q.sample.curvature;
                else{result.minimum=std::min(result.minimum,q.sample.curvature);result.maximum=std::max(result.maximum,q.sample.curvature);}
                result.maxDistance=std::max(result.maxDistance,distance);++result.count;
            }else ++result.rejectedVisibility;}
        if(emit)std::cout<<"BF03_ADAPT_Q,"<<name<<','<<velocity<<','<<panelIndex<<','<<half<<','<<ordinal<<','<<panel.component<<','<<panel.graph<<','
            <<begin<<','<<end<<','<<q.sample.geometryIdentity.value()<<','<<actual.step<<','<<actual.begin<<','<<actual.end<<','<<static_cast<int>(actual.position)<<','
            <<f.nativeCellIndex<<','<<f.point.x<<','<<f.point.y<<','<<f.normal.x<<','<<f.normal.y<<','<<q.weight<<','<<q.physicalWeight<<','
            <<q.sample.curvature<<','<<distance<<','<<visible<<','<<kernel<<','<<weight<<'\n';
        ++ordinal;
    }
    check(ordinal==8,"subarc lost original GL8 node rule");return result;
}
void combineAudit(RebuiltRule& total,const RebuiltRule& next)
{
    total.weight+=next.weight;total.moment+=next.moment;total.numerator+=next.numerator;total.rejectedVisibility+=next.rejectedVisibility;
    if(next.count){if(total.count==0){total.minimum=next.minimum;total.maximum=next.maximum;}
        else{total.minimum=std::min(total.minimum,next.minimum);total.maximum=std::max(total.maximum,next.maximum);}}
    total.count+=next.count;total.maxDistance=std::max(total.maxDistance,next.maxDistance);
}
RebuiltRule rebuildAccepted(const HydroBalancedCsfLoad& load,const Reference& ref,int velocity,const std::string& name,bool emit)
{
    const auto& report=load.supportReports()[velocity];const auto& panels=load.quadraturePanels(velocity);
    check(report.acceptedPanels==panels.size(),"accepted panel ledger count mismatch");
    check(report.surfaceQueries<=load.curvatureConfig().maxSurfaceQueriesPerFace&&report.geometryQueries<=load.curvatureConfig().maxGeometryQueriesPerFace&&report.shadowWallQueries<=load.curvatureConfig().maxShadowWallQueriesPerFace,
        "adaptive query ledger exceeded explicit work caps");
    RebuiltRule total;long double weightError=0,momentError=0;
    std::map<std::pair<std::size_t,std::size_t>,std::vector<std::array<double,2>>> partitions;
    const auto& pub=load.color().phaseGeometry().publication();
    for(std::size_t i=0;i<panels.size();++i){const auto& panel=panels[i];
        check(panel.begin>=0&&panel.end<=1&&panel.begin<panel.end&&panel.depth>=0&&panel.depth<=load.curvatureConfig().maxQuadratureDepth,
            "invalid retained original subarc panel");
        check(panel.depth<=report.deepestPanel,"panel depth exceeds visited-depth ledger");
        check(panel.coarseSubdivisions>0&&panel.coarseSubdivisions<=load.curvatureConfig().maxSurfaceQueriesPerFace/24,"invalid composite packet count");
        RebuiltRule coarse,fine;
        for(std::size_t j=0;j<panel.coarseSubdivisions;++j){const double a=panel.begin+(panel.end-panel.begin)*(double(j)/panel.coarseSubdivisions),
            b=panel.begin+(panel.end-panel.begin)*(double(j+1)/panel.coarseSubdivisions);
            combineAudit(coarse,rebuildSubarc(pub,ref,report,panel,a,b,name,velocity,i,-1,false));}
        for(std::size_t j=0;j<2*panel.coarseSubdivisions;++j){const double a=panel.begin+(panel.end-panel.begin)*(double(j)/(2*panel.coarseSubdivisions)),
            b=panel.begin+(panel.end-panel.begin)*(double(j+1)/(2*panel.coarseSubdivisions));
            combineAudit(fine,rebuildSubarc(pub,ref,report,panel,a,b,name,velocity,i,static_cast<int>(j),emit));}
        const auto ew=std::abs(fine.weight-coarse.weight),em=std::abs(fine.moment-coarse.moment);
        const auto arithmeticScale=std::max(fine.weight+coarse.weight,1e-300L);
        close(panel.weightErrorEstimate,ew,256*std::numeric_limits<double>::epsilon()*arithmeticScale,"parent coarse/fine weight error ledger");
        close(panel.anchoredMomentErrorEstimate,em,256*std::numeric_limits<double>::epsilon()*arithmeticScale*std::max(1.,report.curvatureScale),
            "parent coarse/fine anchored-moment error ledger");
        weightError+=ew;momentError+=em;combineAudit(total,fine);partitions[{panel.component,panel.graph}].push_back({panel.begin,panel.end});
        if(emit)std::cout<<"BF03_PANEL,"<<name<<','<<velocity<<','<<i<<','<<panel.component<<','<<panel.graph<<','<<panel.begin<<','<<panel.end<<','<<panel.depth<<','
            <<coarse.weight<<','<<fine.weight<<','<<coarse.moment<<','<<fine.moment<<','<<panel.weightErrorEstimate<<','<<panel.anchoredMomentErrorEstimate<<','<<panel.coarseSubdivisions<<','<<panel.wholeGraphEstimate<<'\n';
    }
    for(auto& entry:partitions){auto& intervals=entry.second;std::sort(intervals.begin(),intervals.end());double end=0;
        for(const auto& interval:intervals){check(interval[0]==end,"accepted candidate graph has gap or overlap");end=interval[1];}
        check(end==1,"accepted candidate graph does not cover original [0,1]");}
    check(report.surfaceQueries>=24*panels.size(),"accepted panels exceed coarse/two-fine surface query accounting");
    const auto ledgerScale=std::max(total.weight,1e-300L);
    close(report.weightErrorEstimate,weightError,512*std::numeric_limits<double>::epsilon()*ledgerScale,"aggregate weight error ledger");
    close(report.anchoredMomentErrorEstimate,momentError,512*std::numeric_limits<double>::epsilon()*ledgerScale*std::max(1.,report.curvatureScale),"aggregate moment error ledger");
    if(total.weight>0){const auto estimate=(momentError+std::abs(total.moment/total.weight)*weightError)/(total.weight-weightError);
        close(report.curvatureErrorEstimate,estimate,1e-12*std::max(1.,report.curvatureScale),"normalized curvature error propagation");
        const auto scale=std::max(report.curvatureScale,1/report.radius),tolerance=load.curvatureConfig().relativeQuadratureTolerance;
        check(report.weightErrorEstimate<=.01L*report.totalWeight*(1+1e-12),"adaptive denominator conditioning failed");
        check(report.curvatureErrorEstimate<=tolerance*scale*(1+1e-12),"adaptive curvature estimate failed original budget");}
    return total;
}
struct Continuous{long double weight=0,mean=0;};
std::vector<Shape::Arc> analyticNodes(const Shape& shape,int n)
{std::vector<Shape::Arc> q;q.reserve(n);for(int i=0;i<n;++i){auto a=shape.arc(shape.end()*(i+.5L)/n);a.weight*=shape.end()/n;q.push_back(a);}return q;}
Continuous continuousAverage(const std::vector<Shape::Arc>& nodes,Vector2D center,double radius)
{
    Continuous out;long double numerator=0;for(const auto& q:nodes){const auto distance=std::hypot(static_cast<long double>(q.point.x)-center.x,
        static_cast<long double>(q.point.y)-center.y);if(distance>=radius)continue;
        const auto ratio=distance/radius,base=(1-ratio)*(1+ratio),kernel=base*base*base;
        check(kernel>0&&std::isfinite(kernel),"unrepresentable positive continuous compact-C2 kernel");
        const auto w=q.weight*kernel;out.weight+=w;numerator+=w*q.curvature;}
    if(out.weight>0)out.mean=numerator/out.weight;return out;
}
void missingSupportDiagnostic(const Reference& ref,const HydroRegularizedColor& color,const std::vector<Node>& nodes,const std::string& name)
{
    for(std::size_t u=0;u<ref.area.size();++u){const double gradient=ref.gradient(static_cast<int>(u),color.colors());if(gradient==0)continue;
        const int a=ref.minus[u],b=ref.plus[u],p=a>=0?a:b,component=ref.component[p];const auto center=ref.centers[u];
        const double radius=color.config().epsilon+2*std::max(a<0?0:ref.diagonal[a],b<0?0:ref.diagonal[b]);
        double closest=std::numeric_limits<double>::infinity();std::size_t visible=0,local=0;
        for(const auto& q:nodes)if(ref.component[q.pressure]==component&&ref.visible(center,q.point,component)){
            const auto distance=std::hypot(q.point.x-center.x,q.point.y-center.y);closest=std::min(closest,distance);++visible;if(distance<radius)++local;}
        if(local)continue;
        std::cout<<"BF03_MISSING_FACE,"<<name<<','<<color.config().epsilon<<','<<u<<','<<center.x<<','<<center.y<<','<<radius<<','<<gradient<<','
            <<(a<0||b<0)<<','<<a<<','<<b<<','<<(a<0?0:color.colors()[a])<<','<<(b<0?0:color.colors()[b])<<','
            <<(a<0?0:color.cellReports()[a].transitionVolumeUpper)<<','<<(b<0?0:color.cellReports()[b].transitionVolumeUpper)<<','
            <<closest<<','<<visible<<','<<local<<'\n';
    }
}
HydroBalancedCsfLoad diagnosedProjection(const QuadtreeMesh& mesh,const HydroPressureLoad& pressure,const HydroRegularizedColor& color,
    const Reference& ref,const std::vector<Node>& nodes,const std::string& name,double sigma=.7,HydroCsfCurvatureConfig cfg={})
{
    try{return HydroBalancedCsfLoad::capillary(mesh,pressure,color,sigma,cfg);}
    catch(const std::exception&){missingSupportDiagnostic(ref,color,nodes,name);throw;}
}
void inspect(const QuadtreeMesh& mesh,const HydroPhaseGeometry& phase,const HydroRegularizedColor& color,const std::string& name,
    const Shape* analytic=nullptr,int n=0,int amr=0,bool output=true)
{
    ++scenes;const auto pressure=HydroPressureLoad::build(phase.metrics().divergence());const Reference ref(mesh,phase);
    const auto nodes=sourceNodes(phase,ref,name,output);const auto sharpBefore=phase.liquidFractions();
    const auto load=diagnosedProjection(mesh,pressure,color,ref,nodes,name);load.validateFor(mesh,pressure,color);
    check(load.isCapillary(),"projected load channel lost");rejection([&]{load.constantCurvature();},"no prescribed");
    std::vector<Shape::Arc> coarse,fine;if(analytic){coarse=analyticNodes(*analytic,4096);fine=analyticNodes(*analytic,8192);}
    long double continuousError=0,continuousScale=0,sharpError=0,sharpScale=0,forceError=0,forceScale=0,sharpForceError=0,sharpForceScale=0;
    long double referenceCurvChange=0,referenceForceChange=0;
    std::size_t defined=0;double maximumContinuous=0,referenceChange=0;
    for(int u=0;u<phase.layout().velocityCount();++u){const auto& r=load.supportReports()[u];const auto center=ref.centers[u];
        const int a=ref.minus[u],b=ref.plus[u],p=a>=0?a:b,component=ref.component[p];
        const double radius=color.config().epsilon+2*std::max(a<0?0:ref.diagonal[a],b<0?0:ref.diagonal[b]);
        close(r.center.x,center.x,0,"actual native face center x");close(r.center.y,center.y,0,"actual native face center y");
        close(r.radius,radius,1e-14*radius,"fixed physical compact kernel radius");check(r.pressureComponent==component,"pressure component mismatch");
        const auto rebuilt=rebuildAccepted(load,ref,u,name,output);
        const auto sum=rebuilt.weight,numerator=rebuilt.numerator;const auto min=rebuilt.minimum,max=rebuilt.maximum,maxDistance=rebuilt.maxDistance;
        const auto count=rebuilt.count,rejectedVisibility=rebuilt.rejectedVisibility;
        check(r.supportCount==count,"independent selected/visible support count");
        check(r.visibilityRejected>=rejectedVisibility,"final-node visibility rejections exceed complete trial ledger");
        close(r.totalWeight,sum,1e-12*std::max(sum,1e-300L),"independent positive physical support measure");
        const double gradient=ref.gradient(u,color.colors());close(load.colorGradient()[u],gradient,1e-12*std::max(ref.area[u],1e-300),"independent same-C color incidence");
        if(gradient==0)check(r.status==HydroCsfSupportStatus::ZeroGradient&&!r.curvatureDefined&&count==0&&load.quadraturePanels(u).empty()&&r.surfaceQueries==0,
            "exact zero g pretended adaptive curvature computation");
        double expected=0,mean=0;long double continuum=0,sharp=0;
        if(count){++defined;mean=static_cast<double>(numerator/sum);check(r.curvatureDefined&&r.status==HydroCsfSupportStatus::Projected,"supported curvature undefined");
            peakMean=std::max(peakMean,rel(load.curvatureAverages()[u]-mean,std::max(std::abs(mean),1.)));
            check(rel(load.curvatureAverages()[u]-mean,std::max(std::abs(mean),1.))<=1e-12,"independent accepted-subarc positive curvature mean");
            close(r.curvatureMin,min,0,"support curvature minimum");close(r.curvatureMax,max,0,"support curvature maximum");
            close(r.maxSupportDistance,maxDistance,1e-14*radius,"support maximum physical distance");expected=.7*mean*gradient;
            if(analytic){const auto low=continuousAverage(coarse,center,radius),high=continuousAverage(fine,center,radius);
                check(low.weight>0&&high.weight>0,"analytic continuous compact kernel has no support");continuum=(4*high.mean-low.mean)/3;
                sharp=analytic->sharpCurvature(center);const long double weight=std::abs(gradient);
                continuousError+=weight*std::pow(load.curvatureAverages()[u]-continuum,2);continuousScale+=weight*continuum*continuum;
                sharpError+=weight*std::pow(load.curvatureAverages()[u]-sharp,2);sharpScale+=weight*sharp*sharp;
                const long double error=.7*gradient*(load.curvatureAverages()[u]-continuum),target=.7*gradient*continuum;
                forceError+=error*error;forceScale+=target*target;sharpForceError+=std::pow(.7*gradient*(load.curvatureAverages()[u]-sharp),2);
                sharpForceScale+=std::pow(.7*gradient*sharp,2);
                referenceCurvChange+=weight*std::pow(high.mean-low.mean,2);
                referenceForceChange+=std::pow(.7*gradient*(high.mean-low.mean),2);
                maximumContinuous=std::max(maximumContinuous,std::abs(static_cast<double>(load.curvatureAverages()[u]-continuum)));
                referenceChange=std::max(referenceChange,std::abs(static_cast<double>(high.mean-low.mean)));
                if(analytic->kind=="circle")close(load.curvatureAverages()[u],(analytic->rz?2:1)/analytic->a,1e-10/std::min(analytic->a,analytic->b),"literal geometric circle/sphere curvature");}
        }else{check(!r.curvatureDefined,"no-support curvature invented");
            if(gradient==0)check(r.status==HydroCsfSupportStatus::ZeroGradient,"zero-gradient support status");
            else check(r.status==HydroCsfSupportStatus::PureBoundaryReaction&&(a<0||b<0)&&color.cellReports()[p].transitionVolumeUpper==0,"unsupported force hidden as pure reaction");}
        peakLoad=std::max(peakLoad,rel(load.momentumRhsLoad()[u]-expected,std::max(std::abs(expected),ref.area[u])));
        check(rel(load.momentumRhsLoad()[u]-expected,std::max(std::abs(expected),ref.area[u]))<=1e-12,"positive actual S=sigma*kappa*Ccolor");
        if(ref.area[u]==0){check(load.momentumRhsLoad()[u]==0,"RZ axis force");++axisWitness;}
        if(r.crossesAmr)++amrWitness;if(r.visibilityRejected)++visibilityWitness;
        if(output)std::cout<<"BF03_FACE,"<<name<<','<<n<<','<<amr<<','<<color.config().epsilon<<','<<u<<','<<a<<','<<b<<','<<ref.area[u]<<','
            <<center.x<<','<<center.y<<','<<radius<<','<<r.pressureComponent<<','<<static_cast<int>(r.status)<<','<<r.curvatureDefined<<','
            <<(r.curvatureDefined?load.curvatureAverages()[u]:std::numeric_limits<double>::quiet_NaN())<<','<<r.supportCount<<','<<r.totalWeight<<','
            <<r.curvatureMin<<','<<r.curvatureMax<<','<<r.maxSupportDistance<<','<<r.candidateQueries<<','<<r.visibilityRejected<<','<<r.visibilityCellQueries<<','
            <<r.boundary<<','<<r.nearBoundary<<','<<r.nearAxis<<','<<r.crossesAmr<<','<<r.subResolutionGradient<<','<<r.gradientRoundoffScale<<','
            <<r.colorGradientErrorEstimate<<','<<gradient<<','<<load.momentumRhsLoad()[u]<<','<<mean<<','<<continuum<<','<<sharp<<','
            <<r.surfaceQueries<<','<<r.geometryQueries<<','<<r.acceptedPanels<<','<<r.deepestPanel<<','<<r.weightErrorEstimate<<','<<r.anchoredMomentErrorEstimate<<','
            <<r.curvatureErrorEstimate<<','<<r.curvatureScale<<','<<r.quadratureAnchor<<','<<r.shadowWallQueries<<','<<r.quadratureAnchorDefined<<','<<r.trialPositiveGraphsWithoutFineSupport<<'\n';}
    check(phase.liquidFractions()==sharpBefore,"curvature projection changed sharp material");
    std::cout<<"BF03_ERROR,"<<name<<','<<n<<','<<amr<<','<<color.config().epsilon<<','<<defined<<','
        <<(continuousScale>0?std::sqrt(continuousError/continuousScale):0)<<','<<(forceScale>0?std::sqrt(forceError/forceScale):0)<<','
        <<(sharpScale>0?std::sqrt(sharpError/sharpScale):0)<<','<<(sharpForceScale>0?std::sqrt(sharpForceError/sharpForceScale):0)<<','
        <<maximumContinuous<<','<<referenceChange<<','<<load.curvatureReport().candidateQueries<<','<<load.curvatureReport().visibilityRejected<<','
        <<(continuousScale>0?std::sqrt(referenceCurvChange/continuousScale):0)<<','
        <<(forceScale>0?std::sqrt(referenceForceChange/forceScale):0)<<'\n';
}
void kernels()
{
    for(double value:{0.,1.,-3.,1e120,1e-280})for(const auto& weights:{std::vector<double>{1,2,0,4},std::vector<double>{1e308,1e-308,0,1e308}})
        check(HydroBalancedCsfLoad::positiveCurvatureAverage(weights,std::vector<double>(weights.size(),value))==value,"constant positive-average reproduction not exact");
    close(HydroBalancedCsfLoad::positiveCurvatureAverage({1,2,3},{-2,0,4}),10.L/6,1e-15,"independent positive mean");
    close(HydroBalancedCsfLoad::positiveCurvatureAverage({1e308,1e308},{1e120,-1e120}),0,0,"normalized weight overflow avoided");
    rejection([]{HydroBalancedCsfLoad::positiveCurvatureAverage({},{});},"dimensions");
    rejection([]{HydroBalancedCsfLoad::positiveCurvatureAverage({1},{1,2});},"dimensions");
    rejection([]{HydroBalancedCsfLoad::positiveCurvatureAverage({0,0},{1,2});},"positive");
    rejection([]{HydroBalancedCsfLoad::positiveCurvatureAverage({-1,1},{1,2});},"input");
    rejection([]{HydroBalancedCsfLoad::positiveCurvatureAverage({1},{std::numeric_limits<double>::quiet_NaN()});},"input");
    rejection([]{HydroBalancedCsfLoad::positiveCurvatureAverage({1e308,1e-308,1e-308},{0,1e308,-5e307});},"unrepresentable");
    rejection([]{HydroBalancedCsfLoad::positiveCurvatureAverage({1,1e-320,1e-320},{1e-308,1e308,-1e308});},"unrepresentable");
    const double join=.875,gap=std::nextafter(join,std::numeric_limits<double>::infinity());
    const std::vector<Box> touching{{{{0,-1},{join,1}}},{{{join,-1},{3,1}}}},withHole{{{{0,-1},{join,1}}},{{{gap,-1},{3,1}}}};
    for(bool reverse:{false,true}){const Vector2D a=reverse?Vector2D{3,0}:Vector2D{0,0},b=reverse?Vector2D{0,0}:Vector2D{3,0};
        check(HydroBalancedCsfLoad::selectedSegmentVisible(a,b,touching),"closed touching cells not visible");
        check(!HydroBalancedCsfLoad::selectedSegmentVisible(a,b,withHole),"positive-width hole bridged by rounded ray fractions");}
    check(HydroBalancedCsfLoad::selectedSegmentVisible({join,0},{join,0},touching),"boundary zero-length segment lost");
    check(!HydroBalancedCsfLoad::selectedSegmentVisible({-1,0},{0,0},touching),"uncovered segment accepted");
    std::cout<<"BF03_KERNEL,PASS\n";
}
template<class T>void appendBits(std::vector<unsigned char>& output,const T& value)
{const auto* first=reinterpret_cast<const unsigned char*>(&value);output.insert(output.end(),first,first+sizeof(value));}
void appendFrame(std::vector<unsigned char>& bits,const LevelSetHermiteGeometry::Frame& f)
{
    for(double value:{f.point.x,f.point.y,f.normal.x,f.normal.y,f.curvature,f.levelSetResidual,f.stationarityResidual,
        f.cellLocalPoint.x,f.cellLocalPoint.y,f.cellChartPoint.x.offset,f.cellChartPoint.y.offset})appendBits(bits,value);
    for(int value:{f.iterations,f.attemptedSeeds,f.convergedSeeds,f.cellI,f.cellJ,f.nativeCellIndex})appendBits(bits,value);
    for(bool value:{f.hasCellLocalPoint,f.hasCellChartPoint,f.cellChartPoint.x.fromUpper,f.cellChartPoint.y.fromUpper})appendBits(bits,value);
}
std::vector<unsigned char> quadratureBits(const std::vector<HermiteGalsSurfaceQuadrature>& nodes)
{
    std::vector<unsigned char> bits;appendBits(bits,nodes.size());
    for(const auto& q:nodes){appendBits(bits,q.sample.geometryIdentity.value());appendBits(bits,q.sample.timeLayer.step);
        appendBits(bits,q.sample.timeLayer.begin);appendBits(bits,q.sample.timeLayer.end);appendBits(bits,q.sample.timeLayer.position);
        appendFrame(bits,q.sample.frame);appendBits(bits,q.sample.azimuthalCurvature);appendBits(bits,q.sample.curvature);
        appendBits(bits,q.weight);appendBits(bits,q.physicalWeight);}return bits;
}
std::vector<unsigned char> authorityBits(const HermiteGalsPublication& pub)
{
    std::vector<unsigned char> bits;const auto& state=pub.state();const auto& p=state.primary();const auto& c=state.closure();
    appendBits(bits,state.identity().value());const auto clock=state.timeLayer();appendBits(bits,clock.step);appendBits(bits,clock.begin);appendBits(bits,clock.end);appendBits(bits,clock.position);
    appendBits(bits,p.measure);appendBits(bits,p.transport);appendBits(bits,p.vertices.size());appendBits(bits,p.edges.size());
    for(const auto& v:p.vertices){appendBits(bits,v.vertex);appendBits(bits,v.phi);appendBits(bits,v.dx);appendBits(bits,v.dy);}
    for(const auto& e:p.edges){appendBits(bits,e.vertices[0]);appendBits(bits,e.vertices[1]);appendBits(bits,e.phi);appendBits(bits,e.dn);}
    for(int cell:p.geometryCells)appendBits(bits,cell);for(const auto& m:p.phaseVolumes){appendBits(bits,m.cell);appendBits(bits,m.liquid);}
    for(double value:c.vertexMixedDerivatives)appendBits(bits,value);for(double value:c.cellBubbleCoefficients)appendBits(bits,value);
    const auto& r=pub.report();for(std::size_t value:{r.components,r.graphs,r.quadratureNodes,r.halfFaces,r.axisEndpoints})appendBits(bits,value);
    appendBits(bits,r.maxRelativeInventoryResidual);appendBits(bits,r.maxRelativeInterfaceResidual);appendBits(bits,r.maxUnitNormalDefect);appendBits(bits,r.maxAbsCurvature);
    const auto& v=pub.volumes().report();for(std::size_t value:{v.liquidCells,v.gasCells,v.cutCells,v.components,v.graphs,v.faceRootIncidences})appendBits(bits,value);
    for(const long double value:{v.totalCapacity,v.totalRawLiquid,v.totalInitialLiquid,v.relativeBoundaryClosure,v.relativeMomentEstimate,
        v.maxRelativeCapacityDefect,v.maxRelativePureDefect,v.totalSourceInventory,v.geometricInventoryDifference,v.maxCellRelativeInventoryDifference})appendBits(bits,value);
    return bits;
}
void subarcContract(const HydroPhaseGeometry& phase)
{
    const auto& pub=phase.publication();const auto before=authorityBits(pub);const auto* source=&pub.state().geometry().snapshot();
    for(std::size_t c=0;c<pub.topology().curves.size();++c)for(std::size_t g=0;g<pub.topology().curves[c].graphs().size();++g){
        const auto original=quadratureBits(pub.quadrature(c,g));const auto* cached=&pub.quadrature(c,g);
        check(original==quadratureBits(pub.quadrature(c,g,0,1)),"full-span subarc changed any cached field bits");
        const auto& chain=pub.topology().curves[c];const auto native=chain.quadrature(static_cast<int>(g),.125,.875);
        const auto actual=pub.quadrature(c,g,.125,.875);check(actual.size()==native.size(),"publication changed original interval GL rule");
        for(std::size_t i=0;i<actual.size();++i){std::vector<unsigned char> left,right;appendFrame(left,actual[i].sample.frame);appendFrame(right,native[i].frame);
            check(left==right,"subarc frame differs from original native root chart");
            close(actual[i].weight,native[i].weight,0,"subarc ds changed native Jacobian");close(actual[i].physicalWeight,native[i].physicalWeight,0,"subarc physical weight changed native Jacobian");
            check(actual[i].sample.geometryIdentity==phase.geometryIdentity(),"subarc detached source authority");
            const auto clock=actual[i].sample.timeLayer,expected=phase.timeLayer();
            check(clock.step==expected.step&&clock.begin==expected.begin&&clock.end==expected.end&&clock.position==expected.position,"subarc detached full clock");}
        check(original==quadratureBits(pub.quadrature(c,g))&&cached==&pub.quadrature(c,g),"subarc changed original GL8 cache owner/data");
    }
    check(before==authorityBits(pub)&&source==&pub.state().geometry().snapshot(),"subarc queries changed accepted state/geometry/report");
    rejection([&]{pub.quadrature(0,0,-.1,.5);},"interval");rejection([&]{pub.quadrature(0,0,.5,1.1);},"interval");
    rejection([&]{pub.quadrature(0,0,.75,.25);},"interval");
    rejection([&]{pub.quadrature(0,0,std::numeric_limits<double>::quiet_NaN(),1);},"interval");
    rejection([&]{pub.quadrature(0,pub.topology().curves[0].graphs().size(),.25,.75);},"graph index");
    check(before==authorityBits(pub),"subarc rejection mutated source/report");std::cout<<"BF03_SUBARC,PASS\n";
}
void runCase(const std::string& kind,bool rz,int amr,int n)
{
    check(n>=4&&n<=64&&amr>=0&&amr<=2,"invalid case mesh request");auto mesh=makeMesh(n,amr);const auto shape=shapeFor(kind,rz);
    const auto phase=HydroPhaseGeometry::build(*mesh,metrics(*mesh,rz),publish(*mesh,shape.profile(),rz,shape.planes()));
    const auto name=kind+(rz?"-rz":"-xy");
    for(double ratio:{1.5,2.}){const auto color=colorFor(*mesh,phase,ratio/n);inspect(*mesh,phase,color,name,&shape,n,amr);}
}
void basic()
{
    // Revisit the original failed inputs under the uniformly revised local
    // two-diagonal kernel. The version marker keeps the original V1 evidence.
    for(bool rz:{false,true}){auto mesh=makeMesh(4,1);const auto shape=shapeFor("circle",rz);
        const auto phase=HydroPhaseGeometry::build(*mesh,metrics(*mesh,rz),publish(*mesh,shape.profile(),rz));
        inspect(*mesh,phase,colorFor(*mesh,phase,.375),rz?"basic-sphere":"basic-circle",&shape,4,1);}
    {auto mesh=makeMesh(8);const auto shape=shapeFor("circle",false);
        const auto phase=HydroPhaseGeometry::build(*mesh,metrics(*mesh,false),publish(*mesh,shape.profile(),false));
        inspect(*mesh,phase,colorFor(*mesh,phase,.1875),"original-n8-circle",&shape,8,0);}
}
void domains()
{
    for(bool rz:{false,true})for(double sign:{-1.,1.}){auto mesh=makeMesh(4);const Profile f=[sign](Vector2D){return Jet{sign,0,0,0};};
        const auto phase=HydroPhaseGeometry::build(*mesh,metrics(*mesh,rz),publish(*mesh,f,rz));const auto color=colorFor(*mesh,phase,.375);
        inspect(*mesh,phase,color,sign<0?"pure-liquid":"pure-gas");const auto pressure=HydroPressureLoad::build(phase.metrics().divergence());
        const auto projected=HydroBalancedCsfLoad::capillary(*mesh,pressure,color,.7),constant=HydroBalancedCsfLoad::constantCurvature(*mesh,pressure,color,.7,2);
        bool pureWall=false;for(std::size_t u=0;u<projected.colorGradient().size();++u){check(projected.momentumRhsLoad()[u]==0,"pure variable CSF not physically zero");
            if(projected.supportReports()[u].status==HydroCsfSupportStatus::PureBoundaryReaction){pureWall=true;check(constant.momentumRhsLoad()[u]!=0,"constant full wall stress lost");}}
        if(sign<0)check(pureWall,"pure-liquid boundary reaction witness missing");}
    auto mesh=makeMesh(8);Shape shape=shapeFor("circle",false);shape.cx=.23;shape.cy=.5;shape.a=shape.b=.12;
    const auto phase=HydroPhaseGeometry::build(*mesh,metrics(*mesh,false,1),publish(*mesh,shape.profile(),false));
    const auto shadowed=colorFor(*mesh,phase,.125);
    rejection([&]{inspect(*mesh,phase,shadowed,"selected-hole-shadow-unsupported");},"without local selected/visible curvature support");
    shape.cx=shape.cy=.5;shape.a=shape.b=.30;
    const auto visiblePhase=HydroPhaseGeometry::build(*mesh,metrics(*mesh,false,1),publish(*mesh,shape.profile(),false));
    inspect(*mesh,visiblePhase,colorFor(*mesh,visiblePhase,.125),"selected-hole-surrounding");
    check(visibilityWitness>0,"actual selected hole visibility witness missing");
}
void provenance()
{
    auto mesh=makeMesh(8);const auto shape=shapeFor("circle",false);
    const auto phase=HydroPhaseGeometry::build(*mesh,metrics(*mesh,false),publish(*mesh,shape.profile(),false));
    subarcContract(phase);
    {const auto radial=shapeFor("circle",true);const auto radialPhase=HydroPhaseGeometry::build(*mesh,metrics(*mesh,true),publish(*mesh,radial.profile(),true));
        subarcContract(radialPhase);}
    const auto color=colorFor(*mesh,phase,.1875);const auto pressure=HydroPressureLoad::build(phase.metrics().divergence());
    const auto load=HydroBalancedCsfLoad::capillary(*mesh,pressure,color,.7),copy=load;load.validateFor(*mesh,pressure,color);copy.validateFor(*mesh,pressure,color);
    check(&load.supportReports()==&copy.supportReports()&&&load.color().colors()==&color.colors(),"immutable projection owner copied data");
    const auto zero=HydroBalancedCsfLoad::capillary(*mesh,pressure,color,0);check(zero.curvatureReport().projectionSkippedSigmaZero,"sigma zero pretended projection");
    for(std::size_t u=0;u<zero.momentumRhsLoad().size();++u){check(zero.momentumRhsLoad()[u]==0&&!std::signbit(zero.momentumRhsLoad()[u]),"sigma zero not exact positive zero");
        check(zero.supportReports()[u].status==HydroCsfSupportStatus::SigmaZero&&!zero.supportReports()[u].curvatureDefined&&zero.quadraturePanels(u).empty(),"sigma zero invented curvature");}
    rejection([&]{zero.constantCurvature();},"no prescribed");
    const auto constant=HydroBalancedCsfLoad::constantCurvature(*mesh,pressure,color,.7,2);
    rejection([&]{constant.curvatureAverages();},"no projected");rejection([&]{constant.supportReports();},"no curvature support");
    rejection([&]{load.validateFor(*mesh,pressure,HermiteGalsIdentity{},color.timeLayer());},"identity");
    for(auto clock:{HermiteGalsTimeLayer{1,0,0,Position::Initial},HermiteGalsTimeLayer{0,0,0,Position::End},HermiteGalsTimeLayer{0,0,1,Position::Initial}})
        rejection([&]{load.validateFor(*mesh,pressure,color.geometryIdentity(),clock);},"time");
    auto bad=color.timeLayer();bad.end=std::numeric_limits<double>::quiet_NaN();rejection([&]{load.validateFor(*mesh,pressure,color.geometryIdentity(),bad);},"time");
    const auto otherPressure=HydroPressureLoad::build(phase.metrics().divergence());rejection([&]{load.validateFor(*mesh,otherPressure,color);},"pressure owner");
    const auto wrongD=HydroPressureLoad::build(HydroIntegralDivergence::build(*mesh,phase.layout()));
    rejection([&]{HydroBalancedCsfLoad::capillary(*mesh,wrongD,color,.7);},"divergence owner");
    const auto otherColor=colorFor(*mesh,phase,.25);rejection([&]{load.validateFor(*mesh,pressure,otherColor);},"color owner");
    for(double sigma:{-1.,std::numeric_limits<double>::infinity(),std::numeric_limits<double>::quiet_NaN()})
        rejection([&]{HydroBalancedCsfLoad::capillary(*mesh,pressure,color,sigma);},"surface tension");
    auto cfg=HydroCsfCurvatureConfig{};cfg.maxCandidateQueriesPerFace=0;
    rejection([&]{HydroBalancedCsfLoad::capillary(*mesh,pressure,color,0,cfg);},"config");cfg.maxCandidateQueriesPerFace=1;
    rejection([&]{HydroBalancedCsfLoad::capillary(*mesh,pressure,color,.7,cfg);},"support query limit");
    cfg=HydroCsfCurvatureConfig{};cfg.maxVisibilityCellQueriesPerFace=0;
    rejection([&]{HydroBalancedCsfLoad::capillary(*mesh,pressure,color,0,cfg);},"config");
    // The full rectangle deliberately needs no visibility leaf queries. A
    // real nonconvex selected domain must exercise the separate leaf budget.
    {auto holeMesh=makeMesh(8);auto holeShape=shapeFor("circle",false);holeShape.cx=holeShape.cy=.5;holeShape.a=holeShape.b=.30;
        const auto holePhase=HydroPhaseGeometry::build(*holeMesh,metrics(*holeMesh,false,1),publish(*holeMesh,holeShape.profile(),false));
        const auto holeColor=colorFor(*holeMesh,holePhase,.125);const auto holePressure=HydroPressureLoad::build(holePhase.metrics().divergence());
        cfg=HydroCsfCurvatureConfig{};cfg.maxVisibilityCellQueriesPerFace=1;
        rejection([&]{HydroBalancedCsfLoad::capillary(*holeMesh,holePressure,holeColor,.7,cfg);},"visibility cell query limit");
        cfg=HydroCsfCurvatureConfig{};cfg.maxShadowWallQueriesPerFace=0;
        rejection([&]{HydroBalancedCsfLoad::capillary(*holeMesh,holePressure,holeColor,0,cfg);},"config");
        cfg.maxShadowWallQueriesPerFace=1;
        rejection([&]{HydroBalancedCsfLoad::capillary(*holeMesh,holePressure,holeColor,.7,cfg);},"shadow wall query limit");}
    cfg=HydroCsfCurvatureConfig{};cfg.maxSurfaceQueriesPerFace=0;
    rejection([&]{HydroBalancedCsfLoad::capillary(*mesh,pressure,color,0,cfg);},"config");cfg.maxSurfaceQueriesPerFace=1;
    rejection([&]{HydroBalancedCsfLoad::capillary(*mesh,pressure,color,.7,cfg);},"surface query limit");
    cfg=HydroCsfCurvatureConfig{};cfg.maxGeometryQueriesPerFace=0;
    rejection([&]{HydroBalancedCsfLoad::capillary(*mesh,pressure,color,0,cfg);},"config");cfg.maxGeometryQueriesPerFace=1;
    rejection([&]{HydroBalancedCsfLoad::capillary(*mesh,pressure,color,.7,cfg);},"geometry query limit");
    cfg=HydroCsfCurvatureConfig{};cfg.maxQuadratureDepth=0;
    rejection([&]{HydroBalancedCsfLoad::capillary(*mesh,pressure,color,0,cfg);},"config");cfg.maxQuadratureDepth=41;
    rejection([&]{HydroBalancedCsfLoad::capillary(*mesh,pressure,color,0,cfg);},"config");
    cfg=HydroCsfCurvatureConfig{};cfg.maxQuadratureDepth=1;cfg.relativeQuadratureTolerance=1e-16;
    rejection([&]{HydroBalancedCsfLoad::capillary(*mesh,pressure,color,.7,cfg);},"quadrature depth");
    for(double tolerance:{0.,1.,std::numeric_limits<double>::infinity(),std::numeric_limits<double>::quiet_NaN()}){
        cfg=HydroCsfCurvatureConfig{};cfg.relativeQuadratureTolerance=tolerance;
        rejection([&]{HydroBalancedCsfLoad::capillary(*mesh,pressure,color,0,cfg);},"config");}
    const Profile plane=[](Vector2D p){return Jet{p.x-.6,1,0,0};};
    const auto missingPhase=HydroPhaseGeometry::build(*mesh,metrics(*mesh,false,2),publish(*mesh,plane,false,{{1,0},{1,1}}));
    const auto missingColor=colorFor(*mesh,missingPhase,.1875);const auto missingPressure=HydroPressureLoad::build(missingPhase.metrics().divergence());
    rejection([&]{HydroBalancedCsfLoad::capillary(*mesh,missingPressure,missingColor,.7);},"without local selected/visible curvature support");
    const auto harmless=HydroBalancedCsfLoad::capillary(*mesh,missingPressure,missingColor,0);
    check(harmless.curvatureReport().projectionSkippedSigmaZero,"sigma zero missing-support bypass dishonest");
    auto target=makeMesh(8,1);rejection([&]{load.validateFor(*target,pressure,color);},"mesh");
    const auto transfer=transferHermiteGalsFirstJetMesh(*mesh,phase.publication(),*target,HermiteGalsMeshTransferKind::ExactRefinement);
    const auto nextPhase=HydroPhaseGeometry::build(*target,metrics(*target,false),transfer.publication());const auto nextColor=colorFor(*target,nextPhase,.1875);
    inspect(*target,nextPhase,nextColor,"actual-regrid");
    const auto nextPressure=HydroPressureLoad::build(nextPhase.metrics().divergence());const auto next=HydroBalancedCsfLoad::capillary(*target,nextPressure,nextColor,.7);
    rejection([&]{next.validateFor(*target,nextPressure,color.geometryIdentity(),nextColor.timeLayer());},"identity");
}
void extraProvenance()
{
    auto mesh=makeMesh(8);const auto shape=shapeFor("circle",false);
    HermiteGalsPublicationConfig config;config.relativeInterfaceResidualTolerance=1e-10;config.unitNormalTolerance=1e-13;
    const auto custom=HydroPhaseGeometry::build(*mesh,metrics(*mesh,false),publish(*mesh,shape.profile(),false,{},config));
    subarcContract(custom);std::cout<<"BF03_CUSTOM_CONFIG,"<<config.relativeInterfaceResidualTolerance<<','<<config.unitNormalTolerance<<",PASS\n";
    const Profile separate=[](Vector2D p){constexpr double radius=.12;const double scale=1/(radius*radius);
        const double x0=p.x-.18,y0=p.y-.48,x1=p.x-.82,y1=p.y-.52,
            a=(x0*x0+y0*y0)*scale-1,b=(x1*x1+y1*y1)*scale-1,
            ax=2*x0*scale,ay=2*y0*scale,bx=2*x1*scale,by=2*y1*scale;
        return Jet{a*b,ax*b+a*bx,ay*b+a*by,ax*by+ay*bx};};
    const auto disconnected=HydroPhaseGeometry::build(*mesh,metrics(*mesh,false,3),publish(*mesh,separate,false));
    const Reference ref(*mesh,disconnected);check(ref.components==2,"actual disconnected pressure components missing");
    check(disconnected.publication().topology().curves.size()==2,"actual disconnected surface components missing");
    inspect(*mesh,disconnected,colorFor(*mesh,disconnected,.1875),"actual-disconnected");
    std::cout<<"BF03_DISCONNECTED,"<<ref.components<<",PASS\n";
    const Profile plane=[](Vector2D p){return Jet{p.y-.5,0,1,0};};
    const auto tangent=HydroPhaseGeometry::build(*mesh,metrics(*mesh,false),publish(*mesh,plane,false,{{0,0},{0,1}}));
    const double epsilon=.4375-2*std::hypot(.125,.125);
    const auto planeColor=colorFor(*mesh,tangent,epsilon);
    const auto planeLoad=HydroBalancedCsfLoad::capillary(*mesh,HydroPressureLoad::build(tangent.metrics().divergence()),planeColor,.7);
    bool witnessed=false;for(const auto& support:planeLoad.supportReports())if(support.center.x==.4375&&support.center.y==.5){
        check(support.radius==.4375&&support.curvatureDefined,"exact radial endpoint-touch witness missing");witnessed=true;}
    check(witnessed,"endpoint-touch native face absent");
    inspect(*mesh,tangent,planeColor,"kernel-endpoint-touch");
    std::cout<<"BF03_KERNEL_ENDPOINT_TOUCH,PASS\n";
}
}
int main(int argc,char** argv)
{
    try{std::cout<<std::setprecision(17);
    std::cout<<"BF03_KERNEL_VERSION,same-source-radial-adaptive-v7\n";
        std::cout<<"BF03_Q_HEADER,shape,index,surface_component,graph,native_cell,pressure,pressure_component,x,y,nx,ny,ds_weight,physical_weight,curvature,azimuthal_curvature,interface_residual\n";
        std::cout<<"BF03_FACE_HEADER,shape,n,amr,epsilon,velocity,minus,plus,area,x,y,radius,component,status,defined,kappa,count,weight,min,max,max_distance,candidates,visibility_rejected,visibility_cells,boundary,near_boundary,near_axis,amr_crossing,subresolution,roundoff_scale,color_error_estimate,Cc,S,q_reference,continuous_reference,sharp_reference,surface_queries,geometry_queries,accepted_panels,deepest_panel,weight_error_estimate,anchored_moment_error_estimate,curvature_error_estimate,curvature_scale,quadrature_anchor,shadow_wall_queries,quadrature_anchor_defined,trial_positive_graphs_without_fine_support\n";
        std::cout<<"BF03_ADAPT_Q_HEADER,shape,velocity,panel,half,node,surface_component,graph,begin,end,geometry_identity,step,clock_begin,clock_end,clock_position,native_cell,x,y,nx,ny,ds_weight,physical_weight,curvature,distance,visible,kernel,positive_weight\n";
        std::cout<<"BF03_PANEL_HEADER,shape,velocity,panel,surface_component,graph,begin,end,depth,coarse_weight,fine_weight,coarse_anchor_moment,fine_anchor_moment,weight_error_estimate,anchored_moment_error_estimate,coarse_subdivisions,whole_graph_estimate\n";
        std::cout<<"BF03_ERROR_HEADER,shape,n,amr,epsilon,defined,continuous_curvature_relative,continuous_force_relative,sharp_curvature_relative,sharp_force_relative,max_continuous_error,reference_4096_8192_change,candidates,visibility_rejected,reference_curvature_change_relative,reference_force_change_relative\n";
        std::cout<<"BF03_MISSING_FACE_HEADER,shape,epsilon,velocity,x,y,R,Cc,boundary,minus,plus,cminus,cplus,minus_transition_upper,plus_transition_upper,min_visible_same_component_q_distance,visible_q_count,local_q_count\n";
        if(argc==6&&std::string(argv[1])=="--case")runCase(argv[2],std::stoi(argv[3])!=0,std::stoi(argv[4]),std::stoi(argv[5]));
        else if(argc==2&&std::string(argv[1])=="--kernels")kernels();
        else if(argc==2&&std::string(argv[1])=="--provenance")provenance();
        else if(argc==2&&std::string(argv[1])=="--extra-provenance")extraProvenance();
        else if(argc==2&&std::string(argv[1])=="--domains")domains();
        else{check(argc==1,"unknown BF03 contract arguments");kernels();basic();domains();provenance();check(amrWitness>0&&axisWitness>0&&visibilityWitness>0,"coverage witness missing");}
        std::cout<<"BF03_CONTRACT,"<<scenes<<','<<guards<<','<<peakMean<<','<<peakLoad<<",PASS\n";return 0;
    }catch(const std::exception& e){std::cerr<<"BF03_CURVATURE_FAIL,"<<e.what()<<'\n';return 1;}
}
