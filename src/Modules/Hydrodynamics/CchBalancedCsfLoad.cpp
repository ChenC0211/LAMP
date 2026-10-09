#include "Modules/Hydrodynamics/CchBalancedCsfLoad.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Geometry/Detail/ExactFaceTimeRational.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>

namespace Quantum::Hydrodynamics
{
namespace
{
void require(bool value,const char* reason)
{if(!value)throw std::invalid_argument(std::string("[CchBalancedCsfLoad] ")+reason);}
struct Box {Vector2D lo,hi;};
namespace Exact=Quantum::Geometry::Detail::ExactExpansion;
using RayParameter=Quantum::Geometry::Detail::FaceTime::FaceTimeRational;
int compare(const RayParameter& a,const RayParameter& b)
{return Exact::sign(Exact::subtract(Exact::multiply(a.n,b.d),Exact::multiply(b.n,a.d)));}
struct RayLine
{
    Vector2D begin;
    std::array<Exact::Expansion,2> delta;
    RayLine(Vector2D a,Vector2D b):begin(a),delta{{Exact::subtract(Exact::fromDouble(b.x),Exact::fromDouble(a.x)),
        Exact::subtract(Exact::fromDouble(b.y),Exact::fromDouble(a.y))}}{}
    bool clip(Box box,RayParameter& low,RayParameter& high) const
    {
        for(int axis=0;axis<2;++axis){const double start=axis?begin.y:begin.x,lo=axis?box.lo.y:box.lo.x,hi=axis?box.hi.y:box.hi.x;
            if(Exact::sign(delta[axis])==0){if(start<lo||start>hi)return false;}
            else{RayParameter a(Exact::subtract(Exact::fromDouble(lo),Exact::fromDouble(start)),delta[axis]),
                    b(Exact::subtract(Exact::fromDouble(hi),Exact::fromDouble(start)),delta[axis]);
                if(compare(a,b)>0)std::swap(a,b);if(compare(a,low)>0)low=std::move(a);if(compare(b,high)<0)high=std::move(b);}
            if(compare(low,high)>0)return false;}
        return true;
    }
};
using RayIntervals=std::vector<std::pair<RayParameter,RayParameter>>;
bool covered(RayIntervals& intervals)
{
    std::sort(intervals.begin(),intervals.end(),[](const auto& a,const auto& b){const int order=compare(a.first,b.first);
        return order<0||(order==0&&compare(a.second,b.second)<0);});
    RayParameter end(0),one(1);
    for(const auto& span:intervals){if(compare(span.first,end)>0)return false;if(compare(span.second,end)>0)end=span.second;
        if(compare(end,one)>=0)return true;}
    return false;
}
Box surrounding(Vector2D p,double radius)
{
    const double inf=std::numeric_limits<double>::infinity();
    return {{std::nextafter(p.x-radius,-inf),std::nextafter(p.y-radius,-inf)},
            {std::nextafter(p.x+radius,inf),std::nextafter(p.y+radius,inf)}};
}
// Bins only accelerate the fixed physical support. They do not truncate it
// to a ring count or alter either the native mesh or the source quadrature.
struct Buckets
{
    using Key=std::pair<long long,long long>;
    Vector2D origin{0,0};double width=1;
    std::map<Key,std::vector<int>> data;
    Key key(Vector2D p) const
    {
        const double x=std::floor((p.x-origin.x)/width),y=std::floor((p.y-origin.y)/width);
        require(std::isfinite(x)&&std::isfinite(y)&&std::abs(x)<1e15&&std::abs(y)<1e15,"unrepresentable CSF spatial index");
        return {static_cast<long long>(x),static_cast<long long>(y)};
    }
    template<class F>void visit(Box b,F&& action) const
    {
        const auto lo=key(b.lo),hi=key(b.hi);
        require(hi.first>=lo.first&&hi.second>=lo.second&&hi.first-lo.first<=8&&hi.second-lo.second<=8,"unresolved CSF local index extent");
        for(auto x=lo.first;x<=hi.first;++x)for(auto y=lo.second;y<=hi.second;++y){
            const auto it=data.find({x,y});if(it!=data.end())for(int id:it->second)action(id);}
    }
    void insert(Box b,int id)
    {
        const auto lo=key(b.lo),hi=key(b.hi);
        require(hi.first-lo.first<=2&&hi.second-lo.second<=2,"unresolved CSF index insertion");
        for(auto x=lo.first;x<=hi.first;++x)for(auto y=lo.second;y<=hi.second;++y)data[{x,y}].push_back(id);
    }
};
struct SurfacePoint {Vector2D point;double weight,curvature;int pressure,level;};
struct ProjectionIndex
{
    std::vector<Box> cells,walls;
    std::vector<double> diagonal;
    std::vector<int> component,level;
    std::vector<NativeLeafFace> faces;
    std::vector<SurfacePoint> points;
    Buckets cellBins,wallBins,pointBins;
    std::vector<int> scratch;
    std::vector<std::size_t> cellVisitStamp;
    std::size_t visitEpoch=0;
    RayIntervals intervals;
    Box selectedBounds;bool convexRectangle=false;
    ProjectionIndex(const QuadtreeMesh& mesh,const CchRegularizedColor& color,CchCsfCurvatureReport& report)
    {
        const auto& phase=color.phaseGeometry();const auto& l=phase.layout();const int np=l.pressureCount();
        require(np>0,"empty CSF pressure domain");cellVisitStamp.resize(np,0);
        const auto& scalar=phase.publication().state().geometry().snapshot();
        Vector2D origin{std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity()};double maxDiagonal=0;
        for(int p=0;p<np;++p){const int cell=l.nativeCell(p);const auto b=scalar.nativeCellBounds(cell);
            cells.push_back({b[0],b[1]});const auto h=b[1]-b[0];const double d=std::hypot(h.x,h.y);
            require(std::isfinite(d)&&d>0,"unrepresentable CSF leaf diameter");diagonal.push_back(d);maxDiagonal=std::max(maxDiagonal,d);
            origin.x=std::min(origin.x,b[0].x);origin.y=std::min(origin.y,b[0].y);level.push_back(mesh.cells(mesh.cellIndexToLabel[cell]).level());}
        const double width=color.config().epsilon+maxDiagonal;
        require(std::isfinite(width)&&width>0,"unrepresentable CSF support radius");
        for(auto* bins:{&cellBins,&wallBins,&pointBins}){bins->origin=origin;bins->width=width;}
        for(int p=0;p<np;++p)cellBins.insert(cells[p],p);
        selectedBounds={origin,origin};auto area=Exact::fromDouble(0);
        for(const auto& box:cells){selectedBounds.hi.x=std::max(selectedBounds.hi.x,box.hi.x);selectedBounds.hi.y=std::max(selectedBounds.hi.y,box.hi.y);
            area=Exact::add(area,Exact::multiply(Exact::subtract(Exact::fromDouble(box.hi.x),Exact::fromDouble(box.lo.x)),
                Exact::subtract(Exact::fromDouble(box.hi.y),Exact::fromDouble(box.lo.y))));}
        // Native selected leaves are disjoint. Exact area equality proves
        // this particular union is a complete closed rectangle, allowing
        // a convex-domain shortcut. Arbitrary public boxes use exact slabs.
        convexRectangle=Exact::sign(Exact::subtract(area,Exact::multiply(
            Exact::subtract(Exact::fromDouble(selectedBounds.hi.x),Exact::fromDouble(origin.x)),
            Exact::subtract(Exact::fromDouble(selectedBounds.hi.y),Exact::fromDouble(origin.y)))))==0;
        std::vector<int> parent(np);std::iota(parent.begin(),parent.end(),0);
        const auto root=[&](int p){int r=p;while(parent[r]!=r)r=parent[r];while(parent[p]!=p){const int next=parent[p];parent[p]=r;p=next;}return r;};
        for(int u=0;u<l.velocityCount();++u){const auto f=l.domain().face(mesh,l.nativeFace(u));faces.push_back(f);
            const int a=f.negativeCell<0?-1:l.pressureDof(f.negativeCell),b=f.positiveCell<0?-1:l.pressureDof(f.positiveCell);
            if(a>=0&&b>=0&&f.area>0)parent[root(b)]=root(a);
            if(a<0||b<0){const auto& e=f.horizontal?mesh.hEdges(f.edgeLabel):mesh.vEdges(f.edgeLabel);
                const auto x=mesh.verts(e.edgeToVert[0]).coord(),y=mesh.verts(e.edgeToVert[1]).coord();
                walls.push_back({{std::min(x.x,y.x),std::min(x.y,y.y)},{std::max(x.x,y.x),std::max(x.y,y.y)}});
                wallBins.insert(walls.back(),static_cast<int>(walls.size()-1));}}
        std::map<int,int> numbered;
        for(int p=0;p<np;++p){const int r=root(p);if(!numbered.count(r))numbered[r]=static_cast<int>(numbered.size());component.push_back(numbered[r]);}
        const auto& publication=phase.publication();
        const auto clock=color.timeLayer();
        for(std::size_t c=0;c<publication.topology().curves.size();++c)for(std::size_t g=0;g<publication.topology().curves[c].graphs().size();++g)
            for(const auto& q:publication.quadrature(c,g)){
                ++report.sourceNodes;const auto& t=q.sample.timeLayer;
                require(q.sample.geometryIdentity==color.geometryIdentity()&&t.step==clock.step&&t.begin==clock.begin&&t.end==clock.end&&t.position==clock.position,
                    "different surface geometry/full clock");
                require(std::isfinite(q.physicalWeight)&&q.physicalWeight>0&&std::isfinite(q.sample.curvature)&&
                    std::isfinite(q.sample.frame.point.x)&&std::isfinite(q.sample.frame.point.y),"invalid original CSF surface sample");
                const int p=l.pressureDof(q.sample.frame.nativeCellIndex);if(p<0)continue;
                points.push_back({q.sample.frame.point,q.physicalWeight,q.sample.curvature,p,level[p]});
                pointBins.insert({points.back().point,points.back().point},static_cast<int>(points.size()-1));++report.selectedNodes;}
    }
    bool nearWall(Vector2D p,double radius) const
    {
        bool found=false;wallBins.visit(surrounding(p,radius),[&](int id){
            const auto& b=walls[id];const double dx=std::max({0.,b.lo.x-p.x,p.x-b.hi.x}),dy=std::max({0.,b.lo.y-p.y,p.y-b.hi.y});
            found=found||std::hypot(dx,dy)<radius;});return found;
    }
    bool visible(Vector2D a,Vector2D b,int pressureComponent,CchCsfFaceSupportReport& report,std::size_t queryLimit)
    {
        const auto inside=[&](Vector2D p){return p.x>=selectedBounds.lo.x&&p.x<=selectedBounds.hi.x&&p.y>=selectedBounds.lo.y&&p.y<=selectedBounds.hi.y;};
        if(convexRectangle&&inside(a)&&inside(b))return true;
        scratch.clear();intervals.clear();
        if(++visitEpoch==0){std::fill(cellVisitStamp.begin(),cellVisitStamp.end(),0);visitEpoch=1;}
        cellBins.visit({{std::min(a.x,b.x),std::min(a.y,b.y)},{std::max(a.x,b.x),std::max(a.y,b.y)}},[&](int id){
            if(cellVisitStamp[id]==visitEpoch)return;
            require(report.visibilityCellQueries<queryLimit,"CSF visibility cell query limit");++report.visibilityCellQueries;
            cellVisitStamp[id]=visitEpoch;scratch.push_back(id);});
        const RayLine ray(a,b);
        for(int id:scratch){
            if(component[id]!=pressureComponent)continue;RayParameter low(0),high(1);
            if(ray.clip(cells[id],low,high))intervals.emplace_back(std::move(low),std::move(high));}
        return covered(intervals);
    }
};
}
struct CchBalancedCsfLoad::Data
{
    Data(CchPressureLoad p,CchRegularizedColor c,double s,double k)
        :pressure(std::move(p)),color(std::move(c)),sigma(s),curvature(k){}
    CchPressureLoad pressure;
    CchRegularizedColor color;
    double sigma,curvature;
    std::vector<double> gradient,rhs;
    bool capillary=false;
    CchCsfCurvatureConfig config;
    std::vector<double> averages;
    std::vector<CchCsfFaceSupportReport> support;
    CchCsfCurvatureReport report;
};
CchBalancedCsfLoad::CchBalancedCsfLoad(std::shared_ptr<const Data> d):data_(std::move(d)){}
CchBalancedCsfLoad CchBalancedCsfLoad::constantCurvature(const QuadtreeMesh& mesh,
    CchPressureLoad pressure,CchRegularizedColor color,double sigma,double curvature)
{
    require(std::isfinite(sigma)&&sigma>=0,"invalid surface tension");
    require(std::isfinite(curvature),"nonfinite constant curvature");
    require(pressure.layout().sameIdentity(color.phaseGeometry().layout()),"different pressure/color layout");
    require(&pressure.divergence().matrix()==&color.phaseGeometry().metrics().divergence().matrix(),
        "different authoritative divergence owner");
    color.validateFor(mesh,color.geometryIdentity(),color.timeLayer());
    const double coefficient=sigma*curvature;
    require(std::isfinite(coefficient)&&(coefficient!=0||sigma==0||curvature==0),
        "unrepresentable constant CSF coefficient");
    auto d=std::make_shared<Data>(std::move(pressure),std::move(color),sigma,curvature);
    const int count=d->pressure.layout().velocityCount();
    // Reuse the only authoritative gradient. No new stencil, area factor,
    // density factor, sharp alpha substitution, or virtual coarse-face split.
    d->gradient=d->pressure.apply(mesh,CchNativeValues(d->pressure.layout(),std::vector<double>(count,0),d->color.colors()));
    d->rhs.assign(count,0);
    if(sigma!=0&&curvature!=0)for(int u=0;u<count;++u){
        d->rhs[u]=coefficient*d->gradient[u];
        require(std::isfinite(d->rhs[u]),"nonfinite constant CSF momentum RHS");}
    return CchBalancedCsfLoad(std::move(d));
}
double CchBalancedCsfLoad::positiveCurvatureAverage(const std::vector<double>& weights,const std::vector<double>& values)
{
    require(!weights.empty()&&weights.size()==values.size(),"invalid curvature average dimensions");
    double ws=0,ks=0;std::size_t anchor=weights.size();bool positive=false,negative=false,constant=true;
    for(std::size_t i=0;i<weights.size();++i){require(std::isfinite(weights[i])&&weights[i]>=0&&std::isfinite(values[i]),"invalid positive curvature average input");
        if(weights[i]>0){ws=std::max(ws,weights[i]);ks=std::max(ks,std::abs(values[i]));if(anchor==weights.size())anchor=i;
            constant=constant&&values[i]==values[anchor];
            positive=positive||values[i]>0;negative=negative||values[i]<0;}}
    require(ws>0,"no positive curvature support");if(constant)return values[anchor];
    const long double base=static_cast<long double>(values[anchor])/ks;long double total=0,offset=0;
    for(std::size_t i=0;i<weights.size();++i)if(weights[i]>0){const long double w=static_cast<long double>(weights[i])/ws,
            k=static_cast<long double>(values[i])/ks;
        require(w>0&&(k!=0||values[i]==0),"unrepresentable normalized curvature average input");
        const long double delta=k-base,term=w*delta;
        require(term!=0||delta==0,"unrepresentable normalized curvature average moment");total+=w;offset+=term;}
    const double mean=static_cast<double>((base+offset/total)*ks);
    require(std::isfinite(mean)&&(mean!=0||positive==negative),"unrepresentable positive curvature average");return mean;
}
bool CchBalancedCsfLoad::selectedSegmentVisible(Vector2D a,Vector2D b,const std::vector<std::array<Vector2D,2>>& boxes)
{
    Exact::validateEnvironment();const RayLine ray(a,b);RayIntervals spans;
    for(const auto& box:boxes){require(std::isfinite(box[0].x)&&std::isfinite(box[0].y)&&std::isfinite(box[1].x)&&std::isfinite(box[1].y)&&
            box[0].x<=box[1].x&&box[0].y<=box[1].y,"invalid selected visibility box");
        RayParameter low(0),high(1);if(ray.clip({box[0],box[1]},low,high))spans.emplace_back(std::move(low),std::move(high));}
    return covered(spans);
}
CchBalancedCsfLoad CchBalancedCsfLoad::capillary(const QuadtreeMesh& mesh,CchPressureLoad pressure,
    CchRegularizedColor color,double sigma,CchCsfCurvatureConfig cfg)
{
    require(std::isfinite(sigma)&&sigma>=0,"invalid surface tension");
    require(cfg.maxCandidateQueriesPerFace>0&&cfg.maxVisibilityCellQueriesPerFace>0&&std::isfinite(cfg.relativeGradientRoundoff)&&cfg.relativeGradientRoundoff>0&&cfg.relativeGradientRoundoff<1,
        "invalid CSF curvature config");
    require(pressure.layout().sameIdentity(color.phaseGeometry().layout()),"different pressure/color layout");
    require(&pressure.divergence().matrix()==&color.phaseGeometry().metrics().divergence().matrix(),"different authoritative divergence owner");
    color.validateFor(mesh,color.geometryIdentity(),color.timeLayer());
    auto d=std::make_shared<Data>(std::move(pressure),std::move(color),sigma,0);d->capillary=true;d->config=cfg;
    const auto& l=d->pressure.layout();const int nu=l.velocityCount();
    d->gradient=d->pressure.apply(mesh,CchNativeValues(l,std::vector<double>(nu,0),d->color.colors()));
    d->rhs.assign(nu,0);d->averages.assign(nu,0);d->support.resize(nu);
    if(sigma==0){d->report.projectionSkippedSigmaZero=true;for(auto& r:d->support)r.status=CchCsfSupportStatus::SigmaZero;
        return CchBalancedCsfLoad(std::move(d));}
    ProjectionIndex index(mesh,d->color,d->report);std::vector<double> weights,curvatures;
    const auto& cells=d->color.cellReports();
    for(int u=0;u<nu;++u){auto& r=d->support[u];const auto& face=index.faces[u];
        const int a=face.negativeCell<0?-1:l.pressureDof(face.negativeCell),b=face.positiveCell<0?-1:l.pressureDof(face.positiveCell);
        const int p=a>=0?a:b;require(p>=0,"CSF face without native pressure cell");r.pressureComponent=index.component[p];
        r.center=face.center;r.boundary=a<0||b<0;
        const double diameter=std::max(a<0?0:index.diagonal[a],b<0?0:index.diagonal[b]);
        r.radius=d->color.config().epsilon+diameter;r.nearAxis=l.domain().axisymmetric()&&r.center.x<r.radius;
        r.nearBoundary=index.nearWall(r.center,r.radius);
        r.crossesAmr=a>=0&&b>=0&&index.level[a]!=index.level[b];
        long double scale=0,estimatedError=0;
        for(CchPressureLoad::Matrix::InnerIterator it(d->pressure.matrix(),u);it;++it){
            scale+=std::abs(it.value()*d->color.colors()[it.col()]);
            estimatedError+=std::abs(it.value())*cells[it.col()].absoluteErrorEstimate/d->color.phaseGeometry().metrics().cellVolumes()[it.col()];}
        r.gradientRoundoffScale=static_cast<double>(cfg.relativeGradientRoundoff*scale);
        r.colorGradientErrorEstimate=static_cast<double>(estimatedError);
        r.subResolutionGradient=d->gradient[u]!=0&&std::abs(d->gradient[u])<=r.gradientRoundoffScale;
        weights.clear();curvatures.clear();
        const Box supportBox=surrounding(r.center,r.radius);
        index.pointBins.visit(supportBox,[&](int id){
            require(r.candidateQueries<cfg.maxCandidateQueriesPerFace,"CSF support query limit");++r.candidateQueries;
            const auto& q=index.points[id];if(index.component[q.pressure]!=r.pressureComponent)return;
            const double distance=std::hypot(q.point.x-r.center.x,q.point.y-r.center.y);if(distance>=r.radius)return;
            if(!index.visible(r.center,q.point,r.pressureComponent,r,cfg.maxVisibilityCellQueriesPerFace)){++r.visibilityRejected;return;}
            const double hat=1-distance/r.radius;if(hat==0)return;
            const double w=q.weight*hat;require(std::isfinite(w)&&w>0,"unrepresentable positive CSF weight");
            weights.push_back(w);curvatures.push_back(q.curvature);r.totalWeight+=w;
            if(r.supportCount==0)r.curvatureMin=r.curvatureMax=q.curvature;
            else{r.curvatureMin=std::min(r.curvatureMin,q.curvature);r.curvatureMax=std::max(r.curvatureMax,q.curvature);}
            ++r.supportCount;r.maxSupportDistance=std::max(r.maxSupportDistance,distance);
            r.crossesAmr=r.crossesAmr||q.level!=index.level[p];});
        d->report.candidateQueries+=r.candidateQueries;d->report.visibilityRejected+=r.visibilityRejected;
        d->report.nearBoundaryFaces+=r.nearBoundary;d->report.nearAxisFaces+=r.nearAxis;d->report.amrFaces+=r.crossesAmr;
        if(r.supportCount){require(std::isfinite(r.totalWeight)&&r.totalWeight>0,"unrepresentable CSF support measure");
            d->averages[u]=positiveCurvatureAverage(weights,curvatures);r.curvatureDefined=true;r.status=CchCsfSupportStatus::Projected;++d->report.projectedFaces;
            const double coefficient=sigma*d->averages[u];require(std::isfinite(coefficient)&&(coefficient!=0||d->averages[u]==0),"unrepresentable projected CSF coefficient");
            d->rhs[u]=coefficient*d->gradient[u];
            require(std::isfinite(d->rhs[u])&&(d->rhs[u]!=0||coefficient==0||d->gradient[u]==0),"unrepresentable projected CSF momentum RHS");}
        else if(d->gradient[u]==0){r.status=CchCsfSupportStatus::ZeroGradient;}
        else if(r.boundary&&cells[p].transitionVolumeUpper==0){
            // C1 on a certified band-free prescribed wall is a pressure
            // reaction term, not evidence of a physical interface here.
            r.status=CchCsfSupportStatus::PureBoundaryReaction;++d->report.pureBoundaryRows;}
        else throw std::invalid_argument("[CchBalancedCsfLoad] unsupported active CSF face without local selected/visible curvature support, velocity="+std::to_string(u));
    }
    return CchBalancedCsfLoad(std::move(d));
}
const CchPressureLoad& CchBalancedCsfLoad::pressureLoad() const noexcept{return data_->pressure;}
const CchRegularizedColor& CchBalancedCsfLoad::color() const noexcept{return data_->color;}
const CchDofLayout& CchBalancedCsfLoad::layout() const noexcept{return data_->pressure.layout();}
double CchBalancedCsfLoad::surfaceTension() const noexcept{return data_->sigma;}
double CchBalancedCsfLoad::constantCurvature() const{require(!data_->capillary,"projected load has no prescribed constant curvature");return data_->curvature;}
bool CchBalancedCsfLoad::isCapillary() const noexcept{return data_->capillary;}
const CchCsfCurvatureConfig& CchBalancedCsfLoad::curvatureConfig() const{require(data_->capillary,"constant load has no curvature projection config");return data_->config;}
const std::vector<double>& CchBalancedCsfLoad::curvatureAverages() const{require(data_->capillary,"constant load has no projected curvature");return data_->averages;}
const std::vector<CchCsfFaceSupportReport>& CchBalancedCsfLoad::supportReports() const{require(data_->capillary,"constant load has no curvature support reports");return data_->support;}
const CchCsfCurvatureReport& CchBalancedCsfLoad::curvatureReport() const{require(data_->capillary,"constant load has no curvature projection report");return data_->report;}
const std::vector<double>& CchBalancedCsfLoad::colorGradient() const noexcept{return data_->gradient;}
const std::vector<double>& CchBalancedCsfLoad::momentumRhsLoad() const noexcept{return data_->rhs;}
Interface::HermiteGalsIdentity CchBalancedCsfLoad::geometryIdentity() const{return data_->color.geometryIdentity();}
Interface::HermiteGalsTimeLayer CchBalancedCsfLoad::timeLayer() const{return data_->color.timeLayer();}
void CchBalancedCsfLoad::validateFor(const QuadtreeMesh& mesh,const CchPressureLoad& pressure,
    Interface::HermiteGalsIdentity geometry,Interface::HermiteGalsTimeLayer clock) const
{
    require(&pressure.matrix()==&data_->pressure.matrix(),"different authoritative pressure owner");
    data_->color.validateFor(mesh,geometry,clock);
}
void CchBalancedCsfLoad::validateFor(const QuadtreeMesh& mesh,const CchPressureLoad& pressure,const CchRegularizedColor& color) const
{
    validateFor(mesh,pressure,color.geometryIdentity(),color.timeLayer());
    require(&color.colors()==&data_->color.colors(),"different authoritative color owner/epsilon");
    data_->color.validateFor(mesh,color.phaseGeometry());
}
}
