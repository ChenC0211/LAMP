#include "Modules/Interface/HermiteGalsPhaseQuadrature.h"

#include "Core/Geometry/LevelSetMappedSegment.h"
#include "Core/Geometry/LevelSetHermiteGeometry.h"
#include "Core/Mesh/CartesianMeasure.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>

namespace Quantum::Interface
{
namespace
{
constexpr std::array<double,8> gx{{-.96028985649753623168,-.79666647741362673959,
    -.52553240991632898582,-.18343464249564980494,.18343464249564980494,
    .52553240991632898582,.79666647741362673959,.96028985649753623168}};
constexpr std::array<double,8> gw{{.10122853629037625915,.22238103445337447054,
    .31370664587788728734,.36268378337836198297,.36268378337836198297,
    .31370664587788728734,.22238103445337447054,.10122853629037625915}};
constexpr std::array<double,4> iy{{-.86113631159405257522,-.33998104358485626480,
    .33998104358485626480,.86113631159405257522}};
constexpr std::array<double,4> iw{{.34785484513745385737,.65214515486254614263,
    .65214515486254614263,.34785484513745385737}};
using Moments=std::array<long double,90>; // 2 phases * 3 inverse-r weights * 15 P4 monomials.
void require(bool ok,const std::string& reason)
{if(!ok)throw std::invalid_argument("[Interface/HermiteGalsPhaseQuadrature] "+reason);}
struct Rule
{
    std::vector<HermiteGalsPhaseQuadraturePoint> points;
    Moments moments{};
    long double error = 0;
};
void append(Rule& to,Rule&& from)
{
    for(std::size_t i=0;i<to.moments.size();++i)to.moments[i]+=from.moments[i];
    to.error+=from.error;
    to.points.insert(to.points.end(),from.points.begin(),from.points.end());
}
class CellRule
{
public:
    CellRule(const LevelSetHermiteGeometry::Snapshot& scalar,
        std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> owner,
        const HermiteGalsPhaseVolumeCell& source,const HermiteGalsPhaseQuadratureConfig& config,
        bool radial) : scalar_(scalar),owner_(std::move(owner)),source_(source),config_(config),radial_(radial),
        bounds_(scalar.nativeCellBounds(source.cell)) {}
    std::size_t queries() const {return queries_;}
    Rule rectangle() // Exact polynomial moments; off-axis reciprocal moments use adaptive outer GL8.
    {
        if(radial_ && bounds_[0].x>0)return adaptive(bounds_[0].x,bounds_[1].x,config_.relativeMomentTolerance*source_.capacity/4,0,0);
        Rule result;
        for(std::size_t i=0;i<iy.size();++i) {
            const double x=bounds_[0].x+(bounds_[1].x-bounds_[0].x)*(1+iy[i])/2;
            const long double weight=(static_cast<long double>(bounds_[1].x)-bounds_[0].x)*iw[i]/2;
            slice(result,x,weight);
        }
        return result;
    }
    Rule adaptive(double lo,double hi,long double localBudget,int depth,int minimumDepth)
    {
        const double middle=lo+(hi-lo)/2;
        require(middle>lo && middle<hi,"unrepresentable outer subdivision at cell="+std::to_string(source_.cell));
        const auto coarse=outer(lo,hi);
        auto fine=outer(lo,middle);append(fine,outer(middle,hi));
        long double error=0;
        for(std::size_t k=0;k<fine.moments.size();++k)error=std::max(error,std::abs(fine.moments[k]-coarse.moments[k]));
        require(std::isfinite(error),"nonfinite moment estimate");
        if(depth>=minimumDepth && error<=localBudget){fine.error=error;return fine;}
        if(depth>=config_.maximumDepth) {
            std::ostringstream message;message.precision(17);
            message<<"outer quadrature depth limit at cell="<<source_.cell<<" interval=["<<lo<<','<<hi
                <<"] estimate="<<error<<" budget="<<localBudget;
            require(false,message.str());
        }
        auto result=adaptive(lo,middle,localBudget,depth+1,minimumDepth);
        append(result,adaptive(middle,hi,localBudget,depth+1,minimumDepth));return result;
    }
private:
    void point(Rule& result,double x,long double y,long double weight,bool liquid)
    {
        const Vector2D p{x,static_cast<double>(y)};
        require(p.x>=bounds_[0].x && p.x<=bounds_[1].x && p.y>=bounds_[0].y && p.y<=bounds_[1].y,
            "quadrature point outside its native cell");
        require(scalar_.nativeQuarticBoxSign(source_.cell,p,p)==(liquid?-1:1),
            "unresolved phase point at cell="+std::to_string(source_.cell));
        const double physical=static_cast<double>(weight*(radial_?CartesianMeasure::radialWeight(x):1.L));
        require(std::isfinite(physical) && physical>0,"nonpositive/unrepresentable phase weight");
        require(result.points.size()<config_.maximumPoints,"phase point work limit");
        result.points.push_back({source_.cell,p,physical,liquid});
        // Normalize monomials in the ORIGINAL cell, not in refined slices.
        const long double nx=(static_cast<long double>(p.x)-bounds_[0].x)/(bounds_[1].x-bounds_[0].x);
        const long double ny=(static_cast<long double>(p.y)-bounds_[0].y)/(bounds_[1].y-bounds_[0].y);
        std::array<long double,5> xp{{1,0,0,0,0}},yp{{1,0,0,0,0}};
        for(int i=1;i<5;++i){xp[i]=xp[i-1]*nx;yp[i]=yp[i-1]*ny;}
        const int count=radial_ && bounds_[0].x>0?3:1;
        const long double ratio=count==3?(static_cast<long double>(bounds_[0].x)+bounds_[1].x)/(2*p.x):1;
        long double factor=1;
        for(int inverse=0;inverse<count;++inverse) {
            int monomial=0;
            for(int degree=0;degree<=4;++degree)for(int i=0;i<=degree;++i)
                result.moments[(liquid?0:45)+15*inverse+monomial++]+=physical*factor*xp[i]*yp[degree-i];
            factor*=ratio;
        }
    }
    void slice(Rule& result,double x,long double outerWeight)
    {
        require(++queries_<=config_.maximumSliceQueriesPerCell,"slice work limit at cell="+std::to_string(source_.cell));
        const auto interval=[&](bool liquid,bool fromEnd,long double offset,long double length) {
            const long double height=static_cast<long double>(bounds_[1].y)-bounds_[0].y;
            for(std::size_t i=0;i<iy.size();++i) {
                const long double distance=height*(offset+length*((1.L+iy[i])/2));
                const long double y=fromEnd?bounds_[1].y-distance:bounds_[0].y+distance;
                point(result,x,y,outerWeight*height*length*iw[i]/2,liquid);
            }
        };
        if(source_.certifiedPhase!=0){interval(source_.certifiedPhase<0,false,0,1);return;}
        const LevelSetMappedSegment segment(owner_,{x,bounds_[0].y},{x,bounds_[1].y});
        for(const auto& range:segment.ranges())interval(range.negative,range.fromEnd,range.offset,range.length);
    }
    Rule outer(double lo,double hi)
    {
        Rule result;
        for(std::size_t i=0;i<gx.size();++i) {
            // Endpoint-relative embedding avoids 1-t when an interval is short.
            const bool upper=gx[i]>0;const long double fraction=(1.L-std::abs(gx[i]))/2;
            const long double width=static_cast<long double>(hi)-lo;
            const double x=static_cast<double>(upper?hi-width*fraction:lo+width*fraction);
            require(x>lo && x<hi,"unrepresentable outer interior at cell="+std::to_string(source_.cell));
            slice(result,x,width*gw[i]/2);
        }
        return result;
    }
    const LevelSetHermiteGeometry::Snapshot& scalar_;
    std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> owner_;
    const HermiteGalsPhaseVolumeCell& source_;
    const HermiteGalsPhaseQuadratureConfig& config_;
    bool radial_;
    std::array<Vector2D,2> bounds_;
    std::size_t queries_=0;
};
}
struct HermiteGalsPhaseQuadrature::Data
{
    explicit Data(HermiteGalsPublication p) : publication(std::move(p)) {}
    HermiteGalsPublication publication;
    std::vector<HermiteGalsPhaseQuadraturePoint> points;
    std::vector<HermiteGalsPhaseQuadratureCell> cells;
    HermiteGalsPhaseQuadratureReport report;
};
HermiteGalsPhaseQuadrature::HermiteGalsPhaseQuadrature(std::shared_ptr<const Data> d):data_(std::move(d)){}
HermiteGalsPhaseQuadrature HermiteGalsPhaseQuadrature::build(const QuadtreeMesh& mesh,HermiteGalsPublication publication,
    const HermiteGalsPhaseQuadratureConfig& config)
{
    publication.validateMesh(mesh);
    require(publication.state().primary().transport==HermiteGalsTransport::FirstJet,
        "requires FirstJet publication");
    require(std::isfinite(config.relativeMomentTolerance) && config.relativeMomentTolerance>0 && config.relativeMomentTolerance<=1e-6 &&
        std::isfinite(config.relativeVolumeTolerance) && config.relativeVolumeTolerance>0 && config.relativeVolumeTolerance<=1e-6,
        "invalid quadrature tolerance");
    require(config.maximumDepth>=0 && config.maximumDepth<=40 && config.maximumSliceQueriesPerCell>0 && config.maximumPoints>0,
        "invalid quadrature work limit");
    auto data=std::make_shared<Data>(std::move(publication));
    const auto geometry=std::make_shared<const HermiteGalsGeometryView>(data->publication.state().geometry());
    const std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> scalar(geometry,&geometry->snapshot());
    std::vector<std::vector<double>> seeds(scalar->nativeCellCount());
    for(const auto& curve:data->publication.topology().curves)for(std::size_t g=0;g<curve.graphs().size();++g) {
        for(double fraction:{0.,.5,1.}) {
            const auto frame=curve.frameAt(static_cast<int>(g),fraction);
            seeds[frame.nativeCellIndex].push_back(frame.point.x); // Numerical subdivisions only; never a phase test.
        }
        if(!curve.graphs()[g].coordinateIsX)for(int half=0;half<2;++half) {
            double lo=.5*half,hi=.5*(half+1);
            auto first=curve.frameAt(static_cast<int>(g),lo),last=curve.frameAt(static_cast<int>(g),hi);
            // Refine a bracketed x-projection turning point of this SAME
            // regular y graph. These are numerical integration breakpoints,
            // not a phase/topology certificate or a replacement curve.
            if(!((first.normal.y<0 && last.normal.y>0) || (first.normal.y>0 && last.normal.y<0)))continue;
            for(int iteration=0;iteration<40;++iteration) {
                const double middle=lo+(hi-lo)/2;
                require(middle>lo && middle<hi,"unrepresentable projection refinement");
                const auto frame=curve.frameAt(static_cast<int>(g),middle);
                if(frame.normal.y==0){first=last=frame;break;}
                if((frame.normal.y>0)==(first.normal.y>0)){lo=middle;first=frame;}
                else{hi=middle;last=frame;}
            }
            seeds[curve.graphs()[g].nativeCellIndex].push_back(std::min(first.point.x,last.point.x));
        }
    }
    for(const auto& source:data->publication.volumes().cells()) {
        const auto bounds=scalar->nativeCellBounds(source.cell);
        std::vector<double> cuts{bounds[0].x,bounds[1].x};
        for(double x:seeds[source.cell])if(x>bounds[0].x && x<bounds[1].x)cuts.push_back(x);
        std::sort(cuts.begin(),cuts.end());cuts.erase(std::unique(cuts.begin(),cuts.end()),cuts.end());
        const long double targetLiquid=source.certifiedPhase<0?source.capacity:source.certifiedPhase>0?0:source.rawLiquid;
        const long double targetGas=source.capacity-targetLiquid;
        CellRule integrator(*scalar,scalar,source,config,data->publication.domain().axisymmetric());
        Rule rule;bool accepted=false;long double lastVolumeError=0,lastMomentError=0;
        for(int retry=0;retry<5 && !accepted;++retry) {
            rule=Rule{};
            if(source.certifiedPhase!=0)rule=integrator.rectangle();
            else for(std::size_t i=1;i<cuts.size();++i) {
                const long double portion=(static_cast<long double>(cuts[i])-cuts[i-1])/(bounds[1].x-bounds[0].x);
                try {append(rule,integrator.adaptive(cuts[i-1],cuts[i],config.relativeMomentTolerance*source.capacity*portion/
                    (4*std::pow(10.L,retry)),0,retry));}
                catch(const std::exception& e) {
                    std::ostringstream message;message.precision(17);message<<e.what()<<" retry="<<retry
                        <<" previous volume difference="<<lastVolumeError<<" previous moment estimate="<<lastMomentError;
                    require(false,message.str());
                }
            }
            const long double volumeError=std::max(std::abs(rule.moments[0]-targetLiquid),std::abs(rule.moments[45]-targetGas))/source.capacity;
            lastVolumeError=volumeError;lastMomentError=rule.error/source.capacity;
            accepted=std::isfinite(volumeError) && volumeError<=config.relativeVolumeTolerance &&
                rule.error<=config.relativeMomentTolerance*source.capacity &&
                (targetLiquid==0?rule.moments[0]==0:rule.moments[0]>0) && (targetGas==0?rule.moments[45]==0:rule.moments[45]>0);
            if(accepted) {
                data->cells.push_back({source.cell,rule.moments[0],rule.moments[45],volumeError,rule.error/source.capacity});
                data->report.maxRelativeVolumeDifference=std::max(data->report.maxRelativeVolumeDifference,volumeError);
                data->report.maxRelativeMomentEstimate=std::max(data->report.maxRelativeMomentEstimate,rule.error/source.capacity);
            }
        }
        require(accepted,"phase volume/moment agreement failed at cell="+std::to_string(source.cell));
        require(rule.points.size()<=config.maximumPoints-data->points.size(),"total phase point work limit");
        data->points.insert(data->points.end(),rule.points.begin(),rule.points.end());
        data->report.sliceQueries+=integrator.queries();
    }
    data->report.points=data->points.size();return HermiteGalsPhaseQuadrature(std::move(data));
}
const HermiteGalsPublication& HermiteGalsPhaseQuadrature::publication() const noexcept{return data_->publication;}
HermiteGalsIdentity HermiteGalsPhaseQuadrature::geometryIdentity() const{return publication().state().identity();}
HermiteGalsTimeLayer HermiteGalsPhaseQuadrature::timeLayer() const{return publication().state().timeLayer();}
const std::vector<HermiteGalsPhaseQuadraturePoint>& HermiteGalsPhaseQuadrature::points() const noexcept{return data_->points;}
const std::vector<HermiteGalsPhaseQuadratureCell>& HermiteGalsPhaseQuadrature::cells() const noexcept{return data_->cells;}
const HermiteGalsPhaseQuadratureReport& HermiteGalsPhaseQuadrature::report() const noexcept{return data_->report;}
void HermiteGalsPhaseQuadrature::validateMesh(const QuadtreeMesh& mesh) const{publication().validateMesh(mesh);}
}
