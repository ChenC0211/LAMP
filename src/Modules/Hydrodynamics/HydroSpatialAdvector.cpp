#include "Modules/Hydrodynamics/HydroSpatialAdvector.h"
#include "Modules/Hydrodynamics/Detail/HydroRegularAdvectorGeometry.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/CartesianMeasure.h"
#include "Core/Mesh/UniformLeafGrid.h"
#include <algorithm>
#include <atomic>
#include <cfenv>
#include <cmath>
#include <map>
#include <limits>
#include <optional>
#include <mutex>
#include <sstream>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#if defined(_M_X64) || defined(__SSE2__)
#include <xmmintrin.h>
#endif

namespace LAMP::Hydrodynamics
{
namespace
{
void require(bool ok, const char* why)
{ if (!ok) throw std::invalid_argument(std::string("[HydroSpatialAdvector] ")+why); }
void fp()
{
    require(std::fegetround()==FE_TONEAREST,"round-to-nearest required; mode not changed");
#if defined(_M_X64) || defined(__SSE2__)
    require((_mm_getcsr() & ((1U<<15)|(1U<<6)))==0,"gradual underflow required; mode not changed");
#endif
}
long double c(Vector2D p,int axis) { return axis==0?p.x:p.y; }
Vector2D point(int axis,double plane,double t) { return axis==0?Vector2D{plane,t}:Vector2D{t,plane}; }
void sortUnique(std::vector<double>& v) { std::sort(v.begin(),v.end());v.erase(std::unique(v.begin(),v.end()),v.end()); }
bool sameTime(Interface::HermiteGalsTimeLayer a,Interface::HermiteGalsTimeLayer b)
{ a.time();b.time();return a.step==b.step && a.begin==b.begin && a.end==b.end && a.position==b.position; }
std::atomic<std::uint64_t> nextIdentity{1};
struct ModalBubbleJet {double v=0,d=0,dd=0;};
ModalBubbleJet modalBubble(double coordinate,double lower,double upper,bool evenAxis)
{
    const double h=upper-lower,s=evenAxis?coordinate*coordinate/(h*h):(coordinate-lower)/h;
    const double ds=evenAxis?2*coordinate/(h*h):1/h,dds=evenAxis?2/(h*h):0;
    const double g=4*s*(1-s),dg=4-8*s,value=g*g*g,first=3*g*g*dg,second=6*g*dg*dg-24*g*g;
    return {value,first*ds,second*ds*ds+first*dds};
}
}
struct HydroAdvectorReconstruction::Preparation
{
    explicit Preparation(HydroMomentumGeometry g):geometry(std::move(g)){}
    HydroMomentumGeometry geometry;
    std::vector<HydroAdvectorAtom> atoms;
    std::vector<HydroAdvectorPiece> pieces;
    std::vector<std::vector<int>> primal,dual;
    HydroAdvectorPreparationReport report;
    std::shared_ptr<const Detail::HydroRegularAdvectorGeometry> regular;
};
struct HydroAdvectorReconstruction::Data
{
    explicit Data(std::shared_ptr<const Preparation> p):preparation(std::move(p)),
        geometry(preparation->geometry),report(preparation->report),regular(preparation->regular){}
    std::shared_ptr<const Preparation> preparation;
    HydroMomentumGeometry geometry;
    HydroAdvectorPreparationReport report;
    std::shared_ptr<const Detail::HydroRegularAdvectorGeometry> regular;
    std::optional<HydroPhaseGeometry> phase;
    std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> scalar;
    std::vector<HydroTangentialMode> modes;
    std::vector<int> modeByPressure;
};
struct HydroSpatialAdvector::Data
{
    Data(HydroAdvectorReconstruction r,HydroNativeValues v,Interface::HermiteGalsIdentity g,
        Interface::HermiteGalsTimeLayer t,HydroAdvectorIdentity i)
        : reconstruction(std::move(r)),values(std::move(v)),geometry(g),clock(t),identity(i){}
    HydroAdvectorReconstruction reconstruction;
    HydroNativeValues values;
    Interface::HermiteGalsIdentity geometry;
    Interface::HermiteGalsTimeLayer clock;
    HydroAdvectorIdentity identity;
    std::optional<HydroAdvectorBoundary> boundary;
    Detail::HydroRegularCoefficients coefficients;
    mutable std::once_flag curvePreparation;
    mutable std::vector<Detail::HydroRegularPiecePolynomial> curvePolynomials;
    // The virtual grid is a geometric chart, not a request to convert every
    // quarter piece whenever a curve touches one of them. Cache only actual
    // curve queries; old native sources retain their original batch path.
    struct VirtualCurvePolynomial
    {
        std::once_flag ready;
        std::optional<Detail::HydroRegularPiecePolynomial> value;
    };
    mutable std::mutex virtualCurveMutex;
    mutable std::map<int,std::shared_ptr<VirtualCurvePolynomial>> virtualCurvePolynomials;
    const Detail::HydroRegularPiecePolynomial& curvePolynomial(int piece) const
    {
        const auto& regular=reconstruction.regularGeometry();const auto& pieces=reconstruction.pieces();
        if(!regular.usesVirtualCompatibleTensor()){
            std::call_once(curvePreparation,[&]{curvePolynomials=regular.piecePolynomials(coefficients,pieces);});
            return curvePolynomials[piece];
        }
        std::shared_ptr<VirtualCurvePolynomial> cached;
        {std::lock_guard<std::mutex> lock(virtualCurveMutex);auto& slot=virtualCurvePolynomials[piece];
            if(!slot)slot=std::make_shared<VirtualCurvePolynomial>();cached=slot;}
        std::call_once(cached->ready,[&]{auto selected=regular.piecePolynomials(coefficients,{pieces[piece]});cached->value=std::move(selected.front());});
        return *cached->value;
    }
    std::vector<long double> divergence;
    std::vector<HydroNormalPolynomial> polynomials;
    std::vector<double> tangential;
};
struct HydroAdvectorBoundary::Data
{
    Data(HydroMomentumGeometry g,std::vector<int> v,bool complete)
        :geometry(std::move(g)),velocities(std::move(v)),completeLines(complete){}
    HydroMomentumGeometry geometry;
    std::vector<int> velocities;
    bool completeLines=true;
    // Preparation retains geometry only and never this boundary, so the
    // accepted boundary can keep it across steps without an ownership cycle.
    mutable std::mutex reconstructionMutex;
    mutable std::shared_ptr<const HydroAdvectorReconstruction::Preparation> reconstructionPreparation;
    mutable std::shared_ptr<const HydroAdvectorReconstruction::Preparation> unboundReconstructionPreparation;
};
HydroAdvectorBoundary::HydroAdvectorBoundary(std::shared_ptr<const Data> d):data_(std::move(d)){}
HydroAdvectorBoundary HydroAdvectorBoundary::impermeable(const QuadtreeMesh& mesh,HydroMomentumGeometry geometry,const std::vector<int>& dofs)
{
    fp();geometry.validateMesh(mesh);auto ordered=dofs;std::sort(ordered.begin(),ordered.end());
    require(std::adjacent_find(ordered.begin(),ordered.end())==ordered.end(),"duplicate impermeable face");
    const auto& l=geometry.layout();for(int u:ordered){require(u>=0 && u<l.velocityCount(),"impermeable velocity index out of range");
        require(l.domain().face(mesh,l.nativeFace(u)).outwardSign!=0,"impermeable law requires selected-domain boundary face");}
    // A tensor normal trace can be zero on a whole boundary line; one
    // isolated zero interval mean cannot constrain that entire interval.
    // Compute the geometry-only qualification here once, never per stencil.
    std::map<std::pair<int,double>,std::array<std::size_t,2>> lines;
    for(int u=0;u<l.velocityCount();++u){const auto face=l.domain().face(mesh,l.nativeFace(u));
        if(face.outwardSign==0||geometry.metrics().faceMeasures()[u]==0)continue;
        const int normal=face.horizontal?1:0;auto& count=lines[{normal,static_cast<double>(c(face.center,normal))}];++count[0];
        if(std::binary_search(ordered.begin(),ordered.end(),u))++count[1];
    }
    bool complete=true;for(const auto& line:lines)
        complete=complete&&(line.second[1]==0||line.second[1]==line.second[0]);
    return HydroAdvectorBoundary(std::make_shared<Data>(std::move(geometry),std::move(ordered),complete));
}
const HydroMomentumGeometry& HydroAdvectorBoundary::momentumGeometry() const noexcept{return data_->geometry;}
const std::vector<int>& HydroAdvectorBoundary::velocityDofs() const noexcept{return data_->velocities;}
bool HydroAdvectorBoundary::completeBoundaryLines() const noexcept{return data_->completeLines;}
bool HydroAdvectorBoundary::sameIdentity(const HydroAdvectorBoundary& other) const noexcept{return data_==other.data_;}
void HydroAdvectorBoundary::validateFor(const QuadtreeMesh& mesh,const HydroMomentumGeometry& expected) const
{fp();momentumGeometry().validateMesh(mesh);require(momentumGeometry().sameIdentity(expected),"different impermeable momentum geometry instance");}
HydroAdvectorReconstruction::HydroAdvectorReconstruction(std::shared_ptr<const Data> d):data_(std::move(d)){}
const HydroMomentumGeometry& HydroAdvectorReconstruction::momentumGeometry() const noexcept{return data_->geometry;}
const std::vector<HydroAdvectorAtom>& HydroAdvectorReconstruction::atoms() const noexcept{return data_->preparation->atoms;}
const std::vector<HydroAdvectorPiece>& HydroAdvectorReconstruction::pieces() const noexcept{return data_->preparation->pieces;}
const std::vector<std::vector<int>>& HydroAdvectorReconstruction::primalAtoms() const noexcept{return data_->preparation->primal;}
const std::vector<std::vector<int>>& HydroAdvectorReconstruction::dualAtoms() const noexcept{return data_->preparation->dual;}
const HydroAdvectorPreparationReport& HydroAdvectorReconstruction::report() const noexcept{return data_->report;}
const Detail::HydroRegularAdvectorGeometry& HydroAdvectorReconstruction::regularGeometry() const noexcept{return *data_->regular;}
bool HydroAdvectorReconstruction::sameIdentity(const HydroAdvectorReconstruction& r) const noexcept{return data_==r.data_;}
const HydroPhaseGeometry* HydroAdvectorReconstruction::phaseGeometry() const noexcept
{return data_->phase?&*data_->phase:nullptr;}
const std::vector<HydroTangentialMode>& HydroAdvectorReconstruction::tangentialModes() const noexcept{return data_->modes;}
int HydroAdvectorReconstruction::tangentialModeForPressure(int p) const
{
    require(p>=0&&p<momentumGeometry().layout().pressureCount(),"modal pressure index out of range");
    return data_->modeByPressure.empty()?-1:data_->modeByPressure[p];
}
HydroAdvectorReconstruction HydroAdvectorReconstruction::withTangentialModes(const QuadtreeMesh& mesh,HydroPhaseGeometry phase) const
{
    fp();require(!phaseGeometry(),"tangential source is already phase bound");
    phase.validateFor(mesh,phase.geometryIdentity(),phase.timeLayer());
    require(phase.layout().sameIdentity(momentumGeometry().layout()),"different tangential phase layout");
    require(&phase.metrics().divergence().matrix()==&momentumGeometry().metrics().divergence().matrix(),
        "different tangential divergence owner");
    auto d=std::make_shared<Data>(*data_);d->phase=std::move(phase);
    const auto& geometry=d->phase->publication().state().geometry();
    d->scalar=std::shared_ptr<const LevelSetHermiteGeometry::Snapshot>(
        std::make_shared<const Interface::HermiteGalsGeometryView>(geometry),&geometry.snapshot());
    const auto& layout=momentumGeometry().layout();const int np=layout.pressureCount();
    std::vector<long double> measures(np),scales(np);d->modeByPressure.assign(np,-1);
    const auto& publication=d->phase->publication();
    for(std::size_t curve=0;curve<publication.topology().curves.size();++curve)
        for(std::size_t graph=0;graph<publication.topology().curves[curve].graphs().size();++graph)
            for(const auto& q:publication.quadrature(curve,graph)){
                const auto& frame=q.sample.frame;const int p=layout.pressureDof(frame.nativeCellIndex);
                if(p<0||q.physicalWeight==0)continue;
                const auto jet=d->scalar->scalarJetInNativeCell(frame.nativeCellIndex,frame.point);
                const double norm=std::hypot(jet.dx,jet.dy);
                require(std::isfinite(norm)&&norm>0,"unresolved tangential gradient scale");
                measures[p]+=q.physicalWeight;scales[p]+=q.physicalWeight*norm;
            }
    for(int p=0;p<np;++p){const auto& volume=d->phase->cellVolumes()[p];
        if(!(volume.liquid>0&&volume.gas>0))continue;
        require(measures[p]>0,"cut cell has no tangential interface quadrature");
        const auto& cell=regularGeometry().cells()[p];const double scale=static_cast<double>(scales[p]/measures[p]);
        require(std::isfinite(scale)&&scale>0,"invalid tangential gradient scale");
        d->modeByPressure[p]=static_cast<int>(d->modes.size());
        d->modes.push_back({p,layout.nativeCell(p),cell.lower,cell.upper,scale});
    }
    return HydroAdvectorReconstruction(std::move(d));
}
LevelSetBoundaryTrace::VelocityJet HydroAdvectorReconstruction::tangentialModeJet(int index,Vector2D point,bool liquid) const
{
    fp();require(index>=0&&static_cast<std::size_t>(index)<data_->modes.size(),"tangential mode index out of range");
    const auto& mode=data_->modes[index];
    if(!(std::isfinite(point.x)&&std::isfinite(point.y)&&point.x>=mode.lower.x&&point.x<=mode.upper.x
        &&point.y>=mode.lower.y&&point.y<=mode.upper.y)){
        std::ostringstream message;message.precision(17);
        message<<"[HydroSpatialAdvector] point outside tangential mode cell: mode="<<index<<", pressure="<<mode.pressure
            <<", point=("<<point.x<<','<<point.y<<"), box=("<<mode.lower.x<<','<<mode.lower.y<<")-("
            <<mode.upper.x<<','<<mode.upper.y<<')';throw std::invalid_argument(message.str());}
    if(liquid)return {};
    const auto h=mode.upper-mode.lower;const bool rz=momentumGeometry().layout().domain().axisymmetric();
    const auto x=modalBubble(point.x,mode.lower.x,mode.upper.x,rz&&mode.lower.x==0);
    const auto y=modalBubble(point.y,mode.lower.y,mode.upper.y,false);
    const double b=x.v*y.v,bx=x.d*y.v,by=x.v*y.d,bxx=x.dd*y.v,bxy=x.d*y.d,byy=x.v*y.dd;
    const auto phi=data_->scalar->scalarJetInNativeCell(mode.nativeCell,point);
    const double scale=mode.gradientScale*(rz?h.x:1);
    const double a=phi.value*b/scale,ax=(phi.dx*b+phi.value*bx)/scale,ay=(phi.dy*b+phi.value*by)/scale;
    const double axx=(phi.dxx*b+2*phi.dx*bx+phi.value*bxx)/scale;
    const double axy=(phi.dxy*b+phi.dx*by+phi.dy*bx+phi.value*bxy)/scale;
    const double ayy=(phi.dyy*b+2*phi.dy*by+phi.value*byy)/scale;
    LevelSetBoundaryTrace::VelocityJet result;
    if(rz){const double r=point.x;result.value={r*ay,-2*a-r*ax};
        result.derivative={ay+r*axy,r*ayy,-3*ax-r*axx,-2*ay-r*axy};}
    else{result.value={ay,-ax};result.derivative={axy,ayy,-axx,-axy};}
    require(std::isfinite(result.value.x)&&std::isfinite(result.value.y)&&std::isfinite(result.derivative.xx)
        &&std::isfinite(result.derivative.xy)&&std::isfinite(result.derivative.yx)&&std::isfinite(result.derivative.yy),
        "nonfinite tangential mode jet");return result;
}
long double HydroAdvectorReconstruction::tangentialModeStreamfunction(int index,Vector2D point,bool liquid) const
{
    require(index>=0&&static_cast<std::size_t>(index)<data_->modes.size(),"tangential mode index out of range");
    const auto& mode=data_->modes[index];
    require(point.x>=mode.lower.x&&point.x<=mode.upper.x&&point.y>=mode.lower.y&&point.y<=mode.upper.y,
        "streamfunction point outside tangential cell");
    if(liquid)return 0;
    const bool rz=momentumGeometry().layout().domain().axisymmetric();
    const auto x=modalBubble(point.x,mode.lower.x,mode.upper.x,rz&&mode.lower.x==0);
    const auto y=modalBubble(point.y,mode.lower.y,mode.upper.y,false);
    const auto phi=data_->scalar->scalarJetInNativeCell(mode.nativeCell,point);
    const long double primitive=static_cast<long double>(phi.value)*x.v*y.v/mode.gradientScale;
    const long double value=primitive*(rz?point.x*point.x/(mode.upper.x-mode.lower.x):1.);
    require(std::isfinite(value),"nonfinite tangential streamfunction");return value;
}

HydroAdvectorReconstruction HydroAdvectorReconstruction::prepare(const QuadtreeMesh& mesh,HydroMomentumGeometry geometry)
{return prepareBound(mesh,std::move(geometry),nullptr,nullptr);}
HydroAdvectorReconstruction HydroAdvectorReconstruction::prepare(const QuadtreeMesh& mesh,HydroMomentumGeometry geometry,
    const HydroAdvectorBoundary& boundary)
{return prepareBound(mesh,std::move(geometry),&boundary,&boundary);}
HydroAdvectorReconstruction HydroAdvectorReconstruction::prepareCached(const QuadtreeMesh& mesh,HydroMomentumGeometry geometry,
    const HydroAdvectorBoundary& boundary)
{return prepareBound(mesh,std::move(geometry),&boundary,&boundary);}
HydroAdvectorReconstruction HydroAdvectorReconstruction::prepareUnboundCached(const QuadtreeMesh& mesh,HydroMomentumGeometry geometry,
    const HydroAdvectorBoundary& cacheOwner)
{return prepareBound(mesh,std::move(geometry),nullptr,&cacheOwner);}
HydroAdvectorReconstruction HydroAdvectorReconstruction::prepareBound(const QuadtreeMesh& mesh,HydroMomentumGeometry geometry,
    const HydroAdvectorBoundary* boundary,const HydroAdvectorBoundary* cacheOwner)
{
    fp();geometry.validateMesh(mesh);
    if(boundary&&!boundary->momentumGeometry().sameIdentity(geometry))boundary->validateFor(mesh,geometry);
    if(cacheOwner&&cacheOwner!=boundary&&!cacheOwner->momentumGeometry().sameIdentity(geometry))cacheOwner->validateFor(mesh,geometry);
    std::unique_lock<std::mutex> cacheLock;
    if(cacheOwner){
        cacheLock=std::unique_lock<std::mutex>(cacheOwner->data_->reconstructionMutex);
        const auto& retained=boundary?cacheOwner->data_->reconstructionPreparation:cacheOwner->data_->unboundReconstructionPreparation;
        if(retained){
            require(retained->geometry.sameIdentity(geometry),"retained preparation has a different momentum owner");
            return HydroAdvectorReconstruction(std::make_shared<Data>(retained));
        }
    }
    auto d=std::make_shared<Preparation>(std::move(geometry));
    d->regular=std::make_shared<const Detail::HydroRegularAdvectorGeometry>(Detail::HydroRegularAdvectorGeometry::prepare(mesh,d->geometry,boundary));
    const auto& l=d->geometry.layout();const bool rz=l.domain().axisymmetric();
    const auto& cells=d->regular->cells();const auto& faces=d->regular->faces();
    const auto& prepared=d->regular->report();
    d->report.svdCount=prepared.svdCount;d->report.weights=prepared.weights;d->report.rankDeficient=prepared.rankDeficient;
    for(const auto& row:d->regular->rows()) {
        if(row.priorityRank<row.priorityDimension)++d->report.tangentUnidentified;
        d->report.maxTangentDeficit=std::max(d->report.maxTangentDeficit,row.priorityDeficit);
    }
    for(std::size_t p=0;p<cells.size();++p) {
        const auto& knots=cells[p].supportKnots;
        for(std::size_t x=1;x<knots[0].size();++x)for(std::size_t y=1;y<knots[1].size();++y)
            d->pieces.push_back({static_cast<int>(p),{knots[0][x-1],knots[1][y-1]},{knots[0][x],knots[1][y]}});
    }
    d->report.cellPieces=d->pieces.size();
    const auto owner=[&](int u){const auto& f=faces[u];return f.negativePressure>=0?f.negativePressure:f.positivePressure;};
    struct Input{int primal=-1,dual=-1,sourcePrimal=-1,sourceCell=-1;double lo=0,hi=0;
        int normal=0;double plane=0;std::size_t ordinal=0;};
    std::vector<Input> inputs;
    require(d->geometry.segments().size()<=inputs.max_size()-static_cast<std::size_t>(l.velocityCount()),"atom incidence overflow");
    inputs.reserve(static_cast<std::size_t>(l.velocityCount())+d->geometry.segments().size());
    const auto append=[&](int normal,Vector2D a,Vector2D b,int primal,int dual,int sourcePrimal,int sourceCell) {
        const double plane=static_cast<double>(c(a,normal));
        const double lo=static_cast<double>(c(a,1-normal)),hi=static_cast<double>(c(b,1-normal));
        require(lo<hi,"nonpositive actual segment");inputs.push_back({primal,dual,sourcePrimal,sourceCell,lo,hi,normal,plane,inputs.size()});
    };
    d->primal.resize(l.velocityCount());d->dual.resize(d->geometry.segments().size());
    for(int u=0;u<l.velocityCount();++u){const auto& f=faces[u];append(f.normalAxis,f.lower,f.upper,u,-1,u,-1);}
    for(std::size_t i=0;i<d->geometry.segments().size();++i){const auto& s=d->geometry.segments()[i];append(s.normalAxis,s.begin,s.end,-1,static_cast<int>(i),s.primalVelocity,s.ownerPressure);}
    // Equal numeric planes retain the first encountered key, including its zero sign.
    // Recreate each old input sequence before applying the original lo-only sort.
    std::sort(inputs.begin(),inputs.end(),[](const Input& a,const Input& b){
        const auto ak=std::make_pair(a.normal,a.plane),bk=std::make_pair(b.normal,b.plane);
        return ak<bk||(!(bk<ak)&&a.ordinal<b.ordinal);});
    for(std::size_t begin=0;begin<inputs.size();){std::size_t end=begin+1;
        const int normal=inputs[begin].normal;const double plane=inputs[begin].plane;
        while(end<inputs.size()&&inputs[end].normal==normal&&inputs[end].plane==plane)++end;
        std::vector<double> knots;
        for(auto i=begin;i<end;++i){const auto& v=inputs[i];knots.push_back(v.lo);knots.push_back(v.hi);
            const int pressure=v.sourcePrimal>=0?owner(v.sourcePrimal):v.sourceCell;
            for(double knot:cells[pressure].supportKnots[1-normal])if(v.lo<knot && knot<v.hi)knots.push_back(knot);}
        sortUnique(knots);
        std::sort(inputs.begin()+begin,inputs.begin()+end,[](const Input& a,const Input& b){return a.lo<b.lo;});
        std::vector<std::size_t> active;std::size_t next=begin;
        for(std::size_t k=1;k<knots.size();++k){const double lo=knots[k-1],hi=knots[k];
            active.erase(std::remove_if(active.begin(),active.end(),[&](std::size_t i){return inputs[i].hi<=lo;}),active.end());
            while(next<end && inputs[next].lo<=lo){active.push_back(next++);}if(active.empty())continue;
            const auto& origin=inputs[active.front()];for(auto i:active){const auto& v=inputs[i];require(v.hi>=hi && v.sourcePrimal==origin.sourcePrimal && v.sourceCell==origin.sourceCell,"ambiguous canonical atom origin");}
            HydroAdvectorAtom atom;atom.normalAxis=normal;atom.begin=point(normal,plane,lo);atom.end=point(normal,plane,hi);atom.primalVelocity=origin.sourcePrimal;atom.ownerPressure=origin.sourceCell;
            atom.area=(static_cast<long double>(hi)-lo)*(rz?CartesianMeasure::twoPi*(static_cast<long double>(atom.begin.x)+atom.end.x)/2:1.L);
            require(std::isfinite(atom.area) && atom.area>=0,"invalid atom area");const int id=static_cast<int>(d->atoms.size());d->atoms.push_back(atom);
            for(auto i:active){const auto& v=inputs[i];(v.primal>=0?d->primal[v.primal]:d->dual[v.dual]).push_back(id);++d->report.viewEntries;}
        }
        begin=end;
    }
    d->report.atoms=d->atoms.size();
    // Complete the original output allocation before publishing any cache slot.
    // A failed prepare therefore cannot leave a newly seeded bare preparation.
    auto result=HydroAdvectorReconstruction(std::make_shared<Data>(std::move(d)));
    if(cacheOwner){
        auto& retained=boundary?cacheOwner->data_->reconstructionPreparation:cacheOwner->data_->unboundReconstructionPreparation;
        retained=result.data_->preparation;
        // prepare() uses the boundary only to qualify the compatible/virtual
        // chart. For complete boundary lines its bare geometry is exactly the
        // null-boundary chart; every wall/base mode is already prepared. The
        // actual mask and velocity coefficients remain late-bound in apply().
        // Partial-line owners keep the two original, distinct cache slots.
        if(cacheOwner->completeBoundaryLines()){
            auto& other=boundary?cacheOwner->data_->unboundReconstructionPreparation:cacheOwner->data_->reconstructionPreparation;
            other=retained;
        }
    }
    return result;
}
HydroSpatialAdvector HydroAdvectorReconstruction::evaluate(const QuadtreeMesh& mesh,const HydroNativeValues& values,
    Interface::HermiteGalsIdentity geometry,Interface::HermiteGalsTimeLayer clock) const
{return evaluateBound(mesh,values,geometry,clock,nullptr,nullptr);}
HydroSpatialAdvector HydroAdvectorReconstruction::evaluate(const QuadtreeMesh& mesh,const HydroNativeValues& values,
    Interface::HermiteGalsIdentity geometry,Interface::HermiteGalsTimeLayer clock,const HydroAdvectorBoundary& boundary) const
{return evaluateBound(mesh,values,geometry,clock,&boundary,nullptr);}
HydroSpatialAdvector HydroAdvectorReconstruction::evaluate(const QuadtreeMesh& mesh,const HydroNativeValues& values,
    const std::vector<double>& tangential,Interface::HermiteGalsIdentity geometry,Interface::HermiteGalsTimeLayer clock,
    const HydroAdvectorBoundary* boundary) const
{return evaluateBound(mesh,values,geometry,clock,boundary,&tangential);}
HydroSpatialAdvector HydroAdvectorReconstruction::evaluateBound(const QuadtreeMesh& mesh,const HydroNativeValues& values,
    Interface::HermiteGalsIdentity geometry,Interface::HermiteGalsTimeLayer clock,const HydroAdvectorBoundary* boundary,
    const std::vector<double>* tangential) const
{
    fp();values.validateFor(momentumGeometry().layout(),mesh);clock.time();require(clock.step>=0 && geometry.value()!=0,"unbound provenance or negative step");
    require(tangential||data_->modes.empty(),"explicit tangential velocity state required");
    if(tangential){require(tangential->size()==data_->modes.size(),"tangential state size mismatch");
        for(double q:*tangential)require(std::isfinite(q),"nonfinite tangential velocity state");}
    if(data_->phase)data_->phase->validateFor(mesh,geometry,clock);
    if(boundary){boundary->validateFor(mesh,momentumGeometry());for(int u:boundary->velocityDofs())require(values.normalVelocity()[u]==0,"nonzero impermeable native mean");}
    const auto divergence=momentumGeometry().metrics().apply(mesh,values);
    auto d=std::make_shared<HydroSpatialAdvector::Data>(*this,values,geometry,clock,HydroAdvectorIdentity(nextIdentity.fetch_add(1)));
    if(boundary)d->boundary=*boundary;
    if(tangential)d->tangential=*tangential;
    const auto& regular=regularGeometry();
    d->coefficients=regular.apply(mesh,values,boundary);
    d->divergence.assign(divergence.begin(),divergence.end());
    for(const auto& atom:data_->preparation->atoms) {
        int p=atom.ownerPressure;
        if(atom.primalVelocity>=0){const auto& f=regular.faces()[atom.primalVelocity];p=f.negativePressure>=0?f.negativePressure:f.positivePressure;}
        d->polynomials.push_back(regular.normalPolynomial(d->coefficients,p,atom.normalAxis,atom.begin,atom.end));
    }
    return HydroSpatialAdvector(std::move(d));
}
HydroSpatialAdvector HydroAdvectorReconstruction::interpolate(const QuadtreeMesh& mesh,
    const HydroSpatialAdvector& previous,const HydroSpatialAdvector& arrival,
    Interface::HermiteGalsIdentity geometry,Interface::HermiteGalsTimeLayer clock) const
{
    fp();const auto& g=momentumGeometry();
    previous.validateFor(mesh,g,previous.geometryIdentity(),previous.timeLayer());
    arrival.validateFor(mesh,g,arrival.geometryIdentity(),arrival.timeLayer());
    require(previous.reconstruction().sameIdentity(*this)&&arrival.reconstruction().sameIdentity(*this),"different endpoint reconstruction");
    require((!previous.boundary()&&!arrival.boundary())||(previous.boundary()&&arrival.boundary()
        &&previous.boundary()->sameIdentity(*arrival.boundary())),"different endpoint boundary law");
    using P=Interface::HermiteGalsTimeLayer::Position;
    const auto a=previous.timeLayer(),b=arrival.timeLayer();const double time=clock.time();
    require((a.position==P::Initial||a.position==P::End)&&b.position==P::End&&b.step>a.step&&b.step-a.step==1
        &&b.begin==a.time()&&b.end>b.begin,"invalid endpoint slab");
    require(geometry.value()!=0&&clock.step==b.step&&clock.begin==b.begin&&clock.end<=b.end
        &&clock.position!=P::Initial&&time>a.time()&&time<=b.time(),"dense geometry clock outside endpoint slab");
    const double theta=(time-a.time())/(b.time()-a.time());
    const auto mix=[&](double x,double y){return time==b.time()?y:(1-theta)*x+theta*y;};
    auto normal=previous.values().normalVelocity();
    for(std::size_t i=0;i<normal.size();++i)normal[i]=mix(normal[i],arrival.values().normalVelocity()[i]);
    auto d=std::make_shared<HydroSpatialAdvector::Data>(*this,HydroNativeValues(g.layout(),normal,std::vector<double>(g.layout().pressureCount())),
        geometry,clock,HydroAdvectorIdentity(nextIdentity.fetch_add(1)));
    d->boundary=previous.data_->boundary;d->coefficients=previous.data_->coefficients;
    if(!d->coefficients.compatibleNativeMeans.empty())d->coefficients.compatibleNativeMeans=normal;
    const auto& other=arrival.data_->coefficients;
    // This is the SAME actual representation: general cell Hermite data or
    // compatible quarter-cell power patches. Interpolate every retained patch.
    require(d->coefficients.cellDofs.size()==other.cellDofs.size(),"different endpoint coefficient representation");
    for(std::size_t p=0;p<d->coefficients.cellDofs.size();++p)for(int c=0;c<2;++c)for(int k=0;k<25;++k)
        d->coefficients.cellDofs[p][c][k]=mix(d->coefficients.cellDofs[p][c][k],other.cellDofs[p][c][k]);
    for(std::size_t i=0;i<d->coefficients.faceCorrections.size();++i)
        d->coefficients.faceCorrections[i]=mix(d->coefficients.faceCorrections[i],other.faceCorrections[i]);
    d->polynomials=previous.atomPolynomials();
    for(std::size_t i=0;i<d->polynomials.size();++i)for(int k=0;k<5;++k){
        auto& p=d->polynomials[i];const auto& q=arrival.atomPolynomials()[i];
        p.beginCoefficients[k]=(1-theta)*p.beginCoefficients[k]+theta*q.beginCoefficients[k];
        p.endCoefficients[k]=(1-theta)*p.endCoefficients[k]+theta*q.endCoefficients[k];
    }
    d->divergence=previous.cellDivergence();
    for(std::size_t i=0;i<d->divergence.size();++i)d->divergence[i]=(1-theta)*d->divergence[i]+theta*arrival.cellDivergence()[i];
    d->tangential=previous.tangentialVelocity();
    for(std::size_t i=0;i<d->tangential.size();++i)d->tangential[i]=mix(d->tangential[i],arrival.tangentialVelocity()[i]);
    return HydroSpatialAdvector(std::move(d));
}
long double HydroNormalPolynomial::value(bool fromEnd,long double offset) const
{
    fp();require(std::isfinite(offset) && offset>=0 && offset<=length,"invalid polynomial chart offset");const auto& a=fromEnd?endCoefficients:beginCoefficients;
    long double v=a.back();for(std::size_t k=a.size()-1;k>0;--k)v=a[k-1]+offset*v;
    require(std::isfinite(v),"nonfinite polynomial value");return v;
}
long double HydroNormalPolynomial::integral(bool fromEnd,long double offset,long double width,int sign) const
{
    fp();require(std::isfinite(offset) && std::isfinite(width) && offset>=0 && width>0 && offset<=length && width<=length-offset && (sign==1 || sign==-1),"invalid shifted interval or normal sign");
    const auto& a=fromEnd?endCoefficients:beginCoefficients;
    // Translate the complete quartic before integrating the independent width.
    // No rounded global endpoint is constructed, even at a multiple zero.
    constexpr int choose[5][5]={{1,0,0,0,0},{1,1,0,0,0},{1,2,1,0,0},{1,3,3,1,0},{1,4,6,4,1}};
    std::array<long double,5> shifted{};
    for(std::size_t k=0;k<a.size();++k){long double v=choose[4][k]*a[4];
        for(int j=3;j>=static_cast<int>(k);--j)v=choose[j][k]*a[j]+offset*v;shifted[k]=v;}
    long double result=0,firstMoment=0;
    for(std::size_t k=0;k<shifted.size();++k){auto term=shifted[k];for(std::size_t j=0;j<=k;++j)term*=width;
        result+=term/(k+1);if(radialMeasure && radialTangent)firstMoment+=term*width/(k+2);}
    if(radialMeasure){const auto radius=(fromEnd?radiusEnd:radiusBegin)+(radialTangent?(fromEnd?-offset:offset):0);
        result=CartesianMeasure::twoPi*(radius*result+(radialTangent?(fromEnd?-1.L:1.L)*firstMoment:0));}
    require(std::isfinite(result),"nonfinite shifted integral");return sign*result;
}
HydroSpatialAdvector::HydroSpatialAdvector(std::shared_ptr<const Data> d):data_(std::move(d)){}
const HydroAdvectorReconstruction& HydroSpatialAdvector::reconstruction() const noexcept{return data_->reconstruction;}
const HydroMomentumGeometry& HydroSpatialAdvector::momentumGeometry() const noexcept{return reconstruction().momentumGeometry();}
const HydroNativeValues& HydroSpatialAdvector::values() const noexcept{return data_->values;}
const std::vector<double>& HydroSpatialAdvector::tangentialVelocity() const noexcept{return data_->tangential;}
Interface::HermiteGalsIdentity HydroSpatialAdvector::geometryIdentity() const noexcept{return data_->geometry;}
Interface::HermiteGalsTimeLayer HydroSpatialAdvector::timeLayer() const noexcept{return data_->clock;}
HydroAdvectorIdentity HydroSpatialAdvector::identity() const noexcept{return data_->identity;}
const HydroAdvectorBoundary* HydroSpatialAdvector::boundary() const noexcept{return data_->boundary?&*data_->boundary:nullptr;}
const std::vector<HydroNormalPolynomial>& HydroSpatialAdvector::atomPolynomials() const
{
    for(double q:data_->tangential)require(q==0,"physical tangential phase flux coefficients required");
    return data_->polynomials;
}
const std::vector<HydroNormalPolynomial>& HydroSpatialAdvector::regularAtomPolynomials() const noexcept{return data_->polynomials;}
long double HydroSpatialAdvector::tangentialAtomRate(int index,const LevelSetMappedSegment::NativeQuarticTrace& trace) const
{
    const auto& source=reconstruction();
    require(index>=0&&static_cast<std::size_t>(index)<source.atoms().size(),"tangential atom index out of range");
    const auto& atom=source.atoms()[index];
    require(trace.begin().x==atom.begin.x&&trace.begin().y==atom.begin.y&&trace.end().x==atom.end.x
        &&trace.end().y==atom.end.y&&trace.normalAxis()==atom.normalAxis,"different modal physical trace");
    if(atom.primalVelocity>=0)return 0; // Mode and first jet vanish on native cell boundaries.
    const int mode=source.tangentialModeForPressure(atom.ownerPressure);
    if(mode<0||data_->tangential[mode]==0)return 0;
    const auto* phase=source.phaseGeometry();require(phase!=nullptr,"modal atom has no phase source");
    const auto& scalar=phase->publication().state().geometry().snapshot();
    require(trace.usesSource(scalar),"different tangential atom Q4 source");
    const auto endpoint=[&](bool end)->long double {
        for(std::size_t event=0;event<trace.events().size();++event)
            if(trace.events()[event].normalOrder==0&&trace.rootAtInputEndpoint(event,end))return 0;
        for(const auto& zero:trace.coincidentRanges())
            if(zero.interval.offset==0&&(zero.interval.fromEnd==end||zero.interval.length==1))return 0;
        for(const auto& range:trace.ordinaryRanges())
            if(range.offset==0&&(range.fromEnd==end||range.length==1))
                return source.tangentialModeStreamfunction(mode,end?atom.end:atom.begin,range.negative);
        throw std::runtime_error("[HydroSpatialAdvector] unresolved modal atom endpoint phase");
    };
    const long double rate=data_->tangential[mode]*(atom.normalAxis==0?1.L:-1.L)*
        (momentumGeometry().layout().domain().axisymmetric()?CartesianMeasure::twoPi:1.L)*(endpoint(true)-endpoint(false));
    require(std::isfinite(rate),"nonfinite tangential atom rate");return rate;
}
const std::vector<long double>& HydroSpatialAdvector::cellDivergence() const noexcept{return data_->divergence;}
Vector2D HydroSpatialAdvector::cellVelocity(int p,Vector2D point) const
{return cellVelocityJet(p,point).jet.value;}
HydroAdvectorCellJet HydroSpatialAdvector::cellVelocityJet(int p,Vector2D point) const
{
    fp();const auto& cells=reconstruction().regularGeometry().cells();
    require(p>=0 && static_cast<std::size_t>(p)<cells.size(),"pressure index out of range");
    const auto& cell=cells[p];
    require(std::isfinite(point.x)&&std::isfinite(point.y)&&point.x>=cell.lower.x&&point.x<=cell.upper.x
        &&point.y>=cell.lower.y&&point.y<=cell.upper.y,"point outside actual cell");
    return reconstruction().regularGeometry().cellJet(data_->coefficients,p,point);
}
HydroAdvectorCellJet HydroSpatialAdvector::cellPhysicalVelocityJet(int p,Vector2D point,bool liquid) const
{
    auto result=cellVelocityJet(p,point);const int mode=reconstruction().tangentialModeForPressure(p);
    if(mode>=0){const auto jet=reconstruction().tangentialModeJet(mode,point,liquid);const double q=data_->tangential[mode];
        result.jet.value=result.jet.value+jet.value*q;
        result.jet.derivative.xx+=q*jet.derivative.xx;result.jet.derivative.xy+=q*jet.derivative.xy;
        result.jet.derivative.yx+=q*jet.derivative.yx;result.jet.derivative.yy+=q*jet.derivative.yy;}
    return result;
}
HydroAdvectorCellJet HydroSpatialAdvector::pieceVelocityJet(int index,Vector2D point) const
{
    fp();const auto& pieces=reconstruction().pieces();
    require(index>=0 && static_cast<std::size_t>(index)<pieces.size(),"piece index out of range");
    const auto& piece=pieces[index];
    require(point.x>=piece.lower.x && point.x<=piece.upper.x && point.y>=piece.lower.y && point.y<=piece.upper.y,"point outside selected piece");
    auto result=cellVelocityJet(piece.pressure,point);
    // Value and first derivatives are common at every support boundary.
    result.lower=piece.lower;result.upper=piece.upper;return result;
}
HydroAdvectorCellJetPair HydroSpatialAdvector::pieceVelocityJetPair(
    const HydroSpatialAdvector& arrival,int index,Vector2D point) const
{
    fp();require(reconstruction().sameIdentity(arrival.reconstruction()),
        "different reconstruction for paired piece jet");
    const auto& pieces=reconstruction().pieces();
    require(index>=0&&static_cast<std::size_t>(index)<pieces.size(),"piece index out of range");
    const auto& piece=pieces[index];
    require(point.x>=piece.lower.x&&point.x<=piece.upper.x&&point.y>=piece.lower.y&&point.y<=piece.upper.y,
        "point outside selected piece");
    const auto result=reconstruction().regularGeometry().pieceJetPair(data_->coefficients,arrival.data_->coefficients,piece,point);
    fp();return result;
}
HydroAdvectorCurveIntegral HydroSpatialAdvector::pieceCurveIntegral(int piece,Vector2D r0,Vector2D r1,Vector2D r2,
    LevelSetBoundaryTrace::Jacobian P0,LevelSetBoundaryTrace::Jacobian P1,
    LevelSetBoundaryTrace::Jacobian P2,double lo,double hi) const
{
    fp();const auto& pieces=reconstruction().pieces();
    require(piece>=0&&static_cast<std::size_t>(piece)<pieces.size(),"piece index out of range");
    const auto& regular=reconstruction().regularGeometry();
    const auto result=regular.curveIntegral(data_->curvePolynomial(piece),r0,r1,r2,P0,P1,P2,lo,hi);
    fp();return result;
}
HydroAdvectorCurveIntegralPair HydroSpatialAdvector::pieceCurveIntegralPair(
    const HydroSpatialAdvector& arrival,int piece,Vector2D r0,Vector2D r1,Vector2D r2,
    LevelSetBoundaryTrace::Jacobian P0,LevelSetBoundaryTrace::Jacobian P1,
    LevelSetBoundaryTrace::Jacobian P2,double lo,double hi) const
{
    fp();require(reconstruction().sameIdentity(arrival.reconstruction()),
        "different reconstruction for paired curve");
    const auto& pieces=reconstruction().pieces();
    require(piece>=0&&static_cast<std::size_t>(piece)<pieces.size(),"piece index out of range");
    const auto& regular=reconstruction().regularGeometry();
    const auto result=regular.curveIntegralPair(data_->curvePolynomial(piece),
        arrival.data_->curvePolynomial(piece),r0,r1,r2,P0,P1,P2,lo,hi);
    fp();return result;
}
Vector2D HydroSpatialAdvector::nativeVelocity(int u) const
{
    fp();const auto& faces=reconstruction().regularGeometry().faces();
    require(u>=0 && static_cast<std::size_t>(u)<faces.size(),"velocity index out of range");
    const auto& f=faces[u];const int p=f.negativePressure>=0?f.negativePressure:f.positivePressure;
    return cellVelocity(p,(f.lower+f.upper)*.5);
}
struct HydroAdvectorTimeSlab::Data
{
    Data(Interface::HermiteGalsPublication source,HydroAdvectorReconstruction r,HydroSpatialAdvector p,HydroSpatialAdvector a)
        : original(std::move(source)),reconstruction(std::move(r)),previous(std::move(p)),arrival(std::move(a)) {}
    Interface::HermiteGalsPublication original;
    HydroAdvectorReconstruction reconstruction;
    HydroSpatialAdvector previous,arrival;
    int nx=0,ny=0;
    std::vector<int> gridCells;
    std::vector<std::array<Vector2D,2>> bounds;
    std::vector<double> lowerX,upperX,lowerY,upperY;
    std::vector<std::size_t> atomCellOffsets;
    std::vector<int> atomCells;
};
HydroAdvectorTimeSlab::HydroAdvectorTimeSlab(std::shared_ptr<const Data> data):data_(std::move(data)) {}
HydroAdvectorTimeSlab HydroAdvectorTimeSlab::prepare(const QuadtreeMesh& mesh,Interface::HermiteGalsPublication original,
    HydroAdvectorReconstruction reconstruction,const HydroSpatialAdvector& previous,const HydroSpatialAdvector& arrival)
{
    fp();original.validateMesh(mesh);const auto& g=reconstruction.momentumGeometry();
    previous.validateFor(mesh,g,previous.geometryIdentity(),previous.timeLayer());
    arrival.validateFor(mesh,g,arrival.geometryIdentity(),arrival.timeLayer());
    require(previous.reconstruction().sameIdentity(reconstruction)&&arrival.reconstruction().sameIdentity(reconstruction),"different endpoint reconstruction");
    require((!previous.boundary()&&!arrival.boundary())||(previous.boundary()&&arrival.boundary()
        &&previous.boundary()->sameIdentity(*arrival.boundary())),"different endpoint boundary law");
    using P=Interface::HermiteGalsTimeLayer::Position;const auto a=previous.timeLayer(),b=arrival.timeLayer();
    require((a.position==P::Initial||a.position==P::End)&&b.position==P::End&&b.step>a.step&&b.step-a.step==1
        &&b.begin==a.time()&&b.end>b.begin,"invalid endpoint slab");
    require(previous.geometryIdentity()==original.state().identity()&&sameTime(a,original.state().timeLayer()),
        "different original publication");
    require(!g.layout().domain().axisymmetric()&&original.state().primary().measure==Interface::HermiteGalsMeasure::CartesianArea,
        "local time slab requires uniform XY measure");
    // Preserve interpolate's original single-polynomial port rejection. All
    // endpoint modal coefficients must be zero, so their temporal mix is zero.
    (void)previous.atomPolynomials();(void)arrival.atomPolynomials();
    const UniformLeafGrid grid(mesh,CartesianGridScope::UniformLeaves);
    auto d=std::make_shared<Data>(std::move(original),std::move(reconstruction),previous,arrival);
    d->nx=grid.nx();d->ny=grid.ny();d->gridCells.resize(grid.cells().size());d->bounds.resize(mesh.cellIndexNum());
    const double infinity=std::numeric_limits<double>::infinity();
    d->lowerX.assign(d->nx,infinity);d->upperX.assign(d->nx,-infinity);
    d->lowerY.assign(d->ny,infinity);d->upperY.assign(d->ny,-infinity);
    for(int j=0;j<d->ny;++j)for(int i=0;i<d->nx;++i){const auto& cell=mesh.cells(grid.cells()[j*d->nx+i]);const int c=cell.index();
        const auto a=mesh.verts(cell.cellToVert[0]).coord(),b=mesh.verts(cell.cellToVert[2]).coord();
        const Vector2D lo{std::min(a.x,b.x),std::min(a.y,b.y)},hi{std::max(a.x,b.x),std::max(a.y,b.y)};
        d->gridCells[j*d->nx+i]=c;d->bounds[c]={{lo,hi}};
        d->lowerX[i]=std::min(d->lowerX[i],lo.x);d->upperX[i]=std::max(d->upperX[i],hi.x);
        d->lowerY[j]=std::min(d->lowerY[j],lo.y);d->upperY[j]=std::max(d->upperY[j],hi.y);}
    // Exact original closed-box incidence, prepared once rather than per time node.
    d->atomCellOffsets.push_back(0);
    for(const auto& atom:d->reconstruction.atoms()) {
        require((atom.normalAxis==0||atom.normalAxis==1)&&std::isfinite(atom.begin.x)&&std::isfinite(atom.begin.y)
            &&std::isfinite(atom.end.x)&&std::isfinite(atom.end.y),"invalid actual atom geometry");
        std::vector<int> cells;
        const Vector2D lo{std::min(atom.begin.x,atom.end.x),std::min(atom.begin.y,atom.end.y)},
            hi{std::max(atom.begin.x,atom.end.x),std::max(atom.begin.y,atom.end.y)};
        const int firstX=static_cast<int>(std::lower_bound(d->upperX.begin(),d->upperX.end(),lo.x)-d->upperX.begin()),
            endX=static_cast<int>(std::upper_bound(d->lowerX.begin(),d->lowerX.end(),hi.x)-d->lowerX.begin()),
            firstY=static_cast<int>(std::lower_bound(d->upperY.begin(),d->upperY.end(),lo.y)-d->upperY.begin()),
            endY=static_cast<int>(std::upper_bound(d->lowerY.begin(),d->lowerY.end(),hi.y)-d->lowerY.begin());
        bool found=false;for(int j=firstY;j<endY;++j)for(int i=firstX;i<endX;++i){const int c=d->gridCells[j*d->nx+i];
            const auto& bounds=d->bounds[c];if(hi.x>=bounds[0].x&&lo.x<=bounds[1].x&&hi.y>=bounds[0].y&&lo.y<=bounds[1].y){
                cells.push_back(c);found=true;}}
        require(found,"actual atom has no native support");
        std::sort(cells.begin(),cells.end());
        d->atomCells.insert(d->atomCells.end(),cells.begin(),cells.end());
        d->atomCellOffsets.push_back(d->atomCells.size());
    }
    fp();return HydroAdvectorTimeSlab(std::move(d));
}
const HydroAdvectorReconstruction& HydroAdvectorTimeSlab::reconstruction() const noexcept{return data_->reconstruction;}
Interface::HermiteGalsIdentity HydroAdvectorTimeSlab::sourceIdentity() const noexcept{return data_->original.state().identity();}
Interface::HermiteGalsIdentity HydroAdvectorTimeSlab::arrivalGeometryIdentity() const noexcept{return data_->arrival.geometryIdentity();}
HydroAdvectorIdentity HydroAdvectorTimeSlab::previousIdentity() const noexcept{return data_->previous.identity();}
HydroAdvectorIdentity HydroAdvectorTimeSlab::arrivalIdentity() const noexcept{return data_->arrival.identity();}
Interface::HermiteGalsTimeLayer HydroAdvectorTimeSlab::endpointTimeLayer() const noexcept{return data_->arrival.timeLayer();}
void HydroAdvectorTimeSlab::validateSample(Interface::HermiteGalsIdentity geometry,Interface::HermiteGalsIdentity source,
    Interface::HermiteGalsIdentity arrival,Interface::HermiteGalsTimeLayer clock) const
{
    fp();require(source==sourceIdentity(),"different original publication");
    require(arrival==arrivalGeometryIdentity(),"different arrival velocity geometry identity");
    const auto a=data_->previous.timeLayer(),b=data_->arrival.timeLayer();const double time=clock.time();
    require(geometry.value()!=0&&clock.step==b.step&&clock.begin==b.begin&&clock.end<=b.end
        &&clock.position!=Interface::HermiteGalsTimeLayer::Position::Initial&&time>a.time()&&time<=b.time(),
        "dense geometry clock outside endpoint slab");
}
HydroNormalPolynomial HydroAdvectorTimeSlab::atomPolynomial(int atom,Interface::HermiteGalsTimeLayer clock) const
{
    fp();validateSample(sourceIdentity(),sourceIdentity(),arrivalGeometryIdentity(),clock);
    require(atom>=0&&static_cast<std::size_t>(atom)<reconstruction().atoms().size(),"invalid atom");
    const auto a=data_->previous.timeLayer(),b=data_->arrival.timeLayer();const double theta=(clock.time()-a.time())/(b.time()-a.time());
    auto p=data_->previous.regularAtomPolynomials()[atom];const auto& q=data_->arrival.regularAtomPolynomials()[atom];
    for(int k=0;k<5;++k){p.beginCoefficients[k]=(1-theta)*p.beginCoefficients[k]+theta*q.beginCoefficients[k];
        p.endCoefficients[k]=(1-theta)*p.endCoefficients[k]+theta*q.endCoefficients[k];}
    return p;
}
std::size_t HydroAdvectorTimeSlab::nativeCellCount() const noexcept{return data_->bounds.size();}
std::array<Vector2D,2> HydroAdvectorTimeSlab::nativeCellBounds(int cell) const
{fp();require(cell>=0&&static_cast<std::size_t>(cell)<data_->bounds.size(),"native cell index out of range");return data_->bounds[cell];}
HydroAdvectorAtomSupport HydroAdvectorTimeSlab::prepareAtomSupport(const std::vector<int>& atomIds) const
{
    fp();require(!atomIds.empty(),"empty requested atoms");
    std::vector<unsigned char> seenAtoms(reconstruction().atoms().size(),0),seenCells(data_->bounds.size(),0);
    std::vector<int> cells;cells.reserve(std::min(data_->bounds.size(),atomIds.size()*4));
    for(int id:atomIds) {
        require(id>=0&&static_cast<std::size_t>(id)<seenAtoms.size(),"invalid atom");
        require(!seenAtoms[id],"duplicate atom");seenAtoms[id]=1;
        for(std::size_t k=data_->atomCellOffsets[id];k<data_->atomCellOffsets[id+1];++k) {
            const int cell=data_->atomCells[k];if(!seenCells[cell]){seenCells[cell]=1;cells.push_back(cell);}
        }
    }
    std::sort(cells.begin(),cells.end());fp();
    return HydroAdvectorAtomSupport(data_,atomIds,std::move(cells));
}
void HydroAdvectorTimeSlab::validateAtomSupport(const HydroAdvectorAtomSupport& support) const
{fp();require(support.owner_.get()==data_.get(),"different prepared atom support owner");}
std::vector<int> HydroAdvectorTimeSlab::nativeCells(const std::vector<int>& atomIds) const
{return prepareAtomSupport(atomIds).nativeCells();}

void HydroSpatialAdvector::validateFor(const QuadtreeMesh& mesh,const HydroMomentumGeometry& expected,
    Interface::HermiteGalsIdentity geometry,Interface::HermiteGalsTimeLayer clock) const
{
    fp();momentumGeometry().validateMesh(mesh);require(momentumGeometry().sameIdentity(expected),"different momentum geometry instance");
    require(geometry==geometryIdentity(),"different velocity geometry identity");require(sameTime(clock,timeLayer()),"different full time layer");
}
}
