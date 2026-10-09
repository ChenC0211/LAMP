#include "Modules/Hydrodynamics/CchSpatialAdvector.h"
#include "Modules/Hydrodynamics/Detail/CchRegularAdvectorGeometry.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/CartesianMeasure.h"
#include <algorithm>
#include <atomic>
#include <cfenv>
#include <cmath>
#include <map>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#if defined(_M_X64) || defined(__SSE2__)
#include <xmmintrin.h>
#endif

namespace Quantum::Hydrodynamics
{
namespace
{
void require(bool ok, const char* why)
{ if (!ok) throw std::invalid_argument(std::string("[CchSpatialAdvector] ")+why); }
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
struct CchAdvectorReconstruction::Data
{
    explicit Data(CchMomentumGeometry g):geometry(std::move(g)){}
    CchMomentumGeometry geometry;
    std::vector<CchAdvectorAtom> atoms;
    std::vector<CchAdvectorPiece> pieces;
    std::vector<std::vector<int>> primal,dual;
    CchAdvectorPreparationReport report;
    std::shared_ptr<const Detail::CchRegularAdvectorGeometry> regular;
    std::optional<CchPhaseGeometry> phase;
    std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> scalar;
    std::vector<CchTangentialMode> modes;
    std::vector<int> modeByPressure;
};
struct CchSpatialAdvector::Data
{
    Data(CchAdvectorReconstruction r,CchNativeValues v,Interface::HermiteGalsIdentity g,
        Interface::HermiteGalsTimeLayer t,CchAdvectorIdentity i)
        : reconstruction(std::move(r)),values(std::move(v)),geometry(g),clock(t),identity(i){}
    CchAdvectorReconstruction reconstruction;
    CchNativeValues values;
    Interface::HermiteGalsIdentity geometry;
    Interface::HermiteGalsTimeLayer clock;
    CchAdvectorIdentity identity;
    std::optional<CchAdvectorBoundary> boundary;
    Detail::CchRegularCoefficients coefficients;
    std::vector<long double> divergence;
    std::vector<CchNormalPolynomial> polynomials;
    std::vector<double> tangential;
};
struct CchAdvectorBoundary::Data
{
    Data(CchMomentumGeometry g,std::vector<int> v):geometry(std::move(g)),velocities(std::move(v)){}
    CchMomentumGeometry geometry;
    std::vector<int> velocities;
};
CchAdvectorBoundary::CchAdvectorBoundary(std::shared_ptr<const Data> d):data_(std::move(d)){}
CchAdvectorBoundary CchAdvectorBoundary::impermeable(const QuadtreeMesh& mesh,CchMomentumGeometry geometry,const std::vector<int>& dofs)
{
    fp();geometry.validateMesh(mesh);auto ordered=dofs;std::sort(ordered.begin(),ordered.end());
    require(std::adjacent_find(ordered.begin(),ordered.end())==ordered.end(),"duplicate impermeable face");
    const auto& l=geometry.layout();for(int u:ordered){require(u>=0 && u<l.velocityCount(),"impermeable velocity index out of range");
        require(l.domain().face(mesh,l.nativeFace(u)).outwardSign!=0,"impermeable law requires selected-domain boundary face");}
    return CchAdvectorBoundary(std::make_shared<Data>(std::move(geometry),std::move(ordered)));
}
const CchMomentumGeometry& CchAdvectorBoundary::momentumGeometry() const noexcept{return data_->geometry;}
const std::vector<int>& CchAdvectorBoundary::velocityDofs() const noexcept{return data_->velocities;}
bool CchAdvectorBoundary::sameIdentity(const CchAdvectorBoundary& other) const noexcept{return data_==other.data_;}
void CchAdvectorBoundary::validateFor(const QuadtreeMesh& mesh,const CchMomentumGeometry& expected) const
{fp();momentumGeometry().validateMesh(mesh);require(momentumGeometry().sameIdentity(expected),"different impermeable momentum geometry instance");}
CchAdvectorReconstruction::CchAdvectorReconstruction(std::shared_ptr<const Data> d):data_(std::move(d)){}
const CchMomentumGeometry& CchAdvectorReconstruction::momentumGeometry() const noexcept{return data_->geometry;}
const std::vector<CchAdvectorAtom>& CchAdvectorReconstruction::atoms() const noexcept{return data_->atoms;}
const std::vector<CchAdvectorPiece>& CchAdvectorReconstruction::pieces() const noexcept{return data_->pieces;}
const std::vector<std::vector<int>>& CchAdvectorReconstruction::primalAtoms() const noexcept{return data_->primal;}
const std::vector<std::vector<int>>& CchAdvectorReconstruction::dualAtoms() const noexcept{return data_->dual;}
const CchAdvectorPreparationReport& CchAdvectorReconstruction::report() const noexcept{return data_->report;}
const Detail::CchRegularAdvectorGeometry& CchAdvectorReconstruction::regularGeometry() const noexcept{return *data_->regular;}
bool CchAdvectorReconstruction::sameIdentity(const CchAdvectorReconstruction& r) const noexcept{return data_==r.data_;}
const CchPhaseGeometry* CchAdvectorReconstruction::phaseGeometry() const noexcept
{return data_->phase?&*data_->phase:nullptr;}
const std::vector<CchTangentialMode>& CchAdvectorReconstruction::tangentialModes() const noexcept{return data_->modes;}
int CchAdvectorReconstruction::tangentialModeForPressure(int p) const
{
    require(p>=0&&p<momentumGeometry().layout().pressureCount(),"modal pressure index out of range");
    return data_->modeByPressure.empty()?-1:data_->modeByPressure[p];
}
CchAdvectorReconstruction CchAdvectorReconstruction::withTangentialModes(const QuadtreeMesh& mesh,CchPhaseGeometry phase) const
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
    return CchAdvectorReconstruction(std::move(d));
}
LevelSetBoundaryTrace::VelocityJet CchAdvectorReconstruction::tangentialModeJet(int index,Vector2D point,bool liquid) const
{
    fp();require(index>=0&&static_cast<std::size_t>(index)<data_->modes.size(),"tangential mode index out of range");
    const auto& mode=data_->modes[index];
    if(!(std::isfinite(point.x)&&std::isfinite(point.y)&&point.x>=mode.lower.x&&point.x<=mode.upper.x
        &&point.y>=mode.lower.y&&point.y<=mode.upper.y)){
        std::ostringstream message;message.precision(17);
        message<<"[CchSpatialAdvector] point outside tangential mode cell: mode="<<index<<", pressure="<<mode.pressure
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
long double CchAdvectorReconstruction::tangentialModeStreamfunction(int index,Vector2D point,bool liquid) const
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

CchAdvectorReconstruction CchAdvectorReconstruction::prepare(const QuadtreeMesh& mesh,CchMomentumGeometry geometry)
{
    fp();geometry.validateMesh(mesh);auto d=std::make_shared<Data>(std::move(geometry));
    d->regular=std::make_shared<const Detail::CchRegularAdvectorGeometry>(Detail::CchRegularAdvectorGeometry::prepare(mesh,d->geometry));
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
    struct Input{int primal=-1,dual=-1,sourcePrimal=-1,sourceCell=-1;double lo=0,hi=0;};
    struct Line{std::vector<Input> inputs;std::vector<double> knots;};
    std::map<std::pair<int,double>,Line> lines;
    const auto append=[&](int normal,Vector2D a,Vector2D b,int primal,int dual,int sourcePrimal,int sourceCell) {
        auto& line=lines[{normal,static_cast<double>(c(a,normal))}];
        const double lo=static_cast<double>(c(a,1-normal)),hi=static_cast<double>(c(b,1-normal));
        require(lo<hi,"nonpositive actual segment");line.inputs.push_back({primal,dual,sourcePrimal,sourceCell,lo,hi});line.knots.push_back(lo);line.knots.push_back(hi);
        const int pressure=sourcePrimal>=0?owner(sourcePrimal):sourceCell;
        for(double knot:cells[pressure].supportKnots[1-normal])if(lo<knot && knot<hi)line.knots.push_back(knot);
    };
    d->primal.resize(l.velocityCount());d->dual.resize(d->geometry.segments().size());
    for(int u=0;u<l.velocityCount();++u){const auto& f=faces[u];append(f.normalAxis,f.lower,f.upper,u,-1,u,-1);}
    for(std::size_t i=0;i<d->geometry.segments().size();++i){const auto& s=d->geometry.segments()[i];append(s.normalAxis,s.begin,s.end,-1,static_cast<int>(i),s.primalVelocity,s.ownerPressure);}
    for(auto& entry:lines){auto& line=entry.second;sortUnique(line.knots);
        std::sort(line.inputs.begin(),line.inputs.end(),[](const Input& a,const Input& b){return a.lo<b.lo;});std::vector<std::size_t> active;std::size_t next=0;
        for(std::size_t k=1;k<line.knots.size();++k){const double lo=line.knots[k-1],hi=line.knots[k];
            active.erase(std::remove_if(active.begin(),active.end(),[&](std::size_t i){return line.inputs[i].hi<=lo;}),active.end());
            while(next<line.inputs.size() && line.inputs[next].lo<=lo){active.push_back(next++);}if(active.empty())continue;
            const auto& origin=line.inputs[active.front()];for(auto i:active){const auto& v=line.inputs[i];require(v.hi>=hi && v.sourcePrimal==origin.sourcePrimal && v.sourceCell==origin.sourceCell,"ambiguous canonical atom origin");}
            const auto [normal,plane]=entry.first;CchAdvectorAtom atom;atom.normalAxis=normal;atom.begin=point(normal,plane,lo);atom.end=point(normal,plane,hi);atom.primalVelocity=origin.sourcePrimal;atom.ownerPressure=origin.sourceCell;
            atom.area=(static_cast<long double>(hi)-lo)*(rz?CartesianMeasure::twoPi*(static_cast<long double>(atom.begin.x)+atom.end.x)/2:1.L);
            require(std::isfinite(atom.area) && atom.area>=0,"invalid atom area");const int id=static_cast<int>(d->atoms.size());d->atoms.push_back(atom);
            for(auto i:active){const auto& v=line.inputs[i];(v.primal>=0?d->primal[v.primal]:d->dual[v.dual]).push_back(id);++d->report.viewEntries;}
        }
    }
    d->report.atoms=d->atoms.size();return CchAdvectorReconstruction(std::move(d));
}
CchSpatialAdvector CchAdvectorReconstruction::evaluate(const QuadtreeMesh& mesh,const CchNativeValues& values,
    Interface::HermiteGalsIdentity geometry,Interface::HermiteGalsTimeLayer clock) const
{return evaluateBound(mesh,values,geometry,clock,nullptr,nullptr);}
CchSpatialAdvector CchAdvectorReconstruction::evaluate(const QuadtreeMesh& mesh,const CchNativeValues& values,
    Interface::HermiteGalsIdentity geometry,Interface::HermiteGalsTimeLayer clock,const CchAdvectorBoundary& boundary) const
{return evaluateBound(mesh,values,geometry,clock,&boundary,nullptr);}
CchSpatialAdvector CchAdvectorReconstruction::evaluate(const QuadtreeMesh& mesh,const CchNativeValues& values,
    const std::vector<double>& tangential,Interface::HermiteGalsIdentity geometry,Interface::HermiteGalsTimeLayer clock,
    const CchAdvectorBoundary* boundary) const
{return evaluateBound(mesh,values,geometry,clock,boundary,&tangential);}
CchSpatialAdvector CchAdvectorReconstruction::evaluateBound(const QuadtreeMesh& mesh,const CchNativeValues& values,
    Interface::HermiteGalsIdentity geometry,Interface::HermiteGalsTimeLayer clock,const CchAdvectorBoundary* boundary,
    const std::vector<double>* tangential) const
{
    fp();values.validateFor(momentumGeometry().layout(),mesh);clock.time();require(clock.step>=0 && geometry.value()!=0,"unbound provenance or negative step");
    require(tangential||data_->modes.empty(),"explicit tangential velocity state required");
    if(tangential){require(tangential->size()==data_->modes.size(),"tangential state size mismatch");
        for(double q:*tangential)require(std::isfinite(q),"nonfinite tangential velocity state");}
    if(data_->phase)data_->phase->validateFor(mesh,geometry,clock);
    if(boundary){boundary->validateFor(mesh,momentumGeometry());for(int u:boundary->velocityDofs())require(values.normalVelocity()[u]==0,"nonzero impermeable native mean");}
    const auto divergence=momentumGeometry().metrics().apply(mesh,values);
    auto d=std::make_shared<CchSpatialAdvector::Data>(*this,values,geometry,clock,CchAdvectorIdentity(nextIdentity.fetch_add(1)));
    if(boundary)d->boundary=*boundary;
    if(tangential)d->tangential=*tangential;
    const auto& regular=regularGeometry();
    d->coefficients=regular.apply(mesh,values,boundary);
    d->divergence.assign(divergence.begin(),divergence.end());
    for(const auto& atom:data_->atoms) {
        int p=atom.ownerPressure;
        if(atom.primalVelocity>=0){const auto& f=regular.faces()[atom.primalVelocity];p=f.negativePressure>=0?f.negativePressure:f.positivePressure;}
        d->polynomials.push_back(regular.normalPolynomial(d->coefficients,p,atom.normalAxis,atom.begin,atom.end));
    }
    return CchSpatialAdvector(std::move(d));
}
CchSpatialAdvector CchAdvectorReconstruction::interpolate(const QuadtreeMesh& mesh,
    const CchSpatialAdvector& previous,const CchSpatialAdvector& arrival,
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
    auto d=std::make_shared<CchSpatialAdvector::Data>(*this,CchNativeValues(g.layout(),normal,std::vector<double>(g.layout().pressureCount())),
        geometry,clock,CchAdvectorIdentity(nextIdentity.fetch_add(1)));
    d->boundary=previous.data_->boundary;d->coefficients=previous.data_->coefficients;
    const auto& other=arrival.data_->coefficients;
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
    return CchSpatialAdvector(std::move(d));
}
long double CchNormalPolynomial::value(bool fromEnd,long double offset) const
{
    fp();require(std::isfinite(offset) && offset>=0 && offset<=length,"invalid polynomial chart offset");const auto& a=fromEnd?endCoefficients:beginCoefficients;
    long double v=a.back();for(std::size_t k=a.size()-1;k>0;--k)v=a[k-1]+offset*v;
    require(std::isfinite(v),"nonfinite polynomial value");return v;
}
long double CchNormalPolynomial::integral(bool fromEnd,long double offset,long double width,int sign) const
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
CchSpatialAdvector::CchSpatialAdvector(std::shared_ptr<const Data> d):data_(std::move(d)){}
const CchAdvectorReconstruction& CchSpatialAdvector::reconstruction() const noexcept{return data_->reconstruction;}
const CchMomentumGeometry& CchSpatialAdvector::momentumGeometry() const noexcept{return reconstruction().momentumGeometry();}
const CchNativeValues& CchSpatialAdvector::values() const noexcept{return data_->values;}
const std::vector<double>& CchSpatialAdvector::tangentialVelocity() const noexcept{return data_->tangential;}
Interface::HermiteGalsIdentity CchSpatialAdvector::geometryIdentity() const noexcept{return data_->geometry;}
Interface::HermiteGalsTimeLayer CchSpatialAdvector::timeLayer() const noexcept{return data_->clock;}
CchAdvectorIdentity CchSpatialAdvector::identity() const noexcept{return data_->identity;}
const CchAdvectorBoundary* CchSpatialAdvector::boundary() const noexcept{return data_->boundary?&*data_->boundary:nullptr;}
const std::vector<CchNormalPolynomial>& CchSpatialAdvector::atomPolynomials() const
{
    for(double q:data_->tangential)require(q==0,"physical tangential phase flux coefficients required");
    return data_->polynomials;
}
const std::vector<CchNormalPolynomial>& CchSpatialAdvector::regularAtomPolynomials() const noexcept{return data_->polynomials;}
long double CchSpatialAdvector::tangentialAtomRate(int index,const LevelSetMappedSegment::NativeQuarticTrace& trace) const
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
        throw std::runtime_error("[CchSpatialAdvector] unresolved modal atom endpoint phase");
    };
    const long double rate=data_->tangential[mode]*(atom.normalAxis==0?1.L:-1.L)*
        (momentumGeometry().layout().domain().axisymmetric()?CartesianMeasure::twoPi:1.L)*(endpoint(true)-endpoint(false));
    require(std::isfinite(rate),"nonfinite tangential atom rate");return rate;
}
const std::vector<long double>& CchSpatialAdvector::cellDivergence() const noexcept{return data_->divergence;}
Vector2D CchSpatialAdvector::cellVelocity(int p,Vector2D point) const
{return cellVelocityJet(p,point).jet.value;}
CchAdvectorCellJet CchSpatialAdvector::cellVelocityJet(int p,Vector2D point) const
{
    fp();const auto& cells=reconstruction().regularGeometry().cells();
    require(p>=0 && static_cast<std::size_t>(p)<cells.size(),"pressure index out of range");
    const auto& cell=cells[p];
    require(std::isfinite(point.x)&&std::isfinite(point.y)&&point.x>=cell.lower.x&&point.x<=cell.upper.x
        &&point.y>=cell.lower.y&&point.y<=cell.upper.y,"point outside actual cell");
    return reconstruction().regularGeometry().cellJet(data_->coefficients,p,point);
}
CchAdvectorCellJet CchSpatialAdvector::cellPhysicalVelocityJet(int p,Vector2D point,bool liquid) const
{
    auto result=cellVelocityJet(p,point);const int mode=reconstruction().tangentialModeForPressure(p);
    if(mode>=0){const auto jet=reconstruction().tangentialModeJet(mode,point,liquid);const double q=data_->tangential[mode];
        result.jet.value=result.jet.value+jet.value*q;
        result.jet.derivative.xx+=q*jet.derivative.xx;result.jet.derivative.xy+=q*jet.derivative.xy;
        result.jet.derivative.yx+=q*jet.derivative.yx;result.jet.derivative.yy+=q*jet.derivative.yy;}
    return result;
}
CchAdvectorCellJet CchSpatialAdvector::pieceVelocityJet(int index,Vector2D point) const
{
    fp();const auto& pieces=reconstruction().pieces();
    require(index>=0 && static_cast<std::size_t>(index)<pieces.size(),"piece index out of range");
    const auto& piece=pieces[index];
    require(point.x>=piece.lower.x && point.x<=piece.upper.x && point.y>=piece.lower.y && point.y<=piece.upper.y,"point outside selected piece");
    auto result=cellVelocityJet(piece.pressure,point);
    // Value and first derivatives are common at every support boundary.
    result.lower=piece.lower;result.upper=piece.upper;return result;
}
Vector2D CchSpatialAdvector::nativeVelocity(int u) const
{
    fp();const auto& faces=reconstruction().regularGeometry().faces();
    require(u>=0 && static_cast<std::size_t>(u)<faces.size(),"velocity index out of range");
    const auto& f=faces[u];const int p=f.negativePressure>=0?f.negativePressure:f.positivePressure;
    return cellVelocity(p,(f.lower+f.upper)*.5);
}
void CchSpatialAdvector::validateFor(const QuadtreeMesh& mesh,const CchMomentumGeometry& expected,
    Interface::HermiteGalsIdentity geometry,Interface::HermiteGalsTimeLayer clock) const
{
    fp();momentumGeometry().validateMesh(mesh);require(momentumGeometry().sameIdentity(expected),"different momentum geometry instance");
    require(geometry==geometryIdentity(),"different velocity geometry identity");require(sameTime(clock,timeLayer()),"different full time layer");
}
}
