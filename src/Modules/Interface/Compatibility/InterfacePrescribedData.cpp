#include "Modules/Interface/Compatibility/Detail/InterfaceStateInternal.h"
#include "Modules/Interface/Compatibility/InterfacePrescribedData.h"
#include "Modules/Interface/Compatibility/InterfaceMotion.h"
#include "Core/Geometry/LevelSetMappedSegment.h"
#include "Core/Mesh/NativeLeafDomain.h"
#include "Core/Mesh/QuadtreeMesh.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

namespace InterfacePrescribedData
{
namespace {
long double vortexSine(long double p)
{
    const auto half=std::round(2*p),local=p-half/2;
    return (static_cast<int>(half)%2?-1:1)*std::sin(2*std::acos(-1.L)*local);
}
// Keep the local distance to the analytic zero before adding a global
// coordinate. On Windows long double has double precision; forming the
// midpoint near 1/2 first can erase most of a very short phase subinterval.
long double vortexSineOffset(long double anchor,long double length,long double offset)
{
    const auto half=std::round(2*anchor);
    const auto local=std::fma(length,offset,anchor-half/2);
    return (static_cast<int>(half)%2?-1:1)*std::sin(2*std::acos(-1.L)*local);
}
std::shared_ptr<const InterfaceNormalFluxView> slipReversalFlux(const QuadtreeMesh& mesh,const NativeLeafDomain& domain){
    InterfaceNormalFluxView::Queries q;
    q.chart=[&mesh,&domain](int row){
        const auto f=domain.face(mesh,row);const auto& edge=f.horizontal?mesh.hEdges(f.edgeLabel):mesh.vEdges(f.edgeLabel);
        auto a=mesh.verts(edge.edgeToVert[0]).coord(),b=mesh.verts(edge.edgeToVert[1]).coord();
        if(f.horizontal?a.x>b.x:a.y>b.y)std::swap(a,b);
        return InterfaceNormalFluxView::FaceChart{a,f.center,b};
    };
    const auto integral=[chart=q.chart](int row,int half,bool fromEnd,long double offset,long double width){
        const auto f=chart(row);const auto a=half?f.midpoint:f.begin,b=half?f.end:f.midpoint;
        const bool horizontal=a.y==b.y;const double normal=horizontal?a.y:a.x;
        if(normal==-.5||normal==.5)return 0.L;
        const long double pi=std::acos(-1.L),length=horizontal?static_cast<long double>(b.x)-a.x:static_cast<long double>(b.y)-a.y;
        const long double anchor=horizontal?(fromEnd?b.x:a.x):(fromEnd?b.y:a.y);
        const long double center=std::fma(fromEnd?-length:length,offset+.5L*width,anchor);
        return (horizontal?1:-1)*3*std::cos(pi*normal)*std::sin(pi*center)*std::sin(.5L*pi*length*width)/pi;
    };
    q.integral=[integral](int row,int half,double a,double b,double s,double t){
        const long double span=static_cast<long double>(b)-a,width=span*(static_cast<long double>(t)-s);
        const bool fromEnd=a+b>1;
        return integral(row,half,fromEnd,fromEnd?std::fma(span,1-static_cast<long double>(t),1-static_cast<long double>(b)):std::fma(span,s,a),width);
    };
    q.phaseIntegral=[integral](int row,int half,const LevelSetMappedSegment& phase){
        long double value=0;for(const auto& r:phase.ranges())if(r.negative)value+=integral(row,half,r.fromEnd,r.offset,r.length);return value;
    };
    q.rangeIntegral=[integral](int row,int half,double a,double b,bool fromEnd,long double offset,long double length){
        const long double span=static_cast<long double>(b)-a;
        return integral(row,half,fromEnd,std::fma(span,offset,fromEnd?1-static_cast<long double>(b):static_cast<long double>(a)),span*length);
    };
    q.density=[chart=q.chart](int row,int half,bool fromEnd,long double offset){
        const auto f=chart(row);const auto a=half?f.midpoint:f.begin,b=half?f.end:f.midpoint;
        const bool horizontal=a.y==b.y;const double normal=horizontal?a.y:a.x;
        if(normal==-.5||normal==.5)return 0.L;
        const long double pi=std::acos(-1.L),length=horizontal?static_cast<long double>(b.x)-a.x:static_cast<long double>(b.y)-a.y;
        const long double anchor=horizontal?(fromEnd?b.x:a.x):(fromEnd?b.y:a.y);
        return (horizontal?1:-1)*1.5L*length*std::cos(pi*normal)*std::sin(pi*std::fma(fromEnd?-length:length,offset,anchor));
    };
    q.densityMagnitudeBound=[chart=q.chart](int row,int half){const auto f=chart(row);return 1.5L*(half?f.end-f.midpoint:f.midpoint-f.begin).mod();};
    q.roots=[chart=q.chart](int row,int half){
        const auto f=chart(row);const auto a=half?f.midpoint:f.begin,b=half?f.end:f.midpoint;
        const bool horizontal=a.y==b.y;const double normal=horizontal?a.y:a.x;
        Quantum::Geometry::UnitPolynomialRoots result;
        if(normal==-.5||normal==.5){result.identicallyZero=true;return result;}
        const double begin=horizontal?a.x:a.y,end=horizontal?b.x:b.y;
        if(begin<=0&&end>=0){const double root=-begin/(end-begin);
            result.roots.push_back({root,root==0||root==1?root:std::nextafter(root,0.),root==0||root==1?root:std::nextafter(root,1.)});}
        return result;
    };
    return std::make_shared<const InterfaceNormalFluxView>(std::move(q));
}
}
ReversibleVortex::ReversibleVortex(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,double period,Profile profile):period_(period),profile_(profile)
{
    domain.validateMeshIdentity(mesh);
    if(profile_==Profile::SlipReversal){
        if(domain.axisymmetric()||mesh.minCoord.x!=-.5||mesh.minCoord.y!=-.5||mesh.maxCoord.x!=.5||mesh.maxCoord.y!=.5||!std::isfinite(period)||period<=0)
            throw std::invalid_argument("slip reversal requires the centered unit XY square and positive period");
        spatial_=slipReversalFlux(mesh,domain);return;
    }
    if(domain.axisymmetric() || mesh.minCoord.x!=0 || mesh.minCoord.y!=0 || mesh.maxCoord.x!=1 || mesh.maxCoord.y!=1
        || !std::isfinite(period) || period<=0) throw std::invalid_argument("reversible vortex requires a unit XY square and positive period");
    InterfaceNormalFluxView::Queries q;
    q.chart=[&mesh,&domain](int row) {
        const auto f=domain.face(mesh,row);const auto& edge=f.horizontal?mesh.hEdges(f.edgeLabel):mesh.vEdges(f.edgeLabel);
        auto a=mesh.verts(edge.edgeToVert[0]).coord(),b=mesh.verts(edge.edgeToVert[1]).coord();
        if(f.horizontal?a.x>b.x:a.y>b.y) std::swap(a,b);
        return InterfaceNormalFluxView::FaceChart{a,f.center,b};
    };
    const auto integral=[chart=q.chart](int row,int half,bool fromEnd,long double offset,long double width) {
        const auto f=chart(row);const auto a=half?f.midpoint:f.begin,b=half?f.end:f.midpoint;
        const bool horizontal=a.y==b.y;const double normal=horizontal?a.y:a.x;
        const long double length=horizontal?static_cast<long double>(b.x)-a.x:static_cast<long double>(b.y)-a.y;
        const long double anchor=horizontal?(fromEnd?b.x:a.x):(fromEnd?b.y:a.y);
        if(normal==0 || normal==1) return 0.L;
        const auto transverse=std::sin(std::acos(-1.L)*std::min(normal,1-normal));
        return (horizontal?-1:1)*transverse*transverse*vortexSineOffset(anchor,fromEnd?-length:length,offset+.5L*width)
            *std::sin(std::acos(-1.L)*length*width)/std::acos(-1.L);
    };
    q.integral=[integral](int r,int half,double a,double b,double s,double t) {
        const long double span=static_cast<long double>(b)-a,width=span*(static_cast<long double>(t)-s);
        const bool fromEnd=a+b>1;
        return integral(r,half,fromEnd,fromEnd?std::fma(span,1-static_cast<long double>(t),1-static_cast<long double>(b))
            :std::fma(span,s,a),width);
    };
    q.phaseIntegral=[integral](int r,int half,const LevelSetMappedSegment& phase) {
        long double sum=0;for(const auto& range:phase.ranges()) if(range.negative)
            sum+=integral(r,half,range.fromEnd,range.offset,range.length);
        return sum;
    };
    q.rangeIntegral=[integral](int r,int half,double a,double b,bool fromEnd,long double offset,long double length) {
        const long double span=static_cast<long double>(b)-a;
        return integral(r,half,fromEnd,std::fma(span,offset,fromEnd?1-static_cast<long double>(b):static_cast<long double>(a)),span*length);
    };
    q.density=[chart=q.chart](int r,int half,bool fromEnd,long double offset) {
        const auto f=chart(r);const auto a=half?f.midpoint:f.begin,b=half?f.end:f.midpoint;
        const bool horizontal=a.y==b.y;const double normal=horizontal?a.y:a.x;
        const long double length=horizontal?static_cast<long double>(b.x)-a.x:static_cast<long double>(b.y)-a.y;
        const long double anchor=horizontal?(fromEnd?b.x:a.x):(fromEnd?b.y:a.y);
        if(normal==0 || normal==1) return 0.L;
        const auto transverse=std::sin(std::acos(-1.L)*std::min(normal,1-normal));
        return (horizontal?-1:1)*transverse*transverse*vortexSineOffset(anchor,fromEnd?-length:length,offset)*length;
    };
    q.densityMagnitudeBound=[chart=q.chart](int r,int half) {
        const auto f=chart(r);const auto a=half?f.midpoint:f.begin,b=half?f.end:f.midpoint;
        // Each sine factor has magnitude at most one.
        return static_cast<long double>((b-a).mod());
    };
    q.roots=[chart=q.chart](int r,int half) {
        const auto f=chart(r);const auto a=half?f.midpoint:f.begin,b=half?f.end:f.midpoint;
        const bool horizontal=a.y==b.y;const auto normal=horizontal?a.y:a.x;
        Quantum::Geometry::UnitPolynomialRoots result;
        if(normal==0 || normal==1) {result.identicallyZero=true;return result;}
        const long double first=horizontal?a.x:a.y,last=horizontal?b.x:b.y;
        for(long double zero:{0.L,.5L,1.L}) if(zero>=first && zero<=last) {
            const double point=static_cast<double>((zero-first)/(last-first));
            if(zero==first || zero==last) result.roots.push_back({point,point,point});
            else result.roots.push_back({point,std::nextafter(point,0.),std::nextafter(point,1.)});
        }
        return result;
    };
    spatial_=std::make_shared<const InterfaceNormalFluxView>(std::move(q));
}
double ReversibleVortex::amplitude(double time) const
{
    if(!std::isfinite(time)) throw std::invalid_argument("nonfinite prescribed vortex time");
    if(profile_==Profile::SlipReversal)return time==0||time==period_*.5||time==period_?0:std::sin(2*std::acos(-1.)*time/period_);
    return time==period_*.5?0:std::cos(std::acos(-1.)*time/period_);
}
Vector2D ReversibleVortex::sample(Vector2D p,double time) const
{
    if(profile_==Profile::SlipReversal){
        const double pi=std::acos(-1.),cx=p.x==-.5||p.x==.5?0:std::cos(pi*p.x),cy=p.y==-.5||p.y==.5?0:std::cos(pi*p.y);
        return Vector2D{-1.5*cx*std::sin(pi*p.y),1.5*std::sin(pi*p.x)*cy}*amplitude(time);
    }
    if(p.x==0 || p.x==1 || p.y==0 || p.y==1) return {0,0};
    const auto sx=std::sin(std::acos(-1.L)*std::min(p.x,1-p.x));
    const auto sy=std::sin(std::acos(-1.L)*std::min(p.y,1-p.y));
    return Vector2D{static_cast<double>(sx*sx*vortexSine(p.y)),static_cast<double>(-sy*sy*vortexSine(p.x))}*amplitude(time);
}
InterfaceNormalFluxView ReversibleVortex::normalFlux(double time) const
{return InterfaceNormalFluxView::scaled(spatial_,amplitude(time));}
LevelSetBoundaryTrace::VelocityJet ReversibleVortex::sampleJet(Vector2D p,double time) const
{
    if(profile_==Profile::SlipReversal){
        const double pi=std::acos(-1.),factor=1.5*pi*amplitude(time);
        const double cx=p.x==-.5||p.x==.5?0:std::cos(pi*p.x),cy=p.y==-.5||p.y==.5?0:std::cos(pi*p.y);
        const double diagonal=factor*std::sin(pi*p.x)*std::sin(pi*p.y),offDiagonal=factor*cx*cy;
        return {sample(p,time),{diagonal,-offDiagonal,offDiagonal,-diagonal}};
    }
    const auto pi=std::acos(-1.L);const double factor=amplitude(time);
    const auto sx=std::sin(pi*std::min(p.x,1-p.x)),sy=std::sin(pi*std::min(p.y,1-p.y));
    const auto xy=pi*vortexSine(p.x)*vortexSine(p.y)*factor;
    return {sample(p,time),{static_cast<double>(xy),static_cast<double>(2*pi*sx*sx*std::cos(2*pi*p.y)*factor),
        static_cast<double>(-2*pi*sy*sy*std::cos(2*pi*p.x)*factor),static_cast<double>(-xy)}};
}

AxialColumn::AxialColumn(double radius,double peakSpeed) : radius_(radius),speed_(peakSpeed)
{
    if(!std::isfinite(radius_) || radius_<=0 || !std::isfinite(speed_))
        throw std::invalid_argument("prescribed axial column requires positive radius and finite speed");
}
Vector2D AxialColumn::sample(const Vector2D& point) const
{
    return {0,point.x>=radius_?0:speed_*std::pow(1-point.x*point.x/(radius_*radius_),5)};
}

InterfaceNormalFluxView AxialColumn::normalFlux(const QuadtreeMesh& mesh,const NativeLeafDomain& domain) const
{
    domain.validateMeshIdentity(mesh);
    if (!domain.axisymmetric()) throw std::invalid_argument("axial column flux requires full RZ measure");
    const double radius=radius_, speed=speed_;
    InterfaceNormalFluxView::Queries q;
    q.chart=[&](int row) {
        const auto face=domain.face(mesh,row);
        const auto& edge=face.horizontal?mesh.hEdges(face.edgeLabel):mesh.vEdges(face.edgeLabel);
        auto a=mesh.verts(edge.edgeToVert[0]).coord(),b=mesh.verts(edge.edgeToVert[1]).coord();
        if(face.horizontal?a.x>b.x:a.y>b.y)std::swap(a,b);
        return InterfaceNormalFluxView::FaceChart{a,face.center,b};
    };
    const auto ends=[chart=q.chart](int row,int half) {
        if(half!=0 && half!=1) throw std::out_of_range("prescribed flux requires half face 0 or 1");
        const auto face=chart(row);
        return std::array<Vector2D,2>{{half==0?face.begin:face.midpoint,half==0?face.midpoint:face.end}};
    };
    // Exact antiderivative in w=1-r^2/R^2, factored to retain small ranges.
    // Integral = pi*A*R^2*(w0^6-w1^6)/6; no primitive subtraction.
    const auto integral=[radius,speed](long double r0,long double r1,long double width) {
        if(r0>=radius && r1>=radius)return 0.0L;
        const long double w0=1-r0*r0/(radius*radius),w1=1-r1*r1/(radius*radius);
        long double powers=0;
        for(int k=0;k<6;++k)powers+=std::pow(w0,5-k)*std::pow(w1,k);
        return std::acos(-1.0L)*speed/6*width*(r0+r1)*powers;
    };
    q.integral=[ends,integral](int row,int half,double a,double b,double begin,double end) {
        const auto edge=ends(row,half);
        const long double dr=static_cast<long double>(edge[1].x)-edge[0].x;
        const long double span=static_cast<long double>(b)-a;
        const long double r0=edge[0].x+dr*(a+span*begin),r1=edge[0].x+dr*(a+span*end);
        return integral(r0,r1,dr*span*(static_cast<long double>(end)-begin));
    };
    q.phaseIntegral=[ends,integral](int row,int half,const LevelSetMappedSegment& segment) {
        const auto edge=ends(row,half);
        const long double dr=static_cast<long double>(edge[1].x)-edge[0].x;
        long double total=0;
        for(const auto& range:segment.ranges())if(range.negative)
        {
            const long double start=range.fromEnd?edge[1].x-dr*range.offset:edge[0].x+dr*range.offset;
            const long double width=dr*range.length;
            total+=range.fromEnd?integral(start-width,start,width):integral(start,start+width,width);
        }
        return total;
    };
    q.density=[ends,radius,speed](int row,int half,bool fromEnd,long double offset) {
        const auto edge=ends(row,half);
        const long double dr=static_cast<long double>(edge[1].x)-edge[0].x;
        const long double r=fromEnd?edge[1].x-dr*offset:edge[0].x+dr*offset;
        return r>=radius?0.0L:2*std::acos(-1.0L)*speed*dr*r*std::pow(1-r*r/(radius*radius),5);
    };
    q.densityMagnitudeBound=[ends,radius,speed](int row,int half){
        const auto edge=ends(row,half);const long double lo=edge[0].x,hi=edge[1].x,width=hi-lo;
        if(width==0||lo>=radius||speed==0)return 0.L;
        // r(1-r^2/R^2)^5 has its sole interior maximum at R/sqrt(11).
        const long double r=std::max(lo,std::min(hi,static_cast<long double>(radius)/std::sqrt(11.L)));
        const long double value=2*std::acos(-1.L)*std::abs(speed)*width*r*std::pow(1-r*r/(radius*radius),5);
        return value*(1+32*std::numeric_limits<double>::epsilon());
    };
    q.roots=[ends,radius,speed](int row,int half) {
        const auto edge=ends(row,half);const double a=edge[0].x,b=edge[1].x;
        Quantum::Geometry::UnitPolynomialRoots roots;
        if (speed==0 || a==b || a>=radius) {roots.identicallyZero=true;return roots;}
        // The actual profile has a zero interval beyond R. It is not a
        // polynomial on a straddling half-face; require a resolved native
        // break instead of inventing one isolated root for that interval.
        if (a<0 || b>radius) throw std::invalid_argument("axial column root chart straddles its radial support");
        if (a==0) roots.roots.push_back({0,0,0});
        if (b==radius) roots.roots.push_back({1,1,1});
        return roots;
    };
    // Copies retain ONE spatial authority, enabling the existing exact
    // separable time interpolation without expanding this degree-11 trace
    // into the lower-degree generic polynomial representation.
    return InterfaceNormalFluxView::scaled(std::make_shared<const InterfaceNormalFluxView>(std::move(q)),1);
}

void initializeParabolicCap(const QuadtreeMesh& mesh,const InterfaceTransportDomain& transport,
    Vector2D pin,double height,InterfaceState& destination)
{
    if (destination.stage() != InterfaceState::Stage::Empty)
        throw std::invalid_argument("prescribed interface requires a fresh state");
    auto& state = InterfaceStateAccess::get(destination);
    const auto& domain=transport.geometry();
    domain.validateMeshIdentity(mesh);
    if(!domain.uniformLeaves() || !domain.axisymmetric() || mesh.minCoord.x!=0
        || !std::isfinite(pin.x) || !std::isfinite(pin.y) || pin.x<=0
        || !std::isfinite(height) || height<=0 || state.geometry() || state.panels()
        || state.freeLevelSet.size()!=mesh.vertIndexNum())
        throw std::invalid_argument("parabolic cap requires uniform RZ, a positive radius/height and fresh state");
    InterfacePreparation preparation(state);
    state.meshIdentity = domain.meshIdentity();
    Quantum::Geometry::ActiveHermitePatch patch;
    patch.activeCells.assign(mesh.cellIndexNum(),0);
    const double nan=std::numeric_limits<double>::quiet_NaN();
    patch.nodes.assign(mesh.vertIndexNum(),{nan,nan,nan,nan});
    state.activeVertex.assign(mesh.vertIndexNum(),0);
    for(int c:domain.cellIndexByRow())
    {
        patch.activeCells[c]=1;
        for(int label:mesh.cells(mesh.cellIndexToLabel[c]).cellToVert)
            state.activeVertex[mesh.verts(label).index()]=1;
    }
    int pinned=-1;
    for(int i=0;i<mesh.vertIndexNum();++i)if(state.activeVertex[i])
    {
        const auto p=mesh.verts(mesh.vertIndexToLabel[i]).coord();
        const double value=p.y-pin.y-height*(1-p.x*p.x/(pin.x*pin.x));
        patch.nodes[i]={value,2*height*p.x/(pin.x*pin.x),1,0};
        state.freeLevelSet.value[i]=value;
        if(p.x==pin.x && p.y==pin.y)pinned=i;
    }
    if(pinned<0)throw std::invalid_argument("pinned point is not an active mesh vertex");
    int wet=-1;
    for(const auto face:domain.faces(mesh))
    {
        if(transport.boundary(face.nativeRow)!=InterfaceTransportDomain::Boundary::StationaryWall)continue;
        const auto& edge=face.horizontal?mesh.hEdges(face.edgeLabel):mesh.vEdges(face.edgeLabel);
        const int a=mesh.verts(edge.edgeToVert[0]).index(),b=mesh.verts(edge.edgeToVert[1]).index();
        const int other=a==pinned?b:b==pinned?a:-1;
        if(other<0 || patch.nodes[other].value>=0)continue;
        if(face.horizontal || wet>=0)throw std::invalid_argument("parabolic cap requires one vertical wet wall at pin");
        wet=face.negativeCell>=0?face.negativeCell:face.positiveCell;
    }
    if(wet<0)throw std::invalid_argument("parabolic cap has no wet wall at pin");
    patch.cornerSectors.push_back({wet,pinned,{0,2*height/pin.x,0,0},-1});
    state.pinnedCoordinates={pin};state.pinnedVertexIndices={pinned};
    state.scalarPatch=patch;
    auto source=std::make_shared<const LevelSetHermiteGeometry::Snapshot>(mesh,patch,true,CartesianGridScope::UniformLeaves);
    auto panels=std::make_unique<LevelSetInterfacePanels>(mesh,source,true,state.pinnedCoordinates,CartesianGridScope::UniformLeaves);
    InterfaceMotion::validatePinnedToAxisChain(*panels,pin);
    state.retainGeometry(std::move(source),std::move(panels));
    preparation.complete(InterfaceState::Stage::Ready);
}
}
