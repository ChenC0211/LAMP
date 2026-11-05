#include "Modules/Hydrodynamics/HydroKepConvection.h"
#include "Modules/Hydrodynamics/Detail/HydroFiniteContraction.h"
#include "Core/Geometry/Detail/FiniteVolumeClosure.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <Eigen/SVD>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseQR>
#include <Eigen/Jacobi>
#include "Modules/Hydrodynamics/HydroSpatialAdvector.h"
#include "Modules/Hydrodynamics/Detail/HydroRegularAdvectorGeometry.h"
#include <algorithm>
#include <array>
#include <cstdint>
#include <cfenv>
#include <cerrno>
#include <cstring>
#include <chrono>
#include <atomic>
#include <exception>
#include <type_traits>
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
#include <omp.h>
#include <float.h>
#endif
#if defined(_M_X64) || defined(__SSE2__)
#include <xmmintrin.h>
#endif
#include <iomanip>
#include <sstream>
#include <cmath>
#include <limits>
#include <map>
#include <mutex>
#include <numeric>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <tuple>

namespace LAMP::Hydrodynamics
{
HydroKineticHodge::Matrix HydroKepConvection::modalBulkTransportMatrix(const QuadtreeMesh& mesh,
    const HydroKineticHodge& hodge,const HydroMaterialFlux& flux)
{
    hodge.validateFor(mesh,hodge.geometryIdentity(),hodge.timeLayer());
    flux.validateFor(mesh,hodge.materialMass().geometry(),hodge.geometryIdentity(),flux.advectorIdentity(),hodge.timeLayer());
    const auto& source=hodge.reconstruction();
    if(!source.sameIdentity(flux.phaseRate().spatialAdvector().reconstruction())
        ||flux.liquidDensity()!=hodge.materialMass().liquidDensity()||flux.gasDensity()!=hodge.materialMass().gasDensity())
        throw std::invalid_argument("[HydroKepConvection] different modal source/material owner");
    const auto& regular=source.regularGeometry();const int nu=hodge.layout().velocityCount();
    const int n=nu+static_cast<int>(source.tangentialModes().size());
    std::vector<std::optional<Detail::HydroRegularCellStencil>> stencils(hodge.layout().pressureCount());
    std::map<std::pair<int,int>,long double> entries;
    for(const auto& q:hodge.strainQuadrature()){
        const auto f=flux.physicalFlux(q.pressure,q.point,q.liquid);
        if(f.x==0&&f.y==0)continue;
        if(!stencils[q.pressure])stencils[q.pressure]=regular.cellStencil(mesh,q.pressure);
        struct Row{int column;Vector2D value,transport;};std::vector<Row> rows;
        for(const auto& w:regular.cellJetWeights(*stencils[q.pressure],q.point))
            rows.push_back({w.velocityDof,{w.jet[0],w.jet[1]},
                {f.x*w.jet[2]+f.y*w.jet[3],f.x*w.jet[4]+f.y*w.jet[5]}});
        const int mode=source.tangentialModeForPressure(q.pressure);
        if(mode>=0){const auto j=source.tangentialModeJet(mode,q.point,q.liquid);
            rows.push_back({nu+mode,j.value,{f.x*j.derivative.xx+f.y*j.derivative.xy,
                f.x*j.derivative.yx+f.y*j.derivative.yy}});}
        for(const auto& a:rows)for(const auto& b:rows)
            entries[{a.column,b.column}]+=q.physicalWeight*(static_cast<long double>(a.value.x)*b.transport.x
                +static_cast<long double>(a.value.y)*b.transport.y);
    }
    std::vector<Eigen::Triplet<double,Eigen::Index>> triplets;
    for(const auto& e:entries){const double value=static_cast<double>(e.second);
        if(!std::isfinite(value))throw std::runtime_error("[HydroKepConvection] nonfinite modal bulk transport entry");
        if(value!=0)triplets.emplace_back(e.first.first,e.first.second,value);}
    HydroKineticHodge::Matrix result(n,n);result.setFromTriplets(triplets.begin(),triplets.end());result.makeCompressed();return result;
}
namespace
{
struct OpenApplyPlan;
constexpr long double tolerance=256.L*std::numeric_limits<double>::epsilon();
void require(bool ok,const char* why)
{if(!ok)throw std::invalid_argument(std::string("[HydroKepConvection] ")+why);}
long double finite(long double x,const char* why)
{if(Detail::hydroFiniteRangeEnabled()&&!std::isfinite(x))throw Detail::HydroFiniteRange{};require(std::isfinite(x),why);return x;}
struct ProductDiagnosticContext
{
    const char* operation="unbound";
    std::array<std::uint64_t,3> vofSources{};
    LAMP::TimeLayer clock;
};
thread_local const ProductDiagnosticContext* activeProductDiagnostic=nullptr;
struct ProductDiagnosticScope
{
    ProductDiagnosticContext value;
    const ProductDiagnosticContext* previous;
    ProductDiagnosticScope(const HydroMaterialSources& sources,LAMP::TimeLayer clock,const char* operation)
        :previous(activeProductDiagnostic)
    {
        value.operation=operation;value.clock=clock;
        for(std::size_t i=0;i<value.vofSources.size();++i)value.vofSources[i]=sources.vofIds[i].value();
        activeProductDiagnostic=&value;
    }
    ~ProductDiagnosticScope(){activeProductDiagnostic=previous;}
    ProductDiagnosticScope(const ProductDiagnosticScope&)=delete;
    ProductDiagnosticScope& operator=(const ProductDiagnosticScope&)=delete;
};
long double checkedProduct(long double a,long double b,const char* expressionA,const char* expressionB,
    const char* caller,int line)
{
    // The original multiply, finite guard and nonzero criterion are unchanged.
    const auto raw=a*b;Detail::hydroFiniteProductRange(a,b,raw);
    const auto v=finite(raw,"nonfinite convection product");
    if(!(a==0 || b==0 || v!=0)){
        int exponentA=0,exponentB=0;
        const auto mantissaA=std::frexp(a,&exponentA),mantissaB=std::frexp(b,&exponentB);
        std::ostringstream detail;detail<<std::scientific<<std::setprecision(std::numeric_limits<long double>::max_digits10)
            <<"[HydroKepConvection] unrepresentable convection product"
            <<" caller="<<caller<<" line="<<line<<" expression_a=["<<expressionA<<"] expression_b=["<<expressionB<<']'
            <<" a="<<a<<" b="<<b<<" result="<<v
            <<" a_class="<<std::fpclassify(a)<<" b_class="<<std::fpclassify(b)<<" result_class="<<std::fpclassify(v)
            <<" mantissa_a="<<mantissaA<<" exponent_a="<<exponentA
            <<" mantissa_b="<<mantissaB<<" exponent_b="<<exponentB<<" exponent_sum="<<(exponentA+exponentB)
            <<" digits="<<std::numeric_limits<long double>::digits<<" min_exponent="<<std::numeric_limits<long double>::min_exponent
            <<" min_normal="<<std::numeric_limits<long double>::min()<<" denorm_min="<<std::numeric_limits<long double>::denorm_min()
            <<" sizeof_long_double="<<sizeof(long double);
        if(activeProductDiagnostic){const auto& context=*activeProductDiagnostic;const auto& clock=context.clock;
            detail<<" context="<<context.operation<<" vof_source_begin="<<context.vofSources[0]
                <<" vof_source_midpoint="<<context.vofSources[1]<<" vof_source_end="<<context.vofSources[2]
                <<" step="<<clock.step<<" begin="<<clock.begin<<" end="<<clock.end
                <<" position="<<static_cast<int>(clock.position);
        }else detail<<" context=unbound";
        detail<<" a_hex="<<std::hexfloat<<a<<" b_hex="<<b<<" result_hex="<<v
            <<" mantissa_a_hex="<<mantissaA<<" mantissa_b_hex="<<mantissaB;
        throw std::invalid_argument(detail.str());
    }
    return v;
}
// Translation-unit-local diagnostic arguments only. All existing product
// operands/call sites retain their arithmetic and evaluation as before.
#define product(a,b) checkedProduct((a),(b),#a,#b,__func__,__LINE__)
long double identity(long double difference,long double scale,const char* why)
{
    finite(difference,why);finite(scale,why);require(scale>=0 && (scale>0 || difference==0),why);
    const auto error=scale>0?std::abs(difference)/scale:0;require(error<=tolerance,why);return error;
}
int component(const HydroDofLayout& l,int u){return l.normalAxis(u)==HydroDofLayout::NormalAxis::First?0:1;}
bool internal(const HydroMomentumSegment& g){return g.negativeVelocity>=0 && g.positiveVelocity>=0;}
int owner(const HydroMomentumSegment& g){return g.negativeVelocity>=0?g.negativeVelocity:g.positiveVelocity;}
Vector2D center(const HydroMomentumSegment& g,bool rz)
{
    auto p=g.begin+(g.end-g.begin)*.5;
    if(rz && g.normalAxis==1){
        const long double a=g.begin.x,h=static_cast<long double>(g.end.x)-a,d=a+h/2;
        require(d>0,"invalid radial segment centroid");
        p.x=static_cast<double>(a+h*(a/(2*d)+h/(3*d)));
    }
    require(std::isfinite(p.x) && std::isfinite(p.y),"nonfinite segment centroid");return p;
}
void insert(std::vector<int>& values,int v)
{if(std::find(values.begin(),values.end(),v)==values.end())values.push_back(v);}
std::vector<HydroConvectionWeight> reconstruct(const HydroMomentumSegment& g,Vector2D point,
    const HydroMomentumGeometry& geometry,const std::vector<NativeLeafFace>& faces,
    const std::vector<std::vector<int>>& cellFaces,const std::vector<std::vector<int>>& neighbors)
{
    const auto& l=geometry.layout();std::vector<int> roots,patch,support;
    for(int u:{g.negativeVelocity,g.positiveVelocity})for(const auto& c:geometry.volumes()[u].overlaps)insert(roots,c.pressure);
    patch=roots;for(int p:roots)for(int n:neighbors[p])insert(patch,n);
    for(int p:patch)for(int u:cellFaces[p])if(l.normalAxis(u)==l.normalAxis(g.negativeVelocity))insert(support,u);
    std::sort(support.begin(),support.end());require(!support.empty(),"empty common trace support");
    const auto& a=geometry.volumes()[g.negativeVelocity];const auto& b=geometry.volumes()[g.positiveVelocity];
    const Vector2D width{std::max(a.upper.x-a.lower.x,b.upper.x-b.lower.x),std::max(a.upper.y-a.lower.y,b.upper.y-b.lower.y)};
    require(width.x>0 && width.y>0,"invalid common trace scale");
    Eigen::MatrixXd constraints(3,support.size());Eigen::VectorXd w=Eigen::VectorXd::Zero(support.size());
    for(Eigen::Index j=0;j<constraints.cols();++j){const int u=support[static_cast<std::size_t>(j)];
        constraints(0,j)=1;constraints(1,j)=(faces[u].center.x-point.x)/width.x;constraints(2,j)=(faces[u].center.y-point.y)/width.y;
        if(u==g.negativeVelocity || u==g.positiveVelocity)w[j]=.5;
    }
    require(constraints.allFinite(),"nonfinite common trace constraints");
    Eigen::Vector3d target;target<<1,0,0;const Eigen::Vector3d residual=target-constraints*w;
    // Preserve ordinary center averaging exactly where it already reproduces
    // the physical common point. Else solve only this identifiable functional.
    if(residual.cwiseAbs().maxCoeff()>8*std::numeric_limits<double>::epsilon()){
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(constraints,Eigen::ComputeThinU|Eigen::ComputeThinV);svd.setThreshold(1e-12);
        w+=svd.solve(residual);
    }
    require(w.allFinite(),"nonfinite common trace weights");
    const Eigen::Vector3d difference=constraints*w-target;
    for(int k=0;k<3;++k){long double gross=std::abs(target[k]);for(Eigen::Index j=0;j<w.size();++j)gross+=std::abs(constraints(k,j)*w[j]);
        identity(difference[k],gross,"unidentifiable affine common trace");}
    std::vector<HydroConvectionWeight> result;
    for(Eigen::Index j=0;j<w.size();++j)if(w[j]!=0)result.push_back({support[static_cast<std::size_t>(j)],w[j]});
    return result;
}
// Shared sparse native pairs, with block constraints on the FINAL operator.
// Local 1/6-row SVD blocks provide a minimum-change patch preconditioner;
// the transient multipliers are discarded after preparing the shared pairs.
struct MomentEntry {std::size_t edge;std::array<double,6> a{};};
struct MomentBlock
{
    int offset=0,count=0;std::vector<MomentEntry> entries;Eigen::MatrixXd whitening;
};
std::array<double,6> polynomial(Vector2D a,Vector2D b,double h)
{
    const double x=(b.x-a.x)/h,y=(b.y-a.y)/h;
    return {1,x,y,x*x,x*y,y*y};
}
void accumulate(double value,double& sum,double& compensation)
{const double term=value-compensation,next=sum+term;compensation=(next-sum)-term;sum=next;}
long double dot(const std::vector<double>& a,const std::vector<double>& b)
{double value=0,compensation=0;for(std::size_t i=0;i<a.size();++i)accumulate(a[i]*b[i],value,compensation);return finite(value,"nonfinite joint moment dot product");}
double maximum(const std::vector<double>& a)
{double value=0;for(double x:a)value=std::max(value,std::abs(x));return value;}
struct NativeTopology
{
    bool paired=false;Vector2D period{};
    std::vector<int> partner,representative,count;
    std::vector<std::vector<int>> roots;
};
Vector2D displacement(Vector2D a,Vector2D b,const NativeTopology& topology)
{
    auto d=b-a;
    if(topology.period.x>0)d.x-=topology.period.x*std::round(d.x/topology.period.x);
    if(topology.period.y>0)d.y-=topology.period.y*std::round(d.y/topology.period.y);
    return d;
}
std::array<double,6> polynomial(Vector2D a,Vector2D b,double h,const NativeTopology& topology)
{return polynomial({},displacement(a,b,topology),h);}
NativeTopology topologyFor(const HydroDualMaterial& m,const std::vector<NativeLeafFace>& faces,
    std::vector<std::vector<int>>& neighbors,const std::vector<HydroConvectionPeriodicPair>* pairs)
{
    NativeTopology t;t.paired=pairs!=nullptr;const int n=m.layout().velocityCount();const auto& segments=m.geometry().segments();
    t.partner.assign(segments.size(),-1);t.representative.resize(n);std::iota(t.representative.begin(),t.representative.end(),0);
    const auto root=[&](int i){while(t.representative[i]!=i)i=t.representative[i];return i;};
    if(pairs)for(const auto& p:*pairs){const int g=p.firstSegment,h=p.secondSegment;
        require(g>=0 && h>=0 && g!=h && static_cast<std::size_t>(g)<segments.size() && static_cast<std::size_t>(h)<segments.size(),"invalid periodic segment indices");
        const auto& a=segments[g];const auto& b=segments[h];require(!internal(a) && !internal(b) && t.partner[g]<0 && t.partner[h]<0,"repeated or internal periodic segment");
        const int i=owner(a),j=owner(b);require(component(m.layout(),i)==component(m.layout(),j) && a.normalAxis==b.normalAxis && (a.negativeVelocity>=0)!=(b.negativeVelocity>=0),"incompatible periodic segment orientation");
        require(std::isfinite(p.translation.x) && std::isfinite(p.translation.y),"nonfinite periodic translation");
        const double length=a.normalAxis==0?std::abs(p.translation.x):std::abs(p.translation.y);
        require(length>0 && (a.normalAxis==0?p.translation.y==0:p.translation.x==0),"periodic translation must be normal to segments");
        require(!m.layout().domain().axisymmetric() || a.normalAxis==1,"radial periodicity changes RZ measure");
        double& declared=a.normalAxis==0?t.period.x:t.period.y;if(declared>0)identity(length-declared,length+declared,"different periodic lengths");else declared=length;
        const auto begin=b.begin+p.translation,end=b.end+p.translation;
        identity(std::hypot(begin.x-a.begin.x,begin.y-a.begin.y),length,"periodic segment begin mismatch");
        identity(std::hypot(end.x-a.end.x,end.y-a.end.y),length,"periodic segment end mismatch");
        identity(a.area-b.area,a.area+b.area,"periodic segment measure mismatch");
        const auto fa=(a.negativeVelocity>=0?1:-1)*m.faceRates()[g],fb=(b.negativeVelocity>=0?1:-1)*m.faceRates()[h];
        identity(fa+fb,std::abs(fa)+std::abs(fb),"periodic material rates mismatch");t.partner[g]=h;t.partner[h]=g;
        if(component(m.layout(),i)==a.normalAxis){
            const auto d=displacement(faces[i].center,faces[j].center,t);identity(std::hypot(d.x,d.y),length,"periodic normal samples do not coincide");
            const int ri=root(i),rj=root(j);t.representative[std::max(ri,rj)]=std::min(ri,rj);
            // These actual primal boundary faces close pressure adjacency.
            for(const auto& x:m.geometry().volumes()[i].overlaps)for(const auto& y:m.geometry().volumes()[j].overlaps){insert(neighbors[x.pressure],y.pressure);insert(neighbors[y.pressure],x.pressure);}
        }
    }
    if(pairs)for(std::size_t g=0;g<segments.size();++g)if(!internal(segments[g]) && m.faceRates()[g]!=0)
        require(t.partner[g]>=0,"unpaired nonzero material boundary rate");
    t.roots.resize(n);t.count.resize(n);
    for(int i=0;i<n;++i){const int r=root(i);t.representative[i]=r;++t.count[r];for(const auto& p:m.geometry().volumes()[i].overlaps)insert(t.roots[r],p.pressure);}
    return t;
}
// Translation-unit-private and limited to one enclosing public call. Callers
// first check the live mesh before reusing the sealed immutable owner/binding.
// No mutable mesh pointer or state value is retained.
class CheckedVofMaterialBatch
{
    friend class LAMP::Hydrodynamics::HydroKepConvection;
    friend class LAMP::Hydrodynamics::HydroConvectionBoundary;
    const void* dual_=nullptr;
    std::shared_ptr<const NativeMeshIdentity> binding_;
    explicit CheckedVofMaterialBatch(const HydroDualMaterial& material)
    {
        if(material.material().model()==HydroPhaseGeometry::Model::GeometricVof){
            dual_=&material.cells();binding_=material.layout().domain().meshIdentity();}
    }
public:
    static CheckedVofMaterialBatch validate(const QuadtreeMesh& mesh,const HydroDualMaterial& material,
        HydroMaterialSources source,LAMP::TimeLayer clock)
    {material.validateFor(mesh,source,clock);return CheckedVofMaterialBatch(material);}
    bool sameBinding(const NativeLeafDomain& domain) const noexcept
    {return dual_&&binding_==domain.meshIdentity();}
    bool matches(const HydroDualMaterial& material) const noexcept
    {
        // Exact immutable owner identity, not equal hashes/numeric snapshots.
        if(!dual_||dual_!=&material.cells()||!sameBinding(material.layout().domain())
            ||material.material().model()!=HydroPhaseGeometry::Model::GeometricVof)return false;
        const auto& interval=material.material();
        return sameBinding(material.geometry().layout().domain())
            &&sameBinding(interval.beginGeometry().vofPublication().domain())
            &&sameBinding(interval.midpointGeometry().vofPublication().domain())
            &&sameBinding(interval.endGeometry().vofPublication().domain());
    }
};
bool wallExecutionQualified() noexcept
{
    // Preserve both additional public HydroAdvectorBoundary FP predicates.
    // Failure calls its original validator for the original error ordering.
    if(std::fegetround()!=FE_TONEAREST)return false;
#if defined(_M_X64) || defined(__SSE2__)
    if((_mm_getcsr()&((1U<<15)|(1U<<6)))!=0)return false;
#endif
    return true;
}
void validateWallsAfterMaterial(const QuadtreeMesh& mesh,const HydroDualMaterial& material,
    const HydroAdvectorBoundary& walls,const CheckedVofMaterialBatch* checked)
{
    if(checked&&checked->matches(material)&&wallExecutionQualified()
        &&walls.momentumGeometry().sameIdentity(material.geometry())
        &&checked->sameBinding(walls.momentumGeometry().layout().domain()))return;
    walls.validateFor(mesh,material.geometry());
}
void validateConvectionBoundaryAfterMaterial(const QuadtreeMesh& mesh,const HydroDualMaterial& material,
    const HydroConvectionBoundary& boundary,HydroMaterialSources source,LAMP::TimeLayer clock,
    const CheckedVofMaterialBatch* checked)
{
    if(checked&&checked->matches(material)&&checked->matches(boundary.material()))return;
    boundary.validateFor(mesh,material,source,clock);
}
void validateVolumeBoundaryAfterMaterial(const QuadtreeMesh& mesh,const HydroDualMaterial& material,
    const Interface::GeometricVofAdvectorSlab::Boundary& boundary,const Interface::GeometricVofPublication& source,
    const CheckedVofMaterialBatch* checked)
{
    // This actual retained slab/boundary was qualified by the same material
    // Transport validator. Different bindings keep the standalone deep path.
    if(checked&&checked->matches(material)&&checked->sameBinding(boundary.domain()))return;
    boundary.validateFor(mesh,source);
}
const Interface::GeometricVofTransport& validateVofNativeSource(const QuadtreeMesh& mesh,
    const HydroDualMaterial& material,const CheckedVofMaterialBatch* checked=nullptr)
{
    const auto& interval=material.material();const auto& volume=interval.vofTransport();const auto& slab=volume.advector();
    if(!checked||!checked->matches(material))
        volume.validateFor(mesh,interval.beginGeometry().vofPublication(),slab,interval.liquidDensity(),interval.gasDensity());
    require(volume.source().geometryKind()==Interface::GeometricVofGeometryKind::LocalPlic,
        "typed VOF native pairing requires the actual LocalPlic source");
    const auto clock=interval.timeLayer(),mid=volume.midpoint().timeLayer(),end=slab.endTime();
    require(clock.step==mid.step&&clock.begin==mid.begin&&clock.end==mid.end&&clock.position==mid.position
        &&clock.position==Interface::GeometricVofTimeLayer::Position::Midpoint
        &&end.step==clock.step&&end.begin==clock.begin&&end.end==clock.end
        &&end.position==Interface::GeometricVofTimeLayer::Position::End,"different VOF native source interval");
    if(!checked||!checked->matches(material)){
        interval.midpointGeometry().validateVofFor(mesh,volume.midpoint().identity(),mid);
        interval.endGeometry().validateVofFor(mesh,volume.endpoint().identity(),end);
    }
    require(volume.sourceIdentity()==interval.beginGeometry().vofIdentity(),"different VOF native source geometry");
    require(material.layout().domain().meshIdentity()==volume.source().domain().meshIdentity()
        &&static_cast<std::size_t>(material.layout().pressureCount())==volume.source().domain().cellIndexByRow().size()
        &&material.layout().domain().cellIndexByRow()==volume.source().domain().cellIndexByRow()
        &&material.layout().domain().cellRowByIndex()==volume.source().domain().cellRowByIndex(),
        "VOF physical KEP requires the same actual selected native domain");
    return volume;
}
const Interface::GeometricVofTransport& validateVofOpenSource(const QuadtreeMesh& mesh,
    const HydroDualMaterial& material,const HydroAdvectorBoundary& walls,const CheckedVofMaterialBatch* checked=nullptr)
{
    const auto& volume=validateVofNativeSource(mesh,material,checked);const auto& slab=volume.advector();
    const auto* boundary=slab.boundary();require(boundary&&!boundary->allClosed(),"VOF prescribed flow requires actual F1 reservoir owner");
    validateVolumeBoundaryAfterMaterial(mesh,material,*boundary,volume.source(),checked);
    validateWallsAfterMaterial(mesh,material,walls,checked);
    using Mode=Interface::GeometricVofAdvectorSlab::Boundary::Mode;std::vector<int> actualClosed;
    const auto& layout=material.layout();
    for(int u=0;u<layout.velocityCount();++u){const int row=layout.nativeFace(u);const auto face=layout.domain().face(mesh,row);
        if(face.negativeCell>=0&&face.positiveCell>=0)continue;
        if(boundary->rules()[row]==Mode::Closed){actualClosed.push_back(u);
            require(slab.beginNormal()[row]==0&&slab.endNormal()[row]==0,"nonzero actual F1 Closed endpoint normal");}
    }
    require(actualClosed==walls.velocityDofs(),"VOF Closed wall ports differ from actual F1 boundary owner");return volume;
}
std::vector<int> impermeableRows(const QuadtreeMesh& mesh,const HydroDualMaterial& m,
    const HydroAdvectorBoundary& boundary,const CheckedVofMaterialBatch* checked=nullptr)
{
    const auto& l=m.layout();
    const bool vof=m.material().model()==HydroPhaseGeometry::Model::GeometricVof;
    require(vof||l.domain().uniformLeaves(),"impermeable FirstJet native pairing requires uniform XY or RZ leaves");
    require(vof||l.pressureCount()==mesh.cellIndexNum(),"impermeable FirstJet native moments require the complete native domain");
    validateWallsAfterMaterial(mesh,m,boundary,checked);
    std::vector<int> fixed;
    if(vof){
        require(boundary.completeBoundaryLines(),"VOF impermeable source requires complete physical wall lines");
        const auto& slab=validateVofNativeSource(mesh,m,checked).advector();
        for(int u=0;u<l.velocityCount();++u){const auto row=l.nativeFace(u);const auto face=l.domain().face(mesh,row);
            if(face.negativeCell<0||face.positiveCell<0){fixed.push_back(u);
                require(slab.beginNormal()[row]==0&&slab.endNormal()[row]==0,"nonzero VOF impermeable endpoint native mean");}}
    }else{
        const auto& actual=m.material().midpointRate().spatialAdvector();
        actual.validateFor(mesh,m.geometry(),actual.geometryIdentity(),m.timeLayer());
        require(actual.boundary()&&actual.boundary()->sameIdentity(boundary),"different actual impermeable source boundary");
        require(actual.reconstruction().tangentialModes().empty()&&actual.tangentialVelocity().empty(),
            "impermeable native moments require the base HYDRO source");
        for(int u=0;u<l.velocityCount();++u){const auto face=l.domain().face(mesh,l.nativeFace(u));
            if(face.negativeCell<0||face.positiveCell<0){fixed.push_back(u);
                require(actual.values().normalVelocity()[u]==0,"nonzero actual impermeable native mean");}}
    }
    require(!fixed.empty()&&fixed==boundary.velocityDofs(),"impermeable boundary is not the complete native normal wall set");
    const auto& segments=m.geometry().segments();
    for(std::size_t g=0;g<segments.size();++g)if(!internal(segments[g]))
        require(m.faceRates()[g]==0&&m.faceIntegrals()[g]==0,"nonzero impermeable exterior material flux");
    return fixed;
}
void validatePhysicalSourceBoundary(const QuadtreeMesh& mesh,const HydroDualMaterial& material,
    const HydroAdvectorBoundary& walls,bool closed,const CheckedVofMaterialBatch* checked=nullptr)
{
    if(material.material().model()==HydroPhaseGeometry::Model::GeometricVof){
        if(closed)(void)impermeableRows(mesh,material,walls,checked);
        else (void)validateVofOpenSource(mesh,material,walls,checked);
        return;
    }
    walls.validateFor(mesh,material.geometry());
    const auto& actual=material.material().midpointRate().spatialAdvector();
    require(actual.boundary()&&actual.boundary()->sameIdentity(walls),"different retained physical source boundary");
}
std::vector<int> prescribedFlowRows(const QuadtreeMesh& mesh,const HydroDualMaterial& m,
    const HydroAdvectorBoundary& walls,const HydroConvectionBoundary& inlet,
    const std::vector<int>& freeOutflow,bool vof)
{
    const auto& l=m.layout();
    if(vof){require(m.material().model()==HydroPhaseGeometry::Model::GeometricVof&&
            static_cast<std::size_t>(l.pressureCount())==m.material().vofTransport().source().domain().cellIndexByRow().size(),
            "VOF prescribed flow requires the actual selected native domain");
        require(inlet.hasVofEndpointVelocity(),"VOF inflow requires typed actual endpoint full vectors");}
    else {
        require(m.material().model()!=HydroPhaseGeometry::Model::GeometricVof,"LocalPlic KEP supports only the typed closed physical boundary");
        require(l.domain().uniformLeaves()&&l.pressureCount()==mesh.cellIndexNum(),
            "prescribed flow requires the complete FirstJet uniform XY/RZ domain");
    }
    walls.validateFor(mesh,m.geometry());inlet.validateFor(mesh,m,m.sources(),m.timeLayer());
    validatePhysicalSourceBoundary(mesh,m,walls,false);
    const auto normal=[&](int u,bool end){
        if(vof){const auto& slab=m.material().vofTransport().advector();const auto row=l.nativeFace(u);
            return end?slab.endNormal()[row]:slab.beginNormal()[row];}
        return m.material().midpointRate().spatialAdvector().values().normalVelocity()[u];
    };
    if(!vof){const auto& actual=m.material().midpointRate().spatialAdvector();
        actual.validateFor(mesh,m.geometry(),actual.geometryIdentity(),m.timeLayer());
        require(actual.reconstruction().tangentialModes().empty()&&actual.tangentialVelocity().empty(),
            "prescribed flow requires the base HYDRO source");}
    for(int u:walls.velocityDofs())require(normal(u,false)==0&&normal(u,true)==0,"nonzero actual prescribed-flow wall mean");
    std::vector<int> fixed;
    for(int u=0;u<l.velocityCount();++u){const auto f=l.domain().face(mesh,l.nativeFace(u));
        if(f.negativeCell<0 || f.positiveCell<0)fixed.push_back(u);}
    require(!fixed.empty(),"prescribed flow has no physical boundary faces");
    std::vector<bool> freeMask(l.velocityCount(),false);
    for(int u:freeOutflow){
        require(u>=0 && u<l.velocityCount() && !freeMask[u],
            "duplicate or invalid free outflow normal");
        const auto f=l.domain().face(mesh,l.nativeFace(u));
        require(std::binary_search(fixed.begin(),fixed.end(),u) && f.outwardSign!=0 && f.area>0,
            "free outflow normal is internal or zero-area axis");
        require(!std::binary_search(walls.velocityDofs().begin(),walls.velocityDofs().end(),u),
            "free outflow normal overlaps a prescribed wall");
        require(f.outwardSign*normal(u,false)>=0&&f.outwardSign*normal(u,true)>=0,
            "free outflow endpoint native normal points inward");
        freeMask[u]=true;
    }
    for(const auto& g:m.geometry().segments())if(!internal(g))
        require(g.primalVelocity>=0 && std::binary_search(fixed.begin(),fixed.end(),g.primalVelocity),
            "prescribed-flow exterior segment has no physical native port");
    const auto& segments=m.geometry().segments();
    for(std::size_t g=0;g<segments.size();++g){const auto& segment=segments[g];
        if(internal(segment)){require(inlet.traces()[g]==0,"nonzero unused internal inlet trace");continue;}
        const auto q=finite((segment.negativeVelocity>=0?1.L:-1.L)*m.faceRates()[g],
            "nonfinite prescribed-flow exterior mass rate");
        if(std::binary_search(walls.velocityDofs().begin(),walls.velocityDofs().end(),segment.primalVelocity))
            require(m.faceRates()[g]==0 && m.faceIntegrals()[g]==0,
                "nonzero material flux through prescribed-flow wall");
        if(freeMask[segment.primalVelocity]){
            const auto integral=finite((segment.negativeVelocity>=0?1.L:-1.L)*m.faceIntegrals()[g],
                "nonfinite free outflow exterior material integral");
            require(q>=0 && integral>=0,"incoming material flux through a free outflow");
        }
        if(q>=0)require(inlet.traces()[g]==0,"inlet trace supplied on a zero-flow or outflow segment");
    }
    fixed.erase(std::remove_if(fixed.begin(),fixed.end(),[&](int u){return freeMask[u];}),fixed.end());
    return fixed;
}
// The source already owns each native normal polynomial and its true XY/RZ
// measure. Consume those exact atoms, rather than fitting an outlet field or
// equating its point value with the native face mean.
void certifyFreeOutflowMeans(const HydroSpatialAdvector& actual,const std::vector<int>& freeOutflow,
    HydroConvectionMomentReport& report)
{
    const auto& reconstruction=actual.reconstruction();
    const auto& atoms=reconstruction.atoms();const auto& polynomials=actual.regularAtomPolynomials();
    for(int u:freeOutflow){
        long double rate=0,area=0,gross=0;
        require(!reconstruction.primalAtoms()[u].empty(),"free outflow has no actual source atoms");
        for(int id:reconstruction.primalAtoms()[u]){
            const auto& atom=atoms[id];const auto& polynomial=polynomials[id];
            const auto value=finite(polynomial.integral(false,0,polynomial.length),
                "nonfinite free outflow source mean");
            rate=finite(rate+value,"nonfinite free outflow source mean sum");
            area=finite(area+atom.area,"nonfinite free outflow source measure");
            gross=finite(gross+std::abs(value),"nonfinite free outflow source mean scale");
        }
        require(area>0,"free outflow has zero actual source measure");
        const auto target=product(area,actual.values().normalVelocity()[u]);
        const auto defect=identity(rate-target,gross+std::abs(target),
            "actual free outflow source mean differs from original native normal");
        report.maxNativeMeanDefect=std::max(report.maxNativeMeanDefect,defect);
    }
}
void certifyVofFreeOutflowMeans(const QuadtreeMesh&,const HydroDualMaterial&,
    const std::vector<int>&,HydroConvectionMomentReport&,const CheckedVofMaterialBatch* = nullptr);
// Uniform closed XY/RZ MAC: retain the physical two-point dual-face operator.
// F already includes the native physical measure; no Cartesian rescaling is applied.
// The canonical material integral remains authoritative; only its existing
// roundoff certificate is checked, never repaired by a second flux or target.
std::vector<HydroConvectionCoupling> physicalCenteredPairs(const HydroDualMaterial& m,
    const std::vector<bool>& fixed,HydroConvectionMomentReport& report,
    const std::vector<long double>* externalOutflow=nullptr,
    const Detail::HydroOpenConvectionSupport* support=nullptr,std::vector<long double>* sharedCoefficients=nullptr,
    std::vector<long double>* canonicalOutflow=nullptr)
{
    const int n=m.layout().velocityCount();require(fixed.size()==static_cast<std::size_t>(n),"different physical fixed-row dimensions");
    require(!externalOutflow || externalOutflow->size()==static_cast<std::size_t>(n),
        "different prescribed-flow exterior dimensions");
    report.physicalPairing=1;report.fixedRows=static_cast<int>(std::count(fixed.begin(),fixed.end(),true));
    report.freeRows=report.boundedRows=report.interiorRows=n-report.fixedRows;
    std::map<std::pair<int,int>,long double> pairs;std::vector<long double> gross(n),rowSum(n),rowNorm(n);
    if(support)sharedCoefficients->assign(support->pairs.size(),0);
    const auto& segments=m.geometry().segments();
    for(std::size_t g=0;g<segments.size();++g){const auto& segment=segments[g];const auto f=m.faceRates()[g];
        if(!internal(segment)){
            if(externalOutflow)gross[owner(segment)]+=std::abs(f);
            else require(f==0 && m.faceIntegrals()[g]==0,"nonzero physical closed exterior flux");
            continue;
        }
        ++report.physicalSegments;int i=segment.negativeVelocity,j=segment.positiveVelocity;
        require(i!=j && component(m.layout(),i)==component(m.layout(),j),"different physical pair components");
        gross[i]+=std::abs(f);gross[j]+=std::abs(f);auto coefficient=product(.5L,f);
        if(i>j){std::swap(i,j);coefficient=-coefficient;}
        auto& sum=support?(*sharedCoefficients)[support->segmentPairSlot[g]]:pairs[{i,j}];
        sum=finite(sum+coefficient,"nonfinite physical shared pair");
    }
    std::vector<HydroConvectionCoupling> result;
    const auto append=[&](int i,int j,long double k){if(k==0)return;result.push_back({i,j,k});
        rowSum[i]+=k;rowSum[j]-=k;rowNorm[i]+=std::abs(k);rowNorm[j]+=std::abs(k);};
    if(support){for(std::size_t id=0;id<support->pairs.size();++id){const auto& pair=support->pairs[id];append(pair.first,pair.second,(*sharedCoefficients)[id]);}}
    else for(const auto& entry:pairs)append(entry.first.first,entry.first.second,entry.second);
    report.physicalPairs=static_cast<int>(result.size());
    for(int i=0;i<n;++i){const auto canonical=finite(m.cells()[i].outwardIntegral/m.deltaTime(),"nonfinite physical canonical outflow");
        if(canonicalOutflow)(*canonicalOutflow)[i]=canonical;
        const auto internalOutflow=externalOutflow
            ? finite(canonical-(*externalOutflow)[i],"nonfinite physical internal outflow") : canonical;
        const auto difference=finite(rowSum[i]-.5L*internalOutflow,"nonfinite physical constant defect");
        const auto error=identity(difference,gross[i],"physical pair canonical constant certificate failed");
        report.constantError=std::max(report.constantError,error);
        if(gross[i]==0){++report.zeroGrossRows;require(rowNorm[i]==0 && internalOutflow==0,"nonzero physical zero-gross row");}
        if(fixed[i]){report.maxReactionCoefficientL1=std::max(report.maxReactionCoefficientL1,rowNorm[i]);
            report.maxReactionConstantRawError=std::max(report.maxReactionConstantRawError,std::abs(difference));}
        else {require(rowNorm[i]<=2*gross[i]+tolerance*gross[i],"physical free coefficient bound failed");
            if(gross[i]>0)report.coefficientRatio=std::max(report.coefficientRatio,rowNorm[i]/(2*gross[i]));}
    }
    return result;
}
std::vector<HydroConvectionCoupling> constrainOpaqueMoments(const HydroDualMaterial& m,
    const std::vector<NativeLeafFace>& faces,const std::vector<std::vector<int>>& cellFaces,
    const std::vector<std::vector<int>>& neighbors,const std::vector<Vector2D>& points,
    const std::vector<std::vector<HydroConvectionWeight>>& weights,HydroConvectionMomentReport& report)
{
    const auto& l=m.layout();const int n=l.velocityCount();const auto& segments=m.geometry().segments();
    long double rateScale=0;for(auto f:m.faceRates())rateScale=std::max(rateScale,std::abs(f));
    if(rateScale==0)return {};
    using Pair=std::pair<int,int>;std::map<Pair,long double> original;
    const auto add=[&](int i,int j,long double value){if(i==j)return;if(i>j){std::swap(i,j);value=-value;}original[{i,j}]+=value;};
    std::vector<double> rowScale(n),h(n);std::vector<bool> exterior(n);
    std::vector<std::array<double,6>> rhs(n);
    for(int i=0;i<n;++i){const auto& box=m.geometry().volumes()[i];h[i]=std::max(box.upper.x-box.lower.x,box.upper.y-box.lower.y);require(h[i]>0,"invalid native moment scale");}
    for(std::size_t g=0;g<segments.size();++g){const auto& s=segments[g];const long double f=m.faceRates()[g]/rateScale;
        for(int i:{s.negativeVelocity,s.positiveVelocity})if(i>=0)rowScale[i]+=static_cast<double>(std::abs(f));
        if(internal(s))for(const auto& w:weights[g]){add(s.negativeVelocity,w.velocity,.5L*f*w.weight);add(s.positiveVelocity,w.velocity,-.5L*f*w.weight);}
        else exterior[owner(s)]=true;
    }
    for(const auto& edge:original){rowScale[edge.first.first]+=static_cast<double>(std::abs(edge.second));rowScale[edge.first.second]+=static_cast<double>(std::abs(edge.second));}
    for(std::size_t g=0;g<segments.size();++g){const auto& s=segments[g];const double f=static_cast<double>(m.faceRates()[g]/rateScale);
        for(int i:{s.negativeVelocity,s.positiveVelocity})if(i>=0 && !exterior[i]){const auto p=polynomial(faces[i].center,points[g],h[i]);const double sign=i==s.negativeVelocity?1:-1;
            rhs[i][1]+=sign*f*p[1];rhs[i][2]+=sign*f*p[2];}
    }
    for(const auto& edge:original){const int i=edge.first.first,j=edge.first.second;const double k=static_cast<double>(edge.second);
        for(int row:{i,j})if(!exterior[row]){const auto p=polynomial(faces[row].center,faces[row==i?j:i].center,h[row]);const double sign=row==i?1:-1;
            for(int c=1;c<6;++c)rhs[row][c]-=sign*k*p[c];}}
    // All pairs in the same geometry-only one-ring patch, shared globally.
    std::map<Pair,std::size_t> ids;
    for(int u=0;u<n;++u){std::vector<int> roots,patch,support;
        for(const auto& p:m.geometry().volumes()[u].overlaps)insert(roots,p.pressure);
        patch=roots;for(int p:roots)for(int q:neighbors[p])insert(patch,q);
        for(int p:patch)for(int v:cellFaces[p])if(component(l,v)==component(l,u))insert(support,v);
        std::sort(support.begin(),support.end());
        for(std::size_t a=0;a<support.size();++a)for(std::size_t b=a+1;b<support.size();++b){const int i=support[a],j=support[b];
            if(rowScale[i]>0 && rowScale[j]>0)ids.emplace(Pair{i,j},0);}
    }
    std::vector<Pair> pairs;std::vector<double> pairScale;
    for(auto& entry:ids){entry.second=pairs.size();pairs.push_back(entry.first);pairScale.push_back(std::sqrt(rowScale[entry.first.first]*rowScale[entry.first.second]));}
    std::vector<MomentBlock> blocks(n);std::vector<double> target;
    for(int i=0;i<n;++i){auto& b=blocks[i];b.offset=static_cast<int>(target.size());b.count=rowScale[i]>0?(exterior[i]?1:6):0;
        if(b.count==6)++report.interiorRows;for(int c=0;c<b.count;++c)target.push_back(rhs[i][c]/rowScale[i]);}
    report.constraints=static_cast<int>(target.size());
    for(std::size_t id=0;id<pairs.size();++id){const auto [i,j]=pairs[id];
        for(int row:{i,j}){MomentEntry e{id,polynomial(faces[row].center,faces[row==i?j:i].center,h[row])};const double scale=(row==i?1:-1)*pairScale[id]/rowScale[row];
            for(double& a:e.a)a*=scale;blocks[row].entries.push_back(e);}}
    const auto originalBlocks=blocks;const auto originalTarget=target;
    for(auto& b:blocks)if(b.count){Eigen::MatrixXd a(b.count,b.entries.size());for(Eigen::Index j=0;j<a.cols();++j)for(int c=0;c<b.count;++c)a(c,j)=b.entries[static_cast<std::size_t>(j)].a[c];
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(a,Eigen::ComputeThinU|Eigen::ComputeThinV);svd.setThreshold(1e-12);
        require(svd.rank()==b.count,"unidentifiable native moment patch");
        const Eigen::MatrixXd whitening=svd.singularValues().cwiseInverse().asDiagonal()*svd.matrixU().transpose();b.whitening=whitening;
        const Eigen::MatrixXd whitened=whitening*a;Eigen::VectorXd rhsBlock(b.count);for(int c=0;c<b.count;++c)rhsBlock[c]=target[b.offset+c];
        const Eigen::VectorXd whitenedTarget=whitening*rhsBlock;
        require(whitened.allFinite() && whitenedTarget.allFinite(),"nonfinite native moment preconditioner");
        for(int c=0;c<b.count;++c){target[b.offset+c]=whitenedTarget[c];for(std::size_t j=0;j<b.entries.size();++j)b.entries[j].a[c]=whitened(c,static_cast<Eigen::Index>(j));}}
    const auto transpose=[&](const std::vector<double>& y){std::vector<double> x(pairs.size()),compensation(x.size());for(const auto& b:blocks)for(const auto& e:b.entries){double value=0,local=0;for(int c=0;c<b.count;++c)accumulate(e.a[c]*y[b.offset+c],value,local);accumulate(value,x[e.edge],compensation[e.edge]);}return x;};
    const auto multiply=[&](const std::vector<double>& x){std::vector<double> y(target.size());for(const auto& b:blocks)for(int c=0;c<b.count;++c){double compensation=0;for(const auto& e:b.entries)accumulate(e.a[c]*x[e.edge],y[b.offset+c],compensation);}return y;};
    const auto precondition=[](const std::vector<double>& r){return r;};
    const auto residual=[&](const std::vector<double>& x){auto r=originalTarget;for(const auto& b:originalBlocks)for(int c=0;c<b.count;++c){double actual=0,compensation=0;for(const auto& e:b.entries)accumulate(e.a[c]*x[e.edge],actual,compensation);r[b.offset+c]-=actual;}return r;};
    const auto residualScale=[&](const std::vector<double>& x){std::vector<double> scale(originalTarget.size(),1);for(const auto& b:originalBlocks)for(int c=0;c<b.count;++c)for(const auto& e:b.entries)scale[b.offset+c]+=std::abs(e.a[c]*x[e.edge]);return scale;};
    const auto whiten=[&](const std::vector<double>& r){std::vector<double> y(r.size());for(const auto& b:blocks)for(int c=0;c<b.count;++c){double compensation=0;for(int k=0;k<b.count;++k)accumulate(b.whitening(c,k)*r[b.offset+k],y[b.offset+c],compensation);}return y;};
    std::vector<double> solution(pairs.size()),solutionCompensation(solution.size()),r=target,z=precondition(r),direction=z;
    long double rz=dot(r,z);const double stopping=8*std::numeric_limits<double>::epsilon()*std::max(1.,maximum(target));
    int iterations=0;const auto maximumIterations=std::max<std::size_t>(10000,4*target.size());
    while(static_cast<std::size_t>(iterations)<maximumIterations){
        if(maximum(r)<=stopping){const auto actual=residual(solution),scale=residualScale(solution);double error=0;for(std::size_t i=0;i<actual.size();++i)error=std::max(error,std::abs(actual[i])/scale[i]);if(error<=8*std::numeric_limits<double>::epsilon())break;
            r=whiten(actual);direction=precondition(r);rz=dot(r,direction);}
        const auto step=transpose(direction),a=multiply(step);const auto denominator=dot(direction,a);
        require(denominator>0 && rz>0,"incompatible native spatial/skew moments");const double alpha=static_cast<double>(rz/denominator);require(std::isfinite(alpha),"nonfinite native moment step");
        // Accumulate the minimum-norm pair solution directly: forming A^T
        // from a large final multiplier loses the small resolved moments.
        for(std::size_t i=0;i<solution.size();++i)accumulate(alpha*step[i],solution[i],solutionCompensation[i]);
        for(std::size_t i=0;i<r.size();++i)r[i]-=alpha*a[i];++iterations;
        // Retain conjugacy across the long-wavelength patch modes. Restart
        // only for explicit residual refinement, never on an arbitrary timer.
        const bool restart=maximum(r)<=stopping;
        if(restart)r=whiten(residual(solution));
        z=precondition(r);const auto next=dot(r,z);const double beta=restart?0:static_cast<double>(next/rz);
        for(std::size_t i=0;i<direction.size();++i)direction[i]=z[i]+beta*direction[i];rz=next;
    }
    report.iterations=iterations;
    const auto scale=residualScale(solution);
    for(const auto& b:originalBlocks)for(int c=0;c<b.count;++c){double actual=0,compensation=0;for(const auto& e:b.entries)accumulate(e.a[c]*solution[e.edge],actual,compensation);const auto error=std::abs(actual-originalTarget[b.offset+c]);
        const auto normalized=error/scale[b.offset+c];
        if(c==0)report.constantError=std::max(report.constantError,static_cast<long double>(normalized));else if(c<3)report.affineError=std::max(report.affineError,static_cast<long double>(normalized));else report.quadraticError=std::max(report.quadraticError,static_cast<long double>(normalized));}
    if(std::max({report.constantError,report.affineError,report.quadraticError})>tolerance)
        throw std::invalid_argument("[HydroKepConvection] native moment solve did not converge: iterations="+std::to_string(iterations)+
            " scaled residual="+std::to_string(std::max({report.constantError,report.affineError,report.quadraticError})/tolerance)+" tolerances");
    std::vector<HydroConvectionCoupling> result;
    for(std::size_t id=0;id<pairs.size();++id)if(solution[id]!=0){const auto value=product(rateScale,static_cast<long double>(pairScale[id])*solution[id]);result.push_back({pairs[id].first,pairs[id].second,value});}
    return result;
}
// A closed final skew operator: exact C/P1, Taylor P2 slabs and bounded
// coefficient norm. ADMM variables are temporary shared-pair algebra;
// no prepared failure can select the opaque operator.
// A shared-pair convex solve. Row normalization is numerical conditioning,
// never a change of physical moment or coefficient budgets.
using Sparse=Eigen::SparseMatrix<double>;
using Triplet=Eigen::Triplet<double>;
struct TypedBoundedProjectionState
{
    Eigen::VectorXd x;
    std::vector<int> endpoint;
    bool converged=false;
    int peakWorkingRows=0,peakFactorCapacity=0;
    std::size_t peakThinFactorBytes=0; // Owned thin-Q/R arrays incl. growth overlap, not total process/SparseQR memory.
};
// Goldfarb-Idnani dual projection of the SAME Euclidean seed objective.
// Working inequalities have nonnegative multipliers: a violated new row
// either reaches its bound or removes an old inequality at dual zero.
// SparseQR applies Q only to vectors; every original dependent hard row and
// every inactive inequality remains in the final physical certificate.
TypedBoundedProjectionState typedBoundedProjectionComponent(const Sparse& raw,const Eigen::VectorXd& lower,
    const Eigen::VectorXd& upper,const Eigen::VectorXd& seed,const Eigen::VectorXd& scale,
    HydroConvectionMomentReport& report,int pivotBudget)
{
    std::vector<int> hardIndex(raw.rows(),-1),hardRows,upperId(raw.rows(),-1),lowerId(raw.rows(),-1),signedRows,signedDirections;
    for(int i=0;i<raw.rows();++i){if(lower[i]==upper[i]){hardIndex[i]=static_cast<int>(hardRows.size());hardRows.push_back(i);}
        else{if(std::isfinite(upper[i])){upperId[i]=static_cast<int>(signedRows.size());signedRows.push_back(i);signedDirections.push_back(1);}
            if(std::isfinite(lower[i])){lowerId[i]=static_cast<int>(signedRows.size());signedRows.push_back(i);signedDirections.push_back(-1);}}}
    require(!hardRows.empty(),"typed bounded problem has no hard moment rows");
    std::vector<Triplet> hardEntries,inequalityEntries;
    for(int j=0;j<raw.outerSize();++j)for(Sparse::InnerIterator e(raw,j);e;++e){const int i=static_cast<int>(e.row());const double coefficient=scale[i]*e.value();
        if(hardIndex[i]>=0)hardEntries.emplace_back(hardIndex[i],j,coefficient);
        else{if(upperId[i]>=0)inequalityEntries.emplace_back(upperId[i],j,coefficient);
            if(lowerId[i]>=0)inequalityEntries.emplace_back(lowerId[i],j,-coefficient);}}
    Sparse hard(hardRows.size(),raw.cols()),inequality(signedRows.size(),raw.cols());
    hard.setFromTriplets(hardEntries.begin(),hardEntries.end());inequality.setFromTriplets(inequalityEntries.begin(),inequalityEntries.end());
    Eigen::VectorXd allHardTarget(hard.rows());for(int k=0;k<hard.rows();++k)allHardTarget[k]=scale[hardRows[k]]*lower[hardRows[k]];
    TypedBoundedProjectionState state;state.endpoint.assign(raw.rows(),0);for(int i:hardRows)state.endpoint[i]=2;
    int hardRank=0,count=0,capacity=0;Eigen::MatrixXd thin,triangular;Eigen::VectorXd baseTarget;std::vector<Triplet> baseEntries;
    const auto factorBytes=[&](int n)->std::size_t{return sizeof(double)*(static_cast<std::size_t>(raw.cols())*n+static_cast<std::size_t>(n)*n);};
    {
        // H^T*P=Q*R. Selecting P's first rank hard rows in exactly that
        // order makes the working-base column order equal to R11's order.
        const Sparse hardTranspose=hard.transpose();Eigen::SparseQR<Sparse,Eigen::COLAMDOrdering<int>> hardQr;hardQr.compute(hardTranspose);
        require(hardQr.info()==Eigen::Success,"typed hard moment QR failed");
        hardRank=static_cast<int>(hardQr.rank());require(hardRank>0,"unidentifiable typed hard moment rows");
        count=capacity=hardRank;baseTarget.resize(hardRank);std::vector<int> independent(hard.rows(),-1),hardOrder(hardRank);
        for(int k=0;k<hardRank;++k){const int original=hardQr.colsPermutation().indices()[k];
            require(original>=0 && original<hard.rows() && independent[original]<0,"invalid typed hard basis permutation");
            independent[original]=k;hardOrder[k]=original;baseTarget[k]=allHardTarget[original];}
        for(int k=0;k<hardRank;++k)require(baseTarget[k]==allHardTarget[hardOrder[k]],"different ordered typed hard targets");
        for(int j=0;j<hard.outerSize();++j)for(Sparse::InnerIterator e(hard,j);e;++e)
            if(independent[e.row()]>=0)baseEntries.emplace_back(independent[e.row()],j,e.value());
        thin.resize(raw.cols(),capacity);Eigen::VectorXd axis=Eigen::VectorXd::Zero(raw.cols());
        // Materialize only working-rank columns, never a full variable-square Q.
        for(int j=0;j<count;++j){if(j)axis[j-1]=0;axis[j]=1;const Eigen::VectorXd column=hardQr.matrixQ()*axis;thin.col(j)=column;}
        triangular=Eigen::MatrixXd(hardQr.matrixR().topLeftCorner(count,count));
        require(thin.allFinite() && triangular.allFinite(),"nonfinite typed thin hard factors");
    }
    state.peakWorkingRows=count;state.peakFactorCapacity=capacity;state.peakThinFactorBytes=factorBytes(capacity);
    const auto reserve=[&](int requested){
        if(requested<=capacity)return;
        const int next=std::min(static_cast<int>(raw.cols()),((requested+63)/64)*64);
        require(next>=requested,"typed thin factor exceeds variable rank");
        state.peakThinFactorBytes=std::max(state.peakThinFactorBytes,factorBytes(capacity)+factorBytes(next));
        Eigen::MatrixXd enlargedThin(raw.cols(),next),enlargedTriangular=Eigen::MatrixXd::Zero(next,next);
        enlargedThin.leftCols(count)=thin.leftCols(count);enlargedTriangular.topLeftCorner(count,count)=triangular.topLeftCorner(count,count);
        thin.swap(enlargedThin);triangular.swap(enlargedTriangular);capacity=next;state.peakFactorCapacity=std::max(state.peakFactorCapacity,capacity);
    };
    Eigen::VectorXd bounds(signedRows.size());for(int k=0;k<bounds.size();++k){const int i=signedRows[k];
        bounds[k]=scale[i]*(signedDirections[k]>0?upper[i]:-lower[i]);}
    std::vector<int> active;std::vector<bool> included(signedRows.size(),false);
    Sparse base(hardRank,raw.cols());base.setFromTriplets(baseEntries.begin(),baseEntries.end());
    const Eigen::VectorXd initialResidual=baseTarget-base*seed;
    const Eigen::VectorXd initialCorrection=triangular.topLeftCorner(count,count).transpose().triangularView<Eigen::Lower>().solve(initialResidual);
    state.x=seed+thin.leftCols(count)*initialCorrection;
    Eigen::VectorXd multipliers=-triangular.topLeftCorner(count,count).triangularView<Eigen::Upper>().solve(initialCorrection);
    require(state.x.allFinite() && multipliers.allFinite(),"nonfinite typed GI seed projection");
    const Sparse rowVectors=inequality.transpose();int pivots=0;
    while(pivots<pivotBudget){
        const Eigen::VectorXd violation=inequality*state.x-bounds;int next=-1;double largest=tolerance;
        for(int i=0;i<violation.size();++i)if(!included[i] && violation[i]>largest){largest=violation[i];next=i;}
        if(next<0){state.converged=true;break;}
        double newMultiplier=0;
        while(pivots<pivotBudget){
            ++pivots;++report.iterations;Eigen::VectorXd row=Eigen::VectorXd::Zero(raw.cols());
            for(Sparse::InnerIterator e(rowVectors,next);e;++e)row[e.row()]=e.value();
            // Two-pass reorthogonalization preserves both the new-column
            // reconstruction a=Q*h+p and the Euclidean nullspace direction.
            Eigen::VectorXd coordinates=thin.leftCols(count).transpose()*row;
            Eigen::VectorXd primalDirection=row-thin.leftCols(count)*coordinates;
            const Eigen::VectorXd correction=thin.leftCols(count).transpose()*primalDirection;
            coordinates+=correction;primalDirection-=thin.leftCols(count)*correction;
            const Eigen::VectorXd dualDirection=triangular.topLeftCorner(count,count).triangularView<Eigen::Upper>().solve(coordinates);
            require(dualDirection.allFinite() && primalDirection.allFinite(),"nonfinite typed GI incremental direction");
            const double tangentSquared=primalDirection.squaredNorm();
            double dualStep=std::numeric_limits<double>::infinity();int remove=-1;
            for(int i=hardRank;i<multipliers.size();++i)if(dualDirection[i]>0){const double step=std::max(0.,multipliers[i])/dualDirection[i];
                if(step<dualStep){dualStep=step;remove=i-hardRank;}}
            const double roundoff=std::numeric_limits<double>::epsilon()*std::numeric_limits<double>::epsilon()*row.squaredNorm();
            const double primalStep=tangentSquared>roundoff?std::max(0.,row.dot(state.x)-bounds[next])/tangentSquared:std::numeric_limits<double>::infinity();
            const double step=std::min(primalStep,dualStep);
            require(std::isfinite(step) && step>=0,"incompatible typed GI inequality without dual exit");
            state.x-=step*primalDirection;multipliers-=step*dualDirection;newMultiplier+=step;
            require(state.x.allFinite() && multipliers.allFinite() && std::isfinite(newMultiplier),"nonfinite typed GI pivot");
            if(primalStep<=dualStep){
                require(count<raw.cols(),"dependent typed GI append");const double norm=std::sqrt(tangentSquared);
                require(norm>0 && std::isfinite(norm),"nonfinite typed GI append norm");reserve(count+1);
                thin.col(count)=primalDirection/norm;triangular.block(0,count,count,1)=coordinates;
                triangular.block(count,0,1,count).setZero();triangular(count,count)=norm;++count;
                state.peakWorkingRows=std::max(state.peakWorkingRows,count);included[next]=true;active.push_back(next);
                multipliers.conservativeResize(multipliers.size()+1);multipliers[multipliers.size()-1]=newMultiplier;break;
            }
            require(remove>=0,"missing typed GI dual-zero inequality");const int erased=hardRank+remove;const int previousCount=count;
            // Delete the working constraint COLUMN, then chase its upper
            // Hessenberg bulge. makeGivens gives G^T*[p,q]=[r,0]:
            // R_new=G^T*R and Q_new=Q*G preserve the remaining W^T=Q*R.
            if(erased+1<previousCount)triangular.block(0,erased,previousCount,previousCount-erased-1)
                =triangular.block(0,erased+1,previousCount,previousCount-erased-1).eval();
            for(int i=erased;i<previousCount-1;++i){Eigen::JacobiRotation<double> rotation;rotation.makeGivens(triangular(i,i),triangular(i+1,i));
                triangular.topLeftCorner(previousCount,previousCount-1).applyOnTheLeft(i,i+1,rotation.adjoint());
                thin.leftCols(previousCount).applyOnTheRight(i,i+1,rotation);triangular(i+1,i)=0;}
            --count;included[active[remove]]=false;active.erase(active.begin()+remove);
            Eigen::VectorXd retained(multipliers.size()-1);retained.head(erased)=multipliers.head(erased);
            retained.tail(retained.size()-erased)=multipliers.tail(retained.size()-erased);multipliers=std::move(retained);
            // Keep x and the pending new multiplier: this is a factor update,
            // never an intermediate reprojection that changes the GI step.
        }
    }
    if(!state.converged){
        const Eigen::VectorXd violation=inequality*state.x-bounds;state.converged=true;
        for(int i=0;i<violation.size();++i)if(!included[i] && violation[i]>tolerance){state.converged=false;break;}
    }
    require(state.converged,"typed bounded GI reached the shared 20000 pivot work limit");
    require(count==hardRank+static_cast<int>(active.size()),"different typed incremental working dimensions");
    // Anchor the terminal SAME working set with one fresh original SparseQR.
    // This removes factor-update drift without adding any pivot/refinement
    // budget or changing seed, constraints, pending GI state, or operator.
    std::vector<Triplet> finalEntries=baseEntries;std::vector<int> position(signedRows.size(),-1);
    Eigen::VectorXd target(count);target.head(hardRank)=baseTarget;
    for(std::size_t k=0;k<active.size();++k){position[active[k]]=static_cast<int>(k);target[hardRank+k]=bounds[active[k]];}
    for(int j=0;j<inequality.outerSize();++j)for(Sparse::InnerIterator e(inequality,j);e;++e)
        if(position[e.row()]>=0)finalEntries.emplace_back(hardRank+position[e.row()],j,e.value());
    Sparse working(count,raw.cols());working.setFromTriplets(finalEntries.begin(),finalEntries.end());
    const Sparse transpose=working.transpose();Eigen::SparseQR<Sparse,Eigen::COLAMDOrdering<int>> finalQr;finalQr.compute(transpose);
    require(finalQr.info()==Eigen::Success && finalQr.rank()==working.rows(),"dependent typed GI terminal working constraints");
    Sparse leading=finalQr.matrixR().topLeftCorner(count,count);leading.makeCompressed();
    const Eigen::VectorXd residual=target-working*seed,permuted=finalQr.colsPermutation().transpose()*residual;
    const Eigen::VectorXd terminalCorrection=leading.transpose().triangularView<Eigen::Lower>().solve(permuted);
    require(terminalCorrection.allFinite(),"nonfinite typed GI terminal seed projection");
    Eigen::VectorXd terminalCoordinates=Eigen::VectorXd::Zero(raw.cols());terminalCoordinates.head(count)=terminalCorrection;
    state.x=seed+finalQr.matrixQ()*terminalCoordinates;
    const Eigen::VectorXd dual=-leading.triangularView<Eigen::Upper>().solve(terminalCorrection);multipliers=finalQr.colsPermutation()*dual;
    require(state.x.allFinite() && multipliers.allFinite(),"nonfinite typed GI terminal projection");
    const double dualScale=std::max(1.,multipliers.cwiseAbs().maxCoeff());
    for(int i=hardRank;i<multipliers.size();++i)require(multipliers[i]>=-tolerance*dualScale,"typed GI active dual sign certificate failed");
    for(int k:active)state.endpoint[signedRows[k]]=signedDirections[k];
    return state;
}
// The frozen pair-coefficient QP is a product of its raw variable-row
// connected components: each original row touches exactly one component,
// and the Euclidean seed objective has no cross-variable Hessian terms.
// Only GI preparation is decomposed; the caller retains ONE global
// 32-pass refinement and every original physical/facet/free-L1 gate.
TypedBoundedProjectionState typedBoundedProjection(const Sparse& raw,const Eigen::VectorXd& lower,
    const Eigen::VectorXd& upper,const Eigen::VectorXd& seed,const Eigen::VectorXd& scale,
    HydroConvectionMomentReport& report)
{
    std::vector<int> parent(raw.cols()),first(raw.rows(),-1);std::iota(parent.begin(),parent.end(),0);
    const auto root=[&](int i){while(parent[i]!=i){parent[i]=parent[parent[i]];i=parent[i];}return i;};
    for(int j=0;j<raw.outerSize();++j)for(Sparse::InnerIterator e(raw,j);e;++e)if(e.value()!=0){
        if(first[e.row()]<0)first[e.row()]=j;
        else{const int a=root(first[e.row()]),b=root(j);if(a!=b)parent[b]=a;}}
    std::map<int,int> ids;std::vector<std::vector<int>> columns,rows;
    std::vector<int> columnComponent(raw.cols()),localColumn(raw.cols()),rowComponent(raw.rows()),localRow(raw.rows());
    for(int j=0;j<raw.cols();++j){const int representative=root(j);const auto inserted=ids.emplace(representative,static_cast<int>(columns.size()));
        if(inserted.second){columns.emplace_back();rows.emplace_back();}
        const int component=inserted.first->second;columnComponent[j]=component;localColumn[j]=static_cast<int>(columns[component].size());columns[component].push_back(j);}
    for(int i=0;i<raw.rows();++i){require(first[i]>=0,"zero typed component constraint row");
        const int component=columnComponent[first[i]];rowComponent[i]=component;localRow[i]=static_cast<int>(rows[component].size());rows[component].push_back(i);}
    if(columns.size()==1)return typedBoundedProjectionComponent(raw,lower,upper,seed,scale,report,20000);
    std::vector<std::vector<Triplet>> entries(columns.size());
    for(int j=0;j<raw.outerSize();++j)for(Sparse::InnerIterator e(raw,j);e;++e)if(e.value()!=0){const int component=columnComponent[j];
        require(rowComponent[e.row()]==component,"typed constraint crosses independent components");
        entries[component].emplace_back(localRow[e.row()],localColumn[j],e.value());}
    TypedBoundedProjectionState result;result.x=seed;result.endpoint.assign(raw.rows(),0);result.converged=true;
    // This budget is local to this QP, independent of report.iterations
    // accumulated by previous facet solves or other prepared operators.
    int spent=0;
    for(std::size_t c=0;c<columns.size();++c){if(rows[c].empty())continue; // Unconstrained variables retain their unique seed minimum.
        Sparse local(rows[c].size(),columns[c].size());local.setFromTriplets(entries[c].begin(),entries[c].end());
        Eigen::VectorXd lo(rows[c].size()),hi(rows[c].size()),conditioning(rows[c].size()),origin(columns[c].size());
        for(std::size_t i=0;i<rows[c].size();++i){lo[i]=lower[rows[c][i]];hi[i]=upper[rows[c][i]];conditioning[i]=scale[rows[c][i]];}
        for(std::size_t j=0;j<columns[c].size();++j)origin[j]=seed[columns[c][j]];
        const int before=report.iterations;const auto state=typedBoundedProjectionComponent(local,lo,hi,origin,conditioning,report,20000-spent);
        const int used=report.iterations-before;require(used>=0 && used<=20000-spent,"typed component exceeded the shared GI pivot budget");spent+=used;
        require(state.converged,"typed independent component did not converge");
        result.peakWorkingRows=std::max(result.peakWorkingRows,state.peakWorkingRows);
        result.peakFactorCapacity=std::max(result.peakFactorCapacity,state.peakFactorCapacity);
        result.peakThinFactorBytes=std::max(result.peakThinFactorBytes,state.peakThinFactorBytes);
        for(std::size_t j=0;j<columns[c].size();++j)result.x[columns[c][j]]=state.x[j];
        for(std::size_t i=0;i<rows[c].size();++i)result.endpoint[rows[c][i]]=state.endpoint[i];
    }
    require(result.x.allFinite(),"nonfinite assembled typed component coefficients");return result;
}
void certifyFinalBoundedRows(const Sparse& raw,const Eigen::VectorXd& lower,
    const Eigen::VectorXd& upper,const Eigen::VectorXd& x,const std::vector<int>* reactionConstantRows)
{
    // Certify the complete moment/slab/facet QP after the full L1 bounds
    // stabilize. Intermediate facet solves do not publish an operator.
    require(x.allFinite(),"nonfinite final bounded coefficients");
    const Eigen::VectorXd original=raw*x;
    for(int i=0;i<raw.rows();++i){
        if(reactionConstantRows && std::binary_search(reactionConstantRows->begin(),reactionConstantRows->end(),i)){
            require(lower[i]==upper[i],"reaction constant row is not an equality");
            long double arithmetic=0;for(int j=0;j<raw.outerSize();++j)arithmetic+=std::abs(raw.coeff(i,j)*x[j]);
            const auto error=finite(std::abs(original[i]-lower[i]),"nonfinite reaction constant error");
            const auto physicalScale=finite(arithmetic+2*std::abs(lower[i]),"nonfinite reaction constant arithmetic");
            require(physicalScale>0 || error==0,"nonzero reaction constant error at zero physical scale");
            require(error<=tolerance*physicalScale,"physical reaction constant moment certificate failed");
        }else if(!(original[i]>=lower[i]-tolerance && original[i]<=upper[i]+tolerance)){
            long double arithmetic=0;
            for(int j=0;j<raw.outerSize();++j)arithmetic+=std::abs(static_cast<long double>(raw.coeff(i,j))*x[j]);
            throw std::invalid_argument("[HydroKepConvection] physical free moment/slab certificate failed after refinement: row="+std::to_string(i)+
                " kind="+(lower[i]==upper[i]?std::string("hard"):std::string("slab"))+
                " lower_gap/tolerance="+std::to_string((lower[i]-original[i])/tolerance)+
                " upper_gap/tolerance="+std::to_string((original[i]-upper[i])/tolerance)+
                " raw_arithmetic="+std::to_string(arithmetic));}
    }
}
Eigen::VectorXd boundedProjection(const Sparse& raw,const Eigen::VectorXd& lower,
    const Eigen::VectorXd& upper,const Eigen::VectorXd& seed,HydroConvectionMomentReport& report,bool certifyPhysical=false,const std::vector<int>* reactionConstantRows=nullptr)
{
    Eigen::VectorXd scale=Eigen::VectorXd::Zero(raw.rows());
    for(int j=0;j<raw.outerSize();++j)for(Sparse::InnerIterator e(raw,j);e;++e)scale[e.row()]+=e.value()*e.value();
    for(int i=0;i<scale.size();++i){require(scale[i]>0,"unidentifiable bounded moment row");scale[i]=std::sqrt((lower[i]==upper[i]?1000.:1.)/scale[i]);}
    const Sparse a=scale.asDiagonal()*raw;const Eigen::VectorXd lo=scale.cwiseProduct(lower),hi=scale.cwiseProduct(upper);
    const auto clip=[&](const Eigen::VectorXd& x)->Eigen::VectorXd{return x.cwiseMax(lo).cwiseMin(hi);};
    Eigen::VectorXd x=seed,z,y;bool converged=false;std::vector<int> endpoint;
    if(certifyPhysical){auto state=typedBoundedProjection(raw,lower,upper,seed,scale,report);x=std::move(state.x);endpoint=std::move(state.endpoint);converged=state.converged;}
    else{
    const Sparse transpose=a.transpose();const double rho=10;
    Sparse normal=a*transpose;for(int i=0;i<normal.rows();++i)normal.coeffRef(i,i)+=1/rho;normal.makeCompressed();
    Eigen::SimplicialLDLT<Sparse> factor;factor.compute(normal);require(factor.info()==Eigen::Success,"bounded moment factorization failed");
    const Eigen::VectorXd origin=a*seed;z=clip(origin);y=Eigen::VectorXd::Zero(a.rows());

    for(int iteration=0;iteration<20000;++iteration){
        const Eigen::VectorXd multiplier=factor.solve(z-y-origin);require(factor.info()==Eigen::Success && multiplier.allFinite(),"nonfinite bounded moment solve");
        x=seed+transpose*multiplier;const Eigen::VectorXd value=a*x,previous=z,relaxed=1.6*value-.6*previous;
        z=clip(relaxed+y);y+=relaxed-z;++report.iterations;
        if(iteration%25==0){const double primal=(value-z).lpNorm<Eigen::Infinity>(),dual=(rho*transpose*(z-previous)).lpNorm<Eigen::Infinity>();
            if(std::max(primal,dual)<=1e-10){converged=true;break;}}
    }
    }
    if(!converged && !certifyPhysical){
        const Sparse transpose=a.transpose();const double rho=10;
        const Eigen::VectorXd value=a*x,rawValue=raw*x;
        double violation=0;
        for(int i=0;i<raw.rows();++i)violation=std::max(violation,
            std::max(lower[i]-rawValue[i],rawValue[i]-upper[i]));
        throw std::invalid_argument("[HydroKepConvection] bounded moment ADMM did not converge: primal/stopping="
            +std::to_string((value-z).lpNorm<Eigen::Infinity>()/1e-10)
            +" rawViolation/tolerance="+std::to_string(violation/tolerance));
    }
    // Refine the identified active constraints jointly. Typed GI uses
    // implicit sparse QR minimum-norm corrections. The historical periodic
    // branch retains its regularized Schur factor and residual refinement.
    const Eigen::VectorXd value=a*x;std::vector<int> active(a.rows(),-1);std::vector<double> targets;
    for(int i=0;i<a.rows();++i){double target=0;bool use=false;
        if(lo[i]==hi[i]){use=true;target=lo[i];}
        else if(certifyPhysical && endpoint[i]!=0){use=true;target=endpoint[i]<0?lo[i]:hi[i];}
        else if(!certifyPhysical && std::isfinite(lo[i]) && std::abs(value[i]-lo[i])<1e-8){use=true;target=lo[i];}
        else if(!certifyPhysical && std::isfinite(hi[i]) && std::abs(value[i]-hi[i])<1e-8){use=true;target=hi[i];}
        if(use){active[i]=static_cast<int>(targets.size());targets.push_back(target);}}
    std::vector<Triplet> entries;for(int j=0;j<a.outerSize();++j)for(Sparse::InnerIterator e(a,j);e;++e)if(active[e.row()]>=0)entries.emplace_back(active[e.row()],j,e.value());
    Sparse hard(targets.size(),a.cols());hard.setFromTriplets(entries.begin(),entries.end());const Sparse ht=hard.transpose();
    const Eigen::Map<const Eigen::VectorXd> target(targets.data(),targets.size());
    if(certifyPhysical){
        std::vector<Triplet> rawEntries;
        for(int j=0;j<raw.outerSize();++j)for(Sparse::InnerIterator e(raw,j);e;++e)
            if(active[e.row()]>=0)rawEntries.emplace_back(active[e.row()],j,e.value());
        Sparse rawActive(targets.size(),raw.cols());rawActive.setFromTriplets(rawEntries.begin(),rawEntries.end());
        Eigen::VectorXd physicalTarget(targets.size()),activeScale(targets.size());std::vector<bool> reactionActive(targets.size(),false);
        for(int i=0;i<raw.rows();++i)if(active[i]>=0){const int row=active[i];
            physicalTarget[row]=(lo[i]==hi[i] || endpoint[i]<0)?lower[i]:upper[i];activeScale[row]=scale[i];
            reactionActive[row]=reactionConstantRows && std::binary_search(reactionConstantRows->begin(),reactionConstantRows->end(),i);}
        Eigen::SparseQR<Sparse,Eigen::COLAMDOrdering<int>> polish;polish.compute(ht);
        require(polish.info()==Eigen::Success && polish.rank()>0,"typed active moment QR failed");
        const int rank=static_cast<int>(polish.rank());Sparse leading=polish.matrixR().topLeftCorner(rank,rank);leading.makeCompressed();
        // Keep the same scaled factor/rowspace, but resolve its original
        // physical residual. Reaction rows retain their arithmetic-relative scale.
        for(int pass=0;pass<32;++pass){const Eigen::VectorXd rawResidual=physicalTarget-rawActive*x;bool resolved=true;
            for(int row=0;row<rawActive.rows();++row){const auto error=std::abs(rawResidual[row]);
                if(reactionActive[row]){long double arithmetic=0;
                    for(int j=0;j<rawActive.outerSize();++j)arithmetic+=std::abs(rawActive.coeff(row,j)*x[j]);
                    const auto physicalScale=finite(arithmetic+2*std::abs(physicalTarget[row]),"nonfinite active reaction arithmetic");
                    require(physicalScale>0 || error==0,"nonzero active reaction residual at zero physical scale");
                    if(error>8*std::numeric_limits<double>::epsilon()*physicalScale)resolved=false;
                }else if(error>8*std::numeric_limits<double>::epsilon())resolved=false;}
            if(resolved)break;
            const Eigen::VectorXd residual=activeScale.cwiseProduct(rawResidual);
            const Eigen::VectorXd permuted=polish.colsPermutation().transpose()*residual;
            const Eigen::VectorXd correction=leading.transpose().triangularView<Eigen::Lower>().solve(permuted.head(rank));
            require(correction.allFinite(),"nonfinite typed active moment refinement");
            Eigen::VectorXd coordinates=Eigen::VectorXd::Zero(raw.cols());coordinates.head(rank)=correction;
            x+=polish.matrixQ()*coordinates;++report.projectionCalls;}
    }else{
        Sparse schur=hard*ht;for(int i=0;i<schur.rows();++i)schur.coeffRef(i,i)+=1e-9;schur.makeCompressed();
        Eigen::SimplicialLDLT<Sparse> polish;polish.compute(schur);require(polish.info()==Eigen::Success,"bounded active moment factorization failed");
        for(int pass=0;pass<32;++pass){const Eigen::VectorXd residual=target-hard*x;
            if(residual.lpNorm<Eigen::Infinity>()<=8*std::numeric_limits<double>::epsilon())break;
            const Eigen::VectorXd dual=polish.solve(residual);require(dual.allFinite(),"nonfinite bounded active refinement");x+=ht*dual;++report.projectionCalls;}
    }
    require(x.allFinite(),"nonfinite final bounded coefficients");

    return x;
}
struct NativeFaceResponse
{
    std::array<long double,2> mean{};
    std::array<long double,3> quadratic{}; // xx, xy, yy advective responses about this face center.
};
// Native U observes a true face mean. The typed P1 and P2 targets consume
// that same functional of the retained actual source, including tangential
// values. Positive full-degree quadrature splits every original support knot.
NativeFaceResponse nativeFaceMean(const HydroSpatialAdvector& actual,int velocity,
    HydroConvectionMomentReport& report)
{
    const auto& regular=actual.reconstruction().regularGeometry();const auto& face=regular.faces()[velocity];
    require(!actual.momentumGeometry().layout().domain().axisymmetric() && face.area>0,
        "native mean P1 requires a positive XY face measure");
    const int normal=face.normalAxis,tangent=1-normal;
    const int pressure=face.negativePressure>=0?face.negativePressure:face.positivePressure;
    require(pressure>=0,"native mean P1 lacks an actual source owner");
    const auto coordinate=[&](Vector2D p){return tangent==0?p.x:p.y;};
    const double begin=coordinate(face.lower),end=coordinate(face.upper);
    const long double length=static_cast<long double>(end)-begin;
    require(length>0 && std::isfinite(length),"invalid native mean P1 face length");
    std::vector<double> cuts{begin,end};
    for(double knot:regular.cells()[pressure].supportKnots[tangent])if(begin<knot && knot<end)cuts.push_back(knot);
    std::sort(cuts.begin(),cuts.end());cuts.erase(std::unique(cuts.begin(),cuts.end()),cuts.end());
    constexpr double nodes[5]={.0469100770306680036,.230765344947158454,.5,.769234655052841546,.953089922969331996};
    constexpr double weights[5]={.118463442528094544,.239314335249683234,.284444444444444444,.239314335249683234,.118463442528094544};
    NativeFaceResponse response;auto& mean=response.mean;auto& quadratic=response.quadratic;std::array<long double,2> gross{};
    for(std::size_t part=1;part<cuts.size();++part){const double lower=cuts[part-1],upper=cuts[part];
        const long double fraction=(static_cast<long double>(upper)-lower)/length;
        require(fraction>0 && std::isfinite(fraction),"invalid native mean P1 support partition");
        for(int q=0;q<5;++q){auto point=face.lower;(tangent==0?point.x:point.y)=lower+(upper-lower)*nodes[q];
            const auto value=actual.cellVelocity(pressure,point);const long double weight=weights[q]*fraction;
            require(weight>0 && std::isfinite(weight),"invalid native mean P1 quadrature weight");
            const double components[2]={value.x,value.y};
            for(int c=0;c<2;++c){mean[c]+=product(weight,components[c]);gross[c]+=product(weight,std::abs(components[c]));}
            const auto center=(face.lower+face.upper)*.5;
            const long double dx=static_cast<long double>(point.x)-center.x,dy=static_cast<long double>(point.y)-center.y;
            quadratic[0]+=product(weight,2*product(value.x,dx));
            quadratic[1]+=product(weight,product(value.x,dy)+product(value.y,dx));
            quadratic[2]+=product(weight,2*product(value.y,dy));}}
    for(auto value:mean)finite(value,"nonfinite native mean P1 value");
    for(auto value:quadratic)finite(value,"nonfinite native mean P2 response");
    const auto native=static_cast<long double>(actual.values().normalVelocity()[velocity]);
    const auto defect=identity(mean[normal]-native,gross[normal]+std::abs(native),
        "actual native face mean differs from the original normal DOF");
    report.maxNativeMeanDefect=std::max(report.maxNativeMeanDefect,defect);
    const auto point=actual.nativeVelocity(velocity);
    report.maxMeanPointDifference=std::max({report.maxMeanPointDifference,std::abs(mean[0]-point.x),std::abs(mean[1]-point.y)});
    return response;
}
std::vector<HydroConvectionCoupling> constrainBoundedMoments(const HydroDualMaterial& m,
    const std::vector<NativeLeafFace>& faces,const std::vector<std::vector<int>>& cellFaces,
    const std::vector<std::vector<int>>& neighbors,const std::vector<Vector2D>&,
    const std::vector<std::vector<HydroConvectionWeight>>& weights,const NativeTopology& topology,HydroConvectionMomentReport& report,
    const std::vector<bool>* fixed=nullptr,std::vector<HydroConvectionCoupling>* finalPairs=nullptr)
{
    const int n=m.layout().velocityCount();const auto& segments=m.geometry().segments();const auto& actual=m.material().midpointRate().spatialAdvector();
    if(fixed){require(fixed->size()==static_cast<std::size_t>(n) && finalPairs,"different impermeable fixed-row dimensions");
        report.fixedRows=static_cast<int>(std::count(fixed->begin(),fixed->end(),true));report.freeRows=n-report.fixedRows;
        report.boundedRows=report.freeRows;report.interiorRows=report.freeRows;}
    std::vector<NativeFaceResponse> means;
    if(fixed){means.resize(n);for(int i=0;i<n;++i)if(!(*fixed)[i]){means[i]=nativeFaceMean(actual,i,report);++report.meanP1Rows;++report.meanP2Rows;}}
    long double rateScale=0;for(auto f:m.faceRates())rateScale=std::max(rateScale,std::abs(f));
    if(rateScale==0){if(fixed){report.zeroGrossRows=n;
            for(int i=0;i<n;++i){require(m.cells()[i].outwardIntegral==0,"zero material rates have a nonzero authoritative outflow integral");
                if(!(*fixed)[i]){const auto& b=means[i].mean;
                    require(b[0]==0 && b[1]==0,"zero material flux cannot support a nonzero free native-mean moment");
                    for(auto q:means[i].quadratic)require(q==0,"zero material flux cannot support a nonzero free native-quadratic response");}}}
        return {};}
    using Pair=std::pair<int,int>;std::map<Pair,long double> original,central;
    const auto add=[&](std::map<Pair,long double>& map,int i,int j,long double k){i=topology.representative[i];j=topology.representative[j];if(i==j)return;if(i>j){std::swap(i,j);k=-k;}map[{i,j}]+=k;};
    std::vector<double> h(n),gross(n);std::vector<std::array<double,6>> target(n);std::vector<int> row(n,-1),offset(n,-1),width(n,6);
    Vector2D lower=faces.front().center,upper=lower;int count=0;
    for(int i=0;i<n;++i){const int r=topology.representative[i];const auto& box=m.geometry().volumes()[i];h[r]=std::max(h[r],std::max(box.upper.x-box.lower.x,box.upper.y-box.lower.y));
        lower.x=std::min(lower.x,faces[i].center.x);lower.y=std::min(lower.y,faces[i].center.y);upper.x=std::max(upper.x,faces[i].center.x);upper.y=std::max(upper.y,faces[i].center.y);
        if(!fixed || !(*fixed)[i]){const auto b=actual.nativeVelocity(i);const auto mass=m.midpointMass().masses()[i]/rateScale;
            if(fixed){target[r][1]+=static_cast<double>(product(mass,means[i].mean[0]));target[r][2]+=static_cast<double>(product(mass,means[i].mean[1]));
                for(int c=0;c<3;++c)target[r][3+c]+=static_cast<double>(product(mass,means[i].quadratic[c]));}
            else{target[r][1]+=static_cast<double>(mass*b.x);target[r][2]+=static_cast<double>(mass*b.y);}}}
    const double length=std::max({upper.x-lower.x,upper.y-lower.y,topology.period.x,topology.period.y});require(length>0,"invalid bounded domain scale");
    for(std::size_t g=0;g<segments.size();++g){const auto& s=segments[g];const long double f=m.faceRates()[g]/rateScale;
        for(int i:{s.negativeVelocity,s.positiveVelocity})if(i>=0){const int r=topology.representative[i];gross[r]+=static_cast<double>(std::abs(f));target[r][0]+=.5*(i==s.negativeVelocity?1:-1)*static_cast<double>(f);}
        if(internal(s)){add(central,s.negativeVelocity,s.positiveVelocity,.5L*f);for(const auto& w:weights[g]){add(original,s.negativeVelocity,w.velocity,.5L*f*w.weight);add(original,s.positiveVelocity,w.velocity,-.5L*f*w.weight);}}
        else if(static_cast<int>(g)<topology.partner[g]){const auto k=.5L*(s.negativeVelocity>=0?1:-1)*f;const int j=owner(segments[topology.partner[g]]);add(original,owner(s),j,k);add(central,owner(s),j,k);}}
    // The finite-M/uR identity uses the authoritative integrated material
    // defect. Its outward part must follow the same sum-before-divide path;
    // summing already-divided face rates is retained only as a diagnostic.
    if(fixed)for(int i=0;i<n;++i){const auto outward=finite(m.cells()[i].outwardIntegral/m.deltaTime(),"nonfinite canonical impermeable outflow");
        target[i][0]=static_cast<double>((.5L*outward)/rateScale);}
    std::map<Pair,int> ids;
    std::set<Pair> allowedPairs;
    for(int u=0;u<n;++u)if(topology.representative[u]==u){const auto& roots=topology.roots[u];auto patch=roots;std::vector<int> support;
        for(int p:roots)for(int q:neighbors[p])insert(patch,q);for(int p:patch)for(int v:cellFaces[p])if(component(m.layout(),v)==component(m.layout(),u))insert(support,topology.representative[v]);
        std::sort(support.begin(),support.end());for(std::size_t a=0;a<support.size();++a)for(std::size_t b=a+1;b<support.size();++b){const Pair pair{support[a],support[b]};
            if(fixed)allowedPairs.insert(pair);else ids.emplace(pair,0);}}
    std::vector<double> numericScale=gross;
    if(fixed){
        std::vector<int> parent(n);std::iota(parent.begin(),parent.end(),0);
        const auto root=[&](int i){while(parent[i]!=i)i=parent[i];return i;};
        const auto admissible=[&](int i){return (*fixed)[i] || gross[i]>0;};
        for(const auto& e:allowedPairs)if(admissible(e.first)&&admissible(e.second))parent[root(e.second)]=root(e.first);
        std::vector<double> active(n);
        for(int i=0;i<n;++i)if(admissible(i))active[root(i)]=std::max(active[root(i)],gross[i]);
        for(int i=0;i<n;++i)if((*fixed)[i] && gross[i]==0){numericScale[i]=active[root(i)];if(numericScale[i]>0)++report.zeroGrossReactionRows;}
        for(const auto& e:allowedPairs)if(numericScale[e.first]>0 && numericScale[e.second]>0)ids.emplace(e,0);
    }
    int equationCount=0;
    for(int i=0;i<n;++i)if(topology.representative[i]==i){require(h[i]>0,"zero-scale bounded native row");
        if(fixed && gross[i]==0){++report.zeroGrossRows;
            require(m.cells()[i].outwardIntegral==0,"zero-gross row has a nonzero authoritative outflow integral");
            if(!(*fixed)[i]){const auto& b=means[i].mean;
                require(b[0]==0 && b[1]==0,"zero-gross free row cannot support a nonzero native-mean moment");
                for(auto q:means[i].quadratic)require(q==0,"zero-gross free row cannot support a nonzero native-quadratic response");}}
        if(fixed && numericScale[i]==0)continue;
        require(numericScale[i]>0,"zero-scale bounded native row");row[i]=count++;offset[i]=equationCount;
        if(fixed && (*fixed)[i])width[i]=1;equationCount+=width[i];target[i][1]/=h[i];target[i][2]/=h[i];
        if(fixed && !(*fixed)[i])for(int c=3;c<6;++c)target[i][c]/=h[i]*h[i];}
    std::vector<Pair> pairs;std::vector<double> pairScale;std::vector<std::vector<std::pair<int,double>>> incident(n);std::vector<Triplet> entries;
    for(auto& e:ids){const int k=static_cast<int>(pairs.size());e.second=k;pairs.push_back(e.first);const auto [i,j]=e.first;const double scale=std::sqrt(numericScale[i]*numericScale[j]);pairScale.push_back(scale);
        for(int r:{i,j}){const auto p=polynomial(faces[r].center,faces[r==i?j:i].center,h[r],topology);const double sign=r==i?1:-1;incident[r].push_back({k,scale/numericScale[r]});
            report.maximumSupportRatio=std::max(report.maximumSupportRatio,static_cast<long double>(std::hypot(p[1],p[2])));
            for(int c=0;c<width[r];++c)if(p[c]!=0)entries.emplace_back(offset[r]+c,k,sign*scale/numericScale[r]*p[c]);}}
    for(const auto& e:original)require(e.second==0 || (fixed?allowedPairs.count(e.first)==1:ids.count(e.first)==1),"original pair outside declared one-ring");
    for(int i=0;i<n;++i)if(row[i]>=0){const int localCount=width[i]==1?1:3;Eigen::MatrixXd local(localCount,incident[i].size());
        for(std::size_t j=0;j<incident[i].size();++j){const auto& edge=pairs[incident[i][j].first];const int other=edge.first==i?edge.second:edge.first;const auto p=polynomial(faces[i].center,faces[other].center,h[i],topology);for(int c=0;c<localCount;++c)local(c,j)=p[c];}
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(local);svd.setThreshold(1e-12);require(svd.rank()==localCount,"unidentifiable bounded affine patch");}
    Eigen::VectorXd seed(pairs.size());for(int j=0;j<seed.size();++j)seed[j]=static_cast<double>(central[pairs[j]])/pairScale[j];
    std::vector<double> lo,hi;for(int i=0;i<n;++i)if(row[i]>=0)for(int c=0;c<width[i];++c){
        if(fixed && c>=3){const double center=target[i][c]/numericScale[i],halfwidth=h[i]/length;
            lo.push_back(center-halfwidth);hi.push_back(center+halfwidth);
            report.maxQuadraticTargetRatio=std::max(report.maxQuadraticTargetRatio,static_cast<long double>(std::abs(center)/halfwidth));}
        else{lo.push_back(c<3?target[i][c]/numericScale[i]:-h[i]/length);hi.push_back(c<3?target[i][c]/numericScale[i]:h[i]/length);}}
    Sparse matrix(lo.size(),pairs.size());matrix.setFromTriplets(entries.begin(),entries.end());Eigen::VectorXd solution=seed;
    // L1 balls are imposed by their supporting half-spaces. Add only violated
    // facets, always solving the same final-pair QP; no alternate operator.
    std::vector<int> reactionConstants;
    if(fixed)for(int i=0;i<n;++i)if((*fixed)[i] && row[i]>=0)reactionConstants.push_back(offset[i]);
    bool bounded=false;
    for(int pass=0;pass<32;++pass){
        const Eigen::Map<const Eigen::VectorXd> lowerBound(lo.data(),lo.size()),upperBound(hi.data(),hi.size());
        const Eigen::VectorXd value=matrix*solution;bool feasible=true;
        for(int i=0;i<value.size();++i)if(value[i]<lo[i]-tolerance || value[i]>hi[i]+tolerance)feasible=false;
        if(!feasible)solution=boundedProjection(matrix,lowerBound,upperBound,seed,report,fixed!=nullptr,fixed?&reactionConstants:nullptr);
        bounded=true;for(int i=0;i<n;++i)if(row[i]>=0 && (!fixed || !(*fixed)[i])){double norm=0;for(const auto& e:incident[i])norm+=e.second*std::abs(solution[e.first]);
            if(norm>2+tolerance){bounded=false;const int r=static_cast<int>(lo.size());lo.push_back(-std::numeric_limits<double>::infinity());hi.push_back(2);
                for(const auto& e:incident[i])entries.emplace_back(r,e.first,std::copysign(e.second,solution[e.first]));}}
        if(bounded)break;matrix.resize(lo.size(),pairs.size());matrix.setFromTriplets(entries.begin(),entries.end());}
    require(bounded,"bounded coefficient facets did not converge");
    if(fixed){const Eigen::Map<const Eigen::VectorXd> lowerBound(lo.data(),lo.size()),upperBound(hi.data(),hi.size());
        certifyFinalBoundedRows(matrix,lowerBound,upperBound,solution,&reactionConstants);}
    const Eigen::VectorXd value=matrix*solution;if(!fixed)report.boundedRows=report.interiorRows=count;report.constraints=static_cast<int>(lo.size());
    for(int i=0;i<n;++i)if(row[i]>=0){double norm=0;for(const auto& e:incident[i])norm+=e.second*std::abs(solution[e.first]);
        if(fixed && (*fixed)[i])report.maxReactionCoefficientL1=std::max(report.maxReactionCoefficientL1,product(rateScale,numericScale[i]*norm));
        else report.coefficientRatio=std::max(report.coefficientRatio,static_cast<long double>(norm/2));
        for(int c=0;c<width[i];++c){const int r=offset[i]+c;if(c<3){const auto error=std::abs(value[r]-lo[r]);if(c==0){if(fixed && (*fixed)[i])report.maxReactionConstantRawError=std::max(report.maxReactionConstantRawError,static_cast<long double>(error));const double physicalScale=norm+2*std::abs(lo[r]);
                    const long double relative=fixed && (*fixed)[i]?(physicalScale>0?error/physicalScale:(error==0?0:std::numeric_limits<long double>::infinity())):error;
                    report.constantError=std::max(report.constantError,relative);}else report.affineError=std::max(report.affineError,static_cast<long double>(error));}
            else if(fixed){const double center=target[i][c]/numericScale[i],halfwidth=h[i]/length,residual=std::abs(value[r]-center);
                report.quadraticBudgetRatio=std::max(report.quadraticBudgetRatio,static_cast<long double>(residual/halfwidth));
                report.quadraticBoundError=std::max(report.quadraticBoundError,static_cast<long double>(std::max(0.,residual-halfwidth)));}
            else{report.quadraticBudgetRatio=std::max(report.quadraticBudgetRatio,static_cast<long double>(std::abs(value[r])/hi[r]));report.quadraticBoundError=std::max(report.quadraticBoundError,static_cast<long double>(std::max(0.,std::abs(value[r])-hi[r])));}}}
    require(std::max({report.constantError,report.affineError,report.quadraticBoundError})<=tolerance,"final bounded native moments failed");
    if(fixed)require(report.coefficientRatio<=1+.5L*tolerance,"final free native coefficient bound failed");
    std::vector<HydroConvectionCoupling> result;for(int j=0;j<solution.size();++j){const auto final=static_cast<long double>(solution[j])*pairScale[j];
        if(finalPairs && final!=0)finalPairs->push_back({pairs[j].first,pairs[j].second,product(rateScale,final)});
        const auto correction=final-original[pairs[j]];if(correction!=0)result.push_back({pairs[j].first,pairs[j].second,product(rateScale,correction)});}
    // A zero-gross FREE row has the exact L1 bound zero. Its original centered
    // scatter can contain third-row entries; cancel those symbolically and
    // apply the retained final K directly in the typed impermeable path.
    if(fixed)for(const auto& e:original)if(e.second!=0 && (numericScale[e.first.first]==0 || numericScale[e.first.second]==0))
        result.push_back({e.first.first,e.first.second,product(rateScale,-e.second)});
    return result;
}

}
struct DissipationPatch
{
    std::vector<int> support;std::vector<double> squareRootMassWeights;
    Eigen::MatrixXd image;
};
struct DissipationConnection
{
    int first=-1,second=-1,firstSegment=-1,secondSegment=-1;
};
namespace
{
void validateDissipation(const HydroConvectionDissipationConfig& c)
{
    require(std::isfinite(c.strength) && c.strength>=0 && c.strength<=1,
        "dissipation strength must be finite and in [0,1]");
    require(std::isfinite(c.activation) && std::isfinite(c.saturation) &&
        c.activation>=0 && c.activation<c.saturation,
        "dissipation thresholds must be finite with 0 <= activation < saturation");
}
void prepareDissipation(const HydroDualMaterial& m,const NativeTopology& topology,
    const std::vector<NativeLeafFace>& faces,const std::vector<std::vector<int>>& cellFaces,
    const std::vector<std::vector<int>>& neighbors,std::vector<DissipationPatch>& patches,
    std::vector<DissipationConnection>& connections)
{
    const int n=m.layout().velocityCount();std::vector<long double> mass(n);
    std::vector<std::vector<int>> supports(n);
    for(int i=0;i<n;++i)mass[topology.representative[i]]+=m.midpointMass().masses()[i];
    for(int i=0;i<n;++i)if(topology.representative[i]==i){
        auto cells=topology.roots[i];
        for(int c:topology.roots[i])for(int j:neighbors[c])insert(cells,j);
        for(int c:cells)for(int j:cellFaces[c])if(component(m.layout(),i)==component(m.layout(),j))
            insert(supports[i],topology.representative[j]);
    }
    const auto& segments=m.geometry().segments();
    for(std::size_t g=0;g<segments.size();++g){
        const auto& s=segments[g];int i=-1,j=-1,other=-1;
        if(internal(s)){i=s.negativeVelocity;j=s.positiveVelocity;}
        else if(topology.partner[g]>static_cast<int>(g)){other=topology.partner[g];i=owner(s);j=owner(segments[other]);}
        else continue;
        const int ri=topology.representative[i],rj=topology.representative[j];
        if(ri==rj)continue; // Identical periodic samples have exactly zero jump.
        DissipationPatch p;p.support=supports[ri];for(int k:supports[rj])insert(p.support,k);
        std::sort(p.support.begin(),p.support.end());require(!p.support.empty(),"empty dissipation sensor support");
        const auto& box=m.geometry().volumes()[i];const double h=std::max(box.upper.x-box.lower.x,box.upper.y-box.lower.y);
        const auto point=faces[i].center+displacement(faces[i].center,faces[j].center,topology)*.5;
        long double largest=0,total=0;for(int j:p.support)largest=std::max(largest,mass[j]);
        require(std::isfinite(largest) && largest>0 && h>0,"invalid dissipation sensor mass or scale");
        for(int j:p.support)total+=mass[j]/largest;
        Eigen::MatrixXd a(p.support.size(),6);
        for(std::size_t k=0;k<p.support.size();++k){
            const double w=static_cast<double>((mass[p.support[k]]/largest)/total);
            require(w>0 && std::isfinite(w),"unrepresentable dissipation sensor mass weight");p.squareRootMassWeights.push_back(std::sqrt(w));
            const auto q=polynomial(point,faces[p.support[k]].center,h,topology);
            for(int c=0;c<6;++c)a(k,c)=p.squareRootMassWeights.back()*q[c];
        }
        // A sampling projection for the sensor only. No constraint or change
        // to the final convection operator, material flux, or actual B/J.
        require(a.allFinite(),"nonfinite dissipation sensor geometry");
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(a,Eigen::ComputeThinU|Eigen::ComputeThinV);svd.setThreshold(1e-12);
        p.image=svd.matrixU().leftCols(svd.rank());require(p.image.allFinite(),"nonfinite dissipation sensor projection");
        patches.push_back(std::move(p));connections.push_back({i,j,static_cast<int>(g),other});
    }
}
long double sensorValue(const DissipationPatch& p,const std::vector<double>& u,std::vector<long double>& y)
{
    const int root=p.support.front();long double scale=0;
    for(int j:p.support)scale=std::max(scale,std::abs(.5L*u[j]-.5L*u[root]));
    if(scale==0)return 0;
    std::array<long double,6> coefficient{};long double mean=0;y.resize(p.support.size());
    for(std::size_t k=0;k<p.support.size();++k){
        y[k]=p.squareRootMassWeights[k]*(.5L*u[p.support[k]]-.5L*u[root])/scale;
        mean+=p.squareRootMassWeights[k]*y[k];
        for(int c=0;c<p.image.cols();++c)coefficient[c]+=p.image(k,c)*y[k];
    }
    long double variance=0,residual=0;
    for(std::size_t k=0;k<p.support.size();++k){
        const auto deviation=y[k]-p.squareRootMassWeights[k]*mean;variance+=deviation*deviation;
        auto r=y[k];for(int c=0;c<p.image.cols();++c)r-=p.image(k,c)*coefficient[c];residual+=r*r;
    }
    return variance>0?finite(std::sqrt(residual/variance),"nonfinite Favre dissipation sensor"):0;
}
}
struct HydroConvectionBoundary::Data
{
    Data(HydroDualMaterial m,std::vector<double> t):material(std::move(m)),traces(std::move(t)){}
    HydroDualMaterial material;std::vector<double> traces;
    std::optional<std::array<std::vector<Vector2D>,2>> vofEndpointVelocity;
    std::optional<std::vector<HydroVofBoundaryMomentumFlux>> vofOpenMomentum;
    std::vector<int> freePressurePorts;
    struct ProductionReceipt {std::vector<int> closedVelocityDofs;long double maximumNativeMeanDefect=0;
        std::vector<std::size_t> exteriorSegments;std::size_t exteriorBytes=0;};
    std::optional<ProductionReceipt> productionReceipt;
};
HydroConvectionBoundary::HydroConvectionBoundary(std::shared_ptr<const Data> d):data_(std::move(d)){}
HydroConvectionBoundary HydroConvectionBoundary::bind(const QuadtreeMesh& mesh,HydroDualMaterial m,
    const std::vector<double>& traces,HydroMaterialSources sources,LAMP::TimeLayer clock)
{
    m.validateFor(mesh,sources,clock);const auto& segments=m.geometry().segments();
    require(traces.size()==segments.size(),"boundary trace dimensions do not match shared segments");
    for(std::size_t id=0;id<segments.size();++id){const auto& g=segments[id];require(std::isfinite(traces[id]),"nonfinite boundary trace");
        if(internal(g))require(traces[id]==0,"nonzero unused internal boundary trace");
        else if(m.layout().domain().axisymmetric() && g.normalAxis==0 && g.begin.x==0 && component(m.layout(),owner(g))==0)
            require(traces[id]==0,"nonzero radial axis boundary trace");
    }
    return HydroConvectionBoundary(std::make_shared<Data>(std::move(m),traces));
}
HydroConvectionBoundary HydroConvectionBoundary::bindVofEndpointVelocity(const QuadtreeMesh& mesh,HydroDualMaterial m,
    const std::vector<Vector2D>& begin,const std::vector<Vector2D>& end,HydroMaterialSources sources,LAMP::TimeLayer clock)
{
    m.validateFor(mesh,sources,clock);require(m.material().model()==HydroPhaseGeometry::Model::GeometricVof,"VOF endpoint vectors require actual VOF material");
    const auto& volume=validateVofNativeSource(mesh,m);const auto& slab=volume.advector();const auto* boundary=slab.boundary();
    require(boundary&&!boundary->allClosed(),"VOF endpoint vectors require retained actual F1 boundary");boundary->validateFor(mesh,volume.source());
    const auto& segments=m.geometry().segments();require(begin.size()==segments.size()&&end.size()==segments.size(),"different actual endpoint full-vector dimensions");
    const auto zero=[](Vector2D v){return v.x==0&&v.y==0;};const auto coordinate=[](Vector2D v,int axis){return axis==0?v.x:v.y;};
    std::vector<double> traces(segments.size());
    for(std::size_t id=0;id<segments.size();++id){const auto& g=segments[id];
        require(std::isfinite(begin[id].x)&&std::isfinite(begin[id].y)&&std::isfinite(end[id].x)&&std::isfinite(end[id].y),"nonfinite prescribed endpoint full vector");
        const auto q=internal(g)?0.L:(g.negativeVelocity>=0?1.L:-1.L)*m.faceRates()[id];
        if(internal(g)||q>=0){require(zero(begin[id])&&zero(end[id]),"nonzero unused internal/zero/outflow endpoint vector");continue;}
        require(g.primalVelocity>=0,"incoming dual segment lacks actual primal native port");const int row=m.layout().nativeFace(g.primalVelocity);
        const auto face=m.layout().domain().face(mesh,row);require((face.negativeCell<0||face.positiveCell<0)&&face.area>0,"incoming endpoint vector has no positive-area exterior native face");
        require((face.horizontal?1:0)==g.normalAxis,"incoming segment/native normal axis differs");
        using Mode=Interface::GeometricVofAdvectorSlab::Boundary::Mode;const auto mode=boundary->rules()[row];
        require(mode==Mode::PureLiquidReservoir||mode==Mode::PureGasReservoir,"incoming vector requires actual F1 pure reservoir");
        require(coordinate(begin[id],g.normalAxis)==slab.beginNormal()[row]&&coordinate(end[id],g.normalAxis)==slab.endNormal()[row],
            "prescribed full-vector normal differs from actual F1 endpoint native mean");
        const int transported=component(m.layout(),owner(g));using Dyadic=Geometry::Detail::FiniteVolume::Dyadic;
        Dyadic midpoint;midpoint.add((Dyadic::product)({.5L,coordinate(begin[id],transported)}));
        midpoint.add((Dyadic::product)({.5L,coordinate(end[id],transported)}));
        const double stored=Geometry::Detail::FiniteVolume::roundDouble(midpoint).value;
        require(std::isfinite(stored)&&(midpoint.compare(Dyadic{})==0||stored!=0),"unrepresentable prescribed midpoint component");
        traces[id]=stored;
    }
    auto data=std::make_shared<Data>(std::move(m),std::move(traces));data->vofEndpointVelocity=std::array<std::vector<Vector2D>,2>{{begin,end}};
    return HydroConvectionBoundary(std::move(data));
}
bool HydroConvectionBoundary::hasVofEndpointVelocity() const noexcept{return data_&&data_->vofEndpointVelocity.has_value();}
const std::vector<Vector2D>& HydroConvectionBoundary::beginVelocity() const
{require(hasVofEndpointVelocity(),"boundary has no actual VOF endpoint full vectors");return (*data_->vofEndpointVelocity)[0];}
const std::vector<Vector2D>& HydroConvectionBoundary::endVelocity() const
{require(hasVofEndpointVelocity(),"boundary has no actual VOF endpoint full vectors");return (*data_->vofEndpointVelocity)[1];}
const HydroDualMaterial& HydroConvectionBoundary::material() const noexcept{return data_->material;}
const std::vector<double>& HydroConvectionBoundary::traces() const noexcept{return data_->traces;}
void HydroConvectionBoundary::validateFor(const QuadtreeMesh& mesh,const HydroDualMaterial& expected,
    HydroMaterialSources sources,LAMP::TimeLayer clock) const
{
    require(&material().faceRates()==&expected.faceRates(),"different prepared material owner");
    material().validateFor(mesh,sources,clock);
}
namespace {
using OpenProducts=Geometry::Detail::FiniteVolume::Binary64ExactProducts;
struct OpenInteriorProducts
{
    std::vector<OpenProducts::Factor> diagonal,pairs;
    std::size_t bytes=0;
};
}
struct HydroKepConvection::Data
{
    explicit Data(HydroDualMaterial m):material(std::move(m)){}
    HydroDualMaterial material;std::vector<Vector2D> points;std::vector<std::vector<HydroConvectionWeight>> weights;
    std::vector<HydroConvectionCoupling> couplings;HydroConvectionMomentReport momentReport;
    std::optional<HydroAdvectorBoundary> impermeable,prescribedWalls;
    std::optional<HydroConvectionBoundary> prescribedInlet;
    bool vofPrescribed=false,vofOpenMomentum=false,productionPrepared=false;
    std::vector<Detail::HydroFiniteContraction> openOutgoing,openIncomingContinuation;
    mutable std::shared_ptr<const OpenApplyPlan> openApplyPlan;
    mutable bool openApplyPlanAttempted=false;
    mutable std::mutex openApplyPlanMutex;
    std::shared_ptr<const Detail::HydroOpenConvectionSupport> openSupport;
    std::vector<long double> openPairCoefficients;
    std::shared_ptr<const OpenInteriorProducts> openInteriorProducts;
    std::vector<int> openBoundarySlot;
    int boundarySlot(int i) const noexcept {return openSupport?openSupport->boundaryRowSlot[i]:openBoundarySlot[i];}
    std::vector<long double> externalOutflow,externalAbsolute,boundaryLoad;
    std::vector<int> fixed;
    std::vector<HydroConvectionCoupling> finalPairs;
    std::vector<HydroConvectionPeriodicPair> periodicPairs;NativeTopology topology;std::vector<long double> outflow;
    HydroConvectionDissipationConfig dissipation;bool sensorDissipation=false;
    std::vector<DissipationPatch> sensorPatches;std::vector<DissipationConnection> dissipationConnections;
};
HydroKepConvection::HydroKepConvection(std::shared_ptr<const Data> d):data_(std::move(d)){}
HydroKepConvection HydroKepConvection::build(const QuadtreeMesh& mesh,HydroDualMaterial m)
{return buildWithPorts(mesh,std::move(m),nullptr);}
HydroKepConvection HydroKepConvection::build(const QuadtreeMesh& mesh,HydroDualMaterial m,const std::vector<HydroConvectionPeriodicPair>& pairs)
{return buildWithPorts(mesh,std::move(m),&pairs);}
HydroKepConvection HydroKepConvection::buildControlled(const QuadtreeMesh& mesh,HydroDualMaterial m,HydroConvectionDissipationConfig config)
{validateDissipation(config);return buildWithPorts(mesh,std::move(m),nullptr,&config);}
HydroKepConvection HydroKepConvection::buildControlled(const QuadtreeMesh& mesh,HydroDualMaterial m,const std::vector<HydroConvectionPeriodicPair>& pairs,HydroConvectionDissipationConfig config)
{validateDissipation(config);return buildWithPorts(mesh,std::move(m),&pairs,&config);}
HydroKepConvection HydroKepConvection::buildImpermeable(const QuadtreeMesh& mesh,HydroDualMaterial m,const HydroAdvectorBoundary& boundary)
{return buildWithPorts(mesh,std::move(m),nullptr,nullptr,&boundary);}
HydroKepConvection HydroKepConvection::buildImpermeableControlled(const QuadtreeMesh& mesh,HydroDualMaterial m,const HydroAdvectorBoundary& boundary,HydroConvectionDissipationConfig config)
{validateDissipation(config);return buildWithPorts(mesh,std::move(m),nullptr,&config,&boundary);}
HydroKepConvection HydroKepConvection::buildPrescribedFlow(const QuadtreeMesh& mesh,HydroDualMaterial m,
    const HydroAdvectorBoundary& walls,const HydroConvectionBoundary& inlet)
{return buildPrescribedFlow(mesh,std::move(m),walls,inlet,{});}
HydroKepConvection HydroKepConvection::buildPrescribedFlow(const QuadtreeMesh& mesh,HydroDualMaterial m,
    const HydroAdvectorBoundary& walls,const HydroConvectionBoundary& inlet,
    const std::vector<int>& freeOutflow)
{return buildPrescribedFlowImpl(mesh,std::move(m),walls,inlet,freeOutflow,false);}
HydroKepConvection HydroKepConvection::buildVofPrescribedFlow(const QuadtreeMesh& mesh,HydroDualMaterial m,
    const HydroAdvectorBoundary& walls,const HydroConvectionBoundary& inlet,const std::vector<int>& freeOutflow)
{return buildPrescribedFlowImpl(mesh,std::move(m),walls,inlet,freeOutflow,true);}
HydroKepConvection HydroKepConvection::buildPrescribedFlowImpl(const QuadtreeMesh& mesh,HydroDualMaterial m,
    const HydroAdvectorBoundary& walls,const HydroConvectionBoundary& inlet,const std::vector<int>& freeOutflow,bool vof)
{
    m.validateFor(mesh,m.sources(),m.timeLayer());
    ProductDiagnosticScope productDiagnostic(m.sources(),m.timeLayer(),vof?"buildVofPrescribedFlow":"buildPrescribedFlow");
    require(!inlet.hasVofOpenMomentum(),"ambient momentum requires explicit buildVofOpenFlow");
    auto fixed=prescribedFlowRows(mesh,m,walls,inlet,freeOutflow,vof);
    auto d=std::make_shared<Data>(std::move(m));const auto& l=d->material.layout();
    d->prescribedWalls=walls;d->prescribedInlet=inlet;d->fixed=std::move(fixed);
    d->vofPrescribed=vof;std::map<int,Detail::HydroFiniteContraction> vofBoundaryLoad;
    d->externalOutflow.resize(l.velocityCount());d->externalAbsolute.resize(l.velocityCount());
    d->boundaryLoad.resize(l.velocityCount());d->outflow.resize(l.velocityCount());
    std::vector<NativeLeafFace> faces;std::vector<std::vector<int>> neighbors(l.pressureCount());
    for(int u=0;u<l.velocityCount();++u){const auto f=l.domain().face(mesh,l.nativeFace(u));faces.push_back(f);
        const int a=f.negativeCell<0?-1:l.pressureDof(f.negativeCell),b=f.positiveCell<0?-1:l.pressureDof(f.positiveCell);
        if(a>=0 && b>=0){insert(neighbors[a],b);insert(neighbors[b],a);}}
    const auto& segments=d->material.geometry().segments();
    for(std::size_t g=0;g<segments.size();++g){const auto& segment=segments[g];const auto f=d->material.faceRates()[g];
        d->points.push_back(center(segment,l.domain().axisymmetric()));
        if(internal(segment))d->weights.push_back({{segment.negativeVelocity,.5},{segment.positiveVelocity,.5}});
        else {d->weights.emplace_back();const int i=owner(segment);
            const auto q=finite((segment.negativeVelocity>=0?1.L:-1.L)*f,"nonfinite prescribed-flow exterior rate");
            d->externalOutflow[i]=finite(d->externalOutflow[i]+q,"nonfinite prescribed-flow exterior sum");
            d->externalAbsolute[i]=finite(d->externalAbsolute[i]+std::abs(q),"nonfinite prescribed-flow exterior work");
            if(q<0){if(vof)vofBoundaryLoad[i].addProduct({-q,inlet.traces()[g]});
                else d->boundaryLoad[i]=finite(d->boundaryLoad[i]-product(q,inlet.traces()[g]),
                    "nonfinite prescribed-flow inlet momentum load");}}
        if(segment.negativeVelocity>=0)d->outflow[segment.negativeVelocity]+=f;
        if(segment.positiveVelocity>=0)d->outflow[segment.positiveVelocity]-=f;
    }
    if(vof)for(const auto& value:vofBoundaryLoad)d->boundaryLoad[value.first]=value.second.rounded();
    d->topology=topologyFor(d->material,faces,neighbors,nullptr);
    std::vector<bool> fixedMask(l.velocityCount(),false);for(int u:d->fixed)fixedMask[u]=true;
    d->finalPairs=physicalCenteredPairs(d->material,fixedMask,d->momentReport,&d->externalOutflow);
    if(vof)certifyVofFreeOutflowMeans(mesh,d->material,freeOutflow,d->momentReport);
    else certifyFreeOutflowMeans(d->material.material().midpointRate().spatialAdvector(),freeOutflow,d->momentReport);
    for(int i=0;i<l.velocityCount();++i){const auto canonical=finite(d->material.cells()[i].outwardIntegral/d->material.deltaTime(),
            "nonfinite canonical prescribed-flow outflow");
        const auto difference=finite(d->outflow[i]-canonical,"nonfinite material rate/integral rounding difference");
        auto& report=d->momentReport;report.rateSumMinusIntegralRateSigned+=difference;
        report.rateSumMinusIntegralRateL1+=std::abs(difference);
        report.maxAbsRateSumMinusIntegralRate=std::max(report.maxAbsRateSumMinusIntegralRate,std::abs(difference));
        d->outflow[i]=canonical;
    }
    return HydroKepConvection(std::move(d));
}
HydroKepConvection HydroKepConvection::buildWithPorts(const QuadtreeMesh& mesh,HydroDualMaterial m,const std::vector<HydroConvectionPeriodicPair>* pairs,const HydroConvectionDissipationConfig* config,const HydroAdvectorBoundary* boundary)
{
    const auto checked=CheckedVofMaterialBatch::validate(mesh,m,m.sources(),m.timeLayer());
    ProductDiagnosticScope productDiagnostic(m.sources(),m.timeLayer(),"buildWithPorts");
    std::vector<int> fixed;
    require(m.material().model()!=HydroPhaseGeometry::Model::GeometricVof||(boundary&&!pairs),
        "volume-transport KEP requires the typed complete impermeable path; opaque/periodic preparation is not supported");
    if(boundary){require(!pairs,"impermeable source cannot declare periodic aliases");fixed=impermeableRows(mesh,m,*boundary,&checked);}
    auto d=std::make_shared<Data>(std::move(m));const auto& l=d->material.layout();
    if(boundary){d->impermeable=*boundary;d->fixed=std::move(fixed);}
    std::vector<NativeLeafFace> faces;std::vector<std::vector<int>> cellFaces(l.pressureCount()),neighbors(l.pressureCount());
    for(int u=0;u<l.velocityCount();++u){const auto f=l.domain().face(mesh,l.nativeFace(u));faces.push_back(f);
        const int a=f.negativeCell<0?-1:l.pressureDof(f.negativeCell),b=f.positiveCell<0?-1:l.pressureDof(f.positiveCell);
        if(a>=0)cellFaces[a].push_back(u);if(b>=0)cellFaces[b].push_back(u);
        if(a>=0 && b>=0){insert(neighbors[a],b);insert(neighbors[b],a);}
    }
    for(const auto& g:d->material.geometry().segments()){
        const auto p=center(g,l.domain().axisymmetric());d->points.push_back(p);
        if(boundary && internal(g))d->weights.push_back({{g.negativeVelocity,.5},{g.positiveVelocity,.5}});
        else d->weights.push_back(internal(g)?reconstruct(g,p,d->material.geometry(),faces,cellFaces,neighbors):std::vector<HydroConvectionWeight>{});
    }
    d->topology=topologyFor(d->material,faces,neighbors,pairs);if(pairs)d->periodicPairs=*pairs;
    if(boundary){std::vector<bool> fixedMask(l.velocityCount(),false);for(int u:d->fixed)fixedMask[u]=true;
        d->finalPairs=physicalCenteredPairs(d->material,fixedMask,d->momentReport);}
    else d->couplings=d->topology.paired?constrainBoundedMoments(d->material,faces,cellFaces,neighbors,d->points,d->weights,d->topology,d->momentReport):
        constrainOpaqueMoments(d->material,faces,cellFaces,neighbors,d->points,d->weights,d->momentReport);
    d->outflow.resize(l.velocityCount());for(std::size_t g=0;g<d->material.faceRates().size();++g){const auto& s=d->material.geometry().segments()[g];const auto f=d->material.faceRates()[g];if(s.negativeVelocity>=0)d->outflow[s.negativeVelocity]+=f;if(s.positiveVelocity>=0)d->outflow[s.positiveVelocity]-=f;}
    if(boundary)for(int i=0;i<l.velocityCount();++i){const auto canonical=finite(d->material.cells()[i].outwardIntegral/d->material.deltaTime(),"nonfinite canonical impermeable outflow");
        const auto difference=finite(d->outflow[i]-canonical,"nonfinite material rate/integral rounding difference");
        auto& report=d->momentReport;report.rateSumMinusIntegralRateSigned+=difference;report.rateSumMinusIntegralRateL1+=std::abs(difference);
        report.maxAbsRateSumMinusIntegralRate=std::max(report.maxAbsRateSumMinusIntegralRate,std::abs(difference));d->outflow[i]=canonical;}
    if(config && config->strength>0){d->sensorDissipation=true;d->dissipation=*config;prepareDissipation(d->material,d->topology,faces,cellFaces,neighbors,d->sensorPatches,d->dissipationConnections);}
    return HydroKepConvection(std::move(d));
}
const HydroDualMaterial& HydroKepConvection::material() const noexcept{return data_->material;}
const HydroDofLayout& HydroKepConvection::layout() const noexcept{return material().layout();}
const std::vector<Vector2D>& HydroKepConvection::points() const noexcept{return data_->points;}
const std::vector<std::vector<HydroConvectionWeight>>& HydroKepConvection::weights() const noexcept{return data_->weights;}
const std::vector<HydroConvectionCoupling>& HydroKepConvection::couplings() const noexcept{return data_->couplings;}
const HydroConvectionMomentReport& HydroKepConvection::momentReport() const noexcept{return data_->momentReport;}
bool HydroKepConvection::isImpermeable() const noexcept{return data_->impermeable.has_value();}
bool HydroKepConvection::isPrescribedFlow() const noexcept{return data_->prescribedInlet.has_value();}
const std::vector<int>& HydroKepConvection::fixedVelocityDofs() const noexcept{return data_->fixed;}
const std::vector<HydroConvectionPeriodicPair>& HydroKepConvection::periodicPairs() const noexcept{return data_->periodicPairs;}
const std::vector<int>& HydroKepConvection::representatives() const noexcept{return data_->topology.representative;}
HydroKineticHodge::Matrix HydroKepConvection::frozenMomentumOperator(const QuadtreeMesh& mesh,
    HydroMaterialSources source,LAMP::TimeLayer clock) const
{
    const auto checked=CheckedVofMaterialBatch::validate(mesh,material(),source,clock);
    require(data_->impermeable.has_value() || data_->prescribedInlet.has_value(),
        "frozen momentum requires typed impermeable source or typed prescribed-flow source");
    require(!data_->sensorDissipation,"frozen momentum rejects sensor dissipation");
    const auto& walls=data_->impermeable?*data_->impermeable:*data_->prescribedWalls;
    validatePhysicalSourceBoundary(mesh,material(),walls,data_->impermeable.has_value(),&checked);

    if(data_->prescribedInlet)validateConvectionBoundaryAfterMaterial(mesh,material(),*data_->prescribedInlet,source,clock,&checked);
    return frozenMomentumOperatorAfterValidation();
}
HydroKineticHodge::Matrix HydroKepConvection::frozenMomentumOperatorAfterValidation() const
{
    const auto entry=[](long double value){
        finite(value,"nonfinite frozen momentum coefficient");
        const double stored=static_cast<double>(value);
        require(std::isfinite(stored) && (value==0 || stored!=0),
            "unrepresentable frozen momentum coefficient");return stored;};
    const int n=layout().velocityCount();
    std::vector<Eigen::Triplet<double,Eigen::Index>> entries;
    entries.reserve(static_cast<std::size_t>(n)+2*data_->finalPairs.size());
    const auto dt=material().deltaTime();
    for(int i=0;i<n;++i){
        long double diagonalValue;
        if(data_->vofOpenMomentum){Detail::HydroFiniteContraction diagonal;
            diagonal.addProduct({.5L,data_->outflow[i]});diagonal.addProduct({-.5L,data_->externalOutflow[i]});
            const int slot=data_->boundarySlot(i);if(slot>=0){diagonal.add(data_->openOutgoing[slot]);diagonal.add(data_->openIncomingContinuation[slot],-1);}
            if(!material().usesConservativeVofMomentum())diagonal.addQuotientProduct({-1.L,material().cells()[i].integratedDefect},dt);
            diagonalValue=diagonal.rounded();}
        else {const auto consistency=material().usesConservativeVofMomentum()?0.L:finite(material().cells()[i].integratedDefect/dt,
                "nonfinite frozen material consistency rate");
            diagonalValue=data_->prescribedInlet
                ? .5L*data_->outflow[i]-consistency+.5L*data_->externalAbsolute[i]
                : .5L*data_->outflow[i]-consistency;}
        const double diagonal=entry(diagonalValue);
        if(diagonal!=0)entries.emplace_back(i,i,diagonal);
    }
    for(const auto& pair:data_->finalPairs){const double coefficient=entry(pair.coefficient);
        if(coefficient!=0){entries.emplace_back(pair.first,pair.second,coefficient);
            entries.emplace_back(pair.second,pair.first,-coefficient);}}
    HydroKineticHodge::Matrix result(n,n);
    result.setFromTriplets(entries.begin(),entries.end());result.makeCompressed();return result;
}
std::vector<long double> HydroKepConvection::frozenBoundaryMomentumLoad(const QuadtreeMesh& mesh,
    HydroMaterialSources source,LAMP::TimeLayer clock) const
{
    const auto checked=CheckedVofMaterialBatch::validate(mesh,material(),source,clock);
    require(data_->impermeable.has_value() || data_->prescribedInlet.has_value(),
        "frozen boundary momentum requires typed physical source");
    require(!data_->sensorDissipation,"frozen boundary momentum rejects sensor dissipation");
    const auto& walls=data_->impermeable?*data_->impermeable:*data_->prescribedWalls;
    validatePhysicalSourceBoundary(mesh,material(),walls,data_->impermeable.has_value(),&checked);

    if(!data_->prescribedInlet)return std::vector<long double>(layout().velocityCount(),0.L);
    validateConvectionBoundaryAfterMaterial(mesh,material(),*data_->prescribedInlet,source,clock,&checked);
    return data_->boundaryLoad;
}
std::pair<HydroKineticHodge::Matrix,std::vector<long double>> HydroKepConvection::prepareFrozenMomentum(
    const QuadtreeMesh& mesh,const HydroDualMaterial& validatedMaterial,
    HydroMaterialSources source,LAMP::TimeLayer clock) const
{
    const auto expected=material().sources();const auto expectedClock=material().timeLayer();
    const bool sameSource=source.model==expected.model&&source.begin==expected.begin&&source.midpoint==expected.midpoint
        &&source.end==expected.end&&source.advector==expected.advector&&source.vofStamp&&expected.vofStamp
        &&source.vofStamp->sameInstance(*expected.vofStamp)&&source.vofIds[0]==expected.vofIds[0]
        &&source.vofIds[1]==expected.vofIds[1]&&source.vofIds[2]==expected.vofIds[2];
    const bool sameClock=clock.step==expectedClock.step&&clock.begin==expectedClock.begin
        &&clock.end==expectedClock.end&&clock.position==expectedClock.position;
    if(data_->productionPrepared&&sameSource&&sameClock&&wallExecutionQualified()
        &&&material().cells()==&validatedMaterial.cells()
        &&material().geometry().sameIdentity(validatedMaterial.geometry())
        &&layout().sameIdentity(validatedMaterial.layout())){
        // The dedicated production factory sealed these immutable wall/inlet
        // contents. Frozen just validated this exact material and live mesh;
        // no callback or mutable mesh operation occurs before this private use.
        require(data_->impermeable.has_value()||data_->prescribedInlet.has_value(),
            "frozen momentum requires typed impermeable source or typed prescribed-flow source");
        require(!data_->sensorDissipation,"frozen momentum rejects sensor dissipation");
        auto matrix=frozenMomentumOperatorAfterValidation();
        return {std::move(matrix),data_->boundaryLoad};
    }
    // Foreign/unqualified and legacy owners retain their original validation
    // order. Only the immutable boundary-load copy moves before the A action.
    auto matrix=frozenMomentumOperator(mesh,source,clock);
    auto boundary=frozenBoundaryMomentumLoad(mesh,source,clock);
    return {std::move(matrix),std::move(boundary)};
}
struct HydroKepConvectionResult::Data
{
    Data(HydroDualMaterial m,HydroConvectionBoundary b):material(std::move(m)),boundary(std::move(b)){}
    HydroDualMaterial material;HydroConvectionBoundary boundary;std::vector<long double> center,rates,loads,dissipation;HydroKepConvectionReport report;
    bool centerAvailable=true,powerAvailable=true,componentAvailable=true;
    std::optional<HydroKepOpenBoundaryBudget> boundaryBudget;
    std::vector<HydroConvectionDissipationSample> dissipationSamples;
    HydroKepOpenApplyReport openApply;
};
HydroKepConvectionResult::HydroKepConvectionResult(std::shared_ptr<const Data> d):data_(std::move(d)){}
HydroKepConvectionResult HydroKepConvection::apply(const QuadtreeMesh& mesh,const HydroNativeValues& values,
    const HydroConvectionBoundary& boundary,HydroMaterialSources source,LAMP::TimeLayer clock) const
{
    const auto checked=CheckedVofMaterialBatch::validate(mesh,material(),source,clock);values.validateFor(layout(),mesh);
    validateConvectionBoundaryAfterMaterial(mesh,material(),boundary,source,clock,&checked);
    require(!data_->vofOpenMomentum,"ambient momentum requires explicit applyVofOpen");
    ProductDiagnosticScope productDiagnostic(source,clock,"apply");
    if(data_->impermeable)validatePhysicalSourceBoundary(mesh,material(),*data_->impermeable,true,&checked);

    if(data_->prescribedInlet){validateWallsAfterMaterial(mesh,material(),*data_->prescribedWalls,&checked);
        validateConvectionBoundaryAfterMaterial(mesh,material(),*data_->prescribedInlet,source,clock,&checked);
        validatePhysicalSourceBoundary(mesh,material(),*data_->prescribedWalls,false,&checked);
        for(int u:data_->prescribedWalls->velocityDofs())
            require(values.normalVelocity()[u]==0,"nonzero applied prescribed-flow wall mean");}
    const bool finiteVofClosed=material().material().model()==HydroPhaseGeometry::Model::GeometricVof
        && data_->impermeable.has_value() && !data_->topology.paired && data_->sensorPatches.empty();
    const bool finiteVofPhysical=finiteVofClosed||(material().material().model()==HydroPhaseGeometry::Model::GeometricVof
        &&data_->vofPrescribed&&!data_->topology.paired&&data_->sensorPatches.empty());
    const bool conservative=material().usesConservativeVofMomentum();
    const auto ordinary=[&]()->HydroKepConvectionResult{
    const auto& u=values.normalVelocity();
    const auto impulse=conservative?std::vector<long double>(u.size()):
        (checked.matches(material())?material().consistencyImpulseAfterValidation(values)
        :material().consistencyImpulse(mesh,values,source,clock));
    if(finiteVofPhysical)for(std::size_t i=0;i<u.size();++i){
        // These original nested operands are otherwise invisible to product().
        Detail::hydroFiniteProductRange(.5L,u[i],.5L*u[i]);
        Detail::hydroFiniteProductRange(.5L,data_->outflow[i],.5L*data_->outflow[i]);
        Detail::hydroFiniteQuotientRange(impulse[i],material().deltaTime(),impulse[i]/material().deltaTime());
        if(conservative)Detail::hydroFiniteQuotientRange(material().cells()[i].integratedDefect,
            material().deltaTime(),material().cells()[i].defectRate);
        const auto change=material().cells()[i].massChange;
        Detail::hydroFiniteQuotientRange(change,material().deltaTime(),change/material().deltaTime());
    }
    if(data_->topology.paired){
        for(int i=0;i<layout().velocityCount();++i){const int r=representatives()[i];if(i!=r)identity(u[i]-static_cast<long double>(u[r]),std::abs(static_cast<long double>(u[i]))+std::abs(static_cast<long double>(u[r])),"different periodic alias velocities");}
        for(const auto& p:periodicPairs()){const int i=owner(material().geometry().segments()[p.firstSegment]),j=owner(material().geometry().segments()[p.secondSegment]);const auto trace=.5L*u[i]+.5L*u[j];
            for(int g:{p.firstSegment,p.secondSegment})identity(boundary.traces()[g]-trace,std::abs(trace)+std::abs(static_cast<long double>(boundary.traces()[g])),"different periodic port trace");}
    }
    auto d=std::make_shared<HydroKepConvectionResult::Data>(material(),boundary);d->rates.resize(u.size());d->loads.resize(u.size());
    std::array<long double,2> momentumScale{},powerScale{};const auto& segments=material().geometry().segments();
    for(std::size_t id=0;id<segments.size();++id){const auto& g=segments[id];const auto f=material().faceRates()[id];long double trace=0;
        if(internal(g))for(const auto& w:weights()[id])trace+=product(w.weight,u[w.velocity]);else trace=boundary.traces()[id];
        if(data_->prescribedInlet && !internal(g)){const auto q=(g.negativeVelocity>=0?1.L:-1.L)*f;
            const double expected=q<0?data_->prescribedInlet->traces()[id]:q>0?u[owner(g)]:0.;
            require(boundary.traces()[id]==expected,"different prescribed-flow upwind boundary trace");}
        const auto raw=product(f,finite(trace,"nonfinite common velocity trace"));d->center.push_back(raw);
        if(data_->impermeable || (data_->prescribedInlet && internal(g)))continue;
        if(internal(g)){
            const int i=g.negativeVelocity,j=g.positiveVelocity,a=component(layout(),i);const auto halfJump=.5L*u[i]-.5L*u[j];
            const auto left=product(f,.5L*trace+.5L*u[i]),right=product(f,.5L*trace+.5L*u[j]);d->rates[i]+=left;d->rates[j]-=right;
            momentumScale[a]+=std::abs(left)+std::abs(right);powerScale[a]+=std::abs(product(u[i],left))+std::abs(product(u[j],right));
            const auto factor=product(f,halfJump);for(const auto& w:weights()[id]){const auto correction=product(factor,w.weight);d->rates[w.velocity]-=correction;
                momentumScale[a]+=std::abs(correction);powerScale[a]+=std::abs(product(u[w.velocity],correction));}
        }else{
            const int i=owner(g),a=component(layout(),i);const long double sign=g.negativeVelocity>=0?1:-1;
            const auto flux=sign*raw;d->rates[i]+=flux;d->report.components[a].boundaryMomentumRate+=flux;
            // Multiply by the mass rate before the second velocity, avoiding
            // an unrepresentable u*u intermediate when the physical power fits.
            const auto work=product(u[i],flux-product(sign*f,.5L*u[i]));d->report.components[a].boundaryKineticPower+=work;
            if(data_->prescribedInlet){const auto q=sign*f;
                const auto physical=product(.5L*trace,flux);
                const auto halfJump=.5L*u[i]-.5L*trace;
                const auto dissipation=q<0?finite(2*product(halfJump,product(-q,halfJump)),
                    "nonfinite inlet upwind dissipation"):0.L;
                d->report.components[a].boundaryPhysicalKineticPower+=physical;
                d->report.components[a].inletDissipationPower+=dissipation;
                require(dissipation>=0,"negative inlet upwind dissipation");
                powerScale[a]+=std::abs(physical)+dissipation;
            }
            momentumScale[a]+=std::abs(flux);powerScale[a]+=std::abs(work);
        }
    }
    if(data_->impermeable){
        // Retain the solved FINAL shared K. Direct application makes a
        // certified zero-gross row exactly zero, without subtracting two
        // centered terms or altering the unique material F/R.
        momentumScale.fill(0);powerScale.fill(0);
        for(int i=0;i<layout().velocityCount();++i){const int a=component(layout(),i);const auto value=product(.5L*data_->outflow[i],u[i]);
            d->rates[i]=value;momentumScale[a]+=std::abs(value);powerScale[a]+=std::abs(product(u[i],value));}
        for(const auto& e:data_->finalPairs){const int i=e.first,j=e.second,a=component(layout(),i);
            const auto left=product(e.coefficient,u[j]),right=product(e.coefficient,u[i]);d->rates[i]+=left;d->rates[j]-=right;
            momentumScale[a]+=std::abs(left)+std::abs(right);powerScale[a]+=std::abs(product(u[i],left))+std::abs(product(u[j],right));}
    }
    else if(data_->prescribedInlet){
        for(int i=0;i<layout().velocityCount();++i){const int a=component(layout(),i);
            const auto internalOutflow=finite(data_->outflow[i]-data_->externalOutflow[i],
                "nonfinite prescribed-flow internal outflow");
            const auto value=product(.5L*internalOutflow,u[i]);d->rates[i]+=value;
            momentumScale[a]+=std::abs(value);powerScale[a]+=std::abs(product(u[i],value));}
        for(const auto& e:data_->finalPairs){const int i=e.first,j=e.second,a=component(layout(),i);
            const auto left=product(e.coefficient,u[j]),right=product(e.coefficient,u[i]);
            d->rates[i]+=left;d->rates[j]-=right;
            momentumScale[a]+=std::abs(left)+std::abs(right);
            powerScale[a]+=std::abs(product(u[i],left))+std::abs(product(u[j],right));}
    }
    else for(const auto& e:couplings()){
        const int i=e.first,j=e.second,a=component(layout(),i);const auto left=product(e.coefficient,u[j]),right=product(e.coefficient,u[i]);
        d->rates[i]+=left;d->rates[j]-=right;momentumScale[a]+=std::abs(left)+std::abs(right);
        powerScale[a]+=std::abs(product(u[i],left))+std::abs(product(u[j],right));
    }
    d->dissipation.resize(u.size());
    if(!data_->sensorPatches.empty()){
        std::vector<long double> scratch;
        const auto& config=data_->dissipation;
        for(std::size_t k=0;k<data_->dissipationConnections.size();++k){const auto& e=data_->dissipationConnections[k];
            const int i=e.first,j=e.second,a=component(layout(),i);
            HydroConvectionDissipationSample sample;sample.first=i;sample.second=j;sample.firstSegment=e.firstSegment;sample.secondSegment=e.secondSegment;
            sample.sensor=sensorValue(data_->sensorPatches[k],u,scratch);
            const auto ramp=sample.sensor<=config.activation?0.L:sample.sensor>=config.saturation?1.L:
                (sample.sensor-config.activation)/(config.saturation-config.activation);
            sample.chi=config.strength*ramp*ramp;sample.halfJump=.5L*u[i]-.5L*u[j];
            sample.momentumRate=product(product(sample.chi,std::abs(material().faceRates()[e.firstSegment])),sample.halfJump);
            sample.power=finite(2*product(sample.halfJump,sample.momentumRate),"nonfinite dissipation power");
            require(sample.power>=0,"negative dissipation power");
            d->dissipation[i]+=sample.momentumRate;d->dissipation[j]-=sample.momentumRate;
            d->rates[i]+=sample.momentumRate;d->rates[j]-=sample.momentumRate;
            auto& r=d->report.components[a];r.dissipationPower+=sample.power;r.maximumSensor=std::max(r.maximumSensor,sample.sensor);
            ++r.dissipationConnections;if(sample.chi>0 && material().faceRates()[e.firstSegment]!=0)++r.activeDissipationConnections;
            momentumScale[a]+=2*std::abs(sample.momentumRate);
            powerScale[a]+=std::abs(product(u[i],sample.momentumRate))+std::abs(product(u[j],sample.momentumRate))+sample.power;
            d->dissipationSamples.push_back(sample);
        }
    }
    std::vector<long double> acceleration(u.size()),groupMass(u.size()),groupDissipation(u.size());
    if(data_->topology.paired){
        for(int i=0;i<layout().velocityCount();++i){const int r=representatives()[i];if(data_->topology.count[r]>1){acceleration[r]+=-d->rates[i]+product(u[i],data_->outflow[i]);groupMass[r]+=material().midpointMass().masses()[i];groupDissipation[r]+=d->dissipation[i];}}
        for(int r=0;r<layout().velocityCount();++r)if(data_->topology.count[r]>1)acceleration[r]=finite(acceleration[r]/finite(groupMass[r],"nonfinite periodic class mass"),"nonfinite periodic acceleration");
        for(int i=0;i<layout().velocityCount();++i){const int r=representatives()[i];if(data_->topology.count[r]<=1)continue;const int a=component(layout(),i);
            const auto massLoad=product(material().midpointMass().masses()[i],acceleration[r]),outward=product(u[i],data_->outflow[i]);d->rates[i]=finite(outward-massLoad,"nonfinite periodic native convection");
            d->dissipation[i]=product(material().midpointMass().masses()[i]/groupMass[r],groupDissipation[r]);
            momentumScale[a]+=std::abs(outward)+std::abs(massLoad);powerScale[a]+=std::abs(product(u[i],outward))+std::abs(product(u[i],massLoad));}
    }
    for(int i=0;i<layout().velocityCount();++i){const int a=component(layout(),i);auto& report=d->report.components[a];
        const auto consistency=finite(impulse[i]/material().deltaTime(),"nonfinite consistency rate");
        d->loads[i]=finite(-finite(d->rates[i],"nonfinite compatible convection")+consistency,"nonfinite material momentum load");
        const int r=representatives()[i];if(!conservative && data_->topology.paired && data_->topology.count[r]>1)
            d->loads[i]=finite(product(material().midpointMass().masses()[i],acceleration[r])+product(u[i],material().cells()[i].massChange/material().deltaTime()),"nonfinite periodic material load");
        const auto power=product(u[i],d->rates[i]);
        const auto sourcePower=conservative
            ?product(-.5L*u[i],product(u[i],material().cells()[i].defectRate))
            :product(.5L*u[i],consistency);
        const auto massPower=product(.5L*u[i],product(u[i],material().cells()[i].massChange/material().deltaTime()));
        const auto kinetic=finite(product(u[i],d->loads[i])-massPower,"nonfinite material kinetic rate");
        report.momentumRate+=d->rates[i];report.convectivePower+=power;report.materialKineticPower+=sourcePower;report.kineticRate+=kinetic;
        momentumScale[a]+=std::abs(d->rates[i]);powerScale[a]+=std::abs(power)+std::abs(2*sourcePower)+std::abs(massPower)+std::abs(kinetic);
    }
    for(int a=0;a<2;++a){auto& r=d->report.components[a];
        if(data_->prescribedInlet)identity(r.boundaryKineticPower-r.boundaryPhysicalKineticPower-r.inletDissipationPower,
            powerScale[a],"prescribed-flow boundary kinetic decomposition failed");
        r.momentumIdentity=identity(r.momentumRate-r.boundaryMomentumRate,momentumScale[a],"compatible momentum balance failed");
        r.kineticIdentity=identity(r.kineticRate-r.materialKineticPower+r.boundaryKineticPower+r.dissipationPower,powerScale[a],"nonconservative kinetic identity failed");
    }
    return HydroKepConvectionResult(std::move(d));
    };
    if(!finiteVofPhysical)return ordinary();
    try{Detail::HydroFiniteRangeScope range(true);return ordinary();}
    catch(const Detail::HydroFiniteRange&){
        // No authoritative F, R, mass, phase stock or frozen coefficient is
        // changed. Retain the complete scalar contractions only in this
        // qualified typed physical path. Open VOF includes its actual
        // exterior momentum/work; legacy/opaque/periodic stay strict.
        using Exact=Detail::HydroFiniteContraction;
        const auto& u=values.normalVelocity();const auto dt=material().deltaTime();
        auto d=std::make_shared<HydroKepConvectionResult::Data>(material(),boundary);
        d->rates.resize(u.size());d->loads.resize(u.size());d->dissipation.resize(u.size());
        std::vector<Exact> rates(u.size());std::array<Exact,2> momentumScale,powerScale;
        std::array<Exact,2> boundaryMomentum,boundaryWork,physicalWork,inletDissipation;
        const auto& segments=material().geometry().segments();
        for(std::size_t id=0;id<segments.size();++id){
            Exact traceFlux;const auto f=material().faceRates()[id];
            if(internal(segments[id]))for(const auto& w:weights()[id])traceFlux.addProduct({f,w.weight,u[w.velocity]});
            else {
                if(data_->vofPrescribed){const auto& segment=segments[id];const int i=owner(segment),a=component(layout(),i);
                    const long double q=(segment.negativeVelocity>=0?1.L:-1.L)*f;
                    const double upwind=q<0?data_->prescribedInlet->traces()[id]:q>0?u[i]:0.;
                    require(boundary.traces()[id]==upwind,"different prescribed-flow upwind boundary trace");
                    Exact load,work,physical,dissipation;load.addProduct({q,upwind});rates[i].add(load);boundaryMomentum[a].add(load);
                    work.addProduct({u[i],q,upwind});work.addProduct({-.5L,q,u[i],u[i]});boundaryWork[a].add(work);
                    physical.addProduct({.5L,q,upwind,upwind});physicalWork[a].add(physical);
                    if(q<0){dissipation.addProduct({-.5L,q,u[i],u[i]});dissipation.addProduct({q,u[i],upwind});
                        dissipation.addProduct({-.5L,q,upwind,upwind});
                        require(dissipation.compare(Exact{})>=0,"negative exact inlet upwind dissipation");inletDissipation[a].add(dissipation);}
                    momentumScale[a].addAbsolute(load);powerScale[a].addAbsolute(work);
                    powerScale[a].addAbsolute(physical);powerScale[a].addAbsolute(dissipation);
                }
                traceFlux.addProduct({f,boundary.traces()[id]});
            }
            // The actual stored diagnostic is F*trace; no early quantized trace.
            d->center.push_back(traceFlux.rounded());
        }
        for(std::size_t i=0;i<u.size();++i){const int a=component(layout(),static_cast<int>(i));
            Exact term,work;term.addProduct({.5L,data_->outflow[i],u[i]});
            work.addProduct({.5L,data_->outflow[i],u[i],u[i]});
            if(data_->vofPrescribed){term.addProduct({-.5L,data_->externalOutflow[i],u[i]});
                work.addProduct({-.5L,data_->externalOutflow[i],u[i],u[i]});}
            rates[i].add(term);
            momentumScale[a].addAbsolute(term);powerScale[a].addAbsolute(work);
        }
        for(const auto& e:data_->finalPairs){const int i=e.first,j=e.second,a=component(layout(),i);
            Exact left,right,leftWork,rightWork;
            left.addProduct({e.coefficient,u[j]});right.addProduct({e.coefficient,u[i]});
            rates[i].add(left);rates[j].add(right,-1);
            leftWork.addProduct({u[i],e.coefficient,u[j]});rightWork.addProduct({u[j],e.coefficient,u[i]});
            momentumScale[a].addAbsolute(left);momentumScale[a].addAbsolute(right);
            powerScale[a].addAbsolute(leftWork);powerScale[a].addAbsolute(rightWork);
        }
        std::array<Exact,2> momentum,power,sourcePower,kinetic;
        for(std::size_t i=0;i<u.size();++i){const int a=component(layout(),static_cast<int>(i));const auto& cell=material().cells()[i];
            d->rates[i]=rates[i].rounded();Exact load;load.add(rates[i],-1);
            if(!conservative)load.addQuotientProduct({u[i],cell.integratedDefect},dt);d->loads[i]=load.rounded();
            // Component work contracts the actual stored row values, not an
            // ideal unrounded operator that could hide a failed physical gate.
            Exact rowPower,rowSource,rowMass,rowKinetic;
            rowPower.addProduct({u[i],d->rates[i]});
            rowSource.addQuotientProduct({conservative?-.5L:.5L,u[i],u[i],cell.integratedDefect},dt);
            rowMass.addQuotientProduct({.5L,u[i],u[i],cell.massChange},dt);
            rowKinetic.addProduct({u[i],d->loads[i]});rowKinetic.add(rowMass,-1);
            momentum[a].add(d->rates[i]);power[a].add(rowPower);sourcePower[a].add(rowSource);kinetic[a].add(rowKinetic);
            Exact rowRate;rowRate.add(d->rates[i]);momentumScale[a].addAbsolute(rowRate);
            powerScale[a].addAbsolute(rowPower);powerScale[a].addAbsolute(rowSource);powerScale[a].addAbsolute(rowSource);
            powerScale[a].addAbsolute(rowMass);powerScale[a].addAbsolute(rowKinetic);
        }
        for(int a=0;a<2;++a){auto& r=d->report.components[a];
            r.momentumRate=momentum[a].rounded();r.convectivePower=power[a].rounded();
            r.materialKineticPower=sourcePower[a].rounded();r.kineticRate=kinetic[a].rounded();
            if(data_->vofPrescribed){r.boundaryMomentumRate=boundaryMomentum[a].rounded();r.boundaryKineticPower=boundaryWork[a].rounded();
                r.boundaryPhysicalKineticPower=physicalWork[a].rounded();r.inletDissipationPower=inletDissipation[a].rounded();
                identity(r.boundaryKineticPower-r.boundaryPhysicalKineticPower-r.inletDissipationPower,
                    powerScale[a].rounded(),"prescribed-flow boundary kinetic decomposition failed");}
            // Original identity formulae and tolerance, with the same gross
            // algebra accumulated without prematurely quantizing its terms.
            r.momentumIdentity=identity(r.momentumRate-r.boundaryMomentumRate,momentumScale[a].rounded(),"compatible momentum balance failed");
            r.kineticIdentity=identity(r.kineticRate-r.materialKineticPower+r.boundaryKineticPower+r.dissipationPower,powerScale[a].rounded(),"nonconservative kinetic identity failed");
        }
        return HydroKepConvectionResult(std::move(d));
    }
}
HydroKepConvectionResult HydroKepConvection::applyPaired(const QuadtreeMesh& mesh,const HydroNativeValues& values,HydroMaterialSources source,LAMP::TimeLayer clock) const
{
    require(data_->topology.paired,"paired apply requires declared port topology");material().validateFor(mesh,source,clock);values.validateFor(layout(),mesh);
    std::vector<double> traces(material().geometry().segments().size());for(const auto& p:periodicPairs()){
        const auto& s=material().geometry().segments();const double trace=.5*values.normalVelocity()[owner(s[p.firstSegment])]+.5*values.normalVelocity()[owner(s[p.secondSegment])];traces[p.firstSegment]=traces[p.secondSegment]=trace;}
    return apply(mesh,values,HydroConvectionBoundary::bind(mesh,material(),traces,source,clock),source,clock);
}
const HydroDualMaterial& HydroKepConvectionResult::material() const noexcept{return data_->material;}
const HydroConvectionBoundary& HydroKepConvectionResult::boundary() const noexcept{return data_->boundary;}
bool HydroKepConvectionResult::hasCenterMomentumRates() const noexcept{return data_->centerAvailable;}
const std::vector<long double>& HydroKepConvectionResult::centerMomentumRates() const noexcept{return data_->center;}
const std::vector<long double>& HydroKepConvectionResult::momentumRates() const noexcept{return data_->rates;}
const std::vector<long double>& HydroKepConvectionResult::momentumLoads() const noexcept{return data_->loads;}
const std::vector<long double>& HydroKepConvectionResult::dissipationRates() const noexcept{return data_->dissipation;}
const std::vector<HydroConvectionDissipationSample>& HydroKepConvectionResult::dissipationSamples() const noexcept{return data_->dissipationSamples;}
bool HydroKepConvectionResult::hasPowerReport() const noexcept{return data_->powerAvailable;}
bool HydroKepConvectionResult::hasComponentReport() const noexcept{return data_->componentAvailable;}
const HydroKepOpenBoundaryBudget* HydroKepConvectionResult::openBoundaryBudget() const noexcept{return data_->boundaryBudget?&*data_->boundaryBudget:nullptr;}
const HydroKepConvectionReport& HydroKepConvectionResult::report() const noexcept{return data_->report;}
const HydroKepOpenApplyReport& HydroKepConvectionResult::openApplyReport() const noexcept{return data_->openApply;}
namespace
{
using OpenExact=Detail::HydroFiniteContraction;

using OpenWall=std::chrono::steady_clock;
// Diagnostic conversion does not perturb the numerical environment or errno.
double openSeconds(OpenWall::duration duration) noexcept
{
    const int error=errno;std::fenv_t native{};const bool saved=std::fegetenv(&native)==0;
#if defined(_M_X64) || defined(__SSE2__)
    const unsigned csr=_mm_getcsr();
#endif
    const double result=std::chrono::duration<double>(duration).count();
    if(saved&&std::fesetenv(&native)!=0)std::terminate();
#if defined(_M_X64) || defined(__SSE2__)
    _mm_setcsr(csr);
#endif
    errno=error;return result;
}
struct OpenApplyPlan
{
    struct Pair {std::uint32_t id=0;int sign=1;};
    std::vector<std::uint32_t> boundaryOffsets,boundaryItems,pairOffsets;
    std::vector<Pair> pairItems;
    // Pure indices of every geometrically supported row in each original
    // fixed chunk, independent of current U or whether its value is zero.
    std::vector<std::uint32_t> boundaryChunkOffsets,boundaryChunkRows,pairChunkOffsets,pairChunkRows;
    std::size_t bytes=0;
};
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
#pragma float_control(precise,on,push)
#pragma float_control(except,on)
#pragma fenv_access(on)
struct OpenFloatingState
{
    std::fenv_t native{};unsigned csr=0;bool valid=false;
    bool capture() noexcept {csr=_mm_getcsr();valid=std::fegetenv(&native)==0;return valid;}
    void restore() const noexcept {if(valid&&std::fesetenv(&native)!=0)std::terminate();_mm_setcsr(csr);}
    bool enterFresh(const OpenFloatingState& from) const noexcept
    {auto clean=from.native;clean._Fe_stat=0;if(std::fesetenv(&clean)!=0)return false;_mm_setcsr(from.csr&~0x3fu);return true;}
};
bool openSameControl(const OpenFloatingState& a,const OpenFloatingState& b) noexcept
{auto x=a.native,y=b.native;x._Fe_stat=0;y._Fe_stat=0;
    return std::memcmp(&x,&y,sizeof x)==0&&(a.csr&~0x3fu)==(b.csr&~0x3fu);}
bool openArithmeticQualified(OpenFloatingState& entry) noexcept
{
#if defined(_M_FP_STRICT) && _M_FP_STRICT && !defined(__FAST_MATH__) && (!defined(_M_FP_FAST) || !_M_FP_FAST)
    if(!entry.capture())return false;unsigned control=0;
    if(_controlfp_s(&control,0,0)!=0||(control&_MCW_EM)!=_MCW_EM||(control&_MCW_RC)!=_RC_NEAR
        ||std::fegetround()!=FE_TONEAREST||(entry.csr&0x1f80u)!=0x1f80u||(entry.csr&0xe040u)!=0)return false;
#ifdef _EM_AMBIGUOUS
    if(control&_EM_AMBIGUOUS)return false;
#endif
    return sizeof(double)==8&&sizeof(long double)==8&&std::numeric_limits<double>::is_iec559
        &&std::numeric_limits<double>::digits==53&&std::numeric_limits<long double>::digits==53;
#else
    (void)entry;return false;
#endif
}
bool openQualified(OpenFloatingState& entry) noexcept
{
#if defined(_M_FP_STRICT) && _M_FP_STRICT && !defined(__FAST_MATH__) && (!defined(_M_FP_FAST) || !_M_FP_FAST)
    if(!entry.capture()||omp_in_parallel()||omp_get_max_threads()<2)return false;
    unsigned control=0;
    if(_controlfp_s(&control,0,0)!=0||(control&_MCW_EM)!=_MCW_EM||(control&_MCW_RC)!=_RC_NEAR
        ||std::fegetround()!=FE_TONEAREST||(entry.csr&0x1f80u)!=0x1f80u||(entry.csr&0xe040u)!=0)return false;
#ifdef _EM_AMBIGUOUS
    if(control&_EM_AMBIGUOUS)return false;
#endif
    return sizeof(double)==8&&sizeof(long double)==8&&std::numeric_limits<double>::is_iec559
        &&std::numeric_limits<double>::digits==53&&std::numeric_limits<long double>::digits==53;
#else
    (void)entry;return false;
#endif
}
struct OpenCleanup
{
    OpenFloatingState finalState;int error=0;bool armed=false;
    void arm() noexcept {error=errno;armed=finalState.capture();errno=error;if(!armed)std::terminate();}
    ~OpenCleanup() noexcept {if(armed){finalState.restore();errno=error;}}
};
struct OpenOptionalFailure {};
// Bounded, original-order slices. No whole-domain array of arbitrary integers.
constexpr std::size_t openBatch=1024,openScratchLimit=64u*1024u*1024u;
// One qualified apply owns this dedicated team. The caller (thread zero)
// drives original serial folds; other threads wait at the next phase boundary.
// This is neither a cross-apply pool nor a nested invocation of the old helper.
class OpenCooperativeTeam
{
    struct PhaseWorker {
        unsigned long nativeFlags=0;
        unsigned csrFlags=0;
        bool failed=false;
    };
    const int requested_;
    int actual_=0;
    bool stop_=false;
    std::size_t count_=0;
    const void* context_=nullptr;
    void (*invoke_)(const void*,std::size_t)=nullptr;
    OpenFloatingState prepared_;
    int parentErrno_=0;
    std::vector<PhaseWorker> workers_;

    // Every thread encounters this SAME orphaned barrier region, including
    // the stop handshake. Master and workers do not use different directives.
    // Current MSVC LLVM OpenMP misbinds the global thread ID if this orphaned
    // barrier is inlined into the outlined worker; keep one real helper call.
    __declspec(noinline) static void synchronize() noexcept
    {
#pragma omp barrier
    }
    void phase(int tid,bool threadReady=true) noexcept
    {
        auto& result=workers_[static_cast<std::size_t>(tid)];
        result.nativeFlags=0;result.csrFlags=0;result.failed=false;
        OpenFloatingState original;const int originalErrno=errno;
        const bool captured=original.capture();
        const bool ready=threadReady&&captured&&original.enterFresh(prepared_);
        if(!ready)result.failed=true;else errno=parentErrno_;
        // Use the ACTUAL team, including a runtime-reduced team. Each index
        // belongs to one contiguous range; no arithmetic term is repartitioned.
        const auto team=static_cast<std::size_t>(actual_),id=static_cast<std::size_t>(tid);
        const auto width=count_/team,remainder=count_%team;
        const auto begin=id*width+std::min(id,remainder);
        const auto end=begin+width+(id<remainder?1u:0u);
        for(auto i=begin;i<end;++i){if(!ready)continue;
            try{invoke_(context_,i);}catch(...){result.failed=true;}}
        if(ready){const int numericalError=errno;OpenFloatingState after;const bool got=after.capture();
            if(!got||numericalError!=parentErrno_||!openSameControl(after,prepared_))result.failed=true;
            if(got){result.nativeFlags=after.native._Fe_stat;result.csrFlags=after.csr&0x3fu;}}
        if(captured)original.restore();errno=originalErrno;
    }
    void serve() noexcept
    {
        OpenFloatingState original;const int originalErrno=errno;const bool captured=original.capture();
        for(;;){
            synchronize();
            if(stop_)break;
            phase(omp_get_thread_num(),captured);
            synchronize();
        }
        if(captured)original.restore();errno=originalErrno;
    }
public:
    static constexpr std::size_t workerStorageBytes() noexcept {return sizeof(PhaseWorker);}
    explicit OpenCooperativeTeam(int requested):requested_(requested),workers_(static_cast<std::size_t>(requested)){}
    template<class Kernel>bool rows(std::size_t count,Kernel&& kernel,unsigned& actualWorkers)
    {
        const int parentErrno=errno;if(!prepared_.capture())return false;
        parentErrno_=parentErrno;count_=count;context_=&kernel;
        using Function=std::remove_reference_t<Kernel>;
        invoke_=[](const void* context,std::size_t index){(*static_cast<const Function*>(context))(index);};
        // kernel/context remain alive until EVERY worker has passed completion.
        synchronize();
        phase(0);
        synchronize();
        bool failed=false;unsigned long nativeFlags=0;unsigned csrFlags=0;
        for(int tid=0;tid<actual_;++tid){const auto& worker=workers_[static_cast<std::size_t>(tid)];
            failed=failed||worker.failed;nativeFlags|=worker.nativeFlags;csrFlags|=worker.csrFlags;}
        if(failed||actual_<2)return false;
        auto merged=prepared_.native;merged._Fe_stat|=nativeFlags;
        if(std::fesetenv(&merged)!=0)std::terminate();_mm_setcsr(prepared_.csr|csrFlags);errno=parentErrno;
        actualWorkers=std::max(actualWorkers,static_cast<unsigned>(actual_));return true;
    }
    template<class Body>auto run(Body&& body)->decltype(body())
    {
        const int callerErrno=errno;OpenFloatingState callerState;
        if(!callerState.capture()){errno=callerErrno;throw OpenOptionalFailure{};}errno=callerErrno;
        decltype(body()) result{};std::exception_ptr failure;
        OpenFloatingState finalState;int finalErrno=callerErrno;bool finalReady=false;
        const int requested=requested_;
#pragma omp parallel num_threads(requested) shared(result,failure,finalState,finalErrno,finalReady,callerState,callerErrno)
        {
            if(omp_get_thread_num()==0){
                actual_=omp_get_num_threads();
                // No exception crosses the OpenMP structured block. Workers
                // always receive exactly one stop at their next START barrier.
                try{
                    // The runtime fork must not become a new parent numerical state.
                    callerState.restore();errno=callerErrno;
                    if(actual_<2||actual_>requested_)throw OpenOptionalFailure{};result=body();}
                catch(...){failure=std::current_exception();}
                finalErrno=errno;finalReady=finalState.capture();errno=finalErrno;
                stop_=true;synchronize();
            }else serve();
        }
        // Shield the caller's numerical state from outer-team exit effects.
        // The original OpenCleanup still restores it after scratch destruction.
        if(finalReady)finalState.restore();errno=finalErrno;
        if(!finalReady)throw OpenOptionalFailure{};
        if(failure)std::rethrow_exception(failure);
        return result;
    }
};
#pragma fenv_access(off)
#pragma float_control(pop)
#endif
// Geometry-only indices are immutable and shared by the two different-U calls.
// Failure remains optional: no new build/acceptance failure or FP side effect.
template<class Segments>std::shared_ptr<const OpenApplyPlan> prepareOpenApplyPlan(const Segments& segments,
    const std::vector<HydroConvectionCoupling>& pairs,std::size_t count)
{
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
    const int error=errno;OpenFloatingState entry;if(!entry.capture()){errno=error;return {};}
    std::shared_ptr<OpenApplyPlan> out;
    try{
        const auto cap=std::numeric_limits<std::uint32_t>::max();
        if(count>=cap||segments.size()>=cap||pairs.size()>cap/2)throw OpenOptionalFailure{};
        const auto boundaryChunks=(segments.size()+openBatch-1)/openBatch,pairChunks=(pairs.size()+openBatch-1)/openBatch;
        const std::size_t bytes=sizeof(OpenApplyPlan)+(count+1)*2*sizeof(std::uint32_t)
            +segments.size()*sizeof(std::uint32_t)+pairs.size()*2*sizeof(OpenApplyPlan::Pair)+count*sizeof(std::uint32_t)
            // At most one touched row per original incidence. Both retained
            // chunk offsets/rows and the largest temporary cursor are charged.
            +(boundaryChunks+pairChunks+2+segments.size()+2*pairs.size()
                +std::max(boundaryChunks,pairChunks))*sizeof(std::uint32_t);
        if(bytes>openScratchLimit/4)throw OpenOptionalFailure{};
        out=std::make_shared<OpenApplyPlan>();out->bytes=bytes;
        out->boundaryOffsets.resize(count+1);out->pairOffsets.resize(count+1);
        for(std::size_t id=0;id<segments.size();++id)if(!internal(segments[id]))++out->boundaryOffsets[owner(segments[id])+1];
        for(const auto& pair:pairs){++out->pairOffsets[pair.first+1];++out->pairOffsets[pair.second+1];}
        for(std::size_t i=0;i<count;++i){out->boundaryOffsets[i+1]+=out->boundaryOffsets[i];out->pairOffsets[i+1]+=out->pairOffsets[i];}
        out->boundaryItems.resize(out->boundaryOffsets.back());out->pairItems.resize(out->pairOffsets.back());
        {std::vector<std::uint32_t> cursor(out->boundaryOffsets.begin(),out->boundaryOffsets.end()-1);
            for(std::size_t id=0;id<segments.size();++id)if(!internal(segments[id]))out->boundaryItems[cursor[owner(segments[id])]++]=static_cast<std::uint32_t>(id);}
        {std::vector<std::uint32_t> cursor(out->pairOffsets.begin(),out->pairOffsets.end()-1);
            for(std::size_t id=0;id<pairs.size();++id){const auto& pair=pairs[id];
                out->pairItems[cursor[pair.first]++]={static_cast<std::uint32_t>(id),1};out->pairItems[cursor[pair.second]++]={static_cast<std::uint32_t>(id),-1};}}
        const auto prepareTouched=[&](const auto& offsets,const auto& items,std::size_t chunks,const auto& itemId,
            auto& chunkOffsets,auto& chunkRows){
            chunkOffsets.resize(chunks+1);
            // Original incidence IDs are ascending within each native row.
            // Visiting rows in native order makes every resulting chunk list
            // sorted and unique without sorting, field tests or numeric work.
            for(std::size_t row=0;row<count;++row){std::size_t previous=chunks;
                for(auto at=offsets[row];at<offsets[row+1];++at){const auto chunk=static_cast<std::size_t>(itemId(items[at]))/openBatch;
                    if(chunk>=chunks||(previous!=chunks&&chunk<previous))throw OpenOptionalFailure{};
                    if(chunk!=previous){++chunkOffsets[chunk+1];previous=chunk;}}}
            for(std::size_t chunk=0;chunk<chunks;++chunk)chunkOffsets[chunk+1]+=chunkOffsets[chunk];
            chunkRows.resize(chunkOffsets.back());
            std::vector<std::uint32_t> cursor(chunkOffsets.begin(),chunkOffsets.end()-1);
            for(std::size_t row=0;row<count;++row){std::size_t previous=chunks;
                for(auto at=offsets[row];at<offsets[row+1];++at){const auto chunk=static_cast<std::size_t>(itemId(items[at]))/openBatch;
                    if(chunk!=previous){chunkRows[cursor[chunk]++]=static_cast<std::uint32_t>(row);previous=chunk;}}}
        };
        prepareTouched(out->boundaryOffsets,out->boundaryItems,boundaryChunks,[](std::uint32_t id){return id;},
            out->boundaryChunkOffsets,out->boundaryChunkRows);
        prepareTouched(out->pairOffsets,out->pairItems,pairChunks,[](const OpenApplyPlan::Pair& item){return item.id;},
            out->pairChunkOffsets,out->pairChunkRows);
    }catch(...){out.reset();entry.restore();errno=error;return {};}
    entry.restore();errno=error;return out;
#else
    (void)segments;(void)pairs;(void)count;return {};
#endif
}

std::shared_ptr<const OpenInteriorProducts> prepareOpenInteriorProducts(const HydroDualMaterial& material,
    const std::shared_ptr<const Detail::HydroOpenConvectionSupport>& support,
    const std::vector<long double>& coefficients,const std::vector<long double>& outflow,std::size_t exteriorBytes) noexcept
{
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
    const int error=errno;OpenFloatingState entry;
    const bool qualified=material.usesConservativeVofMomentum()&&support&&openArithmeticQualified(entry);
    if(!qualified){if(entry.valid)entry.restore();errno=error;return {};}
    std::shared_ptr<const OpenInteriorProducts> result;
    try{
        const auto n=outflow.size();constexpr std::size_t reserve=2u*1024u*1024u;
        const auto cap=std::numeric_limits<std::uint32_t>::max();
        if(n>=cap||support->bytes>openScratchLimit-reserve||exteriorBytes>openScratchLimit-reserve-support->bytes)throw OpenOptionalFailure{};
        const auto available=openScratchLimit-reserve-support->bytes-exteriorBytes;
        if(sizeof(OpenInteriorProducts)>available)throw OpenOptionalFailure{};
        constexpr auto rowBytes=sizeof(OpenProducts::Factor)+sizeof(OpenProducts::Value);
        if(n>(available-sizeof(OpenInteriorProducts))/rowBytes)throw OpenOptionalFailure{};
        const auto fixed=sizeof(OpenInteriorProducts)+n*rowBytes;
        if(coefficients.size()>(available-fixed)/sizeof(OpenProducts::Factor))throw OpenOptionalFailure{};
        auto plan=std::make_shared<OpenInteriorProducts>();plan->diagonal.resize(n);plan->pairs.resize(coefficients.size());
        // A row without exterior segments has exactly zero externalOutflow:
        // that array is value-initialized and written only by exterior owners.
        for(std::size_t i=0;i<n;++i)if(support->boundaryOffsets[i]==support->boundaryOffsets[i+1])
            plan->diagonal[i]=OpenProducts::compile(.5,static_cast<double>(outflow[i]));
        for(std::size_t i=0;i<coefficients.size();++i)if(coefficients[i]!=0)
            plan->pairs[i]=OpenProducts::compile(static_cast<double>(coefficients[i]),1.);
        plan->bytes=sizeof(*plan)+plan->diagonal.capacity()*sizeof(OpenProducts::Factor)+plan->pairs.capacity()*sizeof(OpenProducts::Factor);
        if(plan->bytes>available||n>(available-plan->bytes)/sizeof(OpenProducts::Value))throw OpenOptionalFailure{};
        result=std::move(plan);
    }catch(...){result.reset();}
    entry.restore();errno=error;return result;
#else
    (void)material;(void)support;(void)coefficients;(void)outflow;(void)exteriorBytes;return {};
#endif
}

// The public LocalPlic restriction validates the live mesh/source on every
// call. Batch only this immutable checker, preserving every per-port exact
// fold. Any optional error replays the untouched original loop/first error.
void certifyVofFreeOutflowMeans(const QuadtreeMesh& mesh,const HydroDualMaterial& material,
    const std::vector<int>& freeOutflow,HydroConvectionMomentReport& report,
    const CheckedVofMaterialBatch* checked)
{
    const auto& volume=validateVofNativeSource(mesh,material,checked);const auto& slab=volume.advector();
    const auto& layout=material.layout();const auto& reconstruction=material.material().fluxReconstruction();
    const auto literal=[&]{
    for(int u:freeOutflow){const int row=layout.nativeFace(u);const auto face=layout.domain().face(mesh,row);
        std::vector<Interface::GeometricVofNativeSubface> queries;
        for(int id:reconstruction.primalAtoms()[u]){const auto& atom=reconstruction.atoms()[id];
            queries.push_back({row,{atom.begin.x,atom.begin.y},{atom.end.x,atom.end.y}});}
        require(!queries.empty(),"VOF free outflow has no actual source atoms");
        const auto atoms=volume.restrictNativeFacesFull(mesh,queries);Detail::HydroFiniteContraction sum,scale,expected;
        for(const auto& atom:atoms){sum.add(atom.full.total);Detail::HydroFiniteContraction value;value.add(atom.full.total);scale.addAbsolute(value);}
        const auto total=volume.faces()[row].full.total;Detail::HydroFiniteContraction native;native.add(total);scale.addAbsolute(native);
        const auto partition=identity(sum.rounded()-total,scale.rounded(),"VOF free-outflow total differs from original native accepted flux");
        expected.addProduct({.5L,face.area,material.deltaTime(),slab.beginNormal()[row]});
        expected.addProduct({.5L,face.area,material.deltaTime(),slab.endNormal()[row]});
        const auto mean=identity(total-expected.rounded(),std::abs(total)+std::abs(expected.rounded()),"VOF free-outflow flux differs from actual slab native mean");
        report.maxNativeMeanDefect=std::max({report.maxNativeMeanDefect,partition,mean});
    }
    };
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
    const int entryErrno=errno;OpenFloatingState entry;
    const bool qualified=freeOutflow.size()>1&&openQualified(entry);
    if(entry.valid)entry.restore();errno=entryErrno;
    if(qualified){
        try{
            OpenCleanup cleanup;
            struct Port {int row=0;long double area=0;std::size_t begin=0,end=0;};
            using Query=Interface::GeometricVofNativeSubface;
            using Face=Interface::GeometricVofNativeFullFlux;
            // The payload plus a fixed allowance for this serial exact fold
            // is bounded by the existing optional scratch cap. No U/cache.
            constexpr std::size_t allowance=64u*1024u,itemBytes=sizeof(Query)+sizeof(Face);
            if(freeOutflow.size()>(openScratchLimit-allowance)/sizeof(Port))throw OpenOptionalFailure{};
            const auto fixed=allowance+freeOutflow.size()*sizeof(Port);
            const auto queryLimit=(openScratchLimit-fixed)/itemBytes;std::size_t queryCount=0;
            for(int u:freeOutflow){const auto count=reconstruction.primalAtoms()[u].size();
                require(count>0,"VOF free outflow has no actual source atoms");
                if(count>queryLimit-queryCount)throw OpenOptionalFailure{};queryCount+=count;}
            std::vector<Query> queries;queries.reserve(queryCount);
            std::vector<Port> ports;ports.reserve(freeOutflow.size());
            for(int u:freeOutflow){const int row=layout.nativeFace(u);const auto face=layout.domain().face(mesh,row);
                const auto begin=queries.size();
                for(int id:reconstruction.primalAtoms()[u]){const auto& atom=reconstruction.atoms()[id];
                    queries.push_back({row,{atom.begin.x,atom.begin.y},{atom.end.x,atom.end.y}});}
                ports.push_back({row,face.area,begin,queries.size()});}
            const auto atoms=volume.restrictNativeFacesFull(mesh,queries);
            if(atoms.size()!=queryCount)throw OpenOptionalFailure{};
            auto maximum=report.maxNativeMeanDefect;
            for(const auto& port:ports){Detail::HydroFiniteContraction sum,scale,expected;
                for(auto k=port.begin;k<port.end;++k){const auto& atom=atoms[k];sum.add(atom.full.total);
                    Detail::HydroFiniteContraction value;value.add(atom.full.total);scale.addAbsolute(value);}
                const auto total=volume.faces()[port.row].full.total;Detail::HydroFiniteContraction native;native.add(total);scale.addAbsolute(native);
                const auto partition=identity(sum.rounded()-total,scale.rounded(),"VOF free-outflow total differs from original native accepted flux");
                expected.addProduct({.5L,port.area,material.deltaTime(),slab.beginNormal()[port.row]});
                expected.addProduct({.5L,port.area,material.deltaTime(),slab.endNormal()[port.row]});
                const auto mean=identity(total-expected.rounded(),std::abs(total)+std::abs(expected.rounded()),"VOF free-outflow flux differs from actual slab native mean");
                maximum=std::max({maximum,partition,mean});}
            const int finalErrno=errno;OpenFloatingState finalState;
            if(!finalState.capture()||!openSameControl(entry,finalState)||finalErrno!=entryErrno)throw OpenOptionalFailure{};
            errno=finalErrno;report.maxNativeMeanDefect=maximum;cleanup.arm();return;
        }catch(...){entry.restore();errno=entryErrno;}
    }
#endif
    literal();
}

template<class... Factors> void addOpenMass(OpenExact& value,const HydroVofBoundaryMomentumFlux& flux,
    const HydroDualMaterial& material,bool incoming,Factors... factors)
{
    const auto& q=incoming?flux.incomingPhaseIntegral:flux.outgoingPhaseIntegral;
    value.addQuotientProduct({material.material().liquidDensity(),q[0],static_cast<long double>(factors)...},material.deltaTime());
    value.addQuotientProduct({material.material().gasDensity(),q[1],static_cast<long double>(factors)...},material.deltaTime());
}
}
HydroConvectionBoundary HydroConvectionBoundary::bindVofOpenMomentum(const QuadtreeMesh& mesh,HydroDualMaterial m,
    const std::vector<Vector2D>& begin,const std::vector<Vector2D>& end,const std::vector<int>& freePorts,
    const std::vector<Vector2D>& ambient,HydroMaterialSources sources,LAMP::TimeLayer clock,HydroVofPressureOutletMomentumModel model)
{
    require(model==HydroVofPressureOutletMomentumModel::Ambient||model==HydroVofPressureOutletMomentumModel::InteriorContinuationZeroShear,"unknown open momentum model");
    const auto checked=CheckedVofMaterialBatch::validate(mesh,m,sources,clock);
    const auto& volume=validateVofNativeSource(mesh,m,&checked);
    const auto& slab=volume.advector();const auto* law=slab.boundary();
    require(law&&!law->allClosed(),"open momentum requires actual F1 boundary");
    validateVolumeBoundaryAfterMaterial(mesh,m,*law,volume.source(),&checked);
    const auto& layout=m.layout();const auto& segments=m.geometry().segments();
    require(begin.size()==segments.size()&&end.size()==segments.size(),"open momentum endpoint segment dimensions");
    require(ambient.size()==static_cast<std::size_t>(layout.domain().nativeFaceCount()),"open momentum ambient native-row dimensions");
    const auto zero=[](Vector2D q){return q.x==0&&q.y==0;};const auto finiteVector=[](Vector2D q){return std::isfinite(q.x)&&std::isfinite(q.y);};
    std::vector<bool> free(layout.velocityCount()),freeRows(ambient.size());
    using Mode=Interface::GeometricVofAdvectorSlab::Boundary::Mode;
    for(int u:freePorts){require(u>=0&&u<layout.velocityCount()&&!free[u],"duplicate or invalid ambient free pressure port");
        const int row=layout.nativeFace(u);const auto face=layout.domain().face(mesh,row);
        require(face.outwardSign!=0&&face.area>0,"ambient pressure port is internal or axis");
        require(law->rules()[row]==Mode::PureLiquidReservoir||law->rules()[row]==Mode::PureGasReservoir,
            "ambient pressure port requires declared F1 pure reservoir");free[u]=true;freeRows[row]=true;}
    require(!freePorts.empty(),"explicit open momentum requires physical free pressure ports");
    for(std::size_t row=0;row<ambient.size();++row){require(finiteVector(ambient[row]),"nonfinite ambient momentum");
        if(!freeRows[row]||model==HydroVofPressureOutletMomentumModel::InteriorContinuationZeroShear)
            require(zero(ambient[row]),"ambient momentum supplied outside ambient pressure ports");}
    std::vector<HydroVofBoundaryMomentumFlux> fluxes(segments.size());
    std::vector<Interface::GeometricVofNativeSubface> queries;std::vector<std::size_t> ids;
    for(std::size_t id=0;id<segments.size();++id){const auto& g=segments[id];
        require(finiteVector(begin[id])&&finiteVector(end[id]),"nonfinite prescribed open momentum endpoint");
        if(internal(g)){require(zero(begin[id])&&zero(end[id]),"nonzero unused internal open momentum");continue;}
        require(g.primalVelocity>=0&&g.primalVelocity<layout.velocityCount(),"open momentum segment lacks native port");
        const int row=layout.nativeFace(g.primalVelocity);const auto face=layout.domain().face(mesh,row);
        const int sign=g.negativeVelocity>=0?1:-1;require(face.outwardSign==sign&&(face.horizontal?1:0)==g.normalAxis,
            "open momentum segment/native exterior orientation differs");
        auto& f=fluxes[id];f.nativeRow=row;f.velocityDof=owner(g);
        f.canonicalOutwardMassRate=sign*m.faceRates()[id];
        require(std::isfinite(f.canonicalOutwardMassRate),"nonfinite canonical open momentum mass rate");
        if(face.area==0){require(law->rules()[row]==Mode::Closed&&slab.beginNormal()[row]==0&&slab.endNormal()[row]==0
                &&m.faceRates()[id]==0&&m.faceIntegrals()[id]==0&&zero(begin[id])&&zero(end[id]),"nonzero actual axis open momentum");continue;}
        queries.push_back({row,{g.begin.x,g.begin.y},{g.end.x,g.end.y}});ids.push_back(id);
    }
    const auto actual=volume.restrictNativeFacesFull(mesh,queries);require(actual.size()==ids.size(),"open momentum actual gross query count");
    for(std::size_t k=0;k<ids.size();++k){const auto id=ids[k];const auto& g=segments[id];auto& f=fluxes[id];
        const auto& a=actual[k].full;const int sign=g.negativeVelocity>=0?1:-1;
        f.outgoingPhaseIntegral=sign>0?std::array<long double,2>{a.positiveLiquid,a.positiveGas}:std::array<long double,2>{a.negativeLiquid,a.negativeGas};
        f.incomingPhaseIntegral=sign>0?std::array<long double,2>{a.negativeLiquid,a.negativeGas}:std::array<long double,2>{a.positiveLiquid,a.positiveGas};
        for(auto q:f.outgoingPhaseIntegral)require(std::isfinite(q)&&q>=0,"invalid actual outgoing phase gross");
        for(auto q:f.incomingPhaseIntegral)require(std::isfinite(q)&&q>=0,"invalid actual incoming phase gross");
        OpenExact outgoing,incoming;addOpenMass(outgoing,f,m,false,1);addOpenMass(incoming,f,m,true,1);
        f.outgoingMassRate=outgoing.rounded();f.incomingMassRate=incoming.rounded();
        const bool hasIncoming=f.incomingPhaseIntegral[0]>0||f.incomingPhaseIntegral[1]>0;
        if(free[g.primalVelocity]){require(zero(begin[id])&&zero(end[id]),"free pressure ambient must not prescribe endpoint normal");
            f.incomingModel=model;const auto b=ambient[f.nativeRow];
            f.incomingMomentum=component(layout,f.velocityDof)==0?b.x:b.y;}
        else if(hasIncoming){require(begin[id].x==end[id].x&&begin[id].y==end[id].y,"open momentum first scope requires steady prescribed vectors");
            require((g.normalAxis==0?begin[id].x:begin[id].y)==slab.beginNormal()[f.nativeRow]
                &&(g.normalAxis==0?end[id].x:end[id].y)==slab.endNormal()[f.nativeRow],"prescribed open momentum normal differs from actual slab");
            const auto mode=law->rules()[f.nativeRow];require(mode==Mode::PureLiquidReservoir||mode==Mode::PureGasReservoir,"incoming open momentum lacks pure reservoir");
            f.incomingMomentum=component(layout,f.velocityDof)==0?begin[id].x:begin[id].y;
        }else require(zero(begin[id])&&zero(end[id]),"nonzero unused outgoing open momentum vector");
        if(law->rules()[f.nativeRow]==Mode::Closed)require(outgoing.compare(OpenExact{})==0&&incoming.compare(OpenExact{})==0,
            "nonzero phase gross through Closed open momentum wall");
        OpenExact difference,scale;difference.add(f.canonicalOutwardMassRate);difference.add(outgoing,-1);difference.add(incoming);
        scale.addAbsolute(outgoing);scale.addAbsolute(incoming);OpenExact canonical;canonical.add(f.canonicalOutwardMassRate);scale.addAbsolute(canonical);
        identity(difference.rounded(),scale.rounded(),"canonical open momentum mass rate differs from actual phase gross budget");
    }
    std::vector<double> traces;traces.reserve(fluxes.size());for(const auto& f:fluxes)traces.push_back(f.incomingMomentum);
    auto data=std::make_shared<Data>(std::move(m),std::move(traces));data->vofEndpointVelocity=std::array<std::vector<Vector2D>,2>{{begin,end}};
    data->vofOpenMomentum=std::move(fluxes);data->freePressurePorts=freePorts;return HydroConvectionBoundary(std::move(data));
}
HydroConvectionBoundary HydroConvectionBoundary::bindVofOpenProductionMomentum(const QuadtreeMesh& mesh,HydroDualMaterial m,
    const std::vector<int>& freePorts,const std::vector<Vector2D>& ambient,
    const std::function<Vector2D(int,LAMP::TimeLayer)>& endpointVelocity,
    HydroMaterialSources sources,LAMP::TimeLayer clock,HydroVofPressureOutletMomentumModel model)
{
    require(model==HydroVofPressureOutletMomentumModel::Ambient||model==HydroVofPressureOutletMomentumModel::InteriorContinuationZeroShear,"unknown open momentum model");
    const auto actualSources=m.sources();const auto actualClock=m.timeLayer();
    // Dual::prepare publishes only immutable material/geometry after its actual
    // construction gates. Match every supplied tag before deferring the LIVE
    // mesh scan to the original public Full restriction below.
    bool deferFaceChecks=m.material().model()==HydroPhaseGeometry::Model::GeometricVof
        &&m.material().vofTransport().source().geometryKind()==Interface::GeometricVofGeometryKind::LocalPlic
        &&sources.model==actualSources.model&&sources.begin==actualSources.begin&&sources.midpoint==actualSources.midpoint
        &&sources.end==actualSources.end&&sources.advector==actualSources.advector&&sources.vofStamp&&actualSources.vofStamp
        &&sources.vofStamp->sameInstance(m.material().vofTransport())&&sources.vofIds==actualSources.vofIds
        &&clock.step==actualClock.step&&clock.begin==actualClock.begin&&clock.end==actualClock.end
        &&clock.position==actualClock.position&&wallExecutionQualified();
    std::optional<CheckedVofMaterialBatch> checked;
    if(!deferFaceChecks)checked=CheckedVofMaterialBatch::validate(mesh,m,sources,clock);
    const auto& volume=deferFaceChecks?m.material().vofTransport():validateVofNativeSource(mesh,m,&*checked);
    const auto& slab=volume.advector();const auto* law=slab.boundary();
    // Bounded invalid inputs keep the original live validator/error path and
    // actual face checks; the immutable preflight is never a new rejection.
    const auto validateOriginal=[&]{
        if(!checked)checked=CheckedVofMaterialBatch::validate(mesh,m,sources,clock);
        (void)validateVofNativeSource(mesh,m,&*checked);
        if(law)validateVolumeBoundaryAfterMaterial(mesh,m,*law,volume.source(),&*checked);
        deferFaceChecks=false;
    };
    const auto requireInput=[&](bool ok,const char* why){if(!ok){validateOriginal();require(false,why);}};
    requireInput(law&&!law->allClosed(),"open momentum requires actual F1 boundary");
    if(checked)validateVolumeBoundaryAfterMaterial(mesh,m,*law,volume.source(),&*checked);
    const auto& layout=m.layout();const auto& segments=m.geometry().segments();
    const auto& faceMeasures=m.geometry().metrics().faceMeasures();
    requireInput(ambient.size()==static_cast<std::size_t>(layout.domain().nativeFaceCount()),"open momentum ambient native-row dimensions");
    const auto zero=[](Vector2D q){return q.x==0&&q.y==0;};const auto finiteVector=[](Vector2D q){return std::isfinite(q.x)&&std::isfinite(q.y);};
    // Snapshot the actual ports before arbitrary user callbacks can mutate an
    // externally captured vector. Needed ambient components are copied below.
    auto ownedPorts=freePorts;std::vector<bool> free(layout.velocityCount()),freeRows(ambient.size());
    using Mode=Interface::GeometricVofAdvectorSlab::Boundary::Mode;
    for(int u:ownedPorts){requireInput(u>=0&&u<layout.velocityCount()&&!free[u],"duplicate or invalid ambient free pressure port");
        const int row=layout.nativeFace(u);
        if(deferFaceChecks){
            // acrossPrimalFaces always appends the exterior native-u segment
            // to this same volume's sides. Otherwise use the actual slow face.
            bool exterior=false;
            for(const auto& side:m.geometry().volumes()[u].sides)for(int id:side){const auto& g=segments[id];
                if(!internal(g)&&g.primalVelocity==u)exterior=true;}
            if(!exterior||!(faceMeasures[u]>0))validateOriginal();
        }
        if(!deferFaceChecks){const auto face=layout.domain().face(mesh,row);
            require(face.outwardSign!=0&&face.area>0,"ambient pressure port is internal or axis");}
        requireInput(law->rules()[row]==Mode::PureLiquidReservoir||law->rules()[row]==Mode::PureGasReservoir,
            "ambient pressure port requires declared F1 pure reservoir");free[u]=true;freeRows[row]=true;}
    requireInput(!ownedPorts.empty(),"explicit open momentum requires physical free pressure ports");
    for(std::size_t row=0;row<ambient.size();++row){requireInput(finiteVector(ambient[row]),"nonfinite ambient momentum");
        if(!freeRows[row]||model==HydroVofPressureOutletMomentumModel::InteriorContinuationZeroShear)
            requireInput(zero(ambient[row]),"ambient momentum supplied outside ambient pressure ports");}
    std::vector<HydroVofBoundaryMomentumFlux> fluxes(segments.size());
    using Query=Interface::GeometricVofNativeSubface;
    using Key=std::tuple<int,std::uint64_t,std::uint64_t,std::uint64_t,std::uint64_t>;
    std::vector<Query> queries;std::map<Key,std::size_t> querySlots;
    // Every query endpoint here is an exact extension of a native Vector2D double.
    const auto bits=[](long double value){const double native=static_cast<double>(value);std::uint64_t word=0;std::memcpy(&word,&native,sizeof word);return word;};
    const auto append=[&](const Query& q){const Key key{q.nativeRow,bits(q.begin.x),bits(q.begin.y),bits(q.end.x),bits(q.end.y)};
        const auto found=querySlots.find(key);if(found!=querySlots.end())return found->second;
        const auto index=queries.size();queries.push_back(q);querySlots.emplace(key,index);return index;};
    const auto closedZero=[&](int row,std::size_t id){
        if(!deferFaceChecks||row<0||static_cast<std::size_t>(row)>=volume.faces().size()
            ||law->rules()[row]!=Mode::Closed||slab.beginNormal()[row]!=0||slab.endNormal()[row]!=0
            ||m.faceRates()[id]!=0||m.faceIntegrals()[id]!=0)return false;
        const auto& parent=volume.faces()[row];if(parent.nativeRow!=row||(parent.negativeCell<0&&parent.positiveCell<0))return false;
        const auto& a=parent.full;
        return a.total==0&&a.positiveTotal==0&&a.negativeTotal==0&&a.rawLiquid==0&&a.liquid==0&&a.gas==0&&a.mass==0
            &&a.positiveLiquid==0&&a.negativeLiquid==0&&a.positiveGas==0&&a.negativeGas==0
            &&a.exhaustionLiquidCorrection==0&&a.exhaustionGasCorrection==0&&a.limiter==1;
    };
    std::vector<std::size_t> ids,dualQueries,exteriorIds;
    for(std::size_t id=0;id<segments.size();++id){const auto& g=segments[id];if(internal(g))continue;
        exteriorIds.push_back(id);
        requireInput(g.primalVelocity>=0&&g.primalVelocity<layout.velocityCount(),"open momentum segment lacks native port");
        const int row=layout.nativeFace(g.primalVelocity);const int sign=g.negativeVelocity>=0?1:-1;
        bool zeroArea=faceMeasures[g.primalVelocity]==0;
        if(!deferFaceChecks){const auto face=layout.domain().face(mesh,row);
            require(face.outwardSign==sign&&(face.horizontal?1:0)==g.normalAxis,"open momentum segment/native exterior orientation differs");zeroArea=face.area==0;}
        auto& f=fluxes[id];f.nativeRow=row;f.velocityDof=owner(g);f.canonicalOutwardMassRate=sign*m.faceRates()[id];
        require(std::isfinite(f.canonicalOutwardMassRate),"nonfinite canonical open momentum mass rate");
        if(free[g.primalVelocity]){f.incomingModel=model;const auto b=ambient[row];f.incomingMomentum=component(layout,f.velocityDof)==0?b.x:b.y;}
        if(zeroArea){requireInput(law->rules()[row]==Mode::Closed&&slab.beginNormal()[row]==0&&slab.endNormal()[row]==0
            &&m.faceRates()[id]==0&&m.faceIntegrals()[id]==0,"nonzero actual axis open momentum");continue;}
        // Accepted LocalPlic P0 restrictions of this entire zero measure are
        // zero. Metadata/canonical signed zero are retained above; true live
        // face checks below still cover every exterior segment.
        if(closedZero(row,id))continue;
        ids.push_back(id);dualQueries.push_back(append({row,{g.begin.x,g.begin.y},{g.end.x,g.end.y}}));
    }
    struct Port {int row=0;long double area=0;std::size_t begin=0,end=0;};
    std::vector<Port> ports;std::vector<std::size_t> atomQueries;const auto& reconstruction=m.material().fluxReconstruction();
    for(int u:ownedPorts){const int row=layout.nativeFace(u);const auto begin=atomQueries.size();
        for(int id:reconstruction.primalAtoms()[u]){const auto& atom=reconstruction.atoms()[id];
            atomQueries.push_back(append({row,{atom.begin.x,atom.begin.y},{atom.end.x,atom.end.y}}));}
        requireInput(atomQueries.size()>begin,"VOF free outflow has no actual source atoms");
        const auto area=deferFaceChecks?0.L:static_cast<long double>(layout.domain().face(mesh,row).area);
        ports.push_back({row,area,begin,atomQueries.size()});}
    // This original public operation performs the first live mesh scan on the
    // qualified route, before it reads any actual face for the immutable union.
    const auto actual=volume.restrictNativeFacesFull(mesh,queries);require(actual.size()==queries.size(),"open momentum actual gross query count");
    if(deferFaceChecks){
        checked=CheckedVofMaterialBatch(m);
        (void)validateVofNativeSource(mesh,m,&*checked); // Includes BOTH selected-domain index arrays.
        validateVolumeBoundaryAfterMaterial(mesh,m,*law,volume.source(),&*checked);
        for(auto& port:ports){const auto face=layout.domain().face(mesh,port.row);
            require(face.outwardSign!=0&&face.area>0,"ambient pressure port is internal or axis");port.area=face.area;}
        for(auto id:exteriorIds){const auto& g=segments[id];
            const int row=layout.nativeFace(g.primalVelocity);const auto face=layout.domain().face(mesh,row);const int sign=g.negativeVelocity>=0?1:-1;
            require(face.outwardSign==sign&&(face.horizontal?1:0)==g.normalAxis,"open momentum segment/native exterior orientation differs");
            if(face.area==0)require(law->rules()[row]==Mode::Closed&&slab.beginNormal()[row]==0&&slab.endNormal()[row]==0
                &&m.faceRates()[id]==0&&m.faceIntegrals()[id]==0,"nonzero actual axis open momentum");
        }
    }
    for(std::size_t k=0;k<ids.size();++k){const auto id=ids[k];const auto& g=segments[id];auto& f=fluxes[id];const auto& a=actual[dualQueries[k]].full;
        const int sign=g.negativeVelocity>=0?1:-1;
        f.outgoingPhaseIntegral=sign>0?std::array<long double,2>{a.positiveLiquid,a.positiveGas}:std::array<long double,2>{a.negativeLiquid,a.negativeGas};
        f.incomingPhaseIntegral=sign>0?std::array<long double,2>{a.negativeLiquid,a.negativeGas}:std::array<long double,2>{a.positiveLiquid,a.positiveGas};
        for(auto q:f.outgoingPhaseIntegral)require(std::isfinite(q)&&q>=0,"invalid actual outgoing phase gross");
        for(auto q:f.incomingPhaseIntegral)require(std::isfinite(q)&&q>=0,"invalid actual incoming phase gross");
    }
    // No optional replay encloses these callbacks. Preserve the old Step's
    // first-segment/native-row source->End order and use the actual clocks.
    std::vector<Vector2D> first(segments.size()),last(segments.size()),sourceVector(ambient.size()),endVector(ambient.size());
    std::vector<bool> queried(ambient.size());bool invokedEndpoint=false;
    for(auto id:ids){const auto& g=segments[id];const auto& f=fluxes[id];if(free[g.primalVelocity])continue;
        if(!(f.incomingPhaseIntegral[0]>0||f.incomingPhaseIntegral[1]>0))continue;const int row=f.nativeRow;
        if(!queried[row]){require(bool(endpointVelocity),"missing actual VOF endpoint momentum query");
            invokedEndpoint=true;sourceVector[row]=endpointVelocity(row,volume.sourceTime());endVector[row]=endpointVelocity(row,volume.endTime());queried[row]=true;}
        first[id]=sourceVector[row];last[id]=endVector[row];
    }
    // Content ownership cannot authenticate a mesh or FP mode changed by a callback.
    if(invokedEndpoint)layout.validateMesh(mesh);
    require(wallExecutionQualified(),"production open momentum callback changed floating execution mode");
    const auto expected=m.sources();const auto expectedClock=m.timeLayer();
    require(sources.model==expected.model&&sources.begin==expected.begin&&sources.midpoint==expected.midpoint&&sources.end==expected.end
        &&sources.advector==expected.advector&&sources.vofStamp&&expected.vofStamp&&sources.vofStamp->sameInstance(*expected.vofStamp)
        &&sources.vofIds==expected.vofIds&&clock.step==expectedClock.step&&clock.begin==expectedClock.begin&&clock.end==expectedClock.end
        &&clock.position==expectedClock.position,"different actual production open momentum source/clock");
    for(auto id:ids){const auto& g=segments[id];auto& f=fluxes[id];
        require(finiteVector(first[id])&&finiteVector(last[id]),"nonfinite prescribed open momentum endpoint");
        OpenExact outgoing,incoming;addOpenMass(outgoing,f,m,false,1);addOpenMass(incoming,f,m,true,1);
        f.outgoingMassRate=outgoing.rounded();f.incomingMassRate=incoming.rounded();
        const bool hasIncoming=f.incomingPhaseIntegral[0]>0||f.incomingPhaseIntegral[1]>0;
        if(!free[g.primalVelocity]&&hasIncoming){require(first[id].x==last[id].x&&first[id].y==last[id].y,"open momentum first scope requires steady prescribed vectors");
            require((g.normalAxis==0?first[id].x:first[id].y)==slab.beginNormal()[f.nativeRow]
                &&(g.normalAxis==0?last[id].x:last[id].y)==slab.endNormal()[f.nativeRow],"prescribed open momentum normal differs from actual slab");
            const auto mode=law->rules()[f.nativeRow];require(mode==Mode::PureLiquidReservoir||mode==Mode::PureGasReservoir,"incoming open momentum lacks pure reservoir");
            f.incomingMomentum=component(layout,f.velocityDof)==0?first[id].x:first[id].y;}
        if(law->rules()[f.nativeRow]==Mode::Closed)require(outgoing.compare(OpenExact{})==0&&incoming.compare(OpenExact{})==0,"nonzero phase gross through Closed open momentum wall");
        OpenExact difference,scale;difference.add(f.canonicalOutwardMassRate);difference.add(outgoing,-1);difference.add(incoming);
        scale.addAbsolute(outgoing);scale.addAbsolute(incoming);OpenExact canonical;canonical.add(f.canonicalOutwardMassRate);scale.addAbsolute(canonical);
        identity(difference.rounded(),scale.rounded(),"canonical open momentum mass rate differs from actual phase gross budget");
    }
    Data::ProductionReceipt receipt;
    for(const auto& port:ports){OpenExact sum,scale,expectedIntegral;
        for(auto k=port.begin;k<port.end;++k){const auto total=actual[atomQueries[k]].full.total;sum.add(total);OpenExact term;term.add(total);scale.addAbsolute(term);}
        const auto total=volume.faces()[port.row].full.total;OpenExact native;native.add(total);scale.addAbsolute(native);
        const auto partition=identity(sum.rounded()-total,scale.rounded(),"VOF free-outflow total differs from original native accepted flux");
        expectedIntegral.addProduct({.5L,port.area,m.deltaTime(),slab.beginNormal()[port.row]});
        expectedIntegral.addProduct({.5L,port.area,m.deltaTime(),slab.endNormal()[port.row]});
        const auto mean=identity(total-expectedIntegral.rounded(),std::abs(total)+std::abs(expectedIntegral.rounded()),"VOF free-outflow flux differs from actual slab native mean");
        receipt.maximumNativeMeanDefect=std::max({receipt.maximumNativeMeanDefect,partition,mean});
    }
    for(int u=0;u<layout.velocityCount();++u){const int row=layout.nativeFace(u);const auto face=layout.domain().face(mesh,row);
        if(face.negativeCell>=0&&face.positiveCell>=0)continue;if(law->rules()[row]==Mode::Closed){receipt.closedVelocityDofs.push_back(u);
            require(slab.beginNormal()[row]==0&&slab.endNormal()[row]==0,"nonzero actual F1 Closed endpoint normal");}}
    // Reuse the original global segment order without another table or sort.
    // Leave the existing16MiB geometry-support ceiling and2MiB workspace
    // reserve available inside the SAME64MiB optional budget.
    constexpr std::size_t supportCeiling=16u*1024u*1024u,reserve=2u*1024u*1024u;
    receipt.exteriorBytes=sizeof(receipt.exteriorSegments)+sizeof(receipt.exteriorBytes);
    if(exteriorIds.capacity()<=(openScratchLimit-supportCeiling-reserve-receipt.exteriorBytes)/sizeof(std::size_t)){
        receipt.exteriorSegments=std::move(exteriorIds);receipt.exteriorBytes+=receipt.exteriorSegments.capacity()*sizeof(std::size_t);}
    std::vector<double> traces;traces.reserve(fluxes.size());for(const auto& f:fluxes)traces.push_back(f.incomingMomentum);
    auto data=std::make_shared<Data>(std::move(m),std::move(traces));
    data->vofEndpointVelocity=std::array<std::vector<Vector2D>,2>{{std::move(first),std::move(last)}};
    data->vofOpenMomentum=std::move(fluxes);data->freePressurePorts=std::move(ownedPorts);data->productionReceipt=std::move(receipt);
    return HydroConvectionBoundary(std::move(data));
}
bool HydroConvectionBoundary::hasVofOpenMomentum() const noexcept{return data_&&data_->vofOpenMomentum.has_value();}
const std::vector<HydroVofBoundaryMomentumFlux>& HydroConvectionBoundary::vofOpenMomentumFluxes() const
{require(hasVofOpenMomentum(),"no typed open momentum receipt");return *data_->vofOpenMomentum;}
const std::vector<int>& HydroConvectionBoundary::freePressureOutletVelocityDofs() const
{require(hasVofOpenMomentum(),"no typed open momentum pressure ports");return data_->freePressurePorts;}
HydroKepConvection HydroKepConvection::buildVofOpenFlow(const QuadtreeMesh& mesh,HydroDualMaterial m,
    const HydroAdvectorBoundary& walls,const HydroConvectionBoundary& incoming)
{return buildVofOpenImpl(mesh,std::move(m),walls,incoming,true);}
HydroKepConvection HydroKepConvection::buildVofOpenProductionFlow(const QuadtreeMesh& mesh,HydroDualMaterial m,
    const HydroAdvectorBoundary& walls,const HydroConvectionBoundary& incoming)
{return buildVofOpenImpl(mesh,std::move(m),walls,incoming,false);}
HydroKepConvection HydroKepConvection::buildVofOpenImpl(const QuadtreeMesh& mesh,HydroDualMaterial m,
    const HydroAdvectorBoundary& walls,const HydroConvectionBoundary& incoming,bool includeCenterTables)
{
    const bool prepared=!includeCenterTables&&incoming.hasVofOpenMomentum()&&incoming.data_->productionReceipt
        &&&incoming.material().cells()==&m.cells()&&wallExecutionQualified();
    if(prepared)m.layout().validateMesh(mesh);
    const auto checked=prepared?CheckedVofMaterialBatch(m):CheckedVofMaterialBatch::validate(mesh,m,m.sources(),m.timeLayer());
    validateConvectionBoundaryAfterMaterial(mesh,m,incoming,m.sources(),m.timeLayer(),&checked);
    require(incoming.hasVofOpenMomentum(),"open KEP requires explicit typed ambient momentum");
    if(prepared){validateWallsAfterMaterial(mesh,m,walls,&checked);
        require(incoming.data_->productionReceipt->closedVelocityDofs==walls.velocityDofs(),"VOF Closed wall ports differ from actual F1 boundary owner");}
    else validateVofOpenSource(mesh,m,walls,&checked);
    const auto& l=m.layout();std::vector<bool> free(l.velocityCount());
    for(int u:incoming.freePressureOutletVelocityDofs())free[u]=true;
    auto d=std::make_shared<Data>(std::move(m));d->prescribedWalls=walls;d->prescribedInlet=incoming;d->vofPrescribed=true;d->vofOpenMomentum=true;
    const int count=l.velocityCount();d->externalOutflow.resize(count);if(includeCenterTables)d->externalAbsolute.resize(count);d->boundaryLoad.resize(count);d->outflow.resize(count);
    d->openSupport=d->material.geometry().openConvectionSupport(mesh);
    std::size_t boundaryCount=0;
    if(d->openSupport){boundaryCount=d->openSupport->boundaryRows.size();
        for(int u:d->openSupport->exteriorVelocityDofs)if(!free[u])d->fixed.push_back(u);}
    else {d->openBoundarySlot.assign(count,-1);
        for(const auto& g:d->material.geometry().segments())if(!internal(g)&&d->openBoundarySlot[owner(g)]<0)d->openBoundarySlot[owner(g)]=0;
        for(int u=0;u<count;++u){if(d->openBoundarySlot[u]>=0)d->openBoundarySlot[u]=static_cast<int>(boundaryCount++);
            if(l.domain().face(mesh,l.nativeFace(u)).outwardSign!=0&&!free[u])d->fixed.push_back(u);}}
    d->openOutgoing.resize(boundaryCount);d->openIncomingContinuation.resize(boundaryCount);
    std::vector<OpenExact> affine(boundaryCount);
    for(int u:walls.velocityDofs())require(!free[u],"ambient pressure port overlaps Closed wall");
    const auto& segments=d->material.geometry().segments();const auto& gross=incoming.vofOpenMomentumFluxes();
    require(gross.size()==segments.size(),"different open momentum gross segment table");
    const auto boundarySegment=[&](std::size_t id){const auto& g=segments[id];const int i=owner(g),slot=d->boundarySlot(i);const auto& q=gross[id];
        require(q.velocityDof==i&&q.nativeRow==l.nativeFace(g.primalVelocity),"different actual open momentum segment incidence");
        d->externalOutflow[i]=finite(d->externalOutflow[i]+q.canonicalOutwardMassRate,"nonfinite canonical exterior rate sum");
        addOpenMass(d->openOutgoing[slot],q,d->material,false,1);
        if(q.incomingModel==HydroVofPressureOutletMomentumModel::InteriorContinuationZeroShear)
            addOpenMass(d->openIncomingContinuation[slot],q,d->material,true,1);
        else addOpenMass(affine[slot],q,d->material,true,q.incomingMomentum);
        if(includeCenterTables){OpenExact absRate;addOpenMass(absRate,q,d->material,false,1);addOpenMass(absRate,q,d->material,true,1);
            d->externalAbsolute[i]=finite(d->externalAbsolute[i]+absRate.rounded(),"nonfinite actual exterior gross rate sum");}
    };
    if(prepared&&d->openSupport&&!incoming.data_->productionReceipt->exteriorSegments.empty()){
        // The retained binder ids preserve the ORIGINAL GLOBAL error order.
        for(auto id:incoming.data_->productionReceipt->exteriorSegments)boundarySegment(id);
        for(int i:d->openSupport->boundaryRows)d->boundaryLoad[i]=affine[d->boundarySlot(i)].rounded();
    }else {
        for(std::size_t id=0;id<segments.size();++id){const auto& g=segments[id];const auto f=d->material.faceRates()[id];
            if(includeCenterTables){d->points.push_back(center(g,l.domain().axisymmetric()));
                if(internal(g))d->weights.push_back({{g.negativeVelocity,.5},{g.positiveVelocity,.5}});else d->weights.emplace_back();}
            if(!internal(g))boundarySegment(id);
            if(includeCenterTables){if(g.negativeVelocity>=0)d->outflow[g.negativeVelocity]+=f;if(g.positiveVelocity>=0)d->outflow[g.positiveVelocity]-=f;}
        }
        for(int i=0;i<count;++i){const int slot=d->boundarySlot(i);if(slot>=0)d->boundaryLoad[i]=affine[slot].rounded();}
    }
    // Open physical pairing has no periodic roots, partners or moment patches.
    d->topology.representative.resize(count);std::iota(d->topology.representative.begin(),d->topology.representative.end(),0);
    std::vector<bool> fixed(count);for(int u:d->fixed)fixed[u]=true;
    d->finalPairs=physicalCenteredPairs(d->material,fixed,d->momentReport,&d->externalOutflow,d->openSupport.get(),&d->openPairCoefficients,prepared?&d->outflow:nullptr);
    if(prepared)d->momentReport.maxNativeMeanDefect=incoming.data_->productionReceipt->maximumNativeMeanDefect;
    else certifyVofFreeOutflowMeans(mesh,d->material,incoming.freePressureOutletVelocityDofs(),d->momentReport,&checked);
    if(!prepared)for(int i=0;i<count;++i){const auto canonical=finite(d->material.cells()[i].outwardIntegral/d->material.deltaTime(),"nonfinite canonical open KEP outflow");
        if(includeCenterTables){const auto difference=finite(d->outflow[i]-canonical,"nonfinite actual rate/integral diagnostic");d->momentReport.rateSumMinusIntegralRateSigned+=difference;
            d->momentReport.rateSumMinusIntegralRateL1+=std::abs(difference);d->momentReport.maxAbsRateSumMinusIntegralRate=std::max(d->momentReport.maxAbsRateSumMinusIntegralRate,std::abs(difference));}d->outflow[i]=canonical;}
    if(!includeCenterTables){const auto unavailable=std::numeric_limits<long double>::quiet_NaN();
        d->momentReport.rateSumMinusIntegralRateSigned=d->momentReport.rateSumMinusIntegralRateL1=d->momentReport.maxAbsRateSumMinusIntegralRate=unavailable;}
    if(!includeCenterTables)d->openInteriorProducts=prepareOpenInteriorProducts(d->material,d->openSupport,d->openPairCoefficients,d->outflow,
        incoming.data_->productionReceipt?incoming.data_->productionReceipt->exteriorBytes:0u);
    d->productionPrepared=!includeCenterTables;
    return HydroKepConvection(std::move(d));
}
bool HydroKepConvection::isVofOpenFlow() const noexcept{return data_->vofOpenMomentum;}
HydroKepConvectionResult HydroKepConvection::applyVofOpen(const QuadtreeMesh& mesh,const HydroNativeValues& values,
    HydroMaterialSources source,LAMP::TimeLayer clock) const
{return applyVofOpenImpl(mesh,values,source,clock,true,true);}
HydroKepConvectionResult HydroKepConvection::applyVofOpenWithoutCenterDiagnostic(const QuadtreeMesh& mesh,const HydroNativeValues& values,
    HydroMaterialSources source,LAMP::TimeLayer clock) const
{return applyVofOpenImpl(mesh,values,source,clock,false,true);}
HydroKepConvectionResult HydroKepConvection::applyVofOpenMomentumLoads(const QuadtreeMesh& mesh,const HydroNativeValues& values,
    HydroMaterialSources source,LAMP::TimeLayer clock) const
{return applyVofOpenImpl(mesh,values,source,clock,false,false);}
HydroKepConvectionResult HydroKepConvection::applyVofOpenImpl(const QuadtreeMesh& mesh,const HydroNativeValues& values,
    HydroMaterialSources source,LAMP::TimeLayer clock,bool includeCenterDiagnostic,bool includePowerReport) const
{
    const auto start=OpenWall::now();
    require(data_->vofOpenMomentum&&data_->prescribedInlet,"applyVofOpen requires actual open momentum operator");
    const auto checked=CheckedVofMaterialBatch::validate(mesh,material(),source,clock);values.validateFor(layout(),mesh);const auto& incoming=*data_->prescribedInlet;
    validateConvectionBoundaryAfterMaterial(mesh,material(),incoming,source,clock,&checked);
    validateVofOpenSource(mesh,material(),*data_->prescribedWalls,&checked);
    const auto& u=values.normalVelocity();for(int i:data_->prescribedWalls->velocityDofs())require(u[i]==0,"nonzero actual Closed normal in open KEP apply");
    const auto& gross=incoming.vofOpenMomentumFluxes();const auto& segments=material().geometry().segments();const auto dt=material().deltaTime();
    HydroKepOpenApplyReport receipt;receipt.available=true;receipt.nativeRows=u.size();receipt.segments=segments.size();receipt.pairs=data_->finalPairs.size();receipt.timeLayer=clock;
    for(std::size_t i=0;i<3;++i)receipt.vofSources[i]=source.vofIds[i].value();
    receipt.validationSeconds=openSeconds(OpenWall::now()-start);
    const auto markPowerUnavailable=[](HydroKepComponentReport& report){
        const auto unavailable=std::numeric_limits<long double>::quiet_NaN();
        report.convectivePower=report.boundaryKineticPower=report.materialKineticPower=report.kineticRate=unavailable;
        report.boundaryPhysicalKineticPower=report.inletDissipationPower=report.boundaryFluxRoundingPower=unavailable;
        report.boundaryDiagonalPower=report.kineticIdentity=unavailable;
    };
    const auto executeSerial=[&]{auto mark=OpenWall::now();
    struct BoundaryPower {OpenExact work,physical,dissipation,representation;};
    struct PairPower {OpenExact lwork,rwork;};
    struct RowPower {OpenExact rp,rs,rm,rk;};
    auto d=std::make_shared<HydroKepConvectionResult::Data>(material(),incoming);
    d->centerAvailable=includeCenterDiagnostic;d->powerAvailable=includePowerReport;
    d->rates.resize(u.size());d->loads.resize(u.size());d->dissipation.resize(u.size());
    std::vector<OpenExact> rates(u.size());std::array<OpenExact,2> momentumScale,powerScale,boundaryMomentum,boundaryWork,physicalWork,dissipationWork,representationWork,diagonalPower;
    receipt.preparationSeconds=openSeconds(OpenWall::now()-mark);mark=OpenWall::now();
    for(std::size_t id=0;id<segments.size();++id){const auto& g=segments[id];
        std::optional<OpenExact> diagnostic;if(includeCenterDiagnostic)diagnostic.emplace();
        if(internal(g)){if(diagnostic){diagnostic->addProduct({material().faceRates()[id],.5L,u[g.negativeVelocity]});diagnostic->addProduct({material().faceRates()[id],.5L,u[g.positiveVelocity]});}}
        else {const auto& q=gross[id];const int i=owner(g),a=component(layout(),i);
            const auto b=q.incomingModel==HydroVofPressureOutletMomentumModel::InteriorContinuationZeroShear?u[i]:q.incomingMomentum;
            OpenExact load;std::optional<BoundaryPower> power;if(includePowerReport)power.emplace();OpenExact eta;
            addOpenMass(load,q,material(),false,u[i]);addOpenMass(load,q,material(),true,-b);
            rates[i].add(load);boundaryMomentum[a].add(load);momentumScale[a].addAbsolute(load);
            if(power){auto& work=power->work;auto& physical=power->physical;auto& dissipation=power->dissipation;auto& representation=power->representation;
                addOpenMass(diagonalPower[a],q,material(),false,u[i],u[i]);
                if(q.incomingModel==HydroVofPressureOutletMomentumModel::InteriorContinuationZeroShear)
                    addOpenMass(diagonalPower[a],q,material(),true,-1,u[i],u[i]);
                addOpenMass(work,q,material(),false,u[i],u[i]);addOpenMass(work,q,material(),true,-1,u[i],b);
                work.addProduct({-.5L,q.canonicalOutwardMassRate,u[i],u[i]});boundaryWork[a].add(work);
                addOpenMass(physical,q,material(),false,.5L,u[i],u[i]);addOpenMass(physical,q,material(),true,-.5L,b,b);physicalWork[a].add(physical);
                addOpenMass(dissipation,q,material(),true,.5L,u[i],u[i]);addOpenMass(dissipation,q,material(),true,-1,u[i],b);
                addOpenMass(dissipation,q,material(),true,.5L,b,b);require(dissipation.compare(OpenExact{})>=0,"negative exact gross inlet dissipation");dissipationWork[a].add(dissipation);
                representation.addProduct({-.5L,q.canonicalOutwardMassRate,u[i],u[i]});addOpenMass(representation,q,material(),false,.5L,u[i],u[i]);addOpenMass(representation,q,material(),true,-.5L,u[i],u[i]);representationWork[a].add(representation);
            }
            eta.add(q.canonicalOutwardMassRate);addOpenMass(eta,q,material(),false,-1);addOpenMass(eta,q,material(),true,1);
            d->report.components[a].maximumBoundaryCanonicalGrossRateDefect=std::max(d->report.components[a].maximumBoundaryCanonicalGrossRateDefect,std::abs(eta.rounded()));
            if(power){powerScale[a].addAbsolute(power->work);powerScale[a].addAbsolute(power->physical);powerScale[a].addAbsolute(power->dissipation);powerScale[a].addAbsolute(power->representation);}
            if(diagnostic)diagnostic->add(load,g.negativeVelocity>=0?1:-1);
        }
        if(diagnostic)d->center.push_back(diagnostic->rounded());
    }
    receipt.segmentSeconds=openSeconds(OpenWall::now()-mark);mark=OpenWall::now();
    for(std::size_t i=0;i<u.size();++i){const int a=component(layout(),static_cast<int>(i));OpenExact term;
        std::optional<OpenExact> work;if(includePowerReport)work.emplace();
        term.addProduct({.5L,data_->outflow[i],u[i]});term.addProduct({-.5L,data_->externalOutflow[i],u[i]});rates[i].add(term);
        if(work){work->addProduct({.5L,data_->outflow[i],u[i],u[i]});work->addProduct({-.5L,data_->externalOutflow[i],u[i],u[i]});
            diagonalPower[a].addProduct({-.5L,data_->externalOutflow[i],u[i],u[i]});}
        momentumScale[a].addAbsolute(term);if(work)powerScale[a].addAbsolute(*work);}
    receipt.diagonalSeconds=openSeconds(OpenWall::now()-mark);mark=OpenWall::now();
    for(const auto& e:data_->finalPairs){const int i=e.first,j=e.second,a=component(layout(),i);OpenExact left,right;
        std::optional<PairPower> power;if(includePowerReport)power.emplace();
        left.addProduct({e.coefficient,u[j]});right.addProduct({e.coefficient,u[i]});rates[i].add(left);rates[j].add(right,-1);
        if(power){power->lwork.addProduct({u[i],e.coefficient,u[j]});power->rwork.addProduct({u[j],e.coefficient,u[i]});}
        momentumScale[a].addAbsolute(left);momentumScale[a].addAbsolute(right);
        if(power){powerScale[a].addAbsolute(power->lwork);powerScale[a].addAbsolute(power->rwork);}}
    receipt.pairSeconds=openSeconds(OpenWall::now()-mark);mark=OpenWall::now();
    std::array<OpenExact,2> momentum,power,sourcePower,kinetic;
    for(std::size_t i=0;i<u.size();++i){const int a=component(layout(),static_cast<int>(i));const auto& cell=material().cells()[i];d->rates[i]=rates[i].rounded();
        OpenExact load;load.add(rates[i],-1);
        if(!material().usesConservativeVofMomentum())load.addQuotientProduct({u[i],cell.integratedDefect},dt);d->loads[i]=load.rounded();
        std::optional<RowPower> rowPower;if(includePowerReport)rowPower.emplace();
        if(rowPower){auto& rp=rowPower->rp;auto& rs=rowPower->rs;auto& rm=rowPower->rm;auto& rk=rowPower->rk;
            rp.addProduct({u[i],d->rates[i]});rs.addQuotientProduct({material().usesConservativeVofMomentum()?-.5L:.5L,u[i],u[i],cell.integratedDefect},dt);rm.addQuotientProduct({.5L,u[i],u[i],cell.massChange},dt);
            rk.addProduct({u[i],d->loads[i]});rk.add(rm,-1);}
        momentum[a].add(d->rates[i]);if(rowPower){power[a].add(rowPower->rp);sourcePower[a].add(rowPower->rs);kinetic[a].add(rowPower->rk);}
        OpenExact rate;rate.add(d->rates[i]);momentumScale[a].addAbsolute(rate);
        if(rowPower){powerScale[a].addAbsolute(rowPower->rp);powerScale[a].addAbsolute(rowPower->rs);powerScale[a].addAbsolute(rowPower->rs);powerScale[a].addAbsolute(rowPower->rm);powerScale[a].addAbsolute(rowPower->rk);}}
    receipt.rowSeconds=openSeconds(OpenWall::now()-mark);mark=OpenWall::now();
    for(int a=0;a<2;++a){auto& r=d->report.components[a];
        if(includePowerReport){r.momentumRate=momentum[a].rounded();r.convectivePower=power[a].rounded();r.materialKineticPower=sourcePower[a].rounded();r.kineticRate=kinetic[a].rounded();
            r.boundaryMomentumRate=boundaryMomentum[a].rounded();r.boundaryKineticPower=boundaryWork[a].rounded();r.boundaryPhysicalKineticPower=physicalWork[a].rounded();r.inletDissipationPower=dissipationWork[a].rounded();r.boundaryFluxRoundingPower=representationWork[a].rounded();r.boundaryDiagonalPower=diagonalPower[a].rounded();
            identity(r.boundaryKineticPower-r.boundaryPhysicalKineticPower-r.inletDissipationPower-r.boundaryFluxRoundingPower,powerScale[a].rounded(),"gross boundary physical/dissipation/representation decomposition failed");
            r.momentumIdentity=identity(r.momentumRate-r.boundaryMomentumRate,momentumScale[a].rounded(),"compatible open momentum balance failed");
            r.kineticIdentity=identity(r.kineticRate-r.materialKineticPower+r.boundaryKineticPower+r.dissipationPower,powerScale[a].rounded(),"nonconservative open kinetic identity failed");}
        else{r.momentumRate=momentum[a].rounded();r.boundaryMomentumRate=boundaryMomentum[a].rounded();
            r.momentumIdentity=identity(r.momentumRate-r.boundaryMomentumRate,momentumScale[a].rounded(),"compatible open momentum balance failed");
            markPowerUnavailable(r);}}
    receipt.ledgerSeconds=openSeconds(OpenWall::now()-mark);
    return d;
    };

#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
    const int entryErrno=errno;OpenFloatingState entry;const bool qualified=openQualified(entry);
    receipt.parallelRequested=!omp_in_parallel()&&omp_get_max_threads()>1;
    if(entry.valid)entry.restore();errno=entryErrno;
    std::shared_ptr<const OpenApplyPlan> applyPlan;
    if(qualified){try{std::lock_guard<std::mutex> lock(data_->openApplyPlanMutex);
        if(!data_->openApplyPlanAttempted){data_->openApplyPlanAttempted=true;
            data_->openApplyPlan=prepareOpenApplyPlan(segments,data_->finalPairs,u.size());}
        applyPlan=data_->openApplyPlan;}catch(...){}
        entry.restore();errno=entryErrno;}
    receipt.parallelEligible=qualified&&static_cast<bool>(applyPlan);
    if(receipt.parallelEligible){
        const auto optionalStart=OpenWall::now();
        try{
            const auto executeParallel=[&]{
                // Cleanup is declared first and armed before leaving: optional
                // storage destruction cannot leak allocator FP/errno effects.
                OpenCleanup cleanup;
                struct BoundaryPower {OpenExact work,physical,dissipation,representation,diagonal;};
                struct BoundaryWork {OpenExact load;std::optional<BoundaryPower> power;};
                struct SegmentWork {long double center=0,eta=0;std::optional<BoundaryWork> boundary;};
                struct DiagonalPower {OpenExact work,diagonal;};
                struct DiagonalWork {OpenExact term;std::optional<DiagonalPower> power;};
                struct PairPower {OpenExact lwork,rwork;};
                struct PairWork {OpenExact left,right;std::optional<PairPower> power;};
                struct RowPower {OpenExact rp,rs,rm,rk;};
                struct RowWork {std::optional<RowPower> power;OpenExact rate;};
                const auto& plan=*applyPlan;
                // Four binary64 factors and one original dt have <11000
                // exponent-span bits; uint32 sums add at most32 carry bits.
                // 8KiB per retained rational plus128KiB/worker conservatively
                // covers bigint capacity/temporaries, not just pointer sizes.
                // Completed stages merge and clear their channels before reuse.
                // At most the original sixteen exact partials live per worker.
                struct LedgerPartials {
                    std::array<OpenExact,16> ledger;
                    std::array<long double,2> eta{};
                };
                constexpr std::size_t rationalBytes=8192;
                // Retained stage partials are additional to the original
                // per-worker temporary budget, including all 16 Data objects.
                constexpr std::size_t workerBytes=128u*1024u+sizeof(LedgerPartials)+16*rationalBytes;
                // sizeof charges inline optional metadata even when disengaged.
                // Only the requested exact payloads exist in each exclusive stage.
                const auto stageBytes=std::max({sizeof(SegmentWork)+(includePowerReport?6u:1u)*rationalBytes,
                    sizeof(DiagonalWork)+(includePowerReport?3u:1u)*rationalBytes,
                    sizeof(PairWork)+(includePowerReport?4u:2u)*rationalBytes,
                    sizeof(RowWork)+(includePowerReport?5u:1u)*rationalBytes});
                const auto fixed=plan.bytes+openBatch*stageBytes
                    +(incoming.data_->productionReceipt?incoming.data_->productionReceipt->exteriorBytes:0u);
                if(fixed>openScratchLimit||openScratchLimit-fixed<2*workerBytes)throw OpenOptionalFailure{};
                const auto teamCount=std::min<std::size_t>({static_cast<std::size_t>(omp_get_max_threads()),
                    std::max(u.size(),segments.size()),(openScratchLimit-fixed)/workerBytes});
                if(teamCount<2)throw OpenOptionalFailure{};
                const int team=static_cast<int>(teamCount);
                receipt.scratchBytes=fixed+teamCount*workerBytes;
                auto sw=std::make_unique<SegmentWork[]>(openBatch);
                std::unique_ptr<DiagonalWork[]> dw;std::unique_ptr<PairWork[]> pw;std::unique_ptr<RowWork[]> rw;
                // Only optional allocator/destructor effects are neutralized.
                // No original numerical work or reduction is in this scope.
                const auto nextStageBuffer=[&](auto& previous,auto& next){
                    const int error=errno;OpenFloatingState numerical;
                    if(!numerical.capture())throw OpenOptionalFailure{};
                    previous.reset();
                    using Work=typename std::remove_reference_t<decltype(next)>::element_type;
                    next=std::make_unique<Work[]>(openBatch);
                    if(includePowerReport)for(std::size_t k=0;k<openBatch;++k)next[k].power.emplace();
                    numerical.restore();errno=error;
                };
                auto partials=std::make_unique<LedgerPartials[]>(teamCount);
                auto d=std::make_shared<HydroKepConvectionResult::Data>(material(),incoming);
                d->centerAvailable=includeCenterDiagnostic;d->powerAvailable=includePowerReport;
                d->rates.resize(u.size());d->loads.resize(u.size());d->dissipation.resize(u.size());
                if(includeCenterDiagnostic)d->center.resize(segments.size());
                std::vector<OpenExact> rates(u.size());
                std::array<OpenExact,2> momentumScale,powerScale,boundaryMomentum,boundaryWork,physicalWork,dissipationWork,representationWork,diagonalPower;
                OpenCooperativeTeam cooperative(team);
                if(errno!=entryErrno)throw OpenOptionalFailure{};
                entry.restore();errno=entryErrno;
                auto mark=OpenWall::now();receipt.preparationSeconds=openSeconds(mark-optionalStart);
                return cooperative.run([&]{
                const auto rows=[&](std::size_t count,auto&& kernel,unsigned& workers){
                    if(!cooperative.rows(count,std::forward<decltype(kernel)>(kernel),workers))throw OpenOptionalFailure{};};
                for(std::size_t begin=0;begin<segments.size();begin+=openBatch){const auto count=std::min(openBatch,segments.size()-begin),end=begin+count;
                    rows(count,[&](std::size_t k){const auto id=begin+k;const auto& g=segments[id];auto& w=sw[k];w.boundary.reset();w.eta=0;
                        std::optional<OpenExact> diagnostic;if(includeCenterDiagnostic)diagnostic.emplace();
                        if(internal(g)){if(diagnostic){diagnostic->addProduct({material().faceRates()[id],.5L,u[g.negativeVelocity]});diagnostic->addProduct({material().faceRates()[id],.5L,u[g.positiveVelocity]});}}
                        else {w.boundary.emplace();auto& b=*w.boundary;if(includePowerReport)b.power.emplace();const auto& q=gross[id];const int i=owner(g);
                            const auto trace=q.incomingModel==HydroVofPressureOutletMomentumModel::InteriorContinuationZeroShear?u[i]:q.incomingMomentum;
                            addOpenMass(b.load,q,material(),false,u[i]);addOpenMass(b.load,q,material(),true,-trace);
                            if(b.power){auto& p=*b.power;
                                addOpenMass(p.diagonal,q,material(),false,u[i],u[i]);
                                if(q.incomingModel==HydroVofPressureOutletMomentumModel::InteriorContinuationZeroShear)addOpenMass(p.diagonal,q,material(),true,-1,u[i],u[i]);
                                addOpenMass(p.work,q,material(),false,u[i],u[i]);addOpenMass(p.work,q,material(),true,-1,u[i],trace);
                                p.work.addProduct({-.5L,q.canonicalOutwardMassRate,u[i],u[i]});
                                addOpenMass(p.physical,q,material(),false,.5L,u[i],u[i]);addOpenMass(p.physical,q,material(),true,-.5L,trace,trace);
                                addOpenMass(p.dissipation,q,material(),true,.5L,u[i],u[i]);addOpenMass(p.dissipation,q,material(),true,-1,u[i],trace);
                                addOpenMass(p.dissipation,q,material(),true,.5L,trace,trace);require(p.dissipation.compare(OpenExact{})>=0,"negative exact gross inlet dissipation");
                                p.representation.addProduct({-.5L,q.canonicalOutwardMassRate,u[i],u[i]});addOpenMass(p.representation,q,material(),false,.5L,u[i],u[i]);addOpenMass(p.representation,q,material(),true,-.5L,u[i],u[i]);}
                            OpenExact eta;eta.add(q.canonicalOutwardMassRate);addOpenMass(eta,q,material(),false,-1);addOpenMass(eta,q,material(),true,1);w.eta=std::abs(eta.rounded());
                            if(diagnostic)diagnostic->add(b.load,g.negativeVelocity>=0?1:-1);}
                        if(diagnostic){w.center=diagnostic->rounded();d->center[id]=w.center;}
                        if(w.boundary){const auto& b=*w.boundary;const int a=component(layout(),owner(g));
                            auto& partial=partials[static_cast<std::size_t>(omp_get_thread_num())];auto& sum=partial.ledger;
                            sum[a].add(b.load);sum[2+a].addAbsolute(b.load);
                            if(b.power){const auto& p=*b.power;sum[4+a].add(p.diagonal);
                                sum[6+a].add(p.work);sum[8+a].add(p.physical);sum[10+a].add(p.dissipation);sum[12+a].add(p.representation);}
                            partial.eta[a]=std::max(partial.eta[a],w.eta);
                            if(b.power){const auto& p=*b.power;sum[14+a].addAbsolute(p.work);sum[14+a].addAbsolute(p.physical);
                                sum[14+a].addAbsolute(p.dissipation);sum[14+a].addAbsolute(p.representation);}}},receipt.contributionWorkers);
                    const auto touchedBegin=plan.boundaryChunkOffsets[begin/openBatch],touchedEnd=plan.boundaryChunkOffsets[begin/openBatch+1];
                    // A geometry-empty gather has no kernel or numerical operands.
                    if(touchedBegin!=touchedEnd)rows(touchedEnd-touchedBegin,[&](std::size_t local){const auto i=plan.boundaryChunkRows[touchedBegin+local];
                        auto first=plan.boundaryItems.begin()+plan.boundaryOffsets[i],last=plan.boundaryItems.begin()+plan.boundaryOffsets[i+1];
                        for(auto at=std::lower_bound(first,last,begin);at!=last&&*at<end;++at)rates[i].add(sw[*at-begin].boundary->load);},receipt.gatherWorkers);
                }
                // Exact worker partials are never rounded. Channels and eta
                // are merged only from the actual team, then cleared for reuse.
                rows(includePowerReport?16:4,[&](std::size_t channel){const int a=static_cast<int>(channel%2),ledger=static_cast<int>(channel/2);
                    OpenExact* sum=nullptr;
                    switch(ledger){case 0:sum=&boundaryMomentum[a];break;case 1:sum=&momentumScale[a];break;
                    case 2:sum=&diagonalPower[a];break;case 3:sum=&boundaryWork[a];break;case 4:sum=&physicalWork[a];break;
                    case 5:sum=&dissipationWork[a];break;case 6:sum=&representationWork[a];break;case 7:sum=&powerScale[a];break;}
                    for(int tid=0;tid<omp_get_num_threads();++tid){auto& partial=partials[static_cast<std::size_t>(tid)];
                        sum->add(partial.ledger[channel]);partial.ledger[channel].reset();
                        if(channel<2){auto& maximum=d->report.components[a].maximumBoundaryCanonicalGrossRateDefect;
                            maximum=std::max(maximum,partial.eta[a]);partial.eta[a]=0;}}},receipt.gatherWorkers);
                receipt.segmentSeconds=openSeconds(OpenWall::now()-mark);mark=OpenWall::now();
                nextStageBuffer(sw,dw);
                for(std::size_t begin=0;begin<u.size();begin+=openBatch){const auto count=std::min(openBatch,u.size()-begin);
                    rows(count,[&](std::size_t k){const auto i=begin+k;auto& w=dw[k];w.term.reset();
                        if(w.power){w.power->work.reset();w.power->diagonal.reset();}
                        w.term.addProduct({.5L,data_->outflow[i],u[i]});w.term.addProduct({-.5L,data_->externalOutflow[i],u[i]});rates[i].add(w.term);
                        if(w.power){auto& p=*w.power;
                            p.work.addProduct({.5L,data_->outflow[i],u[i],u[i]});p.work.addProduct({-.5L,data_->externalOutflow[i],u[i],u[i]});
                            p.diagonal.addProduct({-.5L,data_->externalOutflow[i],u[i],u[i]});}
                        const int a=component(layout(),static_cast<int>(i));auto& sum=partials[static_cast<std::size_t>(omp_get_thread_num())].ledger;
                        if(w.power){sum[a].add(w.power->diagonal);sum[2+a].addAbsolute(w.term);sum[4+a].addAbsolute(w.power->work);}
                        else sum[a].addAbsolute(w.term);},receipt.contributionWorkers);
                }
                rows(includePowerReport?6:2,[&](std::size_t channel){const int a=static_cast<int>(channel%2),ledger=static_cast<int>(channel/2);
                    OpenExact* sum=includePowerReport?(ledger==0?&diagonalPower[a]:ledger==1?&momentumScale[a]:&powerScale[a]):&momentumScale[a];
                    for(int tid=0;tid<omp_get_num_threads();++tid){auto& partial=partials[static_cast<std::size_t>(tid)].ledger[channel];
                        sum->add(partial);partial.reset();}},receipt.gatherWorkers);
                receipt.diagonalSeconds=openSeconds(OpenWall::now()-mark);mark=OpenWall::now();
                nextStageBuffer(dw,pw);
                for(std::size_t begin=0;begin<data_->finalPairs.size();begin+=openBatch){const auto count=std::min(openBatch,data_->finalPairs.size()-begin),end=begin+count;
                    rows(count,[&](std::size_t k){const auto& e=data_->finalPairs[begin+k];auto& w=pw[k];w.left.reset();w.right.reset();
                        if(w.power){w.power->lwork.reset();w.power->rwork.reset();}
                        w.left.addProduct({e.coefficient,u[e.second]});w.right.addProduct({e.coefficient,u[e.first]});
                        if(w.power){w.power->lwork.addProduct({u[e.first],e.coefficient,u[e.second]});w.power->rwork.addProduct({u[e.second],e.coefficient,u[e.first]});}
                        const int a=component(layout(),e.first);auto& sum=partials[static_cast<std::size_t>(omp_get_thread_num())].ledger;
                        sum[a].addAbsolute(w.left);sum[a].addAbsolute(w.right);
                        if(w.power){sum[2+a].addAbsolute(w.power->lwork);sum[2+a].addAbsolute(w.power->rwork);}},receipt.contributionWorkers);
                    const auto touchedBegin=plan.pairChunkOffsets[begin/openBatch],touchedEnd=plan.pairChunkOffsets[begin/openBatch+1];
                    if(touchedBegin!=touchedEnd)rows(touchedEnd-touchedBegin,[&](std::size_t local){const auto i=plan.pairChunkRows[touchedBegin+local];
                        auto first=plan.pairItems.begin()+plan.pairOffsets[i],last=plan.pairItems.begin()+plan.pairOffsets[i+1];
                        auto at=std::lower_bound(first,last,begin,[](const auto& e,std::size_t id){return e.id<id;});
                        for(;at!=last&&at->id<end;++at){const auto& w=pw[at->id-begin];rates[i].add(at->sign>0?w.left:w.right,at->sign);}},receipt.gatherWorkers);
                }
                // Exact reassociation only: partials never undergo native RNE.
                // Runtime-reduced teams merge only the slots actually written.
                rows(includePowerReport?4:2,[&](std::size_t channel){const int a=static_cast<int>(channel%2),ledger=static_cast<int>(channel/2);
                    auto& sum=ledger==0?momentumScale[a]:powerScale[a];
                    for(int tid=0;tid<omp_get_num_threads();++tid){auto& partial=partials[static_cast<std::size_t>(tid)].ledger[channel];
                        sum.add(partial);partial.reset();}},receipt.gatherWorkers);
                receipt.pairSeconds=openSeconds(OpenWall::now()-mark);mark=OpenWall::now();
                nextStageBuffer(pw,rw);
                std::array<OpenExact,2> momentum,power,sourcePower,kinetic;
                for(std::size_t begin=0;begin<u.size();begin+=openBatch){const auto count=std::min(openBatch,u.size()-begin);
                    rows(count,[&](std::size_t k){const auto i=begin+k;const auto& cell=material().cells()[i];auto& w=rw[k];
                        if(w.power){w.power->rp.reset();w.power->rs.reset();w.power->rm.reset();w.power->rk.reset();}
                        w.rate.reset();d->rates[i]=rates[i].rounded();
                        OpenExact load;load.add(rates[i],-1);if(!material().usesConservativeVofMomentum())load.addQuotientProduct({u[i],cell.integratedDefect},dt);d->loads[i]=load.rounded();
                        if(w.power){auto& p=*w.power;
                            p.rp.addProduct({u[i],d->rates[i]});p.rs.addQuotientProduct({material().usesConservativeVofMomentum()?-.5L:.5L,u[i],u[i],cell.integratedDefect},dt);
                            p.rm.addQuotientProduct({.5L,u[i],u[i],cell.massChange},dt);p.rk.addProduct({u[i],d->loads[i]});p.rk.add(p.rm,-1);}
                        w.rate.add(d->rates[i]);
                        const int a=component(layout(),static_cast<int>(i));auto& sum=partials[static_cast<std::size_t>(omp_get_thread_num())].ledger;
                        sum[a].add(d->rates[i]);
                        if(w.power){const auto& p=*w.power;sum[2+a].add(p.rp);sum[4+a].add(p.rs);sum[6+a].add(p.rk);
                            sum[8+a].addAbsolute(w.rate);sum[10+a].addAbsolute(p.rp);sum[10+a].addAbsolute(p.rs);
                            sum[10+a].addAbsolute(p.rs);sum[10+a].addAbsolute(p.rm);sum[10+a].addAbsolute(p.rk);}
                        else sum[2+a].addAbsolute(w.rate);},receipt.contributionWorkers);
                }
                rows(includePowerReport?12:4,[&](std::size_t channel){const int a=static_cast<int>(channel%2),ledger=static_cast<int>(channel/2);
                    OpenExact* sum=nullptr;
                    if(includePowerReport){switch(ledger){case 0:sum=&momentum[a];break;case 1:sum=&power[a];break;case 2:sum=&sourcePower[a];break;
                        case 3:sum=&kinetic[a];break;case 4:sum=&momentumScale[a];break;case 5:sum=&powerScale[a];break;}}
                    else sum=ledger==0?&momentum[a]:&momentumScale[a];
                    for(int tid=0;tid<omp_get_num_threads();++tid){auto& partial=partials[static_cast<std::size_t>(tid)].ledger[channel];
                        sum->add(partial);partial.reset();}},receipt.gatherWorkers);
                receipt.rowSeconds=openSeconds(OpenWall::now()-mark);mark=OpenWall::now();
                for(int a=0;a<2;++a){auto& r=d->report.components[a];
                    if(includePowerReport){r.momentumRate=momentum[a].rounded();r.convectivePower=power[a].rounded();r.materialKineticPower=sourcePower[a].rounded();r.kineticRate=kinetic[a].rounded();
                        r.boundaryMomentumRate=boundaryMomentum[a].rounded();r.boundaryKineticPower=boundaryWork[a].rounded();r.boundaryPhysicalKineticPower=physicalWork[a].rounded();r.inletDissipationPower=dissipationWork[a].rounded();r.boundaryFluxRoundingPower=representationWork[a].rounded();r.boundaryDiagonalPower=diagonalPower[a].rounded();
                        identity(r.boundaryKineticPower-r.boundaryPhysicalKineticPower-r.inletDissipationPower-r.boundaryFluxRoundingPower,powerScale[a].rounded(),"gross boundary physical/dissipation/representation decomposition failed");
                        r.momentumIdentity=identity(r.momentumRate-r.boundaryMomentumRate,momentumScale[a].rounded(),"compatible open momentum balance failed");
                        r.kineticIdentity=identity(r.kineticRate-r.materialKineticPower+r.boundaryKineticPower+r.dissipationPower,powerScale[a].rounded(),"nonconservative open kinetic identity failed");}
                    else{r.momentumRate=momentum[a].rounded();r.boundaryMomentumRate=boundaryMomentum[a].rounded();
                        r.momentumIdentity=identity(r.momentumRate-r.boundaryMomentumRate,momentumScale[a].rounded(),"compatible open momentum balance failed");
                        markPowerUnavailable(r);}}
                receipt.ledgerSeconds=openSeconds(OpenWall::now()-mark);cleanup.arm();return d;
                });
            };
            auto d=executeParallel();receipt.parallelUsed=true;receipt.totalSeconds=openSeconds(OpenWall::now()-start);d->openApply=receipt;
            return HydroKepConvectionResult(std::move(d));
        }catch(...){
            // Original-order literal replay chooses the original first error;
            // no partial candidate or optional worker exception escapes.
            const double failed=openSeconds(OpenWall::now()-optionalStart);
            const bool requested=receipt.parallelRequested,eligible=receipt.parallelEligible;
            receipt.preparationSeconds=receipt.segmentSeconds=receipt.diagonalSeconds=receipt.pairSeconds=receipt.rowSeconds=receipt.ledgerSeconds=0;
            receipt.contributionWorkers=receipt.gatherWorkers=0;receipt.scratchBytes=0;receipt.parallelUsed=false;
            receipt.parallelRequested=requested;receipt.parallelEligible=eligible;receipt.parallelFallback=true;receipt.fallbackSeconds=failed;
            entry.restore();errno=entryErrno;
        }
    }
#endif
    auto d=executeSerial();receipt.totalSeconds=openSeconds(OpenWall::now()-start);d->openApply=receipt;
    return HydroKepConvectionResult(std::move(d));
}

HydroKepConvectionResult HydroKepConvection::applyVofOpenProductionLoads(const QuadtreeMesh& mesh,const HydroNativeValues& values,
    HydroMaterialSources source,LAMP::TimeLayer clock) const
{return applyVofOpenProductionImpl(mesh,values,source,clock,false);}
HydroKepConvectionResult HydroKepConvection::applyVofOpenBoundaryBudget(const QuadtreeMesh& mesh,const HydroNativeValues& values,
    HydroMaterialSources source,LAMP::TimeLayer clock) const
{return applyVofOpenProductionImpl(mesh,values,source,clock,true);}
HydroKepConvectionResult HydroKepConvection::applyVofOpenProductionImpl(const QuadtreeMesh& mesh,const HydroNativeValues& values,
    HydroMaterialSources source,LAMP::TimeLayer clock,bool includeBoundaryBudget) const
{
    const auto start=OpenWall::now();
    require(data_->vofOpenMomentum&&data_->prescribedInlet,"applyVofOpen requires actual open momentum operator");
    const auto& incoming=*data_->prescribedInlet;
    const auto expected=material().sources();const auto expectedClock=material().timeLayer();
    const bool sameSource=source.model==expected.model&&source.begin==expected.begin&&source.midpoint==expected.midpoint
        &&source.end==expected.end&&source.advector==expected.advector&&source.vofStamp&&expected.vofStamp
        &&source.vofStamp->sameInstance(*expected.vofStamp)&&source.vofIds[0]==expected.vofIds[0]
        &&source.vofIds[1]==expected.vofIds[1]&&source.vofIds[2]==expected.vofIds[2];
    const bool sameClock=clock.step==expectedClock.step&&clock.begin==expectedClock.begin&&clock.end==expectedClock.end&&clock.position==expectedClock.position;
    if(data_->productionPrepared&&sameSource&&sameClock&&wallExecutionQualified()){
        // Only immutable contents were certified by the dedicated factory.
        // This public request still checks the actual current mutable mesh.
        values.validateFor(layout(),mesh);
    }else {
        const auto checked=CheckedVofMaterialBatch::validate(mesh,material(),source,clock);values.validateFor(layout(),mesh);
        validateConvectionBoundaryAfterMaterial(mesh,material(),incoming,source,clock,&checked);
        validateVofOpenSource(mesh,material(),*data_->prescribedWalls,&checked);
    }
    const auto& u=values.normalVelocity();for(int i:data_->prescribedWalls->velocityDofs())require(u[i]==0,"nonzero actual Closed normal in open KEP apply");
    const auto& gross=incoming.vofOpenMomentumFluxes();const auto& segments=material().geometry().segments();const auto dt=material().deltaTime();
    HydroKepOpenApplyReport receipt;receipt.available=true;receipt.nativeRows=u.size();receipt.segments=segments.size();receipt.pairs=data_->finalPairs.size();receipt.timeLayer=clock;
    for(std::size_t i=0;i<3;++i)receipt.vofSources[i]=source.vofIds[i].value();
    receipt.validationSeconds=openSeconds(OpenWall::now()-start);
    const auto commonPreparationStart=OpenWall::now();
    bool arithmeticQualified=false;
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
    {const int error=errno;OpenFloatingState state;arithmeticQualified=openArithmeticQualified(state);
        if(state.valid)state.restore();errno=error;}
#endif
    // Execute the original environment probe once on BOTH literal and team
    // requests before optional setup can be shielded; retain its native flags.
    if(arithmeticQualified){OpenExact qualifiedEnvironment;(void)qualifiedEnvironment;}
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
    OpenCleanup productsCleanup; // Restores successful numerics after fresh-buffer destruction only.
#endif
    std::unique_ptr<OpenProducts::Value[]> productValues;
    const auto retainedProducts=data_->openInteriorProducts?data_->openInteriorProducts->bytes:0u;
    const auto retainedSupport=data_->openSupport?data_->openSupport->bytes:0u;
    const auto retainedExterior=incoming.data_->productionReceipt?incoming.data_->productionReceipt->exteriorBytes:0u;
    bool useProducts=arithmeticQualified&&bool(data_->openInteriorProducts);
    std::size_t freshProductBytes=0;
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
    if(useProducts){const int error=errno;OpenFloatingState state;
        if(!state.capture()){useProducts=false;errno=error;}
        else {try{
                constexpr auto reserve=2u*1024u*1024u;
                if(retainedSupport>openScratchLimit-reserve||retainedExterior>openScratchLimit-reserve-retainedSupport
                    ||retainedProducts>openScratchLimit-reserve-retainedSupport-retainedExterior)throw OpenOptionalFailure{};
                const auto available=openScratchLimit-reserve-retainedSupport-retainedExterior-retainedProducts;
                if(sizeof(productValues)>available||u.size()>(available-sizeof(productValues))/sizeof(OpenProducts::Value))throw OpenOptionalFailure{};
                productValues=std::make_unique<OpenProducts::Value[]>(u.size());freshProductBytes=sizeof(productValues)+u.size()*sizeof(OpenProducts::Value);
            }catch(...){productValues.reset();useProducts=false;}
            state.restore();errno=error;}}
#endif
    const auto fixedProductStorage=retainedSupport+retainedProducts+retainedExterior+freshProductBytes;
    const auto commonPreparationDuration=OpenWall::now()-commonPreparationStart;
    struct BudgetWork {std::array<std::array<OpenExact,6>,2> value;OpenExact p,b,v,c,s,d,zero;};
    struct RowAccumulator {OpenExact sum,negative;std::optional<BudgetWork> budget;
        std::array<OpenProducts::Factor,64> productFactors;std::array<std::uint32_t,64> productRows;};
    const auto reset=[&](OpenExact& value){if(arithmeticQualified)value.resetQualifiedStorage();else value.reset();};
    const auto round=[&](const OpenExact& value){return arithmeticQualified?value.roundedQualified():value.rounded();};
    const auto pairTerm=[&](OpenExact& value,long double a,long double b){
        if(arithmeticQualified)value.addBinary64Product(a,b,0);else value.addProduct({a,b});};
    const auto diagonalTerms=[&](OpenExact& value,std::size_t i){
        if(arithmeticQualified){value.addBinary64Product(data_->outflow[i],u[i],-1);value.addBinary64Product(-data_->externalOutflow[i],u[i],-1);}
        else {value.addProduct({.5L,data_->outflow[i],u[i]});value.addProduct({-.5L,data_->externalOutflow[i],u[i]});}};
    const auto makeResult=[&]{auto out=std::make_shared<HydroKepConvectionResult::Data>(material(),incoming);
        out->centerAvailable=out->powerAvailable=out->componentAvailable=false;
        out->rates.resize(u.size());out->loads.resize(u.size());out->dissipation.resize(u.size());
        if(includeBoundaryBudget)out->boundaryBudget.emplace();const auto absent=std::numeric_limits<long double>::quiet_NaN();
        for(auto& q:out->report.components){q.momentumRate=q.boundaryMomentumRate=q.convectivePower=absent;
            q.boundaryKineticPower=q.materialKineticPower=q.kineticRate=absent;
            q.boundaryPhysicalKineticPower=q.inletDissipationPower=q.boundaryFluxRoundingPower=absent;
            q.maximumBoundaryCanonicalGrossRateDefect=q.boundaryDiagonalPower=q.momentumIdentity=q.kineticIdentity=absent;}
        return out;};
    // p=.5*A*U^2, b=.5*B*ub^2, v=.5*B*U^2, c=B*U*ub,
    // s=.5*canonical*U^2. A/B retain both exact phase products/original dt.
    const auto budgetSegment=[&](BudgetWork& w,std::size_t id){const auto& q=gross[id];const int i=owner(segments[id]),a=component(layout(),i);
        const bool continuation=q.incomingModel==HydroVofPressureOutletMomentumModel::InteriorContinuationZeroShear;
        const auto ub=continuation?u[i]:q.incomingMomentum;auto& value=w.value[a];
        for(auto* term:{&w.p,&w.b,&w.v,&w.c,&w.s,&w.d})reset(*term);
        addOpenMass(w.p,q,material(),false,.5L,u[i],u[i]);
        addOpenMass(w.v,q,material(),true,.5L,u[i],u[i]);
        if(continuation){w.b.add(w.v);w.c.add(w.v);w.c.add(w.v);}
        else {addOpenMass(w.b,q,material(),true,.5L,ub,ub);addOpenMass(w.c,q,material(),true,1,u[i],ub);}
        w.s.addProduct({.5L,q.canonicalOutwardMassRate,u[i],u[i]});
        value[1].add(w.p);value[1].add(w.p);value[1].add(w.c,-1);value[1].add(w.s,-1);
        value[2].add(w.p);value[2].add(w.b,-1);
        w.d.add(w.v);w.d.add(w.c,-1);w.d.add(w.b);
        require(w.d.compare(w.zero)>=0,"negative exact gross inlet dissipation");value[3].add(w.d);
        value[4].add(w.p);value[4].add(w.v,-1);value[4].add(w.s,-1);
        value[5].add(w.p);value[5].add(w.p);if(continuation){value[5].add(w.v,-1);value[5].add(w.v,-1);}
    };
    const auto storeRow=[&](std::size_t i,OpenExact& sum,RowAccumulator& scratch,HydroKepConvectionResult::Data& out){
        const auto rate=round(sum);out.rates[i]=rate;
        if(material().usesConservativeVofMomentum()){
            // Symmetric RNE, including nonzero exact values below the quantum.
            // The retained negative accumulator is still exact zero here.
            out.loads[i]=rate!=0?-rate:sum.compare(scratch.negative)==0?0.L:-rate;
        }else {scratch.negative.add(sum,-1);scratch.negative.addQuotientProduct({u[i],material().cells()[i].integratedDefect},dt);
            out.loads[i]=round(scratch.negative);}
    };
    const auto row=[&](std::size_t i,RowAccumulator& scratch,HydroKepConvectionResult::Data& out){
        const auto& plan=*data_->openSupport;
        const auto degree=plan.pairOffsets[i+1]-plan.pairOffsets[i];
        if(useProducts&&plan.boundaryOffsets[i]==plan.boundaryOffsets[i+1]&&degree<scratch.productFactors.size()){
            const auto& compiled=*data_->openInteriorProducts;std::size_t count=0;
            scratch.productFactors[count]=compiled.diagonal[i];scratch.productRows[count++]=static_cast<std::uint32_t>(i);
            for(auto at=plan.pairOffsets[i];at<plan.pairOffsets[i+1];++at){const auto& item=plan.pairItems[at];
                if(data_->openPairCoefficients[item.id]==0)continue;const auto& pair=plan.pairs[item.id];
                scratch.productFactors[count]=item.sign>0?compiled.pairs[item.id]:compiled.pairs[item.id].negated();
                scratch.productRows[count++]=static_cast<std::uint32_t>(item.sign>0?pair.second:pair.first);}
            const auto exact=OpenProducts::sum(scratch.productFactors.data(),scratch.productRows.data(),count,productValues.get(),u.size());
            const auto rate=exact.sum.roundedDouble();out.rates[i]=rate;
            out.loads[i]=rate!=0?-rate:exact.sum.compare(Geometry::Detail::FiniteVolume::Dyadic{})==0?0.L:-rate;return;
        }
        reset(scratch.sum);reset(scratch.negative);auto& sum=scratch.sum;auto* budget=scratch.budget?&*scratch.budget:nullptr;
        for(auto at=plan.boundaryOffsets[i];at<plan.boundaryOffsets[i+1];++at){const auto id=plan.boundaryItems[at];const auto& q=gross[id];
            const auto ub=q.incomingModel==HydroVofPressureOutletMomentumModel::InteriorContinuationZeroShear?u[i]:q.incomingMomentum;
            addOpenMass(sum,q,material(),false,u[i]);addOpenMass(sum,q,material(),true,-ub);if(budget)budgetSegment(*budget,id);}
        if(budget&&plan.boundaryOffsets[i]!=plan.boundaryOffsets[i+1]){const int a=component(layout(),static_cast<int>(i));
            budget->value[a][0].add(sum); // Same actual boundary-only momentum, before diagonal/pairs.
            budget->value[a][5].addProduct({-.5L,data_->externalOutflow[i],u[i],u[i]});}
        diagonalTerms(sum,i);
        for(auto at=plan.pairOffsets[i];at<plan.pairOffsets[i+1];++at){const auto& item=plan.pairItems[at];const auto k=data_->openPairCoefficients[item.id];
            if(k==0)continue;const auto& pair=plan.pairs[item.id];pairTerm(sum,item.sign>0?k:-k,u[item.sign>0?pair.second:pair.first]);}
        storeRow(i,sum,scratch,out);
    };
    const auto finishBudget=[&](RowAccumulator* workers,std::size_t count,HydroKepConvectionResult::Data& out){
        if(!includeBoundaryBudget)return;auto& first=*workers[0].budget;
        for(int a=0;a<2;++a)for(int field=0;field<6;++field)
            for(std::size_t worker=1;worker<count;++worker)first.value[a][field].add(workers[worker].budget->value[a][field]);
        for(int a=0;a<2;++a){auto& target=out.boundaryBudget->components[a];const auto& value=first.value[a];
            target.boundaryMomentumRate=round(value[0]);target.boundaryKineticPower=round(value[1]);
            target.boundaryPhysicalKineticPower=round(value[2]);target.inletDissipationPower=round(value[3]);
            target.boundaryFluxRoundingPower=round(value[4]);target.boundaryDiagonalPower=round(value[5]);}
    };
    const auto literal=[&]{auto mark=OpenWall::now();auto out=makeResult();RowAccumulator scratch;
        if(includeBoundaryBudget)scratch.budget.emplace();
        // Without support, the legacy domain Exact scatter has no complete
        // stock report. Charge only the known new retained index there.
        receipt.scratchBytes=data_->openSupport?fixedProductStorage+sizeof(RowAccumulator)+(2u+(includeBoundaryBudget?19u:0u))*8192u+128u*1024u+(useProducts?128u*1024u:0u):retainedExterior;
        receipt.preparationSeconds=openSeconds(commonPreparationDuration+(OpenWall::now()-mark));mark=OpenWall::now();
        if(useProducts)for(std::size_t i=0;i<u.size();++i)productValues[i]=OpenProducts::decode(u[i]);
        if(data_->openSupport)for(std::size_t i=0;i<u.size();++i)row(i,scratch,*out);
        else {
            // The original linear exact scatter remains the bounded-route-failure
            // fallback; only this path retains a domain-sized Exact vector.
            std::vector<OpenExact> rates(u.size());
            for(std::size_t id=0;id<segments.size();++id)if(!internal(segments[id])){const int i=owner(segments[id]);const auto& q=gross[id];
                const auto ub=q.incomingModel==HydroVofPressureOutletMomentumModel::InteriorContinuationZeroShear?u[i]:q.incomingMomentum;
                addOpenMass(rates[i],q,material(),false,u[i]);addOpenMass(rates[i],q,material(),true,-ub);if(scratch.budget)budgetSegment(*scratch.budget,id);}
            for(std::size_t i=0;i<u.size();++i){if(scratch.budget&&data_->boundarySlot(static_cast<int>(i))>=0){const int a=component(layout(),static_cast<int>(i));
                    scratch.budget->value[a][0].add(rates[i]);scratch.budget->value[a][5].addProduct({-.5L,data_->externalOutflow[i],u[i],u[i]});}
                diagonalTerms(rates[i],i);}
            for(const auto& pair:data_->finalPairs){pairTerm(rates[pair.first],pair.coefficient,u[pair.second]);pairTerm(rates[pair.second],-pair.coefficient,u[pair.first]);}
            for(std::size_t i=0;i<u.size();++i){reset(scratch.negative);storeRow(i,rates[i],scratch,*out);}
        }
        finishBudget(&scratch,1,*out);receipt.rowSeconds=openSeconds(OpenWall::now()-mark);return out;};
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
    const int entryErrno=errno;OpenFloatingState entry;const bool qualified=openQualified(entry);
    receipt.parallelRequested=!omp_in_parallel()&&omp_get_max_threads()>1;
    if(entry.valid)entry.restore();errno=entryErrno;
    receipt.parallelEligible=qualified&&static_cast<bool>(data_->openSupport);
    if(receipt.parallelEligible){const auto optionalStart=OpenWall::now();
        try{const auto parallel=[&]{OpenCleanup cleanup;
            constexpr std::size_t rationalBytes=8192,limit=64u*1024u*1024u;
            // sizeof includes the disengaged optional's inline metadata. The
            // budget request constructs12 sums+7 exact temporaries per worker.
            const auto exactCount=2u+(includeBoundaryBudget?19u:0u);
            const auto workerBytes=sizeof(RowAccumulator)+OpenCooperativeTeam::workerStorageBytes()+exactCount*rationalBytes+128u*1024u+(useProducts?128u*1024u:0u);
            const auto fixed=fixedProductStorage+sizeof(OpenCooperativeTeam);
            if(fixed>limit||limit-fixed<2*workerBytes)throw OpenOptionalFailure{};
            const auto count=std::min<std::size_t>({static_cast<std::size_t>(omp_get_max_threads()),u.size(),(limit-fixed)/workerBytes});
            if(count<2)throw OpenOptionalFailure{};
            auto scratch=std::make_unique<RowAccumulator[]>(count);if(includeBoundaryBudget)for(std::size_t i=0;i<count;++i)scratch[i].budget.emplace();
            auto out=makeResult();OpenCooperativeTeam team(static_cast<int>(count));
            if(errno!=entryErrno)throw OpenOptionalFailure{};entry.restore();errno=entryErrno;receipt.scratchBytes=fixed+count*workerBytes;
            auto mark=OpenWall::now();receipt.preparationSeconds=openSeconds(commonPreparationDuration+(mark-optionalStart));
            return team.run([&]{
                if(useProducts&&!team.rows(u.size(),[&](std::size_t i){productValues[i]=OpenProducts::decode(u[i]);},receipt.contributionWorkers))throw OpenOptionalFailure{};
                if(!team.rows(u.size(),[&](std::size_t i){row(i,scratch[static_cast<std::size_t>(omp_get_thread_num())],*out);},receipt.contributionWorkers))throw OpenOptionalFailure{};
                finishBudget(scratch.get(),static_cast<std::size_t>(omp_get_num_threads()),*out);
                // The numerical row interval includes fused boundary budgets,
                // the one exact partial merge and the requested final6 RNs.
                receipt.rowSeconds=openSeconds(OpenWall::now()-mark);cleanup.arm();return out;});};
            auto out=parallel();receipt.parallelUsed=true;receipt.totalSeconds=openSeconds(OpenWall::now()-start);out->openApply=receipt;
            productsCleanup.arm();return HydroKepConvectionResult(std::move(out));
        }catch(...){const double failed=openSeconds(OpenWall::now()-optionalStart);
            receipt.preparationSeconds=receipt.rowSeconds=0;receipt.contributionWorkers=receipt.gatherWorkers=0;
            receipt.scratchBytes=0;receipt.parallelUsed=false;receipt.parallelFallback=true;receipt.fallbackSeconds=failed;useProducts=false;entry.restore();errno=entryErrno;}
    }
    try{
        const auto execute=[&]{
            if(!useProducts)return literal();const int error=errno;OpenFloatingState state;
            if(!state.capture()){useProducts=false;errno=error;return literal();}
            const auto attempt=OpenWall::now();
            try{return literal();}catch(...){receipt.fallbackSeconds+=openSeconds(OpenWall::now()-attempt);useProducts=false;
                receipt.preparationSeconds=receipt.rowSeconds=0;state.restore();errno=error;return literal();}
        };
        auto out=execute();receipt.totalSeconds=openSeconds(OpenWall::now()-start);out->openApply=receipt;productsCleanup.arm();
        return HydroKepConvectionResult(std::move(out));
    }catch(...){productsCleanup.arm();throw;}
#else
    auto out=literal();receipt.totalSeconds=openSeconds(OpenWall::now()-start);out->openApply=receipt;
    return HydroKepConvectionResult(std::move(out));
#endif
}

}
#undef product
