#include "Modules/Hydrodynamics/HydroVofMeshTransfer.h"
#include "Modules/Hydrodynamics/Detail/HydroRegularAdvectorGeometry.h"
#include "Core/Mesh/CartesianLeafLocator.h"
#include "Core/Mesh/NativeLeafOverlap.h"
#include "Core/Mesh/NativeMeshIdentity.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Geometry/Detail/FiniteVolumeClosure.h"
#include <algorithm>
#include <chrono>
#include <array>
#include <cfenv>
#include <cerrno>
#include <exception>
#if defined(_M_X64) || defined(__SSE2__)
#include <xmmintrin.h>
#endif
#if defined(_MSC_VER) && defined(_M_X64)
#include <float.h>
#endif
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

namespace LAMP::Hydrodynamics
{
namespace
{
using Clock=std::chrono::steady_clock;
// Passive observations only: never enclose a physical call or its exception.
class PassiveRemapSampler
{
    struct FloatingRestore {
        std::fenv_t native{};int error=errno;
#if defined(_M_X64) || defined(__SSE2__)
        unsigned csr=_mm_getcsr();
#endif
        bool captured=std::fegetenv(&native)==0;
        ~FloatingRestore() noexcept {
            if(captured&&std::fesetenv(&native)!=0)std::terminate();
#if defined(_M_X64) || defined(__SSE2__)
            _mm_setcsr(csr);
#endif
            errno=error;
        }
    };
public:
    struct Stamp {Clock::time_point value{};bool available=false;};
    bool available=true;
    Stamp mark() noexcept {
        FloatingRestore restore;
        if(!restore.captured){available=false;return {};}
        return {Clock::now(),true};
    }
    double seconds(Clock::time_point begin) noexcept {
        FloatingRestore restore;
        if(!restore.captured){available=false;return 0.;}
        const double value=std::chrono::duration<double>(Clock::now()-begin).count();
        if(!std::isfinite(value)||value<0){available=false;return 0.;}return value;
    }
    double seconds(Stamp begin) noexcept {
        if(!begin.available){available=false;return 0.;}return seconds(begin.value);
    }
    double unattributed(double total,const std::array<double,10>& children) noexcept {
        FloatingRestore restore;
        if(!restore.captured||!available){available=false;return 0.;}
        double value=total;for(double child:children)value-=child;
        if(!std::isfinite(value)){available=false;return 0.;}return value;
    }
};
void require(bool ok,const char* why){if(!ok)throw std::invalid_argument(std::string("[HydroVofMeshTransfer] ")+why);}
struct Sum{long double value=0,low=0;void add(long double x){const auto y=x-low,t=value+y;low=(t-value)-y;value=t;}};
double coordinate(Vector2D p,int axis){return axis==0?p.x:p.y;}
Vector2D at(int normal,double fixed,double tangent){return normal==0?Vector2D{fixed,tangent}:Vector2D{tangent,fixed};}
double elapsed(Clock::time_point begin){return std::chrono::duration<double>(Clock::now()-begin).count();}
double relative(double difference,double scale)
{
    require(std::isfinite(difference)&&std::isfinite(scale)&&scale>=0,"nonfinite remap divergence");
    require(scale>0||difference==0,"nonzero remap divergence at zero flux scale");return scale>0?std::abs(difference)/scale:0;
}
void divergence(const QuadtreeMesh& mesh,const HydroPhysicalDivergence& metrics,const HydroNativeValues& values,
    double& maximum,long double& l1)
{
    const auto actual=metrics.divergence().apply(mesh,values);const auto& d=metrics.divergence().matrix();
    for(int p=0;p<values.layout().pressureCount();++p){double scale=0;
        for(HydroIntegralDivergence::Matrix::InnerIterator it(d,p);it;++it)scale+=std::abs(it.value()*values.normalVelocity()[it.col()]);
        maximum=std::max(maximum,relative(actual[p],scale));l1+=std::abs(static_cast<long double>(actual[p]));}
}
HydroVofRemapNativeBudget budget(const HydroMomentumMass& mass,const std::vector<double>& velocity)
{
    HydroVofRemapNativeBudget result;
    for(int u=0;u<mass.layout().velocityCount();++u){const int component=mass.layout().normalAxis(u)==HydroDofLayout::NormalAxis::First?0:1;
        const long double m=mass.masses()[u],v=velocity[u];result.mass[component]+=m;result.momentum[component]+=m*v;result.kinetic+=.5L*m*v*v;}
    return result;
}
namespace FV=LAMP::Geometry::Detail::FiniteVolume;
using Mode=Interface::GeometricVofAdvectorSlab::Boundary::Mode;
struct FaceInterval {int normal;double fixed,lower,upper;};
FaceInterval interval(const QuadtreeMesh& mesh,const NativeLeafFace& face)
{
    const int normal=face.horizontal?1:0,tangent=1-normal;
    const auto& edge=face.horizontal?mesh.hEdges(face.edgeLabel):mesh.vEdges(face.edgeLabel);
    const auto a=mesh.verts(edge.edgeToVert[0]).coord(),b=mesh.verts(edge.edgeToVert[1]).coord();
    return {normal,coordinate(face.center,normal),std::min(coordinate(a,tangent),coordinate(b,tangent)),
        std::max(coordinate(a,tangent),coordinate(b,tangent))};
}
void samePhysicalDomain(const QuadtreeMesh& a,const NativeLeafDomain& ad,const QuadtreeMesh& b,const NativeLeafDomain& bd)
{
    ad.validateMeshIdentity(a);bd.validateMeshIdentity(b);
    require(ad.axisymmetric()==bd.axisymmetric()&&a.minCoord.x==b.minCoord.x&&a.minCoord.y==b.minCoord.y
        &&a.maxCoord.x==b.maxCoord.x&&a.maxCoord.y==b.maxCoord.y,"different actual physical remap domains");
    // A selected mask is a physical domain, not just a total volume. Reuse the
    // native overlap owner to certify every active source/target leaf and the
    // empty inactive groups before any boundary polynomial is consumed.
    if(ad.cellIndexByRow().size()!=static_cast<std::size_t>(a.cellIndexNum())||
        bd.cellIndexByRow().size()!=static_cast<std::size_t>(b.cellIndexNum()))
        (void)Mesh::NativeLeafOverlap::build(a,ad,b,bd);
}
double normalMean(const HydroTwoPhaseState& source,const QuadtreeMesh& targetMesh,const NativeLeafDomain& targetDomain,
    const NativeLeafFace& face,const Detail::HydroRegularAdvectorGeometry& regular,
    const Detail::HydroRegularCoefficients& coefficients,const CartesianLeafLocator& locator,HydroVofMeshTransferReport& report)
{
        require(face.area>0,"nonpositive physical interior face area");const int normal=face.horizontal?1:0,tangent=1-normal;
        const auto& edge=face.horizontal?targetMesh.hEdges(face.edgeLabel):targetMesh.vEdges(face.edgeLabel);
        const auto a=targetMesh.verts(edge.edgeToVert[0]).coord(),b=targetMesh.verts(edge.edgeToVert[1]).coord();
        const auto lo=std::min(coordinate(a,tangent),coordinate(b,tangent)),hi=std::max(coordinate(a,tangent),coordinate(b,tangent));
        const auto fixed=coordinate(face.center,normal);const auto begin=at(normal,fixed,lo),end=at(normal,fixed,hi);
        const auto candidates=locator.intersectingBox(begin,end);report.sourceLeafCandidates+=candidates.size();std::vector<double> knots{lo,hi};
        for(int native:candidates){const int p=source.values().layout().pressureDof(native);if(p<0)continue;
            const auto& cell=regular.cells()[p];if(fixed<coordinate(cell.lower,normal)||fixed>coordinate(cell.upper,normal))continue;
            const auto left=std::max(lo,coordinate(cell.lower,tangent)),right=std::min(hi,coordinate(cell.upper,tangent));if(!(left<right))continue;
            knots.push_back(left);knots.push_back(right);for(double knot:cell.supportKnots[tangent])if(left<knot&&knot<right)knots.push_back(knot);}
        std::sort(knots.begin(),knots.end());knots.erase(std::unique(knots.begin(),knots.end()),knots.end());Sum flux;
        for(std::size_t i=1;i<knots.size();++i){const auto left=knots[i-1],right=knots[i];if(!(left<right))continue;
            const auto center=at(normal,fixed,left+(right-left)/2);const auto incident=locator.incident(center);
            require(incident.count>0,"new native face leaves the original domain");int p=-1;
            // A wall/port span may also touch an inactive solid leaf. Retain
            // the original locator order among actual selected incident cells.
            for(int j=0;j<incident.count&&p<0;++j)p=source.values().layout().pressureDof(incident.indices[j]);
            require(p>=0,"new face span has no actual original pressure leaf");
            const auto polynomial=regular.normalPolynomial(coefficients,p,normal,at(normal,fixed,left),at(normal,fixed,right));
            // This bool is fromEnd (chart), NOT a measure selector. integral()
            // uses the actual radialMeasure/radialTangent and already applies
            // 2*pi*r in RZ. Never apply another radial weight to this flux.
            require(polynomial.radialMeasure==targetDomain.axisymmetric(),"different actual polynomial physical measure");
            flux.add(polynomial.integral(false,0,polynomial.length));++report.faceSpans;++report.polynomialIntegrals;}
        const auto value=flux.value/face.area;const auto mean=static_cast<double>(value);
        require(std::isfinite(value)&&std::isfinite(mean)&&(value==0||mean!=0),"unrepresentable actual face-mean remap");
        return mean;
}
// Roles are topology/physical data, independent of the projected solution.
// 0=interior, 1=wall, 2=fixed reservoir normal, 3=free pressure outlet.
std::vector<int> roles(const QuadtreeMesh& mesh,const HydroDofLayout& layout,const HydroVofPrescribedFlow& flow)
{
    require(flow.normalVelocity.size()==static_cast<std::size_t>(layout.velocityCount())&&
        flow.rulesByNativeRow.size()==static_cast<std::size_t>(layout.domain().nativeFaceCount()),
        "different target physical boundary dimensions/data");
    if(flow.pressureDrivenReservoirs){
        require(flow.prescribedVelocityDofs.empty()&&!flow.endpointVelocity&&flow.pressureOutlets.size()>=2
            &&flow.pressureOutletBackflowModel.has_value(),"invalid pressure-driven remap boundary data");
        bool distinct=false;for(const auto& outlet:flow.pressureOutlets)
            distinct=distinct||outlet.pressure!=flow.pressureOutlets.front().pressure;
        require(distinct,"pressure-driven remap requires distinct physical pressures");
        for(double value:flow.normalVelocity)require(value==0,"pressure-driven remap has a fixed normal payload");
    }else require(bool(flow.endpointVelocity),"target fixed flow has no endpoint velocity data");
    std::vector<int> result(layout.velocityCount());
    for(int u=0;u<layout.velocityCount();++u){const auto f=layout.domain().face(mesh,layout.nativeFace(u));
        require(std::isfinite(flow.normalVelocity[u]),"nonfinite target normal payload");
        result[u]=f.outwardSign!=0?1:0;
        if(f.outwardSign==0)require(flow.rulesByNativeRow[layout.nativeFace(u)]==Mode::Interior,"interior target row is not Interior");}
    for(int u:flow.prescribedVelocityDofs){require(u>=0&&u<layout.velocityCount()&&result[u]==1,"invalid/duplicate target fixed port");
        const auto f=layout.domain().face(mesh,layout.nativeFace(u));
        require(f.area>0&&flow.normalVelocity[u]!=0,"target fixed port is axis or zero normal");result[u]=2;}
    for(const auto& o:flow.pressureOutlets){const int u=o.velocityDof;
        require(u>=0&&u<layout.velocityCount()&&result[u]==1&&std::isfinite(o.pressure),"invalid/overlapping target pressure outlet");
        require(layout.domain().face(mesh,layout.nativeFace(u)).area>0&&flow.normalVelocity[u]==0,"target free outlet has fixed payload/zero area");result[u]=3;}
    if(!flow.pressureDrivenReservoirs)require(!flow.prescribedVelocityDofs.empty(),"target has no fixed physical inlet");
    for(int u=0;u<layout.velocityCount();++u){const auto mode=flow.rulesByNativeRow[layout.nativeFace(u)];
        if(result[u]==2||result[u]==3)require(mode==Mode::PureLiquidReservoir||mode==Mode::PureGasReservoir,"target open port has no phase reservoir");
        else require(flow.normalVelocity[u]==0,"target payload outside fixed port");
        if(result[u]==1)require(mode==Mode::Closed,"target unlisted wall is not Closed");}
    require(std::is_sorted(flow.stationaryNoSlipVelocityDofs.begin(),flow.stationaryNoSlipVelocityDofs.end())
        &&std::adjacent_find(flow.stationaryNoSlipVelocityDofs.begin(),flow.stationaryNoSlipVelocityDofs.end())
            ==flow.stationaryNoSlipVelocityDofs.end(),"remap stationary no-slip wall rows must be sorted and unique");
    for(int u:flow.stationaryNoSlipVelocityDofs)
        require(u>=0&&u<layout.velocityCount()&&result[u]==1
            &&layout.domain().face(mesh,layout.nativeFace(u)).area>0
            &&flow.rulesByNativeRow[layout.nativeFace(u)]==Mode::Closed,
            "remap stationary no-slip requires a positive-area Closed wall, not port/interior/axis");
    return result;
}
double outletPressure(const HydroVofPrescribedFlow& flow,int u)
{
    const auto it=std::find_if(flow.pressureOutlets.begin(),flow.pressureOutlets.end(),[&](const auto& o){return o.velocityDof==u;});
    require(it!=flow.pressureOutlets.end(),"missing actual physical outlet pressure");return it->pressure;
}
Vector2D ambientMomentum(const HydroVofPrescribedFlow& flow,int row,LAMP::TimeLayer fullMid)
{
    require(bool(flow.pressureOutletBackflowMomentum),"missing actual ambient momentum query on remap");
    NativeMeshIdentity::ValidationBatch::Pause callback;
    const auto first=flow.pressureOutletBackflowMomentum(row,fullMid),repeat=flow.pressureOutletBackflowMomentum(row,fullMid);
    require(std::isfinite(first.x)&&std::isfinite(first.y)&&first.x==repeat.x&&first.y==repeat.y,
        "nonfinite/nonrepeatable actual ambient momentum on remap");return first;
}
// roles() already queried every source face in this same immutable batch.
// On the callback-free strict path, repeating an interior face can only OR
// the same sticky FP flags; it has no later numerical or ownership effect.
bool reuseBoundarySourceRoles(const HydroVofPrescribedFlow& source,const HydroVofPrescribedFlow& target) noexcept
{
#if defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE) && defined(_M_FP_STRICT) && _M_FP_STRICT \
    && !defined(__FAST_MATH__) && (!defined(_M_FP_FAST) || !_M_FP_FAST)
    using Model=HydroVofPressureOutletMomentumModel;
    if(!source.pressureDrivenReservoirs||!target.pressureDrivenReservoirs
        ||source.endpointVelocity||target.endpointVelocity
        ||source.pressureOutletBackflowMomentum||target.pressureOutletBackflowMomentum
        ||source.pressureOutletBackflowModel!=Model::InteriorContinuationZeroShear
        ||target.pressureOutletBackflowModel!=Model::InteriorContinuationZeroShear)return false;
    const int error=errno;unsigned control=0;const unsigned csr=_mm_getcsr();
    bool qualified=_controlfp_s(&control,0,0)==0&&(control&_MCW_EM)==_MCW_EM
        &&(control&_MCW_RC)==_RC_NEAR&&std::fegetround()==FE_TONEAREST
        &&(csr&0x1f80u)==0x1f80u&&(csr&0xe040u)==0;
#ifdef _EM_AMBIGUOUS
    qualified=qualified&&(control&_EM_AMBIGUOUS)==0;
#endif
    errno=error;return qualified&&sizeof(double)==8&&sizeof(long double)==8
        &&std::numeric_limits<double>::is_iec559&&std::numeric_limits<double>::digits==53
        &&std::numeric_limits<long double>::digits==53;
#else
    (void)source;(void)target;return false;
#endif
}
std::vector<int> authenticateBoundary(const QuadtreeMesh& sourceMesh,const HydroTwoPhaseState& source,
    const QuadtreeMesh& targetMesh,const HydroPhaseGeometry& phase,const HydroVofPrescribedFlow& targetFlow,
    const std::vector<double>& seed,HydroVofMeshTransferReport& report)
{
    const auto& old=*source.vofPrescribedFlow();const auto& sl=source.values().layout();const auto& tl=phase.layout();
    samePhysicalDomain(sourceMesh,sl.domain(),targetMesh,tl.domain());
    require(old.pressureDrivenReservoirs==targetFlow.pressureDrivenReservoirs,
        "target changes the actual pressure-driven reservoir mode");
    require(old.pressureOutletBackflowModel==targetFlow.pressureOutletBackflowModel,
        "target changes the actual pressure-outlet backflow model");
    bool ambientBackflow=false;auto fullMid=source.timeLayer();
    if(old.pressureOutletBackflowModel){using Model=HydroVofPressureOutletMomentumModel;
        const auto model=*old.pressureOutletBackflowModel;
        require(model==Model::Ambient||model==Model::InteriorContinuationZeroShear,"unknown remap backflow model");
        require(!old.pressureOutlets.empty()&&!targetFlow.pressureOutlets.empty(),"remap backflow has no actual pressure ports");
        ambientBackflow=model==Model::Ambient;
        if(ambientBackflow){
            require(bool(old.pressureOutletBackflowMomentum)&&bool(targetFlow.pressureOutletBackflowMomentum),"ambient remap requires actual source/target queries");
            require(fullMid.position==LAMP::TimeLayer::Position::End,"ambient remap requires the actual accepted End interval");
            fullMid.position=LAMP::TimeLayer::Position::Midpoint;fullMid.time();
        }else require(!old.pressureOutletBackflowMomentum&&!targetFlow.pressureOutletBackflowMomentum,"interior continuation cannot consume ambient remap data");
    }else require(!old.pressureOutletBackflowMomentum&&!targetFlow.pressureOutletBackflowMomentum,"legacy remap cannot consume an unused ambient query");
    const auto boundary=Interface::GeometricVofAdvectorSlab::Boundary::prepare(targetMesh,phase.vofPublication(),targetFlow.rulesByNativeRow);
    require(!boundary.allClosed(),"target physical flow is all Closed");const auto sr=roles(sourceMesh,sl,old),tr=roles(targetMesh,tl,targetFlow);
    const bool skipInterior=reuseBoundarySourceRoles(old,targetFlow);
    std::vector<FV::Dyadic> covered(sl.velocityCount());
    for(int u=0;u<tl.velocityCount();++u){const auto tf=tl.domain().face(targetMesh,tl.nativeFace(u));if(tf.outwardSign==0)continue;
        const auto ti=interval(targetMesh,tf);FV::Dyadic length;Sum tangential,tangentWeight;
        bool hasTangent=false,constantTangent=true;double firstTangent=0;
        for(int v=0;v<sl.velocityCount();++v){if(skipInterior&&sr[v]==0)continue;
            const auto sf=sl.domain().face(sourceMesh,sl.nativeFace(v));if(sf.outwardSign==0)continue;
            const auto si=interval(sourceMesh,sf);if(si.normal!=ti.normal||si.fixed!=ti.fixed||sf.outwardSign!=tf.outwardSign)continue;
            const double left=std::max(si.lower,ti.lower),right=std::min(si.upper,ti.upper);if(!(left<right))continue;
            FV::Dyadic span(right);span.add(left,-1);length.add(span);covered[v].add(span);++report.boundaryCoverageSpans;
            require(sr[v]==tr[u]&&old.rulesByNativeRow[sl.nativeFace(v)]==targetFlow.rulesByNativeRow[tl.nativeFace(u)],
                "target boundary merges/changes source physical role or phase reservoir");
            require(std::binary_search(old.stationaryNoSlipVelocityDofs.begin(),old.stationaryNoSlipVelocityDofs.end(),v)
                ==std::binary_search(targetFlow.stationaryNoSlipVelocityDofs.begin(),targetFlow.stationaryNoSlipVelocityDofs.end(),u),
                "target boundary adds/removes a physical stationary no-slip wall span");
            if(tr[u]==3){require(outletPressure(old,v)==outletPressure(targetFlow,u),"target merges/changes physical outlet pressure");
                if(ambientBackflow){
                    const auto from=ambientMomentum(old,sl.nativeFace(v),fullMid),to=ambientMomentum(targetFlow,tl.nativeFace(u),fullMid);
                    require(from.x==to.x&&from.y==to.y,"target changes/merges distinct actual P0 ambient momentum vectors");
                }}
            if(tr[u]==2){require((sf.outwardSign*old.normalVelocity[v]<0)==(tf.outwardSign*targetFlow.normalVelocity[u]<0),
                    "target fixed port merges opposite physical flow directions");
                if(tf.outwardSign*targetFlow.normalVelocity[u]<0){const auto q=[&]{
                    NativeMeshIdentity::ValidationBatch::Pause callback;
                    return old.endpointVelocity(sl.nativeFace(v),source.timeLayer());}();
                    require(std::isfinite(q.x)&&std::isfinite(q.y),"nonfinite actual source inlet vector");
                    const auto area=CartesianMeasure::face(right-left,ti.normal==1?(left+(right-left)/2):ti.fixed,tl.domain().axisymmetric());
                    const double value=coordinate(q,1-ti.normal);
                    if(!hasTangent){firstTangent=value;hasTangent=true;}else constantTangent=constantTangent&&value==firstTangent;
                    tangentWeight.add(area);tangential.add(static_cast<long double>(area)*value);}}
        }
        FV::Dyadic whole(ti.upper);whole.add(ti.lower,-1);require(length.compare(whole)==0,"target boundary lacks exact source native coverage");
        if(tr[u]==2){++report.fixedNormals;require(targetFlow.normalVelocity[u]==seed[u],"target fixed normal differs from actual source polynomial integral");
            if(tf.outwardSign*targetFlow.normalVelocity[u]<0){const auto q=[&]{
                    NativeMeshIdentity::ValidationBatch::Pause callback;
                    return targetFlow.endpointVelocity(tl.nativeFace(u),source.timeLayer());}();
                const auto repeat=[&]{NativeMeshIdentity::ValidationBatch::Pause callback;
                    return targetFlow.endpointVelocity(tl.nativeFace(u),source.timeLayer());}();
                require(std::isfinite(q.x)&&std::isfinite(q.y)&&q.x==repeat.x&&q.y==repeat.y,"nonfinite/nonrepeatable target inlet vector");
                require(hasTangent&&tangentWeight.value>0,"missing actual source inlet tangent measure");
                // A constant declared P0 value has this exact weighted mean;
                // do not alter it through nonadditive rounded RZ subareas.
                const double tangent=constantTangent?firstTangent:static_cast<double>(tangential.value/tangentWeight.value);
                require(std::isfinite(tangent)&&coordinate(q,ti.normal)==seed[u]&&coordinate(q,1-ti.normal)==tangent,
                    "target inlet vector differs from same-End source P0 boundary mean");}}
        else if(tr[u]==1){++report.wallNormals;require(seed[u]==0,"nonzero actual source normal on target wall");
            if(tf.area==0){++report.axisNormals;require(targetFlow.rulesByNativeRow[tl.nativeFace(u)]==Mode::Closed,"target axis is not Closed");}}
        else ++report.pressureOutletNormals;
    }
    for(int v=0;v<sl.velocityCount();++v){const auto sf=sl.domain().face(sourceMesh,sl.nativeFace(v));if(sf.outwardSign==0)continue;
        const auto si=interval(sourceMesh,sf);FV::Dyadic whole(si.upper);whole.add(si.lower,-1);
        require(covered[v].compare(whole)==0,"source physical boundary loses exact target native coverage");}
    return tr;
}
void openProjectionEnergy(const HydroMomentumMass& mass,const HydroPressureLoad& pressure,const std::vector<int>& free,
    const std::vector<double>& b,const HydroNativeValues& raw,const HydroNativeValues& values,
    const std::vector<double>& lambda,HydroVofMeshTransferReport& report)
{
    // Full contractions in exact native-input arithmetic, once rounded only
    // after each complete report quantity. No residual is repaired/subtracted.
    std::vector<bool> isFree(b.size());for(int u:free)isFree[u]=true;
    FV::Dyadic delta,correction,wb,wrhs,freeWork,continuity,scale;
    const auto& c=pressure.matrix();
    for(int u=0;u<mass.layout().velocityCount();++u){const long double m=mass.masses()[u],v=raw.normalVelocity()[u],x=values.normalVelocity()[u];
        delta.add(FV::Dyadic::product({.5L,m,x,x}));delta.add(FV::Dyadic::product({.5L,m,v,v}),-1);
        correction.add(FV::Dyadic::product({.5L,m,x,x}));correction.add(FV::Dyadic::product({.5L,m,v,v}));
        correction.add(FV::Dyadic::product({m,x,v}),-1);
        FV::Dyadic rhsWork=FV::Dyadic::product({x,b[u]});rhsWork.add(FV::Dyadic::product({x,m,v}),-1);wrhs.add(rhsWork);
        FV::Dyadic rowWork=FV::Dyadic::product({x,m,x});rowWork.add(FV::Dyadic::product({x,b[u]}),-1);
        for(HydroPressureLoad::Matrix::InnerIterator it(c,u);it;++it)rowWork.add(FV::Dyadic::product({x,it.value(),lambda[it.col()]}));
        if(isFree[u])freeWork.add(rowWork);else wb.add(rowWork);
    }
    // A separate physical scale, independent of the equation residual.
    scale=FV::Dyadic{};
    for(int u=0;u<mass.layout().velocityCount();++u){const long double m=mass.masses()[u],v=raw.normalVelocity()[u],x=values.normalVelocity()[u];
        scale.add(FV::Dyadic::product({m,x,x}).absolute());scale.add(FV::Dyadic::product({m,v,v}).absolute());
        scale.add(FV::Dyadic::product({m,x,v}).absolute());scale.add(FV::Dyadic::product({x,b[u]}).absolute());
        for(HydroPressureLoad::Matrix::InnerIterator it(c,u);it;++it)scale.add(FV::Dyadic::product({x,it.value(),lambda[it.col()]}).absolute());}
    const auto& d=pressure.divergence().matrix();
    for(int p=0;p<static_cast<int>(lambda.size());++p)for(HydroIntegralDivergence::Matrix::InnerIterator it(d,p);it;++it)
        continuity.add(FV::Dyadic::product({lambda[p],it.value(),values.normalVelocity()[it.col()]}));
    FV::Dyadic identity=delta;identity.add(correction);identity.add(wb,-1);identity.add(wrhs,-1);
    report.projectionKineticChange=delta.rounded();report.projectionCorrectionEnergy=correction.rounded();
    report.projectionBoundaryWork=wb.rounded();report.projectionRhsRoundingWork=wrhs.rounded();
    report.projectionFreeResidualWork=freeWork.rounded();report.projectionContinuityWork=continuity.rounded();
    report.projectionEnergyIdentityResidual=identity.rounded();report.projectionEnergyIdentityScale=scale.rounded();
    require(scale.compare(FV::Dyadic{})>0||identity.compare(FV::Dyadic{})==0,"nonzero open projection energy at zero physical scale");
    report.projectionEnergyIdentityRelativeResidual=report.projectionEnergyIdentityScale>0?
        std::abs(report.projectionEnergyIdentityResidual)/report.projectionEnergyIdentityScale:0;
    require(report.projectionEnergyIdentityRelativeResidual<=2e-11L,"open projection constraint-energy identity failed");
}

}
struct HydroVofMeshTransfer::Data
{
    Data(HydroTwoPhaseState s,Interface::GeometricVofMeshTransfer t,HydroNativeValues r,std::vector<double> multiplier,
        HydroTwoPhaseState c,HydroVofMeshTransferReport report,std::vector<double> momentumRhs)
        :source(std::move(s)),interfaceTransfer(std::move(t)),raw(std::move(r)),multiplier(std::move(multiplier)),candidate(std::move(c)),report(std::move(report)),momentumRhs(std::move(momentumRhs)){}
    HydroTwoPhaseState source;Interface::GeometricVofMeshTransfer interfaceTransfer;HydroNativeValues raw;
    std::vector<double> multiplier;HydroTwoPhaseState candidate;HydroVofMeshTransferReport report;std::vector<double> momentumRhs;
};
HydroVofMeshTransfer::HydroVofMeshTransfer(std::shared_ptr<const Data> d):data_(std::move(d)){}
HydroVofMeshTransfer HydroVofMeshTransfer::prepare(const QuadtreeMesh& sourceMesh,const HydroTwoPhaseState& source,
    const QuadtreeMesh& targetMesh,Interface::GeometricVofMeshTransfer transfer,double rhoL,double rhoG,const HydroVofMeshTransferConfig& config)
{return prepareImpl(sourceMesh,source,targetMesh,std::move(transfer),rhoL,rhoG,config,nullptr);}
HydroVofMeshTransfer HydroVofMeshTransfer::preparePrescribedFlow(const QuadtreeMesh& sourceMesh,const HydroTwoPhaseState& source,
    const QuadtreeMesh& targetMesh,Interface::GeometricVofMeshTransfer transfer,double rhoL,double rhoG,
    HydroVofPrescribedFlow targetFlow,const HydroVofMeshTransferConfig& config)
{return prepareImpl(sourceMesh,source,targetMesh,std::move(transfer),rhoL,rhoG,config,&targetFlow);}
std::vector<double> HydroVofMeshTransfer::boundaryNormalMeans(const QuadtreeMesh& sourceMesh,const HydroTwoPhaseState& source,
    const QuadtreeMesh& targetMesh,const NativeLeafDomain& targetDomain)
{
    {NativeMeshIdentity::ValidationBatch::Pause state;source.validateFor(sourceMesh);}
    NativeMeshIdentity::ValidationBatch sourceValidation(sourceMesh),targetValidation(targetMesh);
    require(source.model()==HydroPhaseGeometry::Model::GeometricVof&&source.vofPrescribedFlow(),
        "boundary means require actual open VOF source");
    require(source.timeLayer().position==LAMP::TimeLayer::Position::End,"boundary means require accepted End");
    samePhysicalDomain(sourceMesh,source.values().layout().domain(),targetMesh,targetDomain);
    const auto& regular=source.reconstruction().regularGeometry();
    const auto coefficients=regular.apply(sourceMesh,source.values(),&source.boundary());const CartesianLeafLocator locator(sourceMesh);
    std::vector<double> result(targetDomain.nativeFaceCount());HydroVofMeshTransferReport report;
    for(int row=0;row<targetDomain.nativeFaceCount();++row){const auto f=targetDomain.face(targetMesh,row);if(f.outwardSign==0||f.area==0)continue;
        result[row]=normalMean(source,targetMesh,targetDomain,f,regular,coefficients,locator,report);}
    return result;
}
HydroVofMeshTransfer HydroVofMeshTransfer::prepareImpl(const QuadtreeMesh& sourceMesh,const HydroTwoPhaseState& source,
    const QuadtreeMesh& targetMesh,Interface::GeometricVofMeshTransfer transfer,double rhoL,double rhoG,
    const HydroVofMeshTransferConfig& config,const HydroVofPrescribedFlow* targetFlow)
{
    PassiveRemapSampler sampler;HydroVofRemapPreparationTimings preparation;
    const auto start=Clock::now();
    {NativeMeshIdentity::ValidationBatch::Pause state;source.validateFor(sourceMesh);}
    // Each immutable mesh has its own proof chain; opaque boundary queries
    // below suspend and clear both before entering caller-owned code.
    NativeMeshIdentity::ValidationBatch sourceValidation(sourceMesh),targetValidation(targetMesh);
    require(source.model()==HydroPhaseGeometry::Model::GeometricVof&&!source.prescribedFlow(),"requires a genuine closed VOF accepted source");
    if(targetFlow)require(source.vofPrescribedFlow()!=nullptr,"explicit open remap requires actual open source");
    else require(!source.vofPrescribedFlow(),"open VOF remap requires an explicit physical boundary projection");
    const auto clock=source.timeLayer();require(clock.position==LAMP::TimeLayer::Position::End,"remap requires the same accepted End clock");
    require(config.velocityTransfer==HydroVofMeshVelocityTransfer::RegularFluxThenNativeMassProjection,"invalid preselected velocity remap law");
    require(std::isfinite(rhoL)&&rhoL>0&&std::isfinite(rhoG)&&rhoG>0,"invalid physical remap densities");
    require(transfer.sourceIdentity()==source.phaseGeometry().vofIdentity(),"Interface transfer has a different actual accepted source");
    const auto& overlaps=transfer.nativeOverlaps();overlaps.validateFor(sourceMesh,targetMesh);
    require(source.phaseGeometry().vofPublication().geometryKind()==Interface::GeometricVofGeometryKind::LocalPlic,
        "requires the actual accepted LocalPlic source");
    const auto& target=transfer.candidate();target.validateFor(targetMesh,target.identity(),clock);
    require(target.geometryKind()==Interface::GeometricVofGeometryKind::LocalPlic,"requires LocalPlic transfer candidate");
    require(std::any_of(overlaps.pieces().begin(),overlaps.pieces().end(),[](const auto& p){
        return p.kind!=Mesh::NativeLeafOverlapPiece::Kind::Unchanged;}),"no-op tree must skip fluid remap");
    require(target.domain().axisymmetric()==source.phaseGeometry().layout().domain().axisymmetric(),"different source/target physical measure");
    preparation.inputValidationSeconds=sampler.seconds(start);const auto layoutStarted=sampler.mark();
    const auto layout=HydroDofLayout::build(targetMesh,target.domain());
    const auto metrics=HydroPhysicalDivergence::build(targetMesh,HydroIntegralDivergence::build(targetMesh,layout));
    preparation.targetLayoutMetricSeconds=sampler.seconds(layoutStarted);const auto geometryStarted=sampler.mark();
    const auto geometry=HydroMomentumGeometry::buildWithRemappedRegularPreparation(sourceMesh,source.reconstruction().regularGeometry(),targetMesh,metrics,overlaps);const auto phase=HydroPhaseGeometry::build(targetMesh,metrics,target);
    preparation.targetMomentumPhaseSeconds=sampler.seconds(geometryStarted);
    using Density=HydroMomentumMass::DensityRepresentation;
    require(config.momentumDensityRepresentation==Density::SharpGeometry||config.momentumDensityRepresentation==Density::CellAverageProjection,
        "invalid selected remap density representation");
    const auto makeMass=[&](const QuadtreeMesh& mesh,const HydroMomentumGeometry& g,const HydroPhaseGeometry& p){
        return config.momentumDensityRepresentation==Density::CellAverageProjection
            ?HydroMomentumMass::buildCellAverageProjection(mesh,g,p,rhoL,rhoG)
            :HydroMomentumMass::build(mesh,g,p,rhoL,rhoG);};
    const auto targetMassStarted=sampler.mark();
    const auto mass=makeMass(targetMesh,geometry,phase);
    preparation.targetMassSeconds=sampler.seconds(targetMassStarted);const auto sourceMassStarted=sampler.mark();
    const auto sourceMass=[&]{if(const auto* retained=source.retainedMass()){
        require(retained->liquidDensity()==rhoL&&retained->gasDensity()==rhoG
            &&retained->representation()==config.momentumDensityRepresentation,"remap cannot change accepted densities/representation");
        retained->validateVofFor(sourceMesh,source.phaseGeometry().vofIdentity(),clock);return *retained;}
        return makeMass(sourceMesh,source.momentumGeometry(),source.phaseGeometry());}();
    preparation.sourceMassSeconds=sampler.seconds(sourceMassStarted);
    HydroVofMeshTransferReport report;report.velocityTransfer=config.velocityTransfer;
    report.momentumDensityRepresentation=config.momentumDensityRepresentation;report.prescribedFlow=targetFlow!=nullptr;
    const auto fluxStart=Clock::now();const auto& regular=source.reconstruction().regularGeometry();
    const auto coefficients=regular.apply(sourceMesh,source.values(),&source.boundary());const CartesianLeafLocator locator(sourceMesh);
    std::vector<double> velocity(layout.velocityCount()),pressure(layout.pressureCount());
    for(int u=0;u<layout.velocityCount();++u){
        const auto face=layout.domain().face(targetMesh,layout.nativeFace(u));
        if(!targetFlow&&(face.negativeCell<0||face.positiveCell<0)){++report.wallNormals;if(face.area==0)++report.axisNormals;continue;}
        if(face.area==0){require(face.outwardSign!=0,"zero-area interior face");continue;}
        velocity[u]=normalMean(source,targetMesh,layout.domain(),face,regular,coefficients,locator,report);
    }
    if(targetFlow)authenticateBoundary(sourceMesh,source,targetMesh,phase,*targetFlow,velocity,report);
    // Physical pressure is a cell scalar: preserve its overlap-volume average
    // before resetting the one complete-domain volume gauge. The remap saddle
    // multiplier below never replaces this physical pressure.
    for(int p=0;p<layout.pressureCount();++p){const int targetCell=layout.nativeCell(p);Sum integral;
        for(std::size_t id:overlaps.byTarget()[targetCell]){const auto& piece=overlaps.pieces()[id];const int oldP=source.values().layout().pressureDof(piece.sourceCell);
            require(oldP>=0,"pressure overlap has no original native source");integral.add(piece.capacity*source.values().pressure()[oldP]);}
        pressure[p]=static_cast<double>(integral.value/metrics.cellVolumes()[p]);require(std::isfinite(pressure[p]),"unrepresentable physical pressure average");}
    Sum pressureIntegral,totalVolume;for(int p=0;p<layout.pressureCount();++p){pressureIntegral.add(static_cast<long double>(metrics.cellVolumes()[p])*pressure[p]);totalVolume.add(metrics.cellVolumes()[p]);}
    require(totalVolume.value>0,"empty pressure-gauge volume");report.pressureIntegralBeforeGauge=pressureIntegral.value;
    report.pressureGaugeApplied=!targetFlow||targetFlow->pressureOutlets.empty();
    if(report.pressureGaugeApplied){
    report.pressureGaugeShift=pressureIntegral.value/totalVolume.value;
    for(double& value:pressure){value=static_cast<double>(value-report.pressureGaugeShift);require(std::isfinite(value),"nonfinite physical pressure gauge");}
    Sum newMean,newGross;for(int p=0;p<layout.pressureCount();++p){const auto value=static_cast<long double>(metrics.cellVolumes()[p])*pressure[p];newMean.add(value);newGross.add(std::abs(value));}
    require(newGross.value>0||newMean.value==0,"nonzero physical pressure mean at zero scale");
    report.pressureGaugeRelativeResidual=newGross.value>0?std::abs(newMean.value)/newGross.value:0;
    require(report.pressureGaugeRelativeResidual<=1e-11L,"physical pressure volume gauge failed");
    }
    const HydroNativeValues raw(layout,velocity,pressure);divergence(targetMesh,metrics,raw,report.maximumRawRelativeDivergence,report.rawDivergenceL1);
    report.rawFluxSeconds=elapsed(fluxStart);
    const auto projectionStart=Clock::now();HydroSaddleVelocityBlock::Matrix diagonal(layout.velocityCount(),layout.velocityCount());std::vector<double> rhs(velocity.size()),prescribed(velocity.size());
    for(int u=0;u<layout.velocityCount();++u){diagonal.insert(u,u)=mass.masses()[u];const auto value=static_cast<long double>(mass.masses()[u])*velocity[u];rhs[u]=static_cast<double>(value);
        require(std::isfinite(rhs[u])&&(value==0||rhs[u]!=0),"unrepresentable original mass-projection RHS");}
    diagonal.makeCompressed();
    std::vector<HydroSaddlePressureOutlet> projectionOutlets;
    if(targetFlow){for(int u:targetFlow->prescribedVelocityDofs)prescribed[u]=targetFlow->normalVelocity[u];
        for(const auto& outlet:targetFlow->pressureOutlets)projectionOutlets.push_back({outlet.velocityDof,0});}
    const auto pressureLoad=HydroPressureLoad::build(metrics.divergence());
    const auto saddle=projectionOutlets.empty()?HydroSinglePhaseSaddle::buildClosed(pressureLoad,HydroSaddleVelocityBlock::build(layout,diagonal),config.saddle):
        HydroSinglePhaseSaddle::buildWithPressureOutlets(pressureLoad,HydroSaddleVelocityBlock::build(layout,diagonal),projectionOutlets,config.saddle);
    const auto projected=saddle.solve(targetMesh,HydroSaddleRhs(layout,rhs,prescribed));report.remapSolves=1;report.projectionSaddle=projected.report();
    report.projectionSeconds=elapsed(projectionStart);
    const auto postProjectionStarted=sampler.mark();
    const HydroNativeValues values(layout,projected.values().normalVelocity(),pressure);
    divergence(targetMesh,metrics,values,report.maximumFinalRelativeDivergence,report.finalDivergenceL1);
    require(report.maximumFinalRelativeDivergence<=1e-11,"projected native remap is not divergence free");
    preparation.postProjectionValidationSeconds=sampler.seconds(postProjectionStarted);const auto targetStateStarted=sampler.mark();
    // Reuse the actual target geometry/M, retaining one endpoint mass owner.
    const auto candidate=[&]{NativeMeshIdentity::ValidationBatch::Pause state;
        return targetFlow?HydroTwoPhaseState::createVofPrescribedFlowRemapped(targetMesh,phase,values,mass,*targetFlow):
            HydroTwoPhaseState::createVofRemapped(targetMesh,phase,values,mass);}();
    preparation.targetStateSeconds=sampler.seconds(targetStateStarted);const auto budgetStarted=sampler.mark();
    const auto nextClock=candidate.timeLayer();
    require(candidate.phaseGeometry().vofIdentity()==target.identity()&&nextClock.step==clock.step
        &&nextClock.begin==clock.begin&&nextClock.end==clock.end&&nextClock.position==clock.position,
        "different actual accepted remap source/clock");
    report.historyReset=candidate.velocityHistory()==nullptr;require(report.historyReset,"remap retained incompatible old-layout velocity history");
    report.source=budget(sourceMass,source.values().normalVelocity());report.rawTarget=budget(mass,raw.normalVelocity());report.target=budget(mass,values.normalVelocity());
    for(int c=0;c<2;++c){report.momentumChange[c]=report.target.momentum[c]-report.source.momentum[c];report.projectionMomentumChange[c]=report.target.momentum[c]-report.rawTarget.momentum[c];}
    report.kineticChange=report.target.kinetic-report.source.kinetic;report.projectionKineticChange=report.target.kinetic-report.rawTarget.kinetic;
    for(int u=0;u<layout.velocityCount();++u){const long double difference=values.normalVelocity()[u]-static_cast<long double>(raw.normalVelocity()[u]);report.projectionCorrectionEnergy+=.5L*mass.masses()[u]*difference*difference;}
    report.projectionEnergyIdentityResidual=report.projectionKineticChange+report.projectionCorrectionEnergy;
    if(targetFlow)openProjectionEnergy(mass,pressureLoad,saddle.freeVelocityDofs(),rhs,raw,values,projected.values().pressure(),report);
    preparation.budgetValidationSeconds=sampler.seconds(budgetStarted);
    report.totalSeconds=elapsed(start);
    preparation.unattributedSeconds=sampler.unattributed(report.totalSeconds,{{
        preparation.inputValidationSeconds,preparation.targetLayoutMetricSeconds,preparation.targetMomentumPhaseSeconds,
        preparation.targetMassSeconds,preparation.sourceMassSeconds,report.rawFluxSeconds,report.projectionSeconds,
        preparation.postProjectionValidationSeconds,preparation.targetStateSeconds,preparation.budgetValidationSeconds}});
    preparation.available=sampler.available;report.preparation=preparation;
    return HydroVofMeshTransfer(std::make_shared<Data>(source,std::move(transfer),raw,projected.values().pressure(),candidate,std::move(report),std::move(rhs)));
}
const std::vector<double>& HydroVofMeshTransfer::projectionMomentumRhs()const noexcept{return data_->momentumRhs;}
const HydroTwoPhaseState& HydroVofMeshTransfer::sourceState()const noexcept{return data_->source;}
const Interface::GeometricVofMeshTransfer& HydroVofMeshTransfer::interfaceTransfer()const noexcept{return data_->interfaceTransfer;}
const HydroNativeValues& HydroVofMeshTransfer::rawValues()const noexcept{return data_->raw;}
const std::vector<double>& HydroVofMeshTransfer::projectionMultiplier()const noexcept{return data_->multiplier;}
const HydroTwoPhaseState& HydroVofMeshTransfer::candidateState()const noexcept{return data_->candidate;}
const HydroVofMeshTransferReport& HydroVofMeshTransfer::report()const noexcept{return data_->report;}
}
