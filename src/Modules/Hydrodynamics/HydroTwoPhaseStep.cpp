#include "Modules/Hydrodynamics/HydroTwoPhaseStep.h"
#include "Modules/Hydrodynamics/HydroFrozenMidpoint.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/UniformLeafGrid.h"
#include "Modules/Interface/HermiteGalsAmrClosure.h"
#include "Modules/Hydrodynamics/Detail/HydroRegularAdvectorGeometry.h"
#include "Core/Geometry/Detail/ExactExpansion.h"
#include "Core/IO/BinaryArchive.h"
#include "Core/Mesh/FaceSegmentGraph.h"
#include "Core/Mesh/NativeMeshIdentity.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <chrono>
#include <cstring>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>

namespace LAMP::Hydrodynamics
{
namespace
{
void require(bool ok,const char* why)
{if(!ok)throw std::invalid_argument(std::string("[HydroTwoPhaseStep] ")+why);}
double relative(double error,double scale)
{
    require(std::isfinite(error)&&std::isfinite(scale)&&error>=0&&scale>=0,"nonfinite update norm");
    if(scale==0){require(error==0,"nonzero update at exactly zero field scale");return 0;}
    const double result=error/scale;require(std::isfinite(result),"nonfinite relative update");return result;
}
double difference(const std::vector<double>& a,const std::vector<double>& b)
{
    require(a.size()==b.size(),"different compared field sizes");double error=0,na=0,nb=0;
    for(std::size_t i=0;i<a.size();++i){require(std::isfinite(a[i])&&std::isfinite(b[i]),"nonfinite compared field");error=std::hypot(error,a[i]-b[i]);na=std::hypot(na,a[i]);nb=std::hypot(nb,b[i]);}
    return relative(error,std::max(na,nb));
}
std::vector<double> geometryValues(const Interface::HermiteGalsPublication& p)
{
    const auto& primary=p.state().primary();const auto& closure=p.state().closure();std::vector<double> values;
    values.reserve(3*primary.vertices.size()+2*primary.edges.size()+closure.vertexMixedDerivatives.size()+closure.cellBubbleCoefficients.size());
    for(const auto& v:primary.vertices){values.push_back(v.phi);values.push_back(v.dx);values.push_back(v.dy);}
    for(const auto& e:primary.edges){values.push_back(e.phi);values.push_back(e.dn);}
    values.insert(values.end(),closure.vertexMixedDerivatives.begin(),closure.vertexMixedDerivatives.end());
    values.insert(values.end(),closure.cellBubbleCoefficients.begin(),closure.cellBubbleCoefficients.end());return values;
}
double geometryDifference(const Interface::HermiteGalsPublication& a,const Interface::HermiteGalsPublication& b)
{
    const auto& x=a.state().primary();const auto& y=b.state().primary();
    require(x.transport==y.transport&&x.measure==y.measure&&x.geometryCells==y.geometryCells&&x.vertices.size()==y.vertices.size()&&x.edges.size()==y.edges.size(),"different FirstJet primary layouts");
    for(std::size_t i=0;i<x.vertices.size();++i)require(x.vertices[i].vertex==y.vertices[i].vertex,"different FirstJet vertex layout");
    for(std::size_t i=0;i<x.edges.size();++i)require(x.edges[i].vertices==y.edges[i].vertices,"different FirstJet edge layout");
    return difference(geometryValues(a),geometryValues(b));
}
void supportedGeometry(const QuadtreeMesh& mesh,const HydroPhaseGeometry& phase)
{
    require(phase.model()==HydroPhaseGeometry::Model::FirstJet,"FirstJet state factory requires its genuine model");
    const auto& l=phase.layout();l.validateMesh(mesh);phase.validateFor(mesh,phase.geometryIdentity(),phase.timeLayer());const auto& pub=phase.publication();
    require(l.domain().uniformLeaves(),"two-phase step requires uniform XY or RZ leaves");
    require(l.pressureCount()==mesh.cellIndexNum()&&pub.domain().cellIndexByRow()==l.domain().cellIndexByRow(),"two-phase step requires the complete native domain");
    const auto& p=pub.state().primary();
    const auto measure=l.domain().axisymmetric()?Interface::HermiteGalsMeasure::AxisymmetricVolume
        :Interface::HermiteGalsMeasure::CartesianArea;
    require(p.transport==Interface::HermiteGalsTransport::FirstJet&&p.phaseVolumes.empty()&&p.measure==measure,
        "two-phase step requires pure FirstJet with the native XY/RZ measure and no inventory");
    require(p.geometryCells.size()==static_cast<std::size_t>(mesh.cellIndexNum()),"FirstJet scalar support is not complete");
    const auto clock=phase.timeLayer();clock.time();using P=Interface::HermiteGalsTimeLayer::Position;
    require(clock.position==P::Initial||clock.position==P::End,"accepted state must be Initial or End");
}
template<class CheckLayout,class CheckPhase>
void supportedVofGeometryImpl(const QuadtreeMesh& mesh,const HydroPhaseGeometry& phase,
    CheckLayout&& checkLayout,CheckPhase&& checkPhase)
{
    require(phase.model()==HydroPhaseGeometry::Model::GeometricVof,"VOF state factory requires its genuine model");
    const auto& layout=phase.layout();checkLayout();
    checkPhase();
    const auto& source=phase.vofPublication().domain();
    require(static_cast<std::size_t>(layout.pressureCount())==source.cellIndexByRow().size()
        &&source.cellIndexByRow()==layout.domain().cellIndexByRow()
        &&source.cellRowByIndex()==layout.domain().cellRowByIndex()
        &&source.meshIdentity()==layout.domain().meshIdentity(),"VOF state requires the actual selected same-source native domain");
    const auto clock=phase.timeLayer();clock.time();using P=Interface::GeometricVofTimeLayer::Position;
    require(clock.step>=0&&(clock.position==P::Initial||clock.position==P::End),"VOF accepted state must be Initial or End");
}
void supportedVofGeometry(const QuadtreeMesh& mesh,const HydroPhaseGeometry& phase)
{
    supportedVofGeometryImpl(mesh,phase,[&]{phase.layout().validateMesh(mesh);},
        [&]{phase.validateVofFor(mesh,phase.vofIdentity(),phase.timeLayer());});
}
std::vector<int> closedFaces(const QuadtreeMesh& mesh,const HydroDofLayout& l)
{
    std::vector<int> walls;for(int u=0;u<l.velocityCount();++u){const auto f=l.domain().face(mesh,l.nativeFace(u));if(f.negativeCell<0||f.positiveCell<0)walls.push_back(u);}require(!walls.empty(),"closed domain has no boundary ports");return walls;
}
template<class CheckValues,class ApplyDivergence>
void acceptedValuesImpl(const QuadtreeMesh&,const HydroPhaseGeometry& phase,const HydroNativeValues& values,
    const std::vector<int>& walls,bool pressureAnchored,CheckValues&& checkValues,ApplyDivergence&& applyDivergence)
{
    checkValues();for(int u:walls)require(values.normalVelocity()[u]==0,"nonzero closed boundary normal velocity");
    const auto& d=phase.metrics().divergence().matrix();const auto flow=applyDivergence();
    for(int p=0;p<phase.layout().pressureCount();++p){double scale=0;for(HydroIntegralDivergence::Matrix::InnerIterator it(d,p);it;++it)scale+=std::abs(it.value()*values.normalVelocity()[it.col()]);require(relative(std::abs(flow[p]),scale)<=1e-11,"accepted native velocity is not divergence free");}
    // The supported complete uniform rectangle is one connected pressure
    // component. A true traction outlet anchors it; do not shift that pressure.
    if(pressureAnchored)return;
    long double mean=0,gross=0;for(int p=0;p<phase.layout().pressureCount();++p){const auto t=static_cast<long double>(phase.metrics().cellVolumes()[p])*values.pressure()[p];mean+=t;gross+=std::abs(t);}require(relative(static_cast<double>(std::abs(mean)),static_cast<double>(gross))<=1e-11,"accepted pressure gauge is not mean zero");
}
void acceptedValues(const QuadtreeMesh& mesh,const HydroPhaseGeometry& phase,const HydroNativeValues& values,
    const std::vector<int>& walls,bool pressureAnchored=false)
{
    acceptedValuesImpl(mesh,phase,values,walls,pressureAnchored,[&]{values.validateFor(phase.layout(),mesh);},
        [&]{return phase.metrics().divergence().apply(mesh,values);});
}
std::vector<int> prescribedWalls(const QuadtreeMesh& mesh,const HydroPhaseGeometry& phase,
    const HydroNativeValues& values,const HydroTwoPhasePrescribedFlow& flow)
{
    const auto& l=phase.layout();values.validateFor(l,mesh);
    require(flow.normalVelocity.size()==static_cast<std::size_t>(l.velocityCount()),"different prescribed-flow normal dimensions");
    require(!flow.inletVelocityDofs.empty()&&(!flow.outletVelocityDofs.empty()||!flow.pressureOutlets.empty())
        &&bool(flow.inletVelocity)&&bool(flow.inletFirstJet),"prescribed flow requires inlet/outlet and complete inlet data");
    std::vector<int> direction(l.velocityCount());
    for(int sign:{-1,1})for(int u:sign<0?flow.inletVelocityDofs:flow.outletVelocityDofs){
        require(u>=0&&u<l.velocityCount()&&direction[u]==0,"duplicate or invalid prescribed-flow port");
        const auto f=l.domain().face(mesh,l.nativeFace(u));
        require(f.outwardSign!=0&&phase.metrics().faceMeasures()[u]>0,"prescribed-flow port is internal or zero-area axis");
        require(std::isfinite(flow.normalVelocity[u])&&sign*f.outwardSign*flow.normalVelocity[u]>0,"prescribed-flow port has wrong direction");
        require(values.normalVelocity()[u]==flow.normalVelocity[u],"accepted prescribed-flow normal mean differs");direction[u]=sign;
    }
    for(const auto& outlet:flow.pressureOutlets){
        const int u=outlet.velocityDof;
        require(u>=0&&u<l.velocityCount()&&direction[u]==0,"duplicate, overlapping or invalid pressure outlet");
        const auto f=l.domain().face(mesh,l.nativeFace(u));
        require(f.outwardSign!=0&&phase.metrics().faceMeasures()[u]>0,"pressure outlet is internal or zero-area axis");
        require(std::isfinite(outlet.pressure),"nonfinite pressure outlet traction");
        require(flow.normalVelocity[u]==0,"pressure outlet prescribed normal payload must be zero");
        require(f.outwardSign*values.normalVelocity()[u]>0,"accepted pressure outlet reverses or has zero flow");
        direction[u]=2;
    }
    std::vector<int> walls;long double net=0,gross=0;
    for(int u=0;u<l.velocityCount();++u){require(std::isfinite(flow.normalVelocity[u]),"nonfinite prescribed-flow normal mean");
        if(direction[u]==0)require(flow.normalVelocity[u]==0,"prescribed-flow payload on a wall or free face");
        const auto f=l.domain().face(mesh,l.nativeFace(u));
        if(f.outwardSign!=0){if(direction[u]==0)walls.push_back(u);
            const auto normal=flow.pressureOutlets.empty()?flow.normalVelocity[u]:values.normalVelocity()[u];
            const auto q=static_cast<long double>(f.outwardSign)*phase.metrics().faceMeasures()[u]*normal;net+=q;gross+=std::abs(q);}
    }
    require(relative(static_cast<double>(std::abs(net)),static_cast<double>(gross))<=1e-12,"prescribed-flow boundary flux compatibility failed");
    return walls;
}
// VOF physical ports never reuse a FirstJet arrival/extension callback.
std::vector<int> vofPrescribedWalls(const QuadtreeMesh& mesh,const HydroPhaseGeometry& phase,
    const HydroNativeValues& values,const HydroVofPrescribedFlow& flow)
{
    using Mode=Interface::GeometricVofAdvectorSlab::Boundary::Mode;
    const auto& l=phase.layout();values.validateFor(l,mesh);
    require(flow.normalVelocity.size()==static_cast<std::size_t>(l.velocityCount()),"different VOF prescribed normal dimensions");
    if(flow.pressureDrivenReservoirs){
        require(flow.prescribedVelocityDofs.empty()&&!bool(flow.endpointVelocity),
            "pressure-driven VOF reservoirs cannot prescribe endpoint velocity");
        require(flow.pressureOutlets.size()>=2,"pressure-driven VOF reservoirs require multiple free pressure ports");
        require(flow.pressureOutletBackflowModel.has_value(),
            "pressure-driven VOF reservoirs require an explicit backflow momentum model");
    }else require(!flow.prescribedVelocityDofs.empty()&&bool(flow.endpointVelocity),"VOF prescribed flow requires real exterior velocity data");
    const auto sourceBoundary=Interface::GeometricVofAdvectorSlab::Boundary::prepare(mesh,phase.vofPublication(),flow.rulesByNativeRow);
    require(!sourceBoundary.allClosed(),"VOF prescribed flow has no open reservoir");
    const bool openMomentum=flow.pressureOutletBackflowModel.has_value();
    require(openMomentum||!bool(flow.pressureOutletBackflowMomentum),
        "VOF backflow callback requires an explicit momentum model");
    if(openMomentum){
        require(*flow.pressureOutletBackflowModel==HydroVofPressureOutletMomentumModel::Ambient||
            *flow.pressureOutletBackflowModel==HydroVofPressureOutletMomentumModel::InteriorContinuationZeroShear,
            "unknown VOF pressure-outlet momentum model");
        require(!flow.pressureOutlets.empty(),"VOF open momentum requires actual pressure ports");
        if(*flow.pressureOutletBackflowModel==HydroVofPressureOutletMomentumModel::Ambient)
            require(bool(flow.pressureOutletBackflowMomentum),"VOF ambient backflow momentum requires an explicit query");
        else require(!bool(flow.pressureOutletBackflowMomentum),
            "VOF interior continuation cannot consume ambient callback data");
    }
    std::vector<int> direction(l.velocityCount());bool incoming=false,outgoing=false;
    bool distinctPressure=false;
    const auto reservoir=[&](int u){const auto mode=flow.rulesByNativeRow[l.nativeFace(u)];
        return mode==Mode::PureLiquidReservoir||mode==Mode::PureGasReservoir;};
    for(const int u:flow.prescribedVelocityDofs){
        require(u>=0&&u<l.velocityCount()&&direction[u]==0,"duplicate or invalid VOF prescribed port");
        const auto face=l.domain().face(mesh,l.nativeFace(u));
        require(face.outwardSign!=0&&phase.metrics().faceMeasures()[u]>0&&reservoir(u),"VOF prescribed port requires positive-area exterior reservoir");
        const auto normal=flow.normalVelocity[u];require(std::isfinite(normal)&&normal!=0,"invalid steady VOF prescribed normal");
        require(values.normalVelocity()[u]==normal,"accepted VOF prescribed normal differs from steady source");
        const auto outward=face.outwardSign*normal;incoming=incoming||outward<0;outgoing=outgoing||outward>0;
        if(outward<0){const auto vector=flow.endpointVelocity(l.nativeFace(u),phase.timeLayer());
            require(std::isfinite(vector.x)&&std::isfinite(vector.y),"nonfinite VOF exterior velocity vector");
            require((face.horizontal?vector.y:vector.x)==normal,"VOF exterior vector normal differs from native mean");}
        direction[u]=1;
    }
    for(const auto& port:flow.pressureOutlets){const int u=port.velocityDof;
        require(u>=0&&u<l.velocityCount()&&direction[u]==0,"duplicate or overlapping VOF pressure outlet");
        const auto face=l.domain().face(mesh,l.nativeFace(u));
        require(face.outwardSign!=0&&phase.metrics().faceMeasures()[u]>0&&reservoir(u)&&std::isfinite(port.pressure),
            "VOF pressure outlet requires positive-area exterior reservoir");
        require(flow.normalVelocity[u]==0,"free VOF pressure outlet has prescribed normal payload");
        if(flow.pressureDrivenReservoirs)distinctPressure=distinctPressure||port.pressure!=flow.pressureOutlets.front().pressure;
        if(openMomentum){
            require(std::isfinite(values.normalVelocity()[u]),"nonfinite accepted VOF free pressure normal");
            // Free-port direction belongs to the actual accepted velocity, not
            // to the reservoir's phase rule or the sign of its pressure data.
            if(flow.pressureDrivenReservoirs)incoming=incoming||face.outwardSign*values.normalVelocity()[u]<0;
            outgoing=outgoing||face.outwardSign*values.normalVelocity()[u]>0;
        }else{
            require(face.outwardSign*values.normalVelocity()[u]>0,"accepted VOF pressure outlet reverses or has zero flow");
            outgoing=true;
        }
        direction[u]=2;
    }
    if(flow.pressureDrivenReservoirs){
        require(distinctPressure,"pressure-driven VOF reservoirs require distinct physical pressures");
        // A resting accepted state has no incoming/outgoing flux yet. The
        // original Saddle traction and actual VOF net/gross transport solve
        // those directions; no artificial seed velocity or phase stock enters.
    }else require(incoming&&outgoing,"VOF prescribed flow requires actual incoming and outgoing ports");
    std::vector<int> walls;long double net=0,gross=0;
    for(int u=0;u<l.velocityCount();++u){require(std::isfinite(flow.normalVelocity[u]),"nonfinite VOF prescribed normal payload");
        if(direction[u]==0)require(flow.normalVelocity[u]==0,"VOF prescribed normal payload outside a prescribed port");
        const auto face=l.domain().face(mesh,l.nativeFace(u));
        if(face.outwardSign==0)continue;
        if(direction[u]==0){require(flow.rulesByNativeRow[l.nativeFace(u)]==Mode::Closed,"unlisted VOF exterior is not Closed");walls.push_back(u);}
        const auto q=static_cast<long double>(face.outwardSign)*phase.metrics().faceMeasures()[u]*values.normalVelocity()[u];net+=q;gross+=std::abs(q);
    }
    require(std::is_sorted(flow.stationaryNoSlipVelocityDofs.begin(),flow.stationaryNoSlipVelocityDofs.end())
        &&std::adjacent_find(flow.stationaryNoSlipVelocityDofs.begin(),flow.stationaryNoSlipVelocityDofs.end())
            ==flow.stationaryNoSlipVelocityDofs.end(),"VOF stationary no-slip wall rows must be sorted and unique");
    for(const int u:flow.stationaryNoSlipVelocityDofs){
        require(u>=0&&u<l.velocityCount()&&std::binary_search(walls.begin(),walls.end(),u)
            &&phase.metrics().faceMeasures()[u]>0,"VOF stationary no-slip requires a positive-area Closed wall, not port/interior/axis");
    }
    require(relative(static_cast<double>(std::abs(net)),static_cast<double>(gross))<=1e-12,"VOF physical boundary flux compatibility failed");
    return walls;
}
// Both endpoint quartics have one strict sign over the ENTIRE physical atom.
// Their linear temporal interpolation then cannot hide counterflow in net F.
void prescribedDirections(const QuadtreeMesh& mesh,const HydroCharacteristicSource& source,
    const HydroTwoPhasePrescribedFlow& flow)
{
    const auto& r=source.geometry().reconstruction();const auto& l=r.momentumGeometry().layout();
    std::vector<int> direction(l.velocityCount());for(int u:flow.inletVelocityDofs)direction[u]=-1;for(int u:flow.outletVelocityDofs)direction[u]=1;
    for(const auto& outlet:flow.pressureOutlets)direction[outlet.velocityDof]=1;
    constexpr long double choose[5][5]={{1,0,0,0,0},{1,1,0,0,0},{1,2,1,0,0},{1,3,3,1,0},{1,4,6,4,1}};
    for(int u:closedFaces(mesh,l)){const auto f=l.domain().face(mesh,l.nativeFace(u));
        for(int id:r.primalAtoms()[u])for(const auto* endpoint:{&source.previous(),&source.arrival()}){
            const auto& polynomial=endpoint->atomPolynomials()[id];const auto& c=polynomial.beginCoefficients;
            if(direction[u]==0){for(auto a:c)require(a==0,"nonzero actual impermeable source polynomial");continue;}
            std::array<long double,5> power{};long double lengthPower=1;
            for(int k=0;k<5;++k){power[k]=c[k]*lengthPower;lengthPower*=polynomial.length;}
            for(int j=0;j<5;++j){long double value=0;for(int k=0;k<=j;++k)value+=choose[j][k]/choose[4][k]*power[k];
                require(std::isfinite(value)&&direction[u]*f.outwardSign*value>0,"prescribed-flow source reverses or cannot certify its atom direction");}
        }
    }
}
HydroConvectionBoundary inletBoundary(const QuadtreeMesh& mesh,const HydroDualMaterial& dual,
    const HydroTwoPhasePrescribedFlow& flow)
{
    const auto& l=dual.layout();const auto& segments=dual.geometry().segments();std::vector<double> traces(segments.size());
    for(std::size_t g=0;g<segments.size();++g){const auto& s=segments[g];if(s.negativeVelocity>=0&&s.positiveVelocity>=0)continue;
        const auto q=(s.negativeVelocity>=0?1:-1)*dual.faceRates()[g];if(q>=0)continue;
        require(s.primalVelocity>=0&&std::find(flow.inletVelocityDofs.begin(),flow.inletVelocityDofs.end(),s.primalVelocity)!=flow.inletVelocityDofs.end(),"incoming dual flux lacks a declared inlet");
        const auto point=s.begin+(s.end-s.begin)*.5;const auto side=s.normalAxis==0
            ?(s.negativeVelocity>=0?Boundary::east:Boundary::west):(s.negativeVelocity>=0?Boundary::north:Boundary::south);
        const auto v=flow.inletVelocity(point,dual.timeLayer(),side);
        LAMP::Geometry::Detail::ExactExpansion::validateActualEnvironment();
        require(std::isfinite(v.x)&&std::isfinite(v.y),"nonfinite prescribed inlet velocity");
        const auto normal=s.normalAxis==0?v.x:v.y;const auto given=flow.normalVelocity[s.primalVelocity];
        require(relative(std::abs(normal-given),std::abs(normal)+std::abs(given))<=1e-11,"inlet velocity disagrees with prescribed native normal mean");
        const int u=s.negativeVelocity>=0?s.negativeVelocity:s.positiveVelocity;
        traces[g]=l.normalAxis(u)==HydroDofLayout::NormalAxis::First?v.x:v.y;
    }
    return HydroConvectionBoundary::bind(mesh,dual,traces,dual.sources(),dual.timeLayer());
}
HydroConvectionBoundary prescribedBoundary(const QuadtreeMesh& mesh,const HydroDualMaterial& dual,
    const HydroConvectionBoundary& inlet,const HydroNativeValues& values)
{
    inlet.validateFor(mesh,dual,dual.sources(),dual.timeLayer());auto traces=inlet.traces();
    const auto& segments=dual.geometry().segments();
    for(std::size_t g=0;g<segments.size();++g){const auto& s=segments[g];if(s.negativeVelocity>=0&&s.positiveVelocity>=0)continue;
        const auto q=(s.negativeVelocity>=0?1:-1)*dual.faceRates()[g];
        if(q>0){const int u=s.negativeVelocity>=0?s.negativeVelocity:s.positiveVelocity;traces[g]=values.normalVelocity()[u];}
        else if(q==0)traces[g]=0;
    }
    return HydroConvectionBoundary::bind(mesh,dual,traces,dual.sources(),dual.timeLayer());
}
HydroConvectionBoundary convectionBoundary(const QuadtreeMesh& mesh,const HydroDualMaterial& dual,const HydroSpatialAdvector& actual)
{
    const auto& l=dual.layout();const auto& segments=dual.geometry().segments();std::vector<double> traces(segments.size());
    for(std::size_t i=0;i<segments.size();++i){const auto& s=segments[i];if(s.negativeVelocity>=0&&s.positiveVelocity>=0)continue;const int u=s.negativeVelocity>=0?s.negativeVelocity:s.positiveVelocity;require(u>=0,"exterior dual segment lacks native velocity owner");int p=s.ownerPressure;
        if(p<0){require(s.primalVelocity>=0,"exterior dual segment lacks physical pressure owner");const auto f=l.domain().face(mesh,l.nativeFace(s.primalVelocity));const int native=f.negativeCell>=0?f.negativeCell:f.positiveCell;p=l.pressureDof(native);}require(p>=0,"exterior trace pressure is unselected");
        const auto point=s.begin+(s.end-s.begin)*.5;const auto jet=actual.cellVelocityJet(p,point).jet.value;traces[i]=l.normalAxis(u)==HydroDofLayout::NormalAxis::First?jet.x:jet.y;
    }return HydroConvectionBoundary::bind(mesh,dual,traces,dual.sources(),dual.timeLayer());
}
HydroNativeValues midpoint(const HydroDofLayout& l,const HydroNativeValues& a,const HydroNativeValues& b)
{
    std::vector<double> u(a.normalVelocity().size());for(std::size_t i=0;i<u.size();++i)u[i]=.5*a.normalVelocity()[i]+.5*b.normalVelocity()[i];return {l,u,b.pressure()};
}

long double interfaceMeasure(const Interface::HermiteGalsPublication& publication)
{
    long double area=0;
    for(std::size_t c=0;c<publication.topology().curves.size();++c)
        for(std::size_t g=0;g<publication.topology().curves[c].graphs().size();++g)
            for(const auto& q:publication.quadrature(c,g)){
                require(std::isfinite(q.physicalWeight)&&q.physicalWeight>0,
                    "nonpositive physical interface quadrature weight");
                area+=static_cast<long double>(q.physicalWeight);
            }
    require(std::isfinite(area)&&area>=0,"nonfinite physical interface measure");return area;
}

double interfaceMismatch(const QuadtreeMesh& mesh,const Interface::HermiteGalsPublication& publication,
    const HydroSpatialAdvector& predicted,const HydroSpatialAdvector& solved,double dt)
{
    const auto& layout=predicted.values().layout();
    const auto spacing=layout.domain().spacing();const double h=std::min(spacing.x,spacing.y);
    require(std::isfinite(h)&&h>0,"invalid interface mismatch grid scale");double peak=0;
    for(std::size_t c=0;c<publication.topology().curves.size();++c)
        for(std::size_t g=0;g<publication.topology().curves[c].graphs().size();++g)
            for(const auto& q:publication.quadrature(c,g)){
                const auto& frame=q.sample.frame;const int row=layout.pressureDof(frame.nativeCellIndex);
                require(row>=0,"interface mismatch sample is outside the native layout");
                const auto a=predicted.cellVelocityJet(row,frame.point).jet.value;
                const auto b=solved.cellVelocityJet(row,frame.point).jet.value;
                const double value=std::abs((b.x-a.x)*frame.normal.x+(b.y-a.y)*frame.normal.y);
                require(std::isfinite(value),"nonfinite interface velocity mismatch");peak=std::max(peak,value);
            }
    const double eta=(dt/h)*peak;require(std::isfinite(eta),"nonfinite interface mismatch");return eta;
}
double correctionRelative(const HydroNativeValues& previous,const HydroNativeValues& solved,
    const HydroMomentumMass& mass,double referenceVelocity)
{
    const auto& a=previous.normalVelocity();const auto& b=solved.normalVelocity();
    require(a.size()==b.size()&&a.size()==mass.masses().size(),"different correction metric sizes");
    long double error=0,norm=0,total=0;
    for(std::size_t i=0;i<a.size();++i){const long double m=mass.masses()[i],d=static_cast<long double>(b[i])-a[i];
        require(std::isfinite(m)&&m>0,"nonpositive correction metric mass");error+=m*d*d;norm+=m*b[i]*b[i];total+=m;}
    const long double scale=std::max(std::sqrt(norm),static_cast<long double>(referenceVelocity)*std::sqrt(total));
    require(std::isfinite(error)&&std::isfinite(scale)&&scale>0,"nonfinite correction metric");
    const double result=static_cast<double>(std::sqrt(error)/scale);require(std::isfinite(result),"nonfinite correction relative norm");return result;
}
bool recoverableGeometry(const std::string& message)
{
    // Explicit numerical failures only; source/config/observer errors propagate.
    for(const char* text:{"spacetime refinement depth exhausted","spacetime geometry sample budget exhausted",
        "unresolved spacetime center","unresolved spacetime panel","unresolved atom-local cumulative error",
        "interface residual exceeds","nonunit/nonfinite interface normal","nonfinite meridional curvature",
        "nonfinite physical curvature","unresolved interface topology"})
        if(message.find(text)!=std::string::npos)return true;
    return false;
}
struct Trial
{
    Trial(HydroNativeValues guess,HydroCharacteristicSource source,HydroMaterialTransport transport,HydroDualMaterial dual,
        HydroKepConvection convection,HydroKepConvectionResult loads,HydroPhaseViscosity viscosity,HydroSinglePhaseSaddle saddle,
        HydroSaddleResult solve,HydroSaddleVelocityBlock::Matrix matrix,std::vector<double> rhs,
        std::vector<long double> bodyLoads,std::optional<HydroBalancedCsfLoad> capillary,
        std::optional<HydroInterfaceNormalLoad> sharp,std::vector<double> sharpRhs,HydroFrozenMidpoint frozen)
        :guess(std::move(guess)),source(std::move(source)),transport(std::move(transport)),dual(std::move(dual)),convection(std::move(convection)),loads(std::move(loads)),viscosity(std::move(viscosity)),saddle(std::move(saddle)),solve(std::move(solve)),matrix(std::move(matrix)),rhs(std::move(rhs)),bodyLoads(std::move(bodyLoads)),capillary(std::move(capillary)),sharp(std::move(sharp)),sharpRhs(std::move(sharpRhs)),frozen(std::move(frozen)){}
    HydroNativeValues guess;HydroCharacteristicSource source;HydroMaterialTransport transport;HydroDualMaterial dual;
    HydroKepConvection convection;HydroKepConvectionResult loads;HydroPhaseViscosity viscosity;
    HydroSinglePhaseSaddle saddle;HydroSaddleResult solve;HydroSaddleVelocityBlock::Matrix matrix;std::vector<double> rhs;
    std::vector<long double> bodyLoads;
    std::optional<HydroBalancedCsfLoad> capillary;
    std::optional<HydroInterfaceNormalLoad> sharp;
    std::vector<double> sharpRhs; // Derived sign conversion of the retained sharp balance owner.
    HydroFrozenMidpoint frozen; // The actual shared NS preparation, not another solver path.
};
}

struct HydroVofSteadyVelocity::Data
{
    Data(HydroDofLayout layout,std::vector<Vector2D> values,std::vector<bool> covered,ClockRole role)
        :layout(std::move(layout)),values(std::move(values)),covered(std::move(covered)),role(role){}
    HydroDofLayout layout;
    std::vector<Vector2D> values;
    std::vector<bool> covered;
    ClockRole role;
};
HydroVofSteadyVelocity::HydroVofSteadyVelocity(std::shared_ptr<const Data> data):data_(std::move(data)){}
HydroVofSteadyVelocity HydroVofSteadyVelocity::prepare(const QuadtreeMesh& mesh,const HydroDofLayout& layout,
    std::vector<Vector2D> values,std::vector<bool> covered,ClockRole role)
{
    require(role==ClockRole::Endpoint||role==ClockRole::FullMidpoint,"unknown steady VOF velocity clock role");
    layout.validateMesh(mesh);const auto& domain=layout.domain();
    require(domain.meshIdentity()!=nullptr,"steady VOF velocity has no native mesh identity");
    const auto count=static_cast<std::size_t>(domain.nativeFaceCount());
    require(values.size()==count&&covered.size()==count,"steady VOF velocity requires complete native-row arrays");
    for(std::size_t row=0;row<count;++row){
        require(std::isfinite(values[row].x)&&std::isfinite(values[row].y),"nonfinite steady VOF native velocity vector");
        if(!covered[row])continue;
        const int velocity=layout.velocityDof(static_cast<int>(row));
        require(velocity>=0&&layout.nativeFace(velocity)==static_cast<int>(row),
            "covered steady VOF native row has no inverse Hydro velocity DOF");
        const auto face=domain.face(mesh,static_cast<int>(row));
        require((face.negativeCell>=0)!=(face.positiveCell>=0),"steady VOF velocity coverage requires a physical boundary row");
        require(face.outwardSign!=0&&std::isfinite(face.area)&&face.area>0,
            "steady VOF velocity coverage requires a positive-area physical boundary");
    }
    return HydroVofSteadyVelocity(std::make_shared<Data>(layout,std::move(values),std::move(covered),role));
}
Vector2D HydroVofSteadyVelocity::operator()(int nativeRow,LAMP::TimeLayer clock) const
{
    require(bool(data_),"consumed steady VOF velocity");clock.time();
    using Position=LAMP::TimeLayer::Position;
    require(data_->role==ClockRole::Endpoint?(clock.position==Position::Initial||clock.position==Position::End):
        data_->role==ClockRole::FullMidpoint&&clock.position==Position::Midpoint,
        "steady VOF velocity query has a different clock role");
    require(nativeRow>=0&&static_cast<std::size_t>(nativeRow)<data_->covered.size()&&data_->covered[nativeRow],
        "steady VOF velocity queried outside its covered native rows");
    return data_->values[nativeRow];
}
bool HydroVofSteadyVelocity::matches(const HydroDofLayout& layout,ClockRole role) const noexcept
{return data_&&data_->role==role&&data_->layout.sameIdentity(layout);}
void HydroVofSteadyVelocity::validateFor(const QuadtreeMesh& mesh,const HydroDofLayout& layout,ClockRole role) const
{
    require(matches(layout,role),"steady VOF velocity has a different layout or clock role");
    layout.validateMesh(mesh);
}
HydroVofSteadyVelocity HydroVofSteadyVelocity::reboundFor(const QuadtreeMesh& mesh,const HydroDofLayout& layout,ClockRole role) const
{
    require(bool(data_)&&data_->role==role,"steady VOF velocity cannot change clock role while rebinding");
    // An explicit new layout may name the SAME unchanged native mesh. Both
    // live identities must validate; this is never a changed-tree field remap.
    data_->layout.validateMesh(mesh);layout.validateMesh(mesh);
    const auto& old=data_->layout;
    require(old.domain().axisymmetric()==layout.domain().axisymmetric()&&
        old.domain().nativeFaceCount()==layout.domain().nativeFaceCount()&&
        old.domain().cellRowByIndex()==layout.domain().cellRowByIndex()&&
        old.domain().cellIndexByRow()==layout.domain().cellIndexByRow()&&old.velocityCount()==layout.velocityCount(),
        "steady VOF velocity rebind has a different native selection");
    for(int velocity=0;velocity<old.velocityCount();++velocity)
        require(old.nativeFace(velocity)==layout.nativeFace(velocity)&&old.normalAxis(velocity)==layout.normalAxis(velocity),
            "steady VOF velocity rebind has a different native velocity mapping");
    if(matches(layout,role))return *this;
    return prepare(mesh,layout,data_->values,data_->covered,role);
}

using VelocityHistory=HydroTwoPhaseVelocityHistory;
struct HydroTwoPhaseState::Data
{
    Data(HydroPhaseGeometry phase,HydroNativeValues values,HydroMomentumGeometry geometry,HydroAdvectorReconstruction reconstruction,
        HydroAdvectorBoundary boundary,std::optional<HydroCharacteristicGeometry> characteristics,std::optional<VelocityHistory> history={},
        std::optional<HydroMomentumMass> endMass={},std::optional<HydroTwoPhasePrescribedFlow> flow={},
        std::optional<HydroVofPrescribedFlow> vofFlow={})
        :phase(std::move(phase)),values(std::move(values)),geometry(std::move(geometry)),reconstruction(std::move(reconstruction)),boundary(std::move(boundary)),characteristics(std::move(characteristics)),history(std::move(history)),endMass(std::move(endMass)),flow(std::move(flow)),vofFlow(std::move(vofFlow)){}
    HydroPhaseGeometry phase;HydroNativeValues values;HydroMomentumGeometry geometry;HydroAdvectorReconstruction reconstruction;
    HydroAdvectorBoundary boundary;std::optional<HydroCharacteristicGeometry> characteristics;
    std::optional<VelocityHistory> history; // One actual accepted predecessor, never a state-handle chain.
    std::optional<HydroMomentumMass> endMass; // This state's actual accepted endpoint mass, with its source and clock.
    std::optional<HydroTwoPhasePrescribedFlow> flow;
    std::optional<HydroVofPrescribedFlow> vofFlow; // Distinct actual VOF physical ports, no FirstJet.
    // Only restored states retain one genuine, history-free predecessor.
    // Normal advancement never retains a previous State handle or a chain.
    std::shared_ptr<const HydroTwoPhaseState> restartPredecessor;
};
HydroTwoPhaseState::HydroTwoPhaseState(std::shared_ptr<const Data> data):data_(std::move(data)){}
HydroTwoPhaseState HydroTwoPhaseState::create(const QuadtreeMesh& mesh,HydroPhaseGeometry phase,const HydroNativeValues& values)
{
    supportedGeometry(mesh,phase);const auto walls=closedFaces(mesh,phase.layout());acceptedValues(mesh,phase,values,walls);
    const auto geometry=HydroMomentumGeometry::build(mesh,phase.metrics());
    const auto boundary=HydroAdvectorBoundary::impermeable(mesh,geometry,walls);
    const auto reconstruction=HydroAdvectorReconstruction::prepare(mesh,geometry,boundary);const auto characteristics=HydroCharacteristicGeometry::prepare(mesh,reconstruction);
    return HydroTwoPhaseState(std::make_shared<Data>(std::move(phase),values,geometry,reconstruction,boundary,characteristics));
}
HydroTwoPhaseState HydroTwoPhaseState::createVof(const QuadtreeMesh& mesh,HydroPhaseGeometry phase,const HydroNativeValues& values)
{
    supportedVofGeometry(mesh,phase);const auto walls=closedFaces(mesh,phase.layout());acceptedValues(mesh,phase,values,walls);
    const auto geometry=HydroMomentumGeometry::build(mesh,phase.metrics());
    const auto boundary=HydroAdvectorBoundary::impermeable(mesh,geometry,walls);
    const auto reconstruction=HydroAdvectorReconstruction::prepare(mesh,geometry,boundary);
    HydroTwoPhaseState result(std::make_shared<Data>(std::move(phase),values,geometry,reconstruction,boundary,std::nullopt));
    result.validateFor(mesh);return result;
}
HydroTwoPhaseState HydroTwoPhaseState::createVofRemapped(const QuadtreeMesh& mesh,HydroPhaseGeometry phase,
    const HydroNativeValues& values,HydroMomentumMass retained)
{
    supportedVofGeometry(mesh,phase);
    require(phase.timeLayer().position==LAMP::TimeLayer::Position::End,"VOF remap requires the same accepted End");
    retained.validateVofFor(mesh,phase.vofIdentity(),phase.timeLayer());
    const auto geometry=retained.geometry();
    require(retained.layout().sameIdentity(phase.layout())
        && &geometry.metrics().divergence().matrix()==&phase.metrics().divergence().matrix(),
        "different remapped VOF phase/mass layout or divergence");
    const auto walls=closedFaces(mesh,phase.layout());
    const auto boundary=HydroAdvectorBoundary::impermeable(mesh,geometry,walls);
    const auto reconstruction=HydroAdvectorReconstruction::prepare(mesh,geometry,boundary);
    // No previous-layout field is retained as a velocity predictor.
    HydroTwoPhaseState result(std::make_shared<Data>(std::move(phase),values,geometry,reconstruction,boundary,
        std::nullopt,std::nullopt,std::move(retained)));
    result.validateFor(mesh);return result;
}
HydroTwoPhaseState HydroTwoPhaseState::createVofAccepted(const QuadtreeMesh& mesh,HydroPhaseGeometry phase,
    const HydroNativeValues& values,const HydroTwoPhaseState& previous,std::optional<HydroMomentumMass> retained)
{
    supportedVofGeometry(mesh,phase);previous.validateFor(mesh);
    require(previous.model()==HydroPhaseGeometry::Model::GeometricVof&&!previous.prescribedFlow(),"VOF history requires a genuine closed VOF predecessor");
    require(phase.layout().sameIdentity(previous.values().layout()),"different accepted VOF history layout; explicit transfer/reset required");
    require(&phase.metrics().divergence().matrix()==&previous.phaseGeometry().metrics().divergence().matrix(),"different accepted VOF history divergence");
    const auto clock=phase.timeLayer(),old=previous.timeLayer();
    require(clock.position==Interface::GeometricVofTimeLayer::Position::End&&old.step<std::numeric_limits<std::int64_t>::max()
        &&clock.step==old.step+1&&clock.begin==old.time(),"noncontiguous VOF accepted history clock");
    if(previous.vofPrescribedFlow()){
        const auto walls=vofPrescribedWalls(mesh,phase,values,*previous.vofPrescribedFlow());
        require(walls==previous.boundary().velocityDofs(),"different accepted VOF physical wall partition");
        acceptedValues(mesh,phase,values,walls,!previous.vofPrescribedFlow()->pressureOutlets.empty());
    }else acceptedValues(mesh,phase,values,previous.boundary().velocityDofs());
    VelocityHistory history{previous.values(),old,{},HydroPhaseGeometry::Model::GeometricVof,previous.phaseGeometry().vofIdentity()};
    HydroTwoPhaseState result(std::make_shared<Data>(std::move(phase),values,previous.momentumGeometry(),previous.reconstruction(),previous.boundary(),
        std::nullopt,std::move(history),std::move(retained),std::nullopt,previous.data_->vofFlow));
    result.validateFor(mesh);return result;
}
HydroTwoPhaseState HydroTwoPhaseState::createPrescribedFlow(const QuadtreeMesh& mesh,HydroPhaseGeometry phase,
    const HydroNativeValues& values,HydroTwoPhasePrescribedFlow flow)
{
    supportedGeometry(mesh,phase);const auto walls=prescribedWalls(mesh,phase,values,flow);acceptedValues(mesh,phase,values,walls,!flow.pressureOutlets.empty());
    const auto geometry=HydroMomentumGeometry::build(mesh,phase.metrics());
    const auto boundary=HydroAdvectorBoundary::impermeable(mesh,geometry,walls);
    const auto reconstruction=HydroAdvectorReconstruction::prepare(mesh,geometry,boundary);const auto characteristics=HydroCharacteristicGeometry::prepare(mesh,reconstruction);
    return HydroTwoPhaseState(std::make_shared<Data>(std::move(phase),values,geometry,reconstruction,boundary,characteristics,
        std::nullopt,std::nullopt,std::move(flow)));
}
HydroTwoPhaseState HydroTwoPhaseState::createVofPrescribedFlow(const QuadtreeMesh& mesh,HydroPhaseGeometry phase,
    const HydroNativeValues& values,HydroVofPrescribedFlow flow)
{
    supportedVofGeometry(mesh,phase);const auto walls=vofPrescribedWalls(mesh,phase,values,flow);
    acceptedValues(mesh,phase,values,walls,!flow.pressureOutlets.empty());
    const auto geometry=HydroMomentumGeometry::build(mesh,phase.metrics());
    const auto boundary=HydroAdvectorBoundary::impermeable(mesh,geometry,walls);
    const auto reconstruction=HydroAdvectorReconstruction::prepare(mesh,geometry,boundary);
    HydroTwoPhaseState result(std::make_shared<Data>(std::move(phase),values,geometry,reconstruction,boundary,
        std::nullopt,std::nullopt,std::nullopt,std::nullopt,std::move(flow)));
    result.validateFor(mesh);return result;
}
HydroTwoPhaseState HydroTwoPhaseState::createVofPrescribedFlowRemapped(const QuadtreeMesh& mesh,HydroPhaseGeometry phase,
    const HydroNativeValues& values,HydroMomentumMass retained,HydroVofPrescribedFlow flow)
{
    supportedVofGeometry(mesh,phase);
    require(phase.timeLayer().position==LAMP::TimeLayer::Position::End,"VOF open remap requires the same accepted End");
    retained.validateVofFor(mesh,phase.vofIdentity(),phase.timeLayer());const auto geometry=retained.geometry();
    require(retained.layout().sameIdentity(phase.layout())
        && &geometry.metrics().divergence().matrix()==&phase.metrics().divergence().matrix(),
        "different remapped open VOF phase/mass layout or divergence");
    const auto walls=vofPrescribedWalls(mesh,phase,values,flow);
    acceptedValues(mesh,phase,values,walls,!flow.pressureOutlets.empty());
    const auto boundary=HydroAdvectorBoundary::impermeable(mesh,geometry,walls);
    const auto reconstruction=HydroAdvectorReconstruction::prepare(mesh,geometry,boundary);
    HydroTwoPhaseState result(std::make_shared<Data>(std::move(phase),values,geometry,reconstruction,boundary,
        std::nullopt,std::nullopt,std::move(retained),std::nullopt,std::move(flow)));
    result.validateFor(mesh);return result;
}
const HydroVofPrescribedFlow* HydroTwoPhaseState::vofPrescribedFlow() const noexcept
{return data_->vofFlow?&*data_->vofFlow:nullptr;}
const Interface::GeometricVofAdvectorSlab::Boundary* HydroTwoPhaseState::vofTransportBoundary() const noexcept
{
    // This State does not retain a publication-bound transport Boundary.
    // Callers prepare the genuine current-owner boundary from vofFlow rules.
    return nullptr;
}
const HydroTwoPhasePrescribedFlow* HydroTwoPhaseState::prescribedFlow() const noexcept
{return data_->flow?&*data_->flow:nullptr;}
const HydroPhaseGeometry& HydroTwoPhaseState::phaseGeometry() const noexcept{return data_->phase;}
HydroPhaseGeometry::Model HydroTwoPhaseState::model() const noexcept{return phaseGeometry().model();}
const Interface::HermiteGalsPublication& HydroTwoPhaseState::publication() const
{require(model()==HydroPhaseGeometry::Model::FirstJet,"FirstJet publication port requires a FirstJet state");return data_->phase.publication();}
const HydroTwoPhaseVelocityHistory* HydroTwoPhaseState::velocityHistory() const noexcept{return data_->history?&*data_->history:nullptr;}
const HydroMomentumMass* HydroTwoPhaseState::retainedMass() const noexcept{return data_->endMass?&*data_->endMass:nullptr;}
const HydroNativeValues& HydroTwoPhaseState::values() const noexcept{return data_->values;}
std::vector<Vector2D> HydroTwoPhaseState::sampleCellCenterVelocities(const QuadtreeMesh& mesh) const
{
    validateFor(mesh);
    if(model()==HydroPhaseGeometry::Model::GeometricVof
        &&phaseGeometry().vofPublication().geometryKind()==Interface::GeometricVofGeometryKind::LocalPlic)
    {
        // Same native P0 traces as LocalPlic transport. Sum actual hanging
        // subfaces into each whole side, then evaluate its RT0/Piola field.
        // This diagnostic never replaces or modifies the native face DOFs.
        const auto& layout=values().layout();const auto& domain=layout.domain();
        std::vector<std::array<long double,4>> flux(layout.pressureCount()),area(layout.pressureCount());
        for(int u=0;u<layout.velocityCount();++u)
        {
            const auto f=domain.face(mesh,layout.nativeFace(u));
            const int axis=layout.normalAxis(u)==HydroDofLayout::NormalAxis::First?0:1;
            const long double a=phaseGeometry().metrics().faceMeasures()[u],q=a*values().normalVelocity()[u];
            for(int native:{f.negativeCell,f.positiveCell})if(native>=0)
            {
                const int p=layout.pressureDof(native);if(p<0)continue;
                const int side=2*axis+(native==f.negativeCell?1:0);
                flux[p][side]+=q;area[p][side]+=a;
            }
        }
        std::vector<Vector2D> result(layout.pressureCount());
        for(int p=0;p<layout.pressureCount();++p)
        {
            const auto& a=area[p];const auto& q=flux[p];
            require(a[1]>0&&a[2]>0&&a[3]>0&&(a[0]>0||domain.axisymmetric()),"incomplete VOF native side support");
            long double first=0;
            if(domain.axisymmetric())
            {
                const auto& box=phaseGeometry().vofPublication().plicCells()[layout.nativeCell(p)].rectangle;
                const long double lo=box.lower.x,hi=box.upper.x,r=(lo+hi)/2;
                require(r>0&&hi>lo,"invalid VOF radial cell support");
                const auto xi=(r*r-lo*lo)/(hi*hi-lo*lo);
                first=((1-xi)*q[0]+xi*q[1])/((a[0]+a[1])/2);
            }
            else first=(q[0]/a[0]+q[1]/a[1])/2;
            const auto second=(q[2]/a[2]+q[3]/a[3])/2;
            require(std::isfinite(first)&&std::isfinite(second),"nonfinite VOF RT0 cell sample");
            result[p]={static_cast<double>(first),static_cast<double>(second)};
        }
        return result;
    }
    const auto& regular=reconstruction().regularGeometry();
    const auto coefficients=regular.apply(mesh,values(),&boundary());
    const auto& layout=values().layout();std::vector<Vector2D> result(layout.pressureCount());
    for(int p=0;p<layout.pressureCount();++p)
    {
        const auto native=layout.nativeCell(p);
        const auto point=mesh.cells(mesh.cellIndexToLabel[native]).centerCoord();
        result[p]=regular.cellJet(coefficients,p,point).jet.value;
    }
    return result;
}
const HydroMomentumGeometry& HydroTwoPhaseState::momentumGeometry() const noexcept{return data_->geometry;}
const HydroAdvectorReconstruction& HydroTwoPhaseState::reconstruction() const noexcept{return data_->reconstruction;}
const HydroAdvectorBoundary& HydroTwoPhaseState::boundary() const noexcept{return data_->boundary;}
const HydroCharacteristicGeometry& HydroTwoPhaseState::characteristicGeometry() const
{require(model()==HydroPhaseGeometry::Model::FirstJet&&data_->characteristics.has_value(),"FirstJet characteristic port requires a FirstJet state and characteristic geometry");return *data_->characteristics;}
Interface::HermiteGalsIdentity HydroTwoPhaseState::geometryIdentity() const{return data_->phase.geometryIdentity();}
Interface::HermiteGalsTimeLayer HydroTwoPhaseState::timeLayer() const{return data_->phase.timeLayer();}
bool HydroTwoPhaseState::sameIdentity(const HydroTwoPhaseState& other) const noexcept{return data_==other.data_;}
void HydroTwoPhaseState::validateVofPhaseInBatch(const QuadtreeMesh& mesh,const HydroPhaseGeometry& phase,
    Interface::GeometricVofIdentity owner,LAMP::TimeLayer time,Detail::HydroStateMeshValidationBatch& batch)
{
    const auto& source=phase.vofPublication();
    const auto& domain=source.domain(); // Original checked()/consumed check FIRST.
    time.time();const auto actual=source.timeLayer();
    // Only a repeated pure mesh scan is omitted. Failed/foreign intrinsic
    // predicates use the original public entry and its original first error.
    if(batch.sameBinding(mesh,domain)&&owner==source.identity()&&owner.value()!=0
        &&time.step==actual.step&&time.begin==actual.begin&&time.end==actual.end&&time.position==actual.position)return;
    phase.validateVofFor(mesh,owner,time);
    batch.rememberAfterValidation(mesh,domain);
}
void HydroTwoPhaseState::supportedVofGeometryInBatch(const QuadtreeMesh& mesh,const HydroPhaseGeometry& phase,
    Detail::HydroStateMeshValidationBatch& batch)
{
    supportedVofGeometryImpl(mesh,phase,[&]{batch.validateLayout(mesh,phase.layout());},
        [&]{validateVofPhaseInBatch(mesh,phase,phase.vofIdentity(),phase.timeLayer(),batch);});
}
void HydroTwoPhaseState::acceptedValuesInBatch(const QuadtreeMesh& mesh,const HydroPhaseGeometry& phase,
    const HydroNativeValues& values,const std::vector<int>& walls,bool pressureAnchored,
    Detail::HydroStateMeshValidationBatch& batch)
{
    acceptedValuesImpl(mesh,phase,values,walls,pressureAnchored,
        [&]{batch.validateValues(mesh,values,phase.layout());},
        [&]{return phase.metrics().divergence().applyInBatch(mesh,values,batch);});
}
void HydroTwoPhaseState::validateFor(const QuadtreeMesh& mesh) const
{
    Detail::HydroStateMeshValidationBatch batch;
    validateForInBatch(mesh,batch);
}
void HydroTwoPhaseState::validateForInBatch(const QuadtreeMesh& mesh,Detail::HydroStateMeshValidationBatch& batch) const
{
    const bool vof=model()==HydroPhaseGeometry::Model::GeometricVof;
    if(vof){supportedVofGeometryInBatch(mesh,phaseGeometry(),batch);require(!data_->characteristics&&!data_->flow,"closed VOF state has no FirstJet characteristic or flow owner");}
    else{require(!data_->vofFlow,"FirstJet state has a VOF physical flow owner");supportedGeometry(mesh,phaseGeometry());characteristicGeometry().validateFor(mesh);}
    if(vof)batch.validateLayout(mesh,momentumGeometry().layout());else momentumGeometry().validateMesh(mesh);
    boundary().validateFor(mesh,momentumGeometry());
    require(reconstruction().momentumGeometry().sameIdentity(momentumGeometry()),"different retained reconstruction geometry");
    const bool anchored=data_->flow&&!data_->flow->pressureOutlets.empty();
    if(vof&&data_->vofFlow){
        // The ENTIRE helper can call endpointVelocity. No identity proof or
        // presumed-safe domain/getter reference crosses that user callback.
        batch.clear();
        require(vofPrescribedWalls(mesh,phaseGeometry(),values(),*data_->vofFlow)==boundary().velocityDofs(),"different VOF physical wall owner");
        acceptedValuesInBatch(mesh,phaseGeometry(),values(),boundary().velocityDofs(),!data_->vofFlow->pressureOutlets.empty(),batch);
    }else{
        if(vof)acceptedValuesInBatch(mesh,phaseGeometry(),values(),boundary().velocityDofs(),anchored,batch);
        else acceptedValues(mesh,phaseGeometry(),values(),boundary().velocityDofs(),anchored);
        if(vof)require(boundary().completeBoundaryLines()&&closedFaces(mesh,phaseGeometry().layout())==boundary().velocityDofs(),"closed VOF state requires the complete closed wall set");
    }
    if(data_->flow)require(prescribedWalls(mesh,phaseGeometry(),values(),*data_->flow)==boundary().velocityDofs(),"different prescribed-flow wall owner");
    if(data_->history){const auto& previous=*data_->history;if(vof)batch.validateValues(mesh,previous.values,values().layout());else previous.values.validateFor(values().layout(),mesh);
        require(previous.model==model(),"different accepted velocity history model");
        if(vof)require(previous.source.value()==0&&previous.vofSource.value()!=0,"VOF history has a fabricated or missing source identity");
        else require(previous.source.value()!=0,"FirstJet history has a fabricated or missing source identity");
        const auto clock=timeLayer();require(clock.position==Interface::HermiteGalsTimeLayer::Position::End
            &&clock.step>0&&previous.clock.step==clock.step-1&&previous.clock.time()==clock.begin
            &&(previous.clock.position==Interface::HermiteGalsTimeLayer::Position::Initial||previous.clock.position==Interface::HermiteGalsTimeLayer::Position::End),
            "noncontiguous accepted velocity history");}
    if(data_->endMass){const auto& mass=*data_->endMass;
        require(mass.model()==model()&&mass.layout().sameIdentity(values().layout())&&mass.geometry().sameIdentity(momentumGeometry()),"different accepted retained mass model/layout/geometry");
        require(&mass.geometry().metrics().divergence().matrix()==&phaseGeometry().metrics().divergence().matrix(),"different accepted retained mass divergence");
        if(vof)validateVofPhaseInBatch(mesh,mass.phaseGeometry(),phaseGeometry().vofIdentity(),timeLayer(),batch);
        else mass.validateFor(mesh,geometryIdentity(),timeLayer());
    }
}
void HydroTwoPhaseState::validateVofFor(const QuadtreeMesh& mesh,Interface::GeometricVofIdentity identity,Interface::GeometricVofTimeLayer clock) const
{
    require(model()==HydroPhaseGeometry::Model::GeometricVof,"VOF state authentication requires its genuine model");
    phaseGeometry().validateVofFor(mesh,identity,clock); // Original complete live entry.
    Detail::HydroStateMeshValidationBatch batch;
    batch.rememberAfterValidation(mesh,phaseGeometry().vofPublication().domain());
    validateForInBatch(mesh,batch);
}
namespace
{
namespace Checkpoint=LAMP::IO::Checkpoint;
bool checkpointBitsEqual(const std::vector<double>& a,const std::vector<double>& b)
{return a.size()==b.size()&&(a.empty()||std::memcmp(a.data(),b.data(),a.size()*sizeof(double))==0);}
std::string checkpointFlow(const HydroVofPrescribedFlow* flow)
{
    Checkpoint::Writer out("LAMP.HydroVofFlow.1");out.boolean(flow!=nullptr);
    if(flow){
        std::vector<int> rules;rules.reserve(flow->rulesByNativeRow.size());
        for(const auto rule:flow->rulesByNativeRow)rules.push_back(static_cast<int>(rule));
        out.values(rules);out.values(flow->prescribedVelocityDofs);out.values(flow->normalVelocity);
        out.count(flow->pressureOutlets.size());
        for(const auto& port:flow->pressureOutlets){out.scalar(port.velocityDof);out.scalar(port.pressure);}
        out.boolean(flow->pressureOutletBackflowModel.has_value());
        if(flow->pressureOutletBackflowModel)out.scalar(static_cast<int>(*flow->pressureOutletBackflowModel));
        out.boolean(flow->pressureDrivenReservoirs);out.values(flow->stationaryNoSlipVelocityDofs);
        // Functions are supplied by the unchanged App configuration, never
        // deserialized code. Their declared capability and physical port data
        // must still match; the App authenticates the complete configuration.
        out.boolean(bool(flow->endpointVelocity));out.boolean(bool(flow->pressureOutletBackflowMomentum));
    }
    return out.finish();
}
std::vector<int> checkpointVelocityFaces(const HydroDofLayout& layout)
{
    std::vector<int> faces;faces.reserve(layout.velocityCount());
    for(int u=0;u<layout.velocityCount();++u)faces.push_back(layout.nativeFace(u));return faces;
}
std::vector<int> checkpointVelocityAxes(const HydroDofLayout& layout)
{
    std::vector<int> axes;axes.reserve(layout.velocityCount());
    for(int u=0;u<layout.velocityCount();++u)axes.push_back(static_cast<int>(layout.normalAxis(u)));return axes;
}
std::string checkpointMass(const HydroMomentumMass& mass)
{
    require(mass.model()==HydroPhaseGeometry::Model::GeometricVof
        &&mass.representation()==HydroMomentumMass::DensityRepresentation::CellAverageProjection,
        "VOF checkpoint supports the production cell-average momentum mass");
    Checkpoint::Writer out("LAMP.HydroVofMass.1");out.scalar(static_cast<int>(mass.representation()));
    out.scalar(mass.liquidDensity());out.scalar(mass.gasDensity());out.values(mass.masses());
    const auto measure=[&](const HydroPhaseMeasure& q){out.scalar(q.liquid);out.scalar(q.gas);};
    out.count(mass.representedPhaseVolumes().size());for(const auto& q:mass.representedPhaseVolumes())measure(q);
    out.count(mass.representedOverlapPhaseVolumes().size());
    for(const auto& pieces:mass.representedOverlapPhaseVolumes()){
        out.count(pieces.size());for(const auto& q:pieces)measure(q);
    }
    out.scalar(mass.layout().pressureCount());
    for(int p=0;p<mass.layout().pressureCount();++p)for(int component=0;component<2;++component)
        for(int phase=0;phase<2;++phase)out.boolean(mass.usesFiniteP0Partition(p,component,phase));
    return out.finish();
}
}
std::string HydroTwoPhaseState::checkpoint(const QuadtreeMesh& mesh,const HydroTwoPhaseState* previous) const
{
    require(model()==HydroPhaseGeometry::Model::GeometricVof&&!prescribedFlow(),
        "checkpoint requires a genuine VOF state");validateFor(mesh);
    const auto& publication=phaseGeometry().vofPublication();const auto& layout=values().layout();
    require(publication.domain().cellIndexByRow().size()==static_cast<std::size_t>(mesh.cellIndexNum()),
        "VOF checkpoint currently requires the complete native rectangle");
    const auto* history=velocityHistory();
    if(history){
        if(!previous)previous=data_->restartPredecessor.get();
        require(previous!=nullptr,"VOF checkpoint with velocity history requires its genuine predecessor");
        previous->validateVofFor(mesh,history->vofSource,history->clock);
        require(history->model==HydroPhaseGeometry::Model::GeometricVof
            &&history->values.layout().sameIdentity(previous->values().layout())
            &&checkpointBitsEqual(history->values.normalVelocity(),previous->values().normalVelocity())
            &&checkpointBitsEqual(history->values.pressure(),previous->values().pressure()),
            "VOF checkpoint predecessor differs from the actual velocity history");
        require(checkpointFlow(vofPrescribedFlow())==checkpointFlow(previous->vofPrescribedFlow()),
            "VOF checkpoint predecessor has different physical ports");
    }
    Checkpoint::Writer out("LAMP.HydroVofState.1");const NativeMeshIdentity identity(mesh);
    out.values(identity.checkpointWords());out.boolean(publication.domain().axisymmetric());
    out.values(publication.domain().cellIndexByRow());out.values(checkpointVelocityFaces(layout));out.values(checkpointVelocityAxes(layout));
    out.text(checkpointFlow(vofPrescribedFlow()));out.text(publication.checkpoint());
    out.values(values().normalVelocity());out.values(values().pressure());out.boolean(history!=nullptr);
    if(history){out.text(previous->phaseGeometry().vofPublication().checkpoint());
        out.values(previous->values().normalVelocity());out.values(previous->values().pressure());}
    out.boolean(retainedMass()!=nullptr);if(retainedMass())out.text(checkpointMass(*retainedMass()));
    return out.finish();
}
HydroTwoPhaseState HydroTwoPhaseState::restoreCheckpoint(const QuadtreeMesh& mesh,const std::string& bytes,
    std::optional<HydroVofPrescribedFlow> configuredFlow)
{
    Checkpoint::Reader in(bytes,"LAMP.HydroVofState.1");const NativeMeshIdentity actualMesh(mesh);
    require(in.values<std::uint64_t>()==actualMesh.checkpointWords(),"VOF checkpoint native mesh geometry/connectivity differs");
    const bool radial=in.boolean();const auto nativeCells=in.values<int>(),nativeFaces=in.values<int>(),axes=in.values<int>();
    require(in.text()==checkpointFlow(configuredFlow?&*configuredFlow:nullptr),
        "VOF checkpoint physical port configuration differs");
    const auto publicationBytes=in.text();const auto velocity=in.values<double>(),pressure=in.values<double>();
    const bool hasHistory=in.boolean();std::string predecessorBytes;std::vector<double> oldVelocity,oldPressure;
    if(hasHistory){predecessorBytes=in.text();oldVelocity=in.values<double>();oldPressure=in.values<double>();}
    const bool hasMass=in.boolean();const auto massBytes=hasMass?in.text():std::string{};in.finish();
    FaceSegmentGraph graph;graph.build(mesh);const auto domain=NativeLeafDomain::build(mesh,graph,radial);
    const auto publication=Interface::GeometricVofPublication::restoreCheckpoint(mesh,domain,publicationBytes);
    const auto layout=HydroDofLayout::build(mesh,publication.domain());
    require(nativeCells==publication.domain().cellIndexByRow()&&nativeFaces==checkpointVelocityFaces(layout)
        &&axes==checkpointVelocityAxes(layout),"VOF checkpoint native field ordering differs");
    const auto metrics=HydroPhysicalDivergence::build(mesh,HydroIntegralDivergence::build(mesh,layout));
    const auto phase=HydroPhaseGeometry::build(mesh,metrics,publication);const HydroNativeValues restoredValues(layout,velocity,pressure);
    std::optional<HydroTwoPhaseState> predecessor;
    if(hasHistory){
        const auto oldPublication=Interface::GeometricVofPublication::restoreCheckpoint(mesh,domain,predecessorBytes);
        const auto oldPhase=HydroPhaseGeometry::build(mesh,metrics,oldPublication);
        const HydroNativeValues oldValues(layout,oldVelocity,oldPressure);
        predecessor=configuredFlow?createVofPrescribedFlow(mesh,oldPhase,oldValues,*configuredFlow):createVof(mesh,oldPhase,oldValues);
    }
    auto restored=predecessor?createVofAccepted(mesh,phase,restoredValues,*predecessor):
        configuredFlow?createVofPrescribedFlow(mesh,phase,restoredValues,*configuredFlow):createVof(mesh,phase,restoredValues);
    auto data=std::make_shared<Data>(*restored.data_);
    if(hasMass){
        Checkpoint::Reader material(massBytes,"LAMP.HydroVofMass.1");
        const auto representation=material.scalar<int>();
        require(representation==static_cast<int>(HydroMomentumMass::DensityRepresentation::CellAverageProjection),
            "VOF checkpoint has an unsupported momentum density representation");
        const double liquidDensity=material.scalar<double>(),gasDensity=material.scalar<double>();
        auto mass=HydroMomentumMass::buildCellAverageProjection(mesh,restored.momentumGeometry(),phase,liquidDensity,gasDensity);
        // The production P0 mass is a deterministic function of BOTH stocks
        // and the same native overlap geometry. Rebuilding is allowed only
        // when every saved M, represented phase moment and finite flag agrees
        // exactly; no averaging or tolerance can replace its canonical values.
        require(checkpointMass(mass)==massBytes,"restored momentum mass differs from the exact checkpoint moment payload");
        data->endMass=std::move(mass);
    }
    if(predecessor)data->restartPredecessor=std::make_shared<HydroTwoPhaseState>(*predecessor);
    restored.data_=std::move(data);restored.validateFor(mesh);return restored;
}
struct HydroTwoPhaseStep::Data
{
    Data(HydroTwoPhaseState begin,HydroTwoPhaseState candidate,Trial trial,HydroNativeValues midpoint,HydroKepConvectionResult reassessed,HydroTwoPhaseStepReport report)
        :begin(std::move(begin)),candidate(std::move(candidate)),trial(std::move(trial)),midpoint(std::move(midpoint)),reassessed(std::move(reassessed)),report(std::move(report)){}
    HydroTwoPhaseState begin,candidate;Trial trial;HydroNativeValues midpoint;HydroKepConvectionResult reassessed;HydroTwoPhaseStepReport report;
};
HydroTwoPhaseStep::HydroTwoPhaseStep(std::shared_ptr<const Data> data):data_(std::move(data)){}
HydroTwoPhaseStep HydroTwoPhaseStep::prepare(const QuadtreeMesh& mesh,const HydroTwoPhaseState& accepted,
    Interface::HermiteGalsTimeLayer clock,const HydroTwoPhaseStepConfig& config)
{
    require(accepted.model()==HydroPhaseGeometry::Model::FirstJet,"FirstJet spacetime preparation requires a FirstJet state");
    accepted.validateFor(mesh);using P=Interface::HermiteGalsTimeLayer::Position;const auto oldClock=accepted.timeLayer();const double dt=clock.end-clock.begin;
    require(oldClock.step<std::numeric_limits<std::int64_t>::max()&&clock.step==oldClock.step+1&&clock.position==P::End&&clock.begin==oldClock.time()&&std::isfinite(dt)&&dt>0,"invalid next End clock");clock.time();
    require(std::isfinite(config.liquidDensity)&&config.liquidDensity>0&&std::isfinite(config.gasDensity)&&config.gasDensity>0,"invalid two-phase densities");
    using Density=HydroMomentumMass::DensityRepresentation;
    require(config.momentumDensityRepresentation==Density::SharpGeometry
        ||config.momentumDensityRepresentation==Density::CellAverageProjection,"invalid momentum density representation");
    require(std::isfinite(config.liquidViscosity)&&config.liquidViscosity>=0&&std::isfinite(config.gasViscosity)&&config.gasViscosity>=0,"invalid phase viscosities");
    require(std::isfinite(config.bodyAcceleration.x)&&std::isfinite(config.bodyAcceleration.y),"invalid body acceleration");
    const bool bodyActive=config.bodyAcceleration.x!=0||config.bodyAcceleration.y!=0;
    require(std::isfinite(config.surfaceTension)&&config.surfaceTension>=0,"invalid surface tension");require(config.maximumIterations>=1&&config.maximumIterations<=12,"corrector iteration budget must be 1..12");
    require(config.capillaryModel==HydroCapillaryModel::BalancedCsf||config.capillaryModel==HydroCapillaryModel::SharpNormalStress,
        "invalid capillary model");
    require(config.csfCurvature.forceMethod==HydroCsfForceMethod::MotionAdjoint
        ||config.csfCurvature.forceMethod==HydroCsfForceMethod::HeightFunction
        ||config.csfCurvature.forceMethod==HydroCsfForceMethod::HeightFunctionIncrement,"invalid CSF force method");
    if(config.capillaryModel==HydroCapillaryModel::BalancedCsf&&config.csfCurvature.forceMethod==HydroCsfForceMethod::MotionAdjoint)
        require(std::isfinite(config.csfCurvature.relativeQuadratureTolerance)&&config.csfCurvature.relativeQuadratureTolerance>0&&config.csfCurvature.relativeQuadratureTolerance<1,
            "invalid CSF curvature integration tolerance");
    using Gradient=Interface::HermiteGalsGradientUpdate;
    require(config.publication.gradientUpdate==Gradient::CharacteristicPullback
        ||config.publication.gradientUpdate==Gradient::CompatibleIncrement,"invalid gradient update");
    if(config.publication.gradientUpdate==Gradient::CompatibleIncrement)
        require(config.capillaryModel==HydroCapillaryModel::BalancedCsf
            &&(config.csfCurvature.forceMethod==HydroCsfForceMethod::HeightFunction
                ||config.csfCurvature.forceMethod==HydroCsfForceMethod::HeightFunctionIncrement),
            "compatible gradient increment requires height-function or height-function-increment force");
    if(config.csfCurvature.forceMethod==HydroCsfForceMethod::HeightFunctionIncrement)
        require(config.publication.gradientUpdate==Gradient::CompatibleIncrement,
            "height-function-increment requires compatible gradient increment");
    if(config.csfCurvature.forceMethod==HydroCsfForceMethod::HeightFunction
        ||config.csfCurvature.forceMethod==HydroCsfForceMethod::HeightFunctionIncrement){
        require(config.capillaryModel==HydroCapillaryModel::BalancedCsf,"height-function requires balancedCsf force model");
        require(!config.interfaceLoad,"height-function external interface traction is unavailable");
        require(!accepted.prescribedFlow(),"height-function requires homogeneous closed boundary");
        require(!accepted.phaseGeometry().layout().domain().axisymmetric(),"height-function requires XY geometry");
        require(accepted.phaseGeometry().layout().domain().uniformLeaves()
            &&accepted.reconstruction().regularGeometry().usesCompatibleTensor(),"height-function requires complete uniform XY motion support");
        const UniformLeafGrid heightGrid(mesh,CartesianGridScope::UniformLeaves);
        require(heightGrid.nx()>=3&&heightGrid.ny()>=3,"height-function requires at least three uniform XY cells in each direction");
        require(accepted.phaseGeometry().layout().pressureCount()==mesh.cellIndexNum()
            &&accepted.publication().state().primary().geometryCells.size()==static_cast<std::size_t>(mesh.cellIndexNum()),
            "height-function requires complete native geometry and pressure domain");
        require(config.publication.scalarRepresentation==Interface::HermiteGalsScalarRepresentation::SingleQuartic
            &&accepted.publication().state().geometry().snapshot().reconstruction()==LevelSetHermiteGeometry::Reconstruction::EnrichedQuartic,
            "height-function requires SingleQuartic scalar representation");
    }
    require(std::isfinite(config.velocityUpdateTolerance)&&config.velocityUpdateTolerance>0&&config.velocityUpdateTolerance<1
        &&std::isfinite(config.geometryUpdateTolerance)&&config.geometryUpdateTolerance>0&&config.geometryUpdateTolerance<1
        &&std::isfinite(config.momentumPowerTolerance)&&config.momentumPowerTolerance>0&&config.momentumPowerTolerance<1,
        "invalid corrector update tolerances");
    using Scheme=HydroTwoPhaseTimeScheme;
    require(config.timeScheme==Scheme::FullyCoupledMidpoint||config.timeScheme==Scheme::ExtrapolatedMidpoint
        ||config.timeScheme==Scheme::AdaptiveLinearMidpoint,"invalid two-phase time scheme");
    const bool linear=config.timeScheme!=Scheme::FullyCoupledMidpoint;
    if(config.interfaceLoad)require(linear&&config.capillaryModel==HydroCapillaryModel::BalancedCsf,
        "external interface load requires the AC balancedCsf path; Sharp/Picard are unsupported");
    if(linear){
        require(!config.dissipation,"linear midpoint does not support nonlinear sensor dissipation");
        require(std::isfinite(config.interfaceMismatchTolerance)&&config.interfaceMismatchTolerance>0
            &&std::isfinite(config.maximumCorrectionRelative)&&config.maximumCorrectionRelative>0
            &&std::isfinite(config.correctionReferenceVelocity)&&config.correctionReferenceVelocity>0,
            "invalid linear midpoint controls");
    }
    if(config.timeScheme==Scheme::AdaptiveLinearMidpoint)
        require(std::isfinite(config.capillaryCorrectionCourant)&&config.capillaryCorrectionCourant>0
            &&std::isfinite(config.maximumCapillaryCourant)&&config.maximumCapillaryCourant>0
            &&config.capillaryCorrectionCourant<=config.maximumCapillaryCourant,"invalid adaptive capillary controls");
    const auto& geometry=accepted.momentumGeometry();const auto& l=geometry.layout();const auto& reconstruction=accepted.reconstruction();const auto& before=accepted.values();
    HydroTwoPhaseStepReport report;
    report.hasExternalInterfaceLoad=static_cast<bool>(config.interfaceLoad);
    report.capillaryPowerAvailable=!report.hasExternalInterfaceLoad;
    if(!report.capillaryPowerAvailable){report.capillaryPower=std::numeric_limits<long double>::quiet_NaN();
        report.capillaryEnergyDefect=std::numeric_limits<long double>::quiet_NaN();}
    report.deltaTime=dt;report.surfaceTension=config.surfaceTension;report.momentumDensityRepresentation=config.momentumDensityRepresentation;report.capillaryModel=config.capillaryModel;report.timeScheme=config.timeScheme;
    if(config.timeScheme==Scheme::AdaptiveLinearMidpoint&&config.surfaceTension>0
        &&!accepted.publication().topology().curves.empty()){
        const auto spacing=l.domain().spacing();const double h=std::min(spacing.x,spacing.y);
        require(std::isfinite(h)&&h>0,"invalid capillary grid scale");
        // C_sigma=dt/sqrt(.5*(rho_l+rho_g)*h^3/(pi*sigma)).
        // Log evaluation avoids intermediate density sums/products overflowing.
        const double rhoScale=std::max(config.liquidDensity,config.gasDensity);
        const double logMeanDensity=std::log(rhoScale)+std::log(.5*(config.liquidDensity/rhoScale+config.gasDensity/rhoScale));
        const double courant=std::exp(std::log(dt)+.5*(std::log(std::acos(-1.))+std::log(config.surfaceTension)
            -logMeanDensity-3*std::log(h)));
        if(!std::isfinite(courant)||courant>config.maximumCapillaryCourant){
            std::ostringstream message;message<<"[HydroTwoPhaseStep] capillary Courant "<<courant
                <<" exceeds configured bound "<<config.maximumCapillaryCourant;
            throw HydroTwoPhaseStepRetry(HydroTwoPhaseRetryReason::Capillary,message.str());
        }
        report.capillaryCourant=courant;
        report.capillaryForcedCorrection=courant>=config.capillaryCorrectionCourant;
    }
    require(!l.domain().axisymmetric()||config.bodyAcceleration.x==0,"RZ body acceleration must be axial");
    const auto* flow=accepted.prescribedFlow();
    const bool pressureActive=flow&&!flow->pressureOutlets.empty();
    std::vector<int> freeOutflow;
    if(pressureActive)for(const auto& outlet:flow->pressureOutlets)freeOutflow.push_back(outlet.velocityDof);
    if(flow){require(!config.dissipation,"prescribed-flow step rejects nonlinear sensor dissipation");
        require(config.spaceTime.firstJetBoundaryTreatment==Interface::HermiteGalsBoundaryTreatment::InflowHistory,
            "prescribed-flow step requires explicit inlet-only FirstJet boundary treatment");}
    const auto previous=reconstruction.evaluate(mesh,before,accepted.geometryIdentity(),oldClock,accepted.boundary());const auto pressure=HydroPressureLoad::build(geometry.metrics().divergence());
    const auto makeMass=[&](const HydroPhaseGeometry& phase){
        return config.momentumDensityRepresentation==Density::SharpGeometry
            ?HydroMomentumMass::build(mesh,geometry,phase,config.liquidDensity,config.gasDensity)
            :HydroMomentumMass::buildCellAverageProjection(mesh,geometry,phase,config.liquidDensity,config.gasDensity);
    };
    // The previous accepted endpoint is this transaction's exact begin source.
    const auto beginMass=[&]{
        const auto& retained=accepted.data_->endMass;
        if(!retained||retained->liquidDensity()!=config.liquidDensity||retained->gasDensity()!=config.gasDensity
            ||retained->representation()!=config.momentumDensityRepresentation)
            return makeMass(accepted.phaseGeometry());
        retained->validateFor(mesh,accepted.geometryIdentity(),oldClock);
        require(retained->geometry().sameIdentity(geometry),"different retained begin momentum geometry");
        require(retained->layout().sameIdentity(l),"different retained begin layout");
        require(&retained->geometry().metrics().divergence().matrix()==&accepted.phaseGeometry().metrics().divergence().matrix(),
            "different retained begin D_I");
        require(&retained->phaseGeometry().publication().state().geometry().snapshot()==&accepted.publication().state().geometry().snapshot(),
            "different retained begin scalar source");
        LAMP::Geometry::Detail::ExactExpansion::validateActualEnvironment();
        return *retained;
    }();
    std::vector<double> guess=before.normalVelocity();std::optional<Trial> preceding;
    report.usedVelocityHistory=linear&&accepted.data_->history.has_value()&&!config.resetVelocityHistory;
    report.startupCorrection=linear&&!report.usedVelocityHistory;
    if(report.usedVelocityHistory){
        const auto& history=*accepted.data_->history;const double oldDt=oldClock.end-oldClock.begin;
        require(std::isfinite(oldDt)&&oldDt>0,"invalid accepted history timestep");
        for(std::size_t i=0;i<guess.size();++i){guess[i]+=dt/oldDt*(before.normalVelocity()[i]-history.values.normalVelocity()[i]);
            require(std::isfinite(guess[i]),"nonfinite extrapolated endpoint velocity");}
    }
    const int stageLimit=linear?2:config.maximumIterations;
    for(int iteration=1;iteration<=stageLimit;++iteration){
        try {
            std::array<double,8> stageSeconds{};auto stageStart=std::chrono::steady_clock::now();
            const auto mark=[&](int stage){const auto now=std::chrono::steady_clock::now();stageSeconds[stage]=std::chrono::duration<double>(now-stageStart).count();stageStart=now;};
            HydroNativeValues endpoint(l,guess,before.pressure());const auto arrival=reconstruction.evaluate(mesh,endpoint,accepted.geometryIdentity(),clock,accepted.boundary());
            const auto source=accepted.characteristicGeometry().linear(mesh,previous,arrival);
            if(flow)prescribedDirections(mesh,source,*flow);
            mark(0);
            const auto transport=[&]{try{
                return HydroMaterialTransport::prepareSpaceTime(mesh,accepted.publication(),source,config.liquidDensity,config.gasDensity,config.spaceTime,flow?flow->inletFirstJet:Interface::HermiteGalsBoundaryFirstJetQuery{},config.publication,preceding?&preceding->transport:nullptr);
            }catch(const std::exception& e){if(linear&&recoverableGeometry(e.what()))
                    throw HydroTwoPhaseStepRetry(HydroTwoPhaseRetryReason::Geometry,e.what());throw;}}();
            mark(1);
            const auto dual=HydroDualMaterial::prepare(mesh,transport,beginMass,makeMass(transport.midpointRate().geometry()),makeMass(transport.endGeometry()));
            mark(2);
            const auto& advector=transport.midpointRate().spatialAdvector();
            std::optional<HydroConvectionBoundary> inlet;
            if(flow)inlet.emplace(inletBoundary(mesh,dual,*flow));
            const auto convection=flow?HydroKepConvection::buildPrescribedFlow(mesh,dual,accepted.boundary(),*inlet,freeOutflow)
                :config.dissipation?HydroKepConvection::buildImpermeableControlled(mesh,dual,accepted.boundary(),*config.dissipation)
                :HydroKepConvection::buildImpermeable(mesh,dual,accepted.boundary());
            const auto cb=flow?prescribedBoundary(mesh,dual,*inlet,advector.values()):convectionBoundary(mesh,dual,advector);
            const auto loads=convection.apply(mesh,advector.values(),cb,dual.sources(),dual.timeLayer());
            mark(3);
            const auto& midpointPhase=transport.midpointRate().geometry();
            std::optional<HydroBalancedCsfLoad> capillary;
            std::optional<HydroInterfaceNormalLoad> sharp;
            std::vector<double> sharpRhs;
            const std::vector<double>* capillaryRhs=nullptr;
            double csfColorSeconds=0,csfProjectionSeconds=0;
            std::optional<Interface::HermiteGalsSurfaceVectorView> external;
            if(config.interfaceLoad){
                external.emplace(config.interfaceLoad(mesh,transport));
                external->validateFor(mesh,midpointPhase.publication(),transport.timeLayer());
            }
            if(config.surfaceTension>0||external){
                if(config.capillaryModel==HydroCapillaryModel::BalancedCsf){
                    const auto projectionStart=std::chrono::steady_clock::now();
                    // Select by the prepared source capability; construction
                    // failures are never caught as a switch to another load.
                    const bool projected=!flow && reconstruction.regularGeometry().usesCompatibleTensor()
                        && midpointPhase.publication().state().primary().geometryCells.size()==
                            static_cast<std::size_t>(mesh.cellIndexNum());
                    if(config.csfCurvature.forceMethod==HydroCsfForceMethod::HeightFunctionIncrement){
                        require(projected&&!external,"height-function-increment requires complete uniform XY closed motion without external traction");
                        const auto mixed=Interface::HermiteGalsAmrClosureSensitivity::buildUniformFirstJet(
                            mesh,midpointPhase.publication().state(),config.publication.volume.symmetryPlanes);
                        capillary=HydroBalancedCsfLoad::heightFunctionIncrementCapillary(mesh,pressure,midpointPhase,
                            config.surfaceTension,advector,&mixed,config.csfCurvature);
                        capillary->validateFor(mesh,pressure,dual.sources().midpoint,dual.timeLayer(),advector,&mixed);
                    }else if(config.csfCurvature.forceMethod==HydroCsfForceMethod::HeightFunction){
                        require(projected&&!external,"height-function requires complete uniform XY closed motion without external traction");
                        capillary=HydroBalancedCsfLoad::heightFunctionCapillary(mesh,pressure,midpointPhase,
                            config.surfaceTension,config.csfCurvature);
                        capillary->validateFor(mesh,pressure,dual.sources().midpoint,dual.timeLayer());
                    }else{
                        require(!midpointPhase.publication().state().geometry().snapshot().usesC2Macro() || projected,
                            "C2 macro CSF requires complete uniform projected interface motion");
                        if(projected){
                            const auto mixed=Interface::HermiteGalsAmrClosureSensitivity::buildUniformFirstJet(
                                mesh,midpointPhase.publication().state(),config.publication.volume.symmetryPlanes);
                            if(external)capillary=HydroBalancedCsfLoad::projectedMotionInterfaceLoad(mesh,pressure,midpointPhase,
                                config.surfaceTension,*external,advector,&mixed,config.csfCurvature);
                            else capillary=HydroBalancedCsfLoad::projectedMotionCapillary(mesh,pressure,midpointPhase,
                                config.surfaceTension,advector,&mixed,config.csfCurvature);
                            capillary->validateFor(mesh,pressure,dual.sources().midpoint,dual.timeLayer(),advector,&mixed);
                        }else{
                            if(external)capillary=HydroBalancedCsfLoad::adjointInterfaceLoad(mesh,pressure,midpointPhase,
                                config.surfaceTension,*external,advector.reconstruction(),advector.boundary(),config.csfCurvature);
                            else capillary=HydroBalancedCsfLoad::adjointCapillary(mesh,pressure,midpointPhase,
                                config.surfaceTension,advector.reconstruction(),advector.boundary(),config.csfCurvature);
                            capillary->validateFor(mesh,pressure,dual.sources().midpoint,dual.timeLayer(),
                                advector.reconstruction(),advector.boundary());
                        }
                        if(external)capillary->validateFor(mesh,pressure,dual.sources().midpoint,dual.timeLayer(),*external);
                    }
                    csfProjectionSeconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-projectionStart).count();
                    require(&capillary->phaseGeometry().cellVolumes()==&midpointPhase.cellVolumes(),
                        "different retained CSF midpoint phase owner");
                    capillaryRhs=&capillary->momentumRhsLoad();
                }else{
                    sharp=HydroInterfaceNormalLoad::capillary(mesh,pressure,midpointPhase,
                        config.surfaceTension,advector.reconstruction(),advector.boundary());
                    sharp->validateFor(mesh,pressure,dual.sources().midpoint,dual.timeLayer(),
                        advector.reconstruction(),advector.boundary());
                    require(&sharp->geometry().cellVolumes()==&midpointPhase.cellVolumes(),
                        "different retained sharp midpoint phase owner");
                    sharpRhs=sharp->balanceLoads();
                    for(double& value:sharpRhs)value=-value;
                    capillaryRhs=&sharpRhs;
                }
                require(capillaryRhs->size()==static_cast<std::size_t>(l.velocityCount()),
                    "different capillary load dimensions");
            }

            mark(4);
            const auto viscosity=config.liquidViscosity==0&&config.gasViscosity==0
                ?HydroPhaseViscosity::zero(mesh,midpointPhase)
                :HydroPhaseViscosity::build(mesh,midpointPhase,
                    Interface::HermiteGalsPhaseQuadrature::build(mesh,midpointPhase.publication()),config.liquidViscosity,config.gasViscosity);
            mark(5);
            HydroFrozenMidpointForces forces;forces.bodyAcceleration=config.bodyAcceleration;
            if(capillary)forces.interfaceLoad=*capillary;else if(sharp)forces.interfaceLoad=*sharp;
            const auto physicalBoundary=flow?HydroFrozenMidpointBoundary::prescribed(accepted.boundary(),flow->normalVelocity,flow->pressureOutlets)
                :HydroFrozenMidpointBoundary::closed(accepted.boundary());
            const auto frozen=[&]{try{
                return HydroFrozenMidpoint::prepare(mesh,dual,convection,loads,viscosity,pressure,before,std::move(forces),physicalBoundary,
                    linear?HydroFrozenMidpointScheme::FrozenConvection:HydroFrozenMidpointScheme::LaggedConvection,dual.sources(),dual.timeLayer());
            }catch(const std::exception& e){if(linear&&std::string(e.what()).find("velocity block is not positive definite")!=std::string::npos)
                    throw HydroTwoPhaseStepRetry(HydroTwoPhaseRetryReason::VelocityBlock,e.what());throw;}}();
            const auto& k=frozen.viscosityMatrix();const auto& physicalA=frozen.physicalVelocityMatrix();const auto& a=frozen.velocityMatrix();
            const Eigen::Map<const Eigen::VectorXd> u0(before.normalVelocity().data(),l.velocityCount());const auto& ku0=frozen.beginViscousAction();
            const auto& bodyLoads=frozen.bodyLoads();const auto& rhs=frozen.integralMomentumRhs();
            const auto& saddle=frozen.saddle();const auto& solved=frozen.solveResult();const auto& traction=saddle.pressureOutletLoads();
            mark(6);
            Trial trial(endpoint,source,transport,dual,convection,loads,viscosity,saddle,solved,a,rhs,bodyLoads,capillary,sharp,sharpRhs,frozen);
            HydroTwoPhaseIterationReport ir;ir.iteration=iteration;ir.timeScheme=config.timeScheme;ir.saddle=solved.report();ir.velocityUpdateRelative=difference(endpoint.normalVelocity(),solved.values().normalVelocity());
            if(preceding){ir.geometryCompared=true;ir.midpointGeometryUpdateRelative=geometryDifference(preceding->transport.midpointRate().geometry().publication(),transport.midpointRate().geometry().publication());ir.endGeometryUpdateRelative=geometryDifference(preceding->transport.endGeometry().publication(),transport.endGeometry().publication());ir.geometryUpdateRelative=std::max(ir.midpointGeometryUpdateRelative,ir.endGeometryUpdateRelative);}report.history.push_back(ir);
            const auto midpointSolved=midpoint(l,before,solved.values());const auto solveAdvector=reconstruction.evaluate(mesh,midpointSolved,midpointPhase.geometryIdentity(),midpointPhase.timeLayer(),accepted.boundary());
            const auto reassessed=convection.apply(mesh,midpointSolved,flow?prescribedBoundary(mesh,dual,*inlet,midpointSolved):convectionBoundary(mesh,dual,solveAdvector),dual.sources(),dual.timeLayer());
            const Eigen::Map<const Eigen::VectorXd> u1(solved.values().normalVelocity().data(),l.velocityCount());
            const Eigen::VectorXd au=a*u1,physicalAu=physicalA*u1;const auto cp=pressure.apply(mesh,solved.values());
            double linearError=0,linearScale=0;
            for(int u:saddle.freeVelocityDofs()){
                const double appliedRhs=pressureActive?rhs[u]+traction[u]:rhs[u];
                linearError=std::hypot(linearError,au[u]+cp[u]-appliedRhs);
                linearScale=std::hypot(linearScale,std::abs(au[u])+std::abs(cp[u])+std::abs(rhs[u])+(pressureActive?std::abs(traction[u]):0));}
            report.linearMomentumDefect=relative(linearError,linearScale);
            if(linear){
                const auto eta=interfaceMismatch(mesh,midpointPhase.publication(),advector,solveAdvector,dt);
                report.finalInterfaceMismatch=eta;if(iteration==1)report.interfaceMismatch=eta;
                report.history.back().interfaceMismatch=eta;
                if(iteration==2){
                    report.correctionRelative=correctionRelative(midpoint(l,before,preceding->solve.values()),midpointSolved,dual.midpointMass(),config.correctionReferenceVelocity);
                    report.corrected=true;
                }
            }
            double error=0,scale=0;for(int u:saddle.freeVelocityDofs()){
                long double value=static_cast<long double>(dual.beginMass().masses()[u])*u0[u]/dt-.5L*ku0[u]+reassessed.momentumLoads()[u];
                if(capillaryRhs)value+=(*capillaryRhs)[u];if(bodyActive)value+=bodyLoads[u];if(pressureActive)value+=traction[u];const double corrected=static_cast<double>(value);
                error=std::hypot(error,physicalAu[u]+cp[u]-corrected);scale=std::hypot(scale,std::abs(physicalAu[u])+std::abs(cp[u])+std::abs(corrected));}
            report.frozenOwnerMomentumDefect=relative(error,scale);
            const Eigen::Map<const Eigen::VectorXd> um(midpointSolved.normalVelocity().data(),l.velocityCount());const Eigen::VectorXd kum=k*um;
            std::vector<bool> fixedNormal(l.velocityCount(),true);
            for(int u:saddle.freeVelocityDofs())fixedNormal[u]=false;
            long double powerResidual=0,powerWork=0,boundaryReactionPower=0;
            for(int u=0;u<l.velocityCount();++u){
                const long double v=um[u],m0=dual.beginMass().masses()[u],m1=dual.endMass().masses()[u];
                const long double rate=(m1*u1[u]-m0*u0[u])/dt,viscous=kum[u],pressureForce=cp[u],load=reassessed.momentumLoads()[u];
                const long double force=capillaryRhs?(*capillaryRhs)[u]:0;
                const long double outletForce=pressureActive?traction[u]:0;
                if(bodyActive){powerResidual+=v*(rate+viscous+pressureForce-load-force-outletForce-bodyLoads[u]);
                    if(flow&&fixedNormal[u])boundaryReactionPower+=v*(rate+viscous+pressureForce-load-force-outletForce-bodyLoads[u]);
                }else{powerResidual+=v*(rate+viscous+pressureForce-load-force-outletForce);
                    if(flow&&fixedNormal[u])boundaryReactionPower+=v*(rate+viscous+pressureForce-load-force-outletForce);}
                powerWork+=std::abs(v*rate)+std::abs(v*viscous)+std::abs(v*pressureForce)+std::abs(v*load)+std::abs(v*force)+std::abs(v*outletForce);
                if(bodyActive)powerWork+=std::abs(v*bodyLoads[u]);
            }
            if(flow){powerResidual-=boundaryReactionPower;powerWork+=std::abs(boundaryReactionPower);}
            require(std::isfinite(powerResidual)&&std::isfinite(powerWork)&&powerWork>=0,"nonfinite physical momentum power work");
            if(powerWork==0){require(powerResidual==0,"nonzero physical momentum power with zero work");report.momentumPowerRelativeDefect=0;}
            else {const auto defect=std::abs(powerResidual)/powerWork;report.momentumPowerRelativeDefect=static_cast<double>(defect);
                require(std::isfinite(report.momentumPowerRelativeDefect)&&(defect==0||report.momentumPowerRelativeDefect>0),"unrepresentable physical momentum power defect");}
            mark(7);auto& completed=report.history.back();completed.stageSeconds=stageSeconds;
            completed.csfColorSeconds=csfColorSeconds;completed.csfProjectionSeconds=csfProjectionSeconds;
            const auto& timeReport=transport.spaceTimeReport();
            completed.samplingCost=timeReport.samplingCost;
            completed.geometrySamples=timeReport.geometrySamples;completed.acceptedPanels=timeReport.acceptedPanels;
            completed.initialPlanPanels=timeReport.initialPlanPanels;completed.deepestPanel=timeReport.deepestPanel;
            completed.endpointSupportPanels=timeReport.endpointSupportPanels;
            completed.reusedTimePlan=timeReport.reusedTimePlan;completed.timeLayer=transport.timeLayer();
            completed.frozenOwnerMomentumDefect=report.frozenOwnerMomentumDefect;
            completed.momentumPowerRelativeDefect=report.momentumPowerRelativeDefect;completed.linearMomentumDefect=report.linearMomentumDefect;
            if(config.onIteration)config.onIteration(completed);
            if(linear&&iteration==2&&report.correctionRelative>config.maximumCorrectionRelative)
                throw HydroTwoPhaseStepRetry(HydroTwoPhaseRetryReason::Correction,"[HydroTwoPhaseStep] linear midpoint correction exceeds configured bound");
            const bool chooseLinear=linear&&(iteration==2||(!report.startupCorrection&&!report.capillaryForcedCorrection
                &&(config.timeScheme==Scheme::ExtrapolatedMidpoint||report.interfaceMismatch<=config.interfaceMismatchTolerance)));
            if(chooseLinear){require(report.linearMomentumDefect<=1e-11,"linear midpoint frozen equation residual failed");
                require(report.frozenOwnerMomentumDefect<=1e-11&&report.momentumPowerRelativeDefect<=config.momentumPowerTolerance,
                    "linear midpoint physical momentum power gate failed");}
            if(chooseLinear||(!linear&&ir.geometryCompared&&ir.velocityUpdateRelative<=config.velocityUpdateTolerance&&ir.geometryUpdateRelative<=config.geometryUpdateTolerance&&report.frozenOwnerMomentumDefect<=1e-11&&report.momentumPowerRelativeDefect<=config.momentumPowerTolerance)){
                report.iterations=iteration;report.velocityUpdateRelative=ir.velocityUpdateRelative;report.geometryUpdateRelative=ir.geometryUpdateRelative;report.advectorGuessSolveRelative=ir.velocityUpdateRelative;
                const auto impulse=dual.consistencyImpulse(mesh,midpointSolved,dual.sources(),dual.timeLayer());
                long double timeScale=0;for(int u=0;u<l.velocityCount();++u){const long double b0=dual.beginMass().masses()[u],b1=dual.endMass().masses()[u],v0=u0[u],v1=u1[u],v=um[u],dm=b1-b0,dv=v1-v0;const long double e0=.5L*b0*v0*v0,e1=.5L*b1*v1*v1,inertia=v*(b1*dv+dm*v0),cross=.5L*v0*v1*dm,change=.5L*dm*v1*v1+.5L*b0*dv*(v1+v0);
                    report.kineticBegin+=e0;report.kineticEnd+=e1;report.kineticChange+=change;report.inertialImpulseWork+=inertia;report.temporalMassCrossTerm+=cross;timeScale+=std::abs(inertia)+std::abs(change)+std::abs(cross);
                    report.pressurePower+=v*cp[u];if(pressureActive)report.pressureOutletPower+=v*traction[u];report.viscousPower+=v*kum[u];report.convectionPower+=v*reassessed.momentumRates()[u];report.consistencyPower+=v*impulse[u]/dt;report.numericalDissipationPower+=v*reassessed.dissipationRates()[u];
                    if(capillaryRhs){report.interfacePower+=v*(*capillaryRhs)[u];
                        if(report.capillaryPowerAvailable)report.capillaryPower+=v*(*capillaryRhs)[u];}
                    if(bodyActive)report.bodyPower+=v*bodyLoads[u];
                    long double corrected=static_cast<long double>(dual.beginMass().masses()[u])*u0[u]/dt-.5L*ku0[u]+reassessed.momentumLoads()[u];if(capillaryRhs)corrected+=(*capillaryRhs)[u];
                    if(bodyActive)corrected+=bodyLoads[u];
                    auto physicalReaction=flow
                        ? (b1*v1-b0*v0)/dt+kum[u]+cp[u]-reassessed.momentumLoads()[u]-(capillaryRhs?(*capillaryRhs)[u]:0)
                        : physicalAu[u]+cp[u]-corrected;
                    if(flow&&bodyActive)physicalReaction-=bodyLoads[u];
                    if(pressureActive)physicalReaction-=traction[u];
                    if(fixedNormal[u])report.wallReactionPower+=v*(flow?physicalReaction:static_cast<long double>(au[u]+cp[u]-rhs[u]));
                    report.momentumPowerDefect+=v*physicalReaction;
                }report.inertialPower=report.inertialImpulseWork/dt;report.temporalIdentityDefect=relative(static_cast<double>(std::abs(report.inertialImpulseWork-report.kineticChange-report.temporalMassCrossTerm)),static_cast<double>(timeScale));
                if(flow)report.momentumPowerDefect-=report.wallReactionPower;
                if(bodyActive)require(std::isfinite(report.bodyPower),"nonfinite reported body power");
                for(int axis=0;axis<2;++axis){const auto& component=reassessed.report().components[axis];
                    report.boundaryMomentumRate[axis]=component.boundaryMomentumRate;
                    report.boundaryKineticPower+=component.boundaryKineticPower;
                    report.boundaryPhysicalKineticPower+=component.boundaryPhysicalKineticPower;
                    report.inletDissipationPower+=component.inletDissipationPower;}
                for(const auto& component:loads.report().components)report.sourceNumericalDissipationPower+=component.dissipationPower;
                if(capillaryRhs){
                    report.interfaceMeasureBegin=interfaceMeasure(transport.beginGeometry().publication());
                    report.interfaceMeasureMidpoint=interfaceMeasure(midpointPhase.publication());
                    report.interfaceMeasureEnd=interfaceMeasure(transport.endGeometry().publication());
                    report.surfaceEnergyBegin=static_cast<long double>(config.surfaceTension)*report.interfaceMeasureBegin;
                    report.surfaceEnergyEnd=static_cast<long double>(config.surfaceTension)*report.interfaceMeasureEnd;
                    report.surfaceEnergyChange=report.surfaceEnergyEnd-report.surfaceEnergyBegin;
                    if(report.capillaryPowerAvailable)report.capillaryEnergyDefect=report.surfaceEnergyChange+static_cast<long double>(dt)*report.capillaryPower;
                    require(std::isfinite(report.surfaceEnergyBegin)&&std::isfinite(report.surfaceEnergyEnd)&&
                        std::isfinite(report.surfaceEnergyChange)&&std::isfinite(report.interfacePower)
                        &&(!report.capillaryPowerAvailable||std::isfinite(report.capillaryEnergyDefect)),
                        "nonfinite reported capillary energy");
                }

                auto next=std::make_shared<HydroTwoPhaseState::Data>(transport.endGeometry(),solved.values(),geometry,reconstruction,accepted.boundary(),accepted.characteristicGeometry(),VelocityHistory{before,oldClock,accepted.geometryIdentity()},dual.endMass(),accepted.data_->flow);HydroTwoPhaseState state(next);state.validateFor(mesh);
                if(pressureActive){
                    // Certify the actual solved endpoint as well as the retained
                    // material guess. This checks direction; it creates no new F.
                    const auto solvedEndpoint=reconstruction.evaluate(mesh,state.values(),state.geometryIdentity(),clock,accepted.boundary());
                    prescribedDirections(mesh,accepted.characteristicGeometry().linear(mesh,previous,solvedEndpoint),*flow);
                }
                return HydroTwoPhaseStep(std::make_shared<Data>(accepted,std::move(state),std::move(trial),midpointSolved,reassessed,std::move(report)));
            }guess=solved.values().normalVelocity();preceding=std::move(trial);
        }catch(const HydroTwoPhaseStepRetry& e){std::ostringstream msg;msg<<"[HydroTwoPhaseStep] stage="<<iteration<<": "<<e.what();throw HydroTwoPhaseStepRetry(e.reason(),msg.str());
        }catch(const std::exception& e){std::ostringstream msg;msg<<"[HydroTwoPhaseStep] iteration="<<iteration<<": "<<e.what();throw std::runtime_error(msg.str());}
    }
    std::ostringstream msg;msg<<"[HydroTwoPhaseStep] corrector did not converge in "<<config.maximumIterations<<" iterations";if(!report.history.empty()){const auto& last=report.history.back();msg<<", velocity="<<last.velocityUpdateRelative<<", geometry="<<last.geometryUpdateRelative<<", frozen-owner momentum="<<report.frozenOwnerMomentumDefect<<", physical power="<<report.momentumPowerRelativeDefect;}throw std::runtime_error(msg.str());
}
const HydroTwoPhaseState& HydroTwoPhaseStep::beginState() const noexcept{return data_->begin;}
const HydroTwoPhaseState& HydroTwoPhaseStep::candidateState() const noexcept{return data_->candidate;}
const HydroNativeValues& HydroTwoPhaseStep::endpointGuess() const noexcept{return data_->trial.guess;}
const HydroNativeValues& HydroTwoPhaseStep::solvedMidpointValues() const noexcept{return data_->midpoint;}
const HydroCharacteristicSource& HydroTwoPhaseStep::characteristicSource() const noexcept{return data_->trial.source;}
const HydroMaterialTransport& HydroTwoPhaseStep::materialTransport() const noexcept{return data_->trial.transport;}
const HydroDualMaterial& HydroTwoPhaseStep::dualMaterial() const noexcept{return data_->trial.dual;}
const HydroKepConvection& HydroTwoPhaseStep::convection() const noexcept{return data_->trial.convection;}
const HydroKepConvectionResult& HydroTwoPhaseStep::convectionResult() const noexcept{return data_->trial.loads;}
const HydroKepConvectionResult& HydroTwoPhaseStep::reassessedConvectionResult() const noexcept{return data_->reassessed;}
bool HydroTwoPhaseStep::hasCapillaryLoad() const noexcept{return data_->trial.capillary.has_value()||data_->trial.sharp.has_value();}
HydroCapillaryModel HydroTwoPhaseStep::capillaryModel() const noexcept{return data_->report.capillaryModel;}
const std::vector<double>& HydroTwoPhaseStep::capillaryMomentumRhsLoad() const
{require(hasCapillaryLoad(),"step has no interface load");return data_->trial.capillary?data_->trial.capillary->momentumRhsLoad():
    data_->trial.sharpRhs;}
const HydroBalancedCsfLoad& HydroTwoPhaseStep::capillaryLoad() const
{require(hasCapillaryLoad(),"step has no interface load");require(data_->trial.capillary.has_value(),"sharpNormalStress step has no balancedCsf load");return *data_->trial.capillary;}
const HydroInterfaceNormalLoad& HydroTwoPhaseStep::sharpNormalStressLoad() const
{require(hasCapillaryLoad(),"step has no interface load");require(data_->trial.sharp.has_value(),"balancedCsf step has no sharpNormalStress load");return *data_->trial.sharp;}
const HydroPhaseViscosity& HydroTwoPhaseStep::viscosity() const noexcept{return data_->trial.viscosity;}
const HydroSinglePhaseSaddle& HydroTwoPhaseStep::saddle() const noexcept{return data_->trial.saddle;}
const HydroSaddleResult& HydroTwoPhaseStep::solveResult() const noexcept{return data_->trial.solve;}
const HydroTwoPhaseStep::Matrix& HydroTwoPhaseStep::velocityMatrix() const noexcept{return data_->trial.matrix;}
const std::vector<double>& HydroTwoPhaseStep::integralMomentumRhs() const noexcept{return data_->trial.rhs;}
const std::vector<long double>& HydroTwoPhaseStep::bodyLoads() const noexcept{return data_->trial.bodyLoads;}
const HydroTwoPhaseStepReport& HydroTwoPhaseStep::report() const noexcept{return data_->report;}
void HydroTwoPhaseStep::commit(const QuadtreeMesh& mesh,HydroTwoPhaseState& accepted) const
{require(accepted.sameIdentity(beginState()),"stale candidate cannot replace a different accepted handle");accepted.validateFor(mesh);candidateState().validateFor(mesh);accepted.data_=data_->candidate.data_;}
}


