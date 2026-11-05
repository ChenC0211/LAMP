#include "Apps/Electrohydrodynamics/ElectrohydrodynamicStep.h"
#include "Modules/Electric/ElectricMaxwellTrace.h"
#include "Modules/Interface/HermiteGalsFirstJetStep.h"
#include "Core/IO/DictionaryMesh.h"
#include "Core/Mesh/FaceSegmentGraph.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <Eigen/LU>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
namespace C=LAMP::Coupling;
namespace H=LAMP::Hydrodynamics;
namespace E=LAMP::Electric;
namespace I=LAMP::Interface;
using Clock=I::HermiteGalsTimeLayer;
using Position=Clock::Position;
constexpr long double pi=3.1415926535897932384626433832795L;

void check(bool value,const std::string& reason) {
    if(!value)throw std::runtime_error(reason);
}
bool sameDouble(double a,double b) {return std::memcmp(&a,&b,sizeof(double))==0;}
bool sameClock(Clock a,Clock b) {
    return a.step==b.step&&sameDouble(a.begin,b.begin)&&sameDouble(a.end,b.end)&&a.position==b.position;
}
bool sameBits(const std::vector<double>& a,const std::vector<double>& b) {
    return a.size()==b.size()&&(a.empty()||std::memcmp(a.data(),b.data(),a.size()*sizeof(double))==0);
}
double relative(long double gap,long double scale) {
    check(std::isfinite(gap)&&std::isfinite(scale)&&scale>=0,"nonfinite independent residual");
    return scale>0?double(std::abs(gap)/scale):(gap==0?0:std::numeric_limits<double>::infinity());
}
void scaled(long double actual,long double expected,long double scale,double tolerance,const char* reason) {
    check(relative(actual-expected,scale)<=tolerance,reason);
}
template<class F> void rejected(F action,const std::string& marker) {
    try {action();}catch(const std::exception& e) {
        check(marker.empty()||std::string(e.what()).find(marker)!=std::string::npos,
            "unexpected rejection: "+std::string(e.what()));
        std::cout<<"EJ1_REJECT,"<<marker<<','<<e.what()<<'\n';return;
    }
    throw std::runtime_error("missing joint rejection: "+marker);
}
long double faceArea(const QuadtreeMesh& mesh,const H::HydroDofLayout& layout,int u) {
    const auto f=layout.domain().face(mesh,layout.nativeFace(u));
    const auto& edge=f.horizontal?mesh.hEdges(f.edgeLabel):mesh.vEdges(f.edgeLabel);
    auto a=mesh.verts(edge.edgeToVert[0]).coord(),b=mesh.verts(edge.edgeToVert[1]).coord();
    if(f.horizontal?a.x>b.x:a.y>b.y)std::swap(a,b);
    return f.horizontal?pi*(static_cast<long double>(b.x)*b.x-static_cast<long double>(a.x)*a.x)
        :2*pi*a.x*(b.y-a.y);
}
H::HydroPhysicalDivergence metrics(const QuadtreeMesh& mesh) {
    FaceSegmentGraph graph;graph.build(mesh);
    const auto layout=H::HydroDofLayout::build(mesh,NativeLeafDomain::build(mesh,graph,true));
    return H::HydroPhysicalDivergence::build(mesh,H::HydroIntegralDivergence::build(mesh,layout));
}
H::HydroNativeValues initialVelocity(const QuadtreeMesh& mesh,const H::HydroDofLayout& layout) {
    const auto psi=[](Vector2D p) {
        const double r2=p.x*p.x,z=p.y*(1-p.y);
        return 16*r2*(1-r2)*(1-r2)*z*z*(1+.7*p.y);
    };
    std::vector<double> velocity(layout.velocityCount());
    for(int u=0;u<layout.velocityCount();++u) {
        const auto f=layout.domain().face(mesh,layout.nativeFace(u));
        const auto& edge=f.horizontal?mesh.hEdges(f.edgeLabel):mesh.vEdges(f.edgeLabel);
        auto a=mesh.verts(edge.edgeToVert[0]).coord(),b=mesh.verts(edge.edgeToVert[1]).coord();
        if(f.horizontal?a.x>b.x:a.y>b.y)std::swap(a,b);
        const auto area=faceArea(mesh,layout,u);
        velocity[u]=area>0?double(2*pi*(f.horizontal?psi(a)-psi(b):psi(b)-psi(a))/area):0;
    }
    return {layout,velocity,std::vector<double>(layout.pressureCount())};
}

// Capture the unique accepted owners and every FirstJet/closure/native payload.
// The external Field can be checked separately during the intentional edit.
struct AcceptedCopy {
    H::HydroTwoPhaseState hydro;
    I::HermiteGalsPrimaryState primary;
    I::HermiteGalsClosure closure;
    std::vector<double> velocity,pressure,charge,previousCharge;
    Clock clock,representationClock;
    std::optional<Clock> previousClock;
    const I::HermiteGalsPrimaryState* primaryOwner;
    const LevelSetHermiteGeometry::Snapshot* geometryOwner;
    const E::SurfaceChargeField* fieldOwner;
    const E::ElectricFixedGeometryChargeSolver* solverOwner;
    explicit AcceptedCopy(const C::ElectrohydrodynamicState& state)
        :hydro(state.hydro()),primary(state.publication().state().primary()),
        closure(state.publication().state().closure()),
        velocity(state.hydro().values().normalVelocity()),pressure(state.hydro().values().pressure()),
        charge(state.electric().integratedCharge()),
        previousCharge(state.electric().previousIntegratedCharge()?*state.electric().previousIntegratedCharge():std::vector<double>{}),
        clock(state.timeLayer()),representationClock(state.electric().geometrySourceTimeLayer()),
        previousClock(state.electric().previousTimeLayer()?std::optional<Clock>(*state.electric().previousTimeLayer()):std::nullopt),
        primaryOwner(&state.publication().state().primary()),
        geometryOwner(&state.publication().state().geometry().snapshot()),
        fieldOwner(&state.electric().field(0)),solverOwner(&state.electric().solver()) {}
};
void unchanged(const C::ElectrohydrodynamicState& state,const AcceptedCopy& old,bool checkField=true) {
    check(state.hydro().sameIdentity(old.hydro)&&sameClock(state.timeLayer(),old.clock)
        &&sameClock(state.electric().timeLayer(),old.clock)
        &&sameClock(state.electric().geometrySourceTimeLayer(),old.representationClock),
        "failed joint preparation/commit changed accepted owner or complete clock");
    const auto& p=state.publication().state().primary();
    const auto& c=state.publication().state().closure();
    check(&p==old.primaryOwner&&&state.publication().state().geometry().snapshot()==old.geometryOwner
        &&&state.electric().solver()==old.solverOwner,"failed joint operation replaced accepted geometry/solver owner");
    check(sameBits(state.hydro().values().normalVelocity(),old.velocity)
        &&sameBits(state.hydro().values().pressure(),old.pressure)
        &&sameBits(state.electric().integratedCharge(),old.charge),
        "failed joint operation changed u/p/Q receipts");
    check(bool(state.electric().previousTimeLayer())==old.previousClock.has_value()
        &&bool(state.electric().previousIntegratedCharge())==old.previousClock.has_value()
        &&(!old.previousClock||(sameClock(*state.electric().previousTimeLayer(),*old.previousClock)
            &&sameBits(*state.electric().previousIntegratedCharge(),old.previousCharge))),
        "failed joint operation changed true accepted charge history");
    check(p.measure==old.primary.measure&&p.transport==old.primary.transport
        &&p.geometryCells==old.primary.geometryCells&&p.phaseVolumes.empty()&&old.primary.phaseVolumes.empty()
        &&p.vertices.size()==old.primary.vertices.size()&&p.edges.size()==old.primary.edges.size(),
        "failed joint operation changed pure FirstJet topology");
    for(std::size_t v=0;v<p.vertices.size();++v) {
        const auto& a=p.vertices[v];const auto& b=old.primary.vertices[v];
        check(a.vertex==b.vertex&&sameDouble(a.phi,b.phi)&&sameDouble(a.dx,b.dx)&&sameDouble(a.dy,b.dy),"failed joint operation changed vertex FirstJet");
    }
    for(std::size_t e=0;e<p.edges.size();++e) {
        const auto& a=p.edges[e];const auto& b=old.primary.edges[e];
        check(a.vertices==b.vertices&&sameDouble(a.phi,b.phi)&&sameDouble(a.dn,b.dn),"failed joint operation changed edge FirstJet");
    }
    check(sameBits(c.vertexMixedDerivatives,old.closure.vertexMixedDerivatives)
        &&sameBits(c.cellBubbleCoefficients,old.closure.cellBubbleCoefficients),"failed joint operation changed closure");
    if(checkField)check(&state.electric().field(0)==old.fieldOwner
        &&sameBits(state.electric().field(0).transport().panelCharge(),old.charge),"failed joint operation changed Field owner/Q");
}

// Original Q/V equations solved directly, independently of the matrix-free Schur
// implementation. Coefficients come from physical nodes, not report receipts.
void electricBlock(const E::ElectricMovingGeometryChargeStep& moving,const E::ElectricPhaseProblem& problem,
    const std::array<double,2>& conductivity) {
    const auto& step=moving.step();const auto& old=moving.mappedIntegratedCharge();
    const auto* previous=moving.mappedPreviousIntegratedCharge();
    const double h=step.timeLayer().time()-moving.sourceTimeLayer().time();
    const double d=moving.previousSourceTimeLayer()?moving.sourceTimeLayer().time()-moving.previousSourceTimeLayer()->time():0;
    check(h>0&&(!previous||d>0),"invalid original Q/V physical time nodes");
    const double a0=previous?(2*h+d)/(h+d):1,a1=previous?-(h+d)/d:-1,a2=previous?h*h/(d*(h+d)):0;
    const auto& report=step.report();
    scaled(report.dt,h,h,1e-14,"stage increment differs from actual physical increment");
    if(previous)scaled(report.previousDt,d,d,1e-14,"stage used midpoint history increment");
    check(report.temporalOrder==(previous?2:1),"stage lost End history or invented startup history");
    scaled(report.a0,a0,std::abs(a0),1e-14,"original physical-node a0");
    scaled(report.a1,a1,std::abs(a1),1e-14,"original physical-node a1");
    scaled(report.a2,a2,std::max(1.,std::abs(a2)),1e-14,"original physical-node a2");
    const auto& op=step.current().momentOperator();
    const auto zero=E::solveElectricPhaseElectrostatics(step.potential().geometry(),problem,op.apply(std::vector<double>(op.panelCount(),0)));
    const auto current=zero.preparePanelCurrent(op,conductivity);
    const Eigen::MatrixXd physical=op.interfaceInjectionMatrix()*op.projectionMatrix();
    const int nu=int(zero.skeletonSolution().size()),nq=int(op.panelCount());
    Eigen::MatrixXd injection=Eigen::MatrixXd::Zero(nu,nq);int offset=0;
    for(std::size_t face=0;face<zero.geometry().faces().size();++face) {
        const int side=zero.geometry().faces()[face].boundarySide;
        const bool prescribed=side>=0&&problem.boundaries[std::size_t(side)].type==E::ElectricCartesianBoundary::Type::Dirichlet;
        if(!prescribed) {check(offset+3<=nu,"original free electric layout");
            injection.middleRows(offset,3)=physical.middleRows(3*Eigen::Index(face),3);offset+=3;}
    }
    check(offset==nu,"original charge injection missed free modes");
    Eigen::MatrixXd block=Eigen::MatrixXd::Zero(nu+nq,nu+nq);
    block.topLeftCorner(nu,nu)=Eigen::MatrixXd(zero.skeletonMatrix());
    block.topRightCorner(nu,nq)=-injection;
    block.bottomLeftCorner(nq,nu)=-h*Eigen::MatrixXd(current.matrix());
    block.bottomRightCorner(nq,nq)=a0*Eigen::MatrixXd::Identity(nq,nq);
    Eigen::VectorXd rhs(nu+nq);rhs.head(nu)=Eigen::Map<const Eigen::VectorXd>(zero.skeletonRhs().data(),nu);
    for(int q=0;q<nq;++q)rhs[nu+q]=-a1*old[std::size_t(q)]-(previous?a2*(*previous)[std::size_t(q)]:0)+h*current.offset()[std::size_t(q)];
    const Eigen::FullPivLU<Eigen::MatrixXd> direct(block);
    check(direct.isInvertible(),"original Q/V block is singular");
    const Eigen::VectorXd reference=direct.solve(rhs);
    check((block*reference-rhs).norm()/(rhs.norm()+block.norm()*reference.norm())<=1e-11,"original Q/V reference residual");
    const Eigen::VectorXd referenceCurrent=Eigen::MatrixXd(current.matrix())*reference.head(nu)
        +Eigen::Map<const Eigen::VectorXd>(current.offset().data(),nq);
    double qGap=0,vGap=0,iGap=0,panelGap=0;long double total=0,totalScale=0;
    for(int q=0;q<nq;++q) {
        const auto k=std::size_t(q);
        qGap=std::max(qGap,std::abs(step.integratedCharge()[k]-reference[nu+q])/(1+std::abs(reference[nu+q])));
        iGap=std::max(iGap,std::abs(step.current().panelCurrent()[k]-referenceCurrent[q])/(1+std::abs(referenceCurrent[q])));
        const long double history=previous?a2*(*previous)[k]:0;
        const long double residual=a0*step.integratedCharge()[k]+a1*old[k]+history-h*step.current().panelCurrent()[k];
        const long double scale=std::abs(a0*step.integratedCharge()[k])+std::abs(a1*old[k])+std::abs(history)+std::abs(h*step.current().panelCurrent()[k]);
        panelGap=std::max(panelGap,relative(residual,scale));total+=residual;totalScale+=scale;
    }
    for(int v=0;v<nu;++v)vGap=std::max(vGap,std::abs(step.potential().skeletonSolution()[std::size_t(v)]-reference[v])/(1+std::abs(reference[v])));
    check(qGap<=1e-9&&vGap<=1e-9&&iGap<=1e-9&&panelGap<=1e-10&&relative(total,totalScale)<=1e-10,
        "actual stage/End Q/V/I differs from original coupled equations");
    check(step.current().conductivityByPhase()==conductivity&&step.current().epsilonByPhase()==problem.epsilonByPhase,
        "actual electric candidate reused foreign material");
    std::cout<<"EJ1_QV_BLOCK,"<<step.timeLayer().step<<','<<int(step.timeLayer().position)<<','<<h<<','<<d
        <<','<<qGap<<','<<vGap<<','<<iGap<<','<<panelGap<<','<<relative(total,totalScale)<<'\n';
}
void coeval(const QuadtreeMesh& mesh,const E::ElectricMovingGeometryChargeStep& moving,Clock clock) {
    const auto& step=moving.step();const auto& pub=step.potential().geometry().publication();
    const auto& op=step.current().momentOperator();const auto maxwell=E::ElectricMaxwellTrace::prepare(step.potential());
    check(sameClock(step.timeLayer(),clock)&&sameClock(pub.state().timeLayer(),clock)
        &&sameClock(step.geometrySourceTimeLayer(),clock)&&sameClock(moving.geometrySourceTimeLayer(),clock)
        &&sameClock(step.potential().timeLayer(),clock)&&sameClock(step.current().timeLayer(),clock)
        &&sameClock(op.timeLayer(),clock)&&sameClock(maxwell.timeLayer(),clock)
        &&sameClock(maxwell.geometrySourceTimeLayer(),clock),"stage/End Q/V/I/Maxwell full clock mismatch");
    check(op.geometry().geometryIdentity()==pub.state().identity()
        &&moving.field(0).geometry().usesSource(pub.state().geometry().snapshot())
        &&sameBits(moving.field(0).transport().panelCharge(),step.integratedCharge()),
        "stage/End field/P_Q uses foreign publication or Q");
    maxwell.validateSource(mesh,pub,clock);
    long double tangent=0;
    for(auto face:maxwell.interfaceFaces())for(const auto& q:maxwell.quadrature(face)) {
        scaled(q.tangentialTraction,q.chargeDensity*q.tangentialField,
            std::max(1.,std::abs(q.tangentialTraction)),1e-12,"actual q_h Et traction identity");
        tangent+=q.physicalWeight*std::abs(q.tangentialTraction);
    }
    check(tangent>0,"signed real Q lost tangential Maxwell stress");
}

void characteristic(const QuadtreeMesh& mesh,const H::HydroTwoPhaseStep& step,double dt) {
    const auto& source=step.characteristicSource();const auto& curve=step.beginState().publication().topology().curves.front();
    check(!curve.graphs().empty(),"missing source native interface");
    H::HydroCharacteristicConfig config;config.maxStep=.00025;
    std::vector<Vector2D> axisPoints;std::optional<Vector2D> interior;
    for(std::size_t g=0;g<curve.graphs().size();++g) {
        for(double fraction:{0.,1.}) {
            const auto p=curve.frameAt(int(g),fraction).point;
            if(p.x==0&&std::none_of(axisPoints.begin(),axisPoints.end(),[&](Vector2D q){return q.y==p.y;}))axisPoints.push_back(p);
        }
        for(double fraction:{.37,.23,.61}) {
            const auto p=curve.frameAt(int(g),fraction).point;
            if(!interior&&p.x>0&&source.geometry().locatePieces(p).size()==1)interior=p;
        }
    }
    check(axisPoints.size()==2&&interior.has_value(),"actual material map lacks two poles or an interior chart witness");
    const std::array<Vector2D,3> points{{axisPoints[0],*interior,axisPoints[1]}};
    int poles=0;double movement=0,roundtrip=0;
    // The order-three forward/backward maps have fourth-order local error.
    // This short physical window uses a fixed dt^4 bound plus FP roundoff.
    const double bound=64*std::pow(dt,4)+4096*std::numeric_limits<double>::epsilon();
    for(auto p:points) {
        const auto forward=source.traceForward(p,step.beginState().timeLayer().time(),step.candidateState().timeLayer().time(),config);
        check(!forward.boundary()&&forward.time()==step.candidateState().timeLayer().time(),"actual forward material trace stopped early");
        if(p.x==0) {++poles;check(forward.point().x==0,"actual material pole left the radial axis");}
        const auto backward=source.traceBackward(forward.point(),forward.time(),step.beginState().timeLayer().time(),config);
        check(!backward.boundary()&&backward.time()==step.beginState().timeLayer().time(),"actual roundtrip trace stopped early");
        roundtrip=std::max(roundtrip,std::hypot(backward.point().x-p.x,backward.point().y-p.y));
        movement=std::max(movement,std::hypot(forward.point().x-p.x,forward.point().y-p.y));
    }
    check(poles==2&&movement>0&&roundtrip<=bound,"actual material map pole/motion/physical roundtrip witness");
    const auto mapped=source.traceForward(*interior,step.beginState().timeLayer().time(),step.candidateState().timeLayer().time(),config);
    check(mapped.jacobian().has_value(),"actual interior forward map lacks its differentiated kernel");
    constexpr double epsilon=1e-6;std::array<Vector2D,2> derivative;
    for(int axis=0;axis<2;++axis) {
        auto plus=*interior,minus=*interior;
        if(axis==0){plus.x+=epsilon;minus.x-=epsilon;}else{plus.y+=epsilon;minus.y-=epsilon;}
        const auto a=source.traceForward(plus,step.beginState().timeLayer().time(),step.candidateState().timeLayer().time(),config);
        const auto b=source.traceForward(minus,step.beginState().timeLayer().time(),step.candidateState().timeLayer().time(),config);
        check(!a.boundary()&&!b.boundary(),"central forward-map witness hit a physical boundary");
        derivative[axis]={(a.point().x-b.point().x)/(2*epsilon),(a.point().y-b.point().y)/(2*epsilon)};
    }
    const auto j=*mapped.jacobian();
    const double jGap=std::hypot(std::hypot(derivative[0].x-j.xx,derivative[1].x-j.xy),
        std::hypot(derivative[0].y-j.yx,derivative[1].y-j.yy));
    const double jScale=std::max(1.,std::hypot(std::hypot(j.xx,j.xy),std::hypot(j.yx,j.yy)));
    const double fdBound=16*epsilon*epsilon+1024*std::numeric_limits<double>::epsilon()/epsilon;
    check(jGap/jScale<=fdBound,"forward Jacobian differs from central differences of the same actual map");
    (void)mesh;
    std::cout<<"EJ1_MATERIAL_MAP,"<<movement<<','<<roundtrip<<','<<bound<<",axis_poles="<<poles
        <<",forward_J_relative="<<jGap/jScale<<",forward_J_bound="<<fdBound<<'\n';
}
void momentum(const QuadtreeMesh& mesh,const H::HydroTwoPhaseStep& step) {
    const auto& dual=step.dualMaterial();const auto& layout=dual.layout();
    const auto& m0=dual.beginMass().masses();const auto& m1=dual.endMass().masses();
    const auto& u0=step.beginState().values().normalVelocity();const auto& u1=step.candidateState().values().normalVelocity();
    const auto& bar=step.solvedMidpointValues().normalVelocity();const auto& pressure=step.candidateState().values().pressure();
    const auto& force=step.capillaryMomentumRhsLoad();const auto& physical=step.reassessedConvectionResult();
    const auto frozen=step.convection().frozenMomentumOperator(mesh,dual.sources(),dual.timeLayer());
    const Eigen::VectorXd hu0=frozen*Eigen::Map<const Eigen::VectorXd>(u0.data(),Eigen::Index(u0.size()));
    const auto& report=step.report();const double dt=report.deltaTime;
    const auto& k=step.viscosity().viscosity().matrix();
    check(k.nonZeros()==0&&!step.viscosity().hasQuadrature()&&report.viscousPower==0,"joint step lost exact unsampled mu=0");
    long double rhsGap=0,rhsScale=0,gap=0,scale=0,power=0,powerScale=0,loadPower=0,loadScale=0,tangent=0;
    long double kinetic0=0,kinetic1=0,inertial=0,cross=0,temporalScale=0;
    for(int u=0;u<layout.velocityCount();++u) {
        const auto f=layout.domain().face(mesh,layout.nativeFace(u));
        const long double pm=f.negativeCell<0?0:pressure[layout.pressureDof(f.negativeCell)];
        const long double pp=f.positiveCell<0?0:pressure[layout.pressureDof(f.positiveCell)];
        const long double cp=faceArea(mesh,layout,u)*(pp-pm);
        const long double rate=(static_cast<long double>(m1[u])*u1[u]-static_cast<long double>(m0[u])*u0[u])/dt;
        const long double residual=rate+cp-physical.momentumLoads()[u]-force[u]-step.bodyLoads()[u];
        const long double work=(std::abs(static_cast<long double>(m1[u])*u1[u])+std::abs(static_cast<long double>(m0[u])*u0[u]))/dt
            +std::abs(cp)+std::abs(physical.momentumLoads()[u])+std::abs(force[u])+std::abs(step.bodyLoads()[u]);
        const long double rhs=static_cast<long double>(m0[u])*u0[u]/dt-.5L*hu0[u]+force[u]+step.bodyLoads()[u];
        rhsGap=std::hypot(rhsGap,step.integralMomentumRhs()[u]-rhs);rhsScale=std::hypot(rhsScale,std::abs(rhs)+std::abs(force[u]));
        if(f.negativeCell>=0&&f.positiveCell>=0) {gap=std::hypot(gap,residual);scale=std::hypot(scale,work);}
        power+=bar[u]*residual;powerScale+=std::abs(bar[u]*rate)+std::abs(bar[u]*cp)
            +std::abs(bar[u]*physical.momentumLoads()[u])+std::abs(bar[u]*force[u])+std::abs(bar[u]*step.bodyLoads()[u]);
        loadPower+=bar[u]*force[u];loadScale+=std::abs(bar[u]*force[u]);
        tangent+=std::abs(step.capillaryLoad().tangentialRhsLoad()[u]);
        kinetic0+=.5L*m0[u]*u0[u]*u0[u];kinetic1+=.5L*m1[u]*u1[u]*u1[u];
        inertial+=bar[u]*(static_cast<long double>(m1[u])*u1[u]-static_cast<long double>(m0[u])*u0[u]);
        cross+=.5L*u0[u]*u1[u]*(m1[u]-m0[u]);
        temporalScale+=std::abs(bar[u]*(static_cast<long double>(m1[u])*u1[u]-static_cast<long double>(m0[u])*u0[u]))
            +.5L*(std::abs(static_cast<long double>(m0[u])*u0[u]*u0[u])+std::abs(static_cast<long double>(m1[u])*u1[u]*u1[u]));
        scaled(bar[u],.5L*(u0[u]+u1[u]),std::abs(.5L*u0[u])+std::abs(.5L*u1[u]),1e-12,"actual solved midpoint velocity");
    }
    check(relative(rhsGap,rhsScale)<=1e-12&&relative(gap,scale)<=1e-11&&relative(power,powerScale)<=1e-11,
        "actual full electric/capillary load is absent from native momentum or physical power");
    scaled(report.interfacePower,loadPower,loadScale,1e-12,"complete interface power receipt");
    scaled(report.momentumPowerRelativeDefect*powerScale,std::abs(power),powerScale,1e-12,"actual physical power residual receipt");
    scaled(report.kineticBegin,kinetic0,kinetic0,1e-12,"true variable-density initial kinetic energy");
    scaled(report.kineticEnd,kinetic1,kinetic1,1e-12,"true variable-density End kinetic energy");
    scaled(inertial-(kinetic1-kinetic0),cross,temporalScale,1e-12,"true variable-mass temporal identity");
    check(report.hasExternalInterfaceLoad&&!report.capillaryPowerAvailable&&std::isnan(report.capillaryPower)
        &&std::isnan(report.capillaryEnergyDefect)&&std::isfinite(report.interfacePower)
        &&std::isfinite(report.surfaceEnergyBegin)&&std::isfinite(report.surfaceEnergyEnd)&&tangent>0&&loadScale>0,
        "mixed load misreported total work as pure capillary work or erased tangent force");
    check(u1!=u0&&step.beginState().publication().state().closure().cellBubbleCoefficients
        !=step.candidateState().publication().state().closure().cellBubbleCoefficients,"joint window froze real velocity/FirstJet geometry");
    for(int axis:step.beginState().phaseGeometry().metrics().axisVelocityDofs())
        check(u1[axis]==0&&step.materialTransport().faceIntegrals()[axis]==0,"joint step leaked true radial-axis velocity/material flux");
    std::cout<<"EJ1_MOMENTUM,"<<report.surfaceTension<<','<<relative(gap,scale)<<','<<relative(power,powerScale)
        <<",interface_power="<<double(loadPower)<<",tangent_l1="<<double(tangent)<<'\n';
}
}
int main()try {
    std::cout<<std::unitbuf<<std::setprecision(17);
    DictionaryMesh meshInput;meshInput.Nx=meshInput.Ny=4;meshInput.maxInitialIterations=0;meshInput.maxCellNum=512;
    QuadtreeMesh mesh(meshInput);
    const auto publication=I::initializeHermiteGalsFirstJet(mesh,true,[](Vector2D p) {
        const double z=p.y-.51;return I::HermiteGalsBoundaryFirstJet{p.x*p.x+z*z-.27*.27,{2*p.x,2*z}};
    });
    const auto phase=H::HydroPhaseGeometry::build(mesh,metrics(mesh),publication);
    const auto hydro=H::HydroTwoPhaseState::create(mesh,phase,initialVelocity(mesh,phase.layout()));
    const auto electricGeometry=E::ElectricPhaseGeometry::build(mesh,publication);
    E::SurfaceChargeField initialField(mesh,publication.topology().curves.front(),0);
    auto q0=initialField.transport().panelCharge();bool positive=false,negative=false;
    for(std::size_t q=0;q<q0.size();++q) {
        q0[q]=initialField.transport().panelMeasure()[q]*(std::sin(1.7*(q+1))+.2+.01*q);
        positive=positive||q0[q]>0;negative=negative||q0[q]<0;
    }
    check(positive&&negative,"real RZ fixture lacks signed nonuniform integrated charge");initialField.assignIntegratedCharge(q0);
    C::ElectrohydrodynamicStepConfig config;
    config.hydro.liquidDensity=1000;config.hydro.gasDensity=1;config.hydro.liquidViscosity=config.hydro.gasViscosity=0;
    config.hydro.surfaceTension=1;config.hydro.spaceTime.initialPanels=1;config.hydro.spaceTime.characteristic.maxStep=.00025;
    config.electricProblem.epsilonByPhase={{5,1}};config.conductivityByPhase={{.7,.03}};
    for(int side=0;side<4;++side) {
        config.electricProblem.boundaries[side].type=side==0?E::ElectricCartesianBoundary::Type::Neumann:E::ElectricCartesianBoundary::Type::Dirichlet;
        config.electricProblem.boundaries[side].value=side==0?std::function<double(Vector2D)>([](Vector2D){return 0.;})
            :std::function<double(Vector2D)>([](Vector2D p){return .05*p.y;});
    }
    const auto moments=E::ElectricSurfaceChargeMomentOperator::prepare(electricGeometry,{&initialField},electricGeometry.timeLayer());
    const auto potential=E::solveElectricPhaseElectrostatics(electricGeometry,config.electricProblem,moments.apply(q0));
    const auto solver=E::ElectricFixedGeometryChargeSolver::prepare(potential,moments,config.conductivityByPhase);
    auto electric=E::ElectricSurfaceChargeTimeState::initialize(solver,{&initialField},electricGeometry.timeLayer());
    auto state=C::ElectrohydrodynamicState::create(mesh,hydro,std::move(electric));
    check(config.hydro.timeScheme==H::HydroTwoPhaseTimeScheme::AdaptiveLinearMidpoint
        &&config.hydro.interfaceMismatchTolerance==.01&&config.hydro.maximumCorrectionRelative==.25
        &&config.hydro.momentumPowerTolerance==1e-11,"joint contract changed adopted AC gates");
    const Clock firstClock{1,0,.00025,Position::End};const AcceptedCopy initialCopy(state);
    for(const auto failurePosition:{Position::Midpoint,Position::End}) {
        auto failing=config;int callbacks=0,midpoints=0;
        failing.onElectricCandidate=[&](const E::ElectricMovingGeometryChargeStep& candidate) {
            ++callbacks;unchanged(state,initialCopy);
            if(candidate.step().timeLayer().position==Position::Midpoint)++midpoints;
            if(candidate.step().timeLayer().position==failurePosition)throw std::runtime_error("EJ1_OBSERVER_FAILURE");
        };
        rejected([&]{C::ElectrohydrodynamicStep::prepare(mesh,state,firstClock,failing);},"EJ1_OBSERVER_FAILURE");
        unchanged(state,initialCopy);
        check(callbacks>0&&(failurePosition!=Position::End||midpoints==2),"End failure did not follow real startup A+C candidates");
    }
    int completedHydroStages=0;auto failingHydro=config;
    failingHydro.hydro.onIteration=[&](const H::HydroTwoPhaseIterationReport& report) {
        ++completedHydroStages;check(report.saddle.freeVelocityCount>0&&report.geometrySamples>=2,
            "late joint failure did not reach actual material/saddle");unchanged(state,initialCopy);
        throw std::runtime_error("EJ1_COMPLETED_HYDRO_STAGE_FAILURE");
    };
    rejected([&]{C::ElectrohydrodynamicStep::prepare(mesh,state,firstClock,failingHydro);},"EJ1_COMPLETED_HYDRO_STAGE_FAILURE");
    check(completedHydroStages==1,"late joint stage failure was skipped/repeated");unchanged(state,initialCopy);

    for(int acceptedStep=1;acceptedStep<=2;++acceptedStep) {
        const AcceptedCopy before(state);const double dt=acceptedStep==1?.00025:.0003;
        const Clock end{acceptedStep,before.clock.time(),before.clock.time()+dt,Position::End};
        const Clock mid{end.step,end.begin,end.end,Position::Midpoint};
        auto observed=config;if(acceptedStep==2)observed.hydro.surfaceTension=0;
        int callbacks=0;observed.onElectricCandidate=[&](const E::ElectricMovingGeometryChargeStep&) {
            ++callbacks;unchanged(state,before);
        };
        const auto step=C::ElectrohydrodynamicStep::prepare(mesh,state,end,observed);
        unchanged(state,before);
        const auto& h=step.hydroStep();const auto& stages=step.electricStages();
        check(!stages.empty()&&stages.size()==std::size_t(h.report().iterations)
            &&callbacks==int(stages.size())+1&&h.report().history.size()==stages.size(),"actual Hydro/electric A+C stage inventory mismatch");
        check(h.report().startupCorrection==(acceptedStep==1)&&h.report().usedVelocityHistory==(acceptedStep==2)
            &&(acceptedStep!=1||stages.size()==2),"joint AC startup or accepted velocity history was lost");
        for(const auto& stage:stages) {
            coeval(mesh,stage,mid);electricBlock(stage,observed.electricProblem,observed.conductivityByPhase);
            check(sameClock(stage.sourceTimeLayer(),before.clock)
                &&bool(stage.previousSourceTimeLayer())==before.previousClock.has_value()
                &&(!before.previousClock||sameClock(*stage.previousSourceTimeLayer(),*before.previousClock)),
                "stage advanced or retagged physical accepted End history");
        }
        const auto& last=stages.back();const auto& load=h.capillaryLoad();
        check(last.step().potential().geometry().geometryIdentity()==h.materialTransport().midpointRate().geometry().geometryIdentity()
            &&sameClock(load.timeLayer(),mid)&&load.hasExternalSurfaceVector(),"final actual stage geometry/load is not coeval");
        const auto maxwell=E::ElectricMaxwellTrace::prepare(last.step().potential());
        const auto& midPub=h.materialTransport().midpointRate().geometry().publication();
        for(std::size_t g=0;g<midPub.topology().curves.front().graphs().size();++g) {
            // Authenticate at an actual accepted electric/publication quadrature node.
            const double pieceEnd=std::ldexp(1.,-last.step().potential().geometry().report().acceptedRefinement);
            const auto nodes=midPub.quadrature(0,g,0,pieceEnd);
            check(!nodes.empty(),"actual accepted electric graph quadrature is empty");
            const auto& q=nodes.front().sample;
            const auto a=load.externalSurfaceVector().sample(0,g,q),b=maxwell.sample(0,g,q).tractionJump;
            check(std::hypot(a.x-b.x,a.y-b.y)<=1e-12*(1+std::hypot(b.x,b.y)),"Hydro load retained another electric candidate");
        }
        const auto& electricEnd=step.electricEnd();coeval(mesh,electricEnd,end);
        electricBlock(electricEnd,observed.electricProblem,observed.conductivityByPhase);
        check(electricEnd.step().potential().geometry().geometryIdentity()==h.candidateState().geometryIdentity()
            &&sameClock(electricEnd.sourceTimeLayer(),before.clock),"End potential detached actual final Hydro/FirstJet candidate");
        momentum(mesh,h);characteristic(mesh,h,dt);
        if(acceptedStep==1) {
            auto edited=q0;edited.front()+=initialField.transport().panelMeasure().front()*.125;
            initialField.assignIntegratedCharge(edited);
            rejected([&]{step.commit(mesh,state);},"");
            unchanged(state,before,false);
            check(sameBits(initialField.transport().panelCharge(),edited),"failed joint commit overwrote external Field edit");
            initialField.assignIntegratedCharge(q0);unchanged(state,before);
        }
        step.commit(mesh,state);
        check(state.hydro().sameIdentity(h.candidateState())&&sameClock(state.timeLayer(),end)
            &&sameClock(state.electric().timeLayer(),end)&&sameClock(state.electric().geometrySourceTimeLayer(),end)
            &&sameBits(state.electric().integratedCharge(),electricEnd.step().integratedCharge())
            &&state.electric().previousTimeLayer()&&sameClock(*state.electric().previousTimeLayer(),before.clock)
            &&state.electric().previousIntegratedCharge()
            &&sameBits(*state.electric().previousIntegratedCharge(),electricEnd.mappedIntegratedCharge())
            &&&state.publication().state().geometry().snapshot()==&electricEnd.step().potential().geometry().publication().state().geometry().snapshot(),
            "joint acceptance did not publish one actual End geometry/u/p/Q/history");
        check(!sameBits(state.electric().integratedCharge(),last.step().integratedCharge())
            &&sameBits(initialField.transport().panelCharge(),q0),"accepted physical Q reused stage Q or overwrote old Field");
        const AcceptedCopy committed(state);rejected([&]{step.commit(mesh,state);},"");unchanged(state,committed);
        std::cout<<"EJ1_ACCEPT,RZ,n=4,rho=1000:1,mu=0,step="<<acceptedStep<<",dt="<<dt
            <<",sigma="<<observed.hydro.surfaceTension<<",stages="<<stages.size()<<",physical_history=End_only,PASS\n";
    }
    std::cout<<"EJ1_JOINT_STEP,PASS,accepted=2,observer_failures=2,late_Hydro_failure=1,external_edit_rejection=1,stale_commits=2\n";
    return 0;
}catch(const std::exception& e) {std::cerr<<"EJ1_JOINT_STEP,FAIL,"<<e.what()<<'\n';return 1;}
