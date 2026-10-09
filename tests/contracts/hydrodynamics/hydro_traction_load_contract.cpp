#include "Modules/Hydrodynamics/HydroFittedInterfaceTractionLoad.h"
#include "Modules/Hydrodynamics/HydroFittedDomain.h"
#include "Modules/Hydrodynamics/HydroFittedCapillary.h"
#include "Modules/Hydrodynamics/HydroFittedForce.h"
#include "Modules/Interface/Compatibility/InterfacePrescribedData.h"
#include "Modules/Interface/Compatibility/InterfaceState.h"
#include "Core/IO/DictionaryMesh.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/FaceSegmentGraph.h"
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace {
using namespace Quantum::Hydrodynamics;
void check(bool ok,const char* reason){if(!ok)throw std::runtime_error(reason);}
template<class F>void rejects(F action){
    try{action();}catch(const std::invalid_argument&){return;}
    throw std::runtime_error("invalid physical traction/source was accepted");
}
void run(bool native){
    DictionaryMesh input;input.Nx=32;input.Ny=64;input.maxCoord={2.,4.};input.maxCellNum=8192;input.maxInitialIterations=0;
    QuadtreeMesh mesh(input);FaceSegmentGraph graph;
    graph.build(mesh,[](const QuadtreeMesh&,const Cell& cell){const auto p=cell.centerCoord();return p.x<.5||p.y>2.;});
    using Type=HydroFittedDomain::BoundaryType;
    std::vector<HydroFittedDomain::Boundary> boundaries(graph.segments().size());
    for(std::size_t i=0;i<boundaries.size();++i){
        const auto& face=graph.segments()[i];
        const bool active=face.L>=0&&(mesh.cells(face.L).centerCoord().x<.5||mesh.cells(face.L).centerCoord().y>2.);
        if(active)continue;
        boundaries[i].type=face.L>=0?Type::NoSlipWall:face.side==FaceSide::West?Type::Symmetry:Type::PressureOutlet;
    }
    const auto domain=HydroFittedDomain::build(mesh,graph,boundaries,true);
    using Boundary=InterfaceTransportDomain::Boundary;
    InterfaceTransportDomain transport(domain.geometry(),[&](int row){
        const auto type=domain.faces()[row].boundary.type;
        return type==Type::Interior?Boundary::Interior:type==Type::NoSlipWall?Boundary::StationaryWall
            :type==Type::Symmetry?Boundary::Symmetry:Boundary::SinglePhasePort;
    });
    InterfaceState initial(mesh);InterfacePrescribedData::initializeParabolicCap(mesh,transport,{.5,2.},.5,initial);
    auto snapshot=initial.geometry();
    std::unique_ptr<LevelSetInterfacePanels> nativePanels;
    const auto* panels=initial.panels();
    if(native){
        snapshot=std::make_shared<const LevelSetHermiteGeometry::Snapshot>(
            LevelSetHermiteGeometry::Snapshot::fromNativeLeaves(mesh,std::get<Quantum::Geometry::ActiveHermitePatch>(initial.scalarPatch()),true));
        nativePanels=std::make_unique<LevelSetInterfacePanels>(LevelSetInterfacePanels::discoverNativeGraphs(mesh,snapshot,true,{{.5,2.}}));
        panels=nativePanels.get();
    }
    const auto geometry=native?buildNativeFittedInterfaceLoadGeometry(mesh,domain,snapshot,*panels,3.*.0625,
        FittedTransverseGradientClosure::ActiveHermiteWallGradient)
        :buildFittedInterfaceLoadGeometry(mesh,domain,snapshot,*panels,3.*.0625,FittedTransverseGradientClosure::ActiveHermiteWallGradient);
    const auto& faces=geometry.faces();check(!faces.empty(),"no physical load faces");
    std::vector<Vector2D> traction(faces.size(),{2.,-.75});
    const auto load=buildFittedInterfaceTractionLoad(mesh,domain,geometry,geometry.traceGeometry(),traction);
    Eigen::VectorXd h(domain.cellIndexByRow().size());
    for(std::size_t i=0;i<domain.cellIndexByRow().size();++i)h[i]=geometry.indicator[domain.cellIndexByRow()[i]];
    const Eigen::VectorXd pressure=domain.pressureGradientLoad()*h;
    HydroFittedCapillaryState capillary(mesh);
    constexpr double sigma=.375;
    assembleFittedCapillaryForce(mesh,domain,geometry,sigma,capillary);
    long double capX=0.,capY=0.;
    for(const auto& face:faces) {
        const int row=face.edgeIndex+(face.horizontal?mesh.vEdgeIndexNum():0);
        const auto& receiving=domain.faces()[row];
        // Analytic RZ curvature of z=2.5-2*r^2, including its r=0 limit.
        const double r=face.frame.point.x,square=1.+16.*r*r;
        const double curvature=4./std::pow(square,1.5)+4./std::sqrt(square);
        const double expected=-sigma*curvature*(pressure[receiving.velocityRow]/receiving.dualVolume);
        const double actual=face.horizontal?capillary.forceY.value[face.edgeIndex]:capillary.forceX.value[face.edgeIndex];
        check(std::abs(actual-expected)<2e-9,"capillary load differs from analytic curvature times actual pressure gradient");
        (face.horizontal?capY:capX)+=static_cast<long double>(expected)*face.dualMeasure;
    }
    check(std::abs(capillary.integralX-static_cast<double>(capX))<2e-9
        &&std::abs(capillary.integralY-static_cast<double>(capY))<2e-9,"capillary RZ load integrals disagree");
    if(native) {
        const auto rebuiltCap=buildNativeFittedCapillaryForce(mesh,domain,snapshot,*panels,geometry.halfWidth,sigma);
        for(int i=0;i<mesh.vEdgeIndexNum();++i)check(rebuiltCap.forceX.value[i]==capillary.forceX.value[i],"native capillary radial path differs");
        for(int i=0;i<mesh.hEdgeIndexNum();++i)check(rebuiltCap.forceY.value[i]==capillary.forceY.value[i],"native capillary axial path differs");
    } else rejects([&]{buildNativeFittedCapillaryForce(mesh,domain,snapshot,*panels,geometry.halfWidth,sigma);});
    // The completed load owns its values and provenance across workspace reuse.
    rejects([&]{HydroFittedForce::combine(mesh,domain,capillary,load.forceX,load.forceY);});
    capillary.source=HydroForceSource::fitted(domain,.125,HydroForceEvolution::TimeDependent,snapshot);
    auto combined=HydroFittedForce::combine(mesh,domain,capillary,load.forceX,load.forceY,2.);
    const auto retained=combined;
    const auto density=retained.density();
    check(density.source==capillary.source && density.evaluationTime==.125
        && density.faceX==combined.density().faceX && density.faceY==combined.density().faceY,
        "completed load lost its source or copied its owned arrays");
    for(const auto& face:faces) {
        const int i=face.edgeIndex;
        check(face.horizontal ? retained.y()[i]==capillary.forceY.value[i]+2.*load.forceY.value[i]
            : retained.x()[i]==capillary.forceX.value[i]+2.*load.forceX.value[i],
            "completed load differs from capillary plus converted traction");
    }
    int freeX=-1,fixedX=-1;
    for(const auto& face:domain.faces()) if(!face.horizontal) {
        if(face.velocityRow>=0)freeX=face.nativeRow; else fixedX=face.nativeRow;
    }
    check(freeX>=0 && fixedX>=0,"force rejection fixture needs free and constrained rows");
    const auto& first=faces.front();
    const auto retainedValue=first.horizontal?retained.y()[first.edgeIndex]:retained.x()[first.edgeIndex];
    check(retainedValue!=0.,"retained load fixture needs a nonzero value");
    auto external=load.forceX;
    external.value[freeX]=std::numeric_limits<double>::infinity();
    rejects([&]{HydroFittedForce::combine(mesh,domain,capillary,external,load.forceY);});
    external.value[freeX]=std::numeric_limits<double>::max();
    rejects([&]{HydroFittedForce::combine(mesh,domain,capillary,external,load.forceY,2.);});
    external.value[freeX]=std::numeric_limits<double>::denorm_min();
    rejects([&]{HydroFittedForce::combine(mesh,domain,capillary,external,load.forceY,.25);});
    external=load.forceX;external.value[fixedX]=1.;
    rejects([&]{HydroFittedForce::combine(mesh,domain,capillary,external,load.forceY);});
    external.value.reset();
    rejects([&]{HydroFittedForce::combine(mesh,domain,capillary,external,load.forceY);});
    rejects([&]{HydroFittedForce::combine(mesh,domain,capillary,load.forceX,load.forceY,0.);});
    rejects([&]{assembleFittedCapillaryForce(mesh,domain,geometry,-1.,capillary);});
    assembleFittedCapillaryForce(mesh,domain,geometry,0.,capillary);
    combined=HydroFittedForce{};
    check(!combined.source() && (first.horizontal?retained.y()[first.edgeIndex]:retained.x()[first.edgeIndex])==retainedValue
        && retained.source()==density.source,
        "reusing component/owner workspace changed a retained completed load");
    for(int i=0;i<mesh.vEdgeIndexNum();++i)check(capillary.forceX.value[i]==0.,"zero tension retained old radial load");
    for(int i=0;i<mesh.hEdgeIndexNum();++i)check(capillary.forceY.value[i]==0.,"zero tension retained old axial load");
    long double ix=0.,iy=0.;int transverseOnly=0;
    for(const auto& face:faces){
        const int row=face.edgeIndex+(face.horizontal?mesh.vEdgeIndexNum():0);
        const auto& receiving=domain.faces()[row];
        const double normal=2.*pressure[receiving.velocityRow]/receiving.dualVolume;
        const double tangent=-.75*face.transverseGradient;
        const double expected=normal+(face.horizontal?tangent:-tangent);
        const double actual=face.horizontal?load.forceY.value[face.edgeIndex]:load.forceX.value[face.edgeIndex];
        check(std::abs(actual-expected)<1e-11,"load differs from pressure-compatible normal plus rotated transverse components");
        if(face.gradient==0.&&face.transverseGradient!=0.){++transverseOnly;check(actual!=0.,"transverse-only force was dropped");}
        (face.horizontal?iy:ix)+=static_cast<long double>(expected)*face.dualMeasure;
    }
    check(transverseOnly>0,"transverse-only support was not exercised");
    check(std::abs(load.totalForceX-static_cast<double>(ix))<1e-11&&std::abs(load.totalForceY-static_cast<double>(iy))<1e-11,
        "physical RZ load integrals disagree");
    const auto copiedDomain=domain;
    const auto copied=buildFittedInterfaceTractionLoad(mesh,copiedDomain,geometry,geometry.traceGeometry(),traction);
    check(copied.totalForceX==load.totalForceX,"copying the immutable domain changed the load");
    const auto foreign=std::make_shared<const InterfaceTraceGeometry>(snapshot,geometry.traceGeometry()->points());
    rejects([&]{buildFittedInterfaceTractionLoad(mesh,domain,geometry,foreign,traction);});
    const auto rebuilt=HydroFittedDomain::build(mesh,graph,boundaries,true);
    capillary.source=density.source;
    rejects([&]{HydroFittedForce::combine(mesh,rebuilt,capillary,load.forceX,load.forceY);});
    rejects([&]{assembleFittedCapillaryForce(mesh,rebuilt,geometry,sigma,capillary);});
    rejects([&]{buildFittedInterfaceTractionLoad(mesh,rebuilt,geometry,geometry.traceGeometry(),traction);});
    auto bad=traction;bad.pop_back();
    rejects([&]{buildFittedInterfaceTractionLoad(mesh,domain,geometry,geometry.traceGeometry(),bad);});
    bad=traction;bad[0].y=std::numeric_limits<double>::quiet_NaN();
    rejects([&]{buildFittedInterfaceTractionLoad(mesh,domain,geometry,geometry.traceGeometry(),bad);});
    const auto repeat=buildFittedInterfaceTractionLoad(mesh,domain,geometry,geometry.traceGeometry(),traction);
    check(repeat.totalForceX==load.totalForceX&&repeat.totalForceY==load.totalForceY,"rejected load modified retained data");
    auto& vertex=mesh.verts(mesh.vertIndexToLabel[0]);vertex.setCoord(vertex.coord()+Vector2D{.01,0.});
    rejects([&]{buildFittedInterfaceTractionLoad(mesh,domain,geometry,geometry.traceGeometry(),traction);});
}
}
int main()try{
    run(false);run(true);
    std::cout<<"Hydro-only normal/transverse support, native branch, RZ integrals, owned load combination and source rejection passed\n";
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
