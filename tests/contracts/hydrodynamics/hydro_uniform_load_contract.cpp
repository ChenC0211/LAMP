#include "Modules/Hydrodynamics/HydroUniformLoad.h"
#include "Modules/Hydrodynamics/HydroInterfaceForce.h"
#include "Modules/Hydrodynamics/HydroInterfaceLoadGeometry.h"
#include "Core/Mesh/NativeMeshIdentity.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/FaceSegmentGraph.h"
#include "Core/Mesh/CartesianMeasure.h"
#include "Core/IO/DictionaryMesh.h"
#include "Core/Field/VertField.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <type_traits>

namespace {
using namespace Quantum::Hydrodynamics;
void check(bool ok, const char* reason) { if (!ok) throw std::runtime_error(reason); }
template<class F> void rejects(F action) {
    bool rejected=false;
    try { action(); } catch(const std::exception&) { rejected=true; }
    check(rejected,"invalid load/source was accepted");
}
void assemblyContract(QuadtreeMesh& mesh) {
    using Assembly=HydroInterfaceForce::Assembly;
    static_assert(!std::is_copy_constructible_v<Assembly> && !std::is_move_constructible_v<Assembly>);
    HydroInterfaceForce::Report report; std::string stage="caller-stage";
    auto fresh=[&] {return std::make_unique<Assembly>(mesh,nullptr,-1,false,false,
        HydroInterfaceForce::Options{},report,stage);};
    auto empty=fresh(); rejects([&]{empty->faceX();}); rejects([&]{empty->faceY();});
    const auto noLoad=empty->publish(0.);
    check(!noLoad.cellX && !noLoad.cellY && !noLoad.faceX && !noLoad.faceY && !noLoad.source,
        "empty assembly manufactured force fields or source identity");
    check(stage=="caller-stage","empty assembly changed the failure stage");
    rejects([&]{empty->publish(0.);});
    CellField<double> x(mesh,2.),y(mesh,-3.);
    VEdgeField<double> fx(mesh,.25); HEdgeField<double> fy(mesh,-.5);
    auto combined=fresh();
    rejects([&]{combined->addCouplingFaces(fx,fy,true);});
    combined->addCouplingCells(x,y);
    rejects([&]{combined->publish(.001);});
    rejects([&]{combined->addCouplingCells(x,y);});
    combined->addCouplingFaces(fx,fy,true);
    rejects([&]{combined->addCouplingFaces(fx,fy,true);});
    rejects([&]{combined->publish(std::numeric_limits<double>::quiet_NaN());});
    const auto load=combined->publish(.125);
    check(load.cellX && load.cellY && load.faceX && load.faceY && load.evaluationTime==.125 && !load.source,
        "load publication lost additive sources or invented provenance");
    x.value[0]=99.;y.value[0]=99.;fx.value[0]=99.;fy.value[0]=99.;
    for(int i=0;i<mesh.cellIndexNum();++i) check(load.cellX->value[i]==2. && load.cellY->value[i]==-3.,"combined cells alias provider storage");
    for(int i=0;i<mesh.vEdgeIndexNum();++i) check(load.faceX->value[i]==.25,"combined x faces alias provider storage");
    for(int i=0;i<mesh.hEdgeIndexNum();++i) check(load.faceY->value[i]==-.5,"combined y faces alias provider storage");
    rejects([&]{combined->addCouplingCells(x,y);}); rejects([&]{combined->publish(.125);});
    auto cellOnly=fresh(); cellOnly->addCouplingCells(x,y);
    fx.initialize(mesh,0.);fy.initialize(mesh,-0.);cellOnly->addCouplingFaces(fx,fy,false);
    const auto cells=cellOnly->publish(0.);
    check(cells.cellX && cells.cellY && !cells.faceX && !cells.faceY,"cell-only provider enabled face loading");
    for(bool horizontal:{false,true}) {
        auto invalid=fresh(); invalid->addCouplingCells(x,y);
        fx.initialize(mesh,0.);fy.initialize(mesh,0.);
        if(horizontal) fy.value[fy.size()-1]=1.;else fx.value[fx.size()-1]=1.;
        rejects([&]{invalid->addCouplingFaces(fx,fy,false);});
        rejects([&]{invalid->publish(0.);});
        fx.initialize(mesh,0.);fy.initialize(mesh,0.);
        rejects([&]{invalid->addCouplingFaces(fx,fy,true);});
    }
    auto badCell=fresh();y.value[y.size()-1]=std::numeric_limits<double>::infinity();
    rejects([&]{badCell->addCouplingCells(x,y);});rejects([&]{badCell->publish(0.);});
    rejects([&]{badCell->faceX();});rejects([&]{badCell->faceY();});
    y.initialize(mesh,0.);rejects([&]{badCell->addCouplingCells(x,y);});
    auto badSize=fresh();CellField<double> missing;
    rejects([&]{badSize->addCouplingCells(missing,y);});rejects([&]{badSize->publish(0.);});
    auto badFace=fresh();badFace->addCouplingCells(x,y);fy.value[0]=std::numeric_limits<double>::quiet_NaN();
    rejects([&]{badFace->addCouplingFaces(fx,fy,true);});rejects([&]{badFace->publish(0.);});
    auto badFaceSize=fresh();badFaceSize->addCouplingCells(x,y);VEdgeField<double> missingFace;
    rejects([&]{badFaceSize->addCouplingFaces(missingFace,fy,true);});rejects([&]{badFaceSize->publish(0.);});
    HydroInterfaceForce::Options active;active.surfaceTension=.001;
    rejects([&]{Assembly missingGeometry(mesh,nullptr,-1,false,false,active,report,stage);});
}
void scalarGradient(QuadtreeMesh& mesh) {
    NativeMeshIdentity source(mesh);
    FaceSegmentGraph graph; graph.build(mesh);
    std::vector<double> scalar(mesh.cellIndexNum()), coefficient(mesh.vEdgeIndexNum()+mesh.hEdgeIndexNum(),0.);
    for(int i=0;i<mesh.cellIndexNum();++i) {
        const auto p=mesh.cells(mesh.cellIndexToLabel[i]).centerCoord();
        scalar[i]=p.x+2.*p.y;
    }
    int boundaryRow=-1;
    for(const auto& f:graph.segments()) {
        const int row=f.isHorizontal?mesh.vEdgeIndexNum()+mesh.hEdges(f.edgeLabel).index():mesh.vEdges(f.edgeLabel).index();
        if(f.boundary==Boundary::non) coefficient[row]=3.; else boundaryRow=row;
    }
    const auto load=buildUniformScalarGradientLoad(mesh,source,scalar,coefficient);
    for(const auto& f:graph.segments()) {
        const double value=f.isHorizontal?load.faceForceY.value[mesh.hEdges(f.edgeLabel).index()]
            :load.faceForceX.value[mesh.vEdges(f.edgeLabel).index()];
        const double expected=f.boundary==Boundary::non?(f.isHorizontal?6.:3.):0.;
        check(value==expected,"linear scalar gradient did not load the actual free face");
    }
    check(std::abs(load.totalForceX-3.*15./16.)<1e-13 && std::abs(load.totalForceY-6.*15./16.)<1e-13,
        "linear-gradient integral includes constrained boundary velocity rows");
    auto wrong=coefficient; wrong[boundaryRow]=1.;
    rejects([&]{buildUniformScalarGradientLoad(mesh,source,scalar,wrong);});
    wrong=coefficient;wrong.pop_back();
    rejects([&]{buildUniformScalarGradientLoad(mesh,source,scalar,wrong);});
    wrong=coefficient;wrong.front()=std::numeric_limits<double>::quiet_NaN();
    rejects([&]{buildUniformScalarGradientLoad(mesh,source,scalar,wrong);});
    auto& vertex=mesh.verts(mesh.vertIndexToLabel[0]);
    const auto saved=vertex.coord();vertex.setCoord(saved+Vector2D{.01,0.});
    rejects([&]{buildUniformScalarGradientLoad(mesh,source,scalar,coefficient);});
    vertex.setCoord(saved);
    std::fill(scalar.begin(),scalar.end(),7.);
    const auto zero=buildUniformScalarGradientLoad(mesh,source,scalar,coefficient);
    check(zero.maxForce==0. && zero.totalForceX==0. && zero.totalForceY==0.,"constant scalar produced a load");
}
void interfaceTraction(QuadtreeMesh& mesh, bool axisymmetric, bool tangential) {
    VertField<double> scalar(mesh,0.);
    for(int i=0;i<mesh.vertIndexNum();++i) {
        const auto p=mesh.verts(mesh.vertIndexToLabel[i]).coord();
        scalar.value[i]=std::hypot(p.x-.5,p.y-.5)-.23;
    }
    HydroInterfaceLoadGeometry::Options options;
    options.layer=HydroInterfaceLoadGeometry::Layer::Regularized;
    options.halfWidth=.08;options.includeTangentialLoads=tangential;options.axisymmetric=axisymmetric;
    HydroInterfaceLoadGeometry geometry(mesh,scalar,options);
    HydroInterfaceForce::Report zeroReport;std::string zeroStage="unchanged";
    HydroInterfaceForce::Assembly zero(mesh,nullptr,-1,axisymmetric,false,{},zeroReport,zeroStage,&geometry);
    const auto zeroLoad=zero.publish(0.);
    check(zeroLoad.cellX && zeroLoad.cellY && zeroLoad.faceX && zeroLoad.faceY,
        "explicit Hermite zero load changed its field presence");
    check(zeroStage=="unchanged" && zeroReport.nodeCount==0,"zero Hermite load built a reference contour");
    for(int i=0;i<mesh.cellIndexNum();++i)check(zeroLoad.cellX->value[i]==0. && zeroLoad.cellY->value[i]==0.,"nonzero Hermite zero cell load");
    for(int i=0;i<mesh.vEdgeIndexNum();++i)check(zeroLoad.faceX->value[i]==0.,"nonzero Hermite zero face load");
    std::vector<Vector2D> traction(geometry.faces().size(),{2.,tangential?-.75:0.});
    const auto query=geometry.traceGeometry();
    const auto load=buildUniformInterfaceTractionLoad(mesh,scalar,geometry,query,traction,tangential);
    long double integralX=0.,integralY=0.;
    bool transverseOnly=false;
    for(const auto& f:geometry.faces()) {
        // Independent pressure-gradient reference from retained cell H.
        const double gradient=(geometry.loadIndicator()[f.highCellIndex]-geometry.loadIndicator()[f.lowCellIndex])/.0625;
        const double expected=2.*gradient+(tangential?(f.horizontal?-1.:1.)*.75*f.transverseGradient:0.);
        const double actual=f.horizontal?load.faceForceY.value[f.edgeIndex]:load.faceForceX.value[f.edgeIndex];
        check(std::abs(actual-expected)<1e-12,"interface traction differs from pressure-gradient reference");
        (f.horizontal?integralY:integralX)+=static_cast<long double>(expected)*f.dualMeasure;
        if(f.gradient==0. && f.transverseGradient!=0.) { transverseOnly=true;check(actual!=0.,"transverse-only load was lost"); }
    }
    check(!tangential||transverseOnly,"test did not exercise transverse-only support");
    check(std::abs(load.totalForceX-static_cast<double>(integralX))<1e-11
        && std::abs(load.totalForceY-static_cast<double>(integralY))<1e-11,"cell diagnostics disagree with dual-load integral");
    HydroInterfaceLoadGeometry other(mesh,scalar,options);
    rejects([&]{buildUniformInterfaceTractionLoad(mesh,scalar,geometry,other.traceGeometry(),traction,tangential);});
    auto wrong=traction;wrong.pop_back();
    rejects([&]{buildUniformInterfaceTractionLoad(mesh,scalar,geometry,query,wrong,tangential);});
    wrong=traction;wrong.front().x=std::numeric_limits<double>::infinity();
    rejects([&]{buildUniformInterfaceTractionLoad(mesh,scalar,geometry,query,wrong,tangential);});
    if(tangential) rejects([&]{buildUniformInterfaceTractionLoad(mesh,scalar,geometry,query,traction,false);});
    const double saved=scalar.value[0];scalar.value[0]+=.01;
    rejects([&]{buildUniformInterfaceTractionLoad(mesh,scalar,geometry,query,traction,tangential);});
    scalar.value[0]=saved;
    HydroInterfaceLoadGeometry moved(std::move(geometry));
    const auto repeated=buildUniformInterfaceTractionLoad(mesh,scalar,moved,query,traction,tangential);
    check(repeated.totalForceX==load.totalForceX,"moving load geometry changed source association");
    std::shared_ptr<const InterfaceTraceGeometry> retained;
    { HydroInterfaceLoadGeometry temporary(mesh,scalar,options);retained=temporary.traceGeometry(); }
    check(std::isfinite(retained->source()->value(retained->points().front().frame.point)),"query source did not retain its lifetime");
}
}
int main() try {
    DictionaryMesh input;input.Nx=input.Ny=16;input.maxCellNum=1024;input.maxInitialIterations=0;
    QuadtreeMesh mesh(input);
    assemblyContract(mesh);
    scalarGradient(mesh);
    for(bool rz:{false,true})for(bool tangent:{false,true})interfaceTraction(mesh,rz,tangent);
    std::cout<<"Hydro uniform loads: analytic gradient, XY/RZ normal/tangential traction, source and lifetime checks passed\n";
    return 0;
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}
