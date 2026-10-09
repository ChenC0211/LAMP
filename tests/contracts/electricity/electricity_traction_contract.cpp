#include "Modules/Electricity/FittedSharpDielectricElectrostatics.h"
#include "Modules/Electricity/SurfaceChargeSource.h"
#include "Core/Geometry/LevelSetInterfacePanels.h"
#include "Core/IO/DictionaryMesh.h"
#include "Core/Field/VertField.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <array>
#include <chrono>
#include <iomanip>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace {
using Snapshot=LevelSetHermiteGeometry::Snapshot;
void check(bool ok,const char* message) { if(!ok) throw std::runtime_error(message); }
template<class F> void rejects(F action) {
    try { action(); } catch(const std::exception&) { return; }
    throw std::runtime_error("invalid electric trace request accepted");
}
template<class F> auto singleThread(F action) {
#ifdef _OPENMP
    struct Restore {int count=omp_get_max_threads();~Restore(){omp_set_num_threads(count);}} restore;
    omp_set_num_threads(1);
#endif
    return action();
}
void sameTraces(const FittedSharpDielectricTractions& a,const FittedSharpDielectricTractions& b) {
    check(a.geometry==b.geometry && a.interfaceTraces.size()==b.interfaceTraces.size(),"parallel trace association changed");
    const auto values=[](const auto& t) {return std::array<double,14>{t.interfacePotential,t.tangentialElectricField,
        t.normalDisplacement,t.tangentialFieldDerivative,t.normalDisplacementDerivative,t.tangentialFieldSecondDerivative,
        t.normalDisplacementSecondDerivative,t.normalTraction,t.condition,t.reproduction,t.positiveNormalDisplacement,
        t.surfaceChargeDensity,t.tangentialTraction,t.modelRelativeResidual};};
    for(std::size_t i=0;i<a.interfaceTraces.size();++i) {
        const auto& x=a.interfaceTraces[i];const auto& y=b.interfaceTraces[i];
        check(values(x)==values(y) && x.sampleCellIndices==y.sampleCellIndices && x.axisCapExteriorModel==y.axisCapExteriorModel,
            "parallel Maxwell fit changed an ordered physical value or stencil");
    }
    const auto diagnostics=[](const auto& d) {return std::array<double,8>{double(d.samplePoints),double(d.localTraceFits),
        d.maxTraceCondition,d.maxTraceReproductionError,d.maxNormalTraction,d.maxTangentialTraction,
        double(d.axisCapExteriorSamples),d.maxAxisCapRelativeResidual};};
    check(diagnostics(a.diagnostics)==diagnostics(b.diagnostics),"parallel Maxwell diagnostics changed");
}
double derivative(const Snapshot& source,const LevelSetHermiteGeometry::Frame& frame) {
    if(source.usesNativeLeaves())
        return source.curvatureTangentialDerivativeAtNativeCellChart(frame.nativeCellIndex,frame.cellChartPoint);
    return frame.hasCellChartPoint
        ?source.curvatureTangentialDerivativeAtCellChart(frame.cellI,frame.cellJ,frame.cellChartPoint)
        :source.curvatureTangentialDerivativeAtCellLocal(frame.cellI,frame.cellJ,frame.cellLocalPoint);
}
void run(bool axisymmetric,bool native) {
    DictionaryMesh input; input.Nx=input.Ny=16;input.maxCellNum=1024;input.maxInitialIterations=0;
    QuadtreeMesh mesh(input);FaceSegmentGraph graph;graph.build(mesh);
    Quantum::Geometry::ActiveHermitePatch patch;
    VertField<double> scalar(mesh,0.);
    patch.activeCells.assign(mesh.cellIndexNum(),1);patch.nodes.resize(mesh.vertIndexNum());
    for(int i=0;i<mesh.vertIndexNum();++i) {
        const auto p=mesh.verts(mesh.vertIndexToLabel[i]).coord();
        const double x=p.x-(axisymmetric?0.:.5),y=p.y-.5;
        patch.nodes[i]={x*x+y*y-.23*.23,2*x,2*y,0.};
        scalar.value[i]=patch.nodes[i].value;
    }
    auto snapshot=std::make_shared<const Snapshot>(native?Snapshot::fromNativeLeaves(mesh,patch,axisymmetric)
        :Snapshot(mesh,scalar,axisymmetric,CartesianGridScope::UniformLeaves));
    auto panels=native?std::make_unique<LevelSetInterfacePanels>(
        LevelSetInterfacePanels::discoverNativeGraphs(mesh,snapshot,axisymmetric))
        :std::make_unique<LevelSetInterfacePanels>(mesh,snapshot,axisymmetric,std::vector<Vector2D>{},CartesianGridScope::UniformLeaves);
    // Exact harmonic potential Phi=3-.2*y in BOTH media, epsilon-=2, epsilon+=1.
    // q=(epsilon+-epsilon-)*E.n=-.2*n_y supplies its consistent local derivatives.
    SurfaceChargeSource charge;charge.axisymmetric=axisymmetric;
    charge.cellIntegratedCharge.assign(mesh.cellIndexNum(),0.);
    charge.evaluate=[snapshot](const auto& frame) {
#ifdef _OPENMP
        check(!omp_in_parallel(),"Maxwell preparation called user charge code in a worker region");
#endif
        const double k=frame.curvature,ks=derivative(*snapshot,frame);
        return SurfaceChargeJet{-.2*frame.normal.y,-.2*k*frame.normal.x,
            -.2*(ks*frame.normal.x-k*k*frame.normal.y)};
    };
    for(const auto& panel:panels->panels())
        for(const auto& q:panel.quadratureNodes)
            charge.cellIntegratedCharge[panel.cellIndex]+=-.2*q.frame.normal.y*q.physicalWeight;
    using Solver=PolarizationPoissonFVSolverCellCentered;
    std::vector<Solver::FittedBoundary> boundary(graph.segments().size());
    for(std::size_t i=0;i<boundary.size();++i) {
        const auto& segment=graph.segments()[i]; if(segment.L>=0)continue;
        boundary[i].type=axisymmetric&&segment.side==FaceSide::West
            ?Solver::FittedBoundaryType::HomogeneousNeumann:Solver::FittedBoundaryType::Dirichlet;
        if(boundary[i].type==Solver::FittedBoundaryType::Dirichlet)boundary[i].value=3.-.2*segment.faceCenter.y;
    }
    DielectricConfig config;config.model="leakyDielectric";config.interfaceDiscretization="sharpHermiteP3";
    config.axisymmetric=axisymmetric;config.permittivityNegative=2.;config.permittivityPositive=1.;
    config.linearResidualTolerance=1e-10;
    CellField<double> guess(mesh,0.);
    FittedSharpDielectricSolver solver;
    const auto solved=native?solver.solveNative(mesh,graph,boundary,snapshot,*panels,config,guess,charge)
        :solver.solve(mesh,graph,boundary,snapshot,*panels,CartesianGridScope::UniformLeaves,config,guess,charge);
    check(solved.report.success&&solved.report.usedSuperLU,"charged analytic field did not use SuperLU");
    for(int c=0;c<mesh.cellIndexNum();++c)
        check(std::abs(solved.potential[c]-(3.-.2*mesh.cells(mesh.cellIndexToLabel[c]).centerCoord().y))<1e-10,
            "charged field disagrees with exact harmonic potential");
    std::vector<InterfaceTraceGeometry::Point> points;
    for(const auto& panel:panels->panels()) {
        const auto projected=panels->projectedFrame(panel.frame.point);
        points.push_back({projected.frame,derivative(*snapshot,projected.frame),projected.kktResidual});
    }
    // Several receiving faces may share a projection. Retain their ordering
    // and use enough requests to exercise the production OpenMP threshold.
    const auto originalPoints=points;
    while(points.size()<96)points.insert(points.end(),originalPoints.begin(),originalPoints.end());
    const auto geometry=std::make_shared<const InterfaceTraceGeometry>(snapshot,points);
    const auto reconstruct=[&](std::shared_ptr<const InterfaceTraceGeometry> request) {
        return native?reconstructNativeFittedSharpDielectricTractions(mesh,request,*panels,config,solved.potential,charge)
            :reconstructFittedSharpDielectricTractions(mesh,request,*panels,CartesianGridScope::UniformLeaves,config,solved.potential,charge);
    };
    const auto begin=std::chrono::steady_clock::now();
    const auto serial=singleThread([&]{return reconstruct(geometry);});
    const auto middle=std::chrono::steady_clock::now();
    const auto result=reconstruct(geometry);
    const auto end=std::chrono::steady_clock::now();
    sameTraces(serial,result);
    std::cout<<"Parallel Maxwell axis="<<axisymmetric<<" native="<<native<<" samples="<<points.size()
        <<" serial_s="<<std::chrono::duration<double>(middle-begin).count()
        <<" parallel_s="<<std::chrono::duration<double>(end-middle).count()<<" exact=1\n";
    check(result.geometry==geometry&&result.interfaceTraces.size()==points.size(),"lost physical sample association");
    const auto tractions=result.tractionComponents();
    double tangentPeak=0.;
    for(std::size_t i=0;i<points.size();++i) {
        const auto n=points[i].frame.normal;
        const auto& trace=result.interfaceTraces[i];
        check(std::abs(trace.normalDisplacement-.4*n.y)<1e-9
            &&std::abs(trace.positiveNormalDisplacement-.2*n.y)<1e-9
            &&std::abs(trace.tangentialElectricField-.2*n.x)<1e-9,"two-sided electric trace is incorrect");
        // Cartesian stress jump is diag(+.02,-.02) for this exact field.
        check(std::abs(tractions[i].x-.02*(n.x*n.x-n.y*n.y))<1e-9
            &&std::abs(tractions[i].y+.04*n.x*n.y)<1e-9,"charged normal/tangential Maxwell stress is incorrect");
        tangentPeak=std::max(tangentPeak,std::abs(tractions[i].y));
    }
    check(tangentPeak>.01,"test did not exercise nonzero tangential traction");
    if(native && axisymmetric) {
        // Independent harmonic-field patch check for the explicit axis-cap
        // closure. The ordinary native path above stays the reference; this
        // manufactured declaration is a model check, not a breakup event.
        FittedAxisCapTractionModel model{snapshot,{}};
        for(const auto& panel:panels->panels())if(snapshot->nativeCellBounds(panel.cellIndex)[0].x==0
            && std::find(model.cells.begin(),model.cells.end(),panel.cellIndex)==model.cells.end())
            model.cells.push_back(panel.cellIndex);
        const auto modeled=reconstructNativeFittedSharpDielectricTractions(mesh,geometry,*panels,config,
            solved.potential,charge,&model);
        sameTraces(singleThread([&]{return reconstructNativeFittedSharpDielectricTractions(mesh,geometry,*panels,
            config,solved.potential,charge,&model);}),modeled);
        check(modeled.diagnostics.axisCapExteriorSamples>0 && modeled.diagnostics.maxAxisCapRelativeResidual<1e-10,
            "axis-cap harmonic control was not reproduced");
        for(std::size_t i=0;i<points.size();++i) {
            const auto& a=modeled.interfaceTraces[i];const auto& b=result.interfaceTraces[i];
            check(std::abs(a.normalTraction-b.normalTraction)<1e-9
                && std::abs(a.tangentialTraction-b.tangentialTraction)<1e-9
                && std::abs(a.tangentialFieldDerivative-b.tangentialFieldDerivative)<1e-8
                && std::abs(a.normalDisplacementDerivative-b.normalDisplacementDerivative)<1e-8
                && std::abs(a.tangentialFieldSecondDerivative-b.tangentialFieldSecondDerivative)<1e-7
                && std::abs(a.normalDisplacementSecondDerivative-b.normalDisplacementSecondDerivative)<1e-7,
                "axis-cap model disagrees with the analytic field or its retained charge transmission derivatives");
        }
        auto stale=model;stale.source.reset();
        rejects([&]{reconstructNativeFittedSharpDielectricTractions(mesh,geometry,*panels,config,solved.potential,charge,&stale);});
        std::cout<<"Axis-cap exterior model: analytic harmonic field/Maxwell/derivatives, declared scope and stale-source rejection passed\n";
    }
    auto invalid=points;invalid[0].curvatureTangentialDerivative+=1.;
    rejects([&]{reconstruct(std::make_shared<const InterfaceTraceGeometry>(snapshot,invalid));});
    invalid=points;invalid[0].frame.normal.x+=.1;
    rejects([&]{reconstruct(std::make_shared<const InterfaceTraceGeometry>(snapshot,invalid));});
    rejects([&]{reconstruct({});});
    invalid=points;invalid[2].frame.normal.x+=.1;invalid[7].frame.normal.x+=.2;
    const auto bad=std::make_shared<const InterfaceTraceGeometry>(snapshot,invalid);
    const auto failure=[&]() {try {reconstruct(bad);}catch(const std::exception& e){return std::string(e.what());}
        throw std::runtime_error("parallel trace preparation accepted invalid frames");};
    const auto serialFailure=singleThread(failure),parallelFailure=failure();
    check(serialFailure==parallelFailure && serialFailure.find("sample-2:")!=std::string::npos,
        "parallel Maxwell failure lost original first-query order/context");
    auto foreign=std::make_shared<const Snapshot>(native?Snapshot::fromNativeLeaves(mesh,patch,axisymmetric)
        :Snapshot(mesh,scalar,axisymmetric,CartesianGridScope::UniformLeaves));
    rejects([&]{reconstruct(std::make_shared<const InterfaceTraceGeometry>(foreign,points));});
}
void chargedSphereAccuracy() {
    // An independently known curved, non-polynomial exterior field. A sphere
    // carrying uniform q has V-=qR, V+=qR^2/rho for epsilon+=1, E-=0 and
    // Maxwell normal jump q^2/2. This tests the actual FV solve plus traction
    // reconstruction; the old affine harmonic test cannot measure this error.
    constexpr double radius=.47,chargeDensity=.001;
    double errors[2]{};
    for(int level=0;level<2;++level) {
        const int n=32<<level;
        DictionaryMesh input;input.minCoord={0,-2};input.maxCoord={2,2};input.Nx=n;input.Ny=2*n;
        input.maxCellNum=4*n*n;input.maxInitialIterations=0;QuadtreeMesh mesh(input);
        FaceSegmentGraph graph;graph.build(mesh);
        Quantum::Geometry::ActiveHermitePatch patch;
        patch.activeCells.assign(mesh.cellIndexNum(),1);patch.nodes.resize(mesh.vertIndexNum());
        for(int v=0;v<mesh.vertIndexNum();++v) {
            const auto p=mesh.verts(mesh.vertIndexToLabel[v]).coord();
            patch.nodes[v]={p.x*p.x+p.y*p.y-radius*radius,2*p.x,2*p.y,0};
        }
        auto snapshot=std::make_shared<const Snapshot>(Snapshot::fromNativeLeaves(mesh,patch,true));
        LevelSetInterfacePanels panels(LevelSetInterfacePanels::discoverNativeGraphs(mesh,snapshot,true));
        SurfaceChargeSource charge;charge.axisymmetric=true;
        charge.cellIntegratedCharge.assign(mesh.cellIndexNum(),0);
        charge.evaluate=[chargeDensity](const auto&){return SurfaceChargeJet{chargeDensity,0,0};};
        for(const auto& panel:panels.panels())for(const auto& q:panel.quadratureNodes)
            charge.cellIntegratedCharge[panel.cellIndex]+=chargeDensity*q.physicalWeight;
        const auto exactPotential=[radius,chargeDensity](Vector2D p) {
            return chargeDensity*radius*radius/std::max(radius,std::hypot(p.x,p.y));
        };
        using Solver=PolarizationPoissonFVSolverCellCentered;
        std::vector<Solver::FittedBoundary> boundary(graph.segments().size());
        for(std::size_t i=0;i<boundary.size();++i) {
            const auto& face=graph.segments()[i];if(face.L>=0)continue;
            boundary[i].type=face.side==FaceSide::West?Solver::FittedBoundaryType::HomogeneousNeumann
                :Solver::FittedBoundaryType::Dirichlet;
            if(boundary[i].type==Solver::FittedBoundaryType::Dirichlet)boundary[i].value=exactPotential(face.faceCenter);
        }
        DielectricConfig config;config.model="leakyDielectric";config.interfaceDiscretization="sharpHermiteP3";
        config.axisymmetric=true;config.permittivityNegative=2;config.permittivityPositive=1;config.linearResidualTolerance=1e-10;
        CellField<double> initial(mesh,0.);
        FittedSharpDielectricSolver solver;
        const auto solved=solver.solveNative(mesh,graph,boundary,snapshot,panels,config,initial,charge);
        check(solved.report.success && solved.report.usedSuperLU,"charged sphere did not complete the actual FV solve");
        FittedAxisCapTractionModel model{snapshot,{}};
        std::vector<InterfaceTraceGeometry::Point> points;std::vector<double> weights;
        for(const auto& panel:panels.panels()) {
            if(snapshot->nativeCellBounds(panel.cellIndex)[0].x==0
                && std::find(model.cells.begin(),model.cells.end(),panel.cellIndex)==model.cells.end())model.cells.push_back(panel.cellIndex);
            for(const auto& q:panel.quadratureNodes) {
                points.push_back({q.frame,derivative(*snapshot,q.frame),0});weights.push_back(q.physicalWeight);
            }
        }
        for(bool end:{false,true}) {
            const auto frame=panels.nativeChain().frameAt(end?static_cast<int>(panels.panels().size())-1:0,end?1:0);
            check(frame.point.x==0,"charged sphere lost its actual axis pole");
            points.push_back({frame,derivative(*snapshot,frame),0});weights.push_back(0);
        }
        const auto queries=std::make_shared<const InterfaceTraceGeometry>(snapshot,points);
        const auto loads=reconstructNativeFittedSharpDielectricTractions(mesh,queries,panels,config,solved.potential,charge,&model);
        const double reference=.5*chargeDensity*chargeDensity;
        double capError=0,capWeight=0,ordinaryError=0,ordinaryWeight=0,poleError=0,tangentialError=0,potentialError=0;
        for(std::size_t i=0;i<points.size();++i) {
            const auto& trace=loads.interfaceTraces[i];const double error=trace.normalTraction/reference-1;
            if(points[i].frame.point.x==0)poleError=std::max(poleError,std::abs(error));
            if(trace.axisCapExteriorModel){capError+=weights[i]*error*error;capWeight+=weights[i];}
            else {ordinaryError+=weights[i]*error*error;ordinaryWeight+=weights[i];}
            tangentialError=std::max(tangentialError,std::abs(trace.tangentialTraction)/reference);
        }
        for(int c=0;c<mesh.cellIndexNum();++c)potentialError=std::max(potentialError,
            std::abs(solved.potential[c]-exactPotential(mesh.cells(mesh.cellIndexToLabel[c]).centerCoord()))/(chargeDensity*radius));
        check(capWeight>0 && ordinaryWeight>0,"charged sphere did not exercise both actual traction providers");
        errors[level]=std::sqrt(capError/capWeight);
        std::cout<<std::setprecision(17)<<"CHARGED_SPHERE N="<<n<<" radius_over_h="<<radius/mesh.cellSize.x
            <<" cap_traction_relative_L2="<<errors[level]<<" ordinary_traction_relative_L2="<<std::sqrt(ordinaryError/ordinaryWeight)
            <<" pole_traction_relative_Linf="<<poleError<<" tangential_scaled_Linf="<<tangentialError
            <<" potential_scaled_Linf="<<potentialError<<" model_fit_residual="<<loads.diagnostics.maxAxisCapRelativeResidual<<'\n';
    }
    check(errors[1]<.01 && errors[1]<errors[0],"resolved charged sphere axis-model accuracy/trend failed");
}
}
void nativeElectricSamplingContract();
int main(int argc,char** argv) try {
    if(argc==2 && std::string(argv[1])=="--native-sampling-study"){nativeElectricSamplingContract();return 0;}
    if(argc==2 && std::string(argv[1])=="--charged-sphere-study"){chargedSphereAccuracy();return 0;}
    nativeElectricSamplingContract();
    for(bool axis:{false,true})for(bool native:{false,true})run(axis,native);
    chargedSphereAccuracy();
    std::cout<<"Electric-only charged XY/RZ and native P3 traces match analytic Maxwell stress; invalid sources rejected\n";
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
