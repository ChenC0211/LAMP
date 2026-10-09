#include "Modules/Electricity/FittedChargeEvolution.h"
#include "Core/IO/DictionaryMesh.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Geometry/ActiveQuarticHermiteTransfer.h"
#include "Core/Geometry/ActiveHermitePatch.h"
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <map>
#include <iomanip>

namespace {
using namespace Quantum::Electricity;
static_assert(std::is_same_v<decltype(std::declval<const FittedChargeState&>().charge()),std::shared_ptr<const SurfaceChargeField>>);
static_assert(!std::is_assignable_v<decltype(std::declval<const FittedChargeState&>().electric().potential[0]),double>);
static_assert(!std::is_assignable_v<decltype(std::declval<const FittedChargeState&>().previousMaterialCharge()[0]),double>);
void check(bool value, const char* message) { if (!value) throw std::runtime_error(message); }
template<class F> void rejects(F action) {
    bool rejected=false;
    try { action(); } catch (const std::exception&) { rejected=true; }
    check(rejected,"invalid/incomplete fitted charge candidate was accepted");
}
std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> sphere(const QuadtreeMesh& mesh) {
    VertField<double> scalar(mesh,0.);
    for (int i=0;i<mesh.vertIndexNum();++i) {
        const auto p=mesh.verts(mesh.vertIndexToLabel[i]).coord();
        scalar.value[i]=p.x*p.x+(p.y-.5)*(p.y-.5)-.23*.23;
    }
    return std::make_shared<const LevelSetHermiteGeometry::Snapshot>(
        mesh,scalar,true,CartesianGridScope::UniformLeaves);
}
void multipleNativeSurfaces() {
    DictionaryMesh input;input.minCoord={0,0};input.maxCoord={1,1};input.Nx=32;input.Ny=64;
    input.maxCellNum=5000;input.maxInitialIterations=0;QuadtreeMesh mesh(input);
    constexpr double radius=.17;
    const std::array<double,2> centers{{.25,.75}},densities{{.001,-.0006}};
    Quantum::Geometry::ActiveHermitePatch patch;
    patch.activeCells.assign(mesh.cellIndexNum(),1);patch.nodes.resize(mesh.vertIndexNum());
    for(int v=0;v<mesh.vertIndexNum();++v) {
        const auto p=mesh.verts(mesh.vertIndexToLabel[v]).coord();
        const double z=p.y-centers[p.y>.5?1:0];
        patch.nodes[v]={p.x*p.x+z*z-radius*radius,2*p.x,p.y==.5?0:2*z,0};
    }
    auto geometry=std::make_shared<const LevelSetHermiteGeometry::Snapshot>(
        LevelSetHermiteGeometry::Snapshot::fromNativeLeaves(mesh,patch,true));
    auto components=LevelSetInterfacePanels::discoverNativeComponents(mesh,geometry,true);
    check(components.curves.size()==2,"two charged spheres lost their separate curves");
    std::vector<std::shared_ptr<SurfaceChargeField>> fields;
    std::vector<const LevelSetInterfacePanels*> surfaces;
    std::vector<SurfaceChargeSource> sources;std::vector<int> owner(mesh.cellIndexNum(),-1);
    SurfaceChargeSource charge;charge.axisymmetric=true;charge.cellIntegratedCharge.assign(mesh.cellIndexNum(),0);
    for(auto& chain:components.curves) {
        const int side=chain.frameAt(0,.5).point.y>.5?1:0;
        fields.push_back(std::make_shared<SurfaceChargeField>(mesh,std::move(chain),densities[side]));
        const int index=static_cast<int>(fields.size())-1;
        surfaces.push_back(&fields.back()->geometry());sources.push_back(fields.back()->source());
        for(const auto& panel:surfaces.back()->panels()) {
            check(owner[panel.cellIndex]<0 || owner[panel.cellIndex]==index,"sphere fixture shares a physical cut cell");
            owner[panel.cellIndex]=index;
        }
        for(std::size_t c=0;c<charge.cellIntegratedCharge.size();++c)
            charge.cellIntegratedCharge[c]+=sources.back().cellIntegratedCharge[c];
    }
    charge.evaluate=[sources,owner](const auto& frame) {
        const int side=owner.at(frame.nativeCellIndex);check(side>=0,"global field asked an unowned charge frame");
        return sources.at(side).evaluate(frame);
    };
    const auto exact=[&](Vector2D p) {
        double value=0;for(int s=0;s<2;++s)
            value+=densities[s]*radius*radius/std::max(radius,std::hypot(p.x,p.y-centers[s]));
        return value;
    };
    FaceSegmentGraph graph;graph.build(mesh);using Solver=PolarizationPoissonFVSolverCellCentered;
    std::vector<Solver::FittedBoundary> boundary(graph.segments().size());
    for(std::size_t i=0;i<boundary.size();++i) {
        const auto& face=graph.segments()[i];if(face.L>=0)continue;
        boundary[i].type=face.side==FaceSide::West?Solver::FittedBoundaryType::HomogeneousNeumann:Solver::FittedBoundaryType::Dirichlet;
        if(boundary[i].type==Solver::FittedBoundaryType::Dirichlet)boundary[i].value=exact(face.faceCenter);
    }
    DielectricConfig config;config.model="leakyDielectric";config.interfaceDiscretization="sharpHermiteP3";
    config.axisymmetric=true;config.permittivityNegative=config.permittivityPositive=1;
    config.lowerPotential=config.upperPotential=0;config.linearResidualTolerance=1e-10;
    const DielectricConductivity conductivity{0,.1};CellField<double> initial(mesh,0.);
    FittedSharpDielectricSolver solver;
    rejects([&]{solver.solveNative(mesh,graph,boundary,geometry,
        std::vector<const LevelSetInterfacePanels*>{surfaces[0]},config,initial,charge,&conductivity);});
    rejects([&]{solver.solveNative(mesh,graph,boundary,geometry,
        std::vector<const LevelSetInterfacePanels*>{surfaces[0],surfaces[0]},config,initial,charge,&conductivity);});
    const auto both=solver.solveNative(mesh,graph,boundary,geometry,surfaces,config,initial,charge,&conductivity);
    check(both.report.success && both.report.usedSuperLU && both.report.conductionCurrent,
        "two surfaces did not complete one global charged/current solve");
    double error=0,scale=0;
    for(int c=0;c<mesh.cellIndexNum();++c) {
        const double reference=exact(mesh.cells(mesh.cellIndexToLabel[c]).centerCoord());
        error=std::max(error,std::abs(both.potential[c]-reference));scale=std::max(scale,std::abs(reference));
    }
    check(error/scale<.03,"two-sphere global potential differs from the analytic superposition");
    const auto repeated=solver.solveNative(mesh,graph,boundary,geometry,surfaces,config,initial,charge,&conductivity);
    check(repeated.geometryAssemblyReused,"complete multi-interface geometry was not reused on an unchanged solve");
    std::reverse(surfaces.begin(),surfaces.end());
    const auto reordered=solver.solveNative(mesh,graph,boundary,geometry,surfaces,config,initial,charge,&conductivity);
    double orderError=0,reuseError=0;
    for(int c=0;c<mesh.cellIndexNum();++c) {
        orderError=std::max(orderError,std::abs(both.potential[c]-reordered.potential[c]));
        reuseError=std::max(reuseError,std::abs(both.potential[c]-repeated.potential[c]));
    }
    check(orderError/scale<1e-10 && reuseError/scale<1e-10,"component ordering or geometry reuse changed the global solution");
    double stressError=0,stressScale=0;std::size_t stressPoints=0;
    for(const auto& field:fields) {
        const auto& panels=field->geometry();
        const int side=panels.panels().front().frame.point.y>.5?1:0,other=1-side;
        std::vector<InterfaceTraceGeometry::Point> points;std::vector<double> weights;
        FittedAxisCapTractionModel model{geometry,{}};
        for(const auto& panel:panels.panels()) {
            if(geometry->nativeCellBounds(panel.cellIndex)[0].x==0
                && std::find(model.cells.begin(),model.cells.end(),panel.cellIndex)==model.cells.end())model.cells.push_back(panel.cellIndex);
            for(const auto& q:panel.quadratureNodes) {
                points.push_back({q.frame,geometry->curvatureTangentialDerivativeAtNativeCellChart(
                    q.frame.nativeCellIndex,q.frame.cellChartPoint),0});weights.push_back(q.physicalWeight);
            }
        }
        const auto queries=std::make_shared<const InterfaceTraceGeometry>(geometry,points);
        // Both loads consume exactly the same globally solved potential and
        // combined charge. Only their physical interface query sets differ.
        const auto loads=reconstructNativeFittedSharpDielectricTractions(mesh,queries,panels,config,
            both.potential,charge,&model).tractionComponents();
        for(std::size_t k=0;k<points.size();++k) {
            const auto& frame=points[k].frame;const auto p=frame.point,n=frame.normal;
            const double z=p.y-centers[other],distance=std::hypot(p.x,z);
            const double magnitude=densities[other]*radius*radius/(distance*distance*distance);
            const double er=magnitude*p.x,ez=magnitude*z;
            const double normal=densities[side]*(er*n.x+ez*n.y)+.5*densities[side]*densities[side];
            const double tangent=densities[side]*(-er*n.y+ez*n.x);
            stressError+=weights[k]*(std::pow(loads[k].x-normal,2)+std::pow(loads[k].y-tangent,2));
            stressScale+=weights[k]*(normal*normal+tangent*tangent);++stressPoints;
        }
    }
    const double relativeStress=std::sqrt(stressError/stressScale);
    check(relativeStress<.05,"child Maxwell loads disagree with the interacting-shell analytic stresses");
    FittedChargeProblem problem{graph,boundary,{config,conductivity,1e-8}};
    std::vector<FittedChargeComponentInitial> initialComponents;
    for(std::size_t s=0;s<fields.size();++s) {
        const int side=fields[s]->geometry().panels().front().frame.point.y>.5?1:0;
        initialComponents.push_back({static_cast<std::uint64_t>(side+2),&fields[s]->geometry(),densities[side]});
    }
    auto state=initializeFittedChargeState(mesh,problem,geometry,initialComponents,initial);
    check(state.componentCount()==2,"module did not retain both child charge states");
    rejects([&]{state.charge();});rejects([&]{state.previousMaterialCharge();});
    const auto retained=state;const auto retainedSource=state.source();
    const auto retainedQ0=state.charge(0)->transport().panelCharge(),retainedQ1=state.charge(1)->transport().panelCharge();
    const auto retainedFrame=state.charge(0)->geometry().panels().front().frame;
    const auto retainedDensity=retainedSource.evaluate(retainedFrame).density;
    double maximumChildBalance=0,maximumGlobalBalance=0;
    for(int step=1;step<=2;++step) {
        std::vector<FittedChargeComponentMotion> motions;
        for(std::size_t s=0;s<state.componentCount();++s)
            motions.push_back({state.surfaceId(s),{geometry,&state.charge(s)->geometry(),{},
                [](Vector2D p){return p;},.001,state.time(),step*.001}});
        // Reordering is not a surface identity change.
        std::reverse(motions.begin(),motions.end());
        if(step==1) {
            auto bad=motions;bad.back().surfaceId=bad.front().surfaceId;
            rejects([&]{prepareFittedChargeCollectionStep(mesh,problem,state,bad);});
            bad=motions;bad.back().motion.forwardMaterialPoint=[](Vector2D)->Vector2D{throw std::runtime_error("second child failure");};
            rejects([&]{prepareFittedChargeCollectionStep(mesh,problem,state,bad);});
            check(state.completedStep()==0 && state.charge(0)->transport().panelCharge()==retainedQ0
                && state.charge(1)->transport().panelCharge()==retainedQ1,"failed child partially committed the collection");
        }
        auto next=prepareFittedChargeCollectionStep(mesh,problem,state,motions);
        for(std::size_t s=0;s<next.state.componentCount();++s) {
            std::size_t old=0;while(state.surfaceId(old)!=next.state.surfaceId(s))++old;
            long double expected=0,actual=0,reference=0;
            for(double q:state.charge(old)->transport().panelCharge()){expected+=q;reference+=std::abs(q);}
            const auto& cells=state.charge(old)->geometry().nativeChain().cellToGraphs();
            for(std::size_t c=0;c<cells.size();++c)if(!cells[c].empty())
                expected+=.001L*state.electric().report.conductionCurrent->cellIntegratedCurrent[c];
            for(double q:next.state.charge(s)->transport().panelCharge())actual+=q;
            maximumChildBalance=std::max(maximumChildBalance,static_cast<double>(std::abs(actual-expected)/reference));
        }
        maximumGlobalBalance=std::max(maximumGlobalBalance,next.report.relativeChargeClosure);
        check(next.state.completedStep()==step && next.state.time()==step*.001,"collection did not advance one common clock");
        state=std::move(next.state);
    }
    check(maximumChildBalance<1e-12,"old current was duplicated or assigned to the wrong child");
    check(retained.completedStep()==0 && retained.charge(0)->transport().panelCharge()==retainedQ0
        && retained.charge(1)->transport().panelCharge()==retainedQ1
        && retainedSource.evaluate(retainedFrame).density==retainedDensity,"new charge time level changed retained state/source");
    std::cout<<std::setprecision(17)<<"MULTIPLE_CHARGE surfaces=2 global_Poisson=1 analytic_relative_max="<<error/scale
        <<" component_order_error="<<orderError/scale<<" reuse_error="<<reuseError/scale
        <<" current_segments="<<both.conductionCurrentBySegment.size()<<" missing_and_duplicate_rejected=1"
        <<" Maxwell_points="<<stressPoints<<" coupled_stress_relative_L2="<<relativeStress
        <<" material_steps=2 child_current_relative_balance="<<maximumChildBalance
        <<" global_electrode_relative_balance="<<maximumGlobalBalance<<" atomic_failure_and_retained_source=1"<<'\n';
}
void steepArcMapping() {
    DictionaryMesh input;input.Nx=input.Ny=16;input.maxCellNum=1024;input.maxInitialIterations=0;
    QuadtreeMesh mesh(input);std::size_t refined=0;
    // The nearly tangent endpoint also needs a well-conditioned chart when
    // one globally signed partial is much smaller than the other. Keep the
    // analytic whole-circle/sphere measures and material-charge checks.
    for(bool axis:{false,true})for(double gap:{.007,.001,.000001}) {
        constexpr double radius=.04,density=.001;
        const double center=.5-gap,radialCenter=axis?0:.5;
        Quantum::Geometry::ActiveQuarticHermitePatch patch;
        patch.activeCells.assign(mesh.cellIndexNum(),1);patch.nodes.resize(mesh.vertIndexNum());
        patch.cellCenters.resize(mesh.cellIndexNum());
        const auto value=[&](Vector2D p){return (p.x-radialCenter)*(p.x-radialCenter)+(p.y-center)*(p.y-center)-radius*radius;};
        for(int v=0;v<mesh.vertIndexNum();++v) {
            const auto p=mesh.verts(mesh.vertIndexToLabel[v]).coord();const double z=p.y-center,r=p.x-radialCenter;
            patch.nodes[v]={r*r+z*z-radius*radius,2*r,2*z,0};
        }
        constexpr int pairs[4][2]={{0,1},{3,2},{0,3},{1,2}};
        std::map<std::array<int,2>,bool> sides;
        for(int c=0;c<mesh.cellIndexNum();++c) {
            const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);patch.cellCenters[c]=value(cell.centerCoord());
            for(int side=0;side<4;++side)sides.emplace(std::array<int,2>{mesh.verts(cell.cellToVert[pairs[side][0]]).index(),
                mesh.verts(cell.cellToVert[pairs[side][1]]).index()},side<2);
        }
        for(const auto& side:sides) {
            const auto a=mesh.verts(mesh.vertIndexToLabel[side.first[0]]).coord();
            const auto b=mesh.verts(mesh.vertIndexToLabel[side.first[1]]).coord();const Vector2D p{(a.x+b.x)*.5,(a.y+b.y)*.5};
            patch.edgeMidpoints.push_back({side.first,value(p),side.second?2*(p.y-center):2*(p.x-radialCenter)});
        }
        auto geometry=std::make_shared<const LevelSetHermiteGeometry::Snapshot>(
            LevelSetHermiteGeometry::Snapshot::fromNativeLeaves(mesh,patch,axis));
        auto components=LevelSetInterfacePanels::discoverNativeComponents(mesh,geometry,axis);
        check(components.curves.size()==1,"steep sphere did not retain one complete meridian");
        SurfaceChargeField field(mesh,std::move(components.curves[0]),density);
        const auto& chain=field.geometry().nativeChain();const auto source=field.source();
        const auto before=field.transport().panelCharge();
        double length=0,area=0,maxEstimate=0;std::size_t intervals=0;
        for(const auto& panel:field.geometry().panels()) {
            length+=panel.measure;area+=panel.physicalMeasure;
            refined+=panel.measureQuadratureIntervals>1;
            intervals=std::max(intervals,panel.measureQuadratureIntervals);
            for(double e:panel.measureQuadratureRelativeError)maxEstimate=std::max(maxEstimate,e);
            double previousLength=0,previousArea=0;
            for(const auto& q:panel.quadratureNodes) {
                check(q.localArclength>previousLength && q.localPhysicalMeasure>previousArea
                    && q.localArclength<panel.measure && q.localPhysicalMeasure<panel.physicalMeasure,
                    "native cumulative measures are not strictly within their panel");
                previousLength=q.localArclength;previousArea=q.localPhysicalMeasure;
            }
        }
        const double lengthError=std::abs(length/((axis?1:2)*std::acos(-1.)*radius)-1);
        const double areaError=std::abs(area/((axis?4*radius:2)*std::acos(-1.)*radius)-1);
        check(lengthError<1e-8 && areaError<1e-8 && maxEstimate<=1e-10,
            "adaptive native sphere measures disagree with the independent analytic integrals");
        for(int g=0;g<static_cast<int>(chain.graphs().size());++g)for(int k=0;k<=64;++k) {
            const auto q=source.evaluate(chain.frameAt(g,k/64.));
            check(q.density==density && q.firstDerivative==0 && q.secondDerivative==0,
                "steep arc mapping changed a constant charge jet");
        }
        const auto& report=field.measureInterpolation();
        check(before==field.transport().panelCharge(),"measure interpolation changed panel-integrated charge");
        auto varying=before;
        for(std::size_t i=0;i<varying.size();++i)varying[i]*=1+.1*i;
        field.assignIntegratedCharge(varying);const auto nonuniform=field.source();double momentError=0;
        for(std::size_t i=0;i<varying.size();++i) {
            double integral=0;
            for(const auto& q:field.geometry().panels()[i].quadratureNodes)
                integral+=nonuniform.evaluate(q.frame).density*q.physicalWeight;
            momentError=std::max(momentError,std::abs(integral-varying[i])/std::abs(varying[i]));
        }
        check(momentError<1e-12,"rational map changed the saved charge reconstruction moments");
        if(axis) {
            const auto remap=field.remapMoving(mesh,geometry,{},[](Vector2D p){return p;},{varying},&field.geometry());
            double error=0;
            for(std::size_t i=0;i<varying.size();++i)
                error=std::max(error,std::abs(remap.field->transport().panelCharge()[i]-varying[i])/std::abs(varying[i]));
            std::cout<<std::setprecision(17)<<"STEEP_ARC identity_panel_relative_error="<<error
                <<" charge_closure="<<remap.report.remapConservationResidual
                <<" projection="<<remap.report.maxProjectionDistance<<'\n';
            check(error<1e-12 && remap.history.at(0)==remap.field->transport().panelCharge()
                && field.transport().panelCharge()==varying,
                "rational measure map changed identity material charge/history beyond roundoff");
        }
        std::cout<<std::setprecision(17)<<"STEEP_ARC axis="<<axis<<" gap="<<gap<<" panels="<<chain.graphs().size()
            <<" rational_arc="<<report.rationalArcIntervals<<" rational_area="<<report.rationalAreaIntervals
            <<" max_measure_intervals="<<intervals<<" max_measure_error_estimate="<<maxEstimate
            <<" saved_length_relative_error="<<lengthError
            <<" saved_area_relative_error="<<areaError
            <<" charge_moment_relative_error="<<momentError<<'\n';
    }
    check(refined>0,"steep native graph did not exercise adaptive surface integration");
}
FittedChargeState initial(const QuadtreeMesh& mesh, const FittedChargeProblem& problem,
    std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> geometry) {
    const LevelSetInterfacePanels panels(mesh,geometry,true,{},CartesianGridScope::UniformLeaves);
    CellField<double> guess(mesh,0.);
    auto invalidSeed=guess; invalidSeed.value.reset();
    rejects([&]{initializeFittedChargeState(mesh,problem,geometry,panels,{},invalidSeed);});
    auto invalidLabels=guess; invalidLabels.setLabel(0,-1);
    rejects([&]{initializeFittedChargeState(mesh,problem,geometry,panels,{},invalidLabels);});
    for(bool flux:{true,false}) {
        auto unused=problem;
        if(flux) unused.parameters.field.boundaryFluxTolerance=1e-9;
        else unused.parameters.field.chargeClosureTolerance=1e-9;
        rejects([&]{initializeFittedChargeState(mesh,unused,geometry,panels,{},guess);});
    }
    auto state=initializeFittedChargeState(mesh,problem,geometry,panels,{},guess);
    check(state.electric().report.success && state.electric().report.usedSuperLU,"real initial P3/SuperLU solve failed");
    for(int i=0;i<mesh.cellIndexNum();++i)
        check(std::abs(state.electric().potential[i]-.01*mesh.cells(mesh.cellIndexToLabel[i]).centerCoord().y)<1e-12,
            "uniform-permittivity initial field differs from analytic potential");
    return state;
}
void run(const QuadtreeMesh& mesh,const FittedChargeProblem& problem,bool implicit) {
    const auto geometry=sphere(mesh);
    auto state=initial(mesh,problem,geometry);
    // Identity material motion still remaps both histories and rebuilds the
    // charged field; no mock current or relaxation callback supplies physics.
    FittedChargeMotion motion{geometry,&state.charge()->geometry(),{},[](Vector2D p){return p;},.001,0.,.001};
    SurfaceChargeRelaxationOptions options;
    options.underRelaxation=.1; options.maxIterations=500;
    options.densityChangeTolerance=1e-10; options.relativeEquationTolerance=1e-10;
    auto oldQ=state.charge()->transport().panelCharge();
    std::vector<double> oldPotential(mesh.cellIndexNum());
    for(int i=0;i<mesh.cellIndexNum();++i) oldPotential[i]=state.electric().potential[i];
    const auto unchanged=[&] {
        check(oldQ==state.charge()->transport().panelCharge(),"failure changed accepted charge");
        check(state.previousMaterialCharge().empty() && state.previousChargeDt()==0.,"failure published history");
        for(int i=0;i<mesh.cellIndexNum();++i)
            check(oldPotential[i]==state.electric().potential[i],"failure changed accepted potential");
    };
    auto bad=motion; bad.dt=0.;
    rejects([&]{prepareFittedChargeStep(mesh,problem,state,bad);}); unchanged();
    bad=motion; bad.sourceTime=.001;
    rejects([&]{prepareFittedChargeStep(mesh,problem,state,bad);}); unchanged();
    bad=motion; bad.targetTime=.002;
    rejects([&]{prepareFittedChargeStep(mesh,problem,state,bad);}); unchanged();
    bad=motion; bad.geometry=sphere(mesh);
    rejects([&]{prepareFittedChargeStep(mesh,problem,state,bad);}); unchanged();
    bad=motion; bad.forwardMaterialPoint=[](Vector2D)->Vector2D{throw std::runtime_error("injected map failure");};
    rejects([&]{prepareFittedChargeStep(mesh,problem,state,bad);}); unchanged();
    auto shortIteration=options; shortIteration.maxIterations=1;
    rejects([&]{prepareFittedChargeStep(mesh,problem,state,motion,&shortIteration);}); unchanged();
    rejects([&]{prepareFittedChargeStep(mesh,problem,state,motion,&options,
        [](int,int,bool,bool){throw std::runtime_error("injected observer failure");});}); unchanged();

    for(int step=1;step<=3;++step) {
        const auto retained=state;
        const auto retainedQ=retained.charge()->transport().panelCharge();
        const auto retainedHistory=retained.previousMaterialCharge();
        std::vector<double> retainedPotential(mesh.cellIndexNum());
        for(int i=0;i<mesh.cellIndexNum();++i) retainedPotential[i]=retained.electric().potential[i];
        const auto source=retained.source();
        const auto frame=retained.charge()->geometry().projectedFrame(retained.charge()->geometry().panels().front().frame.point).frame;
        const auto jet=source.evaluate(frame);
        motion.panels=&state.charge()->geometry();
        motion.sourceTime=state.time(); motion.targetTime=step*motion.dt;
        int observations=0;
        auto next=prepareFittedChargeStep(mesh,problem,state,motion,implicit ? &options : nullptr,
            [&](int count,int limit,bool,bool){check(count==++observations && count<=limit+1,"invalid observation count");});
        check(next.chargeStep.transport.success && next.chargeStep.remapped,"material update incomplete");
        check(next.state.electric().report.success && next.state.electric().report.usedSuperLU,"P3 field missing");
        check(next.report.relativeChargeClosure<=problem.parameters.closureTolerance,"electrode closure failed");
        if(implicit) {
            check(next.report.relaxation && next.report.relaxation->success,"implicit iteration incomplete");
            check(next.report.relaxation->temporalOrder==(step==1?1:2),"history time order changed");
            check(next.state.previousMaterialCharge()==state.charge()->transport().panelCharge(),"identity history remap changed integrals");
            check(next.report.electricFactorReuses>0 && next.report.electricGeometryReuses>0,"production caches not exercised");
            check(observations>1,"implicit field was not iterated");
        } else {
            check(!next.report.relaxation && observations==0,"explicit path performed hidden iteration");
            const double residual=next.chargeStep.transport.totalChargeAfter-next.chargeStep.transport.totalChargeBefore
                +motion.dt*state.electric().report.conductionCurrent->netBoundaryCurrent;
            check(residual==next.report.chargeClosureResidual,"explicit closure did not use old current");
        }
        double peak=0.; for(double q:next.state.charge()->transport().density()) peak=std::max(peak,std::abs(q));
        check(peak>0. && std::isfinite(peak),"zero/invalid charge does not test conduction");
        state=std::move(next.state);
        check(state.completedStep()==step && state.time()==motion.targetTime,"candidate time publication is incorrect");
        check(retained.completedStep()==step-1 && retained.charge()->transport().panelCharge()==retainedQ
            && retained.previousMaterialCharge()==retainedHistory,"commit changed retained charge/history");
        for(int i=0;i<mesh.cellIndexNum();++i)
            check(retained.electric().potential[i]==retainedPotential[i],"commit changed retained field");
        const auto after=source.evaluate(frame);
        check(after.density==jet.density && after.firstDerivative==jet.firstDerivative
            && after.secondDerivative==jet.secondDerivative,"commit changed retained charge source");
    }
    const auto source=state.source();
    const auto frame=state.charge()->geometry().projectedFrame(state.charge()->geometry().panels().front().frame.point).frame;
    const auto expected=source.evaluate(frame);
    const auto potential=state.electric().potential;
    const double first=potential[0];
    auto scratch=potential.copyForMeshTransfer(); scratch.value[0]+=1.;
    check(potential[0]==first,"AMR interpolation copy aliases the immutable potential");
    state={};
    check(source.evaluate(frame).density==expected.density && potential[0]==first,
        "destroying state invalidated retained source/potential");
}
void meshTransfer(const QuadtreeMesh& mesh,const FittedChargeProblem& problem) {
    using Patch=Quantum::Geometry::ActiveQuarticHermitePatch;
    using Snapshot=LevelSetHermiteGeometry::Snapshot;
    Patch patch; patch.activeCells.assign(mesh.cellIndexNum(),1);
    patch.nodes.resize(mesh.vertIndexNum()); patch.cellCenters.resize(mesh.cellIndexNum());
    const auto value=[](Vector2D p){return p.x*p.x+(p.y-.5)*(p.y-.5)-.23*.23;};
    for(int v=0;v<mesh.vertIndexNum();++v) {
        const auto p=mesh.verts(mesh.vertIndexToLabel[v]).coord();
        patch.nodes[v]={value(p),2*p.x,2*(p.y-.5),0.};
    }
    constexpr int pairs[4][2]={{0,1},{3,2},{0,3},{1,2}};
    std::map<std::array<int,2>,bool> sides;
    for(int c=0;c<mesh.cellIndexNum();++c) {
        const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);
        patch.cellCenters[c]=value(cell.centerCoord());
        for(int side=0;side<4;++side)
            sides.emplace(std::array<int,2>{mesh.verts(cell.cellToVert[pairs[side][0]]).index(),
                mesh.verts(cell.cellToVert[pairs[side][1]]).index()},side<2);
    }
    for(const auto& side:sides) {
        const auto a=mesh.verts(mesh.vertIndexToLabel[side.first[0]]).coord();
        const auto b=mesh.verts(mesh.vertIndexToLabel[side.first[1]]).coord();
        const Vector2D p{(a.x+b.x)*.5,(a.y+b.y)*.5};
        patch.edgeMidpoints.push_back({side.first,value(p),side.second?2*(p.y-.5):2*p.x});
    }
    const auto geometry=std::make_shared<const Snapshot>(Snapshot::fromNativeLeaves(mesh,patch,true));
    const LevelSetInterfacePanels panels(LevelSetInterfacePanels::discoverNativeGraphs(mesh,geometry,true));
    CellField<double> guess(mesh,0.);
    auto invalidSeed=guess; invalidSeed.value.reset();
    rejects([&]{initializeFittedChargeState(mesh,problem,geometry,panels,{},invalidSeed);});
    auto invalidLabels=guess; invalidLabels.setLabel(0,-1);
    rejects([&]{initializeFittedChargeState(mesh,problem,geometry,panels,{},invalidLabels);});
    auto state=initializeFittedChargeState(mesh,problem,geometry,panels,{},guess);
    SurfaceChargeRelaxationOptions options; options.underRelaxation=.1; options.maxIterations=500;
    options.densityChangeTolerance=options.relativeEquationTolerance=1e-10;
    FittedChargeMotion motion{geometry,&state.charge()->geometry(),{},[](Vector2D p){return p;},.001,0.,.001};
    state=prepareFittedChargeStep(mesh,problem,state,motion,&options).state;
    motion={geometry,&state.charge()->geometry(),{},[](Vector2D p){return p;},.001,state.time(),.002};
    state=prepareFittedChargeStep(mesh,problem,state,motion,&options).state;
    bool nonzeroHistory=false;
    for(double q:state.previousMaterialCharge()) nonzeroHistory=nonzeroHistory || q!=0.;
    check(nonzeroHistory,"mesh transfer must preserve a nonzero previous material layer");
    auto target=mesh.clone();
    DictionaryMesh::RefineBox box; box.min={.75,0.}; box.max={1.,1.}; box.level=1;
    target->applyRefineBoxes({box});
    check(target->cellIndexNum()>mesh.cellIndexNum(),"mesh-transfer test did not change actual topology");
    std::cout<<"Native Q4 charge/history transfer: "<<mesh.cellIndexNum()<<" -> "<<target->cellIndexNum()<<" cells\n";
    const auto mapped=ActiveQuarticHermiteTransfer::prepare(mesh,patch,*target,
        std::vector<std::uint8_t>(target->cellIndexNum(),1),true);
    const auto nextGeometry=std::make_shared<const Snapshot>(Snapshot::fromNativeLeaves(*target,mapped.patch,true));
    const auto& cellMap=mapped.sourceToTargetCell;
    auto badMap=cellMap; badMap[state.charge()->geometry().panels().front().cellIndex]=-1;
    rejects([&]{prepareFittedChargeMeshTransfer(mesh,state,*target,nextGeometry,badMap);});
    const auto directTransfer=prepareFittedChargeMeshTransfer(mesh,state,*target,nextGeometry,cellMap);
    const auto preparedPanels=state.charge()->geometry().rebindNative(
        LevelSetInterfacePanels::discoverNativeGraphs(*target,nextGeometry,true),cellMap);
    rejects([&]{prepareFittedChargeMeshTransfer(mesh,state,*target,nextGeometry,cellMap,&panels);});
    const auto transfer=prepareFittedChargeMeshTransfer(mesh,state,*target,nextGeometry,cellMap,&preparedPanels);
    check(transfer.charge()->transport().panelCharge()==directTransfer.charge()->transport().panelCharge(),
        "prebuilt Interface geometry changed transferred charge");
    FaceSegmentGraph graph; graph.build(*target);
    using Solver=PolarizationPoissonFVSolverCellCentered;
    std::vector<Solver::FittedBoundary> boundaries(graph.segments().size());
    for(std::size_t i=0;i<boundaries.size();++i) {
        const auto& face=graph.segments()[i]; if(face.L>=0) continue;
        boundaries[i].type=Solver::FittedBoundaryType::HomogeneousNeumann;
        if(face.boundary==Boundary::south || face.boundary==Boundary::north) {
            boundaries[i].type=Solver::FittedBoundaryType::Dirichlet;
            boundaries[i].value=face.boundary==Boundary::north ? .01 : 0.;
        }
    }
    FittedChargeProblem nextProblem{graph,boundaries,problem.parameters};
    CellField<double> seed(*target,0.);
    auto invalid=nextProblem; invalid.parameters.conductivity.negative=std::numeric_limits<double>::quiet_NaN();
    rejects([&]{solveFittedChargeMeshTransfer(*target,invalid,transfer,nextGeometry,seed);});
    const auto rebound=solveFittedChargeMeshTransfer(*target,nextProblem,transfer,nextGeometry,seed);
    check(rebound.time()==state.time() && rebound.completedStep()==state.completedStep()
        && rebound.previousChargeDt()==state.previousChargeDt()
        && rebound.previousMaterialCharge()==state.previousMaterialCharge()
        && rebound.charge()->transport().panelCharge()==state.charge()->transport().panelCharge(),
        "zero-time native Q4 rebind changed charge/history/time");
    check(rebound.electric().potential.source()==nextGeometry
        && rebound.charge()->geometry().usesSource(*nextGeometry),"rebind did not publish a coherent target source");
    motion={nextGeometry,&rebound.charge()->geometry(),{},[](Vector2D p){return p;},.001,rebound.time(),.003};
    const auto advanced=prepareFittedChargeStep(*target,nextProblem,rebound,motion,&options);
    check(advanced.report.relaxation->temporalOrder==2 && advanced.state.completedStep()==3,
        "mesh rebind lost the two-level material history");
}
}
void axisTransmissionContract();
int main(int argc,char** argv) try {
    if(argc==2 && std::string(argv[1])=="--axis-transmission-study"){axisTransmissionContract();return 0;}
    if(argc==2 && std::string(argv[1])=="--multiple-surfaces-study"){multipleNativeSurfaces();return 0;}
    if(argc==2 && std::string(argv[1])=="--arc-map-study"){steepArcMapping();return 0;}
    axisTransmissionContract();
    multipleNativeSurfaces();
    steepArcMapping();
    DictionaryMesh input; input.Nx=input.Ny=16; input.maxCellNum=1024; input.maxInitialIterations=0;
    QuadtreeMesh mesh(input); FaceSegmentGraph graph; graph.build(mesh);
    using Solver=PolarizationPoissonFVSolverCellCentered;
    std::vector<Solver::FittedBoundary> boundaries(graph.segments().size());
    for(std::size_t i=0;i<boundaries.size();++i) {
        const auto& face=graph.segments()[i];
        if(face.L>=0) continue;
        boundaries[i].type=Solver::FittedBoundaryType::HomogeneousNeumann;
        if(face.boundary==Boundary::south || face.boundary==Boundary::north) {
            boundaries[i].type=Solver::FittedBoundaryType::Dirichlet;
            boundaries[i].value=face.boundary==Boundary::north ? .01 : 0.;
        }
    }
    FittedChargeParameters p; p.field.model="leakyDielectric"; p.field.interfaceDiscretization="sharpHermiteP3";
    p.field.axisymmetric=true; p.field.permittivityNegative=p.field.permittivityPositive=1.;
    p.field.linearResidualTolerance=1e-10; p.conductivity={.1,.01}; p.closureTolerance=1.6e-9;
    const FittedChargeProblem problem{graph,boundaries,p};
    rejects([&]{preflightFittedChargeStep(p,3.);});
    auto bad=p; bad.conductivity.negative=std::numeric_limits<double>::quiet_NaN();
    rejects([&]{preflightFittedChargeStep(bad,.001);});
    run(mesh,problem,false); run(mesh,problem,true); meshTransfer(mesh,problem);
    std::cout<<"Fitted Electricity alone: immutable field/charge/history, owned source lifetime, explicit and BE/two-level charge, native Q4 mesh rebind passed\n";
    return 0;
} catch(const std::exception& error) { std::cerr<<error.what()<<'\n'; return 1; }
