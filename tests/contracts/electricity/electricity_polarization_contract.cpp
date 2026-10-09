#include "Modules/Electricity/FittedDielectricPolarization.h"
#include "Modules/Electricity/FittedCellHarmonic.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/IO/DictionaryMesh.h"
#include "Core/Geometry/LevelSetHermiteGeometry.h"
#include "Core/Field/VertField.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <type_traits>
#include <limits>
#include <cstdlib>

namespace {
using namespace Quantum::Electricity;
static_assert(!std::is_assignable_v<decltype(std::declval<const FittedCellHarmonicState&>().potential()[0]),double>);
static_assert(!std::is_assignable_v<decltype(std::declval<const FittedCellHarmonicState&>().epsilonByCell()[0]),double>);
void check(bool ok, const char* reason) { if (!ok) throw std::runtime_error(reason); }
template<class F> void rejects(F action) {
    try { action(); } catch (const std::invalid_argument&) { return; }
    throw std::runtime_error("invalid electric source accepted");
}
void run(bool axisymmetric, bool uniform, bool refined=false) {
    DictionaryMesh input; input.Nx=input.Ny=8; input.maxCellNum=256; input.maxInitialIterations=0;
    QuadtreeMesh mesh(input);
    if(refined) {
        DictionaryMesh::RefineBox box;box.min={0.,.5};box.max={1.,1.};box.level=1;
        mesh.applyRefineBoxes({box});check(mesh.cellIndexNum()==160,"layered test did not create real coarse/fine leaves");
    }
    FaceSegmentGraph graph; graph.build(mesh);
    const auto domain=NativeLeafDomain::build(mesh,graph,axisymmetric);
    const double upper=uniform?1.:2., tangent=(axisymmetric||refined)?0.:.25;
    const auto potential=[&](Vector2D p) { return (p.y<.5?p.y:.5+(p.y-.5)/upper)+tangent*p.x; };
    std::vector<double> epsilon(mesh.cellIndexNum()), coefficient(graph.segments().size());
    for(int k=0;k<mesh.cellIndexNum();++k)
        epsilon[k]=mesh.cells(mesh.cellIndexToLabel[k]).centerCoord().y<.5?1.:upper;
    using Solver=PolarizationPoissonFVSolverCellCentered;
    std::vector<Solver::FittedBoundary> boundary(coefficient.size());
    for(std::size_t s=0;s<coefficient.size();++s) {
        const auto& face=graph.segments()[s]; const int k=mesh.cells(face.K).index();
        coefficient[s]=epsilon[k];
        if(face.L>=0) {
            const int other=mesh.cells(face.L).index();
            const auto a=mesh.cells(face.K).centerCoord(),b=mesh.cells(face.L).centerCoord();
            const double dk=std::abs(face.isHorizontal?face.faceCenter.y-a.y:face.faceCenter.x-a.x);
            const double dl=std::abs(face.isHorizontal?b.y-face.faceCenter.y:b.x-face.faceCenter.x);
            coefficient[s]=(dk+dl)/(dk/epsilon[k]+dl/epsilon[other]);
        } else if(axisymmetric && face.side==FaceSide::West) {
            boundary[s].type=Solver::FittedBoundaryType::HomogeneousNeumann;
        } else {
            boundary[s].type=Solver::FittedBoundaryType::Dirichlet;
            boundary[s].value=potential(face.faceCenter);
        }
    }
    CellField<double> rhs(mesh,0.), phi(mesh,0.);
    Solver solver; Solver::Options options; options.requireSuperLU=true; options.useCG=false; options.useIC=false;
    options.axisymmetric=axisymmetric; options.relTol=1e-11;
    Solver::SolveReport report;
    const bool solved=solver.solveFittedBoundary(mesh,graph,rhs,phi,coefficient,boundary,options,&report);
    if (!solved || !report.success || !report.usedSuperLU)
        throw std::runtime_error("layered electric SuperLU solve failed: "+report.message);
    const FittedCellHarmonicControls controls{1e-11,1e-10,3.,0.,3.*potential({1.,1.})};
    const auto prepared=prepareFittedCellHarmonic(mesh,graph,boundary,epsilon,axisymmetric,0.);
    check(prepared.stage()==FittedCellHarmonicState::Stage::Prepared && !prepared.report().success,
        "unsolved binding was published as a solved field");
    rejects([&]{prepared.faceValues(mesh);});
    auto state=prepareFittedCellHarmonic(mesh,graph,boundary,epsilon,axisymmetric,0.,nullptr,&controls);
    check(state.stage()==FittedCellHarmonicState::Stage::Solved && state.report().success && state.report().usedSuperLU,
        "module did not perform the actual cell-harmonic SuperLU solve");
    check(state.epsilonBySegment()==coefficient,"module harmonic coefficients differ from the layered reference");
    const auto result=buildFittedDielectricPolarization(mesh,domain,graph,state.potential(),epsilon,coefficient,boundary,1e-10);
    for(int k=0;k<mesh.cellIndexNum();++k) {
        check(state.potential()[k]==phi.value[k],"module field differs from the original direct solver");
        check(std::abs(phi.value[k]-potential(mesh.cells(mesh.cellIndexToLabel[k]).centerCoord()))<1e-12,
            "potential disagrees with layered analytic solution");
        check(std::abs(result.fieldX[k]+tangent)<1e-11
            && std::abs(result.fieldY[k]+1./epsilon[k])<1e-11,"reconstructed electric field is incorrect");
    }
    int loaded=0;
    for(const auto& face:domain.faces(mesh)) {
        const bool jump=face.negativeCell>=0 && face.positiveCell>=0
            && epsilon[face.negativeCell]!=epsilon[face.positiveCell];
        if(jump) {
            ++loaded;
            check(std::abs(result.displacementByNative[face.nativeRow]+1.)<1e-11,"normal displacement is incorrect");
            // Independent Maxwell stress jump divided by epsilonHigh-epsilonLow.
            const double stressJump=.5*(1./upper-1.)-.5*tangent*tangent*(upper-1.);
            check(std::abs(result.coefficientByNative[face.nativeRow]-stressJump/(upper-1.))<1e-11,
                "polarization coefficient disagrees with Maxwell stress jump");
        } else check(result.coefficientByNative[face.nativeRow]==0.,"uniform/boundary face has a coefficient");
    }
    check(loaded==(uniform?0:refined?16:8),"wrong number of material-jump faces");
    if(refined) {
        int split=0;
        const auto& offsets=graph.cellOffsets();
        for(int k=0;k<mesh.cellIndexNum();++k)split+=offsets[k+1]-offsets[k]>4;
        check(split==8,"coarse-side split incidence was not exercised");
        std::cout<<"Layered "<<(axisymmetric?"RZ":"XY")<<": 160 cells, 8 split owners, "<<loaded<<" material faces\n";
    }
    const auto faceValues=state.faceValues(mesh);
    for(std::size_t j=0;j<graph.segments().size();++j) {
        const auto& face=graph.segments()[j]; const int k=mesh.cells(face.K).index();
        const double area=axisymmetric?6.283185307179586476925286766559*face.faceCenter.x*face.Af:face.Af;
        const double expected=(-epsilon[k]*tangent*face.n.x-face.n.y)*area;
        check(std::abs(faceValues[j].chargeCode-expected)<1e-11,"face displacement differs from analytic flux");
    }
    auto badEpsilon=epsilon; badEpsilon[0]=0.;
    rejects([&]{updateFittedCellHarmonic(mesh,state,badEpsilon,controls,.001);});
    auto badBounds=controls; badBounds.upperPotential=.01;
    rejects([&]{updateFittedCellHarmonic(mesh,state,epsilon,badBounds,.001);});
    auto changed=epsilon;
    for(int k=0;k<mesh.cellIndexNum();++k) if(mesh.cells(mesh.cellIndexToLabel[k]).centerCoord().y>=.5) changed[k]*=3.;
    VertField<double> scalar(mesh,0.);
    for(int v=0;v<mesh.vertIndexNum();++v) scalar.value[v]=mesh.verts(mesh.vertIndexToLabel[v]).coord().y-.5;
    const auto makeSource=[&](bool axis) {
        if(!refined)return std::make_shared<const LevelSetHermiteGeometry::Snapshot>(mesh,scalar,axis,CartesianGridScope::UniformLeaves);
        Quantum::Geometry::ActiveHermitePatch patch;patch.activeCells.assign(mesh.cellIndexNum(),1);patch.nodes.resize(mesh.vertIndexNum());
        for(int v=0;v<mesh.vertIndexNum();++v)patch.nodes[v]={scalar.value[v],0.,1.,0.};
        return std::make_shared<const LevelSetHermiteGeometry::Snapshot>(LevelSetHermiteGeometry::Snapshot::fromNativeLeaves(mesh,patch,axis));
    };
    auto source=makeSource(axisymmetric);
    const auto updated=updateFittedCellHarmonic(mesh,state,changed,controls,.001,source);
    check(updated.potential().source()==source,"solved field lost its actual retained scalar");
    const auto wrongMeasure=makeSource(!axisymmetric);
    rejects([&]{updateFittedCellHarmonic(mesh,state,changed,controls,.001,wrongMeasure);});
    source.reset();scalar.value[0]+=1.;
    check(updated.potential().source()->value({.5,.25})==-.25,"input mutation changed the retained source scalar");
    check(&updated.graph()==&state.graph() && &updated.boundaries()==&state.boundaries(),"update copied the fixed problem");
    check(updated.time()==.001 && state.time()==0. && state.epsilonByCell()==epsilon,"update changed retained materials/time");
    double difference=0.;
    for(int k=0;k<mesh.cellIndexNum();++k) {
        check(state.potential()[k]==phi.value[k],"update/failure changed the retained field");
        difference=std::max(difference,std::abs(updated.potential()[k]-state.potential()[k]));
    }
    check(difference>1e-4,"material update did not exercise a different field");
    rejects([&]{updateFittedCellHarmonic(mesh,updated,epsilon,controls,0.);});
    const auto retainedPotential=state.potential(); const double first=retainedPotential[0];
    state={};
    check(retainedPotential[0]==first,"discarding state invalidated its owned potential");
    auto wrong=coefficient; wrong.front()*=2.;
    rejects([&]{buildFittedDielectricPolarization(mesh,domain,graph,phi,epsilon,wrong,boundary,1e-10);});
    phi.value[0]+=.1;
    rejects([&]{buildFittedDielectricPolarization(mesh,domain,graph,phi,epsilon,coefficient,boundary,1e-10);});
    phi.value[0]-=.1;
    auto& vertex=mesh.verts(mesh.vertIndexToLabel[0]);
    vertex.setCoord(vertex.coord()+Vector2D{.01,0.});
    rejects([&]{buildFittedDielectricPolarization(mesh,domain,graph,phi,epsilon,coefficient,boundary,1e-10);});
    rejects([&]{updated.faceValues(mesh);});
}
void inactivePayload() {
    DictionaryMesh input;input.Nx=input.Ny=8;input.maxCellNum=256;input.maxInitialIterations=0;
    QuadtreeMesh mesh(input);FaceSegmentGraph graph;
    graph.build(mesh,[](const QuadtreeMesh&,const Cell& c){return c.centerCoord().x>=.25;});
    using Solver=PolarizationPoissonFVSolverCellCentered;
    std::vector<Solver::FittedBoundary> boundary(graph.segments().size());
    std::vector<double> epsilon(mesh.cellIndexNum(),0.);
    for(int k=0;k<mesh.cellIndexNum();++k) if(mesh.cells(mesh.cellIndexToLabel[k]).centerCoord().x>=.25) epsilon[k]=1.;
    for(std::size_t j=0;j<boundary.size();++j) {
        const auto& face=graph.segments()[j];
        if(face.L<0 || epsilon[mesh.cells(face.L).index()]==0.) {
            boundary[j].type=Solver::FittedBoundaryType::Dirichlet;boundary[j].value=2.;
        }
    }
    const FittedCellHarmonicControls controls{1e-11,1e-10,1.,0.,8.};
    const auto state=prepareFittedCellHarmonic(mesh,graph,boundary,epsilon,true,7.,nullptr,&controls);
    for(int k=0;k<mesh.cellIndexNum();++k)
        check(epsilon[k]==0.?state.potential()[k]==7.:std::abs(state.potential()[k]-2.)<1e-12,
            "inactive payload was lost or substituted for a fitted wall voltage");
    CellField<double> invalidLabels(mesh,0.);invalidLabels.setLabel(0,-1);
    const FittedPotentialInput stale(invalidLabels);
    rejects([&]{prepareFittedCellHarmonic(mesh,graph,boundary,epsilon,true,7.,&stale,&controls);});
    CellField<double> missingValues(mesh,0.);missingValues.value.reset();
    const FittedPotentialInput missing(missingValues);
    rejects([&]{prepareFittedCellHarmonic(mesh,graph,boundary,epsilon,true,7.,&missing,&controls);});
    auto bad=boundary;bad.front().value=std::numeric_limits<double>::quiet_NaN();
    rejects([&]{prepareFittedCellHarmonic(mesh,graph,bad,epsilon,true,7.);});
}
void directFactorIdentity(bool superLU) {
    using Solver=PolarizationPoissonFVSolverCellCentered;
    DictionaryMesh input; input.Nx=input.Ny=8; input.maxCellNum=256; input.maxInitialIterations=0;
    QuadtreeMesh mesh(input); FaceSegmentGraph graph; graph.build(mesh);
    PicDomainBC bc; bc.top.type=bc.bottom.type=PicBCType::Neumann;
    CellField<double> rhs(mesh,2.), phi(mesh,0.);
    std::vector<double> kappa(graph.segments().size(),1.);
    Solver solver; Solver::Options opt; opt.useCG=false; opt.useIC=false; opt.useSuperLU=superLU;
    opt.requireSuperLU=superLU; opt.freezeDirectMatrixUntilPatternChange=true;
    const auto solve=[&](bool reused,double coefficient,double source,double voltage) {
        std::fill(kappa.begin(),kappa.end(),coefficient);
        for(int k=0;k<mesh.cellIndexNum();++k)rhs.value[k]=source;
        bc.left.value=bc.right.value=voltage;
        Solver::SolveReport report;
        check(solver.solve(mesh,graph,bc,rhs,phi,&kappa,nullptr,&report,opt),"cached direct solve failed");
        check(report.success && report.usedSuperLU==superLU && report.directMatrixFrozen==reused,
              "direct factors reused for a changed matrix or lost for an unchanged matrix");
        // Independent discrete 1D parabola, including the half-cell Dirichlet closure.
        for(int k=0;k<mesh.cellIndexNum();++k) {
            const double x=mesh.cells(mesh.cellIndexToLabel[k]).centerCoord().x;
            const double expected=voltage+source/(2.*coefficient)*(x*(1.-x)+1./256.);
            if(std::abs(phi.value[k]-expected)>=1e-12)
                std::cerr<<"backend="<<superLU<<" kappa="<<coefficient<<" rhs="<<source
                         <<" cell="<<k<<" x="<<x<<" actual="<<phi.value[k]<<" expected="<<expected<<'\n';
            check(std::abs(phi.value[k]-expected)<1e-12,"cached solve disagrees with discrete Poisson solution");
        }
    };
    solve(false,1.,2.,0.);
    solve(true,1.,6.,0.);  // RHS-only change retains the factors.
    solve(false,4.,6.,0.); // Same sparsity, different coefficients must refactorize.
    solve(true,4.,6.,1.);  // Dirichlet values only change the RHS.
    solve(false,1.,2.,0.);
    solver.clearDirectCache();
    solve(false,1.,2.,0.);
    solve(true,1.,2.,0.);
    // Spatial coefficient changes also preserve sparsity but must invalidate factors.
    kappa[kappa.size()/2]*=3.;
    Solver fresh; CellField<double> expected(mesh,0.); Solver::SolveReport report;
    check(solver.solve(mesh,graph,bc,rhs,phi,&kappa,nullptr,&report,opt) && !report.directMatrixFrozen,
          "local coefficient change reused stale direct factors");
    check(fresh.solve(mesh,graph,bc,rhs,expected,&kappa,nullptr,nullptr,opt),"fresh direct reference failed");
    for(int k=0;k<mesh.cellIndexNum();++k)
        check(std::abs(phi.value[k]-expected.value[k])<1e-12,"cached and fresh variable coefficient solves disagree");
}
}
int main() try {
#ifdef _MSC_VER
    _putenv_s("QUANTUM_POL_SUPERLU_LOCAL","0");
#else
    setenv("QUANTUM_POL_SUPERLU_LOCAL","0",1);
#endif
    for(bool superLU:{false,true}) directFactorIdentity(superLU);
    for(bool axis:{false,true}) for(bool uniform:{false,true}) run(axis,uniform);
    for(bool axis:{false,true})for(bool uniform:{false,true})run(axis,uniform,true);
    inactivePayload();
    std::cout<<"Electric-only cell-harmonic solve/update, analytic XY/RZ fields/flux/Maxwell coefficients, immutable state and inactive payload passed\n";
    return 0;
} catch(const std::exception& e) { std::cerr<<e.what()<<'\n'; return 1; }
